!> @file

!> Module to handle DEC-MP2 gradient construction, see the mp2grad structure for details.
!> The theory is given in JCP 137, 114102 (2012). 
!> \author Kasper Kristensen

!> Module to handle DEC-MP2 gradient construction, see the mp2grad structure for details.
module rimp2_gradient_module

!  use fundamental
  use precision
  use tensor_type_def_module
  use dec_typedef_module
  use lsparameters
#ifdef VAR_MPI
  use lsmpi_type
  use infpar_module
#endif
  ! DEC DEPENDENCIES (within deccc directory)   
  ! *****************************************!
  use ri_util_module
#ifdef VAR_MPI
  use decmpi_module
#endif
  public :: RIMP2_gradient_driver
  private

contains
  subroutine RIMP2_gradient_driver(MyFragment,ThetaOCC,RIMP2grad,natoms,&
       & nbasis,noccEOS,nvirtAOS,CoccEOS,CvirtAOS,dopair_occ)
    implicit none
    integer,intent(in) :: natoms,nbasis,noccEOS,nvirtAOS
    !> Atomic fragment
    type(decfrag),intent(inout) :: MyFragment
    !> The MO Theta intermediate (intent in but due to bcast inout)
    real(realk),intent(inout)   :: ThetaOcc(nvirtAOS*(noccEOS*i8*noccEOS)*nvirtAOS)
    !> MP2 gradient
    real(realk) :: RIMP2grad(3,natoms)
    !> Occupied MO coef
    real(realk) :: CoccEOS(nbasis,noccEOS)
    !> Vitual MO coef
    real(realk) :: CvirtAOS(nbasis,nvirtAOS)
    !> Do the Pair of occupied indexes 
    logical :: dopair_occ(noccEOS,noccEOS)
    !
    logical :: wakeslave,CollaborateWithSlaves,master,FORCEPRINT
    integer :: LUPRI,mynum,numnodes,nAtomsAux,nBasis2
    integer :: nBasisAux,J,I,Jeos
    integer(kind=long) :: nsize
    LUPRI = DECinfo%output
    FORCEPRINT = .FALSE.
#ifdef VAR_MPI
    master= (infpar%lg_mynum == infpar%master)
    if(infpar%lg_nodtot.GT.1) then
       wakeslave=.true.
       CollaborateWithSlaves=.true.
    else
       wakeslave=.false.
       CollaborateWithSlaves=.false.
    end if
    StartUpSlavesRIMP2grad: if(wakeslave .and. master) then
       call ls_mpibcast(DECRIMP2GRAD,infpar%master,infpar%lg_comm)
       call mpi_communicate_MyFragment(MyFragment)
       nsize = nvirtAOS*(noccEOS*i8*nvirtAOS)*noccEOS
       call ls_mpibcast(ThetaOcc,nsize,infpar%master,infpar%lg_comm)
       call ls_mpibcast(dopair_occ,noccEOS,noccEOS,infpar%master,infpar%lg_comm)
    endif StartUpSlavesRIMP2grad
    IF(.NOT.master) LUPRI = 6 !standard Output
    mynum = infpar%lg_mynum
    numnodes = infpar%lg_nodtot
#else
    master=.true.
    wakeslave=.false.
    CollaborateWithSlaves=.false.
    mynum = 0
    numnodes = 1
#endif    
    if(natoms.NE.MyFragment%natoms)call lsquit('Error in RIMP2 natoms1 dim mismatch',-1)
    call getMolecularDimensions(MyFragment%mylsitem%SETTING%MOLECULE(1)%p,nAtomsAux,nBasis2,nBasisAux)
    if(nbasis.NE.nbasis2)call lsquit('Error in RIMP2 nbasis dim mismatch',-1)
    if(natoms.NE.natomsAux)call lsquit('Error in RIMP2 natoms dim mismatch',-1)

    call Build_RIMP2grad(MyFragment%myLSitem,master,nbasis,nbasisAux,LUPRI,FORCEPRINT,&
         & CollaborateWithSlaves,CvirtAOS,nvirtAOS,CoccEOS,noccEOS,mynum,numnodes,nAtomsAux,&
         & natoms,ThetaOcc,RIMP2grad,dopair_occ)

  end subroutine RIMP2_gradient_driver

end module rimp2_gradient_module

#ifdef VAR_MPI
subroutine RIMP2_gradient_slave()
  use precision
  use tensor_type_def_module
  use dec_typedef_module
  use lsparameters
  use lsmpi_type
  use infpar_module
  
  ! DEC DEPENDENCIES (within deccc directory)   
  ! *****************************************!
  use rimp2_gradient_module
  use array2_simple_operations
  use decmpi_module
  implicit none
  type(decfrag) :: MyFragment
  !local variables
  integer :: natoms,nbasis,noccEOS,nvirtAOS
  type(array2) :: CvirtAOS, CoccEOS
  real(realk),pointer :: ThetaOcc(:),RIMP2grad(:,:)
  logical,pointer :: dopair_occ(:,:)
  integer(kind=long) :: nsize
  call mpi_communicate_MyFragment(MyFragment)
  noccEOS = MyFragment%noccEOS   ! Number of occupied functions in EOS 
  nvirtAOS = MyFragment%nvirtAOS ! Number of virtual functions in AOS 
  nbasis = MyFragment%nbasis     ! Number of basis functions in fragment
  natoms = MyFragment%natoms     ! Number of atoms in fragment

  nsize = nvirtAOS*(noccEOS*i8*noccEOS)*nvirtAOS
  call mem_alloc(ThetaOcc,nsize)
  !communicate ThetaOcc - assume this fits on all nodes 
  call ls_mpibcast(ThetaOcc,nsize,infpar%master,infpar%lg_comm)
  call mem_alloc(dopair_occ,noccEOS,noccEOS)
  !communicate dopair_occ - assume this fits on all nodes 
  call ls_mpibcast(dopair_occ,noccEOS,noccEOS,infpar%master,infpar%lg_comm)

  ! Get virtual MO coefficient matrix in array2 form
  CvirtAOS = array2_init([nbasis,nvirtAOS],MyFragment%Cv(1:nbasis,1:nvirtAOS))  
  ! Get occupied MO coefficient matrix for EOS orbitals in array2 form
  CoccEOS = array2_init_plain([MyFragment%nbasis,MyFragment%noccEOS])
  call extract_occupied_EOS_MO_indices(CoccEOS,MyFragment)

  call mem_alloc(RIMP2GRAD,3,natoms)

  call RIMP2_gradient_driver(MyFragment,ThetaOCC,RIMP2grad,natoms,&
       nbasis,noccEOS,nvirtAOS,CoccEOS%val,CvirtAOS%val,dopair_occ)

  call mem_dealloc(RIMP2GRAD)
  call mem_dealloc(ThetaOcc)
  call mem_dealloc(dopair_occ)
  call array2_free(CvirtAOS)
  call array2_free(CoccEOS) 

end subroutine RIMP2_gradient_slave
#endif

