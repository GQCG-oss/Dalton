!> @file

!> Module to handle DEC-MP2 gradient construction, see the mp2grad structure for details.
!> The theory is given in JCP 137, 114102 (2012). 
!> \author Kasper Kristensen

!> Module to handle DEC-MP2 gradient construction, see the mp2grad structure for details.
module mp2_gradient_module

  use fundamental
  use precision
  use lstiming!,only: lstimer
  use typedeftype!,only: lsitem
  use matrix_module!,only:matrix,matrixp
  use files!,only:lsopen,lsclose
  use matrix_util!,only:util_MO_to_AO_different_trans
  use matrix_operations!,only: mat_dotproduct, mat_init, mat_zero, &
!       & mat_free, mat_set_from_full, mat_trans, mat_daxpy, mat_mul,&
!       & mat_add, mat_assign
  use matrix_operations_aux
  use memory_handling!,only: mem_alloc,mem_dealloc, mem_turnonthread_memory,&
!       & init_threadmemvar, collect_thread_memory, mem_turnoffthread_memory
  use dec_typedef_module
  use IntegralInterfaceMOD!,only: ii_get_twoelectron_gradient, ii_get_reorthonormalization, &
!       & ii_get_oneelectron_gradient, ii_get_nn_gradient


  ! DEC DEPENDENCIES (within deccc directory)   
  ! *****************************************
  use dec_fragment_utils
  use crop_tools_module
  use array2_simple_operations!,only: array2_init, array2_transpose, array2_free,&
 !      & extract_occupied_EOS_MO_indices
  use array4_simple_operations!,only: array4_init, array4_contract1, array4_free,&
 !      & array4_reorder, array4_contract3
  use ccintegrals!,only:dec_fock_transformation
  use rimp2_gradient_module

  public :: init_fullmp2grad,free_fullmp2grad,single_calculate_mp2gradient_driver,&
       &pair_calculate_mp2gradient_driver,read_gradient_and_energies_for_restart, &
       &write_gradient_and_energies_for_restart,free_mp2grad,get_mp2gradient_main,&
       &dec_get_error_for_geoopt,update_full_mp2gradient,nullify_mp2dens,nullify_mp2grad,&
       &init_mp2grad
  private

contains


  !> \brief Initiate full mp2 gradient structure.
  !> \author Kasper Kristensen
  !> \date October 2011
  subroutine init_fullmp2grad(MyMolecule,fullgrad)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Full MP2 gradient structure
    type(fullmp2grad),intent(inout) :: fullgrad

    ! Dimensions
    ! **********

    ! Number of occupied orbitals
    fullgrad%nocc = MyMolecule%nocc
    ! Number of virtual orbitals
    fullgrad%nvirt = MyMolecule%nvirt
    ! Number of basis functions
    fullgrad%nbasis = MyMolecule%nbasis
    ! Number of atoms
    fullgrad%natoms = MyMolecule%natoms

    ! Zero energies
    fullgrad%Ecorr = 0E0_realk
    fullgrad%EHF = 0E0_realk
    fullgrad%Etot = 0E0_realk


    ! Initialize two-dimensional arrays and set elements to zero
    ! **********************************************************
    call mem_alloc(fullgrad%rho,fullgrad%nbasis,fullgrad%nbasis)
    fullgrad%rho = 0.0_realk
    call mem_alloc(fullgrad%Phi,fullgrad%nbasis,fullgrad%nbasis)
    fullgrad%Phi = 0.0_realk
    call mem_alloc(fullgrad%Ltheta,3,fullgrad%natoms)
    fullgrad%Ltheta=0.0_realk
    call mem_alloc(fullgrad%mp2gradient,3,fullgrad%natoms)
    fullgrad%mp2gradient=0.0_realk

  end subroutine init_fullmp2grad



  !> \brief Free full mp2 gradient structure.
  !> \author Kasper Kristensen
  !> \date October 2011
  subroutine free_fullmp2grad(fullgrad)

    implicit none
    !> Full MP2 gradient structure
    type(fullmp2grad),intent(inout) :: fullgrad

    call mem_dealloc(fullgrad%rho)
    call mem_dealloc(fullgrad%Phi)
    call mem_dealloc(fullgrad%Ltheta)
    call mem_dealloc(fullgrad%mp2gradient)

  end subroutine free_fullmp2grad


  !> \brief Nullify pointers in fragment gradient structure
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine nullify_mp2grad(grad)
    implicit none
    !> MP2 gradient for fragment
    type(mp2grad), intent(inout) :: grad

    nullify(grad%atoms_idx)
    nullify(grad%Phioo)
    nullify(grad%Phivv)
    nullify(grad%Ltheta)
    nullify(grad%PhiAO)

  end subroutine nullify_mp2grad


  !> \brief Nullify pointers in fragment density structure
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine nullify_mp2dens(dens)
    implicit none
    !> MP2 density for fragment
    type(mp2dens), intent(inout) :: dens

    nullify(dens%EOSatoms)
    nullify(dens%basis_idx)
    nullify(dens%Y)
    nullify(dens%X)
    nullify(dens%rho)
    nullify(dens%Phivo)
    nullify(dens%Phiov)

  end subroutine nullify_mp2dens



  !> \brief Init MP2 gradient structure for single fragment, assumes that the fragment structure (decfrag)
  !> has already been initiated!
  !> \author Kasper Kristensen
  !> \date October 2011
  subroutine init_mp2grad(fragment,grad)

    implicit none
    !> Fragment (single or pair) for which to initialize MP2 gradient info
    type(decfrag),intent(inout) :: fragment
    !> MP2 gradient for fragment
    type(mp2grad), intent(inout) :: grad
    integer :: i,atom1,atom2
    logical :: is_pair


    ! Check if it is a pair fragment
    ! ******************************
    atom1 = fragment%EOSatoms(1)
    if(fragment%nEOSatoms==2) then  ! pair fragment
       is_pair=.true.
       write(DECinfo%output,*) 'Initiating MP2 gradient structure for pair fragment ', &
            & fragment%EOSatoms
       atom2 = fragment%EOSatoms(2)
    else
       is_pair=.false.
       write(DECinfo%output,*) 'Initiating MP2 gradient structure for single fragment ', fragment%EOSatoms(1)
       atom2 = 0
    end if


    ! Initiate information used for fragment contribution to MP2 density
    ! ******************************************************************
    if(is_pair) then
       call pair_init_mp2dens(fragment,grad%dens,atom1,atom2)
    else
       call single_init_mp2dens(fragment,grad%dens)
    end if



    ! Remaining information, not used for MP2 density (see mp2grad structure)
    ! ***********************************************************************

    ! Number of atoms in fragment (atomic extent)
    grad%natoms = fragment%natoms

    ! Atomic indices for atoms in atomic extent
    if(associated(grad%atoms_idx)) then
       call mem_dealloc(grad%atoms_idx)
       nullify(grad%atoms_idx)
    end if
    call mem_alloc(grad%atoms_idx,grad%natoms)
    do i=1,grad%natoms
       grad%atoms_idx(i)=fragment%atoms_idx(i)
    end do

    ! Occ-occ block of Phi matrix - dimension: valence,(core+valence)  for frozen core
    if(associated(grad%Phioo)) then
       call mem_dealloc(grad%Phioo)
       nullify(grad%Phioo)
    end if
    call mem_alloc(grad%Phioo,grad%dens%nocc,grad%dens%nocctot)
    grad%Phioo=0E0_realk

    ! Virt-virt block of Phi matrix
    if(associated(grad%Phivv)) then
       call mem_dealloc(grad%Phivv)
       nullify(grad%Phivv)
    end if
    call mem_alloc(grad%Phivv,grad%dens%nvirt,grad%dens%nvirt)
    grad%Phivv=0E0_realk

    ! Total Phi matrix in AO basis
    if(associated(grad%PhiAO)) then
       call mem_dealloc(grad%PhiAO)
       nullify(grad%PhiAO)
    end if
    call mem_alloc(grad%PhiAO,grad%dens%nbasis,grad%dens%nbasis)
    grad%PhiAO=0E0_realk

    ! Lheta fragment contribution to molecular MP2 gradient
    if(associated(grad%Ltheta)) then
       call mem_dealloc(grad%Ltheta)
       nullify(grad%Ltheta)
    end if
    call mem_alloc(grad%Ltheta,3,grad%natoms)
    grad%Ltheta=0E0_realk


  end subroutine init_mp2grad



  !> \brief Free MP2 gradient structure
  !> \author Kasper Kristensen
  !> \date October 2011
  subroutine free_mp2grad(grad)

    implicit none
    type(mp2grad), intent(inout) :: grad

    ! Free MP2 density structure
    call free_mp2dens(grad%dens)
    grad%natoms=0

    ! Deallocate vectors/arrays
    if(associated(grad%atoms_idx)) call mem_dealloc(grad%atoms_idx)
    if(associated(grad%Phioo)) call mem_dealloc(grad%Phioo)
    if(associated(grad%Phivv)) call mem_dealloc(grad%Phivv)
    if(associated(grad%PhiAO)) call mem_dealloc(grad%PhiAO)
    if(associated(grad%Ltheta)) call mem_dealloc(grad%Ltheta)

  end subroutine free_mp2grad




  !> \brief Driver for calculating contributions to MP2 gradient for single fragments.
  !> \author Kasper Kristensen
  !> \date October 2011
  subroutine single_calculate_mp2gradient_driver(MyFragment,t2occ,t2virt,&
       & VOOO,VOVV,VOVOocc,VOVOvirt,grad)


    implicit none
    !> Atomic fragment
    type(decfrag),intent(inout) :: MyFragment
    !> t2 amplitudes t_{IJ}^{CD}, only for EOS orbitals using occupied partitioning, order:  (C,I,D,J)
    type(tensor),intent(inout) :: t2occ  ! ordered as (C,I,J,D) at output
    !> t2 amplitudes t_{KL}^{AB}, only for EOS orbitals using virtual partitioning, order: (A,K,B,L)
    type(tensor),intent(in) :: t2virt
    !> (C I | J L) integrals stored as (C,I,J,L)    [see index conventions in mp2.f90]
    type(tensor),intent(in) :: VOOO
    !> (B K | A C) integrals stored as (B,K,A,C)    [see index conventions in mp2.f90]
    type(tensor),intent(in) :: VOVV
    !> (C I | D J) integrals stored as (C,I,D,J)   [using occ partitioning]
    type(tensor),intent(inout) :: VOVOocc   ! ordered as (C,I,J,D) at output
    !> (A K | B L) integrals stored as (A,K,B,L)   [using virt partitioning]
    type(tensor),intent(in)    :: VOVOvirt
    type(array4) :: t2occ_arr4
    type(array4) :: t2virt_arr4
    type(array4) :: VOOO_arr4
    type(array4) :: VOVV_arr4
    type(array4) :: VOVOocc_arr4
    type(array4) :: VOVOvirt_arr4
    !> MP2 gradient structure
    type(mp2grad),intent(inout) :: grad
    type(array4) :: ThetaOCC, ThetaVIRT
    real(realk) :: tcpu,twall, tcpu1,tcpu2, twall1,twall2
    logical :: all_ok

    call LSTIMER('START',tcpu,twall,DECinfo%output)
    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    ! Sanity check
    all_ok=.false.
    if( t2occ%itype == TT_DENSE .and. t2virt%itype == TT_DENSE .and. VOOO%itype == TT_DENSE .and.&
         & VOVV%itype == TT_DENSE .and. VOVOocc%itype == TT_DENSE .and. VOVOvirt%itype == TT_DENSE ) then
       all_ok=.true.
    end if
    ! Unrelaxed density, we don't use VOOO and VOVV so don't check those
    if(DECinfo%unrelaxed) then
       if( t2occ%itype == TT_DENSE .and. t2virt%itype == TT_DENSE .and.&
            & VOVOocc%itype == TT_DENSE .and. VOVOvirt%itype == TT_DENSE ) then
          all_ok=.true.
       end if
    end if

    if(.not. all_ok) then
       call lsquit("ERROR(single_calculate_mp2gradient_driver) a PDM version needs to be implemented",-1)
    endif

    write(DECinfo%output,*) 'Calculating MP2 gradient for fragment', MyFragment%EOSatoms(1)

    ! Init MP2 gradient structure
    call init_mp2grad(MyFragment,grad)

    !use the trick for arrays to just associate them
    t2occ_arr4%val          => t2occ%elm4
    t2occ_arr4%dims         =  t2occ%dims
    t2occ_arr4%nelements    =  t2occ%nelms
    t2virt_arr4%val         => t2virt%elm4
    t2virt_arr4%dims        =  t2virt%dims
    t2virt_arr4%nelements   =  t2virt%nelms
    VOVOocc_arr4%val        => VOVOocc%elm4
    VOVOocc_arr4%dims       =  VOVOocc%dims
    VOVOocc_arr4%nelements  =  VOVOocc%nelms
    VOVOvirt_arr4%val       => VOVOvirt%elm4
    VOVOvirt_arr4%dims      =  VOVOvirt%dims
    VOVOvirt_arr4%nelements =  VOVOvirt%nelms
    if(.not. DECinfo%unrelaxed) then
       VOOO_arr4%val           => VOOO%elm4
       VOOO_arr4%dims          =  VOOO%dims
       VOOO_arr4%nelements     =  VOOO%nelms
       VOVV_arr4%val           => VOVV%elm4
       VOVV_arr4%dims          =  VOVV%dims
       VOVV_arr4%nelements     =  VOVV%nelms
    end if

    ! Get Theta arrays for occ and virt EOS
    ! *************************************
    call construct_theta_array(t2occ_arr4,ThetaOCC)
    call construct_theta_array(t2virt_arr4,ThetaVIRT)

    ! Reorder t2occ, ThetaOCC, and VOVOocc for easy contractions
    ! **********************************************************
    ! Theta(C,I,D,J) --> Theta(C,I,J,D)
    call array4_reorder(ThetaOCC,[1,2,4,3])
    ! t2(C,I,D,J) --> t2(C,I,J,D)
    call tensor_reorder(t2occ,[1,2,4,3])
    t2occ_arr4%dims = t2occ%dims
    ! VOVOocc(C,I,D,J) --> VOVOocc(C,I,J,D)
    call tensor_reorder(VOVOocc,[1,2,4,3])
    VOVOocc_arr4%dims =  VOVOocc%dims

    ! Calculate MP2 density contributions to gradient
    ! ***********************************************
    call single_calculate_mp2density(MyFragment,t2occ_arr4,t2virt_arr4,ThetaOCC,&
         &ThetaVIRT,VOOO_arr4,VOVV_arr4,grad%dens)

    ! Calculate remaining contributions to gradient (not used for MP2 density)
    ! ************************************************************************
    if(DECinfo%gradient) then
       call single_calculate_mp2gradient(MyFragment,ThetaOCC,ThetaVIRT,VOVOocc_arr4,VOVOvirt_arr4,grad)
    end if

    ! Free stuff
    call array4_free(ThetaOCC)
    call array4_free(ThetaVIRT)

    t2occ_arr4%val          => null()
    t2occ_arr4%dims         =  0
    t2occ_arr4%nelements    =  0
    t2virt_arr4%val         => null()
    t2virt_arr4%dims        =  0
    t2virt_arr4%nelements   =  0
    VOVOocc_arr4%val        => null()
    VOVOocc_arr4%dims       =  0
    VOVOocc_arr4%nelements  =  0
    VOVOvirt_arr4%val       => null()
    VOVOvirt_arr4%dims      =  0
    VOVOvirt_arr4%nelements =  0
    if(.not. DECinfo%unrelaxed) then
       VOVV_arr4%val           => null()
       VOVV_arr4%dims          =  0
       VOVV_arr4%nelements     =  0
       VOOO_arr4%val           => null()
       VOOO_arr4%dims          =  0
       VOOO_arr4%nelements     =  0
    end if

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)
    call LSTIMER('SINGLE MP2DENS',tcpu,twall,DECinfo%output)

  end subroutine single_calculate_mp2gradient_driver



  !> \brief Calculate contributions to MP2 gradient from single fragment.
  !> NOTE: Only the MP2 gradient data not contained in the mp2dens sub-structure are calculated here, i.e.,
  !> Phioo, Phivv, and Ltheta -- see mp2 grad structure.
  !> (The quantities in the mp2dens substructure are calculated using single_calculate_mp2density.)
  !> \author Kasper Kristensen
  !> \date October 2011
  subroutine single_calculate_mp2gradient(MyFragment,ThetaOCC,ThetaVIRT,VOVOocc,VOVOvirt,grad)

    implicit none

    !> Atomic fragment
    type(decfrag),intent(inout) :: MyFragment
    !> Theta array Theta_{IJ}^{CD}, only for EOS orbitals using occupied partitioning, order:  (C,I,J,D)
    type(array4),intent(in) :: ThetaOCC
    !> Theta array Theta_{KL}^{AB}, only for EOS orbitals using virtual partitioning, order:  (A,K,B,L)
    type(array4),intent(in) :: ThetaVIRT
    !> (C I | D J) integrals stored as (C,I,J,D)   [using occ partitioning]
    type(array4),intent(inout) :: VOVOocc
    !> (A K | B L) integrals stored as (A,K,B,L)   [using virt partitioning]
    type(array4),intent(in) :: VOVOvirt
    !> MP2 gradient contribution from given single fragment
    type(mp2grad),intent(inout) :: grad
    type(array2) :: Phioo, Phivv
    integer :: noccEOS,nvirtEOS,noccAOS,nvirtAOS,nocctot
    logical :: something_wrong
    real(realk) :: tcpu,twall


    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Init stuff
    noccEOS = Thetaocc%dims(2) ! number of occupied EOS orbitals
    noccAOS = Thetavirt%dims(2) ! number of occupied AOS orbitals
    nvirtEOS = Thetavirt%dims(1) ! number of virtual EOS orbitals
    nvirtAOS = Thetaocc%dims(1) ! number of virtual AOS orbitals
    ! number of occupied core+valence orbitals (only different from noccAOS for frozen approx)
    nocctot = VOVOvirt%dims(4)
    ! Just in case, zero matrices in grad structure
    grad%Phioo = 0e0_realk
    grad%Phivv = 0e0_realk
    grad%Ltheta = 0e0_realk

    ! Sanity check 1: Gradient input structure is correct
    something_wrong=.false.
    if(noccAOS/=grad%dens%nocc) something_wrong=.true.
    if(nvirtAOS/=grad%dens%nvirt) something_wrong=.true.
    if(something_wrong) then
       write(DECinfo%output,*) 'Grad structure: nocc =   ', grad%dens%nocc
       write(DECinfo%output,*) 'Grad structure: nvirt =  ', grad%dens%nvirt
       write(DECinfo%output,*) 'Theta, occ AOS dimension ', noccAOS
       write(DECinfo%output,*) 'Theta, virt AOS dimension', nvirtAOS
       call lsquit('single_calculate_mp2gradient: &
            & Something wrong with input dimensions!',DECinfo%output)
    end if

    ! Sanity check 2: Array dimensions match for last two indices of Theta and VOVO
    if(noccEOS /= Thetaocc%dims(3)) something_wrong=.true.
    if(nvirtAOS /= Thetaocc%dims(4)) something_wrong=.true.
    if(nvirtEOS /= Thetavirt%dims(3)) something_wrong=.true.
    if(noccAOS /= Thetavirt%dims(4)) something_wrong=.true.
    if(noccEOS /= VOVOocc%dims(3)) something_wrong=.true.
    if(nvirtAOS /= VOVOocc%dims(4)) something_wrong=.true.
    if(nvirtEOS /= VOVOvirt%dims(3)) something_wrong=.true.
    if(.not. DECinfo%frozencore) then
       if(noccAOS /= VOVOvirt%dims(4)) something_wrong=.true.
    end if
    if(something_wrong) then
       write(DECinfo%output,*) 'EOS: nocc, nvirt', noccEOS, nvirtEOS
       write(DECinfo%output,*) 'AOS: nocc, nvirt', noccAOS, nvirtAOS
       write(DECinfo%output,*) 'VOVOocc dims    ', VOVOocc%dims
       write(DECinfo%output,*) 'Thetaocc dims   ', Thetaocc%dims
       write(DECinfo%output,*) 'VOVOvirt dims   ', VOVOvirt%dims
       write(DECinfo%output,*) 'Thetavirt dims  ', Thetavirt%dims
       call lsquit('single_calculate_mp2gradient: &
            & Something wrong with Theta and/or integral dimensions!',DECinfo%output)
    end if



    ! ****************************************************************************************
    !                  Calculate contribution to virt-virt block of Phi matrix               *
    ! ****************************************************************************************

    ! Phivv(D,E) = sum_{CIJ} Theta(C,I,J,D) g(C,I,J,E)
    Phivv = array2_init([nvirtAOS,nvirtAOS])

    call array4_contract3(ThetaOCC,VOVOocc,Phivv)
    grad%Phivv(1:nvirtAOS,1:nvirtAOS) = Phivv%val(1:nvirtAOS,1:nvirtAOS)

    call array2_free(Phivv)
    call LSTIMER('PHIVV MATRIX',tcpu,twall,DECinfo%output)



    ! ****************************************************************************************
    !                  Calculate contribution to occ-occ block of Phi matrix                 *
    ! ****************************************************************************************

    ! Phioo(L,M) = sum_{ABK} Theta(B,K,A,L) g(B,K,A,M)
    Phioo = array2_init([noccAOS,nocctot])
    call array4_contract3(ThetaVIRT,VOVOvirt,Phioo)
    grad%Phioo(1:noccAOS,1:nocctot) = Phioo%val(1:noccAOS,1:nocctot)
    call array2_free(Phioo)

    call LSTIMER('PHIOO MATRIX',tcpu,twall,DECinfo%output)


    ! ****************************************************************************************
    !                          Calculate Phi matrix in AO basis                              *
    ! ****************************************************************************************
    call fragment_Phi_matrix_in_AO_basis_wrapper(MyFragment,grad)
 

    ! ****************************************************************************************
    !                  Calculate contribution to Ltheta conponent of MP2 grad                *
    ! ****************************************************************************************
    call single_fragment_Ltheta_contribution(ThetaOCC,MyFragment,grad)


  end subroutine single_calculate_mp2gradient


  !> \brief Get Ltheta contribution to MP2 gradient for single fragment:
  !> Ltheta = 1/2 * sum_{ijab} Theta_{ij}^{ab} gx_{ij}^{ab}   (i and j belong to central atom in MyFragment)
  !> Theta_{ij}^{ab} = 8*t2_{ij}^{ab} - 4*t2_{ij}^{ba}        (see construct_theta_array)
  !> gx: Two-electron integrals differentiated w.r.t. nuclear coordinates.
  !> \author Kasper Kristensen
  !> \date October 2011
  subroutine single_fragment_Ltheta_contribution(ThetaOCC,MyFragment,grad)

    implicit none
    !> Theta array Theta_{IJ}^{CD}, only for EOS orbitals using occupied partitioning, order:  (C,I,J,D)
    type(array4),intent(in) :: ThetaOCC
    !> Atomic fragment
    type(decfrag),intent(inout) :: MyFragment
    !> Structure containing MP2 gradient info, including Ltheta.
    type(mp2grad), intent(inout) :: grad
    type(array4) :: ThetaAO, dens4index
    type(array2) :: CvirtAOS, CoccEOS
    integer :: noccEOS, nvirtAOS, nbasis, matdim,i,j,natoms
    real(realk) :: tcpu,twall,tKcpu1,tKwall1
    logical :: something_wrong
    logical,pointer :: dopair_occ(:,:)

    ! If only MP2 density is requested, skip expensive Ltheta calculation,
    ! set Ltheta (which is not used) to zero for consistency.
    if(.not. DECinfo%gradient) then
       Grad%Ltheta = 0.0_realk
       return
    end if

    call LSTIMER('START',tcpu,twall,DECinfo%output)


    ! Initialize stuff
    ! ****************
    ! Number of occupied EOS orbitals in fragment
    noccEOS = ThetaOCC%dims(2)
    ! Number of virtual AOS orbitals in fragment
    nvirtAOS = ThetaOCC%dims(1)
    ! Number of basis functions in fragment
    nbasis = MyFragment%nbasis
    ! Number of atoms in fragment
    natoms = MyFragment%natoms
    ! Number of AO matrices to use as input in exchange gradient routine
    matdim = noccEOS*noccEOS

    ! Sanity check 1: Consistency of fragment and LSitem
    ! **************************************************
    if(natoms /= MyFragment%MyLsitem%setting%molecule(1)%p%natoms) then
       write(DECinfo%output,*) 'Number of atoms in fragment = ', natoms
       write(DECinfo%output,*) 'Number of atoms in LsItem   = ', MyFragment%MyLsitem%setting%molecule(1)%p%natoms
       call lsquit('single_fragment_Ltheta_contribution:&
            & The number of atoms in the fragment does not match the&
            & number of atoms in LsItem!',-1)
    end if

    ! Sanity check 2: Dimensions in fragment and gradient structures
    ! **************************************************************
    something_wrong=.false.
    if(grad%dens%nocc /= MyFragment%noccAOS) something_wrong=.true.
    if(grad%dens%nvirt /= MyFragment%nvirtAOS) something_wrong=.true.
    if(grad%natoms /= MyFragment%natoms) something_wrong=.true.
    if(something_wrong) then
       write(DECinfo%output,*) 'Gradient: nocc   = ', grad%dens%nocc
       write(DECinfo%output,*) 'Gradient: nvirt  = ', grad%dens%nvirt
       write(DECinfo%output,*) 'Gradient: natoms = ', grad%natoms
       write(DECinfo%output,*) 'Fragment: nocc   = ', MyFragment%noccAOS
       write(DECinfo%output,*) 'Fragment: nvirt  = ', MyFragment%nvirtAOS
       write(DECinfo%output,*) 'Fragment: natoms = ', MyFragment%natoms
       call lsquit('Something wrong in single_fragment_Ltheta_contribution:&
            & Dimensions in fragment does not match dimensions in gradient structure', DECinfo%output)
    end if



    ! Prepare arrays for constructing Ltheta gradient
    ! ***********************************************

    ! Get virtual MO coefficient matrix in array2 form
    CvirtAOS = array2_init([nbasis,nvirtAOS],&
         &MyFragment%Cv(1:nbasis,1:nvirtAOS))

    ! Get occupied MO coefficient matrix for EOS orbitals in array2 form
    CoccEOS = array2_init_plain([MyFragment%nbasis,MyFragment%noccEOS])
    call extract_occupied_EOS_MO_indices(CoccEOS,MyFragment)
    
    ! Calculate Ltheta contribution
    ! *****************************

    call LSTIMER('START',tKcpu1,tKwall1,DECinfo%output)

    IF(DECinfo%ccModel .EQ. MODEL_RIMP2)THEN
       call mem_alloc(dopair_occ,MyFragment%noccEOS,MyFragment%noccEOS)
       dopair_occ = .TRUE.
       call RIMP2_gradient_driver(MyFragment,ThetaOCC%val,Grad%Ltheta,natoms,&
            & nbasis,MyFragment%noccEOS,nvirtAOS,CoccEOS%val,CvirtAOS%val,&
            & dopair_occ)
       call mem_dealloc(dopair_occ)
       call array2_free(CvirtAOS)
       call array2_free(CoccEOS)
       ! Timing
       call LSTIMER('LTHETA RI GRAD',tKcpu1,tKwall1,DECinfo%output)
    ELSE
       ! Transform virtual Theta indices to AO indices
       call transform_virtual_Theta_indices_to_AO(ThetaOCC, CvirtAOS, ThetaAO)       
       ! Get 4-index density: dens4index(mu,nu,i,j) = Cocc(mu,i)*Cocc(nu,j)
       ! for occupied EOS indices i and j.
       call get_4_dimensional_HFdensity(CoccEOS,dens4index)      
       ! Done with MO coefficient matrices
       call array2_free(CvirtAOS)
       call array2_free(CoccEOS)
       call II_get_K_gradientfull(Grad%Ltheta(1:3,1:natoms),&
            & dens4index%val,ThetaAO%val,nbasis,matdim,matdim,&
            & MyFragment%MyLsitem%setting,DECinfo%output,DECinfo%output)
       ! Timing
       call LSTIMER('LTHETA K GRAD',tKcpu1,tKwall1,DECinfo%output)
       ! By definition:
       !
       ! Ltheta = 1/2 * sum_{ijab} Theta_{ij}^{ab} gx_{ij}^{ab}
       !
       ! Due to conventions in II_get_K_gradientfull we effectively get
       ! -1/4 sum_{ijab} Theta_{ij}^{ab} gx_{ij}^{ab}
       ! using the procedure above.
       ! Therefore we multiply the final result by (-2)
       Grad%Ltheta(1:3,1:natoms) = -2E0_realk*Grad%Ltheta(1:3,1:natoms)
       ! Free stuff
       call array4_free(ThetaAO)
       call array4_free(dens4index)
    ENDIF

    call LSTIMER('LTHETA TOTAL',tcpu,twall,DECinfo%output)


  end subroutine single_fragment_Ltheta_contribution



  !> \brief Driver for calculating contributions to MP2 gradient for pair fragments.
  !> \author Kasper Kristensen
  !> \date October 2011
  subroutine pair_calculate_mp2gradient_driver(Fragment1,Fragment2,PairFragment,&
       & t2occ,t2virt,VOOO,VOVV,VOVOocc,VOVOvirt,grad)


    implicit none
    !> Fragment 1 in the pair fragment
    type(decfrag),intent(in) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(decfrag),intent(in) :: Fragment2
    !> Pair fragment
    type(decfrag),intent(inout) :: pairfragment
    !> t2 amplitudes t_{IJ}^{CD}, only for EOS orbitals using occupied partitioning, order:  (C,I,D,J)
    type(tensor),intent(in) :: t2occ
    !> t2 amplitudes t_{KL}^{AB}, only for EOS orbitals using virtual partitioning, order: (A,K,B,L)
    type(tensor),intent(in) :: t2virt
    !> (C I | J L) integrals stored as (C,I,J,L)    [see index conventions in mp2.f90]
    type(tensor),intent(in) :: VOOO
    !> (B K | A C) integrals stored as (B,K,A,C)    [see index conventions in mp2.f90]
    type(tensor),intent(in) :: VOVV
    !> (C I | D J) integrals stored as (C,I,D,J)   [using occ partitioning]
    type(tensor),intent(inout) :: VOVOocc
    !> (A K | B L) integrals stored as (A,K,B,L)   [using virt partitioning]
    type(tensor),intent(in) :: VOVOvirt
    type(array4) :: t2occ_arr4
    type(array4) :: t2virt_arr4
    type(array4) :: VOOO_arr4
    type(array4) :: VOVV_arr4
    type(array4) :: VOVOocc_arr4
    type(array4) :: VOVOvirt_arr4
    type(array4) :: ThetaOCC, ThetaVIRT
    !> MP2 gradient structure for pair
    type(mp2grad),intent(inout) :: grad
    real(realk) :: tcpu,twall, tcpu1,tcpu2, twall1,twall2
    logical :: all_ok

    call LSTIMER('START',tcpu,twall,DECinfo%output)
    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    ! Sanity check
    all_ok=.false.
    if( t2occ%itype == TT_DENSE .and. t2virt%itype == TT_DENSE .and. VOOO%itype == TT_DENSE .and.&
         & VOVV%itype == TT_DENSE .and. VOVOocc%itype == TT_DENSE .and. VOVOvirt%itype == TT_DENSE ) then
       all_ok=.true.
    end if
    ! Unrelaxed density, we don't use VOOO and VOVV so don't check those
    if(DECinfo%unrelaxed) then
       if( t2occ%itype == TT_DENSE .and. t2virt%itype == TT_DENSE .and.&
            & VOVOocc%itype == TT_DENSE .and. VOVOvirt%itype == TT_DENSE ) then
          all_ok=.true.
       end if
    end if

    if(.not. all_ok) then
       call lsquit("ERROR(pair_calculate_mp2gradient_driver) a PDM version needs to be implemented",-1)
    endif


    call LSTIMER('START',tcpu,twall,DECinfo%output)
    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    write(DECinfo%output,*) 'Calculating MP2 gradient for pair fragment', &
         & PairFragment%EOSatoms

    ! Init MP2 gradient structure
    call init_mp2grad(PairFragment,grad)

    !use the trick for arrays to just associate them
    t2occ_arr4%val          => t2occ%elm4
    t2occ_arr4%dims         =  t2occ%dims
    t2occ_arr4%nelements    =  t2occ%nelms
    t2virt_arr4%val         => t2virt%elm4
    t2virt_arr4%dims        =  t2virt%dims
    t2virt_arr4%nelements   =  t2virt%nelms
    VOVOocc_arr4%val        => VOVOocc%elm4
    VOVOocc_arr4%dims       =  VOVOocc%dims
    VOVOocc_arr4%nelements  =  VOVOocc%nelms
    VOVOvirt_arr4%val       => VOVOvirt%elm4
    VOVOvirt_arr4%dims      =  VOVOvirt%dims
    VOVOvirt_arr4%nelements =  VOVOvirt%nelms
    if(.not. DECinfo%unrelaxed) then
       VOOO_arr4%val           => VOOO%elm4
       VOOO_arr4%dims          =  VOOO%dims
       VOOO_arr4%nelements     =  VOOO%nelms
       VOVV_arr4%val           => VOVV%elm4
       VOVV_arr4%dims          =  VOVV%dims
       VOVV_arr4%nelements     =  VOVV%nelms
    end if

    ! Get Theta arrays for occ and virt EOS
    ! *************************************
    call construct_theta_array(t2occ_arr4,ThetaOCC)

    call construct_theta_array(t2virt_arr4,ThetaVIRT)
    ! NOTE: It is not useful to reorder these arrays here as in done in single_calculate_mp2gradient_driver.

    ! Calculate MP2 density contributions to gradient
    ! ***********************************************
    call pair_calculate_mp2density(fragment1,fragment2,PairFragment,t2occ_arr4,t2virt_arr4,&
         & ThetaOCC, ThetaVIRT, VOOO_arr4,VOVV_arr4,grad%dens)

    ! Calculate remaining contributions to gradient (not used for MP2 density)
    ! ************************************************************************
    if(DECinfo%gradient) then
       call pair_calculate_mp2gradient(fragment1,fragment2,pairfragment,&
            & ThetaOCC,ThetaVIRT,VOVOocc_arr4,VOVOvirt_arr4,grad)
    end if

    ! Free stuff
    call array4_free(ThetaOCC)
    call array4_free(ThetaVIRT)

    t2occ_arr4%val          => null()
    t2occ_arr4%dims         =  0
    t2occ_arr4%nelements    =  0
    t2virt_arr4%val         => null()
    t2virt_arr4%dims        =  0
    t2virt_arr4%nelements   =  0
    VOVOocc_arr4%val        => null()
    VOVOocc_arr4%dims       =  0
    VOVOocc_arr4%nelements  =  0
    VOVOvirt_arr4%val       => null()
    VOVOvirt_arr4%dims      =  0
    VOVOvirt_arr4%nelements =  0
    if(.not. DECinfo%unrelaxed) then
       VOVV_arr4%val           => null()
       VOVV_arr4%dims          =  0
       VOVV_arr4%nelements     =  0
       VOOO_arr4%val           => null()
       VOOO_arr4%dims          =  0
       VOOO_arr4%nelements     =  0
    end if


    call LSTIMER('START',tcpu2,twall2,DECinfo%output)
    call LSTIMER('PAIR MP2DENS',tcpu,twall,DECinfo%output)

  end subroutine pair_calculate_mp2gradient_driver



  !> \brief Calculate contributions to MP2 gradient from single fragment.
  !> NOTE: Only the MP2 gradient data not contained in the mp2dens sub-structure are calculated here, i.e.,
  !> Phioo, Phivv, and Ltheta -- see mp2 grad structure.
  !> (The quantities in the mp2dens substructure are calculated using single_calculate_mp2density.)
  !> \author Kasper Kristensen
  !> \date October 2011
  subroutine pair_calculate_mp2gradient(fragment1,fragment2,pairfragment,&
       & ThetaOCC,ThetaVIRT,VOVOocc,VOVOvirt,grad)

    implicit none

    !> Fragment 1 in the pair fragment
    type(decfrag),intent(in) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(decfrag),intent(in) :: Fragment2
    !> Pair fragment
    type(decfrag),intent(inout) :: pairfragment
    !> Theta array, only for EOS orbitals using occupied partitioning, order:  (C,I,D,J)
    type(array4),intent(inout) :: ThetaOCC       ! ordered as (C,I,J,D) at output
    !> Theta array, only for EOS orbitals using virtual partitioning, order:  (A,K,B,L)
    type(array4),intent(in) :: ThetaVIRT
    !> (C I | D J) integrals stored as (C,I,D,J)   [using occ partitioning]
    type(array4),intent(inout) :: VOVOocc
    !> (A K | B L) integrals stored as (A,K,B,L)   [using virt partitioning]
    type(array4),intent(in) :: VOVOvirt
    !> MP2 gradient contribution from given single fragment
    type(mp2grad),intent(inout) :: grad
    real(realk),pointer :: Phivv_tmp(:,:), Phioo_tmp(:,:)
    integer :: noccEOS,nvirtEOS,noccAOS,nvirtAOS,nocctot
    logical :: something_wrong
    real(realk) :: tcpu,twall
    integer :: i,j,k,l,m,a,b,c,d,e
    logical, pointer :: dopair_occ(:,:),dopair_virt(:,:)


    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Init stuff
    noccEOS = Thetaocc%dims(2) ! number of occupied EOS orbitals
    noccAOS = Thetavirt%dims(2) ! number of occupied AOS orbitals
    nvirtEOS = Thetavirt%dims(1) ! number of virtual EOS orbitals
    nvirtAOS = Thetaocc%dims(1) ! number of virtual AOS orbitals
    ! number of occupied core+valence orbitals (only different from noccAOS for frozen core approx)
    nocctot = VOVOvirt%dims(4)
    ! Just in case, zero matrices in grad structure
    grad%Phioo = 0e0_realk
    grad%Phivv = 0e0_realk
    grad%Ltheta = 0e0_realk

    ! Which "interaction pairs" to include for occ and virt space (avoid double counting)
    call mem_alloc(dopair_occ,noccEOS,noccEOS)
    call mem_alloc(dopair_virt,nvirtEOS,nvirtEOS)
    call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)
    call which_pairs_virt(Fragment1,Fragment2,PairFragment,dopair_virt)


    ! Sanity check 1: Gradient input structure is correct
    something_wrong=.false.
    if(noccAOS/=grad%dens%nocc) something_wrong=.true.
    if(nvirtAOS/=grad%dens%nvirt) something_wrong=.true.
    if(something_wrong) then
       write(DECinfo%output,*) 'Grad structure: nocc =   ', grad%dens%nocc
       write(DECinfo%output,*) 'Grad structure: nvirt =  ', grad%dens%nvirt
       write(DECinfo%output,*) 'Theta, occ AOS dimension ', noccAOS
       write(DECinfo%output,*) 'Theta, virt AOS dimension', nvirtAOS
       call lsquit('pair_calculate_mp2density: &
            & Something wrong with input dimensions!',DECinfo%output)
    end if

    ! Sanity check 2: Array dimensions match for last two indices of Theta and VOVO
    if(nvirtAOS /= Thetaocc%dims(3)) something_wrong=.true.
    if(noccEOS /= Thetaocc%dims(4)) something_wrong=.true.
    if(nvirtEOS /= Thetavirt%dims(3)) something_wrong=.true.
    if(noccAOS /= Thetavirt%dims(4)) something_wrong=.true.
    if(nvirtAOS /= VOVOocc%dims(3)) something_wrong=.true.
    if(noccEOS /= VOVOocc%dims(4)) something_wrong=.true.
    if(nvirtEOS /= VOVOvirt%dims(3)) something_wrong=.true.
    if(.not. DECinfo%frozencore) then
       if(noccAOS /= VOVOvirt%dims(4)) something_wrong=.true.
    end if
    if(something_wrong) then
       write(DECinfo%output,*) 'EOS: nocc, nvirt', noccEOS, nvirtEOS
       write(DECinfo%output,*) 'AOS: nocc, nvirt', noccAOS, nvirtAOS
       write(DECinfo%output,*) 'VOVOocc dims    ', VOVOocc%dims
       write(DECinfo%output,*) 'Thetaocc dims   ', Thetaocc%dims
       write(DECinfo%output,*) 'VOVOvirt dims   ', VOVOvirt%dims
       write(DECinfo%output,*) 'Thetavirt dims  ', Thetavirt%dims
       call lsquit('pair_calculate_mp2gradient: &
            & Something wrong with Theta and/or integral dimensions!',DECinfo%output)
    end if




    ! ****************************************************************************************
    !                  Calculate contribution to virt-virt block of Phi matrix               *
    ! ****************************************************************************************


    ! Use occupied partitioning scheme
    ! ! Phivv(D,E) = sum_{CIJ} Theta(C,I,D,J) g(C,I,E,J)
    ! -- where we only include the contributions if I and J belong to different atoms!
    call LSTIMER('START',tcpu,twall,DECinfo%output)

    call mem_TurnONThread_Memory()
    !$OMP PARALLEL DEFAULT(shared) PRIVATE(Phivv_tmp,i,j,c,d,e)
    call init_threadmemvar()
    call mem_alloc(Phivv_tmp,nvirtAOS,nvirtAOS)
    Phivv_tmp = 0E0_realk

    !$OMP DO SCHEDULE(dynamic,1)

    do i=1,noccEOS
       do j=1,noccEOS

          ! Only update for "interaction orbital pairs" - see which_pairs_occ
          if(dopair_occ(i,j)) then !PhivvDoPair

             do e=1,nvirtAOS
                do d=1,nvirtAOS
                   do c=1,nvirtAOS
                      Phivv_tmp(d,e) = Phivv_tmp(d,e) + ThetaOCC%val(c,i,d,j)*VOVOocc%val(c,i,e,j)
                   end do
                end do
             end do

          end if

       end do
    end do
    !$OMP END DO NOWAIT

    ! Total Phivv matrix is found by summing all thread contributions
    !$OMP CRITICAL
    grad%Phivv = grad%Phivv + Phivv_tmp
    !$OMP END CRITICAL

    call mem_dealloc(Phivv_tmp)
    call collect_thread_memory()
    !$OMP END PARALLEL
    call mem_TurnOffThread_Memory()
    call LSTIMER('PHIVV MATRIX',tcpu,twall,DECinfo%output)




    ! ****************************************************************************************
    !                  Calculate contribution to occ-occ block of Phi matrix                 *
    ! ****************************************************************************************

    ! Phioo(L,M) = sum_{ABK} Theta(A,K,B,L) g(A,K,B,M)
    ! -- where we only include the contributions if A and B belong to different atoms!

    Call mem_TurnONThread_Memory()
    !$OMP PARALLEL DEFAULT(shared) PRIVATE(Phioo_tmp,a,b,k,l,m)
    call init_threadmemvar()
    call mem_alloc(Phioo_tmp,noccAOS,nocctot)
    Phioo_tmp = 0E0_realk

    !$OMP DO SCHEDULE(dynamic,1)


    do a=1,nvirtEOS
       do b=1,nvirtEOS

          ! Only update for "interaction orbital pairs" - see which_pairs_virt
          if(dopair_virt(a,b)) then !PhiooDoPair

             do l=1,noccAOS
                do m=1,nocctot
                   do k=1,noccAOS
                      Phioo_tmp(l,m) = Phioo_tmp(l,m) + ThetaVIRT%val(a,k,b,l)*VOVOvirt%val(a,k,b,m)
                   end do
                end do
             end do

          end if

       end do
    end do
    !$OMP END DO NOWAIT

    ! Total Phioo matrix is found by summing all thread contributions
    !$OMP CRITICAL
    grad%Phioo = grad%Phioo + Phioo_tmp
    !$OMP END CRITICAL

    call mem_dealloc(Phioo_tmp)
    call collect_thread_memory()
    !$OMP END PARALLEL
    call mem_TurnOffThread_Memory()

    call LSTIMER('PHIOO MATRIX',tcpu,twall,DECinfo%output)



    ! ****************************************************************************************
    !                          Calculate Phi matrix in AO basis                              *
    ! ****************************************************************************************
    call fragment_Phi_matrix_in_AO_basis_wrapper(PairFragment,grad)


    ! ****************************************************************************************
    !                  Calculate contribution to Ltheta conponent of MP2 grad                *
    ! ****************************************************************************************

    ! We need to be consistent with Theta ordering assumed in pair_fragment_Ltheta_contribution.
    ! Reorder: Theta(C,I,D,J) --> Theta(C,I,J,D)
    call array4_reorder(ThetaOCC,[1,2,4,3])
    call pair_fragment_Ltheta_contribution(ThetaOCC,PairFragment,noccEOS,dopair_occ,grad)
    call mem_dealloc(dopair_occ)
    call mem_dealloc(dopair_virt)

  end subroutine pair_calculate_mp2gradient



  !> \brief Get Ltheta contribution to MP2 gradient for pair fragment:
  !> Ltheta = 1/2 * sum_{ijab} Theta_{ij}^{ab} gx_{ij}^{ab}   (i,j belong to different atoms in pair fragment)
  !> Theta_{ij}^{ab} = 8*t2_{ij}^{ab} - 4*t2_{ij}^{ba}        (see construct_theta_array)
  !> gx: Two-electron integrals differentiated w.r.t. nuclear coordinates.
  !> \author Kasper Kristensen
  !> \date October 2011
  subroutine pair_fragment_Ltheta_contribution(ThetaOCC,PairFragment,noccEOS,dopair_occ,grad)

    implicit none
    !> Theta array Theta_{IJ}^{CD}, only for EOS orbitals using occupied partitioning, order:  (C,I,J,D)
    type(array4),intent(in) :: ThetaOCC
    !> Pair fragment
    type(decfrag),intent(inout) :: PairFragment
    !> Number of occupied EOS orbitals
    integer,intent(in) :: noccEOS
    !> Which "interaction orbital pairs" to include (see which_pairs_occ)
    logical, dimension(noccEOS,noccEOS),intent(in) :: dopair_occ
    !> Structure containing MP2 gradient info, including Ltheta.
    type(mp2grad), intent(inout) :: grad
    type(array4) :: ThetaAO, dens4index
    type(array2) :: CvirtAOS, CoccEOS
    integer ::nvirtAOS, nbasis, matdim,ij,i,j,natoms_frag
    real(realk),pointer :: dens4index_pair(:,:,:), ThetaAO_pair(:,:,:)
    real(realk) :: tcpu,twall,tKcpu1,tKwall1
    logical :: something_wrong


    ! If only MP2 density is requested, skip expensive Ltheta calculation,
    ! set Ltheta (which is not used) to zero for consistency.
    if(.not. DECinfo%gradient) then
       Grad%Ltheta = 0.0_realk
       return
    end if

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Initialize stuff
    ! ****************
    ! Number of virtual AOS orbitals in fragment
    nvirtAOS = ThetaOCC%dims(1)
    ! Number of basis functions in fragment
    nbasis = PairFragment%nbasis
    ! Number of atoms in fragment
    natoms_frag = PairFragment%natoms
    ! Number of occupied pairs to consider
    matdim = count(dopair_occ)


    ! Sanity check 1: Consistency of fragment and LSitem
    ! **************************************************
    if(natoms_frag /= PairFragment%MyLsitem%setting%molecule(1)%p%natoms) then
       write(DECinfo%output,*) 'Number of atoms in fragment = ', natoms_frag
       write(DECinfo%output,*) 'Number of atoms in LsItem   = ', &
            & PairFragment%MyLsitem%setting%molecule(1)%p%natoms
       call lsquit('pair_fragment_Ltheta_contribution:&
            & The number of atoms in the fragment does not match the&
            & number of atoms in LsItem!',-1)
    end if

    ! Sanity check 2: Dimensions in fragment and gradient structures
    ! **************************************************************
    something_wrong=.false.
    if(grad%dens%nocc /= PairFragment%noccAOS) something_wrong=.true.
    if(grad%dens%nvirt /= PairFragment%nvirtAOS) something_wrong=.true.
    if(grad%natoms /= PairFragment%natoms) something_wrong=.true.
    if(something_wrong) then
       write(DECinfo%output,*) 'Gradient: nocc   = ', grad%dens%nocc
       write(DECinfo%output,*) 'Gradient: nvirt  = ', grad%dens%nvirt
       write(DECinfo%output,*) 'Gradient: natoms = ', grad%natoms
       write(DECinfo%output,*) 'Fragment: nocc   = ', PairFragment%noccAOS
       write(DECinfo%output,*) 'Fragment: nvirt  = ', PairFragment%nvirtAOS
       write(DECinfo%output,*) 'Fragment: natoms = ', PairFragment%natoms
       call lsquit('Something wrong in pair_fragment_Ltheta_contribution:&
            & Dimensions in fragment does not match dimensions in gradient structure', DECinfo%output)
    end if


    ! Prepare arrays for constructing Ltheta gradient
    ! *************************************************

    ! Get virtual MO coefficient matrix in array2 form
    CvirtAOS = array2_init([nbasis,nvirtAOS],&
         &PairFragment%Cv(1:nbasis,1:nvirtAOS))

    ! Get occupied MO coefficient matrix for EOS orbitals in array2 form
    CoccEOS = array2_init_plain([PairFragment%nbasis,PairFragment%noccEOS])
    call extract_occupied_EOS_MO_indices(CoccEOS,PairFragment)

    IF(DECinfo%ccModel .EQ. MODEL_RIMP2)THEN
       call RIMP2_gradient_driver(PairFragment,ThetaOCC%val,Grad%Ltheta,natoms_frag,&
            & nbasis,PairFragment%noccEOS,nvirtAOS,CoccEOS%val,CvirtAOS%val,&
            & dopair_occ)
       call array2_free(CvirtAOS)
       call array2_free(CoccEOS)
       ! Timing
       call LSTIMER('LTHETA RI GRAD',tKcpu1,tKwall1,DECinfo%output)
    ELSE
       ! Transform virtual Theta indices to AO indices
       ! (required for calling gradient exchange routine).
       call transform_virtual_Theta_indices_to_AO(ThetaOCC, CvirtAOS, ThetaAO)
       
       ! Get 4-index density: dens4index(mu,nu,i,j) = Cocc(mu,i)*Cocc(nu,j)
       ! for occupied EOS indices i and j.
       call get_4_dimensional_HFdensity(CoccEOS,dens4index)
       
       ! Done with MO coefficient matrices
       call array2_free(CvirtAOS)
       call array2_free(CoccEOS)
       
       ! Extract occupied pair EOS indices where "i" and "j" are assigned to different atoms
       ! ***********************************************************************************

       ! ThetaAO and 4-index density arrays contain all (i,j) combinations where
       ! i and j are assigned to atom P OR atom Q. We now extract the ones where
       !
       ! (i assigned to P,   j assigned to Q)           or
       ! (j assigned to P,   i assigned to Q)          
       !
       ! The are matdim possible (i,j) combinations satisfying this.
       call mem_alloc(dens4index_pair,nbasis,nbasis,matdim)
       call mem_alloc(thetaAO_pair,nbasis,nbasis,matdim)
       
       ij = 0
       do j=1,noccEOS
          do i=1,noccEOS
             
             if(dopair_occ(i,j)) then ! Extract orbital pair information
                ij = ij+1
                
                ! Theta array
                thetaAO_pair(:,:,ij) = ThetaAO%val(:,:,i,j)
                
                ! Density array
                dens4index_pair(:,:,ij) = dens4index%val(:,:,i,j)
                
             end if
             
          enddo
       enddo
       
       ! Sanity check
       if(ij /= matdim) then
          write(DECinfo%output,*) 'ij     = ', ij
          write(DECinfo%output,*) 'matdim = ', matdim
          call lsquit('pair_fragment_Ltheta_contribution: ij counter is different from matdim', DECinfo%output)
       end if
       
       ! Done with 4 dimensional arrays in type array4 form
       call array4_free(ThetaAO)
       call array4_free(dens4index)
       
       
       
       ! Calculate Ltheta contribution
       ! *****************************
       call LSTIMER('START',tKcpu1,tKwall1,DECinfo%output)
       call II_get_K_gradientfull(Grad%Ltheta,&
            & dens4index_pair,ThetaAO_pair,nbasis,matdim,matdim,&
            & PairFragment%MyLsitem%setting,DECinfo%output,DECinfo%output)
       call LSTIMER('LTHETA K GRAD',tKcpu1,tKwall1,DECinfo%output)
       
       
       ! By definition:
       !
       ! Ltheta = 1/2 * sum_{ijab} Theta_{ij}^{ab} gx_{ij}^{ab}
       !
       ! Due to conventions in II_get_K_gradient we effectively get
       ! -1/4 sum_{ijab} Theta_{ij}^{ab} gx_{ij}^{ab}
       ! using the procedure above.
       ! Therefore we multiply the final result by (-2)
       Grad%Ltheta(1:3,1:natoms_frag) = -2E0_realk*Grad%Ltheta(1:3,1:natoms_frag)
       ! Free stuff
       call mem_dealloc(dens4index_pair)
       call mem_dealloc(thetaAO_pair)
    ENDIF
    call LSTIMER('LTHETA TOTAL',tcpu,twall,DECinfo%output)

  end subroutine pair_fragment_Ltheta_contribution




  !> \brief Transform virtual indices of Theta array from MO to AO basis
  !> ThetaAO_{IJ}^{alpha beta} =
  !> sum_{CD} Cvirt_{alpha C} Cvirt_{beta D} ThetaMO_{IJ}^{CD}
  !> ThetaAO is also initialized here.
  !> Input order: (C,I,J,D)
  !> Output order: (beta,alpha,J,I)  ( equivalent to (alpha,beta,I,J) )
  !> where D-->beta and C--> alpha.
  !> \author Kasper Kristensen
  !> \date October 2011
  subroutine transform_virtual_Theta_indices_to_AO(ThetaMO, Cvirt, ThetaAO)

    implicit none

    !> Theta MO array Theta_{IJ}^{CD}, only for EOS orbitals using occupied partitioning, order:  (C,I,J,D)
    type(array4),intent(in) :: ThetaMO
    !> Theta AO array Theta_{IJ}^{alpha beta}, order: (beta,alpha,J,I)
    type(array4),intent(inout) :: ThetaAO
    type(array2), intent(inout) :: Cvirt
    type(array4) :: tmp1
    integer :: nocc, nvirt, nbasis

    ! Initialize stuff
    nbasis = Cvirt%dims(1)
    nocc = ThetaMO%dims(2)
    nvirt = Cvirt%dims(2)



    ! *****************************************************************************
    !            Transform virtual Theta indices from MO to AO basis              !
    ! *****************************************************************************

    ! Transpose Cvirt to be able to use array4_contract1
    ! (sum over the first index) to transform Theta from MO to AO basis.
    call array2_transpose(Cvirt)

    ! --------------------------------------------------------------------------------
    !        STEP 1: TRANSFORM INDEX 1 OF ORIGINAL THETA ARRAY IN MO BASIS           !
    ! --------------------------------------------------------------------------------

    ! Dimension of ThetaMO is (nvirt, nocc, nvirt, nocc)
    ! Dimension of tmp1 where we transform the (virtual) first index is:
    ! (nbasis, nocc, nvirt, nocc). I.e. we do the transformation:
    ! (C,I,J,D) --> (alpha,I,J,D)
    tmp1 = array4_init([nbasis,nocc,nocc,nvirt])
    call array4_contract1(ThetaMO,Cvirt,tmp1,.true.)

    ! --------------------------------------------------------------------------------
    !            STEP 2: TRANSFORM INDEX 3 OF ORIGINAL THETA ARRAY IN MO BASIS       !
    ! --------------------------------------------------------------------------------

    ! Reorder tmp1: (alpha,I,J,D) --> (D,alpha,J,I)
    call array4_reorder(tmp1,[4,1,3,2])
    ThetaAO = array4_init([nbasis,nbasis,nocc,nocc])
    ! Contract: (D,alpha,J,I) --> (beta,alpha,J,I)
    call array4_contract1(tmp1,Cvirt,ThetaAO,.true.)
    call array4_free(tmp1)
    ! Transpose Cvirt back
    call array2_transpose(Cvirt)


  end subroutine transform_virtual_Theta_indices_to_AO



  !> \brief Get 4 index density for HF orbitals:
  !> dens4index(mu,nu,i,j) = Cocc(mu,i)*Cocc(nu,j)
  !> \author Kasper Kristensen
  !> \date October 2011
  subroutine get_4_dimensional_HFdensity(Cocc,dens4index)

    implicit none

    type(array4) :: dens4index
    type(array2), intent(inout) :: Cocc
    integer :: nocc, nbasis
    integer :: mu,nu,i,j

    ! Initialize stuff
    nbasis = Cocc%dims(1)
    nocc = Cocc%dims(2)

    dens4index = array4_init([nbasis,nbasis,nocc,nocc])

    ! Simple loop to calculate 4 index density,
    do j=1,nocc
       do i=1,nocc
          do nu=1,nbasis
             do mu=1,nbasis
                dens4index%val(mu,nu,i,j) = Cocc%val(mu,i)*Cocc%val(nu,j)
             end do
          end do
       end do
    end do


  end subroutine get_4_dimensional_HFdensity






  !> Calculate MP2 gradient from fragment contributions.
  !> \author Kasper Kristensen
  !> \date October 2011
  subroutine get_mp2gradient_main(MyMolecule,mylsitem,D,mp2gradient,fullgrad)


    implicit none

    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS setting info
    type(lsitem), intent(inout) :: mylsitem
    !> HF density matrix
    type(matrix),intent(in) :: D
    !> MP2 gradient
    real(realk), intent(inout) :: mp2gradient(3,MyMolecule%natoms)
    !> Full MP2 gradient structure
    type(fullmp2grad),intent(inout) :: fullgrad
    integer :: nvirt,nocc,nbasis,ncore,nval
    type(matrix) :: rho,DMP2,Phivo,Phiov,Phioo
    type(matrix) :: C,F,W
    logical :: full_equation
    real(realk), dimension(3) :: HFdipole, MP2dipole
    real(realk) :: tcpu,twall, tcpu1,tcpu2, twall1,twall2, TrSrho
    real(realk),pointer :: basis(:,:)
    type(matrix) :: Phi

    call LSTIMER('START',tcpu,twall,DECinfo%output)
    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    if(MyMolecule%mem_distributed)then
       call lsquit("ERROR(get_mp2gradient_main) not implemented for distributed matrices in molecule",-1)
    endif


    ! Easy reference to molecule info
    ! *******************************
    nbasis = MyMolecule%nbasis
    nocc = MyMolecule%nocc
    nvirt = MyMolecule%nvirt
    ncore = MyMolecule%ncore
    nval = MyMolecule%nval

    ! Get MP2 gradient matrices in type(matrix) form
    if(DECinfo%frozencore) then
       call convert_mp2gradient_matrices_to_typematrix(mylsitem,MyMolecule,fullgrad,&
            & rho,Phi,Phivo, Phiov,Phioo=Phioo)
    else
       call convert_mp2gradient_matrices_to_typematrix(mylsitem,MyMolecule,fullgrad,&
            & rho,Phi,Phivo, Phiov)
    end if

    ! Get MP2 density
    if(DECinfo%frozencore) then
       call get_full_mp2density(MyMolecule,mylsitem,D,&
            & Phivo, Phiov, DMP2,rho,Phioo=Phioo)
    else
       call get_full_mp2density(MyMolecule,mylsitem,D,&
            & Phivo, Phiov, DMP2,rho)
    end if
!!$    write(DECinfo%output,*) 'Final rho matrix in AO basis'
!!$    call mat_print(rho, 1,rho%nrow, 1, rho%ncol,DECinfo%output)
    call mat_free(Phivo)
    call mat_free(Phiov)
    if(DECinfo%frozencore) call mat_free(Phioo)


    GradientCalc: if(DECinfo%gradient) then ! only for gradient, simple MP2 density

       ! Get MP2 reorthonormalization matrix
       ! ***********************************
       ! MO coefficient matrix and Fock matrix in type matrix form
       call mat_init(C,nbasis,nbasis)
       call mat_init(F,nbasis,nbasis)
       call mem_alloc(basis,nbasis,nbasis)
       basis(1:nbasis,1:nocc) = MyMolecule%Co%elm2(1:nbasis,1:nocc)
       basis(1:nbasis,nocc+1:nbasis) = MyMolecule%Cv%elm2(1:nbasis,1:nvirt)
       call mat_set_from_full(basis(1:nbasis,1:nbasis), 1E0_realk, C)
       call mem_dealloc(basis)
       call mat_set_from_full(MyMolecule%fock%elm2(1:nbasis,1:nbasis), 1E0_realk, F)

       ! Reorthonormalization matrix W
       call util_get_symm_part(rho)
       call get_mp2_reorthonormalization_matrix(F,D,Phi,rho,C,MyLsitem,W)
       call mat_free(C)


       ! Finally, we can determine the MP2 gradient using the calculated matrices
       ! ************************************************************************
       call calculate_MP2_gradient(D,F,rho,W,DMP2,MyLsitem,FullGrad)
       mp2gradient = fullgrad%mp2gradient


       ! Free stuff
       ! **********
       call mat_free(F)
       call mat_free(W)

    else
       write(DECinfo%output,*) 'Only MP2 density requested, skipping MP2 gradient part'

    end if GradientCalc


    call mat_free(DMP2)
    call mat_free(rho)
    call mat_free(Phi)

    call LSTIMER('GET MP2GRAD',tcpu,twall,DECinfo%output)
    call LSTIMER('START',tcpu2,twall2,DECinfo%output)

  end subroutine get_mp2gradient_main



  !> \brief Calculate MP2 gradient by contracting input matrices
  !> against differentiated integrals.
  !> \author Kasper Kristensen
  !> \date November 2010
  subroutine calculate_MP2_gradient(D,F,rho,W,DMP2,MyLsitem,FullGrad)

    implicit none
    !> HF Density matrix in AO basis
    type(matrix), intent(in) :: D
    !> Fock matrix in AO basis
    type(matrix), intent(in) :: F
    !> rho matrix (see dec_get_rho_matrix_in_AO_basis)
    type(matrix), intent(in) :: rho
    !> Reorthogonalization matrix W (see dec_get_gradient_reorthonormalization_matrix)
    type(matrix), intent(in) :: W
    !> DMP2 = 2*D + rho
    type(matrix), intent(in) :: DMP2
    !> Full MP2 gradient matrices (Ltheta is used, Etot is calculated)
    type(FullMP2grad), intent(inout) :: Fullgrad
    !> LSDALTON INFO
    type(lsitem), intent(inout) :: MyLsitem
    ! These are the contributions to the TOTAL MP2 gradient
    real(realk), pointer :: gradient_nuc(:,:), gradient_1el(:,:),&
         &gradient_fock2(:,:), gradient_reort(:,:), gradient_tot(:,:),&
         &gradient_Ltheta(:,:)
    real(realk) :: twall, tcpu,ttotcpu,ttotwall,gradnorm,gradrms
    integer :: nbasis,i,natoms,j,counter,nST
    character(len=80) :: FileName
    character(len=6)  :: MODELSTRING


    call LSTIMER('START',tcpu,twall,DECinfo%output)
    call LSTIMER('START',ttotcpu,ttotwall,DECinfo%output)


    ! Initialize stuff
    ! ****************
    nbasis = D%nrow
    natoms = FullGrad%natoms

    ! Total MP2 gradient contributions
    call mem_alloc(gradient_nuc,3,natoms)
    call mem_alloc(gradient_1el,3,natoms)
    call mem_alloc(gradient_fock2,3,natoms)
    call mem_alloc(gradient_reort,3,natoms)
    call mem_alloc(gradient_Ltheta,3,natoms)
    call mem_alloc(gradient_tot,3,natoms)
    gradient_nuc=0E0_realk
    gradient_1el=0E0_realk
    gradient_fock2=0E0_realk
    gradient_reort=0E0_realk
    gradient_Ltheta=0E0_realk
    gradient_tot=0E0_realk


    ! ************************************************************************************
    !                 Calculate purely nuclear contribution to gradient                  !
    ! ************************************************************************************
    call II_get_nn_gradient(gradient_nuc,MyLsitem%setting,DECinfo%output,DECinfo%output)
    call LSTIMER('GRAD: NUC.',tcpu,twall,DECinfo%output)


    ! ************************************************************************************
    !              Calculate reorthonormalization contribution - Tr(Sx*W)                !
    ! ************************************************************************************
    call get_reorthonormalization_gradient(natoms,gradient_reort,MyLsitem,W)
    call LSTIMER('GRAD: REORT.',tcpu,twall,DECinfo%output)


    ! ************************************************************************************
    !               Calculate one-electron contribution Tr( 2*hx*D + hx*rho )
    ! ************************************************************************************
    call get_one_electron_gradient(natoms,gradient_1el,MyLsitem,DMP2)
    call LSTIMER('GRAD: ONE-EL',tcpu,twall,DECinfo%output)


    ! ************************************************************************************
    !                   Calculate coulomb and exchange contributions                     !
    !                         Tr { [D+rho] * [2*Jx(D) - Kx(D)] }                         !
    ! ************************************************************************************
    call get_Fock2_derivative_gradient(natoms, gradient_fock2, MyLsitem, D, rho)
    call LSTIMER('GRAD: COU+EXC',tcpu,twall,DECinfo%output)


    ! ************************************************************************************
    !                          Ltheta contributions contributions                        !
    ! ************************************************************************************
    ! Just copy Ltheta from FullGrad structure
    do j=1,natoms
       do i=1,3
          gradient_Ltheta(i,j) = FullGrad%Ltheta(i,j)
       end do
    end do



    ! ************************************************************************************
    !                                Add contributions                                   !
    ! ************************************************************************************

    ! TOTAL MP2 gradient (including Hartree-Fock contribution)
    ! 1-electron contribution             (gradient_1el)
    ! reorthonormalization contribution   (gradient_reort)
    ! 2-electron contribution             (gradient_fock2 + gradient_ltheta)
    ! Purely nuclear contribution         (gradient_nuc)
    gradient_tot = gradient_1el + gradient_reort &
         & + gradient_fock2 + gradient_nuc + gradient_ltheta
    ! Set gradient in fullmp2grad structure
    fullgrad%mp2gradient = gradient_tot

    ! Gradient norm
    gradnorm=0.0_realk
    do j=1,natoms
       do i=1,3
          gradnorm = gradnorm + gradient_tot(i,j)**2
       end do
    end do
    gradnorm = sqrt(gradnorm)
    ! RMS for gradient
    gradrms = gradnorm/sqrt(real(3*natoms-1))


    ! ********************************************************************
    !                          Also calculate MP2 energy                 *
    ! ********************************************************************
    ! Total MP2 energy: Hartree-Fock + correlation
    FullGrad%Etot = fullgrad%EHF + FullGrad%Ecorr



    ! ************************************************************************************
    !                                       Print gradient                               !
    ! ************************************************************************************
    !MODIFY FOR MODEL
    SELECT CASE(DECinfo%ccModel)
    CASE(MODEL_MP2)
       nST = 3; MODELSTRING(1:nST)='MP2'
    CASE(MODEL_RIMP2)
       nST = 6; MODELSTRING(1:nST)='RI-MP2'
    CASE DEFAULT
       WRITE (DECinfo%output,'(A)') ' Unknown Gradient Model'
       CALL lsQUIT('Unknown Gradient Model in calculate_MP2_gradient',DECinfo%output)
    END SELECT

    write(DECinfo%output,'(5X,a)') '*****************************************************&
         &***************************'
    write(DECinfo%output,'(5X,a,17X,a,a,18X,a)') '*',MODELSTRING(1:nST),' MOLECULAR GRADIENT FROM DEC CALCULATION','*'
    write(DECinfo%output,'(5X,a)') '*****************************************************&
         &***************************'
    write(DECinfo%output,*)
    write(DECinfo%output,'(5X,a,a,a,g20.8)') 'DEC-',MODELSTRING(1:nST),' GRADIENT NORM: ', gradnorm
    write(DECinfo%output,'(5X,a,a,a,g20.8)') 'DEC-',MODELSTRING(1:nST),' GRADIENT RMS : ', gradrms

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(25X,a,a,a)') 'TOTAL ',MODELSTRING(1:nST),' MOLECULAR GRADIENT'
    write(DECinfo%output,'(25X,a)') "****************************"
    call print_gradient(natoms,gradient_tot,MyLsitem)
    write(DECinfo%output,*)
    write(DECinfo%output,'(25X,a)') 'Individual contributions'
    write(DECinfo%output,'(25X,a)') '************************'

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(25X,a)') 'Nuclear repulsion gradient'
    write(DECinfo%output,'(25X,a)') "--------------------------"
    call print_gradient(natoms,gradient_nuc,MyLsitem)

    write(DECinfo%output,'(24X,a)') 'One-electron integral gradient'
    write(DECinfo%output,'(24X,a)') "------------------------------"
    call print_gradient(natoms,gradient_1el,MyLsitem)

    write(DECinfo%output,'(24X,a)') 'Reorthonormalization gradient'
    write(DECinfo%output,'(24X,a)') "-----------------------------"
    call print_gradient(natoms,gradient_reort,MyLsitem)

    write(DECinfo%output,'(23X,a)') 'Two-electron integral gradient 1'
    write(DECinfo%output,'(20X,a)') "(Fock matrix derivative contribution)"
    write(DECinfo%output,'(20X,a)') "-------------------------------------"
    call print_gradient(natoms,gradient_fock2,MyLsitem)

    write(DECinfo%output,'(23X,a)') 'Two-electron integral gradient 2'
    write(DECinfo%output,'(29X,a)') "(Ltheta contribution)"
    write(DECinfo%output,'(23X,a)') "--------------------------------"
    call print_gradient(natoms,gradient_ltheta,MyLsitem)


    ! Free stuff
    ! **********
    call mem_dealloc(gradient_nuc)
    call mem_dealloc(gradient_1el)
    call mem_dealloc(gradient_fock2)
    call mem_dealloc(gradient_reort)
    call mem_dealloc(gradient_Ltheta)
    call mem_dealloc(gradient_tot)

    call LSTIMER('GRAD: ALL CONT',ttotcpu,ttotwall,DECinfo%output)


  end subroutine calculate_MP2_gradient



  !> \brief Update full mp2 gradient structure with fragment contribution.
  !> \author Kasper Kristensen
  !> \date September 2011
  subroutine update_full_mp2gradient(fraggrad,fullgrad)

    implicit none
    !> Full MP2 gradient structure
    type(fullmp2grad),intent(inout) :: fullgrad
    !> Fragment MP2 gradient contribution
    type(mp2grad),intent(inout) :: fraggrad
    integer :: i,j,ix,jx,a,b,ax,bx


    ! Update full unrelaxed corr. density matrix rho in AO basis
    ! ----------------------------------------------------------
    do j=1,fraggrad%dens%nbasis
       jx=fraggrad%dens%basis_idx(j)
       do i=1,fraggrad%dens%nbasis
          ix=fraggrad%dens%basis_idx(i)
          fullgrad%rho(ix,jx) = fullgrad%rho(ix,jx) + fraggrad%dens%rho(i,j)
       end do
    end do


    ! Update full Phi matrix in AO basis
    ! ----------------------------------
    do j=1,fraggrad%dens%nbasis
       jx=fraggrad%dens%basis_idx(j)
       do i=1,fraggrad%dens%nbasis
          ix=fraggrad%dens%basis_idx(i)
          fullgrad%Phi(ix,jx) = fullgrad%Phi(ix,jx) + fraggrad%PhiAO(i,j)
       end do
    end do


    ! Update correlation energy
    ! -------------------------
    fullgrad%Ecorr = fullgrad%Ecorr + fraggrad%dens%energy


    ! Loop over elements in fragment Ltheta vector and update full Ltheta vector
    ! --------------------------------------------------------------------------
    do j=1,fraggrad%natoms
       jx=fraggrad%atoms_idx(j)
       ! Update full matrix Ltheta element by fragment (i,j) contribution
       fullgrad%Ltheta(1,jx) = fullgrad%Ltheta(1,jx) + fraggrad%Ltheta(1,j)
       fullgrad%Ltheta(2,jx) = fullgrad%Ltheta(2,jx) + fraggrad%Ltheta(2,j)
       fullgrad%Ltheta(3,jx) = fullgrad%Ltheta(3,jx) + fraggrad%Ltheta(3,j)
    end do


  end subroutine update_full_mp2gradient


  !> \brief Get reorthonormalization matrix W for molecular MP2 gradient:
  !> W = 2*D*F*D + rho*F*Dext + D*G(rho)*Dext + 0.5*Phi_AO
  !> where
  !> D: HF density matrix
  !> F: Fock matrix
  !> rho: Effective "MP2 density" (without Hartree Fock density)
  !> Dext = C * C^T: Extended density matrix where the full set of orbitals are used
  !> (in contrast to D where only the occupied orbitals are used).
  !> G(rho) = 2*Coulomb(rho) - Exchange(rho)
  !> Phi_AO: Phi matrix in AO basis
  !> ALL MATRICES IN AO BASIS
  !> \author Kasper Kristensen
  subroutine get_mp2_reorthonormalization_matrix(F,D,Phi,rho,C,MyLsitem,W)


    implicit none
    !> Fock matrix in AO basis
    type(matrix), intent(in) :: F
    !> Density matrix in AO basis
    type(matrix), intent(in) :: D
    !> Phi matrix in AO basis (see subroutine dec_get_Phi_matrix_in_AO_basis)
    type(matrix), intent(in) :: Phi
    !> rho matrix in AO basis (see subroutine dec_get_rho_matrix_in_AO_basis)
    type(matrix), intent(in) :: rho
    !> MO coefficients
    type(matrix), intent(in) :: C
    !> LSDALTON info
    type(lsitem), intent(inout) :: MyLsitem
    !> Reorthonormalization gradient for MP2
    type(matrix), intent(inout) :: W
    type(matrix) :: Dext,tmp1,tmp2
    integer :: nbasis
    real(realk) :: tcpu,twall


    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Initialize stuff
    nbasis = D%nrow
    call mat_init(W,nbasis,nbasis)
    call mat_zero(W)



    ! *******************************************************************
    ! *                           2D*F*D term                           *
    ! *******************************************************************

    ! DFD
    call mat_init(tmp1,nbasis,nbasis)
    call mat_mul(D,F,'n','n',1E0_realk,0E0_realk,tmp1)
    call mat_mul(tmp1,D,'n','n',2E0_realk,0E0_realk,W)


    ! *******************************************************************
    ! *                    rho*F(D)*Dext term                           *
    ! *******************************************************************


    ! Get extended density matrix
    ! '''''''''''''''''''''''''''
    ! For the extended density matrix we use the full set of MO coefficients:
    ! (in contrast to the "normal" density matrix where only occupied MO coefficients are used):
    !
    ! Dext = C C^T
    !
    ! where C = (Cocc|Cvirt).

    ! Note: It is necessary to send in both C and a copy of C (tmp1).
    call mat_assign(tmp1,C)

    call mat_init(Dext,nbasis,nbasis)
    call mat_mul(C,tmp1,'n','t',1E0_realk,0E0_realk,Dext)


    ! Get rho F Dext term
    ! '''''''''''''''''''
    ! tmp1 = rho F
    call mat_mul(rho,F,'n','n',1E0_realk,0E0_realk,tmp1)
    call mat_init(tmp2,nbasis,nbasis)
    ! tmp2 = tmp1 Dext = rho F Dext
    call mat_mul(tmp1,Dext,'n','n',1E0_realk,0E0_realk,tmp2)
    call mat_free(tmp1)
    ! W -> W + tmp2 = 2*DFD +  rho F Dext
    call mat_daxpy(1E0_realk,tmp2, W)
    call mat_free(tmp2)



    ! ***************************************************************
    ! *                        D*G(rho)*Dext                        *
    ! ***************************************************************

    ! 1. Construct tmp1 = G(rho) = 2*Coulomb(rho) - Exchange(rho)
    call mat_init(tmp1,nbasis,nbasis)
    call dec_fock_transformation(tmp1,rho,MyLsitem,.true.)

    ! 2. Calculate D*G(rho)*Dext
    call mat_init(tmp2,nbasis,nbasis)
    call mat_mul(D,tmp1,'n','n',1E0_realk,0E0_realk,tmp2)  ! tmp2 = D*G(rho)
    call mat_mul(tmp2,Dext,'n','n',1E0_realk,0E0_realk,tmp1) ! tmp1 = D*G(rho)*Dext
    call mat_free(tmp2)

    ! W -> W + tmp1 = (2*DFD +  rho F Dext) + D*G(rho)*Dext
    call mat_daxpy(1E0_realk,tmp1, W)

    ! Done with Dext
    call mat_free(tmp1)
    call mat_free(Dext)



    ! *******************************************************************
    ! *                   Phi matrix contribution                       *
    ! *******************************************************************
    ! W -> W + 0.5*Phi = (2*DFD +  rho F Dext + D*G(rho)*Dext) + 0.5*Phi
    call mat_daxpy(0.5E0_realk,Phi, W)


    call LSTIMER('W MATRIX',tcpu,twall,DECinfo%output)


  end subroutine get_mp2_reorthonormalization_matrix



  !> \brief Calculate reorthonormalization contribution to gradient:
  !> - Tr(Sx W)   (total MP2 gradient: Hartree-Fock+correlation)
  !> where Sx is the differentiated overlap matrix, and
  !> W is defined in get_mp2_reorthonormalization_matrix.
  !> \author Kasper Kristensen
  !> \date October 2011
  subroutine get_reorthonormalization_gradient(natoms,gradient,MyLsitem,W)


    implicit none
    !> Number of atoms in the molecule
    integer, intent(in) :: natoms
    !> MP2 reortonormalization gradient
    real(realk), dimension(3,natoms), intent(inout) :: gradient
    !> LSDALTON info
    type(lsitem), intent(inout) :: MyLsitem 
    !> W matrix (see get_mp2_reorthonormalization_matrix)
    type(matrix), target, intent(in) :: W
    type(matrixp) :: Wp(1)

    ! Init stuff
    gradient = 0E0_realk
    Wp(1)%p => W

    ! MP2 reorthonormalization gradient
    call II_get_reorthoNormalization(gradient,Wp,1,&
         MyLsitem%setting,DECinfo%output,DECinfo%output)

    ! MP2 reorthonormalization contribution: -Tr(Sx W)
    gradient = -1.0E0_realk*gradient

  end subroutine get_reorthonormalization_gradient



  !> \brief Get one-electron contribution to MP2 gradient: Tr[ (2*D+rho)hx ]
  !> hx are differentiated one-electron integrals, D is the Hartree-Fock density, and
  !> rho is the MP2 correlation density matrix, see dec_get_rho_matrix_in_AO_basis.
  !> \author Kasper Kristensen
  !> \date October 2011
  subroutine get_one_electron_gradient(natoms,gradient_1el,MyLsitem,DMP2)


    implicit none
    !> Number of atoms in the molecule
    integer, intent(in) :: natoms
    !> One-electron contribution to MP2 gradient
    real(realk), dimension(3,natoms), intent(inout) :: gradient_1el
    !> LS item info
    type(lsitem), intent(inout) :: MyLsitem 
    !> Total MP2 density matrix: HF + MP2 correlation contribution
    type(matrix), target, intent(in) :: DMP2
    type(matrixp) :: DMP2p(1)

    gradient_1el(:,:) = 0E0_realk

    ! Associate pointers
    DMP2p(1)%p => DMP2

    ! 1-electron contribution to MP2 gradient
    call II_get_oneElectron_gradient(gradient_1el,DMP2p,1,&
         MyLsitem%setting,DECinfo%output,DECinfo%output)

  end subroutine get_one_electron_gradient



  !> \brief Get contribution to MP2 gradient from
  !> two-electron part of differentiated Fock matrix: Tr [ (D + rho) Gx(D) ]
  !> where Gx(D) = 2Jx(D) - Kx(D)
  !> is the two-electron part of the Fock matrix with differentiated integrals, i.e.
  !> Jx(D) and Kx(D) are coulomb and exchange contributions
  !> with differentiated two-electron integrals.
  !> The correlation density rho is defined in subroutine dec_get_rho_matrix_in_AO_basis.
  !> \author Kasper Kristensen
  !> \date October 2011
  subroutine get_Fock2_derivative_gradient(natoms, gradient_fock2, MyLsitem, D, rho)


    implicit none
    !> Number of atoms in the molecule
    integer, intent(in) :: natoms
    !> Contribution to gradient from differentiated Fock matrix
    real(realk), dimension(3,natoms), intent(inout) :: gradient_fock2
    !> LSitem info
    type(lsitem), intent(inout) :: MyLsitem
    !> Hartree-Fock density matrix
    type(matrix), target, intent(in) :: D
    !> Correlation density matrix
    type(matrix), intent(in) :: rho
    type(matrix), target :: Dtot
    type(matrixp) :: Dtotp(1),Dp(1)
    real(realk), dimension(3,natoms) :: CorrGradient_fock2
    integer :: nbasis



    ! Initialize stuff
    gradient_fock2(:,:) = 0E0_realk
    nbasis = D%nrow
    call mat_init(Dtot,nbasis,nbasis)
    call mat_add(1E0_realk,D,1E0_realk,rho,Dtot)

    ! Associate pointer
    Dtotp(1)%p => Dtot
    Dp(1)%p => D

    ! Coulomb+exchange contributions to MP2 gradient
    call II_get_twoElectron_gradient(gradient_fock2,natoms,Dp,Dtotp, &
         1,1,MyLsitem%setting,DECinfo%output,DECinfo%output)

    ! Due to conventions in II_get_twoElectron_gradient we have to multiply
    ! the gradient contribution by 4 (only works for closed-shell systems)
    gradient_fock2 = 4E0_realk*gradient_fock2
    call mat_free(Dtot)

  end subroutine get_Fock2_derivative_gradient


  !> \brief Print molecular gradient
  !> \author Kasper Kristensen
  !> \date November 2010
  subroutine print_gradient(natoms,grad,MyLsitem)
    implicit none
    !> Number of atoms in the molecule
    integer :: natoms
    !> Gradient to be printed
    real(realk), dimension(3,natoms) :: grad
    !> LSDALTON INFO
    type(lsitem), intent(inout) :: MyLsitem
    integer :: i

    write(DECinfo%output,*)
    do i=1,natoms
       write(DECinfo%output,'(4X,a4,3F20.10)') MyLsitem%input%molecule%atom(i)%name, &
            grad(1:3,i)
    enddo
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)

  end subroutine print_gradient



  !> \brief Write intermediates for full MP2 gradient structure and
  !> fragment energies to file mp2grad.info for easy restart.
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine write_gradient_and_energies_for_restart(natoms,FragEnergies,jobs,fullgrad)

    implicit none
    !> Number of atoms in the molecule
    integer,intent(in) :: natoms
    !> Fragment energies (see decfrag type def)
    real(realk),dimension(natoms,natoms,ndecenergies),intent(in) :: FragEnergies
    !> Job list of fragments
    type(joblist),intent(in) :: jobs
    !> Full MP2 gradient structure
    type(fullmp2grad),intent(in) :: fullgrad
    character(len=40) :: FileName
    integer :: funit,i
    logical :: file_exist

    ! Init stuff
    funit = -1
    FileName='mp2grad.info'

    ! backup exisiting file
    inquire(file=FileName,exist=file_exist)
    if(file_exist) then
#ifdef SYS_AIX
       call rename('mp2grad.info\0','mp2grad.backup\0')
#else
       call rename('mp2grad.info','mp2grad.backup')
#endif
    end if

    ! Create a new file mp2grad.info
    call lsopen(funit,FileName,'NEW','UNFORMATTED')


    ! Write info to file
    ! ******************
    
    ! Job list and energies
    call basic_write_jobs_and_fragment_energies_for_restart(natoms,FragEnergies,jobs,funit,filename)

    ! Gradient stuff
    write(funit) fullgrad%nocc
    write(funit) fullgrad%nvirt
    write(funit) fullgrad%natoms
    write(funit) fullgrad%EHF
    write(funit) fullgrad%Ecorr
    write(funit) fullgrad%Etot
    write(funit) fullgrad%rho
    write(funit) fullgrad%Phi
    write(funit) fullgrad%Ltheta
    write(funit) fullgrad%mp2gradient
    call lsclose(funit,'KEEP')

  end subroutine write_gradient_and_energies_for_restart



  !> \brief Read intermediates for full MP2 gradient structure and
  !> fragment energies from file mp2grad.info for easy restart.
  !> \author Kasper Kristensen
  !> \date June 2012
  subroutine read_gradient_and_energies_for_restart(nfrags,FragEnergies,jobs,fullgrad)

    implicit none
    !> Number of fragments
    integer,intent(in) :: nfrags
    !> Fragment energies (see decfrag type def)
    real(realk),dimension(nfrags,nfrags,ndecenergies),intent(inout) :: FragEnergies
    !> Job list of fragments
    type(joblist),intent(inout) :: jobs
    !> Full MP2 gradient structure
    !> (the values are changed here, but it is assumed that the matrices have been initialized)
    type(fullmp2grad),intent(inout) :: fullgrad
    character(len=40) :: FileName
    integer :: funit
    logical :: file_exist

    ! Init stuff
    funit = -1
    FileName='mp2grad.info'

    ! Check if file exists
    inquire(file=FileName,exist=file_exist)
    if(.not. file_exist) then
       write(DECinfo%output,*) 'WARNING: Restart file: ', FileName
       write(DECinfo%output,*) 'does not exist! I therefore calculate it from scratch...'
       return
    end if

    ! Open file mp2grad.info
    call lsopen(funit,FileName,'OLD','UNFORMATTED')


    ! Read info from file
    ! *******************

    ! Job and energy info
    call basic_read_jobs_and_fragment_energies_for_restart(nfrags,FragEnergies,jobs,funit,FileName)

    ! Gradient stuff
    if(DECinfo%convert64to32) then
       call read_64bit_to_32bit(funit,fullgrad%nocc)
       call read_64bit_to_32bit(funit,fullgrad%nvirt)
       call read_64bit_to_32bit(funit,fullgrad%natoms)
    elseif(DECinfo%convert32to64) then
       call read_32bit_to_64bit(funit,fullgrad%nocc)
       call read_32bit_to_64bit(funit,fullgrad%nvirt)
       call read_32bit_to_64bit(funit,fullgrad%natoms)
    else
       read(funit) fullgrad%nocc
       read(funit) fullgrad%nvirt
       read(funit) fullgrad%natoms
    end if

    read(funit) fullgrad%EHF
    read(funit) fullgrad%Ecorr
    read(funit) fullgrad%Etot
    read(funit) fullgrad%rho
    read(funit) fullgrad%Phi
    read(funit) fullgrad%Ltheta
    read(funit) fullgrad%mp2gradient
    call lsclose(funit,'KEEP')


  end subroutine read_gradient_and_energies_for_restart

  
  !> \brief Convert rho and Phi matrices in fullmp2grad structure to type(matrix) format.
  !> Matrices are also initialized here.
  !> \author Kasper Kristensen
  !> \date December 2012
  subroutine convert_mp2gradient_matrices_to_typematrix(mylsitem,MyMolecule,fullgrad,rho,Phi,&
       & Phivo, Phiov,Phioo)
    implicit none

    !> LS setting info
    type(lsitem), intent(inout) :: mylsitem
    !> Full molecular info
    type(fullmolecule) :: MyMolecule
    !> MP2 gradient structure
    type(fullmp2grad),intent(in) :: fullgrad
    !> rho matrix in AO basis
    type(matrix),intent(inout) :: rho
    !> Phi matrix in AO basis
    type(matrix),intent(inout) :: Phi
    !> Virt-occ blocks of Phi matrix
    type(matrix),intent(inout) :: Phivo
    !> Occ-virt blocks of Phi matrix
    type(matrix),intent(inout) :: Phiov
    !> Occ-occ blocks of Phi matrix (only for frozen core)
    type(matrix),intent(inout),optional :: Phioo


    ! rho and Phi in AO basis
    call mat_init(rho,fullgrad%nbasis,fullgrad%nbasis)
    call mat_set_from_full(fullgrad%rho, 1.0_realk, rho)
    call mat_init(Phi,fullgrad%nbasis,fullgrad%nbasis)
    call mat_set_from_full(fullgrad%Phi, 1.0_realk, Phi)

    ! Phi in MO basis
    if(DECinfo%frozencore) then
       if(.not. present(Phioo) ) then
          call lsquit('convert_mp2gradient_matrices_to_typematrix: &
               & Phioo must be present for frozen core!',-1)
       end if
       call get_Phi_MO_blocks(mylsitem,MyMolecule,fullgrad,Phivo,Phiov,Phioo=Phioo)
    else
       call get_Phi_MO_blocks(mylsitem,MyMolecule,fullgrad,Phivo,Phiov)
    end if

  end subroutine convert_mp2gradient_matrices_to_typematrix


  !> \brief Get intrinsic DEC energy error for geometry optimization
  !>
  !> UNDER INVESTIGATION!!!
  !> 
  !> It is not yet clear what the best strategy is here.
  !> We could simply take the estimated error as it is ( DECinfo%EerrFactor = 1)
  !> or scale it ( DECinfo%EerrFactor /= 1)
  !> or compare it to the error at the previous geometry (code currently commented out).
  !> \author Kasper Kristensen
  !> \date December 2012
  subroutine dec_get_error_for_geoopt(Eerr)
    implicit none
    !> Input: Estimated intrinsic DEC energy error
    !> Output: Absolute difference between intrinsic energy error from this
    !> and the previous geometry.
    real(realk),intent(inout) :: Eerr
    real(realk) :: Eerrsave

    ! Scale intrinsic error 
    Eerr = DECinfo%EerrFactor*Eerr
    Eerrsave = Eerr

    ! Save existing energy error in DECinfo%EerrOLD
    DECinfo%EerrOLD = Eerrsave

!!$    ! Energy error returned to optimizer is difference between error at this 
!!$    ! geometry and the previous geometry.
!!$    if( DECinfo%ncalc(DECinfo%FOTlevel)==0 ) then  
!!$       ! This is the very first gradient calculation for the current FOT level.
!!$       ! Therefore, we don't have an error at a different geometry to compare against,
!!$       ! and we simply set error to zero to avoid artefacts for the dynamic geometry optimizer
!!$       Eerr = 0.0_realk
!!$    else
!!$       ! Set error equal to difference in errors between this and the previous geometry.
!!$       Eerr = abs(Eerr - DECinfo%EerrOLD)
!!$    end if
!!$
!!$    write(DECinfo%output,'(1X,a)') 'DEC STABILITY'
!!$    write(DECinfo%output,'(1X,a,g20.10)') 'DEC STABILITY: Intrinsic absolute error   :', Eerrsave
!!$    if( DECinfo%ncalc(DECinfo%FOTlevel)/=0 ) then
!!$       write(DECinfo%output,'(1X,a,g20.10)') 'DEC STABILITY: Intrinsic error difference :', Eerr
!!$    end if

  end subroutine dec_get_error_for_geoopt



  ! MP2 DENSITY STUFF - SUBSET OF GRADIENT
  ! ======================================


  !> \brief Init MP2 density structure for single fragment, assumes that the fragment structure (decfrag)
  !> has already been initiated!
  !> \author Kasper Kristensen
  !> \date September 2011
  subroutine single_init_mp2dens(fragment,dens)

    implicit none
    !> Fragment with information corresponding to density
    type(decfrag),intent(inout) :: fragment
    !> MP2 density for fragment
    type(mp2dens), intent(inout) :: dens
    integer :: i

    write(DECinfo%output,*) 'Initiating MP2 density matrix...'


    ! Central atom and orbital space sizes
    ! ************************************
    dens%centralatom = fragment%EOSatoms(1)
    dens%centralatom2 = 0    ! only used for pairs
    dens%nbasis = fragment%nbasis
    dens%nvirt = fragment%nvirtAOS
    dens%nocc = fragment%noccAOS
    dens%nocctot = fragment%nocctot
    dens%energy = fragment%energies(FRAGMODEL_OCCMP2)
    dens%pairdist = 0E0_realk
    dens%nEOSatoms = fragment%nEOSatoms

    call mem_alloc(dens%EOSatoms,dens%nEOSatoms)
    do i=1,dens%nEOSatoms
       dens%EOSatoms(i) = fragment%EOSatoms(i)
    end do

    ! Sanity check
    ! ************
    if( (dens%nvirt==0) .or. (dens%nocc==0) ) then
       write(DECinfo%output,*) 'nvirt,nocc =', dens%nvirt,dens%nocc
       call lsquit('single_init_mp2dens: Number of orbitals is zero, it seems that the fragment &
            & has not been initialized before initiating the density',-1)
    end if


    ! Atomic orbital indices
    ! **********************
    call mem_alloc(dens%basis_idx,dens%nbasis)
    dens%basis_idx = fragment%basis_idx

    ! Initiate array for virt-virt block of density (Y)
    ! *************************************************
    ! Dimension of Y is (nvirt,nvirt)
    call mem_alloc(dens%Y,dens%nvirt,dens%nvirt)
    dens%Y=0E0_realk

    ! Initiate array for occ-occ block of density (X)
    ! ***********************************************
    ! Dimension of X is (nocc,nocc)
    call mem_alloc(dens%X,dens%nocc,dens%nocc)
    dens%X=0E0_realk

    ! Initiate array for unrelaxed corr dens in AO basis
    ! **************************************************
    call mem_alloc(dens%rho,dens%nbasis,dens%nbasis)
    dens%rho=0E0_realk


    ! Initiate array for virt-occ block of Phi matrix
    ! ***********************************************
    ! Dimension of Phivo is (nvirt,nocctot)
    call mem_alloc(dens%Phivo,dens%nvirt,dens%nocctot)
    dens%Phivo=0E0_realk

    ! Initiate array for occ-virt block of Phi matrix
    ! ***********************************************
    ! Dimension of Phiov is (nocc.nvirt)
    call mem_alloc(dens%Phiov,dens%nocc,dens%nvirt)
    dens%Phiov=0E0_realk

  end subroutine single_init_mp2dens



  !> \brief Init MP2 density structure for pair fragment, assumes that the fragment
  !> has already been initiated!
  !> \author Kasper Kristensen
  !> \date August 2011
  subroutine pair_init_mp2dens(pairfragment,dens,atom1,atom2)

    implicit none
    !> Pair fragment with information corresponding to density
    type(decfrag),intent(inout) :: pairfragment
    !> MP2 density for pair fragment
    type(mp2dens), intent(inout) :: dens
    !> First atom in pair
    integer, intent(in) :: atom1
    !> Second atom in pair
    integer, intent(in) :: atom2

    ! Simply init as in corresponding routine for single fragment,
    ! then set the atomic site information and pair distance for the pair.
    call single_init_mp2dens(pairfragment,dens)
    dens%centralatom = atom1
    dens%centralatom2 = atom2
    dens%pairdist = pairfragment%pairdist

  end subroutine pair_init_mp2dens



  !> \brief Free MP2 density structure
  !> \author Kasper Kristensen
  !> \date August 2011
  subroutine free_mp2dens(dens)

    implicit none
    type(mp2dens), intent(inout) :: dens

    ! Set stuff to zero
    dens%nocc=0
    dens%nocctot=0
    dens%nvirt=0
    dens%centralatom=0
    dens%centralatom2=0
    dens%energy=0E0_realk
    dens%pairdist=0E0_realk

    ! Deallocate vectors/arrays
    if(associated(dens%EOSatoms)) call mem_dealloc(dens%EOSatoms)
    if(associated(dens%basis_idx)) call mem_dealloc(dens%basis_idx)
    if(associated(dens%Y)) call mem_dealloc(dens%Y)
    if(associated(dens%X)) call mem_dealloc(dens%X)
    if(associated(dens%rho)) call mem_dealloc(dens%rho)
    if(associated(dens%Phivo)) call mem_dealloc(dens%Phivo)
    if(associated(dens%Phiov)) call mem_dealloc(dens%Phiov)

  end subroutine free_mp2dens


  !> \brief Calculate contributions to MP2 density from t2 amplitudes for single fragment.
  !> See mp2dens structure for equations.
  !> \author Kasper Kristensen
  !> \date August 2011
  subroutine single_calculate_mp2density(MyFragment,t2occ,t2virt,ThetaOCC,ThetaVIRT,VOOO,VOVV,dens)


    implicit none

    !> Atomic fragment
    type(decfrag),intent(inout) :: MyFragment
    !> t2 amplitudes t_{IJ}^{CD}, only for EOS orbitals using occupied partitioning, order:  (C,I,J,D)
    type(array4),intent(in) :: t2occ
    !> t2 amplitudes t_{KL}^{AB}, only for EOS orbitals using virtual partitioning, order: (A,K,B,L)
    type(array4),intent(in) :: t2virt
    !> Theta array, only for EOS orbitals using occupied partitioning, order:  (C,I,J,D)
    type(array4),intent(in) :: ThetaOCC
    !> Theta array, only for EOS orbitals using virtual partitioning, order:  (A,K,B,L)
    type(array4),intent(in) :: ThetaVIRT
    !> (C I | J L) integrals stored as (C,I,J,L)    [see index conventions in mp2.f90]
    type(array4),intent(in) :: VOOO
    !> (B K | A C) integrals stored as (B,K,A,C)    [see index conventions in mp2.f90]
    type(array4),intent(in) :: VOVV
    !> MP2 density matrix contribution from given single fragment
    type(mp2dens),intent(inout) :: dens
    type(array2) :: Phivo, Phiov, X, Y
    integer :: noccEOS,nvirtEOS,noccAOS,nvirtAOS,nocctot
    logical :: something_wrong
    real(realk) :: tcpu,twall

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Init stuff
    noccEOS = t2occ%dims(2) ! number of occupied EOS orbitals
    noccAOS = t2virt%dims(2) ! number of occupied AOS orbitals
    nvirtEOS = t2virt%dims(1) ! number of virtual EOS orbitals
    nvirtAOS = t2occ%dims(1) ! number of virtual AOS orbitals

    ! number of occupied core+valence orbitals (only different from noccAOS for frozen approx)
    if(DECinfo%unrelaxed) then
       nocctot = dens%nocctot   ! not used for unrelaxed density
    else
       nocctot = VOOO%dims(4)   
    end if

    ! Just in case, zero matrices in dens structure
    dens%X = 0e0_realk
    dens%Y = 0e0_realk
    dens%Phivo = 0e0_realk
    dens%Phiov = 0e0_realk



    ! Sanity check 1: Density input structure is correct
    something_wrong=.false.
    if(noccAOS/=dens%nocc) something_wrong=.true.
    if(nvirtAOS/=dens%nvirt) something_wrong=.true.
    if(something_wrong) then
       write(DECinfo%output,*) 'dens%nocc', dens%nocc
       write(DECinfo%output,*) 'dens%nvirt', dens%nvirt
       write(DECinfo%output,*) 'Amplitudes, occ AOS dimension ', noccAOS
       write(DECinfo%output,*) 'Amplitudes, virt AOS dimension', nvirtAOS
       call lsquit('single_calculate_mp2density: &
            & Something wrong with input dimensions!',DECinfo%output)
    end if

    ! Sanity check 2: Array dimensions match for last two indices of t2 and Theta
    if(noccEOS /= t2occ%dims(3)) something_wrong=.true.
    if(nvirtAOS /= t2occ%dims(4)) something_wrong=.true.
    if(nvirtEOS /= t2virt%dims(3)) something_wrong=.true.
    if(noccAOS /= t2virt%dims(4)) something_wrong=.true.
    if(noccEOS /= Thetaocc%dims(3)) something_wrong=.true.
    if(nvirtAOS /= Thetaocc%dims(4)) something_wrong=.true.
    if(nvirtEOS /= Thetavirt%dims(3)) something_wrong=.true.
    if(noccAOS /= Thetavirt%dims(4)) something_wrong=.true.
    if(something_wrong) then
       write(DECinfo%output,*) 'EOS: nocc, nvirt', noccEOS, nvirtEOS
       write(DECinfo%output,*) 'AOS: nocc, nvirt', noccAOS, nvirtAOS
       write(DECinfo%output,*) 't2occ dims     ', t2occ%dims
       write(DECinfo%output,*) 'Thetaocc dims  ', Thetaocc%dims
       write(DECinfo%output,*) 't2virt dims    ', t2virt%dims
       write(DECinfo%output,*) 'Thetavirt dims ', Thetavirt%dims
       call lsquit('single_calculate_mp2density: &
            & Something wrong with amplitude and Theta dimensions!',DECinfo%output)
    end if


    ! ******************************************************************************************
    !                Calculate contribution to virt-virt block of density matrix (Y)           *
    ! ******************************************************************************************

    ! Use occupied partitioning scheme
    ! Y_{ab} = sum_{cij} t_{ji}^{ca} * mult_{ji}^{cb}           [see the mp2dens structure]
    !        = 1/2* sum_{cij} t_{ji}^{ca} * Theta_{ji}^{cb}     [see construct_theta_array]
    Y = array2_init([nvirtAOS,nvirtAOS])

    call array4_contract3(t2occ,ThetaOCC,Y)
    dens%Y(1:nvirtAOS,1:nvirtAOS) = 0.5e0_realk*Y%val(1:nvirtAOS,1:nvirtAOS)

    call array2_free(Y)
    call LSTIMER('Y MATRIX',tcpu,twall,DECinfo%output)



    ! ****************************************************************************************
    !                Calculate contribution to occ-occ block of density matrix (X)           *
    ! ****************************************************************************************

    ! Use virtual partitioning scheme
    ! X_{ij} = sum_{abk} t_{ki}^{ba} * mult_{kj}^{ba}          [see the mp2dens structure]
    !        = 1/2* sum_{abk} t_{ki}^{ba} * Theta_{kj}^{ba}    [see construct_theta_array]
    X = array2_init([noccAOS,noccAOS])

    call array4_contract3(t2virt,ThetaVIRT,X)
    dens%X(1:noccAOS,1:noccAOS) = 0.5e0_realk*X%val(1:noccAOS,1:noccAOS)

    call array2_free(X)
    call LSTIMER('X MATRIX',tcpu,twall,DECinfo%output)

    ! X and Y matrices collected and transformed to AO basis (unrelaxed corr. density)
    call get_unrelaxed_corrdens_in_AO_basis(dens%nocc,dens%nvirt,dens%nbasis,MyFragment%Co,&
         & MyFragment%Cv,dens%X,dens%Y,dens%rho)

    UNRELAXED: if(DECinfo%unrelaxed) then
       ! Skip Phi matrix
       dens%Phivo=0.0E0_realk
       dens%Phiov=0.0E0_realk

    else

       ! ****************************************************************************************
       !                  Calculate contribution to virt-occ block of Phi matrix                *
       ! ****************************************************************************************


       ! Phivo(D,L) = sum_{CIJ} Theta(C,I,J,D) g(C,I,J,L)     [see the mp2dens structure]
       ! Note: L is both core+valence, also for frozen core approximation.
       Phivo = array2_init([nvirtAOS,nocctot])
       call array4_contract3(ThetaOCC,VOOO,Phivo)
       dens%Phivo(1:nvirtAOS,1:nocctot) = Phivo%val(1:nvirtAOS,1:nocctot)

       call array2_free(Phivo)
       call LSTIMER('PHIVO MATRIX',tcpu,twall,DECinfo%output)




       ! ****************************************************************************************
       !                  Calculate contribution to occ-virt block of Phi matrix                *
       ! ****************************************************************************************


       ! Phiov(L,C) = sum_{ABK} Theta(B,K,A,L) g(B,K,A,C)      [see the mp2dens structure]
       Phiov = array2_init([noccAOS,nvirtAOS])
       call array4_contract3(ThetaVIRT,VOVV,Phiov)
       dens%Phiov(1:noccAOS,1:nvirtAOS) = Phiov%val(1:noccAOS,1:nvirtAOS)

       call array2_free(Phiov)

       call LSTIMER('PHIOV MATRIX',tcpu,twall,DECinfo%output)

    end if UNRELAXED


  end subroutine single_calculate_mp2density



  !> \brief Calculate contributions to MP2 density from t2 amplitudes for pair fragment.
  !> See mp2dens structure for equations.
  !> \author Kasper Kristensen
  !> \date August 2011
  subroutine pair_calculate_mp2density(fragment1,fragment2,PairFragment,t2occ,t2virt,&
       & ThetaOCC, ThetaVIRT, VOOO,VOVV,dens)


    implicit none
    !> Fragment 1 in the pair fragment
    type(decfrag),intent(in) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(decfrag),intent(in) :: Fragment2
    !> Pair fragment
    type(decfrag),intent(inout) :: pairfragment
    !> t2 amplitudes, only for EOS orbitals using occupied partitioning, order:  (A,I,B,J)
    type(array4),intent(in) :: t2occ
    !> t2 amplitudes, only for EOS orbitals using virtual partitioning, order:  (A,I,B,J)
    type(array4),intent(in) :: t2virt
    !> Theta array, only for EOS orbitals using occupied partitioning, order:  (A,I,B,J)
    type(array4),intent(inout) :: ThetaOCC
    !> Theta array, only for EOS orbitals using virtual partitioning, order:  (A,I,B,J)
    type(array4),intent(in) :: ThetaVIRT
    !> (C I | J L) integrals stored as (C,I,J,L)    [see index conventions in mp2.f90]
    type(array4),intent(in) :: VOOO
    !> (B K | A C) integrals stored as (B,K,A,C)    [see index conventions in mp2.f90]
    type(array4),intent(in) :: VOVV
    !> MP2 density matrix contribution from given pair fragment
    type(mp2dens),intent(inout) :: dens
    logical, pointer :: dopair_occ(:,:), dopair_virt(:,:)
    integer :: i,j,a,b,c,k,d,l,noccEOS,nvirtEOS,noccAOS,nvirtAOS,natoms,nocctot
    logical :: something_wrong
    real(realk) :: tcpu, twall
    real(realk),pointer :: Y_tmp(:,:),X_tmp(:,:), Phivo_tmp(:,:), Phiov_tmp(:,:)

    ! Init stuff
    noccEOS = t2occ%dims(2) ! number of occupied EOS orbitals
    noccAOS = t2virt%dims(2) ! number of occupied AOS orbitals
    nvirtEOS = t2virt%dims(1) ! number of virtual EOS orbitals
    nvirtAOS = t2occ%dims(1) ! number of virtual AOS orbitals

    ! number of occupied core+valence orbitals (only different from noccAOS for frozen approx)
    if(DECinfo%unrelaxed) then
       nocctot = dens%nocctot   ! not used for unrelaxed density
    else
       nocctot = VOOO%dims(4)   
    end if

    ! Which "interaction pairs" to include for occ and virt space (avoid double counting)
    call mem_alloc(dopair_occ,noccEOS,noccEOS)
    call mem_alloc(dopair_virt,nvirtEOS,nvirtEOS)
    call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)
    call which_pairs_virt(Fragment1,Fragment2,PairFragment,dopair_virt)

    ! Just in case, zero matrices in dens structure
    dens%X = 0e0_realk
    dens%Y = 0e0_realk
    dens%Phivo = 0e0_realk
    dens%Phiov = 0e0_realk


    ! Sanity check
    something_wrong=.false.
    if(noccAOS/=dens%nocc) something_wrong=.true.
    if(nvirtAOS/=dens%nvirt) something_wrong=.true.
    if(nocctot/=dens%nocctot) something_wrong=.true.
    if(something_wrong) then
       write(DECinfo%output,*) 'dens%nocc', dens%nocc
       write(DECinfo%output,*) 'dens%nvirt', dens%nvirt
       write(DECinfo%output,*) 'dens%nocctot', dens%nocctot
       write(DECinfo%output,*) 'Amplitudes, occ AOS dimension ', noccAOS
       write(DECinfo%output,*) 'Amplitudes, virt AOS dimension', nvirtAOS
       call lsquit('pair_calculate_mp2density: &
            & Something wrong with input dimensions!',DECinfo%output)
    end if



    ! ******************************************************************************************
    !                Calculate contribution to virt-virt block of density matrix (Y)           *
    ! ******************************************************************************************

    ! Use occupied partitioning scheme
    ! Y_{ab} = sum_{cij} t_{ji}^{ca} * mult_{ji}^{cb}           [see the mp2dens structure]
    !        = 1/2* sum_{cij} t_{ji}^{ca} * Theta_{ji}^{cb}     [see construct_theta_array]
    ! -- where we only include the contributions if i and j belong to different fragments!
    call LSTIMER('START',tcpu,twall,DECinfo%output)

call mem_TurnONThread_Memory()
!$OMP PARALLEL DEFAULT(shared) PRIVATE(Y_tmp,i,j,a,b,c)
call init_threadmemvar()

    call mem_alloc(Y_tmp,nvirtAOS,nvirtAOS)
    Y_tmp = 0E0_realk

!$OMP DO SCHEDULE(dynamic,1)

    do i=1,noccEOS
       do j=1,noccEOS

          ! Only update for "interaction orbital pairs" - see which_pairs_occ
          if(dopair_occ(i,j)) then !YDoPair

             do b=1,nvirtAOS
                do c=1,nvirtAOS
                   do a=1,nvirtAOS
                      Y_tmp(a,b) = Y_tmp(a,b) + t2occ%val(c,j,a,i)*ThetaOCC%val(c,j,b,i)
                   end do
                end do
             end do

          end if

       end do
    end do
!$OMP END DO NOWAIT

! Total Y matrix is found by summing all thread contributions
!$OMP CRITICAL
dens%Y = dens%Y + Y_tmp
!$OMP END CRITICAL

call mem_dealloc(Y_tmp)
call collect_thread_memory()
!$OMP END PARALLEL
call mem_TurnOffThread_Memory()

! Multiply by 1/2 due to convention for Theta array [see above]
dens%Y = 0.5E0_realk*dens%Y
    call LSTIMER('Y MATRIX',tcpu,twall,DECinfo%output)



    ! ****************************************************************************************
    !                Calculate contribution to occ-occ block of density matrix (X)           *
    ! ****************************************************************************************

    ! Use virtual partitioning scheme
    ! X_{ij} = sum_{abk} t_{ki}^{ba} * mult_{kj}^{ba}          [see the mp2dens structure]
    !        = 1/2* sum_{abk} t_{ki}^{ba} * Theta_{kj}^{ba}    [see construct_theta_array]
    ! -- where we only include the contributions if A and B belong to different fragments!

Call mem_TurnONThread_Memory()
!$OMP PARALLEL DEFAULT(shared) PRIVATE(X_tmp,i,j,k,a,b)
call init_threadmemvar()
    call mem_alloc(X_tmp,noccAOS,noccAOS)
    X_tmp = 0E0_realk

!$OMP DO SCHEDULE(dynamic,1)


    do a=1,nvirtEOS
       do b=1,nvirtEOS

          ! Only update for "interaction orbital pairs" - see which_pairs_virt
          if(dopair_virt(a,b)) then !XDoPair

             do j=1,noccAOS
                do k=1,noccAOS
                   do i=1,noccAOS
                      X_tmp(i,j) = X_tmp(i,j) + t2virt%val(b,k,a,i)*ThetaVIRT%val(b,k,a,j)
                   end do
                end do
             end do

          end if

       end do
    end do
!$OMP END DO NOWAIT

! Total X matrix is found by summing all thread contributions
!$OMP CRITICAL
dens%X = dens%X + X_tmp
!$OMP END CRITICAL

call mem_dealloc(X_tmp)
call collect_thread_memory()
!$OMP END PARALLEL
call mem_TurnOffThread_Memory()

! Multiply by 1/2 due to convention for Theta array [see above]
    dens%X = 0.5E0_realk*dens%X
    call LSTIMER('X MATRIX',tcpu,twall,DECinfo%output)




    ! ****************************************************************************************
    !                  Calculate contribution to virt-occ block of Phi matrix                *
    ! ****************************************************************************************


    ! Phivo(d,l) = sum_{cij} Theta(c,i,d,j) g(c,i,j,l)        [see the mp2dens structure]
    ! -- where we only include the contributions if I and J belong to different fragments!

! Skip for unrelaxed density
UNRELAXED1: if(DECinfo%unrelaxed) then
dens%Phivo = 0.0E0_realk

else

call mem_TurnONThread_Memory()
!$OMP PARALLEL DEFAULT(shared) PRIVATE(Phivo_tmp,i,j,c,l,d)
call init_threadmemvar()
    call mem_alloc(Phivo_tmp,nvirtAOS,nocctot)
    Phivo_tmp = 0E0_realk

!$OMP DO SCHEDULE(dynamic,1)


    do i=1,noccEOS
       do j=1,noccEOS

          if(dopair_occ(i,j)) then !PhivoDoPair

             do d=1,nvirtAOS
                do l=1,nocctot
                   do c=1,nvirtAOS
                      Phivo_tmp(d,l) = Phivo_tmp(d,l) + ThetaOCC%val(c,i,d,j)*VOOO%val(c,i,j,l)
                   end do
                end do
             end do

          end if

       end do
    end do
!$OMP END DO NOWAIT

! Total Phivo matrix is found by summing all thread contributions
!$OMP CRITICAL
dens%Phivo = dens%Phivo + Phivo_tmp
!$OMP END CRITICAL

    call mem_dealloc(Phivo_tmp)
call collect_thread_memory()
!$OMP END PARALLEL
call mem_TurnOffThread_Memory()

    call LSTIMER('PHIVO MATRIX',tcpu,twall,DECinfo%output)

end if UNRELAXED1



    ! ****************************************************************************************
    !                  Calculate contribution to occ-virt block of Phi matrix                *
    ! ****************************************************************************************

    ! Phiov(l,c) = sum_{abk} Theta(b,k,a,l) g(b,k,a,c)      [see the mp2dens structure])
    ! -- where we only include the contributions if A and B belong to different fragments!

! Skip for unrelaxed density
UNRELAXED2: if(DECinfo%unrelaxed) then
dens%Phiov = 0.0E0_realk

else

call mem_TurnONThread_Memory()
!$OMP PARALLEL DEFAULT(shared) PRIVATE(Phiov_tmp,a,b,c,k,l)
call init_threadmemvar()
    call mem_alloc(Phiov_tmp,noccAOS,nvirtAOS)
    Phiov_tmp = 0E0_realk

!$OMP DO SCHEDULE(dynamic,1)

    do a=1,nvirtEOS
       do b=1,nvirtEOS

          if(dopair_virt(a,b)) then ! PhiovDoPair

             do l=1,noccAOS
                do c=1,nvirtAOS
                   do k=1,noccAOS
                      Phiov_tmp(l,c) = Phiov_tmp(l,c) + ThetaVIRT%val(b,k,a,l)*VOVV%val(b,k,a,c)
                   end do
                end do
             end do

          end if

       end do
    end do
!$OMP END DO NOWAIT

! Total Phiov matrix is found by summing all thread contributions
!$OMP CRITICAL
dens%Phiov = dens%Phiov + Phiov_tmp
!$OMP END CRITICAL

    call mem_dealloc(Phiov_tmp)
call collect_thread_memory()
!$OMP END PARALLEL
call mem_TurnOffThread_Memory()

end if UNRELAXED2

    call mem_dealloc(dopair_occ)
    call mem_dealloc(dopair_virt)


    ! X and Y matrices collected and transformed to AO basis (unrelaxed corr. density)
    call get_unrelaxed_corrdens_in_AO_basis(dens%nocc,dens%nvirt,dens%nbasis,PairFragment%Co,&
         & PairFragment%Cv,dens%X,dens%Y,dens%rho)

    call LSTIMER('PHIOV MATRIX',tcpu,twall,DECinfo%output)

  end subroutine pair_calculate_mp2density





  !> \brief Save total MP2 density and MP2 correlation densities to files MP2.dens and MP2corr.dens,
  !> respectively. See routine get_full_mp2density for details.
  !> \author Kasper Kristensen
  !> \date September 2011
  subroutine save_mp2density_matrices_to_file(DMP2,DMP2_scaled)


    implicit none

    !> Total MP2 density matrix (including HF contribution) -- NOT changed at output
    type(matrix),intent(inout) :: DMP2
    !> DMP2 scaled to give correct number of electrons
    type(matrix),intent(inout) :: DMP2_scaled
    character(len=80) :: FileName
    integer :: wunit
    logical :: file_exist,OnMaster

    write(DECinfo%output,*) 'Writing full molecular MP2 density to file....'

    ! Save densities to file
    ! **********************

    ! Important note: Before saving to file we multiply by 0.5.
    ! In this way the MP2.dens and MP2corr.dens files use the same convention
    ! as the Hartree-Fock density matrix file "dens.restart", i.e.,
    ! Tr( S dens.restart) = 0.5 * Nelectrons
    ! Tr( S MP2.dens.restart) = 0.5 * Nelectrons
    ! However, at the end we scale back by a factor 2, such that the matrices are unchanged at output
    call mat_scal(0.5E0_realk,DMP2)
    call mat_scal(0.5E0_realk,DMP2_scaled)


    ! MP2 density
    ! -----------
    FileName='MP2.dens'
    inquire(file=FileName,exist=file_exist)
    if(file_exist) then
       write(DECinfo%output,'(a)') 'warning :: overwriting MP2 density file ',FileName
       wunit=-1
       call lsopen(wunit,FileName,'OLD','UNFORMATTED')
       call lsclose(wunit,'DELETE')
    end if

    ! Open file and write data
    wunit=-1
    call lsopen(wunit,FileName,'NEW','UNFORMATTED')
    OnMaster = .TRUE.
    call mat_write_to_disk(wunit,DMP2,OnMaster)
    call lsclose(wunit,'KEEP')



    ! MP2 density scaled
    ! ------------------
    FileName='MP2_scaled.dens'
    inquire(file=FileName,exist=file_exist)
    if(file_exist) then
       write(DECinfo%output,'(a)') 'warning :: overwriting scaled MP2 density file ',FileName
       wunit=-1
       call lsopen(wunit,FileName,'OLD','UNFORMATTED')
       call lsclose(wunit,'DELETE')
    end if

    ! Open file and write data
    wunit=-1
    call lsopen(wunit,FileName,'NEW','UNFORMATTED')
    OnMaster = .TRUE.
    call mat_write_to_disk(wunit,DMP2_scaled,OnMaster)
    call lsclose(wunit,'KEEP')


    ! Scale MP2 densities back such that they are unchanged at output
    call mat_scal(2.0E0_realk,DMP2)
    call mat_scal(2.0E0_realk,DMP2_scaled)



  end subroutine save_mp2density_matrices_to_file




  !> \brief Calculate MP2 density matrices in the AO basis (also initialized here):
  !> 1. Full MP2 density matrix including Hartree-Fock contribution
  !> 2. Correlation density matrix (Dcorr) -- difference between MP2 and HF densities matrices.
  !> \author Kasper Kristensen
  !> \date August 2011
  subroutine get_full_mp2density(MyMolecule,mylsitem,DHF,&
       & Phivo_MO, Phiov_MO, DMP2,rho, Phioo)

    implicit none

    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS setting info
    type(lsitem), intent(inout) :: mylsitem
    !> HF density matrix
    type(matrix),intent(in) :: DHF
    !> Phivo matrix for full molecule (see mp2dens structure)
    type(matrix),intent(in) :: Phivo_MO
    !> Phiov matrix for full molecule (see mp2dens structure)
    type(matrix),intent(in) :: Phiov_MO
    !> MP2 density matrix in AO basis (HF + corr dens), initialized here
    type(matrix),intent(inout) :: DMP2
    !> Relaxed MP2 correlation density matrix in AO basis 
    !> (at input should be unrelaxed density matrix containing only occ-occ and virt-virt blocks,
    !>  while at output it is a relaxed correlation density where all blocks are nonzero)
    type(matrix), intent(inout) :: rho
    !> Phioo matrix required for frozen core (not changed here but intent(inout) for practical reasons)
    type(matrix),intent(inout), optional :: Phioo
    integer :: nvirt,nocc,nbasis,nelectrons,i,ncore,nval
    type(matrix) :: Cocc,Cvirt,S,RHS,kappabarUO,kappabarVC,DMP2_scaled
    real(realk) :: nel_HF, nel_MP2,scale
    logical :: full_equation
    real(realk), dimension(3) :: HFdipole, MP2dipole
    real(realk) :: tcpu,twall, tcpu1,tcpu2, twall1,twall2

    call LSTIMER('START',tcpu,twall,DECinfo%output)
    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    if(MyMolecule%mem_distributed)then
       call lsquit("ERROR(get_full_mp2density) not implemented for distributed matrices in molecule",-1)
    endif

    ! Frozen core sanity check
    if(DECinfo%frozencore) then
       if(.not. present(Phioo)) then
          call lsquit('get_full_mp2density: Phioo matrix required for frozen core!',-1)
       end if
    end if


    ! Easy reference to molecule info
    ! *******************************
    nbasis = MyMolecule%nbasis
    nocc = MyMolecule%nocc
    nvirt = MyMolecule%nvirt
    nelectrons = MyMolecule%nelectrons
    ncore = MyMolecule%ncore
    nval = MyMolecule%nval

    if(DECinfo%frozencore) then
       call mat_init(kappabarVC,nval,ncore)
       call get_kappabar_frozen_core(MyLsitem,MyMolecule,Phioo,kappabarVC)
    end if


    ! Get MO coefficients in matrix form
    ! **********************************
    call mat_init(Cocc,nbasis,nocc)
    call mat_init(Cvirt,nbasis,nvirt)
    call mat_set_from_full(MyMolecule%Co%elm2(1:nbasis,1:nocc), 1E0_realk, Cocc)
    call mat_set_from_full(MyMolecule%Cv%elm2(1:nbasis,1:nvirt), 1E0_realk, Cvirt)


    ! Get RHS matrix for kappabar orbital rotation multiplier equation (dimension nvirt,nocc)
    ! ***************************************************************************************
    UNRELAXED1: if(.not. DECinfo%unrelaxed) then
       if(DECinfo%frozencore) then
          call get_kappabar_RHS(MyLsitem,rho,Phivo_MO,Phiov_MO,Cocc,Cvirt,RHS,kappabarVC=kappabarVC)
       else
          call get_kappabar_RHS(MyLsitem,rho,Phivo_MO,Phiov_MO,Cocc,Cvirt,RHS)
       end if
    end if UNRELAXED1


    ! Solve kappabar (orbital-rotation) multiplier equation
    ! **************************************************
    ! 1. Solve simplified equation without coulomb+exchange transformations
    ! This gives a VERY cheap and reasonable starting guess for the full equation
    full_equation = .false.
    call mat_init(kappabarUO,nvirt,nocc)
    call mat_zero(kappabarUO)
    UNRELAXED2: if(.not. DECinfo%unrelaxed) then
       call dec_solve_kappabar_equation(RHS,full_equation,&
            &MyLsitem,MyMolecule,Cocc,Cvirt,kappabarUO,.false.)

       ! 2. Solve full kappabar multiplier equation using solution from simplified
       !    equation as starting guess.
       full_equation = .true.
       call dec_solve_kappabar_equation(RHS,full_equation,&
            &MyLsitem,MyMolecule,Cocc,Cvirt,kappabarUO,.false.)
       call mat_free(RHS)
    end if UNRELAXED2

    ! Construct correlation density matrix (rho) in AO basis
    ! ******************************************************
    if(DECinfo%frozencore .and. (.not. DECinfo%unrelaxed) ) then
       call dec_get_rho_matrix_in_AO_basis(kappabarUO,MyMolecule,rho,&
            & kappabarVCmat=kappabarVC)
    else
       call dec_get_rho_matrix_in_AO_basis(kappabarUO,MyMolecule,rho)
    end if
    ! Done with kappabar
    call mat_free(Cocc)
    call mat_free(Cvirt)
    call mat_free(kappabarUO)
    !    call mat_print(rho, 1, nbasis, 1,nbasis, DECinfo%output)



    ! Full MP2 density matrix = Hartree-Fock + correlation contribution
    ! *****************************************************************
    ! DMP2 = 2*DHF + rho  (the factor two is due to double occupation -> only for closed-shell)
    call mat_init(DMP2,nbasis,nbasis)
    call mat_add(2E0_realk,DHF,1E0_realk,rho,DMP2)


    ! Sanity checks: The density matrices integrates to the correct number of electrons
    ! *********************************************************************************
    ! Get overlap matrix
    call mat_init(S,nbasis,nbasis)
    call mat_zero(S)
    call II_get_overlap(DECinfo%output,DECinfo%output,mylsitem%setting,S)

    ! Number of electrons calculated for HF state
    Nel_HF= 2E0_realk*mat_dotproduct(DHF,S)

    ! Number of electrons calculated for MP2 state
    Nel_MP2= mat_dotproduct(DMP2,S)

    scale = real(nelectrons)/Nel_MP2
    call mat_init(DMP2_scaled,nbasis,nbasis)
    call mat_copy(1.0_realk,DMP2, DMP2_scaled)
    call mat_scal(scale,DMP2_scaled)

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '*****************************************************************'
    write(DECinfo%output,*) '*                     MP2 DENSITY SANITY CHECK                  *'
    write(DECinfo%output,*) '*****************************************************************'
    write(DECinfo%output,'(1X,a,g20.10)') 'HF state : Number of electrons = ', Nel_HF
    write(DECinfo%output,'(1X,a,g20.10)') 'MP2 state: Number of electrons = ', Nel_MP2
    write(DECinfo%output,'(1X,a,I12)')   'Number of electrons from input =', nelectrons
    write(DECinfo%output,*)
    write(DECinfo%output,*)


    ! Also calculate and print HF and MP2 electric dipole moments
    ! ***********************************************************
    call get_HF_and_MP2_dipole_moments(mylsitem,DHF,DMP2,HFdipole,MP2dipole)
    call print_HF_and_MP2_dipoles(HFdipole, MP2dipole)

    ! Free stuff
    call mat_free(S)

    if(DECinfo%frozencore) then
       call mat_free(kappabarVC)
    end if

    call save_mp2density_matrices_to_file(DMP2,DMP2_scaled)
    call mat_free(DMP2_scaled)

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)
    call LSTIMER('FULL MP2DENS',tcpu,twall,DECinfo%output)

  end subroutine get_full_mp2density




  !> \brief Print HF and MP2 dipole info to LSDALTON.OUT
  !> \author Kasper Kristensen
  !> \date October 2011
  subroutine print_HF_and_MP2_dipoles(HFdipole, MP2dipole)


    implicit none

    !> Hartree-Fock dipole moment
    real(realk), dimension(3),intent(in) :: HFdipole
    !> MP2 dipole moment
    real(realk), dimension(3),intent(in) :: MP2dipole
    real(realk) :: au_to_debye, au_to_SI, normHFdipole,normMP2dipole
    integer :: i

    normHFdipole=0E0_realk
    normMP2dipole=0E0_realk
    do i=1,3
       normHFdipole = normHFdipole + HFdipole(i)**2
       normMP2dipole = normMP2dipole + MP2dipole(i)**2
    end do
    normHFdipole = sqrt(normHFdipole)
    normMP2dipole = sqrt(normMP2dipole)

    au_to_debye=2.54175
    au_to_SI=8.47835
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    if(DECinfo%unrelaxed) then
       write(DECinfo%output,*) '    UNRELAXED DIPOLE MOMENTS FOR HARTREE-FOCK AND MP2 '
       write(DECinfo%output,*) '    ================================================= '
    else
       write(DECinfo%output,*) '    DIPOLE MOMENTS AT THE HARTREE-FOCK AND MP2 LEVELS OF THEORY '
       write(DECinfo%output,*) '    =========================================================== '
    end if
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(6X,A)') '                 HF: Permanent dipole moment'
    write(DECinfo%output,'(6X,A)') '                 ---------------------------'
    write(DECinfo%output,'(12X,A)') '   au              Debye           10**-30 C m'
    write(DECinfo%output,'(5X,3g18.6)') normHFDipole, au_to_debye*normHFDipole, au_to_SI*normHFDipole
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(6X,A)') '                 HF: Dipole moment components'
    write(DECinfo%output,'(6X,A)') '                 ----------------------------'
    write(DECinfo%output,'(12X,A)') '   au              Debye           10**-30 C m'
    write(DECinfo%output,'(1X,A,3X,3g18.6)') 'x', HFdipole(1), &
         & au_to_debye*HFdipole(1), au_to_SI*HFdipole(1)
    write(DECinfo%output,'(1X,A,3X,3g18.6)') 'y', HFdipole(2), &
         & au_to_debye*HFdipole(2), au_to_SI*HFdipole(2)
    write(DECinfo%output,'(1X,A,3X,3g18.6)') 'z', HFdipole(3), &
         & au_to_debye*HFdipole(3), au_to_SI*HFdipole(3)

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)

    write(DECinfo%output,'(6X,A)') '                 MP2: Permanent dipole moment'
    write(DECinfo%output,'(6X,A)') '                 ----------------------------'
    write(DECinfo%output,'(12X,A)') '   au              Debye           10**-30 C m'
    write(DECinfo%output,'(5X,3g18.6)') normMP2Dipole, au_to_debye*normMP2Dipole, au_to_SI*normMP2Dipole
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(6X,A)') '                 MP2: Dipole moment components'
    write(DECinfo%output,'(6X,A)') '                 -----------------------------'
    write(DECinfo%output,'(12X,A)') '   au              Debye           10**-30 C m'
    write(DECinfo%output,'(1X,A,3X,3g18.6)') 'x', MP2dipole(1), &
         & au_to_debye*MP2dipole(1), au_to_SI*MP2dipole(1)
    write(DECinfo%output,'(1X,A,3X,3g18.6)') 'y', MP2dipole(2), &
         & au_to_debye*MP2dipole(2), au_to_SI*MP2dipole(2)
    write(DECinfo%output,'(1X,A,3X,3g18.6)') 'z', MP2dipole(3), &
         & au_to_debye*MP2dipole(3), au_to_SI*MP2dipole(3)

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)

  end subroutine print_HF_and_MP2_dipoles



  !> \brief Construct Theta array required for MP2 density.
  !> Equations are given inside the subroutine.
  !> NOTE: The input EOS amplitudes and the output EOS Theta array are assumed to be kept in memory!
  !> \author Kasper Kristensen
  !> \date September 2011
  subroutine construct_theta_array(t2,Theta)

    implicit none
    !> MP2 amplitudes for EOS
    type(array4) :: t2
    !> Output Theta array (see equation inside subroutine)
    type(array4) :: Theta
    integer :: i,a,j,b
    real(realk) :: tcpu,twall


    ! The Theta array is given by:
    !
    ! Theta_{ij}^{ab} = 4*t2_{ij}^{ab} - 2*t2_{ij}^{ba} + mult_{ij}^{ab}
    !
    ! where the t2 are the MP2 amplitudes and mult are the MP2 multipliers.
    ! The multipliers can simply be calculated from the amplitudes:
    !
    ! mult_{ij}^{ab} = 4*t2_{ij}^{ab} - 2*t2_{ij}^{ba}
    !
    ! Hence, the Theta array is simply given by:
    !
    ! Theta_{ij}^{ab} = 8*t2_{ij}^{ab} - 4*t2_{ij}^{ba}
    !
    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Init Theta array
    Theta = array4_init(t2%dims)


    do j=1,t2%dims(4)
       do b=1,t2%dims(3)
          do i=1,t2%dims(2)
             do a=1,t2%dims(1)
                Theta%val(a,i,b,j) = 8E0_realk*t2%val(a,i,b,j) - 4E0_realk*t2%val(b,i,a,j)
             end do
          end do
       end do
    end do

    call LSTIMER('GET THETA',tcpu,twall,DECinfo%output)

  end subroutine construct_theta_array



  !> \brief Calculate RHS matrix for orbital-rotation multiplier (kappabar):
  !> RHS_{ai} = Phivo_{ai} - Phiov_{ia} + 2G_{ai}(X-Y)
  !> Frozen core: RHS_{ai} = Phivo_{ai} - Phiov_{ia} + 2G_{ai}(X-Y-kappabarVCsym)
  !> The RHS matrix is calculated in the MO basis with dimensions (nvirt,nocc)
  !> See the mp2dens structure for more details.
  !> kappabarVCsym is a symmetrized matrix for valence-core multiplier.
  !> \author Kasper Kristensen
  !> \date November 2010
  subroutine get_kappabar_RHS(MyLsitem,rho_unrelaxed,Phivo,Phiov,Cocc,Cvirt,RHS,kappabarVC)

    implicit none
    !> LSDALTON info
    type(lsitem), intent(inout) :: MyLsitem
    !> Unrelaxed correlation density matrix in AO basis
    !> containing only X and Y blocks, see mp2dens structure
    type(matrix), intent(in) :: rho_unrelaxed
    !> Phivo matrix in MO basis, see mp2dens structure
    type(matrix), intent(in) :: Phivo
    !> Phiov matrix in MO basis, see mp2dens structure
    type(matrix), intent(in) :: Phiov
    !> Occupied MOs in type matrix form (not changed though intent(inout))
    type(matrix),intent(inout) :: Cocc
    !> Virtual MOs in type matrix form (not changed though intent(inout))
    type(matrix),intent(inout) :: Cvirt
    !> RHS matrix for kappabar equation
    type(matrix), intent(inout) :: RHS
    !> Kappabar frozen core multipliers for valence,core block
    type(matrix),intent(in),optional :: kappabarVC
    type(matrix) :: M,FockM_AO,FockM, Phiov_trans, VCAOsym,VCAO, VCAOT, Ccore,Cval
    integer :: nocc,nvirt,nbasis,ncore,nval,i
    real(realk) :: scaling_factor
    real(realk),pointer :: tmp(:,:)
    real(realk) :: tcpu,twall

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Sanity check
    if(DECinfo%frozencore) then
       if(.not. present(kappabarVC)) then
          call lsquit('get_kappabar_RHS: Frozen core but core-valence multipliers not present!',-1)
       end if
    end if

    ! Initialize stuff
    ! ****************
    nocc  = Cocc%ncol     ! Number of occupied orbitals
    nvirt = Cvirt%ncol    ! Number of virtual orbitals
    nbasis = Cocc%nrow    ! Number of basis functions


    ! Calculate matrix M = -2*rho_unrelaxed
    ! *************************************
    ! This is the matrix M in Eq. (A16) in JCP 137, 114102 (2012)
    call mat_init(M,nbasis,nbasis)
    call mat_assign(M,rho_unrelaxed)
    scaling_factor = -2.0_realk
    call mat_scal(scaling_factor,M)


    ! For frozen core, add symmetrized kappabar multiplier for core-valence splitting
    if(DECinfo%frozencore) then
       ncore = kappabarVC%ncol  ! #core orbitals
       nval = kappabarVC%nrow   ! #valence orbitals

       ! Get core MO coefficients
       call mem_alloc(tmp,nbasis,ncore)
       call mat_retrieve_block(Cocc,tmp,nbasis,ncore,1,1)
       call mat_init(Ccore,nbasis,ncore)
       call mat_set_from_full(tmp,1.0_realk,Ccore)
       call mem_dealloc(tmp)

       ! Get valence MO coefficients
       call mem_alloc(tmp,nbasis,nval)
       call mat_retrieve_block(Cocc,tmp,nbasis,nval,1,ncore+1)
       call mat_init(Cval,nbasis,nval)
       call mat_set_from_full(tmp,1.0_realk,Cval)
       call mem_dealloc(tmp)

       ! Transform kappabar(valence,core) to AO basis
       call mat_init(VCAO,nbasis,nbasis)
       call util_MO_to_AO_different_trans(Cval,kappabarVC,Ccore,VCAO)
       call mat_free(Cval)
       call mat_free(Ccore)

       ! Symmtrize kappabar: KappabarSYM_{mu,nu} = kappabar(mu,nu) + kappabar(nu,mu)
       call mat_init(VCAOT,nbasis,nbasis)
       call mat_init(VCAOsym,nbasis,nbasis)
       call mat_trans(VCAO,VCAOT)
       call mat_add(1E0_realk, VCAO, 1E0_realk, VCAOT, VCAOsym)

       ! Update M --> M = 2*X - 2*Y - 2*kappabarVCsym
       call mat_daxpy(-2.0_realk, VCAOsym, M)
       call mat_free(VCAOsym)
       call mat_free(VCAOT)
       call mat_free(VCAO)

    end if


    ! Fock transformation (coulomb+exchange) on symmetric M matrix (FockM_AO)
    ! ***********************************************************************

    call mat_init(FockM_AO,nbasis,nbasis)
    call util_get_symm_part(M)
    call dec_fock_transformation(FockM_AO,M,MyLsitem,.true.)
    call mat_free(M)



    ! Get Fock transformed matrix in the mixed MO basis (virt,occ): FockM_MO = Cvirt^T FockM_AO Cocc
    ! **********************************************************************************************

    call mat_init(FockM,nvirt,nocc)
    call util_AO_to_MO_different_trans(Cvirt,FockM_AO,Cocc,FockM)
    call mat_free(FockM_AO)



    ! Transpose Phiov matrix to get same dimensions (nvirt,nocc) as the other matrices
    ! ********************************************************************************
    call mat_init(Phiov_trans,nvirt,nocc)
    call mat_trans(Phiov,Phiov_trans)



    ! Calculate RHS matrix
    ! ********************

    ! Now we have the following in the MO basis of dimension (nvirt,nocc):
    !
    ! Phivo_{ai}   (from input)
    ! Phiov_trans_{ai} = Phiov_{ia}
    ! FockM_{ai} = G_{ai}(M)          [G transformation corresponds to both coulomb and exchange]
    !
    ! Thus, we can now calculate the RHS matrix (dimension: nvirt,nocc):
    ! RHS_{ai} =  Phivo_{ai} - Phiov_{ia} + G_{ai}(M)

    call mat_init(RHS,nvirt,nocc)

    ! 1. RHS_{ai} = Phivo_{ai} - Phiov_{ia}
    call mat_add(1E0_realk,Phivo,-1E0_realk,Phiov_trans,RHS)

    ! 2. RHS_{ai} = Phivo_{ai} - Phiov_{ia} + G_{ai}(M)
    call mat_daxpy(1E0_realk,FockM,RHS)


    ! Free remaining matrices
    call mat_free(FockM)
    call mat_free(Phiov_trans)

    call LSTIMER('KAPPABAR RHS',tcpu,twall,DECinfo%output)


  end subroutine get_kappabar_RHS



  !> \brief Solve kappabar (orbital-rotation) multiplier equation.
  !> Note:
  !> 1. The kappabar_final input matrix is used as starting guess.
  !> 2. The coulomb+exchange terms are omitted if full_equation=.false.
  !> 3. For frozen core we have full_equation=.false.
  !>    and the virtual/occupied Fock matrix transformations in the
  !>    virt-occ kappabar multiplier equation is replaced by valence/core transformations.
  !> \author Kasper Kritensen
  subroutine dec_solve_kappabar_equation(RHS,full_equation,&
       &MyLsitem,MyMolecule,Cocc,Cvirt,kappabar_final,fc)


    implicit none
    !> Right-hand side matrix for kappabar multiplier equation
    type(matrix), intent(in) :: RHS
    !> Solve full kappabar equation (true) or reduced equation (false).
    logical, intent(in) :: full_equation
    !> LSDALTON info
    type(lsitem), intent(inout) :: MyLsitem
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Occupied MO coefficients
    type(matrix), intent(in) :: Cocc
    !> Virtual MO coefficients
    type(matrix), intent(in) :: Cvirt
    !> Final kappabar solution vector. The kappabar_final input is used as starting guess.
    type(matrix), intent(inout) :: kappabar_final
    !> Solve frozen core equation rather than (virt,occ) kappabar equation
    logical,intent(in) :: fc
    integer :: nocc, nvirt,ncore,nval
    type(matrix), pointer :: kappabar(:), residual(:)
    type(matrix) :: residual_opt, residual_prec,kappabar_opt
    type(matrix) :: ppfock, qqfock,prec
    reaL(realk) :: tcpu,twall
    integer, dimension(2) :: kappabar_dims
    real(realk), pointer :: B(:,:),c(:), Fockcorecore(:,:), Fockvalval(:,:),precfull(:,:)
    logical :: crop_ok,break_iterations
    real(realk) :: prev_norm, one_norm_total, two_norm_total
    real(realk) :: convthr
    integer :: iter, last_iter, i,j,dim1,dim2,a

    call LSTIMER('START',tcpu,twall,DECinfo%output)
    if(MyMolecule%mem_distributed)then
       call lsquit("ERROR(dec_solve_kappabar_equation) not implemented for distributed matrices in molecule",-1)
    endif


    ! Dimensions for RHS (and thus for kappabar)
    dim1 = kappabar_final%nrow
    dim2 = kappabar_final%ncol

    ! Sanity checks
    ! *************
    if( (dim1 /= RHS%nrow) .or. (dim2/=RHS%ncol) ) then
       write(DECinfo%output,*) 'Kappabar dims', dim1,dim2
       write(DECinfo%output,*) 'RHS   dims', RHS%nrow, RHS%ncol
       call lsquit('dec_solve_kappabar_equation: Kappabar and RHS dimension mismatch',-1)
    end if

    if(full_equation .and. fc) then
       call lsquit('dec_solve_kappabar_equation: Not implemented for full equation AND frozen core!',-1)
    end if

    if(DECinfo%ccModel /= MODEL_MP2 .and. DECinfo%ccModel /= MODEL_RIMP2) then
       print *, 'CC model: ', DECinfo%ccModel
       call lsquit('dec_solve_kappabar_equation:&
            & kappabar multiplier equation is only implemented for RIMP2 and MP2',-1)
    end if

    write(DECinfo%output,*)
    if(full_equation) then
       write(DECinfo%output,*) 'Starting kappabar multiplier solver for full equation...'
    else
       if(fc) then
          write(DECinfo%output,*) 'Starting frozen core kappabar multiplier solver for zeroth order equation...'
       else
          write(DECinfo%output,*) 'Starting kappabar multiplier solver for zeroth order equation...'
       end if
    endif
    write(DECinfo%output,*)


    ! Initialize stuff
    ! ****************
    nocc=MyMolecule%nocc
    nvirt=MyMolecule%nvirt
    ncore = MyMolecule%ncore
    nval = MyMolecule%nval

    ! The convergence threshold
    convthr = DECinfo%kappaTHR
    if(fc) then
       write(DECinfo%output,'(a,g16.6)') 'Frozen core kappabar multiplier equation: Conv. thr. set to', convthr
    else
       write(DECinfo%output,'(a,g16.6)') 'Virt-occ kappabar multiplier equation: Conv. thr. set to', convthr
    end if

    ! CROP/DIIS matrices
    call mem_alloc(B,DECinfo%kappaMaxIter,DECinfo%kappaMaxIter)
    call mem_alloc(c,DECinfo%kappaMaxIter)

    ! Kappabar vectors and residuals
    call mem_alloc(kappabar,DECinfo%kappaMaxIter)
    call mem_alloc(residual,DECinfo%kappaMaxIter)

    ! iterate
    break_iterations = .false.
    crop_ok = .false.
    prev_norm = 1.0E6_realk


    ! Transformation matrices
    ! ***********************

    ! For kappabar(virt,occ) and kappabar(valence,core) we need equivalent transformation matrices.
    ! For kappabar(virt,occ) we need virt-virt and occ-occ Fock matrix blocks.
    ! For kappabar(valence,core) we simple replace virt-virt and occ-occ Fock matrix
    ! blocks by valence,valence and core-core blocks, respectively.
    ! We therefore simply set ppfock and qqfock differently according to the fc
    ! input and the same code can then be used for kappabar(virt,occ) and kappabar(valence,core).
   
    if(fc) then  ! kappabar(valence,core)
       
       ! Get core-core and val-val blocks
       call mem_alloc(Fockcorecore,ncore,ncore)
       call mem_alloc(Fockvalval,nval,nval)
       do j=1,ncore
          do i=1,ncore
             Fockcorecore(i,j) = MyMolecule%oofock%elm2(i,j)
          end do
       end do
       do j=1,nval
          do i=1,nval
             Fockvalval(i,j) = MyMolecule%oofock%elm2(i+ncore,j+ncore)
          end do
       end do

       ! Preconditioning matrix
       call mem_alloc(precfull,nval,ncore)
       do j=1,ncore
          do i=1,nval
             precfull(i,j) = Fockvalval(i,i) - Fockcorecore(j,j)
          end do
       end do

       ! Convert to matrix format and store in ppfock and qqfock
       call mat_init(ppfock,ncore,ncore)
       call mat_set_from_full(Fockcorecore, 1E0_realk,ppfock)
       call mat_init(qqfock,nval,nval)
       call mat_set_from_full(Fockvalval, 1E0_realk,qqfock)
       call mem_dealloc(Fockcorecore)
       call mem_dealloc(Fockvalval)
       call mat_init(prec,nval,ncore)
       call mat_set_from_full(precfull, 1.0_realk, prec)

    else

       ! Occ-occ Fock block (both core and valence)
       call mat_init(ppfock,nocc,nocc)
       call mat_set_from_full(MyMolecule%oofock%elm2(1:nocc,1:nocc), 1E0_realk,ppfock)

       ! Virt-virt Fock matrix block
       call mat_init(qqfock,nvirt,nvirt)
       call mat_set_from_full(MyMolecule%vvfock%elm2(1:nvirt,1:nvirt), 1E0_realk,qqfock)

       ! Preconditioning matrix
       call mem_alloc(precfull,nvirt,nocc)
       do i=1,nocc
          do a=1,nvirt
             precfull(a,i) = MyMolecule%vvfock%elm2(a,a) - MyMolecule%oofock%elm2(i,i)
          end do
       end do
       call mat_init(prec,nvirt,nocc)
       call mat_set_from_full(precfull, 1.0_realk, prec)

    end if
    call mem_dealloc(precfull)


    write(DECinfo%output,*)
    write(DECinfo%output,*) 'Starting kappabar iterations'
    write(DECinfo%output,*) '-------------------------'
    write(DECinfo%output,'(1X,a)')  '***  Iteration     Residual norm'


    KappabarIteration : do iter=1,DECinfo%kappaMaxIter


       ! remove old vectors
       RemoveOldVectors : if(iter > DECinfo%kappaMaxDIIS) then
          call mat_free(kappabar(iter-DECinfo%kappaMaxDIIS))
          call mat_free(residual(iter-DECinfo%kappaMaxDIIS))
       end if RemoveOldVectors

       ! get new vectors
       GetGuessVectors : if(iter == 1) then
          call mat_init(kappabar(iter),dim1,dim2)
          ! Use input kappabar (might be zero) as starting guess
          call mat_assign(kappabar(iter),kappabar_final)
       end if GetGuessVectors

       call mat_init(residual(iter),dim1,dim2)
       ! Get current residual
       if(full_equation) then ! Full equation to be solved
          call dec_get_kappabar_residual(kappabar(iter), RHS, &
               & ppfock, qqfock, Cocc, Cvirt, MyLsitem,residual(iter))
       else ! zeroth order equation to be solved
          call get_zeroth_order_kappabar_residual(kappabar(iter), RHS, &
               & ppfock, qqfock,residual(iter))
       endif

       KappabarDebug : if(DECinfo%kappa_driver_debug) then

          write(DECinfo%output,*) 'Kappabar iteration number', iter
          write(DECinfo%output,'(a,f16.10)') ' debug :: kappabar norm    ',&
               & mat_sqnorm2(kappabar(iter))
          write(DECinfo%output,'(a,f16.10)') ' debug :: residual norm    ',&
               & mat_sqnorm2(residual(iter))
          write(DECinfo%output,'(a,f16.10)') ' debug :: ppfock norm    ',&
               & mat_sqnorm2(ppfock)
          write(DECinfo%output,'(a,f16.10)') ' debug :: qqfock norm    ',&
               & mat_sqnorm2(qqfock)
          write(DECinfo%output,'(a,f16.10)') ' debug :: RHS norm    ',&
               & mat_sqnorm2(RHS)

       end if KappabarDebug

       ! calculate crop/diis matrix
       B=0.0E0_realk; c=0.0E0_realk
       do i=iter,max(iter-DECinfo%kappaMaxDIIS+1,1),-1
          do j=iter,i,-1
             if(DECinfo%kappa_use_preconditioner_in_b) then
                call mat_init(residual_prec,dim1,dim2)
                call dec_precondition_kappabar(residual(j),prec,residual_prec)
                B(i,j) = mat_dotproduct(residual(i),residual_prec)
                call mat_free(residual_prec)
             else
                B(i,j) = mat_dotproduct(residual(i),residual(j))
             end if
             B(j,i) = B(i,j)
          end do
       end do


       ! solve crop/diis equation
       call CalculateDIIScoefficients(DECinfo%kappaMaxDIIS,DECinfo%kappaMaxIter,iter,B,c, &
            DECinfo%kappa_driver_debug)

       ! mixing to get optimal
       call mat_init(kappabar_opt,dim1,dim2)
       call mat_init(residual_opt,dim1,dim2)
       call mat_zero(kappabar_opt)
       call mat_zero(residual_opt)
       do i=iter,max(iter-DECinfo%kappaMaxDIIS+1,1),-1
          call mat_daxpy(c(i),kappabar(i),kappabar_opt)
          call mat_daxpy(c(i),residual(i),residual_opt)
       end do

       ! if crop, put the optimal in place of trial (not for diis)
       if(DECinfo%use_crop) then
          call mat_assign(kappabar(iter),kappabar_opt)
          call mat_assign(residual(iter),residual_opt)
       end if

       ! check for convergence
       one_norm_total = mat_sqnorm2(residual(iter))  ! norm squared
       two_norm_total = sqrt(one_norm_total)   ! norm itself

       write(DECinfo%output,'(1X,a,2X,i4,5X,g18.8)')  '*** ',iter, two_norm_total

       ! simple crop diagnostics
       if(two_norm_total < prev_norm) then
          crop_ok=.true.
       else
          crop_ok=.false.
          write(DECinfo%output,'(a)') ' warning :: total norm was smaller in previous iteration !!! '
       end if
       prev_norm=two_norm_total

       ! check if this is the last iteration
       if( (iter == DECinfo%kappaMaxIter) .or. (two_norm_total < convthr) ) then
          break_iterations=.true.
       end if

       ! Save final kappabar
       if(break_iterations) then
         call mat_assign(kappabar_final,kappabar_opt)
       endif

       ! generate next trial vector only if this is not the last iteration
       if(.not.break_iterations) then
          if(DECinfo%kappa_use_preconditioner) then
             call mat_init(residual_prec,dim1,dim2)
             call dec_precondition_kappabar(residual_opt,prec,residual_prec)
             ! set kappabar(n+1) = kappabar_opt(n) + C^{-1}*residual(n)
             ! (Eq. 5.50 in Marcin's thesis)
             call mat_init(kappabar(iter+1),dim1,dim2)
             call mat_add(1E0_realk,kappabar_opt,1E0_realk,residual_prec,kappabar(iter+1))
             call mat_free(residual_prec)
          else
             call mat_init(kappabar(iter+1),dim1,dim2)
             call mat_add(1E0_realk,kappabar_opt,1E0_realk,residual_opt,kappabar(iter+1))
          end if
       end if

       ! delete optimals
       call mat_free(kappabar_opt)
       call mat_free(residual_opt)

       last_iter = iter
       if(break_iterations) exit

    end do KappabarIteration



    if(break_iterations) then
       if(full_equation) then
          write(DECinfo%output,'(A41,I5,A12)') &
               & 'Full kappabar multiplier equation solved in', last_iter, ' iterations!'
       else
          write(DECinfo%output,'(A49,I5,A12)') &
               & 'Zeroth order kappabar multiplier equation solved in ',&
               last_iter, ' iterations!'
       endif
       write(DECinfo%output,*) 'Residual norm:', two_norm_total
    endif


    ! deallocate stuff
    call mem_dealloc(B)
    call mem_dealloc(c)

    ! remove rest of the amplitudes and residuals
    do i=last_iter,max(last_iter-DECinfo%kappaMaxDIIS+1,1),-1
       call mat_free(kappabar(i))
       call mat_free(residual(i))
    end do
    call mat_free(ppfock)
    call mat_free(qqfock)
    call mat_free(prec)

    call mem_dealloc(kappabar)
    call mem_dealloc(residual)

    if(fc) then
       call LSTIMER('FC KAPPA EQ.',tcpu,twall,DECinfo%output)
    else
       if(full_equation) then
          call LSTIMER('FULL KAPPA EQ.',tcpu,twall,DECinfo%output)
       else
          call LSTIMER('RED. KAPPA EQ.',tcpu,twall,DECinfo%output)
       end if
    end if
 
  end subroutine dec_solve_kappabar_equation




  !> \brief Preconditioning for kappabar equation (dimension nvirt x nocc) using
  !> prec_{a,i} = residual(a,i) / [Fock(a,a) - Fock(i,i)]
  !> (Always works for frozen core equation, then virtual->valence and occupied->core).
  !> All matrices are expressed in the MO basis.
  !> \author Kasper Kristensen
  !> \date November 2010
  subroutine dec_precondition_kappabar(residual,prec,residual_prec)

    implicit none
    !> Original residual
    type(matrix), intent(in) :: residual
    !> Preconditioner: prec(a,i) = Fock(a,a) - Fock(i,i)
    type(matrix),intent(in) :: prec
    !> Preconditioned residual
    type(matrix), intent(inout) :: residual_prec
    real(realk) :: mu

    mu=0.0_realk
    ! Copy residual into output vector
    call mat_assign(residual_prec,residual)
    ! "Divide" by preconditioning matrix
    call mat_hdiv(residual_prec,prec,mu)

  end subroutine dec_precondition_kappabar




  !> \brief Get residual for kappabar multiplier equation,
  !> Residual = RHS - [E2 kappabar]
  !> where
  !> [E2 kappabar] = E2_fock kappabar + 4*[2*J(kappabar) - K(kappabar)]
  !> where E2fock kappabar is the part of the E2 transformation
  !> converning the Fock matrix,
  !> and J(kappabar) and K(kappabar) are Coulomb and exchange transformations.
  !> (see subroutines dec_get_E2fock_kappabar and
  !>  dec_get_coulomb_exchange_on_kappabar for details).
  !> \author Kasper Kristensen
  !> \date November 2010
  subroutine dec_get_kappabar_residual(kappabar, RHS, &
       & ppfock, qqfock, Cocc, Cvirt, MyLsitem,residual)

    !> Current solution matrix to kappabar multiplier equation (MO basis)
    type(matrix), intent(in) :: kappabar
    !> RHS matrix for kappabar multiplier equation (MO basis)
    type(matrix), intent(in) :: RHS
    !> Occ-occ block of Fock matrix in MO basis
    type(matrix), intent(in) :: ppfock
    !> Virt-virt block of Fock matrix in MO basis
    type(matrix), intent(in) :: qqfock
    !> Occupied MO coefficients
    type(matrix), intent(in) :: Cocc
    !> Virtual MO coefficients
    type(matrix), intent(in) :: Cvirt
    !> LSDALTON info
    type(lsitem), intent(inout) :: MyLsitem
    !> Residual kappabar multiplier equation
    type(matrix),intent(inout) :: residual
    type(matrix) :: E2_kappabar
    integer :: nocc,nvirt


    ! Initialize
    ! **********
    nvirt=kappabar%nrow
    nocc=kappabar%ncol

    ! Calculate E2 transformation on kappabar (only Fock matrix contribution)
    ! ********************************************************************
    call mat_init(E2_kappabar,nvirt,nocc)
    call dec_get_E2_kappabar(kappabar,ppfock,qqfock,Cocc,Cvirt,MyLsitem,E2_kappabar)


    ! Subtract from RHS matrix to get residual
    ! ****************************************
    call mat_add(1E0_realk,RHS,-1E0_realk,E2_kappabar,residual)

    ! Free stuff
    ! **********
    call mat_free(E2_kappabar)

  end subroutine dec_get_kappabar_residual




  !> \brief Get linear transformation for electronic Hessian on kappabar vector:
  !> [E2 kappabar] = E2_fock kappabar + 4*[2*J(kappabar) - K(kappabar)]
  !> where E2fock kappabar is the part of the E2 transformation
  !> converning the Fock matrix,
  !> and J(kappabar) and K(kappabar) are Coulomb and exchange transformations.
  !> (see subroutines dec_get_E2fock_kappabar and
  !>  dec_get_coulomb_exchange_on_kappabar for details).
  !> \author Kasper Kristensen
  !> \date November 2010
  subroutine dec_get_E2_kappabar(kappabar, ppfock, qqfock, &
       &Cocc, Cvirt, Mylsitem, E2_kappabar)

    implicit none
    !> Current solution matrix to kappabar multiplier equation (MO basis)
    type(matrix), intent(in) :: kappabar
    !> Occ-occ block of Fock matrix in MO basis
    type(matrix), intent(in) :: ppfock
    !> Virt-virt block of Fock matrix in MO basis
    type(matrix), intent(in) :: qqfock
    !> Occupied MO coefficients
    type(matrix), intent(in) :: Cocc
    !> Virtual MO coefficients
    type(matrix), intent(in) :: Cvirt
    !> LSDALTON info
    type(lsitem), intent(inout) :: MyLsitem
    !> E2-transformed kappabar matrix
    type(matrix), intent(inout) :: E2_kappabar
    type(matrix) :: E2fock_kappabar, Gkappabar
    integer :: nvirt, nocc


    ! Initialize
    ! **********
    nvirt=kappabar%nrow
    nocc=kappabar%ncol

    ! Get Fock matrix transformation of kappabar matrix
    ! **********************************************
    call mat_init(E2fock_kappabar,nvirt,nocc)
    ! (E2_fock kappabar)_{ai} = 2 * [ sum_{b} F_{ab} kappabar_{bi} - sum_{j} kappabar_{aj} F_{ji} ]
    call dec_get_E2fock_kappabar(kappabar, ppfock, qqfock, E2fock_kappabar)


    ! Get coulomb+exchange (G) transformations on kappabar matrix
    ! ********************************************************
    call mat_init(Gkappabar,nvirt,nocc)
    call dec_get_coulomb_exchange_on_kappabar(kappabar,ppfock,qqfock,Cocc,Cvirt,&
         &MyLsitem,Gkappabar)


    ! Total E2 transformation: E2_kappabar = E2fock_kappabar + Gkappabar
    ! *********************************************************
    call mat_add(1E0_realk, E2fock_kappabar, 1E0_realk, Gkappabar, E2_kappabar)
    ! Free matrices
    ! *************
    call mat_free(E2fock_kappabar)
    call mat_free(Gkappabar)


  end subroutine dec_get_E2_kappabar




  !> \brief Get linear transformation for dominating (Fock) terms in electronic Hessian:
  !> (E[2]_fock kappabar)_{ai} = 2 * [ sum_{b} F_{ab} kappabar_{bi} - sum_{j} kappabar_{aj} F_{ji} ]
  !> \author Kasper Kristensen
  !> \date November 2010
  subroutine dec_get_E2fock_kappabar(kappabar, ppfock, qqfock, E2fock_kappabar)

    implicit none
    !> E2-transformed kappabar matrix
    type(matrix), intent(inout) :: E2fock_kappabar
    !> Current solution matrix to kappabar multiplier equation (MO basis)
    type(matrix), intent(in) :: kappabar
    !> Occ-occ block of Fock matrix in MO basis
    type(matrix), intent(in) :: ppfock
    !> Virt-virt block of Fock matrix in MO basis
    type(matrix), intent(in) :: qqfock
    type(matrix) :: tmp1,tmp2
    integer :: nvirt, nocc


    ! Initialize
    ! **********
    nvirt=kappabar%nrow
    nocc=kappabar%ncol


    ! tmp1_{ai} = 2*sum_{b} F_{ab} kappabar_{bi}
    ! *************************************
    call mat_init(tmp1,nvirt,nocc)
    call mat_mul(qqfock,kappabar,'n','n',2E0_realk,0E0_realk,tmp1)


    ! tmp2_{ai} = - 2*sum_{j} kappabar_{aj} F_{ji}
    ! *************************************
    call mat_init(tmp2,nvirt,nocc)
    call mat_mul(kappabar,ppfock,'n','n',-2E0_realk,0E0_realk,tmp2)


    ! E[2]_fock kappabar = tmp1 + tmp2
    ! *********************************
    call mat_add(1E0_realk,tmp1,1E0_realk,tmp2,E2fock_kappabar)


    ! Free matrices
    ! *************
    call mat_free(tmp1)
    call mat_free(tmp2)


  end subroutine dec_get_E2fock_kappabar


  !> \brief Calculate coulomb and exchange transformations on kappabar matrix.
  !> MO basis:
  !> G_MO(kappabar) = Cvirt^T G_AO(kappabar) Cocc
  !> G_AO(kappabar) = 4*J(kappabar_sym) - 2*K(kappabar_sym)
  !> where kappabar_sym is a symmetrized kappabar matrix in the AO basis
  !> kappabar_sym = kappabar_AO + kappabar_AO^T
  !> and the Coulomb and exchange transformations in the AO basis are
  !> J_{mu nu}(X) = sum_{alpha beta} (mu nu | alpha beta) X_{alpha beta}
  !> K_{mu nu}(X) = sum_{alpha beta} (mu alpha | beta nu) X_{alpha beta}
  !> \author Kasper Kristensen
  !> \date November 2010
  subroutine dec_get_coulomb_exchange_on_kappabar(kappabar,ppfock,qqfock,Cocc,Cvirt,&
       &MyLsitem,Gkappabar_MO)


    implicit none
    !> Coulomb+exchange on kappabar: G(kappabar) = 4*J(kappabar_sym) - 2*K(kappabar_sym)
    type(matrix), intent(inout) :: Gkappabar_MO
    !> Kappabar multiplier matrix
    type(matrix),intent(in) :: kappabar
    !> Occ-occ block of Fock matrix (MO basis)
    type(matrix),intent(in) :: ppfock
    !> Virt-virt block of Fock matrix (MO basis)
    type(matrix),intent(in) :: qqfock
    !> Occupied MO coeffiecients
    type(matrix),intent(in) :: Cocc
    !> Virtual MO coeffiecients
    type(matrix),intent(in) :: Cvirt
    !> LSDALTON info
    type(lsitem), intent(inout) :: MyLsitem
    type(matrix) :: GkappabarAO, kappabar_AO,kappabar_sym,TMP
    integer :: nbasis



    ! Number of basis functions
    nbasis = Cocc%nrow



    ! -------------------------------------------------------------------------------
    !                            Transform kappabar to AO basis                        !
    ! -------------------------------------------------------------------------------
    ! kappabar_AO = Cvirt kappabar_MO Cocc^T
    call mat_init(kappabar_AO,nbasis,nbasis)
    call util_MO_to_AO_different_trans(Cvirt,kappabar,Cocc,kappabar_AO)



    ! -------------------------------------------------------------------------------
    !                Construct symmetrized kappabar_AO matrix kappabar_sym                !
    ! -------------------------------------------------------------------------------

    call mat_init(TMP,nbasis,nbasis)
    call mat_trans(kappabar_AO,TMP)

    call mat_init(kappabar_sym,nbasis,nbasis)

    ! kappabar_sym = kappabar_AO + kappabar_AO^T
    call mat_add(1E0_realk, kappabar_AO, 1E0_realk, TMP, kappabar_sym)

    ! Done with kappabar_AO and TMP
    call mat_free(kappabar_AO)
    call mat_free(TMP)

    ! -------------------------------------------------------------------------------
    !                  Carry out Fock transformation on kappabar_sym                   !
    ! -------------------------------------------------------------------------------

    ! GkappabarAO = 2*J(kappabar_sym) - K(kappabar_sym)
    call mat_init(GkappabarAO,nbasis,nbasis)
    call util_get_symm_part(kappabar_sym)
    call dec_fock_transformation(GkappabarAO,kappabar_sym,MyLsitem,.true.)

    ! Done with kappabar_sym
    call mat_free(kappabar_sym)

    ! The total coulomb and exchange contributions are (p. 96 in notes)
    ! 4 Coulomb(kappabar_sym) - 2 Exchange(kappabar_sym).
    ! Full_fock_transformation calculates
    ! 2 Coulomb(kappabar_sym) - 1 Exchange(kappabar_sym)
    ! and we therefore multiply GkappabarAO by two.
    call mat_scal(2E0_realk,GkappabarAO)



    ! -------------------------------------------------------------------------------
    !                              Convert AO to MO                                 !
    ! -------------------------------------------------------------------------------

    ! The Gkappabar_MO (virtupied,occupied) in the MO basis can be found
    ! from a similairty transformation (notes p. 96)
    ! Gkappabar_MO = Cvirt^T * Gkappabar_AO * Cocc
    call util_AO_to_MO_different_trans(Cvirt,GkappabarAO,Cocc,Gkappabar_MO)

    ! Done with GkappabarAO
    call mat_free(GkappabarAO)


  end subroutine dec_get_coulomb_exchange_on_kappabar



  !> \brief Get zeroth order residual for kappabar multiplier equation,
  !> meaning the we omit the computationally expensive coulomb+exchange part
  !> of the E2 transformation,
  !> Residual = RHS - (E[2]_fock kappabar)
  !> where
  !> (E[2]_fock * kappabar)_{ai} = 2 * [ sum_{b} F_{ab} kappabar_{bi} - sum_{j} F_{ij} kappabar_{aj} ]
  !> \author Kasper Kristensen
  !> \date November 2010
  subroutine get_zeroth_order_kappabar_residual(kappabar, RHS, &
       & ppfock, qqfock,residual)

    !> Residual kappabar multiplier equation
    type(matrix),intent(inout) :: residual
    !> Current solution matrix to kappabar multiplier equation (MO basis)
    type(matrix), intent(in) :: kappabar
    !> RHS matrix for kappabar multiplier equation (MO basis)
    type(matrix), intent(in) :: RHS
    !> Occ-occ block of Fock matrix in MO basis
    type(matrix), intent(in) :: ppfock
    !> Virt-virt block of Fock matrix in MO basis
    type(matrix), intent(in) :: qqfock
    type(matrix) :: E2fock_kappabar
    integer :: nocc,nvirt


    ! Initialize
    ! **********
    nvirt=kappabar%nrow
    nocc=kappabar%ncol

    ! Calculate E2 transformation on kappabar (only Fock matrix contribution)
    ! ********************************************************************
    call mat_init(E2fock_kappabar,nvirt,nocc)
    call dec_get_E2fock_kappabar(kappabar,ppfock,qqfock,E2fock_kappabar)


    ! Subtract from RHS matrix to get residual
    ! ****************************************
    call mat_add(1E0_realk,RHS,-1E0_realk,E2fock_kappabar,residual)


    ! Free stuff
    ! **********
    call mat_free(E2fock_kappabar)

  end subroutine get_zeroth_order_kappabar_residual




  !> \brief Construct MP2 correlation density matrix (rho) in AO basis.
  !> See expression for rho inside subroutine.
  !> \author Kasper Kristensen
  !> \date November 2010
  subroutine dec_get_rho_matrix_in_AO_basis(kappabarmat,MyMolecule,rhoAOmat,kappabarVCmat)

    implicit none
    !> Orbital-rotation multiplers
    type(matrix), intent(in) :: kappabarmat
    !> Full molecule structure
    type(fullmolecule),intent(in) :: MyMolecule
    !> Correlation density matrix (rho) in AO basis:
    !> At input: Unrelaxed correlation density with only occ-occ and virt-virt blocks set (see mp2dens type def.)
    !> At output: Relaxed correlation density with also occ-virt and virt-occ blocks set.
    type(matrix), intent(inout) :: rhoAOmat
    !> Kappabar frozen core multipliers for valence,core block in MO basis
    type(matrix),intent(in),optional :: kappabarVCmat
    type(matrix) :: rhoAOmat_relaxation
    integer :: nocc, nvirt, nbasis,nval,ncore
    integer :: i,j,a,b,ax,bx,start,ix
    real(realk) :: tcpu,twall
    real(realk),pointer :: kappabar(:,:),kappabarVC(:,:),rhoMO_relaxation(:,:),rhoAO_relaxation(:,:),C(:,:),tmp(:,:)

    if(DECinfo%unrelaxed) then
       ! Skipping calculation of relaxation part!
       return
    end if

    if(MyMolecule%mem_distributed)then
       call lsquit("ERROR(dec_get_rho_matrix_in_AO_basis) not implemented for distributed matrices in molecule",-1)
    endif

    call LSTIMER('START',tcpu,twall,DECinfo%output)


    if(DECinfo%frozencore) then
       if(.not. present(kappabarVCmat)) then
          call lsquit('dec_get_rho_matrix_in_AO_basis: Valence-core multipliers not present!',-1)
       end if
    end if


    ! Initialize stuff
    ! ****************
    nocc = MyMolecule%nocc
    nvirt = MyMolecule%nvirt
    nbasis = MyMolecule%nbasis
    ncore = MyMolecule%ncore
    nval = MyMolecule%nval
    if(DECinfo%frozencore) then
       start = ncore+1
    else
       start=1
    end if


    ! Convert components of rho matrix to full arrays
    ! ***********************************************
    ! I know this is not pretty or very efficient but it is negligible work in the big picture...
    call mem_alloc(kappabar,kappabarmat%nrow,kappabarmat%ncol)
    call mat_to_full(kappabarmat,1.0_realk,kappabar)
    if(DECinfo%frozencore) then
       call mem_alloc(kappabarVC,kappabarVCmat%nrow,kappabarVCmat%ncol)
       call mat_to_full(kappabarVCmat,1.0_realk,kappabarVC)       
    end if
    ! Relaxation components (occ-virt and valence-corr for frozen core) for density matrix
    call mem_alloc(rhoMO_relaxation,nbasis,nbasis)
    rhoMO_relaxation = 0.0_realk


    ! rho in MO basis
    ! ***************
    ! Occ-occ block   : rho_{ij} = - X_{ij}          ! This is already set in the input matrix rhoAOmat
    ! Virt-virt block : rho_{ab} = Y_{ab}            ! This is already set in the input matrix rhoAOmat
    ! Virt-occ block  : rho_{ai} = kappabar_{ai}     ! This needs to be set here
    ! Occ-virt block  : rho_{ia} = kappabar_{ai}     ! This needs to be set here
    !
    ! Matrix structure:
    ! '''''''''''''''''
    ! Occ-Occ    Occ-virt
    ! Occ-Virt   Virt-virt

    ! Special case: Frozen core (letting capital/small letter denote core/valence)
    ! ----------------------------------------------------------------------------
    !
    !     CORE        VALENCE     VIRTUAL
    !      0          kappabar_{Ii}  kappabar_{Ia}
    !  kappabar_{iI}     -X_{ij}     kappabar_{ia}
    !  kappabar_{aI}     kappabar_{ai}  Y_{ab}
    ! 



    ! Virt-occ and occ-virt blocks
    ! ****************************

    do a=1,nvirt
       ax = a+nocc
       do i=1,nocc
          rhoMO_relaxation(ax,i) = kappabar(a,i)
          rhoMO_relaxation(i,ax) = kappabar(a,i)
       enddo
    enddo


    ! Frozen core: core-valence and valence-core blocks
    ! *************************************************
    ! for simplity a/i denote valence/core indices here

    if(DECinfo%frozencore) then
       do a=1,nval
          ax = a+ncore
          do i=1,ncore
             rhoMO_relaxation(ax,i) = kappabarVC(a,i)
             rhoMO_relaxation(i,ax) = kappabarVC(a,i)
          enddo
       enddo
    end if


    call mem_dealloc(kappabar)
    if(DECinfo%frozencore) then
       call mem_dealloc(kappabarVC)
    end if


    ! Transform rho from MO to AO basis: rhoAO = C rhoMO C^T  (just relaxation part)
    ! ******************************************************************************

    ! Collect occ and virt MO coefficients
    call mem_alloc(C,nbasis,nbasis)
    do i=1,nocc
       C(:,i) = MyMolecule%Co%elm2(:,i)
    end do
    do i=1,nvirt
       C(:,i+nocc) = MyMolecule%Cv%elm2(:,i)
    end do

    ! tmp = C rhoMO_relaxation
    call mem_alloc(tmp,nbasis,nbasis)
    call dec_simple_dgemm(nbasis,nbasis,nbasis,C,rhoMO_relaxation,tmp,'n','n')

    ! rhoAO_relaxation = C rhoMO_relaxation C^T = tmp C^T
    call mem_alloc(rhoAO_relaxation,nbasis,nbasis)
    call dec_simple_dgemm(nbasis,nbasis,nbasis,tmp,C,rhoAO_relaxation,'n','t')
    call mem_dealloc(C)
    call mem_dealloc(tmp)
    call mem_dealloc(rhoMO_relaxation)

    ! Convert relaxation part to type(matrix)
    call mat_init(rhoAOmat_relaxation,nbasis,nbasis)
    call mat_set_from_full(rhoAO_relaxation, 1.0_realk, rhoAOmat_relaxation)
    call mem_dealloc(rhoAO_relaxation)

    ! Add relaxation contribution to unrelaxed input correlation density
    call mat_daxpy(1E0_realk,rhoAOmat_relaxation, rhoAOmat)
    call mat_free(rhoAOmat_relaxation)

    call LSTIMER('RHO MATRIX',tcpu,twall,DECinfo%output)


  end subroutine dec_get_rho_matrix_in_AO_basis



  !> \brief Calculate electric dipole moment at the HF and MP2 levels of theory.
  !> \author Kasper Kristensen
  !> \date August 2011
  subroutine get_HF_and_MP2_dipole_moments(mylsitem,DHF,DMP2,HFdipole,MP2dipole)


    implicit none

    !> LS setting info
    type(lsitem), intent(inout) :: mylsitem
    !> HF density matrix, (corresponds to dens.restart, thus there is NO factor 2)
    type(matrix),intent(in) :: DHF
    !> MP2 density matrix (incl HF contribution, which is multiplied by two)
    type(matrix),intent(in) :: DMP2
    type(matrix),target :: DipoleIntegral(3)
    real(realk), intent(inout) :: HFdipole(3)
    real(realk), intent(inout) :: MP2dipole(3)
    real(realk),dimension(3) :: NucDipole, HF_el_dipole, MP2_el_dipole
    integer :: nbasis,i
    character(len=7) :: string




    ! Calculate pure nuclear contribution to dipole (of course the same for HF and MP2)
    ! *********************************************************************************
    call II_get_nucdip(mylsitem%setting,NucDipole)


    ! Calculate electronic dipole matrices
    ! ************************************
    string(1:7) = 'DIPLEN '
    nbasis = DHF%nrow
    do i=1,3
       call mat_init(DipoleIntegral(i),nbasis,nbasis)
    end do
    call II_get_integral(DECinfo%output, DECinfo%output ,mylsitem%setting, DipoleIntegral,3,string)


    ! Electronic contribution to HF dipole
    ! ************************************
    ! Dot product of dipole integrals and HF density (factor two, due to double occupation for closed-shell)
    do i=1,3
       HF_el_dipole(i) = 2E0_realk*mat_dotproduct(DHF,DipoleIntegral(i))
    end do


    ! Electronic contribution to MP2 dipole
    ! ************************************
    ! Dot product of dipole integrals and MP2 density
    do i=1,3
       MP2_el_dipole(i) = mat_dotproduct(DMP2,DipoleIntegral(i))
    end do


    ! Total dipoles
    ! *************
    ! "Sum" of nuclear and electronic contributions. However, the
    ! dipole by definition is MINUS the derivative of the energy/Lagrangian with respect
    ! to an electric field component. Our density matrices are associated with an HF or an MP2 energy.
    ! Therefore, the electronic contribution should be multiplied by (-1),
    ! while the nuclear contribution has the correct sign when calculated using II_get_nucdip.

    ! HF total dipole
    do i=1,3
       HFdipole(i) = -HF_el_dipole(i) + NucDipole(i)
    end do

    ! MP2 total dipole
    do i=1,3
       MP2dipole(i) = -MP2_el_dipole(i) + NucDipole(i)
    end do



    ! Free stuff
    ! **********
    do i=1,3
       call mat_free(DipoleIntegral(i))
    end do

  end subroutine get_HF_and_MP2_dipole_moments



  !> Get kappabar multipliers (dimension nvalence,ncore) which formally ensure
  !> that core/valence block of Fock matrix is 0.
  !> \author Kasper Kristensen
  !> \date December 2012
  subroutine get_kappabar_frozen_core(MyLsitem,MyMolecule,Phioo,kappabarVC)
    implicit none
   !> LS setting info
    type(lsitem), intent(inout) :: mylsitem
    !> Molecule info
    type(fullmolecule),intent(in) :: MyMolecule
    !> RHS for kappabar equation (only valence,occ block of Phioo is used)
    !> Not changed here altough intent(inout) for practical reasons!
    type(matrix),intent(inout) :: Phioo
    !> Kappabar frozen core multipliers for valence,core block
    type(matrix),intent(inout) :: kappabarVC
    real(realk),pointer :: tmp(:,:)
    integer :: nval,ncore,nocc
    type(matrix) :: RHS,dummy1,dummy2

    nocc = MyMolecule%nocc
    ncore = MyMolecule%ncore
    nval = MyMolecule%nval

    ! Get valence,core block of Phioo and store in matrix structure
    call mem_alloc(tmp,nval,ncore)
    call mat_retrieve_block(Phioo,tmp,nval,ncore,ncore+1,1)
    call mat_init(RHS,nval,ncore)
    call mat_set_from_full(tmp,1.0_realk,RHS)
    call mem_dealloc(tmp)

    ! Solve equation
    ! **************
    call mat_zero(kappabarVC)      ! Use zero matrix as starting guess
    call dec_solve_kappabar_equation(RHS,.false.,&
         &MyLsitem,MyMolecule,dummy1,dummy2,kappabarVC,.true.)
    call mat_free(RHS)

  end subroutine get_kappabar_frozen_core

  
  !> \brief Calculate MP2 correlation density matrix in AO basis
  !> from X and Y matrices in MO basis (see type mp2dens).
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine get_unrelaxed_corrdens_in_AO_basis(nocc,nvirt,nbasis,Cocc,Cvirt,X,Y,rho)
    implicit none
    !> Number of occupied orbitals, virtual orbitals, and basis functions
    integer,intent(in) :: nocc,nvirt,nbasis
    !> Occupied MO coefficients
    real(realk),intent(in) :: Cocc(nbasis,nocc)
    !> Virtual MO coefficients
    real(realk),intent(in) :: Cvirt(nbasis,nvirt)
    !> (Minus) occ-occ block of correlation density matrix (see type mp2dens)
    real(realk),intent(in) :: X(nocc,nocc)
    !> Virt-virt block of correlation density matrix (see type mp2dens)
    real(realk),intent(in) :: Y(nvirt,nvirt)
    !> Unrelaxed correlation density in AO basis (see rho in type mp2dens WITHOUT the kappa elements)
    real(realk),intent(inout) :: rho(nbasis,nbasis)
    real(realk),pointer :: XAO(:,:),YAO(:,:)

    ! Get X matrix in AO basis
    call mem_alloc(XAO,nbasis,nbasis)
    call dec_simple_basis_transform2(nbasis,nocc,Cocc,X,XAO)

    ! Get Y matrix in AO basis
    call mem_alloc(YAO,nbasis,nbasis)
    call dec_simple_basis_transform2(nbasis,nvirt,Cvirt,Y,YAO)

    ! rho = Y-X (see type mp2dens)
    rho = YAO-XAO

    call mem_dealloc(XAO)
    call mem_dealloc(YAO)

  end subroutine get_unrelaxed_corrdens_in_AO_basis



  !> \brief Get fragment Phi matrix in AO basis from MO blocks,
  !> wrapper to hide frozen core stuff behind the curtain.
  !> Assumes that MO blocks are stored in grad%Phioo, grad%Phivv, grad%dens%Phivo, grad%dens%Phiov
  !> at input, and Phi AO matrix is stored in grad%PhiAO matrix at output.
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine fragment_Phi_matrix_in_AO_basis_wrapper(fragment,grad)
    implicit none

    !> Fragment under consideration
    type(decfrag),intent(in) :: fragment
    !> MP2 gradient structure for fragment
    type(mp2grad), intent(inout) :: grad

    ! Frozen core requires special treatment of core orbitals
    if(DECinfo%frozencore) then
       call fragment_Phi_matrix_in_AO_basis_fc(fragment%ncore,fragment%noccAOS,&
            & fragment%nvirtAOS,fragment%nbasis,fragment%CoreMO,fragment%Co,&
            & fragment%Cv,grad)
    else
       call fragment_Phi_matrix_in_AO_basis_standard(fragment%noccAOS,fragment%nvirtAOS,&
            & fragment%nbasis,fragment%Co,fragment%Cv,grad)
    end if

  end subroutine fragment_Phi_matrix_in_AO_basis_wrapper



  !> \brief Get fragment Phi matrix in AO basis from MO blocks.
  !> Assumes that MO blocks are stored in grad%Phioo, grad%Phivv, grad%dens%Phivo, grad%dens%Phiov
  !> at input, and Phi AO matrix is stored in grad%PhiAO matrix at output.
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine fragment_Phi_matrix_in_AO_basis_fc(ncore,nval,nvirt,nbasis,Ccore,Cval,Cvirt,grad)
    implicit none
    !> Number of core orbitals, valence orbitals, 
    !> virtual orbitals, and basis functions
    integer,intent(in) :: ncore,nval,nvirt,nbasis
    !> Occupied core MO coefficients
    real(realk),intent(in) :: Ccore(nbasis,ncore)
    !> Occupied valence MO coefficients
    real(realk),intent(in) :: Cval(nbasis,nval)
    !> Virtual MO coefficients
    real(realk),intent(in) :: Cvirt(nbasis,nvirt)
    !> MP2 gradient structure for fragment
    type(mp2grad), intent(inout) :: grad
    real(realk),pointer :: tmp(:,:),Cocc(:,:)
    integer :: nb2,i,nocc
    real(realk) :: alpha

    ! For frozen core we need to consider both occupied valence
    ! orbitals separately (e.g. Phiov) and the total set of occupied 
    ! MOs (core+valence) (e.g. for Phivo), see single_init_mp2dens and init_mp2grad.

    ! Total number of occupied orbitals (core+valence)
    nocc = ncore+nval
    call mem_alloc(Cocc,nbasis,nocc)

    ! Copy core MO coefficients
    do i=1,ncore
       Cocc(:,i) = Ccore(:,i)
    end do

    ! Copy valence MO coefficients
    do i=1,nval
       Cocc(:,i+ncore) = Cval(:,i)
    end do

    ! Get Phioo matrix in AO basis 
    ! First dimension is valence, second dimension is core+valence, see init_mp2grad
    call dec_diff_basis_transform2(nbasis,nval,nocc,Cval,Cocc,grad%Phioo,grad%PhiAO)

    call mem_alloc(tmp,nbasis,nbasis)
    nb2 = nbasis*nbasis
    alpha=1.0_realk

    ! Add Phivv matrix in AO basis
    call dec_simple_basis_transform2(nbasis,nvirt,Cvirt,grad%Phivv,tmp)
    call daxpy(nb2,alpha,tmp,1,grad%PhiAO,1)

    ! Add Phivo matrix in AO basis  (second dimension is core+valence, see single_init_mp2dens)
    call dec_diff_basis_transform2(nbasis,nvirt,nocc,Cvirt,&
         & Cocc,grad%dens%Phivo,tmp)
    call daxpy(nb2,alpha,tmp,1,grad%PhiAO,1)

    ! Add Phiov matrix in AO basis (first dimension is just valence, see single_init_mp2dens)
    call dec_diff_basis_transform2(nbasis,nval,nvirt,Cval,&
         & Cvirt,grad%dens%Phiov,tmp)
    call daxpy(nb2,alpha,tmp,1,grad%PhiAO,1)

    call mem_dealloc(tmp)
    call mem_dealloc(Cocc)

  end subroutine fragment_Phi_matrix_in_AO_basis_fc



  !> \brief Get fragment Phi matrix in AO basis from MO blocks.
  !> Assumes that MO blocks are stored in grad%Phioo, grad%Phivv, grad%dens%Phivo, grad%dens%Phiov
  !> at input, and Phi AO matrix is stored in grad%PhiAO matrix at output.
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine fragment_Phi_matrix_in_AO_basis_standard(nocc,nvirt,nbasis,Cocc,Cvirt,grad)
    implicit none
    !> Number of occupied orbitals, virtual orbitals, and basis functions
    integer,intent(in) :: nocc,nvirt,nbasis
    !> Occupied MO coefficients
    real(realk),intent(in) :: Cocc(nbasis,nocc)
    !> Virtual MO coefficients
    real(realk),intent(in) :: Cvirt(nbasis,nvirt)
    !> MP2 gradient structure for fragment
    type(mp2grad), intent(inout) :: grad
    real(realk),pointer :: tmp(:,:)
    integer :: nb2
    real(realk) :: alpha

    ! Get Phioo matrix in AO basis
    call dec_simple_basis_transform2(nbasis,nocc,Cocc,grad%Phioo,grad%PhiAO)

    call mem_alloc(tmp,nbasis,nbasis)
    nb2 = nbasis*nbasis
    alpha=1.0_realk

    ! Add Phivv matrix in AO basis
    call dec_simple_basis_transform2(nbasis,nvirt,Cvirt,grad%Phivv,tmp)
    call daxpy(nb2,alpha,tmp,1,grad%PhiAO,1)

    ! Add Phivo matrix in AO basis
    call dec_diff_basis_transform2(nbasis,nvirt,nocc,Cvirt,&
         & Cocc,grad%dens%Phivo,tmp)
    call daxpy(nb2,alpha,tmp,1,grad%PhiAO,1)

    ! Add Phiov matrix in AO basis
    call dec_diff_basis_transform2(nbasis,nocc,nvirt,Cocc,&
         & Cvirt,grad%dens%Phiov,tmp)
    call daxpy(nb2,alpha,tmp,1,grad%PhiAO,1)

    call mem_dealloc(tmp)

  end subroutine fragment_Phi_matrix_in_AO_basis_standard

  
  !> \brief Calculate occ-virt and virt-occ blocks of Phi
  !> matrix in MO basis from Phi matrix in AO basis.
  !> For frozen core approx occ-occ block is also required.
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine get_Phi_MO_blocks_from_AO_simple(nbasis,nocc,nvirt,S,&
       & Cocc,Cvirt,PhiAO,Phivo,Phiov,Phioo)

    implicit none
    !> Number of occupied orbitals, virtual orbitals, and basis functions
    integer,intent(in) :: nocc,nvirt,nbasis
    !> Overlap matrix
    real(realk),dimension(nbasis,nbasis),intent(in) :: S
    !> Occupied MO coefficients
    real(realk),dimension(nbasis,nocc),intent(in) :: Cocc
    !> Virtual MO coefficients
    real(realk),dimension(nbasis,nvirt),intent(in) :: Cvirt
    !> Phi matrix in AO basis
    real(realk),dimension(nbasis,nbasis),intent(in) :: PhiAO
    !> Virt-occ block of Phi matrix in MO basis
    real(realk),dimension(nvirt,nocc),intent(inout) :: Phivo
    !> Occ-virt block of Phi matrix in MO basis
    real(realk),dimension(nocc,nvirt),intent(inout) :: Phiov
    !> Occ-occ block of Phi matrix (only for frozen core)
    real(realk),dimension(nocc,nocc),intent(inout) :: Phioo
    real(realk),pointer :: S_PhiAO_S(:,:)


    ! Phi in AO and MO bases are related as follows:
    ! (see also fragment_Phi_matrix_in_AO_basis)
    ! PhiAO = C PhiMO C^T
    ! PhiMO = C^T S PhiAO S C    (*)

    ! Now we want to calculate blocks of PhiMO so the procedure is:
    ! 1. Calculate the intermediate S PhiAO S
    ! 2. Insert the relevant MO coefficients in (*) to get
    !    the relevant MO blocks of Phi.



    ! Calculate S PhiAO S
    ! *******************

    call mem_alloc(S_PhiAO_S,nbasis,nbasis)
    call dec_simple_basis_transform1(nbasis,nbasis,S,PhiAO,S_PhiAO_S)
    ! (we can use dec_simple_basis_transform1 for something which is 
    !  not its original purpose because S is symmetric)
    

    ! Transform to MO basis
    ! *********************

    ! Phivo = Cvirt^T (S PhiAO S) Cocc
    call dec_diff_basis_transform1(nbasis,nvirt,nocc,Cvirt,Cocc,S_PhiAO_S,Phivo)

    ! Phiov = Cocc^T (S PhiAO S) Cvirt
    call dec_diff_basis_transform1(nbasis,nocc,nvirt,Cocc,Cvirt,S_PhiAO_S,Phiov)

    ! Phioo = Cocc^T (S PhiAO S) Cocc
    call dec_simple_basis_transform1(nbasis,nocc,Cocc,S_PhiAO_S,Phioo)

    call mem_dealloc(S_PhiAO_S)

  end subroutine get_Phi_MO_blocks_from_AO_simple



  !> \brief Calculate occ-virt and virt-occ blocks of Phi
  !> matrix in MO basis from full molecular structure and
  !> full molecular gradient structure.
  !> For frozen core approx occ-occ block is also required.
  !> Matrices are also intialized here!
  !> \author Kasper Kristesen
  !> \date March 2013
  subroutine get_Phi_MO_blocks(mylsitem,MyMolecule,fullgrad,Phivo,Phiov,Phioo)

    implicit none
    !> LS setting info
    type(lsitem), intent(inout) :: mylsitem
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Full MP2 gradient structure
    type(fullmp2grad),intent(in) :: fullgrad
    !> Virt-occ blocks of Phi matrix
    type(matrix),intent(inout) :: Phivo
    !> Occ-virt blocks of Phi matrix
    type(matrix),intent(inout) :: Phiov
    !> Occ-occ blocks of Phi matrix (only for frozen core)
    type(matrix),intent(inout),optional :: Phioo
    integer :: nbasis,nocc,nvirt
    real(realk),pointer :: Phivo_simple(:,:), Phiov_simple(:,:), &
         & Phioo_simple(:,:),S(:,:)

    ! Sanity check
    if(DECinfo%frozencore) then
       if(.not. present(Phioo)) then
          call lsquit('get_Phi_MO_blocks : &
               & Phioo must be present for frozen core!',-1)
       end if
    end if


    nbasis = MyMolecule%nbasis
    nocc = MyMolecule%nocc
    nvirt = MyMolecule%nvirt

    ! AO overlap matrix
    call mem_alloc(S,nbasis,nbasis)
    call II_get_mixed_overlap_full(DECinfo%output,DECinfo%output,mylsitem%SETTING,&
         & S,nbasis,nbasis,AORdefault,AORdefault)


    ! Get MO blocks for Phi matrix in simple fortran form
    call mem_alloc(Phivo_simple,nvirt,nocc)
    call mem_alloc(Phiov_simple,nocc,nvirt)
    call mem_alloc(Phioo_simple,nocc,nocc)
    call get_Phi_MO_blocks_from_AO_simple(nbasis,nocc,nvirt,S,&
         & MyMolecule%Co%elm2,MyMolecule%Cv%elm2,fullgrad%Phi,&
         & Phivo_simple,Phiov_simple,Phioo_simple)
    call mem_dealloc(S)

    ! Init and convert to type(matrix) form
    call mat_init(Phivo,nvirt,nocc)
    call mat_set_from_full(Phivo_simple, 1.0_realk, Phivo)
    call mem_dealloc(Phivo_simple)

    call mat_init(Phiov,nocc,nvirt)
    call mat_set_from_full(Phiov_simple, 1.0_realk, Phiov)
    call mem_dealloc(Phiov_simple)

    if(DECinfo%frozencore) then
       call mat_init(Phioo,nocc,nocc)
       call mat_set_from_full(Phioo_simple, 1.0_realk, Phioo)
    end if
    call mem_dealloc(Phioo_simple)

  end subroutine get_Phi_MO_blocks


end module mp2_gradient_module
