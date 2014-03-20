!> @file
!> DEC-CCSD(T) routines
!> \brief: ccsd(t) module
!> \author: Janus Juul Eriksen
!> \date: 2012-2013, Aarhus
module ccsdpt_module

#ifdef VAR_MPI
      use infpar_module
      use lsmpi_type
#endif
  use precision
  use dec_typedef_module
  use memory_handling
  use lstiming!, only: lstimer
  use screen_mod!, only: DECscreenITEM
  use BUILDAOBATCH
  use typedeftype!, only: Lsitem,lssetting
  use screen_mod!,only: free_decscreen, DECSCREENITEM
  use IntegralInterfaceDEC!, only: II_precalc_DECScreenMat,&
!       & II_getBatchOrbitalScreen, II_GET_DECPACKED4CENTER_J_ERI
  use IntegralInterfaceMOD
  use Fundamental, only: bohr_to_angstrom
  use tensor_interface_module


  ! DEC DEPENDENCIES (within deccc directory)  
  ! *****************************************
#ifdef VAR_MPI
      use decmpi_module
#endif
  use dec_fragment_utils
  use array2_simple_operations
  use array3_simple_operations
  use array4_simple_operations
  
  public :: ccsdpt_driver,ccsdpt_energy_e4_frag,ccsdpt_energy_e5_frag,&
       & ccsdpt_energy_e4_pair, ccsdpt_energy_e5_pair,&
       & ccsdpt_energy_e4_full, print_e4_full, ccsdpt_energy_e5_full,&
       & print_e5_full
  private

contains

#ifdef MOD_UNRELEASED

  !> \brief: driver routine for dec-ccsd(t)
  !> \author: Janus Juul Eriksen
  !> \date: july 2012
  subroutine ccsdpt_driver(nocc,nvirt,nbasis,ppfock,qqfock,Co,Cv,mylsitem,ccsd_doubles,&
                         & ccsdpt_singles,ccsdpt_doubles)

    implicit none

    !> nocc, nvirt, and nbasis for fragment or full molecule
    integer, intent(in) :: nocc, nvirt, nbasis
    !> ppfock and qqfock for fragment or full molecule
    real(realk), intent(in) :: ppfock(nocc,nocc), qqfock(nvirt,nvirt)
    !> mo coefficents for occ and virt space for fragment or full molecule
    real(realk), intent(in) :: Co(nbasis,nocc), Cv(nbasis,nvirt)
    !> mylsitem for fragment or full molecule
    type(lsitem), intent(inout) :: mylsitem
    !> ccsd doubles amplitudes
    type(array4), intent(inout) :: ccsd_doubles
    !> 2-el integrals
    type(array4) :: jaik ! integrals (AI|JK) in the order (J,A,I,K)
    type(array4) :: abij ! integrals (AI|BJ) in the order (A,B,I,J)
    ! cbai is of type DENSE, if this is a serial calculation, and TILED_DIST,
    ! if this is a parallel calculation
    type(array) :: cbai ! integrals (AI|BC) in the order (C,B,A,I)
#ifdef VAR_MPI
    integer :: nodtotal
    real(realk) :: jaik_norm, abij_norm, cbai_norm
    real(realk) :: ccsdpt_doubles_norm, ccsdpt_doubles_2_norm, ccsdpt_singles_norm
#endif
    !> orbital energies
    real(realk), pointer :: eivalocc(:), eivalvirt(:)
    !> MOs and unitary transformation matrices
    type(array2) :: C_can_occ, C_can_virt, Uocc, Uvirt
    !> dimensions
    integer, dimension(2) :: occdims, virtdims, virtoccdims,occAO,virtAO
    integer, dimension(3) :: dims_aaa
    integer, dimension(4) :: dims_iaai, dims_aaii
    !> input for the actual triples computation
    type(array4) :: ccsdpt_doubles_2
    type(array4),intent(inout) :: ccsdpt_doubles
    type(array2),intent(inout) :: ccsdpt_singles

    ! init dimensions
    occdims     = [nocc,nocc]
    virtdims    = [nvirt,nvirt]
    virtoccdims = [nvirt,nocc]
    dims_iaai   = [nocc,nvirt,nvirt,nocc]
    dims_aaii   = [nvirt,nvirt,nocc,nocc]
    dims_aaa    = [nvirt,nvirt,nvirt]
    occAO       = [nbasis,nocc]
    virtAO      = [nbasis,nvirt]

#ifdef VAR_MPI

    nodtotal = infpar%lg_nodtot

    ! bcast the JOB specifier and distribute data to all the slaves within local group
    waking_the_slaves: if ((nodtotal .gt. 1) .and. (infpar%lg_mynum .eq. infpar%master)) then

       ! slaves are in lsmpi_slave routine (or corresponding dec_mpi_slave) and are now awaken
       call ls_mpibcast(CCSDPTSLAVE,infpar%master,infpar%lg_comm)

       ! distribute ccsd doubles and fragment or full molecule quantities to the slaves
       call mpi_communicate_ccsdpt_calcdata(nocc,nvirt,nbasis,ppfock,qqfock,Co,Cv,ccsd_doubles%val,mylsitem)

    end if waking_the_slaves

#endif

    ! *************************************
    ! get arrays for transforming integrals
    ! *************************************
    ! C_can_occ, C_can_virt:  MO coefficients for canonical basis
    ! Uocc, Uvirt: unitary transformation matrices for canonical --> local basis (and vice versa)
    ! note: Uocc and Uvirt have indices (local,canonical)

    call mem_alloc(eivalocc,nocc)
    call mem_alloc(eivalvirt,nvirt)
    Uocc       = array2_init(occdims)
    Uvirt      = array2_init(virtdims)
    C_can_occ  = array2_init(occAO)
    C_can_virt = array2_init(virtAO)
    call get_canonical_integral_transformation_matrices(nocc,nvirt,nbasis,ppfock,qqfock,Co,Cv,&
                         & C_can_occ%val,C_can_virt%val,Uocc%val,Uvirt%val,eivalocc,eivalvirt)

    ! ***************************************************
    ! get vo³, v²o², and v³o integrals in proper sequence
    ! ***************************************************
    ! note: the integrals are calculated in canonical basis

    call get_CCSDpT_integrals(mylsitem,nbasis,nocc,nvirt,C_can_occ%val,C_can_virt%val,jaik,abij,cbai)

#ifdef VAR_MPI

    print *,'proc no. ',infpar%lg_mynum,'after get_CCSDpT_integrals'
    call array4_print_norm_nrm(jaik,jaik_norm)
    call array4_print_norm_nrm(abij,abij_norm)
    call array_print_norm_nrm(cbai,cbai_norm)
    print *,'proc no. ',infpar%lg_mynum,'jaik_norm = ',jaik_norm
    print *,'proc no. ',infpar%lg_mynum,'abij_norm = ',abij_norm
    print *,'proc no. ',infpar%lg_mynum,'cbai_norm = ',cbai_norm

#endif

    ! release occ and virt canonical MOs
    call array2_free(C_can_occ)
    call array2_free(C_can_virt)

    ! ***************************************************
    ! transform ccsd doubles amplitudes to diagonal basis
    ! ***************************************************

    call ccsdpt_local_can_trans(ccsd_doubles,nocc,nvirt,Uocc,Uvirt)

#ifdef VAR_MPI

    print *,'proc no. ',infpar%lg_mynum,'after ccsdpt_local_can_trans'
    call array4_print_norm_nrm(jaik,jaik_norm)
    call array4_print_norm_nrm(abij,abij_norm)
    call array_print_norm_nrm(cbai,cbai_norm)
    print *,'proc no. ',infpar%lg_mynum,'jaik_norm = ',jaik_norm
    print *,'proc no. ',infpar%lg_mynum,'abij_norm = ',abij_norm
    print *,'proc no. ',infpar%lg_mynum,'cbai_norm = ',cbai_norm

#endif

    ! Now we transpose the unitary transformation matrices as we will need these in the transformation
    ! of the ^{ccsd}T^{ab}_{ij}, ^{*}T^{a}_{i}, and ^{*}T^{ab}_{ij} amplitudes from canonical to local basis
    ! later on
    call array2_transpose(Uocc)
    call array2_transpose(Uvirt)

    ! ********************************
    ! begin actual triples calculation
    ! ********************************

    ! in all comments in the below, we employ the notation of eqs. (14.6.60) [with (i,j,k)/(a,c,d)]
    ! and (14.6.64).

    ! objective is three-fold:
    ! 1) calculate triples amplitudes, collect in array3 structures, trip_*** [canonical basis]
    ! 2) calculate ^{*}T^{a}_{i} and ^{*}T^{ab}_{ij} amplitudes in array2 and array4 structures, 
    !    ccsdpt_singles and ccsdpt_doubles [canonical basis]
    !    here: ccsdpt_doubles_2 is a temp array towards the generation of ccsdpt_doubles
    ! 3) transform ccsd_doubles, ccsdpt_singles and ccsdpt_doubles into local basis [local basis]

    ! *****************************************************
    ! ***************** trip generation *******************
    ! *****************************************************

    ! init ccsdpt_doubles_2 array4 structure.
    ! we merge ccsdpt_doubles and ccsdpt_doubles_2 at the end into ccsdpt_doubles. 
    ! we have dimensioned ccsdpt_doubles as dims_aaii and ccsdpt_doubles_2 as dims_iaai 
    ! in order to load in data consecutive in memory inside ccsdpt_contract_21 
    ! and ccsdpt_contract_22, respectively.
    ccsdpt_doubles_2 = array4_init_standard(dims_iaai)

    ! initially, reorder ccsd_doubles
    call array4_reorder(ccsd_doubles,[3,1,4,2]) ! ccsd_doubles(a,i,b,j) --> ccsd_doubles(b,a,j,i)

    ! if cbai is tiled distributed, then put the three tiles into
    ! an array structure, cbai_pdm. here, initialize the array structure.

    !**********************************************************
    ! here: the main ijk-loop: this is where the magic happens!
    !**********************************************************

#ifdef VAR_MPI

    ! the parallel version of the ijk-loop
    call ijk_loop_par(nocc,nvirt,jaik,abij,cbai,ccsd_doubles,&
                    & ccsdpt_doubles%val,ccsdpt_doubles_2%val,ccsdpt_singles%val,eivalocc,eivalvirt,nodtotal)

#ifdef VAR_MPI

    print *,'proc no. ',infpar%lg_mynum,'after ijk_loop_par'
    call array4_print_norm_nrm(jaik,jaik_norm)
    call array4_print_norm_nrm(abij,abij_norm)
    call array_print_norm_nrm(cbai,cbai_norm)
    call array4_print_norm_nrm(ccsdpt_doubles,ccsdpt_doubles_norm)
    call array4_print_norm_nrm(ccsdpt_doubles_2,ccsdpt_doubles_2_norm)
    call array2_print_norm_nrm(ccsdpt_singles,ccsdpt_singles_norm)
    print *,'proc no. ',infpar%lg_mynum,'jaik_norm = ',jaik_norm
    print *,'proc no. ',infpar%lg_mynum,'abij_norm = ',abij_norm
    print *,'proc no. ',infpar%lg_mynum,'cbai_norm = ',cbai_norm
    print *,'proc no. ',infpar%lg_mynum,'ccsdpt_doubles_norm = ',ccsdpt_doubles_norm
    print *,'proc no. ',infpar%lg_mynum,'ccsdpt_doubles_2_norm = ',ccsdpt_doubles_2_norm
    print *,'proc no. ',infpar%lg_mynum,'ccsdpt_singles_norm = ',ccsdpt_singles_norm

#endif

#else

    ! the serial version of the ijk-loop
    call ijk_loop_ser(nocc,nvirt,jaik,abij,cbai,ccsd_doubles,&
                    & ccsdpt_doubles%val,ccsdpt_doubles_2%val,ccsdpt_singles%val,eivalocc,eivalvirt)

#endif

    ! *******************************************
    ! *********** done w/ main loop *************
    ! *******************************************

#ifdef VAR_MPI

    ! here, synchronize all procs
    call lsmpi_barrier(infpar%lg_comm)

    ! reduce singles and doubles arrays into that residing on the master
    reducing_to_master: if (nodtotal .gt. 1) then

       call lsmpi_local_reduction(ccsdpt_singles%val,nocc,nvirt,infpar%master)
       call lsmpi_local_reduction(ccsdpt_doubles%val,nvirt,nocc,nvirt,nocc,infpar%master)
       call lsmpi_local_reduction(ccsdpt_doubles_2%val,nvirt,nocc,nvirt,nocc,infpar%master)

    end if reducing_to_master

    if (infpar%lg_mynum .eq. infpar%master) then

       print *,'proc no. ',infpar%lg_mynum,'after lsmpi_local_reduction'
       call array4_print_norm_nrm(ccsdpt_doubles,ccsdpt_doubles_norm)
       call array4_print_norm_nrm(ccsdpt_doubles_2,ccsdpt_doubles_2_norm)
       call array2_print_norm_nrm(ccsdpt_singles,ccsdpt_singles_norm)
       print *,'proc no. ',infpar%lg_mynum,'ccsdpt_doubles_norm = ',ccsdpt_doubles_norm
       print *,'proc no. ',infpar%lg_mynum,'ccsdpt_doubles_2_norm = ',ccsdpt_doubles_2_norm
       print *,'proc no. ',infpar%lg_mynum,'ccsdpt_singles_norm = ',ccsdpt_singles_norm

    end if

    ! release stuff located on slaves
    releasing_the_slaves: if ((nodtotal .gt. 1) .and. (infpar%lg_mynum .ne. infpar%master)) then

       ! release stuff initialized herein
       call array2_free(Uocc)
       call array2_free(Uvirt)
       call array4_free(ccsdpt_doubles_2) 
       call mem_dealloc(eivalocc)
       call mem_dealloc(eivalvirt)
       call array4_free(abij)
       call array_free(cbai)
       call array4_free(jaik)

       ! now, release the slaves  
       return

    end if releasing_the_slaves

#endif

    ! now everything resides on the master...

    ! collect ccsdpt_doubles and ccsdpt_doubles_2 into ccsdpt_doubles array4 structure
    ! ccsdpt_doubles(a,b,i,j) = ccsdpt_doubles(a,b,i,j) + ccsdpt_doubles_2(j,a,b,i) (*)
    ! (*) here, ccsdpt_doubles_2 is simultaneously reordered as (j,a,b,i) --> (a,b,i,j)
    call array_reorder_4d(1.0E0_realk,ccsdpt_doubles_2%val,ccsdpt_doubles_2%dims(1),&
                               &ccsdpt_doubles_2%dims(2),ccsdpt_doubles_2%dims(3),ccsdpt_doubles_2%dims(4),&
                               &[2,3,4,1],1.0E0_realk,ccsdpt_doubles%val)

    ! release ccsdpt_doubles_2 array4 structure
    call array4_free(ccsdpt_doubles_2)

    ! *************************************************
    ! ***** do canonical --> local transformation *****
    ! *************************************************

    call ccsdpt_can_local_trans(ccsd_doubles,ccsdpt_singles,ccsdpt_doubles,nocc,nvirt,Uocc,Uvirt)

    ! now, release Uocc and Uvirt
    call array2_free(Uocc)
    call array2_free(Uvirt)

    ! clean up
    call mem_dealloc(eivalocc)
    call mem_dealloc(eivalvirt)
    call array4_free(abij)
    call array_free(cbai)
    call array4_free(jaik)

    ! **************************************************************
    ! *** do final reordering of amplitudes and clean the dishes ***
    ! **************************************************************

    ! reorder ccsdpt_doubles and ccsd_doubles back to (a,b,i,j) sequence
    call array4_reorder(ccsdpt_doubles,[3,4,1,2])
    call array4_reorder(ccsd_doubles,[4,3,2,1])

  end subroutine ccsdpt_driver


  !> \brief: main ijk-loop (parallel version)
  !> \author: Janus Juul Eriksen
  !> \date: january 2014
  subroutine ijk_loop_par(nocc,nvirt,ovoo,vvoo,vvvo,ccsd_doubles,&
                        & ccsdpt_doubles,ccsdpt_doubles_2,ccsdpt_singles,eivalocc,eivalvirt,nodtotal)

    implicit none

    !> nocc,nvirt
    integer, intent(in)      :: nocc,nvirt
    !> 2-el integrals
    type(array4), intent(inout) :: ovoo ! integrals (AI|JK) in the order (J,A,I,K)
    type(array4), intent(inout) :: vvoo ! integrals (AI|BJ) in the order (A,B,I,J)
    type(array), intent(inout)  :: vvvo ! integrals (AI|BC) in the order (C,B,A,I)
    type(array)                 :: vvvo_pdm ! v^3 tiles from cbai, 1 == i, 2 == j, 3 == k
    !> ccsd doubles amplitudes
    type(array4), intent(inout) :: ccsd_doubles
    ! o*v^2 portions of ccsd_doubles
    real(realk), pointer, dimension(:,:,:) :: ccsd_doubles_portions_i,ccsd_doubles_portions_j,ccsd_doubles_portions_k
    !> triples amplitudes and 3d work array
    real(realk), pointer, dimension(:,:,:) :: trip_tmp, trip_ampl
    !> ccsd(t) intermediates 
    real(realk), dimension(nvirt,nvirt,nocc,nocc) :: ccsdpt_doubles
    real(realk), dimension(nocc,nvirt,nvirt,nocc) :: ccsdpt_doubles_2
    real(realk), dimension(nvirt,nocc) :: ccsdpt_singles
    !> orbital energies
    real(realk), intent(inout)  :: eivalocc(nocc), eivalvirt(nvirt)
    !> job distribution
    integer :: b_size,njobs,nodtotal,ij,ij_count,i_old,j_old
    integer, pointer :: ij_array(:),jobs(:)
    !> loop integers
    integer :: i,j,k,idx,tuple_type

    ! init ccsd_doubles_help_arrays
    call mem_alloc(ccsd_doubles_portions_i,nocc,nvirt,nvirt)
    call mem_alloc(ccsd_doubles_portions_j,nocc,nvirt,nvirt)
    call mem_alloc(ccsd_doubles_portions_k,nocc,nvirt,nvirt)

    ! init triples tuples structure
    call mem_alloc(trip_ampl,nvirt,nvirt,nvirt)
    ! init 3d wrk array
    call mem_alloc(trip_tmp,nvirt,nvirt,nvirt)

    ! init pdm work array for vvvo integrals
    vvvo_pdm = array_init([nvirt,nvirt,nvirt,3],4)

    ! create job distribution list
    ! first, determine common batch size from number of tasks and nodes
    ! in the ij matrix, njobs is the number of elements in the lower triangular matrix
    ! always an even number [ n(n+1) is always an even number ]
    njobs = int((nocc**2 + nocc)/2)
    b_size = int(njobs/nodtotal)

    ! ij_array stores all jobs for composite ij indices in descending order
    call mem_alloc(ij_array,njobs)
    ! init list (one more than b_size since mod(njobs,nodtotal) is not necessearily zero
    call mem_alloc(jobs,b_size + 1)

    ! create ij_array
    call create_ij_array_ccsdpt(njobs,nocc,ij_array)
    ! fill the list
    call job_distrib_ccsdpt(b_size,njobs,ij_array,jobs)

    ! release ij_array
    call mem_dealloc(ij_array)

    ! now follows the main loop

    ! a note on the mpi scheme.
    ! since we (in a dec picture) often have many nodes compared to nocc, we explicitly collapse the i- and j-loop.
    ! by doing this, we are guaranteed that all nodes participate.
    ! the composite index ij is incremented in the collapsed loop, and we may calculate i and j from ij.

    ! a note on the openacc implementation.
    ! we use pointers to alias specific memory blocks of the vvoo and ovoo integrals.
    ! this is because we thus avoid to allocate more memory on the device for these integral blocks than what is
    ! actually needed within one loop cycle AND we avoid doing explicit copies into memory buffers on the host.

    ! init ij and i_old/j_old
    ij = 0
    i_old = 0
    j_old = 0

 ijrun_par: do ij_count = 1,b_size + 1

               ! get value of ij from job disttribution list
               ij = jobs(ij_count)
    
               ! no more jobs to be done? otherwise leave the loop
               if (ij .lt. 0) exit

               ! calculate i and j from composite ij value
               call calc_i_and_j(ij,nocc,i,j)
    
               ! has the i and j index changed?
               if (i .eq. i_old) then
    
                  ! get the j'th v^3 tile only
                  call array_get_tile(vvvo,j,vvvo_pdm%elm1(nvirt**3+1:2*nvirt**3),nvirt**3)
    
                  ! store portion of ccsd_doubles (the j'th index) to avoid unnecessary reorderings
                  call array_reorder_3d(1.0E0_realk,ccsd_doubles%val(:,:,:,j),nvirt,nvirt,&
                          & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_j)
    
                  ! store j index
                  j_old = j
    
               else if (j .eq. j_old) then
    
                  ! get the i'th v^3 tile only
                  call array_get_tile(vvvo,i,vvvo_pdm%elm1(1:nvirt**3),nvirt**3)
    
                  ! store portion of ccsd_doubles (the i'th index) to avoid unnecessary reorderings
                  call array_reorder_3d(1.0E0_realk,ccsd_doubles%val(:,:,:,i),nvirt,nvirt,&
                          & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_i)
    
                  ! store i index
                  i_old = i
    
               else
    
                  ! get the i'th and j'th v^3 tile
                  call array_get_tile(vvvo,i,vvvo_pdm%elm1(1:nvirt**3),nvirt**3)
                  call array_get_tile(vvvo,j,vvvo_pdm%elm1(nvirt**3+1:2*nvirt**3),nvirt**3)
    
                  ! store portion of ccsd_doubles (the i'th index) to avoid unnecessary reorderings
                  call array_reorder_3d(1.0E0_realk,ccsd_doubles%val(:,:,:,i),nvirt,nvirt,&
                          & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_i)
    
                  ! store portion of ccsd_doubles (the j'th index) to avoid unnecessary reorderings
                  call array_reorder_3d(1.0E0_realk,ccsd_doubles%val(:,:,:,j),nvirt,nvirt,&
                          & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_j)
    
                  ! store i and j indices
                  i_old = i
                  j_old = j
    
               end if

        krun_par: do k=1,j
 
                     ! get the k'th tile
                     call array_get_tile(vvvo,k,vvvo_pdm%elm1(2*nvirt**3+1:3*nvirt**3),nvirt**3)
     
                     ! store portion of ccsd_doubles (the k'th index) to avoid unnecessary reorderings
                     call array_reorder_3d(1.0E0_realk,ccsd_doubles%val(:,:,:,k),nvirt,nvirt,&
                             & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_k)
     
                     ! select type of tuple
                     tuple_type = -1
                     ! i == j == k
                     ! this always gives zero contribution
                     ! i == j > k
                     if ((i .eq. j) .and. (j .gt. k) .and. (i .gt. k)) tuple_type = 1
                     ! i > j == k
                     if ((i .gt. j) .and. (j .eq. k) .and. (i .gt. k)) tuple_type = 2
                     ! i > j > k
                     if ((i .gt. j) .and. (j .gt. k) .and. (i .gt. k)) tuple_type = 3
     
                     ! generate tuple(s)
                     TypeOfTuple_par: select case(tuple_type)
     
                     case(1)

                        call trip_generator_case1(i,k,nocc,nvirt,ccsd_doubles%val(:,:,i,i),ccsd_doubles%val(:,:,i,k),&
                                                & ccsd_doubles%val(:,:,k,i),ccsd_doubles_portions_i,&
                                                & ccsd_doubles_portions_k,&
                                                & vvvo_pdm%elm4(:,:,:,1),vvvo_pdm%elm4(:,:,:,3),&
                                                & ovoo%val(:,:,i,i),ovoo%val(:,:,i,k),ovoo%val(:,:,k,i),&
                                                & trip_tmp,trip_ampl)
     
                        ! generate triples amplitudes from trip arrays
     
                        call trip_denom(i,i,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl)
     
                        ! now do the contractions
     
                        call ccsdpt_driver_case1(i,k,nocc,nvirt,vvoo%val(:,:,i,i),vvoo%val(:,:,i,k),vvoo%val(:,:,k,i),&
                                             & ovoo%val(:,:,i,i),ovoo%val(:,:,i,k),ovoo%val(:,:,k,i),&
                                             & vvvo_pdm%elm4(:,:,:,1),vvvo_pdm%elm4(:,:,:,3),&
                                             & ccsdpt_singles,ccsdpt_doubles,ccsdpt_doubles_2,trip_tmp,trip_ampl)
     
                     case(2)

                        call trip_generator_case2(i,j,nocc,nvirt,ccsd_doubles%val(:,:,i,j),ccsd_doubles%val(:,:,j,i),&
                                                & ccsd_doubles%val(:,:,j,j),ccsd_doubles_portions_i,&
                                                & ccsd_doubles_portions_j,&
                                                & vvvo_pdm%elm4(:,:,:,1),vvvo_pdm%elm4(:,:,:,2),&
                                                & ovoo%val(:,:,i,j),ovoo%val(:,:,j,i),ovoo%val(:,:,j,j),&
                                                & trip_tmp,trip_ampl)

                        ! generate triples amplitudes from trip arrays
     
                        call trip_denom(i,j,j,nocc,nvirt,eivalocc,eivalvirt,trip_ampl)
     
                        ! now do the contractions
     
                        call ccsdpt_driver_case2(i,j,nocc,nvirt,vvoo%val(:,:,i,j),vvoo%val(:,:,j,i),vvoo%val(:,:,j,j),&
                                             & ovoo%val(:,:,i,j),ovoo%val(:,:,j,i),ovoo%val(:,:,j,j),&
                                             & vvvo_pdm%elm4(:,:,:,1),vvvo_pdm%elm4(:,:,:,2),&
                                             & ccsdpt_singles,ccsdpt_doubles,ccsdpt_doubles_2,trip_tmp,trip_ampl)
     
                     case(3)

                        call trip_generator_case3(i,j,k,nocc,nvirt,ccsd_doubles%val(:,:,i,j),ccsd_doubles%val(:,:,i,k),&
                                                & ccsd_doubles%val(:,:,j,i),ccsd_doubles%val(:,:,j,k),&
                                                & ccsd_doubles%val(:,:,k,i),ccsd_doubles%val(:,:,k,j),&
                                                & ccsd_doubles_portions_i,ccsd_doubles_portions_j,&
                                                & ccsd_doubles_portions_k,vvvo_pdm%elm4(:,:,:,1),&
                                                & vvvo_pdm%elm4(:,:,:,2),vvvo_pdm%elm4(:,:,:,3),&
                                                & ovoo%val(:,:,i,j),ovoo%val(:,:,i,k),ovoo%val(:,:,j,i),&
                                                & ovoo%val(:,:,j,k),ovoo%val(:,:,k,i),ovoo%val(:,:,k,j),&
                                                & trip_tmp,trip_ampl)

                        ! generate triples amplitudes from trip arrays
     
                        call trip_denom(i,j,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl)
     
                        ! now do the contractions

                        call ccsdpt_driver_case3(i,j,k,nocc,nvirt,vvoo%val(:,:,i,j),vvoo%val(:,:,i,k),&
                                             & vvoo%val(:,:,j,i),vvoo%val(:,:,j,k),vvoo%val(:,:,k,i),&
                                             & vvoo%val(:,:,k,j),ovoo%val(:,:,i,j),ovoo%val(:,:,i,k),&
                                             & ovoo%val(:,:,j,i),ovoo%val(:,:,j,k),ovoo%val(:,:,k,i),ovoo%val(:,:,k,j),&
                                             & vvvo_pdm%elm4(:,:,:,1),vvvo_pdm%elm4(:,:,:,2),vvvo_pdm%elm4(:,:,:,3),&
                                             & ccsdpt_singles,ccsdpt_doubles,ccsdpt_doubles_2,trip_tmp,trip_ampl)
     
                     end select TypeOfTuple_par
  
                  end do krun_par
  
            end do ijrun_par

    ! release ccsd_doubles_help_arrays
    call mem_dealloc(ccsd_doubles_portions_i)
    call mem_dealloc(ccsd_doubles_portions_j)
    call mem_dealloc(ccsd_doubles_portions_k)

    ! release pdm work array and job list
    call array_free(vvvo_pdm)
    call mem_dealloc(jobs)

    ! release triples ampl structures
    call mem_dealloc(trip_ampl)
    call mem_dealloc(trip_tmp)

  end subroutine ijk_loop_par


  !> \brief: main ijk-loop (serial version)
  !> \author: Janus Juul Eriksen
  !> \date: january 2014
  subroutine ijk_loop_ser(nocc,nvirt,ovoo,vvoo,vvvo,ccsd_doubles,&
                        & ccsdpt_doubles,ccsdpt_doubles_2,ccsdpt_singles,eivalocc,eivalvirt)

    implicit none

    !> nocc,nvirt
    integer, intent(in)      :: nocc,nvirt
    !> 2-el integrals
    type(array4), intent(inout) :: ovoo ! integrals (AI|JK) in the order (J,A,I,K)
    type(array4), intent(inout) :: vvoo ! integrals (AI|BJ) in the order (A,B,I,J)
    type(array), intent(inout)  :: vvvo ! integrals (AI|BC) in the order (C,B,A,I)
    !> ccsd doubles amplitudes
    type(array4), intent(inout) :: ccsd_doubles
    ! o*v^2 portions of ccsd_doubles
    real(realk), pointer, dimension(:,:,:) :: ccsd_doubles_portions_i,ccsd_doubles_portions_j,ccsd_doubles_portions_k
    !> triples amplitudes and 3d work array
    real(realk), pointer, dimension(:,:,:) :: trip_tmp, trip_ampl
    !> ccsd(t) intermediates
    real(realk), dimension(nvirt,nvirt,nocc,nocc) :: ccsdpt_doubles
    real(realk), dimension(nocc,nvirt,nvirt,nocc) :: ccsdpt_doubles_2
    real(realk), dimension(nvirt,nocc) :: ccsdpt_singles
    !> orbital energiesi
    real(realk), intent(inout)  :: eivalocc(nocc), eivalvirt(nvirt)
    !> loop integers
    integer :: i,j,k,idx,tuple_type

    ! init ccsd_doubles_help_arrays
    call mem_alloc(ccsd_doubles_portions_i,nocc,nvirt,nvirt)
    call mem_alloc(ccsd_doubles_portions_j,nocc,nvirt,nvirt)
    call mem_alloc(ccsd_doubles_portions_k,nocc,nvirt,nvirt)

    ! init triples tuples structure
    call mem_alloc(trip_ampl,nvirt,nvirt,nvirt)
    ! init 3d wrk array
    call mem_alloc(trip_tmp,nvirt,nvirt,nvirt)

    ! a note on the openacc implementation.
    ! we use pointers to alias specific memory blocks of the vvoo and ovoo integrals.
    ! this is because we thus avoid to allocate more memory on the device for these integral blocks than what is
    ! actually needed within one loop cycle AND we avoid doing explicit copies into memory buffers on the host.

 irun_ser: do i=1,nocc

          ! store portion of ccsd_doubles (the i'th index) to avoid unnecessary reorderings
          call array_reorder_3d(1.0E0_realk,ccsd_doubles%val(:,:,:,i),nvirt,nvirt,&
                  & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_i)

    jrun_ser: do j=1,i

             ! store portion of ccsd_doubles (the j'th index) to avoid unnecessary reorderings
             call array_reorder_3d(1.0E0_realk,ccsd_doubles%val(:,:,:,j),nvirt,nvirt,&
                     & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_j)

       krun_ser: do k=1,j

                    ! store portion of ccsd_doubles (the k'th index) to avoid unnecessary reorderings
                    call array_reorder_3d(1.0E0_realk,ccsd_doubles%val(:,:,:,k),nvirt,nvirt,&
                            & nocc,[3,2,1],0.0E0_realk,ccsd_doubles_portions_k)
    
                    ! select type of tuple
                    tuple_type = -1
                    ! i == j == k
                    ! this always gives zero contribution
                    ! i == j > k
                    if ((i .eq. j) .and. (j .gt. k) .and. (i .gt. k)) tuple_type = 1
                    ! i > j == k
                    if ((i .gt. j) .and. (j .eq. k) .and. (i .gt. k)) tuple_type = 2
                    ! i > j > k
                    if ((i .gt. j) .and. (j .gt. k) .and. (i .gt. k)) tuple_type = 3
    
                    ! generate tuple(s)
                    TypeOfTuple_ser: select case(tuple_type)
    
                    case(1)

                       call trip_generator_case1(i,k,nocc,nvirt,ccsd_doubles%val(:,:,i,i),ccsd_doubles%val(:,:,i,k),&
                                               & ccsd_doubles%val(:,:,k,i),ccsd_doubles_portions_i,&
                                               & ccsd_doubles_portions_k,&
                                               & vvvo%elm4(:,:,:,i),vvvo%elm4(:,:,:,k),&
                                               & ovoo%val(:,:,i,i),ovoo%val(:,:,i,k),ovoo%val(:,:,k,i),&
                                               & trip_tmp,trip_ampl)
 
                       ! generate triples amplitudes from trip arrays
    
                       call trip_denom(i,i,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl)
    
                       ! now do the contractions
    
                       call ccsdpt_driver_case1(i,k,nocc,nvirt,vvoo%val(:,:,i,i),vvoo%val(:,:,i,k),vvoo%val(:,:,k,i),&
                                            & ovoo%val(:,:,i,i),ovoo%val(:,:,i,k),ovoo%val(:,:,k,i),&
                                            & vvvo%elm4(:,:,:,i),vvvo%elm4(:,:,:,k),&
                                            & ccsdpt_singles,ccsdpt_doubles,ccsdpt_doubles_2,trip_tmp,trip_ampl)
    
                    case(2)

                       call trip_generator_case2(i,j,nocc,nvirt,ccsd_doubles%val(:,:,i,j),ccsd_doubles%val(:,:,j,i),&
                                               & ccsd_doubles%val(:,:,j,j),ccsd_doubles_portions_i,&
                                               & ccsd_doubles_portions_j,&
                                               & vvvo%elm4(:,:,:,i),vvvo%elm4(:,:,:,j),&
                                               & ovoo%val(:,:,i,j),ovoo%val(:,:,j,i),ovoo%val(:,:,j,j),&
                                               & trip_tmp,trip_ampl)

                       ! generate triples amplitudes from trip arrays
    
                       call trip_denom(i,j,j,nocc,nvirt,eivalocc,eivalvirt,trip_ampl)
    
                       ! now do the contractions

                       call ccsdpt_driver_case2(i,j,nocc,nvirt,vvoo%val(:,:,i,j),vvoo%val(:,:,j,i),vvoo%val(:,:,j,j),&
                                            & ovoo%val(:,:,i,j),ovoo%val(:,:,j,i),ovoo%val(:,:,j,j),&
                                            & vvvo%elm4(:,:,:,i),vvvo%elm4(:,:,:,j),&
                                            & ccsdpt_singles,ccsdpt_doubles,ccsdpt_doubles_2,trip_tmp,trip_ampl)
    
                    case(3)

                       call trip_generator_case3(i,j,k,nocc,nvirt,ccsd_doubles%val(:,:,i,j),ccsd_doubles%val(:,:,i,k),&
                                               & ccsd_doubles%val(:,:,j,i),ccsd_doubles%val(:,:,j,k),&
                                               & ccsd_doubles%val(:,:,k,i),ccsd_doubles%val(:,:,k,j),&
                                               & ccsd_doubles_portions_i,ccsd_doubles_portions_j,&
                                               & ccsd_doubles_portions_k,vvvo%elm4(:,:,:,i),&
                                               & vvvo%elm4(:,:,:,j),vvvo%elm4(:,:,:,k),&
                                               & ovoo%val(:,:,i,j),ovoo%val(:,:,i,k),ovoo%val(:,:,j,i),&
                                               & ovoo%val(:,:,j,k),ovoo%val(:,:,k,i),ovoo%val(:,:,k,j),&
                                               & trip_tmp,trip_ampl)

                       ! generate triples amplitudes from trip arrays
    
                       call trip_denom(i,j,k,nocc,nvirt,eivalocc,eivalvirt,trip_ampl)
    
                       ! now do the contractions

                       call ccsdpt_driver_case3(i,j,k,nocc,nvirt,vvoo%val(:,:,i,j),vvoo%val(:,:,i,k),vvoo%val(:,:,j,i),&
                                            & vvoo%val(:,:,j,k),vvoo%val(:,:,k,i),vvoo%val(:,:,k,j),ovoo%val(:,:,i,j),&
                                            & ovoo%val(:,:,i,k),ovoo%val(:,:,j,i),ovoo%val(:,:,j,k),ovoo%val(:,:,k,i),&
                                            & ovoo%val(:,:,k,j),vvvo%elm4(:,:,:,i),vvvo%elm4(:,:,:,j),vvvo%elm4(:,:,:,k),&
                                            & ccsdpt_singles,ccsdpt_doubles,ccsdpt_doubles_2,trip_tmp,trip_ampl)

                    end select TypeOfTuple_ser

                 end do krun_ser
    
              end do jrun_ser
    
           end do irun_ser

    ! release ccsd_doubles_help_arrays
    call mem_dealloc(ccsd_doubles_portions_i)
    call mem_dealloc(ccsd_doubles_portions_j)
    call mem_dealloc(ccsd_doubles_portions_k)

    ! release triples ampl structures
    call mem_dealloc(trip_ampl)
    call mem_dealloc(trip_tmp)

  end subroutine ijk_loop_ser


  !> \brief: create ij_array for ccsd(t)
  !> \author: Janus Juul Eriksen
  !> \date: july 2013
  subroutine create_ij_array_ccsdpt(njobs,no,ij_array)

    implicit none

    !> njobs and nocc
    integer, intent(in) :: njobs,no
    !> ij_array
    integer, dimension(njobs), intent(inout) :: ij_array
    !> integers
    integer :: counter,offset,fill_1,fill_2

    ! since i .ge. j, the composite ij indices will make up a lower triangular matrix.
    ! for each ij, k (where j .ge. k) jobs have to be carried out.
    ! thus, the largest jobs for a given i-value will be those that have the largest j-value,
    ! and the largest jobs will thus be those for which the ij index appears near the diagonal.
    ! as the value of j specifies how large a given job is, we fill up the ij_array with jobs
    ! for j-values in descending order.

    ! the below is the lower triangular part of the ij (5*5) matrix written in row-major order

    ! ||  1               ||
    ! ||  2   3           ||
    ! ||  4   5  6        ||
    ! ||  7   8  9 10     ||
    ! ||  11 12 13 14 15  ||

    ! examples of ij --> i,j conversion
    ! - ij index 15 corresponds to (i,j)=(5,5) and thus to k=1,2,3,4,5
    ! - ij index 9  corresponds to (i,j)=(4,3) and thus to k=1,2,3
    ! - ij index 11  corresponds to (i,j)=(5,1) and thus to k=1

    ! we want ij_array to look like this
    ! (15 , 14 , 10 , 13 , 9 , 6 , 12 , 8 , 5 , 3 , 11 , 7 , 4 , 2 , 1)

    ! counter specifies the index of ij_array
    counter = 1

    do fill_1 = 0,no-1

       ! zero the offset
       offset = 0

       if (fill_1 .eq. 0) then

          ! this is largest possible job, i.e., the (no,no)-th entry in the ij matrix
          ij_array(counter) = njobs
          ! increment counter
          counter = counter + 1

       else

          do fill_2 = 0,fill_1

             if (fill_2 .eq. 0) then

                ! this is the largest i-value, for which we have to do k number of jobs, 
                ! that is, we are at the no'th row essentially moving from right towards left.
                ij_array(counter) = njobs - fill_1
                ! increment counter
                counter = counter + 1

             else

                ! we loop through the i-values keeping the j-value (and k-range) fixed
                ! we thus loop from i == no up towards the diagonal of the lower triangular matrix
                offset = offset + (no - fill_2)
                ! 'njobs - fill_1' gives the current column, while 'offset' moves us up through the rows
                ! while staying below or on the diagonal.(still row-major numbering)
                ij_array(counter) = njobs - fill_1 - offset
                ! increment counter
                counter = counter + 1

             end if

          end do

       end if

    end do

  end subroutine create_ij_array_ccsdpt


  !> \brief: make job distribution list for ccsd(t)
  !> \author: Janus Juul Eriksen
  !> \date: july 2013
  subroutine job_distrib_ccsdpt(b_size,njobs,ij_array,jobs)

    implicit none

    !> batch size (without remainder contribution) and njobs 
    integer, intent(in) :: b_size,njobs
    !> ij_array
    integer, dimension(njobs), intent(inout) :: ij_array
    !> jobs array
    integer, dimension(b_size+1), intent(inout) :: jobs
    !> integers
    integer :: nodtotal,fill,fill_sum

#ifdef VAR_MPI

    nodtotal = infpar%lg_nodtot

    ! fill the jobs array with values of ij stored in ij_array
    ! there are (nocc**2 + nocc)/2 jobs in total (njobs)

    ! the below algorithm distributes the jobs evenly among the nodes.

    do fill = 0,b_size

       fill_sum = infpar%lg_mynum + 1 + fill*nodtotal

       if (fill_sum .le. njobs) then

          jobs(fill + 1) = ij_array(infpar%lg_mynum + 1 + fill*nodtotal) 

       else

          ! fill jobs array with negative number such that this number won't appear for any value of ij
          jobs(fill + 1) = -1

       end if

    end do

#endif

  end subroutine job_distrib_ccsdpt


  !> \brief: determine i and j from ij
  !> \author: Janus Juul Eriksen
  !> \date: july 2013
  subroutine calc_i_and_j(ij,no,i,j)

    implicit none

    !> composite ij index
    integer, intent(in) :: ij,no
    !> i and j
    integer, intent(inout) :: i,j
    !> integers
    integer :: gauss_sum,gauss_sum_old,series

    ! in a N x N lower triangular matrix, there is a total of (N**2 + N)/2 elements.
    ! in column 1, there are N rows, in column 2, there are (N-1) rows, ...,  in
    ! column (N-1), there are 2 rows, and in column N, there are 1 row.
    ! this is a gauss sum of 1 + 2 + 3 + ... + N-2 + N-1 + N
    ! for a given value of i, the value of ij can thus at max be (i**2 + i)/2 (gauss_sum).
    ! if gauss_sum for a given i (series) is smaller than ij, we loop.
    ! when gauss_sum is greater than ij, we use the value of i for the present loop cycle
    ! and calculate the value of j from the present ij and previous gauss_sum values.
    ! when gauss_sum is equal to ij, we are on the diagonal and i == j (== series).

    do series = 1,no

       gauss_sum = int((series**2 + series)/2)

       if (gauss_sum .lt. ij) then

          gauss_sum_old = gauss_sum

          cycle

       else if (gauss_sum .eq. ij) then

          j = series
          i = series

          exit

       else

          j = ij - gauss_sum_old
          i = series

          exit

       end if

    end do

  end subroutine calc_i_and_j


  !> \brief: generator for triples amplitudes, case(1)
  !> \author: Janus Juul Eriksen
  !> \date: february 2014
  subroutine trip_generator_case1(oindex1,oindex3,no,nv,ccsd_doubles_11,ccsd_doubles_13,&
                                & ccsd_doubles_31,ccsd_doubles_portions_1,ccsd_doubles_portions_3,&
                                & cbai_tile_1,cbai_tile_3,jaik_tile_11,&
                                & jaik_tile_13,jaik_tile_31,trip_tmp,trip_ampl)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex3,no,nv
    !> nv**2 tiles of ccsd_doubles
    real(realk), dimension(nv,nv) :: ccsd_doubles_11, ccsd_doubles_13, ccsd_doubles_31
    !> no*nv**2 tiles of ccsd_doubles
    real(realk), dimension(no,nv,nv) :: ccsd_doubles_portions_1,ccsd_doubles_portions_3
    !> tiles of jaik 2-el integrals
    real(realk), dimension(no,nv) :: jaik_tile_11, jaik_tile_13, jaik_tile_31
    !> tiles of cbai 2-el integrals
    real(realk), dimension(nv,nv,nv), intent(inout) :: cbai_tile_1, cbai_tile_3
    !> triples amplitude and work array
    real(realk), dimension(nv,nv,nv) :: trip_tmp, trip_ampl

    ! iik,iki
    call trip_amplitudes_virt(oindex1,oindex1,oindex3,no,nv,ccsd_doubles_11,&
                            & cbai_tile_3,trip_tmp)
    call trip_amplitudes_occ(oindex1,oindex3,oindex1,no,nv,ccsd_doubles_portions_1,&
                            & jaik_tile_13,trip_tmp)
    trip_ampl = trip_tmp
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,1,3],1.0E0_realk,trip_ampl)

    ! kii,iik
    call trip_amplitudes_virt(oindex3,oindex1,oindex1,no,nv,ccsd_doubles_13,&
                            & cbai_tile_1,trip_tmp)
    call trip_amplitudes_occ(oindex1,oindex1,oindex3,no,nv,ccsd_doubles_portions_1,&
                            & jaik_tile_31,trip_tmp)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,3,1],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,2,1],1.0E0_realk,trip_ampl)

    ! iki.kii
    call trip_amplitudes_virt(oindex1,oindex3,oindex1,no,nv,ccsd_doubles_31,&
                            & cbai_tile_1,trip_tmp)
    call trip_amplitudes_occ(oindex3,oindex1,oindex1,no,nv,ccsd_doubles_portions_3,&
                            & jaik_tile_11,trip_tmp)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [1,3,2],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,1,2],1.0E0_realk,trip_ampl)

  end subroutine trip_generator_case1


  !> \brief: generator for triples amplitudes, case(2)
  !> \author: Janus Juul Eriksen
  !> \date: february 2014
  subroutine trip_generator_case2(oindex1,oindex2,no,nv,ccsd_doubles_12,ccsd_doubles_21,&
                                & ccsd_doubles_22,ccsd_doubles_portions_1,ccsd_doubles_portions_2,&
                                & cbai_tile_1,cbai_tile_2,jaik_tile_12,&
                                & jaik_tile_21,jaik_tile_22,trip_tmp,trip_ampl)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex2,no,nv
    !> nv**2 tiles of ccsd_doubles
    real(realk), dimension(nv,nv) :: ccsd_doubles_12, ccsd_doubles_21, ccsd_doubles_22
    !> no*nv**2 tiles of ccsd_doubles
    real(realk), dimension(no,nv,nv) :: ccsd_doubles_portions_1,ccsd_doubles_portions_2
    !> tiles of jaik 2-el integrals
    real(realk), dimension(no,nv) :: jaik_tile_12, jaik_tile_21, jaik_tile_22
    !> tiles of cbai 2-el integrals
    real(realk), dimension(nv,nv,nv), intent(inout) :: cbai_tile_1, cbai_tile_2
    !> triples amplitude and work array
    real(realk), dimension(nv,nv,nv) :: trip_tmp, trip_ampl

    ! ijj.jji
    call trip_amplitudes_virt(oindex1,oindex2,oindex2,no,nv,ccsd_doubles_21,&
                            & cbai_tile_2,trip_tmp)
    call trip_amplitudes_occ(oindex2,oindex2,oindex1,no,nv,ccsd_doubles_portions_2,&
                            & jaik_tile_12,trip_tmp)
    trip_ampl = trip_tmp
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [1,3,2],1.0E0_realk,trip_ampl)
 
    ! jij,ijj
    call trip_amplitudes_virt(oindex2,oindex1,oindex2,no,nv,ccsd_doubles_12,&
                            & cbai_tile_2,trip_tmp)
    call trip_amplitudes_occ(oindex1,oindex2,oindex2,no,nv,ccsd_doubles_portions_1,&
                            & jaik_tile_22,trip_tmp)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,1,3],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,3,1],1.0E0_realk,trip_ampl)
 
    ! jji,jij
    call trip_amplitudes_virt(oindex2,oindex2,oindex1,no,nv,ccsd_doubles_22,&
                            & cbai_tile_1,trip_tmp)
    call trip_amplitudes_occ(oindex2,oindex1,oindex2,no,nv,ccsd_doubles_portions_2,&
                            & jaik_tile_21,trip_tmp)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,1,2],1.0E0_realk,trip_ampl)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,2,1],1.0E0_realk,trip_ampl)

  end subroutine trip_generator_case2


  !> \brief: generator for triples amplitudes, case(3)
  !> \author: Janus Juul Eriksen
  !> \date: february 2014
  subroutine trip_generator_case3(oindex1,oindex2,oindex3,no,nv,ccsd_doubles_12,ccsd_doubles_13,&
                                & ccsd_doubles_21,ccsd_doubles_23,ccsd_doubles_31,ccsd_doubles_32,&
                                & ccsd_doubles_portions_1,ccsd_doubles_portions_2,&
                                & ccsd_doubles_portions_3,cbai_tile_1,cbai_tile_2,cbai_tile_3,&
                                & jaik_tile_12,jaik_tile_13,jaik_tile_21,jaik_tile_23,jaik_tile_31,&
                                & jaik_tile_32,trip_tmp,trip_ampl)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex2,oindex3,no,nv
    !> nv**2 tiles of ccsd_doubles
    real(realk), dimension(nv,nv) :: ccsd_doubles_12, ccsd_doubles_13, ccsd_doubles_21
    real(realk), dimension(nv,nv) :: ccsd_doubles_23, ccsd_doubles_31, ccsd_doubles_32
    !> no*nv**2 tiles of ccsd_doubles
    real(realk), dimension(no,nv,nv) :: ccsd_doubles_portions_1,ccsd_doubles_portions_2,ccsd_doubles_portions_3
    !> tiles of jaik 2-el integrals
    real(realk), dimension(no,nv) :: jaik_tile_12, jaik_tile_13, jaik_tile_21
    real(realk), dimension(no,nv) :: jaik_tile_23, jaik_tile_31, jaik_tile_32
    !> tiles of cbai 2-el integrals 
    real(realk), dimension(nv,nv,nv), intent(inout) :: cbai_tile_1, cbai_tile_2, cbai_tile_3
    !> triples amplitude and work array
    real(realk), dimension(nv,nv,nv) :: trip_tmp, trip_ampl

    ! ijk.jki
    call trip_amplitudes_virt(oindex1,oindex2,oindex3,no,nv,ccsd_doubles_21,&
                            & cbai_tile_3,trip_ampl)
    call trip_amplitudes_occ(oindex2,oindex3,oindex1,no,nv,ccsd_doubles_portions_2,&
                            & jaik_tile_13,trip_ampl)

    ! kij,ijk
    call trip_amplitudes_virt(oindex3,oindex1,oindex2,no,nv,ccsd_doubles_13,&
                            & cbai_tile_2,trip_tmp)
    call trip_amplitudes_occ(oindex1,oindex2,oindex3,no,nv,ccsd_doubles_portions_1,&
                            & jaik_tile_32,trip_tmp)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,3,1],1.0E0_realk,trip_ampl)

    ! jki,kij
    call trip_amplitudes_virt(oindex2,oindex3,oindex1,no,nv,ccsd_doubles_32,&
                            & cbai_tile_1,trip_tmp)
    call trip_amplitudes_occ(oindex3,oindex1,oindex2,no,nv,ccsd_doubles_portions_3,&
                            & jaik_tile_21,trip_tmp)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,1,2],1.0E0_realk,trip_ampl)

    ! ikj,kji
    call trip_amplitudes_virt(oindex1,oindex3,oindex2,no,nv,ccsd_doubles_31,&
                            & cbai_tile_2,trip_tmp)
    call trip_amplitudes_occ(oindex3,oindex2,oindex1,no,nv,ccsd_doubles_portions_3,&
                            & jaik_tile_12,trip_tmp)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [1,3,2],1.0E0_realk,trip_ampl)
 
    ! jik,ikj
    call trip_amplitudes_virt(oindex2,oindex1,oindex3,no,nv,ccsd_doubles_12,&
                            & cbai_tile_3,trip_tmp)
    call trip_amplitudes_occ(oindex1,oindex3,oindex2,no,nv,ccsd_doubles_portions_1,&
                            & jaik_tile_23,trip_tmp)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [2,1,3],1.0E0_realk,trip_ampl)
 
    ! kji,jik
    call trip_amplitudes_virt(oindex3,oindex2,oindex1,no,nv,ccsd_doubles_23,&
                            & cbai_tile_1,trip_tmp)
    call trip_amplitudes_occ(oindex2,oindex1,oindex3,no,nv,ccsd_doubles_portions_2,&
                            & jaik_tile_31,trip_tmp)
    call array_reorder_3d(1.0E0_realk,trip_tmp,nv,nv,nv,&
                        & [3,2,1],1.0E0_realk,trip_ampl)

  end subroutine trip_generator_case3


  !> \brief: driver routine for contractions in case(1) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: march 2013
  subroutine ccsdpt_driver_case1(oindex1,oindex3,no,nv,abij_tile_12,abij_tile_13,abij_tile_31,&
                            & jaik_tile_12,jaik_tile_13,jaik_tile_31,&
                            & int_virt_tile_o1,int_virt_tile_o3,&
                            & ccsdpt_singles,ccsdpt_doubles,ccsdpt_doubles_2,wrk_3d,trip)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex3,no,nv
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(nv,nv,no,no) :: ccsdpt_doubles
    real(realk), dimension(no,nv,nv,no) :: ccsdpt_doubles_2
    real(realk), dimension(nv,no) :: ccsdpt_singles
    !> aibj and jaik 2-el integrals
    !> tiles of jaik integrals (AI|JK) in the order (J,A,I,K)
    real(realk), dimension(no,nv) :: jaik_tile_12, jaik_tile_13, jaik_tile_31
    !> tiles of abij integrals (AI|BJ) in the order (A,B,I,J)
    real(realk), dimension(nv,nv) :: abij_tile_12, abij_tile_13, abij_tile_31
    !> tiles of cbai 2-el integrals determined by incomming oindex1,oindex3
    real(realk), dimension(nv,nv,nv), intent(inout) :: int_virt_tile_o1
    real(realk), dimension(nv,nv,nv), intent(inout) :: int_virt_tile_o3
    !> triples amplitude and work array
    real(realk), dimension(nv,nv,nv) :: trip, wrk_3d
    !> loop integer
    integer :: idx

    ! before the calls to the contractions in ccsdpt_contract_211/212 and ccsdpt_contract_221/222,
    ! we implicitly do a [2,3,1] reordering. in order to minimize the number of reorderings needed to be
    ! performed, and in order to take optimal advantage of the symmetry of the amplitudes, we carry out
    ! the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : iik --312--> kii --312--> iki
    ! in 211/212 : kii ........ iki ........ iik
    ! in 221/222 : iki ........ iik ........ kii

    do idx = 1,3

       if (idx .eq. 1) then

          ! calculate contribution to ccsdpt_singles:

          call ccsdpt_contract_11(oindex1,oindex1,oindex3,nv,no,abij_tile_13,abij_tile_31,ccsdpt_singles,&
                       & trip,.false.)
          call ccsdpt_contract_12(oindex1,oindex1,oindex3,nv,no,abij_tile_12,abij_tile_12,ccsdpt_singles,&
                       & trip,.false.)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_211(oindex3,oindex1,oindex1,nv,no,&
                           & ccsdpt_doubles,wrk_3d,trip,int_virt_tile_o1,.false.)
          call ccsdpt_contract_212(oindex3,oindex1,oindex1,nv,no,&
                           & ccsdpt_doubles,wrk_3d,trip,int_virt_tile_o3)

          ! now do occ part:

          call ccsdpt_contract_221(oindex1,oindex3,oindex1,no,nv,jaik_tile_31,jaik_tile_13,&
                           & ccsdpt_doubles_2,trip,.true.)

       else if (idx .eq. 2) then

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array
         call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,wrk_3d)

          call ccsdpt_contract_11(oindex3,oindex1,oindex1,nv,no,abij_tile_12,abij_tile_12,ccsdpt_singles,&
                       & wrk_3d,.true.)
          call ccsdpt_contract_12(oindex3,oindex1,oindex1,nv,no,abij_tile_13,abij_tile_31,ccsdpt_singles,&
                       & wrk_3d,.true.)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_211(oindex1,oindex3,oindex1,nv,no,&
                           & ccsdpt_doubles,trip,wrk_3d,int_virt_tile_o1,.true.)

          ! now do occ part:

          call ccsdpt_contract_221(oindex1,oindex1,oindex3,no,nv,jaik_tile_13,jaik_tile_31,&
                           & ccsdpt_doubles_2,wrk_3d,.false.)
          call ccsdpt_contract_222(oindex1,oindex1,oindex3,no,nv,jaik_tile_12,&
                           & ccsdpt_doubles_2,wrk_3d)

       else if (idx .eq. 3) then

          ! iki: this case is redundant since both the coulumb and the exchange contributions
          ! will be contructed from the ampl_iki trip amplitudes and therefore end up
          ! canceling each other when added to ccsdpt_singles

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array
          call array_reorder_3d(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,trip)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_211(oindex1,oindex1,oindex3,nv,no,&
                           & ccsdpt_doubles,wrk_3d,trip,int_virt_tile_o3,.false.)
          call ccsdpt_contract_212(oindex1,oindex1,oindex3,nv,no,&
                           & ccsdpt_doubles,wrk_3d,trip,int_virt_tile_o1)

          ! now do occ part:

          call ccsdpt_contract_221(oindex3,oindex1,oindex1,no,nv,jaik_tile_12,jaik_tile_12,&
                           & ccsdpt_doubles_2,trip,.false.)
          call ccsdpt_contract_222(oindex3,oindex1,oindex1,no,nv,jaik_tile_31,&
                           & ccsdpt_doubles_2,trip)

       end if

    end do

  end subroutine ccsdpt_driver_case1


  !> \brief: driver routine for contractions in case(2) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: march 2013
  subroutine ccsdpt_driver_case2(oindex1,oindex2,no,nv,abij_tile_12,abij_tile_21,abij_tile_23,&
                            & jaik_tile_12,jaik_tile_21,jaik_tile_23,&
                            & int_virt_tile_o1,int_virt_tile_o2,&
                            & ccsdpt_singles,ccsdpt_doubles,ccsdpt_doubles_2,wrk_3d,trip)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex2,no,nv
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(nv,nv,no,no) :: ccsdpt_doubles
    real(realk), dimension(no,nv,nv,no) :: ccsdpt_doubles_2
    real(realk), dimension(nv,no) :: ccsdpt_singles
    !> aibj and jaik 2-el integrals
    !> tiles of jaik integrals (AI|JK) in the order (J,A,I,K)
    real(realk), dimension(no,nv) :: jaik_tile_12, jaik_tile_21, jaik_tile_23
    !> tiles of abij integrals (AI|BJ) in the order (A,B,I,J)
    real(realk), dimension(nv,nv) :: abij_tile_12, abij_tile_21, abij_tile_23
    !> tiles of cbai 2-el integrals determined by incomming oindex1,oindex2
    real(realk), dimension(nv,nv,nv), intent(inout) :: int_virt_tile_o1
    real(realk), dimension(nv,nv,nv), intent(inout) :: int_virt_tile_o2
    !> triples amplitude and work array
    real(realk), dimension(nv,nv,nv) :: trip, wrk_3d
    !> loop integer
    integer :: idx

    ! before the calls to the contractions in ccsdpt_contract_211/212 and ccsdpt_contract_221/222,
    ! we implicitly do a [2,3,1] reordering. in order to minimize the number of reorderings needed to be
    ! performed, and in order to take optimal advantage of the symmetry of the amplitudes, we carry out
    ! the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : ijj --312--> jij --312--> jji
    ! in 211/212 : jij ........ jji ........ ijj
    ! in 221/222 : jji ........ ijj ........ jij

    do idx = 1,3

       if (idx .eq. 1) then
  
          ! calculate contributions to ccsdpt_singles:
 
          call ccsdpt_contract_11(oindex1,oindex2,oindex2,nv,no,abij_tile_23,abij_tile_23,ccsdpt_singles,&
                       & trip,.true.)
          call ccsdpt_contract_12(oindex1,oindex2,oindex2,nv,no,abij_tile_21,abij_tile_12,ccsdpt_singles,&
                       & trip,.true.)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_211(oindex2,oindex1,oindex2,nv,no,&
                           & ccsdpt_doubles,wrk_3d,trip,int_virt_tile_o2,.true.)

          ! now do occ part:

          call ccsdpt_contract_221(oindex2,oindex2,oindex1,no,nv,jaik_tile_21,jaik_tile_12,&
                           & ccsdpt_doubles_2,trip,.false.)
          call ccsdpt_contract_222(oindex2,oindex2,oindex1,no,nv,jaik_tile_23,&
                           & ccsdpt_doubles_2,trip)

       else if (idx .eq. 2) then
   
          ! this case is redundant since both the coulumb and the exchange contributions
          ! will be contructed from the ampl_jij trip amplitudes and therefore end up
          ! canceling each other when added to T_star

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,wrk_3d)
 
          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_211(oindex2,oindex2,oindex1,nv,no,&
                           & ccsdpt_doubles,trip,wrk_3d,int_virt_tile_o1,.false.)
          call ccsdpt_contract_212(oindex2,oindex2,oindex1,nv,no,&
                           & ccsdpt_doubles,trip,wrk_3d,int_virt_tile_o2)

          ! now do occ part:

          call ccsdpt_contract_221(oindex1,oindex2,oindex2,no,nv,jaik_tile_23,jaik_tile_23,&
                           & ccsdpt_doubles_2,wrk_3d,.false.)
          call ccsdpt_contract_222(oindex1,oindex2,oindex2,no,nv,jaik_tile_12,&
                           & ccsdpt_doubles_2,wrk_3d)

       else if (idx .eq. 3) then

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array
          call array_reorder_3d(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,trip)

          ! calculate contributions to ccsdpt_singles:
   
          call ccsdpt_contract_11(oindex2,oindex2,oindex1,nv,no,abij_tile_21,abij_tile_12,ccsdpt_singles,&
                       & trip,.false.)
          call ccsdpt_contract_12(oindex2,oindex2,oindex1,nv,no,abij_tile_23,abij_tile_23,ccsdpt_singles,&
                       & trip,.false.)
   
          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_211(oindex1,oindex2,oindex2,nv,no,&
                           & ccsdpt_doubles,wrk_3d,trip,int_virt_tile_o2,.false.)
          call ccsdpt_contract_212(oindex1,oindex2,oindex2,nv,no,&
                           & ccsdpt_doubles,wrk_3d,trip,int_virt_tile_o1)

          ! now do occ part:

          call ccsdpt_contract_221(oindex2,oindex1,oindex2,no,nv,jaik_tile_12,jaik_tile_21,&
                           & ccsdpt_doubles_2,trip,.true.)

       end if

    end do

  end subroutine ccsdpt_driver_case2


  !> \brief: driver routine for contractions in case(3) of ccsdpt_driver
  !> \author: Janus Juul Eriksen
  !> \date: march 2013
  subroutine ccsdpt_driver_case3(oindex1,oindex2,oindex3,no,nv,&
                            & abij_tile_12, abij_tile_13, abij_tile_21,&
                            & abij_tile_23, abij_tile_31, abij_tile_32,&
                            & jaik_tile_12, jaik_tile_13, jaik_tile_21,&
                            & jaik_tile_23, jaik_tile_31, jaik_tile_32,&
                            & int_virt_tile_o1,int_virt_tile_o2,int_virt_tile_o3,&
                            & ccsdpt_singles,ccsdpt_doubles,ccsdpt_doubles_2,wrk_3d,trip)

    implicit none

    !> i, j, k, nocc, and nvirt
    integer, intent(in) :: oindex1,oindex2,oindex3,no,nv
    !> ccsd(t) singles and doubles amplitudes
    real(realk), dimension(nv,nv,no,no) :: ccsdpt_doubles
    real(realk), dimension(no,nv,nv,no) :: ccsdpt_doubles_2
    real(realk), dimension(nv,no) :: ccsdpt_singles
    !> aibj and jaik 2-el integrals
    !> tiles of jaik integrals (AI|JK) in the order (J,A,I,K)
    real(realk), dimension(no,nv) :: jaik_tile_12, jaik_tile_13, jaik_tile_21
    real(realk), dimension(no,nv) :: jaik_tile_23, jaik_tile_31, jaik_tile_32
    !> tiles of abij integrals (AI|BJ) in the order (A,B,I,J)
    real(realk), dimension(nv,nv) :: abij_tile_12, abij_tile_13, abij_tile_21
    real(realk), dimension(nv,nv) :: abij_tile_23, abij_tile_31, abij_tile_32
    !> tiles of cbai 2-el integrals determined by incomming oindex1,oindex2,oindex3
    real(realk), dimension(nv,nv,nv), intent(inout) :: int_virt_tile_o1
    real(realk), dimension(nv,nv,nv), intent(inout) :: int_virt_tile_o2
    real(realk), dimension(nv,nv,nv), intent(inout) :: int_virt_tile_o3
    !> triples amplitude and work array
    real(realk), dimension(nv,nv,nv) :: trip, wrk_3d
    !> loop integer
    integer :: idx

    ! before the calls to the contractions in ccsdpt_contract_211/212 and ccsdpt_contract_221/222,
    ! we implicitly do a [2,3,1] reordering. in order to minimize the number of reorderings needed to be
    ! performed, and in order to take optimal advantage of the symmetry of the amplitudes, we carry out
    ! the amplitudes in accordance to the following scheme
    !
    ! in 11/12   : ijk --312--> kij --312--> jki --213--> kji --312--> ikj --312--> jik
    ! in 211/212 : kij ........ jki ........ ijk ........ ikj ........ jik ........ kji
    ! in 221/222 : jki ........ ijk ........ kij ........ jik ........ kji ........ ikj

    do idx = 1,6

       if (idx .eq. 1) then

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_11(oindex1,oindex2,oindex3,nv,no,abij_tile_23,abij_tile_32,ccsdpt_singles,&
                       & trip,.false.)
          call ccsdpt_contract_12(oindex1,oindex2,oindex3,nv,no,abij_tile_21,abij_tile_12,ccsdpt_singles,&
                       & trip,.false.)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_211(oindex3,oindex1,oindex2,nv,no,&
                           & ccsdpt_doubles,wrk_3d,trip,int_virt_tile_o2,.false.)
          call ccsdpt_contract_212(oindex3,oindex1,oindex2,nv,no,&
                           & ccsdpt_doubles,wrk_3d,trip,int_virt_tile_o3)

          ! now do occ part:

          call ccsdpt_contract_221(oindex2,oindex3,oindex1,no,nv,jaik_tile_31,jaik_tile_13,&
                           & ccsdpt_doubles_2,trip,.false.)
          call ccsdpt_contract_222(oindex2,oindex3,oindex1,no,nv,jaik_tile_23,&
                           & ccsdpt_doubles_2,trip)

       else if (idx .eq. 2) then

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,wrk_3d)

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_11(oindex3,oindex1,oindex2,nv,no,abij_tile_12,abij_tile_21,ccsdpt_singles,&
                       & wrk_3d,.false.)
          call ccsdpt_contract_12(oindex3,oindex1,oindex2,nv,no,abij_tile_13,abij_tile_31,ccsdpt_singles,&
                       & wrk_3d,.false.)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_211(oindex2,oindex3,oindex1,nv,no,&
                           & ccsdpt_doubles,trip,wrk_3d,int_virt_tile_o1,.false.)
          call ccsdpt_contract_212(oindex2,oindex3,oindex1,nv,no,&
                           & ccsdpt_doubles,trip,wrk_3d,int_virt_tile_o2)

          ! now do occ part:

          call ccsdpt_contract_221(oindex1,oindex2,oindex3,no,nv,jaik_tile_23,jaik_tile_32,&
                           & ccsdpt_doubles_2,wrk_3d,.false.)
          call ccsdpt_contract_222(oindex1,oindex2,oindex3,no,nv,jaik_tile_12,&
                           & ccsdpt_doubles_2,wrk_3d)

       else if (idx .eq. 3) then

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array
          call array_reorder_3d(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,trip)

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_11(oindex2,oindex3,oindex1,nv,no,abij_tile_31,abij_tile_13,ccsdpt_singles,&
                       & trip,.false.)
          call ccsdpt_contract_12(oindex2,oindex3,oindex1,nv,no,abij_tile_32,abij_tile_23,ccsdpt_singles,&
                       & trip,.false.)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_211(oindex1,oindex2,oindex3,nv,no,&
                           & ccsdpt_doubles,wrk_3d,trip,int_virt_tile_o3,.false.)
          call ccsdpt_contract_212(oindex1,oindex2,oindex3,nv,no,&
                           & ccsdpt_doubles,wrk_3d,trip,int_virt_tile_o1)

          ! now do occ part:

          call ccsdpt_contract_221(oindex3,oindex1,oindex2,no,nv,jaik_tile_12,jaik_tile_21,&
                           & ccsdpt_doubles_2,trip,.false.)
          call ccsdpt_contract_222(oindex3,oindex1,oindex2,no,nv,jaik_tile_31,&
                           & ccsdpt_doubles_2,trip)

       else if (idx .eq. 4) then

          !*** special reordering (see note above)
          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[2,1,3],0.0E0_realk,wrk_3d)

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_11(oindex3,oindex2,oindex1,nv,no,abij_tile_21,abij_tile_12,ccsdpt_singles,&
                       & wrk_3d,.false.)
          call ccsdpt_contract_12(oindex3,oindex2,oindex1,nv,no,abij_tile_23,abij_tile_32,ccsdpt_singles,&
                       & wrk_3d,.false.)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_211(oindex1,oindex3,oindex2,nv,no,&
                           & ccsdpt_doubles,trip,wrk_3d,int_virt_tile_o2,.false.)
          call ccsdpt_contract_212(oindex1,oindex3,oindex2,nv,no,&
                           & ccsdpt_doubles,trip,wrk_3d,int_virt_tile_o1)

          ! now do occ part:

          call ccsdpt_contract_221(oindex2,oindex1,oindex3,no,nv,jaik_tile_13,jaik_tile_31,&
                           & ccsdpt_doubles_2,wrk_3d,.false.)
          call ccsdpt_contract_222(oindex2,oindex1,oindex3,no,nv,jaik_tile_21,&
                           & ccsdpt_doubles_2,wrk_3d)

       else if (idx .eq. 5) then

          ! initially, reorder wrk_3d - after reordering, trip holds the triples ampls
          ! and wrk_3d is a 3d work array
          call array_reorder_3d(1.0E0_realk,wrk_3d,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,trip)

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_11(oindex1,oindex3,oindex2,nv,no,abij_tile_32,abij_tile_23,ccsdpt_singles,&
                       & trip,.false.)
          call ccsdpt_contract_12(oindex1,oindex3,oindex2,nv,no,abij_tile_31,abij_tile_13,ccsdpt_singles,&
                       & trip,.false.)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_211(oindex2,oindex1,oindex3,nv,no,&
                           & ccsdpt_doubles,wrk_3d,trip,int_virt_tile_o3,.false.)
          call ccsdpt_contract_212(oindex2,oindex1,oindex3,nv,no,&
                           & ccsdpt_doubles,wrk_3d,trip,int_virt_tile_o2)

          ! now do occ part:

          call ccsdpt_contract_221(oindex3,oindex2,oindex1,no,nv,jaik_tile_21,jaik_tile_12,&
                           & ccsdpt_doubles_2,trip,.false.)
          call ccsdpt_contract_222(oindex3,oindex2,oindex1,no,nv,jaik_tile_32,&
                           & ccsdpt_doubles_2,trip)

       else if (idx .eq. 6) then

          ! initially, reorder trip - after reordering, wrk_3d holds the triples ampls
          ! and trip is a 3d work array
          call array_reorder_3d(1.0E0_realk,trip,nv,nv,&
                           & nv,[3,1,2],0.0E0_realk,wrk_3d)

          ! calculate contributions to ccsdpt_singles:

          call ccsdpt_contract_11(oindex2,oindex1,oindex3,nv,no,abij_tile_13,abij_tile_31,ccsdpt_singles,&
                       & wrk_3d,.false.)
          call ccsdpt_contract_12(oindex2,oindex1,oindex3,nv,no,abij_tile_12,abij_tile_21,ccsdpt_singles,&
                       & wrk_3d,.false.)

          ! calculate contributions to ccsdpt_doubles (virt part):

          call ccsdpt_contract_211(oindex3,oindex2,oindex1,nv,no,&
                           & ccsdpt_doubles,trip,wrk_3d,int_virt_tile_o1,.false.)
          call ccsdpt_contract_212(oindex3,oindex2,oindex1,nv,no,&
                           & ccsdpt_doubles,trip,wrk_3d,int_virt_tile_o3)

          ! now do occ part:

          call ccsdpt_contract_221(oindex1,oindex3,oindex2,no,nv,jaik_tile_32,jaik_tile_23,&
                           & ccsdpt_doubles_2,wrk_3d,.false.)
          call ccsdpt_contract_222(oindex1,oindex3,oindex2,no,nv,jaik_tile_13,&
                           & ccsdpt_doubles_2,wrk_3d)

       end if

    end do

  end subroutine ccsdpt_driver_case3


  !> \brief: transform ccsd doubles from local to canonical basis
  !> \author: Janus Juul Eriksen
  !> \date: september 2012
  !> \param: ccsd_t2, no and nv are nocc and nvirt, respectively, and U_occ and U_virt
  !          are unitary matrices from local --> canonical basis
  subroutine ccsdpt_local_can_trans(ccsd_t2,no,nv,U_occ,U_virt)

    implicit none
    !> ccsd doubles
    type(array4), intent(inout) :: ccsd_t2
    !> unitary transformation matrices
    type(array2), intent(inout) :: U_occ, U_virt
    !> integers
    integer, intent(in) :: no, nv
    !> temp array4 structures
    type(array4) :: tmp1, tmp2

    ! (a,i,b,j) are local basis indices and (A,I,B,J) refer to the canonical basis.
    ! we want to carry out the transformation:
    ! T^{AB}_{IJ} = sum_{aibj} U_{aA} U_{iI} U_{bB} U_{jJ} T^{ab}_{ij}

    ! 1. Init temporary arrays, dims = aiai
    tmp1 = array4_init_standard([nv,no,nv,no])

    ! 2. 1st index: doub_ampl(a,i,b,j) --> tmp1(A,i,b,j)
    call array4_contract1(ccsd_t2,U_virt,tmp1,.true.)
    call array4_free(ccsd_t2)

    ! 3. 2nd index: tmp1(A,i,b,j) --> tmp1(i,A,b,j) --> tmp2(I,A,b,j)
    call array4_reorder(tmp1,[2,1,3,4])
    tmp2 = array4_init_standard([no,nv,nv,no])
    call array4_contract1(tmp1,U_occ,tmp2,.true.)
    call array4_free(tmp1)

    ! 4. 3rd index: tmp2(I,A,b,j) --> tmp2(j,b,A,I) --> tmp1(J,b,A,I)
    call array4_reorder(tmp2,[4,3,2,1])
    tmp1 = array4_init_standard([no,nv,nv,no])
    call array4_contract1(tmp2,U_occ,tmp1,.true.)
    call array4_free(tmp2)

    ! 5. 4th index: tmp1(J,b,A,I) --> tmp1(b,J,A,I) --> doub_ampl(B,J,A,I) = doub_ampl(A,I,B,J)
    call array4_reorder(tmp1,[2,1,3,4])
    ccsd_t2 = array4_init_standard([nv,no,nv,no])
    call array4_contract1(tmp1,U_virt,ccsd_t2,.true.)
    call array4_free(tmp1)

  end subroutine ccsdpt_local_can_trans


  !> \brief: transform ccsd_doubles, ccsdpt_singles and ccsdpt_doubles from canonical to local basis
  !> \author: Janus Juul Eriksen
  !> \date: september 2012
  !> \param: ccsd_t2, ccsdpt_t1, ccsdpt_t2, no and nv are nocc and nvirt, respectively, 
  !<         and U_occ and U_virt are unitary matrices from canonical --> local basis
  subroutine ccsdpt_can_local_trans(ccsd_t2,ccsdpt_t1,ccsdpt_t2,no,nv,U_occ,U_virt)

    implicit none
    !> ccsdpt_singles
    type(array2), intent(inout) :: ccsdpt_t1
    !> ccsd_doubles and ccsdpt_doubles
    type(array4), intent(inout) :: ccsd_t2, ccsdpt_t2
    !> unitary transformation matrices
    type(array2), intent(inout) :: U_occ, U_virt
    !> integers
    integer, intent(in) :: no, nv
    !> temp array2 and array4 structures
    type(array2) :: tmp0
    type(array4) :: tmp1, tmp2, tmp3, tmp4

    ! (a,i,b,j) are local basis indices and (A,I,B,J) refer to the canonical basis.

    ! 1a. init temporary array4s, tmp1 and tmp3
    tmp1 = array4_init_standard([nv,nv,no,no])
    tmp3 = array4_init_standard([nv,nv,no,no])
    ! 1b. init temporary array2, tmp0
    tmp0 = array2_init_plain([nv,no])

    ! 2. 1st index:
    ! ccsdpt_t2(A,B,I,J) --> tmp1(a,B,I,J)
    ! ccsd_t2(B,A,J,I) --> tmp3(b,A,J,I)
    call array4_contract1(ccsdpt_t2,U_virt,tmp1,.true.)
    call array4_contract1(ccsd_t2,U_virt,tmp3,.true.)
    ! ccsdpt_t1(A,I) --> tmp0(a,I)
    call array2_matmul(U_virt,ccsdpt_t1,tmp0,'t','n',1.0E0_realk,0.0E0_realk)

    ! free ccsdpt_doubles and ccsd_doubles
    call array4_free(ccsdpt_t2)
    call array4_free(ccsd_t2)

    ! 3. 2nd index:
    ! tmp1(a,B,I,J) --> tmp1(B,a,I,J) --> tmp2(b,a,I,J)
    ! tmp3(b,A,J,I) --> tmp3(A,b,J,I) --> tmp4(a,b,J,I) 
    ! tmp0(a,I) --> ccsdpt_t1(a,i)
    call array4_reorder(tmp1,[2,1,3,4])
    call array4_reorder(tmp3,[2,1,3,4])

    ! init temporary array4s, tmp2 and tmp4
    tmp2 = array4_init_standard([nv,nv,no,no])
    tmp4 = array4_init_standard([nv,nv,no,no])

    ! transformation time - ccsdpt_doubles and ccsd_doubles case
    call array4_contract1(tmp1,U_virt,tmp2,.true.)
    call array4_contract1(tmp3,U_virt,tmp4,.true.)
    ! ccsdpt_singles case
    call array2_matmul(tmp0,U_occ,ccsdpt_t1,'n','n',1.0E0_realk,0.0E0_realk)

    ! free tmp1 and tmp3
    call array4_free(tmp1)
    call array4_free(tmp3)
    ! free tmp0
    call array2_free(tmp0)

    ! 4. 3rd index:
    ! tmp2(b,a,I,J) --> tmp2(J,I,a,b) --> tmp1(j,I,a,b)
    ! tmp4(a,b,J,I) --> tmp4(I,J,b,a) --> tmp3(i,J,b,a)
    call array4_reorder(tmp2,[4,3,2,1])
    call array4_reorder(tmp4,[4,3,2,1])

    ! init temporary array4s, tmp1 and tmp3, once again
    tmp1 = array4_init_standard([no,no,nv,nv])
    tmp3 = array4_init_standard([no,no,nv,nv])

    ! transformation time
    call array4_contract1(tmp2,U_occ,tmp1,.true.)
    call array4_contract1(tmp4,U_occ,tmp3,.true.)

    ! free tmp2 and tmp4
    call array4_free(tmp2)
    call array4_free(tmp4)

    ! 5. 4th index:
    ! tmp1(j,I,a,b) --> tmp1(I,j,a,b) --> ccsdpt_doubles(i,j,a,b)
    ! tmp3(i,J,b,a) --> tmp3(J,i,b,a) --> ccsd_doubles(j,i,b,a)
    call array4_reorder(tmp1,[2,1,3,4])
    call array4_reorder(tmp3,[2,1,3,4])

    ! init ccsdpt_t2 and ccsd_t2 array4s once again
    ccsdpt_t2 = array4_init_standard([no,no,nv,nv])
    ccsd_t2 = array4_init_standard([no,no,nv,nv])

    ! transformation time
    call array4_contract1(tmp1,U_occ,ccsdpt_t2,.true.)
    call array4_contract1(tmp3,U_occ,ccsd_t2,.true.)

    ! free tmp1 and tmp3
    call array4_free(tmp1)
    call array4_free(tmp3)

  end subroutine ccsdpt_can_local_trans


  !> \brief: create VIRTUAL part of a triples amplitude ([a,b,c] tuple) for a fixed [i,j,k] tuple, that is, t^{***}_{ijk}
  !          saved as an array3 structure (amplitudes)
  !> \author: Janus Juul Eriksen
  !> \date: july 2012
  !> \param: oindex1, oindex2, and oindex3 are the three occupied indices of the outer loop in the ccsd(t) driver
  !> \param: no and nv are nocc and nvirt, respectively
  !> \param: doub_ampl are ccsd ampltidues, t^{ab}_{ij}
  !> \param: int_virt is a v^3 part of cbai of driver routine
  !> \param: trip holds the triples tuple [a,b,c], that is, of the size (virt)³ kept in memory
  subroutine trip_amplitudes_virt(oindex1,oindex2,oindex3,no,nv,doub_ampl_v2,int_virt_tile,trip)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, no, nv
    real(realk), dimension(nv,nv,nv), intent(in) :: int_virt_tile
    real(realk), dimension(nv,nv), intent(in) :: doub_ampl_v2
    real(realk), dimension(nv,nv,nv), intent(inout) :: trip
    !> temporary quantities
    integer :: nv2

    nv2 = nv**2

    ! important: notation adapted from eq. (14.6.60) of MEST. herein, 'e' is the running index of the
    ! equation in the book (therein: 'd')

    ! NOTE: incoming array structures are ordered according to:
    ! canTaibj(a,b,j,i) - in doub_ampl_v2 we have (a,b)
    ! int_virt_tile(c,b,a,i) - only (c,a,b)

    ! ***************************************************
    ! ** contraction time (over the virtual index 'e') **
    ! ***************************************************

    ! do v^4o^3 contraction

!$acc host_data use_device(doub_ampl_v2,int_virt_tile,trip)
    call dgemm('t','n',nv,nv2,nv,1.0E0_realk,doub_ampl_v2,nv,int_virt_tile,nv,&
                   & 0.0E0_realk,trip,nv)
!$acc end host_data

  end subroutine trip_amplitudes_virt

  !> \brief: create OCCUPIED part of a triples amplitude ([a,b,c] tuple) for a fixed [i,j,k] tuple, that is, t^{***}_{ijk}
  !          saved as an array3 structure (amplitudes)
  !> \author: Janus Juul Eriksen
  !> \date: july 2012
  !> \param: oindex1, oindex2, and oindex3 are the three occupied indices of the outer loop in the ccsd(t) driver
  !> \param: no and nv are nocc and nvirt, respectively
  !> \param: doub_ampl are ccsd ampltidues, t^{ab}_{ij}
  !> \param: int_occ is a ov part of jaik of driver routine
  !> \param: trip holds the triples tuple [c,a,b], that is, of the size (virt)³ kept in memory
  subroutine trip_amplitudes_occ(oindex1,oindex2,oindex3,no,nv,doub_ampl_ov2,int_occ_portion,trip)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, no, nv
    real(realk), dimension(no,nv), intent(in) :: int_occ_portion
    real(realk), dimension(no,nv,nv), intent(in) :: doub_ampl_ov2
    real(realk), dimension(nv,nv,nv), intent(inout) :: trip
    !> temporary quantities
    integer :: nv2

    nv2 = nv**2

    ! important: notation adapted from eq. (14.6.60) of MEST. herein, 'm' is the running index of the
    ! equation in the book (therein: 'l')

    ! NOTE: incoming array structures are ordered according to:
    ! canTaibj(a,b,j,i) - in doub_ampl_ov2 we have (j,a,b)
    ! canAIJK(j,a,i,k) - in int_occ_portion we have only (j,a)

    ! ****************************************************
    ! ** contraction time (over the occupied index 'm') **
    ! ****************************************************

    ! do v^3o^4 contraction

!$acc host_data use_device(int_occ_portion,doub_ampl_ov2,trip)
    call dgemm('t','n',nv,nv2,no,-1.0E0_realk,int_occ_portion,no,doub_ampl_ov2,no,&
                   & 1.0E0_realk,trip,nv)
!$acc end host_data

  end subroutine trip_amplitudes_occ


  !> \author: Janus Juul Eriksen
  !> \date: july 2012
  !> \param: oindex1, oindex2, and oindex3 are the three occupied indices of the outer loop in the ccsd(t) driver
  !> \param: no and nv are nocc and nvirt, respectively
  !> \param: eigenocc and eigenvirt are vectors containing occupied and virtual orbital energies, respectively
  !> \param: amplitudes are the final triples amplitude tuple [a,b,c], that is, of the size (virt)³ kept in memory
  subroutine trip_denom(oindex1,oindex2,oindex3,no,nv,eigenocc,eigenvirt,trip)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, no, nv
    real(realk), intent(inout) :: eigenocc(no), eigenvirt(nv)
    real(realk), dimension(nv,nv,nv), intent(inout) :: trip
    !> temporary quantities
    integer :: trip_type, a, b, c
    real(realk) :: e_orb, e_orb_occ

    ! at first, calculate the sum of the three participating occupied orbital energies, as this
    ! is a constant for the three incomming occupied indices

    e_orb_occ = eigenocc(oindex1) + eigenocc(oindex2) + eigenocc(oindex3)

#ifdef VAR_OPENACC
!$acc parallel present(trip,eigenvirt,eigenocc)
!$acc loop gang
#else
!$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(a,b,c,e_orb),SHARED(nv,trip,eigenvirt,e_orb_occ)
#endif
 arun_0: do a=1,nv
#ifdef VAR_OPENACC
!$acc loop worker
#endif
    brun_0: do b=1,nv
#ifdef VAR_OPENACC
!$acc loop vector
#endif
       crun_0: do c=1,nv

                  trip(c,b,a) = trip(c,b,a) / (e_orb_occ - eigenvirt(a) - eigenvirt(b) - eigenvirt(c))

               end do crun_0
#ifdef VAR_OPENACC
!$acc end loop
#endif
            end do brun_0
#ifdef VAR_OPENACC
!$acc end loop
#endif
         end do arun_0
#ifdef VAR_OPENACC
!$acc end loop
!$acc end parallel
#else
!$OMP END PARALLEL DO
#endif

  end subroutine trip_denom


  !> brief: do the first of the two contraction over 'cdkl' (here: 'c' and 'd', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.63) of MEST
  !> author: Janus Juul Eriksen
  !> date: august 2012
  !> param: oindex1-oindex3 are outside loop indices of driver routine. int_normal is abij of driver.
  !> nv is nvirt and T_star is ccsdpt_singles of driver. trip_ampl is the triples amplitude array.
  subroutine ccsdpt_contract_11(oindex1,oindex2,oindex3,nv,no,int_normal_23,int_normal_32,T_star,trip_ampl,special)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, nv, no
    real(realk), dimension(nv,no) :: T_star
    real(realk), dimension(nv,nv) :: int_normal_23, int_normal_32
    real(realk), dimension(nv,nv,nv) :: trip_ampl
    logical, intent(in) :: special
    !> temporary quantities
    integer :: contraction_type, nv2

    nv2 = nv**2

    ! NOTE: incoming array4 structures are ordered according to:
    ! canAIBJ(c,d,k,l) (MEST nomenclature)
    ! T_ast_0(a,i)

    ! determine which type of contraction is to be performed
    contraction_type = -1
    ! is this a special contraction, i.e., can we handle 211 and 212 contractions in one go?
    if (special) contraction_type = 0
    ! otherwise, do the default 211 contraction
    if (.not. special) contraction_type = 1

    ! contraction time (here: over virtual indices 'c' and 'd') with "coulumb minus exchange"
    ! version of canAIBJ (2 * canAIJK(c,k,d,l) - canAIBC(c,l,d,k))

    TypeofContraction_11: select case(contraction_type)

    case(0)

       ! here, the coulumb and exchange parts will be equal and we thus only need to contract with the coulumb part. 

       ! now contract coulumb term over both indices
       call dgemm('n','n',nv,1,nv2,&
                & 1.0E0_realk,trip_ampl,nv,int_normal_23,&
                & nv2,1.0E0_realk,T_star(1,oindex1),nv)

    case(1)

       ! now contract coulumb term over both indices
       call dgemm('n','n',nv,1,nv2,&
                & 2.0E0_realk,trip_ampl,nv,int_normal_23,&
                & nv2,1.0E0_realk,T_star(1,oindex1),nv)

       ! now contract exchange term over both indices
       call dgemm('n','n',nv,1,nv2,&
                & -1.0E0_realk,trip_ampl,nv,int_normal_32,&
                & nv2,1.0E0_realk,T_star(1,oindex1),nv)

    end select TypeofContraction_11

  end subroutine ccsdpt_contract_11

  !> brief: do the second of the two contraction over 'cdkl' (here: 'c' and 'd', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.63) of MEST
  !> author: Janus Juul Eriksen
  !> date: august 2012
  !> param: oindex1-oindex3 are outside loop indices of driver routine. int_normal is abij of driver.
  !> nv is nvirt and T_star is ccsdpt_singles of driver. trip_ampl is the triples amplitude array.
  subroutine ccsdpt_contract_12(oindex1,oindex2,oindex3,nv,no,int_normal_21,int_normal_12,T_star,trip_ampl,special)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, nv, no
    real(realk), dimension(nv,no) :: T_star
    real(realk), dimension(nv,nv) :: int_normal_21, int_normal_12
    real(realk), dimension(nv,nv,nv) :: trip_ampl
    logical, intent(in) :: special
    !> temporary quantities
    integer :: contraction_type, nv2

    nv2 = nv**2

    ! NOTE: incoming array4 structures are ordered according to:
    ! canAIBJ(c,d,k,l) (MEST nomenclature)
    ! T_ast_0(a,i)

    ! contraction time (here: over virtual indices 'c' and 'd') with "coulumb minus exchange"
    ! version of canAIBJ (2 * canAIJK(c,k,d,l) - canAIBC(c,l,d,k))

    ! determine which type of contraction is to be performed
    contraction_type = -1
    ! is this a special contraction, i.e., can we handle 211 and 212 contractions in one go?
    if (special) contraction_type = 0
    ! otherwise, do the default 211 contraction
    if (.not. special) contraction_type = 1

    ! contraction time (here: over virtual indices 'c' and 'd') with "coulumb minus exchange"
    ! version of canAIBJ (2 * canAIJK(c,k,d,l) - canAIBC(c,l,d,k))

    TypeofContraction_12: select case(contraction_type)

    case(0)

       ! here, the coulumb and exchange parts will be equal
       ! and we thus only need to contract with (-1)*coulumb part. 

       ! now contract coulumb term over both indices
       call dgemm('n','n',nv,1,nv2,&
                & -1.0E0_realk,trip_ampl,nv,int_normal_21,&
                & nv2,1.0E0_realk,T_star(1,oindex3),nv)

    case(1)

       ! now contract coulumb term over both indices
       call dgemm('n','n',nv,1,nv2,&
                & -2.0E0_realk,trip_ampl,nv,int_normal_21,&
                & nv2,1.0E0_realk,T_star(1,oindex3),nv)

       ! now contract exchange term over both indices
       call dgemm('n','n',nv,1,nv2,&
                & 1.0E0_realk,trip_ampl,nv,int_normal_12,&
                & nv2,1.0E0_realk,T_star(1,oindex3),nv)

    end select TypeofContraction_12

  end subroutine ccsdpt_contract_12


  !> brief: do the first of the two contractions over 'cdk' (here: 'cd', 'k' is the summation in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: august 2012
  !> param: oindex1-oindex3 are outside loop indices of driver routine.
  !> nv is nvirt and T_star is T_ast_1 of driver. trip_ampl is the triples amplitude array.
  !> int_virt_tile is a v^3 tile determined by driver occ index
  !> tmp_g is a 3d work array
  subroutine ccsdpt_contract_211(oindex1,oindex2,oindex3,nv,no,&
       & T_star,tmp_g,trip_ampl,int_virt_tile,special)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, nv, no
    real(realk), dimension(nv,nv,no,no) :: T_star
    real(realk), dimension(nv,nv,nv) :: tmp_g,trip_ampl
    real(realk), dimension(nv,nv,nv), intent(in) :: int_virt_tile
    logical, intent(in) :: special
    !> temporary quantities
    integer :: contraction_type, nv2

    nv2 = nv**2

    ! NOTE: incoming array(4) structures are ordered according to:
    ! int_virt_tile(c,b,a,x)
    ! T_ast_1(a,b,i,j)

    ! determine which type of contraction is to be performed
    contraction_type = -1
    ! is this a special contraction, i.e., can we handle 211 and 212 contractions in one go?
    if (special) contraction_type = 0
    ! otherwise, do the default 211 contraction
    if (.not. special) contraction_type = 1

    ! contraction time (here: over virtual indices 'c' and 'd') with "coulumb minus exchange"
    ! version of canAIBC (2 * canAIBC(b,c,k,d) - canAIBC(b,d,k,c)) and (if special) canAIBC(b,c,k,d)

    TypeofContraction_211: select case(contraction_type)

    case(0)

       ! note: here we collect contract over L_{dkbc} and g_{dkbc} in one go.

       ! reorder to obtain coulumb term, tmp_g(c,d,b)
       call array_reorder_3d(1.0E0_realk,int_virt_tile,nv,nv,nv,[1,3,2],0.0E0_realk,tmp_g)

       ! now contract coulumb term over 2 first indices
       call dgemm('t','n',nv,nv,nv2, &
            1.0E0_realk,trip_ampl,nv2,tmp_g,nv2,1.0E0_realk,T_star(1,1,oindex1,oindex2),nv)

       ! reorder to obtain exchange term, tmp_g(d,c,b)
       call array_reorder_3d(1.0E0_realk,int_virt_tile,nv,nv,nv,[3,1,2],0.0E0_realk,tmp_g)

       ! now contract exchange term over 2 first indices2
       call dgemm('t','n',nv,nv,nv2, &
            -1.0E0_realk,trip_ampl,nv2,tmp_g,nv2,1.0E0_realk,T_star(1,1,oindex1,oindex2),nv)

    case(1)

       ! note: here we contract over L_{dkbc}.

       ! reorder to obtain coulumb term, tmp_g(c,d,b)
       call array_reorder_3d(1.0E0_realk,int_virt_tile,nv,nv,nv,[1,3,2],0.0E0_realk,tmp_g)

       ! now contract coulumb term over 2 first indices
       call dgemm('t','n',nv,nv,nv2, &
            2.0E0_realk,trip_ampl,nv2,tmp_g,nv2,1.0E0_realk,T_star(1,1,oindex1,oindex2),nv)

       ! reorder to obtain exchange term, tmp_g(d,c,b)
       call array_reorder_3d(1.0E0_realk,int_virt_tile,nv,nv,nv,[3,1,2],0.0E0_realk,tmp_g)

       ! now contract exchange term over 2 first indices
       call dgemm('t','n',nv,nv,nv2, &
            -1.0E0_realk,trip_ampl,nv2,tmp_g,nv2,1.0E0_realk,T_star(1,1,oindex1,oindex2),nv)

    end select TypeofContraction_211

  end subroutine ccsdpt_contract_211

  !> brief: do the second of the two contractions over 'cdk' (here: 'cd', 'k' is the summation in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: august 2012
  !> param: oindex1-oindex3 are outside loop indices of driver routine.
  !> nv is nvirt and T_star is T_ast_1 of driver. trip_ampl is the triples amplitude array.
  !> int_virt_tile is a v^3 tile determined by driver occ index
  !> tmp_g is a 3d work array
  subroutine ccsdpt_contract_212(oindex1,oindex2,oindex3,nv,no,&
       & T_star,tmp_g,trip_ampl,int_virt_tile)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, nv, no
    real(realk), dimension(nv,nv,no,no) :: T_star
    real(realk), dimension(nv,nv,nv), intent(in) :: int_virt_tile
    real(realk), dimension(nv,nv,nv) :: tmp_g,trip_ampl
    !> temporary quantities
    integer :: nv2

    nv2 = nv**2

    ! NOTE: incoming array(4) structures are ordered according to:
    ! int_virt_tile(c,b,a,x)
    ! T_ast_1(a,b,i,j)

    ! contraction time (here: over virtual indices 'c' and 'd') with canAIBC(b,c,k,d)

    ! reorder to obtain tmp_g(c,d,b)
    call array_reorder_3d(1.0E0_realk,int_virt_tile,nv,nv,nv,[1,3,2],0.0E0_realk,tmp_g)

    ! now contract coulumb term over 2 first indices
    call dgemm('t','n',nv,nv,nv2,&
         & -1.0E0_realk,trip_ampl,nv2,tmp_g,nv2,1.0E0_realk,T_star(1,1,oindex3,oindex2),nv)

  end subroutine ccsdpt_contract_212


  !> brief: do the first of the two contractions over 'ckl' (here: 'c', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: august 2012
  !> param: oindex1-oindex3 are outside loop indices of driver routine.
  !> nv is nvirt and T_star is T_ast_2 of driver. trip_ampl is the triples amplitud array.
  subroutine ccsdpt_contract_221(oindex1,oindex2,oindex3,no,nv,&
                               & int_occ_23,int_occ_32,T_star,trip_ampl,special)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, no, nv
    real(realk), dimension(no,nv,nv,no) :: T_star
    real(realk), dimension(no,nv) :: int_occ_23, int_occ_32
    real(realk), dimension(nv,nv,nv) :: trip_ampl
    logical, intent(in) :: special
    !> temporary quantities
    integer :: contraction_type, nv2

    nv2 = nv**2

    ! NOTE: incoming array4 structures are ordered according to:
    ! canAIJK(j,c,l,k) (MEST nomenclature)
    ! T_ast_2(j,a,b,i)

    ! determine which type of contraction is to be performed
    contraction_type = -1
    ! is this a special contraction, i.e., can we handle 221 and 222 contractions in one go?
    if (special) contraction_type = 0
    ! otherwise, do the default 221 contraction
    if (.not. special) contraction_type = 1

    ! contraction time (here: over virtual index 'c') with "coulumb minus exchange"
    ! version of canAIBC (2 * canAIJK(k,j,l,c) - canAIBC(l,j,k,c)) and (if special) canAIBC(k,j,l,c)

    TypeofContraction_221: select case(contraction_type)

    case(0)

       ! now contract coulumb term over first index.
       ! for this special case, we only have to subtract one coulumb term
       call dgemm('n','n',no,nv2,nv,-1.0E0_realk,int_occ_32,no,&
                      & trip_ampl,nv,1.0E0_realk,T_star(1,1,1,oindex1),no)

       ! now contract exchange term over first index
       call dgemm('n','n',no,nv2,nv,1.0E0_realk,int_occ_23,no,&
                      & trip_ampl,nv,1.0E0_realk,T_star(1,1,1,oindex1),no)

    case(1)
 
       ! now contract coulumb term over first index
       call dgemm('n','n',no,nv2,nv,-2.0E0_realk,int_occ_32,no,&
                      & trip_ampl,nv,1.0E0_realk,T_star(1,1,1,oindex1),no)

       ! now contract exchange term over first index
       call dgemm('n','n',no,nv2,nv,1.0E0_realk,int_occ_23,no,&
                      & trip_ampl,nv,1.0E0_realk,T_star(1,1,1,oindex1),no)

    end select TypeofContraction_221

  end subroutine ccsdpt_contract_221


  !> brief: do the second of the two contractions over 'ckl' (here: 'c', 'k' and 'l' are summations in driver routine)
  !         in eq. (14.6.64) of MEST
  !> author: Janus Juul Eriksen
  !> date: august 2012
  !> param: oindex1-oindex3 are outside loop indices of driver routine.
  !> nv is nvirt and T_star is T_ast_2 of driver. trip_ampl is the triples amplitud array.
  subroutine ccsdpt_contract_222(oindex1,oindex2,oindex3,no,nv,int_occ_12,T_star,trip_ampl)

    implicit none
    !> input
    integer, intent(in) :: oindex1, oindex2, oindex3, no, nv
    real(realk), dimension(no,nv,nv,no) :: T_star
    real(realk), dimension(no,nv) :: int_occ_12
    real(realk), dimension(nv,nv,nv) :: trip_ampl
    !> integer
    integer :: nv2

    nv2 = nv**2

    ! NOTE: incoming array4 structures are ordered according to:
    ! canAIJK(j,c,l,k) (MEST nomenclature)
    ! T_ast_2(j,a,b,i)

    ! contraction time (here: over virtual index 'c') with canAIBC(k,j,l,c)

    ! contract coulumb term over first index
    call dgemm('n','n',no,nv2,nv,1.0E0_realk,int_occ_12,no,&
                   & trip_ampl,nv,1.0E0_realk,T_star(1,1,1,oindex3),no)

  end subroutine ccsdpt_contract_222


  !> \brief: calculate E[5] contribution to single fragment ccsd(t) energy correction
  !> \author: Janus Eriksen
  !> \date: september 2012
  subroutine ccsdpt_energy_e5_frag(MyFragment,ccsd_singles,ccsdpt_singles)

    implicit none

    !> fragment info
    type(decfrag), intent(inout) :: MyFragment
    ! ccsd and ccsd(t) singles amplitudes
    type(array2), intent(inout) :: ccsd_singles, ccsdpt_singles
    !> integers
    integer :: nocc_eos, nvirt_eos, i,a, i_eos, a_eos
    !> temp energy
    real(realk) :: energy_tmp, ccsdpt_e5

    ! init dimensions
    nocc_eos = MyFragment%noccEOS
    nvirt_eos = MyFragment%nunoccEOS

    ! ***********************
    !   do E[5] energy part
    ! ***********************

    ! init energy reals to be on the safe side.
    ! note: OccEnergyPT and VirtEnergyPT have been initialized in the e4 routine.
    MyFragment%energies(FRAGMODEL_OCCpT5) = 0.0E0_realk
    MyFragment%energies(FRAGMODEL_VIRTpT5) = 0.0E0_realk

    ! init temp energy
    ccsdpt_e5 = 0.0E0_realk

                    !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,i_eos,a,a_eos,energy_tmp),&
                    !$OMP SHARED(ccsd_singles,ccsdpt_singles,nocc_eos,nvirt_eos,MyFragment),&
                    !$OMP REDUCTION(+:ccsdpt_e5)
  ido_frag_singles: do i=1,nocc_eos
                    i_eos = MyFragment%idxo(i)
     ado_frag_singles: do a=1,nvirt_eos
                       a_eos = MyFragment%idxu(a)

                          energy_tmp = ccsd_singles%val(a_eos,i_eos) * ccsdpt_singles%val(a_eos,i_eos)
                          ccsdpt_e5 = ccsdpt_e5 + energy_tmp

                       end do ado_frag_singles
                    end do ido_frag_singles
                    !$OMP END PARALLEL DO

    MyFragment%energies(FRAGMODEL_OCCpT5) = 2.0E0_realk * ccsdpt_e5

    ! insert into occ. part. scheme part
    MyFragment%energies(FRAGMODEL_OCCpT) = MyFragment%energies(FRAGMODEL_OCCpT) + MyFragment%energies(FRAGMODEL_OCCpT5)

    ! *********************************
    ! do unoccupied partitioning scheme
    ! *********************************

    ! singles contribution is the same as in occupied partitioning scheme
    MyFragment%energies(FRAGMODEL_VIRTpT) = MyFragment%energies(FRAGMODEL_VIRTpT) + MyFragment%energies(FRAGMODEL_OCCpT5)
    ! insert into virt_e5 part
    MyFragment%energies(FRAGMODEL_VIRTpT5) = MyFragment%energies(FRAGMODEL_VIRTpT5) + MyFragment%energies(FRAGMODEL_OCCpT5)

    ! ******************************
    !   done with E[5] energy part
    ! ******************************

  end subroutine ccsdpt_energy_e5_frag 


  !> \brief: calculate E[5] contribution to pair fragment ccsd(t) energy correction
  !> \author: Janus Eriksen
  !> \date: september 2012
  subroutine ccsdpt_energy_e5_pair(Fragment1,Fragment2,PairFragment,ccsd_singles,ccsdpt_singles)

    implicit none

    !> fragment # 1 in the pair fragment
    type(decfrag),intent(inout) :: Fragment1
    !> fragment # 2 in the pair fragment
    type(decfrag),intent(inout) :: Fragment2
    !> fragment info
    type(decfrag), intent(inout) :: PairFragment
    ! ccsd and ccsd(t) singles amplitudes
    type(array2), intent(inout) :: ccsd_singles, ccsdpt_singles
    !> integers
    integer :: nocc_eos, nvirt_eos, i, a, idx, adx, AtomI, AtomA, i_eos, a_eos
    !> logicals to avoid double counting
    logical :: occ_in_frag_1, virt_in_frag_1, occ_in_frag_2, virt_in_frag_2

    ! init dimensions
    nocc_eos = PairFragment%noccEOS
    nvirt_eos = PairFragment%nunoccEOS

    ! ***********************
    !   do E[5] energy part
    ! ***********************

    ! init energy reals to be on the safe side.
    ! note: OccEnergyPT and VirtEnergyPT have been initialized in the e4 routine.
    PairFragment%energies(FRAGMODEL_OCCpT5) = 0.0E0_realk
    PairFragment%energies(FRAGMODEL_VIRTpT5) = 0.0E0_realk

  ido_pair_singles: do i=1,nocc_eos
                    i_eos = PairFragment%idxo(i)
                    AtomI = PairFragment%occAOSorb(i_eos)%CentralAtom
     ado_pair_singles: do a=1,nvirt_eos
                       a_eos = PairFragment%idxu(a)
                       AtomA = PairFragment%unoccAOSorb(a_eos)%CentralAtom

                       ! occ in frag # 1, virt in frag # 2

                          occ_in_frag_1 = .false.
                          do idx = 1,Fragment1%nEOSatoms
                             if (Fragment1%EOSatoms(idx) .eq. AtomI) then
                                occ_in_frag_1 = .true.
                             end if
                          end do

                          virt_in_frag_2 = .false.
                          do adx = 1,Fragment2%nEOSatoms
                             if (Fragment2%EOSatoms(adx) .eq. AtomA) then
                                virt_in_frag_2 = .true.
                             end if
                          end do

                          if (occ_in_frag_1 .and. virt_in_frag_2) then
                             PairFragment%energies(FRAGMODEL_OCCpT5) = PairFragment%energies(FRAGMODEL_OCCpT5) &
                               & + 2.0E0_realk * ccsd_singles%val(a_eos,i_eos) &
                               & * ccsdpt_singles%val(a_eos,i_eos)
                          end if

                       ! virt in frag # 1, occ in frag # 2
                            
                          occ_in_frag_2 = .false. 
                          do idx = 1,Fragment2%nEOSatoms
                             if (Fragment2%EOSatoms(idx) .eq. AtomI) then
                                occ_in_frag_2 = .true.
                             end if
                          end do 
                             
                          virt_in_frag_1 = .false.
                          do adx = 1,Fragment1%nEOSatoms
                             if (Fragment1%EOSatoms(adx) .eq. AtomA) then
                                virt_in_frag_1 = .true.
                             end if
                          end do
                             
                          if (occ_in_frag_2 .and. virt_in_frag_1) then
                             PairFragment%energies(FRAGMODEL_OCCpT5) = PairFragment%energies(FRAGMODEL_OCCpT5) &
                               & + 2.0E0_realk * ccsd_singles%val(a_eos,i_eos) &
                               & * ccsdpt_singles%val(a_eos,i_eos)
                          end if

                       ! sanity checks

                          if (.not. (occ_in_frag_1 .or. occ_in_frag_2)) then
                             call lsquit('Problem in evaluation of E[5] contr. &
                                  & to pair fragment energy: occ orbital neither in frag 1 or 2',DECinfo%output)        
                          end if
                          if (.not. (virt_in_frag_1 .or. virt_in_frag_2)) then
                             call lsquit('Problem in evaluation of E[5] contr. &
                                  & to pair fragment energy: virt orbital neither in frag 1 or 2',DECinfo%output)
                          end if
                          if (occ_in_frag_1 .and. occ_in_frag_2) then
                             call lsquit('Problem in evaluation of E[5] contr. &
                                  & to pair fragment energy: occ orbital both in frag 1 or 2',DECinfo%output)
                          end if
                          if (virt_in_frag_1 .and. virt_in_frag_2) then
                             call lsquit('Problem in evaluation of E[5] contr. &
                                  & to pair fragment energy: virt orbital both in frag 1 or 2',DECinfo%output)
                          end if

                       end do ado_pair_singles
                    end do ido_pair_singles

    ! insert into occ. part. scheme part
    PairFragment%energies(FRAGMODEL_OCCpT) = PairFragment%energies(FRAGMODEL_OCCpT) + PairFragment%energies(FRAGMODEL_OCCpT5)

    ! *********************************
    ! do unoccupied partitioning scheme
    ! *********************************

    ! singles contribution is the same as in occupied partitioning scheme
    PairFragment%energies(FRAGMODEL_VIRTpT) = PairFragment%energies(FRAGMODEL_VIRTpT) + PairFragment%energies(FRAGMODEL_OCCpT5)
    ! insert into virt_e5 part
    PairFragment%energies(FRAGMODEL_VIRTpT5) = PairFragment%energies(FRAGMODEL_VIRTpT5) + PairFragment%energies(FRAGMODEL_OCCpT5)

    ! ******************************
    !   done with E[5] energy part
    ! ******************************

  end subroutine ccsdpt_energy_e5_pair


  !> \brief: calculate E[4] contribution to single fragment ccsd(t) energy correction
  !> \author: Janus Eriksen
  !> \date: september 2012
  subroutine ccsdpt_energy_e4_frag(MyFragment,ccsd_doubles,ccsdpt_doubles,&
                             & occ_contribs,virt_contribs,fragopt_pT)

    implicit none

    !> fragment info
    type(decfrag), intent(inout) :: MyFragment
    ! ccsd and ccsd(t) doubles amplitudes
    type(array4), intent(inout) :: ccsd_doubles, ccsdpt_doubles
    !> is this called from inside the ccsd(t) fragment optimization routine?
    logical, optional, intent(in) :: fragopt_pT
    !> incomming orbital contribution vectors
    real(realk), intent(inout) :: occ_contribs(MyFragment%noccAOS), virt_contribs(MyFragment%nunoccAOS)
    !> integers
    integer :: nocc_eos, nocc_aos, nvirt_eos, nvirt_aos, i,j,a,b, i_eos, j_eos, a_eos, b_eos
    !> energy reals
    real(realk) :: energy_tmp, energy_res_cou, energy_res_exc

    ! init dimensions
    nocc_eos = MyFragment%noccEOS
    nvirt_eos = MyFragment%nunoccEOS
    nocc_aos = MyFragment%noccAOS
    nvirt_aos = MyFragment%nunoccAOS

    ! **************************************************************
    ! ************** do energy for single fragment *****************
    ! **************************************************************

    ! ***********************
    !   do E[4] energy part
    ! ***********************

    ! init energy reals to be on the safe side
    ! note: OccEnergyPT and VirtEnergyPT is also initialized from in here
    !       as this (e4) routine is called before the e5 routine
    MyFragment%energies(FRAGMODEL_OCCpT) = 0.0E0_realk
    MyFragment%energies(FRAGMODEL_VIRTpT) = 0.0E0_realk
    MyFragment%energies(FRAGMODEL_OCCpT4) = 0.0E0_realk
    MyFragment%energies(FRAGMODEL_VIRTpT4) = 0.0E0_realk

    ! *******************************
    ! do occupied partitioning scheme
    ! *******************************

    energy_res_cou = 0.0E0_realk
    energy_res_exc = 0.0E0_realk

                        !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,i_eos,j,j_eos,a,b,energy_tmp),&
                        !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nocc_eos,nvirt_aos,MyFragment),&
                        !$OMP REDUCTION(+:energy_res_cou),REDUCTION(+:virt_contribs)
  jdo_frag_doubles_cou: do j=1,nocc_eos
                        j_eos = MyFragment%idxo(j)
     ido_frag_doubles_cou: do i=1,nocc_eos
                           i_eos = MyFragment%idxo(i)

                              do b=1,nvirt_aos
                                 do a=1,nvirt_aos

                                    energy_tmp = 4.0E0_realk * ccsd_doubles%val(a,b,i_eos,j_eos) &
                                                   & * ccsdpt_doubles%val(a,b,i_eos,j_eos)
                                    energy_res_cou = energy_res_cou + energy_tmp

                                    ! update contribution from aos orbital a
                                    virt_contribs(a) = virt_contribs(a) + energy_tmp

                                    ! update contribution from aos orbital b 
                                    ! (only if different from aos orbital a to avoid double counting)
                                    if (a .ne. b) virt_contribs(b) = virt_contribs(b) + energy_tmp

                                 end do
                              end do

                           end do ido_frag_doubles_cou
                        end do jdo_frag_doubles_cou
                        !$OMP END PARALLEL DO

    ! reorder from (a,b,i,j) to (a,b,j,i)
    call array4_reorder(ccsd_doubles,[1,2,4,3])

                        !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,i_eos,j,j_eos,a,b,energy_tmp),&
                        !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nocc_eos,nvirt_aos,MyFragment),&
                        !$OMP REDUCTION(+:energy_res_exc),REDUCTION(+:virt_contribs)
  jdo_frag_doubles_exc: do j=1,nocc_eos
                        j_eos = MyFragment%idxo(j)
     ido_frag_doubles_exc: do i=1,nocc_eos
                           i_eos = MyFragment%idxo(i)

                              do b=1,nvirt_aos
                                 do a=1,nvirt_aos

                                    energy_tmp = 2.0E0_realk * ccsd_doubles%val(a,b,i_eos,j_eos) &
                                                   & * ccsdpt_doubles%val(a,b,i_eos,j_eos)
                                    energy_res_exc = energy_res_exc - energy_tmp

                                    ! update contribution from aos orbital a
                                    virt_contribs(a) = virt_contribs(a) - energy_tmp

                                    ! update contribution from aos orbital b 
                                    ! (only if different from aos orbital a to avoid double counting)
                                    if (a .ne. b) virt_contribs(b) = virt_contribs(b) - energy_tmp

                                 end do
                              end do

                           end do ido_frag_doubles_exc
                        end do jdo_frag_doubles_exc
                        !$OMP END PARALLEL DO

    !get total fourth--order energy contribution
    MyFragment%energies(FRAGMODEL_OCCpT4) = energy_res_cou + energy_res_exc

    ! insert into occ. part. scheme part
    MyFragment%energies(FRAGMODEL_OCCpT) = MyFragment%energies(FRAGMODEL_OCCpT) + MyFragment%energies(FRAGMODEL_OCCpT4)

    ! *********************************
    ! do unoccupied partitioning scheme
    ! *********************************

    ! initially, reorder ccsd_doubles and ccsdpt_doubles
    ! ccsd_doubles from from (a,b,j,i) sequence to (j,i,a,b) sequence
    ! ccsdpt_doubles from from (a,b,i,j) sequence to (i,j,a,b) sequence
    call array4_reorder(ccsd_doubles,[3,4,1,2])
    call array4_reorder(ccsdpt_doubles,[3,4,1,2])

    energy_res_cou = 0.0E0_realk
    energy_res_exc = 0.0E0_realk

                        !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(a,a_eos,b,b_eos,i,j,energy_tmp),&
                        !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nvirt_eos,nocc_aos,MyFragment),&
                        !$OMP REDUCTION(+:energy_res_exc),REDUCTION(+:occ_contribs)
  bdo_frag_doubles_exc: do b=1,nvirt_eos
                        b_eos = MyFragment%idxu(b)
     ado_frag_doubles_exc: do a=1,nvirt_eos
                           a_eos = MyFragment%idxu(a)

                              do j=1,nocc_aos
                                 do i=1,nocc_aos

                                    energy_tmp = 2.0E0_realk * ccsd_doubles%val(i,j,a_eos,b_eos) &
                                                   & * ccsdpt_doubles%val(i,j,a_eos,b_eos)
                                    energy_res_exc = energy_res_exc - energy_tmp

                                    ! update contribution from aos orbital i
                                    occ_contribs(i) = occ_contribs(i) - energy_tmp

                                    ! update contribution from aos orbital j 
                                    ! (only if different from aos orbital i to avoid double counting)
                                    if (i .ne. j) occ_contribs(j) = occ_contribs(j) - energy_tmp

                                 end do
                              end do

                           end do ado_frag_doubles_exc
                        end do bdo_frag_doubles_exc
                        !$OMP END PARALLEL DO

    ! reorder form (j,i,a,b) to (i,j,a,b)
    call array4_reorder(ccsd_doubles,[2,1,3,4])

                        !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(a,a_eos,b,b_eos,i,j,energy_tmp),&
                        !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nvirt_eos,nocc_aos,MyFragment),&
                        !$OMP REDUCTION(+:energy_res_cou),REDUCTION(+:occ_contribs)
  bdo_frag_doubles_cou: do b=1,nvirt_eos
                        b_eos = MyFragment%idxu(b)
     ado_frag_doubles_cou: do a=1,nvirt_eos
                           a_eos = MyFragment%idxu(a)

                              do j=1,nocc_aos
                                 do i=1,nocc_aos

                                    energy_tmp = 4.0E0_realk * ccsd_doubles%val(i,j,a_eos,b_eos) &
                                                   & * ccsdpt_doubles%val(i,j,a_eos,b_eos)
                                    energy_res_cou = energy_res_cou + energy_tmp

                                    ! update contribution from aos orbital i
                                    occ_contribs(i) = occ_contribs(i) + energy_tmp

                                    ! update contribution from aos orbital j 
                                    ! (only if different from aos orbital i to avoid double counting)
                                    if (i .ne. j) occ_contribs(j) = occ_contribs(j) + energy_tmp

                                 end do
                              end do

                           end do ado_frag_doubles_cou
                        end do bdo_frag_doubles_cou
                        !$OMP END PARALLEL DO

    !get total fourth--order energy contribution
    MyFragment%energies(FRAGMODEL_VIRTpT4) = energy_res_cou + energy_res_exc

    ! insert into virt. part. scheme part
    MyFragment%energies(FRAGMODEL_VIRTpT) = MyFragment%energies(FRAGMODEL_VIRTpT) + MyFragment%energies(FRAGMODEL_VIRTpT4)

    ! ******************************
    !   done with E[4] energy part
    ! ******************************

    ! ************************************************************************
    !   as we need to reuse the ccsd doubles in the fragment optimization,
    !   we here reorder back into (a,i,b,j) sequence IF fragopt_pT == .true. 
    ! ************************************************************************

    if (present(fragopt_pT)) then

       ! reorder from (i,j,a,b) to (a,i,b,j)
       if (fragopt_pT) call array4_reorder(ccsd_doubles,[3,1,4,2])

    end if

    ! *******************************************************************
    ! ************** done w/ energy for single fragment *****************
    ! *******************************************************************

  end subroutine ccsdpt_energy_e4_frag


  !> \brief: calculate E[4] contribution to pair fragment ccsd(t) energy correction
  !> \author: Janus Eriksen
  !> \date: september 2012
  subroutine ccsdpt_energy_e4_pair(Fragment1,Fragment2,PairFragment,ccsd_doubles,ccsdpt_doubles)

    implicit none

    !> fragment # 1 in the pair fragment
    type(decfrag),intent(inout) :: Fragment1
    !> fragment # 2 in the pair fragment
    type(decfrag),intent(inout) :: Fragment2
    !> pair fragment info
    type(decfrag), intent(inout) :: PairFragment
    ! ccsd and ccsd(t) doubles amplitudes
    type(array4), intent(inout) :: ccsd_doubles, ccsdpt_doubles
    ! logical pointers for keeping hold of which pairs are to be handled
    logical, pointer :: dopair_occ(:,:), dopair_virt(:,:)
    !> integers
    integer :: nocc_eos, nocc_aos, nvirt_eos, nvirt_aos, i,j,a,b, i_eos, j_eos, a_eos, b_eos
    !> temporary energy arrays
    type(array2) :: energy_interm_cou, energy_interm_exc, energy_interm_ccsdpt
    !> energy reals
    real(realk) :: energy_tmp, energy_res_cou, energy_res_exc  

    ! init dimensions
    nocc_eos = PairFragment%noccEOS
    nvirt_eos = PairFragment%nunoccEOS
    nocc_aos = PairFragment%noccAOS
    nvirt_aos = PairFragment%nunoccAOS

    ! which pairs are to be included for occ and unocc space (avoid double counting)
    call mem_alloc(dopair_occ,nocc_eos,nocc_eos)
    call mem_alloc(dopair_virt,nvirt_eos,nvirt_eos)
    call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)
    call which_pairs_unocc(Fragment1,Fragment2,PairFragment,dopair_virt)

    ! *************************************************************
    ! ************** do energy for pair fragments *****************
    ! *************************************************************

    ! ***********************
    !   do E[4] energy part
    ! ***********************

    ! init energy reals to be on the safe side
    ! note: OccEnergyPT and VirtEnergyPT is also initialized from in here
    !       as this (e4) routine is called before the e5 routine
    PairFragment%energies(FRAGMODEL_OCCpT) = 0.0E0_realk
    PairFragment%energies(FRAGMODEL_VIRTpT) = 0.0E0_realk
    PairFragment%energies(FRAGMODEL_OCCpT4) = 0.0E0_realk
    PairFragment%energies(FRAGMODEL_VIRTpT4) = 0.0E0_realk

    ! *******************************
    ! do occupied partitioning scheme
    ! *******************************

    energy_res_cou = 0.0E0_realk
    energy_res_exc = 0.0E0_realk

                        !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,i_eos,j,j_eos,a,b,energy_tmp),&
                        !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nocc_eos,nvirt_aos,&
                        !$OMP PairFragment,dopair_occ),REDUCTION(+:energy_res_cou)
  jdo_pair_doubles_cou: do j=1,nocc_eos
                        j_eos = PairFragment%idxo(j)
     ido_pair_doubles_cou: do i=1,nocc_eos
                           i_eos = PairFragment%idxo(i)

                              if (.not. dopair_occ(i,j)) cycle ido_pair_doubles_cou

                              do b=1,nvirt_aos
                                 do a=1,nvirt_aos

                                    energy_tmp = ccsd_doubles%val(a,b,i_eos,j_eos) &
                                               & * ccsdpt_doubles%val(a,b,i_eos,j_eos)
                                    energy_res_cou = energy_res_cou + energy_tmp

                                 end do
                              end do

                           end do ido_pair_doubles_cou
                        end do jdo_pair_doubles_cou
                        !$OMP END PARALLEL DO

    ! reorder from (a,b,i,j) to (a,b,j,i)
    call array4_reorder(ccsd_doubles,[1,2,4,3])

                        !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,i_eos,j,j_eos,a,b,energy_tmp),&
                        !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nocc_eos,nvirt_aos,&
                        !$OMP PairFragment,dopair_occ),REDUCTION(+:energy_res_exc)
  jdo_pair_doubles_exc: do j=1,nocc_eos
                        j_eos = PairFragment%idxo(j)
     ido_pair_doubles_exc: do i=1,nocc_eos
                           i_eos = PairFragment%idxo(i)

                              if (.not. dopair_occ(i,j)) cycle ido_pair_doubles_exc

                              do b=1,nvirt_aos
                                 do a=1,nvirt_aos

                                    energy_tmp = ccsd_doubles%val(a,b,i_eos,j_eos) &
                                               & * ccsdpt_doubles%val(a,b,i_eos,j_eos)
                                    energy_res_exc = energy_res_exc + energy_tmp

                                 end do
                              end do

                           end do ido_pair_doubles_exc
                        end do jdo_pair_doubles_exc
                        !$OMP END PARALLEL DO

    ! get total fourth--order energy contribution
    PairFragment%energies(FRAGMODEL_OCCpT4) = 4.0E0_realk * energy_res_cou - 2.0E0_realk * energy_res_exc

    ! insert into occ. part. scheme part
    PairFragment%energies(FRAGMODEL_OCCpT) = PairFragment%energies(FRAGMODEL_OCCpT) + PairFragment%energies(FRAGMODEL_OCCpT4)

    ! *********************************
    ! do unoccupied partitioning scheme
    ! *********************************

    ! initially, reorder ccsd_doubles and ccsdpt_doubles
    ! ccsd_doubles from from (a,b,j,i) sequence to (j,i,a,b) sequence
    ! ccsdpt_doubles from from (a,b,i,j) sequence to (i,j,a,b) sequence
    call array4_reorder(ccsd_doubles,[3,4,1,2])
    call array4_reorder(ccsdpt_doubles,[3,4,1,2])

    energy_res_cou = 0.0E0_realk
    energy_res_exc = 0.0E0_realk

                        !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(a,a_eos,b,b_eos,i,j,energy_tmp),&
                        !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nvirt_eos,nocc_aos,&
                        !$OMP PairFragment,dopair_virt),REDUCTION(+:energy_res_exc)
  bdo_pair_doubles_exc: do b=1,nvirt_eos
                        b_eos = PairFragment%idxu(b)
     ado_pair_doubles_exc: do a=1,nvirt_eos
                           a_eos = PairFragment%idxu(a)

                              if (.not. dopair_virt(a,b)) cycle ado_pair_doubles_exc
    
                              do j=1,nocc_aos
                                 do i=1,nocc_aos

                                    energy_tmp = ccsd_doubles%val(i,j,a_eos,b_eos) &
                                               & * ccsdpt_doubles%val(i,j,a_eos,b_eos)
                                    energy_res_exc = energy_res_exc + energy_tmp

                                 end do
                              end do
    
                           end do ado_pair_doubles_exc
                        end do bdo_pair_doubles_exc
                        !$OMP END PARALLEL DO

    ! reorder form (j,i,a,b) to (i,j,a,b)
    call array4_reorder(ccsd_doubles,[2,1,3,4])

                        !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(a,a_eos,b,b_eos,i,j,energy_tmp),&
                        !$OMP SHARED(ccsd_doubles,ccsdpt_doubles,nvirt_eos,nocc_aos,&
                        !$OMP PairFragment,dopair_virt),REDUCTION(+:energy_res_cou)
  bdo_pair_doubles_cou: do b=1,nvirt_eos
                        b_eos = PairFragment%idxu(b)
     ado_pair_doubles_cou: do a=1,nvirt_eos
                           a_eos = PairFragment%idxu(a)

                              if (.not. dopair_virt(a,b)) cycle ado_pair_doubles_cou

                              do j=1,nocc_aos
                                 do i=1,nocc_aos

                                    energy_tmp = ccsd_doubles%val(i,j,a_eos,b_eos) &
                                               & * ccsdpt_doubles%val(i,j,a_eos,b_eos)
                                    energy_res_cou = energy_res_cou + energy_tmp

                                 end do
                              end do

                           end do ado_pair_doubles_cou
                        end do bdo_pair_doubles_cou
                        !$OMP END PARALLEL DO

    ! get total fourth--order energy contribution
    PairFragment%energies(FRAGMODEL_VIRTpT4) = 4.0E0_realk * energy_res_cou - 2.0E0_realk * energy_res_exc

    ! insert into virt. part. scheme part
    PairFragment%energies(FRAGMODEL_VIRTpT) = PairFragment%energies(FRAGMODEL_VIRTpT) + PairFragment%energies(FRAGMODEL_VIRTpT4)

    ! ******************************
    !   done with E[4] energy part
    ! ******************************

    ! now release logical pair arrays
    call mem_dealloc(dopair_occ)
    call mem_dealloc(dopair_virt)

    ! ******************************************************************
    ! ************** done w/ energy for pair fragments *****************
    ! ******************************************************************

  end subroutine ccsdpt_energy_e4_pair

  !> \brief: calculate E[4] contribution to ccsd(t) energy correction for full molecule calculation
  !> \author: Janus Juul Eriksen
  !> \date: February 2013
  subroutine ccsdpt_energy_e4_full(nocc,nvirt,natoms,offset,ccsd_doubles,ccsdpt_doubles,occ_orbitals,&
                           & eccsdpt_matrix_cou,eccsdpt_matrix_exc,ccsdpt_e4)

    implicit none

    !> ccsd and ccsd(t) doubles amplitudes
    type(array4), intent(inout) :: ccsd_doubles, ccsdpt_doubles
    !> dimensions
    integer, intent(in) :: nocc, nvirt, natoms, offset
    !> occupied orbital information
    type(decorbital), dimension(nocc+offset), intent(inout) :: occ_orbitals
    !> etot
    real(realk), intent(inout) :: ccsdpt_e4
    real(realk), dimension(natoms,natoms), intent(inout) :: eccsdpt_matrix_cou, eccsdpt_matrix_exc
    !> integers
    integer :: i,j,a,b,atomI,atomJ
    !> energy reals
    real(realk) :: energy_tmp, energy_res_cou, energy_res_exc

    ! *************************************************************
    ! ************** do energy for full molecule ******************
    ! *************************************************************

    ! ***********************
    !   do E[4] energy part
    ! ***********************

    energy_res_cou = 0.0E0_realk
    energy_res_exc = 0.0E0_realk
    ccsdpt_e4 = 0.0E0_realk

    ! ***note: we only run over nval (which might be equal to nocc_tot if frozencore = .false.)
    ! so we only assign orbitals for the space in which the core orbitals (the offset) are omited

    !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,atomJ,a,b,energy_tmp),REDUCTION(+:energy_res_cou),&
    !$OMP REDUCTION(+:eccsdpt_matrix_cou),SHARED(ccsd_doubles,ccsdpt_doubles,nocc,nvirt,occ_orbitals,offset)
    do j=1,nocc
    atomJ = occ_orbitals(j+offset)%CentralAtom
       do i=1,nocc
       atomI = occ_orbitals(i+offset)%CentralAtom

          do b=1,nvirt
             do a=1,nvirt

                energy_tmp = ccsd_doubles%val(a,b,i,j) * ccsdpt_doubles%val(a,b,i,j)
                eccsdpt_matrix_cou(AtomI,AtomJ) = eccsdpt_matrix_cou(AtomI,AtomJ) + energy_tmp
                energy_res_cou = energy_res_cou + energy_tmp

             end do
          end do

       end do
    end do
    !$OMP END PARALLEL DO

    ! reorder from (a,b,i,j) to (a,b,j,i)
    call array4_reorder(ccsd_doubles,[1,2,4,3])

    !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,atomI,j,atomJ,a,b,energy_tmp),REDUCTION(+:energy_res_exc),&
    !$OMP REDUCTION(+:eccsdpt_matrix_exc),SHARED(ccsd_doubles,ccsdpt_doubles,nocc,nvirt,occ_orbitals,offset)
    do j=1,nocc
    atomJ = occ_orbitals(j+offset)%CentralAtom
       do i=1,nocc
       atomI = occ_orbitals(i+offset)%CentralAtom

          do b=1,nvirt
             do a=1,nvirt

                energy_tmp = ccsd_doubles%val(a,b,i,j) * ccsdpt_doubles%val(a,b,i,j)
                eccsdpt_matrix_exc(AtomI,AtomJ) = eccsdpt_matrix_exc(AtomI,AtomJ) + energy_tmp
                energy_res_exc = energy_res_exc + energy_tmp

             end do
          end do

       end do
    end do
    !$OMP END PARALLEL DO

    ! get total fourth--order energy contribution
    eccsdpt_matrix_cou = 4.0E0_realk * eccsdpt_matrix_cou - 2.0E0_realk * eccsdpt_matrix_exc
    ccsdpt_e4 = 4.0E0_realk * energy_res_cou - 2.0E0_realk * energy_res_exc

    ! for the e4 pair fragment energy matrix,
    ! we only consider pairs IJ where J>I; thus, move contributions

    do AtomJ=1,natoms
       do AtomI=AtomJ+1,natoms

          eccsdpt_matrix_cou(AtomI,AtomJ) = eccsdpt_matrix_cou(AtomI,AtomJ) &
                                              & + eccsdpt_matrix_cou(AtomJ,AtomI)
       end do
    end do



    ! ******************************************************************
    ! ************** done w/ energy for full molecule ******************
    ! ******************************************************************

  end subroutine ccsdpt_energy_e4_full

  !> \brief: print out E[4] fragment and pair interaction contribution to 
  !>         ccsd(t) energy correction for full molecule calculation
  !> \author: Janus Juul Eriksen
  !> \date: February 2013
  subroutine print_e4_full(natoms,e4_matrix,orbitals_assigned,distancetable)

    implicit none

    !> number of atoms in molecule
    integer, intent(in) :: natoms
    !> matrices containing E[4] energies and interatomic distances
    real(realk), dimension(natoms,natoms), intent(in) :: e4_matrix, distancetable
    !> vector handling how the orbitals are assigned?
    logical, dimension(natoms), intent(inout) :: orbitals_assigned
    !> loop counters
    integer :: i,j
!    use the one in lsutil/fundamental.f90
!    real(realk), parameter :: bohr_to_angstrom = 0.5291772083E0_realk

    ! print out fragment energies

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') '***************************************************************'
    write(DECinfo%output,'(1X,a)') '*                         E[4] energies                       *'
    write(DECinfo%output,'(1X,a)') '***************************************************************'
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(8X,a)') '-- Atomic fragment energies (fourth--order E[4])'
    write(DECinfo%output,'(8X,a)') '------    --------------------'
    write(DECinfo%output,'(8X,a)') ' Atom            Energy '
    write(DECinfo%output,'(8X,a)') '------    --------------------'
    write(DECinfo%output,*)

    do i=1,natoms

       if (orbitals_assigned(i)) then

          write(DECinfo%output,'(1X,a,i6,4X,g20.10)') '#SING#', i, e4_matrix(i,i)

       end if

    end do

    ! now print out pair interaction energies

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(8X,a)') '-- Pair interaction energies (fourth--order E[4])     '
    write(DECinfo%output,'(8X,a)') '------    ------    ----------    --------------------'
    write(DECinfo%output,'(8X,a)') '   P         Q        R(Ang)              E(PQ)       '
    write(DECinfo%output,'(8X,a)') '------    ------    ----------    --------------------'
    write(DECinfo%output,*)

    do j=1,natoms
       do i=j+1,natoms

          ! write increments only if pair interaction energy is nonzero
          if( orbitals_assigned(i) .and. orbitals_assigned(j) ) then

             write(DECinfo%output,'(1X,a,i6,4X,i6,4X,g10.4,4X,g20.10)') '#PAIR#',j,i,&
                  & bohr_to_angstrom*distancetable(i,j), e4_matrix(i,j)

          end if

       end do
    end do


  end subroutine print_e4_full



  !> \brief: calculate E[5] contribution to ccsd(t) energy correction for full molecule calculation
  !> \author: Janus Juul Eriksen
  !> \date: February 2013
  subroutine ccsdpt_energy_e5_full(nocc,nvirt,natoms,offset,ccsd_singles,ccsdpt_singles,&
                             & occ_orbitals,unocc_orbitals,e5_matrix,ccsdpt_e5)

    implicit none

    !> ccsd and ccsd(t) singles amplitudes
    type(array2), intent(inout) :: ccsd_singles, ccsdpt_singles
    !> dimensions
    integer, intent(in) :: nocc, nvirt, natoms, offset
    !> occupied orbital information
    type(decorbital), dimension(nocc+offset), intent(inout) :: occ_orbitals
    !> virtual orbital information
    type(decorbital), dimension(nvirt), intent(inout) :: unocc_orbitals
    !> etot
    real(realk), intent(inout) :: ccsdpt_e5
    real(realk), dimension(natoms,natoms), intent(inout) :: e5_matrix
    !> integers
    integer :: i,a,AtomI,AtomA
    !> tmp energy real
    real(realk) :: energy_tmp

    ! ***********************
    !   do E[5] energy part
    ! ***********************

    ccsdpt_e5 = 0.0E0_realk

    !$OMP PARALLEL DO DEFAULT(NONE),PRIVATE(i,a,energy_tmp,AtomI,AtomA),&
    !$OMP SHARED(ccsd_singles,ccsdpt_singles,nocc,nvirt,offset,occ_orbitals,unocc_orbitals),&
    !$OMP REDUCTION(+:ccsdpt_e5),REDUCTION(+:e5_matrix)
    do i=1,nocc
    AtomI = occ_orbitals(i+offset)%CentralAtom
       do a=1,nvirt
       AtomA = unocc_orbitals(a)%CentralAtom

           energy_tmp = ccsd_singles%val(a,i) * ccsdpt_singles%val(a,i)
           e5_matrix(AtomA,AtomI) = e5_matrix(AtomA,AtomI) + energy_tmp
           ccsdpt_e5 = ccsdpt_e5 + energy_tmp

       end do
    end do
    !$OMP END PARALLEL DO

    ! get total fifth-order energy correction
    e5_matrix = 2.0E0_realk * e5_matrix
    ccsdpt_e5 = 2.0E0_realk * ccsdpt_e5

    ! ******************************
    !   done with E[5] energy part
    ! ******************************

  end subroutine ccsdpt_energy_e5_full

  !> \brief: print out fifth-order pair interaction energies for full molecule calculation 
  !> \author: Janus Juul Eriksen
  !> \date: February 2013
  subroutine print_e5_full(natoms,e5_matrix,orbitals_assigned,distancetable)

    implicit none

    !> number of atoms in molecule
    integer, intent(in) :: natoms
    !> matrices containing E[4] energies and interatomic distances
    real(realk), dimension(natoms,natoms), intent(in) :: e5_matrix, distancetable
    !> vector handling how the orbitals are assigned?
    logical, dimension(natoms), intent(inout) :: orbitals_assigned
    !> loop counters
    integer :: i,a
!    real(realk), parameter :: bohr_to_angstrom = 0.5291772083E0_realk

    ! print out fragment energies

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') '***************************************************************'
    write(DECinfo%output,'(1X,a)') '*                         E[5] energies                       *'
    write(DECinfo%output,'(1X,a)') '***************************************************************'
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(8X,a)') '-- Pair fragment energies (fifth--order E[5])          '
    write(DECinfo%output,'(9X,a)') '-------    ------    ----------    --------------------'
    write(DECinfo%output,'(9X,a)') 'P(virt)    Q(occ)      R(Ang)           E(PQ)          '
    write(DECinfo%output,'(9X,a)') '-------    ------    ----------    --------------------'
    write(DECinfo%output,*)

    ! the total singles energy must result from an unrestricted summation over all occ and virt indices
    ! as we are only interested in general orbital interactions and hence not the nature (occ/virt)
    ! of the individual orbitals

    do i=1,natoms
       do a=1,natoms

          ! write increments only if pair interaction energy is nonzero
          if( orbitals_assigned(i) .and. orbitals_assigned(a) ) then

             if (i .eq. a) then
                write(DECinfo%output,'(1X,a,i7,4X,i6,4X,g10.4,4X,g20.10)') '#PAIR#',a,i,&
                     &0.000, e5_matrix(a,i)

             else

                write(DECinfo%output,'(1X,a,i7,4X,i6,4X,g10.4,4X,g20.10)') '#PAIR#',a,i,&
                     &bohr_to_angstrom*distancetable(a,i), e5_matrix(a,i)

             end if

          end if

       end do
    end do

  end subroutine print_e5_full


  !> \brief Get MO integrals for CCSD(T) (in canonical basis), see integral storing order below.
  !> \author Janus Eriksen and Kasper Kristensen
  !> \date September-October 2012
  subroutine get_CCSDpT_integrals(MyLsitem,nbasis,nocc,nvirt,Cocc,Cvirt,JAIK,ABIJ,CBAI)

    implicit none

    !> Integral info
    type(lsitem), intent(inout) :: mylsitem
    !> Number of basis functions
    integer,intent(in) :: nbasis
    !> Number of occupied orbitals
    integer,intent(in) :: nocc
    !> Number of virtual orbitals
    integer,intent(in) :: nvirt
    !> Occupied MO coefficients
    real(realk), dimension(nbasis,nocc),intent(in) :: Cocc
    !> Virtual MO coefficients
    real(realk), dimension(nbasis,nvirt),intent(in) :: Cvirt
    ! JIAK: Integrals (AI|JK) in the order (J,A,I,K)
    type(array4), intent(inout) :: JAIK
    ! ABIJ: Integrals (AI|BJ) in the order (A,B,I,J)
    type(array4), intent(inout) :: ABIJ
    ! CBAI: Integrals (AI|BC) in the order (C,B,A,I)
    type(array), intent(inout) :: CBAI
    integer :: gammadim, alphadim,iorb
    integer :: alphaB,gammaB,dimAlpha,dimGamma,idx
    real(realk),pointer :: tmp1(:),tmp2(:),tmp3(:)
    integer(kind=long) :: size1,size2,size3
    integer :: GammaStart, GammaEnd, AlphaStart, AlphaEnd,m,k,n,i,dims(4),order(4)
    logical :: FullRHS,doscreen
    real(realk) :: tcpu, twall
    real(realk),pointer :: CoccT(:,:), CvirtT(:,:)
    type(array4) :: JAIB
    integer :: MaxActualDimAlpha,nbatchesAlpha
    integer :: MaxActualDimGamma,nbatchesGamma
    type(batchtoorb), pointer :: batch2orbAlpha(:)
    type(batchtoorb), pointer :: batch2orbGamma(:)
    integer, pointer :: orb2batchAlpha(:), batchdimAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
    integer, pointer :: orb2batchGamma(:), batchdimGamma(:), batchsizeGamma(:), batchindexGamma(:)
    TYPE(DECscreenITEM)   :: DecScreen
    ! distribution stuff needed for mpi parallelization
    integer, pointer :: distribution(:)
    Character            :: intSpec(5)
    integer :: myload

    ! Lots of timings
    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Integral screening?
    doscreen = mylsitem%setting%scheme%cs_screen.OR.&
         & mylsitem%setting%scheme%ps_screen

    ! allocate arrays to update during integral loop 
    ! **********************************************
    
    ! note 1: this must be done before call to get_optimal_batch_sizes_ccsdpt_integrals
    ! note 2: these integrals will be reordered into the output structures

    ! JAIK: Integrals (AI|KJ) in the order (J,A,I,K)
    dims = [nocc,nvirt,nocc,nocc]
    JAIK = array4_init_standard(dims)

    ! JAIB: Integrals (AI|BJ) in the order (J,A,I,B)
    dims = [nocc,nvirt,nocc,nvirt]
    JAIB = array4_init_standard(dims)

    ! CBAI: Integrals (AB|IC) in the order (C,B,A,I)
    dims = [nvirt,nvirt,nvirt,nocc]

#ifdef VAR_MPI

    CBAI = array_init(dims,4,TILED_DIST,ALL_ACCESS,[nvirt,nvirt,nvirt,1])
    call array_zero_tiled_dist(CBAI)

#else

    CBAI = array_init(dims,4)
    call array_zero(CBAI)

#endif

    ! For efficiency when calling dgemm, save transposed matrices
    call mem_alloc(CoccT,nocc,nbasis)
    call mem_alloc(CvirtT,nvirt,nbasis)
    call mat_transpose(nbasis,nocc,1.0E0_realk,Cocc,0.0E0_realk,CoccT)
    call mat_transpose(nbasis,nvirt,1.0E0_realk,Cvirt,0.0E0_realk,CvirtT)

    ! Determine optimal batchsizes and corresponding sizes of arrays
    call get_optimal_batch_sizes_ccsdpt_integrals(mylsitem,nbasis,nocc,nvirt,alphadim,gammadim,&
         & size1,size2,size3)


    ! ************************************************
    ! * Determine batch information for Gamma batch  *
    ! ************************************************

    ! Orbital to batch information
    ! ----------------------------
    call mem_alloc(orb2batchGamma,nbasis)
    call build_batchesofAOS(DECinfo%output,mylsitem%setting,gammadim,&
         & nbasis,MaxActualDimGamma,batchsizeGamma,batchdimGamma,batchindexGamma,&
         & nbatchesGamma,orb2BatchGamma,'R')

    write(DECinfo%output,*) 'BATCH: Number of Gamma batches   = ', nbatchesGamma

    ! Translate batchindex to orbital index
    ! -------------------------------------
    call mem_alloc(batch2orbGamma,nbatchesGamma)

    do idx=1,nbatchesGamma

       call mem_alloc(batch2orbGamma(idx)%orbindex,batchdimGamma(idx) )
       batch2orbGamma(idx)%orbindex = 0
       batch2orbGamma(idx)%norbindex = 0

    end do

    do iorb=1,nbasis

       idx = orb2batchGamma(iorb)
       batch2orbGamma(idx)%norbindex = batch2orbGamma(idx)%norbindex+1
       K = batch2orbGamma(idx)%norbindex
       batch2orbGamma(idx)%orbindex(K) = iorb

    end do

    ! ************************************************
    ! * Determine batch information for Alpha batch  *
    ! ************************************************

    ! Orbital to batch information
    ! ----------------------------
    call mem_alloc(orb2batchAlpha,nbasis)
    call build_batchesofAOS(DECinfo%output,mylsitem%setting,alphadim,&
         & nbasis,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,&
         & nbatchesAlpha,orb2BatchAlpha,'R')

    write(DECinfo%output,*) 'BATCH: Number of Alpha batches   = ', nbatchesAlpha

    ! Translate batchindex to orbital index
    ! -------------------------------------
    call mem_alloc(batch2orbAlpha,nbatchesAlpha)

    do idx=1,nbatchesAlpha

       call mem_alloc(batch2orbAlpha(idx)%orbindex,batchdimAlpha(idx) )
       batch2orbAlpha(idx)%orbindex = 0
       batch2orbAlpha(idx)%norbindex = 0

    end do

    do iorb=1,nbasis

       idx = orb2batchAlpha(iorb)
       batch2orbAlpha(idx)%norbindex = batch2orbAlpha(idx)%norbindex+1
       K = batch2orbAlpha(idx)%norbindex
       batch2orbAlpha(idx)%orbindex(K) = iorb

    end do

    ! Set integral info
    ! *****************
    INTSPEC(1)='R' !R = Regular Basis set on the 1th center 
    INTSPEC(2)='R' !R = Regular Basis set on the 2th center 
    INTSPEC(3)='R' !R = Regular Basis set on the 3th center 
    INTSPEC(4)='R' !R = Regular Basis set on the 4th center 
    INTSPEC(5)='C' !C = Coulomb operator
    call II_precalc_DECScreenMat(DecScreen,DECinfo%output,6,mylsitem%setting,&
            & nbatchesAlpha,nbatchesGamma,INTSPEC)

    if (doscreen) then

       call II_getBatchOrbitalScreen(DecScreen,mylsitem%setting,&
            & nbasis,nbatchesAlpha,nbatchesGamma,&
            & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
            & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)

    end if

    FullRHS = (nbatchesGamma .eq. 1) .and. (nbatchesAlpha .eq. 1)


    ! Allocate array for AO integrals
    ! *******************************
    call mem_alloc(tmp1,size1)
    call mem_alloc(tmp2,size2)
    call mem_alloc(tmp3,size3)

#ifdef VAR_MPI

    ! alloc distribution array
    nullify(distribution)
    call mem_alloc(distribution,nbatchesGamma*nbatchesAlpha)

    ! init distribution
    distribution = 0
    myload = 0
    call distribute_mpi_jobs(distribution,nbatchesAlpha,nbatchesGamma,&
    &batchdimAlpha,batchdimGamma,myload,infpar%lg_nodtot,infpar%lg_mynum)

#endif

    ! Start looping over gamma and alpha batches and calculate integrals
    ! ******************************************************************

    BatchGamma: do gammaB = 1,nbatchesGamma  ! AO batches
       dimGamma = batchdimGamma(gammaB)                           ! Dimension of gamma batch
       GammaStart = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
       GammaEnd = batch2orbGamma(gammaB)%orbindex(dimGamma)       ! Last index in gamma batch


       BatchAlpha: do alphaB = 1,nbatchesAlpha  ! AO batches
          dimAlpha = batchdimAlpha(alphaB)                                ! Dimension of alpha batch
          AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)                 ! First index in alpha batch
          AlphaEnd = batch2orbAlpha(alphaB)%orbindex(dimAlpha)            ! Last index in alpha batch

#ifdef VAR_MPI

          ! distribute tasks
          if (distribution((alphaB-1)*nbatchesGamma+gammaB) .ne. infpar%lg_mynum) then

             cycle BatchAlpha

          end if

!          write (DECinfo%output, '("Rank(T) ",I3," starting job (",I3,"/",I3,",",I3,"/",I3,")")'),infpar%lg_mynum,alphaB,&
!                          &nbatchesAlpha,gammaB,nbatchesGamma

#endif

          if (doscreen) mylsitem%setting%LST_GAB_LHS => DECSCREEN%masterGabLHS
          if (doscreen) mylsitem%setting%LST_GAB_RHS => DECSCREEN%batchGab(alphaB,gammaB)%p


          ! Get (beta delta | alphaB gammaB) integrals using (beta,delta,alphaB,gammaB) ordering
          ! ************************************************************************************
          call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output, &
               & mylsitem%setting,tmp1,batchindexAlpha(alphaB),batchindexGamma(gammaB),&
               & batchsizeAlpha(alphaB),batchsizeGamma(gammaB),nbasis,nbasis,dimAlpha,dimGamma,&
               & FullRHS,INTSPEC)

          ! tmp2(delta,alphaB,gammaB;A) = sum_{beta} [tmp1(beta;delta,alphaB,gammaB)]^T Cvirt(beta,A)
          m = nbasis*dimGamma*dimAlpha
          k = nbasis
          n = nvirt
!          call dec_simple_dgemm(m,k,n,tmp1,CvirtT,tmp2,'T','T')
          call dgemm('T','N',m,n,k,1.0E0_realk,tmp1,k,Cvirt,k,0.0E0_realk,tmp2,m)

          ! tmp3(B;alphaB,gammaB,A) = sum_{delta} CvirtT(B,delta) tmp2(delta;alphaB,gammaB,A)
          m = nvirt
          k = nbasis
          n = dimAlpha*dimGamma*nvirt
!          call dec_simple_dgemm(m,k,n,CvirtT,tmp2,tmp3,'N','N')
          call dgemm('N','N',m,n,k,1.0E0_realk,CvirtT,m,tmp2,k,0.0E0_realk,tmp3,m)

          ! tmp1(I;,alphaB,gammaB,A) = sum_{delta} CoccT(I,delta) tmp2(delta,alphaB,gammaB,A)
          m = nocc
          k = nbasis
          n = dimAlpha*dimGamma*nvirt
!          call dec_simple_dgemm(m,k,n,CoccT,tmp2,tmp1,'N','N')
          call dgemm('N','N',m,n,k,1.0E0_realk,CoccT,m,tmp2,k,0.0E0_realk,tmp1,m)

          ! Reorder: tmp1(I,alphaB;gammaB,A) --> tmp2(gammaB,A;I,alphaB)
          m = nocc*dimAlpha
          n = dimGamma*nvirt
          call mat_transpose(m,n,1.0E0_realk,tmp1,0.0E0_realk,tmp2)

          ! tmp1(J;A,I,alphaB) = sum_{gamma in gammaB} CoccT(J,gamma) tmp2(gamma,A,I,alphaB)
          m = nocc
          k = dimGamma
          n = nvirt*nocc*dimAlpha
!          call dec_simple_dgemm(m,k,n,CoccT(1:nocc,GammaStart:GammaEnd),tmp2,tmp1,'N','N')
          call dgemm('N','N',m,n,k,1.0E0_realk,CoccT(1:nocc,GammaStart:GammaEnd),m,tmp2,k,0.0E0_realk,tmp1,m)

          ! JAIK(J,A,I;K) += sum_{alpha in alphaB} tmp1(J,A,I,alpha) Cocc(alpha,K)
          m = nvirt*nocc**2
          k = dimAlpha
          n = nocc
!          call dec_simple_dgemm_update(m,k,n,tmp1,&
!                                     & CoccT(1:nocc,AlphaStart:AlphaEnd),JAIK%val,'N','T')
!          call dgemm('N','N',m,n,k,1.0E0_realk,tmp1,m,Cocc(AlphaStart:AlphaEnd,1:nocc),k,1.0E0_realk,JAIK%val,m)
          call dgemm('N','N',m,n,k,1.0E0_realk,tmp1,m,Cocc(AlphaStart:,:),nbasis-AlphaStart+1,1.0E0_realk,JAIK%val,m)

          ! JAIB(J,A,I;B) += sum_{alpha in alphaB} tmp1(J,A,I,alpha) Cvirt(alpha,B)
          m = nvirt*nocc**2
          k = dimAlpha
          n = nvirt
!          call dec_simple_dgemm_update(m,k,n,tmp1,&
!                                     & CvirtT(1:nvirt,AlphaStart:AlphaEnd),JAIB%val,'N','T')
!          call dgemm('N','N',m,n,k,1.0E0_realk,tmp1,m,Cvirt(AlphaStart:AlphaEnd,1:nvirt),k,1.0E0_realk,JAIB%val,m)
          call dgemm('N','N',m,n,k,1.0E0_realk,tmp1,m,Cvirt(AlphaStart:,:),nbasis-AlphaStart+1,1.0E0_realk,JAIB%val,m)

          ! Reorder: tmp3(B,alphaB;gammaB,A) --> tmp1(gammaB,A;B,alphaB)
          m = nvirt*dimAlpha
          n = dimGamma*nvirt
          call mat_transpose(m,n,1.0E0_realk,tmp3,0.0E0_realk,tmp1)

          ! tmp3(C;A,B,alphaB) = sum_{gamma in gammaB} CvirtT(C,gamma) tmp1(gamma,A,B,alphaB)
          m = nvirt
          k = dimGamma
          n = dimAlpha*nvirt**2
!          call dec_simple_dgemm(m,k,n,CvirtT(1:nvirt,GammaStart:GammaEnd),tmp1,tmp3,'N','N')
          call dgemm('N','N',m,n,k,1.0E0_realk,CvirtT(:,GammaStart:),m,tmp1,k,0.0E0_realk,tmp3,m)

          ! reorder tmp1 and do CBAI(B,A,C,I) += sum_{i in IB} tmp1(B,A,C,i)
          m = nvirt**3
          k = dimAlpha
          n = 1

#ifdef VAR_MPI

          do i=1,nocc

             ! tmp1(C,A,B,i) = sum_{alpha in alphaB} tmp3(C,A,B,alpha) Cocc(alpha,i)
!             call dec_simple_dgemm(m,k,n,tmp3,CoccT(i,AlphaStart:AlphaEnd),tmp1,'N','T')
             call dgemm('N','N',m,n,k,1.0E0_realk,tmp3,m,Cocc(AlphaStart:,i),nbasis-AlphaStart+1,0.0E0_realk,tmp1,m)

             ! *** tmp1 corresponds to (AB|iC) in Mulliken notation. Noting that the v³o integrals
             ! are normally written as g_{AIBC}, we may also write this Mulliken integral (with substitution
             ! of dummy indices A=B, B=C, and C=A) as (BC|IA). In order to align with the CBAI order of
             ! ccsd(t) driver routine, we reorder as:
             ! (BC|IA) --> (CB|AI), i.e., tmp1(C,A,B,i) = ABCI(A,B,C,i) (norm. notat.) --> 
             !                                            tmp1(C,B,A,i) (norm. notat.) = tmp1(B,A,C,i) (notat. herein)
             ! 
             ! next, we accumulate
             ! CBAI(B,A,C,I) += sum_{i in IB} tmp1(B,A,C,i)

             call array_reorder_3d(1.0E0_realk,tmp1,nvirt,nvirt,nvirt,[3,2,1],0.0E0_realk,tmp2)

             call array_accumulate_tile(CBAI,i,tmp2,nvirt**3)

          end do

#else

          do i=1,nocc

             ! for description, see mpi section above
!             call dec_simple_dgemm(m,k,n,tmp3,CoccT(i,AlphaStart:AlphaEnd),tmp1,'N','T')
             call dgemm('N','N',m,n,k,1.0E0_realk,tmp3,m,Cocc(AlphaStart:,i),nbasis-AlphaStart+1,0.0E0_realk,tmp1,m)

             call array_reorder_3d(1.0E0_realk,tmp1,nvirt,nvirt,nvirt,[3,2,1],1.0E0_realk,CBAI%elm4(:,:,:,i))

          end do

#endif

       end do BatchAlpha
    end do BatchGamma

#ifdef VAR_MPI

    if (infpar%lg_nodtot .gt. 1) then

       ! now, reduce o^2v^2 and o^3v integrals onto master
       call lsmpi_allreduce(JAIB%val,nocc,nvirt,nocc,nvirt,infpar%lg_comm)
       call lsmpi_allreduce(JAIK%val,nocc,nvirt,nocc,nocc, infpar%lg_comm) 

    end if

    ! dealloc distribution array
    call mem_dealloc(distribution)

#endif

    ! free stuff
    ! **********
    call mem_dealloc(tmp1)
    call mem_dealloc(tmp2)
    call mem_dealloc(tmp3)
    call free_decscreen(DECSCREEN)
    call mem_dealloc(CoccT)
    call mem_dealloc(CvirtT)
    call mem_dealloc(orb2batchGamma)
    call mem_dealloc(batchdimGamma)
    call mem_dealloc(batchsizeGamma)
    call mem_dealloc(batchindexGamma)
    do idx=1,nbatchesGamma
       call mem_dealloc(batch2orbGamma(idx)%orbindex)
    end do
    call mem_dealloc(batch2orbGamma)
    call mem_dealloc(orb2batchAlpha)
    call mem_dealloc(batchdimAlpha)
    call mem_dealloc(batchsizeAlpha)
    call mem_dealloc(batchindexAlpha)
    do idx=1,nbatchesAlpha
       call mem_dealloc(batch2orbAlpha(idx)%orbindex)
       batch2orbAlpha(idx)%orbindex => null()
    end do
    call mem_dealloc(batch2orbAlpha)
    nullify(mylsitem%setting%LST_GAB_LHS)
    nullify(mylsitem%setting%LST_GAB_RHS)

    ! finally, reorder JAIB to final output
    ! *********************************************

    ! ** JAIB corresponds to (AI|BJ) in Mulliken notation. Noting that the v²o² integrals
    ! are normally written as g_{AIBJ} = g_{BJAI}, we may also write this Mulliken integral (with substitution
    ! of dummy indices A=B and B=A) as (BI|AJ). In order to align with the ABIJ order of
    ! ccsd(t) driver routine, we reorder as:
    ! (BI|AJ) --> (AB|IJ), i.e., JAIB(J,A,I,B) = JBIA(J,B,I,A) (norm. notat.) --> 
    !                                            BAIJ(B,A,I,J) (norm. notat.) =
    !                                            ABIJ(A,B,I,J) (notat. herein)

    order = [2,4,3,1]
    dims = [nvirt,nvirt,nocc,nocc]
    ABIJ = array4_init_standard(dims)
    
    call array_reorder_4d(1.0E0_realk,JAIB%val,JAIB%dims(1),JAIB%dims(2),&
         & JAIB%dims(3),JAIB%dims(4),order,0.0E0_realk,ABIJ%val)
    
    call array4_free(JAIB)

    call LSTIMER('CCSD(T) INT',tcpu,twall,DECinfo%output)

  end subroutine get_CCSDpT_integrals

  !> \brief Get optimal batch sizes to be used in get_CCSDpT_integrals
  !> using the available memory.
  !> \author Kasper Kristensen & Janus Eriksen
  !> \date September 2011, rev. October 2012
  subroutine get_optimal_batch_sizes_ccsdpt_integrals(mylsitem,nbasis,nocc,nvirt,alphadim,gammadim,&
     & size1,size2,size3)

    implicit none
  
    !> Integral info
    type(lsitem), intent(inout) :: mylsitem
    !> Number of AO basis functions
    integer,intent(in) :: nbasis
    !> Number of occupied (AOS) orbitals
    integer,intent(in) :: nocc
    !> Number of virt (AOS) orbitals
    integer,intent(in) :: nvirt
    !> Max size for AO alpha batch
    integer,intent(inout) :: alphadim
    !> Max size for AO gamma batch
    integer,intent(inout) :: gammadim
    !> Dimension of temporary array 1
    integer(kind=long),intent(inout) :: size1
    !> Dimension of temporary array 2
    integer(kind=long),intent(inout) :: size2
    !> Dimension of temporary array 3
    integer(kind=long),intent(inout) :: size3
    !> memory reals
    real(realk) :: MemoryNeeded, MemoryAvailable
    integer :: MaxAObatch, MinAOBatch, AlphaOpt, GammaOpt,alpha,gamma


    ! Memory currently available
    ! **************************
    call get_currently_available_memory(MemoryAvailable)
    ! Note: We multiply by 85 % to be on the safe side!
    MemoryAvailable = 0.85*MemoryAvailable
  
  
  
    ! Maximum and minimum possible batch sizes
    ! ****************************************
  
    ! The largest possible AO batch is the number of basis functions
    MaxAObatch = nbasis
  
    ! The smallest possible AO batch depends on the basis set
    ! (More precisely, if all batches are made as small as possible, then the
    !  call below determines the largest of these small batches).
    call determine_maxBatchOrbitalsize(DECinfo%output,mylsitem%setting,MinAObatch,'R')
  
  
    ! Initialize batch sizes to be the minimum possible and then start increasing sizes below
    AlphaDim=MinAObatch
    GammaDim=MinAObatch
  
  
    ! Gamma batch size
    ! =================================
    GammaLoop: do gamma = MaxAObatch,MinAOBatch,-1
  
       call get_max_arraysizes_for_ccsdpt_integrals(alphaDim,gamma,nbasis,nocc,nvirt,&
            & size1,size2,size3,MemoryNeeded)
  
       if(MemoryNeeded < MemoryAvailable .or. (gamma==minAObatch) ) then
          GammaOpt = gamma
          exit
       end if
  
    end do GammaLoop
  
    ! If gamma batch size was set manually we use that value instead
    if(DECinfo%ccsdGbatch/=0) then
       write(DECinfo%output,*) 'Gamma batch size was set manually, use that value instead!'
       GammaOpt=DECinfo%ccsdGbatch
    end if
  
    ! The optimal gamma batch size is GammaOpt.
    ! We now find the maximum possible gamma batch size smaller than or equal to GammaOpt
    ! and store this number in gammadim.
    call determine_MaxOrbitals(DECinfo%output,mylsitem%setting,GammaOpt,gammadim,'R')
  
  
    ! Largest possible alpha batch size
    ! =================================
    AlphaLoop: do alpha = MaxAObatch,MinAOBatch,-1
  
       call get_max_arraysizes_for_ccsdpt_integrals(alpha,gammadim,nbasis,nocc,nvirt,&
            & size1,size2,size3,MemoryNeeded)
  
       if(MemoryNeeded < MemoryAvailable .or. (alpha==minAObatch) ) then
          AlphaOpt = alpha
          exit
       end if
  
    end do AlphaLoop
  
    ! If alpha batch size was set manually we use that value instead
    if(DECinfo%ccsdAbatch/=0) then
       write(DECinfo%output,*) 'Alpha batch size was set manually, use that value instead!'
       AlphaOpt=DECinfo%ccsdAbatch
    end if
  
    ! The optimal alpha batch size is AlphaOpt.
    ! We now find the maximum possible alpha batch size smaller than or equal to AlphaOpt
    ! and store this number in alphadim.
    call determine_MaxOrbitals(DECinfo%output,mylsitem%setting,AlphaOpt,alphadim,'R')
  
  
    ! Print out and sanity check
    ! ==========================
 
#ifdef VAR_MPI

    if (infpar%lg_mynum .ne. infpar%master) goto 666

#endif

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '======================================================================='
    write(DECinfo%output,*) '                     CCSD(T) INTEGRALS: MEMORY SUMMARY                 '
    write(DECinfo%output,*) '======================================================================='
    write(DECinfo%output,*)
    write(DECinfo%output,*) 'To be on the safe side we use only 85% of the estimated available memory'
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a,g10.3)') '85% of available memory (GB)            =', MemoryAvailable
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a,i8)')    'Number of atomic basis functions        =', nbasis
    write(DECinfo%output,'(1X,a,i8)')    'Number of occupied orbitals             =', nocc
    write(DECinfo%output,'(1X,a,i8)')    'Number of virtual  orbitals             =', nvirt
    write(DECinfo%output,'(1X,a,i8)')    'Maximum alpha batch dimension           =', alphadim
    write(DECinfo%output,'(1X,a,i8)')    'Maximum gamma batch dimension           =', gammadim
    write(DECinfo%output,'(1X,a,g14.3)') 'Size of tmp array 1                     =', size1*realk*1.0E-9
    write(DECinfo%output,'(1X,a,g14.3)') 'Size of tmp array 2                     =', size2*realk*1.0E-9
    write(DECinfo%output,'(1X,a,g14.3)') 'Size of tmp array 3                     =', size3*realk*1.0E-9
    write(DECinfo%output,*)
  
#ifdef VAR_MPI

666 continue

#endif

    ! Sanity check
    call get_max_arraysizes_for_ccsdpt_integrals(alphadim,gammadim,nbasis,nocc,nvirt,&
         & size1,size2,size3,MemoryNeeded)  
    if(MemoryNeeded > MemoryAvailable) then
       write(DECinfo%output,*) 'Requested/available memory: ', MemoryNeeded, MemoryAvailable
       call lsquit('CCSD(T) integrals: Insufficient memory!',-1)
    end if


  end subroutine get_optimal_batch_sizes_ccsdpt_integrals



  !> \brief Get sizes of temporary arrays used in CCSD(T) integral routine (get_CCSDpT_integrals)
  !> with the chosen AO batch sizes.
  !> NOTE: If get_CCSDpT_integrals is modified, this routine must be changed accordingly!
  !> \author Kasper Kristensen & Janus Eriksen
  !> \date September 2011, rev. October 2012
  subroutine get_max_arraysizes_for_ccsdpt_integrals(alphadim,gammadim,nbasis,nocc,nvirt,&
                     & size1,size2,size3,mem)
    implicit none
    !> Max size for AO alpha batch
    integer,intent(in) :: alphadim
    !> Max size for AO gamma batch
    integer,intent(in) :: gammadim
    !> Number of AO basis functions
    integer,intent(in) :: nbasis
    !> Number of occupied (AOS) orbitals
    integer,intent(in) :: nocc
    !> Number of virt (AOS) orbitals
    integer,intent(in) :: nvirt
    !> Dimension of temporary array 1
    integer(kind=long),intent(inout) :: size1
    !> Dimension of temporary array 2
    integer(kind=long),intent(inout) :: size2
    !> Dimension of temporary array 3
    integer(kind=long),intent(inout) :: size3
    !> Tot size of temporary arrays (in GB)
    real(realk), intent(inout) :: mem
    real(realk) :: GB
    integer(kind=long) :: tmpI
    GB = 1.000E-9_realk ! 1 GB
    ! Array sizes needed in get_CCSDpT_integrals are checked and the largest one is found
  
    ! Tmp array 1 (five candidates)
    size1 = i8*alphadim*gammadim*nbasis*nbasis
    tmpI = i8*nvirt**2*gammadim*alphadim
    size1 = max(size1,tmpI)
    tmpI = i8*nvirt*nocc*gammadim*alphadim
    size1 = max(size1,tmpI)
    tmpI = i8*nvirt*nocc**2*alphadim
    size1 = max(size1,tmpI)
    tmpI = i8*nvirt**3
    size1 = max(size1,tmpI)
  
    ! tmp array 2 (three candidates)
    size2 = i8*alphadim*gammadim*nbasis*nvirt
    tmpI = alphadim*gammadim*nvirt*nocc
    size2 = max(size2,tmpI)
    tmpI = i8*nvirt**3
    size2 = max(size2,tmpI)
  
    ! Tmp array3 (two candidates)
    size3 = i8*alphadim*gammadim*nvirt**2
    tmpI = i8*alphadim*nvirt**3
    size3 = max(size3,tmpI)
  
    ! Size = size1+size2+size3,  convert to GB
    mem = realk*GB*(size1+size2+size3)


  end subroutine get_max_arraysizes_for_ccsdpt_integrals
!endif mod_unreleased
#endif

  subroutine dummy_ccsdpt_routine()

  end subroutine dummy_ccsdpt_routine

end module ccsdpt_module

#ifdef MOD_UNRELEASED

  !> \brief slaves enter here from lsmpi_slave (or dec_lsmpi_slave) and need to get to work 
  !> \author Janus Juul Eriksen
  !> \date x-mas 2012
#ifdef VAR_MPI

  subroutine ccsdpt_slave()

  use infpar_module
  use lsmpi_type
  use decmpi_module

  use precision
  use dec_typedef_module
  use memory_handling
  use lstiming, only: lstimer
  use typedeftype, only: Lsitem,lssetting

  ! DEC DEPENDENCIES (within deccc directory)  
  ! *****************************************
  use array2_simple_operations, only: array2_init_plain,array2_free 
  use array4_simple_operations, only: array4_init_standard,array4_free
  use atomic_fragment_operations
  use ccsdpt_module, only: ccsdpt_driver

    implicit none
    integer :: nocc, nvirt,nbasis
    real(realk), pointer :: ppfock(:,:), qqfock(:,:), Co(:,:), Cv(:,:)
    type(array2) :: ccsdpt_t1
    type(array4) :: ccsd_t2, ccsdpt_t2
    type(lsitem) :: mylsitem

    ! call ccsd(t) data routine in order to receive data from master
    call mpi_communicate_ccsdpt_calcdata(nocc,nvirt,nbasis,ppfock,qqfock,Co,Cv,ccsd_t2%val,mylsitem)

    ! init and receive ppfock
    call mem_alloc(ppfock,nocc,nocc)
    call ls_mpibcast(ppfock,nocc,nocc,infpar%master,infpar%lg_comm)

    ! init and receive qqfock
    call mem_alloc(qqfock,nvirt,nvirt)
    call ls_mpibcast(qqfock,nvirt,nvirt,infpar%master,infpar%lg_comm)

    ! init and receive Co
    call mem_alloc(Co,nbasis,nocc)
    call ls_mpibcast(Co,nbasis,nocc,infpar%master,infpar%lg_comm)

    ! init and receive Cv
    call mem_alloc(Cv,nbasis,nvirt)
    call ls_mpibcast(Cv,nbasis,nvirt,infpar%master,infpar%lg_comm)

    ! init and receive ccsd_doubles array4 structure
    ccsd_t2 = array4_init([nvirt,nocc,nvirt,nocc])
    call ls_mpibcast(ccsd_t2%val,nvirt,nocc,nvirt,nocc,infpar%master,infpar%lg_comm)

    ! init ccsd(t) singles and ccsd(t) doubles
    ccsdpt_t1 = array2_init_plain([nvirt,nocc])
    ccsdpt_t2 = array4_init_standard([nvirt,nvirt,nocc,nocc])

    ! now enter the ccsd(t) driver routine
    call ccsdpt_driver(nocc,nvirt,nbasis,ppfock,qqfock,Co,Cv,mylsitem,ccsd_t2,&
                         & ccsdpt_t1,ccsdpt_t2)

    ! now, release all amplitude arrays, both ccsd and ccsd(t)
    call array2_free(ccsdpt_t1)
    call array4_free(ccsd_t2)
    call array4_free(ccsdpt_t2)

    ! finally, release fragment or full molecule quantities
    call ls_free(mylsitem)
    call mem_dealloc(ppfock)
    call mem_dealloc(qqfock)
    call mem_dealloc(Co)
    call mem_dealloc(Cv)

  end subroutine ccsdpt_slave

#endif
!endif mod_unreleased
#endif
