!dalton_copyright_start
!
!
!dalton_copyright_end

module lucita_ci_task_interface

  use lucita_integral_density_interface

#ifdef VAR_MPI
#ifdef USE_MPI_MOD_F90
  use mpi
  implicit none
#else
  implicit none
#include "mpif.h"
#endif
#else
  implicit none
#endif

  public CI_task_list_interface
  public create_CI_task_list

  private

#ifdef VAR_MPI
  integer(kind=MPI_INTEGER_KIND) :: my_MPI_REAL8 = MPI_REAL8
  integer(kind=MPI_INTEGER_KIND) :: ierr_mpi, len_mpi, root_comm=0
#endif

contains

!**********************************************************************

  subroutine CI_task_list_interface(ci_task_list,                    &
                                    cref,                            &
                                    hc,                              &
                                    resolution_mat,                  &
                                    int1_or_rho1,                    &
                                    int2_or_rho2,                    & 
                                    block_list,                      &
                                    par_dist_block_list,             &
                                    proclist,                        &
                                    grouplist,                       &
                                    fh_array,                        &
                                    rcctos,                          &
                                    nbatch,                          &
                                    nblock,                          &
                                    print_lvl)
!-------------------------------------------------------------------------------
!
! purpose: interface routine to all individual CI tasks defined in a proper order 
!          in the list "ci_task_list".
!
!-------------------------------------------------------------------------------
      use lucita_mcscf_ci_cfg
      character (len=12), intent(in)   :: ci_task_list(*)
!------------- optional input depending on MCSCF/CI run ------------------------
      real(8),           intent(inout) :: cref(*)
      real(8),           intent(inout) :: hc(*)
      real(8),           intent(inout) :: resolution_mat(*)
      real(8),           intent(inout) :: int1_or_rho1(*)
      real(8),           intent(inout) :: int2_or_rho2(*)
!------------- end of optional input depending on MCSCF/CI run -----------------
      integer,            intent(in)   :: nbatch
      integer,            intent(in)   :: nblock
      integer,            intent(in)   :: block_list(nblock)
      integer,            intent(in)   :: par_dist_block_list(nblock)
      integer,            intent(in)   :: proclist(*)
      integer,            intent(in)   :: grouplist(*)
      integer,            intent(in)   :: fh_array(*)
      integer,            intent(in)   :: rcctos(nblock)
      integer,            intent(in)   :: print_lvl
!-------------------------------------------------------------------------------
      integer                          :: ci_task_ticket
!-------------------------------------------------------------------------------

!     initialize ci_task_ticket counter
      ci_task_ticket = 0

      select case(ci_task_list(ci_task_ticket+1))

        case ('return CIdia', 'perform Davi', 'return sigma', 'return rotVC', 'return densM', 'ijkl resort ', &
              'fci dump    ')

!         fill integral indices pointers and possibly read in integrals
!         -------------------------------------------------------------
          call lucita_pointer_integral_driver(ci_task_list(ci_task_ticket+1),          &
                                              int1_or_rho1,                            &
                                              int2_or_rho2,                            &
                                              mcscf_ci_update_ijkl,                    &
                                              mcscf_orbital_trial_vector,              &
                                              mcscf_ci_trial_vector,                   &
                                              integrals_from_mcscf_env,                &
                                              print_lvl)

        case default ! 'report CIspc', 'report CIana'
          
!         print *, ' skipped integral pointer construction and vector block allocation'

      end select

!          
      do ! infinte loop over CI task tickets

        ci_task_ticket = ci_task_ticket + 1

        select case(ci_task_list(ci_task_ticket))
          

          case ('report CIspc')
            exit ! return because there is nothing else to do
          case ('return CIdia')
            call calculate_CI_diagonal(cref,print_lvl,block_list,par_dist_block_list)
          case ('perform Davi')
            call davidson_ci_driver(block_list,par_dist_block_list,proclist,                      &
                                    grouplist,fh_array,rcctos,nblock,print_lvl,                   &
                                    cref,hc,resolution_mat,docisrdft_mc2lu)
          case ('return sigma')
            call return_sigma_vector(print_lvl,nbatch,block_list,par_dist_block_list,             &
                                     rcctos,grouplist,proclist,                                   &
                                     cref,hc,resolution_mat,int1_or_rho1,int2_or_rho2)
          case ('return rotVC')
            call return_bvec_transformed_2_new_mo_basis(cref,hc,resolution_mat,int1_or_rho1,      &
                                                        par_dist_block_list,block_list,           &
                                                        rcctos,grouplist,proclist)
          case ('report CIana')
            call report_CI_vector_analysis(cref,print_lvl)
          case ('return densM')
            call return_Xp_density_matrix(print_lvl,nbatch,block_list,par_dist_block_list,        &
                                          rcctos,grouplist,proclist,                              &
                                          cref,hc,resolution_mat,int1_or_rho1,int2_or_rho2)
          case default ! standard loop exit criterion
            exit
        end select

      end do

  end subroutine CI_task_list_interface
!**********************************************************************

  subroutine create_CI_task_list(lucita_ci_run_id,                  &
                                 report_ci_analysis,                &
                                 return_density_matrices,           &
                                 ci_task_list,                      &
                                 max_ci_tasks)
!-------------------------------------------------------------------------------
!
! purpose: translate external (MCSCF) CI run ID's in LUCITA internal 
!          CI tasks. add additional CI tasks as requested by input
!          parameters.
!
!-------------------------------------------------------------------------------
      character (len=12), intent(in)    :: lucita_ci_run_id
      integer           , intent(in)    :: report_ci_analysis
      integer           , intent(in)    :: return_density_matrices
      integer           , intent(in)    :: max_ci_tasks
      character (len=12), intent(inout) :: ci_task_list(*)
!-------------------------------------------------------------------------------
      integer                           :: ci_task_ticket
!-------------------------------------------------------------------------------

!     initialize
      do ci_task_ticket = 1, max_ci_tasks
        ci_task_list(ci_task_ticket) = 'undefined   '
      end do

      ci_task_ticket = 1

      select case(lucita_ci_run_id)

        case('return CIdim') ! calculate # of determinants per symmetry irrep
          ci_task_list(ci_task_ticket) = 'report CIspc'

        case('return CIdia') ! calculate the diagonal part of the CI Hamiltonian matrix
          ci_task_list(ci_task_ticket) = 'return CIdia'

        case('sigma vec   ') ! calculate a sigma vector
          ci_task_list(ci_task_ticket) = 'return sigma'

        case('Xp-density m') ! calculate the 1-/2-particle density matrix
          ci_task_list(ci_task_ticket) = 'return densM'

        case('analyze Cvec') ! analyze CI vector
          ci_task_list(ci_task_ticket) = 'report CIana'

        case('rotate  Cvec') ! rotate CI vector
          ci_task_list(ci_task_ticket) = 'return rotVC'

        case('srdft   ci  ') ! srdft ci
          ci_task_list(ci_task_ticket) = 'return CIdia'
          ci_task_ticket               = ci_task_ticket + 1
          ci_task_list(ci_task_ticket) = 'perform Davi'

          if(report_ci_analysis > 0)then
            ci_task_ticket               = ci_task_ticket + 1
            ci_task_list(ci_task_ticket) = 'report CIana'
          end if

          if(return_density_matrices > 0)then
            ci_task_ticket               = ci_task_ticket + 1
            ci_task_list(ci_task_ticket) = 'return densM'
          end if

        case('standard ci ', 'initial ci  ') ! perform Davidson CI run

          ci_task_list(ci_task_ticket) = 'return CIdia'
          ci_task_ticket               = ci_task_ticket + 1
          ci_task_list(ci_task_ticket) = 'perform Davi'

          if(report_ci_analysis > 0)then
            ci_task_ticket               = ci_task_ticket + 1
            ci_task_list(ci_task_ticket) = 'report CIana'
          end if

          if(return_density_matrices > 0)then
            ci_task_ticket               = ci_task_ticket + 1
            ci_task_list(ci_task_ticket) = 'return densM'
          end if

        case('ijkl resort ') ! resort integrals to lucita format (prior to a large-scale CI)
          ci_task_list(ci_task_ticket) = 'ijkl resort '

        case('fci dump    ') ! resort integrals to lucita format (prior to a large-scale CI)
          ci_task_list(ci_task_ticket) = 'fci dump    '

        case default
         call quit('error in create_CI_task_list: undefined CI run id.')
      end select

  end subroutine create_CI_task_list
!**********************************************************************

  subroutine calculate_CI_diagonal(cref,                         &
                                   print_lvl,                    &
                                   block_list,                   &
                                   par_dist_block_list)

  use lucita_energy_types
  use file_type_module
  use ttss_block_module

#ifdef VAR_MPI
  use par_mcci_io
#endif

#include "parluci.h"
#include "mxpdim.inc"
#include "crun.inc"
#include "oper.inc"
#include "clunit.inc"
!-------------------------------------------------------------------------------
      real(8), intent(inout)   :: cref(*)
      integer, intent(in)      :: print_lvl
      integer, intent(in)      :: block_list(num_blocks)
      integer, intent(in)      :: par_dist_block_list(num_blocks)
!-------------------------------------------------------------------------------
#ifdef VAR_MPI
#include "maxorb.h"
#include "infpar.h"
      integer(MPI_INTEGER_KIND):: my_STATUS(MPI_STATUS_SIZE)
      integer(MPI_OFFSET_KIND) :: offset
      integer                  :: block_length
      integer                  :: isblk
      integer                  :: iproc
      integer                  :: my_MPI_COMM_WORLD = MPI_COMM_WORLD
#endif
!-------------------------------------------------------------------------------

!     construct the diagonal part of the Hamiltonian
      if(idiag == 2)  ludia = lusc1
      if(icistr >= 2) call rewino(ludia)

      call gasdiat(cref,ludia,ecore_orig-ecore,icistr,i12,                         &
                   ttss_info%ttss_block_type,num_blocks,ttss_info%ttss_vec_split,  &
                   par_dist_block_list)

      if(nocsf == 1 .and. icistr == 1)then
        call rewino(ludia)
        call todsc_luci(cref,l_combi,l_combi,ludia)
      end if

      if(luci_nmproc > 1)then
#ifdef VAR_MPI
        rewind ludia

        call isetvc(file_info%iluxlist(1,file_info%current_file_nr_diag),1,         &
                    file_info%max_list_length)

!       collect H diagonal 
        call mcci_cp_vcd_mpi_2_seq_io_interface(cref,ludia,                         &
                                                file_info%fh_lu(file_info%          &
                                                current_file_nr_diag),              &
                                                file_info%file_offsets(             &
                                                file_info%current_file_nr_diag),    &
                                                file_info%iluxlist(1,               &
                                                file_info%current_file_nr_diag),    &
                                                par_dist_block_list,                &
                                                block_list,                         &
                                                my_MPI_COMM_WORLD,                  &
                                                num_blocks,1,1,2)
#endif
      end if

#ifdef LUCI_DEBUG
!     debug printing (parallel run)
      if(print_lvl .ge. 10)then
#ifdef VAR_MPI
        if(luci_nmproc > 1)then
          offset = file_info%file_offsets(file_info%current_file_nr_diag)
          do isblk = 1, num_blocks
            call get_block_proc(par_dist_block_list,isblk,iproc)
            if(luci_myproc == iproc)then
              call get_block_length(block_list,isblk,block_length)
              call mpi_file_read_at(file_info%fh_lu(file_info%current_file_nr_diag),&
                                    offset,cref,block_length,my_MPI_REAL8,my_STATUS,ierr_mpi)
              write(luwrt,*) ' # of diagonal elements ==> ',block_length
              call wrtmatmn(cref,1,block_length,1,block_length,luwrt)
              offset = offset + block_length
            end if
          end do
        end if
#endif
      end if
#endif

  end subroutine calculate_CI_diagonal
!**********************************************************************

  subroutine report_CI_vector_analysis(cref,print_lvl)

  use file_type_module
  use ttss_block_module

#include "parluci.h"
#include "mxpdim.inc"
#include "crun.inc"
#include "cstate.inc"
#include "clunit.inc"
!-------------------------------------------------------------------------------
      real(8), intent(inout)   :: cref(*)
      integer, intent(in)      :: print_lvl
!-------------------------------------------------------------------------------
      integer                  :: eigen_state_id
      integer                  :: imzero, iampack
      integer                  :: luc_vector_file
      integer                  :: lusc_vector_file
!-------------------------------------------------------------------------------

      if(luci_myproc == luci_master)then

        luc_vector_file = luc
!       MCSCF - CI task run: file handle set in MCSCF/CI interface routine
        if(file_info%current_file_fh_seqf(1) > 0) luc_vector_file = file_info%current_file_fh_seqf(1)

        call rewino(luc_vector_file)
        
        do eigen_state_id = 1, nroot

          if(icistr == 1)then
            call frmdsc_luci(cref,l_combi,l_combi,luc_vector_file,imzero,iampack)
          else
            lusc_vector_file = luc_vector_file
            if(nroot > 1)then
              call rewino(lusc1)
              call copvcd(luc_vector_file,lusc1,cref,0,-1)
              lusc_vector_file = lusc1
            end if
          end if
!
          write(luwrt,'(/a)')    '   ******************************************************'
          write(luwrt,'(a,i3)')  '   Analysis of leading coefficients for eigen state = ',eigen_state_id
          write(luwrt,'(a)')     '   ******************************************************'


          call rewino(lusc_vector_file)
          call gasana(l0block,num_blocks,ttss_info%ttss_vec_split,ttss_info%ttss_block_type,lusc_vector_file,icistr)
        end do
      end if

  end subroutine report_CI_vector_analysis
!**********************************************************************

  subroutine return_Xp_density_matrix(print_lvl,                    &      
                                      nbatch,                       &
                                      block_list,                   &
                                      par_dist_block_list,          &
                                      rcctos,                       &
                                      grouplist,                    &
                                      proclist,                     &
                                      cref,                         &
                                      hc,                           &
                                      resolution_mat,               &
                                      int1_or_rho1,                 &
                                      int2_or_rho2)
  use file_type_module
  use lucita_cfg
  use lucita_mcscf_ci_cfg
  use lucita_mcscf_srdftci_cfg
#ifdef VAR_MPI
  use par_mcci_io
#endif


#include "priunit.h"
#include "parluci.h"
#include "mxpdim.inc"
#include "crun.inc"
#include "cstate.inc"
#include "cprnt.inc"
#include "clunit.inc"
#include "oper.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "cicisp.inc"
#include "cands.inc"
#include "glbbas.inc"
!-------------------------------------------------------------------------------
      integer, intent(in)              :: print_lvl
      integer, intent(in)              :: nbatch
      integer, intent(in)              :: block_list(num_blocks)
      integer, intent(in)              :: par_dist_block_list(num_blocks)
      integer, intent(in)              :: rcctos(num_blocks)
      integer, intent(in)              :: grouplist(luci_nmproc)
      integer, intent(in)              :: proclist(luci_nmproc)
!-------------------------------------------------------------------------------
      real(8),           intent(inout) :: cref(*)
      real(8),           intent(inout) :: hc(*)
      real(8),           intent(inout) :: resolution_mat(*)
      real(8),           intent(inout) :: int1_or_rho1(*)
      real(8),           intent(inout) :: int2_or_rho2(*)
!-------------------------------------------------------------------------------
      real(8)                          :: work
#include "wrkspc.inc"
      real(8)                          :: test_energy
      real(8)                          :: exps2
      real(8)                          :: cv_dummy, hcv_dummy
      integer, parameter               :: isigden = 2 ! density matrix switch for sigden_ci
      integer                          :: eigen_state_id
      integer                          :: imzero, iampack
      integer                          :: lusc_vector_file
      integer                          :: luhc_vector_file
      integer                          :: luc_internal
      integer                          :: luhc_internal
      integer                          :: twopart_densdim
      integer                          :: rhotype
      integer                          :: rhotype_spin1
      integer                          :: idum
      integer(8)                       :: kdum
      integer(8)                       :: k_scratchsbatch = 1 ! to avoid compiler warnings
      integer(8)                       :: k_dens2_scratch
      integer(8)                       :: k_scratch1
      integer(8)                       :: k_scratch2
      integer(8)                       :: k_scratch3
      integer(8)                       :: k_scratch4
#ifdef VAR_MPI
      integer, allocatable             :: blocks_per_batch(:)
      integer, allocatable             :: batch_length(:)
      integer, allocatable             :: block_offset_batch(:)
      integer, allocatable             :: block_info_batch(:,:)
      integer                          :: iatp, ibtp
      integer                          :: nbatch_par, nblock_par
      integer                          :: iloop, nbatch_par_max
      integer(kind=MPI_INTEGER_KIND)   :: mynew_comm_mpi, my_luci_master
      integer(kind=MPI_OFFSET_KIND)    :: my_lu4_off_tmp
      integer                          :: my_MPI_COMM_WORLD = MPI_COMM_WORLD
#endif
!-------------------------------------------------------------------------------

!     set level of particle-density matrix calculation
      i12 = idensi

!     initialize s**2 expectation value
      exps2 = 0.0d0

!     set rhotype
      rhotype = 1
!     remember: lucita-internal routines for natorbs expects an unpacked density matrix
      if(lucita_cfg_natural_orb_occ_nr) rhotype = 0
      if(lucita_cfg_transition_densm)   rhotype = 3

!     MCSCF - CI task run: file handle set in MCSCF/CI interface routine
      luc_internal  = luc
      luhc_internal = luhc
      if(file_info%current_file_fh_seqf(1) > 0) luc_internal  = file_info%current_file_fh_seqf(1)
      if(file_info%current_file_fh_seqf(2) > 0) luhc_internal = file_info%current_file_fh_seqf(2)

#ifdef VAR_MPI
      my_luci_master = luci_master
      if(luci_nmproc > 1)then

        if(icsm /= issm) & 
        call quit('*** return_Xp_density_matrix: property/response density matrix not parallelized yet. ***')

!       sequential --> MPI I/O
!       ----------------------
        allocate(blocks_per_batch(mxntts))
        allocate(batch_length(mxntts))
        allocate(block_offset_batch(mxntts))
        allocate(block_info_batch(8,mxntts))
        blocks_per_batch   = 0
        batch_length       = 0
        block_offset_batch = 0
        block_info_batch   = 0
        iatp = 1
        ibtp = 2
        call Z_BLKFO_partitioning_parallel(icspc,icsm,iatp,ibtp,             &
                                           blocks_per_batch,batch_length,    &
                                           block_offset_batch,               &
                                           block_info_batch,nbatch_par,      &
                                           nblock_par,par_dist_block_list)

        if(luci_myproc == luci_master .and. rhotype == 1)then
          idum = 0
          call memman(idum,idum,'MARK  ',idum,'Xpden0')
          nbatch_par_max = 0
          do iloop = 1, nbatch_par
            nbatch_par_max = max(nbatch_par_max,batch_length(iloop))
          end do
!         write(lupri,*) 'allocate scratch space for sigma vec: ',nbatch_par_max
          call memman(k_scratchsbatch,nbatch_par_max,'ADDL  ',2,'SBATCH')
          call dzero(work(k_scratchsbatch),nbatch_par_max)
        end if
         

        if(.not.file_info%file_type_mc)then ! not within MCSCF

!         hardwired to ILU1 + ILU2 for CI
          file_info%current_file_nr_active1 = 2
          file_info%current_file_nr_active2 = 2
          call rewino(luc_internal)

          call izero(file_info%iluxlist(1,file_info%current_file_nr_active1),  &
                     file_info%max_list_length)
          call izero(file_info%iluxlist(1,file_info%current_file_nr_active2),  &
                     file_info%max_list_length)

!         step 1: the rhs vector
          call mcci_cp_vcd_mpi_2_seq_io_interface(cref,                                  &
                                                  luc_internal,                          &
                                                  file_info%fh_lu(file_info%             &
                                                  current_file_nr_active1),              &
                                                  file_info%file_offsets(                &
                                                  file_info%current_file_nr_active1),    &
                                                  file_info%iluxlist(1,                  &
                                                  file_info%current_file_nr_active1),    &
                                                  par_dist_block_list,                   &
                                                  block_list,                            &
                                                  my_MPI_COMM_WORLD,                     &
                                                  NUM_BLOCKS,NROOT,1,1)
        end if

        if(lucita_cfg_transition_densm)then

          if(.not.file_info%file_type_mc)then ! not within MCSCF

!           step 2: the lhs vector
            file_info%current_file_nr_active2 = 3
            call rewino(luhc_internal)

            call mcci_cp_vcd_mpi_2_seq_io_interface(hc,                                  &
                                                    luhc_internal,                       &
                                                    file_info%fh_lu(file_info%           &
                                                    current_file_nr_active2),            &
                                                    file_info%file_offsets(              &
                                                    file_info%current_file_nr_active2),  &
                                                    file_info%iluxlist(1,                &
                                                    file_info%current_file_nr_active2),  &
                                                    par_dist_block_list,                 &
                                                    block_list,                          &
                                                    my_MPI_COMM_WORLD,                   &
                                                    NUM_BLOCKS,NROOT,1,1)
          end if

        else

          call icopy(file_info%max_list_length,file_info%iluxlist(1,file_info%current_file_nr_active1),1, &
                                               file_info%iluxlist(1,file_info%current_file_nr_active2),1)
        end if

        luhc_vector_file = file_info%fh_lu(file_info%current_file_nr_active2)
!       save and temporarily re-set the file offset ("ilu4" is used inside sigden_ci)
        my_lu4_off_tmp   = my_lu4_off
        my_lu4_off       = file_info%file_offsets(file_info%current_file_nr_active2)

        lusc_vector_file = file_info%fh_lu(file_info%current_file_nr_bvec)
      end if
#endif

      call rewino(luc_internal)
      call rewino(luhc_internal)

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      iprden = 100
      write(luwrt,'(/a)') ' IPRDEN raised explicitly in return_Xp_density_matrix '
#endif

      do eigen_state_id = 1, nroot

        write(luwrt,'(/a)')    '   **************************************************************'
        write(luwrt,'(a,i3)')  '   calculating 1-/2-particle density matrix for eigen state = ', eigen_state_id
        write(luwrt,'(a)')     '   **************************************************************'

        if(luci_nmproc == 1)then
          if(icistr == 1)then
            call frmdsc_luci(cref,l_combi,l_combi,luc_internal,imzero,iampack)
            if(lucita_cfg_transition_densm)then
!             transition density
              call frmdsc_luci(hc,l_combi,l_combi,luhc_internal,imzero,iampack)
            else
              call dcopy(l_combi,cref,1,hc,1)
            end if
          else
            if(lucita_cfg_transition_densm)then
!             transition density
              lusc_vector_file = luc_internal
              luhc_vector_file = luhc_internal
            else
              call rewino(lusc1)
              call copvcd(luc_internal,lusc1,cref,0,-1)

              call copvcd(lusc1,lusc2,cref,1,-1)
              call rewino(lusc1)
              call rewino(lusc2)
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
              CALL REWINE(lusc1,-1)
              WRITE(LUWRT,*) '  final solution vector for root ==> ',eigen_state_id
              CALL WRTVCD(cref,LUsC1,0,-1)
              CALL REWINE(LUsC1,-1)
#undef LUCI_DEBUG
#endif
              lusc_vector_file = lusc1
              luhc_vector_file = lusc2
            end if
          end if
        else
#ifdef VAR_MPI
!         MPI I/O --> MPI I/O node-master collection file
!         -----------------------------------------------
          mynew_comm_mpi = mynew_comm
          call mpi_barrier(mynew_comm_mpi,ierr_mpi) ! probably not needed if collection of densities is in action
          call izero(file_info%ilublist,file_info%max_list_length_bvec)

          call mcci_cp_vcd_batch(file_info%fh_lu(file_info%current_file_nr_active1),       &    
                                 file_info%fh_lu(file_info%current_file_nr_bvec),          &
                                 cref,nbatch_par,blocks_per_batch,                         &
                                 batch_length,block_offset_batch,block_info_batch,         &
                                 block_list,                                               &
                                 file_info%file_offsets(file_info%current_file_nr_active1),&
                                 file_info%file_offsets(file_info%current_file_nr_bvec),   &
                                 file_info%iluxlist(1,file_info%current_file_nr_active1),  &
                                 file_info%ilublist,                                       &
                                 my_vec2_ioff,my_act_blk2,eigen_state_id-1)
!         offset in MPI I/O file for lhs-vector
          jvec_sf = eigen_state_id - 1
#endif
        end if

!
!       stefan - jan 2012: the following construct is a consequence of the physical memory
!                          identity of cref and hc on the master node in case of a parallel mcscf
!                          and a regular density matrix run
!
        if(luci_nmproc > 1 .and. luci_myproc == luci_master .and.  rhotype == 1)then
          call sigden_ci(cref,work(k_scratchsbatch),resolution_mat,lusc_vector_file,luhc_vector_file,&
                         cv_dummy,hcv_dummy,isigden,rhotype,exps2                                    &
#ifdef VAR_MPI
                        ,file_info%ilublist,file_info%iluxlist(1,file_info%current_file_nr_active2), &
                         block_list,par_dist_block_list,                                             &
                         rcctos,grouplist,proclist                                                   &
#endif
                        )
        else
          call sigden_ci(cref,hc,resolution_mat,lusc_vector_file,luhc_vector_file,                   &
                         cv_dummy,hcv_dummy,isigden,rhotype,exps2                                    &
#ifdef VAR_MPI
                        ,file_info%ilublist,file_info%iluxlist(1,file_info%current_file_nr_active2), &
                         block_list,par_dist_block_list,                                             &
                         rcctos,grouplist,proclist                                                   &
#endif
                        )
        end if


        if(luci_nmproc > 1)then
#ifdef VAR_MPI
!         collect density matrices
!         ------------------------
          if(luci_myproc == luci_master)then
            len_mpi = nacob**2
            call mpi_reduce(mpi_in_place,work(krho1),len_mpi,my_mpi_real8,           &
                            mpi_sum,my_luci_master,mpi_comm_world,ierr_mpi)
            if(ispnden > 0)then
              call mpi_reduce(mpi_in_place,work(ksrho1),len_mpi,my_mpi_real8,        &
                              mpi_sum,my_luci_master,mpi_comm_world,ierr_mpi)
              call mpi_reduce(mpi_in_place,work(ksrho1a),len_mpi,my_mpi_real8,       &
                              mpi_sum,my_luci_master,mpi_comm_world,ierr_mpi)
              call mpi_reduce(mpi_in_place,work(ksrho1b),len_mpi,my_mpi_real8,       &
                              mpi_sum,my_luci_master,mpi_comm_world,ierr_mpi)
            end if
            if(i12 > 1)then
              len_mpi = nacob**2*(nacob**2+1)/2
              call mpi_reduce(mpi_in_place,work(krho2),len_mpi,                      &
                              my_mpi_real8,mpi_sum,my_luci_master,mpi_comm_world,ierr_mpi)
              len_mpi = 1
              call mpi_reduce(mpi_in_place,exps2,len_mpi,                            &
                              my_mpi_real8,mpi_sum,my_luci_master,mpi_comm_world,ierr_mpi)
            end if
          else
            len_mpi = nacob**2
            call mpi_reduce(work(krho1),mpi_in_place,len_mpi,my_mpi_real8,           &
                            mpi_sum,my_luci_master,mpi_comm_world,ierr_mpi)
            if(ispnden > 0)then
              call mpi_reduce(work(ksrho1),mpi_in_place,len_mpi,my_mpi_real8,        &
                              mpi_sum,my_luci_master,mpi_comm_world,ierr_mpi)
              call mpi_reduce(work(ksrho1a),mpi_in_place,len_mpi,my_mpi_real8,       &
                              mpi_sum,my_luci_master,mpi_comm_world,ierr_mpi)
              call mpi_reduce(work(ksrho1b),mpi_in_place,len_mpi,my_mpi_real8,       &
                              mpi_sum,my_luci_master,mpi_comm_world,ierr_mpi)
            end if
            if(i12 > 1)then
              len_mpi = nacob**2*(nacob**2+1)/2
              call mpi_reduce(work(krho2),mpi_in_place,len_mpi,                      &
                              my_mpi_real8,mpi_sum,my_luci_master,mpi_comm_world,ierr_mpi)
              len_mpi = 1
              call mpi_reduce(exps2,mpi_in_place,len_mpi,                            &
                              my_mpi_real8,mpi_sum,my_luci_master,mpi_comm_world,ierr_mpi)
            end if
          end if
#endif
        end if ! luci_nmproc > 1

        if(i12 > 1)then
          write(luwrt,'(/a      )')  '   ------------------------------------------------'
          write(luwrt,'(a,1f10.3)')  '   expectation value of operator <S**2> =', exps2
          write(luwrt,'(a/      )')  '   ------------------------------------------------'
        end if

!       export 1-/2-particle density matrix to MCSCF format
!       ---------------------------------------------------
!
!#ifdef largeFCInotneeded
        idum = 0
        call memman(idum,idum,'MARK  ',idum,'Xpden1')
        if(i12  > 1) call memman(k_dens2_scratch,nacob**4,'ADDL  ',2,'PVfull')
        if(i12 == 1) call memman(k_dens2_scratch,       0,'ADDL  ',2,'PVfull')

!       set rhotype for spin-densities after reordering
        rhotype_spin1 = 1
        if(lucita_cfg_natural_orb_occ_nr) rhotype_spin1 = rhotype

!       activate for MCSCF/improved CI
        if(integrals_from_mcscf_env)then

          !> reoder first spin-densities in order to not overwrite the outgoing 1p density matrix in int1_or_rho1
          if(ispnden == 1)then
            call lucita_spinputdens_1p(work(kSRHO1a),work(krho2),int1_or_rho1,int2_or_rho2,          &
                                       work(k_dens2_scratch),i12,isigden,rhotype_spin1,eigen_state_id,1)
            call lucita_spinputdens_1p(work(kSRHO1b),work(krho2),int1_or_rho1,int2_or_rho2,          &
                                       work(k_dens2_scratch),i12,isigden,rhotype_spin1,eigen_state_id,2)
          end if
          call lucita_putdens_generic(work(krho1),work(krho2),int1_or_rho1,int2_or_rho2,             &
                                      work(k_dens2_scratch),i12,isigden,rhotype,eigen_state_id)

        else
          twopart_densdim = 0
          if(i12 > 1) twopart_densdim = (nacob*(nacob+1)/2)**2
          call memman(k_scratch1     ,nacob**2           ,'ADDL  ',2,'scrat1')
          call memman(k_scratch2     ,twopart_densdim    ,'ADDL  ',2,'scrat2')

          call lucita_putdens_generic(work(krho1),work(krho2),work(k_scratch1),work(k_scratch2),     &
                                      work(k_dens2_scratch),i12,isigden,rhotype,eigen_state_id)
          if(ispnden >= 1)then
            call lucita_spinputdens_1p(work(ksrho1a),work(krho2),work(k_scratch1),work(k_scratch2),  &
                                       work(k_dens2_scratch),  1,isigden,rhotype_spin1,eigen_state_id,1)
            call lucita_spinputdens_1p(work(ksrho1b),work(krho2),work(k_scratch1),work(k_scratch2),  &
                                       work(k_dens2_scratch),  1,isigden,rhotype_spin1,eigen_state_id,2)
          end if
        end if
!
        call memman(kdum ,idum,'FLUSM ',2,'Xpden1')
!#endif

!       natural orbital occupation numbers
!       ----------------------------------
        if(lucita_cfg_natural_orb_occ_nr)then
          idum = 0
          call memman(idum,idum,'MARK  ',idum,'Xpden2')
          call memman(k_scratch1,nacob**2         ,'ADDL  ',2,'scrat1')
          call memman(k_scratch2,nacob**2         ,'ADDL  ',2,'scrat2')
          call memman(k_scratch3,nacob            ,'ADDL  ',2,'scrat3')
          call memman(k_scratch4,nacob*(nacob+1)/2,'ADDL  ',2,'scrat4')
!         if (ispnden == 1 .and. lucita_cfg_is_spin_multiplett /= 1) then
          if (ispnden >= 1)then
            write(luwrt,'(/a)') ' Natural spin-orbital occupation numbers for alpha spin-orbitals'
            call lnatorb(work(ksrho1a),nsmob,ntoobs,nacobs,ninobs,                                 &
                         ireost,work(k_scratch1),work(k_scratch2),work(k_scratch3),nacob,          &
                         work(k_scratch4),1)
            write(luwrt,'(/a)') ' Natural spin-orbital occupation numbers for beta spin-orbitals'
            call lnatorb(work(ksrho1b),nsmob,ntoobs,nacobs,ninobs,                                 &
                         ireost,work(k_scratch1),work(k_scratch2),work(k_scratch3),nacob,          &
                         work(k_scratch4),1)
          end if

          write(luwrt,'(/a)') ' Natural orbital occupation numbers'
          call lnatorb(work(krho1),nsmob,ntoobs,nacobs,ninobs,                                     &
                       ireost,work(k_scratch1),work(k_scratch2),work(k_scratch3),nacob,            &
                       work(k_scratch4),1)

          call memman(kdum ,idum,'FLUSM ',2,'Xpden2')
        end if

        if(i12 > 1 .and. .not. srdft_ci_with_lucita) call en_from_dens(test_energy,i12)

      end do ! loop over eigen states
#ifdef LUCI_DEBUG
          write(luwrt,'(/a)') ' IPRDEN lowered explicitly in return_Xp_density_matrix '
          iprden = 1
#endif

#ifdef VAR_MPI
      if(luci_nmproc > 1)then
        if(luci_myproc == luci_master .and. rhotype == 1)then
          call memman(kdum ,idum,'FLUSM ',2,'Xpden0')
        end if
!       reset file off-set and free memory
        if(.not.lucita_cfg_transition_densm) my_lu4_off = my_lu4_off_tmp
        deallocate(block_info_batch)
        deallocate(block_offset_batch)
        deallocate(batch_length)
        deallocate(blocks_per_batch)
      end if
#endif

  end subroutine return_Xp_density_matrix
!**********************************************************************

  subroutine return_sigma_vector(print_lvl,                    &      
                                 nbatch,                       &
                                 block_list,                   &
                                 par_dist_block_list,          &
                                 rcctos,                       &
                                 grouplist,                    &
                                 proclist,                     &
                                 cref,                         &
                                 hc,                           &
                                 resolution_mat,               &
                                 int1_or_rho1,                 &
                                 int2_or_rho2)
  use file_type_module
  use lucita_cfg
  use lucita_mcscf_ci_cfg
#ifdef VAR_MPI
  use par_mcci_io
#endif


#include "parluci.h"
#include "mxpdim.inc"
#include "crun.inc"
#include "cstate.inc"
#include "clunit.inc"
#include "oper.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "cicisp.inc"
#include "cands.inc"
#include "glbbas.inc"
!-------------------------------------------------------------------------------
      integer, intent(in)              :: print_lvl
      integer, intent(in)              :: nbatch
      integer, intent(in)              :: block_list(num_blocks)
      integer, intent(in)              :: par_dist_block_list(num_blocks)
      integer, intent(in)              :: rcctos(num_blocks)
      integer, intent(in)              :: grouplist(luci_nmproc)
      integer, intent(in)              :: proclist(luci_nmproc)
!-------------------------------------------------------------------------------
      real(8),           intent(inout) :: cref(*)
      real(8),           intent(inout) :: hc(*)
      real(8),           intent(inout) :: resolution_mat(*)
      real(8),           intent(inout) :: int1_or_rho1(*)
      real(8),           intent(inout) :: int2_or_rho2(*)
!-------------------------------------------------------------------------------
      real(8)                          :: cv_dummy, hcv_dummy,exps2
      integer                          :: eigen_state_id
      integer                          :: imzero, iampack
      integer                          :: lusc_vector_file
      integer                          :: luhc_vector_file
      integer                          :: luc_internal
      integer                          :: luhc_internal
      integer, parameter               :: isigden = 1 ! sigma vector switch for sigden_ci
#ifdef VAR_MPI
      integer, allocatable             :: blocks_per_batch(:)
      integer, allocatable             :: batch_length(:)
      integer, allocatable             :: block_offset_batch(:)
      integer, allocatable             :: block_info_batch(:,:)
      integer                          :: iatp, ibtp
      integer                          :: nbatch_par, nblock_par
      integer(kind=MPI_INTEGER_KIND)   :: mynew_comm_mpi
      integer(kind=mpi_offset_kind)    :: my_lu4_off_tmp
      integer                          :: my_MPI_COMM_WORLD = MPI_COMM_WORLD
#endif
!-------------------------------------------------------------------------------

!     MCSCF - CI task run: file handle set in MCSCF/CI interface routine
      luc_internal  = luc
      luhc_internal = luhc
      if(file_info%current_file_fh_seqf(1) > 0) luc_internal  = file_info%current_file_fh_seqf(1)

#ifdef VAR_MPI
      if(luci_nmproc > 1)then

        if(icsm /= issm) & 
        call quit('*** return_sigma_vector: property/response E2b calculation not parallelized yet. ***')

!       sequential --> MPI I/O
!       ----------------------

        allocate(blocks_per_batch(mxntts))
        allocate(batch_length(mxntts))
        allocate(block_offset_batch(mxntts))
        allocate(block_info_batch(8,mxntts))
        blocks_per_batch   = 0
        batch_length       = 0
        block_offset_batch = 0
        block_info_batch   = 0
        iatp = 1
        ibtp = 2
        call Z_BLKFO_partitioning_parallel(icspc,icsm,iatp,ibtp,             &
                                           blocks_per_batch,batch_length,    &
                                           block_offset_batch,               &
                                           block_info_batch,nbatch_par,      &
                                           nblock_par,par_dist_block_list)

!       set HC vector file target
        file_info%current_file_nr_active2 = 3

        if(.not.file_info%file_type_mc)then ! not within MCSCF

!         hardwired to ILU1 for CI
          file_info%current_file_nr_active1 = 2

          call izero(file_info%iluxlist(1,file_info%current_file_nr_active1),  &
                     file_info%max_list_length)
          call izero(file_info%iluxlist(1,file_info%current_file_nr_active2),  &
                     file_info%max_list_length)

          call rewino(luc_internal)
!         step 1: the rhs vector
          call mcci_cp_vcd_mpi_2_seq_io_interface(cref,                                  &
                                                  luc_internal,                          &
                                                  file_info%fh_lu(file_info%             &
                                                  current_file_nr_active1),              &
                                                  file_info%file_offsets(                &
                                                  file_info%current_file_nr_active1),    &
                                                  file_info%iluxlist(1,                  &
                                                  file_info%current_file_nr_active1),    &
                                                  par_dist_block_list,                   &
                                                  block_list,                            &
                                                  my_MPI_COMM_WORLD,                     &
                                                  NUM_BLOCKS,NROOT,1,1)
        end if ! not within MCSCF
      end if
#endif

      call rewino(luhc_internal)
      call rewino(luc_internal)

      do eigen_state_id = 1, nroot

        write(luwrt,'(/a)')    '   **********************************************'
        write(luwrt,'(a,i3)')  '   calculating E2b lhs vector for b state = ', eigen_state_id
        write(luwrt,'(a)')     '   **********************************************'

        if(luci_nmproc == 1)then
          if(icistr == 1)then
            call frmdsc_luci(cref,l_combi,l_combi,luc_internal,imzero,iampack)
          else
            call rewino(lusc1)
            call copvcd(luc_internal,lusc1,cref,0,-1)
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
            call wrtvcd(cref,lusc1,-1,-1)
#endif
            call rewino(lusc1)
          end if
          lusc_vector_file = lusc1
          luhc_vector_file = luhc_internal
        else
#ifdef VAR_MPI
!         MPI I/O --> MPI I/O node-master collection file
!         -----------------------------------------------
          mynew_comm_mpi = mynew_comm
          call mpi_barrier(mynew_comm_mpi,ierr_mpi) 
          call izero(file_info%ilublist,file_info%max_list_length_bvec)

          call mcci_cp_vcd_batch(file_info%fh_lu(file_info%current_file_nr_active1),       &    
                                 file_info%fh_lu(file_info%current_file_nr_bvec),          &
                                 hc,nbatch_par,blocks_per_batch,                           &
                                 batch_length,block_offset_batch,block_info_batch,         &
                                 block_list,                                               &
                                 file_info%file_offsets(file_info%current_file_nr_active1),&
                                 file_info%file_offsets(file_info%current_file_nr_bvec),   &
                                 file_info%iluxlist(1,file_info%current_file_nr_active1),  &
                                 file_info%ilublist,                                       &
                                 my_vec2_ioff,my_act_blk2,eigen_state_id-1)

          lusc_vector_file = file_info%fh_lu(file_info%current_file_nr_bvec)
          luhc_vector_file = file_info%fh_lu(file_info%current_file_nr_active2)

!         save and temporarily re-set the file offset ("ilu4" is used inside sigden_ci)
          my_lu4_off_tmp   = my_lu4_off
          my_lu4_off       = file_info%file_offsets(file_info%current_file_nr_active2)

!         offset in MPI I/O file for lhs-vector
          jvec_sf = eigen_state_id - 1
#endif
        end if
!
        call sigden_ci(cref,hc,resolution_mat,lusc_vector_file,luhc_vector_file,                  &
                       cv_dummy,hcv_dummy,isigden,-1,exps2                                        &
#ifdef VAR_MPI
                      ,file_info%ilublist,file_info%iluxlist(1,file_info%current_file_nr_active2),&
                       block_list,par_dist_block_list,                                            &
                       rcctos,grouplist,proclist                                                  &
#endif
                        )


      end do ! loop over eigen states

      if(luci_nmproc > 1)then
#ifdef VAR_MPI
        if(.not.file_info%file_type_mc)then ! not within MCSCF
!         collect e2b vector(s)
!         --------------------
          call rewino(luhc_internal)
  
!         the lhs vector
          call mcci_cp_vcd_mpi_2_seq_io_interface(hc,luhc_internal,luhc_vector_file,   &
                                                  file_info%file_offsets(              &
                                                  file_info%current_file_nr_active2),  &
                                                  file_info%iluxlist(1,                &
                                                  file_info%current_file_nr_active2),  &
                                                  par_dist_block_list,                 &
                                                  block_list,                          &
                                                  my_MPI_COMM_WORLD,                   &
                                                  num_blocks,nroot,1,2)

        end if ! not within MCSCF

        my_lu4_off = my_lu4_off_tmp

        deallocate(block_info_batch)
        deallocate(block_offset_batch)
        deallocate(batch_length)
        deallocate(blocks_per_batch)
#endif
      end if ! luci_nmproc > 1

  end subroutine return_sigma_vector
!**********************************************************************

  subroutine return_bvec_transformed_2_new_mo_basis(vec1,vec2,c2,mo2mo_mat, &
                                                    par_dist_block_list,    &
                                                    block_list,rcctos,      &
                                                    grouplist,proclist)
!
! purpose: master routine for transforming a CI vector to new orbital basis
!
! Jeppe Olsen, January 98
! re-factored and parallelized by Stefan Knecht, November 2011
!
!-------------------------------------------------------------------------------
 use file_type_module
!-------------------------------------------------------------------------------
#include "mxpdim.inc"
#include "lucinp.inc"
#include "clunit.inc"
#include "cgas.inc"
#include "csm.inc"
#include "crun.inc"
#include "cstate.inc"
#include "glbbas.inc"
#include "orbinp.inc"
#include "cicisp.inc"
#include "cands.inc"
#include "parluci.h"
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),&
                    ADSXA(MXPOBS,2*MXPOBS),                     &
                    SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
      real(8), intent(inout) :: vec1(*)
      real(8), intent(inout) :: vec2(*)
      real(8), intent(inout) :: c2(*)
      real(8), intent(inout) :: mo2mo_mat(*)
      integer, intent(in)    :: par_dist_block_list(*)
      integer, intent(in)    :: block_list(*)
      integer, intent(in)    :: proclist(*)
      integer, intent(in)    :: grouplist(*)
      integer, intent(in)    :: rcctos(*)
!-------------------------------------------------------------------------------
      real(8)                :: work
#include "wrkspc.inc"
      integer                :: my_in_fh, my_out_fh, my_sc1_fh, my_sc2_fh
      integer                :: my_BVC_fh
      integer                :: len_ilu1
      integer                :: len_ilu2
      integer                :: len_iluc
      integer                :: iatp, ibtp
      integer                :: nbatch, nblock
      integer                :: lu_ref, lu_refout
      integer                :: luc_internal, luhc_internal
      integer(kind=8)        :: my_in_off
      integer(kind=8)        :: my_out_off
      integer(kind=8)        :: my_scr_off
      integer(kind=8)        :: my_BVC_off

      integer, allocatable   :: blocks_per_batch(:)
      integer, allocatable   :: batch_length(:)
      integer, allocatable   :: block_offset_batch(:)
      integer, allocatable   :: block_info_batch(:,:)
      integer, allocatable   :: blocktype(:)
      integer(8)             :: my_lu4_off_tmp
!-------------------------------------------------------------------------------

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      WRITE(luwrt,*) ' Welcome to return_bvec_transformed_2_new_mo_basis'
      WRITE(luwrt,*) ' ================================================='
      call flshfo(luwrt)
#endif

!     MCSCF - CI task run: file handle set in MCSCF/CI interface routine
      luc_internal  = luc
      luhc_internal = luhc
      if(file_info%current_file_fh_seqf(1) > 0) luc_internal  = file_info%current_file_fh_seqf(1)

      len_ilu1   = 0
      len_ilu2   = 0
      len_iluc   = 0
      my_in_off  = 0
      my_out_off = 0
      my_scr_off = 0
      my_BVC_off = 0
      my_BVC_fh  = 0
#ifdef VAR_MPI
      file_info%current_file_nr_active2 = 3

      my_in_fh   = file_info%fh_lu(file_info%current_file_nr_active1)
      my_out_fh  = file_info%fh_lu(file_info%current_file_nr_active1)
      my_sc1_fh  = file_info%fh_lu(file_info%current_file_nr_active2)
      my_sc2_fh  = file_info%fh_lu(file_info%current_file_nr_active2)
      my_BVC_fh  = file_info%fh_lu(file_info%current_file_nr_bvec)

      my_in_off  = file_info%file_offsets(            &
                   file_info%current_file_nr_active1)
      my_out_off = file_info%file_offsets(            &
                   file_info%current_file_nr_active1)
      my_scr_off = file_info%file_offsets(            &
                   file_info%current_file_nr_active2)
      my_BVC_off = file_info%file_offsets(            &
                   file_info%current_file_nr_bvec)

      len_ilu1   = file_info%max_list_length
      len_ilu2   = file_info%max_list_length
      len_iluc   = file_info%max_list_length_bvec

!     save and temporarily re-set the file offset ("ilu4" is used inside sigden_ci)
      my_lu4_off_tmp   = my_lu4_off
      my_lu4_off       = file_info%file_offsets(file_info%current_file_nr_active2)

!     update co-workers with single-orbital transformation matrix
      len_mpi = nacob**2
      call mpi_bcast(mo2mo_mat,len_mpi,my_mpi_real8,root_comm,mpi_comm_world,ierr_mpi)

#else
      my_in_fh   = lusc1
      my_out_fh  = luhc_internal
      my_sc1_fh  = lusc2
      my_sc2_fh  = lusc3

      call copvcd(luc_internal,my_in_fh,vec1,1,-1)
      call rewine(my_in_fh,-1)
      CALL rewine(my_out_fh,-1)
#endif
      lu_ref     = luc_internal
      lu_refout  = luhc_internal

!     set up block and batch structure of vector

      allocate(blocks_per_batch(mxntts))
      allocate(batch_length(mxntts))
      allocate(block_offset_batch(mxntts))
      allocate(block_info_batch(8,mxntts))
      allocate(blocktype(nsmst))

      blocks_per_batch   = 0
      batch_length       = 0
      block_offset_batch = 0
      block_info_batch   = 0
      blocktype          = 0

      IATP     = 1
      IBTP     = 2
      CALL Z_BLKFO_partitioning_parallel(ICSPC,ICSM,iatp,ibtp,   &
                                         blocks_per_batch,       &
                                         batch_length,           &
                                         block_offset_batch,     &
                                         block_info_batch,       &
                                         NBATCH,NBLOCK,          &
                                         par_dist_block_list)

      CALL ZBLTP_IDC(ISMOST(1,ICSM),NSMST,IDC,blocktype)

!     transform CI vector
      call tracid(mo2mo_mat,work(kint1),lu_ref,lu_refout,        &
                  my_in_fh,my_out_fh,                            &
                  my_sc1_fh,my_sc2_fh,                           &
                  my_BVC_fh,                                     &
                  vec1,vec2,c2,                                  &
                  NBATCH,NBLOCK,blocks_per_batch,                &
                  batch_length,block_offset_batch,               &
                  block_info_batch,blocktype,                    &
                  par_dist_block_list,block_list,                &
                  rcctos,grouplist,proclist,                     &
                  file_info%iluxlist(1,file_info%                &
                  current_file_nr_active1),                      &
                  file_info%iluxlist(1,file_info%                &
                  current_file_nr_active2),                      &
                  file_info%ilublist,                            &
                  len_ilu1,len_ilu2,len_iluc,                    &
                  my_in_off,my_out_off,my_scr_off,               &
                  my_BVC_off)

#ifdef VAR_MPI
!     reset
      my_lu4_off = my_lu4_off_tmp
#endif

      deallocate(blocks_per_batch)
      deallocate(batch_length)
      deallocate(block_offset_batch)
      deallocate(block_info_batch)
      deallocate(blocktype)

  end subroutine return_bvec_transformed_2_new_mo_basis
!**********************************************************************

end module
