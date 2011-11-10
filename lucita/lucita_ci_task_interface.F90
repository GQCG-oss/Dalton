!dalton_copyright_start
!
!
!dalton_copyright_end

module lucita_ci_task_interface

  use lucita_integral_density_interface

  implicit none

#if defined (VAR_INT64)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

  public CI_task_list_interface
  public create_CI_task_list

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

        case ('return CIdia', 'perform Davi', 'return sigma', 'return rotVC', 'return densM')

!         fill integral indices pointers and possibly read in integrals
!         -------------------------------------------------------------
          call lucita_pointer_integral_driver(ci_task_list(ci_task_ticket+1),          &
                                              int1_or_rho1,                            &
                                              int2_or_rho2,                            &
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
                                    cref,hc,resolution_mat)
          case ('return sigma')
            call return_sigma_vector(print_lvl,nbatch,block_list,par_dist_block_list,             &
                                     rcctos,grouplist,proclist,                                   &
                                     cref,hc,resolution_mat,int1_or_rho1,int2_or_rho2)
          case ('return rotVC')
            call traci_ctl(cref,hc,resolution_mat,int1_or_rho1,par_dist_block_list)
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

        case default
         call quit('error in create_CI_task_list: undefined CI run id.')
      end select

  end subroutine create_CI_task_list
!**********************************************************************

  subroutine calculate_CI_diagonal(cref,                         &
                                   print_lvl,                    &
                                   block_list,                   &
                                   par_dist_block_list)

  use lucita_vector_partitioning_pointer
  use lucita_energy_types

#ifdef VAR_MPI
  use par_mcci_io
#ifndef VAR_USE_MPIF
  use mpi
#else
#include "mpif.h"
#endif
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
      real(8)                  :: work
#include "wrkspc.inc"
#ifdef VAR_MPI
#include "maxorb.h"
#include "infpar.h"
      dimension                :: my_STATUS(MPI_STATUS_SIZE)
      integer(MPI_OFFSET_KIND) :: offset
      integer                  :: block_length
      integer                  :: ierr
      integer                  :: isblk
      integer                  :: iproc
      integer, allocatable     :: dialist(:)
#endif
!-------------------------------------------------------------------------------

!     construct the diagonal part of the Hamiltonian
      if(idiag == 2)  ludia = lusc1
      if(icistr >= 2) call rewino(ludia)

      call gasdiat(cref,ludia,ecore_orig-ecore,icistr,i12,       &
                   work(klcbltp),num_blocks,work(klcibt),        &
                   par_dist_block_list)

      if(nocsf == 1 .and. icistr == 1)then
        call rewino(ludia)
        call todsc_luci(cref,l_combi,l_combi,ludia)
      end if

      if(luci_nmproc > 1)then
#ifdef VAR_MPI
        rewind ludia

        allocate(dialist(iall_lu1))
        dialist = 1

!       collect H diagonal 
        call mcci_cp_vcd_mpi_2_seq_io_interface(cref,ludia,idia,            &
                                                my_dia_off,                 &
                                                dialist,                    &
                                                par_dist_block_list,        &
                                                block_list,                 &
                                                MPI_COMM_WORLD,             &
                                                num_blocks,1,1,2)
       deallocate(dialist)
#endif
      end if

#ifdef LUCI_DEBUG
!     debug printing (parallel run)
      if(print_lvl .ge. 10)then
#ifdef VAR_MPI
        if(luci_nmproc > 1)then
          offset = my_dia_off
          do isblk = 1, num_blocks
            call get_block_proc(par_dist_block_list,isblk,iproc)
            if(luci_myproc == iproc)then
              call get_block_length(block_list,isblk,block_length)
              call mpi_file_read_at(idia,offset,cref,block_length,   &
                                    MPI_REAL8,my_STATUS,ierr)
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

  use lucita_vector_partitioning_pointer

#include "parluci.h"
#include "mxpdim.inc"
#include "crun.inc"
#include "cstate.inc"
#include "clunit.inc"
!-------------------------------------------------------------------------------
      real(8), intent(inout)   :: cref(*)
      integer, intent(in)      :: print_lvl
!-------------------------------------------------------------------------------
      real(8)                  :: work
#include "wrkspc.inc"
      integer                  :: eigen_state_id
      integer                  :: imzero, iampack
      integer                  :: lusc_vector_file
!-------------------------------------------------------------------------------

      call rewino(luc)

      if(luci_myproc == luci_master)then
        do eigen_state_id = 1, nroot

          if(icistr == 1)then
            call frmdsc_luci(cref,l_combi,l_combi,luc,imzero,iampack)
          else
            lusc_vector_file = luc
            if(nroot > 1)then
              call rewino(lusc1)
              call copvcd(luc,lusc1,cref,0,-1)
              lusc_vector_file = lusc1
            end if
          end if
!
          write(luwrt,'(/a)')    '   ******************************************************'
          write(luwrt,'(a,i3)')  '   Analysis of leading coefficients for eigen state = ',eigen_state_id
          write(luwrt,'(a)')     '   ******************************************************'


          call rewino(lusc_vector_file)
          call gasana(cref,num_blocks,work(klcibt),work(klcbltp),lusc_vector_file,icistr)
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
  use lucita_file_list_pointer
  use lucita_cfg
  use lucita_mcscf_ci_cfg
#ifdef VAR_MPI
  use par_mcci_io
#ifndef VAR_USE_MPIF
  use mpi
#else
#include "mpif.h"
#endif
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
      real(8)                          :: exps2
      real(8)                          :: cv_dummy, hcv_dummy
      integer, parameter               :: isigden = 2 ! density matrix switch for sigden_ci
      integer                          :: eigen_state_id
      integer                          :: imzero, iampack
      integer                          :: lusc_vector_file
      integer                          :: luhc_vector_file
      integer                          :: twopart_densdim
      integer                          :: rhotype
      integer                          :: idum
      integer(8)                       :: kdum
      integer(8)                       :: k_dens2_scratch
      integer(8)                       :: k_scratch1
      integer(8)                       :: k_scratch2
      integer(8)                       :: k_scratch3
      integer(8)                       :: k_scratch4
#ifdef VAR_MPI
      integer, allocatable             :: lu1list(:)
      integer, allocatable             :: lu2list(:)
      integer, allocatable             :: luclist(:)
      integer, allocatable             :: blocks_per_batch(:)
      integer, allocatable             :: batch_length(:)
      integer, allocatable             :: block_offset_batch(:)
      integer, allocatable             :: block_info_batch(:,:)
      integer                          :: iatp, ibtp
      integer                          :: nbatch_par, nblock_par
      integer                          :: ierr
      integer(kind=MPI_OFFSET_KIND)    :: my_lu4_off_tmp
#endif
!-------------------------------------------------------------------------------

!     set level of particle-density matrix calculation
      i12 = idensi

!     set rhotype
      rhotype = 1
      if(lucita_cfg_transition_densm) rhotype = 3

#ifdef VAR_MPI
      if(luci_nmproc > 1)then

        if(icsm /= issm) & 
        call quit('*** return_Xp_density_matrix: property/response density matrix not parallelized yet. ***')

!       sequential --> MPI I/O
!       ----------------------

        allocate(lu1list(iall_lu1))
        allocate(lu2list(iall_lu2))
        allocate(luclist(iall_luc))
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
        lu1list  = 0
        lu2list  = 0

        call rewino(luc)
!       step 1: the rhs vector
        call mcci_cp_vcd_mpi_2_seq_io_interface(cref,LUC,ILU1,               &
                                                MY_LU1_OFF,                  &
                                                lu1list,                     &
                                                par_dist_block_list,         &
                                                block_list,                  &
                                                MPI_COMM_WORLD,              &
                                                NUM_BLOCKS,NROOT,1,1)

        if(lucita_cfg_transition_densm)then
          call rewino(luhc)
!         step 2: the lhs vector
          call mcci_cp_vcd_mpi_2_seq_io_interface(hc,LUHC,ILU2,              &
                                                  MY_LU2_OFF,                &
                                                  lu2list,                   &
                                                  par_dist_block_list,       &
                                                  block_list,                &
                                                  MPI_COMM_WORLD,            &
                                                  NUM_BLOCKS,NROOT,1,1)
          luhc_vector_file = ilu2
        else
          luhc_vector_file = ilu1
          iall_lu2         = iall_lu1
          call icopy(iall_lu2,lu1list,1,lu2list,1)
!         save and temporarily re-set the file offset for ilu2 to ilu1 (ilu2 has the same offset as ilu4 which is used inside sigden_ci)
          my_lu4_off_tmp   = my_lu4_off
          my_lu4_off       = my_lu1_off
        end if
        lusc_vector_file = iluc 
      end if
#endif

      call rewino(luhc)
      call rewino(luc)

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
            call frmdsc_luci(cref,l_combi,l_combi,luc,imzero,iampack)
            if(lucita_cfg_transition_densm)then
!             transition density
              call frmdsc_luci(hc,l_combi,l_combi,luhc,imzero,iampack)
            else
              call dcopy(l_combi,cref,1,hc,1)
            end if
          else
            if(lucita_cfg_transition_densm)then
!             transition density
              lusc_vector_file = luc
              luhc_vector_file = luhc
            else
              call rewino(lusc1)
              call copvcd(luc,lusc1,cref,0,-1)
              call copvcd(lusc1,luhc,cref,1,-1)
              call rewino(lusc1)
              call rewino(luhc)
              lusc_vector_file = lusc1
              luhc_vector_file = luhc
            end if
          end if
        else
#ifdef VAR_MPI
!         MPI I/O --> MPI I/O node-master collection file
!         -----------------------------------------------
          call mpi_barrier(mynew_comm,ierr) ! probably not needed if collection of densities is in action
          luclist = 0
          call mcci_cp_vcd_batch(ilu1,iluc,cref,nbatch_par,blocks_per_batch,               &
                                 batch_length,block_offset_batch,block_info_batch,         &
                                 block_list,my_lu1_off,my_luc_off,lu1list,luclist,         &
                                 my_vec2_ioff,my_act_blk2,eigen_state_id-1)
!         offset in MPI I/O file for lhs-vector
          jvec_sf = eigen_state_id - 1
#endif
        end if
!
        call sigden_ci(cref,hc,resolution_mat,lusc_vector_file,luhc_vector_file,           &
                       cv_dummy,hcv_dummy,isigden,rhotype                                  &
#ifdef VAR_MPI
                      ,luclist,lu2list,block_list,par_dist_block_list,                     &
                       rcctos,grouplist,proclist                                           &
#endif
                        )
        if(luci_nmproc > 1)then
#ifdef VAR_MPI
!         collect density matrices
!         ------------------------
          if(luci_myproc == luci_master)then
            call mpi_reduce(mpi_in_place,work(krho1),nacob**2,mpi_real8,             &
                            mpi_sum,luci_master,mpi_comm_world,ierr)
            if(i12 > 1)then
              call mpi_reduce(mpi_in_place,work(krho2),nacob**2*(nacob**2+1)/2,      &
                              mpi_real8,mpi_sum,luci_master,mpi_comm_world,ierr)
            end if
          else
            call mpi_reduce(work(krho1),mpi_in_place,nacob**2,mpi_real8,             &
                            mpi_sum,luci_master,mpi_comm_world,ierr)
            if(i12 > 1)then
              call mpi_reduce(work(krho2),mpi_in_place,nacob**2*(nacob**2+1)/2,      &
                              mpi_real8,mpi_sum,luci_master,mpi_comm_world,ierr)
            end if
          end if
#endif
        end if ! luci_nmproc > 1

!       export 1-/2-particle density matrix to MCSCF format
!       ---------------------------------------------------
!
        idum = 0
        call memman(idum,idum,'MARK  ',idum,'Xpden1')
        call memman(k_dens2_scratch,nacob**4,'ADDL  ',2,'PVfull')

!       activate for MCSCF/improved CI
        if(integrals_from_mcscf_env)then
          call lucita_putdens_generic(work(krho1),work(krho2),int1_or_rho1,int2_or_rho2,             &
                                      work(k_dens2_scratch),i12,isigden,rhotype,eigen_state_id)
        else
          twopart_densdim = 0
          if(i12 > 1) twopart_densdim = (nacob*(nacob+1)/2)**2
          call memman(k_scratch1     ,(nacob*(nacob+1)/2),'ADDL  ',2,'scrat1')
          call memman(k_scratch2     ,twopart_densdim    ,'ADDL  ',2,'scrat2')

          call lucita_putdens_generic(work(krho1),work(krho2),work(k_scratch1),work(k_scratch2),     &
                                      work(k_dens2_scratch),i12,isigden,rhotype,eigen_state_id)
        end if
!
        call memman(kdum ,idum,'FLUSM ',2,'Xpden1')

!       natural orbital occupation numbers
!       ----------------------------------
        if(lucita_cfg_natural_orb_occ_nr)then
          idum = 0
          call memman(idum,idum,'MARK  ',idum,'Xpden2')
          call memman(k_scratch1,nacob**2         ,'ADDL  ',2,'scrat1')
          call memman(k_scratch2,nacob**2         ,'ADDL  ',2,'scrat2')
          call memman(k_scratch3,nacob            ,'ADDL  ',2,'scrat3')
          call memman(k_scratch4,nacob*(nacob+1)/2,'ADDL  ',2,'scrat4')
#ifdef LUCI_DEBUG
          write(luwrt,'(/a)') ' IPRDEN lowered explicitly in return_Xp_density_matrix '
          iprden = 1
#endif

          call lnatorb(work(krho1),nsmob,ntoobs,nacobs,ninobs,                                     &
                       ireost,work(k_scratch1),work(k_scratch2),work(k_scratch3),nacob,            &
                       work(k_scratch4),1)

          call memman(kdum ,idum,'FLUSM ',2,'Xpden2')
        end if

      end do ! loop over eigen states

#ifdef VAR_MPI
      if(luci_nmproc > 1)then
!       reset file off-set and free memory
        if(.not.lucita_cfg_transition_densm) my_lu4_off = my_lu4_off_tmp
        deallocate(block_info_batch)
        deallocate(block_offset_batch)
        deallocate(batch_length)
        deallocate(blocks_per_batch)
        deallocate(luclist)
        deallocate(lu2list)
        deallocate(lu1list)
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
  use lucita_file_list_pointer
  use lucita_cfg
  use lucita_mcscf_ci_cfg
#ifdef VAR_MPI
  use par_mcci_io
#ifndef VAR_USE_MPIF
  use mpi
#else
#include "mpif.h"
#endif
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
      real(8)                          :: cv_dummy, hcv_dummy
      integer                          :: eigen_state_id
      integer                          :: imzero, iampack
      integer                          :: lusc_vector_file
      integer                          :: luhc_vector_file
      integer, parameter               :: isigden = 1 ! sigma vector switch for sigden_ci
#ifdef VAR_MPI
      integer, allocatable             :: lu1list(:)
      integer, allocatable             :: lu2list(:)
      integer, allocatable             :: luclist(:)
      integer, allocatable             :: blocks_per_batch(:)
      integer, allocatable             :: batch_length(:)
      integer, allocatable             :: block_offset_batch(:)
      integer, allocatable             :: block_info_batch(:,:)
      integer                          :: iatp, ibtp
      integer                          :: nbatch_par, nblock_par
      integer                          :: ierr
#endif
!-------------------------------------------------------------------------------

#ifdef VAR_MPI
      if(luci_nmproc > 1)then

        if(icsm /= issm) & 
        call quit('*** return_sigma_vector: property/response E2b calculation not parallelized yet. ***')

!       sequential --> MPI I/O
!       ----------------------

        allocate(lu1list(iall_lu1))
        allocate(lu2list(iall_lu2))
        allocate(luclist(iall_luc))
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
        lu1list = 0
        lu2list = 0

        call rewino(luc)
!       step 1: the rhs vector
        call mcci_cp_vcd_mpi_2_seq_io_interface(cref,LUC,ILU1,               &
                                                MY_LU1_OFF,                  &
                                                lu1list,                     &
                                                par_dist_block_list,         &
                                                block_list,                  &
                                                MPI_COMM_WORLD,              &
                                                NUM_BLOCKS,NROOT,1,1)
      end if
#endif

      call rewino(luhc)
      call rewino(luc)

      do eigen_state_id = 1, nroot

        write(luwrt,'(/a)')    '   **********************************************'
        write(luwrt,'(a,i3)')  '   calculating E2b lhs vector for b state = ', eigen_state_id
        write(luwrt,'(a)')     '   **********************************************'

        if(luci_nmproc == 1)then
          if(icistr == 1)then
            call frmdsc_luci(cref,l_combi,l_combi,luc,imzero,iampack)
          else
            call rewino(lusc1)
            call copvcd(luc,lusc1,cref,0,-1)
            call rewino(lusc1)
          end if
          lusc_vector_file = lusc1
          luhc_vector_file = luhc
        else
#ifdef VAR_MPI
!         MPI I/O --> MPI I/O node-master collection file
!         -----------------------------------------------
          call mpi_barrier(mynew_comm,ierr) ! probably not needed if collection of densities is in action
          luclist = 0
          call mcci_cp_vcd_batch(ilu1,iluc,hc,nbatch_par,blocks_per_batch,                 &
                                 batch_length,block_offset_batch,block_info_batch,         &
                                 block_list,my_lu1_off,my_luc_off,lu1list,luclist,         &
                                 my_vec2_ioff,my_act_blk2,eigen_state_id-1)

          lusc_vector_file = iluc
          luhc_vector_file = ilu2

!         offset in MPI I/O file for lhs-vector
          jvec_sf = eigen_state_id - 1
#endif
        end if
!
        call sigden_ci(cref,hc,resolution_mat,lusc_vector_file,luhc_vector_file,           &
                       cv_dummy,hcv_dummy,isigden,-1                                       &
#ifdef VAR_MPI
                      ,luclist,lu2list,block_list,par_dist_block_list,                     &
                       rcctos,grouplist,proclist                                           &
#endif
                        )

      end do ! loop over eigen states

      if(luci_nmproc > 1)then
#ifdef VAR_MPI
!       collect e2b vector(s)
!       --------------------
        call rewino(luhc)

!       the lhs vector
        call mcci_cp_vcd_mpi_2_seq_io_interface(hc,luhc,luhc_vector_file,  &
                                                my_lu2_off,                &
                                                lu2list,                   &
                                                par_dist_block_list,       &
                                                block_list,                &
                                                MPI_COMM_WORLD,            &
                                                num_blocks,nroot,1,2)

        deallocate(block_info_batch)
        deallocate(block_offset_batch)
        deallocate(batch_length)
        deallocate(blocks_per_batch)
        deallocate(luclist)
        deallocate(lu2list)
        deallocate(lu1list)
#endif
      end if ! luci_nmproc > 1

  end subroutine return_sigma_vector
!**********************************************************************

end module
