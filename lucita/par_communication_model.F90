!dalton_copyright_start
!
!
!dalton_copyright_end

#ifdef VAR_MPI
module communication_model

! stefan: - this module provides all necessary functionality
!           to setup a communication model in parallel mcscf/ci 
!           calculations.
!
!           written by sknecht, may 2007 for DIRAC MCSCF/KR-CI/LUCITA
!           adapted for DALTON by sknecht, november 2010.
  use mpi

#if defined (VAR_INT64)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

  implicit none

  public setup_communication_model

  private

  save

  integer, private                       :: istat(MPI_STATUS_SIZE)
  integer, private                       :: ierr

contains 
 
  subroutine setup_communication_model(nr_file_groups,                 &
                                       io_mode,                        &
                                       shared_memory_mode,             &
                                       shared_memory_lvl_ijkl,         &
                                       shared_memory_lvl_cvec,         &
                                       process_list_glb,               &
                                       process_list_shared_mem_glb,    &
                                       my_process_id_glb,              &
                                       nr_of_process_glb,              &
                                       communicator_glb,               &
                                       joff)
!******************************************************************************
!
!    purpose:  initialize communication model for parallel CI/MCSCF runs:
!                - provides: 
!                           a. communication handles and process-ids in each group
!                           b. arrays with group listings
!                - requires: 
!                           a. global (starting) communication handle
!                           a. number of processors
!                           b. information concerning shared-memory mode(s)
!
!    for a detailed description of the communication model see references:
!    S. Knecht, H. J. Aa. Jensen, and T. Fleig
!       JCP, 128, 014108 (2008)
!       JCP, 132, 014108 (2008)
!
!*******************************************************************************
     integer, intent(in )   :: nr_of_process_glb
     integer, intent(in )   :: communicator_glb
     integer, intent(in )   :: my_process_id_glb
     integer, intent(out)   :: process_list_glb(nr_of_process_glb)
     integer, intent(out)   :: process_list_shared_mem_glb(nr_of_process_glb)
     integer, intent(out)   :: nr_file_groups
     integer, intent(out)   :: io_mode
     integer, intent(out)   :: shared_memory_lvl_ijkl
     integer, intent(out)   :: shared_memory_lvl_cvec
     logical, intent(out)   :: shared_memory_mode
     integer, intent(out)   :: joff
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!     initialize number of process groups sharing, e.g. a c-vector file
      nr_file_groups              =  0
!     parallel i/o mode in general (1 = MPI-I/O)
      io_mode                     =  1
!     set default values (to be consistent with Dirac/KR-CI where these values are dynamic.) 
      shared_memory_lvl_ijkl      = -1
      shared_memory_lvl_cvec      = -1
      process_list_shared_mem_glb =  1 
      shared_memory_mode          = .false.

!     1. determine # process groups and store group-IDs in process_list_glb
      call set_file_groups(process_list_glb,          &
                           nr_file_groups,            &
                           my_process_id_glb,         &
                           nr_of_process_glb,         &
                           communicator_glb)

      joff = -1
!     2. setup communicators and process-id for the various communication levels
!        a. intra-node
!        b. inter-node
!        c. (shared-memory)
!     call set_communication_levels()

  end subroutine setup_communication_model
!*******************************************************************************

  subroutine set_file_groups(process_list_glb,        &
                             nr_file_groups,          &
                             my_process_id_glb,       &
                             nr_of_process_glb,       &
                             communicator_glb)
!*******************************************************************************
     integer, intent(out)   :: process_list_glb(nr_of_process_glb)
     integer, intent(out)   :: nr_file_groups  
     integer, intent(in )   :: nr_of_process_glb
     integer, intent(in )   :: my_process_id_glb
     integer, intent(in )   :: communicator_glb
!-------------------------------------------------------------------------------
     integer                          :: process_name_length  
     integer                          :: local_counter_file_groups  
     integer                          :: current_process_id  
     integer                          :: proc_id  
     integer                          :: finished_loop  
     character (len=100)              :: process_name
     character (len=100)              :: scr_process_name
     integer,             allocatable :: scr_arr_name_length(:)
     character (len=100), allocatable :: scr_arr_process_name(:)
!-------------------------------------------------------------------------------
                                               
      process_list_glb = -1

!#ifdef SYS_LOCALDISKS                          
                                               
!     find system-dependent unique process name
      call mpixprocname(process_name,process_name_length)

      allocate(scr_arr_name_length(nr_of_process_glb))
      scr_arr_name_length = 0

!     1. gather all name length
      call mpi_allgather(process_name_length,1,my_MPI_INTEGER,         &
     &                   scr_arr_name_length,1,my_MPI_INTEGER,         &
     &                   communicator_glb,ierr)

!     2. collect all names in temporary storage
      allocate(scr_arr_process_name(nr_of_process_glb*100))
      scr_arr_process_name(my_process_id_glb+1) = process_name

      do proc_id = 1, nr_of_process_glb

         scr_process_name = process_name
         call mpi_bcast(scr_process_name,                              &
                        scr_arr_name_length(proc_id),                  &
                        MPI_CHARACTER,                                 &
                        proc_id-1,                                     &
                        communicator_glb,                              &
                        ierr)

         if(my_process_id_glb /= proc_id -1)then
           scr_arr_process_name(proc_id) = scr_process_name(1:scr_arr_name_length(proc_id))
         end if

      end do

 
!     3. find all processors on the same deck and reorder (if necessary) 
!     to get the processors as close as possible, starting with the master (id == 0)
!     NOTE: reordering is currently NOT done
      current_process_id        = 1
      local_counter_file_groups = 0

      do ! indefinite loop - we check at the end if we meet an exit condition

        finished_loop = 1

        local_counter_file_groups = local_counter_file_groups + 1
        process_name = scr_arr_process_name(current_process_id)         &
                       (1:scr_arr_name_length(current_process_id))
 
        do proc_id = 1, nr_of_process_glb
 
          scr_process_name(1:scr_arr_name_length(proc_id)) =            &
          scr_arr_process_name(proc_id)                                 &
          (1:scr_arr_name_length(proc_id))
 
          if(scr_process_name(1:scr_arr_name_length(proc_id)) ==        &
             process_name(1:scr_arr_name_length(current_process_id)))   &
             process_list_glb(proc_id) = current_process_id
 
        end do
!
!       check if we are done and set finished_loop: done = 1, else 0
!       search for the next lowest cpu building the 'local group master'
!       which is current_process_id
!
        do proc_id = 1, nr_of_process_glb
          if(finished_loop /= 0 .and. process_list_glb(proc_id) == -1)then 
            finished_loop = 0
            current_process_id = proc_id
          end if
        end do
        if(finished_loop == 1) exit

      end do
 
!     report the count of intra-groups formed by all processes
      if(my_process_id_glb == 0)then
        write(2,'(/a,i4,1x,a)') '  *** final output from the group generator:', &
        local_counter_file_groups,' intra-node groups are constructed. ***'
      end if

      nr_file_groups = local_counter_file_groups
 
!     final consistency check
      if(nr_file_groups == 1) process_list_glb = 1

      deallocate(scr_arr_process_name)
      deallocate(scr_arr_name_length)

!#else
      nr_file_groups   = 1
      process_list_glb = 1
!#endif

  end subroutine set_file_groups
!*******************************************************************************

! subroutine set_communication_levels(IGROUPLIST,IPROCLIST,IPROCLIST_SM)
!
!     IMPLICIT REAL*8 (A-H,O-Z)
!
!include "infpar.h"
!include "mpif.h"
!     DIMENSION ISTAT(MPI_STATUS_SIZE)
!include "parluci.h"
!
!     DIMENSION IPROCLIST(LUCI_NMPROC),IGROUPLIST(LUCI_NMPROC)
!     DIMENSION IPROCLIST_SM(LUCI_NMPROC)
!     INTEGER IKEY, ICOLOR, JKEY, JCOLOR, KCOLOR, KKEY
!     INTEGER LCOLOR, LKEY
!     ITEST = 00
!     ITEST = 10
!
!     'intra-node' communicator MYNEW_COMM ( I/O communicator )
!
!     IKEY = LUCI_MYPROC
!     ICOLOR = IPROCLIST(LUCI_MYPROC+1)
!
!     CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,ICOLOR,IKEY,MYNEW_COMM,IERR)
!
!     collect useful information about each group,
!     store on common block
!
!     NEWCOMM_PROC = 0
!     MYNEW_ID = 0
!     CALL MPI_COMM_SIZE(MYNEW_COMM,NEWCOMM_PROC,IERR)
!     CALL MPI_COMM_RANK(MYNEW_COMM,MYNEW_ID,IERR)
!
!     T-communicator MYNEW_COMM_SM ( 1st "shared memory" communicator)
!
!     KKEY = LUCI_MYPROC
!     IF( IT_SHL .eq. - 1 .or. IT_SHL .eq. - 2 )THEN
!       KCOLOR = IPROCLIST_SM( LUCI_MYPROC + 1 )
!     ELSE IF( IT_SHL .eq. 0 )THEN
!       KCOLOR = IPROCLIST_SM( LUCI_MYPROC + 1 )
!     ELSE IF( IT_SHL .eq. 1 )THEN
!       KCOLOR = 1
!     END IF
!
!     CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,KCOLOR,KKEY,
!    &                    MYNEW_COMM_SM,IERR)
!
!       collect useful information about each group,
!       store on common block
!
!     NEWCOMM_PROC_SM = 0
!     MYNEW_ID_SM     = 0
!     CALL MPI_COMM_SIZE(MYNEW_COMM_SM,NEWCOMM_PROC_SM,IERR)
!     CALL MPI_COMM_RANK(MYNEW_COMM_SM,MYNEW_ID_SM,IERR)
!
!     C-communicator MYNEW_COMM_SM_C ( 2nd "shared memory" communicator)
!
!     LKEY = LUCI_MYPROC
!     IF( IC_SHL .eq. - 1 )THEN
!       LCOLOR = IPROCLIST_SM( LUCI_MYPROC + 1 )
!     ELSE IF( IC_SHL .eq. 0 )THEN
!       LCOLOR = IPROCLIST_SM( LUCI_MYPROC + 1 )
!     ELSE IF( IC_SHL .eq. 1 )THEN
!       LCOLOR = 1
!     END IF
!
!     CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,LCOLOR,LKEY,
!    &                    MYNEW_COMM_SM_C,IERR)
!
!       collect useful information about each group,
!       store on common block
!
!     NEWCOMM_PROC_SM_C = 0
!     MYNEW_ID_SM_C     = 0
!     CALL MPI_COMM_SIZE(MYNEW_COMM_SM_C,NEWCOMM_PROC_SM_C,IERR)
!     CALL MPI_COMM_RANK(MYNEW_COMM_SM_C,MYNEW_ID_SM_C,IERR)
!
!     set node-master (might be MASTER of all CPU's)
!     N_MASTER_SM_C = 0
!
!     inter-node communicator ICOMM
!
!     IF( MYNEW_ID .eq. 0 ) THEN
!       JKEY = IPROCLIST(LUCI_MYPROC+1)
!       JCOLOR = 2
!     ELSE
!       JKEY = IPROCLIST(LUCI_MYPROC+1)
!       JCOLOR = 3
!     END IF
!
!     CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,JCOLOR,JKEY,ICOMM,IERR)
!
!     collect again ...     
!
!     ICOMM_ID = 0
!     ICOMM_SIZE = 0
!     CALL MPI_COMM_RANK(ICOMM,ICOMM_ID,IERR)
!     CALL MPI_COMM_SIZE(ICOMM,ICOMM_SIZE,IERR)
!
!     set-up and store group information 
!
!     CALL SET_GROUP_TABLE_REL(IGROUPLIST,IPROCLIST,ICOLOR)
!
!     store personal group number
!
!     MY_GROUPN = ICOLOR
!     IF ( ITEST .ge. 10 )THEN
!       WRITE(LUWRT,*) ' '
!       WRITE(LUWRT,*) ' OUTPUT FROM GROUP_CONSTRUCTOR_REL'
!       WRITE(LUWRT,*) ' '
!       WRITE(LUWRT,*) ' size of MYNEW_COMM     :',NEWCOMM_PROC
!       WRITE(LUWRT,*) ' size of MYNEW_COMM_SM  :',NEWCOMM_PROC_SM
!       WRITE(LUWRT,*) ' size of MYNEW_COMM_SM_C:',NEWCOMM_PROC_SM_C
!       WRITE(LUWRT,*) ' size of ICOMM          :',ICOMM_SIZE
!     END IF
!
!     END
!**********************************************************************
!                                                                     *
! LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
!           parallelization by Stefan Knecht                          *
!                                                                     *
!**********************************************************************
!     SUBROUTINE SET_GROUP_TABLE_REL(IGROUPLIST,IPROCLIST,ICOLOR)
!
!     IMPLICIT REAL*8 (A-H,O-Z)
!
!include "infpar.h"
!include "mpif.h"
!     DIMENSION ISTAT(MPI_STATUS_SIZE)
!include "parluci.h"
!
!     INPUT
!     DIMENSION IPROCLIST(LUCI_NMPROC)
!     OUTPUT
!     DIMENSION IGROUPLIST(LUCI_NMPROC)
!
!     INUMB = 1
!     DO IPROC = 1, LUCI_NMPROC
!
!       IF( IPROCLIST(IPROC) .eq. ICOLOR ) THEN
!         put in the process tag ( master is 0 )
!         IGROUPLIST(INUMB) = IPROC - 1
!         INUMB = INUMB + 1
!       END IF
!
!     END DO
!
!     IF( INUMB-1 .gt. NEWCOMM_PROC ) THEN
!       WRITE(LUWRT,*) 'Error in SET_GROUP_TABLE_REL: more CPUs assigned
!    & as included in this group!'
!       WRITE(LUWRT,*) 'assigned CPUs, group size',INUMB-1,NEWCOMM_PROC
!       CALL Abend2('Error detected in gp/gplupar.F: 
!    &               SET_GROUP_TABLE_REL')
!     END IF
!
!     END
#else 
module dummy_comm_model
#endif
end module
