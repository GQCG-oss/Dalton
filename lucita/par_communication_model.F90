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

  implicit none

  public setup_communication_model

  private

  save

  integer, private                       :: istat(MPI_STATUS_SIZE)
  integer, private                       :: ierr

contains 
 
#ifdef VAR_MPI
C
C     parallel I/O mode for parallel runs / executable
C
      IIOMOD = 1
C     initialize number of process groups sharing a c-vector file
      NFLGRPS     = 0
C
#ifdef SYS_LOCALDISKS
      IF( IIOMOD .eq. 1 ) THEN
C
         CALL MPI_GET_PROCESSOR_NAME(MACHINENAME,NAMELENGTH,IERR)
C
         CALL FIND_GROUP_OF_PROCS(MACHINENAME,NAMELENGTH,
     &                            WORK(KPROCLIST),WORK(KPROCLENGTH))
C
         IF( NFLGRPS .eq. 1 ) THEN
            CALL ISETVC(WORK(KPROCLIST),1,LUCI_NMPROC)
         END IF
      END IF
#else
C     in case of a parallel filesystem, we should use it
      IF( IIOMOD .eq. 1 ) THEN
        NFLGRPS = 1
        CALL ISETVC(WORK(KPROCLIST),1,LUCI_NMPROC)
      END IF
#endif
C
C     compatibilty with LUCIAREL
C
      SHARED_M = .FALSE.
      IT_SHL   = - 1
      IC_SHL   = - 1
      CALL ISETVC(WORK(KGROUPLIST2),1,LUCI_NMPROC)

  subroutine setup_communication_model(nr_file_groups,                 &
                                       io_mode,                        &
                                       shared_memory_mode,             &
                                       shared_memory_lvl_ijkl,         &
                                       shared_memory_lvl_cvec,         &
                                       process_list_glb,               &
                                       process_list_shared_mem_glb,    &
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
!     2. setup communicators and process-id for the various communication levels
!        a. intra-node
!        b. inter-node
!        c. (shared-memory)
      call set_communication_levels()

  end subroutine setup_communication_model
!*******************************************************************************

  subroutine set_communication_levels(IGROUPLIST,IPROCLIST,IPROCLIST_SM)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
C
      DIMENSION IPROCLIST(LUCI_NMPROC),IGROUPLIST(LUCI_NMPROC)
      DIMENSION IPROCLIST_SM(LUCI_NMPROC)
      INTEGER IKEY, ICOLOR, JKEY, JCOLOR, KCOLOR, KKEY
      INTEGER LCOLOR, LKEY
      ITEST = 00
C     ITEST = 10
C
C     'intra-node' communicator MYNEW_COMM ( I/O communicator )
C
      IKEY = LUCI_MYPROC
      ICOLOR = IPROCLIST(LUCI_MYPROC+1)
C
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,ICOLOR,IKEY,MYNEW_COMM,IERR)
C
C     collect useful information about each group,
C     store on common block
C
      NEWCOMM_PROC = 0
      MYNEW_ID = 0
      CALL MPI_COMM_SIZE(MYNEW_COMM,NEWCOMM_PROC,IERR)
      CALL MPI_COMM_RANK(MYNEW_COMM,MYNEW_ID,IERR)
C
C     T-communicator MYNEW_COMM_SM ( 1st "shared memory" communicator)
C
      KKEY = LUCI_MYPROC
      IF( IT_SHL .eq. - 1 .or. IT_SHL .eq. - 2 )THEN
        KCOLOR = IPROCLIST_SM( LUCI_MYPROC + 1 )
      ELSE IF( IT_SHL .eq. 0 )THEN
        KCOLOR = IPROCLIST_SM( LUCI_MYPROC + 1 )
      ELSE IF( IT_SHL .eq. 1 )THEN
        KCOLOR = 1
      END IF
C
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,KCOLOR,KKEY,
     &                    MYNEW_COMM_SM,IERR)
C
C       collect useful information about each group,
C       store on common block
C
      NEWCOMM_PROC_SM = 0
      MYNEW_ID_SM     = 0
      CALL MPI_COMM_SIZE(MYNEW_COMM_SM,NEWCOMM_PROC_SM,IERR)
      CALL MPI_COMM_RANK(MYNEW_COMM_SM,MYNEW_ID_SM,IERR)
C
C     C-communicator MYNEW_COMM_SM_C ( 2nd "shared memory" communicator)
C
      LKEY = LUCI_MYPROC
      IF( IC_SHL .eq. - 1 )THEN
        LCOLOR = IPROCLIST_SM( LUCI_MYPROC + 1 )
      ELSE IF( IC_SHL .eq. 0 )THEN
        LCOLOR = IPROCLIST_SM( LUCI_MYPROC + 1 )
      ELSE IF( IC_SHL .eq. 1 )THEN
        LCOLOR = 1
      END IF
C
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,LCOLOR,LKEY,
     &                    MYNEW_COMM_SM_C,IERR)
C
C       collect useful information about each group,
C       store on common block
C
      NEWCOMM_PROC_SM_C = 0
      MYNEW_ID_SM_C     = 0
      CALL MPI_COMM_SIZE(MYNEW_COMM_SM_C,NEWCOMM_PROC_SM_C,IERR)
      CALL MPI_COMM_RANK(MYNEW_COMM_SM_C,MYNEW_ID_SM_C,IERR)
C
C     set node-master (might be MASTER of all CPU's)
      N_MASTER_SM_C = 0
C
C     inter-node communicator ICOMM
C
      IF( MYNEW_ID .eq. 0 ) THEN
        JKEY = IPROCLIST(LUCI_MYPROC+1)
        JCOLOR = 2
      ELSE
        JKEY = IPROCLIST(LUCI_MYPROC+1)
        JCOLOR = 3
      END IF
C
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,JCOLOR,JKEY,ICOMM,IERR)
C
C     collect again ...     
C
      ICOMM_ID = 0
      ICOMM_SIZE = 0
      CALL MPI_COMM_RANK(ICOMM,ICOMM_ID,IERR)
      CALL MPI_COMM_SIZE(ICOMM,ICOMM_SIZE,IERR)
C
C     set-up and store group information 
C
      CALL SET_GROUP_TABLE_REL(IGROUPLIST,IPROCLIST,ICOLOR)
C
C     store personal group number
C
      MY_GROUPN = ICOLOR
      IF ( ITEST .ge. 10 )THEN
        WRITE(LUWRT,*) ' '
        WRITE(LUWRT,*) ' OUTPUT FROM GROUP_CONSTRUCTOR_REL'
        WRITE(LUWRT,*) ' '
        WRITE(LUWRT,*) ' size of MYNEW_COMM     :',NEWCOMM_PROC
        WRITE(LUWRT,*) ' size of MYNEW_COMM_SM  :',NEWCOMM_PROC_SM
        WRITE(LUWRT,*) ' size of MYNEW_COMM_SM_C:',NEWCOMM_PROC_SM_C
        WRITE(LUWRT,*) ' size of ICOMM          :',ICOMM_SIZE
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SET_GROUP_TABLE_REL(IGROUPLIST,IPROCLIST,ICOLOR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
C
C     INPUT
      DIMENSION IPROCLIST(LUCI_NMPROC)
C     OUTPUT
      DIMENSION IGROUPLIST(LUCI_NMPROC)
C
      INUMB = 1
      DO IPROC = 1, LUCI_NMPROC
C
        IF( IPROCLIST(IPROC) .eq. ICOLOR ) THEN
C         put in the process tag ( master is 0 )
          IGROUPLIST(INUMB) = IPROC - 1
          INUMB = INUMB + 1
        END IF
C
      END DO
C
      IF( INUMB-1 .gt. NEWCOMM_PROC ) THEN
        WRITE(LUWRT,*) 'Error in SET_GROUP_TABLE_REL: more CPUs assigned
     & as included in this group!'
        WRITE(LUWRT,*) 'assigned CPUs, group size',INUMB-1,NEWCOMM_PROC
        CALL Abend2('Error detected in gp/gplupar.F: 
     &               SET_GROUP_TABLE_REL')
      END IF
C
      END
