!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!**********************************************************************

      SUBROUTINE LUCITA(run_id,                                         &
                        cref,                                           &
                        hc,                                             &
                        int1_or_rho1,                                   &
                        int2_or_rho2,                                   &
                        WRK_DALTON,                                     &
                        LWRK_DALTON)
!
!     GASCI program
!
!     Written by Jeppe Olsen , winter of 1991
!                          GAS version in action summer of 95
! 
!     General restructuring and parallel version: Stefan Knecht spring 2007 - ???
!
!     This routine is originally called from dalton sirius.
!     In the case of a parallel run, the relevant nodes will be 
!     woken up and sent to sleep.
!
      use lucita_mcscf_ci_cfg
      use lucita_cfg
!     parallel lucita
#ifdef VAR_MPI
      use sync_coworkers
#endif
      implicit none

!*******************************************************************************
      integer,           intent(inout) :: LWRK_DALTON
      real(8),           intent(inout) :: WRK_DALTON(LWRK_DALTON)
      character(len=12), intent(in)    :: run_id
!------------- optional input depending on MCSCF/CI run ------------------------
      real(8),           intent(inout) :: cref(*)
      real(8),           intent(inout) :: hc(*)
      real(8),           intent(inout) :: int1_or_rho1(*)
      real(8),           intent(inout) :: int2_or_rho2(*)
!------------- end of optional input depending on MCSCF/CI run -----------------
!
!     stefan: 'optional' is meant such that all parameters have to appear in the 
!             calling list to LUCITA (since we call it from outside a .F90 file) 
!             but it is not guaranteed that these are actually allocated 
!             (memory taken from the Dalton WORK array; if it had been for f90 allocation 
!             in Dalton i could have placed all the stuff in a module... ;))
!             Odense, July 2011.
!             
!

#ifdef VAR_MPI
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#endif
#include "maxorb.h"
#include "infpar.h"
#include "parluci.h"
#include "priunit.h"

      integer :: k1, k_resolution_mat, kfree, kfrsav, lfree
!------------------------------------------------------------------------------

      kfrsav = 1
      kfree  = kfrsav
      lfree  = lwrk_dalton

!     set LUCITA internal master and co-worker IDs
      LUCI_MASTER = MASTER
      LUCI_MYPROC = MYNUM
!     add the master node to obtain the total number of active processes
      LUCI_NMPROC = NODTOT + 1

!     set the ci-task id
      lucita_ci_run_id = run_id
 
#ifdef VAR_MPI

      IF (LUCI_NMPROC .GT. 1) then 
!       summon the co-workers, who are waiting in the general menu routine
        call lucita_start_cw

!       control variables for synchronizing co-workers
!       ----------------------------------------------
!       b. synchronize co-workers with ci/mcscf input data: dynamic and static (the latter only if required)
        sync_ctrl_array(1) = .true.
!       c. synchronize MC/CI array lengths with co-workers
        sync_ctrl_array(2) = .true.
!       a. synchronize ci-run id with co-workers
        sync_ctrl_array(6) = .true.

!       start synchronization process
        call sync_coworkers_cfg()
      end if
#endif
 
!     set marker on incoming work space 
      call memget('REAL',k1,1,wrk_dalton,kfree,lfree)
!     allocate resolution matrix (may have length 0 depending on the lucita_ci_run_id)
      call memget('REAL',k_resolution_mat,len_resolution_mat_mc2lu,     &
                  wrk_dalton,kfree,lfree)

!     enter the LUCITA driver (master becomes now primus inter pares -- "first among equals")
!     ---------------------------------------------------------------------------------------
      call lucita_driver(wrk_dalton(kfree),                             &
                         lfree,                                         &
                         cref,                                          &
                         hc,                                            &
                         wrk_dalton(k_resolution_mat),                  &
                         int1_or_rho1,                                  &
                         int2_or_rho2)

!     release marker on incoming work space
      call memrel('lucita.done',wrk_dalton,kfrsav,kfrsav,kfree,lfree)
      
#ifdef LUCI_DEBUG
!     print info on matrix multiplier
      call pr_matml_stat
#endif
 
      end
!**********************************************************************
 
      SUBROUTINE LUCITA_DRIVER(WRK_DALTON,                              &
                               LWRK_DALTON,                             &
                               cref,                                    &
                               hc,                                      &
                               resolution_mat,                          &
                               int1_or_rho1,                            &
                               int2_or_rho2)
!-----------------------------------------------------------------------
!
!     purpose: LUCITA main driver providing branching points for
!              individual MCSCF/CI tasks
!
!-----------------------------------------------------------------------
!
!     dependencies on F90 modules
!     common modules
      use lucita_cfg
      use lucita_setup
      use lucita_file_list_pointer
      use lucita_vector_partitioning_pointer
      use lucita_ci_task_interface
      use lucita_energy_types
!     parallel lucita
#ifdef VAR_MPI
      use parallel_setup
      use communication_model
      use file_io_model
      use sync_coworkers
#endif

#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "infpar.h"
#ifdef VAR_MPI
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
#include "files.inc"
! parameters for dimensioning
#include "mxpdim.inc"
! file numbers
#include "clunit.inc"
#include "units.inc"
! print flags
#include "cprnt.inc"
#include "lucinp.inc"
#include "cstate.inc"
! contains enviro
#include "crun.inc"
#include "cicisp.inc"
#include "oper.inc"
#include "cgas.inc"
#include "cands.inc"
#include "glbbas.inc"
! memory
#include "wrkspc.inc"
!------------------------------------------------------------------------------
      integer, intent(inout) :: lwrk_dalton
      real(8), intent(inout) :: wrk_dalton(lwrk_dalton)
!------------- optional input depending on MCSCF/CI run ------------------------
      real(8),           intent(inout) :: cref(*)
      real(8),           intent(inout) :: hc(*)
      real(8),           intent(inout) :: resolution_mat(*)
      real(8),           intent(inout) :: int1_or_rho1(*)
      real(8),           intent(inout) :: int2_or_rho2(*)
!------------- end of optional input depending on MCSCF/CI run -----------------
      integer                :: files_to_close
      integer                :: fh_offset
      integer, allocatable   :: proclist(:)
      integer, allocatable   :: grouplist(:)
      integer, allocatable   :: rcctos(:)
      integer, allocatable   :: fh_array(:)
      integer, allocatable   :: par_dist_block_list(:)
      integer, allocatable   :: block_list(:)
!------------------------------------------------------------------------------

      CALL QENTER('LUCITA_MASTER')

      write(lupri,'(/a)')                                               &
      '   ------------------------------------------------'
      write(lupri,'(a,a12)')                                            &
      '   LUCITA module started with run-id = ',lucita_ci_run_id
      write(lupri,'(a/)')                                               &
      '   ------------------------------------------------'

!     set LUCITA internal control and task-individual common block information
!     ------------------------------------------------------------------------
      call setup_lucita_cb_interface()

!     print header
      if(lucita_ci_run_id == 'standard ci ') call hello_dalton_lucita

!     allocate process-id and grouplist-id arrays
      allocate(proclist(luci_nmproc))
      allocate(grouplist(luci_nmproc))
      proclist  = -1
      grouplist = -1

!     set LUCITA internal orbital and string common block information 
!     ---------------------------------------------------------------
      call setup_lucita_orbital_string_cb(lucita_cfg_initialize_cb,     &
                                          iprorb)
!     initialize LUCITA work space
!     ----------------------------
      call setup_lucita_initialize_work_space(wrk_dalton,               &
                                              work,                     &
                                              lwrk_dalton,              &
                                              mxpwrd,                   &
                                              luci_nmproc,              &
                                              nr_files,                 &
                                              luwrt)
 
!     set LUCITA internal work space pointers 
!     ---------------------------------------
!     part 1 - generate string and integral pointers on work space + possibly read in integrals.
      call setup_lucita_pointer_strings_work_space(icspc,icsm,iprorb)

!     information about largest vector block and block structure ==> memory for partitioning of C vector
!     --------------------------------------------------------------------------------------------------
      iatp   = 1
      ibtp   = 2
      call z_blkfo(icspc,icsm,iatp,ibtp,klclbt,klclebt,                 &
                   klci1bt,klcibt,klcbltp,nbatch,nblock)

!     allocate arrays for:
!      a. block length
!      b. block distribution among cpu's (parallel run)
!      c. the c to sigma connections     (parallel run)
!      d. storing file handles           (parallel run)
      allocate(         block_list(nblock))
      allocate(par_dist_block_list(nblock))
      allocate(             rcctos(nblock))
      allocate(           fh_array(nr_files))

!     initialize
      call icopy(nblock,work(klci1bt),1,block_list,1)
      par_dist_block_list = -2
      rcctos              =  0
      fh_array            =  0

      if(luci_nmproc > 1)then
#ifdef VAR_MPI

!       setup parallel model:
!         1. communication model
!         2. file i/o model
!         3. static load balancing (== block distribution among cpu's)
!         4. file offset calculation
!         5. file list pointers creation
!       --------------------------------------------------------------
        call lucita_setup_parallel_model(block_list,par_dist_block_list,&
                                         grouplist,proclist,nblock,     &
                                         rcctos,kiluclist,              &
                                         kilu1list,kilu2list,kilu3list, &
                                         kilu4list,kilu5list,kilu6list, &
                                         kilu7list,fh_array,mxciv,nroot)
#endif
      else
        call setup_lucita_par_dist_in_seq(par_dist_block_list,          &
                                          block_list,nblock)
      end if

!     ----------------------------------------------------------------------------
!      branching point for MCSCF/CI tasks (controlled by entries in ci_task_list)
!     ----------------------------------------------------------------------------
      call CI_task_list_interface(ci_task_list,                         &
                                  cref,                                 &
                                  hc,                                   &
                                  resolution_mat,                       &
                                  int1_or_rho1,                         &
                                  int2_or_rho2,                         &
                                  block_list,                           &
                                  par_dist_block_list,proclist,         &
                                  grouplist,fh_array,rcctos,nbatch,     &
                                  nblock,iprorb)
!     ----------------------------------------------------------------------------------
!      end of branching point for MCSCF/CI tasks (controlled by entries in ci_task_list)
!     ----------------------------------------------------------------------------------

#ifdef VAR_MPI
      if(luci_nmproc > 1)then
!       close parallel model:
!         1. parallel file(s) / file handle(s)
!         2. communication model
        call lucita_close_parallel_model(nr_files,lucita_ci_run_id,     &
                                         fh_array,mynew_comm,icomm,     &
                                         mynew_comm_sm, mynew_comm_sm_c)
      end if
#endif

!     free memory
      deallocate(fh_array)
      deallocate(rcctos)
      deallocate(par_dist_block_list)
      deallocate(block_list)
      deallocate(grouplist)
      deallocate(proclist)
     
      write(lupri,'(/a)')                                               &
      '   ----------------------------------------------------------'
      write(lupri,'(a,a12,a)') '   LUCITA run, run-id = ',              &
            lucita_ci_run_id(1:12),' finished with no errors.'
      write(lupri,'(a )')                                               &
      '   ----------------------------------------------------------'

      call qexit('LUCITA_MASTER')

      END
!**********************************************************************
#ifdef VAR_MPI

      subroutine lucita_start_cw
!
!     wake-up the co-workers
!**********************************************************************
      implicit none
#include "parluci.h"
      integer :: iprtyp
      integer :: idummy
!
!     iprtyp = 42 for parallel LUCITA
!
      iprtyp =  42
      idummy = - 1
      call mpixbcast(iprtyp,1,'INTEGER',luci_master)
      call mpixbcast(idummy,1,'INTEGER',luci_master)
!
      end
!**********************************************************************

      subroutine lucita_coworker_main(work_dalton,lwork_dalton)
!
!     LUCITA routine (DALTON interface) for the co-workers
!
!     adapted version of the co-workers LUCINOD in Dirac
!     LUCITA version by Stefan Knecht, Feb. 06
!
!***********************************************************************
      use sync_coworkers
      use lucita_mcscf_ci_cfg
      use lucita_cfg
      use lucita_setup
      implicit none
#include "priunit.h"
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
!----------------------------------------------------------------------
      integer, intent(in)    :: lwork_dalton
      real(8), intent(inout) :: work_dalton(lwork_dalton)
      character(len=20)      :: lucitabasf
      character(len=24)      :: lucifiln
      integer                :: outfile_node
      integer                :: k1, kfree, kfrsav, lfree
      integer                :: kcref, khc, kresolution_mat 
      integer                :: kint1_or_rho1, kint2_or_rho2
      integer                :: lupri_save, lufil, print_lvl
!----------------------------------------------------------------------

      call qenter('lucita_coworker_main')

      kfrsav = 1
      kfree  = kfrsav
      lfree  = lwork_dalton

!     arrange for the MPI stuff and correct node number
!     to the total number of running invocations.
      LUCI_MASTER = MASTER
      LUCI_MYPROC = MYNUM
!     Add the master node, NODTOT = number of co-workers
      LUCI_NMPROC = NODTOT + 1
!
#ifdef LUCI_DEBUG
      print *, '*** co-worker',LUCI_MYPROC,'has arrived. ***'
      print *, '*** co-worker: priunit =',lupri,'***'
#endif

!     control variables for synchronizing co-workers
!     ----------------------------------------------
!     b. synchronize co-workers with ci/mcscf input data: dynamic and static (the latter only if required)
      sync_ctrl_array(1) = .true.
!     c. synchronize MC/CI array lengths with co-workers
      sync_ctrl_array(2) = .true.
!     a. synchronize ci-run id with co-workers
      sync_ctrl_array(6) = .true.

!     start synchronization process
      call sync_coworkers_cfg()

!     create a node-unique filename as output file. Important on
!     shared file systems. Otherwise all the output gets mingled in one
!     file. You don't really want to do this.
!
      lucitabasf   = "lucita-coworkers.out"
      lupri_save   = lupri
      outfile_node = 6
      lupri        = outfile_node
!
      IF (LUCI_MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (LUCIFILN,'(A20,A1,I1)') LUCITABASF,'.',LUCI_MYPROC
         LUFIL=22
      ELSE IF (LUCI_MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (LUCIFILN,'(A20,A1,I2)') LUCITABASF,'.',LUCI_MYPROC
         LUFIL=23
      ELSE IF (LUCI_MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (LUCIFILN,'(A20,A1,I3)') LUCITABASF,'.',LUCI_MYPROC
         LUFIL=24
      ELSE
         call quit("luci_nmproc.gt.1000!                                &
                    extend the lucita_coworker_main module")
      ENDIF
!
!     open the local input file and the node specific output file.
!     every access to the local stdout handle then automatically writes
!     to the corresponding output file.
      open(outfile_node,file=lucifiln(1:lufil))

!     set marker on incoming work space 
      call memget('REAL',k1,1,work_dalton,kfree,lfree)

      print_lvl = 0
      call setup_lucita_inc_wrkspc_alloc_cw(lucita_ci_run_id,           &
                                            work_dalton,                &
                                            kfree,                      &
                                            lfree,                      &
                                            kcref,                      &
                                            khc,                        &
                                            kresolution_mat,            &
                                            kint1_or_rho1,              &
                                            kint2_or_rho2,              &
                                            len_cref_mc2lu,             &
                                            len_hc_mc2lu,               &
                                            len_resolution_mat_mc2lu,   &
                                            len_int1_or_rho1_mc2lu,     &
                                            len_int2_or_rho2_mc2lu,     &
                                            print_lvl)

!     enter the LUCITA driver -- joining with the master
!     --------------------------------------------------
      call lucita_driver(work_dalton(kfree),lfree,                      &
                         work_dalton(kcref),                            &
                         work_dalton(khc),                              &
                         work_dalton(kresolution_mat),                  &
                         work_dalton(kint1_or_rho1),                    &
                         work_dalton(kint2_or_rho2))

!     release marker on incoming work space
      call memrel('lucita.done',work_dalton,kfrsav,kfrsav,kfree,lfree)

      close(outfile_node,status='KEEP')

      call qexit('lucita_coworker_main')

      END
#endif    /* ifdef VAR_MPI */
