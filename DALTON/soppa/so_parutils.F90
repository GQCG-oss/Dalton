#ifdef VAR_MPI
module so_parutils
!
   use so_info, only: sop_dp
! SOPPA parallel/mpi utilities
!
! This module defines some parameters and some subroutines, which are
! used for distributing work for SOPPA calculations over mpi.
!
! The common sin, sadly a lot of common-blocks still require this
#include "implicit.h"
!
! Various other common-blocks that set occasionally used parameters
#include "mpif.h"
#include "maxorb.h"
#include "maxash.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
#include "iprtyp.h"
#include "maxaqn.h"
#include "chrnos.h"
#include "ibtpar.h"
!
! These common-blocks are needed across this module
!
#include "infpar.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "soppinf.h"
! The module should be used with "only:", if this anyone wants to
! use it differently, we should make private the default (so as not to export
! common blocks)
!   private
!
! Maybe get this from iso_fortran_env, when supported in most compilers?
!

   integer, parameter :: real8 = sop_dp !kind(1.0D0)
!  Flags to be send to slaves to tell them, what work to do
   integer, parameter :: parsoppa_release_slave = 0,  &! Leave
     &                   parsoppa_do_eres  = 1,       &! Call eres routine
     &                   parsoppa_update_amplitudes = 2
!
! SOPPA communicator (needed if not all nodes participate in soppa)
   integer(mpi_integer_kind) :: soppa_comm_active    ! communicator
   integer                   :: soppa_num_active = 0 ! number of nodes

! In order to work with both ERI and HERMITE direct...
   integer            :: soppa_nint   ! number of calls to integral program
!
! Make the defines in infpar a fortran parameter (nicer IMO).
! also, maybe change the name?
   integer(mpi_integer_kind), parameter :: my_mpi_integer = my_MPI_INTEGER
   integer(mpi_integer_kind), parameter :: my_mpi_logical = my_MPI_LOGICAL
   ! Actually every call to MPI functions should be explicitly typed,
   ! but it is a pain to write 1_mpi_integer_kind everywhere...
   integer(mpi_integer_kind), parameter :: one_mpi = 1, zero_mpi = 0
   integer(mpi_integer_kind), parameter :: sop_master = 0
#undef my_MPI_INTEGER

   private soppa_update_common, stupid_isao_bcast_routine

contains


   subroutine soppa_update_common()
!
!  Subroutine that broad-casts various common-blocks from master to slaves
!  in parallel soppa calculations.
!  This routine is defined in order to take these quite verbose
!  calls out of the main flow of the soppa_nodedriver/par_so_eres
!  routines.
!
!  It must be called by the master and all slaves or not at all.
!
!  Rasmus Faber 13/7 - 2015
!  Using code of F. Beyer
!
!  Eventually we may want to replace as much of the broad-casting of easily
!  recalculatable information with calls to the proper initiation routines
!
!Import all the common-blocks, which we don't want have
!poluting the name spaces of routines, that import the module
!
      use dyn_iadrpk

#include "priunit.h"
! These include files depend on previous include files
!#include "infpar.h"
#include "eribuf.h"
#include "eritap.h"

! The rest...
#include "inftap.h"
!#include "ccorb.h"
#include "infind.h"
#include "blocks.h"
#include "ccsdinp.h"
!#include "ccsdsym.h"
#include "ccsdio.h"
#include "distcl.h"
#include "cbieri.h"
! Parameters + stuff initialized in so_init, depends on ccorb.h
!#include "soppinf.h"
#include "aobtch.h"
#include "odclss.h"
#include "ccom.h"
#include "ericom.h"
#include "eridst.h"
#include "erithr.h"
#include "erimem.h"
#include "odbtch.h"
#include "nuclei.h"
#include "symmet.h"
#include "r12int.h"
#include "hertop.h"
#include "cbirea.h"
#include "erisel.h"
#include "symsq.h"
#include "gnrinf.h"
#include "ccpack.h"
#include "ccinftap.h"

      integer(mpi_integer_kind) :: bytesize, ierr, count_mpi 
         !
         !The infinite list of getbytespan -- mpi_bcast starts here
         !
      call getbytespan(lbuf, eribufLAST, bytesize)
      call mpi_bcast(lbuf, bytesize, mpi_byte, sop_master, mpi_comm_world, ierr)

      call getbytespan(luaorc, eritapLAST, bytesize)
      call mpi_bcast(luaorc, bytesize, mpi_byte, sop_master, mpi_comm_world, ierr)

      call getbytespan(nsym, ccorbLAST, bytesize)
      call mpi_bcast(nsym, bytesize, mpi_byte, sop_master, mpi_comm_world, ierr)

      call getbytespan(irow, infindLAST, bytesize)
      call mpi_bcast(irow, bytesize, mpi_byte, sop_master, mpi_comm_world, ierr)

      call getbytespan(centsh, blocksLAST, bytesize)
      call mpi_bcast(centsh, bytesize, mpi_byte,sop_master, mpi_comm_world, ierr)

      call getbytespan(skip, ccsdgninpLAST, bytesize)
      call mpi_bcast(skip, bytesize, mpi_byte, sop_master, mpi_comm_world, ierr)

      call getbytespan(ccs, ccmodelsLAST, bytesize)
      call mpi_bcast(ccs, bytesize, mpi_byte, sop_master, mpi_comm_world, ierr)

      call getbytespan(etmp, etmpLAST, bytesize)
      call mpi_bcast(etmp, bytesize, mpi_byte, sop_master, mpi_comm_world, ierr)

      call getbytespan(nckijmax, ccsdmaxLAST, bytesize)
      call mpi_bcast(nckijmax, bytesize, mpi_byte,sop_master,mpi_comm_world,ierr)

      call getbytespan(nt1amx, ccsdsymLAST, bytesize)
      call mpi_bcast(nt1amx, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call getbytespan(it2del, ccsdioLAST, bytesize)
      call mpi_bcast(it2del, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call getbytespan(mxcall, distclLAST, bytesize)
      call mpi_bcast(mxcall, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call getbytespan(runeri, cbieriLAST, bytesize)
      call mpi_bcast(runeri, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call getbytespan(sotime, soppinfLAST, bytesize)
      call mpi_bcast(sotime, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call getbytespan(nexci2, soppexcLAST, bytesize)
      call mpi_bcast(nexci2, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call getbytespan(soorwc, rwinfLAST, bytesize)
      call mpi_bcast(soorwc, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call getbytespan(expbt, aobtchLAST, bytesize)
      call mpi_bcast(expbt, bytesize, mpi_byte, sop_master, mpi_comm_world, ierr)

      call getbytespan(nitcl, odclssLAST, bytesize)
      call mpi_bcast(nitcl, bytesize, mpi_byte, sop_master, mpi_comm_world, ierr)

      call getbytespan(thrs, ccomLAST, bytesize)
      call mpi_bcast(thrs, bytesize, mpi_byte, sop_master, mpi_comm_world, ierr)

      call getbytespan(gtotyp, ccomcLAST, bytesize)
      call mpi_bcast(gtotyp, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call getbytespan(scrmab, ericomLAST, bytesize)
      call mpi_bcast(scrmab, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call getbytespan(thrsh, erithrLAST, bytesize)
      call mpi_bcast(thrsh, bytesize, mpi_byte, sop_master, mpi_comm_world, ierr)

      call getbytespan(memadd, erimemLAST, bytesize)
      call mpi_bcast(memadd, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call getbytespan(nodbch, odbtchLAST, bytesize)
      call mpi_bcast(nodbch, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call getbytespan(charge, nucleiLAST, bytesize)
      call mpi_bcast(charge, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call getbytespan(ndistr, eridstLAST, bytesize)
      call mpi_bcast(ndistr, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call getbytespan(namn, nuclecLAST, bytesize)
      call mpi_bcast(namn, bytesize, mpi_byte, sop_master, mpi_comm_world, ierr)

      call getbytespan(maxrep, symmtiLAST, bytesize)
      call mpi_bcast(maxrep, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call getbytespan(fmult, symmtrLAST, bytesize)
      call mpi_bcast(fmult, bytesize, mpi_byte, sop_master, mpi_comm_world, ierr)

      call getbytespan(gamac, comr12LAST, bytesize)
      call mpi_bcast(gamac, bytesize, mpi_byte, sop_master, mpi_comm_world, ierr)

      call getbytespan(mbas1, cmmmulLAST, bytesize)
      call mpi_bcast(mbas1, bytesize, mpi_byte, sop_master, mpi_comm_world, ierr)

      call getbytespan(jtop, hertopLAST, bytesize)
      call mpi_bcast(jtop, bytesize, mpi_byte, sop_master, mpi_comm_world, ierr)

      call getbytespan(zcmval, cbireaLAST, bytesize)
      call mpi_bcast(zcmval, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call getbytespan(mulnam, cbirea_cLAST, bytesize)
      call mpi_bcast(mulnam, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call getbytespan(nmulbs, cmmbasLAST, bytesize)
      call mpi_bcast(nmulbs, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call getbytespan(i2bst, symsqLAST, bytesize)
      call mpi_bcast(i2bst, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call getbytespan(thrpckint, ccpackLAST, bytesize)
      call mpi_bcast(thrpckint,bytesize,mpi_byte,sop_master,mpi_comm_world, ierr)

      call getbytespan(luiajb, cc_tapLAST, bytesize)
      call mpi_bcast(luiajb,bytesize,mpi_byte,sop_master,mpi_comm_world, ierr)

      call getbytespan(gradml, gnrinfLAST, bytesize)
      call mpi_bcast(gradml, bytesize, mpi_byte, sop_master, mpi_comm_world,ierr)

      call stupid_isao_bcast_routine()

      ! Common blocks with explicitly set sizes, no need to calculate bytespan
      ! Though we should probably still insert the actual parameter governing
      ! their length, eg. lbasdir+12 = 612
      count_mpi = lbasdir+12
      call mpi_bcast(basdir, count_mpi, mpi_character, sop_master,  &
     &               mpi_comm_world, ierr)
      call mpi_bcast(fnvajkl, 10_mpi_integer_kind, mpi_character, sop_master, mpi_comm_world,ierr)
      call mpi_bcast(vclthr, 24_mpi_integer_kind, mpi_byte, sop_master, mpi_comm_world, ierr)

      ! The slaves need to create the iadrpk array by
      ! calling the module functions rather than with a bcast
      ! iadrpk_dim is initialized in get_iadrpk, no need to send it
!      call mpi_bcast(iadrpk_dim, 1, my_mpi_integer, 0,mpi_comm_world,ierr)
      if (.not. allocated(iadrpk) ) then
         call get_iadrpk(lupri, nsym, muld2h, nbas,           &
     &                   nbast, i2bst, iaodis, iaodpk)
      endif

      if (herdir) then
         soppa_nint = maxshl
      else ! ERIDI
         soppa_nint = mxcall
      endif

!      CALL DZERO(WORK(KAIJ),LAIJ)
!      CALL DZERO(WORK(KAAB),LAAB)

      return

   end subroutine soppa_update_common



   subroutine soppa_nodedriver (work, lwork, iprint)
!
! This is an internal nodedriver routine for parallel AO-based
! SOPPA calculations. The slaves go here from dalton_nodriver, when
! they are activated for SOPPA work.
! Upon entry they are updated with at lot of stuff, which is needed
! across SOPPA routines.
! The slaves then wait to be told to call a particular SOPPA routine.
! (Currently only SO_ERES)
!
! Once the SOPPA part of the calculation is over they can be released
! back to dalton_nodedriver
!
! Rasmus Faber 13/7 - 2015
!
!#include "parsoppa.h"
!
      use so_info, only: sop_models
      implicit none
!
      integer, intent(in) :: iprint, lwork
      real(real8), intent(inout) :: work(*)

      logical :: update_common_blocks
      integer :: soppa_work_kind, maxnumjobs
! Use proper integerkind for mpi-arguments, irrespective of
! what flags are use for compilation of MPI and Dalton
! Not that it will change things right away, but perhaps one day
      integer(mpi_integer_kind) :: ierr, numprocs, mycolor, count_mpi, mynum_mpi
! Lengths of arrays
      integer :: lt2am, lfockd, ldensij, ldensab, ldensai, lworkf, &
                 LAssignedIndices
! Pointers to arrays
      integer :: kt2am, kfockd, kdensij, kdensab, kdensai, kend, &
                 kAssignedIndices
! Other integers
      character(len=5) :: model
! Some info, that we need in each pass
      integer :: nnewtr, noldtr, isymtr, nit, idtype
      integer :: imod, imodel
! Need to ensure that the six above variables are stored
! consecutively, so we can recieve them with a single bcast.
! This is the only purpose of info_array, only address it as
! part of communication!
      integer :: info_array(6)
      equivalence (info_array(1), isymtr), (info_array(2), nit), &
                  (info_array(3), nnewtr), (info_array(4), noldtr), &
                  (info_array(5), idtype), (info_array(6), imodel)
      !
      ! Recieve the method on which to work
      !
      call mpi_bcast( imod, one_mpi, my_mpi_integer, sop_master, mpi_comm_world, &
                      ierr )
      ! Do we need to do an allreduce to check that all
      ! Processes agree not to update the common-blocks?
      ! For now just let master tell the slaves
      call mpi_bcast(update_common_blocks, one_mpi, my_mpi_logical,     &
                     sop_master, mpi_comm_world, ierr)

      if ( update_common_blocks ) then
         ! Recieve all the system and symmetry information
         ! from the master here.
         ! Master MUST also call this routine
         call soppa_update_common()
      endif
      !  Set print-level for slave
      iprsop = iprint
      !
      ! Set up communicator for soppa. Work will later be separated into
      ! integral distributions, which will be distributed among the nodes
      ! If there are more processes than integral distributions,
      ! we don't include the rest in the soppa communicator.
      !
      !numprocs = nodtot + 1 ! nodtot from infpar.h
      call mpi_comm_size( mpi_comm_world, numprocs, ierr)
      if ( numprocs .le. soppa_nint ) then
         !
         ! Usual case:
         ! Just keep mpi_comm_world
         soppa_num_active  = numprocs
         soppa_comm_active = mpi_comm_world
      else
         ! This should happen so rarely, that we don't really need
         ! to worry about it. Should be done in some intelligent
         ! manner though
         soppa_num_active  = soppa_nint
         if ( mynum .ge. soppa_num_active ) then
            mycolor = MPI_UNDEFINED
         else
            mycolor = 0
         endif
         mynum_mpi = mynum
         call mpi_comm_split( mpi_comm_world, mycolor, mynum_mpi, &
                              soppa_comm_active, ierr)
      endif

      !
      ! Allocate the work-array
      !
      ! This is allocations, which persist across iterations
      ! We need at least DENSIJ, DENSAB, FOCKD, DENSAI, and
      ! the mp2 - amplitudes
      !
      kend = 1 ! Start from one, or some suitably alligned location?

      !
      ! These allocations mirror exactly those of so_excit1
      !
      lfockd = norbt   ! from ccorb.h
      kfockd = kend
      kend   = kfockd + lfockd

      ! following only if not RPA...
      if ( imod .gt. 0 ) then
         lt2am   = nt2amx             ! from ccsdsym.h
         ldensij = nijden(1)          ! from soppinf.h
         ldensab = nabden(1)          ! from soppinf.h
         ldensai = naiden(1)          ! from soppinf.h

         kt2am   = kend
         kdensij = kt2am   + lt2am    ! densij and densab is not
         kdensab = kdensij + ldensij  ! currently used on slaves
         kdensai = kdensab + ldensab
         kend    = kdensai + ldensai
         !
         ! Zero densai (To mirror initialization in so_excit1)
         call dzero( work(kdensai), ldensai )
      else
         ! For RPA initialize the addresses as too large, to ensure a crash
         ! if they are for some reason accessed anyway
         kt2am   = huge(lwork)
         kdensij = huge(lwork)
         kdensab = huge(lwork)
         kdensai = huge(lwork)

         lt2am   = 0
         ldensij = 0
         ldensab = 0
         ldensai = 0
      endif

!
!     Allocation of of space for load-balancing
!
      maxnumjobs = soppa_nint - min(soppa_nint, numprocs) + 1
      lAssignedIndices = (maxnumjobs + 1) / irat
      kAssignedIndices = kend
      kend = kAssignedIndices + lAssignedIndices

      lworkf = lwork - kend
      if ( kend .gt. lwork ) call stopit('SOPPA_NODEDRIVER', '2',   &
     &                                    kend, lwork )

!
!     Bugfix: need to call er2ini here, so that it does not overwrite
!             configured values later
      call er2ini

!
! Go to an infinite loop... While we are here, the master
! broadcasts job descriptions to the slaves.
!
      do
            ! Recieve work from master
         call mpi_bcast ( soppa_work_kind, one_mpi, my_mpi_integer, sop_master,       &
     &               mpi_comm_world, ierr )

            ! Act according to the job type recieved
         select case ( soppa_work_kind )
            !
         case (parsoppa_release_slave)
            !
            ! Slaves no longer needed in parallel soppa
            !-------------------------------------------
            ! Free any SOPPA communicators
            if ( (soppa_comm_active .ne.  mpi_comm_world).and.         &
                 (mynum .lt. soppa_num_active)                ) then
               call mpi_comm_free( soppa_comm_active, ierr )
            endif
            soppa_num_active = 0
            !
            ! Return to dalton_nodedriver
            return
            !
         case (parsoppa_do_eres)
            !
            ! Calculate linear tranformation of trial-vectors
            !-------------------------------------------------
            !
            ! We need to communicate NOLDTR, NNEWTR, ISYMTR and NIT.
            ! These are joined together in the info_array.
            ! MODEL has allready have been communicated
            !
            call mpi_bcast( info_array(1), 6_mpi_integer_kind, my_mpi_integer, &
                            sop_master, mpi_comm_world, ierr)

            ! Inactive processes do nothing
            if ( soppa_comm_active .eq. MPI_COMM_NULL ) cycle
            model = sop_models(imodel)
            !
            ! All should now run through same ERES routine
            !
            call so_eres( model, noldtr, nnewtr,          &! General info
                  work(kdensij), ldensij,                 &! Densij
                  work(kdensab), ldensab,                 &! Densab
                  work(kt2am),lt2am,                      &! t2 amplitudes
                  work(kfockd), lfockd,                   &! Fockd
                  work(kdensai), ldensai,                 &! Densai
                  nit, isymtr,                            &! Info
                  idtype,                                 &! dynamic or static?
                  work(kassignedindices),maxnumjobs,      &! Load-balancing space
     &            work(kend), lworkf )                     ! Work-array
            !

         case (parsoppa_update_amplitudes)
            !
            !  The master have read/calculated new amplitudes
            !
            if (imod.eq.0) call quit('Error in AOSOPPA nodedriver:'// &
               ' Slave have not reserved memory for amplitudes')

            count_mpi = lt2am
            call mpi_bcast( work(kt2am), count_mpi, mpi_real8, sop_master,        &
     &                      mpi_comm_world, ierr )

         case default
            !
            ! Anything else
            !--------------
            call quit('Slave recieved invalid job-description'//     &
                            ' in AOSOPPA nodedriver.' )
         endselect

      enddo

      ! We should never reach here...

   end subroutine soppa_nodedriver

   subroutine soppa_initialize_slaves( update_common_blocks, rpa_only )
!    -----------------------------------------------------------
!     This subroutine tells the slaves that hang in
!     dalton_nodedriver to enter the soppa node-driver and sends
!     information, with doesn't change between soppa iterations.
!    -----------------------------------------------------------
      implicit none
#include "iprtyp.h"
#include "distcl.h"
!
! Arguments
      logical, intent(in)        :: update_common_blocks, rpa_only
!
! Locals
      integer(mpi_integer_kind)  :: ierr
      integer                    :: numprocs
      integer                    :: imodel

      if (rpa_only) then
         imodel = 0
      else
         imodel = 1
      end if
!      if (nodtot .eq. 0 ) return
         !
         ! Send the slave to the soppa_nodedriver
      call mpixbcast( PARA_SO_ERES, 1, 'INTEGE', 0)
         ! Set slave print-level
      call mpixbcast( 0, 1, 'INTEGE', 0)
         !
         ! Here the slaves enter soppa_nodedriver
         !
         ! Set the method
      call mpi_bcast( imodel, one_mpi, my_mpi_integer, sop_master, mpi_comm_world, &
                      ierr )
         !
         ! Send the various common blocks if needed
      call mpi_bcast( update_common_blocks, one_mpi, my_mpi_logical,    &
     &                   sop_master, mpi_comm_world, ierr )
!
      if ( update_common_blocks ) then

         call soppa_update_common()
      endif
      !
      !  Setup communicator for SOPPA, see comment in soppa_nodedriver
      !
      numprocs = nodtot + 1 ! from infpar.h
      if ( numprocs .le. soppa_nint ) then ! mxcall from distcl.h
         soppa_num_active  = numprocs
         soppa_comm_active = mpi_comm_world
      else
         soppa_num_active  = soppa_nint
         call mpi_comm_split( mpi_comm_world, zero_mpi, zero_mpi,  &
                              soppa_comm_active, ierr)
      endif
      !
      return

   end subroutine soppa_initialize_slaves

   subroutine soppa_update_amplitudes(t2mp, lt2mp)
!
! Send the amplitudes to the slaves. The slaves must be
! waiting in soppa_nodedriver when this is called.
!
      implicit none
      integer, intent(in) :: lt2mp
      real(sop_dp), intent(in) :: t2mp(lt2mp)
!
      integer(mpi_integer_kind) :: ierr, count_mpi
!
! Tell slaves that the amplitudes will be send
      call mpixbcast(parsoppa_update_amplitudes, 1, 'INTEGE', 0)
!
! Send the actual amplitudes
      count_mpi = lt2mp
      call mpi_bcast(t2mp, count_mpi, mpi_real8, sop_master, mpi_comm_world, ierr)
      return
   end subroutine soppa_update_amplitudes

   subroutine soppa_release_slaves()
!
! Release the slaves and send them back to the main
! node-driver. This routine should be called my the master,
! when we no longer need the slaves for soppa work
!
! The slaves must be in the SOPPA node-driver, when this routine is
! called
      implicit none
      integer(mpi_integer_kind) :: ierr
!
! Send release signal
      call mpixbcast(parsoppa_release_slave, 1, 'INTEGE', 0)
!
! Deallocate the communicator
!
      if (soppa_comm_active .ne. mpi_comm_world )  then
!         write (*,*) 'Freeing soppa comm'
!         write (*,*) soppa_comm_active, mpi_comm_world
         call mpi_comm_free( soppa_comm_active , ierr )
      endif
!
      return
   end subroutine soppa_release_slaves

   subroutine stupid_isao_bcast_routine()
! For the isoa
#include "ccisao.h"
!
!  This is a simple (but stupid) work-around to the problem that
!  the array ISAO exist in both infind.h and ccisao.h and that it seem
!  pretty random which version is needed by which routine
!  (the content appear to be the same).
!  A better solution would maybe be to remove isao from common block
!  in infind.h and explicitly include ccisao.h there...
!  Or even just replace ISAO with the function
!
!  pure function isao ( aoindex ) result(isym)
!  #include ! something that sets nbas(8), nsym
!  integer, intent(in) :: aoindex
!  integer :: i, isym, popsum
!  popsum = 0
!  do isym = 1, nsym
!     popsum = nbas(isym) + popsum
!     if ( aoindex .le. popsum ) return
!  enddo
!  ! Some error statement, here would be in order
!  end function isao
!
      integer(mpi_integer_kind)  :: bytesize, ierr

!    This size SHOULD be MXCORB x "integer size"
!    But since this is a cc block, someone may change that to MXCORB_CC
      call getbytespan(isao, ccisaolast, bytesize )
      call mpi_bcast(isao, bytesize, mpi_byte, sop_master, mpi_comm_world,   &
     &               ierr)

   end subroutine stupid_isao_bcast_routine


endmodule so_parutils
!
! This routine, though currently only used insides the module, sits here
! because fortran only recently got a void type (TYPE(*)), and it will be
! a pain to write overloaded interfaces for all possible combinations.
! /* deck getbytespan */
subroutine getbytespan(firstvar, lastvar, bytespan)
! Frederik Beyer, March 2014.
!
! This subroutine calculates the memory span in bytes between two variables.
!
! This is used for easy updating of common block in parallel calculations.
! The former approaches relied on counting the number of occurrences
! of variables of a particular type and then transferring a common in
! a series of mpi_bcasts; one for every datatype in the common block.
! This approach utilizes the contiguous storage of variables in a common
! block to calculate the span in bytes from the first variable
! in the common block and up to, but not including, the last variable.
!
! Point the firstvar to the first variable in the common block
! and point lastvar to the last variable in the common block.
!
! Lastvar should have a name like <commonblockname>LAST
! ie. CCSDGNINPLAST (see include/ccsdinp.h for an example)
!
! These <name>LAST variables are only there to facilitate easy common block
! transfers. They are never explicitly needed in a calculation.
! They should always be of type int.
!
! Getbytespan calculates the total amount of bytes needed for an MPI_BCAST
! to transfer the whole block in one go (with datatype = mpi_byte),
! including the first but excluding the <name>last variable.
!
! Example of use:
! to update the common block /eribuf/:
!
!      call getbytespan(lbuf, eribufLAST, bytesize)
!      call mpi_bcast(lbuf, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      implicit none
#include "mpif.h"
#include "priunit.h"

      integer(mpi_integer_kind), intent(out) :: bytespan
      integer,intent(in)   :: firstvar, lastvar
!RF This should be mpi_address_kind
      integer(mpi_address_kind) :: firstmem, lastmem
      integer(mpi_integer_kind) :: ierr
      integer              :: totaltransfer = 0
      integer, parameter   :: approxeager = 4096 ! Approximate upper limit for eager protocol transfers. Implementation dependent.
      integer              :: numeagersends=0, numrendezsends=0
      logical, parameter   :: debug = .false.


      call mpi_get_address(firstvar, firstmem, ierr)
      call mpi_get_address(lastvar, lastmem, ierr)
      bytespan = lastmem - firstmem
      if (bytespan.lt.1) then
          call quit('subroutine getbytespan calculated a non-sensical ', &
     &    'size for common block transfer.')
      endif


      if (debug) then
         totaltransfer = bytespan + totaltransfer
         if (bytespan.lt.approxeager) then
            numeagersends   = numeagersends + 1
         else
            numrendezsends  = numrendezsends + 1
         endif

         write(lupri, '(a, i8)')  &
     &   "Running amount of Bytes transferred via getbytespan: "  &
     &   ,totaltransfer
         write(lupri, '(a, i8)') "Estimated transfers using the ", &
     &               "eager protocol: ", numeagersends
         write(lupri, '(a, i8)') "Estimated transfers using the ", &
     &                "rendezvous protocol: ", numrendezsends
      endif

      return
end subroutine

#else
! Dummy subroutune
subroutine so_parutils()
   call quit('MPI dummy subroutine called')
endsubroutine
#endif
