#ifdef VAR_MPI
module so_parutils
!
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

   integer, parameter :: real8 = kind(1.0D0)
!  Flags to be send to slaves to tell them, what work to do
   integer, parameter :: parsoppa_release_slave = 0,   &! Leave
     &                   parsoppa_do_eres  = 1  ! Call eres routine
!
! SOPPA communicator (needed if not all nodes participate in soppa)
   integer(mpi_integer_kind) :: soppa_comm_active    ! communicator
   integer                   :: soppa_num_active     ! number of nodes

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
!#include "parsoppa.h"
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

      integer(mpi_integer_kind) :: bytesize, ierr
         !
         !The infinite list of getbytespan -- mpi_bcast starts here
         !
      call getbytespan(lbuf, eribufLAST, bytesize)
      call mpi_bcast(lbuf, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(luaorc, eritapLAST, bytesize)
      call mpi_bcast(luaorc, bytesize, mpi_byte,0, mpi_comm_world, ierr)

      call getbytespan(nsym, ccorbLAST, bytesize)
      call mpi_bcast(nsym, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(irow, infindLAST, bytesize)
      call mpi_bcast(irow, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(centsh, blocksLAST, bytesize)
      call mpi_bcast(centsh, bytesize, mpi_byte,0, mpi_comm_world, ierr)

      call getbytespan(skip, ccsdgninpLAST, bytesize)
      call mpi_bcast(skip, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(ccs, ccmodelsLAST, bytesize)
      call mpi_bcast(ccs, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(etmp, etmpLAST, bytesize)
      call mpi_bcast(etmp, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(nckijmax, ccsdmaxLAST, bytesize)
      call mpi_bcast(nckijmax, bytesize, mpi_byte,0,mpi_comm_world,ierr)

      call getbytespan(nt1amx, ccsdsymLAST, bytesize)
      call mpi_bcast(nt1amx, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(it2del, ccsdioLAST, bytesize)
      call mpi_bcast(it2del, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(mxcall, distclLAST, bytesize)
      call mpi_bcast(mxcall, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(runeri, cbieriLAST, bytesize)
      call mpi_bcast(runeri, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(sotime, soppinfLAST, bytesize)
      call mpi_bcast(sotime, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(nexci2, soppexcLAST, bytesize)
      call mpi_bcast(nexci2, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(soorwc, rwinfLAST, bytesize)
      call mpi_bcast(soorwc, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(expbt, aobtchLAST, bytesize)
      call mpi_bcast(expbt, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(nitcl, odclssLAST, bytesize)
      call mpi_bcast(nitcl, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(thrs, ccomLAST, bytesize)
      call mpi_bcast(thrs, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(gtotyp, ccomcLAST, bytesize)
      call mpi_bcast(gtotyp, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(scrmab, ericomLAST, bytesize)
      call mpi_bcast(scrmab, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(thrsh, erithrLAST, bytesize)
      call mpi_bcast(thrsh, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(memadd, erimemLAST, bytesize)
      call mpi_bcast(memadd, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(nodbch, odbtchLAST, bytesize)
      call mpi_bcast(nodbch, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(charge, nucleiLAST, bytesize)
      call mpi_bcast(charge, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(ndistr, eridstLAST, bytesize)
      call mpi_bcast(ndistr, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(namn, nuclecLAST, bytesize)
      call mpi_bcast(namn, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(maxrep, symmtiLAST, bytesize)
      call mpi_bcast(maxrep, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(fmult, symmtrLAST, bytesize)
      call mpi_bcast(fmult, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(gamac, comr12LAST, bytesize)
      call mpi_bcast(gamac, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(mbas1, cmmmulLAST, bytesize)
      call mpi_bcast(mbas1, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(jtop, hertopLAST, bytesize)
      call mpi_bcast(jtop, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(zcmval, cbireaLAST, bytesize)
      call mpi_bcast(zcmval, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(mulnam, cbirea_cLAST, bytesize)
      call mpi_bcast(mulnam, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(nmulbs, cmmbasLAST, bytesize)
      call mpi_bcast(nmulbs, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(i2bst, symsqLAST, bytesize)
      call mpi_bcast(i2bst, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(thrpckint, ccpackLAST, bytesize)
      call mpi_bcast(thrpckint,bytesize,mpi_byte,0,mpi_comm_world, ierr)

      call getbytespan(luiajb, cc_tapLAST, bytesize)
      call mpi_bcast(luiajb,bytesize,mpi_byte,0,mpi_comm_world, ierr)

      call getbytespan(gradml, gnrinfLAST, bytesize)
      call mpi_bcast(gradml, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call stupid_isao_bcast_routine()

      ! Common blocks with explicitly set sizes, no need to calculate bytespan
      ! Though we should probably still insert the actual parameter governing
      ! their length, eg. lbasdir+12 = 612
      call mpi_bcast(basdir,lbasdir+12, mpi_character, 0,  &
     &               mpi_comm_world, ierr)
      call mpi_bcast(fnvajkl, 10, mpi_character, 0, mpi_comm_world,ierr)
      call mpi_bcast(vclthr, 24, mpi_byte, 0, mpi_comm_world, ierr)

      ! The slaves need to create the iadrpk array by
      ! calling the module functions rather than with a bcast
      ! iadrpk_dim is initialized in get_iadrpk, no need to send it
!      call mpi_bcast(iadrpk_dim, 1, mpi_integer, 0,mpi_comm_world,ierr)
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
! Need MXCALL
!#include "distcl.h"
!
      integer, intent(in) :: iprint, lwork
      real(real8)         :: work(*)
      logical :: update_common_blocks
      integer :: soppa_work_kind
! Use proper integerkind for mpi-arguments, irrespective of
! what flags are use for compilation of MPI and Dalton
! Not that it will change things right away, but perhaps one day
      integer(mpi_integer_kind) :: ierr, numprocs, mycolor
! Lengths of arrays
      integer :: lt2am, lfockd, ldensij, ldensab, ldensai
! Pointers to arrays
      integer :: kt2am, kfockd, kdensij, kdensab, kdensai
! Other integers
      character(len=5) :: model
! Some info, that we need in each pass
      integer :: nnewtr, noldtr, isymtr, nit
! Need to ensure that the four above variables are stored
! consecutively, so we can recieve them with a single bcast.
! This is the only purpose of info_array, only address it as
! part of communication!
      integer :: info_array(4)
      equivalence (info_array(1), isymtr), (info_array(2), nit), &
                  (info_array(3), nnewtr), (info_array(4), noldtr)
      !
      ! Recieve the method on which to work
      !
      call mpi_bcast( model, 5, mpi_character, 0, mpi_comm_world, &
                      ierr )
      ! Do we need to do an allreduce to check that all
      ! Processes agree not to update the common-blocks?
      ! For now just let master tell the slaves
      call mpi_bcast(update_common_blocks, 1, my_mpi_logical,     &
                     0, mpi_comm_world, ierr)

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
      numprocs = nodtot + 1 ! nodtot from infpar.h
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
         call mpi_comm_split( mpi_comm_world, mycolor, mynum, &
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
      if ( model .ne. 'AORPA' ) then
         lt2am   = nt2amx             ! from ccsdsym.h
         if (triplet) lt2am = 2*nt2amx
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
         ! For RPA initialize the addresses as negative, to ensure a crash
         ! if they are for some reason accessed anyway
         kt2am   = -1
         kdensij = -1
         kdendab = -1
         kdensai = -1
      endif

!
!     Allocation of of space for load-balancing
!
      call mpi_comm_size( mpi_comm_world, numprocs, ierr)
      maxnumjobs = soppa_nint - min(soppa_nint, numprocs) + 1
      lAssignedIndices = (maxnumjobs + 1) / irat
      kAssignedIndices = kend
      kend = kAssignedIndices + lAssignedIndices

      lworkf = lwork - kend
      if ( kend .gt. lwork ) call stopit('SOPPA_NODEDRIVER', '2',   &
     &                                    kend, lwork )

      ! Recieve the stuff, which is already known.
      if ( model .ne. 'AORPA' ) then
         ! The MP2 amplitudes
         call mpi_bcast( work(kt2am), lt2am, mpi_real8, 0,        &
     &                   mpi_comm_world, ierr )
         ! Densab and Densij could be added here, but is currently
         ! not used by the slaves...
      endif
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
         call mpi_bcast ( soppa_work_kind, 1, my_mpi_integer, 0,       &
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
            call mpi_bcast( info_array(1), 4, my_mpi_integer, 0,     &
                            mpi_comm_world, ierr)

            ! Inactive processes do nothing
            if ( soppa_comm_active .eq. MPI_COMM_NULL ) cycle
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
                  work(kassignedindices),maxnumjobs,      &! Load-balancing space
     &            work(kend), lworkf )                     ! Work-array
            !
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

   subroutine soppa_initialize_slaves( update_common_blocks,         &
     &             t2mp, lt2mp, model )
!    -----------------------------------------------------------
!     This subroutine tells the slaves that hang in
!     dalton_nodedriver to enter the soppa node-driver and sends
!     information, with doesn't change between soppa iterations.
!    -----------------------------------------------------------
#include "iprtyp.h"
#include "distcl.h"
!
! Arguments
      logical, intent(in)        :: update_common_blocks
      integer, intent(in)        :: lt2mp
      real(real8)                :: t2mp(lt2mp)
      character(len=5),intent(in):: model
!
! Locals
      integer(mpi_integer_kind)  :: ierr
      integer                    :: numprocs

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
      call mpi_bcast( model, 5, mpi_character, 0, mpi_comm_world, &
                      ierr )
         !
         ! Send the various common blocks if needed
      call mpi_bcast( update_common_blocks, 1, my_mpi_logical,    &
     &                   0, mpi_comm_world, ierr )
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
      ! Send allready known stuff such as the
      !
      if ( model .ne. 'AORPA' ) then
         ! The MP2 (or CC T2) amplitudes
         call mpi_bcast(t2mp, lt2mp, mpi_real8, 0,mpi_comm_world,ierr)
         ! Densab and Densai could be added here...
      endif



      return

   end subroutine soppa_initialize_slaves

   subroutine soppa_release_slaves()
!
! Release the slaves and send them back to the main
! node-driver. This routine should be called my the master,
! when we no longer need the slaves for soppa work
!
! The slaves must be in the SOPPA node-driver, when this routine is
! called
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
   endsubroutine soppa_release_slaves

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
      call mpi_bcast(isao, bytesize, mpi_byte, master, mpi_comm_world,   &
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

#include "mpif.h"
#include "priunit.h"

      integer(mpi_integer_kind), intent(out) :: bytespan
      integer,intent(in)   :: firstvar, lastvar
!RF This should be mpi_address_kind
      integer(mpi_address_kind) :: firstmem, lastmem
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

