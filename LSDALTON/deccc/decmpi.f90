!> @file
!> DEC MPI handling
!> \author Kasper Kristensen
!> \date March 2012
#ifdef VAR_LSMPI
module decmpi_module

  use precision
  use typedeftype!, only: lsitem
  use memory_handling!,only: mem_alloc,mem_dealloc
  use dec_typedef_module
  use infpar_module
  use lsmpi_type
  use lsmpi_op,only: mpicopy_lsitem
  use io!, only: io_init
  use Matrix_module!,only: matrix


  ! DEC DEPENDENCIES (within deccc directory)  
  ! *****************************************
  use array_memory_manager!, only: arr_msg_len
  use dec_pdm_module

contains



  !> \brief Send three fragment energies (for occupied, virtual, and Lagrangian schemes)
  !> from given sender (typically a slave) to given receiver (typically the master).
  !> For (T) also pass (T) fragment energies.
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine mpi_send_recv_fragmentenergy(comm,MySender,MyReceiver,fragenergy,job)

    implicit none

    !> Communicator
    integer(kind=ls_mpik),intent(in) :: comm
    !> Rank of sender WITHIN comm group
    integer(kind=ls_mpik),intent(in) :: MySender
    !> Rank of receiver WITHIN comm group
    integer(kind=ls_mpik),intent(in) :: MyReceiver
    !> Fragment energies, see "energy" in ccatom type def
    real(realk),dimension(ndecenergies),intent(inout) :: fragenergy
    !> Fragment job info (job list containing a single job)
    type(joblist),intent(inout) :: job
    integer(kind=ls_mpik) :: mynum,master
    master=0


    ! Sanity check 1: Meaningful rank
    call get_rank_for_comm(comm,mynum)  ! rank number for given communicator
    if( (mynum/=MySender) .and. (mynum/=MyReceiver) ) then
       print '(a,3i6)', 'Rank,sender,receiver',mynum,MySender,MyReceiver
       call lsquit('mpi_send_recv_fragmentenergy: Rank number is neither sender nor receiver!',-1)
    end if

    ! Init buffer
    call ls_mpiInitBuffer(master,LSMPISENDRECV,comm,sender=MySender, receiver=MyReceiver)

    ! Sender: Copy fragment info into buffer
    ! Receiver: Read fragment info from buffer
    call ls_mpi_buffer(fragenergy,ndecenergies,master)
    call mpicopy_fragment_joblist(job)

    ! Finalize buffer
    call ls_mpiFinalizeBuffer(master,LSMPISENDRECV,comm,sender=MySender, receiver=MyReceiver)

  end subroutine mpi_send_recv_fragmentenergy


  !> \brief Send fragment MP2 gradient information from given sender (typically a slave)
  !> to given receiver (typically the master).
  !> The three fragment energies (for occupied, virtual, and Lagrangian schemes) are also communicated.
  !> \author Kasper Kristensen
  !> \date June 2012
  subroutine mpi_send_recv_mp2grad(comm,MySender,MyReceiver,grad,fragenergy,job,only_update)

    implicit none

    !> Communicator
    integer(kind=ls_mpik),intent(in) :: comm
    !> Rank of sender WITHIN comm group
    integer(kind=ls_mpik),intent(in) :: MySender
    !> Rank of receiver WITHIN comm group
    integer(kind=ls_mpik),intent(in) :: MyReceiver
    !> Fragment energies, see "energy" in ccatom type def
    real(realk),dimension(ndecenergies),intent(inout) :: fragenergy
    !> MP2 fragment gradient matrix information
    type(mp2grad),intent(inout) :: grad
    !> Fragment job info (job list containing a single job)
    type(joblist),intent(inout) :: job
    !> Only MPI copy information needed for update_full_mp2gradient
    !> (this will minimize MPI communication for global master)
    logical,intent(in) :: only_update
    integer(kind=ls_mpik) :: mynum,master
    master=0

    ! Sanity check 1: Meaningful rank
    call get_rank_for_comm(comm,mynum)  ! rank number for given communicator
    if( (mynum/=MySender) .and. (mynum/=MyReceiver) ) then
       print '(a,3i6)', 'Rank,sender,receiver',mynum,MySender,MyReceiver
       call lsquit('mpi_send_recv_mp2grad: Rank number is neither sender nor receiver!',-1)
    end if

    ! Init buffer
    call ls_mpiInitBuffer(master,LSMPISENDRECV,comm,sender=MySender, receiver=MyReceiver)

    ! Sender: Copy fragment info into buffer
    ! Receiver: Read fragment info from buffer
    call mpicopy_fragmentgradient(grad,only_update)
    call ls_mpi_buffer(fragenergy,ndecenergies,master)
    call mpicopy_fragment_joblist(job)

    ! Finalize buffer
    call ls_mpiFinalizeBuffer(master,LSMPISENDRECV,comm,sender=MySender, receiver=MyReceiver)

  end subroutine mpi_send_recv_mp2grad



  !> \brief Send fragment from given sender to given receiver.
  !> NOTE: Only the information OUTSIDE the expensive box in the ccatom
  !> type definition can be sent this way. In general, the information outside
  !> the expensive box combined with full molecular information is enough
  !> to detemine the stuff outside the expensive box.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine mpi_send_recv_single_fragment(MyFragment,comm,MySender,MyReceiver,job)

    implicit none

    !> Fragment under consideration
    type(ccatom),intent(inout) :: MyFragment
    !> Communicator
    integer(kind=ls_mpik),intent(in) :: comm
    !> Rank of sender WITHIN comm group
    integer(kind=ls_mpik),intent(in) :: MySender
    !> Rank of receiver WITHIN comm group
    integer(kind=ls_mpik),intent(in) :: MyReceiver
    !> Fragment job info (job list containing a single job)
    type(joblist),intent(inout) :: job
    integer(kind=ls_mpik) :: mynum,master

    master = 0

    ! Sanity check 1: Meaningful rank
    call get_rank_for_comm(comm,mynum)  ! rank number for given communicator
    if( (mynum/=MySender) .and. (mynum/=MyReceiver) ) then
       print '(a,3i6)', 'Rank,sender,receiver',mynum,MySender,MyReceiver
       call lsquit('mpi_send_recv_single_fragment: Rank number is neither sender nor receiver!',-1)
    end if

    ! Sanity check 2: Basis info in "Expensive box" in ccatom should not be sent
    if(MyFragment%BasisInfoIsSet .and. mynum==MySender) then
       call lsquit('mpi_send_recv_single_fragment: Not implemented for &
            & sending/receiving fragment basis info!',-1)
    end if

    ! Init buffer
    call ls_mpiInitBuffer(master,LSMPISENDRECV,comm,sender=MySender, receiver=MyReceiver)

    ! Sender: Copy fragment info into buffer
    ! Receiver: Read fragment info from buffer
    call mpicopy_fragment(MyFragment,comm,.false.)

    ! Copy fragment statistics
    call mpicopy_fragment_joblist(job)

    ! Finalize buffer
    call ls_mpiFinalizeBuffer(master,LSMPISENDRECV,comm,sender=MySender, receiver=MyReceiver)

  end subroutine mpi_send_recv_single_fragment



  !> \brief Send fragments from given sender to given receiver.
  !> NOTE: Only the information OUTSIDE the expensive box in the ccatom
  !> type definition can be sent this way. In general, the information outside
  !> the expensive box combined with full molecular information is enough
  !> to detemine the stuff outside the expensive box.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine mpi_send_recv_many_fragments(natoms,whichfrags,Fragments,comm,MySender,MyReceiver)

    implicit none


    !> Total number of atoms (natoms is smaller than or equal to actual number
    !>                        of fragments to be send/received).
    integer,intent(in) :: natoms
    !> Logical vector telling which fragments should be sent/received
    !> E.g. if whichfrags(3)=true and whichfrags(4)=false then fragment for atom 3
    !> should be sent/received while no fragment for atom 4 should be communicated.
    logical,dimension(natoms),intent(inout) :: whichfrags
    !> Fragments to be sent/received. Only those entries in Fragments specified
    !> by whichfrags will be touched!
    type(ccatom),dimension(natoms),intent(inout) :: Fragments
    !> Communicator
    integer(kind=ls_mpik),intent(in) :: comm
    !> Rank of sender WITHIN comm group
    integer(kind=ls_mpik),intent(in) :: MySender
    !> Rank of receiver WITHIN comm group
    integer(kind=ls_mpik),intent(in) :: MyReceiver
    integer :: atom
    integer(kind=ls_mpik) :: mynum, master
    master=0

    ! Sanity check: Meaningful rank
    call get_rank_for_comm(comm,mynum)  ! rank number for given communicator
    if( (mynum/=MySender) .and. (mynum/=MyReceiver) ) then
       print '(a,3i6)', 'Rank,sender,receiver',mynum,MySender,MyReceiver
       call lsquit('mpi_send_recv_many_fragments: Rank number is neither sender nor receiver!',-1)
    end if

    ! Send logical list of which fragments were done by MySender
    call ls_mpisendrecv(whichfrags,natoms,comm,MySender,MyReceiver)


    ! Init buffer
    call ls_mPiinitbuffer(master,LSMPISENDRECV,comm,sender=MySender, receiver=MyReceiver)

    do atom=1,natoms

       CommunicateFragment: if(whichfrags(atom)) then

          ! Sanity check: Basis info in "Expensive box" in ccatom should not be sent
          if(Fragments(atom)%BasisInfoIsSet) then
             call lsquit('mpi_send_recv_many_fragments: Not implemented for &
                  & sending/receiving fragment basis info!',-1)
          end if

          ! Sender: Copy fragment info into buffer
          ! Receiver: Read fragment info from buffer into Fragments(atom)
          call mpicopy_fragment(Fragments(atom),comm,.false.)

          ! Check that the correct fragment number was communicated
          if(Fragments(atom)%atomic_number /= atom) then
             print '(a,2i8)', 'Counter atom, actual atom ', atom,Fragments(atom)%atomic_number
             call lsquit('mpi_send_recv_many_fragments: Inconsistency between &
                  & counter and actual fragment',-1)
          end if

       end if CommunicateFragment

    end do

    ! Finalize buffer
    call ls_mpiFinalizeBuffer(master,LSMPISENDRECV,comm,&
         & sender=MySender, receiver=MyReceiver)


  end subroutine mpi_send_recv_many_fragments



  !> \brief Bcast fragments from master to slave (within communicator group).
  !> NOTE: Only the information OUTSIDE the expensive box in the ccatom
  !> type definition can be sent this way. In general, the information outside
  !> the expensive box combined with full molecular information is enough
  !> to detemine the stuff outside the expensive box.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine mpi_bcast_many_fragments(natoms,whichfrags,Fragments,comm)

    implicit none


    !> Total number of atoms (natoms is smaller than or equal to actual number
    !>                        of fragments to be bcasted).
    integer,intent(in) :: natoms
    !> Logical vector telling which fragments should be bcasted
    !> E.g. if whichfrags(3)=true and whichfrags(4)=false then fragment for atom 3
    !> should be bcasted while no fragment for atom 4 should be communicated.
    logical,dimension(natoms),intent(in) :: whichfrags
    !> Fragments to be bcasted. Only those entries in Fragments specified
    !> by whichfrags will be touched!
    type(ccatom),dimension(natoms),intent(inout) :: Fragments
    !> Communicator
    integer(kind=ls_mpik),intent(in) :: comm
    integer :: atom
    integer(kind=ls_mpik) :: mynum,master
    master = 0

    call get_rank_for_comm(comm,mynum)  ! rank number for given communicator

    ! Init buffer
    call ls_mpiInitBuffer(master,LSMPIBROADCAST,comm)

    do atom=1,natoms

       CommunicateFragment: if(whichfrags(atom)) then

          ! Sanity check: Basis info in "Expensive box" in ccatom should not be sent
          if(Fragments(atom)%BasisInfoIsSet .and. mynum==0) then
             call lsquit('mpi_bcast_many_fragments: Not implemented for &
                  & sending/receiving fragment basis info!',-1)
          end if

          ! Master: Copy fragment info into buffer
          ! Slave: Read fragment info from buffer into Fragments(atom)
          call mpicopy_fragment(Fragments(atom),comm,.false.)

          ! Check that the correct fragment number was communicated
          if(Fragments(atom)%atomic_number /= atom) then
             print '(a,2i8)', 'Counter atom, actual atom ', atom,Fragments(atom)%atomic_number
             call lsquit('mpi_bcast_many_fragments: Inconsistency between &
                  & counter and actual fragment',-1)
          end if

       end if CommunicateFragment

    end do

    ! Finalize buffer
    call ls_mpiFinalizeBuffer(master,LSMPIBROADCAST,comm)

  end subroutine mpi_bcast_many_fragments


  !> \brief MPI communcation where the information in the ccatom type
  !> is send (for local master) or received (for local slaves).
  !> In addition, other information required for calculating MP2 integrals
  !> and amplitudes is communicated.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine mpi_communicate_mp2_int_and_amp(MyFragment,bat,first_order_integrals,DoBasis)

    implicit none

    !> Fragment under consideration
    type(ccatom),intent(inout) :: MyFragment
    !> Batch information
    type(mp2_batch_construction),intent(inout) :: bat
    !> Do first order integrals?
    logical,intent(inout) :: first_order_integrals
    !> Communicate basis (expensive box in ccatom)
    logical,intent(in) :: DoBasis
    integer(kind=ls_mpik) :: master
    master = 0

    ! Init MPI buffer which eventually will contain all fragment info
    ! ***************************************************************
    ! MASTER: Prepare for writing to buffer
    ! SLAVE: Receive buffer
    call ls_mpiInitBuffer(master,LSMPIBROADCAST,infpar%lg_comm)


    ! Buffer handling
    ! ***************
    ! MASTER: Put information info buffer
    ! SLAVE: Put buffer information into ccatom structure

    ! Fragment information
    call mpicopy_fragment(MyFragment,infpar%lg_comm,DoBasis)

    ! Logical determining whether to calculate integrals for first order properties
    call ls_mpi_buffer(first_order_integrals,master)

    ! Batch information
    call ls_mpi_buffer(bat%MaxAllowedDimAlpha,master)
    call ls_mpi_buffer(bat%MaxAllowedDimGamma,master)
    call ls_mpi_buffer(bat%virtbatch,master)

    ! Finalize MPI buffer
    ! *******************
    ! MASTER: Send stuff to slaves and deallocate temp. buffers
    ! SLAVE: Deallocate buffer etc.
    call ls_mpiFinalizeBuffer(master,LSMPIBROADCAST,infpar%lg_comm)

    ! Copy 64 bit integers in their own bcast to minimize buffer complications...
    call ls_mpibcast(bat%size1,4,master,infpar%lg_comm)
    call ls_mpibcast(bat%size2,4,master,infpar%lg_comm)
    call ls_mpibcast(bat%size3,4,master,infpar%lg_comm)


  end subroutine mpi_communicate_mp2_int_and_amp


  !> \brief MPI communcation of full molecular information required by slaves
  !> to start fragment calculation.
  !> Assumes that basic dimension information (number of atoms and number of orbitals)
  !> have already been communicated.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine mpi_dec_fullinfo_master_to_slaves(natoms,nocc,nunocc,DistanceTable,&
       & OccOrbitals, UnoccOrbitals, MyMolecule, MyLsitem)

    implicit none

    !> Number of atoms in the molecule
    integer,intent(in) :: natoms
    !> Number of occupied orbitals in the molecule
    integer,intent(in) :: nocc
    !> Number of unoccupied orbitals in the molecule
    integer,intent(in) :: nunocc
    !> Distance table with inter-atomic distances
    !> Intent(in) for main master, intent(out) for local masters
    real(realk),intent(inout) :: DistanceTable(natoms,natoms)
    !> Occupied orbitals in DEC format
    !> Intent(in) for main master, intent(out) for local masters
    type(ccorbital),intent(inout) :: OccOrbitals(nocc)
    !> Unoccupied orbitals in DEC format
    !> Intent(in) for main master, intent(out) for local masters
    type(ccorbital),intent(inout) :: UnoccOrbitals(nunocc)
    !> Full molecule structure (MO coeffecients, Fock matrix etc.)
    !> Intent(in) for main master, intent(out) for local masters
    type(fullmolecule),intent(inout) :: MyMolecule
    !> Lsitem used for integral code
    !> Intent(in) for main master, intent(out) for local masters
    type(lsitem),intent(inout) :: MyLsitem
    integer(kind=ls_mpik) :: ierr, master
    integer :: i
    logical :: gm
    master=0

    ! Everything except the full molecule structure is passed using the LS buffer framework

    IF(infpar%mynum==0)THEN ! Global master
       gm=.true.
    else ! Slave
       gm=.false.
    end IF


    ! Init MPI buffer
    ! ***************
    ! Main master:  Prepare for writing to buffer
    ! Local master: Receive buffer
    call ls_mpiInitBuffer(master,LSMPIBROADCAST,MPI_COMM_LSDALTON)


    ! Buffer handling
    ! ***************
    ! Main master:  Put information info buffer
    ! Local master: Put buffer information into appropriate arrays



    ! DEC Occupied orbitals
    ! ---------------------

    ! Integers and reals inside ccorbital sub-type
    do i=1,nocc
       call ls_mpi_buffer(OccOrbitals(i)%orbitalnumber,master)
       call ls_mpi_buffer(OccOrbitals(i)%centralatom,master)
       call ls_mpi_buffer(OccOrbitals(i)%numberofatoms,master)

       ! Pointer for atoms for orbital "i"
       if(.not. gm) then ! allocate for local masters (not global master)
          Nullify(OccOrbitals(i)%atoms)
          call mem_alloc(OccOrbitals(i)%atoms,OccOrbitals(i)%numberofatoms)
       end if
       ! Buffer handling for both master and slave
       call ls_mpi_buffer(OccOrbitals(i)%atoms, &
            & OccOrbitals(i)%numberofatoms,master)
    end do


    ! DEC Unoccupied orbitals
    ! -----------------------

    ! Integers and reals inside ccorbital sub-type
    do i=1,nunocc
       call ls_mpi_buffer(UnoccOrbitals(i)%orbitalnumber,master)
       call ls_mpi_buffer(UnoccOrbitals(i)%centralatom,master)
       call ls_mpi_buffer(UnoccOrbitals(i)%numberofatoms,master)

       ! Pointer for atoms for orbital "i"
       if(.not. gm) then ! allocate for local masters (not global master)
          Nullify(UnoccOrbitals(i)%atoms)
          call mem_alloc(UnoccOrbitals(i)%atoms,UnoccOrbitals(i)%numberofatoms)
       end if
       ! Buffer handling for both master and slave
       call ls_mpi_buffer(UnoccOrbitals(i)%atoms, &
            & UnoccOrbitals(i)%numberofatoms,master)
    end do


    ! Distance table
    ! --------------
    call ls_mpi_buffer(DistanceTable,natoms,natoms,master)

    ! Integral lsitem
    ! ---------------
    call mpicopy_lsitem(MyLsitem,MPI_COMM_LSDALTON)


    ! Finalize MPI buffer
    ! *******************
    ! Main master:  Send stuff to local masters and deallocate temp. buffers
    ! Local master: Deallocate buffer etc.
    call ls_mpiFinalizeBuffer(master,LSMPIBROADCAST,MPI_COMM_LSDALTON)


    ! Full molecule bcasting
    ! **********************
    call mpi_bcast_fullmolecule(MyMolecule)


  end subroutine mpi_dec_fullinfo_master_to_slaves


  !> \brief MPI Bcast fullmolecule structure from global master to all slaves.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine mpi_bcast_fullmolecule(MyMolecule)

    implicit none

    !> Full molecule structure (MO coeffecients, Fock matrix etc.)
    !> Intent(in) for main master, intent(out) for local masters
    type(fullmolecule),intent(inout) :: MyMolecule
    logical :: gm
    integer(kind=ls_mpik) :: master
    master = 0

    IF(infpar%mynum==0)THEN ! Global master
       gm=.true.
    else ! Local master
       gm=.false.
    end IF

    ! Simple integers
    ! ---------------
    call ls_mpibcast(MyMolecule%nelectrons,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%natoms,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%nbasis,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%nauxbasis,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%numocc,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%ncore,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%nval,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%numvirt,master,MPI_COMM_LSDALTON)

    ! Allocate pointers if local master
    if(.not. gm) then
       call mem_alloc(MyMolecule%atom_size,MyMolecule%natoms)
       call mem_alloc(MyMolecule%atom_start,MyMolecule%natoms)
       call mem_alloc(MyMolecule%atom_end,MyMolecule%natoms)
       call mem_alloc(MyMolecule%ypo,MyMolecule%nbasis,MyMolecule%numocc)
       call mem_alloc(MyMolecule%ypv,MyMolecule%nbasis,MyMolecule%numvirt)
       call mem_alloc(MyMolecule%fock,MyMolecule%nbasis,MyMolecule%nbasis)
       call mem_alloc(MyMolecule%overlap,MyMolecule%nbasis,MyMolecule%nbasis)
       call mem_alloc(MyMolecule%ppfock,MyMolecule%numocc,MyMolecule%numocc)
       call mem_alloc(MyMolecule%qqfock,MyMolecule%numvirt,MyMolecule%numvirt)
    end if


    ! Integer pointers
    ! ----------------
    call ls_mpibcast(MyMolecule%atom_size,MyMolecule%natoms,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%atom_start,MyMolecule%natoms,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%atom_end,MyMolecule%natoms,master,MPI_COMM_LSDALTON)


    ! Real pointers
    ! -------------
    call ls_mpibcast(MyMolecule%ypo,MyMolecule%nbasis,MyMolecule%numocc,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%ypv,MyMolecule%nbasis,MyMolecule%numvirt,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%fock,MyMolecule%nbasis,MyMolecule%nbasis,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%overlap,MyMolecule%nbasis,MyMolecule%nbasis,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%ppfock,MyMolecule%numocc,MyMolecule%numocc,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%qqfock,MyMolecule%numvirt,MyMolecule%numvirt,master,MPI_COMM_LSDALTON)


  end subroutine mpi_bcast_fullmolecule



  !> \brief MPI copy information in ccatom type.
  !> Note 1: The information is copied into the global buffers here,
  !> and it is assumed that ls_mpiInitBuffer (ls_mpiFinalizeBuffer)
  !> is called before (after) calling this subroutine.
  !> Note 2: If the global logical AddToBuffer is FALSE then
  !> pointers in fragment type are allocated here and the approriate information
  !> is read from the buffer. If AddToBuffer is TRUE,
  !> then it is assumed that the pointers in MyFragment have already been allocated
  !> and contains the relevant information to write to buffer.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine mpicopy_fragment(MyFragment,comm,DoBasis)

    implicit none

    !> Fragment under consideration
    type(ccatom),intent(inout) :: MyFragment
    !> Communicator
    integer(kind=ls_mpik),intent(in) :: comm
    !> Copy basis (expensive box in ccatom)
    logical,intent(in) :: DoBasis
    integer :: i
    integer(kind=ls_mpik) :: master
    master = 0


    ! Integers that are not pointers
    ! ------------------------------

    CALL ls_mpi_buffer(MyFragment%atomic_number,master)
    CALL ls_mpi_buffer(MyFragment%noccEOS,master)
    CALL ls_mpi_buffer(MyFragment%nunoccEOS,master)
    CALL ls_mpi_buffer(MyFragment%noccAOS,master)
    CALL ls_mpi_buffer(MyFragment%ncore,master)
    CALL ls_mpi_buffer(MyFragment%nocctot,master)
    CALL ls_mpi_buffer(MyFragment%nunoccAOS,master)
    CALL ls_mpi_buffer(MyFragment%REDnoccAOS,master)
    CALL ls_mpi_buffer(MyFragment%REDnunoccAOS,master)
    CALL ls_mpi_buffer(MyFragment%nEOSatoms,master)
    CALL ls_mpi_buffer(MyFragment%ntasks,master)
    CALL ls_mpi_buffer(MyFragment%t1dims,2,master)
    CALL ls_mpi_buffer(MyFragment%noccFA,master)
    CALL ls_mpi_buffer(MyFragment%nunoccFA,master)
    CALL ls_mpi_buffer(MyFragment%number_atoms,master)
    CALL ls_mpi_buffer(MyFragment%number_basis,master)


    ! Logicals that are not pointers
    ! ------------------------------
    CALL ls_mpi_buffer(MyFragment%BasisInfoIsSet,master)
    CALL ls_mpi_buffer(MyFragment%t1_stored,master)
    CALL ls_mpi_buffer(MyFragment%CDset,master)
    CALL ls_mpi_buffer(MyFragment%FAtransSet,master)
    CALL ls_mpi_buffer(MyFragment%fragmentadapted,master)
    CALL ls_mpi_buffer(MyFragment%pairfrag,master)

    ! Reals that are not pointers
    ! ---------------------------
    ! (Note: This - and much else MPI stuff - must be modified for single precision to work!)
    call ls_mpi_buffer(MyFragment%energies,ndecenergies,master)
    call ls_mpi_buffer(MyFragment%pairdist,master)
    call ls_mpi_buffer(MyFragment%EoccFOP,master)
    call ls_mpi_buffer(MyFragment%EvirtFOP,master)
    call ls_mpi_buffer(MyFragment%LagFOP,master)
    CALL ls_mpi_buffer(MyFragment%flops_slaves,master)
    call ls_mpi_buffer(MyFragment%slavetime,master)
    call ls_mpi_buffer(MyFragment%RejectThr,2,master)


    ! Integer pointers
    ! ----------------
    ! Nullify and allocate stuff for receiver (global addtobuffer is false)
    if(.not. AddToBuffer) then
       nullify(MyFragment%occEOSidx)
       call mem_alloc(MyFragment%occEOSidx,MyFragment%noccEOS)
       nullify(MyFragment%unoccEOSidx)
       call mem_alloc(MyFragment%unoccEOSidx,MyFragment%nunoccEOS)
       nullify(MyFragment%occAOSidx)
       call mem_alloc(MyFragment%occAOSidx,MyFragment%noccAOS)
       nullify(MyFragment%unoccAOSidx)
       call mem_alloc(MyFragment%unoccAOSidx,MyFragment%nunoccAOS)
       nullify(MyFragment%coreidx)
       if(MyFragment%ncore>0) then
          call mem_alloc(MyFragment%coreidx,MyFragment%ncore)
       end if
       nullify(MyFragment%REDoccAOSidx)
       call mem_alloc(MyFragment%REDoccAOSidx,MyFragment%REDnoccAOS)
       nullify(MyFragment%REDunoccAOSidx)
       call mem_alloc(MyFragment%REDunoccAOSidx,MyFragment%REDnunoccAOS)
       nullify(MyFragment%idxo)
       call mem_alloc(MyFragment%idxo,MyFragment%noccEOS)
       nullify(MyFragment%idxu)
       call mem_alloc(MyFragment%idxu,MyFragment%nunoccEOS)
       nullify(MyFragment%EOSatoms)
       call mem_alloc(MyFragment%EOSatoms,MyFragment%nEOSatoms)
       if(MyFragment%t1_stored) then ! only used for CC singles effects
          nullify(MyFragment%t1_occidx)
          call mem_alloc(MyFragment%t1_occidx,MyFragment%t1dims(2) )
          nullify(MyFragment%t1_virtidx)
          call mem_alloc(MyFragment%t1_virtidx,MyFragment%t1dims(1) )
       end if
       nullify(MyFragment%atoms_idx)
       call mem_alloc(MyFragment%atoms_idx,MyFragment%number_atoms)
       nullify(MyFragment%basis_idx)
       call mem_alloc(MyFragment%basis_idx,MyFragment%number_basis)
    end if


    ! Buffer handling
    call ls_mpi_buffer(MyFragment%occEOSidx,MyFragment%noccEOS,master)
    call ls_mpi_buffer(MyFragment%unoccEOSidx,MyFragment%nunoccEOS,master)
    call ls_mpi_buffer(MyFragment%occAOSidx,MyFragment%noccAOS,master)
    call ls_mpi_buffer(MyFragment%unoccAOSidx,MyFragment%nunoccAOS,master)
    if(MyFragment%ncore>0) then
       call ls_mpi_buffer(MyFragment%coreidx,MyFragment%ncore,master)
    end if
    call ls_mpi_buffer(MyFragment%REDoccAOSidx,&
         & MyFragment%REDnoccAOS,master)
    call ls_mpi_buffer(MyFragment%REDunoccAOSidx,&
         MyFragment%REDnunoccAOS,master)
    call ls_mpi_buffer(MyFragment%idxo,MyFragment%noccEOS,master)
    call ls_mpi_buffer(MyFragment%idxu,MyFragment%nunoccEOS,master)
    call ls_mpi_buffer(MyFragment%EOSatoms,MyFragment%nEOSatoms,master)
    call ls_mpi_buffer(MyFragment%atoms_idx,MyFragment%number_atoms,master)
    call ls_mpi_buffer(MyFragment%basis_idx,MyFragment%number_basis,master)
    if(MyFragment%t1_stored) then ! only used for CC singles effects
       call ls_mpi_buffer(MyFragment%t1_occidx,MyFragment%t1dims(2),master)
       call ls_mpi_buffer(MyFragment%t1_virtidx,MyFragment%t1dims(1),master)
    end if


    ! Real pointers
    ! -------------
    ! (Note: This and much else MPI stuff must be modified for single precision to work!)
    if(.not. AddToBuffer) then
       nullify(MyFragment%OccContribs)
       call mem_alloc(MyFragment%OccContribs,MyFragment%noccAOS)
       nullify(MyFragment%VirtContribs)
       call mem_alloc(MyFragment%VirtContribs,MyFragment%nunoccAOS)

       if(MyFragment%CDset) then
          nullify(MyFragment%OccMat)
          call mem_alloc(MyFragment%OccMat,MyFragment%noccAOS,MyFragment%noccAOS)
          nullify(MyFragment%VirtMat)
          call mem_alloc(MyFragment%VirtMat,MyFragment%nunoccAOS,MyFragment%nunoccAOS)
       end if

       if(MyFragment%FAtransSet) then
          nullify(MyFragment%CoccFA)
          call mem_alloc(MyFragment%CoccFA,MyFragment%number_basis,MyFragment%noccFA)
          nullify(MyFragment%CunoccFA)
          call mem_alloc(MyFragment%CunoccFA,MyFragment%number_basis,MyFragment%nunoccFA)
       end if

       if (MyFragment%t1_stored) then ! only used for CC singles effects
          nullify(MyFragment%t1)
          call mem_alloc(MyFragment%t1,MyFragment%t1dims(1),MyFragment%t1dims(2) )
       end if
    end if

    call ls_mpi_buffer(MyFragment%OccContribs,MyFragment%noccAOS,master)
    call ls_mpi_buffer(MyFragment%VirtContribs,MyFragment%nunoccAOS,master)
    if(MyFragment%CDset) then
       call ls_mpi_buffer(MyFragment%OccMat,MyFragment%noccAOS,MyFragment%noccAOS,master)
       call ls_mpi_buffer(MyFragment%VirtMat,MyFragment%nunoccAOS,MyFragment%nunoccAOS,master)
    end if
    if(MyFragment%FAtransSet) then
       call ls_mpi_buffer(MyFragment%CoccFA,MyFragment%number_basis,MyFragment%noccFA,master)
       call ls_mpi_buffer(MyFragment%CunoccFA,MyFragment%number_basis,MyFragment%nunoccFA,master)
    end if
    if (MyFragment%t1_stored) then ! only used for CC singles effects
       call ls_mpi_buffer(MyFragment%t1,MyFragment%t1dims(1),MyFragment%t1dims(2),master)
    end if


    ! CCORBITAL types
    ! ---------------
    ! Allocate ccorbitals
    if(.not. AddToBuffer) then
       nullify(MyFragment%occAOSorb)
       call mem_alloc(MyFragment%occAOSorb,MyFragment%noccAOS)
       nullify(MyFragment%unoccAOSorb)
       call mem_alloc(MyFragment%unoccAOSorb,MyFragment%nunoccAOS)
    end if

    ! Integers and reals inside ccorbital sub-type
    do i=1,MyFragment%noccAOS ! occ orbitals
       call ls_mpi_buffer(MyFragment%occAOSorb(i)%orbitalnumber,master)
       call ls_mpi_buffer(MyFragment%occAOSorb(i)%centralatom,master)
       call ls_mpi_buffer(MyFragment%occAOSorb(i)%numberofatoms,master)

       ! Pointers inside ccorbital sub-type (which again is inside ccatom type)
       if(.not. AddToBuffer) then
          Nullify(MyFragment%occAOSorb(i)%atoms)
          call mem_alloc(MyFragment%occAOSorb(i)%atoms,&
               & MyFragment%occAOSorb(i)%numberofatoms)
       end if
       ! Buffer handling
       call ls_mpi_buffer(MyFragment%occAOSorb(i)%atoms, &
            & MyFragment%occAOSorb(i)%numberofatoms,master)
    end do

    do i=1,MyFragment%nunoccAOS ! unocc orbitals
       call ls_mpi_buffer(MyFragment%unoccAOSorb(i)%orbitalnumber,master)
       call ls_mpi_buffer(MyFragment%unoccAOSorb(i)%centralatom,master)
       call ls_mpi_buffer(MyFragment%unoccAOSorb(i)%numberofatoms,master)
       if(.not. AddToBuffer) then ! allocate for slave
          Nullify(MyFragment%unoccAOSorb(i)%atoms)
          call mem_alloc(MyFragment%unoccAOSorb(i)%atoms,&
               & MyFragment%unoccAOSorb(i)%numberofatoms)
       end if
       call ls_mpi_buffer(MyFragment%unoccAOSorb(i)%atoms, &
            & MyFragment%unoccAOSorb(i)%numberofatoms,master)
    end do



    ! ===========================================================================
    !                              EXPENSIVE BOX
    ! ===========================================================================
    ! This is the information inside the "expensive" box in the ccatom type definition.
    ! This stuff is only copied if MyFragment%BasisInfoIsSet is true.
    ExpensiveBox: if(DoBasis) then

       ! Real pointers
       if(.not. AddToBuffer) then
          call mem_alloc(MyFragment%S,MyFragment%number_basis,MyFragment%number_basis)
          call mem_alloc(MyFragment%ypo,MyFragment%number_basis,MyFragment%noccAOS)
          call mem_alloc(MyFragment%ypv,MyFragment%number_basis,MyFragment%nunoccAOS)
          call mem_alloc(MyFragment%coreMO,MyFragment%number_basis,MyFragment%ncore)
          call mem_alloc(MyFragment%fock,MyFragment%number_basis,MyFragment%number_basis)
          call mem_alloc(MyFragment%ppfock,MyFragment%noccAOS,MyFragment%noccAOS)
          call mem_alloc(MyFragment%ccfock,MyFragment%ncore,MyFragment%ncore)
          call mem_alloc(MyFragment%qqfock,MyFragment%nunoccAOS,MyFragment%nunoccAOS)
       end if
       call ls_mpi_buffer(MyFragment%S,MyFragment%number_basis,MyFragment%number_basis,master)
       call ls_mpi_buffer(MyFragment%ypo,MyFragment%number_basis,MyFragment%noccAOS,master)
       call ls_mpi_buffer(MyFragment%ypv,MyFragment%number_basis,MyFragment%nunoccAOS,master)
       call ls_mpi_buffer(MyFragment%coreMO,MyFragment%number_basis,MyFragment%ncore,master)
       call ls_mpi_buffer(MyFragment%fock,MyFragment%number_basis,MyFragment%number_basis,master)
       call ls_mpi_buffer(MyFragment%ppfock,MyFragment%noccAOS,MyFragment%noccAOS,master)
       call ls_mpi_buffer(MyFragment%ccfock,MyFragment%ncore,MyFragment%ncore,master)
       call ls_mpi_buffer(MyFragment%qqfock,MyFragment%nunoccAOS,&
            &MyFragment%nunoccAOS,master)

       ! INTEGRAL LSITEM
       ! '''''''''''''''
       call mpicopy_lsitem(MyFragment%MyLsitem,comm)

       ! Basis info has now been set
       MyFragment%BasisInfoIsSet=.true.

    else

       ! Basis info was not set
       MyFragment%BasisInfoIsSet=.false.

    end if ExpensiveBox


  End subroutine mpicopy_fragment

  subroutine share_E2_with_slaves(ppf,qqf,t2,xo,yv,Gbi,Had,no,nv,nb,omega2,s)
    implicit none
    real(realk),pointer :: xo(:),yv(:),Gbi(:),Had(:)
    real(realk), intent(inout) :: ppf(:),qqf(:)
    integer :: no,nv,nb,s
    type(array),intent(inout) :: t2,omega2
    integer :: oaddr(infpar%lg_nodtot)
    integer :: taddr(infpar%lg_nodtot)
    logical :: master

    master=(infpar%lg_mynum==infpar%master)
    !print *, infpar%lg_mynum,"comming",master
    if(master) call ls_mpibcast(CCSDSLV4E2,infpar%master,infpar%lg_comm)
    
    call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
    call ls_mpi_buffer(no,infpar%master)
    call ls_mpi_buffer(nv,infpar%master)
    call ls_mpi_buffer(nb,infpar%master)
    call ls_mpi_buffer(s,infpar%master)
    if(master)oaddr=omega2%addr_p_arr
    call ls_mpi_buffer(oaddr,infpar%lg_nodtot,infpar%master)
    if(master)taddr=t2%addr_p_arr
    call ls_mpi_buffer(taddr,infpar%lg_nodtot,infpar%master)
    if(.not.master)then
      call mem_alloc(xo,nb*no)
      call mem_alloc(yv,nb*nv)
      call mem_alloc(Gbi,nb*no)
      call mem_alloc(Had,nv*nb)
    endif
    call ls_mpi_buffer(xo,nb*no,infpar%master)
    call ls_mpi_buffer(yv,nb*nv,infpar%master)
    call ls_mpi_buffer(Gbi,nb*no,infpar%master)
    call ls_mpi_buffer(Had,nv*nb,infpar%master)
    call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
    if(master)then
      call ls_mpibcast(ppf,no*no,infpar%master,infpar%lg_comm)
      call ls_mpibcast(qqf,nv*nv,infpar%master,infpar%lg_comm)
    else
      t2     = get_arr_from_parr(taddr(infpar%lg_mynum+1))
      omega2 = get_arr_from_parr(oaddr(infpar%lg_mynum+1))
    endif
  end subroutine share_E2_with_slaves



  !> \brief MPI communcation where CCSD and CC2 data is transferred
  !> \author Patrick Ettenhuber
  !> \date March 2012
  subroutine mpi_communicate_ccsd_calcdata(om2,t2,govov,xo,xv,yo,yv,MyLsItem,nbas,nvirt,nocc,iter,s,p)
    implicit none
    type(mp2_batch_construction) :: bat
    integer            :: nbas,nocc,nvirt,ierr,iter,s
    !real(realk)        :: t2(:),govov(:)
    type(array),intent(inout) :: t2,govov,om2
    real(realk)        :: xo(:),xv(:),yo(:),yv(:)
    type(lsitem)       :: MyLsItem
    real(realk)        :: norm
    integer(kind=long) :: nelms
    integer :: i,n4,k
    integer :: gaddr(infpar%lg_nodtot)
    integer :: taddr(infpar%lg_nodtot)
    integer :: oaddr(infpar%lg_nodtot)
    logical :: p
    character(ARR_MSG_LEN) :: msg
    logical :: master
    master=(infpar%lg_mynum==infpar%master)
   !communicate mylsitem and integers
    call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
    call ls_mpi_buffer(DECinfo%ccModel,infpar%master)
    call ls_mpi_buffer(nbas,infpar%master)
    call ls_mpi_buffer(nocc,infpar%master)
    call ls_mpi_buffer(nvirt,infpar%master)
    call ls_mpi_buffer(iter,infpar%master)
    call ls_mpi_buffer(s,infpar%master)
    call ls_mpi_buffer(p,infpar%master)
    if(p)then
      if(master)gaddr=govov%addr_p_arr
      call ls_mpi_buffer(gaddr,infpar%lg_nodtot,infpar%master)
      if(master)taddr=t2%addr_p_arr
      call ls_mpi_buffer(taddr,infpar%lg_nodtot,infpar%master)
      if(master)oaddr=om2%addr_p_arr
      call ls_mpi_buffer(oaddr,infpar%lg_nodtot,infpar%master)
    endif
    call mpicopy_lsitem(MyLsItem,infpar%lg_comm)
    call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)

    !communicate rest of the quantities, master here, slaves back in the slave
    !routine, due to crappy pointer/non-pointer issues (->allocations)
    if(master)then

      !split messages in 2GB parts, compare to counterpart in
      !ccsd_data_preparation
      k=250000000

      nelms = nbas*nocc
      call ls_mpibcast_chunks(xo,nelms,infpar%master,infpar%lg_comm,k)
      call ls_mpibcast_chunks(yo,nelms,infpar%master,infpar%lg_comm,k)

      nelms = nbas*nvirt
      call ls_mpibcast_chunks(xv,nelms,infpar%master,infpar%lg_comm,k)
      call ls_mpibcast_chunks(yv,nelms,infpar%master,infpar%lg_comm,k)

      nelms = nvirt*nvirt*nocc*nocc
      call ls_mpibcast_chunks(t2%elm1,nelms,infpar%master,infpar%lg_comm,k)
      if(iter/=1.and.(s==0.or.s==4))then
        call ls_mpibcast_chunks(govov%elm1,nelms,infpar%master,infpar%lg_comm,k)
      endif
    else
      if(p)then
        govov = get_arr_from_parr(gaddr(infpar%lg_mynum+1))
        t2    = get_arr_from_parr(taddr(infpar%lg_mynum+1))
        om2   = get_arr_from_parr(oaddr(infpar%lg_mynum+1))
      endif
    endif
  end subroutine mpi_communicate_ccsd_calcdata

  !> \brief mpi communcation where ccsd(t) data is transferred
  !> \author Janus Juul Eriksen
  !> \date February 2013
  subroutine mpi_communicate_ccsdpt_calcdata(nocc,nvirt,nbasis,ppfock,qqfock,ypo,ypv,ccsd_t2,mylsitem)

    implicit none

    integer            :: nocc,nvirt,nbasis,ierr
    real(realk)        :: ppfock(:,:),qqfock(:,:),ypo(:,:),ypv(:,:),ccsd_t2(:,:,:,:)
    type(lsitem)       :: mylsitem

    ! communicate mylsitem and integers
    call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
    call ls_mpi_buffer(DECinfo%ccModel,infpar%master)
    call ls_mpi_buffer(DECinfo%memory,infpar%master)
    call ls_mpi_buffer(nbasis,infpar%master)
    call ls_mpi_buffer(nocc,infpar%master)
    call ls_mpi_buffer(nvirt,infpar%master)
    call mpicopy_lsitem(mylsitem,infpar%lg_comm)
    call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)

    ! communicate rest of the quantities, master here, slaves back in the slave
    ! routine, due to crappy pointer/non-pointer issues (->allocations)
    if (infpar%lg_mynum .eq. infpar%master) then

       call ls_mpibcast(ppfock,nocc,nocc,infpar%master,infpar%lg_comm)
 
       call ls_mpibcast(qqfock,nvirt,nvirt,infpar%master,infpar%lg_comm)

       call ls_mpibcast(ypo,nbasis,nocc,infpar%master,infpar%lg_comm)

       call ls_mpibcast(ypv,nbasis,nvirt,infpar%master,infpar%lg_comm)
 
       call ls_mpibcast(ccsd_t2,nvirt,nocc,nvirt,nocc,infpar%master,infpar%lg_comm)

    endif

  end subroutine mpi_communicate_ccsdpt_calcdata


  !> \brief get a suitable job distribution in mpi calculations
  !> \author Patrick Ettenhuber
  !> \date March 2012
  subroutine distribute_mpi_jobs(mpi_task_distribution,nbA,nbG,batchdimAlpha,&
             & batchdimGamma,myload,ccsdscheme,no,nv,nb,b2oa,b2og)
    implicit none

    integer, intent(inout) :: myload
    integer,intent(in):: nbA,nbG,batchdimAlpha(nbA),batchdimGamma(nbG)
    integer,intent(out)::mpi_task_distribution(nbA*nbG)
    integer,intent(in),optional :: ccsdscheme,no,nv,nb
    type(batchtoorb),intent(in),optional :: b2oa(nbA),b2og(nbG)
    integer(kind=8),pointer :: workloads(:),jobsize_per_node(:)
    integer(kind=8) :: myload64
    integer(kind=8),pointer :: node_rank_by_job(:),easytrace(:)
    logical,pointer :: touched(:)
    logical :: all_touched
    type(traceback), pointer::trace(:)
    integer :: i,j,k,l,swap,counter,ialpha,igamma,actual_node
    integer(kind=8) :: fa,fg,la,lg

    myload64=0
    nullify(workloads)
    nullify(trace)
    nullify(touched)
    nullify(easytrace)
    nullify(jobsize_per_node)
    nullify(node_rank_by_job)
    call  mem_alloc(touched,nbA*nbG)
    call  mem_alloc(workloads,nbA*nbG)
    call  mem_alloc(easytrace,nbA*nbG)
    call  mem_alloc(jobsize_per_node,infpar%lg_nodtot)
    call  mem_alloc(node_rank_by_job,infpar%lg_nodtot)
    call  mem_alloc(trace,nbA*nbG)
    !some needed initializations for the matrices
    touched = .false.
    jobsize_per_node = 0
    do k=1, infpar%lg_nodtot
      node_rank_by_job(k) = k-1
    enddo

    if(DECinfo%PL>0) then
       print '(a,4i6)', 'LOCAL,GLOBAL,sizeL,sizeG ', &
            & infpar%lg_mynum, infpar%mynum, infpar%lg_nodtot, infpar%nodtot
    end if

    if (.false.) then
     ! simple modulo distribtution --> tured off for the moment
     ! fill the jobsize_per_node matrix with the values from mpi_task_distribution
      do i=1,nbG
        do j=1,nbA
          actual_node=MODULO(j-1+(i-1)*MODULO(nbA,infpar%lg_nodtot),infpar%lg_nodtot)
          mpi_task_distribution((j-1)*nbG+i)=actual_node
          jobsize_per_node(actual_node+1)=jobsize_per_node(actual_node+1) + batchdimGamma(i)*batchdimAlpha(j)
          touched((j-1)*nbg+i) = .true.
        enddo
      enddo


    else
      ! here the actual routine is located -> jobs are assigned according to their sizes

      ! Calculate the sizes of the individual jobs while labeling them with a uniqe number saved
      ! in easytrace. The trace variable contains the the indices of the batches. I think in this
      ! way it is easier to see what is done instead of recalculating the indices from the combined
      ! index
      counter = 0
      do j=1,nbA
        do i=1,nbG
         la=batchdimAlpha(j)
         lg=batchdimGamma(i)
         workloads((j-1)*nbG+i)   = la*lg

         !in CCSD the workloads are not equally distributed, and a more accurate
         !estimation of the work has to be done
         if(present(ccsdscheme))then
           fa=b2oA(j)%orbindex(1)
           fg=b2oG(i)%orbindex(1)
           !account for integral work more "correctly"
           workloads((j-1)*nbG+i) = workloads((j-1)*nbG+i) * nb * nb * 2
           !account for Kobayashi-Term
           if(fa<=fg+lg-1 )then
             workloads((j-1)*nbG+i) = workloads((j-1)*nbG+i) + (nv**2*no**2*la*lg)/4
           endif
         endif

         easytrace((j-1)*nbG+i)   = counter
         trace((j-1)*nbG+i)%ident = counter
         trace((j-1)*nbG+i)%na    = j
         trace((j-1)*nbG+i)%ng    = i
         counter = counter + 1
        enddo
      enddo

      ! Sort the jobs according to their size, and keep track of indices
      do i=1,nbG*nbA
        do j=i+1,nbg*nbA
          if( workloads(j) > workloads(i) )then

            swap=workloads(j)
            workloads(j)=workloads(i)
            workloads(i)=swap

            swap=easytrace(j)
            easytrace(j)=easytrace(i)
            easytrace(i)=swap
          endif
        enddo
      enddo

      ! assign new jobs and thereby keep track of the sizes, assign ranks according
      ! to workload after each round
      do i=1,nbG*nbA-MODULO(nbG*nbA,infpar%lg_nodtot),infpar%lg_nodtot
        !distribute the work according to previous load
        do k=1, infpar%lg_nodtot
          jobsize_per_node(k)=jobsize_per_node(k) + workloads(i+node_rank_by_job(k))
          ialpha = trace(easytrace(i+node_rank_by_job(k))+1)%na
          igamma = trace(easytrace(i+node_rank_by_job(k))+1)%ng
          mpi_task_distribution((ialpha-1)*nbg+igamma) = k-1

          ! keep track of which elements were done if an error occurs, quit, since
          ! the result would
          if(.not. touched((ialpha-1)*nbg+igamma)) then
            touched((ialpha-1)*nbg+igamma) = .true.
          else
            if (infpar%lg_mynum==0)then
              call LSQUIT('Element had already been assigned in distribute_mpi_jobs',-1)
            endif
          endif
        enddo

        ! get the new node ranks by counting how many of the nodes have smaller jobs
        ! --> is an initial rank, because it is ordered the other way round as the jobs
        do k=1, infpar%lg_nodtot
          counter = 0
          do l=1,infpar%lg_nodtot
            if (jobsize_per_node(l) < jobsize_per_node(k)) counter=counter+1
          enddo
          node_rank_by_job(k) = counter
        enddo
        ! if some have the same amount of work, just assign jobs uniquely
        do k=1, infpar%lg_nodtot
          do l=k+1,infpar%lg_nodtot
            if (node_rank_by_job(k)==node_rank_by_job(l))then
              node_rank_by_job(l) = node_rank_by_job(l) +1
            endif
          enddo
        enddo
      enddo

      !distribute the rest of the jobs according to load
      i = nbG*nbA-MODULO(nbG*nbA,infpar%lg_nodtot)
      do k=1, infpar%lg_nodtot
        actual_node= node_rank_by_job(k) + 1
        if(MODULO(nbG*nbA,infpar%lg_nodtot) >= actual_node)then
          jobsize_per_node(k)=jobsize_per_node(k) + workloads(i+actual_node)
          ialpha = trace(easytrace(i+actual_node)+1)%na
          igamma = trace(easytrace(i+actual_node)+1)%ng
          mpi_task_distribution((ialpha-1)*nbg+igamma) = k-1
          if(.not. touched((ialpha-1)*nbg+igamma)) then
            touched((ialpha-1)*nbg+igamma) = .true.
          else
            if (infpar%lg_mynum==0)then
              call LSQUIT('Element had already been assigned in distribute_mpi_jobs',-1)
            endif
          endif
        endif
      enddo
    endif

    ! some sanity checks, that the sizes calculated in the routine do not differ from the
    ! jobsizes calculated with the distribution matrix, and that all items were redistributed
    all_touched = .true.
    do j=1,nbA
      do i=1,nbG
        if (mpi_task_distribution((j-1)*nbG+i)==infpar%lg_mynum)then
           la=batchdimAlpha(j)
           lg=batchdimGamma(i)
           if(present(ccsdscheme))then
             fa=b2oA(j)%orbindex(1)
             fg=b2oG(i)%orbindex(1)
             !account for integral work more "correctly"
             myload64 = myload64 + la*lg*nb*nb*2
             !account for Kobayashi-Term
             if(fa<=fg+lg-1 )then
               myload64 = myload64 + (nv**2*no**2*la*lg)/4
             endif
           else
             myload64 = myload64 + la*lg
           endif
        endif
        if (.not. touched((j-1)*nbg+i)) all_touched = .false.
      enddo
    enddo

    if ((myload64 /= jobsize_per_node(infpar%lg_mynum+1)).or. .not. all_touched)then
      print *, "something is wrong ... quitting",myload64,jobsize_per_node(infpar%lg_mynum+1),all_touched
      write(DECinfo%output,*) 'workloads per node:'
      write(DECinfo%output,*) node_rank_by_job
      write(DECinfo%output,*) jobsize_per_node
      write(DECinfo%output,*) mpi_task_distribution
      print *, 'workloads per node:'
      print *, node_rank_by_job
      print *, jobsize_per_node
      print *, mpi_task_distribution
      print *, touched
      call LSQUIT('Problem in distribute_mpi_jobs (64bit 32bit??), (jobsize or unused elements)',-1)
    endif
    myload=myload64

    !deallocating and nullifying pointers
    call  mem_dealloc(workloads)
    call  mem_dealloc(touched)
    call  mem_dealloc(easytrace)
    call  mem_dealloc(jobsize_per_node)
    call  mem_dealloc(node_rank_by_job)
    call mem_dealloc(trace)

    nullify(workloads)
    nullify(easytrace)
    nullify(trace)
    nullify(touched)
    nullify(jobsize_per_node)
    nullify(node_rank_by_job)

  end subroutine distribute_mpi_jobs


  !> \brief Bcast fragment joblist from master to slaves
  !> for fragment job list AFTER fragment optimization is done.
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine bcast_post_fragopt_joblist(jobs,comm)

    implicit none
    !> Job list (send from master to slaves via communicator)
    type(joblist),intent(inout) :: jobs
    !> Communicator
    integer(kind=ls_mpik),intent(in) :: comm
    integer(kind=ls_mpik) :: master
    master=0


    ! Init buffer
    call ls_mpiInitBuffer(master,LSMPIBROADCAST,comm)

    ! Master: Copy job info into buffer
    ! Slave: Read job info from buffer into Fragments(atom)
    call mpicopy_fragment_joblist(jobs)

    ! Finalize buffer
    call ls_mpiFinalizeBuffer(master,LSMPIBROADCAST,comm)


  end subroutine bcast_post_fragopt_joblist



  !> \brief MPI copy fragment joblist from master to slaves using MPI buffer
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine mpicopy_fragment_joblist(jobs)

    implicit none
    !> Job list (send from master to slaves via communicator)
    type(joblist),intent(inout) :: jobs
    integer(kind=ls_mpik) :: master
    master=0


    ! Master: Copy job info into buffer
    ! Slave: Read job info from buffer into Fragments(atom)


    ! Number of jobs
    call ls_mpi_buffer(jobs%njobs,master)

    ! Allocate pointers for slave
    if(.not. AddToBuffer) then   ! this is a slave
       call mem_alloc(jobs%atom1,jobs%njobs)
       call mem_alloc(jobs%atom2,jobs%njobs)
       call mem_alloc(jobs%jobsize,jobs%njobs)
       call mem_alloc(jobs%jobsdone,jobs%njobs)
       call mem_alloc(jobs%nslaves,jobs%njobs)
       call mem_alloc(jobs%nocc,jobs%njobs)
       call mem_alloc(jobs%nvirt,jobs%njobs)
       call mem_alloc(jobs%nbasis,jobs%njobs)
       call mem_alloc(jobs%ntasks,jobs%njobs)
       call mem_alloc(jobs%flops,jobs%njobs)
       call mem_alloc(jobs%LMtime,jobs%njobs)
       call mem_alloc(jobs%load,jobs%njobs)
    end if

    ! Buffer handling for pointers
    call ls_mpi_buffer(jobs%atom1,jobs%njobs,master)
    call ls_mpi_buffer(jobs%atom2,jobs%njobs,master)
    call ls_mpi_buffer(jobs%jobsize,jobs%njobs,master)
    call ls_mpi_buffer(jobs%jobsdone,jobs%njobs,master)
    call ls_mpi_buffer(jobs%nslaves,jobs%njobs,master)
    call ls_mpi_buffer(jobs%nocc,jobs%njobs,master)
    call ls_mpi_buffer(jobs%nvirt,jobs%njobs,master)
    call ls_mpi_buffer(jobs%nbasis,jobs%njobs,master)
    call ls_mpi_buffer(jobs%ntasks,jobs%njobs,master)
    call ls_mpi_buffer(jobs%flops,jobs%njobs,master)
    call ls_mpi_buffer(jobs%LMtime,jobs%njobs,master)
    call ls_mpi_buffer(jobs%load,jobs%njobs,master)

  end subroutine mpicopy_fragment_joblist


  !> \brief Divide local group in two smaller groups
  !> each of which is half the original group.
  !> The global parameters infpar%lg_mynum, infpar%lg_nodtot, and infpar%lg_comm
  !> will be redefined here.
  !> It is important to ensure that ALL MEMBERS of the local group call
  !> this routine at the same time - if not, some ranks will be waiting for a signal
  !> which never comes...
  !> In other words, synchronization of all nodes in the local groups must be done OUTSIDE this
  !> routine because their communication channel is redefined in here!
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine dec_half_local_group
    implicit none
    integer(kind=ls_mpik) :: ngroups
    integer(kind=ls_mpik) :: groupdims(2)

    ! Sanity check: Exit if local group is only of size 1
    if(infpar%lg_nodtot==1) return

    ! Number of ranks in new group 1:
    ! Half the number of ranks in current local group (plus one if uneven)
    groupdims(1) = ceiling(real(infpar%lg_nodtot)/2.0)

    ! Number of ranks in new group 2: The rest
    groupdims(2) = infpar%lg_nodtot - groupdims(1)

    ! Create new local groups:
    ! Current values of infpar%lg_mynum, infpar%lg_nodtot, and infpar%lg_comm
    ! will be overwritten.
    ngroups=2
    call divide_local_mpi_group(ngroups,groupdims)

  end subroutine dec_half_local_group


  !> \brief MPI copy information in mp2dens type.
  !> The same practical comments (note 1 and note 2) in mpicopy_fragment
  !> applies here.
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine mpicopy_fragmentdensity(dens,only_update)

    implicit none

    !> MP2 fragment density matrix information
    type(mp2dens),intent(inout) :: dens
    !> Only MPI copy information needed for update_full_mp2gradient
    !> (this will minimize MPI communication for global master)
    logical,intent(in) :: only_update


    ! INFORMATION NEEDED BY update_full_mp2gradient
    ! *********************************************
    CALL ls_mpi_buffer(dens%nbasis,infpar%master)
    if(.not. AddToBuffer) then
       nullify(dens%basis_idx)
       call mem_alloc(dens%basis_idx,dens%nbasis)
       nullify(dens%rho)
       call mem_alloc(dens%rho,dens%nbasis,dens%nbasis)
    end if
    ! Buffer handling
    call ls_mpi_buffer(dens%basis_idx,dens%nbasis,infpar%master)
    call ls_mpi_buffer(dens%rho,dens%nbasis,dens%nbasis,infpar%master)
    call ls_mpi_buffer(dens%energy,infpar%master)

    if(only_update) then  ! We do not want to sent the remaining (redundant) information
       return
    end if


    ! Integers that are not pointers
    CALL ls_mpi_buffer(dens%CentralAtom,infpar%master)
    CALL ls_mpi_buffer(dens%CentralAtom2,infpar%master)
    CALL ls_mpi_buffer(dens%nvirt,infpar%master)
    CALL ls_mpi_buffer(dens%nocc,infpar%master)
    CALL ls_mpi_buffer(dens%nocctot,infpar%master)
    CALL ls_mpi_buffer(dens%nEOSatoms,infpar%master)

    ! Reals that are not pointers
    call ls_mpi_buffer(dens%pairdist,infpar%master)


    ! Integer pointers
    ! ----------------
    ! Nullify and allocate stuff for receiver (global addtobuffer is false)
    if(.not. AddToBuffer) then
       nullify(dens%EOSatoms)
       call mem_alloc(dens%EOSatoms,dens%nEOSatoms)
    end if
    ! Buffer handling
    call ls_mpi_buffer(dens%EOSatoms,dens%nEOSatoms,infpar%master)


    ! Real pointers
    ! -------------
    if(.not. AddToBuffer) then
       nullify(dens%Y)
       call mem_alloc(dens%Y,dens%nvirt,dens%nvirt)
       nullify(dens%X)
       call mem_alloc(dens%X,dens%nocc,dens%nocc)
       nullify(dens%Phivo)
       call mem_alloc(dens%Phivo,dens%nvirt,dens%nocctot)
       nullify(dens%Phiov)
       call mem_alloc(dens%Phiov,dens%nocc,dens%nvirt)
    end if
    ! Buffer handling
    call ls_mpi_buffer(dens%Y,dens%nvirt,dens%nvirt,infpar%master)
    call ls_mpi_buffer(dens%X,dens%nocc,dens%nocc,infpar%master)
    call ls_mpi_buffer(dens%Phivo,dens%nvirt,dens%nocctot,infpar%master)
    call ls_mpi_buffer(dens%Phiov,dens%nocc,dens%nvirt,infpar%master)


  end subroutine mpicopy_fragmentdensity




  !> \brief MPI copy information in mp2grad type.
  !> The same practical comments (note 1 and note 2) in mpicopy_fragment
  !> applies here.
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine mpicopy_fragmentgradient(grad,only_update)

    implicit none

    !> MP2 fragment gradient information
    type(mp2grad),intent(inout) :: grad
    !> Only MPI copy information needed for update_full_mp2gradient
    !> (this will minimize MPI communication for global master)
    logical,intent(in) :: only_update
    logical :: only_update_check


    ! Sanity check
    only_update_check = only_update
    call ls_mpi_buffer(only_update_check,infpar%master)
    if(only_update_check .neqv. only_update) then
       call lsquit('mpicopy_fragmentgradient: Mismatch in only_update input!',-1)
    end if


    ! Copy density part of gradient structure
    ! ***************************************
    call mpicopy_fragmentdensity(grad%dens,only_update)


    ! Information needed by update_full_mp2gradient
    ! *********************************************
    call ls_mpi_buffer(grad%natoms,infpar%master)


    if(.not. AddToBuffer) then
       nullify(grad%atoms_idx)
       call mem_alloc(grad%atoms_idx,grad%natoms)
       nullify(grad%PhiAO)
       call mem_alloc(grad%PhiAO,grad%dens%nbasis,grad%dens%nbasis)
       nullify(grad%Ltheta)
       call mem_alloc(grad%Ltheta,3,grad%natoms)
    end if
    call ls_mpi_buffer(grad%atoms_idx, grad%natoms, infpar%master)
    call ls_mpi_buffer(grad%PhiAO, grad%dens%nbasis, grad%dens%nbasis, infpar%master)
    call ls_mpi_buffer(grad%Ltheta, 3, grad%natoms, infpar%master)

    if(only_update) then  ! do not sent remaining (redundant) information
       return
    end if



    ! Information not needed by update_full_mp2gradient
    ! *************************************************

    ! Pointers
    ! --------
    ! Nullify and allocate stuff for receiver (global addtobuffer is false)
    if(.not. AddToBuffer) then
       nullify(grad%Phioo)
       call mem_alloc(grad%Phioo,grad%dens%nocc,grad%dens%nocctot)
       nullify(grad%Phivv)
       call mem_alloc(grad%Phivv,grad%dens%nvirt,grad%dens%nvirt)
    end if
    ! Buffer handling
    call ls_mpi_buffer(grad%Phioo, grad%dens%nocc,grad%dens%nocctot, infpar%master)
    call ls_mpi_buffer(grad%Phivv, grad%dens%nvirt,grad%dens%nvirt, infpar%master)

  end subroutine mpicopy_fragmentgradient



  !> \brief Print MPI fragment statistics (see type joblist)
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine print_MPI_fragment_statistics(jobs,mastertime,string)

    implicit none
    !> Fragment job list
    type(joblist),intent(in) :: jobs
    !> Time measured by global master when fragments were calculated
    real(realk),intent(in) :: mastertime
    !> String to print describing job list
    character(*),intent(in) :: string
    real(realk) :: Gflops, avflop, minflop, maxflop,tmp,totflops,tottime_actual,localloss
    real(realk) :: globalloss, tottime_ideal, slavetime, localuse
    integer :: i, minidx, maxidx,N


    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') '*********************************************************************'
    write(DECinfo%output,'(1X,2a)') ' MPI fragment statistics: ', string
    write(DECinfo%output,'(1X,a)') '*********************************************************************'

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) 'Printed information'
    write(DECinfo%output,*) '-------------------'
    write(DECinfo%output,*) 'Job    : Fragment job number'
    write(DECinfo%output,*) '#occ   : Number of occupied AOS orbitals in fragment'
    write(DECinfo%output,*) '#virt  : Number of virtual AOS orbitals in fragment'
    write(DECinfo%output,*) '#basis : Number of basis functions in fragment'
    write(DECinfo%output,*) 'slotsiz: Number of slaves used for fragment incl. local master (slot size)'
    write(DECinfo%output,*) '#tasks : Number of integral tasks for fragment (nalpha*ngamma)'
    write(DECinfo%output,*) 'GFLOPS : Accumulated GFLOPS for fragment (local master AND local slaves)'
    write(DECinfo%output,*) 'Time(s): Time (in seconds) used by local master (NOT local slaves)'
    write(DECinfo%output,*) 'Load   : Load distribution measure (ideally 1.0, smaller in practice): '
    write(DECinfo%output,*) '           (Sum of effective work times for ALL MPI processes in slot)'
    write(DECinfo%output,*) '         / (slotsiz * [Local master time] )'
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(5X,a,4X,a,3X,a,2X,a,1X,a,2X,a,5X,a,5X,a,5X,a)') 'Job', '#occ', &
         & '#virt', '#basis', 'slotsiz', '#tasks', 'GFLOPS', 'Time(s)', 'Load'
    avflop=0.0E0_realk
    totflops=0.0E0_realk
    tottime_actual = 0.0E0_realk
    slavetime= 0.0_realk

    minflop = huge(1.0)
    maxflop=tiny(1.0)
    minidx=0
    maxidx=0
    N=0
    do i=1,jobs%njobs
       ! If nocc is zero, the job was not done and we do not print it
       if(jobs%nocc(i)==0) cycle
       N=N+1

       ! Giga flops for fragment
       Gflops = jobs%flops(i)*1.0e-9
       ! Update total number of flops
       totflops = totflops + Gflops
       ! Update total time used by ALL nodes (including dead time by local slaves)
       tottime_actual = tottime_actual + jobs%LMtime(i)*jobs%nslaves(i)
       ! Effective slave time (WITHOUT dead time by slaves)
       slavetime = slavetime + jobs%load(i)*jobs%nslaves(i)*jobs%LMtime(i)

       write(DECinfo%output,'(6i8,3X,3g11.3,a)') i, jobs%nocc(i), jobs%nvirt(i), jobs%nbasis(i),&
            & jobs%nslaves(i), jobs%ntasks(i), Gflops, jobs%LMtime(i), jobs%load(i), 'STAT'


       ! Accumulated Gflops per sec
       tmp = Gflops/jobs%LMtime(i)
       ! Gflops per sec per node
       tmp = tmp/real(jobs%nslaves(i))

       ! Update min and max flop per node per sec
       if(minflop>tmp) then
          minflop = tmp
          minidx=i
       end if
       if(maxflop<tmp) then
          maxflop = tmp
          maxidx=i
       end if

    end do


    ! Ideal total time used by all slaves:
    ! Assuming all infpar%nodtot-1 slaves were doing something useful all the time 
    ! while the master were sending out jobs.
    tottime_ideal = (infpar%nodtot-1)*mastertime
    ! If this time is zero (debugging case of restart case - then don't print summary

    if(tottime_ideal > 1.0e-6_realk) then
       ! Average % loss in local groups
       ! Time lost in local groups is the actual local time minus the effective slave time
       localloss = ( ( tottime_actual - slavetime ) / tottime_ideal )*100.0_realk

       ! Average % global loss (idle time at the end when some slots wait for others)
       globalloss = ( 1 - (tottime_actual / tottime_ideal) )*100.0_realk

       ! Average flops per node per sec
       avflop = totflops/tottime_actual
       write(DECinfo%output,*) '-----------------------------------------------------------------------------'
       write(DECinfo%output,'(1X,a,i8)')    'Number of jobs done     = ', N
       write(DECinfo%output,'(1X,a,g12.3)') 'Time-to-solution (s)    = ', mastertime
       write(DECinfo%output,'(1X,a,g12.3)') 'TOTAL time(s) all nodes = ', tottime_actual
#ifdef VAR_PAPI
       write(DECinfo%output,'(1X,a,g12.3)') 'TOTAL Gflops  = ', totflops
       write(DECinfo%output,'(1X,a,g12.3)') 'AVERAGE Gflops/s per MPI process= ', avflop
       write(DECinfo%output,'(1X,a,g12.3,a,i8)') 'MINIMUM Gflops/s per MPI process = ', &
            & minflop, ' for job ', minidx
       write(DECinfo%output,'(1X,a,g12.3,a,i8)') 'MAXIMUM Gflops/s per MPI process = ', &
            & maxflop, ' for job ', maxidx
#endif
       write(DECinfo%output,'(1X,a,g12.3)') 'Local MPI loss (%)  = ', localloss
       write(DECinfo%output,'(1X,a,g12.3)') 'Global MPI loss (%) = ', globalloss
       write(DECinfo%output,'(1X,a,g12.3)') 'Total MPI loss (%)  = ', localloss+globalloss
       write(DECinfo%output,*) '-----------------------------------------------------------------------------'

    end if
    write(DECinfo%output,*)
    write(DECinfo%output,*)

  end subroutine print_MPI_fragment_statistics


  subroutine mpicopy_dec_settings(decitem,comm)
    implicit none
    integer(kind=ls_mpik) :: comm
    type(decsettings) :: DECitem

    call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,comm)
    call mpicopy_dec_settings1(decitem,comm)
    call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,comm)

  end subroutine mpicopy_dec_settings

  ! USE HIGHLY DISCOURAGED!!!
  subroutine mpicopy_dec_settings1(decitem,comm)
    implicit none
    integer(kind=ls_mpik) :: comm,master
    type(decsettings) :: DECitem
    master = 0
    call ls_mpi_buffer(DECitem%doDEC,Master)
    call ls_mpi_buffer(DECitem%doHF,Master)
    call ls_mpi_buffer(DECitem%memory,Master)
    call ls_mpi_buffer(DECitem%full_molecular_cc,Master)
    call ls_mpi_buffer(DECitem%FullDEC,Master)
    call ls_mpi_buffer(DECitem%simulate_full,Master)
    call ls_mpi_buffer(DECitem%simulate_natoms,Master)
    call ls_mpi_buffer(DECitem%SkipReadIn,Master)
    call ls_mpi_buffer(DECitem%SinglesPolari,Master)
    call ls_mpi_buffer(DECitem%singlesthr,Master)
    call ls_mpi_buffer(DECitem%nsubfrag,Master)
    call ls_mpi_buffer(DECitem%nsubfrag2,Master)
    call ls_mpi_buffer(DECitem%subfrag,Master)
    call ls_mpi_buffer(DECitem%restart,Master)
    call ls_mpi_buffer(DECitem%TimeBackup,Master)
    call ls_mpi_buffer(DECitem%fragmentation_debug,Master)
    call ls_mpi_buffer(DECitem%dec_driver_debug,Master)
    call ls_mpi_buffer(DECitem%cc_driver_debug,Master)
    call ls_mpi_buffer(DECitem%ccsd_old,Master)
    call ls_mpi_buffer(DECitem%manual_batchsizes,Master)
    call ls_mpi_buffer(DECitem%ccsdAbatch,Master)
    call ls_mpi_buffer(DECitem%ccsdGbatch,Master)
    call ls_mpi_buffer(DECitem%hack,Master)
    call ls_mpi_buffer(DECitem%hack2,Master)
    call ls_mpi_buffer(DECitem%mpidebug,Master)
    call ls_mpi_buffer(DECitem%output,Master)
    call ls_mpi_buffer(DECitem%mulliken_threshold,Master)
    call ls_mpi_buffer(DECitem%simple_mulliken_threshold,Master)
    call ls_mpi_buffer(DECitem%approximated_norm_threshold,Master)
    call ls_mpi_buffer(DECitem%check_lcm_orbitals,Master)
    call ls_mpi_buffer(DECitem%use_canonical,Master)
    call ls_mpi_buffer(DECitem%reassignHatoms,Master)
    call ls_mpi_buffer(DECitem%mulliken,Master)
    call ls_mpi_buffer(DECitem%BoughtonPulay,Master)
    call ls_mpi_buffer(DECitem%FitOrbitals,Master)
    call ls_mpi_buffer(DECitem%simple_orbital_threshold,Master)
    call ls_mpi_buffer(DECitem%MaxIter,Master)
    call ls_mpi_buffer(DECitem%FOT,Master)
    call ls_mpi_buffer(DECitem%PL,Master)
    call ls_mpi_buffer(DECitem%SkipCC,Master)
    call ls_mpi_buffer(DECitem%purifyMOs,Master)
    call ls_mpi_buffer(DECitem%precondition_with_full,Master)
    call ls_mpi_buffer(DECitem%InclFullMolecule,Master)
    call ls_mpi_buffer(DECitem%HybridScheme,Master)
    call ls_mpi_buffer(DECitem%LagStepSize,Master)
    call ls_mpi_buffer(DECitem%pair_distance_threshold,Master)
    call ls_mpi_buffer(DECitem%PairReductionDistance,Master)
    call ls_mpi_buffer(DECitem%fragment_to_do,Master)
    call ls_mpi_buffer(DECitem%fragment_a_to_pair,Master)
    call ls_mpi_buffer(DECitem%fragment_b_to_pair,Master)
    call ls_mpi_buffer(DECitem%fullmolecule_memory,Master)
    call ls_mpi_buffer(DECitem%ccsd_expl,Master)
    DECitem%cc_models(1)='MP2     '
    DECitem%cc_models(2)='CC2     '
    DECitem%cc_models(3)='CCSD    '
    DECitem%cc_models(4)='CCSD(T) '
    call ls_mpi_buffer(DECitem%t2_restart,Master)
    call ls_mpi_buffer(DECitem%simulate_eri,Master)
    call ls_mpi_buffer(DECitem%fock_with_ri,Master)
    call ls_mpi_buffer(DECitem%ccMaxIter,Master)
    call ls_mpi_buffer(DECitem%ccMaxDIIS,Master)
    call ls_mpi_buffer(DECitem%ccModel,Master)
    call ls_mpi_buffer(DECitem%F12,Master)
    call ls_mpi_buffer(DECitem%ccConvergenceThreshold,Master)
    call ls_mpi_buffer(DECitem%CCthrSpecified,Master)
    call ls_mpi_buffer(DECitem%use_singles,Master)
    call ls_mpi_buffer(DECitem%use_preconditioner,Master)
    call ls_mpi_buffer(DECitem%use_preconditioner_in_b,Master)
    call ls_mpi_buffer(DECitem%use_crop,Master)
    call ls_mpi_buffer(DECitem%show_time,Master)
    call ls_mpi_buffer(DECitem%timing,Master)
    call ls_mpi_buffer(DECitem%show_memory,Master)
    call ls_mpi_buffer(DECitem%skip_full_ao,Master)
    call ls_mpi_buffer(DECitem%array4OnFile,Master)
    call ls_mpi_buffer(DECitem%array4OnFile_specified,Master)
    call ls_mpi_buffer(DECitem%AObasedCC,Master)
    call ls_mpi_buffer(DECitem%zero_threshold,Master)
    call ls_mpi_buffer(DECitem%integral_time_cpu,Master)
    call ls_mpi_buffer(DECitem%integral_time_wall,Master)
    call ls_mpi_buffer(DECitem%MOintegral_time_cpu,Master)
    call ls_mpi_buffer(DECitem%MOintegral_time_wall,Master)
    call ls_mpi_buffer(DECitem%solver_time_cpu,Master)
    call ls_mpi_buffer(DECitem%solver_time_wall,Master)
    call ls_mpi_buffer(DECitem%energy_time_cpu,Master)
    call ls_mpi_buffer(DECitem%energy_time_wall,Master)
    call ls_mpi_buffer(DECitem%density_time_cpu,Master)
    call ls_mpi_buffer(DECitem%density_time_wall,Master)
    call ls_mpi_buffer(DECitem%gradient_time_cpu,Master)
    call ls_mpi_buffer(DECitem%gradient_time_wall,Master)
    call ls_mpi_buffer(DECitem%trans_time_cpu,Master)
    call ls_mpi_buffer(DECitem%trans_time_wall,Master)
    call ls_mpi_buffer(DECitem%reorder_time_cpu,Master)
    call ls_mpi_buffer(DECitem%reorder_time_wall,Master)
    call ls_mpi_buffer(DECitem%memallo_time_cpu,Master)
    call ls_mpi_buffer(DECitem%memallo_time_wall,Master)
    call ls_mpi_buffer(DECitem%MPIgroupsize,Master)
    call ls_mpi_buffer(DECitem%first_order,Master)
    call ls_mpi_buffer(DECitem%MP2density,Master)
    call ls_mpi_buffer(DECitem%gradient,Master)
    call ls_mpi_buffer(DECitem%kappa_use_preconditioner,Master)
    call ls_mpi_buffer(DECitem%kappa_use_preconditioner_in_b,Master)
    call ls_mpi_buffer(DECitem%kappaMaxDIIS,Master)
    call ls_mpi_buffer(DECitem%kappaMaxIter,Master)
    call ls_mpi_buffer(DECitem%kappa_driver_debug,Master)
    call ls_mpi_buffer(DECitem%kappaTHR,Master)
    call ls_mpi_buffer(DECitem%SaveFragFile,Master)
  end subroutine mpicopy_dec_settings1

end module decmpi_module

#else
module decmpi_module

contains

!Added to avoid "has no symbols" linking warning
subroutine decmpi_module_void()
end subroutine decmpi_module_void

end module decmpi_module

#endif
