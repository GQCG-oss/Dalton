!> @file
!> DEC MPI handling
!> \author Kasper Kristensen
!> \date March 2012
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
  use background_buffer_module
  use tensor_basic_module
  use tensor_interface_module
!
  use DEC_settings_mod
  use dec_fragment_utils
  use array4_simple_operations
  use array2_simple_operations
  use CABS_operations

  private
#ifdef VAR_MPI
  public :: mpi_send_recv_fragmentenergy, mpi_send_recv_mp2grad,&
       & mpi_send_recv_single_fragment, mpi_send_recv_many_fragments,&
       & mpi_bcast_many_fragments, mpi_communicate_mp2_int_and_amp,&
       & mpi_communicate_MyFragment, decmpi_bcast_f12_info,&
       & mpi_dec_fullinfo_master_to_slaves, &
       & mpi_bcast_fullmolecule, mpicopy_fragment,&
       & mpicopy_fragmentAOStype,&
       & buffercopy_PNOSpaceInfo_struct,&
       & share_E2_with_slaves, mpi_communicate_ccsd_calcdata,&
       & mpi_communicate_ccsdpt_calcdata, get_mo_ccsd_joblist,&
       & distribute_mpi_jobs, bcast_dec_fragment_joblist,&
       & mpicopy_fragment_joblist, dec_half_local_group,&
       & mpicopy_fragmentdensity, mpicopy_fragmentgradient,&
       & print_MPI_fragment_statistics,&
       & mpibcast_dec_settings, mpi_communicate_get_gmo_data,&
       & mpi_communicate_moccsd_data,&
       & mpicopy_dec_settings, rpa_residual_communicate_data,&
       & rpa_fock_communicate_data, &
       & mpi_dec_fullinfo_master_to_slaves_precursor,&
       & wake_slaves_for_simple_mo,&
       & get_slaves_to_simple_par_mp2_res,&
       & check_job, &
       & mpi_communicate_MyFragment_f12
#else
public :: check_job
#endif


contains
 
#ifdef VAR_MPI

  !> \brief Send three fragment energies (for occupied, virtual, and Lagrangian schemes)
  !> from given sender (typically a slave) to given receiver (typically the master).
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
    !> Fragment energies, see "energy" in decfrag type def
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
     !> Fragment energies, see "energy" in decfrag type def
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
  !> NOTE: Only the information OUTSIDE the expensive box in the decfrag
  !> type definition can be sent this way. In general, the information outside
  !> the expensive box combined with full molecular information is enough
  !> to detemine the stuff outside the expensive box.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine mpi_send_recv_single_fragment(MyFragment,comm,MySender,MyReceiver,job)

    implicit none

    !> Fragment under consideration
    type(decfrag),intent(inout) :: MyFragment
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

    ! Sanity check 2: Basis info in "Expensive box" in decfrag should not be sent
    if(mynum==MySender) then
       if(MyFragment%BasisInfoIsSet) then
          call lsquit('mpi_send_recv_single_fragment: Not implemented for &
               & sending/receiving fragment basis info!',-1)
       end if
    endif
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
  !> NOTE: Only the information OUTSIDE the expensive box in the decfrag
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
    type(decfrag),dimension(natoms),intent(inout) :: Fragments
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

          ! Sanity check: Basis info in "Expensive box" in decfrag should not be sent
          if(Fragments(atom)%BasisInfoIsSet) then
             call lsquit('mpi_send_recv_many_fragments: Not implemented for &
                  & sending/receiving fragment basis info!',-1)
          end if

          ! Sender: Copy fragment info into buffer
          ! Receiver: Read fragment info from buffer into Fragments(atom)
          call mpicopy_fragment(Fragments(atom),comm,.false.)

          ! Check that the correct fragment number was communicated
          if(Fragments(atom)%EOSatoms(1) /= atom) then
             print '(a,2i8)', 'Counter atom, actual atom ', atom,Fragments(atom)%EOSatoms(1)
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
  !> NOTE: Only the information OUTSIDE the expensive box in the decfrag
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
    type(decfrag),dimension(natoms),intent(inout) :: Fragments
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

          ! Sanity check: Basis info in "Expensive box" in decfrag should not be sent
          if(mynum==0) then
             if(Fragments(atom)%BasisInfoIsSet) then
                call lsquit('mpi_bcast_many_fragments: Not implemented for &
                     & sending/receiving fragment basis info!',-1)
             end if
          endif
          ! Master: Copy fragment info into buffer
          ! Slave: Read fragment info from buffer into Fragments(atom)
          call mpicopy_fragment(Fragments(atom),comm,.false.)

          ! Check that the correct fragment number was communicated
          if(Fragments(atom)%EOSatoms(1) /= atom) then
             print '(a,2i8)', 'Counter atom, actual atom ', atom,Fragments(atom)%EOSatoms(1)
             call lsquit('mpi_bcast_many_fragments: Inconsistency between &
                  & counter and actual fragment',-1)
          end if

       end if CommunicateFragment

    end do

    ! Finalize buffer
    call ls_mpiFinalizeBuffer(master,LSMPIBROADCAST,comm)

  end subroutine mpi_bcast_many_fragments


  !> \brief MPI communcation where the information in the decfrag type
  !> is send (for local master) or received (for local slaves).
  !> In addition, other information required for calculating MP2 integrals
  !> and amplitudes is communicated.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine mpi_communicate_mp2_int_and_amp(MyFragment,bat,first_order_integrals,DoBasis)

    implicit none

    !> Fragment under consideration
    type(decfrag),intent(inout) :: MyFragment
    !> Batch information
    type(mp2_batch_construction),intent(inout) :: bat
    !> Do first order integrals?
    logical,intent(inout) :: first_order_integrals
    !> Communicate basis (expensive box in decfrag)
    logical,intent(in) :: DoBasis
    integer(kind=ls_mpik) :: master

    call time_start_phase( PHASE_COMM )

    master = 0

    ! Init MPI buffer which eventually will contain all fragment info
    ! ***************************************************************
    ! MASTER: Prepare for writing to buffer
    ! SLAVE: Receive buffer
    call ls_mpiInitBuffer(master,LSMPIBROADCAST,infpar%lg_comm)


    ! Buffer handling
    ! ***************
    ! MASTER: Put information info buffer
    ! SLAVE: Put buffer information into decfrag structure

    ! Fragment information
    call mpicopy_fragment(MyFragment,infpar%lg_comm,DoBasis)

    ! Logical determining whether to calculate integrals for first order properties
    call ls_mpi_buffer(first_order_integrals,master)

    ! Batch information
    call ls_mpi_buffer(bat%MaxAllowedDimAlpha,master)
    call ls_mpi_buffer(bat%MaxAllowedDimGamma,master)
    call ls_mpi_buffer(bat%virtbatch,master)

    call ls_mpi_buffer(bat%size1,4,master)
    call ls_mpi_buffer(bat%size2,4,master)
    call ls_mpi_buffer(bat%size3,4,master)

    ! Finalize MPI buffer
    ! *******************
    ! MASTER: Send stuff to slaves and deallocate temp. buffers
    ! SLAVE: Deallocate buffer etc.
    call ls_mpiFinalizeBuffer(master,LSMPIBROADCAST,infpar%lg_comm)

    call time_start_phase( PHASE_WORK )
  end subroutine mpi_communicate_mp2_int_and_amp

  !> \brief MPI communcation of MyFragment
  !> \author Thomas Kjaergaard
  !> \date March 2015
  subroutine mpi_communicate_MyFragment(MyFragment)
    implicit none
    !> Fragment under consideration
    type(decfrag),intent(inout) :: MyFragment
    !local variables
    integer(kind=ls_mpik) :: master
    logical :: DoBasis
    DoBasis = .TRUE.
    
    call time_start_phase( PHASE_COMM )
    master = 0
    ! Init MPI buffer which eventually will contain all fragment info
    ! ***************************************************************
    ! MASTER: Prepare for writing to buffer
    ! SLAVE: Receive buffer
    call ls_mpiInitBuffer(master,LSMPIBROADCAST,infpar%lg_comm)

    ! Buffer handling
    ! ***************
    ! MASTER: Put information info buffer
    ! SLAVE: Put buffer information into decfrag structure
    ! Fragment information
    call mpicopy_fragment(MyFragment,infpar%lg_comm,DoBasis)

    ! Finalize MPI buffer
    ! *******************
    ! MASTER: Send stuff to slaves and deallocate temp. buffers
    ! SLAVE: Deallocate buffer etc.
   call ls_mpiFinalizeBuffer(master,LSMPIBROADCAST,infpar%lg_comm)

    call time_start_phase( PHASE_WORK )
  end subroutine mpi_communicate_MyFragment 

  subroutine mpi_communicate_MyFragment_f12(MyFragment,Taibj,fragcase,dopair)
    implicit none
    !> Fragment under consideration
    type(decfrag),intent(inout) :: MyFragment
    !> Case MODEL
    integer :: fragcase
    !> dopair_occ array
    logical,intent(inout) :: dopair

    real(realk),pointer,intent(in) :: Taibj(:,:,:,:)

    !local variables
    integer(kind=ls_mpik) :: master
    integer :: nrow,ncol
    logical :: DoBasis
    DoBasis = .TRUE.
    call time_start_phase( PHASE_COMM )
    master = 0
    ! Init MPI buffer which eventually will contain all fragment info
    ! ***************************************************************
    ! MASTER: Prepare for writing to buffer
    ! SLAVE: Receive buffer
    call ls_mpiInitBuffer(master,LSMPIBROADCAST,infpar%lg_comm)

    ! Buffer handling
    ! ***************
    ! MASTER: Put information info buffer
    ! SLAVE: Put buffer information into decfrag structure
    ! Fragment information
    call mpicopy_fragment(MyFragment,infpar%lg_comm,DoBasis)
 
    if(infpar%lg_mynum .ne. master) then
       call mem_alloc(MyFragment%Ccabs,MyFragment%ncabsAO,MyFragment%ncabsMO)
       call mem_alloc(MyFragment%Cri,MyFragment%ncabsAO,MyFragment%ncabsAO)
       call mem_alloc(MyFragment%Fcp,MyFragment%ncabsMO,MyFragment%nvirtAOS+MyFragment%nocctot)
       call mem_alloc(MyFragment%hJir,MyFragment%noccEOS,MyFragment%ncabsAO)
       call mem_alloc(MyFragment%Krs,MyFragment%ncabsAO,MyFragment%ncabsAO)                                                    
       call mem_alloc(MyFragment%Frs,MyFragment%ncabsAO,MyFragment%ncabsAO)                                                    
       call mem_alloc(MyFragment%Frm,MyFragment%ncabsAO,MyFragment%nocctot)
       call mem_alloc(MyFragment%Fcp,MyFragment%ncabsMO,MyFragment%nvirtAOS+MyFragment%nocctot)   
       call mem_alloc(Taibj,MyFragment%nvirtAOS,MyFragment%noccEOS,MyFragment%nvirtAOS,MyFragment%noccEOS)

    endif

    call ls_mpi_buffer(MyFragment%Ccabs,MyFragment%ncabsAO,MyFragment%ncabsMO,master)
    call ls_mpi_buffer(MyFragment%Cri,MyFragment%ncabsAO,MyFragment%ncabsAO,master)
    call ls_mpi_buffer(MyFragment%Fcp,MyFragment%ncabsMO,MyFragment%nvirtAOS+MyFragment%nocctot,master) 
    call ls_mpi_buffer(MyFragment%hJir,MyFragment%noccEOS,MyFragment%ncabsAO,master) 
    call ls_mpi_buffer(MyFragment%Krs,MyFragment%ncabsAO,MyFragment%ncabsAO,master)                                             
    call ls_mpi_buffer(MyFragment%Frs,MyFragment%ncabsAO,MyFragment%ncabsAO,master)                                             
    call ls_mpi_buffer(MyFragment%Frm,MyFragment%ncabsAO,MyFragment%nocctot,master)
    call ls_mpi_buffer(MyFragment%Fcp,MyFragment%ncabsMO,MyFragment%nvirtAOS+MyFragment%nocctot,master)     
    call ls_mpi_buffer(Taibj,MyFragment%nvirtAOS,MyFragment%noccEOS,MyFragment%nvirtAOS,MyFragment%noccEOS,master)


    call ls_mpi_buffer(fragcase,master) 
    call ls_mpi_buffer(dopair,master)

    ! Finalize MPI buffer
    ! *******************
    ! MASTER: Send stuff to slaves and deallocate temp. buffers
    ! SLAVE: Deallocate buffer etc.
    call ls_mpiFinalizeBuffer(master,LSMPIBROADCAST,infpar%lg_comm)
    call time_start_phase( PHASE_WORK )
 end subroutine mpi_communicate_MyFragment_f12 


  !> \brief MPI communcation where the information in the decfrag type
  !> is send (for local master) or received (for local slaves).
  !> In addition, other information required for calculating F12 integrals
  !> and amplitudes is communicated.
  !> \author Yang Min Wang
  !> \date 2014
  subroutine decmpi_bcast_f12_info(MyFragment, Taibj, Tai, case,PairFrag,DoBasis, Fragment1, Fragment2)
    implicit none
    !> Atomic fragment to be determined (Single or Pair fragment)
    type(decfrag) :: MyFragment
    !> t2EOS amplitudes stored in the order T(a,i,b,j)
    real(realk),pointer :: Taibj(:,:,:,:) 
    !> t1EOS amplitudes stored in the order T(a,i)
    real(realk),pointer :: Tai(:,:) 
    !> Case MODEL
    integer :: case
    !> Pair or Atomic Calc
    logical :: PairFrag
    !> Bcast the Expensive Box
    logical :: DoBasis
    !> Fragment 1 in the pair fragment
    type(decfrag),optional :: Fragment1
    !> Fragment 2 in the pair fragment
    type(decfrag),optional :: Fragment2
    !local variables
    integer(kind=ls_mpik) :: master
    logical :: Slave,DoBasisPair
    integer :: Tain1,Tain2,n1,n2,n3,n4
    master = 0
    Slave = infpar%lg_mynum .NE. master
    call time_start_phase( PHASE_COMM )
    ! Init MPI buffer which eventually will contain all fragment info
    ! ***************************************************************
    ! MASTER: Prepare for writing to buffer
    ! SLAVE: Receive buffer
    call ls_mpiInitBuffer(master,LSMPIBROADCAST,infpar%lg_comm)

    ! Buffer handling
    ! ***************
    ! MASTER: Put information info buffer
    ! SLAVE: Put buffer information into decfrag structure

    ! Fragment information
    call ls_mpi_buffer(DoBasis,master)
    call mpicopy_fragment(MyFragment,infpar%lg_comm,DoBasis)

    call ls_mpi_buffer(PairFrag,master)
    IF(PairFrag)THEN
       DoBasisPair = .FALSE. !this info is not needed for these fragments
       call mpicopy_fragment(Fragment1,infpar%lg_comm,DoBasisPair)
       call mpicopy_fragment(Fragment2,infpar%lg_comm,DoBasisPair)
    ENDIF
    call ls_mpi_buffer(case,master)

    ! Finalize MPI buffer
    ! *******************
    ! MASTER: Send stuff to slaves and deallocate temp. buffers
    ! SLAVE: Deallocate buffer etc.
    IF(.NOT.Slave)THEN
       n1 = size(Taibj,1)
       n2 = size(Taibj,2)
       n3 = size(Taibj,3)
       n4 = size(Taibj,4)
       IF(case.EQ.MODEL_CCSD)THEN
          Tain1 = size(Tai,1)
          Tain2 = size(Tai,2)
       ENDIF
    ENDIF
    call ls_mpi_buffer(n1,master)
    call ls_mpi_buffer(n2,master)
    call ls_mpi_buffer(n3,master)
    call ls_mpi_buffer(n4,master)
    IF(case.EQ.MODEL_CCSD)THEN
       call ls_mpi_buffer(Tain1,master)
       call ls_mpi_buffer(Tain2,master)
    ENDIF
    call ls_mpiFinalizeBuffer(master,LSMPIBROADCAST,infpar%lg_comm)

    IF(Slave)call mem_alloc(Taibj,n1,n2,n3,n4)
    call ls_mpibcast(Taibj,n1,n2,n3,n4,master,infpar%lg_comm)
    IF(case.EQ.MODEL_CCSD)THEN
       IF(Slave)call mem_alloc(Tai,Tain1,Tain2)
       call ls_mpibcast(Tai,Tain1,Tain2,master,infpar%lg_comm)
    ENDIF

    call time_start_phase( PHASE_WORK )
  end subroutine decmpi_bcast_f12_info


  !> \brief MPI communcation of full molecular information required by slaves
  !> to start fragment calculation.
  !> Assumes that basic dimension information (number of atoms and number of orbitals)
  !> have already been communicated.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine mpi_dec_fullinfo_master_to_slaves(nocc,nvirt,&
       & OccOrbitals, virtOrbitals, MyMolecule, MyLsitem)

    implicit none

    !> Number of occupied orbitals in the molecule
    integer,intent(in) :: nocc
    !> Number of virtupied orbitals in the molecule
    integer,intent(in) :: nvirt
    !> Occupied orbitals in DEC format
    !> Intent(in) for main master, intent(out) for local masters
    type(decorbital),intent(inout) :: OccOrbitals(nocc)
    !> virtupied orbitals in DEC format
    !> Intent(in) for main master, intent(out) for local masters
    type(decorbital),intent(inout) :: virtOrbitals(nvirt)
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

    ! Integers and reals inside decorbital sub-type
    do i=1,nocc
       call ls_mpi_buffer(OccOrbitals(i)%orbitalnumber,master)
       call ls_mpi_buffer(OccOrbitals(i)%centralatom,master)
       call ls_mpi_buffer(OccOrbitals(i)%secondaryatom,master)
       call ls_mpi_buffer(OccOrbitals(i)%numberofaos,master)

       ! Pointer for Atomic Orbitals (AOs) for orbital "i"
       if(.not. gm) then ! allocate for local masters (not global master)
          Nullify(OccOrbitals(i)%aos)
          call mem_alloc(OccOrbitals(i)%aos,OccOrbitals(i)%numberofaos)
       end if
       ! Buffer handling for both master and slave
       call ls_mpi_buffer(OccOrbitals(i)%aos, &
            & OccOrbitals(i)%numberofaos,master)
    end do


    ! DEC virtupied orbitals
    ! -----------------------

    ! Integers and reals inside decorbital sub-type
    do i=1,nvirt
       call ls_mpi_buffer(virtOrbitals(i)%orbitalnumber,master)
       call ls_mpi_buffer(virtOrbitals(i)%centralatom,master)
       call ls_mpi_buffer(virtOrbitals(i)%secondaryatom,master)
       call ls_mpi_buffer(virtOrbitals(i)%numberofaos,master)

       ! Pointer for Atomic Orbitals (AOs) for orbital "i"
       if(.not. gm) then ! allocate for local masters (not global master)
          Nullify(virtOrbitals(i)%aos)
          call mem_alloc(virtOrbitals(i)%aos,virtOrbitals(i)%numberofaos)
       end if
       ! Buffer handling for both master and slave
       call ls_mpi_buffer(virtOrbitals(i)%aos, &
            & virtOrbitals(i)%numberofaos,master)
    end do


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
    integer :: natoms2
    integer :: Co_addr(infpar%nodtot)
    integer :: Cv_addr(infpar%nodtot)
    integer :: fock_addr(infpar%nodtot)
    integer :: oofock_addr(infpar%nodtot)
    integer :: vvfock_addr(infpar%nodtot)

    master = 0

    IF(infpar%mynum==0)THEN ! Global master
       gm=.true.
    else ! global slave
       gm=.false.
    end IF

    ! Simple integers
    ! ---------------
    call ls_mpibcast(MyMolecule%nelectrons,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%natoms,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%nfrags,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%nbasis,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%nMO,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%nauxbasis,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%nocc,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%ncore,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%nval,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%nvirt,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%nCabsAO,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%nCabsAOOnly,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%nCabsMO,master,MPI_COMM_LSDALTON)

    !simple logicals
    call ls_mpibcast(MyMolecule%mem_distributed,master,MPI_COMM_LSDALTON)

    ! Simple reals
    call ls_mpibcast(MyMolecule%Edisp,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%Ect,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%Esub,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%EF12singles,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%EF12NLSV1,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%EF12NLSB1,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%EF12NLSX1,master,MPI_COMM_LSDALTON)

    ! Allocate pointers if global slave
    if(.not. gm) then
       call mem_alloc(MyMolecule%atom_size,MyMolecule%natoms)
       call mem_alloc(MyMolecule%atom_start,MyMolecule%natoms)
       call mem_alloc(MyMolecule%atom_end,MyMolecule%natoms)
       call mem_alloc(MyMolecule%bas_start,MyMolecule%nbasis)
       call mem_alloc(MyMolecule%bas_end,MyMolecule%nbasis)
       IF(DECINFO%F12)THEN
          call init_cabs()
          call mem_alloc(MyMolecule%atom_cabssize,MyMolecule%natoms)
          call mem_alloc(MyMolecule%atom_cabsstart,MyMolecule%natoms)
       ENDIF
       if(.not.MyMolecule%mem_distributed)then
          call tensor_init(MyMolecule%Co,[MyMolecule%nbasis,MyMolecule%nocc],2)
          call tensor_init(MyMolecule%Cv,[MyMolecule%nbasis,MyMolecule%nvirt],2)
          if(.not. DECinfo%noaofock) then
             call tensor_init(MyMolecule%fock,[MyMolecule%nbasis,MyMolecule%nbasis],2)
          end if
          call tensor_init(MyMolecule%oofock,[MyMolecule%nocc,MyMolecule%nocc],2)
          call tensor_init(MyMolecule%vvfock,[MyMolecule%nvirt,MyMolecule%nvirt],2)
       endif

       call mem_alloc(MyMolecule%carmomocc,3,MyMolecule%nocc)
       call mem_alloc(MyMolecule%carmomvirt,3,MyMolecule%nvirt)
       call mem_alloc(MyMolecule%AtomCenters,3,MyMolecule%natoms)
       call mem_alloc(MyMolecule%PhantomAtom,MyMolecule%natoms)
       IF(DECinfo%F12)THEN          
          IF(.NOT.DECinfo%full_molecular_cc)THEN
             call mem_alloc(MyMolecule%Fij,MyMolecule%nocc,MyMolecule%nocc)
             call mem_alloc(MyMolecule%hJir,MyMolecule%nocc,MyMolecule%nCabsAO)
             call mem_alloc(MyMolecule%Krs,MyMolecule%nCabsAO,MyMolecule%nCabsAO)
             call mem_alloc(MyMolecule%Frs,MyMolecule%nCabsAO,MyMolecule%nCabsAO)
             !HACK NOT MyMolecule%nvirt,MyMolecule%nCabsMO
             call mem_alloc(MyMolecule%Fac,MyMolecule%nvirt,MyMolecule%nCabsAO)
             !Warning MyMolecule%Frm is allocated with noccfull ?????
             call mem_alloc(MyMolecule%Frm,MyMolecule%nCabsAO,MyMolecule%nocc)
             !HACK NOT MyMolecule%nCabsMO,MyMolecule%nbasis
             call mem_alloc(MyMolecule%Fcp,MyMolecule%nCabsAO,MyMolecule%nbasis)
          ENDIF
       ENDIF
       call mem_alloc(MyMolecule%SubSystemIndex,MyMolecule%nAtoms)
       call mem_alloc(MyMolecule%DistanceTable,MyMolecule%nfrags,MyMolecule%nfrags)
       call mem_alloc(MyMolecule%ccmodel,MyMolecule%nfrags,MyMolecule%nfrags)
       call mem_alloc(MyMolecule%PairFOTlevel,MyMolecule%nfrags,MyMolecule%nfrags)
       if(DECinfo%use_abs_overlap)then
          call mem_alloc(MyMolecule%ov_abs_overlap,MyMolecule%nocc,MyMolecule%nvirt)
       endif
    end if


    ! Integer pointers
    ! ----------------
    call ls_mpibcast(MyMolecule%atom_size,MyMolecule%natoms,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%atom_start,MyMolecule%natoms,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%atom_end,MyMolecule%natoms,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%bas_start,MyMolecule%nbasis,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%bas_end,MyMolecule%nbasis,master,MPI_COMM_LSDALTON)
    IF(DECINFO%F12)THEN
       call ls_mpibcast(MyMolecule%atom_cabssize,MyMolecule%natoms,master,MPI_COMM_LSDALTON)
       call ls_mpibcast(MyMolecule%atom_cabsstart,MyMolecule%natoms,master,MPI_COMM_LSDALTON)
    ENDIF

    ! Real pointers
    ! -------------
    if(MyMolecule%mem_distributed)then

       !master get adresses
       if(gm)then
          Co_addr=MyMolecule%Co%addr_p_arr
          Cv_addr=MyMolecule%Cv%addr_p_arr
          if(.not. DECinfo%noaofock) then
             fock_addr=MyMolecule%fock%addr_p_arr
          endif
          oofock_addr=MyMolecule%oofock%addr_p_arr
          vvfock_addr=MyMolecule%vvfock%addr_p_arr
       endif

       !broadcast addresses
       call ls_mpibcast(Co_addr,infpar%lg_nodtot,infpar%master,MPI_COMM_LSDALTON)
       call ls_mpibcast(Cv_addr,infpar%lg_nodtot,infpar%master,MPI_COMM_LSDALTON)
       if(.not. DECinfo%noaofock) then
          call ls_mpibcast(fock_addr,infpar%lg_nodtot,infpar%master,MPI_COMM_LSDALTON)
       endif
       call ls_mpibcast(oofock_addr,infpar%lg_nodtot,infpar%master,MPI_COMM_LSDALTON)
       call ls_mpibcast(vvfock_addr,infpar%lg_nodtot,infpar%master,MPI_COMM_LSDALTON)

       !slaves set matrices
       if(.not.gm)then
          MyMolecule%Co     = get_tensor_from_parr(Co_addr(infpar%lg_mynum+1))
          MyMolecule%Cv     = get_tensor_from_parr(Cv_addr(infpar%lg_mynum+1))
          if(.not. DECinfo%noaofock) then
             MyMolecule%fock   = get_tensor_from_parr(fock_addr(infpar%lg_mynum+1))
          endif
          MyMolecule%oofock = get_tensor_from_parr(oofock_addr(infpar%lg_mynum+1))
          MyMolecule%vvfock = get_tensor_from_parr(vvfock_addr(infpar%lg_mynum+1))
       endif

    else

       call ls_mpibcast(MyMolecule%Co%elm2,MyMolecule%nbasis,MyMolecule%nocc,master,MPI_COMM_LSDALTON)
       call ls_mpibcast(MyMolecule%Cv%elm2,MyMolecule%nbasis,MyMolecule%nvirt,master,MPI_COMM_LSDALTON)
       if(.not. DECinfo%noaofock) then
          call ls_mpibcast(MyMolecule%fock%elm2,MyMolecule%nbasis,MyMolecule%nbasis,master,MPI_COMM_LSDALTON)
       end if
       call ls_mpibcast(MyMolecule%oofock%elm2,MyMolecule%nocc,MyMolecule%nocc,master,MPI_COMM_LSDALTON)
       call ls_mpibcast(MyMolecule%vvfock%elm2,MyMolecule%nvirt,MyMolecule%nvirt,master,MPI_COMM_LSDALTON)

    endif

    call ls_mpibcast(MyMolecule%carmomocc,3,MyMolecule%nocc,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%carmomvirt,3,MyMolecule%nvirt,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%AtomCenters,3,MyMolecule%natoms,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%PhantomAtom,MyMolecule%natoms,master,MPI_COMM_LSDALTON)
    IF(DECinfo%F12)THEN
       IF(.NOT.DECinfo%full_molecular_cc)THEN
          call ls_mpibcast(MyMolecule%Fij,MyMolecule%nocc,MyMolecule%nocc,master,MPI_COMM_LSDALTON)
          call ls_mpibcast(MyMolecule%hJir,MyMolecule%nocc,MyMolecule%nCabsAO,master,MPI_COMM_LSDALTON)
          call ls_mpibcast(MyMolecule%Krs,MyMolecule%nCabsAO,MyMolecule%nCabsAO,master,MPI_COMM_LSDALTON)
          call ls_mpibcast(MyMolecule%Frs,MyMolecule%nCabsAO,MyMolecule%nCabsAO,master,MPI_COMM_LSDALTON)
          !HACK NOT MyMolecule%nvirt,MyMolecule%nCabsMO
          call ls_mpibcast(MyMolecule%Fac,MyMolecule%nvirt,MyMolecule%nCabsAO,master,MPI_COMM_LSDALTON)
          !Warning MyMolecule%Frm is allocated with noccfull ?????
          call ls_mpibcast(MyMolecule%Frm,MyMolecule%nCabsAO,MyMolecule%nocc,master,MPI_COMM_LSDALTON)
          !HACK NOT MyMolecule%nCabsMO,MyMolecule%nbasis
          call ls_mpibcast(MyMolecule%Fcp,MyMolecule%nCabsAO,MyMolecule%nbasis,master,MPI_COMM_LSDALTON)
       ENDIF
    ENDIF
    call ls_mpibcast(MyMolecule%SubSystemIndex,MyMolecule%natoms,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%DistanceTable,MyMolecule%nfrags,MyMolecule%nfrags,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%ccmodel,MyMolecule%nfrags,MyMolecule%nfrags,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(MyMolecule%PairFOTlevel,MyMolecule%nfrags,MyMolecule%nfrags,&
         & master,MPI_COMM_LSDALTON)
    if(DECinfo%use_abs_overlap)then
       call ls_mpibcast(MyMolecule%ov_abs_overlap,MyMolecule%nocc,MyMolecule%nvirt,master,MPI_COMM_LSDALTON)
    endif

  end subroutine mpi_bcast_fullmolecule



  !> \brief MPI copy information in decfrag type.
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
    type(decfrag),intent(inout) :: MyFragment
    !> Communicator
    integer(kind=ls_mpik),intent(in) :: comm
    !> Copy basis (expensive box in decfrag)
    logical,intent(in) :: DoBasis
    integer :: i
    integer(kind=ls_mpik) :: master
    master = 0

    ! Integer pointers for some dimensions  
    if(.not. AddToBuffer) then
       call fragment_init_dimension_pointers(MyFragment)
    end if

    ! Integers that are not pointers
    ! ------------------------------

    CALL ls_mpi_buffer(MyFragment%noccEOS,master)
    CALL ls_mpi_buffer(MyFragment%nvirtEOS,master)
    CALL ls_mpi_buffer(MyFragment%ncore,master)
    CALL ls_mpi_buffer(MyFragment%nocctot,master)
    CALL ls_mpi_buffer(MyFragment%nEOSatoms,master)
    CALL ls_mpi_buffer(MyFragment%ntasks,master)
    CALL ls_mpi_buffer(MyFragment%t1dims,2,master)
    CALL ls_mpi_buffer(MyFragment%noccFA,master)
    CALL ls_mpi_buffer(MyFragment%nvirtFA,master)
    CALL ls_mpi_buffer(MyFragment%natoms,master)
    CALL ls_mpi_buffer(MyFragment%nbasis,master)
    CALL ls_mpi_buffer(MyFragment%nCabsAO,master)
    CALL ls_mpi_buffer(MyFragment%nCabsAOOnly,master)
    CALL ls_mpi_buffer(MyFragment%nCabsMO,master)
    CALL ls_mpi_buffer(MyFragment%ccmodel,master)
    CALL ls_mpi_buffer(MyFragment%noccLOC,master)
    CALL ls_mpi_buffer(MyFragment%nvirtLOC,master)
    CALL ls_mpi_buffer(MyFragment%nspaces,master)


    ! Logicals that are not pointers
    ! ------------------------------
    CALL ls_mpi_buffer(MyFragment%BasisInfoIsSet,master)
    CALL ls_mpi_buffer(MyFragment%t1_stored,master)
    CALL ls_mpi_buffer(MyFragment%CDset,master)
    CALL ls_mpi_buffer(MyFragment%FAset,master)
    CALL ls_mpi_buffer(MyFragment%PNOset,master)
    CALL ls_mpi_buffer(MyFragment%fragmentadapted,master)
    CALL ls_mpi_buffer(MyFragment%pairfrag,master)
    CALL ls_mpi_buffer(MyFragment%isopt,master)

    ! Reals that are not pointers
    ! ---------------------------
    ! (Note: This - and much else MPI stuff - must be modified for single precision to work!)
    call ls_mpi_buffer(MyFragment%energies,ndecenergies,master)
    call ls_mpi_buffer(MyFragment%pairdist,master)
    call ls_mpi_buffer(MyFragment%Eocc_err,master)
    call ls_mpi_buffer(MyFragment%Evir_err,master)
    call ls_mpi_buffer(MyFragment%Elag_err,master)
    CALL ls_mpi_buffer(MyFragment%flops_slaves,master)
    CALL ls_mpi_buffer(MyFragment%gpu_flops_slaves,master)
    call ls_mpi_buffer(MyFragment%slavetime_work,ndecmodels,master)
    call ls_mpi_buffer(MyFragment%slavetime_comm,ndecmodels,master)
    call ls_mpi_buffer(MyFragment%slavetime_idle,ndecmodels,master)
    call ls_mpi_buffer(MyFragment%RejectThr,2,master)
    call ls_mpi_buffer(MyFragment%DmaxAE,master)
    call ls_mpi_buffer(MyFragment%DmaxAOS,master)
    call ls_mpi_buffer(MyFragment%DaveAE,master)
    call ls_mpi_buffer(MyFragment%DaveAOS,master)
    call ls_mpi_buffer(MyFragment%DsdvAE,master)
    call ls_mpi_buffer(MyFragment%DsdvAOS,master)
    call ls_mpi_buffer(MyFragment%RmaxAE,master)
    call ls_mpi_buffer(MyFragment%RmaxAOS,master)
    call ls_mpi_buffer(MyFragment%RaveAE,master)
    call ls_mpi_buffer(MyFragment%RaveAOS,master)
    call ls_mpi_buffer(MyFragment%RsdvAE,master)
    call ls_mpi_buffer(MyFragment%RsdvAOS,master)

    ! Integer pointers
    ! ----------------
    ! Nullify and allocate stuff for receiver (global addtobuffer is false)
    if(.not. AddToBuffer) then
       nullify(MyFragment%occEOSidx)
       call mem_alloc(MyFragment%occEOSidx,MyFragment%noccEOS)
       nullify(MyFragment%virtEOSidx)
       call mem_alloc(MyFragment%virtEOSidx,MyFragment%nvirtEOS)
       nullify(MyFragment%occAOSidx)
       call mem_alloc(MyFragment%occAOSidx,MyFragment%noccLOC)
       nullify(MyFragment%virtAOSidx)
       call mem_alloc(MyFragment%virtAOSidx,MyFragment%nvirtLOC)
       nullify(MyFragment%coreidx)
       if(MyFragment%ncore>0) then
          call mem_alloc(MyFragment%coreidx,MyFragment%ncore)
       end if
       nullify(MyFragment%idxo)
       call mem_alloc(MyFragment%idxo,MyFragment%noccEOS)
       nullify(MyFragment%idxu)
       call mem_alloc(MyFragment%idxu,MyFragment%nvirtEOS)
       nullify(MyFragment%EOSatoms)
       call mem_alloc(MyFragment%EOSatoms,MyFragment%nEOSatoms)
       if(MyFragment%t1_stored) then ! only used for CC singles effects
          nullify(MyFragment%t1_occidx)
          call mem_alloc(MyFragment%t1_occidx,MyFragment%t1dims(2) )
          nullify(MyFragment%t1_virtidx)
          call mem_alloc(MyFragment%t1_virtidx,MyFragment%t1dims(1) )
       end if
       nullify(MyFragment%atoms_idx)
       call mem_alloc(MyFragment%atoms_idx,MyFragment%natoms)
       nullify(MyFragment%basis_idx)
       call mem_alloc(MyFragment%basis_idx,MyFragment%nbasis)
       nullify(MyFragment%cabsbasis_idx)
       IF(decinfo%F12)THEN
          call mem_alloc(MyFragment%cabsbasis_idx,MyFragment%nCabsAOOnly)
       ENDIF
    end if


    ! Buffer handling
    call ls_mpi_buffer(MyFragment%occEOSidx,MyFragment%noccEOS,master)
    call ls_mpi_buffer(MyFragment%virtEOSidx,MyFragment%nvirtEOS,master)
    call ls_mpi_buffer(MyFragment%occAOSidx,MyFragment%noccLOC,master)
    call ls_mpi_buffer(MyFragment%virtAOSidx,MyFragment%nvirtLOC,master)
    if(MyFragment%ncore>0) then
       call ls_mpi_buffer(MyFragment%coreidx,MyFragment%ncore,master)
    end if
    call ls_mpi_buffer(MyFragment%idxo,MyFragment%noccEOS,master)
    call ls_mpi_buffer(MyFragment%idxu,MyFragment%nvirtEOS,master)
    call ls_mpi_buffer(MyFragment%EOSatoms,MyFragment%nEOSatoms,master)
    call ls_mpi_buffer(MyFragment%atoms_idx,MyFragment%natoms,master)
    call ls_mpi_buffer(MyFragment%basis_idx,MyFragment%nbasis,master)
    IF(decinfo%F12)THEN
       call ls_mpi_buffer(MyFragment%cabsbasis_idx,MyFragment%nCabsAOOnly,master)
    ENDIF
    if(MyFragment%t1_stored) then ! only used for CC singles effects
       call ls_mpi_buffer(MyFragment%t1_occidx,MyFragment%t1dims(2),master)
       call ls_mpi_buffer(MyFragment%t1_virtidx,MyFragment%t1dims(1),master)
    end if


    ! Real pointers
    ! -------------
    ! (Note: This and much else MPI stuff must be modified for single precision to work!)
    if(.not. AddToBuffer) then
       nullify(MyFragment%OccContribs)
       call mem_alloc(MyFragment%OccContribs,MyFragment%noccLOC)
       nullify(MyFragment%VirtContribs)
       call mem_alloc(MyFragment%VirtContribs,MyFragment%nvirtLOC)

       if(MyFragment%CDset) then
          nullify(MyFragment%OccMat)
          call mem_alloc(MyFragment%OccMat,MyFragment%noccLOC,MyFragment%noccLOC)
          nullify(MyFragment%VirtMat)
          call mem_alloc(MyFragment%VirtMat,MyFragment%nvirtLOC,MyFragment%nvirtLOC)
       end if

       if(MyFragment%FAset) then
          nullify(MyFragment%CoFA)
          call mem_alloc(MyFragment%CoFA,MyFragment%nbasis,MyFragment%noccFA)
          nullify(MyFragment%CvFA)
          call mem_alloc(MyFragment%CvFA,MyFragment%nbasis,MyFragment%nvirtFA)
          !nullify(MyFragment%ppfockFA)
          !call mem_alloc(MyFragment%ppfockFA,MyFragment%noccFA,MyFragment%noccFA)
          !nullify(MyFragment%qqfockFA)
          !call mem_alloc(MyFragment%qqfockFA,MyFragment%nvirtFA,MyFragment%nvirtFA)
          if(.not. MyFragment%pairfrag) then
             nullify(MyFragment%CDocceival)
             call mem_alloc(MyFragment%CDocceival,MyFragment%noccFA)
             nullify(MyFragment%CDvirteival)
             call mem_alloc(MyFragment%CDvirteival,MyFragment%nvirtFA)
          end if
       end if

       if (MyFragment%t1_stored) then ! only used for CC singles effects
          nullify(MyFragment%t1)
          call mem_alloc(MyFragment%t1,MyFragment%t1dims(1),MyFragment%t1dims(2) )
       end if
    end if

    call ls_mpi_buffer(MyFragment%OccContribs,MyFragment%noccLOC,master)
    call ls_mpi_buffer(MyFragment%VirtContribs,MyFragment%nvirtLOC,master)
    if(MyFragment%CDset) then
       call ls_mpi_buffer(MyFragment%OccMat,MyFragment%noccLOC,MyFragment%noccLOC,master)
       call ls_mpi_buffer(MyFragment%VirtMat,MyFragment%nvirtLOC,MyFragment%nvirtLOC,master)
    end if
    if(MyFragment%FAset) then
       call ls_mpi_buffer(MyFragment%CoFA,MyFragment%nbasis,MyFragment%noccFA,master)
       call ls_mpi_buffer(MyFragment%CvFA,MyFragment%nbasis,MyFragment%nvirtFA,master)
       !call ls_mpi_buffer(MyFragment%ppfockFA,MyFragment%noccFA,MyFragment%noccFA,master)
       !call ls_mpi_buffer(MyFragment%qqfockFA,MyFragment%nvirtFA,MyFragment%nvirtFA,master)
       if(.not. MyFragment%pairfrag) then
          call ls_mpi_buffer(MyFragment%CDocceival,MyFragment%noccFA,master)
          call ls_mpi_buffer(MyFragment%CDvirteival,MyFragment%nvirtFA,master)
       end if
    end if
    if (MyFragment%t1_stored) then ! only used for CC singles effects
       call ls_mpi_buffer(MyFragment%t1,MyFragment%t1dims(1),MyFragment%t1dims(2),master)
    end if


    ! Reduced fragment types
    ! ----------------------
    if(.not. AddToBuffer) then
       call mem_alloc(MyFragment%REDfrags,DECinfo%nFRAGSred)
    end if
    do i=1,DECinfo%nFRAGSred
       call mpicopy_fragmentAOStype(MyFragment%REDfrags(i),comm)
    end do


    ! CCORBITAL types
    ! ---------------
    ! Allocate decorbitals
    if(.not. AddToBuffer) then
       nullify(MyFragment%occAOSorb)
       call mem_alloc(MyFragment%occAOSorb,MyFragment%noccLOC)
       nullify(MyFragment%virtAOSorb)
       call mem_alloc(MyFragment%virtAOSorb,MyFragment%nvirtLOC)
    end if

    ! Integers and reals inside decorbital sub-type
    do i=1,MyFragment%noccLOC ! occ orbitals
       call ls_mpi_buffer(MyFragment%occAOSorb(i)%orbitalnumber,master)
       call ls_mpi_buffer(MyFragment%occAOSorb(i)%centralatom,master)
       call ls_mpi_buffer(MyFragment%occAOSorb(i)%secondaryatom,master)
       call ls_mpi_buffer(MyFragment%occAOSorb(i)%numberofaos,master)

       ! Pointers inside decorbital sub-type (which again is inside decfrag type)
       if(.not. AddToBuffer) then
          Nullify(MyFragment%occAOSorb(i)%aos)
          call mem_alloc(MyFragment%occAOSorb(i)%aos,&
               & MyFragment%occAOSorb(i)%numberofaos)
       end if
       ! Buffer handling
       call ls_mpi_buffer(MyFragment%occAOSorb(i)%aos, &
            & MyFragment%occAOSorb(i)%numberofaos,master)
    end do

    do i=1,MyFragment%nvirtLOC ! virt orbitals
       call ls_mpi_buffer(MyFragment%virtAOSorb(i)%orbitalnumber,master)
       call ls_mpi_buffer(MyFragment%virtAOSorb(i)%centralatom,master)
       call ls_mpi_buffer(MyFragment%virtAOSorb(i)%secondaryatom,master)
       call ls_mpi_buffer(MyFragment%virtAOSorb(i)%numberofaos,master)
       if(.not. AddToBuffer) then ! allocate for slave
          Nullify(MyFragment%virtAOSorb(i)%aos)
          call mem_alloc(MyFragment%virtAOSorb(i)%aos,&
               & MyFragment%virtAOSorb(i)%numberofaos)
       end if
       call ls_mpi_buffer(MyFragment%virtAOSorb(i)%aos, &
            & MyFragment%virtAOSorb(i)%numberofaos,master)
    end do



    ! ===========================================================================
    !                              EXPENSIVE BOX
    ! ===========================================================================
    ! This is the information inside the "expensive" box in the decfrag type definition.
    ! This stuff is only copied if MyFragment%BasisInfoIsSet is true.
    ExpensiveBox: if(DoBasis) then

       ! Real pointers
       if(.not. AddToBuffer) then
          call mem_alloc(MyFragment%S,MyFragment%nbasis,MyFragment%nbasis)
          call mem_alloc(MyFragment%CoLOC,MyFragment%nbasis,MyFragment%noccLOC)
          call mem_alloc(MyFragment%CvLOC,MyFragment%nbasis,MyFragment%nvirtLOC)
          call mem_alloc(MyFragment%coreMO,MyFragment%nbasis,MyFragment%ncore)
          if(.not. DECinfo%noaofock) then
             call mem_alloc(MyFragment%fock,MyFragment%nbasis,MyFragment%nbasis)
          end if
          call mem_alloc(MyFragment%ccfock,MyFragment%ncore,MyFragment%ncore)
          call mem_alloc(MyFragment%ppfockLOC,MyFragment%noccLOC,MyFragment%noccLOC)
          call mem_alloc(MyFragment%qqfockLOC,MyFragment%nvirtLOC,MyFragment%nvirtLOC)
          if(MyFragment%FAset) then
            call mem_alloc(MyFragment%ppfockFA,MyFragment%noccFA,MyFragment%noccFA)
            call mem_alloc(MyFragment%qqfockFA,MyFragment%nvirtFA,MyFragment%nvirtFA)
          endif
       end if
       call ls_mpi_buffer(MyFragment%S,MyFragment%nbasis,MyFragment%nbasis,master)
       call ls_mpi_buffer(MyFragment%CoLOC,MyFragment%nbasis,MyFragment%noccLOC,master)
       call ls_mpi_buffer(MyFragment%CvLOC,MyFragment%nbasis,MyFragment%nvirtLOC,master)
       call ls_mpi_buffer(MyFragment%coreMO,MyFragment%nbasis,MyFragment%ncore,master)
       if(.not. DECinfo%noaofock) then
          call ls_mpi_buffer(MyFragment%fock,MyFragment%nbasis,MyFragment%nbasis,master)
       end if
       call ls_mpi_buffer(MyFragment%ccfock,MyFragment%ncore,MyFragment%ncore,master)
       call ls_mpi_buffer(MyFragment%ppfockLOC,MyFragment%noccLOC,MyFragment%noccLOC,master)
       call ls_mpi_buffer(MyFragment%qqfockLOC,MyFragment%nvirtLOC,MyFragment%nvirtLOC,master)
       if(MyFragment%FAset) then
         call ls_mpi_buffer(MyFragment%ppfockFA,MyFragment%noccFA,MyFragment%noccFA,master)
         call ls_mpi_buffer(MyFragment%qqfockFA,MyFragment%nvirtFA,MyFragment%nvirtFA,master)
       endif

       if(MyFragment%PNOset) then

         if( .not. AddToBuffer )then

           nullify(MyFragment%CLocPNO)
           call mem_alloc( MyFragment%CLocPNO, MyFragment%nspaces )

         endif

         do i = 1, MyFragment%nspaces
           call buffercopy_PNOSpaceInfo_struct(MyFragment%CLocPNO(i),master)
         enddo
       endif

       ! INTEGRAL LSITEM
       ! '''''''''''''''
       call mpicopy_lsitem(MyFragment%MyLsitem,comm)

       ! Basis info has now been set
       MyFragment%BasisInfoIsSet=.true.

    else

       ! Basis info was not set
       MyFragment%BasisInfoIsSet=.false.

    end if ExpensiveBox


    if(.not. AddToBuffer) then
       
       ! Point to FO data
       if(MyFragment%fragmentadapted.and.DoBasis) then
          call fragment_basis_point_to_FOs(MyFragment)
       else
          ! Point to local orbital data
          call fragment_basis_point_to_LOs(MyFragment)
       end if
    end if

  End subroutine mpicopy_fragment


  !> \brief MPI copy fragment AOS type
  !> \author Kasper Kristensen
  !> \date November 2014
  subroutine mpicopy_fragmentAOStype(MyFragmentAOS,comm)

    implicit none

    !> FragmentAOS information
    type(fragmentAOS),intent(inout) :: MyFragmentAOS
    !> Communicator
    integer(kind=ls_mpik),intent(in) :: comm
    integer(kind=ls_mpik) :: master
    master = 0

    CALL ls_mpi_buffer(MyFragmentAOS%noccAOS,master)
    CALL ls_mpi_buffer(MyFragmentAOS%nvirtAOS,master)
    CALL ls_mpi_buffer(MyFragmentAOS%FOT,master)

    ! Nullify and allocate stuff for receiver (global addtobuffer is false)
    if(.not. AddToBuffer) then
       call mem_alloc(MyFragmentAOS%occAOSidx,MyFragmentAOS%noccAOS)
       call mem_alloc(MyFragmentAOS%virtAOSidx,MyFragmentAOS%nvirtAOS)
    end if
    call ls_mpi_buffer(MyFragmentAOS%occAOSidx,MyFragmentAOS%noccAOS,master)
    call ls_mpi_buffer(MyFragmentAOS%virtAOSidx,MyFragmentAOS%nvirtAOS,master)

  end subroutine mpicopy_fragmentAOStype


  subroutine buffercopy_PNOSpaceInfo_struct(inf,master)
    implicit none
    type(PNOSpaceInfo),intent(inout) :: inf
    integer(kind=ls_mpik),intent(in) :: master
    call ls_mpi_buffer(inf%n,master)
    call ls_mpi_buffer(inf%ns1,master)
    call ls_mpi_buffer(inf%ns2,master)
    call ls_mpi_buffer(inf%rpd,master)
    call ls_mpi_buffer(inf%red1,master)
    call ls_mpi_buffer(inf%red2,master)
    call ls_mpi_buffer(inf%allocd,master)
    call ls_mpi_buffer(inf%s_associated,master)

    if(inf%allocd)then

      !allocate the arrays correctly
      if( .not. AddToBuffer )then

        call mem_alloc(inf%iaos,inf%n)

        if (inf%s_associated) then
          call mem_alloc(inf%s1,inf%ns1,inf%red1)
          call mem_alloc(inf%s2,inf%red2,inf%ns2)
          call mem_alloc(inf%d,inf%red1,inf%red2)
        else
          call mem_alloc(inf%d,inf%ns1,inf%ns2)
        endif
      endif
    
      call ls_mpi_buffer(inf%iaos,inf%n,master)
      if (inf%s_associated) then
        call ls_mpi_buffer(inf%s1,inf%ns1,inf%red1,master)
        call ls_mpi_buffer(inf%s2,inf%red2,inf%ns2,master)
        call ls_mpi_buffer(inf%d,inf%red1,inf%red2,master)
      else
        call ls_mpi_buffer(inf%d,inf%ns1,inf%ns2,master)
      endif

    endif
   
  end subroutine buffercopy_PNOSpaceInfo_struct

  subroutine share_E2_with_slaves(ccmodel,ppf,qqf,t2,xo,yv,Gbi,Had,no,nv,nb,omega2,s,lo)
    implicit none
    integer,intent(inout) :: ccmodel
    real(realk),pointer :: xo(:),yv(:),Gbi(:),Had(:)
    real(realk), intent(inout) :: ppf(:),qqf(:)
    integer :: no,nv,nb,s
    type(tensor),intent(inout) :: t2,omega2
    logical :: lo
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
    call ls_mpi_buffer(lo,infpar%master)
    call ls_mpi_buffer(ccmodel,infpar%master)
    if(master)oaddr=omega2%addr_p_arr
    call ls_mpi_buffer(oaddr,infpar%lg_nodtot,infpar%master)
    if(master)taddr=t2%addr_p_arr
    call ls_mpi_buffer(taddr,infpar%lg_nodtot,infpar%master)
    if(.not.master)then
       if(mem_is_background_buf_init())then
          call mem_pseudo_alloc(xo,i8*nb*no)
          call mem_pseudo_alloc(yv,i8*nb*nv)
          call mem_pseudo_alloc(Gbi,i8*nb*no)
          call mem_pseudo_alloc(Had,i8*nv*nb)
       else
          call mem_alloc(xo,nb*no)
          call mem_alloc(yv,nb*nv)
          call mem_alloc(Gbi,nb*no)
          call mem_alloc(Had,nv*nb)
       endif
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
      t2     = get_tensor_from_parr(taddr(infpar%lg_mynum+1))
      omega2 = get_tensor_from_parr(oaddr(infpar%lg_mynum+1))
    endif
  end subroutine share_E2_with_slaves



  !> \brief MPI communcation where CCSD and CC2 data is transferred
  !> \author Patrick Ettenhuber
  !> \date March 2012
  subroutine mpi_communicate_ccsd_calcdata(ccmodel,om2,t2,govov,xo,xv,yo,yv,MyLsItem,nbas,nvirt,nocc,iter,loc,res_nr)
     implicit none
     integer,intent(inout) :: ccmodel,res_nr
     type(mp2_batch_construction) :: bat
     integer            :: nbas,nocc,nvirt,ierr,iter
     !real(realk)        :: t2(:),govov(:)
     type(tensor),intent(inout) :: t2,govov,om2
     real(realk)        :: xo(:),xv(:),yo(:),yv(:)
     type(lsitem)       :: MyLsItem
     real(realk)        :: norm
     integer(kind=long) :: nelms
     integer :: i,n4,k
     integer :: gaddr(infpar%lg_nodtot)
     integer :: taddr(infpar%lg_nodtot)
     integer :: oaddr(infpar%lg_nodtot)
     logical :: loc
     character(tensor_MSG_LEN) :: msg
     logical :: master
     master=(infpar%lg_mynum==infpar%master)

     !communicate mylsitem and integers
     call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
     !    call ls_mpi_buffer(DECinfo%ccModel,infpar%master)
     call ls_mpi_buffer(nbas,infpar%master)
     call ls_mpi_buffer(nocc,infpar%master)
     call ls_mpi_buffer(nvirt,infpar%master)
     call ls_mpi_buffer(iter,infpar%master)
     call ls_mpi_buffer(loc,infpar%master)
     call ls_mpi_buffer(ccmodel,infpar%master)
     call ls_mpi_buffer(res_nr,infpar%master)
     if(.not.loc)then
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

        nelms = nbas*nocc
        call ls_mpibcast(xo,nelms,infpar%master,infpar%lg_comm)
        call ls_mpibcast(yo,nelms,infpar%master,infpar%lg_comm)

        nelms = nbas*nvirt
        call ls_mpibcast(xv,nelms,infpar%master,infpar%lg_comm)
        call ls_mpibcast(yv,nelms,infpar%master,infpar%lg_comm)

     else
        if(.not.loc)then
           govov = get_tensor_from_parr(gaddr(infpar%lg_mynum+1))
           t2    = get_tensor_from_parr(taddr(infpar%lg_mynum+1))
           om2   = get_tensor_from_parr(oaddr(infpar%lg_mynum+1))
        endif
     endif
  end subroutine mpi_communicate_ccsd_calcdata

  !> \brief mpi communcation where ccsd(t) data is transferred
  !> \author Janus Juul Eriksen
  !> \date February 2013
  subroutine mpi_communicate_ccsdpt_calcdata(nocc,nvirt,nbasis,vovo,ccsd_t2,mylsitem,print_frags,abc,ccsd_t1)

    implicit none

    integer                             :: nocc,nvirt,nbasis,ierr
    type(tensor), intent(inout)         :: vovo,ccsd_t2
    type(lsitem)                        :: mylsitem
    logical                             :: print_frags,abc
    integer,dimension(infpar%lg_nodtot) :: vovo_addr,t2_addr
    type(tensor), intent(inout), optional :: ccsd_t1

    ! communicate mylsitem and integers
    call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
    call ls_mpi_buffer(DECinfo%memory,infpar%master)
    call ls_mpi_buffer(nbasis,infpar%master)
    call ls_mpi_buffer(nocc,infpar%master)
    call ls_mpi_buffer(nvirt,infpar%master)
    call ls_mpi_buffer(print_frags,infpar%master)
    call ls_mpi_buffer(abc,infpar%master)
    if (infpar%lg_mynum .eq. infpar%master) t2_addr = ccsd_t2%addr_p_arr
    call ls_mpi_buffer(t2_addr,infpar%lg_nodtot,infpar%master)
    if (infpar%lg_mynum .eq. infpar%master) vovo_addr = vovo%addr_p_arr
    call ls_mpi_buffer(vovo_addr,infpar%lg_nodtot,infpar%master)
    call mpicopy_lsitem(mylsitem,infpar%lg_comm)
    call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)

    if (infpar%lg_mynum .ne. infpar%master) then

       ccsd_t2 = get_tensor_from_parr(t2_addr(infpar%lg_mynum+1))
       vovo = get_tensor_from_parr(vovo_addr(infpar%lg_mynum+1))

    endif

    ccsd_t2%access_type = AT_ALL_ACCESS
    vovo%access_type = AT_ALL_ACCESS

    if ((infpar%lg_mynum .eq. infpar%master) .and. (.not. print_frags) .and. present(ccsd_t1)) then

       call ls_mpibcast(ccsd_t1%elm1,nvirt*nocc,infpar%master,infpar%lg_comm)

    endif

  end subroutine mpi_communicate_ccsdpt_calcdata

#ifdef MOD_UNRELEASED
  !> Purpose: Get job list to have a good load balance in the
  !           main loop of the MO-CCSD residual calculations
  !
  !> Author:  Pablo Baudin
  !> Date:    January 2014
  subroutine get_mo_ccsd_joblist(MOinfo, joblist, pgmo_diag, pgmo_up)

    implicit none

    type(MObatchInfo), intent(in) :: MOinfo
    type(tensor), intent(in) :: pgmo_diag, pgmo_up
    integer, intent(inout) :: joblist(:)

    integer, pointer :: workloads(:), easytrace(:,:), work_in_node(:)
    integer :: swapar(3)
    integer :: nnod, njob, swap, i, j, next_nod, ijob, dpos, didx, dwidx

    nnod = infpar%lg_nodtot
    njob = MOinfo%Nbatch

    call mem_alloc(workloads,njob)
    call mem_alloc(easytrace,njob,3)
    call mem_alloc(work_in_node,nnod)
 
    workloads=MOinfo%dimTot
 
    ! Associate a integer three-vector to each job:
    ! the first number is the tile index
    ! the second number is 0 for pgmo_diag array
    ! and 1 for pgmo_up array.
    ! and the third number is the global job index
    do i = 1, njob
      easytrace(i,1:2) = MOinfo%tileInd(i,:)
      easytrace(i,3) = i
    end do

    ! Sort the jobs according to their size, and keep track of indices
    do i=1,njob
      do j=i+1,njob
        if( workloads(j) > workloads(i) )then

          swap=workloads(j)
          workloads(j)=workloads(i)
          workloads(i)=swap

          swapar=easytrace(j,:)
          easytrace(j,:)=easytrace(i,:)
          easytrace(i,:)=swapar
        endif
      enddo
    enddo

    ! Associate a rank node to each job:
    work_in_node = 0
    ! for the nnod first jobs, each node treat the batch that stand on its memory.
    do ijob=1,min(nnod,njob)
      if (easytrace(ijob,2)==0) then 
        ! tile in pgmo_diag array
        call get_residence_of_tile(pgmo_diag,easytrace(ijob,1),next_nod,dpos,didx,dwidx)
      else 
        ! tile in pgmo_up array
        call get_residence_of_tile(pgmo_up,easytrace(ijob,1),next_nod,dpos,didx,dwidx)
      end if
      next_nod = next_nod + 1
      ! Update joblist and workload
      joblist(ijob) = next_nod 
      work_in_node(next_nod) = work_in_node(next_nod) + workloads(ijob)
    end do

    ! If more jobs then attribute node depending on workload:
    if (nnod<njob) then
      do ijob=nnod+1,njob

        ! get node with smallest work:
        next_nod = 1
        do i=2, nnod
          if (work_in_node(i)<work_in_node(next_nod)) next_nod = i
        end do

        ! Update joblist and workload
        joblist(ijob) = next_nod
        work_in_node(next_nod) = work_in_node(next_nod) + workloads(ijob)
      end do
    end if
 
    ! go back to initial order:
    do i=1,njob
      do j=i+1,njob
        if( easytrace(j,3) < easytrace(i,3) )then
          swap=joblist(j)
          joblist(j)=joblist(i)
          joblist(i)=swap

          swapar=easytrace(j,:)
          easytrace(j,:)=easytrace(i,:)
          easytrace(i,:)=swapar
        endif
      enddo
    enddo

    call mem_dealloc(workloads)
    call mem_dealloc(easytrace)
    call mem_dealloc(work_in_node)

  end subroutine get_mo_ccsd_joblist
#endif

  !> \brief get a suitable job distribution in mpi calculations
  !> \author Patrick Ettenhuber
  !> \date March 2012
  subroutine distribute_mpi_jobs(mpi_task_distribution,nbA,nbG,batchdimAlpha,&
             & batchdimGamma,myload,lg_nnod,lg_me,ccsdscheme,no,nv,nb,b2oa,b2og)
    implicit none

    integer, intent(inout) :: myload
    integer,intent(in):: nbA,nbG,batchdimAlpha(nbA),batchdimGamma(nbG)
    integer(kind=ls_mpik),intent(in)::lg_nnod,lg_me
    integer,intent(out)::mpi_task_distribution(nbA*nbG)
    integer,intent(in),optional :: ccsdscheme,no,nv,nb
    type(batchtoorb),intent(in),optional :: b2oa(nbA),b2og(nbG)
    integer(kind=8),pointer :: workloads(:),jobsize_per_node(:)
    integer(kind=8) :: myload64,swap
    integer(kind=8),pointer :: node_rank_by_job(:),easytrace(:)
    logical,pointer :: touched(:)
    logical :: all_touched
    type(traceback), pointer::trace(:)
    integer :: i,j,k,l,counter,ialpha,igamma,actual_node
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
    call  mem_alloc(jobsize_per_node,lg_nnod)
    call  mem_alloc(node_rank_by_job,lg_nnod)
    call  mem_alloc(trace,nbA*nbG)
    !some needed initializations for the matrices
    touched = .false.
    jobsize_per_node = 0
    do k=1, lg_nnod
      node_rank_by_job(k) = k-1
    enddo

    if (.false.) then
     ! simple modulo distribtution --> tured off for the moment
     ! fill the jobsize_per_node matrix with the values from mpi_task_distribution
      do i=1,nbG
        do j=1,nbA
          actual_node=MODULO(j-1+(i-1)*MODULO(nbA,lg_nnod),lg_nnod)
          mpi_task_distribution((j-1)*nbG+i)=actual_node
          jobsize_per_node(actual_node+1)=jobsize_per_node(actual_node+1) & 
          &+ int(i8*batchdimGamma(i)*batchdimAlpha(j),kind=8)
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
         workloads((j-1)*nbG+i)   = int((i8*la)*lg,kind=8)

         !in CCSD the workloads are not equally distributed, and a more accurate
         !estimation of the work has to be done
         if(present(ccsdscheme))then
           fa=b2oA(j)%orbindex(1)
           fg=b2oG(i)%orbindex(1)
           !account for integral work more "correctly"
           workloads((j-1)*nbG+i) = workloads((j-1)*nbG+i) * nb * nb * 2
           !account for Kobayashi-Term
           if(fa<=fg+lg-1 )then
             workloads((j-1)*nbG+i) = workloads((j-1)*nbG+i) + int(((i8*nv**2)*no**2*la*lg)/4,kind=8)
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
      do i=1,nbG*nbA-MODULO(nbG*nbA,lg_nnod),lg_nnod
        !distribute the work according to previous load
        do k=1, lg_nnod
          jobsize_per_node(k)=jobsize_per_node(k) + workloads(i+node_rank_by_job(k))
          ialpha = trace(easytrace(i+node_rank_by_job(k))+1)%na
          igamma = trace(easytrace(i+node_rank_by_job(k))+1)%ng
          mpi_task_distribution((ialpha-1)*nbg+igamma) = k-1

          ! keep track of which elements were done if an error occurs, quit, since
          ! the result would
          if(.not. touched((ialpha-1)*nbg+igamma)) then
            touched((ialpha-1)*nbg+igamma) = .true.
          else
            if (lg_me==0)then
              call LSQUIT('Element had already been assigned in distribute_mpi_jobs',-1)
            endif
          endif
        enddo

        ! get the new node ranks by counting how many of the nodes have smaller jobs
        ! --> is an initial rank, because it is ordered the other way round as the jobs
        do k=1, lg_nnod
          counter = 0
          do l=1,lg_nnod
            if (jobsize_per_node(l) < jobsize_per_node(k)) counter=counter+1
          enddo
          node_rank_by_job(k) = counter
        enddo
        ! if some have the same amount of work, just assign jobs uniquely
        do k=1, lg_nnod
          do l=k+1, lg_nnod
            if (node_rank_by_job(k)==node_rank_by_job(l))then
              node_rank_by_job(l) = node_rank_by_job(l) +1
            endif
          enddo
        enddo
      enddo

      !distribute the rest of the jobs according to load
      i = nbG*nbA-MODULO(nbG*nbA,lg_nnod)
      do k=1, lg_nnod
        actual_node= node_rank_by_job(k) + 1
        if(MODULO(nbG*nbA,lg_nnod) >= actual_node)then
          jobsize_per_node(k)=jobsize_per_node(k) + workloads(i+actual_node)
          ialpha = trace(easytrace(i+actual_node)+1)%na
          igamma = trace(easytrace(i+actual_node)+1)%ng
          mpi_task_distribution((ialpha-1)*nbg+igamma) = k-1
          if(.not. touched((ialpha-1)*nbg+igamma)) then
            touched((ialpha-1)*nbg+igamma) = .true.
          else
            if (lg_me==0)then
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
        if (mpi_task_distribution((j-1)*nbG+i)==lg_me)then
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

    if ((myload64 /= jobsize_per_node(lg_me+1)).or. .not. all_touched)then
      print *, "something is wrong ... quitting",myload64,jobsize_per_node(lg_me+1),all_touched
      write(DECinfo%output,*) 'workloads per node:'
      write(DECinfo%output,*) node_rank_by_job
      write(DECinfo%output,*) jobsize_per_node
      write(DECinfo%output,*) mpi_task_distribution
      print *, 'workloads per node:'
      print *,workloads
      print *, node_rank_by_job
      print *, jobsize_per_node
      print *, mpi_task_distribution
      print *, touched
      print *, lg_me,kind(lg_me),kind(jobsize_per_node(lg_me+1)),kind(myload64)
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
  subroutine bcast_dec_fragment_joblist(jobs,comm)

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


  end subroutine bcast_dec_fragment_joblist



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
       call mem_alloc(jobs%dofragopt,jobs%njobs)
       call mem_alloc(jobs%esti,jobs%njobs)
       call mem_alloc(jobs%nslaves,jobs%njobs)
       call mem_alloc(jobs%nocc,jobs%njobs)
       call mem_alloc(jobs%nvirt,jobs%njobs)
       call mem_alloc(jobs%nbasis,jobs%njobs)
       call mem_alloc(jobs%ntasks,jobs%njobs)
       call mem_alloc(jobs%flops,jobs%njobs)
       call mem_alloc(jobs%gpu_flops,jobs%njobs)
       call mem_alloc(jobs%LMtime,jobs%njobs)
       call mem_alloc(jobs%workt,jobs%njobs)
       call mem_alloc(jobs%commt,jobs%njobs)
       call mem_alloc(jobs%idlet,jobs%njobs)
       call mem_alloc(jobs%comm_gl_master_time,jobs%njobs)
    end if

    ! Buffer handling for pointers
    call ls_mpi_buffer(jobs%atom1,jobs%njobs,master)
    call ls_mpi_buffer(jobs%atom2,jobs%njobs,master)
    call ls_mpi_buffer(jobs%jobsize,jobs%njobs,master)
    call ls_mpi_buffer(jobs%jobsdone,jobs%njobs,master)
    call ls_mpi_buffer(jobs%dofragopt,jobs%njobs,master)
    call ls_mpi_buffer(jobs%esti,jobs%njobs,master)
    call ls_mpi_buffer(jobs%nslaves,jobs%njobs,master)
    call ls_mpi_buffer(jobs%nocc,jobs%njobs,master)
    call ls_mpi_buffer(jobs%nvirt,jobs%njobs,master)
    call ls_mpi_buffer(jobs%nbasis,jobs%njobs,master)
    call ls_mpi_buffer(jobs%ntasks,jobs%njobs,master)
    call ls_mpi_buffer(jobs%flops,jobs%njobs,master)
    call ls_mpi_buffer(jobs%gpu_flops,jobs%njobs,master)
    call ls_mpi_buffer(jobs%LMtime,jobs%njobs,master)
    call ls_mpi_buffer(jobs%workt,jobs%njobs,master)
    call ls_mpi_buffer(jobs%commt,jobs%njobs,master)
    call ls_mpi_buffer(jobs%idlet,jobs%njobs,master)
    call ls_mpi_buffer(jobs%comm_gl_master_time,jobs%njobs,master)

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
  !> Here some of the PDM stuff has to be redefined according to the new groups,
  !> because the PDM memory allocation has to happen in the smaller groups now
  !> \author Kasper Kristensen, modified by Patrick Ettenhuber
  !> \date May 2012
  subroutine dec_half_local_group(print_)
    implicit none
    logical, intent(in), optional :: print_
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
    call divide_local_mpi_group(ngroups,groupdims,print_=print_)
    call new_group_reset_persistent_array

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
    type(joblist),intent(inout) :: jobs
    !> Time measured by global master when fragments were calculated
    real(realk),intent(in) :: mastertime
    !> String to print describing job list
    character(*),intent(in) :: string
    real(realk) :: Gflops, avflop, minflop, maxflop,tmp,totflops,tottime_actual,localloss
    real(realk) :: gpuGflops,totgpuflops
    real(realk) :: globalloss, tottime_ideal, slavetime, localuse
    integer :: i, minidx, maxidx,N
    logical :: Fragoptjobs

    if(jobs%njobs<1) then
       write(DECinfo%output,*) 'MPI fragment statistics: No jobs to print!'
       return
    end if


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
    write(DECinfo%output,*) '#vir   : Number of virtual AOS orbitals in fragment'
    write(DECinfo%output,*) '#bas   : Number of basis functions in fragment'
    write(DECinfo%output,*) 'GrpSiz : Number of slaves used for fragment incl. local master (Groupe Size)'
    write(DECinfo%output,*) '#tasks : Number of integral tasks for fragment (nalpha*ngamma)'
    write(DECinfo%output,*) 'GFLOPS : Accumulated GFLOPS for fragment (local master AND local slaves)'
    write(DECinfo%output,*) 'GGFLOPS: Accumulated GPU GFLOPS for fragment (local master AND local slaves)'
    write(DECinfo%output,*) 'Time(s): Time (in seconds) used by local master (NOT local slaves)'
    write(DECinfo%output,*) 'CommT  : Time (in seconds) used by local master to communicate data to the'
    write(DECinfo%output,*) '         global master (correspond to job n+1)'
    write(DECinfo%output,*) 'Load   : Load distribution measure (ideally 1.0, smaller in practice)'
    write(DECinfo%output,*) '      Load1 = (Sum of working and communication times for ALL MPI processes in slot)'
    write(DECinfo%output,*) '           / (slotsize * [Local master time] )'
    write(DECinfo%output,*) '      Load2 = (Sum of working times for ALL MPI processes in slot)'
    write(DECinfo%output,*) '           / (slotsize * [Local master time] )'
    write(DECinfo%output,*)
    write(DECinfo%output,*) 'Note: The Load given below is set to -1 for cases where it is not implemented'
    write(DECinfo%output,*) '      Similarly, GFLOPS is set to -1 if you have not linked to the PAPI library'
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(4X,a,X,a,X,a,X,a,X,a,X,a,3X,a,4X,a,4X,a,6X,a,4X,a,2X,a)') 'Job', '#occ', &
         & '#vir', '#bas', 'GrpSiz', '#tasks', 'GFLOPS', 'GGFLOPS', 'Time(s)', 'CommT', 'Load1', 'Load2'

    avflop         = 0.0E0_realk
    totflops       = 0.0E0_realk
    totgpuflops    = 0.0E0_realk
    tottime_actual = 0.0E0_realk
    slavetime      = 0.0E0_realk

    minflop        = huge(minflop)
    maxflop        = tiny(maxflop)
    minidx         = 0
    maxidx         = 0
    N              = 0
    Fragoptjobs    = any(jobs%dofragopt)

    do i=1,jobs%njobs
       ! If nocc is zero, the job was not done and we do not print it
       if(jobs%nocc(i)==0) cycle
       N=N+1

       ! Giga flops for fragment
#ifdef VAR_PAPI
       Gflops = jobs%flops(i)*1.0e-9_realk
#else
       Gflops=-1.0_realk
#endif
       gpuGflops = jobs%gpu_flops(i)*1.0e-9_realk
       ! Update total number of flops
       totflops = totflops + Gflops
       totgpuflops = totgpuflops + gpuGflops
       ! Update total time used by ALL nodes (including dead time by local slaves)
       tottime_actual = tottime_actual + jobs%LMtime(i)*jobs%nslaves(i)
       ! Effective slave time (WITHOUT dead time by slaves)
       slavetime = slavetime + jobs%workt(i) + jobs%commt(i)

       if (Fragoptjobs .and. .not. jobs%dofragopt(i)) then
           write(DECinfo%output,&
              &'("------------------------------------- Fragment optimization done -------------------------------------")')
           Fragoptjobs = .false.
       endif

       write(DECinfo%output,'(6i8,3X,5g11.3,a)') i, jobs%nocc(i), jobs%nvirt(i), jobs%nbasis(i),&
            & jobs%nslaves(i), jobs%ntasks(i), Gflops, gpuGflops, jobs%LMtime(i), &
            &(jobs%workt(i)+jobs%commt(i))/(jobs%LMtime(i)*jobs%nslaves(i)), &
            &(jobs%workt(i))/(jobs%LMtime(i)*jobs%nslaves(i)), 'STAT'

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

       write(DECinfo%output,*) '-----------------------------------------------------------------------------'
       write(DECinfo%output,'(1X,a,i8)')    'Number of jobs done     = ', N
       write(DECinfo%output,'(1X,a,g12.3)') 'Time-to-solution (s)    = ', mastertime
       write(DECinfo%output,'(1X,a,g12.3)') 'TOTAL time(s) all nodes = ', tottime_actual
#ifdef VAR_PAPI
       ! Average flops per node per sec
       avflop = totflops/tottime_actual
       write(DECinfo%output,'(1X,a,g12.3)') 'TOTAL Gflops     = ', totflops
       write(DECinfo%output,'(1X,a,g12.3)') 'TOTAL GPU Gflops = ', totgpuflops
       write(DECinfo%output,'(1X,a,g12.3)') 'AVERAGE Gflops/s per MPI process= ', avflop
       write(DECinfo%output,'(1X,a,g12.3,a,i8)') 'MINIMUM Gflops/s per MPI process = ', &
            & minflop, ' for job ', minidx
       write(DECinfo%output,'(1X,a,g12.3,a,i8)') 'MAXIMUM Gflops/s per MPI process = ', &
            & maxflop, ' for job ', maxidx
#endif
    write(DECinfo%output,'(1X,a,g12.3)') 'Global MPI loss (%)     = ', globalloss
       !if(.not. any(jobs%dofragopt) ) then
          ! Only print local loss when it is actually implemented
          write(DECinfo%output,'(1X,a,g12.3)') 'Local MPI loss (%)      = ', localloss
          write(DECinfo%output,'(1X,a,g12.3)') 'Total MPI loss (%)      = ', localloss+globalloss
       !end if
    write(DECinfo%output,*) '-----------------------------------------------------------------------------'

    end if
    write(DECinfo%output,*)
    write(DECinfo%output,*)

  end subroutine print_MPI_fragment_statistics


  !> \brief Bcast DEC setting structure
  !> \author Kasper Kristensen
  !> \date June 2013
  subroutine mpibcast_dec_settings(DECitem,comm)
    implicit none
    integer(kind=ls_mpik) :: comm
    type(decsettings) :: DECitem

    call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,comm)
    call mpicopy_dec_settings(DECitem)
    call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,comm)

  end subroutine mpibcast_dec_settings

#ifdef MOD_UNRELEASED
  !> Purpose: Communicate data to the slaves needed to get MO integral.
  !           get_packed_gmo routine.
  !
  !> Author:  Pablo Baudin
  !> Date:    December 2013
  subroutine mpi_communicate_get_gmo_data(mo_ccsd,MyLsItem,Co,Cv, &
             & pgmo_diag,pgmo_up,nbas,nocc,nvir,nbatch)

    implicit none
     
    !> number of orbitals:
    integer :: nbas, nocc, nvir
    !> Number of MO batches
    integer :: nbatch
    !> SCF transformation matrices:
    real(realk), pointer  :: Co(:,:), Cv(:,:)
    !> performed MO-based CCSD calculation ?
    logical :: mo_ccsd
    !> array with gmo on output:
    type(tensor) :: pgmo_diag, pgmo_up
    !> LS item information
    type(lsitem) :: MyLsItem

    integer :: pgmo_diag_addr(infpar%lg_nodtot)   
    integer :: pgmo_up_addr(infpar%lg_nodtot)   
    integer :: ntot
    logical :: master

    master = (infpar%lg_mynum == infpar%master)

    call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
    call ls_mpi_buffer(mo_ccsd,infpar%master)
    call ls_mpi_buffer(nbas,infpar%master)
    call ls_mpi_buffer(nocc,infpar%master)
    call ls_mpi_buffer(nvir,infpar%master)
    call ls_mpi_buffer(nbatch,infpar%master)
    if(.not.master)then
      call mem_alloc(Co,nbas,nocc)
      call mem_alloc(Cv,nbas,nvir)
    endif
    call ls_mpi_buffer(Co,nbas,nocc,infpar%master)
    call ls_mpi_buffer(Cv,nbas,nvir,infpar%master)

    if (master) pgmo_diag_addr=pgmo_diag%addr_p_arr
    call ls_mpi_buffer(pgmo_diag_addr,infpar%lg_nodtot,infpar%master)

    if (nbatch>1) then 
      if (master) pgmo_up_addr=pgmo_up%addr_p_arr
      call ls_mpi_buffer(pgmo_up_addr,infpar%lg_nodtot,infpar%master)
    end if

    call mpicopy_lsitem(MyLsItem,infpar%lg_comm)
    call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)

    if (.not.master) then
      pgmo_diag = get_tensor_from_parr(pgmo_diag_addr(infpar%lg_mynum+1))
      if (nbatch>1) pgmo_up   = get_tensor_from_parr(pgmo_up_addr(infpar%lg_mynum+1))
    endif

  end subroutine mpi_communicate_get_gmo_data


  !> Purpose: Communicate data to the slaves needed to calculate 
  !           MO-CCSD residual.
  !
  !> Author:  Pablo Baudin
  !> Date:    January 2014
  subroutine mpi_communicate_moccsd_data(ccmodel,pgmo_diag,pgmo_up,t1,t2,om2, &
             & govov,nbas,nocc,nvir,iter,MOinfo,MyLsItem,loc)

    implicit none
     
    !> CC model
    integer,intent(inout) :: ccmodel
    integer :: nbas, nocc, nvir, iter
    !> MO pack integrals; amplitudes and residuals:
    type(tensor) :: pgmo_diag, pgmo_up
    type(tensor) :: govov
    type(tensor) :: t1
    type(tensor) :: t2
    type(tensor) :: om2
     
    !> LS item with information needed for integrals
    type(lsitem) :: MyLsItem
     
    !> Batches info:
    type(MObatchInfo) :: MOinfo
    logical :: loc

    integer :: pgmo_diag_addr(infpar%lg_nodtot)   
    integer :: pgmo_up_addr(infpar%lg_nodtot)   
    integer :: t1addr(infpar%lg_nodtot), t2addr(infpar%lg_nodtot)
    integer :: gaddr(infpar%lg_nodtot), oaddr(infpar%lg_nodtot)
    integer :: ntot, k
    integer(kind=long) :: nelms
    logical :: master

    master = (infpar%lg_mynum == infpar%master)

    call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
    call ls_mpi_buffer(nbas,infpar%master)
    call ls_mpi_buffer(nocc,infpar%master)
    call ls_mpi_buffer(nvir,infpar%master)
    call ls_mpi_buffer(iter,infpar%master)
    call ls_mpi_buffer(MOinfo%nbatch,infpar%master)
    call ls_mpi_buffer(loc,infpar%master)
    call ls_mpi_buffer(ccmodel,infpar%master)
    if (.not.master) then
      call mem_alloc(MOinfo%dimInd1,MOinfo%nbatch)
      call mem_alloc(MOinfo%dimInd2,MOinfo%nbatch)
      call mem_alloc(MOinfo%StartInd1,MOinfo%nbatch)
      call mem_alloc(MOinfo%StartInd2,MOinfo%nbatch)
      call mem_alloc(MOinfo%dimTot,MOinfo%nbatch)
      call mem_alloc(MOinfo%tileInd,MOinfo%nbatch,2)
    end if
    call ls_mpi_buffer(MOinfo%dimInd1,MOinfo%nbatch,infpar%master)
    call ls_mpi_buffer(MOinfo%dimInd2,MOinfo%nbatch,infpar%master)
    call ls_mpi_buffer(MOinfo%StartInd1,MOinfo%nbatch,infpar%master)
    call ls_mpi_buffer(MOinfo%StartInd2,MOinfo%nbatch,infpar%master)
    call ls_mpi_buffer(MOinfo%dimTot,MOinfo%nbatch,infpar%master)
    call ls_mpi_buffer(MOinfo%tileInd,MOinfo%nbatch,2,infpar%master)

    if(.not.loc)then
      if(master)t1addr=t1%addr_p_arr
      call ls_mpi_buffer(t1addr,infpar%lg_nodtot,infpar%master)
      if(master)gaddr=govov%addr_p_arr
      call ls_mpi_buffer(gaddr,infpar%lg_nodtot,infpar%master)
      if(master)t2addr=t2%addr_p_arr
      call ls_mpi_buffer(t2addr,infpar%lg_nodtot,infpar%master)
      if(master)oaddr=om2%addr_p_arr
      call ls_mpi_buffer(oaddr,infpar%lg_nodtot,infpar%master)
    endif

    if (master) pgmo_diag_addr=pgmo_diag%addr_p_arr
    call ls_mpi_buffer(pgmo_diag_addr,infpar%lg_nodtot,infpar%master)

    if (MOinfo%nbatch>1) then 
      if (master) pgmo_up_addr=pgmo_up%addr_p_arr
      call ls_mpi_buffer(pgmo_up_addr,infpar%lg_nodtot,infpar%master)
    end if
 
    call mpicopy_lsitem(MyLsItem,infpar%lg_comm)
    call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)

    !communicate rest of the quantities, master here, slaves back in the slave
    !routine, due to crappy pointer/non-pointer issues (->allocations)
    if(master)then

      !split messages in 2GB parts, compare to counterpart in
      !ccsd_data_preparation

      nelms = int(i8*nvir*nvir*nocc*nocc,kind=8)
      call ls_mpibcast(t2%elm1,nelms,infpar%master,infpar%lg_comm)
      call ls_mpibcast(govov%elm1,nelms,infpar%master,infpar%lg_comm)
    else
      if(.not.loc)then
        t1    = get_tensor_from_parr(t1addr(infpar%lg_mynum+1))
        govov = get_tensor_from_parr(gaddr(infpar%lg_mynum+1))
        t2    = get_tensor_from_parr(t2addr(infpar%lg_mynum+1))
        om2   = get_tensor_from_parr(oaddr(infpar%lg_mynum+1))
      endif
      pgmo_diag = get_tensor_from_parr(pgmo_diag_addr(infpar%lg_mynum+1))
      if (MOinfo%nbatch>1) pgmo_up = get_tensor_from_parr(pgmo_up_addr(infpar%lg_mynum+1))
    endif

  end subroutine mpi_communicate_moccsd_data
#endif

  !> \brief Copy DEC setting structure to buffer (master)
  !> or read from buffer (slave)
  !> \author Kasper Kristensen
  !> \date June 2013
  subroutine mpicopy_dec_settings(DECitem)
    implicit none
    type(decsettings) :: DECitem
    integer(kind=ls_mpik) :: master
    master = 0

    call ls_mpi_buffer(DECitem%SNOOP,Master)
    call ls_mpi_buffer(DECitem%SNOOPjustHF,Master)
    call ls_mpi_buffer(DECitem%SNOOPmaxiter,Master)
    call ls_mpi_buffer(DECitem%SNOOPthr,Master)
    call ls_mpi_buffer(DECitem%SNOOPmaxdiis,Master)
    call ls_mpi_buffer(DECitem%SNOOPdebug,Master)
    call ls_mpi_buffer(DECitem%SNOOPsamespace,Master)
    call ls_mpi_buffer(DECitem%SNOOPlocalize,Master)
    call ls_mpi_buffer(DECitem%SNOOPrestart,Master)
    call ls_mpi_buffer(DECitem%SNOOPonesub,Master)
    call ls_mpi_buffer(DECitem%CCexci,Master)
    call ls_mpi_buffer(DECitem%JacobianNumEival,Master)
    call ls_mpi_buffer(DECitem%JacobianLHTR,Master)
    call ls_mpi_buffer(DECitem%JacobianThr,Master)
    call ls_mpi_buffer(DECitem%JacobianMaxSubspace,Master)
    call ls_mpi_buffer(DECitem%JacobianInitialSubspace,Master)
    call ls_mpi_buffer(DECitem%JacobianMaxIter,Master)
    call ls_mpi_buffer(DECitem%JacobianPrecond,Master)
    call ls_mpi_buffer(DECitem%SinglesEW1,Master)
    call ls_mpi_buffer(DECitem%LW1,Master)
    call ls_mpi_buffer(DECitem%P_EOM_MBPT2,Master)
    call ls_mpi_buffer(DECitem%doDEC,Master)
    call ls_mpi_buffer(DECitem%DECCO,Master)
    call ls_mpi_buffer(DECitem%DECNP,Master)
    call ls_mpi_buffer(DECitem%frozencore,Master)
    call ls_mpi_buffer(DECitem%full_molecular_cc,Master)
    call ls_mpi_buffer(DECitem%use_canonical,Master)
    call ls_mpi_buffer(DECitem%simulate_full,Master)
    call ls_mpi_buffer(DECitem%simulate_natoms,Master)
    call ls_mpi_buffer(DECitem%InclFullMolecule,Master)
    call ls_mpi_buffer(DECitem%use_system_memory_info,Master)
    call dec_set_model_names(DECitem)
    call ls_mpi_buffer(DECitem%ccModel,Master)
    call ls_mpi_buffer(DECitem%use_singles,Master)
    call ls_mpi_buffer(DECitem%gcbasis,Master)
    call ls_mpi_buffer(DECitem%HFrestart,Master)
    call ls_mpi_buffer(DECitem%DECrestart,Master)
    call ls_mpi_buffer(DECitem%EnforceRestart,Master)
    call ls_mpi_buffer(DECitem%TimeBackup,Master)
    call ls_mpi_buffer(DECitem%read_dec_orbitals,Master)
    call ls_mpi_buffer(DECitem%only_generate_DECorbs,Master)
    call ls_mpi_buffer(DECitem%memory,Master)
    call ls_mpi_buffer(DECitem%memory_defined,Master)
    call ls_mpi_buffer(DECitem%fullmolecule_memory,Master)    
    call ls_mpi_buffer(DECitem%bg_memory,Master)
    call ls_mpi_buffer(DECitem%array4OnFile,Master)
    call ls_mpi_buffer(DECitem%array4OnFile_specified,Master)
    call ls_mpi_buffer(DECitem%SinglesPolari,Master)
    call ls_mpi_buffer(DECitem%singlesthr,Master)
    call ls_mpi_buffer(DECitem%convert64to32,Master)
    call ls_mpi_buffer(DECitem%convert32to64,Master)
    call ls_mpi_buffer(DECitem%CCSDnosaferun,Master)
    call ls_mpi_buffer(DECitem%solver_par,Master)
    call ls_mpi_buffer(DECitem%force_scheme,Master)
    call ls_mpi_buffer(DECitem%dyn_load,Master)
    call ls_mpi_buffer(DECitem%print_frags,Master)
    call ls_mpi_buffer(DECitem%abc,Master)
    call ls_mpi_buffer(DECitem%ijk_tile_size,Master)
    call ls_mpi_buffer(DECitem%abc_tile_size,Master)
    call ls_mpi_buffer(DECitem%ijk_nbuffs,Master)
    call ls_mpi_buffer(DECitem%abc_nbuffs,Master)
    call ls_mpi_buffer(DECitem%acc_sync,Master)
    call ls_mpi_buffer(DECitem%pt_single_prec,Master)
    call ls_mpi_buffer(DECitem%pt_hack,Master)
    call ls_mpi_buffer(DECitem%pt_hack2,Master)
    call ls_mpi_buffer(DECitem%CCSDno_restart,Master)
    call ls_mpi_buffer(DECitem%CCSD_NO_DEBUG_COMM,Master)
    call ls_mpi_buffer(DECitem%spawn_comm_proc,Master)
    call ls_mpi_buffer(DECitem%CCSDpreventcanonical,Master)
    call ls_mpi_buffer(DECitem%NO_MO_CCSD,Master)
    call ls_mpi_buffer(DECitem%CCDhack,Master)
    call ls_mpi_buffer(DECitem%noPNOtrafo,Master)
    call ls_mpi_buffer(DECitem%noPNOtrunc,Master)
    call ls_mpi_buffer(DECitem%noFAtrafo,Master)
    call ls_mpi_buffer(DECitem%noFAtrunc,Master)
    call ls_mpi_buffer(DECitem%simplePNOthr,Master)
    call ls_mpi_buffer(DECitem%use_pnos,Master)
    call ls_mpi_buffer(DECitem%EOSPNOthr,Master)
    call ls_mpi_buffer(DECitem%noPNOoverlaptrunc,Master)
    call ls_mpi_buffer(DECitem%PNOoverlapthr,Master)
    call ls_mpi_buffer(DECitem%PNOtriangular,Master)
    call ls_mpi_buffer(DECitem%pno_S_on_the_fly,Master)
    call ls_mpi_buffer(DECitem%CCSDmultipliers,Master)
    call ls_mpi_buffer(DECitem%simple_multipler_residual,Master)
    call ls_mpi_buffer(DECitem%CRASHCALC,Master)
    call ls_mpi_buffer(DECitem%CRASHESTI,Master)
    call ls_mpi_buffer(DECitem%cc_driver_debug,Master)
    call ls_mpi_buffer(DECitem%use_bg_buffer,Master)
    call ls_mpi_buffer(DECitem%cc_solver_tile_mem,Master)
    call ls_mpi_buffer(DECitem%en_mem,Master)
    call ls_mpi_buffer(DECitem%ccsolver_overwrite_prec,Master)
    call ls_mpi_buffer(DECitem%precondition_with_full,Master)
    call ls_mpi_buffer(DECitem%ccsd_expl,Master)
    call ls_mpi_buffer(DECitem%ccMaxIter,Master)
    call ls_mpi_buffer(DECitem%ccMaxDIIS,Master)
    call ls_mpi_buffer(DECitem%ccConvergenceThreshold,Master)
    call ls_mpi_buffer(DECitem%CCthrSpecified,Master)
    call ls_mpi_buffer(DECitem%use_preconditioner,Master)
    call ls_mpi_buffer(DECitem%ccsolverskip,Master)
    call ls_mpi_buffer(DECitem%use_preconditioner_in_b,Master)
    call ls_mpi_buffer(DECitem%use_crop,Master)

    call ls_mpi_buffer(DECitem%F12,Master)
    call ls_mpi_buffer(DECitem%F12fragopt,Master)
    call ls_mpi_buffer(DECitem%F12DEBUG,Master)
    call ls_mpi_buffer(DECitem%F12Ccoupling,Master)
    call ls_mpi_buffer(DECitem%F12singles,Master)
    call ls_mpi_buffer(DECitem%F12singlesMaxIter,Master)
    call ls_mpi_buffer(DECitem%F12singlesThr,Master)
    call ls_mpi_buffer(DECitem%F12singlesMaxDIIS,Master)
    call ls_mpi_buffer(DECitem%NaturalLinearScalingF12Terms,Master)
    call ls_mpi_buffer(DECitem%NaturalLinearScalingF12TermsB1,Master)
    call ls_mpi_buffer(DECitem%NaturalLinearScalingF12TermsV1,Master)
    call ls_mpi_buffer(DECitem%NaturalLinearScalingF12TermsX1,Master)

    call ls_mpi_buffer(DECitem%PureHydrogenDebug,Master)
    call ls_mpi_buffer(DECitem%StressTest,Master)
    call ls_mpi_buffer(DECitem%AtomicExtent,Master)
    call ls_mpi_buffer(DECitem%MPMP2,Master)
    !RIMP2 settings
    call ls_mpi_buffer(DECitem%AuxAtomicExtent,Master)
    call ls_mpi_buffer(DECitem%NAF,Master)
    call ls_mpi_buffer(DECitem%NAFthreshold,Master)
    call ls_mpi_buffer(DECitem%RIMPSubGroupSize,Master)
    call ls_mpi_buffer(DECitem%RIMP2PDMTENSOR,Master)
    call ls_mpi_buffer(DECinfo%RIMP2ForcePDMCalpha,Master)
    call ls_mpi_buffer(DECinfo%RIMP2_tiling,Master)
    call ls_mpi_buffer(DECinfo%RIMP2_lowdin,Master)
    call ls_mpi_buffer(DECinfo%RIMP2_Laplace,Master)
    call ls_mpi_buffer(DECinfo%RIMP2_deactivateopenmp,Master)
    call ls_mpi_buffer(DECitem%DFTreference,Master)
    call ls_mpi_buffer(DECitem%mpisplit,Master)
    call ls_mpi_buffer(DECitem%rimpisplit,Master)
    call ls_mpi_buffer(DECitem%MPIgroupsize,Master)
    call ls_mpi_buffer(DECitem%manual_batchsizes,Master)
    call ls_mpi_buffer(DECitem%ccsdAbatch,Master)
    call ls_mpi_buffer(DECitem%ccsdGbatch,Master)
    call ls_mpi_buffer(DECitem%manual_occbatchsizes,Master)
    call ls_mpi_buffer(DECitem%batchOccI,Master)
    call ls_mpi_buffer(DECitem%batchOccJ,Master)
    call ls_mpi_buffer(DECitem%InteractionEnergy,Master)
    call ls_mpi_buffer(DECitem%PrintInteractionEnergy,Master)
    call ls_mpi_buffer(DECitem%hack,Master)
    call ls_mpi_buffer(DECitem%hack2,Master)
    call ls_mpi_buffer(DECitem%test_fully_distributed_integrals,Master)
    call ls_mpi_buffer(DECitem%SkipReadIn,Master)
    call ls_mpi_buffer(DECitem%tensor_test,Master)
    call ls_mpi_buffer(DECitem%tensor_segmenting_scheme,Master)
    call ls_mpi_buffer(DECitem%reorder_test,Master)
    call ls_mpi_buffer(DECitem%check_lcm_orbitals,Master)
    call ls_mpi_buffer(DECitem%check_Occ_SubSystemLocality,Master)
    call ls_mpi_buffer(DECitem%force_Occ_SubSystemLocality,Master)
    call ls_mpi_buffer(DECitem%PL,Master)
    call ls_mpi_buffer(DECitem%MemDebugPrint,Master)
    call ls_mpi_buffer(DECitem%print_small_calc,Master)
    call ls_mpi_buffer(DECitem%skipfull,Master)
    call ls_mpi_buffer(DECitem%distribute_fullmolecule,Master)
    call ls_mpi_buffer(DECitem%force_distribution,Master)
    call ls_mpi_buffer(DECitem%output,Master)
    call ls_mpi_buffer(DECitem%AbsorbHatoms,Master)
    call ls_mpi_buffer(DECitem%FitOrbitals,Master)
    call ls_mpi_buffer(DECitem%simple_orbital_threshold,Master)
    call ls_mpi_buffer(DECitem%purifyMOs,Master)
    call ls_mpi_buffer(DECitem%fragadapt,Master)
    call ls_mpi_buffer(DECitem%mulliken,Master)
    call ls_mpi_buffer(DECitem%distance,Master)
    call ls_mpi_buffer(DECitem%FOT,Master)
    call ls_mpi_buffer(DECitem%GeoFOTs,nFOTs,Master)
    call ls_mpi_buffer(DECitem%MaxIter,Master)
    call ls_mpi_buffer(DECitem%FOTlevel,Master)
    call ls_mpi_buffer(DECitem%Frag_Exp_Scheme,Master)
    call ls_mpi_buffer(DECitem%Frag_RedOcc_Scheme,Master)
    call ls_mpi_buffer(DECitem%Frag_RedVir_Scheme,Master)
    call ls_mpi_buffer(DECitem%Frag_Init_Size,Master)
    call ls_mpi_buffer(DECitem%Frag_Exp_Size,Master)
    call ls_mpi_buffer(DECitem%Frag_red1_thr,Master)
    call ls_mpi_buffer(DECitem%Frag_red2_thr,Master)
    call ls_mpi_buffer(DECitem%Frag_red_occ,Master)
    call ls_mpi_buffer(DECitem%Frag_red_virt,Master)
    call ls_mpi_buffer(DECitem%fragopt_exp_model,Master)
    call ls_mpi_buffer(DECitem%fragopt_red_model,Master)
    call ls_mpi_buffer(DECitem%OnlyOccPart,Master)
    call ls_mpi_buffer(DECitem%OnlyVirtPart,Master)
    call ls_mpi_buffer(DECitem%all_init_radius,Master)
    call ls_mpi_buffer(DECitem%occ_init_radius,Master)
    call ls_mpi_buffer(DECitem%vir_init_radius,Master)
    call ls_mpi_buffer(DECitem%RepeatAF,Master)
    call ls_mpi_buffer(DECitem%CorrDensScheme,Master)
    call ls_mpi_buffer(DECitem%use_abs_overlap,Master)
    call ls_mpi_buffer(DECitem%pairestimateignore,Master)
    call ls_mpi_buffer(DECitem%pair_distance_threshold,Master)
    call ls_mpi_buffer(DECitem%PairMinDist,Master)
    call ls_mpi_buffer(DECitem%checkpairs,Master)
    call ls_mpi_buffer(DECitem%pairFOthr,Master)
    call ls_mpi_buffer(DECitem%PairEstimate,Master)
    call ls_mpi_buffer(DECitem%EstimateInitRadius,Master)
    call ls_mpi_buffer(DECitem%EstimateInitAtom,Master)
    call ls_mpi_buffer(DECitem%PairEstimateModel,Master)
    call ls_mpi_buffer(DECitem%nFRAGSred,Master)
    call ls_mpi_buffer(DECitem%FOTscaling,Master)
    call ls_mpi_buffer(DECitem%first_order,Master)
    call ls_mpi_buffer(DECitem%density,Master)
    call ls_mpi_buffer(DECitem%unrelaxed,Master)
    call ls_mpi_buffer(DECitem%gradient,Master)
    call ls_mpi_buffer(DECitem%kappa_use_preconditioner,Master)
    call ls_mpi_buffer(DECitem%kappa_use_preconditioner_in_b,Master)
    call ls_mpi_buffer(DECitem%kappaMaxDIIS,Master)
    call ls_mpi_buffer(DECitem%kappaMaxIter,Master)
    call ls_mpi_buffer(DECitem%kappa_driver_debug,Master)
    call ls_mpi_buffer(DECitem%kappaTHR,Master)
    call ls_mpi_buffer(DECitem%SOS,Master)
    call ls_mpi_buffer(DECitem%FracOfOrbSpace_red,Master)
    call ls_mpi_buffer(DECitem%ncalc,nFOTs,Master)
    call ls_mpi_buffer(DECitem%EerrFactor,Master)
    call ls_mpi_buffer(DECitem%EerrOLD,Master)
    call ls_mpi_buffer(DECitem%no_pairs,Master)
    call ls_mpi_buffer(DECitem%only_pair_frag_jobs,Master)
    call ls_mpi_buffer(DECitem%only_n_frag_jobs,Master)
    if(DECitem%only_n_frag_jobs>0)then
       if(.not. AddToBuffer)then
          call mem_alloc(DECitem%frag_job_nr,DECitem%only_n_frag_jobs)
       endif
       call ls_mpi_buffer(DECitem%frag_job_nr,DECitem%only_n_frag_jobs,Master)
    endif
    call ls_mpi_buffer(DECitem%UseIchor,Master)
    call ls_mpi_buffer(DECitem%IntegralThreshold,Master)
    call ls_mpi_buffer(DECitem%noaofock,Master)
    call ls_mpi_buffer(DECitem%THCNOPRUN,Master)
    call ls_mpi_buffer(DECitem%THCDUMP,Master)
    call ls_mpi_buffer(DECitem%THCradint,Master)
    call ls_mpi_buffer(DECitem%THC_MIN_RAD_PT,Master)
    call ls_mpi_buffer(DECitem%THCangint,Master)
    call ls_mpi_buffer(DECitem%THCHRDNES,Master)
    call ls_mpi_buffer(DECitem%THCTURBO,Master)
    call ls_mpi_buffer(DECitem%THCRADIALGRID,Master)
    call ls_mpi_buffer(DECitem%THCZdependenMaxAng,Master)
    call ls_mpi_buffer(DECitem%THCPARTITIONING,Master)

  end subroutine mpicopy_dec_settings

!  subroutine rpa_res_communicate_data(gmo,t2,omega2,nvirt,nocc)
!    implicit none
!    !real(realk),intent(inout),pointer :: gmo(:)
!    !type(array4), intent(inout) :: omega2
!    type(tensor), intent(inout) :: gmo
!    type(tensor), intent(inout) :: omega2
!    !type(array4),intent(inout)         :: t2
!    type(tensor),intent(inout)         :: t2
!    !real(realk),intent(inout)         :: t2(:,:,:,:)
!    integer,intent(inout)             :: nvirt,nocc
!    logical :: master
!    integer :: addr1(infpar%lg_nodtot)
!    integer :: addr2(infpar%lg_nodtot)
!    integer :: addr3(infpar%lg_nodtot)
!
!
!    master = (infpar%lg_mynum == infpar%master)
!
!    if(master) then
!   !   write(*,*)'Johannes addr in comm', omega2%addr_p_arr
!      call ls_mpibcast(RPAGETRESIDUAL,infpar%master,infpar%lg_comm)
!      addr1 = omega2%addr_p_arr
!      addr2 = gmo%addr_p_arr
!      addr3 = t2%addr_p_arr
!    endif
!
!    call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
!    call ls_mpi_buffer(nvirt,infpar%master)
!    call ls_mpi_buffer(nocc,infpar%master)
!    call ls_mpi_buffer(addr1,infpar%lg_nodtot,infpar%master)
!    call ls_mpi_buffer(addr2,infpar%lg_nodtot,infpar%master)
!    call ls_mpi_buffer(addr3,infpar%lg_nodtot,infpar%master)
!
!
!    if(.not.master)then
!     ! call mem_alloc(gmo,nvirt*nocc*nocc*nvirt)
!      omega2 = get_tensor_from_parr(addr1(infpar%lg_mynum+1))
!      gmo     = get_tensor_from_parr(addr2(infpar%lg_mynum+1))
!      t2     = get_tensor_from_parr(addr3(infpar%lg_mynum+1))
!      !t2=array4_init([nvirt,nocc,nvirt,nocc])
!      !omega2=array4_init([nvirt,nocc,nvirt,nocc])
!      !omega2=tensor_ainit([nvirt,nvirt,nocc,nocc],4,atype='TDAR')
!      !call tensor_ainit(gmo,[nvirt,nvirt,nocc,nocc],4,local =.true.,atype='TDAR')
!      !call tensor_ainit(t2, [nvirt,nvirt,nocc,nocc],4,local =.true.,atype='TDAR')
!    endif
!    !call ls_mpi_buffer(gmo,nvirt*nocc*nocc*nvirt,infpar%master)
!    !call ls_mpibcast(t2,nvirt,nvirt,nocc,nocc,infpar%master,infpar%lg_comm)
!    call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
!
!    !print*,' inside rpa_res_comm',infpar%lg_mynum
!
!    
!    !print*,' after gmo',infpar%lg_mynum
!
!    !print*,' after omega2',infpar%lg_mynum
!
!
!    !call ls_mpibcast(t2%elm4,nvirt,nvirt,nocc,nocc,infpar%master,infpar%lg_comm)
!    !print*,' after t2',infpar%lg_mynum
!    !call ls_mpibcast(gmo%elm1,nvirt*nvirt*nocc*nocc,infpar%master,infpar%lg_comm)
!
!
!    !call ls_mpibcast(omega2%elm4,nvirt,nvirt,nocc,nocc,infpar%master,infpar%lg_comm)
!
!    !print*,' after finalize',infpar%lg_mynum
!    !call ls_mpibcast(t2%elm1,nvirt*nvirt*nocc*nocc,infpar%master,infpar%lg_comm)
!
!
!  end subroutine rpa_res_communicate_data

  subroutine rpa_residual_communicate_data(t2,omega2,iajb,oof,vvf,no,nv)
    implicit none
    type(tensor),intent(inout) :: omega2,iajb,t2,oof,vvf
    logical :: master
    integer :: addr1(infpar%lg_nodtot)
    integer :: addr2(infpar%lg_nodtot)
    integer :: addr3(infpar%lg_nodtot)
    integer :: addr4(infpar%lg_nodtot)
    integer :: addr5(infpar%lg_nodtot)
    !integer :: iter
    !type(array4), intent(inout) :: omega2
    !type(tensor), intent(inout) :: omega2
    !type(tensor),intent(inout)         :: t2
    !type(tensor), intent(inout) :: pfock,qfock
    integer,intent(inout)             :: nv,no
    !integer :: addr1(infpar%lg_nodtot)
    !integer :: addr2(infpar%lg_nodtot)
    !integer :: addr3(infpar%lg_nodtot)
    !integer :: addr4(infpar%lg_nodtot)
    !!real(realk), intent(inout) :: pfock(no,no),qfock(nv,nv)
    !logical :: master


    master = (infpar%lg_mynum == infpar%master)


    if(master) then
      call ls_mpibcast(RPAGETFOCK,infpar%master,infpar%lg_comm)
      addr1 = omega2%addr_p_arr
      addr2 = iajb%addr_p_arr
      addr3 = t2%addr_p_arr
      addr4 = oof%addr_p_arr
      addr5 = vvf%addr_p_arr
    endif


    call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
    call ls_mpi_buffer(nv,infpar%master)
    call ls_mpi_buffer(no,infpar%master)
    call ls_mpi_buffer(addr1,infpar%lg_nodtot,infpar%master)
    call ls_mpi_buffer(addr2,infpar%lg_nodtot,infpar%master)
    call ls_mpi_buffer(addr3,infpar%lg_nodtot,infpar%master)
    call ls_mpi_buffer(addr4,infpar%lg_nodtot,infpar%master)
    call ls_mpi_buffer(addr5,infpar%lg_nodtot,infpar%master)
    call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)

    if(.not.master)then
      omega2 = get_tensor_from_parr(addr1(infpar%lg_mynum+1))
      iajb   = get_tensor_from_parr(addr2(infpar%lg_mynum+1))
      t2     = get_tensor_from_parr(addr3(infpar%lg_mynum+1))
      oof    = get_tensor_from_parr(addr4(infpar%lg_mynum+1))
      vvf    = get_tensor_from_parr(addr5(infpar%lg_mynum+1))
    endif
    !call ls_mpibcast(t2%elm1,nv*nv*no*no,infpar%master,infpar%lg_comm)


  end subroutine rpa_residual_communicate_data




  subroutine rpa_fock_communicate_data(t2,omega2,iajb,oof,vvf,no,nv)
    implicit none
    type(tensor),intent(inout) :: omega2,iajb,t2,oof,vvf
    logical :: master
    integer :: addr1(infpar%lg_nodtot)
    integer :: addr2(infpar%lg_nodtot)
    integer :: addr3(infpar%lg_nodtot)
    integer :: addr4(infpar%lg_nodtot)
    integer :: addr5(infpar%lg_nodtot)
    !integer :: iter
    !type(array4), intent(inout) :: omega2
    !type(tensor), intent(inout) :: omega2
    !type(tensor),intent(inout)         :: t2
    !type(tensor), intent(inout) :: pfock,qfock
    integer,intent(inout)             :: nv,no
    !integer :: addr1(infpar%lg_nodtot)
    !integer :: addr2(infpar%lg_nodtot)
    !integer :: addr3(infpar%lg_nodtot)
    !integer :: addr4(infpar%lg_nodtot)
    !!real(realk), intent(inout) :: pfock(no,no),qfock(nv,nv)
    !logical :: master


    master = (infpar%lg_mynum == infpar%master)


    if(master) then
      call ls_mpibcast(RPAGETFOCK,infpar%master,infpar%lg_comm)
      addr1 = omega2%addr_p_arr
      addr2 = iajb%addr_p_arr
      addr3 = t2%addr_p_arr
      addr4 = oof%addr_p_arr
      addr5 = vvf%addr_p_arr
    endif


    call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
    call ls_mpi_buffer(nv,infpar%master)
    call ls_mpi_buffer(no,infpar%master)
    call ls_mpi_buffer(addr1,infpar%lg_nodtot,infpar%master)
    call ls_mpi_buffer(addr2,infpar%lg_nodtot,infpar%master)
    call ls_mpi_buffer(addr3,infpar%lg_nodtot,infpar%master)
    call ls_mpi_buffer(addr4,infpar%lg_nodtot,infpar%master)
    call ls_mpi_buffer(addr5,infpar%lg_nodtot,infpar%master)
    call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)

    if(.not.master)then
      omega2 = get_tensor_from_parr(addr1(infpar%lg_mynum+1))
      iajb   = get_tensor_from_parr(addr2(infpar%lg_mynum+1))
      t2     = get_tensor_from_parr(addr3(infpar%lg_mynum+1))
      oof    = get_tensor_from_parr(addr4(infpar%lg_mynum+1))
      vvf    = get_tensor_from_parr(addr5(infpar%lg_mynum+1))
    endif
    !call ls_mpibcast(t2%elm1,nv*nv*no*no,infpar%master,infpar%lg_comm)


  end subroutine rpa_fock_communicate_data



  !> \brief bcast very basic information from master to slaves 
  !> (information which for practical reasons cannot be packed into mpi_dec_fullinfo_master_to_slaves)
  subroutine mpi_dec_fullinfo_master_to_slaves_precursor(esti,nocc,nvirt,master)
    implicit none
    !> Is this an estimated calculation (no FOT optimization)?
    logical,intent(inout) :: esti
    !> Number of occ/virt orbitals in molecule
    integer,intent(inout) :: nocc,nvirt
    !> Master node number
    integer(kind=ls_mpik),intent(in) :: master

    call ls_mpibcast(esti,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(nocc,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(nvirt,master,MPI_COMM_LSDALTON)

  end subroutine mpi_dec_fullinfo_master_to_slaves_precursor

  subroutine wake_slaves_for_simple_mo(integral,trafo1,trafo2,trafo3,trafo4,mylsitem,c)
     implicit none
     type(tensor),intent(inout)   :: integral
     type(tensor),intent(inout)   :: trafo1,trafo2,trafo3,trafo4
     type(lsitem), intent(inout) :: mylsitem
     logical, intent(inout) :: c
     integer :: addr1(infpar%lg_nodtot)
     integer :: addr2(infpar%lg_nodtot)
     integer :: addr3(infpar%lg_nodtot)
     integer :: addr4(infpar%lg_nodtot)
     integer :: addr5(infpar%lg_nodtot)
     logical :: master


     master = (infpar%lg_mynum == infpar%master)

     if(master) call ls_mpibcast(MO_INTEGRAL_SIMPLE,infpar%master,infpar%lg_comm)

     if(master)then
        addr1 = trafo1%addr_p_arr
        addr2 = trafo2%addr_p_arr
        addr3 = trafo3%addr_p_arr
        addr4 = trafo4%addr_p_arr
        addr5 = integral%addr_p_arr
     endif

     call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
     call ls_mpi_buffer(addr1,infpar%lg_nodtot,infpar%master)
     call ls_mpi_buffer(addr2,infpar%lg_nodtot,infpar%master)
     call ls_mpi_buffer(addr3,infpar%lg_nodtot,infpar%master)
     call ls_mpi_buffer(addr4,infpar%lg_nodtot,infpar%master)
     call ls_mpi_buffer(addr5,infpar%lg_nodtot,infpar%master)
     call mpicopy_lsitem(MyLsItem,infpar%lg_comm)
     call ls_mpi_buffer(c,infpar%master)
     call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)

     if(.not.master)then
        trafo1 = get_tensor_from_parr(addr1(infpar%lg_mynum+1))
        trafo2 = get_tensor_from_parr(addr2(infpar%lg_mynum+1))
        trafo3 = get_tensor_from_parr(addr3(infpar%lg_mynum+1))
        trafo4 = get_tensor_from_parr(addr4(infpar%lg_mynum+1))
        integral = get_tensor_from_parr(addr5(infpar%lg_mynum+1))
     endif


  end subroutine wake_slaves_for_simple_mo


  subroutine get_slaves_to_simple_par_mp2_res(omega2,iajb,t2,oof,vvf,iter)
     implicit none
     type(tensor),intent(inout) :: omega2,iajb,t2,oof,vvf
     logical :: master
     integer :: addr1(infpar%lg_nodtot)
     integer :: addr2(infpar%lg_nodtot)
     integer :: addr3(infpar%lg_nodtot)
     integer :: addr4(infpar%lg_nodtot)
     integer :: addr5(infpar%lg_nodtot)
     integer :: iter


     master = (infpar%lg_mynum == infpar%master)

     if(master) call ls_mpibcast(SIMPLE_MP2_PAR,infpar%master,infpar%lg_comm)

     if(master)then
        addr1 = omega2%addr_p_arr
        addr2 = iajb%addr_p_arr
        addr3 = t2%addr_p_arr
        addr4 = oof%addr_p_arr
        addr5 = vvf%addr_p_arr
     endif

     call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
     call ls_mpi_buffer(addr1,infpar%lg_nodtot,infpar%master)
     call ls_mpi_buffer(addr2,infpar%lg_nodtot,infpar%master)
     call ls_mpi_buffer(addr3,infpar%lg_nodtot,infpar%master)
     call ls_mpi_buffer(addr4,infpar%lg_nodtot,infpar%master)
     call ls_mpi_buffer(addr5,infpar%lg_nodtot,infpar%master)
     call ls_mpi_buffer(iter,infpar%master)
     call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)

     if(.not.master)then
        omega2 = get_tensor_from_parr(addr1(infpar%lg_mynum+1))
        iajb   = get_tensor_from_parr(addr2(infpar%lg_mynum+1))
        t2     = get_tensor_from_parr(addr3(infpar%lg_mynum+1))
        oof    = get_tensor_from_parr(addr4(infpar%lg_mynum+1))
        vvf    = get_tensor_from_parr(addr5(infpar%lg_mynum+1))
     endif

  end subroutine get_slaves_to_simple_par_mp2_res

#endif

  !> \brief this routine is used to decide wheter a job has to be carried out on
  !a given node and enables a dynamic load balance for the nodes in an integral
  !direct loop
  !> \author Patrick Ettenhuber
  !> \date unknown
  subroutine check_job(batch,fr,dyn,a,g,nA,nG,static,win,prnt)
     implicit none
     !> combined index in the alpha-gamma counter regime
     integer, intent(inout) :: batch
     !> logical to check whether this routine is called for the first time-this
     !is important for the dynamic load balancing
     logical, intent(inout) :: fr
     !> is dynamic load balancing requested, is print requested
     logical,intent(in) :: dyn,prnt
     !> number of alpha and gamma batches
     integer,intent(in) :: nA,nG
     !> current alpha and gamma counters
     integer,intent(inout) :: a,g
     ! if not dynamic lload balancing this array contains the ranks of the nodes
     ! that have to do a particular job, with gamma as the fast counter!!
     integer, pointer, intent(inout) :: static(:)
     ! this is the window handle used in the dynamic load balancing scheme
     integer(kind=ls_mpik) :: win
     !> internal variables
     real(realk) :: mpi_buf
     integer :: GammaAlpha(2)
     integer :: one
#ifdef VAR_MPI
     one  = 1

     !GET MPI batch number
     if(.not.dyn)then

        !Check the static array for my_rank and return
        batch=batch+1
        FindJobLoop: do while(batch<=nA*nG)

           call get_midx(batch,GammaAlpha,[nA,nG],2)

           a = GammaAlpha(1)
           g = GammaAlpha(2)

           !check whether this job has been assigned to me
           if(static((a-1)*nG+g)/=infpar%lg_mynum)then
              batch=batch+1
           else

              exit FindJobLoop

           endif

        enddo FindJobLoop

     else

        if(fr)then

           !in the first call we just take the rank as the job to be carried out
           fr   = .false.
           batch = infpar%lg_mynum + 1

        else

           ! in all subsequent calls we use lsmpi_get_acc to get the new job
           ! identifier
           call time_start_phase( PHASE_COMM )
#ifdef VAR_HAVE_MPI3
           call lsmpi_get_acc(one,batch,infpar%master,1,win)
           call lsmpi_win_flush(win,rank = infpar%master, local=.true.)
#else
           call lsmpi_win_lock(infpar%master,win,'e')
           call lsmpi_get_acc(one,batch,infpar%master,1,win)
           call lsmpi_win_unlock(infpar%master,win)
#endif
           call time_start_phase( PHASE_WORK )
        endif
     endif

#else
     !Non--MPI job counting, only increment
     batch = batch + 1
#endif

     ! No more jobs to be done exit
     if(batch > nG*nA )return

     call get_midx(batch,GammaAlpha,[nA,nG],2)

     a = GammaAlpha(1)
     g = GammaAlpha(2)

#ifdef VAR_MPI
     if(prnt) write (*, '("Rank ",I3," starting job (",I3,"/",I3,",",I3,"/",I3,")")') infpar%mynum,&
        &a,na,g,ng
#else
     if(prnt) write (*, '("starting job (",I3,"/",I3,",",I3,"/",I3,")")')a,&
        &na,g,ng
#endif
     call lsmpi_poke()

  end subroutine check_job

end module decmpi_module



#ifdef VAR_MPI
subroutine set_dec_settings_on_slaves()
   use infpar_module
   use lsmpi_type
   use lsparameters
   use dec_typedef_module
   use decmpi_module, only:mpibcast_dec_settings
   implicit none
   if(infpar%mynum == infpar%master) call ls_mpibcast(DEC_SETTING_TO_SLAVES,infpar%master,MPI_COMM_LSDALTON)
   call mpibcast_dec_settings(DECinfo,MPI_COMM_LSDALTON)
end subroutine set_dec_settings_on_slaves

#endif
