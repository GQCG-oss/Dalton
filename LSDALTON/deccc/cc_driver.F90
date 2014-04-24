!> @file
!> Main CC driver (mp2, cc2 and ccsd are so far implemented)
!> \author Marcin Ziolkowski @ AU 2009,2010 teozio(at)gmail.com

!> General coupled-cluster solver for both full molecular and dec
!> calculations. This should depend only on the integral program and classes
!> related to storage of two- and four-dimensional arrays. All other parameters
!> should be passed as parameters.
module ccdriver

use precision
use lstiming!, only: lstimer
use typedeftype!,only:lsitem
use typedef
use files!,only:lsopen,lsclose
use memory_handling
use dec_typedef_module
use integralinterfaceMod
use integralparameters
#ifdef VAR_MPI
use infpar_module
#endif
use tensor_interface_module


! DEC DEPENDENCIES (within deccc directory)   
! *****************************************
use dec_fragment_utils!,only: get_density_from_occ_orbitals
use crop_tools_module
use array2_simple_operations
use array4_simple_operations
use mp2_module!,only: get_VOVO_integrals
use atomic_fragment_operations
use ccintegrals!,only:get_full_eri,getL_simple_from_gmo,&
!       & get_gmo_simple,get_h1
use cc_debug_routines_module
use ccsd_module!,only: getDoublesResidualMP2_simple, &
!       & getDoublesResidualCCSD_simple,getDoublesResidualCCSD_simple2, &
!       & precondition_doubles,get_ccsd_residual_integral_driven,&
!       & get_ccsd_residual_integral_driven_oldarray_wrapper
#ifdef MOD_UNRELEASED
use ccsdpt_module
!endif mod_unreleased
#endif
use orbital_operations
use rpa_module



public :: ccsolver, ccsolver_par, fragment_ccsolver, ccsolver_justenergy,&
   & mp2_solver
private

contains



!!> \brief Get coupled-cluster energy by calling general ccsolver.
!!> \author Kasper Kristensen
!function ccsolver_justenergy(ccmodel,MyMolecule,Co_f,Cv_f,fock_f,nbasis,nocc,nvirt, &
!      &mylsitem,ccPrintLevel,fragment_job,ppfock_f,qqfock_f) result(ccenergy)
!
!   implicit none
!
!   !> CC model
!   integer,intent(inout) :: ccmodel
!   !> full molecule information
!   type(fullmolecule), intent(in) :: MyMolecule
!   !> Number of occupied orbitals in full molecule/fragment AOS
!   !> (only number of valence orbitals for frozen core)
!   integer, intent(in) :: nocc
!   !> Number of virtual orbitals in full molecule/fragment AOS
!   integer, intent(in) :: nvirt
!   !> Number of basis functions in full molecule/atomic extent
!   integer, intent(in) :: nbasis
!   !> Fock matrix in AO basis for fragment or full molecule
!   real(realk), dimension(nbasis,nbasis), intent(in) :: fock_f
!   !> Occupied MO coefficients  for fragment/full molecule
!   real(realk), dimension(nbasis,nocc), intent(inout) :: Co_f
!   !> Virtual MO coefficients  for fragment/full molecule
!   real(realk), dimension(nbasis,nvirt), intent(inout) :: Cv_f
!   !> Occ-occ block of Fock matrix in MO basis
!   real(realk), dimension(nocc,nocc), intent(inout) :: ppfock_f
!   !> Virt-virt block of Fock matrix in MO basis
!   real(realk), dimension(nvirt,nvirt), intent(inout) :: qqfock_f
!   !> Is this a fragment job (true) or a full molecular calculation (false)
!   logical, intent(in) :: fragment_job
!   !> LS item information
!   type(lsitem), intent(inout) :: mylsitem
!   !> How much to print? ( ccPrintLevel>0 --> print info stuff)
!   integer, intent(in) :: ccPrintLevel
!   !> Coupled cluster energy for fragment/full molecule
!   real(realk) :: ccenergy!,ccsdpt_e4,ccsdpt_e5,ccsdpt_tot
!   type(array4) :: t2_final,VOVO, mp2_amp!,ccsdpt_t2
!   type(array2) :: t1_final!,ccsdpt_t1
!   logical :: local
!
!   local=.true.
!#ifdef VAR_MPI
!   if(infpar%lg_nodtot>1.or.DECinfo%hack2)local=.false.
!#endif
!
!
!   if(DECinfo%CCDEBUG)then
!
!      if(DECinfo%use_pnos)then
!
!         !GET MP2 AMPLITUDES TO CONSTRUCT PNOS
!         call get_VOVO_integrals( mylsitem, nbasis, nocc, nvirt, Cv_f, Co_f, VOVO )
!         call mp2_solver( nocc, nvirt, ppfock_f, qqfock_f, VOVO, mp2_amp )
!         call array4_free( VOVO )
!
!         !CALL THE SOLVER WITH PNO ARGUMENT
!         call ccsolver_debug(ccmodel,Co_f,Cv_f,fock_f,nbasis,nocc,nvirt, &
!            & mylsitem,ccPrintLevel,fragment_job,ppfock_f,qqfock_f,ccenergy, &
!            & t1_final,t2_final,VOVO,.false.,SOLVE_AMPLITUDES,m2=mp2_amp,use_pnos=DECinfo%use_pnos)
!
!         !FREE MP2 AMPLITUDES
!         call array4_free( mp2_amp )
!
!      else
!
!         !CALL DEBUG SOLVER WITHOUT PNOS
!         call ccsolver_debug(ccmodel,Co_f,Cv_f,fock_f,nbasis,nocc,nvirt, &
!            & mylsitem,ccPrintLevel,fragment_job,ppfock_f,qqfock_f,ccenergy, &
!            & t1_final,t2_final,VOVO,.false.,SOLVE_AMPLITUDES)
!
!      endif
!   else
!
!      ! CALL PRODUCTION SOLVER
!      call ccsolver_par(ccmodel,Co_f,Cv_f,fock_f,nbasis,nocc,nvirt, &
!         & mylsitem,ccPrintLevel,fragment_job,ppfock_f,qqfock_f,ccenergy, &
!         & t1_final,t2_final,VOVO,.false.,local)
!
!   endif
!
!
!   ! Print fragment energies (currently only for occupied partitioning scheme)
!   if(DECinfo%full_molecular_cc .and. DECinfo%full_print_frag_energies) then
!      call fragment_energies_in_fulL_ccsd_occ(MyMolecule,mylsitem,t1_final,t2_final,VOVO) 
!   end if
!
!   ! Free arrays
!   call array2_free(t1_final)
!   call array4_free(t2_final)
!   call array4_free(VOVO)
!
!end function ccsolver_justenergy


!!> \brief Calculate and print CCSD (or CCD) fragment energies from amplitudes
!!> for the full molecular system. Only intended for testing purposes.
!!> Only for occupied partitioning scheme at this stage.
!!> \author Kasper Kristensen
!subroutine fragment_energies_in_fulL_ccsd_occ(MyMolecule,mylsitem,t1,t2,VOVO) 
!
!   implicit none
!
!   !> full molecule information
!   type(fullmolecule), intent(in) :: MyMolecule
!   !> LS item information
!   type(lsitem), intent(inout) :: mylsitem
!   !> Singles amplitudes (not changed but need to be inout)
!   type(array2),intent(inout) :: t1
!   !> Doubles amplitudes (not changed but need to be inout)
!   type(array4),intent(inout) :: t2
!   !> Integrals (ai|bj) stored as (a,i,b,j) (not changed but need to be inout)
!   type(array4),intent(inout) :: VOVO
!   integer :: natoms,nocc_tot,ncore,i,p,pdx,nvirt,nocc
!   type(decorbital), pointer :: occorbitals(:)
!   type(decorbital), pointer :: unoccorbitals(:)
!   logical, pointer :: orbitals_assigned(:)
!   type(array2) :: ccsd_mat_tot,ccsd_mat_tmp
!   real(realk) :: interactionE
!
!
!#ifdef MOD_UNRELEASED
!   natoms = MyMolecule%natoms
!   ! Note: For frozen core approx: nocc_tot = nocc + ncore,  nocc=#valence orbitals
!   !       Without frozen core approx: nocc_tot = nocc
!   nocc_tot = MyMolecule%nocc
!   if(DECinfo%frozencore) then
!      ncore = MyMolecule%ncore
!      nocc = MyMolecule%nval
!   else
!      ncore = 0
!      nocc = nocc_tot
!   endif
!   nvirt = MyMolecule%nunocc
!
!   ! -- Analyze basis and create orbitals
!   call mem_alloc(occorbitals,nocc_tot)
!   call mem_alloc(unoccorbitals,nvirt)
!   call GenerateOrbitals_driver(MyMolecule,mylsitem,nocc_tot,nvirt,natoms, &
!      & occorbitals,unoccorbitals)
!
!   ! Orbital assignment
!   call mem_alloc(orbitals_assigned,natoms)
!   orbitals_assigned=.false.
!   do p=1,nocc_tot
!      pdx = occorbitals(p)%centralatom
!      orbitals_assigned(pdx) = .true.
!   end do
!   do p=1,nvirt
!      pdx = unoccorbitals(p)%centralatom
!      orbitals_assigned(pdx) = .true.
!   end do
!
!   ! reorder VOVO integrals from (a,i,b,j) to (a,b,i,j)
!   call array4_reorder(VOVO,[1,3,2,4])
!   ! reorder doubles amplitudes from (a,i,b,j) to (a,b,i,j)
!   call array4_reorder(t2,[1,3,2,4])
!
!   ! Calculate and print out ccsd fragment and pair interaction energies
!   ccsd_mat_tot = array2_init([natoms,natoms])
!   ccsd_mat_tmp = array2_init([natoms,natoms])
!   call ccsd_energy_full_occ(nocc,nvirt,natoms,ncore,t2,t1,VOVO,occorbitals,&
!      & ccsd_mat_tot%val,ccsd_mat_tmp%val)
!   call print_ccsd_full_occ(natoms,ccsd_mat_tot%val,orbitals_assigned,mymolecule%distancetable)
!   IF(DECinfo%InteractionEnergy)THEN
!      call lsquit('InteractionEnergy not implemented in full ccsd',-1)
!   ENDIF
!   IF(DECinfo%PrintInteractionEnergy)THEN
!      call add_dec_interactionenergies(natoms,ccsd_mat_tot%val,orbitals_assigned,&
!         & interactionE,MyMolecule%SubSystemIndex,2)
!      write(DECinfo%output,'(1X,a,g20.10)') 'CCSD Interaction correlation energy : ',interactionE
!   ENDIF
!
!   ! Delete orbitals 
!   do i=1,nocc_tot
!      call orbital_free(OccOrbitals(i))
!   end do
!   do i=1,nvirt
!      call orbital_free(UnoccOrbitals(i))
!   end do
!   call mem_dealloc(occorbitals)
!   call mem_dealloc(unoccorbitals)
!   call mem_dealloc(orbitals_assigned)
!   call array2_free(ccsd_mat_tot)
!   call array2_free(ccsd_mat_tmp)
!
!   ! reorder VOVO integrals back: from (a,b,i,j) to (a,i,b,j)
!   call array4_reorder(VOVO,[1,3,2,4])
!   ! reorder doubles amplitudes back: from (a,b,i,j) to (a,i,b,j)
!   call array4_reorder(t2,[1,3,2,4])
!#endif
!
!end subroutine fragment_energies_in_fulL_ccsd_occ




!> \brief get ccsd(t) corrections for full molecule.
!> \author Janus Juul Eriksen, modified by Patrick Ettenhuber and TK
!> \date February 2013
function ccsolver_justenergy(ccmodel,MyMolecule,nbasis,nocc,nvirt,mylsitem,&
      & ccPrintLevel,fragment_job,Co_fc,ppfock_fc) result(ccenergy)

   implicit none

   !> CC model
   integer,intent(inout) :: ccmodel
   !> full molecule information
   type(fullmolecule), intent(in) :: MyMolecule
   !> Number of occupied orbitals in full molecule/fragment AOS
   integer, intent(in) :: nocc
   !> Number of virtual orbitals in full molecule/fragment AOS
   integer, intent(in) :: nvirt
   !> Number of basis functions in full molecule/atomic extent
   integer, intent(in) :: nbasis
   !> Is this a fragment job (true) or a full molecular calculation (false)
   logical, intent(in) :: fragment_job
   !> LS item information
   type(lsitem), intent(inout) :: mylsitem
   !> How much to print? ( ccPrintLevel>0 --> print info stuff)
   integer, intent(in) :: ccPrintLevel
   !> Occupied MO coefficients  for fragment/full molecule (only used for Frozen core)
   real(realk), dimension(nbasis,nocc), intent(in),optional :: Co_fc
   !> Occ-occ block of Fock matrix in MO basis (only used for frozen core)
   real(realk), dimension(nocc,nocc), intent(in),optional :: ppfock_fc
   !> Coupled cluster energy for full molecule
   real(realk) :: ccenergy,ccsdpt_e4,ccsdpt_e5,ccsdpt_tot
   type(array4) :: t2_final,ccsdpt_t2,VOVO
   type(array2) :: t1_final,ccsdpt_t1,ccsd_mat_tot,ccsd_mat_tmp,e4_mat_tot,e4_mat_tmp,e5_mat_tot
   integer :: natoms,ncore,nocc_tot,p,pdx,i
   type(decorbital), pointer :: occ_orbitals(:)
   type(decorbital), pointer :: unocc_orbitals(:)
   logical, pointer :: orbitals_assigned(:)
   logical :: local
   real(realk) :: InteractionECCSD,InteractionECCSDPT4,InteractionECCSDPT5

   real(realk) :: time_CCSD_work, time_CCSD_comm, time_CCSD_idle
   real(realk) :: time_pT_work, time_pT_comm, time_pT_idle

   call time_start_phase(PHASE_WORK, swwork = time_CCSD_work, swcomm = time_CCSD_comm, swidle = time_CCSD_idle) 

   ccenergy   = 0.0E0_realk
   ccsdpt_e4  = 0.0E0_realk
   ccsdpt_e5  = 0.0E0_realk
   ccsdpt_tot = 0.0E0_realk

   local=.true.
#ifdef VAR_MPI
   if(infpar%lg_nodtot>1)local=.false.
#endif

   ! is this a frozen core calculation or not?
   if (DECinfo%frozencore) then

      ncore = MyMolecule%ncore

      if(.not. present(Co_fc)) then
         call lsquit('ccsolver_justenergy_pt: Occ MOs not present for frozencore!',-1)
      end if

      if(.not. present(ppfock_fc)) then
         call lsquit('ccsolver_justenergy_pt: Occ-occ Fock matrix not present for frozencore!',-1)
      end if

      if (DECinfo%CCDEBUG) then
         call ccsolver_debug(ccmodel,Co_fc,MyMolecule%Cv,MyMolecule%fock,nbasis,nocc,nvirt,&
            & mylsitem,ccPrintLevel,fragment_job,ppfock_fc,MyMolecule%qqfock,ccenergy,&
            & t1_final,t2_final,VOVO,.false.,SOLVE_AMPLITUDES)
      else
         call ccsolver_par(ccmodel,Co_fc,MyMolecule%Cv,MyMolecule%fock,nbasis,nocc,nvirt,&
            & mylsitem,ccPrintLevel,fragment_job,ppfock_fc,MyMolecule%qqfock,ccenergy,&
            & t1_final,t2_final,VOVO,.false.,local)
      endif

   else
      ncore = 0

      if (DECinfo%CCDEBUG) then
         call ccsolver_debug(ccmodel,MyMolecule%Co,MyMolecule%Cv,MyMolecule%fock,nbasis,nocc,nvirt,&
            & mylsitem,ccPrintLevel,fragment_job,MyMolecule%ppfock,MyMolecule%qqfock,ccenergy,&
            & t1_final,t2_final,VOVO,.false.,SOLVE_AMPLITUDES)
      else
         call ccsolver_par(ccmodel,MyMolecule%Co,MyMolecule%Cv,MyMolecule%fock,nbasis,nocc,nvirt,&
            & mylsitem,ccPrintLevel,fragment_job,MyMolecule%ppfock,MyMolecule%qqfock,ccenergy,&
            & t1_final,t2_final,VOVO,.false.,local)
      endif

   end if

   if(DECinfo%PL>1)then
      call time_start_phase(PHASE_WORK, dwwork = time_CCSD_work, dwcomm = time_CCSD_comm, dwidle = time_CCSD_idle, &
         &swwork = time_pT_work, swcomm = time_pT_comm, swidle = time_pT_idle, &
         &labeldwwork = 'MASTER WORK CC solver: ',&
         &labeldwcomm = 'MASTER COMM CC solver: ',&
         &labeldwidle = 'MASTER IDLE CC solver: ') 
   endif

   natoms = MyMolecule%natoms
   nocc_tot = MyMolecule%nocc

   if(ccmodel == MODEL_CCSDpT)then

      ccsdpt_t1 = array2_init([nvirt,nocc])
      ccsdpt_t2 = array4_init([nvirt,nvirt,nocc,nocc])

      if(DECinfo%frozencore) then
         call ccsdpt_driver(nocc,nvirt,nbasis,ppfock_fc,MyMolecule%qqfock,Co_fc,MyMolecule%Cv,mylsitem,t2_final,&
            & ccsdpt_t1,ccsdpt_t2)
      else
         call ccsdpt_driver(nocc,nvirt,nbasis,MyMolecule%ppfock,MyMolecule%qqfock,MyMolecule%Co,&
            & MyMolecule%Cv,mylsitem,t2_final,ccsdpt_t1,ccsdpt_t2)
      end if

      if(DECinfo%PL>1)then
         call time_start_phase(PHASE_WORK,dwwork = time_pT_work, dwcomm = time_pT_comm, dwidle = time_pT_idle, &
            &labeldwwork = 'MASTER WORK pT: ',&
            &labeldwcomm = 'MASTER COMM pT: ',&
            &labeldwidle = 'MASTER IDLE pT: ') 
      endif

   else

      call array4_reorder(t2_final,[1,3,2,4])

   endif

   ! as we want to  print out fragment and pair interaction fourth-order energy contributions,
   ! then for locality analysis purposes we need occ_orbitals and
   ! unocc_orbitals (adapted from fragment_energy.f90)

   ! -- Analyze basis and create orbitals
   call mem_alloc(occ_orbitals,nocc_tot)
   call mem_alloc(unocc_orbitals,nvirt)
   call GenerateOrbitals_driver(MyMolecule,mylsitem,nocc_tot,nvirt,natoms, &
      & occ_orbitals,unocc_orbitals)

   ! Orbital assignment
   call mem_alloc(orbitals_assigned,natoms)
   orbitals_assigned=.false.
   do p=1,nocc_tot
      pdx = occ_orbitals(p)%centralatom
      orbitals_assigned(pdx) = .true.
   end do
   do p=1,nvirt
      pdx = unocc_orbitals(p)%centralatom
      orbitals_assigned(pdx) = .true.
   end do

   ! reorder VOVO integrals from (a,i,b,j) to (a,b,i,j)
   call array4_reorder(VOVO,[1,3,2,4])

   ! print out ccsd fragment and pair interaction energies
   ccsd_mat_tot = array2_init([natoms,natoms])
   ccsd_mat_tmp = array2_init([natoms,natoms])

   call ccsd_energy_full_occ(nocc,nvirt,natoms,ncore,t2_final,t1_final,VOVO,occ_orbitals,&
      & ccsd_mat_tot%val,ccsd_mat_tmp%val)

   call print_ccsd_full_occ(natoms,ccsd_mat_tot%val,orbitals_assigned,mymolecule%distancetable)


   if(DECinfo%InteractionEnergy)then
      call lsquit('InteractionEnergy not implemented in full ccsd_pt',-1)
   endif
   if(DECinfo%PrintInteractionEnergy)then
      call add_dec_interactionenergies(natoms,ccsd_mat_tot%val,orbitals_assigned,&
         & interactionECCSD,MyMolecule%SubSystemIndex,2)
   endif

   ! release ccsd stuff
   call array2_free(ccsd_mat_tot)
   call array2_free(ccsd_mat_tmp)

   ! free integrals
   call array4_free(VOVO)

   if(ccmodel == MODEL_CCSDpT)then
      ! now we calculate fourth-order (which are printed out in print_e4_full) and fifth-order energies
      e4_mat_tot = array2_init([natoms,natoms])
      e4_mat_tmp = array2_init([natoms,natoms])
      e5_mat_tot = array2_init([natoms,natoms])

      call ccsdpt_energy_e4_full(nocc,nvirt,natoms,ncore,t2_final,ccsdpt_t2,occ_orbitals,&
         & e4_mat_tot%val,e4_mat_tmp%val,ccsdpt_e4)

      call ccsdpt_energy_e5_full(nocc,nvirt,natoms,ncore,t1_final,ccsdpt_t1,&
         & occ_orbitals,unocc_orbitals,e5_mat_tot%val,ccsdpt_e5)

      ! print out the fourth- and fifth-order fragment and pair interactin energies
      call print_e4_full(natoms,e4_mat_tot%val,orbitals_assigned,mymolecule%distancetable)

      call print_e5_full(natoms,e5_mat_tot%val,orbitals_assigned,mymolecule%distancetable)


      if(DECinfo%PrintInteractionEnergy)then
         call add_dec_interactionenergies(natoms,e4_mat_tot%val,orbitals_assigned,&
            & InteractionECCSDPT4,MyMolecule%SubSystemIndex)
         write(DECinfo%output,'(A)') ' '
         write(DECinfo%output,'(A)') ' '
         write(DECinfo%output,'(A)') 'Interaction Energies '
         write(DECinfo%output,'(1X,a,g20.10)') 'CCSD Interaction correlation energy                  :',&
            & interactionECCSD
         write(DECinfo%output,'(1X,a,g20.10)') '(fourth order) CCSD(T) Interaction correlation energy:',&
            & InteractionECCSDPT4
         call add_dec_interactionenergies(natoms,e5_mat_tot%val,orbitals_assigned,&
            & InteractionECCSDPT5,MyMolecule%SubSystemIndex)
         write(DECinfo%output,'(1X,a,g20.10)') '(fifth order) CCSD(T) Interaction correlation energy :',&
            & InteractionECCSDPT5
         write(DECinfo%output,'(1X,a,g20.10)') 'Total CCSD(T) Interaction correlation energy         :',&
            & InteractionECCSDPT4 + InteractionECCSDPT5 + InteractionECCSD
         write(DECinfo%output,'(A)') ' '
      endif

      ! release stuff
      call array2_free(e4_mat_tot)
      call array2_free(e4_mat_tmp)
      call array2_free(e5_mat_tot)
   endif



   do i=1,nocc_tot
      call orbital_free(occ_orbitals(i))
   end do
   call mem_dealloc(occ_orbitals)
   do i=1,nvirt
      call orbital_free(unocc_orbitals(i))
   end do
   call mem_dealloc(unocc_orbitals)
   call mem_dealloc(orbitals_assigned)

   ! sum up energies
   ccsdpt_tot = ccsdpt_e4 + ccsdpt_e5

   !MODIFY FOR NEW MODEL
   write(DECinfo%output,*)
   write(DECinfo%output,*)
   write(DECinfo%output,'(1X,a)') '*****************************************************************************'
   if(ccmodel == MODEL_CCSDpT )then
      write(DECinfo%output,'(1X,a)') '*                      Full CCSD(T) calculation is done !                   *'
   else if (ccmodel == MODEL_CCSD ) then
      write(DECinfo%output,'(1X,a)') '*                      Full CCSD calculation is done !                      *'
   else if (ccmodel == MODEL_CC2 ) then
      write(DECinfo%output,'(1X,a)') '*                      Full CC2 calculation is done !                       *'
   else if (ccmodel == MODEL_MP2 ) then
      write(DECinfo%output,'(1X,a)') '*                      Full MP2 calculation is done !                       *'
   else if (ccmodel == MODEL_RPA ) then
      write(DECinfo%output,'(1X,a)') '*                      Full RPA calculation is done !                       *'
   else
      call lsquit("ERROR(ccsolver_justenergy)model not recognized",-1)
   endif
   write(DECinfo%output,'(1X,a)') '*****************************************************************************'
   write(DECinfo%output,*)
   write(DECinfo%output,*)
   write(DECinfo%output,'(1X,a)')   '-------------------------------------------------------------'
   write(DECinfo%output,*)
   write(DECinfo%output,'(1X,a,g20.10)') 'Total CC solver correlation energy           =', ccenergy
   write(DECinfo%output,*)

   ! now update ccenergy with ccsd(t) correction
   ccenergy = ccenergy + ccsdpt_tot

   if(ccmodel == MODEL_CCSDpT)then
      write(DECinfo%output,'(1X,a)')   '-------------------------------------------------------------'
      write(DECinfo%output,*)
      write(DECinfo%output,'(1X,a,g20.10)') 'The E4 doubles and triples contribution =', ccsdpt_e4
      write(DECinfo%output,*)
      write(DECinfo%output,'(1X,a)')   '-------------------------------------------------------------'
      write(DECinfo%output,*)
      write(DECinfo%output,'(1X,a,g20.10)') 'The E5 singles and triples contribution =', ccsdpt_e5
      write(DECinfo%output,*)
      write(DECinfo%output,'(1X,a)')   '-------------------------------------------------------------'
      write(DECinfo%output,*)
      write(DECinfo%output,'(1X,a,g20.10)') 'Total CCSD(T) energy contribution       =', ccsdpt_tot
      write(DECinfo%output,*)
      write(DECinfo%output,'(1X,a)')   '-------------------------------------------------------------'
      write(DECinfo%output,*)
      write(DECinfo%output,'(1X,a,g20.10)') 'Total CCSD(T) correlation energy        =', ccenergy
      write(DECinfo%output,*)
      write(DECinfo%output,'(1X,a)')   '-------------------------------------------------------------'
      write(DECinfo%output,*)
      write(DECinfo%output,*)
   endif

   if( ccmodel /= MODEL_MP2 .and. ccmodel /= MODEL_RPA ) then
      ! free amplitude arrays
      call array2_free(t1_final)
   endif

   call array4_free(t2_final)

   if(ccmodel == MODEL_CCSDpT)then
      call array2_free(ccsdpt_t1)
      call array4_free(ccsdpt_t2)
   endif

end function ccsolver_justenergy

!> \brief For a given fragment, calculate singles and doubles amplitudes and
!> two-electron integrals (a i | bj ) required for CC energy.
!> Intended to be used for CC2 and CCSD (and NOT for MP2).
!> \author Kasper Kristensen
!> \date January 2012
subroutine fragment_ccsolver(MyFragment,t1,t2,VOVO)

   implicit none

   !> Fragment info (only t1 information in MyFragment may be changed here)
   type(decfrag), intent(inout) :: MyFragment
   !> Singles amplitudes t1(a,i)
   type(array2),intent(inout) :: t1
   !> Doubles amplitudes t2(a,i,b,j)
   type(array4),intent(inout) :: t2
   !> Two electron integrals (a i | b j) stored as (a,i,b,j)
   type(array4),intent(inout) :: VOVO

   !INTERNAL PARAMETERS
   type(array4) :: mp2_amp
   integer :: dims(2)
   real(realk) :: ccenergy
   logical :: local

   ! Sanity check: This routine is not intended for MP2
   if(MyFragment%ccmodel == MODEL_MP2) then
      call lsquit('fragment_ccsolver cannot be used for MP2!',&
      & DECinfo%output)
   end if

   local=.true.
#ifdef VAR_MPI
   if(infpar%lg_nodtot>1)local=.false.
#endif

   ! If MyFragment%t1_stored is TRUE, then we reuse the singles amplitudes
   ! from previous fragment calculations to describe long-range
   ! singles effects.
   ! In this case the fragment t1 amplitudes are stored in MyFragment%t1
   if(MyFragment%t1_stored) then
      dims(1) = MyFragment%nunoccAOS
      dims(2) = MyFragment%noccAOS
      t1 = array2_init(dims,MyFragment%t1)
   end if

   if(DECinfo%CCDEBUG)then
      if(DECinfo%use_pnos)then

         !GET MP2 AMPLITUDES TO CONSTRUCT PNOS
         call get_VOVO_integrals( myfragment%mylsitem, myfragment%nbasis, &
            &myfragment%noccAOS, myfragment%nunoccAOS, myfragment%Cv, myfragment%Co, VOVO )
         call mp2_solver( myfragment%noccAOS, myfragment%nunoccAOS, myfragment%ppfock,&
            & myfragment%qqfock, VOVO, mp2_amp )
         call array4_free( VOVO )

         !GET THE CORRELATION DENSITY FOR THE CENTRAL ATOM
         call mem_alloc(MyFragment%occmat,MyFragment%noccAOS,MyFragment%noccAOS)
         call mem_alloc(MyFragment%virtmat,MyFragment%nunoccAOS,MyFragment%nunoccAOS)
         call calculate_corrdens_EOS(mp2_amp,MyFragment) 
         MyFragment%CDset=.true.

         !CALL THE SOLVER WITH PNO ARGUMENT
         call ccsolver_debug(MyFragment%ccmodel,myfragment%Co,myfragment%Cv,&
            & myfragment%fock, myfragment%nbasis,myfragment%noccAOS,&
            & myfragment%nunoccAOS,myfragment%mylsitem,DECinfo%PL,&
            & .true.,myfragment%ppfock,myfragment%qqfock,ccenergy,t1,t2,VOVO,&
            &MyFragment%t1_stored,SOLVE_AMPLITUDES,m2=mp2_amp,use_pnos=DECinfo%use_pnos, fraginfo=myfragment)

         call array4_free(mp2_amp)

      else

         call ccsolver_debug(MyFragment%ccmodel,myfragment%Co,myfragment%Cv,&
            & myfragment%fock, myfragment%nbasis,myfragment%noccAOS,&
            & myfragment%nunoccAOS,myfragment%mylsitem,DECinfo%PL,&
            & .true.,myfragment%ppfock,myfragment%qqfock,ccenergy,t1,t2,VOVO,MyFragment%t1_stored,SOLVE_AMPLITUDES)

      endif
   else
      call ccsolver_par(MyFragment%ccmodel,myfragment%Co,myfragment%Cv,&
         & myfragment%fock, myfragment%nbasis,myfragment%noccAOS,&
         & myfragment%nunoccAOS,myfragment%mylsitem,DECinfo%PL,&
         & .true.,myfragment%ppfock,myfragment%qqfock,ccenergy,t1,t2,VOVO,MyFragment%t1_stored,local)
   endif

   ! Save singles amplitudes in fragment structure
   if(DECinfo%SinglesPolari) then
      call save_fragment_t1_AOSAOSamplitudes(MyFragment,t1)
   end if

end subroutine fragment_ccsolver



!> \brief For a given fragment, calculate combined doubles+singles amplitudes:
!> u(a,i,b,j) = t2(a,i,b,j) + t1(a,i)*t1(b,j)
!> and two-electron integrals (a i | bj ).
!> The EOS indices are extracted for both the occupied and
!> the virtual spaces
!> Intended to be used for CC2 and CCSD (and NOT for MP2).
!> \author Kasper Kristensen
!> \date January 2012
subroutine CCfragment_get_EOS_amplitudes_and_integrals(MyFragment,&
      & uocc,uvirt,VOVOocc,VOVOvirt)

   implicit none

   !> Fragment info (only t1 information in MyFragment may be changed here)
   type(decfrag), intent(inout) :: MyFragment
   !> Combined doubles+singles amplitudes for occupied partitioning
   type(array4),intent(inout) :: uocc
   !> Combined doubles+singles amplitudes for virtual partitioning
   type(array4),intent(inout) :: uvirt
   !> Two electron integrals (a i | b j) stored as (a,i,b,j) for occ part.
   !> (a,b: AOS orbitals;   i,j: EOS orbitals)
   type(array4),intent(inout) :: VOVOocc
   !> Two electron integrals (a i | b j) stored as (a,i,b,j) for virt part.
   !> (a,b: EOS orbitals;   i,j: AOS orbitals)
   type(array4),intent(inout) :: VOVOvirt
   type(array2) :: t1
   type(array4) :: t2,u,VOVO


   ! Solve CC equation to calculate amplitudes and integrals
   ! *******************************************************
   ! Here all output indices in t1,t2, and VOVO are AOS indices.
   call fragment_ccsolver(MyFragment,t1,t2,VOVO)


   ! Extract EOS indices for integrals
   ! *********************************
   call array4_extract_eos_indices_both_schemes(VOVO, &
      & VOVOocc, VOVOvirt, MyFragment)
   call array4_free(VOVO)


   ! Calculate combined single+doubles amplitudes
   ! ********************************************
   ! u(a,i,b,j) = t2(a,i,b,j) + t1(a,i)*t1(b,j)
   call get_combined_SingleDouble_amplitudes(t1,t2,u)
   call array2_free(t1)
   call array4_free(t2)


   ! Extract EOS indices for amplitudes
   ! **********************************
   call array4_extract_eos_indices_both_schemes(u, &
      & uocc, uvirt, MyFragment)
   call array4_free(u)


end subroutine CCfragment_get_EOS_amplitudes_and_integrals



!> \brief Solve MP2 equation:
!> RHS_{bjai} =  - sum_{c} t_{bjci} F_{ca}
!>               - sum_{c} t_{aicj} F_{cb}
!>               + sum_{k} t_{bjak} F_{ki}
!>               + sum_{k} t_{aibk} F_{kj}
!> It is assumed that RHS_{bjai} = R_{aibj} !
!> \author Kasper Kristensen
!> \date February 2011
subroutine mp2_solver(nocc,nvirt,ppfock,qqfock,RHS,t2)

   implicit none
   !> Number of occupied orbitals (dimension of ppfock)
   integer, intent(in) :: nocc
   !> Number of unoccupied orbitals (dimension of qqfock)
   integer, intent(in) :: nvirt
   !> Occupied-occupied block of Fock matrix
   real(realk) :: ppfock(nocc,nocc)
   !> Unoccupied-unoccupied block of Fock matrix
   real(realk) :: qqfock(nvirt,nvirt)
   !> RHS array
   type(array4), intent(in) :: RHS
   !> Solution array
   type(array4), intent(inout) :: t2
   real(realk) :: tcpu1,twall1,tcpu2,twall2

   call LSTIMER('START',tcpu1,twall1,DECinfo%output)

   if(DECinfo%array4OnFile) then ! RHS and t2 values are stored on file
      call mp2_solver_file(nocc,nvirt,ppfock,qqfock,RHS,t2)
   else ! RHS and t2 values are stored in memory
      call mp2_solver_mem(nocc,nvirt,ppfock,qqfock,RHS,t2)
   end if

   call LSTIMER('START',tcpu2,twall2,DECinfo%output)


end subroutine mp2_solver



!> \brief Solve MP2 equation when RHS and t2 are values are stored in memory.
!> See mp2_solver for details about the equation.
!> \author Kasper Kristensen
!> \date February 2011
subroutine mp2_solver_mem(nocc,nvirt,ppfock,qqfock,RHS,t2)

   implicit none
   !> Number of occupied orbitals (dimension of ppfock)
   integer, intent(in) :: nocc
   !> Number of unoccupied orbitals (dimension of qqfock)
   integer, intent(in) :: nvirt
   !> Occupied-occupied block of Fock matrix
   real(realk) :: ppfock(nocc,nocc)
   !> Unoccupied-unoccupied block of Fock matrix
   real(realk) :: qqfock(nvirt,nvirt)
   !> RHS array
   type(array4), intent(in) :: RHS
   !> Solution array
   type(array4), intent(inout) :: t2
   type(array4) :: tmp1,tmp2
   real(realk),pointer :: Cocc_data(:,:), Cvirt_data(:,:), Socc(:,:), Svirt(:,:)
   real(realk),pointer :: EVocc(:), EVvirt(:)
   type(array2) :: Cocc, Cvirt
   integer :: I,J,A,B
   real(realk) :: tcpu, twall, deltaF
   integer :: dims(4), occdims(2), virtdims(2)
   ! real(realk) :: test


   ! Strategy for solving MP2 equation:
   ! 1. Find basis where Fock matrix is diagonal
   ! 2. Transform 2-electron integrals to diagonal basis
   ! 3. In diagonal basis the solution is trivial and the amplitudes are found.
   ! 4. Transform amplitudes in diagonal basis back to LCM basis.

   call LSTIMER('START',tcpu,twall,DECinfo%output)


   write(DECinfo%output,*)
   write(DECinfo%output,*) 'Entering MP2 solver - store array values in memory'
   write(DECinfo%output,*)


   ! Sanity checks
   ! *************
   ! Check that nvirt /= 0.
   if(nvirt<1 .or. nocc<1) then
      write(DECinfo%output,*) 'Number of occupied orbitals = ', nocc
      write(DECinfo%output,*) 'Number of unoccupied orbitals = ', nvirt
      call lsquit('Error in mp2_solver: Number of orbitals is smaller than one!', DECinfo%output)
   endif




   ! Initialize stuff
   ! ****************
   dims = [nvirt,nocc,nvirt,nocc]
   occdims = [nocc,nocc]
   virtdims = [nvirt,nvirt]




   ! 1. Solve Fock eigenvalue problem - each block separately
   ! ********************************************************

   ! OCCUPIED-OCCUPIED BLOCK
   ! '''''''''''''''''''''''

   ! Eigenvectors
   call mem_alloc(Cocc_data,nocc,nocc)

   ! Eigenvalues
   call mem_alloc(EVocc,nocc)

   ! The overlap matrix is simply the unit matrix because
   ! the LCM/MLM basis is orthogonal.
   call mem_alloc(Socc,nocc,nocc)
   Socc=0.0E0_realk
   do i=1,nocc
      Socc(i,i) = 1E0_realk
   end do

   ! Solve eigenvalue problem
   call solve_eigenvalue_problem(nocc,ppfock,Socc,EVocc,Cocc_data)

   ! For later, it is convenient to keep eigenvectors in array2 form
   Cocc = array2_init(occdims,Cocc_data)

   ! Done with some matrices
   call mem_dealloc(Cocc_data)
   call mem_dealloc(Socc)




   ! VIRTUAL-VIRTUAL BLOCK
   ! '''''''''''''''''''''

   ! Eigenvectors
   call mem_alloc(Cvirt_data,nvirt,nvirt)

   ! Eigenvalues
   call mem_alloc(EVvirt,nvirt)

   ! Unit overlap for virtual space
   call mem_alloc(Svirt,nvirt,nvirt)
   Svirt=0.0E0_realk
   do i=1,nvirt
      Svirt(i,i) = 1E0_realk
   end do

   ! Solve eigenvalue problem
   call solve_eigenvalue_problem(nvirt,qqfock,Svirt,EVvirt,Cvirt_data)

   ! For later, it is convenient to keep eigenvectors in array2 form
   Cvirt = array2_init(virtdims,Cvirt_data)

   ! Done with some matrices
   call mem_dealloc(Cvirt_data)
   call mem_dealloc(Svirt)


   call LSTIMER('SOLVE: EIVAL',tcpu,twall,DECinfo%output)



   ! 2. Transform two-electron integrals to diagonal basis
   ! *****************************************************

   ! Using notation that (a,i,b,j) are LCM indices
   ! and (A,I,B,J) are indices in the diagonal basis,
   ! we want to carry out the transformations:
   ! RHS_{AIBJ} = sum_{aibj} C_{aA} C_{iI} C_{bB} C_{jJ} RHS_{aibj} (*)

   ! 1. Init temporary array
   tmp1= array4_init(dims)

   ! 2. Index A: RHS(a,i,b,j) --> tmp1(A,i,b,j)
   call array4_contract1(RHS,Cvirt,tmp1,.true.)

   ! 3. Index I: tmp1(A,i,b,j) --> tmp1(i,A,b,j) --> tmp2(I,A,b,j)
   call array4_reorder(tmp1,[2,1,3,4])
   tmp2= array4_init([nocc,nvirt,nvirt,nocc])

   call array4_contract1(tmp1,Cocc,tmp2,.true.)
   call array4_free(tmp1)

   ! 4. Index J: tmp2(I,A,b,j) --> tmp2(j,b,A,I) --> tmp1(J,b,A,I)
   call array4_reorder(tmp2,[4,3,2,1])
   tmp1= array4_init([nocc,nvirt,nvirt,nocc])
   call array4_contract1(tmp2,Cocc,tmp1,.true.)
   call array4_free(tmp2)

   ! 5. Index B: tmp1(J,b,A,I) --> tmp1(b,J,A,I) --> t2(B,J,A,I)
   call array4_reorder(tmp1,[2,1,3,4])
   t2 = array4_init(dims)
   call array4_contract1(tmp1,Cvirt,t2,.true.)
   call array4_free(tmp1)

   ! Now t2 contains the two-electron integrals (AI|BJ) = (BJ|AI) in the order (B,J,A,I).
   ! [Due to the symmetry it is not necessary to reorder t2 back to (A,I,B,J)].


   call LSTIMER('SOLVE: TRANS 1',tcpu,twall,DECinfo%output)



   ! 3. Solve MP2 equation in the diagonal basis
   ! *******************************************

   ! In the diagonal basis the solution to the MP2 equation is trivial!
   ! The equation:
   ! RHS_{bjai} =  - sum_{c} t_{bjci} F_{ca}
   !               - sum_{c} t_{aicj} F_{cb}
   !               + sum_{k} t_{bjak} F_{ki}
   !               + sum_{k} t_{aibk} F_{kj}
   !
   ! simply becomes (using t_{BJAI} = t_{AIBJ}):
   !
   ! RHS_{BJAI} =  - t_{BJAI} F_{AA}
   !               - t_{AIBJ} F_{BB}
   !               + t_{BJAI} F_{II}
   !               + t_{AIBJ} F_{JJ}
   !            =  t_{BJAI} [ F_{II} + F_{JJ} - F_{AA} - F_{BB} ]
   !
   ! In other words:
   !
   ! t_{BJAI} = RHS_{BJAI} / [ F_{II} + F_{JJ} - F_{AA} - F_{BB} ] (**)
   !

   ! Recalling that currently t_{BJAI} = RHS_{BJAI} and that the
   ! Fock matrix in the diagonal basis are the eigenvalues EVocc and EVvirt -
   ! we simply modify the t2 elements as follows:

   !    test=0E0_realk
   do I=1,nocc
      do A=1,nvirt
         do J=1,nocc
            do B=1,nvirt

               ! deltaF = F_{II} + F_{JJ} - F_{AA} - F_{BB}
               deltaF = EVocc(I) + EVocc(J) - EVvirt(A) - EVvirt(B)

               ! Sanity check
               if( abs(deltaF) < 1e-9 ) then
                  write(DECinfo%output,*) 'WARNING: SMALL NUMBERS OCCURING IN SOLVER!!!'
                  write(DECinfo%output,*) 'WARNING: SOLVER MAY BE UNSTABLE!!!'
                  write(DECinfo%output,*) 'Delta epsilon value = ', deltaF
               end if

               ! Canonical energy check
               ! test = test &
               ! & + (2E0_realk*t2%val(B,J,A,I) - t2%val(B,I,A,J))*t2%val(B,J,A,I)/deltaF

               ! t2 according to (**)
               t2%val(B,J,A,I) = t2%val(B,J,A,I)/deltaF

            end do
         end do
      end do
   end do

   call LSTIMER('SOLVE: CALC T2',tcpu,twall,DECinfo%output)




   ! 4. Transform amplitudes back to LCM/MLM basis
   ! *********************************************

   ! Since LCM/MLM and the diagonal basis are connected by a unitary
   ! transformation this basically corresponds to repeating step 2
   ! above with the transposed transformation matrices:
   !
   ! t_{aibj} = sum_{AIBJ} C_{Aa} C_{Ia} C_{Bb} C_{Jj} t_{AIBJ}

   ! 1. Transpose transformation matrices
   call array2_transpose(Cocc)
   call array2_transpose(Cvirt)

   ! 2. Index A: t(A,I,B,J) --> tmp1(a,I,B,J)
   tmp1= array4_init(dims)
   call array4_contract1(t2,Cvirt,tmp1,.true.)
   call array4_free(t2)

   ! 3. Index I: tmp1(a,I,B,J) --> tmp1(I,a,B,J) --> tmp2(i,a,B,J)
   call array4_reorder(tmp1,[2,1,3,4])
   tmp2= array4_init([nocc,nvirt,nvirt,nocc])
   call array4_contract1(tmp1,Cocc,tmp2,.true.)
   call array4_free(tmp1)

   ! 4. Index J: tmp2(i,a,B,J) --> tmp2(J,B,a,i) --> tmp1(j,B,a,i)
   call array4_reorder(tmp2,[4,3,2,1])
   tmp1 = array4_init([nocc,nvirt,nvirt,nocc])
   call array4_contract1(tmp2,Cocc,tmp1,.true.)
   call array4_free(tmp2)

   ! 5. Index B: tmp1(j,B,a,i) --> tmp1(B,j,a,i) --> t2(b,j,a,i)
   call array4_reorder(tmp1,[2,1,3,4])
   t2 = array4_init(dims)
   call array4_contract1(tmp1,Cvirt,t2,.true.)
   call array4_free(tmp1)

   call LSTIMER('SOLVE: TRANS 2',tcpu,twall,DECinfo%output)


   ! Clean up
   ! ********
   call mem_dealloc(EVocc)
   call mem_dealloc(EVvirt)
   call array2_free(Cocc)
   call array2_free(Cvirt)



end subroutine mp2_solver_mem



!> \brief Solve MP2 equation when RHS and t2 are values are stored on file.
!> See mp2_solver for details about the equation.
!> \author Kasper Kristensen
!> \date February 2011
subroutine mp2_solver_file(nocc,nvirt,ppfock,qqfock,RHS,t2)

   implicit none
   !> Number of occupied orbitals (dimension of ppfock)
   integer, intent(in) :: nocc
   !> Number of unoccupied orbitals (dimension of qqfock)
   integer, intent(in) :: nvirt
   !> Occupied-occupied block of Fock matrix
   real(realk) :: ppfock(nocc,nocc)
   !> Unoccupied-unoccupied block of Fock matrix
   real(realk) :: qqfock(nvirt,nvirt)
   !> RHS array, storing type 1 (see array4_init_file)
   type(array4), intent(in) :: RHS
   !> Solution array, storing type 1 (see array4_init_file)
   type(array4), intent(inout) :: t2
   type(array4) :: tmp1,tmp2
   type(array4) :: RHSaib,RHSbaj,RHStmp,t2tmp
   real(realk),pointer :: Cocc_data(:,:), Cvirt_data(:,:), Socc(:,:), Svirt(:,:)
   real(realk),pointer :: EVocc(:), EVvirt(:)
   type(array2) :: Cocc, Cvirt,CoccT,CvirtT
   integer :: I,J,A,B
   real(realk) :: tcpu, twall, deltaF
   integer :: occdims(2), virtdims(2)



   ! Strategy for solving MP2 equation:
   ! 1. Find basis where Fock matrix is diagonal
   ! 2. Transform 2-electron integrals to diagonal basis
   ! 3. In diagonal basis the solution is trivial and the amplitudes are found.
   ! 4. Transform amplitudes in diagonal basis back to LCM basis.

   ! Steps 2-4 necessarily overlap when we store array values on file.
   call LSTIMER('START',tcpu,twall,DECinfo%output)


   write(DECinfo%output,*)
   write(DECinfo%output,*) 'Entering MP2 solver - store array values on file'
   write(DECinfo%output,*)


   ! Sanity checks
   ! *************
   ! Check that nvirt /= 0.
   if(nvirt<1 .or. nocc<1) then
      write(DECinfo%output,*) 'Number of occupied orbitals = ', nocc
      write(DECinfo%output,*) 'Number of unoccupied orbitals = ', nvirt
      call lsquit('Error in mp2_solver: Number of orbitals is smaller than one!', DECinfo%output)
   endif


   ! Initialize stuff
   ! ****************
   occdims = [nocc,nocc]
   virtdims = [nvirt,nvirt]



   ! 1. Solve Fock eigenvalue problem - each block separately
   ! ********************************************************

   ! OCCUPIED-OCCUPIED BLOCK
   ! '''''''''''''''''''''''

   ! Eigenvectors
   call mem_alloc(Cocc_data,nocc,nocc)
   ! Eigenvalues
   call mem_alloc(EVocc,nocc)

   ! The overlap matrix is simply the unit matrix because
   ! the LCM/MLM basis is orthogonal.
   call mem_alloc(Socc,nocc,nocc)
   Socc=0.0E0_realk
   do i=1,nocc
      Socc(i,i) = 1E0_realk
   end do

   ! Solve eigenvalue problem
   call solve_eigenvalue_problem(nocc,ppfock,Socc,EVocc,Cocc_data)

   ! For later, it is convenient to keep eigenvectors in array2 form
   ! and also to have the transposed matrices available
   Cocc = array2_init(occdims,Cocc_data)
   CoccT = array2_init(occdims)
   call array2_copy(CoccT,Cocc)
   call array2_transpose(CoccT)

   ! Done with some matrices
   call mem_dealloc(Cocc_data)
   call mem_dealloc(Socc)




   ! VIRTUAL-VIRTUAL BLOCK
   ! '''''''''''''''''''''

   ! Eigenvectors
   call mem_alloc(Cvirt_data,nvirt,nvirt)

   ! Eigenvalues
   call mem_alloc(EVvirt,nvirt)

   ! Unit overlap for virtual space
   call mem_alloc(Svirt,nvirt,nvirt)
   Svirt=0.0E0_realk
   do i=1,nvirt
      Svirt(i,i) = 1E0_realk
   end do

   ! Solve eigenvalue problem
   call solve_eigenvalue_problem(nvirt,qqfock,Svirt,EVvirt,Cvirt_data)

   ! For later, it is convenient to keep eigenvectors in array2 form
   ! and also to have the transposed matrices available
   Cvirt = array2_init(virtdims,Cvirt_data)
   CvirtT = array2_init(virtdims)
   call array2_copy(CvirtT,Cvirt)
   call array2_transpose(CvirtT)


   ! Done with some matrices
   call mem_dealloc(Cvirt_data)
   call mem_dealloc(Svirt)
   call LSTIMER('SOLVE: EIVAL',tcpu,twall,DECinfo%output)



   ! Transform three indices to diagonal basis (aib-->AIB)
   ! *****************************************************

   ! Using the notation that (a,i,b,j) are LCM indices
   ! and (A,I,B,J) are indices in the diagonal basis,
   ! we want to carry out the transformations:
   ! RHS_{AIBJ} = sum_{aibj} C_{aA} C_{iI} C_{bB} C_{jJ} RHS_{aibj}

   ! Since we cannot have full four-dimensional arrays in memory,
   ! we do this in steps. First, for a fixed "j" we transform the other indices:
   !
   ! RHS_{AIBj} = sum_{aib} C_{aA} C_{iI} C_{bB} RHS_{aibj}

   ! Temporary RHS where three indices are transformed, stored on file.
   RHStmp = array4_init_file([nvirt,nvirt,nocc,nocc],2,.false.)
   ! Temporary array for keeping RHS_{aibj} for a fixed j (last dimension is one).
   RHSaib = array4_init([nvirt,nocc,nvirt,1])
   call array4_open_file(RHS)
   call array4_open_file(RHStmp)


   do j=1,nocc

      ! 1. Read in RHS_{aibj} for fixed j
      call array4_read_file_type1(RHS,j,&
         & RHSaib%val(1:nvirt,1:nocc,1:nvirt,1),nvirt,nocc,nvirt)

      ! 2. Index A: RHS(a,i,b,j) --> tmp1(A,i,b,j)
      tmp1= array4_init([nvirt,nocc,nvirt,1])
      call array4_contract1(RHSaib,Cvirt,tmp1,.true.)

      ! 3. Index I: tmp1(A,i,b,j) --> tmp1(i,A,b,j) --> tmp2(I,A,b,j)
      call array4_reorder(tmp1,[2,1,3,4])
      tmp2= array4_init([nocc,nvirt,nvirt,1])
      call array4_contract1(tmp1,Cocc,tmp2,.true.)
      call array4_free(tmp1)

      ! 4. Index B: tmp2(I,A,b,j) --> tmp2(b,A,I,j) --> tmp1(B,A,I,j)
      call array4_reorder(tmp2,[3,2,1,4])
      tmp1= array4_init([nvirt,nvirt,nocc,1])
      call array4_contract1(tmp2,Cvirt,tmp1,.true.)
      call array4_free(tmp2)

      ! 5. Write to file referenced by temporary RHS array (storing type 2)
      do I=1,nocc
         call array4_write_file_type2(RHStmp,I,j,tmp1%val(1:nvirt,1:nvirt,I,1),nvirt,nvirt)
      end do
      call array4_free(tmp1)

   end do
   call array4_close_file(RHS,'KEEP')
   call array4_free(RHSaib)
   call LSTIMER('SOLVE: STEP 1',tcpu,twall,DECinfo%output)
   ! Now the file referenced by RHStmp contains RHS_{AIBj},
   ! stored in the order (B,A,I,j) using storing type 2.



   ! Transform j-->J; Solve equation in diag basis; Back transform ABJ-->abj
   ! ***********************************************************************

   ! Temporary array
   t2tmp = array4_init_file([nvirt,nvirt,nocc,nocc],2,.false.)
   call array4_open_file(t2tmp)


   I_loop: do I=1,nocc


      ! Read in RHS_{AIBj} for fixed I
      ! ------------------------------
      RHSbaj=array4_init([nvirt,nvirt,nocc,1])
      do j=1,nocc
         call array4_read_file_type2(RHStmp,I,j,&
            & RHSbaj%val(1:nvirt,1:nvirt,j,1), nvirt,nvirt)
      end do
      ! Now RHSbaj contains elements RHS_{AIBj} for a fixed I
      ! stored in the order (B,A,j,I)


      ! Transform final index j-->J
      ! ---------------------------
      ! Reorder: RHSbaj(B,A,j,I) --> RHSbaj(j,B,A,I)
      call array4_reorder(RHSbaj,[3,1,2,4])


      ! Transform: RHSbaj(j,B,A,I) --> tmp1(J,B,A,I)
      tmp1=array4_init([nocc,nvirt,nvirt,1])
      call array4_contract1(RHSbaj,Cocc,tmp1,.true.)
      call array4_free(RHSbaj)


      ! Now tmp1 contains the RHS_{BJAI} in the diagonal basis
      ! stored in the order (J,B,A,I) [fixed I].


      ! Divide by deltaF (solve equation in diagonal basis)
      ! ---------------------------------------------------
      ! In the diagonal basis the t2 solution vector is simply:
      ! t2_{BJAI} = RHS_{BJAI} / [ F_{II} + F_{JJ} - F_{AA} - F_{BB} ] (**)
      ! [See (**) in subroutine mp2_solver_mem]
      do A=1,nvirt
         do B=1,nvirt
            do J=1,nocc

               ! deltaF = F_{II} + F_{JJ} - F_{AA} - F_{BB}
               deltaF = EVocc(I) + EVocc(J) - EVvirt(A) - EVvirt(B)

               ! Sanity check
               if( abs(deltaF) < 1e-9 ) then
                  write(DECinfo%output,*) 'WARNING: SMALL NUMBERS OCCURING IN SOLVER!!!'
                  write(DECinfo%output,*) 'WARNING: SOLVER MAY BE UNSTABLE!!!'
                  write(DECinfo%output,*) 'Delta epsilon value = ', deltaF
               end if

               ! Amplitude element according to (**)
               tmp1%val(J,B,A,1) = tmp1%val(J,B,A,1)/deltaF

            end do
         end do
      end do

      ! Now tmp1 contains the t2_{BJAI} solution vector in the diagonal basis
      ! stored in the order (J,B,A,I) [fixed I], and we "just"
      ! need to transform back to the original basis.


      ! Transform back three indices: JBA --> jba
      ! -----------------------------------------
      ! Note: Use Transposed transformation matrices to back transform

      ! 1. Index j: tmp1(J,B,A,I) -->  tmp2(j,B,A,I)
      tmp2= array4_init([nocc,nvirt,nvirt,1])
      call array4_contract1(tmp1,CoccT,tmp2,.true.)
      call array4_free(tmp1)

      ! 2. Index b: tmp2(j,B,A,I) --> tmp2(B,A,j,I) --> tmp1(b,A,j,I)
      call array4_reorder(tmp2,[2,3,1,4])
      tmp1= array4_init([nvirt,nvirt,nocc,1])
      call array4_contract1(tmp2,CvirtT,tmp1,.true.)
      call array4_free(tmp2)

      ! 3. Index a: tmp1(b,A,j,I) --> tmp1(A,b,j,I) --> tmp2(a,b,j,I)
      call array4_reorder(tmp1,[2,1,3,4])
      tmp2= array4_init([nvirt,nvirt,nocc,1])
      call array4_contract1(tmp1,CvirtT,tmp2,.true.)
      call array4_free(tmp1)


      ! Write solution vector t2_{aIbj} to file using storing type 2, order: (a,b,j,I)
      ! ------------------------------------------------------------------------------

      ! Note: Now we only need to transform the last "I" index of t2 back.
      ! To do this we need to save on file such that we later can read in the full set of "I's"
      do j=1,nocc
         call array4_write_file_type2(t2tmp,j,I,&
            & tmp2%val(1:nvirt,1:nvirt,j,1), nvirt,nvirt )
      end do
      call array4_free(tmp2)

   end do I_loop

   call array4_close_file(RHStmp,'DELETE')
   call array4_free(RHStmp)
   call LSTIMER('SOLVE: STEP 2',tcpu,twall,DECinfo%output)
   ! Now the file assosicated with t2tmp contains the final amplitudes,
   ! except that the I index must be transformed back.



   ! Transform final t2 index (I-->i) and write solution vector to file
   ! ******************************************************************

   ! Final t2 solution vector
   t2 = array4_init_file([nvirt,nocc,nvirt,nocc],1,.false.)
   call array4_open_file(t2)

   do j=1,nocc

      tmp1 = array4_init([nvirt,nvirt,nocc,1])
      do I=1,nocc
         call array4_read_file_type2(t2tmp,j,I,tmp1%val(1:nvirt,1:nvirt,I,1),nvirt,nvirt)
      end do
      ! tmp1 now contains amplitudes t2_{aIbj} in the order (a,b,I,j) for fixed j.


      ! Transform final index I-->i
      ! ---------------------------

      ! tmp1(a,b,I,j) --> tmp1(I,a,b,j) --> tmp2(i,a,b,j)
      call array4_reorder(tmp1,[3,1,2,4])
      tmp2 = array4_init([nocc,nvirt,nvirt,1])
      call array4_contract1(tmp1,CoccT,tmp2,.true.)
      call array4_free(tmp1)
      ! Reorder to final storing order: tmp2(i,a,b,j) --> tmp2(a,i,b,j)
      call array4_reorder(tmp2,[2,1,3,4])


      ! Write t2_{aibj} to file for each j
      ! ----------------------------------
      call array4_write_file_type1(t2,j,&
         & tmp2%val(1:nvirt,1:nocc,1:nvirt,1),nvirt,nocc,nvirt)
      call array4_free(tmp2)

   end do

   call array4_close_file(t2tmp,'DELETE')
   call array4_free(t2tmp)
   call array4_close_file(t2,'KEEP')
   call LSTIMER('SOLVE: STEP 3',tcpu,twall,DECinfo%output)


   ! Clean up
   ! ********
   call mem_dealloc(EVocc)
   call mem_dealloc(EVvirt)
   call array2_free(Cocc)
   call array2_free(Cvirt)
   call array2_free(CoccT)
   call array2_free(CvirtT)



end subroutine mp2_solver_file


!> \brief adaption of the ccsolver routine, rebuild for the 
! use of parallel distributed memory. This solver is acutally a bit
! complicated in structure if used in an MPI framework. Most of the work
! happens hidden in the type(array) structure. It is highly recommended to
! begin implementing new features with setting local=.true. at the beginning
! and running without .spawn_comm_procs in the **CC input section. On
! INPUT:
! Co_f,Cv_f : the occupied and virtual orbital transformation coefficients
! fock_f      : the ao fock matrix
! nb,no,nv    : number of atomic, occupied and virtual orbitals respectively
! mylsitem    : the typical lsitem structure
! ccPrintLevel: print level (might be removed due to DECinfo%PL)
! fragment_job: specify whether it is a fragment job or a full calc
! ppfock_f    : occ occ fock matrix
! qqfock_f    : virt virt fock matrix
! longrange_singles : the longrange singles correction for the fock matrix on input
! local       : boolean which steers whether everything should be treated locally
!
! OUTPUT:
! ccenergy    : output correlation energy
! t1_final    : final singles amplitudes, will be allocated in the solver, output
! t2_final    : final doubles amplitudes, will be allocated in the solver, output
! VOVO        : mo-integral WITHOUT T1 trafo on output
!> \author Patrick Ettenhuber (heavily adapted version from Marcin)
subroutine ccsolver_par(ccmodel,Co_f,Cv_f,fock_f,nb,no,nv, &
      & mylsitem,ccPrintLevel,fragment_job,ppfock_f,qqfock_f,ccenergy, &
      & t1_final,t2_final,VOVO,longrange_singles,local)

   implicit none

   !> CC model
   integer,intent(inout) :: ccmodel
   !> Number of occupied orbitals in full molecule/fragment AOS
   integer, intent(in)                       :: no
   !> Number of virtual orbitals in full molecule/fragment AOS
   integer, intent(in)                       :: nv
   !> Number of basis functions in full molecule/atomic extent
   integer, intent(in)                       :: nb
   !> Fock matrix in AO basis for fragment or full molecule
   real(realk), dimension(nb,nb), intent(in) :: fock_f
   !> Occupied MO coefficients for fragment/full molecule
   real(realk), dimension(nb,no), intent(in) :: Co_f
   !> Virtual MO coefficients for fragment/full molecule
   real(realk), dimension(nb,nv), intent(in) :: Cv_f
   !> Occ-occ block of Fock matrix in MO basis
   real(realk), dimension(no,no), intent(in) :: ppfock_f
   !> Virt-virt block of Fock matrix in MO basis
   real(realk), dimension(nv,nv), intent(in) :: qqfock_f
   real(realk),pointer                       :: dens(:,:)
   !> Is this a fragment job (true) or a full molecular calculation (false)
   logical, intent(in)                       :: fragment_job
   !> LS item information
   type(lsitem), intent(inout)               :: mylsitem
   !> How much to print? ( ccPrintLevel>0 --> print info stuff)
   integer, intent(in)                       :: ccPrintLevel
   !> Coupled cluster energy for fragment/full molecule
   real(realk),intent(inout)                 :: ccenergy
   !> Final singles amplitudes
   type(array2),intent(inout)                :: t1_final
   !> Final doubles amplitudes
   type(array4),intent(inout)                :: t2_final
   !> Two electron integrals (a i | b j) stored as (a,i,b,j)
   type(array4),intent(inout)                :: VOVO
   !> Include long-range singles effects using singles amplitudes
   !> from previous fragment calculations.
   !> IMPORTANT: If this it TRUE, then the singles amplitudes for the fragment
   !> (from previous calculations) must be stored in t1_final at input!
   logical,intent(in)                        :: longrange_singles
   logical,intent(in)                        :: local
   !
   !> Do an MO-based CCSD calculation?
   logical :: mo_ccsd
   !> full set of MO integrals (non-T1-transformed)
   type(array) :: pgmo_diag, pgmo_up
   type(MObatchInfo) :: MOinfo
   !
   !work stuff
   real(realk),pointer :: Co_d(:,:),Cv_d(:,:),Co2_d(:,:), Cv2_d(:,:),focc(:),fvirt(:)
   real(realk),pointer :: ppfock_d(:,:),qqfock_d(:,:),Uocc(:,:),Uvirt(:,:)
   integer, dimension(2) :: occ_dims, virt_dims, ao2_dims, ampl2_dims
   integer, dimension(4) :: ampl4_dims
   type(array)  :: fock,Co,Cv,Co2,Cv2
   type(array)  :: ppfock,qqfock,pqfock,qpfock
   type(array)  :: ifock,delta_fock
   type(array)  :: iajb
   type(array), pointer :: t2(:),omega2(:)
   type(array), pointer :: t1(:),omega1(:)
   type(array) :: omega1_opt, t1_opt, omega1_prec
   type(array) :: omega2_opt, t2_opt, omega2_prec, u
   type(array) :: xo,yo,xv,yv,h1
   type(array) :: Lmo
   real(realk)             :: test_norm,two_norm_total, one_norm_total,one_norm1, one_norm2, prev_norm
   real(realk), pointer    :: B(:,:),c(:)
   integer                 :: iter,last_iter,i,j,k,l
   logical                 :: crop_ok,break_iterations,saferun
   type(ri)                :: l_ao
   type(array)             :: ppfock_prec, qqfock_prec, qpfock_prec
   real(realk)             :: tcpu, twall, ttotend_cpu, ttotend_wall, ttotstart_cpu, ttotstart_wall
   real(realk)             :: iter_cpu,iter_wall
   integer                 :: nnodes
   character(3), parameter :: safefilet11 = 't11'
   character(3), parameter :: safefilet12 = 't12'
   character(3), parameter :: safefilet1f = 't1f'
   character(3), parameter :: safefilet21 = 't21'
   character(3), parameter :: safefilet22 = 't22'
   character(3), parameter :: safefilet2f = 't2f'
   real(realk) :: time_work, time_comm, time_idle, time_fock_mat, time_prec1
   real(realk) :: time_start_guess,time_copy_opt,time_crop_mat,time_energy,time_iter 
   real(realk) :: time_main,time_mixing,time_mo_ints,time_new_guess,time_norm,time_residual
   real(realk) :: time_solve_crop,time_t1_trafo,time_write,time_finalize
   !SOME DUMMIES FOR TESTING
   type(array4) :: iajb_a4
   type(array)            :: tmp
   character(ARR_MSG_LEN) :: msg
   integer                :: ii, jj, aa, bb, old_iter
   logical                :: restart, w_cp, restart_from_converged, collective
   character(4)           :: atype

   time_work        = 0.0E0_realk
   time_comm        = 0.0E0_realk
   time_idle        = 0.0E0_realk
   time_fock_mat    = 0.0E0_realk
   time_prec1       = 0.0E0_realk
   time_start_guess = 0.0E0_realk
   time_copy_opt    = 0.0E0_realk
   time_crop_mat    = 0.0E0_realk
   time_energy      = 0.0E0_realk
   time_iter        = 0.0E0_realk
   time_main        = 0.0E0_realk
   time_mixing      = 0.0E0_realk
   time_mo_ints     = 0.0E0_realk
   time_new_guess   = 0.0E0_realk
   time_norm        = 0.0E0_realk
   time_residual    = 0.0E0_realk
   time_solve_crop  = 0.0E0_realk
   time_t1_trafo    = 0.0E0_realk
   time_write       = 0.0E0_realk
   time_finalize    = 0.0E0_realk

   collective       = .true.

   call time_start_phase(PHASE_WORK, twall = ttotstart_wall, tcpu = ttotstart_cpu )

   !Set defaults
   restart     = .false.
   w_cp        = .false.
   saferun     = (.not.DECinfo%CCSDnosaferun)


   nnodes      = 1

#ifdef VAR_MPI
   nnodes      = infpar%lg_nodtot
   w_cp        = DECinfo%spawn_comm_proc

   call time_start_phase(PHASE_COMM, at = time_work)

   if ( w_cp ) call lspdm_start_up_comm_procs

   call time_start_phase(PHASE_WORK, at = time_comm)

#ifndef COMPILER_UNDERSTANDS_FORTRAN_2003
   call lsquit("ERROR(ccsolver_par):Your compiler does not support certain&
      & features needed to run that part of the code. Use a compiler supporting&
      & Fortran 2003 features",-1)
#endif

#endif


   ! Sanity check 1: Number of orbitals
   if( (nv < 1) .or. (no < 1) ) then
      write(DECinfo%output,*) 'Number of occupied orbitals = ', no
      write(DECinfo%output,*) 'Number of virtual  orbitals = ', nv
      call lsquit('ccsolver: Empty occupied or virtual space!',DECinfo%output)
   endif

   ! Sanity check 2: Singles amplitudes initiated appropriately
   if(longrange_singles) then
      if(.not. associated(t1_final%val)) then
         call lsquit('ccsolver: Long range singles corrections requested, &
            & but t1_final does not contain existing amplitudes!',DECinfo%output)
      end if
   end if

   if(ccmodel == MODEL_MP2)then
      atype = 'REAR'
   else
      atype = 'LDAR'
   endif

   ! go to a (pseudo) canonical basis
   call mem_alloc( focc,     no     )
   call mem_alloc( fvirt,    nv     )
   call mem_alloc( Co_d,     nb, no )
   call mem_alloc( Cv_d,     nb, nv )
   call mem_alloc( Co2_d,    nb, no )
   call mem_alloc( Cv2_d,    nb, nv )
   call mem_alloc( ppfock_d, no, no )
   call mem_alloc( qqfock_d, nv, nv )
   call mem_alloc( Uocc,     no, no )
   call mem_alloc( Uvirt,    nv, nv )

   if(DECinfo%CCSDpreventcanonical)then
      !no diagonalization
      Co_d   = Co_f
      Cv_d   = Cv_f
      ppfock_d = ppfock_f
      qqfock_d = qqfock_f
      Uocc=0.0E0_realk
      Uvirt=0.0E0_realk
      do ii=1,no
         Uocc(ii,ii) = 1.0E0_realk
      enddo
      do aa=1,nv
         Uvirt(aa,aa) = 1.0E0_realk
      enddo
   else
      call get_canonical_integral_transformation_matrices(no,nv,nb,ppfock_f,qqfock_f,Co_f,Cv_f,&
         & Co_d,Cv_d,Uocc,Uvirt,focc,fvirt)
      ppfock_d = 0.0E0_realk
      qqfock_d = 0.0E0_realk
      do ii=1,no
         ppfock_d(ii,ii) = focc(ii)
      enddo
      do aa=1,nv
         qqfock_d(aa,aa) = fvirt(aa)
      enddo
   endif

   call mem_dealloc( focc  )
   call mem_dealloc( fvirt )

   ! Copy MO coeffcients. It is very convenient to store them twice to handle transformation
   ! (including transposed MO matrices) efficiently. 
   Co2_d = Co_d
   Cv2_d = Cv_d


   ! dimension vectors
   occ_dims   = [nb,no]
   virt_dims  = [nb,nv]
   ao2_dims   = [nb,nb]
   ampl4_dims = [nv,nv,no,no]
   ampl2_dims = [nv,no]

   ! create transformation matrices in array form
   Co   = array_minit( occ_dims, 2, local=local, atype=atype )
   Cv   = array_minit( virt_dims,2, local=local, atype=atype )
   Co2  = array_minit( occ_dims, 2, local=local, atype=atype )
   Cv2  = array_minit( virt_dims,2, local=local, atype=atype )
   fock = array_minit( ao2_dims, 2, local=local, atype=atype )

   call array_convert( Co_d,   Co   )
   call array_convert( Cv_d,   Cv   )
   call array_convert( Co2_d,  Co2  )
   call array_convert( Cv2_d,  Cv2  )
   call array_convert( fock_f, fock )

   call mem_dealloc( Co_d )
   call mem_dealloc( Cv_d )
   call mem_dealloc( Co2_d )
   call mem_dealloc( Cv2_d )

   ! Get Fock matrix correction (for fragment and/or frozen core)
   ! ************************************************************
   ! Full molecule/frozen core: The correction corresponds to difference between actual Fock matrix
   !                            and Fock matrix where the density is made from only valence orbitals.
   ! Fragment: The correction correspond to the difference between actual Fock matrix
   !           and Fock matrix calculated from a "fragment density" determined from
   !           fragment's occupied molecular orbitals (which for frozen core includes only valence
   !           orbitals).
   call time_start_phase(PHASE_work, at = time_work, twall = time_fock_mat)

   ! Density corresponding to input MOs
   call mem_alloc(dens,nb,nb)
   call get_density_from_occ_orbitals(nb,no,Co%elm2,dens)

   ifock=array_minit( ao2_dims, 2, local=local, atype='LDAR' )

   if(fragment_job) then ! fragment: calculate correction


      if(longrange_singles) then
         ! Get Fock matrix using singles amplitudes from previous
         ! fragment calculation, thereby effectively including long-range
         ! correlated polarization effects
         call Get_AOt1Fock(mylsitem,t1_final,ifock,no,nv,nb,Co,Co2,Cv2)
      else
         ! Fock matrix for fragment from density made from input MOs
         call get_fock_matrix_for_dec(nb,dens,mylsitem,ifock,.true.)
         !print *,"DEBUGGGING: zero fragment iFOck instead of calculating it"
         !call array_zero(ifock)
      end if

      ! Long range Fock correction:
      !delta_fock = getFockCorrection(fock,ifock)
      delta_fock = array_minit( ao2_dims, 2, local=local, atype='LDAR' )

      call array_cp_data(fock,delta_fock)
      call array_add(delta_fock,-1.0E0_realk,ifock)

   else 
      ! Full molecule: deltaF = F(Dcore) for frozen core (0 otherwise)
      if(DECinfo%frozencore) then
         ! Fock matrix from input MOs
         call get_fock_matrix_for_dec(nb,dens,mylsitem,ifock,.true.)
         !print *,"DEBUGGGING: zero iFOck instead of calculating it"
         !call array_zero(ifock)
         ! Correction to actual Fock matrix
         delta_fock=array_minit( ao2_dims, 2, local=local, atype='LDAR' )
         call array_cp_data(fock,delta_fock)
         call array_add(delta_fock,-1.0E0_realk,ifock)
      else
         delta_fock=array_minit( ao2_dims,2,local=local, atype='LDAR' )
         call array_zero(delta_fock)
      end if
   end if



   if(DECinfo%PL>1)call time_start_phase(PHASE_work, at = time_work, ttot = time_fock_mat, &
      &twall = time_prec1 , labelttot = 'CCSOL: AO FOCK MATRIX :', output = DECinfo%output)



   ! get fock matrices, used in Preconditioning and MP2

   ppfock_prec = array_minit( [no,no], 2, local=local, atype='REPD' )
   qqfock_prec = array_minit( [nv,nv], 2, local=local, atype='REPD' )
   qpfock_prec = array_minit( [nv,no], 2, local=local, atype='REPD' )

   call array_change_atype_to_rep( ppfock_prec, local )
   call array_change_atype_to_rep( qqfock_prec, local )

   if(DECinfo%precondition_with_full) then
      call array_convert( ppfock_d, ppfock_prec )
      call array_convert( qqfock_d, qqfock_prec )
      call array_zero(qpfock_prec)
   else
      tmp = array_minit( [nb,no], 2, local=local, atype='LDAR' )
      call array_contract_outer_indices_rl(1.0E0_realk,fock,Co2,0.0E0_realk,tmp)
      call array_contract_outer_indices_ll(1.0E0_realk,Co,tmp,0.0E0_realk,ppfock_prec)
      call array_free(tmp)

      tmp = array_minit( [nb,nv], 2, local=local, atype='LDAR'  )
      call array_contract_outer_indices_rl(1.0E0_realk,fock,Cv2,0.0E0_realk,tmp)
      call array_contract_outer_indices_ll(1.0E0_realk,Cv,tmp,0.0E0_realk,qqfock_prec)
      call array_free(tmp)

      tmp = array_minit( [nb,no], 2, local=local, atype='LDAR'  )
      call array_contract_outer_indices_rl(1.0E0_realk,fock,Co2,0.0E0_realk,tmp)
      call array_contract_outer_indices_ll(1.0E0_realk,Cv,tmp,0.0E0_realk,qpfock_prec)
      call array_free(tmp)
   end if

   if( ccmodel /= MODEL_MP2 .and. ccmodel /= MODEL_RPA )then
      call array_change_atype_to_d( ppfock_prec )
      call array_change_atype_to_d( qqfock_prec )
   endif


   call array_free(ifock)
   call mem_dealloc(dens)


   if(DECinfo%PL>1)call time_start_phase(PHASE_work, at = time_work, ttot = time_prec1,&
      &labelttot = 'CCSOL: PRECOND. INIT. :', output = DECinfo%output)

   call mem_dealloc( ppfock_d )
   call mem_dealloc( qqfock_d )

   ! allocate things
   if(DECinfo%use_singles) then
      call mem_alloc( t1,     DECinfo%ccMaxIter )
      call mem_alloc( omega1, DECinfo%ccMaxIter )
      ppfock=array_minit( [no,no], 2, local=local, atype='LDAR' )
      pqfock=array_minit( [no,nv], 2, local=local, atype='LDAR' )
      qpfock=array_minit( [nv,no], 2, local=local, atype='LDAR' )
      qqfock=array_minit( [nv,nv], 2, local=local, atype='LDAR' )
   end if
   call mem_alloc(t2,DECinfo%ccMaxIter)
   call mem_alloc(omega2,DECinfo%ccMaxIter)

   ! initialize T1 matrices and fock transformed matrices for CC pp,pq,qp,qq
   if(CCmodel /= MODEL_MP2 .and. ccmodel /= MODEL_RPA) then
      xo = array_minit( occ_dims, 2, local=local, atype='LDAR' )
      yo = array_minit( occ_dims, 2, local=local, atype='LDAR' )
      xv = array_minit( virt_dims,2, local=local, atype='LDAR' )
      yv = array_minit( virt_dims,2, local=local, atype='LDAR' )
   end if

   if(CCmodel == MODEL_RPA) then
     iajb=array_minit( [no,nv,no,nv], 4, local=.false., atype='TDAR' )
   else
     iajb=array_minit( [no,nv,no,nv], 4, local=local, atype='TDAR' )
   endif
   call array_zero(iajb)

  ! if(ccmodel == MODEL_RPA) call lsquit("ERROR(ccsolver_par): model not implemented: RPA",-1)


   call mem_alloc( B, DECinfo%ccMaxIter, DECinfo%ccMaxIter )
   call mem_alloc( c, DECinfo%ccMaxIter                    )

   call time_start_phase( PHASE_work, at = time_work, twall = time_start_guess )




   ! get guess amplitude vectors in the first iteration --> zero if no
   ! restart, else the t*.restart files are read
   two_norm_total = DECinfo%ccConvergenceThreshold + 1.0E0_realk
   if(DECinfo%use_singles)then

      t1(1) = array_minit( ampl2_dims, 2, local=local, atype='REPD' )
      t2(1) = array_minit( ampl4_dims, 4, local=local, atype='TDAR' )

      call get_guess_vectors(restart,old_iter,two_norm_total,ccenergy,t2(1),iajb,Co,Cv,&
         & ppfock_prec,qqfock_prec,qpfock_prec, mylsitem, local, safefilet21,safefilet22, safefilet2f, &
         & t1(1),safefilet11,safefilet12, safefilet1f  )
   else

      !if MP2, just zero the array, and keep it in PDM all the time
         atype = 'TDAR'
         t2(1) = array_minit( ampl4_dims, 4, local=local, atype=atype )
      if(ccmodel == MODEL_MP2 )then
         old_iter = 0
      else
         call get_guess_vectors(restart,old_iter,two_norm_total,ccenergy,t2(1),iajb,Co,Cv,&
            & ppfock_prec,qqfock_prec,qpfock_prec,mylsitem,local,safefilet21,safefilet22, safefilet2f )
      endif
   endif
   restart_from_converged = (two_norm_total < DECinfo%ccConvergenceThreshold)

   call array_free(qpfock_prec)


   if(DECinfo%PL>1)call time_start_phase( PHASE_WORK, at = time_work, ttot = time_start_guess,&
      &labelttot = 'CCSOL: STARTING GUESS :', output = DECinfo%output, twall = time_main  )



   ! title
   Call print_ccjob_header(ccmodel,ccPrintLevel,fragment_job,&
      &.false.,nb,no,nv,DECinfo%ccMaxDIIS,restart,restart_from_converged,old_iter)


   If_not_converged: if(.not.restart_from_converged)then

      mo_ccsd = .false.
      if (DECinfo%MOCCSD) mo_ccsd = .true.

      INTEGRAL : if(ccmodel == MODEL_MP2) then

         call get_mo_integral_par( iajb, Co, Cv, Co, Cv, mylsitem, local, collective )
         call get_mp2_starting_guess( iajb, t2(1), ppfock_prec, qqfock_prec, local )

      else

         !============================================================================!
         !                          MO-CCSD initialization                            !
         !____________________________________________________________________________!
         ! Check if there is enough memory to performed an MO-CCSD calculation.
         !   YES: get full set of t1 free gmo and pack them
         !   NO:  returns mo_ccsd == .false. and switch to standard CCSD.
         if (mo_ccsd) then
            if(DECinfo%PL>1)call time_start_phase( PHASE_work, at = time_work, twall = time_mo_ints ) 

            call get_t1_free_gmo(mo_ccsd,mylsitem,Co%elm2,Cv2%elm2,iajb,pgmo_diag,pgmo_up, &
               & nb,no,nv,CCmodel,MOinfo)

            if(DECinfo%PL>1)call time_start_phase( PHASE_work, at = time_work, ttot = time_mo_ints,&
               &labelttot = 'CCSOL: INIT MO INTS   :', output = DECinfo%output )

         end if

         !JOHANNES MODEL_RPA
         if (ccmodel == MODEL_RPA) then
            if(DECinfo%PL>1)call time_start_phase( PHASE_work, at = time_work, twall = time_mo_ints ) 
            write(*,*) 'getting gmo'

            !call  wrapper_to_get_real_t1_free_gmo(nb,no,nv,Co%elm2,Cv2%elm2,&
            !  & iajb,ccmodel,mylsitem)
            call get_t1_free_gmo(mo_ccsd,mylsitem,Co%elm2,Cv2%elm2,iajb,pgmo_diag,pgmo_up, &
               & nb,no,nv,CCmodel,MOinfo)
            write(*,*) 'got gmo'

            if(DECinfo%PL>1)call time_start_phase( PHASE_work, at = time_work, ttot = time_mo_ints,&
               &labelttot = 'CCSOL: INIT MO INTS   :', output = DECinfo%output )

         end if

      end if INTEGRAL



      ! readme : the iteration sequence is universal and may be used for all
      !          iterative cc models (linear or non-linear) and is
      !          semi-independent on the storage of vectors (allocation and
      !          deallocation, etc)


      ! iterate
      break_iterations = .false.
      crop_ok          = .false.
      prev_norm        = 1.0E6_realk


      CCIteration : do iter=1,DECinfo%ccMaxIter

         if(DECinfo%PL>1)call time_start_phase( PHASE_work, at = time_work, twall = time_iter ) 

         ! remove old vectors
         RemoveOldVectors : if(iter > DECinfo%ccMaxDIIS) then

            if(DECinfo%cc_driver_debug) then
               write(DECinfo%output,'(a,i4)') ' debug :: vector to delete : ',iter-DECinfo%ccMaxDIIS
            end if

            if(DECinfo%use_singles) then
               call array_free( t1(iter-DECinfo%ccMaxDIIS)     )
               Call array_free( omega1(iter-DECinfo%ccMaxDIIS) )

            end if
            call array_free(t2(iter-DECinfo%ccMaxDIIS))
            call array_free(omega2(iter-DECinfo%ccMaxDIIS))

         end if RemoveOldVectors


         ! Initialize residual vectors
         if(DECinfo%use_singles)then
            omega1(iter) = array_minit( ampl2_dims, 2 , local=local, atype='LDAR' )
            call array_zero(omega1(iter))
         endif
         omega2(iter) = array_minit( ampl4_dims, 4, local=local, atype='TDAR' )
         call array_zero(omega2(iter))


         if(DECinfo%PL>1)call time_start_phase( PHASE_work, at = time_work, twall = time_t1_trafo ) 
         ! get singles
         T1Related : if(DECinfo%use_singles) then

            ! synchronize singles data on slaves
            call array_sync_replicated(t1(iter))

            ! get the T1 transformation matrices
            call array_cp_data(Cv2,yv)
            call array_cp_data(Cv,xv)
            call array_contract_outer_indices_rr(-1.0E0_realk,Co,t1(iter),1.0E0_realk,xv)

            call array_cp_data(Co2,yo)
            call array_cp_data(Co,xo)
            call array_contract_outer_indices_rl(1.0E0_realk,Cv2,t1(iter),1.0E0_realk,yo)

         end if T1Related

         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_t1_trafo, &
            &labelttot= 'CCIT: T1 TRAFO        :', output = DECinfo%output, twall = time_residual ) 

         ! readme : get residuals, so far this solver supports only singles and doubles
         !          amplitudes (extension to higher iterative model is trivial); differences
         !          are mainly based on the set of residual vectors to be evaluated

         ! MODIFY FOR NEW MODEL
         ! If you implement a new model, please insert call to your own residual routine here!
         SelectCoupledClusterModel : select case( CCmodel )
         case( MODEL_MP2 )

            call get_simple_parallel_mp2_residual(omega2(iter),iajb,t2(iter),ppfock_prec,qqfock_prec,iter,local)

         case( MODEL_CC2, MODEL_CCSD, MODEL_CCSDpT ) !CC2 or  CCSD or CCSD(T)


            call ccsd_residual_wrapper(ccmodel,w_cp,delta_fock,omega2(iter),t2(iter),&
               & fock,iajb,no,nv,ppfock,qqfock,pqfock,qpfock,xo,xv,yo,yv,nb,&
               & MyLsItem,omega1(iter),t1(iter),pgmo_diag,pgmo_up,MOinfo,&
               & mo_ccsd,iter,local,restart)

         case( MODEL_RPA )
           
            !msg = 'Norm of gmo'
            !call print_norm(iajb,msg)

            !msg = 'Norm of t2'
            !call print_norm(t2(iter),msg)

            !msg = 'Norm of omega2'
            !call print_norm(omega2(iter),msg)

           call RPA_residual_par(Omega2(iter),t2(iter),iajb,ppfock_prec,qqfock_prec,no,nv)

            !msg = 'Norm of omega2 after computations'
            !call print_norm(omega2(iter),msg)

            !call lsquit("ERROR(ccsolver_par):no RPA implemented",DECinfo%output)
         case default
            call lsquit("ERROR(ccsolver_par):wrong choice of ccmodel",DECinfo%output)
         end select SelectCoupledClusterModel

         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_residual, &
            &labelttot= 'CCIT: RESIDUAL        :', output = DECinfo%output, twall = time_crop_mat ) 

         ! calculate crop/diis matrix
         B=0.0E0_realk
         do i=iter,max(iter-DECinfo%ccMaxDIIS+1,1),-1
            do j=iter,i,-1
               if(DECinfo%use_singles) then
                  if(DECinfo%use_preconditioner_in_b) then
                     omega1_prec = precondition_singles( omega1(j), ppfock_prec, qqfock_prec        )
                     omega2_prec = precondition_doubles( omega2(j), ppfock_prec, qqfock_prec, local )
                     B(i,j) =          array_ddot( omega1(i), omega1_prec ) 
                     B(i,j) = B(i,j) + array_ddot( omega2(i), omega2_prec )

                     call array_free( omega1_prec )
                     call array_free( omega2_prec )
                  else
                     B(i,j) =          array_ddot( omega1(i), omega1(j) ) 
                     B(i,j) = B(i,j) + array_ddot( omega2(i), omega2(j) )
                  end if
               else
                  ! just doubles
                  if(DECinfo%use_preconditioner_in_b) then
                     omega2_prec = precondition_doubles(omega2(j),ppfock_prec,qqfock_prec,local)
                     B(i,j) = array_ddot( omega2(i), omega2_prec )
                     call array_free( omega2_prec )
                  else
                     B(i,j) = array_ddot( omega2(i), omega2(j) )
                  end if
               end if
               B(j,i) = B(i,j)
            end do
         end do

         !msg="DIIS mat, new"
         !call print_norm(B,DECinfo%ccMaxIter*DECinfo%ccMaxIter,msg)

         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_crop_mat, &
            &labelttot= 'CCIT: CROP MATRIX     :', output = DECinfo%output, twall = time_solve_crop ) 
         ! solve crop/diis equation
         c=0.0E0_realk
         call CalculateDIIScoefficients(DECinfo%ccMaxDIIS,DECinfo%ccMaxIter,iter,B,c, &
            DECinfo%cc_driver_debug)


         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_solve_crop, &
            &labelttot= 'CCIT: SOLVE CROP      :', output = DECinfo%output, twall = time_mixing ) 


         ! mixing omega to get optimal
         if(DECinfo%use_singles) then
            t1_opt     = array_minit( ampl2_dims, 2 , local=local, atype='LDAR')
            omega1_opt = array_minit( ampl2_dims, 2 , local=local, atype='LDAR')
            call array_zero(t1_opt    )
            call array_zero(omega1_opt)
         end if

         omega2_opt  = array_minit( ampl4_dims, 4, local=local, atype='TDAR' )
         t2_opt      = array_minit( ampl4_dims, 4, local=local, atype='TDAR' )
         call array_zero( omega2_opt )
         call array_zero( t2_opt     )

         do i=iter,max(iter-DECinfo%ccMaxDIIS+1,1),-1

            ! mix singles
            if(DECinfo%use_singles) then
               call array_add( omega1_opt, c(i), omega1(i) )
               call array_add( t1_opt,     c(i), t1(i)     )
            end if

            ! mix doubles
            call array_add( omega2_opt, c(i), omega2(i) )
            call array_add( t2_opt,     c(i), t2(i)     )

         end do


         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_mixing, &
            &labelttot= 'CCIT: MIXING          :', output = DECinfo%output, twall = time_copy_opt ) 

         ! if crop, put the optimal in place of trial (not for diis)
         if(DECinfo%use_crop) then
            if(DECinfo%use_singles) then
               call array_cp_data( omega1_opt, omega1(iter) )
               call array_cp_data( t1_opt,     t1(iter)     )
            end if
            call array_cp_data( omega2_opt, omega2(iter) )
            call array_cp_data( t2_opt,     t2(iter)     )
         end if

         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_copy_opt, &
            &labelttot= 'CCIT: COPY OPTIMALS   :', output = DECinfo%output, twall = time_norm ) 

         ! check for the convergence
         one_norm1 = 0.0E0_realk
         one_norm2 = 0.0E0_realk
         if(DECinfo%use_singles)call print_norm(omega1(iter),one_norm1,.true.)
         call print_norm(omega2(iter),one_norm2,.true.)
         one_norm_total = one_norm1 + one_norm2
         two_norm_total = sqrt(one_norm_total)
         test_norm      = two_norm_total

         !intentionally crash the calculation prematurely
         if(iter==5.and.DECinfo%CRASHCALC.and.DECinfo%full_molecular_cc)then
            print*,'Calculation was intentionally crashed due to keyword .CRASHCALC'
            print*,'This keyword is only used for debug and testing purposes'
            print*,'We want to be able to test the .RESTART keyword'
            print*,'In the CC case only quit prematurely, then this keyword is even more handy'
            WRITE(DECinfo%output,*)'Calculation was intentionally crashed due to keyword .CRASHCALC'
            WRITE(DECinfo%output,*)'This keyword is only used for debug and testing purposes'
            WRITE(DECinfo%output,*)'We want to be able to test the .RESTART keyword'
            print*,"SETTING TEST_NORM TO QUIT"
            test_norm=0.9*DECinfo%ccConvergenceThreshold
         endif

         ! simple crop diagnostics
         if(two_norm_total < prev_norm) then
            crop_ok=.true.
         else
            crop_ok=.false.
            write(DECinfo%output,'(a)') ' warning :: total norm was smaller in previous iteration !!! '
         end if
         prev_norm=two_norm_total
         ! check if this is the last iteration
         if(iter == DECinfo%ccMaxIter .or. test_norm < DECinfo%ccConvergenceThreshold) &
            break_iterations=.true.


         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_norm, &
            &labelttot= 'CCIT: GET NORMS       :', output = DECinfo%output, twall = time_energy ) 


         ! calculate the correlation energy and fragment energy
         ! MODIFY FOR NEW MODEL
         ! If you implement a new model, please insert call to energy routine here,
         ! or insert a call to get_cc_energy if your model uses the standard CC energy expression.
         EnergyForCCmodel: select case(CCmodel)
         case( MODEL_MP2 )

            ccenergy = get_mp2_energy(t2(iter),iajb,no,nv)

         case( MODEL_CC2, MODEL_CCSD, MODEL_CCSDpT )

            ! CC2, CCSD, or CCSD(T) (for (T) calculate CCSD contribution here)
            ccenergy = get_cc_energy(t1(iter),t2(iter),iajb,no,nv)

         case(MODEL_RPA)
           ccenergy = get_RPA_energy_arrnew(t2(iter),iajb,no,nv)
           if(DECinfo%SOS) then
             ccenergy =ccenergy+get_SOSEX_cont_arrnew(t2(iter),iajb,no,nv)
           endif


         case default

            call lsquit("ERROR(ccsolver_par):energy expression for your model&
               & not yet implemented",-1)

         end select EnergyForCCmodel

         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_energy, &
            &labelttot= 'CCIT: GET ENERGY      :', output = DECinfo%output, twall = time_new_guess) 

         ! generate next trial vector if this is not the last iteration
         if(.not.break_iterations) then
            if(DECinfo%use_preconditioner) then
               if(DECinfo%use_singles) then
                  omega1_prec = precondition_singles(omega1_opt,ppfock_prec,qqfock_prec)
                  t1(iter+1) = array_minit( ampl2_dims, 2, local=local, atype='REPD' )
                  call array_cp_data(t1_opt,t1(iter+1))
                  call array_add(t1(iter+1),1.0E0_realk,omega1_prec)
                  call array_free(omega1_prec)
               end if
               omega2_prec = precondition_doubles(omega2_opt,ppfock_prec,qqfock_prec,local)
               t2(iter+1) = array_minit( ampl4_dims, 4, local=local, atype='TDAR' )
               call array_cp_data(t2_opt,t2(iter+1))
               call array_add(t2(iter+1),1.0E0_realk,omega2_prec)
               call array_free(omega2_prec)
            else
               if(DECinfo%use_singles)then
                  t1(iter+1) = array_minit( ampl2_dims, 2, local=local, atype='REPD' )
                  call array_cp_data(t1_opt,t1(iter+1))
                  call array_add(t1(iter+1),1.0E0_realk,omega1_opt)
               endif
               t2(iter+1) = array_minit( ampl4_dims, 4, local=local, atype='TDAR' )
               call array_cp_data(t2_opt,t2(iter+1))
               call array_add(t2(iter+1),1.0E0_realk,omega2_opt)
            end if
            !if .not.DECinfo%CCSDnosaferun option is set, make sure data is in dense
            if(saferun)then

               if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, twall = time_write ) 

               if(DECinfo%use_singles)then
                  call save_current_guess(local,iter+old_iter,two_norm_total,ccenergy,t2(iter+1),safefilet21,&
                     &safefilet22,t1(iter+1),safefilet11,safefilet12)
               else
                  call save_current_guess(local,iter+old_iter,two_norm_total,ccenergy,t2(iter+1),safefilet21,&
                     &safefilet22)
               endif

               if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_write, &
                  &labelttot= ' -NG: NEW GUESS WRITE :', output = DECinfo%output) 

            endif
         end if

         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_new_guess, &
            &labelttot= 'CCIT: NEW GUESS VEC.  :', output = DECinfo%output ) 


         ! delete optimals
         if(DECinfo%use_singles) then
            call array_free(t1_opt)
            call array_free(omega1_opt)
         end if
         call array_free(t2_opt)
         call array_free(omega2_opt)

         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_iter, &
            &labelttot= 'CCIT: ITERATION       :', output = DECinfo%output) 

#ifdef __GNUC__
         call flush(DECinfo%output)
#endif

         call print_ccjob_iterinfo(iter+old_iter,two_norm_total,ccenergy,.false.,fragment_job)

         last_iter = iter
         if(break_iterations) exit

      end do CCIteration

   else

      call print_ccjob_iterinfo(old_iter,two_norm_total,ccenergy,.false.,fragment_job)

      break_iterations = .true.
      last_iter        = 1

   endif If_not_converged



   call time_start_phase( PHASE_work, at = time_work, ttot = time_main, labelttot = 'CCSOL: MAIN LOOP      :', &
      & output = DECinfo%output, twall = time_finalize )


   ! Free memory and save final amplitudes
   ! *************************************

   ! remove rest of the singles amplitudes and residuals
   do i=last_iter,max(last_iter-DECinfo%ccMaxDIIS+1,1),-1


      ! Save final double amplitudes (to file if saferun)
      if(i==last_iter) then
         ! Save two-electron integrals in the order (virt,occ,virt,occ)
         t2_final = array4_init([nv,nv,no,no])
         call array_convert(t2(last_iter),t2_final%val)
         call array4_reorder(t2_final,[1,3,2,4])
         !if(.not.restart_from_converged)then
         !   call array_cp_tiled2dense(t2(last_iter),.true.)
         !endif
         !call array_reorder_4d(1.0E0_realk,t2(last_iter)%elm1,nv,nv,no,no,[1,3,2,4],0.0E0_realk,t2_final%val)

         if(DECinfo%use_singles) then
            if(.not.longrange_singles) then ! intitialize and copy, else just copy
               t1_final = array2_init(ampl2_dims)
            end if
            call dcopy(int(t1(i)%nelms),t1(last_iter)%elm1,1,t1_final%val,1)
         endif

         !SAFE THE FINAL AMPLITUDES, NOT YET REORDERED
         if(saferun.and..not.restart_from_converged)then
            if(DECinfo%use_singles)then
               call save_current_guess(local,i+old_iter,two_norm_total,ccenergy,&
                  &t2(last_iter),safefilet2f,safefilet2f,t1(last_iter),safefilet1f,safefilet1f)
            else
               call save_current_guess(local,i+old_iter,two_norm_total,ccenergy,&
                  &t2(last_iter),safefilet2f,safefilet2f)
            endif
         endif

         !call array_change_itype_to_td(t2(last_iter),local)
      end if

      ! Free doubles residuals
      if(.not.restart_from_converged)call array_free(omega2(i))
      ! Free doubles amplitudes
      call array_free(t2(i))


      ! Free singles amplitudes and residuals
      if(DECinfo%use_singles) then
         call array_free( t1(i)     )
         if(.not.restart_from_converged)call array_free( omega1(i) )
      end if

   end do
   
   !in the call to get_combined_SingleDouble_amplitudes
   !t1 is used, for RPA use_singles = .false. 
   !it should be moved out of this routine
   !but I get a wrong energy if it is moved out of the
   !routine. Johannes 
   if(CCmodel == MODEL_RPA)then
     t1_final = array2_init([nv,no])
   endif

   call time_start_phase(PHASE_WORK,at = time_work, twall = ttotend_wall, tcpu = ttotend_cpu )

   ! Write finalization message
   call print_ccjob_summary(break_iterations,.false.,fragment_job,last_iter+old_iter,&
      &ccenergy,ttotend_wall,ttotstart_wall,ttotend_cpu,ttotstart_cpu,t1_final,t2_final)

   ! Save two-electron integrals in the order (virt,occ,virt,occ)
   VOVO = array4_init([no,nv,no,nv])
   call array_convert(iajb,VOVO%val)
   call array_free(iajb)
   call array4_reorder(VOVO,[2,1,4,3])
   call print_norm(VOVO%val,i8*no*no*nv*nv)

   ! deallocate stuff
   if(DECinfo%use_singles) then
      call mem_dealloc(t1)
      call mem_dealloc(omega1)
   end if

   call mem_dealloc(t2)
   call mem_dealloc(omega2)

   call mem_dealloc(B)
   call mem_dealloc(c)


   ! remove fock correction
   call array_free(delta_fock)


   call array_free(ppfock_prec)
   call array_free(qqfock_prec)

   if(DECinfo%use_singles) then
      !call array2_free(h1)
      call array_free(xo)
      call array_free(yo)
      call array_free(xv)
      call array_free(yv)
      call array_free(pqfock)
      call array_free(qpfock)
   end if

   
   if(ccmodel /= MODEL_MP2 .and. ccmodel /= MODEL_RPA)then
     call array_free(ppfock)
     call array_free(qqfock)
   endif

   call array_free(Co)
   call array_free(Co2)
   call array_free(Cv)
   call array_free(Cv2)
   call array_free(fock)

   ! free memory from MO-based CCSD
   if(.not. restart_from_converged)then
      if (mo_ccsd) then
         if (pgmo_diag%dims(2)>1) call array_free(pgmo_up)
         call array_free(pgmo_diag)
         call mem_dealloc(MOinfo%dimInd1)
         call mem_dealloc(MOinfo%dimInd2)
         call mem_dealloc(MOinfo%StartInd1)
         call mem_dealloc(MOinfo%StartInd2)
         call mem_dealloc(MOinfo%dimTot)
         call mem_dealloc(MOinfo%tileInd)
      end if
   endif

   !transform back to original basis   
   if(DECinfo%use_singles)then
      call ccsolver_can_local_trans(no,nv,nb,Uocc,Uvirt,vovo=t2_final%val,vo=t1_final%val)
      call ccsolver_can_local_trans(no,nv,nb,Uocc,Uvirt,vovo=VOVO%val)
   else
      call ccsolver_can_local_trans(no,nv,nb,Uocc,Uvirt,vovo=t2_final%val)
      call ccsolver_can_local_trans(no,nv,nb,Uocc,Uvirt,vovo=VOVO%val)
   endif

   call mem_dealloc(Uocc)
   call mem_dealloc(Uvirt)

#ifdef VAR_MPI
   if ( w_cp ) call lspdm_shut_down_comm_procs
   !print *,"ALL DONE"
   !call sleep(3)
   !stop 0
#endif

#ifdef MOD_UNRELEASED
   if( .not. fragment_job .and. DECinfo%PL>2 )then
      call array_print_mem_info(DECinfo%output,.true.,.false.)
   endif
#endif

   if(DECinfo%PL>1)call time_start_phase( PHASE_work, at = time_work, ttot = time_finalize, &
      &labelttot = 'CCSOL: FINALIZATION   :', output = DECinfo%output )

end subroutine ccsolver_par


!> \brief should be a general subroutine to get the guess amplitudes when
!starting up a CCSD or CC2 calculation and checks for files which contain
!amplitudes of a previous calculation. If none are found the usual zero guess
!is returned
!> \author Patrick Ettenhuber
!> \date December 2012
subroutine get_guess_vectors(restart,iter_start,norm,energy,t2,iajb,Co,Cv,oof,vvf,vof,mylsitem,local,&
   & safefilet21,safefilet22,safefilet2f, t1,safefilet11,safefilet12,safefilet1f)
   implicit none
   logical,intent(out) :: restart
   real(realk), intent(out) :: norm,energy
   !> contains the guess doubles amplitudes on output
   type(array), intent(inout) :: t2,iajb,Co,Cv,oof,vvf,vof
   logical, intent(in) :: local
   !> integral info
   type(lsitem), intent(inout) :: mylsitem
   !> the filenames to check for valid doubles amplitudes
   character(3),intent(in) :: safefilet21,safefilet22,safefilet2f
   !> contains the singles amplitudes on output
   type(array),intent(inout),optional :: t1
   !> the filenames to check for valid singles amplitudes
   character(3),intent(in), optional :: safefilet11,safefilet12,safefilet1f
   integer :: no,nv
   integer :: fu_t11,fu_t12,fu_t21,fu_t22,fu_t1,fu_t2,fu_t2f,fu_t1f
   logical :: file_exists11,file_exists12,file_exists1f,file_exists21,file_exists22,file_exists2f
   logical :: file_status11,file_status12,file_status1f,file_status21,file_status22,file_status2f
   logical :: readfile1, readfile2
   integer :: saved_iter11,saved_iter12,saved_iter1f,saved_iter21,saved_iter22,saved_iter2f,iter_start
   integer :: saved_nel11,saved_nel12,saved_nel21,saved_nel22,saved_nel1f,saved_nel2f
   logical :: all_singles, fin1_exists, fin2_exists
   character(11) :: fullname11, fullname12, fin1, fullname21, fullname22,fin2
   character(ARR_MSG_LEN) :: msg
   integer :: a,i

   all_singles=present(t1).and.present(safefilet11).and.present(safefilet12).and.present(safefilet1f)

   !print *,"CHECK INPUT",safefilet11,safefilet12,all_singles,DECinfo%use_singles
   fu_t11=111
   fu_t12=112
   fu_t1f=113
   fu_t21=121
   fu_t22=122
   fu_t2f=123
   nv=t2%dims(1)
   no=t2%dims(3)


   iter_start=0
   !check for safe files of the amplitudes in the current directory and read
   !them if they exist and ok
   readfile1 = .false.
   readfile2 = .false.

   if(DECinfo%DECrestart)then

      !CHECK IF THERE ARE CONVERGED AMPLITUDES AVAILABLE
      if(DECinfo%use_singles.and.all_singles)then
         fin1=safefilet1f//'.restart'
         INQUIRE(FILE=fin1,EXIST=file_exists1f)
         if(file_exists1f)then
            file_status1f=.true.
            OPEN(fu_t1f,FILE=fin1,STATUS='OLD',FORM='UNFORMATTED')
            READ(fu_t1f)saved_iter1f
            READ(fu_t1f)saved_nel1f
            if(saved_nel1f/=no*nv)then
               print *,"WARNING(ccsolver_par):wrong dimensions in final singles amplitudes file,&
                  & checking previous iterations"
            else
               fu_t1     = fu_t1f
               readfile1 = .true.
            endif
         endif
      endif

      fin2=safefilet2f//'.restart'
      INQUIRE(FILE=fin2,EXIST=file_exists2f)
      if(file_exists2f)then
         file_status2f=.true.
         OPEN(fu_t2f,FILE=fin2,STATUS='OLD',FORM='UNFORMATTED')
         READ(fu_t2f)saved_iter2f
         READ(fu_t2f)saved_nel2f
         if(saved_nel2f/=no**2*nv**2)then
            print *,"WARNING(ccsolver_par):wrong dimensions in final doubles amplitudes file,&
               & checking previous iterations"
         else
            fu_t2=fu_t2f
            readfile2=.true.
         endif
         iter_start = saved_iter2f
      endif

      !THEN CHECK IF THERE ARE AMPLITUDES FROM OTHER ITERATIONS AVAILALBE
      if(DECinfo%use_singles.and.all_singles.and..not.readfile1)then
         fullname11=safefilet11//'.restart'
         fullname12=safefilet12//'.restart'

         file_status11=.false.
         INQUIRE(FILE=fullname11,EXIST=file_exists11)
         if(file_exists11)then
            file_status11=.true.
            OPEN(fu_t11,FILE=fullname11,STATUS='OLD',FORM='UNFORMATTED')
            READ(fu_t11)saved_iter11
            READ(fu_t11)saved_nel11
            if(saved_nel11/=no*nv)then
               call lsquit("ERROR(ccsolver_par):wrong dimensions in amplitude &
                  &file",DECinfo%output)
            endif
         endif

         file_status12=.false.
         INQUIRE(FILE=fullname12,EXIST=file_exists12)
         if(file_exists12)then
            file_status12=.true.
            OPEN(fu_t12,FILE=fullname12,STATUS='OLD',FORM='UNFORMATTED')
            READ(fu_t12)saved_iter12
            READ(fu_t12)saved_nel12
            if(saved_nel12/=no*nv)then
               call lsquit("ERROR(ccsolver_par):wrong dimensions in amplitude &
                  &file",DECinfo%output)
            endif
         endif

         !CHECK WHICH IS THE PREFERRED FILE TO READ
         if(file_status11.and.file_status12)then
            if(saved_iter11>saved_iter12)then
               fu_t1=fu_t11
               CLOSE(fu_t12)
               readfile1=.true.
            else
               fu_t1=fu_t12
               CLOSE(fu_t11)
               readfile1=.true.
            endif
         else if(file_status11)then
            fu_t1=fu_t11
            readfile1=.true.
         else if(file_status12)then
            fu_t1=fu_t12
            readfile1=.true.
         endif  
      endif

      if(.not.readfile2)then

         fullname21=safefilet21//'.restart'
         fullname22=safefilet22//'.restart'

         file_status21=.false.
         INQUIRE(FILE=fullname21,EXIST=file_exists21)
         if(file_exists21)then
            file_status21=.true.
            OPEN(fu_t21,FILE=fullname21,STATUS='OLD',FORM='UNFORMATTED')
            READ(fu_t21)saved_iter21
            READ(fu_t21)saved_nel21
            if(saved_nel21/=no*no*nv*nv)then
               call lsquit("ERROR(ccsolver_par):wrong dimensions in amplitude &
                  &file",DECinfo%output)
            endif
         endif

         file_status22=.false.
         INQUIRE(FILE=fullname22,EXIST=file_exists22)
         if(file_exists22)then
            file_status22=.true.
            OPEN(fu_t22,FILE=fullname22,STATUS='OLD',FORM='UNFORMATTED')
            READ(fu_t22)saved_iter22
            READ(fu_t22)saved_nel22
            if(saved_nel22/=no*no*nv*nv)then
               call lsquit("ERROR(ccsolver_par):wrong dimensions in amplitude &
                  &file",DECinfo%output)
            endif
         endif

         !CHECK WHICH IS THE PREFERRED FILE TO READ
         if(file_status21.and.file_status22)then
            if(saved_iter21>saved_iter22)then
               iter_start=saved_iter21
               fu_t2=fu_t21
               CLOSE(fu_t22)
               readfile2=.true.
            else
               iter_start=saved_iter22
               fu_t2=fu_t22
               CLOSE(fu_t21)
               readfile2=.true.
            endif
            WRITE(DECinfo%output,'("RESTARTING CC CALCULATION WITH TRIAL-VECS FROM: ",I3)')iter_start
         else if(file_status21)then
            iter_start=saved_iter21
            fu_t2=fu_t21
            WRITE(DECinfo%output,'("RESTARTING CC CALCULATION WITH TRIAL-VECS FROM: ",I3)')iter_start
            readfile2=.true.
         else if(file_status22)then
            iter_start=saved_iter22
            fu_t2=fu_t22
            WRITE(DECinfo%output,'("RESTARTING CC CALCULATION WITH TRIAL-VECS FROM: ",I3)')iter_start
            readfile2=.true.
         else
            iter_start=0
         endif  
      endif

   endif


   if(readfile1)then
      READ(fu_t1) t1%elm1
      READ(fu_t1) norm
      READ(fu_t1) energy
      CLOSE(fu_t1)
      restart = .true.
   else
      !do a=1,nv
      !   do i=1,no

      !      t1%elm2(a,i) = vof%elm2(a,i)/( oof%elm2(i,i) - vvf%elm2(a,a) ) 

      !   end do
      !end do
      if(DECinfo%use_singles) call array_zero(t1)
   endif

   if(readfile2)then
      ! allocate dense part of t2 array:
      if (.not.local) call memory_allocate_array_dense(t2)
      READ(fu_t2) t2%elm1
      READ(fu_t2) norm
      READ(fu_t2) energy
      CLOSE(fu_t2)
      ! mv dense part to tiles:
      if (.not.local) call array_mv_dense2tiled(t2,.false.)
      restart = .true.
   else
      !call get_mo_integral_par( iajb, Co, Cv, Co, Cv, mylsitem, local, .true.)
      !call get_mp2_starting_guess( iajb, t2, oof, vvf, local )
      call array_zero(t2)
      !call array_zero(iajb)
   endif
end subroutine get_guess_vectors

!> \brief Subroutine to save the current guess amplitudes for the next
!iteration
!> \author Patrick Ettenhuber
!> \date Dezember 2012
subroutine save_current_guess(local,iter,res_norm,energy,t2,safefilet21,safefilet22,&
      &t1,safefilet11,safefilet12)
   implicit none
   logical, intent(in) :: local
   !> iteration number
   integer,intent(in) :: iter
   !> write the corresponding residual norm into the file
   real(realk), intent(in)    :: res_norm,energy
   !> doubles guess amplitudes for the next iteration
   type(array), intent(inout) :: t2
   !> alternating filenames for the doubles amplitudes
   character(3),intent(in)    :: safefilet21,safefilet22
   !> singles guess amplitudes for the next iteration
   type(array), intent(in), optional :: t1
   !> alternating filenames for the singles amplitudes
   character(3),intent(in), optional :: safefilet11,safefilet12
   integer :: fu_t21,fu_t22
   integer :: fu_t11,fu_t12
   logical :: file_status11,file_status12,file_status21,file_status22
   logical :: all_singles
   character(ARR_MSG_LEN) :: msg
#ifdef SYS_AIX
   character(12) :: fullname11,  fullname12,  fullname21,  fullname22
   character(12) :: fullname11D, fullname12D, fullname21D, fullname22D
#else
   character(11) :: fullname11, fullname12, fullname21, fullname22
   character(11) :: fullname11D, fullname12D, fullname21D, fullname22D
#endif
   ! cp doubles from tile to dense part: (only if t2%itype/=DENSE)
   if (.not.local) call array_cp_tiled2dense(t2,.false.)

   all_singles=present(t1).and.present(safefilet11).and.present(safefilet12)
   fu_t11=111
   fu_t12=112
   fu_t21=121
   fu_t22=122
   if(DECinfo%use_singles.and.all_singles)then
      !msg="singles norm save"
      !call print_norm(t1,msg)
#ifdef SYS_AIX
      fullname11=safefilet11//'.writing\0'
      fullname12=safefilet12//'.writing\0'
#else
      fullname11=safefilet11//'.writing'
      fullname12=safefilet12//'.writing'
#endif

      if(mod(iter,2)==1)then
         file_status11=.false. 
         OPEN(fu_t11,FILE=fullname11,STATUS='REPLACE',FORM='UNFORMATTED')
         WRITE(fu_t11)iter
         WRITE(fu_t11)t1%nelms
         WRITE(fu_t11)t1%elm1
         WRITE(fu_t11)res_norm
         WRITE(fu_t11)energy
         file_status11=.true.
         WRITE(fu_t11)file_status11
         ENDFILE(fu_t11)
         CLOSE(fu_t11)
#ifdef SYS_AIX
         fullname11D=safefilet11//'.restart\0'
         fullname11=safefilet11//'\0'
#else
         fullname11D=safefilet11//'.restart'
#endif
         if(file_status11)call rename(fullname11,fullname11D)

      else if(mod(iter,2)==0)then
         file_status12=.false. 
         OPEN(fu_t12,FILE=fullname12,STATUS='REPLACE',FORM='UNFORMATTED')
         WRITE(fu_t12)iter
         WRITE(fu_t12)t1%nelms
         WRITE(fu_t12)t1%elm1
         WRITE(fu_t12)res_norm
         WRITE(fu_t12)energy
         file_status12=.true.
         WRITE(fu_t12)file_status12
         ENDFILE(fu_t12)
         CLOSE(fu_t12)
#ifdef SYS_AIX
         fullname12D=safefilet12//'.restart\0'
         fullname12=safefilet12//'\0'
#else
         fullname12D=safefilet12//'.restart'
#endif
         if(file_status12)call rename(fullname12,fullname12D)

      else
         call lsquit("ERROR(ccdriver_par):impossible iteration&
            &number)",DECinfo%output)
      endif
   endif

   fullname21=safefilet21//'.writing'
   fullname22=safefilet22//'.writing'
   if(mod(iter,2)==1)then
      file_status21=.false. 
      OPEN(fu_t21,FILE=fullname21,STATUS='REPLACE',FORM='UNFORMATTED')
      WRITE(fu_t21)iter
      WRITE(fu_t21)t2%nelms
      WRITE(fu_t21)t2%elm1
      WRITE(fu_t21)res_norm
      WRITE(fu_t21)energy
      file_status21=.true. 
      WRITE(fu_t22)file_status21
      ENDFILE(fu_t21)
      CLOSE(fu_t21)
#ifdef SYS_AIX
      fullname21D=safefilet21//'.restart\0'
      fullname21=safefilet21//'\0'
#else
      fullname21D=safefilet21//'.restart'
#endif
      if(file_status21)call rename(fullname21,fullname21D)

   else if(mod(iter,2)==0)then
      file_status22=.false. 
      OPEN(fu_t22,FILE=fullname22,STATUS='REPLACE',FORM='UNFORMATTED')
      WRITE(fu_t22)iter
      WRITE(fu_t22)t2%nelms
      WRITE(fu_t22)t2%elm1
      WRITE(fu_t22)res_norm
      WRITE(fu_t22)energy
      file_status22=.true.
      WRITE(fu_t22)file_status22
      ENDFILE(fu_t22)
      CLOSE(fu_t22)
#ifdef SYS_AIX
      fullname22D=safefilet22//'.restart\0'
      fullname22=safefilet22//'\0'
#else
      fullname22D=safefilet22//'.restart'
#endif
      if(file_status22)call rename(fullname22,fullname22D)
   else
      call lsquit("ERROR(ccdriver_par):impossible iteration&
         &number)",DECinfo%output)
   endif
   ! deallocate dense part of doubles:
   if (.not.local) call memory_deallocate_array_dense(t2)
end subroutine save_current_guess

!> Purpose: Wrapper for the RPA model: get MO integrals (non-T1 transformed)
!
!> Author:  Pablo Baudin
!> Date:    January 2014
subroutine wrapper_to_get_real_t1_free_gmo(nb,no,nv,Co,Cv,govov,ccmodel,mylsitem)

  implicit none

  integer, intent(in) :: nb, no, nv
  real(realk), pointer, intent(in) :: Co(:,:), Cv(:,:)
  type(array), intent(inout) :: govov
  integer, intent(in) :: ccmodel
  !> LS item with information needed for integrals
  type(lsitem), intent(inout) :: MyLsItem
     
  ! dummy arguments:
  type(array) :: pgmo_diag, pgmo_up
  type(MObatchInfo) :: MOinfo
  logical :: mo_ccsd  
  
  call get_t1_free_gmo(mo_ccsd,mylsitem,Co,Cv,govov,pgmo_diag,pgmo_up, &
    & nb,no,nv,CCmodel,MOinfo)
 
  end subroutine wrapper_to_get_real_t1_free_gmo

end module ccdriver
