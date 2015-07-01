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
use lsparameters
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
use cc_response_tools_module
use ccsd_module!,only: getDoublesResidualMP2_simple, &
!       & getDoublesResidualCCSD_simple,getDoublesResidualCCSD_simple2, &
!       & precondition_doubles,get_ccsd_residual_integral_driven,&
!       & get_ccsd_residual_integral_driven_oldtensor_wrapper
use pno_ccsd_module
use snoop_tools_module
#ifdef MOD_UNRELEASED
use ccsdpt_module
!endif mod_unreleased
#endif
use orbital_operations
use rpa_module



public :: ccsolver, fragment_ccsolver, ccsolver_justenergy,&
   & mp2_solver, SOLVE_AMPLITUDES,SOLVE_AMPLITUDES_PNO, SOLVE_MULTIPLIERS,&
   & ccsolver_energy_multipliers
private

interface mp2_solver
   module procedure mp2_solver_frag, mp2_solver_mol
end interface mp2_solver


contains

#ifdef MOD_UNRELEASED
   subroutine ccsolver_energy_multipliers(ccmodel,Co,Cv,fock,nb,no,nv, &
        &mylsitem,ccPrintLevel,fragment_job,oof,vvf,ccenergy)

     implicit none

     !> CC model
     integer,intent(in) :: ccmodel
     !> Number of occupied orbitals in full molecule/fragment AOS
     integer, intent(in) :: no
     !> Number of virtual orbitals in full molecule/fragment AOS
     integer, intent(in) :: nv
     !> Number of basis functions in full molecule/atomic extent
     integer, intent(in) :: nb
     !> Fock matrix in AO basis for fragment or full molecule
     real(realk), dimension(nb,nb), intent(inout) :: fock
     !> Occupied MO coefficients  for fragment/full molecule
     real(realk), dimension(nb,no), intent(inout) :: Co
     !> Virtual MO coefficients  for fragment/full molecule
     real(realk), dimension(nb,nv), intent(inout) :: Cv
     !> Occ-occ block of Fock matrix in MO basis
     real(realk), dimension(no,no), intent(inout) :: oof
     !> Virt-virt block of Fock matrix in MO basis
     real(realk), dimension(nv,nv), intent(inout) :: vvf
     !> Is this a fragment job (true) or a full molecular calculation (false)
     logical, intent(in) :: fragment_job
     !> LS item information
     type(lsitem), intent(inout) :: mylsitem
     !> How much to print? ( ccPrintLevel>0 --> print info stuff)
     integer, intent(in) :: ccPrintLevel
     !> Coupled cluster energy for fragment/full molecule
     real(realk),intent(inout) :: ccenergy!,ccsdpt_e4,ccsdpt_e5,ccsdpt_tot
     type(array4) :: t2inal,VOVO!,ccsdpt_t2
     type(array2) :: t1inal!,ccsdpt_t1
     !> stuff needed for pair analysis
     type(array2) :: ccsd_mat_tot,ccsd_mat_tmp
     integer :: natoms,ncore,no_tot,p,pdx,i
     type(decorbital), pointer :: occ_orbitals(:)
     type(decorbital), pointer :: uno_orbitals(:)
     logical, pointer :: orbitals_assigned(:)
     type(array4) :: mult2
     type(array2) :: mult1
     type(tensor) :: t1f,t2f,m1f,m2f,aibj
     integer :: solver_ccmodel
     logical :: loc

     solver_ccmodel = ccmodel
     if(ccmodel == MODEL_CCSDpT)solver_ccmodel = MODEL_CCSD
     loc=.true.
#ifdef VAR_MPI
     loc = (infpar%lg_nodtot == 1)
#endif

     call ccsolver_job(ccmodel,Co,Cv,fock,nb,no,nv,mylsitem,ccPrintLevel,&
        &oof,vvf,ccenergy,aibj,.false.,loc,t1f,t2f,m1f,m2f)

     !call one_el_unrel_dens(nb,no,nv,t1f,t2f,m1f,m2f) 

     call tensor_free(t1f)
     call tensor_free(t2f)
     call tensor_free(m1f)
     call tensor_free(m2f)
     call tensor_free(aibj)


   end subroutine ccsolver_energy_multipliers
#endif

   subroutine one_el_unrel_dens(nb,no,nv,t1f,t2f,m1f,m2f)

      implicit none
      integer, intent(in) :: nb,no,nv
      type(tensor), intent(inout) :: t1f,t2f
      type(tensor), intent(inout) :: m1f,m2f
      integer a,b,i,j
      integer, dimension(2)       :: dims

      type(array4)                :: tbar, temp1, temp2, temp3
      type(array2)                :: Xij, Yba, Eai, temp
      type(array2) :: Dsd
      
      dims = [nb,nb]

      !Initialize SD density
      Dsd  = array2_init(dims)

      !1. virt-occ block
      !*****************
      !D^sd_ai = z^a_i
      do a=no+1,nb
         do i=1,no
             Dsd%val(a,i) = m1f%elm2(a-no,i)
         end do
      end do
 
      !2. virt-virt block
      !******************
      !D^sd_ab = Sum_ckl (t^cb_kl * z^ca_kl) = Y_ba
 
      !container for amps voov shape
      temp1 = array4_init([nv,no,no,nv])
      temp2 = array4_init([nv,no,no,nv])
 
      !sort amps(ckbl) -> bckl
      call array_reorder_4d(1.0E0_realk,m2f%elm4,nv,no,nv,no,[1,2,4,3],0.0E0_realk,temp1%val)
      call array_reorder_4d(1.0E0_realk,t2f%elm4,nv,no,nv,no,[1,2,4,3],0.0E0_realk,temp2%val)
 
      !form Y_ba intermediate
      Yba = array2_init([nv,nv])
      call array4_contract3(temp1,temp2,Yba)
 
      !clean up!
      call array4_free(temp1)
      call array4_free(temp2)
 
      !D^sd_ab = Yba
      do a=no+1,nb
         do b=no+1,nb
             Dsd%val(a,b) = 1.0*Yba%val(a-no,b-no)
         end do
      end do
      call array2_free(Yba)

      !3. occ-occ block
      !****************
      !D^sd_ij = 2*delta_ij - Sum_cdk (t^cd_ki * z^cd_kj) = 2*delta_ij - X_ij
 
      !container for amps vovo shape
      temp1 = array4_init([nv,no,nv,no])
      temp2 = array4_init([nv,no,nv,no])
      call array_reorder_4d(1.0E0_realk,m2f%elm4,nv,no,nv,no,[1,2,3,4],0.0E0_realk,temp1%val)
      call array_reorder_4d(1.0E0_realk,t2f%elm4,nv,no,nv,no,[1,2,3,4],0.0E0_realk,temp2%val)
 
      !form X intermediate
      Xij = array2_init([no,no])
      call array4_contract3(temp2,temp1,Xij)
 
      !clean up!
      call array4_free(temp1)
      call array4_free(temp2)
 
      !add to the SD density
      do i=1,no
         do j=1,no
            if(i==j)then
               Dsd%val(i,j) = Dsd%val(i,j) + 2 - Xij%val(i,j)
            else
               Dsd%val(i,j) = Dsd%val(i,j) - Xij%val(i,j)
            endif
         end do
      end do
      call array2_free(Xij)

      !4. occ-virt block
      !*****************
      !D^sd_ia = Sum_bj z^b_j*(2t^ab_ij - t^ba_ij) = Sum_bj z^b_j * tbar^ab_ij
 
      !intermediate container vovo shape
      tbar = array4_init([nv,no,nv,no])
 
      !Construct the tbar intermediate
      do j=1,t2f%dims(4)
        do b=1,t2f%dims(3)
           do i=1,t2f%dims(2)
              do a=1,t2f%dims(1)
                 tbar%val(a,i,b,j) = tbar%val(a,i,b,j) + 2.0*t2f%elm4(a,i,b,j)&
                                      & - t2f%elm4(b,i,a,j)
              end do
           end do
        end do
      end do
 
      !form eai intermediate
      Eai = array2_init([nv,no])
 
      !container for amps
      temp = array2_init([nv,no],m1f%elm2)
      call array4_contract_array2(tbar,temp,Eai)
      call array2_free(temp)
      call array4_free(tbar)

      !add to the SD density
      do a=no+1,nb
         do i=1,no
             Dsd%val(i,a) = Eai%val(a-no,i)
         end do
      end do

      !print *, ""
      !print *, "Dsd in MO DEBUG Print: "
      !print *, ""
      !print *,no,nv
      !write(*, *) ''
      !do i=1,Dsd%dims(1)
      !  do j=1,Dsd%dims(2)
      !    write(*,'(f16.10)',advance='no') Dsd%val(i,j)
      !  end do
      !  write(*, *) ''
      !end do
      !print *, ""
  
   end subroutine one_el_unrel_dens  

   subroutine ccsolver_job(ccmodel,Co,Cv,fock,nb,no,nv,myls,ccPr,oof,vvf,e,vovo,lrt1,loc,t1f,t2f,m1f,m2f,frag)
      implicit none
      integer, intent(in) :: ccmodel, nb,no,nv, ccPr
      real(realk), intent(out) :: e
      type(lsitem), intent(inout) :: myls
      logical, intent(in)    :: lrt1
      logical, intent(inout) :: loc
      real(realk), intent(inout) :: Co(:,:), Cv(:,:), oof(:,:),  vvf(:,:), fock(:,:)
      type(tensor), intent(inout) :: vovo,t1f,t2f
      type(tensor), intent(inout), optional :: m1f,m2f
      type(decfrag), intent(inout), optional :: frag
      logical :: fj
      integer :: os,vs, solver_job, solver_ccmodel
      type(tensor) :: mp2_amp
      type(array4) :: t2a4,VOVOa4, m2a4, mp2a4
      type(array2) :: t1a2, m1a2

      solver_ccmodel = ccmodel
      if(ccmodel == MODEL_CCSDpT) solver_ccmodel = MODEL_CCSD
      fj = present(frag)
      call get_symm_tensor_segmenting_simple(no,nv,os,vs)


      if(DECinfo%use_pnos)then

         call ccsolver(MODEL_MP2,Co,Cv,fock,nb,no,nv,myls,ccPr,&
            &e,vovo,lrt1,loc,SOLVE_AMPLITUDES,p4=mp2_amp,frag=frag)

         if(fj)then
            !GET THE MP2 CORRELATION DENSITY FOR THE CENTRAL ATOM
            call calculate_MP2corrdens_frag(mp2_amp,frag) 
         endif

         solver_job = SOLVE_AMPLITUDES_PNO
      else
         solver_job = SOLVE_AMPLITUDES
      endif

      call ccsolver(solver_ccmodel,Co,Cv,fock,nb,no,nv,myls,ccPr,e,&
         & vovo,lrt1,loc,solver_job, p2=t1f, p4=t2f, m4=mp2_amp, vovo_supplied=DECinfo%use_pnos, frag=frag )

      if(DECinfo%use_pnos) call tensor_free( mp2_amp )

      if(DECinfo%CCSDmultipliers)then
         call ccsolver(solver_ccmodel,Co,Cv,fock,nb,no,nv,myls,ccPr,e,&
            & vovo,lrt1,loc,SOLVE_MULTIPLIERS, p2=m1f, p4=m2f, m2=t1f, m4=t2f,vovo_supplied=.true., frag=frag )
      endif

      ! Calculate CCSD eigenvalues if requested
      if(DECinfo%CCeival) then
         call ccsd_eigenvalue_solver(nb,no,nv,fock,Co,Cv,myls,t1f,t2f)
      end if

   end subroutine ccsolver_job

!> \brief get ccsd(t) corrections for full molecule.
!> \author Janus Juul Eriksen, modified by Patrick Ettenhuber and TK
!> \date February 2013
function ccsolver_justenergy(ccmodel,MyMolecule,nbasis,nocc,nvirt,mylsitem,&
      & ccPrintLevel,fragment_job,Co_fc,ppfock_fc) result(ccenergy)

   implicit none

   !> CC model
   integer,intent(inout) :: ccmodel
   !> full molecule information
   type(fullmolecule), intent(inout) :: MyMolecule
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
   real(realk), dimension(nbasis,nocc), intent(in),optional,target :: Co_fc
   !> Occ-occ block of Fock matrix in MO basis (only used for frozen core)
   real(realk), dimension(nocc,nocc), intent(in),optional,target :: ppfock_fc
   !type(array2) :: t1_final_arr2
   !type(array4) :: t2_final_arr4, VOVO_arr4
   type(tensor) :: t2_final,ccsdpt_t2,VOVO, mp2_amp
   type(tensor) :: t1_final,ccsdpt_t1,e4_mat_tot,e4_mat_tmp,e5_mat_tot
   type(tensor) :: t2f_local, VOVO_local
   type(tensor) :: m1_final, m2_final
   integer :: natoms,nfrags,ncore,nocc_tot,p,pdx,i
   type(decorbital), pointer :: occ_orbitals(:)
   type(decorbital), pointer :: virt_orbitals(:)
   logical, pointer :: orbitals_assigned(:)
   logical :: local,print_frags,abc
   ! Fragment and total energies as listed in decfrag type def "energies"
   real(realk), pointer :: FragEnergies(:,:,:), FragEnergies_tmp(:,:), ccenergies(:)
   real(realk) :: ccenergy
   integer :: nenergies, cc_sol_o, pT_4_o, pT_5_o, pT_full_o
   integer :: cc_sol_v, pT_4_v, pT_5_v, pT_full_v
   real(realk), pointer :: Co(:,:), Cv(:,:), oof(:,:),  vvf(:,:), fock(:,:)
   !> local variables 
   character(len=30) :: CorrEnergyString
   integer :: iCorrLen, nsingle, npair, njobs
   real(realk) :: norm_pt_1,norm_pt_2
   logical :: bg

   real(realk) :: time_CCSD_work, time_CCSD_comm, time_CCSD_idle
   real(realk) :: time_pT_work, time_pT_comm, time_pT_idle

   call time_start_phase(PHASE_WORK, swwork = time_CCSD_work, swcomm = time_CCSD_comm, swidle = time_CCSD_idle) 

   local= .true.
   bg   = mem_is_background_buf_init()
#ifdef VAR_MPI
   if(infpar%lg_nodtot>1)local=.false.
#endif
   Cv   => MyMolecule%Cv%elm2
   fock => MyMolecule%fock%elm2
   vvf  => MyMolecule%vvfock%elm2

   ! is this a frozen core calculation or not?
   if (DECinfo%frozencore) then

      ncore = MyMolecule%ncore

      if(.not. present(Co_fc)) then
         call lsquit('ccsolver_justenergy_pt: Occ MOs not present for frozencore!',-1)
      end if

      if(.not. present(ppfock_fc)) then
         call lsquit('ccsolver_justenergy_pt: Occ-occ Fock matrix not present for frozencore!',-1)
      end if

      Co  => Co_fc
      oof => ppfock_fc


   else
      ncore = 0
      Co => MyMolecule%Co%elm2
      oof => MyMolecule%oofock%elm2
   end if

#ifdef MOD_UNRELEASED

   ! nenergies is set to 8: a CC solver model plus pT corrections, 
   ! (4th order, 5th order and both), for both virt and occ partionings
   nenergies = 8
   call mem_alloc(ccenergies,nenergies)
   ccenergies = 0.0E0_realk
   cc_sol_o  = 1
   cc_sol_v  = 2
   pT_full_o = 3
   pT_full_v = 4
   pT_4_o    = 5
   pT_4_v    = 6
   pT_5_o    = 7
   pT_5_v    = 8

   print_frags = DECinfo%print_frags
   abc = DECinfo%abc
   
   if (print_frags) then ! should we print fragment energies?

      if (.not. DECinfo%pt_hack) then

         call ccsolver_job(ccmodel,Co,Cv,fock,nbasis,nocc,nvirt,mylsitem,ccPrintLevel, &
            & oof,vvf,ccenergies(cc_sol_o),VOVO,.false.,local,t1_final,t2_final, &
            & m1f=m1_final, m2f=m2_final)
   
         if(DECinfo%CCSDmultipliers)then
            call tensor_free(m2_final)
            if(DECinfo%use_singles)then
               call tensor_free(m1_final)
            endif
         endif
      
         if(DECinfo%PL>1)then
            call time_start_phase(PHASE_WORK, dwwork = time_CCSD_work, dwcomm = time_CCSD_comm, dwidle = time_CCSD_idle, &
               &swwork = time_pT_work, swcomm = time_pT_comm, swidle = time_pT_idle, &
               &labeldwwork = 'MASTER WORK CC solver: ',&
               &labeldwcomm = 'MASTER COMM CC solver: ',&
               &labeldwidle = 'MASTER IDLE CC solver: ') 
         endif
      
         natoms   = MyMolecule%natoms
         nfrags   = MyMolecule%nfrags
         nocc_tot = MyMolecule%nocc

      endif

      if(ccmodel == MODEL_CCSDpT)then
 
         if (abc) then

            call tensor_init(ccsdpt_t1 , [nocc,nvirt],2, bg=bg)
            call tensor_init(ccsdpt_t2 , [nocc,nocc,nvirt,nvirt],4, bg=bg)

         else
 
            call tensor_init(ccsdpt_t1, [nvirt,nocc],2,bg=bg)
            call tensor_init(ccsdpt_t2, [nvirt,nvirt,nocc,nocc],4,bg=bg)
   
         endif

         call ccsdpt_driver(nocc,nvirt,nbasis,oof,vvf,Co,Cv,mylsitem,VOVO,t2_final,&
            & ccsdpt_t1,print_frags,abc,ccsdpt_doubles=ccsdpt_t2)

         if (DECinfo%pt_hack) then

            call tensor_print_norm_nrm(ccsdpt_t1,norm_pt_1)
            call tensor_print_norm_nrm(ccsdpt_t2,norm_pt_2)
            call tensor_free(ccsdpt_t1)
            call tensor_free(ccsdpt_t2)
            call mem_dealloc(ccenergies)
            ccenergy = (norm_pt_1 + norm_pt_2) * 1.0E-3_realk ! arbitrary value - only for testing...

            return

         endif
  
         ! now, reorder amplitude and integral arrays
         if (abc) then

            call tensor_reorder(ccsdpt_t1,[2,1]) ! order (i,a) --> (a,i)
            call tensor_reorder(ccsdpt_t2,[3,4,1,2]) ! order (i,j,a,b) --> (a,b,i,j)
     
         endif
 
         if(DECinfo%PL>1)then
            call time_start_phase(PHASE_WORK,dwwork = time_pT_work, dwcomm = time_pT_comm, dwidle = time_pT_idle, &
               &labeldwwork = 'MASTER WORK pT: ',&
               &labeldwcomm = 'MASTER COMM pT: ',&
               &labeldwidle = 'MASTER IDLE pT: ') 
         endif
   
      endif
  
      !FIXME: all the following should be implemented in PDM
      !THIS IS JUST A WORKAROUND, ccsolver_par gives PDM tensors if more than
      !one node is used
      call tensor_init(VOVO_local, VOVO%dims,     4, bg=bg)
      call tensor_init(t2f_local,  t2_final%dims, 4, bg=bg)
      call tensor_cp_data( VOVO,     VOVO_local  )
      call tensor_cp_data( t2_final, t2f_local   )
      call tensor_reorder(t2f_local,[1,3,2,4]) ! t2 in the order (a,b,i,j)
      call tensor_reorder(VOVO_local,[1,3,2,4]) ! vovo integrals in the order (a,b,i,j)

      ! as we want to  print out fragment and pair interaction fourth-order energy contributions,
      ! then for locality analysis purposes we need occ_orbitals and
      ! virt_orbitals (adapted from fragment_energy.f90)
   
      ! -- Analyze basis and create orbitals
      call mem_alloc(occ_orbitals,nocc_tot)
      call mem_alloc(virt_orbitals,nvirt)
      call GenerateOrbitals_driver(MyMolecule,mylsitem,nocc_tot,nvirt,natoms, &
         & occ_orbitals,virt_orbitals)
   
      ! Orbital assignment
      call mem_alloc(orbitals_assigned,nfrags)
      orbitals_assigned=.false.
      if (DECinfo%onlyoccpart) then
         do p=1,nocc_tot
            pdx = occ_orbitals(p)%centralatom
            orbitals_assigned(pdx) = .true.
         end do
      else if (DECinfo%onlyvirtpart) then
         do p=1,nvirt
            pdx = virt_orbitals(p)%centralatom
            orbitals_assigned(pdx) = .true.
         end do
      else
         do p=1,nocc_tot
            pdx = occ_orbitals(p)%centralatom
            orbitals_assigned(pdx) = .true.
         end do
         do p=1,nvirt
            pdx = virt_orbitals(p)%centralatom
            orbitals_assigned(pdx) = .true.
         end do
      end if

      ! Summary print out
      nsingle = count(orbitals_assigned)
      npair   = nsingle*(nsingle-1)/2
      njobs   = nsingle + npair
      write(DECinfo%output,'(//,1X,a,i10)') 'FULL JOB SUMMARY: Number of single jobs = ', nsingle
      write(DECinfo%output,'(1X,a,i10)') 'FULL JOB SUMMARY: Number of pair jobs   = ', npair
      write(DECinfo%output,'(1X,a,i10,//)') 'FULL JOB SUMMARY: Total number of jobs  = ', njobs
   

      ! Calculate single and pair fragments energies:
      call mem_alloc(FragEnergies,nfrags,nfrags,nenergies)
      call mem_alloc(FragEnergies_tmp,nfrags,nfrags)
      FragEnergies     = 0.0E0_realk
      FragEnergies_tmp = 0.0E0_realk
   
      if (DECinfo%DECNP) then
         call solver_decnp_full(nocc,nvirt,nfrags,ncore,t2f_local,t1_final, &
            & VOVO_local,occ_orbitals,virt_orbitals,FragEnergies(:,:,cc_sol_o), &
            & FragEnergies(:,:,cc_sol_v),FragEnergies_tmp)
      else
         call solver_energy_full(nocc,nvirt,nfrags,ncore,t2f_local,t1_final, &
            & VOVO_local,occ_orbitals,virt_orbitals,FragEnergies(:,:,cc_sol_o), &
            & FragEnergies(:,:,cc_sol_v),FragEnergies_tmp)
      end if

   
      if(ccmodel == MODEL_CCSDpT)then
         if (DECinfo%DECNP) then
            ! now we calculate fourth-order and fifth-order energies
            call ccsdpt_decnp_e4_full(nocc,nvirt,nfrags,ncore,t2f_local,ccsdpt_t2,occ_orbitals,&
               & virt_orbitals,FragEnergies(:,:,pT_4_o),FragEnergies(:,:,pT_4_v), &
               & FragEnergies_tmp,ccenergies(pT_4_o))

            call ccsdpt_decnp_e5_full(nocc,nvirt,nfrags,ncore,t1_final,ccsdpt_t1,&
               & occ_orbitals,virt_orbitals,FragEnergies(:,:,pT_5_o),FragEnergies(:,:,pT_5_v), &
               & ccenergies(pT_5_o))
         else
            ! now we calculate fourth-order and fifth-order energies
            call ccsdpt_energy_e4_full(nocc,nvirt,nfrags,ncore,t2f_local,ccsdpt_t2,occ_orbitals,&
               & virt_orbitals,FragEnergies(:,:,pT_4_o),FragEnergies(:,:,pT_4_v), &
               & FragEnergies_tmp,ccenergies(pT_4_o))

            call ccsdpt_energy_e5_full(nocc,nvirt,nfrags,ncore,t1_final,ccsdpt_t1,&
               & occ_orbitals,virt_orbitals,FragEnergies(:,:,pT_5_o),FragEnergies(:,:,pT_5_v), &
               & ccenergies(pT_5_o))
         end if
   
         ! calculate total (T) contributions:
         ccenergies(pT_full_o) = ccenergies(pT_4_o)+ccenergies(pT_5_o)
         FragEnergies(:,:,pT_full_o:pT_full_v) = FragEnergies(:,:,pT_4_o:pT_4_v) &
            & + FragEnergies(:,:,pT_5_o:pT_5_v)

      endif
   
      ! Print all fragment energies
      call print_fragment_energies_full(nfrags,FragEnergies,ccenergies,orbitals_assigned,&
         & mymolecule%DistanceTable)
   
      do i=1,nocc_tot
         call orbital_free(occ_orbitals(i))
      end do
   
      call mem_dealloc(occ_orbitals)
   
      do i=1,nvirt
         call orbital_free(virt_orbitals(i))
      end do
      call mem_dealloc(virt_orbitals)
      call mem_dealloc(orbitals_assigned)
   
      ! release stuff
      call mem_dealloc(FragEnergies)
      call mem_dealloc(FragEnergies_tmp)

      call tensor_free(t2f_local)
      call tensor_free(VOVO_local)
   
   else ! we do not print fragment energies

      if (.not. DECinfo%pt_hack) then

         call ccsolver_job(ccmodel,Co,Cv,fock,nbasis,nocc,nvirt,mylsitem,ccPrintLevel, &
            & oof,vvf,ccenergies(cc_sol_o),VOVO,.false.,local,t1_final,t2_final, &
            & m1f=m1_final, m2f=m2_final)

         if(DECinfo%CCSDmultipliers)then
            call tensor_free(m2_final)
            if(DECinfo%use_singles)then
               call tensor_free(m1_final)
            endif
         endif
   
         if(DECinfo%PL>1)then
            call time_start_phase(PHASE_WORK, dwwork = time_CCSD_work, dwcomm = time_CCSD_comm, dwidle = time_CCSD_idle, &
               &swwork = time_pT_work, swcomm = time_pT_comm, swidle = time_pT_idle, &
               &labeldwwork = 'MASTER WORK CC solver: ',&
               &labeldwcomm = 'MASTER COMM CC solver: ',&
               &labeldwidle = 'MASTER IDLE CC solver: ')
         endif
   
         ! If there are two subsystems, calculate dispersion, charge transfer and subsystem energy contributions 
         if(mylsitem%input%molecule%nSubSystems==2) then
            !THIS IS JUST A WORKAROUND, ccsolver_par gives PDM tensors if more than
            !one node is used  FIXME
            if(VOVO%itype==TT_DENSE .and. t2_final%itype==TT_DENSE) then
               call SNOOP_partition_energy(VOVO,t1_final,t2_final,mylsitem,MyMolecule)
            else
               call tensor_init( VOVO_local, VOVO%dims,    4, bg=bg )
               call tensor_init( t2f_local, t2_final%dims, 4, bg=bg )
               call tensor_cp_data( VOVO,     VOVO_local )
               call tensor_cp_data( t2_final, t2f_local  )
               call SNOOP_partition_energy(VOVO_local,t1_final,t2f_local,mylsitem,MyMolecule)
               call tensor_free(t2f_local)
               call tensor_free(VOVO_local)
            end if
         end if

      else

         call tensor_init(t1_final,[nvirt,nocc],2)
         call tensor_random(t1_final)
         call tensor_scale(t1_final,1.0E-1_realk)

      endif
 
      if(ccmodel == MODEL_CCSDpT)then

         if (abc) then
            call tensor_init(ccsdpt_t1,[nocc,nvirt],2,bg=bg)
         else
            call tensor_init(ccsdpt_t1,[nvirt,nocc],2,bg=bg)
         endif

         call ccsdpt_driver(nocc,nvirt,nbasis,oof,vvf,Co,Cv,mylsitem,VOVO,t2_final,ccsdpt_t1,print_frags, &
            & abc,e4=ccenergies(pT_4_o))

         if (abc) call tensor_reorder(ccsdpt_t1,[2,1])
         call ccsdpt_energy_e5_ddot(nocc,nvirt,ccsdpt_t1%elm1,t1_final%elm1,ccenergies(pT_5_o))

         ! sum up energies
         ccenergies(pT_full_o) = ccenergies(pT_4_o) + ccenergies(pT_5_o)

         if(DECinfo%PL>1)then
            call time_start_phase(PHASE_WORK,dwwork = time_pT_work, dwcomm = time_pT_comm, dwidle = time_pT_idle, &
               &labeldwwork = 'MASTER WORK pT: ',&
               &labeldwcomm = 'MASTER COMM pT: ',&
               &labeldwidle = 'MASTER IDLE pT: ')
         endif

!         ! free integrals
!         call tensor_free(vovo)

      endif

      CorrEnergyString = 'correlation energy            '
      iCorrLen = 18
      write(DECinfo%output,*)
      write(DECinfo%output,*)
      write(DECinfo%output,*)
      if(ccmodel == MODEL_RPA)then
         write(DECinfo%output,'(1X,a,a,a,g20.10)') 'RPA ', &
            & CorrEnergyString(1:iCorrLen),' : ', ccenergies(cc_sol_o)
      elseif(ccmodel == MODEL_SOSEX)then
         write(DECinfo%output,'(1X,a,a,a,g20.10)') 'SOSEX ', &
            & CorrEnergyString(1:iCorrLen),' : ', ccenergies(cc_sol_o)
      else if(ccmodel == MODEL_MP2)then
         write(DECinfo%output,'(1X,a,a,a,g20.10)') 'MP2 ', &
            & CorrEnergyString(1:iCorrLen),' : ', ccenergies(cc_sol_o)
      else if(ccmodel == MODEL_CC2)then
         write(DECinfo%output,'(1X,a,a,a,g20.10)') 'CC2 ', &
            & CorrEnergyString(1:iCorrLen),' : ', ccenergies(cc_sol_o)
      else if(ccmodel == MODEL_CCSD)then
         write(DECinfo%output,'(1X,a,a,a,g20.10)') 'CCSD ', &
            & CorrEnergyString(1:iCorrLen),' : ', ccenergies(cc_sol_o)
      else if(ccmodel == MODEL_CCSDpT)then
         write(DECinfo%output,'(1X,a,a,a,g20.10)') 'CCSD ', &
            & CorrEnergyString(1:iCorrLen),' : ', ccenergies(cc_sol_o)
         write(DECinfo%output,'(1X,a,g20.10)') '(T) correlation energy  : ', &
            & ccenergies(pT_full_o)
         write(DECinfo%output,'(1X,a,g20.10)') '(T) 4th order energy    : ', &
            & ccenergies(pT_4_o)
         write(DECinfo%output,'(1X,a,g20.10)') '(T) 5th order energy    : ', &
            & ccenergies(pT_5_o)
         write(DECinfo%output,*)
         write(DECinfo%output,'(1X,a,a,a,g20.10)') 'Total CCSD(T) ', &
            & CorrEnergyString(1:iCorrLen),' : ', ccenergies(cc_sol_o)+ccenergies(pT_full_o)
      end if
      write(DECinfo%output,*)

   endif

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
   else if (ccmodel == MODEL_RIMP2) then
      write(DECinfo%output,'(1X,a)') '*                      Full RI-MP2 calculation is done !                    *'
   else if (ccmodel == MODEL_LSTHCRIMP2) then
      write(DECinfo%output,'(1X,a)') '*                      Full LS-THC-RI-MP2 calculation is done !             *'
   else if (ccmodel == MODEL_RPA ) then
      write(DECinfo%output,'(1X,a)') '*                      Full dRPA calculation is done !                       *'
   else if (ccmodel == MODEL_SOSEX ) then
      write(DECinfo%output,'(1X,a)') '*                      Full SOSEX calculation is done !                       *'
   else
      call lsquit("ERROR(ccsolver_justenergy)model not recognized",-1)
   endif
   write(DECinfo%output,'(1X,a)') '*****************************************************************************'
   write(DECinfo%output,*)


   ! now update ccenergy with ccsd(t) correction
   ccenergy = ccenergies(cc_sol_o) + ccenergies(pT_full_o)


   if(ccmodel == MODEL_CCSDpT)then
      if (DECinfo%print_frags) then
         call tensor_free(ccsdpt_t2)
         call tensor_free(ccsdpt_t1)
      else
         call tensor_free(ccsdpt_t1)
      endif
   endif

   call mem_dealloc(ccenergies)
!else mod unreleased
#else

   call ccsolver_job(ccmodel,Co,Cv,fock,nbasis,nocc,nvirt,mylsitem,ccPrintLevel,oof,vvf,ccenergy,&
      & VOVO,.false.,local,t1_final,t2_final)
!endif mod unreleased
#endif

   if( ccmodel /= MODEL_MP2 .and. ccmodel /= MODEL_RPA &
     &.and. ccmodel /= MODEL_SOSEX ) then
      ! free amplitude arrays
      call tensor_free(t1_final)
   endif
   if (.not. DECinfo%pt_hack) then
      call tensor_free(t2_final)
      call tensor_free(VOVO)
   endif

end function ccsolver_justenergy

!> \brief For a given fragment, calculate singles and doubles amplitudes and
!> two-electron integrals (a i | bj ) required for CC energy.
!> Intended to be used for CC2 and CCSD (and NOT for MP2).
!> \author Kasper Kristensen, heavily modifed by PE
!> \date January 2012
subroutine fragment_ccsolver(MyFragment,t1,t2,VOVO,m1,m2)

   implicit none

   !> Fragment info (only t1 information in MyFragment may be changed here)
   type(decfrag), intent(inout) :: MyFragment
   !> Singles amplitudes t1(a,i)
   type(tensor),intent(inout) :: t1
   !> Doubles amplitudes t2(a,i,b,j)
   type(tensor),intent(inout) :: t2
   !> Two electron integrals (a i | b j) stored as (a,i,b,j)
   type(tensor),intent(inout) :: VOVO
   !> Singles multipliers m1(a,i)
   type(tensor),intent(inout), optional :: m1
   !> Doubles multipliers m2(a,i,b,j)
   type(tensor),intent(inout), optional :: m2

   !INTERNAL PARAMETERS
   type(tensor) :: mp2_amp
   integer :: dims(2)
   real(realk) :: ccenergy
   logical :: local

   ! Sanity check: This routine is not intended for MP2 (except for DECNP)
   if(MyFragment%ccmodel == MODEL_MP2 .and. (.not.DECinfo%DECNP)) then
      call lsquit('fragment_ccsolver cannot be used for MP2!',&
      & DECinfo%output)
   end if

   local=.true.
#ifdef VAR_MPI
   if(infpar%lg_nodtot>1) local = .false.
#endif


   ! If MyFragment%t1_stored is TRUE, then we reuse the singles amplitudes
   ! from previous fragment calculations to describe long-range
   ! singles effects.
   ! In this case the fragment t1 amplitudes are stored in MyFragment%t1
   if(MyFragment%t1_stored) then
      dims(1) = MyFragment%nvirtAOS
      dims(2) = MyFragment%noccAOS
      call tensor_init(t1,dims,2)
      call tensor_convert(MyFragment%t1,t1)
   end if

   call ccsolver_job(myfragment%ccmodel,myfragment%Co,myfragment%Cv,myfragment%fock,myfragment%nbasis,myfragment%noccAOS,&
      & myfragment%nvirtAOS,myfragment%mylsitem,DECinfo%PL,myfragment%ppfock,myfragment%qqfock,ccenergy,&
      & VOVO,myfragment%t1_stored,local,t1,t2,m1f=m1,m2f=m2,frag=MyFragment)


   ! Save singles amplitudes in fragment structure
   if(DECinfo%SinglesPolari) then
      call save_fragment_t1_AOSAOSamplitudes(MyFragment,t1%elm2)
   end if


end subroutine fragment_ccsolver



!> \brief Solve MP2 equation:
!> RHS_{bjai} =  - sum_{c} t_{bjci} F_{ca}
!>               - sum_{c} t_{aicj} F_{cb}
!>               + sum_{k} t_{bjak} F_{ki}
!>               + sum_{k} t_{aibk} F_{kj}
!> It is assumed that RHS_{bjai} = R_{aibj} !
!> \author Patrick Ettenhuber adapted from Kasper Kristensen
!> \date October 2014
subroutine mp2_solver_frag(frag,RHS,t2,rhs_input,mp2_energy)

   implicit none
   !> fragment
   type(decfrag), intent(inout) :: frag 
   !> RHS array
   type(tensor), intent(inout)  :: RHS
   !> Solution array
   type(tensor), intent(inout)  :: t2
   !> logical to specify whether to use the supplied RHS or calculate it in the
   !solver
   logical, intent(in)          :: rhs_input
   !> output mp2 energy
   real(realk), intent(out), optional   :: mp2_energy
   real(realk)  :: tcpu1,twall1,tcpu2,twall2
   real(realk) :: e
   integer :: Ncore, no, nv, nb
   logical :: local

   call LSTIMER('START',tcpu1,twall1,DECinfo%output)

   local = .true.
#ifdef VAR_MPI
   local = (infpar%lg_nodtot==1)
#endif

   !if(DECinfo%array4OnFile) then ! RHS and t2 values are stored on file
   !   call mp2_solver_file(nocc,nvirt,ppfock,qqfock,RHS,t2)
   !else ! RHS and t2 values are stored in memory
   !   call mp2_solver_mem(nocc,nvirt,ppfock,qqfock,RHS,t2)
   !end if

   nb    = frag%nbasis
   nv    = frag%nvirtAOS
   no    = frag%noccAOS
   Ncore = frag%ncore

   call ccsolver(MODEL_MP2,frag%Co,frag%Cv,frag%fock,nb,&
      & no,nv, frag%mylsitem,DECinfo%PL,e,RHS,.false.,local,SOLVE_AMPLITUDES,frag=frag,&
      & vovo_supplied = rhs_input, p4=t2)

   if(present(mp2_energy)) mp2_energy = e

   call LSTIMER('START',tcpu2,twall2,DECinfo%output)

end subroutine mp2_solver_frag

!> \brief Solve MP2 equation:
!> RHS_{bjai} =  - sum_{c} t_{bjci} F_{ca}
!>               - sum_{c} t_{aicj} F_{cb}
!>               + sum_{k} t_{bjak} F_{ki}
!>               + sum_{k} t_{aibk} F_{kj}
!> It is assumed that RHS_{bjai} = R_{aibj} !
!> \author Patrick Ettenhuber adapted from Kasper Kristensen
!> \date October 2014
subroutine mp2_solver_mol(mol,mls,RHS,t2,rhs_input,mp2_energy)

   implicit none
   !> Full molecular information
   type(fullmolecule), intent(in) :: mol 
   !corresponding lsitem info
   type(lsitem), intent(inout)  :: mls
   !> RHS array
   type(tensor), intent(inout)  :: RHS
   !> Solution array
   type(tensor), intent(inout)  :: t2
   !> logical to specify whether to use the supplied RHS or calculate it in the
   !solver
   logical, intent(in)          :: rhs_input
   !> output mp2 energy
   real(realk), intent(out), optional   :: mp2_energy

   !internal variables
   real(realk)  :: tcpu1,twall1,tcpu2,twall2
   real(realk)  :: e
   real(realk), pointer :: Co(:,:),oof(:,:)
   integer :: Ncore, no, nv, nb, i, j
   logical :: local

   call LSTIMER('START',tcpu1,twall1,DECinfo%output)

   if( mol%mem_distributed )then
      call lsquit("ERROR(mp2_solver_mol) not implemented for distributed mol",-1)
   endif

   local = .true.
#ifdef VAR_MPI
   local = (infpar%lg_nodtot==1)
#endif

   !if(DECinfo%array4OnFile) then ! RHS and t2 values are stored on file
   !   call mp2_solver_file(nocc,nvirt,ppfock,qqfock,RHS,t2)
   !else ! RHS and t2 values are stored in memory
   !   call mp2_solver_mem(nocc,nvirt,ppfock,qqfock,RHS,t2)
   !end if

   nb    = mol%nbasis
   nv    = mol%nvirt
   Ncore = mol%ncore

   if(DECinfo%frozencore) then

      no = mol%nval

      ! Only copy valence orbitals

      call mem_alloc(Co,nb,no)
      call mem_alloc(oof,no,no)

      do i=1,no
         Co(:,i) = mol%Co%elm2(:,i+Ncore)
      end do

      ! Fock valence
      do j=1,no
         do i=1,no
            oof(i,j) = mol%oofock%elm2(i+Ncore,j+Ncore)
         end do
      end do

   else
      ! No frozen core, simply copy elements for all occupied orbitals
      no = mol%nocc

      call mem_alloc(Co,nb,no)
      call mem_alloc(oof,no,no)

      Co  = mol%Co%elm2
      oof = mol%oofock%elm2

   end if

   call ccsolver(MODEL_MP2,Co,mol%Cv%elm2,mol%fock%elm2,nb,no,nv,mls,DECinfo%PL,&
      & e,RHS,.false.,local,SOLVE_AMPLITUDES, vovo_supplied=rhs_input, p4=t2)

   if(present(mp2_energy)) mp2_energy = e

   call mem_dealloc(Co)
   call mem_dealloc(oof)

   call LSTIMER('START',tcpu2,twall2,DECinfo%output)
end subroutine mp2_solver_mol



!!> \brief Solve MP2 equation when RHS and t2 are values are stored in memory.
!!> See mp2_solver for details about the equation.
!!> \author Kasper Kristensen
!!> \date February 2011
!subroutine mp2_solver_mem(nocc,nvirt,ppfock,qqfock,RHS,t2)
!
!   implicit none
!   !> Number of occupied orbitals (dimension of ppfock)
!   integer, intent(in) :: nocc
!   !> Number of virtupied orbitals (dimension of qqfock)
!   integer, intent(in) :: nvirt
!   !> Occupied-occupied block of Fock matrix
!   real(realk) :: ppfock(nocc,nocc)
!   !> virtupied-virtupied block of Fock matrix
!   real(realk) :: qqfock(nvirt,nvirt)
!   !> RHS array
!   type(array4), intent(in) :: RHS
!   !> Solution array
!   type(array4), intent(inout) :: t2
!   type(array4) :: tmp1,tmp2
!   real(realk),pointer :: Cocc_data(:,:), Cvirt_data(:,:), Socc(:,:), Svirt(:,:)
!   real(realk),pointer :: EVocc(:), EVvirt(:)
!   type(array2) :: Cocc, Cvirt
!   integer :: I,J,A,B
!   real(realk) :: tcpu, twall, deltaF
!   integer :: dims(4), occdims(2), virtdims(2)
!   ! real(realk) :: test
!
!
!   ! Strategy for solving MP2 equation:
!   ! 1. Find basis where Fock matrix is diagonal
!   ! 2. Transform 2-electron integrals to diagonal basis
!   ! 3. In diagonal basis the solution is trivial and the amplitudes are found.
!   ! 4. Transform amplitudes in diagonal basis back to LCM basis.
!
!   call LSTIMER('START',tcpu,twall,DECinfo%output)
!
!
!   write(DECinfo%output,*)
!   write(DECinfo%output,*) 'Entering MP2 solver - store array values in memory'
!   write(DECinfo%output,*)
!
!
!   ! Sanity checks
!   ! *************
!   ! Check that nvirt /= 0.
!   if(nvirt<1 .or. nocc<1) then
!      write(DECinfo%output,*) 'Number of occupied orbitals = ', nocc
!      write(DECinfo%output,*) 'Number of virtupied orbitals = ', nvirt
!      call lsquit('Error in mp2_solver: Number of orbitals is smaller than one!', DECinfo%output)
!   endif
!
!
!
!
!   ! Initialize stuff
!   ! ****************
!   dims = [nvirt,nocc,nvirt,nocc]
!   occdims = [nocc,nocc]
!   virtdims = [nvirt,nvirt]
!
!
!
!
!   ! 1. Solve Fock eigenvalue problem - each block separately
!   ! ********************************************************
!
!   ! OCCUPIED-OCCUPIED BLOCK
!   ! '''''''''''''''''''''''
!
!   ! Eigenvectors
!   call mem_alloc(Cocc_data,nocc,nocc)
!
!   ! Eigenvalues
!   call mem_alloc(EVocc,nocc)
!
!   ! The overlap matrix is simply the unit matrix because
!   ! the LCM/MLM basis is orthogonal.
!   call mem_alloc(Socc,nocc,nocc)
!   Socc=0.0E0_realk
!   do i=1,nocc
!      Socc(i,i) = 1E0_realk
!   end do
!
!   ! Solve eigenvalue problem
!   call solve_eigenvalue_problem(nocc,ppfock,Socc,EVocc,Cocc_data)
!
!   ! For later, it is convenient to keep eigenvectors in array2 form
!   Cocc = array2_init(occdims,Cocc_data)
!
!   ! Done with some matrices
!   call mem_dealloc(Cocc_data)
!   call mem_dealloc(Socc)
!
!
!
!
!   ! VIRTUAL-VIRTUAL BLOCK
!   ! '''''''''''''''''''''
!
!   ! Eigenvectors
!   call mem_alloc(Cvirt_data,nvirt,nvirt)
!
!   ! Eigenvalues
!   call mem_alloc(EVvirt,nvirt)
!
!   ! Unit overlap for virtual space
!   call mem_alloc(Svirt,nvirt,nvirt)
!   Svirt=0.0E0_realk
!   do i=1,nvirt
!      Svirt(i,i) = 1E0_realk
!   end do
!
!   ! Solve eigenvalue problem
!   call solve_eigenvalue_problem(nvirt,qqfock,Svirt,EVvirt,Cvirt_data)
!
!   ! For later, it is convenient to keep eigenvectors in array2 form
!   Cvirt = array2_init(virtdims,Cvirt_data)
!
!   ! Done with some matrices
!   call mem_dealloc(Cvirt_data)
!   call mem_dealloc(Svirt)
!
!
!   call LSTIMER('SOLVE: EIVAL',tcpu,twall,DECinfo%output)
!
!
!
!   ! 2. Transform two-electron integrals to diagonal basis
!   ! *****************************************************
!
!   ! Using notation that (a,i,b,j) are LCM indices
!   ! and (A,I,B,J) are indices in the diagonal basis,
!   ! we want to carry out the transformations:
!   ! RHS_{AIBJ} = sum_{aibj} C_{aA} C_{iI} C_{bB} C_{jJ} RHS_{aibj} (*)
!
!   ! 1. Init temporary array
!   tmp1= array4_init(dims)
!
!   ! 2. Index A: RHS(a,i,b,j) --> tmp1(A,i,b,j)
!   call array4_contract1(RHS,Cvirt,tmp1,.true.)
!
!   ! 3. Index I: tmp1(A,i,b,j) --> tmp1(i,A,b,j) --> tmp2(I,A,b,j)
!   call array4_reorder(tmp1,[2,1,3,4])
!   tmp2= array4_init([nocc,nvirt,nvirt,nocc])
!
!   call array4_contract1(tmp1,Cocc,tmp2,.true.)
!   call array4_free(tmp1)
!
!   ! 4. Index J: tmp2(I,A,b,j) --> tmp2(j,b,A,I) --> tmp1(J,b,A,I)
!   call array4_reorder(tmp2,[4,3,2,1])
!   tmp1= array4_init([nocc,nvirt,nvirt,nocc])
!   call array4_contract1(tmp2,Cocc,tmp1,.true.)
!   call array4_free(tmp2)
!
!   ! 5. Index B: tmp1(J,b,A,I) --> tmp1(b,J,A,I) --> t2(B,J,A,I)
!   call array4_reorder(tmp1,[2,1,3,4])
!   t2 = array4_init(dims)
!   call array4_contract1(tmp1,Cvirt,t2,.true.)
!   call array4_free(tmp1)
!
!   ! Now t2 contains the two-electron integrals (AI|BJ) = (BJ|AI) in the order (B,J,A,I).
!   ! [Due to the symmetry it is not necessary to reorder t2 back to (A,I,B,J)].
!
!
!   call LSTIMER('SOLVE: TRANS 1',tcpu,twall,DECinfo%output)
!
!
!
!   ! 3. Solve MP2 equation in the diagonal basis
!   ! *******************************************
!
!   ! In the diagonal basis the solution to the MP2 equation is trivial!
!   ! The equation:
!   ! RHS_{bjai} =  - sum_{c} t_{bjci} F_{ca}
!   !               - sum_{c} t_{aicj} F_{cb}
!   !               + sum_{k} t_{bjak} F_{ki}
!   !               + sum_{k} t_{aibk} F_{kj}
!   !
!   ! simply becomes (using t_{BJAI} = t_{AIBJ}):
!   !
!   ! RHS_{BJAI} =  - t_{BJAI} F_{AA}
!   !               - t_{AIBJ} F_{BB}
!   !               + t_{BJAI} F_{II}
!   !               + t_{AIBJ} F_{JJ}
!   !            =  t_{BJAI} [ F_{II} + F_{JJ} - F_{AA} - F_{BB} ]
!   !
!   ! In other words:
!   !
!   ! t_{BJAI} = RHS_{BJAI} / [ F_{II} + F_{JJ} - F_{AA} - F_{BB} ] (**)
!   !
!
!   ! Recalling that currently t_{BJAI} = RHS_{BJAI} and that the
!   ! Fock matrix in the diagonal basis are the eigenvalues EVocc and EVvirt -
!   ! we simply modify the t2 elements as follows:
!
!   !    test=0E0_realk
!   do I=1,nocc
!      do A=1,nvirt
!         do J=1,nocc
!            do B=1,nvirt
!
!               ! deltaF = F_{II} + F_{JJ} - F_{AA} - F_{BB}
!               deltaF = EVocc(I) + EVocc(J) - EVvirt(A) - EVvirt(B)
!
!               ! Sanity check
!               if( abs(deltaF) < 1e-9_realk ) then
!                  write(DECinfo%output,*) 'WARNING: SMALL NUMBERS OCCURING IN SOLVER!!!'
!                  write(DECinfo%output,*) 'WARNING: SOLVER MAY BE UNSTABLE!!!'
!                  write(DECinfo%output,*) 'Delta epsilon value = ', deltaF
!               end if
!
!               ! Canonical energy check
!               ! test = test &
!               ! & + (2E0_realk*t2%val(B,J,A,I) - t2%val(B,I,A,J))*t2%val(B,J,A,I)/deltaF
!
!               ! t2 according to (**)
!               t2%val(B,J,A,I) = t2%val(B,J,A,I)/deltaF
!
!            end do
!         end do
!      end do
!   end do
!
!   call LSTIMER('SOLVE: CALC T2',tcpu,twall,DECinfo%output)
!
!
!
!
!   ! 4. Transform amplitudes back to LCM/MLM basis
!   ! *********************************************
!
!   ! Since LCM/MLM and the diagonal basis are connected by a unitary
!   ! transformation this basically corresponds to repeating step 2
!   ! above with the transposed transformation matrices:
!   !
!   ! t_{aibj} = sum_{AIBJ} C_{Aa} C_{Ia} C_{Bb} C_{Jj} t_{AIBJ}
!
!   ! 1. Transpose transformation matrices
!   call array2_transpose(Cocc)
!   call array2_transpose(Cvirt)
!
!   ! 2. Index A: t(A,I,B,J) --> tmp1(a,I,B,J)
!   tmp1= array4_init(dims)
!   call array4_contract1(t2,Cvirt,tmp1,.true.)
!   call array4_free(t2)
!
!   ! 3. Index I: tmp1(a,I,B,J) --> tmp1(I,a,B,J) --> tmp2(i,a,B,J)
!   call array4_reorder(tmp1,[2,1,3,4])
!   tmp2= array4_init([nocc,nvirt,nvirt,nocc])
!   call array4_contract1(tmp1,Cocc,tmp2,.true.)
!   call array4_free(tmp1)
!
!   ! 4. Index J: tmp2(i,a,B,J) --> tmp2(J,B,a,i) --> tmp1(j,B,a,i)
!   call array4_reorder(tmp2,[4,3,2,1])
!   tmp1 = array4_init([nocc,nvirt,nvirt,nocc])
!   call array4_contract1(tmp2,Cocc,tmp1,.true.)
!   call array4_free(tmp2)
!
!   ! 5. Index B: tmp1(j,B,a,i) --> tmp1(B,j,a,i) --> t2(b,j,a,i)
!   call array4_reorder(tmp1,[2,1,3,4])
!   t2 = array4_init(dims)
!   call array4_contract1(tmp1,Cvirt,t2,.true.)
!   call array4_free(tmp1)
!
!   call LSTIMER('SOLVE: TRANS 2',tcpu,twall,DECinfo%output)
!
!
!   ! Clean up
!   ! ********
!   call mem_dealloc(EVocc)
!   call mem_dealloc(EVvirt)
!   call array2_free(Cocc)
!   call array2_free(Cvirt)
!
!
!
!end subroutine mp2_solver_mem
!
!
!
!!> \brief Solve MP2 equation when RHS and t2 are values are stored on file.
!!> See mp2_solver for details about the equation.
!!> \author Kasper Kristensen
!!> \date February 2011
!subroutine mp2_solver_file(nocc,nvirt,ppfock,qqfock,RHS,t2)
!
!   implicit none
!   !> Number of occupied orbitals (dimension of ppfock)
!   integer, intent(in) :: nocc
!   !> Number of virtupied orbitals (dimension of qqfock)
!   integer, intent(in) :: nvirt
!   !> Occupied-occupied block of Fock matrix
!   real(realk) :: ppfock(nocc,nocc)
!   !> virtupied-virtupied block of Fock matrix
!   real(realk) :: qqfock(nvirt,nvirt)
!   !> RHS array, storing type 1 (see array4_init_file)
!   type(array4), intent(in) :: RHS
!   !> Solution array, storing type 1 (see array4_init_file)
!   type(array4), intent(inout) :: t2
!   type(array4) :: tmp1,tmp2
!   type(array4) :: RHSaib,RHSbaj,RHStmp,t2tmp
!   real(realk),pointer :: Cocc_data(:,:), Cvirt_data(:,:), Socc(:,:), Svirt(:,:)
!   real(realk),pointer :: EVocc(:), EVvirt(:)
!   type(array2) :: Cocc, Cvirt,CoccT,CvirtT
!   integer :: I,J,A,B
!   real(realk) :: tcpu, twall, deltaF
!   integer :: occdims(2), virtdims(2)
!
!
!
!   ! Strategy for solving MP2 equation:
!   ! 1. Find basis where Fock matrix is diagonal
!   ! 2. Transform 2-electron integrals to diagonal basis
!   ! 3. In diagonal basis the solution is trivial and the amplitudes are found.
!   ! 4. Transform amplitudes in diagonal basis back to LCM basis.
!
!   ! Steps 2-4 necessarily overlap when we store array values on file.
!   call LSTIMER('START',tcpu,twall,DECinfo%output)
!
!
!   write(DECinfo%output,*)
!   write(DECinfo%output,*) 'Entering MP2 solver - store array values on file'
!   write(DECinfo%output,*)
!
!
!   ! Sanity checks
!   ! *************
!   ! Check that nvirt /= 0.
!   if(nvirt<1 .or. nocc<1) then
!      write(DECinfo%output,*) 'Number of occupied orbitals = ', nocc
!      write(DECinfo%output,*) 'Number of virtupied orbitals = ', nvirt
!      call lsquit('Error in mp2_solver: Number of orbitals is smaller than one!', DECinfo%output)
!   endif
!
!
!   ! Initialize stuff
!   ! ****************
!   occdims = [nocc,nocc]
!   virtdims = [nvirt,nvirt]
!
!
!
!   ! 1. Solve Fock eigenvalue problem - each block separately
!   ! ********************************************************
!
!   ! OCCUPIED-OCCUPIED BLOCK
!   ! '''''''''''''''''''''''
!
!   ! Eigenvectors
!   call mem_alloc(Cocc_data,nocc,nocc)
!   ! Eigenvalues
!   call mem_alloc(EVocc,nocc)
!
!   ! The overlap matrix is simply the unit matrix because
!   ! the LCM/MLM basis is orthogonal.
!   call mem_alloc(Socc,nocc,nocc)
!   Socc=0.0E0_realk
!   do i=1,nocc
!      Socc(i,i) = 1E0_realk
!   end do
!
!   ! Solve eigenvalue problem
!   call solve_eigenvalue_problem(nocc,ppfock,Socc,EVocc,Cocc_data)
!
!   ! For later, it is convenient to keep eigenvectors in array2 form
!   ! and also to have the transposed matrices available
!   Cocc = array2_init(occdims,Cocc_data)
!   CoccT = array2_init(occdims)
!   call array2_copy(CoccT,Cocc)
!   call array2_transpose(CoccT)
!
!   ! Done with some matrices
!   call mem_dealloc(Cocc_data)
!   call mem_dealloc(Socc)
!
!
!
!
!   ! VIRTUAL-VIRTUAL BLOCK
!   ! '''''''''''''''''''''
!
!   ! Eigenvectors
!   call mem_alloc(Cvirt_data,nvirt,nvirt)
!
!   ! Eigenvalues
!   call mem_alloc(EVvirt,nvirt)
!
!   ! Unit overlap for virtual space
!   call mem_alloc(Svirt,nvirt,nvirt)
!   Svirt=0.0E0_realk
!   do i=1,nvirt
!      Svirt(i,i) = 1E0_realk
!   end do
!
!   ! Solve eigenvalue problem
!   call solve_eigenvalue_problem(nvirt,qqfock,Svirt,EVvirt,Cvirt_data)
!
!   ! For later, it is convenient to keep eigenvectors in array2 form
!   ! and also to have the transposed matrices available
!   Cvirt = array2_init(virtdims,Cvirt_data)
!   CvirtT = array2_init(virtdims)
!   call array2_copy(CvirtT,Cvirt)
!   call array2_transpose(CvirtT)
!
!
!   ! Done with some matrices
!   call mem_dealloc(Cvirt_data)
!   call mem_dealloc(Svirt)
!   call LSTIMER('SOLVE: EIVAL',tcpu,twall,DECinfo%output)
!
!
!
!   ! Transform three indices to diagonal basis (aib-->AIB)
!   ! *****************************************************
!
!   ! Using the notation that (a,i,b,j) are LCM indices
!   ! and (A,I,B,J) are indices in the diagonal basis,
!   ! we want to carry out the transformations:
!   ! RHS_{AIBJ} = sum_{aibj} C_{aA} C_{iI} C_{bB} C_{jJ} RHS_{aibj}
!
!   ! Since we cannot have full four-dimensional arrays in memory,
!   ! we do this in steps. First, for a fixed "j" we transform the other indices:
!   !
!   ! RHS_{AIBj} = sum_{aib} C_{aA} C_{iI} C_{bB} RHS_{aibj}
!
!   ! Temporary RHS where three indices are transformed, stored on file.
!   RHStmp = array4_init_file([nvirt,nvirt,nocc,nocc],2,.false.)
!   ! Temporary array for keeping RHS_{aibj} for a fixed j (last dimension is one).
!   RHSaib = array4_init([nvirt,nocc,nvirt,1])
!   call array4_open_file(RHS)
!   call array4_open_file(RHStmp)
!
!
!   do j=1,nocc
!
!      ! 1. Read in RHS_{aibj} for fixed j
!      call array4_read_file_type1(RHS,j,&
!         & RHSaib%val(1:nvirt,1:nocc,1:nvirt,1),nvirt,nocc,nvirt)
!
!      ! 2. Index A: RHS(a,i,b,j) --> tmp1(A,i,b,j)
!      tmp1= array4_init([nvirt,nocc,nvirt,1])
!      call array4_contract1(RHSaib,Cvirt,tmp1,.true.)
!
!      ! 3. Index I: tmp1(A,i,b,j) --> tmp1(i,A,b,j) --> tmp2(I,A,b,j)
!      call array4_reorder(tmp1,[2,1,3,4])
!      tmp2= array4_init([nocc,nvirt,nvirt,1])
!      call array4_contract1(tmp1,Cocc,tmp2,.true.)
!      call array4_free(tmp1)
!
!      ! 4. Index B: tmp2(I,A,b,j) --> tmp2(b,A,I,j) --> tmp1(B,A,I,j)
!      call array4_reorder(tmp2,[3,2,1,4])
!      tmp1= array4_init([nvirt,nvirt,nocc,1])
!      call array4_contract1(tmp2,Cvirt,tmp1,.true.)
!      call array4_free(tmp2)
!
!      ! 5. Write to file referenced by temporary RHS array (storing type 2)
!      do I=1,nocc
!         call array4_write_file_type2(RHStmp,I,j,tmp1%val(1:nvirt,1:nvirt,I,1),nvirt,nvirt)
!      end do
!      call array4_free(tmp1)
!
!   end do
!   call array4_close_file(RHS,'KEEP')
!   call array4_free(RHSaib)
!   call LSTIMER('SOLVE: STEP 1',tcpu,twall,DECinfo%output)
!   ! Now the file referenced by RHStmp contains RHS_{AIBj},
!   ! stored in the order (B,A,I,j) using storing type 2.
!
!
!
!   ! Transform j-->J; Solve equation in diag basis; Back transform ABJ-->abj
!   ! ***********************************************************************
!
!   ! Temporary array
!   t2tmp = array4_init_file([nvirt,nvirt,nocc,nocc],2,.false.)
!   call array4_open_file(t2tmp)
!
!
!   I_loop: do I=1,nocc
!
!
!      ! Read in RHS_{AIBj} for fixed I
!      ! ------------------------------
!      RHSbaj=array4_init([nvirt,nvirt,nocc,1])
!      do j=1,nocc
!         call array4_read_file_type2(RHStmp,I,j,&
!            & RHSbaj%val(1:nvirt,1:nvirt,j,1), nvirt,nvirt)
!      end do
!      ! Now RHSbaj contains elements RHS_{AIBj} for a fixed I
!      ! stored in the order (B,A,j,I)
!
!
!      ! Transform final index j-->J
!      ! ---------------------------
!      ! Reorder: RHSbaj(B,A,j,I) --> RHSbaj(j,B,A,I)
!      call array4_reorder(RHSbaj,[3,1,2,4])
!
!
!      ! Transform: RHSbaj(j,B,A,I) --> tmp1(J,B,A,I)
!      tmp1=array4_init([nocc,nvirt,nvirt,1])
!      call array4_contract1(RHSbaj,Cocc,tmp1,.true.)
!      call array4_free(RHSbaj)
!
!
!      ! Now tmp1 contains the RHS_{BJAI} in the diagonal basis
!      ! stored in the order (J,B,A,I) [fixed I].
!
!
!      ! Divide by deltaF (solve equation in diagonal basis)
!      ! ---------------------------------------------------
!      ! In the diagonal basis the t2 solution vector is simply:
!      ! t2_{BJAI} = RHS_{BJAI} / [ F_{II} + F_{JJ} - F_{AA} - F_{BB} ] (**)
!      ! [See (**) in subroutine mp2_solver_mem]
!      do A=1,nvirt
!         do B=1,nvirt
!            do J=1,nocc
!
!               ! deltaF = F_{II} + F_{JJ} - F_{AA} - F_{BB}
!               deltaF = EVocc(I) + EVocc(J) - EVvirt(A) - EVvirt(B)
!
!               ! Sanity check
!               if( abs(deltaF) < 1e-9_realk ) then
!                  write(DECinfo%output,*) 'WARNING: SMALL NUMBERS OCCURING IN SOLVER!!!'
!                  write(DECinfo%output,*) 'WARNING: SOLVER MAY BE UNSTABLE!!!'
!                  write(DECinfo%output,*) 'Delta epsilon value = ', deltaF
!               end if
!
!               ! Amplitude element according to (**)
!               tmp1%val(J,B,A,1) = tmp1%val(J,B,A,1)/deltaF
!
!            end do
!         end do
!      end do
!
!      ! Now tmp1 contains the t2_{BJAI} solution vector in the diagonal basis
!      ! stored in the order (J,B,A,I) [fixed I], and we "just"
!      ! need to transform back to the original basis.
!
!
!      ! Transform back three indices: JBA --> jba
!      ! -----------------------------------------
!      ! Note: Use Transposed transformation matrices to back transform
!
!      ! 1. Index j: tmp1(J,B,A,I) -->  tmp2(j,B,A,I)
!      tmp2= array4_init([nocc,nvirt,nvirt,1])
!      call array4_contract1(tmp1,CoccT,tmp2,.true.)
!      call array4_free(tmp1)
!
!      ! 2. Index b: tmp2(j,B,A,I) --> tmp2(B,A,j,I) --> tmp1(b,A,j,I)
!      call array4_reorder(tmp2,[2,3,1,4])
!      tmp1= array4_init([nvirt,nvirt,nocc,1])
!      call array4_contract1(tmp2,CvirtT,tmp1,.true.)
!      call array4_free(tmp2)
!
!      ! 3. Index a: tmp1(b,A,j,I) --> tmp1(A,b,j,I) --> tmp2(a,b,j,I)
!      call array4_reorder(tmp1,[2,1,3,4])
!      tmp2= array4_init([nvirt,nvirt,nocc,1])
!      call array4_contract1(tmp1,CvirtT,tmp2,.true.)
!      call array4_free(tmp1)
!
!
!      ! Write solution vector t2_{aIbj} to file using storing type 2, order: (a,b,j,I)
!      ! ------------------------------------------------------------------------------
!
!      ! Note: Now we only need to transform the last "I" index of t2 back.
!      ! To do this we need to save on file such that we later can read in the full set of "I's"
!      do j=1,nocc
!         call array4_write_file_type2(t2tmp,j,I,&
!            & tmp2%val(1:nvirt,1:nvirt,j,1), nvirt,nvirt )
!      end do
!      call array4_free(tmp2)
!
!   end do I_loop
!
!   call array4_close_file(RHStmp,'DELETE')
!   call array4_free(RHStmp)
!   call LSTIMER('SOLVE: STEP 2',tcpu,twall,DECinfo%output)
!   ! Now the file assosicated with t2tmp contains the final amplitudes,
!   ! except that the I index must be transformed back.
!
!
!
!   ! Transform final t2 index (I-->i) and write solution vector to file
!   ! ******************************************************************
!
!   ! Final t2 solution vector
!   t2 = array4_init_file([nvirt,nocc,nvirt,nocc],1,.false.)
!   call array4_open_file(t2)
!
!   do j=1,nocc
!
!      tmp1 = array4_init([nvirt,nvirt,nocc,1])
!      do I=1,nocc
!         call array4_read_file_type2(t2tmp,j,I,tmp1%val(1:nvirt,1:nvirt,I,1),nvirt,nvirt)
!      end do
!      ! tmp1 now contains amplitudes t2_{aIbj} in the order (a,b,I,j) for fixed j.
!
!
!      ! Transform final index I-->i
!      ! ---------------------------
!
!      ! tmp1(a,b,I,j) --> tmp1(I,a,b,j) --> tmp2(i,a,b,j)
!      call array4_reorder(tmp1,[3,1,2,4])
!      tmp2 = array4_init([nocc,nvirt,nvirt,1])
!      call array4_contract1(tmp1,CoccT,tmp2,.true.)
!      call array4_free(tmp1)
!      ! Reorder to final storing order: tmp2(i,a,b,j) --> tmp2(a,i,b,j)
!      call array4_reorder(tmp2,[2,1,3,4])
!
!
!      ! Write t2_{aibj} to file for each j
!      ! ----------------------------------
!      call array4_write_file_type1(t2,j,&
!         & tmp2%val(1:nvirt,1:nocc,1:nvirt,1),nvirt,nocc,nvirt)
!      call array4_free(tmp2)
!
!   end do
!
!   call array4_close_file(t2tmp,'DELETE')
!   call array4_free(t2tmp)
!   call array4_close_file(t2,'KEEP')
!   call LSTIMER('SOLVE: STEP 3',tcpu,twall,DECinfo%output)
!
!
!   ! Clean up
!   ! ********
!   call mem_dealloc(EVocc)
!   call mem_dealloc(EVvirt)
!   call array2_free(Cocc)
!   call array2_free(Cvirt)
!   call array2_free(CoccT)
!   call array2_free(CvirtT)
!
!
!
!end subroutine mp2_solver_file


!> \brief adaption of the ccsolver routine, rebuild for the 
! use of parallel distributed memory. This solver is acutally a bit
! complicated in structure if used in an MPI framework. Most of the work
! happens hidden in the type(tensor) structure. It is highly recommended to
! begin implementing new features with setting local=.true. at the beginning
! and running without .spawn_comm_procs in the **CC input section. On
! INPUT:
! Co_f,Cv_f : the occupied and virtual orbital transformation coefficients
! fock_f      : the ao fock matrix
! nb,no,nv    : number of atomic, occupied and virtual orbitals respectively
! mylsitem    : the typical lsitem structure
! ccPrintLevel: print level (might be removed due to DECinfo%PL)
! fragment_job: specify whether it is a fragment job or a full calc
! longrange_singles : the longrange singles correction for the fock matrix on input
! local       : boolean which steers whether everything should be treated locally
! m1,m2 ...   : additional parameter sets needed for special calcultions (pno ccsd, multipliers)
!
! INPUT AND/OR OUTPUT:
! VOVO        : mo-integral WITHOUT T1 trafo on output
! vovo_supplied : does the VOVO contain the desired integrals or is it an empty containter, t
!
! OUTPUT:
! ccenergy    : output correlation energy
! p1,p2 ...   : parameter sets p with mode number containing the final parameters
!> \author Patrick Ettenhuber (heavily adapted version from Marcin)
subroutine ccsolver(ccmodel,Co_f,Cv_f,fock_f,nb,no,nv, &
      & mylsitem,ccPrintLevel,ccenergy,VOVO,longrange_singles,local,JOB, &
      & frag,vovo_supplied,p1,p2,p3,p4,m1,m2,m3,m4)

   implicit none

   !> CC model
   integer,intent(in) :: ccmodel
   !> Number of occupied orbitals in full molecule/fragment AOS
   integer, intent(in)                       :: no
   !> Number of virtual orbitals in full molecule/fragment AOS
   integer, intent(in)                       :: nv
   !> Number of basis functions in full molecule/atomic extent
   integer, intent(in)                       :: nb
   !> Fock matrix in AO basis for fragment or full molecule
   real(realk), target, intent(in) :: fock_f(nb,nb)
   !> Occupied MO coefficients for fragment/full molecule
   real(realk), target, intent(in) :: Co_f(nb,no)
   !> Virtual MO coefficients for fragment/full molecule
   real(realk), target, intent(in) :: Cv_f(nb,nv)
   !> LS item information
   type(lsitem), intent(inout)               :: mylsitem
   !> How much to print? ( ccPrintLevel>0 --> print info stuff)
   integer, intent(in)                       :: ccPrintLevel
   !> Coupled cluster energy for fragment/full molecule
   real(realk),intent(inout)                 :: ccenergy
   !> Two electron integrals (a i | b j) stored as (a,i,b,j)
   type(tensor),intent(inout)                :: VOVO
   !> Include long-range singles effects using singles amplitudes
   !> from previous fragment calculations.
   !> IMPORTANT: If this it TRUE, then the singles amplitudes for the fragment
   !> (from previous calculations) must be stored in t1_final at input!
   logical,intent(in)                        :: longrange_singles
   logical,intent(inout)                     :: local
   integer,intent(in)                        :: JOB
   !additional information and parameter set
   type(tensor), intent(inout), optional     :: m1,m2,m3,m4
   !final parameters 
   type(tensor), intent(inout), optional     :: p1,p2,p3,p4
   type(decfrag), intent(inout), optional    :: frag
   logical, intent(in), optional             :: vovo_supplied
   !
   !> Occ-occ block of Fock matrix in MO basis
   type(tensor)        :: oo_fock
   !> Virt-virt block of Fock matrix in MO basis
   type(tensor)        :: vv_fock

   !> Do an MO-based CCSD calculation?
   logical :: mo_ccsd
   !> full set of MO integrals (non-T1-transformed)
   type(tensor) :: pgmo_diag, pgmo_up
   type(MObatchInfo) :: MOinfo
   !
   !work stuff
   real(realk),pointer :: Co_d(:,:),Cv_d(:,:),focc(:),fvirt(:)
   real(realk),pointer :: Co_b1(:,:),Co_b2(:,:),Cv_b1(:,:),Cv_b2(:,:)
   real(realk),pointer :: fock_(:,:),oofock_d(:,:),vvfock_d(:,:)
   real(realk) :: ccenergy_check
   integer, dimension(2) :: occ_dims, virt_dims, ao2_dims, ampl2_dims, ord2
   integer, dimension(4) :: ampl4_dims
   type(tensor)  :: fock,Co,Cv,Uo,Uv
   type(tensor)  :: oofock,vvfock,ovfock,vofock, t1fock
   type(tensor)  :: iajb
   type(tensor), pointer :: t2(:),omega2(:)
   type(tensor), pointer :: t1(:),omega1(:)
   type(tensor) :: omega1_opt, t1_opt, omega1_prec
   type(tensor) :: omega2_opt, t2_opt, omega2_prec, u
   type(tensor) :: xo,yo,xv,yv,h1
   type(tensor) :: Lmo, Co_orig, Cv_orig
   real(realk)             :: test_norm,two_norm_total, one_norm_total,one_norm1, one_norm2, prev_norm
   real(realk), pointer    :: B(:,:),c(:)
   integer                 :: iter,last_iter,i,j,k,l,iter_idx,next,nSS,ib,ie
   logical                 :: crop_ok,break_iterations,saferun, use_pseudo_diag_basis, use_pnos, get_multipliers
   type(ri)                :: l_ao
   type(tensor)            :: oofock_prec, vvfock_prec, vofock_prec, ovfock_prec
   real(realk)             :: tcpu, twall, ttotend_cpu, ttotend_wall, ttotstart_cpu, ttotstart_wall
   real(realk)             :: iter_cpu,iter_wall
   integer                 :: nnodes, t2idx, t1idx
   logical                 :: fragment_job
   type(PNOSpaceInfo), pointer :: pno_cv(:), pno_S(:)
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
   real(realk) :: t2fnorm2, t1fnorm2
   !SOME DUMMIES FOR TESTING
   type(array4) :: iajb_a4
   type(tensor)            :: tmp1,tmp2,tmp3
   character(tensor_MSG_LEN) :: msg
   integer(kind=8)        :: o2v2
   real(realk)            :: mem_o2v2, MemFree
   integer                :: ii, jj, aa, bb, cc, old_iter, nspaces, os, vs, counter
   logical                :: restart, restart_from_converged,collective,use_singles,vovo_avail
   logical                :: trafo_vovo, trafo_m4, trafo_m2,bg_was_init, bg
   logical                :: diag_oo_block, diag_vv_block, prec, prec_not1, prec_in_b
   character(4)           :: atype
   character              :: intspec(5)

   call time_start_phase(PHASE_WORK, twall = ttotstart_wall, tcpu = ttotstart_cpu )

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
   fragment_job     = present(frag)
   
   o2v2             = (i8*nv**2)*no**2
   mem_o2v2         = (8.0E0_realk*o2v2)/(1.024E3_realk**3)

   ! Set integral info
   ! *****************
   !R = Regular Basis set on the 1th center 
   !R = Regular Basis set on the 2th center 
   !R = Regular Basis set on the 3th center 
   !R = Regular Basis set on the 4th center 
   !C = Coulomb operator
   !E = Long-Range Erf operator
   if (mylsitem%setting%scheme%CAM) then
      intspec = ['R','R','R','R','E']
   else
      intspec = ['R','R','R','R','C']
   endif


   call get_currently_available_memory(MemFree)

   !Set defaults
   restart          = .false.
   saferun          = (.not.DECinfo%CCSDnosaferun.or.(DECinfo%only_n_frag_jobs>0))
   prec             = .true.
   prec_not1        = .false.

   if( saferun .and. (2.0E0_realk*mem_o2v2 > 0.5E0_realk*MemFree) )then
      print *,"WARNING(ccsolver_par): detected high memory requirements, I will"
      print *,"therefore not save any amplitudes to file. This requires a makeover"

      saferun = .false.
   endif

   nnodes      = 1

#ifdef VAR_MPI
   nnodes      = infpar%lg_nodtot
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
      if(.not.present(p2))then
         call lsquit('ccsolver: Long range singles corrections requested, &
            & but p2 not present!',DECinfo%output)
      endif
      if(.not. p2%initialized) then
         call lsquit('ccsolver: Long range singles corrections requested, &
            & but t1_final does not contain existing amplitudes!',DECinfo%output)
      end if
   end if

   !Settings for the models
   !MODIFY FOR NEW MODEL
   ModelSpecificSettings: select case(ccmodel)
   case( MODEL_MP2 )

      use_singles = .false.
      atype = 'REAR'

   case( MODEL_RIMP2,MODEL_LSTHCRIMP2 )

      call lsquit("ERROR(ccsolver_par): RI-MP2 is not implemented iteratively",-1)

      use_singles = .false.
      atype = 'REAR'

   case( MODEL_CC2, MODEL_CCSD )

      if(.not.present(p2))then
         call lsquit("ERROR(ccsolver_par): singles parameters need to be present in p2",-1)
      endif
      use_singles = .true.
      atype = 'LDAR'

   case(MODEL_RPA)

      use_singles = .false.
      atype = 'LDAR'

   case(MODEL_SOSEX)

      use_singles = .false.
      atype = 'LDAR'

   case default

      call lsquit("ERROR(ccsolver_par): requested model not yet implemented",-1)

   end select ModelSpecificSettings

   call ccdriver_set_tensor_segments_and_alloc_workspace(MyLsitem,nb,no,nv,os,vs,local,saferun,use_singles,JOB,ccmodel,bg_was_init)
   bg = mem_is_background_buf_init()

   ! dimension vectors
   occ_dims   = [nb,no]
   virt_dims  = [nb,nv]
   ao2_dims   = [nb,nb]
   ampl4_dims = [nv,nv,no,no]
   ampl2_dims = [nv,no]

   ! initialize T1 matrices and fock transformed matrices for CC pp,pq,qp,qq
   if(CCmodel /= MODEL_MP2 .and. ccmodel /= MODEL_RPA&
     &.and. ccmodel /= MODEL_SOSEX) then
      call tensor_minit(xo, occ_dims, 2, local=local, atype='LDAR', bg=bg )
      call tensor_minit(yo, occ_dims, 2, local=local, atype='LDAR', bg=bg )
      call tensor_minit(xv, virt_dims,2, local=local, atype='LDAR', bg=bg )
      call tensor_minit(yv, virt_dims,2, local=local, atype='LDAR', bg=bg )
   end if

   ! allocate things
   if(use_singles) then
      call mem_alloc( t1,     DECinfo%ccMaxIter )
      call mem_alloc( omega1, DECinfo%ccMaxIter )

      if( bg )then
         do i=1,DECinfo%ccMaxDIIS
            call tensor_minit(t1(i), ampl2_dims, 2, local=local, atype='REPD', bg=bg )
         enddo
      else
         call tensor_minit(t1(1), ampl2_dims, 2, local=local, atype='REPD' )
      endif
   end if

   call mem_alloc( t2,     DECinfo%ccMaxDIIS )
   call mem_alloc( omega2, DECinfo%ccMaxDIIS )

   if( bg )then
      do i=1,DECinfo%ccMaxDIIS
         call tensor_minit(t2(i), ampl4_dims, 4, local=local, atype='TDAR', tdims=[vs,vs,os,os], bg=bg )
      enddo
   else
      call tensor_minit(t2(1), ampl4_dims, 4, local=local, atype='TDAR', tdims=[vs,vs,os,os] )
   endif

   ! Fock matices and density corresponding to input MOs
   call tensor_minit( fock, ao2_dims, 2, local=local, atype=atype, bg=bg )
   call tensor_convert( fock_f,  fock )

   call tensor_minit( t1fock, ao2_dims, 2, local=local, atype=atype , bg=bg )

   !set job related
   select case(JOB)
   case(SOLVE_AMPLITUDES)

      use_pnos        = .false.
      get_multipliers = .false.
      prec_in_b       = .true.
      trafo_m4        = .false.
      trafo_m2        = .false.

      Co_b1 => Co_f
      Co_b2 => Co_f

      Cv_b1 => Cv_f
      Cv_b2 => Cv_f

      fock_ => fock_f

   case(SOLVE_AMPLITUDES_PNO)

      use_pnos        = .true.
      get_multipliers = .false.
      prec_in_b       = .true.

      if(.not.present(m4)) then
         call lsquit('ccsolver: When using PNOs make sure MP2 amplitudes are &
            & in m4',DECinfo%output)
      end if

      trafo_m4 = .true.
      trafo_m2 = .false.

      Co_b1 => Co_f
      Co_b2 => Co_f

      Cv_b1 => Cv_f
      Cv_b2 => Cv_f

      fock_ => fock_f

   case(SOLVE_MULTIPLIERS)

      use_pnos        = .false.
      get_multipliers = .true.
      prec_in_b       = .true.

      if(.not.present(m4)) then
         call lsquit('ccsolver: When solving for the multipliers, make sure the&
         & doubles amplitudes are given in m4',DECinfo%output)
      end if
      trafo_m4 = .true.

      if(use_singles)then
         if(.not.present(m2)) then
            call lsquit('ccsolver: When solving for the multipliers, make sure the&
               & singles amplitudes are given in m2',DECinfo%output)
         end if
         trafo_m2 = .true.
      else
         trafo_m2 = .false.
      endif

      call tensor_minit(Co_orig, occ_dims, 2, local=local, atype='LDAR',bg=bg )
      call tensor_minit(Cv_orig, virt_dims,2, local=local, atype='LDAR',bg=bg )

      Co_b1 => Co_f
      Cv_b1 => Cv_f
      call tensor_convert( Co_b1, Co_orig )
      call tensor_convert( Cv_b1, Cv_orig )

      call get_t1_matrices(MyLsitem,m2,Co_orig,Cv_orig,xo,yo,xv,yv,fock,t1fock,.false.)

      call tensor_free(Cv_orig)
      call tensor_free(Co_orig)

      Co_b1 => xo%elm2
      Co_b2 => yo%elm2

      Cv_b1 => xv%elm2
      Cv_b2 => yv%elm2

      fock_ => t1fock%elm2


      !SOLVE IN FOCK DIAG BASIS
      Co_b1 => Co_f
      Co_b2 => Co_f

      Cv_b1 => Cv_f
      Cv_b2 => Cv_f

      fock_ => fock_f
   case default
      call lsquit("ERROR(ccsolver_par):unknown jobid",-1)
   end select
 
   if( DECinfo%ccsolver_overwrite_prec )then
      prec      = DECinfo%use_preconditioner
      prec_not1 = DECinfo%precondition_with_full
      prec_in_b = DECinfo%use_preconditioner_in_b
   endif

   vovo_avail = .false.
   trafo_vovo = .false.
   if(present(vovo_supplied)) vovo_avail = vovo_supplied

   if(vovo_avail)then

      if(.not.VOVO%initialized)then
         call lsquit("ERROR(ccsolver_par): input specified that vovo was&
         & supplied, but vovo is not initialized",-1)
      endif

      if(VOVO%itype == TT_TILED_DIST .and.(VOVO%tdim(1)/=vs .or.  VOVO%tdim(2)/=os))then
         call lsquit("ERROR(ccsolver_par): distribution of VOVO not fit for solver",-1)
      endif

      trafo_vovo = .true.

   endif

   ! go to a (pseudo) canonical basis
   ! --------------------------------

   call mem_alloc( focc,     no     )
   call mem_alloc( fvirt,    nv     )
   call mem_alloc( Co_d,     nb, no )
   call mem_alloc( Cv_d,     nb, nv )
   call mem_alloc( oofock_d, no, no )
   call mem_alloc( vvfock_d, nv, nv )
   call tensor_minit(Uo,   [no,no],  2, local=local, atype="TDPD", tdims = [os,os] )
   call tensor_minit(Uv,   [nv,nv],  2, local=local, atype="TDPD", tdims = [vs,vs] )
   call tensor_minit( oo_fock, [no,no], 2, bg=bg )
   call tensor_minit( vv_fock, [nv,nv], 2, bg=bg )

   !misuse Co_d and Cv_d as workspaces
   call dgemm('n','n',nb,no,nb,1.0E0_realk,fock_,nb,Co_b2,nb,0.0E0_realk,Co_d,nb)
   call dgemm('t','n',no,no,nb,1.0E0_realk,Co_b1,nb,Co_d,nb,0.0E0_realk,oo_fock%elm1,no)

   call dgemm('n','n',nb,nv,nb,1.0E0_realk,fock_,nb,Cv_b2,nb,0.0E0_realk,Cv_d,nb)
   call dgemm('t','n',nv,nv,nb,1.0E0_realk,Cv_b1,nb,Cv_d,nb,0.0E0_realk,vv_fock%elm1,nv)

   Co_b1 => null()
   Co_b2 => null()

   Cv_b1 => null()
   Cv_b2 => null()

   use_pseudo_diag_basis = .not.(DECinfo%CCSDpreventcanonical)

   diag_oo_block = use_pseudo_diag_basis.and..not.use_pnos
   diag_vv_block = use_pseudo_diag_basis

   if(use_pseudo_diag_basis)then
      call get_canonical_integral_transformation_matrices(no,nv,nb,oo_fock%elm1,vv_fock%elm1,Co_f,Cv_f,&
         & Co_d,Cv_d,Uo%elm2,Uv%elm2,focc,fvirt)
   endif

   !OO BLOCK
   if( diag_oo_block )then
      oofock_d = 0.0E0_realk
      do ii=1,no
         oofock_d(ii,ii) = focc(ii)
      enddo
   else
      Co_d     = Co_f
      oofock_d = oo_fock%elm2
      Uo%elm2  = 0.0E0_realk
      do ii=1,no
         Uo%elm2(ii,ii) = 1.0E0_realk
      enddo
   endif

   !VV BLOCK
   if( diag_vv_block )then
      vvfock_d = 0.0E0_realk
      do aa=1,nv
         vvfock_d(aa,aa) = fvirt(aa)
      enddo
   else
      Cv_d     = Cv_f
      vvfock_d = vv_fock%elm2
      Uv%elm2  = 0.0E0_realk
      do aa=1,nv
         Uv%elm2(aa,aa) = 1.0E0_realk
      enddo
   endif

   call tensor_free(vv_fock)
   call tensor_free(oo_fock)

   call mem_dealloc( focc  )
   call mem_dealloc( fvirt )
   call tensor_mv_dense2tiled(Uo,.not.local,dealloc_local = .false.)
   call tensor_mv_dense2tiled(Uv,.not.local,dealloc_local = .false.)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !transform input to the diagonal basis!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if(trafo_m2)   call local_can_trans(no,nv,nb,Uo%elm2,Uv%elm2,vo=m2%elm1)
   if(trafo_m4)   call tensor_transform_basis([Uo,Uv],2,[ m4   ], [[2,1,2,1]], [[1,1,1,1]], 4, 1)
   if(trafo_vovo) call tensor_transform_basis([Uo,Uv],2,[ VOVO ], [[2,1,2,1]], [[1,1,1,1]], 4, 1)

   if(fragment_job.and.use_pnos)then
      call local_can_trans(no,nv,nb,Uo%elm2,Uv%elm2,vv=frag%VirtMat,oo=frag%OccMat)
   endif
   if(JOB==SOLVE_MULTIPLIERS)then
      call local_can_trans(no,nv,nb,Uo%elm2,Uv%elm2,bo=xo%elm1,bv=xv%elm1)
      call local_can_trans(no,nv,nb,Uo%elm2,Uv%elm2,bo=yo%elm1,bv=yv%elm1)
   endif

   ! create transformation matrices in array form
   call tensor_minit( Co , occ_dims, 2, local=local, atype="REAR", bg=bg )
   call tensor_minit( Cv , virt_dims,2, local=local, atype="REAR", bg=bg )

   call tensor_convert( Co_d, Co )
   call tensor_convert( Cv_d, Cv )

   call mem_dealloc( Co_d )
   call mem_dealloc( Cv_d )

   if(DECinfo%PL>1)call time_start_phase(PHASE_work, at = time_work, ttot = time_fock_mat, &
      &twall = time_prec1 , labelttot = 'CCSOL: AO FOCK MATRIX :', output = DECinfo%output)


   ! get fock matrices, used in Preconditioning and MP2
   ! --------------------------------------------------

   call tensor_minit(oofock_prec, [no,no], 2, local=local, atype='REPD', bg=bg )
   call tensor_minit(vvfock_prec, [nv,nv], 2, local=local, atype='REPD', bg=bg )
   call tensor_minit(vofock_prec, [nv,no], 2, local=local, atype='REPD' )
   if( JOB == SOLVE_MULTIPLIERS )then
      call tensor_minit(ovfock_prec, [no,nv], 2, local=local, atype='REPD' )
   endif

   call tensor_change_atype_to_rep( oofock_prec, local )
   call tensor_change_atype_to_rep( vvfock_prec, local )

   if( prec_not1 ) then
      call tensor_convert( oofock_d, oofock_prec )
      call tensor_convert( vvfock_d, vvfock_prec )
      call tensor_zero(vofock_prec)
   else
      call tensor_minit(tmp1, [nb,no], 2, local=local, atype='LDAR', bg=bg )
      ord2 = [1,2]
      call tensor_contract(1.0E0_realk,fock,Co,[2],[1],1,0.0E0_realk,tmp1,ord2)
      call tensor_contract(1.0E0_realk,Co,tmp1,[1],[1],1,0.0E0_realk,oofock_prec,ord2)
      call tensor_free(tmp1)

      call tensor_minit(tmp1, [nb,nv], 2, local=local, atype='LDAR', bg=bg )
      call tensor_contract(1.0E0_realk,fock,Cv,[2],[1],1,0.0E0_realk,tmp1,ord2)
      call tensor_contract(1.0E0_realk,Cv,tmp1,[1],[1],1,0.0E0_realk,vvfock_prec,ord2)
      call tensor_free(tmp1)

      call tensor_minit(tmp1, [nb,no], 2, local=local, atype='LDAR', bg=bg )
      call tensor_contract(1.0E0_realk,fock,Co,[2],[1],1,0.0E0_realk,tmp1,ord2)
      call tensor_contract(1.0E0_realk,Cv,tmp1,[1],[1],1,0.0E0_realk,vofock_prec,ord2)
      call tensor_free(tmp1)

      if(JOB == SOLVE_MULTIPLIERS )then
         !For the starting guess
         call tensor_minit(tmp1, [nb,nv], 2, local=local, atype='LDAR', bg=bg )
         call tensor_contract(1.0E0_realk,t1fock,Cv,[2],[1],1,0.0E0_realk,tmp1,ord2)
         call tensor_contract(-2.0E0_realk,Co,tmp1,[1],[1],1,0.0E0_realk,ovfock_prec,ord2)
         call tensor_free(tmp1)
      endif

   end if

   if( ccmodel /= MODEL_MP2 .and. ccmodel /= MODEL_RPA &
    &.and. ccmodel /= MODEL_SOSEX )then
      call tensor_change_atype_to_d( oofock_prec )
      call tensor_change_atype_to_d( vvfock_prec )
   endif

   if(DECinfo%PL>1)call time_start_phase(PHASE_work, at = time_work, ttot = time_prec1,&
      &labelttot = 'CCSOL: PRECOND. INIT. :', output = DECinfo%output)

   call mem_dealloc( oofock_d )
   call mem_dealloc( vvfock_d )


   !Get OVOV integrals from non-T1 transformed VOVO integrals or generate
   !----------------------------------------------------------------------
   call tensor_minit(iajb, [no,nv,no,nv], 4, local=local, atype='TDAR', tdims=[os,vs,os,vs], bg=bg )

   if(.not.vovo_avail)then
      call tensor_zero(iajb)
      call get_mo_integral_par( iajb, Co, Cv, Co, Cv, mylsitem, intspec, local, collective )
   else
      call tensor_cp_data(VOVO, iajb, order = [2,1,4,3])
      !call tensor_free(VOVO)
   endif

   call mem_alloc( B, DECinfo%ccMaxDIIS, DECinfo%ccMaxDIIS )
   call mem_alloc( c, DECinfo%ccMaxDIIS                    )

   call time_start_phase( PHASE_work, at = time_work, twall = time_start_guess )

   ! get guess amplitude vectors in the first iteration --> generate or read t*.restart
   ! ----------------------------------------------------------------------------------
   two_norm_total = DECinfo%ccConvergenceThreshold + 1.0E0_realk
   if(use_singles)then

      call get_guess_vectors(ccmodel,JOB,prec,restart,old_iter,nb,two_norm_total,ccenergy,t2(1),iajb,Co,Cv,Uo%elm2,Uv%elm2,&
         & oofock_prec,vvfock_prec,vofock_prec, ovfock_prec, mylsitem, local, safefilet21,safefilet22, safefilet2f, &
         & t1(1),safefilet11,safefilet12, safefilet1f  )

   else

      call get_guess_vectors(ccmodel,JOB,prec,restart,old_iter,nb,two_norm_total,ccenergy,t2(1),iajb,Co,Cv,Uo%elm2,Uv%elm2,&
         & oofock_prec,vvfock_prec,vofock_prec,ovfock_prec,mylsitem,local,safefilet21,safefilet22, safefilet2f )

   endif
   restart_from_converged = (two_norm_total < DECinfo%ccConvergenceThreshold)

   call tensor_free(vofock_prec)
   if( JOB == SOLVE_MULTIPLIERS )then
      call tensor_free(ovfock_prec)
   endif

   if(DECinfo%PL>1)call time_start_phase( PHASE_WORK, at = time_work, ttot = time_start_guess,&
      &labelttot = 'CCSOL: STARTING GUESS :', output = DECinfo%output, twall = time_main  )

   ! Print Job Header
   ! ----------------
   Call print_ccjob_header(ccmodel,ccPrintLevel,fragment_job,&
      &(JOB==SOLVE_MULTIPLIERS),nb,no,nv,DECinfo%ccMaxDIIS,restart,restart_from_converged,old_iter)


   If_not_converged: if(.not.restart_from_converged)then


      !set special quantities, if the mos are precalculated set mos, if a PNO
      !calculation is requested, set the pno space information
      call ccsolver_get_special(ccmodel,mylsitem,no,nv,nb,use_pnos,mo_ccsd,Co,Cv,pgmo_diag,pgmo_up,&
         &MOinfo,nspaces,pno_cv,pno_S,m4=m4,frag=frag)

      break_iterations = .false.
      crop_ok          = .false.
      prev_norm        = huge(prev_norm)

      if(bg)then
         if(use_singles)then
            do i=1,DECinfo%ccMaxDIIS
               call tensor_minit(omega1(i), ampl2_dims, 2 , local=local, atype='LDAR',bg=bg )
            enddo
         endif
         do i=1,DECinfo%ccMaxDIIS
            call tensor_minit(omega2(i), ampl4_dims, 4, local=local,&
               & atype='TDAR', tdims=[vs,vs,os,os], bg=bg )
         enddo
      endif

      ! Iterate the equations
      ! ---------------------
      CCIteration : do iter=1,DECinfo%ccMaxIter

         if(DECinfo%PL>1)call time_start_phase( PHASE_work, at = time_work, twall = time_iter ) 

         iter_idx = get_iter_idx(iter    )
         next     = get_iter_idx(iter + 1)
         nSS      = min(iter,DECinfo%ccMaxDIIS)

         ! Initialize residual vectors
         ! ---------------------------
         if(use_singles)then
            if(iter <= DECinfo%ccMaxDIIS .and. .not. bg)then
               call tensor_minit(omega1(iter_idx), ampl2_dims, 2 , local=local, atype='LDAR' )
            endif
            call tensor_zero(omega1(iter_idx))
         endif
         if(iter <= DECinfo%ccMaxDIIS .and. .not. bg)then
            call tensor_minit(omega2(iter_idx), ampl4_dims, 4, local=local,&
               & atype='TDAR', tdims=[vs,vs,os,os] )
         endif
         call tensor_zero(omega2(iter_idx))


         if(DECinfo%PL>1)call time_start_phase( PHASE_work, at = time_work, twall = time_t1_trafo ) 

         ! Get special t1 matrices, else just the normal fock and trafo matrices
         ! ---------------------------------------------------------------------
         T1Related : if(use_singles) then
            select case( JOB )
            case( SOLVE_AMPLITUDES, SOLVE_AMPLITUDES_PNO)
               call get_t1_matrices(MyLsitem,t1(iter_idx),Co,Cv,xo,yo,xv,yv,fock,t1fock,.true.)
            end select
         else
            t1fock%elm2=fock%elm2
         end if T1Related

         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_t1_trafo, &
            &labelttot= 'CCIT: T1 TRAFO        :', output = DECinfo%output, twall = time_residual ) 

         !Get the residual r = Ax - b for any of the implemented models
         !-------------------------------------------------------------
         call ccsolver_get_residual(ccmodel,JOB,omega2,t2,fock,t1fock,iajb,no,nv,oofock_prec,vvfock_prec,xo,xv,yo,yv,nb,MyLsItem,&
            &omega1,t1,pgmo_diag,pgmo_up,MOinfo,mo_ccsd,pno_cv,pno_s,nspaces,iter,local,use_pnos,restart,frag=frag,m2=m2,m4=m4)


         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_residual, &
            &labelttot= 'CCIT: RESIDUAL        :', output = DECinfo%output, twall = time_crop_mat ) 

         !calculate crop matrix
         !---------------------
         call ccsolver_calculate_crop_matrix(B,nSS,omega2,omega1,oofock_prec,vvfock_prec,ccmodel,local,prec_in_b)

         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_crop_mat, &
            &labelttot= 'CCIT: CROP MATRIX     :', output = DECinfo%output, twall = time_solve_crop ) 

         ! solve crop/diis equation
         !-------------------------
         call CalculateDIIScoefficientsII(nSS,B,c,DECinfo%PL>3)

         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_solve_crop, &
            &labelttot= 'CCIT: SOLVE CROP      :', output = DECinfo%output, twall = time_mixing ) 


         ! Calculate linear combinations of the residual and parameter vectors
         ! -------------------------------------------------------------------
         if(use_singles) then
            call tensor_minit(omega1_opt, ampl2_dims, 2 , local=local, atype='LDAR',bg=bg )
            call tensor_minit(t1_opt    , ampl2_dims, 2 , local=local, atype='LDAR',bg=bg )
            call tensor_zero(omega1_opt)
            call tensor_zero(t1_opt    )
         end if

         call tensor_minit(omega2_opt, ampl4_dims, 4, local=local, atype='TDAR', tdims=[vs,vs,os,os], bg=bg )
         call tensor_minit(t2_opt    , ampl4_dims, 4, local=local, atype='TDAR', tdims=[vs,vs,os,os], bg=bg )

         call tensor_zero( omega2_opt )
         call tensor_zero( t2_opt     )

         do i=1,nSS

            ! mix singles
            if(use_singles) then
               call tensor_add( omega1_opt, c(i), omega1(i) )
               call tensor_add( t1_opt,     c(i), t1(i)     )
            end if

            ! mix doubles
            call tensor_add( omega2_opt, c(i), omega2(i) )
            call tensor_add( t2_opt,     c(i), t2(i)     )

         end do

         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_mixing, &
            &labelttot= 'CCIT: MIXING          :', output = DECinfo%output, twall = time_norm ) 

         ! Check for the convergence
         ! -------------------------
         one_norm1 = 0.0E0_realk
         one_norm2 = 0.0E0_realk
         if(use_singles)call print_norm(omega1(iter_idx),one_norm1,.true.)
         call print_norm(omega2(iter_idx),one_norm2,.true.)
         one_norm_total = one_norm1 + one_norm2
         two_norm_total = sqrt(one_norm_total)
         test_norm      = two_norm_total

         ! intentionally quit prematurely
         ! ------------------------------
         if(iter==5.and.((DECinfo%CRASHCALC.and.DECinfo%full_molecular_cc).or.DECinfo%ccsolverskip))then
            if(.not.DECinfo%ccsolverskip)then
               print*,'Calculation was intentionally crashed due to keyword .CRASHCALC'
               print*,'This keyword is only used for debug and testing purposes'
               print*,'We want to be able to test the .RESTART keyword'
               print*,'In the CC case only quit prematurely, then this keyword is even more handy'
               WRITE(DECinfo%output,*)'Calculation was intentionally crashed due to keyword .CRASHCALC'
               WRITE(DECinfo%output,*)'This keyword is only used for debug and testing purposes'
               WRITE(DECinfo%output,*)'We want to be able to test the .RESTART keyword'
               print*,"SETTING TEST_NORM TO QUIT"
            endif
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


         if( JOB == SOLVE_AMPLITUDES.or.JOB == SOLVE_AMPLITUDES_PNO )then

            ! calculate the correlation energy and fragment energy
            ! ----------------------------------------------------


            ! MODIFY FOR NEW MODEL
            ! If you implement a new model, please insert call to energy routine here,
            ! or insert a call to get_cc_energy if your model uses the standard CC energy expression.

            EnergyForCCmodel: select case(CCmodel)
            case( MODEL_MP2 )

               ccenergy = get_mp2_energy(t2(iter_idx),iajb,no,nv)

            case( MODEL_RIMP2,MODEL_LSTHCRIMP2)

               ccenergy = get_mp2_energy(t2(iter_idx),iajb,no,nv)

            case( MODEL_CC2, MODEL_CCSD )

               ! CC2, CCSD, or CCSD(T) (for (T) calculate CCSD contribution here)
               ccenergy = get_cc_energy(t1(iter_idx),t2(iter_idx),iajb,no,nv)

            case(MODEL_RPA)

               ccenergy = get_RPA_energy_arrnew(t2(iter_idx),iajb,no,nv)

               if(DECinfo%SOS) then
                  ccenergy =ccenergy+get_SOSEX_cont_arrnew(t2(iter_idx),iajb,no,nv)
               endif

            case(MODEL_SOSEX)
               ccenergy = get_RPA_energy_arrnew(t2(iter_idx),iajb,no,nv)
               ccenergy =ccenergy+get_SOSEX_cont_arrnew(t2(iter_idx),iajb,no,nv)

            case default
               ! MODEL RIMP2 defaults here since it shoule not use this solver

               call lsquit("ERROR(ccsolver_par):energy expression for your model&
                  & not yet implemented",-1)

            end select EnergyForCCmodel

         endif


         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_energy, &
            &labelttot= 'CCIT: GET ENERGY      :', output = DECinfo%output, twall = time_copy_opt) 


         ! if crop, put the optimal in place of trial (not for diis)
         ! ---------------------------------------------------------

         if(DECinfo%use_crop) then
            if(use_singles) then
               call tensor_cp_data( omega1_opt, omega1(iter_idx) )
               call tensor_cp_data( t1_opt,     t1(iter_idx)     )
            end if
            call tensor_cp_data( omega2_opt, omega2(iter_idx) )
            call tensor_cp_data( t2_opt,     t2(iter_idx)     )
         end if

         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_copy_opt, &
            &labelttot= 'CCIT: COPY OPTIMALS   :', output = DECinfo%output, twall = time_new_guess ) 


         ! generate next trial vector if this is not the last iteration
         ! ------------------------------------------------------------
         if(.not.break_iterations) then
            if( prec ) then
               if(use_singles) then
                  omega1_prec = precondition_singles(omega1_opt,oofock_prec,vvfock_prec)
                  if(iter+1<=DECinfo%ccMaxDIIS.and. .not. bg)then
                     call tensor_minit(t1(next),  ampl2_dims, 2, local=local, atype='REPD' )
                  endif
                  call tensor_cp_data(t1_opt,t1(next))
                  call tensor_add(t1(next),1.0E0_realk,omega1_prec)
                  call tensor_free(omega1_prec)
               end if
               omega2_prec = precondition_doubles(omega2_opt,oofock_prec,vvfock_prec,local)

               if(iter+1<=DECinfo%ccMaxDIIS.and..not.bg)then
                  call tensor_minit(t2(next), ampl4_dims, 4, local=local, atype='TDAR', tdims=[vs,vs,os,os] )
               endif
               call tensor_cp_data(t2_opt,t2(next))
               call tensor_add(t2(next),1.0E0_realk,omega2_prec)
               call tensor_free(omega2_prec)

            else
               if(use_singles)then
                  if(iter+1<=DECinfo%ccMaxDIIS.and..not.bg)then
                     call tensor_minit(t1(next), ampl2_dims, 2, local=local, atype='REPD' )
                  endif
                  call tensor_cp_data(t1_opt,t1(next))
                  call tensor_add(t1(next),1.0E0_realk,omega1_opt)
               endif
               if(iter+1<=DECinfo%ccMaxDIIS.and..not.bg)then
                  call tensor_minit(t2(next), ampl4_dims, 4, local=local, atype='TDAR', tdims=[vs,vs,os,os] )
               endif
               call tensor_cp_data(t2_opt,t2(next))
               call tensor_add(t2(next),1.0E0_realk,omega2_opt)
            end if
         end if


         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_new_guess, &
            &labelttot= 'CCIT: NEW GUESS VEC.  :', output = DECinfo%output ) 


         ! delete optimals
         ! ---------------
         call tensor_free(t2_opt)
         call tensor_free(omega2_opt)
         if(use_singles) then
            call tensor_free(t1_opt)
            call tensor_free(omega1_opt)
         end if


         ! Checkpoint the calculation
         ! --------------------------
         if(saferun.and..not.break_iterations)then

            if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, twall = time_write ) 

            if(use_singles)then
               call save_current_guess(local,iter+old_iter,nb,two_norm_total,ccenergy,Uo%elm2,Uv%elm2,t2(next),safefilet21,&
                  &safefilet22,t1(next),safefilet11,safefilet12)                  
            else                                                                     
               call save_current_guess(local,iter+old_iter,nb,two_norm_total,ccenergy,Uo%elm2,Uv%elm2,t2(next),safefilet21,&
                  &safefilet22)
            endif

            if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_write, &
               &labelttot= ' CCIT: NEW GUESS WRITE :', output = DECinfo%output) 

         endif


         if(DECinfo%PL>1) call time_start_phase( PHASE_work, at = time_work, ttot = time_iter, &
            &labelttot= 'CCIT: ITERATION       :', output = DECinfo%output) 

         call print_ccjob_iterinfo(iter+old_iter,two_norm_total,ccenergy,(JOB==SOLVE_MULTIPLIERS),fragment_job)

         last_iter = iter
         if(break_iterations) exit

      end do CCIteration

      !Free residuals
      if(bg)then
         ib = DECinfo%ccMaxDIIS
         ie = 1
      else
         ib = last_iter
         ie = max(last_iter-DECinfo%ccMaxDIIS+1,1)
      endif
      do i=ib,ie,-1
         iter_idx = get_iter_idx(i)
         call tensor_free(omega2(iter_idx))
      end do
      if(use_singles) then
         do i=ib,ie,-1
            iter_idx = get_iter_idx(i)
            call tensor_free( omega1(iter_idx) )
         end do
      endif

   else

      call print_ccjob_iterinfo(old_iter,two_norm_total,ccenergy,(JOB==SOLVE_MULTIPLIERS),fragment_job)

      break_iterations = .true.
      last_iter        = 1

   endif If_not_converged


   call time_start_phase( PHASE_work, at = time_work, ttot = time_main, labelttot = 'CCSOL: MAIN LOOP      :', &
      & output = DECinfo%output, twall = time_finalize )


   ! Free memory and save final amplitudes
   ! *************************************
   ! Save two-electron integrals in the order (virt,occ,virt,occ), save the used RHS or restore the old rhs
   if(.not.vovo_avail)then
      call tensor_minit( tmp2, [nv,no,nv,no], 4, local=local, tdims=[vs,os,vs,os],atype = "TDAR" )
      call tensor_cp_data(iajb, tmp2, order = [2,1,4,3] )
   endif
   call tensor_free(iajb)

   ! remove preconditioning matrices
   call tensor_free(oofock_prec)
   call tensor_free(vvfock_prec)

   call tensor_free(Cv)
   call tensor_free(Co)

   call tensor_free(t1fock)
   call tensor_free(fock)

   ! remove rest of the singles amplitudes and residuals

   iter_idx = get_iter_idx(last_iter)

   ! Save final double amplitudes (to file if saferun)
   call tensor_minit( tmp1, [nv,no,nv,no], 4 , local=local, tdims = [vs,os,vs,os], atype = "TDAR" )
   call tensor_cp_data(t2(iter_idx), tmp1, order = [1,3,2,4] )

   if(use_singles) then
      call tensor_minit(tmp3,[nv,no],2)
      call tensor_cp_data(t1(iter_idx), tmp3 )
   endif

   !SAFE THE FINAL AMPLITUDES, NOT YET REORDERED
   if(saferun.and..not.restart_from_converged.and..not.DECinfo%CRASHCALC)then
      if(use_singles)then
         call save_current_guess(local,last_iter+old_iter,nb,two_norm_total,ccenergy,Uo%elm2,Uv%elm2,&
            &t2(iter_idx),safefilet2f,safefilet2f,t1(iter_idx),safefilet1f,safefilet1f)
      else
         call save_current_guess(local,last_iter+old_iter,nb,two_norm_total,ccenergy,Uo%elm2,Uv%elm2,&
            &t2(iter_idx),safefilet2f,safefilet2f)
      endif
   endif

   !Free workspace
   if(bg)then
      ib = DECinfo%ccMaxDIIS
      ie = 1
   else
      ib = last_iter
      ie = max(last_iter-DECinfo%ccMaxDIIS+1,1)
   endif
   ! Free doubles amplitudes
   do i=ib,ie,-1
      iter_idx = get_iter_idx(i)
      call tensor_free(t2(iter_idx))
   end do
   ! Free singles amplitudes and residuals
   if(use_singles) then
      do i=ib,ie,-1
         iter_idx = get_iter_idx(i)
         call tensor_free( t1(iter_idx) )
      end do
   end if

   if(use_singles) then
      call mem_dealloc(t1)
      call mem_dealloc(omega1)
   end if
   call mem_dealloc(t2)
   call mem_dealloc(omega2)

   call mem_dealloc(B)
   call mem_dealloc(c)


   if(use_singles) then
      !call array2_free(h1)
      call tensor_free(yv)
      call tensor_free(xv)
      call tensor_free(yo)
      call tensor_free(xo)
   end if

   call time_start_phase(PHASE_WORK,at = time_work, twall = ttotend_wall, tcpu = ttotend_cpu )

   call print_norm(tmp1,t2fnorm2,.true.)
   if(use_singles)then
      call print_norm(tmp3,t1fnorm2,.true.)
   endif

   call print_ccjob_summary(break_iterations,(JOB==SOLVE_MULTIPLIERS),fragment_job,&
      &last_iter+old_iter,use_singles,ccenergy,ttotend_wall,&
      &ttotstart_wall,ttotend_cpu,ttotstart_cpu,t1fnorm2,t2fnorm2,nm1=t1fnorm2,nm2=t2fnorm2)
   
   call ccdriver_dealloc_workspace(saferun,local,bg_was_init)

   call ccsolver_free_special(pgmo_up,pgmo_diag,MOinfo,restart_from_converged,&
      &mo_ccsd,pno_cv,pno_S,nspaces,use_pnos,frag)


   !initialize and copy data to output arrays
   if(.not.vovo_avail)then
      call tensor_minit( VOVO, [nv,no,nv,no], 4, local=local, tdims=[vs,os,vs,os],&
         &atype = "TDAR", fo=tmp2%offset, bg=bg_was_init )
      call tensor_cp_data(tmp2, VOVO )
      call tensor_free(tmp2)
   endif

   call tensor_minit( p4, [nv,no,nv,no], 4 , local=local, tdims = [vs,os,vs,os], atype = "TDAR", fo=tmp1%offset, bg=bg_was_init)
   call tensor_cp_data( tmp1, p4 )
   call tensor_free(tmp1)

   if(use_singles) then
      if(.not.longrange_singles) then 
         call tensor_minit(p2,[nv,no],2, bg=bg_was_init)
      end if
      call tensor_cp_data(tmp3, p2 )
      call tensor_free(tmp3)
   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !transform back to original basis!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call tensor_transform_basis( [Uo,Uv], 2, [p4,VOVO], [[2,1,2,1],[2,1,2,1]], [[2,2,2,2],[2,2,2,2]], 4, 2,bg=bg_was_init)

   if(use_singles)then
      !FIXME: all can local trans and local can trans should be replaced by
      !tensor_transform_basis
      call can_local_trans(no,nv,nb,Uo%elm2,Uv%elm2,vo=p2%elm1)
   endif

   if(trafo_m2) call can_local_trans(no,nv,nb,Uo%elm2,Uv%elm2,vo=m2%elm1 )
   if(trafo_m4) call tensor_transform_basis([Uo,Uv],2, [m4], [[2,1,2,1]], [[2,2,2,2]], 4, 1, bg=bg_was_init)

   if(fragment_job.and.use_pnos)then
      call can_local_trans(no,nv,nb,Uo%elm2,Uv%elm2,oo=frag%OccMat,vv=frag%VirtMat)
   endif

   call tensor_free(Uo)
   call tensor_free(Uv)


   if(DECinfo%PL>1)call time_start_phase( PHASE_work, at = time_work, ttot = time_finalize, &
      &labelttot = 'CCSOL: FINALIZATION   :', output = DECinfo%output )

end subroutine ccsolver

subroutine ccdriver_set_tensor_segments_and_alloc_workspace(MyLsitem,nb,no,nv,os,vs,local,saferun,&
     & use_singles,JOB,ccmodel,bg_was_init)
   implicit none
   type(lsitem), intent(inout) :: MyLsitem
   integer, intent(in)  :: nb, no, nv, JOB, ccmodel
   integer, intent(out) :: os,vs
   logical, intent(in)  :: local,saferun,use_singles
   logical, intent(inout) :: bg_was_init
   integer :: ntpm(4),nt, sch, buf
   integer(kind=ls_mpik) :: nnod
   real(realk) :: MemFree
   integer(kind=8) :: mem41,mem31,mem21
   integer(kind=8) :: mem42,mem32,mem22
   integer(kind=8) ::mem_out, mem_in, w1,w2, mem_o, mem_i, mem_in_basic,mem_out_basic
   integer(kind=8) :: bytes, Freebytes, bytes_to_alloc
   integer(kind=8) :: nelms_ccdriver, nelms_int0,nelms_int1,nelms_int2,nelms_int3, nelms_res_in, nelms_res_out
   logical :: use_bg
   integer(kind=8) :: elm
   integer :: MinAO,iAO,s,nbuffs, bs

   nnod = 1
#ifdef VAR_MPI
   nnod = infpar%lg_nodtot
#endif

   !get tensor segments
   call get_symm_tensor_segmenting_simple(no,nv,os,vs)

   ! allocate the buffer in the background
   use_bg = DECinfo%use_bg_buffer.and..not.saferun
   bg_was_init = mem_is_background_buf_init()
   if(use_bg)then

      if( .not. bg_was_init )then
         !get free memory GB -> bytes conversion, use 80% of the free memory
         call get_currently_available_memory(MemFree)
         Freebytes = int(MemFree * (1024.0E0_realk**3) * 0.8E0_realk,kind=8)

         !CCSOLVER OUTSIDE BG BUFFER
         !**************************
         !estimate memory use in solver while the residual is build, the other
         ! FOck matrix, MO trafo matrices (bv,bo), Lambda matrices, basis trafo matrices, fock-blocks
         nelms_ccdriver = 1_long*nb*no + 1_long*nb*nv + 1_long*no**2 + 1_long*nv**2 + 1_long*no*nv
         !singles vectors and residuals, and t1 fock matrix
         if( use_singles )then
            elm = i8*nb**2
            nelms_ccdriver = nelms_ccdriver + elm
         endif
         mem_out_basic = nelms_ccdriver * 8
         mem_out = nelms_ccdriver * 8

         !CCSOLVER INSIDE BG BUFFER
         !*************************
         ! residual-, trial- and integral vectors, additional two tiles for buffering
         !singles vectors and residuals, and t1 fock matrix
         nelms_ccdriver = i8*nb**2 + 2_long*nb*no + 2_long*nb*nv + 1_long*no**2 + 1_long*nv**2 + 1_long*no*nv
         if( use_singles )then
            elm = (i8*nv)*no * (2*DECinfo%ccMaxDIIS )
            nelms_ccdriver = nelms_ccdriver + elm
         endif
         if( local )then
            elm = (i8*nv*nv)*no*no * (2*DECinfo%ccMaxDIIS + 1)
         else
            !memory requirements are negligible 
            call tensor_get_ntpm([nv,nv,no,no],[vs,vs,os,os],4,ntpm,ntiles=nt)
            nt = ceiling(float(nt)/float(nnod))

            elm = nt * (i8*vs**2) * os**2 * (2*DECinfo%ccMaxDIIS + 3 )  
         endif
         nelms_ccdriver = nelms_ccdriver + elm
         mem_in_basic = nelms_ccdriver * 8


         !estimate memory use in residual without workspaces
         if(DECinfo%useIchor)then
            iAO = 4 
            call determine_MinimumAllowedAObatchSize(MyLsItem%setting,iAO,'R',MinAO)
         else
            call determine_maxBatchOrbitalsize(DECinfo%output,MyLsItem%setting,MinAO,'R')
         endif

         !find buffer size for integral calculation
         nbuffs = ceiling(float(nv)/float(vs))
         bs     = get_split_scheme_0(nb)
         s      = 0
         w1     = get_work_array_size(1,nv,no,nv,no,vs,os,vs,os,nb,bs,MinAO,MinAO,s,nbuffs,MyLsItem%setting)
         w2     = get_work_array_size(2,nv,no,nv,no,vs,os,vs,os,nb,bs,MinAO,MinAO,s,nbuffs,MyLsItem%setting)
         nelms_int0 = w1 + w2 + (i8*no*nv)*no*nv

         s = 1
         nelms_int1 = w1 + w2


         s = 2
         w1 = get_work_array_size(1,nv,no,nv,no,vs,os,vs,os,nb,bs,MinAO,MinAO,s,nbuffs,MyLsItem%setting)
         w2 = get_work_array_size(2,nv,no,nv,no,vs,os,vs,os,nb,bs,MinAO,MinAO,s,nbuffs,MyLsItem%setting)
         nelms_int2 = w1 + w2 + (nbuffs*i8*os*vs)*os*vs

         s      = 3
         nbuffs = get_nbuffs_scheme_0()
         w1 = get_work_array_size(1,nv,no,nv,no,vs,os,vs,os,nb,bs,MinAO,MinAO,s,nbuffs,MyLsItem%setting)
         w2 = get_work_array_size(2,nv,no,nv,no,vs,os,vs,os,nb,bs,MinAO,MinAO,s,nbuffs,MyLsItem%setting)
         nelms_int3 = w1 + w2 

         mem_in = min(nelms_int0,nelms_int1,nelms_int2,nelms_int3) * 8

         nelms_res_in  = 0
         nelms_res_out = 0

         if(JOB == SOLVE_AMPLITUDES .and. ccmodel == MODEL_CCSD)then     

            ! Find memory requirements in CCSD residual WITHOUT the space needed for the matrices, that use the background buffer
            mem41=int(get_min_mem_req(no,os,nv,vs,nb,0,MinAO,MinAO,0,5,4,.false.,MyLsItem%setting,'RRRRC')*1024.0E0_realk**3,kind=8)
            mem31=int(get_min_mem_req(no,os,nv,vs,nb,0,MinAO,MinAO,0,5,3,.false.,MyLsItem%setting,'RRRRC')*1024.0E0_realk**3,kind=8)
            mem21=int(get_min_mem_req(no,os,nv,vs,nb,0,MinAO,MinAO,0,5,2,.false.,MyLsItem%setting,'RRRRC')*1024.0E0_realk**3,kind=8)
            !find minimum requirements for all that IS allocated on the BG buffer
            mem42=int(get_min_mem_req(no,os,nv,vs,nb,0,MinAO,MinAO,0,8,4,.false.,MyLsItem%setting,'RRRRC')*1024.0E0_realk**3,kind=8)
            mem32=int(get_min_mem_req(no,os,nv,vs,nb,0,MinAO,MinAO,0,8,3,.false.,MyLsItem%setting,'RRRRC')*1024.0E0_realk**3,kind=8)
            mem22=int(get_min_mem_req(no,os,nv,vs,nb,0,MinAO,MinAO,0,8,2,.false.,MyLsItem%setting,'RRRRC')*1024.0E0_realk**3,kind=8)

            if(DECinfo%PL>2)then
               write( *,'("INFO(ccdriver_set_tensor_segments_and_alloc_workspace): Free         : ",g9.2," GB")')&
                  &Freebytes/1024.0E0_realk**3
               write( *,'("INFO(ccdriver_set_tensor_segments_and_alloc_workspace): scheme 4 out : ",g9.2," GB, in :",g9.2," GB")')&
                  &(mem41)/1024.0E0_realk**3,(mem42)/1024.0E0_realk**3
               write( *,'("INFO(ccdriver_set_tensor_segments_and_alloc_workspace): scheme 3 out : ",g9.2," GB, in :",g9.2," GB")')&
                  &(mem31)/1024.0E0_realk**3,(mem32)/1024.0E0_realk**3
               write( *,'("INFO(ccdriver_set_tensor_segments_and_alloc_workspace): scheme 2 out : ",g9.2," GB, in :",g9.2," GB")')&
                  &(mem21)/1024.0E0_realk**3,(mem22)/1024.0E0_realk**3
            endif

            if(mem41+mem42<=Freebytes)then
               sch   = 4
            else if(mem31+mem32<=Freebytes)then
               sch   = 3
            else if(mem21+mem22<=Freebytes)then
               sch   = 2
            else
               call lsquit("ERROR(ccdriver_set_tensor_segments_and_alloc_workspace):&
                  & calculation is too big for the available memory",-1)
            endif

            if(DECinfo%force_scheme)then
               sch=DECinfo%en_mem
            endif

            if(DECinfo%PL>2)then
               if(DECinfo%force_scheme)then
                  write( *,'("INFO(ccdriver_set_tensor_segments_and_alloc_workspace): Force scheme :",I2)')sch
               else
                  write( *,'("INFO(ccdriver_set_tensor_segments_and_alloc_workspace): Found scheme :",I2)')sch
               endif
            endif


            !simple increase buffer size
            BufferAdaption: do buf = MinAO, nb
               mem_o=int(get_min_mem_req(no,os,nv,vs,nb,0,buf,2*buf,0,5,sch,.false.,&
                  &MyLsItem%setting,'RRRRC')*1024.0E0_realk**3,kind=8)
               mem_i=int(get_min_mem_req(no,os,nv,vs,nb,0,buf,2*buf,0,8,sch,.false.,&
                  &MyLsItem%setting,'RRRRC')*1024.0E0_realk**3,kind=8)
               if( mem_o + mem_i > Freebytes)then
                  mem_o=int(get_min_mem_req(no,os,nv,vs,nb,0,buf-1,2*(buf-1),0,5,&
                     &sch,.false.,MyLsItem%setting,'RRRRC')*1024.0E0_realk**3,kind=8)
                  mem_i=int(get_min_mem_req(no,os,nv,vs,nb,0,buf-1,2*(buf-1),0,8,&
                     &sch,.false.,MyLsItem%setting,'RRRRC')*1024.0E0_realk**3,kind=8)

                  exit BufferAdaption

               endif
            enddo BufferAdaption

            nelms_res_out = ceiling(float(mem_o)/8.0)
            nelms_res_in  = ceiling(float(mem_i)/8.0)

         else if(JOB == SOLVE_AMPLITUDES .and. ccmodel == MODEL_MP2)then     


            call lsquit("ERROR(ccdriver_set_tensor_segments_and_alloc_workspace):&
               & background buffering not implemented for MP2 and amplitudes calculation",-1)

         else

            call lsquit("ERROR(ccdriver_set_tensor_segments_and_alloc_workspace):&
               & background buffering not implemented for your model/job combination",-1)

         endif

         mem_out = mem_out + nelms_res_out*8
         mem_in  = max(mem_in,nelms_res_in*8)

         bytes = mem_out

         !we may use 80% of the difference between bytes and Freebytes now as the background buffer or 2*mem_in
         bytes_to_alloc = max(min( (Freebytes - bytes) * 8 / 10, 5 * mem_in +mem_in_basic, &
            &int(DECinfo%bg_memory * 1024.0E0**3,kind=8)),mem_in+mem_in_basic)

         if(DECinfo%PL>1)then

            write (DECinfo%output,*) "cc driver info:"
            write (DECinfo%output,*) "---------------"
            write (DECinfo%output,'("Free memory is                   : ",g9.2," GB")')Freebytes/1024.0E0_realk**3
            write (DECinfo%output,'("Requirement for the solver   (o) : ",g9.2," GB")')mem_out_basic/1024.0E0_realk**3
            write (DECinfo%output,'("Requirement for the solver   (i) : ",g9.2," GB")')mem_in_basic /1024.0E0_realk**3
            write (DECinfo%output,'("Requirement for int trafo    (i) : ",g9.2," GB")')&
               &(min(nelms_int0,nelms_int1,nelms_int2,nelms_int3)*8)/1024.0E0_realk**3
            write (DECinfo%output,'("Requirement for the residual (o) : ",g9.2," GB")')(nelms_res_out*8)/1024.0E0_realk**3
            write (DECinfo%output,'("Requirement for the residual (i) : ",g9.2," GB")')(nelms_res_in*8)/1024.0E0_realk**3

            write (DECinfo%output,'("Min heap memory requirement (=o) : ",g9.2," GB")') mem_out/1024.0E0_realk**3
            write (DECinfo%output,'("Min buffer requirement      (=i) : ",g9.2," GB")') (mem_in+mem_in_basic)/1024.0E0_realk**3

            write (DECinfo%output,'("Requesting                       : ",g9.2," GB to be allocated in BG buffer")')&
               &bytes_to_alloc/1024.0E0_realk**3

#ifdef COMPILER_UNDERSTANDS_FORTRAN_2003
           flush(DECinfo%output)
#endif

         endif

         if(bytes_to_alloc <= 0) then
            print *,Freebytes/1024.0**3,"GB free, need to alloc",bytes_to_alloc/1024.0**3,"GB"
            call lsquit("ERROR(ccdriver_set_tensor_segments_and_alloc_workspace):&
               & not enough space",-1)
         endif
      endif

      if( bg_was_init )then
         if( .not. local )then
#ifdef VAR_MPI
            call lspdm_init_global_buffer(.true.)
#endif
         endif
      else
         if( local )then
            call mem_init_background_alloc(bytes_to_alloc)
#ifdef VAR_MPI
         else
            call mem_init_background_alloc_all_nodes(infpar%lg_comm,bytes_to_alloc)
            call lspdm_init_global_buffer(.true.)
#endif
         endif
      endif

   endif

end subroutine ccdriver_set_tensor_segments_and_alloc_workspace
subroutine ccdriver_dealloc_workspace(saferun,local,bg_was_init)
   implicit none
   logical, intent(in)  :: saferun,local,bg_was_init
   integer :: ntpm(4),nt
   integer(kind=ls_mpik) :: nnod
   real(realk) :: bytes, Freebytes,bytes_to_alloc,mem4,mem3,mem2
   logical :: use_bg
   integer(kind=8) :: elm
   integer :: MinAO,iAO

   nnod = 1
#ifdef VAR_MPI
   nnod = infpar%lg_nodtot
#endif

   ! deallocate the buffer in the background
   use_bg = DECinfo%use_bg_buffer.and..not.saferun

   if(use_bg)then

      if (bg_was_init)then
#ifdef VAR_MPI
         if(.not. local)then
            call lspdm_free_global_buffer(.true.)
         endif
#endif
      else
         if( local )then
            call mem_free_background_alloc()
#ifdef VAR_MPI
         else
            call mem_free_background_alloc_all_nodes(infpar%lg_comm)
            call lspdm_free_global_buffer(.true.)
#endif
         endif
      endif
   endif

end subroutine ccdriver_dealloc_workspace


function get_iter_idx(iter) result(iter_idx)
   implicit none
   integer, intent(in) :: iter
   integer :: iter_idx
   iter_idx = mod(iter-1, DECinfo%ccMaxDIIS)+1
end function get_iter_idx

subroutine ccsolver_get_residual(ccmodel,JOB,omega2,t2,&
      & fock,t1fock,iajb,no,nv,ppfock_prec,qqfock_prec,xo,xv,yo,yv,nb,&
      & MyLsItem,omega1,t1,pgmo_diag,pgmo_up,MOinfo,mo_ccsd,&
      & pno_cv,pno_s,nspaces, iter,local,use_pnos,restart,frag,m1,m2,m3,m4)
   implicit none
   integer,intent(in) :: ccmodel, JOB, nb,no,nv, nspaces, iter
   type(lsitem),intent(inout) :: MyLsItem
   type(tensor) :: omega2(:),t2(:),fock,t1fock,iajb,xo,xv,yo,yv
   type(tensor) :: omega1(:),t1(:),pgmo_diag,pgmo_up,ppfock_prec,qqfock_prec
   type(MObatchInfo),intent(inout) :: MOinfo
   type(tensor), intent(inout), optional :: m1,m2,m3,m4
   logical, intent(in) :: mo_ccsd, use_pnos, local
   logical, intent(inout) :: restart
   type(PNOSpaceInfo), intent(in), pointer :: pno_cv(:), pno_S(:)
   type(decfrag), intent(inout), optional    :: frag
   integer :: use_i
   !FIXME: remove these by implementing a massively parallel version of the
   !multipliers residual
   type(tensor) :: o2,tl2,ml4, ppfock, qqfock,pqfock,qpfock
   ! readme : get residuals, so far this solver supports only singles and doubles
   !          amplitudes (extension to higher iterative model is trivial); differences
   !          are mainly based on the set of residual vectors to be evaluated

   use_i = get_iter_idx(iter)

   ! MODIFY FOR NEW MODEL
   ! If you implement a new model, please insert call to your own residual routine here!
   SelectCoupledClusterModel : select case( ccmodel )
   case( MODEL_MP2 )

      select case(JOB)
      case( SOLVE_AMPLITUDES )
         call get_simple_parallel_mp2_residual(omega2(use_i),&
            &iajb,t2(use_i),ppfock_prec,qqfock_prec,iter,local)
      case default
         call lsquit("ERROR(ccsolver_get_residual): job not implemented for MP2",-1)
      end select

   case( MODEL_RIMP2,MODEL_LSTHCRIMP2 )

      call lsquit('ERROR(get_residual): RI-MP2 has no residual - non iterative',-1)

   case( MODEL_CC2, MODEL_CCSD ) !CC2 or  CCSD

      select case(JOB)
      case( SOLVE_AMPLITUDES, SOLVE_AMPLITUDES_PNO )

         call ccsd_residual_wrapper(ccmodel,omega2(use_i),t2(use_i),&
            & fock,t1fock,iajb,no,nv,xo,xv,yo,yv,nb,&
            & MyLsItem,omega1(use_i),t1(use_i),pgmo_diag,pgmo_up,MOinfo,mo_ccsd,&
            & pno_cv,pno_s,nspaces, iter,local,use_pnos,restart,frag=frag)


#ifdef MOD_UNRELEASED
      case( SOLVE_MULTIPLIERS )

         if (ccmodel == MODEL_CC2)then
            call lsquit("ERROR(ccsolver_get_residual): CC2 multipliers not implemented",-1)
         endif

         if(.not. local)then
            print *,"WARNING(ccsolver_get_residual): the CCSD multiplier residual is not yet parallelized"
         endif

         !FIXME: remove this dirty workaround due to the noddy implementation
         call tensor_minit(o2,[nv,no,nv,no],4)
         call tensor_zero(o2)
         call tensor_minit(tl2,[nv,no,nv,no],4)
         call tensor_cp_data(t2(use_i),tl2,order=[1,3,2,4])
         call tensor_minit(ml4,[nv,no,nv,no],4)
         call tensor_cp_data(m4,ml4)

         if(DECinfo%simple_multipler_residual)then
            call get_ccsd_multipliers_simple(omega1(use_i)%elm2,o2%elm4,m2%elm2&
               &,ml4%elm4,t1(use_i)%elm2,tl2%elm4,t1fock%elm2,xo%elm2,yo%elm2,xv%elm2,yv%elm2&
               &,no,nv,nb,MyLsItem)
         else
            print *,"The multipliers integral direct scheme is not implemented yet! Quitting..."
            stop 0
         endif

         call tensor_cp_data(o2,omega2(use_i),order=[1,3,2,4])
         call tensor_free(o2)
         call tensor_free(tl2)
         call tensor_free(ml4)
#endif

      case default
         call lsquit("ERROR(ccsolver_get_residual): job not implemented for CC2, CCSD or CCSD(T)",-1)
      end select

   case( MODEL_RPA,MODEL_SOSEX )

      call RPA_residual_par(omega2(use_i),t2(use_i),iajb,ppfock_prec,qqfock_prec,no,nv,local)

   case default

      call lsquit("ERROR(ccsolver_par):wrong choice of ccmodel",DECinfo%output)

   end select SelectCoupledClusterModel
end subroutine ccsolver_get_residual

subroutine ccsolver_calculate_crop_matrix(B,nSS,omega2,omega1,ppfock_prec,qqfock_prec,ccmodel,local,prec)
   implicit none
   real(realk), intent(out) :: B(nSS,nSS)
   integer, intent(in)      :: nSS
   type(tensor), intent(inout) :: omega2(:),omega1(:)
   type(tensor), intent(inout) :: ppfock_prec,qqfock_prec
   integer, intent(in)         :: ccmodel
   logical, intent(in)         :: local, prec
   !internal
   type(tensor) :: omega2_prec,omega1_prec
   integer :: i,j
   logical :: bg

   bg = mem_is_background_buf_init()

   ! calculate crop/diis matrix
   B=0.0E0_realk

   ! MODIFY FOR NEW MODEL
   ! If you implement a new model, please insert call to your own residual routine here!
   SelectCoupledClusterModel : select case( CCmodel )
   case( MODEL_RPA, MODEL_MP2,MODEL_SOSEX ) 

      do i=1,nSS
         do j=1,i
            ! just doubles
            if( prec ) then
               omega2_prec = precondition_doubles(omega2(j),ppfock_prec,qqfock_prec,local,bg=bg)
               B(i,j) = tensor_ddot( omega2(i), omega2_prec )
               call tensor_free( omega2_prec )
            else
               B(i,j) = tensor_ddot( omega2(i), omega2(j) )
            end if
            B(j,i) = B(i,j)
         end do
      end do
   case( MODEL_RIMP2,MODEL_LSTHCRIMP2 ) 
      call lsquit('RIMP2 not implemented in ccsolver_calculate_crop_matrix',-1)
   case( MODEL_CC2, MODEL_CCSD ) 

      do i=1,nSS
         do j=1,i
            if( prec ) then
               omega1_prec = precondition_singles( omega1(j), ppfock_prec, qqfock_prec )
               omega2_prec = precondition_doubles( omega2(j), ppfock_prec, qqfock_prec, local, bg=bg )
               B(i,j) =          tensor_ddot( omega1(i), omega1_prec ) 
               B(i,j) = B(i,j) + tensor_ddot( omega2(i), omega2_prec )

               call tensor_free( omega1_prec )
               call tensor_free( omega2_prec )
            else
               B(i,j) =          tensor_ddot( omega1(i), omega1(j) ) 
               B(i,j) = B(i,j) + tensor_ddot( omega2(i), omega2(j) )
            end if
            B(j,i) = B(i,j)
         end do
      end do

   case default
      call lsquit("ERROR(ccsolver_calculate_crop_matrix): ccmodel not known",-1)
   end select SelectCoupledClusterModel

end subroutine ccsolver_calculate_crop_matrix

subroutine ccsolver_get_special(ccmodel,mylsitem,no,nv,nb,use_pnos,mo_ccsd,Co,Cv,&
      &pgmo_diag,pgmo_up,MOinfo,nspaces,pno_cv,pno_S,m4,frag)
   implicit none
   integer, intent(in) :: ccmodel, no,nv,nb
   type(lsitem), intent(inout) :: mylsitem
   integer, intent(out) :: nspaces
   logical, intent(in) :: use_pnos
   logical, intent(out) :: mo_ccsd
   type(tensor), intent(inout) :: Co, Cv, pgmo_up, pgmo_diag
   type(MObatchInfo),intent(inout) :: MOinfo
   type(PNOSpaceInfo), intent(inout), pointer :: pno_cv(:), pno_S(:)
   type(tensor), optional, intent(inout) :: m4
   type(decfrag), optional :: frag
   type(tensor) :: tmp
   logical :: fragment_job
   real(realk) :: time_mo_ints, time_pno_setup

   fragment_job = present(frag)

   mo_ccsd = .true.
   if (DECinfo%NO_MO_CCSD.or.(no+nv>MAX_ORB_MOCCSD).or.use_pnos.or.(ccmodel==MODEL_MP2) &
      & .or. (ccmodel==MODEL_RPA).or.(ccmodel==MODEL_SOSEX) ) mo_ccsd = .false.

   if (DECinfo%force_scheme) then
      if (DECinfo%en_mem<5) then
         DECinfo%NO_MO_CCSD = .true.
         mo_ccsd            = .false.
      else if (DECinfo%en_mem>=5) then 
         mo_ccsd            = .true.
         if (DECinfo%NO_MO_CCSD) call lsquit('ERROR(CCSD): Inconsistent input, CCSD schemes &
            & 5 and 6 require the MO based algorithm. (Remove NO_MO_CCSD keyword)', DECinfo%output)
      end if
   end if


#ifdef MOD_UNRELEASED
   !============================================================================!
   !                          MO-CCSD initialization                            !
   !____________________________________________________________________________!
   ! Check if there is enough memory to performed an MO-CCSD calculation.
   !   YES: get full set of t1 free gmo and pack them
   !   NO:  returns mo_ccsd == .false. and switch to standard CCSD.
   if (mo_ccsd.and.(.not.ccmodel==MODEL_RPA).and.(.not.ccmodel==MODEL_MP2) &
     & .and. (.not.ccmodel==MODEL_SOSEX)) then
      if(DECinfo%PL>1)call time_start_phase( PHASE_work, twall = time_mo_ints ) 

      call get_t1_free_gmo(mo_ccsd,mylsitem,Co%elm2,Cv%elm2,pgmo_diag,pgmo_up, &
         & nb,no,nv,MOinfo)

      if(DECinfo%PL>1)call time_start_phase( PHASE_work, ttot = time_mo_ints,&
         &labelttot = 'CCSOL: INIT MO INTS   :', output = DECinfo%output )

   end if
#endif

   nspaces = 0
   set_pno_info:if(use_pnos)then

      if(DECinfo%PL>1)call time_start_phase( PHASE_work, twall = time_pno_setup ) 

      !FIXME: do the PNO construction in MPI parallel
      call tensor_init(tmp,[nv,no,nv,no],4)
      call tensor_convert(m4, tmp%elm1 )

      !GET THE PNO TRANSFORMATION MATRICES
      if(fragment_job)then

         !COUNT PAIRS OUTSIDE EOS
         nspaces = ( no - frag%noccEOS ) * ( no - frag%noccEOS + 1) / 2 &
            !COUNT PAIRS WITH 1 IDX IN EOS          !EOS
         &+ frag%noccEOS * ( no - frag%noccEOS ) + 1

         frag%nspaces = nspaces

         call mem_alloc( frag%CLocPNO, nspaces )
         call get_pno_trafo_matrices(no,nv,nb,tmp%elm1,frag%CLocPNO,frag%nspaces,f=frag)
         pno_cv => frag%CLocPNO

      else
         !ALL PAIRS
         nspaces = no * ( no + 1 ) / 2
         call mem_alloc( pno_cv, nspaces )
         call get_pno_trafo_matrices(no,nv,nb,tmp%elm1,pno_cv,nspaces,f=frag)

      endif

      call tensor_free(tmp)

      !GET THE OVERLAP BETWEEN THE PNO SPACES
      call mem_alloc( pno_S , nspaces * (nspaces - 1)/2 )   
      !Get all the overlap matrices necessary
      call get_pno_overlap_matrices(no,nv,pno_cv,pno_S,nspaces,.true.)

      if(DECinfo%PL>1)call time_start_phase( PHASE_work, ttot = time_pno_setup,&
         &labelttot = 'CCSOL: PNO SPACE SETUP:', output = DECinfo%output )
   endif set_pno_info

end subroutine ccsolver_get_special

subroutine ccsolver_free_special(pgmo_up,pgmo_diag,MOinfo,restart,mo_ccsd,pno_cv,pno_S,nspaces,use_pnos,frag)
   implicit none
   type(tensor), intent(inout) :: pgmo_up,pgmo_diag
   type(MObatchInfo),intent(inout) :: MOinfo
   type(PNOSpaceInfo), intent(inout), pointer :: pno_cv(:), pno_S(:)
   logical, intent(in) :: restart, mo_ccsd, use_pnos
   integer, intent(in) :: nspaces
   type(decfrag), optional :: frag
   integer :: i,j,cc
   logical :: fragment_job
#ifdef MOD_UNRELEASED

   fragment_job = present(frag)

   ! free memory from MO-based CCSD
   if(.not. restart)then
      if (mo_ccsd) then
         if (pgmo_diag%dims(2)>1) call tensor_free(pgmo_up)
         call tensor_free(pgmo_diag)
         call mem_dealloc(MOinfo%dimInd1)
         call mem_dealloc(MOinfo%dimInd2)
         call mem_dealloc(MOinfo%StartInd1)
         call mem_dealloc(MOinfo%StartInd2)
         call mem_dealloc(MOinfo%dimTot)
         call mem_dealloc(MOinfo%tileInd)
      end if
   endif

   !Free PNO information
   if(use_pnos)then

      if(.not.fragment_job)then

         do i = 1, nspaces

            if( pno_cv(i)%allocd )then
               call free_PNOSpaceInfo(pno_cv(i))
            endif

            do j = 1, i - 1
               cc = (j - i + 1) + i*(i-1)/2
               if( pno_S(cc)%allocd )  call free_PNOSpaceInfo( pno_S(cc) )
            enddo
         enddo

         call mem_dealloc( pno_cv )

      else

         do i = 1, nspaces

            if( frag%CLocPNO(i)%allocd )then
               call free_PNOSpaceInfo( frag%CLocPNO(i) )
            endif

            do j = 1, i - 1
               cc = (j - i + 1) + i*(i-1)/2
               if( pno_S(cc)%allocd )  call free_PNOSpaceInfo( pno_S(cc) )
            enddo
         enddo

         call mem_dealloc( frag%CLocPNO )
         pno_cv => null()

      endif

      call mem_dealloc( pno_S )

   endif


   if( .not. fragment_job .and. DECinfo%PL>2 )then
      call tensor_print_mem_info(DECinfo%output,.true.,.false.)
   endif
#endif
end subroutine ccsolver_free_special

!> \brief should be a general subroutine to get the guess amplitudes when
!starting up a CCSD or CC2 calculation and checks for files which contain
!amplitudes of a previous calculation. If none are found the usual zero guess
!is returned
!> \author Patrick Ettenhuber
!> \date December 2012
subroutine get_guess_vectors(ccmodel,JOB,prec,restart,iter_start,nb,norm,energy,t2,iajb,Co,Cv,Uo,Uv,oof,vvf,vof,ovf,mylsitem,local,&
   & safefilet21,safefilet22,safefilet2f, t1,safefilet11,safefilet12,safefilet1f)
   implicit none
   integer, intent(in) :: nb,ccmodel,JOB
   logical,intent(out) :: restart
   real(realk), intent(inout) :: norm,energy,Uo(:,:),Uv(:,:)
   !> contains the guess doubles amplitudes on output
   type(tensor), intent(inout) :: t2,iajb,Co,Cv,oof,vvf,vof,ovf
   logical, intent(in) :: local,prec
   !> integral info
   type(lsitem), intent(inout) :: mylsitem
   !> the filenames to check for valid doubles amplitudes
   character(3),intent(in) :: safefilet21,safefilet22,safefilet2f
   !> contains the singles amplitudes on output
   type(tensor),intent(inout),optional :: t1
   !> the filenames to check for valid singles amplitudes
   character(3),intent(in), optional :: safefilet11,safefilet12,safefilet1f
   integer, intent(out) :: iter_start
   integer :: no,nv
   integer(8) :: fu_t11,fu_t12,fu_t21,fu_t22,fu_t1,fu_t2,fu_t2f,fu_t1f
   logical :: file_exists11,file_exists12,file_exists1f,file_exists21,file_exists22,file_exists2f
   logical(8) :: file_status11,file_status12,file_status1f,file_status21,file_status22,file_status2f
   logical(8) :: readfile1, readfile2
   integer(8) :: saved_iter11,saved_iter12,saved_iter1f,saved_iter21,saved_iter22,saved_iter2f
   integer :: saved_nel11,saved_nel12,saved_nel21,saved_nel22,saved_nel1f,saved_nel2f
   logical :: all_singles, fin1_exists, fin2_exists
   character(11) :: fullname11, fullname12, fin1, fullname21, fullname22,fin2
   character(tensor_MSG_LEN) :: msg
   integer :: a,i
   logical :: use_singles, get_multipliers, use_bg

   all_singles=present(t1).and.present(safefilet11).and.present(safefilet12).and.present(safefilet1f)
   get_multipliers = (JOB == SOLVE_MULTIPLIERS)
   use_bg = (mem_is_background_buf_init()) .and. (mem_get_bg_buf_free() > t2%nelms)

   !MODIFY FOR NEW MODEL THAT IS USING THE CC DRIVER
   ! set model specifics here
   select case(ccmodel)
   case(MODEL_MP2)
      use_singles = .false.
   case(MODEL_RIMP2)
      use_singles = .false.
   case(MODEL_LSTHCRIMP2)
      use_singles = .false.
   case(MODEL_CC2)
      use_singles = .true.
   case(MODEL_CCSD)
      use_singles = .true.
   case(MODEL_RPA)
      use_singles = .false.
   case(MODEL_SOSEX)
      use_singles = .false.
   case default
      call lsquit("ERROR(get_guess_vectors) unknown model",-1)
   end select

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
      if(use_singles.and.all_singles)then
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
         iter_start = int(saved_iter2f)
      endif

      !THEN CHECK IF THERE ARE AMPLITUDES FROM OTHER ITERATIONS AVAILALBE
      if(use_singles.and.all_singles.and..not.readfile1)then
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
               iter_start=int(saved_iter21)
               fu_t2=fu_t21
               CLOSE(fu_t22)
               readfile2=.true.
            else
               iter_start=int(saved_iter22)
               fu_t2=fu_t22
               CLOSE(fu_t21)
               readfile2=.true.
            endif
            WRITE(DECinfo%output,'("RESTARTING CC CALCULATION WITH TRIAL-VECS FROM: ",I3)')iter_start
         else if(file_status21)then
            iter_start=int(saved_iter21)
            fu_t2=fu_t21
            WRITE(DECinfo%output,'("RESTARTING CC CALCULATION WITH TRIAL-VECS FROM: ",I3)')iter_start
            readfile2=.true.
         else if(file_status22)then
            iter_start=int(saved_iter22)
            fu_t2=fu_t22
            WRITE(DECinfo%output,'("RESTARTING CC CALCULATION WITH TRIAL-VECS FROM: ",I3)')iter_start
            readfile2=.true.
         else
            iter_start=0
         endif  
      endif

   endif


   if(use_singles)then
      if(readfile1)then
         READ(fu_t1) t1%elm1
         READ(fu_t1) norm
         READ(fu_t1) energy
         CLOSE(fu_t1)
         restart = .true.
         call local_can_trans(no,nv,nb,Uo,Uv,vo=t1%elm1)
      else
         !This was an attempt to start with something else than 0 for CCSD, but
         !it does not seem worth it, even if the start is the amplitudes of the
         !second iteration
         if(get_multipliers)then
            select case(ccmodel)
            case(MODEL_CCSD)
               call array_reorder_2d(-2.0E0_realk,ovf%elm2,no,nv,[2,1],0.0E0_realk,t1%elm2)
            case default
               call tensor_zero(t1)
            end select
         else
            select case(ccmodel)
            case(MODEL_CC2,MODEL_CCSD)
               call tensor_cp_data(vof,t1)
            case default
               call tensor_zero(t1)
            end select

            if(prec)then
               do a=1,nv
                  do i=1,no

                     t1%elm2(a,i) = t1%elm2(a,i)/( oof%elm2(i,i) - vvf%elm2(a,a) ) 

                  end do
               end do
            endif

         endif
      endif
   endif

   if(readfile2)then
      ! allocate dense part of t2 array:
      if (.not.local) call tensor_allocate_dense(t2,bg=use_bg)
      READ(fu_t2) t2%elm1
      READ(fu_t2) norm
      READ(fu_t2) energy
      CLOSE(fu_t2)
      ! mv dense part to tiles:
      call local_can_trans(no,nv,nb,Uo,Uv,vvoo=t2%elm1)
      if (.not.local) call tensor_mv_dense2tiled(t2,.false.)
      restart = .true.
   else
      if(get_multipliers)then
         select case(ccmodel)
         case(MODEL_CCSD)
            call get_starting_guess( iajb, t2, oof, vvf, local, "CCSD_LAG_RHS", prec )
         case default
            call tensor_zero(t2)
         end select
      else
         select case(ccmodel)
         case(MODEL_MP2,MODEL_CC2,MODEL_CCSD,MODEL_CCSDpT)
            call get_starting_guess( iajb, t2, oof, vvf, local, "MP2AMP", prec )
         case default
            call tensor_zero(t2)
         end select
      endif
   endif
end subroutine get_guess_vectors

!> \brief Subroutine to save the current guess amplitudes for the next
!iteration
!> \author Patrick Ettenhuber
!> \date Dezember 2012
subroutine save_current_guess(local,iter,nb,res_norm,energy,Uo,Uv,t2,safefilet21,safefilet22,&
      &t1,safefilet11,safefilet12)
   implicit none
   logical, intent(in) :: local
   !> iteration number
   integer,intent(in) :: iter,nb
   !> write the corresponding residual norm into the file
   real(realk), intent(in)    :: res_norm,energy,Uo(:,:),Uv(:,:)
   !> doubles guess amplitudes for the next iteration
   type(tensor), intent(inout) :: t2
   !> alternating filenames for the doubles amplitudes
   character(3),intent(in)    :: safefilet21,safefilet22
   !> singles guess amplitudes for the next iteration
   type(tensor), intent(inout), optional :: t1
   !> alternating filenames for the singles amplitudes
   character(3),intent(in), optional :: safefilet11,safefilet12
   integer :: fu_t21,fu_t22
   integer :: fu_t11,fu_t12
   integer :: no, nv
   logical(8) :: file_status11,file_status12,file_status21,file_status22
   logical :: all_singles
   character(tensor_MSG_LEN) :: msg
#ifdef SYS_AIX
   character(12) :: fullname11,  fullname12,  fullname21,  fullname22
   character(12) :: fullname11D, fullname12D, fullname21D, fullname22D
#else
   character(11) :: fullname11, fullname12, fullname21, fullname22
   character(11) :: fullname11D, fullname12D, fullname21D, fullname22D
#endif

   nv = t2%dims(1) 
   no = t2%dims(3) 

   ! cp doubles from tile to dense part: (only if t2%itype/=TT_DENSE)
   if (.not.local) call tensor_cp_tiled2dense(t2,.false.)

   call can_local_trans(no,nv,nb,Uo,Uv,vvoo=t2%elm1)

   all_singles=present(t1).and.present(safefilet11).and.present(safefilet12)
   fu_t11=111
   fu_t12=112
   fu_t21=121
   fu_t22=122

   if(DECinfo%use_singles.and.all_singles)then

      call can_local_trans(no,nv,nb,Uo,Uv,vo=t1%elm1)

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
         WRITE(fu_t11)int(iter,kind=8)
         WRITE(fu_t11)int(t1%nelms,kind=8)
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
         WRITE(fu_t12)int(iter,kind=8)
         WRITE(fu_t12)int(t1%nelms,kind=8)
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

      call local_can_trans(no,nv,nb,Uo,Uv,vo=t1%elm1)

   endif

   fullname21=safefilet21//'.writing'
   fullname22=safefilet22//'.writing'
   if(mod(iter,2)==1)then
      file_status21=.false. 
      OPEN(fu_t21,FILE=fullname21,STATUS='REPLACE',FORM='UNFORMATTED')
      WRITE(fu_t21)int(iter,kind=8)
      WRITE(fu_t21)int(t2%nelms,kind=8)
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
      WRITE(fu_t22)int(iter,kind=8)
      WRITE(fu_t22)int(t2%nelms,kind=8)
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

   call local_can_trans(no,nv,nb,Uo,Uv,vvoo=t2%elm1)

   ! deallocate dense part of doubles:
   if (.not.local) call memory_deallocate_tensor_dense(t2)
end subroutine save_current_guess


end module ccdriver
