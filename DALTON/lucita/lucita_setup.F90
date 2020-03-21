!dalton_copyright_start
!
!
!dalton_copyright_end

module lucita_setup

! stefan: interface to lucita setup routines which need to be called for each 
!         mcscf/ci runtype

  implicit none

  public setup_lucita_orbital_string_cb
  public setup_lucita_pointer_strings_work_space
  public setup_lucita_initialize_work_space
  public setup_lucita_cb_interface
  public setup_lucita_check_input
  public setup_lucita_mcci_wrkspc_dimensions
  public setup_lucita_par_dist_in_seq
#ifdef VAR_MPI
  public setup_lucita_inc_wrkspc_alloc_cw
#endif

  private

#include "priunit.h"

contains

  subroutine setup_lucita_orbital_string_cb(set_common_blocks,              &
                                            print_level)
!*******************************************************************************
!
!    purpose:  interface routine to LUCITA setup routines which 
!              will create the mandatory orbital, symmetry and string
!              information for LUCITA ci/mcscf calculations.
!
!*******************************************************************************
    logical, intent(inout) :: set_common_blocks
    integer, intent(in)    :: print_level
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      select case(set_common_blocks)
        case (.true.)
!         translate shell spaces into orbital spaces + create lucita/dalton reorder and index arrays
          call orbinf(print_level)
!         calculate the number of string types
          call strtyp_gas(print_level)
!         divide orbital spaces into inactive/active/secondary
          call gasspc()
!         setup symmetry information 
          call syminf(print_level)
!         calculate the number of 1e- and 2e-integrals
          call intdim(print_level)
          call flshfo(lupri)
        case (.false.)
!         do nothing...
      end select

  end subroutine setup_lucita_orbital_string_cb
!*******************************************************************************

  subroutine setup_lucita_pointer_strings_work_space(ci_space,              &
                                                     current_symmetry,      &
                                                     print_level)
!*******************************************************************************
!
!    purpose:  interface routine to LUCITA setup routines which 
!              will create the mandatory work space pointers, string
!              information for LUCITA ci/mcscf calculations on work space blocks.
!
!*******************************************************************************
    integer, intent(in)    :: ci_space
    integer, intent(in)    :: current_symmetry
    integer, intent(in)    :: print_level
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!     allocate local and global arrays ==> set work space pointers
      call lucita_alloc()
!     internal string information
      call strinf_gas(print_level)
!     internal subspaces
      call lcispc(ci_space,                &
                  current_symmetry,        &
                  print_level)

  end subroutine setup_lucita_pointer_strings_work_space
!*******************************************************************************

  subroutine setup_lucita_initialize_work_space(work_space_in,                 &
                                                work_space_out,                &
                                                work_space_in_length,          & 
                                                work_space_out_length,         &
                                                luci_nmproc,                   &
                                                nr_files,                      &
                                                print_unit)
!*******************************************************************************
!
!    purpose:  interface routine to LUCITA setup routines which 
!              will create the mandatory work space pointers, string
!              information for LUCITA ci/mcscf calculations on work space blocks
!              as well as read in integrals (if applicable).
!
!*******************************************************************************
    real(8), intent(inout) :: work_space_in(*)
    real(8), intent(inout) :: work_space_out(*)
    integer, intent(in   ) :: work_space_in_length
    integer, intent(  out) :: work_space_out_length
    integer, intent(in)    :: luci_nmproc
    integer, intent(in)    :: nr_files
    integer, intent(in)    :: print_unit
!-------------------------------------------------------------------------------
    integer, parameter     :: max_num_ttss_blocks = 20000
    integer, parameter     :: max_subspace_dim    = 300
    integer                :: k1
    integer                :: kbase_lucita
    integer                :: dummy
    integer(8)             :: offset_2_work_space_out
    integer(8)             :: k8base_lucita
    integer(8)             :: work_space_out_length_scratch
    logical, save          :: first_time_print = .true.
!-------------------------------------------------------------------------------

!     compute offset to LUCITA work space
      call compoff(work_space_in,work_space_out,offset_2_work_space_out)

#ifdef integer_test
      print *,'K_OFFSET      =',offset_2_work_space_out
      ! integer   :: kbase_lucita
      ! integer*8 :: k8base_lucita
      KBASE_LUCITA = offset_2_work_space_out
      ! HJAaJ Aug 2008:
      ! MEMMAN in LUCITA is using default integer type
      ! while K_OFFSET is always INTEGER*8.
      ! Therefore we use implicit conversion to KBASE_LUCITA.
      ! Now we test if KBASE_LUCITA is OK:
      K8BASE_LUCITA = KBASE_LUCITA
      print *,'KBASE_LUCITA  =',KBASE_LUCITA
      print *,'K8BASE_LUCITA =',K8BASE_LUCITA
      IF(K8BASE_LUCITA /= offset_2_work_space_out)THEN
         WRITE (print_unit,'(/A/A,2I20//A)')                                   &
         'FATAL ERROR in LUCITA: offset is too big for i*4',                   &
         'K_OFFSET and KBASE_LUCITA:',offset_2_work_space_out, kbase_lucita,   &
         '(See code in lucita/lucita_setup.F90)'
         call quit('FATAL ERROR in LUCITA: offset is too big for i*4')
      END IF
#endif
! 
!     subtract dynamically allocated memory
!     the parameter values are based on an educated guess from my experience - stefan nov 2010
      work_space_out_length_scratch = work_space_in_length                                        &
                                    - 2 * luci_nmproc - nr_files -  8 * max_num_ttss_blocks       &
                                    - (max_subspace_dim*(max_subspace_dim+1)/2)                   &
                                    - 2 * max_subspace_dim**2

      work_space_out_length         = work_space_out_length_scratch

      if(first_time_print)then
        write(print_unit,'(/A,I18)') '          dimension of LUCITA work space : ', work_space_out_length
        first_time_print = .false.
      end if

!     initialize LUCITA internal memory
      call memman(offset_2_work_space_out,work_space_out_length_scratch,'INI   ',dummy,'DUMMY')

  end subroutine setup_lucita_initialize_work_space
!*******************************************************************************

  subroutine setup_lucita_cb_interface()
!
!    purpose: interface routine for general purpose lucita common blocks
!             the ci_run_id controls the redefinition of some variables 
!             according to:
!             - CI calculation in general
!             - sigma-vector run
!             - density-matrix evaluation (1-/2-particle density matrix)
!             - CI vector analysis
!             - calculation of H_diagonal
!
!
!*******************************************************************************
     use lucita_cfg
     use lucita_energy_types
     use lucita_ci_task_interface, only: create_CI_task_list
#ifdef VAR_MPI
     use sync_coworkers
#endif
!-------------------------------------------------------------------------------
#include "priunit.h"
#include "mxpdim.inc"
#include "lucinp.inc"
#include "cstate.inc"
#include "crun.inc"
#include "cprnt.inc"
#include "cgas.inc"
#include "cands.inc"
#include "comjep.inc"
#include "oper.inc"
#include "orbinp.inc"
#include "intform.inc"
#include "cc_exc.inc"
! common block for parallel calculations
#include "parluci.h"
! common blocks with orbital/energy information outside the LUCITA world (== within DALTON)
#include "inforb.h"
#include "maxorb.h"
#include "infinp.h"
!*******************************************************************************
      integer :: i, j, ntoob_local
      integer :: mxelr1, mxelr2, nimx, nimn, nemn, nemx, namx, namn
      integer :: ispin1, ispin2 ! these variables are not active yet...
!-------------------------------------------------------------------------------
      
!     set internal print units (luwrt on common block in parluci.h)
      luwrt     = lupri
!     initialize point group D2H and its subgroups
      pntgrp = 1
!     default number of point group irreps in lucita
      nirrep = lucita_cfg_nr_ptg_irreps
      nsmob  = nirrep
!     orbital environment sync with Dalton
      naos_env(1:nirrep) = naos_lucita(1:nirrep)
      nmos_env(1:nirrep) = nmos_lucita(1:nirrep)
!     nuclear repulsion energy
      ecore_orig = lucita_cfg_core_energy
!     initialize variables used for other point groups
      maxml  = -1
      maxl   = -1
!     internal CI space
      intspc = 1 !???
!     # of active e-
      luci_nactel = lucita_cfg_nr_active_e
!     inactive shells
      call izero(ninash,nirrep)
!     inactive shells - 2
      call izero(ninobs,mxpobs)
!     Core shells ( = RAS0 shells) only relevant for extspace > 0
      call izero(nrs0sh,nirrep)
!     ras 1, 2, 3 and 4
      mxr4tp = 1
      mxer4  = 0
      call izero(nrssh(1,1),nirrep)
      call izero(nrssh(1,2),nirrep)
      call izero(nrssh(1,3),nirrep)
      call izero(nrs4sh(1,mxr4tp),nirrep)
!     ?
      mnrs10 = -1
      mxrs30 = -1
      mnrs1r = mnrs10
      mxrs3r = mxrs30
!     selection of internal configurations - no
      intsel = 0
      mults  = lucita_cfg_is_spin_multiplett
      ms2    = mults - 1 ! MS(MAX) = S set by default
      ispin1 = lucita_cfg_spin1_2el_op
      ispin2 = lucita_cfg_spin2_2el_op
!     reference symmetry
      irefsm_c  = lucita_cfg_csym
      irefsm_hc = lucita_cfg_hcsym
      irefsm    = irefsm_c 
!     transfer to common block cands
      icsm   = irefsm_c
      issm   = irefsm_hc
      icspc  = 1
      isspc  = 1
!     initialize variables on comjep.inc
      mxacj    = 0
      mxacij   = 0
      mxaadst  = 0
!     roots to be obtained
      nroot  = lucita_cfg_nr_roots
      do i = 1, nroot
        iroot(i) = i
      END DO
!     state of interest (if particular convergence to a given state is desired)
      c_state_of_interest = lucita_cfg_icstate
!     standard davidson
      idiag  = 1
!     no perturbative Hamiltonian
      MXP1    = 0
      MXP2    = 0
      MXQ     = 0
      IPERTOP = 0
!     max # of Davidson iterations
      MAXIT  = lucita_cfg_max_nr_dav_ci_iter
!     restart ci
      irestr = lucita_cfg_restart_ci
!     integrals imported from SIRIUS/DALTON
      INTIMP      = 5
      ENVIRO(1:6) = 'DALTON'
!     in core option for integrals - default: in core
      incore = 1
!     no deleted shells
      call izero(ndelsh,nirrep)
!     no spin/ml combination factors
      pssign = 0.0D0
      plsign = 0.0D0
      if(pssign.eq.0.0d0.and.plsign.eq.0.0d0)then
        IDC = 1
      else if(pssign.ne.0.0d0.and.plsign.eq.0.0d0)then
        IDC = 2
      else if(pssign.eq.0.0d0.and.plsign.ne.0.0d0)then
        IDC = 3
      else if(pssign.ne.0.0d0.and.plsign.ne.0.0d0)then
        IDC = 4
      end if
      plus_combi = lucita_cfg_plus_combi ! plus combination of start vectors
!     default print levels - zero for the time being - should be
!     controlled via input keyword
      IPRSTR  = 0
      IPRCIX  = IPRSTR
      IPRORB  = IPRSTR
      IPRDIA  = IPRSTR
      IPRXT   = IPRSTR
      IPRRSP  = IPRSTR
      IPRDEN  = IPRSTR
      IPROCC  = IPRSTR
      IPRNCIV = 0
      if(lucita_cfg_analyze_cvec) iprnciv = 1
!     set timing printouts for parallel calculations
      timing_par = .false.
      if(lucita_cfg_timing_par) timing_par = .true.
!     max davidson subspace
      mxciv = max(6*nroot, lucita_cfg_max_dav_subspace_dim) 
      if(idiag == 2) mxciv = max(2*nroot, lucita_cfg_max_dav_subspace_dim)
!     storage mode for vectors
!     complete vector + block version in core
!     icistr = 1
!     default is three type-type-symmetry blocks
      icistr = 3
!     no csf - determinants == default
      NOCSF  = 1
!     read in integrals
      NOINT  = 0
!     NOINT  = 1 ! for densities/analysis ...
!     do not dum integrals
      IDMPIN = 0
!     define dimension of resolution matrices
!SK12 mxinka = 100000
      mxinka = 1000
!     use CJKAIB matrices as intermediate matrices in alpha-beta-loop
      icjkaib = 1
!     use minimal operatioon count method for alpha-alpha and beta-beta
!     - default: do not use it
      MOCAA = 0
!     use minimal operatioon count method for alpha-beta
!     - default: do not use it
      MOCAB = 0
!     initial CI in reference space
!     - default: no
      INIREF = 0
!     restart from reference CI expansion
      IRESTRF = 0
!     initialize core energy
      ECORE = 0.0D0
!     do not use perturbation theory for zero order space
      IPERT = 0
      NPERT = 0
!     no approximate Hamiltonian in reference space
      IAPRREF = 0
      MNRS1RE = MNRS1R
      MXRS3RE = MXRS3R
!     no approximate Hamiltonian in zero order space
      IAPRZER = 0
      MNRS1ZE = MNRS10
      MXRS3ZE = MXRS30
!     default # of orbital space occupations
      NCISPC  = lucita_cfg_nr_calc_sequences
      NCMBSPC = NCISPC
!     general sequencer : default is just straight sequence of CI with default number of iterations
!     thus no sequence of calculations: ci, pert, cc,...
      DO i = 1, NCMBSPC
        LCMBSPC(i)   = 1
        ICMBSPC(1,i) = i
        NSEQCI(i)    = 1
        CSEQCI(1,i)  = 'CI'
        ISEQCI(1,i)  = MAXIT
      END DO
!     energy convergence of CI
!     THRES_E = 1.0D-10
      THRES_E     = lucita_cfg_convergence_c
      THRES_E_aux = lucita_cfg_auxiliary_conv_c
!     threshold on energy and wave function corrections in truncation of classes (NOT active)
      E_THRE      = 0.0d0
      C_THRE      = 0.0d0
      E_CONV      = 0.0d0
      C_CONV      = 0.0d0
!     do not call EXTENDED KOOPMANS' THEOREM ROUTINE
      IEXTKOP = 1
!     do not save first order wave function correction to disc (if icistr == 1)
      IC1DSC = 0
!     no perturbation in subspaces spaces
      IH0SPC = 0
!     reference root for Perturbation theory
      IRFROOT = NROOT
!     default handling of degenrences of initial CI vectors
!     - no action
      INIDEG = 0
!     do not use modified Hamilton operator in CI optimization
      XLAMBDA = 1.0D0
!     length of smallest block for batch of C an Sigma vectors
      LCSBLK = max(100000,lucita_cfg_max_batch_size)
!     final orbitals
      IFINMO = 1 ! natural orbitals (if at all)
!     default: no class selection unless we do teraci (which is disabled... )
      ICLSSEL = 0
      if(idiag == 2) ICLSSEL = 1

!     perturbation expansion of EKT, default is no
      IPTEKT = 0
!     root used to define zero order Hamiltonian
      IH0ROOT = NROOT
!     no restart in CI calculation 2
      IRST2 =  0
!     skip initial evaluation of energy from CI calc 2
!     ISKIPEI = 1
!     symmetry of X, Y and Z - zero because properties + transition properties are disabled
      ITRAPRP = 0
      NPROP   = 0
      DO i = 1, 3
        IXYZSYM(i) = 0
      END DO
!     no CI response
      IRESPONS = 0
      NRESP    = 0
      N_AVE_OP = 0
!     max # of iterations in linear equations (ci response)
      MXITLE      = 20
!     not root homing
      IROOTHOMING = 0
!     calculation of density matrices
      idensi  = lucita_cfg_density_calc_lvl
!     calculate spin densities
      ispnden =  lucita_cfg_spindensity_calc_lvl
      !> set level of density calculation equal to spin-density level to avoid any inconsistencies...
      if(ispnden > idensi)then 
        idensi                      = ispnden
        lucita_cfg_density_calc_lvl = ispnden
      end if
!     no particle-hole simplification in use for compatibility with densities
!     jeppe + stefan: april 2011
      iuse_ph = 0
      if(idensi > 0) IUSE_PH = 0
!     allow the sigma routine to take advice
      IADVICE = 0
      if(IUSE_PH > 0) IADVICE = 1
!     do not transform CI vectors to alternative orbital representation
      ITRACI    = 0
      ITRACI_CR = 'COMP    '
      ITRACI_CN = 'NATURAL '
!     do no separate strings into active and passive parts
      IUSE_PA = 0
!     no perturbation expansion of Fock matrix
      IPTFOCK = 0
!     no restart of CC with CI coefficients
      I_RESTRT_CC = 0
!     use relaxed densities for properties
!     default: do not use them
      IRELAX= 0
      if(IRELAX > 0) IDENSI = 2
!     external space calculation - default: no
      EXTSPC = 0
      MXHR0  = 0
!     setting of parallel I/O model
      IIOMOD    = 0
      IF(LUCI_NMPROC .GT. 1) IIOMOD = 1
!     parallel distribution routine
      idistroute = lucipar_cfg_ttss_dist_strategy
!     Truncate residual vectors before creating new trial vector?
!     (14-jun-07, hjaaj)
      trunc_fac = lucita_cfg_accepted_truncation
!     memory reduction multiplier ...
      ismemfac  = lucipar_cfg_mem_reduction_multp 
!     LUCITA is used in real calculations: no complex part
      irc_save = 1
!     normal integral accessed
      IH1FORM  = 1
      IH2FORM  = 1
      I_RES_AB = 0
!     CI not CC
      ICC_EXC  = 0
!     default: complete operator (1e- + 2e-)
      i12 = lucita_cfg_el_operator_level

      if(lucita_cfg_initialize_cb)then

!       generalized active space concept, define orbital spaces
!       -------------------------------------------------------
        ngas    = lucita_cfg_nr_gas_spaces 

        if(lucita_cfg_ci_type(1:6) == 'RASCI ')then

          ngas = 3
!         check for atypical RAS (no RAS3):
          i = 0
          do j = 1, nirrep
            i = i + ngsh_lucita(3,j)
          end do
          if(i <= 0 ) ngas = 2 
        end if

        do i = 1, ngas
          do j = 1, nirrep
            NGSSH(j,i) = ngsh_lucita(i,j)
          end do
        end do


!       check for maximum number of orbitals per space and symmetry
        do i = 1, NGAS 
          do j = 1, NIRREP
            if(NGSSH(j,i) > MXTSOB)then
              print *, ' too many orbitals per space and symmetry!'
              print *, ' my maximum is  ',MXTSOB
              print *, ' please redefine your active spaces or if you are doing RASCI please use the GA setup.'
              call quit('*** error in setup_lucita_cb_interface: too many orbitals per space and symmetry!')
            end if
          end do
        end do

!       set orbital space occupations
        select case(lucita_cfg_ci_type(1:6))
!         ras
          case('RASCI ')

            if(luci_nactel <= 0) call quit('*** # active e- <= 0... ***')

            nimx = lucita_cfg_max_e_ras1
            nimn = lucita_cfg_min_e_ras1
            nemn = LUCI_NACTEL
            nemx = LUCI_NACTEL
            namx = nemx - lucita_cfg_min_e_ras3
            namn = namx - lucita_cfg_max_e_ras3
            igsoccx(1,1,1) = nimn
            igsoccx(1,2,1) = nimx
            igsoccx(2,1,1) = namn
            igsoccx(2,2,1) = namx
!           write(lupri,*) 'nimn, nimx, namn, namx, nemx, nemn', &
!           nimn, nimx, namn, namx, nemx, nemn
            if(ngas > 2)then
              igsoccx(3,1,1) = nemn
              igsoccx(3,2,1) = nemx
            end if
            if(namx > luci_nactel ) call quit('*** reconsider your RAS setup - it is wrong...  ***')
!         gas
          case('GASCI ')
            do i = 1, ngas
              do j = 1, 2
                igsoccx(i,j,1) = ngso_lucita(i,j)
              end do
            end do
        end select

        if(LUCI_NACTEL /= igsoccx(ngas,2,1)) then
          write(lupri,*) 'Number of active electrons does not match total number of electrons in active spaces.'
          call quit('*** error in setup_lucita_cb_interface: Number of active'//                  &
                    ' electrons does not match total number of electrons in active spaces.')
        end if

        ntoob_local = 0
        do i = 1, ngas
          do j = 1, nirrep
            ntoob_local = ntoob_local + NGSSH(j,i)
          end do
        end do

        if(LUCI_NACTEL > 2*ntoob_local)then
          write(lupri,*) 'Number of active electrons exceeds the orbital space:',LUCI_NACTEL,'>',2*NTOOB_local
          write(lupri,*) 'Consider Pauli´s famous principle and restart!'
          call quit('*** error in setup_lucita_cb_interface: Number of active electrons exceeds the orbital space.')
        end if

!       set file handles for reading/writing, etc.
        call set_file_handles(irefsm)


      end if ! lucita_cfg_initialize_cb check

!     cross-check input and print LUCITA settings
      call setup_lucita_check_input(lucita_ci_run_id,igsoccx,ngssh,ngas)

!     create the ci_task_list based on the CI run-id and input parameters
      call create_CI_task_list(lucita_ci_run_id,                  &
                               iprnciv,                           &
                               idensi,                            &
                               ci_task_list,                      &
                               max_ci_tasks)

!     we are done :) - return

 end subroutine setup_lucita_cb_interface
!*******************************************************************************

  subroutine setup_lucita_check_input(ci_run_id,igsoccx,ngssh,ngas)
!*******************************************************************************
!
!    purpose:  cross-check the mandatory LUCITA input variables to tally with  
!              a reasonable setting.
!
!*******************************************************************************
    use lucita_cfg
#include "priunit.h"
#include "parluci.h"
#include "mxpdim.inc"
!-------------------------------------------------------------------------------
    character (len=12), intent(in) :: ci_run_id
    integer           , intent(in) :: igsoccx(mxpngas,2,mxpici)
    integer           , intent(in) :: ngssh(mxpirr,mxpngas)
    integer           , intent(in) :: ngas
!-------------------------------------------------------------------------------
    integer                        :: print_lvl

!
!     CI type for LUCIA (no default)
      if(lucita_cfg_ci_type(1:4).eq.'none')then
         write(lupri,*) ' Keyword for type of CI calculation missing. '
         write(lupri,*) ' This keyword is mandatory. '
         call quit(' *** error in setup_lucita_check_input: keyword .CITYPE not specified.')      
      else
        if(lucita_cfg_ci_type(1:6).ne.'GASCI '.and.  lucita_cfg_ci_type(1:6).ne.'RASCI ')then
           write(lupri,'(//A//2A/A)')                                                             &
           ' Type of CI calculation not properly specified.',                                     &
           ' You have chosen: ',lucita_cfg_ci_type,                                               &
           ' Allowed types  :  GASCI and RASCI'
          call quit(' *** error in setup_lucita_check_input: wrong input to keyword .CITYPE specified.')
        end if
      end if

!     # of active e-
      if(lucita_cfg_nr_active_e < 0)then
        write(lupri,'(//A/A/)') ' Number of active electrons NACTEL must be',                     &
               ' specified for this type of calculation. Quitting.'
          call quit(' *** error in setup_lucita_check_input: # electrons .NACTEL must be specified.')
      end if

!     spin multiplicity
      if(lucita_cfg_is_spin_multiplett < 0)then
        write(lupri,*) 'Spin multiplicity .MULTIP is a MANDATORY keyword.'
        write(lupri,*) 'Specify and restart.'
          call quit(' *** error in setup_lucita_check_input: spin multiplicity .MULTIP must be specified in any case.')
      else
!       check consistency of # of e- and spin multiplicity
!       case a: even - even
        if(mod(lucita_cfg_nr_active_e,2) == 0 .and. mod(lucita_cfg_is_spin_multiplett,2) == 0) then
          write(lupri,*) 'Illegal spin multiplicity given: # of e- even, spin mult. even.'
          call quit(' *** error in setup_lucita_check_input: spin multiplicity is impossible.')
!       case b: odd - odd
        else if(mod(lucita_cfg_nr_active_e,2) > 0 .and. mod(lucita_cfg_is_spin_multiplett,2) > 0) then
          write(lupri,*) 'Illegal spin multiplicity given: # of e- odd, spin mult. odd.'
          call quit(' *** error in setup_lucita_check_input: spin multiplicity is impossible.')
!       case c: odd - even
        else if(mod(lucita_cfg_nr_active_e,2) > 0 .and. mod(lucita_cfg_is_spin_multiplett,2) == 0) then
!         case c.1: spin mult. < 2 or > (# active e- + 1)
          if(lucita_cfg_is_spin_multiplett < 2 .or. lucita_cfg_is_spin_multiplett > (lucita_cfg_nr_active_e + 1))then
            write(lupri,'(//A/A/)') ' Illegal spin multiplicity given.',                          &
                                    ' Compare with number of active electrons.'
            call quit(' *** error in setup_lucita_check_input: spin multiplicity is impossible.')
          end if
!       case d: even - odd
        else if(mod(lucita_cfg_nr_active_e,2) == 0 .and. mod(lucita_cfg_is_spin_multiplett,2) > 0)then
!         case d.1: spin mult. < 1 or > (# active e- + 1)
          if(lucita_cfg_is_spin_multiplett < 1 .or. lucita_cfg_is_spin_multiplett > (lucita_cfg_nr_active_e+1))then
            write(lupri,'(//A/A/)') ' Illegal spin multiplicity given.',                          &
                                    ' Compare with number of active electrons.'
            call quit(' *** error in dalton_lucita: spin multiplicity is impossible.')
          end if
        end if
      end if

!     global print parameter
      if(lucita_cfg_global_print_lvl > 4)then
        write(lupri,*) 'Invalid print flag. Check .PRINTG.'
        write(lupri,*) 'lucita_cfg_global_print_lvl = ', lucita_cfg_global_print_lvl
        call quit(' *** error in setup_lucita_check_input: global print level .PRINTG too high. 0 <= range <= 4')
      end if

!     inactive shells
      if(.not.lucita_cfg_inactive_shell_set)then
        write(lupri,*) 'Number of inactive orbitals per sym'
        write(lupri,*) 'has to be specified. This is mandatory'
        write(lupri,*) 'in a GASCI/RASCI calculation.'
        call quit(' *** error in setup_lucita_check_input: inactive orbitals .INACTI not specified for RASCI calculation.')
      end if
!
!     Orbital distribution in GAS spaces (no defaults if GASCI)
      if(lucita_cfg_nr_gas_spaces < 1)then
        if(lucita_cfg_ci_type(1:6).eq.'GASCI ')then
          write(lupri,*) 'GASCI type requires .GASSHE to be specified.'
          write(lupri,*) 'Else, I do not know what to do.'
          call quit(' *** error in setup_lucita_check_input: GASCI run but .GASSHE not specified.')
        end if
      end if

!     cumulative min. and max. numbers of electrons in GAS spaces
      if(lucita_cfg_ci_type(1:6).eq.'GASCI ')then
        if(.not.lucita_cfg_minmax_occ_gas_set)then
          write(lupri,*) 'GASCI type requires .GASSPC to be specified.'
          write(lupri,*) 'Else, I do not know what to do.'
          call quit(' *** error in setup_lucita_check_input: GASCI run but .GASSPC not specified.')
        end if
      end if

!     sequence of CI spaces - default: 1 (all other values are not supported for
!     LUCITA in Dalton/Dirac)
      if(lucita_cfg_nr_calc_sequences > 1)then
        write(lupri,*) ' # of sequence of CI spaces not supported ==> ', lucita_cfg_nr_calc_sequences
        write(lupri,*) ' LUCITA in Dalton/Dirac supports only one CI space.'
        call quit(' *** error in setup_lucita_check_input: # of sequence of CI spaces not supported.')
      end if

!     parallel distribution routine
      if(luci_nmproc > 1)then
        if(lucipar_cfg_ttss_dist_strategy > 2 .or. lucipar_cfg_ttss_dist_strategy < 1)then
          write(lupri,*) 'Value for keyword .DISTRT incorrect ==> ',lucipar_cfg_ttss_dist_strategy
          write(lupri,*) 'input will be ignored'
          lucipar_cfg_ttss_dist_strategy = 2   
        end if
      end if

      print_lvl = 0
      if(ci_run_id == 'standard ci ' .or. ci_run_id == 'initial ci  ') print_lvl = 1

      if(print_lvl > 0)then
        call print_lucita_run_setup(lucita_cfg_run_title,                   &
                                    lucita_cfg_ci_type,                     &
                                    lucita_cfg_nr_roots,                    &
                                    lucita_cfg_ptg_symmetry,                &
                                    lucita_cfg_nr_active_e,                 &
                                    lucita_cfg_is_spin_multiplett,          &
                                    lucita_cfg_global_print_lvl,            &
                                    lucipar_cfg_ttss_dist_strategy,         &
                                    lucita_cfg_accepted_truncation,         &
                                    lucita_cfg_density_calc_lvl,            &
                                    luci_nmproc,                            &
                                    mxpngas,                                &
                                    mxpici,                                 &
                                    mxpirr,                                 &
                                    igsoccx,                                &
                                    ngssh,                                  &
                                    lucita_cfg_nr_ptg_irreps,               &
                                    ngas,                                   &
                                    lupri)
      end if

  end subroutine setup_lucita_check_input

!*******************************************************************************
  subroutine print_lucita_run_setup(lucita_cfg_run_title,                   &
                                    lucita_cfg_ci_type,                     &
                                    lucita_cfg_nr_roots,                    &
                                    lucita_cfg_ptg_symmetry,                &
                                    lucita_cfg_nr_active_e,                 &
                                    lucita_cfg_is_spin_multiplett,          &
                                    lucita_cfg_global_print_lvl,            &
                                    lucipar_cfg_ttss_dist_strategy,         &
                                    lucita_cfg_accepted_truncation,         &
                                    lucita_cfg_density_calc_lvl,            &
                                    luci_nmproc,                            &
                                    mxpngas,                                &
                                    mxpici,                                 &
                                    mxpirr,                                 &
                                    igsoccx,                                &
                                    ngssh,                                  &
                                    lucita_cfg_nr_ptg_irreps,               &
                                    ngas,                                   &
                                    print_unit)
!*******************************************************************************
!
!    purpose:  print LUCITA input settings.
!
!*******************************************************************************
    character (len=72), intent(in)  :: lucita_cfg_run_title
    character (len=72), intent(in)  :: lucita_cfg_ci_type
    integer,            intent(in)  :: lucita_cfg_nr_roots
    integer,            intent(in)  :: lucita_cfg_ptg_symmetry
    integer,            intent(in)  :: lucita_cfg_nr_active_e
    integer,            intent(in)  :: lucita_cfg_is_spin_multiplett
    integer,            intent(in)  :: lucita_cfg_global_print_lvl
    integer,            intent(in)  :: lucipar_cfg_ttss_dist_strategy
    integer,            intent(in)  :: lucita_cfg_density_calc_lvl
    integer,            intent(in)  :: luci_nmproc
    integer,            intent(in)  :: mxpngas
    integer,            intent(in)  :: mxpici 
    integer,            intent(in)  :: mxpirr 
    integer,            intent(in)  :: igsoccx(mxpngas,2,mxpici)
    integer,            intent(in)  :: ngssh(mxpirr,mxpngas)
    integer,            intent(in)  :: lucita_cfg_nr_ptg_irreps
    integer,            intent(in)  :: ngas
    integer,            intent(in)  :: print_unit
    real(8),            intent(in)  :: lucita_cfg_accepted_truncation
!-------------------------------------------------------------------------------
    integer                         :: i
    integer                         :: j
    integer                         :: igas
!-------------------------------------------------------------------------------

!     Title
      write(print_unit,'(/2X,60A1)') ('-',i=1,60)
      write(print_unit,'(2A)') '   title : ',lucita_cfg_run_title
      write(print_unit,'(2X,60A1)') ('-',i=1,60)
!     type of CI calculation
      write(print_unit,'(a,A6/)') '   Type of calculation .................. ', lucita_cfg_ci_type(1:6)
!     number of roots to be treated
      write(print_unit,'(a,I3)')  '   Number of roots to be obtained ....... ', lucita_cfg_nr_roots
!     state symmetry
      write(print_unit,'(a,I3)')  '   Calculation carried out in irrep ..... ', lucita_cfg_ptg_symmetry
!     number of active electrons
      write(print_unit,'(a,I3)')  '   Number of active electrons ........... ', lucita_cfg_nr_active_e
!     spin multiplicity
      write(print_unit,'(a,I3)')  '   Spin multiplicity .................... ', lucita_cfg_is_spin_multiplett
!     LUCITA global print parameter
      write(print_unit,'(a,i3)')  '   Global print level is ................ ', lucita_cfg_global_print_lvl
!     parallel distribution routine
      if(LUCI_NMPROC > 1)then
        write(print_unit,'(a,I3)')'   Parallel distribution routine ........ ', lucipar_cfg_ttss_dist_strategy
      end if
!     truncation factor
      if(lucita_cfg_accepted_truncation /= 1.0d-10)then
        write(print_unit,'(/a,1P,D10.2)') '   Truncation Factor .................... ',lucita_cfg_accepted_truncation
      end if
!     density matrix level
      if(lucita_cfg_density_calc_lvl >= 1)then
        write(print_unit,'(a,I3)') '   Density matrices level ............... ',lucita_cfg_density_calc_lvl
      end if

!     GAS spaces and occupation constraints
      write(print_unit,'(/A/A/A/)')   &
         ' ******************************************************************',   &
         ' Generalized active space: shell spaces and occupation constraints ',   &
         ' ******************************************************************'

      select case(lucita_cfg_nr_ptg_irreps)
      case(1)
        write(print_unit,'(a,1i4,a)')       &
          '                 Irrep:',(i,i = 1,lucita_cfg_nr_ptg_irreps),    '        Min. occ    Max. occ '
        write(print_unit,'(10a)')           &
          '                 ========',('====',i = 1,lucita_cfg_nr_ptg_irreps), '      ========    ======== '
        do igas = 1, ngas
          write(print_unit,'(a,i2,a,1i4,i13,i12)') &
          '        GAS',igas,'          ',            &
          (ngssh(i,igas),i = 1, lucita_cfg_nr_ptg_irreps),igsoccx(igas,1,1),igsoccx(igas,2,1)
        end do
      case(2)
        write(print_unit,'(/a,2i4,a)')      &
          '                 Irrep:',(i,i = 1,lucita_cfg_nr_ptg_irreps),    '        Min. occ    Max. occ '
        write(print_unit,'(10a)')           &
          '                 ========',('====',i = 1,lucita_cfg_nr_ptg_irreps), '      ========    ======== '
        do igas = 1, ngas
          write(print_unit,'(a,i2,a,2i4,i13,i12)') &
          '        GAS',igas,'          ',            &
          (ngssh(i,igas),i = 1, lucita_cfg_nr_ptg_irreps),igsoccx(igas,1,1),igsoccx(igas,2,1)
        end do
      case(4)
        write(print_unit,'(/a,4i4,a)')      &
          '                 Irrep:',(i,i = 1,lucita_cfg_nr_ptg_irreps),    '        Min. occ    Max. occ '
        write(print_unit,'(10a)')           &
          '                 ========',('====',i = 1,lucita_cfg_nr_ptg_irreps), '      ========    ======== '
        do igas = 1, ngas
          write(print_unit,'(a,i2,a,4i4,i13,i12)') &
          '        GAS',igas,'          ',            &
          (ngssh(i,igas),i = 1, lucita_cfg_nr_ptg_irreps),igsoccx(igas,1,1),igsoccx(igas,2,1)
        end do
      case(8)
        write(print_unit,'(/a,8i4,a)')      &
          '                 Irrep:',(i,i = 1,lucita_cfg_nr_ptg_irreps),    '        Min. occ    Max. occ '
        write(print_unit,'(10a)')           &
          '                 ========',('====',i = 1,lucita_cfg_nr_ptg_irreps), '      ========    ======== '
        do igas = 1, ngas
          write(print_unit,'(a,i2,a,8i4,i13,i12)') &
          '        GAS',igas,'          ',            &
          (ngssh(i,igas),i = 1, lucita_cfg_nr_ptg_irreps),igsoccx(igas,1,1),igsoccx(igas,2,1)
        end do
      end select


  end subroutine print_lucita_run_setup
!*******************************************************************************

  subroutine setup_lucita_par_dist_in_seq(par_dist_block_list,   &
                                          block_list,            &
                                          nblock)
!*******************************************************************************
!
!    purpose:  set parallel distribution list to master-only list.
!              this can be useful in sequential calculations. 
!
!*******************************************************************************
    integer, intent(in)           :: nblock
    integer, intent(in)           :: block_list(nblock)
    integer, intent(out)          :: par_dist_block_list(nblock)
!-------------------------------------------------------------------------------
    integer                       :: i
!-------------------------------------------------------------------------------

      i = 1

      do while (i <= nblock)
        if(block_list(i) > 0) par_dist_block_list(i) = 0 
        i = i + 1
      end do

  end subroutine setup_lucita_par_dist_in_seq
!*******************************************************************************

  subroutine setup_lucita_mcci_wrkspc_dimensions(ci_run_id,print_lvl)
!*******************************************************************************
!
!    purpose:  set dimensions of incoming matrices to LUCITA 
!              from the MCSCF environment.  
!              note that the allocation is CI-run dependent.
!
!*******************************************************************************
  use lucita_mcscf_ci_cfg
  use lucita_cfg
#include "priunit.h"
#include "mxpdim.inc"
#include "crun.inc"
#include "orbinp.inc"
#include "parluci.h"
    character(len=12), intent(in) :: ci_run_id
    integer, intent(in)           :: print_lvl
!-------------------------------------------------------------------------------
    integer                       :: tmp_length_p1
    integer                       :: tmp_length_p2
!-------------------------------------------------------------------------------

!     initialize internal length of matrices
      len_cref_mc2lu             = 0
      len_hc_mc2lu               = 0
      len_resolution_mat_mc2lu   = 0
      len_int1_or_rho1_mc2lu     = 0
      len_int2_or_rho2_mc2lu     = 0
      tmp_length_p1              = 0
      tmp_length_p2              = 0

      if(integrals_from_mcscf_env)then
!       1-p density matrix / 1-e integrals
        tmp_length_p1 = nacob**2
!       2-p density matrix / 2-e integrals
        if(lucita_cfg_el_operator_level > 1) tmp_length_p2 = (nacob*(nacob+1)/2)**2
      end if

      select case(ci_run_id)
        case('return CIdim', 'ijkl resort ', 'fci dump    ') ! calculate # of determinants per symmetry irrep / resort integrals to lucita format
!         nothing needs to be allocated
        case('xc vector   ') ! exchange CI/MCSCF vector
          len_cref_mc2lu              =  l_combi
        case('return CIdia', 'analyze Cvec') ! calculate the diagonal part of the CI Hamiltonian matrix , analyze CI vector
          len_cref_mc2lu              =  lblock
        case('sigma vec   ', 'Xp-density m', 'rotate  Cvec', 'standard ci ', 'initial ci  ', 'srdft   ci  ')
!      calculate sigma vector; 1-/2-particle density matrix; rotate CI vector; perform Davidson CI run
          len_cref_mc2lu              =  lblock
          len_hc_mc2lu                =  lblock
          len_resolution_mat_mc2lu    =  max(lblock,2*lsingle_resolution_block)
          len_int1_or_rho1_mc2lu      =  tmp_length_p1
          len_int2_or_rho2_mc2lu      =  tmp_length_p2
        case default
          print *, ' unknown CI run id: ',ci_run_id,' no memory allocated for the coworkers.'
      end select

      if(print_lvl > +1)then
        write(lupri,'(/a     )') '  dimensions of CI matrix / integral arrays:'
        write(lupri,'( a     )') '  -----------------------------------------'
        write(lupri,'( a,i15 )') '  ci matrix #1     ==> ',len_cref_mc2lu
        write(lupri,'( a,i15 )') '  ci matrix #2     ==> ',len_hc_mc2lu
        write(lupri,'( a,i15 )') '  resolution block ==> ',len_resolution_mat_mc2lu
        write(lupri,'( a,i15 )') '  1-el ij          ==> ',len_int1_or_rho1_mc2lu
        write(lupri,'( a,i15/)') '  2-el ijkl        ==> ',len_int2_or_rho2_mc2lu
        call flshfo(lupri)
      end if

  end subroutine setup_lucita_mcci_wrkspc_dimensions
!*******************************************************************************

#ifdef VAR_MPI

  subroutine setup_lucita_inc_wrkspc_alloc_cw(lucita_ci_run_id,             &
                                              work_dalton,                  &
                                              kfree_pointer,                &
                                              lfree,                        &
                                              kcref_pointer,                &
                                              khc_pointer,                  &
                                              kresolution_mat_pointer,      &
                                              kint1_or_rho1_pointer,        &
                                              kint2_or_rho2_pointer,        &
                                              len_cref,                     &
                                              len_hc,                       &
                                              len_resolution_mat,           &
                                              len_int1_or_rho1,             &
                                              len_int2_or_rho2,             &
                                              print_lvl)
!*******************************************************************************
!
!    purpose:  allocate matrices for co-workers which are incoming to LUCITA 
!              from the Dalton sirius environment.  
!              note that the allocation is CI-run dependent.
!
!*******************************************************************************
    character(len=12), intent(in) :: lucita_ci_run_id
    real(8), intent(inout)        :: work_dalton(*)
    integer, intent(inout)        :: kfree_pointer
    integer, intent(inout)        :: kcref_pointer
    integer, intent(inout)        :: khc_pointer
    integer, intent(inout)        :: kresolution_mat_pointer
    integer, intent(inout)        :: kint1_or_rho1_pointer
    integer, intent(inout)        :: kint2_or_rho2_pointer
    integer, intent(inout)        :: lfree
    integer, intent(in)           :: len_cref
    integer, intent(in)           :: len_hc
    integer, intent(in)           :: len_resolution_mat
    integer, intent(in)           :: len_int1_or_rho1
    integer, intent(in)           :: len_int2_or_rho2
    integer, intent(in)           :: print_lvl
!-------------------------------------------------------------------------------
    integer                       :: len_cref_internal
    integer                       :: len_hc_internal
    integer                       :: len_resolution_mat_internal
    integer                       :: len_int1_or_rho1_internal
    integer                       :: len_int2_or_rho2_internal
!-------------------------------------------------------------------------------

!     initialize internal length of matrices
      len_cref_internal           =  0
      len_hc_internal             =  0
      len_resolution_mat_internal =  0
      len_int1_or_rho1_internal   =  0
      len_int2_or_rho2_internal   =  0

!     initialize pointers
      kcref_pointer               = -1
      khc_pointer                 = -1
      kresolution_mat_pointer     = -1
      kint1_or_rho1_pointer       = -1
      kint1_or_rho1_pointer       = -1

      select case(lucita_ci_run_id)
        case('return CIdim', 'ijkl resort ', 'fci dump    ') ! calculate # of determinants per symmetry irrep / resort integrals to lucita format
!         nothing needs to be allocated
        case('return CIdia') ! calculate the diagonal part of the CI Hamiltonian matrix
          len_cref_internal           =  len_cref
          len_int1_or_rho1_internal   =  len_int1_or_rho1
          len_int2_or_rho2_internal   =  len_int2_or_rho2
        case('sigma vec   ',      'Xp-density m',               'rotate  Cvec',       'standard ci ', 'initial ci  ') 
!      calculate sigma vector; 1-/2-particle density matrix;   rotate CI vector;          perform Davidson CI run
          len_cref_internal           =  len_cref
          len_hc_internal             =  len_hc
          len_resolution_mat_internal =  len_resolution_mat
          len_int1_or_rho1_internal   =  len_int1_or_rho1
          len_int2_or_rho2_internal   =  len_int2_or_rho2
!            analyze CI vector; exchange CI/MCSCF vector
        case('analyze Cvec',    'xc vector   ')
          len_cref_internal           =  len_cref
        case default
          print *, ' unknown CI run id: ',lucita_ci_run_id,' no memory allocated for the coworkers.'
      end select

      call memget('REAL',kcref_pointer,          len_cref_internal,          work_dalton,kfree_pointer,lfree)
      call memget('REAL',khc_pointer,            len_hc_internal,            work_dalton,kfree_pointer,lfree)
      call memget('REAL',kresolution_mat_pointer,len_resolution_mat_internal,work_dalton,kfree_pointer,lfree)
      call memget('REAL',kint1_or_rho1_pointer,  len_int1_or_rho1_internal,  work_dalton,kfree_pointer,lfree)
      call memget('REAL',kint2_or_rho2_pointer,  len_int2_or_rho2_internal,  work_dalton,kfree_pointer,lfree)

      if(print_lvl > 2)then
        print *, ' memory allocated for the coworkers; pointer values are:'
        print *, ' -------------------------------------------------------'
        print *, ' kcref_pointer           ==> ',kcref_pointer
        print *, ' khc_pointer             ==> ',khc_pointer
        print *, ' kresolution_mat_pointer ==> ',kresolution_mat_pointer
        print *, ' kint1_or_rho1_pointer   ==> ',kint1_or_rho1_pointer
        print *, ' kint2_or_rho2_pointer   ==> ',kint2_or_rho2_pointer
        print *, ' -------------------------------------------------------'

        print *, ' memory allocated for the coworkers; length  values are:'
        print *, ' -------------------------------------------------------'
        print *, ' kcref_pointer           ==> ',len_cref_internal
        print *, ' khc_pointer             ==> ',len_hc_internal
        print *, ' kresolution_mat_pointer ==> ',len_resolution_mat_internal
        print *, ' kint1_or_rho1_pointer   ==> ',len_int1_or_rho1_internal
        print *, ' kint2_or_rho2_pointer   ==> ',len_int2_or_rho2_internal
      end if

  end subroutine setup_lucita_inc_wrkspc_alloc_cw
!*******************************************************************************
#endif /* VAR_MPI */

end module
