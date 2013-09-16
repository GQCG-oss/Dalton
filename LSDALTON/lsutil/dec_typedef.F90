!> @file
!> The module contains structures used in DEC
!> \author Kasper Kristensen

!> The module contains structures used in DEC (and printing routines)
module dec_typedef_module

  use precision
  use,intrinsic :: iso_c_binding, only:c_ptr
  use TYPEDEFTYPE, only: lsitem
  use Matrix_module, only: matrix
  !Could someone please rename ri to something less generic. TK!!
  private
  public :: DECinfo, ndecenergies,DECsettings,array2,array3,array4,ccorbital,ri,&
       & fullmolecule,ccatom,FullMP2grad,mp2dens,mp2grad,&
       & mp2_batch_construction,mypointer,joblist,traceback,batchTOorb,&
       & SPgridbox
  ! IMPORTANT: Number of possible energies to calculate using the DEC scheme
  ! MUST BE UPDATED EVERYTIME SOMEONE ADDS A NEW MODEL TO THE DEC SCHEME!!!!
  ! MODIFY FOR NEW MODEL
  ! MODIFY FOR NEW CORRECTION
  integer, parameter :: ndecenergies = 14

  !> \author Kasper Kristensen
  !> \date June 2010
  !> \brief Contains settings for DEC calculation, see default settings in dec_set_default_config.
  type DECsettings

     ! ****************************************************************************************
     !                 !!!!!!!!!!!! VERY VERY IMPORTANT !!!!!!!!
     !
     !  IF YOU REMOVE/ADD MEMBERS OF/TO THIS STRUCTURE, REMEMBER TO MODIFY mpicopy_dec_settings
     !  IN decmpi.f90 ACCORDINGLY!!!!!!
     ! 
     ! *****************************************************************************************



     !> MAIN SETTINGS DEFINING DEC CALCULATION
     !> **************************************

     !> Run DEC calculation at all?
     logical :: doDEC
     !> Frozen core calculation?
     logical :: frozencore
     !> Full molecular job
     logical :: full_molecular_cc ! full molecular cc
     !> Enforce canonical orbitals in calculation 
     !> (only meaningful for full_molecular_cc or simulate_full)
     logical :: use_canonical
     !> Default: Full calculation uses canonical orbitals, while DEC calculation uses local orbitals.
     !> This can be overruled by inpiut keywords (see config_dec_input).
     !> If the default choice was overruled "user_defined_orbitals" is set to true.
     logical :: user_defined_orbitals
     !> Simulate full molecular calculation in DEC mode  (debug)
     logical :: simulate_full
     !> How many atoms to use in simulation mode   (debug)
     integer :: simulate_natoms
     !> Includes the whole molecule in all fragments (debug)
     logical :: InclFullMolecule
     !> Label to print CC models
     character(len=8), dimension(10) :: cc_models
     !> Requested CC model
     integer :: ccModel ! 1 - MP2, 2 - CC2, 3 - CCSD, 4 - CCSD(T), 5 - RPA
     !> Use singles
     logical :: use_singles


     !> Restart options
     !> ***************
     !> Restart calculation if some of the fragments were calculated in a previous calculation
     logical :: restart
     !> Creating files for restart: Time (in seconds) passing before backing up restart files
     real(realk) :: TimeBackup
     !> Read DEC orbital file DECOrbitals.info from file (default: Existing file is overwritten)
     logical :: read_dec_orbitals


     !> Memory stuff
     !> ************
     !> Memory available for DEC calculation
     real(realk) :: memory
     !> Memory defined by input? (If not, use system call).
     logical :: memory_defined
     ! Memory use for full molecule structure
     real(realk) :: fullmolecule_memory
     !> Save array4 on file instead of in memory (debug)
     logical :: array4OnFile
     logical :: array4OnFile_specified


     !> Singles polarization (currenly turned off)
     !> ******************************************
     !> Construct full molecular singles amplitudes from fragment calculations?
     logical :: SinglesPolari
     !> Relative difference between singles amplitudes to accept
     !> without invoking an additional set of fragment calculations
     real(realk) :: singlesthr


     !> 32 vs. 64 bit issues
     !> ********************
     !> DEC files (lcm_orbitals.u, fock.restart, dens.restart, overlapmatrix, DECorbitals.info) 
     !> are in 64 (or 32) bit integers but the program was compiled with 32 (or 64) 
     !> bit integers so these files need to be converted during the read-in.
     logical :: convert64to32
     logical :: convert32to64


     !> CCSD residual/solver settings
     !> *****************************
     !> save next guess amplitudes of CCSD in each iteration on disk
     logical :: CCSDnosaferun
     !> Use parallel CCSD solver
     logical :: solver_par
     !> forcing one or the other scheme in get_coubles residual integral_driven
     logical :: force_scheme
     logical :: dyn_load
     !> Use the cc debug routines from the file cc_debug_routines.F90
     logical :: CCDEBUG
     !> skip reading the old amplitudes from disk
     logical :: CCSDno_restart
     !> if mpich is used CCSD has some special treats that can be used
     logical :: CCSD_MPICH
     !> prevent canonicalization in the ccsolver
     logical :: CCSDpreventcanonical
     !> chose left-transformations to be carried out
     logical :: CCSDmultipliers
     !> do not update the singles residual
     logical :: CCDhack
     !> Debug CC driver
     logical :: cc_driver_debug
     !> Integer specifying which scheme to use in CCSD calculations (debug)
     integer :: en_mem
     !> Use full molecular Fock matrix to precondition
     logical :: precondition_with_full
     !> Use devel version of doubles
     logical :: ccsd_expl
     !> Max number of iterations in CC driver
     integer :: ccMaxIter
     !> Max number of vectors in the subspace
     integer :: ccMaxDIIS
     !> CC convergence threshold
     real(realk) :: ccConvergenceThreshold
     !> Was CC convergence threshold specified?
     logical :: CCthrSpecified
     !> Use preconditioner
     logical :: use_preconditioner
     !> Use preconditioner in B matrix
     logical :: use_preconditioner_in_b
     !> Use CROP (if false we use DIIS)
     logical :: use_crop
     !> Simulate full ERI using RI arrays 
     !> (obsolete for the moment, Patrick will remove when cleaning the CC solver)
     logical :: simulate_eri
     !> Construct Fock matrix from RI integrals (obsolete for the moment)
     !> (obsolete for the moment, Patrick will remove when cleaning the CC solver)
     logical :: fock_with_ri

     !> F12 settings
     !> ************
     !> Use F12 correction
     logical :: F12

     !> F12 debug settings
     !> ************
     !> Use F12 correction
     logical :: F12DEBUG

     !> MPI settings
     !> ************
     !> Factor determining when MPI groups should split
     integer :: mpisplit
     !> Manually set starting group size for local MPI group
     integer(kind=ls_mpik) :: MPIgroupsize

     !> Integral batching
     !> *****************
     !> Set integral batch sizes manually
     logical :: manual_batchsizes
     !> Sizes of alpha and gamma batches defined manually
     integer :: ccsdAbatch,ccsdGbatch


     !> General debug and simple tests
     !> ******************************
     !> General HACK parameters, to be used for easy debugging
     logical :: hack
     logical :: hack2
     !> Full calculation where individual pair and single energies are calculated in ONE energy calc.
     !> Only implemented for MP2 and only for debugging purposes.
     logical :: mp2energydebug
     !> Skip the read-in of molecular info files dens.restart, fock.restart, lcm_orbitals.u
     logical :: SkipReadIn
     !> test the array structure
     logical :: array_test
     !> test the array reorderings
     logical :: reorder_test
     !> Check that LCM orbitals are correct
     logical :: check_lcm_orbitals
     !> Debug print level
     integer :: PL
     !> only do fragment part of density or gradient calculation 
     logical :: SkipFull 
     ! --

     !> Output options 
     !> **************
     !> File unit for LSDALTON.OUT
     integer :: output
     ! --

     !> DEC Orbital treatment
     !> *********************
     !> Absorb H atoms into heavy atoms during orbital assignment
     logical :: AbsorbHatoms
     !> Fit orbital coefficients in fragment (default: true)
     logical :: FitOrbitals
     !> Threshold for simple Lowdin procedure for determining atomic extent
     real(realk) :: simple_orbital_threshold
     !> Purify fitted MO coefficients (projection + orthogonalization)
     logical :: PurifyMOs
     !> Use fragment-adapted orbitals for fragment calculations
     logical :: FragAdapt
     !> Has simple orbital threshold been defined manually in input (true),
     !> or should simple orbital threshold be adapted to FOT 
     !> as descripted under FOTlevel (false)?
     logical :: simple_orbital_threshold_set     
     !> Use Boughton-Pulay criteria for generating orbitals rather than simple Lowdin charge criteria
     logical :: BoughtonPulay
     !> Simple Mulliken charge threshold (only for Boughton-Pulay procedure)
     real(realk) :: mulliken_threshold
     !> Simple Mulliken charge criteria 
     logical :: simple_mulliken_threshold
     !> Norm error in approximated (fitted orbitals)
     real(realk) :: approximated_norm_threshold
     !> Use Mulliken population analysis to assign orbitals (default: Lowdin, only for Boughton-Pulay)
     logical :: mulliken
     !> Use Distance criteria to determine central atom
     logical :: Distance
     ! --


     !> Fragment optimization
     !> *********************
     !> Fragment optimization threshold
     real(realk) :: FOT
     !> Max number of iterations for expanding fragment
     integer :: MaxIter
     !> FOT level defining precision of calculation, see set_input_for_fot_level
     integer :: FOTlevel
     !> Max accepted FOT level
     integer :: maxFOTlevel
     !> Use occupied/virtual hybrid partitioning scheme
     !> This is a temporary solution for models where the Lagrangian has not yet been implemented.
     !> The Lagrangian fragment energy is then simply defined as the average value of the two energies
     !> for the occupied and virtual partitioning schemes.
     logical :: HybridScheme
     !> Number of atoms to include in fragment expansion
     integer :: FragmentExpansionSize
     !> Use MP2 energies for expansion part of fragment optimization
     logical :: fragopt_exp_mp2
     !> Use MP2 energies for reduction part of fragment optimization
     logical :: fragopt_red_mp2
     !> Only consider occupied partitioning
     logical :: OnlyOccPart
     !> Repeat atomic fragment calculations after fragment optimization?
     ! (this is necessary e.g. for gradient calculations).
     logical :: RepeatAF
     ! --  

     !> Pair fragments
     !> **************
     !> Distance cutoff for pair fragments
     real(realk) :: pair_distance_threshold
     !> Pair cutoff set manually (will overwrite default pair cutoff defined by FOTlevel)
     logical :: paircut_set
     !> Pair distance beyond which reduced fragments are used
     real(realk) :: PairReductionDistance
     !> When pair regression fit is performed, pair distances smaller than PairMinDist are ignored
     real(realk) :: PairMinDist
     !> Skip pair analysis (debug)
     logical :: checkpairs
     ! --


     ! First order properties
     ! **********************
     !> Do first order properties (MP2 density, electric dipole, mp2 gradient)
     logical :: first_order
     !> MP2 density matrix (and not gradient)
     logical :: MP2density    
     !> Calculate MP2 gradient  (density is then also calculated as a subset of the calculation)
     logical :: gradient
     !> Use preconditioner for kappa multiplier equation
     logical :: kappa_use_preconditioner
     !> Use preconditioner for kappa multiplier equation
     logical :: kappa_use_preconditioner_in_b
     !> Number of vectors to save in kappa multiplier equation
     integer :: kappaMaxDIIS
     !> Maximum number of iteration in kappa multiplier equation
     integer :: kappaMaxIter
     !> Debug print in kappa multiplier equation
     logical :: kappa_driver_debug
     !> Residual threshold for kappa orbital rotation multiplier equation
     real(realk) :: kappaTHR


     !> Geometry optimization
     !> *********************
     !> Book keeping of the number of DEC calculations for each FOT level
     !> (only relevant for geometry optimizations)
     integer,dimension(8) :: ncalc
     !> Factor multiply intrinsic energy error by before returning error to geometry optimizer
     real(realk) :: EerrFactor
     !> Old energy error (used only for geometry opt)
     real(realk) :: EerrOLD


  end type DECSETTINGS



  type array2

     integer, dimension(2) :: dims
     real(realk), pointer :: val(:,:) => null()

  end type array2

  type array3

     !> Dimensions
     integer, dimension(3) :: dims
     !> Current order
     integer, dimension(3) :: order
     !> Data
     real(realk), pointer :: val(:,:,:) => null()

  end type array3

  type array4

     !> Dimensions
     integer, dimension(4) :: dims
     !> Current order
     integer, dimension(4) :: order
     !> Data
     real(realk), pointer :: val(:,:,:,:) => null()
     !> File unit counter
     integer :: FUnit
     !> File name
     character(len=80) :: FileName

     ! Information only used when array4s are stored on file (array4OnFile= .true.)
     ! ****************************************************************************

     !> Address counter
     integer(kind=long) :: address_counter

     !> Storing type
     integer :: storing_type
     ! Storing_type = 2: The values are stored on file as: val(:,:,n1,n2)
     ! Storing_type = 3: The values are stored on file as: val(:,n1,n2,n3)


     !> Number of elements stored in each address on file
     integer(kind=long) :: nelements
     ! For storing_type=2, nelements=dims(1)*dims(2)
     ! For storing_type=3, nelements=dims(1)

     !> List of addresses for the values stored on file
     ! E.g. If one uses storing_type=2 and wants to read in the values
     ! val(:,:,n1,n2), then the corresponding address on file is given by
     ! address(1,1,n1,n2)
     ! where it is understood that all elements are read in for the first two elements
     ! because storing_type=2.
     ! (address is a 4-dimensional array to keep it general and extendable to other cases)
     integer(kind=long), pointer :: address(:,:,:,:) => null()


  end type array4


  type ccorbital

     !> Number of the orbital in full molecular basis
     integer :: orbitalnumber
     !> Cental atom in the population
     integer :: centralatom
     !> Number of significant atoms
     integer :: numberofatoms

     !> List of significant atoms
     integer, pointer :: atoms(:) => null()

  end type ccorbital


  !> Three dimensional array
  type ri

     !> Dimensions
     integer, dimension(3) :: dims
     !> Data
     real(realk), pointer :: val(:,:,:) => null()
     !> File unit
     integer :: FUnit
     !> File name
     character(len=80) :: FileName

  end type ri



  !> All information about full molecule and HF calculation
  type fullmolecule

     !> Number of electrons
     integer :: nelectrons
     !> Number of atoms
     integer :: natoms
     !> Number of basis functions
     integer :: nbasis
     !> Number of auxiliary basis functions
     integer :: nauxbasis
     !> Number of occupied orbitals (core + valence)
     integer :: numocc
     !> Number of core orbitals
     integer :: ncore
     !> Number of valence orbitals (numocc-ncore)
     integer :: nval
     !> Number of unoccupied orbitals
     integer :: numvirt

     !> Number of basis functions on atoms
     integer, pointer :: atom_size(:) => null()
     !> Index of the first basis function for an atom
     integer, pointer :: atom_start(:) => null()
     !> Index of the last basis function for an atom
     integer, pointer :: atom_end(:) => null()

     !> Occupied MO coefficients (mu,i)
     real(realk), pointer :: ypo(:,:) => null()
     !> Virtual MO coefficients (mu,a)
     real(realk), pointer :: ypv(:,:) => null()
     !> CABS MO coefficients (mu,x)
     real(realk), pointer :: cabsMOs(:,:) => null()

     !> Fock matrix (AO basis)
     real(realk), pointer :: fock(:,:) => null()
     !> Overlap matrix (AO basis)
     real(realk), pointer :: overlap(:,:) => null()

     !> Occ-occ block of Fock matrix in MO basis
     real(realk), pointer :: ppfock(:,:) => null()
     !> Virt-virt block of Fock matrix in MO basis
     real(realk), pointer :: qqfock(:,:) => null()
     !> carmom coord for occ
     real(realk), pointer :: carmomocc(:,:) => null()
     !> carmom coord for virt
     real(realk), pointer :: carmomvirt(:,:) => null()
     !> atomic centers
     real(realk), pointer :: AtomCenters(:,:) => null()

  end type fullmolecule


  !> Atomic fragment / Atomic pair fragment
  type ccatom

     !> Number of atom in full molecule
     integer :: atomic_number=0
     !> Number of occupied EOS orbitals 
     integer :: noccEOS=0
     !> Number of unoccupied EOS orbitals 
     integer :: nunoccEOS=0
     !> Number of occupied AOS orbitals (for frozen core approx this is only the valence orbitals)
     integer :: noccAOS=0
     !> Number of core orbitals in AOS
     integer :: ncore=0
     !> Total number of orbitals (core+valence) in AOS (noccAOS + ncore)
     integer :: nocctot=0
     !> Total number of unoccupied orbitals (AOS)
     integer :: nunoccAOS=0

     !> Pair fragment?
     logical :: pairfrag

     !> Occupied orbital EOS indices 
     integer, pointer :: occEOSidx(:) => null()
     !> Unoccupied orbital EOS indices 
     integer, pointer :: unoccEOSidx(:) => null()
     !> Occupied AOS orbital indices (only valence orbitals for frozen core approx)
     integer, pointer :: occAOSidx(:) => null()
     !> Unoccupied AOS orbital indices 
     integer, pointer :: unoccAOSidx(:) => null()
     !> Core orbitals indices (only used for frozen core approx, 
     !> otherwise there are included in the occAOSidx list).
     integer,pointer :: coreidx(:) => null()


     ! Special info for reduced fragment of lower accuracy
     !****************************************************
     !> Number of occupied AOS orbitals 
     integer :: REDnoccAOS=0
     !> Total number of unoccupied orbitals (AOS)
     integer :: REDnunoccAOS=0
     !> Occupied orbital indices (AOS) for reduced fragment
     integer, pointer :: REDoccAOSidx(:) => null()
     !> All unoccupied orbital indices (AOS)
     integer, pointer :: REDunoccAOSidx(:) => null()

     !> Indices of occupied EOS in AOS basis
     integer, pointer :: idxo(:) => null()
     !> Indices of unoccupied EOS in AOS basis
     integer, pointer :: idxu(:) => null()

     ! MODIFY FOR NEW MODEL
     ! MODIFY FOR NEW CORRECTION
     !> DEC fragment energies stored in the following manner:
     !> 1. MP2 Lagrangian partitioning scheme
     !> 2. MP2 occupied partitioning scheme
     !> 3. MP2 virtual partitioning scheme
     !> 4. CC2 occupied partitioning scheme
     !> 5. CC2 virtual partitioning scheme
     !> 6. CCSD occupied partitioning scheme
     !> 7. CCSD virtual partitioning scheme
     !> 8. (T) contribution, occupied partitioning scheme
     !> 9. (T) contribution, virtual partitioning scheme
     !> 10. Fourth order (T) contribution, occupied partitioning scheme
     !> 11. Fourth order (T) contribution, virtual partitioning scheme
     !> 12. Fifth order (T) contribution, occupied partitioning scheme
     !> 13. Fifth order (T) contribution, virtual partitioning scheme
     !> 14. MP2-F12 energy correction

     real(realk),dimension(ndecenergies) :: energies
     ! Note 1: Only the energies requested for the model in question are calculated!
     ! Note 2: Obviously you need to change the the global integer "ndecenergies"
     !         at the top of this file if you add new models!!!


     !> The energy definitions below are only used for fragment optimization (FOP)
     !> These are (in general) identical to the corresponding energies saved in "energies".
     !> However, for fragment optimization it is very convenient to have direct access to the energies
     !> without thinking about which CC model we are using...
     !> Energy using occupied partitioning scheme
     real(realk) :: EoccFOP
     !> Energy using virtual partitioning scheme
     real(realk) :: EvirtFOP
     !> Lagrangian energy 
     !> ( = 0.5*OccEnergy + 0.5*VirtEnergy for models where Lagrangian has not been implemented)
     real(realk) :: LagFOP
  
     !> Contributions to the fragment Lagrangian energy from each individual
     !  occupied or virtual orbital.
     real(realk),pointer :: OccContribs(:) => null()
     real(realk),pointer :: VirtContribs(:) => null()

     !> Number of EOS atoms (1 for atomic fragment, 2 for pair fragment)
     integer :: nEOSatoms
     !> List of EOS atoms
     integer, pointer :: EOSatoms(:) => null()


     !> Information used only when the ccatom is a pair fragment
     !> ********************************************************
     !> Distance between single fragments used to generate pair
     real(realk) :: pairdist

     !> Total occupied orbital space (orbital type)
     type(ccorbital), pointer :: occAOSorb(:) => null()
     !> Total unoccupied orbital space (orbital type)
     type(ccorbital), pointer :: unoccAOSorb(:) => null()

     !> Number of atoms (atomic extent)
     integer :: number_atoms=0
     !> Number of basis functions
     integer :: number_basis=0
     !> Atomic indices
     integer, pointer :: atoms_idx(:) => null()
     !> Corresponding basis function indices
     integer,pointer :: basis_idx(:) => null()

     !> Has the information inside the expensive box below been initialized or not?
     logical :: BasisInfoIsSet

     ! ===========================================================================
     !                       IMPORTANT: EXPENSIVE BOX
     ! ===========================================================================
     ! The information inside this "expensive box" is what takes the time when a
     ! fragment is initialized (MO coefficients, integral input etc.)
     ! When the fragment is initialized using atomic_fragment_init_orbital_specific
     ! with DoBasis=.false. then this information is NOT SET!
     ! In this way the basic fragment information (everything ouside the expensive box)
     ! can be obtained in a very cheap manner, which is convenient for
     ! planning of a large number of fragment calculations.
     ! ---------------------------------------------------------------------------

     !> AO overlap matrix for fragment
     real(realk),pointer :: S(:,:) => null()

     !> Occupied MO coefficients (only valence space for frozen core approx)
     real(realk), pointer :: ypo(:,:) => null()
     !> Virtual MO coefficients
     real(realk), pointer :: ypv(:,:) => null()
     !> Cabs MO coefficients
     real(realk),pointer :: cabsMOs(:,:) => null()     
     !> Core MO coefficients 
     real(realk),pointer :: CoreMO(:,:) => null()

     !> AO Fock matrix
     real(realk), pointer :: fock(:,:) => null()
     !> Occ-occ block of Fock matrix in MO basis  (only valence space for frozen core approx)
     real(realk), pointer :: ppfock(:,:) => null()
     !> Virt-virt block of Fock matrix in MO basis
     real(realk), pointer :: qqfock(:,:) => null()
     !> Core-core block of Fock matrix in MO basis  (subset of ppfock when frozen core is NOT used)
     real(realk), pointer :: ccfock(:,:) => null()

     !> Integral program input
     type(lsitem) :: mylsitem

     ! End of EXPENSIVE BOX
     ! ==============================================================

     
     ! Information used for fragment-adapted orbitals
     ! *******************************************
     !> Correlation density matrices in local AOS basis
     real(realk), pointer :: OccMat(:,:) => null()  ! occ AOS-EOS
     real(realk), pointer :: VirtMat(:,:) => null()  ! virt AOS-EOS
     !> Threshold to use for throwing away fragment-adapted occupied (1) or virtual (2) orbitals
     real(realk) :: RejectThr(2)
     !> Control of whether corr dens matrices have been set (true) or simply initialized (false)
     logical :: CDset
     !> Is this a fragment-adapted fragment?
     logical :: fragmentadapted
     !> Number of occ orbitals for fragment-adapted orbitals 
     integer :: noccFA
     !> Number of unocc orbitals for fragment-adapted orbitals 
     integer :: nunoccFA
     !> Transformation between AO basis and fragment-adapted basis
     !> Index 1: Local,   Index 2: Fragment-adapted
     !> Has fragment-adapted MO coeff been set (not done by default fragment initialization)?
     logical :: FAset
     real(realk),pointer :: CoccFA(:,:) => null()     ! dimension: number_basis,noccFA
     real(realk),pointer :: CunoccFA(:,:) => null()   ! dimension: number_basis,nunoccFA


     !> Information used only for the CC2 and CCSD models to describe
     !> long-range effects described by singles amplitudes properly.
     !> *************************************************************
     !> Are t1 amplitudes stored in the fragment structure?
     logical :: t1_stored
     !> Dimensions of t1 amplitudes (virtual,occupied - can be either EOS or AOS)
     integer,dimension(2) :: t1dims
     !> t1 amplitudes
     real(realk),pointer :: t1(:,:) => null()
     !> Indices for occupied fragment indices in full list of orbitals
     integer,pointer :: t1_occidx(:) => null()
     !> Indices for virtual fragment indices in full list of orbitals
     integer,pointer :: t1_virtidx(:) => null()


     ! FLOP ACCOUNTING
     ! ***************
     ! MPI: Sum of flop counts for local slaves (NOT local master, only local slaves!)
     real(realk) :: flops_slaves
     ! Number of integral tasks
     integer :: ntasks

     ! INTEGRAL TIME ACCOUNTING
     ! ************************
     ! MPI: Time(s) used by local slaves
     real(realk) :: slavetime


  end type ccatom


  !> MP2 gradient matrices for full molecule.
  !> \author Kasper Kristensen
  !> \date October 2010
  type FullMP2grad
     !> Number of occupied orbitals in full molecule
     integer :: nocc
     !> Number of virtual orbitals in full molecule
     integer :: nvirt
     !> Number of basis functions in full molecule
     integer :: nbasis
     !> Number of atoms in full molecule
     integer :: natoms
     !> Hartree-Fock energy
     real(realk) :: EHF
     !> MP2 correlation energy
     real(realk) :: Ecorr
     !> Total MP2 energy: EHF + Ecorr
     real(realk) :: Etot
     !> MP2 correlation density matrix in AO basis (see type mp2dens)
     ! (before kappa-bar equation is solved we only have the occ-occ and virt-virt blocks,
     ! corresponding to the unrelaxed correlation density matrix)
     real(realk),pointer :: rho(:,:)
     !> Phi matrix in AO basis.
     real(realk),pointer :: Phi(:,:)
     !> In MO basis Phi is:
     !> Phivv_{ab} = sum_{cij} Theta_{cjai} (cj|bi)
     !> Phivo_{ab} = sum_{cij} Theta_{cjai} (cj|ki)
     !> Phioo_{ij} = sum_{abk} Theta_{bkai} (bk|aj)
     !> Phiov_{ic} = sum_{abk} Theta_{bkai} (bk|ac)

     !> Ltheta = sum_{aibj} Theta_{aibj} (ai|bj)^x
     real(realk),pointer :: Ltheta(:,:) 
     !> Total MP2 molecular gradient
     real(realk),pointer :: mp2gradient(:,:)
  end type FullMP2grad





  !> MP2 density matrix information for a given atomic fragment or pair fragment
  type mp2dens


     ! ************************************************************************************
     !                   MP2 correlation density matrix rho in MO basis:                  !
     !                                                                                    !
     !                               rho_{ij} = - X_{ij}                                  !
     !                               rho_{ab} = Y_{ab}                                    !
     !                               rho_{ai} = kappa_{ai}                                !
     !                               rho_{ia} = kappa_{ai}                                !
     !                                                                                    !
     ! ************************************************************************************


     ! X_{ij} = sum_{abk} t_{ki}^{ba} * mult_{kj}^{ba}
     ! Y_{ab} = sum_{cij} t_{ji}^{ca} * mult_{ji}^{cb}
     !
     ! where the multipliers can be determined simply from the amplitudes:
     ! mult_{ij}^{ab} = 4*t_{ij}^{ab} - 2*t_{ij}^{ba}
     !
     ! The determination of the kappa orbital rotation multipliers requires the solution
     ! of the full molecular orbital rotation equation. To determine the RHS matrix
     ! for this equation we need to calculate the virt-occ and occ-virt blocks of
     ! the Phi matrix (using the index convention in cc_ao_contractions):
     ! Phivo_{dl} = sum_{cij} Theta_{ij}^{cd} (ci|jl)
     ! Phiov_{lc} = sum_{abk} Theta_{kl}^{ba} (bk|ac)
     !
     ! The Phi matrix is also constructed based on fragment calculations.
     ! When X,Y,Phivo, and Phiov have been determined, the RHS for the kappa equation is:
     !
     ! RHS_{ai} = Phiov_{ia} - Phivo_{ai} + G_{ai}(M)
     !
     ! G is a Fock (Coulomb+exchange) transformation, and M is determined from X and Y:
     ! M = Y + Y^T - X - X^T
     !
     ! For the construction of M it is implicitly understood that X and Y have been
     ! transformed to the AO basis.
     ! Once RHS has been determined, kappa is found by solving:
     !
     ! E2[kappa] = RHS
     !
     ! where E2 is the Hessian transformation, see dec_solve_kappa_equation.


     ! Single fragment:
     ! ij in X_{ij} belongs to CentralAtom
     ! ab in X_{ij} belongs to CentralAtom
     !
     ! Pair fragment:
     ! i belongs to CentralAtom and j belongs to CentralAtom2 - or vice versa
     ! a belongs to CentralAtom and b belongs to CentralAtom2 - or vice versa


     !> Central atom for fragment
     integer :: CentralAtom
     !> Second central atom - only used for pair fragments
     integer :: CentralAtom2
     !> Number of basis functions in fragment
     integer :: nbasis
     !> Number of virtual AOS orbitals in fragment
     integer :: nvirt
     !> Number of occupied AOS orbitals in fragment (only valence for frozen core)
     integer :: nocc
     !> Number of occupied core+valence AOS orbitals (only different from nocc for frozen core)
     integer :: nocctot
     !> Fragment energy (for single fragment or pair fragment)
     real(realk) :: energy
     !> Only pair frags: Distance between (super) fragments in pair (zero for single fragments)
     real(realk) :: pairdist

     !> Number of EOS atoms (1 for atomic fragment, 2 for pair fragment)
     integer :: nEOSatoms
     !> List of  EOS atoms
     integer, pointer :: EOSatoms(:) => null()

     !> Indices for atomic basis functions in the list of basis functions for full molecule
     integer,pointer :: basis_idx(:) => null()

     !> Fragment component of virt-virt block of MP2 density matrix (MO basis)
     real(realk), pointer :: Y(:,:) => null()
     !> Fragment component of occ-occ block of MP2 density matrix (MO basis)
     real(realk), pointer :: X(:,:) => null()

     !> X and Y components of the MP2 correlation density matrix transformed to AO basis
     real(realk),pointer :: rho(:,:) => null()

     !> Virt-occ component of Phi matrix (needed to construct RHS for kappa-bar multiplier equation)
     !> Note: Even for frozen core the occupied index refers to both core+valence!
     real(realk), pointer :: Phivo(:,:) => null()
     !> Occ-virt component of Phi matrix (needed to construct RHS for kappa-bar multiplier equation)
     real(realk), pointer :: Phiov(:,:) => null()

  end type mp2dens



  !> Structure for fragment contribution to MP2 gradient
  type mp2grad

     ! Note: Many of the matrices needed for the MP2 gradient are also required for the MP2 density.
     ! Hence, the MP2 density structure is a subset of the MP2 gradient structure.

     !> Fragment components for MP2 correlation density matrix
     type(mp2dens) :: dens

     !> Number of atoms used to describe MOs in fragment (atomic extent)
     integer :: natoms

     !> Atomic indices for atoms in atomic extent
     integer, pointer :: atoms_idx(:) => null()

     !> Occ-occ component of Phi matrix (needed to construct reorthonormalization matrix)
     real(realk), pointer :: Phioo(:,:) => null()
     !> Virt-virt component of Phi matrix (needed to construct reorthonormalization matrix)
     real(realk), pointer :: Phivv(:,:) => null()

     !> Ltheta contribution to gradient -- dimension: (3,natoms)
     real(realk), pointer :: Ltheta(:,:) => null()

     !> Phi matrix in AO matrix (all compontents occ-occ, virt-virt, occ-virt, virt-occ)
     real(realk),pointer :: PhiAO(:,:) => null()

  end type mp2grad

  !> Batch sizes used for MP2 integral/amplitude calculation
  !> (See get_optimal_batch_sizes_for_mp2_integrals for details)
  type mp2_batch_construction
     !> Maximum allowed size of alpha batch
     integer :: MaxAllowedDimAlpha
     !> Maximum allowed size of gamma batch
     integer :: MaxAllowedDimGamma
     !> Maximum allowed size of virtual batch
     integer :: virtbatch
     !> Sizes of the four temporary arrays in step 1 of integral/amplitude scheme (AO integral part)
     integer(kind=long),dimension(4) :: size1
     !> Sizes of the four temporary arrays in step 2 of integral/amplitude scheme (virtual batch part)
     integer(kind=long),dimension(4) :: size2
     !> Sizes of the four temporary arrays in step 3 of integral/amplitude scheme (after integral loop)
     integer(kind=long),dimension(4) :: size3

  end type mp2_batch_construction

  ! Simple structure for pointer, which points to some limited chunk of larger array
  type mypointer
     !> Start index in larger array
     integer(kind=long) :: start
     !> End index in larger array
     integer(kind=long) :: end
     !> Number of elements = end-start+1
     integer(kind=long) :: N
     !> Pointer
     real(realk),pointer :: p(:) => null()
  end type mypointer


  !> Job list of fragment calculations, both single and pair
  !> Ideally they are listed in order of size with the largest jobs first.
  !> Also includes MPI performance statistics for each job.
  type joblist
     ! Number of superfragment jobs
     integer :: njobs

     ! All pointers below has the dimension njobs
     ! ------------------------------------------

     ! Atom 1 in super fragment (dimension: njobs)
     integer,pointer :: atom1(:) 
     ! Atom 2 in super fragment (dimension: njobs)   (NOTE: atom2=0 for single fragments)
     integer,pointer :: atom2(:) 
     ! Size of job (dimension: njobs)
     integer,pointer :: jobsize(:) 
     ! Is a given job done (true) or not (false) (dimension: njobs)
     logical,pointer :: jobsdone(:) 

     ! MPI statistics

     !> Number of nodes in MPI slot (local master + local slaves)
     integer,pointer:: nslaves(:)
     !> Number of occupied orbitals for given fragment (AOS)
     integer,pointer :: nocc(:)
     !> Number of virtual orbitals for given fragment (AOS)
     integer,pointer :: nvirt(:)
     !> Number of basis functions for given fragment
     integer,pointer :: nbasis(:)
     !> Number of MPI tasks used for integral/transformation (nalpha*ngamma)
     integer,pointer :: ntasks(:)
     !> FLOP count for all local nodes (local master + local slaves)
     real(realk),pointer :: flops(:)
     !> Time used for local master
     real(realk),pointer :: LMtime(:)
     !> Measure of load distribution:
     !> { (total times for nodes) / (time for local master) } / number of nodes
     real(realk),pointer :: load(:)
  end type joblist

  !> Bookkeeping when distributing DEC MPI jobs.
   TYPE traceback
      INTEGER :: na,ng,ident
   END TYPE traceback
    

   !> Integral batch handling
   TYPE batchTOorb
     INTEGER,pointer :: orbindex(:)
     INTEGER :: norbindex
  END TYPE batchTOorb

  !> \brief Grid box handling for analyzing orbitals in specific parts of space
  !> for single precision real grid points
  !> \author Kasper Kristensen
  !> \date November 2012
  type SPgridbox
     !> Center of grid box
     real(4) :: center(3)
     !> Distance between neighbouring points in grid box
     real(4) :: delta
     !> Number of grid points in each x,y,z direction measured from center
     !> (is is assumed that all sides of the grid box have the same length)
     integer :: n
     !> Number of grid points in each direction, nd = 2n+1
     !> (Somewhat redundant since it is given by n, but nice to have direct access to)
     integer :: nd
     !> Grid point values for (x,y,z) coordinates, e.g. entry(1,1,2n+1) is the
     !> point with minimum x and y values and maximum z value.
     real(4), pointer :: val(:,:,:)
  end type SPgridbox
  

  !> Information about DEC calculation
  !> We keep it as a global parameter for now.
  type(DECsettings) :: DECinfo

end module dec_typedef_module
