!> @file
!> The module contains structures used in DEC
!> \author Kasper Kristensen

!> The module contains structures used in DEC (and printing routines)
module dec_typedef_module

  use precision
  use,intrinsic :: iso_c_binding, only:c_ptr
  use TYPEDEFTYPE, only: lsitem
  use Matrix_module, only: matrix
  use tensor_type_def_module, only: tensor
  !Could someone please rename ri to something less generic. TK!!
  !  private
  !  public :: DECinfo, ndecenergies,DECsettings,array2,array3,array4,decorbital,ri,&
  !       & fullmolecule,decfrag,FullMP2grad,mp2dens,mp2grad,&
  !       & mp2_batch_construction,mypointer,joblist,traceback,batchTOorb,&
  !       & SPgridbox,MODEL_MP2,MODEL_CC2,MODEL_CCSD,MODEL_CCSDpT,MODEL_RPA,MODEL_NONE



  ! ***************************************************************************************
  !                         PARAMETERS DEFINING DEC MODELS
  ! ***************************************************************************************
  ! Do never ever use hardcoded values of these numbers inside the DEC routines!


  ! Overall CC model: MODIFY FOR NEW MODEL!
  ! ---------------------------------------
  !> how many real models in total are there, disregard MODEL_NONE
  integer,parameter :: ndecmodels   = 8
  !> Number of different fragment energies
  integer,parameter :: MODEL_NONE       = 0
  integer,parameter :: MODEL_MP2        = 1
  integer,parameter :: MODEL_CC2        = 2
  integer,parameter :: MODEL_CCSD       = 3
  integer,parameter :: MODEL_CCSDpT     = 4
  integer,parameter :: MODEL_RPA        = 5
  integer,parameter :: MODEL_RIMP2      = 6
  integer,parameter :: MODEL_SOSEX      = 7
  integer,parameter :: MODEL_LSTHCRIMP2 = 8

  ! Number of possible FOTs to consider in geometry optimization
  integer,parameter :: nFOTs=8


  ! DEC fragment energies: MODIFY FOR NEW MODEL & MODIFY FOR NEW CORRECTION
  ! -----------------------------------------------------------------------
  ! Given a CC model, there are typically more than one DEC fragment energy to consider.
  ! Parameters defining the fragment energies are given here.

  !> Number of different fragment energies
  integer, parameter :: ndecenergies = 26
  !> Numbers for storing of fragment energies in the decfrag%energies array
  integer,parameter :: FRAGMODEL_LAGMP2   = 1   ! MP2 Lagrangian partitioning scheme
  integer,parameter :: FRAGMODEL_OCCMP2   = 2   ! MP2 occupied partitioning scheme
  integer,parameter :: FRAGMODEL_VIRTMP2  = 3   ! MP2 virtual partitioning scheme
  integer,parameter :: FRAGMODEL_LAGRPA   = 4   ! MP2 Lagrangian partitioning scheme
  integer,parameter :: FRAGMODEL_OCCRPA   = 5   ! RPA occupied partitioning scheme
  integer,parameter :: FRAGMODEL_VIRTRPA  = 6   ! RPA virtual partitioning scheme
  integer,parameter :: FRAGMODEL_OCCCC2   = 7   ! CC2 occupied partitioning scheme
  integer,parameter :: FRAGMODEL_VIRTCC2  = 8   ! CC2 virtual partitioning scheme
  integer,parameter :: FRAGMODEL_OCCCCSD  = 9   ! CCSD occupied partitioning scheme
  integer,parameter :: FRAGMODEL_VIRTCCSD = 10  ! CCSD virtual partitioning scheme
  integer,parameter :: FRAGMODEL_OCCpT    = 11  ! (T) contribution, occupied partitioning scheme
  integer,parameter :: FRAGMODEL_VIRTpT   = 12  ! (T) contribution, virtual partitioning scheme
  integer,parameter :: FRAGMODEL_OCCpT4   = 13  ! Fourth order (T) contribution, occ partitioning scheme
  integer,parameter :: FRAGMODEL_VIRTpT4  = 14  ! Fourth order (T) contribution, virt partitioning scheme
  integer,parameter :: FRAGMODEL_OCCpT5   = 15  ! Fifth order (T) contribution, occ partitioning scheme
  integer,parameter :: FRAGMODEL_VIRTpT5  = 16  ! Fifth order (T) contribution, virt partitioning scheme
  integer,parameter :: FRAGMODEL_MP2f12   = 17  ! MP2-F12 energy correction
  integer,parameter :: FRAGMODEL_CCSDf12  = 18  ! CCSD-F12 energy correction
  integer,parameter :: FRAGMODEL_LAGRIMP2 = 19  ! RI-MP2 Lagrangian partitioning scheme
  integer,parameter :: FRAGMODEL_OCCRIMP2 = 20  ! RI-MP2 occupied partitioning scheme
  integer,parameter :: FRAGMODEL_VIRTRIMP2= 21  ! RI-MP2 virtual partitioning scheme
  integer,parameter :: FRAGMODEL_OCCSOS   = 22  ! SOSEX occupied partitioning scheme
  integer,parameter :: FRAGMODEL_VIRTSOS  = 23  ! SOSEX virtual partitioning scheme
  integer,parameter :: FRAGMODEL_LAGLSTHCRIMP2  = 24 ! LS-THC-RI-MP2 Lagrangian partitioning scheme
  integer,parameter :: FRAGMODEL_OCCLSTHCRIMP2  = 25 ! LS-THC-RI-MP2 occupied partitioning scheme
  integer,parameter :: FRAGMODEL_VIRTLSTHCRIMP2 = 26 ! LS-THC-RI-MP2 virtual partitioning scheme
  integer,parameter :: FRAGMODEL_RIMP2f12 = 27  ! RI-MP2F12 energy correction

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


     ! SNOOP
     ! =====
     !> Do a SNOOP calculation rather than DEC? (later SNOOP and DEC will be somewhat merged)
     logical :: SNOOP
     !> Skip CC calculation in SNOOP and just do HF
     logical :: SNOOPjustHF
     !> Maximum number of iterations in SNOOP HF calculations
     integer :: SNOOPmaxiter
     !> Convergence threshold in SNOOP HF calculations
     real(realk) :: SNOOPthr
     !> Maximum number of DIIS vectors stored in RH/DIIS scheme in SNOOP
     integer :: SNOOPmaxdiis
     !> Debug prints for SNOOP RH/DIIS
     logical :: SNOOPdebug
     !> Impose orthogonality constrant for occupied subsystem orbitals in SNOOP
     logical :: SNOOPort
     !> Use "same" orbital spaces for monomer calculation as for full calculation,
     !> as defined by natural connection
     logical :: SNOOPsamespace
     !> Localize SNOOP subsystem orbitals (cannot be used in connection with SNOOPsamespace)
     logical :: SNOOPlocalize


     ! CC response (no DEC so far)
     ! ===========================
     !> Calculate CC Jacobian eigenvalues
     logical :: CCeival
     !> Number of Jacobian eigenvalues to determine
     integer :: JacobianNumEival
     !> Use Jacobian left-transformations when determining CC eigenvalues 
     !> (false by default such that we use right-transformations)
     logical :: JacobianLHTR
     !> Convergence threshold when solving CC eigenvalue equation
     real(realk) :: JacobianThr
     !> Maximum dimension of subspace when solving Jacobian eigenvalue equation
     integer :: JacobianMaxSubspace
     !> Size of initial subspace when solving Jacobian eigenvalue equation
     integer :: JacobianInitialSubspace
     !> Maximum number of iterations when solving Jacobian eigenvalue equation
     integer :: JacobianMaxIter
     !> Use preconditioning for Jacobian eigenvalue problem
     logical :: JacobianPrecond


     !> MAIN SETTINGS DEFINING DEC CALCULATION
     !> **************************************

     !> Run DEC calculation at all?
     logical :: doDEC
     !> Frozen core calculation?
     logical :: frozencore
     !> Full molecular job
     logical :: full_molecular_cc ! full molecular cc
     !> Print fragment energies for full molecular cc
     logical :: print_frags
     !> Enforce canonical orbitals in calculation 
     logical :: use_canonical
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
     !> is the density and other matrices in the grand-canonical basis?
     logical :: gcbasis
     !> DEC-CC orbital-based (DECCO)
     logical :: DECCO
     !> DEC-CC orbital-based (DECNP)
     logical :: DECNP



     !> Restart options
     !> ***************
     !> Use HF info generated in previous run (does not necessarily require DECrestart to be true)
     logical :: HFrestart
     !> Restart DEC calculation using fragment info files (requires HFrestart to be true)
     logical :: DECrestart
     !> Enforce restart in spite of inconsistencies in restart file - only for advanced users!
     logical :: EnforceRestart
     !> Creating files for restart: Time (in seconds) passing before backing up restart files
     real(realk) :: TimeBackup
     !> Read DEC orbital file DECOrbitals.info from file (default: Existing file is overwritten)
     logical :: read_dec_orbitals

     
     !> Integral Stuff
     !> ************
     !> Screening threshold for Integral evaluation
     real(realk) :: IntegralThreshold
     !> Use Ichor Integral Code
     logical :: UseIchor


     !> Memory stuff
     !> ************
     !> Memory available for DEC calculation
     real(realk) :: memory
     !> Memory defined by input? (If not, use system call).
     logical :: memory_defined
     !> use system information to determine available memory during a dec
     !calculation
     logical :: use_system_memory_info
     !> Memory available for Background buffer in DEC calculation
     real(realk) :: bg_memory

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
     !> DEC files (lcm_orbitals.u, fock.restart, dens.restart, DECorbitals.info) 
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
     !> skip reading the old amplitudes from disk
     logical :: CCSDno_restart
     !> if mpich is used CCSD has some special treats that can be used
     logical :: CCSD_NO_DEBUG_COMM
     !> prevent canonicalization in the ccsolver
     logical :: CCSDpreventcanonical
     !> chose left-transformations to be carried out
     logical :: CCSDmultipliers
     !> use debug multiplier residual
     logical :: simple_multipler_residual
     !> use pnos in dec
     logical :: use_pnos
     !> override the transformation to the PNOs by putting unit matrices as
     !transformation matrices
     logical :: noPNOtrafo, noPNOtrunc, pno_S_on_the_fly
     logical :: noFAtrafo, noFAtrunc
     !> defines a simple cutoff threshold for constructing the PNOs from the
     !correlation density
     real(realk) :: simplePNOthr
     !> cutoff value for the overlap between different PNO spaces
     logical :: noPNOoverlaptrunc
     real(realk) :: PNOoverlapthr
     !> this defines the PNO threshold used for the EOS adapted space
     real(realk) :: EOSPNOthr
     !> use triangular counting in th occupied indices
     logical :: PNOtriangular
     !> Prevent using MO-based algorithm to solve the CCSD equations
     logical :: NO_MO_CCSD
     !> do not update the singles residual
     logical :: CCDhack
     !> Crash Calc Debug keyword - to test restart option
     logical :: CRASHCALC
     !> Crash Calc Debug keyword - to test restart option
     logical :: CRASHESTI
     !> Debug CC driver
     logical :: cc_driver_debug
     !> Use Background buffer
     logical :: use_bg_buffer
     !> Integer specifying which scheme to use in CCSD calculations (debug)
     integer :: en_mem
     !overwrite standard preconditioning settings in solver
     logical :: ccsolver_overwrite_prec
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
     !> CCSOLVER skip iterations and return starting guess
     logical :: ccsolverskip
     !> Use preconditioner in B matrix
     logical :: use_preconditioner_in_b
     !> Use CROP (if false we use DIIS)
     logical :: use_crop
     !> logial to set whether special communication processes should be spawned
     logical :: spawn_comm_proc
     !> set tilesize in ccsolver in GB
     real(realk) :: cc_solver_tile_mem
     !> select tensor segmenting scheme
     integer :: tensor_segmenting_scheme

     !> ccsd(T) settings
     !> ****************
     !> logical for abc scheme
     logical :: abc
     !> force a specific tile size for use with ijk scheme
     integer :: ijk_tile_size
     !> force a specific tile size for use with abc scheme
     integer :: abc_tile_size
     !> number of mpi buffers in ccsdpt ijk loop to prefetch tiles
     integer :: ijk_nbuffs
     !> number of mpi buffers in ccsdpt abc loop to prefetch tiles
     integer :: abc_nbuffs
     !> do we want to do gpu computations synchronous?
     logical :: acc_sync
     !> (T) hack variable - used for omitting CCSD
     logical :: pt_hack
     !> (T) hack variable - used for omitting (t) integrals (may be used in combintion with pt_hack)
     logical :: pt_hack2

     !> F12 settings
     !> ************
     !> Use F12 correction
     logical :: F12
     !> Do F12 also for fragment optimization
     logical :: F12fragopt
     !> Do C coupling in F12 scheme
     logical :: F12Ccoupling

     !> F12 debug settings
     !> ******************
     !> Use F12 correction
     logical :: F12DEBUG

     logical :: SOS

     !> Debug keyword to specify pure hydrogen atoms
     logical :: PUREHYDROGENdebug

     !> Calculate the Interaction Energy (Ref to article)
     logical :: InteractionEnergy
     !> Print the Interaction Energy (Ref to article)
     logical :: PrintInteractionEnergy

     !> Stress Test 
     logical :: StressTest
     !> Kohn-Sham Reference
     logical :: DFTreference

     !> Atomic Extent - include all atomic orbitals of atoms included
     logical :: AtomicExtent

     !> Massively parallel MP2 (Full Molecular canonical MP2)
     logical :: MPMP2
    
     !> RIMP2 settings
     !> ************
     !> Auxiliary Atomic Extent - for now include ALL atomic orbitals in RI
     logical :: AuxAtomicExtent
     !> Use Natural Auxiliary Functions (NAF)
     logical :: NAF
     !> Natural Auxiliary Functions Threshold
     real(realk) :: NAFthreshold
     !> Hardcode the Group size used in RIMP(build_CalphaMO)
     integer :: RIMPSubGroupSize
     !> Use Tensor Framework to Construct Calpha
     logical :: RIMP2PDMTENSOR
     !> Force the use of the code that distribute the Calpha
     logical :: RIMP2ForcePDMCalpha
     !> Force tiling in Step 5 of RIMP2 code
     logical :: RIMP2_tiling
     !> Use lowdin decomposition
     logical :: RIMP2_lowdin
     !> MPI group is split if #nodes > O*V/RIMPIsplit
     integer :: RIMPIsplit

     !> MPI settings
     !> ************
     !> Factor determining when MPI groups should split
     integer :: mpisplit
     !> Manually set starting group size for local MPI group
     integer(kind=ls_mpik) :: MPIgroupsize
     !> set whether to distribute the data in the full molecule structure
     logical :: distribute_fullmolecule,force_distribution

     !> Integral batching
     !> *****************
     !> Set integral batch sizes manually
     logical :: manual_batchsizes
     !> Sizes of alpha and gamma batches defined manually
     integer :: ccsdAbatch,ccsdGbatch
     !> test integral scheme, fully distributed, get_mo_integrals
     logical :: test_fully_distributed_integrals

     !> MP2 occupied batching
     !> *********************
     !> Set batch sizes manually
     logical :: manual_occbatchsizes
     !> Sizes of I and J occupied batches defined manually
     integer :: batchOccI,batchOccJ

     !> General debug and simple tests
     !> ******************************
     !> General HACK parameters, to be used for easy debugging
     logical :: hack
     logical :: hack2
     !> Skip the read-in of molecular info files dens.restart, fock.restart, lcm_orbitals.u
     logical :: SkipReadIn
     !> test the array structure
     logical :: tensor_test
     !> test the array reorderings
     logical :: reorder_test
     !> Check that LCM orbitals are correct
     logical :: check_lcm_orbitals
     !> Check the the Occupied Subsystem locality
     logical :: check_Occ_SubSystemLocality
     !> Enforce a zero Occupied Subsystem locality
     logical :: force_Occ_SubSystemLocality
     !> Debug print level
     integer :: PL
     !> reduce the output if a big calculation is done
     logical :: print_small_calc
     !> only do fragment part of density or gradient calculation 
     logical :: SkipFull 
     !> set fraction of extended orbital space to reduce to in the binary search
     real(realk) :: FracOfOrbSpace_red
     ! --

     !> Output options 
     !> **************
     !> File unit for LSDALTON.OUT
     integer :: output
     ! --

     !> DEC Orbital treatment
     !> *********************
     !> DEC is quitting after generating DEC orbitals to file
     logical :: only_generate_DECorbs
     !> Absorb H atoms into heavy atoms during orbital assignment
     logical :: AbsorbHatoms
     !> Fit orbital coefficients in fragment (default: true)
     logical :: FitOrbitals
     !> Threshold for simple Lowdin procedure for determining atomic extent
     real(realk) :: simple_orbital_threshold
     !> Purify fitted MO coefficients (projection + orthogonalization)
     logical :: PurifyMOs
     !> Use Mulliken population analysis to assign orbitals (default: Lowdin, only for Boughton-Pulay)
     logical :: mulliken
     !> Use Distance criteria to determine central atom
     logical :: Distance
     ! --


     !> Fragment optimization
     !> *********************
     !> Fragment optimization threshold
     real(realk) :: FOT
     !> Fragment optimization thresholds of decreasing magnitude (increasing accuracy)
     !> to possibly be used in geometry optimization
     real(realk) :: GeoFOTs(nFOTs)
     !> Max number of iterations for expanding fragment
     integer :: MaxIter
     !> FOT level (used for geometry opt.)
     integer :: FOTlevel
     !> Which Fragment Expansion Scheme should be used
     integer :: Frag_Exp_Scheme
     !> Which Fragment Reduction Scheme should be used for the occ space
     integer :: Frag_RedOcc_Scheme
     !> Which Fragment Reduction Scheme should be used for the vir space
     integer :: Frag_RedVir_Scheme
     !> Number of atoms to include in initial fragment
     integer :: Frag_Init_Size
     !> Number of atoms to include in fragment expansion
     integer :: Frag_Exp_Size
     !> Threshold in reduction in case of unbalanced red:
     real(realk) :: Frag_red1_thr
     real(realk) :: Frag_red2_thr
     !> space to reduce first in FOP:
     logical :: Frag_red_occ
     logical :: Frag_red_virt
     !> Model to use for fragment expansion
     integer :: fragopt_exp_model
     !> Model to use for fragment reduction
     integer :: fragopt_red_model
     !> Only consider occupied partitioning
     logical :: OnlyOccPart
     !> Only consider virtual partitioning
     logical :: OnlyVirtPart
     !> Fragment initialization radius WITHOUT OPTIMIZING THE FRAGMENT AFTERWARDS
     real(realk) :: all_init_radius
     real(realk) :: vir_init_radius
     real(realk) :: occ_init_radius
     !> Repeat atomic fragment calculations after fragment optimization?
     ! (this is necessary e.g. for gradient calculations).
     logical :: RepeatAF
     !> How to construct correlation density defining fragment-adapted orbitals?
     !> THIS IS WORK IN PROGRESS AND SHOULD BE MODIFIED!!!
     !> For now, for atomic site P: 
     !> CorrDensScheme=1:   Only EOS ampllitudes enter corr. dens . (ij \in P)
     !> CorrDensScheme=2:   EOS and EOS-coupling amplitudes         (i\in P, j \in [P] or vice versa)
     !> CorrDensScheme=3:   All AOS ampllitudes                     (ij \in [P])
     !> Note that scheme 2 is only meaningful for occupied partitioning scheme.
     integer :: CorrDensScheme
     ! --  
     logical :: use_abs_overlap
     !> Use fragment-adapted orbitals for fragment calculations
     logical :: FragAdapt
     !> avoid all pair calculations:
     logical         :: no_pairs
     !> Hack to only do fragment optimization
     integer         :: only_n_frag_jobs
     integer,pointer :: frag_job_nr(:)
     !> Use hack to specify only pair fragment jobs
     logical         :: only_pair_frag_jobs


     !> Pair fragments
     !> **************
     !> Distance cutoff for pair fragments
     real(realk) :: pair_distance_threshold
     !> When pair regression fit is performed, pair distances smaller than PairMinDist are ignored
     real(realk) :: PairMinDist
     !> Skip pair analysis (debug)
     logical :: checkpairs
     !> Fragment-adapted threshold for throwing away orbitals in atomic fragments
     !> that constitute pair fragment (currently on a testing basis)
     real(realk) :: pairFOthr
     !> Estimate pair interaction energies using simple estimates?
     logical :: PairEstimate
     !> Carry out pair estimate, but anyway run all pairs.
     logical :: PairEstimateIgnore
     !> initiation radius of the estimated fragments
     real(realk) :: EstimateINITradius
     !> number of average atoms that will be included in the estimated fragments
     integer :: EstimateInitAtom
     !> Which model to use for pair estimates
     integer :: PairEstimateModel
     !> Number of reduced fragments (increased FOT) to used for pair calculations 
     !> (NOT including fragment for main FOT)
     integer :: nFRAGSred
     !> Factor to scale FOT by for reduced fragments
     integer :: FOTscaling
     ! --


     ! First order properties
     ! **********************
     !> Do first order properties (density, electric dipole, gradient)
     logical :: first_order
     !> density matrix (and not gradient)
     logical :: density    
     !> Unrelaxed density matrix
     logical :: unrelaxed
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
     integer,dimension(nFOTs) :: ncalc
     !> Factor multiply intrinsic energy error by before returning error to geometry optimizer
     real(realk) :: EerrFactor
     !> Old energy error (used only for geometry opt)
     real(realk) :: EerrOLD


     !> Stripped down keywords - reduce memory usage to a minimum
     !> *********************************************************
     !> Should only be used by someone who knows what they are
     !> doing...
     !> Do not store AO Fock matrix in full molecule structure.
     logical :: noaofock

     !> THC keywords
     logical     :: THCNOPRUN
     logical     :: THCDUMP
     real(realk) :: THCradint
     integer     :: THC_MIN_RAD_PT
     integer     :: THCangint
     integer     :: THCHRDNES
     integer     :: THCTURBO
     integer     :: THCRADIALGRID
     logical     :: THCZdependenMaxAng
     integer     :: THCPARTITIONING
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


  type decorbital

     !> Number of the orbital in full molecular basis
     integer :: orbitalnumber
     !> Central atom to which orbital is assigned
     integer :: centralatom
     !> Number of significant atoms
!     integer :: numberofatoms
     !> Secondary central atom for orbital. For example, if virtual orbital "a" is assigned to
     !> a hydrogen atom "H_A" for which there are no occupied orbitals assigned, the secondary
     !> central atom for orbital "a" will be the atom closest to "H_A" which has a nonzero number
     !> of occupied AND virtual orbitals in the original assignment.
     integer :: secondaryatom

     !> List of significant atoms
!     integer, pointer :: atoms(:) => null()
     !> Number of significant Atomic Orbitals
     integer :: numberofaos
     !> List of significant Atomic Orbitals
     integer, pointer :: aos(:) => null()

  end type decorbital


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
  !> REMEMBER TO UPDATE mpi_bcast_fullmolecule IF YOU MODIFY THIS!!!
  type fullmolecule

     !> Number of electrons
     integer :: nelectrons
     !> Number of atoms
     integer :: natoms
     !> Number of basis functions
     integer :: nbasis
     !> Number of MOs (usually equal to nbasis but can be different for subsystems in SNOOP)
     integer :: nMO
     !> Number of auxiliary basis functions
     integer :: nauxbasis
     !> Number of occupied orbitals (core + valence)
     integer :: nocc
     !> Number of core orbitals
     integer :: ncore
     !> Number of valence orbitals (nocc-ncore)
     integer :: nval
     !> Number of virtual (virtupied) orbitals
     integer :: nvirt
     !> Number of cabs AO orbitals
     integer :: nCabsAO
     !> Number of cabs MO orbitals
     integer :: nCabsMO
     !> Number of possible fragments
     integer :: nfrags


     !> Number of basis functions on atoms
     integer, pointer :: atom_size(:) => null()
     !> Index of the first basis function for an atom
     integer, pointer :: atom_start(:) => null()
     !> Index of the last basis function for an atom
     integer, pointer :: atom_end(:) => null()
     !> Index of the first basis function for the angmom
     integer, pointer :: bas_start(:) => null()
     !> Index of the last basis function for the angmom
     integer, pointer :: bas_end(:) => null()
     !> Number of CABS basis functions on atoms
     integer, pointer :: atom_cabssize(:) => null()
     !> Index of the first CABS basis function for an atom
     integer, pointer :: atom_cabsstart(:) => null()

     !> logical that saves whether the tensors are in PDM or dense
     logical :: mem_distributed

     !> Occupied MO coefficients (mu,i)
     type(tensor) :: Co
     !> Virtual MO coefficients (mu,a)
     type(tensor) :: Cv
     !> CABS MO coefficients (mu,x)
     !     real(realk), pointer :: Ccabs(:,:) => null()
     !> RI MO coefficients 
     !     real(realk), pointer :: Cri(:,:) => null() 

     !> Fock matrix (AO basis)
     type(tensor) :: fock
     !> Occ-occ block of Fock matrix in MO basis
     type(tensor) :: oofock
     !> Virt-virt block of Fock matrix in MO basis
     type(tensor) :: vvfock

     !> Abs overlap information
     real(realk), pointer :: ov_abs_overlap(:,:) => null()
     !> carmom coord for occ
     real(realk), pointer :: carmomocc(:,:) => null()
     !> carmom coord for virt
     real(realk), pointer :: carmomvirt(:,:) => null()
     !> atomic centers
     real(realk), pointer :: AtomCenters(:,:) => null()

     !> Which atoms are phantom atoms (only basis functions)
     Logical, pointer :: PhantomAtom(:) => null()


     !> Occ-Occ Fock matrix in MO basis (change)
     real(realk), pointer :: Fij(:,:) => null()

     !> Occ-CABS (one-electron + coulomb matrix) in MO basis
     real(realk), pointer :: hJir(:,:) => null() 
     !> Cabs ri-Cabs ri exchange matrix in MO basis
     real(realk), pointer :: Krs(:,:) => null() 
     !> Cabs ri-Cabs ri Fock matrix in MO basis
     real(realk), pointer :: Frs(:,:) => null() 
     !> Virt-Cabs Fock matrix in MO basis
     real(realk), pointer :: Fac(:,:) => null() 
     !> Cabs ri-Occ Fock matrix in MO basis  
     real(realk), pointer :: Frm(:,:) => null()
     !> Cabs-(Occ+virt) Fock matrix in MO basis
     real(realk), pointer :: Fcp(:,:) => null()
     !> Cabs-Cabs  Fock matrix in MO basis
     !real(realk), pointer :: Fcd(:,:) => null()
     
     integer,pointer :: SubSystemIndex(:) => null()

     !> Pair distance table giving interatomic distances
     real(realk),pointer :: DistanceTable(:,:) => null()
     !> Table describing which model should be used for given fragment calculation:
     !> model=MODEL_NONE:  Skip fragment (only relevant for pairs)
     !> model=MODEL_MP2 :  Do MP2
     !> etc., see MODEL_* definitions below
     integer,pointer :: ccmodel(:,:) => null()
     !> FOT level to use for each pair calculation
     !>  0: Use input FOT
     !>  n>0: Use AOS information from fragment%REDfrags(n)
     integer,pointer :: PairFOTlevel(:,:) => null()

     !> Partitioning of energy into dispersion, charge transfer,
     !> and internal subsystem excitations 
     !> (see SNOOP_partition_energy).
     real(realk) :: Edisp
     real(realk) :: Ect
     real(realk) :: Esub

  end type fullmolecule


  !> Atomic fragment / Atomic pair fragment
  !> IMPORTANT: IF YOU MODIFY THIS STRUCTURE, REMEMBER TO CHANGE mpicopy_fragment ACCORDINGLY!!!
  type decfrag

     !> Number of occupied EOS orbitals 
     integer :: noccEOS=0
     !> Number of virtupied EOS orbitals 
     integer :: nvirtEOS=0
     !> Number of occupied AOS orbitals (for frozen core approx this is only the valence orbitals)
     integer,pointer :: noccAOS
     !> Number of core orbitals in AOS
     integer :: ncore=0
     !> Total number of orbitals (core+valence) in AOS (noccAOS + ncore)
     integer :: nocctot=0
     !> Total number of virtupied orbitals (AOS)
     integer,pointer :: nvirtAOS

     !> Has fragment been optimized
     logical :: isopt

     !> Pair fragment?
     logical :: pairfrag

     !> CC model to use for fragment (see MODEL_* in this file)
     integer :: ccmodel

     !> Occupied orbital EOS indices in the full basis 
     integer, pointer :: occEOSidx(:) => null()
     !> Unoccupied orbital EOS indices in the full basis 
     integer, pointer :: virtEOSidx(:) => null()
     !> Occupied AOS orbital indices (only valence orbitals for frozen core approx)
     integer, pointer :: occAOSidx(:) => null()
     !> Unoccupied AOS orbital indices in the full basis  
     integer, pointer :: virtAOSidx(:) => null()
     !> Core orbitals indices (only used for frozen core approx, 
     !> otherwise these are included in the occAOSidx list).
     integer,pointer :: coreidx(:) => null()
     !> Indices of occupied EOS in AOS basis
     integer, pointer :: idxo(:) => null()
     !> Indices of virtupied EOS in AOS basis
     integer, pointer :: idxu(:) => null()

     ! Reduced fragments (dimension DECinfo%nFOTred)
     type(fragmentAOS),pointer :: REDfrags(:)

     !> DEC fragment energies are stored in the energies array
     !> according to the global integers "FRAGMODEL_*" defined below.
     ! Note 1: Only the energies requested for the model in question are calculated!
     ! Note 2: Obviously you need to change the the global integer "ndecenergies"
     !         at the top of this file if you add new models!!!
     real(realk),dimension(ndecenergies) :: energies

     !> energy error estimates
     real(realk) :: Eocc_err
     real(realk) :: Evir_err
     real(realk) :: Elag_err

     !> Contributions to the fragment Lagrangian energy from each individual
     !  occupied or virtual orbital.
     real(realk),pointer :: OccContribs(:) => null()
     real(realk),pointer :: VirtContribs(:) => null()

     !> Number of EOS atoms (1 for atomic fragment, 2 for pair fragment)
     integer :: nEOSatoms
     !> List of EOS atoms
     integer, pointer :: EOSatoms(:) => null()


     !> Information used only when the decfrag is a pair fragment
     !> ********************************************************
     !> Distance between atomic fragments used to generate pair
     real(realk) :: pairdist
     
     !> Information about fragment size always set, this is the maximum distance
     !between any two atoms in the fragment
     real(realk) :: RmaxAE,RmaxAOS,RaveAE,RaveAOS,RsdvAE,RsdvAOS
     real(realk) :: DmaxAE,DmaxAOS,DaveAE,DaveAOS,DsdvAE,DsdvAOS

     ! NOTE!!! occAOSorb and virtAOSorb are ILL-DEFINED when fragmentadapted=.true. !!!!

     !> Total occupied orbital space (orbital type)
     type(decorbital), pointer :: occAOSorb(:) => null()
     !> Total virtupied orbital space (orbital type)
     type(decorbital), pointer :: virtAOSorb(:) => null()

     !> Number of atoms (atomic extent)
     integer :: natoms=0
     !> Number of basis functions
     integer :: nbasis=0
     !> Number of CABS basis functions
     integer :: ncabsAO=0
     !> Atomic indices
     integer, pointer :: atoms_idx(:) => null()
     !> Corresponding basis function indices
     integer,pointer :: basis_idx(:) => null()
     !> Corresponding CABS basis function indices
     integer,pointer :: cabsbasis_idx(:) => null()

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

     ! Note: Co and Cv will point to CoLOC and CvLOC if local orbitals are used
     !>      (or whatever the input orbitals are)   OR
     !>      Co and Cv will point to CoFA and CvFA (when FO=.true.)

     !> Occupied MO coefficients (only valence space for frozen core approx)
     real(realk), pointer :: Co(:,:) => null()
     !> Virtual MO coefficients
     real(realk), pointer :: Cv(:,:) => null()
     !> Cabs MO coefficients
     real(realk),pointer :: Ccabs(:,:) => null()     
     !> Core MO coefficients 
     real(realk),pointer :: CoreMO(:,:) => null()
     !> RI Mo coefficients
     real(realk),pointer :: Cri(:,:) => null()


     !> AO Fock matrix
     real(realk), pointer :: fock(:,:) => null()
     !> Occ-occ block of Fock matrix in MO basis  (only valence space for frozen core approx)
     real(realk), pointer :: ppfock(:,:) => null()
     !> Virt-virt block of Fock matrix in MO basis
     real(realk), pointer :: qqfock(:,:) => null()
     !> Core-core block of Fock matrix in MO basis  (subset of ppfock when frozen core is NOT used)
     real(realk), pointer :: ccfock(:,:) => null()


     !> Occ-Occ Fock matrix in MO basis
     real(realk), pointer :: Fij(:,:) => null()

     !> Occ-CABS (one-electron + coulomb matrix) in MO basis
     real(realk), pointer :: hJir(:,:) => null() 
     !> Cabs ri-Cabs ri exchange matrix in MO basis
     real(realk), pointer :: Krs(:,:) => null() 
     !> Cabs ri-Cabs ri Fock matrix in MO basis
     real(realk), pointer :: Frs(:,:) => null() 
     !> Virt-Cabs Fock matrix in MO basis
     real(realk), pointer :: Fac(:,:) => null() 
     !> Cabs ri-Occ Fock matrix in MO basis  
     real(realk), pointer :: Frm(:,:) => null()
     !> Cabs-(Occ+virt) Fock matrix in MO basis
     real(realk), pointer :: Fcp(:,:) => null()
    
     ! Information for local orbitals
     ! ******************************
     !> Local occupied MO coefficients
     real(realk), pointer :: CoLOC(:,:) => null()
     !> Local virtual MO coefficients
     real(realk), pointer :: CvLOC(:,:) => null()
     !> Occ-occ block of Fock matrix in local MO basis  (only valence space for frozen core approx)
     real(realk), pointer :: ppfockLOC(:,:) => null()
     !> Virt-virt block of Fock matrix in local MO basis
     real(realk), pointer :: qqfockLOC(:,:) => null()

     ! Information used for fragment-adapted orbitals
     ! **********************************************
     !> Correlation density matrices in local AOS basis
     real(realk), pointer :: OccMat(:,:) => null()  ! occ AOS-EOS
     real(realk), pointer :: VirtMat(:,:) => null()  ! virt AOS-EOS
     !> Threshold to use for throwing away fragment-adapted occupied (1) or virtual (2) orbitals
     !> For pair fragment PQ this refers to the threshold used to throw away FOs from
     !> the underlying atomic fragments P and Q.
     real(realk) :: RejectThr(2)
     !> Control of whether corr dens matrices have been set (true) or simply initialized (false)
     logical :: CDset
     !> Is this a fragment-adapted fragment?
     logical :: fragmentadapted
     !> Number of occ orbitals for fragment-adapted orbitals 
     integer,pointer :: noccFA
     !> Number of virt orbitals for fragment-adapted orbitals 
     integer,pointer :: nvirtFA
     !> Transformation between AO basis and fragment-adapted basis
     !> Index 1: Local,   Index 2: Fragment-adapted
     !> Has fragment-adapted MO coeff been set (not done by default fragment initialization)?
     logical :: FAset
     !> Occupied FA coeff
     real(realk),pointer :: CoFA(:,:) => null()     ! dimension: nbasis,noccFA
     !> Virtual FA coeff
     real(realk),pointer :: CvFA(:,:) => null()   ! dimension: nbasis,nvirtFA
     !> Eigenvalues for correlation density matrices 
     !> --> only set for atomic fragments (pairfrag=.false.) and when FAset=.true.
     real(realk),pointer :: CDocceival(:) => null()    ! dimension noccFA
     real(realk),pointer :: CDvirteival(:) => null()  ! dimension nvirtFA
     !> Occ-occ block of Fock matrix in FO basis  (only valence space for frozen core approx)
     real(realk), pointer :: ppfockFA(:,:) => null()
     !> Virt-virt block of Fock matrix in FO basis
     real(realk), pointer :: qqfockFA(:,:) => null()

     ! Information used for pair-natural orbitals (this only applies to virtual
     ! oribtals, the occupied orbitals will be kept in the local basis)
     ! ******************************************
     !> use PNO information?
     logical :: PNOset
     !> collection of transformation matrices from LO space to PNO space
     type(PNOSpaceInfo), pointer :: CLocPNO(:)
     !> number of spaces to consider
     integer :: nspaces


     !> Integral program input
     type(lsitem) :: mylsitem

     ! End of EXPENSIVE BOX
     ! ==============================================================


     ! Information for local orbitals
     ! ******************************
     !> Number of local occupied orbitals in fragment
     integer,pointer :: noccLOC
     !> Number of local virtupied orbitals in fragment
     integer,pointer :: nvirtLOC


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
     ! MPI: Sum of GPU flop counts for local slaves (NOT local master, only local slaves!)
     real(realk) :: gpu_flops_slaves
     ! Number of integral tasks
     integer :: ntasks

     ! INTEGRAL TIME ACCOUNTING
     ! ************************
     ! MPI: Time(s) used by local slaves
     real(realk),dimension(ndecmodels) :: slavetime_work
     real(realk),dimension(ndecmodels) :: slavetime_comm
     real(realk),dimension(ndecmodels) :: slavetime_idle


  end type decfrag


  !> MP2 gradient matrices for full molecule.
  !> \author Kasper Kristensen
  !> \date October 2010
  type FullMP2grad
     !> Number of occupied orbitals in full molecule
     integer :: nocc
     !> Number of virtupied orbitals in full molecule
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


     ! atomic fragment:
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
     !> Number of virtupied AOS orbitals in fragment
     integer :: nvirt
     !> Number of occupied AOS orbitals in fragment (only valence for frozen core)
     integer :: nocc
     !> Number of occupied core+valence AOS orbitals (only different from nocc for frozen core)
     integer :: nocctot
     !> Fragment energy (for atomic fragment or pair fragment)
     real(realk) :: energy
     !> Only pair frags: Distance between fragments in pair (zero for atomic fragments)
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
     ! Number of fragment jobs
     integer :: njobs

     ! All pointers below has the dimension njobs
     ! ------------------------------------------

     ! Atom 1 in fragment (dimension: njobs)
     integer,pointer :: atom1(:) 
     ! Atom 2 in fragment (dimension: njobs)   (NOTE: atom2=0 for atomic fragments)
     integer,pointer :: atom2(:) 
     ! Size of job (dimension: njobs)
     integer,pointer :: jobsize(:) 
     ! Is a given job done (true) or not (false) (dimension: njobs)
     logical,pointer :: jobsdone(:) 
     ! Does job require fragment optimization?
     logical,pointer :: dofragopt(:)
     ! Does job use estimated fragments?
     logical,pointer :: esti(:)
     !> Number of occupied orbitals for given fragment (AOS)
     integer,pointer :: nocc(:)
     !> Number of virtupied orbitals for given fragment (AOS)
     integer,pointer :: nvirt(:)
     !> Number of basis functions for given fragment
     integer,pointer :: nbasis(:)

     ! MPI statistics

     !> Number of nodes in MPI slot (local master + local slaves)
     integer,pointer:: nslaves(:)
     !> Number of MPI tasks used for integral/transformation (nalpha*ngamma)
     integer,pointer :: ntasks(:)
     !> FLOP count for all local nodes (local master + local slaves)
     real(realk),pointer :: flops(:)
     !> GPU FLOP count for all local nodes (local master + local slaves)
     real(realk),pointer :: gpu_flops(:)
     !> Time used for local master
     real(realk),pointer :: LMtime(:)
     !> Measure of load distribution:
     !> ( work and communication times for nodes) / {(time for local master) * (number of nodes) }
     !> ( work times for nodes) / {(time for local master) * (number of nodes) }
     real(realk),pointer :: commt(:)
     real(realk),pointer :: workt(:)
     real(realk),pointer :: idlet(:)
  end type joblist

  !> Bookkeeping when distributing DEC MPI jobs.
  type traceback
     integer :: na,ng,ident
  end type traceback

  type int_batch
     integer :: nbatches,max_dim
     integer, pointer :: orb2batch(:)
     integer, pointer :: batchdim(:)
     integer, pointer :: batchsize(:)
     integer, pointer :: batchindex(:)
     type(batchtoorb), pointer :: batch2orb(:)
  end type int_batch

  !> Integral batch handling
  type batchTOorb
     integer,pointer :: orbindex(:)
     integer :: norbindex
  end type batchTOorb

  !> MO Integral batch info:
  type MObatchInfo

     !> number of batches:
     integer :: nbatch
     !> dimension of each of the nbatch1:
     integer, pointer :: dimInd1(:) 
     !> dimension of each of the nbatch2:
     integer, pointer :: dimInd2(:)
     !> MO index corresponding to the starting point of each batch:
     integer, pointer :: StartInd1(:) 
     integer, pointer :: StartInd2(:) 
     !> Total dimension of the batch
     integer, pointer :: dimTot(:)
     !> Tile index for pdm arrays
     integer, pointer :: tileInd(:,:)

  end type MObatchInfo

  !> AO Integral batch info:
  type DecAObatchinfo
     integer :: Dim       ! Dimension of DEC batch of AO batches 
     integer :: OrbStart  ! First orbital index in DEC batch
     integer :: OrbEnd    ! Last orbital index in DEC batch
     integer :: AOStart   ! First AO batch index in DEC batch
     integer :: AOEnd     ! Last AO batch index in DEC batch
  end type DecAObatchinfo

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


  !> Space information specifically designed to keep PNO spaces nicely together
  !> IMPORTANT:
  !> if you modify this structure, also modify bufferadd_PNOSpaceInfo_struct in decmpi.F90
  type PNOSpaceInfo

     integer              :: rpd                    ! corresponding dimension of the (restricted) pair space [i(<=)j], read reduced pair dimension

     integer              :: n                      ! number of occ orbitals in the corresponding space
     integer, pointer     :: iaos(:) => null()      ! orbital index in the aos space

     logical              :: s_associated           ! indicate whether s matrices are associated, i.e. a SVD  d = s1 d_new  s2^T  
     integer              :: ns1,ns2,red1,red2      ! dimensions, depending on s_associated
     real(realk), pointer :: d(:,:)  => null()      ! density matrix or overlap matrix. if s_associated d = d_new, either (ns1,ns2) or (red1,red2)
     real(realk), pointer :: s1(:,:) => null()      ! the left unit matrix reduced to the kernel dimensions (ns1,red1)
     real(realk), pointer :: s2(:,:) => null()      ! the right unit matrix reduced to the kernel dimensions (red2,ns2)

     logical              :: allocd,contributes     ! logical to show the allocation and contribution status
     logical              :: is_FA_space            ! save whether this refers to FA space, only important for trafo mats, not for overlap
     logical              :: PS                     ! save wheter it is a triangular pair space
  end type PNOSpaceInfo

  type pno_query_info
     integer(kind=8)          :: n_arrays
     integer(kind=8), pointer :: size_array(:)
  end type pno_query_info

  ! Information for fragment AOS (intended to be used for fragments of reduced FOT)
  ! Remember to modify mpicopy_fragmentAOStype if you add/remove someting here!
  type fragmentAOS
     !> Number of occupied and virtupied orbitals in fragment AOS
     integer          :: noccAOS, nvirtAOS
     !> Occupied and virtupied AOS indices
     integer, pointer :: occAOSidx(:), virtAOSidx(:)
     !> FOT corresponding to orbital spaces
     real(realk) :: FOT
  end type fragmentAOS

end module dec_typedef_module
