# WARNING if you add a new set of sources that contain .F90 files 
# remember to add it to the collection of all free fortran sources
# LSDALTON_FREE_FORTRAN_SOURCES  
# at the end of the file

set(DEC_C_SOURCES
    deccc/crayio.c
    )
set(DFTFUNC_SOURCES
    dft/fun-example.c
    dft/fun-b97-1.c
    dft/fun-becke.c
    dft/fun-cam-b3lyp.c
    dft/fun-camx.c
    dft/fun-cam_compx.c
    dft/fun-gga.c
    dft/fun-hcth120.c
    dft/fun-hcth147.c
    dft/fun-hcth407.c
    dft/fun-hcth93.c
    dft/fun-kt.c
    dft/fun-lb94.c
    dft/fun-lyp.c
    dft/fun-optx.c
    dft/fun-p86c.c
    dft/fun-pbec.c
    dft/fun-pbex.c
    dft/fun-pw86x.c
    dft/fun-pw91c.c
    dft/fun-pz81.c
    dft/fun-slater.c
    dft/fun-vwn.c
    dft/functionals.c
    dft/general.c
    )
set(DFTFUNC_F_SOURCES
    dft/II_util.F90
    )
set(FMM_C_SOURCES
    mm/mm_proc_selector.c
    )
set(LSDALTONMAIN_FORTRAN_SOURCES
    lsdaltonsrc/numerical_derivatives.F90
    lsdaltonsrc/LSDALTON.F90
    lsdaltonsrc/LS_optimizer.F90
    lsdaltonsrc/LSdynamics.F90
    lsdaltonsrc/Energy_and_deriv.F90
    lsdaltonsrc/lsmpiSlave.F90
    lsdaltonsrc/init_lsdalton.F90
    lsdaltonsrc/configuration.F90
    lsdaltonsrc/LSlib.F90
    lsdaltonsrc/Profile.F90
    )
set(DDYNAM_SOURCES
    ddynam/LSinput.F90
    ddynam/Fock_mat_dyn.F90
    ddynam/Modules.F90
    ddynam/Temperature.F90
    ddynam/TimeRev_prop.F90
    ddynam/dyn_utilities.F90
    )
set(DEC_SOURCES
    deccc/array4_memory.F90
    deccc/array3_memory.F90
    deccc/ccarray3_simple.F90
    deccc/CABS.F90
    deccc/mp2.F90
    deccc/ccsdpt.F90
    deccc/cc_crop.F90
    deccc/ccsd.F90
    deccc/rpa.F90
    deccc/f12_integrals.F90
    deccc/cc_driver.F90
    deccc/cc_integrals.F90
    deccc/ccarray2_simple.F90
    deccc/ccarray4_simple.F90
    deccc/ccorbital.F90
    deccc/ccri_simple.F90
    deccc/dec_atom.F90
    deccc/dec_driver.F90
    deccc/dec_driver_slave.F90
    deccc/dec_main.F90
    deccc/dec_settings.F90
    deccc/dec_utils.F90
    deccc/full_driver_f12contractions.F90
    deccc/fullmolecule.F90
    deccc/mp2_gradient.F90
    deccc/fragment_energy.F90
    deccc/full_driver.F90
    deccc/decmpi.F90
    deccc/decmpiSlave.F90
    )	
set(GEOOPT_SOURCES
    geomopt/LSopt-input.F90
    geomopt/ls_opt.F90
    geomopt/ls_opt2.F90
    geomopt/ls_redint.F90
    )
set(LINEARS_SOURCES	
    linears/configurationType.F90
    linears/ChargeLoc.F90
    linears/SCFLOOP.F90
    linears/trustradius.F90
    linears/LSDALTON_RESPONSE.F90
    linears/LSDALTON_RESPONSE_type.F90
    linears/arh_debug.F90
    linears/arh_driver.F90
    linears/average_util.F90
    linears/davidson_settings.F90
    linears/davidson_solver.F90
    linears/debug.F90
    linears/dens_subspace.F90
    linears/densopt.F90
    linears/diag.F90
    linears/diis.F90
    linears/dsm.F90
    linears/dsm_xterm.F90
    linears/ecdata_module.F90
    linears/extra-output.F90
    linears/fock-eval.F90
    linears/kurtosis.F90
    linears/leastchange.F90
    linears/levelshift.F90
    linears/minimize.F90
    linears/optimlocNONMOD.F90
    linears/OrbLoc_input.F90
    linears/localization_util.F90
    linears/localization_orbspread.F90
    linears/localization_charge.F90
    linears/prop_contribs.F90
    linears/response_driver.F90
    linears/response_prop.F90
    linears/rsp_equations.F90
    linears/scfopt-typedef.F90
    linears/soeo-redspace.F90
    linears/soeo-loop.F90
    linears/soeo-debug.F90
    linears/soeo-matop.F90
    linears/soeo-typedef.F90
    linears/soeo-util.F90
    linears/starting_guess.F90
    linears/statistics.F90
    linears/trilevel.F90
    linears/plt_driver.F90
    )

set(RSPSOLVER_SOURCES	
    responsesolver/rsp_cmplx_sym.F90
    responsesolver/rsp_complex.F90
    responsesolver/rsp_solver.F90
    responsesolver/rsp_std_sym.F90
    )

set(SOLVERUTIL_SOURCES	
    SolverUtilities/arh_density.F90
    SolverUtilities/dalton_interface.F90
    SolverUtilities/dd_utilities.F90
    SolverUtilities/dd_utilities_unres.F90
    SolverUtilities/decomp.F90
    SolverUtilities/lsdalton_modules.F90
    SolverUtilities/queue-operations.F90
    SolverUtilities/queue-typedef.F90
    SolverUtilities/rsp_precond.F90
    SolverUtilities/rsp_utilities.F90
    )

set(RSP_PROPERTIES_SOURCES	
    rsp_properties/molecular_hessian.F90
    rsp_properties/test_molHessian.F90
    )

set(PBC_FORTRAN_SOURCES
    pbc2/pbc_compare.F90
    pbc2/pbc-multipole.F90
    pbc2/pbc-matop.F90
    pbc2/pbc-harmonics.F90
    pbc2/pbc-msc.F90
    pbc2/pbc-krsp-op.F90
    pbc2/pbc_kscf.F90
    pbc2/pbc_mmit.F90
    pbc2/pbc_eigsolv.F90
    pbc2/pbc-data.F90
    pbc2/pbc-ffdata.F90
    pbc2/pbc_int.F90
    )

set(FMM_SOURCES
    mm/mm_T_contractors.F90
    mm/mm_T_pair_builder.F90
    mm/mm_T_pair_tests.F90
    mm/mm_T_worker.F90
    mm/mm_Vff_driver.F90
    mm/mm_Vff_processor.F90
    mm/mm_W_contractors.F90
    mm/mm_W_pair_builder.F90
    mm/mm_W_worker.F90
    mm/mm_box_builder.F90
    mm/mm_box_processor.F90
    mm/mm_buffers.F90
    mm/mm_densfit.F90
    mm/mm_driver.F90
    mm/mm_global_paras.F90
    mm/mm_grid_searcher.F90
    mm/mm_input.F90
    mm/mm_interface.F90
    mm/mm_map_builder.F90
    mm/mm_memory.F90
    mm/mm_memory_manager.F90
    mm/mm_overlap_tree.F90
    mm/mm_sort.F90
    mm/mm_stats.F90
    mm/mm_tree_T_buffer.F90
    mm/cbifmm_f90.F90
    )   
set(INTERESTLIB_SOURCES
    interest/src/module_interest.f90
    interest/src/module_interest_hrr.f90
    interest/src/module_interest_osr.f90
    )
set(LSINT_SOURCES
    LSint/dft_gridLL.F90
    LSint/BuildBasis.F90
    LSint/BuildMolFile.F90
    LSint/II_Fragment.F90
    LSint/gridgeneration.F90
    LSint/gridgeneration_boxify.F90
    LSint/II_XC_interface.F90
    LSint/II_absval_int.F90
    LSint/II_dft_int.F90
    LSint/II_dft_ksm.F90
    LSint/II_dft_ksm_worker.F90
    LSint/Integral_interface.F90
    LSint/Integral_interfaceDEC.F90
    LSint/Integral_interfaceDF.F90
    LSint/Integral_interface_test.F90
    LSint/MBIE.F90
    LSint/Molecule.F90
    LSint/ODbatches.F90
    LSint/ReadMolFile.F90
    LSint/ThermiteDistribute.F90
    LSint/ThermiteDistributeK.F90
    LSint/ThermiteDistributeDEC.F90
    LSint/ThermiteDriver.F90
    LSint/ThermiteIntegrals.F90
    LSint/ThermiteOverlapDistribution.F90
    LSint/ThermiteProp.F90
    LSint/ThermiteMem.F90
    LSint/dalton-input.F90
    LSint/linsolv_df.F90
    LSint/ls_IntegralInterface.F90
    LSint/pari.F90
    LSint/lsmpi.F90
    )
#####################################################
#WARNING: READ ME BEFORE ADDING FILES TO LSUTIL
# lsutil is split into several bunches
# The LSUTIL_PRECISION_SOURCES is the first to be compiled
# and contains only lsutil/ls_precision.F90
# the next be compiled is LSUTIL_MATRIXM_SOURCES
# which contains the matrix type definition. 
# then LSUTIL_COMMON_SOURCES which should only contain 
# type definitions and files which only depend on 
# lsutil/ls_precision.F90 and lsutil/matrix_module.F90
# then LSUTIL_MATRIXO_SOURCES and LSUTIL_MATRIXO_C_SOURCES and
# LSUTIL_MATRIXU_SOURCES is compiled. Then
# LSUTIL_TYPE_SOURCES which despite the name
# contains alot of files which deals wiht Operations of 
# types. Then at the end some Utillities are gathered in
# LSUTILLIB_SOURCES
############################################################

set(LSUTIL_PRECISION_SOURCES
    lsutil/ls_precision.F90
    lsutil/ptr_assoc.F90
    )

set(LSUTIL_MATRIXM_SOURCES
    lsutil/matrix_module.F90    
    )

set(LSUTIL_COMMON_SOURCES
    lsutil/rsp-typedef.F90
    lsutil/tensor_type_def.F90
    lsutil/response_prop_type.F90
    lsutil/ls_IOType.F90
    lsutil/AOType.F90
    lsutil/dec_typedef.F90
    lsutil/MoleculeType.F90
    lsutil/ODType.F90
    lsutil/lstensorType.F90
    lsutil/dftType.F90
    lsutil/BasisinfoType.F90
    lsutil/IntegralOutputType.F90 
    lsutil/f12.F90
    lsutil/IntegralType.F90
    lsutil/TYPE-DEF.F90
    lsutil/memory.F90
    lsutil/gridgeneration_memory.F90
    lsutil/Time.F90
    lsutil/common_utilities.F90
    lsutil/ls_parameters.F90
    lsutil/par_mod.F90
    lsutil/lsmpiType.F90
    lsutil/TestMPI.F90
    lsutil/LSmatrixType.F90
    lsutil/lsmatop_dense.F90
    lsutil/file-operations.F90
    lsutil/grid_utilities.F90
    lsutil/mat3dim.F90
    lsutil/papi.F90
    lsutil/ls_math.F90
    lsutil/SphCartMatrices.F90
    lsutil/OverlapDistributionType.F90
    )

set(LSUTIL_MATRIXO_SOURCES
    lsutil/matop_csr.F90
    lsutil/matop_dense.F90
    lsutil/matop_dense_unrest.F90
    lsutil/matop_scalapack.F90    
    )
set(LSUTIL_MATRIXO_C_SOURCES
    lsutil/matop_csr_aux.c
    lsutil/myPAPI_set_inherit.c
    )

set(LSUTIL_MATRIXU_SOURCES
    lsutil/mat-operations-aux.F90
    lsutil/mat-operations-essential.F90
    lsutil/matrix_utilities.F90
    lsutil/matrix_defop_backend.F90
    lsutil/matrix_defop_lowlevel.F90
    lsutil/matrix_defop.F90
    )

set(LSUTIL_TYPE_SOURCES
    lsutil/AO_operations.F90
    lsutil/Molecule_operations.F90
    lsutil/OD_operations.F90
    lsutil/lstensor_Mem.F90
    lsutil/lstensor_operations.F90
    lsutil/Integral_input.F90
    lsutil/ls_IO.F90
    lsutil/ls_screen.F90
    lsutil/dft_operations.F90
    lsutil/dftMem.F90
    lsutil/Basisinfo_operations.F90
    lsutil/IntegralOutput_operations.F90 
    lsutil/TYPE-OP.F90
    lsutil/pbc_lattice_type.F90
    lsutil/Build_AOBATCH.F90
    lsutil/lspdm_basic.F90
    lsutil/tensor_basic.F90
    lsutil/lspdm_tensor_operations.F90
    lsutil/manual_reorderings.F90
    lsutil/manual_utils.F90
    )

set(LSUTILLIB_SOURCES
    lsutil/lowdin.F90
    lsutil/AtomicSparse.F90
    lsutil/fundamental.F90
    lsutil/ks-settings.F90
    lsutil/lsmpi-operations.F90
    lsutil/lsutilities.F90
    lsutil/ls_screenMPI.F90
    lsutil/LSoptType.F90
    lsutil/ddynType.F90
    lsutil/ProfileType.F90
    lsutil/pbc_lattice_vectors.F90
    lsutil/tensor_interface.F90
    lsutil/lspdm_slave.F90
    )   

set(LSLIB_SOURCES
    lsdaltonsrc/LSlib.F90
    lsdaltonsrc/LSlib_tester.F90
    )

set(LSDALTON_FIXED_FORTRAN_SOURCES
    pdpack/jacobi.F
    pdpack/eispack.F
    pdpack/linpack.F
    pdpack/cholesky.F
    pdpack/invroutines.F
    )
set(LSDALTON_OWN_BLAS_SOURCES
    pdpack/blas_wrapper.F
    pdpack/gp_dblas1.F
    pdpack/gp_dblas2.F
    pdpack/gp_dblas3.F
    pdpack/gp_sblas.F
    pdpack/gp_qblas.F
    pdpack/gp_zblas.F
    )
set(LSDALTON_OWN_LAPACK_SOURCES
    pdpack/gp_dlapack.F
    pdpack/gp_zlapack.F
    )
#collection of all free fortran sources
set(LSDALTON_FREE_FORTRAN_SOURCES
    ${DFTFUNC_F_SOURCES}
    ${LSDALTONMAIN_FORTRAN_SOURCES}
    ${DDYNAM_SOURCES}
    ${DEC_SOURCES}
    ${GEOOPT_SOURCES}
    ${LINEARS_SOURCES}	
    ${RSPSOLVER_SOURCES}	
    ${SOLVERUTIL_SOURCES}	
    ${RSP_PROPERTIES_SOURCES}	
    ${PBC_FORTRAN_SOURCES}
    ${FMM_SOURCES}
    ${INTERESTLIB_SOURCES}
    ${LSINT_SOURCES}
    ${LSUTIL_PRECISION_SOURCES}
    ${LSUTIL_MATRIXM_SOURCES}
    ${LSUTIL_COMMON_SOURCES}
    ${LSUTIL_MATRIXO_SOURCES}
    ${LSUTIL_MATRIXO_C_SOURCES}
    ${LSUTIL_MATRIXU_SOURCES}
    ${LSUTIL_TYPE_SOURCES}
    ${LSUTILLIB_SOURCES}
    ${LSLIB_SOURCES}
 )