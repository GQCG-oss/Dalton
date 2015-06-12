# WARNING if you add a new set of sources that contain .F90 files 
# remember to add it to the collection of all free fortran sources
# LSDALTON_FREE_FORTRAN_SOURCES  
# at the end of the file

set(CUDA_GPU_INTERFACE_SOURCES
    LSDALTON/cuda/gpu_interfaces.F90
    )
set(DFTFUNC_SOURCES
    LSDALTON/dft/fun-example.c
    LSDALTON/dft/fun-b97-1.c
    LSDALTON/dft/fun-becke.c
    LSDALTON/dft/fun-cam-b3lyp.c
    LSDALTON/dft/fun-camx.c
    LSDALTON/dft/fun-cam_compx.c
    LSDALTON/dft/fun-gga.c
    LSDALTON/dft/fun-hcth120.c
    LSDALTON/dft/fun-hcth147.c
    LSDALTON/dft/fun-hcth407.c
    LSDALTON/dft/fun-hcth93.c
    LSDALTON/dft/fun-kt.c
    LSDALTON/dft/fun-lb94.c
    LSDALTON/dft/fun-lyp.c
    LSDALTON/dft/fun-optx.c
    LSDALTON/dft/fun-p86c.c
    LSDALTON/dft/fun-pbec.c
    LSDALTON/dft/fun-pbex.c
    LSDALTON/dft/fun-pw86x.c
    LSDALTON/dft/fun-pw91c.c
    LSDALTON/dft/fun-pz81.c
    LSDALTON/dft/fun-slater.c
    LSDALTON/dft/fun-vwn.c
    LSDALTON/dft/fun-revpbex.c
    LSDALTON/dft/fun-rpbex.c
    LSDALTON/dft/fun-mpbex.c
    LSDALTON/dft/fun-pw91x.c
    LSDALTON/dft/fun-g96.c
    LSDALTON/dft/fun-lg93.c
    LSDALTON/dft/functionals.c
    LSDALTON/dft/general.c
    )
set(DFTFUNC_F_SOURCES
    LSDALTON/dft/II_util.F90
    )
set(FMM_C_SOURCES
    LSDALTON/mm/mm_proc_selector.c
    )
set(LSDALTONMAIN_FORTRAN_SOURCES
    LSDALTON/lsdaltonsrc/numerical_derivatives.F90
    LSDALTON/lsdaltonsrc/LSDALTON.F90
    LSDALTON/lsdaltonsrc/LS_optimizer.F90
    LSDALTON/lsdaltonsrc/LSdynamics.F90
    LSDALTON/lsdaltonsrc/Energy_and_deriv.F90
    LSDALTON/lsdaltonsrc/lsmpiSlave.F90
    LSDALTON/lsdaltonsrc/init_lsdalton.F90
    LSDALTON/lsdaltonsrc/configuration.F90
    LSDALTON/lsdaltonsrc/LSlib.F90
    LSDALTON/lsdaltonsrc/LSlibState.F90
    LSDALTON/lsdaltonsrc/Profile.F90
    LSDALTON/lsdaltonsrc/IchorTesting.F90
    LSDALTON/lsdaltonsrc/IchorProfile.F90
    LSDALTON/lsdaltonsrc/InteractionEnergy.F90
    LSDALTON/lsdaltonsrc/ADMMbasisOpt.F90
    )
set(DDYNAM_SOURCES
    LSDALTON/ddynam/LSinput.F90
    LSDALTON/ddynam/Fock_mat_dyn.F90
    LSDALTON/ddynam/Modules.F90
    LSDALTON/ddynam/Temperature.F90
    LSDALTON/ddynam/TimeRev_prop.F90
    LSDALTON/ddynam/dyn_utilities.F90
    )
set(DEC_SOURCES
    LSDALTON/deccc/array4_memory.F90
    LSDALTON/deccc/array3_memory.F90
    LSDALTON/deccc/ccarray3_simple.F90
    LSDALTON/deccc/CABS.F90
    LSDALTON/deccc/mp2.F90
    LSDALTON/deccc/rimp2.F90
    LSDALTON/deccc/ri_util.F90
    LSDALTON/deccc/ls_thc_rimp2.F90
    LSDALTON/deccc/ccsdpt.F90
    LSDALTON/deccc/crop_tools.F90
    LSDALTON/deccc/cc_tools.F90
    LSDALTON/deccc/cc_response_tools.F90
    LSDALTON/deccc/ccsd.F90
    LSDALTON/deccc/pno_ccsd.F90
    LSDALTON/deccc/rpa.F90
    LSDALTON/deccc/f12_integrals.F90
    LSDALTON/deccc/f12_routines.F90
    LSDALTON/deccc/cc_driver.F90
    LSDALTON/deccc/cc_integrals.F90
    LSDALTON/deccc/ccarray2_simple.F90
    LSDALTON/deccc/ccarray4_simple.F90
    LSDALTON/deccc/ccorbital.F90
    LSDALTON/deccc/dec_atom.F90
    LSDALTON/deccc/dec_driver.F90
    LSDALTON/deccc/dec_driver_slave.F90
    LSDALTON/deccc/dec_main.F90
    LSDALTON/deccc/dec_settings.F90
    LSDALTON/deccc/dec_utils.F90
    LSDALTON/deccc/dec_tools.F90
    LSDALTON/deccc/full_driver_f12contractions.F90
    LSDALTON/deccc/fullmolecule.F90
    LSDALTON/deccc/mp2_gradient.F90
    LSDALTON/deccc/ccsd_gradient.F90
    LSDALTON/deccc/fragment_energy.F90
    LSDALTON/deccc/full_driver.F90
    LSDALTON/deccc/full_rimp2.F90
    LSDALTON/deccc/full_rimp2f12.F90
    LSDALTON/deccc/full_mp2.F90
    LSDALTON/deccc/full_ls_thc_rimp2.F90
    LSDALTON/deccc/snoop_main.F90
    LSDALTON/deccc/snoop_tools.F90
    LSDALTON/deccc/decmpi.F90
    LSDALTON/deccc/decmpiSlave.F90
    )
set(GEOOPT_SOURCES
    LSDALTON/geomopt/LSopt-input.F90
    LSDALTON/geomopt/ls_opt.F90
    LSDALTON/geomopt/ls_opt2.F90
    LSDALTON/geomopt/ls_redint.F90
    LSDALTON/geomopt/q_to_x.F90
    LSDALTON/geomopt/dqdx.cpp           
    )
set(LSDALTON_PCM_SOURCES	
    LSDALTON/pcm/ls_pcm_config.F90
    LSDALTON/pcm/ls_pcm_utils.F90
    LSDALTON/pcm/ls_pcm_write.F90
    LSDALTON/pcm/ls_pcm_integrals.F90
    LSDALTON/pcm/ls_pcm_scf.F90
    )
set(LINEARS_SOURCES	
    LSDALTON/linears/configurationType.F90
    LSDALTON/linears/LocTypes.F90
    LSDALTON/linears/ChargePrec.F90
    LSDALTON/linears/ChargeLoc.F90
    LSDALTON/linears/SCFLOOP.F90
    LSDALTON/linears/trustradius.F90
    LSDALTON/linears/LSDALTON_RESPONSE.F90
    LSDALTON/linears/LSDALTON_RESPONSE_type.F90
    LSDALTON/linears/arh_debug.F90
    LSDALTON/linears/arh_driver.F90
    LSDALTON/linears/average_util.F90
    LSDALTON/linears/davidson_settings.F90
    LSDALTON/linears/davidson_solver.F90
    LSDALTON/linears/debug.F90
    LSDALTON/linears/dens_subspace.F90
    LSDALTON/linears/densopt.F90
    LSDALTON/linears/diag.F90
    LSDALTON/linears/diis.F90
    LSDALTON/linears/dsm.F90
    LSDALTON/linears/dsm_xterm.F90
    LSDALTON/linears/ecdata_module.F90
    LSDALTON/linears/extra-output.F90
    LSDALTON/linears/fock-eval.F90
    LSDALTON/linears/kurtosis.F90
    LSDALTON/linears/leastchange.F90
    LSDALTON/linears/levelshift.F90
    LSDALTON/linears/minimize.F90
    LSDALTON/linears/optimlocModf.F90
    LSDALTON/linears/orbspread_hess_prec.F90
    LSDALTON/linears/localityMeasure.F90
    LSDALTON/linears/OrbLoc_input.F90
    LSDALTON/linears/localization_util.F90
    LSDALTON/linears/localization_orbspread.F90
    LSDALTON/linears/orbspread_util.F90
    LSDALTON/linears/localization_charge.F90
    LSDALTON/linears/response_driver.F90
    LSDALTON/linears/response_prop.F90
    LSDALTON/linears/scfopt-typedef.F90
    LSDALTON/linears/soeo-redspace.F90
    LSDALTON/linears/soeo-loop.F90
    LSDALTON/linears/soeo-debug.F90
    LSDALTON/linears/soeo-matop.F90
    LSDALTON/linears/soeo-typedef.F90
    LSDALTON/linears/soeo-util.F90
    LSDALTON/linears/starting_guess.F90
    LSDALTON/linears/statistics.F90
    LSDALTON/linears/trilevel.F90
    LSDALTON/linears/plt_driver.F90
    )
set(RSPSOLVER_SOURCES	
    LSDALTON/responsesolver/rsp_cmplx_sym.F90
    LSDALTON/responsesolver/rsp_cmplx_sym_new.F90
    LSDALTON/responsesolver/rsp_complex.F90
    LSDALTON/responsesolver/rsp_solver.F90
    LSDALTON/responsesolver/rsp_std_sym.F90
    )
set(SOLVERUTIL_SOURCES	
    LSDALTON/SolverUtilities/arh_density.F90
    LSDALTON/SolverUtilities/dalton_interface.F90
    LSDALTON/SolverUtilities/dd_utilities.F90
    LSDALTON/SolverUtilities/dd_utilities_unres.F90
    LSDALTON/SolverUtilities/decomp.F90
    LSDALTON/SolverUtilities/lsdalton_modules.F90
    LSDALTON/SolverUtilities/queue-operations.F90
    LSDALTON/SolverUtilities/queue-typedef.F90
    LSDALTON/SolverUtilities/rsp_precond.F90
    LSDALTON/SolverUtilities/rsp_utilities.F90
    )
set(RSP_PROPERTIES_SOURCES	
    LSDALTON/rsp_properties/response_prop_noOpenRSP.F90
    LSDALTON/rsp_properties/molecular_hessian.F90
    LSDALTON/rsp_properties/test_molHessian.F90
    )
set(PBC_FORTRAN_SOURCES
    LSDALTON/pbc2/pbc_compare.F90
    LSDALTON/pbc2/pbc-matop.F90
    LSDALTON/pbc2/pbc-msc.F90
    LSDALTON/pbc2/pbc-krsp-op.F90
    LSDALTON/pbc2/pbc-nfinteractions.F90
    LSDALTON/pbc2/pbc-ffinteractions.F90
    LSDALTON/pbc2/pbc-scf.F90
    LSDALTON/pbc2/pbc-data.F90
    LSDALTON/pbc2/pbc-ffdata.F90
    LSDALTON/pbc2/pbcmain.F90
    LSDALTON/pbc2/pbc-timings.F90
    )
set(FMM_SOURCES
    LSDALTON/mm/mm_T_contractors.F90
    LSDALTON/mm/mm_T_pair_builder.F90
    LSDALTON/mm/mm_T_pair_tests.F90
    LSDALTON/mm/mm_T_worker.F90
    LSDALTON/mm/mm_Vff_driver.F90
    LSDALTON/mm/mm_Vff_processor.F90
    LSDALTON/mm/mm_W_contractors.F90
    LSDALTON/mm/mm_W_pair_builder.F90
    LSDALTON/mm/mm_W_worker.F90
    LSDALTON/mm/mm_box_builder.F90
    LSDALTON/mm/mm_box_processor.F90
    LSDALTON/mm/mm_buffers.F90
    LSDALTON/mm/mm_densfit.F90
    LSDALTON/mm/mm_driver.F90
    LSDALTON/mm/mm_global_paras.F90
    LSDALTON/mm/mm_grid_searcher.F90
    LSDALTON/mm/mm_input.F90
    LSDALTON/mm/mm_interface.F90
    LSDALTON/mm/mm_map_builder.F90
    LSDALTON/mm/mm_memory.F90
    LSDALTON/mm/mm_memory_manager.F90
    LSDALTON/mm/mm_overlap_tree.F90
    LSDALTON/mm/mm_sort.F90
    LSDALTON/mm/mm_stats.F90
    LSDALTON/mm/mm_tree_T_buffer.F90
    LSDALTON/mm/cbifmm.F90
    )
set(INTERESTLIB_SOURCES
    LSDALTON/interest/src/module_interest.f90
    LSDALTON/interest/src/module_interest_hrr.f90
    LSDALTON/interest/src/module_interest_osr.f90
    )
set(ICHORINT_SOURCES
    LSDALTON/IchorIntegrals/MainIchorInterface.F90
    LSDALTON/IchorIntegrals/IchorPresicion.F90
    LSDALTON/IchorIntegrals/IchorParameters.F90
    LSDALTON/IchorIntegrals/IchorInputInfo.F90
    LSDALTON/IchorIntegrals/IchorCommon.F90
    LSDALTON/IchorIntegrals/IchorMem.F90
    LSDALTON/IchorIntegrals/IchorGammaTabulation.F90
    LSDALTON/IchorIntegrals/IchorGaussianGeminal.F90
    LSDALTON/IchorIntegrals/IchorBatchTools.F90
    LSDALTON/IchorIntegrals/IchorSaveGab.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_general.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_Gen.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_SegQ.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_SegP.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_Seg.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_Seg1Prim.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_general.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_Gen.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_SegQ.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_SegP.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_Seg.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_Seg1Prim.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_Gen2.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_SegQ2.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_SegP2.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_Seg2.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_Seg1Prim2.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_Gen2.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_SegQ2.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_SegP2.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_Seg2.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_Seg1Prim2.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_GenSize.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_SegQSize.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_SegPSize.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_SegSize.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_Seg1PrimSize.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_GenSize.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_SegQSize.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_SegPSize.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_SegSize.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_Seg1PrimSize.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_McM_seg_seg_SSSS.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_McM_WTUV_general.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_McM_Ecoeff_specR.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_McM_Ecoeff_specL.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_McM_Ecoeff_specL2.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_McM_Ecoeff_specL3.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_McM_Ecoeff_specL4.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_McM_Ecoeff_general.F90
    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_McM_general.F90
    LSDALTON/IchorIntegrals/IchorEriTools.F90
    LSDALTON/IchorIntegrals/IchorEriDist.F90
    LSDALTON/IchorIntegrals/IchorEriLink.F90
    LSDALTON/IchorIntegrals/IchorEri.F90
    LSDALTON/IchorIntegrals/IchorEri_GabIntegral_McM_general.F90
    LSDALTON/IchorIntegrals/IchorEri_GabIntegral_OBS_Gen.F90
    LSDALTON/IchorIntegrals/IchorEri_GabIntegral_OBS_Seg.F90
    LSDALTON/IchorIntegrals/IchorEri_GabIntegral_OBS_general.F90
    LSDALTON/IchorIntegrals/IchorEri_GabIntegral_Prim.F90
    LSDALTON/IchorIntegrals/IchorGabLoop.F90
    LSDALTON/IchorIntegrals/IchorGab.F90
    LSDALTON/IchorIntegrals/AGC_CPU_BuildRJ000Gen.F90
    LSDALTON/IchorIntegrals/AGC_CPU_BuildRJ000Seg1Prim.F90
    LSDALTON/IchorIntegrals/AGC_GPU_BuildRJ000Gen.F90
    LSDALTON/IchorIntegrals/AGC_GPU_BuildRJ000Seg1Prim.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPAGen.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPBGen.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPCGen.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPDGen.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPASegQ.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPBSegQ.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPCSegQ.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPDSegQ.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPASegP.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPBSegP.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPCSegP.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPDSegP.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPASeg.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPBSeg.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPCSeg.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPDSeg.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPASeg1Prim.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPBSeg1Prim.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPCSeg1Prim.F90
    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPDSeg1Prim.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPAGen.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPBGen.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPCGen.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPDGen.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPASegQ.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPBSegQ.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPCSegQ.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPDSegQ.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPASegP.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPBSegP.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPCSegP.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPDSegP.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPASeg.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPBSeg.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPCSeg.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPDSeg.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPASeg1Prim.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPBSeg1Prim.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPCSeg1Prim.F90
    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPDSeg1Prim.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceParam.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCGen1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCGen2.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDGen1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDGen2.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoAGen.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoAGen.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoCGen1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoDGen1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoBGen.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoBGen.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSegQ1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSegQ2.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDSegQ1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDSegQ2.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoASegQ.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoASegQ.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoCSegQ1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoDSegQ1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoBSegQ.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoBSegQ.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSegP1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSegP2.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDSegP1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDSegP2.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoASegP.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoASegP.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoCSegP1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoDSegP1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoBSegP.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoBSegP.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSeg1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSeg2.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDSeg1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDSeg2.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoASeg.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoASeg.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoCSeg1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoDSeg1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoBSeg.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoBSeg.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSeg1Prim1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSeg1Prim2.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDSeg1Prim1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDSeg1Prim2.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoASeg1Prim.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoASeg1Prim.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoCSeg1Prim1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoDSeg1Prim1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoBSeg1Prim.F90
    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoBSeg1Prim.F90    
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceParam.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCGen1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCGen2.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDGen1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDGen2.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoAGen.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoAGen.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoCGen1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoDGen1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoBGen.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoBGen.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSegQ1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSegQ2.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDSegQ1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDSegQ2.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoASegQ.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoASegQ.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoCSegQ1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoDSegQ1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoBSegQ.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoBSegQ.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSegP1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSegP2.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDSegP1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDSegP2.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoASegP.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoASegP.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoCSegP1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoDSegP1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoBSegP.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoBSegP.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSeg1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSeg2.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDSeg1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDSeg2.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoASeg.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoASeg.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoCSeg1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoDSeg1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoBSeg.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoBSeg.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSeg1Prim1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSeg1Prim2.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDSeg1Prim1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDSeg1Prim2.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoASeg1Prim.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoASeg1Prim.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoCSeg1Prim1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoDSeg1Prim1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoBSeg1Prim.F90
    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoBSeg1Prim.F90
    LSDALTON/IchorIntegrals/AGC_GPU_PrimContractSegQ.F90
    LSDALTON/IchorIntegrals/AGC_GPU_PrimContractSegP.F90
    LSDALTON/IchorIntegrals/AGC_GPU_PrimContractGen.F90
    LSDALTON/IchorIntegrals/AGC_CPU_PrimContractSegQ.F90
    LSDALTON/IchorIntegrals/AGC_CPU_PrimContractSegP.F90
    LSDALTON/IchorIntegrals/AGC_CPU_PrimContractGen.F90
    LSDALTON/IchorIntegrals/AGC_CPU_HorizontalRecurrencePAtoB.F90
    LSDALTON/IchorIntegrals/AGC_CPU_HorizontalRecurrencePBtoA.F90
    LSDALTON/IchorIntegrals/AGC_CPU_HorizontalRecurrenceQCtoD.F90
    LSDALTON/IchorIntegrals/AGC_CPU_HorizontalRecurrenceQDtoC.F90
    LSDALTON/IchorIntegrals/AGC_GPU_HorizontalRecurrencePAtoB.F90
    LSDALTON/IchorIntegrals/AGC_GPU_HorizontalRecurrencePBtoA.F90
    LSDALTON/IchorIntegrals/AGC_GPU_HorizontalRecurrenceQCtoD.F90
    LSDALTON/IchorIntegrals/AGC_GPU_HorizontalRecurrenceQDtoC.F90
    LSDALTON/IchorIntegrals/AGC_CPU_SphContractOBS1.F90
    LSDALTON/IchorIntegrals/AGC_CPU_SphContractOBS2.F90
    LSDALTON/IchorIntegrals/AGC_GPU_SphContractOBS1.F90
    LSDALTON/IchorIntegrals/AGC_GPU_SphContractOBS2.F90
    )
#    LSDALTON/IchorIntegrals/AGC_CPU_BuildRJ000Gen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_BuildRJ000Seg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPAGen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPBGen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPCGen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPDGen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPASegQ.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPBSegQ.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPCSegQ.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPDSegQ.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPASegP.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPBSegP.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPCSegP.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPDSegP.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPASeg.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPBSeg.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPCSeg.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPDSeg.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPASeg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPBSeg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPCSeg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_VerticalRecurrenceQPDSeg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_BuildRJ000Gen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_BuildRJ000Seg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPAGen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPBGen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPCGen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPDGen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPASegQ.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPBSegQ.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPCSegQ.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPDSegQ.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPASegP.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPBSegP.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPCSegP.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPDSegP.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPASeg.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPBSeg.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPCSeg.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPDSeg.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPASeg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPBSeg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPCSeg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_VerticalRecurrenceQPDSeg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCGen1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCGen2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDGen1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDGen2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoAGen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoAGen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoCGen1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoDGen1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoBGen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoBGen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSegQ1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSegQ2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSegQ3.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDSegQ1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDSegQ2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoASegQ.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoASegQ.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoCSegQ1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoDSegQ1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoBSegQ.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoBSegQ.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSegP1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSegP2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSegP3.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDSegP1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDSegP2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoASegP.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoASegP.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoCSegP1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoDSegP1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoBSegP.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoBSegP.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSeg1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSeg2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSeg3.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDSeg1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDSeg2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoASeg.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoASeg.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoCSeg1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoDSeg1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoBSeg.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoBSeg.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSeg1Prim1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoCSeg1Prim2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDSeg1Prim1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceAtoDSeg1Prim2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoASeg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoASeg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoCSeg1Prim1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceBtoDSeg1Prim1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceDtoBSeg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_TransferRecurrenceCtoBSeg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCGen1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCGen2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDGen1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDGen2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoAGen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoAGen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoCGen1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoDGen1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoBGen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoBGen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSegQ1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSegQ2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSegQ3.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDSegQ1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDSegQ2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoASegQ.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoASegQ.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoCSegQ1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoDSegQ1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoBSegQ.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoBSegQ.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSegP1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSegP2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSegP3.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDSegP1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDSegP2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoASegP.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoASegP.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoCSegP1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoDSegP1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoBSegP.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoBSegP.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSeg1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSeg2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSeg3.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDSeg1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDSeg2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoASeg.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoASeg.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoCSeg1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoDSeg1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoBSeg.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoBSeg.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSeg1Prim1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoCSeg1Prim2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDSeg1Prim1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceAtoDSeg1Prim2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoASeg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoASeg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoCSeg1Prim1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceBtoDSeg1Prim1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceDtoBSeg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_TransferRecurrenceCtoBSeg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_HorizontalRecurrencePAtoB.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_HorizontalRecurrencePBtoA.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_HorizontalRecurrencePAtoB.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_HorizontalRecurrencePBtoA.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_HorizontalRecurrenceQCtoD.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_HorizontalRecurrenceQDtoC.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_HorizontalRecurrenceQCtoD.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_HorizontalRecurrenceQDtoC.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_SphContractOBS1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_SphContractOBS2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_SphContractOBS1.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_SphContractOBS2.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_PrimContractSegQ.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_PrimContractSegP.SP.F90
#    LSDALTON/IchorIntegrals/AGC_GPU_PrimContractGen.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_PrimContractSegQ.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_PrimContractSegP.SP.F90
#    LSDALTON/IchorIntegrals/AGC_CPU_PrimContractGen.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_McM_general.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_McM_Ecoeff_general.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_McM_Ecoeff_specR.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_McM_Ecoeff_specL.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_McM_Ecoeff_specL2.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_McM_WTUV_general.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_Gen.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_SegQ.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_SegP.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_Seg.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_Seg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_Gen.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_SegQ.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_SegP.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_Seg.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_Seg1Prim.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_Gen2.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_SegQ2.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_SegP2.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_Seg2.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_Seg1Prim2.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_Gen2.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_SegQ2.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_SegP2.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_Seg2.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_Seg1Prim2.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_GenSize.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_SegQSize.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_SegPSize.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_SegSize.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_CPU_OBS_Seg1PrimSize.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_GenSize.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_SegQSize.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_SegPSize.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_SegSize.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_GPU_OBS_Seg1PrimSize.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_GabIntegral_OBS_Gen.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_GabIntegral_OBS_Seg.SP.F90
#    LSDALTON/IchorIntegrals/IchorEri_CoulombIntegral_McM_SSSS.F90
set(LSINT_SOURCES
    LSDALTON/LSint/IchorInterface.F90
    LSDALTON/LSint/dft_gridLL.F90
    LSDALTON/LSint/BuildBasis.F90
    LSDALTON/LSint/BuildMolFile.F90
    LSDALTON/LSint/gridgeneration.F90
    LSDALTON/LSint/gridgeneration_boxify.F90
    LSDALTON/LSint/THC_grid.F90
    LSDALTON/LSint/THC_util.F90
    LSDALTON/LSint/II_XC_interface.F90
    LSDALTON/LSint/II_absval_int.F90
    LSDALTON/LSint/II_dft_int.F90
    LSDALTON/LSint/II_dft_ksm.F90
    LSDALTON/LSint/II_dft_ksm_worker.F90
    LSDALTON/LSint/Integral_interface.F90
    LSDALTON/LSint/Integral_interfaceDEC.F90
    LSDALTON/LSint/Integral_interfaceDF.F90
    LSDALTON/LSint/Integral_interface_test.F90
    LSDALTON/LSint/MBIE.F90
    LSDALTON/LSint/Molecule.F90
    LSDALTON/LSint/ODbatches.F90
    LSDALTON/LSint/ReadMolFile.F90
    LSDALTON/LSint/ThermiteDistribute.F90
    LSDALTON/LSint/ThermiteDistribute2.F90
    LSDALTON/LSint/ThermiteDistributeGen.F90
    LSDALTON/LSint/ThermiteDistributeK.F90
    LSDALTON/LSint/ThermiteDistributeK2.F90
    LSDALTON/LSint/ThermiteDistributeDEC.F90
    LSDALTON/LSint/ThermiteDriver.F90
    LSDALTON/LSint/ThermiteIntegrals.F90
    LSDALTON/LSint/ThermiteIntTransform.F90
    LSDALTON/LSint/ThermiteOverlapDistribution.F90
    LSDALTON/LSint/ThermiteProp.F90
    LSDALTON/LSint/ThermiteMem.F90
    LSDALTON/LSint/dalton-input.F90
    LSDALTON/LSint/linsolv_df.F90
    LSDALTON/LSint/ls_IntegralInterface.F90
    LSDALTON/LSint/pari.F90
    LSDALTON/LSint/lsmpi.F90
    LSDALTON/LSint/HODItest.F90
    LSDALTON/LSint/II_dft_dftd.F90
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
    LSDALTON/lsutil/ls_precision.F90
    LSDALTON/lsutil/lsmpi_mod.F90
    )
set(LSUTIL_MATRIXM_SOURCES
    LSDALTON/lsutil/matrix_module.F90
    )
set(LSUTIL_COMMON_C_SOURCES
    LSDALTON/lsutil/myPAPI_set_inherit.c
    LSDALTON/lsutil/crayio.c
    )
set(LSUTIL_COMMON_SOURCES
    LSDALTON/lsutil/crayio_util.F90
    LSDALTON/lsutil/rsp-typedef.F90
    LSDALTON/lsutil/tensor_type_def.F90
    LSDALTON/lsutil/response_prop_type.F90
    LSDALTON/lsutil/ls_IOType.F90
    LSDALTON/lsutil/AOType.F90
    LSDALTON/lsutil/dec_typedef.F90
    LSDALTON/lsutil/MoleculeType.F90
    LSDALTON/lsutil/ODType.F90
    LSDALTON/lsutil/lstensorType.F90
    LSDALTON/lsutil/dftType.F90
    LSDALTON/lsutil/BasisinfoType.F90
    LSDALTON/lsutil/IntegralOutputType.F90 
    LSDALTON/lsutil/f12.F90
    LSDALTON/lsutil/IntegralType.F90
    LSDALTON/lsutil/TYPE-DEF.F90
    LSDALTON/lsutil/background_buffer.F90
    LSDALTON/lsutil/memory.F90
    LSDALTON/lsutil/MemoryLeakTool.F90
    LSDALTON/lsutil/gridgeneration_memory.F90
    LSDALTON/lsutil/Time.F90
    LSDALTON/lsutil/common_utilities.F90
    LSDALTON/lsutil/ls_parameters.F90
    LSDALTON/lsutil/par_mod.F90
    LSDALTON/lsutil/lsmpiType.F90
    LSDALTON/lsutil/TestMPI.F90
    LSDALTON/lsutil/LSmatrixType.F90
    LSDALTON/lsutil/lsmatop_dense.F90
    LSDALTON/lsutil/file-operations.F90
    LSDALTON/lsutil/grid_utilities.F90
    LSDALTON/lsutil/mat3dim.F90
    LSDALTON/lsutil/papi.F90
    LSDALTON/lsutil/ls_math.F90
    LSDALTON/lsutil/SphCartMatrices.F90
    LSDALTON/lsutil/OverlapDistributionType.F90
    LSDALTON/lsutil/pbc_lattice_type.F90
    )
set(LSUTIL_TENSOR_SOURCES
    LSDALTON/lsutil/dec_workarounds.F90
    LSDALTON/lsutil/tensor_interface.F90
    LSDALTON/lsutil/lspdm_tensor_operations.F90
    LSDALTON/lsutil/tensor_algebra_dil.F90
    LSDALTON/lsutil/tensor_basic.F90
    LSDALTON/lsutil/lspdm_basic.F90
    )
set(LSUTIL_MATRIXO_SOURCES
    LSDALTON/lsutil/matop_csr.F90
    LSDALTON/lsutil/matop_dense.F90
    LSDALTON/lsutil/matop_dense_unrest.F90
    LSDALTON/lsutil/matop_scalapack.F90
    LSDALTON/lsutil/matop_pdm.F90
    )
set(LSUTIL_MATRIXO_C_SOURCES
    LSDALTON/lsutil/matop_csr_aux.c
    LSDALTON/lsutil/myPAPI_set_inherit.c
    )
set(LSUTIL_MATRIXU_SOURCES
    LSDALTON/lsutil/mat-operations-aux.F90
    LSDALTON/lsutil/mat-operations-essential.F90
    LSDALTON/lsutil/matrix_utilities.F90
    )
set(LSUTIL_TYPE_SOURCES
    LSDALTON/lsutil/AO_operations.F90
    LSDALTON/lsutil/Molecule_operations.F90
    LSDALTON/lsutil/OD_operations.F90
    LSDALTON/lsutil/lstensor_Mem.F90
    LSDALTON/lsutil/lstensor_operations.F90
    LSDALTON/lsutil/Integral_input.F90
    LSDALTON/lsutil/ls_IO.F90
    LSDALTON/lsutil/ls_screen.F90
    LSDALTON/lsutil/dft_operations.F90
    LSDALTON/lsutil/dftMem.F90
    LSDALTON/lsutil/Basisinfo_operations.F90
    LSDALTON/lsutil/IntegralOutput_operations.F90 
    LSDALTON/lsutil/TYPE-OP.F90
    LSDALTON/lsutil/GCtrans.F90
    LSDALTON/lsutil/Build_AOBATCH.F90
    )
set(LSUTILLIB_SOURCES
    LSDALTON/lsutil/lowdin.F90
    LSDALTON/lsutil/AtomicSparse.F90
    LSDALTON/lsutil/fundamental.F90
    LSDALTON/lsutil/ks-settings.F90
    LSDALTON/lsutil/lsmpi-operations.F90
    LSDALTON/lsutil/lsutilities.F90
    LSDALTON/lsutil/pseudoinverse.F90
    LSDALTON/lsutil/LSoptType.F90
    LSDALTON/lsutil/ddynType.F90
    LSDALTON/lsutil/ProfileType.F90
    LSDALTON/lsutil/pbc_lattice_vectors.F90
    LSDALTON/lsutil/lspdm_slave.F90
    )
set(LSLIB_SOURCES
    LSDALTON/lsdaltonsrc/LSlib.F90
    LSDALTON/lsdaltonsrc/LSlibState.F90
    LSDALTON/lsdaltonsrc/LSlib_tester_module.F90
    LSDALTON/lsdaltonsrc/LSlib_tester.F90
    )
set(LSDALTON_FIXED_FORTRAN_SOURCES
    LSDALTON/pdpack/jacobi.F
    LSDALTON/pdpack/eispack.F
    LSDALTON/pdpack/linpack.F
    LSDALTON/pdpack/cholesky.F
    LSDALTON/pdpack/invroutines.F
    )
set(LSDALTON_OWN_BLAS_SOURCES
    LSDALTON/pdpack/blas_wrapper.F
    LSDALTON/pdpack/gp_dblas1.F
    LSDALTON/pdpack/gp_dblas2.F
    LSDALTON/pdpack/gp_dblas3.F
    LSDALTON/pdpack/gp_sblas.F
    LSDALTON/pdpack/gp_qblas.F
    LSDALTON/pdpack/gp_zblas.F
    )
set(LSDALTON_OWN_LAPACK_SOURCES
    LSDALTON/pdpack/gp_dlapack.F
    LSDALTON/pdpack/gp_zlapack.F
    )
# collect all free fortran sources
set(LSDALTON_FREE_FORTRAN_SOURCES
    ${CUDA_GPU_INTERFACE_SOURCES}
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
    ${ICHORINT_SOURCES}
    ${LSUTIL_PRECISION_SOURCES}
    ${LSUTIL_MATRIXM_SOURCES}
    ${LSUTIL_COMMON_SOURCES}
    ${LSUTIL_MATRIXO_SOURCES}
    ${LSUTIL_MATRIXO_C_SOURCES}
    ${LSUTIL_MATRIXU_SOURCES}
    ${LSUTIL_TYPE_SOURCES}
    ${LSUTILLIB_SOURCES}
    ${LSLIB_SOURCES}
    ${ICHORLIB_SOURCES}
 )
if(ENABLE_PCMSOLVER)
    set(LSDALTON_FREE_FORTRAN_SOURCES
        ${LSDALTON_FREE_FORTRAN_SOURCES} 
        ${LSDALTON_PCM_SOURCES})
endif()
