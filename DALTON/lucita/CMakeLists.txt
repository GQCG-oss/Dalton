add_library(lucita "")

target_sources(
  lucita
  PRIVATE
    ciinfo.F
    dalton_interface.F
    density.F
    diagonal.F
    diagonalize.F
    fock.F
    io.F
    memory.F
    par_lucia.F
    par_solver.F
    program.F
    property.F
    sigma.F
    strings.F
    transform.F
    utils.F
    vecmat.F
    dalton_mpi.F90
    lucita_cfg.F90
    lucita_mcscf_ci_cfg.F90
    lucita_orbitalspaces.F90
    par_communication_model.F90
    par_io.F90
    par_io_fh.F90
    par_sync_cw.F90
  )

target_link_libraries(
  lucita
  PRIVATE
    gp
  )
