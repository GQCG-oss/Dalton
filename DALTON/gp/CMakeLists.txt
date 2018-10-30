add_library(gp "")

target_sources(
  gp
  PRIVATE
    gpdahle.F
    gphjj.F
    gpjeppe.F
    gpjkp.F
    gpjt.F
    gplucia.F
    gplucia_io.F
    gplupar.F
    gpsaue.F
    gpsimen.F
    gptrygve.F
    graphic.F
    intf7790.F
    mempkg.F
    mpimacro.F
    pvmmacro.F
    qpack.F
    vcraypack.F
    one_sided_communication_wrappers.F90
    memory_parallel.F90
  )

if(NOT ENABLE_MPI)
  target_sources(
    gp
    PRIVATE
      mpi_dummyc.c
      mpi_dummy.F
    )
endif()
