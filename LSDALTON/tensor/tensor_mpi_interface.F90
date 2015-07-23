module tensor_mpi_interface_module
   use tensor_parameters_and_counters
   use tensor_mpi_binding

   !NOTE THAT ONLY MPI_IN_PLACE HAS BEEN IMPLEMENTED FOR THE MPI COLLECTIVES
   public :: tensor_mpi_bcast
   public :: tensor_mpi_reduce
   public :: tensor_mpi_allreduce
   private

   interface tensor_mpi_bcast
      module procedure
   end interface tensor_mpi_bcast

   interface tensor_mpi_reduce
      module procedure
   end interface tensor_mpi_reduce

   interface tensor_mpi_allreduce
      module procedure
   end interface tensor_mpi_allreduce

   contains


end module tensor_mpi_interface_module
