#define REALK 8
#define _REALK _8
#define MPI_REALK MPI_DOUBLE_PRECISION

#if defined (VAR_INT64)
#define MPI_INTEGERK MPI_INTEGER8
#define MPI_LOGICALK MPI_INTEGER8
#else
#define MPI_INTEGERK MPI_INTEGER4
#define MPI_LOGICALK MPI_INTEGER4
#endif

#define STDOUT 6

#define MAX_LEN_STR 80

#if defined(VAR_MPI)
#define MANAGER 0
#define REQUEST_WORK 0
#define NO_MORE_WORK 0
#endif
