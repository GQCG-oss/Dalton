#define REALK 8
#define _REALK _8
#define MPI_REALK MPI_DOUBLE_PRECISION

#define STDOUT 6

#define MAX_LEN_STR 80

#if defined(PROG_DIRAC)
#define NUM_COMPONENTS 2
#define LARGE_COMP 1
#define SMALL_COMP 2
#else
#define NUM_COMPONENTS 1
#define LARGE_COMP 1
#endif

#if defined(VAR_MPI)
#define MANAGER 0
#define REQUEST_WORK 0
#define NO_MORE_WORK 0
#endif
