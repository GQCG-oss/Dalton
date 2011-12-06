/* mpi_dummyc.c: just a stub routines to make stuff work. We cannot
 * include mpi.h because this file is compiled even in sequential code
 * and the user might not have this file.  Instead, we try to fake
 * the needed types.  Comments and improvement suggestions to
 * pawsa@theochem.kth.se.  Flames to /dev/null. */
typedef int  MPI_Comm;
typedef int  MPI_Datatype;
typedef int  MPI_Op;
typedef struct MPI_Status_ MPI_Status;
#include <stdlib.h>
#define R return 0
int MPI_Barrier(MPI_Comm c)                                             { /* No OP */ R; }
int MPI_Bcast(void* d, int n, MPI_Datatype t, int m, MPI_Comm c)        { /* No OP */ R; }
int MPI_Comm_size(MPI_Comm c, int *sz)                                  { *sz = 1;    R; }
int MPI_Comm_rank(MPI_Comm c, int *r)                                   { *r  = 0;    R; }
int MPI_Send(void* d, int n, MPI_Datatype t, int a, int b, MPI_Comm c)  { /* No OP */ R; }
int MPI_Recv(void* d, int n, MPI_Datatype t, int a, int b, MPI_Comm c,
	     MPI_Status *stat) { /* No OP */ R; }
int MPI_Reduce(void* s, void*d, int n, MPI_Datatype t, MPI_Op o, int r,
	       MPI_Comm c) { abort(); }
