!> Module for FLOP counting using PAPI library.
!> extended to include GPU info
!> FYI
!> ***
!Compilation on grendel, currently only installed on node s57n41
!icc -I/home/brano/bin/papi-4.2.1/include myPAPI_set_inherit.c -c
!ifort -openmp -fpp papitest1_inherit.f90 myPAPI_set_inherit.o -I/home/brano/bin/papi-4.2.1/include -Wl,-rpath=/home/brano/bin/papi-4.2.1/lib -L/home/brano/bin/papi-4.2.1/lib -lpapi -mkl=parallel -o papitest_inherit.x

!Compilation on jaguar
!module load PrgEnv-intel
!module load papi
!cc -c myPAPI_set_inherit.c
!ftn -fpp1 papitest_inherit.f90 myPAPI_set_inherit.o -o papitest_inherit.x
module papi_module
use precision
!> EPAPI vent set which is initiated at the very beginning of LSDALTON.
integer,save :: eventset

contains

#ifdef VAR_PAPI

  !> \brief Init papi FLOP counting.
  ! Note: This is done ONCE and then one can count FLOPS for many different processes.
  !> \author Branislav Jansik (modified for LSDALTON by Kasper Kristensen)
  subroutine mypapi_init(es)
    implicit none
#include <f90papi.h>
    integer :: es
    integer :: event, retval

    !init PAPI library
    retval= PAPI_VER_CURRENT
    call PAPIf_library_init(retval)

    !set event to count floating point istructions
    event = PAPI_FP_INS
    call PAPIf_query_event(event, retval)

    !instrument thread
    es = PAPI_NULL

    call PAPIf_create_eventset(es, retval)
    call myPAPIf_set_inherit(es, retval)
    call PAPIf_add_event( es, event, retval )

  end subroutine mypapi_init

  !Wrapper to PAPIf_start
  subroutine mypapi_start(es)  
    integer :: es
    integer :: event, retval
    
    event = PAPI_FP_INS
    call PAPIf_query_event(event, retval)
    !instrument thread
    es = PAPI_NULL
    call PAPIf_create_eventset(es, retval)
    call myPAPIf_set_inherit(es, retval)
    call PAPIf_add_event( es, event, retval )
  end subroutine mypapi_start

  !Wrapper to PAPIf_stop
  subroutine mypapi_stop(es,flops)  
    integer,intent(inout) :: es
    integer(8),intent(inout) :: flops
    !local variables
    integer :: retval    
    flops=0 ! zero flops (this is probably redundant)
    call PAPIf_stop(es,flops,retval)
  end subroutine mypapi_stop


  !> Simple example of FLOP counting using PAPI
  !> This assumes that mypapi_init was called at the very beginning of the lsdalton
  !> subroutine to initiate the eventset stored in the global integer "eventset".
  !> \author Branislav Jansik (modified for LSDALTON by Kasper Kristensen)
  subroutine papi_example(lupriIN)
    implicit none
    integer,optional :: lupriIN
    integer(8) :: flops
    integer :: retval, i,lupri
    integer, parameter :: n=500, n1=50
    real(8)     :: a(n,n),b(n,n),c(n,n)
    integer, external :: omp_get_max_threads
    integer :: eventset2
    lupri = 6
    IF(present(lupriIN)) lupri = lupriIN

    ! initialize stuff to be used in this example (not PAPI related)
    ! **************************************************************
    call RANDOM_SEED()
    call random_number(a)
    call random_number(b)
    a=2.0
    b=3.0
    c=0.0
    i=omp_get_max_threads()
    write (lupri,*) 'Using ', i, ' threads'
    write (lupri,*) 'Matrix(n,n)    n=', n
    write (lupri,*) 'Matrix(n1,n1) n1=', n


    ! start PAPI counters
    ! *******************
    call mypapi_start(eventset2)  

    call PAPIf_start(eventset, retval)

    ! Do some operations - dgemm in this simple example
    ! *************************************************
    call dgemm('n','n',n,n,n,1.0E0_realk,a,n,b,n,0.0E0_realk,c,n)

    ! Stop and read PAPI counters
    ! ****************************
    flops=0 ! zero flops (this is probably redundant)
    call PAPIf_stop(eventset,flops,retval)
    write(lupri,*) 'FLOPS for first  dgemm = ', flops


    ! One can now have as many new Papi sections as requested. E.g. we could do:
    call PAPIf_start(eventset, retval)
    call dgemm('n','n',n1,n1,n1,1.0E0_realk,a,n1,b,n1,0.0E0_realk,c,n1)

    flops=0
    call PAPIf_stop(eventset,flops,retval)
    write(lupri,*) 'FLOPS for second dgemm = ', flops

    flops=0
    call mypapi_stop(eventset2,flops)  
    write(lupri,*) 'FLOPS for both dgemm = ', flops

  end subroutine papi_example
#endif

end module papi_module
