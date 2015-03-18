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
!> module variable to count the FLOPs done on the GPU
real(realk),save :: FLOPonGPU

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

  !> Simple example of FLOP counting using PAPI
  !> This assumes that mypapi_init was called at the very beginning of the lsdalton
  !> subroutine to initiate the eventset stored in the global integer "eventset".
  !> \author Branislav Jansik (modified for LSDALTON by Kasper Kristensen)
  subroutine papi_example

    implicit none
    integer(8) :: flops
    integer :: retval, i
    integer, parameter :: n=500, n1=50
    real(8)     :: a(n,n),b(n,n),c(n,n)
    integer, external :: omp_get_max_threads


    ! initialize stuff to be used in this example (not PAPI related)
    ! **************************************************************
    call RANDOM_SEED()
    call random_number(a)
    call random_number(b)
    a=2.0
    b=3.0
    c=0.0
    i=omp_get_max_threads()
    write (*,*) 'Using ', i, ' threads'
    write (*,*) 'Matrix ', n


    ! start PAPI counters
    ! *******************
    call PAPIf_start(eventset, retval)

    ! Do some operations - dgemm in this simple example
    ! *************************************************
    call dgemm('n','n',n,n,n,1.0E0_realk,a,n,b,n,0.0E0_realk,c,n)

    ! Stop and read PAPI counters
    ! ****************************
    flops=0 ! zero flops (this is probably redundant)
    call PAPIf_stop(eventset,flops,retval)
    write(*,*) 'FLOPS for first  dgemm = ', flops


    ! One can now have as many new Papi sections as requested. E.g. we could do:
    call PAPIf_start(eventset, retval)
    call dgemm('n','n',n1,n1,n1,1.0E0_realk,a,n1,b,n1,0.0E0_realk,c,n1)
    flops=0
    call PAPIf_stop(eventset,flops,retval)
    write(*,*) 'FLOPS for second dgemm = ', flops


  end subroutine papi_example
#endif

  subroutine AddFLOP_FLOPonGPUaccouting(inputFLOPonGPU)
    implicit none
    real(realk),intent(in) :: inputFLOPonGPU
    FLOPonGPU = FLOPonGPU + inputFLOPonGPU
  end subroutine AddFLOP_FLOPonGPUaccouting

  subroutine addDGEMM_FLOPonGPUaccouting(M,N,K,beta)
    implicit none
    integer,intent(in) :: M,N,K
    real(realk),intent(in) :: beta
    IF(ABS(beta).GT.1.0E-10_realk)THEN
       FLOPonGPU = FLOPonGPU + 2*M*N*K
    ELSE
       FLOPonGPU = FLOPonGPU + M*N*(2*K-1)
    ENDIF
  end subroutine AddDGEMM_FLOPonGPUaccouting
  
  subroutine init_FLOPonGPUaccouting()
    implicit none
    FLOPonGPU = 0.0E0_realk
    
  end subroutine Init_FLOPonGPUaccouting

  subroutine extract_FLOPonGPUaccouting(outFLOPonGPU)
    implicit none
    real(realk),intent(inout) :: outFLOPonGPU
    outFLOPonGPU = FLOPonGPU
  end subroutine Extract_FLOPonGPUaccouting

end module papi_module
