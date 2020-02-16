#ifdef COMMENT
! -- infpar.h --
!     my_MPI_INTEGER is used both in .c and .F routines in MPI calls
!        so we can handle "-i8" compilations on 32-bit machines,
!        using VAR_INT64 /Jan 2007 hjaaj
!
!    MAXNOD is hardcoded max number of nodes (slaves). If you increase it
!    beyond 9999 you must also allow node file names with more than 4 digits
!    in GPOPEN etc. in gp/gptrygve.F /Jan 2017 hjaaj
!
!    infpar.h depends on maxorb.h
!    use #include "maxorb.h" before including infpar. - FBeyer 20140302
!
!    Note: my_MPI_LOGICAL points to MPI_INTEGER4/8 because in Fortran
!          the default size of logical and integer is the same.
#endif
#if defined (VAR_INT64)
#define my_MPI_INTEGER MPI_INTEGER8
#define my_MPI_LOGICAL MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#define my_MPI_LOGICAL MPI_INTEGER4
#endif

#if defined(__CVERSION__)
#define MAXNOD 9999
#define NPARI  7
extern struct common_infpar {
#if defined (VAR_INT64)
    long iprpar, ntask, ndegdi, master, mynum, mytid;
    long nodtot, nfmat, mtottk, parher, debug, pario;
    long timing, slave, rma_model;
#else
    int  iprpar, ntask, ndegdi, master, mynum, mytid;
    int  nodtot, nfmat, mtottk, parher, debug, pario;
    int  rma_model, timing, slave;
#endif
    char nodnam[MAXNOD][20], myname[20];
} daltoninfpar_;
#else
! File: infpar.h for Dalton; special information for parallel calculations
!
!     Parameters NPARI must be updated after changes (for communication to co-workers)
!
!     NOTE: Integers  (IPRPAR,...,MASTER,...,MYTID)
!           Logicals  (TIMING,SLAVE)
!           Character (NODNAM,MYNAME) should NOT be sent to slaves
!     THUS: NPARI is length from NODTOT,...,PARIO,rma_model
!
      INTEGER   MAXNOD, NPARI
      PARAMETER ( MAXNOD = 9999, NPARI = 7 )
      INTEGER IPRPAR, NTASK, NDEGDI, MASTER, MYNUM, MYTID
      INTEGER NODTOT, NFMAT, MTOTTK
      LOGICAL PARHER, PARIO, INFPAR_DEBUG, TIMING, SLAVE, rma_model
      CHARACTER*20   NODNAM(0:MAXNOD), MYNAME
      COMMON /DALTONINFPAR/                                              &
     &        IPRPAR, NTASK, NDEGDI, MASTER, MYNUM, MYTID                &
     &       ,NODTOT, NFMAT, MTOTTK, PARHER, INFPAR_DEBUG, PARIO         &
     &       ,rma_model, TIMING, SLAVE , NODNAM, MYNAME

! -- end of infpar.h --
#endif
