#ifdef COMMENT
! -- infpar.h --
!     my_MPI_INTEGER is used both in .c and .F routines in MPI calls
!        so we can handle "-i8" compilations on 32-bit machines,
!        using INT_STAR8 /Jan-2007 hjaaj
#endif
#if defined (VAR_INT64)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

#if defined(__CVERSION__)
#define MAXNOD 200
#define MAXCL2 10000
#define MAXTSK (MXSHEL*(MXSHEL+1)/2)
#define NPARI  ((MAXNOD + 1) + 7)
extern struct common_infpar {
#if defined (VAR_INT64)
    long iprpar, ntask, ncode, ndegdi, master, mynum, mytid, nodtot;
    long nodeid[MAXNOD+1], nfmat, mtottk;
    long parher, slave, timing, debug, pario;
#else
    int  iprpar, ntask, ncode, ndegdi, master, mynum, mytid, nodtot;
    int  nodeid[MAXNOD+1], nfmat, mtottk;
    int  parher, slave, timing, debug, pario;
#endif
    char nodnam[MAXNOD][20], myname[20];
} daltoninfpar_;
#else
C File: infpar.h for Dalton; special information for parallel calculations
C
C     Parameters NPARI must be updated after changes (for parallelization)
C
C     NOTE: Integers  (MASTER,...)
C           Logicals  (TIMING,...)
C           Character (NODNAM,...) should NOT be sent to slaves
C
      INTEGER MAXNOD, MAXCL2
      PARAMETER (MAXNOD = 200, MAXCL2 = 10000)
      PARAMETER (MAXTSK=MXSHEL*(MXSHEL+1)/2, NPARI = (MAXNOD + 1) + 7)
      INTEGER IPRPAR, NTASK, NCODE, NDEGDI, MASTER, MYNUM, MYTID, NODTOT
      INTEGER NODEID(0:MAXNOD), NFMAT, MTOTTK
      LOGICAL PARHER, SLAVE, TIMING, DEBUG, PARIO
      CHARACTER*20   NODNAM(0:MAXNOD), MYNAME
      COMMON /DALTONINFPAR/                                              &
     &        IPRPAR, NTASK, NCODE, NDEGDI, MASTER, MYNUM, MYTID, NODTOT &
     &       ,NODEID, NFMAT, MTOTTK                                      &
     &       ,PARHER, SLAVE, TIMING, DEBUG, PARIO                        &
     &       ,NODNAM, MYNAME

#if defined (INT_STAR8)
!     integer array ISTAT contains MPI_SOURCE information.
!     Proper use of ISTAT on 64-bit machines in 
!     combination with INT_STAR8 requires explicit declaration 
!     as INTEGER*4 /March-2007 sk 
      INTEGER*4 ISTAT
#endif

C -- end of infpar.h --
#endif
