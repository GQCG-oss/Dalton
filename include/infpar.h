#if defined(__CVERSION__)
#define MAXNOD 128
#define MAXCL2 10000
#define MAXTSK (MXSHEL*(MXSHEL+1)/2)
#define NPARI  ((MAXNOD + 1) + 7)
extern struct common_infpar {
    int nodtot, nodeid[MAXNOD+1], ncode,  iprpar, mtottk, ntask,  
	nfmat,  ndegdi, master, mynum,  mytid,  timing,   slave,
	debug, pario;
    char nodnam[MAXNOD][20], myname[20];
} infpar_;
#else
C File: infpar.h
C
C     Parameters NPARI must be updated after changes (for parallelization)
C
C     NOTE: Integers  (MASTER,...)
C           Logicals  (TIMING,...)
C           Character (NODNAM,...) should NOT be sent to slaves
C
      PARAMETER (MAXNOD=128, MAXCL2=10000, MAXTSK=MXSHEL*(MXSHEL+1)/2)
      PARAMETER (NPARI = (MAXNOD + 1) + 7)
      LOGICAL        SLAVE, TIMING, DEBUG, PARIO
      CHARACTER*20   NODNAM, MYNAME
      DIMENSION NODEID(0:MAXNOD), NODNAM(0:MAXNOD)
      COMMON /INFPAR/ NODTOT, NODEID, NCODE, IPRPAR, MTOTTK, NTASK,  
     &                NFMAT,  NDEGDI, 
     &                MASTER, MYNUM,  MYTID,
     &                SLAVE,  TIMING, DEBUG, PARIO,
     &                NODNAM, MYNAME
C -- end of infpar.h --
#endif
