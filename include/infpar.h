C
C     Parameters NPARI must be updated after changes (for parallelization)
C
C     NOTE: Integers  (MASTER,...)
C           Logicals  (TIMING,...)
C           Character (NODNAM,...) should NOT be sent to slaves
C
      PARAMETER (MAXNOD=100, MAXCL2=10000, MAXTSK=MXSHEL*(MXSHEL+1)/2)
      PARAMETER (NPARI = (MAXNOD + 1) + 7)
      LOGICAL SLAVE, TIMING, DEBUG, PARIO
      CHARACTER*20 NODNAM, MYNAME
      DIMENSION NODEID(0:MAXNOD), NODNAM(0:MAXNOD)
      COMMON /INFPAR/ NODTOT, NODEID, NCODE, IPRPAR, MTOTTK, NTASK,  
     &                NFMAT,  NDEGDI, 
     &                MASTER, MYNUM,  MYTID,
     &                TIMING, SLAVE,  DEBUG,
     &                NODNAM, MYNAME
      COMMON /PARIOC/ PARIO
