!     FILE: trinfo.h
      PARAMETER ( NSYMP = 36 )
      INTEGER*8      LTEST ! for comparison with LRUPQ, LUSPQ, LTUPQ because
                           ! total number of integrals can be greater than 2G
      COMMON/TRINFO/ ISP,ISQ,ISR,ISS,NBP,NBQ,NBR,NBS,NBPQ,NBRS,IRRST,
     &               NOCP,NOCQ,NOCR,NOCS,NPQ,LADX,LRUPQ,LUSPQ,LTUPQ,
     &               NOP,NOQ,NOR,NOS,LMOP,LMOQ,LMOR,LMOS,
     &               LMOP2,LMOQ2,LMOR2,LMOS2,ITP,ITQ,ITR,ITS,
     &               JMOP2,JMOQ2,JMOR2,JMOS2,LX2X3,
     &               IAD1S,IAD2S,IAD3S,IPQMX1,IPQMX2,IPQMX3
!  -- end of trinfo.h --
