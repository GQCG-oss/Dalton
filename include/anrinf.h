C
C$Id: anrinf.h,v 1.1.1.1 2001-02-08 13:33:28 hjj Exp $
C
C    ANRINF : space nedded in NRLIN = LNRWA + nsim*LNRWB
C             NSIDE=0, symmetric; =1 from left side; =2 from right side
C             NREFS = # of reference vectors to othogonalize against
      COMMON /ANRINF/ THRNR, MAXNR, IPRNR, NRTYPA,NRTYPB,NRTOT, NSIDE,
     *                NREFS, LNRWA, LNRWB, KNRLIN(20),   NNRLIN,LNRLIN
