C
C$Id: leinf.h,v 1.1.1.1 2001-02-08 13:33:28 hjj Exp $
C
C     LEINF : space nedded in LELIN = LLEWA + nsim*LLEWB
C             NSIDE=0, symmetric; =1 from left side; =2 from right side
C             NCREF = # of reference vectors to othogonalize against
      COMMON /LEINF / THRLE, MAXLE, IPRLE, LETYPA,LETYPB,LETOT, NSIDE,
     *                NCREF, LULEA, LULEB, LULEC,
     *                LLEWA, LLEWB, KLELIN(20),   NLELIN,LLELIN
      EQUIVALENCE (LETYP2,LETYPA), (LETYP1,LETYPB)
