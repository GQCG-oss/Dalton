C
C$Id: indqr.h,v 1.1.1.1 2001-02-08 13:33:26 hjj Exp $
C
      PARAMETER ( MXEXQR = 60 , MXLRQR = 60 )
      CHARACTER*8 QRLBL ,TRLBL
      COMMON /INDQR/ JEXQR(MXEXQR),ISEXQR(MXEXQR),
     *               ILRQR(8),NLRQR(8),ITRQR(8),NTRQR(8),IEXQR(8),
     *               NEXQR(8),ISYMTR(MXLRQR),
     *               TRFREQ(MXLRQR),QRFREQ(MXLRQR),
     *               ISYMQR(MXLRQR),EXCITA(8,MXEXQR,2),
     *               NEXLBL,NEXLB2,NLRLBL,NTRLBL,
     *               TRLBL(MXLRQR),QRLBL(MXLRQR)
