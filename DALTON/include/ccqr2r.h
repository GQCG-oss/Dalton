      INTEGER MXQR2O, MXQR2ST
      PARAMETER ( MXQR2O = 60,MXQR2ST = 1000)

      INTEGER NQR2OP, NXQR2ST, NSEQR2

      INTEGER IAQR2OP, IBQR2OP, IQR2STI, IQR2STF, ISEQR2S

      LOGICAL QR22N1,SELQR2,XOSCST,XVELST

      COMMON /INQR2SD/ IQR2STI(MXQR2ST),IQR2STF(MXQR2ST),
     *                 ISEQR2(MXQR2ST,4),
     *                 IAQR2OP(MXQR2O),IBQR2OP(MXQR2O),
     *                 NQR2OP,NSEQR2,NXQR2ST,QR22N1,SELQR2,
     *                 XOSCST,XVELST


