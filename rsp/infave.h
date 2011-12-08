C (notes by HJAaJ Sep 2003)
C infave.h has information for rsp/expone.F:
C       N1AVTOT 1-electron expectation values of LBL1AV(1:N1AVTOT)
C       N2AVTOT 2-electron expectation values of LBL2AV(1:2,1:N2AVTOT)
C
      PARAMETER ( MAX1AV = MAXLBL, MAX2AV = 10 )
      CHARACTER*8 LBL1AV, LBL2AV
      COMMON /INFAV/ N1AVTOT,N2AVTOT
      COMMON /CHRAV/ LBL1AV(MAX1AV), LBL2AV(2,MAX2AV)
