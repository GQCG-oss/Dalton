C
C$Id: huckel.h,v 1.1.1.1 2001-02-08 13:33:28 hjj Exp $
C
      LOGICAL ADDSTO
C
C     Worst case scenario for HUCEXC; a user running STO-3G (Hueckel basis
C     half the size of the total basis set size
C
      COMMON /HUCKEL/ HUCCNT, HUCEXC(MXSHEL), ADDSTO, NHUCAO(8), 
     &                NHUCBA, IHUCPT(MXSHEL)

