      LOGICAL DOHUCKEL
C
C     Worst case scenario for HUCEXC; a user running STO-3G (Hueckel basis
C     half the size of the total basis set size
C
      COMMON /HUCKEL/ HUCCNT, HUCEXC(MXSHEL), DOHUCKEL, NHUCAO(8),
     &                NHUCBA, IHUCPT(MXSHEL)
