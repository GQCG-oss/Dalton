C added by Bin Gao, Dec. 5, 2010
C...  maximum number of AO-blocks
C...  here, AO-blocks means the orbitals with the same exponents
C...  (contraction coefficients may be different)
      integer MAX_NBLOCK
      parameter (MAX_NBLOCK=200)
C...  number of AO-blocks
      integer NUM_BLOCK
C...  start index of shell in AO-blocks
      integer IDX_SHELL(MAX_NBLOCK)
C...  symmetry integral pointers for AO-blocks
      integer INDFA(8,MXAQN,MXCONT), INDFB(8,MXAQN,MXCONT),
     &        ISOFRA(8), ISOFRB(8)
C...  common block used by gen1int_wrapper.F
      common /GENINT_INFO/ NUM_BLOCK, IDX_SHELL, INDFA, INDFB,
     &                     ISOFRA, ISOFRB
