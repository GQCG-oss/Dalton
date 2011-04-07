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
C...  maximum number of centers in the integral
      integer MAX_NCENT
      parameter (MAX_NCENT=4)
C...  indices of atomic centers
      integer idx_cent(MAX_NCENT)
C...  order of derivatives of the corresponding atomic centers
      integer order_cent(MAX_NCENT)
C...  allowed maximum order of geometric derivatives
      integer MAX_ORDER_GEO
      parameter (MAX_ORDER_GEO=30)
C...  allowed maximum order of magnetic derivatives
      integer MAX_ORDER_MAG
      parameter (MAX_ORDER_MAG=10)
C...  common block used by gen1int_wrapper.F
      common /GENINT_INFO/ NUM_BLOCK, IDX_SHELL, idx_cent, order_cent
