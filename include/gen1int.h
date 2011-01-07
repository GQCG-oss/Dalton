C Common block of Gen1Int
C added by Bin Gao, Dec. 5, 2010
C...  if running Gen1Int
      logical RUN_GENINT
C...  maximum number of AO-blocks
      integer MAX_NBLOCK
      parameter (MAX_NBLOCK=200)
C...  number of AO-blocks
      integer NUM_BLOCK
C...  start index of shell in AO-blocks
      integer IDX_SHELL(MAX_NBLOCK)
C...
      common /GENINT_INFO/ RUN_GENINT, NUM_BLOCK, IDX_SHELL
