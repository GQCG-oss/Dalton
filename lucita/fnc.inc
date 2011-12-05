*----------------------------------------------------------------------*
* Translation functions for AIX I/O.                                   *
*----------------------------------------------------------------------*
*     ExtAdr(t)=(t+min_Block_Length-1)/min_Block_Length
      ExtAdr(t)=1+(t-1)/min_Block_Length
      IntAdr(t)=min_Block_Length*t
      BufLen(t)=ItoB*t
