C
C$Id: lbmxsq.h,v 1.1.1.1 2001-02-08 13:33:27 hjj Exp $
C
      PARAMETER (LBMXSQ = 4094)
C =4096-2; 4094*12+4 (16 for Cray) will then fit in n*4kb disk buffers
C NOTE: LBMXSQ must be even
