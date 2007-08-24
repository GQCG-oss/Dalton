C File: memcb1.h   (for mempkg.F package)
      COMMON /MEMCB1/ LUWMEM, LUEMEM, NWNMEM
#ifdef VAR_INT64
Chjaaj mar 2005: standard size for Integer*8
C     (1) tag, (2) length, (3) tag
      PARAMETER ( LENID = 3)
#else
Chjaaj mar 2005: standard size for Integer*4
C     (1) tag + length, (2) tag+length
      PARAMETER ( LENID = 2)
#endif
Chjaaj mar 2005: make LENID .gt. standard size for extra debugging
Cdbg  PARAMETER ( LENID = 100)
cdbg  remember to cpp define DBG_LENID if LENID .gt. 2
C
      REAL*8        WMEMID(2), WMEMCK(3)
#ifdef VAR_INT64
      INTEGER*8     MEMID(2),  MEMCK(3)
#else
      INTEGER*4     MEMID(4),  MEMCK(6)
#endif
      EQUIVALENCE  (MEMID, WMEMID), (MEMCK, WMEMCK)
      SAVE          MEMID
      DATA          MEMID(1) /1234567890/
C -- end of memcb1.h
