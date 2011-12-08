C File: memcb1.h   (for mempkg.F package)

      COMMON /MEMCB1/ WIDENT_KFREE, LUWMEM, LUEMEM, NWNMEM

Chjaaj sep 2008: new standard (MEMGET2)
C     (1) ident, (2) length, (3) tag (here to catch VEC(0) errors)
#ifdef DBG_LENID
Chjaaj mar 2005: make LENID .gt. standard size for extra debugging
      PARAMETER ( LENID = 100)
#else
      PARAMETER ( LENID = 3)
#endif

      REAL*8        WLREAL, WLZERO
      INTEGER*8      LREAL,  LZERO
      EQUIVALENCE  (WLREAL, LREAL), (WLZERO, LZERO)
      CHARACTER*8   IDENT8

      REAL*8        WMEMID, WMEMCK
      INTEGER*8      MEMID,  MEMCK
      EQUIVALENCE  (WMEMID, MEMID), (WMEMCK, MEMCK)

      SAVE          LZERO,    MEMID
      DATA          LZERO/0/, MEMID /1234567890/
C -- end of memcb1.h
