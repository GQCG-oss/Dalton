C
C$Id: nosync.h,v 1.1.1.1 2001-02-08 13:33:27 hjj Exp $
C
#if defined (SYS_ALLIANT)
CVD$ NOSYNC
#endif
#if !defined (SYS_ALLIANT)
Cparallelization note: no synchronization problems
#include <ivdep.h>
#endif
