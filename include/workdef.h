C
C$Id: workdef.h,v 1.1.1.1 2001-02-08 13:33:24 hjj Exp $
C
#include <linesep.h>
#if defined (VAR_MOTECC)
#if defined (SYS_XA)
      PARAMETER (LWORK = 4 000 000)
#else
      PARAMETER (LWORK = 1 000 000)
#endif
#else
#if !defined (VAR_TESTIBM)
      PARAMETER (LWORK = 4 000 000)
#else
C  880626: LWORK reduced for Georgia 9370 IBM
      PARAMETER (LWORK = 500 000)
#endif
#endif
