#include <linesep.h>
C     IRAT  = (real word length) / (integer word length)
C     IRAT2 = (real word length) / (half-integer word length)
C             if available and used, otherwise IRAT2 = IRAT
C     LRAT  = (real word length) / (logical word length)
      INTEGER IRAT, IRAT2, LRAT
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_NEC) || defined (SYS_T90)
#if defined (VAR_STAR4)
      PARAMETER (IRAT = 1, IRAT2 = 2, LRAT = 1)
#else
      PARAMETER (IRAT = 1, IRAT2 = 1, LRAT = 1)
#endif
#else
#if defined (VAR_STAR2)
      PARAMETER (IRAT = 2, IRAT2 = 4, LRAT = 2)
#else
      PARAMETER (IRAT = 2, IRAT2 = 2, LRAT = 2)
#endif
#endif
