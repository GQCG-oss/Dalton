#include <linesep.h>
C     make DUMMY write protected (on most computers)
#if defined (SYS_CRAY)
      REAL DUMMY
#else
      DOUBLE PRECISION DUMMY
#endif
C
      INTEGER IDUMMY
