C
C$Id: dummy_dec.h,v 1.1 2001-02-12 18:24:07 vebjornb Exp $
C
#include <linesep.h>
C     make DUMMY write protected (on most computers)
#if defined (SYS_CRAY)
      REAL DUMMY
#else
      DOUBLE PRECISION DUMMY
#endif
C
      INTEGER IDUMMY
