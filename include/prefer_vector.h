C
C$Id: prefer_vector.h,v 1.1.1.1 2001-02-08 13:33:27 hjj Exp $
C
#if defined (SYS_CONVEX) || defined (SYS_HAL)
C$DIR PREFER_VECTOR
#else
Cvectorization note --> prefer vectorization over this loop
#endif
