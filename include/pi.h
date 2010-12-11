#if defined (__CVERSION__)
  real pi = 3.14159265358979323846;
#define pi2 pi*pi
  real sqrtpi = 1.77245385090551602730;
  real r2pi52 = 5.91496717279561287782;
#else
C File: pi.h
C
      REAL*8     PI, PI2, SQRTPI, R2PI52
      PARAMETER (PI     = 3.14159 26535 89793 23846 D00, PI2 = PI*PI,
     &           SQRTPI = 1.77245 38509 05516 02730 D00,
     &           R2PI52 = 5.91496 71727 95612 87782 D00)
C     R2PI52 = sqrt(2 * sqrt(PI^5) ) -- used in calc. of 2-el. integrals
C -- end of pi.h --
#endif
