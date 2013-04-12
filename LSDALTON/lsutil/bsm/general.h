#ifndef GENERAL
#define GENERAL
#include "Failure.h"
#include "IndexF.h"
#include "InputF.h"
#include "DimensionF.h"
#include "IOF.h"

#if defined(VAR_INT64)
/* #include <stdint.h> */
typedef int64_t integer;
/* typedef long integer */
#else
typedef int integer;
#endif

namespace bsm
{

  typedef double real;
  //typedef float real;
  const real ZERO=0;
  const real ONE=1;
  const integer PREC=10;
}

#endif
