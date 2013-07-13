#ifndef DIMENSIONF
#define DIMENSIONF
#include "Failure.h"
namespace bsm
{
  class DimensionF:public Failure
    {
    public:
      DimensionF()
	:Failure()
	{}
      explicit DimensionF(const std::string& msg)
	:Failure(msg)
	{}
    };
}

#endif
