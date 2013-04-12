#ifndef IOF_
#define IOF_
#include "Failure.h"
namespace bsm
{
  class IOF:public Failure
    {
    public:
      IOF()
	:Failure()
	{}
      explicit IOF(const std::string& msg)
	:Failure(msg)
	{}
    };
}

#endif
