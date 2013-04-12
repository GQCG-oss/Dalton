#ifndef INDEXF
#define INDEXF
#include "Failure.h"
namespace bsm
{
  class IndexF:public Failure
    {
    public:
      IndexF()
	:Failure()
	{}
      explicit IndexF(const std::string& msg)
	:Failure(msg)
	{}
    };
}

#endif
