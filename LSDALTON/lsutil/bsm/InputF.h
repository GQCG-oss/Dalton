#ifndef INPUTF
#define INPUTF
#include "Failure.h"
namespace bsm
{
  class InputF:public Failure
    {
    public:
      InputF()
	:Failure()
	{}
      explicit InputF(const std::string& msg)
	:Failure(msg)
	{}
    };
}

#endif
