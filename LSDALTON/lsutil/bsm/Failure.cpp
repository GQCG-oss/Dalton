#include "Failure.h"
namespace bsm
{
  Failure::Failure (const std::string& msg)
    :message(msg)
  {}
  Failure::~Failure(){}
  const std::string Failure::what() const
  {return message;}
}
