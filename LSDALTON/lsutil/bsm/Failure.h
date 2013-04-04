#ifndef FAILURE
#define FAILURE
#include <string>
namespace bsm
{
  class Failure
    {
      std::string message;
    public:
      Failure()
	:message("No failure information available")
	{}
      explicit Failure (const std::string& msg);
      virtual ~Failure();
      virtual const std::string what() const;
    };
}
#endif
