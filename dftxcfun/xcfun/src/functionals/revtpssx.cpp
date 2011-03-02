#include "functional.h"
#include "constants.h"
#include "revtpssx_eps.h"

template<class num>
static num revtpssx(const densvars<num> &d)
{
  return  revtpssx_eps::revtpssx_eps(d);
}      

FUNCTIONAL(XC_REVTPSSX) = {
  "Reviewed TPSS exchange functional",
  "Reviewed TPSS exchange functional.\n"
  "J.P. Perdew, A. Ruzsinszky, G. I. Csonka, L. A. Constantin, J. Sun,\n"
  "Workhorse Semilocal Density Functional\n"
  "for Condensed Matter Physics and Quantum Chemistry\n" 
  "Phys. Rev. Lett. 103 (2009) 026403\n"
  "Implemented by Andrea Debnarova\n",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(revtpssx)
};

