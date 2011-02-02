#include "functional.h"
#include "constants.h"

template<class num>
static num ENERGY_FUNCTION(XC_OPTX)(const densvars<num> &d)
{
  const parameter a1 = 1.05151, a2 = 1.43169, gamma = 0.006;
  num g_xa2 = gamma*d.gaa*pow(d.a,-8.0/3.0);
  num g_xb2 = gamma*d.gbb*pow(d.b,-8.0/3.0);
  return
    -(d.a_43*(a1*xc_constants::c_slater + a2*pow(g_xa2,2)*pow(1+g_xa2,-2)))
    -(d.b_43*(a1*xc_constants::c_slater + a2*pow(g_xb2,2)*pow(1+g_xb2,-2)));
}


NEW_GGA_FUNCTIONAL(XC_OPTX);
SHORT_DESCRIPTION(XC_OPTX) = "OPTX Handy & Cohen exchange";
LONG_DESCRIPTION(XC_OPTX) =
	     "OPTX Handy & Cohen exchange GGA exchange functional\n"
	     "Implemented by Ulf Ekstrom\n";
NO_TEST(XC_OPTX);

