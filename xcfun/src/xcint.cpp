#include <cstdio>
#include <cstdlib>
#include "xcint.h"

functional_data xcint_funs[XC_NR_FUNCTIONALS]; 
parameter_data xcint_params[XC_NR_PARAMETERS_AND_FUNCTIONALS];

template<int FUN>
void xcint_functional_setup_helper()
{
  fundat_db<FUN>::d.symbol = fundat_db<FUN>::symbol;
  fundat_db<FUN>::d.id = (enum xc_functional_id)FUN;
  xcint_funs[FUN] = fundat_db<FUN>::d;
  xcint_functional_setup_helper<FUN+1>();
}

template<> void xcint_functional_setup_helper<XC_NR_FUNCTIONALS>() {  }

template<int P>
void xcint_parameter_setup_helper()
{
  pardat_db<P>::d.symbol = pardat_db<P>::symbol;
  xcint_params[P] = pardat_db<P>::d;
  xcint_parameter_setup_helper<P+1>();
}

template<> 
void xcint_parameter_setup_helper<XC_NR_PARAMETERS_AND_FUNCTIONALS>() {}

vars_data xcint_vars[XC_NR_VARS] =
  {
    {"XC_A", 1, XC_DENSITY},
    {"XC_N", 1, XC_DENSITY},
    {"XC_A_B", 2, XC_DENSITY},
    {"XC_N_S", 2, XC_DENSITY},
    {"XC_A_GAA", 2, XC_DENSITY | XC_GRADIENT},
    {"XC_N_GNN", 2, XC_DENSITY | XC_GRADIENT},
    {"XC_A_B_GAA_GAB_GBB", 5, XC_DENSITY | XC_GRADIENT},
    {"XC_N_S_GNN_GNS_GSS", 5, XC_DENSITY | XC_GRADIENT},
    {"XC_A_B_GAA_GAB_GBB_LAPA_LAPB",7, XC_DENSITY | XC_GRADIENT | XC_LAPLACIAN},
    {"XC_A_B_GAA_GAB_GBB_TAUA_TAUB",7, XC_DENSITY | XC_GRADIENT | XC_KINETIC},
    {"XC_N_S_GNN_GNS_GSS_LAPN_LAPS",7, XC_DENSITY | XC_GRADIENT | XC_LAPLACIAN},
    {"XC_N_S_GNN_GNS_GSS_TAUN_TAUS",7, XC_DENSITY | XC_GRADIENT | XC_KINETIC},
    {"XC_A_2ND_TAYLOR", 10, XC_DENSITY | XC_GRADIENT | XC_LAPLACIAN},
    {"XC_A_B_2ND_TAYLOR", 20, XC_DENSITY | XC_GRADIENT | XC_LAPLACIAN},
    {"XC_N_2ND_TAYLOR", 10, XC_DENSITY | XC_GRADIENT | XC_LAPLACIAN},
    {"XC_N_S_2ND_TAYLOR", 20, XC_DENSITY | XC_GRADIENT | XC_LAPLACIAN},
  };


void xcint_assure_setup()
{
  static bool is_setup = false;
  if (!is_setup)
    {
      xcint_functional_setup_helper<0>();      
      xcint_parameter_setup_helper<XC_NR_FUNCTIONALS>();
      is_setup = true;
    }
}

void xcint_die(const char *message, int code)
{
  fprintf(stderr,"XCFun fatal error %i: ",code);
  fprintf(stderr,"%s",message);
  fprintf(stderr,"\n");
  exit(-1);
}
