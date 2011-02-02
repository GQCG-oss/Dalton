
#ifndef XC_NO_REGULARIZATION
#error Implement regularization properly, what about the non-constant terms when setting something to 0?
#endif

// Variables for expressing functionals, these are redundant because
// different functionals have different needs.
template<class T>
struct densvars
{
  // Fills all density variables that can be filled from vars. Length of d
  // depends on vars.
  densvars(xc_functional_obj *parent, const T *d)
  {
    this->parent = parent;
    switch (parent->vars)
      {
      case XC_A_GAA:
	gaa = d[1];
	gab = 0;
	gbb = 0;
	gnn  = gaa;
	gss  = gaa;
	gns  = gaa;
      case XC_A:
	b = 0;
	n = a;
	s = a;
	break;
      case XC_A_B_GAA_GAB_GBB:
	gaa = d[2];
	gab = d[3];
	gbb = d[4];
	gnn  = gaa + 2*gab + gbb; 
	gss  = gaa - 2*gab + gbb;
	gns  = gaa - gbb;
      case XC_A_B:
	a = d[0];
	b = d[1];
	n = a+b;
	s = a-b;
	break;
      default:
	xcint_die("Illegal vars value in densvars()",parent->vars);	
      }
    zeta = s/n;
    r_s = pow(3.0/(n*4.0*M_PI),1.0/3.0); // (3/4pi)^1/3*n^(-1/3) !check
    n_m13 = pow(n,-1.0/3.0);
    a_43 = pow(a,4.0/3.0);
    b_43 = pow(b,4.0/3.0);
  }

  const xc_functional_obj *parent;
  double get_param(enum xc_parameter p) const
  {
    return parent->settings[p];
  }
  
  T a, b, gaa, gab, gbb;
  /* na+nb, na-nb, (grad n)^2, (grad n).(grad s), (grad s)^2 */
  T n, s, gnn, gns, gss; 

  T tau, taua, taub; // Kinetic energy densities.

  T lapa,lapb; // Density Laplacians

  T zeta; //s/n
  T r_s; // (3/4pi)^1/3*n^(-1/3)
  T n_m13; // pow(n,-1.0/3.0)
  T a_43, b_43; // pow(a,4.0/3.0), pow(b,4.0/3.0)
};
