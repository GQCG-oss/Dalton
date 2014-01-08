/* 
   Copyright 2012 Ulf Ekstrom, University of Oslo
   <uekstrom@gmail.com> 
   This code can be copied and modified under the LGPL 
   version 3 or later, see http://www.gnu.org/copyleft/lesser.html
*/

#include <cmath>
#include <cstdlib>
#include <cassert>
#include "taylor.hpp"

template<class T>
T bond_length(const T x[])
{
  return sqrt((x[0]-x[3])*(x[0]-x[3]) +
	      (x[1]-x[4])*(x[1]-x[4]) +
	      (x[2]-x[5])*(x[2]-x[5]));
}

template<class T>
T bond_angle(const T x[])
{
  T r12 = bond_length(x);
  T r23 = bond_length(x+3);
  T dot = 0;
  for (int i=0;i<3;i++)
    dot += (x[i] - x[3+i])*(x[6+i] - x[3+i]);
  return acos(dot/(r12*r23));
}

template<class T>
T dihedral(const T x[])
{
  T b1[3], b2[3], b3[3];
  for (int i=0;i<3;i++)
    {
      b1[i] = (x[3+i] - x[i]);
      b2[i] = (x[6+i] - x[3+i]);
      b3[i] = (x[9+i] - x[6+i]);
    }
  T yy = 
      b1[0]*(b2[1]*b3[2] - b2[2]*b3[1])
    + b1[1]*(b2[2]*b3[0] - b2[0]*b3[2])
    + b1[2]*(b2[0]*b3[1] - b2[1]*b3[0]);
  yy *= sqrt(b2[0]*b2[0] + b2[1]*b2[1] + b2[2]*b2[2]);
  T xx =       
      (b1[1]*b2[2] - b1[2]*b2[1])*(b2[1]*b3[2] - b2[2]*b3[1])
    + (b1[2]*b2[0] - b1[0]*b2[2])*(b2[2]*b3[0] - b2[0]*b3[2])
    + (b1[0]*b2[1] - b1[1]*b2[0])*(b2[0]*b3[1] - b2[1]*b3[0]);
  T res = atan(yy/xx);
  if (yy > 0 and xx < 0)
    res += M_PI;
  else if (yy < 0 and xx < 0)
    res -= M_PI;
  return res;
}

template<int N>
static void dispatch(int type,
		     double *q,
		     const double *x)
{
  const taylor<double,1,N> *in;
  taylor<double,1,N> *out;
  in = reinterpret_cast<const taylor<double,1,N> *>(x);
  out = reinterpret_cast<taylor<double,1,N> *>(q);
  switch(type)
    {
    case 1:
      *out = bond_length(in);
      break;
    case 2:
      *out = bond_angle(in);
      break;
    case 3:
      *out = dihedral(in);
      break;
    default:
      assert(type >= 1 && "Invalid coordiate type");      
      assert(type <= 3 && "Invalid coordiate type");      
    }
}

extern "C"
void dqdx(int type,
	  int order,
	  double q[],
	  const double x[])
{
  switch(order)
    {
    case 0:
      return dispatch<0>(type,q,x);
    case 1:
      return dispatch<1>(type,q,x);
    case 2:
      return dispatch<2>(type,q,x);
    case 3:
      return dispatch<3>(type,q,x);
    case 4:
      return dispatch<4>(type,q,x);
    case 5:
      return dispatch<5>(type,q,x);
    case 6:
      return dispatch<6>(type,q,x);
    case 7:
      return dispatch<7>(type,q,x);
    case 8:
      return dispatch<8>(type,q,x);
    case 9:
      return dispatch<9>(type,q,x);
    case 10:
      return dispatch<10>(type,q,x);
    case 11:
      return dispatch<11>(type,q,x);
    case 12:
      return dispatch<12>(type,q,x);
    case 13:
      return dispatch<13>(type,q,x);
    case 14:
      return dispatch<14>(type,q,x);
    case 15:
      return dispatch<15>(type,q,x);
    case 16:
      return dispatch<16>(type,q,x);
    case 17:
      return dispatch<17>(type,q,x);
    case 18:
      return dispatch<18>(type,q,x);
    case 19:
      return dispatch<19>(type,q,x);
    case 20:
      return dispatch<20>(type,q,x);
    default:
      assert(order >= 0 && "Invalid order in dqdx");
      assert(order <= 20 && "Invalid order in dqdx");
    };
}

#ifndef NO_FORTRAN

extern "C"
void dqdx_(int *type, int *order,
	   double *q,
	   double *x)
{
  dqdx(*type,*order,q,x);
}

#endif


#ifdef DQDX_STANDALONE

#include <iostream>
using namespace std;

void exit_on_error(istream &src)
{
  if (src.fail())
    {
      cerr << "Error reading input, quitting.\n";
      exit(-1);
    }
}

template<int N>
void read_eval_print(istream &src)
{
  int type;
  taylor<double, 1, N> x[12], res;
  src >> type;
  exit_on_error(src);
  if (type == 1) // bond length
    {
      for (int i=0;i<6;i++)
	for (int j=0;j<=N;j++)
	  {
	    src >> x[i][j];
	    exit_on_error(src);
	  }
      res = bond_length(x);
    }
  else if (type == 2)
    {
      for (int i=0;i<9;i++)
	for (int j=0;j<=N;j++)
	  {
	    src >> x[i][j];
	    exit_on_error(src);
	  }
      res = bond_angle(x);
    }
  else if (type == 3)
    {
      for (int i=0;i<12;i++)
	for (int j=0;j<=N;j++)
	  {
	    src >> x[i][j];
	    exit_on_error(src);
	  }
      res = dihedral(x);
    }
  else
    {
      cerr << "Error, invalid coordinate type " << type << " , quitting.\n";
      exit(-2);
    }
  for (int i=0;i<=N;i++)
    cout << res[i] << " ";
  cout << endl;
}


int main(int argc, const char **argv)
{
  if (argc > 1)
    {
      cout << "This program reads from standard input and outputs to standard output." << endl;
      cout << "Input:" << endl;
      cout << "ORDER TYPE X10 X11 .. X1ORDER X20 .. X2ORDER" << endl;
      cout << "Where ORDER is the order of derivative, TYPE is one of" << endl;
      cout << "1 - Bond length (X1..X6)" << endl;
      cout << "2 - Bond angle (X1..X9)" << endl;
      cout << "3 - Dihedral angle (X1..X12)" << endl;
      cout << "After TYPE follows ORDER coefficients of each cartesian coordinate" << endl;
      cout << "For example if ORDER = 1 and TYPE = 1 the input would be" << endl;
      cout << "1 1 X10 X11 X20 X21 X30 X31 X40 X41 X50 X51 X60 X61" << endl;
      cout << "where X20 is the constant term of the y coordinate of the first atom," << endl;
      cout << "X41 is the linear term of the x coordinate of the second atom, etc." << endl;
      return EXIT_SUCCESS;
    }
  else 
    {
      cout.precision(16);
      while (true)
	{
	  int n = -1;
	  cin >> n;
	  if (cin.fail())
	    {
	      if (cin.eof())
		return EXIT_SUCCESS;
	      else
		cerr << "Input error reading derivative order\n";
	      return EXIT_FAILURE;
	    }
	  switch (n)
	    {
	    case 0:
	      read_eval_print<0>(cin);
	      break;
	    case 1:
	      read_eval_print<1>(cin);
	      break;
	    case 2:
	      read_eval_print<2>(cin);
	      break;
	    case 3:
	      read_eval_print<3>(cin);
	      break;
	    case 4:
	      read_eval_print<4>(cin);
	      break;
	    case 5:
	      read_eval_print<5>(cin);
	      break;
	    case 6:
	      read_eval_print<6>(cin);
	      break;
	    case 7:
	      read_eval_print<7>(cin);
	      break;
	    case 8:
	      read_eval_print<8>(cin);
	      break;
	    case 9:
	      read_eval_print<9>(cin);
	      break;
	    case 10:
	      read_eval_print<10>(cin);
	      break;
	    case 11:
	      read_eval_print<11>(cin);
	      break;
	    case 12:
	      read_eval_print<12>(cin);
	      break;
	    default:
	      cerr << "Requested order " << n << " too high, quitting\n";
	      return EXIT_FAILURE;
	    }
	}
    }
}

#endif
