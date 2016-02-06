#if defined(VAR_MPI)
#include "mpi.h"
#endif

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <map>
#include <queue>
#include "consts.hpp"
#include "ababps.hpp"
#include "cache64.hpp"

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
    #define omp_get_num_threads() 1
#endif

using namespace std;
using namespace LibAngularVPD;
using namespace LibIGammaVPD;

struct point {
    double x;
    double y;
    double z;
};

struct vector3 {
    double x;
    double y;
    double z;
};

SHTerm LibAngularVPD::CartList[LMAX+2][MMAX];
SphHar LibAngularVPD::SHList  [LMAX+2][2*LMAX+1];

static double AngularIntegral(int nx, int ny, int nz) {
    //gamma[n+1/2]
    const double hGamma[] =
    {1.77245385090551603, 0.886226925452758014, 1.32934038817913702,
    3.32335097044784255, 11.6317283965674489, 52.3427777845535202,
    287.885277815044361, 1871.25430579778835, 14034.4072934834126,
    119292.461994609007, 1.13327838894878557e6,
     1.18994230839622485e7, 1.36843365465565857e8,
     1.71054206831957322e9, 2.30923179223142384e10,
     3.34838609873556457e11, 5.18999845304012508e12,
     8.56349744751620639e13, 1.49861205331533612e15,
     2.77243229863337182e16, 5.40624298233507504e17,
     1.10827981137869038e19, 2.38280159446418433e20,
     5.36130358754441473e21, 1.25990634307293746e23,
     3.08677054052869678e24, 7.87126487834817680e25,
     2.08588519276226685e27, 5.73618428009623384e28,
     1.63481251982742664e30, 4.82269693349090860e31,
     1.47092256471472712e33};

   if (nx%2 || ny%2 || nz%2) return 0;
   return 2*(hGamma[nx/2]*hGamma[ny/2]*hGamma[nz/2])/hGamma[(nx+ny+nz+2)/2];
}


static inline void LenzUpdate(double & A, double & Ap, const double a, const double b) {
    double tmp = A;
    A = b*A + a*Ap;
    Ap = tmp;
}

//Boys' function
//Table[Integrate[t^(2n) exp[-x t^2],{x,0,1}],{n,0,...}]
static inline void calcF(double * F, int n, double x) {
    double expx = exp(-x);

    //compute gamma(n + 1/2, x) / x^(n+1/2) and the lower ones by downward recursion
    if (x < double(n) + 1.5) {
        double ac  = 1./(double(n)+0.5);
        double sum = ac;

        //loop until machine precision convergence
        for (int i=1; (sum+ac)!=sum; ++i) {
            ac *= (x)/(double(n+i)+0.5);
            sum +=ac;
        }

        F[n] = 0.5 * sum;

        //downward recursion
        for (int i=n-1; i>=0; --i)
            F[i] = (2*x*F[i+1] + 1.)/double(2*i+1);

        for (int i=n; i>=0; --i)
            F[i] *= expx;
    }

    //compute F0[x] from continued fraction expansion and the rest by upward recursion
    else {
        double Ap = 1;
        double A  = 0; //b0=0

        double Bp = 0;
        double B  = 1;

        LenzUpdate(A,Ap,1.,x);
        LenzUpdate(B,Bp,1.,x);

        //OPT: optimal loop unrolling
        //OPT: function inlining
        for (int i=0; i<10; ++i) {
            LenzUpdate(A,Ap,double(i)+0.5,1.);
            LenzUpdate(B,Bp,double(i)+0.5,1.);

            LenzUpdate(A,Ap,double(i+1),x);
            LenzUpdate(B,Bp,double(i+1),x);
        }

        F[0] = sqrt(PI/(4*x)) -0.5 * expx * (A/B);

        //upward recursion
        for (int i=0; i<n; ++i)
            F[i+1] = (double(2*i+1)*F[i] - expx) / (2*x);
    }
}


namespace LibIGammaVPD {

    GammaFunction::GammaFunction() {
        gamma_table = NULL;
        gamma_table_vals = NULL;
    }

    GammaFunction::~GammaFunction() {
        free(gamma_table);
        free(gamma_table_vals);
    }

    //nueva funcion gamma
    void GammaFunction::CalcGammas(double * F, int n, double x) const{
        //busca el valor mas cercano de F
        //cout << vg_max << " " << vg_step << endl;
        if (x>vg_max-vg_step) {
            F[0] = 0.5*sqrt(PI/x);

            double ix = 0.5/x;
            for (int i=1; i<=n; ++i) {
                F[i] = F[i-1]*ix*double(2*i-1);
            }
        }
        else {
            double p = ivg_step*(x-vg_min);
            int pos = int(p+0.5);
            double x0 = vg_min + vg_step*double(pos);
            double Ax = x0-x;

            double Axn[NEXP];
            Axn[0] = 1;
            for (int i=1; i<NEXP; ++i) {
                Axn[i] = (Ax * Axn[i-1])/double(i);
            }

            {
                double sum = 0;
                for (int j=0; j<NEXP; ++j)
                    sum += gamma_table[pos][n+j+1] * Axn[j];
                F[n] = sum;
            }

            if (n>0)
            {
                double expx = 1;

                for (int i=1; i<NEXP; ++i) {
                    expx += Axn[i];
                }

                expx *= gamma_table[pos][0];

                for (int i=n-1; i>=0; --i)
                    F[i] = (2*x*F[i+1] + expx)/double(2*i+1);
            }


        }
    }


    void GammaFunction::InitIncompleteGammaTable() {
        if (gamma_table!=NULL) return;

        //init arrays and array positions
        posix_memalign((void**)(&gamma_table)     , CACHE_LINE_SIZE, NVIG     *sizeof(double*));
        posix_memalign((void**)(&gamma_table_vals), CACHE_LINE_SIZE, NVIG*NGD2*sizeof(double));

        gamma_table      = new double*[NVIG];
        gamma_table_vals = new double[NGD*NVIG];

        for (int i=0; i<NVIG; ++i) gamma_table[i] = gamma_table_vals + i*NGD;

        double Av = (vg_max-vg_min) / double (NVIG);
        for (int i=0; i<NVIG; ++i) {
            double v = vg_min + i*Av;

            gamma_table[i][0] = exp(-v);
            calcF(gamma_table[i]+1, 4*LMAX+1+NEXP, v);
        }
    }

    void GammaFunction::InitIncompleteGammaTable(int L) {
        if (gamma_table!=NULL) return;

        //init arrays and array positions
        posix_memalign((void**)(&gamma_table)     , CACHE_LINE_SIZE, NVIG     *sizeof(double*));
        posix_memalign((void**)(&gamma_table_vals), CACHE_LINE_SIZE, NVIG*NGD2*sizeof(double));

        for (int i=0; i<NVIG; ++i) gamma_table[i] = gamma_table_vals + i*NGD2;

        double Av = (vg_max-vg_min) / double (NVIG);
        for (int i=0; i<NVIG; ++i) {
            double v = vg_min + i*Av;

            double vv[32];
            calcF(vv, L+7, v);

            //first value is the exponential
            gamma_table[i][0] = exp(-v);
            //rest of values are the derivatives up to order 6
            for (int j=0; j<7; ++j) gamma_table[i][j+1] = vv[L+j];
        }
    }


    GammaFunction IncompleteGamma;
    GammaFunction IncompleteGammas[4*LMAX+3];


    void InitIncompleteGammas() {
        IncompleteGamma.InitIncompleteGammaTable();

        for (int l=0; l<=4*LMAX+2; ++l) IncompleteGammas[l].InitIncompleteGammaTable(l);
    }



}




void LibAngularVPD::SphHar::SetNPS(int n) {
    nps = n;
}

void LibAngularVPD::SphHar::SetT(int t, double c, int x, int y, int z) {
    T[t].cN = c;
    T[t].nx = x;
    T[t].ny = y;
    T[t].nz = z;

    //find the term
    //int l = x+y+z;
    int mm = nmC[l];

    for (int m=0; m<mm; ++m) {
        if ( (CartList[l][m].nx == x) && (CartList[l][m].ny == y) && (CartList[l][m].nz == z) ) {
            T[t].nc = m;
            return;
        }
    }

}


//initialize solid spherical harmonics tables
void LibAngularVPD::InitSHList() {

    for (int l=0; l<=LMAX; ++l) {
        for (int m=0; m<=2*l; ++m){
            SHList[l][m].l = l;
            SHList[l][m].m = m;
        }
    }

    /*
    //SP; ignore at the moment
    {
        for (int m=0; m<4; ++m){
            SHList[LSP][m].l = LSP;
            SHList[LSP][m].m = m;
        }
    }
    */

    //s
    {
        SHList[0][0].SetNPS(1);
        SHList[0][0].SetT(0, 1, 0,0,0);
    }

    //p

    /*
    if (LMAX > 0) {
        SHList[1][0].SetNPS(1); //y
        SHList[1][0].SetT(0, 1, 0,1,0);

        SHList[1][1].SetNPS(1); //z
        SHList[1][1].SetT(0, 1, 0,0,1);

        SHList[1][2].SetNPS(1); //x
        SHList[1][2].SetT(0, 1, 1,0,0);
    }
    */


    //sp
    /*
    if (LMAX > 0) {
        SHList[LSP][0].SetNPS(1);
        SHList[LSP][0].SetT(0, 1, 0,0,0);

        SHList[LSP][1].SetNPS(1); //y
        SHList[LSP][1].SetT(0, 1, 0,1,0);

        SHList[LSP][2].SetNPS(1); //z
        SHList[LSP][2].SetT(0, 1, 0,0,1);

        SHList[LSP][3].SetNPS(1); //x
        SHList[LSP][3].SetT(0, 1, 1,0,0);
    }
    */

    //p not in normal order!
    if (LMAX > 0) {
        SHList[1][0].SetNPS(1); //x
        SHList[1][0].SetT(0, 1, 1,0,0);

        SHList[1][1].SetNPS(1); //y
        SHList[1][1].SetT(0, 1, 0,1,0);

        SHList[1][2].SetNPS(1); //z
        SHList[1][2].SetT(0, 1, 0,0,1);
    }


    //d
    if (LMAX > 1) {
        SHList[2][0].SetNPS(1); //2i xy
        SHList[2][0].SetT(0, 1, 1,1,0);

        SHList[2][1].SetNPS(1); //i yz
        SHList[2][1].SetT(0, 1, 0,1,1);

        SHList[2][2].SetNPS(3); //2z^2-x^2-y^2 = 3z^2 - r^2
        SHList[2][2].SetT(0,  2, 0,0,2);
        SHList[2][2].SetT(1, -1, 2,0,0);
        SHList[2][2].SetT(2, -1, 0,2,0);

        SHList[2][3].SetNPS(1); //xz
        SHList[2][3].SetT(0, 1, 1,0,1);

        SHList[2][4].SetNPS(2); //x^2-y^2
        SHList[2][4].SetT(0,  1, 2,0,0);
        SHList[2][4].SetT(1, -1, 0,2,0);
    }

    //f
    if (LMAX > 2) {
        // i (3x^2 - y^2) y
        SHList[3][0].SetNPS(2);
        SHList[3][0].SetT(0, -1, 0,3,0);
        SHList[3][0].SetT(1,  3, 2,1,0);

        //2i z xy
        SHList[3][1].SetNPS(1);
        SHList[3][1].SetT(0, 1, 1,1,1);

        // i (4z^2-x^2-y^2) y
        SHList[3][2].SetNPS(3);
        SHList[3][2].SetT(0,  4, 0,1,2);
        SHList[3][2].SetT(1, -1, 2,1,0);
        SHList[3][2].SetT(2, -1, 0,3,0);

        // (2z^2-3x^2-3y^2) z m=0
        SHList[3][3].SetNPS(3);
        SHList[3][3].SetT(0,  2, 0,0,3);
        SHList[3][3].SetT(1, -3, 2,0,1);
        SHList[3][3].SetT(2, -3, 0,2,1);

        // x(4z^2-x^2-y^2)
        SHList[3][4].SetNPS(3);
        SHList[3][4].SetT(0,  4, 1,0,2);
        SHList[3][4].SetT(1, -1, 3,0,0);
        SHList[3][4].SetT(2, -1, 1,2,0);

        // z(x^2-y^2)
        SHList[3][5].SetNPS(2);
        SHList[3][5].SetT(0,  1, 2,0,1);
        SHList[3][5].SetT(1, -1, 0,2,1);

        // x^3 - 3xy^2
        SHList[3][6].SetNPS(2);
        SHList[3][6].SetT(0,  1, 3,0,0);
        SHList[3][6].SetT(1, -3, 1,2,0);
    }

    //g
    if (LMAX > 3) {
        //i xy (x^2 - y^2)
        SHList[4][0].SetNPS(2);
        SHList[4][0].SetT(0,  1, 3,1,0);
        SHList[4][0].SetT(1, -1, 1,3,0);

        //i z y (3x^2 - y^2)
        SHList[4][1].SetNPS(2);
        SHList[4][1].SetT(0, 3, 2,1,1);
        SHList[4][1].SetT(1,-1, 0,3,1);

        //i x y (7z^2 - r^2)
        SHList[4][2].SetNPS(3);
        SHList[4][2].SetT(0,  6, 1,1,2);
        SHList[4][2].SetT(1, -1, 3,1,0);
        SHList[4][2].SetT(2, -1, 1,3,0);

        //i y z (7z^2 -3r^2)
        SHList[4][3].SetNPS(3);
        SHList[4][3].SetT(0,  4, 0,1,3);
        SHList[4][3].SetT(1, -3, 2,1,1);
        SHList[4][3].SetT(2, -3, 0,3,1);

        //3x^4 + 6x^2y^2 + 3y^4 - 24x^2z^2 - 24y^2z^2 + 8z^4 = 35z^4 - 30z^2 r^2 + 3r^4
        //= 35z^4 - 30z^2 r^2 + 3r^4
        //m=0
        SHList[4][4].SetNPS(6);
        SHList[4][4].SetT(0,   3, 4,0,0);
        SHList[4][4].SetT(1,   6, 2,2,0);
        SHList[4][4].SetT(2,   3, 0,4,0);
        SHList[4][4].SetT(3, -24, 2,0,2);
        SHList[4][4].SetT(4, -24, 0,2,2);
        SHList[4][4].SetT(5,   8, 0,0,4);

        //x z (7z^2 - 3r^2)
        SHList[4][5].SetNPS(3);
        SHList[4][5].SetT(0,  4, 1,0,3);
        SHList[4][5].SetT(1, -3, 3,0,1);
        SHList[4][5].SetT(2, -3, 1,2,1);

        //(x^2-y^2)*(7z^2-r^2) = -x^4 + y^4 - 6y^2 z^2 + 6x^2 z^2
        SHList[4][6].SetNPS(4);
        SHList[4][6].SetT(0, -1, 4,0,0);
        SHList[4][6].SetT(1,  1, 0,4,0);
        SHList[4][6].SetT(2, -6, 0,2,2);
        SHList[4][6].SetT(3,  6, 2,0,2);

        //x z (x^2 - 3y^2)
        SHList[4][7].SetNPS(2);
        SHList[4][7].SetT(0,  1, 3,0,1);
        SHList[4][7].SetT(1, -3, 1,2,1);

        //x^4 - 6x^2y^2 + y^4
        SHList[4][8].SetNPS(3);
        SHList[4][8].SetT(0,  1, 4,0,0);
        SHList[4][8].SetT(1,  1, 0,4,0);
        SHList[4][8].SetT(2, -6, 2,2,0);
    }

    //h
    if (LMAX > 4) {
        //i (5 x^4 y - 10 x^2 y^3 + y^5)
        SHList[5][0].SetNPS(3);
        SHList[5][0].SetT(0,   1, 0,5,0);
        SHList[5][0].SetT(1, -10, 2,3,0);
        SHList[5][0].SetT(2,   5, 4,1,0);

        //i z (x^3 y - x y^3)
        SHList[5][1].SetNPS(2);
        SHList[5][1].SetT(0,  1, 3,1,1);
        SHList[5][1].SetT(1, -1, 1,3,1);

        //i (y^5 - 2 y^3 x^2 - 3 y x^4 - 8 y^3 z^2 + 24 y x^2 z^2)
        SHList[5][2].SetNPS(5);
        SHList[5][2].SetT(0,  1, 0,5,0);
        SHList[5][2].SetT(1, -2, 2,3,0);
        SHList[5][2].SetT(2, -3, 4,1,0);
        SHList[5][2].SetT(3, -8, 0,3,2);
        SHList[5][2].SetT(4, 24, 2,1,2);

        //i (-x^3 yz - x y^3 z + 2xyz^3)
        SHList[5][3].SetNPS(3);
        SHList[5][3].SetT(0, -1, 3,1,1);
        SHList[5][3].SetT(1, -1, 1,3,1);
        SHList[5][3].SetT(2,  2, 1,1,3);

        //i (x^4 y + 2 x^2 y^3 + y^5 - 12 x^2 y z^2 - 12 y^3 z^2 + 8 y z^4)
        SHList[5][4].SetNPS(6);
        SHList[5][4].SetT(0,   1, 0,5,0);
        SHList[5][4].SetT(1,   2, 2,3,0);
        SHList[5][4].SetT(2,   1, 4,1,0);
        SHList[5][4].SetT(3, -12, 0,3,2);
        SHList[5][4].SetT(4, -12, 2,1,2);
        SHList[5][4].SetT(5,   8, 0,1,4);

        //15 x^4 z + 30 x^2 y^2 z + 15 y^4 z - 40 x^2 z^3 - 40 y^2 z^3 + 8 z^5
        //63 z^5 - 70 z^3 r^2 + 15 z r^4
        // m=0
        SHList[5][5].SetNPS(6);
        SHList[5][5].SetT(0,   8, 0,0,5);
        SHList[5][5].SetT(1, -40, 2,0,3);
        SHList[5][5].SetT(2, -40, 0,2,3);
        SHList[5][5].SetT(3,  15, 4,0,1);
        SHList[5][5].SetT(4,  30, 2,2,1);
        SHList[5][5].SetT(5,  15, 0,4,1);

        //x^5 + 2 x^3 y^2 + x y^4 - 12 x^3 z^2 - 12 x y^2 z^2 + 8 x z^4
        SHList[5][6].SetNPS(6);
        SHList[5][6].SetT(0,   1, 5,0,0);
        SHList[5][6].SetT(1,   2, 3,2,0);
        SHList[5][6].SetT(2,   1, 1,4,0);
        SHList[5][6].SetT(3, -12, 3,0,2);
        SHList[5][6].SetT(4, -12, 1,2,2);
        SHList[5][6].SetT(5,   8, 1,0,4);

        // (-x^4 z + y^4 z + 2 x^2 y^3 - 2 y^2 x^3)
        SHList[5][7].SetNPS(4);
        SHList[5][7].SetT(0, -1, 4,0,1);
        SHList[5][7].SetT(1,  1, 0,4,1);
        SHList[5][7].SetT(2,  2, 2,0,3);
        SHList[5][7].SetT(3, -2, 0,2,3);

        //3 x y^4 + 2 x^3 y^2 - x^5 - 24 x y^2 z^2 + 8 x^3 z^2
        SHList[5][8].SetNPS(5);
        SHList[5][8].SetT(0, -1, 5,0,0);
        SHList[5][8].SetT(1,  2, 3,2,0);
        SHList[5][8].SetT(2,  3, 1,4,0);
        SHList[5][8].SetT(3,  8, 3,0,2);
        SHList[5][8].SetT(4,-24, 1,2,2);

        //z (x^4 - 6 x^2 y^2 + y^4)
        SHList[5][9].SetNPS(3);
        SHList[5][9].SetT(0,  1, 4,0,1);
        SHList[5][9].SetT(1,  1, 0,4,1);
        SHList[5][9].SetT(2, -6, 2,2,1);

        //x^5 - 10 x^3 y^2 + 5 x y^4
        SHList[5][10].SetNPS(3);
        SHList[5][10].SetT(0,   1, 5,0,0);
        SHList[5][10].SetT(1, -10, 3,2,0);
        SHList[5][10].SetT(2,   5, 1,4,0);
    }


    //normalization
    for (int l=0; l<=LMAX+1; ++l) {
        for (int m=0; m<nmS[l]; ++m) {
            double sum = 0;

            for (int p=0; p<SHList[l][m].nps; ++p)
                for (int q=0; q<SHList[l][m].nps; ++q)
                    sum += SHList[l][m].T[p].cN * SHList[l][m].T[q].cN *
                    AngularIntegral( SHList[l][m].T[p].nx + SHList[l][m].T[q].nx ,
                                     SHList[l][m].T[p].ny + SHList[l][m].T[q].ny ,
                                     SHList[l][m].T[p].nz + SHList[l][m].T[q].nz);

            sum = 1/sqrt(sum);

            for (int p=0; p<SHList[l][m].nps; ++p)
                SHList[l][m].T[p].cN *= sum;

        }
    }
}


//initialize cartesian tables
void LibAngularVPD::InitCartList () {

    for (int l=0; l<=LMAX; ++l) {

        int t = 0;

        for (int a=0; a<=l; ++a) {
            for (int b=0; b<=l-a; ++b) {
                CartList[l][t].nx = l-a-b;
                CartList[l][t].ny = b;
                CartList[l][t].nz = a;
                CartList[l][t].cN = 1./sqrt( AngularIntegral(2*(l-a-b), 2*b, 2*a)); // hGamma[i+1]/(2*hGamma[a] * hGamma[b] * hGamma[i-a-b]) );

                ++t;
            }
        }
    }

}





/*
    Overloaded MPI short-hand
*/

inline void mpi_bcast(int & val) {
  #if defined(VAR_MPI)
     MPI_Bcast (&val, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  #endif
}

inline void mpi_bcast(int * ival, int len) {
  #if defined(VAR_MPI)
     MPI_Bcast (ival, len, MPI_INTEGER, 0, MPI_COMM_WORLD);
  #endif
}

inline void mpi_bcast(double & val) {
  #if defined(VAR_MPI)
    MPI_Bcast (&val, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif
}


inline void mpi_bcast(double * pval, int len) {
  #if defined(VAR_MPI)
    MPI_Bcast (pval, len, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif
}

inline int mpi_nnodes() {
  #if defined(VAR_MPI)
    int numprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    return numprocs;
  #else
    return 1;
  #endif
}

inline int mpi_myrank() {
  #if defined(VAR_MPI)
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
  #else
    return 0;
  #endif
}


void Bcast_Potential(point * NN, double * QQ, double * Rq, vector3 * DD, double * Rd, int Ncenters) {
    mpi_bcast((double*)NN, 3*Ncenters);
    mpi_bcast(QQ, Ncenters);
    mpi_bcast(Rq, Ncenters);
    mpi_bcast((double*)DD, 3*Ncenters);
    mpi_bcast(Rd, Ncenters);
}

void Bcast_Basis (int * nhkt, int * nuco, int * nrco, int * jstrt, double * cent, double * priccf, double * priexp) {

    mpi_bcast(nhkt,  MXSHEL);
    mpi_bcast(nuco,  MXSHEL);
    mpi_bcast(nrco,  MXSHEL);
    mpi_bcast(jstrt, MXSHEL);

    mpi_bcast(cent, 6*MXSHEL);
    mpi_bcast(priccf, MXPRIM*MXCONT);
    mpi_bcast(priexp, MXPRIM);
}

void ReduceV(double * V, double * VV, int lenV2) {
    #if defined(VAR_MPI)
    int ierr;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(VV, V, lenV2, MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    #else
    for (int i=0; i<lenV2; ++i) V[i] = VV[i];
    #endif
}






struct aCenter {
    point P;
    int   n;
    int   pos;
};

struct aCenterDuple {
    double R2;
    int posA;
    int posB;
};

//this is just an implementation of quicksort
void FlashSort(aCenterDuple * keys, int N) {

    //for short lists default to selection sort
    if (N<32) {
        for (int i=0; i<N; ++i)  {
            double best = keys[i].R2;
            int   k   = i;

            for (int j=i+1; j<N; ++j)
                if (keys[j].R2<best) {best = keys[j].R2; k = j;}

            swap(keys[i], keys[k]);
        }
        return;
    }

    //median-of-three pivot
    double kpiv; {
        double & k1 = keys[0].R2;
        double & k2 = keys[N/2].R2;
        double & k3 = keys[N-1].R2;


        if      (k1<=k2 && k2<=k3) kpiv = k2;
        else if (k3<=k2 && k2<=k1) kpiv = k2;
        else if (k1<=k3 && k3<=k2) kpiv = k3;
        else if (k2<=k3 && k3<=k1) kpiv = k3;
        else                       kpiv = k1;
        //else if (k2<k1 && k1<k3) kpiv = k1;
        //else if (k3<k1 && k1<k2) kpiv = k1;
    }


    //count number of elements lower than the pivot
    int nl=0;
    int ng=0;
    for (int i=0; i<N; ++i) {
        if      (keys[i].R2 < kpiv) ++nl;
        else if (keys[i].R2 > kpiv) ++ng;
    }

    int i,j;

    i=0;
    j=nl;

    for (; i<nl; ++i) {
        //if the element does not belong in this sublist
        if (keys[i].R2>=kpiv) {
            for (; j<N; ++j) {
                //found a replacement
                if (keys[j].R2<kpiv) {
                    swap(keys[i], keys[j]);
                    break;
                }
            }
        }
    }

    j=N-ng;

    for (i=nl; i<N-ng; ++i) {
        //if the element does not belong in this sublist
        if (keys[i].R2>kpiv) {
            for (; j<N; ++j) {
                //found a replacement
                if (keys[j].R2==kpiv) {
                    swap(keys[i], keys[j]);
                    break;
                }
            }
        }
    }

    /*

    int i,j,k;

    i=0;
    j=nl;

    while (1) {
        //find the first useful place in the first half of the list
        while (i<nl && keys[i].R2< kpiv) ++i;
        if (i>=nl) break;
        while (j<N  && keys[j].R2>=kpiv) ++j;
        if (j>=N) break;

        swap(keys[i]    , keys[j]);
    }

    j = nl;
    k = N-ng;

    while (1) {
        //find the first useful place in the first half of the list
        while (j<N-ng && keys[j].R2==kpiv) ++j;
        if (j>=N-ng) break;
        while (k<N    && keys[k].R2> kpiv) ++k;
        if (k>=N) break;

        swap(keys[j]    , keys[k]);
    }
    */


    //FlashSort both sublists (skipping the pivot)
    FlashSort(keys         ,  nl);
    FlashSort(keys + (N-ng),  ng);
}



const int DPC = DOUBLES_PER_CACHE_LINE;

//shorthand
struct Xpotential {

    //atomic centers
    point * centers;

    //charge and "radius" of the gaussian
    double * Q;
    double * Rq;

    //dipole and "radius" of the gaussian
    vector3 * D;
    double * Rd;

    //number of centers to take into consideration
    int ncenters;
};

template <int L> struct Triple64 {
    cacheline64 V[L+1][L+1][L+1];

    Triple64() {
        for (int ex=0; ex<=L; ++ex)
            for (int ey=0; ey<=L-ex; ++ey)
                for (int ez=0; ez<=L-ex-ey; ++ez)
                    V[ex][ey][ez] = 0;
    }

    inline cacheline64 & operator()(int i, int j, int k) {
        return V[i][j][k];
    }

    inline cacheline64 operator()(int i, int j, int k) const {
        return V[i][j][k];
    }

};

template <int L> struct TripleX64 {
    cacheline64 V[L+1][L+1][L+1];

    inline cacheline64 & operator()(int i, int j, int k) {
        return V[i][j][k];
    }

    inline cacheline64 operator()(int i, int j, int k) const {
        return V[i][j][k];
    }

};


//computing the one center damped potential integral over Hermite polynomials
template <int L> void RIntegralX(TripleX64<L> & Rtuv, double k, double Rd, const cacheline64 * P, const point & N) {

    cacheline64 Fn[L+1];

    cacheline64 Rx, Ry, Rz;

    Rx.set(N.x);
    Ry.set(N.y);
    Rz.set(N.z);

    Rx -= P[0];
    Ry -= P[1];
    Rz -= P[2];


    cacheline64 R2;
    R2 = Rx*Rx + Ry*Ry + Rz*Rz;

    double ** gamma = IncompleteGamma.gamma_table;
    double *  gamma2 = IncompleteGammas[L].gamma_table_vals;

    //compute the auxilliary integrals
    {
        double kw2 = k / (1+k*Rd*Rd);
        cacheline64 kw2R2 = R2*kw2;

        for (int d=0; d<DPC; ++d) {
            if (kw2R2(d)>LibIGammaVPD::vg_max-LibIGammaVPD::vg_step) {
                double iR2 = 1/R2(d);

                Fn[0](d) = 0.5*sqrt(PI*iR2);

                for (int i=1; i<=L; ++i) {
                    Fn[i](d) = Fn[i-1](d)*iR2*double(2*i-1);
                }
            }
            else {
                double p = LibIGammaVPD::ivg_step*(kw2R2(d)-LibIGammaVPD::vg_min);
                int pos = int(p+0.5);
                double x0 = LibIGammaVPD::vg_min + LibIGammaVPD::vg_step*double(pos);
                double Ax1 = x0-kw2R2(d);
                double Ax2 = 0.5 * Ax1*Ax1;
                double Ax3 = 0.33333333333333333333 * Ax1*Ax2;

                Fn[L](d) = (gamma[pos][L+1] + gamma[pos][L+2] * Ax1 + gamma[pos][L+3] * Ax2 + gamma[pos][L+4] * Ax3);

                if (L>0) {
                    double expx = gamma[pos][0] * (1+Ax1+Ax2+Ax3);

                    for (int i=L-1; i>=0; --i)
                        Fn[i](d) = (2*kw2R2(d)*Fn[i+1](d) + expx)*(1./double(2*i+1));
                }

                double kw2n = sqrt(kw2);

                for (int n=0; n<=L; ++n) {
                    Fn[n](d) *= kw2n;
                    kw2n     *= 2*kw2;
                }
            }
        }
    }


    //scale R^{n}_{000}
    double ac = (2*PI)/(k*sqrt(k));
    for (int n=0; n<=L; ++n) Fn[n] *= ac;



    //use 1-center RRs to compute R_{tuv} reusing intermediate memory
    Rtuv(0,0,0) = Fn[L];
    for (int n=L-1; n>=0; --n) {
        for (int tuv=L-n; tuv>=1; --tuv) {
            for (int uv=0; uv<=tuv; ++uv) {
                for (int v=0; v<=uv; ++v) {
                    int t = tuv - uv;
                    int u =  uv -  v;

                    if(t==0) {
                        if(u==0) {
                            if(v==0) {
                                //(CANNOT HAPPEN)
                            }
                            else if (v==1) Rtuv(t,u,v) = Rz * Rtuv(t,u,v-1);
                            else           Rtuv(t,u,v) = Rz * Rtuv(t,u,v-1) -Rtuv(t,u,v-2)*(v-1);
                        }
                        else if (u==1) Rtuv(t,u,v) = Ry * Rtuv(t,u-1,v);
                        else           Rtuv(t,u,v) = Ry * Rtuv(t,u-1,v) -Rtuv(t,u-2,v)*(u-1);
                    }
                    else if (t==1) Rtuv(t,u,v) = Rx * Rtuv(t-1,u,v);
                    else           Rtuv(t,u,v) = Rx * Rtuv(t-1,u,v) -Rtuv(t-2,u,v)*(t-1);
                }
            }
        }

        //tuv = 0
        Rtuv(0,0,0) = Fn[n];
    }
}



// very inefficient implementation of the RRs (CTE)
template<int L>
static inline cacheline64 R2E(int ex, int ey, int ez, const cacheline64 (&PA)[3],     double iz, int px, int py, int pz, Triple64<L> & R) {

    if (px<0 || py<0 || pz<0) {
        cacheline64 zero;
        zero = 0;
        return zero;
    }

    if (ex>0) {
        return R2E(ex-1,ey,ez, PA,  iz, px-1,py,pz, R) * px    +
                R2E(ex-1,ey,ez, PA,  iz, px  ,py,pz, R) * PA[0] +
                R2E(ex-1,ey,ez, PA,  iz, px+1,py,pz, R) * iz;
    }
    if (ey>0) {
        return R2E(ex,ey-1,ez, PA,  iz, px,py-1,pz, R) * py +
                R2E(ex,ey-1,ez, PA,  iz, px,py  ,pz, R) * PA[1] +
                R2E(ex,ey-1,ez, PA,  iz, px,py+1,pz, R) * iz;
    }
    if (ez>0) {
        return R2E(ex,ey,ez-1, PA,  iz, px,py,pz-1, R) * pz +
                R2E(ex,ey,ez-1, PA,  iz, px,py,pz  , R) * PA[2] +
                R2E(ex,ey,ez-1, PA,  iz, px,py,pz+1, R) * iz;
    }

    return R(px,py,pz);
    //return R[px][py][pz];
}

//horizontal recurrence relation
static inline cacheline64 HRR(int ax, int ay, int az,   int bx, int by, int bz,  const cacheline64 (&BA)[3],  Triple64<2*LMAX> & R) {
    if (bx<0 || by<0 || bz<0) {
        cacheline64 zero;
        zero = 0;
        return zero;
    }

    if (bx>0) {
        return HRR(ax+1,ay,az,  bx-1,by,bz, BA, R) -
                HRR(ax,  ay,az,  bx-1,by,bz, BA, R) * BA[0];
    }
    if (by>0) {
        return HRR(ax,ay+1,az,  bx,by-1,bz, BA, R) -
                HRR(ax,ay,  az,  bx,by-1,bz, BA, R) * BA[1];
    }
    if (bz>0) {
        return HRR(ax,ay,az+1,  bx,by,bz-1, BA, R) -
                HRR(ax,ay,az,    bx,by,bz-1, BA, R) * BA[2];
    }

    return R(ax,ay,az);
    //return R[ax][ay][az];
}



typedef void (*CXPH)
                        (const cacheline64 * AA, int Ja, double ka, const double * wa,
                         const cacheline64 * BB, int Jb, double kb, const double * wb,
                         const Xpotential & Xpot,
                         const cacheline64 & ABr2,
                         cacheline64 * EEE
                         );
//the actual function
template <bool CMONO, bool CDIP, int Ltot>
void CalcExternalPotentialH
                        (const cacheline64 * AA, int Ja, double k1, const double * wa,
                         const cacheline64 * BB, int Jb, double k2, const double * wb,
                         const Xpotential & Xpot,
                         const cacheline64 & ABr2,
                         cacheline64 * EEE
                         ) {


    Triple64<Ltot> R;

    double k12, ik12, f1, f2;

    k12  = k1+k2;
    ik12 = 1/k12;
    f1   = k1*ik12;
    f2   = k2*ik12;

    cacheline64 P[3]; {
        P[0] = AA[0]*f1 + BB[0]*f2;
        P[1] = AA[1]*f1 + BB[1]*f2;
        P[2] = AA[2]*f1 + BB[2]*f2;
    }

    //accumulate contributions
    //************************
    TripleX64<Ltot> TRQ;

    // monopole contributions
    if (CMONO)
    for (int i=0; i<Xpot.ncenters; ++i) {

        //add the damped point charge contributions
        RIntegralX<Ltot>(TRQ, k12, Xpot.Rq[i], P, Xpot.centers[i]);

        for (int px=0; px<=Ltot; ++px) {
            for (int py=0; py<=Ltot-px; ++py) {
                for (int pz=0; pz<=Ltot-px-py; ++pz) {
                    R(px,py,pz) -= TRQ(px,py,pz) * Xpot.Q[i]; // + (Z charge) x - (e- charge) = -
                }
            }
        }
    }

    TripleX64<Ltot+1> TRD;

    if (CDIP)
    // dipole contributions
    for (int i=0; i<Xpot.ncenters; ++i) {
        //add the dipole moment contributions
        RIntegralX<Ltot+1>(TRD, k12, Xpot.Rd[i], P, Xpot.centers[i]);

        for (int px=0; px<=Ltot; ++px) {
            for (int py=0; py<=Ltot-px; ++py) {
                for (int pz=0; pz<=Ltot-px-py; ++pz) {
                    // - (e- charge) x - (r to p|q hermite decoupling) = +
                    R(px,py,pz) += TRD(px+1,py,pz) * Xpot.D[i].x;
                    R(px,py,pz) += TRD(px,py+1,pz) * Xpot.D[i].y;
                    R(px,py,pz) += TRD(px,py,pz+1) * Xpot.D[i].z;
                }
            }
        }
    }

    cacheline64 PA[3]; {
        PA[0] = P[0] - AA[0];
        PA[1] = P[1] - AA[1];
        PA[2] = P[2] - AA[2];
    }


    cacheline64 kR2; kR2 = ABr2 * (k1*k2*ik12);
    cacheline64 wab;
    for (int i=0; i<DPC; ++i)
        wab(i) = exp(-kR2(i));

    int i3=0;
    for (int ex=0; ex<=Ltot; ++ex)
        for (int ey=0; ey<=Ltot-ex; ++ey)
            for (int ez=0; ez<=Ltot-ex-ey; ++ez,++i3) {
                cacheline64 r2e;
                r2e = R2E(ex,ey,ez,PA, 0.5*ik12, 0,0,0, R);

                for (int ja=0; ja<Ja; ++ja) {
                    for (int jb=0; jb<Jb; ++jb) {
                        cacheline64 w12; w12 = wab *  (wa[ja] * wb[jb]);
                        EEE[i3*Ja*Jb + ja*Jb + jb] += w12 * r2e;
                    }
                }
            }


}


//check again with Mathematica!!!!
double CalcVps(const Xpotential & Xpot) {

    const double PI54  = 34.986836655249725693;
    const double SQRT2 =  1.4142135623730951;

    double vQ = 0;
    double vD = 0;

    // add contribution from gaussian charges
    for (int n=0; n<Xpot.ncenters; ++n) {
        double ik = (Xpot.Rq[n]*Xpot.Rq[n]);
        double sik = fabs(Xpot.Rq[n]);
        double N2 = 1./(PI*ik);
        double N  = N2*sqrt(N2);
        vQ += fabs(Xpot.Q[n]) * sqrt( (PI54/SQRT2)* sik*ik*ik * (N*N) );
    }

    //cout << "Q self-interaction: " << vQ << endl;

    // add contribution from gaussian dipoles
    for (int n=0; n<Xpot.ncenters; ++n) {
        double ik = (Xpot.Rd[n]*Xpot.Rd[n]);
        double sik = fabs(Xpot.Rd[n]);

        double D2 = Xpot.D[n].x*Xpot.D[n].x + Xpot.D[n].y*Xpot.D[n].y + Xpot.D[n].z*Xpot.D[n].z;

        double N2 = 1./(PI*ik);
        double N  = N2*sqrt(N2);
        //double N2 = (PI/2)*ik;
        //double N5 = (N2*sqrt(N2)*0.25*ik);

        vD += sqrt( D2 * (PI54*SQRT2/3)* sik*ik  * (N*N) );
    }

    //cout << "D self-interaction: " << vD << endl;

    return (vQ+vD)*(vQ+vD);
}


// modified quicksort;
// performs a check on the pivot and ignores upper list if not passed;
// reduces complexity from O(N^2 log N) to O(N log N)
void FlashSort(aCenterDuple * keys, int & N,   K2PS Kps, const eBasis & A, const eBasis & B, const eBasis2 & AB, double thresh) {

    //for short lists default to selection sort
    if (N<32) {
        for (int i=0; i<N; ++i)  {
            double best = keys[i].R2;
            int   k   = i;

            for (int j=i+1; j<N; ++j)
                if (keys[j].R2<best) {best = keys[j].R2; k = j;}

            swap(keys[i], keys[k]);
        }
        return;
    }

    //median-of-three pivot
    double kpiv; {
        double & k1 = keys[0].R2;
        double & k2 = keys[N/2].R2;
        double & k3 = keys[N-1].R2;


        if      (k1<=k2 && k2<=k3) kpiv = k2;
        else if (k3<=k2 && k2<=k1) kpiv = k2;
        else if (k1<=k3 && k3<=k2) kpiv = k3;
        else if (k2<=k3 && k3<=k1) kpiv = k3;
        else                       kpiv = k1;
    }


    //count number of elements lower than the pivot
    int nl=0;
    int ng=0;
    for (int i=0; i<N; ++i) {
        if      (keys[i].R2 < kpiv) ++nl;
        else if (keys[i].R2 > kpiv) ++ng;
    }

    int i,j;

    i=0;
    j=nl;

    for (; i<nl; ++i) {
        //if the element does not belong in this sublist
        if (keys[i].R2>=kpiv) {
            for (; j<N; ++j) {
                //found a replacement
                if (keys[j].R2<kpiv) {
                    swap(keys[i], keys[j]);
                    break;
                }
            }
        }
    }

    j=N-ng;

    for (i=nl; i<N-ng; ++i) {
        //if the element does not belong in this sublist
        if (keys[i].R2>kpiv) {
            for (; j<N; ++j) {
                //found a replacement
                if (keys[j].R2==kpiv) {
                    swap(keys[i], keys[j]);
                    break;
                }
            }
        }
    }

    /*

    int i,j,k;

    i=0;
    j=nl;

    while (1) {
        //find the first useful place in the first half of the list
        while (i<nl && keys[i].R2< kpiv) ++i;
        if (i>=nl) break;
        while (j<N  && keys[j].R2>=kpiv) ++j;
        if (j>=N) break;

        swap(keys[i]    , keys[j]);
    }

    j = nl;
    k = N-ng;

    while (1) {
        //find the first useful place in the first half of the list
        while (j<N-ng && keys[j].R2==kpiv) ++j;
        if (j>=N-ng) break;
        while (k<N    && keys[k].R2> kpiv) ++k;
        if (k>=N) break;

        swap(keys[j]    , keys[k]);
    }
    */



    //FlashSort both sublists (skipping the pivot)
    int np = N-nl-ng;

    {
        //check if it is necessary to compute at all
        double R =  sqrt(kpiv);
        int nABk2 = Kps(A, B, AB, R, thresh);

        if (nABk2==-1) {
            FlashSort(keys         ,  nl,    Kps, A, B, AB, thresh);
            N = nl; //return the length of the sorted sublist
        }
        else {
            FlashSort(keys         ,  nl,    Kps, A, B, AB, thresh);
            FlashSort(keys + (N-ng),  ng,    Kps, A, B, AB, thresh);
            N = nl+np+ng;  //return the length of the sorted sublist
        }
    }

}


//can be made smaller
struct XPtile {

    // a matrix with N>65536 would be 32Gb, which should never happen in DALTON
    // on the other hand, keeping this struct's size to a minimum should improve performance notably
    unsigned short int posA[DPC]; //16 bytes
    unsigned short int posB[DPC]; //16 bytes

    const eBasis * Abasis; // 8 bytes
    const eBasis * Bbasis; // 8 bytes

    nKinteracting * nKi;  // 8 bytes

    short int nUsed;  // <=DPC
    short int Pstart; // where to start counting primitives

    //these are the important ones for sorting
    short int Ltot; // La+Lb
    short int nP2;  // number of primitives


    //to establish a well-ordered sequence
    bool operator<(const XPtile & rhs) const {
        if (Ltot!=rhs.Ltot) return (Ltot<rhs.Ltot); // begin with high L
        if (nP2 !=rhs.nP2)  return (nP2 <rhs.nP2);  // then highest K2

        //if (Pstart !=rhs.Pstart)  return (Pstart >rhs.Pstart);  // fix this so that the set is sorted by decreasing number of primitives


        //everything else doesn't matter
        if (posA[0] !=rhs.posA[0])  return (posA[0] <rhs.posA[0]);
        if (posB[0] !=rhs.posB[0])  return (posB[0] <rhs.posB[0]);

        return false;
    }
};



//data and memory needed to evaluate one batch
struct XPint {

    unsigned char aa[maxK*maxK];
    unsigned char bb[maxK*maxK];

    cacheline64 AA[3];
    cacheline64 BB[3];
    cacheline64 AB[3];
    cacheline64 ABr2;

    CXPH Kernel;
    cacheline64 ** EEE;

    const eBasis * Abasis;
    const eBasis * Bbasis;

    int posA[DPC];
    int posB[DPC];
    int nUsed;

    int sEEE;

    int nK2;
    volatile int iK2;
    volatile int nThreads;


    XPint() {
        EEE = NULL;
        Abasis = Bbasis = NULL;
    }

    ~XPint() {
        for (int t=0; t<omp_get_max_threads(); ++t) delete[] EEE[t];
    }


    void Init(const XPtile & Unit, const map<int, aCenter> & Cpos) {

        Abasis = Unit.Abasis;
        Bbasis = Unit.Bbasis;

        int Ja = Abasis->J;
        int La = Abasis->L;
        int Jb = Bbasis->J;
        int Lb = Bbasis->L;
        const int Ltot = La+Lb;

        nUsed = Unit.nUsed;
        for (int d=0; d<DPC; ++d) posA[d] = Unit.posA[d];
        for (int d=0; d<DPC; ++d) posB[d] = Unit.posB[d];


        const int P3[] = {1, 4, 10, 20, 35, 56, 84, 120, 165, 220, 286}; //triangular piramidal numbers
        sEEE = Ja*Jb*P3[Ltot];

        EEE = new cacheline64*[omp_get_max_threads()];
        for (int t=0; t<omp_get_max_threads(); ++t) EEE[t] = NULL;


        //gather center information
        //-------------------------
        for (int d=0; d<DPC; ++d) {

            map<int, aCenter>::const_iterator ita, itb;

            ita = Cpos.find(Unit.posA[d]);
            itb = Cpos.find(Unit.posB[d]);

            //CHECK IF EVERYTHING OK!!!!!
            point A = ita->second.P;
            point B = itb->second.P;

            AA[0](d) = A.x;
            AA[1](d) = A.y;
            AA[2](d) = A.z;

            BB[0](d) = B.x;
            BB[1](d) = B.y;
            BB[2](d) = B.z;
        }


        AB[0] = BB[0]-AA[0];
        AB[1] = BB[1]-AA[1];
        AB[2] = BB[2]-AA[2];

        ABr2 = AB[0]*AB[0] + AB[1]*AB[1] + AB[2]*AB[2];


        const bool CMONO = true;
        const bool CDIP  = true;

        switch (Ltot) {
            case  0: Kernel = CalcExternalPotentialH <CMONO,CDIP, 0>; break;
            case  1: Kernel = CalcExternalPotentialH <CMONO,CDIP, 1>; break;
            case  2: Kernel = CalcExternalPotentialH <CMONO,CDIP, 2>; break;
            case  3: Kernel = CalcExternalPotentialH <CMONO,CDIP, 3>; break;
            case  4: Kernel = CalcExternalPotentialH <CMONO,CDIP, 4>; break;
            case  5: Kernel = CalcExternalPotentialH <CMONO,CDIP, 5>; break;
            case  6: Kernel = CalcExternalPotentialH <CMONO,CDIP, 6>; break;
            case  7: Kernel = CalcExternalPotentialH <CMONO,CDIP, 7>; break;
            case  8: Kernel = CalcExternalPotentialH <CMONO,CDIP, 8>; break;
            case  9: Kernel = CalcExternalPotentialH <CMONO,CDIP, 9>; break;
            case 10: Kernel = CalcExternalPotentialH <CMONO,CDIP,10>; break;
        }


        int ab=0;
        int ab0 = Unit.Pstart;
        int abf = Unit.Pstart+Unit.nP2;
        nK2 = 0;

        for (int b=0; b<Unit.nKi->nKb; ++b) {
            for (int a=0; a<Unit.nKi->nKa[b]; ++a) {
                if (ab>=ab0 && ab<abf) {
                    aa[nK2] = a;
                    bb[nK2] = b;
                    ++nK2;
                }

                ++ab;
            }
        }

        iK2 = 0;
        nThreads = 0;
    }

    //more threads can join and "help"
    void Compute(const Xpotential & Xpot) {
        int t = omp_get_thread_num();
        ++nThreads;

        EEE[t] = new cacheline64[sEEE];
        for (int i=0; i<DPC*sEEE; ++i) ((double*)EEE[t])[i] = 0.;

        while (iK2<nK2) {
            int ab;

            #pragma omp critical
            {
                ab = iK2;
                ++iK2;
            }

            int a = aa[ab];
            int b = bb[ab];
            Kernel (AA, Abasis->J, Abasis->k[a], Abasis->w[a],
                    BB, Bbasis->J, Bbasis->k[b], Bbasis->w[b],
                    Xpot,
                    ABr2, EEE[t]);
        }
    }

    void Transform(double * V) {

        int Ja = Abasis->J;
        int La = Abasis->L;
        int Jb = Bbasis->J;
        int Lb = Bbasis->L;
        const int Ltot = La+Lb;

        //sum if there are multiple contributions
        cacheline64 * ES = NULL;

        for (int t=0; t<omp_get_max_threads(); ++t) {
            if (EEE[t]!=NULL) {
                if (ES==NULL)
                    ES = EEE[t];
                //add other contributions
                else {
                    for (int i=0; i<sEEE; ++i)
                        ES[i] += EEE[t][i];
                }
            }
        }

        // transform to spherical
        for (int ja=0; ja<Ja; ++ja) {
            for (int jb=0; jb<Jb; ++jb) {

                //copy to a struct with better layout for RRs
                Triple64<2*LMAX> EE;

                int i3=0;
                for (int ex=0; ex<=Ltot; ++ex)
                    for (int ey=0; ey<=Ltot-ex; ++ey)
                        for (int ez=0; ez<=Ltot-ex-ey; ++ez,++i3) {
                            EE(ex,ey,ez) = ES[i3*Ja*Jb + ja*Jb + jb];
                        }


                //apply HRR
                cacheline64 VV[MMAX][MMAX];

                for (int ma=0; ma<nmC[La]; ++ma) {
                    int ax = CartList[La][ma].nx;
                    int ay = CartList[La][ma].ny;
                    int az = CartList[La][ma].nz;

                    for (int mb=0; mb<nmC[Lb]; ++mb) {
                        int bx = CartList[Lb][mb].nx;
                        int by = CartList[Lb][mb].ny;
                        int bz = CartList[Lb][mb].nz;

                        VV[ma][mb] = HRR(ax, ay, az,   bx, by, bz,  AB, EE);
                    }
                }


                //convert to spherical
                for (int a=0; a<nmS[La]; ++a) {
                    for (int b=0; b<nmS[Lb]; ++b) {

                        cacheline64 sumV;
                        sumV = 0;
                        for (int i=0; i<SHList[La][a].nps; ++i) {
                            for (int j=0; j<SHList[Lb][b].nps; ++j) {
                                const SHTerm & sha = SHList[La][a].T[i];
                                const SHTerm & shb = SHList[Lb][b].T[j];

                                sumV += VV[sha.nc][shb.nc] * (sha.cN * shb.cN);
                            }
                        }

                        //write to matrix
                        for (int d=0; d<nUsed; ++d) {
                            //compute address in triangular matrix
                            int i = posA[d] + ja*nmS[La] + a;
                            int j = posB[d] + jb*nmS[Lb] + b;

                            if (posA[d]==posB[d] && j>i) continue;

                            //compute address in triangular matrix
                            int ij;
                            if (i>j) ij = (i*i+i)/2 + j;
                            else     ij = (j*j+j)/2 + i;

                            V[ij] += sumV(d);
                        }

                    }
                }

            }
        }

    }

};

struct MPIprocess {
    priority_queue<XPtile> Work;
    int ID;

    void AddWork(const XPtile & rhs) {
        Work.push(rhs);
    }
};

//stores information about the elemental basis block
struct EBblock {
    eBasis  * Abasis;
    eBasis  * Bbasis;
    eBasis2 * ABbasis;
    aCenterDuple * pList;
    K2PS Kps;
    int maxt2;

    EBblock() {
        Abasis = Bbasis = NULL;
        ABbasis = NULL;
        pList = NULL;
        maxt2 = 0;
    }

    bool operator<(const EBblock & rhs) const {
        if (Abasis!=rhs.Abasis) return (Abasis<rhs.Abasis);
        if (Bbasis!=rhs.Bbasis) return (Bbasis<rhs.Bbasis);
        return false;
    }

    ~EBblock()  {
        //delete[] pList;
    }
};

struct XPevaluator {

    int MaxEB;
    int lenV;
    map <eBasis, deque<aCenter> > cLists; //generate map

    eBasis * Lbasis;
    eBasis2 ** L2basis;

    map<int, aCenter> Cpos; //centers for every position in the matrix

    set<EBblock> EBlists;


    XPevaluator() {
        Lbasis  = NULL;
        L2basis = NULL;
        MaxEB = 0;
        lenV  = 0;
    }

    void Free() {
        delete[] Lbasis;
        for (int i=0; i<MaxEB; ++i) delete[] L2basis[i];
        delete[] L2basis;

        set<EBblock>::const_iterator itEB;

        for (itEB=EBlists.begin(); itEB!=EBlists.end(); ++itEB) {
            const EBblock & Tblock = *itEB;
            delete[] Tblock.pList;
        }
    }

    ~XPevaluator() {
        Free();
    }

    // count elemental basis and classify functions
    void Init(int Kmax, const int * nhkt, const int * nuco, const int * nrco, const int * jstrt, const double * cent, const double * priccf, const double * priexp) {

        static const double RN[] = {3.54490770181103205433, 2.04665341589297697680, 0.91529123286376899328, 0.34594756847932032850, 0.11531585615977344279, 0.03476903884186927621};

        //make sure everything is empty
        Free();
        cLists.clear();
        Cpos.clear();
        EBlists.clear();

        //initialize

        lenV = 0;

        for (int shellA=0; shellA<Kmax; shellA+=nrco[shellA]) {

            eBasis basis; {
                basis.L = nhkt[shellA]-1;
                basis.J = nrco[shellA];
                basis.K = nuco[shellA];

                int primA = jstrt[shellA];

                for (int k=0; k<basis.K; ++k)
                    basis.k[k] = priexp[primA + k];

                for (int k=0; k<basis.K; ++k)
                    for (int j=0; j<basis.J; ++j)
                        basis.w[k][j] = priccf[j*MXPRIM + primA + k] * RN[basis.L];

                basis.sort();
            }

            aCenter center; {
                center.P.x = cent[shellA];
                center.P.y = cent[shellA +   MXSHEL];
                center.P.z = cent[shellA + 2*MXSHEL];
                center.n   = shellA;
                center.pos = lenV;
            }

            cLists[basis].push_back(center);

            lenV += basis.J * nmS[basis.L];
        }

        MaxEB = cLists.size();



        map <eBasis, deque<aCenter> >::const_iterator it1, it2;

        // make a hard copy of the elemental basis,
        // so that local functions and products can refer to it
        // ****************************************************

        int NEbasis = MaxEB;

        Lbasis = new eBasis[NEbasis];
        int neb = 0;
        for (it1=cLists.begin(); it1!=cLists.end(); ++it1, ++neb)
            Lbasis[neb] = it1->first;

        // create the product basis
        // ************************
        L2basis = new eBasis2*[NEbasis];

        for (int i=0; i<NEbasis; ++i) {
            L2basis[i] = new eBasis2[NEbasis];
            for (int j=0; j<NEbasis; ++j)
                L2basis[i][j].set(Lbasis[i], Lbasis[j]);
        }


        // make a list of position and centers
        // +++++++++++++++++++++++++++++++++++
        for (it1=cLists.begin(); it1!=cLists.end(); ++it1) {

            const deque<aCenter> & Alist = it1->second;
            deque<aCenter>::const_iterator q1;

            for (q1=Alist.begin(); q1!=Alist.end(); ++q1) {
                const aCenter & C = *q1;
                Cpos[C.pos] = C;
            }
        }


        //generate interaction lists
        //**************************

        //different functions (can be in same atom)
        int iit1, iit2;
        for (it1=cLists.begin(),iit1=0; it1!=cLists.end(); ++it1, ++iit1) {
            for (it2=cLists.begin(),iit2=0; it2!=it1; ++it2, ++iit2) {
                const eBasis & Abasis = it1->first;
                const eBasis & Bbasis = it2->first;
                const deque<aCenter> & Alist = it1->second;
                const deque<aCenter> & Blist = it2->second;

                eBasis2 & AB2basis = L2basis[iit1][iit2];

                int La = Abasis.L;
                int Ja = Abasis.J;
                int Ka = Abasis.K;

                int Lb = Bbasis.L;
                int Jb = Bbasis.J;
                int Kb = Bbasis.K;

                //specialized funciton for prescreening
                K2PS Kps; {

                    switch (La) {

                        case  0: {
                            switch(Lb) {
                                case  0: Kps = ABABps[0][0]; break;
                                case  1: Kps = ABABps[0][1]; break;
                                case  2: Kps = ABABps[0][2]; break;
                                case  3: Kps = ABABps[0][3]; break;
                                case  4: Kps = ABABps[0][4]; break;
                                case  5: Kps = ABABps[0][5]; break;
                            }
                            break;
                        }

                        case  1: {
                            switch(Lb) {
                                case  0: Kps = ABABps[1][0]; break;
                                case  1: Kps = ABABps[1][1]; break;
                                case  2: Kps = ABABps[1][2]; break;
                                case  3: Kps = ABABps[1][3]; break;
                                case  4: Kps = ABABps[1][4]; break;
                                case  5: Kps = ABABps[1][5]; break;
                            }
                            break;
                        }

                        case  2: {
                            switch(Lb) {
                                case  0: Kps = ABABps[2][0]; break;
                                case  1: Kps = ABABps[2][1]; break;
                                case  2: Kps = ABABps[2][2]; break;
                                case  3: Kps = ABABps[2][3]; break;
                                case  4: Kps = ABABps[2][4]; break;
                                case  5: Kps = ABABps[2][5]; break;
                            }
                            break;
                        }

                        case  3: {
                            switch(Lb) {
                                case  0: Kps = ABABps[3][0]; break;
                                case  1: Kps = ABABps[3][1]; break;
                                case  2: Kps = ABABps[3][2]; break;
                                case  3: Kps = ABABps[3][3]; break;
                                case  4: Kps = ABABps[3][4]; break;
                                case  5: Kps = ABABps[3][5]; break;
                            }
                            break;
                        }

                        case  4: {
                            switch(Lb) {
                                case  0: Kps = ABABps[4][0]; break;
                                case  1: Kps = ABABps[4][1]; break;
                                case  2: Kps = ABABps[4][2]; break;
                                case  3: Kps = ABABps[4][3]; break;
                                case  4: Kps = ABABps[4][4]; break;
                                case  5: Kps = ABABps[4][5]; break;
                            }
                            break;
                        }

                        case  5: {
                            switch(Lb) {
                                case  0: Kps = ABABps[5][0]; break;
                                case  1: Kps = ABABps[5][1]; break;
                                case  2: Kps = ABABps[5][2]; break;
                                case  3: Kps = ABABps[5][3]; break;
                                case  4: Kps = ABABps[5][4]; break;
                                case  5: Kps = ABABps[5][5]; break;
                            }
                            break;
                        }

                    }
                }


                // generate a list of all pairs and their distance
                // ===============================================

                int maxt2 = it1->second.size() * it2->second.size();
                aCenterDuple * pList = new aCenterDuple[maxt2];

                int imaxt2 = 0;

                deque<aCenter>::const_iterator q1, q2;

                for (q1=Alist.begin(); q1!=Alist.end(); ++q1) {
                    for (q2=Blist.begin(); q2!=Blist.end(); ++q2) {

                        vector3 rAB;
                        rAB.x = q2->P.x - q1->P.x;
                        rAB.y = q2->P.y - q1->P.y;
                        rAB.z = q2->P.z - q1->P.z;

                        double R2 = rAB.x*rAB.x + rAB.y*rAB.y + rAB.z*rAB.z;

                        //pList.insert()
                        pList[imaxt2].R2 =  R2;
                        pList[imaxt2].posA = q1->pos;
                        pList[imaxt2].posB = q2->pos;

                        ++imaxt2;
                    }
                }

                // sort the pairs by distance
                // ==========================
                FlashSort (pList, maxt2);

                EBblock Tblock;
                Tblock.Abasis  = &Lbasis[iit1];
                Tblock.Bbasis  = &Lbasis[iit2];
                Tblock.ABbasis = &L2basis[iit1][iit2];
                Tblock.pList   = pList;
                Tblock.maxt2   = maxt2;
                Tblock.Kps     = Kps;
                EBlists.insert(Tblock);
            }
        }

        //same elemental function
        for (it1=cLists.begin(), iit1=0; it1!=cLists.end(); ++it1, ++iit1) {

            const eBasis & Abasis = it1->first;
            const deque<aCenter> & Alist = it1->second;

            int La = Abasis.L;
            int Ja = Abasis.J;
            int Ka = Abasis.K;

            eBasis2 & AB2basis = L2basis[iit1][iit1];

            //specialized funciton for prescreening
            K2PS Kps; {
                switch (La) {
                    case  0: Kps = ABABps[0][0]; break;
                    case  1: Kps = ABABps[1][1]; break;
                    case  2: Kps = ABABps[2][2]; break;
                    case  3: Kps = ABABps[3][3]; break;
                    case  4: Kps = ABABps[4][4]; break;
                    case  5: Kps = ABABps[5][5]; break;
                }
            }

            // generate a list of all pairs and their distance
            // ===============================================
            int maxt2 = (it1->second.size() * it1->second.size() + it1->second.size())/2;
            aCenterDuple * pList = new aCenterDuple[maxt2];

            int imaxt2 = 0;

            deque<aCenter>::const_iterator q1, q2;

            for (q1=Alist.begin(); q1!=Alist.end(); ++q1) {
                for (q2=Alist.begin(); q2!=q1; ++q2) {

                    vector3 rAB;
                    rAB.x = q2->P.x - q1->P.x;
                    rAB.y = q2->P.y - q1->P.y;
                    rAB.z = q2->P.z - q1->P.z;

                    double R2 = rAB.x*rAB.x + rAB.y*rAB.y + rAB.z*rAB.z;

                    pList[imaxt2].R2 =  R2;
                    pList[imaxt2].posA = q1->pos;
                    pList[imaxt2].posB = q2->pos;

                    ++imaxt2;
                }

                //same basis
                {
                    pList[imaxt2].R2 =  0;
                    pList[imaxt2].posA = q1->pos;
                    pList[imaxt2].posB = q1->pos;

                    ++imaxt2;
                }
            }

            // sort the pairs by distance
            // ==========================
            //FlashSort (pList, maxt2,Kps,Abasis,Abasis,AB2basis,thr2V);
            FlashSort (pList, maxt2);


            EBblock Tblock;
            Tblock.Abasis  = &Lbasis[iit1];
            Tblock.Bbasis  = &Lbasis[iit1];
            Tblock.ABbasis = &L2basis[iit1][iit1];
            Tblock.pList   = pList;
            Tblock.maxt2   = maxt2;
            Tblock.Kps     = Kps;
            EBlists.insert(Tblock);
        }


    }


    void ComputeVintegrals(double * V, const Xpotential & Xpot) const {

        //initialize V
        for (int ij=0; ij<(lenV*lenV+lenV)/2; ++ij) V[ij] = 0;

        //skip intergal blocks with predicted norm below threshold
        const double threshold = 1e-15;
        const double thr2 = threshold*threshold;

        //compute the Cauchy-Schwarz parameter for the potential
        double Vps = CalcVps(Xpot);
        //cout << "Cauchy-Schwarz parameter for potential: " << Vps << endl;
        //cout << "Threads: " << omp_get_max_threads() << endl;

        const double thr2V = thr2/Vps;

        // keep track of total number of primitive
        // products of each angular momentum
        // ***************************************
        int np2l[2*LMAX+1];
        for (int i=0; i<=2*LMAX; ++i) np2l[i] = 0;


        // generate the tiles
        // ==================

        //worst case scenario, we have (N^2)/2DPC interacting elemental functions
        set<XPtile> setXPtiles[2*LMAX+1];

        set<EBblock>::const_iterator itEB;

        for (itEB=EBlists.begin(); itEB!=EBlists.end(); ++itEB) {
            const EBblock & Tblock = *itEB;

            const eBasis & Abasis    = *Tblock.Abasis;
            const eBasis & Bbasis    = *Tblock.Bbasis;
            const eBasis2 & AB2basis = *Tblock.ABbasis;

            const aCenterDuple * pList = Tblock.pList;
            int maxt2 = Tblock.maxt2;
            K2PS Kps  = Tblock.Kps;


            int La = Abasis.L;
            int Ja = Abasis.J;
            int Ka = Abasis.K;

            int Lb = Bbasis.L;
            int Jb = Bbasis.J;
            int Kb = Bbasis.K;


            //#pragma omp parallel for schedule(dynamic)
            for (int q=0; q<maxt2; q+=DPC) {
                int    nABk2 = Kps(Abasis, Bbasis, AB2basis, sqrt(pList[q].R2), thr2V);
                if (nABk2 == -1) break;

                np2l[La+Lb] += nABk2+1;

                XPtile Unit; {
                    Unit.Abasis = &Abasis;
                    Unit.Bbasis = &Bbasis;

                    Unit.nKi = &AB2basis.nKi[nABk2];
                    Unit.Ltot = La+Lb;
                    Unit.nP2  = nABk2+1;
                    Unit.nUsed = min(maxt2-q,DPC);
                    Unit.Pstart = 0;

                    for (int d=0; d<DPC; ++d) {
                        int p = q+d;
                        if (p>=maxt2) p=maxt2-1; //fill the eoa with copiesof the last element
                        Unit.posA[d] = pList[p].posA;
                        Unit.posB[d] = pList[p].posB;
                    }
                }

                #pragma omp critical
                setXPtiles[La+Lb].insert(Unit);
            }
        }


        // EVALUATE THE TILES !!!
        // **********************
/*
        static const int NMPI = 8; //simulate a number of MPI processes

        MPIprocess MPI_procs[NMPI];
        for (int id=0; id<NMPI; ++id) MPI_procs[id].ID = id;


        //first divide work across nodes
        //******************************
        for (int Lab=0; Lab<2*LMAX+1; ++Lab) {
            if (setXPtiles[Lab].size()==0) continue; //skip if empty

            //how many primitives must be evaluated for this angular momentum
            int nNP[NMPI];

            int nNodePrim = np2l[Lab]/NMPI;
            for (int i=0; i<NMPI;           ++i)   nNP[i] = nNodePrim;
            for (int i=0; i<np2l[Lab]%NMPI; ++i) ++nNP[i];

            //algorithm 1: tends to cluster many integrals with low K2 on the nodes with higher ID
            if (0) {

                set<XPtile>::const_iterator XPit, XPend;
                XPit  = setXPtiles[Lab].begin();
                XPend = setXPtiles[Lab].end();


                int ID_MPI = 0;

                //allocate the primitives across nodes
                for (;XPit!=XPend; ++XPit) {
                    XPtile Unit;
                    Unit = *XPit;

                    while (Unit.nP2>0) {
                        if (nNP[ID_MPI] > Unit.nP2) {
                            MPI_procs[ID_MPI].AddWork(Unit);
                            nNP[ID_MPI]-=Unit.nP2;
                            Unit.nP2=0;
                        }
                        else {
                            int AP2 = nNP[ID_MPI];
                            int nP2 = Unit.nP2-AP2;
                            Unit.nP2 = AP2;
                            MPI_procs[ID_MPI].AddWork(Unit);

                            Unit.Pstart += AP2;
                            Unit.nP2     = nP2;

                            nNP[ID_MPI]-=AP2;
                            ++ID_MPI;
                        }
                    }
                }
            }

            //algorithm 2: looks much better balanced;
            //also it might be possible to avoid internode communication
            else {
                set<XPtile>::iterator XPit;

                int ID_MPI = 0;

                //allocate the primitives across nodes
                while (setXPtiles[Lab].size()>0) {
                    //XPit  = setXPtiles[Lab].end();
                    //--XPit;
                    XPit = setXPtiles[Lab].begin();
                    XPtile Unit;
                    Unit = *XPit;
                    setXPtiles[Lab].erase(XPit);


                    while (nNP[ID_MPI]==0) {
                        ++ID_MPI;
                        ID_MPI %= NMPI;
                    }


                    if (nNP[ID_MPI] >= Unit.nP2) {
                        MPI_procs[ID_MPI].AddWork(Unit);
                        nNP[ID_MPI]-=Unit.nP2;
                    }
                    else {
                        int AP2 = nNP[ID_MPI];
                        int nP2 = Unit.nP2-AP2;
                        Unit.nP2 = AP2;
                        MPI_procs[ID_MPI].AddWork(Unit);

                        //
                        Unit.Pstart += AP2;
                        Unit.nP2     = nP2;
                        setXPtiles[Lab].insert(Unit);

                        nNP[ID_MPI]-=AP2;
                    }


                    ++ID_MPI;
                    ID_MPI %= NMPI;
                }
            }
        }

*/


        int NMPI   = mpi_nnodes();
        int MYRANK = mpi_myrank();

        MPIprocess NodeWork;


        //first divide work across nodes
        //******************************
        for (int Lab=0; Lab<2*LMAX+1; ++Lab) {
            if (setXPtiles[Lab].size()==0) continue; //skip if empty

            //how many primitives must be evaluated for this angular momentum
            int nNP[NMPI];

            int nNodePrim = np2l[Lab]/NMPI;
            for (int i=0; i<NMPI;           ++i)   nNP[i] = nNodePrim;
            for (int i=0; i<np2l[Lab]%NMPI; ++i) ++nNP[i];


            //algorithm 2: looks much better balanced;
            //also it might be possible to avoid internode communication
            {
                set<XPtile>::iterator XPit;

                int ID_MPI = 0;

                //allocate the primitives across nodes
                while (setXPtiles[Lab].size()>0) {
                    //XPit  = setXPtiles[Lab].end();
                    //--XPit;
                    XPit = setXPtiles[Lab].begin();
                    XPtile Unit;
                    Unit = *XPit;
                    setXPtiles[Lab].erase(XPit);


                    while (nNP[ID_MPI]==0) {
                        ++ID_MPI;
                        ID_MPI %= NMPI;
                    }


                    if (nNP[ID_MPI] >= Unit.nP2) {
                        if (MYRANK==ID_MPI)
                            NodeWork.AddWork(Unit);
                        nNP[ID_MPI]-=Unit.nP2;
                    }
                    else {
                        int AP2 = nNP[ID_MPI];
                        int nP2 = Unit.nP2-AP2;
                        Unit.nP2 = AP2;

                        if (MYRANK==ID_MPI)
                            NodeWork.AddWork(Unit);

                        //
                        Unit.Pstart += AP2;
                        Unit.nP2     = nP2;
                        setXPtiles[Lab].insert(Unit);

                        nNP[ID_MPI]-=AP2;
                    }


                    ++ID_MPI;
                    ID_MPI %= NMPI;
                }
            }
        }




        //execute every node
        //******************
        //for (int ID_MPI=0; ID_MPI<NMPI; ++ID_MPI)
        {

            //priority_queue<XPtile> & qNode = MPI_procs[ID_MPI].Work;
            //cout << "MPI process " << ID_MPI << endl;

            priority_queue<XPtile> & qNode = NodeWork.Work;


            set<XPint*> qReady;

            const int ogmt = omp_get_max_threads();

            XPint * volatile XPTlist[ogmt];
            for (int t=0; t<ogmt; ++t) XPTlist[t] = NULL;


            #pragma omp parallel
            while(!qNode.empty()) {

                XPtile Unit;

                bool done = false;

                #pragma omp critical
                {
                    if (qNode.empty()) {
                        done = true;
                    }
                    else {
                        Unit = qNode.top();
                        qNode.pop();

                        //cout << omp_get_thread_num() << " initializing block; " << Unit.Ltot << " " << Unit.nP2 << endl; cout.flush();
                    }
                }


                if (done) {

                    XPint  * XPhelp = NULL;

                    do {

                        //should try first with computations being executed only by one thread, then two, and so on
                        //array<set<*>>

                        //#pragma omp critical
                        {
                        /*
                            set<XPint*>::iterator qit;
                            set<XPint*>::iterator qmax = qReady.end();
                            float pmax = 0;

                            for (qit = qReady.begin(); qit!=qReady.end(); ++qit) {
                                float p = ((*qit)->nK2 - (*qit)->iK2) / (*qit)->nThreads;
                                if (p>pmax) {
                                    pmax=p;
                                    qmax = qit;
                                }
                            }
                            */

                            for (int t=0; t<ogmt; ++t)
                                if (XPTlist[t]!=NULL) XPhelp = XPTlist[t];


                            //if (qmax!=qReady.end()) XPhelp = (*qmax);

                            //if (qReady.begin()!=qReady.end()) XPhelp = (*qReady.begin());

                            cout << "thread " << omp_get_thread_num() << " helping" << endl;
                        }

                        if (XPhelp!=NULL) XPhelp->Compute   (Xpot);
                    } while (XPhelp!=NULL);

                    continue;
                }


                //create and initialize object
                XPint  * XPlocal = new XPint;
                XPlocal->Init      (Unit, Cpos);

                //ready to receive help
                /*
                #pragma omp critical
                qReady.insert(XPlocal);
                */

                //XPTlist[omp_get_thread_num()] = XPlocal;

                //compute all integrals
                XPlocal->Compute   (Xpot);

                //ready to receive help
                /*
                #pragma omp critical
                qReady.erase(XPlocal);
                */

                //XPTlist[omp_get_thread_num()] = NULL;

                XPlocal->Transform (V);
                delete XPlocal;
            }

        }

    }


};



//headers for benchmarking
#include <fstream>
#include <sys/time.h>

extern "C" {
    static bool InitLibs = false; //don't export this symbol

    XPevaluator XPeval;


    void init_vpot(int & Kmax, int * nhkt, int * nuco, int * nrco, int * jstrt,  double * cent, double * priccf, double * priexp) {

        if (!InitLibs) {
            //init solid spherical harmonic polynomial expansions
            LibAngularVPD::InitCartList();
            LibAngularVPD::InitSHList();

            //compute coefficients of the incomplete gamma for the evaluation by taylor expqansion
            LibIGammaVPD::InitIncompleteGammas();

            //initialize primitive density-pair prescreener
            ABABinit();

            InitLibs = true;
        }

        mpi_bcast(Kmax);

        Bcast_Basis (nhkt, nuco, nrco, jstrt, cent, priccf, priexp);
        XPeval.Init(Kmax, nhkt, nuco, nrco, jstrt, cent, priccf, priexp);
    }

    void init_vpotdamped_slave_() {

        int Kmax;

        //initialize
        //==========
        int * nhkt = new int[MXSHEL];
        int * nuco = new int[MXSHEL];
        int * nrco = new int[MXSHEL];
        int * jstrt = new int[MXSHEL];

        double * cent   = new double[6*MXSHEL];
        double * priccf = new double[MXPRIM*MXCONT];
        double * priexp = new double[MXPRIM];

        init_vpot(Kmax, nhkt, nuco, nrco, jstrt, cent, priccf, priexp);

        //clean the mess
        //--------------
        delete[] nhkt;
        delete[] nuco;
        delete[] nrco;
        delete[] jstrt;

        delete[] cent;
        delete[] priccf;
        delete[] priexp;
    }

    void init_vpotdamped_(const int * kmax, int * nhkt, int * nuco, int * nrco, int * jstrt, double * cent, double * priccf, double * priexp) {

        int Kmax     = *kmax;

        init_vpot(Kmax, nhkt, nuco, nrco, jstrt, cent, priccf, priexp);
    }


    void ipotdamped_slave_() {


        // set to zero to remove warning about uninitialized Ncenters below
        // is anyway overwritten by mpi_bcast
        int Ncenters = 0;

        mpi_bcast(Ncenters);

        point   * NN = new point  [Ncenters]; //centers
        vector3 * DD = new vector3[Ncenters]; //dipole moments

        double * QQ = new double[Ncenters];
        double * Rq = new double[Ncenters];
        double * Rd = new double[Ncenters];

        Bcast_Potential(NN, QQ, Rq, DD, Rd, Ncenters);

        Xpotential Xpot; {
            Xpot.centers  = NN;
            Xpot.Q        = QQ;
            Xpot.Rq       = Rq;
            Xpot.D        = DD;
            Xpot.Rd       = Rd;
            Xpot.ncenters = Ncenters;
        }


        //compute the integrals
        //*********************
        int lenV2 = (XPeval.lenV*XPeval.lenV+XPeval.lenV)/2;
        double * VV = new double[lenV2];

        XPeval.ComputeVintegrals (VV, Xpot);

        ReduceV(NULL, VV, lenV2);


        //clean the mess
        //--------------
        delete[] NN;
        delete[] DD;
        delete[] QQ;
        delete[] Rq;
        delete[] Rd;

        delete[] VV;
    }

    void ipotdamped_(double * V,
                     const double * centers,
                     double * QQ,  double * Rq,
                     double *  D,  double * Rd,
                     const int * ncenters) {   // external potential parameters


        //scatter the potential data across nodes
        //***************************************
        int Ncenters = *ncenters;
        mpi_bcast(Ncenters);

        point   * NN = new point  [Ncenters]; //centers
        vector3 * DD = new vector3[Ncenters]; //dipole moments

        //copy coordinates from SoA to AoS
        for (int i=0; i<Ncenters; ++i) {
            NN[i].x = centers[3*i];
            NN[i].y = centers[3*i + 1];
            NN[i].z = centers[3*i + 2];
        }

        //copy dipoles from SoA to AoS
        for (int i=0; i<Ncenters; ++i) {
            DD[i].x = D[3*i];
            DD[i].y = D[3*i + 1];
            DD[i].z = D[3*i + 2];
        }

        Bcast_Potential(NN, QQ, Rq, DD, Rd, Ncenters);

        Xpotential Xpot; {
            Xpot.centers  = NN;
            Xpot.Q        = QQ;
            Xpot.Rq       = Rq;
            Xpot.D        = DD;
            Xpot.Rd       = Rd;
            Xpot.ncenters = Ncenters;
        }


        const int USECSPERSEC = 1000000;
        timeval timer[2];

        gettimeofday(&timer[0], NULL);

        int lenV2 = (XPeval.lenV*XPeval.lenV+XPeval.lenV)/2;
        double * VV = new double[lenV2];

        XPeval.ComputeVintegrals (VV, Xpot);

        ReduceV(V, VV, lenV2);

        gettimeofday(&timer[1], NULL);

        int secs (timer[1].tv_sec  - timer[0].tv_sec);
        int usecs(timer[1].tv_usec - timer[0].tv_usec);

        double At = double(secs) + double(usecs)/double(USECSPERSEC);

        //cout << "Time evaluating QMCMM integrals: " << At << endl;

        delete[] NN;
        delete[] DD;
        delete[] VV;
    }


    void vpotdamped_slave_() {

        init_vpotdamped_slave_();
        ipotdamped_slave_();
    }

    // O(N*M)
    // ipotdamped
    void vpotdamped_(const int * kmax, int * nhkt, int * nuco, int * nrco, int * jstrt, double * cent, double * priccf, double * priexp,        // shell parameters
                 double * V,                                                                                                                // output (potential matrix)
                 const double * centers,
                 double * QQ,  double * Rq,
                 double *  D,  double * Rd,
                 const int * ncenters) {   // external potential parameters

        init_vpotdamped_ (kmax, nhkt, nuco, nrco, jstrt, cent, priccf, priexp);
        ipotdamped_      (V,   centers, QQ, Rq, D, Rd, ncenters);
    }

}

