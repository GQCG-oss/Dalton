#include <cmath>
#include <cstdlib>
#include <iostream>
#include "ababps.hpp"
using namespace std;
using namespace LibAngularVPD;
using namespace LibIGammaVPD;

K2PS ABABps[6][6];

eBasis::eBasis() {
    for (int i=0; i<maxK; ++i)
        k[i] = 0;
    for (int i=0; i<maxK; ++i)
        for (int j=0; j<maxJ; ++j)
            w[i][j] = 0;
    L = J = K = 0;
}

 //sort the primtives in ascending order
void eBasis::sort() {

    for (int i=0; i<K; ++i) {
        double kmin = k[i];
        int    imin = i;

        for (int j=i+1; j<K; ++j) {
            if (k[j] < kmin) {
                imin = j;
                kmin = k[j];
            }
        }

        swap(k[i], k[imin]);

        for (int j=0; j<J; ++j)
            swap(w[i][j], w[imin][j]);
    }

}

//for sorting
bool eBasis::operator<(const eBasis & rhs) const {
    if (L!=rhs.L) return (L<rhs.L);
    if (J!=rhs.J) return (J<rhs.J);
    if (K!=rhs.K) return (K<rhs.K);

    for (int i=0; i<K; ++i)
        if (k[i]!=rhs.k[i]) return (k[i]<rhs.k[i]);

    for (int i=0; i<K; ++i)
        for (int j=0; j<J; ++j)
            if (w[i][j]!=rhs.w[i][j]) return (w[i][j]<rhs.w[i][j]);

    return false;
}


void eBasis2::set(const eBasis & Abasis, const eBasis & Bbasis) {
    La = Abasis.L;
    Lb = Bbasis.L;
    Ka = Abasis.K;
    Kb = Bbasis.K;
    Ja = Abasis.J;
    Jb = Bbasis.J;

    K2 = Ka*Kb;

    nKi = new nKinteracting[Ka*Kb+1];

    int nKt=0;

    for (int p=0; p<Ka; ++p) {
        for (int q=0; q<Kb; ++q) {
            nKi[nKt].ikappa = 1/Abasis.k[p] + 1/Bbasis.k[q];
            nKi[nKt].p  = p;
            nKi[nKt].q  = q;
            ++nKt;
        }
    }

    //sort them
    for (int i=0; i<nKt; ++i) {
        double kmax = nKi[i].ikappa;
        int    imax = i;

        for (int j=i+1; j<nKt; ++j) {
            if (nKi[j].ikappa > kmax) {
                imax = j;
                kmax = nKi[j].ikappa;
            }
        }

        swap(nKi[i], nKi[imax]);
    }


    int nKa[maxK];

    int nKb = 0;
    for (int p=0; p<maxK; ++p) nKa[p] = 0;


    for (int i=0; i<nKt; ++i) {
        int ka = nKi[i].p;
        int kb = nKi[i].q;

        if (nKa[kb] == 0) ++nKb;
        ++nKa[kb];

        //Psets[i].nK2 = i+1;
        //Psets[i].maxD2 = pps[i].ik;

        for (int p=0; p<maxK; ++p) nKi[i].nKa[p] = nKa[p];
        nKi[i].nKb   = nKb;
    }

    /*
    //print
    for (int i=0; i<nKt; ++i) {
        cout << "product " << i << endl;
        cout << nKi[i].p << " " << nKi[i].q << "   " << nKi[i].ikappa << endl;
        cout << "Maxb" << int(nKi[i].nKb) << endl;
        for (int p=0; p<maxK; ++p) cout << int(nKi[i].nKa[p]) << " "; cout << endl;
        cout << endl;
    }
    */


}



typedef double TripleXX [4*LMAX+1][4*LMAX+1][4*LMAX+1];

typedef void (*RintXX) (TripleXX & , double, double, double);

template <int L> void RIntegralXX(TripleXX & Rtuv, double k1, double k2, const double R) {

    double Fn[L+1];

    double R2 = R*R;

    double ** gamma = IncompleteGamma.gamma_table;
    double *  gamma2 = IncompleteGammas[L].gamma_table_vals;

    double k = (k1*k2)/(k1+k2);

    //calcula las 'F's
    {
        double Sn[L+1];

        double kR2 = k*R2;

        if (kR2>LibIGammaVPD::vg_max-LibIGammaVPD::vg_step) {
            double iR2 = 1/R2;
            Sn[0] = 0.5*sqrt(PI*iR2);

            for (int i=1; i<=L; ++i) {
                Sn[i] = Sn[i-1]*iR2*double(2*i-1);
            }
        }
        else {
            double p = LibIGammaVPD::ivg_step*(kR2-LibIGammaVPD::vg_min);
            int pos = int(p+0.5);
            double x0 = LibIGammaVPD::vg_min + LibIGammaVPD::vg_step*double(pos);
            double Ax1 = x0-kR2;
            double Ax2 = 0.5 * Ax1*Ax1;
            double Ax3 = 0.33333333333333333333 * Ax1*Ax2;

            Sn[L] = (gamma[pos][L+1] + gamma[pos][L+2] * Ax1 + gamma[pos][L+3] * Ax2 + gamma[pos][L+4] * Ax3);

            if (L>0) {
                double expx = gamma[pos][0] * (1+Ax1+Ax2+Ax3);

                for (int i=L-1; i>=0; --i)
                    Sn[i] = (2*kR2*Sn[i+1] + expx)*(1./double(2*i+1));
            }

            double kn   = sqrt(k);

            for (int n=0; n<=L; ++n) {
                Sn[n] *= kn;
                kn    *= 2*k;
            }
        }

        for (int n=0; n<=L; ++n) Fn[n] = Sn[n];
    }


    //calcula R^{n}_{000}
    //double ac = (2*PI)/(k*sqrt(k));
    double ac = 2*PI*PI*sqrt(PI) / (k1*k2*sqrt(k1*k2));

    for (int n=0; n<=L; ++n) Fn[n] *= ac;

    //calcula R_{tuv} reaprovechando la memoria intermedia
    //(podria tener un impacto negativo debido a la dependencia de variables o positivo por usar menos registros/cache)
    Rtuv[0][0][0] = Fn[L];
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
                            else if (v==1) Rtuv[t][u][v] = R * Rtuv[t][u][v-1];
                            else           Rtuv[t][u][v] = R * Rtuv[t][u][v-1] -(v-1)*Rtuv[t][u][v-2];
                        }
                        else if (u==1) Rtuv[t][u][v] = 0;
                        else           Rtuv[t][u][v] = -(u-1)*Rtuv[t][u-2][v];
                    }
                    else if (t==1) Rtuv[t][u][v] = 0;
                    else           Rtuv[t][u][v] = -(t-1)*Rtuv[t-2][u][v];
                }
            }
        }

        //tuv = 0
        Rtuv[0][0][0] = Fn[n];
    }

}


//for screening pairs
template <int La, int Lb> double CalcGTO2ps(const eBasis & A, const eBasis & B, double R) {

    const int Lab = La+Lb;

    RintXX V1Cint; {
        switch (Lab) {
            case  0: V1Cint = RIntegralXX< 0>; break;
            case  1: V1Cint = RIntegralXX< 2>; break;
            case  2: V1Cint = RIntegralXX< 4>; break;
            case  3: V1Cint = RIntegralXX< 6>; break;
            case  4: V1Cint = RIntegralXX< 8>; break;
            case  5: V1Cint = RIntegralXX<10>; break;
            case  6: V1Cint = RIntegralXX<12>; break;
            case  7: V1Cint = RIntegralXX<14>; break;
            case  8: V1Cint = RIntegralXX<16>; break;
            case  9: V1Cint = RIntegralXX<18>; break;
            case 10: V1Cint = RIntegralXX<20>; break;
        }
    }


    double ps = 0;

    for (int a1=0; a1<A.K; ++a1) {
        for (int b1=0; b1<B.K; ++b1) {
            for (int a2=0; a2<A.K; ++a2) {
                for (int b2=0; b2<B.K; ++b2) {

                    double ka = A.k[a1];
                    double kb = B.k[b1];
                    double kc = A.k[a2];
                    double kd = B.k[b2];

                    double Kab = ka*kb/(ka+kb);
                    double Kcd = kc*kd/(kc+kd);

                    double kab = ka+kb;
                    double kcd = kc+kd;


                    //add up contributions from all GC functions
                    double c2a = 0;

                    for (int ja=0; ja<A.J; ++ja)
                        c2a += A.w[a1][ja] * A.w[a2][ja];

                    double c2b = 0;

                    for (int jb=0; jb<B.J; ++jb)
                        c2b += B.w[b1][jb] * B.w[b2][jb];


                    double c2 = c2a*c2b * exp(-(Kab+Kcd)*R*R);


                    //compute the trace of (ab|ab)
                    double tr = 0;

                    {
                        // A is in 0, B is in R
                        double P =  R * kb/(ka+kb);
                        double Q =  R * kd/(kc+kd);


                        //expand x^n in terms of hermite functions
                        double HexpAB[Lab+1][Lab+1]; {

                            for (int a=0; a<=Lab; ++a)
                                for (int p=0; p<=Lab; ++p)
                                    HexpAB[a][p] = 0;

                            HexpAB[0][0] = 1;
                            HexpAB[1][1] = -(2*kab);

                            for (int i=2; i<=Lab; ++i) {
                                for (int p=1; p<=i; ++p)
                                    HexpAB[i][p] = -2*kab*HexpAB[i-1][p-1];
                                for (int p=0; p<i-1; ++p)
                                    HexpAB[i][p] += (p+1) * HexpAB[i-1][p+1];
                            }
                        }

                        double HexpCD[Lab+1][Lab+1]; {

                            for (int a=0; a<=Lab; ++a)
                                for (int p=0; p<=Lab; ++p)
                                    HexpCD[a][p] = 0;

                            HexpCD[0][0] = 1;
                            HexpCD[1][1] = -(2*kcd);

                            for (int i=2; i<=Lab; ++i) {
                                for (int p=1; p<=i; ++p)
                                    HexpCD[i][p] = -2*kcd*HexpCD[i-1][p-1];
                                for (int p=0; p<i-1; ++p)
                                    HexpCD[i][p] += (p+1) * HexpCD[i-1][p+1];
                            }
                        }


                        double iHexpAB[Lab+1][Lab+1]; {

                            for (int a=0; a<=Lab; ++a)
                                for (int p=0; p<=Lab; ++p)
                                    iHexpAB[a][p] = 0;

                            for (int a=0; a<=Lab; ++a) {
                                double hpol[Lab+1];
                                for (int p=0; p<=Lab; ++p) hpol[p] = 0; hpol[a] = 1;

                                for (int p=a; p>=0; --p) {
                                    double r = hpol[p]/HexpAB[p][p];

                                    iHexpAB[a][p] = r;

                                    for (int i=0; i<=p; ++i) hpol[i] -= r*HexpAB[p][i];
                                }
                            }
                        }

                        double iHexpCD[Lab+1][Lab+1]; {

                            for (int a=0; a<=Lab; ++a)
                                for (int p=0; p<=Lab; ++p)
                                    iHexpCD[a][p] = 0;

                            for (int a=0; a<=Lab; ++a) {
                                double hpol[Lab+1];
                                for (int p=0; p<=Lab; ++p) hpol[p] = 0; hpol[a] = 1;

                                for (int p=a; p>=0; --p) {
                                    double r = hpol[p]/HexpCD[p][p];

                                    iHexpCD[a][p] = r;

                                    for (int i=0; i<=p; ++i) hpol[i] -= r*HexpCD[p][i];
                                }
                            }
                        }


                        double PexpAB[La+1][Lb+1][2*Lab+1]; {

                            for (int a=0; a<=La; ++a)
                                for (int b=0; b<=Lb; ++b)
                                    for (int p=0; p<=Lab; ++p)
                                        PexpAB[a][b][p] = 0;

                            PexpAB[0][0][0] = 1;
                            for (int a=1; a<=La; ++a) {
                                for (int p=0; p<a; ++p) PexpAB[a][0][p]    = -P * PexpAB[a-1][0][p];
                                for (int p=0; p<a; ++p) PexpAB[a][0][p+1] +=     PexpAB[a-1][0][p];

                            }

                            for (int b=1; b<=Lb; ++b) {
                                for (int a=0; a<=La; ++a) {
                                    for (int p=0; p<a+b; ++p) PexpAB[a][b][p]    = (R-P) * PexpAB[a][b-1][p];
                                    for (int p=0; p<a+b; ++p) PexpAB[a][b][p+1] +=         PexpAB[a][b-1][p];
                                }
                            }
                        }

                        double PexpCD[La+1][Lb+1][2*Lab+1]; {

                            for (int a=0; a<=La; ++a)
                                for (int b=0; b<=Lb; ++b)
                                    for (int p=0; p<=Lab; ++p)
                                        PexpCD[a][b][p] = 0;

                            PexpCD[0][0][0] = 1;
                            for (int a=1; a<=La; ++a) {
                                for (int p=0; p<a; ++p) PexpCD[a][0][p]    = -Q * PexpCD[a-1][0][p];
                                for (int p=0; p<a; ++p) PexpCD[a][0][p+1] +=     PexpCD[a-1][0][p];

                            }

                            for (int b=1; b<=Lb; ++b) {
                                for (int a=0; a<=La; ++a) {
                                    for (int p=0; p<a+b; ++p) PexpCD[a][b][p]    = (R-Q) * PexpCD[a][b-1][p];
                                    for (int p=0; p<a+b; ++p) PexpCD[a][b][p+1] +=         PexpCD[a][b-1][p];
                                }
                            }
                        }


                        double YYexpAB[2*La+1][2*Lb+1][Lab+1][Lab+1][Lab+1]; {
                            //expand
                            for (int ma=0; ma<2*La+1; ++ma) {
                                for (int mb=0; mb<2*Lb+1; ++mb) {
                                    for (int xx=0; xx<=Lab; ++xx)
                                    for (int yy=0; yy<=Lab; ++yy)
                                    for (int zz=0; zz<=Lab; ++zz)
                                        YYexpAB[ma][mb][xx][yy][zz] = 0;

                                    for (int i=0; i<SHList[La][ma].nps; ++i) {
                                        for (int j=0; j<SHList[Lb][mb].nps; ++j) {
                                            const SHTerm & sha = SHList[La][ma].T[i];
                                            const SHTerm & shb = SHList[Lb][mb].T[j];
                                            int xx = sha.nx + shb.nx;
                                            int yy = sha.ny + shb.ny;
                                            int zz = sha.nz + shb.nz;

                                            for (int zc=0; zc<=zz; ++zc)
                                                YYexpAB[ma][mb][xx][yy][zc] += sha.cN * shb.cN * PexpAB[sha.nz][shb.nz][zc];
                                        }
                                    }

                                }
                            }
                        }

                        double YYexpCD[2*La+1][2*Lb+1][Lab+1][Lab+1][Lab+1]; {
                            //expand
                            for (int ma=0; ma<2*La+1; ++ma) {
                                for (int mb=0; mb<2*Lb+1; ++mb) {
                                    for (int xx=0; xx<=Lab; ++xx)
                                    for (int yy=0; yy<=Lab; ++yy)
                                    for (int zz=0; zz<=Lab; ++zz)
                                        YYexpCD[ma][mb][xx][yy][zz] = 0;

                                    for (int i=0; i<SHList[La][ma].nps; ++i) {
                                        for (int j=0; j<SHList[Lb][mb].nps; ++j) {
                                            const SHTerm & sha = SHList[La][ma].T[i];
                                            const SHTerm & shb = SHList[Lb][mb].T[j];
                                            int xx = sha.nx + shb.nx;
                                            int yy = sha.ny + shb.ny;
                                            int zz = sha.nz + shb.nz;

                                            for (int zc=0; zc<=zz; ++zc)
                                                YYexpCD[ma][mb][xx][yy][zc] += sha.cN * shb.cN * PexpCD[sha.nz][shb.nz][zc];
                                        }
                                    }

                                }
                            }
                        }

                        /*
                        double T6[2*LMAX+1][2*LMAX+1][2*LMAX+1][2*LMAX+1][2*LMAX+1][2*LMAX+1]; {
                            int Lab = Lab;

                            for (int a=0; a<=Lab; ++a)
                                for (int b=0; b<=Lab; ++b)
                                    for (int c=0; c<=Lab; ++c)
                                        for (int p=0; p<=Lab; ++p)
                                            for (int q=0; q<=Lab; ++q)
                                                for (int r=0; r<=Lab; ++r)
                                                    T6[a][b][c][p][q][r] = 0;


                            for (int ma=0; ma<2*La+1; ++ma) {
                                for (int mb=0; mb<2*Lb+1; ++mb) {

                                    for (int a=0; a<=Lab; ++a)
                                        for (int b=0; b<=Lab; ++b)
                                            for (int c=0; c<=Lab; ++c)
                                                for (int p=0; p<=Lab; ++p)
                                                    for (int q=0; q<=Lab; ++q)
                                                        for (int r=0; r<=Lab; ++r)
                                                            T6[a][b][c][p][q][r] += YYexpAB[ma][mb][a][b][c] * YYexpCD[ma][mb][p][q][r];
                                }
                            }
                        }
                        */


                        double T12r[Lab+1][Lab+1][2*Lab+1]; {

                            for (int a=0; a<=Lab; ++a)
                                for (int b=0; b<=Lab; ++b) {

                                    for (int r=0; r<=2*Lab; ++r)
                                        T12r[a][b][r] = 0;

                                    for (int p=0; p<=Lab; ++p)
                                        for (int q=0; q<=Lab; ++q)
                                            if (q%2==0) T12r[a][b][p+q] += iHexpAB[a][p] * iHexpCD[b][q];
                                            else        T12r[a][b][p+q] -= iHexpAB[a][p] * iHexpCD[b][q];
                                }
                        }

                        /*

                        double R3[4*LMAX+1][4*LMAX+1][4*LMAX+1]; {

                            for (int x=0; x<=2*Lab; ++x)
                                for (int y=0; y<=2*Lab; ++y)
                                    for (int z=0; z<=2*Lab; ++z)
                                        R3[x][y][z] = 0;

                            for (int a=0; a<=Lab; ++a)
                                for (int b=0; b<=Lab; ++b)
                                    for (int c=0; c<=Lab; ++c)
                                        for (int p=0; p<=Lab; ++p)
                                            for (int q=0; q<=Lab; ++q)
                                                for (int r=0; r<=Lab; ++r)
                                                    for (int x=0; x<=2*Lab; ++x)
                                                        for (int y=0; y<=2*Lab; ++y)
                                                            for (int z=0; z<=2*Lab; ++z)
                                                                R3[x][y][z] += T12r[a][p][x]*T12r[b][q][y]*T12r[c][r][z] *  T6[a][b][c][p][q][r];

                        }
                        */

                        double R3[2*Lab+1][2*Lab+1][2*Lab+1]; {

                            for (int x=0; x<=2*Lab; ++x)
                                for (int y=0; y<=2*Lab-x; ++y)
                                    for (int z=0; z<=2*Lab-x-y; ++z)
                                        R3[x][y][z] = 0;

                            for (int ma=0; ma<2*La+1; ++ma)
                                for (int mb=0; mb<2*Lb+1; ++mb)
                            for (int a=0; a<=Lab; ++a)
                                for (int b=0; b<=Lab-a; ++b)
                                    for (int c=0; c<=Lab-a-b; ++c)
                                        for (int p=0; p<=Lab; ++p)
                                            for (int q=0; q<=Lab-p; ++q)
                                                for (int r=0; r<=Lab-p-q; ++r)
                                                    for (int x=0; x<=2*Lab; ++x)
                                                        for (int y=0; y<=2*Lab-x; ++y)
                                                            for (int z=0; z<=2*Lab-x-y; ++z)
                                                                R3[x][y][z] += T12r[a][p][x]*T12r[b][q][y]*T12r[c][r][z] *  YYexpAB[ma][mb][a][b][c] * YYexpCD[ma][mb][p][q][r];
                        }


                        //compute the 1C integral
                        TripleXX TXX;
                        V1Cint(TXX, kab, kcd, Q-P);

                        for (int x=0; x<=2*Lab; ++x)
                            for (int y=0; y<=2*Lab-x; ++y)
                                for (int z=0; z<=2*Lab-x-y; ++z)
                                    tr += R3[x][y][z] * TXX[x][y][z];

                    }

                    ps += c2 * tr;
                }
            }
        }
    }


    return ps;
}



//for screening pairs
template <int La, int Lb> int CalcK2ps(const eBasis & A, const eBasis & B, const eBasis2 & AB, double R, double thresh) {

    const int Lab = La+Lb;

    RintXX V1Cint; {
        switch (Lab) {
            case  0: V1Cint = RIntegralXX< 0>; break;
            case  1: V1Cint = RIntegralXX< 2>; break;
            case  2: V1Cint = RIntegralXX< 4>; break;
            case  3: V1Cint = RIntegralXX< 6>; break;
            case  4: V1Cint = RIntegralXX< 8>; break;
            case  5: V1Cint = RIntegralXX<10>; break;
            case  6: V1Cint = RIntegralXX<12>; break;
            case  7: V1Cint = RIntegralXX<14>; break;
            case  8: V1Cint = RIntegralXX<16>; break;
            case  9: V1Cint = RIntegralXX<18>; break;
            case 10: V1Cint = RIntegralXX<20>; break;
        }
    }


    //horrible implementation!!!
    int inK2 = 0;
    int nK2  = AB.K2;
    double ps = 0;

    do {
        ++inK2;

        //add next L
        for (int ab1=0; ab1<inK2; ++ab1) {
            for (int ab2=0; ab2<inK2; ++ab2) {

                //the inner square har already been computed; skip it
                if (ab1<inK2-1 && ab2<inK2-1) continue;

                //cout << ab1 << "," << ab2 << " :: ";

                int a1 = AB.nKi[nK2-ab1-1].p;
                int b1 = AB.nKi[nK2-ab1-1].q;
                int a2 = AB.nKi[nK2-ab2-1].p;
                int b2 = AB.nKi[nK2-ab2-1].q;

                //cout << a1 << " " << b1 << " " << a2 << " " << b2 << endl;



                    double ka = A.k[a1];
                    double kb = B.k[b1];
                    double kc = A.k[a2];
                    double kd = B.k[b2];

                    double Kab = ka*kb/(ka+kb);
                    double Kcd = kc*kd/(kc+kd);

                    double kab = ka+kb;
                    double kcd = kc+kd;


                    //add up contributions from all GC functions
                    double c2a = 0;

                    for (int ja=0; ja<A.J; ++ja)
                        c2a += A.w[a1][ja] * A.w[a2][ja];

                    double c2b = 0;

                    for (int jb=0; jb<B.J; ++jb)
                        c2b += B.w[b1][jb] * B.w[b2][jb];


                    double c2 = c2a*c2b * exp(-(Kab+Kcd)*R*R);


                    //compute the trace of (ab|ab)
                    double tr = 0;

                    {
                        // A is in 0, B is in R
                        double P =  R * kb/(ka+kb);
                        double Q =  R * kd/(kc+kd);


                        //expand x^n in terms of hermite functions
                        double HexpAB[Lab+1][Lab+1]; {

                            for (int a=0; a<=Lab; ++a)
                                for (int p=0; p<=Lab; ++p)
                                    HexpAB[a][p] = 0;

                            HexpAB[0][0] = 1;
                            HexpAB[1][1] = -(2*kab);

                            for (int i=2; i<=Lab; ++i) {
                                for (int p=1; p<=i; ++p)
                                    HexpAB[i][p] = -2*kab*HexpAB[i-1][p-1];
                                for (int p=0; p<i-1; ++p)
                                    HexpAB[i][p] += (p+1) * HexpAB[i-1][p+1];
                            }
                        }

                        double HexpCD[Lab+1][Lab+1]; {

                            for (int a=0; a<=Lab; ++a)
                                for (int p=0; p<=Lab; ++p)
                                    HexpCD[a][p] = 0;

                            HexpCD[0][0] = 1;
                            HexpCD[1][1] = -(2*kcd);

                            for (int i=2; i<=Lab; ++i) {
                                for (int p=1; p<=i; ++p)
                                    HexpCD[i][p] = -2*kcd*HexpCD[i-1][p-1];
                                for (int p=0; p<i-1; ++p)
                                    HexpCD[i][p] += (p+1) * HexpCD[i-1][p+1];
                            }
                        }


                        double iHexpAB[Lab+1][Lab+1]; {

                            for (int a=0; a<=Lab; ++a)
                                for (int p=0; p<=Lab; ++p)
                                    iHexpAB[a][p] = 0;

                            for (int a=0; a<=Lab; ++a) {
                                double hpol[Lab+1];
                                for (int p=0; p<=Lab; ++p) hpol[p] = 0; hpol[a] = 1;

                                for (int p=a; p>=0; --p) {
                                    double r = hpol[p]/HexpAB[p][p];

                                    iHexpAB[a][p] = r;

                                    for (int i=0; i<=p; ++i) hpol[i] -= r*HexpAB[p][i];
                                }
                            }
                        }

                        double iHexpCD[Lab+1][Lab+1]; {

                            for (int a=0; a<=Lab; ++a)
                                for (int p=0; p<=Lab; ++p)
                                    iHexpCD[a][p] = 0;

                            for (int a=0; a<=Lab; ++a) {
                                double hpol[Lab+1];
                                for (int p=0; p<=Lab; ++p) hpol[p] = 0; hpol[a] = 1;

                                for (int p=a; p>=0; --p) {
                                    double r = hpol[p]/HexpCD[p][p];

                                    iHexpCD[a][p] = r;

                                    for (int i=0; i<=p; ++i) hpol[i] -= r*HexpCD[p][i];
                                }
                            }
                        }


                        double PexpAB[La+1][Lb+1][Lab+1]; {

                            for (int a=0; a<=La; ++a)
                                for (int b=0; b<=Lb; ++b)
                                    for (int p=0; p<=Lab; ++p)
                                        PexpAB[a][b][p] = 0;

                            PexpAB[0][0][0] = 1;
                            for (int a=1; a<=La; ++a) {
                                for (int p=0; p<a; ++p) PexpAB[a][0][p]    = -P * PexpAB[a-1][0][p];
                                for (int p=0; p<a; ++p) PexpAB[a][0][p+1] +=      PexpAB[a-1][0][p];

                            }

                            for (int b=1; b<=Lb; ++b) {
                                for (int a=0; a<=La; ++a) {
                                    for (int p=0; p<a+b; ++p) PexpAB[a][b][p]    = (R-P) * PexpAB[a][b-1][p];
                                    for (int p=0; p<a+b; ++p) PexpAB[a][b][p+1] +=         PexpAB[a][b-1][p];
                                }
                            }
                        }

                        double PexpCD[La+1][Lb+1][Lab+1]; {

                            for (int a=0; a<=La; ++a)
                                for (int b=0; b<=Lb; ++b)
                                    for (int p=0; p<=Lab; ++p)
                                        PexpCD[a][b][p] = 0;

                            PexpCD[0][0][0] = 1;
                            for (int a=1; a<=La; ++a) {
                                for (int p=0; p<a; ++p) PexpCD[a][0][p]    = -Q * PexpCD[a-1][0][p];
                                for (int p=0; p<a; ++p) PexpCD[a][0][p+1] +=      PexpCD[a-1][0][p];

                            }

                            for (int b=1; b<=Lb; ++b) {
                                for (int a=0; a<=La; ++a) {
                                    for (int p=0; p<a+b; ++p) PexpCD[a][b][p]    = (R-Q) * PexpCD[a][b-1][p];
                                    for (int p=0; p<a+b; ++p) PexpCD[a][b][p+1] +=         PexpCD[a][b-1][p];
                                }
                            }
                        }


                        double YYexpAB[2*La+1][2*Lb+1][Lab+1][Lab+1][Lab+1]; {
                            //expand
                            for (int ma=0; ma<2*La+1; ++ma) {
                                for (int mb=0; mb<2*Lb+1; ++mb) {
                                    for (int xx=0; xx<=Lab; ++xx)
                                    for (int yy=0; yy<=Lab; ++yy)
                                    for (int zz=0; zz<=Lab; ++zz)
                                        YYexpAB[ma][mb][xx][yy][zz] = 0;

                                    for (int i=0; i<SHList[La][ma].nps; ++i) {
                                        for (int j=0; j<SHList[Lb][mb].nps; ++j) {
                                            const SHTerm & sha = SHList[La][ma].T[i];
                                            const SHTerm & shb = SHList[Lb][mb].T[j];
                                            int xx = sha.nx + shb.nx;
                                            int yy = sha.ny + shb.ny;
                                            int zz = sha.nz + shb.nz;

                                            for (int zc=0; zc<=zz; ++zc)
                                                YYexpAB[ma][mb][xx][yy][zc] += sha.cN * shb.cN * PexpAB[sha.nz][shb.nz][zc];
                                        }
                                    }

                                }
                            }
                        }

                        double YYexpCD[2*La+1][2*Lb+1][Lab+1][Lab+1][Lab+1]; {
                            //expand
                            for (int ma=0; ma<2*La+1; ++ma) {
                                for (int mb=0; mb<2*Lb+1; ++mb) {
                                    for (int xx=0; xx<=Lab; ++xx)
                                    for (int yy=0; yy<=Lab; ++yy)
                                    for (int zz=0; zz<=Lab; ++zz)
                                        YYexpCD[ma][mb][xx][yy][zz] = 0;

                                    for (int i=0; i<SHList[La][ma].nps; ++i) {
                                        for (int j=0; j<SHList[Lb][mb].nps; ++j) {
                                            const SHTerm & sha = SHList[La][ma].T[i];
                                            const SHTerm & shb = SHList[Lb][mb].T[j];
                                            int xx = sha.nx + shb.nx;
                                            int yy = sha.ny + shb.ny;
                                            int zz = sha.nz + shb.nz;

                                            for (int zc=0; zc<=zz; ++zc)
                                                YYexpCD[ma][mb][xx][yy][zc] += sha.cN * shb.cN * PexpCD[sha.nz][shb.nz][zc];
                                        }
                                    }

                                }
                            }
                        }


                        double T12r[Lab+1][Lab+1][2*Lab+1]; {

                            for (int a=0; a<=Lab; ++a)
                                for (int b=0; b<=Lab; ++b) {

                                    for (int r=0; r<=2*Lab; ++r)
                                        T12r[a][b][r] = 0;

                                    for (int p=0; p<=Lab; ++p)
                                        for (int q=0; q<=Lab; ++q)
                                            if (q%2==0) T12r[a][b][p+q] += iHexpAB[a][p] * iHexpCD[b][q];
                                            else        T12r[a][b][p+q] -= iHexpAB[a][p] * iHexpCD[b][q];
                                }
                        }



                        double R3[2*Lab+1][2*Lab+1][2*Lab+1]; {

                            for (int x=0; x<=2*Lab; ++x)
                                for (int y=0; y<=2*Lab-x; ++y)
                                    for (int z=0; z<=2*Lab-x-y; ++z)
                                        R3[x][y][z] = 0;

                            /*
                            for (int x=0; x<=2*Lab; ++x)
                                for (int y=0; y<=2*Lab-x; ++y)
                                    for (int z=0; z<=2*Lab-x-y; ++z)
                                        R3[x][y][z] += T12r[a][p][x]*T12r[b][q][y]*T12r[c][r][z] *  YYexpAB[ma][mb][a][b][c] * YYexpCD[ma][mb][p][q][r];
                            */


                            for (int a=0; a<=Lab; ++a) {
                                for (int p=0; p<=Lab; ++p) {

                                    double Ryz[2*Lab+1][2*Lab+1];

                                    for (int y=0; y<=2*Lab; ++y)
                                        for (int z=0; z<=2*Lab-y; ++z)
                                            Ryz[y][z] = 0;


                                    for (int b=0; b<=Lab-a; ++b) {
                                        for (int q=0; q<=Lab-p; ++q) {

                                            double Rzz[2*Lab+1];
                                            for (int z=0; z<=2*Lab; ++z) Rzz[z] = 0;

                                            for (int c=0; c<=Lab-a-b; ++c) {
                                                for (int r=0; r<=Lab-p-q; ++r) {
                                                    double sum_mn = 0;

                                                    for (int ma=0; ma<2*La+1; ++ma)
                                                        for (int mb=0; mb<2*Lb+1; ++mb)
                                                            sum_mn += YYexpAB[ma][mb][a][b][c] * YYexpCD[ma][mb][p][q][r];

                                                    for (int z=0; z<=2*Lab; ++z)
                                                        Rzz[z] += T12r[c][r][z] * sum_mn;
                                                }
                                            }

                                            for (int y=0; y<=2*Lab; ++y)
                                                for (int z=0; z<=2*Lab-y; ++z)
                                                    Ryz[y][z] += T12r[b][q][y]*Rzz[z];
                                        }
                                    }

                                    for (int x=0; x<=2*Lab; ++x)
                                        for (int y=0; y<=2*Lab-x; ++y)
                                            for (int z=0; z<=2*Lab-x-y; ++z)
                                                R3[x][y][z] += T12r[a][p][x]*Ryz[y][z];

                                }
                            }
                        }


                        //compute the 1C integral
                        TripleXX TXX;
                        V1Cint(TXX, kab, kcd, Q-P);

                        for (int x=0; x<=2*Lab; ++x)
                            for (int y=0; y<=2*Lab-x; ++y)
                                for (int z=0; z<=2*Lab-x-y; ++z)
                                    tr += R3[x][y][z] * TXX[x][y][z];

                    }

                    ps += c2 * tr;

            }
        }

        //cout << endl;
        //cout << "         " << ps << " " << thresh << endl;

        if (ps>thresh) return nK2-inK2; //residual error above threshold; include density
        if (inK2==nK2) return -1;       //all densities below theshold;
    } while (1);

}


void ABABinit() {

    ABABps[0][0] = CalcK2ps<0,0>;
    ABABps[0][1] = CalcK2ps<0,1>;
    ABABps[0][2] = CalcK2ps<0,2>;
    ABABps[0][3] = CalcK2ps<0,3>;
    ABABps[0][4] = CalcK2ps<0,4>;
    ABABps[0][5] = CalcK2ps<0,5>;

    ABABps[1][0] = CalcK2ps<1,0>;
    ABABps[1][1] = CalcK2ps<1,1>;
    ABABps[1][2] = CalcK2ps<1,2>;
    ABABps[1][3] = CalcK2ps<1,3>;
    ABABps[1][4] = CalcK2ps<1,4>;
    ABABps[1][5] = CalcK2ps<1,5>;

    ABABps[2][0] = CalcK2ps<2,0>;
    ABABps[2][1] = CalcK2ps<2,1>;
    ABABps[2][2] = CalcK2ps<2,2>;
    ABABps[2][3] = CalcK2ps<2,3>;
    ABABps[2][4] = CalcK2ps<2,4>;
    ABABps[2][5] = CalcK2ps<2,5>;

    ABABps[3][0] = CalcK2ps<3,0>;
    ABABps[3][1] = CalcK2ps<3,1>;
    ABABps[3][2] = CalcK2ps<3,2>;
    ABABps[3][3] = CalcK2ps<3,3>;
    ABABps[3][4] = CalcK2ps<3,4>;
    ABABps[3][5] = CalcK2ps<3,5>;

    ABABps[4][0] = CalcK2ps<4,0>;
    ABABps[4][1] = CalcK2ps<4,1>;
    ABABps[4][2] = CalcK2ps<4,2>;
    ABABps[4][3] = CalcK2ps<4,3>;
    ABABps[4][4] = CalcK2ps<4,4>;
    ABABps[4][5] = CalcK2ps<4,5>;

    ABABps[5][0] = CalcK2ps<5,0>;
    ABABps[5][1] = CalcK2ps<5,1>;
    ABABps[5][2] = CalcK2ps<5,2>;
    ABABps[5][3] = CalcK2ps<5,3>;
    ABABps[5][4] = CalcK2ps<5,4>;
    ABABps[5][5] = CalcK2ps<5,5>;
}
