#ifndef __ABAB_PS__
#define __ABAB_PS__

#include "consts.hpp"

namespace LibAngularVPD {
    //término en una expansion en cartesianas de un armónico esférico
    struct SHTerm {
        double cN; //peso del armonico normalzado
        char   nc; //numero del termino correspondiente en cartesianas
        char   nx; //potencia en x
        char   ny; //potencia en y
        char   nz; //potencia en z
    };

    //definicion de un armónico esférico (real) en coordenadas cartesianas
    struct SphHar {
        int l; //momento angular
        int m; //m

        int nps;     //numero de primitivas cartesianas para construir el harmonico esferico
        SHTerm T[6]; //terminos individuales
                     //allow only up to 6 at the moment, which is enough for up to h-functions;
                     //this avoids unnecesary memory fragmentation

        //inicia los arrays de terminos
        void SetNPS(int n);
        void SetT(int t, double c, int x, int y, int z);
    };


    extern SHTerm CartList[LMAX+2][MMAX];
    extern SphHar SHList  [LMAX+2][2*LMAX+1];

    // inicia las tablas de esfericos harmonicos y las seminormaliza
    void InitSHList();
    void InitCartList();
}

namespace LibIGammaVPD {

    const double PI3  = 31.0062766802998202;   // = PI^3
    const double PI54 = 34.986836655249725693; // = 2*PI^(5/2)

    const float fPI3  = 31.0062766802998202;   // = PI^3
    const float fPI54 = 34.986836655249725693; // = 2*PI^(5/2)

    // everything is made public since
    // the inlined and specialized algorithms perform better

    // since gamma is evaluated by downward recursion, it is in principle better to have
    // one GammaFunction for each Fm[z] and adjust the number of points in the grid, grid step
    // and critical gamma acordingly; this will eventually be implemented

    const int NEXP = 4;
    const int NVIG = 1024;
    const int NGD = 32; // >= 4*LMAX + 1 + NEXP + 2;
    const int NGD2 = 8;

    const double vg_min  = 0;
    const double vg_max  = 32;
    const double vg_step  = 32./1024.;  //(vg_max-vg_min)/double(NVIG);
    const double ivg_step = 1024./32.;    //double(NVIG)/(vg_max-vg_min);

    const float fvg_min  = 0;
    const float fvg_max  = 32;
    const float fvg_step  = 32./1024.;  //(vg_max-vg_min)/double(NVIG);
    const float fivg_step = 1024./32.;    //double(NVIG)/(vg_max-vg_min);



    class GammaFunction {
      public:

      public:
        double ** gamma_table;
        double *  gamma_table_vals;
      public:
        GammaFunction();
        ~GammaFunction();
        void InitIncompleteGammaTable();
        void InitIncompleteGammaTable(int L);
        void CalcGammas(double * F, int n, double x) const;
    };

    extern GammaFunction IncompleteGamma;
    extern GammaFunction IncompleteGammas[4*LMAX+3];

    void InitIncompleteGammas();
}



struct eBasis {
    //double w[maxJ][maxK];
    double w[maxK][maxJ];
    double k[maxK];
    int    L;
    int    J;
    int    K;


    eBasis(); //nullify everything
    void sort(); //sort the primtives in ascending order
    bool operator<(const eBasis & rhs) const; //for sorting
};

struct nKinteracting {
    unsigned char nKa[maxK];
    double ikappa;
    int p;
    int q;
    unsigned char nKb;

    nKinteracting() {
        for (int k=0; k<maxK; ++k) nKa[k] = 0;
        nKb = 0;
    }
};

struct eBasis2 {
    int La;
    int Ja;
    int Ka;
    int Lb;
    int Jb;
    int Kb;

    int K2;

    nKinteracting * nKi;

    void set(const eBasis & Abasis, const eBasis & Bbasis);
};


typedef double (*GTO2PS) (const eBasis &, const eBasis &, double);
typedef int    (*K2PS)  (const eBasis &, const eBasis &, const eBasis2 & AB, double, double);

extern K2PS ABABps[6][6];

void ABABinit();

#endif
