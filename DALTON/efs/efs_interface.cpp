/*
    Wrapping of Echidna Fock Solver interface to be able to use the library from DALTON

    Jaime Axel Rosal Sandberg, somewhere mid-2013
*/

#if defined(VAR_MPI)
  #include "mpi.h"
#endif

#include <iostream>
#include <sstream>
#include <set>
#include <queue>
#include <map>

#include "libquimera/libquimera.hpp"
#include "libechidna/libechidna.hpp"
#include "math/tensors.hpp"
#include "low/MPIwrap.hpp"

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
    #define omp_get_num_threads() 1
#endif

using namespace std;


class MessagePlotter;
class AtomDefList;

//datatypes and objects necessary to communicate with dalton
namespace EchidnaDaltonInterface {

    const int maxJK = 35; // taken from DALTON
    const int maxJK2 = maxJK*maxJK;

    struct DIgto {
       int angular_momentum;
       int nr_primitives;
       int nr_contractions;

       double primitive_exp    [maxJK];
       double contraction_coef [maxJK][maxJK];
    };

    struct DIpoint {
        double x;
        double y;
        double z;
    };

    struct DIatomtype {
        int nc;
        int ng;
        DIpoint * center;
        DIgto   * gtos;

        DIatomtype() {
            nc = ng = 0;
            center = NULL;
            gtos   = NULL;
        }
    };


    deque<DIatomtype> q_atypes;  // temporal list of objects

    tensor2 tF, tD;               // Fock and Density matrix wrappers for Fortran 2D arrays

    rAtomPrototypes InternalDefs; // name, charge, basis set, etc.
    rAtom * Atoms;                // atom list
    int nAtoms;
    int nbasis;

    int nCycle;

    Echidna FockSolver;           // Fock Solver
};

using namespace EchidnaDaltonInterface;


inline double RadialIntegral(double k, int L) {

    const double hhGamma[] =
        {1.77245385090551603, 0.886226925452758014, 1.32934038817913702, 3.32335097044784255, 11.6317283965674489, 52.3427777845535202,
        287.885277815044361, 1871.25430579778835, 14034.4072934834126, 119292.461994609007, 1.13327838894878557e6,
         1.18994230839622485e7, 1.36843365465565857e8, 1.71054206831957322e9, 2.30923179223142384e10,
         3.34838609873556457e11, 5.18999845304012508e12, 8.56349744751620639e13, 1.49861205331533612e15,
         2.77243229863337182e16, 5.40624298233507504e17, 1.10827981137869038e19, 2.38280159446418433e20,
         5.36130358754441473e21, 1.25990634307293746e23, 3.08677054052869678e24, 7.87126487834817680e25,
         2.08588519276226685e27, 5.73618428009623384e28, 1.63481251982742664e30, 4.82269693349090860e31,
         1.47092256471472712e33};

    double kn = 1;
    for (int i=0; i<2*L+3; ++i)
        kn *= k;
    return hhGamma[L+1] / (2*sqrt(kn));
}

GTO::GTO () {
    SP = false;
}

GTO::~GTO () {
}

void GTO::Normalize() {

    //sort primitives by ascending gaussian exponent
    for (int kk=0; kk<K; ++kk) {
        double minW = k[kk];
        int kmin = kk;

        for (int jj=kk; jj<K; ++jj) if (k[jj]<minW) {minW = k[jj]; kmin = jj;}

        //swap everything
        swap(k[kk], k[kmin]);

        for (int jj=0; jj<J; ++jj)
            swap(N[jj][kk], N[jj][kmin]);

        if (SP)
        for (int jj=0; jj<J; ++jj)
            swap(Np[jj][kk], Np[jj][kmin]);

    }


    // loop over all contracted functions
    // and normalize the radial part
    for (int jj=0; jj<J; ++jj) {

        double cp=0.;
        for (int i=0;i<K;++i) {
            cp += N[jj][i]*N[jj][i]*RadialIntegral(2*k[i],l);

            for (int j=0;j<i;++j) {
                cp += 2*N[jj][i]*N[jj][j]*RadialIntegral(k[i]+k[j],l);
            }
        }

        cp = 1./sqrt(cp);

        for (int i=0; i<K; ++i)
            N[jj][i] *= cp;
    }

    return;
}

void PrintGTO(const GTO & g) {
    cout << "    angular_momentum " << int(g.l) << endl;
    cout << "    nr_primitives    " << int(g.K) << endl;
    cout << "    nr_contractions  " << int(g.J) << endl;
    cout << "    exponents :: contraction coefs  " << endl;

    for(int k=0; k<g.K; ++k) {
        cout << "    " << g.k[k] << "  ::  ";
        for(int j=0; j<g.J; ++j) {
            cout << g.N[j][k] << " ";
        }
        cout << endl;
    }

}


int popcount(int n) {
    int m = 0;

    while (n>0) {
        m += n%2;
        n /=   2;
    }

    return m;
}

void PrintGTO(const DIgto & g) {
    cout << "    angular_momentum " << g.angular_momentum << endl;
    cout << "    nr_primitives    " << g.nr_primitives << endl;
    cout << "    nr_contractions  " << g.nr_contractions << endl;
    cout << "    exponents :: contraction coefs  " << endl;

    for(int i=0; i<g.nr_primitives; ++i) {
        cout << "    " << g.primitive_exp[i] << "  ::  ";
        for(int j=0; j<g.nr_contractions; ++j) {
            cout << g.contraction_coef[i][j] << " ";
        }
        cout << endl;
    }

}


static const long long int CFLAGS[] =
{      1,        2,        4,         8,        16,        32,         64,        128,
     256,      512,     1024,      2048,      4096,      8192,      16384,      32768,
   65536,   131072,   262144,    524288,   1048576,   2097152,    4194304,    8388608,
16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824, 2147483648 };

int CountGTOs (const DIgto & g) {

    const int K = g.nr_primitives;
    const int J = g.nr_contractions;

    int ngtos = 0;
    long long int ocode = 0;

    for (int j=0; j<J; ++j) {

        //generate a binary code for each function depending on the primitives used
        long long int code=0;

        for (int k=0; k<K; ++k)
            if (g.contraction_coef[k][j] != 0.)
                code += CFLAGS[k]; //(1<<k);

        if (code!=ocode) ++ngtos;
        ocode=code;
    }

    return ngtos;
}

void SplitGTOs (const DIgto & g, DIgto * gtos) {

    const int L = g.angular_momentum;
    const int K = g.nr_primitives;
    const int J = g.nr_contractions;

    //generate a binary code for each function depending on the primitives used
    int n = 0;
    int jj = 0;
    long long int ocode = 0;

    for (int j=0; j<J; ++j) {

        long long int code; {
            code=0;

            for (int k=0; k<K; ++k)
                if (g.contraction_coef[k][j] != 0.)
                    code += CFLAGS[k]; //(1<<k);
        }

        //new function
        if (code!=ocode) {
            DIgto & ngto = gtos[n];
            ++n;

            ngto.angular_momentum = L;
            ngto.nr_primitives    = popcount(code);
            ngto.nr_contractions  = 0; //will increase for every function counted

            //copy the exponents
            int kk=0;
            for (int k=0; k<K; ++k) {
                if (g.contraction_coef[k][j] != 0) {
                    ngto.primitive_exp[kk] = g.primitive_exp[k];
                    ++kk;
                }
            }

            jj=0;
        }

        //copy coefficients for the given contraction
        {
            DIgto & ngto = gtos[n-1];

            int kk=0;
            for (int k=0; k<K; ++k) {
                if (g.contraction_coef[k][j] != 0) {
                    ngto.contraction_coef[kk][jj] = g.contraction_coef[k][j];
                    ++kk;
                }
            }

            ++ngto.nr_contractions;
            ++jj;
        }

        ocode=code;
    }

}


void TriSwap(tensor2 & T, int pos) {

    //swap in first dimension
    for (int i=0; i<T.n; ++i) {
        double tmp = T(pos,i);
        T(pos  ,i) = T(pos+1,i);
        T(pos+1,i) = T(pos+2,i);
        T(pos+2,i) = tmp;
    }

    //swap in second dimension
    for (int i=0; i<T.n; ++i) {
        double tmp = T(i,pos);
        T(i,pos  ) = T(i,pos+1);
        T(i,pos+1) = T(i,pos+2);
        T(i,pos+2) = tmp;
    }
}

void iTriSwap(tensor2 & T, int pos) {

    //swap in first dimension
    for (int i=0; i<T.n; ++i) {
        double tmp = T(pos,i);
        T(pos  ,i) = T(pos+2,i);
        T(pos+2,i) = T(pos+1,i);
        T(pos+1,i) = tmp;
    }

    //swap in second dimension
    for (int i=0; i<T.n; ++i) {
        double tmp = T(i,pos);
        T(i,pos  ) = T(i,pos+2);
        T(i,pos+2) = T(i,pos+1);
        T(i,pos+1) = tmp;
    }
}

void AdaptTensor(tensor2 & T) {

    int nb = 0;

    for (int nat=0; nat<nAtoms; ++nat) {
        const rAtomPrototype & rAP = *Atoms[nat].rAP;

        //loop over all basis of given atom
        for (int g=0; g<rAP.basis.N(); ++g) {

            int j = rAP.basis[g].J;
            int l = rAP.basis[g].l;

            // rotate p functions
            if (l==1) {
                for (int jj=0; jj<j; ++jj) {
                    TriSwap(T, nb);
                    nb += 3;
                }
            }
            else
                nb += j * (2*l+1);
        }
    }

}

void iAdaptTensor(tensor2 & T) {

    int nb = 0;

    for (int nat=0; nat<nAtoms; ++nat) {
        const rAtomPrototype & rAP = *Atoms[nat].rAP;

        //loop over all basis of given atom
        for (int g=0; g<rAP.basis.N(); ++g) {

            int j = rAP.basis[g].J;
            int l = rAP.basis[g].l;

            // rotate p functions
            if (l==1) {
                for (int jj=0; jj<j; ++jj) {
                    iTriSwap(T, nb);
                    nb += 3;
                }
            }
            else
                nb += j * (2*l+1);
        }
    }

}

string GetEnv2(const string & var) {
    const char * val = ::getenv(var.c_str());

    if (val== 0) return "";
    else         return val;
}



//functions used by Dalton
extern "C" {

    //functions called by the master MPI node
    void efs_add_atomtype_(DIpoint * center, DIgto * gtos, const int * nc, const int * ng);
    void efs_print_();

    void efs_init_();
    void efs_generate_basis_product_(double * logGDO); // log of the Gaussian Product threshold
    void efs_init_2efock_(double * Xscaling);
    void efs_fock_update_(double * F, double * D, int * ifthresh, int * inmat, int * ityp);
    void efs_add_2efock_(double * F, double * D0, int * ifthresh);

    void efs_print_matrices_(double * F, double * D, int * nbs, int * inmat);


    //functions called by slave threads in eri2par
    void efs_add_atomtype_slave_();
    void efs_init_slave_();
    void efs_generate_basis_product_slave_();
    void efs_init_2efock_slave_();
    void efs_fock_update_slave_();


    //functions exported by the F90 side of the interface
    void efsprint_ (int * length, const char * str);
    void efserror_ (int * length, const char * str);
}

void Print2Dalton(const string & str) {

    int len = str.length();
    efsprint_ (&len, str.c_str());
}

void Error2Dalton(const string & str) {

    int len = str.length();
    efserror_ (&len, str.c_str());
}

void Error2Dalton(int ecode) {
    string str =  " Echidna Fock Solver unexpectedly exited with error code ";

    ostringstream ss;
    ss << ecode;
    str += ss.str();

    Print2Dalton(str);
    Error2Dalton(str);
}


/*
    Overloaded MPI short-hand
*/

inline void mpi_bcast(int & val) {
  #if defined(VAR_MPI)
    //const int one = 1;
    //int ierr;
    //mpi_bcast_ (&nc, &one, &mpi_integer, &mpi_id_master, &mpi_comm, &ierr);

     MPI_Bcast (&val, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
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


/*
    Interface functions for DALTON
*/

void efs_add_atomtype_slave_() {

    //received data
    int     nc;
    int     ng;
    DIgto * gtos;

    //get nc and ng
    mpi_bcast(nc);
    mpi_bcast(ng);


    //other data
    DIatomtype at;

    //copy center information
    at.nc     = nc;
    at.center = new DIpoint[nc];


    for (int i=0; i<nc; ++i) {
        mpi_bcast(at.center[i].x);
        mpi_bcast(at.center[i].y);
        mpi_bcast(at.center[i].z);
    }


    gtos = new DIgto[ng];

    for (int i=0; i<ng; ++i) {
        mpi_bcast ( gtos[i].angular_momentum );
        mpi_bcast ( gtos[i].nr_primitives    );
        mpi_bcast ( gtos[i].nr_contractions  );

        mpi_bcast ( (double*) gtos[i].primitive_exp,      maxJK  );
        mpi_bcast ( (double*) gtos[i].contraction_coef,   maxJK2 );
    }


    //copy GTOs information; split GTO blocks and make a unique GTO array
    at.ng = 0;
    for (int i=0; i<ng; ++i) at.ng += CountGTOs(gtos[i]);

    at.gtos   = new DIgto[at.ng];

    for (int i=0, nng=0; i<ng; ++i) {
        SplitGTOs(gtos[i], at.gtos + nng);
        nng += CountGTOs(gtos[i]);
    }

    q_atypes.push_back(at);

}


void efs_add_atomtype_(DIpoint * center, DIgto * gtos, const int * pnc, const int * png) {

    int nc = *pnc;
    int ng = *png;

  //broadcast everything
  {
    mpi_bcast (nc); //get nc
    mpi_bcast (ng); //get ng

    for (int i=0; i<nc; ++i) {
        mpi_bcast(center[i].x);
        mpi_bcast(center[i].y);
        mpi_bcast(center[i].z);
    }

    for (int i=0; i<ng; ++i) {
        mpi_bcast ( gtos[i].angular_momentum );
        mpi_bcast ( gtos[i].nr_primitives    );
        mpi_bcast ( gtos[i].nr_contractions  );

        mpi_bcast ( (double*) gtos[i].primitive_exp,      maxJK  );
        mpi_bcast ( (double*) gtos[i].contraction_coef,   maxJK2 );
    }
  }

    //process information
    //+++++++++++++++++++

    DIatomtype at;

    //copy center information
    at.nc     = nc;
    at.center = new DIpoint[at.nc];
    for (int i=0; i<nc; ++i) at.center[i] = center[i];

    //copy GTOs information; split GTO blocks and make a unique GTO array
    at.ng = 0;
    for (int i=0; i<ng; ++i) at.ng += CountGTOs(gtos[i]);

    at.gtos   = new DIgto[at.ng];

    for (int i=0, nng=0; i<ng; ++i) {
        SplitGTOs(gtos[i], at.gtos + nng);
        nng += CountGTOs(gtos[i]);
    }

    q_atypes.push_back(at);
}

void efs_print_() {
    deque<DIatomtype>::iterator it;

    int nat = 1;

    for (it=q_atypes.begin(); it!=q_atypes.end(); ++it) {
        const DIatomtype & AT = *it;

        cout << "ATOM TYPE no " << nat << endl;

        cout << "  Centers" << endl;
        for (int i=0; i<AT.nc; ++i)
            cout << "    " << i << "    " << AT.center[i].x << " " << AT.center[i].y << " " << AT.center[i].z << endl;
        cout << endl;

        cout << "  Basis" << endl;
        for (int i=0; i<AT.ng; ++i) {
            cout << "   Shell " << (i+1) << endl;
            PrintGTO(AT.gtos[i]);
            cout << endl;
        }
        cout << endl;
        cout << endl;

        ++nat;
    }

}


void efs_init_slave_() {
    efs_init_();
}


void efs_init_() {

    //init MPI information
    ThisNode.Init();

    //check MPI processes and OpenMP threads work
    /*
    #pragma omp parallel
    {
        cout << "Node " << ThisNode.rank << " init; Thread " << omp_get_thread_num() << endl;
        cout.flush();
    }
    */


    deque<DIatomtype>::iterator it;
    int nat;

    //merge all atoms into a single list
    //==================================
    nAtoms = 0;

    for (it=q_atypes.begin(); it!=q_atypes.end(); ++it)
        nAtoms += it->nc;

    Atoms = new rAtom[nAtoms];

    for (nat=0, it=q_atypes.begin(); it!=q_atypes.end(); ++it) {
        for (int i=0;i<it->nc;++i,++nat) {
            Atoms[nat].c.x = it->center[i].x;
            Atoms[nat].c.y = it->center[i].y;
            Atoms[nat].c.z = it->center[i].z;
        }
    }



    //make the prototypes && link them to each atom in the list
    //=========================================================
    int Natypes = q_atypes.size();

    InternalDefs.prototypes.set(Natypes);

    int id = 0;

    for (nat=0,it=q_atypes.begin(); it!=q_atypes.end(); ++it) {

        //get this information from the DALTON input, if possible
        InternalDefs.prototypes[id].name  = "";
        InternalDefs.prototypes[id].mass  =  0 ;
        InternalDefs.prototypes[id].Z     =  0 ;

        //in DALTON both values are the same
        InternalDefs.prototypes[id].DId   = id;           // Id in the original definition list
        InternalDefs.prototypes[id].RId   = id;           // Id in the reduced element list

        //copy the basis set
        InternalDefs.prototypes[id].basis.set(it->ng);

        for (int i=0; i<it->ng; ++i) {
            GTO & g = InternalDefs.prototypes[id].basis[i];

            g.l = it->gtos[i].angular_momentum;
            g.K = it->gtos[i].nr_primitives;
            g.J = it->gtos[i].nr_contractions;

            for (int k=0; k<g.K; ++k)
                g.k[k] = it->gtos[i].primitive_exp[k];

            for (int j=0; j<g.J; ++j)
                for (int k=0; k<g.K; ++k)
                    g.N[j][k] = it->gtos[i].contraction_coef[k][j];

            g.SP = false;

            //radial normalization && SORT EXPONENTS!!!!
            g.Normalize(); // hope this works
        }


        InternalDefs.prototypes[id].nAtoms = it->nc;

        for (int i=0;i<it->nc;++i) {
            Atoms[nat+i].rAP = &(InternalDefs.prototypes[id]);
        }

        InternalDefs.prototypes[id].Atoms = Atoms + nat;

        nat += it->nc;
        ++id;
    }


    //clean the mess
    for (nat=0,it=q_atypes.begin(); it!=q_atypes.end(); ++it) {
        delete[] it->center;
        delete[] it->gtos;
    }
    q_atypes.clear();


// small check everything is passed OK
/*
    for (int i=0; i<InternalDefs.prototypes.N(); ++i) {
        cout << "Atom " << i << endl;
        const rAtomPrototype & rAP = InternalDefs.prototypes[i];
	for (int b=0; b<rAP.basis.N(); ++b) {
            cout << "   " << b << endl;
            PrintGTO(rAP.basis[b]);
            cout << endl;
        }
    }
    cout << endl;

    for (int i=0; i<nAtoms; ++i) {
        cout << i << " " << Atoms[i].c.x << " " << Atoms[i].c.y << " " << Atoms[i].c.z  << endl;
    }

    char aa;
    cin >> aa;
*/


    try {
        FockSolver.Init (Atoms, nAtoms);
    }
    catch (int ecode) {
        Error2Dalton(ecode);
    }

}



void efs_generate_basis_product_slave_() {
    double nada;
    efs_generate_basis_product_(&nada);
}

void efs_generate_basis_product_(double * GDO) {

    double logGDO;

    if (ThisNode.IsMaster()) {
        logGDO = -log(*GDO);
    }

    mpi_bcast (logGDO);

    //compute atom pairs for prescreening

    //find smallest exponent
    double * iminK = new double[InternalDefs.prototypes.N()];

    for (int i=0; i<InternalDefs.prototypes.N(); ++i) {
        double kmin = InternalDefs.prototypes[i].basis[0].k[0];

        for (int g=0; g<InternalDefs.prototypes[i].basis.N(); ++g) {
            for (int k=0; k<InternalDefs.prototypes[i].basis[g].K; ++k) {
                kmin = min(kmin, InternalDefs.prototypes[i].basis[g].k[k]);
            }
        }

        iminK[i] = 1./kmin;
    }


    //count interactions for each atom
    unsigned short int * Ninter = new unsigned short int[nAtoms];

    for (int at1=0; at1<nAtoms; ++at1) {
        int nint = 0;

        for (int at2=0; at2<nAtoms; ++at2) {
            vector3 r = Atoms[at2].c - Atoms[at1].c;
            double R2 = r*r;
            double iKpq = iminK[Atoms[at1].rAP->DId] + iminK[Atoms[at2].rAP->DId]; // 1/Kpq = 1/Kp + 1/Kq

            if (R2<=(logGDO)*iKpq) ++nint;
        }

        Ninter[at1] = nint;
    }


    //generate a list8.08187 of all atom pairs which have at least pair of GTOs
    r2tensor<int> AtomPairs;

    AtomPairs.set(Ninter, nAtoms);

    for (int at1=0; at1<nAtoms; ++at1) {
        int nint = 0;

        for (int at2=0; at2<nAtoms; ++at2) {
            vector3 r = Atoms[at2].c - Atoms[at1].c;
            double R2 = r*r;
            double iKpq = iminK[Atoms[at1].rAP->DId] + iminK[Atoms[at2].rAP->DId]; // 1/Kpq = 1/Kp + 1/Kq

            if (R2<=(logGDO)*iKpq) {
                AtomPairs[at1][nint] = at2;
                ++nint;
            }
        }
    }

    delete[] Ninter;
    delete[] iminK;

    try {
        FockSolver.ComputeAtomPairs (AtomPairs, logGDO);
    }
    catch (int ecode) {
        Error2Dalton(ecode);
    }

    //cout << "basis product generated! " << ThisNode.rank << endl;
}



void efs_init_2efock_slave_ () {
    double xs;
    efs_init_2efock_(&xs);
}


#define _GNU_SOURCE
#include <fenv.h>

void efs_init_2efock_(double * Xscaling) {

#ifndef SYS_DARWIN
    // radovan: not available on OS X
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

    if (ThisNode.IsMaster()) {

        //                                  print a splash screen
        //                                  *********************
        Print2Dalton("");
        Print2Dalton("");
        Print2Dalton(" ****************************************************************************");
        Print2Dalton(" *                        Echidna Fock Solver                               *");
        Print2Dalton(" *                Author: Jaime Axel Rosal Sandberg                         *");
        Print2Dalton(" *                                                                          *");
        Print2Dalton(" *   The Echidna Fock Solver is experimental code; use at your own risk!    *");
        Print2Dalton(" *         If you use it in your research, please cite the paper            *");
        Print2Dalton(" *             Read the README file for further information                 *");
        Print2Dalton(" ****************************************************************************");
    }


    int one = 1;
    int ierr;

    mpi_bcast (*Xscaling);

        // set memory limits
        string efs_mem    = GetEnv2("EFS_MEM_MB");
        // export EFS_MEM_MB    = 1024
        if (efs_mem    != "") {
            int mem;
            istringstream iss(efs_mem);
            iss     >> mem;
            Echidna::MAXMEM = mem * Mword;
        }

        string efs_buffer = GetEnv2("EFS_BUFFER_MB");
        // export EFS_BUFFER_MB = 1024
        if (efs_buffer != "") {
            int mem;
            istringstream iss(efs_buffer);
            iss  >> mem;
            Echidna::MEMBUFFER = mem * Mword;
        }

        LibQuimera::IC_DIR = GetEnv2("QIC_DIR");


    nCycle = 0;

    //count the total amount of basis functions
    nbasis = 0;

    for (int i=0; i<InternalDefs.prototypes.N(); ++i) {

        //compute the length of the basis for each element
        int nb = 0;

        for (int g=0; g<InternalDefs.prototypes[i].basis.N(); ++g) {
            int j = InternalDefs.prototypes[i].basis[g].J;
            int l = InternalDefs.prototypes[i].basis[g].l;
            nb += j * (2*l+1);
        }

        //add it times the number of atoms for that element
        nbasis += nb * InternalDefs.prototypes[i].nAtoms;
    }

    //set IO tensor wrappers
    tF.n = tF.m = nbasis;
    tF.c = new double*[nbasis];

    tD.n = tD.m = nbasis;
    tD.c = new double*[nbasis];


    try {
        FockSolver.InitFock2e(*Xscaling);
    }
    catch (int ecode) {
        Error2Dalton(ecode);
    }

    //cout << "2efock initialized!" << endl;
}



// All MPI processes enter Echidna through here
// ********************************************

void fock_update(double * F, double * D, double thresh, int nmat) {

    //get the density matrix(ces) information
    mpi_bcast (D , nmat*nbasis*nbasis);


    //Fock and Density matrix wrappers for Fortran 2D arrays
    tensor2 * tF = new tensor2[nmat];
    tensor2 * tD = new tensor2[nmat];

    for (int n=0; n<nmat; ++n) {
        tensor2 & tD2 = tD[n];
        double * D2 = D + n*nbasis*nbasis;

        tD2.n = tD2.m = nbasis;
        tD2.c = new double*[nbasis];

        tD2.c2 = D2;
        for (int i=0; i<nbasis; ++i) tD2.c[i] = D2 + i*nbasis;
    }

    for (int n=0; n<nmat; ++n) {
        tensor2 & tF2 = tF[n];
        double * F2 = F + n*nbasis*nbasis;

        tF2.n = tF2.m = nbasis;
        tF2.c = new double*[nbasis];

        tF2.c2 = F2;
        for (int i=0; i<nbasis; ++i) tF2.c[i] = F2 + i*nbasis;
    }


    // swap px,py,pz DALTON order in the density matrix to p-1, p0, p+1 order used in the library
    for (int n=0; n<nmat; ++n) AdaptTensor(tD[n]);
    for (int n=0; n<nmat; ++n) tF[n] = tD[n];


    //perform the computation
    try {
        FockSolver.FockUpdate(tF, tD, nmat, thresh, thresh, false);
    }
    catch (int ecode) {
        Error2Dalton(ecode);
    }


    // adapt the order to DALTON format and rescale
    if (ThisNode.IsMaster()) {
        for (int n=0; n<nmat; ++n) iAdaptTensor(tD[n]);
        for (int n=0; n<nmat; ++n) iAdaptTensor(tF[n]);
    }


    //try if swapping the results gives correct answers in DALTON
    for (int n=0; n<nmat; ++n) {
        for (int i=0; i<nbasis; ++i)
            for (int j=0; j<i; ++j)
                swap(tF[n](i,j), tF[n](j,i));
    }


    // dispose of the supporting structures properly
    for (int n=0; n<nmat; ++n) tF[n].c2 = NULL;
    for (int n=0; n<nmat; ++n) tD[n].c2 = NULL;

    delete[] tF;
    delete[] tD;

    ++nCycle;
}



void efs_fock_update_slave_() {

    //get information
    int nmat;
    double thresh;

    mpi_bcast(nmat);
    mpi_bcast(thresh);

    //allocate memory
    double * F = new double[nmat*nbasis*nbasis];
    double * D = new double[nmat*nbasis*nbasis];

    for (int i=0; i<nmat*nbasis*nbasis; ++i) F[i] = 0;

    fock_update(F, D, thresh, nmat);

    delete[] F;
    delete[] D;
}

void efs_fock_update_(double * F, double * D0, int * ifthresh, int * inmat, int * ityp) {

    int nmat;
    int ift = *ifthresh;
    double thresh;

    //DALTON uses ifthresh=20 to indicate no prescreening
    if (ift<=0 || ift>=20) thresh = 1e-50;
    else                   thresh = pow(10., -ift); //the screening method can be too harshon exchange integrals with DALTON's settings

    nmat = *inmat;

    /*
    cout << "threshold: " << thresh << endl;
    cout << "number of matrices: " << nmat << endl;
    cout << "nbasis: " << nbasis << endl;
    cout << "total elements: " << nmat*nbasis*nbasis << endl;
    cout << "pointer to F: " << (void*)F << endl;
    cout << "pointer to D: " << (void*)D0 << endl;
    cout << "matrix types: "; for (int n=0; n<nmat; ++n) cout << ityp[n] << " "; cout << endl;
    */

    mpi_bcast (nmat);
    mpi_bcast (thresh);

    double * D = new double[nmat*nbasis*nbasis];
    for (int i=0; i<nmat*nbasis*nbasis; ++i) D[i] = D0[i];

    for (int i=0; i<nmat*nbasis*nbasis; ++i) F[i] = 0;

    //cout << "in " << endl;
    fock_update(F, D, thresh, nmat);
    //cout << "out" << endl << endl;

    delete[] D;
}

//for simple fock contributions
void efs_add_2efock_(double * F, double * D0, int * ifthresh) {
    int nmat = 1;
    int ityp = 13; // sirfck code for symmetric density; request both coulomb and exchange

    efs_fock_update_(F, D0, ifthresh, &nmat, &ityp);
}



void efs_print_matrices_(double * F, double * D,  int * nbs, int * inmat) {

    nbasis = *nbs;

    int nmat;
    nmat = *inmat;


    cout << "number of matrices: " << nmat << endl;
    cout << "nbasis: " << nbasis << endl;
    cout << "total elements: " << nmat*nbasis*nbasis << endl;
    cout << "pointer to F: " << (void*)F << endl;
    cout << "pointer to D: " << (void*)D << endl;
    cout << endl;


    tensor2 * tF = new tensor2[nmat];
    tensor2 * tD = new tensor2[nmat];

    for (int n=0; n<nmat; ++n) {
        tensor2 & tD2 = tD[n];
        double * D2 = D + n*nbasis*nbasis;

        tD2.n = tD2.m = nbasis;
        tD2.c = new double*[nbasis];

        tD2.c2 = D2;
        for (int i=0; i<nbasis; ++i) tD2.c[i] = D2 + i*nbasis;
    }

    for (int n=0; n<nmat; ++n) {
        tensor2 & tF2 = tF[n];
        double * F2 = F + n*nbasis*nbasis;

        tF2.n = tF2.m = nbasis;
        tF2.c = new double*[nbasis];

        tF2.c2 = F2;
        for (int i=0; i<nbasis; ++i) tF2.c[i] = F2 + i*nbasis;
    }


    // print the stuff
    for (int n=0; n<nmat; ++n) {
        cout << "D" << n << endl << tD[n] << endl;
        cout << "F" << n << endl << tF[n] << endl << endl;
    }


    // dispose of the supporting structures properly
    for (int n=0; n<nmat; ++n) tF[n].c2 = NULL;
    for (int n=0; n<nmat; ++n) tD[n].c2 = NULL;

    delete[] tF;
    delete[] tD;
}


