/* -*-mode:c; c-style:bsd; c-basic-offset:4;indent-tabs-mode:nil; -*- */
#if !defined(GRID_GEN_H)
#define GRID_GEN_H 1

#include "general.h"

typedef struct GridGenAtom_     GridGenAtom;
typedef struct GridGenAtomGrid_ GridGenAtomGrid;
typedef struct GridGenMolGrid_  GridGenMolGrid;

typedef real (*GridGeneratingFunc)(real x, real y, real z, void* arg);

struct GridGenAtom_ {
    real x, y, z; /* coordinates of the atom */
    int icent;    /* number of atom in dalton common blocks */
    int Z;        /* its Z number            */
    int mult;     /* number of symmetry equivalent atoms of this type *
                   * instead of having separate atoms, we multiply the *
                   * grid weights by this value. */
};

typedef enum {
    GRID_PART_BECKE,  /* with  Bragg radii correction */
    GRID_PART_BECKE2, /* w/o   Bragg radii correction */
    GRID_PART_SSF
} GridGenPartScheme;

typedef enum {
    GRID_RAD_GC2,
    GRID_RAD_LMG
} GridGenQuad;

extern int dftgrid_adaptive;

struct RhoEvalData {
    DftGrid* grid;
    real* work;
    int*  lwork;
    real* dmat;
    real* dmagao;
};

real rho_grid_func(real x, real y, real z, void* arg);

GridGenAtom* grid_gen_atom_new(int* atom_cnt);
int grid_gen_save(const char* filename, GridGenMolGrid* mgrid);
int grid_gen_generate(const char* filename, int atom_cnt, 
                      const GridGenAtom* atom_arr, real threshold,
                      GridGeneratingFunc generating_function, void* data,
                      int minang, int maxang, real* work, int *lwork);
void grid_gen_set_part_scheme(GridGenPartScheme scheme);

typedef struct DftGridReader_ DftGridReader;

DftGridReader* grid_open(int nbast, real *work, int *lwork);

int
grid_getchunk_blocked(DftGridReader* rawgrid, int maxlen,
                      int *nblocks, int *shlblocks, 
                      real (*coor)[3], real *weight);

#define grid_getchunk_plain(r,m,coor,w) \
       (grid_getchunk_blocked((r),(m),NULL,NULL,(coor),(w)))
void grid_close(DftGridReader *rawgrid);

#endif /* !defined(GRID_GEN_H) */

 
