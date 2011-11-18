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
    GRID_PART_BECKE_CORR,  /* with  Bragg radii correction */
    GRID_PART_BECKE_ORIG,  /* w/o   Bragg radii correction */
    GRID_PART_SSF,
    GRID_PART_BLOCK
} GridGenPartScheme;

typedef enum {
    GRID_RAD_GC2,
    GRID_RAD_LMG
} GridGenQuad;

extern int dftgrid_adaptive;

struct RhoEvalData {
    DftGrid* grid;
    real* work;
    integer*  lwork;
    real* dmat;
    real* dmagao;
};

real rho_grid_func(real x, real y, real z, void* arg);

int grid_gen_save(const char* filename, GridGenMolGrid* mgrid);
int grid_gen_generate(const char* filename, integer atom_cnt, 
                      const GridGenAtom* atom_arr, real threshold,
                      GridGeneratingFunc generating_function, void* data,
                      int minang, int maxang, real* work, integer *lwork);
void grid_gen_set_part_scheme(GridGenPartScheme scheme);

typedef struct DftGridReader_ DftGridReader;

DftGridReader* grid_open(integer nbast, real *dmat, real *work, integer *lwork);
DftGridReader* grid_open_cmo(integer nbast, const real *cmo, 
                             real *work, integer *lwork);

int
grid_getchunk_blocked(DftGridReader* rawgrid, integer maxlen,
                      integer *nblocks, integer (*shlblocks)[2], 
                      real (*coor)[3], real *weight);

#define grid_getchunk_plain(r,m,coor,w) \
       (grid_getchunk_blocked((r),(m),NULL,NULL,(coor),(w)))
void grid_close(DftGridReader *rawgrid);

/* CARTESIAN GRID ROUTINES */
void do_cartesian_grid(int nbast, const real* dmat, DftGridReader* res);

#endif /* !defined(GRID_GEN_H) */

 
