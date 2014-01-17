/* -*-mode:c; c-style:bsd; c-basic-offset:4;indent-tabs-mode:nil; -*- */
#if !defined(GRID_GEN_H)
#define GRID_GEN_H 1

#include "lsdalton_general.h"

typedef struct GridGenAtom_     GridGenAtom;
typedef struct GridGenAtomGrid_ GridGenAtomGrid;
typedef struct GridGenMolGrid_  GridGenMolGrid;

typedef real (*GridGeneratingFunc)(real x, real y, real z, void* arg);

struct GridGenAtom_ {
    real x, y, z; /* coordinates of the atom */
    integer icent;    /* number of atom in dalton common blocks */
    integer Z;        /* its Z number            */
    integer mult;     /* number of symmetry equivalent atoms of this type *
                   * instead of having separate atoms, we multiply the *
                   * grid weights by this value. */
};

typedef enum {
    GRID_PART_BECKE_CORR,  /* Becke scheme with  Bragg radii correction */
    GRID_PART_BECKE_ORIG,  /* Becke scheme w/o   Bragg radii correction */
    GRID_PART_SSF,         /* SSF scheme */
    GRID_PART_BLOCK,       /* blocked Becke scheme with Bragg radii correction*/
    GRID_PART_BLOCK_SSF    /* blocked SSF scheme */
} GridGenPartScheme;

typedef enum {
    GRID_RAD_GC2,
    GRID_RAD_LMG,
    GRID_RAD_TUR
} GridGenQuad;

extern integer dftgrid_adaptive;

struct RhoEvalData {
    DftGrid* grid;
    real* work;
    integer*  lwork;
    real* dmat;
    real* dmagao;
};

void grid_gen_set_part_scheme(GridGenPartScheme scheme);

typedef struct DftGridReader_ DftGridReader;

DftGrid *dft_grid_new();

void dft_grid_destroy(DftGrid *grid);

DftGridReader* II_grid_get_reader(DftGrid *grid, integer nbast, integer natoms,
				  const real *X, 
                                  const real *Y, const real *Z,
				  const integer *Charge,integer *grdone,
				  const integer *NCENT, const integer *NHKT, 
				  const integer *NUCO, integer NHTYP, integer KMAX,
				  integer MXPRIM, const real *PRIEXP, 
				  const integer *NSTART, real *RSHEL, 
                                  integer *it, integer *maxnlen, integer TURBO,
                                  integer GC2, integer lupri);

DftGridReader* II_grid_open(integer nbast, real radint,
                            integer angmin,integer angint, integer hardness, 
                            integer prune, integer natoms, 
                            const real *X, const real *Y, const real *Z, 
                            const integer *Charge, integer *grdone, const integer *NCENT, 
                            const integer *NHKT, const integer *NUCO, integer NHTYP, 
                            integer KMAX,integer MXPRIM, const real *PRIEXP, 
                            const integer *NSTART, real *RSHEL, integer *it, 
                            integer *maxnlen, integer TURBO, integer GC2, integer lupri);

integer
grid_getchunk_blocked(DftGridReader* rawgrid, integer maxlen,
                      integer *nblocks, integer (*shlblocks)[2],
                      real (*coor)[3], real *weight);

#define grid_getchunk_plain(r,m,coor,w) \
       (grid_getchunk_blocked((r),(m),NULL,NULL,(coor),(w)))
void grid_close(DftGridReader *rawgrid);

/* CARTESIAN GRID ROUTINES */
/*void do_cartesian_grid(integer nbast, const real* dmat, DftGridReader* res);
 */
#endif /* !defined(GRID_GEN_H) */

 
