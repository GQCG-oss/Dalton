/* -*-mode:c; c-style:bsd; c-basic-offset:4;indent-tabs-mode:nil; -*- */
/* Zilvinas Rinkevicus: single atom integration with GC2 radial grid.
   Pawel Salek        : all the rest.
   2002.10.10-10.14
*/
/* define VAR_MPI if you want to use parallel grid generation 
 * for later parallel MPI calculation.

/* define USE_PTHREADS if you want to use multithreaded grid
 * generation.  It allows to take an advantage of SMP
 * architectures. It is currently implemented for BECKE, BECK2 and SSF
 * type grids (i.e those that do the work in the preprocessing phase.
 * USE_PTHREADS and VAR_MPI are do not conflict.
 */
#if !defined(SYS_AIX)
/* for an yet unknown reason, MT-grid generation has a terrible
   performance on regatta/AIX.
*/
/* #define USE_PTHREADS */
#endif /* !defined(SYS_AIX) */

#define __CVERSION__
#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas and pow() prototype
 * is not specified. */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif
#include <assert.h>
#include <limits.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>

#include "grid-gen.h"

/* common screening for cubes/cells */
static const real CELL_SIZE = 4.0;
/* the weight threshold:: a grid point will be ignored if has a weight
 * lower than the WEIGHT_THRESHOLD. This should probably depend on the
 * calculation's accuracy. */
static const real WEIGHT_THRESHOLD = 1e-17;


struct point {
    real x,y,z;
};

static int create_cubes(GridGenMolGrid *mg, const char* fname,
                        int point_cnt, real cell_size,
                        void* work, int worksz);

void*
dal_malloc_(size_t sz, const char *place, int line)
{
    void* res = malloc(sz);
    if(!res) {
        fprintf(stderr, "dal_malloc(sz=%d bytes) at %s (line %u) failed.\n",
                sz, place, line);
        exit(1);
    }
    return res;
}

extern void nucbas_(int*, real* , const int*);
extern void radlmg_(real*rad, real* wght, int *nr, real* raderr,
                    const int*maxrad, int *nucorb, real *aa, const int *);
extern void get_no_atoms_(int* atom_cnt);
extern void get_atom_by_icent_(const int* icent, real* charge, int *cnt,
                               int *mult, real *x, real *y, real *z);

/* generate a list of atoms from Dalton common block data, taking into
   account symmetries in the molecule.
   The list of atoms is accessed in two modes:
   - grid is generated only for symmetry-independent atoms...
   - ... but all atoms need to be taken into account when performing
   the space partitioning.
*/
GridGenAtom*
grid_gen_atom_new(int* atom_cnt)
{
    int nat = 0, icent, i, mult;
    real x[8], y[8], z[8], charge;
    GridGenAtom* atoms;

    get_no_atoms_(atom_cnt);
    atoms = dal_malloc(*atom_cnt*sizeof(GridGenAtom));
    
    icent = 0;
    do {
        int cnt;
        icent++;
        get_atom_by_icent_(&icent, &charge, &cnt, &mult, x,y,z);
        for(i=0; i<cnt; i++) {
            atoms[nat].x = x[i];
            atoms[nat].y = y[i];
            atoms[nat].z = z[i];
            atoms[nat].icent = icent-1;
            atoms[nat].Z = charge;
            atoms[nat].mult = mult;
            nat++;
        }
    } while(nat<*atom_cnt);
    return atoms;
}
/* ------------------------------------------------------------------- */
/* Internal data structures used by the grid generator.                */
/* ------------------------------------------------------------------- */
struct GridGenAtomGrid_ {
    int uniq_no;  /* number of the atom, before the symmetry multiplication */
    int Z;        /* frequently used. */
    int* leb_ang; /* leb_gen index associated with each radial point */
    real* rad;    /* radial points array. */
    real* wght;   /* radial points weights */
    int pnt;      /* number of radial points for this atom*/
};  

struct GridGenMolGrid_ {
    int atom_cnt;                  /* number of atoms grid is generated for*/
    const GridGenAtom* atom_coords;/* array with atom data */
    GridGenAtomGrid**   atom_grids;/* array of atom grids */
    real* rij;  /* triangular array of inverse distances between atoms */
    real* aij;  /* triangular array of of.... */
    FILE* fl;   /* the stream grid is being saved to */
    int total_points; /* total number of generated points */
    int off; /* thread number. 0 in serial. */
    int nt;  /* total number of threads. 1 in serial. */
};

/* work data - per thread */
struct GridGenWork_ {
    real* rj;   /* temp array for confocal elliptical coordinates */
    real* p_kg; /* unormalized weights evaluated for given gridpoint/atom */
    real* vec;  /* temporary buffer */
    real* x;  /* x coord */
    real* y;  /* y coord */
    real* z;  /* z coord */
    real* wg; /* weigths */ 
};
typedef struct GridGenWork_ GridGenWork;

void dzero_(real*arr, const int* cnt);
static int verbose=0;
/* include_selected partitioning points to selected partitioning scheme:
   Becke or SSF. The function is supposed to update weights wg.
*/
static void
include_partitioning_becke1(GridGenMolGrid* mg, int catom, int point_cnt,
                            GridGenWork *ggw, int idx, int verbose);
static void
include_partitioning_becke2(GridGenMolGrid* mg, int catom, int point_cnt,
                            GridGenWork *ggw, int idx, int verbose);
static void
include_partitioning_ssf  (GridGenMolGrid* mg, int catom, int point_cnt,
                           GridGenWork *ggw, int idx, int verbose);
static void gen_gc2_quad(GridGenAtomGrid* grid, real thrl, void *quad_data);
static void* gen_lmg_init(void);
static void gen_lmg_quad(GridGenAtomGrid* grid, real thrl, void *quad_data);
static void gen_lmg_free(void *quad_data);
static int block_partition_postprocess(GridGenMolGrid *mg, struct point *c,
                                       int point_cnt, const int *atom_nums,
                                       real (*coor)[3], real *w);


/* Trond Saue:
    The below data gives atomic radii in Angstroms and stems from table I of 
    J.C.Slater: "Atomic Radii in Crystals"
    J.Chem.Phys. 41(1964) 3199-3204
    Values for elements marked with an asterisk has been
    guessed/interpolated
*/
static const real bragg_radii[] = {
/*       H      He*   */
       0.35,  0.35,  
/*       Li     Be     B      C      N      O      F      Ne*  */
       1.45,  1.05,  0.85,  0.70,  0.65,  0.60,  0.50,  0.45,  
/*       Na     Mg     Al     Si     P      S      Cl     Ar*  */
       1.80,  1.50,  1.25,  1.10,  1.00,  1.00,  1.00,  1.00,  
/*      K      Ca     Sc     Ti     V      Cr     Mn     Fe     Co   */
       2.20,  1.80,  1.60,  1.40,  1.35,  1.40,  1.40,  1.40,  1.35,  
/*      Ni     Cu     Zn     Ga     Ge     As     Se     Br     Kr*  */
       1.35,  1.35,  1.35,  1.30,  1.25,  1.15,  1.15,  1.15,  1.10,  
/*      Rb     Sr     Y      Zr     Nb     Mo     Tc     Ru     Rh  */
       2.35,  2.00,  1.80,  1.55,  1.45,  1.45,  1.35,  1.30,  1.35,  
/*      Pd     Ag     Cd     In     Sn     Sb     Te     I      Xe*  */
       1.40,  1.60,  1.55,  1.55,  1.45,  1.45,  1.40,  1.40,  1.40,  
/*      Cs     Ba     La      */
       2.60,  2.15,  1.95,  
/*      Ce     Pr     Nd     Pm     Sm     Eu     Gd  */
       1.85,  1.85,  1.85,  1.85,  1.85,  1.85,  1.80,  
/*      Tb     Dy     Ho     Er     Tm     Yb     Lu  */
       1.75,  1.75,  1.75,  1.75,  1.75,  1.75,  1.75,  
/*      Hf     Ta     W      Re     Os     Ir     Pt     Au     Hg  */
       1.55,  1.45,  1.35,  1.30,  1.30,  1.35,  1.35,  1.35,  1.50,  
/*      Tl     Pb*    Bi     Po     At*    Rn*  */
       1.90,  1.75,  1.60,  1.90,  1.50,  1.50,  
/*      Fr*    Ra     Ac       */
       2.15,  2.15,  1.95,  
/*     rad(U): 1.75 --> 1.37D0  */
/*      Th     Pa     U      Np     Pu     Am     Cm*       */
       1.80,  1.80,  1.37,  1.75,  1.75,  1.75,  1.75,  
/*      Bk*    Cf*    Es*    Fm*    Md*    No*    Lw*  */
       1.75,  1.75,  1.75,  1.75,  1.75,  1.75,  1.75
    };
    

/* ------------------------------------------------------------------- */
/* The grid generator configuration data with defaults.                */
/* ------------------------------------------------------------------- */
/* include_selected_partitioning scales the atomic grid centered at
 * atom no 'catom'.
 */
struct partitioning_scheme_t {
    char *name;
    void (*preprocess)
    (GridGenMolGrid* mg, int catom, int point_cnt,
     GridGenWork *ggw, int idx, int verbose);
    int (*postprocess)
    (GridGenMolGrid* mg, struct point *c, 
     int point_cnt, const int *atom_nums,
     real (*coor)[3], real *w);
};
static struct partitioning_scheme_t partitioning_becke1 =
{ "Becke partitioning with atomic radius correction", 
  include_partitioning_becke1, NULL };
static struct partitioning_scheme_t partitioning_becke2 =
{ "Original Becke partitioning", include_partitioning_becke2 };
static struct partitioning_scheme_t partitioning_ssf =
{ "SSF partitioning", include_partitioning_ssf, NULL };
static struct partitioning_scheme_t partitioning_block =
{ "Blocked partitioning for large molecules/parallel calc:s.",
  NULL, block_partition_postprocess };

static struct partitioning_scheme_t *selected_partitioning =
&partitioning_becke1;

struct radial_scheme_t {
    char *name;
    void *(*quad_init)(void);
    void (*quad_gen)  (GridGenAtomGrid*, real thrl, void *quad_data);
    void (*quad_free) (void *quad_data);
};
static struct radial_scheme_t quad_lmg = { 
    "LMG scheme",
    gen_lmg_init,
    gen_lmg_quad,
    gen_lmg_free
};
static struct radial_scheme_t quad_gc2 = { 
    "Gauss-Chebychev scheme of second kind",
    NULL,
    gen_gc2_quad,
    NULL
};

static struct radial_scheme_t *radial_quad = &quad_lmg;

/* ------------------------------------------------------------------- */
/*              The angular grid generation interface                  */
/* ------------------------------------------------------------------- */

/* routines for generation of Lebedev grid */
void ld0026_(real* x, real* y, real* z, real* w, int *pnt);
void ld0038_(real* x, real* y, real* z, real* w, int *pnt);
void ld0050_(real* x, real* y, real* z, real* w, int *pnt);
void ld0074_(real* x, real* y, real* z, real* w, int *pnt);
void ld0086_(real* x, real* y, real* z, real* w, int *pnt);
void ld0110_(real* x, real* y, real* z, real* w, int *pnt);
void ld0146_(real* x, real* y, real* z, real* w, int *pnt);
void ld0170_(real* x, real* y, real* z, real* w, int *pnt);
void ld0194_(real* x, real* y, real* z, real* w, int *pnt);
void ld0230_(real* x, real* y, real* z, real* w, int *pnt);
void ld0266_(real* x, real* y, real* z, real* w, int *pnt);
void ld0302_(real* x, real* y, real* z, real* w, int *pnt);
void ld0350_(real* x, real* y, real* z, real* w, int *pnt);
void ld0434_(real* x, real* y, real* z, real* w, int *pnt);
void ld0590_(real* x, real* y, real* z, real* w, int *pnt);
void ld0770_(real* x, real* y, real* z, real* w, int *pnt);
void ld0974_(real* x, real* y, real* z, real* w, int *pnt);
void ld1202_(real* x, real* y, real* z, real* w, int *pnt);
void ld1454_(real* x, real* y, real* z, real* w, int *pnt);

struct leb_gen_{
  int point_cnt, poly_type;
  void (*func)(real* x, real* y, real* z, real* w, int *pnt);
};

/* some of the point type values were guessed but since it was only used
   for ANGMIN -> quadrature mapping type, the accuracy is not really relevant.
   Anyone, who would rely on particular mapping here would be crazy.
   The true values can be found in literature.
*/

struct leb_gen_ leb_gen[] = {
 /* {  26,  0, ld0026_}, CONSIDER disconnecting */
 {  38,  0, ld0038_},
 {  50, 9, ld0050_},
 /* {  74, 12, ld0074_},* point type guessed; consider DISCONNECTING */
 {  86, 11, ld0086_},
 { 110, 15, ld0110_},
 { 146, 17, ld0146_},
 { 170, 19, ld0170_},
 { 194, 21, ld0194_},
 /* { 230, 22, ld0230_}, point type guessed; consider DISCONNECTING */
 /* { 266, 22, ld0266_}, point type guessed; consider DISCONNECTING */
 { 302, 23, ld0302_},
 { 350, 29, ld0350_},
 { 434, 31, ld0434_}, 
 { 590, 35, ld0590_},
 { 770, 41, ld0770_},
 { 974, 47, ld0974_},
 {1202, 53, ld1202_},
 {1454, 59, ld1454_},
};  

/* get_leb_from_point returns index to leb_gen array that points to
   the entry having at least iang points.
*/
int get_leb_from_point(int iang);

/* get_leb_from_type returns index to leb_gen array that points to
   the entry contructing polyhedra of at least given order.
*/
int get_leb_from_type(int poly_type);
/* get_leb_from_point returns index to leb_gen array that points to
   the entry having at least iang points.
*/
int
get_leb_from_point(int iang)
{
    int i;
    for(i=ELEMENTS(leb_gen)-1; i>=0; i--)
        if(leb_gen[i].point_cnt<iang) return i;
    return 0;
}

/* get_leb_from_type returns index to leb_gen array that points to
   the entry contructing polyhedra of at least given order.
*/
int
get_leb_from_type(int poly_type)
{
    int i;
    for(i=ELEMENTS(leb_gen)-1; i>=0; i--)
        if(leb_gen[i].poly_type<poly_type)
           return i;
    return 0;
}

/* ------------------------------------------------------------------- */
/* the integrator/grid generator itself                                */
/* ------------------------------------------------------------------- */
void
grid_gen_set_part_scheme(GridGenPartScheme scheme)
{
    switch(scheme) {
    case GRID_PART_BECKE:
        selected_partitioning = &partitioning_becke1;
        break;
    case GRID_PART_BECKE2:
        selected_partitioning = &partitioning_becke2;
        break;
    case GRID_PART_SSF: 
        selected_partitioning = &partitioning_ssf;
        break;
    case GRID_PART_BLOCK:
        selected_partitioning = &partitioning_block;
        break;
    }
}

/* grid_gen_set_rad_quad:
   set radial quadrature.
*/
void
grid_gen_set_rad_quad(GridGenQuad scheme)
{
    switch(scheme) {
    default:
    case GRID_RAD_GC2: radial_quad = &quad_gc2; break;
    case GRID_RAD_LMG: radial_quad = &quad_lmg; break;
    }
}
                                                           
/* max number of points per atom per angular shell */
#define MAX_PT_PER_SHELL (leb_gen[ELEMENTS(leb_gen)-1].point_cnt)

/* include_partitioning_weights:
   scale weights
*/
#define IDX(arr,i,j) arr[(i)*MAX_PT_PER_SHELL +(j)]
static void
gridgen_compute_rjs(GridGenWork *ggw, GridGenMolGrid* mg,
                    int point_cnt, int idx)
{
    int atno, ptno;

    for(atno=0; atno<mg->atom_cnt; atno++) {
        for(ptno=0; ptno<point_cnt; ptno++) {
            real dx = mg->atom_coords[atno].x - ggw->x[idx+ptno];
            real dy = mg->atom_coords[atno].y - ggw->y[idx+ptno];
            real dz = mg->atom_coords[atno].z - ggw->z[idx+ptno];
      
            assert(atno>=0 && atno<mg->atom_cnt);
            assert(ptno>=0 && ptno<MAX_PT_PER_SHELL);
            ggw->IDX(rj,atno,ptno)   = sqrt(dx*dx + dy*dy + dz*dz);
            ggw->IDX(p_kg,atno,ptno) = 1.0; 
        }
    }
}
/* include_partitioning_becke:
   multiply current sphere/shell generated for specified atom by space
   partitioning weights as described by Becke.
   The computed weights are somehow redundant: unnormalized
   weights p_kg are computed for all atoms, althougth only one is 
   needed at this stage. The algorithm should probably be restructured
   to avoid this.
*/
#define HARDNESS1 11
static void
include_partitioning_becke1(GridGenMolGrid* mg, int atom, int point_cnt,
                            GridGenWork *ggw, int idx, int verbose)
{
    int atno, atno2, ptno, h, isign=-1;
    real mu, mu2, g_mu, apasc;
    real xpasc[HARDNESS1], facult[HARDNESS1];
    /* compute confocal ellipical coordinates for the batch of points to
     * be processed. */

    if(mg->atom_cnt==1)
        return;
    facult[0] = 1;
    for(h=1; h<HARDNESS1; h++)
        facult[h] = facult[h-1]*h;

    for(h=HARDNESS1-1; h>=0; h--) {
        isign = -isign;
        xpasc[h] = isign*facult[HARDNESS1-1]/
            ((2*h+1)*facult[h]*facult[HARDNESS1-1-h]);
    }
    xpasc[0] = 1;
    apasc = 0;
    for(h=0; h<HARDNESS1; h++) apasc += xpasc[h];
    apasc = 0.5/apasc;

    gridgen_compute_rjs(ggw, mg, point_cnt, idx);

    if(verbose)printf("Computing cell functions for atom %d\n", atom);
    for(atno=1; atno<mg->atom_cnt; atno++) {
        for(atno2=0; atno2<atno; atno2++) {
            real dist = mg->rij[atno2+(atno*(atno-1))/2];
            for(ptno=0; ptno<point_cnt; ptno++) {
                mu    =(ggw->IDX(rj,atno,ptno)-ggw->IDX(rj,atno2,ptno))*dist;
                mu2 = mu*mu;
                g_mu = 0;
                for(h=0; h<HARDNESS1; h++) {
                    g_mu += xpasc[h]*mu;
                    mu *= mu2;
                }
                ggw->IDX(p_kg,atno,ptno)  *= 0.5-apasc*g_mu;
                ggw->IDX(p_kg,atno2,ptno) *= 0.5+apasc*g_mu;
            }
        }
    }
    /* compute weight normalization factors */
    dzero_(ggw->vec, &point_cnt);
    for(atno=0; atno<mg->atom_cnt; atno++)
        for(ptno=0; ptno<point_cnt; ptno++)
            ggw->vec[ptno] += ggw->IDX(p_kg,atno,ptno); 

    /*
     *   Apply the computed weights
     */
    for(ptno=0; ptno<point_cnt; ptno++)
        ggw->wg[idx+ptno] *= ggw->IDX(p_kg,atom,ptno)/ggw->vec[ptno];
}

#define HARDNESS2 3
static void
include_partitioning_becke2(GridGenMolGrid* mg, int atom, int point_cnt,
                            GridGenWork *ggw, int idx, int verbose)
{
    int atno, atno2, ptno, h;
    real mu, g_mu;
    /* compute confocal ellipical coordinates for the batch of points to
     * be processed. */

    if(mg->atom_cnt==1)
        return;
    gridgen_compute_rjs(ggw, mg, point_cnt, idx);

    if(verbose)printf("Computing cell functions for atom %d\n", atom);
    for(atno=1; atno<mg->atom_cnt; atno++) {
        for(atno2=0; atno2<atno; atno2++) {
            real dist = mg->rij[atno2+(atno*(atno-1))/2];
            real bfac = mg->aij[atno2+(atno*(atno-1))/2];
            for(ptno=0; ptno<point_cnt; ptno++) {
                mu    =(ggw->IDX(rj,atno,ptno)-ggw->IDX(rj,atno2,ptno))*dist;
                mu   += bfac*(1-mu*mu);
                if(mu<-1) mu= -1; /* numerical error correction, needed? */
                if(mu>1)  mu=  1; /* numerical error correction, needed? */
                for(g_mu = mu, h=0; h<HARDNESS2; h++)
                    g_mu = 0.5*g_mu*(3.0-g_mu*g_mu);
                ggw->IDX(p_kg,atno,ptno)  *= 0.5*(1.0-g_mu);
                ggw->IDX(p_kg,atno2,ptno) *= 0.5*(1.0+g_mu);
            }
        }
    }
    /* compute weight normalization factors */
    dzero_(ggw->vec, &point_cnt);
    for(atno=0; atno<mg->atom_cnt; atno++)
        for(ptno=0; ptno<point_cnt; ptno++)
            ggw->vec[ptno] += ggw->IDX(p_kg,atno,ptno); 

    /*
     *   Apply the computed weights
     */
    for(ptno=0; ptno<point_cnt; ptno++)
        ggw->wg[idx+ptno] *= ggw->IDX(p_kg,atom,ptno)/ggw->vec[ptno];
}

/* include_partitioning_ssf:
   multiply current sphere/shell by space partitioning weights as
   described in SSF article.
   FIXME: compute relevant atoms only once.
*/
static void
include_partitioning_ssf(GridGenMolGrid* mg, int atom, int point_cnt,
                         GridGenWork *ggw, int idx, int verbose)
{
    static const real SSF_CUTOFF = 0.64;
    int atno, ptno;
    real mu, g_mu, pr, rx, ry, rz;
    int * relevant_atoms = dal_malloc(sizeof(int)*mg->atom_cnt);
    int atomi=0, ati, ati2, rel_atom_cnt = 0;

    rx = ggw->x[idx] - mg->atom_coords[atom].x;
    ry = ggw->y[idx] - mg->atom_coords[atom].y;
    rz = ggw->z[idx] - mg->atom_coords[atom].z;
    pr = sqrt(rx*rx+ry*ry+rz*rz); /* radius of the sphere */

    /* find relevant atoms for 'atom':
     * separate loops for atoms below and above */

    for(atno=0; atno<atom; atno++) {
        real dist = mg->rij[atno+(atom*(atom-1))/2];
        if(2*pr*dist>1-SSF_CUTOFF)
            relevant_atoms[rel_atom_cnt++] = atno;
    }
    atomi = rel_atom_cnt;
    relevant_atoms[rel_atom_cnt++] = atom;
    for(atno=atom+1; atno<mg->atom_cnt; atno++) {
        real dist = mg->rij[atom+(atno*(atno-1))/2];
        if(2*pr*dist>1-SSF_CUTOFF)
            relevant_atoms[rel_atom_cnt++] = atno;
    }

    if(rel_atom_cnt==1) {
        free(relevant_atoms);
        return;
    }
           
    /* compute confocal ellipical coordinates for the batch of points to
     * be processed. */
    for(ati=0; ati<rel_atom_cnt; ati++) {
        atno = relevant_atoms[ati];
        for(ptno=0; ptno<point_cnt; ptno++) {
            real dx = mg->atom_coords[atno].x-ggw->x[ptno+idx];
            real dy = mg->atom_coords[atno].y-ggw->y[ptno+idx];
            real dz = mg->atom_coords[atno].z-ggw->z[ptno+idx];
      
            assert(atno>=0 && atno<mg->atom_cnt);
            assert(ptno>=0 && ptno<MAX_PT_PER_SHELL);
            ggw->IDX(rj,ati,ptno)  = sqrt(dx*dx + dy*dy + dz*dz);
            ggw->IDX(p_kg,ati,ptno) = 1.0; 
        }
    }

    if(verbose)printf("Computing cell functions for atom %d\n", atom);
    for(ati=1; ati<rel_atom_cnt; ati++) {
        atno = relevant_atoms[ati];
        for(ati2=0; ati2<ati; ati2++) {
            int atno2 = relevant_atoms[ati2];
            real mu2;
            real dist = mg->rij[atno2+(atno*(atno-1))/2];
            for(ptno=0; ptno<point_cnt; ptno++) {
                mu    =(ggw->IDX(rj,ati,ptno)-ggw->IDX(rj,ati2,ptno))*dist;
                if(mu<=-SSF_CUTOFF) 
                    g_mu = -1;
                else if(mu>=SSF_CUTOFF)
                    g_mu = 1;
                else {
                    mu /= SSF_CUTOFF; mu2=mu*mu;
                    g_mu  = 0.0625*mu*(35+mu2*(-35+mu2*(21-5*mu2)));
                }
                ggw->IDX(p_kg,ati,ptno)  *= 0.5*(1.0-g_mu);
                ggw->IDX(p_kg,ati2,ptno) *= 0.5*(1.0+g_mu);
            }
        }
    }
    /* compute weight normalization factors */
    dzero_(ggw->vec, &point_cnt);
    for(ati=0; ati<rel_atom_cnt; ati++)
        for(ptno=0; ptno<point_cnt; ptno++)
            ggw->vec[ptno] += ggw->IDX(p_kg,ati,ptno); 

    /*
     *   Apply the computed weights
     */
    for(ptno=0; ptno<point_cnt; ptno++)
        ggw->wg[idx+ptno] *= ggw->IDX(p_kg,atomi,ptno)/ggw->vec[ptno];

    free(relevant_atoms);
}

/** include_block_partitioning transforms a set of POINT_CNT points.
*/
static int cmpint(const void *a, const void *b)
{ return *((int*)a)-*((int*)b); }
typedef struct {
    real (*coor)[3];
    real *rj;
    real *p_kg;
    real *vec;
    int LDA; /* leading dimension of rj and p_kg */
} GGBlockWork;

static void
block_work_init(GGBlockWork* ggw, real (*coor)[3],
                int atom_cnt, int point_cnt)
{
    ggw->coor = coor;
    ggw->LDA  = point_cnt;
    ggw->rj   = dal_malloc(atom_cnt*point_cnt*sizeof(real));
    ggw->p_kg = dal_malloc(atom_cnt*point_cnt*sizeof(real));
    ggw->vec  = dal_malloc(point_cnt*sizeof(real));
}

static void
block_work_release(GGBlockWork* ggw)
{
    free(ggw->rj);
    free(ggw->p_kg);
    free(ggw->vec);
}

static void
block_compute_rjs(GGBlockWork *ggw, GridGenMolGrid* mg,
                  int rel_at_cnt, int *relevant_atoms,
                  int point_cnt)
{
    int atno, ptno, i;

    for(i=0; i<rel_at_cnt; i++) {
        int atno = relevant_atoms[i];
        for(ptno=0; ptno<point_cnt; ptno++) {
            real dx = mg->atom_coords[atno].x - ggw->coor[ptno][0];
            real dy = mg->atom_coords[atno].y - ggw->coor[ptno][1];
            real dz = mg->atom_coords[atno].z - ggw->coor[ptno][2];

            ggw->rj  [i*ggw->LDA+ptno] = sqrt(dx*dx + dy*dy + dz*dz);
            ggw->p_kg[i*ggw->LDA+ptno] = 1.0;
        }
    }
}

static int
block_partition_postprocess(GridGenMolGrid *mg, struct point *c,
                            int point_cnt, const int *atom_nums,
                            real (*coor)[3], real *w)
{
    int atno, atno2, ptno, h, isign=-1, i, j;
    real mu, mu2, g_mu, apasc;
    real xpasc[HARDNESS1], facult[HARDNESS1];
    /* compute confocal ellipical coordinates for the batch of points to
     * be processed. */
    int *relevant_atoms = dal_malloc(mg->atom_cnt*sizeof(int));
    /* map2r: map from all to relevant_atoms array */
    int *map2r = dal_malloc(mg->atom_cnt*sizeof(int));
    int last_atom = -1, uniq_atoms, curr_atom_pos;
    GGBlockWork ggw;
    int dest;

    /* we fidn first atoms that relevant for this cell.
     * We do it by linear search which will scale as N^2
     * but we can live with that for now.
     */
    uniq_atoms = 0;
    for(i=0; i<mg->atom_cnt; i++) map2r[i] = -1;
    for(atno=0; atno<mg->atom_cnt; atno++) {
        GridGenAtomGrid *ag = mg->atom_grids[atno];
        real dx = mg->atom_coords[atno].x - c->x;
        real dy = mg->atom_coords[atno].y - c->y;
        real dz = mg->atom_coords[atno].z - c->z;
        real dist2 = dx*dx + dy*dy + dz*dz;
        real r = ag->rad[ag->pnt-1]+CELL_SIZE;
        if(r*r>dist2) {
            map2r[atno] = uniq_atoms;
            relevant_atoms[uniq_atoms++] = atno;

        }
    }
    if(uniq_atoms<=1) { /* 0 cannot happen and 1 - no partitioning. */
        free(relevant_atoms);
        free(map2r);
        return point_cnt;
    }
    block_work_init(&ggw, coor, uniq_atoms, point_cnt);
    block_compute_rjs(&ggw, mg, uniq_atoms, relevant_atoms,
                      point_cnt);

    facult[0] = 1;
    for(h=1; h<HARDNESS1; h++)
        facult[h] = facult[h-1]*h;

    for(h=HARDNESS1-1; h>=0; h--) {
        isign = -isign;
        xpasc[h] = isign*facult[HARDNESS1-1]/
            ((2*h+1)*facult[h]*facult[HARDNESS1-1-h]);
    }
    xpasc[0] = 1;
    apasc = 0;
    for(h=0; h<HARDNESS1; h++) apasc += xpasc[h];
    apasc = 0.5/apasc;

    for(i=1; i<uniq_atoms; i++) {
        int atno = relevant_atoms[i];
        for(j=0; j<i; j++) {
            real dist, bfac;
            int atno2 = relevant_atoms[j];
            dist = mg->rij[atno2+(atno*(atno-1))/2];
            bfac = mg->aij[atno2+(atno*(atno-1))/2];
            for(ptno=0; ptno<point_cnt; ptno++) {
                mu =(ggw.rj[ggw.LDA*i+ptno]-
                     ggw.rj[ggw.LDA*j+ptno])*dist;
                mu += bfac*(1-mu*mu);
                mu2 = mu*mu;
                g_mu = 0;
                for(h=0; h<HARDNESS1; h++) {
                    g_mu += xpasc[h]*mu;
                    mu *= mu2;
                }
                ggw.p_kg[ggw.LDA*i+ptno] *= 0.5-apasc*g_mu;
                ggw.p_kg[ggw.LDA*j+ptno] *= 0.5+apasc*g_mu;
            }
        }
    }
    /* compute weight normalization factors */
    dzero_(ggw.vec, &point_cnt);
    for(i=0; i<uniq_atoms; i++) {
        for(ptno=0; ptno<point_cnt; ptno++)
            ggw.vec[ptno] += ggw.p_kg[ggw.LDA*i+ptno];
    }

    /*
     * Apply the computed weights removing at the same time
     * points with low weight.
     */
    for(dest=0, ptno=0; ptno<point_cnt; ptno++) {
        int atom = map2r[atom_nums[ptno]];
        w[dest] =
            w[ptno]*ggw.p_kg[ggw.LDA*atom+ptno]/ggw.vec[ptno];
        coor[dest][0] = coor[ptno][0];
        coor[dest][1] = coor[ptno][1];
        coor[dest][2] = coor[ptno][2];
        if(w[dest] >= WEIGHT_THRESHOLD) dest++;
    }
    free(relevant_atoms);
    free(map2r);
    block_work_release(&ggw);
    return dest;
}

/* ===================================================================
 *             GRID MEMORY HANDLING AND GENERAL GLUE.
 * =================================================================== */

/* atom grid init */
GridGenAtomGrid*
grid_gen_agrid_new(int uniq_no, int Z)
{
    GridGenAtomGrid* grid = dal_malloc(sizeof(GridGenAtomGrid));
    grid->uniq_no = uniq_no;
    grid->Z       = Z;
    grid->pnt     = 0;
    grid->leb_ang = NULL;
    grid->rad     = NULL;
    grid->wght    = NULL;
    return grid;
}

static void
grid_gen_agrid_set_size(GridGenAtomGrid* grid, unsigned cnt)
{
    grid->pnt = cnt;
    assert(grid->pnt>0);
    if(grid->leb_ang) free(grid->leb_ang);
    if(grid->rad)     free(grid->rad);
    if(grid->wght)    free(grid->wght);
    grid->leb_ang = calloc(grid->pnt, sizeof(real));
    grid->rad     = calloc(grid->pnt, sizeof(real));
    grid->wght    = calloc(grid->pnt, sizeof(real));
}


/* destroy atomic grid */
static void
grid_gen_agrid_free(GridGenAtomGrid* grid)
{
    free(grid->leb_ang);
    free(grid->rad);
    free(grid->wght);
    free(grid);
}

static GridGenMolGrid*
mol_grid_new(int atom_cnt, const GridGenAtom* atoms)
{
    int i, j, index;
    GridGenMolGrid* mg = dal_malloc(sizeof(GridGenMolGrid));
    mg->atom_cnt    = atom_cnt;
    mg->atom_coords = atoms;
    mg->atom_grids  = dal_malloc(atom_cnt*sizeof(GridGenAtomGrid*));
    mg->rij  = dal_malloc(sizeof(real)*(atom_cnt*(atom_cnt-1))/2 );
    mg->aij  = dal_malloc(sizeof(real)*(atom_cnt*(atom_cnt-1))/2 );
    index =0;
    for(i=0; i<atom_cnt; i++) {
        mg->atom_grids[i] = grid_gen_agrid_new(atoms[i].icent, atoms[i].Z);
        for(j=0; j<i; j++) {
            real chi = bragg_radii[atoms[i].Z-1]/bragg_radii[atoms[j].Z-1];
            real temp = (chi-1)/(chi+1);
            real dx = atoms[i].x - atoms[j].x;
            real dy = atoms[i].y - atoms[j].y;
            real dz = atoms[i].z - atoms[j].z;
            temp = temp/(temp*temp-1);
            if(temp>0.5) temp = 0.5;
            else if(temp<-0.5) temp = -0.5;
            mg->rij[index] = 1.0/sqrt(dx*dx + dy*dy + dz*dz);
            mg->aij[index++] = temp;
        }
    }
    return mg;
}

void
mol_grid_free(GridGenMolGrid* mg)
{
    int i;
    free(mg->rij);
    free(mg->aij);
    for(i=0; i< mg->atom_cnt; i++)
        grid_gen_agrid_free(mg->atom_grids[i]);
    free(mg->atom_grids);
}

static void
grid_gen_work_init(GridGenWork* ggw, GridGenMolGrid* mg, int mxshells)
{
    ggw->rj   = dal_malloc(sizeof(real)*mg->atom_cnt*MAX_PT_PER_SHELL);
    ggw->p_kg = dal_malloc(sizeof(real)*mg->atom_cnt*MAX_PT_PER_SHELL);
    ggw->vec  = dal_malloc(sizeof(real)*mg->atom_cnt*MAX_PT_PER_SHELL);
    ggw->x    = dal_malloc(MAX_PT_PER_SHELL*mxshells*sizeof(real));
    ggw->y    = dal_malloc(MAX_PT_PER_SHELL*mxshells*sizeof(real));
    ggw->z    = dal_malloc(MAX_PT_PER_SHELL*mxshells*sizeof(real));
    ggw->wg   = dal_malloc(MAX_PT_PER_SHELL*mxshells*sizeof(real));  
}

static void
grid_gen_work_release(GridGenWork* ggw)
{
    free(ggw->rj);
    free(ggw->p_kg);
    free(ggw->vec);
    free(ggw->x);
    free(ggw->y);
    free(ggw->z);
    free(ggw->wg); 
}
            
/* set_radial_grid:
   precompute radial grid for all the atoms in the molecule.
*/
static void
set_radial_grid(GridGenMolGrid* grd, real thrl)
{
    int atom;
    void *data = NULL; 
    if(radial_quad->quad_init)
        data = radial_quad->quad_init();
    for(atom=0; atom<grd->atom_cnt; atom++) {
        radial_quad->quad_gen(grd->atom_grids[atom], thrl, data);
    }
    if(radial_quad->quad_free)
        radial_quad->quad_free(data);
}

/* ===================================================================
 *             RADIAL QUADRATURES
 * the quadratore has to fill in grid->pnt with number of points
 * and set grid->rad.
 * =================================================================== */

/* gc2_rad_cnt:
 * determinates number of radial points to be used for
 * Gauss-Chebyshev quadrature of second kind needed to integrate atom
 * of specified Z number to specified threshold thrl.
 * the ordinary weights are multiplied by r^2
 */
static int
gc2_rad_cnt(int Z, real thrl)
{
    static const int MIN_RAD_PT = 20;
    int ta=1, ri;
    real nr=-5.0*(3*log10(thrl)-ta+8);
    ri = rint(nr); 
    return ri>MIN_RAD_PT ? ri : MIN_RAD_PT;
}

/* gc2_quad:
   Gauss-Chebyshev quadrature of second kind that we use to generate
   the radial grid. Requested number of points is in grid->pnt.
   The grid->rad and grid->wght arrays are filled in.
*/
static void 
gen_gc2_quad(GridGenAtomGrid* grid, real thrl, void *quad_data)
{
    /* constants */
    static const real pi_2 = 2.0/M_PI;  
    static const real sfac = 2.0/3.0;
    const real rfac = 1.0/log(2.0);
    real n_one, n_pi, wfac;
    /* variables */
    real x = 0.0, angl = 0.0, w = 0.0;
    int i;

    grid_gen_agrid_set_size(grid, gc2_rad_cnt(grid->Z,thrl));
    n_one = grid->pnt+1.0;
    n_pi  = M_PI/n_one;
    wfac = 16.0/(3*n_one);
    /* radial points */ 
    for (i=0; i<grid->pnt; i++) {
        real sinangl, sinangl2;
        x = (grid->pnt-1-2*i)/n_one;
        angl = n_pi*(i+1);
        sinangl = sin(angl); 
        sinangl2 = sinangl*sinangl;
        x += pi_2*(1.0+sfac*sinangl2)*cos(angl)*sinangl;
        grid->rad[i] = rfac*log(2.0/(1-x));
        w = wfac*sinangl2*sinangl2;
        grid->wght[i] = w*rfac/(1.0-x)*grid->rad[i]*grid->rad[i];
        /* transformation factor accumulated in weight */
    }
}
/* gen_lmg_quad:
 *  As proposed by Roland Lindh, Per-Aake Malmqvist and Laura
 *  Gagliardi. */

struct lmg_data {
    int  *nucorb;
    real *aa;
    int maxl;
};
void get_maxl_nucind_(int *maxl, int*nucind);

static void*
gen_lmg_init(void)
{
    int nucind;
    struct lmg_data *lmg = dal_malloc(sizeof(struct lmg_data));

    get_maxl_nucind_(&lmg->maxl, &nucind);
    lmg->nucorb = malloc(2*lmg->maxl*nucind*sizeof(int));
    lmg->aa     = calloc(4*lmg->maxl*nucind,sizeof(real));
    if(!lmg->nucorb|| !lmg->aa) {
        fprintf(stderr,"no enough memory. in gen_lmg_init.\n");
        exit(1);
    }
    nucbas_(lmg->nucorb, lmg->aa, &ONEI); 
    return lmg;
}

static void
gen_lmg_quad(GridGenAtomGrid* grid, real thrl, void *quad_data)
{
    static const int MAXRAD = 2000;
    struct lmg_data *lmg = (struct lmg_data*)quad_data;

    grid_gen_agrid_set_size(grid, MAXRAD); 
    radlmg_(grid->rad, grid->wght, &grid->pnt, &thrl, &MAXRAD,
            lmg->nucorb+2*grid->uniq_no*lmg->maxl,
            lmg->aa    +4*grid->uniq_no*lmg->maxl,
            &ZEROI);

    grid->rad  = realloc(grid->rad,  grid->pnt*sizeof(real));
    grid->wght = realloc(grid->wght, grid->pnt*sizeof(real));
}

static void
gen_lmg_free(void *quad_data)
{
    struct lmg_data *lmg = (struct lmg_data*)quad_data;
    
    free(lmg->nucorb);
    free(lmg->aa);
    free(lmg);
}


/* set_ang_fixed:
   computes the angular points for a set of radial points
   obtained from radial integration scheme.
   The only thing altered is grid->leb_ang vector.
*/
static void
set_ang_fixed(GridGenMolGrid* mgrid, real thrl, int minang, int maxang)
{
    static const real BOHR = 0.529177249;
    int atom, i = 0;

    for(atom=0; atom<mgrid->atom_cnt; atom++) {
        GridGenAtomGrid* grid = mgrid->atom_grids[atom];
        real rbragg = bragg_radii[grid->Z-1]/(5.0*BOHR);
        int current_ang = maxang;
        for (i=0; i<grid->pnt; i++) {
            if(grid->rad[i]<rbragg) {
                /* prune */
                int iang = 
                    (double)leb_gen[maxang].point_cnt*grid->rad[i]/rbragg;
                current_ang = get_leb_from_point(iang);
                if(current_ang<minang) current_ang = minang;
            } /* else current_ang = maxang; */
            grid->leb_ang[i] = current_ang;
        }
    }
}

static int
compress_grid(int point_cnt, real* x, real* y, real *z, real* wg)
{
    int i, dest=0;

    for(i=0; i<point_cnt; i++) {
        x [dest] = x [i];
        y [dest] = y [i];
        z [dest] = z [i];
        wg[dest] = wg[i];
        if(fabs(wg[i]) >WEIGHT_THRESHOLD)
            dest++;
/*
        else fort_print("skipping [%7.4f,%7.4f,%7.4f]: %g",
        x[i], y[i], z[i], wg[i]); */
    }
    /* printf("Compression factor: %f\n", 
       (float)(point_cnt-dest)/(float)point_cnt); */
    return dest;
}
#if defined(USE_PTHREADS)
#include <pthread.h>
pthread_mutex_t grid_mutex = PTHREAD_MUTEX_INITIALIZER;
static void
mol_grid_lock(GridGenMolGrid* mgrid)
{
    pthread_mutex_lock(&grid_mutex);
}
static void
mol_grid_unlock(GridGenMolGrid* mgrid)
{
    pthread_mutex_unlock(&grid_mutex);
}
static int
mol_grid_get_num_threads(void)
{ return 4; }
#else
#define mol_grid_lock(mgrid)
#define mol_grid_unlock(mgrid)
#endif


static void*
grid_gen_save_temp_stream(GridGenMolGrid* mgrid)
{
    int atom, ptno, j, idx, cnt;
    int tpt, off = mgrid->off, done;
    GridGenWork ggw;

    mol_grid_unlock(mgrid); /* off has been copied, can unlock */

    j=0;
    for(atom=off; atom<mgrid->atom_cnt; atom+=mgrid->nt) {
        if(mgrid->atom_grids[atom]->pnt>j) 
            j = mgrid->atom_grids[atom]->pnt;
    }
    if(j==0) return NULL; /* this thread has got nothing to do */
    grid_gen_work_init(&ggw, mgrid, j);
    
    for(atom=off, done=0; atom<mgrid->atom_cnt; 
        atom+=mgrid->atom_coords[atom].mult, done++) {
        GridGenAtomGrid* grid = mgrid->atom_grids[atom];
        int mult = mgrid->atom_coords[atom].mult;
        if(done % mgrid->nt != 0) /* do every mt symmetry-independent atom */
            continue; 
        idx = 0;
        for (ptno=0; ptno<grid->pnt; ptno++) {
            real fact = 4*M_PI*grid->wght[ptno]*mult;
            real rad = grid->rad[ptno];
            int ind = grid->leb_ang[ptno]; 
            leb_gen[ind].func(ggw.x+idx, ggw.y+idx, ggw.z+idx, ggw.wg+idx,
                              &tpt);

            for(j=idx; j<idx+leb_gen[ind].point_cnt; j++) {
                ggw.x[j] = ggw.x[j]*rad + mgrid->atom_coords[atom].x;
                ggw.y[j] = ggw.y[j]*rad + mgrid->atom_coords[atom].y;
                ggw.z[j] = ggw.z[j]*rad + mgrid->atom_coords[atom].z;
                ggw.wg[j] *= fact;
            }
            if(selected_partitioning->preprocess)
                selected_partitioning->preprocess
                    (mgrid, atom, leb_gen[ind].point_cnt, &ggw, idx, 0);
            idx += leb_gen[ind].point_cnt;
        }
        /* degeneracy multiplication here */
        cnt = compress_grid(idx, ggw.x, ggw.y, ggw.z, ggw.wg);

        fort_print("Atom: %4d*%d points=%5d compressed from %5d (%3d radial)", 
                   atom+1, mult, cnt, idx, grid->pnt);
        if(cnt>0) {
            int i;
            mol_grid_lock(mgrid);
            fwrite(&cnt,  sizeof(int), 1,  mgrid->fl);
            fwrite(&atom, sizeof(int), 1,  mgrid->fl);
            for(i=0; i<cnt; i++) {
                fwrite(ggw.x+i, sizeof(real), 1, mgrid->fl);
                fwrite(ggw.y+i, sizeof(real), 1, mgrid->fl);
                fwrite(ggw.z+i, sizeof(real), 1, mgrid->fl);
            }
            fwrite(ggw.wg,sizeof(real), cnt, mgrid->fl);
            mgrid->total_points += cnt;
            mol_grid_unlock(mgrid);
        }
    }
    grid_gen_work_release(&ggw);
    return NULL;
}

int
grid_gen_save_temp(const char* filename, GridGenMolGrid* mgrid)
{
    if( (mgrid->fl = fopen(filename,"wb"))==NULL) {
        fprintf(stderr,"ERROR: Cannot open grid file '%s' for writing.\n",
                filename);
        return 0;
    }
    mgrid->total_points = 0;
#if defined(USE_PTHREADS)
    mgrid->nt = mol_grid_get_num_threads();
    { int i;
    pthread_t *ptid = dal_malloc(mgrid->nt*sizeof(pthread_t));
    for(i=0; i<mgrid->nt; i++) {
        mol_grid_lock(mgrid);
        mgrid->off = i;
        pthread_create(&ptid[i], NULL, 
                       (void *(*)(void *))grid_gen_save_temp_stream, mgrid);
    }
    for(i=0; i<mgrid->nt; i++) 
        pthread_join(ptid[i], NULL);
    }
#else
    mgrid->nt = 1; mgrid->off = 0;
    grid_gen_save_temp_stream(mgrid);
#endif
    fclose(mgrid->fl);

    return mgrid->total_points;
}


/* grid_gen_generate:
   returns number of grid points.
*/
int
grid_gen_generate(const char* filename, int atom_cnt, 
                  const GridGenAtom* atom_arr, real threshold,
                  GridGeneratingFunc generating_function, void* arg,
                  int minang, int maxang, real* work, int *lwork)
{
    int res;
    struct tms starttm, endtm; clock_t utm;
    GridGenMolGrid* mgrid =  mol_grid_new(atom_cnt, atom_arr);
    fort_print("Radial Quadrature : %s", radial_quad->name);
    fort_print("Space partitioning: %s", selected_partitioning->name);
    fort_print("Radial integration threshold: %g", threshold);

    times(&starttm);
    set_radial_grid(mgrid, threshold);

    set_ang_fixed(mgrid, threshold, get_leb_from_type(minang), 
                  get_leb_from_type(maxang));

    res = grid_gen_save_temp(filename, mgrid);

    res = create_cubes(mgrid, filename, res, CELL_SIZE,
                       work, *lwork*sizeof(real));
    mol_grid_free(mgrid);

    times(&endtm);
    utm = endtm.tms_utime-starttm.tms_utime;

    fort_print("Number of grid points: %8d Grid generation time: %9.1f s\n", 
               res, utm/(double)sysconf(_SC_CLK_TCK));
    return res;
}

/* =================================================================== */
/* fixed grid generator                                                */
/* =================================================================== */


/* ------------------------------------------------------------------- */
/* the dft grid input routines.                                        */
/* choose different types of grid:                                     */
/* The syntax is:                                                      */
/* ((ADAPT|FIXED) | (GC2|MALMQUIST) | (SSF|BECKE))*                    */
/* Sets                                                                */
/* radial_quad                                                         */
/* selected_partitioning                                               */
/* ------------------------------------------------------------------- */
static int
get_word(const char *line, int st, int max_len)
{
    int res;
    for(res=st; res<max_len && isalnum(line[res]); res++)
        ;
    if(res>=max_len) res = 0;
    return res;
}

static int
skip_spaces(const char *line, int st, int max_len)
{
    int res;
    for(res=st; res<max_len && isspace(line[res]); res++)
        ;
    return res;
}

void
dftgridinput_(const char *line, int line_len)
{
    static const char* keywords[] = {
        "GC2","LMG", "SSF","BECKE","BECKE2", "BLOCK"
    };
    int st, en, i;
    for(st=0; (en=get_word(line, st, line_len)) != 0; 
        st=skip_spaces(line, en, line_len)) {
        for(i=0; i<ELEMENTS(keywords); i++) {
            if(strncasecmp(keywords[i], line+st, en-st)==0)
                break;
        }
        switch(i) {
        case 0: radial_quad = &quad_gc2; break;
        case 1: radial_quad = &quad_lmg; break;
        case 2: selected_partitioning=&partitioning_ssf;break;
        case 3: 
            selected_partitioning = &partitioning_becke1; break;
        case 4: 
            selected_partitioning = &partitioning_becke2; break;
        case 5: 
            selected_partitioning = &partitioning_block; break;
        default: fort_print("GRIDGEN: Unknown .GRID TYPE option ignored.\n%s",
                            line);
            /* FIXME: should I quit here? */
        }
    }
}

/* =================================================================== */
/* linear grid generator                                               */
/* =================================================================== */
/*the bucketing scheme. first shot: nlog n. Load up all grid points,
  assign them [nx,ny,nz] tuples.  sort the tuples to collect the ones
  sharing same cube. Within each cube select one that is closest to
  the center. Save cubes.
*/
#if defined(SYS_DEC)
#if defined(CHAR_BIT)
#undef CHAR_BIT
#endif
/* other architectures define it properly on their own */
#define CHAR_BIT 8
#endif
#define KEY_BITS (sizeof(GridPointKey)*CHAR_BIT)

typedef unsigned GridPointKey;
struct point_key_t {
    GridPointKey key; /* unique identificator of the box */
    unsigned     index;
};


/* assume currently constant amounts of bits per coordinate */
static int
comp_point_key(const void* a, const void* b)
{
    return ((struct point_key_t*)a)->key-((struct point_key_t*)b)->key;
}

void gtexts_(real* r2);
void getblocks_(real *X,real *Y, real *Z, real *CELLSZ, real RSHEL2[],
                int *NBLCNT, int (*IBLCKS)[2]);
static void
load_grid(const char *fname, int point_cnt, real *work, int worksz,
          real (**coor)[3], real **w, int *x_allocated,
          int *atom_idx)
{
    int idx, cnt;
    real *chunk;
    FILE *f;

    if(worksz<4*point_cnt*sizeof(real)) {
        fprintf(stderr, "wrkmem too small (%d, needed %d) trying malloc.\n",
                worksz, 4*point_cnt*sizeof(real));
        *x_allocated = 1;
        chunk = malloc(4*point_cnt*sizeof(real));
        if(!chunk) {
            fprintf(stderr, "loading grid into mem failed, too.\n"
                    "worksz=%d needed size=%d\n", worksz,
                    4*point_cnt*sizeof(real));
            dalton_quit("no mem in load_grid: 2");
        }
    } else chunk = work;
    *coor = (real(*)[3]) chunk;
    *w    = chunk + 3*point_cnt;
     
    f = fopen(fname,"rb");
    if(!f) {
        fprintf(stderr, "internal error, cannot open the grid file %s.\n",
                fname);
        exit(1);
    }
    idx = 0;
    while(fread(&cnt, sizeof(cnt), 1, f)==1) {
        int aidx, i;
        assert(cnt+idx<=point_cnt);
        if(fread(&aidx, sizeof(int), 1, f) != 1) dalton_quit("ERROR GRID1");
        if(atom_idx)
            for(i=0; i<cnt; i++) atom_idx[idx+i] = aidx;
        if(fread(*coor+idx, sizeof(real), 3*cnt, f) != 3*cnt)
            dalton_quit("ERROR GRID1, cnt=%d", cnt);          
        if(fread(*w   +idx, sizeof(real), cnt, f) != cnt)
            dalton_quit("ERROR GRID2, cnt=%d", cnt);
        idx += cnt;
    }
    fclose(f);
}

static void
create_index(real cell_size, real (*coor)[3], int point_cnt,
             struct point_key_t * index,
             struct point *lo, struct point *hi)
{
    real fac;
    int i, bpc = KEY_BITS/3; /* bits per coordinate */
    lo->x = lo->y = lo->z = +1e20;
    hi->x = hi->y = hi->z = -1e20;

    for(i=0; i<point_cnt; i++) {
        if(coor[i][0] < lo->x)      lo->x = coor[i][0];
        else if(coor[i][0] > hi->x) hi->x = coor[i][0];
        if(coor[i][1] < lo->y)      lo->y = coor[i][1];
        else if(coor[i][1] > hi->y) hi->y = coor[i][1];
        if(coor[i][2] < lo->z)      lo->z = coor[i][2];
        else if(coor[i][2] > hi->z) hi->z = coor[i][2];
    }
    if( (1<<bpc) < (hi->x - lo->x)/cell_size) fputs("error: x\n", stderr);
    if( (1<<bpc) < (hi->y - lo->y)/cell_size) fputs("error: y\n", stderr);
    if( (1<<bpc) < (hi->z - lo->z)/cell_size) fputs("error: z\n", stderr);

    fac = 1/cell_size;
    for(i=0; i<point_cnt; i++) {
        int ix = (coor[i][0] - lo->x)*fac;
        int iy = (coor[i][1] - lo->y)*fac;
        int iz = (coor[i][2] - lo->z)*fac;
        if(ix<0) ix = 0; /* correct numerical error, if any */
        if(iy<0) iy = 0; /* correct numerical error, if any */
        if(iz<0) iz = 0; /* correct numerical error, if any */
        index[i].key = (ix<<(bpc*2)) | (iy<<bpc) | iz;
        index[i].index = i;
    }

    qsort(index, point_cnt, sizeof(struct point_key_t), comp_point_key);
}

/* ===================================================================
 * The parallelization section. Depending on the calculation mode, the
 * code uses different initialization, save_batch and finalize
 * actions. */
static void
save_final_batch_local(FILE *f, int cnt, int nblocks, int shlblocks[][2],
                       real *coor, real *w)
{
    if(fwrite(&cnt, sizeof(cnt), 1, f)!=1) {
        fprintf(stderr, "GRIDGEN: 'too short write' error.\n");
        exit(1);
    }
    if(fwrite(&nblocks,  sizeof(nblocks), 1, f) != 1) abort();
    if(fwrite(shlblocks, sizeof(int), nblocks*2, f) != nblocks*2)
        dalton_quit("write error in %s(), point 1", __FUNCTION__);
    if(fwrite(coor, sizeof(real), 3*cnt, f) != 3*cnt)
        dalton_quit("write error in %s(), point 2", __FUNCTION__);
    if(fwrite(w, sizeof(real), cnt, f) != cnt)
        dalton_quit("write error in %s(), point 3", __FUNCTION__);
}

#ifdef VAR_MPI
#include <mpi.h>
static int mynum, nodes, last;
static void
grid_par_init(void) {
    /* Executed by master and all slaves */
    MPI_Comm_rank(MPI_COMM_WORLD, &mynum);
    MPI_Comm_size(MPI_COMM_WORLD, &nodes);
    last = 0;
}
static char*
grid_get_fname(const char *base, int filenum)
{
    if(filenum == 0)
        return strdup(base);
    else {
        char *res = malloc(strlen(base) + 15);
        sprintf(res, "%s.%05d", base, filenum);
        return res;
    }
}

#define M(s) {if((s) != MPI_SUCCESS) printf("MPI comm failed at %d\n", __LINE__); }
static void
save_final_batch(FILE *f, int cnt, int nbl, int shlbl[][2],
                 real *coor, real *w)
{
    /* PARBLLEN must be low multiplicity of dftcom:MAXBLLEN
     * for performance reasons. */
    const int PARBLLEN = 1000;
    /* Executed by master */
    int i;
    for(i=0; i<cnt; i+= PARBLLEN) {
        int bcnt = i+PARBLLEN<cnt ? PARBLLEN : cnt - i;
        last = (last+1) % nodes;
        if( last == 0)
            save_final_batch_local(f, bcnt, nbl, shlbl, coor + i*3, w+i);
        else {
            int arr[2]; arr[0] = bcnt; arr[1] = nbl;
            M(MPI_Send(arr,     2,      MPI_INT,   last, 1, MPI_COMM_WORLD));
            M(MPI_Send(shlbl,   2*nbl,  MPI_INT,   last, 2, MPI_COMM_WORLD));
            M(MPI_Send(coor+i*3,3*bcnt, MPI_DOUBLE,last, 3, MPI_COMM_WORLD));
            M(MPI_Send(w+i,     bcnt,   MPI_DOUBLE,last, 4, MPI_COMM_WORLD));
        }
    }
}

static void
grid_par_shutdown(void)
{
    unsigned i;
    int arr[2]; arr[0] = 0; arr[1] = 0;
    /* Executed by master only */
    for(i=1; i<nodes; i++)
         M(MPI_Send(arr, 2, MPI_INT, i, 1, MPI_COMM_WORLD));
}
static void
grid_par_slave(const char *fname)
{
    int *shlblocks;
    MPI_Status s;
    char *     nm = grid_get_fname(fname, mynum);
    FILE *f;
    int ishlcnt = FSYM(ishell_cnt)();
    real *dt = NULL;
    int   dt_sz = 0;

    if((f=fopen(nm, "wb")) == NULL) dalton_quit("Slave could not save.");
    shlblocks = malloc(ishlcnt*sizeof(int));
    do {
        int arr[2];
        M(MPI_Recv(arr, 2, MPI_INT, 0, 1, MPI_COMM_WORLD, &s));
        if(arr[0] == 0)break; /* End Of Job */

        M(MPI_Recv(shlblocks, 2*arr[1], MPI_INT, 0, 2, MPI_COMM_WORLD, &s));
        if(arr[0]>dt_sz) {
            dt = realloc(dt, 4*arr[0]*sizeof(real));
            dt_sz = arr[0];
        }
        M(MPI_Recv(dt,        3*arr[0], MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &s));
        M(MPI_Recv(dt+3*arr[0], arr[0], MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, &s));
        if(fwrite(arr,       sizeof(int), 2, f) != 2) abort();
        if(fwrite(shlblocks, sizeof(int), arr[1]*2, f) != arr[1]*2)
            dalton_quit("write error in %s(), point 2", __FUNCTION__);
        /* write coords */
        fwrite(dt, sizeof(real), 3*arr[0], f);
        /* write weights */
        fwrite(dt+3*arr[0], sizeof(real), arr[0], f);
    } while(1);
    fclose(f);
    if(dt)
        free(dt);
    free(nm);
    free(shlblocks);
}
  
#else
#define save_final_batch(f,c,b,s,r,w) \
        save_final_batch_local((f),(c),(b),(s),(r),(w)) 
#define grid_get_fname(base,num) strdup(base)

#endif

/* save_final saves the grid, possibly distributing it over to many
 * nodes. NOTE: atom_idx is used only for partitioning schemes that do
 * postprocessing and should not be accessed without prior checking.
 **/
static int
save_final(GridGenMolGrid *mg, const char *fname, int point_cnt,
           real (*coor)[3], real *w, const int *atom_idx,
           struct point_key_t *keys, struct point *lo, real cell_size)
{
    FILE *f;
    int idx, cnt, newcnt;
    real *rshel2;
    int nblocks, (*shlblocks)[2] = malloc(2*ishell_cnt_()*sizeof(int));
    int bpc = KEY_BITS/3; /* bits per coordinate */
    GridPointKey mask = ~((-1)<<bpc);
    int points_saved = 0;
    real *dt = NULL;
    int dt_sz = 0, *atom_nums = NULL;

    if((f = fopen(fname,"wb")) == NULL) {
        fprintf(stderr,"internal error, cannot save sorted grid file.\n");
        exit(1);
    }
    rshel2 = dal_malloc(ishell_cnt_()*sizeof(real));
    gtexts_(rshel2);
    for(idx=0; idx<point_cnt; idx += cnt) {
        int closest = 0, i;
        unsigned tmp;
        struct point c;
        real mindist = 4*cell_size, maxdist = 0;
        GridPointKey key = keys[idx].key;
        real sx = 0, sy = 0, sz = 0;
        c.x = lo->x + (((key >> (bpc*2)) & mask)+0.5)*cell_size;
        c.y = lo->y + (((key >> bpc)     & mask)+0.5)*cell_size;
        c.z = lo->z + (((key)            & mask)+0.5)*cell_size;
        for(cnt=0; idx+cnt<point_cnt &&
                keys[idx+cnt].key == key; cnt++) {
            real dx = coor[keys[idx+cnt].index][0]-c.x;
            real dy = coor[keys[idx+cnt].index][1]-c.y;
            real dz = coor[keys[idx+cnt].index][2]-c.z;
            real dist2 = dx*dx + dy*dy + dz*dz;
            if(dist2<mindist) { mindist=dist2; closest = cnt; }
            if(dist2>maxdist) { maxdist=dist2; }
            sx += dx; sy += dy; sz += dz;
        }

        /* merge data in one block */
        if(dt_sz<cnt) {
            dt_sz = cnt;
            dt = realloc(dt, 4*dt_sz*sizeof(real));
            atom_nums = realloc(atom_nums, dt_sz*sizeof(int));
        }
        for(i=0; i<cnt; i++) {
            int j = keys[idx+i].index;
            dt[i*3+0]    = coor[j][0];
            dt[i*3+1]    = coor[j][1];
            dt[i*3+2]    = coor[j][2];
            dt[cnt*3+i]  = w[j];
            if(atom_idx) atom_nums[i] = atom_idx[j];
        }
        getblocks_(&c.x, &c.y, &c.z, &cell_size, rshel2, &nblocks, shlblocks);
        if(nblocks==0) continue;        

        if(0) {
        printf("box [%5.1f,%5.1f,%5.1f] nt=%5d "
               "%f <r<%f [%5.1f,%5.1f,%5.1f] %d\n",
               c.x, c.y, c.z, cnt, sqrt(mindist), sqrt(maxdist),
               sx/cnt, sy/cnt, sz/cnt, nblocks);
        for(i=0; i<nblocks; i++)
            printf("(%d,%d)",shlblocks[i][0], shlblocks[i][1]);
        puts("");
        }

        if(selected_partitioning->postprocess)
            newcnt = selected_partitioning->postprocess(mg, &c,cnt, atom_nums,
                                                        (real(*)[3])dt,
                                                        dt+3*cnt);
        else
            newcnt = cnt;

        save_final_batch(f, newcnt, nblocks, shlblocks, 
                         dt, dt+3*cnt);
        points_saved += newcnt;
    }
    fclose(f);
    free(rshel2);
    free(shlblocks);
    if(dt) free(dt);
    if(atom_nums) free(atom_nums);
    if(point_cnt != points_saved)
        fort_print("Postprocessing compression: from %d to %d\n",
                    point_cnt, points_saved);
    return points_saved;
}

static int
create_cubes(GridGenMolGrid *mg, const char* fname,
             int point_cnt, real cell_size,
             void* work, int worksz)
{
    real (*coor)[3], *w;
    struct point lo, hi;
    struct point_key_t * keys =
      dal_malloc(point_cnt*sizeof(struct point_key_t));
    int coorw_allocated = 0;
    int *atom_idx;
    int new_point_cnt;

    atom_idx = selected_partitioning->postprocess ?
        dal_malloc(point_cnt*sizeof(int)) : NULL;

    load_grid(fname, point_cnt, work, worksz,
              &coor, &w, &coorw_allocated, atom_idx);
    create_index(cell_size, coor, point_cnt, keys, &lo, &hi);
    new_point_cnt = save_final(mg, fname, point_cnt, coor, w, atom_idx,
                               keys, &lo, cell_size);
    free(keys);
    if(coorw_allocated) free(coor);
    if(atom_idx) free(atom_idx);
    return new_point_cnt;
}

/* =================================================================== */
/* Grid I/O routines.
 *
 * Used by integrator to fetch data. The routines are located here
 * since the grid file format is an internal implementation issue of
 * the grid generator: it may contain not only grid positions and
 * weights but also other metadata.
 * 
 * Two routines are provided. One is designed for blocked approach,
 * the other one is a plain point-by-point grid reader.
 */
/* private definition */
struct DftGridReader_ {
    FILE *f;
};
void get_grid_paras_(int *grdone, real *radint, int *angmin, int *angint);
void set_grid_done_(void);

DftGridReader*
grid_open(int nbast, real *work, int *lwork)
{
    DftGridReader *res = dal_malloc(sizeof(DftGridReader));
    int grdone, angmin, angint;
    real radint;
    char *fname;

    get_grid_paras_(&grdone, &radint, &angmin, &angint);
    if(!grdone) {
        int atom_cnt, pnt_cnt;
        int lwrk = *lwork - nbast;
        GridGenAtom* atoms = grid_gen_atom_new(&atom_cnt);
        struct RhoEvalData dt; /* = { grid, work, lwork, dmat, dmgao}; */
        dt.grid =  NULL; dt.work = work+nbast; dt.lwork = &lwrk;
        dt.dmat =  NULL; dt.dmagao = work;

#ifdef VAR_MPI
        grid_par_init();
        if(mynum == 0) {
            pnt_cnt = grid_gen_generate("DALTON.QUAD", atom_cnt, atoms,
                                        radint, NULL, &dt,
                                        angmin, angint, work, lwork);
            grid_par_shutdown();
        } else 
            grid_par_slave("DALTON.QUAD");
        /* Stop on barrier here so that we know all nodes managed to save
         * their files. */
        MPI_Barrier(MPI_COMM_WORLD);
#else
        pnt_cnt = grid_gen_generate("DALTON.QUAD", atom_cnt, atoms,
                                    radint, NULL, &dt,
                                    angmin, angint, work, lwork);
#endif
        free(atoms);
	set_grid_done_();
    }
    fname = grid_get_fname("DALTON.QUAD", mynum);
    res->f=fopen(fname, "rb");
    free(fname);
    if(res == NULL) {
	perror("DFT quadrature grid file DALTON.QUAD not found.");
	free(res);
	abort();
    }
    return res;
}

/** grid_getchunk_blocked() reads grid data also with screening
    information if only nblocks and shlblocks are provided.
 */
int
grid_getchunk_blocked(DftGridReader* rawgrid, int maxlen,
                      int *nblocks, int *shlblocks, 
                      real (*coor)[3], real *weight)
{
    int sz = 0, i, rc, bl_cnt;
    FILE *f = rawgrid->f;

    if(fread(&sz, sizeof(int), 1, f) <1)
        return -1; /* end of file */
    if(sz>maxlen) {
        fprintf(stderr,
                "grid_getchunk: too long vector length in file: %d > %d\n"
                "Calculation will stop.\n", sz, maxlen);
        dalton_quit("grid_getchunk: too long vector length in file: %d > %d\n"
                "Calculation will stop.\n", sz, maxlen);
        return -1; /* stop this! */
    }

    if(fread(&bl_cnt, sizeof(unsigned), 1, f) <1) {
        puts("OCNT reading error."); return -1;
    }
    if(nblocks) *nblocks = bl_cnt;

    if(shlblocks) {
        rc = fread(shlblocks, sizeof(int), bl_cnt*2, f);
    } else {
        int buf, cnt;
        for(cnt=0; cnt<bl_cnt*2; cnt+=rc) {
            rc=fread(&buf, sizeof(int), 1, f);
            if( rc < 1)
                break;
        }
    }
    if(rc<1) {
        fprintf(stderr,
                "IBLOCKS reading error: failed to read %d blocks.\n",
                *nblocks);
        return -1; 
    }

    if(fread(coor,   sizeof(real), 3*sz, f) < sz) { puts("XYZ");return -1;}
    if(fread(weight, sizeof(real), sz, f) < sz) { puts("W");return -1;}
    return sz;
}


void
grid_close(DftGridReader *rawgrid)
{
    fclose(rawgrid->f);
    free(rawgrid);
}

/* ------------------------------------------------------------------- */
/*                     FORTRAN INTERFACE                               */
/* ------------------------------------------------------------------- */
static DftGridReader *grid = NULL;
void
opnqua_(int *nbast, real *work, int *lwork)
{
    grid = grid_open(*nbast, work, lwork);
}

void
clsqua_(void)
{
    grid_close(grid);
    grid = NULL;
}

void
reaqua_(int *nshell, int *shell_blocks, int *buf_sz, 
        real (*coor)[3], real *weight, int *nlen)
{
    *nlen = grid_getchunk_blocked(grid, *buf_sz, nshell, shell_blocks,
                                  coor, weight);
}

/* ------------------------------------------------------------------- */
/* self test routines of the integrator/grid generator                 */
/* ------------------------------------------------------------------- */

#if defined(GRID_GEN_TEST) && GRID_GEN_TEST == 1
static const TestMolecule* test_molecule = NULL;
typedef struct TestMolecule_ TestMolecule;
struct TestMolecule_ {
    const GridGenAtom* atoms;
    const real* overlaps;
    int atom_cnt;
};


const static GridGenAtom hydrogen1_atoms[] = {
    { 0.0, 0.0, 0.0, 1}  /* single H atom at (0,0,1) */
};
const static real hydrogen1_overlaps[] = { 1 };

const static GridGenAtom hydrogen2_atoms[] = {
    { 0.0, 0.0, 0.0, 1}, /* single H atom at (0,0,1) */
    { 0.0, 0.0, 1.0, 1}  /* single H atom at (0,0,2) */
};

const static real hydrogen2_overlaps[] = {
    0.53726234893, 0.53726234893
};
const static TestMolecule test_molecules[] = { 
    { hydrogen1_atoms, hydrogen1_overlaps, ELEMENTS(hydrogen1_atoms) },
    { hydrogen2_atoms, hydrogen2_overlaps, ELEMENTS(hydrogen2_atoms) }
};

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* Test functions needed for testing the integrator/grid generator     */
/* ------------------------------------------------------------------- */

/* wave function evaluator;
 * does not include overlaps, this should be fixed. */
static void
get_rho_grad(real x, real y, real z, real* rho, real* grad)
{
    const static real COEF=1.99990037505132565043;
    const static real  ALPHA=0.62341234;
    real val = 0, gradx=0, grady=0, gradz=0;
    int i;
    val = 0;
    for(i=0; i<test_molecule->atom_cnt; i++) {
        real dx = x-test_molecule->atoms[i].x, 
            dy  = y-test_molecule->atoms[i].y,
            dz  = z-test_molecule->atoms[i].z;
        real r2 = dx*dx+dy*dy+dz*dz;
        real ao = test_molecule->overlaps[i]*exp(-ALPHA*r2);
        val  += ao;
        if(grad) {
            real fact = 4*test_molecule->overlaps[i]*ALPHA*exp(-ALPHA*r2)/COEF;
            gradx -= dx*fact;       
            grady -= dy*fact;       
            gradz -= dz*fact;       
        }
    }
    val /= COEF;
    *rho = val*2*val;

    if(grad) {
        gradx *= 2*val;
        grady *= 2*val;
        gradz *= 2*val;
        /* printf("grad: (%g,%g,%g)\n", gradx, grady, gradz); */
        *grad = sqrt(gradx*gradx + grady*grady + gradz*gradz);
    }
}

static real
norm3(real x, real y, real z)
{
    const static real  ALPHA=0.62341234;
    real r2 = x*x + y*y+ z*z;
    real norm= exp(-2*ALPHA*r2);
    return norm;
}

static real
rho(real x, real y, real z)
{
    real rho;
    get_rho_grad(x,y,z, &rho, NULL);
    return rho;
}


/* dirac3:
   single atom: -.76151426465514650000
*/
static real
dirac3(real x, real y, real z)
{
    const real PREF= -3.0/4.0*pow(6/M_PI, 1.0/3.0);
    real rhoa;
    get_rho_grad(x,y,z, &rhoa, NULL);
    rhoa *= 0.5;
    return PREF*2*(pow(rhoa,4.0/3.0));
}

static real
becke3(real x, real y, real z)
{
    real rhoa, grada;
    real xa,asha, ea;
    const static real BETA=0.0042;

    get_rho_grad(x,y,z, &rhoa, &grada);
    rhoa *= 0.5; grada *=0.5;
    if(rhoa<1e-10) return 0;

    xa = grada*pow(rhoa,-1.0/3.0)/rhoa;
    asha = asinh(xa);
    ea  = -2*grada*xa*BETA/(1+6*xa*BETA*asha);
    return ea;
}

static real
lyp3(real x, real y, real z)
{
    const real A = 0.04918, B = 0.132, C = 0.2533, D = 0.349;
    const real CF = 0.3*pow(3*M_PI*M_PI,2.0/3.0);

    real rho, rhom13, denom, omega, delta, ret, ngrad2, ngrada2, ngradb2;
    real rho2, t1, t2, t3, t4, t5, t6;
    real rhoa = 0.0, rhob = 0.0, grada = 0.0, gradb = 0.0;

    get_rho_grad(x,y,z, &rhoa, &grada);
    if(rhoa<1e-10) return 0;
    rhob=rhoa;
    gradb=grada;

    rho = rhoa+rhob;
    rho2 = rho*rho;
    ngrad2  = (grada+gradb)*(grada+gradb);
    ngrada2 = grada*grada;
    ngradb2 = gradb*gradb;
    rhom13 = pow(rho,-1.0/3.0);
    denom = 1+D*rhom13;
    omega = exp(-C*rhom13)/denom*pow(rho,-11.0/3.0);
    delta = rhom13*(C + D/denom);
    t1 =  pow(2.0,11.0/3.0)*CF*(pow(rhoa,8.0/3.0) +pow(rhob,8.0/3.0));
    t2 =  (47.0 - 7.0*delta)*ngrad2/18.0;
    t3 = -(2.5 -delta/18.0)*(ngrada2+ngradb2);
    t4 =  (11.0-delta)/9.0*(rhoa*ngrada2 + rhob*ngradb2)/rho;
    t5 = -2.0/3.0*rho2*ngrad2;
    t6 = ((2.0/3.0*rho2-rhoa*rhoa)*ngradb2 +
          (2.0/3.0*rho2-rhob*rhob)*ngrada2);
    ret = -A*(4*rhoa*rhob/(denom*rho)
              +B*omega*(rhoa*rhob*(t1+t2+t3+t4)+t5+t6));
    return ret;
}

/* ------------------------------------------------------------------- */
/* main test driver                                                    */
/* ------------------------------------------------------------------- */
int
main(int argc, char* argv[])
{
    int i, arg, no;
    real THRESHOLD = 1e-3;
    const static struct ff {
        GGFunc func;
        char * name;
    } functions_to_try[] = {
        /*{ norm3, "Norm" }, */
        { dirac3,"Dirac " },
        { becke3,"Becke " },
        { rho,   "charge" }
        /* { lyp3,  "LYP" } */
    };
    include_selected_partitioning = include_partitioning_becke1; 

    if(argc<=1) {
        fprintf(stderr, "No arguments specified.\n"
                "\t-v    : increase verbosity level.\n"
                "\t-t THR: specify new threshold level.\n"
                "\t-b    : use Becke partitioning (default).\n"
                "\t-s    : use SSF partitioning.\n"
                "\tnumber: molecule to run.\n");
        return 1;
    }
    for(arg=1; arg<argc; arg++) {
        if(strcmp(argv[arg], "-v")==0)
            { verbose++; puts("Verbose++"); }
        else if(strcmp(argv[arg], "-t")==0) {
            THRESHOLD = atof(argv[++arg]);
            printf("Threshold set to: %g\n", THRESHOLD);
        } else if(strcmp(argv[arg], "-b")==0) {
            include_selected_partitioning = include_partitioning_becke1; 
            puts("Becke partitioning"); 
        } else if(strcmp(argv[arg], "-s")==0) {
            include_selected_partitioning = include_partitioning_ssf; 
            puts("SSF partitioning"); 
        } else if( (no=atoi(argv[arg]))>0 && 
                  no<= ELEMENTS(test_molecules)) {
            GridGenMolGrid* mgrid;
            printf("Testing molecule: %d\n", no);
            no--;
            test_molecule = &test_molecules[no];
            mgrid = 
                mol_grid_new(test_molecule->atom_cnt, test_molecule->atoms, 
                             THRESHOLD);
            
            for(i=0; i< ELEMENTS(functions_to_try); i++) {
                if(verbose)
                    printf("%5s: Number of radial grid points: %d\n", 
                           functions_to_try[i].name,
                           mgrid->atom_grids[0]->pnt);
                set_radial_grid(mgrid);
                set_ang(mgrid, THRESHOLD, functions_to_try[i].func);
                /* integration */
                printf("%5s: Integrated function         : %20.15f\n",
                       functions_to_try[i].name,
                       integrate(mgrid, functions_to_try[i].func));
            }
            mol_grid_free(mgrid);
        }
    }
    return 0;
}
#endif /* defined(GRID_GEN_TEST) && GRID_GEN_TEST == 1 */
