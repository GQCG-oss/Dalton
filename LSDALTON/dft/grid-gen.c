/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/*    THIS FILE CONTAINS ROUTINES SIMILAR TO ROUTINES IN grid-gen.c    */
/*    BUT THESE ROUTINES CIRCUMVENT CALLS TO ORTRAN77 ROUTINES         */
/*    WHICH USES COMMONBLOCKS, DUE TO THE STRUCATURE OF THE CODE       */
/*    I HAVE BEEN UNABLE TO INCOORPERATE THIS INTO THE STANDARD        */
/*    ROUTINES IN grid-gen.c                                           */ 
/* ------------------------------------------------------------------- */
/*    DECLARATIONS                                                     */

/*
   The code consists of four main parts:
   1. angular grid generation.
   2. radial grid generation.
   3. partitioning.
   4. box-based screening code.

   Grid is generated as follows:
   A. list of atoms is created.

   B. radial grids are created for each atom (FIXME: generate them
   only for grid type just to have nicer code, no performance impact).

   C. radial grid is multiplied with the angular one. Angular can be
   pruned for small R in order to avoid too many grid points with
   small weights.

   C1. if periodic boundary conditions are used: grid points sticking
   out of the master Cell are wrapped back by shifting appropriately
   the positions of the atoms they are associated with.
   
   D. partitioning preprocessing is performed (if any) and grid points
   are dumped to a file (phase A).
   
   E. grid is read back and it is sorted to boxes. All the grid points
   in the box share same screening parameters (active orbitals).

   F. grid point batches are distributed to different processing units
   (thread, mpi processes - if more than one) and postprocessing phase
   of the partitioning is applied.

   G. Final grid is saved.

   Code layout follows this scheme as well: general atom/grid
   manipulation routines, radial grid routines, angular grid routines,
   partitioning routines, phase A routines, boxing routines, phase B
   routines.
*/

/* define VAR_MPI if you want to use parallel grid generation 
 * for later parallel MPI calculation.
 *
 * define USE_PTHREADS if you want to use multithreaded grid
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

#if defined(VAR_OMP)
#include <omp.h>
#endif

#include "grid-gen.h"

/* common screening for cubes/cells */
/* static const real CELL_SIZE = 3.2;*/
static const real CELL_SIZE = 2.0;
/* the weight threshold:: a grid point will be ignored if has a weight
 * lower than the WEIGHT_THRESHOLD. This should probably depend on the
 * calculation's accuracy. */
static const real WEIGHT_THRESHOLD = 1e-15;

/* cutoff value for cutoff profile in SSF weight scheme*/
static const real SSF_CUTOFF=0.64;

enum { GRID_TYPE_STANDARD, GRID_TYPE_CARTESIAN };
    
static integer gridType = GRID_TYPE_STANDARD;

void*
dal_malloc_(size_t sz, const char *place, unsigned line)
{
    void* res = malloc(sz);
    if(!res) {
        fprintf(stderr, "dal_malloc(sz=%u bytes) at %s (line %u) failed.\n",
                (unsigned)sz, place, line);
        exit(1);
    }
    return res;
}

/* ------------------------------------------------------------------- */
/* Internal data structures used by the grid generator.                */
/* ------------------------------------------------------------------- */
struct GridGenAtomGrid_ {
    integer uniq_no;  /* number of the atom, before the symmetry multiplication */
    integer Z;        /* frequently used. */
    integer* leb_ang; /* leb_gen index associated with each radial point */
    real* rad;        /* radial points array. */
    real* wght;       /* radial points weights */
    integer pnt;      /* number of radial points for this atom*/
};  

struct GridPbcData {
    void (*cart2lat)(const real *cart, real *lat);
    void (*lat2cart)(const real *lat, real *cart);
    integer n[3];
    integer maxl; /**< number of layers for one electron operators */
};

struct radial_scheme_t;
struct partitioning_scheme_t;
struct DftGrid_ {
    struct radial_scheme_t       *radial_quad;
    struct partitioning_scheme_t *partitioning;
    unsigned Z_dependent_maxang:1;  /**< whether to use Z-dependent maxang. */
    integer prune;    /* pruning on or off*/
    integer hardness; /*hardness of the partionining function*/
    struct GridPbcData pbc;
    real radint;
    integer angmin, angmax;
    integer regenerate;
    unsigned pbc_enabled:1;
};

struct GridGenMolGrid_ {
    integer atom_cnt;                  /* number of atoms grid is generated for*/
    const GridGenAtom* atom_coords;/* array with atom data */
    GridGenAtomGrid**   atom_grids;/* array of atom grids */
    real* rij;  /* triangular array of inverse distances between atoms */
    real* distij; /* triangular array of distances between atoms */
    real* nn;  /*distance to nearest neighbor for every atom*/
    real* aij;  /* triangular array of of.... */
    FILE* fl;   /* the stream grid is being saved to */
    integer total_points; /* total number of generated points */
    integer off; /* thread number. 0 in serial. */
    integer nt;  /* total number of threads. 1 in serial. */
    DftGrid *conf; /* configuration */
    struct GridPbcData *pbc;
    unsigned verbose:1; /* whether grid generation should be verbose or not */

    /*some constants for the SSF scheme*/
    real ssf_fac0;
    real ssf_fac1;
    real ssf_fac2;
    real ssf_fac3;
    real ssf_fac4;

    /*ANDREAS: the following should go to another structure, preliminary only!*/
    integer finefac;
    real finefacinv;
    real* tabcutoff0;
    real* tabcutoff1;
    real* tabcutoff2;
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
	integer*  atomv; /* after PBC are applied, grid points from single
			  * batch can end up belonging to different atoms
			  * (i.e. same atom but in different cells). */
	/*ssf preprocessing: temp arrays and variables */
	integer* process;        /* logical array saying whether the weight of 
                                  * the point is set to 1.0 and no further 
                                  * processing needed or not*/
	integer* relevant_atoms;
	integer rel_atom_cnt;
	integer atomi;
        integer firstpoint;
};
typedef struct GridGenWork_ GridGenWork;

/* atom grid init */
	static GridGenAtomGrid*
agrid_new(integer uniq_no, integer Z)
{
	GridGenAtomGrid* grid = dal_new(1,GridGenAtomGrid);
	grid->uniq_no = uniq_no;
	grid->Z       = Z;
	grid->pnt     = 0;
	grid->leb_ang = NULL;
	grid->rad     = NULL;
	grid->wght    = NULL;
	return grid;
}

	static void
agrid_set_radial(GridGenAtomGrid* grid, unsigned cnt)
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
agrid_free(GridGenAtomGrid* grid)
{
	free(grid->leb_ang);
	free(grid->rad);
	free(grid->wght);
	free(grid);
}

struct CellBlocksData {
	real cell_center[3];
	const real (*shell_center)[3];
	const real *rshell;  /**< _grid_ cell radii */
	real celldg;   /**< half of the _grid_ cell main diagonal */
	integer kmax;      /**< shell count */
	integer *nblocks;
	integer (*shlblocks)[2];
	integer *s_idx; /**< 3D array (maxl,maxl,maxl) containing indexes to
			  nblocks and shlblocks. */
	integer used_entries;
};
#define I2OFF(i,ml) (i[0]+ml+(2*ml+1)*(i[1]+ml+(2*ml+1)*(i[2]+ml)))

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
	static integer
gc2_rad_cnt(integer Z, real thrl)
{
	static const integer MIN_RAD_PT = 20;
	integer ta, ri;
	if(Z<=2) ta=0;
	else if(Z<=10) ta=1;
	else if(Z<=18) ta=2;
	else if(Z<=36) ta=3;
	else if(Z<=54) ta=4;
	else if(Z<=86) ta=5;
	else ta=6;

	ri = rint( -3.2*(1.8*log10(thrl)-ta*4.0 ) );
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
	integer i;

	agrid_set_radial(grid, gc2_rad_cnt(grid->Z,thrl));
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
	integer  *nucorb;
	real *aa;
	integer maxl;
};

	static void
gen_lmg_free(void *quad_data)
{
	struct lmg_data *lmg = (struct lmg_data*)quad_data;

	free(lmg->nucorb);
	free(lmg->aa);
	free(lmg);
}

/* Treutler-Ahlrichs scheme [JCP 102, 346 (1995)] */
	static integer
turbo_get_no_of_points(integer Z, real thrl)
{
	integer ta, accuracy_correction, z_correction, gridSize;

	if(Z<=2) ta=0;
	else if(Z<=10) ta=1;
	else if(Z<=18) ta=2;
	else if(Z<=36) ta=3;
	else if(Z<=54) ta=4;
	else if(Z<=86) ta=5;
	else ta=6;

	/* thrl = 1e-5 maps to 0, 1e-13 -> 25, following Table III */
	if ((-log10(thrl)-5.0)*3.0 < 0.0) {dalton_quit("Negative threshold not allowed in turbo_get_no_of_points");}
	accuracy_correction = (integer) ( (-log10(thrl)-5.0)*3.0  + 0.5 );
	if(accuracy_correction<0) accuracy_correction = 0;

	z_correction = ta*5;

	static const integer MIN_RAD_PT = 20;

	gridSize = MIN_RAD_PT + accuracy_correction + z_correction;

	return gridSize;
}

	static void 
turbo_quad_generate(GridGenAtomGrid* grid, real thrl, void *quad_data)
{
	/* static const real DALTON_FUDGE_FACTOR = 1e4; */
	static const real zetas[] =  {
		/* H */ 0.8, /* He */ 0.9,
		/* Li */ 1.8, /* Be */ 1.4, /* B */ 1.3,  /* C */ 1.1,
		/* N */ 0.9,  /* O */  0.9, /* F */ 0.9,  /* Ne */ 0.9,
		/* Na */ 1.4, /* Mg */ 1.3, /* Al */ 1.3, /* Si */ 1.2,
		/* P */ 1.1,  /* S */  1.0, /* Cl */ 1.0, /* Ar */ 1.0,
		/* K */ 1.5,  /* Ca */ 1.4, /* Sc */ 1.3, /* Ti */ 1.2, /* V */ 1.2, 
		/* Cr */ 1.2, /* Mn */ 1.2, /* Fe */ 1.2, /* Co */ 1.2, /* Ni */ 1.1,
		/* Cu */ 1.1, /* Zn */ 1.1, /* Ga */ 1.1, /* Ge */ 1.0, /* As */ 0.9,
		/* Se */ 0.9, /* Br */ 0.9, /* Kr */ 0.9
	};
	real piOverN, zeta, rfac;
	static const real a = 1.0;
	integer i;

	zeta = grid->Z >=1 && grid->Z <= (integer)(sizeof(zetas)/sizeof(zetas[0]))
		?  zetas[grid->Z-1] : 0.9;
	rfac = zeta/M_LN2;

	/* radial points */ 
	agrid_set_radial(grid, turbo_get_no_of_points(grid->Z, thrl));
	piOverN = M_PI/grid->pnt;
	for (i=0; i<grid->pnt; i++) {
		real angle = (i+0.5)*piOverN;
		real x = cos(angle);
		real s = sin(angle);
		real w = piOverN * s;
		real aPlusX06 = pow(a+x, 0.6);
		real logAPlus1Over1MinusX = log( (a+1.0)/(1.0-x) );
		real r = rfac*aPlusX06*logAPlus1Over1MinusX;
		real rdiff = rfac*(aPlusX06/(1.0-x) +
				0.6*logAPlus1Over1MinusX/pow(a+x,0.4));
		grid->wght[i] = w*rdiff*r*r;
		grid->rad[i] = r;
	}
}

struct radial_scheme_t {
	char *name;
	void *(*quad_init)(void);
	void (*quad_gen)  (GridGenAtomGrid*, real thrl, void *quad_data);
	void (*quad_free) (void *quad_data);
};

extern void FSYM(ii_nucbas)(integer*, real* , const integer *NCENT,
			    const integer *NHKT, const integer *NUCO,
			    integer *NHTYP, integer *NUCIND, integer *KMAX,
			    integer *MXPRIM, const real *PRIEXP,
			    const integer *NSTART, const integer*);

extern void FSYM(ii_radlmg)(real*rad, real* wght, integer *nr, real* raderr,
                         const integer *maxrad, integer *nucorb, real *aa,
			    integer *nhtyp, const integer *);

static void*
gen_lmg_init(void)
{
	integer nucind;
	struct lmg_data *lmg = dal_new(1, struct lmg_data);
	dalton_quit("gen_lmg_init should not be solved call gen_lmg_init2");
	return lmg;
}

static void*
gen_lmg_init2(const integer *NCENT, const integer *NHKT, const integer *NUCO,
	     integer NHTYP, integer NUCIND,
	     integer KMAX, integer MXPRIM, const real *PRIEXP, const integer *NSTART)
{
    /*    integer nucind;*/
    struct lmg_data *lmg = dal_new(1, struct lmg_data);
    lmg->maxl = NHTYP;
    lmg->nucorb = malloc(2*lmg->maxl*NUCIND*sizeof(integer));
    lmg->aa     = calloc(4*lmg->maxl*NUCIND,sizeof(real));
    if(!lmg->nucorb|| !lmg->aa) {
        fprintf(stderr,"no enough memory. in gen_lmg_init.\n");
        exit(1);
    }
    FSYM(ii_nucbas)(lmg->nucorb, lmg->aa, NCENT,NHKT,NUCO,&NHTYP,
		    &NUCIND,&KMAX,&MXPRIM,PRIEXP,NSTART,&ONEI);
    return lmg;
}

	static void
gen_lmg_quad(GridGenAtomGrid* grid, real thrl, void *quad_data)
{
  dalton_quit("gen_lmg_quad should not be solved call gen_lmg_quad2");
}

static void
gen_lmg_quad2(GridGenAtomGrid* grid, real thrl, void *quad_data,
	     const integer *NCENT, const integer *NHKT, const integer *NUCO,
	     integer NHTYP, integer NUCIND, integer KMAX, integer MXPRIM, const real *PRIEXP,
	     const integer *NSTART)
{
    static const integer MAXRAD = 2000;
    struct lmg_data *lmg = (struct lmg_data*)quad_data;
    size_t new_mem_size;
    integer ptcnt;
    agrid_set_radial(grid, MAXRAD);
    FSYM(ii_radlmg)(grid->rad, grid->wght, &ptcnt, &thrl, &MAXRAD,
                 lmg->nucorb+2*grid->uniq_no*lmg->maxl,
                 lmg->aa    +4*grid->uniq_no*lmg->maxl,
		    &NHTYP,&ZEROI);
    grid->pnt = ptcnt;
    grid->rad  = realloc(grid->rad, grid->pnt * sizeof(real));
    grid->wght = realloc(grid->wght, grid->pnt * sizeof(real));
}

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

static struct radial_scheme_t quad_turbo = { 
	"Treutler-Ahlrichs M4-T2 scheme",
	NULL,
	turbo_quad_generate,
	NULL
};

/* ------------------------------------------------------------------- */
/*              The angular grid generation interface                  */
/* ------------------------------------------------------------------- */

/* routines for generation of Lebedev grid */
void ld0006_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0014_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0026_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0038_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0050_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0074_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0086_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0110_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0146_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0170_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0194_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0230_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0266_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0302_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0350_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0434_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0590_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0770_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0974_(real* x, real* y, real* z, real* w, integer *pnt);
void ld1202_(real* x, real* y, real* z, real* w, integer *pnt);
void ld1454_(real* x, real* y, real* z, real* w, integer *pnt);

/* some of the point type values were guessed but since it was only used
   for ANGMIN -> quadrature mapping type, the accuracy is not really relevant.
   Anyone, who would rely on particular mapping here would be crazy.
   The true values can be found in literature.
 */

static struct leb_gen_{
	integer point_cnt, poly_order;
	void (*func)(real* x, real* y, real* z, real* w, integer *pnt);
} leb_gen[] = {
	{ 14,   5, ld0014_},
	{  38,  9, ld0038_},
	{  50, 11, ld0050_},
	{  86, 15, ld0086_},
	{ 110, 17, ld0110_},
	{ 146, 19, ld0146_},
	{ 170, 21, ld0170_},
	{ 194, 23, ld0194_},
	{ 230, 25, ld0230_},
	{ 266, 27, ld0266_},
	{ 302, 29, ld0302_},
	{ 350, 31, ld0350_},
	{ 434, 35, ld0434_}, 
	{ 590, 41, ld0590_},
	{ 770, 47, ld0770_},
	{ 974, 53, ld0974_},
	{1202, 59, ld1202_},
	{1454, 64, ld1454_},
};  

/* leb_get_from_point returns index to leb_gen array that points to
   the entry having at least iang points.
 */
	static integer
leb_get_from_point(integer iang)
{
	integer i;
	for(i=ELEMENTS(leb_gen)-1; i>=0; i--)
		if(iang>=leb_gen[i].point_cnt) return i;
	return 0;
}

/* leb_get_from_order returns index to leb_gen array that points to
   the entry contructing polyhedra of at least given order.
 */
	static integer
leb_get_from_order(integer poly_order)
{
	integer i;
	for(i=0; i<ELEMENTS(leb_gen); i++)
		if(leb_gen[i].poly_order>=poly_order)
			return i;
	return 0;
}

/* Trond Saue:
   The below data gives atomic radii in Angstroms and stems from table I of 
   J.C.Slater: "Atomic Radii in Crystals"
   J.Chem.Phys. 41(1964) 3199-3204
   Values for elements marked with an asterisk has been
   guessed/interpolated
 */
static const real bragg_radii[] = {
	/*     dummy         */
	0.75,
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
/* max number of points per atom per angular shell */
#define MAX_PT_PER_SHELL (leb_gen[ELEMENTS(leb_gen)-1].point_cnt)

/* include_partitioning_weights:
   scale weights
 */
#define IDX(arr,i,j) arr[(i)*MAX_PT_PER_SHELL +(j)]
	static void
gridgen_compute_rjs(GridGenWork *ggw, GridGenMolGrid* mg,
		integer point_cnt, integer idx)
{
	integer atno, ptno;

	for(atno=0; atno < mg->atom_cnt; atno++) {
		for(ptno=0; ptno < point_cnt; ptno++) {
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

/* include_partitioning_becke_corr:
   multiply current sphere/shell generated for specified atom by space
   partitioning weights as described by Becke.
   The computed weights are somehow redundant: unnormalized
   weights p_kg are computed for all atoms, althougth only one is 
   needed at this stage. The algorithm should probably be restructured
   to avoid this.
 */

/*becke scheme without atomic size correction*/
	static void
becke_orig_preprocess(GridGenMolGrid* mg, integer point_cnt,
		GridGenWork *ggw, integer idx, integer verbose,integer skip)
{
	integer atno, atno2, ptno, h;
	integer hardness=mg->conf->hardness;
	real mua;
        real *mu =dal_new(point_cnt,real);
        real *mu2 =dal_new(point_cnt,real);
	integer pc;

	if(mg->atom_cnt==1) {
	  free(mu);
	  free(mu2);
	  return;
	}
        /* compute confocal ellipical coordinates for the batch of points to
         * be processed. */
	gridgen_compute_rjs(ggw, mg, point_cnt, idx);

	if(verbose)printf("Computing cell functions for atom %d\n", (int)*ggw->atomv);

        /*compute cell functions*/
	for(atno=1; atno<mg->atom_cnt; atno++) {
		if(mg->atom_grids[atno]->pnt == 0)
			continue;
		for(atno2=0; atno2<atno; atno2++) {
			real dist = mg->rij[atno2+(atno*(atno-1))/2];
			if(mg->atom_grids[atno2]->pnt == 0)
				continue;

                        for(ptno=0; ptno<point_cnt; ptno++) {
                           mua  =(ggw->IDX(rj,atno,ptno)-ggw->IDX(rj,atno2,ptno))*dist;
                           if(mua<-1) mua= -1; /* numerical error correction, needed? */
                           if(mua>1)  mua=  1; /* numerical error correction, needed? */
                           mu[ptno]=mua;
                           mu2[ptno]=mua*mua;
                        }
                        for(h=0;h<hardness;h++){
                           for(ptno=0; ptno<point_cnt; ptno++) {
                              mu[ptno]=0.5*mu[ptno]*(3.0-mu2[ptno]);
                           }
                           if(h<hardness-1){
                              for(ptno=0; ptno<point_cnt; ptno++) mu2[ptno]=mu[ptno]*mu[ptno];
                           }
                        }
                        for(ptno=0; ptno<point_cnt; ptno++) {
                           ggw->IDX(p_kg,atno,ptno)  *= 0.5*(1.0-mu[ptno]);
                           ggw->IDX(p_kg,atno2,ptno) *= 0.5*(1.0+mu[ptno]);
                        }
		}
	}

	/* compute weight normalization factors */
	/*ls_dzero_(ggw->vec, &point_cnt);*/
	for(ptno=0; ptno<point_cnt; ptno++) {
         ggw->vec[ptno] = 0.0000000000000000; 
	}
	for(atno=0; atno<mg->atom_cnt; atno++) {
	  /*
	    We exclude centers which carry no grid points
	    these are point charges without basis sets
	    they should not be included in the partitioning
	   */
		if(mg->atom_grids[atno]->pnt == 0)
			continue;
		for(ptno=0; ptno<point_cnt; ptno++)
			ggw->vec[ptno] += ggw->IDX(p_kg,atno,ptno); 
	}

	/*
	 *   Apply the computed weights
	 */
	for(ptno=0; ptno<point_cnt; ptno++) {
		ggw->wg[idx+ptno] *= ggw->IDX(p_kg,ggw->atomv[ptno],ptno)/ggw->vec[ptno];
	}

        free(mu);
        free(mu2);
}

/*becke scheme with atomic size correction*/
static void
becke_corr_preprocess(GridGenMolGrid* mg, integer point_cnt,
                      GridGenWork *ggw, integer idx, integer verbose,integer skip)
{
    integer atno, atno2, ptno, h;
    integer hardness=mg->conf->hardness;
    real mua;
    real *mu =dal_new(point_cnt,real);
    real *mu2 =dal_new(point_cnt,real);
    integer pc;

    if(mg->atom_cnt==1)
        return;

    /* compute confocal ellipical coordinates for the batch of points to
     * be processed. */
    gridgen_compute_rjs(ggw, mg, point_cnt, idx);

    /*compute cell functions*/
    if(verbose)printf("Computing cell functions for atom %d\n", (int)*ggw->atomv);
    for(atno=1; atno<mg->atom_cnt; atno++) {
        if(mg->atom_grids[atno]->pnt == 0)
            continue;
        for(atno2=0; atno2<atno; atno2++) {
            real dist = mg->rij[atno2+(atno*(atno-1))/2];
            real bfac = mg->aij[atno2+(atno*(atno-1))/2];
            if(mg->atom_grids[atno2]->pnt == 0)
                continue;
            for(ptno=0; ptno<point_cnt; ptno++) {
                mua  =(ggw->IDX(rj,atno,ptno)-ggw->IDX(rj,atno2,ptno))*dist;
                mua += bfac*(1-mua*mua);
                if(mua<-1) mua= -1; /* numerical error correction, needed? */
                if(mua>1)  mua=  1; /* numerical error correction, needed? */
                mu[ptno]=mua;
                mu2[ptno]=mua*mua;
            }
            for(h=0;h<hardness;h++){
               for(ptno=0; ptno<point_cnt; ptno++) {
                  mu[ptno]=0.5*mu[ptno]*(3.0-mu2[ptno]);
               }
               if(h<hardness-1){
                  for(ptno=0; ptno<point_cnt; ptno++) mu2[ptno]=mu[ptno]*mu[ptno];
               }
            }
            for(ptno=0; ptno<point_cnt; ptno++) {
               ggw->IDX(p_kg,atno,ptno)  *= 0.5*(1.0-mu[ptno]);
               ggw->IDX(p_kg,atno2,ptno) *= 0.5*(1.0+mu[ptno]);
            }
        }
    }
    /* compute weight normalization factors */
    pc = point_cnt; 
    /*ls_dzero_(ggw->vec, &pc);*/
    for(ptno=0; ptno<point_cnt; ptno++) {
      ggw->vec[ptno] = 0.0000000000000000; 
    }
    for(atno=0; atno<mg->atom_cnt; atno++) {
        if(mg->atom_grids[atno]->pnt == 0)
            continue;
        for(ptno=0; ptno<point_cnt; ptno++)
            ggw->vec[ptno] += ggw->IDX(p_kg,atno,ptno); 
    }
    /*
     *   Apply the computed weights
     */
    for(ptno=0; ptno<point_cnt; ptno++){
        ggw->wg[idx+ptno] *= ggw->IDX(p_kg,ggw->atomv[ptno],ptno)/ggw->vec[ptno];
    }
    
     free(mu);
     free(mu2);
}

/* include_partitioning_ssf:
   multiply current sphere/shell by space partitioning weights as
   described in SSF article. (Chem. Phys. Lett. 1996, 257, 213)
*/
static void
ssf_preprocess(GridGenMolGrid* mg, integer point_cnt,
               GridGenWork *ggw, integer idx, integer verbose,integer skip)
{

    if(skip) return;

    /*some parameters ... */
    real fac0 = mg->ssf_fac0;
    real fac1 = mg->ssf_fac1;
    real fac2 = mg->ssf_fac2;
    real fac3 = mg->ssf_fac3;
    real fac4 = mg->ssf_fac4;
    
    integer atno, ptno;
    integer aidx,atno1i;
    real mu, g_mu, w,v;
    real pr, rx, ry, rz, mu2;
    real *p = dal_new(mg->atom_cnt,real);
    integer *relevant_atoms = dal_new(mg->atom_cnt, integer);
    integer atomi=0, ati, ati2, rel_atom_cnt = 0, atom=ggw->atomv[idx];
    integer pc;

    rx = ggw->x[idx] - mg->atom_coords[atom].x;
    ry = ggw->y[idx] - mg->atom_coords[atom].y;
    rz = ggw->z[idx] - mg->atom_coords[atom].z;
    pr = sqrt(rx*rx+ry*ry+rz*rz); /* radius of the sphere */

    /*
     * determine the relevant atoms:
     *    we make use of the fact that the outermost points come first, 
     *    so that the list of relevant atoms becomes shorter in the course of evaluation
     */

    /*do the outermost points really come first?*/
    assert(mg->atom_grids[atom]->pnt>=2); 
    assert(mg->atom_grids[atom]->rad[0] > mg->atom_grids[atom]->rad[1]);
    /*for the outermost point we use the original list*/
    if(ggw->firstpoint == 0) {
       for(atno=0; atno<atom; atno++) {
           real dist = mg->rij[atno+(atom*(atom-1))/2];
           if(mg->atom_grids[atno]->pnt > 0 && 2*pr*dist>1-SSF_CUTOFF){
               relevant_atoms[rel_atom_cnt++] = atno;}
       }
       atomi = rel_atom_cnt;
       relevant_atoms[rel_atom_cnt++] = atom;
       for(atno=atom+1; atno<mg->atom_cnt; atno++) {
           real dist = mg->rij[atom+(atno*(atno-1))/2];
           if(mg->atom_grids[atno]->pnt > 0 && 2*pr*dist>1-SSF_CUTOFF){
               relevant_atoms[rel_atom_cnt++] = atno;}
       }
       ggw->rel_atom_cnt = rel_atom_cnt;
       ggw->atomi = atomi;
       integer i=0;
       for(i=0;i<rel_atom_cnt;i++) ggw->relevant_atoms[i]=relevant_atoms[i];
       ggw->firstpoint=1;
    }
    /*for the inner points the list of rel atoms in shorter ...*/
    else{
       rel_atom_cnt=0;
       for(ati=0; ati<ggw->rel_atom_cnt; ati++) {
          integer atno = ggw->relevant_atoms[ati];
          if(atno==atom) {
             atomi=rel_atom_cnt;
             relevant_atoms[rel_atom_cnt++] = atom;
             continue;
          }
          else if(atno<atom){atno1i=atno+(atom*(atom-1))/2;}
          else{atno1i=atom+(atno*(atno-1))/2;}
          real dist = mg->rij[atno1i];
          if(mg->atom_grids[atno]->pnt > 0 && 2*pr*dist>1-SSF_CUTOFF){
            relevant_atoms[rel_atom_cnt++] = atno;
          }
       }
       ggw->rel_atom_cnt = rel_atom_cnt;
       ggw->atomi = atomi;
       integer i=0;
       for(i=0;i<rel_atom_cnt;i++) ggw->relevant_atoms[i]=relevant_atoms[i];
    }
    if(rel_atom_cnt==1) {
        free(relevant_atoms);
        return;
    }


    /*       
     * compute confocal ellipical coordinates for the batch of points to
     * be processed. 
     */
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

    if(verbose)printf("Computing cell functions for atom %d\n", (int)atom);

    /* 
     * evaluate the partition functions
     */

    /*VERSION 1: grid point driven 
     *           seems to be the most efficient one, 
     *           due to very efficient "screening" \Andreas Krapp*/
    for(ptno=0; ptno<point_cnt; ptno++) {
       for(ati=0; ati<rel_atom_cnt; ati++)  p[ati]=1.0;
       /*we first sort out all values with zero weight*/
       atno1i=atom*(atom-1)/2;
       for(ati=0; ati<rel_atom_cnt; ati++) {
          if(ati==atomi) continue;
          integer atno2 = relevant_atoms[ati];
          if(ati<atomi) aidx = atno2+atno1i;
          else aidx = atom+atno2*(atno2-1)/2;
          real dist = mg->rij[aidx];
          mu = (pr-ggw->IDX(rj,ati,ptno))*dist;
          if(mu>=SSF_CUTOFF){
             ggw->wg[idx+ptno]=0.0;
             w = 0.0;
             goto next_point;}
       }
       /*now for the rest of points*/
       for(ati=1; ati<rel_atom_cnt; ati++) {
          atno = relevant_atoms[ati];
          atno1i=atno*(atno-1)/2;
          real ria1 = ggw->IDX(rj,ati,ptno);
          for(ati2=0; ati2<ati; ati2++) {
             integer atno2 = relevant_atoms[ati2];
             real dist = mg->rij[atno2+atno1i];
             mu = (ria1-ggw->IDX(rj,ati2,ptno))*dist;
             if(mu<=-SSF_CUTOFF) {
                p[ati2] = 0.0; }
             else if(mu>=SSF_CUTOFF){
                p[ati] = 0.0; }
             else {
                 mu2 = mu*mu;
                 g_mu=mu*(fac1+mu2*(fac2+mu2*(fac3+fac4*mu2)));
                 p[ati]  *= 0.5-g_mu;
                 p[ati2] *= 0.5+g_mu;
             }
          }
       }
       if(p[atomi]==0.0){
          ggw->wg[idx+ptno] = 0.0;
          w = 0.0;
       }
       else{
          /* compute weight normalization factors */
          v = 0.0;
          for(ati=0; ati<rel_atom_cnt; ati++)  v += p[ati];
          w = p[atomi]/v;
          /* Apply the computed weights     */
          ggw->wg[idx+ptno]*=w;
       }
       next_point:;
    }

    free(relevant_atoms);
    free(p);
}

/** include_block_partitioning transforms a set of POINT_CNT points.
*/
typedef struct {
    real (*coor)[4];
    real *rj;
    real *p_kg;
    real *vec;
    integer LDA; /* leading dimension of rj and p_kg */
} GGBlockWork;

static void
block_work_init(GGBlockWork* ggw, real (*coorw)[4],
                integer atom_cnt, integer point_cnt)
{
    ggw->coor = coorw;
    ggw->LDA  = point_cnt;
    ggw->rj   = dal_new(atom_cnt*point_cnt, real);
    ggw->p_kg = dal_new(atom_cnt*point_cnt, real);
    ggw->vec  = dal_new(point_cnt, real);
}

static void
block_work_release(GGBlockWork* ggw)
{
    free(ggw->rj);
    free(ggw->p_kg);
    free(ggw->vec);
}

/*changed loops for speed*/
static void 
block_compute_rjs(GGBlockWork *ggw, GridGenMolGrid* mg,
                  integer rel_at_cnt, integer *relevant_atoms,
                  integer point_cnt)
{                 
    integer ptno, i;  
    real *dx = dal_new(point_cnt, real);
    real *dy = dal_new(point_cnt, real);
    real *dz = dal_new(point_cnt, real);
    for(i=0; i<rel_at_cnt; i++) {
        integer atno = relevant_atoms[i];
        real x = mg->atom_coords[atno].x;
        real y = mg->atom_coords[atno].y; 
        real z = mg->atom_coords[atno].z;
        integer posi = i*ggw->LDA; 
        for(ptno=0; ptno<point_cnt; ptno++) {
            dx[ptno  ] = x - ggw->coor[ptno][0];
            dy[ptno  ] = y - ggw->coor[ptno][1];
            dz[ptno  ] = z - ggw->coor[ptno][2];
        }
        for(ptno=0; ptno<point_cnt; ptno++) {
            real r2x = dx[ptno]*dx[ptno]; 
            real r2y = dy[ptno]*dy[ptno];
            real r2z = dz[ptno]*dz[ptno];
            ggw->p_kg[posi+ptno] = 1.0;
            ggw->rj  [posi+ptno] = sqrt(r2x+r2y+r2z);
        }   
    }          
    free(dx);
    free(dy);
    free(dz);
}

/* SSF scheme for grid weight evaluation 
 * combined with a blockwise handling of grid points 
 * the later be described in JCP, 2004, 171, 2915.
 *
 * the later allows to reduce the list of relevant atom pairs 
 * considerably 
 * \Andreas Krapp, 07/2010
 */
static integer
block_postprocess_ssf(GridGenMolGrid *mg, const real *center,
                  integer point_cnt, const integer *atom_nums,
                  real (*coorw)[4], const integer *nprocess)
{
    /*some parameters ...*/
    real fac0 = mg->ssf_fac0;
    real fac1 = mg->ssf_fac1;
    real fac2 = mg->ssf_fac2;
    real fac3 = mg->ssf_fac3;
    real fac4 = mg->ssf_fac4;
    /**/
    real local_threshold = WEIGHT_THRESHOLD*1e-5;
    integer atno,atno2,ptno,i,j,pc,dest,dest2,rel_pt_cnt;
    integer atom,atomi,aidx,ati,rptno,posi,posj,posat,posipt,posjpt;
    integer uniq_atoms=0;
    real mu2;
    real mu, g_mu;
    real dist, r1, v, w;
    integer *relevant_atoms = dal_new(mg->atom_cnt, integer);
    integer *map2r = dal_new(mg->atom_cnt, integer);
    GGBlockWork ggw;

    /* we find first atoms that relevant for this cell.
     * We do it by linear search which will scale as N^2
     * but we can live with that for now.
     */
    uniq_atoms = 0;
    for(i=0; i<mg->atom_cnt; i++) map2r[i] = -1;
    for(atno=0; atno<mg->atom_cnt; atno++) {
        GridGenAtomGrid *ag = mg->atom_grids[atno];
        real dx = mg->atom_coords[atno].x - center[0];
        real dy = mg->atom_coords[atno].y - center[1];
        real dz = mg->atom_coords[atno].z - center[2];
        real dist2 = dx*dx + dy*dy + dz*dz;
        /* Increase the radius beyond the last point to let the grid
         * time to end.  Adding a scaled-up difference between last
         * two radial points is possibly the best way out. */
        real r;
        if(ag->pnt == 0)
            continue;
        assert(ag->pnt>=2);
        assert(ag->rad[0]> ag->rad[1]);
        /*ANDREAS:should it not be CELL_SIZE*0.5*sqrt(3.0) instead of CELL_SIZE*/
        r = ag->rad[0] + CELL_SIZE;
        if(r*r>dist2) {
            map2r[atno] = uniq_atoms;
            relevant_atoms[uniq_atoms++] = atno;
        }
    }
    for(i=0; i<point_cnt; i++) {
        integer atnoi = atom_nums[i];
        if(map2r[atnoi] == -1) {
            map2r[atnoi] = uniq_atoms;
            relevant_atoms[uniq_atoms++] = atnoi;
            printf("Internal safety check corrected for atom %d - "
                       "please report.", (int)atnoi);
        }
    }

    /*error checking*/
    if(uniq_atoms<=1) { /* 0 cannot happen and 1 - no partitioning. */
        free(relevant_atoms);
        free(map2r);
        return point_cnt;
    }

    /*allocate memory and initialize*/
    block_work_init(&ggw, coorw, uniq_atoms, point_cnt);

    /*compute distance between grid points and relevant atoms*/
    block_compute_rjs(&ggw, mg, uniq_atoms, relevant_atoms,
                      point_cnt);

    /*compute partition functions*/
#if 1
#if 1
   /* THIS IS A VERSION WHICH FIRST GOES OVER INDIVUDAL POINTS 
    * TO SORT OUT WEIGHTS OF 1.0 AND 0.0
    * AND THEN GOES OVER THE REMAINING POINTS 
    */
    dest2=0;
    dest=0;
    rel_pt_cnt=0;
    integer *rel_points = dal_new(point_cnt, integer);
    real *tmp   = dal_new(point_cnt*4,real);
    real *mua   = dal_new(point_cnt,  real);
    real *g_mua = dal_new(point_cnt,  real);

    /* 1. sort out grid points of 1.0 and 0.0 weight*/
    for(ptno=0; ptno<point_cnt; ptno++) {
       if(nprocess[ptno]==1) {
            /*for points with weight 1.0 nothing has to be done, shift them to a temp. array*/
            tmp[dest2  ] = coorw[ptno][0];
            tmp[dest2+1] = coorw[ptno][1];
            tmp[dest2+2] = coorw[ptno][2];
            tmp[dest2+3] = coorw[ptno][3];
            dest2=dest2+4;
            goto next_point1;
       } 
       else{
          /*throw out all points with zero weight*/
          atom = atom_nums[ptno];
          atomi = map2r[atom];
          for(i=0; i<uniq_atoms; i++) {
             if(i==atomi) continue;
             atno2 = relevant_atoms[i];
             if(i<atomi) aidx = atno2+atom *(atom -1)/2;
             else        aidx = atom +atno2*(atno2-1)/2;
             dist = mg->rij[aidx];
             mu = (ggw.rj[ggw.LDA*atomi+ptno]-ggw.rj[ggw.LDA*i+ptno])*dist;
             if(mu>=SSF_CUTOFF) goto next_point1;
          }
       }
       /*the actual point has neither weight of 1.0 nor 0.0 so 
        * add it to the list of relevant points*/
       rel_points[rel_pt_cnt++]=ptno;
       next_point1:;
    }

    /* 2. compute weights for the remaining points*/
    if(rel_pt_cnt > 0){
       /*compute partition functions*/
       for(i=1; i<uniq_atoms; i++) {
          atno = relevant_atoms[i];
          posi = ggw.LDA*i;
          for(j=0; j<i; j++) {
             atno2 = relevant_atoms[j];
             posat = atno2+(atno*(atno-1))/2;
             dist = mg->rij[posat];
             posj = ggw.LDA*j;
             for(rptno=0;rptno<rel_pt_cnt;rptno++){
                ptno=rel_points[rptno];
                posipt = posi+ptno;
                posjpt = posj+ptno;
                mua[rptno] = (ggw.rj[posipt]-ggw.rj[posjpt])*dist;
             }
             for(rptno=0;rptno<rel_pt_cnt;rptno++){
                if     (mua[rptno]<=-SSF_CUTOFF) g_mua[rptno] = -0.5; 
                else if(mua[rptno]>= SSF_CUTOFF) g_mua[rptno] =  0.5; 
                else{
                   mu2 = mua[rptno]*mua[rptno];
                   g_mua[rptno]=mua[rptno]*(fac1+mu2*(fac2+mu2*(fac3+fac4*mu2)));
                }
             }
             for(rptno=0;rptno<rel_pt_cnt;rptno++){
                ggw.p_kg[posi+rptno] *= 0.5-g_mua[rptno];
                ggw.p_kg[posj+rptno] *= 0.5+g_mua[rptno];
             }
          }
       }
       /* compute weight normalization factors */
       pc = rel_pt_cnt; 
       /*ls_dzero_(ggw.vec, &pc);*/
       for(ptno=0; ptno<rel_pt_cnt; ptno++) {
	 ggw.vec[ptno] = 0.0000000000000000; 
       }
       for(i=0; i<uniq_atoms; i++) {
           posi = ggw.LDA*i;
           for(rptno=0;rptno<rel_pt_cnt;rptno++){
              ggw.vec[rptno] += ggw.p_kg[posi+rptno];
           }
       }
       /* Apply the computed weights removing at the same time
        * points with low weight.  */
       for(rptno=0;rptno<rel_pt_cnt;rptno++){
            ptno=rel_points[rptno];
            atomi = map2r[atom_nums[ptno]];
            real factor = ggw.p_kg[ggw.LDA*atomi+rptno]/ggw.vec[rptno];
            real wg = coorw[ptno][3]*factor;
            if(fabs(wg) > local_threshold){
               coorw[dest][0] = coorw[ptno][0];
               coorw[dest][1] = coorw[ptno][1];
               coorw[dest][2] = coorw[ptno][2];
               coorw[dest][3] = wg;
               dest++;
            }
#if 0
            if(factor >= WEIGHT_THRESHOLD){
               coorw[dest][0] = coorw[ptno][0];
               coorw[dest][1] = coorw[ptno][1];
               coorw[dest][2] = coorw[ptno][2];
               coorw[dest][3] = coorw[ptno][3]*factor;
               dest++;}
#endif
       }
    }

    /* 3. now bring the points with weight 1.0 into place*/
    for(i=0;i<dest2;i=i+4){
       coorw[dest][0]=tmp[i  ];
       coorw[dest][1]=tmp[i+1];
       coorw[dest][2]=tmp[i+2];
       coorw[dest][3]=tmp[i+3];
       dest++;
    }

    free(rel_points);
    free(tmp);
    free(mua);
    free(g_mua);
#else
   /* THIS IS A COMPLETELY POINT DRIVEN VERSION */
    real *p = dal_new(mg->atom_cnt, real);
    dest=0;
    for(ptno=0; ptno<point_cnt; ptno++) {
       integer LDBptno=ggw.LDB*ptno;
       atom = atom_nums[ptno];
       atomi = map2r[atom];
       /*all points which have a weight of 1.0 do not have to be touched 
         save directly*/
       if(nprocess[ptno]==1) {
            coorw[dest][0] = coorw[ptno][0];
            coorw[dest][1] = coorw[ptno][1];
            coorw[dest][2] = coorw[ptno][2];
            coorw[dest][3] = coorw[ptno][3];
            dest++;
            w = 1.0;
            goto next_point;
       }
       /*throw out all points with zero weight*/
       for(i=0; i<uniq_atoms; i++) {
          if(i==atomi) continue;
          atno2 = relevant_atoms[i];
          if(i<atomi) aidx = atno2+atom *(atom -1)/2;
          else        aidx = atom +atno2*(atno2-1)/2;
          dist = mg->rij[aidx];
          mu = (ggw.rj[ggw.LDA*atomi+ptno]-ggw.rj[ggw.LDA*i+ptno])*dist;
          if(mu>=SSF_CUTOFF){
             w = 0.0;
             goto next_point;}
       }
       /*now comes the rest of points*/
       for(i=0; i<uniq_atoms; i++)  p[i]=1.0;
       for(i=1; i<uniq_atoms; i++) {
          atno = relevant_atoms[i];
          ati  = atno*(atno-1)/2;
          r1 = ggw.rj[ggw.LDA*i+ptno];
          for(j=0;j<i;j++){
             atno2 = relevant_atoms[j];
             dist = mg->rij[atno2+ati];
             mu = (r1-ggw.rj[ggw.LDA*j+ptno])*dist;
             if(mu<=-SSF_CUTOFF) {
                p[j] = 0.0; }
             else if(mu>=SSF_CUTOFF){
                p[i] = 0.0; }
             else if(p[i] > 0.0 || p[j] > 0.0)  {
                 mu2 = mu*mu;
                 g_mu=mu*(fac1+mu2*(fac2+mu2*(fac3+fac4*mu2)));
                 if(p[i] > 0.0) p[i] *= 0.5-g_mu;
                 if(p[j] > 0.0) p[j] *= 0.5+g_mu;
             }
          }
       }
       if(p[atomi]==0.0){
          w = 0.0;
       }
       else{
          /* compute weight normalization factors */
          v = 0.0;
          for(i=0; i<uniq_atoms; i++)  v += p[i];
          /* Apply the computed weights */
          w = p[atomi]/v;
          if(w >= WEIGHT_THRESHOLD) {
             coorw[dest][0] = coorw[ptno][0];
             coorw[dest][1] = coorw[ptno][1];
             coorw[dest][2] = coorw[ptno][2];
             coorw[dest][3] = coorw[ptno][3]*w;
             dest++;}
       }
       next_point:;
    }
    free(p);
#endif
#else
    /*"standard" version*/
    for(i=1; i<uniq_atoms; i++) {
       atno = relevant_atoms[i];
       posi = ggw.LDA*i;
       for(j=0; j<i; j++) {
          real dist;
          atno2 = relevant_atoms[j];
          posat = atno2+(atno*(atno-1))/2;
          dist = mg->rij[posat];
          posj = ggw.LDA*j;
          for(ptno=0; ptno<point_cnt; ptno++) {
             posipt = posi+ptno;
             posjpt = posj+ptno;
             mu[ptno] = (ggw.rj[posipt]-ggw.rj[posjpt])*dist;
          }
          for(ptno=0; ptno<point_cnt; ptno++) {
             if(mu[ptno]<=-SSF_CUTOFF) 
                 g_mu[ptno] = -0.5;
             else if(mu[ptno]>=SSF_CUTOFF) 
                 g_mu[ptno] = 0.5;
             else {
                 mu2 = mu[ptno]*mu[ptno];
                 g_mu[ptno]=mu[ptno]*(fac1+mu2*(fac2+mu2*(fac3+fac4*mu2)));
		 /* *mu[ptno] /= SSF_CUTOFF; mu2=mu[ptno]*mu[ptno]; */
                /* *g_mu  = 0.0625*mu*(35+mu2*(-35+mu2*(21-5*mu2)));*/
             }
          }
          for(ptno=0; ptno<point_cnt; ptno++) {
             ggw.p_kg[posi+ptno] *= 0.5-g_mu[ptno];
             ggw.p_kg[posj+ptno] *= 0.5+g_mu[ptno];
          }
       }
    } 

    /* compute weight normalization factors */
    pc = point_cnt; 
    /*ls_dzero_(ggw.vec, &pc);*/
    for(ptno=0; ptno<point_cnt; ptno++) {
      ggw.vec[ptno] = 0.0000000000000000; 
    }
    for(i=0; i<uniq_atoms; i++) {
        posi = ggw.LDA*i;
        for(ptno=0; ptno<point_cnt; ptno++)
            ggw.vec[ptno] += ggw.p_kg[posi+ptno];
    }
    /*
     * Apply the computed weights removing at the same time
     * points with low weight.
     */
    for(dest=0, ptno=0; ptno<point_cnt; ptno++) {
            atom = map2r[atom_nums[ptno]];
            real factor = ggw.p_kg[ggw.LDA*atom+ptno]/ggw.vec[ptno];
            if(factor >= WEIGHT_THRESHOLD) {
               coorw[dest][0] = coorw[ptno][0];
               coorw[dest][1] = coorw[ptno][1];
               coorw[dest][2] = coorw[ptno][2];
               coorw[dest][3] = coorw[ptno][3]*factor;
               dest++;
            }
    }
#endif

    free(relevant_atoms);
    free(map2r);

    block_work_release(&ggw);

    return dest;
}

/* becke scheme for getting partition functions 
 * for a block of points, allowing to reduce the number of atoms
 * to be considered in the evaluation, leading to linear scaling
 * see description in JCP, 2004, 171, 2915.
  */
static integer
block_postprocess(GridGenMolGrid *mg, const real *center,
                  integer point_cnt, const integer *atom_nums,
                  real (*coorw)[4], const integer *nprocess)
{
    integer atno, ptno, h, i, j, atno2, a,b;
    integer hardness=mg->conf->hardness;
    real mu3;
    real *mu   = dal_new(point_cnt,real);
    real *mu2  = dal_new(point_cnt,real);
    real *g_mu = dal_new(point_cnt,real);
    /* compute confocal ellipical coordinates for the batch of points to
     * be processed. */
    integer *relevant_atoms = dal_new(mg->atom_cnt, integer);
    /* map2r: map from all to relevant_atoms array */
    integer *map2r = dal_new(mg->atom_cnt, integer);
    integer uniq_atoms;
    GGBlockWork ggw;
    integer dest;
    integer pc;

    /* we find first atoms that relevant for this cell.
     * We do it by linear search which will scale as N^2
     * but we can live with that for now.
     */
    uniq_atoms = 0;
    for(i=0; i<mg->atom_cnt; i++) map2r[i] = -1;
    for(atno=0; atno<mg->atom_cnt; atno++) {
        GridGenAtomGrid *ag = mg->atom_grids[atno];
        real dx = mg->atom_coords[atno].x - center[0];
        real dy = mg->atom_coords[atno].y - center[1];
        real dz = mg->atom_coords[atno].z - center[2];
        real dist2 = dx*dx + dy*dy + dz*dz;
        /* Increase the radius beyond the last point to let the grid
         * time to end.  Adding a scaled-up difference between last
         * two radial points is possibly the best way out. */
        real r;
        if(ag->pnt == 0)
            continue;
        assert(ag->pnt>=2);
        assert(ag->rad[0]> ag->rad[1]);
        /*ANDREAS:should it not be CELL_SIZE*0.5*sqrt(3.0) instead of CELL_SIZE*/
        r = ag->rad[0] + CELL_SIZE;
        if(r*r>dist2) {
	    map2r[atno] = uniq_atoms;
	    relevant_atoms[uniq_atoms++] = atno;
        }
    }
    for(i=0; i<point_cnt; i++) {
        integer atnoi = atom_nums[i];
	if(map2r[atnoi] == -1) {
            map2r[atnoi] = uniq_atoms;
            relevant_atoms[uniq_atoms++] = atnoi;
	    printf("Internal safety check corrected for atom %d - "
		       "please report.", (int)atnoi);
	}
    }
    if(uniq_atoms<=1) { /* 0 cannot happen and 1 - no partitioning. */
        free(relevant_atoms);
        free(map2r);
	free(mu);
	free(mu2);
	free(g_mu);
        return point_cnt;
    }

    /*allocate memory and initialize*/
    block_work_init(&ggw, coorw, uniq_atoms, point_cnt);

    /*compute distance between grid points and relevant atoms*/
    block_compute_rjs(&ggw, mg, uniq_atoms, relevant_atoms,
                      point_cnt);

    /*calculate the partition functions*/
    for(i=1; i<uniq_atoms; i++) {
	    integer atno = relevant_atoms[i];
	    integer posi = ggw.LDA*i;
	    for(j=0; j<i; j++) {
		    real dist, bfac;
		    integer atno2 = relevant_atoms[j];
		    integer posat = atno2+(atno*(atno-1))/2;
		    dist = mg->rij[posat];
		    bfac = mg->aij[posat];
		    integer posj = ggw.LDA*j;
		    for(ptno=0; ptno<point_cnt; ptno++) {
			    integer posipt = posi+ptno;
			    integer posjpt = posj+ptno;

			    real mu_1 = (ggw.rj[posipt]-ggw.rj[posjpt])*dist;
			    mu_1 += bfac*(1-mu_1*mu_1);
			    mu[ptno] = mu_1;
			    mu2[ptno] = mu_1*mu_1;
		    }

                   for(h=0;h<hardness;h++){
                        for(ptno=0; ptno<point_cnt; ptno++) {
                          mu[ptno] = 0.5*mu[ptno]*(3.0-mu2[ptno]);
                        }
                        if(h<hardness-1){ 
                          for(ptno=0; ptno<point_cnt; ptno++) mu2[ptno] = mu[ptno]*mu[ptno];
                        }
                    }
                    for(ptno=0; ptno<point_cnt; ptno++) {
                            ggw.p_kg[posi+ptno] *= 0.5-0.5*mu[ptno];
                            ggw.p_kg[posj+ptno] *= 0.5+0.5*mu[ptno];
                    }
	    }
    }


    /* compute weight normalization factors */
    pc = point_cnt; 
    /* ls_dzero_(ggw.vec, &pc);*/
    for(ptno=0; ptno<point_cnt; ptno++) {
      ggw.vec[ptno] = 0.0000000000000000; 
    }
    for(i=0; i<uniq_atoms; i++) {
	    integer posi = ggw.LDA*i;
	    for(ptno=0; ptno<point_cnt; ptno++)
		    ggw.vec[ptno] += ggw.p_kg[posi+ptno];
    }

    /*
     * Apply the computed weights removing at the same time
     * points with low weight.
     */

    for(dest=0, ptno=0; ptno<point_cnt; ptno++) {
            integer atom = map2r[atom_nums[ptno]];
            real factor = ggw.p_kg[ggw.LDA*atom+ptno]/ggw.vec[ptno];
            if(factor >= WEIGHT_THRESHOLD){
               coorw[dest][0] = coorw[ptno][0];
               coorw[dest][1] = coorw[ptno][1];
               coorw[dest][2] = coorw[ptno][2];
               coorw[dest][3] = coorw[ptno][3]*factor;
               dest++;
            }
    }

    free(relevant_atoms);
    free(map2r);

    free(mu);
    free(mu2);
    free(g_mu);

    block_work_release(&ggw);
    return dest;
}

/* Structure describing chosen partitioning scheme */
struct partitioning_scheme_t {
	char *name;
	void (*preprocess)
		(GridGenMolGrid* mg, integer point_cnt,
		 GridGenWork *ggw, integer idx, integer verbose, integer skip);
	integer (*postprocess)
		(GridGenMolGrid* mg, const real *center, 
		 integer point_cnt, const integer *atom_nums,
		 real (*coor)[4], const integer *nprocess);
};
static struct partitioning_scheme_t partitioning_becke_corr =
{ "Becke partitioning with atomic radius correction", 
	becke_corr_preprocess, NULL };
static struct partitioning_scheme_t partitioning_becke_orig =
{ "Original Becke partitioning", becke_orig_preprocess };
static struct partitioning_scheme_t partitioning_ssf =
{ "SSF partitioning", ssf_preprocess, NULL };
static struct partitioning_scheme_t partitioning_block =
{ "Blocked partitioning for large molecules/parallel calc:s. based on Beckes scheme with atomic radius correction",
	NULL, block_postprocess };
static struct partitioning_scheme_t partitioning_block_ssf =
{ "Blocked SSF partitioning",
        NULL, block_postprocess_ssf };


/* =================================================================== */
/* linear grid generator                                               */
/* =================================================================== */
/*the bucketing scheme. first shot: nlog n. Load up all grid points,
  assign them [nx,ny,nz] tuples.  sort the tuples to collect the ones
  sharing same cube. Within each cube select one that is closest to
  the center. Save cubes.
 */
#if defined(CHAR_BIT)
#undef CHAR_BIT
#endif
#define CHAR_BIT 8
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

/* lo and hi point to three-element vectors */
	static void
boxify_create_index(real cell_size, real (*coor)[4], integer point_cnt,
		struct point_key_t * index,
		real *lo, real *hi)
{
	static const char coordinate_name[] = "xyz";
	real fac;
	integer i, bpc = KEY_BITS/3; /* bits per coordinate */
	lo[0] = lo[1] = lo[2] = +1e20;
	hi[0] = hi[1] = hi[2] = -1e20;

	for(i=0; i<point_cnt; i++) {
		if(coor[i][0] < lo[0])      lo[0] = coor[i][0];
		else if(coor[i][0] > hi[0]) hi[0] = coor[i][0];
		if(coor[i][1] < lo[1])      lo[1] = coor[i][1];
		else if(coor[i][1] > hi[1]) hi[1] = coor[i][1];
		if(coor[i][2] < lo[2])      lo[2] = coor[i][2];
		else if(coor[i][2] > hi[2]) hi[2] = coor[i][2];
	}

	for(i=0; i<3; i++) {
		if( (1<<bpc) < (hi[i] - lo[i])/cell_size)
			fprintf(stderr, "error: %c dimension [%g %g] cell size %f\n",
					coordinate_name[i], hi[i],lo[i], cell_size);
	}

	/* shift lo so that the grid cells stick equally on the all sides */
	for(i=0; i<3; i++) {
		real w = hi[i]-lo[i];
		integer cnt = ceil(w/CELL_SIZE);
		real rw = cnt*CELL_SIZE;
		lo[i] -= (rw-w)*0.5;
	}

	fac = 1/cell_size;
	for(i=0; i<point_cnt; i++) {
		integer ix = (coor[i][0] - lo[0])*fac;
		integer iy = (coor[i][1] - lo[1])*fac;
		integer iz = (coor[i][2] - lo[2])*fac;
		if(ix<0) ix = 0; /* correct numerical error, if any */
		if(iy<0) iy = 0; /* correct numerical error, if any */
		if(iz<0) iz = 0; /* correct numerical error, if any */
		index[i].key = (ix<<(bpc*2)) | (iy<<bpc) | iz;
		index[i].index = i;
	}

	qsort(index, point_cnt, sizeof(struct point_key_t), comp_point_key);
}

/** OrbData is a struct containing information about active shells at
  given chunk of grid points. One has to store additional
  information when PBC are enabled. */
struct OrbData {
	real *rshell;
	integer (*shlblocks)[2];
	real cell_size;
	integer kmax;
	integer nblocks;
	/* PBC specific block */
	integer max_layers;           /**< 0 means noperiodic calculation. */
	real (*shell_centers)[3]; /**< used only in PBC code */
	integer *pbc_blocks;      /**< replaces nblocks which is scalar */
	integer used_blocks;          /**< number of used blocks in shlblocks
					and pbc_blocks */
	integer *pbc_idx; /**< ERIK!? */
	struct GridPbcData *pbc;
};

	static integer
orbdata_write(struct OrbData *od, FILE *f)
{
  if(fwrite(&od->max_layers,  sizeof(od->max_layers), 1, f) != 1)
    return 1;
  if(fwrite(&od->nblocks,  sizeof(integer), 1, f) != 1)
    return 6;
  if(fwrite(od->shlblocks, sizeof(integer), od->nblocks*2, f)!= od->nblocks*2)
    return 7;
  return 0;
}

	static void
ggen_work_init(GridGenWork* ggw, GridGenMolGrid* mg, integer mxshells)
{
	ggw->rj   = dal_new(mg->atom_cnt*MAX_PT_PER_SHELL, real);
	ggw->p_kg = dal_new(mg->atom_cnt*MAX_PT_PER_SHELL, real);
	ggw->vec  = dal_new(mg->atom_cnt*MAX_PT_PER_SHELL, real);
	ggw->x    = dal_new(MAX_PT_PER_SHELL*mxshells, real);
	ggw->y    = dal_new(MAX_PT_PER_SHELL*mxshells, real);
	ggw->z    = dal_new(MAX_PT_PER_SHELL*mxshells, real);
	ggw->wg   = dal_new(MAX_PT_PER_SHELL*mxshells, real);
	ggw->atomv= dal_new(MAX_PT_PER_SHELL*mxshells, integer);
        ggw->relevant_atoms = dal_new(mg->atom_cnt,integer);
        ggw->process = dal_new(MAX_PT_PER_SHELL*mxshells, integer);
}

	static void
ggen_work_release(GridGenWork* ggw)
{
	free(ggw->rj);
	free(ggw->p_kg);
	free(ggw->vec);
	free(ggw->x);
	free(ggw->y);
	free(ggw->z);
	free(ggw->wg); 
	free(ggw->atomv);
        free(ggw->relevant_atoms);
        free(ggw->process);
}

	static GridGenMolGrid*
mgrid_new(DftGrid *conf, integer atom_cnt, const GridGenAtom* atoms)
{
	integer i, j, index, paircnt;
	GridGenMolGrid* mg = dal_new(1, GridGenMolGrid);
	mg->conf        = conf;
	mg->pbc         = conf->pbc_enabled ? &conf->pbc : NULL;
	mg->atom_cnt    = atom_cnt;
	mg->atom_coords = atoms;
	mg->atom_grids  = dal_new(atom_cnt, GridGenAtomGrid*);

	/* careful / ErikT */
	/*mg->total_points = 0;
	  mg->off = -17;
	  mg->nt = -17; */
	/* careful */

	/* be careful with atoms (=no pairs): some systems (AIX) do not 
	 * like 0 allocations */
	paircnt = (atom_cnt*(atom_cnt-1))/2;

	/* careful / ErikT */
	/* printf("mgrid_new(): atom_cnt = %d, paircnt = %d\n",
	   atom_cnt,paircnt); */
	/* careful / ErikT */
	mg->rij  = dal_new((paircnt == 0 ? 1 : paircnt), real);
        mg->distij  = dal_new((paircnt == 0 ? 1 : paircnt), real);
	mg->aij  = dal_new((paircnt == 0 ? 1 : paircnt), real);
	index =0;
	for(i=0; i<atom_cnt; i++) {
		mg->atom_grids[i] = agrid_new(atoms[i].icent, atoms[i].Z);
		for(j=0; j<i; j++) {
			real chi = bragg_radii[atoms[i].Z]/bragg_radii[atoms[j].Z];
			real temp = (chi-1)/(chi+1);
			real dx = atoms[i].x - atoms[j].x;
			real dy = atoms[i].y - atoms[j].y;
			real dz = atoms[i].z - atoms[j].z;
			temp = temp/(temp*temp-1);
			if(temp>0.5) temp = 0.5;
			else if(temp<-0.5) temp = -0.5;
                        real dist = sqrt(dx*dx + dy*dy + dz*dz);
                        mg->distij[index] = dist;
			mg->rij[index] = 1.0/dist;
			mg->aij[index++] = temp;
		}
	}

        /*distance to the nearest neighbor, needed for SSF scheme*/
        mg->nn = dal_new(atom_cnt,real);
        if(conf->partitioning->name == "SSF partitioning" ||
           conf->partitioning->name == "Blocked SSF partitioning"){
           for(i=0; i<atom_cnt; i++) {
              real dist=1.e10;
              for(j=0; j<atom_cnt; j++) {
                 if(i!=j){
                 if(j<i) {index=j+(i*(i-1))/2;}
                 if(j>i) {index=i+(j*(j-1))/2;}
                 if(dist >= mg->distij[index]) dist=mg->distij[index];
                 }
              }
              mg->nn[i]=0.5*(1-SSF_CUTOFF)*dist;
           }
        }
	mg->verbose = 0;

#if 0
	/*compute the tabulated cutoff profile for use in the evaluation of the partition function
     *------------------------------------------------------------------------------------------------
     * the cutoff profile is defined as 
     *
     *   s = (1-p_h)/2
     *
     * with p being a recursively defined polynomial with the recursion order given by the hardness h
     *   p_h = 1.5*p_{h-1}-0.5*p_{h-1}^3
     *   p_0 = 1.5*mu-0.5mu^3 
     *   
     * the first derivative of s wrt mu can be evaluated recursively as
     *   s1 = -0.5\prod_{h=1}^{hardness}{1.5(1-p_{h-1}^2)}
     *
     * the second derivate reads as
     *   s2 = \sum_{h=1}^{hardness}{s1*\frac{2p_{h-1}*\frac{d p_{h-1}}{d mu}}{1-p_{h-1}^2}}
     *
     * /Andreas Krapp
     */
    real fac, fac2, g_mu, s, s1, s2,mu,mu2,g_mu1,g_mu2;
    int h,index2;
    static const real SSF_CUTOFF = 0.64  ; /*must be smaller than 1.0 and positiv*/
    static const integer finefac = 10000;  
    real conv = 1.0/finefac;
    int ipoints = SSF_CUTOFF*finefac*2;

    mg->SSF_CUTOFF = SSF_CUTOFF;
    mg->finefac = finefac;
    mg->finefacinv = conv;
    mg->tabcutoff0 = dal_new(ipoints, real);
    mg->tabcutoff1 = dal_new(ipoints, real);
    mg->tabcutoff2 = dal_new(ipoints, real);

    for(i = 0; i<SSF_CUTOFF*finefac;i++){
#if 1
/*SSF definition*/
          mu = i*conv;
          mu2 = mu*mu;
          g_mu  = 0.0625*mu*(35+mu2*(-35+mu2*(21-5*mu2)));
          g_mu1 = 2.1875*(1-mu2*(3-mu2*(3-mu2)));
          g_mu2 = 13.125*mu*(-1+mu2*(1-mu2));
	  /*       index = SSF_CUTOFF*finefac-1+i;*/
          index = i;
          mg->tabcutoff0[index] = 0.5*(1-g_mu);
          mg->tabcutoff1[index] = -g_mu1;
          mg->tabcutoff2[index] = -g_mu2/3.0;
/*
          if(i>0){ 
             index2 = SSF_CUTOFF*finefac-1-i;
             mg->tabcutoff0[index2] = 0.5*(1+g_mu);
             mg->tabcutoff1[index2] = g_mu1;
             mg->tabcutoff2[index2] = g_mu2/3.0;
          }
          printf("index %i index2 %i counter i %i mu %e val1 %e\n",index,index2,i,mu,mg->tabcutoff0[index]);
*/
#else
       g_mu = i/finefac;
       fac2 = 1.0;
       s1 = 1.0;
       s2 = 0.0;
       for(h=0; h<conf->hardness; h++){          
          fac  = 1.0-g_mu*g_mu;        /*fac=(1-p_{h-1}^2)*/
          s2  += g_mu*fac2/fac;        /*s2_incr = p_{h-1}*1.5*(1-p_{h-2}^2)/(1-p_{h-1}^2*/
          fac2 = 1.5*fac;               
          g_mu = 0.5*g_mu*(2.0+fac);
          s1  *= fac2;
       }
       mg->tabcutoff0[i] = 0.5*(1.0-g_mu);
       mg->tabcutoff1[i] = -0.5*s1/2.0;
       mg->tabcutoff2[i] = s1*s2/6.0;
#endif
   }
   
   /*-----------------------------------------------------------------------------------------------*/
#endif

    /*some constants for the SSF scheme*/
    real cutoff_inv = 1.0/SSF_CUTOFF;
    real cutoffinv2 = cutoff_inv*cutoff_inv;
    real fac0 = 0.5*0.0625*cutoff_inv;
    real fac1 = 35*fac0;
    real fac2 = -35*cutoffinv2*fac0;
    real fac3 =  21*cutoffinv2*cutoffinv2*fac0;
    real fac4 =  -5*cutoffinv2*cutoffinv2*cutoffinv2*fac0;
    mg->ssf_fac0 = fac0;
    mg->ssf_fac1 = fac1;
    mg->ssf_fac2 = fac2;
    mg->ssf_fac3 = fac3;
    mg->ssf_fac4 = fac4;

    return mg;
}

static void
mgrid_free(GridGenMolGrid* mg)
{
    integer i;
    free(mg->rij);
    free(mg->distij);
    free(mg->aij);
    for(i=0; i< mg->atom_cnt; i++)
        agrid_free(mg->atom_grids[i]);
    free(mg->atom_grids);
    free(mg->nn);
#if 0
    free(mg->tabcutoff0);
    free(mg->tabcutoff1);
    free(mg->tabcutoff2);
#endif
    free(mg);
}

/* set_ang_fixed:
   computes the angular points for a set of radial points
   obtained from radial integration scheme.
   The only thing altered is grid->leb_ang vector.
*/
static void
mgrid_set_angular_fixed(GridGenMolGrid* mgrid, real thrl,
                        integer minang, integer maxang)
{
    static const real BOHR = 0.529177249;
    integer atom, i = 0;

    for(atom=0; atom<mgrid->atom_cnt; atom++) {
	integer atom_maxang = maxang;
        GridGenAtomGrid* grid = mgrid->atom_grids[atom];
        real rbragg = bragg_radii[grid->Z]/(2.0*BOHR);

#if 0
        /* the original scheme, which is different from the turbomole scheme.
         * I commented this out, since we only use this when using the turbomole grids 
         * and then we should use the real turbomole scheme for consistency \Andreas Krapp*/

	if(mgrid->conf->Z_dependent_maxang) {
	  if(grid->Z<=2)
	     atom_maxang -= 10;
	  else if(grid->Z<=10)
	     atom_maxang -= 3;
          else if(grid->Z<=18)
	     atom_maxang += 3;
	  else if(grid->Z<=36)
	     atom_maxang += 10;
	  else atom_maxang += 15;
	}
#else
        /* this is the maxang scheme from turbomole*/
        if(mgrid->conf->Z_dependent_maxang) {
          if(grid->Z<=2){
             atom_maxang -= 6;
             if(maxang == 47) atom_maxang -= 6; }
        }
#endif
        for (i=0; i<grid->pnt; i++) {
	  integer iang = leb_get_from_order(atom_maxang);
#if 1
           /*the original pruning, follows Murray, Handy, Laming (MolPhys.1993,78,997)*/
          if(mgrid->conf->prune && grid->rad[i]<rbragg) {
                integer n_points_optimal = 
                  (real)leb_gen[iang].point_cnt*grid->rad[i]/rbragg;
                integer iang1 = leb_get_from_point(n_points_optimal);
                iang = iang1 <iang ? iang1 : iang;
          }
#else
           /*the turbomole pruning as described in Ahlrichs-Treutler paper, but this is not what is done in 
             the turbomole programm, therefore we do not use it as standard */
          if(mgrid->conf->prune) {
               real intpart;
               integer a,b;
               if (modf(grid->pnt/2.0, &intpart) != 0.0) { 
                  a = (integer)(grid->pnt/2.0);}
               else {a = (integer)(grid->pnt/2.0)-1; }
               if (modf(grid->pnt*2.0/3.0, &intpart) != 0.0) {
                  b = (integer)(2*grid->pnt/3.0);}
               else { b = (integer)(2*grid->pnt/3.0)-1;}

               if(i > a) {
                  if (i > b) {
                     iang = leb_get_from_order(5); }
                  else{
                     iang = leb_get_from_order(11);}
               }
          }
#endif
          grid->leb_ang[i] = iang;
        }
    }
}

#if defined(USE_PTHREADS)
#include <pthread.h>
static pthread_mutex_t grid_mutex = PTHREAD_MUTEX_INITIALIZER;
static void
mgrid_lock(GridGenMolGrid* mgrid)
{
    pthread_mutex_lock(&grid_mutex);
}
static void
mgrid_unlock(GridGenMolGrid* mgrid)
{
    pthread_mutex_unlock(&grid_mutex);
}
static integer
mgrid_get_num_threads(void)
{ return 4; }
#else
#define mgrid_lock(mgrid)
#define mgrid_unlock(mgrid)
#endif

static integer
compress_grid(integer point_cnt, real* x, real* y, real *z, real* wg, integer* atomv, integer* process)
{
    integer i, dest=0;
    /* 1e-5 factor below is really a fudge factor... */
    real local_threshold = WEIGHT_THRESHOLD*1e-5;
    for(i=0; i<point_cnt; i++) {
        x [dest] = x [i];
        y [dest] = y [i];
        z [dest] = z [i];
        wg[dest] = wg[i];
        atomv[dest] = atomv[i];
        process[dest] = process[i];
        if(fabs(wg[i]) > local_threshold)
            dest++;
    }

    return dest;
}

static void*
mgrid_compute_coords_worker(GridGenMolGrid* mgrid, integer printlu)
{
#if defined(VAR_OMP)
  integer off = omp_get_thread_num();
#else
  integer off = mgrid->off;
#endif
  
  integer atom, ptno, j, idx, cnt, tpt, done;
  GridGenWork ggw;
  
    mgrid_unlock(mgrid); /* off has been copied, can unlock */

    j=0;
    for (atom=off; atom < mgrid->atom_cnt; atom += mgrid->nt) {
        if (mgrid->atom_grids[atom]->pnt > j) 
            j = mgrid->atom_grids[atom]->pnt;
    }
    
    if (j == 0)
	return NULL; /* this thread has got nothing to do */

    ggen_work_init(&ggw, mgrid, j);

    for(atom = off, done = 0; atom < mgrid->atom_cnt; 
        atom += mgrid->atom_coords[atom].mult, done++) {
        GridGenAtomGrid* grid = mgrid->atom_grids[atom];
        integer mult = mgrid->atom_coords[atom].mult;
        if(done % mgrid->nt != 0) /* do every mt symmetry-independent atom */
            continue; 
        idx = 0;
        integer skip=0;
        ggw.firstpoint=0;
        for (ptno=0; ptno<grid->pnt; ptno++) {
            real fact = 4*M_PI*grid->wght[ptno]*mult;
            real rad = grid->rad[ptno];
            integer ind = grid->leb_ang[ptno]; 
            leb_gen[ind].func(ggw.x+idx, ggw.y+idx, ggw.z+idx, ggw.wg+idx,
                              &tpt);

            /* for SSF partioning: 
             * if the points are in the sphere defined by the SSF condition, 
             * the partition function is 1.0, 
             * and no processing of this point necessary */
	    if(mgrid->conf->partitioning->name == "SSF partitioning" || 
               mgrid->conf->partitioning->name == "Blocked SSF partitioning"){
               if(grid->rad[ptno]<mgrid->nn[atom]) skip=1;}

            for(j=idx; j<idx+leb_gen[ind].point_cnt; j++) {
                ggw.x[j] = ggw.x[j]*rad + mgrid->atom_coords[atom].x;
                ggw.y[j] = ggw.y[j]*rad + mgrid->atom_coords[atom].y;
                ggw.z[j] = ggw.z[j]*rad + mgrid->atom_coords[atom].z;
                ggw.wg[j] *= fact;
		ggw.atomv[j] = atom;
                ggw.process[j] = skip;
            }

            if(mgrid->conf->partitioning->preprocess) {
                mgrid->conf->partitioning->preprocess
                    (mgrid, leb_gen[ind].point_cnt, &ggw, idx, 0, skip);
	    }
            idx += leb_gen[ind].point_cnt;
        }

	/* degeneracy multiplication here */
	cnt = compress_grid(idx, ggw.x, ggw.y, ggw.z, ggw.wg, ggw.atomv, ggw.process);

        if(mgrid->verbose){
#if defined(VAR_OMP)
#pragma omp critical 
#endif
	  if(printlu !=0 ){
            lsfort_print(printlu,"Atom: %4d*%d points=%5d compressed from %5d (%3d radial)", 
                       atom+1, mult, cnt, idx, grid->pnt);
	  }
	  if(printlu ==0 ){
            printf("Atom: %4d*%d points=%5d compressed from %5d (%3d radial)", 
                       (int)atom+1, (int)mult, (int)cnt, (int)idx, (int)grid->pnt);
	  }}
        if(cnt>0) {
            integer i;
            mgrid_lock(mgrid);
#if defined(VAR_OMP)
#pragma omp critical 
#endif
      {
            fwrite(&cnt,  sizeof(integer), 1,  mgrid->fl);
            fwrite(ggw.atomv, sizeof(integer), cnt,  mgrid->fl); 
            fwrite(ggw.process, sizeof(integer), cnt,  mgrid->fl); 
            for(i=0; i<cnt; i++) {
                fwrite(ggw.x+i, sizeof(real), 1, mgrid->fl);
                fwrite(ggw.y+i, sizeof(real), 1, mgrid->fl);
                fwrite(ggw.z+i, sizeof(real), 1, mgrid->fl);
                fwrite(ggw.wg+i,sizeof(real), 1, mgrid->fl);
            }
            mgrid->total_points += cnt;
      }
            mgrid_unlock(mgrid);
        }
    }
    ggen_work_release(&ggw);
    return NULL;
}

static integer
mgrid_compute_coords(GridGenMolGrid* mgrid, const char* filename, integer printlu)
{
    if( (mgrid->fl = fopen(filename,"wb")) == NULL) {
        fprintf(stderr,"ERROR: Cannot open grid file '%s' for writing.\n",
                filename);
        return 0;
    }
    mgrid->total_points = 0;
#if defined(USE_PTHREADS)
    mgrid->nt = mol_grid_get_num_threads();
    { integer i;
    pthread_t *ptid = dal_new(mgrid->nt, pthread_t);
    for(i=0; i<mgrid->nt; i++) {
        mol_grid_lock(mgrid);
        mgrid->off = i;
        pthread_create(&ptid[i], NULL, 
                       (void *(*)(void *))mgrid_compute_coords_worker, mgrid);
    }
    for(i=0; i<mgrid->nt; i++) 
        pthread_join(ptid[i], NULL);
    }
#elif defined(VAR_OMP)
    mgrid->nt = omp_get_max_threads();
#pragma omp parallel shared(mgrid,printlu)
    {
      mgrid_compute_coords_worker(mgrid,printlu);    
    }
#else
    mgrid->nt = 1; mgrid->off = 0;
    mgrid_compute_coords_worker(mgrid,printlu);
#endif
    fclose(mgrid->fl);

    return mgrid->total_points;
}

/* ------------------------------------------------------------------- */
/* the dft grid input routines.                                        */
/* choose different types of grid:                                     */
/* The syntax is:                                                      */
/* ((GC2|LMG|TURBO) | (SSF|BECKE|BECKEORIG|BLOCK|BLOCKSSF))            */
/* Sets                                                                */
/* radial_quad                                                         */
/* selected_partitioning                                               */
/* ------------------------------------------------------------------- */
static integer
get_word(const char *line, integer st, integer max_len)
{
    integer res;
    for(res=st; res<max_len && isalnum(line[res]); res++)
        ;
    if(res>max_len) res = 0;
    return res;
}

static integer
skip_spaces(const char *line, integer st, integer max_len)
{
    integer res;
    for(res=st+1; res<max_len && isspace(line[res]); res++)
        ;
    return res;
}

/* ------------------------------------------------------------------- */
/* the integrator/grid generator itself                                */
/* ------------------------------------------------------------------- */
 /* TheGrid - a singleton object. Because of the hardcoded grid file
    name we can have only one Grid object. */
static DftGrid TheGrid = {
  &quad_lmg, 
  &partitioning_becke_orig,
  0, /* Z-dependent maxang, per default off */ 
  1  /* prune, per default on*/
};


void
grid_gen_set_part_scheme(GridGenPartScheme scheme)
{
    switch(scheme) {
    case GRID_PART_BECKE_CORR: TheGrid.partitioning = &partitioning_becke_corr;
        break;
    case GRID_PART_BECKE_ORIG: TheGrid.partitioning = &partitioning_becke_orig;
        break;
    case GRID_PART_SSF:        TheGrid.partitioning = &partitioning_ssf;
        break;
    case GRID_PART_BLOCK:      TheGrid.partitioning = &partitioning_block;
        break;
    case GRID_PART_BLOCK_SSF:  TheGrid.partitioning = &partitioning_block_ssf;
        break;
    }
}

/* grid_gen_set_rad_quad:
   set radial quadrature.
void
grid_gen_set_rad_quad(GridGenQuad scheme)
{
    switch(scheme) {
    default:
    case GRID_RAD_GC2: TheGrid.radial_quad = &quad_gc2; break;
    case GRID_RAD_LMG: TheGrid.radial_quad = &quad_lmg; break;
    case GRID_RAD_TUR: TheGrid.radial_quad = &quad_turbo; break;
    }
}
*/


void
FSYM(dftgridinput)(const char *line, integer *TURBO, integer line_len)
{
    static const char* keywords[] = {
        "GC2","LMG", "TURBO", "SSF","BECKE","BECKEORIG","BLOCK","BLOCKSSF","CARTESIAN"
    };
    integer st, en, i;
    integer start = isspace(line[0]) 
	? skip_spaces(line, 0, line_len) : 0;
    for(st=start; (en=get_word(line, st, line_len)) != 0 && en>st; 
        st=skip_spaces(line, en, line_len)) {
        for(i=0; i<ELEMENTS(keywords); i++) {
            if(strncasecmp(keywords[i], line+st, en-st)==0)
                break;
        }
	/* printf("line %4d: %s matched i=%d\n", st, line+st, i); */
        switch(i) {
        case 0: TheGrid.radial_quad = &quad_gc2; break;
        case 1: TheGrid.radial_quad = &quad_lmg; break;
        case 2:
	  TheGrid.radial_quad = &quad_turbo;
	  TheGrid.Z_dependent_maxang = 1;
     *TURBO = 1;
	  break;
        case 3: TheGrid.partitioning = &partitioning_ssf;break;
        case 4: TheGrid.partitioning = &partitioning_becke_corr; break;
        case 5: TheGrid.partitioning = &partitioning_becke_orig; break;
        case 6: TheGrid.partitioning = &partitioning_block; break;
        case 7: TheGrid.partitioning = &partitioning_block_ssf; break;
        case 8: 
            gridType = GRID_TYPE_CARTESIAN; break;
        default: printf("GRIDGEN: Unknown .GRID TYPE option ignored.\n%s",
                            line);
            /* FIXME: should I quit here? */
        }
    }
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

static integer
read_coords(FILE *f, integer maxlen, real (*coor)[3], real *weight)
{
    integer sz = -1;
    if(fread(&sz, sizeof(integer), 1, f) <1) {
        return -1; /* end of file */
    }
    if(sz>maxlen) {
        fprintf(stderr,
                "grid_getchunk: too long vector length in file: %d > %d\n"
                "Calculation will stop.\n", (int)sz, (int)maxlen);
        dalton_quit("grid_getchunk: too long vector length in file: %d > %d\n"
                "Calculation will stop.\n", (int)sz, (int)maxlen);
        return -1; /* stop this! */
    }

    if(fread(coor,   sizeof(real), 3*sz, f) < sz) { puts("XYZ");return -1;}
    if(fread(weight, sizeof(real), sz, f) < sz) { puts("W");return -1;}
    return sz;
}

/** grid_getchunk_blocked() reads grid data from the grid file also
    with screening information if only nblocks and shlblocks are
    provided. The data is read to @param coor array that will contain
    the grid point coordinates and to @param weight that will contain
    the associated weights. These arrays must be preallocated and have
    length @param maxlen. Now, the information about the basis
    function shells relevant for this chunk of grid points is returned
    via @param nblocks and @param shlblocks arguments. We do not
    return the active shell numbers each other separately. They can be
    packed instead into blocks of consecutive active shells. Argument
    nblock will be set to the number of these blocks and shlblocks
    entries will be filled with the beginnings (shlblocks[*][0]) and
    ends (shlblocks[*][1]) of the blocks of active shells.
 */
integer
grid_getchunk_blocked(DftGridReader* rawgrid, integer maxlen,
                      integer *nblocks, integer (*shlblocks)[2],
                      real (*coor)[3], real *weight)
{
    integer layers = 0, sz = 0, rc, bl_cnt;
    FILE *f = rawgrid->f;

    if( (sz = read_coords(f, maxlen, coor, weight)) <0)
        return sz;

    if(fread(&layers, sizeof(integer), 1, f) <1)
        return -1; /* end of file */

    if(layers != 0)
        dalton_quit("grid generated for PBC but read from ordinary code");
    /* if(fread(&bl_cnt, sizeof(unsigned), 1, f) <1) { stinne */
    if(fread(&bl_cnt, sizeof(integer), 1, f) <1) {
        puts("OCNT reading error."); return -1;
    }
    if(nblocks) *nblocks = bl_cnt;
    if(shlblocks) {
        rc = fread(shlblocks, sizeof(integer), bl_cnt*2, f);
    } else {
        integer cnt;
	integer buf;
        rc = 0;
        for(cnt=0; cnt<bl_cnt*2; cnt+=rc) {
            rc=fread(&buf, sizeof(integer), 1, f);
            if( rc < 1)
                break;
        }
    }
    if(rc<1) {
        fprintf(stderr,
                "IBLOCKS reading error: failed to read %d blocks.\n",
                (int)*nblocks);
        return -1; 
    }

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
FSYM(clsqua)(void)
{
    grid_close(grid);
    grid = NULL;
}

void
FSYM(reaqua)(integer *nshell, integer (*shell_blocks)[2], integer *buf_sz, 
	real (*coor)[3], real *weight, integer *nlen)
{
    *nlen = grid_getchunk_blocked(grid, *buf_sz, nshell, shell_blocks,
                                  coor, weight);
}

void FSYM(ii_getblocks)(const real *center, real *CELLSZ, real RSHEL2[],
			integer *KMAX, const integer *NCENT, integer *NBAST, integer *NATOM,
			const real *X, const real *Y,const real *Z,
			integer *NBLCNT, integer (*IBLCKS)[2]);

static void
II_boxify_load_grid(const char *fname, integer point_cnt,real (**coorw)[4],
		 integer *x_allocated, integer *atom_idx, integer *process)
{
    integer idx, cnt;
    real *chunk;
    FILE *f;

    *x_allocated = 1;
    chunk = malloc(4*point_cnt*sizeof(real));
    if(!chunk) {
	/*	fprintf(stderr, "loading grid into mem failed, too.\n"
		"worksz=%d needed size=%u\n", worksz,
		(unsigned) (4*point_cnt*sizeof(real)) );*/
	dalton_quit("no mem in load_grid: 2");
    }
    *coorw = (real(*)[4]) chunk;
     
    f = fopen(fname,"rb");
    if(!f) {
        fprintf(stderr, "internal error, cannot open the grid file %s.\n",
                fname);
        exit(1);
    }
    idx = 0;
    while(fread(&cnt, sizeof(cnt), 1, f)==1) {
        integer aidx, i;
        assert(cnt+idx<=point_cnt);
        if(atom_idx) {
	    if(fread(atom_idx+idx, sizeof(integer), cnt, f) != cnt)
		dalton_quit("ERROR GRID1");
	} else {
            for(i=0; i<cnt; i++)
		if(fread(&aidx, sizeof(integer), 1, f) != 1)
		    dalton_quit("ERROR GRID2");
	}

        if(fread(process+idx, sizeof(integer), cnt, f) != cnt)
            dalton_quit("ERROR GRID3");

        if(fread(*coorw+idx, sizeof(real), 4*cnt, f) != 4*cnt)
            dalton_quit("ERROR GRID1, cnt=%d", cnt);
        idx += cnt;
    }
    fclose(f);
}

static void
II_orbdata_destroy(struct OrbData *od)
{
  /*    free(od->rshell);*/
    free(od->shlblocks);
    if(od->shell_centers) free(od->shell_centers);
}

	static char*
grid_get_fname(const char *base, integer filenum)
{
	if(filenum == 0)
		return strdup(base);
	else {
		char *res = malloc(strlen(base) + 15);
		sprintf(res, "%s.%05d", base, (int)filenum);
		return res;
	}
}

static void
II_orbdata_init(struct OrbData *od, struct GridPbcData *pbc,
             real cell_size, integer shell_cnt, real *RSHEL)
{
    od->rshell = RSHEL;
    od->kmax = shell_cnt;
    od->cell_size = cell_size;
    od->shlblocks = dal_malloc(2*shell_cnt*sizeof(integer));
    od->pbc = NULL;
    od->max_layers = 0;
    od->shell_centers = NULL;
    od->pbc_idx       = NULL;
    od->pbc_blocks    = NULL;
}

static integer
II_orbdata_compute_active_orbs(struct OrbData *od, const real center[],
			       integer KMAX, const integer *NCENT, integer NBAST, 
			       integer NATOM, const real *X, const real *Y, 
			       const real *Z)
{
    if(od->max_layers) {
      /*  pbc code */
	dalton_quit("not implemented for pbc");
    } else {
	integer nbl;
        FSYM(ii_getblocks)(center, &od->cell_size, od->rshell,&KMAX,
			   NCENT,&NBAST,&NATOM,X,Y,Z,&nbl, od->shlblocks);
	od->nblocks = nbl;
    }
    return od->nblocks != 0;
}

static void
II_write_final_coords_and_weights(integer cnt, integer *it, integer *maxnlen, real *coorw, FILE *f)
{
    integer i;
    if(cnt <= 0) return;
    *it = *it + 1;
    if(cnt>*maxnlen){
      *maxnlen=cnt;
    }
/* printf("II_write_final_coords_and_weights, it = %d\n", *it); */
    if(fwrite(&cnt, sizeof(cnt), 1, f)!=1) {
        fprintf(stderr, "GRIDGEN: 'too short write' error.\n");
        dalton_quit("GRIDGEN: 'too short write' error.\n");
    }
    /* qsort(coorw, cnt, 4*sizeof(real), comp_weight); */
    for(i=0; i<cnt; i++)
        if(fwrite(coorw+i*4, sizeof(real), 3, f) != 3)
           dalton_quit("write error in %s(), coords", __FUNCTION__);
    for(i=0; i<cnt; i++)
        if(fwrite(coorw+i*4+3, sizeof(real), 1, f) != 1)
           dalton_quit("write error in %s(), weights", __FUNCTION__);
}

static integer
II_boxify_save_batch_local(GridGenMolGrid *mg, FILE *f, integer cnt,
                        struct OrbData *od, const real *center,
                        const integer *atom_nums, real *coorw, integer KMAX,
			const integer *NCENT, integer NBAST, integer NATOM, const real *X,
			   const real *Y, const real *Z, integer *it, integer *maxnlen, const integer *nprocess)
{

    integer rc;
    if (!II_orbdata_compute_active_orbs(od, center,KMAX,NCENT,NBAST,NATOM,X,Y,Z)) {
        return cnt; /* No orbitals active at given grid point - unlikely
                   but possible. */
    }
    if(mg->conf->partitioning->postprocess) {
        integer icnt=0;
        icnt = mg->conf->partitioning->postprocess(mg, center, cnt, atom_nums,
						  (real(*)[4])coorw,nprocess);
         cnt=icnt;
    }
    if(cnt == 0) return cnt;
    II_write_final_coords_and_weights(cnt, it, maxnlen, coorw, f);

    if ( (rc=orbdata_write(od, f)) != 0) {
	dalton_quit("write error in orbdata_write(), point %d", rc);
    }
    return cnt;
}

/* Simen should be if defined VAR_MPI, but II_boxify_save_batch_local did not work directly*/
#ifdef VAR_MPI
static integer
II_boxify_save_batch(GridGenMolGrid *mg, FILE *f, integer cnt,
                  struct OrbData *od, real *center,
                  integer *atom_nums, real *coorw, integer KMAX, const integer *NCENT,
		  integer NBAST, integer NATOM, const real *X, const real *Y,
		     const real *Z, integer *it, integer *maxnlen,integer *nprocess)
{
    const integer PARBLLEN = 1000;
    integer i;
    if (!II_orbdata_compute_active_orbs(od, center,KMAX,NCENT,NBAST,NATOM,X,Y,Z)) {
        return cnt; /* No orbitals active at given grid point - unlikely
                   but possible. */
    }

    integer icnt=0;
    for(i=0; i<cnt; i+= PARBLLEN) {
        integer bcnt = i+PARBLLEN<cnt ? PARBLLEN : cnt - i;
	icnt+=II_boxify_save_batch_local(mg, f, bcnt, od,
				   center, atom_nums+i, coorw + i*4,KMAX,
					 NCENT,NBAST,NATOM,X,Y,Z,it,maxnlen,nprocess+i);
    }
    return icnt;
}
#else
#define II_boxify_save_batch(mg,f,c,od,center,an,r,km,nc,             \
                             nb,na,x,y,z,it,maxnlen,np)		      \
        II_boxify_save_batch_local((mg),(f),(c),(od),(center),(an),   \
				   (r),(km),(nc),(nb),(na),(x),(y),(z),(it),(maxnlen),(np))
#endif


/** II_boxify_save() saves the grid, possibly distributing it over to
 * many nodes. For each chunk of grid points, it generates the list of
 * relevant basis function shells that is saved along with the
 * points. NOTE: atom_idx is used only for partitioning schemes that
 * do postprocessing and should not be accessed without prior
 * checking.
 **/
static integer
II_boxify_save(GridGenMolGrid *mg, const char *fname, integer point_cnt,
	       real (*coorw)[4], const integer *atom_idx,
	       struct point_key_t *keys, real *lo, real cell_size,
	       const real *X, const real *Y, const real *Z,
	       real *RSHEL,integer NBAST, integer KMAX,integer NATOM,
	       const integer *NCENT, integer *it, integer *maxnlen, integer lupri, const integer *process)
{
    FILE *f;
    integer idx, cnt;
    integer shell_cnt = KMAX;
    integer bpc = KEY_BITS/3; /* bits per coordinate */
    GridPointKey mask = ~((-1)<<bpc);
    integer points_saved = 0;
    real *dt = NULL;
    integer dt_sz = 0, *atom_nums = NULL, *nprocess = NULL;
    struct OrbData od;

    II_orbdata_init(&od, mg->pbc, cell_size, shell_cnt, RSHEL);

    if((f = fopen(fname,"wb")) == NULL) {
        fprintf(stderr,"internal error, cannot save sorted grid file.\n");
        exit(1);
    }
    for(idx=0; idx<point_cnt; idx += cnt) {
        integer i;
        real center[3];
        real mindist = 4*cell_size, maxdist = 0;
        GridPointKey key = keys[idx].key;
        real sx = 0, sy = 0, sz = 0;
        center[0] = lo[0] + (((key >> (bpc*2)) & mask)+0.5)*cell_size;
        center[1] = lo[1] + (((key >> bpc)     & mask)+0.5)*cell_size;
        center[2] = lo[2] + (((key)            & mask)+0.5)*cell_size;
        for(cnt=0; idx+cnt<point_cnt &&
                keys[idx+cnt].key == key; cnt++) {
        }

        /* merge data in one block */
        if(dt_sz<cnt) {
            dt_sz = cnt;
            dt = realloc(dt, 4*dt_sz*sizeof(real));
            atom_nums = realloc(atom_nums, dt_sz*sizeof(integer));
            nprocess = realloc(nprocess,dt_sz*sizeof(integer));
        }
        for(i=0; i<cnt; i++) {
            integer j = keys[idx+i].index;
            dt[i*4+0]    = coorw[j][0];
            dt[i*4+1]    = coorw[j][1];
            dt[i*4+2]    = coorw[j][2];
            dt[i*4+3]    = coorw[j][3];
            if(atom_idx) atom_nums[i] = atom_idx[j];
            nprocess[i] = process[j];
        }

        if(0) {
        printf("box [%5.1f,%5.1f,%5.1f] nt=%5d "
               "%f <r<%f [%5.1f,%5.1f,%5.1f]\n",
               center[0], center[1], center[2], (int)cnt,
               sqrt(mindist), sqrt(maxdist),
               sx/cnt, sy/cnt, sz/cnt);
        puts("");
        }

	if (II_orbdata_compute_active_orbs(&od,center,KMAX,NCENT,
					   NBAST,NATOM,X,Y,Z)){
            integer icnt=0;
	    icnt=II_boxify_save_batch(mg, f, cnt, &od, center, atom_nums, dt,
				      KMAX,NCENT,NBAST,NATOM,X,Y,Z,it,maxnlen,nprocess);
            points_saved += icnt; 
        }
    }
    fclose(f);
    II_orbdata_destroy(&od);
    if(dt) free(dt);
    if(atom_nums) free(atom_nums);
    if(nprocess) free(nprocess);
    if(point_cnt != points_saved)
      lsfort_print(lupri,"Postprocessing compression: from %d to %d",
                    point_cnt, points_saved);
    return points_saved;
}

/* mgrid_set_radial:
   precompute radial grid for all the atoms in the molecule.
*/
static void
II_mgrid_set_radial(GridGenMolGrid* grd, struct radial_scheme_t *radial_quad,
		    real thrl, const integer *NCENT,const integer *NHKT,
		    const integer *NUCO, integer NHTYP,	integer NUCIND, integer KMAX,
		    integer MXPRIM, const real *PRIEXP, const integer *NSTART, integer TURBO, 
		    integer GC2)
{

    integer atom;
    void *data = NULL;
    if (radial_quad->quad_init) {
      if(TURBO){
      data = NULL;
      }
      else{
      data = gen_lmg_init2(NCENT,NHKT,NUCO,
                         NHTYP,NUCIND, KMAX,MXPRIM,PRIEXP,NSTART);}
    }

    for (atom = 0; atom < grd->atom_cnt; atom++) {
      if(TURBO){
        turbo_quad_generate(grd->atom_grids[atom], thrl, data);
      }
      else if(GC2){
        gen_gc2_quad(grd->atom_grids[atom], thrl, data);
      }
      else{
        gen_lmg_quad2(grd->atom_grids[atom], thrl, data, NCENT,NHKT,
			      NUCO,NHTYP,NUCIND,KMAX,MXPRIM,PRIEXP,NSTART);
      }
    }
    if (radial_quad->quad_free)
        radial_quad->quad_free(data);
}

static integer
II_mgrid_boxify(GridGenMolGrid *mg, const char* fname,
		integer point_cnt, real cell_size,
		const real *X, const real *Y, const real *Z,
		real *RSHEL,integer NBAST, integer KMAX,integer NATOM,
		const integer *NCENT, integer *it,integer *maxnlen, integer lupri)
{
    real (*coorw)[4];
    real lo[3], hi[3];
    struct point_key_t * keys =
      dal_new(point_cnt, struct point_key_t);
    integer coorw_allocated = 0;
    integer *atom_idx;
    integer *process=dal_new(point_cnt, integer);
    integer new_point_cnt;

    atom_idx = mg->conf->partitioning->postprocess ?
        dal_new(point_cnt, integer) : NULL;

    II_boxify_load_grid(fname, point_cnt,
                     &coorw, &coorw_allocated, atom_idx,process);
    boxify_create_index(cell_size, coorw, point_cnt, keys, lo, hi);
    new_point_cnt = II_boxify_save(mg, fname, point_cnt, coorw, atom_idx,
                                keys, lo, cell_size,X,Y,Z,RSHEL,NBAST,KMAX,
				   NATOM,NCENT,it,maxnlen,lupri,process);
    free(keys);
    if(coorw_allocated) free(coorw);
    if(atom_idx) free(atom_idx);
    if(process)  free(process);
    return new_point_cnt;
}

/* grid_generate:
   returns number of grid points.
*/

integer
II_grid_generate(DftGrid *grid, const char* filename, integer atom_cnt,
		 const GridGenAtom* atom_arr,
		 GridGeneratingFunc generating_function,
		 const integer *NCENT, const integer *NHKT, const integer *NUCO,
		 integer NHTYP, integer KMAX, integer MXPRIM, const real *PRIEXP,
		 const integer *NSTART, const real *X, const real *Y,
		 const real *Z, real *RSHEL,integer nbast, integer *it, integer *maxnlen,integer TURBO,integer GC2, integer lupri)
{
    integer res;
    struct tms starttm, endtm; clock_t utm;
    GridGenMolGrid* mgrid =  mgrid_new(grid, atom_cnt, atom_arr);
    lsfort_print(lupri,"Radial Quadrature : %s", grid->radial_quad->name);
    lsfort_print(lupri,"Space partitioning: %s", grid->partitioning->name);
    lsfort_print(lupri,"Radial integration threshold: %g", grid->radint);
    lsfort_print(lupri,"Angular integration order range: [%i %i] Angular points %u",
	       grid->angmin, grid->angmax,
	       leb_gen[leb_get_from_order(grid->angmax)].point_cnt);
    if (grid->partitioning->name != "SSF partitioning" && grid->partitioning->name != "Blocked SSF partitioning")
        lsfort_print(lupri,"Hardness of the partioning function: %i", grid->hardness);
    if (grid->prune){lsfort_print(lupri,"Grid pruning: on");}
    else {lsfort_print(lupri,"Grid pruning: off");};   

    times(&starttm);
    mgrid->verbose = 1;
    II_mgrid_set_radial(mgrid, grid->radial_quad, grid->radint,NCENT,
			NHKT,NUCO,NHTYP,atom_cnt,KMAX,MXPRIM,PRIEXP,NSTART,TURBO,GC2);

    mgrid_set_angular_fixed(mgrid, grid->radint, grid->angmin,
                            grid->angmax);
    /*    mgrid_set_angular_fixed(mgrid,grid->radint,leb_get_from_order(grid->angmin)
			    ,leb_get_from_order(grid->angmax));
    */
    res = mgrid_compute_coords(mgrid, filename,lupri);

/* stinne: */
/* printf("II_grid_generate: res = %d\n", res);
printf("(int) II_grid_generate: res = %d\n", (int) res);
printf("sizeof = %d\n", sizeof(res)); */
    res = II_mgrid_boxify(mgrid, filename, res, CELL_SIZE,
			  X,Y,Z,RSHEL,nbast,KMAX,atom_cnt,NCENT,it,maxnlen,lupri);
    mgrid_free(mgrid);

    times(&endtm);
    utm = endtm.tms_utime-starttm.tms_utime;

    lsfort_print(lupri,"Number of grid points: %8d Grid generation time: %9.1f s",
               res, utm/(real)sysconf(_SC_CLK_TCK));

    return res;
}

DftGridReader* II_grid_get_reader(DftGrid *grid, integer nbast, integer natoms,
				  const real *X,
				  const real *Y, const real *Z,
				  const integer *Charge,integer *grdone,
				  const integer *NCENT, const integer *NHKT,
				  const integer *NUCO, integer NHTYP, integer KMAX,
				  integer MXPRIM, const real *PRIEXP,
				  const integer *NSTART, real *RSHEL,integer *it,integer *maxnlen, integer TURBO, integer GC2, integer lupri)
{
  DftGridReader *res = dal_new(1, DftGridReader);
  char *fname;
  
  switch(gridType) {
  case GRID_TYPE_STANDARD:
    if(grid->regenerate) {
      integer atom_cnt;
      integer iqq;
      atom_cnt = natoms;
      GridGenAtom* atoms;
      atoms = dal_new(natoms, GridGenAtom);
      for(iqq=0; iqq<natoms; iqq++){
	atoms[iqq].icent = iqq;
	atoms[iqq].Z = Charge[iqq];
	  atoms[iqq].mult = 1;
	  atoms[iqq].x = X[iqq];
	  atoms[iqq].y = Y[iqq];
	  atoms[iqq].z = Z[iqq];
      }
      
      II_grid_generate(grid, "DALTON.QUAD2", atom_cnt, atoms,
		       NULL, NCENT,NHKT,NUCO,NHTYP,KMAX,
		       MXPRIM,PRIEXP,NSTART,X,Y,Z,RSHEL,nbast,it,maxnlen,TURBO,GC2,lupri);
      
      free(atoms);
      *grdone = 1;
      grid->regenerate = 0;
    }
    fname = grid_get_fname("DALTON.QUAD2", 0);
    res->f=fopen(fname, "rb");
    free(fname);
    if(res == NULL) {
      perror("DFT quadrature grid file DALTON.QUAD2 not found.");
      free(res);
      abort();
    }
    return res;
    
  case GRID_TYPE_CARTESIAN:
    perror("Error: cartesian grid not implemented for NEW II ROUTINE.\n");
    free(res);
    abort();
  default:
    perror("Error in II_grid_open: unknown grid type\n");
    free(res);
    abort();
  }
  /* END SWITCH */
  return NULL; /* to keep some compilers quiet */
}

DftGridReader*
II_grid_open(integer nbast, real radint,
	     integer angmin,integer angint, integer hardness, integer prune, integer natoms,
	     const real *X, const real *Y, const real *Z, const integer *Charge,
	     integer *grdone, const integer *NCENT, const integer *NHKT, const integer *NUCO,
	     integer NHTYP, integer KMAX,integer MXPRIM, const real *PRIEXP,
	     const integer *NSTART, real *RSHEL, integer *it, integer *maxnlen,integer TURBO, integer GC2, integer lupri)  
{

    TheGrid.pbc_enabled = 0;
    TheGrid.radint = radint;
    TheGrid.angmin = angmin;
    TheGrid.angmax = angint;
    TheGrid.hardness = hardness;
    TheGrid.prune = prune;
    if(!*grdone) TheGrid.regenerate = !*grdone;

    return II_grid_get_reader(&TheGrid, nbast, natoms,
			      X,Y,Z,Charge,grdone,NCENT,NHKT,NUCO,NHTYP,KMAX,
			      MXPRIM,PRIEXP,NSTART,RSHEL,it,maxnlen,TURBO,GC2,lupri);
}

void
FSYM(ii_opnqua)(integer *nbast,
		real *radint, integer *angmin, integer *angint,
		integer *hardness, integer *prune, integer *natoms,
		const real *X, const real *Y, const real *Z,
		const integer *Charge, integer *grdone,
		const integer *NCENT, const integer *NHKT, const integer *NUCO,
		integer *NHTYP, integer *KMAX, integer *MXPRIM,
		const real *PRIEXP, const integer *NSTART,real *RSHEL,
		integer *it, integer *maxnlen,
		integer *TURBO, integer *GC2, integer *LUPRI)
{

    grid = II_grid_open(*nbast, *radint,
			*angmin,*angint,*hardness,*prune,*natoms, X,Y,Z, Charge, grdone,
			NCENT, NHKT, NUCO, *NHTYP, *KMAX, *MXPRIM, PRIEXP,
			NSTART, RSHEL,it,maxnlen,*TURBO,*GC2,*LUPRI);
}











