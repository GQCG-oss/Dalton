/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* The linearly-scaling quadratic response XC evaluator.
   (c) Pawel Salek, pawsa@theochem.kth.se.
   2004.05.26 - 2005.03
*/

#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>

#define __CVERSION__
#include "general.h"
#include "integrator.h"
#include "functionals.h"

#include "inforb.h"

const static int DFTQR_DEBUG = 0;

void FSYM(deq27)(const real* cmo, const real* ubo, const real* dv, 
                 real* dxcao, real* dxvao, real* wrk, integer* lfrsav);

typedef struct {
    real *res_omega, *res_omY, *res_omZ;
    real *my, *mz, *myzzy; /* matrices we compute expectation values of *
                            * my is [rho,Y], etc. */
    real *yy, *zz, *yzzy;  /* vectors of expectation values */
    real *vy, *vz;         /* transformed orbitals */
    integer spinY, symY; 
    integer spinZ, symZ;
    integer symYZ; /* symmetry of the result = symA */
    unsigned is_gga:1; /* whether we should do the density-gradient
                        * dependent part of the transformation */
} QuadBlData;

static const real MONER = -1.0;
static const real HALFR =  0.5;
/* [D, X] = 
 *     [ A   B ]             [ 0 B ]
 * X = [ C   D ] -> [D, X] = [-C 0 ] 
 * where the blocks sizes are determined by the occupied/virtual split.
 */
static void
commute_d_x(real * RESTRICT kappaY, integer symY,
            real * RESTRICT commuted_mat)
{
    int isym, i, j;
    integer norbt = inforb_.norbt;
#if 0
    int nocc = inforb_.nocct;
    /* occupied columns */
    for(i=0; i<nocc; i++) {
        for(j=0; j<nocc; j++)     commuted_mat[j+i*norbt] = 0;
        for(j=nocc; j<norbt; j++) commuted_mat[j+i*norbt] = -kappaY[j+i*norbt];
    }
    /* virtual columns */
    for(i=nocc; i<norbt; i++) {
        for(j=0; j<nocc; j++)     commuted_mat[j+i*norbt] = kappaY[j+i*norbt];
        for(j=nocc; j<norbt; j++) commuted_mat[j+i*norbt] = 0;
    }
#else
    memset(commuted_mat, 0, norbt*norbt*sizeof(real));
    for(isym=0; isym<inforb_.nsym; isym++) {
        /* the block in question corresponds to (isym,jsym) block */
        int istart= inforb_.iorb[isym];
        int iocc  = inforb_.nocc[isym];
        int iorb  = inforb_.norb[isym];
        int jsym  = inforb_.muld2h[symY-1][isym]-1;
        int jstart= inforb_.iorb[jsym];
        int jocc  = inforb_.nocc[jsym];
        int jorb  = inforb_.norb[jsym];
        real * RESTRICT res = commuted_mat + istart + jstart*norbt;
        real * RESTRICT ky  = kappaY       + istart + jstart*norbt;
        /* occupied columns */
        for(j=0; j<jocc; j++) {
            for(i=0;    i<iocc; i++) res[i+j*norbt] = 0;
            for(i=iocc; i<iorb; i++) res[i+j*norbt] = -ky[i+j*norbt];
        }
        /* virtual columns */
        for(j=jocc; j<jorb; j++) {
            for(i=0;    i<iocc; i++) res[i+j*norbt] = ky[i+j*norbt];
            for(i=iocc; i<iorb; i++) res[i+j*norbt] = 0;
        }
    }
#endif
}

static __inline__ void
commute_matrices(integer sz, const real* a, const real* b,
                 real* c, int addp)
{
    static const real MONER = -1.0;
    const real* firstpref = addp ? &ONER : &ZEROR;
    /* this could be optimized to use symmetry... */
    dgemm_("N", "N", &sz, &sz, &sz, &ONER, 
	   a, &sz, b, &sz, firstpref, c, &sz);
    dgemm_("N", "N", &sz, &sz, &sz, &MONER, 
           b, &sz, a, &sz, &ONER, c, &sz);
}

int FSYM(isetksymop)(const integer *new_ksymop);
static void
qrbl_data_init(QuadBlData *d, real *cmo, integer is_gga, int max_block_len,
	      real *kappaY, integer symY, integer spinY, 
	      real *kappaZ, integer symZ, integer spinZ,
	      real *dmat,
	      real *work, integer *lwork)
{
    real dummy;
    integer isym;
    integer nbast = inforb_.nbast;
    integer norbt = inforb_.norbt;
    real *commuted_mat = work;
    real *work1 = commuted_mat + inforb_.n2orbx;
    real *work2 = work1        + inforb_.n2basx;
    int allocated = inforb_.n2orbx + 2*inforb_.n2basx;
    int comps = is_gga ?  4 : 1; /* number of components */

    d->res_omega = calloc(inforb_.n2basx,sizeof(real));
    d->res_omY   = calloc(inforb_.n2basx,sizeof(real));
    d->res_omZ   = calloc(inforb_.n2basx,sizeof(real));
    d->my    = dal_malloc(inforb_.n2basx*sizeof(real));
    d->mz    = dal_malloc(inforb_.n2basx*sizeof(real));
    d->myzzy = dal_malloc(inforb_.n2basx*sizeof(real));
    d->yy   = dal_malloc(comps*max_block_len*sizeof(real));
    d->zz   = dal_malloc(comps*max_block_len*sizeof(real));
    d->vy   = dal_malloc(comps*max_block_len*nbast*sizeof(real));
    d->vz   = dal_malloc(comps*max_block_len*nbast*sizeof(real));
    d->yzzy = dal_malloc(comps*max_block_len*sizeof(real));
    d->spinY = spinY; d->symY = symY;
    d->spinZ = spinZ; d->symZ = symZ;
    d->symYZ = inforb_.muld2h[symY-1][symZ-1];

    d->is_gga = is_gga != 0; /* make sure lowest bit is set */
    if(*lwork<allocated)
        dalton_quit("no enough mem in %s", __FUNCTION__);
    isym = isetksymop_(&symY);
    FSYM(deq27)(cmo, kappaY, &dummy, d->my,    &dummy, work, lwork);
    isetksymop_(&symZ);
    FSYM(deq27)(cmo, kappaZ, &dummy, d->mz,    &dummy, work, lwork);
    isetksymop_(&isym);
#if 0
    fort_print("cmo");
    outmat_(cmo,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
            &inforb_.nbast, &inforb_.nbast);
    fort_print("dmat");
    outmat_(dmat,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
            &inforb_.nbast, &inforb_.nbast);
    fort_print("input Y");
    outmat_(kappaY,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
            &inforb_.nbast, &inforb_.nbast);
    fort_print("inptu Z");
    outmat_(kappaZ,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
            &inforb_.nbast, &inforb_.nbast);

    fort_print("mY");
    outmat_(d->my,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
            &inforb_.nbast, &inforb_.nbast);
#endif
    /* yzzy = c*(dmo*(y*z+z*y) - (y*dmo*z + z*dmo*y))*c' */
    /* FIXME: There is also place for optimization since large blocks
     * of these the partial results have 0 blocks and probably entire
     * [Y,[Z,D]] could be done in one shot considerity sparsity of D
     * in MO basis. */
    commute_d_x(kappaY, symY, commuted_mat);
    commute_matrices(norbt, commuted_mat, kappaZ, work1, 0);
    commute_d_x(kappaZ, symZ, commuted_mat);
    commute_matrices(norbt, commuted_mat, kappaY, work1, 1);
    if(DFTQR_DEBUG) {
        fort_print("mYZZY (MO), sym: %d", d->symYZ);
        outmat_(work1,&ONEI,&inforb_.norbt,&ONEI,&inforb_.norbt,
                &inforb_.norbt, &inforb_.norbt);
    }

    /* transform work1 to AO basis, (isym,jsym) blocks... */
    memset(d->myzzy, 0, nbast*nbast*sizeof(real));
    for(isym=0; isym<inforb_.nsym; isym++) {
        int ibasi = inforb_.ibas[isym];
        integer nbasi = inforb_.nbas[isym];
        int iorbi = inforb_.iorb[isym];
        integer norbi = inforb_.norb[isym];
        int icmoi = inforb_.icmo[isym];
        int jsym  = inforb_.muld2h[d->symYZ-1][isym]-1;
        int ibasj = inforb_.ibas[jsym];
        integer nbasj = inforb_.nbas[jsym];
        int iorbj = inforb_.iorb[jsym];
        integer norbj = inforb_.norb[jsym];
        int icmoj = inforb_.icmo[jsym];
        if(norbi == 0 || norbj == 0) continue;
        dgemm_("N","N", &nbasi, &norbj, &norbi,
               &HALFR, cmo + icmoi, &nbasi, work1 + iorbi + iorbj*norbt, &norbt,
               &ZEROR, work2, &nbasi);
        dgemm_("N","T", &nbasi, &nbasj, &norbj,
               &ONER, work2, &nbasi, cmo + icmoj, &nbasj,
               &ZEROR, d->myzzy + ibasi + ibasj*nbast, &nbast);
    }
    if(DFTQR_DEBUG) {
        fort_print("mYZZY (AO), sym: %d", d->symYZ);
        outmat_(d->myzzy,&ONEI,&inforb_.nbast,&ONEI,&inforb_.nbast,
                &inforb_.nbast, &inforb_.nbast);
    }
}

static void
qrbl_data_free(QuadBlData *d)
{
    free(d->res_omega); free(d->res_omY); free(d->res_omZ);
    free(d->my); free(d->mz); free(d->myzzy);
    free(d->yy); free(d->zz); free(d->yzzy);
    free(d->vy); free(d->vz);
}

static void
qrbl_eval_rho_vars(DftIntegratorBl* grid, QuadBlData *data,
		  real *tmp, integer bllen, int blstart, int blend)
{
    /* compute vector of transformed densities vt */
    FSYM2(getexp_blocked_lda)(&data->symY, data->my, grid->atv,
                        grid->bas_bl_cnt, grid->basblocks, &grid->shl_bl_cnt,
                        tmp, &bllen, data->yy);
    FSYM2(getexp_blocked_lda)(&data->symZ, data->mz, grid->atv,
                        grid->bas_bl_cnt, grid->basblocks, &grid->shl_bl_cnt,
                        tmp, &bllen, data->zz);
    FSYM2(getexp_blocked_lda)(&data->symYZ, data->myzzy, grid->atv,
                        grid->bas_bl_cnt, grid->basblocks, &grid->shl_bl_cnt,
                        tmp, &bllen, data->yzzy);
}


/* ===================================================================
 *            LDA-optimized code 
 * =================================================================== */
/* add_dft_contribution:
   computes dft contribution to the quadratic response using precomputed
   traces of the operators (trY, trZ, etc). 
*/
static __inline__ real
min(real a, real b)
{ return a>b ? b : a; }
static void
qrbl_add_lda_contribution(DftIntegratorBl* grid, QuadBlData* d,
                          real *tmp, int bllen, int blstart, int blend)
{
    static const real sgn[2]={1.0,-1.0};
    int i, j, k, isym, ibl, jbl;
#if defined(VAR_PGF77) || defined(SYS_DEC)
    real pref2b[DFT_BLLEN], pref2c[DFT_BLLEN], pref3[DFT_BLLEN];
#else
    real pref2b[bllen], pref2c[bllen], pref3[bllen];
#endif
    real * RESTRICT aos = grid->atv;
    int sY = d->spinY, sZ = d->spinZ;
    real *om = d->res_omega, *omY = d->res_omY, *omZ = d->res_omZ;
    FunDensProp dp = {0};

    for(j=blstart; j<blend; j++) {
        ThirdDrv    drvs; /* the functional derivatives */
        real weight = grid->weight[grid->curr_point+j];
        dp.rhoa = dp.rhob = 0.5*grid->r.rho[j];
        dftpot2_(&drvs, weight, &dp, 0, sY != sZ);
        /* FIXME: simplify this */
        /* and first order contribution */
        pref2b[j] = d->yy[j]*drvs.fRR[sY];
        pref2c[j] = d->zz[j]*drvs.fRR[sZ];
    
        /* third order, and part of second order contribution */
        pref3[j] =
            d->yy[j]*d->zz[j]*drvs.fRRR[sY|sZ] +
            d->yzzy[j]       *drvs.fRR[sY^sZ];
    }

    for(isym=0; isym<grid->nsym; isym++) {
        integer (*RESTRICT iblocks)[2] = BASBLOCK(grid,isym);
        int ibl_cnt = grid->bas_bl_cnt[isym];
        
        for(ibl=0; ibl<ibl_cnt; ibl++)
            for(i=iblocks[ibl][0]-1; i<iblocks[ibl][1]; i++) { 
                int ioff = i*bllen;
                for(k=blstart; k<blend; k++) {
                    d->vz[k+ioff] = -pref2b[k]*aos[k+ioff];
                    d->vy[k+ioff] = -pref2c[k]*aos[k+ioff];
                    tmp[k+ioff]   =  pref3[k]*aos[k+ioff];
                }
            }

        /* Compute contributions to om, omY and omZ. We could have
         * special case when all these three matrices have same
         * symmetry but we do not bother since this code is most
         * likely not going to be used a lot and it is best to keep it
         * simple. */
        for(ibl=0; ibl<ibl_cnt; ibl++) {
            for(i=iblocks[ibl][0]-1; i<iblocks[ibl][1]; i++) { 
                int jsym, jbl_cnt, ioff = i*inforb_.nbast;
                real * RESTRICT tmpi = tmp + i*bllen;
                real * RESTRICT vyi = d->vy + i*bllen;
                real * RESTRICT vzi = d->vz + i*bllen;
                integer (*RESTRICT jblocks)[2];

                jsym = inforb_.muld2h[d->symYZ-1][isym]-1;
                jblocks = BASBLOCK(grid,jsym);
                jbl_cnt = grid->bas_bl_cnt[jsym];
                for(jbl=0; jbl<jbl_cnt; jbl++) {
                    int jtop = min(jblocks[jbl][1],i+1);
                    for(j=jblocks[jbl][0]-1; j<jblocks[jbl][1]; j++) { 
                        real * RESTRICT aosj = aos + j*bllen;
                        real s = 0;
                        for(k=blstart; k<blend; k++) s += aosj[k]*tmpi[k];
                        om[j+ioff] += s;
                    }
                }
                jsym = inforb_.muld2h[d->symZ-1][isym]-1;
                jblocks = BASBLOCK(grid,jsym);
                jbl_cnt = grid->bas_bl_cnt[jsym];
                for(jbl=0; jbl<jbl_cnt; jbl++) {
                    int jtop = min(jblocks[jbl][1],i+1);
                    for(j=jblocks[jbl][0]-1; j<jblocks[jbl][1]; j++) { 
                        real * RESTRICT aosj = aos + j*bllen;
                        real s = 0;
                        for(k=blstart; k<blend; k++) s += aosj[k]*vyi[k];
                        omY[j+ioff] += s;
                    }
                }
                jsym = inforb_.muld2h[d->symY-1][isym]-1;
                jblocks = BASBLOCK(grid,jsym);
                jbl_cnt = grid->bas_bl_cnt[jsym];
                for(jbl=0; jbl<jbl_cnt; jbl++) {
                    int jtop = min(jblocks[jbl][1],i+1);
                    for(j=jblocks[jbl][0]-1; j<jblocks[jbl][1]; j++) { 
                        real * RESTRICT aosj = aos + j*bllen;
                        real s = 0;
                        for(k=blstart; k<blend; k++) s += aosj[k]*vzi[k];
                        omZ[j+ioff] += s;
                    }
                }
                /* Start next column i... */
            }
        }
    }
}
static void
qrbl_lda_cb(DftIntegratorBl* grid, real * RESTRICT tmp,
	    int bllen, int blstart, int blend,
	    QuadBlData* data)
{
    qrbl_eval_rho_vars       (grid, data, tmp, bllen, blstart, blend);
    qrbl_add_lda_contribution(grid, data, tmp, bllen, blstart, blend);
}

/* ===================================================================
 *            GGA-capable code 
 * =================================================================== */
static void
qrbl_eval_gga_vars(DftIntegratorBl* grid, QuadBlData *data,
                    real *tmp, integer bllen, int blstart, int blend)
{
    /* compute vector of transformed densities and density gradients vt */
    FSYM2(getexp_blocked_gga)(&data->symY, data->my, grid->atv,
                        grid->bas_bl_cnt, grid->basblocks, &grid->shl_bl_cnt,
                        tmp, &bllen, (double(*)[4])data->yy);
    FSYM2(getexp_blocked_gga)(&data->symZ, data->mz, grid->atv,
                        grid->bas_bl_cnt, grid->basblocks, &grid->shl_bl_cnt,
                        tmp, &bllen, (double(*)[4])data->zz);
    FSYM2(getexp_blocked_gga)(&data->symYZ, data->myzzy, grid->atv,
                        grid->bas_bl_cnt, grid->basblocks, &grid->shl_bl_cnt,
                        tmp, &bllen, (double(*)[4])data->yzzy);
}

static void
qrbl_add_gga_contribution(DftIntegratorBl* grid, QuadBlData* d,
                          real *tmp, int bllen, int blstart, int blend)
{
    static const real sgn[2]={1.0,-1.0};
    int i, j, k, isym, ibl, jbl;
    /* pref3 is the prefactor of the double-commuted term, pref2a is
     * the prefactor of commuted with Z, and pref2b - with commuted
     * with Y. */
#if 1 || defined(VAR_PGF77)
    real pref2b[DFT_BLLEN][4], pref2c[DFT_BLLEN][4], pref3[DFT_BLLEN][4];
#else
    real pref2b[bllen][4], pref2c[bllen][4], pref3[bllen][4];
#endif
    real * RESTRICT aos = grid->atv;
    real * RESTRICT aox = grid->atv+bllen*inforb_.nbast;
    real * RESTRICT aoy = grid->atv+bllen*inforb_.nbast*2;
    real * RESTRICT aoz = grid->atv+bllen*inforb_.nbast*3;
    int sY = d->spinY, sZ = d->spinZ;
    real *om = d->res_omega, *omY = d->res_omY, *omZ = d->res_omZ;
    FunDensProp dp = {0};
    real (*trY)[4]  = (real (*)[4])d->yy;
    real (*trZ)[4]  = (real (*)[4])d->zz;
    real (*yzzy)[4] = (real (*)[4])d->yzzy;
    real (*grad)[3] = grid->g.grad;
    int can_collapse_loops = 
        (d->symY == d ->symZ) && (d->symY == d ->symYZ);

    for(k=blstart; k<blend; k++) {
        real trYsum, trZsum, trYZZYsum, trYtimesZ, a, b, c, d;
        ThirdDrv    drvs; /* the functional derivatives */
        real weight = grid->weight[grid->curr_point+k];
        real ngrad = sqrt(grad[k][0]*grad[k][0]+
                          grad[k][1]*grad[k][1]+
                          grad[k][2]*grad[k][2]);
        dp.rhoa = dp.rhob = 0.5*grid->r.rho[k];
        dp.grada = dp.gradb = 0.5*ngrad;
        dp.gradab = dp.grada*dp.gradb;
        dftpot2_(&drvs, weight, &dp, 1, sY != sZ);
        trYsum = 0.5*(trY[k][1] * grad[k][0] + trY[k][2] * grad[k][1]
            + trY[k][3] * grad[k][2]);
        trZsum = 0.5*(trZ[k][1] * grad[k][0] + trZ[k][2] * grad[k][1]
            + trZ[k][3] * grad[k][2]);
        trYZZYsum = 0.5*(yzzy[k][1] * grad[k][0] + yzzy[k][2] * grad[k][1]
                         + yzzy[k][3] * grad[k][2]);
        trYtimesZ = trY[k][1]*trZ[k][1] + trY[k][2]*trZ[k][2]
            + trY[k][3]*trZ[k][3];
        /* FIXME: simplify this */
        /* and first order contribution */

        /* third order, and part of second order contribution */
        pref3[k][0] = trY[k][0]*trZ[k][0]*drvs.fRRR[sY|sZ]+
            yzzy[k][0]*drvs.fRR[sY^sZ] +
            2*(drvs.fRRZ[sZ][sY]*trY[k][0]*trZsum +
               drvs.fRRZ[sY][sZ]*trZ[k][0]*trYsum) +
               4*trZsum*trYsum*drvs.fRZZ[sY^sZ][sY] +
               2*(trYtimesZ + trYZZYsum)*drvs.fRZ[sY^sZ];

        pref3[k][0] += drvs.fRRG[sY]*trY[k][0]*trZsum*(1+sgn[sZ])+
            drvs.fRRG[sZ]*trZ[k][0]*trYsum*(1+sgn[sY])+
            0.5*drvs.fRG[0]*(1+sgn[sY]*sgn[sZ])*trYZZYsum +
            0.5*drvs.fRG[0]*(sgn[sY]+sgn[sZ])*trYtimesZ;

        pref2b[k][0] = trY[k][0]*drvs.fRR[sY] + 2*trYsum*drvs.fRZ[sY]
            + trYsum*drvs.fRG[sY];

        pref2c[k][0] = trZ[k][0]*drvs.fRR[sZ] + 2*trZsum*drvs.fRZ[sZ]
            + trZsum*drvs.fRG[sZ];


        a = trY[k][0]*(2*drvs.fRZ[sY]+drvs.fRG[sY])+4*trYsum*drvs.fZZ[sY];
        b = 2*drvs.fZ + sgn[sY]*drvs.fG;
        b *= 2;

	pref2b[k][1] = a*grad[k][0] + b*trY[k][1];
	pref2b[k][2] = a*grad[k][1] + b*trY[k][2];
	pref2b[k][3] = a*grad[k][2] + b*trY[k][3];


        a = trZ[k][0]*(2*drvs.fRZ[sZ]+drvs.fRG[sZ])+4*trZsum*drvs.fZZ[sZ];
        b = 2*drvs.fZ + sgn[sZ]*drvs.fG;
        b *= 2;

	pref2c[k][1] = a*grad[k][0] + b*trZ[k][1];
	pref2c[k][2] = a*grad[k][1] + b*trZ[k][2];
	pref2c[k][3] = a*grad[k][2] + b*trZ[k][3];

        /* xyz components of pref3 */
	a = (8*drvs.fZZZ[sY|sZ]*trYsum*trZsum +
             4*(drvs.fRZZ[sY][sZ]*trY[k][0]*trZsum + 
                drvs.fRZZ[sZ][sY]*trZ[k][0]*trYsum)+
             2*drvs.fRRZ[sY^sZ][sY]*trY[k][0]*trZ[k][0]);
        a += drvs.fRRGX[sY][sZ]*trY[k][0]*trZ[k][0]
            + 4*drvs.fZZ[sY^sZ]*trYZZYsum
            + 4*drvs.fZZ[sY^sZ]*trYtimesZ
            + 2*drvs.fRZ[sY^sZ]*yzzy[k][0];
        a += drvs.fRG[sY^sZ]*yzzy[k][0];

        b = 4*drvs.fZZ[sY]*trYsum + 2*drvs.fRZ[sY]*trY[k][0]
            + drvs.fRG[sY]*trY[k][0]*sgn[sZ];
        c = 4*drvs.fZZ[sZ]*trZsum + 2*drvs.fRZ[sZ]*trZ[k][0]
            + drvs.fRG[sZ]*trZ[k][0]*sgn[sY];
        d = 2*drvs.fZ + sgn[sY]*sgn[sZ]*drvs.fG; 
        b *= 2; c *= 2; d *= 2;

        pref3[k][1] = a*grad[k][0] + b*trZ[k][1] + c*trY[k][1] + d*yzzy[k][1];
        pref3[k][2] = a*grad[k][1] + b*trZ[k][2] + c*trY[k][2] + d*yzzy[k][2];
        pref3[k][3] = a*grad[k][2] + b*trZ[k][3] + c*trY[k][3] + d*yzzy[k][3];
    }

    for(isym=0; isym<grid->nsym; isym++) {
        integer (*RESTRICT iblocks)[2] = BASBLOCK(grid,isym);
        int ibl_cnt = grid->bas_bl_cnt[isym];
        for(ibl=0; ibl<ibl_cnt; ibl++)
            for(i=iblocks[ibl][0]-1; i<iblocks[ibl][1]; i++) { 
                real * RESTRICT a0 = aos + i*bllen;
                real * RESTRICT ax = aox + i*bllen;
                real * RESTRICT ay = aoy + i*bllen;
                real * RESTRICT az = aoz + i*bllen;
                real * RESTRICT vyi = d->vy + i*bllen;
                real * RESTRICT vzi = d->vz + i*bllen;
                real * RESTRICT tmi = tmp + i*bllen;
                for(k=blstart; k<blend; k++) {
                    vzi[k] = -(pref2b[k][0]*a0[k] + pref2b[k][1]*ax[k] +
                               pref2b[k][2]*ay[k] + pref2b[k][3]*az[k]); 
                    vyi[k] = -(pref2c[k][0]*a0[k] + pref2c[k][1]*ax[k] +
                               pref2c[k][2]*ay[k] + pref2c[k][3]*az[k]);
                    tmi[k] =  (pref3 [k][0]*a0[k] + pref3 [k][1]*ax[k] +
                               pref3 [k][2]*ay[k] + pref3 [k][3]*az[k]);
                }
            }
        /* Time to distribute computed prefactors. In principle, we
        should have three separate loops because three generated
        matrices may have different symmetries. On the other hand, it
        is profitable to detect the commonly occuring case of same
        symmetry and collapse the loops to reduce the overhead.
        1:special case when all three vectors have same symmetry */
        if(can_collapse_loops) {
            for(ibl=0; ibl<ibl_cnt; ibl++) {
                for(i=iblocks[ibl][0]-1; i<iblocks[ibl][1]; i++) { 
                    int ioff = i*inforb_.nbast;
                    real sum;
                    real * RESTRICT tmpi = tmp  + i*bllen;
                    real * RESTRICT vyi = d->vy + i*bllen;
                    real * RESTRICT vzi = d->vz + i*bllen;
                    int jsym = inforb_.muld2h[d->symYZ-1][isym]-1;
                    integer (*RESTRICT jblocks)[2] = BASBLOCK(grid,jsym);
                    int jbl_cnt = grid->bas_bl_cnt[jsym];
                    for(jbl=0; jbl<jbl_cnt; jbl++) {
                        for(j=jblocks[jbl][0]-1; j<jblocks[jbl][1]; j++) {
                            real * RESTRICT aosj = aos + j*bllen;
                            real oR = 0, oY = 0, oZ = 0;
                            for(k=blstart; k<blend; k++) {
                                oR += aosj[k]*tmpi[k];
                                oY += aosj[k]*vyi [k];
                                oZ += aosj[k]*vzi [k];
                            }
                            om [j+ioff] += oR;
                            omY[j+ioff] += oY;
                            omZ[j+ioff] += oZ;
                        }
                    }
                }
            }
            continue; /* no need to do the general case here */
        }

        /* 2:general case */
        for(ibl=0; ibl<ibl_cnt; ibl++) {
            for(i=iblocks[ibl][0]-1; i<iblocks[ibl][1]; i++) { 
                int ioff = i*inforb_.nbast;
                real sum;
                real * RESTRICT tmpi = tmp  + i*bllen;
                real * RESTRICT vyi = d->vy + i*bllen;
                real * RESTRICT vzi = d->vz + i*bllen;
                int jsym = inforb_.muld2h[d->symYZ-1][isym]-1;
                integer (*RESTRICT jblocks)[2] = BASBLOCK(grid,jsym);
                int jbl_cnt = grid->bas_bl_cnt[jsym];
                for(jbl=0; jbl<jbl_cnt; jbl++) {
                    for(j=jblocks[jbl][0]-1; j<jblocks[jbl][1]; j++) {
                        real * RESTRICT aosj = aos + j*bllen;
                        sum = 0;
                        for(k=blstart; k<blend; k++) sum += aosj[k]*tmpi[k];
                        om [j+ioff] += sum;
                    }
                }
                jsym = inforb_.muld2h[d->symZ-1][isym]-1;
                jblocks = BASBLOCK(grid,jsym);
                jbl_cnt = grid->bas_bl_cnt[jsym];
                for(jbl=0; jbl<jbl_cnt; jbl++) {
                    for(j=jblocks[jbl][0]-1; j<jblocks[jbl][1]; j++) {
                        real * RESTRICT aosj = aos + j*bllen;
                        sum = 0;
                        for(k=blstart; k<blend; k++) sum += aosj[k]*vyi[k];
                        omY[j+ioff] += sum;
                    }
                }
                jsym = inforb_.muld2h[d->symY-1][isym]-1;
                jblocks = BASBLOCK(grid,jsym);
                jbl_cnt = grid->bas_bl_cnt[jsym];
                for(jbl=0; jbl<jbl_cnt; jbl++) {
                    for(j=jblocks[jbl][0]-1; j<jblocks[jbl][1]; j++) {
                        real * RESTRICT aosj = aos + j*bllen;
                        sum = 0;
                        for(k=blstart; k<blend; k++) sum += aosj[k]*vzi[k];
                        omZ[j+ioff] += sum;
                    }
                }
            }
        }
    }
}

static void
qrbl_gga_cb(DftIntegratorBl* grid, real * RESTRICT tmp,
	    int bllen, int blstart, int blend,
	    QuadBlData* data)
{
    qrbl_eval_gga_vars       (grid, data, tmp, bllen, blstart, blend);
    qrbl_add_gga_contribution(grid, data, tmp, bllen, blstart, blend);
}

/* ===================================================================
 *                   Parallel section
 * =================================================================== */

#if defined(VAR_MPI)
#include <mpi.h>
#include <our_extra_mpi.h>
#define MASTER_NO 0
/* dft_qr_faster_slave:
   this is a slave driver. It's task is to allocate memory needed by
   the main property evaluator (dftqrcf_ in this case) and call it.
*/
void
dft_qrbl_slave(real* work, integer* lwork, const integer* iprint)
{
    real* fi    = malloc(inforb_.n2basx*sizeof(real));              /* OUT */
    real *cmo   = malloc(inforb_.norbt*inforb_.nbast*sizeof(real)); /* IN  */
    real *kappaY= malloc(inforb_.n2orbx*sizeof(real));              /* IN  */
    real *kappaZ= malloc(inforb_.n2orbx*sizeof(real));              /* IN  */
    integer addfock, symY, symZ, spinY, spinZ;                      /* IN  */
    FSYM2(dft_qr_respons)(fi, cmo, kappaY, &symY, &spinY, kappaZ, &symZ, &spinZ, 
			  &addfock, work, lwork);
    free(kappaZ);
    free(kappaY);
    free(cmo);
    free(fi);
}

static void
qrbl_sync_slaves(real* cmo, real* kappaY, real* kappaZ, integer* addfock,
		 integer* symY, integer* symZ, integer* spinY, integer* spinZ)
{
    static const SyncData sync_data[] = {
 	{ inforb_.nocc,   8, fortran_MPI_INT },
 	{ &inforb_.nocct, 1, fortran_MPI_INT },
 	{ &inforb_.nvirt, 1, fortran_MPI_INT },
    };
#ifdef C99_COMPILER
    const SyncData data2[] = {
	{ cmo,     inforb_.norbt*inforb_.nbast,MPI_DOUBLE },
	{ kappaY,  inforb_.n2orbx,             MPI_DOUBLE },
	{ kappaZ,  inforb_.n2orbx,             MPI_DOUBLE },
	{ addfock, 1,                          fortran_MPI_INT    },
	{ symY,    1,                          fortran_MPI_INT    },
	{ symZ,    1,                          fortran_MPI_INT    },
	{ spinY,   1,                          fortran_MPI_INT    },
	{ spinZ,   1,                          fortran_MPI_INT    }
    };
#else /* C99_COMPILER */
    /* this is more error-prone but some compilers (HP/UX)... */
    static SyncData data2[] = 
    { {NULL, 0, MPI_DOUBLE}, {NULL, 0, MPI_DOUBLE}, {NULL, 0, MPI_DOUBLE}, 
      {NULL, 0, fortran_MPI_INT   }, {NULL, 0, fortran_MPI_INT   }, {NULL, 0, fortran_MPI_INT   },
      {NULL, 0, fortran_MPI_INT   }, {NULL, 0, fortran_MPI_INT   } };
    data2[0].data = cmo;     data2[0].count = inforb_.norbt*inforb_.nbast;
    data2[1].data = kappaY;  data2[1].count = inforb_.n2orbx;
    data2[2].data = kappaZ;  data2[2].count = inforb_.n2orbx;
    data2[3].data = addfock; data2[3].count = 1;
    data2[4].data = symY;    data2[4].count = 1;
    data2[5].data = symZ;    data2[5].count = 1;
    data2[6].data = spinY;   data2[6].count = 1;
    data2[7].data = spinZ;   data2[7].count = 1;
#endif /* C99_COMPILER */

    mpi_sync_data(sync_data, ELEMENTS(sync_data));
    mpi_sync_data(data2,     ELEMENTS(data2));
}
static __inline__ void
qrbl_collect_info(real* fi, real*work, integer lwork)
{
    int sz = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &sz);
    if(sz<=1) return;

    CHECK_WRKMEM(inforb_.n2orbx,lwork);
    dcopy_(&inforb_.n2orbx, fi, &ONEI, work, &ONEI);
    int n2orbx_int = inforb_.n2orbx;
    MPI_Reduce(work, fi, n2orbx_int, MPI_DOUBLE, MPI_SUM, 
	       MASTER_NO, MPI_COMM_WORLD);
}
#else /* VAR_MPI */
#define qrbl_sync_slaves(cmo,kappaY,kappaZ,addfck,symY,symZ,spinY,spinZ)
#define qrbl_collect_info(fi,work,lwork)
#endif /* VAR_MPI */

/* ===================================================================
 *                   Main driver
 * =================================================================== */
void
FSYM2(dft_qr_respons)(real *fi, real *cmo,
                      real *kappaY, integer *symY, integer *spinY, 
                      real *kappaZ, integer *symZ, integer *spinZ,
                      integer *addfock, real *work, integer *lwork)
{
    static int msg_printed = 0;
    struct tms starttm, endtm; clock_t utm;
    real electrons;
    QuadBlData qr_data; /* quadratic response data */
    int i, j;
    real *dmat;
    real *tmp1, *tmp2;
    DftBlockCallback cb;

    /* WARNING: NO work MAY BE done before syncing slaves! */
    dft_wake_slaves((DFTPropEvalMaster)FSYM2(dft_qr_respons)); /* NO-OP in serial */
    qrbl_sync_slaves(cmo,kappaY,kappaZ,addfock,
		         symY,symZ, spinY,spinZ);
    if(*addfock) {/* call the old, slow version here */
#ifdef VAR_MPI
	fort_print("This was supposed not to be called.");
	dalton_quit("This calculation not reported to work in parallel.");
#endif
        dftqrcf_(fi, cmo, kappaY, symY, spinY, 
                 kappaZ, symZ, spinZ, addfock, 
                 work, lwork);
        return;
    }
    dmat = dal_malloc(inforb_.n2basx*sizeof(real));
    if(!msg_printed) {
        fort_print("DFT-QR computed in a linearly-scaling fashion.\n");
        msg_printed = 1;
    }
    times(&starttm);
    FSYM2(dft_get_ao_dens_mat)(cmo, dmat, work, lwork);
    qrbl_data_init(&qr_data, cmo, selected_func->is_gga(), DFT_BLLEN,
		   kappaY, *symY, *spinY, 
		   kappaZ, *symZ, *spinZ,
		   dmat, work, lwork);
    cb = (DftBlockCallback)
        (qr_data.is_gga ? qrbl_gga_cb : qrbl_lda_cb );
    electrons = dft_integrate_ao_bl(1, dmat, work, lwork, 0, 
                                    cb, &qr_data);
    free(dmat);
    if(DFTQR_DEBUG) {
        fort_print("DFT quadratic contribution (AO)");
        outmat_(qr_data.res_omega, &ONEI, &inforb_.nbast, &ONEI,
                &inforb_.nbast, &inforb_.nbast, &inforb_.nbast);
        outmat_(qr_data.res_omY, &ONEI, &inforb_.nbast, &ONEI,
                &inforb_.nbast, &inforb_.nbast, &inforb_.nbast);
        outmat_(qr_data.res_omZ, &ONEI, &inforb_.nbast, &ONEI,
                &inforb_.nbast, &inforb_.nbast, &inforb_.nbast);
    }
    /* LDA callback computes only half and GGA needs to be symmetrized */
    for(i=0; i<inforb_.nbast; i++) {
        for(j=0; j<i; j++) {
            int ji = j + i*inforb_.nbast;
            int ij = i + j*inforb_.nbast;
            real avg = 0.5*(qr_data.res_omega[ij]+qr_data.res_omega[ji]);
            qr_data.res_omega[ij] = qr_data.res_omega[ji] = avg;
            avg = 0.5*(qr_data.res_omY[ij] + qr_data.res_omY[ji]);
            qr_data.res_omY[ij] = qr_data.res_omY[ji] = avg;
            avg = 0.5*(qr_data.res_omZ[ij] + qr_data.res_omZ[ji]);
            qr_data.res_omZ[ij] = qr_data.res_omZ[ji] = avg;
        }
    }
    /* Transform to MO Fock matrix contribution: fix the symmetry of
     * the result!  */
    tmp1 = dal_malloc(inforb_.n2orbx*sizeof(real));
    tmp2 = dal_malloc(inforb_.n2orbx*sizeof(real));
    dzero_(tmp1, &inforb_.n2orbx);
    FSYM(lrao2mo)(cmo, &qr_data.symYZ, qr_data.res_omega, tmp1, work, lwork);
    /* compute commutators and add them */
    dzero_(tmp2, &inforb_.n2orbx);
    FSYM(lrao2mo)(cmo, &qr_data.symZ, qr_data.res_omY, tmp2, work, lwork);
    commute_matrices(inforb_.norbt, tmp2, kappaY, tmp1, 1);
    dzero_(tmp2, &inforb_.n2orbx);
    FSYM(lrao2mo)(cmo, &qr_data.symY, qr_data.res_omZ, tmp2, work, lwork);
    commute_matrices(inforb_.norbt, tmp2, kappaZ, tmp1, 1);
    qrbl_collect_info(tmp1, work, *lwork);  /* NO-OP in serial */
    daxpy_(&inforb_.n2orbx, &HALFR, tmp1, &ONEI, fi, &ONEI);
    if(DFTQR_DEBUG) {
        fort_print("DFT quadratic contribution (MO)");
        outmat_(fi, &ONEI, &inforb_.norbt, &ONEI, &inforb_.norbt,
                &inforb_.norbt, &inforb_.norbt);
    }
    free(tmp1); free(tmp2);
    qrbl_data_free(&qr_data);
    times(&endtm);
    utm = endtm.tms_utime-starttm.tms_utime;
    fort_print("      Electrons: %f(%9.3g): QR-DFT/b evaluation time: %9.1f s", 
               electrons, (double)(electrons-(int)(electrons+0.5)),
               utm/(double)sysconf(_SC_CLK_TCK));
}
