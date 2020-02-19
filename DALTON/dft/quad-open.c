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


void FSYM(deq27)(const real* cmo, const real* ubo, const real* dv,
                 real* dxcao, real* dxvao, real* wrk, integer* lfrsav);
integer FSYM(isetksymop)(const integer *new_ksymop);

typedef struct {
    real *dmata, *dmatb;       /* Denisty matrices                      */
    real *kappaY, *kappaZ;     /* Perturbation vectors                  */ 
    integer symY, symZ;            /* Symmetry of perturbations             */
    integer spinY, spinZ;          /* Spin rank of perturbations            */
    integer symYZ;                 /* Symmetry of mixed perturbation        */
    real *rhoa_y, *rhob_y;     /* First order perturbed density (AO)    */ 
    real *rhoa_z, *rhob_z;     /* First order perturbed density (AO)    */ 
    real *rhoa_yzzy;           /* Second order perturbed density (AO)   */
    real *rhob_yzzy;           /* Second order perturbed density (AO)   */ 
    real *rhowa_y, *rhowb_y;   /* First order perturbed density (grid)  */ 
    real *rhowa_z, *rhowb_z;   /* First order perturbed density (grid)  */ 
    real *rhowa_yzzy;          /* Second order perturbed density (grid) */
    real *rhowb_yzzy;          /* Second order perturbed density (grid) */ 
    real *res_omega_a;         /* Omega dependent part of qr contr.     */ 
    real *res_omega_b;         /* Omega dependent part of qr contr.     */   
    real *res_omY_a;           /* [Omega,Y] dependent part of qr contr. */ 
    real *res_omY_b;           /* [Omega,Y] dependent part of qr contr. */
    real *res_omZ_a;           /* [Omega,Z] dependent part of qr contr. */
    real *res_omZ_b;           /* [Omega,Z] dependent part of qr contr. */
    real *vya;                 /* VXC contribution at the grid points   */
    real *vyb;                 /* VXC contribution at the grid points   */ 
    real *vza;                 /* VXC contribution at the grid points   */
    real *vzb;                 /* VXC contribution at the grid points   */
    real *tmpa;                /* temporary variable for integration    */
    real *tmpb;                /* temporary variable for integration    */
    real *pref3a;              /* precomputed prefactors for grip points*/
    real *pref3b;              /* precomputed prefactors for grip points*/
    real *prefb2a;             /* precomputed prefactors for grip points*/
    real *prefb2b;             /* precomputed prefactors for grip points*/
    real *prefc2a;             /* precomputed prefactors for grip points*/
    real *prefc2b;             /* precomputed prefactors for grip points*/
} Quad_Open_Data;

static void
quad_open_lda_cb(DftIntegratorBl* grid, real * RESTRICT tmp,
            integer bllen, integer blstart, integer blend,
            Quad_Open_Data* data)
{
   static const real MONER = -1.0;

   integer i, j, k, isym, ibl, jbl;
   real * RESTRICT aos = grid->atv;
   real *om_a = data->res_omega_a;
   real *om_b = data->res_omega_b;
   real *omYa = data->res_omY_a;
   real *omYb = data->res_omY_b;
   real *omZa = data->res_omZ_a;
   real *omZb = data->res_omZ_b;
   real *tmpa = data->tmpa;
   real *tmpb = data->tmpb;
   real *pref3a = data->pref3a;
   real *pref3b = data->pref3b;
   real *prefb2a = data->prefb2a;
   real *prefb2b = data->prefb2b;
   real *prefc2a = data->prefc2a;
   real *prefc2b = data->prefc2b;
   FunDensProp dp = { 0 };
   FunThirdFuncDrv vxc;
   integer mat_size = DFT_BLLEN;

   /* compute vector of transformed densities */
   FSYM2(getexp_blocked_lda)(&data->symY, data->rhoa_y, grid->atv,
                        grid->bas_bl_cnt, grid->basblocks, &grid->shl_bl_cnt,
                        tmp, &bllen, data->rhowa_y);
   FSYM2(getexp_blocked_lda)(&data->symY, data->rhob_y, grid->atv,
                        grid->bas_bl_cnt, grid->basblocks, &grid->shl_bl_cnt,
                        tmp, &bllen, data->rhowb_y);
   FSYM2(getexp_blocked_lda)(&data->symZ, data->rhoa_z, grid->atv,
                        grid->bas_bl_cnt, grid->basblocks, &grid->shl_bl_cnt,
                        tmp, &bllen, data->rhowa_z);
   FSYM2(getexp_blocked_lda)(&data->symZ, data->rhob_z, grid->atv,
                        grid->bas_bl_cnt, grid->basblocks, &grid->shl_bl_cnt,
                        tmp, &bllen, data->rhowb_z);
   FSYM2(getexp_blocked_lda)(&data->symYZ, data->rhoa_yzzy, grid->atv,
                        grid->bas_bl_cnt, grid->basblocks, &grid->shl_bl_cnt,
                        tmp, &bllen, data->rhowa_yzzy);
   FSYM2(getexp_blocked_lda)(&data->symYZ, data->rhob_yzzy, grid->atv,
                        grid->bas_bl_cnt, grid->basblocks, &grid->shl_bl_cnt,
                        tmp, &bllen, data->rhowb_yzzy);
   /* set triplet perturbations signs */ 
   if (data->spinY==1) {
      FSYM(dscal)(&mat_size,&MONER,data->rhowb_y,&ONEI);
      FSYM(dscal)(&mat_size,&MONER,data->rhowb_yzzy,&ONEI);      
   }
   if (data->spinZ==1) {
      FSYM(dscal)(&mat_size,&MONER,data->rhowb_z,&ONEI);
      FSYM(dscal)(&mat_size,&MONER,data->rhowb_yzzy,&ONEI);      
   }
   /* compute derrivatives over grid points */
   for(i=blstart; i<blend; i++) {
     real weight = grid->weight[grid->curr_point+i];
     dp.rhoa = grid->r.ho.a[i];
     dp.rhob = grid->r.ho.b[i];
     drv3_clear(&vxc);
     selected_func->third(&vxc,weight,&dp);
     pref3a[i]  = vxc.df3000*data->rhowa_y[i]*data->rhowa_z[i]+
                  vxc.df2100*data->rhowb_y[i]*data->rhowa_z[i]+
                  vxc.df2100*data->rhowa_y[i]*data->rhowb_z[i]+
                  vxc.df1200*data->rhowb_y[i]*data->rhowb_z[i]+
                  vxc.df2000*data->rhowa_yzzy[i]+
                  vxc.df1100*data->rhowb_yzzy[i];
     pref3b[i]  = vxc.df0300*data->rhowb_y[i]*data->rhowb_z[i]+
                  vxc.df1200*data->rhowb_y[i]*data->rhowa_z[i]+
                  vxc.df1200*data->rhowa_y[i]*data->rhowb_z[i]+
                  vxc.df2100*data->rhowa_y[i]*data->rhowa_z[i]+
                  vxc.df0200*data->rhowb_yzzy[i]+
                  vxc.df1100*data->rhowa_yzzy[i];
     pref3a[i] *= 0.5; pref3b[i] *= 0.5;
     prefb2a[i] = vxc.df2000*data->rhowa_y[i]+vxc.df1100*data->rhowb_y[i];
     prefb2b[i] = vxc.df0200*data->rhowb_y[i]+vxc.df1100*data->rhowa_y[i];
     prefb2a[i] *= 0.5; prefb2b[i] *= 0.5;
     prefc2a[i] = vxc.df2000*data->rhowa_z[i]+vxc.df1100*data->rhowb_z[i];
     prefc2b[i] = vxc.df0200*data->rhowb_z[i]+vxc.df1100*data->rhowa_z[i];
     prefc2a[i] *= 0.5; prefc2b[i] *= 0.5;
   }

   /* main loop over symmetries */ 
   for(isym=0; isym<grid->nsym; isym++) {
     integer (*RESTRICT iblocks)[2] = BASBLOCK(grid,isym);
     integer ibl_cnt = grid->bas_bl_cnt[isym];
        for(ibl=0; ibl<ibl_cnt; ibl++)
            for(i=iblocks[ibl][0]-1; i<iblocks[ibl][1]; i++) {
                integer ioff = i*bllen;
                for(k=blstart; k<blend; k++) {
                    tmpa[k+ioff] = pref3a[k]*aos[k+ioff];
                    tmpb[k+ioff] = pref3b[k]*aos[k+ioff];
                    data->vza[k+ioff] = -prefb2a[k]*aos[k+ioff];
                    data->vzb[k+ioff] = -prefb2b[k]*aos[k+ioff];
                    data->vya[k+ioff] = -prefc2a[k]*aos[k+ioff];
                    data->vyb[k+ioff] = -prefc2b[k]*aos[k+ioff];
                }
            }
        /* Compute contributions to om, omY and omZ */
        for(ibl=0; ibl<ibl_cnt; ibl++) {
            for(i=iblocks[ibl][0]-1; i<iblocks[ibl][1]; i++) {
                integer jsym, jbl_cnt, ioff = i*inforb_.nbast;
                real * RESTRICT tmpai = tmpa + i*bllen;
                real * RESTRICT tmpbi = tmpb + i*bllen;
                real * RESTRICT vyai = data->vya + i*bllen;
                real * RESTRICT vybi = data->vyb + i*bllen;
                real * RESTRICT vzai = data->vza + i*bllen;
                real * RESTRICT vzbi = data->vzb + i*bllen;
                integer (*RESTRICT jblocks)[2];

                jsym = inforb_.muld2h[data->symYZ-1][isym]-1;
                jblocks = BASBLOCK(grid,jsym);
                jbl_cnt = grid->bas_bl_cnt[jsym];
                for(jbl=0; jbl<jbl_cnt; jbl++) {
                    for(j=jblocks[jbl][0]-1; j<jblocks[jbl][1]; j++) {
                        real * RESTRICT aosj = aos + j*bllen;
                        real sa = 0;
                        real sb = 0; 
                        for(k=blstart; k<blend; k++) {
                            sa += aosj[k]*tmpai[k];
                            sb += aosj[k]*tmpbi[k];
                        }  
                        om_a[j+ioff] += sa;
                        om_b[j+ioff] += sb; 
                    }
                }
                
                jsym = inforb_.muld2h[data->symZ-1][isym]-1;
                jblocks = BASBLOCK(grid,jsym);
                jbl_cnt = grid->bas_bl_cnt[jsym];
                for(jbl=0; jbl<jbl_cnt; jbl++) {
                    for(j=jblocks[jbl][0]-1; j<jblocks[jbl][1]; j++) {
                        real * RESTRICT aosj = aos + j*bllen;
                        real sa = 0;
                        real sb = 0;
                        for(k=blstart; k<blend; k++) {
                             sa += aosj[k]*vyai[k];
                             sb += aosj[k]*vybi[k];
                        }
                        omYa[j+ioff] += sa;
                        omYb[j+ioff] += sb;
                    }
                }
                
                jsym = inforb_.muld2h[data->symY-1][isym]-1;
                jblocks = BASBLOCK(grid,jsym);
                jbl_cnt = grid->bas_bl_cnt[jsym];
                for(jbl=0; jbl<jbl_cnt; jbl++) {
                    for(j=jblocks[jbl][0]-1; j<jblocks[jbl][1]; j++) {
                        real * RESTRICT aosj = aos + j*bllen;
                        real sa = 0;
                        real sb = 0;
                        for(k=blstart; k<blend; k++) {
                            sa += aosj[k]*vzai[k];
                            sb += aosj[k]*vzbi[k];
                        }
                        omZa[j+ioff] += sa;
                        omZb[j+ioff] += sb; 
                    }
                }
            } /* Start next column */
        }/* one batch of blocks */
    } /* Loop over symmetries */
}


static void
quad_open_gga_cb(DftIntegratorBl* grid, real * RESTRICT tmp,
            integer bllen, integer blstart, integer blend,
            Quad_Open_Data* data)
{

    static const real MONER = -1.0;

    integer i,j, k, ibl, jbl, isym; 
    real * RESTRICT aos = grid->atv;
    real * RESTRICT aox = grid->atv+bllen*inforb_.nbast;
    real * RESTRICT aoy = grid->atv+bllen*inforb_.nbast*2;
    real * RESTRICT aoz = grid->atv+bllen*inforb_.nbast*3;
    real *om_a = data->res_omega_a;
    real *om_b = data->res_omega_b;
    real *omYa = data->res_omY_a;
    real *omYb = data->res_omY_b;
    real *omZa = data->res_omZ_a;
    real *omZb = data->res_omZ_b;
    real *vya  = data->vya;
    real *vyb  = data->vyb;
    real *vza  = data->vza;
    real *vzb  = data->vzb;
    real *tmpa = data->tmpa;
    real *tmpb = data->tmpb;  
    real (*pref3a)[4]   = (real (*)[4])data->pref3a;
    real (*pref3b)[4]   = (real (*)[4])data->pref3b; 
    real (*prefb2a)[4]  = (real (*)[4])data->prefb2a;
    real (*prefb2b)[4]  = (real (*)[4])data->prefb2b;
    real (*prefc2a)[4]  = (real (*)[4])data->prefc2a;
    real (*prefc2b)[4]  = (real (*)[4])data->prefc2b; 
    FunDensProp dp = { 0 };
    FunThirdFuncDrv vxc;    
    integer mat_size = 4*DFT_BLLEN;
    real zetaYa, zetaYb, zetaYab;
    real zetaZa, zetaZb, zetaZab;
    real zetaYZZYa, zetaYZZYb, zetaYZZYab; 
    real trYZaa, trYZbb, trYZab; 
    real fza2f, fzb2f, fg2f; 
    real fza1fy, fzb1fy, fg1fy; 
    real fza1fz, fzb1fz, fg1fz;   
    real fza0f, fzb0f, fg0f;   
    real fa1y, fb1y, fg1y;
    real fa1z, fb1z, fg1z;
    
    /* compute vectors of transformed densities and their gradients */
    FSYM2(getexp_blocked_gga)(&data->symY, data->rhoa_y, grid->atv,
                              grid->bas_bl_cnt, grid->basblocks, &grid->shl_bl_cnt,
                              tmp, &bllen, (double(*)[4])data->rhowa_y);
    FSYM2(getexp_blocked_gga)(&data->symY, data->rhob_y, grid->atv,
                              grid->bas_bl_cnt, grid->basblocks, &grid->shl_bl_cnt,
                              tmp, &bllen, (double(*)[4])data->rhowb_y);
    FSYM2(getexp_blocked_gga)(&data->symZ, data->rhoa_z, grid->atv,
                              grid->bas_bl_cnt, grid->basblocks, &grid->shl_bl_cnt,
                              tmp, &bllen, (double(*)[4])data->rhowa_z);
    FSYM2(getexp_blocked_gga)(&data->symZ, data->rhob_z, grid->atv,
                              grid->bas_bl_cnt, grid->basblocks, &grid->shl_bl_cnt,
                              tmp, &bllen, (double(*)[4])data->rhowb_z);
    FSYM2(getexp_blocked_gga)(&data->symYZ, data->rhoa_yzzy, grid->atv,
                              grid->bas_bl_cnt, grid->basblocks, &grid->shl_bl_cnt,
                              tmp, &bllen, (double(*)[4])data->rhowa_yzzy);
    FSYM2(getexp_blocked_gga)(&data->symYZ, data->rhob_yzzy, grid->atv,
                              grid->bas_bl_cnt, grid->basblocks, &grid->shl_bl_cnt,
                              tmp, &bllen, (double(*)[4])data->rhowb_yzzy);
    /* set triplet perturbations signs */ 
    if (data->spinY==1) {
      FSYM(dscal)(&mat_size,&MONER,data->rhowb_y,&ONEI);
      FSYM(dscal)(&mat_size,&MONER,data->rhowb_yzzy,&ONEI);      
   }
   if (data->spinZ==1) {
      FSYM(dscal)(&mat_size,&MONER,data->rhowb_z,&ONEI);
      FSYM(dscal)(&mat_size,&MONER,data->rhowb_yzzy,&ONEI);      
   }  
    real (*trYa)[4]  = (real (*)[4])data->rhowa_y;
    real (*trYb)[4]  = (real (*)[4])data->rhowb_y;
    real (*trZa)[4]  = (real (*)[4])data->rhowa_z;
    real (*trZb)[4]  = (real (*)[4])data->rhowb_z;
    real (*trYZZYa)[4]  = (real (*)[4])data->rhowa_yzzy;
    real (*trYZZYb)[4]  = (real (*)[4])data->rhowb_yzzy;

    for(i=blstart; i<blend; i++) {
     real weight = grid->weight[grid->curr_point+i];
     dp.rhoa = grid->r.ho.a[i];
     dp.rhob = grid->r.ho.b[i];
     dp.grada = sqrt(grid->g.rad.a [i][0]*grid->g.rad.a [i][0]+
                     grid->g.rad.a [i][1]*grid->g.rad.a [i][1]+
                     grid->g.rad.a [i][2]*grid->g.rad.a [i][2]);
     dp.gradb = sqrt(grid->g.rad.b [i][0]*grid->g.rad.b [i][0]+
                     grid->g.rad.b [i][1]*grid->g.rad.b [i][1]+
                     grid->g.rad.b [i][2]*grid->g.rad.b [i][2]);
     dp.gradab = grid->g.rad.a [i][0]*grid->g.rad.b [i][0]+
                 grid->g.rad.a [i][1]*grid->g.rad.b [i][1]+
                 grid->g.rad.a [i][2]*grid->g.rad.b [i][2];
     if(dp.grada<1e-20) dp.grada = 1e-20;
     if(dp.gradb<1e-20) dp.gradb = 1e-20;
     /* Clear derrivatives and compute them for current point*/
     drv3_clear(&vxc);
     selected_func->third(&vxc,weight,&dp); 
     /* First order variations of denisty gradient */
     zetaYa  = trYa[i][1]*grid->g.rad.a[i][0]+trYa[i][2]*grid->g.rad.a[i][1]+
               trYa[i][3]*grid->g.rad.a[i][2];
     zetaYb  = trYb[i][1]*grid->g.rad.b[i][0]+trYb[i][2]*grid->g.rad.b[i][1]+
               trYb[i][3]*grid->g.rad.b[i][2];
     zetaZa  = trZa[i][1]*grid->g.rad.a[i][0]+trZa[i][2]*grid->g.rad.a[i][1]+
               trZa[i][3]*grid->g.rad.a[i][2];
     zetaZb  = trZb[i][1]*grid->g.rad.b[i][0]+trZb[i][2]*grid->g.rad.b[i][1]+
               trZb[i][3]*grid->g.rad.b[i][2];
     zetaYa  = zetaYa/dp.grada; zetaYb = zetaYb/dp.gradb;
     zetaZa  = zetaZa/dp.grada; zetaZb = zetaZb/dp.gradb;
     zetaYab = trYa[i][1]*grid->g.rad.b[i][0]+trYa[i][2]*grid->g.rad.b[i][1]+
               trYa[i][3]*grid->g.rad.b[i][2]+trYb[i][1]*grid->g.rad.a[i][0]+
               trYb[i][2]*grid->g.rad.a[i][1]+trYb[i][3]*grid->g.rad.a[i][2];
     zetaZab = trZa[i][1]*grid->g.rad.b[i][0]+trZa[i][2]*grid->g.rad.b[i][1]+
               trZa[i][3]*grid->g.rad.b[i][2]+trZb[i][1]*grid->g.rad.a[i][0]+
               trZb[i][2]*grid->g.rad.a[i][1]+trZb[i][3]*grid->g.rad.a[i][2]; 
     /* Second order variations of density gradient */
     zetaYZZYa  = trYZZYa[i][1]*grid->g.rad.a[i][0]+
                  trYZZYa[i][2]*grid->g.rad.a[i][1]+
                  trYZZYa[i][3]*grid->g.rad.a[i][2];
     zetaYZZYb  = trYZZYb[i][1]*grid->g.rad.b[i][0]+
                  trYZZYb[i][2]*grid->g.rad.b[i][1]+
                  trYZZYb[i][3]*grid->g.rad.b[i][2];
     zetaYZZYa  = zetaYZZYa/dp.grada; 
     zetaYZZYb  = zetaYZZYb/dp.gradb;  
     zetaYZZYab = trYZZYa[i][1]*grid->g.rad.b[i][0]+trYZZYa[i][2]*grid->g.rad.b[i][1]+
                  trYZZYa[i][3]*grid->g.rad.b[i][2]+trYZZYb[i][1]*grid->g.rad.a[i][0]+
                  trYZZYb[i][2]*grid->g.rad.a[i][1]+trYZZYb[i][3]*grid->g.rad.a[i][2];
     trYZaa     = trYa[i][1]*trZa[i][1]+trYa[i][2]*trZa[i][2]+trYa[i][3]*trZa[i][3];
     trYZbb     = trYb[i][1]*trZb[i][1]+trYb[i][2]*trZb[i][2]+trYb[i][3]*trZb[i][3];
     trYZab     = trYa[i][1]*trZb[i][1]+trYa[i][2]*trZb[i][2]+trYa[i][3]*trZb[i][3]+ 
                  trYb[i][1]*trZa[i][1]+trYb[i][2]*trZa[i][2]+trYb[i][3]*trZa[i][3];
     /*****************************************************************/
     /*** Density operator times 2-th order variation of functional ***/
     /*****************************************************************/
     /* FRRR part: */
     pref3a[i][0]  = vxc.df3000*trYa[i][0]*trZa[i][0]+vxc.df2100*trYb[i][0]*trZa[i][0]+
                     vxc.df2100*trYa[i][0]*trZb[i][0]+vxc.df1200*trYb[i][0]*trZb[i][0];
     pref3b[i][0]  = vxc.df0300*trYb[i][0]*trZb[i][0]+vxc.df1200*trYb[i][0]*trZa[i][0]+
                     vxc.df1200*trYa[i][0]*trZb[i][0]+vxc.df2100*trYa[i][0]*trZa[i][0];
     /* FRRZ part: */  
     pref3a[i][0] += vxc.df2010*trYa[i][0]*zetaZa+vxc.df1110*trYb[i][0]*zetaZa+
                     vxc.df2001*trYa[i][0]*zetaZb+vxc.df1101*trYb[i][0]*zetaZb+
                     vxc.df2010*trZa[i][0]*zetaYa+vxc.df1110*trZb[i][0]*zetaYa+
                     vxc.df2001*trZa[i][0]*zetaYb+vxc.df1101*trZb[i][0]*zetaYb;
     pref3b[i][0] += vxc.df0201*trYb[i][0]*zetaZb+vxc.df1101*trYa[i][0]*zetaZb+
                     vxc.df0210*trYb[i][0]*zetaZa+vxc.df1101*trYa[i][0]*zetaZa+
                     vxc.df0201*trZb[i][0]*zetaYb+vxc.df1101*trZa[i][0]*zetaYb+
                     vxc.df0210*trZb[i][0]*zetaYa+vxc.df1101*trZa[i][0]*zetaYa;
     /* FRRG part: */  
     pref3a[i][0] += vxc.df20001*trYa[i][0]*zetaZab+vxc.df11001*trYb[i][0]*zetaZab+
                     vxc.df20001*trZa[i][0]*zetaYab+vxc.df11001*trZb[i][0]*zetaYab;
     pref3b[i][0] += vxc.df02001*trYb[i][0]*zetaZab+vxc.df11001*trYa[i][0]*zetaZab+
                     vxc.df02001*trZb[i][0]*zetaYab+vxc.df11001*trZa[i][0]*zetaYab;
     /* FRZZ part: */  
     pref3a[i][0] += vxc.df1020*zetaYa*zetaZa+vxc.df1011*zetaYa*zetaZb+
                     vxc.df1011*zetaYb*zetaZa+vxc.df1002*zetaYb*zetaZb;
     pref3b[i][0] += vxc.df0102*zetaYb*zetaZb+vxc.df0111*zetaYb*zetaZa+
                     vxc.df0111*zetaYa*zetaZb+vxc.df0120*zetaYa*zetaZa;
     /* FRZG part: not implemented in Dalton */  
     /* FRGG part: not implemented in Dalton */    
     /* FRR  part: */
     pref3a[i][0] += vxc.df2000*trYZZYa[i][0]+vxc.df1100*trYZZYb[i][0];
     pref3b[i][0] += vxc.df0200*trYZZYb[i][0]+vxc.df1100*trYZZYa[i][0];
     /* FRZ  part: */
     pref3a[i][0] += vxc.df1010*zetaYZZYa+vxc.df1001*zetaYZZYb;
     pref3b[i][0] += vxc.df0101*zetaYZZYb+vxc.df0110*zetaYZZYa;
     pref3a[i][0] += vxc.df1010*trYZaa/dp.grada+vxc.df1001*trYZbb/dp.gradb;
     pref3b[i][0] += vxc.df0101*trYZbb/dp.gradb+vxc.df0110*trYZaa/dp.grada;
     pref3a[i][0] += -(vxc.df1010*zetaYa*zetaZa/dp.grada+vxc.df1001*zetaYb*zetaZb/dp.gradb);
     pref3b[i][0] += -(vxc.df0110*zetaYa*zetaZa/dp.grada+vxc.df0101*zetaYb*zetaZb/dp.gradb);
     /* FRG  part: */
     pref3a[i][0] += vxc.df10001*zetaYZZYab+vxc.df10001*trYZab;
     pref3b[i][0] += vxc.df01001*zetaYZZYab+vxc.df01001*trYZab;
     pref3a[i][0] *= 0.5;  pref3b[i][0] *= 0.5;    
     /***************************************************************************/
     /*** Density gradient operator times 2-th order variation of functional ****/
     /***************************************************************************/
     /* FZZZ part: */ 
     fza2f  = vxc.df0030*zetaYa*zetaZa+vxc.df0021*zetaYa*zetaZb+
              vxc.df0021*zetaYb*zetaZa+vxc.df0012*zetaYb*zetaZb;
     fzb2f  = vxc.df0003*zetaYb*zetaZb+vxc.df0012*zetaYb*zetaZa+
              vxc.df0012*zetaYa*zetaZb+vxc.df0021*zetaYa*zetaZa;
     /* FZZR part: */  
     fza2f += vxc.df1020*trYa[i][0]*zetaZa+vxc.df1011*trYa[i][0]*zetaZb+
              vxc.df0120*trYb[i][0]*zetaZa+vxc.df0111*trYb[i][0]*zetaZb+
              vxc.df1020*trZa[i][0]*zetaYa+vxc.df1011*trZa[i][0]*zetaYb+
              vxc.df0120*trZb[i][0]*zetaYa+vxc.df0111*trZb[i][0]*zetaYa;
     fzb2f += vxc.df0102*trYb[i][0]*zetaZb+vxc.df0111*trYb[i][0]*zetaZa+
              vxc.df1002*trYa[i][0]*zetaZb+vxc.df1011*trYa[i][0]*zetaZa+
              vxc.df0102*trZb[i][0]*zetaYb+vxc.df0111*trZb[i][0]*zetaYa+
              vxc.df1002*trZa[i][0]*zetaYb+vxc.df1011*trZa[i][0]*zetaYa;
     /* FZZG part: not implemented in Dalton */  
     /* FZRR part: */
     fza2f += vxc.df2010*trYa[i][0]*trZa[i][0]+vxc.df1110*trYa[i][0]*trZb[i][0]+
              vxc.df1110*trZa[i][0]*trYb[i][0]+vxc.df0210*trZb[i][0]*trYb[i][0];
     fzb2f += vxc.df0201*trYb[i][0]*trZb[i][0]+vxc.df1101*trYa[i][0]*trZb[i][0]+
              vxc.df1101*trZa[i][0]*trYb[i][0]+vxc.df2001*trYa[i][0]*trZa[i][0];
     /* FZRG part: not implemented in Dalton */  
     /* FZZ  part:  */
     fza2f += vxc.df0020*zetaYZZYa+vxc.df0011*zetaYZZYb;
     fzb2f += vxc.df0011*zetaYZZYa+vxc.df0002*zetaYZZYb;
     fza2f += vxc.df0020*trYZaa/dp.grada+vxc.df0011*trYZbb/dp.gradb;
     fzb2f += vxc.df0002*trYZbb/dp.gradb+vxc.df0011*trYZaa/dp.grada;
     // Need additional testting summed in closed shell case with other terms !!!!
     fza2f += -(vxc.df0020*zetaYa*zetaZa/dp.grada+vxc.df0011*zetaYb*zetaZb/dp.gradb);
     fzb2f += -(vxc.df0011*zetaYa*zetaZa/dp.grada+vxc.df0002*zetaYb*zetaZb/dp.gradb);
     /* FZR part:  */
     fza2f += vxc.df1010*trYZZYa[i][0]+vxc.df0110*trYZZYb[i][0];
     fzb2f += vxc.df1001*trYZZYa[i][0]+vxc.df0101*trYZZYb[i][0];
     /* FZG part: not implemented in Dalton */
     /* Scale contribution with density gradient operator prefactor */
     fza2f *= 1.0/dp.grada; fzb2f *= 1.0/dp.gradb;
     /* distribute this contribution: */
     pref3a[i][1] = fza2f*grid->g.rad.a[i][0]; 
     pref3a[i][2] = fza2f*grid->g.rad.a[i][1];
     pref3a[i][3] = fza2f*grid->g.rad.a[i][2];
     pref3b[i][1] = fzb2f*grid->g.rad.b[i][0]; 
     pref3b[i][2] = fzb2f*grid->g.rad.b[i][1];
     pref3b[i][3] = fzb2f*grid->g.rad.b[i][2];
     /*********************************************************************************/
     /*** Mixed density gradient operator times 2-th order variation of functional ****/
     /*********************************************************************************/
     /* FGGG part: not implemented in Dalton */
     /* FGGZ part: not implemented in Dalton */  
     /* FGGR part: not implemented in Dalton */
     /* FGZZ part: not implemented in Dalton */
     /* FGRR part:  */
     fg2f  = vxc.df20001*trYa[i][0]*trZa[i][0]+vxc.df11001*trYa[i][0]*trZb[i][0]+
             vxc.df11001*trYb[i][0]*trZa[i][0]+vxc.df02001*trYb[i][0]*trZb[i][0];
     /* FGG  part: not implemented in Dalton */ 
     /* FGZ  part: not implemented in Dalton */ 
     /* FGR  part:  */
     fg2f +=  vxc.df10001*trYZZYa[i][0]+vxc.df01001*trYZZYb[i][0];  
     /* distribute this contribution: */
     pref3a[i][1] += fg2f*grid->g.rad.b[i][0]; 
     pref3a[i][2] += fg2f*grid->g.rad.b[i][1];
     pref3a[i][3] += fg2f*grid->g.rad.b[i][2];
     pref3b[i][1] += fg2f*grid->g.rad.a[i][0]; 
     pref3b[i][2] += fg2f*grid->g.rad.a[i][1];
     pref3b[i][3] += fg2f*grid->g.rad.a[i][2];
     /**********************************************************************************/
     /*** Density gradient operator variation (Z) times variation of functional (Y) ****/
     /**********************************************************************************/
     /* FZZ  part:  */ 
     fza1fy  = vxc.df0020*zetaYa+vxc.df0011*zetaYb;
     fzb1fy  = vxc.df0002*zetaYb+vxc.df0011*zetaYa;
     /* FZR  part:  */ 
     fza1fy += vxc.df1010*trYa[i][0]+vxc.df0110*trYb[i][0];
     fzb1fy += vxc.df1001*trYa[i][0]+vxc.df0101*trYb[i][0];    
     /* Scale contribution with density gradient operator prefactor */
     fza1fy *= 1.0/dp.grada; fzb1fy *= 1.0/dp.gradb;    
     /* distribute this contribution: */
     pref3a[i][1] += fza1fy*trZa[i][1];
     pref3a[i][2] += fza1fy*trZa[i][2];
     pref3a[i][3] += fza1fy*trZa[i][3];
     pref3b[i][1] += fzb1fy*trZb[i][1];
     pref3b[i][2] += fzb1fy*trZb[i][2];
     pref3b[i][3] += fzb1fy*trZb[i][3];
     /* Second part: scale this contribution */  
     fza1fy *= -zetaZa/dp.grada; 
     fzb1fy *= -zetaZb/dp.gradb;     
     /* distribute this contribution */
     pref3a[i][1] += fza1fy*grid->g.rad.a[i][0]; 
     pref3a[i][2] += fza1fy*grid->g.rad.a[i][1];
     pref3a[i][3] += fza1fy*grid->g.rad.a[i][2];
     pref3b[i][1] += fzb1fy*grid->g.rad.b[i][0]; 
     pref3b[i][2] += fzb1fy*grid->g.rad.b[i][1];
     pref3b[i][3] += fzb1fy*grid->g.rad.b[i][2];
     /**********************************************************************************/
     /*** Density gradient operator variation (Y) times variation of functional (Z) ****/
     /**********************************************************************************/
     /* FZZ  part:  */ 
     fza1fz  = vxc.df0020*zetaZa+vxc.df0011*zetaZb;
     fzb1fz  = vxc.df0002*zetaZb+vxc.df0011*zetaZa;
     /* FZR  part:  */ 
     fza1fz += vxc.df1010*trZa[i][0]+vxc.df0110*trZb[i][0];
     fzb1fz += vxc.df1001*trZa[i][0]+vxc.df0101*trZb[i][0];    
     /* Scale contribution with density gradient operator prefactor */
     fza1fz *= 1.0/dp.grada; fzb1fz *= 1.0/dp.gradb;    
     /* distribute this contribution: */
     pref3a[i][1] += fza1fz*trYa[i][1];
     pref3a[i][2] += fza1fz*trYa[i][2];
     pref3a[i][3] += fza1fz*trYa[i][3];
     pref3b[i][1] += fzb1fz*trYb[i][1];
     pref3b[i][2] += fzb1fz*trYb[i][2];
     pref3b[i][3] += fzb1fz*trYb[i][3];
     /* Second part: scale this contribution */  
     fza1fz *= -zetaYa/dp.grada; 
     fzb1fz *= -zetaYb/dp.gradb;     
     /* distribute this contribution */
     pref3a[i][1] += fza1fz*grid->g.rad.a[i][0]; 
     pref3a[i][2] += fza1fz*grid->g.rad.a[i][1];
     pref3a[i][3] += fza1fz*grid->g.rad.a[i][2];
     pref3b[i][1] += fzb1fz*grid->g.rad.b[i][0]; 
     pref3b[i][2] += fzb1fz*grid->g.rad.b[i][1];
     pref3b[i][3] += fzb1fz*grid->g.rad.b[i][2];
     /****************************************************************************************/
     /*** Mixed density gradient operator variation (Z) times variation of functional (Y) ****/
     /****************************************************************************************/
     /* FGG not implemented in Dalton */
     /* FGZ not implemented in Dalton */ 
     /* FGR  part: */ 
     fg1fy  = vxc.df10001*trYa[i][0]+vxc.df10001*trYb[i][0]; 
     //fg1fy *= 2.0; 
     pref3a[i][1] += fg1fy*trZb[i][1];
     pref3a[i][2] += fg1fy*trZb[i][2];
     pref3a[i][3] += fg1fy*trZb[i][3];
     pref3b[i][1] += fg1fy*trZa[i][1];
     pref3b[i][2] += fg1fy*trZa[i][2];
     pref3b[i][3] += fg1fy*trZa[i][3]; 
     /****************************************************************************************/
     /*** Mixed density gradient operator variation (Y) times variation of functional (Z) ****/
     /****************************************************************************************/
     /* FGG not implemented in Dalton */
     /* FGZ not implemented in Dalton */ 
     /* FGR  part: */ 
     fg1fz  = vxc.df10001*trZa[i][0]+vxc.df10001*trZb[i][0]; 
     pref3a[i][1] += fg1fz*trYb[i][1];
     pref3a[i][2] += fg1fz*trYb[i][2];
     pref3a[i][3] += fg1fz*trYb[i][3];
     pref3b[i][1] += fg1fz*trYa[i][1];
     pref3b[i][2] += fg1fz*trYa[i][2];
     pref3b[i][3] += fg1fz*trYa[i][3]; 
     /**************************************************************/
     /*** Second order variation of density gradient operator    ***/
     /**************************************************************/
     fza0f = vxc.df0010/dp.grada;  
     fzb0f = vxc.df0001/dp.gradb;  
     pref3a[i][1] += fza0f*trYZZYa[i][1];
     pref3a[i][2] += fza0f*trYZZYa[i][2]; 
     pref3a[i][3] += fza0f*trYZZYa[i][3];
     pref3b[i][1] += fzb0f*trYZZYb[i][1];
     pref3b[i][2] += fzb0f*trYZZYb[i][2];
     pref3b[i][3] += fzb0f*trYZZYb[i][3];
     /* second part of this contribution */
     fza0f = -vxc.df0010*zetaYZZYa/(dp.grada*dp.grada);  
     fzb0f = -vxc.df0001*zetaYZZYb/(dp.gradb*dp.gradb);  
     pref3a[i][1] += fza0f*grid->g.rad.a[i][0];
     pref3a[i][2] += fza0f*grid->g.rad.a[i][1];
     pref3a[i][3] += fza0f*grid->g.rad.a[i][2];
     pref3b[i][1] += fzb0f*grid->g.rad.b[i][0];
     pref3b[i][2] += fzb0f*grid->g.rad.b[i][1];
     pref3b[i][3] += fzb0f*grid->g.rad.b[i][2];
     /* third part of this contribution */
     fza0f = -vxc.df0010*trYZaa/(dp.grada*dp.grada*dp.grada); 
     fzb0f = -vxc.df0001*trYZbb/(dp.gradb*dp.gradb*dp.gradb); 
     pref3a[i][1] += fza0f*grid->g.rad.a[i][0];
     pref3a[i][2] += fza0f*grid->g.rad.a[i][1];
     pref3a[i][3] += fza0f*grid->g.rad.a[i][2];
     pref3b[i][1] += fzb0f*grid->g.rad.b[i][0];
     pref3b[i][2] += fzb0f*grid->g.rad.b[i][1];
     pref3b[i][3] += fzb0f*grid->g.rad.b[i][2];
     /* fourth part of this contribution */
     fza0f = -vxc.df0010*zetaYa/(dp.grada*dp.grada); 
     fzb0f = -vxc.df0001*zetaYb/(dp.gradb*dp.gradb); 
     pref3a[i][1] += fza0f*trZa[i][1];
     pref3a[i][2] += fza0f*trZa[i][2];
     pref3a[i][3] += fza0f*trZa[i][3];
     pref3b[i][1] += fzb0f*trZb[i][1];
     pref3b[i][2] += fzb0f*trZb[i][2];
     pref3b[i][3] += fzb0f*trZb[i][3];
     /* fifth part of this contribution */
     fza0f = -vxc.df0010*zetaZa/(dp.grada*dp.grada); 
     fzb0f = -vxc.df0001*zetaZb/(dp.gradb*dp.gradb); 
     pref3a[i][1] += fza0f*trYa[i][1];
     pref3a[i][2] += fza0f*trYa[i][2];
     pref3a[i][3] += fza0f*trYa[i][3];
     pref3b[i][1] += fzb0f*trYb[i][1];
     pref3b[i][2] += fzb0f*trYb[i][2];
     pref3b[i][3] += fzb0f*trYb[i][3];
     /* third part od this contribution with mixed perturbations */
     fza0f = 3.0*vxc.df0010*zetaYa*zetaZa/(dp.grada*dp.grada*dp.grada);  
     fzb0f = 3.0*vxc.df0001*zetaYb*zetaZb/(dp.gradb*dp.gradb*dp.gradb);  
     pref3a[i][1] += fza0f*grid->g.rad.a[i][0];
     pref3a[i][2] += fza0f*grid->g.rad.a[i][1];
     pref3a[i][3] += fza0f*grid->g.rad.a[i][2];
     pref3b[i][1] += fzb0f*grid->g.rad.b[i][0];
     pref3b[i][2] += fzb0f*grid->g.rad.b[i][1];
     pref3b[i][3] += fzb0f*grid->g.rad.b[i][2];
     /********************************************************************/
     /*** Second order variation of mixed density gradient operator    ***/
     /********************************************************************/
     fg0f = vxc.df00001; 
     pref3a[i][1] += fg0f*trYZZYb[i][1];
     pref3a[i][2] += fg0f*trYZZYb[i][2]; 
     pref3a[i][3] += fg0f*trYZZYb[i][3];
     pref3b[i][1] += fg0f*trYZZYa[i][1];
     pref3b[i][2] += fg0f*trYZZYa[i][2];
     pref3b[i][3] += fg0f*trYZZYa[i][3];
     /*************************************/
     /*** Second order (Y) contribution ***/
     /*************************************/
     prefb2a[i][0] = vxc.df2000*trYa[i][0]+vxc.df1100*trYb[i][0]+
                     vxc.df1010*zetaYa+vxc.df1001*zetaYb+
                     vxc.df10001*zetaYab;  
     prefb2b[i][0] = vxc.df0200*trYb[i][0]+vxc.df1100*trYa[i][0]+
                     vxc.df0110*zetaYa+vxc.df0101*zetaYb+ 
                     vxc.df01001*zetaYab;
     prefb2a[i][0] *= 0.5; prefb2b[i][0] *= 0.5; 
     fa1y  = (vxc.df0020*zetaYa+vxc.df0011*zetaYb+vxc.df1010*trYa[i][0]+
              vxc.df0110*trYb[i][0])/dp.grada-vxc.df0010*zetaYa/(dp.grada*dp.grada);
     fb1y  = (vxc.df0002*zetaYb+vxc.df0011*zetaYa+vxc.df1001*trYa[i][0]+
              vxc.df0101*trYb[i][0])/dp.gradb-vxc.df0001*zetaYb/(dp.gradb*dp.gradb);    
     prefb2a[i][1] = fa1y*grid->g.rad.a[i][0];
     prefb2a[i][2] = fa1y*grid->g.rad.a[i][1];
     prefb2a[i][3] = fa1y*grid->g.rad.a[i][2];
     prefb2b[i][1] = fb1y*grid->g.rad.b[i][0];
     prefb2b[i][2] = fb1y*grid->g.rad.b[i][1];
     prefb2b[i][3] = fb1y*grid->g.rad.b[i][2];
     fg1y  = vxc.df10001*trYa[i][0]+vxc.df01001*trYb[i][0]; 
     prefb2a[i][1] += fg1y*grid->g.rad.b[i][0];
     prefb2a[i][2] += fg1y*grid->g.rad.b[i][1];
     prefb2a[i][3] += fg1y*grid->g.rad.b[i][2];
     prefb2b[i][1] += fg1y*grid->g.rad.a[i][0];
     prefb2b[i][2] += fg1y*grid->g.rad.a[i][1];
     prefb2b[i][3] += fg1y*grid->g.rad.a[i][2]; 
     fa1y  = vxc.df0010/dp.grada;
     fb1y  = vxc.df0001/dp.gradb;
     prefb2a[i][1] += fa1y*trYa[i][1];
     prefb2a[i][2] += fa1y*trYa[i][2];
     prefb2a[i][3] += fa1y*trYa[i][3];
     prefb2b[i][1] += fb1y*trYb[i][1];
     prefb2b[i][2] += fb1y*trYb[i][2];
     prefb2b[i][3] += fb1y*trYb[i][3];
     fg1y  = vxc.df00001; 
     prefb2a[i][1] += fg1y*trYb[i][1];
     prefb2a[i][2] += fg1y*trYb[i][2];
     prefb2a[i][3] += fg1y*trYb[i][3];
     prefb2b[i][1] += fg1y*trYa[i][1];
     prefb2b[i][2] += fg1y*trYa[i][2];
     prefb2b[i][3] += fg1y*trYa[i][3];
     /*************************************/
     /*** Second order (Z) contribution ***/
     /*************************************/
     prefc2a[i][0] = vxc.df2000*trZa[i][0]+vxc.df1100*trZb[i][0]+
                     vxc.df1010*zetaZa+vxc.df1001*zetaZb+
                     vxc.df10001*zetaZab;  
     prefc2b[i][0] = vxc.df0200*trZb[i][0]+vxc.df1100*trZa[i][0]+
                     vxc.df0110*zetaZa+vxc.df0101*zetaZb+ 
                     vxc.df01001*zetaZab;
     prefc2a[i][0] *= 0.5; prefc2b[i][0] *= 0.5;
     fa1z  = (vxc.df0020*zetaZa+vxc.df0011*zetaZb+vxc.df1010*trZa[i][0]+
              vxc.df0110*trZb[i][0])/dp.grada-vxc.df0010*zetaZa/(dp.grada*dp.grada);
     fb1z  = (vxc.df0002*zetaZb+vxc.df0011*zetaZa+vxc.df1001*trZa[i][0]+
              vxc.df0101*trZb[i][0])/dp.gradb-vxc.df0001*zetaZb/(dp.gradb*dp.gradb);    
     prefc2a[i][1] = fa1z*grid->g.rad.a[i][0];
     prefc2a[i][2] = fa1z*grid->g.rad.a[i][1];
     prefc2a[i][3] = fa1z*grid->g.rad.a[i][2];
     prefc2b[i][1] = fb1z*grid->g.rad.b[i][0];
     prefc2b[i][2] = fb1z*grid->g.rad.b[i][1];
     prefc2b[i][3] = fb1z*grid->g.rad.b[i][2];
     fg1z  = vxc.df10001*trZa[i][0]+vxc.df01001*trZb[i][0]; 
     prefc2a[i][1] += fg1z*grid->g.rad.b[i][0];
     prefc2a[i][2] += fg1z*grid->g.rad.b[i][1];
     prefc2a[i][3] += fg1z*grid->g.rad.b[i][2];
     prefc2b[i][1] += fg1z*grid->g.rad.a[i][0];
     prefc2b[i][2] += fg1z*grid->g.rad.a[i][1];
     prefc2b[i][3] += fg1z*grid->g.rad.a[i][2]; 
     fa1z  = vxc.df0010/dp.grada;
     fb1z  = vxc.df0001/dp.gradb;
     prefc2a[i][1] += fa1z*trZa[i][1];
     prefc2a[i][2] += fa1z*trZa[i][2];
     prefc2a[i][3] += fa1z*trZa[i][3];
     prefc2b[i][1] += fb1z*trZb[i][1];
     prefc2b[i][2] += fb1z*trZb[i][2];
     prefc2b[i][3] += fb1z*trZb[i][3];
     fg1z  = vxc.df00001; 
     prefc2a[i][1] += fg1z*trZb[i][1];
     prefc2a[i][2] += fg1z*trZb[i][2];
     prefc2a[i][3] += fg1z*trZb[i][3];
     prefc2b[i][1] += fg1z*trZa[i][1];
     prefc2b[i][2] += fg1z*trZa[i][2];
     prefc2b[i][3] += fg1z*trZa[i][3];  
   }   
   /* Compute VXC[3] contribution */
    for(isym=0; isym<grid->nsym; isym++) {
        integer (*RESTRICT iblocks)[2] = BASBLOCK(grid,isym);
        integer ibl_cnt = grid->bas_bl_cnt[isym];
        for(ibl=0; ibl<ibl_cnt; ibl++)
            for(i=iblocks[ibl][0]-1; i<iblocks[ibl][1]; i++) {
                real * RESTRICT a0 = aos + i*bllen;
                real * RESTRICT ax = aox + i*bllen;
                real * RESTRICT ay = aoy + i*bllen;
                real * RESTRICT az = aoz + i*bllen;
                real * RESTRICT vyai  = vya + i*bllen;
                real * RESTRICT vybi  = vyb + i*bllen;
                real * RESTRICT vzai  = vza + i*bllen;
                real * RESTRICT vzbi  = vzb + i*bllen;
                real * RESTRICT tmpai = tmpa + i*bllen;
                real * RESTRICT tmpbi = tmpb + i*bllen;
                for(k=blstart; k<blend; k++) {
		    vzai[k]  = -(a0[k]*prefb2a[k][0]+
		                 ax[k]*prefb2a[k][1]+
                                 ay[k]*prefb2a[k][2]+
			         az[k]*prefb2a[k][3]);  
		    vzbi[k]  = -(a0[k]*prefb2b[k][0]+
		                 ax[k]*prefb2b[k][1]+
                                 ay[k]*prefb2b[k][2]+
			         az[k]*prefb2b[k][3]);  
                    vyai[k]  = -(a0[k]*prefc2a[k][0]+
		                 ax[k]*prefc2a[k][1]+
                                 ay[k]*prefc2a[k][2]+
			         az[k]*prefc2a[k][3]);  
                    vybi[k]  = -(a0[k]*prefc2b[k][0]+
		                 ax[k]*prefc2b[k][1]+
                                 ay[k]*prefc2b[k][2]+
			         az[k]*prefc2b[k][3]);  
		    tmpai[k] = a0[k]*pref3a[k][0]+
                               ax[k]*pref3a[k][1]+
                               ay[k]*pref3a[k][2]+
			       az[k]*pref3a[k][3];  
                    tmpbi[k] = a0[k]*pref3b[k][0]+
                               ax[k]*pref3b[k][1]+
                               ay[k]*pref3b[k][2]+
			       az[k]*pref3b[k][3]; 
                }
            }
	/* distribute VXC[3] contributions */ 
        for(ibl=0; ibl<ibl_cnt; ibl++) {
            for(i=iblocks[ibl][0]-1; i<iblocks[ibl][1]; i++) {
                integer ioff = i*inforb_.nbast;
		real suma;
                real sumb;
                real * RESTRICT tmpai = tmpa  + i*bllen;
                real * RESTRICT tmpbi = tmpb  + i*bllen;  
                real * RESTRICT vyai  = vya + i*bllen;
                real * RESTRICT vybi  = vyb + i*bllen;
                real * RESTRICT vzai  = vza + i*bllen;
                real * RESTRICT vzbi  = vzb + i*bllen;
                integer jsym = inforb_.muld2h[data->symYZ-1][isym]-1;
                integer (*RESTRICT jblocks)[2] = BASBLOCK(grid,jsym);
                integer jbl_cnt = grid->bas_bl_cnt[jsym];
                for(jbl=0; jbl<jbl_cnt; jbl++) {
                    for(j=jblocks[jbl][0]-1; j<jblocks[jbl][1]; j++) {
                        real * RESTRICT aosj = aos + j*bllen;
                        suma = 0; sumb = 0;
                        for(k=blstart; k<blend; k++) { 
                            suma += aosj[k]*tmpai[k];
                            sumb += aosj[k]*tmpbi[k];
                        }
                        om_a[j+ioff] += suma;
                        om_b[j+ioff] += sumb;
                    }
                }
                jsym = inforb_.muld2h[data->symZ-1][isym]-1;
                jblocks = BASBLOCK(grid,jsym);
                jbl_cnt = grid->bas_bl_cnt[jsym];
                for(jbl=0; jbl<jbl_cnt; jbl++) {
                    for(j=jblocks[jbl][0]-1; j<jblocks[jbl][1]; j++) {
                        real * RESTRICT aosj = aos + j*bllen;
                        suma = 0; sumb = 0;
                        for(k=blstart; k<blend; k++) {
                            suma += aosj[k]*vyai[k];
                            sumb += aosj[k]*vybi[k];
                        }
                        omYa[j+ioff] += suma;
                        omYb[j+ioff] += sumb;
                    }
                }
                jsym = inforb_.muld2h[data->symY-1][isym]-1;
                jblocks = BASBLOCK(grid,jsym);
                jbl_cnt = grid->bas_bl_cnt[jsym];
                for(jbl=0; jbl<jbl_cnt; jbl++) {
                    for(j=jblocks[jbl][0]-1; j<jblocks[jbl][1]; j++) {
                        real * RESTRICT aosj = aos + j*bllen;
                        suma = 0; sumb =0;
                        for(k=blstart; k<blend; k++) {
                           suma += aosj[k]*vzai[k];
                           sumb += aosj[k]*vzbi[k];
                        }
                        omZa[j+ioff] += suma;
                        omZb[j+ioff] += sumb;
                    }
		} 
            }
        }  
    } /*Loop over symmetries */   
}

Quad_Open_Data 
quad_open_init(real *cmo, real *kappaY, real *kappaZ, 
               integer *symY, integer *symZ, integer *spinY, integer *spinZ)
{
    Quad_Open_Data data; 

    /* Initialize repsonse data */
    /* FIX ME: change for dal_maloc after testing! */
    data.dmata   = calloc(2*inforb_.n2basx,sizeof(real));
    data.dmatb   = data.dmata+inforb_.n2basx;
    data.kappaY  = kappaY;
    data.kappaZ  = kappaZ;
    data.symY    = *symY;
    data.symZ    = *symZ;
    data.symYZ   = inforb_.muld2h[data.symY-1][data.symZ-1];
    data.spinY   = *spinY;
    data.spinZ   = *spinZ;
    
    /* Allocate one-index transformed densities */
    data.rhoa_y   = calloc(inforb_.n2basx,sizeof(real));
    data.rhob_y   = calloc(inforb_.n2basx,sizeof(real));
    data.rhoa_z   = calloc(inforb_.n2basx,sizeof(real));
    data.rhob_z   = calloc(inforb_.n2basx,sizeof(real));     
    /* Allocate two-index transformed densities */ 
    data.rhoa_yzzy  = calloc(inforb_.n2basx,sizeof(real));
    data.rhob_yzzy  = calloc(inforb_.n2basx,sizeof(real));
    /* Allocate one-index transfromed densities over grid batch */
    if(selected_func->is_gga()) {
        data.rhowa_y    = malloc(4*DFT_BLLEN*sizeof(real));
        data.rhowb_y    = malloc(4*DFT_BLLEN*sizeof(real));
        data.rhowa_z    = malloc(4*DFT_BLLEN*sizeof(real));
        data.rhowb_z    = malloc(4*DFT_BLLEN*sizeof(real)); 
        data.rhowa_yzzy = malloc(4*DFT_BLLEN*sizeof(real));
        data.rhowb_yzzy = malloc(4*DFT_BLLEN*sizeof(real));
        data.vya        = malloc(4*DFT_BLLEN*inforb_.nbast*sizeof(real));
        data.vyb        = malloc(4*DFT_BLLEN*inforb_.nbast*sizeof(real));
        data.vza        = malloc(4*DFT_BLLEN*inforb_.nbast*sizeof(real));
        data.vzb        = malloc(4*DFT_BLLEN*inforb_.nbast*sizeof(real));  
    } else {
        data.rhowa_y    = malloc(DFT_BLLEN*sizeof(real));
        data.rhowb_y    = malloc(DFT_BLLEN*sizeof(real));
        data.rhowa_z    = malloc(DFT_BLLEN*sizeof(real));
        data.rhowb_z    = malloc(DFT_BLLEN*sizeof(real)); 
        data.rhowa_yzzy = malloc(DFT_BLLEN*sizeof(real));
        data.rhowb_yzzy = malloc(DFT_BLLEN*sizeof(real)); 
        data.vya        = malloc(DFT_BLLEN*inforb_.nbast*sizeof(real));
        data.vyb        = malloc(DFT_BLLEN*inforb_.nbast*sizeof(real));
        data.vza        = malloc(DFT_BLLEN*inforb_.nbast*sizeof(real));
        data.vzb        = malloc(DFT_BLLEN*inforb_.nbast*sizeof(real));
    }
    /* Allocate VX[3] contr. components AO basis */
    data.res_omega_a = calloc(inforb_.n2basx,sizeof(real));
    data.res_omega_b = calloc(inforb_.n2basx,sizeof(real));
    data.res_omY_a   = calloc(inforb_.n2basx,sizeof(real));
    data.res_omY_b   = calloc(inforb_.n2basx,sizeof(real));
    data.res_omZ_a   = calloc(inforb_.n2basx,sizeof(real));
    data.res_omZ_b   = calloc(inforb_.n2basx,sizeof(real));
    /* Alocate temp. var. for integration */ 
    if(selected_func->is_gga()) {
        data.tmpa = calloc(4*DFT_BLLEN*inforb_.nbast,sizeof(real));
        data.tmpb = calloc(4*DFT_BLLEN*inforb_.nbast,sizeof(real));
    } else {
        data.tmpa = calloc(DFT_BLLEN*inforb_.nbast,sizeof(real));
        data.tmpb = calloc(DFT_BLLEN*inforb_.nbast,sizeof(real));
    }
    /* Prefactors */
    if(selected_func->is_gga()) {
        data.pref3a  = calloc(4*DFT_BLLEN,sizeof(real));
        data.pref3b  = calloc(4*DFT_BLLEN,sizeof(real));
        data.prefb2a = calloc(4*DFT_BLLEN,sizeof(real));
        data.prefb2b = calloc(4*DFT_BLLEN,sizeof(real));
        data.prefc2a = calloc(4*DFT_BLLEN,sizeof(real));
        data.prefc2b = calloc(4*DFT_BLLEN,sizeof(real)); 
    } else {
        data.pref3a  = calloc(DFT_BLLEN,sizeof(real));
        data.pref3b  = calloc(DFT_BLLEN,sizeof(real));
        data.prefb2a = calloc(DFT_BLLEN,sizeof(real));
        data.prefb2b = calloc(DFT_BLLEN,sizeof(real));
        data.prefc2a = calloc(DFT_BLLEN,sizeof(real));
        data.prefc2b = calloc(DFT_BLLEN,sizeof(real));  
    }
    return data;
}

static void  
quad_open_free(Quad_Open_Data data)
{
    /* Free memory */
    free(data.dmata);
    free(data.rhoa_y);
    free(data.rhob_y);
    free(data.rhoa_z);
    free(data.rhob_z);  
    free(data.rhoa_yzzy);
    free(data.rhob_yzzy);
    free(data.rhowa_y);
    free(data.rhowb_y);
    free(data.rhowa_z);
    free(data.rhowb_z);  
    free(data.rhowa_yzzy);
    free(data.rhowb_yzzy);
    free(data.vya); 
    free(data.vyb);
    free(data.vza); 
    free(data.vzb);
    free(data.res_omega_a);
    free(data.res_omega_b);
    free(data.res_omY_a);
    free(data.res_omY_b);
    free(data.res_omZ_a);
    free(data.res_omZ_b);
    free(data.tmpa); 
    free(data.tmpb);
    free(data.pref3a);
    free(data.pref3b);
    free(data.prefb2a);
    free(data.prefb2b);
    free(data.prefc2a);
    free(data.prefc2b);
}

static void
commute_den_x(integer den, real * RESTRICT kappaY, 
              integer symY, real * RESTRICT com_mat)
{
    integer isym, i, j;
    integer norbt = inforb_.norbt;
    integer iocc, jocc;
   
    real * res, * ky;
    /* if den 1 alpha, else beta density */ 
    for(isym=0; isym<inforb_.nsym; isym++) {
        /* the block in question corresponds to (isym,jsym) block */
        integer istart= inforb_.iorb[isym];
        integer iorb  = inforb_.norb[isym];
        integer jsym  = inforb_.muld2h[symY-1][isym]-1;
        integer jstart= inforb_.iorb[jsym];
        integer jorb  = inforb_.norb[jsym];
        if (den) {
            iocc  = inforb_.nocc[isym];
            jocc  = inforb_.nocc[jsym];            
        } else {
            iocc  = inforb_.nish[isym];
            jocc  = inforb_.nish[jsym];
	}
        res = com_mat  + istart + jstart*norbt;
        ky  = kappaY   + istart + jstart*norbt;
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
}


static __inline__ void
commute_matrices(integer sz, const real* a, const real* b,
                 real* c, integer addp)
{
    static const real MONER = -1.0;
    const real* firstpref = addp ? &ONER : &ZEROR;
    /* this could be optimized to use symmetry... */
    dgemm_("N", "N", &sz, &sz, &sz, &ONER,
           a, &sz, b, &sz, firstpref, c, &sz);
    dgemm_("N", "N", &sz, &sz, &sz, &MONER,
           b, &sz, a, &sz, &ONER, c, &sz);
}


static void 
transform_mat(real *cmo, real *matmo, integer symYZ, real *matao)
{
     static const real HALFR =  0.5;
     real * tmp;
     integer isym, norbt, nbast;
    
     norbt = inforb_.norbt;
     nbast = inforb_.nbast;
     tmp   = calloc(norbt*nbast,sizeof(real));
     for(isym=0; isym<inforb_.nsym; isym++) {
         integer ibasi = inforb_.ibas[isym];
         integer nbasi = inforb_.nbas[isym];
         integer iorbi = inforb_.iorb[isym];
         integer norbi = inforb_.norb[isym];
         integer icmoi = inforb_.icmo[isym];
         integer jsym  = inforb_.muld2h[symYZ-1][isym]-1;
         integer ibasj = inforb_.ibas[jsym];
         integer nbasj = inforb_.nbas[jsym];
         integer iorbj = inforb_.iorb[jsym];
         integer norbj = inforb_.norb[jsym];
         integer icmoj = inforb_.icmo[jsym];
         if(norbi == 0 || norbj == 0) continue;
         dgemm_("N","N", &nbasi, &norbj, &norbi,
                &HALFR, cmo + icmoi, &nbasi, matmo + iorbi + iorbj*norbt, &norbt,
                &ZEROR, tmp, &nbasi);
         dgemm_("N","T", &nbasi, &nbasj, &norbj,
                &ONER, tmp, &nbasi, cmo + icmoj, &nbasj,
                &ZEROR, matao + ibasi + ibasj*nbast, &nbast);
    }
    free(tmp);
}

/***********************************************************************/
/*             Quadratic Response for High Spin Open-Shell             */
/*                          Z. Rinkevicius                             */
/* =================================================================== */
void
FSYM2(dft_qr_ab)(real * fi, real * fo, real *cmo,
                 real *kappaY, integer *symY, integer *spinY, 
                 real *kappaZ, integer *symZ, integer *spinZ,
                 integer *addfock, real *work, integer *lwork, integer *iprint)
{ 
    static const real DP5R  =  0.5;
    static const real MDP5R = -0.5;
    static const real ONER  =  1.0;
    integer i, j, isymsav; 
    real * dv;
    struct tms starttm, endtm; clock_t utm;
    Quad_Open_Data qr_data; 
    real electrons, exp_el = 2.0*inforb_.nisht+inforb_.nasht;
    real *tmp, *tmpa, *tmpb, *kappa_yz;    
    times(&starttm);

    qr_data = quad_open_init(cmo,kappaY,kappaZ,symY,symZ,spinY,spinZ); 
   
    /* Get Density in AO basis */
    FSYM2(dft_get_ao_dens_matab)(cmo,qr_data.dmatb,qr_data.dmata,work,lwork);
    FSYM(daxpy)(&inforb_.n2basx,&DP5R,qr_data.dmatb,&ONEI,qr_data.dmata,&ONEI);
    FSYM(dscal)(&inforb_.n2basx,&DP5R,qr_data.dmatb,&ONEI);
    /* Get active density matrix in MO basis*/
    dv=dal_malloc(inforb_.n2ashx*sizeof(real));
    FSYM(dunit)(dv,&inforb_.nasht);
    /* Compute one index transformed densities */      
    isymsav = FSYM(isetksymop)(symY);
    FSYM(deq27)(cmo,kappaY,dv,qr_data.rhob_y,qr_data.rhoa_y,work,lwork);
    dscal_(&inforb_.n2basx,&DP5R,qr_data.rhob_y,&ONEI);    
    FSYM(daxpy)(&inforb_.n2basx,&ONER,qr_data.rhob_y,&ONEI,qr_data.rhoa_y,&ONEI);
    FSYM(isetksymop)(symZ);
    FSYM(deq27)(cmo,kappaZ,dv,qr_data.rhob_z,qr_data.rhoa_z,work,lwork);
    dscal_(&inforb_.n2basx,&DP5R,qr_data.rhob_z,&ONEI);    
    FSYM(daxpy)(&inforb_.n2basx,&ONER,qr_data.rhob_z,&ONEI,qr_data.rhoa_z,&ONEI);
    FSYM(isetksymop)(&isymsav);
    free(dv);
    /* Compute two-index transformed densities */
    tmp      = calloc(inforb_.n2orbx,sizeof(real));
    kappa_yz = calloc(inforb_.n2basx,sizeof(real));
    commute_den_x(1,kappaY,*symY,tmp);
    commute_matrices(inforb_.norbt,tmp,kappaZ,kappa_yz,0);
    FSYM(dzero)(tmp,&inforb_.n2orbx);
    commute_den_x(1,kappaZ,*symZ,tmp);
    commute_matrices(inforb_.norbt,tmp,kappaY,kappa_yz,1);
    transform_mat(cmo,kappa_yz,qr_data.symYZ,qr_data.rhoa_yzzy);  
    FSYM(dzero)(tmp,&inforb_.n2orbx);
    FSYM(dzero)(kappa_yz,&inforb_.n2orbx);
    commute_den_x(0,kappaY,*symY,tmp);
    commute_matrices(inforb_.norbt,tmp,kappaZ,kappa_yz,0);
    FSYM(dzero)(tmp,&inforb_.n2orbx);
    commute_den_x(0,kappaZ,*symZ,tmp);
    commute_matrices(inforb_.norbt,tmp,kappaY,kappa_yz,1);
    transform_mat(cmo,kappa_yz,qr_data.symYZ,qr_data.rhob_yzzy);  
    free(tmp);
    free(kappa_yz);
    /* Integrate VXC[3] */
    electrons = dft_integrate_ao_bl( 2, qr_data.dmata, work, lwork, iprint,  0,
                                     (DftBlockCallback)(selected_func->is_gga() ? 
                                      quad_open_gga_cb:quad_open_lda_cb),
                                     &qr_data); 
     
    /* Symmetrize integrated contributions */
    for(i=0; i<inforb_.nbast; i++) {
        for(j=0; j<i; j++) {
            integer ji = j + i*inforb_.nbast;
            integer ij = i + j*inforb_.nbast;
            real avg = 0.5*(qr_data.res_omega_a[ij]+qr_data.res_omega_a[ji]);
            qr_data.res_omega_a[ij] = qr_data.res_omega_a[ji] = avg;
            avg = 0.5*(qr_data.res_omega_b[ij]+qr_data.res_omega_b[ji]);
            qr_data.res_omega_b[ij] = qr_data.res_omega_b[ji] = avg;
            avg = 0.5*(qr_data.res_omY_a[ij]+qr_data.res_omY_a[ji]);
            qr_data.res_omY_a[ij] = qr_data.res_omY_a[ji] = avg;
            avg = 0.5*(qr_data.res_omY_b[ij]+qr_data.res_omY_b[ji]);
            qr_data.res_omY_b[ij] = qr_data.res_omY_b[ji] = avg;
            avg = 0.5*(qr_data.res_omZ_a[ij]+qr_data.res_omZ_a[ji]);
            qr_data.res_omZ_a[ij] =  qr_data.res_omZ_a[ji] = avg;
            avg = 0.5*(qr_data.res_omZ_b[ij]+qr_data.res_omZ_b[ji]);
            qr_data.res_omZ_b[ij] = qr_data.res_omZ_b[ji] = avg;
        }
    }
    /* Transform contrbutions to MO basis and add them */
    tmpa = calloc(inforb_.n2orbx,sizeof(real));
    tmpb = calloc(inforb_.n2orbx,sizeof(real));
    FSYM(lrao2mo)(cmo, &qr_data.symYZ, qr_data.res_omega_a, tmpa, work, lwork);
    FSYM(lrao2mo)(cmo, &qr_data.symYZ, qr_data.res_omega_b, tmpb, work, lwork);
    FSYM(daxpy)(&inforb_.n2orbx, &ONER, tmpa, &ONEI, fi, &ONEI);
    FSYM(daxpy)(&inforb_.n2orbx, &MDP5R, tmpa, &ONEI, fo, &ONEI);
    FSYM(daxpy)(&inforb_.n2orbx, &DP5R, tmpb, &ONEI, fo, &ONEI);
    /* [ky,Vxc(Z)] */
    FSYM(dzero)(tmpa, &inforb_.n2orbx);
    FSYM(dzero)(tmpb, &inforb_.n2orbx);
    FSYM(lrao2mo)(cmo, &qr_data.symZ, qr_data.res_omY_a, tmpa, work, lwork);
    commute_matrices(inforb_.norbt, tmpa, kappaY, tmpb, 1);
    FSYM(daxpy)(&inforb_.n2orbx, &ONER, tmpb, &ONEI, fi, &ONEI);
    FSYM(daxpy)(&inforb_.n2orbx, &MDP5R, tmpb, &ONEI, fo, &ONEI);
    FSYM(dzero)(tmpa, &inforb_.n2orbx);
    FSYM(dzero)(tmpb, &inforb_.n2orbx);
    FSYM(lrao2mo)(cmo, &qr_data.symZ, qr_data.res_omY_b, tmpa, work, lwork);
    commute_matrices(inforb_.norbt, tmpa, kappaY, tmpb, 1);
    FSYM(daxpy)(&inforb_.n2orbx, &DP5R, tmpb, &ONEI, fo, &ONEI);
    /* [kz,Vxc(Y)] */
    FSYM(dzero)(tmpa, &inforb_.n2orbx);
    FSYM(dzero)(tmpb, &inforb_.n2orbx);
    FSYM(lrao2mo)(cmo, &qr_data.symY, qr_data.res_omZ_a, tmpa, work, lwork);
    commute_matrices(inforb_.norbt, tmpa, kappaZ, tmpb, 1);
    FSYM(daxpy)(&inforb_.n2orbx, &ONER, tmpb, &ONEI, fi, &ONEI);
    FSYM(daxpy)(&inforb_.n2orbx, &MDP5R, tmpb, &ONEI, fo, &ONEI);
    FSYM(dzero)(tmpa, &inforb_.n2orbx);
    FSYM(dzero)(tmpb, &inforb_.n2orbx);
    FSYM(lrao2mo)(cmo, &qr_data.symY, qr_data.res_omZ_b, tmpa, work, lwork);
    commute_matrices(inforb_.norbt, tmpa, kappaZ, tmpb, 1);
    FSYM(daxpy)(&inforb_.n2orbx, &DP5R, tmpb, &ONEI, fo, &ONEI);
    /* clean up memory */
    free(tmpa);
    free(tmpb);
    quad_open_free(qr_data);
    /* print time and exit */
    if (*iprint>0) {
      times(&endtm);
      utm = endtm.tms_utime-starttm.tms_utime;
      fort_print("Electrons: %f(%9.3f): QR-DFT-AB evaluation time: %9.1f s",
                 electrons,exp_el, utm/(double)sysconf(_SC_CLK_TCK));
    }
}
