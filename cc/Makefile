all : libcc.a

include ../Makefile.config

F_OBJS = herdirect_cc.o  cc_drv.o       ccsd_input.o   ccsd_energy.o \
ccsd_sortao.o  ccsd_rhs.o     cc3_t3int.o    cc3_omeg.o    \
cc_exci.o      cc_lr.o        cc_fop.o       cc_pqi.o      \
cc_exgr.o      cc_rhtr.o      cc_lhtr.o      cceq_solv.o   \
cc_fmat.o      cc_landf.o     cc_landf2.o    cc_landf3.o   \
cc_den.o       cc_d1ao.o      ccsd_diis.o    cc_grad2e.o   \
cc_grad0.o  \
cc_solrsp.o    cc_getgd.o     cc_rdrsp.o     cc_wrrsp.o    \
ccrspsym.o     cclr_util.o    cc_util.o      ccmcd_setup.o \
cc_hyppol.o    cc_2hyp.o      cc_3hyp.o      cc_gmat.o     \
cc_bmat.o      cc_famat.o      ccqr_setup.o   cc_fropro.o   \
cccr_setup.o   cc_mofcon2.o   cc_cphf.o      cc_getgp.o    \
cccr_aa.o      cccr_util.o    cc_etadrv.o    ccexpfck.o    \
cc_cmat.o      cc_fdc.o       cc_qr2rs.o     cc_dmat.o     \
cc_fdd.o       ccrhsvec.o     cc_bamat.o     cc_ftst.o     \
cc_aamat.o     cc_soldrv.o    cc_batst.o     cc_fdba.o     \
cc_4hyp.o      cc_smcal.o     cc_cauchy.o    cc_bfden.o    \
cc_hmat.o      cc_tmcal.o     ccexlr.o                     \
cc_aofock2.o   crayio2.o      ccder1.o       cc_pckutil.o  \
cc_iajb.o      cc_xieta.o     cc_xetst.o     \
cc_fmattst.o   cc_grad2.o     cc_rdaod.o                   \
ccrlxxim.o     cckappasq.o    cc_1idxtran.o  cc_effckmo.o  \
cc_trdrv.o     cc_rhtr3.o     ccsd_rou3.o    ccsd_rhs3.o   \
cc_setdorps.o  cc_sop.o       cc_setfa.o     cc_setxe.o    \
ccfndop.o      cc_get_rmat.o  cc_qmat.o      cc_mmommo.o   \
cc_nuccon.o    cc_f1inact.o   cc_lambdaq.o   iroper2.o     \
cc_gradient.o  cc_fdaa.o      cc_qaos.o      cc_setfb12.o  \
cc_maomao.o    cc_dhfao.o     cc_1idx_d2eff.o cclr_setup.o \
ccximcon.o     ccreo2con.o    cc_get_lambda1.o cc_grad.o   \
cc2_refop.o    cctrbt1.o      cc_21i2.o                    \
cc_fdfb.o      cc_fbtst.o     cc_etadrv1.o                 \
cc_fbta.o      ccfbint1.o     ccfbint2.o     ccfbtaao1.o   \
ccfbtaao2.o    ccfbtaf.o      cc_lamtra.o    cc_aodens3.o  \
cc_ijcb.o      cc_bfif1.o                                  \
cc_21i3.o      cc_fckrlx1.o   cc_hpqta.o     cc_e1aim.o    \
cc_e21con.o    cc_e22con.o    cclt_yps1.o    cc_sort4o.o   \
cc_sortz2.o    cc_zkjam.o     cc_int4o.o     \
ccrhs_ipm1.o   cctrbt2.o      cc_hfr1den.o   \
cc_gim1.o      cc_bfbsort1.o  cc_int3o2.o    \
cc_mcdny.o     ccsd_triple.o  \
cc_slv.o       cc3_omeg3.o    cc3_t3int3.o   ccdiis_solv.o \
cc_foppt.o     cc3_xisd.o \
cc3_lhtr.o     ccetacor.o     cc_den2_pt.o \
cc_den_pt2.o   ccpt_symmback.o     ccpt_resort_den.o \
ccpt_btrden.o  \
ccpt_den.o     ccpt_eta.o     ccpt_fropro.o \
ccpt_grad0.o  \
cc_lhtr_noddy.o ccsdt_xi_noddy.o \
ccsdt_f_noddy.o ccsdt_ints_noddy.o \
ccsdt_eta_noddy.o ccsdt_b_noddy.o ccsdt_ba_noddy.o \
ccsdt_eta_cont.o cc3_fmat.o     ccsdt_fbmat.o  \
cc3_etasd.o cc3_ft3b.o cc3_fmatsd.o cc3_xi_cont.o \
cc3_init.o ccexnorm.o cc_opa.o ccsdt_aden_noddy.o cc3_aden.o \
cc_foppt_sc.o cc_r12.o  cc_mkrv.o cc_rhtr_noddy.o \
cc_2eexp_1.o  cc_grad_1.o   cc_2eexp_2.o   asqr.o \
cc3_gmat.o cc3_bfmat.o  ccsoppa.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > cc.ftnchek
  
# cc_pckcr8.o    cc_pckutil.o   
io : crayio.o

libcc.a: $(F_OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

tags:;          $(RM) -f tags ; etags  *.F

depend: 
	../tools/check_depend $(INCLUDES) > Makefile.depend

include Makefile.depend
