all : libsoppa.a

include ../Makefile.config

F_OBJS =  \
ccsd_tcmepkx.o  rp_lrsolv.o   rp_trial3.o     so_eres.o       rp_rspleq.o   \
so_lrsolv.o     so_read.o     so_residual.o   so_sigda.o      so_trial2.o   \
dc_calc.o       hr_test.o     rp_newtrial.o   rp_trmom.o      so_excit1.o   \
so_memmax.o     so_readset.o  so_res_o.o      so_sort.o       so_trial3.o   \
dc_eres.o       hr_tr1new.o   rp_optvec.o     so_anal.o       so_excout.o   \
so_moenergy.o   so_reave.o    so_res_ovl.o    so_sqait2.o     so_trmom.o    \
dc_r1vec.o      hr_trial1.o   rp_orth_trn.o   so_aodens.o     so_firgp.o    \
so_newtrial.o   so_redev.o    so_res_ovlr.o   so_sres.o       so_twofock.o  \
dc_res_o.o      hr_trial2.o   rp_prop.o       so_bcktr.o      so_fndpos.o   \
so_onefock.o    so_redgp.o    so_res_sym.o    so_t2m1.o       so_write.o    \
dc_rsplex.o     hr_trmom.o    rp_redev.o      so_bextract.o   so_getgp.o    \
so_onefocko.o   so_redle.o    so_res_tcb.o    so_t2mp.o       so_writeset.o \
dc_sres.o       rp_anal.o     rp_redgp.o      so_check.o      so_getmo.o    \
so_onep.o       so_res_a.o    so_res_tcbt.o   so_tbextract.o  so_wrtve.o    \
hr_check.o      rp_charge.o   rp_redle.o      so_close.o      so_gett2.o    \
so_onepmo.o     so_res_b.o    so_ropt.o       so_test1.o      hr_diag.o     \
rp_check.o      rp_redlej.o   so_collect_tm.o so_incred.o     so_open.o     \
so_res_cb.o     so_roptj.o    so_test.o       hr_eres.o       rp_diag.o     \
rp_residual.o   so_densai1.o  so_init.o       so_optvec.o     so_res_cbt.o  \
so_rpprp1.o     so_tfset.o    hr_getgp.o      rp_eres.o       rp_res_o.o    \
so_densai2.o    so_initial.o  so_ordeig.o     so_res_cd.o     so_rspleq.o   \
so_time.o       hr_lrsolv.o   rp_getgp.o      rp_rsplex.o     so_dens.o     \
so_intrp.o      so_orth_trn.o so_res_cdt.o    so_rsplex.o     so_tmltr.o    \
hr_newtrial.o   rp_incred.o   rp_sres.o       so_diag.o       so_lnraba.o   \
so_polar.o      so_res_c.o    so_sdiag1.o     so_tr1new.o     hr_orth_trn.o \
rp_init.o       rp_trial1.o   so_ediag1.o     so_lnrout.o     so_prop.o     \
so_res_cp.o     so_secgp.o    so_tr1rp.o      hr_rsplex.o     rp_lrsoeq.o   \
rp_trial2.o     so_ediag2.o   so_lrsoeq.o     so_prpint.o     so_res_fck.o  \
so_sigai.o      so_trial1.o   rp_residualj.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > soppa.ftnchek

libsoppa.a: $(F_OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean :
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

Makefile.depend depend: ../include/*.h
	../tools/check_depend $(INCLUDES) > Makefile.depend

include Makefile.depend




