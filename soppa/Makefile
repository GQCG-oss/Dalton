all : libsoppa.a

include ../Makefile.config

F_OBJS =  \
so_excit1.o   so_excout.o   so_memmax.o   so_time.o     so_moenergy.o \
rp_rsplex.o   rp_trmom.o    rp_charge.o   rp_init.o     rp_lrsolv.o   \
rp_res_o.o    rp_anal.o     rp_getgp.o    rp_diag.o     rp_trial1.o   \
rp_trial2.o   rp_eres.o     rp_sres.o     rp_redev.o    rp_check.o    \
rp_orth_trn.o rp_residual.o rp_newtrial.o rp_incred.o   rp_optvec.o   \
dc_rsplex.o   dc_calc.o     dc_eres.o     dc_sres.o     dc_r1vec.o    \
dc_res_o.o    \
so_prpint.o   so_init.o     so_gett2.o    so_rsplex.o   so_trmom.o    \
so_onepmo.o   so_dens.o     so_t2mp.o     so_lrsolv.o   so_res_o.o    \
so_anal.o     so_getgp.o    so_getmo.o    so_onep.o     so_diag.o     \
so_trial1.o   so_trial2.o   so_eres.o     so_sres.o     so_redev.o    \
so_check.o    so_res_ovl.o  so_res_ovlr.o so_aodens.o   so_bcktr.o    \
so_twofock.o  so_onefock.o  so_secgp.o    so_firgp.o    so_intrp.o    \
so_ediag1.o   so_ediag2.o   so_sdiag1.o   so_orth_trn.o so_tr1new.o   \
so_tr1rp.o    so_residual.o so_tmltr.o    so_tmltr.o    so_residual.o \
so_newtrial.o so_t2m1.o     so_res_a.o    so_res_tcb.o  so_res_cb.o   \
so_sigai.o    so_sigda.o    so_res_b.o    so_res_c.o    so_res_sym.o  \
so_res_fck.o  so_res_cd.o   so_densai1.o  so_densai2.o  so_incred.o   \
so_ordeig.o   so_ropt.o     so_optvec.o   so_sort.o     so_test.o     \
so_test1.o    so_rpprp1.o   so_bextract.o so_sqait2.o   so_res_cp.o   \
so_initial.o  so_tfset.o    so_collect_tm.o \
so_open.o     so_close.o    so_read.o     so_write.o    so_wrtve.o  \
so_fndpos.o   so_readset.o  so_writeset.o \
ccsd_tcmepkx.o


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




