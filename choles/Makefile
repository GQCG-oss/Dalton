all : libcholes.a

include ../Makefile.config

F_OBJS = \
	cc_choatr.o  cc_choftr.o    cc_cholrim.o  cc_chopt.o \
	cc_chocc2.o  cc_chogd.o     cc_chomp2.o   cc_chopt_dbg.o \
	cc_chodbg.o  cc_choio.o     cc_cia.o \
	cc_chofop.o  cc_cholesky.o  cc_chopden.o  cc_ciadbg.o \
	cc_cyidbg.o  cc_cyi.o       cc_rmrsp.o    cc_decmo.o  \
	cc_wrrspd.o  cc_rdrspd.o    cc_caa.o      nocc_util.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > choles.ftnchek

cc_chopt.o: cc_chopt.F
	$(F90) $(INCLUDES) $(CPPFLAGS) $(SAFEFFLAGS) -c cc_chopt.F

libcholes.a : $(F_OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

Makefile.depend depend: ../include/*.h
	../tools/check_depend $(INCLUDES) > Makefile.depend

include Makefile.depend
