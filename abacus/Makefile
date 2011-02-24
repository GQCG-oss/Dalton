all: libabacus.a

include ../Makefile.config

F_OBJS = 	dalgnr.o     aba2tex.o    abaaat.o     abadip.o\
	abadrv.o     abadtra.o    abaedrv.o    abaenr.o\
	abaesym.o    abaexc.o     abaflt.o     abafsd.o\
	abalnr.o     abalrs.o     abaltra.o    abamag.o\
	abanuc.o     abaopg.o     abaort.o     abaptr.o\
	abarhs.o     abarlx.o     abarspn.o    her2her.o\
	abaset.o     abatro.o     abavib.o     abawalk.o\
	her1ave.o    her1car.o    her1drv.o    her1int.o\
	her1odc.o    her1out.o    her1pro.o    her1sym.o\
	her2ave.o    her2dir.o    her2drv.o    her2el1.o\
	her2el2.o    her2fck.o    her2her.o    her2odc.o\
	her2out.o    her2sup.o    her2sym.o    herbas.o\
	herdirect.o  herdkh.o     herdrv.o     hergam.o     hergeo.o\
	hergp.o      hergroup.o   herrdn.o     hersol.o     huckel.o\
	her2gab.o    hersym.o     hertst.o     herqm3.o\
	abaopt.o     abavrml.o    abaqr.o      abarint.o\
	abander.o    abandsym.o   herfckcon.o  abaop2.o\
	abasop.o     abamolden.o  eckart.o     ecp.o\
	abafmm.o     abanvib.o    aba2r12.o    aba2eth.o\
	abactocd.o   abaqpg.o     abapcm.o     abaqm3.o\
	abaqmmm.o    gen1int_block.o    gen1int_sym.o    gen1int_wrapper.o

F_CPPS = $(F_OBJS:.o=.i) dalton.i 

pre: $(F_CPPS)

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > abacus.ftnchek

libabacus.a: $(F_OBJS)
	$(AR) $(ARFLAGS) $@ $?

main  : dalton.o

slave : herpar.o

linuxmem : linux_mem_allo.o

clean :	
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

Makefile.depend depend: ../include/*.h
	../tools/check_depend $(INCLUDES) > Makefile.depend

include Makefile.depend

her2el2.o: her2el2.F
	$(F90) $(INCLUDES) $(CPPFLAGS) $(SAFE2_FFLAGS) -c $*.F 

hersym.o: hersym.F
	$(F90) $(INCLUDES) $(CPPFLAGS) $(SAFE2_FFLAGS) -c $*.F 
