all : libsirius.a

include ../Makefile.config

F_OBJS = \
	ciauth.o          ciham.o           cioff.o\
	cisig.o           csfstuff.o        dens.o\
	detfi.o           detpr.o           dets0.o\
	dets1.o           drcctl.o          getcix.o\
	getint.o          siraba.o          sirave.o\
	sircan.o          sirci.o\
	sircmo.o          sirctl.o\
	sirdiis.o         sirfck.o          sirgp.o\
	sirgrad.o         sirinp.o          sirlintrn.o\
	sirmp2.o          sirneo.o          sirnexkap.o\
	sirnr.o           siropt.o          sirorb.o\
	sirout.o          sirpop.o          sirsav.o\
	sirset.o          sirtest.o         sirtr1.o\
	sirtra.o          sirupd.o\
	sirsol.o          symchk.o\
	newtrahso.o       sirntra.o         intsort.o \
	r12aux.o          sirloc.o\
	koopro4.o         dalcip.o          dypc.o \
        sirpcm.o          sirsolpcm.o       sirpcmder.o \
        sirief.o          sircav.o          sirqm3.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > sirius.ftnchek

calltree: pre
	ftnchek $(CHEKFLAGS) -calltree=tree -nocheck -mkhtml *.i > sirius.calltree

libsirius.a : $(F_OBJS)
	$(AR) $(ARFLAGS) $@ $?

# sirtra.o has to be compiled without less optimization under HP-UX

#sirtra.o: sirtra.F
#	$(F77) $(INCLUDE_DIR) $(CPP_FLAGS) $(SAFEOPTFLAGS) -c $*.F

clean : 
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

Makefile.depend depend: ../include/*.h
	../tools/check_depend $(INCLUDES) > Makefile.depend

include Makefile.depend
