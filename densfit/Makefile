all: libdensfit.a

include ../Makefile.config

F_OBJS =  df_dirac.o dposv.o denfit_noddy.o #df_dofit.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

ftnchek: pre
	        ftnchek $(CHEKFLAGS) *.i > densfit.ftnchek

libdensfit.a: $(F_OBJS)
	        $(AR) $(ARFLAGS) $@ $?

clean :
	        $(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

Makefile.depend depend: ../include/*.h
	        ../tools/check_depend $(INCLUDES) > Makefile.depend

include Makefile.depend
