
include ../Makefile.config

all : libdensfit.a

F_OBJS =  df_dirac.o df_dofit.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

ftnchek: pre
	        ftnchek $(CHEKFLAGS) *.i > densfit.ftnchek

libdensfit.a: $(F_OBJS)
	        $(AR) $(ARFLAGS) $@ $?

clean :
	        $(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

depend:
	        ../tools/check_depend $(INCLUDES) > Makefile.depend

