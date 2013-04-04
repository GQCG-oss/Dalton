all :	libpdpack.a

include ../Makefile.config

C_OBJS = scatter-io.o
F_OBJS = jacobi.o gp_zblas.o gp_qblas.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > pdpack.ftnchek

libpdpack.a: $(F_OBJS) $(C_OBJS) $(PDPACK_EXTRAS)
	 $(AR) $(ARFLAGS) $@ $?

# eispack.o has to be compiled without optimization under Linux/pgf77

eispack.o: eispack.F
	$(F90) $(SAFEFFLAGS) -c $*.F

clean :
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

Makefile.depend depend:  $(F_OBJS:.o=.F)
	../tools/check_depend $(INCLUDES) -I../dft > Makefile.depend

include Makefile.depend
