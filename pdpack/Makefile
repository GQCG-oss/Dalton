all :	libpdpack.a

include ../Makefile.config

C_OBJS = scatter-io.o
F_OBJS = arhpack.o linextra.o printpkg.o jacobi.o gp_zblas.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > pdpack.ftnchek

libpdpack.a: $(F_OBJS) $(C_OBJS) $(PDPACK_EXTRAS)
	 $(AR) $(ARFLAGS) $@ $?

# eispack.o has to be compiled without optimization under Linux/pgf77

eispack.o: eispack.F
	$(F77) $(SAFEFFLAGS) -c $*.F

clean :
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

depend:
	../tools/check_depend $(INCLUDES) -I../dft > Makefile.depend

include Makefile.depend
