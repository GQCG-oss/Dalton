all :	libpdpack.a

include ../Makefile.config

F_OBJS = arhpack.o fastio_g07.o linextra.o printpkg.o jacobi.o gp_zblas.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > pdpack.ftnchek

libpdpack.a: $(F_OBJS) $(PDPACK_EXTRAS)
	 $(AR) $(ARFLAGS) $@ $?

# eispack.o has to be compiled without optimization under Linux/pgf77

eispack.o: eispack.F
	$(F77) $(SAFEFFLAGS) -c $*.F

clean :
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

depend:
	../tools/check_depend $(INCLUDES) > Makefile.depend

include Makefile.depend
