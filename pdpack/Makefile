include ../Makefile.config

all :	libpdpack.a

OBJS = arhpack.o fastio_g07.o linextra.o printpkg.o jacobi.o

libpdpack.a: $(OBJS) $(PDPACK_EXTRAS)
	 $(AR) $(ARFLAGS) $@ $?

# eispack.o has to be compiled without optimization under Linux/pgf77

eispack.o: eispack.F
	$(F77) $(INCLUDE_DIR) $(CPP_FLAGS) $(SAFEOPTFLAGS) -c $?

clean :
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a

depend:; ../tools/check_depend $(INCLUDE_DIR) > Makefile.depend

######################################################################
# read dependencies:
include Makefile.depend
###################################################################### 
