include ../Makefile.config

all : libamfi.a

F_OBJS = 	amfi.o symtra.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

libamfi.a: $(F_OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean :	
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a

depend:
	../tools/check_depend $(INCLUDE_DIR) > Makefile.depend

include Makefile.depend
