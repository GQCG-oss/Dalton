include ../Makefile.config

all : libamfi.a

F_OBJS = 	amfi.o symtra.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

libamfi.a: $(F_OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean :	
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a

depend:; ../tools/check_depend $(INCLUDE_DIR) > Makefile.depend

######################################################################
# read dependencies using sinclude (silent include)
# which continues if Makefile.depend is not readable (incl. doesn't exist!).
# Makefile.depend can be created with 'make depend'.
-include Makefile.depend
######################################################################
