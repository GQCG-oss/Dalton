include ../Makefile.config

all : libamfi.a

OBJS = 	amfi.o symtra.o

libamfi.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean :	
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a

depend:; ../tools/check_depend $(INCLUDE_DIR) > Makefile.depend

######################################################################
# read dependencies:
include Makefile.depend
###################################################################### 
