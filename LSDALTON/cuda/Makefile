all: libcuda.a

include ../Makefile.config 

INCLUDES := $(INCLUDES) -I../lsutil 

# Keep the file lists sorted for easier modification.

F90_OBJS = gpu_interfaces.o

libcuda.a: $(F90_OBJS) 
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD

Makefile.depend depend: $(F90_OBJS:.o=.F90) 
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -M../lsutil -- $(F90_OBJS:.o=.F90) > Makefile.depend

include Makefile.depend
######################################################################
# read dependencies using sinclude (silent include)
# which continues if Makefile.depend is not readable (incl. doesn't exist!).
# Makefile.depend can be created with 'make depend'.
#-include Makefile.depend
###################################################################### 

%.o: %.mod
