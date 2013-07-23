all: libsolverutilities.a

include ../Makefile.config

INCLUDES := $(INCLUDES) -I../lsutil -I../LSint

# Keep the file lists sorted for easier modification.
# PLEASE KEEP THEM SORTED ALPHABETICALLY!!!

F90_OBJS = \
        arh_density.o \
	dalton_interface.o \
        dd_utilities.o \
        dd_utilities_unres.o \
        decomp.o        \
	lsdalton_modules.o \
	queue-operations.o \
        queue-typedef.o \
        rsp_precond.o \
        rsp_utilities.o 

libsolverutilities.a: $(F90_OBJS) 
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD


Makefile.depend depend: $(F90_OBJS:.o=.F90)
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -M../lsutil -M../LSint -- $(F_OBJS:.o=.F) $(F90_OBJS:.o=.F90) > Makefile.depend

include Makefile.depend
######################################################################
# read dependencies using sinclude (silent include)
# which continues if Makefile.depend is not readable (incl. doesn't exist!).
# Makefile.depend can be created with 'make depend'.
#-include Makefile.depend
###################################################################### 

# Gao: Fortran compiler creates a file .mod for each module name.
# This confuses the gnumake, which thinks that .mod files are
# Modula2 source code files. We must override this built-in rule
# by adding the following line:
%.o: %.mod
