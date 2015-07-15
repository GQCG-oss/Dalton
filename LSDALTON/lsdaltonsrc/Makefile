all: liblsdaltonsrc.a

include ../Makefile.config

INCLUDES := $(INCLUDES) -I../lsutil -I../LSint -I../linears -I../geomopt -I../ddynam -I../deccc -I../pbc2 -I../rsp_properties -I../SolverUtilities -I../xcfun_host -I../IchorIntegrals -I../responsesolver -I../cuda

# Keep the file lists sorted for easier modification.
# PLEASE KEEP THEM SORTED ALPHABETICALLY!!!

F90_OBJS = \
	ADMMbasisOpt.o\
	numerical_derivatives.o\
        LS_optimizer.o \
	LSdynamics.o \
	Energy_and_deriv.o \
	lsmpiSlave.o \
	init_lsdalton.o \
	configuration.o \
	LSlibState.o \
	LSlib.o \
	Profile.o\
	IchorTesting.o\
	IchorProfile.o\
	InteractionEnergy.o\
	LSDALTON.o
SMALLF90_OBJS = \
        lsdalton_wrapper.o


FOBJS = \
	$(SMALLF90_OBJS) \
	$(F90_OBJS)
     

FORTFILES = \
        $(SMALLF90_OBJS:.o=.f90) \
	$(F90_OBJS:.o=.F90)

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

# we should really have a separate target for bsm to properly
# track dependencies but we would not like it to be too troublesome.
#$(C_OBJS) 
liblsdaltonsrc.a: $(FOBJS) 
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD


Makefile.depend depend: $(FORTFILES)
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -M../lsutil -M../LSint -M../linears -M../geomopt -M../deccc -M../SolverUtilities -M../xcfun_host -M../IchorIntegrals -M../responsesolver -M../cuda -- $(FORTFILES) > Makefile.depend


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
