all: liblinears.a

include ../Makefile.config

ifdef VAR_RSP
INCLUDES := $(INCLUDES) -I../lsutil -I../LSint -I../geomopt -I../ddynam -I../deccc -I../pbc2 -I../SolverUtilities -I../responsesolver \
	-I../pcm \
	-I../../external/matrix-defop/src -I../../external/openrsp/src
else
INCLUDES := $(INCLUDES) -I../lsutil -I../LSint -I../geomopt -I../ddynam -I../deccc -I../pbc2 -I../SolverUtilities -I../responsesolver \
	-I../pcm
endif
# Keep the file lists sorted for easier modification.
# PLEASE KEEP THEM SORTED ALPHABETICALLY!!!

#ajt Birefring compiles extremely slowly with pgf, so I disable compiling
#ajt it until I have gotten around to simplifying the code.
#ajt Calls to it from linscf.F90 have been commented out.
#ajt0410 Temporarily added again for Radovan's merging with other versions

F90_OBJS = \
        arh_debug.o \
        arh_driver.o \
        average_util.o \
	LocTypes.o \
	ChargePrec.o\
	ChargeLoc.o \
	configurationType.o \
	davidson_settings.o \
	davidson_solver.o \
	debug.o		\
	densopt.o	\
	dens_subspace.o \
	diag.o		\
	diis.o		\
	dsm.o		\
	dsm_xterm.o     \
        ecdata_module.o \
        extra-output.o \
	fock-eval.o	\
	kurtosis.o \
        leastchange.o \
        levelshift.o \
	SCFLOOP.o\
	trustradius.o\
        LSDALTON_RESPONSE.o        \
        LSDALTON_RESPONSE_type.o        \
	minimize.o \
	optimlocModf.o \
	orbspread_hess_prec.o \
	localityMeasure.o \
	OrbLoc_input.o  \
	localization_charge.o \
	localization_orbspread.o \
	orbspread_util.o \
	localization_util.o \
        response_driver.o \
        response_prop.o \
        scfopt-typedef.o \
	soeo-loop.o \
	soeo-util.o \
	soeo-matop.o \
	soeo-redspace.o \
	soeo-typedef.o \
	soeo-debug.o \
	starting_guess.o \
	statistics.o   	\
	plt_driver.o  \
        trilevel.o 
#	cgop.o
#C_OBJS = memory.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

# we should really have a separate target for bsm to properly
# track dependencies but we would not like it to be too troublesome.
#$(C_OBJS) 
liblinears.a: $(F_OBJS) $(F90_OBJS) 
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD

Makefile.depend depend: $(F_OBJS:.o=.F) $(F90_OBJS:.o=.F90)
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -M../lsutil -M../LSint -M../geomopt -M../deccc -M../SolverUtilities -M../responsesolver -M../pcm -- $(F_OBJS:.o=.F) $(F90_OBJS:.o=.F90) > Makefile.depend

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
