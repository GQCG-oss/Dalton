all: libresponsesolver.a

include ../Makefile.config

INCLUDES := $(INCLUDES) -I../lsutil -I../LSint -I../SolverUtilities

# Keep the file lists sorted for easier modification.
# PLEASE KEEP THEM SORTED ALPHABETICALLY!!!

F90_OBJS = \
        rsp_complex.o \
        rsp_cmplx_sym.o\
        rsp_cmplx_sym_new.o\
        rsp_solver.o    \
        rsp_std_sym.o  

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

libresponsesolver.a: $(F90_OBJS) 
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD


Makefile.depend depend: $(F90_OBJS:.o=.F90)
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -M../lsutil -M../LSint -M../SolverUtilities -- $(F90_OBJS:.o=.F90) > Makefile.depend

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
