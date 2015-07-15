all: librsp_properties.a

include ../Makefile.config

INCLUDES := $(INCLUDES) -I../lsutil -I../LSint -I../linears -I../responsesolver -I../SolverUtilities

# Keep the file lists sorted for easier modification.
# PLEASE KEEP THEM SORTED ALPHABETICALLY!!!

F90_OBJS = \
	molecular_hessian.o\
	nuclei_selected_shielding.o\
	response_prop_noOpenRSP.o\
	test_molHessian.o

# we should really have a separate target for bsm to properly
# track dependencies but we would not like it to be too troublesome.
#$(C_OBJS) 
librsp_properties.a: $(F90_OBJS) 
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD


Makefile.depend depend: $(F90_OBJS:.o=.F90)
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -M../lsutil -M../LSint -M../linears -M../responsesolver -M../SolverUtilities -- $(F90_OBJS:.o=.F90) > Makefile.depend

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
