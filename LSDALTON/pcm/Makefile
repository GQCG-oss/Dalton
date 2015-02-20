all: liblspcm.a

include ../Makefile.config

INCLUDES := $(INCLUDES) -I../lsutil -I../LSint -I../linears

F90_OBJS = ls_pcm_config.o\
           ls_pcm_utils.o\
	   ls_pcm_write.o\
	   ls_pcm_integrals.o\
	   ls_pcm_scf.o

#C_OBJS = memory.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

#$(C_OBJS) 
liblspcm.a: $(F_OBJS) $(F90_OBJS) 
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD

Makefile.depend depend: $(F90_OBJS:.o=.F90)
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -M../lsutil -M../LSint -M../linears -- $(F_OBJS:.o=.F) $(F90_OBJS:.o=.F90) > Makefile.depend

######################################################################
# read dependencies using sinclude (silent include)
# which continues if Makefile.depend is not readable (incl. doesn't exist!).
# Makefile.depend can be created with 'make depend'.
include Makefile.depend
###################################################################### 

# Gao: Fortran compiler creates a file .mod for each module name.
# This confuses the gnumake, which thinks that .mod files are
# Modula2 source code files. We must override this built-in rule
# by adding the following line:
%.o: %.mod
