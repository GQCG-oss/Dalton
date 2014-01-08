all: libgeomopt.a

include ../Makefile.config 

INCLUDES := $(INCLUDES) -I../lsutil -I../LSint -I./dqdx

F90_OBJS = \
        ls_opt.o \
        ls_opt2.o \
        ls_redint.o \
        LSopt-input.o \
        q_to_x.o \

F_CPPS = $(F_OBJS:.o=.i)

CXXFILES =\
	dqdx.o

OBJS = $(F90_OBJS:.F90=.o) $(CXXFILES:.cpp=.o)

pre: $(F_CPPS)

libgeomopt.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD


Makefile.depend depend: $(F90_OBJS:.o=.F90)
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -M../lsutil -M../LSint -- $(F90_OBJS:.o=.F90) > Makefile.depend

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
