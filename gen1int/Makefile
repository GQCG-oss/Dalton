#
# FILE    : Makefile for dalton version of Gen1Int interface
#           Oct. 2, 2011, by Bin Gao
#

all : mod depend libgen1int.a

include ../Makefile.config

# F90_MODS should contain all files with F90 modules in proper order
#          (i.e. if b.F90 calls a module in k.F90, then the proper order is k.F90 b.F90)
F90_MODS = gen1int_shell.F90

F90_MOD_OBJS = $(F90_MODS:.F90=.o)

F90_OBJS = gen1int_dal_init.o gen1int_dal_clean.o gen1int_dal_prop.o \
           gen1int_dal_main.o gen1int_dal_test.o

libgen1int.a : $(F90_OBJS) $(F90_MOD_OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a *.mod

Makefile.depend depend: ../include/*.h
	../tools/check_depend $(INCLUDES) > Makefile.depend

include Makefile.depend

mod:
	$(MAKE) --no-print-directory $(F90_MOD_OBJS)

#############
# the next may possibly be modified to automatically create the F90_MODS list ???? /SK+HJAaJ Nov 2010
# $(MAKE) Makefile.modules  # create list of which file to compile to get the .mod files
# $(MAKE) module_files
#Makefile.modules: $(F90_MODS)
#	../tools/f90dep.pl $(CPPFLAGS) -I $(SOURCEDIR)/include -M $(SOURCEDIR)/modules $(F_OBJS:.o=.F) $(F90_OBJS:.o=.F) $(F90_OBJS:.o=.F90) > Makefile.modules
#############

# Fortran compiler creates a file .mod for each module name.
# This confuses the gnumake, which thinks that .mod files are
# Modula2 source code files. We must override this built-in rule
# by adding the following line:
%.o: %.mod
