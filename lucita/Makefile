#
# FILE    : Makefile for dalton version of lucita
#           hjaaj May 2008 / Nov 2010
#

all : liblucita.a depend

include ../Makefile.config

F77_OBJS = ciinfo.o  density.o  diagonal.o  diagonalize.o  \
           fock.o  dalton_interface.o  io.o  \
           memory.o program.o \
           property.o  sigma.o  strings.o  transform.o  utils.o  \
           vecmat.o  par_lucia.o  par_solver.o lucidal_input.o

F90_OBJS = par_io.o

F90_MODS = par_io.o

F77_CPPS = $(F77_OBJS:.o=.i)

pre:	$(F77_CPPS)

ftnchek:pre
	ftnchek $(CHEKFLAGS) $(F77_CPPS) > ftnchek.out

liblucita.a : $(F77_OBJS) $(F90_OBJS)
	$(AR) $(ARFLAGS) $@ $?


clean : 
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a *.mod

Makefile.depend depend: ../include/*.h *.inc
	../tools/check_depend $(INCLUDES) > Makefile.depend

include Makefile.depend

mod:
	$(MAKE) f90dep  # create list of which file to compile to get the .mod files
	$(MAKE) module_files

module_files: $(F90_MODS)

f77depend:
	makedepend -- $(DEPEND_FLAGS) -- $(F77_OBJS:.o=.F) $(F90_OBJS:.o=.F90)

f90dep: ;\
       ../utils/sfmakedepend.pl -I $(SOURCEDIR)/include -M $(SOURCEDIR)/modules $(F90_OBJS:.o=.F90) $(F90_OBJS:.o=.F) $(F77_OBJS:.o=.F)
 
