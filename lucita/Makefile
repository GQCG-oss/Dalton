#
# FILE    : Makefile for dalton version of lucita
#           hjaaj May 2008
#

all : liblucita.a

include ../Makefile.config

F77_OBJS = ciinfo.o  density.o  diagonal.o  diagonalize.o  \
           fock.o  dalton_interface.o  io.o  \
           lucia_cc.o  lucia_icci.o  memory.o program.o \
           property.o  sigma.o  strings.o  transform.o  utils.o  \
           vecmat.o  par_lucia.o  par_solver.o \
           gplucia.o gplupar.o
# files gp/gplucia.F gp/gplupar.F include/parluci.h include/luctinp.h
# and include/mxmcgas.h in DIRAC are all moved to the lucita/
# in dalton (for the time being). /hjaaj May 2008

F90_OBJS =

F90_MODS =

F77_CPPS = $(F77_OBJS:.o=.i)

pre:	$(F77_CPPS)

ftnchek:pre
	ftnchek $(CHEKFLAGS) $(F77_CPPS) > ftnchek.out

liblucita.a : $(F77_OBJS)
	$(AR) $(ARFLAGS) $@ $?


clean : 
	$(RM) -f *.o *~ *.f *.f90 *.lst *.i

mod:
	$(MAKE) f90dep  # create list of which file to compile to get the .mod files
	$(MAKE) module_files

module_files: $(F90_MODS)

depend:
	makedepend -- $(DEPEND_FLAGS) -- $(F77_OBJS:.o=.F) $(F90_OBJS:.o=.F90)

f90dep: ;\
       ../utils/sfmakedepend.pl -I $(SOURCEDIR)/include -M $(SOURCEDIR)/modules $(F90_OBJS:.o=.F90) $(F90_OBJS:.o=.F) $(F77_OBJS:.o=.F)
 

# DO NOT DELETE THIS LINE -- make depend depends on it.
