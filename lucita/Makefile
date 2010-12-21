#
# FILE    : Makefile for dalton version of lucita
#           hjaaj May 2008 / Nov 2010
#

all : mod depend liblucita.a

include ../Makefile.config

F77_OBJS = ciinfo.o  density.o  diagonal.o  diagonalize.o  \
           fock.o  dalton_interface.o  io.o  \
           memory.o program.o \
           property.o  sigma.o  strings.o  transform.o  utils.o  \
           vecmat.o  par_lucia.o  par_solver.o

# F90_MODS should contain all files with F90 modules in proper order
#          (i.e. if b.F90 calls a module in k.F90, then the proper order is k.F90 b.F90)
F90_MODS = dalton_mpi.F90 lucita_mcscf_ci_cfg.F90 lucita_cfg.F90 lucita_orbitalspaces.F90 par_sync_cw.F90 par_io.F90 par_communication_model.F90 par_io_fh.F90

F90_MOD_OBJS = $(F90_MODS:.F90=.o)

F90_OBJS = 

C_OBJS =

F77_CPPS = $(F77_OBJS:.o=.i)

pre:	$(F77_CPPS)

ftnchek:pre
	ftnchek $(CHEKFLAGS) $(F77_CPPS) > ftnchek.out

liblucita.a : $(F77_OBJS) $(F90_OBJS) $(F90_MOD_OBJS) $(C_OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a *.mod

Makefile.depend depend: ../include/*.h *.inc
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
