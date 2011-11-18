all : mod depend libgp.a

include ../Makefile.config

F_OBJS = gphjj.o    gpjeppe.o  gptrygve.o     mempkg.o \
	qpack.o    graphic.o  vcraypack.o    gpsaue.o \
	mpimacro.o pvmmacro.o gpdahle.o      intf7790.o \
	gpjkp.o gpjt.o gplucia.o gplucia_io.o gplupar.o \
	gpsimen.o $(GP_EXTRAS)

# F90_MODS should contain all files with F90 modules in proper order
# #          (i.e. if b.F90 calls a module in k.F90, then the proper order is k.F90 b.F90)
F90_MODS = character_processing.F90 keyword.F90 dalton_input_processing.F90 input_reader.F90
#
F90_MOD_OBJS = $(F90_MODS:.F90=.o)
#
F90_OBJS =

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS) mpi_dummy.i

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > gp.ftnchek

all: libgp.a mpi_dummy.o mpi_dummyc.o

libgp.a: $(F_OBJS) $(F90_OBJS) $(F90_MOD_OBJS) $(C_OBJS)
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
# # $(MAKE) Makefile.modules  # create list of which file to compile to get the .mod files
# # $(MAKE) module_files
# #Makefile.modules: $(F90_MODS)
# #     ../tools/f90dep.pl $(CPPFLAGS) -I $(SOURCEDIR)/include -M $(SOURCEDIR)/modules $(F_OBJS:.o=.F) $(F90_OBJS:.o=.F) $(F90_OBJS:.o=.F90) > Makefile.modules
# #############
