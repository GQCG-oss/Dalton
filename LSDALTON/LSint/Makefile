all: libLSint.a

include ../Makefile.config

INCLUDES := $(INCLUDES) -I../lsutil -I../mm -I../xcfun_host -I../dft -I../IchorIntegrals

# Keep the file lists sorted for easier modification.

F90_OBJS = \
	dft_gridLL.o\
	BuildMolFile.o\
        BuildBasis.o   \
	dalton-input.o\
	gridgeneration.o\
	gridgeneration_boxify.o\
	HODItest.o\
	MBIE.o\
	Molecule.o\
	IchorInterface.o\
	Integral_interfaceDEC.o\
	Integral_interface.o\
	Integral_interfaceDF.o\
	Integral_interface_test.o\
        II_absval_int.o \
	II_dft_dftd.o \
	II_dft_int.o \
	II_dft_ksm.o \
	II_dft_ksm_worker.o \
	II_XC_interface.o\
        ThermiteDistribute.o \
        ThermiteDistribute2.o \
        ThermiteDistributeGen.o \
        ThermiteDistributeK.o \
        ThermiteDistributeK2.o \
        ThermiteDistributeDEC.o \
        ThermiteDriver.o \
        ThermiteOverlapDistribution.o \
        ThermiteIntegrals.o \
        ThermiteIntTransform.o \
        ThermiteProp.o \
        ThermiteMem.o \
        THC_grid.o\
        THC_util.o\
        ODbatches.o     \
        pari.o     \
        lsmpi.o     \
        ReadMolFile.o  \
	ls_IntegralInterface.o\
	linsolv_df.o 
F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

# we should really have a separate target for bsm to properly
# track dependencies but we would not like it to be too troublesome.
libLSint.a: $(F_OBJS) $(F90_OBJS) 
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD


Makefile.depend depend: $(F_OBJS:.o=.F) $(F90_OBJS:.o=.F90)
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -M../lsutil -M../xcfun_host -M../dft -M../IchorIntegrals -- $(F_OBJS:.o=.F) $(F90_OBJS:.o=.F90) > Makefile.depend

include Makefile.depend
######################################################################
# read dependencies using sinclude (silent include)
# which continues if Makefile.depend is not readable (incl. doesn't exist!).
# Makefile.depend can be created with 'make depend'.
#-include Makefile.depend
###################################################################### 
