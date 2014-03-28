all: libdeccc.a

include ../Makefile.config

INCLUDES := $(INCLUDES) -I../LSint -I../lsutil

F90_OBJS = dec_settings.o \
        dec_workarounds.o \
	decmpi.o \
	decmpiSlave.o \
	CABS.o \
	crop_tools.o \
        array4_memory.o \
        array3_memory.o \
	dec_utils.o \
	ccarray2_simple.o \
	ccarray4_simple.o \
	ccarray3_simple.o \
	fullmolecule.o \
	mp2.o \
	ccsdpt.o \
	ccorbital.o \
	cc_integrals.o \
        cc_debug_routines.o\
	ccsd.o \
	rpa.o \
	f12_integrals.o \
	f12_routines.o \
	wangy_playground.o \
	cc_driver.o \
	dec_atom.o \
	mp2_gradient.o \
	fragment_energy.o \
	dec_driver.o \
	dec_driver_slave.o \
	full_driver.o \
	full_driver_f12contractions.o\
	dec_main.o
C_OBJS = crayio.o



# we should really have a separate target for bsm to properly
# track dependencies but we would not like it to be too troublesome.
libdeccc.a: $(F90_OBJS) $(C_OBJS) 
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD


Makefile.depend depend: $(F90_OBJS:.o=.F90) $(C_OBJS:.o=.c) 
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -M../LSint -M../lsutil -- $(F90_OBJS:.o=.F90) > Makefile.depend

include Makefile.depend
######################################################################
# read dependencies using sinclude (silent include)
# which continues if Makefile.depend is not readable (incl. doesn't exist!).
# Makefile.depend can be created with 'make depend'.
#-include Makefile.depend
###################################################################### 

%.o: %.mod
