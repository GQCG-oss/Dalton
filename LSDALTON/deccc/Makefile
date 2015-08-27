all: libdeccc.a

include ../Makefile.config

INCLUDES := $(INCLUDES) -I../LSint -I../lsutil -I../linears -I../cuda

F90_OBJS = dec_settings.o \
        cc_tools.o \
        cc_response_tools.o \
	decmpi.o \
	decmpiSlave.o \
	CABS.o \
	crop_tools.o \
        array4_memory.o \
        array3_memory.o \
	dec_utils.o \
	dec_tools.o \
	ccarray2_simple.o \
	ccarray4_simple.o \
	ccarray3_simple.o \
	fullmolecule.o \
	mp2.o \
	rimp2.o \
	ri_util.o \
	ccsdpt_driver.o \
	ccsdpt_tools.o \
	ccsdpt_kernels.o \
	ccsdpt_full.o \
	ccsdpt_dec.o \
	ccorbital.o \
	cc_integrals.o \
	ccsd.o \
	pno_ccsd.o \
	rpa.o \
	f12_integrals.o \
	f12_routines.o \
	rif12_integrals.o \
	dec_f12routines.o
	ls_thc_rimp2.o\
	wangy_playground.o \
	cc_driver.o \
	dec_atom.o \
	rimp2_gradient.o \
	mp2_gradient.o \
        ccsd_gradient.o \
	fragment_energy.o \
	dec_driver.o \
	dec_driver_slave.o \
	full_driver.o \
	full_mp2.o \
	full_rimp2.o \
	full_rimp2f12.o \
	full_driver_f12contractions.o\
	full_ls_thc_rimp2.o\
	snoop_main.o\
	snoop_tools.o\
	dec_main.o



# we should really have a separate target for bsm to properly
# track dependencies but we would not like it to be too troublesome.
libdeccc.a: $(F90_OBJS) $(C_OBJS) 
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD


Makefile.depend depend: $(F90_OBJS:.o=.F90) $(C_OBJS:.o=.c) 
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -M../LSint -M../lsutil -M../cuda -- $(F90_OBJS:.o=.F90) > Makefile.depend

include Makefile.depend
######################################################################
# read dependencies using sinclude (silent include)
# which continues if Makefile.depend is not readable (incl. doesn't exist!).
# Makefile.depend can be created with 'make depend'.
#-include Makefile.depend
###################################################################### 

%.o: %.mod
