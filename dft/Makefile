all : libdft.a

include ../Makefile.config

OBJS =  dftdrv.o dftgrd.o  dftden.o dftaos.o dftaux.o \
        dftfun.o  dftest.o dftac.o dftpar.o Lebedev-Laikov.o

libdft.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean :
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a

depend:; ../tools/check_depend $(INCLUDE_DIR) > Makefile.depend

######################################################################
# read dependencies:
include Makefile.depend
######################################################################
