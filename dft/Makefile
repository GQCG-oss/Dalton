include ../Makefile.config

all : libdft.a

OBJS =  dftdrv.o dftgrd.o  dftden.o dftaos.o dftaux.o \
        dftfun.o  dftest.o dftac.o dftpar.o

libdft.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean :
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a

depend:; ../tools/check_depend $(INCLUDE_DIR) > Makefile.depend

######################################################################
# read dependencies:
include Makefile.depend
######################################################################
