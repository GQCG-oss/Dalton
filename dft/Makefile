all : libdft.a

include ../Makefile.config

F_OBJS =  dftdrv.o dftgrd.o  dftden.o dftaos.o dftaux.o \
        dftfun.o  dftest.o dftac.o dftpar.o Lebedev-Laikov.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

libdft.a: $(F_OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean :
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a

depend:
	../tools/check_depend $(INCLUDE_DIR) > Makefile.depend

include Makefile.depend
