all: libgp.a

include ../Makefile.config

F_OBJS = gphjj.o    gpjeppe.o  gptrygve.o     mempkg.o \
	qpack.o    graphic.o  vcraypack.o    gpsaue.o \
	mpimacro.o pvmmacro.o gpdahle.o      intf7790.o \
	gpsimen.o $(GP_EXTRAS)

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > gp.ftnchek

all: libgp.a mpi_dummy.o mpi_dummyc.o

libgp.a: $(F_OBJS) $(C_OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

Makefile.depend depend: ../include/*.h
	../tools/check_depend $(INCLUDES) > Makefile.depend

include Makefile.depend
