include ../Makefile.config

OBJS =	gphjj.o    gpjeppe.o  gptrygve.o     mempkg.o \
	qpack.o    graphic.o  vcraypack.o    gpsaue.o \
	mpimacro.o pvmmacro.o gpdahle.o

all: libgp.a

libgp.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a
