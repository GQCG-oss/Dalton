include ../Makefile.config

all : liberi.a

liberi.a : eri2cft.o        eri2her.o     eri2sym.o \
eri2aob.o     eri2dis.o     eri2odb.o     eri2vec.o \
eri2car1.o    eri2drv.o     eri2out.o     eri2pro.o \
eri2car2.o    eri2fck.o     eri2sph.o     erisph.o
	$(AR) $(ARFLAGS) $@ $?

slave : eri2par.o

node : eriprg.o

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a
