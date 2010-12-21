all: liberi.a

include ../Makefile.config


F_OBJS =  eri2cft.o        eri2her.o     eri2sym.o \
eri2aob.o     eri2dis.o     eri2odb.o     eri2vec.o \
eri2car1.o    eri2drv.o     eri2out.o     eri2pro.o \
eri2car2.o    eri2fck.o     eri2sph.o     eri2ave.o \
erisph.o      eri2r12.o     eri2csi.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > eri.ftnchek

liberi.a : $(F_OBJS)
	$(AR) $(ARFLAGS) $@ $?

slave : eri2par.o

node : eriprg.o

clean : 
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

Makefile.depend depend: ../include/*.h
	../tools/check_depend $(INCLUDES) > Makefile.depend

include Makefile.depend
eri2cft.o: eri2cft.F
	$(F90) $(INCLUDES) $(CPPFLAGS) $(SAFEFFLAGS) -c $*.F 
