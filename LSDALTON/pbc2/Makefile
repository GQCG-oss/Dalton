all:libpbc2.a

include ../Makefile.config

INCLUDES := $(INCLUDES) -I../lsutil -I../LSint

F90_OBJS = \
			  pbc-msc.o\
			  pbc-krsp-op.o\
			  pbc-matop.o\
			  pbc-data.o\
			  pbc-ffdata.o\
			  pbc-nfinteractions.o\
			  pbcmain.o\
			  pbc-ffinteractions.o\
			  pbc-scf.o\
			  pbc_compare.o\
			  pbc-timings.o
F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > pbc2.ftnchek 

libpbc2.a : $(F_OBJS) $(F90_OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean: 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod

Makefile.depend depend:  $(F_OBJS:.o=.F) $(F90_OBJS:.o=.F90)
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -M../lsutil -M../LSint -- $(F90_OBJS:.o=.F90) > Makefile.depend

include Makefile.depend
