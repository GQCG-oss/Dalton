all: libsrdft.a

include ../Makefile.config

FFILES = \
	dftaux_srdft.F \
	dftdrv_srdft.F \
	dftfunjt.F \
	dftfun_srdft.F \
	dftgrd_srdft.F \
	sirmp2_srdft.F \
	sir_srdft.F \
	srdft.F \
	srdftjt.F

# dftpar.F dft_ac.F


# we define functional files separately since they are used also 
# for the fun-tester.
#FUNFILES = \
#	fun-b86.c	\
#  .. cut ..
#       fun-wl90c.c

#CFILES =\
#	$(FUNFILES)	\
#	basisinfo.c	\
#	quad-strict.c

#EXTRA_FILES = fun-tester.c

F_CPPS = $(FFILES:.F=.i)

pre: $(F_CPPS)

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > dft.ftnchek

OBJS = $(FFILES:.F=.o) $(CFILES:.c=.o)

clean :
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

Makefile.depend depend: *.h ../include/*.h
	../tools/check_depend $(INCLUDES) > Makefile.depend

libsrdft.a: Makefile.depend $(OBJS) Makefile
	-rm -f $@
	$(AR) $(ARFLAGS) $@ $(OBJS)

include Makefile.depend
