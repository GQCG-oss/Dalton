all: libdft.a

include ../Makefile.config

INCLUDES := $(INCLUDES) -I../LSint -I../linears -I../lsutil

F77FILES = \
        dft_gridLL.F

F90FILES = \
	II_util.f90

# we define functional files separately since they are used also 
# for the fun-tester.
FUNFILES = \
	fun-becke.c	\
	fun-cam-b3lyp.c \
	fun-example.c	\
	fun-b97-1.c	\
	fun-gga.c	\
	fun-hcth120.c       \
	fun-hcth147.c       \
	fun-hcth407.c       \
	fun-hcth93.c        \
	fun-kt.c	\
	fun-lb94.c	\
	fun-lyp.c	\
	fun-optx.c	\
	fun-p86c.c	\
	fun-pbec.c	\
	fun-pbex.c	\
	fun-pw86x.c	\
	fun-pw91c.c	\
	fun-pz81.c	\
	fun-slater.c	\
	fun-vwn.c

CFILES =\
	$(FUNFILES)	\
	functionals.c	\
	general.c	\
	grid-gen.c	

EXTRA_FILES = fun-tester.c

F_CPPS = $(F77FILES:.F=.i) $(F90FILES:.f90=.i)

pre: $(F_CPPS)

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > dft.ftnchek

OBJS = $(F77FILES:.F=.o) $(F90FILES:.f90=.o) $(CFILES:.c=.o)

clean :
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

Makefile.depend depend:  $(F77FILES) $(F90FILES) $(CFILES)
	../tools/check_depend $(INCLUDES) > Makefile.depend
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -M../LSint -M../linears -M../lsutil -- $(F90FILES) >> Makefile.depend

libdft.a: $(OBJS) Makefile
	-rm $@
	$(AR) $(ARFLAGS) $@ $(OBJS)


fun-tester: fun-tester.c functionals.c $(FUNFILES:.c=.o)
	$(CC) -g -O -I../include -DNO_BACKWARD_COMP $^ -o $@  -lm $(LDLIBS)

fun4-tester: fun4-tester.c functionals.c $(FUNFILES)
	$(CC) -g -O -I../include -DNO_BACKWARD_COMP $^ -o $@  -lm


include Makefile.depend
