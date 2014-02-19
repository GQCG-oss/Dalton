all: libdft.a

include ../Makefile.config

INCLUDES := $(INCLUDES) -I../LSint -I../linears -I../lsutil

F90FILES = \
	II_util.F90

# we define functional files separately since they are used also 
# for the fun-tester.
FUNFILES = \
	fun-becke.c	\
	fun-cam-b3lyp.c \
	fun-camx.c \
	fun-cam_compx.c \
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
	fun-revpbex.c	\
	fun-rpbex.c  	\
	fun-mpbex.c  	\
	fun-pw91x.c  	\
	fun-g96.c    	\
	fun-lg93.c   	\
	fun-vwn.c

CFILES =\
	$(FUNFILES)	\
	functionals.c	\
	general.c	

EXTRA_FILES = fun-tester.c

F_CPPS = $(F90FILES:.F90=.i)

pre: $(F_CPPS)

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > dft.ftnchek

OBJS = $(F90FILES:.F90=.o) $(CFILES:.c=.o)

clean :
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

Makefile.depend depend:  $(F90FILES) $(CFILES)
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
