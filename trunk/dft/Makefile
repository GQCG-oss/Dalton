all: libdft.a

include ../Makefile.config

FFILES = \
        dft_den.F    \
        dft_grid.F   \
        dft_ksm.F    \
        dft_rho.F    \
        dft_aos.F    \
        dft_exp.F    \
        dft_gridLL.F \
        dft_mag.F    \
	dft_mos.F    \
        dft_grad.F   \
        dft_hes.F    \
        dft_aux.F    \
        dft_int.F    \
        dft_inp.F    \
	dso_cb.F

#        dft_ac.F     dft_old.F

# we define functional files separately since they are used also 
# for the fun-tester.
FUNFILES = \
	fun-b86.c	\
	fun-b86mx.c	\
	fun-b97-1.c	\
	fun-b97-2.c	\
	fun-b97-3.c	\
	fun-b97-k.c	\
	fun-b97.c	\
	fun-becke.c	\
	fun-becke35.c	\
	fun-cam-b3lyp.c \
	fun-dk87x.c	\
	fun-example.c	\
	fun-example2.c	\
	fun-example3.c	\
	fun-example4.c	\
	fun-example5.c	\
	fun-example6.c	\
	fun-example7.c	\
	fun-example8.c	\
	fun-example9.c	\
	fun-g96.c	\
	fun-gga.c	\
	fun-hcth93.c	\
	fun-hcth93m.c	\
	fun-hcth120.c	\
	fun-hcth147.c	\
	fun-hcth407.c	\
	fun-hcth407p.c	\
	fun-kt.c	\
	fun-lb94.c	\
	fun-lg93.c	\
	fun-lrc95x.c	\
	fun-lyp.c	\
	fun-lypr.c	\
	fun-mpbex.c	\
	fun-mpwx.c	\
	fun-optx.c	\
	fun-p86c.c	\
	fun-pbec.c	\
	fun-pbex.c	\
	fun-pw86x.c	\
	fun-pw91c.c	\
	fun-pw91nc.c	\
	fun-pw91x.c	\
	fun-pw91x2.c	\
	fun-pw92ac.c	\
	fun-pw92c.c	\
	fun-pz81.c	\
	fun-revpbex.c	\
	fun-rpbex.c	\
	fun-slater.c	\
	fun-vwn.c       \
        fun-wignerc.c   \
        fun-wl90c.c

CFILES =\
	$(FUNFILES)	\
	basisinfo.c	\
	cube-fast.c	\
	eval-dso.c	\
	functionals.c	\
	gen-so-int.c	\
	general.c	\
	grid-gen.c	\
	grid-gen2.c	\
	integrator.c	\
	prop-eval.c	\
	quad-fast.c	\
	quad-faster.c	\
	quad-strict.c

EXTRA_FILES = fun-tester.c

F_CPPS = $(FFILES:.F=.i)

pre: $(F_CPPS)

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > dft.ftnchek

OBJS = $(FFILES:.F=.o) $(CFILES:.c=.o)

clean :
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

Makefile.depend depend: *.h ../include/*.h
	../tools/check_depend $(INCLUDES) > Makefile.depend

libdft.a: $(OBJS) Makefile
	-rm -f $@
	$(AR) $(ARFLAGS) $@ $(OBJS)


fun-tester: fun-tester.c functionals.c $(FUNFILES:.c=.o)
	$(CC) -g -O -I../include -DNO_BACKWARD_COMP $^ -o $@  -lm $(LDLIBS)

fun4-tester: fun4-tester.c functionals.c $(FUNFILES)
	$(CC) -g -O -I../include -DNO_BACKWARD_COMP $^ -o $@  -lm


include Makefile.depend
