all : libdft.a

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
        dft_grad.F   \
        dft_hes.F    \
        dft_aux.F    \
        dft_int.F    \
        dft_inp.F    

#        dft_ac.F     dft_old.F

# we define functional files separately since they are used also 
# for the fun-tester.
FUNFILES = \
	fun-becke.c	\
	fun-cam-b3lyp.c \
	fun-example.c	\
	fun-example2.c	\
	fun-example3.c	\
	fun-example4.c	\
	fun-example5.c	\
	fun-example6.c	\
	fun-example7.c	\
	fun-example8.c	\
	fun-example9.c	\
	fun-gga.c	\
	fun-kt.c	\
	fun-lb94.c	\
	fun-lyp.c	\
	fun-optx.c	\
	fun-p86c.c	\
	fun-pbe.c	\
	fun-pw86x.c	\
	fun-pw91cl.c	\
	fun-pw92.c	\
	fun-pw92pe.c	\
	fun-pwggaIIc.c	\
	fun-pwggaIIc2.c	\
	fun-pwggaIIx.c	\
	fun-pz81.c	\
	fun-slater.c	\
	fun-vwn.c

CFILES =\
	$(FUNFILES)	\
	functionals.c	\
	general.c	\
	grid-gen.c	\
	grid-gen2.c	\
	basisinfo.c	\
	gen-so-int.c	\
	integrator.c	\
	prop-eval.c	\
	quad-fast.c	\
	quad-strict.c

EXTRA_FILES = fun-tester.c

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > dft.ftnchek

OBJS = $(FFILES:.F=.o) $(CFILES:.c=.o)

clean :
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

depend:
	../tools/check_depend $(INCLUDES) > Makefile.depend

libdft.a: $(OBJS) Makefile
	-rm $@
	$(AR) $(ARFLAGS) $@ $(OBJS)


fun-tester: fun-tester.c functionals.c $(FUNFILES:.c=.o)
	$(CC) -g -O -I../include -DNO_BACKWARD_COMP $^ -o $@  -lm

fun4-tester: fun4-tester.c functionals.c $(FUNFILES)
	$(CC) -g -O -I../include -DNO_BACKWARD_COMP $^ -o $@  -lm


include Makefile.depend
