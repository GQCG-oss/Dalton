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
        dft_aux.F    \
        dft_grad.F   \
        dft_int.F    \
        dft_inp.F

#        dft_ac.F     dft_old.F

# we define functional files separately since they are used also 
# for the fun-tester.
FUNFILES = \
	fun-becke.c	\
	fun-dirac.c	\
	fun-gga.c	\
	fun-lb94.c	\
	fun-lyp.c	\
	fun-example.c	\
	fun-example2.c	\
	fun-example3.c	\
	fun-example4.c	\
	fun-example5.c	\
	fun-example6.c	\
	fun-example7.c	\
	fun-example8.c	\
	fun-example9.c	\
	fun-p86c.c	\
	fun-pbe.c	\
	fun-pw86x.c	\
	fun-pw92.c	\
	fun-pwggaIIx.c	\
	fun-pwggaIIc.c	\
	fun-pwggaIIc2.c	\
	fun-pw91cl.c	\
	fun-pz81.c	\
	fun-vwn.c

CFILES =\
	$(FUNFILES)	\
	functionals.c	\
	general.c	\
	grid-gen.c	\
	gen-so-int.c	\
	integrator.c	\
	prop-eval.c	\
	quad-fast.c	\
	quad-strict.c

EXTRA_FILES = fun-tester.c

# this ifeq section below is GNU specific thing. Feel free to remove it.
ifeq ($(CC),gcc)
FFLAGS += -Wall -Wunused -Wuninitialized 
CFLAGS += -Wall -Wunused -Wuninitialized
endif


CCOMPILE = $(CC)  $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)
FCOMPILE = $(F77) $(DEFS) $(INCLUDES) $(CPPFLAGS) $(FFLAGS)

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


fun-tester: fun-tester.c functionals.c $(FUNFILES)
	$(CC) -g -O -I../include -DNO_BACKWARD_COMP $^ -o $@  -lm

ifeq ($(CC),gcc)
# dependency part, should go to a separate, dynamically-created dependency 
# file.
DEP_FILES := $(addprefix .deps/, $(OBJS:.o=.P))
DEPS_MAGIC := $(shell mkdir .deps > /dev/null 2>&1 || :)

sinclude $(DEP_FILES)

# This is  adapted automake output.
# This is GNU specific thing. Feel free to remove it.
%.o: %.c
	@echo '$(CCOMPILE) -c $<'; \
	$(CCOMPILE) -Wp,-MD,.deps/$(*F).pp -c $<
	@-cp .deps/$(*F).pp .deps/$(*F).P; \
	tr ' ' '\012' < .deps/$(*F).pp \
	  | sed -e 's/^\\$$//' -e '/^$$/ d' -e '/:$$/ d' -e 's/$$/ :/' \
	    >> .deps/$(*F).P; \
	rm .deps/$(*F).pp
endif
ifeq ($(F77),g77)
%.o: %.F
	@echo '$(FCOMPILE) -c $<'; \
	$(FCOMPILE) -Wp,-MD,.deps/$(*F).pp -c $<
	@-cp .deps/$(*F).pp .deps/$(*F).P; \
	tr ' ' '\012' < .deps/$(*F).pp \
	  | sed -e 's/^\\$$//' -e '/^$$/ d' -e '/:$$/ d' -e 's/$$/ :/' \
	    >> .deps/$(*F).P; \
	rm .deps/$(*F).pp
endif
