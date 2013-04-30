all: libxcfun_host.a

include ../Makefile.config

INCLUDES := $(INCLUDES) -I../lsutil 

# Keep the file lists sorted for easier modification.

F90_OBJS = \
	xcfun_host.o
F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

# we should really have a separate target for bsm to properly
# track dependencies but we would not like it to be too troublesome.
libxcfun_host.a: $(F_OBJS) $(F90_OBJS) 
#	if grep -q HAVE_BSM ../Makefile.config; then cd bsm && $(MAKE) install prefix=../.. AR="$(AR)"; fi
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD


Makefile.depend depend: $(F90_OBJS:.o=.F90)
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -M../lsutil -- $(F90_OBJS:.o=.F90) > Makefile.depend

include Makefile.depend
######################################################################
# read dependencies using sinclude (silent include)
# which continues if Makefile.depend is not readable (incl. doesn't exist!).
# Makefile.depend can be created with 'make depend'.
#-include Makefile.depend
###################################################################### 
