all: libopenrsp.a

include ../../../LSDALTON/Makefile.config

INCLUDES := $(INCLUDES) -I../../../LSDALTON/lsutil -I../../matrix-defop/src -I../../../LSDALTON/LSint -I../../../LSDALTON/SolverUtilities -I../../../LSDALTON/responsesolver

F90_OBJS = \
	rsp_contribs.o \
	rsp_equations.o

libopenrsp.a: $(F90_OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.a *.mod *.MOD

Makefile.depend depend: $(F90_OBJS:.o=.F90)
	../../../LSDALTON/tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -- $(F90_OBJS:.o=.F90) > Makefile.depend

%.o: %.mod
