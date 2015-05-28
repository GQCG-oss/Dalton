all: libddynam.a

include ../Makefile.config

INCLUDES:=$(INCLUDES) -I../lsutil -I../LSint -I/geomopt -I/deccc

F90_OBJS =  Modules.o\
            LSinput.o\
            dyn_utilities.o\
            Fock_mat_dyn.o\
            Temperature.o\
            TimeRev_prop.o

F_CPPS = $(F90_FILES:.o=.i) 

pre: $(F_CPPS)

libddynam.a: $(F90_OBJS) 
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod

Makefile.depend depend: $(F90_OBJS:.o=.F90)
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -M../lsutil -M../LSint -M../geomopt -M../deccc -- $(F_OBJS:.o=.F) $(F90_OBJS:.o=.F90) > Makefile.depend

######################################################################
# read dependencies:
include Makefile.depend
###################################################################### 

%.o: %.mod

