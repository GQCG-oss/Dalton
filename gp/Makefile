include ../Makefile.config

F_OBJS =	gphjj.o    gpjeppe.o  gptrygve.o     mempkg.o \
	qpack.o    graphic.o  vcraypack.o    gpsaue.o \
	mpimacro.o pvmmacro.o gpdahle.o      intf7790.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

all: libgp.a

libgp.a: $(F_OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a

depend:; ../tools/check_depend $(INCLUDE_DIR) > Makefile.depend

######################################################################
# read dependencies using sinclude (silent include)
# which continues if Makefile.depend is not readable (incl. doesn't exist!).
# Makefile.depend can be created with 'make depend'.
sinclude Makefile.depend
######################################################################

