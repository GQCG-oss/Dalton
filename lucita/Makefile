include ../Makefile.config

all : liblucita.a

F_OBJS = 	par_io.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > lucita.ftnchek

liblucita.a: $(F_OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean :	
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

Makefile.depend depend: ../include/*.h
	../tools/check_depend $(INCLUDES) > Makefile.depend

include Makefile.depend
par_io.o: par_io.f90
	$(F77) $(INCLUDES) $(CPPFLAGS) $(SAFE2_FFLAGS) -c $*.f90
