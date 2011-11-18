include ../Makefile.config

all : libamfi.a

F_OBJS = 	amfi.o symtra.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > amfi.ftnchek

libamfi.a: $(F_OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean :	
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

# ifort has (had?) problems compiling amfi.F
amfi.o: amfi.F
	$(F90) $(INCLUDES) $(CPPFLAGS) $(SAFE_FFLAGS_for_ifort) -c $*.F 

Makefile.depend depend: ../include/*.h
	../tools/check_depend $(INCLUDES) > Makefile.depend

include Makefile.depend
# -- end of Makefile for amfi/
