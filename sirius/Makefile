include ../Makefile.config

all : libsirius.a

OBJS = ciauth.o          ciham.o           cioff.o\
cisig.o           csfstuff.o        dens.o\
detfi.o           detpr.o           dets0.o\
dets1.o           drcctl.o          getcix.o\
getint.o          siraba.o          sirave.o\
sirbkd.o          sircan.o          sirci.o\
sircmo.o          sirctl.o\
sirdiis.o         sirfck.o          sirgp.o\
sirgrad.o         sirinp.o          sirlintrn.o\
sirmp2.o          sirneo.o          sirnexkap.o\
sirnr.o           siropt.o          sirorb.o\
sirout.o          sirpop.o          sirsav.o\
sirset.o          sirtest.o         sirtr1.o\
sirtra.o          sirupd.o\
sirsol.o          symchk.o\
newtrahso.o       sirntra.o         intsort.o

libsirius.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $?

# sirtra.o has to be compiled without less optimization under HP-UX

sirtra.o: sirtra.F
	$(F77) $(INCLUDE_DIR) $(CPP_FLAGS) $(SAFEOPTFLAGS) -c $*.F

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a

depend:; ../tools/check_depend $(INCLUDE_DIR) > Makefile.depend

######################################################################
# read dependencies:
include Makefile.depend
######################################################################
