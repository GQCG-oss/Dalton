include ../Makefile.config

all : libsirius.a

libsirius.a :	ciauth.o          ciham.o           cioff.o\
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
	$(AR) $(ARFLAGS) $@ $?

#sircino.o snglet.o

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a

