include ../Makefile.config

all: librsp.a

F_OBJS = \
	deq27.o		\
	expone.o	\
	fone.o		\
	hsoal2.o	\
	hsoprp.o	\
	hypinp.o	\
	lagrang.o	\
	qfock.o		\
	qrhso.o		\
	rspave.o	\
	rspc8.o		\
	rspc82.o	\
	rspcr1.o	\
	rspcr2.o	\
	rspcr3.o	\
	rspcr4.o	\
	rspcr5.o	\
	rspcr6.o	\
	rspcr7.o	\
	rspcr8.o	\
	rspcr9.o	\
	rspdm.o		\
	rspe2c.o	\
	rspfxd.o	\
	rspmai.o	\
	rspmc1.o	\
	rspmc2.o	\
	rspmc3.o	\
	rspoit.o	\
	rspoli.o	\
	rsporp.o	\
	rspprp.o	\
	rspqrx3.o	\
	rspsol.o	\
	rspsol1.o	\
	rspsol2.o	\
	rspsoppa.o	\
	rspsoppa.o	\
	rsptmo1.o	\
	rsptmo2.o	\
	rsptmo3.o	\
	rsptmo4.o	\
	rsptpa1.o	\
	rsptpa2.o	\
	rsptpa3.o	\
	rspvec.o	\
	absorp.o	\
	x2expl.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

librsp.a: $(F_OBJS)
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

