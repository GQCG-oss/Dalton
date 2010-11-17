#
# FILE    : Makefile for dalton version of lucita
#           hjaaj May 2008
#

all : liblucita.a

include ../Makefile.config

F77_OBJS = ciinfo.o  density.o  diagonal.o  diagonalize.o  \
           fock.o  dalton_interface.o  io.o  \
           memory.o program.o \
           property.o  sigma.o  strings.o  transform.o  utils.o  \
           vecmat.o  par_lucia.o  par_solver.o \
           gplucia.o gplupar.o gplucia_io.o
# files gp/gplucia.F gp/gplupar.F include/parluci.h include/luctinp.h
# and include/mxmcgas.h in DIRAC are all moved to the lucita/
# in dalton (for the time being). /hjaaj May 2008

F90_OBJS = par_io.o

F90_MODS =

F77_CPPS = $(F77_OBJS:.o=.i)

pre:	$(F77_CPPS)

ftnchek:pre
	ftnchek $(CHEKFLAGS) $(F77_CPPS) > ftnchek.out

liblucita.a : $(F77_OBJS)
	$(AR) $(ARFLAGS) $@ $?


clean : 
	$(RM) -f *.o *~ *.mod *.lst *.i

mod:
	$(MAKE) f90dep  # create list of which file to compile to get the .mod files
	$(MAKE) module_files

module_files: $(F90_MODS)

depend:
	makedepend -- $(DEPEND_FLAGS) -- $(F77_OBJS:.o=.F) $(F90_OBJS:.o=.F90)

f90dep: ;\
       ../utils/sfmakedepend.pl -I $(SOURCEDIR)/include -M $(SOURCEDIR)/modules $(F90_OBJS:.o=.F90) $(F90_OBJS:.o=.F) $(F77_OBJS:.o=.F)
 

# DO NOT DELETE THIS LINE -- make depend depends on it.

ciinfo.o: mxpdim.inc cgas.inc orbinp.inc wrkspc.inc lucinp.inc glbbas.inc
ciinfo.o: crun.inc strinp.inc cicisp.inc clunit.inc csm.inc cstate.inc
ciinfo.o: stinf.inc strbas.inc cprnt.inc oper.inc gasstr.inc parluci.h
ciinfo.o: spinfo.inc multd2h.inc machine.inc SysDef.inc irat.inc
density.o: mxpdim.inc orbinp.inc cicisp.inc strbas.inc cstate.inc strinp.inc
density.o: stinf.inc csm.inc wrkspc.inc crun.inc cgas.inc gasstr.inc
density.o: cprnt.inc spinfo.inc glbbas.inc lucinp.inc clunit.inc oper.inc
density.o: multd2h.inc
diagonal.o: mxpdim.inc orbinp.inc cicisp.inc strbas.inc cstate.inc strinp.inc
diagonal.o: stinf.inc csm.inc wrkspc.inc cprnt.inc cgas.inc gasstr.inc
diagonal.o: oper.inc crun.inc glbbas.inc cecore.inc parluci.h lucinp.inc
diagonalize.o: parluci.h files.inc mxpdim.inc orbinp.inc cicisp.inc
diagonalize.o: strbas.inc cstate.inc strinp.inc stinf.inc csm.inc wrkspc.inc
diagonalize.o: gasstr.inc cgas.inc lucinp.inc cprnt.inc glbbas.inc oper.inc
diagonalize.o: spinfo.inc csfbas.inc clunit.inc crun.inc
fock.o: mxpdim.inc wrkspc.inc orbinp.inc glbbas.inc lucinp.inc cgas.inc
fock.o: cicisp.inc clunit.inc csm.inc cstate.inc crun.inc strinp.inc
fock.o: stinf.inc strbas.inc cprnt.inc oper.inc
dalton_interface.o: luctinp.h clunit.inc mxpdim.inc luci.inc units.inc
dalton_interface.o: parluci.h glbbas.inc wrkspc.inc crun.inc cgas.inc
dalton_interface.o: lucinp.inc orbinp.inc multd2h.inc SysDef.inc cc_exc.inc
dalton_interface.o: intform.inc csm.inc orbinf_lucita.inc irat.inc
dalton_interface.o: machine.inc files.inc
io.o: parluci.h files.inc clunit.inc units.inc
memory.o: mxpdim.inc wrkspc.inc lucinp.inc orbinp.inc cstate.inc csm.inc
memory.o: crun.inc cprnt.inc glbbas.inc irat.inc
program.o: mxpdim.inc cicisp.inc wrkspc.inc orbinp.inc clunit.inc csm.inc
program.o: cstate.inc crun.inc strinp.inc stinf.inc strbas.inc glbbas.inc
program.o: cprnt.inc oper.inc gasstr.inc cgas.inc lucinp.inc intform.inc
program.o: comjep.inc cc_exc.inc parluci.h files.inc units.inc luctinp.h
program.o: machine.inc
property.o: mxpdim.inc wrkspc.inc glbbas.inc lucinp.inc orbinp.inc oper.inc
property.o: crun.inc multd2h.inc cprnt.inc cicisp.inc clunit.inc csm.inc
property.o: cstate.inc strinp.inc stinf.inc strbas.inc cgas.inc
sigma.o: mxpdim.inc gasstr.inc orbinp.inc cgas.inc crun.inc parluci.h
sigma.o: cicisp.inc strbas.inc cstate.inc strinp.inc stinf.inc csm.inc
sigma.o: wrkspc.inc lucinp.inc cprnt.inc glbbas.inc oper.inc comjep.inc
strings.o: mxpdim.inc orbinp.inc wrkspc.inc strinp.inc strbas.inc cgas.inc
strings.o: gasstr.inc stinf.inc csm.inc lucinp.inc comjep.inc multd2h.inc
strings.o: cprnt.inc cicisp.inc cstate.inc irat.inc crun.inc
transform.o: mxpdim.inc glbbas.inc wrkspc.inc orbinp.inc lucinp.inc
transform.o: strinp.inc cicisp.inc strbas.inc gasstr.inc crun.inc csm.inc
transform.o: multd2h.inc clunit.inc cstate.inc oper.inc intform.inc
utils.o: units.inc irat.inc
vecmat.o: mxpdim.inc wrkspc.inc glbbas.inc multd2h.inc oper.inc orbinp.inc
vecmat.o: cicisp.inc crun.inc strbas.inc stinf.inc csm.inc cgas.inc
vecmat.o: gasstr.inc cstate.inc strinp.inc cprnt.inc clunit.inc lucinp.inc
vecmat.o: parluci.h
gplucia.o: parluci.h
