include ./Makefile.config
# DALTON_LIBS, MODULES and SUBDIRS must be modified at the same time!

# VAR_RSP
#-L../external/openrsp/src -lopenrsp 
#-L../external/matrix-defop/src -lmatrix-defop 

DALTON_LIBS =             \
-Llsdaltonsrc -llsdaltonsrc\
-Lpcm -llspcm             \
-Llinears -llinears       \
-Lddynam -lddynam         \
-Lgeomopt -lgeomopt       \
-Ldeccc -ldeccc           \
-Lrsp_properties -lrsp_properties  \
-Lresponsesolver -lresponsesolver \
-LSolverUtilities -lsolverutilities \
-Ldft -ldft		  \
-Lpbc2 -lpbc2            \
-LLSint -lLSint           \
-L$(FMMDIR) -lmm          \
-Linterest/lib -linterest \
-Lpdpack -lpdpack         \
-Lxcfun_host -lxcfun_host  \
-Lcuda -lcuda\
-Llsutil -llsutil

# VAR_RSP needs MDO_OBJ ORS_OBJ 
MODULES = LSU_OBJ CUD_OBJ XFH_OBJ PD_OBJ INT_OBJ FMM_OBJ LSI_OBJ PBC_OBJ2 \
	DFT_OBJ SOU_OBJ RSO_OBJ RSP_OBJ DEC_OBJ GEO_OBJ DYN_OBJ LIN_OBJ LSD_OBJ \
	PCM_OBJ

# VAR_RSP needs ../external/matrix-defop/src ../external/openrsp/src
SUBDIRS = lsutil cuda xcfun_host pdpack interest/lib $(FMMDIR) LSint pbc2 dft SolverUtilities \
	responsesolver rsp_properties deccc geomopt ddynam linears pcm lsdaltonsrc
PDLIB = pdpack/pdlib.a

OBJPRG = lsdaltonsrc/lsdalton_wrapper.f90

#
#     Most common build of LSDalton
#
lsdalton.x: $(MODULES)
	@echo "---------------> Linking sequential lsdalton.x ..."
	$(F90) $(FFLAGS) \
	-o $(INSTALLDIR)/lsdalton.x $(OBJPRG) $(DALTON_LIBS) $(LIBS)

#
#       Linux version of the program
#
linux.x: $(MODULES)
	@echo "---------------> Linking sequential lsdalton.x ..."
	$(F90) $(FFLAGS) \
	-o $(INSTALLDIR)/lsdalton.x $(OBJPRG) $(DALTON_LIBS) $(LIBS)
#
#       Build tools: make tool_name.x
#
%.x: $(MODULES)
	@echo "---------------> Linking tool" $*.x
	cd tools && $(MAKE) $*.o
	$(F90) $(FFLAGS) \
	-o $(INSTALLDIR)/$*.x tools/$*.o $(DALTON_LIBS) $(LIBS)
#
#     Linux MPI parallel build (first create sequential build lsdalton.x, then dalpar.x)
#
#linuxparallel.x: $(MODULES)
#	@echo "---------------> Linking parallel dalpar.x ..."
#	$(F90) $(FFLAGS) \
#	-o $(INSTALLDIR)/lsdalpar.x $(OBJPRG) $(OBJSMEM) \
#	$(OBJSLAVE) $(OBJSAMFI) $(DALTON_LIBS) $(LIBS) \
#	$(MPI_LIB_PATH) $(MPI_LIB) 
##
##	MPI parallel build (first create sequential build lsdalton.x, then dalpar.x)
##
#parallel.x: lsdalton.x
#	@echo "---------------> Linking parallel dalpar.x ..."
#	$(F90) $(FFLAGS) \
#	-o $(INSTALLDIR)/lsdalpar.x $(OBJPRG) $(OBJSLAVE) $(OBJSAMFI) \
#	$(DALTON_LIBS) $(LIBS) $(MPI_LIB_PATH) $(MPI_LIB)
##
##	Build with PVMe
##	We will never need MPILIB when using PVM
##       First build master, then slave
##
#dalpvm.x: $(MODULES)
#	$(F90) $(FFLAGS) \
#	-o $(INSTALLDIR)/lsdalpvm.x $(OBJPRG) $(OBJSMEM) \
#	$(OBJSLAVE) $(OBJSAMFI) $(LIBS) \
#	$(DALTON_LIBS) $(PVM_LIB_PATH) \
#	$(PVM_LIBS) $(PVM_INC_PATH)
#	$(F90) $(FFLAGS) -o $(INSTALLDIR)/lsdal_node.x $(OBJNODE) $(OBJSLAVE) \
#	$(OBJSAMFI) $(LIBS) $(DALTON_LIBS) $(PVM_LIB_PATH) \
#	$(PVM_LIBS)  $(PVM_INC_PATH)
##
##	Not tested MPI on Cray yet. Thus no MPILIB, and no OBJSMXM nor
##	OBJSEIS
##
#cray.x: $(MODULES)
#	$(F90) $(FFLAGS) -o $(INSTALLDIR)/lsdalton.x $(OBJPRG) \
#	$(OBJSLAVE) $(OBJSAMFI) $(LIBS) $(DALTON_LIBS)
##
##	This is a proper build for the Cray-T3D
##
#t3d.x: $(MODULES)
#	$(F90) $(FFLAGS) $(MPI_LIB_PATH) $(MPI_LIB) \
#	-o $(INSTALLDIR)/lsdalpar.x $(OBJPRG) $(OBJSLAVE) \
#	$(LIBS) $(OBJSAMFI) $(DALTON_LIBS)
#
#t90.x: $(MODULES)
#	$(F90) $(FFLAGS) $(MPI_LIB_PATH) $(MPI_LIB) \
#	-o $(INSTALLDIR)/lsdalpar.x $(OBJPRG) $(OBJSLAVE) \
#	$(LIBS) $(OBJSAMFI) $(DALTON_LIBS)
#
#int_tester: dft/int_tester.f90 $(MODULES)
#	$(F90) $(F90FLAGS)  -o $@ dft/int_tester.f90 \
#	$(OBJSAMFI) $(OBJSLAVE) $(OBJSMEM) $(DALTON_LIBS) $(LIBS)
#
#
#	Update all dependencies
#
depend  :
	for i in $(SUBDIRS); do ( cd $$i  && touch Makefile.depend && $(MAKE) $@ ); done

#
#	Make it a bit cleaner, remover all .o/.lst/*.f -files
#
clean :
	for i in $(SUBDIRS); do ( cd $$i  && $(MAKE) $@ ); done
	$(RM) -f *~
#
#	We remove the entire source code as well if we do not plan to debug
#
veryclean : 
	$(RM) -rf $(SUBDIRS)

LSD_OBJ	: LSU_OBJ CUD_OBJ LSI_OBJ DFT_OBJ GEO_OBJ DYN_OBJ DEC_OBJ INT_OBJ LIN_OBJ PBC_OBJ2 RSP_OBJ RSO_OBJ SOU_OBJ
	cd lsdaltonsrc && $(MAKE) all

PD_OBJ	:
	cd pdpack && $(MAKE) all

LSU_OBJ :
	cd lsutil && $(MAKE) all

CUD_OBJ : LSU_OBJ
	cd cuda && $(MAKE) all

XFH_OBJ : LSU_OBJ
	cd xcfun_host && $(MAKE) all

DFT_OBJ	: LSU_OBJ
	cd dft && $(MAKE) all

GEO_OBJ : LSU_OBJ LSI_OBJ
	cd geomopt && $(MAKE) all

DYN_OBJ : LSU_OBJ LSI_OBJ GEO_OBJ
	cd ddynam && $(MAKE) all

LSI_OBJ : LSU_OBJ INT_OBJ XFH_OBJ DFT_OBJ 
	cd LSint && $(MAKE) all

DEC_OBJ : LSU_OBJ LSI_OBJ LIN_OBJ CUD_OBJ
	cd deccc && $(MAKE) all

RSP_OBJ : LSU_OBJ LSI_OBJ LIN_OBJ
	cd rsp_properties && $(MAKE) all

RSO_OBJ : LSU_OBJ LSI_OBJ SOU_OBJ
	cd responsesolver && $(MAKE) all

INT_OBJ : 
	cd interest/lib && $(MAKE) all

SOU_OBJ : LSU_OBJ LSI_OBJ
	cd SolverUtilities && $(MAKE) all

PCM_OBJ : LSU_OBJ LSI_OBJ LIN_OBJ
	cd pcm && $(MAKE) all

# VAR_RSP needs 
#MDO_OBJ : LSU_OBJ
#	cd ../external/matrix-defop/src && $(MAKE) all
#
#ORS_OBJ : LSU_OBJ MDO_OBJ LSI_OBJ SOU_OBJ RSO_OBJ
#	cd ../external/openrsp/src && $(MAKE) all
#
# FOR VAR_RSP LIN_OBJ depends on ORS_OBJ MDO_OBJ
LIN_OBJ : LSU_OBJ LSI_OBJ SOU_OBJ RSO_OBJ 
	cd linears && $(MAKE) all

FMM_OBJ : LSU_OBJ
	cd $(FMMDIR) && $(MAKE) all

PBC_OBJ2: LSU_OBJ LSI_OBJ
	cd pbc2 && $(MAKE) all

ftnchek: pre
	ftnchek $(CHEKFLAGS) */*.i > dalton.ftnchek

pre :
	for i in $(SUBDIRS); do ( cd $$i  && $(MAKE) $@ ); done

check : default
	cd test && ./TEST -y linsca
