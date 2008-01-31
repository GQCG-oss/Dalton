include ./Makefile.config

# DALTON_LIBS, MODULES and SUBDIRS must be modified at the same time!
DALTON_LIBS = \
-Labacus -labacus	\
-Lrsp -lrsp 		\
-Lsirius -lsirius	\
-labacus	\
-Leri -leri		\
-Ldensfit -ldensfit	\
-Lcc  -lcc		\
-Ldft -ldft		\
-Lgp -lgp		\
-Lpdpack -lpdpack

MODULES = MAIN_OBJ ABA_OBJ SIR_OBJ RSP_OBJ GP_OBJ SLAVE_OBJ ERI_OBJ \
	DFIT_OBJ PD_OBJ CC_OBJ DFT_OBJ AMFI_OBJ

SUBDIRS = abacus sirius rsp gp cc eri densfit pdpack dft amfi 

OBJSLAVE = abacus/herpar.o eri/eri2par.o

OBJNODE = eri/eriprg.o

OBJSAMFI = amfi/amfi.o amfi/symtra.o

OBJS_MPI_DUMMY = gp/mpi_dummy.o gp/mpi_dummyc.o
#
#     Most common build of Dalton
#
dalton.x: $(MODULES)
	@echo "---------------> Linking sequential dalton.x ..."
	$(F77) $(FFLAGS) \
	-o $(INSTALLDIR)/dalton.x abacus/dalton.o $(OBJSLAVE) \
	$(OBJSAMFI) $(OBJS_MPI_DUMMY) $(DALTON_LIBS) $(LIBS)
#
#       Linux version of the program
#
linux.x: $(MODULES)
	@echo "---------------> Linking sequential dalton.x ..."
	$(F77) $(FFLAGS) \
	-o $(INSTALLDIR)/dalton.x abacus/dalton.o \
	$(OBJSLAVE) $(OBJSAMFI) $(OBJS_MPI_DUMMY) $(DALTON_LIBS) $(LIBS)
#
#     Linux MPI parallel build (first create sequential build dalton.x, then dalpar.x)
#
linuxparallel.x: linux.x
	@echo "---------------> Linking parallel dalpar.x ..."
	$(F77) $(FFLAGS) \
	-o $(INSTALLDIR)/dalpar.x abacus/dalton.o \
	$(OBJSLAVE) $(OBJSAMFI) $(DALTON_LIBS) $(LIBS) \
	$(MPI_LIB_PATH) $(MPI_LIB) 
#
#	MPI parallel build (first create sequential build dalton.x, then dalpar.x)
#
parallel.x: dalton.x
	@echo "---------------> Linking parallel dalpar.x ..."
	$(F77) $(FFLAGS) \
	-o $(INSTALLDIR)/dalpar.x abacus/dalton.o $(OBJSLAVE) $(OBJSAMFI) \
	$(DALTON_LIBS) $(LIBS) $(MPI_LIB_PATH) $(MPI_LIB)
#
#	Build with PVMe
#	We will never need MPILIB when using PVM
#       First build master, then slave
#
dalpvm.x: $(MODULES)
	$(F77) $(FFLAGS) \
	-o $(INSTALLDIR)/dalpvm.x abacus/dalton.o \
	$(OBJSLAVE) $(OBJSAMFI) $(LIBS) \
	$(DALTON_LIBS) $(PVM_LIB_PATH) \
	$(PVM_LIBS) $(PVM_INC_PATH)
	$(F77) $(FFLAGS) -o $(INSTALLDIR)/node.x $(OBJNODE) $(OBJSLAVE) \
	$(OBJSAMFI) $(LIBS) $(DALTON_LIBS) $(PVM_LIB_PATH) \
	$(PVM_LIBS)  $(PVM_INC_PATH)
#
#	Not tested MPI on Cray yet. Thus no MPILIB, and no OBJSMXM nor
#	OBJSEIS
#
cray.x: $(MODULES)
	$(F77) $(FFLAGS) -o $(INSTALLDIR)/dalton.x \
	$(OBJSLAVE) $(OBJSAMFI) $(LIBS) $(DALTON_LIBS)
#
#	This is a proper build for the Cray-T3D
#
t3d.x: $(MODULES)
	$(F77) $(FFLAGS) $(MPI_LIB_PATH) $(MPI_LIB) \
	-o $(INSTALLDIR)/dalpar.x $(OBJSLAVE) \
	$(LIBS) $(OBJSAMFI) $(DALTON_LIBS)

t90.x: $(MODULES)
	$(F77) $(FFLAGS) $(MPI_LIB_PATH) $(MPI_LIB) \
	-o $(INSTALLDIR)/dalpar.x $(IO_OBJS) $(OBJSLAVE) \
	$(LIBS) $(OBJSAMFI) $(DALTON_LIBS)
#
#	Update all dependencies
#
depend :
	for i in $(SUBDIRS); do ( cd $$i  && $(MAKE) $@ ); done
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
	$(RM) -rf abacus sirius rsp gp include eri densfit pdpack cc dft

MAIN_OBJ :
	cd abacus && $(MAKE) main

ABA_OBJ :
	cd abacus && $(MAKE) all

SIR_OBJ :
	cd sirius && $(MAKE) all

RSP_OBJ :
	cd rsp && $(MAKE) all

GP_OBJ :
	cd gp && $(MAKE) all

SLAVE_OBJ :
	cd abacus && $(MAKE) slave
	cd eri && $(MAKE) slave

NODE_OBJ :
	cd eri && $(MAKE) node

ERI_OBJ	:
	cd eri && $(MAKE) all

DFIT_OBJ	:
	cd densfit && $(MAKE) all

CC_OBJ	:
	cd cc && $(MAKE) all

DFT_OBJ	:
	cd dft && $(MAKE) all

IO_OBJ	:
	cd cc && $(MAKE) io

PD_OBJ	:
	cd pdpack && $(MAKE) all

AMFI_OBJ :
	cd amfi && $(MAKE) all

ftnchek: pre
	ftnchek $(CHEKFLAGS) */*.i > dalton.ftnchek

pre :
	for i in $(SUBDIRS); do ( cd $$i  && $(MAKE) $@ ); done
check:
	cd test && ./TEST -y short
