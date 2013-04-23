all: liblsutil.a

include ../Makefile.config

# Keep the file lists sorted for easier modification.

C_OBJS = matop_csr_aux.o \
	myPAPI_set_inherit.o
F90_OBJS = \
	dftType.o\
	gridgeneration_memory.o\
	dft_operations.o\
	dftMem.o\
	ddynType.o\
	BasisinfoType.o\
	Basisinfo_operations.o\
        common_utilities.o \
        file-operations.o \
        grid_utilities.o \
	lsutilities.o\
	lowdin.o\
	ls_precision.o\
	ls_IO.o\
	ls_IOType.o\
	ls_parameters.o\
	ls_screen.o\
	ls_math.o\
	ptr_assoc.o \
	OverlapDistributionType.o\
	ls_screenMPI.o\
	lstensor_operations.o\
	lstensorType.o\
	lstensor_Mem.o\
	LSoptType.o\
	par_mod.o\
	ProfileType.o\
	lsmpiType.o\
	lsmpi-operations.o\
	matop_scalapack.o\
	matop_dense.o\
	matop_csr.o\
	lsmatop_dense.o\
	matop_dense_unrest.o\
	mat-operations-aux.o\
	mat-operations-essential.o\
	matrix_defop.o\
	matrix_defop_backend.o\
        matrix_utilities.o \
	matrix_module.o\
	MoleculeType.o\
	Molecule_operations.o\
        Integral_input.o\
        IntegralOutputType.o\
        IntegralOutput_operations.o\
	IntegralType.o\
	AOType.o\
	ODType.o\
	Build_AOBATCH.o\
	AO_operations.o\
	OD_operations.o\
	LSmatrixType.o\
        memory.o \
	Time.o\
	TestMPI.o\
	ks-settings.o\
	dec_typedef.o\
	SphCartMatrices.o\
	TYPE-DEF.o\
	TYPE-OP.o\
	pbc_lattice_type.o\
	pbc_lattice_vectors.o\
	AtomicSparse.o\
	f12.o\
	mat3dim.o\
	rsp-typedef.o\
	response_prop_type.o\
        fundamental.o\
	papi.o

#       AoBatchType.o\

# we should really have a separate target for bsm to properly
# track dependencies but we would not like it to be too troublesome.
liblsutil.a: $(F90_OBJS) $(C_OBJS)
	if grep -q HAVE_BSM ../Makefile.config; then cd bsm && $(MAKE) install prefix=../.. AR="$(AR)"; fi
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD


Makefile.depend depend: $(F90_OBJS:.o=.f90) $(C_OBJS:.o=.c)
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -- $(F90_OBJS:.o=.f90) $(C_OBJS:.o=.c) > Makefile.depend

include Makefile.depend
######################################################################
# read dependencies using sinclude (silent include)
# which continues if Makefile.depend is not readable (incl. doesn't exist!).
# Makefile.depend can be created with 'make depend'.
#-include Makefile.depend
###################################################################### 

%.o: %.mod
