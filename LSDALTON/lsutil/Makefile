all: liblsutil.a

include ../Makefile.config

reord2d_2_reord.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord3d_1_reord.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord3d_2_reord.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord3d_3_reord.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord3d_1_utils_f2t.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord3d_2_utils_f2t.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord3d_3_utils_f2t.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord3d_1_utils_t2f.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord3d_2_utils_t2f.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord3d_3_utils_t2f.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord4d_1_reord.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord4d_2_reord.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord4d_3_reord.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord4d_4_reord.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord4d_1_utils_f2t.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord4d_2_utils_f2t.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord4d_3_utils_f2t.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord4d_4_utils_f2t.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord4d_1_utils_t2f.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord4d_2_utils_t2f.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord4d_3_utils_t2f.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reord4d_4_utils_t2f.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS}
reorder_frontend.F90: ./autogen/generate_man_reord.py
	python ../../LSDALTON/lsutil/autogen/generate_man_reord.py ${CPPFLAGS} "nocollapse"

# Keep the file lists sorted for easier modification.

C_OBJS = matop_csr_aux.o \
	myPAPI_set_inherit.o
F90_OBJS = \
	crayio_util.o\
	dftType.o\
	gridgeneration_memory.o\
	dft_operations.o\
	dftMem.o\
	ddynType.o\
	background_buffer.o\
	BasisinfoType.o\
	Basisinfo_operations.o\
        common_utilities.o \
        file-operations.o \
        grid_utilities.o \
	gpu_devices.o\
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
	lstensor_operations.o\
        reorder_frontend.o\
	reord2d_2_reord.o\
	reord3d_1_reord.o\
	reord3d_2_reord.o\
	reord3d_3_reord.o\
	reord4d_1_reord.o\
	reord4d_2_reord.o\
	reord4d_3_reord.o\
	reord4d_4_reord.o\
	reord2d_1_utils_f2t.o\
	reord2d_2_utils_f2t.o\
	reord2d_1_utils_t2f.o\
	reord2d_2_utils_t2f.o\
	reord3d_1_utils_f2t.o\
	reord3d_2_utils_f2t.o\
	reord3d_3_utils_f2t.o\
	reord3d_1_utils_t2f.o\
	reord3d_2_utils_t2f.o\
	reord3d_3_utils_t2f.o\
	reord4d_1_utils_f2t.o\
	reord4d_2_utils_f2t.o\
	reord4d_3_utils_f2t.o\
	reord4d_4_utils_f2t.o\
	reord4d_1_utils_t2f.o\
	reord4d_2_utils_t2f.o\
	reord4d_3_utils_t2f.o\
	reord4d_4_utils_t2f.o\
        tensor_type_def.o\
        tensor_interface.o\
        lspdm_slave.o\
        tensor_basic.o\
        lspdm_basic.o\
	lstensorType.o\
	lstensor_Mem.o\
	LSoptType.o\
	par_mod.o\
	pseudoinverse.o\
	ProfileType.o\
	lsmpiType.o\
	lsmpi-operations.o\
	matop_pdm.o\
	matop_scalapack.o\
	matop_dense.o\
	matop_csr.o\
	lsmatop_dense.o\
	matop_dense_unrest.o\
	mat-operations-aux.o\
	mat-operations-essential.o\
	matrix_utilities.o\
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
        MemoryLeakTool.o \
	Time.o\
	TestMPI.o\
	ks-settings.o\
	dec_typedef.o\
	dec_workarounds.o\
	SphCartMatrices.o\
	GCtrans.o\
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
	lspdm_tensor_operations.o\
	tensor_algebra_dil.o\
	tensor_interface.o\
	papi.o
C_OBJS = crayio.o



liblsutil.a: $(F90_OBJS) $(C_OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD reord*.F90


Makefile.depend depend: $(F90_OBJS:.o=.F90) $(C_OBJS:.o=.c)
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -- $(F90_OBJS:.o=.F90) $(C_OBJS:.o=.c) > Makefile.depend

include Makefile.depend
######################################################################
# read dependencies using sinclude (silent include)
# which continues if Makefile.depend is not readable (incl. doesn't exist!).
# Makefile.depend can be created with 'make depend'.
#-include Makefile.depend
###################################################################### 

%.o: %.mod
