objects = ls_precision.o \
	common_utilities.o \
	memory.o \
	matrix_module.o \
	matop_dense.o \
	matop_dense_unrest.o \
	matop_csr.o \
	matop_csr_aux.o \
	mat-operations.o \
	unittest_utilities.o \
	unittest_matrix_dense.o \
	unittest_matrix_reference.o \
	unittest_matrix.o \
	unittest_driver.o \
	linpack.o \
	eispack.o \
	gp_dlapack.o \
	gp_zlapack.o \
	gp_dblas3.o \
	gp_dblas2.o \
	gp_dblas1.o \
	gp_zblas.o \
	lsmpi.o
#
FC = ifort
CC = icc
FLAGS = -fpp -O3 -DUNITTEST -DVAR_IFC -DSYS_LINUX
CFLAGS = -O3
#
unittest: $(objects) 
	$(FC) -o unittest $(FLAGS) $(objects)	 
#
unittest_utilities.o: unittest_utilities.F90
	$(FC) -c $(FLAGS) unittest_utilities.F90
unittest_matrix_reference.o: unittest_matrix_reference.F90
	$(FC) -c $(FLAGS) unittest_matrix_reference.F90
unittest_matrix_dense.o: unittest_matrix_dense.F90
	$(FC) -c $(FLAGS) unittest_matrix_dense.F90
unittest_matrix.o: unittest_matrix.F90
	$(FC) -c $(FLAGS) unittest_matrix.F90
unittest_driver.o: unittest_driver.F90
	$(FC) -c $(FLAGS) unittest_driver.F90
mat-operations.o: ../lsutil/mat-operations.F90
	$(FC) -c $(FLAGS) ../lsutil/mat-operations.F90
matrix_module.o: ../lsutil/matrix_module.F90
	$(FC) -c $(FLAGS) ../lsutil/matrix_module.F90
matop_dense.o: ../lsutil/matop_dense.F90
	$(FC) -c $(FLAGS) ../lsutil/matop_dense.F90
matop_dense_unrest.o: ../lsutil/matop_dense_unrest.F90
	$(FC) -c $(FLAGS) ../lsutil/matop_dense_unrest.F90
matop_csr.o: ../lsutil/matop_csr.F90
	$(FC) -c $(FLAGS) ../lsutil/matop_csr.F90
matop_csr_aux.o: ../lsutil/matop_csr_aux.c
	$(CC) -c $(CFLAGS) ../lsutil/matop_csr_aux.c
ls_precision.o: ../lsutil/ls_precision.F90
	$(FC) -c $(FLAGS) ../lsutil/ls_precision.F90
memory.o: ../lsutil/memory.F90
	$(FC) -c $(FLAGS) ../lsutil/memory.F90
common_utilities.o: ../lsutil/common_utilities.F90
	$(FC) -c $(FLAGS) ../lsutil/common_utilities.F90
linpack.o: ../pdpack/linpack.F
	$(FC) -c $(FLAGS) ../pdpack/linpack.F
eispack.o: ../pdpack/eispack.F
	$(FC) -c $(FLAGS) ../pdpack/eispack.F
gp_dlapack.o: ../pdpack/gp_dlapack.F
	$(FC) -c $(FLAGS) ../pdpack/gp_dlapack.F
gp_zlapack.o: ../pdpack/gp_zlapack.F
	$(FC) -c $(FLAGS) ../pdpack/gp_zlapack.F
gp_dblas3.o: ../pdpack/gp_dblas3.F
	$(FC) -c $(FLAGS) ../pdpack/gp_dblas3.F
gp_dblas2.o: ../pdpack/gp_dblas2.F
	$(FC) -c $(FLAGS) ../pdpack/gp_dblas2.F
gp_dblas1.o: ../pdpack/gp_dblas1.F
	$(FC) -c $(FLAGS) ../pdpack/gp_dblas1.F
gp_zblas.o: ../pdpack/gp_zblas.F
	$(FC) -c $(FLAGS) ../pdpack/gp_zblas.F
lsmpi.o: ../lsutil/lsmpi.F90
	$(FC) -c $(FLAGS) ../lsutil/lsmpi.F90
%.o: %.F90
	$(FC) -c $(FLAGS) $<
clean:
	rm *.mod *.o 
