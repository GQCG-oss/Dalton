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
unittest_utilities.o: unittest_utilities.f90
	$(FC) -c $(FLAGS) unittest_utilities.f90
unittest_matrix_reference.o: unittest_matrix_reference.f90
	$(FC) -c $(FLAGS) unittest_matrix_reference.f90
unittest_matrix_dense.o: unittest_matrix_dense.f90
	$(FC) -c $(FLAGS) unittest_matrix_dense.f90
unittest_matrix.o: unittest_matrix.f90
	$(FC) -c $(FLAGS) unittest_matrix.f90
unittest_driver.o: unittest_driver.f90
	$(FC) -c $(FLAGS) unittest_driver.f90
mat-operations.o: ../lsutil/mat-operations.f90
	$(FC) -c $(FLAGS) ../lsutil/mat-operations.f90
matrix_module.o: ../lsutil/matrix_module.f90
	$(FC) -c $(FLAGS) ../lsutil/matrix_module.f90
matop_dense.o: ../lsutil/matop_dense.f90
	$(FC) -c $(FLAGS) ../lsutil/matop_dense.f90
matop_dense_unrest.o: ../lsutil/matop_dense_unrest.f90
	$(FC) -c $(FLAGS) ../lsutil/matop_dense_unrest.f90
matop_csr.o: ../lsutil/matop_csr.f90
	$(FC) -c $(FLAGS) ../lsutil/matop_csr.f90
matop_csr_aux.o: ../lsutil/matop_csr_aux.c
	$(CC) -c $(CFLAGS) ../lsutil/matop_csr_aux.c
ls_precision.o: ../lsutil/ls_precision.f90
	$(FC) -c $(FLAGS) ../lsutil/ls_precision.f90
memory.o: ../lsutil/memory.f90
	$(FC) -c $(FLAGS) ../lsutil/memory.f90
common_utilities.o: ../lsutil/common_utilities.f90
	$(FC) -c $(FLAGS) ../lsutil/common_utilities.f90
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
lsmpi.o: ../lsutil/lsmpi.f90
	$(FC) -c $(FLAGS) ../lsutil/lsmpi.f90
%.o: %.f90
	$(FC) -c $(FLAGS) $<
clean:
	rm *.mod *.o 
