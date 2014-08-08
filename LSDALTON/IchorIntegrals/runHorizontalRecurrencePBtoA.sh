ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runHorizontalRecurrencePBtoA.x runHorizontalRecurrencePBtoA.f90
./runHorizontalRecurrencePBtoA.x 
ifort -g -check all -debug all -openmp -o runHorizontalRecurrencePBtoAoutput.x runHorizontalRecurrenceCPULHSModBtoA.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runHorizontalRecurrencePBtoAoutput.x runHorizontalRecurrenceGPULHSModBtoA.F90 mathfile.o IchorPresicion.o

cp runHorizontalRecurrenceCPULHSModBtoA.F90 AGC_CPU_HorizontalRecurrencePBtoA.F90
cp runHorizontalRecurrenceGPULHSModBtoA.F90 AGC_GPU_HorizontalRecurrencePBtoA.F90
