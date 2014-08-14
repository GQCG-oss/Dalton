ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -c -g -check all -debug all Strings.f90
ifort -g -check all -debug all -o runHorizontalRecurrencePAtoB.x runHorizontalRecurrencePAtoB.f90 Strings.o
./runHorizontalRecurrencePAtoB.x 
ifort -g -check all -debug all -openmp -o runHorizontalRecurrencePAtoBoutput.x runHorizontalRecurrenceCPULHSModAtoB.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runHorizontalRecurrencePAtoBoutput.x runHorizontalRecurrenceGPULHSModAtoB.F90 mathfile.o IchorPresicion.o

cp runHorizontalRecurrenceCPULHSModAtoB.F90 AGC_CPU_HorizontalRecurrencePAtoB.F90
cp runHorizontalRecurrenceGPULHSModAtoB.F90 AGC_GPU_HorizontalRecurrencePAtoB.F90
