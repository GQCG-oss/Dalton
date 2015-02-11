ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -c -g -check all -debug all Strings.f90
ifort -g -check all -debug all -o runHorizontalRecurrenceQCtoD.x runHorizontalRecurrenceQCtoD.f90 Strings.o
./runHorizontalRecurrenceQCtoD.x 
ifort -g -check all -debug all -openmp -o runHorizontalRecurrenceQCtoDoutput.x runHorizontalRecurrenceCPURHSModCtoD.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runHorizontalRecurrenceQCtoDoutput.x runHorizontalRecurrenceGPURHSModCtoD.F90 mathfile.o IchorPresicion.o

cp runHorizontalRecurrenceCPURHSModCtoD.F90 AGC_CPU_HorizontalRecurrenceQCtoD.F90
cp runHorizontalRecurrenceGPURHSModCtoD.F90 AGC_GPU_HorizontalRecurrenceQCtoD.F90


