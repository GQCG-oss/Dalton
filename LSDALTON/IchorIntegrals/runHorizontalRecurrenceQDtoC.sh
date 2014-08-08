ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runHorizontalRecurrenceQDtoC.x runHorizontalRecurrenceQDtoC.f90
./runHorizontalRecurrenceQDtoC.x 
ifort -g -check all -debug all -openmp -o runHorizontalRecurrenceQDtoCoutput.x runHorizontalRecurrenceCPURHSModDtoC.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runHorizontalRecurrenceQDtoCoutput.x runHorizontalRecurrenceGPURHSModDtoC.F90 mathfile.o IchorPresicion.o

cp runHorizontalRecurrenceCPURHSModDtoC.F90 AGC_CPU_HorizontalRecurrenceQDtoC.F90
cp runHorizontalRecurrenceGPURHSModDtoC.F90 AGC_GPU_HorizontalRecurrenceQDtoC.F90


