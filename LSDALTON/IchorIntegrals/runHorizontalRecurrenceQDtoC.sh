ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runHorizontalRecurrenceQDtoC.x runHorizontalRecurrenceQDtoC.f90
./runHorizontalRecurrenceQDtoC.x >& runHorizontalRecurrenceQDtoCoutput.F90
ifort -g -check all -debug all -openmp -o runHorizontalRecurrenceQDtoCoutput.x runHorizontalRecurrenceQDtoCoutput.F90 mathfile.o IchorPresicion.o
cp runHorizontalRecurrenceQDtoCoutput.F90 AGC_HorizontalRecurrenceQDtoC.F90
