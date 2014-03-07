ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runHorizontalRecurrenceQCtoD.x runHorizontalRecurrenceQCtoD.f90
./runHorizontalRecurrenceQCtoD.x >& runHorizontalRecurrenceQCtoDoutput.F90
ifort -g -check all -debug all -openmp -o runHorizontalRecurrenceQCtoDoutput.x runHorizontalRecurrenceQCtoDoutput.F90 mathfile.o IchorPresicion.o
cp runHorizontalRecurrenceQCtoDoutput.F90 AGC_HorizontalRecurrenceQCtoD.F90


