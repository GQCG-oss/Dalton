ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runVerticalRecurrenceQP.x runVerticalRecurrenceQP.f90
./runVerticalRecurrenceQP.x >& runVerticalRecurrenceQPoutput.F90
ifort -g -check all -debug all -o runVerticalRecurrenceQPoutput.x runVerticalRecurrenceQPoutput.F90 mathfile.o IchorPresicion.o

