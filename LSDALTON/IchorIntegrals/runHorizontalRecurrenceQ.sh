ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runHorizontalRecurrenceQ.x runHorizontalRecurrenceQ.f90
./runHorizontalRecurrenceQ.x >& runHorizontalRecurrenceQoutput.F90
ifort -g -check all -debug all -o runHorizontalRecurrenceQoutput.x runHorizontalRecurrenceQoutput.F90 mathfile.o IchorPresicion.o


