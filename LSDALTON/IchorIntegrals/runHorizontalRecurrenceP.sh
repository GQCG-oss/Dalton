ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runHorizontalRecurrenceP.x runHorizontalRecurrenceP.f90
./runHorizontalRecurrenceP.x >& runHorizontalRecurrencePoutput.F90
ifort -g -check all -debug all -o runHorizontalRecurrencePoutput.x runHorizontalRecurrencePoutput.F90 mathfile.o IchorPresicion.o


