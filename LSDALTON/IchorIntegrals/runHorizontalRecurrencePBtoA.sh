ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runHorizontalRecurrencePBtoA.x runHorizontalRecurrencePBtoA.f90
./runHorizontalRecurrencePBtoA.x >& runHorizontalRecurrencePBtoAoutput.F90
ifort -g -check all -debug all -o runHorizontalRecurrencePBtoAoutput.x runHorizontalRecurrencePBtoAoutput.F90 mathfile.o IchorPresicion.o


