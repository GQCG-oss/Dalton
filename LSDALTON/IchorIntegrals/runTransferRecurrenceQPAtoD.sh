ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runTransferRecurrenceQPAtoD.x runTransferRecurrenceQPAtoD.f90
./runTransferRecurrenceQPAtoD.x >& runTransferRecurrenceQPAtoDoutput.F90
ifort -g -check all -debug all -o runTransferRecurrenceQPAtoDoutput.x runTransferRecurrenceQPAtoDoutput.F90 mathfile.o IchorPresicion.o

