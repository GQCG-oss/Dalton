ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runTransferRecurrenceQPBtoD.x runTransferRecurrenceQPBtoD.f90
./runTransferRecurrenceQPBtoD.x >& runTransferRecurrenceQPBtoDoutput.F90
ifort -g -check all -debug all -o runTransferRecurrenceQPBtoDoutput.x runTransferRecurrenceQPBtoDoutput.F90 mathfile.o IchorPresicion.o

