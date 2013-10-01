ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runTransferRecurrenceQP.x runTransferRecurrenceQP.f90
./runTransferRecurrenceQP.x >& runTransferRecurrenceQPoutput.F90
ifort -g -check all -debug all -o runTransferRecurrenceQPoutput.x runTransferRecurrenceQPoutput.F90 mathfile.o IchorPresicion.o

