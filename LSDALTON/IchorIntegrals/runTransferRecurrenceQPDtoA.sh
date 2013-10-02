ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runTransferRecurrenceQPDtoA.x runTransferRecurrenceQPDtoA.f90
./runTransferRecurrenceQPDtoA.x >& runTransferRecurrenceQPDtoAoutput.F90
ifort -g -check all -debug all -o runTransferRecurrenceQPDtoAoutput.x runTransferRecurrenceQPDtoAoutput.F90 mathfile.o IchorPresicion.o

