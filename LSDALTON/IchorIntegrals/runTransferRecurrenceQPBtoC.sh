ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runTransferRecurrenceQPBtoC.x runTransferRecurrenceQPBtoC.f90
./runTransferRecurrenceQPBtoC.x >& runTransferRecurrenceQPBtoCoutput.F90
ifort -g -check all -debug all -o runTransferRecurrenceQPBtoCoutput.x runTransferRecurrenceQPBtoCoutput.F90 mathfile.o IchorPresicion.o

