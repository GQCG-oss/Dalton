ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runTransferRecurrenceQPDtoB.x runTransferRecurrenceQPDtoB.f90
./runTransferRecurrenceQPDtoB.x >& runTransferRecurrenceQPDtoBoutput.F90
ifort -g -check all -debug all -o runTransferRecurrenceQPDtoBoutput.x runTransferRecurrenceQPDtoBoutput.F90 mathfile.o IchorPresicion.o

