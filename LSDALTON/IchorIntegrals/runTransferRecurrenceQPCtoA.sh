ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -c -g -check all -debug all runTransferRecurrenceQ.f90
ifort -g -check all -debug all -o runTransferRecurrenceQPCtoA.x runTransferRecurrenceQPCtoA.f90 runTransferRecurrenceQ.o
./runTransferRecurrenceQPCtoA.x >& runTransferRecurrenceQPCtoAoutput.F90
ifort -g -check all -debug all -o runTransferRecurrenceQPCtoAoutput.x runTransferRecurrenceQPCtoAoutput.F90 mathfile.o IchorPresicion.o

