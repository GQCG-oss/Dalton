ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runTransferRecurrenceQPCtoB.x runTransferRecurrenceQPCtoB.f90
./runTransferRecurrenceQPCtoB.x >& runTransferRecurrenceQPCtoBoutput.F90
ifort -g -check all -debug all -o runTransferRecurrenceQPCtoBoutput.x runTransferRecurrenceQPCtoBoutput.F90 mathfile.o IchorPresicion.o

