ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runTransferRecurrenceQPAtoC.x runTransferRecurrenceQPAtoC.f90
./runTransferRecurrenceQPAtoC.x >& runTransferRecurrenceQPAtoCoutput.F90
ifort -g -check all -debug all -o runTransferRecurrenceQPAtoCoutput.x runTransferRecurrenceQPAtoCoutput.F90 mathfile.o IchorPresicion.o

