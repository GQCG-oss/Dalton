ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runVerticalRecurrenceQPC.x runVerticalRecurrenceQPC.f90
./runVerticalRecurrenceQPC.x >& runVerticalRecurrenceQPCoutput.F90
ifort -g -check all -debug all -o runVerticalRecurrenceQPCoutput.x runVerticalRecurrenceQPCoutput.F90 mathfile.o IchorPresicion.o

