ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runVerticalRecurrenceQPD.x runVerticalRecurrenceQPD.f90
./runVerticalRecurrenceQPD.x >& runVerticalRecurrenceQPDoutput.F90
ifort -g -check all -debug all -o runVerticalRecurrenceQPDoutput.x runVerticalRecurrenceQPDoutput.F90 mathfile.o IchorPresicion.o

