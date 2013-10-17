ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runVerticalRecurrenceQPB.x runVerticalRecurrenceQPB.f90
./runVerticalRecurrenceQPB.x >& runVerticalRecurrenceQPBoutput.F90
ifort -g -check all -debug all -o runVerticalRecurrenceQPBoutput.x runVerticalRecurrenceQPBoutput.F90 mathfile.o IchorPresicion.o

