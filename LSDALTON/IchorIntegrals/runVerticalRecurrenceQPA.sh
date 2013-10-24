ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runVerticalRecurrenceQPA.x runVerticalRecurrenceQPA.f90
./runVerticalRecurrenceQPA.x >& runVerticalRecurrenceQPAoutput.F90
ifort -g -check all -debug all -o runVerticalRecurrenceQPAoutput.x runVerticalRecurrenceQPAoutput.F90 mathfile.o IchorPresicion.o
#echo 'doing vec reort'
#ifort -O3 -xHost -vec-report4 -o runVerticalRecurrenceQPAoutput.x runVerticalRecurrenceQPAoutput.F90 mathfile.o IchorPresicion.o

