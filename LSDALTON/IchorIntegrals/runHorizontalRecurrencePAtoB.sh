ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runHorizontalRecurrencePAtoB.x runHorizontalRecurrencePAtoB.f90
./runHorizontalRecurrencePAtoB.x >& runHorizontalRecurrencePAtoBoutput.F90
ifort -g -check all -debug all -o runHorizontalRecurrencePAtoBoutput.x runHorizontalRecurrencePAtoBoutput.F90 mathfile.o IchorPresicion.o

cp runHorizontalRecurrencePAtoBoutput.F90 AGC_HorizontalRecurrencePAtoB.F90
