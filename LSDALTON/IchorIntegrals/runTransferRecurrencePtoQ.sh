ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -c -g -check all -debug all Strings.f90
ifort -g -check all -debug all -o runTransferRecurrencePtoQ.x runTransferRecurrencePtoQ.f90 Strings.o
./runTransferRecurrencePtoQ.x #>& runNewTransferRecurrenceQPAtoCoutput.F90
ifort -g -check all -debug all -o runNewTransferRecurrenceQPAtoCoutput.x runNewTransferRecurrenceQPAtoCoutput.F90 mathfile.o IchorPresicion.o

ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPAtoCGenoutput.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPBtoCGenoutput.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPAtoDGenoutput.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPBtoDGenoutput.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPAtoCSegQoutput.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPBtoCSegQoutput.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPAtoDSegQoutput.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPBtoDSegQoutput.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPAtoCSegPoutput.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPBtoCSegPoutput.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPAtoDSegPoutput.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPBtoDSegPoutput.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPAtoCSegoutput.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPBtoCSegoutput.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPAtoDSegoutput.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPBtoDSegoutput.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPAtoCSeg1Primoutput.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPBtoCSeg1Primoutput.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPAtoDSeg1Primoutput.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPBtoDSeg1Primoutput.F90 mathfile.o IchorPresicion.o     

