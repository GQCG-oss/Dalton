ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -c -g -check all -debug all Strings.f90
ifort -g -check all -debug all -o runTransferRecurrenceQtoP.x runTransferRecurrenceQtoP.f90 Strings.o
./runTransferRecurrenceQtoP.x #>& runNewTransferRecurrenceQPAtoCoutput.F90

ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPCtoAGenoutput.F90 mathfile.o IchorPresicion.o     
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPCtoBGenoutput.F90 mathfile.o IchorPresicion.o     
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPDtoAGenoutput.F90 mathfile.o IchorPresicion.o          
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPDtoBGenoutput.F90 mathfile.o IchorPresicion.o          
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPCtoASegQoutput.F90 mathfile.o IchorPresicion.o         
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPCtoBSegQoutput.F90 mathfile.o IchorPresicion.o         
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPDtoASegQoutput.F90 mathfile.o IchorPresicion.o         
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPDtoBSegQoutput.F90 mathfile.o IchorPresicion.o         
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPCtoASegPoutput.F90 mathfile.o IchorPresicion.o         
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPCtoBSegPoutput.F90 mathfile.o IchorPresicion.o         
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPDtoASegPoutput.F90 mathfile.o IchorPresicion.o         
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPDtoBSegPoutput.F90 mathfile.o IchorPresicion.o         
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPCtoASegoutput.F90 mathfile.o IchorPresicion.o          
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPCtoBSegoutput.F90 mathfile.o IchorPresicion.o          
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPDtoASegoutput.F90 mathfile.o IchorPresicion.o          
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPDtoBSegoutput.F90 mathfile.o IchorPresicion.o          
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPCtoASeg1Primoutput.F90 mathfile.o IchorPresicion.o     
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPCtoBSeg1Primoutput.F90 mathfile.o IchorPresicion.o     
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPDtoASeg1Primoutput.F90 mathfile.o IchorPresicion.o     
#ifort -g -check all -debug all -o a.x runNewTransferRecurrenceQPDtoBSeg1Primoutput.F90 mathfile.o IchorPresicion.o     

