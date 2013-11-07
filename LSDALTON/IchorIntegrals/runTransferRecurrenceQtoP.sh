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

cp runNewTransferRecurrenceQPCtoAGenoutput.F90 AGC_TransferRecurrenceCtoAGen.F90 
cp runNewTransferRecurrenceQPCtoASeg1Primoutput.F90 AGC_TransferRecurrenceCtoASeg1Prim.F90 
cp runNewTransferRecurrenceQPCtoASegoutput.F90 AGC_TransferRecurrenceCtoASeg.F90 
cp runNewTransferRecurrenceQPCtoASegPoutput.F90 AGC_TransferRecurrenceCtoASegP.F90 
cp runNewTransferRecurrenceQPCtoASegQoutput.F90 AGC_TransferRecurrenceCtoASegQ.F90 

cp runNewTransferRecurrenceQPCtoBGenoutput.F90 AGC_TransferRecurrenceCtoBGen.F90 
cp runNewTransferRecurrenceQPCtoBSeg1Primoutput.F90 AGC_TransferRecurrenceCtoBSeg1Prim.F90 
cp runNewTransferRecurrenceQPCtoBSegoutput.F90 AGC_TransferRecurrenceCtoBSeg.F90 
cp runNewTransferRecurrenceQPCtoBSegPoutput.F90 AGC_TransferRecurrenceCtoBSegP.F90 
cp runNewTransferRecurrenceQPCtoBSegQoutput.F90 AGC_TransferRecurrenceCtoBSegQ.F90 

cp runNewTransferRecurrenceQPDtoAGenoutput.F90 AGC_TransferRecurrenceDtoAGen.F90 
cp runNewTransferRecurrenceQPDtoASeg1Primoutput.F90 AGC_TransferRecurrenceDtoASeg1Prim.F90 
cp runNewTransferRecurrenceQPDtoASegoutput.F90 AGC_TransferRecurrenceDtoASeg.F90 
cp runNewTransferRecurrenceQPDtoASegPoutput.F90 AGC_TransferRecurrenceDtoASegP.F90 
cp runNewTransferRecurrenceQPDtoASegQoutput.F90 AGC_TransferRecurrenceDtoASegQ.F90 

cp runNewTransferRecurrenceQPDtoBGenoutput.F90 AGC_TransferRecurrenceDtoBGen.F90 
cp runNewTransferRecurrenceQPDtoBSeg1Primoutput.F90 AGC_TransferRecurrenceDtoBSeg1Prim.F90 
cp runNewTransferRecurrenceQPDtoBSegoutput.F90 AGC_TransferRecurrenceDtoBSeg.F90 
cp runNewTransferRecurrenceQPDtoBSegPoutput.F90 AGC_TransferRecurrenceDtoBSegP.F90 
cp runNewTransferRecurrenceQPDtoBSegQoutput.F90 AGC_TransferRecurrenceDtoBSegQ.F90 
