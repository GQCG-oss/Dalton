ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -c -g -check all -debug all Strings.f90
ifort -g -check all -debug all -o runTransferRecurrencePtoQ.x runTransferRecurrencePtoQ.f90 Strings.o
./runTransferRecurrencePtoQ.x #>& runNewTransferRecurrenceQPAtoCoutput.F90
#ifort -g -check all -debug all -o runNewTransferRecurrenceQPAtoCoutput.x runNewTransferRecurrenceQPAtoCoutput.F90 mathfile.o IchorPresicion.o

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


cp runNewTransferRecurrenceQPAtoCGenoutput.F90 AGC_TransferRecurrenceAtoCGen.F90 
cp runNewTransferRecurrenceQPAtoCSeg1Primoutput.F90 AGC_TransferRecurrenceAtoCSeg1Prim.F90 
cp runNewTransferRecurrenceQPAtoCSegoutput.F90 AGC_TransferRecurrenceAtoCSeg.F90 
cp runNewTransferRecurrenceQPAtoCSegPoutput.F90 AGC_TransferRecurrenceAtoCSegP.F90 
cp runNewTransferRecurrenceQPAtoCSegQoutput.F90 AGC_TransferRecurrenceAtoCSegQ.F90 

cp runNewTransferRecurrenceQPAtoDGenoutput.F90 AGC_TransferRecurrenceAtoDGen.F90 
cp runNewTransferRecurrenceQPAtoDSeg1Primoutput.F90 AGC_TransferRecurrenceAtoDSeg1Prim.F90 
cp runNewTransferRecurrenceQPAtoDSegoutput.F90 AGC_TransferRecurrenceAtoDSeg.F90 
cp runNewTransferRecurrenceQPAtoDSegPoutput.F90 AGC_TransferRecurrenceAtoDSegP.F90 
cp runNewTransferRecurrenceQPAtoDSegQoutput.F90 AGC_TransferRecurrenceAtoDSegQ.F90 

cp runNewTransferRecurrenceQPBtoCGenoutput.F90 AGC_TransferRecurrenceBtoCGen.F90 
cp runNewTransferRecurrenceQPBtoCSeg1Primoutput.F90 AGC_TransferRecurrenceBtoCSeg1Prim.F90 
cp runNewTransferRecurrenceQPBtoCSegoutput.F90 AGC_TransferRecurrenceBtoCSeg.F90 
cp runNewTransferRecurrenceQPBtoCSegPoutput.F90 AGC_TransferRecurrenceBtoCSegP.F90 
cp runNewTransferRecurrenceQPBtoCSegQoutput.F90 AGC_TransferRecurrenceBtoCSegQ.F90 

cp runNewTransferRecurrenceQPBtoDGenoutput.F90 AGC_TransferRecurrenceBtoDGen.F90 
cp runNewTransferRecurrenceQPBtoDSeg1Primoutput.F90 AGC_TransferRecurrenceBtoDSeg1Prim.F90 
cp runNewTransferRecurrenceQPBtoDSegoutput.F90 AGC_TransferRecurrenceBtoDSeg.F90 
cp runNewTransferRecurrenceQPBtoDSegPoutput.F90 AGC_TransferRecurrenceBtoDSegP.F90 
cp runNewTransferRecurrenceQPBtoDSegQoutput.F90 AGC_TransferRecurrenceBtoDSegQ.F90 
