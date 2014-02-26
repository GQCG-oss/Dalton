ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -c -g -check all -debug all Strings.f90
ifort -g -check all -debug all -o runTransferRecurrencePtoQ.x runTransferRecurrencePtoQ.f90 Strings.o
echo "execute runTransferRecurrencePtoQ.x"
./runTransferRecurrencePtoQ.x #>& runNewTransferRecurrenceQPAtoCoutput.F90
#ifort -g -check all -debug all -o runNewTransferRecurrenceQPAtoCoutput.x runNewTransferRecurrenceQPAtoCoutput.F90 mathfile.o IchorPresicion.o

echo "compile the files"

ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoCGenoutputCPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoCGenoutputCPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoDGenoutputCPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoDGenoutputCPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoCSegQoutputCPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoCSegQoutputCPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoDSegQoutputCPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoDSegQoutputCPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoCSegPoutputCPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoCSegPoutputCPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoDSegPoutputCPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoDSegPoutputCPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoCSegoutputCPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoCSegoutputCPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoDSegoutputCPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoDSegoutputCPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoCSeg1PrimoutputCPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoCSeg1PrimoutputCPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoDSeg1PrimoutputCPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoDSeg1PrimoutputCPU.F90 mathfile.o IchorPresicion.o     


cp runNewTransferRecurrenceQPAtoCGenoutputCPU.F90 AGC_CPU_TransferRecurrenceAtoCGen.F90 
cp runNewTransferRecurrenceQPAtoCSeg1PrimoutputCPU.F90 AGC_CPU_TransferRecurrenceAtoCSeg1Prim.F90 
cp runNewTransferRecurrenceQPAtoCSegoutputCPU.F90 AGC_CPU_TransferRecurrenceAtoCSeg.F90 
cp runNewTransferRecurrenceQPAtoCSegPoutputCPU.F90 AGC_CPU_TransferRecurrenceAtoCSegP.F90 
cp runNewTransferRecurrenceQPAtoCSegQoutputCPU.F90 AGC_CPU_TransferRecurrenceAtoCSegQ.F90 

cp runNewTransferRecurrenceQPAtoDGenoutputCPU.F90 AGC_CPU_TransferRecurrenceAtoDGen.F90 
cp runNewTransferRecurrenceQPAtoDSeg1PrimoutputCPU.F90 AGC_CPU_TransferRecurrenceAtoDSeg1Prim.F90 
cp runNewTransferRecurrenceQPAtoDSegoutputCPU.F90 AGC_CPU_TransferRecurrenceAtoDSeg.F90 
cp runNewTransferRecurrenceQPAtoDSegPoutputCPU.F90 AGC_CPU_TransferRecurrenceAtoDSegP.F90 
cp runNewTransferRecurrenceQPAtoDSegQoutputCPU.F90 AGC_CPU_TransferRecurrenceAtoDSegQ.F90 

cp runNewTransferRecurrenceQPBtoCGenoutputCPU.F90 AGC_CPU_TransferRecurrenceBtoCGen.F90 
cp runNewTransferRecurrenceQPBtoCSeg1PrimoutputCPU.F90 AGC_CPU_TransferRecurrenceBtoCSeg1Prim.F90 
cp runNewTransferRecurrenceQPBtoCSegoutputCPU.F90 AGC_CPU_TransferRecurrenceBtoCSeg.F90 
cp runNewTransferRecurrenceQPBtoCSegPoutputCPU.F90 AGC_CPU_TransferRecurrenceBtoCSegP.F90 
cp runNewTransferRecurrenceQPBtoCSegQoutputCPU.F90 AGC_CPU_TransferRecurrenceBtoCSegQ.F90 

cp runNewTransferRecurrenceQPBtoDGenoutputCPU.F90 AGC_CPU_TransferRecurrenceBtoDGen.F90 
cp runNewTransferRecurrenceQPBtoDSeg1PrimoutputCPU.F90 AGC_CPU_TransferRecurrenceBtoDSeg1Prim.F90 
cp runNewTransferRecurrenceQPBtoDSegoutputCPU.F90 AGC_CPU_TransferRecurrenceBtoDSeg.F90 
cp runNewTransferRecurrenceQPBtoDSegPoutputCPU.F90 AGC_CPU_TransferRecurrenceBtoDSegP.F90 
cp runNewTransferRecurrenceQPBtoDSegQoutputCPU.F90 AGC_CPU_TransferRecurrenceBtoDSegQ.F90 

