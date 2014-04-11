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

echo "compile the GPU files"

ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoCGenoutputGPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoCGenoutputGPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoDGenoutputGPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoDGenoutputGPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoCSegQoutputGPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoCSegQoutputGPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoDSegQoutputGPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoDSegQoutputGPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoCSegPoutputGPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoCSegPoutputGPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoDSegPoutputGPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoDSegPoutputGPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoCSegoutputGPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoCSegoutputGPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoDSegoutputGPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoDSegoutputGPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoCSeg1PrimoutputGPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoCSeg1PrimoutputGPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPAtoDSeg1PrimoutputGPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPBtoDSeg1PrimoutputGPU.F90 mathfile.o IchorPresicion.o     


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

#GPU
cp runNewTransferRecurrenceQPAtoCGenoutputGPU.F90 AGC_GPU_TransferRecurrenceAtoCGen.F90 
cp runNewTransferRecurrenceQPAtoCSeg1PrimoutputGPU.F90 AGC_GPU_TransferRecurrenceAtoCSeg1Prim.F90 
cp runNewTransferRecurrenceQPAtoCSegoutputGPU.F90 AGC_GPU_TransferRecurrenceAtoCSeg.F90 
cp runNewTransferRecurrenceQPAtoCSegPoutputGPU.F90 AGC_GPU_TransferRecurrenceAtoCSegP.F90 
cp runNewTransferRecurrenceQPAtoCSegQoutputGPU.F90 AGC_GPU_TransferRecurrenceAtoCSegQ.F90 

cp runNewTransferRecurrenceQPAtoDGenoutputGPU.F90 AGC_GPU_TransferRecurrenceAtoDGen.F90 
cp runNewTransferRecurrenceQPAtoDSeg1PrimoutputGPU.F90 AGC_GPU_TransferRecurrenceAtoDSeg1Prim.F90 
cp runNewTransferRecurrenceQPAtoDSegoutputGPU.F90 AGC_GPU_TransferRecurrenceAtoDSeg.F90 
cp runNewTransferRecurrenceQPAtoDSegPoutputGPU.F90 AGC_GPU_TransferRecurrenceAtoDSegP.F90 
cp runNewTransferRecurrenceQPAtoDSegQoutputGPU.F90 AGC_GPU_TransferRecurrenceAtoDSegQ.F90 

cp runNewTransferRecurrenceQPBtoCGenoutputGPU.F90 AGC_GPU_TransferRecurrenceBtoCGen.F90 
cp runNewTransferRecurrenceQPBtoCSeg1PrimoutputGPU.F90 AGC_GPU_TransferRecurrenceBtoCSeg1Prim.F90 
cp runNewTransferRecurrenceQPBtoCSegoutputGPU.F90 AGC_GPU_TransferRecurrenceBtoCSeg.F90 
cp runNewTransferRecurrenceQPBtoCSegPoutputGPU.F90 AGC_GPU_TransferRecurrenceBtoCSegP.F90 
cp runNewTransferRecurrenceQPBtoCSegQoutputGPU.F90 AGC_GPU_TransferRecurrenceBtoCSegQ.F90 

cp runNewTransferRecurrenceQPBtoDGenoutputGPU.F90 AGC_GPU_TransferRecurrenceBtoDGen.F90 
cp runNewTransferRecurrenceQPBtoDSeg1PrimoutputGPU.F90 AGC_GPU_TransferRecurrenceBtoDSeg1Prim.F90 
cp runNewTransferRecurrenceQPBtoDSegoutputGPU.F90 AGC_GPU_TransferRecurrenceBtoDSeg.F90 
cp runNewTransferRecurrenceQPBtoDSegPoutputGPU.F90 AGC_GPU_TransferRecurrenceBtoDSegP.F90 
cp runNewTransferRecurrenceQPBtoDSegQoutputGPU.F90 AGC_GPU_TransferRecurrenceBtoDSegQ.F90 

