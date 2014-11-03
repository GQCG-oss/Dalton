ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -c -g -check all -debug all Strings.f90
ifort -g -check all -debug all -o runTransferRecurrencePtoQ.x runTransferRecurrencePtoQ.f90 Strings.o
echo "execute runTransferRecurrencePtoQ.x"
./runTransferRecurrencePtoQ.x #>& runNewTransferRecurrenceQPAtoCoutput.F90
#ifort -g -check all -debug all -o runNewTransferRecurrenceQPAtoCoutput.x runNewTransferRecurrenceQPAtoCoutput.F90 mathfile.o IchorPresicion.o

echo "compile the files"
ifort -c -g -check all -debug all AGC_CPU_TransferRecurrenceParam.F90
ifort -c -g -check all -debug all AGC_GPU_TransferRecurrenceParam.F90

ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoCGen1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o    
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoCGen2.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o         
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceBtoCGen1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o         
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoDGen1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o  
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoDGen2.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o  
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceBtoDGen1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o         
     
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoCSegQ1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o    
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoCSegQ2.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o    
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoCSegQ3.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o    
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceBtoCSegQ1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoDSegQ1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoDSegQ2.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceBtoDSegQ1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o    
         
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoCSegP1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoCSegP2.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoCSegP3.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceBtoCSegP1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoDSegP1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoDSegP2.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceBtoDSegP1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o 

ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoCSeg1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o  
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoCSeg2.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o  
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoCSeg3.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o  
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceBtoCSeg1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o  
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoDSeg1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o  
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoDSeg2.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o  
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceBtoDSeg1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o  

ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoCSeg1Prim1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o         
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoCSeg1Prim2.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o         
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceBtoCSeg1Prim1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o         
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoDSeg1Prim1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o         
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceAtoDSeg1Prim2.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o         
ifort -g -check all -openmp -debug all -o a.x AGC_CPU_TransferRecurrenceBtoDSeg1Prim1.F90 mathfile.o IchorPresicion.o AGC_CPU_TransferRecurrenceParam.o         


echo "compile the GPU files"

ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoCGen1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o         
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoCGen2.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o         
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceBtoCGen1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o         
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoDGen1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o  
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoDGen2.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o  
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceBtoDGen1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o  
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoCSegQ1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o    
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoCSegQ2.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o    
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoCSegQ3.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o    
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceBtoCSegQ1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoDSegQ1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoDSegQ2.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceBtoDSegQ1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoCSegP1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoCSegP2.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoCSegP3.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceBtoCSegP1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoDSegP1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoDSegP2.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceBtoDSegP1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o 
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoCSeg1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o  
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoCSeg2.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o  
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoCSeg3.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o  
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceBtoCSeg1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o  
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoDSeg1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o  
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoDSeg2.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o  
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceBtoDSeg1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o  
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoCSeg1Prim1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoCSeg1Prim2.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceBtoCSeg1Prim1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoDSeg1Prim1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceAtoDSeg1Prim2.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o
ifort -g -check all -openmp -debug all -o a.x AGC_GPU_TransferRecurrenceBtoDSeg1Prim1.F90 mathfile.o IchorPresicion.o AGC_GPU_TransferRecurrenceParam.o



