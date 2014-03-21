ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -c -g -check all -debug all Strings.f90
ifort -g -check all -debug all -o runTransferRecurrenceQtoP.x runTransferRecurrenceQtoP.f90 Strings.o
./runTransferRecurrenceQtoP.x #>& runNewTransferRecurrenceQPAtoCoutput.F90

ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoAGenoutputCPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoBGenoutputCPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoAGenoutputCPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoBGenoutputCPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoASegQoutputCPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoBSegQoutputCPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoASegQoutputCPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoBSegQoutputCPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoASegPoutputCPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoBSegPoutputCPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoASegPoutputCPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoBSegPoutputCPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoASegoutputCPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoBSegoutputCPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoASegoutputCPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoBSegoutputCPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoASeg1PrimoutputCPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoBSeg1PrimoutputCPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoASeg1PrimoutputCPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoBSeg1PrimoutputCPU.F90 mathfile.o IchorPresicion.o     

ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoAGenoutputGPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoBGenoutputGPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoAGenoutputGPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoBGenoutputGPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoASegQoutputGPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoBSegQoutputGPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoASegQoutputGPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoBSegQoutputGPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoASegPoutputGPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoBSegPoutputGPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoASegPoutputGPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoBSegPoutputGPU.F90 mathfile.o IchorPresicion.o         
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoASegoutputGPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoBSegoutputGPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoASegoutputGPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoBSegoutputGPU.F90 mathfile.o IchorPresicion.o          
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoASeg1PrimoutputGPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPCtoBSeg1PrimoutputGPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoASeg1PrimoutputGPU.F90 mathfile.o IchorPresicion.o     
ifort -g -check all -openmp -debug all -o a.x runNewTransferRecurrenceQPDtoBSeg1PrimoutputGPU.F90 mathfile.o IchorPresicion.o     

cp runNewTransferRecurrenceQPCtoAGenoutputCPU.F90 AGC_CPU_TransferRecurrenceCtoAGen.F90 
cp runNewTransferRecurrenceQPCtoASeg1PrimoutputCPU.F90 AGC_CPU_TransferRecurrenceCtoASeg1Prim.F90 
cp runNewTransferRecurrenceQPCtoASegoutputCPU.F90 AGC_CPU_TransferRecurrenceCtoASeg.F90 
cp runNewTransferRecurrenceQPCtoASegPoutputCPU.F90 AGC_CPU_TransferRecurrenceCtoASegP.F90 
cp runNewTransferRecurrenceQPCtoASegQoutputCPU.F90 AGC_CPU_TransferRecurrenceCtoASegQ.F90 

cp runNewTransferRecurrenceQPCtoBGenoutputCPU.F90 AGC_CPU_TransferRecurrenceCtoBGen.F90 
cp runNewTransferRecurrenceQPCtoBSeg1PrimoutputCPU.F90 AGC_CPU_TransferRecurrenceCtoBSeg1Prim.F90 
cp runNewTransferRecurrenceQPCtoBSegoutputCPU.F90 AGC_CPU_TransferRecurrenceCtoBSeg.F90 
cp runNewTransferRecurrenceQPCtoBSegPoutputCPU.F90 AGC_CPU_TransferRecurrenceCtoBSegP.F90 
cp runNewTransferRecurrenceQPCtoBSegQoutputCPU.F90 AGC_CPU_TransferRecurrenceCtoBSegQ.F90 

cp runNewTransferRecurrenceQPDtoAGenoutputCPU.F90 AGC_CPU_TransferRecurrenceDtoAGen.F90 
cp runNewTransferRecurrenceQPDtoASeg1PrimoutputCPU.F90 AGC_CPU_TransferRecurrenceDtoASeg1Prim.F90 
cp runNewTransferRecurrenceQPDtoASegoutputCPU.F90 AGC_CPU_TransferRecurrenceDtoASeg.F90 
cp runNewTransferRecurrenceQPDtoASegPoutputCPU.F90 AGC_CPU_TransferRecurrenceDtoASegP.F90 
cp runNewTransferRecurrenceQPDtoASegQoutputCPU.F90 AGC_CPU_TransferRecurrenceDtoASegQ.F90 

cp runNewTransferRecurrenceQPDtoBGenoutputCPU.F90 AGC_CPU_TransferRecurrenceDtoBGen.F90 
cp runNewTransferRecurrenceQPDtoBSeg1PrimoutputCPU.F90 AGC_CPU_TransferRecurrenceDtoBSeg1Prim.F90 
cp runNewTransferRecurrenceQPDtoBSegoutputCPU.F90 AGC_CPU_TransferRecurrenceDtoBSeg.F90 
cp runNewTransferRecurrenceQPDtoBSegPoutputCPU.F90 AGC_CPU_TransferRecurrenceDtoBSegP.F90 
cp runNewTransferRecurrenceQPDtoBSegQoutputCPU.F90 AGC_CPU_TransferRecurrenceDtoBSegQ.F90 

#GPU

cp runNewTransferRecurrenceQPCtoAGenoutputGPU.F90 AGC_GPU_TransferRecurrenceCtoAGen.F90 
cp runNewTransferRecurrenceQPCtoASeg1PrimoutputGPU.F90 AGC_GPU_TransferRecurrenceCtoASeg1Prim.F90 
cp runNewTransferRecurrenceQPCtoASegoutputGPU.F90 AGC_GPU_TransferRecurrenceCtoASeg.F90 
cp runNewTransferRecurrenceQPCtoASegPoutputGPU.F90 AGC_GPU_TransferRecurrenceCtoASegP.F90 
cp runNewTransferRecurrenceQPCtoASegQoutputGPU.F90 AGC_GPU_TransferRecurrenceCtoASegQ.F90 

cp runNewTransferRecurrenceQPCtoBGenoutputGPU.F90 AGC_GPU_TransferRecurrenceCtoBGen.F90 
cp runNewTransferRecurrenceQPCtoBSeg1PrimoutputGPU.F90 AGC_GPU_TransferRecurrenceCtoBSeg1Prim.F90 
cp runNewTransferRecurrenceQPCtoBSegoutputGPU.F90 AGC_GPU_TransferRecurrenceCtoBSeg.F90 
cp runNewTransferRecurrenceQPCtoBSegPoutputGPU.F90 AGC_GPU_TransferRecurrenceCtoBSegP.F90 
cp runNewTransferRecurrenceQPCtoBSegQoutputGPU.F90 AGC_GPU_TransferRecurrenceCtoBSegQ.F90 

cp runNewTransferRecurrenceQPDtoAGenoutputGPU.F90 AGC_GPU_TransferRecurrenceDtoAGen.F90 
cp runNewTransferRecurrenceQPDtoASeg1PrimoutputGPU.F90 AGC_GPU_TransferRecurrenceDtoASeg1Prim.F90 
cp runNewTransferRecurrenceQPDtoASegoutputGPU.F90 AGC_GPU_TransferRecurrenceDtoASeg.F90 
cp runNewTransferRecurrenceQPDtoASegPoutputGPU.F90 AGC_GPU_TransferRecurrenceDtoASegP.F90 
cp runNewTransferRecurrenceQPDtoASegQoutputGPU.F90 AGC_GPU_TransferRecurrenceDtoASegQ.F90 

cp runNewTransferRecurrenceQPDtoBGenoutputGPU.F90 AGC_GPU_TransferRecurrenceDtoBGen.F90 
cp runNewTransferRecurrenceQPDtoBSeg1PrimoutputGPU.F90 AGC_GPU_TransferRecurrenceDtoBSeg1Prim.F90 
cp runNewTransferRecurrenceQPDtoBSegoutputGPU.F90 AGC_GPU_TransferRecurrenceDtoBSeg.F90 
cp runNewTransferRecurrenceQPDtoBSegPoutputGPU.F90 AGC_GPU_TransferRecurrenceDtoBSegP.F90 
cp runNewTransferRecurrenceQPDtoBSegQoutputGPU.F90 AGC_GPU_TransferRecurrenceDtoBSegQ.F90 
