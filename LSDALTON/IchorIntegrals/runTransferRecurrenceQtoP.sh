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
