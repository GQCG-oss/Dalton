rm runVerticalRecurrenceQP.x 
ifort -c mathfile.f90
ifort -c IchorPresicion.F90
ifort -c Strings.f90
ifort -o runVerticalRecurrenceQP.x runVerticalRecurrenceQP.f90 Strings.o
./runVerticalRecurrenceQP.x 

ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPAGen.x runVerticalRecurrenceCPUQPAGen.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPBGen.x runVerticalRecurrenceCPUQPBGen.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPCGen.x runVerticalRecurrenceCPUQPCGen.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPDGen.x runVerticalRecurrenceCPUQPDGen.F90 mathfile.o IchorPresicion.o

ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPASegQ.x runVerticalRecurrenceCPUQPASegQ.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPBSegQ.x runVerticalRecurrenceCPUQPBSegQ.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPCSegQ.x runVerticalRecurrenceCPUQPCSegQ.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPDSegQ.x runVerticalRecurrenceCPUQPDSegQ.F90 mathfile.o IchorPresicion.o

ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPASegP.x runVerticalRecurrenceCPUQPASegP.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPBSegP.x runVerticalRecurrenceCPUQPBSegP.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPCSegP.x runVerticalRecurrenceCPUQPCSegP.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPDSegP.x runVerticalRecurrenceCPUQPDSegP.F90 mathfile.o IchorPresicion.o

ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPASeg.x runVerticalRecurrenceCPUQPASeg.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPBSeg.x runVerticalRecurrenceCPUQPBSeg.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPCSeg.x runVerticalRecurrenceCPUQPCSeg.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPDSeg.x runVerticalRecurrenceCPUQPDSeg.F90 mathfile.o IchorPresicion.o

ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPASeg1Prim.x runVerticalRecurrenceCPUQPASeg1Prim.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPBSeg1Prim.x runVerticalRecurrenceCPUQPBSeg1Prim.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPCSeg1Prim.x runVerticalRecurrenceCPUQPCSeg1Prim.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o runVerticalRecurrenceCPUQPDSeg1Prim.x runVerticalRecurrenceCPUQPDSeg1Prim.F90 mathfile.o IchorPresicion.o


ifort -g -check all -debug all -openmp -o BUILDRJ000CPUQPSeg1Prim.x BUILDRJ000CPUQPSeg1Prim.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o BUILDRJ000CPUQPGen.x BUILDRJ000CPUQPGen.F90 mathfile.o IchorPresicion.o


#echo 'doing vec reort'
#ifort -O3 -xHost -vec-report4 -o runVerticalRecurrenceQPAoutput.x runVerticalRecurrenceQPAoutput.F90 mathfile.o IchorPresicion.o
cp runVerticalRecurrenceCPUQPAGen.F90 AGC_CPU_VerticalRecurrenceQPAGen.F90
cp runVerticalRecurrenceCPUQPBGen.F90 AGC_CPU_VerticalRecurrenceQPBGen.F90
cp runVerticalRecurrenceCPUQPCGen.F90 AGC_CPU_VerticalRecurrenceQPCGen.F90
cp runVerticalRecurrenceCPUQPDGen.F90 AGC_CPU_VerticalRecurrenceQPDGen.F90

cp runVerticalRecurrenceCPUQPASegQ.F90 AGC_CPU_VerticalRecurrenceQPASegQ.F90
cp runVerticalRecurrenceCPUQPBSegQ.F90 AGC_CPU_VerticalRecurrenceQPBSegQ.F90
cp runVerticalRecurrenceCPUQPCSegQ.F90 AGC_CPU_VerticalRecurrenceQPCSegQ.F90
cp runVerticalRecurrenceCPUQPDSegQ.F90 AGC_CPU_VerticalRecurrenceQPDSegQ.F90

cp runVerticalRecurrenceCPUQPASegP.F90 AGC_CPU_VerticalRecurrenceQPASegP.F90
cp runVerticalRecurrenceCPUQPBSegP.F90 AGC_CPU_VerticalRecurrenceQPBSegP.F90
cp runVerticalRecurrenceCPUQPCSegP.F90 AGC_CPU_VerticalRecurrenceQPCSegP.F90
cp runVerticalRecurrenceCPUQPDSegP.F90 AGC_CPU_VerticalRecurrenceQPDSegP.F90

cp runVerticalRecurrenceCPUQPASeg.F90 AGC_CPU_VerticalRecurrenceQPASeg.F90
cp runVerticalRecurrenceCPUQPBSeg.F90 AGC_CPU_VerticalRecurrenceQPBSeg.F90
cp runVerticalRecurrenceCPUQPCSeg.F90 AGC_CPU_VerticalRecurrenceQPCSeg.F90
cp runVerticalRecurrenceCPUQPDSeg.F90 AGC_CPU_VerticalRecurrenceQPDSeg.F90

cp runVerticalRecurrenceCPUQPASeg1Prim.F90 AGC_CPU_VerticalRecurrenceQPASeg1Prim.F90
cp runVerticalRecurrenceCPUQPBSeg1Prim.F90 AGC_CPU_VerticalRecurrenceQPBSeg1Prim.F90
cp runVerticalRecurrenceCPUQPCSeg1Prim.F90 AGC_CPU_VerticalRecurrenceQPCSeg1Prim.F90
cp runVerticalRecurrenceCPUQPDSeg1Prim.F90 AGC_CPU_VerticalRecurrenceQPDSeg1Prim.F90

cp BUILDRJ000CPUQPSeg1Prim.F90 AGC_CPU_BuildRJ000Seg1Prim.F90
cp BUILDRJ000CPUQPGen.F90 AGC_CPU_BuildRJ000Gen.F90

#GPU
cp runVerticalRecurrenceGPUQPAGen.F90 AGC_GPU_VerticalRecurrenceQPAGen.F90
cp runVerticalRecurrenceGPUQPBGen.F90 AGC_GPU_VerticalRecurrenceQPBGen.F90
cp runVerticalRecurrenceGPUQPCGen.F90 AGC_GPU_VerticalRecurrenceQPCGen.F90
cp runVerticalRecurrenceGPUQPDGen.F90 AGC_GPU_VerticalRecurrenceQPDGen.F90

cp runVerticalRecurrenceGPUQPASegQ.F90 AGC_GPU_VerticalRecurrenceQPASegQ.F90
cp runVerticalRecurrenceGPUQPBSegQ.F90 AGC_GPU_VerticalRecurrenceQPBSegQ.F90
cp runVerticalRecurrenceGPUQPCSegQ.F90 AGC_GPU_VerticalRecurrenceQPCSegQ.F90
cp runVerticalRecurrenceGPUQPDSegQ.F90 AGC_GPU_VerticalRecurrenceQPDSegQ.F90

cp runVerticalRecurrenceGPUQPASegP.F90 AGC_GPU_VerticalRecurrenceQPASegP.F90
cp runVerticalRecurrenceGPUQPBSegP.F90 AGC_GPU_VerticalRecurrenceQPBSegP.F90
cp runVerticalRecurrenceGPUQPCSegP.F90 AGC_GPU_VerticalRecurrenceQPCSegP.F90
cp runVerticalRecurrenceGPUQPDSegP.F90 AGC_GPU_VerticalRecurrenceQPDSegP.F90

cp runVerticalRecurrenceGPUQPASeg.F90 AGC_GPU_VerticalRecurrenceQPASeg.F90
cp runVerticalRecurrenceGPUQPBSeg.F90 AGC_GPU_VerticalRecurrenceQPBSeg.F90
cp runVerticalRecurrenceGPUQPCSeg.F90 AGC_GPU_VerticalRecurrenceQPCSeg.F90
cp runVerticalRecurrenceGPUQPDSeg.F90 AGC_GPU_VerticalRecurrenceQPDSeg.F90

cp runVerticalRecurrenceGPUQPASeg1Prim.F90 AGC_GPU_VerticalRecurrenceQPASeg1Prim.F90
cp runVerticalRecurrenceGPUQPBSeg1Prim.F90 AGC_GPU_VerticalRecurrenceQPBSeg1Prim.F90
cp runVerticalRecurrenceGPUQPCSeg1Prim.F90 AGC_GPU_VerticalRecurrenceQPCSeg1Prim.F90
cp runVerticalRecurrenceGPUQPDSeg1Prim.F90 AGC_GPU_VerticalRecurrenceQPDSeg1Prim.F90

cp BUILDRJ000GPUQPSeg1Prim.F90 AGC_GPU_BuildRJ000Seg1Prim.F90
cp BUILDRJ000GPUQPGen.F90 AGC_GPU_BuildRJ000Gen.F90

