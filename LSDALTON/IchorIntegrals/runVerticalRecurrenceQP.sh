ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -c -g -check all -debug all Strings.f90
ifort -g -check all -debug all -o runVerticalRecurrenceQP.x runVerticalRecurrenceQP.f90 Strings.o
./runVerticalRecurrenceQP.x 

ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPA.x runVerticalRecurrenceCPUQPA.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPB.x runVerticalRecurrenceCPUQPB.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPC.x runVerticalRecurrenceCPUQPC.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPD.x runVerticalRecurrenceCPUQPD.F90 mathfile.o IchorPresicion.o

ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPASegQ.x runVerticalRecurrenceCPUQPASegQ.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPBSegQ.x runVerticalRecurrenceCPUQPBSegQ.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPCSegQ.x runVerticalRecurrenceCPUQPCSegQ.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPDSegQ.x runVerticalRecurrenceCPUQPDSegQ.F90 mathfile.o IchorPresicion.o

ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPASegP.x runVerticalRecurrenceCPUQPASegP.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPBSegP.x runVerticalRecurrenceCPUQPBSegP.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPCSegP.x runVerticalRecurrenceCPUQPCSegP.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPDSegP.x runVerticalRecurrenceCPUQPDSegP.F90 mathfile.o IchorPresicion.o

ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPASeg.x runVerticalRecurrenceCPUQPASeg.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPBSeg.x runVerticalRecurrenceCPUQPBSeg.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPCSeg.x runVerticalRecurrenceCPUQPCSeg.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPDSeg.x runVerticalRecurrenceCPUQPDSeg.F90 mathfile.o IchorPresicion.o

ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPASeg1Prim.x runVerticalRecurrenceCPUQPASeg1Prim.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPBSeg1Prim.x runVerticalRecurrenceCPUQPBSeg1Prim.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPCSeg1Prim.x runVerticalRecurrenceCPUQPCSeg1Prim.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceCPUQPDSeg1Prim.x runVerticalRecurrenceCPUQPDSeg1Prim.F90 mathfile.o IchorPresicion.o
#echo 'doing vec reort'
#ifort -O3 -xHost -vec-report4 -o runVerticalRecurrenceQPAoutput.x runVerticalRecurrenceQPAoutput.F90 mathfile.o IchorPresicion.o
cp runVerticalRecurrenceCPUQPA.F90 AGC_CPU_VerticalRecurrenceQPA.F90
cp runVerticalRecurrenceCPUQPB.F90 AGC_CPU_VerticalRecurrenceQPB.F90
cp runVerticalRecurrenceCPUQPC.F90 AGC_CPU_VerticalRecurrenceQPC.F90
cp runVerticalRecurrenceCPUQPD.F90 AGC_CPU_VerticalRecurrenceQPD.F90

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



