ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -c -g -check all -debug all Strings.f90
ifort -g -check all -debug all -o runVerticalRecurrenceQP.x runVerticalRecurrenceQP.f90 Strings.o
./runVerticalRecurrenceQP.x 

ifort -g -check all -debug all -o runVerticalRecurrenceQPAoutput2.x runVerticalRecurrenceQPAoutput2.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceQPBoutput2.x runVerticalRecurrenceQPBoutput2.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceQPCoutput2.x runVerticalRecurrenceQPCoutput2.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceQPDoutput2.x runVerticalRecurrenceQPDoutput2.F90 mathfile.o IchorPresicion.o

ifort -g -check all -debug all -o runVerticalRecurrenceQPAoutput2SegQ.x runVerticalRecurrenceQPAoutput2SegQ.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceQPBoutput2SegQ.x runVerticalRecurrenceQPBoutput2SegQ.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceQPCoutput2SegQ.x runVerticalRecurrenceQPCoutput2SegQ.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceQPDoutput2SegQ.x runVerticalRecurrenceQPDoutput2SegQ.F90 mathfile.o IchorPresicion.o

ifort -g -check all -debug all -o runVerticalRecurrenceQPAoutput2SegP.x runVerticalRecurrenceQPAoutput2SegP.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceQPBoutput2SegP.x runVerticalRecurrenceQPBoutput2SegP.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceQPCoutput2SegP.x runVerticalRecurrenceQPCoutput2SegP.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceQPDoutput2SegP.x runVerticalRecurrenceQPDoutput2SegP.F90 mathfile.o IchorPresicion.o

ifort -g -check all -debug all -o runVerticalRecurrenceQPAoutput2Seg.x runVerticalRecurrenceQPAoutput2Seg.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceQPBoutput2Seg.x runVerticalRecurrenceQPBoutput2Seg.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceQPCoutput2Seg.x runVerticalRecurrenceQPCoutput2Seg.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceQPDoutput2Seg.x runVerticalRecurrenceQPDoutput2Seg.F90 mathfile.o IchorPresicion.o

ifort -g -check all -debug all -o runVerticalRecurrenceQPAoutput2Seg1Prim.x runVerticalRecurrenceQPAoutput2Seg1Prim.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceQPBoutput2Seg1Prim.x runVerticalRecurrenceQPBoutput2Seg1Prim.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceQPCoutput2Seg1Prim.x runVerticalRecurrenceQPCoutput2Seg1Prim.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -o runVerticalRecurrenceQPDoutput2Seg1Prim.x runVerticalRecurrenceQPDoutput2Seg1Prim.F90 mathfile.o IchorPresicion.o
#echo 'doing vec reort'
#ifort -O3 -xHost -vec-report4 -o runVerticalRecurrenceQPAoutput.x runVerticalRecurrenceQPAoutput.F90 mathfile.o IchorPresicion.o
cp runVerticalRecurrenceQPAoutput2.F90 AGC_VerticalRecurrenceQPA.F90
cp runVerticalRecurrenceQPBoutput2.F90 AGC_VerticalRecurrenceQPB.F90
cp runVerticalRecurrenceQPCoutput2.F90 AGC_VerticalRecurrenceQPC.F90
cp runVerticalRecurrenceQPDoutput2.F90 AGC_VerticalRecurrenceQPD.F90

cp runVerticalRecurrenceQPAoutput2SegQ.F90 AGC_VerticalRecurrenceQPASegQ.F90
cp runVerticalRecurrenceQPBoutput2SegQ.F90 AGC_VerticalRecurrenceQPBSegQ.F90
cp runVerticalRecurrenceQPCoutput2SegQ.F90 AGC_VerticalRecurrenceQPCSegQ.F90
cp runVerticalRecurrenceQPDoutput2SegQ.F90 AGC_VerticalRecurrenceQPDSegQ.F90

cp runVerticalRecurrenceQPAoutput2SegP.F90 AGC_VerticalRecurrenceQPASegP.F90
cp runVerticalRecurrenceQPBoutput2SegP.F90 AGC_VerticalRecurrenceQPBSegP.F90
cp runVerticalRecurrenceQPCoutput2SegP.F90 AGC_VerticalRecurrenceQPCSegP.F90
cp runVerticalRecurrenceQPDoutput2SegP.F90 AGC_VerticalRecurrenceQPDSegP.F90

cp runVerticalRecurrenceQPAoutput2Seg.F90 AGC_VerticalRecurrenceQPASeg.F90
cp runVerticalRecurrenceQPBoutput2Seg.F90 AGC_VerticalRecurrenceQPBSeg.F90
cp runVerticalRecurrenceQPCoutput2Seg.F90 AGC_VerticalRecurrenceQPCSeg.F90
cp runVerticalRecurrenceQPDoutput2Seg.F90 AGC_VerticalRecurrenceQPDSeg.F90

cp runVerticalRecurrenceQPAoutput2Seg1Prim.F90 AGC_VerticalRecurrenceQPASeg1Prim.F90
cp runVerticalRecurrenceQPBoutput2Seg1Prim.F90 AGC_VerticalRecurrenceQPBSeg1Prim.F90
cp runVerticalRecurrenceQPCoutput2Seg1Prim.F90 AGC_VerticalRecurrenceQPCSeg1Prim.F90
cp runVerticalRecurrenceQPDoutput2Seg1Prim.F90 AGC_VerticalRecurrenceQPDSeg1Prim.F90


