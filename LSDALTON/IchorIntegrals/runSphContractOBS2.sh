ifort -c -g -check all -debug all mathfile.f90
ifort -g -check all -debug all -o runSphContractOBS2.x runSphContractOBS2.f90 mathfile.o
./runSphContractOBS2.x
cp AutoGenCoderunSphContractOBS2_new.f90 AGC_SphContractOBS2.F90

