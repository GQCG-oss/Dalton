ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -g -check all -debug all -o runSphContractOBS2.x runSphContractOBS2.f90 mathfile.o
./runSphContractOBS2.x
ifort -g -check all -debug all -openmp -o AutoGenCoderunSphContractOBS2_new.x AutoGenCoderunSphContractOBS2_new.F90 mathfile.o IchorPresicion.o
cp AutoGenCoderunSphContractOBS2_new.F90 AGC_SphContractOBS2.F90

