ifort -c -g -check all -debug all mathfile.f90
ifort -g -check all -debug all -o runSphContractOBS1.x runSphContractOBS1.f90 mathfile.o
./runSphContractOBS1.x
#cp AutoGenCoderunSphContractOBS1_new.f90 ../LSint/AGC_SphContractOBS1.F90
