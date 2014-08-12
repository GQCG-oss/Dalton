ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -c -g -check all -debug all Strings.f90
ifort -g -check all -debug all -o runSphContractOBS2.x runSphContractOBS2.f90 mathfile.o Strings.o
./runSphContractOBS2.x
ifort -g -check all -debug all -openmp -o AutoGenCoderunSphContractOBS2_new.x AutoGenCoderunSphContractOBS2_CPU_new.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o AutoGenCoderunSphContractOBS2_new.x AutoGenCoderunSphContractOBS2_GPU_new.F90 mathfile.o IchorPresicion.o
cp AutoGenCoderunSphContractOBS2_CPU_new.F90 AGC_CPU_SphContractOBS2.F90
cp AutoGenCoderunSphContractOBS2_GPU_new.F90 AGC_GPU_SphContractOBS2.F90

