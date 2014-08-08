ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all IchorPresicion.F90
ifort -c -g -check all -debug all Strings.f90
ifort -g -check all -debug all -o runSphContractOBS1.x runSphContractOBS1.f90 mathfile.o Strings.o
./runSphContractOBS1.x
#cp AutoGenCoderunSphContractOBS1_new.f90 ../LSint/AGC_SphContractOBS1.F90
ifort -g -check all -debug all -openmp -o AutoGenCoderunSphContractOBS1_new.x AutoGenCoderunSphContractOBS1_CPU_new.F90 mathfile.o IchorPresicion.o
ifort -g -check all -debug all -openmp -o AutoGenCoderunSphContractOBS1_new.x AutoGenCoderunSphContractOBS1_GPU_new.F90 mathfile.o IchorPresicion.o

cp AutoGenCoderunSphContractOBS1_CPU_new.F90 AGC_CPU_SphContractOBS1.F90
cp AutoGenCoderunSphContractOBS1_GPU_new.F90 AGC_GPU_SphContractOBS1.F90

