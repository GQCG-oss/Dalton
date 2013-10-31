ifort -c -g -check all -debug all mathfile.f90
ifort -g -check all -debug all -o runOBSdriver.x runOBSdriver.f90 mathfile.o
./runOBSdriver.x
#mv MAIN_OBS_DRIVER.f90 IchorEri_CoulombIntegral_OBS_general.F90
