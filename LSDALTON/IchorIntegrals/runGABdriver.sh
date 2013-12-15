ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all Strings.f90
ifort -g -check all -debug all -o runGABdriver.x runGABdriver.f90 mathfile.o Strings.o 
./runGABdriver.x
cp GAB_OBS_DRIVER.f90 IchorEri_GabIntegral_OBS_general.F90
cp GAB_OBS_DRIVERGen.f90 IchorEri_GabIntegral_OBS_Gen.F90
cp GAB_OBS_DRIVERSeg.f90 IchorEri_GabIntegral_OBS_Seg.F90

