ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all Strings.f90
ifort -g -check all -debug all -o runOBSdriver.x runOBSdriver.f90 mathfile.o Strings.o 
./runOBSdriver.x
cp MAIN_OBS_DRIVER.f90 IchorEri_CoulombIntegral_OBS_general.F90
cp MAIN_OBS_DRIVERGen.f90 IchorEri_CoulombIntegral_OBS_Gen.F90
cp MAIN_OBS_DRIVERSegQ.f90 IchorEri_CoulombIntegral_OBS_SegQ.F90
cp MAIN_OBS_DRIVERSegP.f90 IchorEri_CoulombIntegral_OBS_SegP.F90
cp MAIN_OBS_DRIVERSeg.f90 IchorEri_CoulombIntegral_OBS_Seg.F90
cp MAIN_OBS_DRIVERSeg1Prim.f90 IchorEri_CoulombIntegral_OBS_Seg1Prim.F90

