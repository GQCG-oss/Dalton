ifort -c -g -check all -debug all mathfile.f90
ifort -c -g -check all -debug all Strings.f90
ifort -g -check all -debug all -o runOBSdriver.x runOBSdriver.f90 mathfile.o Strings.o 
./runOBSdriver.x
cp MAIN_CPU_OBS_DRIVER.f90 IchorEri_CoulombIntegral_CPU_OBS_general.F90
cp MAIN_CPU_OBS_DRIVERGen.f90 IchorEri_CoulombIntegral_CPU_OBS_Gen.F90
cp MAIN_CPU_OBS_DRIVERSegQ.f90 IchorEri_CoulombIntegral_CPU_OBS_SegQ.F90
cp MAIN_CPU_OBS_DRIVERSegP.f90 IchorEri_CoulombIntegral_CPU_OBS_SegP.F90
cp MAIN_CPU_OBS_DRIVERSeg.f90 IchorEri_CoulombIntegral_CPU_OBS_Seg.F90
cp MAIN_CPU_OBS_DRIVERSeg1Prim.f90 IchorEri_CoulombIntegral_CPU_OBS_Seg1Prim.F90

cp MAIN_GPU_OBS_DRIVER.f90 IchorEri_CoulombIntegral_GPU_OBS_general.F90
cp MAIN_GPU_OBS_DRIVERGen.f90 IchorEri_CoulombIntegral_GPU_OBS_Gen.F90
cp MAIN_GPU_OBS_DRIVERSegQ.f90 IchorEri_CoulombIntegral_GPU_OBS_SegQ.F90
cp MAIN_GPU_OBS_DRIVERSegP.f90 IchorEri_CoulombIntegral_GPU_OBS_SegP.F90
cp MAIN_GPU_OBS_DRIVERSeg.f90 IchorEri_CoulombIntegral_GPU_OBS_Seg.F90
cp MAIN_GPU_OBS_DRIVERSeg1Prim.f90 IchorEri_CoulombIntegral_GPU_OBS_Seg1Prim.F90

