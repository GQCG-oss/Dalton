#!/bin/bash 

pwd_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

list="ccsdpt_kernels"

  for i in $list
     do
	 sed "s/real_pt/real_sp/g" $pwd_dir/$i.F90 > $pwd_dir/tmp1
	 sed "s/subroutine /subroutine sp_/g" $pwd_dir/tmp1 > $pwd_dir/tmp2
	 sed "s/module /module sp_/g" $pwd_dir/tmp2 > $pwd_dir/tmp1

	 sed "s/call ijk_loop/call sp_ijk_loop/g" $pwd_dir/tmp1 > $pwd_dir/tmp2
         sed "s/call abc_loop/call sp_abc_loop/g" $pwd_dir/tmp2 > $pwd_dir/tmp1
         sed "s/call ccsdpt_energy_full/call sp_ccsdpt_energy_full/g" $pwd_dir/tmp1 > $pwd_dir/tmp2
         sed "s/call trip_generator/call sp_trip_generator/g" $pwd_dir/tmp2 > $pwd_dir/tmp1
         sed "s/call ptr_init_/call sp_ptr_init_/g" $pwd_dir/tmp1 > $pwd_dir/tmp2
         sed "s/call ptr_final_/call sp_ptr_final_/g" $pwd_dir/tmp2 > $pwd_dir/tmp1
         sed "s/call ptr_aliasing_/call sp_ptr_aliasing_/g" $pwd_dir/tmp1 > $pwd_dir/tmp2
         sed "s/public :: abc/public :: sp_abc/g" $pwd_dir/tmp2 > $pwd_dir/tmp1
         sed "s/public :: ijk/public :: sp_ijk/g" $pwd_dir/tmp1 > $pwd_dir/tmp2

	 mv $pwd_dir/tmp2 $pwd_dir/$i\_sp.F90
         rm $pwd_dir/tmp1
     done

list="ccsdpt_full"

  for i in $list
     do
         sed "s/real_pt/real_sp/g" $pwd_dir/$i.F90 > $pwd_dir/tmp1
         sed "s/subroutine /subroutine sp_/g" $pwd_dir/tmp1 > $pwd_dir/tmp2
         sed "s/module /module sp_/g" $pwd_dir/tmp2 > $pwd_dir/tmp1

         sed "s/call trip_amplitudes_/call sp_trip_amplitudes_/g" $pwd_dir/tmp1 > $pwd_dir/tmp2
         sed "s/call trip_denom/call sp_trip_denom/g" $pwd_dir/tmp2 > $pwd_dir/tmp1
         sed "s/call ccsdpt_contract_/call sp_ccsdpt_contract_/g" $pwd_dir/tmp1 > $pwd_dir/tmp2
         sed "s/call ls_ddot_acc/call ls_sdot_acc/g" $pwd_dir/tmp2 > $pwd_dir/tmp1
         sed "s/call ls_dgemm_acc/call ls_sgemm_acc/g" $pwd_dir/tmp1 > $pwd_dir/tmp2
         sed "s/call array_reorder_3d/call array_reorder_3d_sp/g" $pwd_dir/tmp2 > $pwd_dir/tmp1
         sed "s/call array_reorder_3d_acc/call array_reorder_3d_acc_sp/g" $pwd_dir/tmp1 > $pwd_dir/tmp2
         sed "s/public :: trip/public :: sp_trip/g" $pwd_dir/tmp2 > $pwd_dir/tmp1
         sed "s/public :: ccsdpt/public :: sp_ccsdpt/g" $pwd_dir/tmp1 > $pwd_dir/tmp2

         mv $pwd_dir/tmp2 $pwd_dir/$i\_sp.F90
         rm $pwd_dir/tmp1
     done

