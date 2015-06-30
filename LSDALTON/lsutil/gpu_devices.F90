MODULE gpu_device_handling

#ifdef VAR_OPENACC
#ifdef VAR_PGF90
use openacc, only: acc_get_device_type,ACC_DEVICE_NONE,&
     & ACC_DEVICE_DEFAULT,ACC_DEVICE_HOST,ACC_DEVICE_NOT_HOST,&
     & acc_get_num_devices,acc_set_device_num, acc_init, acc_shutdown, &
     & acc_handle_kind, acc_device_kind, ACC_DEVICE_NVIDIA
!ACC_DEVICE_NVIDIA only guarenteed defined in PGF90
#else
use openacc, only: acc_get_device_type,ACC_DEVICE_NONE,&
     & ACC_DEVICE_DEFAULT,ACC_DEVICE_HOST,ACC_DEVICE_NOT_HOST,&
     & acc_get_num_devices,acc_set_device_num, acc_init, acc_shutdown, & 
     & acc_handle_kind, acc_device_kind
#endif
#endif
use iso_c_binding, only: C_int
private

public Init_GPU_devices, Shutdown_GPU_devices,Set_GPU_devices

#ifdef VAR_OPENACC
integer(kind=acc_device_kind),save :: lsdalton_acc_device_type
integer(kind=C_int),save :: lsdalton_max_acc_device_num
integer(kind=C_int),save :: lsdalton_acc_device_num
#else
integer,save :: lsdalton_acc_device_type
integer,save :: lsdalton_max_acc_device_num
integer,save :: lsdalton_acc_device_num
#endif
CONTAINS

  subroutine Init_GPU_devices()
    implicit none
#ifdef VAR_OPENACC
    lsdalton_acc_device_type = acc_get_device_type()
    print*,'acc_get_device_type = ',lsdalton_acc_device_type
    print*,'================================'
    !4 device types are always supported
    print*,'ACC_DEVICE_NONE     = ',ACC_DEVICE_NONE
    print*,'ACC_DEVICE_DEFAULT  = ',ACC_DEVICE_DEFAULT
    print*,'ACC_DEVICE_HOST     = ',ACC_DEVICE_HOST
    print*,'ACC_DEVICE_NOT_HOST = ',ACC_DEVICE_NOT_HOST
#ifdef VAR_PGF90
    !PGF90 supports
    print*,'ACC_DEVICE_NVIDIA   = ',ACC_DEVICE_NVIDIA
#endif
    print*,'================================'

#ifdef VAR_PGF90
    IF(lsdalton_acc_device_type.EQ.ACC_DEVICE_NVIDIA)print*,'ACC_DEVICE_NVIDIA have been selected'
#endif
    IF(lsdalton_acc_device_type.EQ.ACC_DEVICE_NONE)THEN
       print*,'ACC_DEVICE_NONE have been selected'
       print*,'please use the command'
       print*,'export ACC_DEVICE=NVIDIA'
       print*,'to chose the NVIDIA graphics card or '
       print*,'export ACC_DEVICE=HOST'
       print*,'to chose to run the GPU kernel on the CPU'
    ENDIF
    lsdalton_max_acc_device_num = acc_get_num_devices(lsdalton_acc_device_type)
    print*,'There are ',lsdalton_max_acc_device_num,'devices'
    print*,'The Program only support 1 device at present'
    lsdalton_acc_device_num = 1_acc_handle_kind
    IF(acc_get_num_devices(lsdalton_acc_device_num).GT.1_acc_handle_kind)THEN
       print*,'You can change the number of devices by'
       print*,'export ACC_DEVICE_NUM=1'
       print*,'Setting the selected device number to ',lsdalton_acc_device_num
       call acc_set_device_num(lsdalton_acc_device_num,lsdalton_acc_device_type)     
    ENDIF
    call acc_init(lsdalton_acc_device_type)
#endif
  end subroutine Init_GPU_devices

  subroutine Shutdown_GPU_devices()
    implicit none
#ifdef VAR_OPENACC
    call acc_shutdown(lsdalton_acc_device_type)
#endif
  end subroutine Shutdown_GPU_devices

  subroutine Set_GPU_devices(ngpus)
    implicit none
    integer,intent(in) :: ngpus
#ifdef VAR_OPENACC
    lsdalton_acc_device_num = ngpus
    IF(lsdalton_acc_device_num.LT.lsdalton_max_acc_device_num)THEN
       call acc_set_device_num(lsdalton_acc_device_num,lsdalton_acc_device_type)     
    ELSE
       call lsquit('Set_GPU_devices error',-1)
    ENDIF
#endif
  end subroutine Set_GPU_devices

END MODULE gpu_device_handling

