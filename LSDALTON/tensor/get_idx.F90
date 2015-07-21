module get_idx_mod
  use tensor_parameters_and_counters
  
  public :: get_midx, get_cidx
  private

  interface get_midx
    module procedure get_mode_idx8888,&
                    &get_mode_idx8884,&
                    &get_mode_idx8848,&
                    &get_mode_idx8844,&
                    &get_mode_idx8488,&
                    &get_mode_idx8484,&
                    &get_mode_idx8448,&
                    &get_mode_idx8444,&
                    &get_mode_idx4888,&
                    &get_mode_idx4844,&
                    &get_mode_idx4484,&
                    &get_mode_idx4444
  end interface get_midx

  interface get_cidx
    module procedure get_comp_idx888,&
                    &get_comp_idx884,&
                    &get_comp_idx848,&
                    &get_comp_idx844,&
                    &get_comp_idx488,&
                    &get_comp_idx484,&
                    &get_comp_idx448,&
                    &get_comp_idx444
  end interface get_cidx

  contains

  !> \author Patrick Ettenhuber
  !> \date June 2012
  !> \brief get mode index from composite index
  !> \author Patrick Ettenhuber
  !> \date June 2012
  !> \brief get mode index from composite index
  subroutine get_mode_idx8888(a,inds,dims,modes)
    implicit none
    integer(kind=tensor_long_int),intent(in)    :: a
    integer(kind=tensor_long_int),intent(in)    :: modes
    integer(kind=tensor_long_int),intent(inout) :: inds(modes)
    integer(kind=tensor_long_int),intent(in)    :: dims(modes)
    include "midx.inc"
  end subroutine get_mode_idx8888
  subroutine get_mode_idx8884(a,inds,dims,modes)
    implicit none
    integer(kind=tensor_long_int),intent(in)    :: a
    integer(kind=tensor_standard_int),intent(in)    :: modes
    integer(kind=tensor_long_int),intent(inout) :: inds(modes)
    integer(kind=tensor_long_int),intent(in)    :: dims(modes)
    include "midx.inc"
  end subroutine get_mode_idx8884
  subroutine get_mode_idx8848(a,inds,dims,modes)
    implicit none
    integer(kind=tensor_long_int),intent(in)    :: a
    integer(kind=tensor_long_int),intent(in)    :: modes
    integer(kind=tensor_long_int),intent(inout) :: inds(modes)
    integer(kind=tensor_standard_int),intent(in)    :: dims(modes)
    include "midx.inc"
  end subroutine get_mode_idx8848
  subroutine get_mode_idx8844(a,inds,dims,modes)
    implicit none
    integer(kind=tensor_long_int),intent(in)    :: a
    integer(kind=tensor_standard_int),intent(in)    :: modes
    integer(kind=tensor_long_int),intent(inout) :: inds(modes)
    integer(kind=tensor_standard_int),intent(in)    :: dims(modes)
    include "midx.inc"
  end subroutine get_mode_idx8844
  subroutine get_mode_idx8488(a,inds,dims,modes)
    implicit none
    integer(kind=tensor_long_int),intent(in)    :: a
    integer(kind=tensor_long_int),intent(in)    :: modes
    integer(kind=tensor_standard_int),intent(inout) :: inds(modes)
    integer(kind=tensor_long_int),intent(in)    :: dims(modes)
    include "midx.inc"
  end subroutine get_mode_idx8488
  subroutine get_mode_idx8484(a,inds,dims,modes)
    implicit none
    integer(kind=tensor_long_int),intent(in)    :: a
    integer(kind=tensor_standard_int),intent(in)    :: modes
    integer(kind=tensor_standard_int),intent(inout) :: inds(modes)
    integer(kind=tensor_long_int),intent(in)    :: dims(modes)
    include "midx.inc"
  end subroutine get_mode_idx8484
  subroutine get_mode_idx8448(a,inds,dims,modes)
    implicit none
    integer(kind=tensor_long_int),intent(in)    :: a
    integer(kind=tensor_long_int),intent(in)    :: modes
    integer(kind=tensor_standard_int),intent(inout) :: inds(modes)
    integer(kind=tensor_standard_int),intent(in)    :: dims(modes)
    include "midx.inc"
  end subroutine get_mode_idx8448
  subroutine get_mode_idx8444(a,inds,dims,modes)
    implicit none
    integer(kind=tensor_long_int),intent(in)    :: a
    integer(kind=tensor_standard_int),intent(in)    :: modes
    integer(kind=tensor_standard_int),intent(inout) :: inds(modes)
    integer(kind=tensor_standard_int),intent(in)    :: dims(modes)
    include "midx.inc"
  end subroutine get_mode_idx8444
  subroutine get_mode_idx4888(a,inds,dims,modes)
    implicit none
    integer(kind=tensor_standard_int),intent(in)    :: a
    integer(kind=tensor_long_int),intent(in)    :: modes
    integer(kind=tensor_long_int),intent(inout) :: inds(modes)
    integer(kind=tensor_long_int),intent(in)    :: dims(modes)
    include "midx.inc"
  end subroutine get_mode_idx4888
  subroutine get_mode_idx4844(a,inds,dims,modes)
    implicit none
    integer(kind=tensor_standard_int),intent(in)    :: a
    integer(kind=tensor_standard_int),intent(in)    :: modes
    integer(kind=tensor_long_int),intent(inout) :: inds(modes)
    integer(kind=tensor_standard_int),intent(in)    :: dims(modes)
    include "midx.inc"
  end subroutine get_mode_idx4844
  subroutine get_mode_idx4488(a,inds,dims,modes)
    implicit none
    integer(kind=tensor_standard_int),intent(in)    :: a
    integer(kind=tensor_long_int),intent(in)    :: modes
    integer(kind=tensor_standard_int),intent(inout) :: inds(modes)
    integer(kind=tensor_long_int),intent(in)    :: dims(modes)
    include "midx.inc"
  end subroutine get_mode_idx4488
  subroutine get_mode_idx4484(a,inds,dims,modes)
    implicit none
    integer(kind=tensor_standard_int),intent(in)    :: a
    integer(kind=tensor_standard_int),intent(in)    :: modes
    integer(kind=tensor_long_int),intent(in)    :: dims(modes)
    integer(kind=tensor_standard_int),intent(inout) :: inds(modes)
    include "midx.inc"
  end subroutine get_mode_idx4484
  subroutine get_mode_idx4444(a,inds,dims,modes)
    implicit none
    integer(kind=tensor_standard_int),intent(in)    :: a
    integer(kind=tensor_standard_int),intent(in)    :: modes
    integer(kind=tensor_standard_int),intent(in)    :: dims(modes)
    integer(kind=tensor_standard_int),intent(inout) :: inds(modes)
    include "midx.inc"
  end subroutine get_mode_idx4444

  !> \author Patrick Ettenhuber
  !> \date June 2012
  !> \brief get composite index from mode index
  function get_comp_idx888(inds,dims,modes) result(a)
    implicit none
    integer(kind=tensor_long_int),intent(in) :: modes
    integer(kind=tensor_long_int),intent(in) :: inds(modes)
    integer(kind=tensor_long_int),intent(in) :: dims(modes)
    include "cidx.inc"
  end function get_comp_idx888
  function get_comp_idx884(inds,dims,modes) result(a)
    implicit none
    integer(kind=tensor_standard_int),intent(in) :: modes
    integer(kind=tensor_long_int),intent(in) :: inds(modes)
    integer(kind=tensor_long_int),intent(in) :: dims(modes)
    include "cidx.inc"
  end function get_comp_idx884
  function get_comp_idx848(inds,dims,modes) result(a)
    implicit none
    integer(kind=tensor_long_int),intent(in) :: modes
    integer(kind=tensor_long_int),intent(in) :: inds(modes)
    integer(kind=tensor_standard_int),intent(in) :: dims(modes)
    include "cidx.inc"
  end function get_comp_idx848
  function get_comp_idx844(inds,dims,modes) result(a)
    implicit none
    integer(kind=tensor_standard_int),intent(in) :: modes
    integer(kind=tensor_long_int),intent(in) :: inds(modes)
    integer(kind=tensor_standard_int),intent(in) :: dims(modes)
    include "cidx.inc"
  end function get_comp_idx844
  function get_comp_idx488(inds,dims,modes) result(a)
    implicit none
    integer(kind=tensor_long_int),intent(in) :: modes
    integer(kind=tensor_standard_int),intent(in) :: inds(modes)
    integer(kind=tensor_long_int),intent(in) :: dims(modes)
    include "cidx.inc"
  end function get_comp_idx488
  function get_comp_idx484(inds,dims,modes) result(a)
    implicit none
    integer(kind=tensor_standard_int),intent(in) :: modes
    integer(kind=tensor_standard_int),intent(in) :: inds(modes)
    integer(kind=tensor_long_int),intent(in) :: dims(modes)
    include "cidx.inc"
  end function get_comp_idx484
  function get_comp_idx448(inds,dims,modes) result(a)
    implicit none
    integer(kind=tensor_long_int),intent(in) :: modes
    integer(kind=tensor_standard_int),intent(in) :: inds(modes)
    integer(kind=tensor_standard_int),intent(in) :: dims(modes)
    include "cidx.inc"
  end function get_comp_idx448
  function get_comp_idx444(inds,dims,modes) result(a)
    implicit none
    integer(kind=tensor_standard_int),intent(in) :: modes
    integer(kind=tensor_standard_int),intent(in) :: inds(modes)
    integer(kind=tensor_standard_int),intent(in) :: dims(modes)
    include "cidx.inc"
  end function get_comp_idx444
end module get_idx_mod

