module xcfun_dalton
  use xcfun
  implicit none
  integer :: funobj = -1
contains
  function xcfun_select_by_name(conf_string)
    integer xcfun_select_by_name
    character(*) conf_string
    funobj = xc_new_functional()
    select case(trim(adjustl(conf_string)))
    case ('LDA')
       call xc_set(funobj,XC_SLATERX,1.0D0);
       call xc_set(funobj,XC_VWN5C,1.0D0);
    case ('BLYP')
       call xc_set(funobj,XC_BECKEX,1.0D0);
       call xc_set(funobj,XC_LYPC,1.0D0);
    case ('PBE')
       call xc_set(funobj,XC_PBEX,1.0D0);
       call xc_set(funobj,XC_PBEC,1.0D0);
    case DEFAULT
       print *,'Unknown functional in xcfun_select_by_name()', conf_string
       call xc_free_functional(funobj)
       funobj = -1
       xcfun_select_by_name = -1       
       return
    end select
    xcfun_select_by_name = 0
  end function xcfun_select_by_name
end module
