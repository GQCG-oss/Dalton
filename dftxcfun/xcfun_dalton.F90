module xcfun_dalton
  use xcfun
  implicit none
  integer :: funobj = -1
contains
  function xcfun_select_by_name(conf_string)
    integer xcfun_select_by_name
    character(*) conf_string
    funobj = xc_new_functional()
    if (trim(conf_string).eq.'LDA') then
       call xc_set(funobj,XC_SLATERX,1.0D0);
       call xc_set(funobj,XC_VWN5C,1.0D0);
    else if (trim(conf_string).eq.'BLYP') then
       call xc_set(funobj,XC_BECKEX,1.0D0);
       call xc_set(funobj,XC_LYPC,1.0D0);
    else if (trim(conf_string).eq.'PBE') then
       call xc_set(funobj,XC_PBEX,1.0D0);
       call xc_set(funobj,XC_PBEC,1.0D0);
    else
       print *,'Unknown functional in xcfun_select_by_name()'
       call xc_free_functional(funobj)
       funobj = -1
       xcfun_select_by_name = -1       
       return
    endif
    xcfun_select_by_name = 0
  end function xcfun_select_by_name
end module
