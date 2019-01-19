integer,parameter :: name_len=16
integer,parameter :: attr_len=16
integer,parameter :: data_len=16

!
integer,parameter :: real_kind=8 

! For layout
!CJ: the e3 at the end of the format is to use 3 digits for the exponent
character(len=80),parameter :: real_fmt='(e28.16e3)'
integer, parameter          :: tab_default=3
