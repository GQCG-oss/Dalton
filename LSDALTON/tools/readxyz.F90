program readxyz
use print_moorb_grid_mod
!use typedef
!use matrix_operations
integer :: nocc
type(LsItem) :: ls
character*(80) :: outputfile


 call getarg(1,outputfile)


    write(6,*) 'Reading lsitem...'
    call read_lsitem_from_disk(ls)

    ls%lupri=6
  
    nocc=int(ls%input%molecule%nelectrons/2)
    write(6,*) 'nocc= ', nocc

    call print_xyz(outputfile,ls)

    stop 'readxyz done'

end program readxyz
