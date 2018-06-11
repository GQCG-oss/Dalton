!  Copyright (C) 2018 Andre Severo Pereira Gomes, Christoph Jacob, Lucas Visscher and collaborators
!
!  This file is part of Embed, a program implementing the Frozen Density Embedding (FDE) framework
! 
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at http://mozilla.org/MPL/2.0/.
!

module fde_io

   use fde_cfg
      
   private

   public fde_file_open
   public fde_file_close
   public read_grid

   interface fde_file_open
      module procedure fde_wrapper_file_open
   end interface fde_file_open

   interface fde_file_close
      module procedure fde_wrapper_file_close
   end interface fde_file_close

   interface read_grid
      module procedure read_grid_onecol
      module procedure read_grid_manycol
   end interface read_grid

   real(kind=8) :: threshold = 1.0d-18

   contains

     subroutine fde_wrapper_file_open(name,unit)
        character(len=60), intent(in) :: name
        logical                       :: use_qccode_fileops 
        integer, intent(in)           :: unit

        call fde_get_qccode_fileops(use_qccode_fileops)
        if (use_qccode_fileops) then
            call fde_qccode_file_open(unit,name)
        else
            call fde_plain_fortran_file_open(name,unit)
        end if
     end subroutine 

     subroutine fde_wrapper_file_close(unit)
        integer, intent(in) :: unit
        logical             :: use_qccode_fileops 

        call fde_get_qccode_fileops(use_qccode_fileops)
        if (use_qccode_fileops) then
            call fde_qccode_file_close(unit)
        else
            call fde_plain_fortran_file_close(unit)
        end if
     end subroutine

     subroutine fde_plain_fortran_file_open(name,unit)
        character(len=60), intent(in) :: name
        logical                       :: file_found 
        integer, intent(in)           :: unit

        inquire(FILE=name,EXIST=file_found)
  
        if (file_found) then
            open(unit,                   &
                 FORM   = 'FORMATTED',   &
                 STATUS = 'UNKNOWN',     &
                 ACCESS = 'SEQUENTIAL',  &
                 FILE   =  name)
         else
            call fde_quit("File not found!")
         endif
     end subroutine fde_plain_fortran_file_open 

     subroutine fde_plain_fortran_file_close(unit)
        integer, intent(in) :: unit
        close(unit)
     end subroutine fde_plain_fortran_file_close 


! description : this subroutine will read some operator (i.e. an embedding
!               potential) from file.
!
!               the routine expects that the grid file has the following
!               format:
!
!               line        one : number of grid point
!               lines to to n-1 : quadruplet of numbers, where the
!                                 first three are the (x,y,z) coords
!                                 and the last the quadrature weight;
!                                 optionally, the value of the operator/property 
!                                 (the 5th column) is also read in
!

  SUBROUTINE READ_GRID_onecol(file,points,vc)
    REAL(kind=8),POINTER :: points(:,:)
    REAL(kind=8),pointer,optional :: vc(:)
    INTEGER                 :: file
    INTEGER                 :: i,npoints

    REWIND(file)
    READ(file,*) npoints
    
    if (associated(points)) nullify(points)
    allocate(points(4,npoints))
    
    if (present(vc)) then
       if (associated(vc)) nullify(vc)
       allocate(vc(npoints))
       DO i=1,npoints
          READ(file,*) points(:,i),vc(i)
       END DO

    else
       DO i=1,npoints
          READ(file,*) points(:,i)
       END DO
    end if
  END SUBROUTINE READ_GRID_onecol

! description : this subroutine will read one or more operators/properties (i.e. an embedding
!               potential, a density) from file.
!
!               the routine expects that the grid file has the following
!               format:
!
!               line        one : number of grid point, number of properties on file
!               lines to to n-1 : quadruplet of numbers, where the
!                                 first three are the (x,y,z) coords
!                                 and the fourth the quadrature weight.
!                                 these are followed by whatever properties we
!                                 have in file
!
  SUBROUTINE READ_GRID_manycol(file,points,properties)
    REAL(kind=8),POINTER :: points(:,:), ptmp_many(:,:)
    REAL(kind=8),pointer :: properties(:,:)
    INTEGER              :: file
    INTEGER              :: i,npoints, nprop
   
    REWIND(file)
    READ(file,*) npoints, nprop

    if (associated(points)) nullify(points)
    allocate(points(4,npoints))
    if (associated(properties)) nullify(properties)
    allocate(properties(nprop,npoints))

    DO i=1,npoints
       READ(file,*) points(:,i),properties(:,i)
    END DO

  END SUBROUTINE READ_GRID_manycol
  

end module fde_io
