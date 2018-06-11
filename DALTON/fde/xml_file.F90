!  Copyright (C) 2018 Andre Severo Pereira Gomes, Christoph Jacob, Lucas Visscher and collaborators
!
!  This file is part of Embed, a program implementing the Frozen Density Embedding (FDE) framework
! 
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at http://mozilla.org/MPL/2.0/.
!

module xml_file
  implicit none

  private 
  public open_file, read_var

! Settings. 
! sets the maximum length of variable names.
!  integer :: var_width=132
  
  interface open_file
     module procedure open_file_name, open_file_num
  end interface

  interface read_var
     module procedure read_char
! can be extended for integers and reals
  end interface

! Global variables
!  integer,save :: file_index
        
contains

  integer function open_file_name(filename)
! Opens a file and returns it's identifying number (>0) if the operation was succesful.
!     If the opening was insuccesful a negative or 0 is returned. 
!      0  -   File does not exists.
!     -1  -   To many files opened.
!     -99  -   Error while executing the open statement.
    CHARACTER(LEN=*), INTENT(IN) :: filename
    character(len=1)             :: cur_stat
    INTEGER, parameter :: maxfile=99
    LOGICAL :: ex
    INTEGER :: n
    
    cur_stat=''
    INQUIRE(file=filename, exist=ex)
    IF (.not.(ex)) cur_stat='N'
    
    DO n=1,maxfile
       INQUIRE(n, opened=ex)
       IF (.NOT.(ex)) THEN
          if (cur_stat==' ') then; OPEN(n, file=filename, err=563)   
          else; OPEN(n, file=filename, err=563,STATUS='NEW')   
          end if
          EXIT
       END IF
    END DO
    IF (n>maxfile) THEN
       open_file_name=-1 ! To many files opened
    ELSE
       open_file_name=n ! Open file
    END IF
    return

563 open_file_name=-99

  END FUNCTION OPEN_FILE_NAME

  INTEGER FUNCTION OPEN_FILE_NUM(filename, filenum, formatted, lupri)
!     Opens a file and returns it's identifying number (>0) if the operation was succesful.
!     If the opening was insuccesful a negative or 0 is returned. 
!      0   -   File does not exists.
!     -1   -   To many files opened.
!     -2   -   Number already assigned
!     -3   -   File already opened with a different number.
!     -99  -   Error while executing the open statement.
    CHARACTER(LEN=*), INTENT(IN)  :: filename
    INTEGER, INTENT(IN)           :: filenum
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: formatted
    INTEGER, OPTIONAL, INTENT(IN) :: lupri
    INTEGER, PARAMETER            :: maxfile=99
    LOGICAL                       :: ex
    INTEGER                       :: n

    INQUIRE(file=filename, exist=ex)
    IF (.NOT. ex) THEN
       open_file_num=0
       RETURN
    END IF
    INQUIRE(filenum, exist=ex)
    IF (.NOT. ex) THEN
       open_file_num=-2
       RETURN
    END IF
    OPEN(filenum, file=filename, err=463)   
    open_file_num=filenum
    
    RETURN
463 open_file_num=-99
    
  END FUNCTION OPEN_FILE_NUM

      
  integer function read_char(filenum,var,fixedwidth,skip)
!     TMP adjust for an optional linewith.

    integer           :: filenum
    character(len=*)  :: var
    integer, optional :: fixedwidth
    logical, optional :: skip   

    logical :: is_quoted=.false.,lskip=.false.
!    logical :: is_tag,save=.false.
    integer :: index,ios
    character(len=1) :: chr,prev_chr

    read_char=0
    if (present(skip)) lskip = skip
    var=''; index=1

    select case(prev_chr)
    case('<','>','=')
       var=prev_chr
       prev_chr=''
       return
    case default
       prev_chr=''
    end select

    do 
       read(filenum, '(a)',advance='no',iostat=ios) chr
       if (ios==0) then
          select case(chr)
          case('<','>','=')
             if (.not.(is_quoted)) then
!               if (.not.(chr=='')) prev_chr=chr
               prev_chr=chr
!               var(index:index+1)=chr
               return
             end if             
             var(index:index+1)=chr
             index=index+1
          case(' ')
             if (.not.(is_quoted)) return
             var(index:index+1)=chr
             index=index+1
          case('"')
             is_quoted=.not.(is_quoted)
          case default
             var(index:index+1)=chr
             index=index+1
          end select          
       else; read_char=ios; return
       end if
    end do
  end function read_char

END MODULE XML_FILE


