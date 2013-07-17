!> @file
!> Contains debugging module

!> \brief Debug routines for linsca branch.
!> \author L. Thogersen. Documented by S. Host.
!> \date 2003
module LINSCA_DEBUG
  use precision
  use files, only: lsopen, lsclose
  use matrix_module,only: matrix
  use matrix_operations
  use matrix_operations_aux, only: mat_write_to_disk2, mat_read_from_disk2
  use integralinterfaceMod, only: II_get_overlap, II_get_nucel_mat, II_get_kinetic
  use opttype, only: optItem
  use typedeftype, only: lssetting
!  use typedef
  use lstiming, only: lstimer
contains

  !> \brief Convert files between formatted and unformatted form.
  !> \author S. Host
  !> \date 2006
  !>
  !> This is not a real debug routine. It writes the initial (typically dens.restart) density
  !> to disk in formatted form, so that it can be transferred to another machine, avoiding
  !> little/big-endian problems. Works only with matrix types dense.
  !>
  subroutine debug_convert_density(opt,D)
  implicit none
    !> Contains info about SCF optimization
    type(optItem),intent(inout) :: opt
    !> Density matrix to be converted
    type(Matrix), intent(inout) :: D
    type(Matrix)                :: Dsparse
    real(realk),allocatable     :: Dfull(:,:)
    integer                     :: luconvert, i, j, ndim
    logical                     :: file_exists,OnMaster
    OnMaster=.TRUE.
  ndim = D%nrow
  luconvert = -1
  if (opt%cfg_which_conversion == 1) then !Write density in formatted form in dens.format, 
                                          !and then stop the program
     CALL lsOPEN(luconvert,'dens.format','unknown','FORMATTED')
     call mat_write_to_disk2(luconvert,D)
     write(opt%lupri,*) 'Formatted density successfully written to dens.format, program will stop as requested'
     STOP 'Formatted density successfully written to dens.format' 
  else if (opt%cfg_which_conversion == 2) then !Take formatted density from dens.format, and return it like
                                               !a standard density matrix, and run the program
     INQUIRE(file='dens.format',EXIST=file_exists)
     if (file_exists) then
        CALL lsOPEN(luconvert,'dens.format','old','FORMATTED')
        call mat_read_from_disk2(luconvert,D)
        write(opt%lupri,*) 'Formatted density successfully read from dens.format!'
     else
        WRITE(opt%LUPRI,'(/A)') 'File dens.format must be present when using DEBUG_CONVERT in input'
        CALL lsQUIT(' dens.format not present ',opt%lupri)
     endif
  else
     WRITE(opt%LUPRI,'(/A)') 'Incorrect conversion specification with DEBUG_CONVERT (can be 1-4)'
     CALL lsQUIT('Incorrect use of keyword DEBUG_CONVERT',opt%lupri)
  endif 
  end subroutine debug_convert_density

  !> \brief Test sparse matrix multilication - requested by R. Andersen
  !> \author S. Host
  !> \date May 2010
  !>
  !> Multiply S by it self and print!
  !>
  subroutine sparsetest(setting, lupri)
  implicit none
     type(lssetting)          :: setting
     type(matrix)             :: S, h, T,ra
     integer, intent(in)      :: lupri
     type(matrix)             :: prod1, prod2
     integer                  :: nsize,luerr,j
     real(realk)              :: tstart,tend
     real(realk)              ::mat_trace
     real(realk) :: alpha

     IF(matrix_type .NE. mtype_csr)call lsQUIT('ERROR: we agreed that you should use the keyword .CSR ',lupri)
     luerr = 6
     nsize = setting%molecule(1)%p%nbastREG
     CALL LSTIMER('START ',tstart,tend,lupri)
     call mat_init(S,nsize,nsize)
     CALL II_get_overlap(lupri,luerr,setting,S)
     CALL LSTIMER('*S     ',tstart,tend,lupri)
     
     alpha = 1.0E-10_realk
     print *,"alpha is ", alpha

!     WRITE(lupri,*)'The overlap matrix test in CSR'
!     call mat_print(S,1,S%nrow,1,S%ncol,lupri)

     CALL LSTIMER('START ',tstart,tend,lupri)
     call mat_init(h,nsize,nsize)
     CALL II_get_nucel_mat(lupri,luerr,setting,h)
     CALL LSTIMER('NUCPOT ',tstart,tend,lupri)
!     WRITE(lupri,*)'The nuclear attraction matrix'
!     call mat_print(h,1,h%nrow,1,h%ncol,lupri)

     CALL LSTIMER('START ',tstart,tend,lupri)
     call mat_init(T,nsize,nsize)
     CALL II_get_kinetic(lupri,luerr,setting,T)
     CALL LSTIMER('Kinetic',tstart,tend,lupri)
!    WRITE(lupri,*)'The kinetic energy matrix'
!    call mat_print(T,1,T%nrow,1,T%ncol,lupri)

!     call mat_daxpy(1E0_realk,T,h)
!     CALL mat_free(T)

     nsize = S%nrow
     !WARNING AT THE MOMENT THIS MAT_INIT ROUTINES ONLY SET
     !prod1%nrow = nsize and prod1%ncol = nsize
     !AND NULLIFIES val,col,row arrays
     call mat_init(prod1,nsize,nsize)
     call mat_init(prod2,nsize,nsize)
!     write (lupri,*) 'Overlap matrix:'

     S%val(2) = 1.99
     print *, "Matrix S:" 
     call MAT_PRINT(S, 1, nsize, 1, nsize, lupri)
     
     print *, "Transposing S:"
     call mat_csr_trans(S, prod2)
     call MAT_PRINT(S, 1, nsize, 1, nsize, lupri)
     print *, "Multiplying S*S = prod1:"
     call mat_mul(prod2,prod2,'n','n',1E0_realk,0E0_realk,prod1)
     print *, "prod1:"
     call MAT_PRINT(prod1, 1, nsize, 1, nsize, lupri)
     alpha = 1.02
     print *, "Scaling S:"
     call mat_scal(alpha, S)
     call MAT_PRINT(S, 1, nsize, 1, nsize, lupri)
     mat_trace = mat_tr(S)
     print *,"Trace is ", mat_trace
     call mat_abs_max_elm(prod1, alpha)
     print *,"max abs elm in prod1 is:", alpha
!     call mat_mul(S,S,'n','n',1E0_realk,0E0_realk,prod1)
!     call MAT_PRINT(prod1, 1, nsize, 1, nsize, lupri)
!     call MAT_PRINT(S, 1, nsize, 1, nsize, lupri)

     print *,"Calculating dot product of prod1 and S:"
     alpha = mat_dotproduct(prod1,S)
     print *,"Dot product between prod1 and S is: ", alpha
     call mat_init(ra,nsize,nsize)
     print *,"Copying prod1 to ra:"
     call mat_assign(ra, prod1)
     call mat_assign(prod2, ra)
     call MAT_PRINT(ra, 1, nsize, 1, nsize, lupri)

!     print *, "Adding S+prod1 = ra:"
!     call mat_csr_add(S,prod1,'n','n',1E0_realk,ra)
!     call MAT_PRINT(ra, 1, nsize, 1, nsize, lupri)

     print *,"matrix S:"
     call MAT_PRINT(S, 1, nsize, 1, nsize, lupri)

     print *,"Calculating trace(ra,S) using transpose and dotproduct"
     alpha = mat_trAB(ra,S)
     print *, "trace(ra, S): ", alpha

     print *,"Calculating trace(ra,S) using mmul and trace"
     call mat_mul(ra,S,'n','n',1E0_realk,0E0_realk,prod2)
     mat_trace = mat_tr(prod2)
     print *,"Trace is ", mat_trace

     print *, "Calculating daxpy  (3.0, S, ra)"
     alpha = 3.0
     call mat_daxpy(alpha,S,ra)
     call MAT_PRINT(ra, 1, nsize, 1, nsize, lupri)
     print *, "Calculating precond:"
     call mat_csr_ao_precond(0, alpha, ra, S,prod1, prod2)
     call MAT_PRINT(prod2, 1, nsize,1, nsize, lupri)
!     write (lupri,*) 'Overlap matrix multiplied by itself:'
!     call MAT_PRINT(prod1, 1, nsize, 1, nsize, lupri)   
     
     !Use S*h as test matrix such that it is not symmetric:
!     call mat_mul(S,h,'n','n',1E0_realk,0E0_realk,prod1)
     
!     write (lupri,*) 'Test matrix:'
!     call MAT_PRINT(prod1, 1, nsize, 1, nsize, lupri)   
!     
!     !Multiply test matrix by itself:
!    call mat_mul(prod1,prod1,'n','n',1E0_realk,0E0_realk,prod2)
!     
!     write (lupri,*) 'Test matrix multiplied by itself:'
!     call MAT_PRINT(prod1, 1, nsize, 1, nsize, lupri)   
 
    
!     call mat_free(S)
!     call mat_free(T)
!     call mat_free(h)
!     call mat_free(ra)
!     call mat_free(prod1)
!     call mat_free(prod2)
     call lsquit('Sparsetest exited normally!',lupri)
     
  end subroutine sparsetest

end module LINSCA_DEBUG
