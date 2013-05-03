! Add/subtract plt files to generate new plt file with the same number of gridpoints.
! 
! Manual
! ======
! 1. Compile, e.g.: gfortran AddPltFiles.f90 -o AddPltFiles.x
! 2. Assume you want to generate a plt file output.plt with values:
!    fac1*file1.plt + fac2*file2.plt
!    where file1.plt and file2.plt are existing plt files (of the same size, of course).
!    E.g. use fac1=1.0 and fac2=-1.0 to subtract file2 from file1.
! 3. Run program: ./AddPltFiles.x fac1 file1.plt fac2 file2.plt output.plt
! 

program addplt

  implicit none

  integer(4) :: II,TD,nX,nY,nZ,iunit
  integer(4) :: ngridpoints,reclen,i
  real(4) :: x1,xN,y1,yN,z1,zN, fac1,fac2,tmp
  real(4),allocatable :: values1(:), values2(:), valcom(:)
  character*(80) :: file1,file2,output,charfac1,charfac2

  ! Set stuff according to input
  call getarg(1,charfac1)
  call getarg(2,file1)
  call getarg(3,charfac2)
  call getarg(4,file2)
  call getarg(5,output)
  read (charfac1,'(F20.10)') fac1
  read (charfac2,'(F20.10)') fac2
  II=3; TD=200; IUNIT=111


  ! Determine number of gridpoints
  ! ******************************
  reclen=11+1
  reclen=reclen*4 
  OPEN (IUNIT,FILE=trim(file1),STATUS='OLD',FORM='UNFORMATTED', &
       & ACCESS='DIRECT',RECL=reclen)
  READ(IUNIT,REC=1) II,TD,nZ,nY,nX,Z1,Zn,Y1,Yn,X1,Xn,tmp
  CLOSE(IUNIT,STATUS='KEEP')
  nGRIDPOINTS = nX*nY*nZ

  ! Set record length according to number of gridpoints
  reclen=(nGRIDPOINTS + 11)
  reclen=reclen*4 


  ! Print input information
  ! ***********************
  print *
  print *, '============================'
  print *, '      Input information'
  print *, '============================'
  print *, 'Input File1 = ', trim(file1)
  print *, 'Input File2 = ', trim(file2)
  print *, 'Output File = ', trim(output)
  print *, 'Number of gridpoints = ', ngridpoints
  print *, 'Scaling factors for file 1 and file 2', fac1, fac2
  allocate(values1(ngridpoints))
  allocate(values2(ngridpoints))
  allocate(valcom(ngridpoints))
  values1=0.0_4
  values2=0.0_4
  valcom=0.0_4


  ! File1
  ! *****
  OPEN (IUNIT,FILE=trim(file1),STATUS='OLD',FORM='UNFORMATTED', &
       & ACCESS='DIRECT',RECL=reclen)
  READ(IUNIT,REC=1) II,TD,nZ,nY,nX,Z1,Zn,Y1,Yn,X1,Xn,values1
  CLOSE(IUNIT,STATUS='KEEP')

  ! File2
  ! *****
  OPEN (IUNIT,FILE=trim(file2),STATUS='OLD',FORM='UNFORMATTED', &
       & ACCESS='DIRECT',RECL=reclen)
  READ(IUNIT,REC=1) II,TD,nZ,nY,nX,Z1,Zn,Y1,Yn,X1,Xn,values2
  CLOSE(IUNIT,STATUS='KEEP')

  ! Calculate new plt file according to input factors
  ! *************************************************
  do i=1,ngridpoints
     valcom(i)=fac1*values1(i) + fac2*values2(i)
  end do
  print '(1X,a,3g16.4)', 'Min (file1,file2,output)', minval(values1), minval(values2), minval(valcom)
  print '(1X,a,3g16.4)', 'Max (file1,file2,output)', maxval(values1), maxval(values2), maxval(valcom)

  ! Write new file to output
  ! ************************
  II=3; TD=200; IUNIT=111
  OPEN (IUNIT,FILE=trim(output),STATUS='REPLACE',FORM='UNFORMATTED',&
       & ACCESS='DIRECT',RECL=reclen)
  WRITE(IUNIT,REC=1) II,TD,nZ,nY,nX,Z1,Zn,Y1,Yn,X1,Xn,valcom
  CLOSE(IUNIT,STATUS='KEEP')


end program addplt
