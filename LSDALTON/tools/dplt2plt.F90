program dplt2plt
implicit none
character*(80) :: inputfile
character*(80) :: outputfile
character*(80) :: corb
integer        :: iorb

integer        :: UINP,UOUT
integer        :: nZ,nY,nX,nORBITALS,nGRIDPOINTS,nGRS
real(4)        :: fZ1,fZn,fY1,fYn,fX1,fXn
real(4),allocatable :: moorb(:), buffer(:)
integer        :: i,reclen, nCHUNKS, pos, irec, nLeftover, II,TD
integer        :: CEILING

 call getarg(1,inputfile)
 call getarg(2,outputfile)
 call getarg(3,corb)
 read(corb,*) iorb

 !open input file
 UINP=111; UOUT=112
 OPEN (UINP,FILE=trim(inputfile),STATUS='UNKNOWN',FORM='UNFORMATTED',&
        & ACCESS='DIRECT',RECL=512)
 READ(UINP,REC=1) nZ,nY,nX,fZ1,fZn,fY1,fYn,fX1,fXn,nORBITALS,nGRIDPOINTS,nGRS

 CLOSE(UINP,STATUS='KEEP')

 WRITE(*,*) 'Reading file ', trim(inputfile)
 WRITE(*,*) 'Stored: ', nORBITALS, ' orbitals on ', nGRIDPOINTS,' gridpoints'
 WRITE(*,*) 'Grid:', nZ,' x ',nY,'x',nX, ' points'
 WRITE(*,*) 'Z direction:', fZ1,fZn,  '(',nZ,')'
 WRITE(*,*) 'Y direction:', fY1,fYn,  '(',nY,')'
 WRITE(*,*) 'X direction:', fX1,fXn,  '(',nX,')'
 WRITE(*,*) 'Record length: ',nGRS , ' (4 byte words)'


 !reopen input file with correct record size
 reclen=nGRS
 OPEN (UINP,FILE=trim(inputfile),STATUS='OLD',FORM='UNFORMATTED',&
        & ACCESS='DIRECT',RECL=reclen)

 allocate(moorb(nGRIDPOINTS),buffer(nGRS))

 !read orbital iorb
 nCHUNKS=CEILING(REAL(nGRIDPOINTS,4)/REAL(nGRS,4))

 WRITE(*,*) 'Reading orbital no: ', iorb

 pos=1
 DO i=1,nCHUNKS-1
   irec = iorb + (i-1)*nORBITALS + 1
   READ(UINP,REC=irec) buffer
   moorb(pos:pos+nGRS-1)=buffer
   pos=pos+nGRS
 ENDDO

 !read last chunk
 nLeftover = mod(nGRIDPOINTS,nGRS)
 IF (nLeftover.eq.0) nLeftover = nGRS
 irec = iorb + (nCHUNKS-1)*nORBITALS + 1
 READ(UINP,REC=irec) buffer
 moorb(pos:pos+nLeftover-1)=buffer(1:nLeftover)

 deallocate(buffer)

 CLOSE(UINP,STATUS='KEEP')
 
 !Write plt file
 WRITE(*,*) 'Writing plt file ', trim(outputfile)

 II=3; TD=200
 reclen=(nGRIDPOINTS + 11)

 OPEN (UOUT,FILE=trim(outputfile),STATUS='UNKNOWN',FORM='UNFORMATTED',&
        & ACCESS='DIRECT',RECL=reclen)

 WRITE(UOUT,REC=1) II,TD,nZ,nY,nX,fZ1,fZn,fY1,fYn,fX1,fXn,moorb

 CLOSE(UOUT,STATUS='KEEP')

 deallocate(moorb)

 STOP 'DONE'


end program
