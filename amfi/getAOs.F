      subroutine getAOs(lhigh)                  
      implicit double precision(a-h,o-z) 
cbs   get expansions of atomic orbitals in contracted functions 
      include 'param.h'
      common /nucleus/ charge
      character*12    occtext,occread
      character*18  textnorbmf,textnorbmf2    
      integer closedshells(0:LMAX),openshells(0:LMAX) 
      call getocc(int(charge),closedshells,openshells)
      occtext='OCCUPATION: ' 
      do lrun=0,lhigh
      do irun=1,MxcontL
      do jrun=1,MxcontL
      AOcoeffs(jrun,irun,lrun)=0d0
      enddo
      enddo
      enddo
      write(6,*) 'Orbitals for mean-field' 
      do lrun=0,lhigh 
      write(6,*) 'L= ',lrun 
      write(6,*) 'closed shells ', closedshells(lrun)
      do i=1,closedshells(lrun)
      occup(i,lrun)=2.0
      AOcoeffs(i,i,lrun)=1d0
      enddo 
      i=closedshells(lrun)+1
      occup(i,lrun)=1d0*openshells(lrun)/dfloat(lrun+lrun+1)
      AOcoeffs(i,i,lrun)=1d0
      noccorb(lrun)=i
      write(6,*) 'number of orbitals ',noccorb(lrun)
      do iorbital=1,noccorb(lrun)     
      write(6,'(A,F8.4)') occtext,occup(iorbital,lrun)
      write(6,'(8F10.4)') (AOcoeffs(icont,iorbital,lrun),
     *icont=1,ncontrac(lrun))
      write(6,*) ' '
      enddo
      enddo
      return 
      end 
cbs 
      subroutine getocc(icharge,iclosed,iopen) 
      implicit double precision(a-h,o-z)
      include 'para.h' 
      parameter (ichargemax=36) 
      dimension iclocc(0:Lmax,ichargemax)
      dimension iopocc(0:Lmax,ichargemax)
      dimension iclosed(0:LMAX),iopen(0:LMAX)
      data ((iclocc(i,j),i=0,LMAX),j=1,ichargemax) / 
     & 0 , 0, 0, 0, 0,
     & 1 , 0, 0, 0, 0,
     & 1 , 0, 0, 0, 0,
     & 2 , 0, 0, 0, 0,
     & 2 , 0, 0, 0, 0,
     & 2 , 0, 0, 0, 0,
     & 2 , 0, 0, 0, 0,
     & 2 , 0, 0, 0, 0,
     & 2 , 0, 0, 0, 0,
     & 2 , 0, 0, 0, 0,
     & 2 , 1, 0, 0, 0,
     & 3 , 1, 0, 0, 0,
     & 3 , 1, 0, 0, 0,
     & 3 , 1, 0, 0, 0,
     & 3 , 1, 0, 0, 0,
     & 3 , 1, 0, 0, 0,
     & 3 , 1, 0, 0, 0,
     & 3 , 1, 0, 0, 0,
     & 3 , 2, 0, 0, 0,
     & 4 , 2, 0, 0, 0,
     & 4 , 2, 0, 0, 0,
     & 4 , 2, 0, 0, 0,
     & 4 , 2, 0, 0, 0,
     & 4 , 2, 0, 0, 0,
     & 4 , 2, 0, 0, 0,
     & 4 , 2, 0, 0, 0,
     & 4 , 2, 0, 0, 0,
     & 4 , 2, 0, 0, 0,
     & 4 , 2, 0, 0, 0,
     & 4 , 2, 0, 0, 0,
     & 4 , 2, 1, 0, 0,
     & 4 , 2, 1, 0, 0,
     & 4 , 2, 1, 0, 0,
     & 4 , 2, 1, 0, 0,
     & 4 , 2, 1, 0, 0,
     & 4 , 2, 1, 0, 0/
cbs
      data ((iopocc(i,j),i=0,LMAX),j=1,ichargemax) / 
     & 0 , 0, 0, 0, 0,
     & 0 , 0, 0, 0, 0,
     & 1 , 0, 0, 0, 0,
     & 0 , 0, 0, 0, 0,
     & 0 , 1, 0, 0, 0,
     & 0 , 2, 0, 0, 0,
     & 0 , 3, 0, 0, 0,
     & 0 , 4, 0, 0, 0,
     & 0 , 5, 0, 0, 0,
     & 0 , 6, 0, 0, 0,
     & 1 , 0, 0, 0, 0,
     & 0 , 0, 0, 0, 0,
     & 0 , 1, 0, 0, 0,
     & 0 , 2, 0, 0, 0,
     & 0 , 3, 0, 0, 0,
     & 0 , 4, 0, 0, 0,
     & 0 , 5, 0, 0, 0,
     & 0 , 6, 0, 0, 0,
     & 1 , 0, 0, 0, 0,
     & 0 , 0, 0, 0, 0,
     & 1 , 0, 1, 0, 0,
     & 1 , 0, 2, 0, 0,
     & 1 , 0, 3, 0, 0,
     & 1 , 0, 4, 0, 0,
     & 1 , 0, 5, 0, 0,
     & 1 , 0, 6, 0, 0,
     & 1 , 0, 7, 0, 0,
     & 1 , 0, 8, 0, 0,
     & 1 , 0, 9, 0, 0,
     & 1 , 0,10, 0, 0,
     & 1 , 1, 0, 0, 0,
     & 1 , 2, 0, 0, 0,
     & 1 , 3, 0, 0, 0,
     & 1 , 4, 0, 0, 0,
     & 1 , 5, 0, 0, 0,
     & 1 , 6, 0, 0, 0/
cbs
      if (icharge.gt.ichargemax) then 
      stop 'occupations not implemented'    	
      endif 
      do irun=0,lmax
      iclosed(irun)=iclocc(irun,icharge)
       iopen(irun)=iopocc(irun,icharge) 
      enddo 
      return
      end 
