      program addatom2   
      implicit doUble precision (a-h,o-z) 
      parameter(maxorbs=150)
      parameter(maxorbs3=(maxorbs*maxorbs+maxorbs)/2 )
      parameter(maxsymm=8)                                     
      character*8 xa,ya,za,xa2,ya2,za2
      dimension SCR(maxorbs3) , SCR2(maxorbs3) ,SCR3(maxorbs3) 
      dimension xa(4),ya(4),za(4)
      dimension xa2(4),Ya2(4),za2(4)
      INTEGER*4 SAOsperIR1(8),SAOsperIR2(8),SAOsperIR3(8)  
      dimension newlist1(maxorbs),newlist2(maxorbs) 
c###############################################################################
C
C   addatom:
C
C   Couples two files with spin-orbit integrals to a big file 
c
c   reads: AOPROPER_MF1 AOPROPER_MF2 
c
c   writes: AOPROPER_MF_NEW   
C
C
C
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccC
      IPNT(I,J)=(max(i,j)*max(i,j)-max(i,j))/2 +min(i,j) 
      xa(1)='********'
      ya(1)='********'
      za(1)='********'
      xa(2)='        '
      ya(2)='        '
      Za(2)='        '
      xa(3)='ANTISYMM'
      ya(3)='ANTISYMM'
      Za(3)='ANTISYMM'
      xa(4)='X1SPNORB'
      ya(4)='Y1SPNORB'
      ZA(4)='Z1SPNORB'
c
c
      luprop1=19
      luprop2=20
      luprop3=21
      open(LUPROP1,file='AOPROPER_MF1',FORM='UNFORMATTED',STATUS='OLD')
      open(LUPROP2,file='AOPROPER_MF2',FORM='UNFORMATTED',STATUS='OLD')
      open(LUPROP3,file='AOPROPER_MF_NEW',FORM='UNFORMATTED',
     *STATUS='UNKNOWN')
      rewind luprop1   
      rewind luprop2   
      rewind luprop3   
      read(LUPROP1) Xa2,nsymm1,(SAOsperIR1(I),I=1,nsymm1)
      write(6,*) 'INFO FILE 1: ',nsymm1,(SAOsperIR1(I),I=1,nsymm1)
      if (XA(1).ne.Xa2(1).or.XA(2).ne.Xa2(2).or.XA(3).ne.
     *Xa2(3).or.XA(4).ne.Xa2(4)) then 
      write(6,*) 'wrong header for x '
      endif 
      read(LUPROP2) Xa2,nsymm2,(SAOsperIR2(I),I=1,nsymm2)
      write(6,*) 'INFO FILE 2: ',nsymm2,(SAOsperIR2(I),I=1,nsymm2)
      if (XA(1).ne.Xa2(1).or.XA(2).ne.Xa2(2).or.XA(3).ne.
     *Xa2(3).or.XA(4).ne.Xa2(4)) then 
      write(6,*) 'wrong header for x '
      endif 
cbs   make some checks   
      if (nsymm1.ne.nsymm2) stop 'different numbers  of symmetries'  
      if (nsymm1.gt.maxsymm) stop 'too many symmetries' 
c
c    determine total number of functions  
      norbs1=0
      norbs2=0
      do isymrun=1,nsymm1
      norbs1=norbs1+SAOsperIR1(isymrun) 
      norbs2=norbs2+SAOsperIR2(isymrun) 
      SAOsperIR3(isymrun)=SAOsperIR1(isymrun)+SAOsperIR2(isymrun)
      enddo 
c    
c    check if array is big enough 
c 
      norbs3=norbs1+norbs2 
      if (norbs3.gt.maxorbs) stop 'too many functions' 
      norbs1dr=(norbs1*norbs1+norbs1)/2  
      norbs2dr=(norbs2*norbs2+norbs2)/2  
      norbs3dr=(norbs3*norbs3+norbs3)/2  
      write(LUPROP3) Xa2,nsymm1,(SAOsperIR3(I),I=1,nsymm1)
      write(6,*) 'INFO FILE 3: ',nsymm2,(SAOsperIR3(I),I=1,nsymm2)
C
      call readscr(SCR,NORBS1DR,LUPROP1)
      call readscr(SCR2,NORBS2DR,LUPROP2)
C
      indnew=0
      indold1=0
      indold2=0 
      do isymrun=1,nsymm1
      do iorbs1=1,SAOsperIR1(isymrun)
      indnew=indnew+1
      indold1=indold1+1
      newlist1(indold1)=indnew
      enddo 
      do iorbs2=1,SAOsperIR2(isymrun)
      indnew=indnew+1
      indold2=indold2+1
      newlist2(indold2)=indnew
      enddo 
      enddo 
C
      do irun=1,norbs3dr
      SCR3(irun)=0d0
      enddo 
c 
      iorbrun=0
      do iorb1=1,norbs1  
      do iorb2=1,iorb1
      iorbrun=iorbrun+1   
      SCR3(IPNT(newlist1(iorb1),newlist1(iorb2)))=SCR(iorbrun) 
      enddo   
      enddo   
c 
      iorbrun=0
      do iorb1=1,norbs2  
      do iorb2=1,iorb1
      iorbrun=iorbrun+1   
      SCR3(IPNT(newlist2(iorb1),newlist2(iorb2)))=SCR2(iorbrun) 
      enddo   
      enddo   
C
      call writescR(SCR3,NORBS3DR,LUPROP3)
c   
      read(LUPROP1) Ya2
      if (YA(1).ne.Ya2(1).or.YA(2).ne.Ya2(2).or.YA(3).ne.
     *Ya2(3).or.YA(4).ne.Ya2(4)) then 
      write(6,*) 'wrong header for Y '
      endif 
      read(LUPROP2) Ya2
      if (YA(1).ne.Ya2(1).or.YA(2).ne.Ya2(2).or.YA(3).ne.
     *Ya2(3).or.YA(4).ne.Ya2(4)) then 
      write(6,*) 'wrong header for Y '
      endiF
      write(LUPROP3) Ya2
      call readscr(SCR,NORBS1DR,LUPROP1)
      call readscr(SCR2,NORBS2DR,LUPROP2)
C
      do irun=1,norbs3dr
      SCR3(irun)=0d0
      enddo 
c 
      iorbrun=0
      do iorb1=1,norbs1  
      do iorb2=1,iorb1
      iorbrun=iorbrun+1   
      SCR3(IPNT(newlist1(iorb1),newlist1(iorb2)))=SCR(iorbrun) 
      enddo   
      enddo   
c 
      iorbrun=0
      do iorb1=1,norbs2  
      do iorb2=1,iorb1
      iorbrun=iorbrun+1   
      SCR3(IPNT(newlist2(iorb1),newlist2(iorb2)))=SCR2(iorbrun) 
      enddo   
      enddo   
C
      call writescR(SCR3,NORBS3DR,LUPROP3)
c   
      read(LUPROP1) Za2
      if (zA(1).ne.za2(1).or.zA(2).ne.za2(2).or.zA(3).ne.
     *za2(3).or.zA(4).ne.za2(4)) then 
      write(6,*) 'wrong header for z '
      endif 
      read(LUPROP2) Za2
      if (zA(1).ne.za2(1).or.zA(2).ne.za2(2).or.zA(3).ne.
     *za2(3).or.zA(4).ne.za2(4)) then 
      write(6,*) 'wrong header for z '
      endif 
      write(LUPROP3) za2
      call readscr(SCR,NORBS1Dr,LUPROP1)
      call readscr(SCR2,NORBS2Dr,LUPROP2)
C
      do irun=1,norbs3dr
      SCR3(irun)=0d0
      enddo 
c 
      iorbrun=0
      do iorb1=1,norbs1  
      do iorb2=1,iorb1
      iorbrun=iorbrun+1   
      SCR3(IPNT(newlist1(iorb1),newlist1(iorb2)))=SCR(iorbrun) 
      enddo   
      enddo   
c 
      iorbrun=0
      do iorb1=1,norbs2  
      do iorb2=1,iorb1
      iorbrun=iorbrun+1   
      SCR3(IPNT(newlist2(iorb1),newlist2(iorb2)))=SCR2(iorbrun) 
      enddo   
      enddo   
C
      call writescR(SCR3,NORBS3DR,LUPROP3)
c   
      close(LUPROP1) 
      close(LUPROP2) 
      close(LUPROP3) 
      call system('rm -f fort.19 fort.20 fort.21') 
      stop 
      end 
      subroutine readscr(SCR,norbs3,luprop) 
      double precision scr(norbs3) 
      read(LUPROP) SCR
      return 
      end 
      subroutine writescr(SCR,norbs3,luprop) 
      double precision scr(norbs3) 
      write(LUPROP) SCR
      return 
      end 
      
