      subroutine gencoul(l1,l2,l3,l4,makemean,
     *writeTWO,bonn,breit,sameorb)
      implicit double precision(a-h,o-z)
cbs   SUBROUTINE to generate all required radial 
cbs   integrals for the four angular momenta l1-l4
      include 'param.h' 
      logical*4 makemean,writeTWO,bonn,breit,
     *sameorb 
      dimension                      
     *quot1(MxprimL*MxprimL*MxprimL*MxprimL),
     *quotpow1(mxfrac),
     *quot2(MxprimL*MxprimL*MxprimL*MxprimL),
     *quotpow2(mxfrac) 
      max1=1  !starting values for limits of precalculated 
c             ! powers of function Cfunct(X)
      max2=1
cbs   first of all, this routine determines, for which L 
cbs   values the radial integrals have to be solved
cbs   initialize the number of blocks for the different 
cbs   l-combinations  
c     write(6,*) 'start gencoul l1,l2,l3,l4 ',l1,l2,l3,l4
      if (writeTWO) 
     *write(6,*) 'start gencoul l1,l2,l3,l4 ',l1,l2,l3,l4
cbs   no (ss|ss) contributions 
      if (l1.eq.0.and.l2.eq.0.and.l3.eq.0.and.l4.eq.0) return  ! no integrals for <ss|ss> 
      if (makemean) then 
      	if (l1.eq.1.and.l2.eq.0.and.l3.eq.1.and.l4.eq.0) then ! this should be changed for a version 
c      	                                                       of the code to calculate just some 
c                                                             required blocks with given l-values 
      		nblock=1  ! sp sp are the first, so the first block
      		Lstarter(1,1)=1
      	else 
      		nblock=nblock+1 
cbs   nanz should still have the former value = size of last block 
cbs   but there were strange things on silly 
      nanz=ncontrac(Lvalues(1,nblock-1))*
     *ncontrac(Lvalues(2,nblock-1))*
     *ncontrac(Lvalues(3,nblock-1))*
     *ncontrac(Lvalues(4,nblock-1)) 
c     write(6,'(A,2I6)') 'nanzold,nanz ',nanzold,nanz
     	 	Lstarter(1,nblock)=Lstarter(4,nblock-1)+
     *		nanz*Lblocks(4,nblock-1)
      	endif  
      else 
      	if (l1.eq.0.and.l2.eq.1.and.l3.eq.0.and.l4.eq.1) then ! this should be changed for a version 
c                                                             of the code to calculate just some 
c      	                                                       required blocks with given l-values 
      		nblock=1  ! sp sp are the first, so the first block
      		Lstarter(1,1)=1
      	else 
      		nblock=nblock+1 
cbs   nanz should still have the former value = size of last block 
cbs   but there were strange things on silly 
      nanz=ncontrac(Lvalues(1,nblock-1))*
     *ncontrac(Lvalues(2,nblock-1))*
     *ncontrac(Lvalues(3,nblock-1))*
     *ncontrac(Lvalues(4,nblock-1)) 
c     write(6,'(A,2I6)') 'nanzold,nanz ',nanzold,nanz
      		Lstarter(1,nblock)=Lstarter(4,nblock-1)+
     *		nanz*Lblocks(4,nblock-1)
      	endif  
      endif 
c     write(6,*) 'nblock = ',nblock,' Lstarter ',Lstarter(1,nblock)
cbs   keep track of L-values for later purposes   
      Lvalues(1,nblock)=l1
      Lvalues(2,nblock)=l2
      Lvalues(3,nblock)=l3
      Lvalues(4,nblock)=l4
cbs   now nanz is given the new value
      nanz=ncontrac(l1)*ncontrac(l2)*ncontrac(l3)*ncontrac(l4)
      nprimprod=nprimit(l1)*nprimit(l2)*nprimit(l3)*nprimit(l4)
c
      call initfrac(nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4),
     *quot1,quotpow1,quot2,quotpow2,exponents(1,l1),exponents(1,l2),
     *exponents(1,l3),exponents(1,l4))
cbs   prepare the powers needed for cfunctx
c
c
c     There are seven different CASES of integrals following 
c       (   A  --  C) 
c
c     The structure is the same for all cases, therefore comments can be found only on case A
c
c
c
cbs   ###########################################################################################################
cbs   the (+2) cases          CASE A
cbs   ##########################################################################################################
c     write(6,*) ' the  2  0 cases'
      incl1=1  !  Those increments define the case 
      incl3=1
cbs   determine the possible L-values for the integrals by checking for triangular equation 
c
      call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
c
cbs   returns first and last L-values (Lanf,Lend), for which 
cbs   radial integrals have to be calculated 
      if(Lend-Lanf.ge.0) then 
cbs   if there are blocks 
      	Lblocks(1,nblock)=(Lend-Lanf)/2+1 ! L increases in steps of 2, 
cbs                                       due to parity conservation 
      	Lfirst(1,nblock)=Lanf
      	Llast(1,nblock)=Lend 
      else 
      	Lblocks(1,nblock)=0 
      endif 
      if (Lblocks(1,nblock).gt.0) then    ! integrals have to be calculated 
cbs### check, whether integrals fit on array ################
      if  (Lstarter(1,nblock)+nanz*Lblocks(1,nblock).gt.icont4) then 
      write(6,*) 'end at: ',Lstarter(1,nblock)+nanz*Lblocks(1,nblock) 
      stop 'increase icont4 in para.h'
      endif 
cbs### check, whether integrals fit on array ################
      istart=Lstarter(1,nblock) ! gives the address, where to write the contracted integrals 
c     write(6,*) 'istart = ',istart
cbs   ipow1 and ipow2 are the the numbers of powers in the prefactor
cbs   of the function Cfunct 
      ipow1=2+(l2+l4+Llast(1,nblock))/2
      ipow2=2+(l1+l3+incl1+incl3+Llast(1,nblock))/2
      if (ipow1.gt.max1) then 
cbs   those powers have to be generated... 
      call incrpow(max1,ipow1,quot1,quotpow1,
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
      endif 
      if (ipow2.gt.max2) then 
cbs   those powers have to be generated... 
      call incrpow(max2,ipow2,quot2,quotpow2,
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
      endif 
cbs   now loop over possible L-values 
      do Lrun= Lfirst(1,nblock),Llast(1,nblock),2
      		ipow1=2+(l2+l4+Lrun)/2
      		ipow2=2+(l1+l3+incl1+incl3+Lrun)/2
      if (writetwo) write(6,*) 'Lrun ',lrun
c     write(6,*) 'call buildcoul '
c     in buildcoul the radial integrals are calculated 
      		call buildcoul(l1,l2,l3,l4,incl1,incl3,
     *		Lrun,primints,nprimit(l1),nprimit(l2),nprimit(l3),
     *		nprimit(l4), 
     *		exponents(1,l1),exponents(1,l2),
     *		exponents(1,l3),exponents(1,l4),
     *		powexp(1,1,l3,l1,lrun),powexp(1,1,l4,l2,lrun),
     *		quotpow1(1+nprimprod*(ipow1-1)),
     *		quotpow2(1+nprimprod*(ipow2-1)),writeTWO)
cbs   in the contcas_ routines the integrals are contracted, including exponents as prefactors... 
      		if (bonn.or.breit.or.sameorb) then      
      		call contcasASO(l1,l2,l3,l4,istart,writetwo)
                else 
      		call contcasASO(l1,l2,l3,l4,istart,writetwo)
                call contcasAOO(l1,l2,l3,l4,istart,writetwo)
                endif 
      		istart=istart+nanz  ! start-address for the next block of contracted integrals 
      enddo 
      endif   
cbs   ##########################################################################################################
cbs   the (0) cases         CASE  B
cbs   ##########################################################################################################
c     write(6,*) ' the  0  0 cases'
      incl1=0
      incl3=0
      call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
      if(Lend-Lanf.ge.0) then 
      Lblocks(2,nblock)=(Lend-Lanf)/2+1
      Lfirst(2,nblock)=Lanf
      Llast(2,nblock)=Lend 
      Lblocks(3,nblock)=(Lend-Lanf)/2+1
      Lfirst(3,nblock)=Lanf
      Llast(3,nblock)=Lend 
      else 
      Lblocks(2,nblock)=0 
      Lblocks(3,nblock)=0 
      endif 
      Lstarter(2,nblock)=Lstarter(1,nblock)+
     *nanz*Lblocks(1,nblock)
      Lstarter(3,nblock)=Lstarter(2,nblock)+
     *nanz*Lblocks(2,nblock)
cbs   primitive integrals are the same for type 2 and 3  !!!!!
      if (Lblocks(2,nblock).gt.0) then    
cbs### check, whether integrals fit on array ################
      if  (Lstarter(2,nblock)+2*nanz*Lblocks(2,nblock).gt.icont4) then 
      write(6,*) 'end at: ',Lstarter(2,nblock)+2*nanz*Lblocks(2,nblock) 
      stop 'increase icont4 in para.h'
      endif 
cbs### check, whether integrals fit on array ################
      istart=Lstarter(2,nblock)
c     write(6,*) 'istart = ',istart
      istart2=Lstarter(3,nblock)
c     write(6,*) 'istart2 = ',istart2
      ipow1=2+(l2+l4+Llast(2,nblock))/2
      ipow2=2+(l1+l3+incl1+incl3+Llast(2,nblock))/2
      if (ipow1.gt.max1) then 
      call incrpow(max1,ipow1,quot1,quotpow1,
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
      endif 
      if (ipow2.gt.max2) then 
      call incrpow(max2,ipow2,quot2,quotpow2,
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
      endif 
      do Lrun= Lfirst(2,nblock),Llast(2,nblock),2
      ipow1=2+(l2+l4+Lrun)/2
      ipow2=2+(l1+l3+incl1+incl3+Lrun)/2
      if (writetwo) write(6,*) 'Lrun ',lrun
      call buildcoul(l1,l2,l3,l4,incl1,incl3,
     *Lrun,primints,nprimit(l1),nprimit(l2),nprimit(l3),
     *nprimit(l4), 
     *exponents(1,l1),exponents(1,l2),
     *exponents(1,l3),exponents(1,l4),
     *powexp(1,1,l3,l1,lrun),powexp(1,1,l4,l2,lrun),   
     *quotpow1(1+nprimprod*(ipow1-1)),
     *quotpow2(1+nprimprod*(ipow2-1)),writeTWO)
      if (bonn.or.breit.or.sameorb) then      
      call contcasB1SO(l1,l2,l3,l4,istart,writetwo)
      call contcasB2SO(l1,l2,l3,l4,istart2,writetwo)
      else 
      call contcasB1SO(l1,l2,l3,l4,istart,writetwo)
      call contcasB2SO(l1,l2,l3,l4,istart2,writetwo)
      Call contcasB1OO(l1,l2,l3,l4,istart,writetwo)
      Call contcasB2OO(l1,l2,l3,l4,istart2,writetwo)
      endif   
      istart=istart+nanz
      istart2=istart2+nanz
      enddo 
      endif   
cbs   ##########################################################################################################
cbs   the (-2) cases      CASE C
cbs   ##########################################################################################################
c     write(6,*) ' the -2  0 cases'
      if (l1.eq.0.or.l3.eq.0) then 
      Lblocks(4,nblock)=0
      else 
      incl1=-1
      incl3=-1
      call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
      if(Lend-Lanf.ge.0) then 
      Lblocks(4,nblock)=(Lend-Lanf)/2+1
      Lfirst(4,nblock)=Lanf
      Llast(4,nblock)=Lend 
      else 
      Lblocks(4,nblock)=0 
      endif 
      endif
      Lstarter(4,nblock)=Lstarter(3,nblock)+
     *nanz*Lblocks(3,nblock)
      if (Lblocks(4,nblock).gt.0) then    
cbs### check, whether integrals fit on array ################
      if  (Lstarter(4,nblock)+nanz*Lblocks(4,nblock).gt.icont4) then 
      write(6,*) 'end at: ',Lstarter(4,nblock)+nanz*Lblocks(4,nblock) 
      stop 'increase icont4 in para.h'
      endif 
cbs### check, whether integrals fit on array ################
      istart=Lstarter(4,nblock)
c     write(6,*) 'istart = ',istart
      ipow1=2+(l2+l4+Llast(4,nblock))/2
      ipow2=2+(l1+l3+incl1+incl3+Llast(4,nblock))/2
      if (ipow1.gt.max1) then 
      call incrpow(max1,ipow1,quot1,quotpow1,
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
      endif 
      if (ipow2.gt.max2) then 
      call incrpow(max2,ipow2,quot2,quotpow2,
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
      endif 
      do Lrun= Lfirst(4,nblock),Llast(4,nblock),2
      ipow1=2+(l2+l4+Lrun)/2
      ipow2=2+(l1+l3+incl1+incl3+Lrun)/2
      if (writetwo) write(6,*) 'Lrun ',lrun
      call buildcoul(l1,l2,l3,l4,incl1,incl3,
     *Lrun,primints,nprimit(l1),nprimit(l2),nprimit(l3),
     *nprimit(l4), 
     *exponents(1,l1),exponents(1,l2),
     *exponents(1,l3),exponents(1,l4),
     *powexp(1,1,l3,l1,lrun),powexp(1,1,l4,l2,lrun),
     *quotpow1(1+nprimprod*(ipow1-1)),
     *quotpow2(1+nprimprod*(ipow2-1)),writeTWO)
      if (bonn.or.breit.or.sameorb) then      
      call contcasCSO(l1,l2,l3,l4,istart,writetwo)
      else 
      call contcasCSO(l1,l2,l3,l4,istart,writetwo)
      call contcasCOO(l1,l2,l3,l4,istart,writetwo)
      endif 
      istart=istart+nanz
      enddo 
      endif   
c     write(6,*) 'end gencoul'
c     nanzold=nanz
      return 
      end  
