      Subroutine contandmult(Lhigh,makemean,HERMIT,MOLCAS,AIMP)
      implicit double precision (a-h,o-z) 
      include 'param.h'
      include 'ired.h'
      logical*4 makemean,HERMIT,EX,MOLCAS,AIMP 
      character*8 xa,ya,za
      dimension xa(4),ya(4),za(4)
      common /nucleus/ charge
      double precision normasHERMIT(-Lmax:Lmax,0:Lmax)
      data ((normasHERMIT(ml,l),ml=-l,l),l=0,lmax) /1.0d0,
     &     1.0d0,1.0d0,1.0d0,
     &     1.0d0,1.0d0,3.46410162d0,1.0d0,2.0d0,
     &     4.8989795d0,1.0d0,6.3245553d0,-2.5819889d0,6.3245553d0,
     &     2.0d0,-1.6329932d0,
     &     3.4641016d0,4.89897949d0,9.16515139d0,4.3204938d0,
     &   -3.4156503d0,4.3204938d0,18.330303d0,-1.6329932d0,-6.9282032d0/
      IPNT(I,J)=(J*J-J)/2+I 
c###############################################################################
cbs   subroutine to contract radial one-electron integrals 
cbs   and multiply them with angular factors 
c###############################################################################
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
cbs   clean the arrays for cartesian integrals
C
      length3=(numbalLcarT*numbalLcarT+numbalLcarT)/2
      call dzero(length3,onecartallx)
      call dzero(length3,onecartally)
      call dzero(length3,onecartallz)
      call dzero(length3,onecartallx2)
      call dzero(length3,onecartally2)
      call dzero(length3,onecartallz2)
c
cbs   generate one-electron integrals for all L greater/equal 1
      do L=1,Lhigh   
      	call contone(L,oneoverr3(1,L),onecontr(1,1,-Lmax,1,L),
     *  Lmax,contrarray(iaddtyp3(L)),nprimit(L),ncontrac(L),
     *	MxcontL,dummyone,
     *	onecartx(1,1,1,L),onecartY(1,1,1,L),onecartZ(1,1,1,L),
     *	charge)
      Enddo 
c
ckr
ckr   The previous structure was rather inconvenient. normasHERMIT is now
ckr   a data statement. Kenneth Ruud
ckr      if (HERMIT.and.(.not.MOLCAS)) then
cbs     read factors to normalize like HERMIT, so far implemented up 
cbs     to L=4
cbs     clean up a bit...
ckr      	do Lrun=0,Lhigh 
ckr      	do Mrun=-Lrun,Lrun 
ckr      	normasHERMIT(Mrun,Lrun)=0d0
ckr      	enddo 
ckr      	enddo 
cbs     look for the file with renormalization coefficients 
ckr      	inquire(File='RENORMAL_HERMIT',exist=EX)
ckr      	if (EX) then 
ckr      		open(unit=33,status='UNKNOWN',FIle='RENORMAL_HERMIT') 
ckr                rewind(33)
ckr      		do lrun=0,Lhigh 
ckr      		read(33,*) (normasHERMIT(mrun,Lrun),
ckr     *		mrun=-Lrun,Lrun)
ckr      		enddo
ckr                close(33)
ckr      	else 
ckr      		stop ' RENORMAL_HERMIT not present SORRY' 
ckr      	endif 
ckr      endif 
cbs   ***********************************************************************
cbs   now move all integrals to one big arrays for X,Y,Z
cbs   ***********************************************************************
      print *, 'normas',(normasHERMIT(i,0),i=-Lmax,Lmax)
      print *, 'normas',(normasHERMIT(i,1),i=-Lmax,Lmax)
      do Lrun=1,Lhigh  !loop over L-values (integrals are diagonal in L)
      mrun=0
      do Msec=-Lrun,Lrun    ! cartesian M-values  (Mfirst,Msec) with 
      do Mfirst=-Lrun,Msec  ! Mfirst <= Msec (actually '=' does never appear
c                             as there is no L-component  in Ag
C
       If (HERMIT.and.(.not.MOLCAS)) then 
       factor=normasHERMIT(Mfirst,Lrun)*normasHERMIT(Msec,Lrun)
       else
       factor=1d0
       endif 
c
cbs   determine  if L_X L_Y or L_Z 
        ipowx=ipowxyz(1,mfirst,Lrun)+ipowxyz(1,msec,Lrun)
        ipowy=ipowxyz(2,mfirst,Lrun)+ipowxyz(2,msec,Lrun)
        ipowz=ipowxyz(3,mfirst,Lrun)+ipowxyz(3,msec,Lrun)
c
      	mrun=mrun+1
cbs     now determine the irreducable representations 
      	iredfirst=iredLM(Mfirst,Lrun)
      	iredsec=iredLM(Msec,Lrun)
cbs     check out which IR is the lower one. 
      	if (iredfirst.le.iredsec) then 
        if (HERMIT.and.(.not.MOLCAS)) then 
        incblock=nmbMperIRL(iredseC,Lrun)
        incrfirst=nmbMperIRL(iredfirsT,Lrun)
        elseif (HERMIT.and.MOLCAS) then
        incblock=1
        incrfirst=1
        endif 
cbs     calculate shift to get to the beginning of the block 
      		iredired=shiftIRIR((iredsec*iredsec-iredsec)/2+iredfirst)
     *		+incrlm(Mfirst,Lrun)*itotalperIR(iredsec)+
     *		incrLM(Msec,Lrun)
       if (mod(ipowx,2).eq.0.and.mod(ipowy,2).eq.1.and.
     * mod(ipowz,2).eq.1) then  
      		do icartfirst=1,ncontrac(Lrun) ! loop over functions first index
      		do icartsec=1,ncontrac(Lrun)   ! loop over functions second index
      		onecartallx(iredired+1+(icartsec-1)*iNcblock)=
     *		onecartallx(iredired+1+(icartsec-1)*incblocK)
     *		+onecartx(icartfirst,icartsec,mrun,Lrun)*factor 
      		enddo
cbs             shift pointer by number of functions in IR
      		iredired=iredired+itotalperIR(iredsec)
     *          *incrfirst                              
      		enddo
        endif 
       if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.0.and.
     * mod(ipowz,2).eq.1) then  
      		do icartfirst=1,ncontrac(Lrun) ! loop over functions first index
      		do icartsec=1,ncontrac(Lrun)   ! loop over functions second index
      	        onecartally(iredired+1+(icartsec-1)*iNcblock)=
     *	        onecartally(iredired+1+(icartsec-1)*incblocK)
     *	        +onecarty(icartfirst,icartsec,mrun,Lrun)*factor 
      		enddo
cbs             shift pointer by number of functions in IR
      		iredired=iredired+itotalperIR(iredsec)*
     *          incrfirst                             
      		enddo
        endif 
       if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.1.and.
     * mod(ipowz,2).eq.0) then  
      		do icartfirst=1,ncontrac(Lrun) ! loop over functions first index
      		do icartsec=1,ncontrac(Lrun)   ! loop over functions second index
      	        onecartallz(iredired+1+(icartsec-1)*iNcblock)=
     *	        onecartallz(iredired+1+(icartsec-1)*incblocK)
     *	        +onecartz(icartfirst,icartsec,mrun,Lrun)*factor 
      		enddo
cbs             shift pointer by number of functions in IR
      		iredired=iredired+itotalperIR(iredsec)*
     *          incrfirst                             
      		enddo
        endif 
      	elseif (iredfirst.gt.iredsec) then 
cbs     In this case, indices are exchanged with respect to former 
cbs     symmetry of blocks. Therefore, there will be a minus sign 
        if (HERMIT.and.(.not.MOLCAS)) then 
        incblock=nmbMperIRL(iredfirst,Lrun)
        incrsec=nmbMperIRL(iredsec,Lrun)
        elseif (HERMIT.and.MOLCAS) then
        incblock=1
        incrsec=1                  
        endif 
c
cbs     calculate shift to get to the beginning of the block 
      		iredired=shiftIRIR((iredfirst*iredfirst-iredfirst)/2+
     *		iredsec)+
     *		incrLM(Msec,Lrun)*itotalperIR(iredfirst)+  
     *		incrLM(Mfirst,Lrun)
       if (mod(ipowx,2).eq.0.and.mod(ipowy,2).eq.1.and.
     * mod(ipowz,2).eq.1) then  
      		do icartsec=1,ncontrac(Lrun) !loop over functions second index
      		do icartfirst=1,ncontrac(Lrun) !loop over functions first index
      	        onecartallx(iredired+1+(icartfirsT-1)*
     *          incblock)=
     *	        onecartallx(iredired+1+(icartfirsT-1)*
     *          incblock) 
     *         -onecartx(icartsec,icartfirst,mrun,Lrun)*factor 
      		enddo
cbs             shift pointer by number of functions in IR
      		iredired=iredired+itotalperIR(iredfirst)*
     *          incrsec                      
      		enddo
        endif 
       if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.0.and.
     * mod(ipowz,2).eq.1) then  
      		do icartsec=1,ncontrac(Lrun) !loop over functions second index
      		do icartfirst=1,ncontrac(Lrun) !loop over functions first index
      	        onecartally(iredired+1+(icartfirsT-1)*
     *          incblock)=
     *          onecartally(iredired+1+(icartfirsT-1)*
     *          incblock) 
     *	       -onecarty(icartsec,icartfirst,mrun,Lrun)*factor 
      		enddo
cbs             shift pointer by number of functions in IR
      		iredired=iredired+itotalperIR(iredfirst)*
     *          incrsec                         
      		enddo
        endif 
       if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.1.and.
     * mod(ipowz,2).eq.0) then  
      		do icartsec=1,ncontrac(Lrun) !loop over functions second index
      		do icartfirst=1,ncontrac(Lrun) !loop over functions first index
      	        onecartallz(iredired+1+(icartfirsT-1)*
     *          incblock)=
     *	        onecartallz(iredired+1+(icartfirsT-1)*
     *          incblock) 
     *	       -onecartz(icartsec,icartfirst,mrun,Lrun)*factor   
      		enddo
      		iredired=iredired+itotalperIR(iredfirst)*
     *          incrsec                         
      		enddo
        endif 
      endif 
      enddo
      enddo
      enddo
C
C
      if (HERMIT) then 
cbs   copy integrals on arrays with no symmetry blocking at all 
cbs   which means huge triangular matrices 
      irun=0
      do norb2=1,numballcarT  
      ired2=iredoffunctnew(norb2)
      norbsh2=norb2-shiftIRED(ired2) 
      do norb1=1,norb2    
      ired1=iredoffunctnew(norb1)
      norbsh1=noRb1-shiftIRED(ired1) 
      irun=irun+1 
      iredirEd=shiftIRIR((ired2*ired2-ired2)/2+
     *		ired1)
      if (ired1.ne.ired2) then 
        onecartallx2(irun)=onecartallx(iredired+ 
     *  norbsh2+(norbsH1-1)*itotalperIR(IREd2)) 
        onecartally2(irun)=onecartally(iredired+ 
     *  norBsh2+(norbsH1-1)*itotalperIR(IreD2)) 
        onecartallz2(irun)=onecartallz(iredired+ 
     *  norbsh2+(norbsh1-1)*itotalperIR(IRed2)) 
      else 
        onecartallx2(irun)=onecartallx(iredired+ 
     *  norbsh2*(norbsH2-1)/2+norbsh1)                 
        onecartallY2(irun)=onecartallY(iredired+ 
     *  norbsh2*(norbsH2-1)/2+norbsh1)                 
        onecartallZ2(irun)=onecartallz(iredired+ 
     *  norbsh2*(norbsH2-1)/2+norbsh1)                 
      endif 
      Enddo 
      enddo 
c     write out integrals for test-purposes 
c     write(6,*) 'The X-integrals '
c     ie=0
c     do inumber=1,numballcart  
c     write(6,*) 'Function and IR: ',inumber,iredoffunctnew(inumber)                 
c     ia=ie+1
c     ie=ie+inumber
c     write(6,'(4(E20.14,X))') (onecartallx2(irun),irun=ia,ie)
c     write(6,*) 
c     enddo 
c     write(6,*) 'The Y-integrals '
c     ie=0
c     do inumber=1,numballcarT
c     write(6,*) 'Function and IR: ',inumber,iredoffunctnew(Inumber)                 
c     ia=ie+1
c     ie=ie+inumber
c     write(6,'(4(E20.14,X))') (onecartally2(irun),irun=ia,ie)
c     write(6,*) 
c     enddo 
c     write(6,*) 'The Z-integrals '
c     ie=0
c     do inumber=1,numballcart 
c     Write(6,*) 'Function and IR: ',inumber,iredoffunctnew(Inumber)                 
c     ia=ie+1
c     ie=ie+inumber
c     write(6,'(4(E20.14,X))') (onecartallZ2(irun),irun=ia,ie)
c     write(6,*) 
c     enddo 
      if (.not.AIMP) then 
c     write a hermit-like file   b.s. 4.10.96   
      write(6,*) 'number of orbitals ',numbalLcarT
      write(6,*) 'length of triangular matrix ', length3
       Luprop=19
c             OPEN(LUPROP,STATUS='UNKNOWN',FORM='UNFORMATTED',
c    *        FILE='AOPROPER_MF')
c      rewind LUPROP  
              write(LUPROP)  xa,numbofsym,(nrtofiperIR(I),
     *        i=1,numbofsym),
     *        numballcart,(Loffunction(I),I=1,numballcart), 
     *        (Moffunction(I),I=1,numballcart),
     *        Lhigh,(ncontrac(I),I=0,Lhigh),
     *        ((nmbMperIRL(i,j),i=1,numbofsym),j=0,Lhigh),          
     *        ((iredLM(i,Lrun),i=-Lrun,Lrun),Lrun=0,Lhigh) 
              write(LUPROP) (onecartallx2(irun),irun=1,length3) 
              write(LUPROP)  Ya  
              write(LUPROP) (onecartalLy2(irun),irun=1,length3)
              write(LUPROP)  Za    
              write(LUPROP) (onecartallZ2(irun),irun=1,Length3)
c     close(luprop)
      else
cbs   reorder for AIMP 
      stop ' you should not be here with DALTON'
      write(6,*) 'reorder integrals for AIMP' 
      length3=ikeeporb*(ikeeporb+1)/2
      write(6,*) 'number of orbitals ',ikeeporb        
      write(6,*) 'length of triangular matrix ', length3
      do irun2=1,ikeeporb
      do irun1=1,irun2
      ind2=ikeeplist(irun2)
      ind1=ikeeplist(irun1)
      ipntold=ipnt(ind1,ind2)
      ipntnew=ipnt(irun1,irun2)  
      onecartallx3(ipntnew)=onecartallx2(ipntold) 
      onecartally3(ipntnew)=onecartally2(ipntold) 
      onecartallz3(ipntnew)=onecartallz2(ipntold) 
      enddo  
      enddo  
      write(6,*) 'transfered to new blocks'
       Luprop=19
              OPEN(LUPROP,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *        FILE='AOPROPER_MF')
	      rewind LUPROP  
              write(LUPROP)  xa,numbofsym,(nrtofiperIR(I),
     *        i=1,numbofsym),
     *        ikeeporb,(Loffunction(ikeeplist(i)),i=1,ikeeporb), 
     *        (Moffunction(ikeeplist(i)),I=1,ikeeporb),   
     *        Lhigh,((NContrac(I)-icore(I)),I=0,Lhigh)          
              write(LUPROP) (onecartallx3(irun),irun=1,length3) 
              write(LUPROP)  Ya  
              write(LUPROP) (onecartalLy3(irun),irun=1,length3)
              write(LUPROP)  Za    
              write(LUPROP) (onecartallZ3(irun),irun=1,Length3)
      close(luprop)
      endif 
      Endif 
c
cbs   the following block is actually not needed in the moment, 
cbs   but may-be in the next version connected to the molecule/cas-line 
C
cbs   move now integrals seperately for X,Y,Z and nonvanishing 
cbs   IRED*IRED blocks only 
cbs   
cbs   ONLY BLOCKS with first IRED lower/equal  then second IRED  
cbs   there is no interaction on diagonal blocks  for higher symmetries 
cbs   
cbs         <-------first-------> 
cbs  |     
cbs  |    1
cbs  |   
cbs  s    2    3 
cbs  e
cbs  c    4    5    6  
cbs  o
cbs  n    7    8    9   10 
cbs  d
cbs  |   11   12   13   14   15   
cbs  |
cbs  |
cbs   
cbs   now loop over second index and write out 
cbs   the blocks of integrals with interacting IRs 
cbs   for X,Y and Z. 
cbs   
cbs   
cbs   
cbs   keep track how much was allready moved  
      ishftX=0
      ishftY=0
      ishftZ=0
      if (HERMIT.and.(.not.MOLCAS)) write(6,*) 
     *'Integrals in HERMIT units and normalization' 
      if (HERMIT.and.MOLCAS) write(6,*) 
     *'Integrals in MOLCAS units and normalization' 
      do iredrun=1,numbofsym 
cbs   get  interacting IR
      IRX=IRwithLX(iredrun)    ! IRX=iredrun*IR(LX)
      IRY=IRwithLY(iredrun)    ! IRY=iredrun*IR(LY)
      IRZ=IRwithLZ(iredrun)    ! IRZ=iredrun*IR(LZ)
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
cbs     non-diagonal-block for L_X
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      	if (IRX.lt.iredrun) then ! lower triangle only 
c      write(6,*) 'for X interaction between IRs ',IRX,iredrun
cbs    check where block starts on onecartall 
      		intstartX=shiftIRIR((iredrun*iredrun-iredrun)/2+IRX)     
cbs          determine length of block
      		do irun=1,itotalperIR(iredrun)*itotalperIR(IRX)
      		Xireds(ishftX+irun)=onecartallx(intstartX+irun)
      		enddo
c               write(6,'(4E20.14)') (onecartallx(intstartX+irun),
c    *          irun=1,itotalperIR(iredrun)*itotalperIR(IRX))
cbs    update number of moved integrals 
      		ishftX=ishftX+itotalperIR(iredrun)*itotalperIR(IRX)
c               write(6,*) 'number of integrals ',
c    *          itotalperIR(iredrun)*itotalperIR(IRX)
      	endif 
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
cbs     diagonal-block for L_X
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      	if (IRX.eq.iredrun) then ! lower triangle only 
        ishftXold=ishftX
c      write(6,*) 'for X interaction between IRs ',IRX,iredrun
cbs    check where block starts on onecartall 
      		intstartX=shiftIRIR((iredrun*iredrun-iredrun)/2+IRX)     
cbs          determine length of block
      		do irunlines=1,itotalperIR(iredrun)
                do irun=1,irunlines
      		Xireds(ishftX+irun)=onecartallx(intstartX+
     *          (itotalperIR(iredrun)-1)*irun+irunlines)
      		enddo
                ishftX=ishftX+irunlines
      		enddo
c               write(6,'(4E20.14)') (Xireds(ishftXold+irun),
c    *          irun=1,(itotalperIR(iredrun)*itotalperIR(IRX)+
c    *          itotalperIR(IRX))/2)
      	endif 
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
cbs     non-diagonal-block for L_Y
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      	if (IRY.lt.iredrun) then ! lower triangle only 
c      write(6,*) 'for Y interaction between IRs ',IRY,iredrun
cbs    check where block starts on onecartall 
      		intstartY=shiftIRIR((iredrun*iredrun-iredrun)/2+IRY)
cbs          determine length of block
      		do irun=1,itotalperIR(iredrun)*itotalperIR(IRY)
      		Yireds(ishftY+irun)=onecartally(intstartY+irun)
      		enddo
c               write(6,'(4E20.14)') (onecartally(intstartY+irun),
c    *          irun=1,itotalperIR(iredrun)*itotalperIR(IRY))
cbs    update number of moved integrals 
      		ishftY=ishftY+itotalperIR(iredrun)*itotalperIR(IRY)
c               write(6,*) 'number of integrals ',
c    *          itotalperIR(iredrun)*itotalperIR(IRy)
      	endif 
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
cbs     diagonal-block for L_Y
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      	if (IRY.eq.iredrun) then ! lower triangle only 
        ishftYold=ishftY
c      write(6,*) 'for Y interaction between IRs ',IRY,iredrun
cbs    check where block starts on onecartall 
      		intstartY=shiftIRIR((iredrun*iredrun-iredrun)/2+IRY)     
cbs          determine length of block
      		do irunlines=1,itotalperIR(iredrun)
                do irun=1,irunlines
      		Yireds(ishftX+irun)=onecartally(intstartY+
     *          (itotalperIR(iredrun)-1)*irun+irunlines)
      		enddo
                ishftY=ishftY+irunlines
      		enddo
c               write(6,'(4E20.14)') (Yireds(ishftYold+irun),
c    *          irun=1,(itotalperIR(iredrun)*itotalperIR(IRY)+
c    *          itotalperIR(IRY))/2)
      	endif 
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
cbs     non-diagonal-block for L_Z
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      	if (IRZ.lt.iredrun) then ! lower triangle only 
c               write(6,*) 'for Z interaction between IRs ',IRZ,iredrun
cbs    check where block starts on onecartall 
      		intstartZ=shiftIRIR((iredrun*iredrun-iredrun)/2+IRZ)
c               write(6,*) 'instartz ',intstartz
cbs          determine length of block
      		do irun=1,itotalperIR(iredrun)*itotalperIR(IRZ)
      		Zireds(ishftZ+irun)=onecartallz(intstartZ+irun)
      		enddo
c               write(6,'(4E20.14)') (onecartallz(intstartZ+irun),
c    *          irun=1,itotalperIR(iredrun)*itotalperIR(IRZ))
cbs    update number of moved integrals 
      		ishftZ=ishftZ+itotalperIR(iredrun)*itotalperIR(IRZ)
c               write(6,*) 'number of integrals ',
c    *          itotalperIR(iredrun)*itotalperIR(IRZ)
      	endif 
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
cbs     diagonal-block for L_Z
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      	if (IRZ.eq.iredrun) then ! lower triangle only 
        ishftZold=ishftZ
c      write(6,*) 'for Z interaction between IRs ',IRZ,iredrun
cbs    check where block starts on onecartall 
      		intstartZ=shiftIRIR((iredrun*iredrun-iredrun)/2+IRZ)     
cbs          determine length of block
      		do irunlines=1,itotalperIR(iredrun)
                do irun=1,irunlines
      		Zireds(ishftZ+irun)=onecartallz(intstartZ+
     *          (itotalperIR(iredrun)-1)*irun+irunlines)
      		enddo
                ishftZ=ishftZ+irunlines
      		enddo
c               write(6,'(4E20.14)') (Zireds(ishftZold+irun),
c    *          irun=1,(itotalperIR(iredrun)*itotalperIR(IRZ)+
c    *          itotalperIR(IRZ))/2)
      	endif 
      enddo
cbs######################################################################
cbs    now copy the integrals on the files 
CBs######################################################################
      open(unit=1,status='UNKNOWN',FIle='SOX-integrals') 
      if (makemean) then 
      write(1,'(A30)') 'no-pair X mean-field-integrals'
      else
      write(1,'(A22)') 'no-pair X 1e-integrals'
      endif 
      write(1,'(I5)') ishftX
      write(1,'(4E20.14)') (Xireds(irun),irun=1,ishftX)
      close(1) 
      open(unit=1,status='UNKNOWN',FIle='SOY-integrals') 
      if (makemean) then 
      write(1,'(A30)') 'no-pair Y mean-field-integrals'
      else
      write(1,'(A22)') 'no-pair Y 1E-Integrals'
      endif 
      write(1,'(I5)') ishftY
      write(1,'(4E20.14)') (Yireds(irun),irun=1,ishftY)
      close(1) 
      open(unit=1,status='UNKNOWN',FIle='SOZ-integrals') 
      if (makemean) then 
      write(1,'(A30)') 'no-pair Z mean-field-integrals'
      else
      write(1,'(A22)') 'no-pair Z 1e-integrals'
      endif 
      write(1,'(I5)') ishftZ
      write(1,'(4E20.14)') (Zireds(irun),irun=1,ishftZ)
      close(1) 
cbs   
cbs   that's it!!
cbs   
      return
      end
