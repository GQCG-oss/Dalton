      subroutine angular(Lhigh,keep,keepcart,makemean,bonn,
     *breit,writeTWO,sameorb)
c
cbs   COMBINES THE RADIAL INTEGRALS WITH THE ANGULAR FACTORS 
c
cbs   if keep=.true. then 
cbs   all the integrals will be kept in memory. 
cbs   Perhaps, there will be the option to make the 
cbs   transformation to the cartesian basis-sets 
cbs   everytime, they are required. 
cbs   Therefore, the integrals are kept in memory and 
cbs   can be further transformed, whenever required.   
cbs   in order not to waste to much memory, the atomic  
cbs   integrals are thrown away after each l,l,l,l-block 
      implicit double precision(a-h,o-z)
      include 'param.h' 
      logical*4 keep,keepcart,icheck,mcheckxy,mcheckz,makemean,bonn,   
     *writeTWO,breiT,sameorb,cleaner
      dimension ncont(intblmx),l2block(0:Lmax,0:Lmax,0:Lmax,0:Lmax)
cbs #####################################################################
cbs   some preparation of factors needed later on..                     #
cbs ######################################################################
      ipnt(i,j)=(max(i,j)*max(i,j)-max(i,j))/2+min(i,j) 
      roottwo=dsqrt(2d0) 
cbs   calculate some prefactors that will be needed quite often      
      call prefac(Lmax,preroots,clebsch) 
cbs     clean array for one electron integrals
        iprod=MxcontL*MxcontL*(Lmax+Lmax+1)*(Lmax+1)*Lmax     
        call dzero(iprod,onecartX) 
        call dzero(iprod,onecartY) 
        call dzero(iprod,onecartZ) 
cbs   generate an array with sign for (even/odd) m-values
c     write(6,*) 'set isignM'
      isignM(0)=1
      do I=2,Lmax,2
      isignM(I)=1
      isignM(-I)=1
      enddo
      do I=1,Lmax,2
      isignM(I)=-1
      isignM(-I)=-1
      enddo
cbs #####################################################################
cbs   prefactors preXZ und preY include the factors 1/root(2)
cbs   for the +/- linear combinations of spherical harmonics 
cbs #####################################################################
      do M4=-Lmax,Lmax
      do M3=-Lmax,Lmax
      do M2=-Lmax,Lmax
      do M1=-Lmax,Lmax
	      preXZ(m1,m2,m3,m4)=0.25d0
      enddo
      enddo
      enddo
      enddo
      do M3=-Lmax,Lmax
      do M2=-Lmax,Lmax
      do M1=-Lmax,Lmax
	      preXZ(m1,m2,m3,0)=preXZ(m1,m2,m3,0)*roottwo  
      enddo
      enddo
      enddo
      do M3=-Lmax,Lmax
      do M2=-Lmax,Lmax
      do M1=-Lmax,Lmax
	      preXZ(m1,m2,0,m3)=preXZ(m1,m2,0,m3)*roottwo  
      enddo
      enddo
      enddo
      do M3=-Lmax,Lmax
      do M2=-Lmax,Lmax
      do M1=-Lmax,Lmax
	      preXZ(m1,0,m2,m3)=preXZ(m1,0,m2,m3)*roottwo  
      enddo
      enddo
      enddo
      do M3=-Lmax,Lmax
      do M2=-Lmax,Lmax
      do M1=-Lmax,Lmax
	      preXZ(0,m1,m2,m3)=preXZ(0,m1,m2,m3)*roottwo  
      enddo
      enddo
      enddo
      do M4=-Lmax,Lmax
      do M3=-Lmax,Lmax
      do M2=-Lmax,Lmax
      do M1=-Lmax,Lmax
	      preY(m1,m2,m3,m4)=preXZ(m1,m2,m3,m4) 
      enddo
      enddo
      enddo
      enddo
cbs #####################################################################
cbs   additional (-) signs from the (-i) factors  in the 
cbs   (-) linear combinations   (see tosigX(Y,Z).f)
cbs #####################################################################
cbs   + - - -   =>   minus 
      do M4=-Lmax,-1    
      do M3=-Lmax,-1    
         do M2=-Lmax,-1   
         do M1= 0,Lmax
            preXZ(m1,m2,m3,m4)=-preXZ(m1,m2,m3,m4)
         enddo
         enddo
cbs   - + - -   =>   minus 
         do M2= 0,Lmax
         do M1=-Lmax,-1   
            preXZ(m1,m2,m3,m4)=-preXZ(m1,m2,m3,m4)
         enddo
         enddo
      enddo
      enddo
      do M2= 0,Lmax   
      do M1= 0,Lmax
cbs   + + + -   =>   minus 
         do M4=-Lmax,-1    
         do M3= 0,Lmax    
            preXZ(m1,m2,m3,m4)=-preXZ(m1,m2,m3,m4)
         enddo
         enddo
cbs   + + - +   =>   minus 
         do M4= 0,Lmax    
         do M3=-Lmax,-1    
            preXZ(m1,m2,m3,m4)=-preXZ(m1,m2,m3,m4)
         enddo
         enddo
      enddo
      enddo
cbs   + +  - -  >   - 
      do M4=-Lmax,-1  
      do M3=-Lmax,-1  
      do M2=0,Lmax
      do M1=0,Lmax
	      preY(m1,m2,m3,m4)=-preY(m1,m2,m3,m4) 
      enddo
      enddo
      enddo
      enddo
cbs   - -  + +  >   - 
      do M4=0,Lmax
      do M3=0,Lmax
      do M2=-Lmax,-1  
      do M1=-Lmax,-1  
	      preY(m1,m2,m3,m4)=-preY(m1,m2,m3,m4) 
      enddo
      enddo
      enddo
      enddo
cbs #####################################################################
cbs   some quick decision for interaction    
cbs #####################################################################
c     write(6,*) 'set icheckxy' 
      do M4=0,Lmax
      do M3=0,Lmax
      do M2=0,Lmax
      do M1=0,Lmax
	      icheck=mcheckxy(m1,m2,m3,m4)
	      icheckxy(m1,m2,m3,m4)=icheck 
	      icheckxy(m1,m2,m3,-m4)=icheck 
	      icheckxy(m1,m2,-m3,m4)=icheck 
	      icheckxy(m1,-m2,m3,m4)=icheck 
	      icheckxy(-m1,m2,m3,m4)=icheck 
	      icheckxy(m1,m2,-m3,-m4)=icheck 
	      icheckxy(m1,-m2,m3,-m4)=icheck 
	      icheckxy(m1,-m2,-m3,m4)=icheck 
	      icheckxy(m1,-m2,-m3,-m4)=icheck 
	      icheckxy(-m1,m2,m3,-m4)=icheck 
	      icheckxy(-m1,m2,-m3,m4)=icheck 
	      icheckxy(-m1,m2,-m3,-m4)=icheck 
	      icheckxy(-m1,-m2,m3,m4)=icheck 
	      icheckxy(-m1,-m2,m3,-m4)=icheck 
	      icheckxy(-m1,-m2,-m3,m4)=icheck 
	      icheckxy(-m1,-m2,-m3,-m4)=icheck 
      enddo 
      enddo 
      enddo 
      enddo 
c     write(6,*) 'set icheckz' 
      do M4=0,Lmax
      do M3=0,Lmax
      do M2=0,Lmax
      do M1=0,Lmax
	      icheck=mcheckz(m1,m2,m3,m4)
	      icheckz(m1,m2,m3,m4)=icheck 
	      icheckz(m1,m2,m3,-m4)=icheck 
	      icheckz(m1,m2,-m3,m4)=icheck 
	      icheckz(m1,m2,-m3,-m4)=icheck 
	      icheckz(m1,-m2,m3,m4)=icheck 
	      icheckz(m1,-m2,m3,-m4)=icheck 
	      icheckz(m1,-m2,-m3,m4)=icheck 
	      icheckz(m1,-m2,-m3,-m4)=icheck 
	      icheckz(-m1,m2,m3,m4)=icheck 
	      icheckz(-m1,m2,m3,-m4)=icheck 
	      icheckz(-m1,m2,-m3,m4)=icheck 
	      icheckz(-m1,m2,-m3,-m4)=icheck 
	      icheckz(-m1,-m2,m3,m4)=icheck 
	      icheckz(-m1,-m2,m3,-m4)=icheck 
	      icheckz(-m1,-m2,-m3,m4)=icheck 
	      icheckz(-m1,-m2,-m3,-m4)=icheck 
      enddo 
      enddo 
      enddo 
      enddo 
cbs #####################################################################
cbs   there are at most 16 possible combinations of signs ( 2**4) 
cbs #####################################################################
c     write(6,*) 'set interxyz' 
      do M4=0,Lmax
      do M3=0,Lmax
      do M2=0,Lmax
      do M1=0,Lmax
      do irun=1,16
      interxyz(irun,m1,m2,m3,m4)=0         
      enddo
      enddo
      enddo
      enddo
      enddo
cbs   the following M values are the ones from the cartesian 
cbs   linear combinations. interxyz gives the sign sequence 
cbs   for interacting spherical functions, starting with 
cbs   type 1 (++++) and ending with type 16 (-++-)   
      do M4=0,Lmax
      do M3=0,Lmax
      do M2=0,Lmax
      do M1=0,Lmax
      if (icheckxy(m1,m2,m3,m4).or.icheckz(m1,m2,m3,m4)) then 
          irun=0     
          if (iabs(m1+m2-m3-m4).le.1) then 
          irun=irun+1
          interxyz(irun,m1,m2,m3,m4)=1          ! + + + + 
	          if (m1.gt.0.and.m2.gt.0.and.
     *            m3.gt.0.and.m4.gt.0) then
	          irun=irun+1
	          interxyz(irun,m1,m2,m3,m4)=2  ! - - - - 
	          endif 
          endif 
          if (iabs(m1+m2-m3+m4).le.1) then 
	          if (m4.gt.0) then
                  irun=irun+1
                  interxyz(irun,m1,m2,m3,m4)=3  ! + + + - 
	          endif 
	          if (m1.gt.0.and.m2.gt.0.and.
     *            m3.gt.0) then
	          irun=irun+1
	          interxyz(irun,m1,m2,m3,m4)=4  ! - - - + 
	          endif 
          endif 
          if (iabs(m1+m2+m3-m4).le.1) then 
	          if (m3.gt.0) then
                  irun=irun+1
                  interxyz(irun,m1,m2,m3,m4)=5  ! + + - + 
	          endif 
	          if (m1.gt.0.and.m2.gt.0.and.
     *            m4.gt.0) then
	          irun=irun+1
	          interxyz(irun,m1,m2,m3,m4)=6  ! - - + - 
	          endif 
          endif 
          if (iabs(m1-m2-m3-m4).le.1) then 
	          if (m2.gt.0) then
                  irun=irun+1
                  interxyz(irun,m1,m2,m3,m4)=7  ! + - + + 
	          endif 
	          if (m1.gt.0.and.m3.gt.0.and.
     *            m4.gt.0) then
	          irun=irun+1
	          interxyz(irun,m1,m2,m3,m4)=8  ! - + - - 
	          endif 
          endif 
          if (iabs(-m1+m2-m3-m4).le.1) then 
	          if (m1.gt.0) then
                  irun=irun+1
                  interxyz(irun,m1,m2,m3,m4)=9  ! - + + + 
	          endif 
	          if (m2.gt.0.and.m3.gt.0.and.
     *            m4.gt.0) then
	          irun=irun+1
	          interxyz(irun,m1,m2,m3,m4)=10 ! + - - - 
	          endif 
          endif 
          if (iabs(m1+m2+m3+m4).le.1) then 
	          if (m3.gt.0.and.m4.gt.0) then
                  irun=irun+1
                  interxyz(irun,m1,m2,m3,m4)=11 ! + + - - 
	          endif 
	          if (m1.gt.0.and.m2.gt.0) then  
	          irun=irun+1
	          interxyz(irun,m1,m2,m3,m4)=12 ! - - + +   
	          endif 
          endif 
          if (iabs(m1-m2-m3+m4).le.1) then 
	          if (m2.gt.0.and.m4.gt.0) then
                  irun=irun+1
                  interxyz(irun,m1,m2,m3,m4)=13 ! + - + - 
	          endif 
	          if (m1.gt.0.and.m3.gt.0) then  
	          irun=irun+1
	          interxyz(irun,m1,m2,m3,m4)=14 ! - + - +   
	          endif 
          endif 
          if (iabs(m1-m2+m3-m4).le.1) then 
	          if (m2.gt.0.and.m3.gt.0) then
                  irun=irun+1
                  interxyz(irun,m1,m2,m3,m4)=15 ! + - - + 
	          endif 
	          if (m1.gt.0.and.m4.gt.0) then  
	          irun=irun+1
	          interxyz(irun,m1,m2,m3,m4)=16 ! - + + -   
	          endif 
          endif 
c     write(6,*) 'm1,m2,m3,m4 ',m1,m2,m3,m4,
c    *' interactions: ',(interxyz(I,m1,m2,m3,m4),I=1,irun) 
      endif 
      enddo
      enddo
      enddo
      enddo
cbs #####################################################################
cbs   isgnprod gives the sign due to powers (-1)**M  this are again 
cbs   angular m-values 
cbs #####################################################################
c     write(6,*) 'set isgnprod '
      do M4=-Lmax,Lmax
      if (M4.gt.0) then 
      inter4=isignM(M4) 
      else 
      inter4=1
      endif 
      do M3=-Lmax,Lmax
      if (M3.gt.0) then 
      inter3=inter4*isignM(M3) 
      else
      inter3=inter4 
      endif 
      do M2=-Lmax,Lmax
      if (M2.gt.0) then 
      inter2=inter3*isignM(M2) 
      else
      inter2=inter3 
      endif 
      do M1=-Lmax,Lmax
      if (M1.gt.0) then 
      isgnprod(m1,m2,m3,m4)=inter2*isignM(M1)
      else 
      isgnprod(m1,m2,m3,m4)=inter2
      endif   
      enddo
      enddo
      enddo
      enddo
cbs #####################################################################
cbs   some preparation of factors needed later on..  finished           # 
cbs #####################################################################
c
c
c
cbs   counter for total number of cartesian integrals                   !  set some counters 
      numbcart=0                                                        !
cbs   istartX,istartY,istartZ are the adresses for the first cartesian  !  and starting addresses. 
cbs   integrals                                                         !
      istartX(1)=1                                                      !
      istartY(1)=1                                                       
      istartZ(1)=1                                                       
cbs   same orbit integrals integrals  on carteXSO carteYSO and carteZSO                            
cbs   other orbit integrals  on carteXOO carteYOO and carteZOO                            
      iangfirst(1)=1 ! first block of angular integrals                   
cbs #####################################################################
cbs   loop over all (l,l,l,l) blocks generated in the radial part       # 
cbs #####################################################################
      do lrun4=0,Lmax
      do lrun3=0,Lmax
      do lrun2=0,Lmax
      do lrun1=0,Lmax
      l2block(lrun1,lrun2,lrun3,lrun4)=0
      enddo
      enddo
      enddo
      enddo
      do iblock=1,nblock
        l1=lvalues(1,iblock)
        l2=lvalues(2,iblock)
        l3=lvalues(3,iblock)
        l4=lvalues(4,iblock)
        l2block(l1,l2,l3,l4)=iblock  ! can be used for getting the
cbs     address of a cartesian integral  with given l and m-values
        if (writeTWO) write(6,'(A,4I3)') 'in angular l1-4 ',l1,l2,l3,l4
cbs   local counter for integral adresses  
        ilocalX=istartX(iblock)
        ilocalY=istartY(iblock)
        ilocalZ=istartZ(iblock)
        mblock=0 ! counter of (m,m,m,m)-blocks for (l1,l2,l3,l4)   
cbs     number of contracted integrals for each block 
        ncont(iblock)=ncontrac(l1)*ncontrac(l2)*
     *  ncontrac(l3)*ncontrac(l4) 
cbs     if keep is set to false, the angular integrals are 
cbs     thrown away after each block of l-values 
cbs     which means integrals start at address 1 
        if (.not.keep) iangfirst(iblock)=1
        locstar=iangfirst(iblock) ! local starting adress counter 
        do m1=-l1,l1
        do m2=-l2,l2
        do m3=-l3,l3
        do m4=-l4,l4
        mcombina(1,m1,m2,m3,m4)=0  ! will hold type of integrals (1,2,3) 
        mcombina(2,m1,m2,m3,m4)=0  ! will hold number of block 
        enddo 
        enddo 
        enddo 
        enddo 
c       write(6,*) 'the first radial integrals',
c    *  (cont4(Lstarter(I,iblock)),I=1,4) 
        do m1=-l1,l1
        do m2=-l2,l2
        do m3=-l3,l3
cbs     m4 is more or less fixed by m1-3 
c####################################################################################
c####################################################################################
c########## the L- -type block to be combined with sigma+ ###########################
c####################################################################################
c####################################################################################
        m4=m1+m2-m3+1
        if (iabs(m4).le.l4) then !the  L- -block to be combined with sigma+
cbs     not all m-combinations are needed for the mean-field 
        if ((.not.makemean).or.
     *  (l1.eq.l3.and.l2.eq.l4.and.iabs(m2).eq.iabs(m4)).or.
     *  (l1.eq.l2.and.l3.eq.l4.and.
     *  (iabs(m1).eq.iabs(m2).or.iabs(m3).eq.iabs(m4)))) then 
c       write(6,*) 'L- for m4= ',m4
        mcombina(1,m1,m2,m3,m4)=1    
        mblock=mblock+1
        if (locstar+ncont(iblock).gt.mxangint) then 
        write(6,*) 'not enough space allocated for angular integrals'
        write(6,*) 'increase mxangint to at least ',
     *  locstar+ncont(iblock)        
        stop 
        endif  
cbs mkangLmin = make_angular_integrals_for_L- type operator 
cbs really generates  the angular prefactors and combines them with 
cbs the radial integrals
        call mkangLmin(Lmax,l1,l2,l3,l4,m1,m2,m3,m4,
     *       angintSO(locstar),
     *       angintOO(locstar),
     *       Lfirst(1,iblock),Llast(1,iblock),Lblocks(1,iblock),
     *       ncontrac(l1),ncontrac(l2),ncontrac(l3),ncontrac(l4),
     *       cont4SO(Lstarter(1,iblock)),
     *       cont4SO(Lstarter(2,iblock)),
     *       cont4SO(Lstarter(3,iblock)),
     *       cont4SO(Lstarter(4,iblock)),
     *       cont4OO(Lstarter(1,iblock)),
     *       cont4OO(Lstarter(2,iblock)),
     *       cont4OO(Lstarter(3,iblock)),
     *       cont4OO(Lstarter(4,iblock)),
     *       preroots,clebsch,scratch4,bonn,breit,writeTWO,
     *       sameorb) 
        locstar=locstar+ncont(iblock) ! increase starting address 
        mcombina(2,m1,m2,m3,m4)=mblock  ! set the block number 
c####################################################################################
c####################################################################################
c########## the L+ -type block to be combined with sigma- ###########################
c####################################################################################
c####################################################################################
c
c   these integrals are obtained by changing the signs of the m-values.
c   As the integrals are the same, the pointer points to the same integrals...
c
c
        mcombina(1,-m1,-m2,-m3,-m4)=3
        mcombina(2,-m1,-m2,-m3,-m4)=mblock  
        endif 
        Endif 
cbs     write(6,*) 'take next m-values' 
        enddo 
        enddo 
        enddo 
c####################################################################################
c####################################################################################
c########## the L0 -type block to be combined with sigma0 ###########################
c####################################################################################
c####################################################################################
        do m1=  0,l1
        do m2=-l2,l2
        do m3=-l3,l3
cbs     m4 is more or less fixed by m1-3 
        m4=m1+m2-m3 ! the L0-block to be combined with sigma0 
cbs     not all m-combinations are needed for the mean-field 
        if ((.not.makemean).or.
     *  (l1.eq.l3.and.l2.eq.l4.and.iabs(m2).eq.iabs(m4)).or.
     *  (l1.eq.l2.and.l3.eq.l4.and.
     *  (iabs(m1).eq.iabs(m2).or.iabs(m3).eq.iabs(m4)))) then 
c       
        if (m1.ne.0.or.m2.ne.0.or.m3.ne.0) then !  all m eqal 0 make no sense....
        if (iabs(m4).le.l4)  then    
c       write(6,*) 'L0 for m4= ',m4
        mcombina(1,m1,m2,m3,m4)=2
        mblock=mblock+1
        if (locstar+ncont(iblock).gt.mxangint) then 
        write(6,*) 'not enough space allocated for angular integrals'
        write(6,*) 'increase mxangint to at least ',
     *  locstar+ncont(iblock)        
        stop 
        endif  
        call mkangL0(Lmax,l1,l2,l3,l4,m1,m2,m3,m4,
     *       angintSO(locstar),
     *       angintOO(locstar),
     *       Lfirst(1,iblock),Llast(1,iblock),Lblocks(1,iblock),
     *       ncontrac(l1),ncontrac(l2),ncontrac(l3),ncontrac(l4),
     *       cont4SO(Lstarter(1,iblock)),
     *       cont4SO(Lstarter(2,iblock)),
     *       cont4SO(Lstarter(3,iblock)),
     *       cont4SO(Lstarter(4,iblock)),
     *       cont4OO(Lstarter(1,iblock)),
     *       cont4OO(Lstarter(2,iblock)),
     *       cont4OO(Lstarter(3,iblock)),
     *       cont4OO(Lstarter(4,iblock)),
     *       preroots,clebsch,scratch4,bonn,breit,writeTWO,
     *       sameorb) 
        locstar=locstar+ncont(iblock)
        mcombina(2,m1,m2,m3,m4)=mblock  
        endif 
        endif 
        endif 
        enddo
        enddo
        enddo
        iangfirst(iblock+1)=iangfirst(iblock)+mblock*ncont(iblock)   
cbs  ##################################################################################
cbs  ##################################################################################
cbs     transformation to l,m dependent integrals is finished 
cbs  ##################################################################################
c
c
c
c
cbs  ##################################################################################
cbs     begin transformation to cartesian integrals 
cbs  ##################################################################################
cbs  ##################################################################################
cbs     check out, which combinations of m-values will 
cbs     contribute to cartesian integrals    
        do m1=-l1,l1       !    
        do m2=-l2,l2       ! these indices now run over the real harmonics      
        do m3=-l3,l3       !
        do m4=-l4,l4       !
        mcombcart(1,m1,m2,m3,m4,iblock)=0     ! will hold the type  x=1 y=2 z=3 
        mcombcart(2,m1,m2,m3,m4,iblock)=0     ! will hold the block number
        enddo 
        enddo 
        enddo 
        enddo 
        mblockx=0
        mblocky=0
        mblockz=0
        do m3=-l3,l3            
        do m4=-l4,l4       
cbs     if the l-values are the same : triangular matrix over m-values is sufficient       
        if (l1.eq.l3) then    
        m1upper=m3
        else
        m1upper=l1
        endif 
        if (makemean) m1upper=l1
cbs     if the l-values are the same : triangular matrix over m-values is sufficient       
        if (l2.eq.l4) then 
        m2upper=m4
        else
        m2upper=l2 
        endif 
        if (makemean) m2upper=l2
        do m1=-l1,m1upper   
        If (l1.eq.l3.and.m1.eq.m3) then ! clean real zeros by symmetry to be exactly zero 
cbs     this a problem of the spin-other-orbit integrals, as they are by formula 
cbs     not antisymmetric in the indices for particle 1. 
        cleaner=.true.
        else
        cleaner=.false.
        endif  
        do m2=-l2,m2upper   
cbs     not all m-combinations are needed for the mean-field 
        if ((.not.makemean).or.
     *  (l1.eq.l3.and.l2.eq.l4.and.m2.eq.m4).or.
     *  (l1.eq.l2.and.l3.eq.l4.and.(m1.eq.m2.or.m3.eq.m4))) then 
C
c        write(6,'(A,4I3)') 'M-werte ',m1,m2,m3,m4  
        indx=ipowxyz(1,m1,l1)+ipowxyz(1,m2,l2)+
     *  ipowxyz(1,m3,l3)+ipowxyz(1,m4,l4)
        indy=ipowxyz(2,m1,l1)+ipowxyz(2,m2,l2)+
     *  ipowxyz(2,m3,l3)+ipowxyz(2,m4,l4)
        indz=ipowxyz(3,m1,l1)+ipowxyz(3,m2,l2)+
     *  ipowxyz(3,m3,l3)+ipowxyz(3,m4,l4)
        indx=mod(indx,2)
        indy=mod(indy,2)
        indz=mod(indz,2)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++      SIGMA X      ++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if (indx.eq.0.and.indy.eq.1.and.indz.eq.1.and.      
     *  icheckxy(m1,m2,m3,m4)) then  ! Y*Z ->  transforms like  L_x (B1)
cbs     integrals for sigma_x 
        mblockx=mblockx+1
        mcombcart(1,m1,m2,m3,m4,iblock)=1
        mcombcart(2,m1,m2,m3,m4,iblock)=mblockx                               
c       write(6,'(A21,4I3,A7,I4)') 'Integral for sigma_x ',
c    *  m1,m2,m3,m4,' BLOCK ',mblockx
cbs     check if integrals are still fitting on array..
        if ((ilocalX+ncont(iblock)).gt.mxicart) then   
        write(6,*) 'mxicart is too small. Sorry ..'
        stop 
        endif 
        call tosigX(m1,m2,m3,m4,angintSO(iangfirst(iblock)),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),carteXSO(ilocalX),preXZ,writeTWO,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner) 
c
        if (.not.bonn.and.(.not.breiT)) 
     *  call tosigX(m1,m2,m3,m4,angintOO(iangfirst(iblock)),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),carteXOO(ilocalX),preXZ,writeTWO,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner)
        if (makemean) then ! generate mean-field-contributions
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
             if (l1.eq.l3.and.l2.eq.l4) then 
             if (m2.eq.m4.and.m1.lt.m3.and.
     *       iabs(m1+m3).eq.1.and.l1.ne.0) then 
c            write(6,'(5(A,I2))')            
c    *       'contribution to X-mean-field for L= ',
c    *       l1,' m-values ',m1,' and ',m3,' from L= ',l2,' M= ',
c    *       m2
             call two2mean13(cartexSO(ilocalX),occup(1,l2),
     *       AOcoeffs(1,1,l2),onecartx(1,1,ipnt(m1+l1+1,m3+l3+1),l1),
     *       ncontrac(l1),ncontrac(l2),noccorb(l2))         
             endif 
             endif 
             if (l1.eq.l2.and.l3.eq.l4) then 
             if (m1.eq.m2.and.l3.ne.0.and.l3.ne.l1) then 
             if (m3.lt.m4.and.iabs(m4+m3).eq.1) then 
cbs   for the "Bonn-approach"   exchange cartexOO by cartexSO 
c            write(6,'(5(A,I2))')            
c    *       'contribution to X-mean-field for L= ',
c    *       l3,' m-values ',m3,' and ',m4,' from L= ',l1,' M= ',
c    *       m1
             if (bonn.or.breiT) then 
             call two2mean34a(cartexSO(ilocalX),cartexSO(ilocalX),
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartx(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             else 
             call two2mean34a(cartexSO(ilocalX),cartexOO(ilocalX),
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartx(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             endif 
             endif 
             if (m3.gt.m4.and.iabs(m4+m3).eq.1) then 
cbs   for the "Bonn-approach"   exchange cartexOO by cartexSO 
c            write(6,'(5(A,I2))')            
c    *       'contribution to X-mean-field for L= ',
c    *       l3,' m-values ',m3,' and ',m4,' from L= ',l1,' M= ',
c    *       m1
             if (bonn.or.breiT) then 
             call two2mean34b(cartexSO(ilocalX),cartexSO(ilocalX),
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartx(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             else 
             call two2mean34b(cartexSO(ilocalX),cartexOO(ilocalX),
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartx(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             endif 
             endif 
             elseif(m3.eq.m4.and.l1.ne.0) then 
             if (m1.lt.m2.and.iabs(m1+m2).eq.1) then 
cbs   for the "Bonn-approach"   exchange cartexOO by cartexSO 
c            write(6,'(5(A,I2))')            
c    *       'contribution to X-mean-field for L= ',
c    *       l1,' m-values ',m1,' and ',m2,' from L= ',l3,' M= ',
c    *       m3
             if (bonn.or.breiT) then 
             call two2mean12a(cartexSO(ilocalX),cartexSO(ilocalX),
     *       occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartx(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             else 
             call two2mean12a(cartexSO(ilocalX),cartexOO(ilocalX),
     *       occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartx(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             endif 
             endif 
             if (m1.gt.m2.and.iabs(m1+m2).eq.1) then 
cbs   for the "Bonn-approach"   exchange cartexOO by cartexSO 
c            write(6,'(5(A,I2))')            
c    *       'contribution to X-mean-field for L= ',
c    *       l1,' m-values ',m1,' and ',m2,' from L= ',l3,' M= ',
c    *       m3
             if (bonn.or.breiT) then 
             call two2mean12b(cartexSO(ilocalX),cartexSO(ilocalX),
     *       occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartx(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             else 
             call two2mean12b(cartexSO(ilocalX),cartexOO(ilocalX),
     *       occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartx(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             endif 
             endif 
             endif 
             endif 
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
        endif 
        ilocalX=ilocalX+ncont(iblock)
c       write(6,*) 'finished'
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++      SIGMA Y      ++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        elseif (indx.eq.1.and.indy.eq.0.and.indz.eq.1.and.  
     *  icheckxy(m1,m2,m3,m4)) then  ! X*Z transforms like L_y  (B2) 
cbs     integrals for sigma_y 
        mblocky=mblocky+1
        mcombcart(1,m1,m2,m3,m4,iblock)=2
        mcombcart(2,m1,m2,m3,m4,iblock)=mblocky                               
c       write(6,'(A21,4I3,A7,I4)') 'Integral for sigma_y ',
c    *  m1,m2,m3,m4,' BLOCK ',mblocky
cbs     check if integrals are still fitting on array..
        if ((ilocalY+ncont(iblock)).gt.mxicart) then   
        write(6,*) 'mxicart is too small. Sorry ..'
        stop 
        endif 
        call tosigY(m1,m2,m3,m4,angintSO(iangfirst(iblock)),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),carteYSO(ilocalY),preY,writeTWO,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner) 
c
        if (.not.bonn.and.(.not.breit)) 
     *  call tosigY(m1,m2,m3,m4,angintOO(iangfirst(iblock)),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),carteYOO(ilocalY),preY,writeTWO,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner) 
        if (makemean) then ! generate mean-field-contributions
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
             if (l1.eq.l3.and.l2.eq.l4) then 
             if (m2.eq.m4.and.m1.lt.m3.
     *       and.iabs(m3-m1).eq.1.and.l1.ne.0) then 
c            write(6,'(5(A,I2))')
c    *       'contribution to Y-mean-field for L= ',
c    *       l1,' m-values ',m1,' and ',m3,' from L= ',l2,' M= ',
c    *       m2
             call two2mean13(carteYSO(ilocalY),occup(1,l2),
     *       AOcoeffs(1,1,l2),onecartY(1,1,ipnt(m1+l1+1,m3+l3+1),l1),
     *       ncontrac(l1),ncontrac(l2),noccorb(l2))         
             endif 
             endif 
             if (l1.eq.l2.and.l3.eq.l4) then 
             if (m1.eq.m2.and.l3.ne.0.and.l3.ne.l1) then 
             if (m3.lt.m4.and.iabs(m3-m4).eq.1) then 
cbs   for the "Bonn-approach"   exchange carteYOO by carteYSO 
c            write(6,'(5(A,I2))')
c    *       'contribution to Y-mean-field for L= ',
c    *       l3,' m-values ',m3,' and ',m4,' from L= ',l1,' M= ',
c    *       m1
             if (bonn.or.breiT) then 
             call two2mean34a(carteYSO(ilocalY),carteYSO(ilocalY),
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartY(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             else 
             call two2mean34a(carteYSO(ilocalY),carteYOO(ilocalY),
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartY(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             endif 
             endif 
             if (m3.gt.m4.and.iabs(m3-m4).eq.1) then 
cbs   for the "Bonn-approach"   exchange carteYOO by carteYSO 
c            write(6,'(5(A,I2))')
c    *       'contribution to Y-mean-field for L= ',
c    *       l3,' m-values ',m3,' and ',m4,' from L= ',l1,' M= ',
c    *       m1
             if (bonn.or.breiT) then 
             call two2mean34b(carteYSO(ilocalY),carteYSO(ilocalY),
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartY(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             else 
             call two2mean34b(carteYSO(ilocalY),carteYOO(ilocalY),
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartY(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             endif 
             endif 
             elseif(m3.eq.m4.and.l1.ne.0) then 
             if (m1.lt.m2.and.iabs(m1-m2).eq.1) then 
cbs   for the "Bonn-approach"   exchange carteYOO by carteYSO 
c            write(6,'(5(A,I2))')
c    *       'contribution to Y-mean-field for L= ',
c    *       l1,' m-values ',m1,' and ',m2,' from L= ',l3,' M= ',
c    *       m3
             if (bonn.or.breiT) then 
             call two2mean12a(carteYSO(ilocalY),carteYSO(ilocalY),
     *       occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartY(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             else 
             call two2mean12a(carteYSO(ilocalY),carteYOO(ilocalY),
     *       occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartY(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             endif 
             endif 
             if (m1.gt.m2.anD.Iabs(m1-m2).eq.1) then 
cbs   for the "Bonn-approach"   exchange carteYOO by carteYSO 
c            write(6,'(5(A,I2))')
c    *       'contribution to Y-mean-field for L= ',
c    *       l1,' m-values ',m1,' and ',m2,' from L= ',l3,' M= ',
c    *       m3
             if (bonn.or.breiT) then 
             call two2mean12b(carteYSO(ilocalY),carteYSO(ilocalY),
     *       occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartY(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             else 
             call two2mean12b(carteYSO(ilocalY),carteYOO(ilocalY),
     *       occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartY(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             endif 
             endif 
             endif 
             endif 
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
        endif 
        ilocalY=ilocalY+ncont(iblock)
c       write(6,*) 'finished'
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++      SIGMA Z      ++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        elseif (indx.eq.1.and.indy.eq.1.and.indz.eq.0.and.  
     *  icheckz(m1,m2,m3,m4)) then ! X*Y transforms like L_z  (A2) 
cbs     integrals for sigma_z 
        mblockz=mblockz+1
        mcombcart(1,m1,m2,m3,m4,iblock)=3
        mcombcart(2,m1,m2,m3,m4,iblock)=mblockz                               
c       write(6,'(A21,4I3,A7,I4)') 'Integral for sigma_z ',
c    *  m1,m2,m3,m4,' BLOCK ',mblockz
cbs     check if integrals are still fitting on array..
        if ((ilocalZ+ncont(iblock)).gt.mxicart) then   
        write(6,*) 'mxicart is too small. Sorry ..'
        stop 
        endif 
        call tosigZ(m1,m2,m3,m4,angintSO(iangfirst(iblock)),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),carteZSO(ilocalZ),preXZ,writeTWO,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner) 
c
        if (.not.bonn.and.(.not.breit)) 
     *  call tosigZ(m1,m2,m3,m4,angintOO(iangfirst(iblock)),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),carteZOO(ilocalZ),preXZ,writeTWO,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner) 
        if (makemean) then ! generate mean-field-contributions
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
             if (l1.eq.l3.and.l2.eq.l4) then 
             if (m2.eq.m4.and.m1.lt.m3.
     *       and.m1.eq.-m3.and.l1.ne.0) then 
c            write(6,'(5(A,I2))')
c    *       'contribution to Z-mean-field 13  for L= ',
c    *       l1,' m-values ',m1,' and ',m3,' from L= ',l2,' M= ',
c    *       m2
             call two2mean13(cartezSO(ilocalZ),occup(1,l2),
     *       AOcoeffs(1,1,l2),onecartz(1,1,ipnt(m1+l1+1,m3+l3+1),l1),
     *       ncontrac(l1),ncontrac(l2),noccorb(l2))         
             endif 
             endif 
             if (l1.eq.l2.and.l3.eq.l4) then 
             if (m1.eq.m2.and.l3.ne.0.and.l3.ne.l1) then 
             if (m3.lt.m4.and.m3.eq.-m4) then 
cbs   for the "Bonn-approach"   exchange carteZOO by carteZSO 
c            write(6,'(5(A,I2))')
c    *       'contribution to Z-mean-field 34 for L= ',
c    *       l3,' m-values ',m3,' and ',m4,' from L= ',l1,' M= ',
c    *       m1
             if (bonn.or.breiT) then 
             call two2mean34a(cartezSO(ilocalZ),cartezSO(ilocalZ),
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartz(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             else 
             call two2mean34a(cartezSO(ilocalZ),cartezOO(ilocalZ),
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartz(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             endif 
             endif 
             if (m3.gt.m4.and.m3.eq.-m4) then 
cbs   for the "Bonn-approach"   exchange carteZOO by carteZSO 
c            write(6,'(5(A,I2))')
c    *       'contribution to Z-mean-field 34 for L= ',
c    *       l3,' m-values ',m3,' and ',m4,' from L= ',l1,' M= ',
c    *       m1
             if (bonn.or.breiT) then 
             call two2mean34b(cartezSO(ilocalZ),cartezSO(ilocalZ),
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartz(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             else 
             call two2mean34b(cartezSO(ilocalZ),cartezOO(ilocalZ),
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartz(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             endif 
             endif 
             elseif(m3.eq.m4.and.l1.ne.0) then 
             if (m1.lt.m2.and.m1.eq.-m2) then 
cbs   for the "Bonn-approach"   exchange carteZOO by carteZSO 
c            write(6,'(5(A,I2))')
c    *       'contribution to Z-mean-field 12  for L= ',
c    *       l1,' m-values ',m1,' and ',m2,' from L= ',l3,' M= ',
c    *       m3
             if (bonn.or.breiT) then 
             call two2mean12a(cartezSO(ilocalZ),cartezSO(ilocalZ),
     *       occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartz(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             else 
             call two2mean12a(cartezSO(ilocalZ),cartezOO(ilocalZ),
     *       occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartz(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             endif 
             endif 
             if (m1.gt.m2.and.m1.eq.-m2) then 
cbs   for the "Bonn-approach"   exchange carteZOO by carteZSO 
c            write(6,'(5(A,I2))')
c    *       'contribution to Z-mean-field 12  for L= ',
c    *       l1,' m-values ',m1,' and ',m2,' from L= ',l3,' M= ',
c    *       m3
             if (bonn.or.breiT) then 
             call two2mean12b(cartezSO(ilocalZ),cartezSO(ilocalZ),
     *       occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartz(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             else 
             call two2mean12b(cartezSO(ilocalZ),cartezOO(ilocalZ),
     *       occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartz(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             endif 
             endif 
             endif 
             endif 
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
        endif 
        ilocalZ=ilocalZ+ncont(iblock)
        endif    
c       write(6,*) 'finished'
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        endif ! for check of significance for meanfield.  
        enddo 
        enddo 
        enddo 
        enddo 
        numbcart=numbcart+(mblockx+mblocky+mblockz)*ncont(iblock)
c     write(6,'(A32,4I4,3(A3,I5))') 
c    *'integral blocks for l1,2,l3,l4: ',
c    *l1,l2,l3,l4,' X:',mblockx,' Y:',mblocky,' Z:',mblockz
cbs   just controlling if x and y integrals have the same number of blocks 
      if (mblockx.ne.mblocky) then 
      write(6,*) 
     *'numbers of integrals for sigma_x and sigma_y not equal!'    
      write(6,'(A12,4I3,2(A3,I5))') 
     *'l1,l2,l3,l4 ',l1,l2,l3,l4,' X:',mblockx,' Y:',mblocky  
      write(6,*) ' check the ipowxyz-array'
      stop   
      endif   
cbs   start adresses for the next <ll|ll> block of integrals 
      if (keepcart) then 
      istartX(iblock+1)=istartX(iblock)+mblockX*ncont(iblock)
      istartY(iblock+1)=istartY(iblock)+mblockY*ncont(iblock)
      istartZ(iblock+1)=istartZ(iblock)+mblockZ*ncont(iblock)
      else
      istartX(iblock+1)=istartX(iblock)
      istartY(iblock+1)=istartY(iblock)
      istartZ(iblock+1)=istartZ(iblock)
      endif 
cbs 
cbs   Now take the next  (l,l,l,l)-block   
cbs 
      enddo 
      write(6,*) 'number of cartesian integrals: ',numbcart 
      return 
      end  
