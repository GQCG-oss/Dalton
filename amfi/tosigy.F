      subroutine   tosigY(m1,m2,m3,m4,angint,
     *mcombina,ncontl1,ncontl2,ncontl3,
     *ncontl4,carteY,preY,writeTWO,interxyz,isgnprod,
     *cleaner) 
cbs   this subroutine combines the angular integrals 
cbs   to the integrals for the real-valued linear 
cbs   combinations for the sigma_X part 
cbs   definition of the real-valued linear combinations:
cbs
cbs
cbs   M=0  is the same as   Y(L,0)
cbs
cbs
cbs   M > 0
cbs   
cbs   | L,M,+> = 2**(-0.5) ( (-1)**M Y(L,M) + Y(L,-M)) 
cbs   
cbs   | L,M,-> = -i 2**(-0.5) ( (-1)**M Y(L,M) - Y(L,-M)) ($$$$) 
cbs
cbs
cbs   due to symmetry, there can be only integrals 
cbs   with one or three (sigma_+ and sigma_-)  - combinations 
cbs
      implicit double precision (a-h,o-z)
      include 'para.h'
      parameter (fine=7.29735308D-03) !TO_BE_CHECKED
cbs   at least it's identical with Odd's valuE
      parameter (speed=1d0/fine)
      parameter (speed2=speed*speed)
      dimension mcombina(2,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *angint(ncontl1,ncontl2,ncontl3,ncontl4,*),
cbs  !!!!!!!!!!!changing now to the order of chemists notation!!!!!!!!!!
     *carteY(ncontl1,ncontl3,ncontl2,ncontl4),
     *preY(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *interxyz(*),
     *isgnprod(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *isgnM(-1:1,-1:1,-1:1,-1:1)
      logical*4 writeTWO,cleaner
c     write(6,*) 'begin tosigy '
cbs   cleaning up the integral-array
      irun=ncontl4*ncontl2*ncontl3*ncontl1
      call dzero(irun,carteY)
cbs   set some signs
cbs   isgnM will give an additonal minus-sign if both m-values   
cbs   (cartesian and angular) are negative  see $$$$
      do irun4=-1,1
      do irun3=-1,1
      do irun2=-1,1
      do irun1=-1,1
      isgnM(irun1,irun2,irun3,irun4)=1
      enddo
      enddo
      enddo
      enddo
      if (m1.lt.0) then  
      do irun4=-1,1
      do irun3=-1,1
      do irun2=-1,1
      isgnM(-1,irun2,irun3,irun4)=
     *-isgnM(-1,irun2,irun3,irun4)
      enddo
      enddo
      enddo
      endif 
      if (m2.lt.0) then
      do irun4=-1,1
      do irun3=-1,1
      do irun1=-1,1
      isgnM(irun1,-1,irun3,irun4)=
     *-isgnM(irun1,-1,irun3,irun4)
      enddo
      enddo
      enddo
      endif
      if (m3.lt.0) then 
      do irun4=-1,1
      do irun2=-1,1
      do irun1=-1,1
      isgnM(irun1,irun2,-1,irun4)=
     *-isgnM(irun1,irun2,-1,irun4)
      enddo
      enddo
      enddo
      endif
      if (m4.lt.0) then 
      do irun3=-1,1
      do irun2=-1,1
      do irun1=-1,1
      isgnM(irun1,irun2,irun3,-1)=
     *-isgnM(irun1,irun2,irun3,-1)
      enddo
      enddo
      enddo
      endif
cbs   define absolute m-values
      Mabs1=iabs(m1)
      Mabs2=iabs(m2)
      Mabs3=iabs(m3)
      Mabs4=iabs(m4)
      irun=0
      if (interxyz(1).eq.0) then 
      write(6,*) 'tosigy: no interaction: ',m1,m2,m3,m4
      stop
      endif
      prey1234=preY(m1,m2,m3,m4)
c     write(6,*) 'prey ',prey1234
      do while (interxyz(irun+1).gt.0) 
      irun=irun+1
c     write(6,*) 'tosigy: ',irun,interxyz(irun)
c
cbs
cbs
cbs   This could be done with gotos, but I'm biased to hate those..
cbs
cbs
         if (interxyz(irun).eq.1) then   
         ityp=mcombina(1,Mabs1,Mabs2,Mabs3,Mabs4) 
         if (ityp.ne.1.and.ityp.ne.3) stop 'wrong ityp in tosigY 1'
         iblock=mcombina(2,Mabs1,Mabs2,Mabs3,Mabs4)
         factor=isgnM(1,1,1,1)*prey1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.2) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,-Mabs3,-Mabs4) 
         if (ityp.ne.1.and.ityp.ne.3) stop 'wrong ityp in tosigY 2'
         iblock=mcombina(2,-Mabs1,-Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(-1,-1,-1,-1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,-Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.3) then   
         ityp=mcombina(1,Mabs1,Mabs2,Mabs3,-Mabs4) 
         if (ityp.ne.1.and.ityp.ne.3) stop 'wrong ityp in tosigY 3'
         iblock=mcombina(2,Mabs1,Mabs2,Mabs3,-Mabs4)
         factor=isgnM(1,1,1,-1)*prey1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.4) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,-Mabs3,Mabs4) 
         if (ityp.ne.1.and.ityp.ne.3) stop 'wrong ityp in tosigY 4'
         iblock=mcombina(2,-Mabs1,-Mabs2,-Mabs3,Mabs4)
         factor=isgnM(-1,-1,-1,1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,-Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.5) then   
         ityp=mcombina(1,Mabs1,Mabs2,-Mabs3,Mabs4) 
         if (ityp.ne.1.and.ityp.ne.3) stop 'wrong ityp in tosigY 5'
         iblock=mcombina(2,Mabs1,Mabs2,-Mabs3,Mabs4)
         factor=isgnM(1,1,-1,1)*prey1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,-Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.6) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,Mabs3,-Mabs4) 
         if (ityp.ne.1.and.ityp.ne.3) stop 'wrong ityp in tosigY 6'
         iblock=mcombina(2,-Mabs1,-Mabs2,Mabs3,-Mabs4)
         factor=isgnM(-1,-1,1,-1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.7) then   
         ityp=mcombina(1,Mabs1,-Mabs2,Mabs3,Mabs4) 
         if (ityp.ne.1.and.ityp.ne.3) stop 'wrong ityp in tosigY 7'
         iblock=mcombina(2,Mabs1,-Mabs2,Mabs3,Mabs4)
         factor=isgnM(1,-1,1,1)*prey1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.8) then   
         ityp=mcombina(1,-Mabs1,Mabs2,-Mabs3,-Mabs4) 
         if (ityp.ne.1.and.ityp.ne.3) stop 'wrong ityp in tosigY 8'
         iblock=mcombina(2,-Mabs1,Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(-1,1,-1,-1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,-Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.9) then   
         ityp=mcombina(1,-Mabs1,Mabs2,Mabs3,Mabs4) 
         if (ityp.ne.1.and.ityp.ne.3) stop 'wrong ityp in tosigY 9'
         iblock=mcombina(2,-Mabs1,Mabs2,Mabs3,Mabs4)
         factor=isgnM(-1,1,1,1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.10) then   
         ityp=mcombina(1,Mabs1,-Mabs2,-Mabs3,-Mabs4) 
         if (ityp.ne.1.and.ityp.ne.3) stop 'wrong ityp in tosigY 10'
         iblock=mcombina(2,Mabs1,-Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(1,-1,-1,-1)*prey1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,-Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.11) then  
         ityp=mcombina(1,Mabs1,Mabs2,-Mabs3,-Mabs4) 
         if (ityp.ne.1.and.ityp.ne.3) stop 'wrong ityp in tosigY 11'
         iblock=mcombina(2,Mabs1,Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(1,1,-1,-1)*prey1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,-Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.12) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,Mabs3,Mabs4) 
         if (ityp.ne.1.and.ityp.ne.3) stop 'wrong ityp in tosigY 12'
         iblock=mcombina(2,-Mabs1,-Mabs2,Mabs3,Mabs4)
         factor=isgnM(-1,-1,1,1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.13) then   
         ityp=mcombina(1,Mabs1,-Mabs2,Mabs3,-Mabs4) 
         if (ityp.ne.1.and.ityp.ne.3) stop 'wrong ityp in tosigY 13'
         iblock=mcombina(2,Mabs1,-Mabs2,Mabs3,-Mabs4)
         factor=isgnM(1,-1,1,-1)*prey1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.14) then   
         ityp=mcombina(1,-Mabs1,Mabs2,-Mabs3,Mabs4) 
         if (ityp.ne.1.and.ityp.ne.3) stop 'wrong ityp in tosigY 14'
         iblock=mcombina(2,-Mabs1,Mabs2,-Mabs3,Mabs4)
         factor=isgnM(-1,1,-1,1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,-Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.15) then   
         ityp=mcombina(1,Mabs1,-Mabs2,-Mabs3,Mabs4) 
         if (ityp.ne.1.and.ityp.ne.3) stop 'wrong ityp in tosigY 15'
         iblock=mcombina(2,Mabs1,-Mabs2,-Mabs3,Mabs4)
         factor=isgnM(1,-1,-1,1)*prey1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,-Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.16) then   
         ityp=mcombina(1,-Mabs1,Mabs2,Mabs3,-Mabs4) 
         if (ityp.ne.1.and.ityp.ne.3) stop 'wrong ityp in tosigY 16'
         iblock=mcombina(2,-Mabs1,Mabs2,Mabs3,-Mabs4)
         factor=isgnM(-1,1,1,-1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         endif  
      Enddo 
        if (cleaner) then   
        do irun4=1,ncontl4
        do irun2=1,ncontl2
        do irun1=1,ncontl1
        cartey(irun1,irun1,irun2,irun4)=0d0
        enddo
        enddo
        enddo
        endiF
        if (writeTWO) then
        write(6,*) 'tosigY m1,m2,m3,m4 ',m1,m2,m3,m4
        do irun1=1,ncontl1
        do irun2=1,ncontl2
        do irun3=1,ncontl3
        do irun4=1,ncontl4
        write(6,'(A,4I3,A,E20.14)') 'Indices ',
     *  irun1,irun3,irun2,irun4,'Integral: ',
     *  4d0*speed2*carteY(irun1,irun3,irun2,irun4)
        enddo
        enddo
        enddo
        enddo
        endif
c     write(6,*) 'tosigy:',m1,m2,m3,m4,cartey(1,1,1,1)
c     write(6,*) 'end tosigy '
      return 
      end 
