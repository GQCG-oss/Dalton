      subroutine mkangLmin(Lmax,l1,l2,l3,l4,m1,m2,m3,m4,             
     *angintSO,angintOO,
     *Lfirst,Llast,Lblocks,
     *ncont1,ncont2,ncont3,
     *ncont4,
     *caseaSO,caseb1SO,caseb2SO,casecSO,
     *caseaOO,caseb1OO,caseb2OO,casecOO,
     *preroots,clebsch,dummy,bonn,breit,
     *writeTWO,sameorb)
      implicit double precision (a-h,o-z) 
cbs   subroutine for combining radial intgrls with angular 
cbs   factors for the block with l1,l2,l3,l4,m1,m2,m3m,m4  
cbs   this routine mkangLmin = make angular factors for the L- -part 
cbs   includes both, spin-same and spin-other-orbit parts. 
      double precision LMdepang
      dimension 
     *angintSO(ncont1,ncont2,ncont3,ncont4),
     *angintOO(ncont1,ncont2,ncont3,ncont4),
     *Lfirst(*),Llast(*),Lblocks(*),
cbs   all the arrays with the radial intgrls for 
cbs   this combination of l-values   
     *caseaSO(ncont1*ncont2*ncont3*ncont4,*),  ! (2,0)   intgrls with alpha1*alpha3
     *caseb1SO(ncont1*ncont2*ncont3*ncont4,*), ! (0,0)   intgrls with alpha1
     *caseb2SO(ncont1*ncont2*ncont3*ncont4,*), ! (0,0)   intgrls with alpha3
     *casecSO(ncont1*ncont2*ncont3*ncont4,*),  ! (-2,0)  intgrls with factor 1          
     *caseaOO(ncont1*ncont2*ncont3*ncont4,*),  ! (2,0)   intgrls with alpha1*alpha3
     *caseb1OO(ncont1*ncont2*ncont3*ncont4,*), ! (0,0)   intgrls with alpha1
     *caseb2OO(ncont1*ncont2*ncont3*ncont4,*), ! (0,0)   intgrls with alpha3
     *casecOO(ncont1*ncont2*ncont3*ncont4,*),  ! (-2,0)  intgrls with factor 1          
     *preroots(2,0:Lmax),                    ! some prefactors: dsqrt( (l(+1))/(2l+1)) 
     *clebsch(3,2,-Lmax:Lmax,0:Lmax)         ! some clebsch gordans, that appear regulary 
      dimension dummy(0:*)                                   
      logical*4 writeTWO,bonn,breiT,sameorb 
      root2=dsqrt(2.0d0) 
      root2inv=1d0/root2            
c     write(6,*) 'begin mkangL- ',
c    *l1,l2,l3,l4,m1,m2,m3,m4
cbs  
      ncontall=ncont1*ncont2*ncont3*ncont4
cbs   cheater introduced to correct signs, because they were different from HERMIT
      if (mod(l1+l2+l3+l4,4).eq.2) then
      cheater=1d0
      else
      cheater=-1d0
      endiF
cbs   cleaning up 
      if (bonn.or.breit.or.sameorb) then 
      call dzero(ncontall,angintSO)
      else 
      call dzero(ncontall,angintSO)
      call dzero(ncontall,angintOO)
      endif  
cbs  starting with the same-orbit-contributions 
cbs  first term: ###########################################################################
      factor=-root2inv*preroots(2,l1)*preroots(2,l3)*
     *clebsch(3,2,m1,l1)*
     *clebsch(2,2,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1+1,l2,l3+1,l4,m1+1,m2,m3,m4,cheater)
      if (dummy(L).ne.0d0) then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      Endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs   second term: ###########################################################################
      factor=-root2inv*preroots(1,l1)*preroots(2,l3)*
     *clebsch(3,1,m1,l1)*
     *clebsch(2,2,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1+1,m2,m3,m4,cheater) 
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(3).gt.0) then 
      M=m2-m4
	      if (Lfirst(3).lt.Kfirst) then 
	      do L=Lfirst(3),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1+1,m2,m3,m4,cheater) 
	      enddo 
	      Kfirst=Lfirst(3)
	      endif 
	      if (Llast(3).gt.Klast) then 
	      do L=Klast,Llast(3),2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1+1,m2,m3,m4,cheater)
	      enddo 
	      Klast=Llast(3)
	      endif 
      Lrun=1
      do L=Lfirst(3),Llast(3),2  
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,-(2+4*l1)*
     *   factor*dummy(L),caseb2OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs   third term: ###########################################################################
      factor=-root2inv*preroots(2,l1)*preroots(1,l3)*
     *clebsch(3,2,m1,l1)*
     *clebsch(2,1,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      if (dummy(L).ne.0d0) then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(2).gt.0) then 
      M=m2-m4
	      if (Lfirst(2).lt.Kfirst) then 
	      do L=Lfirst(2),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1+1,m2,
     *                 m3,m4,Cheater)
	      enddo 
	      Kfirst=Lfirst(2)
	      endif 
	      if (Llast(2).gt.Klast) then 
	      do L=Klast,Llast(2),2
	      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
	      enddo 
	      Klast=Llast(2)
	      endif 
      Lrun=1
      do L=Lfirst(2),Llast(2),2  
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
      else   
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,-(2+4*l3)*
     *   factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs   fourth term: ###########################################################################
      factor=-root2inv*preroots(1,l1)*preroots(1,l3)*
     *clebsch(3,1,m1,l1)*
     *clebsch(2,1,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
      if (dummy(L).ne.0d0) then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(2).gt.0) then 
      M=m2-m4
	      if (Lfirst(2).lt.Kfirst) then 
	      do L=Lfirst(2),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
	      enddo 
	      Kfirst=Lfirst(2)
	      endif 
	      if (Llast(2).gt.Klast) then 
	      do L=Klast,Llast(2),2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
	      enddo 
	      Klast=Llast(2)
	      endif 
      Lrun=1
      do L=Lfirst(2),Llast(2),2  
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,-(2+4*l3)*
     *   factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif  
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(3).gt.0) then 
      M=m2-m4
	      if (Lfirst(3).lt.Kfirst) then 
	      do L=Lfirst(3),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
	      enddo 
	      Kfirst=Lfirst(3)
	      endif 
	      if (Llast(3).gt.Klast) then 
	      do L=Klast,Llast(3),2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
	      enddo 
	      Klast=Llast(3)
	      endif 
      Lrun=1
      do L=Lfirst(3),Llast(3),2  
      if (dummy(L).ne.0d0) then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,-(2+4*l1)*
     *   factor*dummy(L),caseb2OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(4).gt.0) then 
      M=m2-m4
	      if (Lfirst(4).lt.Kfirst) then 
	      do L=Lfirst(4),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
	      enddo 
	      Kfirst=Lfirst(4)
	      endif 
	      if (Llast(4).gt.Klast) then 
	      do L=Klast,Llast(4),2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3,m4,cheater)
	      enddo 
	      Klast=Llast(4)   
	      endif 
      Lrun=1
      do L=Lfirst(4),Llast(4),2  
      if (dummy(L).ne.0d0) then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,(4*l1*l3+2*l1+2*l3+1)*factor*dummy(L),
     *   casecSO(1,Lrun),1,angintSO,1)
      else  
         call daxpy(ncontall,(4*l1*l3+2*l1+2*l3+1)*factor*dummy(L),
     *   casecSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   (4*l1*l3+2*l1+2*l3+1)*factor*dummy(L),
     *   casecOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs  fifth term: ###########################################################################
      factor=-root2inv*preroots(2,l1)*preroots(2,l3)*
     *clebsch(2,2,m1,l1)*
     *clebsch(1,2,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1+1,l2,l3+1,l4,m1,m2,m3-1,m4,cheater)
      if (dummy(L).ne.0d0)  then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo
      endif 
      endif 
cbs   sixth  term: ###########################################################################
      factor=-root2inv*preroots(1,l1)*preroots(2,l3)*
     *clebsch(2,1,m1,l1)*
     *clebsch(1,2,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1,m2,m3-1,m4,cheater)
      if (dummy(L).ne.0d0) then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,4*
     *   factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(3).gt.0) then 
      M=m2-m4
	      if (Lfirst(3).lt.Kfirst) then 
	      do L=Lfirst(3),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1,m2,m3-1,m4,cheater)
	      enddo 
	      Kfirst=Lfirst(3)
	      endif 
	      if (Llast(3).gt.Klast) then 
	      do L=Klast,Llast(3),2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1,m2,m3-1,m4,cheater)
	      enddo 
	      Klast=Llast(3)
	      endif 
      Lrun=1
      do L=Lfirst(3),Llast(3),2  
      if (dummy(L).ne.0d0)  then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,-(2+4*l1)*
     *   factor*dummy(L),caseb2OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs   seventh term: ###########################################################################
      factor=-root2inv*preroots(2,l1)*preroots(1,l3)*
     *clebsch(2,2,m1,l1)*
     *clebsch(1,1,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo
      endif 
      if (Lblocks(2).gt.0) then 
      M=m2-m4
	      if (Lfirst(2).lt.Kfirst) then 
	      do L=Lfirst(2),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
	      enddo 
	      Kfirst=Lfirst(2)
	      endif 
	      if (Llast(2).gt.Klast) then 
	      do L=Klast,Llast(2),2
	      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
	      enddo 
	      Klast=Llast(2)
	      endif 
      Lrun=1
      do L=Lfirst(2),Llast(2),2  
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,-(2+4*l3)*
     *   factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs   eigth term: ###########################################################################
      factor=-root2inv*preroots(1,l1)*preroots(1,l3)*
     *clebsch(2,1,m1,l1)*
     *clebsch(1,1,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      Klast=0
      Kfirst=Lmax+Lmax+1 ! just to be sure ..
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Kfirst=Lfirst(1)
      Klast=Llast(1)
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
      if (dummy(L).ne.0d0)  then
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,4*
     *   factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(2).gt.0) then 
      M=m2-m4
	      if (Lfirst(2).lt.Kfirst) then 
	      do L=Lfirst(2),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
	      enddo 
	      Kfirst=Lfirst(2)
	      endif 
	      if (Llast(2).gt.Klast) then 
	      do L=Klast,Llast(2),2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
	      enddo 
	      Klast=Llast(2)
	      endif 
      Lrun=1
      do L=Lfirst(2),Llast(2),2  
      if (dummy(L).ne.0d0) then 
      if (bonn.or.breit.or.sameorb) then 
        call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *  caseb1SO(1,Lrun),1,angintSO,1)
      else 
        call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *  caseb1SO(1,Lrun),1,angintSO,1)
        call daxpy(ncontall,-(2+4*l3)*
     *factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(3).gt.0) then 
      M=m2-m4
	      if (Lfirst(3).lt.Kfirst) then 
	      do L=Lfirst(3),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
	      enddo 
	      Kfirst=Lfirst(3)
	      endif 
	      if (Llast(3).gt.Klast) then 
	      do L=Klast,Llast(3),2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
	      enddo 
	      Klast=Llast(3)
	      endif 
      Lrun=1
      do L=Lfirst(3),Llast(3),2  
      if (dummy(L).ne.0d0)  then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,-(2+4*l1)*
     *   factor*dummy(L),caseb2OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(4).gt.0) then 
      M=m2-m4
	      if (Lfirst(4).lt.Kfirst) then 
	      do L=Lfirst(4),Kfirst,2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
	      enddo 
	      Kfirst=Lfirst(4)
	      endif 
	      if (Llast(4).gt.Klast) then 
	      do L=Klast,Llast(4),2
	      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1,m2,m3-1,m4,cheater)
	      enddo 
	      Klast=Llast(4)   
	      endif 
      Lrun=1
      do L=Lfirst(4),Llast(4),2  
      if (dummy(L).ne.0d0) then 
      if (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,(4*l1*l3+2*l1+2*l3+1)*
     *   factor*dummy(L),
     *   casecSO(1,Lrun),1,angintSO,1)
      else   
         call daxpy(ncontall,(4*l1*l3+2*l1+2*l3+1)*
     *   factor*dummy(L),
     *   casecSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   (4*l1*l3+2*l1+2*l3+1)*factor*dummy(L),
     *   casecOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo
      endif 
      endif 
      if (writeTWO) then
      write(6,'(A,4I3)') 
     *'L-: SO-Integrals in spherical harmonics for m-values '
     *,m1,m2,m3,m4
      do irun1=1,ncont1
      do irun2=1,ncont2
      do irun3=1,ncont3
      do irun4=1,ncont4
      write(6,'(A,4I3,A,E20.14)')
     *'indices ',irun1,irun2,irun3,irun4,' integral ',
     *angintSO(irun1,irun2,irun3,irun4)
      enddo
      enddo
      enddo
      enddo
      if ((.not.bonn).and.(.not.breit).and.(.not.sameorb)) then 
      write(6,'(A,4I3)')
     *'L-: OO-Integrals in spherical harmonics for m-values '
     *,m1,m2,m3,m4
      do irun1=1,ncont1
      do irun2=1,ncont2
      do irun3=1,ncont3
      do irun4=1,ncont4
      write(6,'(A,4I3,A,E20.14)')
     *'indices ',irun1,irun2,irun3,irun4,' integral ',
     *angintOO(irun1,irun2,irun3,irun4)
      enddo
      enddo
      enddo
      enddo
      endif
      endif
      return   
      end
   
