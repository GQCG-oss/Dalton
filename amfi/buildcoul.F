      subroutine buildcoul(l1,l2,l3,l4,! angular momenta of primitives
     *incl1,incl3, ! shifts for different radial integrals
     *Lrun, ! L-value for coulomb integrals 
     *prmints,
     *nprim1,nprim2,nprim3,nprim4,  ! number of primitives
     *expo1,expo2,expo3,expo4, ! arrays with the exponents
     *power13,
     *power24,
     *quotpow1,quotpow2,writeTWO
     *)
cbs ##################################################################
c
cbs  purpose: builds up the coulomb integrals 
cbs  inbetween primitives and multiplies 
cbs  with extra factors to correct the 
cbs  normalization                  
c
cbs ##################################################################
      implicit double precision (a-h,o-z) 
      include 'param.h'
      include 'dofuc.h'
      logical*4 writeTWO    
      common /pifact/ pi,root8ovpi
      dimension expo1(nprim1),
     *expo2(nprim2),
     *expo3(nprim3),
     *expo4(nprim4), ! the exponents
     *prmints(nprim1,nprim2,nprim3,nprim4), ! scratch array for integrals over primitives
     *power13(MxprimL,MxprimL),
     *power24(MxprimL,MxprimL),
     *quotpow1(nprim1,nprim2,nprim3,nprim4),
     *quotpow2(nprim1,nprim2,nprim3,nprim4),
     *fraclist1(0:Lmax+3),fraclist2(0:Lmax+3),fact(MxprimL),
     *frac(MxprimL),cfunctx1(MxprimL),cfunctx2(MxprimL)
c     Write(6,*) 'begin buildcoul' 
c     write(6,'(14I3)') l1,l2,l3,l4,incl1,incl3,
c    *nprim1,nprim2,nprim3,nprim4,Lrun   
cbs ##################################################################
cbs        prepare indices for coulint
cbs ##################################################################
      n1=l1+incl1+1
      n2=l2+1
      n3=l3+incl3+1
      n4=l4+1
c     write(6,*) 'n1,n4 ',n1,n2,n3,n4
      n13=n1+n3
      n24=n2+n4
      index1=N13-Lrun-1
      index2=n24+Lrun
      index3=N24-Lrun-1
      index4=n13+Lrun
c     write(6,*) 'indices ',index1,index2,index3,index4
      do krun=0,(index1-1)/2
      fraclist1(krun)=dffrac(krun+krun+index2-1,krun+krun)*
     *dffrac(1,index2-1)
      enddo
      do krun=0,(index3-1)/2
      fraclist2(krun)=dffrac(krun+krun+index4-1,krun+krun)*
     *dffrac(1,index4-1)
      enddo
cbs ##################################################################
cbs   common factors including double factorials 
cbs ##################################################################
      doff1=dffrac(index1-1,n13-1)*dffrac(n24+Lrun-1,n24-1)
      doff2=dffrac(index3-1,n24-1)*dffrac(n13+Lrun-1,n13-1)
      if (index1.eq.1) then 
                do irun4=1,nprim4
                do irun3=1,nprim3
                if (l2.eq.l4) then          
                limit2=irun4
                else
                limit2=nprim2
                endif 
  	        do irun2=1,limit2 
   	        pow24inv=doff1/power24(irun4,irun2)
                if (l1.eq.l3) then          
                limit1=irun3
                else
                limit1=nprim1
                endif 
        	do irun1=1,limit1 
                prmints(irun1,irun2,irun3,irun4)=    
     *          quotpow1(irun1,irun2,irun3,irun4)*          
     *          dsqrt(0.5d0*(expo1(irun1)+expo3(irun3)))*
     *          power13(irun3,irun1)*pow24inv
      		enddo 
      		enddo 
      		enddo 
      		enddo 
      else 
                do irun4=1,nprim4
                do irun3=1,nprim3
                if (l2.eq.l4) then          
                limit2=irun4
                else
                limit2=nprim2
                endif 
  	        do irun2=1,limit2 
                alpha24inv=1d0/(expo2(irun2)+expo4(irun4))
                pow24inv=doff1/power24(irun4,irun2)
                if (l1.eq.l3) then          
                limit1=irun3
                else
                limit1=nprim1
                endif 
        	do irun1=1,limit1 
                a1324= alpha24inv*(expo1(irun1)+expo3(irun3))
c               write(6,*) 'a1324 ',a1324
                   Cfunctx1(irun1)=fraclist1(0)
                   frac(irun1)=a1324/(1d0+a1324)
c                  write(6,*) 'frac(1) ',frac(1)  
                   fact(irun1)=frac(irun1) 
                enddo
*vocl    loop,repeat(Lmax+3)
                   do k=1,(index1-1)/2
                   do irun1=1,limit1  
                   Cfunctx1(irun1)=Cfunctx1(irun1)+fraclist1(k)
     *            *fact(irun1)
                   enddo
                   do irun1=1,limit1  
                   fact(irun1)=fact(irun1)*frac(irun1)
                   enddo
                   enddo
                do irun1=1,limit1  
                alpha13=0.5d0*(expo1(irun1)+expo3(irun3))
                prmints(irun1,irun2,irun3,irun4)=    
     *          quotpow1(irun1,irun2,irun3,irun4)*       
     *          dsqrt(alpha13)*power13(irun3,irun1)*pow24inv*
     *          Cfunctx1(irun1) 
		enddo
		enddo
		enddo
		enddo
      endif   
      if (index3.eq.1) then 
                do irun4=1,nprim4
                do irun3=1,nprim3
                if (l2.eq.l4) then          
                limit2=irun4
                else
                limit2=nprim2
                endif 
  	        do irun2=1,limit2 
		pow24=doff2*power24(irun4,irun2)*
     *          dsqrt(0.5d0*(expo2(irun2)+expo4(irun4)))
                if (l1.eq.l3) then          
                limit1=irun3
                else
                limit1=nprim1
                endif 
        	do irun1=1,limit1 
                prmints(irun1,irun2,irun3,irun4)=    
     *          prmints(irun1,irun2,irun3,irun4)+    
     *          pow24*quotpow2(irun1,irun2,irun3,irun4)/
     *          power13(irun3,irun1)
                enddo 
                enddo 
                enddo 
                enddo 
c               write(6,*) 'second ',prmints(1,1,1,1)
      else
                do irun4=1,nprim4
                do irun3=1,nprim3
                if (l2.eq.l4) then          
                limit2=irun4
                else
                limit2=nprim2
                endif 
  	        do irun2=1,limit2 
                alpha24=expo2(irun2)+expo4(irun4)
		pow24=doff2*power24(irun4,irun2)*
     *          dsqrt(0.5d0*alpha24)                                  
                if (l1.eq.l3) then          
                limit1=irun3
                else
                limit1=nprim1
                endif 
        	do irun1=1,limit1 
                a2413= alpha24/(expo1(irun1)+expo3(irun3))
                   Cfunctx2(irun1)=fraclist2(0)
                   frac(irun1)=a2413/(1d0+a2413)
                   fact(irun1)=frac(irun1)
        	enddo 
*vocl    loop,repeat(Lmax+3)
                   do k=1,(index3-1)/2
                   do irun1=1,limit1 
                   Cfunctx2(irun1)=Cfunctx2(irun1)+
     *             fraclist2(k)*fact(irun1)
                   enddo
                   do irun1=1,limit1 
                   fact(irun1)=fact(irun1)*frac(irun1)
                   enddo
                   enddo
c               write(6,*) 'cfunctx2 ',cfunctx2(1),
c    *          quotpow2(1,1,1,1)
                do irun1=1,limit1 
                prmints(irun1,irun2,irun3,irun4)=    
     *          prmints(irun1,irun2,irun3,irun4)+    
     *          quotpow2(irun1,irun2,irun3,irun4)*
     *          Cfunctx2(irun1)*
     *          pow24/power13(irun3,irun1)
                enddo 
                enddo 
                enddo 
                enddo 
      endif 
cbs   make some mirroring for identical l-values
cbs   for the case that l1=l3 
      if (l1.eq.l3) then 
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=irun3+1,nprim1
      prmints(irun1,irun2,irun3,irun4)=    
     *prmints(irun3,irun2,irun1,irun4)
      enddo 
      enddo 
      enddo 
      enddo 
      endif  
cbs   for the case that l2=l4 
      if (l2.eq.l4) then 
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=irun4+1,nprim2
      do irun1=1,nprim1
      prmints(irun1,irun2,irun3,irun4)=    
     *prmints(irun1,irun4,irun3,irun2)
      enddo 
      enddo 
      enddo 
      enddo 
      endif 
cbs   some factors which are the same for all cases 
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=1,nprim1
      prmints(irun1,irun2,irun3,irun4)=    
     *prmints(irun1,irun2,irun3,irun4)*
     *coulovlp(irun4,irun2,0,0,l4,l2)*
     *coulovlp(irun3,irun1,incl3,incl1,l3,l1)*
     *root8ovpi
      enddo 
      enddo 
      enddo 
      enddo 
c     write(6,*) 'buildcoul ',l1+incl1,l2,l3+incl3,l4,
c    *prmints(1,1,1,1) 
cbs   
cbs  look for additional factors, as the 
cbs  coulomb integrals are calculated 
cbs  for normalized functions with that 
cbs  specific l 
cbs  
cbs  if l was increased by one, the factor is
cbs  0.5*dsqrt((2l+3)/(exponent))
cbs  if l was decreased by one, the factor is
cbs  2d0*dsqrt(exponent/(2l+1))
cbs  
cbs
cbs   check for first function 
cbs
cbs  
c     write(6,*) 'index 1'
      if (incl1.eq.1) then 
      fact1=0.5d0*dsqrt(dfloat(l1+l1+3))
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=1,nprim1
      factor=fact1/dsqrt(expo1(irun1))                
      prmints(irun1,irun2,irun3,irun4)=
     *prmints(irun1,irun2,irun3,irun4)*factor                                 
      enddo 
      enddo 
      enddo 
      enddo 
      elseif (incl1.eq.-1) then 
      fact1=2d0/dsqrt(dfloat(l1+l1+1))
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=1,nprim1
      factor=fact1*dsqrt(expo1(irun1))                
      prmints(irun1,irun2,irun3,irun4)=
     *prmints(irun1,irun2,irun3,irun4)*factor                                 
      enddo 
      enddo 
      enddo 
      enddo 
      endif
cbs  
cbs
cbs   check for third function 
cbs
cbs  
c     write(6,*) 'index 3'
      if (incl3.eq.1) then 
      fact1=0.5d0*dsqrt(dfloat(l3+l3+3))
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=1,nprim1
      factor=fact1/dsqrt(expo3(irun3))                
      prmints(irun1,irun2,irun3,irun4)=
     *prmints(irun1,irun2,irun3,irun4)*factor                                 
      enddo 
      enddo 
      enddo 
      enddo 
      elseif (incl3.eq.-1) then 
      fact1=2d0/dsqrt(dfloat(l3+l3+1))
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=1,nprim1
      factor=fact1*dsqrt(expo3(irun3))                
      prmints(irun1,irun2,irun3,irun4)=
     *prmints(irun1,irun2,irun3,irun4)*factor                                 
      enddo 
      enddo 
      enddo 
      enddo 
      endif
c     if (writeTWO) then 
c     do irun1=1,nprim1 
c     do irun2=1,nprim2 
c     do irun3=1,nprim3 
c     do irun4=1,nprim4 
c     write(6,'(A,4I3,A,E20.14)') 'Indices ',irun1,irun2,
c    *irun3,irun4,' integral ',prmints(irun1,irun2,irun3,irun4) 
c     enddo
c     enddo
c     enddo
c     enddo
c     endif 
c     Write(6,*) 'end buildcoul' 
      return 
      end 
