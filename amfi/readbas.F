      subroutine readbas(Lhigh,makemean,bonn,breit,writetwo,
     *HERMIT,MOLCAS,symmetry,sameorb,AIMP) 
cbs   suposed to read the maximum of l-values, the number of primitive and contracted 
cbs   functions, the exponents and contraction coefficients
      implicit double precision (a-h,o-z) 
      include 'param.h' 
      include 'ired.h'
      character*76 Header 
      character*6 setBreit                  
      character*5 setbonn,readsame,readaimp,readcore                   
      character*9 write2el                   
      character*7 setHERMIT
      character*4 symmetry   
      character*13 Llimit 
      character*11 inttype   
      character*13 Llimit2 
      character*19 chcharge,chcharge2
      character*30 Nofprim,Nofprim2 
      character*28 addtext,addtext2 
      character*32 Nofcont,Nofcont2 
      character*76 Stars                  
      logical*4 makemean,bonn,breit,writetwo,HERMIT,
     *sameorb,MOLCAS,AIMP 
      common /nucleus/ charge 
      dimension numboffunct(Mxcart),numboftype(Mxcart)
      Integer*4 OUT,ibeginIRED(8),idelpersym(8) 
      IN=INPFILE    
      OUT=6
      stars='*********************************************************** 
     *******************'
      Llimit='MAX. L-VALUE:'
      chcharge=' CHARGE OF NUCLEUS:'
      Nofprim='NUMBER OF PRIMITIVE FUNCTIONS:' 
      Nofcont=' NUMBER OF CONTRACTED FUNCTIONS:' 
      addtext='ADDITIONAL FUNCTIONS in IRS:'
      write(OUT,*) stars                                                
      write(OUT,*) 'ATOMIC NO-PAIR INTEGRAL AND MEAN-FIELD CODE'
      write(OUT,*) ' developed and coded by Bernd Schimmelpfennig'
      write(OUT,*) stars                                                
C     read(IN,'(A76)') header
C     write(OUT,'(A76)') header
C     read(IN,'(A11,A5,A6,A9,A7,A4,A5,A5)') inttype,setbonn,setbreit,
C    *write2el,SetHERMIT,symmetry,readsame,readaimp   
cbs   set keywords to default for DALTON 
      inttype='MEAN-FIELD '
      setbreit='BREIT ' 
      write2el='writetwo ' 
      setHERMIT='HERMIT '
      symmetry='D2H '
      readsame= '     ' 
      readaimp= '     ' 
      if (inttype.eq.'MEAN-FIELD ') then 
      makemean=.true.
      write(OUT,*) stars                                                
      write(OUT,*) '2e-integrals for the mean-field only'
      write(OUT,*) '    mean-field will be generated         ' 
      write(OUT,*) stars                                                
      else 
      makemean=.false.
      write(OUT,*) stars                                                
      write(OUT,*) 'all two-electron-integrals are calculated'
      write(OUT,*) '       for other purposes  '
      write(OUT,*) stars                                                
      endif 
      if (setbonn.eq.'BONN ') then
      write(OUT,*) 'Bonn-approach for spin-other-orbit part'
      bonn=.TRUE.
      else
      bonn=.FALSE. 
      endif 
      if (setbreit.eq.'BREIT ') then
      write(OUT,*) ' Breit-Pauli-Approximation'
      breit=.TRUE.
      else
      breit=.FALSE. 
      write(OUT,*) 'Douglas-Kroll type operators '
      endif 
      if (write2el.eq.'WRITETWO ') then
      write(OUT,*) stars                                                
      write(OUT,*) 'write out two-electron integrals for tests'  
      write(OUT,*) ' generates HUGE output !!!!'  
      write(OUT,*) stars                                                
      writetwo=.TRUE.
      else
      writetwo=.FALSE. 
      endif 
      if (setHERMIT.eq.'HERMIT ') then
      write(OUT,*) stars                                                
      write(OUT,*) 'write out one-electron integrals in HERMIT-style'
      write(OUT,*) '   and with HERMIT normalization ' 
      write(OUT,*) stars                                                
      HERMIT=.TRUE.
      MOLCAS=.FALSE. 
      elseif (setHERMIT.eq.'MOLCAS ') then  
      write(OUT,*) stars                                                
      write(OUT,*) 'write out one-electron integrals in MOLCAS-style'
      write(OUT,*) '   and with MOLCAS normalization ' 
      write(OUT,*) stars                                                
      HERMIT=.TRUE.
      MOLCAS=.TRUE.
      else
      write(OUT,*) stars                                                
      write(OUT,*) 
     *'write out 1e(MF)-integrals assuming normalization to one'
      write(OUT,*) stars                                                
      HERMIT=.FALSE. 
      MOLCAS=.FALSE. 
      endif 
      if (symmetry.eq.'D2H ') then
      write(OUT,*) 'Symmetry is D2H'
      write(OUT,*) 'check whether order of IRs is correct!!!'
      numbofsym=8    
      elseif (symmetry.eq.'C2V ') then 
      write(OUT,*) 'Symmetry is C2V'
      write(OUT,*) 'check whether order of IRs is correct!!!'
      numbofsym=4    
      elseif (symmetry.eq.'D2  ') then 
      write(OUT,*) 'Symmetry is D2 '
      write(OUT,*) 'check whether order of IRs is correct!!!'
      numbofsym=4    
      elseif (symmetry.eq.'C2H ') then 
      write(OUT,*) 'Symmetry is C2H'
      write(OUT,*) 'check whether order of IRs is correct!!!'
      numbofsym=4    
      elseif (symmetry.eq.'CS  ') then 
      write(OUT,*) 'Symmetry is CS'
      write(OUT,*) 'check whether order of IRs is correct!!!'
      numbofsym=2    
      elseif (symmetry.eq.'    ') then 
      write(OUT,*) 'assumed symmetry is D2H' 
      write(OUT,*) 'check whether order of IRs is correct!!!'
      symmetry='D2H '
      numbofsym=8    
      else
      write(OUT,*) 'Symmetry ',symmetry, 'is not supported so far'
      stop 'Sorry, choose another symmetry or include it into code!!'
      endif 
      if (readsame.eq.'SAME ') then 
      write(OUT,*) 'SAME-ORBIT only'
      sameorb=.true.
      else
      write(OUT,*) 'OTHER-ORBIT included' 
      sameorb=.false.
      endif 
      if (readaimp.eq.'AIMP ') then 
      write(OUT,*) 'CORE removed for use with AIMP' 
      AIMP=.true.
      else
      AIMP=.false.
      endif 
cbs   read(IN,'(A28,8I5)') addtext2,(iadditional(I),I=1,numbofsym)
cbs   if (addtext2.ne.addtext) then  
cbs   stop 'wrong keyword for additional functions ' 
cbs   else 
cbs   write(OUT,'(A,8I4)') 'Additonal functions for the symmetries: ',
cbs  *(iadditional(I),I=1,numbofsym)
cbs   endif 
      read(IN,'(A13,I5,A19,F5.2)') Llimit2,Lhigh,chcharge2,charge     
      if (Llimit2.ne.Llimit) stop 'wrong keyword in readbas' 
      if (chcharge2.ne.chcharge) stop 'wrong keyword in readbas' 
      if (Lhigh.gt.Lmax) then 
      write(OUT,*) 'Sorry, so far this code deals only ',
     *'with maximum l-values of ',Lmax
      stop 
      endif  
      write(OUT,*) ' Functions will go up to an L-value of : ',Lhigh
      write(OUT,'(A19,F5.2)') chcharge,charge    
      if (AIMP) then 
      	read(IN,'(A5,10I5)') readcore,(icore(I),I=0,Lmax) 
      	if (readcore.ne.'CORE:') stop 'wrong keyword to read core' 
      write(6,*) 
      	write(6,*) 'CORE to be deleted '  
      	write(6,*) '   L   #orbs.  '
      	do i=0,Lmax
      	write(6,'(2I5)') i,icore(i)
      	enddo 
      write(6,*) 
        else 
        do i=0,Lmax 
        icore(i)=0
        enddo   
      endif 
      
      call initired(symmetry)
      Do iredrun=1,numbofsym
      do Lrun=0,Lhigh  
      nmbMperIRL(iredrun,Lrun)=0 
      enddo 
      enddo 
      do Lrun=0,Lhigh 
      write(OUT,*) 'ANGULAR MOMENTUM ',LRUN  
      read(IN,'(A30,I5,A32,I5)') Nofprim2,nprimit(Lrun),
     *Nofcont2,ncontrac(Lrun) 
      write(OUT,'(A30,I3,A32,I3)') Nofprim2,nprimit(Lrun),
     *Nofcont2,ncontrac(Lrun) 
cbs   check keywords
      if (Nofprim.ne.Nofprim2.or.Nofcont.ne.Nofcont2)  then 
      write(OUT,*) 'wrong keywords read for information on number',
     *' of functions for L= ',Lrun 
      stop 
      endif 
cbs   check maximum numbers
      if (nprimit(Lrun).gt.MxprimL) then 
      write(OUT,*) 'To many primitives for L=',L,
     *' increase MxprimL in para.h or reduce ',
     *'the number of primitives to at most ',MxprimL
      stop 
      endif   
      if (ncontrac(Lrun).gt.MxcontL) then 
      write(OUT,*) 'To many contracted fncts for L=',L,
     *' increase MxcontL in para.h or ',
     *'reduce the number of contracted functions',
     *'to at most ',MxcontL
      stop
      endif   
      if (ncontrac(Lrun).gt.nprimit(Lrun)) then 
      write(OUT,*) 'You have more contracted than ',
     *'uncontracted functions, I don''t believe ',
     *'that. Sorry!! '
      stop 
      endif
c     write(OUT,'(A7,I3,A15,I3,A33,I3,A24)') 'For L= ',Lrun,
c    *' there will be ',
c    *ncontrac(Lrun),' contracted functions, built from ',
c    *nprimit(Lrun),
c    *' uncontracted functions.' 
      if (MOLCAS) then 
cbs   read input in MOLCAS-style
	      do ILINE=1,nprimit(Lrun)
	      read(IN,*) exponents(ILINE,Lrun) 
	      enddo 
	      do ILINE=1,nprimit(Lrun)
	      read(IN,*) (cntscrtch(ILINE,JRUN,Lrun),Jrun=1,
     *        ncontrac(Lrun))                    
	      enddo 
      else 
cbs   old input-type for HERMIT 
cbs   read exponents and contraction coefficients UNFORMATTED.  
cbs   First BLOCK:
cbs   nprimit(Lrun) times: 
cbs   exponent    up to ncffline contraction coefficients  
cbs   empty line 
cbs   if there are more than ncffline contracted functions, there 
cbs   will be another block: 
cbs   nprimit(Lrun) times: 
cbs   up to ncffline contraction coefficients  
CBS   empty line 
cbs   and so on 
      Irest=ncontrac(Lrun)
      ibegin=1       
      iend=ncffline
      do ILINE=1,nprimit(Lrun)
      read(IN,*) exponents(ILINE,Lrun), 
     *(cntscrtch(ILINE,JRUN,Lrun),Jrun=ibegin,min(iend,ncontrac(Lrun)))
      enddo 
      read(IN,'(A76)') header
      Irest=Irest-ncffline           
c
 4711 if (Irest.gt.0)  then 
	      ibegin=ibegin+ncffline
	      iend=iend+ncffline
	      do ILINE=1,nprimit(Lrun)
	      read(IN,*) 
     *        (cntscrtch(ILINE,JRUN,Lrun),Jrun=ibegin,
     *        min(iend,ncontrac(Lrun)))
              enddo 
              read(IN,'(A76)') header
              Irest=Irest-ncffline            
      endif 
      if (Irest.gt.0) goto 4711 
      endif  
c
cbs   
cbs   end of reading for the current L-value 
cbs   
      do  Irun=1,ncontrac(Lrun)
      writE(OUT,*) 'orbital : ',irun   
      write(OUT,'(6(X,E12.6))') 
     *(cntscrtch(I,Irun,Lrun),I=1,nprimit(Lrun))
      enddo 
      write(OUT,*) ' ' 
cbs   setting the numbers of cartesians per IR
      do iredrun=1,numbofsym 
      nfunctions(iredrun,Lrun)=0
      enddo
      do mrun=-Lrun,Lrun
      nfunctions(ipow2ired(ipowxyz(1,mrun,Lrun),ipowxyz(2,mrun,Lrun),
     *Ipowxyz(3,mrun,Lrun)),Lrun)=
     *nfunctions(ipow2ired(ipowxyz(1,mrun,Lrun),ipowxyz(2,mrun,Lrun),
     *ipowxyz(3,mrun,Lrun)),Lrun)+ncontrac(Lrun)
      enddo
      do mrun=-Lrun,Lrun
      nmbMperIRL(ipow2ired(ipowxyz(1,mrun,Lrun),ipowxyz(2,mrun,Lrun),
     *Ipowxyz(3,mrun,Lrun)),lruN)=
     *nmbMperIRL(ipOw2ired(ipowxyz(1,mrun,Lrun),ipowxyz(2,mrun,Lrun),
     *IpowxYz(3,mrun,Lrun)),lruN)+1
      enddo
      write(OUT,*) stars                                                          
      write(OUT,'(A,8I4)') 
     *'Number of functions per IR: ',(nfunctions(iredrun,Lrun),
     *iredrun=1,numbofsym)
      write(OUT,*) stars                                                          
      enddo   ! enddo for loop over L-values 
C     write(OUT,*) 'distribution of M-values'
c     do Lrun=0,Lhigh
c     write(OUT,*) (nmbMperIRL(nsym,Lrun),nsym=1,numbofsym)
c     endDo 
      numbofcart=0
      do lrun=0,Lhigh 
      numbofcart=numbofcart+(Lrun+Lrun+1)*
     *ncontrac(Lrun)
      enddo
      if (numbofcart.gt.Mxcart) 
     *stop 'sorry, maximum of cartesian functions is exceeded'
      do iredrun=1,numbofsym 
      nfunctperIRED(iredrun)=0
      enddo
      do Lrun=0,Lhigh
      do iredrun=1,numbofsym 
      nfunctperIRED(iredrun)=nfunctperIRED(iredrun)+
     *nfunctions(iredrun,Lrun)
      enddo
      enddo
      write(OUT,*) stars                                                          
      write(OUT,'(A,8I3)') 'total number of atomic functions per IRED ',
     *(nfunctperIRED(iredrun),iredrun=1,numbofsym) 
      write(OUT,*) stars                                                          
      isum=0
      do iredrun=1,numbofsym
      itotalperIR(iredrun)=nfunctperIRED(iredrun)+iadditional(iredrun)
      isum=isum+itotalperIR(iredrun)
      enddo 
      numballcart=isuM
      iorbrun=0
      do iredrun=1,numbofsym 
      do inired=1,itotalperIR(iredrun)
      iorbrun=iorbrun+1 
      IREDoffunctnew(Iorbrun)=iredrun 
      enddo 
      enddo 
      if (isum.gt.Mxcart) stop 'sorry to many cartesian functions' 
      write(OUT,*) stars                                                          
      write(OUT,'(A,8I3)') 'including additional functions per IRED ',
     *(itotalperIR(iredrun),iredrun=1,numbofsym) 
      write(OUT,*) stars                                                          
      do iredrun=1,numbofsym  
      ibeginIRED(iredrun)=0                              
      enddo
      do lrun=0,Lhigh
      do mrun=-lrun,lrun
      iredLM(mrun,lrun)=ipow2ired(ipowxyz(1,mrun,Lrun),
     *ipowxyz(2,mrun,Lrun),
     *ipowxyz(3,mrun,Lrun))
      incrLM(mrun,lrun)=ibeginIRED(iredLM(mrun,lrun))
      	if (HERMIT.and.(.not.MOLCAS)) then 
      		ibeginIRED(iredLM(mrun,lrUn))=
     *		ibeginIRED(iredLM(mrun,lrun))+1                
      	elseif (HERMIT.and.MOLCAS) then 
      		ibeginIRED(iredLM(mrun,lrUn))=
     *		ibeginIRED(iredLM(mrun,lrun))+ncontrac(lrun) 
      	endif 
      enddo
      if (HERMIT.and.(.not.MOLCAS)) then
      do iredrun=1,numbofsyM
      ibeginIRED(iredrun)=
     *ibeginIRED(iredrun)+(ncontrac(lrun)-1)*
     *nmbMperIRL(iredrun,lrun)   
      enddo   
      endif 
      enddo   
c     do lrun=0,Lhigh
c     write(OUT,'(A,I4,A,21I3)') 'L= ',lrun,
c    *' shifts inside the IRED',
c    *(incrLM(mrun,lrun),mrun=-lrun,lrun)
c     enddo
      shiftIRED(1)=0
      do iredrun=2,numbofsym  
      shiftIRED(iredrun)=shiftIRED(iredrun-1)
     *                   +itotalperIR(iredrun-1)
      enddo
      shiftIREDold(1)=0
      do iredrun=2,numbofsym  
      shiftIREDold(iredrun)=shiftIREDold(iredrun-1)
     *                   +iadditional(iredrun-1)   
      enddo
c     write(OUT,'(A,8I4)') 'shifts for the IREDs ',
c    *(shiftIRED(iredrun),iredrun=1,numbofsym) 
cbs   test all orbital numbers
c     do lrun=0,Lhigh
c     do mrun=-Lrun,Lrun
c     do irun=1,ncontrac(lrun)
c     write(OUT,*) 'L,M,contr funct, absolute number ',
c    *lrun,mrun,irun,shiftired(iredLM(mrun,lrun))+
c    *incrLM(mrun,Lrun)+irun
c     enddo
c     enddo
c     enddo
      shiftIRIR(1)=0
      irun=1
      do ired1=2,numbofsym 
      do ired2=1,ired1 
      irun=irun+1
      if (ired2.eq.1) then
      shiftIRIR(irun)=shiftIRIR(irun-1)+
     *(itotalperIR(ired1-1)*itotalperIR(ired1-1)+
     *itotalperIR(ired1-1))/2
      else
      shiftIRIR(irun)=shiftIRIR(irun-1)+
     *itotalperIR(ired1)*itotalperIR(ired2-1)
      endif 
c     write(OUT,*) 'ired1,ired2 ',ired1,ired2,
c    *irun,shiftIRIR(irun)
      enddo
      enddo
cbs  
      if (HERMIT.and.(.not.MOLCAS)) then  
      do lrun=0,Lhigh
      do Mrun=-Lrun,Lrun
      ired=iredLM(Mrun,Lrun)
      ishifter=shiftIRED(ired)-shiftiredold(ired)+incrLM(mrun,lrun)
      do icart=1,ncontrac(Lrun)
      moffunction(ishifter+1+(icart-1)*nmbMperIRL(ired,lruN))=Mrun
      Loffunction(ishifter+1+(icart-1)*nmbMperIRL(ired,lruN))=Lrun
      IREDoffunction(ishifter+1+(Icart-1)*nmbMperIRL(ired,lruN))=ired  
      numboftype(ishifter+1+(Icart-1)*nmbMperIRL(ired,lruN))=icart 
      enddo
      enddo
      enddo
      elseif (HERMIT.and.MOLCAS) then  
      do lrun=0,Lhigh
      do Mrun=-Lrun,Lrun
      ired=iredLM(Mrun,Lrun)
      ishifter=shiftIRED(ired)-shiftiredold(ired)+incrLM(mrun,lrun)      
      do icart=1,ncontrac(Lrun)
      moffunction(ishifter+icart)=Mrun
      Loffunction(ishifter+icart)=Lrun
      IREDoffunction(ishifter+Icart)=ired  
      numboftype(ishifter+Icart)=icart 
      enddo
      enddo
      enddo
      endif 
      write(OUT,*) stars
      write(OUT,*) 'SYMMETRY-INFORMATION ON FUNCTIONS '
      write(OUT,*) stars
      do irun = 1, numbofcart
      write(OUT,'(4(A,I3))') 'Number of function: ',
     *irun+shiftIREDold(IREDoffunction(irun)),
     *' IR of function: ',IREDoffunction(irun),
     *' L-value: ',Loffunction(irun),
     *' M-value: ',Moffunction(irun)
      numboffunct(irun)=irun+shiftIREDold(IREDoffunction(irun))
      if (IREDoffunction(irun).ne.IREDoffunction(irun+1)) 
     *write(OUT,*)   
      enddo
      do nsymrun=1,numbofsym
      idelpersym(nsymrun)=0
      enddo 
      do nsymrun=1,numbofsym
      nrtofiperIR(nsymrun)=itotalperIR(nsymrun) 
      enddo
      if (AIMP) then 
cbs   generate list of orbitals to be removed
      ikeeporb=0
      numbprev=0
      do irun=1,numbofcart
4712  if (irun.eq.1.or.(irun.ge.2.and.numboffunct(irun).eq.  
     *numbprev+1)) then 
      Lval=Loffunction(irun)
      number=numboffunct(irun)
      itype=numboftype(irun) 
      if (itype.le.icore(lval)) then
      write(6,777) number,itype,lval
      idelpersym(IREDoffunction(irun))=
     *               idelpersym(IREDoffunction(irun))+1
      numbprev=number
      else 
      ikeeporb=ikeeporb+1
      ikeeplist(ikeeporb)=number
      numbprev=number
      endif 
      else
      ikeeporb=ikeeporb+1
      ikeeplist(ikeeporb)=numbprev+1
      numbprev=numbprev+1     
      goto 4712
      endif 
      enddo 
      ikeeporb=0  
      do nsymrun=1,numbofsym
      nrtofiperIR(nsymrun)=itotalperIR(nsymrun)-idelpersym(nsymrun)  
      enddo
      do nsymrun=1,numbofsym
      ikeeporb=ikeeporb+nrtofiperIR(nsymrun)                        
      enddo
      write(OUT,*) stars                                                          
      write(OUT,'(A,8I3)') '# of funct. per IRED after removing core ',
     *(nrtofiperIR(iredrun),iredrun=1,numbofsym) 
      write(6,*) ikeeporb,' orbitals left after deleting core' 
      endif 
      write(OUT,*) stars
      return 
777   format('ORBITAL NUMBER ',I4,' IS THE ',I3,'TH of L-value ',I3,
     *' IT WILL BE REMOVED !!!')  
      end 
