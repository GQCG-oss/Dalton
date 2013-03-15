!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
#include "single.h"
C
      Subroutine amfi(IN,LUPROP,BREIT,FINITE,EXP_FIN,
     *                WRK,LFREE)            
CBS
CBS IN:  Input file, to be replaced by direct reading from DALTON arrays..
CBS LUPROP: Unit for writing the atomic integrals 
CBS BREIT: FLAG to switch to Breit-Pauli (Douglas-Kroll is the default)
CBS FINITE: flag whether to use a finite nucleus or not ...
CB  EXP_FIN: the finite nucleus exponent (if required)
CBS WRK, KFREE,LFREE standard work-array parameters in DALTON
CBS
CBS
#include "implicit.h"
c###########################################################################
c
c          A M F I 
c
c    Atomic Mean-Field Spin-Orbit Integral Program   
c
c Integral-code to generate the one- and two-electron spin-orbit integrals 
c in the no-pair approximation for an atom.
c 
c basis set is built by atomic functions of the form:
c 
c     f(r,Omega)= r**l Y_(lm) (Omega) 
c
c Allthough the code is created with a lot of care and love for 
c the details, the author doesn't give any warranty for it's 
c correctness. 
c
c B.Schimmelpfennig  Fysikum/Stockholm Summer 1996 
c
c If you use this code, please honour the authors work 
c by citing this work properly. 
c
c The author would like to thank the Deutsche Forschungsgemeinschaft 
c for financing this project by a Forschungsstipendium.
c
c
c   The spatial integrals are expected to be used with a spin part 
c   expressed in Paulis spin-matrices rather than with the Spin-operator
c   itself. So if a factor of two is somehow missing, check whether the 
c   same form of the operator is used.
c
c
c   WARNING !!!   WARNING !!   WARNING !!  WARNING !!   WARNING !!
c 
c   when writing spin-same-orbit and spin-other-oribt with sigma_i:
c 
c   For the spin-other-orbit-integrals particle 1 and 2 are exchanged
c   on the arrays carteXOO,carteYOO,carteZOO!!!!!!!!!
c 
c   The reason is to use most of the same-orbit part again and to 
c   have the same symmetry for the integrals on the arrays.
c
c
c   if the spin-other-orbit-part is used in the formulation with 
c   sigma_j, the particles are of cause not interchanged.
c
c
c
c   (i|HSO_mean|j) = (ij) + 1/2 * sum_M  occ(M) {
c                   2(ij|MM)_same - (iM|jM)_same -2(iM|jM)_other
c                   + (jM|iM)_same +2(jM|iM)_other } 
c
c   in the subroutines some signs are changed  to reorder indices 
c   in the integrals to (iM|jM) or (Mi|Mj) accoding to the way they 
c   were calculated before. 
c
c
c
c   one-particle integrals (really one-particle or mean-field) 
c   are written to files in CONTANDMULT. Look there for information on 
c   the format of files.  
c
c
c  BUGS:  There is still a strange sign-error in the two-electron-integrals 
c  if one applies straight-forward the formulae of the documentation.
c  This problem has been solved by the the cheater...
c
c  Everybody is welcome to find the problem in the formulas ........
c
c  First reasonable results on Thallium (SD with frozen 5D) 14.10.96
c
c
c
c
c
c  Connection to MOLCAS: 
c  How wonderful, they normalize the functions exactly as I do, which 
c  means they use the correct linear combinations.  
c
c  Exponents and coefficients are expected in the MOLCAS-Format      
c  first exponents 
c  coefficients afterwards   
c
c                                           8.5.97   
c                                                           
c  New version for DALTON canibalized from the MOLCAS version september 2000
c                                                           
c###########################################################################
#include "para.h"
      logical keep    ! parameter to decide about keeping angular 
cbs                     ! integrals in memory 
      logical keepcart    ! parameter to decide about keeping cartesian
cbs                         ! integrals in memory 
      logical makemean   ! parameter to decide about generating a meanfield   
      logical bonn       ! if bonn is set, Bonn-approach for spin-other orbit
      logical breit      ! if breit is set, BREIT-PAULI only                 
      logical SAMEORB    ! parameter for same-orbit only 
      logical AIMP       ! parameter to delete CORE for AIMP     
      logical oneonly    ! parameter to use only oneelectron integrals    
      logical FINITE
      character*4  symmetry  
#include "datapow.h"
      common ipowxyz(3,-Lmax:Lmax,0:Lmax)
      dimension WRK(LFREE)
c##########################################################################
cbs  #####################################################################
cbs         version with all angular integrals in memory 
c         keep=.true.
cbs  #####################################################################
cbs         version without  all angular integrals in memory 
          keep=.false. 
cbs  #####################################################################
cbs         version without  all cartesian integrals in memory 
          keepcart=.false. 
cbs  #####################################################################
cbs         version with all cartesian integrals in memory 
c         keepcart=.true.
cbs  #####################################################################
cbs   initialize tables with double facultatives... 
      call inidf
cbs   move some powers of x,y,z to the right place   BEGIN 
cbs   check if Lpowmax is high enough..
      if (Lpowmax.lt.Lmax) then 
      CALL QUIT('increase lpowmax and edit ixyzpow')
      endif 
      jrun=1
      do irun=0,Lmax
      do Mval=-irun,irun
      ipowxyz(1,Mval,irun)=ixyzpow(jrun)
      ipowxyz(2,Mval,irun)=ixyzpow(jrun+1)
      ipowxyz(3,Mval,irun)=ixyzpow(jrun+2)
      jrun=jrun+3
      enddo
      enddo
cbs   move some powers of x,y,z to the right place   END   
      if (FINITE) then 
        ifinite=1
      else
        ifinite=0
      endif 
cbs   read the input 
      call readbas(Lhigh,makemean,bonn,breit,
     *symmetry,sameorb,AIMP,oneonly,ncont4,numballcart,IN,
     *ifinite,EXP_FIN) 
cbs
cbs
 123  if (ifinite.eq.2) call finsub  
cbs
cbs
! Lhigh is the highest l-value in the basis set
      if (makemean.and.(.not.oneonly).and.ifinite.le.1) 
     *call getAOs(Lhigh)
      call genpowers(Lhigh) !generate powers of exponents and overlaps 
cbs   start generating modified contraction coefficients
cbs   generate starting adresses of contraction coefficients  on 
cbs   contrarray 
      call genstar(Lhigh)  
cbs   generate ovlp of normalized primitives 
      call genovlp(Lhigh)    
      do lrun=0,Lhigh
cbs   cont(L) arranges all the contraction coefficients for a given L-value
cbs   and renormalizes them 
      call cont(lrun,breit,ifinite)
      enddo 
cbs                        
cbs        beginning the angular part  
      if (.not.oneonly) then  
CBS   write(6,*) '***************************************************' 
CBS   write(6,*) '********   beginning the 2e-part ******************' 
CBS   write(6,*) '***************************************************' 
cbs                        
cbs  ##################################################################################### 
cbs  ##################################################################################### 
cbs  ##################################################################################### 
cbs
cbs       
      call angular(Lhigh,keep,keepcart,makemean,bonn,breit,
     *sameorb,ifinite,WRK,LFREE) ! subroutine for angular part 
      endif 
      if (ifinite.eq.1) then ! redo everything for finite core
CBS   write(6,*) 'once more the two-electron integrals'
      ifinite=2
      goto 123
      endif 
cbs ########################################################################################
cbs ########################################################################################
cbs ########################################################################################
CBS   write(6,*) '***************************************************' 
CBS   write(6,*) '*******   beginning the 1-electron-part  **********' 
CBS   write(6,*) '***************************************************' 
cbs    the one-electron spin-orbit integrals   
      call gen1overR3(Lhigh)   ! generates the 1/r**3 integrals  for normalized functions 
      call contandmult(Lhigh,makemean,AIMP,oneonly,numballcart,LUPROP,
     *ifinite,WRK,LFREE) ! multiplies radial integrals with l,m-dependent
cbs                             factors and contraction coefficients 
CBS   write(6,*) '***************************************************' 
CBS   write(6,*) '*******   end of  the 1-electron-part    **********' 
CBS   write(6,*) '***************************************************' 
cbs ########################################################################################
cbs ########################################################################################
cbs ########################################################################################
      Return
      end   
      subroutine finsub
cbs
cbs   subroutine to set up parameters for finite nucleus. The s-functions are replaced 
cbs   by just one exponent which models the nucleus.
cbs
#include "implicit.h"
#include "para.h"
#include "amfi_param.h"
      common /nucleus/ charge,Exp_finite    
      noccorb(0)=1  
      do l=1,lmax_occ
      noccorb(l)=0              
      enddo 
      occup(1,0)=-charge 
      nprimit_keep=nprimit(0)
      ncontrac_keep=ncontrac(0)
      nprimit(0)=1
      ncontrac(0)=1
      exponents(1,0)=0.5d0*Exp_finite       
      return 
      end  
      
      
      subroutine angular(Lhigh,keep,keepcart,makemean,bonn,
     *breit,sameorb,ifinite,WRK,LFREE)
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
#include "implicit.h"
#include "priunit.h"
#include "para.h"
#include "amfi_param.h"
      logical keep,keepcart,icheck,mcheckxy,mcheckz,makemean,bonn,   
     *breiT,sameorb,cleaner,NFINI
cbs   NFINI means not finite nucleus  
      dimension l2block(0:Lmax,0:Lmax,0:Lmax,0:Lmax)
      dimension WRK(LFREE)
cbs #####################################################################
cbs   some preparation of factors needed later on..                     #
cbs ######################################################################
      ipnt(i,j)=(max(i,j)*max(i,j)-max(i,j))/2+min(i,j) 
      roottwo=dsqrt(2d0) 
cbs   calculate some prefactors that will be needed quite often      
      call prefac(Lmax,preroots,clebsch) 
        if (ifinite.ne.2) then 
cbs     clean array for one electron integrals
        iprod=MxcontL*MxcontL*(Lmax+Lmax+1)*(Lmax+1)*Lmax     
        call dzero(onecartX,iprod) 
        call dzero(onecartY,iprod) 
        call dzero(onecartZ,iprod) 
        NFINI=.true.   
        else
        NFINI=.false.   
        endif 
cbs   generate an array with sign for (even/odd) m-values
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
      endif 
      enddo
      enddo
      enddo
      enddo
cbs #####################################################################
cbs   isgnprod gives the sign due to powers (-1)**M  this are again 
cbs   angular m-values 
cbs #####################################################################
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
cbs   same orbit integrals integrals  on carteXSO carteYSO and carteSO                            
cbs   other orbit integrals  on carteXOO carteYOO and carteOO                            
      iangfirst=0 ! first block of angular integrals                   
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
cbs   loop over all possible < l1 l2, l3 l4 > blocks
CBS   write(6,'(A)') '   L1   L2   L3   L4' 
      do l1=0,Lhigh   ! improving is probably possible...
      do l2=0,Lhigh
      do l3=0,l1
      do l4=0,l2
cbs   check parity
      if (mod(l1+l2+l3+l4,2).eq.0) then
cbs   check that Lleft and Lright do not always differ by more than one
cbs   a difference of two means two spin flips and is therefore not allowed
      Lleftmax=l1+l2
      Lrightmax=l3+l4
      Lleftmin=iabs(l1-l2)
      Lrightmin=iabs(l3-l4)
      if ((Lrightmin-Lleftmax.le.1.and.Lrightmax-Lleftmin.gt.-1).or.
     *(Lleftmin-Lrightmax.le.1.and.Lleftmax-Lrightmin.gt.-1)) then
cbs   additional check for mean-field
      if ((l1.eq.l3.and.l2.eq.l4).or.(l1.eq.l2.and.l3.eq.l4)) then
      if (l1+l3.ne.0) then
CBS   write(6,'(4I5)') l1,l2,l3,l4
CBS   now I determine the size of the angular integral arrays
        jblock=0
        do m1=-l1,l1
        do m2=-l2,l2
        do m3=-l3,l3           
        m4=m1+m2-m3+1
        if (iabs(m4).le.l4) then
        if ((.not.makemean).or.
     *  (l1.eq.l3.and.l2.eq.l4.and.iabs(m2).eq.iabs(m4)).or.
     *  (l1.eq.l2.and.l3.eq.l4.and.
     *  (iabs(m1).eq.iabs(m2).or.iabs(m3).eq.iabs(m4)))) then   
        jblock=jblock+1
        endif 
        endif 
        enddo 
        enddo 
        enddo 
        do m1=  0,l1
        do m2=-l2,l2
        do m3=-l3,l3 
        m4=m1+m2-m3
        if ((.not.makemean).or.
     *  (l1.eq.l3.and.l2.eq.l4.and.iabs(m2).eq.iabs(m4)).or.
     *  (l1.eq.l2.and.l3.eq.l4.and.
     *  (iabs(m1).eq.iabs(m2).or.iabs(m3).eq.iabs(m4)))) then         
        if (m1.ne.0.or.m2.ne.0.or.m3.ne.0) then !  all m eqal 0 make no sense....
        if (iabs(m4).le.l4)  then                                               
        jblock=jblock+1
        endif 
        endif 
        endif 
        enddo 
        enddo 
        enddo 
CBS   done !!                                                     
cbs     number of contracted integrals for each block 
        ncont=ncontrac(l1)*ncontrac(l2)*
     *  ncontrac(l3)*ncontrac(l4) 
      mxangint=jblock*ncont
cbs   determine the size icont4 for the radial integrals 
      call gencoulDIM(l1,l2,l3,l4,makemean,bonn,breit,
     *sameorb,icont4) 
      IANGSO = 1
      iangOO=iangSO+mxangint
      icartSO=iangOO+mxangint   
      icartOO=icartSO+ncont
      iconSO=icartOO+ncont
      iconOO=iconSO+icont4   
      KLAST = ICONOO + ICONT4
      IF (KLAST .GT. LFREE) CALL STOPIT('AMFI  ','ANGULAR',KLAST,LFREE)
      LLEFT = LFREE - KLAST + 1
      call gencoul(l1,l2,l3,l4,makemean,bonn,breit,
     *sameorb,WRK(iconSO),WRK(iconOO),icont4,
     *WRK(KLAST),LLEFT) ! generates and transforms integrals
        l2block(l1,l2,l3,l4)=1  ! can be used for getting the
cbs   local counter for integral adresses  
        mblock=0 ! counter of (m,m,m,m)-blocks for (l1,l2,l3,l4)   
cbs     if keep is set to false, the angular integrals are 
cbs     thrown away after each block of l-values 
cbs     which means integrals start at address 0 
        if (.not.keep) iangfirst=0
        locstar=iangfirst ! local starting adress counter 
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
        mcombina(1,m1,m2,m3,m4)=1    
        mblock=mblock+1
        if (locstar+ncont.gt.mxangint) then 
        write(LUPRI,*)'not enough space allocated for angular integrals'
        write(LUPRI,*) 'increase mxangint to at least ',
     *  locstar+ncont        
        CALL QUIT('Out of dimensional bounds in AMFI')
        endif  
cbs mkangLmin = make_angular_integrals_for_L- type operator 
cbs really generates  the angular prefactors and combines them with 
cbs the radial integrals
        call mkangLmin(Lmax,l1,l2,l3,l4,m1,m2,m3,m4,
     *       WRK(iangSO+locstar),
     *       WRK(iangOO+locstar),
     *       Lfirst(1),Llast(1),Lblocks(1),
     *       ncontrac(l1),ncontrac(l2),ncontrac(l3),ncontrac(l4),
     *       WRK(iconSO+Lstarter(1)-1),
     *       WRK(iconSO+Lstarter(2)-1),
     *       WRK(iconSO+Lstarter(3)-1),
     *       WRK(iconSO+Lstarter(4)-1),
     *       WRK(iconOO+Lstarter(1)-1),
     *       WRK(iconOO+Lstarter(2)-1),
     *       WRK(iconOO+Lstarter(3)-1),
     *       WRK(iconOO+Lstarter(4)-1),
     *       preroots,clebsch,scratch4,bonn,breit,
     *       sameorb) 
        locstar=locstar+ncont ! increase starting address 
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
        mcombina(1,m1,m2,m3,m4)=2
        mblock=mblock+1
        if (locstar+ncont.gt.mxangint) then 
        write(LUPRI,*)'not enough space allocated for angular integrals'
        write(LUPRI,*) 'increase mxangint to at least ',
     *  locstar+ncont        
        CALL QUIT('Out of dimensional bounds in AMFI')
        endif  
        call mkangL0(Lmax,l1,l2,l3,l4,m1,m2,m3,m4,
     *       WRK(iangSO+locstar),
     *       WRK(iangOO+locstar),
     *       Lfirst(1),Llast(1),Lblocks(1),
     *       ncontrac(l1),ncontrac(l2),ncontrac(l3),ncontrac(l4),
     *       WRK(iconSO+Lstarter(1)-1),
     *       WRK(iconSO+Lstarter(2)-1),
     *       WRK(iconSO+Lstarter(3)-1),
     *       WRK(iconSO+Lstarter(4)-1),
     *       WRK(iconOO+Lstarter(1)-1),
     *       WRK(iconOO+Lstarter(2)-1),
     *       WRK(iconOO+Lstarter(3)-1),
     *       WRK(iconOO+Lstarter(4)-1),
     *       preroots,clebsch,scratch4,bonn,breit,
     *       sameorb) 
        locstar=locstar+ncont
        mcombina(2,m1,m2,m3,m4)=mblock  
        endif 
        endif 
        endif 
        enddo
        enddo
        enddo
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
        mcombcart(1,m1,m2,m3,m4)=0     ! will hold the type  x=1 y=2 z=3 
        mcombcart(2,m1,m2,m3,m4)=0     ! will hold the block number
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
        mcombcart(1,m1,m2,m3,m4)=1
        mcombcart(2,m1,m2,m3,m4)=mblockx                               
        call tosigX(m1,m2,m3,m4,WRK(iangSO+iangfirst),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),WRK(icartSO),preXZ,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner) 
c
        if (.not.bonn.and.(.not.breiT)) 
     *  call tosigX(m1,m2,m3,m4,WRK(iangOO+iangfirst),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),WRK(icartOO),preXZ,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner)
        if (makemean) then ! generate mean-field-contributions
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
             if (l1.eq.l3.and.l2.eq.l4) then 
             if (m2.eq.m4.and.m1.lt.m3.and.
     *       iabs(m1+m3).eq.1.and.l1.ne.0) then 
             call two2mean13(WRK(icartSO),occup(1,l2),
     *       AOcoeffs(1,1,l2),onecartx(1,1,ipnt(m1+l1+1,m3+l3+1),l1),
     *       ncontrac(l1),ncontrac(l2),noccorb(l2))         
             endif 
             endif 
             if (l1.eq.l2.and.l3.eq.l4) then 
             if (m1.eq.m2.and.l3.ne.0.and.l3.ne.l1) then 
             if (m3.lt.m4.and.iabs(m4+m3).eq.1) then 
cbs   for the "Bonn-approach"   exchange cartexOO by cartexSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean34a(WRK(icartSO),
     *       WRK(icartSO),
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartx(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             else 
             if(NFINI) call two2mean34a(WRK(icartSO),
     *       WRK(icartOO),
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartx(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             endif 
             endif 
             if (m3.gt.m4.and.iabs(m4+m3).eq.1) then 
cbs   for the "Bonn-approach"   exchange cartexOO by cartexSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean34b(WRK(icartSO),
     *       WRK(icartSO),       
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartx(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             else 
             if (NFINI) call two2mean34b(WRK(icartSO),    
     *       WRK(icartOO),    
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartx(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             endif 
             endif 
             elseif(m3.eq.m4.and.l1.ne.0) then 
             if (m1.lt.m2.and.iabs(m1+m2).eq.1) then 
cbs   for the "Bonn-approach"   exchange cartexOO by cartexSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean12a(WRK(icartSO),   
     *       WRK(icartSO),occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartx(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             else 
             if (NFINI) call two2mean12a(WRK(icartSO),     
     *       WRK(icartOO),occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartx(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             endif 
             endif 
             if (m1.gt.m2.and.iabs(m1+m2).eq.1) then 
cbs   for the "Bonn-approach"   exchange cartexOO by cartexSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean12b(WRK(icartSO),    
     *       WRK(icartSO),occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartx(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             else 
             if (NFINI) call two2mean12b(WRK(icartSO),    
     *       WRK(icartOO),occup(1,l3),
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++      SIGMA Y      ++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        elseif (indx.eq.1.and.indy.eq.0.and.indz.eq.1.and.  
     *  icheckxy(m1,m2,m3,m4)) then  ! X*Z transforms like L_y  (B2) 
cbs     integrals for sigma_y 
        mblocky=mblocky+1
        mcombcart(1,m1,m2,m3,m4)=2
        mcombcart(2,m1,m2,m3,m4)=mblocky                               
        call tosigY(m1,m2,m3,m4,WRK(iangSO+iangfirst),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),WRK(icartSO),preY,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner) 
c
        if (.not.bonn.and.(.not.breit)) 
     *  call tosigY(m1,m2,m3,m4,WRK(iangOO+iangfirst),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),WRK(icartOO),preY,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner) 
        if (makemean) then ! generate mean-field-contributions
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
             if (l1.eq.l3.and.l2.eq.l4) then 
             if (m2.eq.m4.and.m1.lt.m3.
     *       and.iabs(m3-m1).eq.1.and.l1.ne.0) then 
             call two2mean13(WRK(icartSO),occup(1,l2),
     *       AOcoeffs(1,1,l2),onecartY(1,1,ipnt(m1+l1+1,m3+l3+1),l1),
     *       ncontrac(l1),ncontrac(l2),noccorb(l2))         
             endif 
             endif 
             if (l1.eq.l2.and.l3.eq.l4) then 
             if (m1.eq.m2.and.l3.ne.0.and.l3.ne.l1) then 
             if (m3.lt.m4.and.iabs(m3-m4).eq.1) then 
cbs   for the "Bonn-approach"   exchange carteYOO by carteYSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean34a(WRK(icartSO),
     *       WRK(icartSO),occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartY(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             else 
             if (NFINI) call two2mean34a(WRK(icartSO),
     *       WRK(icartOO),occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartY(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             endif 
             endif 
             if (m3.gt.m4.and.iabs(m3-m4).eq.1) then 
cbs   for the "Bonn-approach"   exchange carteYOO by carteYSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean34b(WRK(icartSO),
     *       WRK(icartSO),occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartY(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             else 
             if (NFINI) call two2mean34b(WRK(icartSO),
     *       WRK(icartOO),occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartY(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             endif 
             endif 
             elseif(m3.eq.m4.and.l1.ne.0) then 
             if (m1.lt.m2.and.iabs(m1-m2).eq.1) then 
cbs   for the "Bonn-approach"   exchange carteOO by carteSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean12a(WRK(icartSO),
     *       WRK(icartSO),occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartY(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             else 
             if (NFINI) call two2mean12a(WRK(icartSO),
     *       WRK(icartOO),occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartY(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             endif 
             endif 
             if (m1.gt.m2.anD.Iabs(m1-m2).eq.1) then 
cbs   for the "Bonn-approach"   exchange carteYOO by carteYSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean12b(WRK(icartSO),
     *       WRK(icartSO),occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartY(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             else 
             if (NFINI) call two2mean12b(WRK(icartSO),
     *       WRK(icartOO),occup(1,l3),
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++      SIGMA Z      ++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        elseif (indx.eq.1.and.indy.eq.1.and.indz.eq.0.and.  
     *  icheckz(m1,m2,m3,m4)) then ! X*Y transforms like L_z  (A2) 
cbs     integrals for sigma_z 
        mblockz=mblockz+1
        mcombcart(1,m1,m2,m3,m4)=3
        mcombcart(2,m1,m2,m3,m4)=mblockz                               
        call tosigZ(m1,m2,m3,m4,WRK(iangSO+iangfirst),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),WRK(icartSO),preXZ,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner) 
c
        if (.not.bonn.and.(.not.breit)) 
     *  call tosigZ(m1,m2,m3,m4,WRK(iangOO+iangfirst),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),WRK(icartOO),preXZ,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner) 
        if (makemean) then ! generate mean-field-contributions
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
             if (l1.eq.l3.and.l2.eq.l4) then 
             if (m2.eq.m4.and.m1.lt.m3.
     *       and.m1.eq.-m3.and.l1.ne.0) then 
             call two2mean13(WRK(icartSO),occup(1,l2),
     *       AOcoeffs(1,1,l2),onecartz(1,1,ipnt(m1+l1+1,m3+l3+1),l1),
     *       ncontrac(l1),ncontrac(l2),noccorb(l2))         
             endif 
             endif 
             if (l1.eq.l2.and.l3.eq.l4) then 
             if (m1.eq.m2.and.l3.ne.0.and.l3.ne.l1) then 
             if (m3.lt.m4.and.m3.eq.-m4) then 
cbs   for the "Bonn-approach"   exchange carteOO by carteSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean34a(WRK(icartSO),
     *       WRK(icartSO),occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartz(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             else 
             if (NFINI) call two2mean34a(WRK(icartSO),
     *       WRK(icartOO),occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartz(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             endif 
             endif 
             if (m3.gt.m4.and.m3.eq.-m4) then 
cbs   for the "Bonn-approach"   exchange carteOO by carteSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean34b(WRK(icartSO),
     *       WRK(icartSO),occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartz(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             else 
             if (NFINI) call two2mean34b(WRK(icartSO),
     *       WRK(icartOO),
     *       occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartz(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)         
             endif 
             endif 
             elseif(m3.eq.m4.and.l1.ne.0) then 
             if (m1.lt.m2.and.m1.eq.-m2) then 
cbs   for the "Bonn-approach"   exchange carteOO by carteSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean12a(WRK(icartSO),
     *       WRK(icartSO),occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartz(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             else 
             if (NFINI) call two2mean12a(WRK(icartSO),
     *       WRK(icartOO),
     *       occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartz(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             endif 
             endif 
             if (m1.gt.m2.and.m1.eq.-m2) then 
cbs   for the "Bonn-approach"   exchange carteOO by carteSO 
             if (bonn.or.breiT) then 
             if (NFINI) call two2mean12b(WRK(icartSO),
     *       WRK(icartSO),
     *       occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartz(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)         
             else 
             if (NFINI) call two2mean12b(WRK(icartSO),
     *       WRK(icartOO),
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
        endif    
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        endif ! for check of significance for meanfield.  
        enddo 
        enddo 
        enddo 
        enddo 
        numbcart=numbcart+(mblockx+mblocky+mblockz)*ncont
cbs   just controlling if x and y integrals have the same number of blocks 
      if (mblockx.ne.mblocky) then 
      write(LUPRI,*) 
     *'numbers of integrals for sigma_x and sigma_y not equal!'    
      write(LUPRI,'(A12,4I3,2(A3,I5))') 
     *'l1,l2,l3,l4 ',l1,l2,l3,l4,' X:',mblockx,' Y:',mblocky  
      write(LUPRI,*) ' check the ipowxyz-array'
      CALL QUIT('Problems with IPOWXYA array in AMFI')
      endif   
cbs   start adresses for the next <ll|ll> block of integrals 
      endif
      endif
      endif
      endif
      enddo
      enddo
      enddo
      enddo
      return 
      end  
      subroutine buildcoul(l1,l2,l3,l4,! angular momenta of primitives
     *incl1,incl3, ! shifts for different radial integrals
     *Lrun, ! L-value for coulomb integrals 
     *prmints,
     *nprim1,nprim2,nprim3,nprim4,  ! number of primitives
     *expo1,expo2,expo3,expo4, ! arrays with the exponents
     *power13,
     *power24,
     *quotpow1,quotpow2
     *)
cbs ##################################################################
c
cbs  purpose: builds up the coulomb integrals 
cbs  inbetween primitives and multiplies 
cbs  with extra factors to correct the 
cbs  normalization                  
c
cbs ##################################################################
#include "implicit.h"
#include "para.h"
#include "amfi_param.h"
#include "dofuc.h"
#include "pi.h"
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
      root8ovpi=dsqrt(8d0/pi)    
cbs ##################################################################
cbs        prepare indices for coulint
cbs ##################################################################
      n1=l1+incl1+1
      n2=l2+1
      n3=l3+incl3+1
      n4=l4+1
      n13=n1+n3
      n24=n2+n4
      index1=N13-Lrun-1
      index2=n24+Lrun
      index3=N24-Lrun-1
      index4=n13+Lrun
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
                   Cfunctx1(irun1)=fraclist1(0)
                   frac(irun1)=a1324/(1d0+a1324)
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
      return 
      end 
      subroutine cartoneX(L,Lmax,onecontr,ncontrac,
     *MxcontL,onecartX)
#include "implicit.h"
      dimension onecontr(MxcontL,MxcontL,-Lmax:Lmax,3),
     *onecartX(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1))  
cbs   arranges the cartesian one-elctron-integrals for X  on a 
cbs   quadratic matrix 
      ipnt(I,J)=(max(i,j)*(max(i,j)-1))/2+min(i,j)
cbs   - + Integrals    m || mprime     mprime=m+1
      do Mprime=2,L
      M=mprime-1 
      iaddr=ipnt(Mprime+L+1,-M+L+1)
      do jcont=1,ncontrac
      do icont=1,ncontrac
      onecartX(icont,jcont,iaddr)=
     *onecartX(icont,jcont,iaddr)
     *-0.25d0*(
     *onecontr(icont,jcont,Mprime,1)+
     *onecontr(icont,jcont,-Mprime,3))
      enddo
      enddo 
      enddo 
cbs   - + Integrals    m || mprime     mprime=m-1
      do Mprime=1,L-1
      M=mprime+1 
      iaddr=ipnt(Mprime+L+1,-M+L+1)
      do jcont=1,ncontrac
      do icont=1,ncontrac
      onecartX(icont,jcont,iaddr)=
     *onecartX(icont,jcont,iaddr)
     *-0.25d0*(
     *onecontr(icont,jcont,Mprime,3)+
     *onecontr(icont,jcont,-Mprime,1))
      enddo
      enddo 
      enddo 
cbs   -1 || 0 integrals 
      pre=dsqrt(0.125d0) 
      iaddr=ipnt(L,L+1)
      do jcont=1,ncontrac
      do icont=1,ncontrac
      onecartX(icont,jcont,iaddr)=
     *onecartX(icont,jcont,iaddr)
     *-pre* (onecontr(icont,jcont,0,3)+
     *onecontr(icont,jcont,0,1) )
      enddo
      enddo 
      return
      end 
      subroutine cartoneY(L,Lmax,onecontr,ncontrac,
     *MxcontL,onecartY)
#include "implicit.h"
      dimension onecontr(MxcontL,MxcontL,-Lmax:Lmax,3),
     *onecartY(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1))  
cbs   arranges the cartesian one-electron integrals for Y  
cbs   on a quadratic matrix 
      ipnt(I,J)=(max(i,j)*(max(i,j)-1))/2+min(i,j)
cbs   + + Integrals    m || mprime     mprime=m+1
      do Mprime=2,L
      M=mprime-1 
      iaddr=ipnt(Mprime+L+1,M+L+1)
      do jcont=1,ncontrac
      do icont=1,ncontrac
      onecartY(icont,jcont,iaddr)=
     *onecartY(icont,jcont,iaddr) 
     *-0.25d0*(
     *onecontr(icont,jcont,Mprime,1)+
     *onecontr(icont,jcont,-Mprime,3))
      enddo
      enddo 
      enddo 
cbs   - - Integrals    m || mprime     mprime=m-1
      do Mprime=1,L-1
      M=mprime+1 
      iaddr=ipnt(-Mprime+L+1,-M+L+1)
      do jcont=1,ncontrac
      do icont=1,ncontrac
      onecartY(icont,jcont,iaddr)=
     *onecartY(icont,jcont,iaddr) 
     *+0.25d0*(
     *onecontr(icont,jcont,Mprime,3)+
     *onecontr(icont,jcont,-Mprime,1))
      enddo
      enddo 
      enddo 
cbs    0 || 1 integrals 
      pre=-dsqrt(0.125d0)      
      iaddr=ipnt(L+1,L+2)
      do jcont=1,ncontrac
      do icont=1,ncontrac
      onecartY(icont,jcont,iaddr)=
     *onecartY(icont,jcont,iaddr) 
     *+pre*    
     *(onecontr(icont,jcont,1,1)+
     *onecontr(icont,jcont,-1,3)) 
      enddo 
      enddo 
      return
      end 
      subroutine cartoneZ(L,Lmax,onecontr,ncontrac,
     *MxcontL,onecartZ)
#include "implicit.h"
      dimension onecontr(MxcontL,MxcontL,-Lmax:Lmax,3),
     *onecartZ(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1))  
cbs   arranges the cartesian one-electron integrals for Z  
cbs   on a quadratic matrix 
      ipnt(I,J)=(max(i,j)*(max(i,j)-1))/2+min(i,j)
cbs   - + Integrals    m || mprime     mprime=m
      do Mprime=1,L
      iaddr=ipnt(Mprime+L+1,-mprime+L+1)
      do jcont=1,ncontrac
      do icont=1,ncontrac
      onecartZ(icont,jcont,iaddr)=
     *onecartZ(icont,jcont,iaddr)+
     *0.5d0*(
     *onecontr(icont,jcont,Mprime,2)-
     *onecontr(icont,jcont,-Mprime,2))
      enddo
      enddo 
      enddo 
      return
      end 
      subroutine chngcont(coeffs,coeffst1,coeffst1a,coeffst2,
     *coeffst2a,ncont,nprims,evec,
     *type1,type2,work,work2,work3,MxprimL,
     *rootOVLP,OVLPinv,exponents)              
c###############################################################################
cbs   purpose: makes out of old contraction coefficients(in normalized functions)
cbs   new coefficients including the kinematical factors
cbs   using the diagonal matrices on type1 and type2 (see subroutine kinemat)
cbs   coeffst1a and coeffst2a additionally include the exponents alpha 
cbs   (that is why ....a). So the exponents in the integrals are moved 
cbs   to the contraction coefficients and not in some way into the primitive 
cbs   integrals. 
cbs
cbs   the different cases for contracted integrals differ later on in the
cbs   choice of different sets of contraction coefficients. 
cbs
c###############################################################################
#include "implicit.h"
      dimension coeffs(nprims,ncont),    ! original contraction coefficients   
     *coeffst1(nprims,ncont),            ! A * contraction coefficients
     *coeffst1a(nprims,ncont),           ! A * alpha*contraction coefficients
     *coeffst2a(nprims,ncont),           ! c*A/(E+m) * contraction coefficients
     *coeffst2(nprims,ncont),            ! c*A/(E+m) * alpha *contraction coefficients    
     *evec(nprims,nprims),
     *work(nprims,nprims) ,
     *work2(nprims,nprims) ,
     *work3(nprims,nprims) ,
     *rootOVLP(MxprimL,*),
     *OVLPinv(MxprimL,*),
     *type1(*),type2(*),
     *exponents(*) 
cbs   
cbs   first new coefficients for type1 (A) 
cbs   generate a transformation matrix on work
cbs   
      do J=1,nprims
      do I=1,nprims
      work(I,J)=0d0
      work2(I,J)=0d0
      work3(I,J)=0d0
      enddo
      enddo
cbs   build up the transformation matrix 
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work(I,J)=work(I,J)+evec(I,K)*type1(K)*evec(J,K)
      enddo
      enddo
      enddo
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work2(I,J)=work2(I,J)+work(I,K)*rootOVLP(K,J)    
      enddo
      enddo
      enddo
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work3(I,J)=work3(I,J)+rootOVLP(I,K)*work2(K,J)    
      enddo
      enddo
      enddo
      do J=1,nprims
      do I=1,nprims
      work(I,J)=0d0
      enddo
      enddo
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work(J,I)=work(J,I)+OVLPinv(I,K)*work3(K,J)       
      enddo
      enddo
      enddo
      do K=1,ncont
      do I=1,nprims
      coeffst1(I,K)=0d0
      enddo
      enddo
cbs   now transform the vectors 
      do K=1,ncont
      do J=1,nprims
      do I=1,nprims
      coeffst1(I,K)=coeffst1(I,K)+work(J,I)*coeffs(J,K)
      enddo 
      enddo 
      enddo 
cbs  
cbs   now with exponent   
cbs    
      do K=1,ncont
      do I=1,nprims
      coeffst1a(I,K)=exponents(I)*coeffst1(I,K) 
      enddo
      enddo
cbs   
cbs   and now the same for the other type  A/(E+m) 
cbs   
      do J=1,nprims
      do I=1,nprims
      work(I,J)=0d0
      work2(I,J)=0d0
      work3(I,J)=0d0
      enddo
      enddo
cbs   build up the transformation matrix 
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work(I,J)=work(I,J)+evec(I,K)*type2(K)*evec(J,K)
      enddo
      enddo
      enddo
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work2(I,J)=work2(I,J)+work(I,K)*rootOVLP(K,J)    
      enddo
      enddo
      enddo
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work3(I,J)=work3(I,J)+rootOVLP(I,K)*work2(K,J)    
      enddo
      enddo
      enddo
      do J=1,nprims
      do I=1,nprims
      work(I,J)=0d0
      enddo
      enddo
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work(J,I)=work(J,I)+OVLPinv(I,K)*work3(K,J)       
      enddo
      enddo
      enddo
      do K=1,ncont
      do I=1,nprims
      coeffst2(I,K)=0d0
      enddo
      enddo
cbs   now transform the vectors 
      do K=1,ncont
      do J=1,nprims
      do I=1,nprims
      coeffst2(I,K)=coeffst2(I,K)+work(J,I)*coeffs(J,K)
      enddo 
      enddo 
      enddo 
cbs  
cbs   now with exponent   
cbs    
      do K=1,ncont
      do I=1,nprims
      coeffst2a(I,K)=exponents(I)*coeffst2(I,K) 
      enddo
      enddo
      return 
      end 
      subroutine cont(L,breit,ifinite)
cbs###########################################################################
cbs   cont prepares all required contraction coefficients for functions 
cbs   with angular momentum L
cbs###########################################################################
#include "implicit.h"
#include "para.h"
#include "amfi_param.h"
      dimension tkintria((MxprimL*MxprimL+MxprimL)/2)
      logical breit,breit_finite
      breit_finite=.true.
cbs   transcon transfers and normalizes contracted functions
cbs   ore more precizely the coefficients  
      call transcon(cntscrtch(1,1,L),MxprimL,
     *MxcontL,normovlp(1,1,L),
     *contrarray(iaddori(L)),nprimit(L),ncontrac(L))
cbs   gentkin generates the matrix of kinetic energy  TKIN 
      call gentkin(L,TKIN,nprimit(L),exponents(1,L),rootOVLPinv(1,1,L))
cbs   kindiag diagonalizes TKIN 
cbs   for finite nucleus 
      if (ifinite.eq.2.and.L.eq.0) then  
      call kindiag(TKIN,TKINTRIA,nprimit(L),evec,eval,breit_finite)  
      else 
      call kindiag(TKIN,TKINTRIA,nprimit(L),evec,eval,breit)
      endif 
cbs   kinemat generates kinematic factors in 
cbs   the basis of eigenvectors   
      call kinemat(L,nprimit(L),eval,type1,type2,Energy)
      incr=nprimit(L)*ncontrac(L)
cbs   chngcont= changecont generates the contraction coeffs
cbs   including kinematic factors and even exponents as factors
      call chngcont(
     *contrarray(iaddori(L)),                  
     *contrarray(iaddtyp1(L)),                  
     *contrarray(iaddtyp2(L)),                  
     *contrarray(iaddtyp3(L)),                  
     *contrarray(iaddtyp4(L)),                  
     *ncontrac(L),nprimit(L),evec,
     *type1,type2,scratch4,scratch4(nprimit(L)*nprimit(L)+1),
     *scratch4(2*nprimit(L)*nprimit(L)+1),MxprimL,
     *rootOVLP(1,1,L),OVLPinv(1,1,L),
     *exponents(1,L))
      return
      end
      Subroutine contandmult(Lhigh,makemean,AIMP,oneonly,numballcart,
     *LUPROP,ifinite,WRK,LWRK)
#include "implicit.h"
#include "para.h"
#include "amfi_param.h"
#include "ired.h"
      logical makemean,AIMP,oneonly  
      character*8 xa,ya,za
      dimension xa(4),ya(4),za(4)
      DIMENSION WRK(LWRK)
      common /nucleus/ charge,Exp_Finite
      double precision normasHERMIT(-Lmax:Lmax,0:Lmax)
      data ((normasHERMIT(ml,l),ml=-lmax,lmax),l=0,lmax) 
     &     /0.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     &     0.0d0,0.0d0,0.0d0,1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,0.0d0,
     &     0.0d0,0.0d0,1.0d0,1.0d0,3.46410162d0,
     &     1.0d0,2.0d0,0.0d0,0.0d0,
     &     0.0d0,4.8989795d0,1.0d0,6.3245553d0,-2.5819889d0,6.3245553d0,
     &     2.0d0,-1.6329932d0,0.0d0,
     &     3.4641016d0,4.89897949d0,9.16515139d0,4.3204938d0,
     &   -3.4156503d0,4.3204938d0,18.330303d0,-1.6329932d0,-6.9282032d0/   
cbs   get back the real number of functions for the finite nucleus   
      if (ifinite.eq.2) ncontrac(0)=ncontrac_keep
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
      xa(4)='X1MNF-SO'
      ya(4)='Y1MNF-SO'
      ZA(4)='Z1MNF-SO'
c
cbs   clean the arrays for cartesian integrals
C
      length3=(numbalLcart*numbalLcart+numbalLcart)/2
      iloca=length3
CBS   print *, 'iloca',iloca
      IOCAX = 1
      iocay=iocax+iloca
      iocaz=iocay+iloca
      iocax2=iocaz+iloca
      iocay2=iocax2+iloca
      iocaz2=iocay2+iloca
      KLAST = IOCAZ2 + ILOCA
      IF (KLAST .GT. LWRK) CALL STOPIT('AMFI  ','CAMUL',KLAST,LFREE)
      call dzero(WRK(iocax),6*length3)
c
c
c
c
cbs   one-electron-integrals:
cbs   1. index: number of first contracted function
cbs   2. index: number of second contracted function
cbs   3. index: pointer(m1,m2)    m1< m2 otherwise change sign of integral
cbs   4. index: L-value
cbs    onecartX(mxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax),
cbs    onecartY(mxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax),
cbs    onecartZ(mxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax)
c
c
c
cbs   generate one-electron integrals for all L greater/equal 1
      if (ifinite.eq.2) charge=0d0 ! nuclear integrals are modelled for finite nucleus somewhere else
      do L=1,Lhigh   
        call contone(L,oneoverr3(1,L),onecontr(1,1,-Lmax,1,L),
     *  Lmax,contrarray(iaddtyp3(L)),nprimit(L),ncontrac(L),
     *  MxcontL,dummyone,
     *  onecartx(1,1,1,L),onecartY(1,1,1,L),onecartZ(1,1,1,L),
     *  charge,oneonly)
      Enddo 
c
cbs   ***********************************************************************
cbs   now move all integrals to one big arrays for X,Y,Z
cbs   ***********************************************************************
      do Lrun=1,Lhigh  !loop over L-values (integrals are diagonal in L)
      mrun=0
      do Msec=-Lrun,Lrun    ! cartesian M-values  (Mfirst,Msec) with 
      do Mfirst=-Lrun,Msec  ! Mfirst <= Msec (actually '=' does never appear
C
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
cbs     calculate shift to get to the beginning of the block 
           iredired=shiftIRIR((iredsec*iredsec-iredsec)/2+iredfirst)
     *       +incrlm(Mfirst,Lrun)*itotalperIR(iredsec)+
     *        incrLM(Msec,Lrun)
       if (mod(ipowx,2).eq.0.and.mod(ipowy,2).eq.1.and.
     * mod(ipowz,2).eq.1) then  
            do icartfirst=1,ncontrac(Lrun) ! loop over functions first index
            do icartsec=1,ncontrac(Lrun)   ! loop over functions second index
CBS                print *, 'iocax',iocax,iredired,icartsec
                WRK(iocax+iredired+(icartsec-1))=
     *          WRK(iocax+iredired+(icartsec-1)) 
     *          +onecartx(icartfirst,icartsec,mrun,Lrun)
                enddo
cbs             shift pointer by number of functions in IR
                iredired=iredired+itotalperIR(iredsec)
                enddo
        endif 
       if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.0.and.
     * mod(ipowz,2).eq.1) then  
                do icartfirst=1,ncontrac(Lrun) ! loop over functions first index
                do icartsec=1,ncontrac(Lrun)   ! loop over functions second index
                WRK(iocay+iredired+(icartsec-1))=
     *          WRK(iocay+iredired+(icartsec-1)) 
     *          +onecarty(icartfirst,icartsec,mrun,Lrun)
                enddo
cbs             shift pointer by number of functions in IR
                iredired=iredired+itotalperIR(iredsec)
                enddo
        endif 
       if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.1.and.
     * mod(ipowz,2).eq.0) then  
                do icartfirst=1,ncontrac(Lrun) ! loop over functions first index
                do icartsec=1,ncontrac(Lrun)   ! loop over functions second index
                WRK(iocaz+iredired+(icartsec-1))=
     *          WRK(iocaz+iredired+(icartsec-1)) 
     *          +onecartz(icartfirst,icartsec,mrun,Lrun)
                enddo
cbs             shift pointer by number of functions in IR
                iredired=iredired+itotalperIR(iredsec)
                enddo
        endif 
        elseif (iredfirst.gt.iredsec) then 
cbs     In this case, indices are exchanged with respect to former 
cbs     symmetry of blocks. Therefore, there will be a minus sign 
c
cbs     calculate shift to get to the beginning of the block 
                iredired=shiftIRIR((iredfirst*iredfirst-iredfirst)/2+
     *          iredsec)+
     *          incrLM(Msec,Lrun)*itotalperIR(iredfirst)+  
     *          incrLM(Mfirst,Lrun)
       if (mod(ipowx,2).eq.0.and.mod(ipowy,2).eq.1.and.
     * mod(ipowz,2).eq.1) then  
                do icartsec=1,ncontrac(Lrun) !loop over functions second index
                do icartfirst=1,ncontrac(Lrun) !loop over functions first index
                WRK(iocax+iredired+(icartfirst-1))=
     *          WRK(iocax+iredired+(icartfirst-1)) 
     *         -onecartx(icartsec,icartfirst,mrun,Lrun)
                enddo
cbs             shift pointer by number of functions in IR
                iredired=iredired+itotalperIR(iredfirst)
                enddo
        endif 
       if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.0.and.
     * mod(ipowz,2).eq.1) then  
                do icartsec=1,ncontrac(Lrun) !loop over functions second index
                do icartfirst=1,ncontrac(Lrun) !loop over functions first index
                WRK(iocay+iredired+(icartfirst-1))=
     *          WRK(iocay+iredired+(icartfirst-1)) 
     *         -onecarty(icartsec,icartfirst,mrun,Lrun)
                enddo
cbs             shift pointer by number of functions in IR
                iredired=iredired+itotalperIR(iredfirst)
                enddo
        endif 
       if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.1.and.
     * mod(ipowz,2).eq.0) then  
                do icartsec=1,ncontrac(Lrun) !loop over functions second index
                do icartfirst=1,ncontrac(Lrun) !loop over functions first index
                WRK(iocaz+iredired+(icartfirst-1))=
     *          WRK(iocaz+iredired+(icartfirst-1)) 
     *         -onecartz(icartsec,icartfirst,mrun,Lrun)
                enddo
                iredired=iredired+itotalperIR(iredfirst)
                enddo
        endif 
      endif 
      enddo
      enddo
      enddo
C
C
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
     *          ired1)
      if (ired1.ne.ired2) then 
        WRK(iocax2+irun-1)=WRK(iocax-1+iredired+norbsh2+
     * (norbsH1-1)*itotalperIR(IREd2)) 
        WRK(iocay2+irun-1)=WRK(iocay-1+iredired+norbsh2+
     * (norbsH1-1)*itotalperIR(IREd2)) 
        WRK(iocaz2+irun-1)=WRK(iocaz-1+iredired+norbsh2+
     * (norbsH1-1)*itotalperIR(IREd2)) 
      else 
       WRK(iocax2+irun-1)=WRK(iocax-1+iredired+norbsh2*
     * (norbsH2-1)/2+norbsh1)                 
       WRK(iocay2+irun-1)=WRK(iocay-1+iredired+norbsh2*
     * (norbsH2-1)/2+norbsh1)                 
       WRK(iocaz2+irun-1)=WRK(iocaz-1+iredired+norbsh2*
     * (norbsH2-1)/2+norbsh1)                 
      endif 
      Enddo 
      enddo 
c     write a hermit-like file   b.s. 4.10.96   
CBS   write(6,*) 'number of orbitals ',numbalLcarT
CBS   write(6,*) 'length of triangular matrix ', length3
              write(LUPROP)  xa,numbofsym,(nrtofiperIR(I),
     *        i=1,numbofsym),
     *        numballcart,(Loffunction(I),I=1,numballcart), 
     *        (Moffunction(I),I=1,numballcart),
     *        Lhigh,(ncontrac(I),I=0,Lhigh)          
              write(LUPROP) (WRK(iocax2+irun),irun=0,length3-1) 
              write(LUPROP)  Ya  
              write(LUPROP) (WRK(iocay2+irun),irun=0,length3-1) 
              write(LUPROP)  Za    
              write(LUPROP) (WRK(iocaz2+irun),irun=0,length3-1) 
cbs   
cbs   that is it!!
cbs   
      return
      end
      subroutine contcasaOO(l1,l2,l3,l4,nstart,primints,
     *scratch1,scratch2,cont4OO)
cbs   contraction for powers (+2)  with alpha1*alpha3
cbs   other-orbit term 
cbs   use averaged integrals by interchanging kinematic factors   
cbs   this is case a in the documentation 
#include "implicit.h"
#include "para.h"
#include "amfi_param.h"
      dimension ncont(4),nprim(4),primints(*),scratch1(*),scratch2(*) 
     *,cont4OO(*)
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
C
C
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp2(l1)), !A *alpha                                       
     *contrarray(iaddtyp3(l2)), !A/E+m                                         
     *contrarray(iaddtyp4(l3)), !A/E+m *alpha     
     *contrarray(iaddtyp1(l4)), !A                                       
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=0.25d0*scratch1(irun)
      enddo 
C
C
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp4(l1)), 
     *contrarray(iaddtyp3(l2)), 
     *contrarray(iaddtyp2(l3)), 
     *contrarray(iaddtyp1(l4)), 
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+0.25d0*
     *scratch1(irun)
      enddo 
C
C
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp2(l1)),
     *contrarray(iaddtyp1(l2)),
     *contrarray(iaddtyp4(l3)),
     *contrarray(iaddtyp3(l4)),
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+0.25d0*
     *scratch1(irun)
      enddo 
C
C
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp4(l1)), 
     *contrarray(iaddtyp1(l2)),
     *contrarray(iaddtyp2(l3)),
     *contrarray(iaddtyp3(l4)),
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+0.25d0*
     *scratch1(irun)
      enddo 
      return 
      end
      subroutine contcasaSO(l1,l2,l3,l4,nstart,primints,
     *scratch1,scratch2,cont4SO)
cbs   contraction for powers (+2)  with alpha1*alpha3
cbs   same orbit term 
cbs   this is case a in the documentation 
#include "implicit.h"
#include "para.h"
#include "amfi_param.h"
      dimension ncont(4),nprim(4),primints(*),scratch1(*),scratch2(*),
     *cont4SO(*) 
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
c     write(6,*) 'scratch1 ',(scratch1(I),I=1,ilength)          
c     write(6,*) 'contraction coeff'
c     write(6,*) (contrarray(iaddtyp4(l1)+I),I=0,nprim(1)-1)
c     write(6,*) (contrarray(iaddtyp1(l2)+I),I=0,nprim(2)-1)
c     write(6,*) (contrarray(iaddtyp4(l3)+I),I=0,nprim(3)-1)
c     write(6,*) (contrarray(iaddtyp1(l4)+I),I=0,nprim(4)-1)
      call contract(
     *contrarray(iaddtyp4(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp4(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
c     write(6,*) 'nstart ',nstart 
      do irun=1,nprod
      cont4SO(nstart+irun-1)=scratch1(irun)
      enddo 
      return 
      end
      subroutine contcasb1OO(l1,l2,l3,l4,nstart,primints,
     *scratch1,scratch2,cont4OO)
cbs   contraction for powers (0)  with alpha1
cbs   this is one of the cases b in the documentation 
cbs   use averaged integrals by interchanging kinematic factors 
#include "implicit.h"
#include "para.h"
#include "amfi_param.h"
      dimension ncont(4),nprim(4),primints(*),scratch1(*),scratch2(*) 
     *,cont4OO(*)
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
C
C
c
cbs   copy primitive integrals to scratch1
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp2(l1)),    
     *contrarray(iaddtyp3(l2)),  
     *contrarray(iaddtyp3(l3)),  
     *contrarray(iaddtyp1(l4)), 
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=0.25d0*scratch1(irun)
      enddo
C
C
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp4(l1)),
     *contrarray(iaddtyp3(l2)),
     *contrarray(iaddtyp1(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
C
C
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp2(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp3(l3)),                                        
     *contrarray(iaddtyp3(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
C
C
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp4(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp1(l3)),                                        
     *contrarray(iaddtyp3(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
      return 
      end
      subroutine contcasb1SO(l1,l2,l3,l4,nstart,primints,
     *scratch1,scratch2,cont4SO)
cbs   contraction for powers (0)  with alpha1
cbs   this is one of the cases b in the documentation 
#include "implicit.h"
#include "para.h"
#include "amfi_param.h"
      dimension ncont(4),nprim(4),primints(*),scratch1(*),scratch2(*),
     *cont4SO(*) 
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
cbs   copy primitive integrals to scratch1
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp4(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp3(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      call dcopy(nprod,scratch1(1),1,cont4SO(nstart),1)
      return 
      end
      subroutine contcasb2OO(l1,l2,l3,l4,nstart,primints,
     *scratch1,scratch2,cont4OO)
cbs   contraction for powers (0)  with alpha3
cbs   this is one of the cases b in the documentation 
cbs   use averaged integrals by interchanging kinematic factors 
#include "implicit.h"
#include "para.h"
#include "amfi_param.h"
      dimension ncont(4),nprim(4),primints(*),scratch1(*),scratch2(*) 
     *,cont4OO(*)
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp1(l1)),                                        
     *contrarray(iaddtyp3(l2)),                                        
     *contrarray(iaddtyp4(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=0.25d0*scratch1(irun)
      enddo
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp3(l1)),                                        
     *contrarray(iaddtyp3(l2)),                                        
     *contrarray(iaddtyp2(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp1(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp4(l3)),                                        
     *contrarray(iaddtyp3(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp3(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp2(l3)),                                        
     *contrarray(iaddtyp3(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
      return 
      end
      subroutine contcasb2SO(l1,l2,l3,l4,nstart,primints,
     *scratch1,scratch2,cont4SO)
cbs   contraction for powers (0)  with alpha3
cbs   this is one of the cases b in the documentation 
#include "implicit.h"
#include "para.h"
#include "amfi_param.h"
      dimension ncont(4),nprim(4),primints(*),scratch1(*),scratch2(*),
     *cont4SO(*) 
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp3(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp4(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      call dcopy(nprod,scratch1(1),1,cont4SO(nstart),1)
      return 
      end
      SUBroutine contcascOO(l1,l2,l3,l4,nstart,primints,
     *scratch1,scratch2,cont4OO)
cbs   contraction for powers (-2)  with factor 1 
cbs   this is case c in the documentation 
cbs   use averaged integrals by interchanging kinematic factors 
#include "implicit.h"
#include "para.h"
#include "amfi_param.h"
      dimension ncont(4),nprim(4),primints(*),scratch1(*),scratch2(*) 
     *,cont4OO(*)
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp1(l1)),                                        
     *contrarray(iaddtyp3(l2)),                                        
     *contrarray(iaddtyp3(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=0.25d0*scratch1(irun)
      enddo
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp3(l1)),                                        
     *contrarray(iaddtyp3(l2)),                                        
     *contrarray(iaddtyp1(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp1(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp3(l3)),                                        
     *contrarray(iaddtyp3(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
c
c
C
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp3(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp1(l3)),                                        
     *contrarray(iaddtyp3(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      do irun=1,nprod
      cont4OO(nstart+irun-1)=cont4OO(nstart+irun-1)+
     *0.25d0*scratch1(irun)
      enddo
      return                           
      end
      subroutine contcascSO(l1,l2,l3,l4,nstart,primints,
     *scratch1,scratch2,cont4SO)
cbs   contraction for powers (-2)  with factor 1 
cbs   this is case c in the documentation 
#include "implicit.h"
#include "para.h"
#include "amfi_param.h"
      dimension ncont(4),nprim(4),primints(*),scratch1(*),scratch2(*),
     *cont4SO(*) 
      ncont(1)=ncontrac(l1)
      ncont(2)=ncontrac(l2)
      ncont(3)=ncontrac(l3)
      ncont(4)=ncontrac(l4)
      nprod=ncont(1)*ncont(2)*ncont(3)*ncont(4)
      nprim(1)=nprimit(l1)
      nprim(2)=nprimit(l2)
      nprim(3)=nprimit(l3)
      nprim(4)=nprimit(l4)
      ilength=nprim(1)*nprim(2)*nprim(3)*nprim(4)
cbs   copy primitive integrals to scratch1
      do IRUN=1,ilength
      scratch1(IRUN)=primints(IRUN)
      enddo 
      call contract(
     *contrarray(iaddtyp3(l1)),                                        
     *contrarray(iaddtyp1(l2)),                                        
     *contrarray(iaddtyp3(l3)),                                        
     *contrarray(iaddtyp1(l4)),                                        
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
     *scratch1,scratch2)
      call dcopy(nprod,scratch1(1),1,cont4SO(nstart),1)
      return 
      end
      subroutine contone(L,oneoverr3,onecontr,Lmax,
     *contcoeff,nprim,ncont,MxcontL,dummy,
     *onecartx,onecartY,onecartZ,charge,oneonly)
cbs   contracts one-electron integrals and multiplies with l,m-dependent 
cbs   factors for L-,L0,L+ 
#include "implicit.h"
      dimension oneoverR3(*),
     *onecontr(MxcontL,MxcontL,-Lmax:Lmax,3),
     *contcoeff(nprim,ncont),dummy(ncont,ncont),
     *onecartx(MxcontL,MxcontL,
     *(Lmax+Lmax+1)*(Lmax+1)),  
     *onecarty(MxcontL,MxcontL,
     *(Lmax+Lmax+1)*(Lmax+1)),  
     *onecartz(MxcontL,MxcontL,
     *(Lmax+Lmax+1)*(Lmax+1))   
      logical oneonly  
      ipnt(I,J)=(max(i,j)*(max(i,j)-1))/2+min(i,j)
cbs   first of all cleaning dummy and onecontr
      do jrun=1,ncont
      do irun=1,ncont
      dummy(irun,jrun)=0d0   
      enddo 
      enddo 
      if (oneonly) then
      iprod=MxcontL*MxcontL*(Lmax+Lmax+1)*(Lmax+1)
      call dzero(onecartx,iprod)
      call dzero(onecarty,iprod)
      call dzero(onecartz,iprod)
      endif 
      iprod=3*(Lmax+lmax+1)*MxcontL*MxcontL
      call dzero(onecontr,iprod)
cbs   contract onto dummy 
      do icont2=1,ncont
      do icont1=1,ncont 
      do iprim2=1,nprim
      do iprim1=1,nprim
      dummy(icont1,icont2)=dummy(icont1,icont2)+
     *contcoeff(iprim1,icont1)*contcoeff(iprim2,icont2)*
     *oneoverR3(ipnt(iprim1,iprim2))
      enddo
      enddo
      enddo 
      enddo 
      do icont2=1,ncont
      do icont1=1,ncont 
      dummy(icont1,icont2)=dummy(icont1,icont2)*charge 
      enddo
      enddo
cbs   start to add l,m dependent factors    
      do M=-L,L 
      factormin=dsqrt(dfloat(L*L-M*M+L+M))
      factor0=dfloat(M)
      factorplus=dsqrt(dfloat(L*L-M*M+L-M))
      do irun=1,ncont 
      do jrun=1,ncont  
      onecontr(irun,jrun,M,1)=dummy(jrun,irun)*factormin  ! L-minus 
      enddo
      enddo
      do irun=1,ncont 
      do jrun=1,ncont  
      onecontr(irun,jrun,M,2)=dummy(jrun,irun)*factor0    ! L-0     
      enddo
      enddo
      do irun=1,ncont 
      do jrun=1,ncont  
      onecontr(irun,jrun,M,3)=dummy(jrun,irun)*factorplus ! L-plus  
      enddo
      enddo
      enddo
cbs   make the final cartesian integrals 
      call cartoneX(L,Lmax,onecontr,ncont,
     *MxcontL,onecartX(1,1,1))
      call cartoneY(L,Lmax,onecontr,ncont,
     *MxcontL,onecartY(1,1,1))
      call cartoneZ(L,Lmax,onecontr,ncont,
     *MxcontL,onecartZ(1,1,1))
      return 
      end
      subroutine contract( coeffs1, coeffs2, coeffs3, coeffs4,
     *  ncont, nprim, arr1, arr2 )
c coeffs1, !(nprim(1),ncont(1)) modified contraction coefficients
c coeffs2, !(nprim(2),ncont(2)) modified contraction coefficients
c coeffs3, !(nprim(3),ncont(3)) modified contraction coefficients
c coeffs4, !(nprim(4),ncont(4)) modified contraction coefficients
c ncont,   ! i-th element is number of contracted functions i. index
c nprim,   ! i-th element is number of primitive functions  i. index
cbs  array one contains at the beginning the uncontracted integrals 
c arr1,  ! array of size (nprim(1)*nprim(2)*nprim(3)*nprim(4))
c arr2   ! array of size (nprim(1)*nprim(2)*nprim(3)*nprim(4))
#include "implicit.h"
      dimension coeffs1(*),coeffs2(*),coeffs3(*),coeffs4(*),   
     *arr1(*),arr2(*),ncont(4),nprim(4),nolds(4),nnew(4)
C
cbs   makes four indextransformations in a row....
cbs   try to find out, which indices should be transformed first...
c
      ratio1=dfloat(nprim(1))/dfloat(ncont(1))
      ratio2=dfloat(nprim(2))/dfloat(ncont(2))
      ratio3=dfloat(nprim(3))/dfloat(ncont(3))
      ratio4=dfloat(nprim(4))/dfloat(ncont(4))
      do IBM=1,4
      nolds(IBM)=nprim(IBM)  
      nnew(IBM)=nprim(IBM)  
      enddo
cbs   determine first, second,third and last index
cbs   determine the first
      xmax=max(ratio1,ratio2,ratio3,ratio4)
      if (xmax.eq.ratio1) then 
      ifirst=1
      ratio1=0
      nnew(ifirst)=ncont(ifirst)
      call trans(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
      else if (xmax.eq.ratio2) then 
      ifirst=2
      ratio2=0
      nnew(ifirst)=ncont(ifirst)
      call trans(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      else if (xmax.eq.ratio3) then 
      ifirst=3
      ratio3=0
      nnew(ifirst)=ncont(ifirst)
      call trans(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      else if (xmax.eq.ratio4) then 
      ifirst=4
      ratio4=0
      nnew(ifirst)=ncont(ifirst)
      call trans(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      endif
      nolds(ifirst)=nnew(ifirst)
cbs   determine the second   
      xmax=max(ratio1,ratio2,ratio3,ratio4)
      if (xmax.eq.ratio1) then 
      isec=1
      ratio1=0
      nnew(isec)=ncont(isec)
      call trans(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      else if (xmax.eq.ratio2) then 
      isec=2
      ratio2=0
      nnew(isec)=ncont(isec)
      call trans(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      else if (xmax.eq.ratio3) then 
      isec=3
      ratio3=0
      nnew(isec)=ncont(isec)
      call trans(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      else if (xmax.eq.ratio4) then 
      isec=4
      ratio4=0
      nnew(isec)=ncont(isec)
      call trans(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      endif
      nolds(isec)=nnew(isec)
cbs   determine the third    
      xmax=max(ratio1,ratio2,ratio3,ratio4)
      if (xmax.eq.ratio1) then 
      ithird=1
      ratio1=0
      nnew(ithird)=ncont(ithird)
      call trans(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      else if (xmax.eq.ratio2) then 
      ithird=2
      ratio2=0
      nnew(ithird)=ncont(ithird)
      call trans(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      else if (xmax.eq.ratio3) then 
      ithird=3
      ratio3=0
      nnew(ithird)=ncont(ithird)
      call trans(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      else if (xmax.eq.ratio4) then 
      ithird=4
      ratio4=0
      nnew(ithird)=ncont(ithird)
      call trans(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)     
      endif
      nolds(ithird)=nnew(ithird)
cbs   determine the last    
      xmax=max(ratio1,ratio2,ratio3,ratio4)
      if (xmax.eq.ratio1) then 
      ifourth=1
      ratio1=0
      nnew(ifourth)=ncont(ifourth)    
      call trans(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      else if (xmax.eq.ratio2) then 
      ifourth=2
      ratio2=0
      nnew(ifourth)=ncont(ifourth)    
      call trans(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      else if (xmax.eq.ratio3) then 
      ifourth=3
      ratio3=0
      nnew(ifourth)=ncont(ifourth)    
      call trans(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      else if (xmax.eq.ratio4) then 
      ifourth=4
      ratio4=0
      nnew(ifourth)=ncont(ifourth)    
      call trans(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),
     *nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)    
      endif
cbs   contracted integrals are now on 
cbs   arr1(ncont1,ncont2,ncont3,ncont4)   
      return 
      end
      double precision function  couple3J(
     *l1, l2, l3, m1, m2, m3)
cbs this routine calculates the coupling of three angular momenta to  zero   
cbs   
cbs   
cbs   Int dOmega i^(l1+l2+l3) Y^l1_m1 (Omega) Y^l2_m2 (Omega) Y^l3_m3 (Omega) =
cbs   sqrt( (2l1+1)(2l2+1)(2l2+3)/ 4Pi)  * 3J(l1,l2,l3,0,0,0) *
cbs   3J(l1,l2,l3,m1,m2,m3)
cbs   
cbs   
#include "implicit.h"
#include "pi.h"
      double precision inv4pi
cbs   (4*PI)**-1 
      inv4pi=0.25d0/pi                                   
cbs   initialize couple3J-coefficient
      couple3J=0d0
cbs   quick check 
      if (m1+m2+m3.ne.0) return 
cbs   double all values for regge3j   
      l1d=l1+l1
      l2d=l2+l2
      l3d=l3+l3
      m1d=m1+m1
      m2d=m2+m2
      m3d=m3+m3
      fac1=dsqrt(dfloat(l1d+1)*dfloat(l2d+1)*dfloat(l3d+1)*inv4pi)
      fac2=regge3j(l1d,l2d,l3d,0,0,0)
      fac3=regge3j(l1d,l2d,l3d,m1d,m2d,m3d)
      couple3J=fac1*fac2*fac3
      return 
      end 
      subroutine daxpint(from,to,fact,ndim1,ndim2,ndim3,ndim4)
#include "implicit.h"
cbs   subroutine similar to daxpy with interchange of two indices 
cbs   change from physicists notation to chemists notaion 
cbs   to(i,j,k,l)=to(i,j,k,l)+fact*from(i,k,j,l) 
      dimension from(ndim1,ndim2,ndim3,ndim4),
     *to(ndim1,ndim3,ndim2,ndim4)
      if (fact.eq.0d0) return 
      do irun4=1,ndim4
      do irun3=1,ndim3
      do irun2=1,ndim2
      do irun1=1,ndim1
      to(irun1,irun3,irun2,irun4)=to(irun1,irun3,irun2,irun4)+
     *fact*from(irun1,irun2,irun3,irun4)  
      enddo
      enddo
      enddo
      enddo
      return   
      end 
      subroutine gen1overR3(Lhigh)
#include "implicit.h"
cbs   generates the radial integrals  for the one electron spin orbit integrals   
cbs   taken the 1/r**3 formula from the documentation and included additional 
cbs   factors for normalization   
#include "para.h"
#include "amfi_param.h"
#include "dofuc.h"
#include "pi.h"
      do L=1,Lhigh    
      icount=0
      do iprim2=1,nprimit(L)
      alpha2=exponents(iprim2,L)
      do iprim1=1,iprim2    
      alpha1=exponents(iprim1,L)
      icount=icount+1
      oneoverR3(icount,L)=dsqrt(2d0/pi)*
     *(df(L+L-2)*2**(L+3)*
     *(alpha1*alpha2)**(0.25d0*
     *(L+L+3)))/((alpha1+alpha2)**L*df(L+L+1))
      enddo
      enddo
      enddo 
      return 
      end   
      subroutine gencoul(l1,l2,l3,l4,makemean,
     *bonn,breit,sameorb,cont4SO,cont4OO,icont4,
     *WRK,LFREE)
#include "implicit.h"
cbs   SUBROUTINE to generate all required radial 
cbs   integrals for the four angular momenta l1-l4
#include "priunit.h"
#include "para.h"
#include "amfi_param.h"
      logical makemean,bonn,breit,sameorb 
      dimension cont4SO(*),cont4OO(*),WRK(LFREE)
      max1=1  !starting values for limits of precalculated 
c             ! powers of function Cfunct(X)
      max2=1
cbs   first of all, this routine determines, for which L 
cbs   values the radial integrals have to be solved
cbs   initialize the number of blocks for the different 
cbs   l-combinations  
cbs   no (ss|ss) contributions 
      if (l1.eq.0.and.l2.eq.0.and.l3.eq.0.and.l4.eq.0) return  ! no integrals for <ss|ss> 
      if (makemean) then 
                nblock=1  ! sp sp are the first, so the first block
                Lstarter(1)=1
      else 
      CALL QUIT('only mean-field with this version')
      endif 
cbs   keep track of L-values for later purposes   
      Lvalues(1)=l1
      Lvalues(2)=l2
      Lvalues(3)=l3
      Lvalues(4)=l4
cbs   now nanz is given the new value
      nanz=ncontrac(l1)*ncontrac(l2)*ncontrac(l3)*ncontrac(l4)
      nprimprod=nprimit(l1)*nprimit(l2)*nprimit(l3)*nprimit(l4)
      IQUOT1 = 1
      iquot2=iquot1+nprimprod
      iquotp1=iquot2+nprimprod
      iquotp2=iquotp1+nprimprod
      iprim=iquotp2+nprimprod
      iscr1=iprim+nprimprod
      iscr2=iscr1+nprimprod
      KLAST = ISCR2 + NPRIMPROD
      IF (KLAST .GT. LFREE) CALL STOPIT('AMFI  ','GENCOU',KLAST,LFREE)
c
      call initfrac(nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4),
     *WRK(iquot1),WRK(iquot2),exponents(1,l1),exponents(1,l2),
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
        Lblocks(1)=(Lend-Lanf)/2+1 ! L increases in steps of 2, 
cbs                                       due to parity conservation 
        Lfirst(1)=Lanf
        Llast(1)=Lend 
      else 
        Lblocks(1)=0 
      endif 
      if (Lblocks(1).gt.0) then    ! integrals have to be calculated 
cbs### check, whether integrals fit on array ################
      if  (Lstarter(1)+nanz*Lblocks(1).gt.icont4) then 
      write(LUPRI,*) 'end at: ',Lstarter(1)+nanz*Lblocks(1) 
      CALL QUIT('increase icont4 in amfi.F')
      endif 
cbs### check, whether integrals fit on array ################
      istart=Lstarter(1) ! gives the address, where to write the contracted integrals 
cbs   ipow1 and ipow2 are the the numbers of powers in the prefactor
cbs   of the function Cfunct 
cbs   now loop over possible L-values 
      do Lrun= Lfirst(1),Llast(1),2
                ipow1=2+(l2+l4+Lrun)/2
                ipow2=2+(l1+l3+incl1+incl3+Lrun)/2
cbs   those powers have to be generated... 
      call getpow(ipow1,WRK(iquot1),WRK(iquotp1),
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
cbs   those powers have to be generated... 
      call getpow(ipow2,WRK(iquot2),WRK(iquotp2),
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
c     in buildcoul the radial integrals are calculated 
                call buildcoul(l1,l2,l3,l4,incl1,incl3,
     *          Lrun,WRK(iprim),nprimit(l1),nprimit(l2),nprimit(l3),
     *          nprimit(l4), 
     *          exponents(1,l1),exponents(1,l2),
     *          exponents(1,l3),exponents(1,l4),
     *          powexp(1,1,l3,l1,lrun),powexp(1,1,l4,l2,lrun),
     *          WRK(iquotp1),WRK(iquotp2))
cbs   in the contcas_ routines the integrals are contracted, including exponents as prefactors... 
                if (bonn.or.breit.or.sameorb) then      
                call contcasASO(l1,l2,l3,l4,istart,WRK(iprim),
     *           WRK(iscr1),WRK(iscr2),cont4SO)
                else 
                call contcasASO(l1,l2,l3,l4,istart,WRK(iprim),
     *           WRK(iscr1),WRK(iscr2),cont4SO)
                call contcasAOO(l1,l2,l3,l4,istart,WRK(iprim),
     *           WRK(iscr1),WRK(iscr2),cont4OO)
                endif 
                istart=istart+nanz  ! start-address for the next block of contracted integrals 
      enddo 
      endif   
cbs   ##########################################################################################################
cbs   the (0) cases         CASE  B
cbs   ##########################################################################################################
      incl1=0
      incl3=0
      call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
      if(Lend-Lanf.ge.0) then 
      Lblocks(2)=(Lend-Lanf)/2+1
      Lfirst(2)=Lanf
      Llast(2)=Lend 
      Lblocks(3)=(Lend-Lanf)/2+1
      Lfirst(3)=Lanf
      Llast(3)=Lend 
      else 
      Lblocks(2)=0 
      Lblocks(3)=0 
      endif 
      Lstarter(2)=Lstarter(1)+
     *nanz*Lblocks(1)
      Lstarter(3)=Lstarter(2)+
     *nanz*Lblocks(2)
cbs   primitive integrals are the same for type 2 and 3  !!!!!
      if (Lblocks(2).gt.0) then    
cbs### check, whether integrals fit on array ################
      if  (Lstarter(2)+2*nanz*Lblocks(2).gt.icont4) then 
      write(LUPRI,*) 'end at: ',Lstarter(2)+2*nanz*Lblocks(2) 
      CALL QUIT('increase icont4 in amfi.F')
      endif 
cbs### check, whether integrals fit on array ################
      istart=Lstarter(2)
      istart2=Lstarter(3)
      do Lrun= Lfirst(2),Llast(2),2
      ipow1=2+(l2+l4+Lrun)/2
      ipow2=2+(l1+l3+incl1+incl3+Lrun)/2
      call getpow(ipow1,WRK(iquot1),WRK(iquotp1),
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
      call getpow(ipow2,WRK(iquot2),WRK(iquotp2),
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
      call buildcoul(l1,l2,l3,l4,incl1,incl3,
     *Lrun,WRK(iprim),nprimit(l1),nprimit(l2),nprimit(l3),
     *nprimit(l4), 
     *exponents(1,l1),exponents(1,l2),
     *exponents(1,l3),exponents(1,l4),
     *powexp(1,1,l3,l1,lrun),powexp(1,1,l4,l2,lrun),   
     *WRK(iquotp1),WRK(iquotp2))
      if (bonn.or.breit.or.sameorb) then      
      call contcasB1SO(l1,l2,l3,l4,istart,WRK(iprim),
     *WRK(iscr1),WRK(iscr2),cont4SO)
      call contcasB2SO(l1,l2,l3,l4,istart2,WRK(iprim),
     *WRK(iscr1),WRK(iscr2),cont4SO)
      else 
      call contcasB1SO(l1,l2,l3,l4,istart,WRK(iprim),
     *WRK(iscr1),WRK(iscr2),cont4SO)
      call contcasB2SO(l1,l2,l3,l4,istart2,WRK(iprim),
     *WRK(iscr1),WRK(iscr2),cont4SO)
      Call contcasB1OO(l1,l2,l3,l4,istart,WRK(iprim),
     *WRK(iscr1),WRK(iscr2),cont4OO)
      Call contcasB2OO(l1,l2,l3,l4,istart2,WRK(iprim),
     *WRK(iscr1),WRK(iscr2),cont4OO)
      endif   
      istart=istart+nanz
      istart2=istart2+nanz
      enddo 
      endif   
cbs   ##########################################################################################################
cbs   the (-2) cases      CASE C
cbs   ##########################################################################################################
      if (l1.eq.0.or.l3.eq.0) then 
      Lblocks(4)=0
      else 
      incl1=-1
      incl3=-1
      call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
      if(Lend-Lanf.ge.0) then 
      Lblocks(4)=(Lend-Lanf)/2+1
      Lfirst(4)=Lanf
      Llast(4)=Lend 
      else 
      Lblocks(4)=0 
      endif 
      endif
      Lstarter(4)=Lstarter(3)+
     *nanz*Lblocks(3)
      if (Lblocks(4).gt.0) then    
cbs### check, whether integrals fit on array ################
      if  (Lstarter(4)+nanz*Lblocks(4).gt.icont4) then 
      write(LUPRI,*) 'end at: ',Lstarter(4)+nanz*Lblocks(4) 
      CALL QUIT('increase icont4 in amfi.F')
      endif 
cbs### check, whether integrals fit on array ################
      istart=Lstarter(4)
      do Lrun= Lfirst(4),Llast(4),2
      ipow1=2+(l2+l4+Lrun)/2
      ipow2=2+(l1+l3+incl1+incl3+Lrun)/2
      call getpow(ipow1,WRK(iquot1),WRK(iquotp1),
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
      call getpow(ipow2,WRK(iquot2),WRK(iquotp2),
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
      call buildcoul(l1,l2,l3,l4,incl1,incl3,
     *Lrun,WRK(iprim),nprimit(l1),nprimit(l2),nprimit(l3),
     *nprimit(l4), 
     *exponents(1,l1),exponents(1,l2),
     *exponents(1,l3),exponents(1,l4),
     *powexp(1,1,l3,l1,lrun),powexp(1,1,l4,l2,lrun),
     *WRK(iquotp1),WRK(iquotp2))
      if (bonn.or.breit.or.sameorb) then      
      call contcasCSO(l1,l2,l3,l4,istart,WRK(iprim),
     *WRK(iscr1),WRK(iscr2),cont4SO)
      else 
      call contcasCSO(l1,l2,l3,l4,istart,WRK(iprim),
     *WRK(iscr1),WRK(iscr2),cont4SO)
      call contcasCOO(l1,l2,l3,l4,istart,WRK(iprim),
     *WRK(iscr1),WRK(iscr2),cont4OO)
      endif 
      istart=istart+nanz
      enddo 
      endif   
      return 
      end  
      subroutine gencoulDIM(l1,l2,l3,l4,makemean,
     *bonn,breit,sameorb,icont4)
#include "implicit.h"
#include "priunit.h"
#include "para.h"
#include "amfi_param.h"
cbs   SUBROUTINE to calculate the dimemsion of the radial integral 
cbs   arrays. BASICALLY GENCOUL WITHOUT EXPLICIT INTEGRAL CALCULATION 
cbs   integrals for the four angular momenta l1-l4
      logical makemean,bonn,breit,sameorb 
      max1=1  !starting values for limits of precalculated 
c             ! powers of function Cfunct(X)
      max2=1
c
      incont4=0
c
cbs   first of all, this routine determines, for which L 
cbs   values the radial integrals have to be solved
cbs   initialize the number of blocks for the different 
cbs   l-combinations  
cbs   no (ss|ss) contributions 
      if (l1.eq.0.and.l2.eq.0.and.l3.eq.0.and.l4.eq.0) return  ! no integrals for <ss|ss> 
      if (makemean) then 
                nblock=1  ! sp sp are the first, so the first block
                Lstarter(1)=1
      else 
      CALL QUIT('only mean-field with this version')
      endif 
cbs   keep track of L-values for later purposes   
      Lvalues(1)=l1
      Lvalues(2)=l2
      Lvalues(3)=l3
      Lvalues(4)=l4
cbs   now nanz is given the new value
      nanz=ncontrac(l1)*ncontrac(l2)*ncontrac(l3)*ncontrac(l4)
      nprimprod=nprimit(l1)*nprimit(l2)*nprimit(l3)*nprimit(l4)
c
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
        Lblocks(1)=(Lend-Lanf)/2+1 ! L increases in steps of 2, 
cbs                                       due to parity conservation 
        Lfirst(1)=Lanf
        Llast(1)=Lend 
      else 
        Lblocks(1)=0 
      endif 
cbs   ##########################################################################################################
cbs   the (0) cases         CASE  B
cbs   ##########################################################################################################
      incl1=0
      incl3=0
      call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
      if(Lend-Lanf.ge.0) then 
      Lblocks(2)=(Lend-Lanf)/2+1
      Lfirst(2)=Lanf
      Llast(2)=Lend 
      Lblocks(3)=(Lend-Lanf)/2+1
      Lfirst(3)=Lanf
      Llast(3)=Lend 
      else 
      Lblocks(2)=0 
      Lblocks(3)=0 
      endif 
      Lstarter(2)=Lstarter(1)+
     *nanz*Lblocks(1)
      Lstarter(3)=Lstarter(2)+
     *nanz*Lblocks(2)
cbs   ##########################################################################################################
cbs   the (-2) cases      CASE C
cbs   ##########################################################################################################
      if (l1.eq.0.or.l3.eq.0) then 
      Lblocks(4)=0
      else 
      incl1=-1
      incl3=-1
      call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
      if(Lend-Lanf.ge.0) then 
      Lblocks(4)=(Lend-Lanf)/2+1
      Lfirst(4)=Lanf
      Llast(4)=Lend 
      else 
      Lblocks(4)=0 
      endif 
      endif
      Lstarter(4)=Lstarter(3)+
     *nanz*Lblocks(3)
c
CBS   now the hole purpose of this routine 
c
      icont4=Lstarter(4)+nanz*Lblocks(4)
      return 
      end  
      subroutine genovlp(Lhigh)
#include "implicit.h"
#include "para.h"
#include "amfi_param.h"
cbs   generates overlap of normalized  primitives. 
      dimension evecinv(MxprimL,MxprimL)
      do L=0,Lhigh 
        do Jrun=1,nprimit(L)
        do Irun=1,nprimit(L)
        normovlp(Irun,Jrun,L)=coulovlp(irun,jrun,0,0,
     *  L,L)   
        enddo 
        enddo 
cbs   invert the matrix, not very elegant, but sufficient
      ipnt=0
      do jrun=1,nprimit(L)
      do irun=1,jrun
      ipnt=ipnt+1
      scratchinv(ipnt)=normovlp(irun,jrun,L)
      enddo
      enddo
      do Jrun=1,nprimit(L)
      do Irun=1,nprimit(L)
      evecinv(Irun,Jrun)=0d0                                
      enddo 
      enddo 
      do Jrun=1,nprimit(L)
      evecinv(jrun,jrun)=1d0
      enddo 
      call jacobi(scratchinv,evecinv,nprimit(L),MxprimL) 
      do irun=1,nprimit(L)
      eval(irun)=dsqrt(scratchinv((irun*irun+irun)/2))
      enddo
cbs   ensure normalization of the vectors.
      do IRUN=1,nprimit(L)
      fact=0d0
      do JRUN=1,nprimit(L)
      fact=fact+evecinv(JRUN,IRUN)*evecinv(JRUN,IRUN)
      enddo
      fact=1d0/dsqrt(fact)
      do JRUN=1,nprimit(L)
      evecinv(JRUN,IRUN)=fact*evecinv(JRUN,IRUN)
      enddo
      enddo  
cbs   now generate rootOVLP           
      do irun=1,nprimit(L)
      do jrun=1,nprimit(L)
      rootOVLP(irun,jrun,l)=0d0
      enddo  
      enddo  
      do jrun=1,nprimit(L)
      do irun=1,nprimit(L)
      do krun=1,nprimit(L)
      rootOVLP(irun,jrun,L)=rootOVLP(irun,jrun,L)+
     *evecinv(irun,krun)*evecinv(jrun,krun)*eval(krun)
      enddo  
      enddo  
      enddo  
cbs   now generate rootOVLPinv           
      do irun=1,nprimit(L)
      eval(irun)=1d0/eval(irun)                            
      enddo
      do irun=1,nprimit(L)
      do jrun=1,nprimit(L)
      rootOVLPinv(irun,jrun,l)=0d0
      enddo  
      enddo  
      do jrun=1,nprimit(L)
      do irun=1,nprimit(L)
      do krun=1,nprimit(L)
      rootOVLPinv(irun,jrun,L)=rootOVLPinv(irun,jrun,L)+
     *evecinv(irun,krun)*evecinv(jrun,krun)*eval(krun)
      enddo  
      enddo  
      enddo  
cbs   now generate OVLPinv           
      do irun=1,nprimit(L)
      eval(irun)=eval(irun)*eval(irun)
      enddo
      do irun=1,nprimit(L)
      do jrun=1,nprimit(L)
      OVLPinv(irun,jrun,l)=0d0
      enddo  
      enddo  
      do jrun=1,nprimit(L)
      do irun=1,nprimit(L)
      do krun=1,nprimit(L)
      OVLPinv(irun,jrun,L)=OVLPinv(irun,jrun,L)+
     *evecinv(irun,krun)*evecinv(jrun,krun)*eval(krun)
      enddo  
      enddo  
      enddo  
      enddo 
      return 
      end   
      subroutine genpowers(Lhigh)   
#include "implicit.h"
#include "para.h"
#include "amfi_param.h"
#include "dofuc.h"
cbs   set some often used powers of exponents 
      do L2=0,Lhigh 
      do L1=0,L2    
      do irun1=1,nprimit(L1)
      do irun2=1,nprimit(L2)
      powexp(irun1,irun2,L1,L2,0)=1d0
      enddo
      enddo
      enddo
      enddo
      do L2=0,Lhigh 
      do L1=0,L2    
      do Lrun=1,(L1+L2+5)
      do irun2=1,nprimit(L2)
      do irun1=1,nprimit(L1)
      fact=dsqrt(0.5d0*(exponents(irun1,L1)+exponents(irun2,L2)))
      powexp(irun1,irun2,L1,L2,Lrun)= powexp(irun1,irun2,L1,L2,Lrun-1)*
     *fact
      enddo
      enddo
      enddo
      enddo
      enddo
cbs   generate coulovlp = overlap for normalized functions, but sometimes
cbs   with shifted l-values
      do l2=0,lhigh
      do incl2=-1,1
         if (l2+incl2.ge.0) then  ! do not lower l for s-functions 
         n2=l2+incl2+1
         df2=1d0/dsqrt(df(n2+n2-1))
         do l1=0,l2
         do incl1=-1,1
         if (l1+incl1.ge.0) then ! do not lower l for s-functions 
         n1=l1+incl1+1
         df1=1d0/dsqrt(df(n1+n1-1))
         df12=df(n1+n2-1)
         do iprim2=1,nprimit(l2)
         fact2=dsqrt(powexp(iprim2,iprim2,l2,l2,n2+n2+1))
         factor=fact2*df1*df2*df12
         do iprim1=1,nprimit(l1)
         fact1=dsqrt(powexp(iprim1,iprim1,l1,l1,n1+n1+1))
         coulovlp(iprim1,iprim2,incl1,incl2,l1,l2)= 
     *   fact1*factor/powexp(iprim1,iprim2,l1,l2,n1+n2+1) 
         enddo
         enddo
         endif
         enddo
         enddo
         endif
      enddo
      enddo
      return 
      end 
 
    
      subroutine genstar(Lhigh) 
#include "implicit.h"
#include "para.h"
#include "amfi_param.h"
cbs   purpose: generate start adresses of contraction coeffs on 
cbs   contrarray for the different L-Blocks 
      istart=1
      do L=0,Lhigh
      inc=nprimit(L)*ncontrac(L)
      iaddori(L)=istart
      istart=istart+inc  
      iaddtyp1(L)=istart 
      istart=istart+inc  
      iaddtyp2(L)=istart 
      istart=istart+inc  
      iaddtyp3(L)=istart 
      istart=istart+inc  
      iaddtyp4(L)=istart 
      istart=istart+inc  
      enddo 
      return 
      end  
      subroutine gentkin(L,TKIN,nprims,exponents,rootOVLPinv)   
#include "implicit.h"
#include "para.h"
cbs   subroutine to generate the kinetic energy 
      dimension TKIN(nprims,nprims),exponents(*),
     *dummy(MxprimL,MxprimL),dummy2(MxprimL,MxprimL),
     *rootOVLPinv(MxprimL,MxprimL)
cbs   one triangular part of the matrix
      do irun2=1,nprims  
      do irun1=1,irun2   
        dummy(irun1,irun2)=
     *  Tkinet(l,exponents(irun1),
     *  exponents(irun2))
      enddo     
      enddo     
cbs   copy to the other triangular part.... 
      do irun2=1,nprims-1 
      do irun1=irun2+1,nprims      
        dummy(irun1,irun2)=dummy(irun2,irun1)                       
      enddo     
      enddo     
cbs   now transform by rootovlp*dummy*rootovlp 
      do jrun=1,nprims
      do irun=1,nprims
        TKIN(irun,jrun)=0d0
        dummy2(irun,jrun)=0d0
      enddo
      enddo
      do irun=1,nprims
      do jrun=1,nprims
      do krun=1,nprims
        dummy2(irun,jrun)=dummy2(irun,jrun)+
     *  dummy(irun,krun)*rootovlpinv(krun,jrun)
      enddo
      enddo
      enddo
      do irun=1,nprims
      do jrun=1,nprims
      do krun=1,nprims
        Tkin(irun,jrun)=Tkin(irun,jrun)+
     *  dummy2(krun,jrun)*rootovlpinv(irun,krun)
      enddo
      enddo
      enddo
      return 
      end  
      subroutine getAOs(lhigh)                  
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "para.h"
#include "amfi_param.h"
cbs   get expansions of atomic orbitals in contracted functions 
      character*12    occtext,occread
      character*18  textnorbmf,textnorbmf2    
      logical EX
      occtext='OCCUPATION: ' 
      textnorbmf='Number of orbitals'
      Inquire(File='AO-expansion',exist=EX) 
      if (.not.EX)  then
CBS   write(6,*) 'get occupations from DATA-block'
      call getAOs2(lhigh)
      return
      endif 
      LUAOEX = -1
      CALL GPOPEN(LUAOEX,'AO-expansion','UNKNOWN',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      write(LUPRI,*) 'Orbitals for mean-field' 
      do lrun=0,lhigh 
      write(LUPRI,'(A3,I3)') 'L= ',lrun 
      read(LUAOEX,'(A18,I3)') textnorbmf2,noccorb(lrun)
      if (textnorbmf.ne.textnorbmf2) 
     *CALL QUIT('wrong keyword for number of orbitals in getAOs')
      write(LUPRI,*) 'number of orbitals ',noccorb(lrun)
      do iorbital=1,noccorb(lrun)     
      read(LUAOEX,'(A12,F5.3)')  occread,occup(iorbital,lrun)
      write(LUPRI,'(A,F8.4)') occtext,occup(iorbital,lrun)
      if (occread.ne.occtext) CALL QUIT('error reading AOs')
      read(LUAOEX,*) (AOcoeffs(icont,iorbital,lrun),
     *icont=1,ncontrac(lrun))
      write(LUPRI,'(8F10.4)') (AOcoeffs(icont,iorbital,lrun),
     *icont=1,ncontrac(lrun))
      write(LUPRI,*) ' '
      read(LUAOEX,*) 
      enddo
      enddo
      call gpclose(LUAOEX,'KEEP')
      return 
      end 
      subroutine getAOs2(lhigh)
#include "implicit.h"
#include "para.h"
#include "amfi_param.h"
cbs   get expansions of atomic orbitals in contracted functions
      common /nucleus/ charge,Exp_finite
      character*12    occtext
      integer closedshells(0:LMAX),openshells(0:LMAX)
      call getocc_ao(int(charge),closedshells,openshells)
      occtext='OCCUPATION: '
      do lrun=0,lhigh
      do irun=1,MxcontL
      do jrun=1,MxcontL
      AOcoeffs(jrun,irun,lrun)=0d0
      enddo
      enddo
      enddo
CBS   write(6,*) 'Orbitals for mean-field'
      do lrun=0,lhigh
CBS   write(6,'(A3,I3)') 'L= ',lrun
      do i=1,closedshells(lrun)
      occup(i,lrun)=2.0
      AOcoeffs(i,i,lrun)=1d0
      enddo
      noccorb(lrun)=closedshells(lrun)
      if (openshells(lrun).gt.0) then
      i=closedshells(lrun)+1
      occup(i,lrun)=1d0*openshells(lrun)/dfloat(lrun+lrun+1)
      AOcoeffs(i,i,lrun)=1d0
      noccorb(lrun)=i
      endif
      if (noccorb(lrun).gt.0) then
CBS   write(6,'(A,I3)') 'number of orbitals ',noccorb(lrun)
CBS   do iorbital=1,noccorb(lrun)
CBS   write(6,'(A,8F8.4)') occtext,(occup(iorbital,lrun),
CBS  *iorbital=1,noccorb(lrun))
CBS   enddo
      endif
      enddo
      return
      end
cbs
      subroutine getocc_ao(icharge,iclosed,iopen)
#include "implicit.h"
#include "priunit.h"
#include "para.h"
      parameter (ichargemax=96)
      dimension iclocc(0:Lmax_occ,0:ichargemax)
      dimension iopocc(0:Lmax_occ,0:ichargemax)
      character*30 occtxt(0:96) 
      character*32 txt  
      data txt/'SO-integrals are calculated for '/
      dimension iclosed(0:LMAX),iopen(0:LMAX)
      data (occtxt(i),i=0,96) /
     *'dummy atom (no integrals)     ',
     *' H: no mean-field             ',
     *'He: 1s^2                      ',
     *'Li: [He]2s^1                  ',
     *'Be: [He]2s^2                  ',
     *' B: [He]2s^2 2p^1             ',
     *' C: [He]2s^2 2p^2             ',
     *' N: [He]2s^2 2p^3             ',
     *' O: [He]2s^2 2p^4             ',
     *' F: [He]2s^2 2p^5             ',
     *'Ne: [He]2s^2 2p^6             ',
     *'Na: [Ne]3s^1                  ',
     *'Mg: [Ne]3s^2                  ',
     *'Al: [Ne]3s^2 3p^1             ',
     *'Si: [Ne]3s^2 3p^2             ',
     *' P: [Ne]3s^2 3p^3             ',
     *' S: [Ne]3s^2 3p^4             ',
     *'Cl: [Ne]3s^2 3p^5             ',
     *'Ar: [Ne]3s^2 3p^6             ',
     *' K: [Ar]4s^1                  ',
     *'Ca: [Ar]4s^2                  ',
     *'Sc: [Ar]4s^2 3d^1             ',
     *'Ti: [Ar]4s^2 3d^2             ',
     *' V: [Ar]4s^2 3d^3             ',
     *'Cr: [Ar]4s^2 3d^4             ',
     *'Mn: [Ar]4s^2 3d^5             ',
     *'Fe: [Ar]4s^2 3d^6             ',
     *'Co: [Ar]4s^2 3d^7             ',
     *'Ni: [Ar]4s^2 3d^8             ',
     *'Cu: [Ar]4s^1 3d^10            ',
     *'Zn: [Ar]4s^2 3d^10            ',
     *'Ga: [Ar]4s^2 3d^10 4p^1       ',
     *'Ge: [Ar]4s^2 3d^10 4p^2       ',
     *'As: [Ar]4s^2 3d^10 4p^3       ',
     *'Se: [Ar]4s^2 3d^10 4p^4       ',
     *'Br: [Ar]4s^2 3d^10 4p^5       ',
     *'Kr: [Ar]4s^2 3d^10 4p^6       ',
     *'Rb: [Kr]5s^1                  ',
     *'Sr: [Kr]5s^2                  ',
     *' Y: [Kr]5s^2 4d^1             ',
     *'Zr: [Kr]5s^2 4d^2             ',
     *'Nb: [Kr]5s^2 4d^3             ',
     *'Mo: [Kr]5s^2 4d^4             ',
     *'Tc: [Kr]5s^2 4d^5             ',
     *'Ru: [Kr]5s^2 4d^6             ',
     *'Rh: [Kr]5s^2 4d^7             ',
     *'Pd: [Kr]5s^2 4d^8             ',
     *'Ag: [Kr]5s^1 4d^10            ',
     *'Cd: [Kr]5s^2 4d^10            ',
     *'In: [Kr]5s^2 4d^10 5p^1       ',
     *'Sn: [Kr]5s^2 4d^10 5p^2       ',
     *'Sb: [Kr]5s^2 4d^10 5p^3       ',
     *'Te: [Kr]5s^2 4d^10 5p^4       ',
     *' I: [Kr]5s^2 4d^10 5p^5       ',
     *'Xe: [Kr]5s^2 4d^10 5p^6       ',
     *'Cs: [Xe]6s^1                  ',
     *'Ba: [Xe]6s^2                  ',
     *'La: [Xe]6s^2 5d^1             ',
     *'Ce: [Xe]6s^2 4f^2             ',
     *'Pr: [Xe]6s^2 4f^3             ',
     *'Nd: [Xe]6s^2 4f^4             ',
     *'Pm: [Xe]6s^2 4f^5             ',
     *'Sm: [Xe]6s^2 4f^6             ',
     *'Eu: [Xe]6s^2 4f^7             ',
     *'Gd: [Xe]6s^2 4f^8             ',
     *'Tb: [Xe]6s^2 4f^9             ',
     *'Dy: [Xe]6s^2 4f^10            ',
     *'Ho: [Xe]6s^2 4f^11            ',
     *'Er: [Xe]6s^2 4f^12            ',
     *'Tm: [Xe]6s^2 4f^13            ',
     *'Yb: [Xe]6s^2 4f^14            ',
     *'Lu: [Xe+4f^14]6s^2 5d^1       ',
     *'Hf: [Xe+4f^14]6s^2 5d^2       ',
     *'Ta: [Xe+4f^14]6s^2 5d^3       ',
     *' W: [Xe+4f^14]6s^2 5d^4       ',
     *'Re: [Xe+4f^14]6s^2 5d^5       ',
     *'Os: [Xe+4f^14]6s^2 5d^6       ',
     *'Ir: [Xe+4f^14]6s^2 5d^7       ',
     *'Pt: [Xe+4f^14]6s^1 5d^9       ', 
     *'Au: [Xe+4f^14]6s^1 5d^10      ',
     *'Hg: [Xe+4f^14]6s^2 5d^10      ',
     *'Tl: [Xe+4f^14+5d^10]6s^2 6p^1 ',
     *'Pb: [Xe+4f^14+5d^10]6s^2 6p^2 ',
     *'Bi: [Xe+4f^14+5d^10]6s^2 6p^3 ',
     *'Po: [Xe+4f^14+5d^10]6s^2 6p^4 ',
     *'At: [Xe+4f^14+5d^10]6s^2 6p^5 ',
     *'Rn: [Xe+4f^14+5d^10]6s^2 6p^6 ',
     *'Fr: [Rn]7s^1                  ',
     *'Ra: [Rn]7s^2                  ',
     *'Ac: [Rn]7s^2 6d^1             ',
     *'Th: [Rn]7s^2 6d^2             ',
     *'Pa: [Rn]7s^2 6d^1 5f^2        ',
     *' U: [Rn]7s^2 6d^1 5f^3        ',
     *'Np: [Rn]7s^2 6d^1 5f^4        ',
     *'Pu: [Rn]7s^2 6d^0 5f^6        ',
     *'Am: [Rn]7s^2 6d^0 5f^7        ',
     *'Cm: [Rn]7s^2 6d^0 5f^8        '/
      data ((iclocc(i,j),i=0,LMAX_occ),j=0,ichargemax) /
     & 0 , 0, 0, 0,       !0
     & 0 , 0, 0, 0,       !1
     & 1 , 0, 0, 0,       !2  
     & 1 , 0, 0, 0,       !3  
     & 2 , 0, 0, 0,       !4  
     & 2 , 0, 0, 0,       !5  
     & 2 , 0, 0, 0,       !6  
     & 2 , 0, 0, 0,       !7  
     & 2 , 0, 0, 0,       !8  
     & 2 , 0, 0, 0,       !9  
     & 2 , 1, 0, 0,       !10 
c 
     & 2 , 1, 0, 0,       !11 
     & 3 , 1, 0, 0,       !12 
     & 3 , 1, 0, 0,       !13 
     & 3 , 1, 0, 0,       !14 
     & 3 , 1, 0, 0,       !15 
     & 3 , 1, 0, 0,       !16 
     & 3 , 1, 0, 0,       !17 
     & 3 , 2, 0, 0,       !18 
     & 3 , 2, 0, 0,       !19 
     & 4 , 2, 0, 0,       !20 
c 
     & 4 , 2, 0, 0,       !21  
     & 4 , 2, 0, 0,       !22  
     & 4 , 2, 0, 0,       !23 
     & 4 , 2, 0, 0,       !24 
     & 4 , 2, 0, 0,       !25 
     & 4 , 2, 0, 0,       !26 
     & 4 , 2, 0, 0,       !27 
     & 4 , 2, 0, 0,       !28 
     & 3 , 2, 1, 0,       !29 
     & 4 , 2, 1, 0,       !30  
c 
     & 4 , 2, 1, 0,       !31 
     & 4 , 2, 1, 0,       !32 
     & 4 , 2, 1, 0,       !33  
     & 4 , 2, 1, 0,       !34 
     & 4 , 2, 1, 0,       !35 
     & 4 , 3, 1, 0,       !36 
     & 4 , 3, 1, 0,       !37 
     & 5 , 3, 1, 0,       !38 
     & 5 , 3, 1, 0,       !39 
     & 5 , 3, 1, 0,       !40 
c 
     & 5 , 3, 1, 0,       !41 
     & 5 , 3, 1, 0,       !42 
     & 5 , 3, 1, 0,       !43 
     & 5 , 3, 1, 0,       !44 
     & 5 , 3, 1, 0,       !45 
     & 5 , 3, 1, 0,       !46 
     & 4 , 3, 2, 0,       !47 
     & 5 , 3, 2, 0,       !48 
c 
     & 5 , 3, 2, 0,       !49  
     & 5 , 3, 2, 0,       !50  
     & 5 , 3, 2, 0,       !51  
     & 5 , 3, 2, 0,       !52  
     & 5 , 3, 2, 0,       !53  
     & 5 , 4, 2, 0,       !54  
     & 5 , 4, 2, 0,       !55  
     & 6 , 4, 2, 0,       !56 
     & 6 , 4, 2, 0,       !57 
     & 6 , 4, 2, 0,       !58 
     & 6 , 4, 2, 0,       !59 
     & 6 , 4, 2, 0,       !60 
c  
     & 6 , 4, 2, 0,       !61  
     & 6 , 4, 2, 0,       !62  
     & 6 , 4, 2, 0,       !63  
     & 6 , 4, 2, 0,       !64  
     & 6 , 4, 2, 0,       !65  
     & 6 , 4, 2, 0,       !66  
     & 6 , 4, 2, 0,       !67  
     & 6 , 4, 2, 0,       !68  
     & 6 , 4, 2, 0,       !69  
     & 6 , 4, 2, 1,       !70  
c   
     & 6 , 4, 2, 1,       !71 
     & 6 , 4, 2, 1,       !72 
     & 6 , 4, 2, 1,       !73 
     & 6 , 4, 2, 1,       !74 
     & 6 , 4, 2, 1,       !75 
     & 6 , 4, 2, 1,       !76 
     & 6 , 4, 2, 1,       !77 
     & 5 , 4, 2, 1,       !78 
     & 5 , 4, 3, 1,       !79 
     & 6 , 4, 3, 1,       !80 
c
     & 6 , 4, 3, 1,       !81  
     & 6 , 4, 3, 1,       !82  
     & 6 , 4, 3, 1,       !83  
     & 6 , 4, 3, 1,       !84  
     & 6 , 4, 3, 1,       !85  
     & 6 , 5, 3, 1,       !86  
     & 6 , 5, 3, 1,       !87  
     & 7 , 5, 3, 1,       !88 
     & 7 , 5, 3, 1,       !89 
     & 7 , 5, 3, 1,       !90 
c    
     & 7 , 5, 3, 1,       !91  
     & 7 , 5, 3, 1,       !92  
     & 7 , 5, 3, 1,       !93  
     & 7 , 5, 3, 1,       !94  
     & 7 , 5, 3, 1,       !95 
     & 7 , 5, 3, 1/       !96 
cbs
      data ((iopocc(i,j),i=0,LMAX_occ),j=0,ichargemax) /
     & 0 , 0, 0, 0,    !0 
c
     & 0 , 0, 0, 0,    ! 1  
     & 0 , 0, 0, 0,    ! 2 
     & 1 , 0, 0, 0,    ! 3 
     & 0 , 0, 0, 0,    ! 4  
     & 0 , 1, 0, 0,    ! 5 
     & 0 , 2, 0, 0,    ! 6 
     & 0 , 3, 0, 0,    ! 7  
     & 0 , 4, 0, 0,    ! 8 
     & 0 , 5, 0, 0,    ! 9 
     & 0 , 0, 0, 0,    ! 10 
c 
     & 1 , 0, 0, 0,    ! 11  
     & 0 , 0, 0, 0,    ! 12   
     & 0 , 1, 0, 0,    ! 13  
     & 0 , 2, 0, 0,    ! 14   
     & 0 , 3, 0, 0,    ! 15   
     & 0 , 4, 0, 0,    ! 16   
     & 0 , 5, 0, 0,    ! 17  
     & 0 , 0, 0, 0,    ! 18   
     & 1 , 0, 0, 0,    ! 19   
     & 0 , 0, 0, 0,    ! 20   
c
     & 0 , 0, 1, 0,    ! 21  
     & 0 , 0, 2, 0,    ! 22   
     & 0 , 0, 3, 0,    ! 23  
     & 0 , 0, 4, 0,    ! 24   
     & 0 , 0, 5, 0,    ! 25   
     & 0 , 0, 6, 0,    ! 26   
     & 0 , 0, 7, 0,    ! 27  
     & 0 , 0, 8, 0,    ! 28  
     & 1 , 0, 0, 0,    ! 29   
     & 0 , 0, 0, 0,    ! 30  
c
     & 0 , 1, 0, 0,    ! 31  
     & 0 , 2, 0, 0,    ! 32  
     & 0 , 3, 0, 0,    ! 33  
     & 0 , 4, 0, 0,    ! 34   
     & 0 , 5, 0, 0,    ! 35  
     & 0 , 0, 0, 0,    ! 36   
     & 1 , 0, 0, 0,    ! 37   
     & 0 , 0, 0, 0,    ! 38  
     & 0 , 0, 1, 0,    ! 39  
     & 0 , 0, 2, 0,    ! 40  
c
     & 0 , 0, 3, 0,    ! 41 
     & 0 , 0, 4, 0,    ! 42  
     & 0 , 0, 5, 0,    ! 43  
     & 0 , 0, 6, 0,    ! 44  
     & 0 , 0, 7, 0,    ! 45   
     & 0 , 0, 8, 0,    ! 46   
     & 1 , 0, 0, 0,    ! 47  
     & 0 , 0, 0, 0,    ! 48  
     & 0 , 1, 0, 0,    ! 49  
     & 0 , 2, 0, 0,    ! 50  
c
     & 0 , 3, 0, 0,    ! 51  
     & 0 , 4, 0, 0,    ! 52  
     & 0 , 5, 0, 0,    ! 53  
     & 0 , 0, 0, 0,    ! 54  
     & 1 , 0, 0, 0,    ! 55  
     & 0 , 0, 0, 0,    ! 56     
     & 0 , 0, 1, 0,    ! 57  
     & 0 , 0, 0, 2,    ! 58  
     & 0 , 0, 0, 3,    ! 59  
     & 0 , 0, 0, 4,    ! 60  
c
     & 0 , 0, 0, 5,    ! 61  
     & 0 , 0, 0, 6,    ! 62  
     & 0 , 0, 0, 7,    ! 63  
     & 0 , 0, 0, 8,    ! 64  
     & 0 , 0, 0, 9,    ! 65  
     & 0 , 0, 0, 10,    ! 66   
     & 0 , 0, 0, 11,    ! 67  
     & 0 , 0, 0, 12,    ! 68  
     & 0 , 0, 0, 13,    ! 69  
     & 0 , 0, 0,  0,    ! 70  
c
     & 0 , 0, 1, 0,    ! 71  
     & 0 , 0, 2, 0,    ! 72  
     & 0 , 0, 3, 0,    ! 73  
     & 0 , 0, 4, 0,    ! 74  
     & 0 , 0, 5, 0,    ! 75  
     & 0 , 0, 6, 0,    ! 76  
     & 0 , 0, 7, 0,    ! 77   
     & 1 , 0, 9, 0,    ! 78  
     & 1 , 0, 0, 0,    ! 79   
     & 0 , 0, 0, 0,    ! 80   
c
     & 0 , 1, 0, 0,    ! 81  
     & 0 , 2, 0, 0,    ! 82  
     & 0 , 3, 0, 0,    ! 83  
     & 0 , 4, 0, 0,    ! 84  
     & 0 , 5, 0, 0,    ! 85  
     & 0 , 0, 0, 0,    ! 86  
     & 1 , 0, 0, 0,    ! 87  
     & 0 , 0, 0, 0,    ! 88   
     & 0 , 0, 1, 0,    ! 89   
     & 0 , 0, 2, 0,    ! 90  
c
     & 0 , 0, 1, 2,    ! 91  
     & 0 , 0, 1, 3,    ! 92  
     & 0 , 0, 1, 4,    ! 93  
     & 0 , 0, 0, 6,    ! 94  
     & 0 , 0, 0, 7,    ! 95  
     & 0 , 0, 0, 8/    ! 96   
cbs
      if (icharge.gt.ichargemax) then
         CALL QUIT('occupations not implemented')
      endif
      write(LUPRI,'(A32,A30)') txt,occtxt(icharge)
      do irun=0,min(lmax,lmax_occ)
         iclosed(irun)=iclocc(irun,icharge)
         iopen(irun)=iopocc(irun,icharge)
      end do
      do irun=min(lmax,lmax_occ)+1,lmax
         iclosed(irun)=0
         iopen(irun)=0
      end do
      return
      end
      double precision function  getCG(
     *j1, j2, j3, m1, m2, m3)
c    *j1,     ! integer  2*j1
c    *j2,     ! integer  2*j2
c    *j3,     ! integer  2*j3
c    *m1,     ! integer  2*m1
c    *m2,     ! integer  2*m2
c    *m3)     ! integer  2*m2
cbs this routine calculates the Clebsch-Gordon-coefficients
cbs by actually calculating the 3j-symbol 
cbs  ---                 ---
cbs  |  j1   j2    |   j3   |         j1+m1+j2-m2     
cbs  |             |        |  =  (-)                 sqrt (2  j3+1) *
cbs  |  m1   m2    |   m3   |
cbs  ---                 ---
cbs
cbs                             ---             ---
cbs                             |  j1   j2   j3   |     
cbs                             |                 |     
cbs                             |  m1   m2  -m3   |
cbs                              ---            ---
#include "implicit.h"
cbs   initialize CG-coefficient
      getCG=0d0
cbs   quick check 
      if (m1+m2.ne.m3) return 
      if (j1.lt.0.or.j2.lt.0.or.j3.lt.0) return   
cbs   check the correct sign    beginning  
      idummy=(j1+j2+m1-m2)/2
      if (mod(idummy,2).eq.0) then 
      isign=1
      else
      isign=-1
      endif 
cbs   check the correct sign    end          
      fac1=dsqrt(dfloat(j3+1))
      fac2=regge3j(j1,j2,j3,m1,m2,-m3)
      getCG=isign*fac1*fac2
      return 
      end 


      
      





      subroutine getLIMIT(l1,l2,l3,l4,Lanf,Lend)
#include "implicit.h"
#include "priunit.h"
cbs   get the minimum and maximum L-values 
cbs   of the the coulomb-potential to interact 
cbs   with l1-l4
      lower1=iabs(l1-l3) 
      lower2=iabs(l2-l4)         
      lupper1=l1+l3 
      lupper2=l2+l4
      Lanf=max(lower1,lower2)
      Lend=min(lupper1,lupper2) 
cbs     check for parity   
      lsum=Lanf+l1+l3
      if (mod(lsum,2).eq.1) Lanf=Lanf+1
      lsum=Lend+l1+l3
      if (mod(lsum,2).eq.1) Lend=Lend-1   
cbs   check the other parity 
      lsum=Lanf+l2+l4
      if (mod(lsum,2).eq.1) then 
      write(LUPRI,*) ' error in getLIMIT: '
      write(LUPRI,*) ' parity inconsistency for '
      write(LUPRI,*) 'l1,l2,l3,l4= ',l1,l2,l3,l4
      CALL QUIT('Inconsistency error in getLIMIT')              
      endif  
      return 
      end 
      subroutine getpow(max,quot,quotpow,
     *nprim1,nprim2,nprim3,nprim4)
cbs   generates some powers of for the prefactors of cfunct(X)
cbs   look out for details there and in initfrac   
#include "implicit.h"
#include "para.h"
      dimension quotpow(nprim1,nprim2,
     *nprim3,nprim4),
     *quot(nprim1,nprim2,nprim3,nprim4) 
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=1,nprim1
      quotpow(irun1,irun2,irun3,irun4)=
     *dsqrt(quot(irun1,irun2,irun3,irun4))
      enddo
      enddo
      enddo
      enddo
      if (max.eq.1) return 
cbs  
      do irun=2,max     
      do irun4=1,nprim4
      do irun3=1,nprim3
      do irun2=1,nprim2
      do irun1=1,nprim1
      quotpow(irun1,irun2,irun3,irun4)=
     *quotpow(irun1,irun2,irun3,irun4)*
     *quot(irun1,irun2,irun3,irun4)
      enddo
      enddo
      enddo
      enddo
      enddo
      return 
      end 
      subroutine inidf
cbs   initializes the df on common block  with double facultatives
#include "implicit.h"
#include "para.h"
#include "amfi_param.h"
#include "dofuc.h"
      df(0)=1.d0
      df(1)=1.d0
      do irun=2,ndfmx
      df(irun)=dfloat(irun)*df(irun-2)
      enddo
      do jbm=0,ndfmx-1  
      do ibm=jbm,ndfmx   
      dffrac(ibm,jbm)=df(ibm)/df(jbm) 
      enddo 
      enddo 
      do jbm=1,ndfmx  
      do ibm=0,jbm-1       
      dffrac(ibm,jbm)=1d0/dffrac(jbm,ibm)   
      enddo 
      enddo 
      return 
      end
      subroutine initfrac(nprimit1,nprimit2,
     *nprimit3,nprimit4,
     *quot1,quot2,expo1,expo2,
     *expo3,expo4)
cbs   initialize some arrays with factors  needed for cfunct(x) 
#include "implicit.h"
      dimension expo1(*),expo2(*),expo3(*),expo4(*), 
     *quot1(nprimit1,nprimit2,nprimit3,nprimit4),
     *quot2(nprimit1,nprimit2,nprimit3,nprimit4)
      do irun4=1,nprimit4 
      do irun3=1,nprimit3 
      do irun2=1,nprimit2 
        sum24=expo2(irun2)+expo4(irun4)
                do irun1=1,nprimit1 
                quot1(irun1,irun2,irun3,irun4)= 
     *          1d0/(1d0+(expo1(irun1)+expo3(irun3))/
     *          sum24)                             
                enddo
      enddo
      enddo
      enddo
      do irun4=1,nprimit4 
      do irun3=1,nprimit3 
      do irun2=1,nprimit2 
        sum24=expo2(irun2)+expo4(irun4)
                do irun1=1,nprimit1 
                quot2(irun1,irun2,irun3,irun4)= 
     *          1d0/(1d0+sum24/
     *          (expo1(irun1)+expo3(irun3)))
                enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine initired
#include "implicit.h"
cbs   initialize all information for ireducible representations 
cbs   later on, it might be useful to have a switch for 
cbs    changing to other orders of IREDs like e.g. in TURBOMOLE
c
c
c   HOW2ADD another symmetry:
c
c   1. add it in readbas.f to be accepted. Add the number of IRs
c
c   2. copy one of the symmetry-blocks in this subroutine and 
c      edit the multiplication-table for the group
c
c   3. assign the right IRs to L_X, L_Y and L_Z 
c
c   that is  all. Good luck!!!
c
#include "priunit.h"
#include "para.h"
#include "ired.h"
      character*3 symmetry 
      symmetry='D2H'  ! MOLCAS-Version
      if (symmetry.eq.'D2H') then 
      mult(2,1)=2
      mult(3,1)=3
      mult(4,1)=4
      mult(5,1)=5
      mult(6,1)=6
      mult(7,1)=7
      mult(8,1)=8
c  
      mult(3,2)=4
      mult(4,2)=3
      mult(5,2)=6
      mult(6,2)=5
      mult(7,2)=8
      mult(8,2)=7
c  
      mult(4,3)=2
      mult(5,3)=7
      mult(6,3)=8
      mult(7,3)=5
      mult(8,3)=6
c  
      mult(5,4)=8
      mult(6,4)=7
      mult(7,4)=6
      mult(8,4)=5
c  
      mult(6,5)=2
      mult(7,5)=3
      mult(8,5)=4
c  
      mult(7,6)=4
      mult(8,6)=3
c  
      mult(8,7)=2
c  
C      
      do ired=1,8
      mult(ired,ired)=1
      enddo
      do irun=2,8
      do jrun=1,irun-1
      mult(jrun,irun)=mult(irun,jrun)
      enddo
      enddo
CBS   write(6,*) 
CBS   write(6,*) 
CBS  *'multiplicitation table (atkins,child and phillips)'
CBS   write(6,*) 
CBS   do ired=1,8
CBS   write(6,'(8I5)') (mult(jred,ired),jred=1,8) 
CBS   write(6,*) 
CBS   enddo
      
c     
      IRLX=4
      IRLY=3
      IRLZ=2
cbs   assume same order of ireds as Atkins Child and Phillips use..
cbs   would lead to an order with 1 to 1, 2 to 2 ...
cbs   however, this is the molecule/ seward order.   
      iredorder(1)=1
      iredorder(2)=4
      iredorder(3)=6
      iredorder(4)=7
      iredorder(5)=8
      iredorder(6)=5
      iredorder(7)=3
      iredorder(8)=2
      do ired=1,8
      iredorderinv(iredorder(ired))=ired
      enddo
      ipow2ired(0,0,0)=iredorder(1)
      ipow2ired(1,1,0)=iredorder(2)
      ipow2ired(1,0,1)=iredorder(3)
      ipow2ired(0,1,1)=iredorder(4)
      ipow2ired(1,1,1)=iredorder(5)
      ipow2ired(0,0,1)=iredorder(6)
      ipow2ired(0,1,0)=iredorder(7)
      ipow2ired(1,0,0)=iredorder(8)
c     write(6,*) 'interacting IRs '
      do ired=1,8
      IRwithLX(ired)=
     *iredorder(mult(IRLX,iredorderinv(ired)))
      IRwithLY(ired)=
     *iredorder(mult(IRLY,iredorderinv(ired)))
      IRwithLZ(ired)=
     *iredorder(mult(IRLZ,iredorderinv(ired)))
c     write(6,*) IRwithLX(ired),IRwithLY(ired),
c    *IRwithLZ(ired)
      enddo
      elseif(symmetry.eq.'C2V') then 
cbs   1. A1 2. A2 3. B1 4. B2
      mult(2,1)=2
      mult(3,1)=3
      mult(4,1)=4
c  
      mult(3,2)=4
      mult(4,2)=3
c  
      mult(4,3)=2
C      
      do ired=1,4
      mult(ired,ired)=1
      enddo
      do irun=2,4
      do jrun=1,irun-1
      mult(jrun,irun)=mult(irun,jrun)
      enddo
      enddo
      write(LUPRI,*) 
      write(LUPRI,*) 
     *'multiplicitation table '
      write(LUPRI,*) 
      do ired=1,4
      write(LUPRI,'(4I5)') (mult(jred,ired),jred=1,4) 
      write(LUPRI,*) 
      enddo
      
c     
      IRLX=4
      IRLY=3
      IRLZ=2
cbs   this is the molecule/ seward order.   
      iredorder(1)=1
      iredorder(2)=4
      iredorder(3)=2
      iredorder(4)=3
      do ired=1,4
      iredorderinv(iredorder(ired))=ired
      enddo
      ipow2ired(0,0,0)=iredorder(1)
      ipow2ired(1,1,0)=iredorder(2)
      ipow2ired(1,0,1)=iredorder(3)
      ipow2ired(0,1,1)=iredorder(4)
      ipow2ired(1,1,1)=iredorder(2)
      ipow2ired(0,0,1)=iredorder(1)
      ipow2ired(0,1,0)=iredorder(4)
      ipow2ired(1,0,0)=iredorder(3)
c     write(6,*) 'interacting IRs '
      do ired=1,4
      IRwithLX(ired)=
     *iredorder(mult(IRLX,iredorderinv(ired)))
      IRwithLY(ired)=
     *iredorder(mult(IRLY,iredorderinv(ired)))
      IRwithLZ(ired)=
     *iredorder(mult(IRLZ,iredorderinv(ired)))
c     write(6,*) IRwithLX(ired),IRwithLY(ired),
c    *IRwithLZ(ired)
      enddo
      elseif(symmetry.eq.'D2 ') then 
cbs   1. A1 2. B1 3. B2 4. B3
      mult(2,1)=2
      mult(3,1)=3
      mult(4,1)=4
c  
      mult(3,2)=4
      mult(4,2)=3
      mult(4,3)=2
C      
      do ired=1,4
      mult(ired,ired)=1
      enddo
      do irun=2,4
      do jrun=1,irun-1
      mult(jrun,irun)=mult(irun,jrun)
      enddo
      enddo
      write(LUPRI,*) 
      write(LUPRI,*) 
     *'multiplicitation table '
      write(LUPRI,*) 
      do ired=1,4
      write(LUPRI,'(4I5)') (mult(jred,ired),jred=1,4) 
      write(LUPRI,*) 
      enddo
      
c     
      IRLX=4
      IRLY=3
      IRLZ=2
      iredorder(1)=1
      iredorder(2)=2
      iredorder(3)=3
      iredorder(4)=4
      do ired=1,4
      iredorderinv(iredorder(ired))=ired
      enddo
      ipow2ired(0,0,0)=iredorder(1)
      ipow2ired(1,1,0)=iredorder(2)
      ipow2ired(1,0,1)=iredorder(3)
      ipow2ired(0,1,1)=iredorder(4)
      ipow2ired(1,1,1)=iredorder(1)
      ipow2ired(0,0,1)=iredorder(2)
      ipow2ired(0,1,0)=iredorder(3)
      ipow2ired(1,0,0)=iredorder(4)
c     write(6,*) 'interacting IRs '
      do ired=1,4
      IRwithLX(ired)=
     *iredorder(mult(IRLX,iredorderinv(ired)))
      IRwithLY(ired)=
     *iredorder(mult(IRLY,iredorderinv(ired)))
      IRwithLZ(ired)=
     *iredorder(mult(IRLZ,iredorderinv(ired)))
c     write(6,*) IRwithLX(ired),IRwithLY(ired),
c    *IRwithLZ(ired)
      enddo
      elseif(symmetry.eq.'C2H') then 
cbs   assume 1.Ag 2.Au 3.Bg 4.Bu 
      mult(2,1)=2
      mult(3,1)=3
      mult(4,1)=4
c  
      mult(3,2)=4
      mult(4,2)=3
c  
      mult(4,3)=2
C      
      do ired=1,4
      mult(ired,ired)=1
      enddo
      do irun=2,4
      do jrun=1,irun-1
      mult(jrun,irun)=mult(irun,jrun)
      enddo
      enddo
      write(LUPRI,*) 
      write(LUPRI,*) 
     *'multiplicitation table '
      write(LUPRI,*) 
      do ired=1,4
      write(LUPRI,'(4I5)') (mult(jred,ired),jred=1,4) 
      write(LUPRI,*) 
      enddo
      
c     
      IRLX=3
      IRLY=3
      IRLZ=1
      iredorder(1)=1
      iredorder(2)=2
      iredorder(3)=3
      iredorder(4)=4
      do ired=1,4
      iredorderinv(iredorder(ired))=ired
      enddo
      ipow2ired(0,0,0)=iredorder(1)
      ipow2ired(1,1,0)=iredorder(1)
      ipow2ired(1,0,1)=iredorder(3)
      ipow2ired(0,1,1)=iredorder(3)
      ipow2ired(1,1,1)=iredorder(2)
      ipow2ired(0,0,1)=iredorder(2)
      ipow2ired(0,1,0)=iredorder(4)
      ipow2ired(1,0,0)=iredorder(4)
c     write(6,*) 'interacting IRs '
      do ired=1,4
      IRwithLX(ired)=
     *iredorder(mult(IRLX,iredorderinv(ired)))
      IRwithLY(ired)=
     *iredorder(mult(IRLY,iredorderinv(ired)))
      IRwithLZ(ired)=
     *iredorder(mult(IRLZ,iredorderinv(ired)))
c     write(6,*) IRwithLX(ired),IRwithLY(ired),
c    *IRwithLZ(ired)
      enddo
      elseif(symmetry.eq.'CS ') then 
      write(LUPRI,*) 'CS in initired '
cbs   assume 1.A' 2.A'                       
      mult(2,1)=2
C      
      do ired=1,2
      mult(ired,ired)=1
      enddo
      do irun=2,2
      do jrun=1,irun-1
      mult(jrun,irun)=mult(irun,jrun)
      enddo
      enddo
      write(LUPRI,*) 
      write(LUPRI,*) 
     *'multiplicitation table '
      write(LUPRI,*) 
      do ired=1,2
      write(LUPRI,'(2I5)') (mult(jred,ired),jred=1,2) 
      write(LUPRI,*) 
      enddo
      
c     
      IRLX=2
      IRLY=2
      IRLZ=1
      iredorder(1)=1
      iredorder(2)=2
      do ired=1,2
      iredorderinv(iredorder(ired))=ired
      enddo
      ipow2ired(0,0,0)=iredorder(1)
      ipow2ired(1,1,0)=iredorder(1)
      ipow2ired(1,0,1)=iredorder(2)
      ipow2ired(0,1,1)=iredorder(2)
      ipow2ired(1,1,1)=iredorder(2)
      ipow2ired(0,0,1)=iredorder(2)
      ipow2ired(0,1,0)=iredorder(1)
      ipow2ired(1,0,0)=iredorder(1)
c     write(6,*) 'interacting IRs '
      do ired=1,2
      IRwithLX(ired)=
     *iredorder(mult(IRLX,iredorderinv(ired)))
      IRwithLY(ired)=
     *iredorder(mult(IRLY,iredorderinv(ired)))
      IRwithLZ(ired)=
     *iredorder(mult(IRLZ,iredorderinv(ired)))
c     write(6,*) IRwithLX(ired),IRwithLY(ired),
c    *IRwithLZ(ired)
      enddo
      endif   
      return 
      end 
      subroutine kindiag(TKIN,TKINTRIA,ndim,evec,eval,breit)
#include "implicit.h"
cbs   determines eigenvectors and -values of TKIN  
      dimension tkin(ndim,ndim),
     *TKINTRIA((ndim*ndim+ndim)/2),eval(ndim),evec(ndim,ndim)
      logical breit
cbs   move symmetric matrix to triangular matrix 
      itria=1
      do irun2=1,ndim
      do irun1=1,irun2 
      TKINTRIA(itria)=TKIN(irun1,irun2)
      itria=itria+1 
      enddo
      enddo
      do irun2=1,ndim
      do irun1=1,ndim
      evec(irun1,irun2)=0d0
      enddo
      enddo
      do irun1=1,ndim
      evec(irun1,irun1)=1d0
      enddo
cbs   now diagonalize  
            CALL jacobi(TKINTRIA,evec,ndim,ndim)    
cbs   get the eigenvalues   
      do irun=1,ndim
      eval(irun)=TKINTRIA((irun*irun+irun)/2)
      enddo
      if (breit) then
      do irun=1,ndim
      eval(irun)=0d0 
      enddo
      endif 
cbs   ensure normalization of the vectors. 
      do IRUN=1,ndim
      fact=0d0
      do JRUN=1,ndim 
      fact=fact+evec(JRUN,IRUN)*evec(JRUN,IRUN) 
      enddo
      fact=1d0/dsqrt(fact)
      do JRUN=1,ndim
      evec(JRUN,IRUN)=fact*evec(JRUN,IRUN)
      enddo 
      enddo
      return   
      end   
      Subroutine kinemat(L,ndim,evtkin,type1,type2,Energy)
#include "implicit.h"
#include "codata.h"
cbs   at least it's identical with Odd's valuE
      parameter (speed2=CVEL*CVEL) 
      parameter (speed4=speed2*speed2) 
cbs   this routine generates the kinematic A-factors=dsqrt((E+mc^2)/(2E))  
cbs   (type1) and   c*A/(E+mc^2) (type2)
cbs   The c in the second kinematic factor comes from Jan Almloef and 
cbs   Odd Gropen in Rev in Comp.Chem. 8(1996)
      dimension evtkin(*),type1(*),type2(*),Energy(*)  
c     E= dsqrt(p**2 c**2 + m**2 c**4) 
c     p**2= 2*m*TKIN    
c     with m = 1 
      do Irun=1,ndim
      if (evtkin(Irun).lt.0) CALL QUIT('strange kinetic energy ')
      Energy(Irun)=(evtkin(Irun)+evtkin(Irun))*speed2+speed4 
      enddo
      do Irun=1,ndim
      Energy(Irun)=dsqrt(energy(irun))
      enddo
      do Irun=1,ndim
!     dsqrt((E+mc^2)/(2E)):
      type1(Irun)=dsqrt(0.5d0*(1d0+speed2/Energy(Irun)))
      enddo
!      c*A/(E+mc^2) 
      do Irun=1,ndim
      type2(Irun)=CVEL*type1(Irun)/(Energy(Irun)+speed2)
      enddo
              do Irun=1,ndim
              type2(Irun)=2*CVEL*type2(Irun)
              enddo                        
      return 
      end
      Double precision function LMdepang(
     *L,M,l1,l2,l3,l4,m1,m2,m3,m4,cheater)
cbs   l1-l4 and m1-m4 are already shifted !!
cbs   purpose: calculates the angular part of the   
cbs   coulomb-type integrals. See documentation for details...
cbs   LMdepang= LM dependent angular factors 
cbs   cheater included for a correcting signs, as there were some 
cbs   signs (only signs!!!!) missing when compared to HERMIT  
cbs                                        B.S.  08.10.96 
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
      LMdepang=0d0
cbs   some quick checks
      if (L.lt.abs(M)) return 
      if (l1.lt.abs(m1)) return 
      if (l2.lt.abs(m2)) return 
      if (l3.lt.abs(m3)) return 
      if (l4.lt.abs(m4)) return 
cbs   prefactor
      fact1=4d0*pi/dfloat(L+L+1)
cbs   determining the sign
      isum=-l3-l1-l4-l2+2*(M+m3+m4)   !???? I am not sure 
      if (mod(isum,4).eq.0) then 
      isign=1
      elseif (iabs(mod(isum,4)).eq.2) then 
      isign=-1
      else 
      write(LUPRI,*) 'L,l1,l2,l3,l4,M,m1,m2,m3,m4'
      write(LUPRI,'(10I3)') L,l1,l2,l3,l4,M,m1,m2,m3,m4
      write(LUPRI,*) 'isum= ',isum,' mod = ',mod(isum,4)
      CALL QUIT('error in lmdepang')
      endif
      fact2=couple3J(L,l3,l1,-M,m3,-m1) 
      fact3=couple3J(L,l4,l2,M,m4,-m2)
C     write(6,*) 'fact2,fact3 ',fact2,fact3
      LMdepang=cheater*dfloat(isign)*fact1*fact2*fact3
      return 
      end 
      logical function mcheckxy(m1,m2,m3,m4)
      integer m1,m2,m3,m4,int12a,int12b,
     *int34a,int34b
cbs   makes a check, if there is an interaction inbetween cartesian functions 
cbs   with m-values m1-m4 
      mcheckxy=.true.   
      int12a=m1+m2
      int12b=-m1+m2
      int34a=m3+m4
      int34b=-m3+m4
cbs   lots of checks 
      if (iabs(int12a+int34a).eq.1) return
      if (iabs(int12a-int34a).eq.1) return
      if (iabs(int12b+int34b).eq.1) return
      if (iabs(int12b-int34b).eq.1) return
      if (iabs(int12a+int34b).eq.1) return
      if (iabs(int12a-int34b).eq.1) return
      if (iabs(int12b+int34a).eq.1) return
      if (iabs(int12b-int34a).eq.1) return
      mcheckxy=.false.
      return 
      end 
      logical function mcheckz(m1,m2,m3,m4)
cbs   makes a check, if there is an interaction inbetween cartesian functions 
cbs   with m-values m1-m4 
      integer m1,m2,m3,m4,int12a,int12b,
     *int34a,int34b
      mcheckz=.true.   
      int12a=m1+m2
      int12b=-m1+m2
      int34a=m3+m4
      int34b=-m3+m4
cbs   lots of checks 
      if (iabs(int12a+int34a).eq.0) return
      if (iabs(int12a-int34a).eq.0) return
      if (iabs(int12b+int34b).eq.0) return
      if (iabs(int12b-int34b).eq.0) return
      if (iabs(int12a+int34b).eq.0) return
      if (iabs(int12a-int34b).eq.0) return
      if (iabs(int12b+int34a).eq.0) return
      if (iabs(int12b-int34a).eq.0) return
      mcheckz=.false.
      return 
      end 
      subroutine mkangL0(Lmax,l1,l2,l3,l4,m1,m2,m3,m4,
     *angintSO,angintOO,
     *Lfirst,Llast,Lblocks,
     *ncont1,ncont2,ncont3,
     *ncont4,
     *caseaSO,caseb1SO,caseb2SO,casecSO,
     *caseaOO,caseb1OO,caseb2OO,casecOO,
     *preroots,clebsch,dummy,bonn,breit,
     *sameorb)
#include "implicit.h"
cbs   subroutine for combining radial integrals with angular 
cbs   factors for the block with l1,l2,l3,l4,m1,m2,m3m,m4  
cbs   this routine mkangL0 = make angular factors for the L0-part 
cbs   includes both, spin-same and spin-other-orbit parts. 
      double precision LMdepang
      dimension 
     *angintSO(ncont1,ncont2,ncont3,ncont4),
     *angintOO(ncont1,ncont2,ncont3,ncont4),
     *Lfirst(*),Llast(*),Lblocks(*),
cbs   all the arrays with the radial integrals for 
cbs   this combination of l-values   
     *caseaSO(ncont1*ncont2*ncont3*ncont4,*),  ! (2,0)   integrals with alpha1*alpha3
     *caseb1SO(ncont1*ncont2*ncont3*ncont4,*), ! (0,0)   integrals with alpha1
     *caseb2SO(ncont1*ncont2*ncont3*ncont4,*), ! (0,0)   integrals with alpha3
     *casecSO(ncont1*ncont2*ncont3*ncont4,*),  ! (-2,0)  integrals with factor 1
     *caseaOO(ncont1*ncont2*ncont3*ncont4,*),  ! (2,0)   integrals with alpha1*alpha3
     *caseb1OO(ncont1*ncont2*ncont3*ncont4,*), ! (0,0)   integrals with alpha1
     *caseb2OO(ncont1*ncont2*ncont3*ncont4,*), ! (0,0)   integrals with alpha3
     *casecOO(ncont1*ncont2*ncont3*ncont4,*),  ! (-2,0)  integrals with factor 1
     *preroots(2,0:Lmax),                    ! some prefactors: dsqrt( (l(+1))/(2l+1))
     *clebsch(3,2,-Lmax:Lmax,0:Lmax)         ! some clebsch gordans, that appear regulary
      dimension dummy(0:*)                                    
      logical bonn,breiT,sameorb 
c     write(6,*) 'begin mkangL0 ',
c    *l1,l2,l3,l4,m1,m2,m3,m4
cbs  
      ncontall=ncont1*ncont2*ncont3*ncont4
cbs   cheater introduced to correct signs, because they were different from HERMIT 
      if (mod(l1+l2+l3+l4,4).eq.2) then                  
      cheater=1d0
      else 
      cheater=-1d0
      endif 
cbs   cleaning up 
      if (bonn.or.breit.or.sameorb) then 
      call dzero(angintSO,ncontall) 
      else 
      call dzero(angintSO,ncontall) 
      call dzero(angintOO,ncontall) 
      endif  
cbs  starting with the same-orbit-contributions 
cbs  first term: ###########################################################################
      factor=-preroots(2,l1)*preroots(2,l3)*
     *clebsch(1,2,m1,l1)*
     *clebsch(1,2,m3,l3)
      if (factor.ne.0d0) then 
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
      M=m2-m4
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1+1,l2,l3+1,l4,m1-1,m2,m3-1,m4,cheater) 
      if (dummy(L).ne.0d0) then
      if (bonn.or.breit.or.sameorb) then 
         Call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),CaseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs   second term: ###########################################################################
      factor=-preroots(1,l1)*preroots(2,l3)*
     *clebsch(1,1,m1,l1)*
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
      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1-1,m2,m3-1,m4,cheater)
      if (dummy(L).ne.0d0) then
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),caseaSO(1,Lrun),1,
     *   angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),caseaSO(1,Lrun),1,
     *   angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,AngintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(3).gt.0) then 
      M=m2-m4
        if (Lfirst(3).lt.Kfirst) then 
        do L=Lfirst(3),Kfirst,2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
        enddo 
        Kfirst=Lfirst(3)
        endif 
        if (Llast(3).gt.Klast) then 
        do L=Klast,Llast(3),2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
        enddo 
        Klast=Llast(3)
        endif 
      Lrun=1
      do L=Lfirst(3),Llast(3),2  
      if (dummy(L).ne.0d0) then
      If (bonn.or.breit.or.sameorb) then 
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
      factor=-preroots(2,l1)*preroots(1,l3)*
     *clebsch(1,2,m1,l1)*
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
      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1-1,m2,
     *m3-1,m4,cheater)
      if (dummy(L).ne.0d0) then 
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),CaseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(2).gt.0) then 
      M=m2-m4
        if (Lfirst(2).lt.Kfirst) then 
        do L=Lfirst(2),Kfirst,2
        dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
        enddo 
        Kfirst=Lfirst(2)
        endif 
        if (Llast(2).gt.Klast) then 
        do L=Klast,Llast(2),2
        dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
        enddo 
        Klast=Llast(2)
        endif 
      Lrun=1
      do L=Lfirst(2),Llast(2),2  
      if (dummy(L).ne.0d0) then 
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   -(2+4*l3)*factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs   fourth term: ###########################################################################
      factor=-preroots(1,l1)*preroots(1,l3)*
     *clebsch(1,1,m1,l1)*
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
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1-1,m2,m3-1,m4,cheater)
      if (dummy(L).ne.0d0) then 
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),caseaSO(1,Lrun),1,
     *   angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),caseaSO(1,Lrun),1,
     *   angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,AngintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(2).gt.0) then 
      M=m2-m4
        if (Lfirst(2).lt.Kfirst) then 
        do L=Lfirst(2),Kfirst,2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
        enddo 
        Kfirst=Lfirst(2)
        endif 
        if (Llast(2).gt.Klast) then 
        do L=Klast,Llast(2),2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
        enddo 
        Klast=Llast(2)
        endif 
      Lrun=1
      do L=Lfirst(2),Llast(2),2  
      if (dummy(L).ne.0d0)  then 
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   -(2+4*l3)*factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(3).gt.0) then 
      M=m2-m4
        if (Lfirst(3).lt.Kfirst) then 
        do L=Lfirst(3),Kfirst,2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
        enddo 
        Kfirst=Lfirst(3)
        endif 
        if (Llast(3).gt.Klast) then 
        do L=Klast,Llast(3),2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
        enddo 
        Klast=Llast(3)
        endif 
      Lrun=1
      do L=Lfirst(3),Llast(3),2  
      if (dummy(L).ne.0d0)  then
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l1)*factor*dummy(L),
     *   caseb2SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   -(2+4*l1)*factor*dummy(L),caseb2OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(4).gt.0) then 
      M=m2-m4
        if (Lfirst(4).lt.Kfirst) then 
        do L=Lfirst(4),Kfirst,2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
        enddo 
        Kfirst=Lfirst(4)
        endif 
        if (Llast(4).gt.Klast) then 
        do L=Klast,Llast(4),2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1-1,m2,
     *        m3-1,m4,cheater)
        enddo 
        Klast=Llast(4)   
        endif 
      Lrun=1
      do L=Lfirst(4),Llast(4),2  
      if (dummy(L).ne.0d0)  then 
      If (bonn.or.breit.or.sameorb) then 
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
      factor=preroots(2,l1)*preroots(2,l3)*
     *clebsch(3,2,m1,l1)*
     *clebsch(3,2,m3,l3)
      if (factor.ne.0d0) then 
      do I=0,Lmax+Lmax+1
      dummy(I)=0d0
      enddo 
cbs   get the L,M dependent coefficients 
      if (Lblocks(1).gt.0) then 
      M=m2-m4
      Lrun=1
      do L=Lfirst(1),Llast(1),2  
      dummy(L)=LMdepang(L,M,l1+1,l2,l3+1,l4,m1+1,m2,m3+1,m4,cheater)
      if (dummy(L).ne.0d0) then 
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),caseaSO(1,Lrun),1,
     *   angintSO,1)
      else  
         call daxpy(ncontall,4*factor*dummy(L),caseaSO(1,Lrun),1,
     *   angintSO,1)
         call daxpy(ncontall,
     *   4*factor*dummy(L),caseaOO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo
      endif 
      endif 
cbs   sixth  term: ###########################################################################
      factor=preroots(1,l1)*preroots(2,l3)*
     *clebsch(3,1,m1,l1)*
     *clebsch(3,2,m3,l3)
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
      dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1+1,m2,m3+1,m4,cheater)
      if (dummy(L).ne.0d0)  then 
      If (bonn.or.breit.or.sameorb) then 
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
        dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
        enddo 
        Kfirst=Lfirst(3)
        endif 
        if (Llast(3).gt.Klast) then 
        do L=Klast,Llast(3),2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3+1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
        enddo 
        Klast=Llast(3)
        endif 
      Lrun=1
      do L=Lfirst(3),Llast(3),2  
      if (dummy(L).ne.0d0) then 
      If (bonn.or.breit.or.sameorb) then 
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
      factor=preroots(2,l1)*preroots(1,l3)*
     *clebsch(3,2,m1,l1)*
     *clebsch(3,1,m3,l3)
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
      dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1+1,m2,m3+1,m4,cheater)
      if (dummy(L).ne.0d0) then 
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,4*factor*dummy(L),
     *   caseaSO(1,Lrun),1,angintSO,1)
         Call daxpy(ncontall,
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
     *        m3+1,m4,cheater)
        enddo 
        Kfirst=Lfirst(2)
        endif 
        if (Llast(2).gt.Klast) then 
        do L=Klast,Llast(2),2
        dummy(L)=LMdepang(L,M,l1+1,l2,l3-1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
        enddo 
        Klast=Llast(2)
        endif 
      Lrun=1
      do L=Lfirst(2),Llast(2),2  
      if (dummy(L).ne.0d0)  then
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
      else 
         Call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
         Call daxpy(ncontall,-(2+4*l3)*
     *   factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      endif 
cbs   eigth term: ###########################################################################
      factor=preroots(1,l1)*preroots(1,l3)*
     *clebsch(3,1,m1,l1)*
     *clebsch(3,1,m3,l3)
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
      dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,m3+1,m4,cheater)
      if (dummy(L).ne.0d0) then
      If (bonn.or.breit.or.sameorb) then 
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
        dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
        enddo 
        Kfirst=Lfirst(2)
        endif 
        if (Llast(2).gt.Klast) then 
        do L=Klast,Llast(2),2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
        enddo 
        Klast=Llast(2)
        endif 
      Lrun=1
      do L=Lfirst(2),Llast(2),2  
      if (dummy(L).ne.0d0)  then
      If (bonn.or.breit.or.sameorb) then 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
      else 
         call daxpy(ncontall,-(2+4*l3)*factor*dummy(L),
     *   caseb1SO(1,Lrun),1,angintSO,1)
         call daxpy(ncontall,
     *   -(2+4*l3)*factor*dummy(L),caseb1OO(1,Lrun),1,angintOO,1)
      endif 
      endif 
      Lrun=Lrun+1
      enddo 
      endif 
      if (Lblocks(3).gt.0) then 
      M=m2-m4
        if (Lfirst(3).lt.Kfirst) then 
        do L=Lfirst(3),Kfirst,2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
        enddo 
        Kfirst=Lfirst(3)
        endif 
        if (Llast(3).gt.Klast) then 
        do L=Klast,Llast(3),2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
        enddo 
        Klast=Llast(3)
        endif 
      Lrun=1
      do L=Lfirst(3),Llast(3),2  
      if (dummy(L).ne.0d0)  then
      If (bonn.or.breit.or.sameorb) then 
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
        dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
        enddo 
        Kfirst=Lfirst(4)
        endif 
        if (Llast(4).gt.Klast) then 
        do L=Klast,Llast(4),2
        dummy(L)=LMdepang(L,M,l1-1,l2,l3-1,l4,m1+1,m2,
     *        m3+1,m4,cheater)
        enddo 
        Klast=Llast(4)   
        endif 
      Lrun=1
      do L=Lfirst(4),Llast(4),2  
      if (dummy(L).ne.0d0) then 
      If (bonn.or.breit.or.sameorb) then 
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
      return   
      end   
      subroutine mkangLmin(Lmax,l1,l2,l3,l4,m1,m2,m3,m4,             
     *angintSO,angintOO,
     *Lfirst,Llast,Lblocks,
     *ncont1,ncont2,ncont3,
     *ncont4,
     *caseaSO,caseb1SO,caseb2SO,casecSO,
     *caseaOO,caseb1OO,caseb2OO,casecOO,
     *preroots,clebsch,dummy,bonn,breit,
     *sameorb)
#include "implicit.h"
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
      logical bonn,breiT,sameorb 
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
      call dzero(angintSO,ncontall)
      else 
      call dzero(angintSO,ncontall)
      call dzero(angintOO,ncontall)
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
      return   
      end
   
      subroutine prefac(Lmax,preroots,clebsch)
#include "implicit.h"
      dimension preroots(2,0:Lmax),
     *clebsch(3,2,-Lmax:Lmax,0:Lmax)
cbs   the roots appearing in front of all 
cbs   the contributions 
c     write(6,*) 'begin of prefac'
      do L=0,Lmax
      fact=1d0/dsqrt(dfloat(L+L+1))
      preroots(1,L)=dsqrt(dfloat(L))*fact 
      preroots(2,L)=dsqrt(dfloat(L+1))*fact 
      enddo
cbs   there are Clebsch-Gordon-Coefficients 
cbs   which always appear:
cbs 
cbs   -----                       ------
cbs  |                                 |
cbs  |  l +/- 1     1        |      l  |
cbs  |                       |         |
cbs  |                       |         |   
cbs  |  m+/-1,0   -1,1,0     |      m  |   
cbs  |                       |         |   
cbs  |                                 |   
cbs   -----                       -----
cbs 
cbs 
cbs  array clebsch (3,2,-Lmax:Lmax,0:Lmax)
cbs  first index    1:  m-1 
cbs                 2:  m 
cbs                 3:  m+1
cbs  second index   1:  l-1
cbs                 2:  l+1
cbs  third index        m 
cbs  fourth index       l 
cbs  
c     write(6,*),'start to generate CGs' 
      do L=0,Lmax
      L2=L+L
      do M=-L,L
c     write(6,*) 'L,M: ',L,M 
      M2=M+M
cbs   getCG calculates CG-coeffecients. In order to avoid fractions,
cbs   e.g. for spins, arguments are doubled values...
      clebsch(1,1,M,L)=
     *getCG(L2-2,2,L2,M2-2,2,M2)
      clebsch(2,1,M,L)=
     *getCG(L2-2,2,L2,M2,0,M2)
      clebsch(3,1,M,L)=
     *getCG(L2-2,2,L2,M2+2,-2,M2)
      clebsch(1,2,M,L)=
     *getCG(L2+2,2,L2,M2-2,2,M2)
      clebsch(2,2,M,L)=
     *getCG(L2+2,2,L2,M2,0,M2)
      clebsch(3,2,M,L)=
     *getCG(L2+2,2,L2,M2+2,-2,M2)
      enddo
      enddo 
      return 
      end
     

      subroutine readbas(Lhigh,makemean,bonn,breit,
     *symmetry,sameorb,AIMP,oneonly,ncont4,numballcart,IN,ifinite,
     *EXP_FIN) 
cbs   suposed to read the maximum of l-values, the number of primitive and contracted 
cbs   functions, the exponents and contraction coefficients
#include "implicit.h"
#include "priunit.h"
#include "para.h"
#include "amfi_param.h"
#include "ired.h"
      character*4 WORD                         
      character*4 symmetry   
      character*13 Llimit 
      character*19 chcharge
      character*30 Nofprim
      character*28 addtext
      character*32 Nofcont
      character*76 Stars                  
      logical makemean,bonn,breit,
     *sameorb,AIMP,oneonly  
      common /nucleus/ charge,Exp_finite 
      Integer ibeginIRED(8),idelpersym(8) 
      dimension INOFT(Mxcart),INOFF(MxCart)
      stars='********************************************************'//
     * '********************'
      Llimit='MAX. L-VALUE:'
      chcharge=' CHARGE OF NUCLEUS:'
      Nofprim='NUMBER OF PRIMITIVE FUNCTIONS:' 
      Nofcont=' NUMBER OF CONTRACTED FUNCTIONS:' 
      addtext='ADDITIONAL FUNCTIONS in IRS:'
CBS   write(LUPRI,*)                                                      
CBS   write(LUPRI,*) 'ATOMIC NO-PAIR SO-MF CODE starts'
CBS   write(LUPRI,*)                                                      
      bonn=.false.
      sameorb=.false.
      aimp=.false.
      oneonly=.false. 
      makemean=.true.
CBS   write(LUPRI,*) stars                                                
CBS   write(LUPRI,*) '2e-integrals for the mean-field only'
CBS   write(LUPRI,*) '    mean-field will be generated         ' 
CBS   write(LUPRI,*) stars                                                
        do i=0,Lmax 
        icore(i)=0
        enddo   
       if (ifinite.eq.1) Exp_finite=EXP_FIN
      if (BONN) then
CBS   write(LUPRI,*) 'Bonn-approach for spin-other-orbit part'
      endif 
      if (BREIT) then
CBS   write(LUPRI,*) ' Breit-Pauli-Approximation'
      else
CBS   write(LUPRI,*) 'Douglas-Kroll type operators '
      endif 
      if (ifinite.eq.0) then 
CBS   write(LUPRI,*) 'Point-nucleus '
      else
CBS   write(LUPRI,*) 'Finite Nucleus' 
      endif  
CBS   write(LUPRI,*) stars                                                
CBS   write(LUPRI,*) 'write out one-electron integrals in MOLCAS-style'
CBS   write(LUPRI,*) '   and with MOLCAS normalization ' 
CBS   write(LUPRI,*) stars                                                
CBS   write(LUPRI,*) stars                                                
CBS   write(LUPRI,*) 
      symmetry='D2H' 
CBS   write(LUPRI,*) 'Symmetry is D2H'
CBS   write(LUPRI,*) 'check whether order of IRs is correct!!!'
      numbofsym=8    
      if (SAMEORB) then 
CBS   write(LUPRI,*) 'SAME-ORBIT only'
      else
CBS   write(LUPRI,*) 'OTHER-ORBIT included' 
      endif 
      if (AIMP) then 
CBS   write(LUPRI,*) 'CORE removed for use with AIMP' 
      endif 
      read(IN,*) charge,Lhigh
      if (Lhigh.gt.Lmax) then 
      write(LUPRI,*) 'Sorry, so far this code deals only ',
     *'with maximum l-values of ',Lmax
      CALL QUIT('Too high angular momentum values in AMFI')            
      endif  
CBS   write(LUPRI,*) ' Functions will go up to an L-value of : ',Lhigh
CBS   write(LUPRI,'(A19,F5.2)') chcharge,charge    
      call initired
      Do iredrun=1,numbofsym
      do Lrun=0,Lhigh  
      nmbMperIRL(iredrun,Lrun)=0 
      enddo 
      enddo 
      do Lrun=0,Lhigh 
CBS   write(LUPRI,*) 'ANGULAR MOMENTUM ',LRUN  
         read(IN,*) nprimit(Lrun),ncontrac(Lrun) 
CBS   write(LUPRI,'(I3,I3)') nprimit(Lrun),ncontrac(Lrun) 
cbs   check keywords
cbs   check maximum numbers
         if (nprimit(Lrun).gt.MxprimL) then 
            write(LUPRI,*) 'Too many primitives for L=',Lrun,
     *           ' increase MxprimL in para.h or reduce ',
     *           'the number of primitives to at least ',MxprimL
            CALL QUIT('Too many primitive functions in AMFI')            
         endif   
         if (ncontrac(Lrun).gt.MxcontL) then 
            write(LUPRI,*) 'Too many contracted fncts for L=',Lrun,
     *           ' increase MxcontL in para.h or ',
     *           'reduce the number of contracted functions',
     *           'to at most ',MxcontL
            CALL QUIT('Too many contracted functions in AMFI')            
         endif   
         if (ncontrac(Lrun).gt.nprimit(Lrun)) then 
            write(LUPRI,*) 'You have more contracted than ',
     *           'uncontracted functions, I don''t believe ',
     *           'that. Sorry!! '
            CALL QUIT('Inconsistent input detected in AMFI')            
         endif
C     write(LUPRI,'(A7,I3,A15,I3,A33,I3,A24)') 'For L= ',Lrun,
C    *' there will be ',
C    *ncontrac(Lrun),' contracted functions, built from ',
C    *nprimit(Lrun),
C    *' uncontracted functions.' 
         do ILINE=1,nprimit(Lrun)
            read(IN,*) exponents(ILINE,Lrun), 
     *           (cntscrtch(ILINE,JRUN,Lrun),
     *           Jrun=1,ncontrac(Lrun))
         enddo 
ckr         read(IN,'(A76)') header
c
cbs   
cbs   end of reading for the current L-value 
cbs   
c     do  Irun=1,ncontrac(Lrun)
c     writE(LUPRI,*) 'orbital : ',irun   
c     write(LUPRI,'(6(X,E12.6))') 
c    *(cntscrtch(I,Irun,Lrun),I=1,nprimit(Lrun))
c     enddo 
c     write(LUPRI,*) ' ' 
cbs   setting the numbers of cartesians per IR
         do iredrun=1,numbofsym 
            nfunctions(iredrun,Lrun)=0
         enddo
         do mrun=-Lrun,Lrun
            nfunctions(ipow2ired(ipowxyz(1,mrun,Lrun),
     *           ipowxyz(2,mrun,Lrun),Ipowxyz(3,mrun,Lrun)),Lrun)=
     *           nfunctions(ipow2ired(ipowxyz(1,mrun,Lrun),
     *           ipowxyz(2,mrun,Lrun),
     *           ipowxyz(3,mrun,Lrun)),Lrun)+ncontrac(Lrun)
         enddo
         do mrun=-Lrun,Lrun
            nmbMperIRL(ipow2ired(ipowxyz(1,mrun,Lrun),
     *           ipowxyz(2,mrun,Lrun),Ipowxyz(3,mrun,Lrun)),lruN)=
     *           nmbMperIRL(ipOw2ired(ipowxyz(1,mrun,Lrun),
     *           ipowxyz(2,mrun,Lrun),IpowxYz(3,mrun,Lrun)),lruN)+1
         enddo
CBS   write(LUPRI,*) stars                                                       
CBS   write(LUPRI,'(A,8I4)') 
CBS  *'Number of functions per IR: ',(nfunctions(iredrun,Lrun),
CBS  *iredrun=1,numbofsym)
CBS   write(LUPRI,*) stars                                                          
      enddo                     ! enddo for loop over L-values 
C     write(LUPRI,*) 'distribution of M-values'
c     do Lrun=0,Lhigh
c     write(LUPRI,*) (nmbMperIRL(nsym,Lrun),nsym=1,numbofsym)
c     endDo 
      numbofcart=0
      do lrun=0,Lhigh 
      numbofcart=numbofcart+(Lrun+Lrun+1)*
     *ncontrac(Lrun)
      enddo
      do iredrun=1,numbofsym 
      nfunctperIRED(iredrun)=0
      enddo
      do Lrun=0,Lhigh
      do iredrun=1,numbofsym 
      nfunctperIRED(iredrun)=nfunctperIRED(iredrun)+
     *nfunctions(iredrun,Lrun)
      enddo
      enddo
CBS   write(LUPRI,*) stars                                                        
CBS   write(LUPRI,'(A,8I3)') 'total number of atomic functions per IRED ',
CBS  *(nfunctperIRED(iredrun),iredrun=1,numbofsym) 
CBS   write(LUPRI,*) stars                                                        
      isum=0
      do iredrun=1,numbofsym
      itotalperIR(iredrun)=nfunctperIRED(iredrun)
      isum=isum+itotalperIR(iredrun)
      enddo 
      numballcart=isum
      iorbrun=0
      do iredrun=1,numbofsym 
      do inired=1,itotalperIR(iredrun)
      iorbrun=iorbrun+1 
      IREDoffunctnew(Iorbrun)=iredrun 
      enddo 
      enddo 
CBS   write(LUPRI,*) stars                                                        
CBS   write(LUPRI,'(A,8I3)') 'including additional functions per IRED ',
CBS  *(itotalperIR(iredrun),iredrun=1,numbofsym) 
CBS   write(LUPRI,*) stars                                                          
      do iredrun=1,numbofsym  
      ibeginIRED(iredrun)=0                              
      enddo
      do lrun=0,Lhigh
      do mrun=-lrun,lrun
      iredLM(mrun,lrun)=ipow2ired(ipowxyz(1,mrun,Lrun),
     *ipowxyz(2,mrun,Lrun),
     *ipowxyz(3,mrun,Lrun))
      incrLM(mrun,lrun)=ibeginIRED(iredLM(mrun,lrun))
      ibeginIRED(iredLM(mrun,lrUn))=
     *ibeginIRED(iredLM(mrun,lrun))+ncontrac(lrun) 
      enddo
      enddo   
c     do lrun=0,Lhigh
c     write(LUPRI,'(A,I4,A,21I3)') 'L= ',lrun,
c    *' shifts inside the IRED',
c    *(incrLM(mrun,lrun),mrun=-lrun,lrun)
c     enddo
      shiftIRED(1)=0
      do iredrun=2,numbofsym  
      shiftIRED(iredrun)=shiftIRED(iredrun-1)
     *                   +itotalperIR(iredrun-1)
      enddo
c     write(LUPRI,'(A,8I4)') 'shifts for the IREDs ',
c    *(shiftIRED(iredrun),iredrun=1,numbofsym) 
cbs   test all orbital numbers
c     do lrun=0,Lhigh
c     do mrun=-Lrun,Lrun
c     do irun=1,ncontrac(lrun)
c     write(LUPRI,*) 'L,M,contr funct, absolute number ',
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
c     write(LUPRI,*) 'ired1,ired2 ',ired1,ired2,
c    *irun,shiftIRIR(irun)
      enddo
      enddo
cbs  
      do lrun=0,Lhigh
      do Mrun=-Lrun,Lrun
      ired=iredLM(Mrun,Lrun)
      ishifter=shiftIRED(ired)+incrLM(mrun,lrun)      
      do icart=1,ncontrac(Lrun)
      moffunction(ishifter+icart)=Mrun
      Loffunction(ishifter+icart)=Lrun
      IREDoffunction(ishifter+Icart)=ired  
      INOFT(ishifter+Icart)=icart     
      enddo
      enddo
      enddo
CBS   write(LUPRI,*) stars
CBS   write(LUPRI,*) 'SYMMETRY-INFORMATION ON FUNCTIONS '
CBS   write(LUPRI,*) stars
      do irun = 1, numbofcart
CBS   write(LUPRI,'(4(A,I3))') 'Number of function: ',
CBS  *irun,
CBS  *' IR of function: ',IREDoffunction(irun),
CBS  *' L-value: ',Loffunction(irun),
CBS  *' M-value: ',Moffunction(irun)
CBS   numboffunct(irun)=irun
      INOFF(irun)=irun
CBS   if (IREDoffunction(irun).ne.IREDoffunction(irun+1)) 
CBS  *write(LUPRI,*)   
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
4712  if (irun.eq.1.or.(irun.ge.2.and.INOFF(irun).eq.  
     *numbprev+1)) then 
      Lval=Loffunction(irun)
      number=INOFF(irun)
      itype=INOFT(irun)  
      if (itype.le.icore(lval)) then
      write(LUPRI,777) number,itype,lval
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
CBS   write(LUPRI,*) stars                                                          
      write(LUPRI,'(A,8I3)')'# of funct. per IRED after removing core ',
     *(nrtofiperIR(iredrun),iredrun=1,numbofsym) 
      write(LUPRI,*) ikeeporb,' orbitals left after deleting core' 
      endif 
CBS   write(LUPRI,*) stars
      nmax=max(6,ncontrac(0))  
      do lrun=1,Lhigh
      nmax=max(nmax,ncontrac(lrun))
      enddo
      ncont4=nmax*nmax*nmax*nmax
      return 
777   format('ORBITAL NUMBER ',I4,' IS THE ',I3,'TH of L-value ',I3,
     *' IT WILL BE REMOVED !!!')  
      end 
      double precision function  regge3j(
     *j1,     ! integer  2*j1
     *j2,     ! integer  2*j2
     *j3,     ! integer  2*j3
     *m1,     ! integer  2*m1
     *m2,     ! integer  2*m2
     *m3)     ! integer  2*m3
cbs   uses magic square of regge (see Lindner pp. 38-39)
cbs 
cbs    ---                                            ---
cbs   |                                                  |
cbs   | -j1+j2+j3     j1-j2+j3         j1+j2-j3          |
cbs   |                                                  |
cbs   |                                                  |
cbs   |  j1-m1        j2-m2            j3-m3             |
cbs   |                                                  |
cbs   |                                                  |
cbs   |  j1+m1        j2+m2            j3+m3             |
cbs   |                                                  |
cbs    ---                                            ---
cbs 
#include "implicit.h"
      dimension MAT(3,3)
      logical testup,testdown
#include "Regge.h"
cbs  facul,   integer array (nprim,0:mxLinRE) prime-expansion of factorials 
cbs  mxLinRE,    integer max. number for facul is given
cbs  nprim,   number of primes for expansion of factorials 
cbs  prim,    integer array with the first nprim prime numbers
cbs  iwork)   integer array of size nprim
      regge3j=0d0
c     write(6,'(A24,6I3)') '3J to be calculated for ',
c    *j1,j2,j3,m1,m2,m3
cbs   quick check  if =/= 0 at all
      icheck=m1+m2+m3 
      if (icheck.ne.0) then 
c     write(6,*) 'sum over m =/= 0'
      return    
      endif 
cbs   check triangular relation (|j1-j2|<= j3 <= j1+j2 )
      imini=iabs(j1-j2)
      imaxi=j1+j2 
      if (j3.lt.imini.or.j3.gt.imaxi) then 
c     write(6,*) 'triangular relation not fulfilled'
      return   
      endif
cbs   quick check  if =/= 0 at all  end 
cbs  
cbs   3J-symbol is not zero by simple rules 
cbs  
cbs   initialize MAT 
      MAT(1,1) =-j1+j2+j3
      MAT(2,1) =j1-m1    
      MAT(3,1) =j1+m1     
      MAT(1,2) =j1-j2+j3  
      MAT(2,2) =j2-m2      
      MAT(3,2) =j2+m2     
      MAT(1,3) =j1+j2-j3  
      MAT(2,3) =j3-m3      
      MAT(3,3) =j3+m3      
      do I=1,3
      do J=1,3
cbs   check for even numbers (2*integer) and positive or zero
      if (mod(MAT(J,I),2).ne.0.or.MAT(J,I).lt.0)  then 
c     write(6,*) 'J,I,MAT(J,I): ',J,I,MAT(J,I)
      return
      endif 
      MAT(J,I)=MAT(J,I)/2
      if (Mat(j,i).gt.mxLinRE)
     *CALL QUIT('increase mxLinRE for regge3j')
      enddo
      enddo
      Isigma=(j1+j2+j3)/2
cbs   check the magic sums
      do I=1,3
      IROW=0
      ICOL=0
      do J=1,3
      IROW=IROW+MAT(I,J)
      ICOL=ICOL+MAT(J,I)
      enddo
      if (IROW.ne.Isigma.or.ICOL.ne.Isigma) then 
c     write(6,*) 'I,IROW,ICOL ',I,IROW,ICOL  
      return
      endif
      enddo
cbs   if j1+j2+j3 is odd: check for equal rows or columns 
      Isign=1
      if (iabs(mod(Isigma,2)).eq.1) then 
      isign=-1
         do I=1,3
         do J=I+1,3
            if (MAT(1,I).eq.MAT(1,J).and. 
     *         MAT(2,I).eq.MAT(2,J).and.
     *         MAT(3,I).eq.MAT(3,J)) return 
            if (MAT(I,1).eq.MAT(J,1).and. 
     *         MAT(I,2).eq.MAT(J,2).and.
     *         MAT(I,3).eq.MAT(J,3)) return 
         enddo
         enddo
      endif 
cbs   look for the lowest element indices: IFIRST,ISECOND
      imini=MAT(1,1) 
      IFIRST=1
      ISECOND=1
      do I=1,3
      do J=1,3 
      if (MAT(J,I).lt.imini) then 
      IFIRST=J
      ISECOND=I
      imini=MAT(J,I)
      endif 
      enddo
      enddo
c     write(6,*) 'Matrix before commuting vectors'
      do ibm=1,3
c     write(6,'(3I5)') (Mat(ibm,j),j=1,3) 
      enddo
      if (IFIRST.ne.1) then  !interchange rows
c     write(6,*) 'IFIRST = ',ifirst
      do I=1,3
      IDUMMY=MAT(1,I) 
      MAT(1,I)=MAT(IFIRST,I)
      MAT(IFIRST,I)=IDUMMY
      enddo
      endif 
      if (ISECOND.ne.1) then  !interchange columns
c     write(6,*) 'ISECOND = ',isecond
      do I=1,3
      IDUMMY=MAT(I,1) 
      MAT(I,1)=MAT(I,ISECOND)
      MAT(I,ISECOND)=IDUMMY
      enddo
      endif 
cbs   lowest element is now on (1,1)
c     write(6,*) 'Matrix after commuting vectors'
c     do ibm=1,3
c     write(6,'(3I5)') (Mat(ibm,j),j=1,3) 
c     enddo
cbs   begin to calculate Sum over s_n
cbs   first the simple cases
      if (Mat(1,1).eq.0) then 
      isum=1 
      elseif (Mat(1,1).eq.1) then 
      isum=Mat(2,3)*Mat(3,2)-Mat(2,2)*Mat(3,3)
      elseif (Mat(1,1).eq.2) then 
      isum=Mat(2,3)*(Mat(2,3)-1)*Mat(3,2)*(Mat(3,2)-1)-
     *2*Mat(2,3)*Mat(3,2)*Mat(2,2)*Mat(3,3)+
     *Mat(2,2)*(Mat(2,2)-1)*Mat(3,3)*(Mat(3,3)-1)
      else !  all the cases with Mat(1,1) >= 3 
        Icoeff=1
        do Ibm=Mat(3,2)-Mat(1,1)+1,Mat(3,2)
          icoeff=icoeff*ibm
        enddo
        do Ibm=Mat(2,3)-Mat(1,1)+1,Mat(2,3)
          icoeff=icoeff*ibm
        enddo
        isum=icoeff
        do Icount=1,MAT(1,1) 
           icoeff=-icoeff*(Mat(1,1)+1-icount)*(Mat(2,2)+1-icount)*
     *           (Mat(3,3)+1-icount)
           Idenom=icount*(Mat(2,3)-Mat(1,1)+icount)*
     *           (Mat(3,2)-Mat(1,1)+icount)
           icoeff=icoeff/Idenom
           isum=isum+icoeff
        enddo
      endif   
cbs  additional sign from interchanging rows or columns
      if (ifirst.ne.1) isum=isum*isign
      if (isecond.ne.1) isum=isum*isign
c     write(6,*) 'isum = ',isum 
cbs       Mat(2,3)+Mat(3,2) 
cbs    (-) 
      if (iabs(mod((Mat(2,3)+Mat(3,2)),2)).eq.1) isum=-isum   
cbs   final factor
      LIMIT=ihigh(max(Mat(1,1),Mat(1,2),Mat(1,3),
     *Mat(2,1),Mat(2,2),Mat(2,3),Mat(3,1),Mat(3,2),
     *Mat(3,3),(Isigma+1)))
      do I=1,LIMIT 
      iwork(I)=facul(I,Mat(1,2))+facul(I,Mat(2,1))+
     *facul(I,Mat(3,1))+facul(I,Mat(1,3))-
     *facul(I,Mat(1,1))-facul(I,Mat(2,2))-
     *facul(I,Mat(3,3))-facul(I,(Isigma+1))-
     *facul(I,Mat(2,3))-facul(I,Mat(3,2))
      enddo
c     write(6,*) 'Iwork: ',(iwork(i),i=1,LIMIT)
      factor=1d0
      iup=1
      idown=1
      testup=.true.
      testdown=.true.
      do I=1,LIMIT   
      do J=1,iwork(I)
      iup=iup*prim(i)
      if (iup.lt.0) testup=.false. !check for Integer overflow
      enddo 
      Enddo 
      up=dfloat(iup) 
      if(.not.testup) then ! if the integers did not run correctly  
        up=1d0
        do I=1,LIMIT   
              do J=1,iwork(I)
              up=up*dfloat(prim(i))
              enddo 
        enddo 
      endif 
      do I=1,LIMIT   
      do J=1,-iwork(I)
      idown=idown*prim(i)
      if (idown.lt.0) testdown=.false. 
      enddo 
      enddo 
      down=dfloat(idown) 
      if(.not.testdown) then 
        down=1d0
        do I=1,LIMIT   
              do J=1,-iwork(I)
              down=down*dfloat(prim(i))
              enddo 
        enddo 
      endif 
c     if (.not.(testup.and.testdown)) then 
c     write(6,*) 'j1,j2,j3,m1,m2,m3 ',j1,j2,j3,m1,m2,m3
c     write(6,*) 'iup,idown ',iup,idown,'up,down ',up,down 
c     endif 
      factor=factor*up/down
cbs   final result
      regge3j=dsqrt(factor)*dfloat(isum)
      return 
      end 
      double precision function Tkinet(l,alpha1,alpha2)
cbs   calculates the matrix element of kinetic energy 
cbs   for primitive normalized functions with the same angular momentum l 
cbs   and exponents alpha1 and alpha2 
cbs   works only, if r**l is assumed for an l-value
cbs   formular obtained from the symmetric expression (d/dr's to (')
cbs   the left and to the right. 
cbs   Overlaps of the different powers are partially crossed out 
cbs   with  the overlap of functions with angular momentum l
cbs   final formula:
cbs   Tkinet=0.5*alpha12 (2l+3) (alpha1*alpha2/alpha12*alpha12)**((2L+7)/4)
cbs   with alpha12=0.5*(alpha1+alpha2)
cbs   as alpha12 has the dimensions 1/length**2, this can not be that bad...
      Implicit double precision (a-h,o-z)
Cbs   alpha12 is the effective exponent 
      Alpha12=0.5d0*(alpha1+alpha2)
      alphpro=alpha1*alpha2   
      ll3=l+l+3        
      ll7=l+l+7
      Tkinet=0.5d0*alpha12*ll3*(alphpro/
     *(alpha12*alpha12))**(0.25*dfloat(ll7))
      return 
      end    
      subroutine   tosigX(m1,m2,m3,m4,angint,
     *mcombina,ncontl1,ncontl2,ncontl3,
     *ncontl4,carteX,preXZ,interxyz,isgnprod,
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
cbs   | L,M,-> = -i 2**(-0.5) ( (-1)**M Y(L,M) - Y(L,-M))  ($$$$) 
cbs
cbs
cbs   due to symmetry, there can be only integrals 
cbs   with indices one or three  (sigma_+ and sigma_-)- combinations 
cbs
#include "implicit.h"
#include "para.h"
#include "priunit.h"
      logical cleaner
      dimension mcombina(2,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *angint(ncontl1,ncontl2,ncontl3,ncontl4,*),
cbs  !!!!!!!!!!!changing now to the order of chemists notation!!!!!!!!!!
     *carteX(ncontl1,ncontl3,ncontl2,ncontl4),
     *preXZ(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *interxyz(*),
     *isgnprod(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),    
     *isgnM(-1:1,-1:1,-1:1,-1:1)
c     write(6,*) ' begin tosigx' 
cbs   cleaning up the integral-array
      irun=ncontl1*ncontl2*ncontl3*ncontl4 
      call dzero(cartex,irun)
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
      write(LUPRI,*) 'tosigx: no interaction: ',m1,m2,m3,m4
      CALL QUIT('Error in TOSIGX in AMFI')            
      endif 
      prexz1234=preXZ(m1,m2,m3,m4)
      do while (interxyz(irun+1).gt.0) 
      irun=irun+1
c     write(6,*) 'tosigx: ',irun,interxyz(irun)
c
cbs
cbs
cbs   This could be done with gotos, but I am biased to hate those..
cbs
cbs
         if (interxyz(irun).eq.1) then   
         ityp=mcombina(1,Mabs1,Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,Mabs3,Mabs4)
         factor=isgnM(1,1,1,1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.2) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,-Mabs1,-Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(-1,-1,-1,-1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,-Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.3) then   
         ityp=mcombina(1,Mabs1,Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,Mabs3,-Mabs4)
         factor=isgnM(1,1,1,-1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.4) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,-Mabs1,-Mabs2,-Mabs3,Mabs4)
         factor=isgnM(-1,-1,-1,1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,-Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.5) then  
         ityp=mcombina(1,Mabs1,Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,-Mabs3,Mabs4)
         factor=isgnM(1,1,-1,1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,-Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.6) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,-Mabs1,-Mabs2,Mabs3,-Mabs4)
         factor=isgnM(-1,-1,1,-1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.7) then   
         ityp=mcombina(1,Mabs1,-Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,Mabs3,Mabs4)
         factor=isgnM(1,-1,1,1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.8) then   
         ityp=mcombina(1,-Mabs1,Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,-Mabs1,Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(-1,1,-1,-1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,-Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.9) then   
         ityp=mcombina(1,-Mabs1,Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,-Mabs1,Mabs2,Mabs3,Mabs4)
         factor=isgnM(-1,1,1,1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.10) then   
         ityp=mcombina(1,Mabs1,-Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(1,-1,-1,-1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,-Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.11) then   
         ityp=mcombina(1,Mabs1,Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(1,1,-1,-1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,-Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.12) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,-Mabs1,-Mabs2,Mabs3,Mabs4)
         factor=isgnM(-1,-1,1,1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.13) then   
         ityp=mcombina(1,Mabs1,-Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,Mabs3,-Mabs4)
         factor=isgnM(1,-1,1,-1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.14) then   
         ityp=mcombina(1,-Mabs1,Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,-Mabs1,Mabs2,-Mabs3,Mabs4)
         factor=isgnM(-1,1,-1,1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,-Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.15) then   
         ityp=mcombina(1,Mabs1,-Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,-Mabs3,Mabs4)
         factor=isgnM(1,-1,-1,1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,-Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.16) then   
         ityp=mcombina(1,-Mabs1,Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,-Mabs1,Mabs2,Mabs3,-Mabs4)
         factor=isgnM(-1,1,1,-1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteX,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
         endif  
       enddo 
        if (cleaner) then   
        do irun4=1,ncontl4
        do irun2=1,ncontl2
        do irun1=1,ncontl1
        cartex(irun1,irun1,irun2,irun4)=0d0
        enddo
        enddo
        enddo
        endif 
      return 
      end 
      subroutine   tosigY(m1,m2,m3,m4,angint,
     *mcombina,ncontl1,ncontl2,ncontl3,
     *ncontl4,carteY,preY,interxyz,isgnprod,
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
#include "implicit.h"
#include "priunit.h"
#include "para.h"
      dimension mcombina(2,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *angint(ncontl1,ncontl2,ncontl3,ncontl4,*),
cbs  !!!!!!!!!!!changing now to the order of chemists notation!!!!!!!!!!
     *carteY(ncontl1,ncontl3,ncontl2,ncontl4),
     *preY(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *interxyz(*),
     *isgnprod(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *isgnM(-1:1,-1:1,-1:1,-1:1)
      logical cleaner
c     write(6,*) 'begin tosigy '
cbs   cleaning up the integral-array
      irun=ncontl4*ncontl2*ncontl3*ncontl1
      call dzero(carteY,irun)
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
      write(LUPRI,*) 'tosigy: no interaction: ',m1,m2,m3,m4
      CALL QUIT('Error in TOSIGY in AMFI')
      endif
      prey1234=preY(m1,m2,m3,m4)
c     write(6,*) 'prey ',prey1234
      do while (interxyz(irun+1).gt.0) 
      irun=irun+1
c     write(6,*) 'tosigy: ',irun,interxyz(irun)
c
cbs
cbs
cbs   This could be done with gotos, but I am biased to hate those..
cbs
cbs
         if (interxyz(irun).eq.1) then   
         ityp=mcombina(1,Mabs1,Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,Mabs3,Mabs4)
         factor=isgnM(1,1,1,1)*prey1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.2) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,-Mabs1,-Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(-1,-1,-1,-1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,-Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.3) then   
         ityp=mcombina(1,Mabs1,Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,Mabs3,-Mabs4)
         factor=isgnM(1,1,1,-1)*prey1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.4) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,-Mabs1,-Mabs2,-Mabs3,Mabs4)
         factor=isgnM(-1,-1,-1,1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,-Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.5) then   
         ityp=mcombina(1,Mabs1,Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,-Mabs3,Mabs4)
         factor=isgnM(1,1,-1,1)*prey1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,-Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.6) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,-Mabs1,-Mabs2,Mabs3,-Mabs4)
         factor=isgnM(-1,-1,1,-1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.7) then   
         ityp=mcombina(1,Mabs1,-Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,Mabs3,Mabs4)
         factor=isgnM(1,-1,1,1)*prey1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.8) then   
         ityp=mcombina(1,-Mabs1,Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,-Mabs1,Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(-1,1,-1,-1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,-Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.9) then   
         ityp=mcombina(1,-Mabs1,Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,-Mabs1,Mabs2,Mabs3,Mabs4)
         factor=isgnM(-1,1,1,1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.10) then   
         ityp=mcombina(1,Mabs1,-Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(1,-1,-1,-1)*prey1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,-Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.11) then  
         ityp=mcombina(1,Mabs1,Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(1,1,-1,-1)*prey1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,-Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.12) then   
         ityp=mcombina(1,-Mabs1,-Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,-Mabs1,-Mabs2,Mabs3,Mabs4)
         factor=isgnM(-1,-1,1,1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.13) then   
         ityp=mcombina(1,Mabs1,-Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,Mabs3,-Mabs4)
         factor=isgnM(1,-1,1,-1)*prey1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,Mabs3,-Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.14) then   
         ityp=mcombina(1,-Mabs1,Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,-Mabs1,Mabs2,-Mabs3,Mabs4)
         factor=isgnM(-1,1,-1,1)*prey1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,-Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.15) then   
         ityp=mcombina(1,Mabs1,-Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,-Mabs3,Mabs4)
         factor=isgnM(1,-1,-1,1)*prey1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,-Mabs3,Mabs4))
         if (ityp.eq.3) factor=-factor    
         call daxpint(angint(1,1,1,1,iblock),carteY,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.16) then   
         ityp=mcombina(1,-Mabs1,Mabs2,Mabs3,-Mabs4) 
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
        endif
      return 
      end 
      subroutine   tosigZ(m1,m2,m3,m4,angint,
     *mcombina,ncontl1,ncontl2,ncontl3,
     *ncontl4,carteZ,preXZ,interxyz,isgnprod,
     *cleaner) 
cbs   this subroutine combines the angular integrals 
cbs   to the integrals for the real-valued linear 
cbs   combinations for the sigma_Z part 
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
cbs   only angular integrals of type 2 (sigma_0) contribute  
cbs
#include "implicit.h"
#include "priunit.h"
#include "para.h"
      dimension mcombina(2,-Lmax:Lmax,-Lmax:Lmax,
     *-Lmax:Lmax,-Lmax:Lmax),
     *angint(ncontl1,ncontl2,ncontl3,ncontl4,*),
cbs  !!!!!!!!!!!changing now to the order of chemists notation!!!!!!!!!!
     *carteZ(ncontl1,ncontl3,ncontl2,ncontl4),
     *preXZ(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *interxyz(*),
     *isgnprod(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *isgnM(-1:1,-1:1,-1:1,-1:1)                          
      logical cleaner
cbs   cleaning up the integral-array
      irun=ncontl4*ncontl2*ncontl3*ncontl1
      call dzero(carteZ,irun)
c     write(6,*) 'begin tosigz'
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
      write(LUPRI,*) 'tosigz: no interaction: ',m1,m2,m3,m4
      CALL QUIT('Error in TOSIGZ in AMFI')
      endif
      prexz1234=preXZ(m1,m2,m3,m4)
      do while (interxyz(irun+1).gt.0) 
      irun=irun+1
c
cbs
cbs
cbs   This could be done with gotos, but I am biased to hate those..
cbs
cbs
         if (interxyz(irun).eq.1) then   
         ityp=mcombina(1,Mabs1,Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,Mabs3,Mabs4)
         factor=isgnM(1,1,1,1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.2) then   
         ityp=mcombina(1,Mabs1,Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,Mabs3,Mabs4)
         factor=-isgnM(-1,-1,-1,-1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,-Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.3) then   
         ityp=mcombina(1,Mabs1,Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,Mabs3,-Mabs4)
         factor=isgnM(1,1,1,-1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.4) then   
         ityp=mcombina(1,Mabs1,Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,Mabs3,-Mabs4)
         factor=-isgnM(-1,-1,-1,1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,-Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.5) then  
         ityp=mcombina(1,Mabs1,Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,-Mabs3,Mabs4)
         factor=isgnM(1,1,-1,1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,-Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.6) then   
         ityp=mcombina(1,Mabs1,Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,-Mabs3,Mabs4)
         factor=-isgnM(-1,-1,1,-1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.7) then   
         ityp=mcombina(1,Mabs1,-Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,Mabs3,Mabs4)
         factor=isgnM(1,-1,1,1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.8) then   
         ityp=mcombina(1,Mabs1,-Mabs2,Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,Mabs3,Mabs4)
         factor=-isgnM(-1,1,-1,-1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,-Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.9) then   
         ityp=mcombina(1,Mabs1,-Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,-Mabs3,-Mabs4)
         factor=-isgnM(-1,1,1,1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.10) then   
         ityp=mcombina(1,Mabs1,-Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(1,-1,-1,-1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,-Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.11) then   
         ityp=mcombina(1,Mabs1,Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,-Mabs3,-Mabs4)
         factor=isgnM(1,1,-1,-1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,Mabs2,-Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.12) then   
         ityp=mcombina(1,Mabs1,Mabs2,-Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,Mabs2,-Mabs3,-Mabs4)
         factor=-isgnM(-1,-1,1,1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,-Mabs2,Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.13) then   
         ityp=mcombina(1,Mabs1,-Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,Mabs3,-Mabs4)
         factor=isgnM(1,-1,1,-1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.14) then   
         ityp=mcombina(1,Mabs1,-Mabs2,Mabs3,-Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,Mabs3,-Mabs4)
         factor=-isgnM(-1,1,-1,1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,-Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.15) then   
         ityp=mcombina(1,Mabs1,-Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,-Mabs3,Mabs4)
         factor=isgnM(1,-1,-1,1)*prexz1234*
     *   dfloat(isgnprod(Mabs1,-Mabs2,-Mabs3,Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         elseif (interxyz(irun).eq.16) then   
         ityp=mcombina(1,Mabs1,-Mabs2,-Mabs3,Mabs4) 
         iblock=mcombina(2,Mabs1,-Mabs2,-Mabs3,Mabs4)
         factor=-isgnM(-1,1,1,-1)*prexz1234*
     *   dfloat(isgnprod(-Mabs1,Mabs2,Mabs3,-Mabs4))
         call daxpint(angint(1,1,1,1,iblock),carteZ,
     *   factor,ncontl1,ncontl2,ncontl3,ncontl4) 
c
         endif  
      enddo 
        if (cleaner) then   
        do irun4=1,ncontl4
        do irun2=1,ncontl2
        do irun1=1,ncontl1
        cartez(irun1,irun1,irun2,irun4)=0d0
        enddo
        enddo
        enddo
        endif
      return 
      end 
      subroutine trans(
cbs   makes the transformation for the ich-th index
     *coeffs, !(nolds(ith),nnew(ith)) modified contraction coefficients
     *idim1,  !  first dimension   
     *idim2,  !  second dimension  
     *ich,    ! index to be changed                         
     *nolds1,nolds2,nolds3,nolds4,  ! old dimensions                    
     *nnew1,nnew2,nnew3,nnew4, ! new dimensions      
     *array1, ! array of size (nolds1,nolds2,nolds3,nolds4)
     *array2  ! array of size (nnew1,nnew2,nnew3,nnew4)
     *)
#include "implicit.h"
      dimension coeffs(idim1,idim2),           
     *array1(nolds1,nolds2,nolds3,nolds4),
     *array2(nnew1,nnew2,nnew3,nnew4)
c     write(6,*) 'begin trans ' ,ich 
c     write(6,'(8I5)') nolds1,nolds2,nolds3,nolds4,
c    *nnew1,nnew2,nnew3,nnew4
      do ind4=1,nnew4
      do ind3=1,nnew3
      do ind2=1,nnew2
      do ind1=1,nnew1
      array2(ind1,ind2,ind3,ind4)=0d0
      enddo
      enddo
      enddo
      enddo
      if (ich.eq.1) then 
      do ind4=1,nnew4
      do ind3=1,nnew3
      do ind2=1,nnew2
      do ind5=1,nnew1  
      do ind1=1,nolds1
      array2(ind5,ind2,ind3,ind4)=array2(ind5,ind2,ind3,ind4)+
     *coeffs(ind1,ind5)*array1(ind1,ind2,ind3,ind4)
      enddo
      enddo
      enddo
      enddo
      enddo
      elseif (ich.eq.2) then 
c     write(6,*) 'transform second index '
      do ind4=1,nnew4
      do ind3=1,nnew3
      do ind5=1,nnew2  
      do ind2=1,nolds2
      coeff=coeffs(ind2,ind5)
      do ind1=1,nnew1
      array2(ind1,ind5,ind3,ind4)=array2(ind1,ind5,ind3,ind4)+
     *coeff*array1(ind1,ind2,ind3,ind4)
      enddo
      enddo
      enddo
      enddo
      enddo
c     write(6,*) 'end  to transform second index '
      elseif (ich.eq.3) then 
      do ind4=1,nnew4
      do ind5=1,nnew3  
      do ind3=1,nolds3
      coeff=coeffs(ind3,ind5)
      do ind2=1,nnew2
      do ind1=1,nnew1
      array2(ind1,ind2,ind5,ind4)=array2(ind1,ind2,ind5,ind4)+
     *coeff*array1(ind1,ind2,ind3,ind4)
      enddo
      enddo
      enddo
      enddo
      enddo
      elseif (ich.eq.4) then 
      do ind5=1,nnew4  
      do ind4=1,nolds4
      coeff=coeffs(ind4,ind5)
      do ind3=1,nnew3
      do ind2=1,nnew2
      do ind1=1,nnew1
      array2(ind1,ind2,ind3,ind5)=array2(ind1,ind2,ind3,ind5)+
     *coeff*array1(ind1,ind2,ind3,ind4)
      enddo
      enddo
      enddo
      enddo
      enddo
      endif  
c     write(6,*) 'end  trans ' 
      return
      end
      subroutine transcon(contold,idim1,idim2,ovlp,contnew,nprim,ncont)
#include "implicit.h"
      dimension contold(idim1,idim2),contnew(nprim,ncont),
     *ovlp(idim1,idim1)
c     write(6,*) 'begin transcon nprim,ncont ',nprim,ncont 
cbs   copy old contraction coefficients in dense form to common block 
      do Jrun=1,ncont
      do Irun=1,nprim
      contnew(Irun,Jrun)=contold(Irun,Jrun)
      enddo 
      enddo 
cbs   ensure normalization 
      do ICONT=1,ncont
        xnorm=0d0   
        do Jrun=1,nprim
        do Irun=1,nprim
        xnorm=xnorm+contnew(Irun,ICONT)*contnew(Jrun,ICONT)
     *  *ovlp(Irun,Jrun)
c       write(6,*) 'Icont,jrun,irun,xnorm ',
c    *  icont,jrun,irun,xnorm 
        enddo 
        enddo 
c       write(6,*) 'ICONT ',ICONT,xnorm   
        xnorm=1d0/dsqrt(xnorm)
cbs   scale with normalization factor 
        do Irun=1,nprim
        contnew(Irun,ICONT)=xnorm*contnew(Irun,ICONT)
        enddo 
      enddo 
c     write(6,*) 'end transcon nprim,ncont ',nprim,ncont 
      return 
      end 
      subroutine two2mean12a(carteSO,carteOO,occup,AOcoeffs,onecart,
     *ncontmf,norbsum,noccorb,sameorb)
#include "implicit.h"
#include "para.h"
      logical sameorb 
      dimension 
     *carteSO(ncontmf,norbsum,ncontmf,norbsum),
     *carteOO(ncontmf,norbsum,ncontmf,norbsum),
     *occup(*),AOcoeffs(MxcontL,*),onecart(MxcontL,MxcontL)
      if (sameorb) THEN 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5d0*coeff  
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)-coeff*
     *carteSO(irun,icartleft,jrun,icartright)
      enddo
      enddo
      enddo
      enddo
      else 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5d0*coeff  
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)-coeff*
     *(carteSO(irun,icartleft,jrun,icartright)+
     *2d0*carteOO(irun,icartleft,jrun,icartright))
      enddo
      enddo
      enddo
      enddo
      endif 
      return 
      end 
      
      subroutine two2mean12b(carteSO,carteOO,occup,AOcoeffs,onecart,
     *ncontmf,norbsum,noccorb,sameorb)
#include "implicit.h"
#include "para.h"
      logical sameorb 
      dimension 
     *carteSO(ncontmf,norbsum,ncontmf,norbsum),
     *carteOO(ncontmf,norbsum,ncontmf,norbsum),
     *occup(*),AOcoeffs(MxcontL,*),onecart(MxcontL,MxcontL)
      if (sameorb) then 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5d0*coeff  
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)+coeff*
     *carteSO(jrun,icartleft,irun,icartright)
      enddo
      enddo
      enddo
      enddo
      else 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5d0*coeff  
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)+coeff*
     *(carteSO(jrun,icartleft,irun,icartright)+
     *2d0*carteOO(jrun,icartleft,irun,icartright))
      enddo
      enddo
      enddo
      enddo
      endif 
      return 
      end 
      
      subroutine two2mean13(carteSO,occup,AOcoeffs,onecart,
     *ncontmf,norbsum,noccorb)
cbs   gives the two first contributions
cbs   < i M | j M >  with Malpha  and Mbeta 
cbs   the other orbit parts cancel    
#include "implicit.h"
#include "para.h"
      dimension carteSO(ncontmf,ncontmf,norbsum,norbsum),
     *occup(*),AOcoeffs(MxcontL,*),onecart(MxcontL,MxcontL)
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb        
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)+coeff*
     *carteSO(irun,jrun,icartleft,icartright)
      enddo
      enddo
      enddo
      enddo
c     write(6,*) 'effective integrals' 
c     do jrun=1,ncontmf
c     write(6,'(4E20.14)') (onecart(irun,jrun),irun=1,ncontmf)
c     enddo
      return 
      end 
      
      subroutine two2mean34a(carteSO,carteOO,occup,AOcoeffs,onecart,
     *ncontmf,norbsum,noccorb,sameorb)
#include "implicit.h"
#include "para.h"
      logical sameorb 
      dimension 
     *carteSO(norbsum,ncontmf,norbsum,ncontmf),
     *carteOO(norbsum,ncontmf,norbsum,ncontmf),
     *occup(*),AOcoeffs(MxcontL,*),onecart(MxcontL,MxcontL)
      if (sameorb) then 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb   
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5d0*coeff
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)+coeff*
     *carteSO(icartleft,irun,icartright,jrun)
      enddo
      enddo
      enddo
      enddo
      else 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb   
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5d0*coeff
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)+coeff*
     *(carteSO(icartleft,irun,icartright,jrun)+
     *2d0*carteOO(icartleft,irun,icartright,jrun))
      enddo
      enddo
      enddo
      enddo
      endif 
      return 
      end 
      
      subroutine two2mean34b(carteSO,carteOO,occup,AOcoeffs,onecart,
     *ncontmf,norbsum,noccorb,sameorb)
#include "implicit.h"
#include "para.h"
      logical sameorb 
      dimension 
     *carteSO(norbsum,ncontmf,norbsum,ncontmf),
     *carteOO(norbsum,ncontmf,norbsum,ncontmf),
     *occup(*),AOcoeffs(MxcontL,*),onecart(MxcontL,MxcontL)
      if (sameorb) then 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb   
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5D0*coeff  
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)-coeff*
     *carteSO(icartleft,jrun,icartright,irun)
      enddo
      enddo
      enddo
      enddo
      else 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb   
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5D0*coeff  
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)-coeff*
     *(carteSO(icartleft,jrun,icartright,irun)+
     *2d0*carteOO(icartleft,jrun,icartright,irun))
      enddo
      enddo
      enddo
      enddo
      endif 
      return 
      end 
