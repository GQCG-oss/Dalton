      Double precision function LMdepang(
     *L,M,l1,l2,l3,l4,m1,m2,m3,m4,cheater)
cbs   l1-l4 and m1-m4 are already shifted !!
cbs   purpose: calculates the angular part of the   
cbs   coulomb-type integrals. See documentation for details...
cbs   LMdepang= LM dependent angular factors 
cbs   cheater included for a correcting signs, as there were some 
cbs   signs (only signs!!!!) missing when compared to HERMIT  
cbs                                        B.S.  08.10.96 
      implicit double precision (a-h,o-z) 
      common /pifact/ pi,root8ovpi
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
      isum=-l3-l1-l4-l2+2*(M+m3+m4)   !???? I'm not sure 
      if (mod(isum,4).eq.0) then 
      isign=1
      elseif (iabs(mod(isum,4)).eq.2) then 
      isign=-1
      else 
      write(6,*) 'L,l1,l2,l3,l4,M,m1,m2,m3,m4'
      write(6,'(10I3)') L,l1,l2,l3,l4,M,m1,m2,m3,m4
      write(6,*) 'isum= ',isum,' mod = ',mod(isum,4)
      stop 'error in lmdepang' 
      endif
      fact2=couple3J(L,l3,l1,-M,m3,-m1) 
      fact3=couple3J(L,l4,l2,M,m4,-m2)
C     write(6,*) 'fact2,fact3 ',fact2,fact3
      LMdepang=cheater*dfloat(isign)*fact1*fact2*fact3
      return 
      end 
