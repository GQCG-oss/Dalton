      double precision function  couple3J(
     *l1,     ! integer  l1
     *l2,     ! integer  l2
     *l3,     ! integer  l3
     *m1,     ! integer  m1
     *m2,     ! integer  m2
     *m3)     ! integer  m3
cbs this routine calculates the coupling of three angular momenta to  zero   
cbs   
cbs   
cbs   Int dOmega i^(l1+l2+l3) Y^l1_m1 (Omega) Y^l2_m2 (Omega) Y^l3_m3 (Omega) =
cbs   sqrt( (2l1+1)(2l2+1)(2l2+3)/ 4Pi)  * 3J(l1,l2,l3,0,0,0) *
cbs   3J(l1,l2,l3,m1,m2,m3)
cbs   
cbs   
      implicit double precision(a-h,o-z)
      double precision inv4pi
      common /pifact/ pi,root8ovpi 
cbs   (4*PI)**-1 
      inv4pi=0.25d0/pi                                   
c     write(6,'(A10,6I5)') 'couple3J: ',l1,l2,l3,m1,m2,m3
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
c     write(6,*) 'fac1,fac2,fac3 in couple3J ',fac1,fac2,fac3 
      couple3J=fac1*fac2*fac3
      return 
      end 
