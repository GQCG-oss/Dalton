      subroutine prefac(Lmax,preroots,clebsch)
      implicit double precision (a-h,o-z)   
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
     

