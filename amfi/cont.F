      subroutine cont(L,breit,HERMIT)
cbs###########################################################################
cbs   cont prepares all required contraction coefficients for functions 
cbs   with angular momentum L
cbs###########################################################################
      implicit double precision (a-h,o-z) 
      include 'param.h'
      dimension tkintria((MxprimL*MxprimL+MxprimL)/2)
      logical*4 breit,HERMIT   
cbs   transcon transfers and normalizes contracted functions
cbs   ore more precizely the coefficients  
c
c     write(6,*) 'L,nprimit(L) ',L,nprimit(L)
c     write(6,*) 'MxprimL,MxcontL,iaddori(L),nprimit(L),ncontrac(L)',
c    *MxprimL,MxcontL,iaddori(L),nprimit(L),ncontrac(L)
c
c     write(6,*) 'call transcon '
      call transcon(cntscrtch(1,1,L),MxprimL,
     *MxcontL,normovlp(1,1,L),
     *contrarray(iaddori(L)),nprimit(L),ncontrac(L))
cbs   gentkin generates the matrix of kinetic energy  TKIN 
c     write(6,*) 'call gentkin  '
      call gentkin(L,TKIN,nprimit(L),exponents(1,L),rootOVLPinv(1,1,L))
cbs   kindiag diagonalizes TKIN 
c     write(6,*) 'call kindiag ' 
      call kindiag(TKIN,TKINTRIA,nprimit(L),evec,eval,breit)  
cbs   kinemat generates kinematic factors in 
cbs   the basis of eigenvectors   
c     write(6,*) 'call kinemat ' 
      call kinemat(L,nprimit(L),eval,type1,type2,Energy,HERMIT)
c     write(6,*) 'after kinemat ' 
      incr=nprimit(L)*ncontrac(L)
cbs   chngcont= changecont generates the contraction coeffs
cbs   including kinematic factors and even exponents as factors
c     write(6,*) 'call chngcont' 
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
