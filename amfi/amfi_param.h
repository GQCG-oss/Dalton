c --- file: amfi_param.h ---
c
c     Parameters and common blocks for the AMFI module written by
c     Bernd Schimmelpfennig (bs)
c
cbs   
cbs   this include files hold a lot of dimensioning parameters and arrays for 
cbs   exponents, contraction coefficients and integrals. 
cbs   I hope, most of the names are selfexplaining...
cbs   
cbs   All those parameters are constructed to be blown up without any problem
cbs   (if there is sufficient memory) except Lmax which is limited to 4 (g-functions)
cbs   
cbs   ###############################################################################
cbs   ################ parameter block ##############################################
cbs   ###############################################################################
      parameter (Lanz=Lmax+1)   ! number of angular momenta                       
cbs   ###############################################################################
cbs   ################ parameter block ##############################################
cbs   ###############################################################################
cbs   overlap of normalized functions
      double precision   normovlp
      common  /normovl/  normovlp(MxprimL,MxprimL,0:Lmax),
     *OVLPinv(MxprimL,MxprimL,0:Lmax),
     *rootOVLP(MxprimL,MxprimL,0:Lmax),
     *rootOVLPinv(MxprimL,MxprimL,0:Lmax),
     *scratchinv((MxprimL*MxprimL+MxprimL)/2) 
cbs   defining a big array with enough space for all modified contraction coefficients
cbs   for each l-value there are five blocks of size (nprimit(l),ncontrac(l)) 
cbs   the original contraction coefficients (for normalized functions) 
cbs   and four modified blocks  depending on different kinematic factors and included exponents
c
cbs   the iadd* arrays hold the start adresses of the the contraction coefficients for each l-value 
      common /contco/ 
     *contrarray((Lmax+1)*5*MxcontL*MxprimL), 
     *iaddori(0:Lmax),iaddtyp1(0:Lmax),iaddtyp2(0:Lmax),
     *iaddtyp3(0:Lmax),iaddtyp4(0:Lmax)
cbs   the exponents 
      common /expo/ exponents(MxprimL,0:Lmax)
cbs   the numbers of primitive and contracted functions for each l-value   
      common /dims/ nprimit(0:Lmax),ncontrac(0:Lmax),nprimit_keep,
     *ncontrac_keep 
cbs   scratch should explain itself ...........
      common /scratch/ 
     *scratch4(4*MxprimL*MxprimL),
     *mcombina(2,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *cntscrtch(MxprimL,MxcontL,0:Lmax) 
      common /contint/ 
     *Lfirst(4),
     *Llast(4),Lblocks(4),Lstarter(4),
     *nblock,Lvalues(4)
cbs   cont4 will keep 4 blocks (label(i)) of structure
cbs      (ncontrac(l1),ncontrac(l2),ncontrac(l3),ncontrac(l4),(Llast(i)-Lfirst(i))/2+1)
cbs                                                               or 0   if no L-value at all 
cbs                                                               = Lblocks(i) 
cbs   for each l1,l2,l3,l4-block 
cbs   Lfirst(i,j) gives the first L-value, for which radial integrals are calculated 
cbs   for type i and l1,l2,l3,l4 - Integral block.
cbs   Llast(i,j) gives the last L-value
cbs   Lblocks gives the number of L-values
cbs   Lstarter gives the adress of each integral block on cont4
cbs   the following block contains a lot stuff for calculating the kinematic factors
      common /diagonalize/ TKIN(MxprimL*MxprimL),evec(MxprimL*MxprimL),
     *eval(MxprimL),Energy(MxprimL),type1(MxprimL),type2(MxprimL)
cbs   
cbs   some factors that appear a lot of times 
cbs   in the angular factors 
      common /prefs/ preroots(2,0:Lmax),clebsch(3,2,-Lmax:Lmax,0:Lmax)
cbs   common block with the cartesian integrals 
      common /cartint/ 
     *mcombcart(2,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *istartX,istartY,istartZ,
     *isignM(-Lmax:Lmax)
cbs   the following array includes information about 
cbs   odd powers of  x y z in the real harmonics 
cbs   this is used to check whether integrals in the 
cbs   cartesian representation appear
      common  ipowxyz(3,-Lmax:Lmax,0:Lmax)
      logical icheckxy,icheckz
cbs   preXYZ holds some prefactors which are used regulary 
      common /preXYZ/ 
     *preXZ(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *preY(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *icheckxy(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),        
     *icheckz(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),        
     *interxyz(16,0:Lmax,0:Lmax,0:Lmax,0:Lmax),       
     *isgnprod(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
c /onepar/: one-particle integrals 
cbs   onecart* arrays for one-electron-integrals: 
cbs   1. index: number of first contracted function
cbs   2. index: number of second contracted function
cbs   3. index: pointer(m1,m2)    m1< m2 otherwise change sign of integral
cbs   4. index: L-value  
      common /onepar/
     *oneoverR3((MxprimL*MxprimL+MxprimL)/2,Lmax),
     *onecontr(mxcontL,MxcontL,-Lmax:Lmax,3,Lmax),
     *dummyone(MxcontL*MxcontL),
     *onecartX(mxcontL,MxcontL,
     *(Lmax+Lmax+1)*(Lmax+1),Lmax),
     *onecartY(mxcontL,MxcontL,
     *(Lmax+Lmax+1)*(Lmax+1),Lmax),
     *onecartZ(mxcontL,MxcontL,
     *(Lmax+Lmax+1)*(Lmax+1),Lmax),
     *onecartallX((Mxcart*Mxcart+Mxcart)/2),
     *onecartallY((Mxcart*Mxcart+Mxcart)/2),
     *onecartallZ((Mxcart*Mxcart+Mxcart)/2),
     *onecartallX2((Mxcart*Mxcart+Mxcart)/2),
     *onecartallY2((Mxcart*Mxcart+Mxcart)/2),
     *onecartallZ2((Mxcart*Mxcart+Mxcart)/2),
     *onecartallX3((Mxcart*Mxcart+Mxcart)/2),
     *onecartallY3((Mxcart*Mxcart+Mxcart)/2),
     *onecartallZ3((Mxcart*Mxcart+Mxcart)/2)
cbs  powexp holds powers of exponents and meam exponents 
cbs  coulovp holds overlap of functions with shifted   
cbs  l-values   
      common /coulpow/ powexp(MxprimL,MxprimL,0:Lmax,
     *0:Lmax,0:(Lmax+Lmax+5)),
     *coulovlp(MxprimL,MxprimL,-1:1,-1:1,
     *0:Lmax,0:Lmax)
cbs   express AOs in contracted functions 
cbs   first index: number of contracted function                                    
cbs   second index: number of AO                                                     
cbs   third index: L-value 
      common /AOincont/ AOcoeffs(MxcontL,MxcontL,0:Lmax),
     *occup(MxcontL,0:Lmax),noccorb(0:Lmax)  
cbs   occupation numbers 
cbs   first index: number of AO                    
cbs   second index: L-value 
      common /corelist/ icore(0:Lmax),ikeeporb,ikeeplist(Mxcart),
     *nrtofiperIR(8)   
c --- end of amfi_param.h ---
