!> @file
!> Calculate canonical MP2 energy using Tensor hyper contraction

module full_ls_thc_rimp2Mod

#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
  use decmpi_module, only: mpi_bcast_fullmolecule
  use lsmpi_op
#endif
  use fundamental
  use precision
  use typedeftype
  use typedef
  use dec_typedef_module
  use matrix_module
  use matrix_operations
  use memory_handling
  !  DEC DEPENDENCIES (within deccc directory)   
  !  *****************************************
  use dec_fragment_utils
  use rimp2_module
  use full_molecule
  use THC_util
  use files
  public :: full_canonical_ls_thc_rimp2
  private

contains
  !> \brief Calculate canonical MP2 energy using Tensor hyper contraction
  !> 
  !> LS-THC-RIMP2
  !> 
  !> \author Thomas Kjaergaard
  !> \date October 2014
subroutine full_canonical_ls_thc_rimp2(MyMolecule,MyLsitem,mp2_energy)
  implicit none
  !> Full molecule info
  type(fullmolecule), intent(inout) :: MyMolecule
  !> Lsitem structure
  type(lsitem), intent(inout) :: mylsitem
  !> Canonical MP2 correlation energy
  real(realk),intent(inout) :: mp2_energy    
  !
  integer :: nbasis,nocc,nvirt,naux,noccfull,mynum,numnodes,nb,nbasisAux
  integer :: nAtoms,lupri,ngrid,Iprint,NBA,M,N,K,I,J,A,B,ALPHA,P,Q,offset
  real(realk) :: TS,TE,epsilon,TMP,TS2,TE2
  real(realk),pointer :: EpsOcc(:),EpsVirt(:),X(:,:),S(:,:),S_inv(:,:),Zpq(:,:)
  real(realk),pointer :: Calpha(:),ABdecomp(:,:),Mmat(:,:),Epq(:,:),XO(:,:)
  real(realk),pointer :: XV(:,:),TZpq(:,:),IntTHC(:,:,:,:),IntRI(:,:,:,:)
  real(realk),pointer :: SC(:,:),Identity(:,:)
  real(realk),pointer :: Tvirt(:,:),Tocc(:,:),BmatPR(:,:,:),AmatPR(:,:,:)
  real(realk),pointer :: ABmatPR(:,:,:),ZCDmat(:,:,:),IntAO(:,:,:,:)
!  integer :: BETA,GAMMA,DELTA
  logical :: master,FORCEPRINT,CollaborateWithSlaves,ABdecompCreate,WriteToDisk
  character :: intspec(5)
  integer,parameter :: nLaplace=10
  real(realk),parameter,dimension(10) :: LaplaceAmp = (/ -0.003431, &
       & -0.023534, -0.088984, -0.275603, -0.757121, -1.906218, -4.485611, &
       & -10.008000, -21.491075, -45.877205 /)
  real(realk),parameter,dimension(10) :: LaplaceW = (/ 0.009348, &
       & 0.035196, 0.107559, 0.293035, 0.729094, 1.690608, 3.709278, &
       & 7.810243, 16.172017, 35.929402 /)

  CALL LSTIMER('START ',TS,TE,DECINFO%OUTPUT)
  CALL LSTIMER('START ',TS2,TE2,DECINFO%OUTPUT)
  mp2_energy = 0.0E0_realk
  WriteToDisk = DECinfo%THCDUMP
  !sanity check
  if(.NOT.DECinfo%use_canonical) then
     call lsquit('Error: full_canonical_ls_thc_rimp2 require canonical Orbitals',-1)
  endif

  ! Init stuff
  ! **********
  nbasis = MyMolecule%nbasis
  nb = nbasis
  nocc   = MyMolecule%nocc
  nvirt  = MyMolecule%nvirt
  noccfull = nocc
  nAtoms = MyMolecule%nAtoms
  LUPRI = DECinfo%output
  nBasisAux = MyMolecule%nauxbasis

  master = .TRUE.
  FORCEPRINT=.FALSE.
  CollaborateWithSlaves = .FALSE.
  ABdecompCreate = .TRUE.
  mynum = 0
  numnodes = 1

  Iprint = 0
  print*,'DECinfo%THCradint',DECinfo%THCradint
  print*,'DECinfo%THCangint',DECinfo%THCangint
  print*,'DECinfo%THCHRDNES',DECinfo%THCHRDNES
  print*,'DECinfo%THCNOPRUN',DECinfo%THCNOPRUN
  print*,'DECinfo%THCTURBO',DECinfo%THCTURBO
  print*,'DECinfo%THCRADIALGRID',DECinfo%THCRADIALGRID
  print*,'DECinfo%THCZdependenMaxAng',DECinfo%THCZdependenMaxAng
  print*,'DECinfo%THCPARTITIONING',DECinfo%THCPARTITIONING
  print*,'DECinfo%THC_MIN_RAD_PT',DECinfo%THC_MIN_RAD_PT
  !Step 1 Build the grid and return the number of gridpoints
  call Get_THC_AO_grid_ngrid(DECinfo%output,iprint,mylsitem%setting,nb,ngrid,&
       & DECinfo%THCradint,DECinfo%THCangint,DECinfo%THCHRDNES,&
       & DECinfo%THCNOPRUN,DECinfo%THCTURBO,DECinfo%THCRADIALGRID,&
       & DECinfo%THCZdependenMaxAng,DECinfo%THCPARTITIONING,&
       & DECinfo%THC_MIN_RAD_PT)
  print*,'THC ngrid',ngrid 
  CALL LSTIMER('THC:Get_THC_AO_grid_ngrid',TS2,TE2,DECINFO%OUTPUT)

  !Step 2a Get the Atomic Orbital function values on the gridpoints
  call mem_alloc(X,ngrid,nbasis)
  !remove grid points that do not contribute to any functions with more than 
  !THC_GRID_THRESHOLD
  call Get_THC_AO_grid_X(DECinfo%output,iprint,mylsitem%setting,nb,ngrid,X)

  CALL LSTIMER('THC:Get_THC_AO_grid_X',TS2,TE2,DECINFO%OUTPUT)

  IF(WriteToDisk)THEN
     !Keep things in the AO basis
     !S(P,Q) = (X(P,mu)*X(Q,mu))(X(P,nu)*X(Q,nu))
     call mem_alloc(S,ngrid,ngrid)
     call Get_THC_grid_overlap(S,X,nbasis,ngrid)
     CALL LSTIMER('THC:OVERLAP2',TS2,TE2,DECINFO%OUTPUT)     
  ELSE
     ! Change X to MO  
     !step 2b XV(P,a)=X(P,mu)*Cv(mu,a)
     call mem_alloc(XV,ngrid,nvirt)
     M = ngrid   !rows of Output Matrix
     N = nvirt   !columns of Output Matrix
     K = nbasis  !summation dimension
     call dgemm('N','N',M,N,K,1.0E0_realk,X,M,MyMolecule%Cv%elm2,K,0.0E0_realk,XV,M)
     
     CALL LSTIMER('THC:DGEMM1',TS2,TE2,DECINFO%OUTPUT)

     !step 2c XO(P,i)=X(P,mu)*Co(mu,i)
     call mem_alloc(XO,ngrid,nocc)
     M = ngrid   !rows of Output Matrix
     N = nocc    !columns of Output Matrix
     K = nbasis  !summation dimension
     call dgemm('N','N',M,N,K,1.0E0_realk,X,M,MyMolecule%Co%elm2,K,0.0E0_realk,XO,M)
     CALL LSTIMER('THC:DGEMM2',TS2,TE2,DECINFO%OUTPUT)
     call mem_dealloc(X)

     !Step 3 Construct the Grid Overlap Matrix Eq. 28 of JCP 137, 224106
     call mem_alloc(S,ngrid,ngrid)    
     call Get_THC_grid_overlap2(S,XO,XV,nocc,nvirt,ngrid)
     CALL LSTIMER('THC:OVERLAP2',TS2,TE2,DECINFO%OUTPUT)
  ENDIF

  !Step 4 Construct the Inverse Grid Overlap Matrix
  call mem_alloc(S_inv,ngrid,ngrid)
  epsilon = 1.0E-10_realk
  call Get_THC_grid_overlap_inv(S,S_inv,ngrid,epsilon)
  CALL LSTIMER('THC:OVERLAP_INV',TS2,TE2,DECINFO%OUTPUT)
  call mem_dealloc(S)

  IF(WriteToDisk)THEN
     !step 5a construct (alpha|mu nu) 
     call mem_alloc(Identity,nbasis,nbasis)
     call ls_dzero(Identity,nbasis*nbasis)
     DO I=1,nbasis
        Identity(I,I)=1.0E0_realk
     ENDDO     
     intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
     intspec(2) = 'R' !Regular AO basis function on center 3
     intspec(3) = 'R' !Regular AO basis function on center 4
     intspec(4) = 'C' !Coulomb Operator
     intspec(5) = 'C' !Coulomb Operator
     call mem_alloc(ABdecomp,nbasisAux,nbasisAux)
     ABdecompCreate = .TRUE.
     call Build_CalphaMO2(Mylsitem,master,nbasis,nbasis,nbasisAux,LUPRI,&
          & FORCEPRINT,CollaborateWithSlaves,Identity,nbasis,&
          & Identity,nbasis,mynum,numnodes,Calpha,NBA,ABdecomp,&
          & ABdecompCreate,intspec,.FALSE.)

     call mem_dealloc(ABdecomp)
     CALL LSTIMER('THC:Build_CalphaMO2',TS2,TE2,DECINFO%OUTPUT)
  ELSE
     !step 5a construct (alpha|ai) 
     intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
     intspec(2) = 'R' !Regular AO basis function on center 3
     intspec(3) = 'R' !Regular AO basis function on center 4
     intspec(4) = 'C' !Coulomb Operator
     intspec(5) = 'C' !Coulomb Operator
     call mem_alloc(ABdecomp,nbasisAux,nbasisAux)
     ABdecompCreate = .TRUE.
     call Build_CalphaMO2(Mylsitem,master,nbasis,nbasis,nbasisAux,LUPRI,&
          & FORCEPRINT,CollaborateWithSlaves,MyMolecule%Cv%elm2,nvirt,&
          & MyMolecule%Co%elm2,nocc,mynum,numnodes,Calpha,NBA,ABdecomp,&
          & ABdecompCreate,intspec,.FALSE.)
     call mem_dealloc(ABdecomp)

     CALL LSTIMER('THC:Build_CalphaMO2',TS2,TE2,DECINFO%OUTPUT)
  ENDIF

  IF(WriteToDisk)THEN
     !step 6 construct M(P,alpha) 
     !M(P,alpha) = (alpha|mu nu)*X(P,mu)*X(P,nu)   !scaling: N*N*Nalpha*Ngrid
     call mem_alloc(Mmat,ngrid,NBA)
     CALL build_THC_MalphaP(Calpha,NBA,nbasis,nbasis,X,ngrid,X,Mmat)
     CALL LSTIMER('THC:build_THC_MalphaP',TS2,TE2,DECINFO%OUTPUT)
  ELSE
     !step 6 construct M(P,alpha) 
     !M(P,alpha) = (alpha|ai)*X(P,a)*X(P,i)   !scaling: O*V*Nalpha*Ngrid
     call mem_alloc(Mmat,ngrid,NBA)
     CALL build_THC_MalphaP(Calpha,NBA,nvirt,nocc,XV,ngrid,XO,Mmat)
     CALL LSTIMER('THC:build_THC_MalphaP',TS2,TE2,DECINFO%OUTPUT)
     call mem_dealloc(Calpha)
  ENDIF

  !step 7 construct Epq
  !E(P,Q) = M(P,alpha)*M(Q,alpha)
  M = ngrid   !rows of Output Matrix
  N = ngrid   !columns of Output Matrix
  K = NBA     !summation dimension
  call mem_alloc(Epq,ngrid,ngrid)
  call dgemm('N','T',M,N,K,1.0E0_realk,Mmat,M,Mmat,M,0.0E0_realk,Epq,M)
  CALL LSTIMER('THC:DGEMM3',TS2,TE2,DECINFO%OUTPUT)
  call mem_dealloc(Mmat)

  !step 8 construct Zpq(S,R) = S_inv(S,P)*Epq(P,Q)*S_inv(Q,R)
  M = ngrid   !rows of Output Matrix
  N = ngrid   !columns of Output Matrix
  K = ngrid   !summation dimension
  !TZpq(S,R) = S_inv(S,P)*Epq(P,Q)
  call mem_alloc(TZpq,ngrid,ngrid)
  call dgemm('N','N',M,N,K,1.0E0_realk,S_inv,M,Epq,K,0.0E0_realk,TZpq,M)
  CALL LSTIMER('THC:DGEMM4',TS2,TE2,DECINFO%OUTPUT)
  call mem_dealloc(Epq)

  !Zpq(S,R) = TZpq(S,Q)*S_inv(Q,R)
  call mem_alloc(Zpq,ngrid,ngrid)
  call dgemm('N','N',M,N,K,1.0E0_realk,TZpq,M,S_inv,K,0.0E0_realk,Zpq,M)
  CALL LSTIMER('THC:DGEMM5',TS2,TE2,DECINFO%OUTPUT)
  call mem_dealloc(TZpq)
  call mem_dealloc(S_inv)

  IF(WriteToDisk)call DumpTHCmatrices(nbasis,ngrid,X,Zpq)

  offset = 0
  call mem_alloc(EpsOcc,nocc)
  !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I) &
  !$OMP SHARED(nocc,MyMolecule,EpsOcc,offset)
  do I=1,nocc
     EpsOcc(I) = MyMolecule%oofock%elm2(I+offset,I+offset)
  enddo
  !$OMP END PARALLEL DO
  call mem_alloc(EpsVirt,nvirt)
  !$OMP PARALLEL DO DEFAULT(none) PRIVATE(A) &
  !$OMP SHARED(nvirt,MyMolecule,EpsVirt)
  do A=1,nvirt
     EpsVirt(A) = MyMolecule%vvfock%elm2(A,A)
  enddo
  !$OMP END PARALLEL DO
  mp2_energy=0.0E0_realk
  CALL LSTIMER('THC:Eps',TS2,TE2,DECINFO%OUTPUT)
  IF(WriteToDisk)THEN
     !verify
     call TestTHCRIAOintegrals(X,ZPQ,Calpha,nbasis,nocc,nvirt,ngrid,NBA,lupri)
     call mem_dealloc(Calpha)

     call mem_alloc(IntAO,nbasis,nbasis,nbasis,nbasis)
     intspec(1) = 'R'
     intspec(2) = 'R'
     intspec(3) = 'R'
     intspec(4) = 'R'
     intspec(5) = 'C'
     call II_get_4center_eri(DECINFO%OUTPUT,DECINFO%OUTPUT,mylsitem%SETTING,&
          & IntAO,nbasis,nbasis,nbasis,nbasis,intspec)     
     call TestTHCAOintegrals(X,ZPQ,IntAO,nbasis,nocc,nvirt,ngrid,lupri)
     call mem_dealloc(IntAO)

     call LS_THC_AO_RIMP2_Ecorr(nocc,nvirt,ngrid,X,Zpq,EpsOcc,EpsVirt,&
          & mp2_energy,MyMolecule%Co%elm2,MyMolecule%Cv%elm2,nbasis)
  ELSE
     call LS_THC_RIMP2_Ecorr(nocc,nvirt,ngrid,XO,XV,Zpq,EpsOcc,EpsVirt,mp2_energy)
  ENDIF
  CALL LSTIMER('THC:Ecorr',TS2,TE2,DECINFO%OUTPUT)
  print*,'LS-THC-RI-MP2 energy=',mp2_energy

!  Coulomb Contribution
!  Coulomb:X(a,P)*X(i,P)*Z(P,Q)*X(b,Q)*X(j,Q)*X(a,R)*X(i,R)*Z(R,S)*X(b,S)*X(j,S)
!  tau(i,l) = exp(-epsilon_I*amp_l)   !l is the laplace points
!  tau(a,l) = exp(epsilon_A*amp_l)   !l is the laplace points
!  A(P,R,l) = X(a,P)*X(a,R)*tau(a,l)
!  B(P,R,l) = X(i,P)*X(i,R)*tau(i,l)
!  C(Q,S,l) = X(b,Q)*X(b,S)*tau(b,l)
!  D(Q,S,l) = X(j,Q)*X(j,S)*tau(j,l)
!  Coulomb: weight_l*(A(P,R)*B(P,R))*Z(P,Q)*(C(Q,S)*D(Q,S))*Z(R,S)
!  AB(P,R,l) = A(P,R,l)*B(P,R,l)
!  CD(Q,S,l) = C(Q,S,l)*D(Q,S,l)
!  Coulomb: AB(P,R,l)*Z(P,Q)*CD(Q,S,l)*Z(R,S)  (Ngrid**3)

!  call mem_alloc(Tvirt,nvirt,nLaplace)
!  call BuildTvirt(Tvirt,nvirt,nLaplace,EpsVirt,LaplaceAmp)
!  call mem_alloc(AmatPR,ngrid,ngrid,nLaplace)
!  call BuildAmat(nLaplace,ngrid,nvirt,XV,Tvirt,AmatPR)
!  call mem_dealloc(Tvirt)

!  call mem_alloc(Tocc,nocc,nLaplace)
!  call BuildTocc(Tocc,nocc,nLaplace,EpsOcc,LaplaceAmp)
!  call mem_alloc(BmatPR,ngrid,ngrid,nLaplace)
!  call BuildAmat(nLaplace,ngrid,nocc,XO,Tocc,BmatPR)
!  call mem_dealloc(Tocc)

!  call mem_alloc(ABmatPR,ngrid,ngrid,nLaplace)
!  call BuildABmat(nLaplace,ngrid,ABmatPR,AmatPR,BmatPR)
!  call mem_dealloc(BmatPR)
!  call mem_dealloc(AmatPR)

  !  Coulomb: AB(P,R,L)*Z(P,Q)*CD(Q,S,L)*Z(R,S)  (Ngrid**3)
!  M = ngrid            !rows of Output Matrix
!  N = ngrid*nLaplace   !columns of Output Matrix
!  K = ngrid            !summation dimension
  !ZCD(P,S,L) = Z(P,Q)*CD(Q,S,L)
!  call mem_alloc(ZCDmat,ngrid,ngrid,nLaplace)
!  call dgemm('N','N',M,N,K,1.0E0_realk,Zpq,M,ABmatPR,K,0.0E0_realk,ZCDmat,M)
!  call LS_THC_RIMP2_LaplaceEcorr(ngrid,nLaplace,Zpq,ABmatPR,ZCDmat,LaplaceW,mp2_energy)
!  call mem_dealloc(ZCDmat)
!  call mem_dealloc(ABmatPR)
!  print*,'LS-THC-RI-MP2 Laplace Coulomb energy=',mp2_energy

!  Exchange:X(a,P)*X(i,P)*Z(P,Q)*X(b,Q)*X(j,Q)*X(b,R)*X(i,R)*Z(R,S)*X(a,S)*X(j,S)

!  A(j,b,Q) = X(j,Q)*X(b,Q)                               (Ngrid*Nocc*Nvirt)
!  B(j,b,P) = A(j,b,Q)*Z(P,Q)                             (Ngrid*2**Nocc*Nvirt)
!  C(j,a,S) = X(j,S)*X(a,S)                               (Ngrid*Nocc*Nvirt)
!  D(j,a,R) = Z(R,S)*C(j,a,S)                             (Ngrid**2*Nocc*Nvirt)
  
!  E(j,R,P) =  X(b,R)*B(j,b,P)                            (Ngrid**2*Nocc*Nvirt)   
!  F(j,R,P) =  X(a,P)*D(j,a,R)                            (Ngrid**2*Nocc*Nvirt)   
!  EF(R,P) = E(j,R,P)*F(j,R,P)                            (Ngrid**2*Nocc)   
!  G(P,R) = X(i,P)*X(i,R)                                 (Ngrid**2*Nocc)   
  
!  Exchange: G(P,R)*EF(R,P)                               (Ngrid**2)

!  MY VERSION
!  A(P,S) = X(a,P)*X(a,S)                                 (Ngrid**2*Nvirt)
!  B(P,R) = X(i,P)*X(i,R)                                 (Ngrid**2*Nocc)
!  C(Q,R) = X(b,Q)*X(b,R)                                 (Ngrid**2*Nvirt)
!  D(Q,S) = X(j,Q)*X(j,S)                                 (Ngrid**2*Nocc)
!  Exchange: A(P,S)*Z(P,Q)*B(P,R)*C(Q,R)*D(Q,S)*Z(R,S)    (Ngrid**4)

  IF(WriteToDisk)THEN
     call mem_dealloc(X)
  ELSE
     call mem_dealloc(XO)
     call mem_dealloc(XV)
  ENDIF
  call mem_dealloc(Zpq)
  call mem_dealloc(EpsOcc)
  call mem_dealloc(EpsVirt)

  CALL LSTIMER('LS_THC_RIMP2 ',TS,TE,DECINFO%OUTPUT)
end subroutine full_canonical_ls_thc_rimp2

subroutine TestTHCRIAOintegrals(X,ZPQ,Calpha,nbasis,nocc,nvirt,ngrid,NBA,lupri)
  implicit none
  integer :: nbasis,nocc,nvirt,ngrid,nba,lupri
  real(realk) :: X(ngrid,nbasis)
  real(realk) :: Zpq(ngrid,ngrid)
  real(realk) :: Calpha(NBA,nbasis,nbasis)
  !
  integer :: ALPHA,BETA,GAMMA,DELTA,P,Q
  real(realk) :: TMP,TMPRI,DIFFsq
  DIFFsq = 0.0E0_realk
  !$OMP PARALLEL DO DEFAULT(none) PRIVATE(ALPHA,BETA,GAMMA,DELTA,P,Q,&
  !$OMP TMP,TMPRI) SHARED(X,ZPQ,nbasis,nocc,nvirt,ngrid,nba,Calpha) REDUCTION(+:DIFFsq)
  DO DELTA=1,nbasis
     DO GAMMA=1,nbasis
        DO BETA=1,nbasis
           DO ALPHA=1,nbasis
              TMP=0.0E0_realk
              DO P=1,ngrid
                 DO Q=1,ngrid
                    TMP=TMP + X(P,ALPHA)*X(P,BETA)*Zpq(P,Q)*X(Q,GAMMA)*X(Q,DELTA)
                 ENDDO
              ENDDO
              TMPRI = 0.0E0_realk
              DO P=1,NBA
                 TMPRI = TMPRI + Calpha(P,ALPHA,BETA)*Calpha(P,GAMMA,DELTA)
              ENDDO
              DIFFSQ = DIFFSQ  + (TMP-TMPRI)**2
              IF(ABS(TMP-TMPRI).GT.1.0E-5_realk)THEN
                 print*,'ALPHA,BETA,GAMMA,DELTA',ALPHA,BETA,GAMMA,DELTA,'RIINT DIFF',ABS(TMP-TMPRI)
                 print*,'TMP',TMP
                 print*,'RI Int',TMPRI
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  WRITE(lupri,*)'DIFF(THC-RI) = ',SQRT(DIFFsq)
end subroutine TestTHCRIAOintegrals

subroutine TestTHCAOintegrals(X,ZPQ,IntAO,nbasis,nocc,nvirt,ngrid,lupri)
  implicit none
  integer :: nbasis,nocc,nvirt,ngrid,lupri
  real(realk) :: X(ngrid,nbasis)
  real(realk) :: Zpq(ngrid,ngrid)
  real(realk) :: IntAO(nbasis,nbasis,nbasis,nbasis)
  !
  integer :: ALPHA,BETA,GAMMA,DELTA,P,Q
  real(realk) :: TMP,DIFFsq
  DIFFsq = 0.0E0_realk
  !$OMP PARALLEL DO DEFAULT(none) PRIVATE(ALPHA,BETA,GAMMA,DELTA,P,Q,&
  !$OMP TMP) SHARED(X,ZPQ,IntAO,nbasis,nocc,nvirt,ngrid) REDUCTION(+:DIFFsq)
  DO DELTA=1,nbasis
     DO GAMMA=1,nbasis
        DO BETA=1,nbasis
           DO ALPHA=1,nbasis
              TMP=0.0E0_realk
              DO P=1,ngrid
                 DO Q=1,ngrid
                    TMP=TMP + X(P,ALPHA)*X(P,BETA)*Zpq(P,Q)*X(Q,GAMMA)*X(Q,DELTA)
                 ENDDO
              ENDDO
              DIFFSQ = DIFFSQ  + (TMP-IntAO(ALPHA,BETA,GAMMA,DELTA))**2
              !              IntAOTHC(ALPHA,BETA,GAMMA,DELTA) = TMP
              IF(ABS(TMP-IntAO(ALPHA,BETA,GAMMA,DELTA)).GT.1.0E-5_realk)THEN
                 print*,'ALPHA,BETA,GAMMA,DELTA',ALPHA,BETA,GAMMA,DELTA,'AOINT DIFF',&
                      & ABS(TMP-IntAO(ALPHA,BETA,GAMMA,DELTA))
                 print*,'TMP',TMP
                 print*,'IntAO(ALPHA,BETA,GAMMA,DELTA)',IntAO(ALPHA,BETA,GAMMA,DELTA)
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  WRITE(lupri,*)'DIFF(THC-AOINT) = ',SQRT(DIFFsq)
end subroutine TestTHCAOintegrals
!  Coulomb: AB(P,R,l)*ZCD(P,S,l)*Z(R,S)  (Ngrid**3)
subroutine LS_THC_RIMP2_LaplaceEcorr(ngrid,nLaplace,Zpq,ABmatPR,ZCDmat,LaplaceW,mp2_energy)
  implicit none
  integer,intent(in) ::  ngrid,nLaplace
  real(realk),intent(in) :: Zpq(ngrid,ngrid),ABmatPR(ngrid,ngrid,nLaplace)
  real(realk),intent(in) :: ZCDmat(ngrid,ngrid,nLaplace)
  real(realk),intent(in) :: LaplaceW(nLaplace)
  real(realk),intent(inout) :: mp2_energy
  !
  integer :: P,R,S,L
  real(realk) :: E,TMPZ
  E = 0.0E0_realk  
  !$OMP PARALLEL DO DEFAULT(none) PRIVATE(P,R,S,L,TMPZ) SHARED(ngrid,&
  !$OMP nLaplace,ABmatPR,ZCDmat,Zpq) REDUCTION(+:E)
  DO L=1,nLaplace
     DO S=1,ngrid
        DO R=1,ngrid
           TMPZ=Zpq(R,S)
           DO P=1,ngrid
              E = E + ABmatPR(P,R,L)*ZCDmat(P,S,L)*TMPZ
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  mp2_energy = E
end subroutine LS_THC_RIMP2_LaplaceEcorr

subroutine BuildABmat(nLaplace,ngrid,ABmatPR,AmatPR,BmatPR)
  implicit none
  integer :: nLaplace,ngrid
  real(realk),intent(in) :: AmatPR(ngrid*ngrid*nLaplace)
  real(realk),intent(in) :: BmatPR(ngrid*ngrid*nLaplace)
  real(realk),intent(inout) :: ABmatPR(ngrid*ngrid*nLaplace)
  !
  integer :: I
  !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I) SHARED(ngrid,&
  !$OMP nLaplace,AmatPR,BmatPR,ABmatPR)
  DO I=1,ngrid*ngrid*nLaplace
     ABmatPR(I) = AmatPR(I)*BmatPR(I)
  ENDDO
  !$OMP END PARALLEL DO
end subroutine BuildABmat

subroutine BuildAmat(nLaplace,ngrid,nvirt,XV,Tvirt,AmatPR)
  implicit none
  integer :: nLaplace,ngrid,nvirt
  real(realk),intent(in) :: XV(nvirt,ngrid),Tvirt(nvirt,nLaplace)
  real(realk),intent(inout) :: AmatPR(ngrid,ngrid,nLaplace)
  !
  integer :: l,R,P,A
  real(realk) :: TMP
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(l,&
  !$OMP R,P,A,TMP) SHARED(nvirt,nLaplace,ngrid,Tvirt,XV,AmatPR)
  do l=1,nLaplace
     do R=1,ngrid
        do P=1,ngrid
           TMP = 0.0E0_realk
           do A=1,nvirt
              TMP = TMP + XV(A,P)*XV(A,R)*Tvirt(A,L)
           enddo
           AmatPR(P,R,l) = TMP
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO
end subroutine BuildAmat

subroutine BuildTvirt(Tvirt,nvirt,nLaplace,EpsVirt,LaplaceAmp)
  implicit none
  integer,intent(in) :: nvirt,nLaplace
  real(realk),intent(in) :: EpsVirt(nvirt),LaplaceAmp(nLaplace)
  real(realk),intent(inout) :: Tvirt(nvirt,nLaplace)
  !
  integer :: L,A
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(L,&
  !$OMP A) SHARED(nvirt,nLaplace,EpsVirt,LaplaceAmp,Tvirt)
  do L=1,nLaplace
     do A=1,nvirt
        Tvirt(A,L) = exp(EpsVirt(A)*LaplaceAmp(L))
     enddo
  enddo
  !$OMP END PARALLEL DO
end subroutine BuildTvirt

subroutine BuildTocc(Tocc,nocc,nLaplace,EpsOcc,LaplaceAmp)
  implicit none
  integer,intent(in) :: nocc,nLaplace
  real(realk),intent(in) :: EpsOcc(nocc),LaplaceAmp(nLaplace)
  real(realk),intent(inout) :: Tocc(nocc,nLaplace)
  !
  integer :: L,I
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(L,&
  !$OMP I) SHARED(nocc,nLaplace,EpsOcc,LaplaceAmp,Tocc)
  do L=1,nLaplace
     do I=1,nocc
        Tocc(I,L) = exp(-EpsOcc(I)*LaplaceAmp(L))
     enddo
  enddo
  !$OMP END PARALLEL DO
end subroutine BuildTocc

subroutine DumpTHCmatrices(nbasis,ngrid,X,Zpq)
  implicit none
  integer,intent(in) :: ngrid,nbasis
  real(realk),intent(in) :: X(ngrid,nbasis)
  real(realk),intent(in) :: Zpq(ngrid,ngrid)
  !local variables
  integer :: LU
  logical :: save_access_stream
  LU = -1
  save_access_stream = access_stream
  access_stream = .TRUE.
  call lsopen(LU,'THC.restart','UNKNOWN','UNFORMATTED')
  WRITE(LU) ngrid
  WRITE(LU) nbasis
  WRITE(LU) X
  WRITE(LU) Zpq
  call lsclose(LU,'KEEP')
  access_stream = save_access_stream
end subroutine DumpTHCmatrices

subroutine LS_THC_RIMP2_Ecorr(nocc,nvirt,ngrid,XO,XV,Zpq,EpsOcc,EpsVirt,Ecorr)
  implicit none
  integer,intent(in) :: nocc,nvirt,ngrid
  real(realk),intent(in) :: XO(ngrid,nocc),XV(ngrid,nvirt)
  real(realk),intent(in) :: Zpq(ngrid,ngrid),EpsOcc(nocc),EpsVirt(nvirt)
  real(realk),intent(inout) :: Ecorr
  !
  integer :: A,B,I,J,P,Q
  real(realk) :: gaibj,gbiaj,TMP,EpsIAB,EcorrJ
  real(realk),pointer :: Z(:,:,:)
  call mem_alloc(Z,ngrid,nvirt,nocc)
  Ecorr=0.0E0_realk
  EcorrJ=0.0E0_realk
!  DO A=1,nvirt
!     DO B=1,nvirt
!        DO I=1,nocc
!           DO J=1,nocc
!              gaibj = 0.0E0_realk
!              DO P=1,ngrid
!                 DO Q=1,ngrid
!                    gaibj = gaibj + XV(P,A)*XO(P,I)*Zpq(P,Q)*XV(Q,B)*XO(Q,J)
!                 ENDDO
!              ENDDO
!              gbiaj = 0.0E0_realk
!              DO P=1,ngrid
!                 DO Q=1,ngrid
!                    gbiaj = gbiaj + XV(P,B)*XO(P,I)*Zpq(P,Q)*XV(Q,A)*XO(Q,J)
!                 ENDDO
!              ENDDO
!              Ecorr = Ecorr + gaibj*(2.0E0_realk*gaibj-gbiaj)/(EpsOcc(I)+EpsOcc(J)-EpsVirt(A)-EpsVirt(B))
!           ENDDO
!        ENDDO
!     ENDDO
!  ENDDO
!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(3) PRIVATE(I,A,Q,P,&
!$OMP TMP) SHARED(XO,XV,Z,Zpq,nocc,nvirt,ngrid)
  DO I=1,nocc
   DO A=1,nvirt
    DO Q=1,ngrid
     TMP = 0.0E0_realk
     DO P=1,ngrid
      TMP = TMP + XV(P,A)*XO(P,I)*Zpq(P,Q)
     ENDDO
     Z(Q,A,I) = TMP
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(3) PRIVATE(I,A,B,J,Q,EpsIAB,&
!$OMP gaibj,gbiaj) SHARED(XO,XV,Z,EpsOcc,EpsVirt,nocc,nvirt,&
!$OMP ngrid) REDUCTION(+:Ecorr,EcorrJ) 
  DO I=1,nocc
   DO A=1,nvirt
    DO B=1,nvirt
     EpsIAB = EpsOcc(I)-EpsVirt(A)-EpsVirt(B)
     DO J=1,nocc
      gaibj = 0.0E0_realk
      DO Q=1,ngrid
         gaibj = gaibj + Z(Q,A,I)*XV(Q,B)*XO(Q,J)
      ENDDO
      gbiaj = 0.0E0_realk
      DO Q=1,ngrid
         gbiaj = gbiaj + Z(Q,B,I)*XV(Q,A)*XO(Q,J)
      ENDDO
      Ecorr = Ecorr + gaibj*(2.0E0_realk*gaibj-gbiaj)/(EpsOcc(J)+EpsIAB)
      EcorrJ = EcorrJ + gaibj*(2.0E0_realk*gaibj)/(EpsOcc(J)+EpsIAB)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO
  print*,'EcorrJ',EcorrJ
  call mem_dealloc(Z)
end subroutine LS_THC_RIMP2_Ecorr

subroutine LS_THC_AO_RIMP2_Ecorr(nocc,nvirt,ngrid,X,Zpq,EpsOcc,EpsVirt,Ecorr,Cocc,Cvirt,nbasis)
  implicit none
  integer,intent(in) :: nocc,nvirt,ngrid,nbasis
  real(realk),intent(in) :: X(ngrid,nbasis),Cvirt(nbasis,nvirt),Cocc(nbasis,nocc)
  real(realk),intent(in) :: Zpq(ngrid,ngrid),EpsOcc(nocc),EpsVirt(nvirt)
  real(realk),intent(inout) :: Ecorr
  !
  integer :: A,B,I,J,P,Q,M,N,K
  real(realk) :: gaibj,gbiaj,TMP,EpsIAB,EcorrJ
  real(realk),pointer :: Z(:,:,:),Z1(:,:,:),Z2(:,:,:),XO(:,:),XV(:,:)
  call mem_alloc(Z,ngrid,nbasis,nbasis)
  Ecorr=0.0E0_realk
  EcorrJ=0.0E0_realk
!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(3) PRIVATE(I,A,Q,P,&
!$OMP TMP) SHARED(X,Z,Zpq,nbasis,ngrid)
  DO I=1,nbasis
   DO A=1,nbasis
    DO Q=1,ngrid
     TMP = 0.0E0_realk
     DO P=1,ngrid
      TMP = TMP + X(P,A)*X(P,I)*Zpq(P,Q)
     ENDDO
     Z(Q,A,I) = TMP
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO
  call mem_alloc(Z1,ngrid,nbasis,nocc)
  M = ngrid*nbasis   !rows of Output Matrix
  N = nocc   !columns of Output Matrix
  K = nbasis  !summation dimension
  call dgemm('N','N',M,N,K,1.0E0_realk,Z,M,Cocc,K,0.0E0_realk,Z1,M)
  call mem_dealloc(Z)
  call mem_alloc(Z2,ngrid,nvirt,nocc)
!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(2) PRIVATE(I,A,Q,P,&
!$OMP TMP) SHARED(X,Z1,Z2,nocc,nvirt,ngrid,nbasis,Cvirt)
  DO I=1,nocc
     DO A=1,nvirt
        DO Q=1,ngrid
           Z2(Q,A,I) = Z1(Q,1,I)*Cvirt(1,A)
        ENDDO
        DO P=2,nbasis
           TMP = Cvirt(P,A)
           DO Q=1,ngrid
              Z2(Q,A,I) = Z2(Q,A,I)+Z1(Q,P,I)*TMP
           ENDDO
        ENDDO
     ENDDO
  ENDDO
!$OMP END PARALLEL DO
  call mem_dealloc(Z1)

  !XO(ngrid,nocc) = X(ngrid,nbasis)*Cocc(nbasis,nocc)
  call mem_alloc(XO,ngrid,nocc)
  M = ngrid   !rows of Output Matrix
  N = nocc    !columns of Output Matrix
  K = nbasis  !summation dimension
  call dgemm('N','N',M,N,K,1.0E0_realk,X,M,Cocc,K,0.0E0_realk,XO,M)

  !XV(ngrid,nvirt) = X(ngrid,nbasis)*Cvirt(nbasis,nvirt)
  call mem_alloc(XV,ngrid,nvirt)
  M = ngrid   !rows of Output Matrix
  N = nvirt   !columns of Output Matrix
  K = nbasis  !summation dimension
  call dgemm('N','N',M,N,K,1.0E0_realk,X,M,Cvirt,K,0.0E0_realk,XV,M)

!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(3) PRIVATE(I,A,B,J,Q,EpsIAB,&
!$OMP gaibj,gbiaj) SHARED(XO,XV,Z,EpsOcc,EpsVirt,nocc,nvirt,&
!$OMP ngrid,Z2) REDUCTION(+:Ecorr,EcorrJ) 
  DO I=1,nocc
   DO A=1,nvirt
    DO B=1,nvirt
     EpsIAB = EpsOcc(I)-EpsVirt(A)-EpsVirt(B)
     DO J=1,nocc
      gaibj = 0.0E0_realk
      DO Q=1,ngrid
         gaibj = gaibj + Z2(Q,A,I)*XV(Q,B)*XO(Q,J)
      ENDDO
      gbiaj = 0.0E0_realk
      DO Q=1,ngrid
         gbiaj = gbiaj + Z2(Q,B,I)*XV(Q,A)*XO(Q,J)
      ENDDO
      Ecorr = Ecorr + gaibj*(2.0E0_realk*gaibj-gbiaj)/(EpsOcc(J)+EpsIAB)
      EcorrJ = EcorrJ + gaibj*(2.0E0_realk*gaibj)/(EpsOcc(J)+EpsIAB)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO
  call mem_dealloc(XO)
  call mem_dealloc(XV)
  print*,'EcorrJ',EcorrJ
  call mem_dealloc(Z2)
end subroutine LS_THC_AO_RIMP2_Ecorr

subroutine build_THC_MalphaP(Calpha,NBA,nvirt,nocc,XV,ngrid,XO,M)
  implicit none
  integer,intent(in) :: NBA,nvirt,nocc,ngrid
  real(realk),intent(in) :: Calpha(NBA,nvirt,nocc),XO(ngrid,nocc)
  real(realk),intent(in) :: XV(ngrid,nvirt)
  real(realk),intent(inout) :: M(ngrid,NBA)
  !
  integer :: I,A,ALPHA,K
  real(realk) :: TI,TMP
!$OMP PARALLEL DO DEFAULT(none) COLLAPSE(2) PRIVATE(I,A,ALPHA,K,&
!$OMP TI,TMP) SHARED(NBA,nvirt,nocc,ngrid,M,XV,XO,Calpha)
  DO ALPHA=1,NBA
     DO K=1,ngrid
        TMP = 0.0E0_realk
        DO I=1,nocc             
           TI = XO(K,I)
           DO A=1,nvirt
              TMP = TMP + Calpha(ALPHA,A,I)*TI*XV(K,A)
           ENDDO
        ENDDO
        M(K,ALPHA) = TMP
     ENDDO
  ENDDO
!$OMP END PARALLEL DO
end subroutine build_THC_MalphaP

end module full_ls_thc_rimp2Mod
