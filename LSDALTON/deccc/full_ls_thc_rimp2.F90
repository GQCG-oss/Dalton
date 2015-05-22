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
  logical :: master,FORCEPRINT,CollaborateWithSlaves,ABdecompCreate,WriteToDisk
  character :: intspec(5)

  CALL LSTIMER('START ',TS,TE,DECINFO%OUTPUT)
  CALL LSTIMER('START ',TS2,TE2,DECINFO%OUTPUT)
  mp2_energy = 0.0E0_realk
  WriteToDisk = .TRUE.  
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
     !Step 3 Construct the Grid Overlap Matrix Eq. 28 of JCP 137, 224106
     call mem_alloc(S,ngrid,ngrid)    
     call Get_THC_grid_overlap2(S,X,X,nbasis,nbasis,ngrid)
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
     call mem_dealloc(Calpha)
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
!     usefull for testing but too slow to do
!     call LS_THC_RIMP2_EcorrAO(nocc,nvirt,ngrid,nbasis,X,Zpq,EpsOcc,EpsVirt,&
!          & MyMolecule%Cv%elm2,MyMolecule%Co%elm2,mp2_energy)     
  ELSE
     call LS_THC_RIMP2_Ecorr(nocc,nvirt,ngrid,XO,XV,Zpq,EpsOcc,EpsVirt,mp2_energy)
  ENDIF
  CALL LSTIMER('THC:Ecorr',TS2,TE2,DECINFO%OUTPUT)
  print*,'LS-THC-RI-MP2 energy=',mp2_energy
  call mem_dealloc(XO)
  call mem_dealloc(XV)
  call mem_dealloc(Zpq)
  call mem_dealloc(EpsOcc)
  call mem_dealloc(EpsVirt)

!  epsilon missing due to lagrange ???? 

!  Coulomb:X(a,P)*X(i,P)*Z(P,Q)*X(b,Q)*X(j,Q)*X(a,R)*X(i,R)*Z(R,S)*X(b,S)*X(j,S)
  
!  A(P,R) = X(a,P)*X(a,R)
!  B(P,R) = X(i,P)*X(i,R)
!  C(Q,S) = X(b,Q)*X(b,S)
!  D(Q,S) = X(j,Q)*X(j,S)
!  Coulomb: (A(P,R)*B(P,R))*Z(P,Q)*(C(Q,S)*D(Q,S))*Z(R,S)/(epsilon_I+epsilon_J-epsilon_A+epsilon_B)
!  AB(P,R) = A(P,R)*B(P,R)
!  CD(Q,S) = C(Q,S)*D(Q,S)
!  Coulomb: AB(P,R)*Z(P,Q)*CD(Q,S)*Z(R,S)  (Ngrid**3)

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


  CALL LSTIMER('LS_THC_RIMP2 ',TS,TE,DECINFO%OUTPUT)
end subroutine full_canonical_ls_thc_rimp2

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

subroutine LS_THC_RIMP2_EcorrAO(nocc,nvirt,ngrid,nbasis,X,Zpq,EpsOcc,EpsVirt,Cv,Co,Ecorr)
  implicit none
  integer,intent(in) :: nocc,nvirt,ngrid,nbasis
  real(realk),intent(in) :: X(ngrid,nbasis)
  real(realk),intent(in) :: Co(nbasis,nocc),Cv(nbasis,nvirt)
  real(realk),intent(in) :: Zpq(ngrid,ngrid),EpsOcc(nocc),EpsVirt(nvirt)
  real(realk),intent(inout) :: Ecorr
  !
  integer :: A,B,I,J,P,Q,ALPHA,BETA,GAMMA,DELTA,M,N,K
  real(realk) :: gaibj,gbiaj,TMP,EpsIAB,TMP2
  real(realk),pointer :: IntAO(:,:,:,:),IntALPHABETAGAMMAJ(:,:,:,:),IntABETAGAMMAJ(:,:,:,:)
  real(realk),pointer :: IntAIGAMMAJ(:,:,:,:),IntMO(:,:,:,:)
  call mem_alloc(IntAO,nbasis,nbasis,nbasis,nbasis)
  !$OMP PARALLEL DO DEFAULT(none) PRIVATE(ALPHA,BETA,GAMMA,DELTA,P,Q,&
  !$OMP TMP) SHARED(X,ZPQ,IntAO,nbasis,nocc,nvirt)
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
              IntAO(ALPHA,BETA,GAMMA,DELTA) = TMP             
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  M=nbasis*nbasis*nbasis
  N=nocc
  K=nbasis
  call mem_alloc(IntALPHABETAGAMMAJ,nbasis,nbasis,nbasis,nocc)
  call DGEMM('N','N',M,N,K,1.0E0_realk,IntAO,M,Co,K,0.0E0_realk,IntALPHABETAGAMMAJ,M)
  call mem_dealloc(IntAO)
  M=nvirt
  N=nbasis*nbasis*nocc
  K=nbasis
  call mem_alloc(IntABETAGAMMAJ,nvirt,nbasis,nbasis,nocc)
  call DGEMM('T','N',M,N,K,1.0E0_realk,Cv,K,IntALPHABETAGAMMAJ,K,0.0E0_realk,IntABETAGAMMAJ,M)
  call mem_dealloc(IntALPHABETAGAMMAJ)
  call mem_alloc(IntAIGAMMAJ,nvirt,nocc,nbasis,nocc)
  !$OMP PARALLEL DO DEFAULT(none) PRIVATE(BETA,GAMMA,A,I,&
  !$OMP TMP) SHARED(IntABETAGAMMAJ,IntAIGAMMAJ,Co,nbasis,nocc,nvirt)
  DO J=1,nocc
     DO A=1,nvirt
        DO GAMMA=1,nbasis
           DO I=1,nocc
              TMP = 0.0E0_realk
              DO BETA=1,nbasis
                 TMP = TMP + IntABETAGAMMAJ(A,BETA,GAMMA,J)*Co(I,BETA)
              ENDDO
              IntAIGAMMAJ(A,I,GAMMA,J)=TMP
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  call mem_dealloc(IntABETAGAMMAJ)
  call mem_alloc(IntMO,nvirt,nocc,nvirt,nocc)
  !$OMP PARALLEL DO DEFAULT(none) PRIVATE(A,I,B,J,GAMMA,&
  !$OMP TMP) SHARED(IntAIGAMMAJ,IntMO,Cv,nbasis,nocc,nvirt)
  DO J=1,nocc
     DO A=1,nvirt
        DO I=1,nocc
           DO B=1,nvirt
              TMP = 0.0E0_realk
              DO GAMMA=1,nbasis
                 TMP = TMP + IntAIGAMMAJ(A,I,GAMMA,J)*Cv(B,GAMMA)
              ENDDO
              IntMO(A,I,B,J)=TMP
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  call mem_dealloc(IntAIGAMMAJ)
  Ecorr=0.0E0_realk
  DO B=1,nvirt
     DO J=1,nocc
        DO I=1,nocc
           DO A=1,nvirt
              EpsIAB = EpsOcc(J)+EpsOcc(I)-EpsVirt(A)-EpsVirt(B)
              Ecorr = Ecorr + IntMO(A,I,B,J)*(2.0E0_realk*IntMO(A,I,B,J)-IntMO(A,J,B,I))/EpsIAB
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  call mem_dealloc(IntMO)
end subroutine LS_THC_RIMP2_EcorrAO

subroutine LS_THC_RIMP2_Ecorr(nocc,nvirt,ngrid,XO,XV,Zpq,EpsOcc,EpsVirt,Ecorr)
  implicit none
  integer,intent(in) :: nocc,nvirt,ngrid
  real(realk),intent(in) :: XO(ngrid,nocc),XV(ngrid,nvirt)
  real(realk),intent(in) :: Zpq(ngrid,ngrid),EpsOcc(nocc),EpsVirt(nvirt)
  real(realk),intent(inout) :: Ecorr
  !
  integer :: A,B,I,J,P,Q
  real(realk) :: gaibj,gbiaj,TMP,EpsIAB
  real(realk),pointer :: Z(:,:,:) 
  call mem_alloc(Z,ngrid,nvirt,nocc)
  Ecorr=0.0E0_realk
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
!$OMP ngrid) REDUCTION(+:Ecorr) 
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
     ENDDO
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO
  call mem_dealloc(Z)
end subroutine LS_THC_RIMP2_Ecorr

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
