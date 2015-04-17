!> @file
!> Full molecular calculation of RI-MP2-F12

module fullrimp2f12 

  use precision
  use typedeftype!,only:lsitem
  use typedef
  use dec_typedef_module
  use matrix_module
  use matrix_operations
  use memory_handling
  use dec_fragment_utils
  use CABS_operations

  use ccintegrals

  !WARNING FOR TESTING
  use f12_routines_module
  use IntegralInterfaceMOD
  use rimp2_module
  
  public :: full_canonical_rimp2_f12

  private

contains
  !> \brief Calculate canonical MP2 energy for full molecular system
  !> keeping full AO integrals in memory. Only for testing.
  !> \author Thomas Kjaergaard
  !> \date 2015
  subroutine full_canonical_rimp2_f12(MyMolecule,MyLsitem,Dmat,mp2f12_energy)
    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> HF density matrix
    type(matrix),intent(in) :: Dmat
    !> MP2-F12 correlation energy
    real(realk),intent(inout) :: mp2f12_energy
    !> Canonical MP2 correlation energy
    real(realk) :: mp2_energy
    !local variables
    integer :: nbasis,nocc,nvirt,ncabsAO,ncabsMO,noccfull
    real(realk) :: E21,Econt(1),ExchangeF12,CoulombF12
!    real(realk) :: CMO_cabs(:,:)        
!    type(matrix) :: mCMO_cabs
    !========================================================
    !WARNING THESE SHOULD NOT BE HERE - ONLY FOR TESTING
    !========================================================
    Real(realk),pointer  :: Vijij(:,:)
    Real(realk),pointer  :: Vjiij(:,:)
    Real(realk),pointer  :: Ripjq(:,:,:,:) !nocc,nbasis,nocc,nbasis
    Real(realk),pointer  :: Gipjq(:,:,:,:) !nocc,nbasis,nocc,nbasis
    Real(realk),pointer  :: Rimjc(:,:,:,:) !nocc,noccfull,nocc,ncabsMO
    Real(realk),pointer  :: Gimjc(:,:,:,:) !nocc,noccfull,nocc,ncabsMO
    Real(realk),pointer  :: Fijkl(:,:,:,:) !nocc,nocc,nocc,nocc
    Real(realk),pointer  :: gao(:,:,:,:) !nbasis,nbasis,nbasis,nbasis
    integer :: i,j,p,q,c,m,mynum,numnodes,nAtoms,lupri
    !========================================================
    ! RI variables
    integer :: nAux,NBA
    real(realk),pointer :: CalphaR(:),CalphaG(:),CalphaF(:)
    real(realk),pointer :: CalphaRocc(:),CalphaGocc(:)
    real(realk),pointer :: Cfull(:,:),ABdecompR(:,:),ABdecompG(:,:)
    real(realk),pointer :: ABdecompF(:,:)
    logical :: Test 
    logical :: master,wakeslaves,ABdecompCreateR,ABdecompCreateG,ABdecompCreateF
    logical :: RIF12,FORCEPRINT,use_bg_buf
    character :: intspec(4)
    type(matrix) :: CMO_CABS,CMO_RI
    if(MyMolecule%mem_distributed)then
       call lsquit("ERROR(full_canonical_rimp2_f12): does not work with PDM type fullmolecule",-1)
    endif

    Test =.TRUE.
    RIF12 = .TRUE.
    lupri = DECinfo%output
#ifdef VAR_TIME
    FORCEPRINT = .TRUE.
#else
    FORCEPRINT = .FALSE.
#endif    
#ifdef VAR_MPI
    master= (infpar%mynum == infpar%master)
    mynum = infpar%mynum
    numnodes = infpar%nodtot
    wakeslaves = infpar%nodtot.GT.1
    if(.NOT.master)lupri = 6
#else
    ! If MPI is not used, consider the single node to be "master"
    master=.true.
    mynum = 0
    numnodes = 1
    wakeslaves = .false.
#endif

    ! Init stuff
    ! **********
    nbasis = MyMolecule%nbasis
    nocc   = MyMolecule%nocc
    nvirt  = MyMolecule%nvirt
    naux   = MyMolecule%nauxbasis
    nAtoms = MyMolecule%nAtoms
    call determine_CABS_nbast(ncabsAO,ncabsMO,mylsitem%setting,DECinfo%output)

    call mat_init(CMO_CABS,nCabsAO,ncabsMO)
    call build_CABS_MO(CMO_CABS,nCabsAO,mylsitem%SETTING,lupri)    
    call mat_init(CMO_RI,nCabsAO,nCabsAO)
    call build_RI_MO(CMO_RI,nCabsAO,mylsitem%SETTING,lupri)

    IF(naux.EQ.0)call lsquit('Error no Aux functions in full_canonical_rimp2_f12',-1)

    noccfull = nocc
    IF(DECinfo%frozencore)call lsquit('RI-MP2-F12 frozen core not implemented',-1)

    IF(RIF12)THEN
       call mem_alloc(Cfull,nbasis,nbasis)
       do J=1,nocc
          do I=1,nbasis
             Cfull(I,J) = MyMolecule%Co%elm2(I,J)
          enddo
       enddo
       do P=1,nvirt
          do I=1,nbasis
             Cfull(I,nocc+P) = MyMolecule%Cv%elm2(I,P)
          enddo
       enddo
    ENDIF

    !==============================================================
    !=                                                            =
    != Step 1:  Fijkl                                             =
    != The Gaussian geminal divided by the Coulomb operator g/r12 =     
    !=                                                            =
    !==============================================================
    call II_get_CoulombEcont(DECinfo%output,DECinfo%output,mylsitem%setting,&
         & [Dmat],Econt,1,GGemCouOperator)
    CoulombF12 = Econt(1) 
    Econt(1) = 0.0E0_realk
    call II_get_exchangeEcont(DECinfo%output,DECinfo%output,mylsitem%setting,&
         & [Dmat],Econt,1,GGemCouOperator)
    ExchangeF12 = Econt(1)       
    E21 = -0.5E0_realk*((5.0E0_realk/4.0E0_realk)*CoulombF12+ExchangeF12*0.5E0_realk)
    WRITE(DECINFO%OUTPUT,*)'E(Fijkl,LS) = ',E21
    mp2f12_energy = E21
    use_bg_buf = .FALSE.
    IF(RIF12)THEN       !Use RI 
       call mem_alloc(ABdecompF,nAux,nAux)
       ABdecompCreateF = .TRUE.
!       call Build_CalphaMO(mylsitem,master,nbasis,nAux,LUPRI,FORCEPRINT,&
!            & wakeslaves,MyMolecule%Co%elm2,nocc,MyMolecule%Co%elm2,nocc,mynum,numnodes,nAtoms,&
!            & CalphaF,NBA,ABdecompF,ABdecompCreateF,GGemCouOperator)
       intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
       intspec(2) = 'R' !Regular AO basis function on center 3
       intspec(3) = 'R' !Regular AO basis function on center 4
       intspec(4) = 'F' !The Gaussian geminal divided by the Coulomb operator g/r12 (GGemCouOperator)
       call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
            & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,MyMolecule%Co%elm2,nocc,&
            & mynum,numnodes,CalphaF,NBA,ABdecompF,ABdecompCreateF,intspec,use_bg_buf)
       ABdecompCreateF = .FALSE.
       call ContractOne4CenterF12IntegralsRI(NBA,nocc,CalphaF,CoulombF12,ExchangeF12)
       E21 = -0.5E0_realk*((5.0E0_realk/2.0E0_realk)*CoulombF12-ExchangeF12*0.5E0_realk)
       mp2f12_energy = mp2f12_energy  + E21
       WRITE(DECINFO%OUTPUT,*)'E(Fijkl,RI) = ',E21       
       WRITE(DECINFO%OUTPUT,*)'E(Fijkl,RI,Coulomb) = ',CoulombF12
       WRITE(DECINFO%OUTPUT,*)'E(Fijkl,RI,Exchange) = ',ExchangeF12
       call mem_dealloc(CalphaF)
       call mem_dealloc(ABdecompF)
    ENDIF


    IF(Test)THEN
       call mem_alloc(Vijij,nocc,nocc)
       call mem_alloc(Vjiij,nocc,nocc)
       call mem_alloc(gao,nbasis,nbasis,nbasis,nbasis)
       call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRF')
       call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
            &                          MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'iiii',gAO,Fijkl)
       call mem_dealloc(gao)
       do j=1,nocc
          do i=1,nocc
             Vijij(i,j) = Fijkl(i,i,j,j)
          enddo
       enddo
       do j=1,nocc
          do i=1,nocc
             Vjiij(i,j) = Fijkl(i,j,j,i)
          enddo
       enddo
       call mem_dealloc(Fijkl)
       E21 = 2.0E0_REALK*mp2f12_E21(Vijij,Vjiij,nocc)
       call mem_dealloc(Vijij)
       call mem_dealloc(Vjiij)
       WRITE(DECINFO%OUTPUT,*)'E(Fijkl,FullIntegral) = ',E21
    ENDIF
    !==========================================================
    !=                                                        =
    !=             Step 2  Ripjq*Gipjq                        =
    != The Coulomb Operator Int multiplied with               =
    != The Gaussian geminal operator g                        =
    != Dim(nocc,nbasis,nocc,nbasis)                           =
    !=                                                        =
    !==========================================================
    !Exchange Ripjq*Gjpiq
    IF(RIF12)THEN       !Use RI 
       call mem_alloc(ABdecompR,nAux,nAux)
       ABdecompCreateR = .TRUE.
       intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
       intspec(2) = 'R' !Regular AO basis function on center 3
       intspec(3) = 'R' !Regular AO basis function on center 4
       intspec(4) = 'C' !The Coulomb Operator
       call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
            & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,Cfull,nbasis,&
            & mynum,numnodes,CalphaR,NBA,ABdecompR,ABdecompCreateR,intspec,use_bg_buf)
       ABdecompCreateR = .FALSE.
       call mem_alloc(ABdecompG,nAux,nAux)
       ABdecompCreateG = .TRUE.
       intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
       intspec(2) = 'R' !Regular AO basis function on center 3
       intspec(3) = 'R' !Regular AO basis function on center 4
       intspec(4) = 'G' !The Gaussian geminal operator g
       call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
            & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,Cfull,nbasis,&
            & mynum,numnodes,CalphaG,NBA,ABdecompG,ABdecompCreateG,intspec,use_bg_buf)
       ABdecompCreateG = .FALSE.
       call ContractTwo4CenterF12IntegralsRI(nBA,nocc,nbasis,CalphaR,CalphaG,&
            & CoulombF12,ExchangeF12)
       E21 = 0.5E0_realk*((5.0E0_realk/2.0E0_realk)*CoulombF12-ExchangeF12*0.5E0_realk)
       mp2f12_energy = mp2f12_energy  + E21
       WRITE(DECINFO%OUTPUT,*)'E(Ripjq*Gjpiq,RI) = ',E21       
       WRITE(DECINFO%OUTPUT,*)'E(Ripjq*Gjpiq,RI,Coulomb) = ',CoulombF12
       WRITE(DECINFO%OUTPUT,*)'E(Ripjq*Gjpiq,RI,Exchange) = ',ExchangeF12
       call mem_dealloc(CalphaR)
       call mem_dealloc(CalphaG)
       call mem_dealloc(Cfull)
       call mem_dealloc(ABdecompR)
       call mem_dealloc(ABdecompG)
    ELSE
       !Something else that is fast????? and small memory
    ENDIF
    
    IF(Test)THEN
       call mem_alloc(Vijij,nocc,nocc)
       call mem_alloc(Vjiij,nocc,nocc)
       call mem_alloc(gao,nbasis,nbasis,nbasis,nbasis)
       call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRC')
       call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
            &                          MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'ipip',gAO,Ripjq)
       call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRG')
       call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
            &                          MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'ipip',gAO,Gipjq)
       call mem_dealloc(gao)
       Vijij = 0.0E0_realk
       Vjiij = 0.0E0_realk
       do q=1,nbasis
          do j=1,nocc
             do p=1,nbasis
                do i=1,nocc
                   Vijij(i,j) = Vijij(i,j) - Ripjq(i,p,j,q)*Gipjq(i,p,j,q)
                enddo
             enddo
          enddo
       enddo
       
       do q=1,nbasis
          do j=1,nocc
             do p=1,nbasis
                do i=1,nocc
                   Vjiij(i,j) = Vjiij(i,j) - Ripjq(i,p,j,q)*Gipjq(j,p,i,q)
                enddo
             enddo
          enddo
       enddo
       call mem_dealloc(Ripjq)
       call mem_dealloc(Gipjq)
       
       E21 = 2.0E0_REALK*mp2f12_E21(Vijij,Vjiij,nocc)
       WRITE(DECINFO%OUTPUT,*)'E(Ripjq*Gipjq,FullIntegral) = ',E21
       call mem_dealloc(Vijij)
       call mem_dealloc(Vjiij)
    ENDIF

    !==========================================================
    !=                                                        =
    !=             Step 3  Rimjc*Gimjc                        =
    != The Coulomb Operator Int multiplied with               =
    != The Gaussian geminal operator g                        =
    != Dim: (nocc,noccfull,nocc,ncabsMO)  need 4 Calphas        =
    !=                                                        =
    !==========================================================
    IF(RIF12)THEN !Use RI 
       call mem_alloc(ABdecompR,nAux,nAux)
       ABdecompCreateR = .TRUE.
       !CalphaRocc(NBA,nocc,nocc)
       intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
       intspec(2) = 'R' !Regular AO basis function on center 3
       intspec(3) = 'R' !CABS AO basis function on center 4
       intspec(4) = 'C' !The Coulomb Operator
       call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
            & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,MyMolecule%Co%elm2,nocc,&
            & mynum,numnodes,CalphaRocc,NBA,ABdecompR,ABdecompCreateR,intspec,use_bg_buf)
       !CalphaR(NBA,nocc,ncabsMO)
       ABdecompCreateR = .FALSE.
       intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
       intspec(2) = 'R' !Regular AO basis function on center 3
       intspec(3) = 'C' !CABS AO basis function on center 4
       intspec(4) = 'C' !The Coulomb Operator
       call Build_CalphaMO2(mylsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
            & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,CMO_CABS%elms,ncabsMO,&
            & mynum,numnodes,CalphaR,NBA,ABdecompR,ABdecompCreateR,intspec,use_bg_buf)
       call mem_dealloc(ABdecompR)
       !CalphaGocc(NBA,nocc,nocc)
       call mem_alloc(ABdecompG,nAux,nAux)
       ABdecompCreateG = .TRUE.
       intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
       intspec(2) = 'R' !Regular AO basis function on center 3
       intspec(3) = 'R' !Regular AO basis function on center 4
       intspec(4) = 'G' !The Gaussian geminal operator g
       call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
            & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,MyMolecule%Co%elm2,nocc,&
            & mynum,numnodes,CalphaGocc,NBA,ABdecompG,ABdecompCreateG,intspec,use_bg_buf)
       !CalphaG(NBA,nocc,ncabsMO)
       ABdecompCreateG = .FALSE.
       intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
       intspec(2) = 'R' !Regular AO basis function on center 3
       intspec(3) = 'C' !CABS AO basis function on center 4
       intspec(4) = 'G' !The Gaussian geminal operator g
       call Build_CalphaMO2(mylsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
            & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,CMO_CABS%elms,ncabsMO,&
            & mynum,numnodes,CalphaG,NBA,ABdecompG,ABdecompCreateG,intspec,use_bg_buf)
       call mem_dealloc(ABdecompG)

       call ContractTwo4CenterF12IntegralsRI2(NBA,nocc,noccfull,ncabsMO,CalphaR,CalphaG,&
            & CalphaRocc,CalphaGocc,CoulombF12,ExchangeF12)
       E21 = 0.5E0_realk*((5.0E0_realk/2.0E0_realk)*CoulombF12-ExchangeF12*0.5E0_realk)
       mp2f12_energy = mp2f12_energy  + E21
       WRITE(DECINFO%OUTPUT,*)'E(Rimjc*Gjmic,RI) = ',E21       
       WRITE(DECINFO%OUTPUT,*)'E(Rimjc*Gjmic,RI,Coulomb) = ',CoulombF12
       WRITE(DECINFO%OUTPUT,*)'E(Rimjc*Gjmic,RI,Exchange) = ',ExchangeF12
       call mem_dealloc(CalphaRocc)
       call mem_dealloc(CalphaGocc)
       call mem_dealloc(CalphaR)
       call mem_dealloc(CalphaG)
    ELSE
       !Something else that is fast????? and small memory
    ENDIF

    IF(Test)THEN
       call mem_alloc(Vijij,nocc,nocc)
       call mem_alloc(Vjiij,nocc,nocc)
       call mem_alloc(gao,nbasis,nbasis,nbasis,ncabsAO)
       call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRCC')
       call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
            & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'imic',gAO,Rimjc)
       call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRCG')
       call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
            & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'imic',gAO,Gimjc)
       call mem_dealloc(gao)

       call mem_alloc(Vijij,nocc,nocc)
       call mem_alloc(Vjiij,nocc,nocc)
       Vijij = 0.0E0_realk
       Vjiij = 0.0E0_realk
       do c=1,ncabsMO
          do j=1,nocc
             do m=1,noccfull
                do i=1,nocc
                   Vijij(i,j) = Vijij(i,j) - Rimjc(i,m,j,c)*Gimjc(i,m,j,c)
                enddo
             enddo
          enddo
       enddo
       do c=1,ncabsMO
          do j=1,nocc
             do m=1,noccfull
                do i=1,nocc
                   Vjiij(i,j) = Vjiij(i,j) - Rimjc(i,m,j,c)*Gimjc(j,m,i,c)
                enddo
             enddo
          enddo
       enddo
       
       E21 = 2.0E0_REALK*mp2f12_E21(Vijij,Vjiij,nocc)
       WRITE(DECINFO%OUTPUT,*)'E(Rimjc*Gimjc,FullIntegral) = ',E21


!!$       print*,'Rimjc(1,1,1,1)=',Rimjc(1,1,1,1)
!!$       print*,'Rimjc(nocc,1,1,1)=',Rimjc(nocc,1,1,1)
!!$       print*,'Rimjc(1,noccfull,1,1)=',Rimjc(1,noccfull,1,1)
!!$       print*,'Rimjc(1,1,nocc,1)=',Rimjc(1,1,nocc,1)
!!$       print*,'Rimjc(1,1,1,ncabsMO)=',Rimjc(1,1,1,ncabsMO)
!!$       print*,'Rimjc(nocc,1,nocc,1)=',Rimjc(nocc,1,nocc,1)
!!$       print*,'Gimjc(1,1,1,1)=',Gimjc(1,1,1,1)
!!$       print*,'Gimjc(nocc,1,1,1)=',Gimjc(nocc,1,1,1)
!!$       print*,'Gimjc(1,noccfull,1,1)=',Gimjc(1,noccfull,1,1)
!!$       print*,'Gimjc(1,1,nocc,1)=',Gimjc(1,1,nocc,1)
!!$       print*,'Gimjc(1,1,1,ncabsMO)=',Gimjc(1,1,1,ncabsMO)
!!$       print*,'Gimjc(nocc,1,nocc,1)=',Gimjc(nocc,1,nocc,1)

       call mem_dealloc(Rimjc)
       call mem_dealloc(Gimjc)

       call mem_dealloc(Vijij)
       call mem_dealloc(Vjiij)
    ENDIF

  end subroutine full_canonical_rimp2_f12

  !> Function for finding the E21 energy  
  function mp2f12_E21(Vijij,Vjiij,nocc) result(energy)
    implicit none
    Integer,intent(IN)     :: nocc
    Real(realk),intent(IN) :: Vijij(nocc,nocc),Vjiij(nocc,nocc)
    Real(realk) :: energy
    !
    Integer     :: i,j
    Real(realk) :: tmp
    
    tmp = 0E0_realk
    DO i=1,nocc
       tmp = tmp + Vijij(i,i)
    ENDDO
    
    energy = -0.5E0_realk*tmp
    tmp = 0E0_realk
    
    DO j=1,nocc
       DO i=j+1,nocc
          tmp = tmp + 5E0_realk * Vijij(i,j) - Vjiij(i,j)
       ENDDO
    ENDDO
    energy = energy - 0.25E0_realk*tmp
  end function mp2f12_E21


  function mp2f12_E21A(Vijij,Vjiij,nocc) result(energy)
    implicit none
    Integer,intent(IN)     :: nocc
    Real(realk),intent(IN) :: Vijij(nocc,nocc),Vjiij(nocc,nocc)
    Real(realk) :: energy
    !
    Integer     :: i,j
    Real(realk) :: tmp
    
    tmp = 0E0_realk
    DO i=1,nocc
       tmp = tmp + Vijij(i,i)
    ENDDO    
    energy = -0.5E0_realk*tmp
  end function mp2f12_E21A

  function mp2f12_E21B(Vijij,Vjiij,nocc) result(energy)
    implicit none
    Integer,intent(IN)     :: nocc
    Real(realk),intent(IN) :: Vijij(nocc,nocc),Vjiij(nocc,nocc)
    Real(realk) :: energy
    !
    Integer     :: i,j
    Real(realk) :: tmp
    
    tmp = 0E0_realk    
    DO j=1,nocc
       DO i=j+1,nocc
          tmp = tmp + 5E0_realk * Vijij(i,j)
       ENDDO
    ENDDO
    energy = - 0.25E0_realk*tmp
  end function mp2f12_E21B

  function mp2f12_E21B2(Vijij,Vjiij,nocc) result(energy)
    implicit none
    Integer,intent(IN)     :: nocc
    Real(realk),intent(IN) :: Vijij(nocc,nocc),Vjiij(nocc,nocc)
    Real(realk) :: energy
    !
    Integer     :: i,j
    Real(realk) :: tmp
    
    tmp = 0E0_realk    
    DO j=1,nocc
       DO i=1,nocc
          tmp = tmp + 5E0_realk * Vijij(i,j)
       ENDDO
    ENDDO
    energy = - 0.25E0_realk*tmp
  end function mp2f12_E21B2

  function mp2f12_E21C(Vijij,Vjiij,nocc) result(energy)
    implicit none
    Integer,intent(IN)     :: nocc
    Real(realk),intent(IN) :: Vijij(nocc,nocc),Vjiij(nocc,nocc)
    Real(realk) :: energy
    !
    Integer     :: i,j
    Real(realk) :: tmp
    
    tmp = 0E0_realk    
    DO j=1,nocc
       DO i=j+1,nocc
          tmp = tmp - Vjiij(i,j)
       ENDDO
    ENDDO
    energy = - 0.25E0_realk*tmp
  end function mp2f12_E21C

  function mp2f12_E21C2(Vijij,Vjiij,nocc) result(energy)
    implicit none
    Integer,intent(IN)     :: nocc
    Real(realk),intent(IN) :: Vijij(nocc,nocc),Vjiij(nocc,nocc)
    Real(realk) :: energy
    !
    Integer     :: i,j
    Real(realk) :: tmp
    
    tmp = 0E0_realk    
    DO j=1,nocc
       DO i=1,nocc
          tmp = tmp - Vjiij(i,j)
       ENDDO
    ENDDO
    energy = - 0.25E0_realk*tmp
  end function mp2f12_E21C2
  
!Perform on GPU 
subroutine ContractOne4CenterF12IntegralsRI(nBA,n,Calpha,EJ,EK)
  implicit none
  integer,intent(in)        :: nBA,n
  real(realk),intent(in)    :: Calpha(nBA,n,n)
  real(realk),intent(inout) :: EJ,EK
  !local variables
  integer :: I,ALPHA,J
  real(realk) :: TMP,TMPV(n),TMPI
  !Exchange Fiijj
  EK = 0.0E0_realk
  EJ = 0.0E0_realk
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(I,J,&
  !$OMP ALPHA) SHARED(Calpha,n,nba) REDUCTION(+:EK,EJ)
  DO I=1,n
     DO J=1,n
        DO ALPHA = 1,NBA
           EJ = EJ + Calpha(ALPHA,I,I)*Calpha(ALPHA,J,J)
           EK = EK + Calpha(ALPHA,I,J)*Calpha(ALPHA,J,I)
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO
end subroutine ContractOne4CenterF12IntegralsRI

subroutine ContractTwo4CenterF12IntegralsRI(nBA,n1,n2,CalphaR,CalphaG,EJ,EK)
  implicit none
  integer,intent(in)        :: nBA,n1,n2
  real(realk),intent(in)    :: CalphaR(nBA,n1,n2),CalphaG(nBA,n1,n2)
  real(realk),intent(inout) :: EJ,EK
  !local variables
  integer :: Q,P,I,J,ALPHA,BETA
  real(realk) :: TMPR,TMPG,TMP1(NBA,NBA),TMP2(NBA,NBA)
  !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
  EK = 0.0E0_realk
  EJ = 0.0E0_realk
  !$OMP PARALLEL DEFAULT(none) PRIVATE(I,J,P,Q,TMPR,TMPG) SHARED(CalphaR,CalphaG,n2,n1,nba) REDUCTION(+:EK,EJ)
  !$OMP DO COLLAPSE(3)
  DO Q=1,n2
     DO P=1,n2
        DO I=1,n1
           DO J=1,n1
              TMPR = 0.0E0_realk
              DO ALPHA = 1,NBA
                 TMPR = TMPR + CalphaR(ALPHA,I,P)*CalphaR(ALPHA,J,Q)
              ENDDO
              TMPG = 0.0E0_realk
              DO BETA = 1,NBA
                 TMPG = TMPG + CalphaG(BETA,J,P)*CalphaG(BETA,I,Q)
              ENDDO
              EK = EK + TMPR*TMPG
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$OMP END DO NOWAIT
  !Coulomb Ripjq*Gipjq Scaling(N*O*Naux**2)
  !FIXME THIS IS NOT OPTIMAL - NON UNIT STRIDE 
  !$OMP DO COLLAPSE(2)
  DO BETA = 1,NBA
     DO ALPHA = 1,NBA
        TMPR = 0.0E0_realk                    
        DO P=1,n2
           DO I=1,n1
              TMPR = TMPR + CalphaR(ALPHA,I,P)*CalphaG(BETA,I,P)
           ENDDO
        ENDDO
        EJ = EJ + TMPR*TMPR
     ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL


end subroutine ContractTwo4CenterF12IntegralsRI

subroutine ContractTwo4CenterF12IntegralsRI2(nBA,n1,n3,n2,CalphaR,CalphaG,&
     & CalphaRocc,CalphaGocc,EJ,EK)
  implicit none
  integer,intent(in)        :: nBA,n1,n2,n3
  real(realk),intent(in)    :: CalphaR(nBA,n1,n2),CalphaG(nBA,n1,n2)
  real(realk),intent(in)    :: CalphaRocc(nBA,n1,n3),CalphaGocc(nBA,n1,n3)
  real(realk),intent(inout) :: EJ,EK
  !local variables
  integer :: M,C,I,J,ALPHA,BETA
  real(realk) :: TMPR,TMPG
  !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
  EK = 0.0E0_realk
  EJ = 0.0E0_realk
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(I,J,M,C,TMPR,&
  !$OMP TMPG) SHARED(CalphaR,CalphaRocc,CalphaG,CalphaGocc,n3,n2,n1,&
  !$OMP nba) REDUCTION(+:EK,EJ)
  DO M=1,n3
     DO C=1,n2
        DO I=1,n1
           DO J=1,n1
              TMPR = 0.0E0_realk
              DO ALPHA = 1,NBA
                 TMPR = TMPR + CalphaRocc(ALPHA,I,M)*CalphaR(ALPHA,J,C)
              ENDDO
              TMPG = 0.0E0_realk
              DO BETA = 1,NBA
                 TMPG = TMPG + CalphaGocc(BETA,J,M)*CalphaG(BETA,I,C)
              ENDDO
              EK = EK + TMPR*TMPG
!              TMPR = 0.0E0_realk
!              DO ALPHA = 1,NBA
!                 TMPR = TMPR + CalphaRocc(ALPHA,I,M)*CalphaR(ALPHA,J,C)
!              ENDDO
              TMPG = 0.0E0_realk
              DO BETA = 1,NBA
                 TMPG = TMPG + CalphaGocc(BETA,I,M)*CalphaG(BETA,J,C)
              ENDDO
              EJ = EJ + TMPR*TMPG
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO

!!$  print*,'n1',n1,'n2',n2
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaRocc(ALPHA,1,1)*CalphaR(ALPHA,1,1)
!!$  ENDDO
!!$  print*,'Rimjc(1,1,1,1)=',TMPR
!!$  
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaRocc(ALPHA,n1,1)*CalphaR(ALPHA,1,1)
!!$  ENDDO
!!$  print*,'Rimjc(n1,1,1,1)=',TMPR
!!$
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaRocc(ALPHA,1,n3)*CalphaR(ALPHA,1,1)
!!$  ENDDO
!!$  print*,'Rimjc(1,n3,1,1)=',TMPR
!!$  
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaRocc(ALPHA,1,1)*CalphaR(ALPHA,n1,1)
!!$  ENDDO
!!$  print*,'Rimjc(1,1,n1,1)=',TMPR
!!$
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaRocc(ALPHA,1,1)*CalphaR(ALPHA,1,n2)
!!$  ENDDO
!!$  print*,'Rimjc(1,1,1,n2)=',TMPR
!!$
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaRocc(ALPHA,n1,1)*CalphaR(ALPHA,n1,1)
!!$  ENDDO
!!$  print*,'Rimjc(n1,1,n1,1)=',TMPR
!!$
!!$
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaGocc(ALPHA,1,1)*CalphaG(ALPHA,1,1)
!!$  ENDDO
!!$  print*,'Gimjc(1,1,1,1)=',TMPR
!!$  
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaGocc(ALPHA,n1,1)*CalphaG(ALPHA,1,1)
!!$  ENDDO
!!$  print*,'Gimjc(n1,1,1,1)=',TMPR
!!$
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaGocc(ALPHA,1,n3)*CalphaG(ALPHA,1,1)
!!$  ENDDO
!!$  print*,'Gimjc(1,n3,1,1)=',TMPR
!!$  
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaGocc(ALPHA,1,1)*CalphaG(ALPHA,n1,1)
!!$  ENDDO
!!$  print*,'Gimjc(1,1,n1,1)=',TMPR
!!$
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaGocc(ALPHA,1,1)*CalphaG(ALPHA,1,n2)
!!$  ENDDO
!!$  print*,'Gimjc(1,1,1,n2)=',TMPR
!!$
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaGocc(ALPHA,n1,1)*CalphaG(ALPHA,n1,1)
!!$  ENDDO
!!$  print*,'Gimjc(n1,1,n1,1)=',TMPR

end subroutine ContractTwo4CenterF12IntegralsRI2

  !Exchange Ripjq*Gjpiq
  !Ripjq*Gjpiq = (AB|OperR|CD)*(EF|OperR|GH)*Docc(A,G)*Docc(C,E)*Dvirt(B,F)*Dvirt(D,H)
!!$  subroutine ContractTwo4CenterF12IntegralsExchange(OperR,OperG,nocc,nbasis,SETTING,INTSPEC)
!!$    implicit none
!!$    integer,intent(in) :: OperR,OperG,nocc,nbasis
!!$    TYPE(LSSETTING),intent(inout)     :: SETTING
!!$    
!!$    NOFAMILY = ls%setting%SCHEME%NOFAMILY
!!$    ls%setting%SCHEME%NOFAMILY = .TRUE.
!!$    !TODO: 
!!$    nullify(batchsize)
!!$    nullify(batchdim)
!!$    nullify(batchindex)
!!$    nullify(orb2batch)
!!$    nullify(batch2orb)
!!$
!!$    doscreen = ls%setting%SCHEME%CS_SCREEN.OR.ls%setting%SCHEME%PS_SCREEN
!!$
!!$    call build_minimalbatchesofAOS(DECinfo%output,setting,nbasis,&
!!$         & batchsize,batchdim,batchindex,nbatches,orb2Batch,INTSPEC(1))
!!$
!!$    call mem_alloc(batch2orb,nbatchesAB)
!!$    do idx=1,nbatchesAB
!!$       call mem_alloc(batch2orb(idx)%orbindex,batchdim(idx) )
!!$       batch2orb(idx)%orbindex = 0
!!$       batch2orb(idx)%norbindex = 0
!!$    end do
!!$    do iorb=1,nbast
!!$       idx = orb2batch(iorb)
!!$       batch2orb(idx)%norbindex = batch2orb(idx)%norbindex+1
!!$       k = batch2orb(idx)%norbindex
!!$       batch2orb(idx)%orbindex(k) = iorb
!!$    end do
!!$
!!$    call mem_alloc(Dbast,nbasis,nbasis)
!!$    call DGEMM
!!$    call mem_alloc(Docc,nbatches,nbatches)
!!$    call ConvertBASTGabToBatchesGab(nbasis,nbatches,setting,Dbast,Docc,lupri,luerr)
!!$    call mem_dealloc(Dbast)
!!$    MaxDocc = MAXVAL(Docc)
!!$    call mem_alloc(MaxDoccV,nbasis)
!!$    DO J=1,nbasis
!!$       MaxDoccV(J) = MAXVAL(MaxDocc(:,J))
!!$    ENDDO
!!$
!!$    call mem_alloc(Dbast,nbasis,nbasis)
!!$    call DGEMM
!!$    call mem_alloc(Dvirt,nbatches,nbatches)
!!$    call ConvertBASTGabToBatchesGab(nbasis,nbatches,setting,Dbast,Dvirt,lupri,luerr)
!!$    call mem_dealloc(Dbast)
!!$    MaxDvirt = MAXVAL(Dvirt)
!!$    call mem_alloc(MaxDvirtV,nbasis)
!!$    DO J=1,nbasis
!!$       MaxDvirtV(J) = MAXVAL(MaxDvirt(:,J))
!!$    ENDDO
!!$
!!$    call mem_alloc(Rscreen,nbatches,nbatches)
!!$    call II_get_2int_BatchScreenMat(DECinfo%output,DECinfo%output,SETTING,&
!!$         & nbatches,Rscreen,nbasis,OperR)
!!$
!!$    WRITE(lupri,*)'Rscreen'
!!$    call ls_output(Rscreen,1,nbatches,1,nbatches,nbatches,nbatches,1,lupri)
!!$
!!$    MaxRscreen = MAXVAL(Rscreen)
!!$
!!$    call mem_alloc(Gscreen,nbatches,nbatches)
!!$    call II_get_2int_BatchScreenMat(DECinfo%output,DECinfo%output,SETTING,&
!!$         & nbatches,Gscreen,nbasis,OperG)
!!$
!!$    WRITE(lupri,*)'Gscreen'
!!$    call ls_output(Gscreen,1,nbatches,1,nbatches,nbatches,nbatches,1,lupri)
!!$
!!$    MaxGscreen = MAXVAL(Gscreen)
!!$
!!$    Threshold_CS = ls%SETTING%SCHEME%THRESHOLD*ls%SETTING%SCHEME%J_THR
!!$    intThreshold = ls%SETTING%SCHEME%THRESHOLD*ls%SETTING%SCHEME%J_THR
!!$    call II_precalc_DECScreenMat(DecScreenR,lupri,luerr,ls%setting,nbatches,&
!!$         & nbatches,INTSPEC,intThreshold)
!!$    IF(doscreen)then
!!$       call II_getBatchOrbitalScreen(DecScreenR,ls%setting,&
!!$            & nbasis,nbatches,nbatches,batchsize,batchsize,batchindex,batchindex,&
!!$            & batchdim,batchdim,INTSPEC,lupri,luerr)
!!$    endif
!!$
!!$    call II_precalc_DECScreenMat(DecScreenG,lupri,luerr,ls%setting,nbatches,&
!!$         & nbatches,INTSPEC,intThreshold)
!!$    IF(doscreen)then
!!$       call II_getBatchOrbitalScreen(DecScreenG,ls%setting,&
!!$            & nbasis,nbatches,nbatches,batchsize,batchsize,batchindex,batchindex,&
!!$            & batchdim,batchdim,INTSPEC,lupri,luerr)
!!$    endif
!!$    FullRHS = .FALSE.
!!$
!!$    !Ripjq*Gjpiq = (AB|OperR|CD)*(EF|OperR|GH)*Docc(A,G)*Docc(C,E)*Dvirt(B,F)*Dvirt(D,H)
!!$    !starting loop over the most sparse of the 2 operators (create list )
!!$    !MPI PARALLIZE THE D,C Loop
!!$    !E = 0
!!$    !do C
!!$    ! do A
!!$    !   IF(RscreenV(A)*Rscreen(C))THEN
!!$    !    done Nac < nbatchA*nbatchC times 
!!$    !    Calc (A,Bfull|OperR|C,Dfull)
!!$    !    construct Rtensor(A,Bfull,C,Dfull)
!!$    !    construct Rtensor(A,Bfull,C,Hfull)=Dvirt(Dfull,Hfull)*Rtensor(A,Bfull,C,Dfull) !N*N*N*Nac
!!$    !    construct Rtensor(A,Ffull,C,Hfull)=Dvirt(Ffull,Bfull)*Rtensor(A,Bfull,C,Hfull)
!!$    !    do G
!!$    !      construct Rtensor(G,Ffull,C,Hfull)=Docc(G,A)*Rtensor(A,Ffull,C,Hfull)
!!$    !      do E On GPU - While CPU does Integral GPU does DGEMM? 
!!$    !        construct Rtensor(G,Ffull,E,Hfull)=Docc(C,E)*Rtensor(G,Ffull,C,Hfull)
!!$    !        maxR = MAXVAL(Rtensor(G,Ffull,E,Hfull)) 
!!$    !        IF(maxR*GscreenV(E)*GscreenV(G))THEN
!!$    !          Calc (E,Ffull|OperG|G,Hfull) !Modified Screening Threshold with maxR
!!$    !          E = E + (E,Ffull|OperG|G,Hfull)*Rtensor(G,Ffull,E,Hfull)
!!$    !        ENDIF
!!$    !      enddo
!!$    !    enddo
!!$    !   ENDIF
!!$    ! enddo
!!$    !enddo
!!$
!!$
!!$    !TODO VERIFY THAT THE SCREENING IS CORRECT - when you do MaxValRjFiD*MaxGscreen*Gscreen(G,H)*Dvirt(D,H)*MaxCMOV(G)*MaxCMO*maxDvirt
!!$    ! it is not correct as you are summing over the elements - but can do something like what they do in Lapalce 
!!$    !MaxValRjFiD*MaxValGjFiH
!!$    !where MaxValGjFiH = MaxGscreen*SUM_G (MaxCMOV(G)*Gscreen(G,H))
!!$    !Psedo Code 
!!$    !Ripjq*Gjpiq = (AB|OperR|CD)*(EF|OperR|GH)*Docc(A,G)*Docc(C,E)*Dvirt(B,F)*Dvirt(D,H)
!!$    !starting loop over the most sparse of the 2 operators (create list )
!!$    !MPI PARALLIZE THE D,C Loop
!!$    !E = 0
!!$    !do D
!!$    ! do C
!!$    !  IF()
!!$    !   done Ncd < nbatchC*nbatchD times 
!!$    !   (std integral AO to MO speed) expect taking CD screening into account 
!!$    !   construct (AfullBfull|OperR|CD)          
!!$    !   construct Rtensor(Afull,Bfull,C,D)       
!!$    !   construct Rtensor(I,Bfull,C,D)=CMO(I,Afull)*Rtensor(Afull,Bfull,C,D) !DGEMM
!!$    !   Reorder   Rtensor(C,D,I,Bfull)<=Rtensor(I,Bfull,C,D)            
!!$    !   construct Rtensor(J,D,I,Bfull)=Rtensor(J,D,I,Bfull) + CMO(J,C)*Rtensor(C,D,I,Bfull) !DGEMM
!!$    !  ENDIF
!!$    ! enddo
!!$    ! Reorder   Rtensor(J,Bfull,I,D)<=Rtensor(J,D,I,Bfull)
!!$    ! construct Rtensor(J,Ffull,I,D)=Rtensor(J,Bfull,I,D)*Dvirt(Bfull,Ffull)
!!$    ! MaxValRjFiD = MAXVAL(Rtensor(J,Ffull,I,D)) !This is not a small number
!!$    ! do H 
!!$    !  do G
!!$    !   IF(MaxValRjFiD*MaxGscreen*Gscreen(G,H).GT.Threshold_CS)THEN 
!!$    !    done Ngh < nbatchG*nbatchH*nbatchD times 
!!$    !    construct Gtensor(EfullFfull|operG|GH) !USING A MODIFIED SCREENING THRESHOLD *MaxValRjFiD 
!!$    !    MAXVAL = Gtensor(Efull,Ffull,G,H)
!!$    !    IF(MAXVAL*MaxValRjFiD.GT.Threshold_CS)THEN 
!!$    !      construct Gtensor(Efull,Ffull,G,H)
!!$    !      construct Gtensor(J,Ffull,G,H)=CMO(J,Efull)*Gtensor(Efull,Ffull,G,H)                !DGEMM
!!$    !      construct Gtensor(J,Ffull,I,H)=Gtensor(J,Ffull,I,H) + CMO(I,G)*Gtensor(J,Ffull,G,H)
!!$    !    ENDIF
!!$    !   ENDIF
!!$    !  enddo
!!$    !  construct Gtensor(J,Ffull,I,D)=Gtensor(J,Ffull,I,H)*Dvirt(D,H) !DGEMM
!!$    ! enddo
!!$    ! E = E + Rtensor(J,Ffull,I,D)*Gtensor(J,Ffull,I,D)
!!$    !enddo
!!$
!!$    BatchD: do D = 1,nbatches
!!$       dimD = batchdim(D)
!!$       
!!$       BatchC: do C = 1,nbatches
!!$          dimC = batchdim(C)           
!!$          
!!$          !Ripjq*Gjpiq = (AB|OperR|CD)*(EF|OperR|GH)*Docc(A,G)*Docc(C,E)*Dvirt(B,F)*Dvirt(D,H)
!!$          IF(MaxRscreen*Rscreen(C,D)*MaxGscreen*MaxGscreen*MaxDvirt(D)*MaxDocc*MaxDoccV(C)*maxDvirt.GT.Threshold_CS)THEN
!!$
!!$             IF(doscreen)ls%setting%LST_GAB_LHS => DECSCREENR%masterGabLHS
!!$             IF(doscreen)ls%setting%LST_GAB_RHS => DECSCREENR%batchGab(C,D)%p
!!$          
!!$             call II_GET_DECPACKED4CENTER_J_ERI(LUPRI,LUERR,ls%SETTING,&
!!$                  & integralsR,batchindex(C),batchindex(D),batchsize(C),batchsize(D),&
!!$                  & nbast,nbast,dimC,dimD,fullRHS,INTSPEC,intThreshold)
!!$
!!$             !          do batch_iD = 1,dimD
!!$             !           iD = batch2orb(D)%orbindex(batch_iD) !Global index
!!$             !           do batch_iC = 1,dimC
!!$             !            iC = batch2orb(C)%orbindex(batch_iC) !Global index
!!$             !Output integrals(1:nbasis,1:nbasis,batch_iC,batch_iD)
!!$             !Reorder to Radcb <= Rabcd
!!$             MaxValRabcd = MAXVAL(integralsR)
!!$             BatchH: do H = 1,nbatches
!!$                dimH = batchdim(H)
!!$                
!!$                BatchG: do G = 1,nbatches
!!$                   dimG = batchdim(G)           
!!$                   
!!$                   !Ripjq*Gjpiq = (AB|OperR|CD)*(EF|OperR|GH)*Docc(A,G)*Docc(C,E)*Dvirt(B,F)*Dvirt(D,H)
!!$                   IF(MaxValRabcd*MaxGscreen*Gscreen(G,H)*Dvirt(D,H)*MaxDoccV(G)*MaxDoccV(C)*maxDvirt.GT.Threshold_CS)THEN
!!$                      
!!$                      IF(doscreen)ls%setting%LST_GAB_LHS => DECSCREENG%masterGabLHS
!!$                      IF(doscreen)ls%setting%LST_GAB_RHS => DECSCREENG%batchGab(G,H)%p
!!$                      
!!$                      call II_GET_DECPACKED4CENTER_J_ERI(LUPRI,LUERR,ls%SETTING,&
!!$                           & integralsG,batchindex(G),batchindex(H),batchsize(G),batchsize(H),&
!!$                           & nbast,nbast,dimG,dimH,fullRHS,INTSPEC,intThreshold)
!!$                      
!!$                      !          do batch_iH = 1,dimH
!!$                      !           iH = batch2orb(H)%orbindex(batch_iH) !Global index
!!$                      !           do batch_iG = 1,dimG
!!$                      !            iG = batch2orb(G)%orbindex(batch_iG) !Global index
!!$                      !Output integrals(1:nbasis,1:nbasis,batch_iG,batch_iH)
!!$                      !Transform to Gabcd = (EF|OperR|GH)*Docc(A,G)*Docc(C,E)*Dvirt(B,F)*Dvirt(D,H)
!!$
!!$                      !Transform1 Gefgd = (EF|OperR|GH)*Dvirt(H,D)
!!$                      !Transform2 Gcfgd = Docc(C,E)*Gefgd
!!$                      !Reorder Ggdcf <= Gcfgd
!!$                      !Transform1 Ggdcb = Ggdcf*Dvirt(B,F)
!!$                      !Transform2 Gadcb = Docc(A,G)*Ggdcb
!!$                      ! E = Radcb*Gadcb
!!$                      
!!$                   ENDIF
!!$                enddo BatchG
!!$             enddo BatchH
!!$          ENDIF
!!$       enddo BatchC
!!$    enddo BatchD
!!$    
!!$
!!$    !Ripjq*Gjpiq = Rpiqj*Gpjqi
!!$    !starting loop over the most sparse of the 2 operators (create list )
!!$    !do A
!!$    ! do B
!!$    !  do C
!!$    !   do D
!!$    !    construct (AB|OperR|CD)    dim: (nbastOnA,nbastOnB,nbastOnC,nbastOnD)
!!$    !    do E 
!!$    !     do F
!!$    !      do G 
!!$    !       do H
!!$    !        construct (EF|operG|GH) dim: (nbastOnE,nbastOnF,nbastOnG,nbastOnH)  
!!$    !        transform to (AD|operG|CB)  dim: (nbastOnA,nbastOnD,nbastOnC,nbastOnB)
!!$    !        E = (AB|operR|CD)*(AD|operG|CB)   scaling: O(N) 
!!$    !       enddo
!!$    !      enddo
!!$    !     enddo
!!$    !    enddo
!!$    !   enddo
!!$    !  enddo
!!$    ! enddo
!!$    !enddo
!!$    ls%setting%SCHEME%NOFAMILY = NOFAMILY
!!$
!!$  end subroutine ContractTwo4CenterF12IntegralsExchange

!!$  print*,'n1',n1,'n2',n2
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaR(ALPHA,1,1)*CalphaR(ALPHA,1,1)
!!$  ENDDO
!!$  print*,'Ripjq(1,1,1,1)=',TMPR
!!$  
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaR(ALPHA,n1,1)*CalphaR(ALPHA,1,1)
!!$  ENDDO
!!$  print*,'Ripjq(n1,1,1,1)=',TMPR
!!$
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaR(ALPHA,1,n2)*CalphaR(ALPHA,1,1)
!!$  ENDDO
!!$  print*,'Ripjq(1,n2,1,1)=',TMPR
!!$  
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaR(ALPHA,1,1)*CalphaR(ALPHA,n1,1)
!!$  ENDDO
!!$  print*,'Ripjq(1,1,n1,1)=',TMPR
!!$
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaR(ALPHA,1,1)*CalphaR(ALPHA,1,n2)
!!$  ENDDO
!!$  print*,'Ripjq(1,1,1,n2)=',TMPR
!!$
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaR(ALPHA,n1,1)*CalphaR(ALPHA,n1,1)
!!$  ENDDO
!!$  print*,'Ripjq(n1,1,n1,1)=',TMPR
!!$
!!$
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaG(ALPHA,1,1)*CalphaG(ALPHA,1,1)
!!$  ENDDO
!!$  print*,'Gipjq(1,1,1,1)=',TMPR
!!$  
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaG(ALPHA,n1,1)*CalphaG(ALPHA,1,1)
!!$  ENDDO
!!$  print*,'Gipjq(n1,1,1,1)=',TMPR
!!$
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaG(ALPHA,1,n2)*CalphaG(ALPHA,1,1)
!!$  ENDDO
!!$  print*,'Gipjq(1,n2,1,1)=',TMPR
!!$  
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaG(ALPHA,1,1)*CalphaG(ALPHA,n1,1)
!!$  ENDDO
!!$  print*,'Gipjq(1,1,n1,1)=',TMPR
!!$
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaG(ALPHA,1,1)*CalphaG(ALPHA,1,n2)
!!$  ENDDO
!!$  print*,'Gipjq(1,1,1,n2)=',TMPR
!!$
!!$  TMPR = 0.0E0_realk
!!$  DO ALPHA = 1,NBA                    
!!$     TMPR = TMPR +  CalphaG(ALPHA,n1,1)*CalphaG(ALPHA,n1,1)
!!$  ENDDO
!!$  print*,'Gipjq(n1,1,n1,1)=',TMPR



end module fullrimp2f12

