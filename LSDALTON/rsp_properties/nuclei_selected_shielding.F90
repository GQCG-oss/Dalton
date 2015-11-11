!> @file
!> Contains nuclei-selected shielding subroutines
!> \author Thomas Kjaergaard, Chandan Kumar
!> \date 2014
module nuclei_selected_shielding_mod
  use precision
  use configurationType, only: configitem
  use TYPEDEFTYPE,   only: LSSETTING,lsitem
  use memory_handling, only: mem_alloc, mem_dealloc
  use files, only: lsopen, lsclose
  use lstiming
  use matrix_module, only: matrix
  use matrix_operations
  use matrix_operations_aux 
  use matrix_util
  use rspsolver
  use decompMod, only: decompItem
  use IntegralInterfaceMOD
  use II_XC_interfaceModule
  use dal_interface
  use rsp_util
  use response_wrapper_type_module, only: ALPHAinputItem
  use ls_Integral_Interface

  !***********************************************************************
  ! Driver routine for linear response (polarizability), NMR shielding and 
  ! excitation energies (quadratic response and hessian in the future)
  !
  ! This driver is to be used as a test code. 
  ! It serves as a way to test for instance new density fitting integrals 
  ! in the context of response. 
  ! It may serve as an easy way to test the response solver for several 
  ! frequencies at a time or several RHS at a time.  
  ! This driver therefor have limited functionality. 
  ! OpenRSP is used as default and have to be used to treat exotic properties. 
  !  
  ! If a new development show alot of promise you can consider making
  ! a new driver using this development throughout or implement it in
  ! OpenRSP. LSDALTON is a platform to develop new methods and as 
  ! such this linear response code serve as such a platform. 
  ! So keep it simpel. 
  !
  ! Be aware that the Code uses the same input types and input reading 
  ! as the OpenRSP drivers and framework. 
  !
  ! Written by Thomas Kjaergard 2014 based on work of Sonia Coriani 
  ! and Stinne Hoest 
  !***********************************************************************

  public :: NMRshieldresponse_RSTNS 

  private

Contains   
!#######################################################
! Compute Nuclei Shiedling tensors. 
! Debugging purpose 
!
!#######################################################

subroutine NMRshieldresponse_RSTNS(ls,molcfg,F,D,S)
  implicit none
  TYPE(lsitem),target     :: ls
  type(rsp_molcfg), intent(inout) :: molcfg
  type(Matrix),intent(in) :: F(1),D(1),S
!  logical :: Dsym
  real(realk),pointer          ::expval(:,:)
  real(realk),pointer          ::NucSpecNMSTtotal(:,:),Prodtotal(:,:),eivalkF(:)
  real(realk)                  :: Factor
  real(realk)                  :: TS,TE
  integer                      :: natoms,icoor,jcoor,lupri,nbast,luerr,k
  Character(len=4),allocatable :: atomName(:)
!  type(Matrix)     :: NDX(3),Sx(3),tempm1,GbJ(3), GbK(3),GbDX(3)
!  type(Matrix) :: GbXc(3), GbDXc(3)
  type(Matrix),pointer ::  RHSk(:),Xk(:)
  type(Matrix) :: ProdA(3),TmpA(3),TmpB(3),tempm1
  integer :: ntrial,nrhs,nsol,nomega,nstart,nAtomsSelected,nnonZero2
  integer :: istart,iend,iAtom,kcoor,offset,iy,Xcoor,nnonZero,nnz
  integer,pointer :: AtomList(:)
  character(len=1)        :: CHRXYZ(-3:3)
  logical :: FoundNMRLabel
  logical,pointer :: NonZero(:)
  real(realk),pointer :: magderivKcont(:,:)
  DATA CHRXYZ /'z','y','x',' ','X','Y','Z'/
  Factor=53.2513539566280 !1e6*alpha^2 

  nbast = D(1)%nrow
  lupri = molcfg%lupri
  luerr = lupri
  natoms = molcfg%natoms
  CALL LSTIMER('START ',TS,TE,LUPRI)

  WRITE(LUPRI,*) '=========================================================='
  WRITE(LUPRI,*) '  NUCLEAR MAGNETIC SHIELDING(NMS) TENSOR RESULTS  (IANS)  '
  WRITE(LUPRI,*) '=========================================================='

  !#############################################################################
  !##                                                                         ## 
  !##   Input analysis                                                        ##
  !##                                                                         ## 
  !#############################################################################

  IF(ls%input%MOLECULE%nSubSystems.NE.0)THEN
     WRITE(LUPRI,'(2X,A38,2X,I7)')'number of Subsystem Labels         :',ls%input%MOLECULE%nSubSystems
     FoundNMRLabel=.FALSE.
     DO Icoor=1,ls%input%MOLECULE%nSubSystems
        WRITE(LUPRI,'(2X,A,I3,A,A)')'Subsystem Label(',Icoor,'):',TRIM(ls%input%MOLECULE%SubsystemLabel(Icoor))
        IF(TRIM(ls%input%MOLECULE%SubsystemLabel(Icoor)).EQ.'NMR')THEN
           FoundNMRLabel=.TRUE.
           nAtomsSelected = 0 
           DO Jcoor=1,ls%input%MOLECULE%nAtoms
              IF(ls%input%MOLECULE%ATOM(Jcoor)%SubSystemIndex.EQ.Icoor)THEN
                 nAtomsSelected = nAtomsSelected + 1 
              ENDIF
           ENDDO
           call mem_alloc(AtomList,nAtomsSelected)
           nAtomsSelected = 0 
           DO Jcoor=1,ls%input%MOLECULE%nAtoms
              IF(ls%input%MOLECULE%ATOM(Jcoor)%SubSystemIndex.EQ.Icoor)THEN
                 nAtomsSelected = nAtomsSelected + 1 
                 AtomList(nAtomsSelected) = Jcoor
              ENDIF
           ENDDO
        ENDIF
     ENDDO
     IF(FoundNMRLabel)THEN
        WRITE(LUPRI,'(2X,A,2X,I7)')'number of Atoms selected:',nAtomsSelected
        WRITE(LUPRI,'(2X,A,2X,I7)')'List of selected Atoms  :'
        DO Icoor=1,nAtomsSelected,10
           write(LUPRI,'(1X,10i7)') AtomList(Icoor:MIN(Icoor+9,nAtomsSelected))
        ENDDO
     ELSE
        CALL LSQUIT('IANS NUCLEAR MAGNETIC SHIELDING(NMS) TENSOR Requires SubSystem=NMR specified',-1)
     ENDIF
  ELSE
     CALL LSQUIT('IANS NUCLEAR MAGNETIC SHIELDING(NMS) TENSOR Requires SubSystems specified',-1)
  ENDIF

  call mem_alloc(NucSpecNMSTtotal,3,3*nAtomsSelected)
  call mem_alloc(Prodtotal,3,3*nAtomsSelected)
  !icoor is to be used for the 3 magnetic coordinates
  !jcoor is to be used for the 3*natoms magnetic moment coordinates

  do icoor=1,3      
     call mat_init(ProdA(icoor),nbast,nbast)
     call mat_init(TmpA(icoor),nbast,nbast)
     call mat_init(TmpB(icoor),nbast,nbast)
  enddo
  call mat_init(tempm1,nbast,nbast)
  
  call II_get_magderivOverlap(TmpA,molcfg%setting,lupri,luerr)
  !TmpA is now Sx=Sb, TmpB will be NDX = -D*SX*D
  do icoor = 1,3 
     call mat_mul(TmpA(icoor),D(1),'n','n',1E0_realk,0.E0_realk,tempm1)     !tempm1 = SX*D
     call mat_mul(D(1),tempm1,'n','n',-1.0E0_realk,0.E0_realk,TmpB(icoor))  !NDX = -D*SX*D
     !Contribution: F(NDX=D0^b)S 
     call mat_mul(TmpB(icoor),S,'n','n',1E0_realk,0E0_realk,tempm1) 
     call mat_mul(F(1),tempm1,'n','n',-1E0_realk,0E0_realk,ProdA(icoor))
     !Contribution: -S(NDX=D0^b)F
     call mat_mul(TmpB(icoor),F(1),'n','n',1E0_realk,0E0_realk,tempm1) 
     call mat_mul(S,tempm1,'n','n',1E0_realk,1E0_realk,ProdA(icoor))
     !Contribution: -FD0SB 
     call mat_mul(D(1),TmpA(icoor), 'n','n',1E0_realk,0E0_realk,tempm1) 
     call mat_mul(F(1),tempm1,'n','n',1E0_realk,1E0_realk,ProdA(icoor)) 
     !Contribution: SBD0F
     call mat_mul(D(1),F(1), 'n','n',1E0_realk,0E0_realk,tempm1) 
     call mat_mul(TmpA(icoor),tempm1,'n','n',-1E0_realk,1E0_realk,ProdA(icoor)) 
  enddo
  call mat_mul(D(1),S, 'n','n',1E0_realk,0E0_realk,tempm1)
  call II_get_prop(LUPRI,LUERR,molcfg%SETTING,TmpA,3,'MAGMOM ')   ! Generate h^b

  !TmpA is now hb; tempm1 = DS
  do icoor=1,3      
     !Contribution:  -hbDS 
     call mat_mul(TmpA(icoor),tempm1,'n','n',1E0_realk,1E0_realk,ProdA(icoor)) 
     !Contribution:  SDhb
     call mat_mul(tempm1,TmpA(icoor),'t','n',-1E0_realk,1E0_realk,ProdA(icoor)) 
  enddo
  call II_get_magderivJ(LUPRI,LUERR,molcfg%SETTING,nbast,D,TmpA)   !Generate J^b
  !TmpA is now GbK
  do icoor=1,3      
     !Contribution: J^b(D)D0S 
     call mat_mul(TmpA(icoor),tempm1,'n','n',1E0_realk,1E0_realk,ProdA(icoor)) 
     !Contribution: -SD0J^b(D)
     call mat_mul(tempm1,TmpA(icoor),'t','n',-1E0_realk,1E0_realk,ProdA(icoor)) 
  enddo
!==============================================================================
  !  This contribution is treated in a special way due to efficiency
!  do icoor=1,3      
!     call mat_zero(TmpA(icoor))
!  enddo
!  call II_get_magderivK(LUPRI,LUERR,molcfg%SETTING,nbast,D,TmpA)   !Add K^b 
!  !TmpA is now GbK
!  do icoor=1,3      
!     !Contribution: K^b(D)D0S      
!     call mat_mul(TmpA(icoor),tempm1,'n','n',1E0_realk,1E0_realk,ProdA(icoor)) 
!     !Contribution: - SD0K^b(D) 
!     call mat_mul(tempm1,TmpA(icoor),'t','n',-1E0_realk,1E0_realk,ProdA(icoor)) 
!  enddo
!==============================================================================
  call di_GET_GbDs(lupri,luerr,TmpB,TmpA,3,molcfg%setting)  ! G(D0^B)
  !TmpA is now GbDX = G(D0^B)
  do icoor=1,3      
     !Contribution: G(D0B)D0S
     call mat_mul(TmpA(icoor),tempm1,'n','n',1E0_realk,1E0_realk,ProdA(icoor)) 
     !Contribution: -SD0G(D0B)
     call mat_mul(tempm1,TmpA(icoor),'t','n',-1E0_realk,1E0_realk,ProdA(icoor)) 
  enddo

  IF(molcfg%setting%do_dft)THEN
     call II_get_xc_linrsp(lupri,luerr,molcfg%setting,nbast,TmpB,D(1),TmpA,3) 
     !TmpA is now GbDXc = Gxc(D0^B)
     do icoor=1,3      
        !Contribution: G^bDOS (DFT XC) - SD0G^b (DFT XC)
        call mat_mul(TmpA(icoor),tempm1,'n','n',1E0_realk,1E0_realk,ProdA(icoor)) 
        call mat_mul(tempm1,TmpA(icoor),'t','n',-1E0_realk,1E0_realk,ProdA(icoor)) 
     enddo
     call II_get_xc_magderiv_kohnsham_mat(LUPRI,LUERR,molcfg%SETTING,nbast,D(1),TmpA)
     !TmpA is now GbXc
     do icoor=1,3      
        !Contribution: G^bDOS (DFT XC) - SD0G^b (DFT XC)
        call mat_mul(TmpA(icoor),tempm1,'n','n',1E0_realk,1E0_realk,ProdA(icoor)) 
        call mat_mul(tempm1,TmpA(icoor),'t','n',-1E0_realk,1E0_realk,ProdA(icoor)) 
     enddo
  ENDIF
  do icoor=1,3
     call mat_free(tmpB(icoor))
  enddo 

  IF(ls%input%DALTON%SolveNMRResponseSimultan)THEN
     !#############################################################################
     !##                                                                         ## 
     !##   Building D^k (Density matrix derivative w.r.t nuclei magnetic moment) ##
     !##     Calc RHS of eq 70 JCP 115 page 10349                                ##
     !##      RHS in this case;                                                  ##
     !##      ! RHSk = -P_anti(hkDS)                                             ##
     !##                                                                         ## 
     !#############################################################################
     ! Generation of hkDS      
     allocate(RHSk(3*nAtomsSelected))   ! RHSk - R.H.S. to solve eq. 70 for D^K 
     !tempm1 = DS
     do jcoor=1,nAtomsSelected
        istart=3*jcoor-2
        iend=istart+2
        iAtom = AtomList(jcoor)
        call II_get_PSO_spec(LUPRI,LUERR,molcfg%SETTING,TmpA,iAtom)
        !TmpA(1:3) is now hk(istart:iend)
        iy = istart+1
        !X coordinate
        call mat_init(RHSk(istart),nbast,nbast)
        call mat_mul(TmpA(1),tempm1,'n','n',1.0E0_realk,0.0E0_realk,RHSk(istart))  
        call mat_mul(tempm1,TmpA(1),'t','n',1.0E0_realk,-1.0E0_realk,RHSk(istart))  
        !Y coordinate
        call mat_init(RHSk(iy),nbast,nbast)
        call mat_mul(TmpA(2),tempm1,'n','n',1.0E0_realk,0.0E0_realk,RHSk(iy))  
        call mat_mul(tempm1,TmpA(2),'t','n',1.0E0_realk,-1.0E0_realk,RHSk(iy))  
        !Z coordinate
        call mat_init(RHSk(iend),nbast,nbast)
        call mat_mul(TmpA(3),tempm1,'n','n',1.0E0_realk,0.0E0_realk,RHSk(iend))  
        call mat_mul(tempm1,TmpA(3),'t','n',1.0E0_realk,-1.0E0_realk,RHSk(iend))  
     enddo
     do icoor=1,3
        call mat_free(tmpA(icoor))
     enddo
     !###########################################################################
     !Solve K([D,Xk]s) = RHSk is now ready                                
     !Solve to get Dk=[D,Xk]_S + D0^k --> Xk (Eq. 70)                    
     !We can use the RSP solver, since K([D,Xk]s) = sigma = -1/2 * E[2]Xk
     !
     !############################################################################     
     call mem_alloc(eivalkF,3*nAtomsSelected)
     eivalkF=0.0E0_realk
     write(lupri,*)'Calling rsp solver for all Xk  '
     allocate(Xk(3*nAtomsSelected))   
     nnonZero = 0
     call mem_alloc(NonZero,3*nAtomsSelected)
     do jcoor=1,3*nAtomsSelected
        call mat_init(Xk(jcoor),nbast,nbast)                           
        call util_scriptPx('T',D(1),S,RHSk(jcoor))     
        if ( mat_dotproduct(RHSk(jcoor),RHSk(jcoor)).LT.1.0d-10) then
           print*,'WARNING RHS(jcoor=',jcoor,') = Zero '
           NonZero(jcoor) = .FALSE.
        else
           NonZero(jcoor) = .TRUE.
           nnonZero = nnonZero + 1
        endif
     enddo
     ntrial = 3*nAtomsSelected !# of trial vectors in a given iteration (number of RHS)
     nrhs = 3*nAtomsSelected   !# of RHS only relevant for linear equations (lineq_x = TRUE)
     nsol = 3*nAtomsSelected   !# of solution (output) vectors
     nomega = 3*nAtomsSelected !If lineq_x, number of laser freqs (input)
     !Otherwise number of excitation energies (output) 
     nstart = 3*nAtomsSelected !Number of start vectors. Only relevant for eigenvalue problem
     !ntrial and nstart seem to be obsolete 
     call rsp_init(ntrial,nrhs,nsol,nomega,nstart)
     if (molcfg%solver%info_rsp_sparsity) then
        do jcoor=1,3*nAtomsSelected
           call mat_report_sparsity(RHSk(jcoor),'NMRshield IANS RHS',nnz,lupri)
           Write(lupri,'(A,I3,A,I9)')'RHSk(',jcoor,') NNZ=',NNZ
        enddo
     endif
     call rsp_solver(molcfg,D(1),S,F(1),.true.,nrhs,RHSk,EIVALKF,Xk)
     do jcoor=1,3*nAtomsSelected
        do icoor = 1,3 
           Prodtotal(icoor,jcoor)=-4E0_realk*factor*mat_trAB(Xk(jcoor),ProdA(icoor))
        enddo
     enddo
     do jcoor=1,3*nAtomsSelected
        call mat_free(RHSk(jcoor))
     enddo
     call mem_dealloc(eivalkF)
     deallocate(RHSk)
     !Contribution Tr(K^b(D)*[D,X]_S) ( the missing term to ProdA: K^b(D)D0S - SD0K^b(D) )
     !First: DSX-XSD  
     nnonZero2 = 0
     do jcoor=1,3*nAtomsSelected
        IF(NonZero(jcoor))THEN !This should only be done for non Zero X's 
           nnonZero2 = nnonZero2 + 1
           call mat_mul(tempm1,Xk(jcoor),'n','n',1E0_realk,0E0_realk,ProdA(1)) 
           call mat_mul(Xk(jcoor),tempm1,'n','t',-1E0_realk,1E0_realk,ProdA(1)) 
           call mat_assign(Xk(nnonZero2),ProdA(1))
        ENDIF
     enddo
     IF(nnonZero.NE.nnonZero2)call lsquit('dimmismatch for nuclear selected MagderivKcont',-1)
     do jcoor=nnonZero+1,3*nAtomsSelected
        call mat_free(Xk(jcoor))
     enddo
     !Xk is now [D,X]_S a contribution to the perturbed density matrix
     call mem_alloc(magderivKcont,3,nnonZero)
     call II_get_magderivKcont(LUPRI,LUERR,molcfg%SETTING,Xk,&
          & nnonZero,magderivKcont,D)

     nnonZero2 = 0
     do jcoor=1,3*nAtomsSelected
        IF(NonZero(jcoor))THEN
           nnonZero2 = nnonZero2 + 1
           do icoor=1,3
              Prodtotal(icoor,jcoor)=Prodtotal(icoor,jcoor)+&
                   & 4E0_realk*factor*magderivKcont(icoor,nnonZero2)
           enddo
        ENDIF
     enddo     
     call mem_dealloc(NonZero)
     call mem_dealloc(magderivKcont)
     do jcoor=1,nnonZero2
        call mat_free(Xk(jcoor))
     enddo
     deallocate(Xk)     
  ELSE
     !#############################################################################
     !##                                                                         ## 
     !##   Building D^k (Density matrix derivative w.r.t nuclei magnetic moment) ##
     !##   Calc RHS of eq 70 JCP 115 page 10349                                  ##
     !##   RHS in this case;    -P_anti(hkDS)                                    ##
     !#############################################################################
     ! Generation of hkDS      
     allocate(RHSk(1))   ! RHSk - R.H.S. to solve eq. 70 for D^K      
     !tempm1 = DS
     allocate(Xk(3*nAtomsSelected))   
     call mem_alloc(eivalkF,1)
     call mat_init(RHSk(1),nbast,nbast)
     call mem_alloc(NonZero,3*nAtomsSelected)
     eivalkF=0.0E0_realk
     nNonZero = 0 
     do jcoor=1,nAtomsSelected
        istart=3*jcoor-2
        iend=istart+2
        iAtom = AtomList(jcoor)
        call II_get_PSO_spec(LUPRI,LUERR,molcfg%SETTING,TmpA,iAtom)
        !TmpA(1:3) is now the hk(istart:iend)
        do Xcoor=1,3
           !-hk*DS + SD*hk
           call mat_mul(TmpA(Xcoor),tempm1,'n','n',1.0E0_realk,0.0E0_realk,RHSk(1))  
           call mat_mul(tempm1,TmpA(Xcoor),'t','n',1.0E0_realk,-1.0E0_realk,RHSk(1))  
           !###########################################################################
           !Solve K([D,Xk]s) = RHSk is now ready                                
           !Solve to get Dk=[D,Xk]_S + D0^k --> Xk (Eq. 70)                    
           !We can use the RSP solver, since K([D,Xk]s) = sigma = -1/2 * E[2]Xk
           !############################################################################     
           write(lupri,*)'Calling rsp solver for Xk  '
           call util_scriptPx('T',D(1),S,RHSk(1))
           if ( mat_dotproduct(RHSk(1),RHSk(1))>1.0d-10) then
              NonZero(Xcoor+(jcoor-1)*3) = .TRUE.
              nNonZero = nNonZero + 1
              ntrial = 1 !# of trial vectors in a given iteration (number of RHS)
              nrhs = 1   !# of RHS only relevant for linear equations (lineq_x = TRUE)
              nsol = 1   !# of solution (output) vectors
              nomega = 1 !If lineq_x, number of laser freqs (input)
              !Otherwise number of excitation energies (output) 
              nstart = 1 !Number of start vectors. Only relevant for eigenvalue problem
              !ntrial and nstart seem to be obsolete 
              call rsp_init(ntrial,nrhs,nsol,nomega,nstart)
              call mat_init(Xk(nNonZero),nbast,nbast)
              if (molcfg%solver%info_rsp_sparsity) then
                 call mat_report_sparsity(RHSk(1),'NMRshield IANS RHS(1)',nnz,lupri)
                 Write(lupri,*)'RHSk(',Xcoor+(jcoor-1)*3,') NNZ=',NNZ
              endif
              call rsp_solver(molcfg,D(1),S,F(1),.true.,nrhs,RHSk(1:1),EIVALKF,Xk(nNonZero:nNonZero))
              do icoor = 1,3 
                 Prodtotal(icoor,Xcoor+(jcoor-1)*3)=-4E0_realk*factor*mat_trAB(Xk(nNonZero),ProdA(icoor))
              enddo
           else
              NonZero(Xcoor+(jcoor-1)*3) = .FALSE.
              write(lupri,*) 'WARNING: RHSk norm is less than threshold'
              write(lupri,*) 'LIN RSP equations NOT solved for this RHS    '
              do icoor = 1,3 
                 Prodtotal(icoor,Xcoor+(jcoor-1)*3)=0E0_realk
              enddo
           end if
        enddo
     enddo
     call mat_free(RHSk(1))
     deallocate(RHSk)
     call mem_dealloc(eivalkF)
     do icoor = 1,3
        call mat_free(TmpA(icoor))
     enddo
     !Contribution Tr(K^b(D)*[D,X]_S) ( the missing term to ProdA: K^b(D)D0S - SD0K^b(D) )
     !First: DSX-XSD  
     do jcoor=1,nnonZero
        call mat_mul(tempm1,Xk(jcoor),'n','n',1E0_realk,0E0_realk,ProdA(1)) 
        call mat_mul(Xk(jcoor),tempm1,'n','t',-1E0_realk,1E0_realk,ProdA(1)) 
        call mat_assign(Xk(jcoor),ProdA(1))
     enddo
     !Xk is now [D,X]_S a contribution to the perturbed density matrix
     call mem_alloc(magderivKcont,3,nnonZero)
     call II_get_magderivKcont(LUPRI,LUERR,molcfg%SETTING,Xk(1:nnonZero),&
          & nnonZero,magderivKcont,D)
     nnonZero2 = 0
     do jcoor=1,3*nAtomsSelected
        IF(NonZero(jcoor))THEN
           nnonZero2 = nnonZero2 + 1
           do icoor=1,3
              Prodtotal(icoor,jcoor)=Prodtotal(icoor,jcoor)+&
                   & 4E0_realk*factor*magderivKcont(icoor,nnonZero2)
           enddo     
        ENDIF
     enddo
     call mem_dealloc(magderivKcont)
     do jcoor=1,nnonZero
        call mat_free(Xk(jcoor))
     enddo
     deallocate(Xk)     
     call mem_dealloc(NonZero)
  ENDIF
  do icoor=1,3
    call mat_free(ProdA(icoor))
  enddo 
  write(lupri,*) '---------------------------------------------------------'
  write(lupri,*) 'NUCLEI MAGNETIC TENSOR  | Response (Xk) test             '
  write(lupri,*) '---------------------------------------------------------'
     
  !#############################################################################
  !     Tr D*h^kb 
  ! 
  !#############################################################################
  !expval(3,3*NATOM) (magnetic X,Y,Z, atomic moment X,Y,Z, for each Atom)

  call mem_alloc(expval,3,3)  
  do jcoor=1,nAtomsSelected
     offset = (jcoor-1)*3
     iAtom = AtomList(jcoor)
     call II_get_nst_spec_expval(LUPRI,LUERR,molcfg%SETTING,expval,D,iAtom)
     do icoor=1,3     ! magnetic koordinate
        do kcoor=1,3  ! magnetic moment koordinate
           NucSpecNMSTtotal(icoor,offset+kcoor) = 2*factor*expval(icoor,kcoor)
        enddo
     enddo
  enddo
  call mem_dealloc(expval)
  !############################
  !Print D*h^kb- diamagnetic part
  !
  !###############################
  allocate(atomname(natoms))
  do jcoor=1,natomsselected  
    iAtom = AtomList(jcoor)
    atomname(jcoor)=molcfg%setting%molecule(1)%p%ATOM(iAtom)%Name
  enddo
  
  WRITE(LUPRI,*) "      D*h^kb - diamagnetic term"
  do jcoor=1,nAtomsSelected 
       WRITE (LUPRI,'(2X,A,I7)') 'Atom Identity:',AtomList(jcoor)
       WRITE (LUPRI,'(2X,A8,f15.8,A8)')  atomname(jcoor), &
        & 1.0E0_realk/3.0E0_realk*(NucSpecNMSTtotal(1,3*jCOOR-2)+NucSpecNMSTtotal(2,3*jCOOR-1) &
        & +NucSpecNMSTtotal(3,3*jCOOR))
  enddo

  !#################################################################################
  !                   Calculate (h^k)*(D0^b)
  !
  !#################################################################################
     
  do icoor=1,3      
     call mat_init(TmpA(icoor),nbast,nbast)
     call mat_init(TmpB(icoor),nbast,nbast)
  enddo
  call II_get_magderivOverlap(TmpB,molcfg%setting,lupri,luerr)
  do icoor=1,3     ! magnetic
     call mat_mul(TmpB(icoor),D(1),'n','n',1E0_realk,0.E0_realk,tempm1)     !tempm1 = SX*D
     call mat_mul(D(1),tempm1,'n','n',-1.0E0_realk,0.E0_realk,TmpA(icoor))       !tempm2 = -D*SX*D
  enddo
  do jcoor=1,nAtomsSelected
     istart=3*jcoor-2
     iend=istart+2
     iAtom = AtomList(jcoor)
     call II_get_PSO_spec(LUPRI,LUERR,molcfg%SETTING,TmpB,iAtom)          
     !TmpA(1:3) is now -D*SX*D and TmpB(1:3) is hk(istart:iend)
     do icoor=1,3     ! magnetic
        !X magnetic moment koordinate        
        NucSpecNMSTtotal(icoor,1+(jcoor-1)*3) = NucSpecNMSTtotal(icoor,1+(jcoor-1)*3) &
             & - 2.0E0_realk*factor*mat_dotproduct(TmpB(1),TmpA(icoor))
        !Y magnetic moment koordinate        
        NucSpecNMSTtotal(icoor,2+(jcoor-1)*3) = NucSpecNMSTtotal(icoor,2+(jcoor-1)*3) &
             & - 2.0E0_realk*factor*mat_dotproduct(TmpB(2),TmpA(icoor))
        !Z magnetic moment koordinate        
        NucSpecNMSTtotal(icoor,3+(jcoor-1)*3) = NucSpecNMSTtotal(icoor,3+(jcoor-1)*3) &
             & - 2.0E0_realk*factor*mat_dotproduct(TmpB(3),TmpA(icoor))
     enddo
  enddo
  do icoor=1,3      
     call mat_free(TmpA(icoor))
     call mat_free(TmpB(icoor))
  enddo
  call mat_free(tempm1)

  !####################################################
  !Print D*h^kb + h^k D0^b
  !
  !####################################################
  WRITE(LUPRI,*) "      D*h^kb + h^k D0^b "
  do jcoor=1,nAtomsSelected 
       WRITE (LUPRI,'(2X,A,I7)') 'Atom Identity:',AtomList(jcoor)
       WRITE (LUPRI,'(2X,A8,f15.8,A8)')  atomname(jcoor), &
        & 1.0E0_realk/3.0E0_realk*(NucSpecNMSTtotal(1,3*jCOOR-2)+NucSpecNMSTtotal(2,3*jCOOR-1) &
        & +NucSpecNMSTtotal(3,3*jCOOR))
  enddo
  
  !#######################################################
  !Print D^(B,k)*F+D^k*G^b*D+D^k*G*D0^b+D^k*h^b
  !
  !######################################################


  WRITE(LUPRI,*) "      D^(B,k)*F+D^k*G^b*D+D^k*G*D0^b+D^k*h^b "
  do jcoor=1,nAtomsSelected 
       WRITE (LUPRI,'(2X,A,I7)') 'Atom Identity:',AtomList(jcoor)
       WRITE (LUPRI,'(2X,A8,f15.8,A8)')  atomname(jcoor), &
        & 1.0E0_realk/3.0E0_realk*(Prodtotal(1,3*jCOOR-2)+Prodtotal(2,3*jCOOR-1) &
        & +Prodtotal(3,3*jCOOR))
  enddo
 !######################################################################################
 !  Adding: Prodtotal+ NucSpecNMSTtotal
 !
 !#####################################################################################
 
  do icoor=1,3     ! magnetic
    do jcoor=1,3*NATOMSselected ! magnetic moment koordinate
       NucSpecNMSTtotal(icoor,jcoor) = NucSpecNMSTtotal(icoor,jcoor) &
            & + Prodtotal(icoor,jcoor)
    enddo
  enddo

!  allocate(atomname(natoms))
!  do jcoor=1,natomsselected  
!    iAtom = AtomList(jcoor)
!    atomname(jcoor)=molcfg%setting%molecule(1)%p%ATOM(iAtom)%Name
!  enddo
 
  Write (lupri,*)   "The nuclear magnetic shielding for selected atoms" 
  do jcoor=1,natomsselected  
    WRITE (LUPRI,'(2X,A,I7)') 'Atom Identity:',AtomList(jcoor)
    WRITE (LUPRI,'(20X,3(A,13X),/)') 'Bx', 'By', 'Bz'
    WRITE (LUPRI,'(2X,A8,3F15.8)')  atomname(jcoor)//CHRXYZ(1),(NucSpecNMSTtotal(K,3*jCOOR-2),K=1,3)
    WRITE (LUPRI,'(2X,A8,3F15.8)')  atomname(jcoor)//CHRXYZ(2),(NucSpecNMSTtotal(K,3*jCOOR-1),K=1,3)
    WRITE (LUPRI,'(2X,A8,3F15.8)')  atomname(jcoor)//CHRXYZ(3),(NucSpecNMSTtotal(K,3*jCOOR),K=1,3)
  enddo

  WRITE(LUPRI,*) "      Absolute chemical shift"
  do jcoor=1,nAtomsSelected 
       WRITE (LUPRI,'(2X,A,I7)') 'Atom Identity:',AtomList(jcoor)
       WRITE (LUPRI,'(2X,A8,f15.8,A8)')  atomname(jcoor), &
        & 1.0E0_realk/3.0E0_realk*(NucSpecNMSTtotal(1,3*jCOOR-2)+NucSpecNMSTtotal(2,3*jCOOR-1) &
        & +NucSpecNMSTtotal(3,3*jCOOR))
  enddo
  call mem_dealloc(NucSpecNMSTtotal)
  call mem_dealloc(Prodtotal)
  deallocate(atomname)
  call mem_dealloc(AtomList)
  WRITE(LUPRI,*) " Done with shielding tensor calculation"

 end subroutine NMRshieldresponse_RSTNS
end module nuclei_selected_shielding_mod 
