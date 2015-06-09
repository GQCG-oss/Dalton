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

  public :: NMRshieldresponse_IANS_TK, NMRshieldresponse_NS, &
            & NMRshieldresponse_RSTNS 

  private

Contains   
           
           
  !####################################################################
  !
  ! Subroutine - NMRshieldresponse_IANS 
  ! It computes Shileding tensor for nuclei magnetic moment response Xk.
  ! However, It compute for all the atoms. HF, and DFT.  
  !
  !#####################################################################


subroutine NMRshieldresponse_IANS_TK(ls,molcfg,F,D,S)
  implicit none
  TYPE(lsitem),target     :: ls
  type(rsp_molcfg), intent(inout) :: molcfg
  type(Matrix),intent(in) :: F(1),D(1),S
  !
  logical :: Dsym
  real(realk),pointer          ::NMST(:,:),expval(:,:),NMST1(:,:),NMST2(:,:),NMST3(:,:),NMST4(:,:)
  real(realk),pointer          :: DSX4(:,:),DSX6(:,:),NucSpecNMSTtotal(:,:)
  real(realk)                  :: Factor,eival(1),eivalk(1)
  integer                      :: natoms,icoor,jcoor,k,lupri,nbast,luerr,n_rhs
  Character(len=4),allocatable :: atomName(:)
  real(realk)                  :: TS,TE
  type(Matrix)     :: Dx(3),NDX(3),Fx(3),Sx(3),tempm1,RHS(3),GbDs(3),Xx(1),Xk(1),GbJ(3), GbK(3),Gm(3),tf1,tf2, tempmx,GbDX(3)
  type(Matrix) :: GbXc(3)
  type(Matrix),pointer ::  RHSk(:),  DXk(:),hk(:),hb(:),GkD(:),tempk(:),sdf(:),DXD(:), DX0(:), DXk2np1A(:,:), DXk2np1B(:,:)
  type(matrix) :: Prod,Prod1
  integer              :: ntrial,nrhs,nsol,nomega,nstart,nAtomsSelected
  integer,pointer :: AtomList(:)
  character(len=1)        :: CHRXYZ(-3:3)
  logical :: FoundNMRLabel
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



  !icoor is to be used for the 3 magnetic coordinates
  !jcoor is to be used for the 3*natoms magnetic moment coordinates

  !#############################################################################
  !##                                                                         ## 
  !##   Building D^k (Density matrix derivative w.r.t nuclei magnetic moment) ##
  !##     Calc RHS of eq 70 JCP 115 page 10349                                ##
  !##      RHS in this case;                                                  ##
  !##      ! RHSk = -P_anti(hkDS)                                             ##
  !##                                                                         ## 
  !#############################################################################
 
  ! Generation of hkDS      
  allocate(hk(3*natoms))     ! Matrix h^k: Generation by 'PSO    ' 
  allocate(RHSk(3*natoms))   ! RHSk - R.H.S. to solve eq. 70 for D^K 
  allocate(Dxk(3*natoms))    ! (D^K Derivative of Density matrix) 
  
  do jcoor = 1,3*NATOMS
     call mat_init(hk(jcoor),nbast,nbast)
  enddo
  call II_get_prop(LUPRI,LUERR,molcfg%SETTING,hk,3*NATOMS,'PSO    ')   !Generation of h^k  

  WRITE(lupri,*)'hk icoor=1 '
  call mat_print(hk(1),1,hk(1)%nrow,1,hk(1)%ncol,lupri)
  
  call mat_init(tempm1,nbast,nbast) 
  call mat_mul(D(1),S,'n','n',1.0E0_realk,0.0E0_realk,tempm1)
  do jcoor = 1,3*NATOMS
     call mat_init(RHSk(jcoor),nbast,nbast)
     !Generation of RHSk = h^kDS  
     call mat_mul(hk(jcoor),tempm1,'n','n',1.0E0_realk,0.0E0_realk,RHSk(jcoor))  
     call mat_mul(tempm1,hk(jcoor),'t','n',1.0E0_realk,-1.0E0_realk,RHSk(jcoor))  
     call mat_free(hk(jcoor))
  enddo  
  call mat_free(tempm1)
   
  !###########################################################################
  !Solve K([D,Xk]s) = RHSk is now ready                                
  !Solve to get Dk=[D,Xk]_S + D0^k --> Xk (Eq. 70)                    
  !We can use the RSP solver, since K([D,Xk]s) = sigma = -1/2 * E[2]Xk
  !
  !############################################################################
     
  call mat_init(Xk(1),nbast,nbast)                           
  allocate(DXk2np1A(3*natoms,3))
  allocate(DXk2np1B(3*natoms,3))
  do jcoor=1,3*natoms
     eivalk(1)=0.0E0_realk
     write(lupri,*)'Calling rsp solver for Xk  '
     call util_scriptPx('T',D(1),S,RHSk(jcoor))
     if ( mat_dotproduct(RHSk(jcoor),RHSk(jcoor))>1.0d-10) then

        ntrial = 1 !# of trial vectors in a given iteration (number of RHS)
        nrhs = 1   !# of RHS only relevant for linear equations (lineq_x = TRUE)
        nsol = 1   !# of solution (output) vectors
        nomega = 1 !If lineq_x, number of laser freqs (input)
                   !Otherwise number of excitation energies (output) 
        nstart = 1 !Number of start vectors. Only relevant for eigenvalue problem
        !ntrial and nstart seem to be obsolete 
        call rsp_init(ntrial,nrhs,nsol,nomega,nstart)
        call rsp_solver(molcfg,D(1),S,F(1),.true.,nrhs,RHSk(jcoor:jcoor),EIVALK,Xk)
     else
        write(lupri,*) 'WARNING: RHSk norm is less than threshold'
        write(lupri,*) 'LIN RSP equations NOT solved for this RHS    '
        call mat_zero(Xk(1))
     end if
     call mat_free(RHSk(jcoor))

     call mat_init(Prod,nbast,nbast)
     call mat_init(Prod1,nbast,nbast)
     call mat_init(SX(1),nbast,nbast)
     call mat_init(SX(2),nbast,nbast)
     call mat_init(SX(3),nbast,nbast)                                     
     call II_get_magderivOverlap(Sx,molcfg%setting,lupri,luerr)
     call mat_init(tempm1,nbast,nbast)
     do icoor = 1,3 
        call mat_mul(SX(icoor),D(1),'n','n',1E0_realk,0.E0_realk,tempm1)     !tempm1 = SX*D
        call mat_init(NDX(icoor),nbast,nbast)        
        call mat_mul(D(1),tempm1,'n','n',-1.0E0_realk,0.E0_realk,NDX(icoor))  !NDX = -D*SX*D
        
        ![D0b,Xk]_S
        call mat_init(DXk2np1A(jcoor,icoor),nbast,nbast)
        call mat_init(DXk2np1B(jcoor,icoor),nbast,nbast)
        call mat_mul(NDX(icoor),S,'n','n',1E0_realk,0E0_realk,Prod)    
        call mat_mul(Prod,Xk(1),'n','n',1E0_realk,0E0_realk,DXk2np1A(jcoor,icoor))  
        call mat_mul(S,NDX(icoor),'n','n',1E0_realk,0E0_realk,Prod)    
        call mat_mul(Xk(1),Prod,'n','n',1E0_realk,0E0_realk,Prod1) 
        call mat_daxpy(-1E0_realk,Prod1,DXk2np1A(jcoor,icoor)) 
        call mat_free(NDX(icoor))
        
        !+[D0,Xk]_Sb
        call mat_mul(D(1),SX(icoor),'n','n',1E0_realk,0E0_realk,Prod)    
        call mat_mul(Prod,Xk(1),'n','n',1E0_realk,0E0_realk,DXk2np1B(jcoor,icoor))  
        call mat_mul(SX(icoor),D(1),'n','n',1E0_realk,0E0_realk,Prod)    
        call mat_mul(Xk(1),Prod,'n','n',1E0_realk,0E0_realk,Prod1) 
        call mat_daxpy(-1E0_realk,Prod1,DXk2np1B(jcoor,icoor)) 
        call mat_free(SX(icoor))                    
     enddo
     call mat_free(tempm1)
     call mat_free(Prod)
     call mat_free(Prod1)
     
     !Make D^k=[D_0,X^k]_s
     call mat_init(DXk(jcoor),nbast,nbast)
     call ABCcommutator(nbast,D(1),Xk(1),S,DXk(jcoor))       !Generate [D_0,X^k]_s 
     call mat_scal(-4.0d0,DXk(jcoor))
  enddo
  call mat_free(Xk(1)) 

  write(lupri,*) '---------------------------------------------------------'
  write(lupri,*) 'NUCLEA MAGNETIC TENSOR (INASHIELD)                                  '
  write(lupri,*) '---------------------------------------------------------'

  call mem_alloc(NucSpecNMSTtotal,3,3*NATOMS)

  !#############################################################################
  !    Tr D*h^kb 
  ! 
  !#############################################################################

  !expval(3,3*NATOM) (magnetic X,Y,Z, atomic moment X,Y,Z, for each Atom)
  call mem_alloc(expval,3,3*NATOMS)
  call II_get_prop_expval(LUPRI,LUERR,molcfg%SETTING,expval,D,1,9*NATOMS,'NST    ')
  ! Printing  (D^*(h^kb) in output file  over (icoor,jcoor) 
  WRITE(lupri,*)  'icoor', 'jcoor', 'X1, |PURE (D)*(h^kb)|  :  2*factor*expval(jcoor,icoor)'  
  do icoor=1,3
     do jcoor=1,3*natoms  ! magnetic moment koordinate
        WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'expval1:',   2*factor*expval(icoor,jcoor)
        NucSpecNMSTtotal(icoor,jcoor) = 2*factor*expval(icoor,jcoor)
     enddo
  enddo
  call mem_dealloc(expval)

  !###########################################################################
  !    Tr D^k*h^b
  !
  !###########################################################################
  
  allocate(hb(3))   ! Allocation to matrix h^b 
  do icoor = 1,3
     call mat_init(hb(icoor),nbast,nbast)
  enddo
   
  call II_get_prop(LUPRI,LUERR,molcfg%SETTING,hb,3,'MAGMOM ')   ! Generate h^b  

  ! Printing  (D^k)*(h^b) in output file  over (icoor,jcoor) 
  WRITE(lupri,*)  'icoor', 'jcoor', 'X1, |PURE (D^k)*(h^b)|  :  - factor*mat_trAB(DXk(jcoor),hb(icoor))'  
  do icoor=1,3
     do jcoor=1,3*natoms ! magnetic moment coordinate
       WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'expval2:',+ factor*mat_trAB(DXk(jcoor),hb(icoor))
       NucSpecNMSTtotal(icoor,jcoor) = NucSpecNMSTtotal(icoor,jcoor) &
            & + factor*mat_trAB(DXk(jcoor),hb(icoor))
     enddo
     call mat_free(hb(icoor))
  enddo

  !###########################################################################
  !    Additional Term for Shielding Tensor :   Tr (D^k)*(G^b(D)) = X5
  !
  !###########################################################################
  do icoor=1,3 
     call mat_init(GbJ(icoor),nbast,nbast)
     call mat_init(GbK(icoor),nbast,nbast)
     call mat_zero(GbK(icoor))
     call mat_init(Gm(icoor),nbast,nbast) 
  enddo
  call II_get_magderivJ(LUPRI,LUERR,molcfg%SETTING,nbast,D,GbJ)   !Generate J^b
  call II_get_magderivK(LUPRI,LUERR,molcfg%SETTING,nbast,D,GbK)   !Generate K^b 
  ! Printing  (D^k)*(G^b(D)) in output file  over (icoor,jcoor) 
  WRITE(lupri,*)  'icoor', 'jcoor', ' X5 |PURE (D^K)*(J^b(D))|  '
  do icoor=1,3
     do jcoor=1,3*natoms
        WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'NMST2J:', factor*mat_trAB(DXk(jcoor),GbJ(icoor)) 
     enddo
  enddo
  WRITE(lupri,*)  'icoor', 'jcoor', ' X5 |PURE (D^K)*(K^b(D))|  '
  do icoor=1,3
     do jcoor=1,3*natoms
        WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'NMST2K:', factor*mat_trAB(DXk(jcoor),GbK(icoor)) 
        !Generate G^b= J^b+K^b
        call mat_add(1.0E0_realk,GbJ(icoor),1.0E0_realk,GbK(icoor),Gm(icoor)) 
        NucSpecNMSTtotal(icoor,jcoor) = NucSpecNMSTtotal(icoor,jcoor) &
             & + factor*mat_trAB(DXk(jcoor),Gm(icoor)) 
     enddo
     call mat_free(GbJ(icoor))
     call mat_free(GbK(icoor))
     call mat_free(Gm(icoor))
  enddo

  if(molcfg%setting%do_dft)then
  
  do icoor=1,3 
     call mat_init(GbXc(icoor),nbast,nbast)
  enddo
  call II_get_xc_magderiv_kohnsham_mat(LUPRI,LUERR,molcfg%SETTING,nbast,D(1),GbXc)
  ! Printing  (D^k)*(G^b(D)) in output file  over (icoor,jcoor) 
  WRITE(lupri,*)  'icoor', 'jcoor', ' X5 |PURE (D^K)*(G^bxc(D))|  '
  do icoor=1,3
     do jcoor=1,3*natoms
        WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'NMST2K:', factor*mat_trAB(DXk(jcoor),GbXc(icoor)) 
       
        NucSpecNMSTtotal(icoor,jcoor) = NucSpecNMSTtotal(icoor,jcoor) &
             & + factor*mat_trAB(DXk(jcoor),GbXc(icoor)) 
     enddo
     call mat_free(GbXc(icoor))
  enddo
  endif
  !#################################################################################
  !      Generate D0X = -D*SX*D
  !
  !#################################################################################
  call mat_init(SX(1),nbast,nbast)
  call mat_init(SX(2),nbast,nbast)
  call mat_init(SX(3),nbast,nbast)
  call II_get_magderivOverlap(Sx,molcfg%setting,lupri,luerr)
  call mat_init(tempm1,nbast,nbast)
  do icoor = 1,3 
     call mat_mul(SX(icoor),D(1),'n','n',1E0_realk,0.E0_realk,tempm1)     !tempm1 = SX*D
     call mat_free(SX(icoor))                    
     call mat_init(NDX(icoor),nbast,nbast)        
     call mat_mul(D(1),tempm1,'n','n',-1.0E0_realk,0.E0_realk,NDX(icoor))  !NDX = -D*SX*D
  enddo
  call mat_free(tempm1)

  !#################################################################################
  !      Debugging the term X6: Tr (D^k)*(G(D0^b)) 
  !
  !#################################################################################

  call mat_init(GbDX(1),nbast,nbast)
  call mat_init(GbDX(2),nbast,nbast)  
  call mat_init(GbDX(3),nbast,nbast)                                  
  call di_GET_GbDs(lupri,luerr,NDX,GbDX,3,molcfg%setting)  ! G(D0^B)
  WRITE(lupri,*)  'icoor', 'jcoor', 'DEBUGGING: X6 |PURE (D^k)*(G(D^b))| &
  & DSX6:  - factor*mat_trAB(DXk(jcoor),GbDX(icoor))'
  do icoor=1,3
    do jcoor=1,3*natoms
       WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'DSX6:',  factor*mat_trAB(DXk(jcoor),GbDX(icoor))
       NucSpecNMSTtotal(icoor,jcoor) = NucSpecNMSTtotal(icoor,jcoor) &
            & + factor*mat_trAB(DXk(jcoor),GbDX(icoor))
    enddo
    call mat_free(GbDX(icoor))
 enddo

 IF(molcfg%setting%do_dft)THEN
 call mat_init(GbDX(1),nbast,nbast)
 call mat_init(GbDX(2),nbast,nbast)  
 call mat_init(GbDX(3),nbast,nbast)                                  
 call II_get_xc_linrsp(lupri,luerr,molcfg%setting,nbast,NDX,D(1),GbDX,3) 
     
 WRITE(lupri,*)  'icoor', 'jcoor', 'DEBUGGING: X6- XC inculded &
  & |PURE (D^k)*(G(D^b))|  DSX6:  - factor*mat_trAB(DXk(jcoor),GbDX(icoor))'
  do icoor=1,3
     do jcoor=1,3*natoms
       WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'DSX6:',  factor*mat_trAB(DXk(jcoor),GbDX(icoor))
       NucSpecNMSTtotal(icoor,jcoor) = NucSpecNMSTtotal(icoor,jcoor) &
            & + factor*mat_trAB(DXk(jcoor),GbDX(icoor))
      enddo
      call mat_free(GbDX(icoor))
  enddo
  ENDIF

 ! #######################################################
 ! Debugging the X4:Tr (h^k)(D_0^B) Equivalent to X4 
 !  
 ! ######################################################
 ! Tr (h^k)(D_0)^B) 
 WRITE(lupri,*)  'icoor', 'jcoor', 'DEBUGGING:X4 |PURE (h^k)*(D0^b)| & 
 & DSX4A -2.0E0_realk*factor*mat_dotproduct(hk(jcoor),NDX(icoor))'
 do jcoor = 1,3*NATOMS
    call mat_init(hk(jcoor),nbast,nbast)
 enddo
 call II_get_prop(LUPRI,LUERR,molcfg%SETTING,hk,3*NATOMS,'PSO    ')   !Generation of h^k  
 do icoor=1,3     ! magnetic
    do jcoor=1,3*NATOMS  ! magnetic moment koordinate
       WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'DSX4A:',-2.0E0_realk*factor*mat_dotproduct(hk(jcoor),NDX(icoor))
       NucSpecNMSTtotal(icoor,jcoor) = NucSpecNMSTtotal(icoor,jcoor) &
            & - 2.0E0_realk*factor*mat_dotproduct(hk(jcoor),NDX(icoor))
    enddo
 enddo
 do icoor = 1,3*NATOMS
    call mat_free(hk(icoor))
 enddo
 do icoor = 1,3
    call mat_free(NDX(icoor))
 enddo

 !#############################################################################
 !     term: Tr(Dkb_2n+1 F)
 !
 !#############################################################################


 WRITE(lupri,*)  'icoor', 'jcoor', 'DEBUGGING: Tr(Dkb_2n+1 F) - 4.0E0_realk*factor*mat_dotproduct(DXk2np1A(jcoor,icoor),F(1))'
 do icoor=1,3     ! magnetic
    do jcoor=1,3*natoms  ! magnetic moment koordinate
       WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'DSX4B:',- 4.0E0_realk*factor*mat_dotproduct(DXk2np1A(jcoor,icoor),F(1))
       NucSpecNMSTtotal(icoor,jcoor) = NucSpecNMSTtotal(icoor,jcoor) &
            & - 4.0E0_realk*factor*mat_dotproduct(DXk2np1A(jcoor,icoor),F(1))
    enddo
 enddo
 WRITE(lupri,*)  'icoor', 'jcoor', 'DEBUGGING: Tr(Dkb_2n+1 F) - 4.0E0_realk*factor*mat_dotproduct(DXk2np1B(jcoor,icoor),F(1))'
 do icoor=1,3     ! magnetic
    do jcoor=1,3*natoms  ! magnetic moment koordinate
       WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'DSX4B:',- 4.0E0_realk*factor*mat_dotproduct(DXk2np1B(jcoor,icoor),F(1))
       NucSpecNMSTtotal(icoor,jcoor) = NucSpecNMSTtotal(icoor,jcoor) &
            & - 4.0E0_realk*factor*mat_dotproduct(DXk2np1B(jcoor,icoor),F(1))
    enddo
 enddo

 do icoor=1,3
   do jcoor = 1,3*NATOMS
    call mat_free(DXk2np1A(jcoor,icoor))
    call mat_free(DXk2np1B(jcoor,icoor))
   enddo
 enddo 
 do jcoor=1,3 
    call mat_free(DXk(jcoor))
 enddo
 deallocate(DXk2np1A)
 deallocate(DXk2np1B)
 ! start  printing numbers 
 allocate(atomname(natoms))
 do jcoor=1,natoms  
    atomname(jcoor)=molcfg%setting%molecule(1)%p%ATOM(jcoor)%Name
 enddo
 
 ! Prinitng for individual term X3 
 Write (lupri,*)   "Preexist term  - X3 :: NMST " 
 do jcoor=1,natoms  
    WRITE (LUPRI,'(20X,3(A,13X),/)') 'Bx', 'By', 'Bz'
    WRITE (LUPRI,'(2X,A8,3F15.8)')  atomname(jcoor)//CHRXYZ(1),(NucSpecNMSTtotal(K,3*jCOOR-2),K=1,3)
    WRITE (LUPRI,'(2X,A8,3F15.8)')  atomname(jcoor)//CHRXYZ(2),(NucSpecNMSTtotal(K,3*jCOOR-1),K=1,3)
    WRITE (LUPRI,'(2X,A8,3F15.8)')  atomname(jcoor)//CHRXYZ(3),(NucSpecNMSTtotal(K,3*jCOOR),K=1,3)
 enddo
  
 WRITE(LUPRI,*) "      Absolute chemical shift"
 WRITE(LUPRI,'(A,I12)') "number of atoms:",natoms
 do jcoor=1,natoms 
       WRITE (LUPRI,'(2X,A8,f15.8,A8)')  atomname(jcoor), &
        & 1.0E0_realk/3.0E0_realk*(NucSpecNMSTtotal(1,3*jCOOR-2)+NucSpecNMSTtotal(2,3*jCOOR-1) &
        & + NucSpecNMSTtotal(3,3*jCOOR))
 enddo
 deallocate(hk)
 deallocate(hb)
 deallocate(RHSk)
 deallocate(DXk)
 deallocate(atomname)
 WRITE(LUPRI,*) " Done with shielding tensor calculation"
 deallocate(NucSpecNMSTtotal)
! call mem_dealloc(AtomList)

end subroutine NMRshieldresponse_IANS_TK


!########################################################
! Subroutine for Nuclei Selected shielding. HF, and DFT  
! Based on Xk response 
!########################################################



subroutine NMRshieldresponse_NS(ls,molcfg,F,D,S)
  implicit none
  TYPE(lsitem),target     :: ls
  type(rsp_molcfg), intent(inout) :: molcfg
  type(Matrix),intent(in) :: F(1),D(1),S
  !
  logical :: Dsym
  real(realk),pointer          ::NMST(:,:),expval(:,:),NMST1(:,:),NMST2(:,:),NMST3(:,:),NMST4(:,:)
  real(realk),pointer          :: DSX4(:,:),DSX6(:,:),NucSpecNMSTtotal(:,:)
  real(realk)                  :: Factor,eival(1),eivalk(1)
  integer                      :: natoms,icoor,jcoor,k,lupri,nbast,luerr,n_rhs
  Character(len=4),allocatable :: atomName(:)
  real(realk)                  :: TS,TE
  type(Matrix)     :: Dx(3),NDX(3),Fx(3),Sx(3),tempm1,RHS(3),GbDs(3),Xx(1),Xk(1),GbJ(3), GbK(3),Gm(3),tf1,tf2, tempmx,GbDX(3)
  type(Matrix) :: GbXc(3)
  type(Matrix),pointer ::  RHSk(:),  DXk(:),hk(:),hb(:),GkD(:),tempk(:),sdf(:),DXD(:), DX0(:), DXk2np1A(:,:), DXk2np1B(:,:)
  type(matrix) :: Prod,Prod1
  integer              :: ntrial,nrhs,nsol,nomega,nstart,nAtomsSelected
  integer              :: istart,iend,iAtom,kcoor,offset
  integer,pointer :: AtomList(:)
  character(len=1)        :: CHRXYZ(-3:3)
  logical :: FoundNMRLabel
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


  !icoor is to be used for the 3 magnetic coordinates
  !jcoor is to be used for the 3*natoms magnetic moment coordinates

  !#############################################################################
  !##                                                                         ## 
  !##   Building D^k (Density matrix derivative w.r.t nuclei magnetic moment) ##
  !##     Calc RHS of eq 70 JCP 115 page 10349                                ##
  !##      RHS in this case;                                                  ##
  !##      ! RHSk = -P_anti(hkDS)                                             ##
  !##                                                                         ## 
  !#############################################################################
  print *, 'test 1' 
  ! Generation of hkDS      
  allocate(hk(3*nAtomsSelected))     ! Matrix h^k: Generation by 'PSO    ' 
  allocate(RHSk(3*nAtomsSelected))   ! RHSk - R.H.S. to solve eq. 70 for D^K 
  allocate(Dxk(3*nAtomsSelected))    ! (D^K Derivative of Density matrix) 
  
  do jcoor = 1,3*nAtomsSelected
     call mat_init(hk(jcoor),nbast,nbast)
  enddo
  print *, 'test 2'
  do jcoor=1,nAtomsSelected
    istart=3*jcoor-2
    iend=istart+2
    iAtom = AtomList(jcoor)
    call II_get_PSO_spec(LUPRI,LUERR,molcfg%SETTING,hk(istart:iend),iAtom)
    !call II_get_prop(LUPRI,LUERR,molcfg%SETTING,hk,3*NATOMS,'PSO    ')   !Generation of h^k  
  enddo

  print *, 'test 3'
  WRITE(lupri,*)'hk icoor=1 '
  call mat_print(hk(1),1,hk(1)%nrow,1,hk(1)%ncol,lupri)
  
  call mat_init(tempm1,nbast,nbast) 
  call mat_mul(D(1),S,'n','n',1.0E0_realk,0.0E0_realk,tempm1)
  do jcoor = 1,3*nAtomsSelected
     call mat_init(RHSk(jcoor),nbast,nbast)
     !Generation of RHSk = h^kDS  
     call mat_mul(hk(jcoor),tempm1,'n','n',1.0E0_realk,0.0E0_realk,RHSk(jcoor))  
     call mat_mul(tempm1,hk(jcoor),'t','n',1.0E0_realk,-1.0E0_realk,RHSk(jcoor))  
     call mat_free(hk(jcoor))
  enddo  
  call mat_free(tempm1)
  print *, 'test 4'
  !###########################################################################
  !Solve K([D,Xk]s) = RHSk is now ready                                
  !Solve to get Dk=[D,Xk]_S + D0^k --> Xk (Eq. 70)                    
  !We can use the RSP solver, since K([D,Xk]s) = sigma = -1/2 * E[2]Xk
  !
  !############################################################################
     
  call mat_init(Xk(1),nbast,nbast)                           
  allocate(DXk2np1A(3*nAtomsSelected,3))
  allocate(DXk2np1B(3*nAtomsSelected,3))
  do jcoor=1,3*nAtomsSelected
     eivalk(1)=0.0E0_realk
     write(lupri,*)'Calling rsp solver for Xk  '
     call util_scriptPx('T',D(1),S,RHSk(jcoor))
     if ( mat_dotproduct(RHSk(jcoor),RHSk(jcoor))>1.0d-10) then

        ntrial = 1 !# of trial vectors in a given iteration (number of RHS)
        nrhs = 1   !# of RHS only relevant for linear equations (lineq_x = TRUE)
        nsol = 1   !# of solution (output) vectors
        nomega = 1 !If lineq_x, number of laser freqs (input)
                   !Otherwise number of excitation energies (output) 
        nstart = 1 !Number of start vectors. Only relevant for eigenvalue problem
        !ntrial and nstart seem to be obsolete 
        call rsp_init(ntrial,nrhs,nsol,nomega,nstart)
        call rsp_solver(molcfg,D(1),S,F(1),.true.,nrhs,RHSk(jcoor:jcoor),EIVALK,Xk)
     else
        write(lupri,*) 'WARNING: RHSk norm is less than threshold'
        write(lupri,*) 'LIN RSP equations NOT solved for this RHS    '
        call mat_zero(Xk(1))
     end if
     call mat_free(RHSk(jcoor))

     call mat_init(Prod,nbast,nbast)
     call mat_init(Prod1,nbast,nbast)
     call mat_init(SX(1),nbast,nbast)
     call mat_init(SX(2),nbast,nbast)
     call mat_init(SX(3),nbast,nbast)                                     
     call II_get_magderivOverlap(Sx,molcfg%setting,lupri,luerr)
     call mat_init(tempm1,nbast,nbast)
     do icoor = 1,3 
        call mat_mul(SX(icoor),D(1),'n','n',1E0_realk,0.E0_realk,tempm1)     !tempm1 = SX*D
        call mat_init(NDX(icoor),nbast,nbast)        
        call mat_mul(D(1),tempm1,'n','n',-1.0E0_realk,0.E0_realk,NDX(icoor))  !NDX = -D*SX*D
        
        ![D0b,Xk]_S
        call mat_init(DXk2np1A(jcoor,icoor),nbast,nbast)
        call mat_init(DXk2np1B(jcoor,icoor),nbast,nbast)
        call mat_mul(NDX(icoor),S,'n','n',1E0_realk,0E0_realk,Prod)    
        call mat_mul(Prod,Xk(1),'n','n',1E0_realk,0E0_realk,DXk2np1A(jcoor,icoor))  
        call mat_mul(S,NDX(icoor),'n','n',1E0_realk,0E0_realk,Prod)    
        call mat_mul(Xk(1),Prod,'n','n',1E0_realk,0E0_realk,Prod1) 
        call mat_daxpy(-1E0_realk,Prod1,DXk2np1A(jcoor,icoor)) 
        call mat_free(NDX(icoor))
        
        !+[D0,Xk]_Sb
        call mat_mul(D(1),SX(icoor),'n','n',1E0_realk,0E0_realk,Prod)    
        call mat_mul(Prod,Xk(1),'n','n',1E0_realk,0E0_realk,DXk2np1B(jcoor,icoor))  
        call mat_mul(SX(icoor),D(1),'n','n',1E0_realk,0E0_realk,Prod)    
        call mat_mul(Xk(1),Prod,'n','n',1E0_realk,0E0_realk,Prod1) 
        call mat_daxpy(-1E0_realk,Prod1,DXk2np1B(jcoor,icoor)) 
        call mat_free(SX(icoor))                    
     enddo
     call mat_free(tempm1)
     call mat_free(Prod)
     call mat_free(Prod1)
     
     !Make D^k=[D_0,X^k]_s
     call mat_init(DXk(jcoor),nbast,nbast)
     call ABCcommutator(nbast,D(1),Xk(1),S,DXk(jcoor))       !Generate [D_0,X^k]_s 
     call mat_scal(-4.0d0,DXk(jcoor))
  enddo
  call mat_free(Xk(1)) 
  print *, 'test 5'
  write(lupri,*) '---------------------------------------------------------'
  write(lupri,*) 'NUCLEA MAGNETIC TENSOR (INASHIELD)                                  '
  write(lupri,*) '---------------------------------------------------------'

  call mem_alloc(NucSpecNMSTtotal,3,3*nAtomsSelected)

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
     ! Printing  (D^*(h^kb) in output file  over (icoor,jcoor) 
     WRITE(lupri,*)  'icoor', 'jcoor', 'X1, |PURE (D)*(h^kb)|  :  2*factor*expval(jcoor,icoor)'  
     do icoor=1,3     ! magnetic koordinate
        do kcoor=1,3  ! magnetic moment koordinate
           WRITE(lupri,*) 'icoor', icoor, 'jcoor', kcoor, 'expval1:',   2*factor*expval(icoor,kcoor)
           NucSpecNMSTtotal(icoor,offset+kcoor) = 2*factor*expval(icoor,kcoor)
        enddo
     enddo
  enddo
  call mem_dealloc(expval)
  print *, 'test 6'
  !###########################################################################
  !    Tr D^k*h^b
  !
  !###########################################################################
  
  allocate(hb(3))   ! Allocation to matrix h^b 
  do icoor = 1,3
     call mat_init(hb(icoor),nbast,nbast)
  enddo
   
  call II_get_prop(LUPRI,LUERR,molcfg%SETTING,hb,3,'MAGMOM ')   ! Generate h^b  

  ! Printing  (D^k)*(h^b) in output file  over (icoor,jcoor) 
  WRITE(lupri,*)  'icoor', 'jcoor', 'X1, |PURE (D^k)*(h^b)|  :  - factor*mat_trAB(DXk(jcoor),hb(icoor))'  
  do icoor=1,3
     do jcoor=1,3*nAtomsSelected ! magnetic moment coordinate
       WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'expval2:',+ factor*mat_trAB(DXk(jcoor),hb(icoor))
       NucSpecNMSTtotal(icoor,jcoor) = NucSpecNMSTtotal(icoor,jcoor) &
            & + factor*mat_trAB(DXk(jcoor),hb(icoor))
     enddo
     call mat_free(hb(icoor))
  enddo
  print *, 'test 7'
  !###########################################################################
  !   Tr (D^k)*(G^b(D)) 
  !
  !###########################################################################
  do icoor=1,3 
     call mat_init(GbJ(icoor),nbast,nbast)
     call mat_init(GbK(icoor),nbast,nbast)
     call mat_zero(GbK(icoor))
     call mat_init(Gm(icoor),nbast,nbast) 
  enddo
  call II_get_magderivJ(LUPRI,LUERR,molcfg%SETTING,nbast,D,GbJ)   !Generate J^b
  call II_get_magderivK(LUPRI,LUERR,molcfg%SETTING,nbast,D,GbK)   !Generate K^b 
  ! Printing  (D^k)*(G^b(D)) in output file  over (icoor,jcoor) 
  WRITE(lupri,*)  'icoor', 'jcoor', ' X5 |PURE (D^K)*(J^b(D))|  '
  do icoor=1,3
     do jcoor=1,3*nAtomsSelected
        WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'NMST2J:', factor*mat_trAB(DXk(jcoor),GbJ(icoor)) 
     enddo
  enddo
  WRITE(lupri,*)  'icoor', 'jcoor', ' X5 |PURE (D^K)*(K^b(D))|  '
  do icoor=1,3
     do jcoor=1,3*nAtomsSelected
        WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'NMST2K:', factor*mat_trAB(DXk(jcoor),GbK(icoor)) 
        !Generate G^b= J^b+K^b
        call mat_add(1.0E0_realk,GbJ(icoor),1.0E0_realk,GbK(icoor),Gm(icoor)) 
        NucSpecNMSTtotal(icoor,jcoor) = NucSpecNMSTtotal(icoor,jcoor) &
             & + factor*mat_trAB(DXk(jcoor),Gm(icoor)) 
     enddo
     call mat_free(GbJ(icoor))
     call mat_free(GbK(icoor))
     call mat_free(Gm(icoor))
  enddo

  if(molcfg%setting%do_dft)then  
     do icoor=1,3 
        call mat_init(GbXc(icoor),nbast,nbast)
     enddo
     call II_get_xc_magderiv_kohnsham_mat(LUPRI,LUERR,molcfg%SETTING,nbast,D(1),GbXc)
     ! Printing  (D^k)*(G^b(D)) in output file  over (icoor,jcoor) 
     WRITE(lupri,*)  'icoor', 'jcoor', ' X5 |PURE (D^K)*(G^bxc(D))|  '
     do icoor=1,3
        do jcoor=1,3*nAtomsSelected
           WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'NMST2K:', factor*mat_trAB(DXk(jcoor),GbXc(icoor)) 
           
           NucSpecNMSTtotal(icoor,jcoor) = NucSpecNMSTtotal(icoor,jcoor) &
                & + factor*mat_trAB(DXk(jcoor),GbXc(icoor)) 
        enddo
        call mat_free(GbXc(icoor))
     enddo
  endif
  print *, 'test 8'
  !#################################################################################
  !      Generate D0X = -D*SX*D
  !
  !#################################################################################
  call mat_init(SX(1),nbast,nbast)
  call mat_init(SX(2),nbast,nbast)
  call mat_init(SX(3),nbast,nbast)
  call II_get_magderivOverlap(Sx,molcfg%setting,lupri,luerr)
  call mat_init(tempm1,nbast,nbast)
  do icoor = 1,3 
     call mat_mul(SX(icoor),D(1),'n','n',1E0_realk,0.E0_realk,tempm1)     !tempm1 = SX*D
     call mat_free(SX(icoor))                    
     call mat_init(NDX(icoor),nbast,nbast)        
     call mat_mul(D(1),tempm1,'n','n',-1.0E0_realk,0.E0_realk,NDX(icoor))  !NDX = -D*SX*D
  enddo
  call mat_free(tempm1)

  !#################################################################################
  !     Tr (D^k)*(G(D0^b)) 
  !
  !#################################################################################

  call mat_init(GbDX(1),nbast,nbast)
  call mat_init(GbDX(2),nbast,nbast)  
  call mat_init(GbDX(3),nbast,nbast)                                  
  call di_GET_GbDs(lupri,luerr,NDX,GbDX,3,molcfg%setting)  ! G(D0^B)
  WRITE(lupri,*)  'icoor', 'jcoor', 'DEBUGGING: X6 |PURE (D^k)*(G(D^b))| &
  & DSX6:  - factor*mat_trAB(DXk(jcoor),GbDX(icoor))'
  do icoor=1,3
    do jcoor=1,3*nAtomsSelected
       WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'DSX6:',  factor*mat_trAB(DXk(jcoor),GbDX(icoor))
       NucSpecNMSTtotal(icoor,jcoor) = NucSpecNMSTtotal(icoor,jcoor) &
            & + factor*mat_trAB(DXk(jcoor),GbDX(icoor))
    enddo
    call mat_free(GbDX(icoor))
 enddo

 IF(molcfg%setting%do_dft)THEN
 call mat_init(GbDX(1),nbast,nbast)
 call mat_init(GbDX(2),nbast,nbast)  
 call mat_init(GbDX(3),nbast,nbast)                                  
 call II_get_xc_linrsp(lupri,luerr,molcfg%setting,nbast,NDX,D(1),GbDX,3) 
     
 WRITE(lupri,*)  'icoor', 'jcoor', 'DEBUGGING: X6- XC inculded &
  & |PURE (D^k)*(G(D^b))|  DSX6:  - factor*mat_trAB(DXk(jcoor),GbDX(icoor))'
  do icoor=1,3
     do jcoor=1,3*nAtomsSelected
       WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'DSX6:',  factor*mat_trAB(DXk(jcoor),GbDX(icoor))
       NucSpecNMSTtotal(icoor,jcoor) = NucSpecNMSTtotal(icoor,jcoor) &
            & + factor*mat_trAB(DXk(jcoor),GbDX(icoor))
      enddo
      call mat_free(GbDX(icoor))
  enddo
  ENDIF

 print *, 'test 9'
 ! #######################################################
 ! Tr (h^k)(D_0^B) 
 !  
 ! ######################################################
 ! Tr (h^k)(D_0)^B) 
 WRITE(lupri,*)  'icoor', 'jcoor', 'DEBUGGING:X4 |PURE (h^k)*(D0^b)| & 
 & DSX4A -2.0E0_realk*factor*mat_dotproduct(hk(jcoor),NDX(icoor))'
 do jcoor = 1,3*nAtomsSelected
    call mat_init(hk(jcoor),nbast,nbast)
 enddo
 print *, 'test 2'
 do jcoor=1,nAtomsSelected
    istart=3*jcoor-2
    iend=istart+2
    iAtom = AtomList(jcoor)
    call II_get_PSO_spec(LUPRI,LUERR,molcfg%SETTING,hk(istart:iend),iAtom)
    !call II_get_prop(LUPRI,LUERR,molcfg%SETTING,hk,3*NATOMS,'PSO    ')   !Generation of h^k  
 enddo 
 do icoor=1,3     ! magnetic
    do jcoor=1,3*NATOMSselected ! magnetic moment koordinate
       WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'DSX4A:',-2.0E0_realk*factor*mat_dotproduct(hk(jcoor),NDX(icoor))
       NucSpecNMSTtotal(icoor,jcoor) = NucSpecNMSTtotal(icoor,jcoor) &
            & - 2.0E0_realk*factor*mat_dotproduct(hk(jcoor),NDX(icoor))
    enddo
 enddo
 do icoor = 1,3*nAtomsSelected
    call mat_free(hk(icoor))
 enddo

 print *, 'test 10'
 !#############################################################################
 !     term: Tr(Dkb_2n+1 F)
 !
 !#############################################################################


 WRITE(lupri,*)  'icoor', 'jcoor', 'DEBUGGING: Tr(Dkb_2n+1 F) - 4.0E0_realk*factor*mat_dotproduct(DXk2np1A(jcoor,icoor),F(1))'
 do icoor=1,3     ! magnetic
    do jcoor=1,3*nAtomsSelected  ! magnetic moment koordinate
       WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'DSX4B:',- 4.0E0_realk*factor*mat_dotproduct(DXk2np1A(jcoor,icoor),F(1))
       NucSpecNMSTtotal(icoor,jcoor) = NucSpecNMSTtotal(icoor,jcoor) &
            & - 4.0E0_realk*factor*mat_dotproduct(DXk2np1A(jcoor,icoor),F(1))
    enddo
 enddo
 WRITE(lupri,*)  'icoor', 'jcoor', 'DEBUGGING: Tr(Dkb_2n+1 F) - 4.0E0_realk*factor*mat_dotproduct(DXk2np1B(jcoor,icoor),F(1))'
 do icoor=1,3     ! magnetic
    do jcoor=1,3*nAtomsSelected  ! magnetic moment koordinate
       WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'DSX4B:',- 4.0E0_realk*factor*mat_dotproduct(DXk2np1B(jcoor,icoor),F(1))
       NucSpecNMSTtotal(icoor,jcoor) = NucSpecNMSTtotal(icoor,jcoor) &
            & - 4.0E0_realk*factor*mat_dotproduct(DXk2np1B(jcoor,icoor),F(1))
    enddo
 enddo
 print *, 'test 11'
 ! start  printing numbers 
 allocate(atomname(natoms))
 do jcoor=1,natomsselected  
    iAtom = AtomList(jcoor)
    atomname(jcoor)=molcfg%setting%molecule(1)%p%ATOM(iAtom)%Name
 enddo
 
 ! Prinitng for individual term X3 
 Write (lupri,*)   "Preexist term  - X3 :: NMST " 
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
 deallocate(hk)
 deallocate(hb)
 deallocate(RHSk)
 deallocate(DXk)
 deallocate(atomname)
 WRITE(LUPRI,*) " Done with shielding tensor calculation"

 call mem_dealloc(AtomList)

end subroutine NMRshieldresponse_NS


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
  type(Matrix)     :: NDX(3),Sx(3),tempm1,GbJ(3), GbK(3),GbDX(3)
  type(Matrix) :: GbXc(3), GbDXc(3)
  type(Matrix),pointer ::  RHSk(:), DXk(:),hk(:),hb(:),Xk(:)
  type(Matrix),pointer :: ProdA(:)
  type(matrix) :: DS,SD
  integer              :: ntrial,nrhs,nsol,nomega,nstart,nAtomsSelected
  integer              :: istart,iend,iAtom,kcoor,offset
  integer,pointer :: AtomList(:)
  character(len=1)        :: CHRXYZ(-3:3)
  logical :: FoundNMRLabel
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


  !icoor is to be used for the 3 magnetic coordinates
  !jcoor is to be used for the 3*natoms magnetic moment coordinates

  !#############################################################################
  !##                                                                         ## 
  !##   Building D^k (Density matrix derivative w.r.t nuclei magnetic moment) ##
  !##     Calc RHS of eq 70 JCP 115 page 10349                                ##
  !##      RHS in this case;                                                  ##
  !##      ! RHSk = -P_anti(hkDS)                                             ##
  !##                                                                         ## 
  !#############################################################################
  ! Generation of hkDS      
  allocate(hk(3*nAtomsSelected))     ! Matrix h^k: Generation by 'PSO    ' 
  allocate(RHSk(3*nAtomsSelected))   ! RHSk - R.H.S. to solve eq. 70 for D^K 
  allocate(Dxk(3*nAtomsSelected))    ! (D^K Derivative of Density matrix) 
  
  do jcoor = 1,3*nAtomsSelected
     call mat_init(hk(jcoor),nbast,nbast)
  enddo
  call mat_init(tempm1,nbast,nbast) 
  call mat_mul(D(1),S,'n','n',1.0E0_realk,0.0E0_realk,tempm1)
  do jcoor=1,nAtomsSelected
     istart=3*jcoor-2
     iend=istart+2
     iAtom = AtomList(jcoor)
     call II_get_PSO_spec(LUPRI,LUERR,molcfg%SETTING,hk(istart:iend),iAtom)
     do icoor=istart,iend
        call mat_init(RHSk(icoor),nbast,nbast)
        call mat_mul(hk(icoor),tempm1,'n','n',1.0E0_realk,0.0E0_realk,RHSk(icoor))  
        call mat_mul(tempm1,hk(icoor),'t','n',1.0E0_realk,-1.0E0_realk,RHSk(icoor))  
        call mat_free(hk(icoor))
     enddo
  enddo
  call mat_free(tempm1)
  !###########################################################################
  !Solve K([D,Xk]s) = RHSk is now ready                                
  !Solve to get Dk=[D,Xk]_S + D0^k --> Xk (Eq. 70)                    
  !We can use the RSP solver, since K([D,Xk]s) = sigma = -1/2 * E[2]Xk
  !
  !############################################################################
     
  call mem_alloc(ProdA,3)
  do icoor=1,3      
     call mat_init(ProdA(icoor),nbast,nbast)
  enddo
  call mem_alloc(NucSpecNMSTtotal,3,3*nAtomsSelected)
  call mem_alloc(Prodtotal,3,3*nAtomsSelected)
  allocate(hb(3))   ! Allocation to matrix h^b 
  call mat_init(GbDX(1),nbast,nbast)
  call mat_init(GbDX(2),nbast,nbast)  
  call mat_init(GbDX(3),nbast,nbast)                                  
  call mat_init(GbDXc(1),nbast,nbast)
  call mat_init(GbDXc(2),nbast,nbast)  
  call mat_init(GbDXc(3),nbast,nbast)                                  
!  call mat_init(Prod,nbast,nbast)
  call mat_init(DS,nbast,nbast)
  call mat_init(SD,nbast,nbast)
  call mat_init(SX(1),nbast,nbast)
  call mat_init(SX(2),nbast,nbast)
  call mat_init(SX(3),nbast,nbast)                                     
  do icoor=1,3
     call mat_init(GbJ(icoor),nbast,nbast)
     call mat_init(GbK(icoor),nbast,nbast)
     call mat_zero(GbK(icoor))
     call mat_init(hb(icoor),nbast,nbast)
  enddo 
  call II_get_magderivOverlap(Sx,molcfg%setting,lupri,luerr)
  call II_get_prop(LUPRI,LUERR,molcfg%SETTING,hb,3,'MAGMOM ')   ! Generate h^b  
  call II_get_magderivJ(LUPRI,LUERR,molcfg%SETTING,nbast,D,GbJ)   !Generate J^b
  call II_get_magderivK(LUPRI,LUERR,molcfg%SETTING,nbast,D,GbK)   !Generate K^b 
  call mat_init(tempm1,nbast,nbast)
  do icoor = 1,3 
     call mat_mul(SX(icoor),D(1),'n','n',1E0_realk,0.E0_realk,tempm1)     !tempm1 = SX*D
     call mat_init(NDX(icoor),nbast,nbast)        
     call mat_mul(D(1),tempm1,'n','n',-1.0E0_realk,0.E0_realk,NDX(icoor))  !NDX = -D*SX*D

     !Contribution To RHS: -FD0SB 
     call mat_mul(D(1),SX(icoor), 'n','n',1E0_realk,0E0_realk,tempm1) 
     call mat_mul(F(1),tempm1,'n','n',1E0_realk,0E0_realk,ProdA(icoor)) 
     !Contribution To RHS: SBD0F
     call mat_mul(D(1),F(1), 'n','n',1E0_realk,0E0_realk,tempm1) 
     call mat_mul(SX(icoor),tempm1,'n','n',-1E0_realk,1E0_realk,ProdA(icoor)) 
  enddo
  call mat_free(tempm1)
  call di_GET_GbDs(lupri,luerr,NDX,GbDX,3,molcfg%setting)  ! G(D0^B)
  IF(molcfg%setting%do_dft)THEN
     call II_get_xc_linrsp(lupri,luerr,molcfg%setting,nbast,NDX,D(1),GbDXc,3) 
  endif          
  do icoor=1,3 
     call mat_init(GbXc(icoor),nbast,nbast)
  enddo
  if(molcfg%setting%do_dft)then  
     call II_get_xc_magderiv_kohnsham_mat(LUPRI,LUERR,molcfg%SETTING,nbast,D(1),GbXc)
  endif
  call mat_mul(D(1),S, 'n','n',1E0_realk,0E0_realk,DS) 
  call mat_mul(S,D(1), 'n','n',1E0_realk,0E0_realk,SD) 
   
  call mat_init(tempm1,nbast,nbast) 
  do icoor=1,3      
     !Contribution To RHS: FD0^bS 
     call mat_mul(NDX(icoor),S,'n','n',1E0_realk,0E0_realk,tempm1) 
     call mat_mul(F(1),tempm1,'n','n',-1E0_realk,1E0_realk,ProdA(icoor))
     !Contribution To RHS: -SDOBF
     call mat_mul(NDX(icoor),F(1),'n','n',1E0_realk,0E0_realk,tempm1) 
     call mat_mul(S,tempm1,'n','n',1E0_realk,1E0_realk,ProdA(icoor))
     !Contribution To RHS: -FD0SB 
!     call mat_mul(D(1),SX(icoor), 'n','n',1E0_realk,0E0_realk,tempm1) 
!     call mat_mul(F(1),tempm1,'n','n',1E0_realk,1E0_realk,ProdA(icoor)) 
     !Contribution To RHS:  SBD0F
!     call mat_mul(D(1),F(1), 'n','n',1E0_realk,0E0_realk,tempm1) 
!     call mat_mul(SX(icoor),tempm1,'n','n',-1E0_realk,1E0_realk,ProdA(icoor)) 
     !Contribution To RHS:  -hbDS 
     call mat_mul(hb(icoor),DS,'n','n',1E0_realk,1E0_realk,ProdA(icoor)) 
     !Contribution To RHS:  SDhb
     call mat_mul(SD,hb(icoor),'n','n',-1E0_realk,1E0_realk,ProdA(icoor)) 
     !Contribution To RHS: -J^b(D)D0S 
     
     call mat_mul(GbJ(icoor),DS,'n','n',1E0_realk,1E0_realk,ProdA(icoor)) 
     !Contribution To RHS:   - SD0J^b(D)
     call mat_mul(SD,GbJ(icoor),'n','n',-1E0_realk,1E0_realk,ProdA(icoor)) 
     
     !Contribution To RHS: -K^b(D)D0S  
     call mat_mul(GbK(icoor),DS,'n','n',1E0_realk,1E0_realk,ProdA(icoor))
     !Contribution To RHS: SD0K^b(D) 
     call mat_mul(SD,Gbk(icoor),'n','n',-1E0_realk,1E0_realk,ProdA(icoor)) 
     !Contribution To RHS: G(D0B)D0S
     call mat_mul(GbDX(icoor),DS,'n','n',1E0_realk,1E0_realk,ProdA(icoor)) 
     !Contribution To RHS: -SD0G(D0B)
     call mat_mul(SD,GbDX(icoor),'n','n',-1E0_realk,1E0_realk,ProdA(icoor)) 
     !Contribution To RHS: G^bDOS (DFT XC) - SD0G^b (DFT XC)
     if(molcfg%setting%do_dft)then  
        call mat_mul(GbXc(icoor),DS,'n','n',1E0_realk,1E0_realk,ProdA(icoor)) 
        call mat_mul(SD,GbXc(icoor),'n','n',-1E0_realk,1E0_realk,ProdA(icoor)) 
     endif
     !??????????????????????????????????????????+
     IF(molcfg%setting%do_dft)THEN
        call mat_mul(GbDXc(icoor),DS,'n','n',1E0_realk,1E0_realk,ProdA(icoor)) 
        call mat_mul(SD,GbDxc(icoor),'n','n',-1E0_realk,1E0_realk,ProdA(icoor)) 
     endif
  enddo
  call mat_free(tempm1)
  call mat_free(DS)
  call mat_free(SD)
  do icoor=1,3
     call mat_free(SX(icoor))
     call mat_free(GbJ(icoor))
     call mat_free(GbK(icoor))
     call mat_free(GbDX(icoor))
     call mat_free(GbXc(icoor))
     call mat_free(GbDXc(icoor))                                  
     call mat_free(hb(icoor))
  enddo 
  IF(ls%input%DALTON%SolveNMRResponseSimultan)THEN
     call mem_alloc(eivalkF,3*nAtomsSelected)
     eivalkF=0.0E0_realk
     write(lupri,*)'Calling rsp solver for all Xk  '
     allocate(Xk(3*nAtomsSelected))   
     do jcoor=1,3*nAtomsSelected
        call mat_init(Xk(jcoor),nbast,nbast)                           
        call util_scriptPx('T',D(1),S,RHSk(jcoor))     
        if ( mat_dotproduct(RHSk(jcoor),RHSk(jcoor)).LT.1.0d-10) then
           print*,'WARNING RHS(jcoor=',jcoor,') = Zero '
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
     call rsp_solver(molcfg,D(1),S,F(1),.true.,nrhs,RHSk,EIVALKF,Xk)
     do jcoor=1,3*nAtomsSelected
        do icoor = 1,3 
           Prodtotal(icoor,jcoor)=-4E0_realk*factor*mat_trAB(Xk(jcoor),ProdA(icoor))
        enddo
     enddo
     print*,'FINAL -4E0_realk*factor*mat_trAB(Xk(jcoor),ProdA(icoor))'
     call ls_output(Prodtotal,1,3,1,3*nAtomsSelected,3,3*nAtomsSelected,1,6)

     do jcoor=1,3*nAtomsSelected
        call mat_free(RHSk(jcoor))
        call mat_free(Xk(jcoor))
     enddo
     deallocate(Xk)
     call mem_dealloc(eivalkF)
  ELSE
     allocate(Xk(1))   
     call mat_init(Xk(1),nbast,nbast)                                
     call mem_alloc(eivalkF,1)
     eivalkF=0.0E0_realk
     do jcoor=1,3*nAtomsSelected
        write(lupri,*)'Calling rsp solver for Xk  '
        call util_scriptPx('T',D(1),S,RHSk(jcoor))
        if ( mat_dotproduct(RHSk(jcoor),RHSk(jcoor))>1.0d-10) then
           ntrial = 1 !# of trial vectors in a given iteration (number of RHS)
           nrhs = 1   !# of RHS only relevant for linear equations (lineq_x = TRUE)
           nsol = 1   !# of solution (output) vectors
           nomega = 1 !If lineq_x, number of laser freqs (input)
           !Otherwise number of excitation energies (output) 
           nstart = 1 !Number of start vectors. Only relevant for eigenvalue problem
           !ntrial and nstart seem to be obsolete 
           call rsp_init(ntrial,nrhs,nsol,nomega,nstart)
           call rsp_solver(molcfg,D(1),S,F(1),.true.,nrhs,RHSk(jcoor:jcoor),EIVALKF,Xk)
        else
           write(lupri,*) 'WARNING: RHSk norm is less than threshold'
           write(lupri,*) 'LIN RSP equations NOT solved for this RHS    '
           call mat_zero(Xk(1))
        end if
        call mat_free(RHSk(jcoor))     
        do icoor = 1,3 
           Prodtotal(icoor,jcoor)=-4E0_realk*factor*mat_trAB(Xk(1),ProdA(icoor))
        enddo
     enddo
     call mat_free(Xk(1)) 
     deallocate(Xk)           
     call mem_dealloc(eivalkF)
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
     ! Printing  (D^*(h^kb) in output file  over (icoor,jcoor) 
!     WRITE(lupri,*)  'icoor', 'jcoor', 'X1, |PURE (D)*(h^kb)|  :  2*factor*expval(jcoor,icoor)'  
     do icoor=1,3     ! magnetic koordinate
        do kcoor=1,3  ! magnetic moment koordinate
!          WRITE(lupri,*) 'icoor', icoor, 'jcoor', kcoor, 'expval1:',   2*factor*expval(icoor,kcoor)
           NucSpecNMSTtotal(icoor,offset+kcoor) = 2*factor*expval(icoor,kcoor)
        enddo
     enddo
  enddo
  call mem_dealloc(expval)
 ! print *, 'test latest 1'
  !#################################################################################
  !(h^k)*(D0^b)
  !
  !#################################################################################
     
  do jcoor = 1,3*nAtomsSelected
     call mat_init(hk(jcoor),nbast,nbast)
  enddo
  !WRITE(lupri,*)  'icoor', 'jcoor', ' (h^k)*(D0^b)| & 
  !& DSX4A -2.0E0_realk*factor*mat_dotproduct(hk(jcoor),NDX(icoor))'
  do jcoor=1,nAtomsSelected
    istart=3*jcoor-2
    iend=istart+2
    iAtom = AtomList(jcoor)
    call II_get_PSO_spec(LUPRI,LUERR,molcfg%SETTING,hk(istart:iend),iAtom)
    !call II_get_prop(LUPRI,LUERR,molcfg%SETTING,hk,3*NATOMS,'PSO    ')   !Generation of h^k  
  enddo 
  do icoor=1,3     ! magnetic
    do jcoor=1,3*NATOMSselected ! magnetic moment koordinate
!      WRITE(lupri,*) 'icoor', icoor, 'jcoor', jcoor, 'DSX4A:',-2.0E0_realk*factor*mat_dotproduct(hk(jcoor),NDX(icoor))
       NucSpecNMSTtotal(icoor,jcoor) = NucSpecNMSTtotal(icoor,jcoor) &
            & - 2.0E0_realk*factor*mat_dotproduct(hk(jcoor),NDX(icoor))
    enddo
  enddo
  do icoor = 1,3*nAtomsSelected
    call mat_free(hk(icoor))
  enddo
  do icoor=1,3 
    call mat_free(NDX(icoor))
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


  allocate(atomname(natoms))
  do jcoor=1,natomsselected  
    iAtom = AtomList(jcoor)
    atomname(jcoor)=molcfg%setting%molecule(1)%p%ATOM(iAtom)%Name
  enddo
 
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
  call mem_dealloc(ProdA)
  call mem_dealloc(Prodtotal)
  deallocate(hk)
  deallocate(hb)
  deallocate(RHSk)
  deallocate(DXk)
  deallocate(atomname)
  call mem_dealloc(AtomList)
  WRITE(LUPRI,*) " Done with shielding tensor calculation"

 end subroutine NMRshieldresponse_RSTNS
end module nuclei_selected_shielding_mod 
