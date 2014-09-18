!> @fileddd
!> Contains property subroutines
!> \author Thomas Kjaergaard
!> \date 2014
module response_noOpenRSP_module
  use precision
  use configurationType, only: configitem
  use TYPEDEFTYPE,   only: LSSETTING,lsitem
  use memory_handling, only: mem_alloc, mem_dealloc
  use files, only: lsopen, lsclose
  use lstiming
  use matrix_module, only: matrix
  use matrix_operations
  use matrix_util
  use rspsolver
  use decompMod, only: decompItem
  use RSPsolver
  use IntegralInterfaceMOD
  use II_XC_interfaceModule
  use dal_interface
  use rsp_util
  public NMRshieldresponse_noOpenRSP, lsdalton_response_noOpenRSP

  private

Contains
subroutine lsdalton_response_noOpenRSP(ls,config,F,D,S)
  implicit none
  TYPE(lsitem),target     :: ls
  type(configItem),target :: config
  TYPE(Matrix)            :: F(1),D(1),S
  type(rsp_molcfg)        :: molcfg
  real(realk)             :: Tstart,Tend,t1,t2 !,ten,tstr,E,gradnrm

  call CPU_TIME(tstart)
  call LSTIMER('START',t1,t2,config%LUPRI)
  IF(config%response%tasks%doResponse)THEN
     if(config%response%RSPSOLVERinput%rsp_mo_precond) then
        call util_save_MOinfo(F(1),S,config%decomp%nocc) 
     endif
     
     call init_rsp_molcfg(molcfg,S,ls%setting%MOLECULE(1)%p%Natoms, &
          & config%decomp%lupri,config%decomp%luerr, &
          & ls%setting,config%decomp,config%response%rspsolverinput)     
     
     if(config%response%tasks%doNMRshield) then
        call NMRshieldresponse_noOpenRSP(molcfg,F,D,S)
     endif

     if(config%response%RSPSOLVERinput%rsp_mo_precond)then
        call util_free_MOstuff()
     endif
  ENDIF
  call LSTIMER('LSDALTON RSP',t1,t2,config%LUPRI,.TRUE.)
  call CPU_TIME(tend)
  WRITE(config%lupri,*) "*****************************************************"
  Write(config%lupri,*) "**     CPU-TIME USED IN LSDALTON RESPONSE: ",tend-tstart,"   **"
  WRITE(config%lupri,*) "*****************************************************"

end subroutine lsdalton_response_noOpenRSP

subroutine NMRshieldresponse_noOpenRSP(molcfg,F,D,S)
  implicit none
  type(rsp_molcfg), intent(inout) :: molcfg
  type(Matrix),intent(in) :: F(1),D(1),S
  !
  real(realk),pointer          :: NMST(:,:),expval(:,:)
  real(realk)                  :: Factor,eival(1)
  integer                      :: natoms,icoor,jcoor,k,lupri,nbast,luerr,n_rhs
  Character(len=4),allocatable :: atomName(:)
  real(realk)                  :: TS,TE
  type(Matrix)                 :: Dx(3),Fx(3),Sx(3),tempm1,RHS(3),GbDs(3),Xx(1)
  character(len=1)        :: CHRXYZ(-3:3)
  DATA CHRXYZ /'z','y','x',' ','X','Y','Z'/

  nbast = D(1)%nrow
  lupri = molcfg%lupri
  luerr = lupri
  natoms = molcfg%natoms
  CALL LSTIMER('START ',TS,TE,LUPRI)

  WRITE(LUPRI,*) '=========================================================='
  WRITE(LUPRI,*) '      NUCLEAR MAGNETIC SHIELDING(NMS) TENSOR RESULTS   ST '
  WRITE(LUPRI,*) '=========================================================='

  !#############################################################################
  !##      Calc RHS of eq 70 JCP 115 page 10349                               ## 
  !##      Note that we have used the notation: FX = hX+GX(D)                 ##
  !##                                                                         ## 
  !##      Using First formula                                                ##
  !##                                                                         ## 
  !##       ! RHS = -P_anti(FxDS + FDSx + FD0xS + G(D0x)DS)                   ##
  !##                                                                         ## 
  !#############################################################################

  call mat_init(SX(1),nbast,nbast)
  call mat_init(SX(2),nbast,nbast)
  call mat_init(SX(3),nbast,nbast)                                     !# Matrices Allocated 3
  call II_get_magderivOverlap(Sx,molcfg%setting,lupri,luerr)

  !RHSX = FDSx
  call mat_init(tempm1,nbast,nbast)                                    !# Matrices Allocated 4
  call mat_mul(F(1),D(1),'n','n',1.0E0_realk,0.0E0_realk,tempm1)
  do icoor = 1,3 
     call mat_init(RHS(icoor),nbast,nbast)
     call mat_mul(tempm1,SX(icoor),'n','n',1.0E0_realk,0.0E0_realk,RHS(icoor))
  enddo                                                                !# Matrices Allocated 7

  !Generate D0X = -D*SX*D
  do icoor = 1,3 
     call mat_mul(SX(icoor),D(1),'n','n',1E0_realk,0.E0_realk,tempm1)     !tempm1 = SX*D
     call mat_free(SX(icoor))                    
     call mat_init(DX(icoor),nbast,nbast)        
     call mat_mul(D(1),tempm1,'n','n',-1.0E0_realk,0.E0_realk,DX(icoor))  !DX = -D*SX*D
  enddo
  call mat_free(tempm1)                                                !# Matrices Allocated 6 (DX,RHS)

  call mat_init(GbDs(1),nbast,nbast)
  call mat_init(GbDs(2),nbast,nbast)  
  call mat_init(GbDs(3),nbast,nbast)                                   !# Matrices Allocated 9 (DX,RHS,GbDs)

  ! Generate G(DX):  The 2-e contribution to sigma vector in RSP
  ! G(DX) = J(DX) + K(DX)
  call di_GET_GbDs(lupri,luerr,DX,GbDs,3,molcfg%setting)
  !RHSX = G(D0X)DS
  call mat_init(tempm1,nbast,nbast)                                    !# Matrices Allocated 10 (DX,RHS,GbDs,tempm1)
  call mat_mul(D(1),S,'n','n',1.0E0_realk,0.0E0_realk,tempm1)
  do icoor = 1,3 
     call mat_mul(GbDs(icoor),tempm1,'n','n',1.0E0_realk,1.0E0_realk,RHS(icoor))
  enddo
  call mat_free(tempm1)

  IF(molcfg%setting%do_dft)THEN
     ! Generate G(DX):  The xc cont to the linear response
     ! Gxc(DX)
     call II_get_xc_linrsp(lupri,luerr,molcfg%setting,nbast,DX,D(1),GbDs,3) 
     ! RHSX = Gxc(D0X)DS
     call mat_init(tempm1,nbast,nbast)
     call mat_mul(D(1),S,'n','n',1.0E0_realk,0.0E0_realk,tempm1)
     do icoor = 1,3 
        call mat_mul(GbDs(icoor),tempm1,'n','n',1.0E0_realk,1.0E0_realk,RHS(icoor))
     enddo
     call mat_free(tempm1)
  ENDIF

  !  Calculate the one electron Magnetic derivative Fock matrix contribution (dh/dBX) to RHS
  call II_get_prop(LUPRI,LUERR,molcfg%SETTING,GbDs,3,'MAGMOM ')
  ! [dh/dBX,S]_D) 
  call mat_init(tempm1,nbast,nbast)
  call mat_mul(D(1),S,'n','n',1.0E0_realk,0.0E0_realk,tempm1)
  do icoor = 1,3 
     call mat_mul(GbDs(icoor),tempm1,'n','n',1.0E0_realk,1.0E0_realk,RHS(icoor))
  enddo
  call mat_free(tempm1)

  !  Calculate the two electron Magnetic derivative Coulomb matrix contribution to RHS
  call II_get_magderivJ(LUPRI,LUERR,molcfg%SETTING,nbast,D,GbDs)
  ! [JX,S]_D) 
  call mat_init(tempm1,nbast,nbast)
  call mat_mul(D(1),S,'n','n',1.0E0_realk,0.0E0_realk,tempm1)
  do icoor = 1,3 
     call mat_mul(GbDs(icoor),tempm1,'n','n',1.0E0_realk,1.0E0_realk,RHS(icoor))
  enddo
  call mat_free(tempm1)

  !  Calculate the two electron Magnetic derivative Exchange matrix contribution to RHS
  call II_get_magderivK(LUPRI,LUERR,molcfg%SETTING,nbast,D,GbDs)
  ! [KX,S]_D) 
  call mat_init(tempm1,nbast,nbast)
  call mat_mul(D(1),S,'n','n',1.0E0_realk,0.0E0_realk,tempm1)
  do icoor = 1,3 
     call mat_mul(GbDs(icoor),tempm1,'n','n',1.0E0_realk,1.0E0_realk,RHS(icoor))
  enddo
  call mat_free(tempm1)

  if(molcfg%setting%do_dft)then
     !  Adding DFT contribution to G^x(D)  
     call II_get_xc_magderiv_kohnsham_mat(LUPRI,LUERR,molcfg%SETTING,nbast,D(1),GbDs)
     ! [FxcX,S]_D) 
     call mat_init(tempm1,nbast,nbast)
     call mat_mul(D(1),S,'n','n',1.0E0_realk,0.0E0_realk,tempm1)
     do icoor = 1,3 
        call mat_mul(GbDs(icoor),tempm1,'n','n',1.0E0_realk,1.0E0_realk,RHS(icoor))
     enddo
     call mat_free(tempm1)
  endif

  call mat_free(GbDs(1))
  call mat_free(GbDs(2))
  call mat_free(GbDs(3))                                        !# Matrices Allocated 6 (DX,RHS)

  do icoor = 1,3 
     call mat_scal(-1.0d0,RHS(icoor))
     !factor 1/2 outside is taken care of in util_get...
     call util_get_symm_part(RHS(icoor))
  enddo

  !#############################################################################
  !##      Solve K([D,Xa]s) = RHS is now ready                                ## 
  !##      Solve to get Da=[D,Xa]_S + D0^a --> Xa (Eq. 70)                    ##
  !## We can use the RSP solver, since K([D,Xa]s) = sigma = -1/2 * E[2]Xa     ##
  !#############################################################################

  do icoor = 1,3 

     eival(1)=0.0E0_realk
     write(lupri,*)'Calling rsp solver for Xa  '
     n_rhs=1
     call mat_init(Xx(1),nbast,nbast)                            !# Matrices Allocated 7 (DX,RHS,Xx)
     if ( mat_dotproduct(RHS(icoor),RHS(icoor))>1.0d-10) then 
        call util_scriptPx('T',D(1),S,RHS(icoor))
        call rsp_init(1,1,1,1,1)
        call rsp_solver(molcfg,D(1),S,F(1),.true.,n_rhs,RHS(icoor:icoor),EIVAL,Xx)
     else
        write(lupri,*) 'WARNING: RHS norm is less than threshold'
        write(lupri,*) 'LIN RSP equations NOT solved for this RHS    '
        call mat_zero(Xx(1))
     end if
     call mat_free(RHS(icoor))

     !############################################################################
     !##      STEP 2: Make D^b=D_0^b+[D_0,X^b]_s                                ## 
     !############################################################################   
     !Generate [D_0,X^b]_s
     call mat_init(tempm1,nbast,nbast)
     call ABCcommutator(nbast,D(1),Xx(1),S,tempm1)
     call mat_free(Xx(1))
     call mat_daxpy(-4.0d0,tempm1,Dx(icoor))
     call mat_free(tempm1)

  enddo !B-field komponen
                                                                 !# Matrices Allocated 3 (DX)
  ! now all 3 DX is known and lies in Dx(1),Dx(2) and Dx(3) 

  !############################################################################
  !##      STEP 3: Make NMST = Tr D^b*h^k + Tr D*h^kb                        ## 
  !############################################################################

  write(lupri,*) '---------------------------------------------------------'
  write(lupri,*) 'NUCLEA MAGNETIC TENSOR                                   '
  write(lupri,*) '---------------------------------------------------------'

  !Now compute the nuclear magnetic shield tensor (ACCORDING TO EQ 83.)
  !NMST = tr D* h^kb + tr D^b*h^k

  !#############################################################################
  !##     Tr D^b*h^k
  !#############################################################################

  !h^k er det det samme som (.PSO) paramagnetic spin-orbit integrals se ref. 61
  !jcp 96:6120, 1992 (PSO abc )

  call mem_alloc(NMST,3*NATOMS,3)
  !expval(3,NATOM,ndmat) X,Y,Z comp for each atom for each B derivate Density Matrix
  call II_get_prop_expval(LUPRI,LUERR,molcfg%SETTING,NMST,Dx,3,3*NATOMS,'PSO    ')  
  do icoor=1,3
     call mat_free(DX(icoor))
  enddo

  !#############################################################################
  !##      Tr D* h^kb
  !#############################################################################

  Factor=53.2513539566280 !1e6*alpha^2 
  !expval(3,3*NATOM) (magnetic X,Y,Z, atomic moment X,Y,Z, for each Atom)
  call mem_alloc(expval,3,3*NATOMS)
  call II_get_prop_expval(LUPRI,LUERR,molcfg%SETTING,expval,D,1,9*NATOMS,'NST    ')
  do icoor=1,3
     do jcoor=1,3*natoms  ! magnetic moment koordinate
        NMST(jcoor,icoor) = 2*factor*(- NMST(jcoor,icoor) + expval(icoor,jcoor))
     enddo
  enddo
  call mem_dealloc(expval)

  allocate(atomname(natoms))
  do jcoor=1,natoms  
     atomname(jcoor)=molcfg%setting%molecule(1)%p%ATOM(jcoor)%Name
  enddo

  WRITE(LUPRI,*) " Total shielding tensor"

  do jcoor=1,natoms  
     WRITE (LUPRI,'(20X,3(A,13X),/)') 'Bx', 'By', 'Bz'
     WRITE (LUPRI,'(2X,A8,3F15.8)')  atomname(jcoor)//CHRXYZ(1),(NMST(3*jCOOR-2,K),K=1,3)
     WRITE (LUPRI,'(2X,A8,3F15.8)')  atomname(jcoor)//CHRXYZ(2),(NMST(3*jCOOR-1,K),K=1,3)
     WRITE (LUPRI,'(2X,A8,3F15.8)')  atomname(jcoor)//CHRXYZ(3),(NMST(3*jCOOR,K),K=1,3)
  enddo

  WRITE(LUPRI,*) "      Absolute chemical shift"
  WRITE(LUPRI,'(A,I12)') "number of atoms:",natoms
  do jcoor=1,natoms  
   WRITE (LUPRI,'(2X,A8,f15.8,A8)')  atomname(jcoor), &
        & 1.0E0_realk/3.0E0_realk*(NMST(3*jCOOR-2,1)+NMST(3*jCOOR-1,2)+NMST(3*jCOOR,3))
  enddo
  deallocate(atomname)
  call mem_dealloc(NMST)
  WRITE(LUPRI,*) " Done with shielding tensor calculation"

end subroutine NMRshieldresponse_noOpenRSP

end module response_noOpenRSP_module
