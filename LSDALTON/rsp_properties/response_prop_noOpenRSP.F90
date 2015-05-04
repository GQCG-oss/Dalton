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
  use nuclei_selected_shielding_mod
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

  public lsdalton_response_noOpenRSP

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
        WRITE(config%LUPRI,*)'NMRshieldresponse_noOpenRSP'
        call LSTIMER('START',t1,t2,config%LUPRI)
        call NMRshieldresponse_noOpenRSP(molcfg,F,D,S)
  !      WRITE(config%LUPRI,*)'NMRshieldresponse_noOpenRSP_TK'
  !      call NMRshieldresponse_noOpenRSP_TK(molcfg,F,D,S)
        call LSTIMER('NMRshield',t1,t2,config%LUPRI)
     endif
     if(config%response%tasks%doNMRshield_selected) then
!        WRITE(config%LUPRI,*)'NMRshieldresponse_IANS'
        call LSTIMER('START',t1,t2,config%LUPRI)     
!        WRITE(config%LUPRI,*)'NMRshieldresponse_IANS_TK'
!        call NMRshieldresponse_IANS_TK(ls,molcfg,F,D,S)
!        WRITE(config%LUPRI,*)'NMRshieldresponse_NS'
!        call NMRshieldresponse_NS(ls,molcfg,F,D,S)
        WRITE(config%LUPRI,*)'NMRshieldresponse_RSTNS'
        call NMRshieldresponse_RSTNS(ls,molcfg,F,D,S)
        call LSTIMER('NMRshield',t1,t2,config%LUPRI)
     endif
     if(config%response%tasks%doALPHA)then
        call alpha_driver(molcfg,F,D,S,config%response%alphainput)
     endif
     if(config%response%tasks%doOPA)then
        !one-photon absorption driver (OPA)
        call excitation_energy_driver(molcfg,F,D,S,config%response%alphainput)
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

!> \brief Driver for calculating polarizabilites for frequencies defined in input.
!> If no frequencies are specified the frequency is set to zero.
!> \author Thomas Kjaergaard
!> \date 2014
subroutine alpha_driver(molcfg,F,D,S,alphainput)
  implicit none
  !> Info on molecule needed by solver and integral programs
  type(rsp_molcfg),intent(inout) :: molcfg
  !> Unperturbed Fock, Density and Overlap matrix
  type(Matrix),intent(in)        :: F(1),D(1),S
  !> Contains alpha input
  type(ALPHAinputItem),intent(inout) :: ALPHAinput
  !Local variables
  type(Matrix)         :: Bgrad(3),DIPLEN(3),XSOL(1)
  type(Matrix)         :: tempm1,tempm2
  real(realk),pointer  :: linrspfunc(:,:,:) !Array for response functions
  logical              :: lineq_x
  integer              :: i,j,ifreq,nfreq,ndim,lupri,luerr
  integer              :: ntrial,nrhs,nsol,nomega,nstart
  ndim = S%nrow
  lupri = molcfg%lupri
  luerr = molcfg%lupri

  IF(alphainput%imag_frequencies_in_input)THEN
     CALL LSQUIT('Error Imaginary frequencies in alpha_driver',-1)
  ENDIF

  nfreq = alphainput%nfreq
  lineq_x = .TRUE.

  call mat_init(DIPLEN(1),ndim,ndim)
  call mat_init(DIPLEN(2),ndim,ndim)
  call mat_init(DIPLEN(3),ndim,ndim)
  call II_get_integral(LUPRI,LUERR,molcfg%SETTING,DIPLEN,3,'DIPLEN ')
  !BDS-SDB
  call mat_init(tempm1,ndim,ndim)
  call mat_mul(D(1),S,'n','n',1.0E0_realk,0.0E0_realk,tempm1)
  DO I = 1,3
     call mat_init(Bgrad(I),ndim,ndim)
     call mat_mul(DIPLEN(I),tempm1,'N','N',1.0E0_realk,0.0E0_realk,Bgrad(I))
     call mat_mul(tempm1,DIPLEN(I),'T','N',-1.0E0_realk,1.0E0_realk,Bgrad(I))
     call mat_free(DIPLEN(I))
     call util_scriptPx('T',D(1),S,Bgrad(I))
  ENDDO
  call mat_free(tempm1)

  call mat_init(XSOL(1),ndim,ndim)
  call mem_alloc(linrspfunc,3,nfreq,3)
  do i = 1, 3
   !FIXME try to solve response equations for more frequencies at a time
   if ( mat_dotproduct(Bgrad(i),Bgrad(i)) > 1.0E-10_realk) then 
      do ifreq = 1, nfreq
         !Solve Linear response equation
         ntrial = 1 !# of trial vectors in a given iteration (number of RHS)
         nrhs = 1   !# of RHS only relevant for linear equations (lineq_x = TRUE)
         nsol = 1   !# of solution (output) vectors
         nomega = 1 !If lineq_x, number of laser freqs (input)
                    !Otherwise number of excitation energies (output) 
         nstart = 1 !Number of start vectors. Only relevant for eigenvalue problem
         !ntrial and nstart seem to be obsolete 
         call rsp_init(ntrial,nrhs,nsol,nomega,nstart)
         call rsp_solver(molcfg,D(1),S,F(1),lineq_x,nrhs,Bgrad(I:I),&
              & alphainput%bfreq(ifreq:ifreq),XSOL)
         ! <<GD_A,GD_B>>_omega
         DO J = 1,3        
            linrspfunc(J,ifreq,I) = -4.E0_realk*mat_dotproduct(Bgrad(J),XSOL(1))
         ENDDO
      enddo
   else
      !XSOL is zero
      do ifreq = 1, nfreq
         DO J = 1,3        
            linrspfunc(J,ifreq,I) = 0.0E0_realk
         ENDDO
      enddo
   endif
  enddo
  call mat_free(XSOL(1))
  do i = 1, 3
     call mat_free(Bgrad(i))
  enddo
  
  ! Print the results
  call print_alpha(linrspfunc,ALPHAinput,nfreq,lupri)
  call mem_dealloc(linrspfunc)
  
end subroutine alpha_driver

!> \brief Driver routine for excitation energy & trans. moments 
!> \author Thomas Kjaergaard
!> \date 2014
subroutine excitation_energy_driver(molcfg,F,D,S,alphainput)
  implicit none
  !> Info on molecule needed by solver and integral programs
  type(rsp_molcfg),intent(inout) :: molcfg
  !> Unperturbed Fock, Density and Overlap matrix
  type(Matrix),intent(in)        :: F(1),D(1),S
  !> Contains alpha input
  type(ALPHAinputItem),intent(inout) :: ALPHAinput
  !Local variables
  type(Matrix)         :: Bgrad(3),DIPLEN(3),GDB(1)
  type(Matrix),pointer :: XSOL(:)
  type(Matrix)         :: tempm1,tempm2
  real(realk),pointer  :: eival(:),ExciMoments(:,:)
  logical              :: lineq_x
  integer              :: i,j,ifreq,nfreq,ndim,lupri,luerr,nexci,lueigvec
  integer              :: ntrial,nrhs,nsol,nomega,nstart

  ! 1. Set number of excitation energies
  ! ''''''''''''''''''''''''''''''''''''
  ! (Defined by .NEXCIT in input)
  nexci = molcfg%decomp%cfg_rsp_nexcit

  ndim = S%nrow
  lupri = molcfg%lupri
  luerr = molcfg%lupri
  lineq_x = .FALSE.

  call mem_alloc(xsol,nexci)
  do i = 1,nexci
     call mat_init(xsol(i),ndim,ndim)
  enddo
  call mem_alloc(eival,nexci)

  ntrial = 1     !# of trial vectors in a given iteration (number of RHS)
  nrhs = 1       !# of RHS only relevant for linear equations (lineq_x = TRUE)
  nsol = nexci   !# of solution (output) vectors
  nomega = nexci !If lineq_x, number of laser freqs (input)
                 !Otherwise number of excitation energies (output) 
  nstart = nexci !Number of start vectors. Only relevant for eigenvalue problem

  !ntrial and nstart seem to be obsolete 
  call rsp_init(ntrial,nrhs,nsol,nomega,nstart)
  call rsp_solver(molcfg,D(1),S,F(1),lineq_x,nexci,GDB,eival,xsol)

!  lueigvec = -1
!  CALL GPOPEN(lueigvec,'rsp_eigenvecs','unknown','SEQUENTIAL','UNFORMATTED',i,j)
!  do i = 1, nexci
!     write(lueigvec) i, EIVAL(i)
!     call mat_write_to_disk(lueigvec,xsol(i))
!  enddo
!  call GPCLOSE(lueigvec,'KEEP') !Might be used for exc state stuff

  call mat_init(DIPLEN(1),ndim,ndim)
  call mat_init(DIPLEN(2),ndim,ndim)
  call mat_init(DIPLEN(3),ndim,ndim)
  call II_get_integral(LUPRI,LUERR,molcfg%SETTING,DIPLEN,3,'DIPLEN ')
  call mat_init(tempm1,ndim,ndim)
  call mat_mul(D(1),S,'n','n',1.0E0_realk,0.0E0_realk,tempm1)
  DO I = 1,3
     call mat_init(Bgrad(I),ndim,ndim)
     call mat_mul(DIPLEN(I),tempm1,'N','N',1.0E0_realk,0.0E0_realk,Bgrad(I))
     call mat_mul(tempm1,DIPLEN(I),'T','N',-1.0E0_realk,1.0E0_realk,Bgrad(I))
     call mat_free(DIPLEN(I))
     call util_scriptPx('T',D(1),S,Bgrad(I))
  ENDDO
  call mat_free(tempm1)
  call mem_alloc(ExciMoments,nexci,3)
  DO J = 1,3        
     do I = 1, nexci     
        ExciMoments(I,J) = -2.E0_realk*mat_dotproduct(Bgrad(J),XSOL(I))
     enddo
     call mat_free(Bgrad(J))
  ENDDO
  do i = 1,nexci
     call mat_free(xsol(i))
  enddo
  call mem_dealloc(xsol)
  call print_excit(ExciMoments,eival,AlphaInput,nexci,lupri)
  call mem_dealloc(eival)
  call mem_dealloc(ExciMoments)
  
end subroutine excitation_energy_driver

subroutine print_excit(ExciMoments,e_excit,AlphaInput,nexci,lupri)
  implicit none
  !> Contains alpha input
  type(ALPHAinputItem),intent(inout) :: AlphaInput
  integer, intent(in) :: nexci,lupri
  real(realk), intent(in) :: ExciMoments(nexci,3),e_excit(nexci)
  integer :: i,j
  real(realk) :: OscillatorStrength(nexci)

  do i=1,nexci
     OscillatorStrength(i) = 0E0_realk
     do j=1,3
        OscillatorStrength(i) = OscillatorStrength(i) + ExciMoments(i,j)**2
     end do
     OscillatorStrength(i) = (2E0_realk/3E0_realk)*e_excit(i)*OscillatorStrength(i)
  enddo

  write(lupri,*) 
  write(lupri,*) 
  write(lupri,'(2X,A)') '********************************************************&
       &**********************'
  write(lupri,'(2X,A)') '*                   ONE-PHOTON ABSORPTION RESULTS (in a.u.)&
       &                  *'
  write(lupri,'(2X,A)') '********************************************************&
       &**********************'
  write(lupri,*)
  write(lupri,*) 
  write(lupri,*) 
  write(lupri,*) '     Excitation              Transition Dipole Moments      &
       &         Oscillator'
  write(lupri,*) '      Energies            x               y               z  &
       &         Strengths'
  write(lupri,*) '===================================================================&
       &============'
  
  do i=1,nexci
     write(lupri,'(5f16.8)') e_excit(i), ExciMoments(i,1:3), OscillatorStrength(i)
  end do
  write(lupri,*) 
  write(lupri,*) 
  write(lupri,'(1X,A)') ' End of excitation energy calculation'
end subroutine print_excit

subroutine print_alpha(alpha,ALPHAinput,nfreq,lupri)
  implicit none
  !> Contains alpha input
  type(ALPHAinputItem),intent(inout) :: ALPHAinput
  integer, intent(in) :: nfreq,lupri
  real(realk), intent(in) :: alpha(3,nfreq,3)
  integer :: i
  write(lupri,*) 
  write(lupri,*) 
  write(lupri,*) 
  write(lupri,'(1X,A)') '*************************************************************'
  write(lupri,'(1X,A)') '*          POLARIZABILITY TENSOR RESULTS (in a.u.)          *'
  write(lupri,'(1X,A)') '*************************************************************'
  write(lupri,*) 
  write(lupri,*)
  write(lupri,*) 
  write(lupri,*) 
  PrintFrequencyLoop: do i=1,nfreq

     write(lupri,'(1X,A,g12.5)')  'Frequency = ', alphainput%bfreq(i)
     write(lupri,*) '=================='
     write(lupri,*) 
     write(lupri,*) "                Ex               Ey               Ez"
     write(lupri,'(1X,A, 3g18.8)') "Ex    ", real(alpha(1,i,:))
     write(lupri,'(1X,A, 3g18.8)') "Ey    ", real(alpha(2,i,:))
     write(lupri,'(1X,A, 3g18.8)') "Ez    ", real(alpha(3,i,:))
     write(lupri,*) 
     write(lupri,'(1X,A,g18.8)') 'Isotropic polarizability = ', &
          & ( alpha(1,i,1) + alpha(2,i,2) + alpha(3,i,3) )/3E0_realk

     write(lupri,*) 
     write(lupri,*) 
     write(lupri,*) 
     write(lupri,*) 
     write(lupri,*) 
     
  end do PrintFrequencyLoop
  write(lupri,'(1X,A)') ' End of polarizability calculation'
end subroutine print_alpha


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
!  type(Matrix),pointer :: ChandanMat(:)
  integer              :: ntrial,nrhs,nsol,nomega,nstart
  character(len=1)        :: CHRXYZ(-3:3)
  DATA CHRXYZ /'z','y','x',' ','X','Y','Z'/
  Factor=53.2513539566280 !1e6*alpha^2 

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
     call mat_mul(tempm1,SX(icoor),'n','n',1.0E0_realk,0.0E0_realk,RHS(icoor))  !RHS = FDSx
  enddo                                                                !# Matrices Allocated 7

  !Generate D0X = -D*SX*D
  do icoor = 1,3 
     call mat_mul(SX(icoor),D(1),'n','n',1E0_realk,0.E0_realk,tempm1)     !tempm1 = SX*D
     call mat_free(SX(icoor))                    
     call mat_init(DX(icoor),nbast,nbast)        
     call mat_mul(D(1),tempm1,'n','n',-1.0E0_realk,0.E0_realk,DX(icoor))  !DX = -D*SX*D
  enddo

  do icoor = 1,3 
     call mat_mul(F(1),Dx(icoor),'n','n',1.0E0_realk,0.0E0_realk,tempm1)
     call mat_mul(tempm1,S,'n','n',1.0E0_realk,1.E0_realk,RHS(icoor))    !RHS = FDSx + FD0xS
  enddo
  call mat_free(tempm1) 


!  INFO
!  call mem_alloc(NMST,3*NATOMS,3)
!  !expval(3,NATOM,ndmat) X,Y,Z comp for each atom for each B derivate Density Matrix
!  call II_get_prop_expval(LUPRI,LUERR,molcfg%SETTING,NMST,Dx,3,3*NATOMS,'PSO    ')  
!  do icoor=1,3
!     do jcoor=1,3*natoms  ! magnetic moment koordinate
!        WRITE(lupri,*)'PURE D0b*hk NMST:',- 2*factor*NMST(jcoor,icoor)
!     enddo
!  enddo
!  call mem_dealloc(NMST)

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

  !Calculate the two electron Magnetic derivative Exchange matrix contribution to RHS
  call mat_zero(GbDs(1))
  call mat_zero(GbDs(2))
  call mat_zero(GbDs(3))
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
     call util_get_symm_part(RHS(icoor)) !Eq. 60 in the paper
  enddo

  !#############################################################################
  !##      Solve K([D,Xa]s) = RHS is now ready                                ## 
  !##      Solve to get Da=[D,Xa]_S + D0^a --> Xa (Eq. 70)                    ##
  !## We can use the RSP solver, since K([D,Xa]s) = sigma = -1/2 * E[2]Xa     ##
  !#############################################################################
  !FIXME try to solve response equations for more RHS at a time
  do icoor = 1,3 

     eival(1)=0.0E0_realk
     write(lupri,*)'Calling rsp solver for Xa  '
     call mat_init(Xx(1),nbast,nbast)                            !# Matrices Allocated 7 (DX,RHS,Xx)
     call util_scriptPx('T',D(1),S,RHS(icoor))
     if ( mat_dotproduct(RHS(icoor),RHS(icoor))>1.0d-10) then 

        ntrial = 1 !# of trial vectors in a given iteration (number of RHS)
        nrhs = 1   !# of RHS only relevant for linear equations (lineq_x = TRUE)
        nsol = 1   !# of solution (output) vectors
        nomega = 1 !If lineq_x, number of laser freqs (input)
                   !Otherwise number of excitation energies (output) 
        nstart = 1 !Number of start vectors. Only relevant for eigenvalue problem
        !ntrial and nstart seem to be obsolete 
        call rsp_init(ntrial,nrhs,nsol,nomega,nstart)
        call rsp_solver(molcfg,D(1),S,F(1),.true.,nrhs,RHS(icoor:icoor),EIVAL,Xx)
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

!     call mat_init(GbDs(icoor),nbast,nbast)
!     call ABCcommutator(nbast,D(1),Xx(1),S,GbDs(icoor))
!     call mat_free(Xx(1))
!     call mat_daxpy(-4.0d0,GbDs(icoor),Dx(icoor))
  enddo !B-field komponen
                                                                 !# Matrices Allocated 3 (DX)
!  call mem_alloc(NMST,3*NATOMS,3)
!  do icoor=1,3
!     !expval(3,NATOM,ndmat) X,Y,Z comp for each atom for each B derivate Density Matrix
!     call II_get_prop_expval(LUPRI,LUERR,molcfg%SETTING,NMST,GbDs,3,3*NATOMS,'PSO    ')  
!     do jcoor=1,3*natoms  ! magnetic moment koordinate
!        WRITE(lupri,*)'PURE [D,Xb]*hk NMST:',- 2*factor*NMST(jcoor,icoor)
!     enddo
!  enddo
!  call mem_dealloc(NMST)
!  call mat_free(GbDs(1))
!  call mat_free(GbDs(2))
!  call mat_free(GbDs(3))

  ! now all 3 DX is known and lies in Dx(1),Dx(2) and Dx(3) 
  



  !############################################################################
  !##      STEP 3: Make NMST = Tr D^b*h^k + Tr D*h^kb                        ## 
  !############################################################################

  write(lupri,*) '---------------------------------------------------------'
  write(lupri,*) 'NUCLEA MAGNETIC TENSOR                                   '
  write(lupri,*) '---------------------------------------------------------'


  !Now compute the nuclear magnetic shield tensor (ACCORDING TO EQ 83.)
  !NMST = tr D* h^kb + tr D^b*h^k +Tr D*h^bk
 
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

  !expval(3,3*NATOM) (magnetic X,Y,Z, atomic moment X,Y,Z, for each Atom)
   call mem_alloc(expval,3,3*NATOMS)
   call II_get_prop_expval(LUPRI,LUERR,molcfg%SETTING,expval,D,1,9*NATOMS,'NST    ')
  do icoor=1,3
     do jcoor=1,3*natoms  ! magnetic moment koordinate
!        WRITE(lupri,*)'NMST:',- 2*factor*NMST(jcoor,icoor), 2*factor*expval(icoor,jcoor),'=',&
!             & 2*factor*(- NMST(jcoor,icoor) + expval(icoor,jcoor))
        NMST(jcoor,icoor) = 2*factor*(- NMST(jcoor,icoor) +expval(icoor,jcoor))
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

  !##################################################################
  !  Subroutine NMRshieldresponse_noOpenRSP_TK
  !  Compute shielding tensor in conventional manner - 'Xb response'. 
  !  This code is for debugging purpose. 
  !
  !
  !##################################################################
subroutine NMRshieldresponse_noOpenRSP_TK(molcfg,F,D,S)
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
  type(Matrix)                 :: D0x(3)
!  type(Matrix),pointer :: ChandanMat(:)
  integer              :: ntrial,nrhs,nsol,nomega,nstart
  character(len=1)        :: CHRXYZ(-3:3)
  DATA CHRXYZ /'z','y','x',' ','X','Y','Z'/
  Factor=53.2513539566280 !1e6*alpha^2 

  nbast = D(1)%nrow
  lupri = molcfg%lupri
  luerr = lupri
  natoms = molcfg%natoms
  CALL LSTIMER('START ',TS,TE,LUPRI)

  allocate(atomname(natoms))    
  do jcoor=1,natoms  
     atomname(jcoor)=molcfg%setting%molecule(1)%p%ATOM(jcoor)%Name
  enddo

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

  call mem_alloc(NMST,3*NATOMS,3)
  !#############################################################################
  !##      Tr D* h^kb
  !#############################################################################

  !expval(3,3*NATOM) (magnetic X,Y,Z, atomic moment X,Y,Z, for each Atom)
  call mem_alloc(expval,3,3*NATOMS)
  call II_get_prop_expval(LUPRI,LUERR,molcfg%SETTING,expval,D,1,9*NATOMS,'NST    ')
  WRITE(lupri,*)'TTT NMST(Tr D* h^kb):'
  do icoor=1,3
     do jcoor=1,3*natoms  ! magnetic moment koordinate
        NMST(jcoor,icoor) = 2*factor*expval(icoor,jcoor)
        WRITE(lupri,*)'NMST(Tr D* h^kb):',NMST(jcoor,icoor)
     enddo
  enddo
  call mem_dealloc(expval)

  call mem_alloc(expval,3*NATOMS,3)
  !Generate D0X = -D*SX*D
  do icoor = 1,3 
     call mat_init(SX(icoor),nbast,nbast)
     call mat_init(DX(icoor),nbast,nbast)        
     call mat_init(D0X(icoor),nbast,nbast)        
     call mat_init(RHS(icoor),nbast,nbast)
  enddo

  call II_get_magderivOverlap(Sx,molcfg%setting,lupri,luerr)
  call mat_init(tempm1,nbast,nbast)       
  do icoor = 1,3 
     call mat_mul(SX(icoor),D(1),'n','n',1E0_realk,0.E0_realk,tempm1)     !tempm1 = SX*D
     call mat_mul(D(1),tempm1,'n','n',-1.0E0_realk,0.E0_realk,DX(icoor))  !DX = -D*SX*D
     call mat_mul(D(1),tempm1,'n','n',-1.0E0_realk,0.E0_realk,D0X(icoor))  !DX = -D*SX*D
  enddo
  call mat_free(tempm1) 
  call II_get_prop_expval(LUPRI,LUERR,molcfg%SETTING,expval,D0x,3,3*NATOMS,'PSO    ')  
  WRITE(lupri,*)'TTT NMST(Tr D0b*h^k):'
  do icoor=1,3
     do jcoor=1,3*natoms  ! magnetic moment koordinate
        WRITE(lupri,*)'NMST(Tr D0b*h^k):',- 2*factor*expval(jcoor,icoor)
        NMST(jcoor,icoor) = NMST(jcoor,icoor) - 2*factor*expval(jcoor,icoor)
     enddo
  enddo

  !RHSX = FD0xS
  call mat_init(tempm1,nbast,nbast)       
  do icoor = 1,3 
     call mat_mul(F(1),D0x(icoor),'n','n',1.0E0_realk,0.0E0_realk,tempm1)
     call mat_mul(tempm1,S,'n','n',1.0E0_realk,0.E0_realk,RHS(icoor))   
       
     WRITE(lupri,*)'Printing RHS'
     call mat_print(RHS(icoor),1,RHS(icoor)%nrow,1,RHS(icoor)%ncol,lupri)
  enddo
  call mat_free(tempm1) 
  call GetXfromRHS(molcfg,RHS,Dx,D,S,F,lupri)
  call II_get_prop_expval(LUPRI,LUERR,molcfg%SETTING,expval,Dx,3,3*NATOMS,'PSO    ')  
  WRITE(lupri,*)'TTT '
  do icoor=1,3
     do jcoor=1,3*natoms  ! magnetic moment koordinate
        WRITE(lupri,*)'NMST(Tr X(FD0bS)*h^k):',- 2*factor*expval(jcoor,icoor)
        NMST(jcoor,icoor) = NMST(jcoor,icoor) - 2*factor*expval(jcoor,icoor) 
     enddo
  enddo

  !RHSX = FDSx
  call mat_init(tempm1,nbast,nbast)       
  call mat_mul(F(1),D(1),'n','n',1.0E0_realk,0.0E0_realk,tempm1)
  do icoor = 1,3 
     call mat_mul(tempm1,SX(icoor),'n','n',1.0E0_realk,0.0E0_realk,RHS(icoor))  !RHS = FDSx
  enddo                                   
  call mat_free(tempm1) 
  call GetXfromRHS(molcfg,RHS,Dx,D,S,F,lupri)
  call II_get_prop_expval(LUPRI,LUERR,molcfg%SETTING,expval,Dx,3,3*NATOMS,'PSO    ')  
  WRITE(lupri,*)'TTT '
  do icoor=1,3
     do jcoor=1,3*natoms  ! magnetic moment koordinate
        WRITE(lupri,*)'NMST(Tr X(FDSb)*h^k):',- 2*factor*expval(jcoor,icoor)
        NMST(jcoor,icoor) = NMST(jcoor,icoor) - 2*factor*expval(jcoor,icoor) 
     enddo
  enddo
  
  call mat_init(GbDs(1),nbast,nbast)
  call mat_init(GbDs(2),nbast,nbast)  
  call mat_init(GbDs(3),nbast,nbast)                                   
  call di_GET_GbDs(lupri,luerr,D0X,GbDs,3,molcfg%setting)
  !RHSX = G(D0b)DS
  call mat_init(tempm1,nbast,nbast)                                    
  call mat_mul(D(1),S,'n','n',1.0E0_realk,0.0E0_realk,tempm1)
  do icoor = 1,3 
     call mat_mul(GbDs(icoor),tempm1,'n','n',1.0E0_realk,0.0E0_realk,RHS(icoor))
  enddo
  call mat_free(tempm1)
  call GetXfromRHS(molcfg,RHS,Dx,D,S,F,lupri)
  call II_get_prop_expval(LUPRI,LUERR,molcfg%SETTING,expval,Dx,3,3*NATOMS,'PSO    ')  
  WRITE(lupri,*)'TTT '
  do icoor=1,3
     do jcoor=1,3*natoms  ! magnetic moment koordinate
        WRITE(lupri,*)'NMST(Tr X(G(D0b)DS)*h^k):',- 2*factor*expval(jcoor,icoor)
        NMST(jcoor,icoor) = NMST(jcoor,icoor) - 2*factor*expval(jcoor,icoor)
     enddo
  enddo

  IF(molcfg%setting%do_dft)THEN
     ! Generate G(DX):  The xc cont to the linear response
     ! Gxc(DX)
     call II_get_xc_linrsp(lupri,luerr,molcfg%setting,nbast,D0X,D(1),GbDs,3) 
     ! RHSX = Gxc(D0X)DS
     call mat_init(tempm1,nbast,nbast)
     call mat_mul(D(1),S,'n','n',1.0E0_realk,0.0E0_realk,tempm1)
     do icoor = 1,3 
        call mat_mul(GbDs(icoor),tempm1,'n','n',1.0E0_realk,0.0E0_realk,RHS(icoor))
     enddo
     call mat_free(tempm1)
     call GetXfromRHS(molcfg,RHS,Dx,D,S,F,lupri)
     call II_get_prop_expval(LUPRI,LUERR,molcfg%SETTING,expval,Dx,3,3*NATOMS,'PSO    ')  
  WRITE(lupri,*)'TTT '
     do icoor=1,3
        do jcoor=1,3*natoms  ! magnetic moment koordinate
           WRITE(lupri,*)'NMST(Tr X(Gxc(D0b)DS)*h^k):',- 2*factor*expval(jcoor,icoor)
           NMST(jcoor,icoor) = NMST(jcoor,icoor) - 2*factor*expval(jcoor,icoor) 
        enddo
     enddo
  ENDIF

  !  Calculate the one electron Magnetic derivative Fock matrix contribution (dh/dBX) to RHS
  call II_get_prop(LUPRI,LUERR,molcfg%SETTING,GbDs,3,'MAGMOM ')
  ! [dh/dBX,S]_D) 
  call mat_init(tempm1,nbast,nbast)
  call mat_mul(D(1),S,'n','n',1.0E0_realk,0.0E0_realk,tempm1)
  do icoor = 1,3 
     call mat_mul(GbDs(icoor),tempm1,'n','n',1.0E0_realk,0.0E0_realk,RHS(icoor))
  enddo
  call mat_free(tempm1)
  call GetXfromRHS(molcfg,RHS,Dx,D,S,F,lupri)
  call II_get_prop_expval(LUPRI,LUERR,molcfg%SETTING,expval,Dx,3,3*NATOMS,'PSO    ')  
  WRITE(lupri,*)'TTT '
  do icoor=1,3
     do jcoor=1,3*natoms  ! magnetic moment koordinate
        WRITE(lupri,*)'NMST(Tr X(hb*DS)*h^k):',- 2*factor*expval(jcoor,icoor)
        NMST(jcoor,icoor) = NMST(jcoor,icoor) - 2*factor*expval(jcoor,icoor)
     enddo
  enddo

  !  Calculate the two electron Magnetic derivative Coulomb matrix contribution to RHS
  call II_get_magderivJ(LUPRI,LUERR,molcfg%SETTING,nbast,D,GbDs)
  ! [JX,S]_D) 
  call mat_init(tempm1,nbast,nbast)
  call mat_mul(D(1),S,'n','n',1.0E0_realk,0.0E0_realk,tempm1)
  do icoor = 1,3 
     call mat_mul(GbDs(icoor),tempm1,'n','n',1.0E0_realk,0.0E0_realk,RHS(icoor))
  enddo
  call mat_free(tempm1)
  call GetXfromRHS(molcfg,RHS,Dx,D,S,F,lupri)
  call II_get_prop_expval(LUPRI,LUERR,molcfg%SETTING,expval,Dx,3,3*NATOMS,'PSO    ')  
  WRITE(lupri,*)'TTT '
  do icoor=1,3
     do jcoor=1,3*natoms  ! magnetic moment koordinate
        WRITE(lupri,*)'NMST(Tr X(Jb(D)*DS)*h^k):',- 2*factor*expval(jcoor,icoor)
        NMST(jcoor,icoor) = NMST(jcoor,icoor) - 2*factor*expval(jcoor,icoor)
     enddo
  enddo

  !Calculate the two electron Magnetic derivative Exchange matrix contribution to RHS
  do icoor=1,3
    call mat_zero(GbDs(icoor))
  enddo
  call II_get_magderivK(LUPRI,LUERR,molcfg%SETTING,nbast,D,GbDs)
  ! [KX,S]_D) 
  call mat_init(tempm1,nbast,nbast)
  call mat_mul(D(1),S,'n','n',1.0E0_realk,0.0E0_realk,tempm1)
  do icoor = 1,3 
     call mat_mul(GbDs(icoor),tempm1,'n','n',1.0E0_realk,0.0E0_realk,RHS(icoor))
  enddo
  call mat_free(tempm1)
  call GetXfromRHS(molcfg,RHS,Dx,D,S,F,lupri)
  call II_get_prop_expval(LUPRI,LUERR,molcfg%SETTING,expval,Dx,3,3*NATOMS,'PSO    ')  
  WRITE(lupri,*)'TTT '
  do icoor=1,3
     do jcoor=1,3*natoms  ! magnetic moment koordinate
        WRITE(lupri,*)'NMST(Tr X(Kb(D)*DS)*h^k):',- 2*factor*expval(jcoor,icoor)
        NMST(jcoor,icoor) = NMST(jcoor,icoor) - 2*factor*expval(jcoor,icoor) 
     enddo
  enddo

  if(molcfg%setting%do_dft)then
     !  Adding DFT contribution to G^x(D)  
     call II_get_xc_magderiv_kohnsham_mat(LUPRI,LUERR,molcfg%SETTING,nbast,D(1),GbDs)
     ! [FxcX,S]_D) 
     call mat_init(tempm1,nbast,nbast)
     call mat_mul(D(1),S,'n','n',1.0E0_realk,0.0E0_realk,tempm1)
     do icoor = 1,3 
        call mat_mul(GbDs(icoor),tempm1,'n','n',1.0E0_realk,0.0E0_realk,RHS(icoor))
     enddo
     call mat_free(tempm1)
     call GetXfromRHS(molcfg,RHS,Dx,D,S,F,lupri)
     call II_get_prop_expval(LUPRI,LUERR,molcfg%SETTING,expval,Dx,3,3*NATOMS,'PSO    ')  
     WRITE(lupri,*)'TTT '
     do icoor=1,3
        do jcoor=1,3*natoms  ! magnetic moment koordinate
           WRITE(lupri,*)'NMST(Tr X(Gxcb(D)*DS)*h^k):',- 2*factor*expval(jcoor,icoor)
           NMST(jcoor,icoor) = NMST(jcoor,icoor) - 2*factor*expval(jcoor,icoor) 
        enddo
     enddo
  endif

  do icoor = 1,3 
     call mat_free(Dx(icoor))
     call mat_free(D0x(icoor))
     call mat_free(GbDs(icoor))
     call mat_free(RHS(icoor))
  enddo
  call mem_dealloc(expval)

  write(lupri,*) '---------------------------------------------------------'
  write(lupri,*) 'FINAL NUCLEA MAGNETIC TENSOR                             '
  write(lupri,*) '---------------------------------------------------------'

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

end subroutine NMRshieldresponse_noOpenRSP_TK

subroutine GetXfromRHS(molcfg,RHS,Dx,D,S,F,lupri)
  implicit none
  type(rsp_molcfg), intent(inout) :: molcfg
  integer :: lupri
  type(matrix),intent(in) :: D(1),S,F(1)
  type(matrix),intent(inout) :: RHS(3)
  type(matrix),intent(inout) :: Dx(3)
  !
  real(realk) :: eival(1)    
  type(matrix) :: Xx(1),tempm1
  integer :: ntrial,nrhs,nsol,nomega,nstart,icoor,nbast
  nbast = S%nrow
  do icoor = 1,3
     call mat_zero(DX(icoor))
     call mat_scal(-1.0d0,RHS(icoor))
     call util_get_symm_part(RHS(icoor)) !Eq. 60 in the paper
  enddo
  
  !#############################################################################
  !##      Solve K([D,Xa]s) = RHS is now ready                                ## 
  !##      Solve to get Da=[D,Xa]_S + D0^a --> Xa (Eq. 70)                    ##
  !## We can use the RSP solver, since K([D,Xa]s) = sigma = -1/2 * E[2]Xa     ##
  !#############################################################################
  !FIXME try to solve response equations for more RHS at a time
  do icoor = 1,3 
     
     eival(1)=0.0E0_realk
     write(lupri,*)'Calling rsp solver for Xa  '
     call mat_init(Xx(1),nbast,nbast)                            !# Matrices Allocated 7 (DX,RHS,Xx)
     call util_scriptPx('T',D(1),S,RHS(icoor))
     if ( mat_dotproduct(RHS(icoor),RHS(icoor))>1.0d-10) then 
        ntrial = 1 !# of trial vectors in a given iteration (number of RHS)
        nrhs = 1   !# of RHS only relevant for linear equations (lineq_x = TRUE)
        nsol = 1   !# of solution (output) vectors
        nomega = 1 !If lineq_x, number of laser freqs (input)
        !Otherwise number of excitation energies (output) 
        nstart = 1 !Number of start vectors. Only relevant for eigenvalue problem
        !ntrial and nstart seem to be obsolete 
        call rsp_init(ntrial,nrhs,nsol,nomega,nstart)
        call rsp_solver(molcfg,D(1),S,F(1),.true.,nrhs,RHS(icoor:icoor),EIVAL,Xx)
     else
        write(lupri,*) 'WARNING: RHS norm is less than threshold'
        write(lupri,*) 'LIN RSP equations NOT solved for this RHS    '
        call mat_zero(Xx(1))
     end if
     
     !############################################################################
     !##      STEP 2: Make D^b=D_0^b+[D_0,X^b]_s                                ## 
     !############################################################################   
     !Generate [D_0,X^b]_s
     
     call mat_init(tempm1,nbast,nbast)
     call ABCcommutator(nbast,D(1),Xx(1),S,tempm1)
     call mat_free(Xx(1))
     call mat_daxpy(-4.0d0,tempm1,Dx(icoor))
     call mat_free(tempm1)
     call mat_zero(RHS(icoor))
  enddo
end subroutine GetXfromRHS


end module response_noOpenRSP_module
