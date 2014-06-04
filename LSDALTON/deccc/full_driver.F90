!> @file
!> Full calculation
!> This file is mainly a playground for new developments, not intended to be included in a release.

module full 
  use fundamental
  use precision
  use typedeftype!,only:lsitem
  use typedef
  use dec_typedef_module
  use matrix_module
  use matrix_operations
  use memory_handling

  ! DEC DEPENDENCIES (within deccc directory)   
  ! *****************************************
  use dec_fragment_utils
  use CABS_operations
#ifdef MOD_UNRELEASED
  use cc_debug_routines_module
  use full_f12contractions
  use f12_routines_module   ! Moved to August 2013 by Yang M. Wang
#endif
  use array4_simple_operations
  use array3_simple_operations
  use array2_simple_operations
  use mp2_module
  !  use orbital_operations
  use full_molecule
  use ccintegrals!,only: get_full_AO_integrals,get_AO_hJ,get_AO_K,get_AO_Fock
  use ccdriver!,only: ccsolver_justenergy, ccsolver
  use fragment_energy_module,only : Full_DECMP2_calculation

  public :: full_driver
  private

contains

  !> \brief Main part for full molecular coupled-cluster calculations.
  subroutine full_driver(MyMolecule,mylsitem,D)

    implicit none
    !> Full molecule structure
    type(fullmolecule), intent(inout) :: MyMolecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> HF density matrix
    type(matrix),intent(in) :: D
    real(realk) :: Ecorr,EHF,Eerr

    write(DECinfo%output,'(/,a)') ' ================================================ '
    write(DECinfo%output,'(a)')   '              Full molecular driver               '
    write(DECinfo%output,'(a,/)') ' ================================================ '

#ifdef VAR_MPI
    call set_dec_settings_on_slaves()
#endif

    ! run cc program
    if(DECinfo%F12) then ! F12 correction
#ifdef MOD_UNRELEASED
       if(DECinfo%ccModel==MODEL_MP2) then
          call full_canonical_mp2_f12(MyMolecule,MyLsitem,D,Ecorr)
       else
          call full_get_ccsd_f12_energy(MyMolecule,MyLsitem,D,Ecorr)
       end if
#else
       call lsquit('f12 not released',-1)
#endif
    else
       !if(DECinfo%ccModel==MODEL_MP2) then

       !   if(DECinfo%use_canonical ) then
       !      ! simple conventional MP2 calculation, only works for canonical orbitals
       !      call full_canonical_mp2_correlation_energy(MyMolecule,mylsitem,Ecorr)
       !   else
       !      ! Call routine which calculates individual fragment contributions and prints them,
       !      ! works both for canonical and local orbitals
       !      call Full_DECMP2_calculation(MyMolecule,mylsitem,Ecorr)
       !   end if

       !else
          call full_cc_dispatch(MyMolecule,mylsitem,Ecorr)          
       !end if
    end if

    ! Get HF energy
    Ehf = get_HF_energy_fullmolecule(MyMolecule,Mylsitem,D)

    ! Print summary
    Eerr = 0.0_realk   ! zero error for full calculation by definition
    call print_total_energy_summary(EHF,Ecorr,Eerr)

  end subroutine full_driver

  !> \brief Dispatch cc job for full molecule
  subroutine full_cc_dispatch(MyMolecule,mylsitem,Ecorr)

    implicit none
    !> FUll molecule info
    type(fullmolecule), intent(inout) :: MyMolecule
    !> Integral info
    type(lsitem), intent(inout) :: mylsitem
    !> Correlation energy
    real(realk),intent(inout) :: Ecorr
    integer :: nocc,nunocc,nbasis,print_level,i,j
    logical :: fragment_job
    real(realk),pointer :: ppfock_fc(:,:), Co_fc(:,:)

    Ecorr = 0.0E0_realk

    if(DECinfo%FrozenCore) then
       nocc = MyMolecule%nval
    else
       nocc = MyMolecule%nocc
    end if
    nunocc = MyMolecule%nunocc
    nbasis = MyMolecule%nbasis

    fragment_job = .false.
    print_level = 2 ! this is not used

    if(DECinfo%frozencore) then
       ! Pass only valence orbitals
       call mem_alloc(ppfock_fc,nocc,nocc)
       do j=1,nocc
          do i=1,nocc
             ppfock_fc(i,j) = MyMolecule%ppfock(MyMolecule%ncore+i,MyMolecule%ncore+j)
          end do
       end do

       ! Frozen core component of MO coeff.
       call mem_alloc(Co_fc,nbasis,nocc)
       do j=1,nocc
          do i=1,nbasis
             Co_fc(i,j) = MyMolecule%Co(i,MyMolecule%ncore+j)
          end do
       end do

       Ecorr = ccsolver_justenergy(DECinfo%ccmodel,MyMolecule,nbasis,nocc,nunocc,&
          & mylsitem,print_level,fragment_job,Co_fc=Co_fc,ppfock_fc=ppfock_fc)

       call mem_dealloc(ppfock_fc)
       call mem_dealloc(Co_fc)

    else


#ifdef MOD_UNRELEASED
       if(DECinfo%CCSDmultipliers)then
          call ccsolver_energy_multipliers(DECinfo%ccmodel,MyMolecule%Co,MyMolecule%Cv,&
             & MyMolecule%fock, nbasis,nocc,nunocc,mylsitem, &
             & print_level,fragment_job,MyMolecule%ppfock,MyMolecule%qqfock,ecorr)
       else
          Ecorr = ccsolver_justenergy(DECinfo%ccmodel,MyMolecule,nbasis,nocc,nunocc,&
             & mylsitem,print_level,fragment_job)
       endif
#else
       Ecorr = ccsolver_justenergy(DECinfo%ccmodel,MyMolecule,nbasis,nocc,nunocc,&
             & mylsitem,print_level,fragment_job)
#endif

    end if

  end subroutine full_cc_dispatch

#ifdef MOD_UNRELEASED
  !> \brief Calculate canonical MP2 energy for full molecular system
  !> keeping full AO integrals in memory. Only for testing.
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine full_canonical_mp2_f12(MyMolecule,MyLsitem,Dmat,mp2f12_energy)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> HF density matrix
    type(matrix),intent(in) :: Dmat
    !> Canonical MP2-F12 correlation energy
    real(realk),intent(inout) :: mp2f12_energy
    !> Canonical MP2 correlation energy
    real(realk) :: mp2_energy
    !> E22 energies
    real(realk) :: X1,X2,X3,X4
    !> E23 energies
    real(realk) :: B1,B2,B3,B4,B5,B6,B7,B8,B9
    
    real(realk),pointer :: gao(:,:,:,:)
    real(realk),pointer :: gmo(:,:,:,:)
    real(realk),pointer :: Ripjq(:,:,:,:)
    real(realk),pointer :: Fijkl(:,:,:,:)
    real(realk),pointer :: Tijkl(:,:,:,:)
    real(realk),pointer :: Rimjc(:,:,:,:)
    real(realk),pointer :: Dijkl(:,:,:,:)
    real(realk),pointer :: Tirjk(:,:,:,:)
    real(realk),pointer :: Tijkr(:,:,:,:)
    real(realk),pointer :: Gipjq(:,:,:,:)
    real(realk),pointer :: Gimjc(:,:,:,:)
    real(realk),pointer :: Girjs(:,:,:,:)
    real(realk),pointer :: Girjm(:,:,:,:)
    real(realk),pointer :: Grimj(:,:,:,:)

    real(realk),pointer :: Gipja(:,:,:,:)
    real(realk),pointer :: Gpiaj(:,:,:,:)

    real(realk),pointer :: Gicjm(:,:,:,:)
    real(realk),pointer :: Gcimj(:,:,:,:)
    real(realk),pointer :: Gcirj(:,:,:,:)
    real(realk),pointer :: Gciaj(:,:,:,:)
    real(realk),pointer :: Giajc(:,:,:,:)
    real(realk),pointer :: Ciajb(:,:,:,:)
    real(realk),pointer :: Cjaib(:,:,:,:)
    real(realk),pointer :: Taibj(:,:,:,:) !amplitudes not integrals

    real(realk),pointer :: Vijkl_term3(:,:,:,:)

    real(realk),pointer :: Vijij(:,:)
    real(realk),pointer :: Vijij_term1(:,:)
    real(realk),pointer :: Vijij_term2(:,:)
    real(realk),pointer :: Vijij_term3(:,:)
    real(realk),pointer :: Vijij_term4(:,:)
    real(realk),pointer :: Vijij_term5(:,:)

    real(realk),pointer :: Vjiij(:,:)    
    real(realk),pointer :: Vjiij_term1(:,:)
    real(realk),pointer :: Vjiij_term2(:,:)
    real(realk),pointer :: Vjiij_term3(:,:)
    real(realk),pointer :: Vjiij_term4(:,:)
    real(realk),pointer :: Vjiij_term5(:,:)

    real(realk),pointer :: Xijkl(:,:,:,:)
    real(realk),pointer :: Xijkl_term1(:,:,:,:)
    real(realk),pointer :: Xijkl_term2(:,:,:,:)
    real(realk),pointer :: Xijkl_term3(:,:,:,:)
    real(realk),pointer :: Xijkl_term4(:,:,:,:)
  
    real(realk),pointer :: Xijij(:,:)
    real(realk),pointer :: Xijij_term1(:,:)
    real(realk),pointer :: Xijij_term2(:,:)
    real(realk),pointer :: Xijij_term3(:,:)
    real(realk),pointer :: Xijij_term4(:,:)

    real(realk),pointer :: Xjiij(:,:)    
    real(realk),pointer :: Xjiij_term1(:,:)
    real(realk),pointer :: Xjiij_term2(:,:)
    real(realk),pointer :: Xjiij_term3(:,:)
    real(realk),pointer :: Xjiij_term4(:,:)

    real(realk),pointer :: Bijij(:,:)
    real(realk),pointer :: Bijij_term1(:,:)
    real(realk),pointer :: Bijij_term2(:,:)
    real(realk),pointer :: Bijij_term3(:,:)
    real(realk),pointer :: Bijij_term4(:,:)
    real(realk),pointer :: Bijij_term5(:,:) 
    real(realk),pointer :: Bijij_term6(:,:)
    real(realk),pointer :: Bijij_term7(:,:)
    real(realk),pointer :: Bijij_term8(:,:)
    real(realk),pointer :: Bijij_term9(:,:)

    real(realk),pointer :: Bjiij(:,:)
    real(realk),pointer :: Bjiij_term1(:,:)
    real(realk),pointer :: Bjiij_term2(:,:)    
    real(realk),pointer :: Bjiij_term3(:,:)
    real(realk),pointer :: Bjiij_term4(:,:)
    real(realk),pointer :: Bjiij_term5(:,:)
    real(realk),pointer :: Bjiij_term6(:,:)
    real(realk),pointer :: Bjiij_term7(:,:)
    real(realk),pointer :: Bjiij_term8(:,:)
    real(realk),pointer :: Bjiij_term9(:,:)

    real(realk),pointer :: Bijij_debug(:,:)
    real(realk),pointer :: Bjiij_debug(:,:)

    integer :: nbasis,ncabs,nocc,nvirt,I,A,B,J,noccfull,ncabsAO
    integer :: l,m,p,q,c,r
    real(realk) :: tmp, V3energy

    real(realk) :: eps
    character :: string(4)
    type(matrix) :: K
    !    type(matrix) :: HJrc
    type(matrix) :: HJir
    !    type(matrix) :: Kcc
    type(matrix) :: Krr
    !    type(matrix) :: Fcc
    type(matrix) :: Frr
    type(matrix) :: Frc
    type(matrix) :: Fpp
    type(matrix) :: Fmm
    type(matrix) :: Frm
    type(matrix) :: Fcp
    type(matrix) :: Fii
    type(matrix) :: Fac
    Real(realk)  :: E21, E21_debug, E22, E22_debug, E23_debug, Gtmp
    type(array4) :: array4Taibj,array4gmo
    !    logical :: fulldriver 
    !    fulldriver = .TRUE.
    !    call init_cabs(fulldriver)

    ! Init stuff
    ! **********
    nbasis = MyMolecule%nbasis
    nocc   = MyMolecule%nocc
    nvirt  = MyMolecule%nunocc
    call determine_CABS_nbast(ncabsAO,ncabs,mylsitem%setting,DECinfo%output)
    noccfull = nocc

    ! Memory check!
    ! ********************
    call full_canonical_mp2_memory_check(nbasis,nocc,nvirt)

    ! Get all F12 Fock Matrices
    ! ********************
    call get_F12_mixed_MO_Matrices(MyLsitem,MyMolecule,Dmat,nbasis,ncabsAO,&
         & nocc,noccfull,nvirt,ncabs,HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp)

    ! Get all AO integrals
    ! ********************
    call mem_alloc(gao,nbasis,nbasis,nbasis,nbasis)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRC')
    ! Transform AO integrals to MO integrals (A I | B J)
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'aiai',gAO,gMO)
    call mem_dealloc(gao)

    call get_4Center_F12_integrals(mylsitem,MyMolecule,nbasis,nocc,noccfull,nvirt,ncabsAO,&
         & Ripjq,Fijkl,Tijkl,Rimjc,Dijkl,Tirjk,Tijkr,Gipjq,Gimjc,Girjs,Girjm,&
         & Grimj,Gipja,Gpiaj,Gicjm,Gcimj,Gcirj,Gciaj,Giajc)

    call mem_alloc(Vijij,nocc,nocc)
    call mem_alloc(Vjiij,nocc,nocc)

    call mp2f12_Vijij(Vijij,Ripjq,Gipjq,Fijkl,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)
    call mp2f12_Vjiij(Vjiij,Ripjq,Gipjq,Fijkl,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)

    call mem_alloc(Ciajb,nocc,nvirt,nocc,nvirt)    
    !   call mem_alloc(Cjaib,nocc,nvirt,nocc,nvirt)
    call mp2f12_Ciajb(Ciajb,Giajc,Fac%elms,nocc,nvirt,ncabs)
    !   call mp2f12_Cjaib(Cjaib,Giajc,Fac%elms,nocc,nvirt,ncabs)
    
    if(DECinfo%use_canonical) then
       !construct canonical T amplitudes
       call mem_alloc(Taibj,nvirt,nocc,nvirt,nocc)
       do J=1,nocc
          do B=1,nvirt
             do I=1,nocc
                do A=1,nvirt
                   ! Difference in orbital energies: eps(I) + eps(J) - eps(A) - eps(B)
                   eps = MyMolecule%ppfock(I,I) + MyMolecule%ppfock(J,J) &
                        & - MyMolecule%qqfock(A,A) - MyMolecule%qqfock(B,B)
                   eps = gmo(A,I,B,J)/eps
                   Taibj(a,i,b,j) = eps
                enddo
             enddo
          enddo
       enddo
       ! Calculate canonical MP2 energy
       ! ******************************
       mp2_energy = 0.0E0_realk
       do J=1,nocc
          do B=1,nvirt
             do I=1,nocc
                do A=1,nvirt

                   ! Difference in orbital energies: eps(I) + eps(J) - eps(A) - eps(B)
                   eps = MyMolecule%ppfock(I,I) + MyMolecule%ppfock(J,J) &
                        & - MyMolecule%qqfock(A,A) - MyMolecule%qqfock(B,B)

                   ! Energy = sum_{AIBJ} (AI|BJ) * [ 2(AI|BJ) - (BI|AJ) ] / (epsI + epsJ - epsA - epsB)
                   mp2_energy = mp2_energy + gmo(A,I,B,J)*(2E0_realk*gmo(A,I,B,J)-gmo(B,I,A,J))/eps

                end do
             end do
          end do
       end do

    else
       !  THIS PIECE OF CODE IS MORE GENERAL AS IT DOES NOT REQUIRE CANONICAL ORBITALS
       !    ! Get full MP2 (as specified in input)

       ! KK: Quick and dirty solution to the fact that the MP2 solver requires array4 format.
       array4gmo = array4_init([nvirt,nocc,nvirt,nocc])
       array4gmo%val=gmo
       call mp2_solver(nocc,nvirt,MyMolecule%ppfock,MyMolecule%qqfock,array4gmo,array4Taibj)
       call array4_free(array4gmo)

       call mem_alloc(Taibj,nvirt,nocc,nvirt,nocc)

       mp2_energy = 0.0E0_realk
       do j=1,nocc
          do b=1,nvirt
             do i=1,nocc
                do a=1,nvirt
                   ! Energy = sum_{ijab} ( Taibj) * (ai | bj)
                   Taibj(a,i,b,j) = array4Taibj%val(a,i,b,j)
                   Gtmp = 2.0E0_realk * gmo(a,i,b,j) - gmo(b,i,a,j)
                   mp2_energy = mp2_energy + Taibj(a,i,b,j) * Gtmp
                end do
             end do
          end do
       end do
    endif
  
    call mp2f12_Vijij_coupling(Vijij,Ciajb,Taibj,nocc,nvirt)
    call mp2f12_Vjiij_coupling(Vjiij,Ciajb,Taibj,nocc,nvirt)

    !> Calculate E21 Energy
    E21 = 2.0E0_REALK*mp2f12_E21(Vijij,Vjiij,nocc)

    if(DECinfo%F12DEBUG) then    
       call mem_alloc(Vijij_term1,nocc,nocc)
       call mem_alloc(Vijij_term2,nocc,nocc)
       call mem_alloc(Vijij_term3,nocc,nocc)
       call mem_alloc(Vijij_term4,nocc,nocc)
       call mem_alloc(Vijij_term5,nocc,nocc)

       call mem_alloc(Vjiij_term1,nocc,nocc)
       call mem_alloc(Vjiij_term2,nocc,nocc)
       call mem_alloc(Vjiij_term3,nocc,nocc)
       call mem_alloc(Vjiij_term4,nocc,nocc)
       call mem_alloc(Vjiij_term5,nocc,nocc)

       call mp2f12_Vijij_term1(Vijij_term1,Fijkl,nocc,noccfull,nbasis,ncabs)
       call mp2f12_Vijij_term2(Vijij_term2,Ripjq,Gipjq,nocc,noccfull,nbasis,ncabs)
       call mp2f12_Vijij_term3(Vijij_term3,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)
       call mp2f12_Vijij_term4(Vijij_term4,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)

       call mp2f12_Vjiij_term1(Vjiij_term1,Fijkl,nocc,noccfull,nbasis,ncabs)
       call mp2f12_Vjiij_term2(Vjiij_term2,Ripjq,Gipjq,nocc,noccfull,nbasis,ncabs)
       call mp2f12_Vjiij_term3(Vjiij_term3,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)
       call mp2f12_Vjiij_term4(Vjiij_term4,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)

       !> Coupling with the C-matrix, only needs to be done once
       call mp2f12_Vijij_term5(Vijij_term5,Ciajb,Taibj,nocc,nvirt)
       call mp2f12_Vjiij_term5(Vjiij_term5,Ciajb,Taibj,nocc,nvirt)
        
       print *, '----------------------------------------'
       print *, '           C - matrix terms             '
       print *, '----------------------------------------'
       print *, 'norm4D(Ciajb): ', norm4D(Ciajb)
       print *, 'norm4D(Giajc): ', norm4D(Giajc)
       print *, 'norm4D(Taibj): ', norm4D(Taibj)
       print *, '----------------------------------------'
       print *, '           V - matrix terms             '
       print *, '----------------------------------------'
       print *,'norm4D(Fijkl): ', norm4D(Fijkl)
       print *,'norm4D(Ripjq): ', norm4D(Ripjq)
       print *,'norm4D(Gipjq): ', norm4D(Gipjq)
       print *, '----------------------------------------'
       print *,'norm4D(Rimjc): ', norm4D(Rimjc)
       print *,'norm4D(Gimjc): ', norm4D(Gimjc)
       print *, '----------------------------------------'
       print *, '           E21  V terms                 '
       print *, '----------------------------------------'
       print *, 'E21_V_term1: ', 2.0E0_REALK*mp2f12_E21(Vijij_term1,Vjiij_term1,nocc)
       print *, 'E21_V_term2: ', 2.0E0_REALK*mp2f12_E21(Vijij_term2,Vjiij_term2,nocc)
       print *, 'E21_V_term3: ', 2.0E0_REALK*mp2f12_E21(Vijij_term3,Vjiij_term3,nocc)
       print *, 'E21_V_term4: ', 2.0E0_REALK*mp2f12_E21(Vijij_term4,Vjiij_term4,nocc)
       print *, 'E21_V_term5: ', 2.0E0_REALK*mp2f12_E21(Vijij_term5,Vjiij_term5,nocc)
       print *, '----------------------------------------'

       E21_debug = 2.0E0_REALK*(mp2f12_E21(Vijij_term1,Vjiij_term1,nocc) + mp2f12_E21(Vijij_term2,Vjiij_term2,nocc) &
            & + mp2f12_E21(Vijij_term3,Vjiij_term3,nocc) + mp2f12_E21(Vijij_term4,Vjiij_term4,nocc) &
            & + mp2f12_E21(Vijij_term5,Vjiij_term5,nocc)) 

       print *, 'E21_Vsum: ', E21_debug
       print *, 'E21_debug: ', 2.0E0_REALK*mp2f12_E21(Vijij,Vjiij,nocc)
    endif

    call mem_dealloc(Vijij)
    call mem_dealloc(Vjiij)   
    call mem_dealloc(Taibj)
    call mem_dealloc(Ciajb)
    
    if(DECinfo%F12DEBUG) then
       call mem_dealloc(Vijij_term1)
       call mem_dealloc(Vijij_term2)
       call mem_dealloc(Vijij_term3)
       call mem_dealloc(Vijij_term4)
       call mem_dealloc(Vijij_term5)

       call mem_dealloc(Vjiij_term1)
       call mem_dealloc(Vjiij_term2)
       call mem_dealloc(Vjiij_term3)
       call mem_dealloc(Vjiij_term4)
       call mem_dealloc(Vjiij_term5)
    endif

    if(DECinfo%F12DEBUG) then
       call mem_alloc(Xijkl_term1,nocc,nocc,nocc,nocc)
       call mem_alloc(Xijkl_term2,nocc,nocc,nocc,nocc)
       call mem_alloc(Xijkl_term3,nocc,nocc,nocc,nocc)
       call mem_alloc(Xijkl_term4,nocc,nocc,nocc,nocc) 

       call mem_alloc(Xijij_term1,nocc,nocc)
       call mem_alloc(Xijij_term2,nocc,nocc)
       call mem_alloc(Xijij_term3,nocc,nocc)
       call mem_alloc(Xijij_term4,nocc,nocc)      

       call mem_alloc(Xjiij_term1,nocc,nocc)
       call mem_alloc(Xjiij_term2,nocc,nocc)
       call mem_alloc(Xjiij_term3,nocc,nocc)
       call mem_alloc(Xjiij_term4,nocc,nocc)

       call mem_alloc(Bijij_term1,nocc,nocc)
       call mem_alloc(Bijij_term2,nocc,nocc)
       call mem_alloc(Bijij_term3,nocc,nocc)
       call mem_alloc(Bijij_term4,nocc,nocc)   
       call mem_alloc(Bijij_term5,nocc,nocc)
       call mem_alloc(Bijij_term6,nocc,nocc)
       call mem_alloc(Bijij_term7,nocc,nocc)
       call mem_alloc(Bijij_term8,nocc,nocc)
       call mem_alloc(Bijij_term9,nocc,nocc)

       call mem_alloc(Bjiij_term1,nocc,nocc)
       call mem_alloc(Bjiij_term2,nocc,nocc)
       call mem_alloc(Bjiij_term3,nocc,nocc)
       call mem_alloc(Bjiij_term4,nocc,nocc)   
       call mem_alloc(Bjiij_term5,nocc,nocc)
       call mem_alloc(Bjiij_term6,nocc,nocc)
       call mem_alloc(Bjiij_term7,nocc,nocc)
       call mem_alloc(Bjiij_term8,nocc,nocc)
       call mem_alloc(Bjiij_term9,nocc,nocc)

       if(DECinfo%use_canonical) then 
          call mp2f12_Xijij_term1(Xijij_term1,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
          call mp2f12_Xijij_term2(Xijij_term2,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
          call mp2f12_Xijij_term3(Xijij_term3,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
          call mp2f12_Xijij_term4(Xijij_term4,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)

          call mp2f12_Xjiij_term1(Xjiij_term1,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
          call mp2f12_Xjiij_term2(Xjiij_term2,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
          call mp2f12_Xjiij_term3(Xjiij_term3,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
          call mp2f12_Xjiij_term4(Xjiij_term4,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)

          print *,'-----------------------------------------'
          print *,'          X - matrix terms               '
          print *,'-----------------------------------------'
          print *,'norm2D(Xijij_term1): ', norm2D(Xijij_term1)
          print *,'norm2D(Xijij_term2): ', norm2D(Xijij_term2)
          print *,'norm2D(Xijij_term3): ', norm2D(Xijij_term3)
          print *,'norm2D(Xijij_term4): ', norm2D(Xijij_term4)
          print *,'-----------------------------------------'

       else !> Non canonical

          call mp2f12_Xijijfull_term1(Xijkl_term1,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
          call mp2f12_Xijijfull_term2(Xijkl_term2,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
          call mp2f12_Xijijfull_term3(Xijkl_term3,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
          call mp2f12_Xijijfull_term4(Xijkl_term4,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)

          print *,'-----------------------------------------'
          print *,'          X - matrix terms               '
          print *,'-----------------------------------------'
          print *,'norm4D(Xijkl_term1): ', norm4D(Xijkl_term1)
          print *,'norm4D(Xijkl_term2): ', norm4D(Xijkl_term2)
          print *,'norm4D(Xijkl_term3): ', norm4D(Xijkl_term3)
          print *,'norm4D(Xijkl_term4): ', norm4D(Xijkl_term4)
          print *,'-----------------------------------------'

       endif

       call mp2f12_Bijij_term1(Bijij_term1,Bjiij_term1,nocc,Dijkl)
       call mp2f12_Bijij_term2(Bijij_term2,Bjiij_term2,nocc,ncabsAO,Tirjk,hJir%elms)
       call mp2f12_Bijij_term3(Bijij_term3,Bjiij_term3,nocc,ncabsAO,Tijkr,hJir%elms)    
       call mp2f12_Bijij_term4(Bijij_term4,Bjiij_term4,nocc,noccfull,ncabsAO,Girjs,Krr%elms)

       call mp2f12_Bijij_term5(Bijij_term5,Bjiij_term5,nocc,noccfull,ncabsAO,Girjm,Grimj,Frr%elms)
       call mp2f12_Bijij_term6(Bijij_term6,Bjiij_term6,nocc,noccfull,ncabsAO,nvirt,nbasis,Gipja,Gpiaj,Fpp%elms)
       call mp2f12_Bijij_term7(Bijij_term7,Bjiij_term7,nocc,noccfull,ncabs,Gicjm,Gcimj,Fmm%elms)
       call mp2f12_Bijij_term8(Bijij_term8,Bjiij_term8,nocc,noccfull,ncabsAO,ncabs,Gicjm,Gcirj,Frm%elms)
       call mp2f12_Bijij_term9(Bijij_term9,Bjiij_term9,nocc,noccfull,nvirt,ncabs,nbasis,Gipja,Gciaj,Fcp%elms)

       print *,'-----------------------------------------'
       print *,'         B - matrix terms                '
       print *,'-----------------------------------------'
       print *, '(B1 Term):'
       print *,'-----------------------------------------'
       print *,'norm4D(Dijkl): ', norm4D(Dijkl)
       print *,'-----------------------------------------'
       print *, '(B2 Term):'
       print *,'-----------------------------------------'
       print *,'norm4D(Tirjk): ', norm4D(Tirjk)
       print *,'-----------------------------------------'
       print *, '(B3 Term):'
       print *,'-----------------------------------------'
       print *,'norm4D(Tijkr): ', norm4D(Tijkr)
       print *,'-----------------------------------------'
       print *, '(B4 Term):'
       print *,'-----------------------------------------'
       print *,'norm4D(Girjs): ', norm4D(Girjs)
       print *,'-----------------------------------------'
       print *, '(B5 Term):'
       print *,'-----------------------------------------'
       print *,'norm4D(Girjm): ', norm4D(Girjm)
       print *,'norm4D(Grimj): ', norm4D(Grimj)
       print *,'-----------------------------------------'
       print *, '(B6 Term):'
       print *,'-----------------------------------------'
       print *,'norm4D(Gipja): ', norm4D(Gipja)
       print *,'norm4D(Gpiaj): ', norm4D(Gpiaj)
       print *,'-----------------------------------------'
       print *, '(B7 Term):'
       print *,'-----------------------------------------'
       print *,'norm4D(Gicjm): ', norm4D(Gicjm)
       print *,'-----------------------------------------'
       print *, '(B8 Term):'
       print *,'-----------------------------------------'
       print *,'norm4D(Gcirj): ', norm4D(Gcirj)
       print *,'-----------------------------------------'
       print *, '(B9 Term):'
       print *,'-----------------------------------------'
       print *,'norm4D(Gciaj): ', norm4D(Gciaj)
       print *,'-----------------------------------------'
       print *,'norm2D(Bijij_term1): ', norm2D(Bijij_term1)
       print *,'norm2D(Bijij_term2): ', norm2D(Bijij_term2)
       print *,'norm2D(Bijij_term3): ', norm2D(Bijij_term3)
       print *,'norm2D(Bijij_term4): ', norm2D(Bijij_term4)
       print *,'norm2D(Bijij_term5): ', norm2D(Bijij_term5)
       print *,'norm2D(Bijij_term6): ', norm2D(Bijij_term6)
       print *,'norm2D(Bijij_term7): ', norm2D(Bijij_term7)
       print *,'norm2D(Bijij_term8): ', norm2D(Bijij_term8)
       print *,'norm2D(Bijij_term9): ', norm2D(Bijij_term9)      
       print *,'-----------------------------------------'
       print *,'full_canonical_mp2_f12: Get all F12 Fock integrals'
       print *,'-----------------------------------------'
       print *, "norm2D(hJir)", norm1D(hJir%elms)
       print *, "norm2D(Krr)", norm1D(Krr%elms)
       print *, "norm2D(Frr)", norm1D(Frr%elms)
       print *, "norm2D(Fac)", norm1D(Fac%elms)
       print *, "norm2D(Fpp)", norm1D(Fpp%elms)
       print *, "norm2D(Fii)", norm1D(Fii%elms)
       print *, "norm2D(Fmm)", norm1D(Fmm%elms)
       print *, "norm2D(Frm)", norm1D(Frm%elms)
       print *, "norm2D(Fcp)", norm1D(Fcp%elms)
       print *,'-----------------------------------------' 

    endif

    if(DECinfo%use_canonical) then    
       call mem_alloc(Xijij,nocc,nocc)
       call mem_alloc(Xjiij,nocc,nocc)

       call mp2f12_Xijij(Xijij,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
       call mp2f12_Xjiij(Xjiij,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)

    else !> Non - canonical
       
       call mem_alloc(Xijkl,nocc,nocc,nocc,nocc)
       call mp2f12_Xijijfull(Xijkl,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)

    endif

    call mem_alloc(Bijij,nocc,nocc)
    call mem_alloc(Bjiij,nocc,nocc)

    if(DECinfo%F12DEBUG) then
       call mem_alloc(Bjiij_debug,nocc,nocc)
       call mem_alloc(Bijij_debug,nocc,nocc)
    endif

    call mp2f12_Bijij(Bijij,Dijkl,Tirjk,Tijkr,Girjs,hJir%elms,Krr%elms,&
         & Frr%elms,Fpp%elms,Fmm%elms,Frm%elms,Fcp%elms,&
         & Girjm,Grimj,Gipja,Gpiaj,Gicjm,Gcimj,Gcirj,Gciaj,&
         & nocc,noccfull,nbasis,ncabsAO,nvirt,ncabs)

    call mp2f12_Bjiij(Bjiij,Dijkl,Tirjk,Tijkr,Girjs,hJir%elms,Krr%elms,&
         & Frr%elms,Fpp%elms,Fmm%elms,Frm%elms,Fcp%elms,&
         & Girjm,Grimj,Gipja,Gpiaj,Gicjm,Gcimj,Gcirj,Gciaj,&
         & nocc,noccfull,nbasis,ncabsAO,nvirt,ncabs)

    if(DECinfo%use_canonical) then

       if(DECinfo%F12DEBUG) then
          !> Setting Bmatrix = 0
          Bijij_debug = 0.0E0_realk
          Bjiij_debug = 0.0E0_realk  
          call submp2f12_EBX(E22_debug,Bijij_debug,Bjiij_debug,Xijij,Xjiij,Fii%elms,nocc)

       else

          call submp2f12_EBX(E22,Bijij,Bjiij,Xijij,Xjiij,Fii%elms,nocc)

       endif

       if(DECinfo%F12DEBUG) then
          print *, '----------------------------------------'
          print *, '          E_22 X term                   '
          print *, '----------------------------------------'
          print *, 'E22_X_term1: ', mp2f12_E22(Xijij_term1,Xjiij_term1,Fii%elms,nocc)
          print *, 'E22_X_term2: ', mp2f12_E22(Xijij_term2,Xjiij_term2,Fii%elms,nocc)
          print *, 'E22_X_term3: ', mp2f12_E22(Xijij_term3,Xjiij_term3,Fii%elms,nocc)
          print *, 'E22_X_term4: ', mp2f12_E22(Xijij_term4,Xjiij_term4,Fii%elms,nocc)
          print *, '----------------------------------------'
          E22_debug = mp2f12_E22(Xijij_term1,Xjiij_term1,Fii%elms,nocc) & 
               & + mp2f12_E22(Xijij_term2,Xjiij_term2,Fii%elms,nocc) &
               & + mp2f12_E22(Xijij_term3,Xjiij_term3,Fii%elms,nocc) + mp2f12_E22(Xijij_term4,Xjiij_term4,Fii%elms,nocc)  
          print *, 'E22_Xsum: ', E22_debug  
          print *, 'E22_debug: ',  mp2f12_E22(Xijij,Xjiij,Fii%elms,nocc)
          print *, '----------------------------------------'
          print *, '          E_23 B term                   '
          print *, '----------------------------------------'
          print *, 'E23_B_term1: ', mp2f12_E23(Bijij_term1,Bjiij_term1,nocc)
          print *, 'E23_B_term2: ', mp2f12_E23(Bijij_term2,Bjiij_term2,nocc)
          print *, 'E23_B_term3: ', mp2f12_E23(Bijij_term3,Bjiij_term3,nocc)
          print *, 'E23_B_term4: ', mp2f12_E23(Bijij_term4,Bjiij_term4,nocc)
          print *, 'E23_B_term5: ', mp2f12_E23(Bijij_term5,Bjiij_term5,nocc)
          print *, 'E23_B_term6: ', mp2f12_E23(Bijij_term6,Bjiij_term6,nocc)
          print *, 'E23_B_term7: ', mp2f12_E23(Bijij_term7,Bjiij_term7,nocc)
          print *, 'E23_B_term8: ', mp2f12_E23(Bijij_term8,Bjiij_term8,nocc)
          print *, 'E23_B_term9: ', mp2f12_E23(Bijij_term9,Bjiij_term9,nocc)   
          print *, '----------------------------------------'
          E23_debug = mp2f12_E23(Bijij_term1,Bjiij_term1,nocc) & 
               & + mp2f12_E23(Bijij_term2,Bjiij_term2,nocc) + mp2f12_E23(Bijij_term3,Bjiij_term3,nocc) &
               & + mp2f12_E23(Bijij_term4,Bjiij_term4,nocc) + mp2f12_E23(Bijij_term5,Bjiij_term5,nocc) &
               & + mp2f12_E23(Bijij_term6,Bjiij_term6,nocc) + mp2f12_E23(Bijij_term7,Bjiij_term7,nocc) &
               & + mp2f12_E23(Bijij_term8,Bjiij_term8,nocc) + mp2f12_E23(Bijij_term9,Bjiij_term9,nocc)
          print *, 'E23_Bsum: ',  E23_debug
          !print *, 'E23_Bsum_debug: ',  mp2f12_E23(Bijij,Bjiij,nocc)
          print *, '----------------------------------------'
       endif
       
    else !> Non - canoical

       call submp2f12_EBXfull(E22,Bijij,Bjiij,Xijkl,Fii%elms,nocc)
       
       if(DECinfo%F12DEBUG) then
          X1 = mp2f12_E22X(Xijkl_term1,Fii%elms,nocc)
          X2 = mp2f12_E22X(Xijkl_term2,Fii%elms,nocc)
          X3 = mp2f12_E22X(Xijkl_term3,Fii%elms,nocc)
          X4 = mp2f12_E22X(Xijkl_term4,Fii%elms,nocc)
          print *, '----------------------------------------'
          print *, '          E_22 X term                   '
          print *, '----------------------------------------'
          print *, 'E22_X_term1: ', X1
          print *, 'E22_X_term2: ', X2
          print *, 'E22_X_term3: ', X3
          print *, 'E22_X_term4: ', X4
          print *, '----------------------------------------'
          E22_debug = X1 + X2 + X3 + X4  
          print *, 'E22_Xsum: ', E22_debug  
          print *, '----------------------------------------'
          print *, '          E_23 B term                   '
          print *, '----------------------------------------'
          print *, 'E23_B_term1: ', mp2f12_E23(Bijij_term1,Bjiij_term1,nocc)
          print *, 'E23_B_term2: ', mp2f12_E23(Bijij_term2,Bjiij_term2,nocc)
          print *, 'E23_B_term3: ', mp2f12_E23(Bijij_term3,Bjiij_term3,nocc)
          print *, 'E23_B_term4: ', mp2f12_E23(Bijij_term4,Bjiij_term4,nocc)
          print *, 'E23_B_term5: ', mp2f12_E23(Bijij_term5,Bjiij_term5,nocc)
          print *, 'E23_B_term6: ', mp2f12_E23(Bijij_term6,Bjiij_term6,nocc)
          print *, 'E23_B_term7: ', mp2f12_E23(Bijij_term7,Bjiij_term7,nocc)
          print *, 'E23_B_term8: ', mp2f12_E23(Bijij_term8,Bjiij_term8,nocc)
          print *, 'E23_B_term9: ', mp2f12_E23(Bijij_term9,Bjiij_term9,nocc)   
          print *, '----------------------------------------'
          E23_debug = mp2f12_E23(Bijij_term1,Bjiij_term1,nocc) & 
               & + mp2f12_E23(Bijij_term2,Bjiij_term2,nocc) + mp2f12_E23(Bijij_term3,Bjiij_term3,nocc) &
               & + mp2f12_E23(Bijij_term4,Bjiij_term4,nocc) + mp2f12_E23(Bijij_term5,Bjiij_term5,nocc) &
               & + mp2f12_E23(Bijij_term6,Bjiij_term6,nocc) + mp2f12_E23(Bijij_term7,Bjiij_term7,nocc) &
               & + mp2f12_E23(Bijij_term8,Bjiij_term8,nocc) + mp2f12_E23(Bijij_term9,Bjiij_term9,nocc)
          print *, 'E23_Bsum: ',  E23_debug
          !print *, 'E23_Bsum_debug: ',  mp2f12_E23(Bijij,Bjiij,nocc)
          print *, '----------------------------------------'

       endif
       
    endif
    call free_F12_mixed_MO_Matrices(HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp)

    if(DECinfo%use_canonical) then
       call mem_dealloc(Xijij)
       call mem_dealloc(Xjiij)

       if(DECinfo%F12DEBUG) then
          call mem_dealloc(Xijij_term1)
          call mem_dealloc(Xijij_term2)
          call mem_dealloc(Xijij_term3)
          call mem_dealloc(Xijij_term4)

          call mem_dealloc(Xjiij_term1)
          call mem_dealloc(Xjiij_term2)
          call mem_dealloc(Xjiij_term3)
          call mem_dealloc(Xjiij_term4)

          call mem_dealloc(Bijij_term1)
          call mem_dealloc(Bijij_term2)
          call mem_dealloc(Bijij_term3)
          call mem_dealloc(Bijij_term4)
          call mem_dealloc(Bijij_term5)
          call mem_dealloc(Bijij_term6)          
          call mem_dealloc(Bijij_term7)
          call mem_dealloc(Bijij_term8)
          call mem_dealloc(Bijij_term9)

          call mem_dealloc(Bjiij_term1)
          call mem_dealloc(Bjiij_term2)
          call mem_dealloc(Bjiij_term3)
          call mem_dealloc(Bjiij_term4)
          call mem_dealloc(Bjiij_term5)
          call mem_dealloc(Bjiij_term6)
          call mem_dealloc(Bjiij_term7)
          call mem_dealloc(Bjiij_term8)
          call mem_dealloc(Bjiij_term9)
       endif

    else
  
       call mem_dealloc(Xijkl)
       
       if(DECinfo%F12DEBUG) then
          call mem_dealloc(Xijkl_term1)
          call mem_dealloc(Xijkl_term2)
          call mem_dealloc(Xijkl_term3)
          call mem_dealloc(Xijkl_term4)   
       endif
    endif
    
    call mem_dealloc(Bijij)
    call mem_dealloc(Bjiij)

    if(DECinfo%F12DEBUG) then
       call mem_dealloc(Bijij_debug) 
       call mem_dealloc(Bjiij_debug)
    endif

    call free_4Center_F12_integrals(&
         & Ripjq,Fijkl,Tijkl,Rimjc,Dijkl,Tirjk,Tijkr,Gipjq,Gimjc,Girjs,Girjm,&
         & Grimj,Gipja,Gpiaj,Gicjm,Gcimj,Gcirj,Gciaj,Giajc)
    call free_cabs()

    if(DECinfo%F12DEBUG) then

       mp2f12_energy = 0.0E0_realk
       mp2f12_energy = mp2_energy+E21_debug+E22_debug+E23_debug

       write(*,'(1X,a,f20.10)') 'TOYCODE: MP2 CORRELATION ENERGY =           ', mp2_energy
       write(*,'(1X,a,f20.10)') 'TOYCODE: F12 E21 CORRECTION TO ENERGY =     ', E21_debug
       write(*,'(1X,a,f20.10)') 'TOYCODE: F12 E22 CORRECTION TO ENERGY =     ', E22_debug
       write(*,'(1X,a,f20.10)') 'TOYCODE: F12 E23 CORRECTION TO ENERGY =     ', E23_debug
       write(*,'(1X,a,f20.10)') 'TOYCODE: F12 E22+E23 CORRECTION TO ENERGY = ', E22_debug + E23_debug
       write(*,'(1X,a)') '-----------------------------------------------------------------'
       write(*,'(1X,a,f20.10)') 'TOYCODE: F12 CORRECTION TO ENERGY = ', E21_debug+E22_debug+E23_debug
       write(*,'(1X,a,f20.10)') 'TOYCODE: MP2-F12 ENERGY =           ', mp2_energy+E21_debug+E22_debug+E23_debug

    else

       mp2f12_energy = 0.0E0_realk
       mp2f12_energy = mp2_energy+E21+E22
       
       write(DECinfo%output,*)  'TOYCODE: MP2 CORRELATION ENERGY =        ', mp2_energy
       write(*,'(1X,a,f20.10)') 'TOYCODE: MP2 CORRELATION ENERGY =        ', mp2_energy
       write(*,'(1X,a,f20.10)') 'TOYCODE: F12 E21 CORRECTION TO ENERGY =  ', E21
       write(DECinfo%output,*)  'TOYCODE: F12 E21 CORRECTION TO ENERGY =  ', E21
       write(*,'(1X,a,f20.10)') 'TOYCODE: F12 E22 CORRECTION TO ENERGY =  ', E22
       write(DECinfo%output,*)  'TOYCODE: F12 E22 CORRECTION TO ENERGY =  ', E22
       write(*,'(1X,a,f20.10)') 'TOYCODE: F12 CORRECTION TO ENERGY =      ', E21+E22
       write(DECinfo%output,*)  'TOYCODE: F12 CORRECTION TO ENERGY =      ', E21+E22       
       ! Total MP2-F12 correlation energy
       ! Getting this energy 
       write(*,'(1X,a)') '----------------------------------------------------'
       write(*,'(1X,a,f20.10)') 'TOYCODE: MP2-F12 CORRELATION ENERGY =    ', mp2f12_energy
       write(DECinfo%output,*)  'TOYCODE: MP2-F12 CORRELATION ENERGY =    ', mp2f12_energy
    endif

    call array4_free(array4Taibj)
    call mem_dealloc(gmo)

  end subroutine full_canonical_mp2_f12
#endif

  !> \brief Memory check for full_canonical_mp2 subroutine
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine full_canonical_mp2_memory_check(nbasis,nocc,nvirt)
    implicit none
    integer,intent(in) :: nbasis
    integer,intent(in) :: nocc
    integer,intent(in) :: nvirt
    real(realk) :: MemRequired,GB

    GB = 1.0E9_realk

    ! Check for integer overflow
    if(nbasis**4 < 1) then
       call lsquit('full_canonical_mp2: Possible integer overflow! &
            & Try to compile with 64 bit integers.',-1)
    end if

    ! Check that arrays fit in memory (hard-coded)
    MemRequired = real(nbasis**4) + real(nocc*(nbasis**3)) &
         & + real(nocc*nvirt*(nbasis**2)) + real(nvirt*nocc*nvirt*nocc)
    MemRequired = MemRequired*realk/GB
    if(MemRequired > DECinfo%memory) then
       print *, 'Mem required (GB) = ', MemRequired
       print *, 'Max memory (GB)   = ', DECinfo%memory
       call lsquit('full_canonical_mp2: Memory exceeded! Try to increase memory &
            & using the .MaxMemory (in GB) keyword in the *DEC section.',-1)
    end if

  end subroutine full_canonical_mp2_memory_check


#ifdef MOD_UNRELEASED
  subroutine submp2f12_EBX(mp2f12_EBX,Bijij,Bjiij,Xijij,Xjiij,Fii,nocc)
    implicit none
    Real(realk)               :: mp2f12_EBX
    Real(realk),intent(INOUT) :: Bijij(nocc,nocc),Bjiij(nocc,nocc)
    Real(realk),intent(IN)    :: Xijij(nocc,nocc),Xjiij(nocc,nocc)
    Real(realk),intent(IN)    :: Fii(nocc,nocc)
    Integer,intent(IN)        :: nocc
    !
    Integer     :: i,j,k
    Real(realk) :: tmp

    DO j=1,nocc
       DO i=1,nocc
          Bijij(i,j) = Bijij(i,j)-(Fii(i,i)+Fii(j,j))*Xijij(i,j)
          Bjiij(i,j) = Bjiij(i,j)-(Fii(i,i)+Fii(j,j))*Xjiij(i,j)
       ENDDO
    ENDDO

    tmp = 0E0_realk
    DO i=1,nocc
       tmp = tmp + Bijij(i,i)
    ENDDO
    mp2f12_EBX = 0.25E0_realk*tmp

    tmp = 0E0_realk
    DO j=1,nocc
       DO i=j+1,nocc
          tmp = tmp + 7E0_realk * Bijij(i,j) + Bjiij(i,j)
       ENDDO
    ENDDO
    mp2f12_EBX = mp2f12_EBX + tmp/16E0_realk
  end subroutine submp2f12_EBX

 !> Function for finding the E22X energy (non-canonical)
  function mp2f12_E22X(Xijkl,Fii,nocc) result(energy)
    implicit none
    integer,intent(IN)  :: nocc
    real(realk), pointer :: Bijij(:,:), Bjiij(:,:)
    !
    Real(realk),intent(IN) :: Xijkl(nocc,nocc,nocc,nocc)
    real(realk),intent(IN) :: Fii(nocc,nocc)
    real(realk) :: energy
    !
    integer     :: i,j,k
    real(realk) :: tmp

    call mem_alloc(Bijij,nocc,nocc)
    call mem_alloc(Bjiij,nocc,nocc)

    Bijij = 0.0E0_realk
    Bjiij = 0.0E0_realk

    DO j=1,nocc
       DO i=1,nocc
          DO k=1,nocc
             Bijij(i,j) = Bijij(i,j)-(Fii(k,i)*Xijkl(i,k,j,j)+Fii(k,j)*Xijkl(i,i,j,k))
             Bjiij(i,j) = Bjiij(i,j)-(Fii(k,i)*Xijkl(i,j,j,k)+Fii(k,j)*Xijkl(i,j,k,i))
          ENDDO
       ENDDO
    ENDDO

    energy = 0E0_realk
    tmp = 0E0_realk
    DO i=1,nocc
       tmp = tmp + Bijij(i,i)
    ENDDO
    energy = 0.25E0_realk*tmp

    tmp = 0E0_realk
    DO j=1,nocc
       DO i=j+1,nocc
          tmp = tmp + 7E0_realk * Bijij(i,j) + Bjiij(i,j)
       ENDDO
    ENDDO
    energy = energy + tmp/16E0_realk    
    
    call mem_dealloc(Bijij)
    call mem_dealloc(Bjiij)

  end function mp2f12_E22X

  subroutine submp2f12_EBXfull(mp2f12_EBX,Bijij,Bjiij,Xijkl,Fii,nocc)
    implicit none
    Real(realk)               :: mp2f12_EBX
    Real(realk),intent(INOUT) :: Bijij(nocc,nocc),Bjiij(nocc,nocc)
    Real(realk),intent(IN)    :: Xijkl(nocc,nocc,nocc,nocc)
    Real(realk),intent(IN)    :: Fii(nocc,nocc)
    Integer,intent(IN)        :: nocc
    !
    Integer     :: i,j,k
    Real(realk) :: tmp

    DO j=1,nocc
       DO i=1,nocc
          DO k=1,nocc
             Bijij(i,j) = Bijij(i,j)-(Fii(k,i)*Xijkl(i,k,j,j)+Fii(k,j)*Xijkl(i,i,j,k))
             Bjiij(i,j) = Bjiij(i,j)-(Fii(k,i)*Xijkl(i,j,j,k)+Fii(k,j)*Xijkl(i,j,k,i))
          ENDDO
       ENDDO
    ENDDO

    tmp = 0E0_realk
    DO i=1,nocc
       tmp = tmp + Bijij(i,i)
    ENDDO
    mp2f12_EBX = 0.25E0_realk*tmp

    tmp = 0E0_realk
    DO j=1,nocc
       DO i=j+1,nocc
          tmp = tmp + 7E0_realk * Bijij(i,j) + Bjiij(i,j)
       ENDDO
    ENDDO
    mp2f12_EBX = mp2f12_EBX + tmp/16E0_realk
  end subroutine submp2f12_EBXfull

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

  !> Function for finding the E22 energy (canonical)
  function mp2f12_E22(Xijij,Xjiij,Fii,nocc) result(energy)
    implicit none
    integer,intent(IN)  :: nocc
    real(realk), pointer :: Bijij(:,:), Bjiij(:,:)
    !
    real(realk),intent(IN) :: Xijij(nocc,nocc), Xjiij(nocc,nocc)
    real(realk),intent(IN) :: Fii(nocc,nocc)
    real(realk) :: energy
    !
    integer     :: i,j
    real(realk) :: tmp

    call mem_alloc(Bijij,nocc,nocc)
    call mem_alloc(Bjiij,nocc,nocc)

    Bijij = 0.0E0_realk
    Bjiij = 0.0E0_realk

    DO j=1,nocc
       DO i=1,nocc
          Bijij(i,j) = -1.0E0_realk*(Fii(i,i)+Fii(j,j))*Xijij(i,j)
          Bjiij(i,j) = -1.0E0_realk*(Fii(i,i)+Fii(j,j))*Xjiij(i,j)
       ENDDO
    ENDDO

    tmp = 0E0_realk
    DO i=1,nocc
       tmp = tmp + Bijij(i,i)
    ENDDO

    energy = 0.25E0_realk*tmp
    tmp = 0E0_realk

    DO j=1,nocc
       DO i=j+1,nocc
          tmp = tmp + 7.0E0_realk * Bijij(i,j) + Bjiij(i,j)
       ENDDO
    ENDDO
    energy = energy + 0.0625E0_realk*tmp

    call mem_dealloc(Bijij)
    call mem_dealloc(Bjiij)

  end function mp2f12_E22

  !> Function for finding the E23 energy for the B-matrix
  function mp2f12_E23(Bijij,Bjiij,nocc) result(energy)
    implicit none
    integer,intent(IN)  :: nocc
    !>
    real(realk),intent(IN) :: Bijij(nocc,nocc), Bjiij(nocc,nocc)
    real(realk) :: energy
    !
    integer     :: i,j
    real(realk) :: tmp

    tmp = 0E0_realk
    DO i=1,nocc
       tmp = tmp + Bijij(i,i)
    ENDDO

    energy = 0.25E0_realk*tmp
    tmp = 0E0_realk

    DO j=1,nocc
       DO i=j+1,nocc
          tmp = tmp + 7.0E0_realk * Bijij(i,j) + Bjiij(i,j)
       ENDDO
    ENDDO
    energy = energy + 0.0625E0_realk*tmp !1/16
  end function mp2f12_E23

#endif

  subroutine get_4Center_MO_integrals(mylsitem,lupri,nbasis,nocc,noccfull,nvirt,&
       & Cocc,Cvirt,inputstring,gAO,gMO)
    implicit none
    integer :: nocc,noccfull,nvirt,nCabsAO,nCabs,nbasis
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    character(len=4) :: inputstring
    integer :: ndim2(4),ndim1(4)
    real(realk),pointer :: gAO(:,:,:,:)
    real(realk),pointer :: gMO(:,:,:,:) ,elms(:)
    type(matrix) :: CMO(4)
    real(realk),dimension(nbasis,nocc),intent(in) :: Cocc
    !> Virtual MO coefficients
    real(realk),dimension(nbasis,nvirt),intent(in) :: Cvirt
    type(matrix) :: CMO_cabs,CMO_ri
    real(realk),pointer :: tmp(:,:,:,:)
    real(realk),pointer :: tmp2(:,:,:,:)
    character :: string(4)
    logical :: doCABS,doRI
    integer :: i,lupri
    string(1) = inputstring(1:1)
    string(2) = inputstring(2:2)
    string(3) = inputstring(3:3)
    string(4) = inputstring(4:4)
    doCABS = .FALSE.
    do i=1,4
       if(string(i).EQ.'c')then !occupied active
          doCABS = .TRUE.
       endif
    enddo
    doRI = .FALSE.
    do i=1,4
       if(string(i).EQ.'r')then !occupied active
          doRI = .TRUE.
       endif
    enddo
    call determine_CABS_nbast(nCabsAO,nCabs,mylsitem%SETTING,lupri)
    IF(doCABS)THEN
       call mat_init(CMO_cabs,nCabsAO,nCabs)
       call build_CABS_MO(CMO_cabs,nCabsAO,mylsitem%SETTING,lupri)
    ENDIF
    IF(doRI)THEN
       call mat_init(CMO_ri,nCabsAO,nCabsAO)
       call build_RI_MO(CMO_ri,nCabsAO,mylsitem%SETTING,lupri)
    ENDIF
    do i=1,4
       if(string(i).EQ.'i')then !occupied active
          ndim1(i) = nbasis
          ndim2(i) = nocc
       elseif(string(i).EQ.'m')then !all occupied
          ndim1(i) = nbasis
          ndim2(i) = noccfull
       elseif(string(i).EQ.'p')then !all occupied + virtual
          ndim1(i) = nbasis
          ndim2(i) = nbasis
       elseif(string(i).EQ.'a')then !virtual
          ndim1(i) = nbasis
          ndim2(i) = nvirt
       elseif(string(i).EQ.'c')then !cabs
          ndim1(i) = ncabsAO
          ndim2(i) = ncabs
       elseif(string(i).EQ.'r')then !ri - MOs
          ndim1(i) = ncabsAO
          ndim2(i) = ncabsAO
       endif
       call mat_init(CMO(i),ndim1(i),ndim2(i))
       if(string(i).EQ.'i')then !occupied active
          call dcopy(ndim2(i)*ndim1(i),Cocc,1,CMO(i)%elms,1)
       elseif(string(i).EQ.'m')then !all occupied
          call dcopy(ndim2(i)*ndim1(i),Cocc,1,CMO(i)%elms,1)
       elseif(string(i).EQ.'p')then !all occupied + virtual
          call dcopy(noccfull*nbasis,Cocc,1,CMO(i)%elms,1)
          call dcopy(nvirt*nbasis,Cvirt,1,CMO(i)%elms(noccfull*nbasis+1),1)
       elseif(string(i).EQ.'a')then !virtual
          call dcopy(ndim2(i)*ndim1(i),Cvirt,1,CMO(i)%elms,1)
       elseif(string(i).EQ.'c')then !cabs
          call dcopy(ndim2(i)*ndim1(i),CMO_cabs%elms,1,CMO(i)%elms,1)
       elseif(string(i).EQ.'r')then !ri - MOs
          call dcopy(ndim2(i)*ndim1(i),CMO_RI%elms,1,CMO(i)%elms,1)
       endif
    enddo
    IF(doCABS)THEN
       call mat_free(CMO_cabs)
    ENDIF
    IF(doRI)THEN
       call mat_free(CMO_ri)
    ENDIF
    call mem_alloc(tmp,ndim2(1),ndim1(2),ndim1(3),ndim1(4))
    call ls_dzero(tmp,ndim2(1)*ndim1(2)*ndim1(3)*ndim1(4))
    call sub1(gao,tmp,CMO(1)%elms,ndim2,ndim1)

    call mem_alloc(tmp2,ndim2(1),ndim2(2),ndim1(3),ndim1(4))
    call ls_dzero(tmp2,ndim2(1)*ndim2(2)*ndim1(3)*ndim1(4))
    call sub2(tmp,tmp2,CMO(2)%elms,ndim2,ndim1)
    call mem_dealloc(tmp)

    call mem_alloc(tmp,ndim2(1),ndim2(2),ndim2(3),ndim1(4))
    call ls_dzero(tmp,ndim2(1)*ndim2(2)*ndim2(3)*ndim1(4))
    call sub3(tmp2,tmp,CMO(3)%elms,ndim2,ndim1)
    call mem_dealloc(tmp2)

    call mem_alloc(gMO,ndim2(1),ndim2(2),ndim2(3),ndim2(4))
    call ls_dzero(gMO,ndim2(1)*ndim2(2)*ndim2(3)*ndim2(4))
    call sub4(tmp,gMO,CMO(4)%elms,ndim2,ndim1)
    call mem_dealloc(tmp)
    do i=1,4
       call mat_free(CMO(i))
    enddo

  contains

    subroutine sub1(gao,tmp,elms,ndim2,ndim1)
      implicit none
      integer :: ndim1(4),ndim2(4)
      real(realk),intent(in) :: elms(ndim1(1),ndim2(1))
      real(realk),intent(in) :: gao(ndim1(1),ndim1(2),ndim1(3),ndim1(4))
      real(realk) :: tmp(ndim2(1),ndim1(2),ndim1(3),ndim1(4))
      integer :: a,b,c,d,p

      do d=1,ndim1(4)
         do c=1,ndim1(3)
            do b=1,ndim1(2)
               do a=1,ndim2(1)

                  do p=1,ndim1(1)
                     tmp(a,b,c,d) = tmp(a,b,c,d) + gao(p,b,c,d)*elms(p,a)
                  end do

               end do
            end do
         end do
      end do
    end subroutine sub1

    subroutine sub2(tmp,tmp2,elms,ndim2,ndim1)
      implicit none
      integer :: ndim1(4),ndim2(4)
      real(realk),intent(in) :: elms(ndim1(2),ndim2(2))
      real(realk),intent(in) :: tmp(ndim2(1),ndim1(2),ndim1(3),ndim1(4))
      real(realk) :: tmp2(ndim2(1),ndim2(2),ndim1(3),ndim1(4))
      integer :: a,b,c,d,p

      do d=1,ndim1(4)
         do c=1,ndim1(3)
            do b=1,ndim2(2)
               do a=1,ndim2(1)

                  do p=1,ndim1(2)
                     tmp2(a,b,c,d) = tmp2(a,b,c,d) + tmp(a,p,c,d)*elms(p,b)
                  end do

               end do
            end do
         end do
      end do
    end subroutine sub2

    subroutine sub3(tmp2,tmp,elms,ndim2,ndim1)
      implicit none
      integer :: ndim1(4),ndim2(4)
      real(realk),intent(in) :: elms(ndim1(3),ndim2(3))
      real(realk),intent(in) :: tmp2(ndim2(1),ndim2(2),ndim1(3),ndim1(4))
      real(realk) :: tmp(ndim2(1),ndim2(2),ndim2(3),ndim1(4))
      integer :: a,b,c,d,p

      do d=1,ndim1(4)
         do c=1,ndim2(3)
            do b=1,ndim2(2)
               do a=1,ndim2(1)

                  do p=1,ndim1(3)
                     tmp(a,b,c,d) = tmp(a,b,c,d) + tmp2(a,b,p,d)*elms(p,c)
                  end do

               end do
            end do
         end do
      end do
    end subroutine sub3

    subroutine sub4(tmp,gMO,elms,ndim2,ndim1)
      implicit none
      integer :: ndim1(4),ndim2(4)
      real(realk),intent(in) :: elms(ndim1(4),ndim2(4))
      real(realk),intent(in) :: tmp(ndim2(1),ndim2(2),ndim2(3),ndim1(4))
      real(realk) :: gMO(ndim2(1),ndim2(2),ndim2(3),ndim2(4))
      integer :: a,b,c,d,p

      do d=1,ndim2(4)
         do c=1,ndim2(3)
            do b=1,ndim2(2)
               do a=1,ndim2(1)

                  do p=1,ndim1(4)
                     gMO(a,b,c,d) = gMO(a,b,c,d) + tmp(a,b,c,p)*elms(p,d)
                  end do

               end do
            end do
         end do
      end do
    end subroutine sub4

  end subroutine get_4Center_MO_integrals

#ifdef MOD_UNRELEASED
  !> \brief Get CCSD-F12 energy, testing code.
  !> \date May 2012
  subroutine full_get_ccsd_f12_energy(MyMolecule,MyLsitem,Dmat,ECCSD_F12)

    implicit none

    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> HF density matrix
    type(matrix),intent(in) :: Dmat
    !> CCSD-F12 correlation energy
    real(realk), intent(inout) :: ECCSD_F12
    real(realk),pointer :: gao(:,:,:,:)
    real(realk),pointer :: AIBJ(:,:,:,:)
    real(realk) :: ECCSD, Ef12, Ttmp, Gtmp,E21,E22
    integer :: ncabs, ncabsAO, nocc,nvirt,nbasis,noccfull, i,j,a,b
    type(array2) :: Tai
    type(array4) :: Taibj
    real(realk),pointer :: Ripjq(:,:,:,:)
    real(realk),pointer :: Fijkl(:,:,:,:)
    real(realk),pointer :: Tijkl(:,:,:,:)
    real(realk),pointer :: Rimjc(:,:,:,:)
    real(realk),pointer :: Dijkl(:,:,:,:)
    real(realk),pointer :: Tirjk(:,:,:,:)
    real(realk),pointer :: Tijkr(:,:,:,:)
    real(realk),pointer :: Gipjq(:,:,:,:)
    real(realk),pointer :: Gimjc(:,:,:,:)
    real(realk),pointer :: Girjs(:,:,:,:)
    real(realk),pointer :: Girjm(:,:,:,:)
    real(realk),pointer :: Grimj(:,:,:,:)
    real(realk),pointer :: Gipja(:,:,:,:)
    real(realk),pointer :: Gpiaj(:,:,:,:)
    real(realk),pointer :: Gicjm(:,:,:,:)
    real(realk),pointer :: Gcimj(:,:,:,:)
    real(realk),pointer :: Gcirj(:,:,:,:)
    real(realk),pointer :: Gciaj(:,:,:,:)
    real(realk),pointer :: Giajc(:,:,:,:)
    real(realk),pointer :: Ciajb(:,:,:,:)
    real(realk),pointer :: Cjaib(:,:,:,:)
    real(realk),pointer :: Tiajb(:,:,:,:)
    real(realk),pointer :: Tjaib(:,:,:,:)
    type(matrix) :: K
    type(matrix) :: HJir
    type(matrix) :: Krr
    type(matrix) :: Frr
    type(matrix) :: Frc
    type(matrix) :: Fpp
    type(matrix) :: Fmm
    type(matrix) :: Frm
    type(matrix) :: Fcp
    type(matrix) :: Fii
    type(matrix) :: Fac

    !   F12 specific
    real(realk),pointer :: Vijij(:,:)
    real(realk),pointer :: Vjiij(:,:)    
    real(realk),pointer :: Xijkl(:,:,:,:)
    real(realk),pointer :: Xijij(:,:)
    real(realk),pointer :: Xjiij(:,:)
    real(realk),pointer :: Bijij(:,:)
    real(realk),pointer :: Bjiij(:,:)

    !   CCSD specific
    real(realk),pointer :: Rapbq(:,:,:,:)
    real(realk),pointer :: Rambc(:,:,:,:)
    real(realk),pointer :: Fiajb(:,:,:,:)
    real(realk),pointer :: Viajb(:,:,:,:)
    real(realk),pointer :: Ripaq(:,:,:,:)
    real(realk),pointer :: Rimac(:,:,:,:)
    real(realk),pointer :: Ramic(:,:,:,:)
    real(realk),pointer :: Fijka(:,:,:,:)
    real(realk),pointer :: Viija(:,:,:)
    real(realk),pointer :: Vijja(:,:,:)
    real(realk),pointer :: Viaji(:,:,:)
    real(realk),pointer :: Viajj(:,:,:)
    !    logical :: fulldriver 
    !    fulldriver = .TRUE.
    !    call init_cabs(fulldriver)

    ! Init dimensions
    nocc = MyMolecule%nocc
    nvirt = MyMolecule%nunocc
    nbasis = MyMolecule%nbasis
    noccfull = nocc
    call determine_CABS_nbast(ncabsAO,ncabs,mylsitem%setting,DECinfo%output)

    ! Get full CCSD singles (Tai) and doubles (Taibj) amplitudes
    call full_get_ccsd_singles_and_doubles(MyMolecule,MyLsitem,Tai,Taibj)

    ! Get all MO mixed matrices
    call get_F12_mixed_MO_Matrices(MyLsitem,MyMolecule,Dmat,nbasis,ncabsAO,&
         & nocc,noccfull,nvirt,ncabs,HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp)

    ! Get all AO integrals in regular basis
    call mem_alloc(gao,nbasis,nbasis,nbasis,nbasis)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRC')

    ! Transform AO integrals to MO integrals (A I | B J)
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'aiai',gAO,AIBJ)
    call mem_dealloc(gao)

    call get_4Center_F12_integrals(mylsitem,MyMolecule,nbasis,nocc,noccfull,nvirt,ncabsAO,&
         & Ripjq,Fijkl,Tijkl,Rimjc,Dijkl,Tirjk,Tijkr,Gipjq,Gimjc,Girjs,Girjm,&
         & Grimj,Gipja,Gpiaj,Gicjm,Gcimj,Gcirj,Gciaj,Giajc)

    !   Rrrrr
    call mem_alloc(gao,nbasis,nbasis,nbasis,nbasis)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRC')
    !   Rapbq
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'apap',gAO,Rapbq)
    !   Ripaq
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'ipap',gAO,Ripaq)
    call mem_dealloc(gao)

    !   Rrrrc
    call mem_alloc(gao,nbasis,nbasis,nbasis,ncabsAO)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRCC')
    !   Rambc
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'amac',gAO,Rambc)
    !   Rimac
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'imac',gAO,Rimac)
    !   Ramic
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'amic',gAO,Ramic)
    call mem_dealloc(gao)

    !   Rrrrc
    call mem_alloc(gao,nbasis,nbasis,nbasis,ncabsAO)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRF')
    !   Fiajb
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'iaia',gAO,Fiajb)
    !   Fijka
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co, MyMolecule%Cv,'iiia',gAO,Fijka)
    call mem_dealloc(gao)

    ! Calculate standard CCSD energy (brainless summation in this test code)
    ECCSD=0.0E0_realk
    do j=1,nocc
       do b=1,nvirt
          do i=1,nocc
             do a=1,nvirt
                ! Energy = sum_{ijab} ( Tai*Tbj + Taibj) * (ai | bj)
                Ttmp = Tai%val(a,i)*Tai%val(b,j) + Taibj%val(a,i,b,j)
                Gtmp = 2.0E0_realk * AIBJ(a,i,b,j) - AIBJ(b,i,a,j)
                ECCSD = ECCSD + Ttmp * Gtmp
             end do
          end do
       end do
    end do

    print *, 'TOYCODE: CCSD CORRELATION ENERGY = ', ECCSD
    write(DECinfo%output,*) 'TOYCODE: CCSD CORRELATION ENERGY = ', ECCSD


    ! INTEGRAL GUYS: DO YOUR MAGIC FOR THE F12 CONTRIBUTION HERE...
    call mem_alloc(Vijij,nocc,nocc)
    call mem_alloc(Vjiij,nocc,nocc)
    call mp2f12_Vijij(Vijij,Ripjq,Gipjq,Fijkl,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)
    call mp2f12_Vjiij(Vjiij,Ripjq,Gipjq,Fijkl,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)

    call mem_alloc(Ciajb,nocc,nvirt,nocc,nvirt)
    call mem_alloc(Cjaib,nocc,nvirt,nocc,nvirt)
    call mp2f12_Ciajb(Ciajb,Giajc,Fac%elms,nocc,nvirt,ncabs)


    call mem_alloc(Viajb,nocc,nvirt,nocc,nvirt)
    call ccsdf12_Viajb(Viajb,Rapbq,Gipjq,Fiajb,Rambc,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)

    call mem_alloc(Viija,nocc,nocc,nvirt)
    call ccsdf12_Viija(Viija,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)

    call mem_alloc(Vijja,nocc,nocc,nvirt)
    call ccsdf12_Vijja(Vijja,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)

    call mem_alloc(Viaji,nocc,nvirt,nocc)
    call ccsdf12_Viaji(Viaji,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)

    call mem_alloc(Viajj,nocc,nvirt,nocc)
    call ccsdf12_Viajj(Viajj,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)

    call ccsdf12_Vijij_coupling(Vijij,Ciajb,Taibj%val,Viajb,Viija,Viajj,Tai%val,nocc,nvirt)
    call ccsdf12_Vjiij_coupling(Vjiij,Ciajb,Taibj%val,Viajb,Vijja,Viaji,Tai%val,nocc,nvirt)

    ! CCSD Specific
    call mem_dealloc(Viija)
    call mem_dealloc(Vijja)
    call mem_dealloc(Viajj)
    call mem_dealloc(Viaji)
    call mem_dealloc(Viajb)

    E21 = 2.0E0_realk*mp2f12_E21(Vijij,Vjiij,nocc)
    print*,'E21',E21

    ! F12 Specific
    call mem_dealloc(Vijij)
    call mem_dealloc(Vjiij)
    call mem_dealloc(Ciajb)
    call mem_dealloc(Cjaib)

    if(DECinfo%use_canonical) then
       call mem_alloc(Xijij,nocc,nocc)
       call mem_alloc(Xjiij,nocc,nocc)
       call mp2f12_Xijij(Xijij,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
       call mp2f12_Xjiij(Xjiij,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
   
    else
       ! non-canonical
       call mem_alloc(Xijkl,nocc,nocc,nocc,nocc)

       call mp2f12_Xijijfull(Xijkl,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
       ! the way you build Xijij frok Xijkl
       !    do j=1,nocc
       !       do i=1,nocc
       !          Xijij(i,j) = Xijkl(i,i,j,j)
       !          Xjiij(i,j) = Xijkl(i,j,j,i)          
       !       enddo
       !    enddo
    endif

    call mem_alloc(Bijij,nocc,nocc)
    call mp2f12_Bijij(Bijij,Dijkl,Tirjk,Tijkr,Girjs,hJir%elms,Krr%elms,&
         & Frr%elms,Fpp%elms,Fmm%elms,Frm%elms,Fcp%elms,&
         & Girjm,Grimj,Gipja,Gpiaj,Gicjm,Gcimj,Gcirj,Gciaj,&
         & nocc,noccfull,nbasis,ncabsAO,nvirt,ncabs)

    call mem_alloc(Bjiij,nocc,nocc)
    call mp2f12_Bjiij(Bjiij,Dijkl,Tirjk,Tijkr,Girjs,hJir%elms,Krr%elms,&
         & Frr%elms,Fpp%elms,Fmm%elms,Frm%elms,Fcp%elms,&
         & Girjm,Grimj,Gipja,Gpiaj,Gicjm,Gcimj,Gcirj,Gciaj,&
         & nocc,noccfull,nbasis,ncabsAO,nvirt,ncabs)

    if(DECinfo%use_canonical) then
       call submp2f12_EBX(E22,Bijij,Bjiij,Xijij,Xjiij,Fii%elms,nocc)

    else
       call submp2f12_EBXfull(E22,Bijij,Bjiij,Xijkl,Fii%elms,nocc)
    endif

    call free_F12_mixed_MO_Matrices(HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp)
    call mem_dealloc(Rapbq)
    call mem_dealloc(Ripaq)
    call mem_dealloc(Rambc)
    call mem_dealloc(Rimac)
    call mem_dealloc(Ramic)
    call mem_dealloc(Fiajb)
    call mem_dealloc(Fijka)
    if(DECinfo%use_canonical) then
       call mem_dealloc(Xijij)
       call mem_dealloc(Xjiij)
    else
       call mem_dealloc(Xijkl)
    endif
    call mem_dealloc(Bijij)
    call mem_dealloc(Bjiij)

    EF12 = E21 + E22

    ! Add contributions
    ECCSD_F12 = ECCSD + EF12
    print *, 'TOYCODE: CCSD-F12 CORRECTION TO ENERGY  = ', EF12
    print *, 'TOYCODE: CCSD-F12 CORRELATION ENERGY = ', ECCSD_F12
    write(DECinfo%output,*) 'TOYCODE: CCSD-F12 CORRECTION TO ENERGY  = ', EF12
    write(DECinfo%output,*) 'TOYCODE: CCSD-F12 CORRELATION ENERGY = ', ECCSD_F12

    call mem_dealloc(AIBJ)
    call array4_free(Taibj)
    call array2_free(Tai)
    call free_4Center_F12_integrals(&
         & Ripjq,Fijkl,Tijkl,Rimjc,Dijkl,Tirjk,Tijkr,Gipjq,Gimjc,Girjs,Girjm,&
         & Grimj,Gipja,Gpiaj,Gicjm,Gcimj,Gcirj,Gciaj,Giajc)
    call free_cabs()

  end subroutine full_get_ccsd_f12_energy
#endif

  !> \brief Get CCSD singles and doubles amplitude for full molecule,
  !> only to be used for debugging purposes.
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine full_get_ccsd_singles_and_doubles(MyMolecule,MyLsitem,Tai, Taibj)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Singles amplitudes
    type(array2),intent(inout) :: Tai
    !> Doubles amplitudes
    type(array4),intent(inout) :: Taibj
    integer :: nocc,nunocc,nbasis,print_level,save_model,startidx,endidx,i,j
    logical :: fragment_job
    real(realk) :: energy
    type(array4) :: VOVO
    real(realk),pointer :: ppfock(:,:)
    logical :: local
    local = .true.
#ifdef VAR_MPI
    local = .false.
#endif


    ! Quick fix to always use CCSD model
    !    save_model=DECinfo%ccmodel
    !    DECinfo%ccmodel=3

    if(DECinfo%FrozenCore) then
       nocc = MyMolecule%nval
    else
       nocc = MyMolecule%nocc
    end if
    nunocc = MyMolecule%nunocc
    nbasis = MyMolecule%nbasis

    fragment_job = .false.
    print_level = 0 ! this is not used


    if(DECinfo%frozencore) then
       ! Pass only valence orbitals
       call mem_alloc(ppfock,nocc,nocc)
       do j=1,nocc
          do i=1,nocc
             ppfock(i,j) = MyMolecule%ppfock(MyMolecule%ncore+i,MyMolecule%ncore+j)
          end do
       end do

       startidx = MyMolecule%ncore+1  
       endidx = MyMolecule%nocc
       call ccsolver_par(DECinfo%ccmodel,MyMolecule%Co(1:nbasis,startidx:endidx),&
            & MyMolecule%Cv,MyMolecule%fock, nbasis,nocc,nunocc,mylsitem,&
            & print_level,fragment_job,&
            & ppfock,MyMolecule%qqfock,energy, Tai, Taibj, VOVO,.false.,local)
       call mem_dealloc(ppfock)

    else

       call ccsolver_par(DECinfo%ccmodel,MyMolecule%Co,MyMolecule%Cv,&
            & MyMolecule%fock, nbasis,nocc,nunocc,mylsitem, print_level, fragment_job,&
            & MyMolecule%ppfock,MyMolecule%qqfock,&
            & energy, Tai, Taibj, VOVO,.false.,local)
    end if

    call array4_free(VOVO)

  end subroutine full_get_ccsd_singles_and_doubles


#ifdef MOD_UNRELEASED
  subroutine get_4Center_F12_integrals(mylsitem,MyMolecule,nbasis,nocc,noccfull,nvirt,ncabsAO,&
       & Ripjq,Fijkl,Tijkl,Rimjc,Dijkl,Tirjk,Tijkr,Gipjq,Gimjc,Girjs,Girjm,&
       & Grimj,Gipja,Gpiaj,Gicjm,Gcimj,Gcirj,Gciaj,Giajc)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    integer :: nbasis,nocc,nvirt,noccfull,ncabsAO
    real(realk),pointer :: Ripjq(:,:,:,:)
    real(realk),pointer :: Fijkl(:,:,:,:)
    real(realk),pointer :: Tijkl(:,:,:,:)
    real(realk),pointer :: Rimjc(:,:,:,:)
    real(realk),pointer :: Dijkl(:,:,:,:)
    real(realk),pointer :: Tirjk(:,:,:,:)
    real(realk),pointer :: Tijkr(:,:,:,:)
    real(realk),pointer :: Gipjq(:,:,:,:)
    real(realk),pointer :: Gimjc(:,:,:,:)
    real(realk),pointer :: Girjs(:,:,:,:)
    real(realk),pointer :: Girjm(:,:,:,:)
    real(realk),pointer :: Grimj(:,:,:,:)
    real(realk),pointer :: Gipja(:,:,:,:)
    real(realk),pointer :: Gpiaj(:,:,:,:)
    real(realk),pointer :: Gicjm(:,:,:,:)
    real(realk),pointer :: Gcimj(:,:,:,:)
    real(realk),pointer :: Gcirj(:,:,:,:)
    real(realk),pointer :: Gciaj(:,:,:,:)
    real(realk),pointer :: Giajc(:,:,:,:)
    !
    real(realk),pointer :: gao(:,:,:,:)

    call mem_alloc(gao,nbasis,nbasis,nbasis,nbasis)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRC')

    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co, MyMolecule%Cv,'ipip',gAO,Ripjq)

    !Calculate the various Gaussian geminal integrals with four regular AO indeces
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co, MyMolecule%Cv,'ipip',gAO,Gipjq)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRF')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co, MyMolecule%Cv,'iiii',gAO,Fijkl)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRD')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                        MyMolecule%Co, MyMolecule%Cv,'iiii',gAO,Dijkl)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRR2')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co, MyMolecule%Cv,'iiii',gAO,Tijkl)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%Co, MyMolecule%Cv,'ipia',gAO,Gipja)


    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%Co, MyMolecule%Cv,'piai',gAO,Gpiaj)


    call mem_dealloc(gao)

    !Calculate the various Gaussian geminal integrals with RRRC
    call mem_alloc(gao,nbasis,nbasis,nbasis,ncabsAO)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRC2')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co, MyMolecule%Cv,'iiir',gAO,Tijkr)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRCC')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co, MyMolecule%Cv,'imic',gAO,Rimjc)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRCG')

    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co, MyMolecule%Cv,'imic',gAO,Gimjc)
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co, MyMolecule%Cv,'iaic',gAO,Giajc)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRC2')

    call mem_dealloc(gao)
    !Calculate the various Gaussian geminal integrals with RCRR
    call mem_alloc(gao,nbasis,ncabsAO,nbasis,nbasis)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RCRR2')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co, MyMolecule%Cv,'irii',gAO,Tirjk)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RCRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%Co, MyMolecule%Cv,'irim',gAO,Girjm)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RCRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%Co, MyMolecule%Cv,'icim',gAO,Gicjm)

    call mem_dealloc(gao)
    !Calculate the various Gaussian geminal integrals with RCRC
    call mem_alloc(gao,nbasis,ncabsAO,nbasis,ncabsAO)

    !Calculate the various Gaussian geminal integrals with four regular AO indeces
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RCRCG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                          MyMolecule%Co, MyMolecule%Cv,'irir',gAO,Girjs)

    call mem_dealloc(gao)
    !Calculate the various Gaussian geminal integrals with CRRR
    call mem_alloc(gao,ncabsAO,nbasis,nbasis,nbasis)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'CRRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%Co, MyMolecule%Cv,'rimi',gAO,Grimj)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'CRRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%Co, MyMolecule%Cv,'cimi',gAO,Gcimj)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'CRRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%Co, MyMolecule%Cv,'ciai',gAO,Gciaj)

    call mem_dealloc(gao)
    !Calculate the various Gaussian geminal integrals with CRCR
    call mem_alloc(gao,ncabsAO,nbasis,ncabsAO,nbasis)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'CRCRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%Co, MyMolecule%Cv,'ciri',gAO,Gcirj)
    call mem_dealloc(gao)

  end subroutine get_4Center_F12_integrals

  subroutine free_4Center_F12_integrals(&
       & Ripjq,Fijkl,Tijkl,Rimjc,Dijkl,Tirjk,Tijkr,Gipjq,Gimjc,Girjs,Girjm,&
       & Grimj,Gipja,Gpiaj,Gicjm,Gcimj,Gcirj,Gciaj,Giajc)
    implicit none
    real(realk),pointer :: Ripjq(:,:,:,:)
    real(realk),pointer :: Fijkl(:,:,:,:)
    real(realk),pointer :: Tijkl(:,:,:,:)
    real(realk),pointer :: Rimjc(:,:,:,:)
    real(realk),pointer :: Dijkl(:,:,:,:)
    real(realk),pointer :: Tirjk(:,:,:,:)
    real(realk),pointer :: Tijkr(:,:,:,:)
    real(realk),pointer :: Gipjq(:,:,:,:)
    real(realk),pointer :: Gimjc(:,:,:,:)
    real(realk),pointer :: Girjs(:,:,:,:)
    real(realk),pointer :: Girjm(:,:,:,:)
    real(realk),pointer :: Grimj(:,:,:,:)
    real(realk),pointer :: Gipja(:,:,:,:)
    real(realk),pointer :: Gpiaj(:,:,:,:)
    real(realk),pointer :: Gicjm(:,:,:,:)
    real(realk),pointer :: Gcimj(:,:,:,:)
    real(realk),pointer :: Gcirj(:,:,:,:)
    real(realk),pointer :: Gciaj(:,:,:,:)
    real(realk),pointer :: Giajc(:,:,:,:)
    call mem_dealloc(Ripjq)
    call mem_dealloc(Fijkl)
    call mem_dealloc(Tijkl)
    call mem_dealloc(Rimjc)
    call mem_dealloc(Dijkl)
    call mem_dealloc(Tirjk)
    call mem_dealloc(Tijkr)
    call mem_dealloc(Gipjq)
    call mem_dealloc(Gimjc)
    call mem_dealloc(Girjs)
    call mem_dealloc(Girjm)
    call mem_dealloc(Grimj)
    call mem_dealloc(Gipja)
    call mem_dealloc(Gpiaj)
    call mem_dealloc(Gicjm)
    call mem_dealloc(Gcimj)
    call mem_dealloc(Gcirj)
    call mem_dealloc(Gciaj)
    call mem_dealloc(Giajc)
  end subroutine free_4Center_F12_integrals

#endif

  !> \brief Full canonical MP2 calculation, not particularly efficient, mainly to be used for
  !> testing.
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine full_canonical_mp2_correlation_energy(MyMolecule,mylsitem,Ecorr)

    implicit none
    !> Molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS Dalton info
    type(lsitem), intent(inout) :: mylsitem
    !> Canonical MP2 correlation energy
    real(realk),intent(inout) :: Ecorr
    real(realk),pointer :: Cocc(:,:), Cunocc(:,:)
    type(array4) :: g
    integer :: nbasis,i,j,a,b,ncore,offset,nocc,nunocc
    real(realk) :: eps
    real(realk), pointer :: ppfock(:,:)

    ! Sanity check
    if(.not. DECinfo%use_canonical) then
       call lsquit('full_canonical_mp2_correlation_energy requires canonical orbitals! &
            & Insert .CANONICAL keyword OR insert .PRINTFRAGS keyword to run test calculation,&
            & where the individual fragment energies are calculated',-1)
    end if

    ! Initialize stuff
    ! ****************

    if(DECinfo%frozencore) then
       ! Frozen core: Only valence orbitals
       nocc = MyMolecule%nval
    else
       nocc = MyMolecule%nocc
    end if

    nunocc = MyMolecule%nunocc
    ncore = MyMolecule%ncore
    nbasis=MyMolecule%nbasis
    call mem_alloc(ppfock,nocc,nocc)
    if(DECinfo%frozencore) then
       ! Only copy valence orbitals into array2 structure
       call mem_alloc(Cocc,nbasis,nocc)
       do i=1,nocc
          Cocc(:,i) = MyMolecule%Co(:,i+Ncore)
       end do

       ! Fock valence
       do j=1,nocc
          do i=1,nocc
             ppfock(i,j) = MyMolecule%ppfock(i+Ncore,j+Ncore)
          end do
       end do
       offset = ncore
    else
       ! No frozen core, simply copy elements for all occupied orbitals
       call mem_alloc(Cocc,nbasis,nocc)
       Cocc=MyMolecule%Co
       ppfock = MyMolecule%ppfock
       offset=0
    end if
    call mem_alloc(Cunocc,nbasis,nunocc)
    Cunocc = MyMolecule%Cv

    ! Get (AI|BJ) integrals stored in the order (A,I,B,J)
    ! ***************************************************
    call get_VOVO_integrals(mylsitem,nbasis,nocc,nunocc,Cunocc,Cocc,g)
    call mem_dealloc(Cocc)
    call mem_dealloc(Cunocc)

    ! Calculate canonical MP2 energy
    ! ******************************
    Ecorr = 0.0_realk
    do J=1,nocc
       do B=1,nunocc
          do I=1,nocc
             do A=1,nunocc
                ! Difference in orbital energies: eps(I) + eps(J) - eps(A) - eps(B)
                eps = MyMolecule%ppfock(I+offset,I+offset) + MyMolecule%ppfock(J+offset,J+offset) &
                     & - MyMolecule%qqfock(A,A) - MyMolecule%qqfock(B,B)

                ! Ecorr = sum_{IJAB} (AI|BJ) * [ 2*(AI|BJ) - (BI|AJ) ] / [eps(I)+eps(J)-eps(A)-eps(B)]
                Ecorr = Ecorr + g%val(A,I,B,J)*(2E0_realk*g%val(A,I,B,J)-g%val(B,I,A,J))/eps
             enddo
          enddo
       enddo
    enddo

    call mem_dealloc(ppfock)
    call array4_free(g)

  end subroutine Full_canonical_mp2_correlation_energy


end module full
