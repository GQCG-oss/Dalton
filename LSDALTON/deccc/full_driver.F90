!> @file
!> Full calculation
!> This file is mainly a playground for new developments, not intended to be included in a release.

module full 

#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
  use decmpi_module, only: mpi_bcast_fullmolecule
  use lsmpi_op
#endif
  use fundamental
  use precision
  use typedeftype!,only:lsitem
  use typedef
  use dec_typedef_module
  use matrix_module
  use matrix_operations
  use memory_handling
  use MemoryLeakToolMod
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
!  use fragment_energy_module,only : Full_DECMP2_calculation

  public :: full_driver, full_canonical_rimp2
  private

contains

  !> \brief Main part for full molecular coupled-cluster calculations.
  subroutine full_driver(MyMolecule,mylsitem,D,EHF,Ecorr)

    implicit none
    !> Full molecule structure
    type(fullmolecule), intent(inout) :: MyMolecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> HF density matrix
    type(matrix),intent(in) :: D
    !> HF Energy 
    real(realk),intent(inout) :: EHF
    !> Correlation Energy 
    real(realk),intent(inout) :: Ecorr
    !local variables
    real(realk) :: Eerr

    write(DECinfo%output,'(/,a)') ' ================================================ '
    write(DECinfo%output,'(a)')   '              Full molecular driver               '
    write(DECinfo%output,'(a,/)') ' ================================================ '

#ifdef VAR_MPI
    call set_dec_settings_on_slaves()
#endif

    !MODIFY FOR NEW MODEL

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
    elseif(DECinfo%ccModel==MODEL_RIMP2)then
!       call lsquit('RIMP2 currently not implemented for **CC ',-1)
       call full_canonical_rimp2(MyMolecule,MyLsitem,Ecorr)       
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
       print *, '  E21  V terms                          '
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
       !print *, 'E21_debug: ', 2.0E0_REALK*mp2f12_E21(Vijij,Vjiij,nocc)
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

       if(DECinfo%use_canonical) then
          call mem_alloc(Xijij_term1,nocc,nocc)
          call mem_alloc(Xijij_term2,nocc,nocc)
          call mem_alloc(Xijij_term3,nocc,nocc)
          call mem_alloc(Xijij_term4,nocc,nocc)  

          call mem_alloc(Xjiij_term1,nocc,nocc)
          call mem_alloc(Xjiij_term2,nocc,nocc)
          call mem_alloc(Xjiij_term3,nocc,nocc)
          call mem_alloc(Xjiij_term4,nocc,nocc)
       endif

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

       else !> Non canonical

          call mp2f12_Xijijfull_term1(Xijkl_term1,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
          call mp2f12_Xijijfull_term2(Xijkl_term2,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
          call mp2f12_Xijijfull_term3(Xijkl_term3,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
          call mp2f12_Xijijfull_term4(Xijkl_term4,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)

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
          print *, '  E_22 X term                           '
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
          print *, 'E22_debug: ', mp2f12_E22(Xijij,Xjiij,Fii%elms,nocc)
          print *, '----------------------------------------'
          print *, '  E_23 B term                           '
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
          print *, 'E23_Bsum_debug: ',  mp2f12_E23(Bijij,Bjiij,nocc)
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

  !> \brief Memory check for full_canonical_rimp2 subroutine
  !> \author Thomas Kjaergaard
  !> \date October 2014
  subroutine full_canonical_rimp2_memory_check(nbasis,nocc,nvirt,nAux,numnodes)
    implicit none
    integer,intent(in) :: nbasis,nocc,nvirt,nAux,numnodes
    real(realk) :: MemRequired,GB

    GB = 1.0E9_realk

    ! Check that arrays fit in memory (hard-coded)
    MemRequired = 2*real(nAux*nocc*nvirt/numnodes) + real(nAux*nAux)
    MemRequired = MemRequired*realk/GB
    if(MemRequired > DECinfo%memory) then
       print *, 'RIMP2 Mem required (GB) = ', MemRequired
       print *, 'RIMP2 Max memory (GB)   = ', DECinfo%memory
       call lsquit('full_canonical_rimp2: Memory exceeded! Try to increase memory &
            & using the .MaxMemory (in GB) keyword in the *DEC section.',-1)
    end if

  end subroutine full_canonical_rimp2_memory_check

  !> \brief Calculate canonical RIMP2 energy for full molecular system
  !> \author Thomas Kjaergaard
  !> \date October 2014
  subroutine full_canonical_rimp2(MyMolecule,MyLsitem,rimp2_energy)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Canonical MP2 correlation energy
    real(realk),intent(inout) :: rimp2_energy    

    integer :: nbasis,nocc,nvirt,naux,noccfull,mynum,numnodes
    logical :: master,wakeslaves
    integer,pointer :: IPVT(:)
    real(realk), pointer   :: work1(:)
    real(realk)            :: RCOND,dummy(2)
    integer(kind=long) :: maxsize
    real(realk) :: tcpuTOT,twallTOT,tcpu_start,twall_start, tcpu_end,twall_end
    real(realk),pointer :: AlphaBeta_inv(:,:),AlphaCD(:,:,:),AlphaCD2(:,:,:),AlphaCD5(:,:,:)
    real(realk),pointer :: TMPAlphaBeta_inv(:,:),Calpha(:,:,:),Calpha2(:,:,:),AlphaCD6(:,:,:)
    real(realk),pointer :: EpsOcc(:),EpsVirt(:)
    integer(kind=ls_mpik)  :: COUNT,TAG,IERR,request,Receiver,sender,J,COUNT2,comm,TAG1,TAG2
    integer ::CurrentWait(2),nAwaitDealloc,iAwaitDealloc,I,NBA,OriginalRanknauxMPI
    integer :: myOriginalRank,node,natoms,MynauxMPI,A
    logical :: useAlphaCD5,useAlphaCD6,MessageRecieved,RoundRobin
    integer(kind=ls_mpik)  :: request1,request2,request5,request6
    integer,pointer :: nbasisauxMPI(:),startAuxMPI(:,:),AtomsMPI(:,:),nAtomsMPI(:),nAuxMPI(:,:)
    integer(KIND=long) :: MaxMemAllocated,MemAllocated
    !use Memory leak tool
    MaxMemAllocated = 0
    MemAllocated = 0
    call set_LeakTool_memvar(MaxMemAllocated,MemAllocated)
    TAG = 1023242411
    TAG1 = 1023242412
    TAG2 = 1023242413

    !sanity check
    if(.NOT.DECinfo%use_canonical) then
       call lsquit('Error: full_canonical_rimp2 require canonical Orbitals',-1)
    endif
    ! Init stuff
    ! **********
    nbasis = MyMolecule%nbasis
    nocc   = MyMolecule%nocc
    nvirt  = MyMolecule%nunocc
    naux   = MyMolecule%nauxbasis
    noccfull = nocc
    nAtoms = MyMolecule%nAtoms

#ifdef VAR_MPI
    comm = MPI_COMM_LSDALTON
    master= (infpar%mynum == infpar%master)
    mynum = infpar%mynum
    numnodes = infpar%nodtot
    wakeslaves = infpar%nodtot.GT.1
#else
    ! If MPI is not used, consider the single node to be "master"
    master=.true.
    mynum = 0
    numnodes = 1
    wakeslaves = .false.
#endif

    ! Memory check!
    ! ********************
    call full_canonical_rimp2_memory_check(nbasis,nocc,nvirt,nAux,numnodes)
    call Test_if_64bit_integer_required(nAux,nAux)

#ifdef VAR_MPI
    ! Master starts up slave
    StartUpSlaves: if(wakeslaves .and. master) then
       ! Wake up slaves to do the job: slaves awoken up with (RIMP2FULL)
       ! and call full_canonical_rimp2_slave which communicate info 
       ! then calls full_canonical_rimp2.
       call ls_mpibcast(RIMP2FULL,infpar%master,comm)
       ! Communicate fragment information to slaves
       call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,comm)
       call mpicopy_lsitem(MyLsitem,comm)
       call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,comm)
       call mpi_bcast_fullmolecule(MyMolecule)
    endif StartUpSlaves
#endif

    call mem_alloc(AlphaBeta_inv,nAux,nAux)
    call mem_LeakTool_alloc(AlphaBeta_inv,LT_RIMP2_AlphaBeta_inv)
    
    IF(master)THEN
       !=====================================================================================
       ! Major Step 1: Master Obtains Overlap (alpha|beta) in Auxiliary Basis 
       !=====================================================================================
       !This part of the Code is NOT MPI/OpenMP parallel - all nodes calculate the full overlap
       !this should naturally be changed      

       call II_get_RI_AlphaBeta_2centerInt(DECinfo%output,DECinfo%output,&
            & AlphaBeta_inv,Mylsitem%setting,nAux)

       !=====================================================================================
       ! Major Step 2: Calculate the inverse (alpha|beta)^(-1) and BCAST
       !=====================================================================================
       ! Warning the inverse is not unique so in order to make sure all slaves have the same
       ! inverse matrix we calculate it on the master a BCAST to slaves

       !Create the inverse AlphaBeta = (alpha|beta)^-1
       call mem_alloc(work1,naux)
       call mem_alloc(IPVT,naux)
       call mem_leaktool_alloc(work1,LT_TMP)
       IPVT = 0 ; RCOND = 0.0E0_realk  
       call DGECO(AlphaBeta_inv,naux,naux,IPVT,RCOND,work1)
       call DGEDI(AlphaBeta_inv,naux,naux,IPVT,dummy,work1,01)
       call mem_leaktool_dealloc(work1,LT_TMP)
       call mem_dealloc(work1)
       call mem_dealloc(IPVT)
#ifdef VAR_MPI
       call ls_mpibcast(AlphaBeta_inv,naux,naux,infpar%master,comm)
#endif
    ELSE
#ifdef VAR_MPI
       call ls_mpibcast(AlphaBeta_inv,naux,naux,infpar%master,comm)
#endif
    ENDIF

    IF(wakeslaves)then 
       !all nodes have info about all nodes 
       call mem_alloc(nbasisauxMPI,numnodes)        !number of Aux basis func assigned to rank
       call mem_alloc(nAtomsMPI,numnodes)           !atoms assign to rank
       call mem_alloc(startAuxMPI,nAtoms,numnodes)  !startindex in full (naux)
       call mem_alloc(nAuxMPI,nAtoms,numnodes)      !nauxBasis functions for each of the nAtomsMPI
       call mem_alloc(AtomsMPI,nAtoms,numnodes)     !identity of atoms in full molecule
       call getRIbasisMPI(Mylsitem%SETTING%MOLECULE(1)%p,nAtoms,numnodes,&
            & nbasisauxMPI,startAuxMPI,AtomsMPI,nAtomsMPI,nAuxMPI)
       MynauxMPI = nbasisauxmpi(mynum+1)
       call mem_dealloc(AtomsMPI) !not used in this subroutine 
    ELSE
       MynauxMPI = naux
    ENDIF

    IF(wakeslaves)then 
       IF(MynauxMPI.GT.0)THEN
          !We wish to build
          !c_(alpha,ai) = (alpha|beta)^-1 (beta|ai)
          !where alpha runs over the Aux basis functions allocated for this rank
          !and beta run over the full set of naux
          call mem_alloc(TMPAlphaBeta_inv,MynauxMPI,naux)
          call mem_LeakTool_alloc(TMPAlphaBeta_inv,LT_RIMP2_TMPAlphaBeta_inv)
          call RIMP2_buildTMPAlphaBeta_inv(TMPAlphaBeta_inv,MynauxMPI,naux,&
               & nAtomsMPI,mynum,startAuxMPI,nAuxMPI,AlphaBeta_inv,numnodes,nAtoms)
          call mem_LeakTool_dealloc(AlphaBeta_inv,LT_RIMP2_AlphaBeta_inv)
          call mem_dealloc(AlphaBeta_inv)
       ELSE
          call mem_LeakTool_dealloc(AlphaBeta_inv,LT_RIMP2_AlphaBeta_inv)
          call mem_dealloc(AlphaBeta_inv)          
       ENDIF
    ENDIF

    IF(MynauxMPI.GT.0)THEN
       !call mem_alloc(AlphaCD,naux,nvirt,nocc)
       !It is very annoying but I allocated AlphaCD inside 
       !II_get_RI_AlphaCD_3centerInt2 due to memory concerns
       !This Part of the Code is MPI/OpenMP parallel and AlphaCD will have the dimensions
       !(MynauxMPI,nvirt,nocc) 
       !nauxMPI is naux divided out on the nodes so roughly nauxMPI=naux/numnodes
       call II_get_RI_AlphaCD_3centerInt2(DECinfo%output,DECinfo%output,&
            & AlphaCD,mylsitem%setting,naux,nbasis,&
            & nvirt,nocc,MyMolecule%Cv,MyMolecule%Co,maxsize,mynum,numnodes)
       call mem_LeakTool_alloc(AlphaCD,LT_AlphaCD)
    ENDIF

    IF(wakeslaves)THEN !START BY SENDING MY OWN PACKAGE alphaCD
#ifdef VAR_MPI
       !A given rank always recieve a package from the same node 
       !rank 0 recieves from rank 1, rank 1 recieves from 2 .. 
       Receiver = MOD(1+mynum,numnodes)
       !A given rank always send to the same node 
       !rank 2 sends to rank 1, rank 1 sends to rank 0 ...
       Sender = MOD(mynum-1+numnodes,numnodes)
       COUNT = MynauxMPI*nocc*nvirt !size of alphaCD
       IF(MynauxMPI.GT.0)THEN !only send package if I have been assigned some basis functions
          call Test_if_64bit_integer_required(MynauxMPI,nocc,nvirt)
          call MPI_ISEND(AlphaCD,COUNT,MPI_DOUBLE_PRECISION,Sender,TAG,comm,request,ierr)
          call MPI_Request_free(request,ierr)
       ENDIF
       IF(MynauxMPI.GT.0)THEN 
          !consider 
          !c_(alpha,ai) = (alpha|beta)^(-1/2) (beta|ai)
          !so that you only need 1 3dim quantity      
          call mem_alloc(Calpha,MynauxMPI,nvirt,nocc)
          call mem_leaktool_alloc(Calpha,LT_Calpha)
          !Use own AlphaCD to obtain part of Calpha
          call RIMP2_buildOwnCalphaFromAlphaCD(nocc,nvirt,mynum,numnodes,&
               & natoms,MynauxMPI,nAtomsMPI,startAuxMPI,nAuxMPI,AlphaCD,&
               & Calpha,TMPAlphaBeta_inv,naux)
       ENDIF
       !To complete construction of  c_(nauxMPI,nvirt,nocc) we need all
       !alphaCD(nauxMPI,nvirt,nocc) contributions from all ranks
       !so we do:
       ! 1. MPI recieve a alphaCD from 'Receiver' 
       !        the first package should already have arrived 
       !        originating from the ISEND immidiately after
       !        II_get_RI_AlphaCD_3centerInt2
       ! 2. Obtain part of Calpha from this contribution
       ! 3. MPI send the recieved alphaCD to 'Sender' 
       ! 4. Repeat untill all contributions have been added

       useAlphaCD5 = .TRUE. 
       useAlphaCD6 = .FALSE.

       CurrentWait(1) = 0
       CurrentWait(2) = 0
       nAwaitDealloc = 0
       DO node=1,numnodes-1 !should recieve numnodes-1 packages 
          !When node=1 the package rank 0 recieves is from rank 1 and was created on rank 1
          !When node=2 the package rank 0 recieves is from rank 1 but was originally created on rank 2
          ! ...
          !myOriginalRank therefore determine the size of NodenauxMPI
          myOriginalRank = MOD(mynum+node,numnodes)
          OriginalRanknauxMPI = nbasisauxMPI(myOriginalRank+1) !dim1 of recieved package
          IF(OriginalRanknauxMPI.GT.0)THEN
             call Test_if_64bit_integer_required(OriginalRanknauxMPI,nocc,nvirt)
             !Step 1 : MPI recieve a alphaCD from 'Receiver' 
             IF(nAwaitDealloc.EQ.2)THEN
                !all buffers are allocated and await to be deallocated once the memory
                !have been recieved by the reciever.
                IF(CurrentWait(1).EQ.5)THEN
                   call MPI_WAIT(request5,status,ierr)
                   call mem_leaktool_dealloc(AlphaCD5,LT_AlphaCD5)
                   call mem_dealloc(AlphaCD5)
                ELSEIF(CurrentWait(1).EQ.6)THEN
                   call MPI_WAIT(request6,status,ierr)
                   call mem_leaktool_dealloc(AlphaCD6,LT_AlphaCD6)
                   call mem_dealloc(AlphaCD6)
                ENDIF
                nAwaitDealloc = 1
                CurrentWait(1) = CurrentWait(2)
                CurrentWait(2) = 0 
             ENDIF
             IF(useAlphaCD5)THEN
                call mem_alloc(AlphaCD5,OriginalRanknauxMPI,nvirt,nocc)
                call mem_leaktool_alloc(AlphaCD5,LT_AlphaCD5)
             ELSEIF(useAlphaCD6)THEN
                call mem_alloc(AlphaCD6,OriginalRanknauxMPI,nvirt,nocc)
                call mem_leaktool_alloc(AlphaCD6,LT_AlphaCD6)
             ENDIF
             COUNT = OriginalRanknauxMPI*nocc*nvirt
             MessageRecieved = .FALSE.
             IF(useAlphaCD5)THEN
                call MPI_RECV(AlphaCD5,COUNT,MPI_DOUBLE_PRECISION,Receiver,TAG,comm,status,ierr)
                call MPI_ISEND(AlphaCD5,COUNT,MPI_DOUBLE_PRECISION,Sender,TAG,comm,request5,ierr)
             ELSEIF(useAlphaCD6)THEN
                call MPI_RECV(AlphaCD6,COUNT,MPI_DOUBLE_PRECISION,Receiver,TAG,comm,status,ierr)
                call MPI_ISEND(AlphaCD6,COUNT,MPI_DOUBLE_PRECISION,Sender,TAG,comm,request6,ierr)
             ENDIF
             IF(MynauxMPI.GT.0)THEN
                !Step 2: Obtain part of Calpha from this contribution
                IF(useAlphaCD5)THEN
                   call RIMP2_buildCalphaContFromAlphaCD(nocc,nvirt,myOriginalRank,numnodes,natoms,&
                        & OriginalRanknauxMPI,MynauxMPI,nAtomsMPI,startAuxMPI,nAuxMPI,&
                        & AlphaCD5,Calpha,TMPAlphaBeta_inv,naux)
                ELSEIF(useAlphaCD6)THEN
                   call RIMP2_buildCalphaContFromAlphaCD(nocc,nvirt,myOriginalRank,numnodes,natoms,&
                        & OriginalRanknauxMPI,MynauxMPI,nAtomsMPI,startAuxMPI,nAuxMPI,&
                        & AlphaCD6,Calpha,TMPAlphaBeta_inv,naux)
                ENDIF
             ENDIF
             !Step 3: MPI send the recieved alphaCD to 'Sender' 
             IF(node.NE.numnodes-1)THEN
                IF(useAlphaCD5)THEN
                   useAlphaCD5 = .FALSE.; useAlphaCD6=.TRUE.
                   nAwaitDealloc = nAwaitDealloc + 1
                   CurrentWait(nAwaitDealloc) = 5
                ELSEIF(useAlphaCD6)THEN
                   useAlphaCD6 = .FALSE.; useAlphaCD5=.TRUE.
                   nAwaitDealloc = nAwaitDealloc + 1
                   CurrentWait(nAwaitDealloc) = 6
                ENDIF
             ELSE
                IF(useAlphaCD5)THEN
                   call mem_leaktool_dealloc(AlphaCD5,LT_AlphaCD5)
                   call mem_dealloc(AlphaCD5)
                ELSEIF(useAlphaCD6)THEN
                   call mem_leaktool_dealloc(AlphaCD6,LT_AlphaCD6)
                   call mem_dealloc(AlphaCD6)
                ENDIF
             ENDIF
          ENDIF
       ENDDO
       IF(MynauxMPI.GT.0)THEN
          call mem_LeakTool_dealloc(TMPAlphaBeta_inv,LT_RIMP2_TMPAlphaBeta_inv)
          call mem_dealloc(TMPAlphaBeta_inv)
       ENDIF
       NBA = MynauxMPI
       IF(nAwaitDealloc.NE.0)THEN
          do iAwaitDealloc=1,nAwaitDealloc
             IF(CurrentWait(iAwaitDealloc).EQ.5)THEN
                call MPI_WAIT(request5,status,ierr)
                call mem_leaktool_dealloc(AlphaCD5,LT_AlphaCD5)
                call mem_dealloc(AlphaCD5)
             ELSEIF(CurrentWait(iAwaitDealloc).EQ.6)THEN
                call MPI_WAIT(request6,status,ierr)
                call mem_leaktool_dealloc(AlphaCD6,LT_AlphaCD6)
                call mem_dealloc(AlphaCD6)
             ENDIF
          enddo
       ENDIF
#endif
    ELSE
       call mem_alloc(Calpha,naux,nvirt,nocc)
       call mem_leaktool_alloc(Calpha,LT_Calpha)
       !Calpha(naux,nvirt,nocc) = AlphaBeta_inv(naux,naux)*AlphaCD(naux,nvirt,nocc)
       call DGEMM('N','N',naux,nvirt*nocc,naux,1.0E0_realk,AlphaBeta_inv,&
            & naux,AlphaCD,naux,0.0E0_realk,Calpha,naux)
       call mem_LeakTool_dealloc(AlphaBeta_inv,LT_RIMP2_AlphaBeta_inv)
       call mem_dealloc(AlphaBeta_inv)
       NBA = naux
    ENDIF

    call mem_alloc(EpsOcc,nocc)
    call mem_leaktool_alloc(EpsOcc,LT_Eps)
    do I=1,nocc
       EpsOcc(I) = MyMolecule%ppfock(I,I)
    enddo
    call mem_alloc(EpsVirt,nvirt)
    call mem_leaktool_alloc(EpsVirt,LT_Eps)
    do A=1,nvirt
       EpsVirt(A) = MyMolecule%qqfock(A,A)
    enddo

    rimp2_energy = 0.0E0_realk

    nullify(AlphaCD2)
    nullify(Calpha2)
    IF(Wakeslaves)THEN
#ifdef VAR_MPI
       IF(MynauxMPI.GT.0)THEN 
          !Energy = sum_{AIBJ} (AI|BJ)_N*[ 2(AI|BJ)_N - (BI|AJ)_N ]/(epsI+epsJ-epsA-epsB)
          call RIMP2_CalcOwnEnergyContribution(nocc,nvirt,EpsOcc,EpsVirt,&
               & NBA,alphaCD,Calpha,rimp2_energy)
          !Energy = sum_{AIBJ} (AI|BJ)_K*[ 2(AI|BJ)_N - (BI|AJ)_N ]/(epsI+epsJ-epsA-epsB)
          COUNT = NBA*nocc*nvirt
          call MPI_ISEND(AlphaCD,COUNT,MPI_DOUBLE_PRECISION,Sender,TAG1,comm,request1,ierr)
!          call MPI_Request_free(request1,ierr)
          call MPI_ISEND(Calpha,COUNT,MPI_DOUBLE_PRECISION,Sender,TAG1,comm,request2,ierr)
!          call MPI_Request_free(request2,ierr)
       ENDIF
       RoundRobin = .FALSE.
       DO node=1,numnodes-1 !should recieve numnodes-1 packages 
          myOriginalRank = MOD(mynum+node,numnodes)         
          OriginalRanknauxMPI = nBasisauxMPI(myOriginalRank+1) !dim1 of recieved package
          IF(OriginalRanknauxMPI.GT.0)THEN
             IF(associated(AlphaCD2))THEN
                call MPI_WAIT(request5,status,ierr)
                call mem_leaktool_dealloc(AlphaCD2,LT_AlphaCD2)
                call mem_dealloc(AlphaCD2)              
                call MPI_WAIT(request6,status,ierr)
                call mem_leaktool_dealloc(Calpha2,LT_Calpha2)
                call mem_dealloc(Calpha2)
                nullify(AlphaCD2)
                nullify(Calpha2)                
             ENDIF
             call mem_alloc(AlphaCD2,OriginalRanknauxMPI,nvirt,nocc)
             call mem_leaktool_alloc(AlphaCD2,LT_AlphaCD2)
             call mem_alloc(Calpha2,OriginalRanknauxMPI,nvirt,nocc)
             call mem_leaktool_alloc(Calpha2,LT_Calpha2)
             COUNT = OriginalRanknauxMPI*nvirt*nocc
             IF(node.EQ.1)THEN
                !recieve from the ISEND above 
                call MPI_RECV(AlphaCD2,COUNT,MPI_DOUBLE_PRECISION,Receiver,TAG1,comm,status,ierr)
                call MPI_RECV(Calpha2,COUNT,MPI_DOUBLE_PRECISION,Receiver,TAG1,comm,status,ierr)
             ELSE
                call MPI_RECV(AlphaCD2,COUNT,MPI_DOUBLE_PRECISION,Receiver,TAG2,comm,status,ierr)
                call MPI_RECV(Calpha2,COUNT,MPI_DOUBLE_PRECISION,Receiver,TAG2,comm,status,ierr)
             ENDIF
             IF(node.NE.numnodes-1)THEN
                RoundRobin = .TRUE.
                call MPI_ISEND(AlphaCD2,COUNT,MPI_DOUBLE_PRECISION,Sender,TAG2,comm,request5,ierr)
                call MPI_ISEND(Calpha2,COUNT,MPI_DOUBLE_PRECISION,Sender,TAG2,comm,request6,ierr)
             ENDIF
             IF(MynauxMPI.GT.0)THEN
                call RIMP2_CalcEnergyContribution(nocc,nvirt,EpsOcc,EpsVirt,&
                     & NBA,alphaCD,Calpha,alphaCD2,Calpha2,OriginalRanknauxMPI,rimp2_energy)
             ENDIF
          ENDIF
       ENDDO
       IF(RoundRobin)THEN
          IF(associated(AlphaCD2))THEN
             call MPI_WAIT(request5,status,ierr)
             call mem_leaktool_dealloc(AlphaCD2,LT_AlphaCD2)
             call mem_dealloc(AlphaCD2)
             call MPI_WAIT(request6,status,ierr)
             call mem_leaktool_dealloc(Calpha2,LT_Calpha2)
             call mem_dealloc(Calpha2)              
          ENDIF
       ELSE
          call mem_leaktool_dealloc(AlphaCD2,LT_AlphaCD2)
          call mem_dealloc(AlphaCD2)
          call mem_leaktool_dealloc(Calpha2,LT_Calpha2)
          call mem_dealloc(Calpha2)              
       ENDIF
       IF(MynauxMPI.GT.0)THEN
          call MPI_WAIT(request1,status,ierr)
          call mem_leaktool_dealloc(AlphaCD,LT_AlphaCD)
          call mem_dealloc(AlphaCD)              
          call MPI_WAIT(request2,status,ierr)
          call mem_leaktool_dealloc(Calpha,LT_Calpha)
          call mem_dealloc(Calpha)              
       ENDIF
       call mem_dealloc(nbasisauxMPI)
       call mem_dealloc(startAuxMPI)
       call mem_dealloc(nAtomsMPI)
       call mem_dealloc(nAuxMPI)
#endif
    ELSE
       !Energy = sum_{AIBJ} (AI|BJ)_N*[ 2(AI|BJ)_N - (BI|AJ)_N ]/(epsI+epsJ-epsA-epsB)
       call RIMP2_CalcOwnEnergyContribution(nocc,nvirt,EpsOcc,EpsVirt,&
            & NBA,alphaCD,Calpha,rimp2_energy)
       call mem_leaktool_dealloc(AlphaCD,LT_AlphaCD)
       call mem_dealloc(AlphaCD)              
       call mem_leaktool_dealloc(Calpha,LT_Calpha)
       call mem_dealloc(Calpha)              
    ENDIF
    call mem_leaktool_dealloc(EpsOcc,LT_Eps)
    call mem_dealloc(EpsOcc)
    call mem_leaktool_dealloc(EpsVirt,LT_Eps)
    call mem_dealloc(EpsVirt)
#ifdef VAR_MPI
    call lsmpi_reduction(rimp2_energy,infpar%master,comm)
#endif    
    IF(MASTER)THEN
       write(DECinfo%output,*)  'RIMP2 CORRELATION ENERGY = ', rimp2_energy
       write(*,'(1X,a,f20.10)') 'RIMP2 CORRELATION ENERGY = ', rimp2_energy
    ENDIF

    write(DECinfo%output,*)  'LEAK TOOL STATISTICS IN full_canonical_rimp2'
    call LeakTools_stat_mem(DECinfo%output)

  end subroutine full_canonical_rimp2

! Calculate canonical MP2 energy
subroutine RIMP2_CalcOwnEnergyContribution(nocc,nvirt,EpsOcc,EpsVirt,NBA,alphaCD,Calpha,rimp2_energy)
  implicit none
  integer,intent(in) :: nocc,nvirt,NBA
  real(realk),intent(in) :: EpsOcc(nocc),EpsVirt(nvirt)
  real(realk),intent(in) :: alphaCD(NBA,nvirt,nocc)
  real(realk),intent(in) :: Calpha(NBA,nvirt,nocc)
  real(realk),intent(inout) :: rimp2_energy
  !
  integer :: J,B,A,I,ALPHA
  real(realk) :: eps,gmoAIBJ,gmoBIAJ,TMP,epsIJB
  Tmp = 0.0E0_realk
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(NONE) &
  !$OMP PRIVATE(J,B,A,I,eps,gmoAIBJ,gmoBIAJ,epsIJB,ALPHA) REDUCTION(+:TMP) &
  !$OMP SHARED(nocc,nvirt,EpsOcc,EpsVirt,NBA,alphaCD,Calpha)
  do J=1,nocc
     do B=1,nvirt
        do I=1,nocc
           epsIJB = EpsOcc(I) + EpsOcc(J) - EpsVirt(B)
           do A=1,nvirt
              ! Difference in orbital energies: eps(I) + eps(J) - eps(A) - eps(B)
              eps = epsIJB - EpsVirt(A)
              gmoAIBJ = 0.0E0_realk
              DO ALPHA = 1,NBA
                 gmoAIBJ = gmoAIBJ + alphaCD(ALPHA,A,I)*Calpha(ALPHA,B,J)
              ENDDO
              gmoBIAJ = 0.0E0_realk                   
              DO ALPHA = 1,NBA
                 gmoBIAJ = gmoBIAJ + alphaCD(ALPHA,B,I)*Calpha(ALPHA,A,J)
              ENDDO
              !Energy = sum_{AIBJ} (AI|BJ)*[ 2(AI|BJ) - (BI|AJ) ]/(epsI + epsJ - epsA - epsB)
              Tmp = Tmp + gmoAIBJ*(2E0_realk*gmoAIBJ-gmoBIAJ)/eps
           end do
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  rimp2_energy = rimp2_energy + Tmp
end subroutine RIMP2_CalcOwnEnergyContribution

! Calculate canonical MP2 energy
!Energy = sum_{AIBJ} (AI|BJ)_K*[ 2(AI|BJ)_N - (BI|AJ)_N ]/(epsI+epsJ-epsA-epsB)
subroutine RIMP2_CalcEnergyContribution(nocc,nvirt,EpsOcc,EpsVirt,NBA,&
     & alphaCD,Calpha,alphaCD2,Calpha2,NBA2,rimp2_energy)
  implicit none
  integer,intent(in) :: nocc,nvirt,NBA,NBA2
  real(realk),intent(in) :: EpsOcc(nocc),EpsVirt(nvirt)
  real(realk),intent(in) :: alphaCD(NBA,nvirt,nocc)
  real(realk),intent(in) :: Calpha(NBA,nvirt,nocc)
  real(realk),intent(in) :: alphaCD2(NBA2,nvirt,nocc)
  real(realk),intent(in) :: Calpha2(NBA2,nvirt,nocc)
  real(realk),intent(inout) :: rimp2_energy
  !
  integer :: J,B,A,I,ALPHA
  real(realk) :: eps,gmoAIBJ,gmoBIAJ,TMP,epsIJB,gmoAIBJ2
  Tmp = 0.0E0_realk
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(NONE) &
  !$OMP PRIVATE(J,B,A,I,eps,gmoAIBJ,gmoBIAJ,epsIJB,gmoAIBJ2,ALPHA) &
  !$OMP REDUCTION(+:TMP) &
  !$OMP SHARED(nocc,nvirt,EpsOcc,EpsVirt,NBA,alphaCD,Calpha,NBA2,alphaCD2,&
  !$OMP Calpha2)
  do J=1,nocc
     do B=1,nvirt
        do I=1,nocc
           epsIJB = EpsOcc(I) + EpsOcc(J) - EpsVirt(B)
           do A=1,nvirt
              ! Difference in orbital energies: eps(I)+eps(J)-eps(A)-eps(B)
              eps = epsIJB - EpsVirt(A)
              gmoAIBJ2 = 0.0E0_realk
              DO ALPHA = 1,NBA2
                 gmoAIBJ2 = gmoAIBJ2 + alphaCD2(ALPHA,A,I)*Calpha2(ALPHA,B,J)
              ENDDO
              gmoAIBJ = 0.0E0_realk
              DO ALPHA = 1,NBA
                 gmoAIBJ = gmoAIBJ + alphaCD(ALPHA,A,I)*Calpha(ALPHA,B,J)
              ENDDO
              gmoBIAJ = 0.0E0_realk                   
              DO ALPHA = 1,NBA
                 gmoBIAJ = gmoBIAJ + alphaCD(ALPHA,B,I)*Calpha(ALPHA,A,J)
              ENDDO
              !Energy = sum_{AIBJ} (AI|BJ)*[ 2(AI|BJ) - (BI|AJ) ]/(epsI + epsJ - epsA - epsB)
              Tmp = Tmp + gmoAIBJ2*(2E0_realk*gmoAIBJ-gmoBIAJ)/eps
           end do
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  rimp2_energy = rimp2_energy + Tmp
end subroutine RIMP2_CalcEnergyContribution

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

  !> ***********************************************************************
  !> ****************************  CCSD F12 ********************************
  !> ***********************************************************************

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
    type(tensor) :: Tai
    type(tensor) :: Taibj
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
    !real(realk),pointer :: Vijij_term1(:,:)

    real(realk),pointer :: Vjiij(:,:)    
    !real(realk),pointer :: Vjiij_term1(:,:)

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
                Ttmp = Tai%elm2(a,i)*Tai%elm2(b,j) + Taibj%elm4(a,i,b,j)
                Gtmp = 2.0E0_realk * AIBJ(a,i,b,j) - AIBJ(b,i,a,j)
                ECCSD = ECCSD + Ttmp * Gtmp
             end do
          end do
       end do
    end do

    write(*,'(1X,a)') '----------------------------------------------------'  
    write(*,'(1X,a,f20.10)') 'TOYCODE: CCSD CORRELATION ENERGY = ', ECCSD
    
    !> DEC Input
    write(DECinfo%output,*) 'TOYCODE: CCSD CORRELATION ENERGY = ', ECCSD

    ! INTEGRAL GUYS: DO YOUR MAGIC FOR THE F12 CONTRIBUTION HERE...
    call mem_alloc(Vijij,nocc,nocc)
    call mem_alloc(Vjiij,nocc,nocc)
    call mp2f12_Vijij(Vijij,Ripjq,Gipjq,Fijkl,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)
    call mp2f12_Vjiij(Vjiij,Ripjq,Gipjq,Fijkl,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)

    call mem_alloc(Ciajb,nocc,nvirt,nocc,nvirt)
    call mem_alloc(Cjaib,nocc,nvirt,nocc,nvirt)
    call mp2f12_Ciajb(Ciajb,Giajc,Fac%elms,nocc,nvirt,ncabs)

    !> --------------------------------------
    !>               Viajb ccsd
    !>---------------------------------------
    !V_ij^ab
    call mem_alloc(Viajb,nocc,nvirt,nocc,nvirt) 

    call ccsdf12_Viajb(Viajb,Rapbq,Gipjq,Fiajb,Rambc,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
   
    if(DECinfo%F12DEBUG) then

!!$       print *, '----------------------------------------'
!!$       print *, '  E21  V terms                          '
!!$       print *, '----------------------------------------'
!!$       print *, 'E21_V_term1: ', 2.0E0_REALK*mp2f12_E21(Vijij_term1, Vjiij_term1,nocc)
!!$       print *, 'E21_V_term2: ', 2.0E0_REALK*mp2f12_E21(Vijij_term2, Vjiij_term2,nocc)
!!$       print *, 'E21_V_term3: ', 2.0E0_REALK*mp2f12_E21(Vijij_term3, Vjiij_term3,nocc)
!!$       print *, 'E21_V_term4: ', 2.0E0_REALK*mp2f12_E21(Vijij_term4, Vjiij_term4,nocc)
!!$       print *, 'E21_V_term5: ', 2.0E0_REALK*mp2f12_E21(Vijij_term5, Vjiij_term5,nocc)
!!$       print *, '----------------------------------------'
!!$
!!$       E21_debug = 2.0E0_REALK*(mp2f12_E21(Vijij_term1,Vjiij_term1,nocc) + mp2f12_E21(Vijij_term2,Vjiij_term2,nocc) &
!!$            & + mp2f12_E21(Vijij_term3,Vjiij_term3,nocc) + mp2f12_E21(Vijij_term4,Vjiij_term4,nocc) &
!!$            & + mp2f12_E21(Vijij_term5,Vjiij_term5,nocc)) 
!!$
!!$       print *, 'E21_Vsum: ', E21_debug
!!$      !print *, 'E21_debug: ', 2.0E0_REALK*mp2f12_E21(Vijij,Vjiij,nocc)

    endif

    !V_ij^ia
    call mem_alloc(Viija,nocc,nocc,nvirt)
    call ccsdf12_Viija(Viija,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)

    !V_ij^ja
    call mem_alloc(Vijja,nocc,nocc,nvirt)
    call ccsdf12_Vijja(Vijja,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)

    !V_ij^ai
    call mem_alloc(Viaji,nocc,nvirt,nocc)
    call ccsdf12_Viaji(Viaji,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)

    !V_ij^aj
    call mem_alloc(Viajj,nocc,nvirt,nocc)
    call ccsdf12_Viajj(Viajj,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)

    ! CCSD coupling V_ij^ij & V_ji^ij 
    call ccsdf12_Vijij_coupling(Vijij,Ciajb,Taibj%elm4,Viajb,Viija,Viajj,Tai%elm2,nocc,nvirt)
    call ccsdf12_Vjiij_coupling(Vjiij,Ciajb,Taibj%elm4,Viajb,Vijja,Viaji,Tai%elm2,nocc,nvirt)

    !call ccsdf12_Vijij_coupling_term1(Vijij_term1,Ciajb,Taibj%val,Viajb,Viija,Viajj,Tai%val,nocc,nvirt)
    !call ccsdf12_Vjiij_coupling_term1(Vjiij_term1,Ciajb,Taibj%val,Viajb,Vijja,Viaji,Tai%val,nocc,nvirt)
    
    ! CCSD Specific
    call mem_dealloc(Viija)
    call mem_dealloc(Vijja)
    call mem_dealloc(Viajj)
    call mem_dealloc(Viaji)
    call mem_dealloc(Viajb)

    E21 = 2.0E0_realk*mp2f12_E21(Vijij,Vjiij,nocc)
    
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

    write(*,'(1X,a)') '----------------------------------------------------'  
    write(*,'(1X,a,f20.10)') 'TOYCODE: CCSD-F12 CORRECTION TO ENERGY  = ', EF12
    write(*,'(1X,a,f20.10)') 'TOYCODE: CCSD-F12 CORRELATION ENERGY    = ', ECCSD_F12
    
    !> Input to DEC
    write(DECinfo%output,*) 'TOYCODE: CCSD-F12 CORRECTION TO ENERGY  = ', EF12
    write(DECinfo%output,*) 'TOYCODE: CCSD-F12 CORRELATION ENERGY = ', ECCSD_F12

    call mem_dealloc(AIBJ)
    call tensor_free(Taibj)
    call tensor_free(Tai)
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
  subroutine full_get_ccsd_singles_and_doubles(MyMolecule,MyLsitem,Tai_local, Taibj_local)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Singles amplitudes
    type(tensor),intent(inout) :: Tai_local
    !> Doubles amplitudes
    type(tensor),intent(inout) :: Taibj_local
    integer :: nocc,nunocc,nbasis,print_level,save_model,startidx,endidx,i,j
    logical :: fragment_job
    real(realk) :: energy
    type(tensor) :: VOVO,Tai,Taibj
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
            & print_level,&
            & ppfock,MyMolecule%qqfock,energy, Tai, Taibj,&
            & VOVO,.false.,local,DECinfo%use_pnos)
       call mem_dealloc(ppfock)

    else

       call ccsolver_par(DECinfo%ccmodel,MyMolecule%Co,MyMolecule%Cv,&
            & MyMolecule%fock, nbasis,nocc,nunocc,mylsitem, print_level, &
            & MyMolecule%ppfock,MyMolecule%qqfock,&
            & energy, Tai, Taibj, VOVO,.false.,local,DECinfo%use_pnos)

    end if

    call tensor_free(VOVO)

    !convert the parallel distributed quantities to local quantities, this is
    !essentially a copying if the tensors are not parallel distributed
    call tensor_minit(Tai_local,[nunocc,nocc],2,atype="LDAR")
    call tensor_add(Tai_local,1.0E0_realk,Tai,a=0.0E0_realk)
    call tensor_minit(Taibj_local,[nunocc,nocc,nunocc,nocc],4,atype="LDAR")
    call tensor_add(Taibj_local,1.0E0_realk,Taibj,a=0.0E0_realk)

    call tensor_free(Taibj)
    call tensor_free(Tai)


  end subroutine full_get_ccsd_singles_and_doubles

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

#ifdef VAR_MPI
subroutine full_canonical_rimp2_slave
  use full,only: full_canonical_rimp2
  use infpar_module !infpar
  use lsmpi_type,only:ls_mpiInitBuffer,ls_mpiFinalizeBuffer,&
       & LSMPIBROADCAST,MPI_COMM_LSDALTON 
  use lsmpi_op,only: mpicopy_lsitem
  use precision
  use typedeftype,only:lsitem
  use lsparameters
  use decmpi_module, only: mpi_bcast_fullmolecule
  use DALTONINFO, only: ls_free
!  use typedef
!  use dec_typedef_module
  ! DEC DEPENDENCIES (within deccc directory)   
  ! *****************************************
!  use dec_fragment_utils
  use full_molecule, only:fullmolecule, molecule_finalize
  implicit none
  !> Full molecule info
  type(fullmolecule) :: MyMolecule
  !> Lsitem structure
  type(lsitem) :: mylsitem
  !> Canonical MP2 correlation energy
  real(realk) :: rimp2_energy    
  
  ! Init MPI buffer
  ! ***************
  ! Main master:  Prepare for writing to buffer
  ! Local master: Receive buffer
  call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,MPI_COMM_LSDALTON)
  ! Integral lsitem
  ! ---------------
  call mpicopy_lsitem(MyLsitem,MPI_COMM_LSDALTON)
  
  call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,MPI_COMM_LSDALTON)
  ! Full molecule bcasting
  ! **********************
  call mpi_bcast_fullmolecule(MyMolecule)
  
  ! Finalize MPI buffer
  ! *******************
  ! Main master:  Send stuff to local masters and deallocate temp. buffers
  ! Local master: Deallocate buffer etc.
  call full_canonical_rimp2(MyMolecule,MyLsitem,rimp2_energy)

  call ls_free(MyLsitem)
  call molecule_finalize(MyMolecule)
  
end subroutine full_canonical_rimp2_slave
#endif

