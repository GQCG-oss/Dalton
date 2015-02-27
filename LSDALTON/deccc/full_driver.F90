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
  !  DEC DEPENDENCIES (within deccc directory)   
  !  *****************************************
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

  public :: full_driver, full_canonical_rimp2, full_canonical_mp2, &
       & full_canonical_mp2B, canonical_mp2B_memreq_test
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
    real(realk) :: Eerr,Edft
    logical :: Success
    write(DECinfo%output,'(/,a)') ' ================================================ '
    write(DECinfo%output,'(a)')   '              Full molecular driver               '
    write(DECinfo%output,'(a,/)') ' ================================================ '

#ifdef VAR_MPI
    call set_dec_settings_on_slaves()
#endif


    ! For SNOOP we might want to skip correlated calculation
    DoCorrelatedCalculation: if(DECinfo%SNOOPjustHF) then

       write(DECinfo%output,*) 'Skipping correlated calculation for SNOOP!'
       Ecorr = 0.0_realk

    else


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
          if(DECinfo%ccModel==MODEL_MP2) then
             if(DECinfo%use_canonical ) then
                !simple conventional MP2 calculation only works for canonical orbitals
                !no amplitudes stored. MP2B requires (nb,nb,nb) can be fully distributed
                call full_canonical_mp2B(MyMolecule,MyLsitem,Ecorr)       
             else
                !Call routine which calculates individual fragment 
                !contributions and prints them,
                !works both for canonical and local orbitals
                !call Full_DECMP2_calculation(MyMolecule,mylsitem,Ecorr)
                call full_cc_dispatch(MyMolecule,mylsitem,Ecorr)          
             end if
          else
             call full_cc_dispatch(MyMolecule,mylsitem,Ecorr)          
          end if
       end if

       ! Get HF energy
       Ehf = get_HF_energy_fullmolecule(MyMolecule,Mylsitem,D)
       !DFT energy
       if(DECinfo%DFTreference) then
         Edft = get_dft_energy_fullmolecule(MyMolecule,Mylsitem,D) 
       endif

       ! Print summary, set the error estimates to zero
       call print_total_energy_summary(EHF,Edft,Ecorr,0.0E0_realk,0.0E0_realk,0.0E0_realk)

    end if DoCorrelatedCalculation

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
    !> MP2-F12 correlation energy
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

    !> Singles correction
    type(matrix) :: Fcd
    type(matrix) :: Fic  
        
    !> Singles correction energy
    real(realk)  :: ES2
     
    real(realk)  :: E21, E21_debug, E22, E22_debug, E23_debug, Gtmp
    type(tensor) :: tensor_Taibj,tensor_gmo
    integer :: vs, os
    logical :: local
    local = .true.
#ifdef VAR_MPI
    local = (infpar%lg_nodtot==1)
#endif
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

    ! Get all F12 Fock Matrices
    ! ********************
    call get_F12_mixed_MO_Matrices(MyLsitem,MyMolecule,Dmat,nbasis,ncabsAO,&
         & nocc,noccfull,nvirt,ncabs,HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp,Fic,Fcd)

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
    
    ! MP2-F12 Singles correction (Yang M. Wang 03.12.2014)
    ! ***************************    
    call get_ES2(ES2,Fic,Fii,Fcd,nocc,ncabs)
   
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
       ! PE: even more quicker and dirtier solution for that the new solver
       ! needs the tensor format
       call mp2_solver(MyMolecule,mylsitem,tensor_gmo,tensor_Taibj,.false.)
       call tensor_free(tensor_gmo)

       call mem_alloc(Taibj,nvirt,nocc,nvirt,nocc)
       call tensor_convert(tensor_Taibj,Taibj)
       call tensor_free(tensor_Taibj)

       mp2_energy = 0.0E0_realk
       do j=1,nocc
          do b=1,nvirt
             do i=1,nocc
                do a=1,nvirt
                   ! Energy = sum_{ijab} ( Taibj) * (ai | bj)
                   !Taibj(a,i,b,j) = array4Taibj%elm4(a,i,b,j)
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
       print *, ' E21 V terms                            '
       print *, '----------------------------------------'
       print *, ' E21_V_term1: ', 2.0E0_REALK*mp2f12_E21(Vijij_term1,Vjiij_term1,nocc)
       print *, ' E21_V_term2: ', 2.0E0_REALK*mp2f12_E21(Vijij_term2,Vjiij_term2,nocc)
       print *, ' E21_V_term3: ', 2.0E0_REALK*mp2f12_E21(Vijij_term3,Vjiij_term3,nocc)
       print *, ' E21_V_term4: ', 2.0E0_REALK*mp2f12_E21(Vijij_term4,Vjiij_term4,nocc)
       print *, ' E21_V_term5: ', 2.0E0_REALK*mp2f12_E21(Vijij_term5,Vjiij_term5,nocc)
       print *, '----------------------------------------'

       E21_debug = 2.0E0_REALK*(mp2f12_E21(Vijij_term1,Vjiij_term1,nocc) + mp2f12_E21(Vijij_term2,Vjiij_term2,nocc) &
            & + mp2f12_E21(Vijij_term3,Vjiij_term3,nocc) + mp2f12_E21(Vijij_term4,Vjiij_term4,nocc) &
            & + mp2f12_E21(Vijij_term5,Vjiij_term5,nocc)) 

       print *, ' E21_Vsum: ', E21_debug
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
 
       if(DECinfo%use_canonical) then
          call mem_alloc(Xijij_term1,nocc,nocc)
          call mem_alloc(Xijij_term2,nocc,nocc)
          call mem_alloc(Xijij_term3,nocc,nocc)
          call mem_alloc(Xijij_term4,nocc,nocc)  

          call mem_alloc(Xjiij_term1,nocc,nocc)
          call mem_alloc(Xjiij_term2,nocc,nocc)
          call mem_alloc(Xjiij_term3,nocc,nocc)
          call mem_alloc(Xjiij_term4,nocc,nocc)
       else
          call mem_alloc(Xijkl_term1,nocc,nocc,nocc,nocc)
          call mem_alloc(Xijkl_term2,nocc,nocc,nocc,nocc)
          call mem_alloc(Xijkl_term3,nocc,nocc,nocc,nocc)
          call mem_alloc(Xijkl_term4,nocc,nocc,nocc,nocc) 
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
          print *, ' E_22 X term                            '
          print *, '----------------------------------------'
          print *, ' E22_X_term1: ', mp2f12_E22(Xijij_term1,Xjiij_term1,Fii%elms,nocc)
          print *, ' E22_X_term2: ', mp2f12_E22(Xijij_term2,Xjiij_term2,Fii%elms,nocc)
          print *, ' E22_X_term3: ', mp2f12_E22(Xijij_term3,Xjiij_term3,Fii%elms,nocc)
          print *, ' E22_X_term4: ', mp2f12_E22(Xijij_term4,Xjiij_term4,Fii%elms,nocc)
          print *, '----------------------------------------'
          E22_debug = mp2f12_E22(Xijij_term1,Xjiij_term1,Fii%elms,nocc) & 
               & + mp2f12_E22(Xijij_term2,Xjiij_term2,Fii%elms,nocc) &
               & + mp2f12_E22(Xijij_term3,Xjiij_term3,Fii%elms,nocc) + mp2f12_E22(Xijij_term4,Xjiij_term4,Fii%elms,nocc)  
          print *, ' E22_Xsum: ', E22_debug  
          print *, ' E22_debug: ', mp2f12_E22(Xijij,Xjiij,Fii%elms,nocc)
          print *, '----------------------------------------'
          print *, ' E_23 B term                           '
          print *, '----------------------------------------'
          print *, ' E23_B_term1: ', mp2f12_E23(Bijij_term1,Bjiij_term1,nocc)
          print *, ' E23_B_term2: ', mp2f12_E23(Bijij_term2,Bjiij_term2,nocc)
          print *, ' E23_B_term3: ', mp2f12_E23(Bijij_term3,Bjiij_term3,nocc)
          print *, ' E23_B_term4: ', mp2f12_E23(Bijij_term4,Bjiij_term4,nocc)
          print *, ' E23_B_term5: ', mp2f12_E23(Bijij_term5,Bjiij_term5,nocc)
          print *, ' E23_B_term6: ', mp2f12_E23(Bijij_term6,Bjiij_term6,nocc)
          print *, ' E23_B_term7: ', mp2f12_E23(Bijij_term7,Bjiij_term7,nocc)
          print *, ' E23_B_term8: ', mp2f12_E23(Bijij_term8,Bjiij_term8,nocc)
          print *, ' E23_B_term9: ', mp2f12_E23(Bijij_term9,Bjiij_term9,nocc)   
          print *, '----------------------------------------'
          E23_debug = mp2f12_E23(Bijij_term1,Bjiij_term1,nocc) & 
               & + mp2f12_E23(Bijij_term2,Bjiij_term2,nocc) + mp2f12_E23(Bijij_term3,Bjiij_term3,nocc) &
               & + mp2f12_E23(Bijij_term4,Bjiij_term4,nocc) + mp2f12_E23(Bijij_term5,Bjiij_term5,nocc) &
               & + mp2f12_E23(Bijij_term6,Bjiij_term6,nocc) + mp2f12_E23(Bijij_term7,Bjiij_term7,nocc) &
               & + mp2f12_E23(Bijij_term8,Bjiij_term8,nocc) + mp2f12_E23(Bijij_term9,Bjiij_term9,nocc)
          print *, ' E23_Bsum: ',  E23_debug
          print *, ' E23_Bsum_debug: ',  mp2f12_E23(Bijij,Bjiij,nocc)
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
          print *, ' E_22 X term                            '
          print *, '----------------------------------------'
          print *, ' E22_X_term1: ', X1
          print *, ' E22_X_term2: ', X2
          print *, ' E22_X_term3: ', X3
          print *, ' E22_X_term4: ', X4
          print *, '----------------------------------------'
          E22_debug = X1 + X2 + X3 + X4  
          print *, 'E22_Xsum: ', E22_debug  
          print *, '----------------------------------------'
          print *, ' E_23 B term                            '
          print *, '----------------------------------------'
          print *, ' E23_B_term1: ', mp2f12_E23(Bijij_term1,Bjiij_term1,nocc)
          print *, ' E23_B_term2: ', mp2f12_E23(Bijij_term2,Bjiij_term2,nocc)
          print *, ' E23_B_term3: ', mp2f12_E23(Bijij_term3,Bjiij_term3,nocc)
          print *, ' E23_B_term4: ', mp2f12_E23(Bijij_term4,Bjiij_term4,nocc)
          print *, ' E23_B_term5: ', mp2f12_E23(Bijij_term5,Bjiij_term5,nocc)
          print *, ' E23_B_term6: ', mp2f12_E23(Bijij_term6,Bjiij_term6,nocc)
          print *, ' E23_B_term7: ', mp2f12_E23(Bijij_term7,Bjiij_term7,nocc)
          print *, ' E23_B_term8: ', mp2f12_E23(Bijij_term8,Bjiij_term8,nocc)
          print *, ' E23_B_term9: ', mp2f12_E23(Bijij_term9,Bjiij_term9,nocc)   
          print *, '----------------------------------------'
          E23_debug = mp2f12_E23(Bijij_term1,Bjiij_term1,nocc) & 
               & + mp2f12_E23(Bijij_term2,Bjiij_term2,nocc) + mp2f12_E23(Bijij_term3,Bjiij_term3,nocc) &
               & + mp2f12_E23(Bijij_term4,Bjiij_term4,nocc) + mp2f12_E23(Bijij_term5,Bjiij_term5,nocc) &
               & + mp2f12_E23(Bijij_term6,Bjiij_term6,nocc) + mp2f12_E23(Bijij_term7,Bjiij_term7,nocc) &
               & + mp2f12_E23(Bijij_term8,Bjiij_term8,nocc) + mp2f12_E23(Bijij_term9,Bjiij_term9,nocc)
          print *, ' E23_Bsum: ',  E23_debug
          print *, ' E23_Bsum_debug: ',  mp2f12_E23(Bijij,Bjiij,nocc)
          print *, '----------------------------------------'

       endif      
    endif

    call free_F12_mixed_MO_Matrices(HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp,Fic,Fcd)

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
    else !non canonical 

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

       write(*,'(1X,a)') '-----------------------------------------------------------------'
       write(*,'(1X,a,f20.10)') 'TOYCODE: MP2 CORRELATION ENERGY =           ', mp2_energy
       write(*,'(1X,a,f20.10)') 'TOYCODE: F12 E21 CORRECTION TO ENERGY =     ', E21_debug
       write(*,'(1X,a,f20.10)') 'TOYCODE: F12 E22 CORRECTION TO ENERGY =     ', E22_debug
       write(*,'(1X,a,f20.10)') 'TOYCODE: F12 E23 CORRECTION TO ENERGY =     ', E23_debug
       write(*,'(1X,a,f20.10)') 'TOYCODE: F12 E22+E23 CORRECTION TO ENERGY = ', E22_debug + E23_debug
       write(*,'(1X,a)') '-----------------------------------------------------------------'
       write(*,'(1X,a,f20.10)') 'TOYCODE: F12 CORRECTION TO ENERGY =         ', E21_debug+E22_debug+E23_debug
       write(*,'(1X,a,f20.10)') 'TOYCODE: F12 ES2 CORRECTION TO ENERGY =     ', ES2
       write(*,'(1X,a,f20.10)') 'TOYCODE: FULL F12 CORRECTION TO ENERGY =    ', E21_debug+E22_debug+E23_debug+ES2
       write(*,'(1X,a)') '-----------------------------------------------------------------'
       write(*,'(1X,a,f20.10)') 'TOYCODE: MP2-F12 ENERGY =                   ', mp2_energy+E21_debug+E22_debug+E23_debug
       write(*,'(1X,a,f20.10)') 'TOYCODE: MP2-F12-ES2 ENERGY =               ', mp2_energy+E21_debug+E22_debug+E23_debug+ES2
    else

       mp2f12_energy = 0.0E0_realk
       mp2f12_energy = mp2_energy+E21+E22

       write(DECinfo%output,*) '----------------------------------------------------------------'
       write(*,'(1X,a,f20.10)') '----------------------------------------------------------------'
       write(DECinfo%output,*) '                 TOYCODE MP2-F12 CALCULATION                    '
       write(*,'(1X,a,f20.10)') '                 TOYCODE MP2-F12 CALCULATION                    '
       write(DECinfo%output,*) '----------------------------------------------------------------'
       write(*,'(1X,a,f20.10)') '----------------------------------------------------------------'
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

    call mem_dealloc(gmo)

  end subroutine full_canonical_mp2_f12
#endif

  !> \brief Memory check for full_canonical_rimp2 subroutine
  !> \author Thomas Kjaergaard
  !> \date October 2014
  subroutine full_canonical_rimp2_memory_check(nbasis,nocc,nvirt,nAux,&
       & numnodes,MinAtomicnAux,MemoryReduced)
    implicit none
    integer,intent(in) :: nbasis,nocc,nvirt,nAux,numnodes,MinAtomicnAux
    real(realk) :: MemRequired,GB,MemStep1,MemStep2
    logical :: MemoryReduced
    GB = 1.0E9_realk

    ! Check that arrays fit in memory (hard-coded)
    MemStep1=2*real(MinAtomicnAux*nbasis*nbasis)+real(nAux*nocc*nvirt)/numnodes
    MemStep2=2*real(nAux*nocc*nvirt)/numnodes + real(nAux*nAux)
    MemRequired = MAX(MemStep1,MemStep2)
!    MemRequired = 2*real(nAux*nocc*nvirt/numnodes) + real(nAux*nAux)
    MemRequired = MemRequired*realk/GB

    if(MemRequired > 0.80E0_realk*DECinfo%memory) then
       print *, 'RIMP2 Have 2 memory intensive steps:'
       print *, 'Step1:'
       print *, '2 arrays of size MinAtomicnAux*nbasis*nbasis=',&
            & MinAtomicnAux*nbasis*nbasis
       print *, '1 array of size (nAux/numnodes)*nocc*nvirt  =',&
            & (nAux/numnodes)*nocc*nvirt
       print *, 'Resulting in =',MemStep1/GB,' GB'
       print *, 'Step1:'
       print *, '2 arrays of size (nAux/numnodes)*nocc*nvirt  =',&
            & (nAux/numnodes)*nocc*nvirt
       print *, '1 array of size nAux*nAux  =',(nAux/numnodes)*nocc*nvirt
       print *, 'Resulting in =',MemStep2/GB,' GB'
       print *, 'With'
       print *, 'MinAtomicnAux=',MinAtomicnAux
       print *, 'nbasis       =',nbasis
       print *, 'nAux         =',nAux
       print *, 'RIMP2 Mem required (GB) = ', MemRequired
       print *, 'We restrict ourselves to 80 % of Max memory (GB)'
       print *, 'RIMP2 Max memory (GB)   = ', DECinfo%memory
       call lsquit('full_canonical_rimp2: Memory exceeded! Try to increase memory &
            & using the .MaxMemory (in GB) keyword in the *DEC section.',-1)
    end if

    IF(MemoryReduced)THEN
       WRITE(DECinfo%output,*)'MemoryReduced RIMP2 scheme chosen as default'
    ELSE
       !determine if MemoryReduced scheme must be used
       MemStep1=2*real(MinAtomicnAux*nbasis*nbasis)
       MemStep1=MemStep1+real(nAux*nocc*nvirt)/numnodes
       MemStep2=4*real(nAux*nocc*nvirt)/numnodes + real(nAux*nAux)
       MemRequired = MAX(MemStep1,MemStep2)
       MemRequired = MemRequired*realk/GB       
       if(MemRequired .GE. 0.80E0_realk*DECinfo%memory) then 
          MemoryReduced = .TRUE.
          WRITE(DECinfo%output,*)'MemoryReduced RIMP2 scheme chosen'
          WRITE(DECinfo%output,*)'MemRequired Step 1:',MemStep1*realk/GB,' GB'
          WRITE(DECinfo%output,*)'MemRequired Step 2:',MemStep2*realk/GB,' GB'
          WRITE(DECinfo%output,*)'DECinfo%memory',DECinfo%memory
       endif
    ENDIF

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
    real(realk),pointer :: AlphaCD(:,:,:),AlphaCD5(:,:,:),AlphaCD6(:,:,:)
    real(realk),pointer :: Calpha(:,:,:),Calpha2(:,:,:),Calpha3(:,:,:)
    real(realk),pointer :: AlphaBeta(:,:),AlphaBeta_minus_sqrt(:,:),TMPAlphaBeta_minus_sqrt(:,:),EpsOcc(:),EpsVirt(:)
    integer(kind=ls_mpik)  :: COUNT,TAG,IERR,request,Receiver,sender,COUNT2,comm,TAG1,TAG2
    integer :: J,CurrentWait(2),nAwaitDealloc,iAwaitDealloc,I,NBA,OriginalRanknauxMPI
    integer :: myOriginalRank,node,natoms,MynauxMPI,A,lupri,MinAtomicnAux,restart_lun
    integer :: noccJstart
    logical :: useAlphaCD5,useAlphaCD6,MessageRecieved,RoundRobin,RoundRobin5,RoundRobin6
    logical :: RoundRobin2,RoundRobin3,useCalpha2,useCalpha3,FORCEPRINT
    integer(kind=ls_mpik)  :: request1,request2,request5,request6,request7,request8
    integer,pointer :: nbasisauxMPI(:),startAuxMPI(:,:),AtomsMPI(:,:),nAtomsMPI(:),nAuxMPI(:,:)
    integer(KIND=long) :: MaxMemAllocated,MemAllocated
    real(realk) :: TS,TE,TS2,TE2,TS3,TE3,CPU1,CPU2,WALL1,WALL2,CPU_MPICOMM,WALL_MPICOMM
    real(realk) :: CPU_MPIWAIT,WALL_MPIWAIT,MemInGBCollected,epsIJ
    logical :: MemoryReduced,AlphaCDAlloced,AlphaCD_Deallocate
    logical(kind=ls_mpik) :: TransferCompleted
    logical :: NotMatSet,file_exists
    real(realk),pointer :: Amat(:,:),Bmat(:,:)
    MemoryReduced = MyLsitem%setting%scheme%ForceRIMP2memReduced
#ifdef VAR_TIME    
    FORCEPRINT = .TRUE.
#endif
    CPU_MPICOMM = 0.0E0_realk
    WALL_MPICOMM = 0.0E0_realk
    CPU_MPIWAIT = 0.0E0_realk
    WALL_MPIWAIT = 0.0E0_realk
    !use Memory leak tool
    CALL LSTIMER('START ',TS,TE,DECINFO%OUTPUT)
    MaxMemAllocated = 0
    MemAllocated = 0
    call set_LeakTool_memvar(MaxMemAllocated,MemAllocated)
    TAG = 1411
    TAG1 = 1412
    TAG2 = 1413

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
    LUPRI = DECinfo%output
#ifdef VAR_MPI
    comm = MPI_COMM_LSDALTON
    master= (infpar%mynum == infpar%master)
    mynum = infpar%mynum
    numnodes = infpar%nodtot
    wakeslaves = infpar%nodtot.GT.1
    IF(.NOT.master)LUPRI = 6 !standard Output
#else
    ! If MPI is not used, consider the single node to be "master"
    master=.true.
    mynum = 0
    numnodes = 1
    wakeslaves = .false.
#endif

    ! Memory check!
    ! ********************
    CALL LSTIMER('START ',TS2,TE2,LUPRI,FORCEPRINT)
!    call getMaxAtomicnAux(Mylsitem%SETTING%MOLECULE(1)%p,MaxAtomicnAux,nAtoms)
    call determine_maxBatchOrbitalsize(DECinfo%output,&
         & Mylsitem%setting,MinAtomicnAux,'D')
    call full_canonical_rimp2_memory_check(nbasis,nocc,nvirt,nAux,&
         & numnodes,MinAtomicnAux,MemoryReduced)
    call Test_if_64bit_integer_required(nAux,nAux)
    CALL LSTIMER('RIMP2: MemCheck ',TS2,TE2,LUPRI,FORCEPRINT)
#ifdef VAR_MPI 
    ! Master starts up slave
    StartUpSlaves: if(wakeslaves .and. master) then
       ! Wake up slaves to do the job: slaves awoken up with (RIMP2FULL)
       ! and call full_canonical_rimp2_slave which communicate info 
       ! then calls full_canonical_rimp2.
       CALL LS_GETTIM(CPU1,WALL1)
       call ls_mpibcast(RIMP2FULL,infpar%master,comm)
       ! Communicate fragment information to slaves
       call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,comm)
       call mpicopy_lsitem(MyLsitem,comm)
       call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,comm)
       call mpi_bcast_fullmolecule(MyMolecule)    
       CALL LS_GETTIM(CPU2,WALL2)
       CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
       WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
    endif StartUpSlaves
#endif
    CALL LSTIMER('RIMP2: WakeSlaves ',TS2,TE2,LUPRI,FORCEPRINT)

    IF(master)THEN
       !=====================================================================================
       ! Major Step 1: Master Obtains Overlap (alpha|beta) in Auxiliary Basis 
       !=====================================================================================
       !This part of the Code is NOT MPI/OpenMP parallel - all nodes calculate the full overlap
       !this should naturally be changed      

       call mem_alloc(AlphaBeta,nAux,nAux)    
       CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
       call II_get_RI_AlphaBeta_2centerInt(lupri,lupri,AlphaBeta,Mylsitem%setting,nAux)
       CALL LSTIMER('AlphaBeta ',TS3,TE3,LUPRI,FORCEPRINT)

       !=====================================================================================
       ! Major Step 2: Calculate the inverse (alpha|beta)^(-1) and BCAST
       !=====================================================================================
       ! Warning the inverse is not unique so in order to make sure all slaves have the same
       ! inverse matrix we calculate it on the master and BCAST to slaves

       !Create the squareroot AlphaBeta = (alpha|beta)^-(1/2)
       ! Given matrix S, computes S^{-1/2}.
       
       call mem_alloc(AlphaBeta_minus_sqrt,nAux,nAux)
       call lowdin_diag_S_minus_sqrt(nAux, AlphaBeta,AlphaBeta_minus_sqrt, lupri)
       call mem_dealloc(AlphaBeta)

       CALL LSTIMER('AlphaBetamSq ',TS3,TE3,LUPRI,FORCEPRINT)
    ELSE
       call mem_alloc(AlphaBeta_minus_sqrt,nAux,nAux)
    ENDIF
#ifdef VAR_MPI
    !The barrier is mostly here to detect the time spent waiting vs the 
    !time spent in communication. 
    CALL LS_GETTIM(CPU1,WALL1)
    call lsmpi_barrier(comm)
    CALL LS_GETTIM(CPU2,WALL2)
    CPU_MPIWAIT = CPU_MPIWAIT + (CPU2-CPU1)
    WALL_MPIWAIT = WALL_MPIWAIT + (WALL2-WALL1)
    CALL LS_GETTIM(CPU1,WALL1)
    call ls_mpibcast(AlphaBeta_minus_sqrt,naux,naux,infpar%master,comm)
    CALL LS_GETTIM(CPU2,WALL2)
    CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
    WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
#endif
    CALL LSTIMER('RIMP2: AlphaBeta ',TS2,TE2,LUPRI,FORCEPRINT)

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
#ifdef VAR_TIME    
       IF(master)THEN
          DO I=1,numnodes
             WRITE(lupri,'(A,I5,A,I8)')'Aux Basis Functions assign to rank ',I,' :',nbasisauxMPI(I)
          ENDDO
       ENDIF
#endif
    ELSE
       MynauxMPI = naux
    ENDIF
    CALL LSTIMER('RIMP2: MPI Basis ',TS2,TE2,LUPRI,FORCEPRINT)

    IF(wakeslaves)then 
       IF(MynauxMPI.GT.0)THEN
          !We wish to build
          !c_(alpha,ai) = (alpha|beta)^(-1/2) (beta|ai)
          !where alpha runs over the Aux basis functions allocated for this rank
          !and beta run over the full set of naux
          call mem_alloc(TMPAlphaBeta_minus_sqrt,MynauxMPI,naux)
          call RIMP2_buildTMPAlphaBeta_inv(TMPAlphaBeta_minus_sqrt,MynauxMPI,naux,&
               & nAtomsMPI,mynum,startAuxMPI,nAuxMPI,AlphaBeta_minus_sqrt,numnodes,nAtoms)
          call mem_dealloc(AlphaBeta_minus_sqrt)
       ELSE
          call mem_dealloc(AlphaBeta_minus_sqrt)          
       ENDIF
    ENDIF
    CALL LSTIMER('RIMP2: MPI AlphaBetaTmp ',TS2,TE2,LUPRI,FORCEPRINT)

    call get_currently_available_memory(MemInGBCollected)
    !maxsize = mem in number of floating point elements
    maxsize = NINT(MemInGBCollected*1.E9_realk) 

    AlphaCDAlloced = .FALSE.
    IF(MynauxMPI.GT.0)THEN
       !call mem_alloc(AlphaCD,naux,nvirt,nocc)
       !It is very annoying but I allocated AlphaCD inside 
       !II_get_RI_AlphaCD_3centerInt2 due to memory concerns
       !This Part of the Code is MPI/OpenMP parallel and AlphaCD will have the dimensions
       !(MynauxMPI,nvirt,nocc) 
       !nauxMPI is naux divided out on the nodes so roughly nauxMPI=naux/numnodes
       call II_get_RI_AlphaCD_3centerInt2(lupri,lupri,&
            & AlphaCD,mylsitem%setting,naux,nbasis,&
            & nvirt,nocc,MyMolecule%Cv,MyMolecule%Co,maxsize,mynum,numnodes)
       call mem_LeakTool_alloc(AlphaCD,LT_AlphaCD)
       AlphaCDAlloced = .TRUE.
    ENDIF
    CALL LSTIMER('RIMP2: AlphaCD ',TS2,TE2,LUPRI,FORCEPRINT)

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
          CALL LS_GETTIM(CPU1,WALL1)
          call MPI_ISEND(AlphaCD,COUNT,MPI_DOUBLE_PRECISION,Sender,TAG,comm,request,ierr)
          CALL LS_GETTIM(CPU2,WALL2)
          CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
          WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
!          call MPI_Request_free(request,ierr)
       ENDIF
       IF(MynauxMPI.GT.0)THEN 
          !consider 
          !c_(alpha,ai) = (alpha|beta)^(-1/2) (beta|ai)
          !so that you only need 1 3dim quantity      
          call mem_alloc(Calpha,MynauxMPI,nvirt,nocc)
          call mem_leaktool_alloc(Calpha,LT_Calpha)
          !Use own AlphaCD to obtain part of Calpha
          CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
          call RIMP2_buildOwnCalphaFromAlphaCD(nocc,nvirt,mynum,numnodes,&
               & natoms,MynauxMPI,nAtomsMPI,startAuxMPI,nAuxMPI,AlphaCD,&
               & Calpha,TMPAlphaBeta_minus_sqrt,naux)
          CALL LSTIMER('OwnCalpha ',TS3,TE3,LUPRI,FORCEPRINT)
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
       RoundRobin5 = .FALSE.
       RoundRobin6 = .FALSE.
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
                !all buffers are allocated and await to be deallocated .OR.
                !Memory Reduced version is used: 
                !once the memory have been recieved by the reciever the 
                !memory can be deallocated. 
                IF(CurrentWait(1).EQ.5)THEN
                   !I need to wait for AlphaCD5 to be received before I can deallocate
                   CALL LS_GETTIM(CPU1,WALL1)
                   call MPI_WAIT(request5,lsmpi_status,ierr)
                   CALL LS_GETTIM(CPU2,WALL2)
                   CPU_MPIWAIT = CPU_MPIWAIT + (CPU2-CPU1)
                   WALL_MPIWAIT = WALL_MPIWAIT + (WALL2-WALL1)

                   call mem_leaktool_dealloc(AlphaCD5,LT_AlphaCD5)
                   call mem_dealloc(AlphaCD5)
                ELSEIF(CurrentWait(1).EQ.6)THEN
                   !I need to wait for AlphaCD5 to be received before I can deallocate
                   CALL LS_GETTIM(CPU1,WALL1)
                   call MPI_WAIT(request6,lsmpi_status,ierr)
                   CALL LS_GETTIM(CPU2,WALL2)
                   CPU_MPIWAIT = CPU_MPIWAIT + (CPU2-CPU1)
                   WALL_MPIWAIT = WALL_MPIWAIT + (WALL2-WALL1)
                   call mem_leaktool_dealloc(AlphaCD6,LT_AlphaCD6)
                   call mem_dealloc(AlphaCD6)
                ENDIF
                IF(nAwaitDealloc.EQ.2)THEN
                   nAwaitDealloc = 1
                   CurrentWait(1) = CurrentWait(2)
                   CurrentWait(2) = 0 
                ELSE
                   nAwaitDealloc = 0
                   CurrentWait(2) = 0 
                   CurrentWait(1) = 0
                ENDIF
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
                !First time (node=1) RECV from the ISEND at line 1068, which sends AlphaCD
                !all other times RECV AlphaCD5
                CALL LS_GETTIM(CPU1,WALL1)
                call MPI_RECV(AlphaCD5,COUNT,MPI_DOUBLE_PRECISION,Receiver,TAG,comm,lsmpi_status,ierr)
                CALL LS_GETTIM(CPU2,WALL2)
                CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
                WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)

                IF(node.NE.numnodes-1)THEN
                   RoundRobin5 = .TRUE.
                   !SEND AlphaCD5 to sender
                   CALL LS_GETTIM(CPU1,WALL1)
                   call MPI_ISEND(AlphaCD5,COUNT,MPI_DOUBLE_PRECISION,Sender,TAG,comm,request5,ierr)
                   CALL LS_GETTIM(CPU2,WALL2)
                   CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
                   WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
                ELSE
                   !Last time (node=numnodes-1) Do not send because the sender is the 
                   !original owner of the block. Since I do not send, I do not need to MPI_WAIT
                   RoundRobin5 = .FALSE.
                ENDIF
             ELSEIF(useAlphaCD6)THEN
                CALL LS_GETTIM(CPU1,WALL1)                
                call MPI_RECV(AlphaCD6,COUNT,MPI_DOUBLE_PRECISION,Receiver,TAG,comm,lsmpi_status,ierr)
                CALL LS_GETTIM(CPU2,WALL2)
                CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
                WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
                IF(node.NE.numnodes-1)THEN
                   RoundRobin6 = .TRUE.
                   !SEND AlphaCD6 to sender
                   CALL LS_GETTIM(CPU1,WALL1)                
                   call MPI_ISEND(AlphaCD6,COUNT,MPI_DOUBLE_PRECISION,Sender,TAG,comm,request6,ierr)
                   CALL LS_GETTIM(CPU2,WALL2)
                   CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
                   WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
                ELSE
                   !Last time (node=numnodes-1) Do not send because the sender is the 
                   !original owner of the block. Since I do not send, I do not need to MPI_WAIT
                   RoundRobin6 = .FALSE.
                ENDIF
             ENDIF
             IF(AlphaCDAlloced)THEN
                AlphaCD_Deallocate = .FALSE.
                IF(MemoryReduced)THEN 
                   !Wait for alphaCD have been received and deallocate it 
                   !to avoid having 4 (nocc,nvirt,nAuxMPI) - reduce to 3.
                   IF(MynauxMPI.GT.0)THEN
                      CALL LS_GETTIM(CPU1,WALL1)
                      call MPI_WAIT(request,lsmpi_status,ierr)
                      CALL LS_GETTIM(CPU2,WALL2)
                      CPU_MPIWAIT = CPU_MPIWAIT + (CPU2-CPU1)
                      WALL_MPIWAIT = WALL_MPIWAIT + (WALL2-WALL1)
                   ENDIF
                   AlphaCD_Deallocate = .TRUE.
                ELSE
                   IF(MynauxMPI.GT.0)THEN 
                      call MPI_TEST(request,TransferCompleted,lsmpi_status,ierr)
                      AlphaCD_Deallocate = TransferCompleted
                   ENDIF
                ENDIF
                IF(AlphaCD_Deallocate)THEN
                   call mem_leaktool_dealloc(alphaCD,LT_alphaCD)
                   call mem_dealloc(AlphaCD) !no longer need this
                   AlphaCDAlloced = .FALSE.      
                ENDIF
             ENDIF
             IF(MynauxMPI.GT.0)THEN
                !Step 2: Obtain part of Calpha from this contribution
                CALL LSTIMER('START ',TS3,TE3,LUPRI,FORCEPRINT)
                IF(useAlphaCD5)THEN
                   call RIMP2_buildCalphaContFromAlphaCD(nocc,nvirt,myOriginalRank,numnodes,natoms,&
                        & OriginalRanknauxMPI,MynauxMPI,nAtomsMPI,startAuxMPI,nAuxMPI,&
                        & AlphaCD5,Calpha,TMPAlphaBeta_minus_sqrt,naux)
                ELSEIF(useAlphaCD6)THEN
                   call RIMP2_buildCalphaContFromAlphaCD(nocc,nvirt,myOriginalRank,numnodes,natoms,&
                        & OriginalRanknauxMPI,MynauxMPI,nAtomsMPI,startAuxMPI,nAuxMPI,&
                        & AlphaCD6,Calpha,TMPAlphaBeta_minus_sqrt,naux)
                ENDIF
                CALL LSTIMER('CalphaOther ',TS3,TE3,LUPRI,FORCEPRINT)
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
                !I can deallocate directly I do not need to MPI_WAIT since I did not send
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
          call mem_dealloc(TMPAlphaBeta_minus_sqrt)
       ENDIF
       NBA = MynauxMPI
       IF(nAwaitDealloc.NE.0)THEN
          do iAwaitDealloc=1,nAwaitDealloc
             IF(CurrentWait(iAwaitDealloc).EQ.5)THEN
                IF(RoundRobin5)THEN
                   CALL LS_GETTIM(CPU1,WALL1)                
                   call MPI_WAIT(request5,lsmpi_status,ierr)
                   CALL LS_GETTIM(CPU2,WALL2)
                   CPU_MPIWAIT = CPU_MPIWAIT + (CPU2-CPU1)
                   WALL_MPIWAIT = WALL_MPIWAIT + (WALL2-WALL1)
                ENDIF
                call mem_leaktool_dealloc(AlphaCD5,LT_AlphaCD5)
                call mem_dealloc(AlphaCD5)
             ELSEIF(CurrentWait(iAwaitDealloc).EQ.6)THEN
                IF(RoundRobin6)THEN
                   CALL LS_GETTIM(CPU1,WALL1)                
                   call MPI_WAIT(request6,lsmpi_status,ierr)
                   CALL LS_GETTIM(CPU2,WALL2)
                   CPU_MPIWAIT = CPU_MPIWAIT + (CPU2-CPU1)
                   WALL_MPIWAIT = WALL_MPIWAIT + (WALL2-WALL1)
                ENDIF
                call mem_leaktool_dealloc(AlphaCD6,LT_AlphaCD6)
                call mem_dealloc(AlphaCD6)
             ENDIF
          enddo
       ENDIF
       IF(AlphaCDAlloced)THEN
          IF(MynauxMPI.GT.0)THEN
             CALL LS_GETTIM(CPU1,WALL1)
             call MPI_WAIT(request,lsmpi_status,ierr)
             CALL LS_GETTIM(CPU2,WALL2)
             CPU_MPIWAIT = CPU_MPIWAIT + (CPU2-CPU1)
             WALL_MPIWAIT = WALL_MPIWAIT + (WALL2-WALL1)
             call mem_leaktool_dealloc(alphaCD,LT_alphaCD)
             call mem_dealloc(AlphaCD) !no longer need this
             AlphaCDAlloced = .FALSE.      
          ENDIF
       ENDIF
#endif
    ELSE
       call mem_alloc(Calpha,naux,nvirt,nocc)
       call mem_leaktool_alloc(Calpha,LT_Calpha)
       !Calpha(naux,nvirt,nocc) = AlphaBeta_minus_sqrt(naux,naux)*AlphaCD(naux,nvirt,nocc)
       call DGEMM('N','N',naux,nvirt*nocc,naux,1.0E0_realk,AlphaBeta_minus_sqrt,&
            & naux,AlphaCD,naux,0.0E0_realk,Calpha,naux)
       call mem_dealloc(AlphaBeta_minus_sqrt)
       NBA = naux
       call mem_leaktool_dealloc(alphaCD,LT_alphaCD)
       call mem_dealloc(AlphaCD) !no longer need this
       AlphaCDAlloced = .FALSE.      
    ENDIF
    CALL LSTIMER('RIMP2: Calpha ',TS2,TE2,LUPRI,FORCEPRINT)


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

    IF(DECinfo%DECrestart)THEN
     !CHECK IF THERE ARE ENERGY CONTRIBUTIONS AVAILABLE
     INQUIRE(FILE='FULLRIMP2.restart',EXIST=file_exists)
     IF(file_exists)THEN
      IF(master)THEN
       WRITE(DECinfo%output,*)'Restart of Full molecular RIMP2 calculation:'
      ENDIF
      restart_lun = -1  !initialization
      call lsopen(restart_lun,'FULLRIMP2.restart','OLD','FORMATTED')
      rewind restart_lun
      read(restart_lun,'(I9)') noccJstart
      IF(noccJstart.EQ.nocc)THEN
         IF(master)WRITE(DECinfo%output,*)'All energies is on file'
         noccJstart = nocc+1
         read(restart_lun,'(F28.16)') rimp2_energy
      ELSEIF(noccJstart.GT.nocc.OR.noccJstart.LT.1)THEN
       IF(master)THEN
        WRITE(DECinfo%output,*)'RIMP2 restart error first integer is wrong. Read:',noccJstart
       ENDIF
       call lsquit('RIMP2 restart error first integer is wrong')             
      ELSE
       noccJstart = noccJstart + 1
       read(restart_lun,'(F28.16)') rimp2_energy
      ENDIF
      call lsclose(restart_lun,'KEEP')
     ELSE
      noccJstart=1
      rimp2_energy = 0.0E0_realk
     ENDIF
    ELSE
     noccJstart=1
     rimp2_energy = 0.0E0_realk
    ENDIF

    IF(Wakeslaves)THEN
#ifdef VAR_MPI
       NotMatSet = .TRUE.
       call mem_alloc(Amat,nvirt,nvirt)
       call mem_alloc(Bmat,nvirt,nvirt)
       do J=noccJstart,nocc
          do I=1,nocc
             epsIJ = EpsOcc(I) + EpsOcc(J)
             IF(MynauxMPI.GT.0)THEN 
                !A(nvirt,nvirt)
                CALL CalcAmat(nocc,nvirt,NBA,Calpha,Amat,I,J)
                CALL CalcBmat(nvirt,EpsIJ,EpsVirt,Amat,Bmat)
             ELSE
                IF(NotMatSet)THEN
                   NotMatSet = .FALSE.
                   call ls_dzero(Amat,nvirt*nvirt) 
                   call ls_dzero(Bmat,nvirt*nvirt) 
                ENDIF
             ENDIF
             !The barrier is mostly here to detect the time spent 
             !waiting vs the time spent in communication. 
             CALL LS_GETTIM(CPU1,WALL1)
             call lsmpi_barrier(comm)
             CALL LS_GETTIM(CPU2,WALL2)
             CPU_MPIWAIT = CPU_MPIWAIT + (CPU2-CPU1)
             WALL_MPIWAIT = WALL_MPIWAIT + (WALL2-WALL1)
             !Reduce A,B
             call lsmpi_reduction(Amat,nvirt,nvirt,infpar%master,comm)
             call lsmpi_reduction(Bmat,nvirt,nvirt,infpar%master,comm)
             CALL LS_GETTIM(CPU1,WALL1)
             CPU_MPICOMM = CPU_MPICOMM + (CPU1-CPU2)
             WALL_MPICOMM = WALL_MPICOMM + (WALL1-WALL2)             
             IF(master)THEN
                call MP2_EnergyContribution(nvirt,Amat,Bmat,rimp2_energy)
             ENDIF
          enddo
          !Write Restart File 
          restart_lun = -1  !initialization
          call lsopen(restart_lun,'FULLRIMP2.restart','UNKNOWN','FORMATTED')
          rewind restart_lun
          write(restart_lun,'(I9)') J
          write(restart_lun,'(F28.16)') rimp2_energy
          call lsclose(restart_lun,'KEEP')
       enddo
       call mem_dealloc(Amat)
       call mem_dealloc(Bmat)
       call mem_dealloc(nbasisauxMPI)
       call mem_dealloc(startAuxMPI)
       call mem_dealloc(nAtomsMPI)
       call mem_dealloc(nAuxMPI)
       IF(MynauxMPI.GT.0)THEN 
          call mem_leaktool_dealloc(Calpha,LT_Calpha)
          call mem_dealloc(Calpha)              
       ENDIF
#endif
    ELSE
       !Energy = sum_{AIBJ} (AI|BJ)_N*[ 2(AI|BJ)_N - (BI|AJ)_N ]/(epsI+epsJ-epsA-epsB)
       call RIMP2_CalcOwnEnergyContribution(nocc,nvirt,EpsOcc,EpsVirt,&
            & NBA,Calpha,rimp2_energy)
       call mem_leaktool_dealloc(Calpha,LT_Calpha)
       call mem_dealloc(Calpha)              
    ENDIF
    CALL LSTIMER('RIMP2: EnergyCont ',TS2,TE2,LUPRI,FORCEPRINT)
    call mem_leaktool_dealloc(EpsOcc,LT_Eps)
    call mem_dealloc(EpsOcc)
    call mem_leaktool_dealloc(EpsVirt,LT_Eps)
    call mem_dealloc(EpsVirt)
    IF(MASTER)THEN
       write(lupri,*)  'RIMP2 CORRELATION ENERGY = ', rimp2_energy
       write(*,'(1X,a,f20.10)') 'RIMP2 CORRELATION ENERGY = ', rimp2_energy
    ENDIF
    write(lupri,*)  'LEAK TOOL STATISTICS IN full_canonical_rimp2'
    call LeakTools_stat_mem(lupri)
    CALL LSTIMER('RIMP2: Finalize ',TS2,TE2,LUPRI,FORCEPRINT)
    CALL LSTIMER('FULL RIMP2 ',TS,TE,DECINFO%OUTPUT)
#ifdef VAR_MPI
    write(lupri,*)'Overall Time spent in MPI Communication and MPI Wait for rank=',infpar%mynum
    CALL ls_TIMTXT('>>>  WALL Time used MPI Communication inc. some Wait',WALL_MPICOMM,lupri)
    CALL ls_TIMTXT('>>>  CPU Time used MPI Communication inc. some Wait',CPU_MPICOMM,lupri)
    CALL ls_TIMTXT('>>>  WALL Time used in MPI Wait',WALL_MPIWAIT,lupri)
    CALL ls_TIMTXT('>>>  CPU Time used in MPI Wait',CPU_MPIWAIT,lupri)
#endif    
    write(lupri,*) ' '
  end subroutine full_canonical_rimp2

  subroutine CalcAmat(nocc,nvirt,NBA,Calpha,Amat,I,J)
    implicit none
    integer,intent(in) :: nocc,nvirt,NBA,I,J
    real(realk),intent(in) :: Calpha(NBA,nvirt,nocc)
    real(realk),intent(inout) :: Amat(nvirt,nvirt)
    !
    integer :: A,B,ALPHA
    real(realk) :: TMP
    !permutational symmetry
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
    !$OMP PRIVATE(B,A,ALPHA,TMP) SHARED(nocc,nvirt,NBA,Calpha,Amat,I,J)
    do B=1,nvirt        
       do A=1,nvirt      
          TMP = 0.0E0_realk
          DO ALPHA = 1,NBA
             TMP = TMP + Calpha(ALPHA,A,I)*Calpha(ALPHA,B,J)
          ENDDO
          Amat(A,B) = TMP
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine CalcAmat

  subroutine CalcBmat(nvirt,EpsIJ,EpsVirt,Amat,Bmat)
    implicit none
    integer,intent(in) :: nvirt
    real(realk),intent(in) :: EpsIJ,Amat(nvirt,nvirt),EpsVirt(nvirt)
    real(realk),intent(inout) :: Bmat(nvirt,nvirt)
    !
    integer :: A,B
    !permutational symmetry ?
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
    !$OMP PRIVATE(B,A) SHARED(nvirt,Amat,Bmat,epsIJ,EpsVirt)
    do B=1,nvirt        
       do A=1,nvirt
          Bmat(A,B) = (2.0E0_realk*Amat(A,B) - Amat(B,A))/(epsIJ-EpsVirt(A)-EpsVirt(B))
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine CalcBmat

  !FIXME USE A DOT - is that faster?
  subroutine MP2_EnergyContribution(nvirt,Amat,Bmat,rimp2_energy)
    implicit none
    integer,intent(in) :: nvirt
    real(realk),intent(in) :: Amat(nvirt*nvirt),Bmat(nvirt*nvirt)
    real(realk),intent(inout) :: rimp2_energy
    !
    integer :: A
    real(realk) :: TMP

    TMP = 0.0E0_realk
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP PRIVATE(A) REDUCTION(+:TMP) SHARED(nvirt,Amat,Bmat)
    DO A=1,nvirt*nvirt
       TMP = TMP + Amat(A)*Bmat(A)      
    ENDDO
    !$OMP END PARALLEL DO
    rimp2_energy = rimp2_energy + TMP
  end subroutine MP2_EnergyContribution

! Calculate canonical MP2 energy
subroutine RIMP2_CalcOwnEnergyContribution(nocc,nvirt,EpsOcc,EpsVirt,NBA,Calpha,rimp2_energy)
  implicit none
  integer,intent(in) :: nocc,nvirt,NBA
  real(realk),intent(in) :: EpsOcc(nocc),EpsVirt(nvirt)
  real(realk),intent(in) :: Calpha(NBA,nvirt,nocc)
  real(realk),intent(inout) :: rimp2_energy
  !
  integer :: J,B,A,I,ALPHA
  real(realk) :: eps,gmoAIBJ,gmoBIAJ,TMP,epsIJB
  Tmp = 0.0E0_realk
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
  !$OMP PRIVATE(J,B,A,I,eps,gmoAIBJ,gmoBIAJ,epsIJB,ALPHA) REDUCTION(+:TMP) &
  !$OMP SHARED(nocc,nvirt,EpsOcc,EpsVirt,NBA,Calpha)
  do J=1,nocc
     do B=1,nvirt        
        epsIJB = EpsOcc(J) + EpsOcc(J) - EpsVirt(B)
        !==================================================================
        ! I=J AND A=B
        !==================================================================
        ! Difference in orbital energies: eps(I) + eps(J) - eps(A) - eps(B)
        eps = epsIJB - EpsVirt(B)
        gmoAIBJ = 0.0E0_realk
        DO ALPHA = 1,NBA
           gmoAIBJ = gmoAIBJ + Calpha(ALPHA,B,J)*Calpha(ALPHA,B,J)
        ENDDO
        !Energy = sum_{AIBJ} (AI|BJ)*[ 2(AI|BJ) - (BI|AJ) ]/(epsI + epsJ - epsA - epsB)
        Tmp = Tmp + gmoAIBJ*gmoAIBJ/eps
        !==================================================================
        ! I=J AND A>B
        !==================================================================
        do A=B+1,nvirt
           ! Difference in orbital energies: eps(I) + eps(J) - eps(A) - eps(B)
           eps = epsIJB - EpsVirt(A)
           gmoAIBJ = 0.0E0_realk
           DO ALPHA = 1,NBA
              gmoAIBJ = gmoAIBJ + Calpha(ALPHA,A,J)*Calpha(ALPHA,B,J)
           ENDDO
           !Energy = sum_{AIBJ} (AI|BJ)*[ 2(AI|BJ) - (BI|AJ) ]/(epsI + epsJ - epsA - epsB)
           Tmp = Tmp + 2.0E0_realk*gmoAIBJ*gmoAIBJ/eps
        end do
        do I=J+1,nocc
           epsIJB = EpsOcc(I) + EpsOcc(J) - EpsVirt(B)
           !==================================================================
           ! I>J AND A=B
           !==================================================================
           ! Difference in orbital energies: eps(I) + eps(J) - eps(A) - eps(B)
           eps = epsIJB - EpsVirt(B)
           gmoAIBJ = 0.0E0_realk
           DO ALPHA = 1,NBA
              gmoAIBJ = gmoAIBJ + Calpha(ALPHA,B,I)*Calpha(ALPHA,B,J)
           ENDDO
           !Energy = sum_{AIBJ} (AI|BJ)*[ 2(AI|BJ) - (BI|AJ) ]/(epsI + epsJ - epsA - epsB)
           Tmp = Tmp + 2.0E0_realk*gmoAIBJ*gmoAIBJ/eps
           !==================================================================
           ! I>J AND A>B
           !==================================================================
           do A=B+1,nvirt
              ! Difference in orbital energies: eps(I) + eps(J) - eps(A) - eps(B)
              eps = epsIJB - EpsVirt(A)
              gmoAIBJ = 0.0E0_realk
              DO ALPHA = 1,NBA
                 gmoAIBJ = gmoAIBJ + Calpha(ALPHA,A,I)*Calpha(ALPHA,B,J)
              ENDDO
              gmoBIAJ = 0.0E0_realk                   
              DO ALPHA = 1,NBA
                 gmoBIAJ = gmoBIAJ + Calpha(ALPHA,B,I)*Calpha(ALPHA,A,J)
              ENDDO
              !Energy = sum_{AIBJ} (AI|BJ)*[ 2(AI|BJ) - (BI|AJ) ]/(epsI + epsJ - epsA - epsB)
!              Tmp = Tmp + 2.0E0_realk*gmoAIBJ*(2E0_realk*gmoAIBJ-gmoBIAJ)/eps
!              Tmp = Tmp + 2.0E0_realk*gmoBIAJ*(2E0_realk*gmoBIAJ-gmoAIBJ)/eps
              Tmp = Tmp + 4.0E0_realk*(gmoAIBJ*gmoAIBJ-gmoAIBJ*gmoBIAJ+gmoBIAJ*gmoBIAJ)/eps
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
     & Calpha,Calpha2,NBA2,rimp2_energy)
  implicit none
  integer,intent(in) :: nocc,nvirt,NBA,NBA2
  real(realk),intent(in) :: EpsOcc(nocc),EpsVirt(nvirt)
  real(realk),intent(in) :: Calpha(NBA,nvirt,nocc)
  real(realk),intent(in) :: Calpha2(NBA2,nvirt,nocc)
  real(realk),intent(inout) :: rimp2_energy
  !
  integer :: J,B,A,I,ALPHA
  real(realk) :: eps,gmoAIBJ,gmoBIAJ,TMP,epsIJB,gmoAIBJ2,gmoBIAJ2
  Tmp = 0.0E0_realk
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
  !$OMP PRIVATE(J,B,A,I,eps,gmoAIBJ,gmoBIAJ,epsIJB,gmoAIBJ2,ALPHA,gmoBIAJ2) &
  !$OMP REDUCTION(+:TMP) &
  !$OMP SHARED(nocc,nvirt,EpsOcc,EpsVirt,NBA,Calpha,NBA2,Calpha2)
  do J=1,nocc
     do B=1,nvirt
        epsIJB = EpsOcc(J) + EpsOcc(J) - EpsVirt(B)
        !================================================================
        ! I = J, A = B       
        !================================================================
        ! Difference in orbital energies: eps(I)+eps(J)-eps(A)-eps(B)
        eps = epsIJB - EpsVirt(B)
        gmoAIBJ2 = 0.0E0_realk
        DO ALPHA = 1,NBA2
           gmoAIBJ2 = gmoAIBJ2 + Calpha2(ALPHA,B,J)*Calpha2(ALPHA,B,J)
        ENDDO
        gmoAIBJ = 0.0E0_realk
        DO ALPHA = 1,NBA
           gmoAIBJ = gmoAIBJ + Calpha(ALPHA,B,J)*Calpha(ALPHA,B,J)
        ENDDO
        !Energy = sum_{AIBJ} (AI|BJ)*[ 2(AI|BJ) - (BI|AJ) ]/(epsI + epsJ - epsA - epsB)
        Tmp = Tmp + gmoAIBJ2*gmoAIBJ/eps
        !================================================================
        ! I = J, A > B       
        !================================================================
        do A=B+1,nvirt
           ! Difference in orbital energies: eps(I)+eps(J)-eps(A)-eps(B)
           eps = epsIJB - EpsVirt(A)
           gmoAIBJ2 = 0.0E0_realk
           DO ALPHA = 1,NBA2
              gmoAIBJ2 = gmoAIBJ2 + Calpha2(ALPHA,A,J)*Calpha2(ALPHA,B,J)
           ENDDO
           gmoAIBJ = 0.0E0_realk
           DO ALPHA = 1,NBA
              gmoAIBJ = gmoAIBJ + Calpha(ALPHA,A,J)*Calpha(ALPHA,B,J)
           ENDDO
           !Energy = sum_{AIBJ} (AI|BJ)*[ 2(AI|BJ) - (BI|AJ) ]/(epsI + epsJ - epsA - epsB)
           Tmp = Tmp + 2.0E0_realk*gmoAIBJ2*gmoAIBJ/eps
        end do
        !================================================================
        ! I > J, A = B       
        !================================================================
        do I=J+1,nocc
           epsIJB = EpsOcc(I) + EpsOcc(J) - EpsVirt(B)
           ! Difference in orbital energies: eps(I)+eps(J)-eps(A)-eps(B)
           eps = epsIJB - EpsVirt(B)
           gmoAIBJ2 = 0.0E0_realk
           DO ALPHA = 1,NBA2
              gmoAIBJ2 = gmoAIBJ2 + Calpha2(ALPHA,B,I)*Calpha2(ALPHA,B,J)
           ENDDO
           gmoAIBJ = 0.0E0_realk
           DO ALPHA = 1,NBA
              gmoAIBJ = gmoAIBJ + Calpha(ALPHA,B,I)*Calpha(ALPHA,B,J)
           ENDDO
           !Energy = sum_{AIBJ} (AI|BJ)*[ 2(AI|BJ) - (BI|AJ) ]/(epsI + epsJ - epsA - epsB)
           Tmp = Tmp + 2.0E0_realk*gmoAIBJ2*gmoAIBJ/eps
           !================================================================
           ! I > J, A > B       
           !================================================================
           do A=B+1,nvirt
              ! Difference in orbital energies: eps(I)+eps(J)-eps(A)-eps(B)
              eps = epsIJB - EpsVirt(A)
              gmoAIBJ2 = 0.0E0_realk
              gmoBIAJ2 = 0.0E0_realk
              DO ALPHA = 1,NBA2
                 gmoAIBJ2 = gmoAIBJ2 + Calpha2(ALPHA,A,I)*Calpha2(ALPHA,B,J)
                 gmoBIAJ2 = gmoBIAJ2 + Calpha2(ALPHA,B,I)*Calpha2(ALPHA,A,J)
              ENDDO
              gmoAIBJ = 0.0E0_realk
              gmoBIAJ = 0.0E0_realk                   
              DO ALPHA = 1,NBA
                 gmoAIBJ = gmoAIBJ + Calpha(ALPHA,A,I)*Calpha(ALPHA,B,J)
                 gmoBIAJ = gmoBIAJ + Calpha(ALPHA,B,I)*Calpha(ALPHA,A,J)
              ENDDO
              !Energy = sum_{AIBJ} (AI|BJ)*[ 2(AI|BJ) - (BI|AJ) ]/(epsI + epsJ - epsA - epsB)
              Tmp = Tmp + 2.0E0_realk*gmoAIBJ2*(2E0_realk*gmoAIBJ-gmoBIAJ)/eps
              Tmp = Tmp + 2.0E0_realk*gmoBIAJ2*(2E0_realk*gmoBIAJ-gmoAIBJ)/eps
           end do
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  rimp2_energy = rimp2_energy + Tmp
end subroutine RIMP2_CalcEnergyContribution

!> \brief Memory check for full_canonical_mp2 subroutine
!> \author Thomas Kjaergaard
!> \date October 2014
subroutine full_canonical_mp2_memory_check(nbasis,nvirt,MinAtomic)
  implicit none
  integer,intent(in) :: nbasis,nvirt,MinAtomic
  real(realk) :: MemRequired,GB,MemStep1,MemStep2
  GB = 1.0E-9_realk
  ! Check that arrays fit in memory (hard-coded)
  MemRequired = real(MinAtomic*MinAtomic*nbasis*(nbasis+nvirt))
  MemRequired = MemRequired*realk*GB  
  if(MemRequired > 0.80E0_realk*DECinfo%memory) then
     call lsquit('full_canonical_mp2: Memory exceeded! ',-1)
  end if
end subroutine full_canonical_mp2_memory_check

!> \brief Calculate canonical MP2 energy for full molecular system
!> \author Thomas Kjaergaard
!> \date October 2014
subroutine full_canonical_mp2(MyMolecule,MyLsitem,mp2_energy)
  implicit none
  !> Full molecule info
  type(fullmolecule), intent(inout) :: MyMolecule
  !> Lsitem structure
  type(lsitem), intent(inout) :: mylsitem
  !> Canonical MP2 correlation energy
  real(realk),intent(inout) :: mp2_energy    
  !
  integer :: nbasis,nocc,nvirt,naux,noccfull,mynum,numnodes
  logical :: master,wakeslaves
  real(realk),pointer :: EpsOcc(:),EpsVirt(:)
  integer :: J,I,node,natoms,A,lupri,restart_lun
  logical :: MessageRecieved,MessageRecievedW,FORCEPRINT,file_exists
  real(realk) :: TS,TE,TS2,TE2,TS3,TE3,CPU1,CPU2,WALL1,WALL2,CPU_MPICOMM,WALL_MPICOMM
  real(realk) :: CPU_MPIWAIT,WALL_MPIWAIT,MemInGBCollected,epsIJ,tmp_mp2_energy
  real(realk) :: CPU3,CPU4,WALL3,WALL4,CPU_AOINT,WALL_AOINT,tmp_mp2_energy2
  real(realk) :: CPU_AOTOMO,WALL_AOTOMO,CPU_ECONT,WALL_ECONT
  real(realk),pointer :: Amat(:,:),Bmat(:,:),tmp1(:,:),tmp2(:,:),tmp3(:,:),tmp4(:,:)
  real(realk),pointer :: tmp5(:,:),tmp6(:,:),tmp7(:,:),CoBatchA(:,:),CoBatchB(:,:)
  real(realk),pointer :: CoBI(:,:),CoBJ(:,:),tmp62(:,:),tmp72(:,:),CoI(:,:),CoJ(:,:)
  real(realk),pointer :: VOVO(:,:),VGVO(:,:),CvA(:,:),CoIG(:,:)
!
  type(DecAObatchinfo),pointer :: AOGammabatchinfo(:)
  type(DecAObatchinfo),pointer :: AOAlphabatchinfo(:)
  type(batchtoorb), pointer :: batch2orbAlpha(:)
  type(batchtoorb), pointer :: batch2orbGamma(:)
  integer, pointer :: orb2batchAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
  integer, pointer :: orb2batchGamma(:), batchsizeGamma(:), batchindexGamma(:)
  integer, pointer :: batchdimAlpha(:), batchdimGamma(:)
  integer :: MaxAllowedDimAlpha,MaxAllowedDimGamma,MaxActualDimAlpha,MaxActualDimGamma
  integer :: nOccBatchDimImax,nOccBatchDimJmax,K,iorb,idx,nbatchesAlpha,nbatchesGamma,nb
  integer :: iAO,nAObatches,AOGammaStart,AOGammaEnd,AOAlphaStart,AOAlphaEnd,iprint,M,N,iB
  integer :: MinAObatch,gammaB,alphaB,nOccBatchesI,nOccBatchesJ,jB,nOccBatchDimI,nOccBatchDimJ
  integer :: nbatchesGammaRestart,nbatchesAlphaRestart,dimGamma,GammaStart,GammaEnd,dimAlpha
  integer :: AlphaStart,AlphaEnd,B,noccRestartI,noccRestartJ,nJobs,startjB 
  integer :: nOccBatchesIrestart,noccIstart
  logical :: MoTrans, NoSymmetry,SameMol,JobDone,JobInfo1Free,FullRHS,doscreen,NotAllMessagesRecieved
  logical :: PermutationalSymmetryIJ
  logical,pointer :: JobsCompleted(:,:)
  integer(kind=8) :: maxsize
  TYPE(DECscreenITEM)   :: DecScreen
  Character            :: intSpec(5)
  integer(kind=ls_mpik)  :: nMPI,TAG,IERR,request,Receiver,sender,comm,TAG1,TAG2
  integer(kind=ls_mpik)  :: request1,request2,masterrank,senderID
#ifdef VAR_MPI
  integer(kind=ls_mpik)  :: mpistatus(MPI_STATUS_SIZE) 
#endif
  integer(kind=4) :: JobInfo1(2)
!  Character(80)        :: FilenameCS,FilenamePS

#ifdef VAR_TIME    
    FORCEPRINT = .TRUE.
#endif

  CALL LSTIMER('START ',TS,TE,DECINFO%OUTPUT)
  mp2_energy = 0.0E0_realk
  
  CPU_AOINT = 0.0E0_realk
  WALL_AOINT = 0.0E0_realk
  CPU_AOTOMO = 0.0E0_realk
  WALL_AOTOMO = 0.0E0_realk
  CPU_ECONT = 0.0E0_realk
  WALL_ECONT = 0.0E0_realk

  CPU_MPICOMM = 0.0E0_realk
  WALL_MPICOMM = 0.0E0_realk
  CPU_MPIWAIT = 0.0E0_realk
  WALL_MPIWAIT = 0.0E0_realk
  TAG = 1411; TAG1 = 1412; TAG2 = 1413
  
  !sanity check
  if(.NOT.DECinfo%use_canonical) then
     call lsquit('Error: full_canonical_mp2 require canonical Orbitals',-1)
  endif
  ! Init stuff
  ! **********
  nbasis = MyMolecule%nbasis
  nb = nbasis
  nocc   = MyMolecule%nocc
  nvirt  = MyMolecule%nunocc
  noccfull = nocc
  nAtoms = MyMolecule%nAtoms
  LUPRI = DECinfo%output

#ifdef VAR_MPI
  comm = MPI_COMM_LSDALTON
  master= (infpar%mynum == infpar%master)
  mynum = infpar%mynum
  numnodes = infpar%nodtot
  wakeslaves = infpar%nodtot.GT.1
  IF(.NOT.master)LUPRI = 6 !standard Output
#else
  ! If MPI is not used, consider the single node to be "master"
  master=.true.
  mynum = 0
  numnodes = 1
  wakeslaves = .false.
#endif

  ! Set integral info
  ! *****************
  INTSPEC(1)='R' !R = Regular Basis set on the 1th center 
  INTSPEC(2)='R' !R = Regular Basis set on the 2th center 
  INTSPEC(3)='R' !R = Regular Basis set on the 3th center 
  INTSPEC(4)='R' !R = Regular Basis set on the 4th center 
  INTSPEC(5)='C' !C = Coulomb operator

  !determine MinAObatch: the minimum allowed AObatch size + number of AO batches
  IF(DECinfo%useIchor)THEN
     !Determine the minimum allowed AObatch size MinAObatch
     !In case of pure Helium atoms in cc-pVDZ ((4s,1p) -> [2s,1p]) MinAObatch = 3 (Px,Py,Pz)
     !In case of pure Carbon atoms in cc-pVDZ ((9s,4p,1d) -> [3s,2p,1d]) MinAObatch=6 (the 2*(Px,Py,Pz))
     !In case of pure Carbon atoms in 6-31G   ((10s,4p) -> [3s,2p]) MinAObatch = 3 (Px,Py,Pz) 
     !'R'  !Specifies that it is the Regular AO basis that should be batched
     iAO = 4 !the center that the batching should occur on (they are all the same in this case)  
     call determine_MinimumAllowedAObatchSize(MyLsItem%setting,iAO,'R',MinAObatch)
     iprint = 0           !print level for Ichor Integral code
     MoTrans = .FALSE.    !Do not transform to MO basis! 
     NoSymmetry = .FALSE. !Use Permutational Symmetry! 
     SameMol = .TRUE.     !Same molecule on all centers of the 4 center 2 electron integral
     !Determine the full number of AO batches - not to be confused with the batches of AOs
     !Required by the MAIN_ICHORERI_DRIVER unless all four dimensions are batched 
     iAO = 1
     call determine_Ichor_nAObatches(mylsitem%setting,iAO,'R',nAObatches,DECinfo%output)
     nullify(AOGammabatchinfo)
     nullify(AOalphabatchinfo)    
  ELSE
     ! The smallest possible AO batch depends on the basis set
     ! (More precisely, if all batches are made as small as possible, then the
     !  call below determines the largest of these small batches).
     call determine_maxBatchOrbitalsize(DECinfo%output,mylsitem%setting,MinAObatch,'R')
     nullify(orb2batchAlpha)
     nullify(batchdimAlpha)
     nullify(batchsizeAlpha)
     nullify(batch2orbAlpha)
     nullify(batchindexAlpha)
     nullify(orb2batchGamma)
     nullify(batchdimGamma)
     nullify(batchsizeGamma)
     nullify(batch2orbGamma)
     nullify(batchindexGamma)
  ENDIF

  doscreen = mylsitem%setting%scheme%cs_screen.OR.&
       & mylsitem%setting%scheme%ps_screen
  
  ! ***************************************************************************************
  !Get optimal values of: MaxAllowedDimAlpha,MaxAllowedDimGamma,nOccBatchDimImax,nOccBatchDimJmax
  ! ***************************************************************************************

  call get_optimal_batch_sizes_for_canonical_mp2(MinAObatch,nbasis,nocc,nvirt,&
       & MaxAllowedDimAlpha,MaxAllowedDimGamma,nOccBatchDimImax,nOccBatchDimJmax,&
       & numnodes,DECinfo%output)

  write(DECinfo%output,*)'nbasis',nbasis
  write(DECinfo%output,*)'nocc  ',nocc
  write(DECinfo%output,*)'nvirt ',nvirt
  write(DECinfo%output,*)'nOccBatchDimImax  ',nOccBatchDimImax
  write(DECinfo%output,*)'nOccBatchDimJmax  ',nOccBatchDimJmax
  write(DECinfo%output,*)'MinAObatch        ',MinAObatch
  write(DECinfo%output,*)'MaxAllowedDimAlpha',MaxAllowedDimAlpha
  write(DECinfo%output,*)'MaxAllowedDimGamma',MaxAllowedDimGamma
  
  nOccBatchesI = nOcc/nOccBatchDimImax
  IF(MOD(nOcc,nOccBatchDimImax).NE.0)nOccBatchesI = nOccBatchesI + 1  
  nOccBatchesJ = nOcc/nOccBatchDimJmax
  IF(MOD(nOcc,nOccBatchDimJmax).NE.0)nOccBatchesJ = nOccBatchesJ + 1

  PermutationalSymmetryIJ = .FALSE.
  IF(nOccbatchesI.EQ.nOccbatchesJ)THEN
     PermutationalSymmetryIJ = .TRUE.
     WRITE(DECinfo%output,*)'Permutational Symmetry exploited'
  ENDIF
  write(DECinfo%output,*)'nOccBatchesI',nOccBatchesI
  write(DECinfo%output,*)'nOccBatchesJ',nOccBatchesJ

#ifdef VAR_MPI 
  ! Master starts up slave
  StartUpSlaves: if(wakeslaves .and. master) then
     ! Wake up slaves to do the job: slaves awoken up with (RIMP2FULL)
     ! and call full_canonical_rimp2_slave which communicate info 
     ! then calls full_canonical_rimp2.
     CALL LS_GETTIM(CPU1,WALL1)
     call ls_mpibcast(CANONMP2FULL,infpar%master,comm)
     ! Communicate fragment information to slaves
     call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,comm)
     call mpicopy_lsitem(MyLsitem,comm)
     call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,comm)
     call mpi_bcast_fullmolecule(MyMolecule)    
     CALL LS_GETTIM(CPU2,WALL2)
     CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
     WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
  endif StartUpSlaves
#endif

  ! ************************************************
  ! * Restart Option                               *
  ! ************************************************
  call mem_alloc(JobsCompleted,nOccBatchesI,nOccBatchesJ)
  JobDone = .FALSE.           
  IF(master)THEN
     IF(DECinfo%DECrestart)THEN
        !CHECK IF THERE ARE ENERGY CONTRIBUTIONS AVAILABLE
        INQUIRE(FILE='FULLMP2.restart',EXIST=file_exists)
        IF(file_exists)THEN
           WRITE(DECinfo%output,*)'Restart of Full canonical molecular MP2 calculation:'
           restart_lun = -1  !initialization
           call lsopen(restart_lun,'FULLMP2.restart','OLD','UNFORMATTED')
           rewind restart_lun
           read(restart_lun) noccRestartI,noccRestartJ
           IF(noccRestartI.NE.nOccBatchesI)THEN
              print*,'noccRestartI,nOccBatchesI',noccRestartI,nOccBatchesI
              call lsquit('nOccBatchesI must be same in Restart for MP2',-1)
           ENDIF
           IF(noccRestartJ.NE.nOccBatchesJ)THEN
              print*,'noccRestartJ,nOccBatchesJ',noccRestartJ,nOccBatchesJ
              call lsquit('nOccBatchesJ must be same in Restart for MP2',-1)
           ENDIF
           read(restart_lun) JobsCompleted
           IF(COUNT(JobsCompleted).EQ.nOccBatchesI*nOccBatchesJ)THEN
              WRITE(DECinfo%output,*)'All MP2 energies is on file JobsCompleted=',JobsCompleted
              JobDone = .TRUE.
           ELSE
              WRITE(DECinfo%output,*)'Restarting from file '
              WRITE(DECinfo%output,*) COUNT(JobsCompleted),' jobs completed out of ',nOccBatchesI*nOccBatchesJ
              JobDone = .FALSE.           
           ENDIF
           read(restart_lun) mp2_energy
           WRITE(DECinfo%output,*)'MP2 Energy Read From File: ',mp2_energy
           call lsclose(restart_lun,'KEEP')
        ELSE
           JobsCompleted = .FALSE.
        ENDIF
     ELSE
        JobDone = .FALSE.           
        JobsCompleted = .FALSE.
     ENDIF
  ELSE
     JobsCompleted = .FALSE.
     JobDone = .FALSE.
  ENDIF

#ifdef VAR_MPI
  IF(DECinfo%DECrestart)THEN
     CALL LS_GETTIM(CPU1,WALL1)
     call ls_mpibcast(JobsCompleted,nOccBatchesI,nOccBatchesJ,infpar%master,comm)
     CALL LS_GETTIM(CPU2,WALL2)
     CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
     WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
     IF(COUNT(JobsCompleted).EQ.nOccBatchesI*nOccBatchesJ)JobDone = .TRUE.
  ENDIF
#endif

  IF(JobDone)THEN
     !do nothing
  ELSE
     call mem_alloc(EpsOcc,nocc)
     !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I) &
     !$OMP SHARED(nocc,MyMolecule,EpsOcc)
     do I=1,nocc
        EpsOcc(I) = MyMolecule%ppfock(I,I)
     enddo
     !$OMP END PARALLEL DO
     call mem_alloc(EpsVirt,nvirt)
     !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I) &
     !$OMP SHARED(nvirt,MyMolecule,EpsVirt)
     do A=1,nvirt
        EpsVirt(A) = MyMolecule%qqfock(A,A)
     enddo
     !$OMP END PARALLEL DO
     
     ! ************************************************
     ! * Determine batch information for Gamma batch  *
     ! * And 
     ! * Determine batch information for Alpha batch  *
     ! ************************************************
     
     IF(DECinfo%useIchor)THEN
        iAO = 2 !Gamma is the 2. Center of the 4 center two electron coulomb integral
        !Determine how many batches of AOS based on the bat%MaxAllowedDimGamma, the requested
        !size of the AO batches. iAO is the center that the batching should occur on. 
        !'R'  !Specifies that it is the Regular AO basis that should be batched 
        call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimGamma,&
             & nbatchesGamma,DECinfo%output)
        call mem_alloc(AOGammabatchinfo,nbatchesGamma)
        !Construct the batches of AOS based on the bat%MaxAllowedDimGamma, the requested
        !size of the AO batches - bat%MaxAllowedDimGamma must be unchanged since the call 
        !to determine_Ichor_nbatchesofAOS
        !MaxActualDimGamma is an output parameter indicating How big the biggest batch was, 
        !So MaxActualDimGamma must be less og equal to bat%MaxAllowedDimGamma
        call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimGamma,&
             & nbatchesGamma,AOGammabatchinfo,MaxActualDimGamma,DECinfo%output)
        
        iAO = 1 !Alpha is the 1. Center of the 4 center two electron coulomb integral
        !Determine how many batches of AOS based on the bat%MaxAllowedDimAlpha, the requested
        !size of the AO batches. iAO is the center that the batching should occur on. 
        !'R'  !Specifies that it is the Regular AO basis that should be batched 
        call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimAlpha,&
             & nbatchesAlpha,DECinfo%output)
        call mem_alloc(AOAlphabatchinfo,nbatchesAlpha)
        !Construct the batches of AOS based on the bat%MaxAllowedDimAlpha, the requested
        !size of the AO batches - bat%MaxAllowedDimAlpha must be unchanged since the call 
        !to determine_Ichor_nbatchesofAOS
        !MaxActualDimAlpha is an output parameter indicating How big the biggest batch was, 
        !So MaxActualDimAlpha must be less og equal to bat%MaxAllowedDimAlpha
        call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimAlpha,&
             & nbatchesAlpha,AOAlphabatchinfo,MaxActualDimAlpha,DECinfo%output)
     ELSE
        ! Orbital to batch information
        ! ----------------------------
        call mem_alloc(orb2batchGamma,nbasis)
        call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimGamma,&
             & nbasis,MaxActualDimGamma,batchsizeGamma,batchdimGamma,batchindexGamma,nbatchesGamma,&
             & orb2BatchGamma,'R')
        ! Translate batchindex to orbital index
        ! -------------------------------------
        call mem_alloc(batch2orbGamma,nbatchesGamma)
        do idx=1,nbatchesGamma
           call mem_alloc(batch2orbGamma(idx)%orbindex,batchdimGamma(idx) )
           batch2orbGamma(idx)%orbindex = 0
           batch2orbGamma(idx)%norbindex = 0
        end do
        do iorb=1,nbasis
           idx = orb2batchGamma(iorb)
           batch2orbGamma(idx)%norbindex = batch2orbGamma(idx)%norbindex+1
           K = batch2orbGamma(idx)%norbindex
           batch2orbGamma(idx)%orbindex(K) = iorb
        end do
        ! Orbital to batch information
        ! ----------------------------
        call mem_alloc(orb2batchAlpha,nbasis)
        call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimAlpha,&
             & nbasis,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,nbatchesAlpha,&
             & orb2BatchAlpha,'R')
        ! Translate batchindex to orbital index
        ! -------------------------------------
        call mem_alloc(batch2orbAlpha,nbatchesAlpha)
        do idx=1,nbatchesAlpha
           call mem_alloc(batch2orbAlpha(idx)%orbindex,batchdimAlpha(idx) )
           batch2orbAlpha(idx)%orbindex = 0
           batch2orbAlpha(idx)%norbindex = 0
        end do
        do iorb=1,nbasis
           idx = orb2batchAlpha(iorb)
           batch2orbAlpha(idx)%norbindex = batch2orbAlpha(idx)%norbindex+1
           K = batch2orbAlpha(idx)%norbindex
           batch2orbAlpha(idx)%orbindex(K) = iorb
        end do
     ENDIF
     
     CALL LS_GETTIM(CPU4,WALL4)
     
     ! ************************************************
     ! * Screening                                    *
     ! ************************************************
     IF(DECinfo%useIchor)THEN
        !Calculate Screening integrals 
        SameMOL = .TRUE. !Specifies same molecule on all centers 
        call SCREEN_ICHORERI_DRIVER(DECinfo%output,iprint,mylsitem%setting,INTSPEC,SameMOL)
     ELSE
        ! This subroutine builds the full screening matrix.
        call II_precalc_DECScreenMat(DECscreen,DECinfo%output,6,mylsitem%setting,&
             & nbatchesAlpha,nbatchesGamma,INTSPEC,DECinfo%IntegralThreshold)
        IF(doscreen)THEN
           call II_getBatchOrbitalScreen(DecScreen,mylsitem%setting,&
                & nb,nbatchesAlpha,nbatchesGamma,batchsizeAlpha,batchsizeGamma,&
                & batchindexAlpha,batchindexGamma,&
                & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
        ENDIF
     ENDIF
     
     CALL LS_GETTIM(CPU3,WALL3)
     CPU_AOINT = CPU_AOINT + (CPU3-CPU4)
     WALL_AOINT = WALL_AOINT + (WALL3-WALL4)

     ! ************************************************
     ! * Main Loop                                    *
     ! ************************************************
     JobInfo1Free = .FALSE.
     FullRHS = (nbatchesGamma.EQ.1).AND.(nbatchesAlpha.EQ.1)
     
     nJobs = 0 
     BatchOccI: do iB = 1,nOccbatchesI
        nOccBatchDimI = nOccBatchDimImax
        IF(MOD(nOcc,nOccBatchDimI).NE.0.AND.iB.EQ.nOccBatchesI)THEN
           !the remainder
           nOccBatchDimI = MOD(nOcc,nOccBatchDimImax)
        ENDIF
        
        !construct CoI(nb,nOccBatchDimI)
        call mem_alloc(CoI,nb,nOccBatchDimI)       
        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I) &
        !$OMP SHARED(nOccBatchDimI,MyMolecule,CoI,iB,nOccBatchDimImax)
        do I=1,nOccBatchDimI
           CoI(:,I) = Mymolecule%Co(:,I+(iB-1)*nOccBatchDimImax) 
        enddo
        !$OMP END PARALLEL DO

        startjB = 1
        IF(PermutationalSymmetryIJ) startjB = iB
        BatchOccJ: do jB = startjB,nOccbatchesJ
           nOccBatchDimJ = nOccBatchDimJmax
           IF(MOD(nOcc,nOccBatchDimJ).NE.0.AND.jB.EQ.nOccBatchesJ)THEN
              !the remainder
              nOccBatchDimJ = MOD(nOcc,nOccBatchDimJmax)
           ENDIF
           IF(JobsCompleted(iB,jB))CYCLE
           
           nJobs = nJobs + 1 
#ifdef VAR_MPI
           ! MPI: Only do this Job if this is a task for this particular rank
           if(MOD(nJobs,numnodes) .NE. mynum) cycle
#endif
           !construct CoJ(nb,nOccBatchDimJ)
           call mem_alloc(CoJ,nb,nOccBatchDimJ)       
           !TODO: OMP Workshare/Loop
           !$OMP PARALLEL DO DEFAULT(none) PRIVATE(J) &
           !$OMP SHARED(nOccBatchDimJ,MyMolecule,CoJ,jB,nOccBatchDimJmax)
           do J=1,nOccBatchDimJ
              CoJ(:,J) = Mymolecule%Co(:,J+(jB-1)*nOccBatchDimJmax) 
           enddo
           !$OMP END PARALLEL DO
           IF(.NOT.FullRHS)THEN
              call mem_alloc(VOVO,nvirt*nOccBatchDimJ,nvirt*nOccBatchDimI)
           ENDIF
           BatchGamma: do gammaB = 1,nbatchesGamma  ! AO batches
              IF(DECinfo%useIchor)THEN
                 dimGamma = AOGammabatchinfo(gammaB)%dim         ! Dimension of gamma batch
                 GammaStart = AOGammabatchinfo(gammaB)%orbstart  ! First orbital index in gamma batch
                 GammaEnd = AOGammabatchinfo(gammaB)%orbEnd      ! Last orbital index in gamma batch
                 AOGammaStart = AOGammabatchinfo(gammaB)%AOstart ! First AO batch index in gamma batch
                 AOGammaEnd = AOGammabatchinfo(gammaB)%AOEnd     ! Last AO batch index in gamma batch
              ELSE
                 dimGamma = batchdimGamma(gammaB)                           ! Dimension of gamma batch
                 GammaStart = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
                 GammaEnd = batch2orbGamma(gammaB)%orbindex(dimGamma)       ! Last index in gamma batch
              ENDIF
              IF(nbatchesAlpha.GT.1)THEN
                 call mem_alloc(VGVO,nvirt*dimGamma,nvirt*nOccBatchDimJ)
              ENDIF
              BatchAlpha: do alphaB = 1,nbatchesAlpha  ! AO batches
                 IF(DECinfo%useIchor)THEN
                    dimAlpha = AOAlphabatchinfo(alphaB)%dim         ! Dimension of alpha batch
                    AlphaStart = AOAlphabatchinfo(alphaB)%orbstart  ! First orbital index in alpha batch
                    AlphaEnd = AOAlphabatchinfo(alphaB)%orbEnd      ! Last orbital index in alpha batch
                    AOAlphaStart = AOAlphabatchinfo(alphaB)%AOstart ! First AO batch index in alpha batch
                    AOAlphaEnd = AOAlphabatchinfo(alphaB)%AOEnd     ! Last AO batch index in alpha batch
                 ELSE
                    dimAlpha = batchdimAlpha(alphaB)                                ! Dimension of alpha batch
                    AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)                 ! First index in alpha batch
                    AlphaEnd = batch2orbAlpha(alphaB)%orbindex(dimAlpha)            ! Last index in alpha batch
                 ENDIF
                 
                 
                 CALL LS_GETTIM(CPU4,WALL4)
                 IF(DECinfo%useIchor)THEN
                    call mem_alloc(tmp1,dimAlpha*dimGamma,nb*nb)                    
                    !(dimAlpha,dimGamma,nb,nb)
                    call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,Mylsitem%setting,dimAlpha,dimGamma,nb,nb,&
                         & tmp1,INTSPEC,FULLRHS,AOAlphaStart,AOAlphaEnd,AOGammaStart,AOGammaEnd,&
                         & 1,nAObatches,1,nAObatches,MoTrans,dimAlpha,dimGamma,nb,nb,NoSymmetry,DECinfo%IntegralThreshold)
                    CALL LS_GETTIM(CPU3,WALL3)
                    CPU_AOINT = CPU_AOINT + (CPU3-CPU4)
                    WALL_AOINT = WALL_AOINT + (WALL3-WALL4)
                 ELSE
                    call mem_alloc(tmp1,dimAlpha*dimGamma,nb*nb)
                    IF(doscreen) mylsitem%setting%LST_GAB_RHS => DECSCREEN%masterGabRHS
                    IF(doscreen) mylsitem%setting%LST_GAB_LHS => DECSCREEN%batchGab(alphaB,gammaB)%p
                    call II_GET_DECPACKED4CENTER_J_ERI2(DECinfo%output,DECinfo%output, &
                         & mylsitem%setting,tmp1,batchindexAlpha(alphaB),batchindexGamma(gammaB),&
                         & batchsizeAlpha(alphaB),batchsizeGamma(gammaB),dimAlpha,&
                         & dimGamma,nbasis,nbasis,FullRHS,INTSPEC,DECinfo%IntegralThreshold)
                    CALL LS_GETTIM(CPU3,WALL3)
                    CPU_AOINT = CPU_AOINT + (CPU3-CPU4)
                    WALL_AOINT = WALL_AOINT + (WALL3-WALL4)
                 ENDIF

                 !tmp2(dimAlpha,dimGamma,nb,nOccBatchDimJ) = tmp1(dimAlpha,dimGamma,nb,nb)*CoJ(nb,nOccBatchDimJ)
                 call mem_alloc(tmp2,dimAlpha*dimGamma,nb*nOccBatchDimJ)
                 M = dimAlpha*dimGamma*nb     !rows of Output Matrix
                 N = nOccBatchDimJ            !columns of Output Matrix
                 K = nb                       !summation dimension
                 call dgemm('N','N',M,N,K,1.0E0_realk,tmp1,M,CoJ,K,0.0E0_realk,tmp2,M)
                 call mem_dealloc(tmp1)
                 
                 !reorder: tmp3(nb,nOccBatchDimJ,dimAlpha,dimGamma) = tmp2(dimAlpha,dimGamma,nb,nOccBatchDimJ)
                 call mem_alloc(tmp3,nb*nOccBatchDimJ,dimAlpha*dimGamma)
                 M = dimAlpha*dimGamma   !row of Input Matrix
                 N = nb*nOccBatchDimJ    !columns of Input Matrix
                 call mat_transpose(M,N,1.0E0_realk,tmp2,0.0E0_realk,tmp3)
                 call mem_dealloc(tmp2)

                 !tmp4(nvirt,nOccBatchDimJ,dimAlpha,dimGamma) = Cv(nb,nvirt)*tmp3(nb,nOccBatchDimJ,dimAlpha,dimGamma)
                 call mem_alloc(tmp4,nvirt*nOccBatchDimJ,dimAlpha*dimGamma)
                 M = nvirt                            !rows of Output Matrix
                 N = nOccBatchDimJ*dimAlpha*dimGamma  !columns of Output Matrix
                 K = nb                               !summation dimension
                 call dgemm('T','N',M,N,K,1.0E0_realk,MyMolecule%Cv,K,tmp3,K,0.0E0_realk,tmp4,M)
                 call mem_dealloc(tmp3)
                 
                 !reorder: tmp5(dimAlpha,dimGamma,nvirt,nOccBatchDimJ) <= tmp4(nvirt,nOccBatchDimJ,dimAlpha,dimGamma)
                 call mem_alloc(tmp5,dimAlpha*dimGamma,nvirt*nOccBatchDimJ)
                 M = nvirt*nOccBatchDimJ   !row of Input Matrix
                 N = dimAlpha*dimGamma     !columns of Input Matrix
                 call mat_transpose(M,N,1.0E0_realk,tmp4,0.0E0_realk,tmp5)
                 call mem_dealloc(tmp4)
                 
                 call mem_alloc(CvA,dimAlpha,nvirt)
                 do B=1,nvirt
                    CvA(1:dimAlpha,B) = Mymolecule%Cv(AlphaStart:AlphaEnd,B)
                 enddo
                 !VGVO(nvirt,dimGamma,nvirt,nOccBatchDimJ) = CvA(dimAlpha,nvirt)*tmp5(dimAlpha,dimGamma,nvirt,nOccBatchDimJ)
                 IF(nbatchesAlpha.EQ.1)THEN
                    call mem_alloc(VGVO,nvirt*dimGamma,nvirt*nOccBatchDimJ)
                 ENDIF
                 M = nvirt                         !rows of Output Matrix
                 N = dimGamma*nvirt*nOccBatchDimJ  !columns of Output Matrix
                 K = dimAlpha                      !summation dimension
                 IF(alphaB.EQ.1)THEN
                    call dgemm('T','N',M,N,K,1.0E0_realk,CvA,K,tmp5,K,0.0E0_realk,VGVO,M) 
                 ELSE
                    call dgemm('T','N',M,N,K,1.0E0_realk,CvA,K,tmp5,K,1.0E0_realk,VGVO,M) !ADD
                 ENDIF
                 call mem_dealloc(tmp5)
                 call mem_dealloc(CvA)

                 CALL LS_GETTIM(CPU4,WALL4)
                 CPU_AOTOMO = CPU_AOTOMO + (CPU4-CPU3)
                 WALL_AOTOMO = WALL_AOTOMO + (WALL4-WALL3)
              ENDDO BatchAlpha
              CALL LS_GETTIM(CPU3,WALL3)
              
              !reorder: tmp7(nvirt,nOccBatchDimJ,nvirt,dimGamma) <= VGVO(nvirt,dimGamma,nvirt,nOccBatchDimJ)
              call mem_alloc(tmp7,nvirt*nOccBatchDimJ,nvirt*dimGamma)
              M = nvirt*dimGamma          !row of Input Matrix
              N = nvirt*nOccBatchDimJ     !columns of Input Matrix
              call mat_transpose(M,N,1.0E0_realk,VGVO,0.0E0_realk,tmp7)
              call mem_dealloc(VGVO)
              
              call mem_alloc(CoIG,dimGamma,nOccBatchDimI)
              !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I) SHARED(nOccBatchDimI,CoI,CoIG,dimGamma,GammaStart,GammaEnd)              
              do I=1,nOccBatchDimI
                 CoIG(1:dimGamma,I) = CoI(GammaStart:GammaEnd,I)
              enddo
              !$OMP END PARALLEL DO
              !VOVO(nvirt,nOccBatchDimJ,nvirt,nOccBatchDimI) = tmp7(nvirt,nOccBatchDimJ,nvirt,dimGamma)*CoIG(dimGamma,nOccBatchDimI)
              IF(FullRHS)THEN
                 call mem_alloc(VOVO,nvirt*nOccBatchDimJ,nvirt*nOccBatchDimI)
              ENDIF
              M = nvirt*nOccBatchDimJ*nvirt     !rows of Output Matrix
              N = nOccBatchDimI                 !columns of Output Matrix
              K = dimGamma                      !summation dimension
              IF(gammaB.EQ.1)THEN
                 call dgemm('N','N',M,N,K,1.0E0_realk,tmp7,M,CoIG,K,0.0E0_realk,VOVO,M) 
              ELSE
                 call dgemm('N','N',M,N,K,1.0E0_realk,tmp7,M,CoIG,K,1.0E0_realk,VOVO,M) !ADD TO VOVO
              ENDIF
              call mem_dealloc(tmp7)     
              
              call mem_dealloc(CoIG)

              CALL LS_GETTIM(CPU4,WALL4)
              CPU_AOTOMO = CPU_AOTOMO + (CPU4-CPU3)
              WALL_AOTOMO = WALL_AOTOMO + (WALL4-WALL3)
           ENDDO BatchGamma
           CALL LS_GETTIM(CPU3,WALL3)
           call mem_dealloc(CoJ)
           call mem_alloc(Amat,nvirt,nvirt)
           call mem_alloc(Bmat,nvirt,nvirt)
           IF(PermutationalSymmetryIJ)THEN
              WRITE(DECinfo%output,*)'Occ Contribution iB=',iB,', jB=',jB
              IF(iB.NE.jB)THEN 
                 tmp_mp2_energy = 0.0E0_realk
                 do I=1,nOccBatchDimI
                    do J=1,nOccBatchDimJ
                       epsIJ = EpsOcc(I+(iB-1)*nOccBatchDimImax) + EpsOcc(J+(jB-1)*nOccBatchDimJmax)
                       CALL CalcAmat2(nOccBatchDimJ,nOccBatchDimI,nvirt,VOVO,Amat,J,I)
                       call CalcBmat(nvirt,EpsIJ,EpsVirt,Amat,Bmat)
                       tmp_mp2_energy2 = 0.0E0_realk
                       call MP2_EnergyContribution(nvirt,Amat,Bmat,tmp_mp2_energy2)
                       WRITE(DECinfo%output,*)'E1(',iB,',',jB,')=',tmp_mp2_energy2
                       tmp_mp2_energy = tmp_mp2_energy + tmp_mp2_energy2
                    enddo
                 enddo
                 !all these contributions appear twice due to permutational symmetry ( I <-> J )
                 tmp_mp2_energy = 2.0E0_realk*tmp_mp2_energy
                 WRITE(DECinfo%output,*)'PermutationalSymmetryIJ E(Triangular)',tmp_mp2_energy
              ELSE !iB = jB same block 
                 tmp_mp2_energy = 0.0E0_realk
                 do I=1,nOccBatchDimI
                    do J=I+1,nOccBatchDimJ
                       epsIJ = EpsOcc(I+(iB-1)*nOccBatchDimImax) + EpsOcc(J+(jB-1)*nOccBatchDimJmax)
                       CALL CalcAmat2(nOccBatchDimJ,nOccBatchDimI,nvirt,VOVO,Amat,J,I)
                       call CalcBmat(nvirt,EpsIJ,EpsVirt,Amat,Bmat)
                       tmp_mp2_energy2 = 0.0E0_realk
                       call MP2_EnergyContribution(nvirt,Amat,Bmat,tmp_mp2_energy2)
                       WRITE(DECinfo%output,*)'E2(',iB,',',jB,')=',tmp_mp2_energy2
                       tmp_mp2_energy = tmp_mp2_energy + tmp_mp2_energy2
                    enddo
                 enddo
                 !all these contributions appear twice due to permutational symmetry ( I <-> J )
                 tmp_mp2_energy = 2.0E0_realk*tmp_mp2_energy
                 WRITE(DECinfo%output,*)'PermutationalSymmetryIJ E(Triangular)',tmp_mp2_energy
                 !all these contributions only appear once since I=J in diagonal (iB,iB) block
                 do I=1,nOccBatchDimI
                    epsIJ = 2.0E0_realk*EpsOcc(I+(iB-1)*nOccBatchDimImax)
                    CALL CalcAmat2(nOccBatchDimI,nOccBatchDimI,nvirt,VOVO,Amat,I,I)
                    call CalcBmat(nvirt,EpsIJ,EpsVirt,Amat,Bmat)
                    tmp_mp2_energy2 = 0.0E0_realk
                    call MP2_EnergyContribution(nvirt,Amat,Bmat,tmp_mp2_energy2)
                    WRITE(DECinfo%output,*)'E3(',iB,',',jB,')=',tmp_mp2_energy2
                    tmp_mp2_energy = tmp_mp2_energy + tmp_mp2_energy2
                 enddo
                 WRITE(DECinfo%output,*)'PermutationalSymmetryIJ E(Triangular+diagonal)',tmp_mp2_energy
              ENDIF
              WRITE(DECinfo%output,*)'canon MP2 energy contribution =',tmp_mp2_energy
           ELSE
              tmp_mp2_energy = 0.0E0_realk
              do I=1,nOccBatchDimI
                 do J=1,nOccBatchDimJ
                    epsIJ = EpsOcc(I+(iB-1)*nOccBatchDimImax) + EpsOcc(J+(jB-1)*nOccBatchDimJmax)                    
                    CALL CalcAmat2(nOccBatchDimJ,nOccBatchDimI,nvirt,VOVO,Amat,J,I)
                    call CalcBmat(nvirt,EpsIJ,EpsVirt,Amat,Bmat)
                    tmp_mp2_energy2 = 0.0E0_realk
                    call MP2_EnergyContribution(nvirt,Amat,Bmat,tmp_mp2_energy2)
                    WRITE(DECinfo%output,*)'E4(',iB,',',jB,')=',tmp_mp2_energy2
                    tmp_mp2_energy = tmp_mp2_energy + tmp_mp2_energy2
                 enddo
              enddo
              WRITE(DECinfo%output,*)'canon MP2 energy contribution =',tmp_mp2_energy
           ENDIF
           call mem_dealloc(Amat)
           call mem_dealloc(Bmat)
           call mem_dealloc(VOVO)
           CALL LS_GETTIM(CPU4,WALL4)
           CPU_ECONT = CPU_ECONT + (CPU4-CPU3)
           WALL_ECONT = WALL_ECONT + (WALL4-WALL3)           
#ifdef VAR_MPI
           IF(master)THEN
              mp2_energy = mp2_energy + tmp_mp2_energy
              JobsCompleted(iB,jB) = .TRUE.       
              IF(PermutationalSymmetryIJ) JobsCompleted(jB,iB) = .TRUE.       
              !test to see if job info have been recieved from slaves
              MessageRecieved = .TRUE.
              MessageRecievedW = .TRUE.
              CALL LS_GETTIM(CPU1,WALL1)
              DO WHILE(MessageRecievedW)
                 call lsmpi_iprobe(comm,MessageRecieved,mpistatus)
                 MessageRecievedW = MessageRecieved
                 IF(MessageRecievedW)THEN
                    !get the sender ID
                    senderID = mpistatus(MPI_SOURCE)
                    nMPI = 2
                    call lsmpi_recv(JobInfo1,nMPI,senderID,TAG1,comm)
                    call lsmpi_recv(tmp_mp2_energy,senderID,TAG2,comm)
                    mp2_energy = mp2_energy + tmp_mp2_energy
                    JobsCompleted(JobInfo1(1),JobInfo1(2)) = .TRUE.
                    IF(PermutationalSymmetryIJ) JobsCompleted(JobInfo1(2),JobInfo1(1)) = .TRUE.
                 ENDIF
              ENDDO
              CALL LS_GETTIM(CPU2,WALL2)
              CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
              WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
           ELSE
              masterrank = 0
              !send info to master
              IF(JobInfo1Free)THEN
                 !wait for master to recieve the first mp2_energy
                 CALL LS_GETTIM(CPU1,WALL1)
                 call lsmpi_wait(request1)
                 call lsmpi_wait(request2)
                 CALL LS_GETTIM(CPU2,WALL2)
                 CPU_MPIWAIT = CPU_MPIWAIT + (CPU2-CPU1)
                 WALL_MPIWAIT = WALL_MPIWAIT + (WALL2-WALL1)
              ENDIF
              nMPI = 2
              JobInfo1(1) = iB
              JobInfo1(2) = jB
              call lsmpi_Isend(JobInfo1,nMPI,masterrank,TAG1,comm,request1)
              mp2_energy = tmp_mp2_energy
              call lsmpi_Isend(mp2_energy,masterrank,TAG2,comm,request2)
              JobInfo1Free = .TRUE.
           ENDIF
#else
           mp2_energy = mp2_energy + tmp_mp2_energy
           JobsCompleted(iB,jB) = .TRUE.
           IF(PermutationalSymmetryIJ) JobsCompleted(jB,iB) = .TRUE.
#endif
           IF(master)THEN
              !Restart File 
              restart_lun = -1  !initialization
              call lsopen(restart_lun,'FULLMP2.restart','UNKNOWN','UNFORMATTED')
              rewind restart_lun
              write(restart_lun) nOccbatchesI,nOccbatchesJ
              write(restart_lun) JobsCompleted
              write(restart_lun) mp2_energy
              call lsclose(restart_lun,'KEEP')
           ENDIF
        enddo BatchOccJ !Batched Occupied J
        call mem_dealloc(CoI)       
     enddo BatchOccI !Batched Occupied I

#ifdef VAR_MPI
     !Wait for all slaves to be finished
     IF(master)THEN
        NotAllMessagesRecieved = COUNT(JobsCompleted).NE.nOccbatchesI*nOccbatchesJ
        DO WHILE(NotAllMessagesRecieved)
           call lsmpi_probe(comm,mpistatus)
           senderID = mpistatus(MPI_SOURCE)
           nMPI = 2
           CALL LS_GETTIM(CPU1,WALL1)
           call lsmpi_recv(JobInfo1,nMPI,senderID,TAG1,comm)
           CALL LS_GETTIM(CPU2,WALL2)
           CPU_MPIWAIT = CPU_MPIWAIT + (CPU2-CPU1)
           WALL_MPIWAIT = WALL_MPIWAIT + (WALL2-WALL1)
           call lsmpi_recv(tmp_mp2_energy,senderID,TAG2,comm)
           CALL LS_GETTIM(CPU1,WALL1)
           CPU_MPICOMM = CPU_MPICOMM + (CPU1-CPU2)
           WALL_MPICOMM = WALL_MPICOMM + (WALL1-WALL2)
           mp2_energy = mp2_energy + tmp_mp2_energy
           IF(JobInfo1(1).NE.JobInfo1(2))mp2_energy = mp2_energy + tmp_mp2_energy
           JobsCompleted(JobInfo1(1),JobInfo1(2)) = .TRUE.
           IF(PermutationalSymmetryIJ) JobsCompleted(JobInfo1(2),JobInfo1(1)) = .TRUE.
           NotAllMessagesRecieved = COUNT(JobsCompleted).NE.nOccbatchesI*nOccbatchesJ
        ENDDO
        !Restart File 
        restart_lun = -1  !initialization
        call lsopen(restart_lun,'FULLMP2.restart','UNKNOWN','UNFORMATTED')
        rewind restart_lun
        write(restart_lun) nOccbatchesI,nOccbatchesJ
        write(restart_lun) JobsCompleted
        write(restart_lun) mp2_energy
        call lsclose(restart_lun,'KEEP')
     ENDIF
#endif
     
     IF(DECinfo%useIchor)THEN
        call FREE_SCREEN_ICHORERI()
        call mem_dealloc(AOGammabatchinfo)
        call mem_dealloc(AOAlphabatchinfo)
     ELSE
        nullify(mylsitem%setting%LST_GAB_LHS)
        nullify(mylsitem%setting%LST_GAB_RHS)
        call free_decscreen(DECSCREEN)
        ! Free gamma batch stuff
        call mem_dealloc(orb2batchGamma)
        call mem_dealloc(batchdimGamma)
        call mem_dealloc(batchsizeGamma)
        call mem_dealloc(batchindexGamma)
        orb2batchGamma => null()
        batchdimGamma => null()
        batchsizeGamma => null()
        batchindexGamma => null()
        do idx=1,nbatchesGamma
           call mem_dealloc(batch2orbGamma(idx)%orbindex)
           batch2orbGamma(idx)%orbindex => null()
        end do
        
        call mem_dealloc(batch2orbGamma)
        batch2orbGamma => null()
        
        ! Free alpha batch stuff
        call mem_dealloc(orb2batchAlpha)
        call mem_dealloc(batchdimAlpha)
        call mem_dealloc(batchsizeAlpha)
        call mem_dealloc(batchindexAlpha)
        orb2batchAlpha => null()
        batchdimAlpha => null()
        batchsizeAlpha => null()
        batchindexAlpha => null()
        do idx=1,nbatchesAlpha
           call mem_dealloc(batch2orbAlpha(idx)%orbindex)
           batch2orbAlpha(idx)%orbindex => null()
        end do
        call mem_dealloc(batch2orbAlpha)
        batch2orbAlpha => null()
     ENDIF
     call mem_dealloc(EpsOcc)
     call mem_dealloc(EpsVirt)
  ENDIF
  call mem_dealloc(JobsCompleted)

  IF(MASTER)THEN
     write(lupri,*)  ''
     write(lupri,*)  'MP2 CORRELATION ENERGY = ', mp2_energy
     write(*,'(1X,a,f20.10)') 'MP2 CORRELATION ENERGY = ', mp2_energy
     write(lupri,*)  ''
  ENDIF
  CALL LSTIMER('FULL CANONICAL MP2 ',TS,TE,DECINFO%OUTPUT,FORCEPRINT)
  CALL ls_TIMTXT('>>>  WALL Time used in MP2 AO integral evaluation ',WALL_AOINT,lupri)
  CALL ls_TIMTXT('>>>  CPU Time used in MP2 AO integral evaluation  ',CPU_AOINT,lupri)
  CALL ls_TIMTXT('>>>  WALL Time used in MP2 AO to MO transformation',WALL_AOTOMO,lupri)
  CALL ls_TIMTXT('>>>  CPU Time used in MP2 AO to MO transformation ',CPU_AOTOMO,lupri)
  CALL ls_TIMTXT('>>>  WALL Time used in MP2 Energy evaluation      ',WALL_ECONT,lupri)
  CALL ls_TIMTXT('>>>  CPU Time used in MP2 Energy evaluation       ',CPU_ECONT,lupri)
#ifdef VAR_MPI
  write(lupri,*)'Overall Time spent in MPI Communication and MPI Wait for rank=',infpar%mynum
  CALL ls_TIMTXT('>>>  WALL Time used MPI Communication inc. some Wait',WALL_MPICOMM,lupri)
  CALL ls_TIMTXT('>>>  CPU Time used MPI Communication inc. some Wait',CPU_MPICOMM,lupri)
  CALL ls_TIMTXT('>>>  WALL Time used in MPI Wait',WALL_MPIWAIT,lupri)
  CALL ls_TIMTXT('>>>  CPU Time used in MPI Wait',CPU_MPIWAIT,lupri)
#endif    
  write(lupri,*) ' '
end subroutine full_canonical_mp2

subroutine get_optimal_batch_sizes_for_canonical_mp2(MinAObatch,nbasis,nocc,nvirt,&
     & MaxAllowedDimAlpha,MaxAllowedDimGamma,nOccBatchDimI,nOccBatchDimJ,numnodes,lupri)
  implicit none
  integer,intent(in) :: MinAObatch,nbasis,nocc,nvirt,lupri,numnodes
  integer,intent(inout) :: MaxAllowedDimAlpha,MaxAllowedDimGamma
  integer,intent(inout) :: nOccBatchDimI,nOccBatchDimJ
  !local variables
  real(realk) :: MemoryAvailable,GB,AG,maxsize
  integer :: iB,jB,iiB,jjB,nb,dimGamma,dimAlpha,AB,iGB,nTasks
  real(realk) :: nOccTMP,nbasisTMP
  logical :: BOTH
  ! Memory currently available
  ! **************************
  call get_currently_available_memory(MemoryAvailable)
  ! Note: We multiply by 85 % to be on the safe side!
  MemoryAvailable = 0.85E0_realk*MemoryAvailable
  GB = 1.000E-9_realk 
  nb = nbasis
  nOccTMP = nocc
  nbasisTMP = nbasis
  !test if full is possible
  if(DECinfo%manual_occbatchsizes)then
   nOccBatchDimI = DECinfo%batchOccI
   nOccBatchDimJ = DECinfo%batchOccJ
   WRITE(DECinfo%output,*)'CANONMP2MEM: Occupied batches chosen in input:',nOccBatchDimI,nOccBatchDimJ
   IF(DECinfo%manual_batchsizes)THEN
      MaxAllowedDimAlpha = DECinfo%ccsdAbatch
      MaxAllowedDimGamma = DECinfo%ccsdGbatch
      WRITE(DECinfo%output,*)'CANONMP2MEM: AO batches chosen in input:',MaxAllowedDimAlpha,MaxAllowedDimGamma
   ELSE
    iB = nOccBatchDimI
    jB = nOccBatchDimJ
    BatchAlpha3: do AB = 1,nbasis
     dimAlpha = MAX(MinAObatch,CEILING(nbasisTMP/AB)) !(nbasis/1,nbasis/2,..)
     BatchGamma3: do iGB = 1,nbasis
      dimGamma = MAX(MinAObatch,CEILING(nbasisTMP/iGB)) !(nbasis/1,nbasis/2,..)
      call memestimateCANONMP2(iB,jB,dimAlpha,dimGamma,nb,nvirt,maxsize)
      MaxAllowedDimAlpha = dimAlpha
      MaxAllowedDimGamma = dimGamma
      IF(maxsize.LT.MemoryAvailable)THEN
         WRITE(DECinfo%output,*)'CANONMP2MEM: mem estimate means a batching of Occupied Index I'
         WRITE(DECinfo%output,*)'Estimated Memory usage is ',maxsize,' GB'
         exit BatchAlpha3
      ENDIF
      IF(dimGamma.EQ.MinAObatch)exit BatchGamma3
     enddo BatchGamma3
     IF(dimAlpha.EQ.MinAObatch)exit BatchAlpha3
    enddo BatchAlpha3
    WRITE(DECinfo%output,*)'CANONMP2MEM: Final AO batches chosen:',MaxAllowedDimAlpha,MaxAllowedDimGamma
   ENDIF
  else
   WRITE(DECinfo%output,*)'Test Full N**4 possible'
   call memestimateCANONMP2(nOcc,nOcc,nb,nb,nb,nvirt,maxsize)
   IF(maxsize.LT.MemoryAvailable)THEN
    MaxAllowedDimAlpha = nb
    MaxAllowedDimGamma = nb
    nOccBatchDimI = nocc
    nOccBatchDimJ = nocc
    WRITE(DECinfo%output,*)'CANONMP2MEM: Final Occupied batches chosen:',nOccBatchDimI,nOccBatchDimJ
    WRITE(DECinfo%output,*)'CANONMP2MEM: Final AO batches chosen:',MaxAllowedDimAlpha,MaxAllowedDimGamma
    WRITE(DECinfo%output,*)'CANONMP2MEM: mem estimate means full N**4 dim are used'
    WRITE(DECinfo%output,*)'Estimated Memory is ',maxsize,' GB'
   ELSE
    WRITE(DECinfo%output,*)'Full scheme would require',maxsize,' GB'
    IF(numnodes.GT.1)THEN
     WRITE(DECinfo%output,*)'Test scheme reducing the size of the nOccBatchDimI since numnodes=',numnodes
     !reduce the size of the nOccBatchDimI
     iB = CEILING(nOccTMP/numnodes)
     !assume minimum AO batches dimGamma = MinAObatch, dimAlpha = MinAObatch
     call memestimateCANONMP2(iB,nOcc,MinAObatch,MinAObatch,nb,nvirt,maxsize)
     IF(maxsize.LT.0.5E0_realk*MemoryAvailable)THEN
      nOccBatchDimI = iB
      nOccBatchDimJ = nocc
      WRITE(DECinfo%output,*)'CANONMP2MEM: Final Occupied batches chosen:',nOccBatchDimI,nOccBatchDimJ
      IF(DECinfo%manual_batchsizes)THEN
       MaxAllowedDimAlpha = DECinfo%ccsdAbatch
       MaxAllowedDimGamma = DECinfo%ccsdGbatch
      ELSE
       BatchAlpha: do AB = 1,nbasis
        dimAlpha = MAX(MinAObatch,CEILING(nbasisTMP/AB)) !(nbasis/1,nbasis/2,..)
        BatchGamma: do iGB = 1,nbasis
         dimGamma = MAX(MinAObatch,CEILING(nbasisTMP/iGB)) !(nbasis/1,nbasis/2,..)
         call memestimateCANONMP2(iB,nOcc,dimAlpha,dimGamma,nb,nvirt,maxsize)
         MaxAllowedDimAlpha = dimAlpha
         MaxAllowedDimGamma = dimGamma
         IF(maxsize.LT.MemoryAvailable)THEN
            WRITE(DECinfo%output,*)'CANONMP2MEM: mem estimate means a batching of Occupied Index I'
            WRITE(DECinfo%output,*)'Estimated Memory usage is ',maxsize,' GB'
            exit BatchAlpha
         ENDIF
         IF(dimGamma.EQ.MinAObatch)exit BatchGamma
        enddo BatchGamma
        IF(dimAlpha.EQ.MinAObatch)exit BatchAlpha
       enddo BatchAlpha
       WRITE(DECinfo%output,*)'CANONMP2MEM: Final AO batches chosen:',MaxAllowedDimAlpha,MaxAllowedDimGamma
      ENDIF
      BOTH = .FALSE.
     ELSE
      BOTH = .TRUE.
     ENDIF
    ELSE
     BOTH = .TRUE.
    ENDIF
    IF(BOTH)THEN
     WRITE(DECinfo%output,*)'Test scheme reducing both the size of the nOccBatchDimI and J'
     !reduce the size of the nOccBatchDimI and nOccBatchDimJ
     OccLoopI: do iiB = 1,nOcc
      iB = nOcc/iiB !(nOcc/2,nOcc/3)
      jB = iB
      nTasks = iiB*iiB 
      IF(nTasks.LT.numnodes)THEN
         print*,'nTasks=',ntasks,'.LT.numnodes=',numnodes,' CYCLE '
         CYCLE OccLoopI
      ENDIF
      WRITE(DECinfo%output,*)'Test Reduction in OccI,OccJ : nOccBatchDimI',IB,'nOccBatchDimJ',JB
      call memestimateCANONMP2(iB,jB,MinAObatch,MinAObatch,nb,nvirt,maxsize)
      WRITE(DECinfo%output,*)'mem for MinAObatch,MinAObatch:',maxsize,'GB'
      IF(maxsize.LT.0.5E0_realk*MemoryAvailable)THEN
       nOccBatchDimI = iB
       nOccBatchDimJ = jB
       WRITE(DECinfo%output,*)'CANONMP2MEM: Final Occupied batches chosen:',nOccBatchDimI,nOccBatchDimJ
       IF(DECinfo%manual_batchsizes)THEN
        MaxAllowedDimAlpha = DECinfo%ccsdAbatch
        MaxAllowedDimGamma = DECinfo%ccsdGbatch
       ELSE
        BatchAlpha2: do AB = 1,nbasis
         dimAlpha = MAX(MinAObatch,CEILING(nbasisTMP/AB)) !(nbasis/1,nbasis/2,..)
         BatchGamma2: do iGB = 1,nbasis
          dimGamma = MAX(MinAObatch,CEILING(nbasisTMP/iGB)) !(nbasis/1,nbasis/2,..)
          WRITE(DECinfo%output,*)'Obtain Alpha Gamma',dimAlpha,dimGamma
          call memestimateCANONMP2(iB,jB,dimAlpha,dimGamma,nb,nvirt,maxsize)
          MaxAllowedDimAlpha = dimAlpha
          MaxAllowedDimGamma = dimGamma
          IF(maxsize.LT.MemoryAvailable)THEN
             WRITE(DECinfo%output,*)'CANONMP2MEM: mem estimate means a batching of Occupied Index I and J'
             WRITE(DECinfo%output,*)'Estimated Memory usage is ',maxsize,' GB'
             exit BatchAlpha2
          ENDIF
          IF(dimGamma.EQ.MinAObatch)exit BatchGamma2
         enddo BatchGamma2
         IF(dimAlpha.EQ.MinAObatch)exit BatchAlpha2
        enddo BatchAlpha2
        WRITE(DECinfo%output,*)'CANONMP2MEM: Final AO batches chosen:',MaxAllowedDimAlpha,MaxAllowedDimGamma
       ENDIF
       EXIT OccLoopI
      ENDIF
      IF(iiB.EQ.nOcc*nOcc)THEN
         WRITE(DECinfo%output,*)
         nOccBatchDimI = 1
         nOccBatchDimJ = 1
         MaxAllowedDimAlpha = MinAObatch
         MaxAllowedDimGamma = MinAObatch
         WRITE(DECinfo%output,*)'CANONMP2MEM: Final Occupied batches chosen:',nOccBatchDimI,nOccBatchDimJ
         WRITE(DECinfo%output,*)'CANONMP2MEM: Final AO batches chosen:',MaxAllowedDimAlpha,MaxAllowedDimGamma
      ENDIF
     enddo OccLoopI
    ENDIF
   ENDIF
  ENDIF
end subroutine get_optimal_batch_sizes_for_canonical_mp2

subroutine memestimateCANONMP2(iB,jB,dimAlpha,dimGamma,nbasis,nvirt,maxsize)
  implicit none
  integer,intent(in) :: iB,jB,dimAlpha,dimGamma,nvirt,nbasis
  real(realk),intent(inout) :: maxsize
  !
  real(realk) :: step1,step2,step3,step4,step5,step6,step7,step8,step9
  real(realk) :: step10,step11,step12,GB,nb
  nb=nbasis
  GB = 1.000E-9_realk 
  !construct Co(nb,nOccBatchDimI)
  step1 = nbasis*iB
  !CoJ(nb,nOccBatchDimJ)       
  step2 = nbasis*iB + nbasis*jB

  IF(dimGamma.EQ.nbasis.AND.dimAlpha.EQ.nbasis)THEN
     step3 = nbasis*iB + nbasis*jB !noallocation of VOVO
  ELSE
     !VOVO(nvirt*nOccBatchDimJ,nvirt*nOccBatchDimI)
     step3 = nbasis*iB + nbasis*jB + nvirt*nvirt*iB*jB
  ENDIF
  IF(dimAlpha.NE.nbasis)THEN
     !VGVO(nvirt*dimGamma,nvirt*nOccBatchDimJ)
     step4 = step3 + jB*nvirt*nvirt*dimGamma
  ELSE
     step4 = step3 !noallocation of VGVO
  ENDIF
  !tmp1(dimAlpha*dimGamma,nb*nb)
  step5 = step4 + dimAlpha*dimGamma*nb*nb
  !tmp2(dimAlpha,dimGamma,nb,nOccBatchDimJ)
  step6 = step4 + dimAlpha*dimGamma*nb*nb + dimAlpha*dimGamma*nb*jB
  !call mem_alloc(tmp3,nb*nOccBatchDimJ,dimAlpha*dimGamma)
  step7 = step4 + dimAlpha*dimGamma*nb*jB + nb*jB*dimAlpha*dimGamma
  !call mem_alloc(tmp4,nvirt*nOccBatchDimJ,dimAlpha*dimGamma)
  step8 = step4 + nb*jB*dimAlpha*dimGamma + nvirt*jB*dimAlpha*dimGamma
  !call mem_alloc(tmp5,dimAlpha*dimGamma,nvirt*nOccBatchDimJ)
  step9 = step4 + nvirt*jB*dimAlpha*dimGamma + dimAlpha*dimGamma*nvirt*jB
  IF(dimAlpha.EQ.nbasis)THEN
     !call mem_alloc(CvA,dimAlpha,nvirt)+VGVO(nvirt*dimGamma,nvirt*nOccBatchDimJ)
     step10 = step4 + dimAlpha*dimGamma*nvirt*jB + dimAlpha*nvirt + jB*nvirt*nvirt*dimGamma
  ELSE
     !call mem_alloc(CvA,dimAlpha,nvirt)
     step10 = step4 + dimAlpha*dimGamma*nvirt*jB + dimAlpha*nvirt
  ENDIF
  !call mem_alloc(tmp7,nvirt*nOccBatchDimJ,nvirt*dimGamma)
  step11 = step3 + nvirt*jB*nvirt*dimGamma 
  IF(dimGamma.EQ.nbasis.AND.dimAlpha.NE.nbasis)THEN
     !call mem_alloc(CoIG,dimGamma,nOccBatchDimI)
     step12 = step3 + nvirt*jB*nvirt*dimGamma + dimGamma*iB
  ELSE
     !alloc VOVO(nvirt*nOccBatchDimJ,nvirt*nOccBatchDimI) + CoIG
     step12 = step3 + nvirt*jB*nvirt*dimGamma + dimGamma*iB + nvirt*nvirt*iB*jB
  ENDIF
  maxsize = MAX(step1,step2,step3,step4,step5,step6,step7,step8,step9,step10,step11,step12)*realk*GB
!  WRITE(DECinfo%output,*)'MemRequired Step  1:',Step1*realk*GB,' GB'
!  WRITE(DECinfo%output,*)'MemRequired Step  2:',Step2*realk*GB,' GB'
!  WRITE(DECinfo%output,*)'MemRequired Step  3:',Step3*realk*GB,' GB'
!  WRITE(DECinfo%output,*)'MemRequired Step  4:',Step4*realk*GB,' GB'
!  WRITE(DECinfo%output,*)'MemRequired Step  5:',Step5*realk*GB,' GB'
!  WRITE(DECinfo%output,*)'MemRequired Step  6:',Step6*realk*GB,' GB'
!  WRITE(DECinfo%output,*)'MemRequired Step  7:',Step7*realk*GB,' GB'
!  WRITE(DECinfo%output,*)'MemRequired Step  8:',Step8*realk*GB,' GB'
!  WRITE(DECinfo%output,*)'MemRequired Step  9:',Step9*realk*GB,' GB'
!  WRITE(DECinfo%output,*)'MemRequired Step 10:',Step10*realk*GB,' GB'
!  WRITE(DECinfo%output,*)'MemRequired Step 11:',Step11*realk*GB,' GB'
!  WRITE(DECinfo%output,*)'MemRequired Step 12:',Step12*realk*GB,' GB'
!  WRITE(DECinfo%output,*)'DECinfo%memory',DECinfo%memory
end subroutine memestimateCANONMP2

!> \brief Calculate canonical MP2 energy for full molecular system
!> \author Thomas Kjaergaard
!> \date October 2014
subroutine full_canonical_mp2B(MyMolecule,MyLsitem,mp2_energy)
  implicit none
  !> Full molecule info
  type(fullmolecule), intent(inout) :: MyMolecule
  !> Lsitem structure
  type(lsitem), intent(inout) :: mylsitem
  !> Canonical MP2 correlation energy
  real(realk),intent(inout) :: mp2_energy    
  !
  integer :: nbasis,nocc,nvirt,naux,noccfull,numnodes
  logical :: master,wakeslaves
  real(realk),pointer :: EpsOcc(:),EpsVirt(:)
  integer :: J,I,node,natoms,A,lupri,restart_lun
  logical :: MessageRecieved,MessageRecievedW,FORCEPRINT,file_exists
  real(realk) :: TS,TE,TS2,TE2,TS3,TE3,CPU1,CPU2,WALL1,WALL2,CPU_MPICOMM,WALL_MPICOMM
  real(realk) :: CPU_MPIWAIT,WALL_MPIWAIT,MemInGBCollected,epsIJ,tmp_mp2_energy
  real(realk) :: CPU3,CPU4,WALL3,WALL4,CPU_AOINT,WALL_AOINT,tmp_mp2_energy2
  real(realk) :: CPU_AOTOMO,WALL_AOTOMO,CPU_ECONT,WALL_ECONT,numnodesR
  real(realk),pointer :: Amat(:,:),Bmat(:,:),tmp1(:,:),tmp2(:,:,:,:),tmp3(:,:),tmp4(:,:)
  real(realk),pointer :: tmp5(:,:),tmp6(:,:),tmp7(:,:),CoBatchA(:,:),CoBatchB(:,:)
  real(realk),pointer :: CoBI(:,:),CoBJ(:,:),tmp62(:,:),tmp72(:,:),CoI(:,:),CoJ(:,:)
  real(realk),pointer :: CAV(:,:),CgammaMPI(:,:),tmp1b(:,:,:,:),CgammaMPI2(:,:)
  real(realk),pointer :: VOVO(:,:),VGVO(:,:),CvA(:,:),CoIG(:,:),VOVO2(:,:)
!
  type(DecAObatchinfo),pointer :: Alphabatchinfo(:),Gammabatchinfo(:)
  type(batchtoorb), pointer :: batch2orbAlpha(:),batch2orbGamma(:)
  integer, pointer :: orb2batchGamma(:), batchsizeGamma(:), batchindexGamma(:),DimAOGamma(:),JobGamma(:)
  integer, pointer :: orb2batchAlpha(:), batchsizeAlpha(:), batchindexAlpha(:),DimAOAlpha(:),JobAlpha(:)
  integer, pointer :: batchdimGamma(:),dimGammaArray(:),GammaIndexArray(:),offsetGamma(:,:)
  integer, pointer :: AOstartGamma(:,:),AOendGamma(:,:),AOdimGamma(:,:),AOdimGammaMPI(:)
  integer, pointer :: batchdimAlpha(:),dimAlphaArray(:),AlphaIndexArray(:),offsetAlpha(:,:)
  integer, pointer :: AOstartAlpha(:,:),AOendAlpha(:,:),AOdimAlpha(:,:),AOdimAlphaMPI(:)
  integer, pointer :: nBlocksGamma(:),nBlocksAlpha(:),nOccBatchDimJrank(:),OccIndexJrank(:,:)
  integer :: inodeLoop,jnodeLoop,nodeLoop,nrownodes,ncolnodes,nBlocks,MAXnBlocksGamma,MAXnBlocksAlpha
  integer :: MaxAllowedDimAlpha,MaxAllowedDimGamma,MaxAllowedDimAlphaMPI,MaxAllowedDimGammaMPI
  integer :: MaxActualDimAlpha,MaxActualDimGamma,dimAOoffsetG,dimAOoffsetA
  integer :: nOccBatchDimImax,K,iorb,nbatchesGamma,nbatchesAlpha,nb,inode,jnode,ibatch,offsetG,offsetA
  integer :: iAO,nAObatches,AOGammaStart,AOGammaEnd,AOAlphaStart,AOAlphaEnd,iprint,M,N,iB
  integer :: MinAObatch,gammaB,alphaB,nOccBatchesI,jB,nOccBatchDimI,nOccBatchDimJ
  integer :: nbatchesGammaRestart,nbatchesAlphaRestart,dimGamma,GammaStart,GammaEnd,dimAlpha
  integer :: AlphaStart,AlphaEnd,B,noccRestartI,noccRestartJ,nJobs,startjB,idxx
  integer :: sqrtnumnodes,gB,idx(1),dimAOoffset,kk,nBlocksG,nBlocksA
  integer :: dimGammaMPI,dimAlphaMPI,ibatchG,ibatchA,dimAlpha2,dimGamma2
  integer :: IMYNUMNBATCHES1,IMYNUMNBATCHES2,nOccBatchDimJmax
  integer :: nOccbatchesIrestart,noccIstart,nbuf1,nbuf2,Ibuf(8)
  logical :: MoTrans, NoSymmetry,SameMol,JobDone,JobInfo1Free,FullRHS,doscreen,NotAllMessagesRecieved
  logical :: PermutationalSymmetryIJ,SetdimGamma,dim1Includezero,dim2Includezero
  logical,pointer :: JobsCompleted(:,:)
  integer(kind=8) :: maxsize
  TYPE(DECscreenITEM)   :: DecScreen
  Character            :: intSpec(5)
  integer(kind=ls_mpik)  :: nMPI,TAG,IERR,request,Receiver,sender,comm,TAG1,TAG2
  integer(kind=ls_mpik)  :: request1,request2,masterrank,senderID,mynum,lsmpinode
!  integer(kind=ls_mpik)  :: Sender
#ifdef VAR_MPI
  integer(kind=ls_mpik)  :: mpistatus(MPI_STATUS_SIZE)
#endif
  integer(kind=4) :: JobInfo1(2)
!  Character(80)        :: FilenameCS,FilenamePS

#ifdef VAR_TIME    
    FORCEPRINT = .TRUE.
#endif

  CALL LSTIMER('START ',TS,TE,DECINFO%OUTPUT)
  mp2_energy = 0.0E0_realk
  
  CPU_AOINT = 0.0E0_realk
  WALL_AOINT = 0.0E0_realk
  CPU_AOTOMO = 0.0E0_realk
  WALL_AOTOMO = 0.0E0_realk
  CPU_ECONT = 0.0E0_realk
  WALL_ECONT = 0.0E0_realk

  CPU_MPICOMM = 0.0E0_realk
  WALL_MPICOMM = 0.0E0_realk
  CPU_MPIWAIT = 0.0E0_realk
  WALL_MPIWAIT = 0.0E0_realk
  TAG = 1411; TAG1 = 1412; TAG2 = 1413
  
  !sanity check
  if(.NOT.DECinfo%use_canonical) then
     call lsquit('Error: full_canonical_mp2B require canonical Orbitals',-1)
  endif
  ! Init stuff
  ! **********
  nbasis = MyMolecule%nbasis
  nb = nbasis
  nocc   = MyMolecule%nocc
  nvirt  = MyMolecule%nunocc
  noccfull = nocc
  nAtoms = MyMolecule%nAtoms
  LUPRI = DECinfo%output

#ifdef VAR_MPI
  comm = MPI_COMM_LSDALTON
  master= (infpar%mynum == infpar%master)
  mynum = infpar%mynum
  numnodes = infpar%nodtot
  wakeslaves = infpar%nodtot.GT.1
  IF(.NOT.master)LUPRI = 6 !standard Output
#else
  ! If MPI is not used, consider the single node to be "master"
  master=.true.
  mynum = 0
  numnodes = 1
  wakeslaves = .false.
#endif

#ifdef VAR_MPI 
  ! Master starts up slave
  StartUpSlaves: if(wakeslaves .and. master) then
     ! Wake up slaves to do the job: slaves awoken up with (RIMP2FULL)
     ! and call full_canonical_rimp2_slave which communicate info 
     ! then calls full_canonical_rimp2.
     CALL LS_GETTIM(CPU1,WALL1)
     call ls_mpibcast(CANONMP2FULL,infpar%master,comm)
     ! Communicate fragment information to slaves
     call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,comm)
     call mpicopy_lsitem(MyLsitem,comm)
     call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,comm)
     call mpi_bcast_fullmolecule(MyMolecule)    
     CALL LS_GETTIM(CPU2,WALL2)
     CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
     WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
  endif StartUpSlaves
#endif

  ! Set integral info
  ! *****************
  INTSPEC(1)='R' !R = Regular Basis set on the 1th center 
  INTSPEC(2)='R' !R = Regular Basis set on the 2th center 
  INTSPEC(3)='R' !R = Regular Basis set on the 3th center 
  INTSPEC(4)='R' !R = Regular Basis set on the 4th center 
  INTSPEC(5)='C' !C = Coulomb operator

  doscreen = mylsitem%setting%scheme%cs_screen.OR.&
       & mylsitem%setting%scheme%ps_screen

  !determine MinAObatch: the minimum allowed AObatch size + number of AO batches
  IF(DECinfo%useIchor)THEN
     !Determine the minimum allowed AObatch size MinAObatch
     !In case of pure Helium atoms in cc-pVDZ ((4s,1p) -> [2s,1p]) MinAObatch = 3 (Px,Py,Pz)
     !In case of pure Carbon atoms in cc-pVDZ ((9s,4p,1d) -> [3s,2p,1d]) MinAObatch=6 (the 2*(Px,Py,Pz))
     !In case of pure Carbon atoms in 6-31G   ((10s,4p) -> [3s,2p]) MinAObatch = 3 (Px,Py,Pz) 
     !'R'  !Specifies that it is the Regular AO basis that should be batched
     iAO = 4 !the center that the batching should occur on (they are all the same in this case)  
     call determine_MinimumAllowedAObatchSize(MyLsItem%setting,iAO,'R',MinAObatch)
     iprint = 0           !print level for Ichor Integral code
     MoTrans = .FALSE.    !Do not transform to MO basis! 
     NoSymmetry = .FALSE. !Use Permutational Symmetry! 
     SameMol = .TRUE.     !Same molecule on all centers of the 4 center 2 electron integral
     !Determine the full number of AO batches - not to be confused with the batches of AOs
     !Required by the MAIN_ICHORERI_DRIVER unless all four dimensions are batched 
     iAO = 1
     call determine_Ichor_nAObatches(mylsitem%setting,iAO,'R',nAObatches,DECinfo%output)
     nullify(Gammabatchinfo)
     nullify(Alphabatchinfo)    
  ELSE
     ! The smallest possible AO batch depends on the basis set
     ! (More precisely, if all batches are made as small as possible, then the
     !  call below determines the largest of these small batches).
     call determine_maxBatchOrbitalsize(DECinfo%output,mylsitem%setting,MinAObatch,'R')
     nullify(orb2batchAlpha)
     nullify(batchdimAlpha)
     nullify(batchsizeAlpha)
     nullify(batch2orbAlpha)
     nullify(batchindexAlpha)
     nullify(orb2batchGamma)
     nullify(batchdimGamma)
     nullify(batchsizeGamma)
     nullify(batch2orbGamma)
     nullify(batchindexGamma)
  ENDIF

  
  ! ***************************************************************************************
  ! Get optimal values of: MaxAllowedDimAlpha,MaxAllowedDimGamma,nOccBatchDimImax,nOccBatchDimJmax
  ! ***************************************************************************************
  call get_optimal_batch_sizes_for_canonical_mp2B(MinAObatch,nbasis,nocc,nvirt,&
       & numnodes,nrownodes,ncolnodes,MaxAllowedDimAlpha,MaxAllowedDimGamma,&
       & MaxAllowedDimAlphaMPI,MaxAllowedDimGammaMPI,nOccBatchDimImax,nOccBatchDimJmax)

#ifdef VAR_MPI
  !use the numbers obtained by master  
  IF(master)THEN
     Ibuf(1) = nrownodes
     Ibuf(2) = ncolnodes
     Ibuf(3) = MaxAllowedDimAlpha
     Ibuf(4) = MaxAllowedDimGamma
     Ibuf(5) = MaxAllowedDimAlphaMPI
     Ibuf(6) = MaxAllowedDimGammaMPI
     Ibuf(7) = nOccBatchDimImax
     Ibuf(8) = nOccBatchDimJmax
  ENDIF
  nbuf1 = 8 
  call ls_mpibcast(Ibuf,nbuf1,infpar%master,comm)
  IF(.NOT.master)THEN
     nrownodes = Ibuf(1) 
     ncolnodes = Ibuf(2)
     MaxAllowedDimAlpha = Ibuf(3)
     MaxAllowedDimGamma = Ibuf(4)
     MaxAllowedDimAlphaMPI = Ibuf(5)
     MaxAllowedDimGammaMPI = Ibuf(6)
     nOccBatchDimImax = Ibuf(7)
     nOccBatchDimJmax = Ibuf(8)
  ENDIF
#endif

  IF(master)THEN
     write(DECinfo%output,*)'nbasis               ',nbasis
     write(DECinfo%output,*)'nocc                 ',nocc
     write(DECinfo%output,*)'nvirt                ',nvirt
     write(DECinfo%output,*)'MinAObatch           ',MinAObatch
     write(DECinfo%output,*)'numnodes             ',numnodes
     write(DECinfo%output,*)'nrownodes            ',nrownodes
     write(DECinfo%output,*)'ncolnodes            ',ncolnodes
     write(DECinfo%output,*)'MaxAllowedDimAlpha   ',MaxAllowedDimAlpha
     write(DECinfo%output,*)'MaxAllowedDimGamma   ',MaxAllowedDimGamma
     write(DECinfo%output,*)'MaxAllowedDimAlphaMPI',MaxAllowedDimAlphaMPI
     write(DECinfo%output,*)'MaxAllowedDimGammaMPI',MaxAllowedDimGammaMPI
     write(DECinfo%output,*)'nOccBatchDimImax     ',nOccBatchDimImax
     write(DECinfo%output,*)'nOccBatchDimJmax     ',nOccBatchDimJmax
  ENDIF

  ! ************************************************
  ! * Determine batch information for Gamma batch  *
  ! * And 
  ! * Determine batch information for Alpha batch  *
  ! ************************************************
  
  IF(DECinfo%useIchor)THEN
     iAO = 2 !Gamma is the 2. Center of the 4 center two electron coulomb integral
     !Determine how many batches of AOS based on the bat%MaxAllowedDimGamma, the requested
     !size of the AO batches. iAO is the center that the batching should occur on. 
     !'R'  !Specifies that it is the Regular AO basis that should be batched 
     call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimGamma,&
          & nbatchesGamma,DECinfo%output)
     call mem_alloc(Gammabatchinfo,nbatchesGamma)
     !Construct the batches of AOS based on the bat%MaxAllowedDimGamma, the requested
     !size of the AO batches - bat%MaxAllowedDimGamma must be unchanged since the call 
     !to determine_Ichor_nbatchesofAOS
     !MaxActualDimGamma is an output parameter indicating How big the biggest batch was, 
     !So MaxActualDimGamma must be less og equal to bat%MaxAllowedDimGamma
     call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimGamma,&
          & nbatchesGamma,Gammabatchinfo,MaxActualDimGamma,DECinfo%output)
     
     iAO = 1 !Alpha is the 1. Center of the 4 center two electron coulomb integral
     !Determine how many batches of AOS based on the bat%MaxAllowedDimAlpha, the requested
     !size of the AO batches. iAO is the center that the batching should occur on. 
     !'R'  !Specifies that it is the Regular AO basis that should be batched 
     call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimAlpha,&
          & nbatchesAlpha,DECinfo%output)
     call mem_alloc(Alphabatchinfo,nbatchesAlpha)
     !Construct the batches of AOS based on the bat%MaxAllowedDimAlpha, the requested
     !size of the AO batches - bat%MaxAllowedDimAlpha must be unchanged since the call 
     !to determine_Ichor_nbatchesofAOS
     !MaxActualDimAlpha is an output parameter indicating How big the biggest batch was, 
     !So MaxActualDimAlpha must be less og equal to bat%MaxAllowedDimAlpha
     call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimAlpha,&
          & nbatchesAlpha,Alphabatchinfo,MaxActualDimAlpha,DECinfo%output)
  ELSE
     ! Orbital to batch information
     ! ----------------------------
     call mem_alloc(orb2batchGamma,nbasis)
     call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimGamma,&
          & nbasis,MaxActualDimGamma,batchsizeGamma,batchdimGamma,batchindexGamma,nbatchesGamma,&
          & orb2BatchGamma,'R')
     ! Translate batchindex to orbital index
     ! -------------------------------------
     call mem_alloc(batch2orbGamma,nbatchesGamma)
     do idxx=1,nbatchesGamma
        call mem_alloc(batch2orbGamma(idxx)%orbindex,batchdimGamma(idxx) )
        batch2orbGamma(idxx)%orbindex = 0
        batch2orbGamma(idxx)%norbindex = 0
     end do
     do iorb=1,nbasis
        idxx = orb2batchGamma(iorb)
        batch2orbGamma(idxx)%norbindex = batch2orbGamma(idxx)%norbindex+1
        K = batch2orbGamma(idxx)%norbindex
        batch2orbGamma(idxx)%orbindex(K) = iorb
     end do
     ! Orbital to batch information
     ! ----------------------------
     call mem_alloc(orb2batchAlpha,nbasis)
     call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimAlpha,&
          & nbasis,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,nbatchesAlpha,&
          & orb2BatchAlpha,'R')
     ! Translate batchindex to orbital index
     ! -------------------------------------
     call mem_alloc(batch2orbAlpha,nbatchesAlpha)
     do idxx=1,nbatchesAlpha
        call mem_alloc(batch2orbAlpha(idxx)%orbindex,batchdimAlpha(idxx) )
        batch2orbAlpha(idxx)%orbindex = 0
        batch2orbAlpha(idxx)%norbindex = 0
     end do
     do iorb=1,nbasis
        idxx = orb2batchAlpha(iorb)
        batch2orbAlpha(idxx)%norbindex = batch2orbAlpha(idxx)%norbindex+1
        K = batch2orbAlpha(idxx)%norbindex
        batch2orbAlpha(idxx)%orbindex(K) = iorb
     end do
  ENDIF

  IF(master)THEN
     write(DECinfo%output,*)'MaxActualDimAlpha    ',MaxActualDimAlpha
     write(DECinfo%output,*)'MaxActualDimGamma    ',MaxActualDimGamma
  ENDIF

  ! ************************************************
  ! * Screening                                    *
  ! ************************************************
  CALL LS_GETTIM(CPU3,WALL3)
  IF(DECinfo%useIchor)THEN
     !Calculate Screening integrals 
     SameMOL = .TRUE. !Specifies same molecule on all centers 
     call SCREEN_ICHORERI_DRIVER(DECinfo%output,iprint,mylsitem%setting,INTSPEC,SameMOL)
  ELSE
     ! This subroutine builds the full screening matrix.
     call II_precalc_DECScreenMat(DECscreen,DECinfo%output,6,mylsitem%setting,&
          & nbatchesAlpha,nbatchesGamma,INTSPEC,DECinfo%IntegralThreshold)
     IF(doscreen)THEN
        call II_getBatchOrbitalScreen(DecScreen,mylsitem%setting,&
             & nb,nbatchesAlpha,nbatchesGamma,batchsizeAlpha,batchsizeGamma,&
             & batchindexAlpha,batchindexGamma,&
             & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
     ENDIF
  ENDIF
  CALL LS_GETTIM(CPU4,WALL4)
  CPU_AOINT = CPU_AOINT + (CPU4-CPU3)
  WALL_AOINT = WALL_AOINT + (WALL4-WALL3)

  !GAMMA SET (ncolnodes)
  call mem_alloc(dimGammaArray,nbatchesGamma)
  call mem_alloc(GammaIndexArray,nbatchesGamma)
  !determine which node have which gammaB batch 
  do gammaB = 1,nbatchesGamma  ! AO batches
     IF(DECinfo%useIchor)THEN
        dimGamma = Gammabatchinfo(gammaB)%dim         ! Dimension of gamma batch
     ELSE
        dimGamma = batchdimGamma(gammaB)           ! Dimension of gamma batch
     ENDIF
     dimGammaArray(gammaB) = dimGamma
     GammaIndexArray(gammaB) = gammaB
  enddo
  call integer_inv_sort_with_tracking(dimGammaArray,GammaIndexArray,nbatchesGamma)
  call mem_alloc(DimAOGamma,ncolnodes)
  call mem_alloc(JobGamma,nbatchesGamma)
  call mem_alloc(nBlocksGamma,ncolnodes)
  DimAOGamma = 0
  nBlocksGamma = 0
  do gB = 1,nbatchesGamma  ! AO batches
     gammaB = GammaIndexArray(gB)
     dimGamma = dimGammaArray(gB)
     idx = MINLOC(DimAOGamma)
     DimAOGamma(idx(1)) = DimAOGamma(idx(1)) + DimGamma
     JobGamma(gammaB) = idx(1)
     nBlocksGamma(idx(1)) = nBlocksGamma(idx(1)) + 1
  enddo
  call mem_dealloc(DimAOGamma)
  call mem_dealloc(dimGammaArray)
  call mem_dealloc(GammaIndexArray)

!  print*,'nBlocksGamma',nBlocksGamma,'MYNUM',MYNUM
  MAXnBlocksGamma = MAXVAL(nBlocksGamma)

  call mem_alloc(offsetGamma,MAXnBlocksGamma,ncolnodes)
  offsetGamma=0
  call mem_alloc(AOstartGamma,MAXnBlocksGamma,ncolnodes)
  AOstartGamma=0
  call mem_alloc(AOendGamma,MAXnBlocksGamma,ncolnodes)
  AOendGamma=0
  call mem_alloc(AOdimGamma,MAXnBlocksGamma,ncolnodes)
  AOdimGamma=0
  call mem_alloc(AOdimGammaMPI,ncolnodes)
  AOdimGammaMPI=0
  
  do jnode = 1,ncolnodes
     dimAOoffsetG = 0 
     nBlocks = 0
     do gammaB = 1,nbatchesGamma
        IF(JobGamma(gammaB).EQ.jnode)THEN 
           nBlocks = nBlocks + 1
           offsetGamma(nBlocks,jnode) = dimAOoffsetG
           IF(DECinfo%useIchor)THEN
              dimGamma = Gammabatchinfo(gammaB)%dim         ! Dimension of gamma batch
              GammaStart = Gammabatchinfo(gammaB)%orbstart  ! First orbital index in gamma batch
              GammaEnd = Gammabatchinfo(gammaB)%orbEnd      ! Last orbital index in gamma batch
           ELSE
              dimGamma = batchdimGamma(gammaB)                ! Dimension of gamma batch
              GammaStart = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
              GammaEnd = batch2orbGamma(gammaB)%orbindex(dimGamma)       ! Last index in gamma batch
           ENDIF
           dimAOoffsetG = dimAOoffsetG + dimGamma
           AOstartGamma(nBlocks,jnode) = GammaStart
           AOendGamma(nBlocks,jnode) = GammaEnd
           AOdimGamma(nBlocks,jnode) = dimGamma
           AOdimGammaMPI(jnode) =  dimAOoffsetG
        endif
     enddo
  enddo
!  do jnode = 1,ncolnodes
!     print*,'AOstartGamma(:,jnode=',jnode,')',AOstartGamma(1:nBlocks,jnode)
!     print*,'AOendGamma(:,jnode=',jnode,')',AOendGamma(1:nBlocks,jnode)
!     print*,'AOdimGamma(:,jnode=',jnode,')',AOdimGamma(1:nBlocks,jnode)
!     print*,'AOdimGammaMPI(jnode=',jnode,')',AOdimGammaMPI(jnode)
!  enddo

  call mem_alloc(dimAlphaArray,nbatchesAlpha)
  call mem_alloc(AlphaIndexArray,nbatchesAlpha)
  !ALPHA SET (nrownodes)
  !determine which node have which alphaB batch 
  do alphaB = 1,nbatchesAlpha  ! AO batches
     IF(DECinfo%useIchor)THEN
        dimAlpha = Alphabatchinfo(alphaB)%dim         ! Dimension of alpha batch
     ELSE
        dimAlpha = batchdimAlpha(alphaB)           ! Dimension of alpha batch
     ENDIF
     dimAlphaArray(alphaB) = dimAlpha
     AlphaIndexArray(alphaB) = alphaB
  enddo
  call integer_inv_sort_with_tracking(dimAlphaArray,AlphaIndexArray,nbatchesAlpha)
  call mem_alloc(DimAOAlpha,nrownodes)
  call mem_alloc(JobAlpha,nbatchesAlpha)
  call mem_alloc(nBlocksAlpha,nrownodes)
  DimAOAlpha = 0
  nBlocksALpha = 0
  do gB = 1,nbatchesAlpha  ! AO batches
     alphaB = AlphaIndexArray(gB)
     dimAlpha = dimAlphaArray(gB)
     idx = MINLOC(DimAOAlpha)
     DimAOAlpha(idx(1)) = DimAOAlpha(idx(1)) + DimAlpha
     JobAlpha(alphaB) = idx(1)
     nBlocksAlpha(idx(1)) = nBlocksAlpha(idx(1)) + 1
  enddo
  call mem_dealloc(DimAOAlpha)
  call mem_dealloc(dimAlphaArray)
  call mem_dealloc(AlphaIndexArray)

!  print*,'nBlocksAlpha',nBlocksAlpha,'MYNUM',MYNUM
  MAXnBlocksAlpha = MAXVAL(nBlocksAlpha)

  call mem_alloc(offsetAlpha,MAXnBlocksAlpha,nrownodes)
  offsetAlpha=0
  call mem_alloc(AOstartAlpha,MAXnBlocksAlpha,nrownodes)
  AOstartAlpha=0
  call mem_alloc(AOendAlpha,MAXnBlocksAlpha,nrownodes)
  AOendAlpha=0
  call mem_alloc(AOdimAlpha,MAXnBlocksAlpha,nrownodes)
  AOdimAlpha=0
  call mem_alloc(AOdimAlphaMPI,nrownodes)
  AOdimAlphaMPI=0
  
  do inode = 1,nrownodes
     dimAOoffsetA = 0 
     nBlocks = 0
     do alphaB = 1,nbatchesAlpha
        IF(JobAlpha(alphaB).EQ.inode)THEN 
           nBlocks = nBlocks + 1
           offsetAlpha(nBlocks,inode) = dimAOoffsetA
           IF(DECinfo%useIchor)THEN
              dimAlpha = Alphabatchinfo(alphaB)%dim         ! Dimension of alpha batch
              AlphaStart = Alphabatchinfo(alphaB)%orbstart  ! First orbital index in alpha batch
              AlphaEnd = Alphabatchinfo(alphaB)%orbEnd      ! Last orbital index in alpha batch
           ELSE
              dimAlpha = batchdimAlpha(alphaB)                ! Dimension of alpha batch
              AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)            ! First index in alpha batch
              AlphaEnd = batch2orbAlpha(alphaB)%orbindex(dimAlpha)       ! Last index in alpha batch
           ENDIF
           dimAOoffsetA = dimAOoffsetA + dimAlpha
           AOstartAlpha(nBlocks,inode) = AlphaStart
           AOendAlpha(nBlocks,inode) = AlphaEnd
           AOdimAlpha(nBlocks,inode) = dimAlpha
           AOdimAlphaMPI(inode) =  dimAOoffsetA
        endif
     enddo
  enddo
!  do inode = 1,nrownodes
!     print*,'AOstartAlpha(:,inode=',inode,')',AOstartAlpha(1:nBlocks,inode)
!     print*,'AOendAlpha(:,inode=',inode,')',AOendAlpha(1:nBlocks,inode)
!     print*,'AOdimAlpha(:,inode=',inode,')',AOdimAlpha(1:nBlocks,inode)
!     print*,'AOdimAlphaMPI(inode=',inode,')',AOdimAlphaMPI(inode)
!  enddo

  inode = mod(mynum,nrownodes)+1   !0 => 1, 1 => 2, 2 => 1, 3 => 2
  jnode = (mynum)/(nrownodes) + 1  !0 => 1, 1 => 1, 2 => 2, 3 => 2
  IF(inode+(jnode-1)*nrownodes.NE.mynum+1)call lsquit('node mismatch in MP2',-1)
     
  dimAlphaMPI = AOdimAlphaMPI(inode)
  dimGammaMPI = AOdimGammaMPI(jnode)

  !building
  !nOccBatchDimJrank(mynum)
  !OccIndexJrank(J,mynum)              
  IF(numnodes.GT.nocc)THEN
     call lsquit('Error more nodes then occupied orbitals: FIXME',-1)
  ENDIF

  dim1Includezero = .TRUE.
  J = numnodes-1
  call mem_alloc(nOccBatchDimJrank,J,dim1Includezero) !include zero nOccBatchDimJrank(0:numnodes-1)
  nOccBatchDimJrank = 0
  do J=1,nocc
     idx = MINLOC(nOccBatchDimJrank(0:numnodes-1))     
     nOccBatchDimJrank(idx(1)-1) = nOccBatchDimJrank(idx(1)-1) + 1
  enddo
  J = numnodes-1
  IF(MAXVAL(nOccBatchDimJrank).GT.nOccBatchDimJmax.OR.MAXVAL(nOccBatchDimJrank).LT.0)THEN
     print*,'nOccBatchDimJrank',nOccBatchDimJrank
     print*,'nOccBatchDimJmax',nOccBatchDimJmax
     print*,'MAXVAL(nOccBatchDimJrank)',MAXVAL(nOccBatchDimJrank)
     call lsquit('miscalc nOccBatchDimJrank full canon mp2b ',-1)
  ENDIF
  dim1Includezero = .FALSE.
  dim2Includezero = .TRUE.
  I = MAXVAL(nOccBatchDimJrank)
  J = numnodes-1
  call mem_alloc(OccIndexJrank,I,J,dim1Includezero,dim2Includezero)   
  nOccBatchDimJrank = 0
  do J=1,nocc
     idx = MINLOC(nOccBatchDimJrank(0:numnodes-1))
     nOccBatchDimJrank(idx(1)-1) = nOccBatchDimJrank(idx(1)-1) + 1
     OccIndexJrank(nOccBatchDimJrank(idx(1)-1),idx(1)-1) = J
  enddo

  nOccBatchesI = nOcc/nOccBatchDimImax
  IF(MOD(nOcc,nOccBatchDimImax).NE.0)nOccBatchesI = nOccBatchesI + 1  
  PermutationalSymmetryIJ = .FALSE.
  write(DECinfo%output,*)'nOccBatchesI         ',nOccBatchesI
  !IF(nrownodes.NE.nOccBatchesI)call lsquit('error in nocbatches in mp2',-1)

  call mem_alloc(EpsOcc,nocc)
  !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I) &
  !$OMP SHARED(nocc,MyMolecule,EpsOcc)
  do I=1,nocc
     EpsOcc(I) = MyMolecule%ppfock(I,I)
  enddo
  !$OMP END PARALLEL DO
  call mem_alloc(EpsVirt,nvirt)
  !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I) &
  !$OMP SHARED(nvirt,MyMolecule,EpsVirt)
  do A=1,nvirt
     EpsVirt(A) = MyMolecule%qqfock(A,A)
  enddo
  !$OMP END PARALLEL DO

  ! ***************************************************************************************
  ! Restart
  ! ***************************************************************************************
  IF(DECinfo%DECrestart)THEN
     !CHECK IF THERE ARE ENERGY CONTRIBUTIONS AVAILABLE
     INQUIRE(FILE='FULLMP2B.restart',EXIST=file_exists)
     IF(file_exists)THEN
        IF(master)THEN
           WRITE(DECinfo%output,*)'Restart of Full molecular MP2(MP2B) calculation:'
        ENDIF
        restart_lun = -1  !initialization
        call lsopen(restart_lun,'FULLMP2B.restart','OLD','FORMATTED')
        rewind restart_lun
        read(restart_lun,'(I9)') nOccbatchesIrestart
        read(restart_lun,'(I9)') noccIstart
        IF(nOccbatchesIrestart.NE.nOccbatchesI)THEN
           print*,'Restart Error: nOccbatchesIrestart=',nOccbatchesIrestart
           print*,'Restart Error: nOccbatchesI       =',nOccbatchesI
           call lsquit('MP2 restart error first integer is wrong')
        ELSE
           IF(noccIstart.EQ.nOccbatchesI)THEN
              IF(master)WRITE(DECinfo%output,*)'All energies is on file'
              noccIstart = nocc+1
              read(restart_lun,'(F28.16)') mp2_energy
           ELSEIF(noccIstart.GT.nOccbatchesI.OR.noccIstart.LT.1)THEN
              IF(master)THEN
                 WRITE(DECinfo%output,*)'MP2 restart error, second integer is wrong. Read:',noccIstart
              ENDIF
              call lsquit('MP2 restart error second integer is wrong')             
           ELSE
              noccIstart = noccIstart + 1
              read(restart_lun,'(F28.16)') mp2_energy
           ENDIF
        ENDIF
        call lsclose(restart_lun,'KEEP')
     ELSE
        noccIstart=1
        mp2_energy = 0.0E0_realk
     ENDIF
  ELSE
     noccIstart=1
     mp2_energy = 0.0E0_realk
  ENDIF
  
  ! ************************************************
  ! * Main Loop                                    *
  ! ************************************************
  JobInfo1Free = .FALSE.
  FullRHS = (nbatchesGamma.EQ.1).AND.(nbatchesAlpha.EQ.1)
  
  BatchOccI: do iB = noccIstart,nOccbatchesI
     nOccBatchDimI = nOccBatchDimImax
     IF(MOD(nOcc,nOccBatchDimI).NE.0.AND.iB.EQ.nOccBatchesI)THEN
        !the remainder
        nOccBatchDimI = MOD(nOcc,nOccBatchDimImax)
     ENDIF

     !construct CoI(nb,nOccBatchDimI)
     call mem_alloc(CoI,nb,nOccBatchDimI)       
     !TODO: OMP Workshare/Loop
     !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I) &
     !$OMP SHARED(nOccBatchDimI,MyMolecule,CoI,iB,nOccBatchDimImax)
     do I=1,nOccBatchDimI
        CoI(:,I) = Mymolecule%Co(:,I+(iB-1)*nOccBatchDimImax) 
     enddo
     !$OMP END PARALLEL DO
     call mem_alloc(tmp2,dimAlphaMPI,dimGammaMPI,nb,nOccBatchDimI)
     !dgemm:   tmp2(dimAlpha,dimGamma,nbasis,noccB)        
     nBlocksG = 0 
     BatchGamma: do gammaB = 1,nbatchesGamma
        nBlocksA = 0 
        BatchAlpha: do alphaB = 1,nbatchesAlpha  ! AO batches
           IF(JobAlpha(alphaB).EQ.inode.AND.JobGamma(gammaB).EQ.jnode)THEN 
              IF(nBlocksA.EQ.0) nBlocksG = nBlocksG + 1
              nBlocksA = nBlocksA + 1

              offsetA = offsetAlpha(nBlocksA,inode)
              offsetG = offsetGamma(nBlocksG,jnode)
              IF(DECinfo%useIchor)THEN
                 dimGamma = Gammabatchinfo(gammaB)%dim         ! Dimension of gamma batch
                 GammaStart = Gammabatchinfo(gammaB)%orbstart  ! First orbital index in gamma batch
                 GammaEnd = Gammabatchinfo(gammaB)%orbEnd      ! Last orbital index in gamma batch
                 AOGammaStart = Gammabatchinfo(gammaB)%AOstart ! First AO batch index in gamma batch
                 AOGammaEnd = Gammabatchinfo(gammaB)%AOEnd     ! Last AO batch index in gamma batch
                 dimAlpha = Alphabatchinfo(alphaB)%dim         ! Dimension of alpha batch
                 AlphaStart = Alphabatchinfo(alphaB)%orbstart  ! First orbital index in alpha batch
                 AlphaEnd = Alphabatchinfo(alphaB)%orbEnd      ! Last orbital index in alpha batch
                 AOAlphaStart = Alphabatchinfo(alphaB)%AOstart ! First AO batch index in alpha batch
                 AOAlphaEnd = Alphabatchinfo(alphaB)%AOEnd     ! Last AO batch index in alpha batch
              ELSE
                 dimGamma = batchdimGamma(gammaB)                      ! Dimension of gamma batch
                 GammaStart = batch2orbGamma(gammaB)%orbindex(1)       ! First index in gamma batch
                 GammaEnd = batch2orbGamma(gammaB)%orbindex(dimGamma)  ! Last index in gamma batch
                 dimAlpha = batchdimAlpha(alphaB)                      ! Dimension of alpha batch
                 AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)       ! First index in alpha batch
                 AlphaEnd = batch2orbAlpha(alphaB)%orbindex(dimAlpha)  ! Last index in alpha batch
              ENDIF

              CALL LS_GETTIM(CPU4,WALL4)
              call mem_alloc(tmp1,dimAlpha*dimGamma,nb*nb)
              IF(DECinfo%useIchor)THEN
                 !(dimAlpha,dimGamma,nb,nb)
                 call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,Mylsitem%setting,dimAlpha,dimGamma,nb,nb,&
                      & tmp1,INTSPEC,FULLRHS,AOAlphaStart,AOAlphaEnd,AOGammaStart,AOGammaEnd,&
                      & 1,nAObatches,1,nAObatches,MoTrans,dimAlpha,dimGamma,nb,nb,NoSymmetry,DECinfo%IntegralThreshold)
                 CALL LS_GETTIM(CPU3,WALL3)
                 CPU_AOINT = CPU_AOINT + (CPU3-CPU4)
                 WALL_AOINT = WALL_AOINT + (WALL3-WALL4)
              ELSE
                 IF(doscreen) mylsitem%setting%LST_GAB_RHS => DECSCREEN%masterGabRHS
                 IF(doscreen) mylsitem%setting%LST_GAB_LHS => DECSCREEN%batchGab(alphaB,gammaB)%p
                 call II_GET_DECPACKED4CENTER_J_ERI2(DECinfo%output,DECinfo%output, &
                      & mylsitem%setting,tmp1,batchindexAlpha(alphaB),batchindexGamma(gammaB),&
                      & batchsizeAlpha(alphaB),batchsizeGamma(gammaB),dimAlpha,&
                      & dimGamma,nbasis,nbasis,FullRHS,INTSPEC,DECinfo%IntegralThreshold)
              ENDIF
              CALL LS_GETTIM(CPU3,WALL3)
              CPU_AOINT = CPU_AOINT + (CPU3-CPU4)
              WALL_AOINT = WALL_AOINT + (WALL3-WALL4)
              
              !tmp1b(dimAlpha,dimGamma,nb,nOccBatchDimI)=tmp1(dimAlpha,dimGamma,nb,nb)*CoI(nb,nOccBatchDimI)
              call mem_alloc(tmp1b,dimAlpha,dimGamma,nb,nOccBatchDimI)
              M = dimAlpha*dimGamma*nb     !rows of Output Matrix
              N = nOccBatchDimI            !columns of Output Matrix
              K = nb                       !summation dimension
              call dgemm('N','N',M,N,K,1.0E0_realk,tmp1,M,CoI,K,0.0E0_realk,tmp1b,M)
              call mem_dealloc(tmp1)
              do j=1,nOccBatchDimI
                 do b=1,nb
                    do i=1,dimGamma
                       do a=1,dimAlpha
                          tmp2(a+offsetA,i+offsetG,b,j) = tmp1b(a,i,b,j)
                       enddo
                    enddo
                 enddo
              enddo
              call mem_dealloc(tmp1b)
              CALL LS_GETTIM(CPU4,WALL4)
              CPU_AOTOMO = CPU_AOTOMO + (CPU4-CPU3)
              WALL_AOTOMO = WALL_AOTOMO + (WALL4-WALL3)
           ENDIF
        enddo BatchAlpha
     enddo BatchGamma
     CALL LS_GETTIM(CPU3,WALL3)
     call mem_dealloc(CoI)

     !reorder: tmp3(nb,nOccBatchDimJ,dimAlphaMPI,dimGammaMPI) = tmp2(dimAlphaMPI,dimGammaMPI,nb,nOccBatchDimI)
     call mem_alloc(tmp3,nb*nOccBatchDimI,dimAlphaMPI*dimGammaMPI)
     M = dimAlphaMPI*dimGammaMPI   !row of Input Matrix
     N = nb*nOccBatchDimI    !columns of Input Matrix
     call mat_transpose(M,N,1.0E0_realk,tmp2,0.0E0_realk,tmp3)
     call mem_dealloc(tmp2)

     !tmp4(nvirt,nOccBatchDimI,dimAlphaMPI,dimGammaMPI) = Cv(nb,nvirt)*tmp3(nb,nOccBatchDimI,dimAlphaMPI,dimGammaMPI)
     call mem_alloc(tmp4,nvirt*nOccBatchDimI,dimAlphaMPI*dimGammaMPI)
     M = nvirt                                 !rows of Output Matrix
     N = nOccBatchDimI*dimAlphaMPI*dimGammaMPI !columns of Output Matrix
     K = nb                                    !summation dimension
     call dgemm('T','N',M,N,K,1.0E0_realk,MyMolecule%Cv,K,tmp3,K,0.0E0_realk,tmp4,M)
     call mem_dealloc(tmp3)
     
     nOccBatchDimJ = nOccBatchDimJrank(mynum)
      
     call mem_alloc(tmp5,nvirt*nOccBatchDimI,dimAlphaMPI*nOccBatchDimJ)
     tmp5 = 0.0E0_realk
     CALL LS_GETTIM(CPU4,WALL4)
     CPU_AOTOMO = CPU_AOTOMO + (CPU4-CPU3)
     WALL_AOTOMO = WALL_AOTOMO + (WALL4-WALL3)

     do jnodeLoop = 1,ncolnodes
      do inodeLoop = 1,nrownodes !alpha 
       nodeLoop = inodeLoop+(jnodeLoop-1)*nrownodes
       IF(nodeLoop.EQ.mynum+1)THEN
#ifdef VAR_MPI
        CALL LS_GETTIM(CPU1,WALL1)
        nbuf1 = nvirt*nOccBatchDimI
        nbuf2 = dimAlphaMPI*dimGammaMPI
        call ls_mpibcast(tmp4,nbuf1,nbuf2,mynum,comm)
        CALL LS_GETTIM(CPU2,WALL2)
        CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
        WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
#endif
        !construct CgammaMPI(dimGammaMPI,nOccBatchDimI)
        IF(jnodeLoop.NE.jnode)call lsquit('jnodeLoop error',-1)
        CALL LS_GETTIM(CPU3,WALL3)
        call mem_alloc(CgammaMPI,dimGammaMPI,nOccBatchDimJ)
        do J=1,nOccBatchDimJ
         JB = OccIndexJrank(J,mynum)
         do kk = 1,nBlocksGamma(jnodeLoop)
            CgammaMPI(offsetGamma(kk,jnodeLoop)+1:offsetGamma(kk,jnodeLoop)+AOdimGamma(kk,jnodeLoop),J)=&
                 &Mymolecule%Co(AOstartGamma(kk,jnodeLoop):AOendGamma(kk,jnodeLoop),JB)
         enddo
        enddo        

        !tmp5(nvirt,noccBI,dimAlpha,noccBJ)=tmp4(nvirt,noccI,dimAlpha,dimGammaMPI)*C(dimGammaMPI,noccBJ)
        M = nvirt*nOccBatchDimI*dimAlphaMPI !rows of Output Matrix
        N = nOccBatchDimJ  !columns of Output Matrix
        K = dimGammaMPI    !summation dimension
        call dgemm('N','N',M,N,K,1.0E0_realk,tmp4,M,CgammaMPI,K,1.0E0_realk,tmp5,M)
        call mem_dealloc(CgammaMPI)
        CALL LS_GETTIM(CPU4,WALL4)
        CPU_AOTOMO = CPU_AOTOMO + (CPU4-CPU3)
        WALL_AOTOMO = WALL_AOTOMO + (WALL4-WALL3)
       ELSE
        dimAlpha2 = AOdimAlphaMPI(inodeLoop)
        dimGamma2 = AOdimGammaMPI(jnodeLoop)
        call mem_alloc(tmp6,nvirt*nOccBatchDimI,dimAlpha2*dimGamma2)              
        !recv
#ifdef VAR_MPI
        CALL LS_GETTIM(CPU1,WALL1)
        nbuf1 = nvirt*nOccBatchDimI
        nbuf2 = dimAlpha2*dimGamma2
        lsmpinode = nodeLoop-1
        call ls_mpibcast(tmp6,nbuf1,nbuf2,lsmpinode,comm)
        CALL LS_GETTIM(CPU2,WALL2)
        CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
        WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)
#endif
        CALL LS_GETTIM(CPU3,WALL3)
        IF(inodeLoop.EQ.inode)THEN
         !same Alpha Batch - receiving a Gamma Batch and contract to OccJ
         IF(dimAlpha2.NE.dimAlphaMPI)call lsquit('dimAlpha2.NE.dimAlphaMPI',-1)
         call mem_alloc(CgammaMPI2,dimGamma2,nOccBatchDimJ)         
         do J=1,nOccBatchDimJ
          JB = OccIndexJrank(J,mynum)
          do kk = 1,nBlocksGamma(jnodeLoop)
           CgammaMPI2(offsetGamma(kk,jnodeLoop)+1:offsetGamma(kk,jnodeLoop)+AOdimGamma(kk,jnodeLoop),J) = &
                & Mymolecule%Co(AOstartGamma(kk,jnodeLoop):AOendGamma(kk,jnodeLoop),JB)
          enddo
         enddo
         !share alpha batches
         !tmp5(nvirt,noccBI,dimAlpha,noccBJ)=tmp6(nvirt,noccBI,dimAlpha,dimGammaMPI)*C(dimGammaMPI,noccBJ)
         M = nvirt*nOccBatchDimI*dimAlphaMPI !rows of Output Matrix
         N = nOccBatchDimJ                   !columns of Output Matrix
         K = dimGamma2                       !summation dimension
         call dgemm('N','N',M,N,K,1.0E0_realk,tmp6,M,CgammaMPI2,K,1.0E0_realk,tmp5,M)
         call mem_dealloc(CgammaMPI2)
        ELSE
           !recv something in BCAST that was not needed - could make new comm for subset.
        ENDIF
        call mem_dealloc(tmp6)
        CALL LS_GETTIM(CPU4,WALL4)
        CPU_AOTOMO = CPU_AOTOMO + (CPU4-CPU3)
        WALL_AOTOMO = WALL_AOTOMO + (WALL4-WALL3)
       ENDIF
      enddo 
     enddo
     CALL LS_GETTIM(CPU3,WALL3)
     call mem_dealloc(tmp4)

     !reorder: tmp7(dimAlpha,noccBJ,nvirt,noccBI) <= tmp5(nvirt,noccBI,dimAlpha,noccBJ)
     call mem_alloc(tmp7,dimAlphaMPI*nOccBatchDimJ,nvirt*nOccBatchDimI)
     M = nvirt*nOccBatchDimI          !row of Input Matrix
     N = dimAlphaMPI*nOccBatchDimJ    !columns of Input Matrix
     call mat_transpose(M,N,1.0E0_realk,tmp5,0.0E0_realk,tmp7)
     call mem_dealloc(tmp5)

     call mem_alloc(VOVO,nvirt*nOccBatchDimJ,nvirt*nOccBatchDimI)
     !Construct CAI(dimAlphaMPI,nvirt) = 
     call mem_alloc(CAV,dimAlphaMPI,nvirt)     
     do a=1,nvirt
        do kk = 1,nBlocksAlpha(inode)
           CAV(offsetAlpha(kk,inode)+1:offsetAlpha(kk,inode)+AOdimAlpha(kk,inode),A) = &
                Mymolecule%Cv(AOstartAlpha(kk,inode):AOendAlpha(kk,inode),A) 
        enddo
     enddo

     !VOVO(nvirt,noccBJ,nvirt,noccBI) = CAV(dimAlpha,nvirt)*tmp7(dimAlpha*nOccBatchDimJ,nvirt*nOccBatchDimI)
     M = nvirt                              !rows of Output Matrix
     N = nOccBatchDimJ*nvirt*nOccBatchDimI  !columns of Output Matrix
     K = dimAlphaMPI                        !summation dimension
     call dgemm('T','N',M,N,K,1.0E0_realk,CAV,K,tmp7,K,0.0E0_realk,VOVO,M)     
     call mem_dealloc(CAV)
     call mem_dealloc(tmp7)
     CALL LS_GETTIM(CPU4,WALL4)
     CPU_AOTOMO = CPU_AOTOMO + (CPU4-CPU3)
     WALL_AOTOMO = WALL_AOTOMO + (WALL4-WALL3)

     CALL LS_GETTIM(CPU1,WALL1)
     !at this point all have a VOVO(nvirt,noccBJ,nvirt,noccBI)
     !inode,jnode have the full VOVO(nvirt,noccBJ,nvirt,noccBI) 
     !with noccBJ determined by jnode but it only contains 
     !contributions from Alpha AO batch(inode)
     !so we now collect the Alpha AO batch(inode) to have a full VOVO(nvirt,noccBJ,nvirt,noccBI) 
     !on node inode=1 
     nbuf1 = nvirt*nOccBatchDimJ  
     nbuf2 = nvirt*nOccBatchDimI  
     do jnodeLoop = 1,ncolnodes !for all jnodes, all noccBJ contributions
      IF(jnodeLoop.EQ.jnode)THEN 
         !I share jnode and noccBJ with the rest (1:nrownodes,jnodeLoop)
         receiver = (1+(jnodeLoop-1)*nrownodes)-1
         IF(inode.EQ.1)THEN !I collect results
            do inodeLoop = 2,nrownodes !all send their contribution to inode=1,jnode=jnodeLoop
               sender = inodeLoop+(jnodeLoop-1)*nrownodes-1
               call mem_alloc(VOVO2,nvirt*nOccBatchDimJ,nvirt*nOccBatchDimI)
#ifdef VAR_MPI
               call ls_mpisendrecv(VOVO2,nbuf1,nbuf2,comm,sender,receiver)               
#endif
               do J=1,nOccBatchDimI*nvirt
                  do I=1,nOccBatchDimJ*nvirt
                     VOVO(I,J) = VOVO(I,J) + VOVO2(I,J)
                  enddo
               enddo
               call mem_dealloc(VOVO2)
            enddo
         ELSE
#ifdef VAR_MPI
            call ls_mpisendrecv(VOVO,nbuf1,nbuf2,comm,mynum,receiver)
#endif
         ENDIF
      ENDIF
     enddo
     CALL LS_GETTIM(CPU2,WALL2)
     CPU_MPICOMM = CPU_MPICOMM + (CPU2-CPU1)
     WALL_MPICOMM = WALL_MPICOMM + (WALL2-WALL1)

     CALL LS_GETTIM(CPU3,WALL3)
     IF(inode.EQ.1)THEN
        call mem_alloc(Amat,nvirt,nvirt)
        call mem_alloc(Bmat,nvirt,nvirt)
        tmp_mp2_energy = 0.0E0_realk

        !VOVO(nvirt,nOccBatchDimJ,nvirt,nOccBatchDimI) (partial cont from Alpha Batch)
        do I=1,nOccBatchDimI
           do J=1,nOccBatchDimJ
              JB = OccIndexJrank(J,mynum)              
              epsIJ = EpsOcc(I+(iB-1)*nOccBatchDimImax) + EpsOcc(JB)
              CALL CalcAmat2(nOccBatchDimJ,nOccBatchDimI,nvirt,VOVO,Amat,J,I)
              call CalcBmat(nvirt,EpsIJ,EpsVirt,Amat,Bmat)
              tmp_mp2_energy2 = 0.0E0_realk
              call MP2_EnergyContribution(nvirt,Amat,Bmat,tmp_mp2_energy2)
              tmp_mp2_energy = tmp_mp2_energy + tmp_mp2_energy2
           enddo
        enddo
!        call lsmpi_reduction(tmp_mp2_energy,infpar%master,comm)
        receiver = 0 !master ?        
        IF(master)THEN 
           do jnodeLoop = 2,ncolnodes !for all jnodes, all noccBJ contributions
              tmp_mp2_energy2 = 0.0E0_realk
              sender = 1+(jnodeLoop-1)*nrownodes-1
#ifdef VAR_MPI
              call ls_mpisendrecv(tmp_mp2_energy2,comm,sender,receiver)
#endif
              tmp_mp2_energy = tmp_mp2_energy + tmp_mp2_energy2
           enddo
        ELSE
           call ls_mpisendrecv(tmp_mp2_energy,comm,mynum,receiver)
        ENDIF
        call mem_dealloc(Amat)
        call mem_dealloc(Bmat)
     ENDIF
     call mem_dealloc(VOVO)
     CALL LS_GETTIM(CPU4,WALL4)
     CPU_ECONT = CPU_ECONT + (CPU4-CPU3)
     WALL_ECONT = WALL_ECONT + (WALL4-WALL3)           
     IF(master)THEN 
        WRITE(DECinfo%output,*)'MP2 Energy(iB=',iB,') = ',tmp_mp2_energy
        mp2_energy = mp2_energy + tmp_mp2_energy        
        !Write Restart File
        restart_lun = -1  !initialization
        call lsopen(restart_lun,'FULLMP2B.restart','UNKNOWN','FORMATTED')
        rewind restart_lun
        write(restart_lun,'(I9)') nOccbatchesI
        write(restart_lun,'(I9)') iB
        write(restart_lun,'(F28.16)') mp2_energy
        call lsclose(restart_lun,'KEEP')
     ENDIF
  enddo BatchOccI
  IF(master)THEN 
     print*,'FINAL MP2 ENERGY',mp2_energy,'MYNUM',MYNUM
  ENDIF

  call mem_dealloc(nOccBatchDimJrank)
  call mem_dealloc(OccIndexJrank)
  call mem_dealloc(JobAlpha)
  call mem_dealloc(JobGamma)
  call mem_dealloc(nBlocksAlpha)
  call mem_dealloc(nBlocksGamma)

  call mem_dealloc(offsetGamma)
  call mem_dealloc(AOstartGamma)
  call mem_dealloc(AOendGamma)
  call mem_dealloc(AOdimGamma)
  call mem_dealloc(AOdimGammaMPI)

  call mem_dealloc(offsetAlpha)
  call mem_dealloc(AOstartAlpha)
  call mem_dealloc(AOendAlpha)
  call mem_dealloc(AOdimAlpha)
  call mem_dealloc(AOdimAlphaMPI)

  IF(DECinfo%useIchor)THEN
     call FREE_SCREEN_ICHORERI()
     call mem_dealloc(Alphabatchinfo)
     call mem_dealloc(Gammabatchinfo)
  ELSE
     nullify(mylsitem%setting%LST_GAB_LHS)
     nullify(mylsitem%setting%LST_GAB_RHS)
     call free_decscreen(DECSCREEN)
     ! Free gamma batch stuff
     call mem_dealloc(orb2batchAlpha)
     call mem_dealloc(batchdimAlpha)
     call mem_dealloc(batchsizeAlpha)
     call mem_dealloc(batchindexAlpha)
     orb2batchAlpha => null()
     batchdimAlpha => null()
     batchsizeAlpha => null()
     batchindexAlpha => null()
     do idxx=1,nbatchesAlpha
        call mem_dealloc(batch2orbAlpha(idxx)%orbindex)
        batch2orbAlpha(idxx)%orbindex => null()
     end do     
     call mem_dealloc(batch2orbAlpha)
     batch2orbAlpha => null()

     call mem_dealloc(orb2batchGamma)
     call mem_dealloc(batchdimGamma)
     call mem_dealloc(batchsizeGamma)
     call mem_dealloc(batchindexGamma)
     orb2batchGamma => null()
     batchdimGamma => null()
     batchsizeGamma => null()
     batchindexGamma => null()
     do idxx=1,nbatchesGamma
        call mem_dealloc(batch2orbGamma(idxx)%orbindex)
        batch2orbGamma(idxx)%orbindex => null()
     end do
     call mem_dealloc(batch2orbGamma)
     batch2orbGamma => null()
  ENDIF
  call mem_dealloc(EpsOcc)
  call mem_dealloc(EpsVirt)
!  call mem_dealloc(JobsCompleted)

  IF(MASTER)THEN
     write(lupri,*)  ''
     write(lupri,*)  'MP2 CORRELATION ENERGY = ', mp2_energy
     write(*,'(1X,a,f20.10)') 'MP2 CORRELATION ENERGY = ', mp2_energy
     write(lupri,*)  ''
  ENDIF
  CALL LSTIMER('FULL CANONICAL MP2 ',TS,TE,DECINFO%OUTPUT,FORCEPRINT)
  CALL ls_TIMTXT('>>>  WALL Time used in MP2 AO integral evaluation ',WALL_AOINT,lupri)
  CALL ls_TIMTXT('>>>  CPU Time used in MP2 AO integral evaluation  ',CPU_AOINT,lupri)
  CALL ls_TIMTXT('>>>  WALL Time used in MP2 AO to MO transformation',WALL_AOTOMO,lupri)
  CALL ls_TIMTXT('>>>  CPU Time used in MP2 AO to MO transformation ',CPU_AOTOMO,lupri)
  CALL ls_TIMTXT('>>>  WALL Time used in MP2 Energy evaluation      ',WALL_ECONT,lupri)
  CALL ls_TIMTXT('>>>  CPU Time used in MP2 Energy evaluation       ',CPU_ECONT,lupri)
#ifdef VAR_MPI
  write(lupri,*)'Overall Time spent in MPI Communication and MPI Wait for rank=',infpar%mynum
  CALL ls_TIMTXT('>>>  WALL Time used MPI Communication inc. some Wait',WALL_MPICOMM,lupri)
  CALL ls_TIMTXT('>>>  CPU Time used MPI Communication inc. some Wait',CPU_MPICOMM,lupri)
  CALL ls_TIMTXT('>>>  WALL Time used in MPI Wait',WALL_MPIWAIT,lupri)
  CALL ls_TIMTXT('>>>  CPU Time used in MPI Wait',CPU_MPIWAIT,lupri)
#endif    
  write(lupri,*) ' '
end subroutine full_canonical_mp2B

subroutine get_optimal_batch_sizes_for_canonical_mp2B(MinAObatch,nbasis,nocc,nvirt,&
     & numnodes,nrownodes,ncolnodes,MaxAllowedDimAlpha,MaxAllowedDimGamma,&
     & MaxAllowedDimAlphaMPI,MaxAllowedDimGammaMPI,nOccBatchDimImax,nOccBatchDimJmax)
  implicit none
  integer,intent(in) :: MinAObatch,nbasis,nocc,nvirt,numnodes
  integer,intent(inout) :: MaxAllowedDimAlpha,MaxAllowedDimGamma
  integer,intent(inout) :: MaxAllowedDimAlphaMPI,MaxAllowedDimGammaMPI
  integer,intent(inout) :: nOccBatchDimImax,nOccBatchDimJmax
  integer,intent(inout) :: nrownodes,ncolnodes
  !local variables
  real(realk) :: MemoryAvailable,GB,AG,maxsize
  integer :: iB,jB,iiB,jjB,nb,dimGamma,dimAlpha,AB,iGB,nTasks,K,tmprow,tmpcol,I,J
  real(realk) :: nbasisR,noccR,nvirtR,numnodesR
  logical :: Success
  integer(kind=ls_mpik) :: nodtot 
  nbasisR = nbasis
  noccR = nocc
  nvirtR = nvirt
  numnodesR = numnodes   

  !The full_canonical_mp2B requires that (nb,nb,nb) can be distributed across all nodes 
  nodtot = INT(numnodes)
  call canonical_mp2B_memreq_test(nbasis,nodtot,Success)

  call get_currently_available_memory(MemoryAvailable)
  ! Note: We multiply by 85 % to be on the safe side!
  MemoryAvailable = 0.85E0_realk*MemoryAvailable
  GB = 1.000E-9_realk 

  !assume you have 
  !numnodes = 144
  !nbasis = 3772
  !nvirt = 3508
  !nocc = 264
  !1. canonical_mp2B_memreq_test already done
  !   nbasis*nbasis*nbasis = 430 GB can be distributed among 144 nodes (3 GB on each)  

  !2.Choose  nOccBatchDimImax as big as possible (nb*dimAlphaMPI*dimGammaMPI*nOccBatchDimImax) need to fit in mem!
  !          Same as (nb*nb*nb*nOccBatchDimImax/numnodes) 
  nOccBatchDimImax = MIN(nocc,FLOOR((MemoryAvailable*numnodes)/(nbasisR*nbasisR*nbasisR*GB))) 
  !nOccBatchDimImax = 10 for this example 
  !This means recalculation of integrals 26 times for this example 
  
  !We divide the numnodes into an array of nodes (inode,jnode)
  !for numnodes=4 
  !mynum=0    means inode=1, jnode=1
  !mynum=1    means inode=2, jnode=1
  !mynum=2    means inode=1, jnode=2
  !mynum=3    means inode=2, jnode=2
  !Chose nrownodes and ncolnodes so that: 
  !numnodes = nrownodes*ncolnodes
  !and as square as possible (but VOVO should still fit in mem)
  !ncolnodes .GE. nrownodes
  ncolnodes = numnodes
  nrownodes = 1
  K=1
  do 
     K=K+1
     IF(numnodes+1.LE.K*K)EXIT
     tmprow = K
     tmpcol = numnodes/K
     IF(tmprow*tmpcol.EQ.numnodes)THEN
        IF(tmprow+tmpcol.LE.ncolnodes+ncolnodes)THEN
           MaxAllowedDimAlphaMPI = CEILING(1.0E0_realk*nbasis/nrownodes) 
           MaxAllowedDimGammaMPI = CEILING(1.0E0_realk*nbasis/ncolnodes) 
           nOccBatchDimJmax = CEILING(1.0E0_realk*nocc/ncolnodes) 
!           print*,'CANONMP2B MinAObatch',MinAObatch
!           print*,'CANONMP2B MinAObatch',MinAObatch
!           print*,'CANONMP2B MaxAllowedDimAlphaMPI',MaxAllowedDimAlphaMPI
!           print*,'CANONMP2B MaxAllowedDimGammaMPI',MaxAllowedDimGammaMPI
!           print*,'CANONMP2B nbasis',nbasis
!           print*,'CANONMP2B nvirt',nvirt
!           print*,'CANONMP2B nOccBatchDimImax',nOccBatchDimImax
!           print*,'CANONMP2B nOccBatchDimJmax',nOccBatchDimJmax
           call memestimateCANONMP2B(MinAObatch,MinAObatch,MaxAllowedDimAlphaMPI,MaxAllowedDimGammaMPI,&
                & nbasis,nvirt,nOccBatchDimImax,nOccBatchDimJmax,maxsize)
!           print*,'TRY tmprow,tmpcol',tmprow,tmpcol,'maxsize',maxsize,'MemoryAvailable',MemoryAvailable
           IF(maxsize.LT.MemoryAvailable)THEN
              nrownodes = tmprow
              ncolnodes = tmpcol
           ENDIF
        ENDIF
     ENDIF
  enddo
  !nOccBatchesI must match number of Alpha batches(1dim) = nrownodes
  MaxAllowedDimAlphaMPI = CEILING(1.0E0_realk*nbasis/nrownodes) !19/2 = 10
  MaxAllowedDimGammaMPI = CEILING(1.0E0_realk*nbasis/ncolnodes) !19/2 = 10
  nOccBatchDimJmax = CEILING(1.0E0_realk*nocc/ncolnodes) 
!  print*,'MinAObatch',MinAObatch
  call memestimateCANONMP2B(MinAObatch,MinAObatch,MaxAllowedDimAlphaMPI,MaxAllowedDimGammaMPI,&
          & nbasis,nvirt,nOccBatchDimImax,nOccBatchDimJmax,maxsize)
  IF(maxsize.GT.MemoryAvailable)THEN
     call lsquit('get_optimal_batch_sizes_for_canonical_mp2B Error MinAoBatches',-1)        
  ENDIF

  !assume 
  MaxAllowedDimGamma = MinAObatch  
  !find MaxAllowedDimAlpha as big as possible
  DO I=MaxAllowedDimAlphaMPI,MinAObatch,-1  
     call memestimateCANONMP2B(I,MinAObatch,MaxAllowedDimAlphaMPI,MaxAllowedDimGammaMPI,nbasis,&
          & nvirt,nOccBatchDimImax,nOccBatchDimJmax,maxsize)
     MaxAllowedDimAlpha = I
     IF(maxsize.LT.MemoryAvailable)EXIT
     IF(I.EQ.MinAObatch)THEN
        call lsquit('get_optimal_batch_sizes_for_canonical_mp2B Error MinAoBatchAlpha',-1)        
     ENDIF
  ENDDO

  DO I=MaxAllowedDimGammaMPI,MinAObatch,-1
     call memestimateCANONMP2B(MaxAllowedDimAlpha,I,MaxAllowedDimAlphaMPI,MaxAllowedDimGammaMPI,nbasis,&
          & nvirt,nOccBatchDimImax,nOccBatchDimJmax,maxsize)
     MaxAllowedDimGamma = I
     IF(maxsize.LT.MemoryAvailable)EXIT
     IF(I.EQ.MinAObatch)THEN
        call lsquit('get_optimal_batch_sizes_for_canonical_mp2B Error MinAoBatchGamma',-1)        
     ENDIF
  ENDDO  

end subroutine get_optimal_batch_sizes_for_canonical_mp2B

subroutine memestimateCANONMP2B(MaxAllowedDimAlpha,MaxAllowedDimGamma,dimAlphaMPI,dimGammaMPI,nb,&
          & nvirt,nOccBatchDimImax,nOccBatchDimJmax,maxsize)
implicit none
integer,intent(in) :: MaxAllowedDimAlpha,MaxAllowedDimGamma,dimAlphaMPI,dimGammaMPI,nb
integer,intent(in) :: nvirt,nOccBatchDimImax,nOccBatchDimJmax
real(realk),intent(inout) :: maxsize
real(realk) :: GB,nOccBatchDimI,nOccBI,nOccBatchDimJ,dimAlpha,dimGamma
nOccBatchDimI = nOccBatchDimImax
nOccBI = nOccBatchDimI
nOccBatchDimJ = nOccBatchDimJmax
dimAlpha=MaxAllowedDimAlpha
dimGamma=MaxAllowedDimGamma
!construct CoI(nb,nOccBatchDimI)
maxsize = nb*nOccBI
!call mem_alloc(tmp2,dimAlphaMPI,dimGammaMPI,nb,nOccBatchDimI)
maxsize = nb*nOccBI+dimAlphaMPI*dimGammaMPI*nb*nOccBatchDimI
!BatchGamma: do gammaB = 1,nbatchesGamma
! BatchAlpha: do alphaB = 1,nbatchesAlpha  ! AO batches
!  call mem_alloc(tmp1,dimAlpha*dimGamma,nb*nb)
maxsize = nb*nOccBI+dimAlphaMPI*dimGammaMPI*nb*nOccBatchDimI+dimAlpha*dimGamma*nb*nb
!  call mem_alloc(tmp1b,dimAlpha,dimGamma,nb,nOccBatchDimI)
maxsize = nb*nOccBI+dimAlphaMPI*dimGammaMPI*nb*nOccBatchDimI+dimAlpha*dimGamma*nb*nb+&
     & dimAlpha*dimGamma*nb*nOccBatchDimI
!  call mem_dealloc(tmp1)
!  call mem_dealloc(tmp1b)
! enddo BatchAlpha
!enddo BatchGamma
!call mem_dealloc(CoI)
!call mem_alloc(tmp3,nb*nOccBatchDimI,dimAlphaMPI*dimGammaMPI)
maxsize = dimAlphaMPI*dimGammaMPI*nb*nOccBatchDimI+nb*nOccBatchDimI*dimAlphaMPI*dimGammaMPI !tmp2+tmp3
!call mem_dealloc(tmp2)
!call mem_alloc(tmp4,nvirt*nOccBatchDimI,dimAlphaMPI*dimGammaMPI)
maxsize = nb*nOccBatchDimI*dimAlphaMPI*dimGammaMPI+nvirt*nOccBatchDimI*dimAlphaMPI*dimGammaMPI !tmp3+tmp4
!call mem_dealloc(tmp3)
!call mem_alloc(tmp5,nvirt*nOccBatchDimI,dimAlphaMPI*nOccBatchDimJ)
maxsize = nvirt*nOccBatchDimI*dimAlphaMPI*dimGammaMPI+nvirt*nOccBatchDimI*dimAlphaMPI*nOccBatchDimJ !tmp5+tmp4
!OPTION1
!call mem_alloc(CgammaMPI,dimGammaMPI,nOccBatchDimJ)
maxsize = nvirt*nOccBatchDimI*dimAlphaMPI*dimGammaMPI+nvirt*nOccBatchDimI*dimAlphaMPI*nOccBatchDimJ+& !tmp5+tmp4
     & dimGammaMPI*nOccBatchDimJ
!call mem_dealloc(CgammaMPI)
!OPTION1
!call mem_alloc(tmp6,nvirt*nOccBatchDimI,dimAlpha2*dimGamma2)              
!call mem_alloc(CgammaMPI2,dimGamma2,nOccBatchDimJ)
maxsize = nvirt*nOccBatchDimI*dimAlphaMPI*dimGammaMPI+nvirt*nOccBatchDimI*dimAlphaMPI*nOccBatchDimJ+& !tmp5+tmp4
     & nvirt*nOccBatchDimI*dimAlphaMPI*dimGammaMPI + dimGammaMPI*nOccBatchDimJ !tmp6
!call mem_dealloc(CgammaMPI2)
!call mem_dealloc(tmp6)
!call mem_dealloc(tmp4)
!call mem_alloc(tmp7,dimAlphaMPI*nOccBatchDimJ,nvirt*nOccBatchDimI)
maxsize = nvirt*nOccBatchDimI*dimAlphaMPI*nOccBatchDimJ & 
     & + dimAlphaMPI*nOccBatchDimJ*nvirt*nOccBatchDimI
!call mem_dealloc(tmp5)
!call mem_alloc(VOVO,nvirt*nOccBatchDimJ,nvirt*nOccBatchDimI)
maxsize = dimAlphaMPI*nOccBatchDimJ*nvirt*nOccBatchDimI+nvirt*nOccBatchDimJ*nvirt*nOccBatchDimI !tmp7+VOVOV
!call mem_alloc(CAV,dimAlphaMPI,nvirt)     
maxsize = dimAlphaMPI*nOccBatchDimJ*nvirt*nOccBatchDimI+nvirt*nOccBatchDimJ*nvirt*nOccBatchDimI+& !tmp7+VOVOV
     & dimAlphaMPI*nvirt
!VOVO(nvirt,noccBJ,nvirt,noccBI) = CAV(dimAlpha,nvirt)*tmp7(dimAlpha*nOccBatchDimJ,nvirt*nOccBatchDimI)
!call mem_dealloc(CAV)
!call mem_dealloc(tmp7)
maxsize = nvirt*nOccBatchDimJ*nvirt*nOccBatchDimI*2 !2 VOVO 
GB = 1.000E-9_realk 
maxsize = maxsize*GB
end subroutine memestimateCANONMP2B

!The full_canonical_mp2B requires that (nb,nb,nb) can be distributed across all nodes 
subroutine canonical_mp2B_memreq_test(nbasis,numnodes,Success)
implicit none
integer(kind=ls_mpik),intent(in) :: numnodes
integer,intent(in) :: nbasis
logical,intent(inout) :: Success
!
real(realk) :: MemoryAvailable,GB,nbasisR,numnodesR
! Memory currently available
! **************************
call get_currently_available_memory(MemoryAvailable)
! Note: We multiply by 85 % to be on the safe side!
MemoryAvailable = 0.85E0_realk*MemoryAvailable
GB = 1.000E-9_realk 
nbasisR = nbasis
numnodesR = numnodes
IF(nbasisR*nbasisR*nbasisR*GB/numnodesR.GT.MemoryAvailable)THEN
   print*,'canonical_mp2B_memreq_test  size(nb*nb*nb/numnodes) =',nbasisR*nbasisR*(nbasisR*GB)/numnodesR,' GB'
   print*,'canonical_mp2B_memreq_test  MemoryAvailable         =',MemoryAvailable,' GB'
   call lsquit('canonical_mp2B_memreq_test failure',-1)
ENDIF
Success = .TRUE.
end subroutine canonical_mp2B_memreq_test

subroutine CalcAmat2(nOccBatchDimI,nOccBatchDimJ,nvirt,tmp7,Amat,I,J)
  implicit none
  integer,intent(in) :: nOccBatchDimI,nOccBatchDimJ,nvirt,I,J
  real(realk),intent(in) :: tmp7(nvirt,nOccBatchDimI,nvirt,nOccBatchDimJ)
  real(realk),intent(inout) :: Amat(nvirt,nvirt)
  !local variables
  integer :: A,B
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
  !$OMP PRIVATE(B,A) SHARED(nvirt,tmp7,Amat,I,J)
  do B=1,nvirt
     do A=1,nvirt
        Amat(A,B) = tmp7(A,I,B,J) 
     enddo
  enddo
  !$OMP END PARALLEL DO
end subroutine CalcAmat2

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

    !>   Singles correction
    type(matrix) :: Fic
    type(matrix) :: Fcd
    real(realk)  :: ES2 
    
    !   F12 specific
    real(realk),pointer :: Vijij(:,:)
    real(realk),pointer :: Vjiij(:,:)    

    real(realk),pointer :: Vijij_mp2f12(:,:)
    real(realk),pointer :: Vjiij_mp2f12(:,:)    

    real(realk),pointer :: Vijij_Viajb_term1(:,:)
    real(realk),pointer :: Vijij_Viajb_term2(:,:)
    real(realk),pointer :: Vijij_Viajb_term3(:,:)
    real(realk),pointer :: Vijij_Viajb_term4(:,:)
    real(realk),pointer :: Vijij_Viajb_term5(:,:)

    real(realk),pointer :: Vjiij_Viajb_term1(:,:)
    real(realk),pointer :: Vjiij_Viajb_term2(:,:)
    real(realk),pointer :: Vjiij_Viajb_term3(:,:)
    real(realk),pointer :: Vjiij_Viajb_term4(:,:)  
    real(realk),pointer :: Vjiij_Viajb_term5(:,:)  

    real(realk),pointer :: Vijij_Viija_term1(:,:)
    real(realk),pointer :: Vijij_Viija_term2(:,:)
    real(realk),pointer :: Vijij_Viija_term3(:,:)
    real(realk),pointer :: Vijij_Viija_term4(:,:)
    real(realk),pointer :: Vijij_Viija_all(:,:)
    
    real(realk),pointer :: Vjiij_Viaji_term1(:,:)
    real(realk),pointer :: Vjiij_Viaji_term2(:,:)
    real(realk),pointer :: Vjiij_Viaji_term3(:,:)
    real(realk),pointer :: Vjiij_Viaji_term4(:,:)
    real(realk),pointer :: Vjiij_Viaji_all(:,:)

    real(realk),pointer :: Vijij_Viajj_term1(:,:)
    real(realk),pointer :: Vijij_Viajj_term2(:,:)
    real(realk),pointer :: Vijij_Viajj_term3(:,:)
    real(realk),pointer :: Vijij_Viajj_term4(:,:)
    real(realk),pointer :: Vijij_Viajj_all(:,:)
    
    real(realk),pointer :: Vjiij_Vijja_term1(:,:)
    real(realk),pointer :: Vjiij_Vijja_term2(:,:)
    real(realk),pointer :: Vjiij_Vijja_term3(:,:)
    real(realk),pointer :: Vjiij_Vijja_term4(:,:)
    real(realk),pointer :: Vjiij_Vijja_all(:,:)
     
    real(realk),pointer :: Xijkl(:,:,:,:)
    real(realk),pointer :: Xijij(:,:)
    real(realk),pointer :: Xjiij(:,:)
    real(realk),pointer :: Bijij(:,:)
    real(realk),pointer :: Bjiij(:,:)

    ! CCSD specific
    real(realk),pointer :: Rapbq(:,:,:,:)
    real(realk),pointer :: Rambc(:,:,:,:)
    real(realk),pointer :: Fiajb(:,:,:,:)

    real(realk),pointer :: Viajb(:,:,:,:)
    real(realk),pointer :: Viajb_term1(:,:,:,:)
    real(realk),pointer :: Viajb_term2(:,:,:,:)
    real(realk),pointer :: Viajb_term3(:,:,:,:)
    real(realk),pointer :: Viajb_term4(:,:,:,:)

    real(realk),pointer :: Ripaq(:,:,:,:)
    real(realk),pointer :: Rimac(:,:,:,:)
    real(realk),pointer :: Ramic(:,:,:,:)
    real(realk),pointer :: Fijka(:,:,:,:)

    real(realk),pointer :: Viija(:,:,:)
    real(realk),pointer :: Viija_term1(:,:,:)
    real(realk),pointer :: Viija_term2(:,:,:)
    real(realk),pointer :: Viija_term3(:,:,:)
    real(realk),pointer :: Viija_term4(:,:,:)
    real(realk),pointer :: Viija_full(:,:,:)

    real(realk),pointer :: Viaji(:,:,:)
    real(realk),pointer :: Viaji_term1(:,:,:)
    real(realk),pointer :: Viaji_term2(:,:,:)
    real(realk),pointer :: Viaji_term3(:,:,:)
    real(realk),pointer :: Viaji_term4(:,:,:)
    real(realk),pointer :: Viaji_full(:,:,:)
    
    real(realk),pointer :: Viajj(:,:,:)
    real(realk),pointer :: Viajj_term1(:,:,:)
    real(realk),pointer :: Viajj_term2(:,:,:)
    real(realk),pointer :: Viajj_term3(:,:,:)
    real(realk),pointer :: Viajj_term4(:,:,:)
    real(realk),pointer :: Viajj_full(:,:,:)

    real(realk),pointer :: Vijja(:,:,:)
    real(realk),pointer :: Vijja_term1(:,:,:)
    real(realk),pointer :: Vijja_term2(:,:,:)
    real(realk),pointer :: Vijja_term3(:,:,:)
    real(realk),pointer :: Vijja_term4(:,:,:)
    real(realk),pointer :: Vijja_full(:,:,:)

    !    logical :: fulldriver 
    !    fulldriver = .TRUE.
    !    call init_cabs(fulldriver)

    real(realk) :: tmp
    real(realk) :: E21_Viajb, E21_Viija, E21_Viajj 

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
         & nocc,noccfull,nvirt,ncabs,HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp,Fic,Fcd)

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

    call mem_alloc(Vijij,nocc,nocc)
    call mem_alloc(Vjiij,nocc,nocc)

    ! INTEGRAL GUYS: DO YOUR MAGIC FOR THE F12 CONTRIBUTION HERE...
    call mp2f12_Vijij(Vijij,Ripjq,Gipjq,Fijkl,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)
    call mp2f12_Vjiij(Vjiij,Ripjq,Gipjq,Fijkl,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)

    call mem_alloc(Ciajb,nocc,nvirt,nocc,nvirt)
    call mem_alloc(Cjaib,nocc,nvirt,nocc,nvirt)

    call mp2f12_Ciajb(Ciajb,Giajc,Fac%elms,nocc,nvirt,ncabs)
   
    !> --------------------------------------
    !>               Viajb ccsd
    !>---------------------------------------
    ! V_ij^ab
    call mem_alloc(Viajb,nocc,nvirt,nocc,nvirt) 
    call ccsdf12_Viajb(Viajb,Rapbq,Gipjq,Fiajb,Rambc,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
    
    if(DECinfo%F12DEBUG) then
       call mem_alloc(Viajb_term1,nocc,nvirt,nocc,nvirt) 
       call mem_alloc(Viajb_term2,nocc,nvirt,nocc,nvirt) 
       call mem_alloc(Viajb_term3,nocc,nvirt,nocc,nvirt) 
       call mem_alloc(Viajb_term4,nocc,nvirt,nocc,nvirt) 

       call ccsdf12_Viajb_term1(Viajb_term1,Rapbq,Gipjq,Fiajb,Rambc,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
       call ccsdf12_Viajb_term2(Viajb_term2,Rapbq,Gipjq,Fiajb,Rambc,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
       call ccsdf12_Viajb_term3(Viajb_term3,Rapbq,Gipjq,Fiajb,Rambc,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
       call ccsdf12_Viajb_term4(Viajb_term4,Rapbq,Gipjq,Fiajb,Rambc,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)

       call mem_alloc(Viija_term1,nocc,nocc,nvirt) 
       call mem_alloc(Viija_term2,nocc,nocc,nvirt) 
       call mem_alloc(Viija_term3,nocc,nocc,nvirt) 
       call mem_alloc(Viija_term4,nocc,nocc,nvirt) 
       call mem_alloc(Viija_full,nocc,nocc,nvirt) 

       call mem_alloc(Viaji_term1,nocc,nvirt,nocc) 
       call mem_alloc(Viaji_term2,nocc,nvirt,nocc) 
       call mem_alloc(Viaji_term3,nocc,nvirt,nocc) 
       call mem_alloc(Viaji_term4,nocc,nvirt,nocc)
       call mem_alloc(Viaji_full,nocc,nvirt,nocc)   
       
       Viija_term1 = 0.0E0_realk
       Viaji_term1 = 0.0E0_realk
       call ccsdf12_Viija0(Viija_term1,Fijka,nocc,nvirt)
       call ccsdf12_Viaji0(Viaji_term1,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)      

       Viija_term2 = 0.0E0_realk
       Viaji_term2 = 0.0E0_realk
       call ccsdf12_Viija1(Viija_term2,Ripaq,Gipjq,nocc,noccfull,nbasis,ncabs,nvirt)
       call ccsdf12_Viaji1(Viaji_term2,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)      

       Viija_term3 = 0.0E0_realk
       Viaji_term3 = 0.0E0_realk     
       call ccsdf12_Viija2(Viija_term3,Rimac,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
       call ccsdf12_Viaji2(Viaji_term3,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)      

       Viija_term4 = 0.0E0_realk
       Viaji_term4 = 0.0E0_realk
       call ccsdf12_Viija3(Viija_term4,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
       call ccsdf12_Viaji3(Viaji_term4,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)      

       Viija_full = 0.0E0_realk
       Viaji_full = 0.0E0_realk
       call ccsdf12_Viija(Viija_full,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
       call ccsdf12_Viaji(Viaji_full,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)

       !> Laster 1/3 of the terms

       call mem_alloc(Viajj_term1,nocc,nvirt,nocc) 
       call mem_alloc(Viajj_term2,nocc,nvirt,nocc) 
       call mem_alloc(Viajj_term3,nocc,nvirt,nocc) 
       call mem_alloc(Viajj_term4,nocc,nvirt,nocc) 
       call mem_alloc(Viajj_full, nocc,nvirt,nocc) 

       call mem_alloc(Vijja_term1,nocc,nocc,nvirt) 
       call mem_alloc(Vijja_term2,nocc,nocc,nvirt) 
       call mem_alloc(Vijja_term3,nocc,nocc,nvirt) 
       call mem_alloc(Vijja_term4,nocc,nocc,nvirt)
       call mem_alloc(Vijja_full, nocc,nocc,nvirt)
       
       call ccsdf12_Viajj0(Viajj_term1,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
       call ccsdf12_Viajj1(Viajj_term2,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
       call ccsdf12_Viajj2(Viajj_term3,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
       call ccsdf12_Viajj3(Viajj_term4,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
       
       call ccsdf12_Vijja0(Vijja_term1,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
       call ccsdf12_Vijja1(Vijja_term2,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
       call ccsdf12_Vijja2(Vijja_term3,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
       call ccsdf12_Vijja3(Vijja_term4,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
       
       Viajj_full = 0.0E0_realk
       Vijja_full = 0.0E0_realk

       call ccsdf12_Viajj(Viajj_full,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)
       call ccsdf12_Vijja(Vijja_full,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)

    endif

    ! V_ij^ia
    call mem_alloc(Viija,nocc,nocc,nvirt)
    call ccsdf12_Viija(Viija,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)

   ! V_ij^ai
    call mem_alloc(Viaji,nocc,nvirt,nocc)
    call ccsdf12_Viaji(Viaji,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)

    ! V_ij^aj
    call mem_alloc(Viajj,nocc,nvirt,nocc)
    call ccsdf12_Viajj(Viajj,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)

    ! V_ij^ja
    call mem_alloc(Vijja,nocc,nocc,nvirt)
    call ccsdf12_Vijja(Vijja,Ripaq,Gipjq,Fijka,Rimac,Ramic,Gimjc,nocc,noccfull,nbasis,ncabs,nvirt)

    call ccsdf12_Vijij_coupling(Vijij,Ciajb,Taibj%elm4,Viajb,Viija,Viajj,Tai%elm2,nocc,nvirt)
    call ccsdf12_Vjiij_coupling(Vjiij,Ciajb,Taibj%elm4,Viajb,Vijja,Viaji,Tai%elm2,nocc,nvirt)

    !> Debug Coupling Terms
    if(DECinfo%F12DEBUG) then

       call mem_alloc(Vijij_Viajb_term1,nocc,nocc)
       call mem_alloc(Vijij_Viajb_term2,nocc,nocc)
       call mem_alloc(Vijij_Viajb_term3,nocc,nocc)
       call mem_alloc(Vijij_Viajb_term4,nocc,nocc)
       call mem_alloc(Vijij_Viajb_term5,nocc,nocc)

       call mem_alloc(Vjiij_Viajb_term1,nocc,nocc)
       call mem_alloc(Vjiij_Viajb_term2,nocc,nocc)
       call mem_alloc(Vjiij_Viajb_term3,nocc,nocc)
       call mem_alloc(Vjiij_Viajb_term4,nocc,nocc)
       call mem_alloc(Vjiij_Viajb_term5,nocc,nocc)
       
       call  ccsdf12_Viajb_coupling(Vijij_Viajb_term1,Vjiij_Viajb_term1,Viajb_term1,Taibj%elm4,nocc,nvirt)
       call  ccsdf12_Viajb_coupling(Vijij_Viajb_term2,Vjiij_Viajb_term2,Viajb_term2,Taibj%elm4,nocc,nvirt)
       call  ccsdf12_Viajb_coupling(Vijij_Viajb_term3,Vjiij_Viajb_term3,Viajb_term3,Taibj%elm4,nocc,nvirt)
       call  ccsdf12_Viajb_coupling(Vijij_Viajb_term4,Vjiij_Viajb_term4,Viajb_term4,Taibj%elm4,nocc,nvirt)
       call  ccsdf12_Viajb_coupling(Vijij_Viajb_term5,Vjiij_Viajb_term5,Ciajb,Taibj%elm4,nocc,nvirt)

!!$       print *, '----------------------------------------'
!!$       print *, '  Taibj                                 '
!!$       print *, '----------------------------------------'
!!$       DO i=1,nocc
!!$          DO j=1,nocc
!!$             DO a=1,nvirt
!!$                DO b=1,nvirt
!!$                      print *, "a b i j value: ", a,b,i,j, Taibj%elm4(a,i,b,j) 
!!$                ENDDO
!!$             ENDDO
!!$          ENDDO
!!$       ENDDO
!!$
!!$       print *, '----------------------------------------'
!!$       print *, '  Tai                                 '
!!$       print *, '----------------------------------------'
!!$       DO i=1,nocc
!!$          DO a=1,nvirt
!!$             print *, "a i  value: ", a,i, Tai%velm2(a,i) 
!!$          ENDDO
!!$       ENDDO

       print *, '----------------------------------------'
       print *, ' E21 CCSD Viajb-terms (DEBUG)          '
       print *, '----------------------------------------'
       print *, ' E21_Viajb_term1: ', 2.0E0_REALK*mp2f12_E21(Vijij_Viajb_term1, Vjiij_Viajb_term1, nocc)
       print *, ' E21_Viajb_term2: ', 2.0E0_REALK*mp2f12_E21(Vijij_Viajb_term2, Vjiij_Viajb_term2, nocc)
       print *, ' E21_Viajb_term3: ', 2.0E0_REALK*mp2f12_E21(Vijij_Viajb_term3, Vjiij_Viajb_term3, nocc)
       print *, ' E21_Viajb_term4: ', 2.0E0_REALK*mp2f12_E21(Vijij_Viajb_term4, Vjiij_Viajb_term4, nocc)
       print *, ' E21_Viajb_term5: ', 2.0E0_REALK*mp2f12_E21(Vijij_Viajb_term5, Vjiij_Viajb_term5, nocc)
       print *, '----------------------------------------'
       print *, ' sum: ', 2.0E0_REALK*(mp2f12_E21(Vijij_Viajb_term1, Vjiij_Viajb_term1, nocc) + & 
            & mp2f12_E21(Vijij_Viajb_term2, Vjiij_Viajb_term2, nocc) + &
            & mp2f12_E21(Vijij_Viajb_term3, Vjiij_Viajb_term3, nocc) + mp2f12_E21(Vijij_Viajb_term4, Vjiij_Viajb_term4, nocc) + &
            & mp2f12_E21(Vijij_Viajb_term5, Vjiij_Viajb_term5, nocc))

       E21_Viajb = 2.0E0_REALK*(mp2f12_E21(Vijij_Viajb_term1, Vjiij_Viajb_term1, nocc) + & 
            & mp2f12_E21(Vijij_Viajb_term2, Vjiij_Viajb_term2, nocc) + &
            & mp2f12_E21(Vijij_Viajb_term3, Vjiij_Viajb_term3, nocc) + mp2f12_E21(Vijij_Viajb_term4, Vjiij_Viajb_term4, nocc) + &
            & mp2f12_E21(Vijij_Viajb_term5, Vjiij_Viajb_term5, nocc))

       call mem_alloc(Vijij_Viija_term1,nocc,nocc)
       call mem_alloc(Vijij_Viija_term2,nocc,nocc)
       call mem_alloc(Vijij_Viija_term3,nocc,nocc)
       call mem_alloc(Vijij_Viija_term4,nocc,nocc)

       call mem_alloc(Vjiij_Viaji_term1,nocc,nocc)
       call mem_alloc(Vjiij_Viaji_term2,nocc,nocc)
       call mem_alloc(Vjiij_Viaji_term3,nocc,nocc)
       call mem_alloc(Vjiij_Viaji_term4,nocc,nocc)

       call  ccsdf12_Viija_coupling(Vijij_Viija_term1,Vjiij_Viaji_term1,Viija_term1,Viaji_term1,Tai%elm2,nocc,nvirt)
       call  ccsdf12_Viija_coupling(Vijij_Viija_term2,Vjiij_Viaji_term2,Viija_term2,Viaji_term2,Tai%elm2,nocc,nvirt)
       call  ccsdf12_Viija_coupling(Vijij_Viija_term3,Vjiij_Viaji_term3,Viija_term3,Viaji_term3,Tai%elm2,nocc,nvirt)
       call  ccsdf12_Viija_coupling(Vijij_Viija_term4,Vjiij_Viaji_term4,Viija_term4,Viaji_term4,Tai%elm2,nocc,nvirt)

       !> Debug Test
       call mem_alloc(Vijij_Viija_all,nocc,nocc)
       call mem_alloc(Vjiij_Viaji_all,nocc,nocc)

       call ccsdf12_Viija_coupling(Vijij_Viija_all,Vjiij_Viaji_all,Viija_full,Viaji_full,Tai%elm2,nocc,nvirt)
     
       print *, '----------------------------------------'
       print *, ' E21 CCSD Viija_terms (DEBUG)          '
       print *, '----------------------------------------'
       print *, ' E21_Viija_term1: ', 2.0E0_REALK*mp2f12_E21(Vijij_Viija_term1, Vjiij_Viaji_term1, nocc)
       print *, ' E21_Viija_term2: ', 2.0E0_REALK*mp2f12_E21(Vijij_Viija_term2, Vjiij_Viaji_term2, nocc)     
       print *, ' E21_Viija_term3: ', 2.0E0_REALK*mp2f12_E21(Vijij_Viija_term3, Vjiij_Viaji_term3, nocc)     
       print *, ' E21_Viija_term4: ', 2.0E0_REALK*mp2f12_E21(Vijij_Viija_term4, Vjiij_Viaji_term4, nocc)     
       print *, '----------------------------------------'
       print *, ' sum: ', 2.0E0_REALK*(mp2f12_E21(Vijij_Viija_term1, Vjiij_Viaji_term1, nocc) + &
            & mp2f12_E21(Vijij_Viija_term2, Vjiij_Viaji_term2, nocc) + &
            & mp2f12_E21(Vijij_Viija_term3, Vjiij_Viaji_term3, nocc) + mp2f12_E21(Vijij_Viija_term4, Vjiij_Viaji_term4, nocc))
       print *, '----------------------------------------'
       print *, ' debug: ', 2.0E0_REALK*mp2f12_E21(Vijij_Viija_all,Vjiij_Viaji_all, nocc)
       
       E21_Viija =  2.0E0_REALK*(mp2f12_E21(Vijij_Viija_term1, Vjiij_Viaji_term1, nocc) + &
            & mp2f12_E21(Vijij_Viija_term2, Vjiij_Viaji_term2, nocc) + &
            & mp2f12_E21(Vijij_Viija_term3, Vjiij_Viaji_term3, nocc) + mp2f12_E21(Vijij_Viija_term4, Vjiij_Viaji_term4, nocc))

       !> Last term check
       call mem_alloc(Vijij_Viajj_term1,nocc,nocc)
       call mem_alloc(Vijij_Viajj_term2,nocc,nocc)
       call mem_alloc(Vijij_Viajj_term3,nocc,nocc)
       call mem_alloc(Vijij_Viajj_term4,nocc,nocc)

       call mem_alloc(Vjiij_Vijja_term1,nocc,nocc)
       call mem_alloc(Vjiij_Vijja_term2,nocc,nocc)
       call mem_alloc(Vjiij_Vijja_term3,nocc,nocc)
       call mem_alloc(Vjiij_Vijja_term4,nocc,nocc)

       call ccsdf12_Viajj_coupling(Vijij_Viajj_term1,Vjiij_Vijja_term1,Viajj_term1,Vijja_term1,Tai%elm2,nocc,nvirt)
       call ccsdf12_Viajj_coupling(Vijij_Viajj_term2,Vjiij_Vijja_term2,Viajj_term2,Vijja_term2,Tai%elm2,nocc,nvirt)
       call ccsdf12_Viajj_coupling(Vijij_Viajj_term3,Vjiij_Vijja_term3,Viajj_term3,Vijja_term3,Tai%elm2,nocc,nvirt)
       call ccsdf12_Viajj_coupling(Vijij_Viajj_term4,Vjiij_Vijja_term4,Viajj_term4,Vijja_term4,Tai%elm2,nocc,nvirt)
     
       !> Debug Test
      call mem_alloc(Vijij_Viajj_all,nocc,nocc)
      call mem_alloc(Vjiij_Vijja_all,nocc,nocc)

       call  ccsdf12_Viajj_coupling(Vijij_Viajj_all,Vjiij_Vijja_all,Viajj_full,Vijja_full,Tai%elm2,nocc,nvirt)

       print *, '----------------------------------------'
       print *, '  E21 CCSD Viajj_terms (DEBUG)          '
       print *, '----------------------------------------'
       print *, 'E21_Viajj_term1: ', 2.0E0_REALK*mp2f12_E21(Vijij_Viajj_term1, Vjiij_Vijja_term1, nocc)
       print *, 'E21_Viajj_term2: ', 2.0E0_REALK*mp2f12_E21(Vijij_Viajj_term2, Vjiij_Vijja_term2, nocc)     
       print *, 'E21_Viajj_term3: ', 2.0E0_REALK*mp2f12_E21(Vijij_Viajj_term3, Vjiij_Vijja_term3, nocc)     
       print *, 'E21_Viajj_term4: ', 2.0E0_REALK*mp2f12_E21(Vijij_Viajj_term4, Vjiij_Vijja_term4, nocc)     
       print *, '----------------------------------------'
       print *, 'sum: ', 2.0E0_REALK*(mp2f12_E21(Vijij_Viajj_term1, Vjiij_Vijja_term1, nocc) + &
            & mp2f12_E21(Vijij_Viajj_term2, Vjiij_Vijja_term2, nocc) + &
            & mp2f12_E21(Vijij_Viajj_term3, Vjiij_Vijja_term3, nocc) + mp2f12_E21(Vijij_Viajj_term4, Vjiij_Vijja_term4, nocc))
       print *, '----------------------------------------'
       print *, 'debug: ', 2.0E0_REALK*mp2f12_E21(Vijij_Viajj_all,Vjiij_Vijja_all, nocc)

       E21_Viajj = 2.0E0_REALK*(mp2f12_E21(Vijij_Viajj_term1, Vjiij_Vijja_term1, nocc) + &
            & mp2f12_E21(Vijij_Viajj_term2, Vjiij_Vijja_term2, nocc) + &
            & mp2f12_E21(Vijij_Viajj_term3, Vjiij_Vijja_term3, nocc) + mp2f12_E21(Vijij_Viajj_term4, Vjiij_Vijja_term4, nocc))

    endif

    ! CCSD Specific
    call mem_dealloc(Viija)
    call mem_dealloc(Vijja)
    call mem_dealloc(Viajj)
    call mem_dealloc(Viaji)
    call mem_dealloc(Viajb)

    if(DECinfo%F12DEBUG) then
       call mem_dealloc(Viajb_term1)
       call mem_dealloc(Viajb_term2)
       call mem_dealloc(Viajb_term3)
       call mem_dealloc(Viajb_term4)

       call mem_dealloc(Viija_term1)
       call mem_dealloc(Viija_term2)
       call mem_dealloc(Viija_term3)
       call mem_dealloc(Viija_term4)
       call mem_dealloc(Viija_full)

       call mem_dealloc(Viaji_term1)
       call mem_dealloc(Viaji_term2)
       call mem_dealloc(Viaji_term3)
       call mem_dealloc(Viaji_term4)
       call mem_dealloc(Viaji_full)
       
       call mem_dealloc(Viajj_term1)
       call mem_dealloc(Viajj_term2)
       call mem_dealloc(Viajj_term3)
       call mem_dealloc(Viajj_term4)
       call mem_dealloc(Viajj_full)
       
       call mem_dealloc(Vijja_term1)
       call mem_dealloc(Vijja_term2)
       call mem_dealloc(Vijja_term3)
       call mem_dealloc(Vijja_term4)
       call mem_dealloc(Vijja_full)
    end if

    !> CCSD Specific MP2-F12 energy
    E21 = 2.0E0_realk*mp2f12_E21(Vijij,Vjiij,nocc)

    ! F12 Specific
    call mem_dealloc(Vijij)
    call mem_dealloc(Vjiij)

    if(DECinfo%F12DEBUG) then

       call mem_dealloc(Vijij_Viajb_term1)
       call mem_dealloc(Vijij_Viajb_term2)
       call mem_dealloc(Vijij_Viajb_term3)
       call mem_dealloc(Vijij_Viajb_term4)   
       call mem_dealloc(Vijij_Viajb_term5)   

       call mem_dealloc(Vjiij_Viajb_term1)
       call mem_dealloc(Vjiij_Viajb_term2)
       call mem_dealloc(Vjiij_Viajb_term3)
       call mem_dealloc(Vjiij_Viajb_term4)
       call mem_dealloc(Vjiij_Viajb_term5)

       call mem_dealloc(Vijij_Viija_term1)
       call mem_dealloc(Vijij_Viija_term2)
       call mem_dealloc(Vijij_Viija_term3)
       call mem_dealloc(Vijij_Viija_term4)   
       call mem_dealloc(Vijij_Viija_all)   

       call mem_dealloc(Vjiij_Viaji_term1)
       call mem_dealloc(Vjiij_Viaji_term2)
       call mem_dealloc(Vjiij_Viaji_term3)
       call mem_dealloc(Vjiij_Viaji_term4) 
       call mem_dealloc(Vjiij_Viaji_all)

       call mem_dealloc(Vijij_Viajj_term1)
       call mem_dealloc(Vijij_Viajj_term2)
       call mem_dealloc(Vijij_Viajj_term3)
       call mem_dealloc(Vijij_Viajj_term4)
       call mem_dealloc(Vijij_Viajj_all)

       call mem_dealloc(Vjiij_Vijja_term1)
       call mem_dealloc(Vjiij_Vijja_term2)
       call mem_dealloc(Vjiij_Vijja_term3)
       call mem_dealloc(Vjiij_Vijja_term4) 
       call mem_dealloc(Vjiij_Vijja_all)

    endif

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

    ! CCSD-F12 Singles Correction Energy
    call get_ES2(ES2,Fic,Fii,Fcd,nocc,ncabs)


    call free_F12_mixed_MO_Matrices(HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp,Fic,Fcd)
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

    if(DECinfo%F12Debug) then
       call mem_alloc(Vijij_mp2f12,nocc,nocc)
       call mem_alloc(Vjiij_mp2f12,nocc,nocc)

       call mp2f12_Vijij(Vijij_mp2f12,Ripjq,Gipjq,Fijkl,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)
       call mp2f12_Vjiij(Vjiij_mp2f12,Ripjq,Gipjq,Fijkl,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)
       
       ! CCSD coupling V_ij^ij & V_ji^ij 
       ! call  ccsdf12_Viajb_coupling(Vijij_mp2f12,Vjiij_mp2f12,Ciajb,Taibj%elm4,nocc,nvirt)

       print *, '----------------------------------------'
       print *, ' CCSD-F12 energy summary                '
       print *, '----------------------------------------'
       write(*,'(1X,a,f25.10)') 'E21_Vijij_MP2_noC: ', 2.0E0_REALK*mp2f12_E21(Vijij_mp2f12, Vjiij_mp2f12, nocc)
       write(*,'(1X,a,f25.10)') 'E22+E23_Vijij_MP2: ', E22
       write(*,'(1X,a,f25.10)') 'E21_Vijij_CCSD:    ', E21_Viajb + E21_Viija + E21_Viajj
       print *, '----------------------------------------'
       write(*,'(1X,a,f25.10)') 'sum:               ', 2.0E0_REALK*mp2f12_E21(Vijij_mp2f12, Vjiij_mp2f12, nocc) + & 
            & E22 + E21_Viajb + E21_Viija + E21_Viajj
        
       call mem_dealloc(Vijij_mp2f12)
       call mem_dealloc(Vjiij_mp2f12)

    endif
    
    call mem_dealloc(Ciajb)
    call mem_dealloc(Cjaib)

    ! Add contributions
    ECCSD_F12 = ECCSD + EF12
   
    write(*,'(1X,a)') '---------------------------------------------------------------------'  
    write(*,'(1X,a,f20.10)') 'TOYCODE: CCSD CORRECTION TO ENERGY      =         ', ECCSD
    write(*,'(1X,a,f20.10)') 'TOYCODE: CCSD-F12 CORRECTION TO ENERGY  =         ', EF12
    write(*,'(1X,a,f20.10)') 'TOYCODE: CCSD-F12 CORRELATION ENERGY    =         ', ECCSD_F12
    write(*,'(1X,a)') '---------------------------------------------------------------------'  
    write(*,'(1X,a,f20.10)') 'TOYCODE: CCSD-F12 SINGLES CORRECTION TO ENERGY  = ', ES2
    write(*,'(1X,a,f20.10)') 'TOYCODE: FULL F12 CORRECTION TO ENERGY          = ', EF12 + ES2
    write(*,'(1X,a,f20.10)') 'TOYCODE: FULL CCSD-F12 CORRECTION TO ENERGY     = ', ECCSD_F12 + ES2
   
    !> Input to DEC
    write(DECinfo%output,*) 'TOYCODE: CCSD CORRECTION TO ENERGY      = ', ECCSD
    write(DECinfo%output,*) 'TOYCODE: CCSD-F12 CORRECTION TO ENERGY  = ', EF12
    write(DECinfo%output,*) 'TOYCODE: CCSD-F12 CORRELATION ENERGY    = ', ECCSD_F12

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
    integer :: solver_ccmodel
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
    solver_ccmodel = DECinfo%ccmodel
    if(DECinfo%ccmodel == MODEL_CCSDpT) solver_ccmodel = MODEL_CCSD


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
       call ccsolver_par(solver_ccmodel,MyMolecule%Co(1:nbasis,startidx:endidx),&
            & MyMolecule%Cv,MyMolecule%fock, nbasis,nocc,nunocc,mylsitem,&
            & print_level,energy,&
            & VOVO,.false.,local,SOLVE_AMPLITUDES,p2=Tai,p4=Taibj)
       call mem_dealloc(ppfock)

    else

       call ccsolver_par(solver_ccmodel,MyMolecule%Co,MyMolecule%Cv,&
            & MyMolecule%fock, nbasis,nocc,nunocc,mylsitem, print_level, &
            & energy,VOVO,.false.,local,SOLVE_AMPLITUDES,p2=Tai,p4=Taibj)

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

subroutine full_canonical_mp2_slave
  use full,only: full_canonical_mp2,full_canonical_mp2B
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
  real(realk) :: mp2_energy    
  logical :: Success
  
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
  call full_canonical_mp2B(MyMolecule,MyLsitem,mp2_energy)
  call ls_free(MyLsitem)
  call molecule_finalize(MyMolecule)
  
end subroutine full_canonical_mp2_slave
#endif

