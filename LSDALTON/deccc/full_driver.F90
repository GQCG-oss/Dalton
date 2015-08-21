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
  use fullrimp2 !,only: full_canonical_rimp2
  use fullrimp2f12 !,only: full_canonical_rimp2_f12
  use fullmp2 
  use full_mp3_module
  use full_ls_thc_rimp2Mod

  public  :: full_driver
  private :: mp2f12_E22X

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
    real(realk) :: Ecorr_rimp2, Ecorr_rimp2f12
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
!When the code is a production code it should be released! TK
          if(DECinfo%ccModel==MODEL_MP2) then
             call full_canonical_mp2_f12(MyMolecule,MyLsitem,D,Ecorr)
          elseif(DECinfo%ccModel==MODEL_RIMP2) then
             call full_canonical_rimp2(MyMolecule,MyLsitem,Ecorr_rimp2)       
             Ecorr_rimp2f12 = Ecorr_rimp2 !in order to provide full_canonical_rimp2_f12 the MP2 energy
             call full_canonical_rimp2_f12(MyMolecule,MyLsitem,D,Ecorr_rimp2f12)
             Ecorr = Ecorr_rimp2 + Ecorr_rimp2f12
             write(DECinfo%output,'(/,a)') ' ================================================ '
             write(DECinfo%output,'(a)')   '                 Energy Summary                   '
             write(DECinfo%output,'(a,/)') ' ================================================ '
             write(*,'(/,a)') ' ================================================ '
             write(*,'(a)')   '                 Energy Summary                   '
             write(*,'(a,/)') ' ================================================ '
             write(*,'(1X,a,f20.10)') 'TOYCODE: RI-MP2 CORRECTION TO ENERGY =    ', Ecorr_rimp2
             write(DECinfo%output,'(1X,a,f20.10)')  'TOYCODE: RI-MP2 CORRECTION TO ENERGY =    ', Ecorr_rimp2
             write(*,'(1X,a,f20.10)') 'TOYCODE: RI-MP2F12 CORRECTION TO ENERGY = ', Ecorr_rimp2f12
             write(DECinfo%output,'(1X,a,f20.10)')  'TOYCODE: RI-MP2F12 CORRECTION TO ENERGY = ', Ecorr_rimp2f12
          else
             call full_get_ccsd_f12_energy(MyMolecule,MyLsitem,D,Ecorr)
          end if
#else
          call lsquit('f12 not released',-1)
#endif
       elseif(DECinfo%ccModel==MODEL_RIMP2)then
          !       call lsquit('RIMP2 currently not implemented for **CC ',-1)
          call full_canonical_rimp2(MyMolecule,MyLsitem,Ecorr)       
       elseif(DECinfo%ccModel==MODEL_LSTHCRIMP2)then
          call full_canonical_ls_thc_rimp2(MyMolecule,MyLsitem,Ecorr) 
       elseif(DECinfo%ccmodel==MODEL_MP3) then
          call full_canonical_mp3(MyMolecule,MyLsitem,Ecorr)
       else
          if(DECinfo%ccModel==MODEL_MP2) then
             if(DECinfo%use_canonical .and. (.not. DECinfo%CCexci) ) then
                ! Do not use canonical MP2 routine when CC eigenvalues are requested
                ! because we need output amplitudes.

                !simple conventional MP2 calculation only works for canonical orbitals
                !no amplitudes stored. MP2B requires (nb,nb,nb) can be fully distributed
                IF(DECinfo%MPMP2)THEN
                   call full_canonical_mpmp2(MyMolecule,MyLsitem,Ecorr)
                ELSE
                   call full_canonical_mp2(MyMolecule,MyLsitem,Ecorr)
                ENDIF
             else
                !Call routine which calculates individual fragment 
                !contributions and prints them,
                !works both for canonical and local orbitals
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
       call print_total_energy_summary(EHF,Edft,Ecorr,MyMolecule%EF12singles,&
            & 0.0E0_realk,0.0E0_realk,0.0E0_realk)

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
    integer :: nocc,nvirt,nbasis,print_level,i,j
    logical :: fragment_job
    real(realk),pointer :: ppfock_fc(:,:), Co_fc(:,:)

    if(MyMolecule%mem_distributed)then
       call lsquit("ERROR(full_cc_dispatch): does not work with distributed full molecule",-1)
    endif

    Ecorr = 0.0E0_realk

    if(DECinfo%FrozenCore) then
       nocc = MyMolecule%nval
    else
       nocc = MyMolecule%nocc
    end if
    nvirt = MyMolecule%nvirt
    nbasis = MyMolecule%nbasis

    fragment_job = .false.
    print_level = 2 ! this is not used

    if(DECinfo%frozencore) then
       ! Pass only valence orbitals
       call mem_alloc(ppfock_fc,nocc,nocc)
       do j=1,nocc
          do i=1,nocc
             ppfock_fc(i,j) = MyMolecule%oofock%elm2(MyMolecule%ncore+i,MyMolecule%ncore+j)
          end do
       end do

       ! Frozen core component of MO coeff.
       call mem_alloc(Co_fc,nbasis,nocc)
       do j=1,nocc
          do i=1,nbasis
             Co_fc(i,j) = MyMolecule%Co%elm2(i,MyMolecule%ncore+j)
          end do
       end do

       Ecorr = ccsolver_justenergy(DECinfo%ccmodel,MyMolecule,nbasis,nocc,nvirt,&
          & mylsitem,print_level,fragment_job,Co_fc=Co_fc,ppfock_fc=ppfock_fc)

       call mem_dealloc(ppfock_fc)
       call mem_dealloc(Co_fc)

    else


#ifdef MOD_UNRELEASED
       if(DECinfo%CCSDmultipliers)then
          call ccsolver_energy_multipliers(DECinfo%ccmodel,MyMolecule%Co%elm2,MyMolecule%Cv%elm2,&
             & MyMolecule%fock%elm2, nbasis,nocc,nvirt,mylsitem, &
             & print_level,fragment_job,MyMolecule%oofock%elm2,MyMolecule%vvfock%elm2,ecorr)
       else
          Ecorr = ccsolver_justenergy(DECinfo%ccmodel,MyMolecule,nbasis,nocc,nvirt,&
             & mylsitem,print_level,fragment_job)
       endif
#else
       Ecorr = ccsolver_justenergy(DECinfo%ccmodel,MyMolecule,nbasis,nocc,nvirt,&
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
    type(matrix) :: Fab
        
    !> Singles correction energy
    real(realk)  :: ES2
     
    real(realk)  :: E21, E21_C, E21_debug, E22, E22_debug, E23_debug, Gtmp
    type(tensor) :: tensor_Taibj,tensor_gmo
    integer :: vs, os,offset
    logical :: local

    !> Additional
    real(realk) :: EK
    real(realk) :: EJ
    EK = 0E0_realk
    EJ = 0E0_realk

    local = .true.

#ifdef VAR_MPI
    local = (infpar%lg_nodtot==1)
#endif
    !    logical :: fulldriver 
    !    fulldriver = .TRUE.
    !    call init_cabs(fulldriver)

    if(MyMolecule%mem_distributed)then
       call lsquit("ERROR(full_canonical_mp2_f12): does not work with PDM type fullmolecule",-1)
    endif

    ! Init stuff
    ! **********
    nbasis = MyMolecule%nbasis
    nvirt  = MyMolecule%nvirt
    call determine_CABS_nbast(ncabsAO,ncabs,mylsitem%setting,DECinfo%output)

    ! Set number of occupied orbitals
    if(DECinfo%frozencore) then
       ! Frozen core: nocc = #valence orbitals
       nocc = MyMolecule%nval
    else
       ! Not frozen core: nocc = total number of occ orbitals
       nocc = MyMolecule%nocc
    end if
    ! noccfull: Always equal to total number of occ orbitals
    noccfull = MyMolecule%nocc

    ! Offset:   Frozen core    : ncore
    !           Not frozen core: 0
    offset = noccfull - nocc

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
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'aiai',gAO,gMO)
    call mem_dealloc(gao)

    call get_4Center_F12_integrals(mylsitem,MyMolecule,nbasis,nocc,noccfull,nvirt,ncabsAO,&
         & Ripjq,Fijkl,Tijkl,Rimjc,Dijkl,Tirjk,Tijkr,Gipjq,Gimjc,Girjs,Girjm,&
         & Grimj,Gipja,Gpiaj,Gicjm,Gcimj,Gcirj,Gciaj,Giajc)

    call mem_alloc(Vijij,nocc,nocc)
    call mem_alloc(Vjiij,nocc,nocc)

    call mp2f12_Vijij(Vijij,Ripjq,Gipjq,Fijkl,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)
    call mp2f12_Vjiij(Vjiij,Ripjq,Gipjq,Fijkl,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)

    call mem_alloc(Ciajb,nocc,nvirt,nocc,nvirt)    
    call mp2f12_Ciajb(Ciajb,Giajc,Fac%elms,nocc,nvirt,ncabs)
    
    ! MP2-F12 Singles correction 
    ! ***************************    
    ES2 = MyMolecule%EF12singles

    !call mat_free(Fab)

    E21 = 0.0E0_realk
    E21_C = 0.0E0_realk
    E21_debug = 0.0E0_realk
    DoCanonical: if(DECinfo%use_canonical) then
       !construct canonical T amplitudes
       call mem_alloc(Taibj,nvirt,nocc,nvirt,nocc)
       do J=1,nocc
          do B=1,nvirt
             do I=1,nocc
                do A=1,nvirt
                   ! Difference in orbital energies: eps(I) + eps(J) - eps(A) - eps(B)
                   eps = MyMolecule%oofock%elm2(I+offset,I+offset) &
                        & + MyMolecule%oofock%elm2(J+offset,J+offset) &
                        & - MyMolecule%vvfock%elm2(A,A) - MyMolecule%vvfock%elm2(B,B)
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
!                   eps = MyMolecule%oofock%elm2(I+offset,I+offset) &
!                        & + MyMolecule%oofock%elm2(J+offset,J+offset) &
!                        & - MyMolecule%vvfock%elm2(A,A) - MyMolecule%vvfock%elm2(B,B)

                   ! Energy = sum_{AIBJ} (AI|BJ) * [ 2(AI|BJ) - (BI|AJ) ] / (epsI + epsJ - epsA - epsB)
                   mp2_energy = mp2_energy + Taibj(a,i,b,j)*(2E0_realk*gmo(A,I,B,J)-gmo(B,I,A,J))

                end do
             end do
          end do
       end do
       IF(DECinfo%F12Ccoupling)THEN
          !overwrite the amplitudes with F12 modified amplitudes which includes the C coupling. TK
          !Build delta T amplitudes with ONLY C coupling 
          tmp = 0.0E0_realk
          do B=1,nvirt
             do J=1,nocc
                do A=1,nvirt
                   do I=1,nocc
                      ! Difference in orbital energies: eps(I) + eps(J) - eps(A) - eps(B)
                      eps = MyMolecule%oofock%elm2(I+offset,I+offset) + MyMolecule%oofock%elm2(J+offset,J+offset) &
                           & - MyMolecule%vvfock%elm2(A,A) - MyMolecule%vvfock%elm2(B,B)
                      tmp = tmp + (7.0E0_realk*Ciajb(I,A,J,B)*Ciajb(I,A,J,B) + 1.0E0_realk*Ciajb(I,A,J,B)*Ciajb(J,A,I,B))/eps
                   enddo
                enddo
             enddo
          enddo
          E21_C = tmp/32.0E0_realk
          print *, "E21_Cterm ...", E21_C
       ENDIF
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

    endif DoCanonical
  
    call mp2f12_Vijij_coupling(Vijij,Ciajb,Taibj,nocc,nvirt)
    call mp2f12_Vjiij_coupling(Vjiij,Ciajb,Taibj,nocc,nvirt)

    !> Calculate E21 Energy
    E21 = E21 + 2.0E0_REALK*mp2f12_E21(Vijij,Vjiij,nocc)

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
       
       call mp2f12_Vijij_term5(Vijij_term5,Ciajb,Taibj,nocc,nvirt)
       call mp2f12_Vjiij_term5(Vjiij_term5,Ciajb,Taibj,nocc,nvirt)
      
       EJ = 0.0E0_realk
       EK = 0.0E0_realk

       DO j=1,nocc
          DO i=1,nocc
             EJ = EJ + Vijij_term5(i,j) 
             EK = EK + Vjiij_term5(i,j)
          ENDDO
       ENDDO

       !print *, "EJ_V5: ", -5.0/4.0*EJ
       !print *, "EK_V5: ", 1.0/4.0*EK
       !print *, "EK_V5 + EJ_V5: ", -1.0*(5.0/4.0*EJ - 1.0/4.0*EK)

 
       E21_debug = E21_debug + 2.0E0_REALK*(mp2f12_E21(Vijij_term1,Vjiij_term1,nocc) + mp2f12_E21(Vijij_term2,Vjiij_term2,nocc) &
                                        & + mp2f12_E21(Vijij_term3,Vjiij_term3,nocc) + mp2f12_E21(Vijij_term4,Vjiij_term4,nocc) &
                                        & + mp2f12_E21(Vijij_term5,Vjiij_term5,nocc)) 

       IF(DECinfo%F12Ccoupling)THEN
          E21_debug = E21_debug + E21_C
       ENDIF

       ! ***********************************************************
       !   Printing Input variables 
       ! ***********************************************************
       print *, "-------------------------------------------------"
       print *, "     F12-integrals.F90                           "
       print *, "-------------------------------------------------"
       print *, "nbasis:  ", nbasis
       print *, "nocc:    ", nocc
       print *, "nvirt:   ", nvirt
       print *, "-------------------------------------------------"
       print *, "noccfull ", noccfull
       print *, "ncabsAO  ", ncabsAO
       print *, "ncabsMO  ", ncabs

       print *, '----------------------------------------'
       print *, ' E21 V terms                            '
       print *, '----------------------------------------'
       write(*,'(1X,a,g25.16)') ' E21_CC_term: ', E21_C
       write(*,'(1X,a,g25.16)') ' E21_V_term1: ', 2.0E0_REALK*mp2f12_E21(Vijij_term1,Vjiij_term1,nocc)
       write(*,'(1X,a,g25.16)') ' E21_V_term2: ', 2.0E0_REALK*mp2f12_E21(Vijij_term2,Vjiij_term2,nocc)
       write(*,'(1X,a,g25.16)') ' E21_V_term3: ', 2.0E0_REALK*mp2f12_E21(Vijij_term3,Vjiij_term3,nocc)
       write(*,'(1X,a,g25.16)') ' E21_V_term4: ', 2.0E0_REALK*mp2f12_E21(Vijij_term4,Vjiij_term4,nocc)
       write(*,'(1X,a,g25.16)') ' E21_V_term5: ', 2.0E0_REALK*mp2f12_E21(Vijij_term5,Vjiij_term5,nocc)
       print *, '----------------------------------------'
       write(*,'(1X,a,g25.16)') ' E21_Vsum:    ', E21_debug
       !write(*,'(1X,a,g25.16)') ' E21_debug:   ', E21
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

       E22 = 0.0E0_realk
       E22_debug = 0.0E0_realk
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
          write(*,'(1X,a,g25.16)') ' E22_X_term1: ', mp2f12_E22(Xijij_term1,Xjiij_term1,Fii%elms,nocc)
          write(*,'(1X,a,g25.16)') ' E22_X_term2: ', mp2f12_E22(Xijij_term2,Xjiij_term2,Fii%elms,nocc)
          write(*,'(1X,a,g25.16)') ' E22_X_term3: ', mp2f12_E22(Xijij_term3,Xjiij_term3,Fii%elms,nocc)
          write(*,'(1X,a,g25.16)') ' E22_X_term4: ', mp2f12_E22(Xijij_term4,Xjiij_term4,Fii%elms,nocc)
          print *, '----------------------------------------'
          E22_debug = mp2f12_E22(Xijij_term1,Xjiij_term1,Fii%elms,nocc) & 
               & + mp2f12_E22(Xijij_term2,Xjiij_term2,Fii%elms,nocc) &
               & + mp2f12_E22(Xijij_term3,Xjiij_term3,Fii%elms,nocc) + mp2f12_E22(Xijij_term4,Xjiij_term4,Fii%elms,nocc)  
          write(*,'(1X,a,g25.16)') ' E22_Xsum:    ', E22_debug  
          !write(*,'(1X,a,g25.16)') ' E22_debug:   ', mp2f12_E22(Xijij,Xjiij,Fii%elms,nocc)
          print *, '----------------------------------------'
          print *, ' E_23 B term                           '
          print *, '----------------------------------------'
          write(*,'(1X,a,g25.16)') ' E23_B_term1: ', mp2f12_E23(Bijij_term1,Bjiij_term1,nocc)
          write(*,'(1X,a,g25.16)') ' E23_B_term2: ', mp2f12_E23(Bijij_term2,Bjiij_term2,nocc)
          write(*,'(1X,a,g25.16)') ' E23_B_term3: ', mp2f12_E23(Bijij_term3,Bjiij_term3,nocc)
          write(*,'(1X,a,g25.16)') ' E23_B_term4: ', mp2f12_E23(Bijij_term4,Bjiij_term4,nocc)
          write(*,'(1X,a,g25.16)') ' E23_B_term5: ', mp2f12_E23(Bijij_term5,Bjiij_term5,nocc)
          write(*,'(1X,a,g25.16)') ' E23_B_term6: ', mp2f12_E23(Bijij_term6,Bjiij_term6,nocc)
          write(*,'(1X,a,g25.16)') ' E23_B_term7: ', mp2f12_E23(Bijij_term7,Bjiij_term7,nocc)
          write(*,'(1X,a,g25.16)') ' E23_B_term8: ', mp2f12_E23(Bijij_term8,Bjiij_term8,nocc)
          write(*,'(1X,a,g25.16)') ' E23_B_term9: ', mp2f12_E23(Bijij_term9,Bjiij_term9,nocc)   
          print *, '----------------------------------------'
          E23_debug = mp2f12_E23(Bijij_term1,Bjiij_term1,nocc) & 
               & + mp2f12_E23(Bijij_term2,Bjiij_term2,nocc) + mp2f12_E23(Bijij_term3,Bjiij_term3,nocc) &
               & + mp2f12_E23(Bijij_term4,Bjiij_term4,nocc) + mp2f12_E23(Bijij_term5,Bjiij_term5,nocc) &
               & + mp2f12_E23(Bijij_term6,Bjiij_term6,nocc) + mp2f12_E23(Bijij_term7,Bjiij_term7,nocc) &
               & + mp2f12_E23(Bijij_term8,Bjiij_term8,nocc) + mp2f12_E23(Bijij_term9,Bjiij_term9,nocc)
          print *, ' E23_Bsum:       ', E23_debug
          !print *, ' E23_Bsum_debug: ', mp2f12_E23(Bijij,Bjiij,nocc)
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
          write(*,'(1X,a,g25.16)') ' E22_X_term1: ', X1
          write(*,'(1X,a,g25.16)') ' E22_X_term2: ', X2
          write(*,'(1X,a,g25.16)') ' E22_X_term3: ', X3
          write(*,'(1X,a,g25.16)') ' E22_X_term4: ', X4
          print *, '----------------------------------------'
          E22_debug = X1 + X2 + X3 + X4  
          write(*,'(1X,a,g25.16)') ' E22_Xsum:    ', E22_debug  
          print *, '----------------------------------------'
          print *, ' E_23 B term                            '
          print *, '----------------------------------------'
          write(*,'(1X,a,g25.16)') ' E23_B_term1: ', mp2f12_E23(Bijij_term1,Bjiij_term1,nocc)
          write(*,'(1X,a,g25.16)') ' E23_B_term2: ', mp2f12_E23(Bijij_term2,Bjiij_term2,nocc)
          write(*,'(1X,a,g25.16)') ' E23_B_term3: ', mp2f12_E23(Bijij_term3,Bjiij_term3,nocc)
          write(*,'(1X,a,g25.16)') ' E23_B_term4: ', mp2f12_E23(Bijij_term4,Bjiij_term4,nocc)
          write(*,'(1X,a,g25.16)') ' E23_B_term5: ', mp2f12_E23(Bijij_term5,Bjiij_term5,nocc)
          write(*,'(1X,a,g25.16)') ' E23_B_term6: ', mp2f12_E23(Bijij_term6,Bjiij_term6,nocc)
          write(*,'(1X,a,g25.16)') ' E23_B_term7: ', mp2f12_E23(Bijij_term7,Bjiij_term7,nocc)
          write(*,'(1X,a,g25.16)') ' E23_B_term8: ', mp2f12_E23(Bijij_term8,Bjiij_term8,nocc)
          write(*,'(1X,a,g25.16)') ' E23_B_term9: ', mp2f12_E23(Bijij_term9,Bjiij_term9,nocc)   
          print *, '----------------------------------------'
      
          E23_debug = mp2f12_E23(Bijij_term1,Bjiij_term1,nocc) & 
               & + mp2f12_E23(Bijij_term2,Bjiij_term2,nocc) + mp2f12_E23(Bijij_term3,Bjiij_term3,nocc) &
               & + mp2f12_E23(Bijij_term4,Bjiij_term4,nocc) + mp2f12_E23(Bijij_term5,Bjiij_term5,nocc) &
               & + mp2f12_E23(Bijij_term6,Bjiij_term6,nocc) + mp2f12_E23(Bijij_term7,Bjiij_term7,nocc) &
               & + mp2f12_E23(Bijij_term8,Bjiij_term8,nocc) + mp2f12_E23(Bijij_term9,Bjiij_term9,nocc)
       
          write(*,'(1X,a,g25.16)') ' E23_Bsum:       ', E23_debug
          write(*,'(1X,a,g25.16)') ' E23_Bsum_debug: ', mp2f12_E23(Bijij,Bjiij,nocc)
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
       mp2f12_energy = mp2_energy+E21+E22+ES2

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
       write(*,'(1X,a,f20.10)') 'TOYCODE: F12 ES2 CORRECTION TO ENERGY =  ', ES2
       write(DECinfo%output,*)  'TOYCODE: F12 ES2 CORRECTION TO ENERGY =  ', ES2
       write(*,'(1X,a,f20.10)') 'TOYCODE: F12 CORRECTION TO ENERGY =      ', E21+E22+ES2
       write(DECinfo%output,*)  'TOYCODE: F12 CORRECTION TO ENERGY =      ', E21+E22+ES2       
       ! Total MP2-F12 correlation energy
       ! Getting this energy 
       write(*,'(1X,a)') '----------------------------------------------------'
       write(*,'(1X,a,f20.10)') 'TOYCODE: MP2-F12 CORRELATION ENERGY =    ', mp2f12_energy
       write(DECinfo%output,*)  'TOYCODE: MP2-F12 CORRELATION ENERGY =    ', mp2f12_energy
    endif

    call mem_dealloc(gmo)

  end subroutine full_canonical_mp2_f12
#endif

#ifdef MOD_UNRELEASED
  subroutine submp2f12_EBX(mp2f12_EBX,Bijij,Bjiij,Xijij,Xjiij,Fii,nocc)
    implicit none
    Real(realk)               :: mp2f12_EBX
    Real(realk),intent(INOUT) :: Bijij(nocc,nocc),Bjiij(nocc,nocc)
    Real(realk),intent(IN)    :: Xijij(nocc,nocc),Xjiij(nocc,nocc)
    Real(realk),intent(IN)    :: Fii(nocc,nocc)
    Integer,intent(IN)        :: nocc
    !
    Integer     :: i,j
    Real(realk) :: tmp

    mp2f12_EBX = 0.0E0_realk
    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,j) SHARED(nocc,Fii,Xijij,Bijij)
    DO j=1,nocc
       DO i=1,nocc
          Bijij(i,j) = Bijij(i,j)-(Fii(i,i)+Fii(j,j))*Xijij(i,j)
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    tmp = 0E0_realk
    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,j) SHARED(nocc,Fii,Xijij,Bijij) REDUCTION(+:tmp)
    DO i=1,nocc
       tmp = tmp + Bijij(i,i)
    ENDDO
    !$OMP END PARALLEL DO 
    mp2f12_EBX = mp2f12_EBX + 0.25E0_realk*tmp

    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,j) SHARED(nocc,Fii,Xjiij,Bjiij)
    DO j=1,nocc
       DO i=1,nocc
          Bjiij(i,j) = Bjiij(i,j)-(Fii(i,i)+Fii(j,j))*Xjiij(i,j)
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    tmp = 0E0_realk
    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,j) SHARED(nocc,Bjiij,Bijij) REDUCTION(+:tmp)
    DO j=1,nocc
       DO i=j+1,nocc
          tmp = tmp + 7E0_realk * Bijij(i,j) + Bjiij(i,j)
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO 
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
    Real(realk) :: energy,EK,EJ 
    !
    Integer     :: i,j
    real(realk) :: tmp

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
    EK = 0.0E0_realk
    EJ = 0.0E0_realk

    DO j=1,nocc
       DO i=1,nocc
          EJ = EJ + Vijij(i,j)
          EK = EK + Vjiij(i,j)
       ENDDO
    ENDDO

    !print *, "EJ: ", 5.0/4.0*EJ
    !print *, "EK: ", 1.0/4.0*EK
    !print *, "EJ + EK: ", 5.0/4.0*EJ + 1.0/4.0*EK 

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

    if(MyMolecule%mem_distributed)then
       call lsquit("ERROR(full_get_ccsd_f12_energy): does not work with PDM fullmolecule",-1)
    endif
    ! Init dimensions
    nocc = MyMolecule%nocc
    nvirt = MyMolecule%nvirt
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
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'aiai',gAO,AIBJ)
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
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'apap',gAO,Rapbq)
    !   Ripaq
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'ipap',gAO,Ripaq)
    call mem_dealloc(gao)

    !   Rrrrc
    call mem_alloc(gao,nbasis,nbasis,nbasis,ncabsAO)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRCC')
    !   Rambc
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'amac',gAO,Rambc)
    !   Rimac
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'imac',gAO,Rimac)
    !   Ramic
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'amic',gAO,Ramic)
    call mem_dealloc(gao)

    !   Rrrrc
    call mem_alloc(gao,nbasis,nbasis,nbasis,ncabsAO)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRF')
    !   Fiajb
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'iaia',gAO,Fiajb)
    !   Fijka
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'iiia',gAO,Fijka)
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
    ES2=MyMolecule%EF12singles

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
    integer :: nocc,nvirt,nbasis,print_level,save_model,startidx,endidx,i,j
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
    if(MyMolecule%mem_distributed)then
       call lsquit("ERROR(full_get_ccsd_singles_and_doubles): does not work with PDM fullmolecule",-1)
    endif

    ! Quick fix to always use CCSD model
    !    save_model=DECinfo%ccmodel
    !    DECinfo%ccmodel=3

    if(DECinfo%FrozenCore) then
       nocc = MyMolecule%nval
    else
       nocc = MyMolecule%nocc
    end if
    nvirt = MyMolecule%nvirt
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
             ppfock(i,j) = MyMolecule%oofock%elm2(MyMolecule%ncore+i,MyMolecule%ncore+j)
          end do
       end do

       startidx = MyMolecule%ncore+1  
       endidx = MyMolecule%nocc
       call ccsolver(solver_ccmodel,MyMolecule%Co%elm2(1:nbasis,startidx:endidx),&
            & MyMolecule%Cv%elm2,MyMolecule%fock%elm2, nbasis,nocc,nvirt,mylsitem,&
            & print_level,energy,&
            & VOVO,.false.,local,SOLVE_AMPLITUDES,p2=Tai,p4=Taibj)
       call mem_dealloc(ppfock)

    else

       call ccsolver(solver_ccmodel,MyMolecule%Co%elm2,MyMolecule%Cv%elm2,&
            & MyMolecule%fock%elm2, nbasis,nocc,nvirt,mylsitem, print_level, &
            & energy,VOVO,.false.,local,SOLVE_AMPLITUDES,p2=Tai,p4=Taibj)

    end if

    call tensor_free(VOVO)

    !convert the parallel distributed quantities to local quantities, this is
    !essentially a copying if the tensors are not parallel distributed
    call tensor_minit(Tai_local,[nvirt,nocc],2,atype="LDAR")
    call tensor_add(Tai_local,1.0E0_realk,Tai,a=0.0E0_realk)
    call tensor_minit(Taibj_local,[nvirt,nocc,nvirt,nocc],4,atype="LDAR")
    call tensor_add(Taibj_local,1.0E0_realk,Taibj,a=0.0E0_realk)

    call tensor_free(Taibj)
    call tensor_free(Tai)


  end subroutine full_get_ccsd_singles_and_doubles

end module full
