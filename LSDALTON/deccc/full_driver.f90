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
  use full_f12contractions
  use array4_simple_operations
  use array3_simple_operations
  use array2_simple_operations
  use mp2_module
!  use orbital_operations
  use full_molecule
  use ccintegrals!,only: get_full_AO_integrals,get_AO_hJ,get_AO_K,get_AO_Fock
  use ccdriver!,only: ccsolver_justenergy, ccsolver

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

    ! run cc program
    if(DECinfo%F12) then ! F12 correction
       if(DECinfo%ccModel==1) then
          call full_canonical_mp2_f12(MyMolecule,MyLsitem,D,Ecorr)
       else
          call full_get_ccsd_f12_energy(MyMolecule,MyLsitem,D,Ecorr)
       end if
    else
       if(DECinfo%ccModel==1) then
          call full_canonical_mp2_correlation_energy(MyMolecule,mylsitem,Ecorr)
       else
          call full_cc_dispatch(MyMolecule,mylsitem,Ecorr)          
       end if
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
    real(realk),pointer :: ppfock_fc(:,:), ypo_fc(:,:)

    Ecorr = 0.0E0_realk

    if(DECinfo%FrozenCore) then
       nocc = MyMolecule%nval
    else
       nocc = MyMolecule%numocc
    end if
    nunocc = MyMolecule%numvirt
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
       call mem_alloc(ypo_fc,nbasis,nocc)
       do j=1,nocc
          do i=1,nbasis
             ypo_fc(i,j) = MyMolecule%ypo(i,MyMolecule%ncore+j)
          end do
       end do       

       if (DECinfo%ccModel .eq. 4) then
          ! ccsd(t) correction
          Ecorr = ccsolver_justenergy_pt(MyMolecule,nbasis,nocc,nunocc,&
               & mylsitem,print_level,fragment_job,ypo_fc=ypo_fc,ppfock_fc=ppfock_fc)
       else
          Ecorr = ccsolver_justenergy(ypo_fc,&
               & MyMolecule%ypv,MyMolecule%fock, nbasis,nocc,nunocc,mylsitem,&
               & print_level,fragment_job,ppfock_fc,MyMolecule%qqfock)
       end if

       call mem_dealloc(ppfock_fc)
       call mem_dealloc(ypo_fc)

    else

       if (Decinfo%ccModel .eq. 4) then

          Ecorr = ccsolver_justenergy_pt(MyMolecule,nbasis,nocc,nunocc,&
               & mylsitem,print_level,fragment_job)

       else

          Ecorr = ccsolver_justenergy(MyMolecule%ypo,MyMolecule%ypv,&
               & MyMolecule%fock, nbasis,nocc,nunocc,mylsitem, &
               & print_level,fragment_job,MyMolecule%ppfock,MyMolecule%qqfock)

       end if

    end if

  end subroutine full_cc_dispatch

  !> \brief Calculate canonical MP2 energy for full molecular system
  !> keeping full AO integrals in memory. Only for testing.
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine full_canonical_mp2_f12(MyMolecule,MyLsitem,Dmat,energy)


    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> HF density matrix
    type(matrix),intent(in) :: Dmat
    !> Canonical MP2 correlation energy
    real(realk),intent(inout) :: energy
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

    real(realk),pointer :: Vjiij(:,:)    
    real(realk),pointer :: Vjiij_term1(:,:)
    real(realk),pointer :: Vjiij_term2(:,:)
    real(realk),pointer :: Vjiij_term3(:,:)
    real(realk),pointer :: Vjiij_term4(:,:)

    real(realk),pointer :: Xijkl(:,:,:,:)
    real(realk),pointer :: Xijij(:,:)
    real(realk),pointer :: Xjiij(:,:)
    real(realk),pointer :: Bijij(:,:)
    real(realk),pointer :: Bjiij(:,:)
    integer :: nbasis,ncabs,nocc,nvirt,I,A,B,J,noccfull,ncabsAO
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
    Real(realk)  :: E21,E22,Gtmp
    type(array2) :: array2Tai
    type(array4) :: array4Taibj

    ! Init stuff
    ! **********
    call init_cabs
    nbasis = MyMolecule%nbasis
    nocc   = MyMolecule%numocc
    nvirt  = MyMolecule%numvirt
    call determine_CABS_nbast(ncabsAO,ncabs,mylsitem%setting,DECinfo%output)
    noccfull = nocc

    ! Memory check
    call full_canonical_mp2_memory_check(nbasis,nocc,nvirt)

    call get_F12_mixed_MO_Matrices(MyLsitem,MyMolecule,Dmat,nbasis,ncabsAO,&
         & nocc,noccfull,nvirt,ncabs,HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp)

    ! Get all AO integrals
    ! ********************

    call mem_alloc(gao,nbasis,nbasis,nbasis,nbasis)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRC')
    ! Transform AO integrals to MO integrals (A I | B J)
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%ypo, MyMolecule%ypv,'aiai',gAO,gMO)
    call mem_dealloc(gao)

    call get_4Center_F12_integrals(mylsitem,MyMolecule,nbasis,nocc,noccfull,nvirt,ncabsAO,&
         & Ripjq,Fijkl,Tijkl,Rimjc,Dijkl,Tirjk,Tijkr,Gipjq,Gimjc,Girjs,Girjm,&
         & Grimj,Gipja,Gpiaj,Gicjm,Gcimj,Gcirj,Gciaj,Giajc)

    call mem_alloc(Vijij,nocc,nocc)
    call mem_alloc(Vijij_term1,nocc,nocc)
    call mem_alloc(Vijij_term2,nocc,nocc)
    call mem_alloc(Vijij_term3,nocc,nocc)
    call mem_alloc(Vijij_term4,nocc,nocc)

    call mem_alloc(Vjiij,nocc,nocc)
    call mem_alloc(Vjiij_term1,nocc,nocc)
    call mem_alloc(Vjiij_term2,nocc,nocc)
    call mem_alloc(Vjiij_term3,nocc,nocc)
    call mem_alloc(Vjiij_term4,nocc,nocc)

    call mp2f12_Vijij(Vijij,Ripjq,Gipjq,Fijkl,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)
    call mp2f12_Vjiij(Vjiij,Ripjq,Gipjq,Fijkl,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)

    call mp2f12_Vijij_term1(Vijij_term1,Fijkl,nocc,noccfull,nbasis,ncabs)
    call mp2f12_Vijij_term2(Vijij_term2,Ripjq,Gipjq,nocc,noccfull,nbasis,ncabs)
    call mp2f12_Vijij_term3(Vijij_term3,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)
    call mp2f12_Vijij_term4(Vijij_term4,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)

    call mp2f12_Vjiij_term1(Vjiij_term1,Fijkl,nocc,noccfull,nbasis,ncabs)
    call mp2f12_Vjiij_term2(Vjiij_term2,Ripjq,Gipjq,nocc,noccfull,nbasis,ncabs)
    call mp2f12_Vjiij_term3(Vjiij_term3,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)
    call mp2f12_Vjiij_term4(Vjiij_term4,Rimjc,Gimjc,nocc,noccfull,nbasis,ncabs)

    print *, 'E_21_V_term1: ', 2.0E0_REALK*mp2f12_EV(Vijij_term1,Vjiij_term1,nocc)
    print *, 'E_21_V_term2: ', 2.0E0_REALK*mp2f12_EV(Vijij_term2,Vjiij_term2,nocc)
    print *, 'E_21_V_term3: ', 2.0E0_REALK*mp2f12_EV(Vijij_term3,Vjiij_term3,nocc)
    print *, 'E_21_V_term4: ', 2.0E0_REALK*mp2f12_EV(Vijij_term4,Vjiij_term4,nocc)
    print *, '----------------------------------------'
    print *, 'E_21_Vsum: ',  2.0E0_REALK*(mp2f12_EV(Vijij_term1,Vjiij_term1,nocc) + mp2f12_EV(Vijij_term2,Vjiij_term2,nocc) &
         & + mp2f12_EV(Vijij_term3,Vjiij_term3,nocc) + mp2f12_EV(Vijij_term4,Vjiij_term4,nocc) ) 

    print *, 'E_21_V: ',  2.0E0_REALK*mp2f12_EV(Vijij,Vjiij,nocc)

    call mem_alloc(Ciajb,nocc,nvirt,nocc,nvirt)
 !   call mem_alloc(Cjaib,nocc,nvirt,nocc,nvirt)
    call mp2f12_Ciajb(Ciajb,Giajc,Fac%elms,nocc,nvirt,ncabs)
!    call mp2f12_Cjaib(Cjaib,Giajc,Fac%elms,nocc,nvirt,ncabs)


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
       energy = 0.0E0_realk
       do J=1,nocc
          do B=1,nvirt
             do I=1,nocc
                do A=1,nvirt
                   
                   ! Difference in orbital energies: eps(I) + eps(J) - eps(A) - eps(B)
                   eps = MyMolecule%ppfock(I,I) + MyMolecule%ppfock(J,J) &
                        & - MyMolecule%qqfock(A,A) - MyMolecule%qqfock(B,B)
                   
                   ! Energy = sum_{AIBJ} (AI|BJ) * [ 2(AI|BJ) - (BI|AJ) ] / (epsI + epsJ - epsA - epsB)
                   energy = energy + gmo(A,I,B,J)*(2E0_realk*gmo(A,I,B,J)-gmo(B,I,A,J))/eps
                   
                end do
             end do
          end do
       end do
       write(DECinfo%output,*) 'TOYCODE: MP2 CORRELATION ENERGY = ', energy
       print *, 'TOYCODE: MP2 CORRELATION ENERGY = ', energy

    else
       !  THIS PIECE OF CODE IS MORE GENERAL AS IT DOES NOT REQUIRE CANONICAL ORBITALS
       !    ! Get full MP2 (as specified in input)
       call full_get_ccsd_singles_and_doubles(MyMolecule,MyLsitem,array2Tai, array4Taibj)
       !Calculate standard MP2 energy (both canonical and noncanonical)

       call mem_alloc(Taibj,nvirt,nocc,nvirt,nocc)

       energy=0.0E0_realk
       do j=1,nocc
          do b=1,nvirt
             do i=1,nocc
                do a=1,nvirt
                   ! Energy = sum_{ijab} ( Taibj) * (ai | bj)
                   Taibj(a,i,b,j) = array4Taibj%val(a,i,b,j)
                   Gtmp = 2.0E0_realk * gmo(a,i,b,j) - gmo(b,i,a,j)
                   energy = energy + Taibj(a,i,b,j) * Gtmp
                end do
             end do
          end do
       end do
       write(DECinfo%output,*) 'TOYCODE: MP2 CORRELATION ENERGY = ', energy
       print *, 'TOYCODE: MP2 CORRELATION ENERGY = ', energy
    end if

    call mp2f12_Vijij_coupling(Vijij,Ciajb,Taibj,nocc,nvirt)
    call mp2f12_Vjiij_coupling(Vjiij,Ciajb,Taibj,nocc,nvirt)

    E21 = 2.0E0_REALK*mp2f12_EV(Vijij,Vjiij,nocc)
!   write(*,*) 'MP2f12 energy term 2 <0|H1|1>',  E21

    call mem_dealloc(Vijij)
    call mem_dealloc(Vijij_term1)
    call mem_dealloc(Vijij_term2)
    call mem_dealloc(Vijij_term3)
    call mem_dealloc(Vijij_term4)

    call mem_dealloc(Vjiij)
    call mem_dealloc(Vjiij_term1)
    call mem_dealloc(Vjiij_term2)
    call mem_dealloc(Vjiij_term3)
    call mem_dealloc(Vjiij_term4)

    call mem_dealloc(Taibj)
    call mem_dealloc(Ciajb)

    if(DECinfo%use_canonical) then
       call mem_alloc(Xijij,nocc,nocc)
       call mem_alloc(Xjiij,nocc,nocc)
       call mp2f12_Xijij(Xijij,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
       call mp2f12_Xjiij(Xjiij,Gipjq,Tijkl,Gimjc,nocc,noccfull,nbasis,ncabs)
    else
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
!   write(*,*) 'MP2f12 energy term <1|H0-E0|1>',E22
    call free_F12_mixed_MO_Matrices(HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp)

    if(DECinfo%use_canonical) then
       call mem_dealloc(Xijij)
       call mem_dealloc(Xjiij)
    else
       call mem_dealloc(Xijkl)
    endif

    call mem_dealloc(Bijij)
    call mem_dealloc(Bjiij)

    call free_4Center_F12_integrals(&
         & Ripjq,Fijkl,Tijkl,Rimjc,Dijkl,Tirjk,Tijkr,Gipjq,Gimjc,Girjs,Girjm,&
         & Grimj,Gipja,Gpiaj,Gicjm,Gcimj,Gcirj,Gciaj,Giajc)
    call free_cabs

!    write(DECinfo%output,*) 'TOYCODE: MP2 CORRELATION ENERGY = ', energy
!    print *, 'TOYCODE: MP2 CORRELATION ENERGY = ', energy

    write(*,*) 'TOYCODE: F12 E21 CORRECTION TO ENERGY = ',E21
    write(DECinfo%output,*) 'TOYCODE: F12 E21 CORRECTION TO ENERGY = ',E21
    write(*,*) 'TOYCODE: F12 E22 CORRECTION TO ENERGY = ',E22
    write(DECinfo%output,*) 'TOYCODE: F12 E22 CORRECTION TO ENERGY = ',E22

    write(*,*) 'TOYCODE: F12 CORRECTION TO ENERGY = ',E21+E22
    write(DECinfo%output,*) 'TOYCODE: F12 CORRECTION TO ENERGY = ',E21+E22

    ! Total MP2-F12 correlation energy
    ! Getting this energy 

    energy = energy + E21 + E22
    print *, 'TOYCODE: MP2-F12 CORRELATION ENERGY = ', energy
    write(DECinfo%output,*) 'TOYCODE: MP2-F12 CORRELATION ENERGY = ', energy

    call array4_free(array4Taibj)
    call array2_free(array2Tai)
    call mem_dealloc(gmo)

  end subroutine full_canonical_mp2_f12


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

  function mp2f12_EV(Vijij,Vjiij,nocc)
  implicit none
  Real(realk) :: mp2f12_EV
  Real(realk),intent(IN) :: Vijij(nocc,nocc),Vjiij(nocc,nocc)
  Integer,intent(IN)     :: nocc
  !
  Integer     :: i,j
  Real(realk) :: tmp

  tmp = 0E0_realk
  DO i=1,nocc
    tmp = tmp + Vijij(i,i)
  ENDDO
  mp2f12_EV = -0.5E0_realk*tmp

  tmp = 0E0_realk
  DO j=1,nocc
    DO i=j+1,nocc
      tmp = tmp + 5E0_realk * Vijij(i,j) - Vjiij(i,j)
    ENDDO
  ENDDO
  mp2f12_EV = mp2f12_EV - 0.25E0_realk*tmp
  end function mp2f12_EV

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
          call dcopy(nvirt*nbasis,Cvirt,1,CMO(i)%elms(noccfull*nbasis+1:nbasis*nbasis),1)
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

  subroutine MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
       & Cocc,Cvirt,inputstring,matAO,matMO)
    implicit none
    !> Lsitem structure
    integer :: nocc,noccfull,nvirt,nCabsAO,nCabs,nbasis
    type(lsitem), intent(inout) :: mylsitem
    integer :: ndim2(2),ndim1(2)
    type(matrix) :: matAO,matMO
    real(realk),pointer :: elms(:)
    type(matrix) :: CMO(2)
    real(realk),dimension(nbasis,nocc),intent(in) :: Cocc
    !> Virtual MO coefficients
    real(realk),dimension(nbasis,nvirt),intent(in) :: Cvirt
    type(matrix) :: CMO_cabs,CMO_ri,tmp
    character(len=2) :: inputstring
    logical :: doCABS,doRI
    integer :: i,lupri
    character :: string(2)
    string(1)=inputstring(1:1) 
    string(2)=inputstring(2:2) 
    lupri=6
    doCABS = .FALSE.
    do i=1,2
       if(string(i).EQ.'c')then !cabs
          doCABS = .TRUE.
       endif
    enddo
    doRI = .FALSE.
    do i=1,2
       if(string(i).EQ.'r')then !RI
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
    do i=1,2
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
          call dcopy(nvirt*nbasis,Cvirt,1,CMO(i)%elms(noccfull*nbasis+1:nbasis*nbasis),1)
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
    call mat_init(tmp,CMO(1)%ncol,matAO%ncol)
    call mat_mul(CMO(1),matAO,'t','n',1E0_realk,0E0_realk,tmp)
    call mat_mul(tmp,CMO(2),'n','n',1E0_realk,0E0_realk,matMO)
    call mat_free(tmp)
    do i=1,2
       call mat_free(CMO(i))
    enddo

  end subroutine MO_transform_AOMatrix

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

    ! Init dimensions
    call init_cabs
    nocc = MyMolecule%numocc
    nvirt = MyMolecule%numvirt
    nbasis = MyMolecule%nbasis
    noccfull = nocc
    call determine_CABS_nbast(ncabsAO,ncabs,mylsitem%setting,DECinfo%output)

    ! Get full CCSD singles (Tai) and doubles (Taibj) amplitudes
    call full_get_ccsd_singles_and_doubles(MyMolecule,MyLsitem,Tai, Taibj)

    ! Get all MO mixed matrices
    call get_F12_mixed_MO_Matrices(MyLsitem,MyMolecule,Dmat,nbasis,ncabsAO,&
         & nocc,noccfull,nvirt,ncabs,HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp)

    ! Get all AO integrals in regular basis
    call mem_alloc(gao,nbasis,nbasis,nbasis,nbasis)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRC')

    ! Transform AO integrals to MO integrals (A I | B J)
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%ypo, MyMolecule%ypv,'aiai',gAO,AIBJ)
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
         & MyMolecule%ypo, MyMolecule%ypv,'apap',gAO,Rapbq)
!   Ripaq
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%ypo, MyMolecule%ypv,'ipap',gAO,Ripaq)
    call mem_dealloc(gao)

!   Rrrrc
    call mem_alloc(gao,nbasis,nbasis,nbasis,ncabsAO)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRCC')
!   Rambc
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%ypo, MyMolecule%ypv,'amac',gAO,Rambc)
!   Rimac
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%ypo, MyMolecule%ypv,'imac',gAO,Rimac)
!   Ramic
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%ypo, MyMolecule%ypv,'amic',gAO,Ramic)
    call mem_dealloc(gao)

!   Rrrrc
    call mem_alloc(gao,nbasis,nbasis,nbasis,ncabsAO)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRF')
!   Fiajb
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%ypo, MyMolecule%ypv,'iaia',gAO,Fiajb)
!   Fijka
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         & MyMolecule%ypo, MyMolecule%ypv,'iiia',gAO,Fijka)
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

    E21 = 2.0E0_realk*mp2f12_EV(Vijij,Vjiij,nocc)
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

    EF12=E21 + E22

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
    call free_cabs

  end subroutine full_get_ccsd_f12_energy


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


    ! Quick fix to always use CCSD model
!    save_model=DECinfo%ccmodel
!    DECinfo%ccmodel=3

    if(DECinfo%FrozenCore) then
       nocc = MyMolecule%nval
    else
       nocc = MyMolecule%numocc
    end if
    nunocc = MyMolecule%numvirt
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
       endidx = MyMolecule%numocc
       call ccsolver(MyMolecule%ypo(1:nbasis,startidx:endidx),&
            & MyMolecule%ypv,MyMolecule%fock, nbasis,nocc,nunocc,mylsitem,&
            & print_level,fragment_job,&
            & ppfock,MyMolecule%qqfock,energy, Tai, Taibj, VOVO,.false.)
       call mem_dealloc(ppfock)

    else

       call ccsolver(MyMolecule%ypo,MyMolecule%ypv,&
            & MyMolecule%fock, nbasis,nocc,nunocc,mylsitem, print_level, fragment_job,&
            & MyMolecule%ppfock,MyMolecule%qqfock,&
            & energy, Tai, Taibj, VOVO,.false.)
    end if

    call array4_free(VOVO)

  end subroutine full_get_ccsd_singles_and_doubles


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
       &                          MyMolecule%ypo, MyMolecule%ypv,'ipip',gAO,Ripjq)

    !Calculate the various Gaussian geminal integrals with four regular AO indeces
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
       &                          MyMolecule%ypo, MyMolecule%ypv,'ipip',gAO,Gipjq)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRF')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
       &                          MyMolecule%ypo, MyMolecule%ypv,'iiii',gAO,Fijkl)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRD')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &                        MyMolecule%ypo, MyMolecule%ypv,'iiii',gAO,Dijkl)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRR2')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
       &                          MyMolecule%ypo, MyMolecule%ypv,'iiii',gAO,Tijkl)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%ypo, MyMolecule%ypv,'ipia',gAO,Gipja)


    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%ypo, MyMolecule%ypv,'piai',gAO,Gpiaj)


    call mem_dealloc(gao)

    !Calculate the various Gaussian geminal integrals with RRRC
    call mem_alloc(gao,nbasis,nbasis,nbasis,ncabsAO)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRC2')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
       &                          MyMolecule%ypo, MyMolecule%ypv,'iiir',gAO,Tijkr)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRCC')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
       &                          MyMolecule%ypo, MyMolecule%ypv,'imic',gAO,Rimjc)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRCG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
       &                          MyMolecule%ypo, MyMolecule%ypv,'imic',gAO,Gimjc)
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
       &                          MyMolecule%ypo, MyMolecule%ypv,'iaic',gAO,Giajc)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRC2')

    call mem_dealloc(gao)
    !Calculate the various Gaussian geminal integrals with RCRR
    call mem_alloc(gao,nbasis,ncabsAO,nbasis,nbasis)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RCRR2')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
       &                          MyMolecule%ypo, MyMolecule%ypv,'irii',gAO,Tirjk)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RCRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%ypo, MyMolecule%ypv,'irim',gAO,Girjm)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RCRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%ypo, MyMolecule%ypv,'icim',gAO,Gicjm)

    call mem_dealloc(gao)
    !Calculate the various Gaussian geminal integrals with RCRC
    call mem_alloc(gao,nbasis,ncabsAO,nbasis,ncabsAO)

    !Calculate the various Gaussian geminal integrals with four regular AO indeces
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RCRCG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
       &                          MyMolecule%ypo, MyMolecule%ypv,'irir',gAO,Girjs)

    call mem_dealloc(gao)
    !Calculate the various Gaussian geminal integrals with CRRR
    call mem_alloc(gao,ncabsAO,nbasis,nbasis,nbasis)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'CRRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%ypo, MyMolecule%ypv,'rimi',gAO,Grimj)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'CRRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%ypo, MyMolecule%ypv,'cimi',gAO,Gcimj)

    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'CRRRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%ypo, MyMolecule%ypv,'ciai',gAO,Gciaj)

    call mem_dealloc(gao)
    !Calculate the various Gaussian geminal integrals with CRCR
    call mem_alloc(gao,ncabsAO,nbasis,ncabsAO,nbasis)
    gao = 0.0E0_realk
    call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'CRCRG')
    call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
         &  MyMolecule%ypo, MyMolecule%ypv,'ciri',gAO,Gcirj)
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

    subroutine get_F12_mixed_MO_Matrices(MyLsitem,MyMolecule,Dmat,nbasis,ncabsAO,&
         & nocc,noccfull,nvirt,ncabs,HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp)

      implicit none
      !> Full molecule info
      type(fullmolecule), intent(in) :: MyMolecule
      !> Lsitem structure
      type(lsitem), intent(inout) :: mylsitem
      integer :: nbasis,nocc,nvirt,noccfull,ncabsAO,ncabs
      type(matrix) :: Dmat,K
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
! Temp
      type(matrix) :: HJrc
      type(matrix) :: Kcc
      type(matrix) :: Fcc

    ! Mixed regular/CABS one-electron and Coulomb matrix (h+J) combination in AO basis
    call mat_init(HJrc,nbasis,ncabsAO)
    call get_AO_hJ(nbasis,ncabsAO,HJrc,Dmat,MyLsitem,'RCRRC')
    call mat_init(HJir,nocc,ncabsAO)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
       & MyMolecule%ypo, MyMolecule%ypv,'ir',HJrc,HJir)
    call mat_free(HJrc)

    ! Mixed CABS/CABS exchange matrix
    call mat_init(Kcc,ncabsAO,ncabsAO)
    call get_AO_K(nbasis,ncabsAO,Kcc,Dmat,MyLsitem,'CCRRC')
    call mat_init(Krr,ncabsAO,ncabsAO)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
       & MyMolecule%ypo, MyMolecule%ypv,'rr',Kcc,Krr)
    call mat_free(Kcc)

    ! Mixed CABS/CABS Fock matrix
    call mat_init(Fcc,ncabsAO,ncabsAO)
    call get_AO_Fock(nbasis,ncabsAO,Fcc,Dmat,MyLsitem,'CCRRC')
    call mat_init(Frr,ncabsAO,ncabsAO)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
       & MyMolecule%ypo, MyMolecule%ypv,'rr',Fcc,Frr)
    call mat_free(Fcc)

    ! Mixed AO/CABS Fock matrix
    call mat_init(Frc,nbasis,ncabsAO)
    call get_AO_Fock(nbasis,ncabsAO,Frc,Dmat,MyLsitem,'RCRRC')
    call mat_init(Fac,nvirt,ncabs)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
       & MyMolecule%ypo, MyMolecule%ypv,'ac',Frc,Fac)
    call mat_free(Frc)

    ! Mixed AO/AO full MO Fock matrix
    call mat_init(Fcc,nbasis,nbasis)
    call get_AO_Fock(nbasis,ncabsAO,Fcc,Dmat,MyLsitem,'RRRRC')
    !Fpp
    call mat_init(Fpp,nbasis,nbasis)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
       & MyMolecule%ypo, MyMolecule%ypv,'pp',Fcc,Fpp)
    !Fii
    call mat_init(Fii,nocc,nocc)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
       & MyMolecule%ypo, MyMolecule%ypv,'ii',Fcc,Fii)
    !Fmm
    call mat_init(Fmm,noccfull,noccfull)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
       & MyMolecule%ypo, MyMolecule%ypv,'mm',Fcc,Fmm)
    call mat_free(Fcc)

    ! Mixed CABS/AO MO Fock matrix
    call mat_init(Fcc,ncabsAO,nbasis)
    call get_AO_Fock(nbasis,ncabsAO,Fcc,Dmat,MyLsitem,'CRRRC')
    !Frm
    call mat_init(Frm,ncabsAO,noccfull)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
       & MyMolecule%ypo, MyMolecule%ypv,'rm',Fcc,Frm)
    !Fcc
    call mat_init(Fcp,ncabs,nbasis)
    call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
       & MyMolecule%ypo, MyMolecule%ypv,'cp',Fcc,Fcp)
    call mat_free(Fcc)
  end subroutine get_F12_mixed_MO_Matrices

    subroutine free_F12_mixed_MO_Matrices(HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp)

      implicit none
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
      call mat_free(HJir)
      call mat_free(Krr)
      call mat_free(Frr)
      call mat_free(Fac)
      call mat_free(Fpp)
      call mat_free(Fii)
      call mat_free(Fmm)
      call mat_free(Frm)
      call mat_free(Fcp)

  end subroutine free_F12_mixed_MO_Matrices




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
    type(array2) :: Cocc, Cunocc
    type(array4) :: g
    integer :: nbasis,i,j,a,b,ncore,offset,nocc,nunocc
    integer, dimension(2) :: occ_dims,unocc_dims
    real(realk) :: eps
    real(realk), pointer :: ppfock(:,:)

    ! Sanity check
    if(.not. DECinfo%use_canonical) then
       call lsquit('full_canonical_mp2_correlation_energy requires canonical orbitals!',-1)
    end if


    ! Initialize stuff
    ! ****************

    if(DECinfo%frozencore) then
       ! Frozen core: Only valence orbitals
       nocc = MyMolecule%nval
    else
       nocc = MyMolecule%numocc
    end if

    nunocc = MyMolecule%numvirt
    ncore = MyMolecule%ncore
    nbasis=MyMolecule%nbasis
    occ_dims = [nbasis,nocc]
    unocc_dims = [nbasis,nunocc]
    call mem_alloc(ppfock,nocc,nocc)
    if(DECinfo%frozencore) then
       ! Only copy valence orbitals into array2 structure
       Cocc=array2_init(occ_dims)
       do i=1,nocc
          Cocc%val(:,i) = MyMolecule%ypo(:,i+Ncore)
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
       Cocc=array2_init(occ_dims,MyMolecule%ypo)
       ppfock = MyMolecule%ppfock
       offset=0
    end if
    Cunocc=array2_init(unocc_dims,MyMolecule%ypv)



    ! Get (AI|BJ) integrals stored in the order (A,I,B,J)
    ! ***************************************************
    call get_VOVO_integrals(mylsitem,nbasis,nocc,nunocc,Cunocc,Cocc,g)
    call array2_free(Cocc)
    call array2_free(Cunocc)



    ! Calculate canonical MP2 energy
    ! ******************************
    Ecorr = 0.0_realk
    do J=1,nocc
       do B=1,nunocc
          do I=1,nocc
             do A=1,nunocc
                ! Difference in orbital energies: eps(I) + eps(J) - eps(A) - eps(B)
                eps = MyMolecule%ppfock(I,I) + MyMolecule%ppfock(J,J) &
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