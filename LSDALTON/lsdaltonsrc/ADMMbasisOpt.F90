module ADMMbasisOptMod
use files
  use precision
  use matrix_module!, only: matrix
  use TYPEDEFTYPE!, only: lsitem
  use configurationType!, only: configitem
  use memory_handling!, only: mem_alloc, mem_dealloc
  use Energy_and_deriv
  use dal_interface
!===========================================================!
!    The main driver for ADMM basis set optimization in LSDALTON  !
!===========================================================!
! Written by Thomas Kjaergaard 2014
!
  private
  public :: ADMMbasisOptSub
CONTAINS
  SUBROUTINE ADMMbasisOptSub(E,config,H1,F,D,S,CMO,ls)
    Implicit none
    !  All these general entities needed to get energy and gradient
    Type(lsitem)                    :: ls 
    Type(Matrix), intent(inout)     :: F(1),D(1),S ! Fock,density,overlap matrices
    Type(Matrix), intent(inout)     :: H1   ! One electron matrix
    Type(Matrix), intent(inout)     :: CMO       ! Orbitals
    Type(ConfigItem), intent(inout) :: Config ! General information
    Real(realk),intent(inout)       :: E(1)   ! Energy
    ! local variables
    Real(realk):: EcontADMM(5),Econt(5),Etotal(1),Etmp
    integer :: lupri,luerr,ndmat,NAtoms,LUADMM
    logical :: ADMM_EXCHANGE,PRINT_EK3,ADMMBASISFILE,ADMM_separateX,CFG_restart
    lupri = config%lupri
    luerr = config%lupri
    ndmat = 1
    !At this point the code have converged the SCF using the ADMM approximation of choice
    !we therefore calculate the energy 
    call di_get_fock_LSDALTON(D,h1,F,ndmat,Etotal,lupri,luerr,ls,EcontADMM)
 
    Econt(2) = EcontADMM(2) !X(D_admm)
    Econt(3) = EcontADMM(3) !k(d_admm)
    Econt(4) = EcontADMM(4) !x(d_admm)

    !step 2 deactivate the ADMM and recalculate SCF energy     
    ADMM_EXCHANGE = ls%setting%scheme%ADMM_EXCHANGE
    ls%setting%scheme%ADMM_EXCHANGE = .FALSE.
    ls%input%dalton%ADMM_EXCHANGE = .FALSE.
    
    NAtoms = config%Molecule%nAtoms
    !make sure to restart to save time
    CFG_restart = config%diag%CFG_restart
    config%diag%CFG_restart = .TRUE.
    call Get_Energy(Etotal,Etmp,config,H1,F,D,S,ls,CMO,Natoms,lupri,luerr)
    config%diag%CFG_restart = CFG_restart

    !reactivate the ADMM and extract the K(D) energy
    ls%setting%scheme%ADMM_EXCHANGE = ADMM_EXCHANGE
    ls%input%dalton%ADMM_EXCHANGE = ADMM_EXCHANGE

    PRINT_EK3 = ls%setting%scheme%PRINT_EK3
    ls%setting%scheme%PRINT_EK3 = .TRUE. 
    ADMMBASISFILE = ls%input%dalton%ADMMBASISFILE
    ls%input%dalton%ADMMBASISFILE = .TRUE.
    call di_get_fock_LSDALTON(D,h1,F,ndmat,Etotal,lupri,luerr,ls,EcontADMM)
    Econt(1) = EcontADMM(1) !K(D)
    !K(D) - k(d_admm) - X(D_admm) + x(d_admm)
    Econt(5) = Econt(1) - Econt(3) - Econt(2) + Econt(4) 


    write(lupri,*) "FINAL ADMMmin: K(D) energy = ",Econt(1)
    write(lupri,*) "FINAL ADMMmin: X(D) energy = ",Econt(2)
    write(lupri,*) "FINAL ADMMmin: k(d) energy = ",Econt(3)
    write(lupri,*) "FINAL ADMMmin: x(d) energy = ",Econt(4)
    Econt(5) = Econt(1) - Econt(3) - Econt(2) + Econt(4)
    write(lupri,*) "===================================================================="
    write(lupri,*) "FINAL ADMMmin: K(D) - k(d) - X(D) + x(d) = ",Econt(5)
    write(lupri,*) "===================================================================="

    write(lupri,*) "ADMMminDATA"


!    CALL LSOPEN(LUADMM,'ADMMmin.dat','UNKNOWN','FORMATTED')    
!    WRITE(LUADMM,'(5F20.13)') Econt(5), Econt(1), Econt(2), Econt(3), Econt(4)                  
    WRITE(LUPRI,'(5F20.13)') Econt(5), Econt(1), Econt(2), Econt(3), Econt(4)                  
!    CALL LSCLOSE(LUADMM,'KEEP')

    ls%input%dalton%ADMMBASISFILE = ADMMBASISFILE
    ls%setting%scheme%PRINT_EK3 = PRINT_EK3
  end SUBROUTINE ADMMbasisOptSub

end module ADMMbasisOptMod




