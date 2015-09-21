module InteractionEnergyMod
use files
  use precision
  use matrix_module!, only: matrix
  use TYPEDEFTYPE!, only: lsitem
  use configurationType!, only: configitem
  use memory_handling!, only: mem_alloc, mem_dealloc
  use molecule_module!, only: print_geometry
  use lstiming!, only: lstimer
  use WRITEMOLEFILE!, only: write_molecule_output
  use Energy_and_deriv
  use dec_typedef_module,only: DECinfo
!===========================================================!
!    The main driver for geometry optimization in LSDALTON  !
!===========================================================!
! Written by Thomas Kjaergaard 2014
!
  private
  public :: InteractionEnergy
CONTAINS
  SUBROUTINE InteractionEnergy(E,config,H1,F,D,S,CMO,ls)
    Implicit none
    !  All these general entities needed to get energy and gradient
    Type(lsitem)                    :: ls 
    Type(Matrix), intent(inout)     :: F(1),D(1),S ! Fock,density,overlap matrices
    Type(Matrix), intent(inout)     :: H1   ! One electron matrix
    Type(Matrix), intent(inout)     :: CMO       ! Orbitals
    Type(ConfigItem), intent(inout) :: Config ! General information
    Real(realk),intent(inout)       :: E(1)   ! Energy
    ! local variables
    Real(realk):: E1(1),E2(1)

    CALL SubSystemEnergy(E1,config,H1,F,D,S,CMO,ls,1)
    IF(config%SameSubSystems)THEN
       E2 = E1
    ELSE
       CALL SubSystemEnergy(E2,config,H1,F,D,S,CMO,ls,2)
    ENDIF
    WRITE(config%lupri,*)''
    WRITE(config%lupri,*)' Summary For Counter Poise Correction'
    WRITE(config%lupri,*)''
    WRITE(config%lupri,'(A23,A,A2,F20.12)')'  Energy for subsystem ',TRIM(ls%input%Molecule%SubSystemLabel(1)),': ',E1(1)
    WRITE(config%lupri,'(A23,A,A2,F20.12)')'  Energy for subsystem ',TRIM(ls%input%Molecule%SubSystemLabel(2)),': ',E2(1)
    WRITE(config%lupri,'(A26,F20.12)')     '  Energy for full System: ',E(1)
    WRITE(config%lupri,'(A)') '---------------------------------------------------------'
    WRITE(config%lupri,'(A26,F20.12)') '  Interaction Energy    : ',E(1) - E1(1) - E2(1)
    WRITE(config%lupri,'(A)') '========================================================='
    WRITE(config%lupri,*)''

  end SUBROUTINE InteractionEnergy

  subroutine SubSystemEnergy(Esub,config,H1,F,D,S,CMO,ls,subIndex)
    Implicit none
    !  All these general entities needed to get energy and gradient
    Type(lsitem)                    :: ls 
    Type(Matrix), intent(inout)     :: F(1),D(1),S ! Fock,density,overlap matrices
    Type(Matrix), intent(inout)     :: H1   ! One electron matrix
    Type(Matrix), intent(inout)     :: CMO       ! Orbitals
    Type(ConfigItem), intent(inout) :: Config ! General information
    Real(realk),intent(inout)       :: Esub(1)   ! Energy
    integer,intent(in) :: subIndex
    ! local variables
    integer :: lupri,luerr,Natoms,nbast,nbast1,nbast2
    real(realk) :: Etmp,nelectrons,nelectronsSub
    logical,pointer :: Phantom(:)
    integer :: I,J,R,jcharge,jtype,norbJ,nOrbI,itype,icharge,iOrb,jOrb
    integer :: restart_lun,dens_lun
    logical :: CFG_restart,CFG_purifyrestart,dens_exsist,OnMaster,gcbasis,OnlyOccPart
    type(matrix) :: D1(1),TMP
    real(realk),pointer :: Dfull(:,:),D1full(:,:)
    Character(len=80) :: SubIndexString(2)

    !Get SubSystem1 Density matrix as part of full D
    !or build new from scratch
    !config%SubSystemDensity
    !The SubSetDensity does seem to be closer to the right 
    !Density - better energy and smaller gradient 
    !but for some reason it converges in more iterations
    !so this have been deactivated

    SubIndexString(1) = ls%input%Molecule%SubSystemLabel(1)
    SubIndexString(2) = ls%input%Molecule%SubSystemLabel(2)
    OnMaster = .TRUE.
    lupri = config%lupri
    luerr = config%luerr
    NAtoms = config%Molecule%nAtoms
    nelectrons = config%Molecule%nelectrons
    
    call mem_alloc(Phantom,NAtoms)
    !Get SubSystem1 
    nelectronsSub = 0 
    Do i = 1,NAtoms
       Phantom(i) = ls%input%Molecule%Atom(i)%Phantom
       IF(ls%input%Molecule%Atom(i)%SubSystemIndex.NE.subIndex)THEN
          ls%input%Molecule%Atom(i)%Phantom = .TRUE.
       ELSE
          nelectronsSub = nelectronsSub + INT(ls%input%Molecule%Atom(i)%Charge)
       ENDIF
    Enddo
    ls%input%Molecule%nelectrons = nelectronsSub
    config%Molecule%nelectrons = nelectronsSub
    config%Integral%nelectrons = nelectronsSub
    config%decomp%nocc = nelectronsSub/2
    config%diag%nocc = nelectronsSub/2

    WRITE(lupri,*)'Subsystem ',TRIM(SubIndexString(SubIndex))
    CALL PRINT_MOLECULEINFO(LUPRI,ls%input%MOLECULE,ls%input%BASIS,ls%input%DALTON%MOLPRINT)
    CALL PRINT_MOLECULE_AND_BASIS(LUPRI,ls%input%MOLECULE,ls%input%BASIS%BINFO(REGBASPARAM))

    IF(config%decomp%nactive.NE.0)call lsquit('counter poise require closed shell',-1)
    
    nbast = D(1)%nrow
    IF(config%SubSystemDensity)THEN
       call mem_alloc(Dfull,nbast,nbast)
       call mem_alloc(D1full,nbast,nbast)
       call mat_to_full(D(1),1.0E0_realk,Dfull)
       call ls_dzero(D1full,nbast*nbast)
       !Get SubSystem1 Density matrix as part of full D
       R = ls%input%basis%binfo(RegBasParam)%labelindex
       nbast2 = 0
       do j=1, nAtoms
          IF(ls%input%MOLECULE%ATOM(j)%pointcharge)CYCLE
          IF(R.EQ.0)THEN
             jcharge = INT(ls%input%MOLECULE%ATOM(j)%charge) 
             jtype = ls%input%basis%binfo(RegBasParam)%chargeindex(jcharge)
          ELSE
             jtype = ls%input%MOLECULE%ATOM(j)%IDtype(1)
          ENDIF
          norbJ = ls%input%basis%binfo(RegBasParam)%ATOMTYPE(jtype)%ToTnorb
          IF(ls%input%Molecule%Atom(j)%SubSystemIndex.EQ.subIndex)THEN
             
             nbast1 = 0
             do i=1, nAtoms
                IF(ls%input%MOLECULE%ATOM(i)%pointcharge)CYCLE
                IF(R.EQ.0)THEN
                   icharge = INT(ls%input%MOLECULE%ATOM(i)%charge) 
                   itype = ls%input%basis%binfo(RegBasParam)%chargeindex(icharge)
                ELSE
                   itype = ls%input%MOLECULE%ATOM(i)%IDtype(1)
                ENDIF
                norbI = ls%input%basis%binfo(RegBasParam)%ATOMTYPE(itype)%ToTnorb
                IF(ls%input%Molecule%Atom(i)%SubSystemIndex.EQ.subIndex)THEN
                   
                   do jOrb=1, nOrbJ
                      do iOrb=1, nOrbI
                         D1full(nbast1+iOrb,nbast2+jOrb)=Dfull(nbast1+iOrb,nbast2+jOrb)
                      enddo
                   enddo
                   nbast1 = nbast1 + nOrbI
                ELSE
                   nbast1 = nbast1 + nOrbI
                ENDIF
             enddo
             nbast2 = nbast2 + nOrbJ          
          ELSE
             nbast2 = nbast2 + nOrbJ          
          ENDIF
       enddo
       call mem_dealloc(Dfull)
       call mat_init(D1(1),nbast,nbast)
       call mat_set_from_full(D1full,1.0E0_realk,D1(1))
       call mem_dealloc(D1full)
    ELSE
       call mat_init(D1(1),nbast,nbast)
    ENDIF

    IF(config%SubSystemDensity)THEN
       INQUIRE(file='dens.restart',EXIST=dens_exsist) 
       IF(dens_exsist)THEN
          restart_lun = -1  !initialization
          call lsopen(restart_lun,'dens.restart','OLD','UNFORMATTED')
          rewind restart_lun
          call mat_init(TMP,nbast,nbast)
          call mat_read_from_disk(restart_lun,TMP,OnMaster)
          call mat_read_info_from_disk(restart_lun,gcbasis)
          
          dens_lun = -1
          call lsopen(dens_lun,'dens.restart.Full','UNKNOWN','UNFORMATTED')
          rewind dens_lun
          call mat_write_to_disk(dens_lun,TMP,OnMaster)
          call mat_write_info_to_disk(dens_lun,gcbasis)
          call mat_free(TMP)
          
          call lsclose(dens_lun,'KEEP')
          call lsclose(restart_lun,'KEEP')
       ENDIF
       
       restart_lun = -1  !initialization
       call lsopen(restart_lun,'dens.restart','OLD','UNFORMATTED')
       rewind restart_lun
       call mat_write_to_disk(restart_lun,D1(1),OnMaster)
       call mat_write_info_to_disk(restart_lun,gcbasis)
       call lsclose(restart_lun,'KEEP')
       
       !Make sure that the get_energy uses McWeeny purification of dens.restart
       CFG_restart = config%diag%CFG_restart
       CFG_purifyrestart = config%diag%CFG_purifyrestart
       config%diag%CFG_restart = .TRUE.
       config%diag%CFG_purifyrestart = .TRUE.
       !Get Energy
    ELSE
       !restart density should not be used - this is the full density
       CFG_restart = config%diag%CFG_restart
       CFG_purifyrestart = config%diag%CFG_purifyrestart
       config%diag%CFG_restart = .FALSE.
       config%diag%CFG_purifyrestart = .FALSE.
    ENDIF

    IF(config%doDEC)call mat_init(CMO,H1%nrow,H1%ncol)
    OnlyOccPart = DECinfo%OnlyOccPart
    IF(config%doDEC.AND.(.not. DECinfo%full_molecular_cc))THEN
     IF(.NOT.DECinfo%OnlyOccPart)THEN
      WRITE(lupri,'(A)')'Note due to the presence of Phantom atoms in the SubSystem calculation'
      WRITE(lupri,'(A)')'The .ONLYOCCPART keyword is activated during the SubSystem calculation'
     ENDIF
     DECinfo%OnlyOccPart  = .TRUE.
    ENDIF
    call Get_Energy(Esub,Etmp,config,H1,F,D1,S,ls,CMO,Natoms,lupri,luerr)
    IF(config%doDEC)call mat_free(CMO)
    DECinfo%OnlyOccPart = OnlyOccPart

    config%diag%CFG_restart = CFG_restart
    config%diag%CFG_purifyrestart = CFG_purifyrestart
    IF(config%SubSystemDensity)THEN
       !Revert the dens.restart
       IF(dens_exsist)THEN
          restart_lun = -1  !initialization
          call lsopen(restart_lun,'dens.restart.Full','OLD','UNFORMATTED')
          rewind restart_lun
          call mat_init(TMP,nbast,nbast)
          call mat_read_from_disk(restart_lun,TMP,OnMaster)
          call mat_read_info_from_disk(restart_lun,gcbasis)
          
          dens_lun = -1
          call lsopen(dens_lun,'dens.restart','UNKNOWN','UNFORMATTED')
          rewind dens_lun
          call mat_write_to_disk(dens_lun,TMP,OnMaster)
          call mat_write_info_to_disk(dens_lun,gcbasis)
          call mat_free(TMP)
          
          call lsclose(dens_lun,'KEEP')
          call lsclose(restart_lun,'KEEP')
       ENDIF
       
       restart_lun = -1  !initialization
       call lsopen(restart_lun,'dens.restart.Sub'//TRIM(SubIndexString(SubIndex)),&
            & 'UNKNOWN','UNFORMATTED')
       rewind restart_lun
       call mat_write_to_disk(restart_lun,D1(1),OnMaster)
       call mat_write_info_to_disk(restart_lun,gcbasis)
       call lsclose(restart_lun,'KEEP')
    ENDIF
    call mat_free(D1(1))

    Do i = 1,NAtoms
       ls%input%Molecule%Atom(i)%Phantom = Phantom(i)
    Enddo
    call mem_dealloc(Phantom)
    ls%input%Molecule%nelectrons = nelectrons
    config%Molecule%nelectrons = nelectrons
    config%Integral%nelectrons = nelectrons
    config%decomp%nocc = nelectrons/2
    config%diag%nocc = nelectrons/2

  end SUBROUTINE SubSystemEnergy
!
end module InteractionEnergyMod




