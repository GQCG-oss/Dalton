!> @file
!> Full molecular calculation of RI-MP2-F12

module fullrimp2f12 

#ifdef MOD_UNRELEASED 

use precision
use typedeftype!,only:lsitem
use typedef
use dec_typedef_module
use matrix_module
use matrix_operations
use memory_handling
use LSTIMING
!use dec_fragment_utils
use CABS_operations
use ccintegrals,only:get_ao_fock
use f12ri_util_module,only: GeneralTwo4CenterF12RICoef1112, &
       & GeneralTwo4CenterF12RICoef1223,F12RIB9,F12RIB9MPI, &
       & F12RIB8,F12RIB8MPI,F12RIB7,F12RIB7MPI

!#ifdef MOD_UNRELEASED
use f12_routines_module,only: &
     & ContractOne4CenterF12IntegralsRI,&
     & ContractOne4CenterF12IntegralsRI2,&
     & ContractOne4CenterF12IntegralsRI2_nc,&
     & ContractOne4CenterF12IntegralsRobustRI,&
     & ContractOne4CenterF12IntegralsRIB23,&
     & ContractTwo4CenterF12IntegralsRI_pf,&
     & ContractTwo4CenterF12IntegralsRIX,&
     & ContractTwo4CenterF12IntegralsRIX_nc,&
     & ContractTwo4CenterF12IntegralsRIX_ncMPI,&
     & ContractTwo4CenterF12IntegralsRI2V3V4,&
     & ContractTwo4CenterF12IntegralsRIC_pf,&
     & ContractTwo4CenterF12IntegralsRIX3X4_nc,&
     & ContractTwo4CenterF12IntegralsRIX3X4_ncMPI,&
     & ContractTwo4CenterF12IntegralsRIB4,&
     & ContractTwo4CenterF12IntegralsRIB4MPI,&
     & ContractTwo4CenterF12IntegralsRIB5,&
     & ContractTwo4CenterF12IntegralsRIB5MPI,&
     & ContractTwo4CenterF12IntegralsRIB6,&
     & ContractTwo4CenterF12IntegralsRIB6MPI,&
!     & ContractTwo4CenterF12IntegralsRIB7,&
!     & ContractTwo4CenterF12IntegralsRIB7MPI,&
!     & ContractTwo4CenterF12IntegralsRIB8,&
!     & ContractTwo4CenterF12IntegralsRIB8MPI,&
!     & ContractTwo4CenterF12IntegralsRIB9,&
!     & ContractTwo4CenterF12IntegralsRIB9MPI,&
     & get_F12_mixed_MO_Matrices, free_f12_mixed_mo_matrices,&
     & MO_transform_AOMatrix

use IntegralInterfaceMOD
use ri_util_module
#ifdef VAR_MPI
use lsmpi_op,only: mpicopy_lsitem
use decmpi_module,only: mpi_bcast_fullmolecule
use lsmpi_type,only:ls_mpiInitBuffer,ls_mpiFinalizeBuffer,&
     & LSMPIBROADCAST,MPI_COMM_LSDALTON 
#endif
!#endif 

public :: full_canonical_rimp2_f12, lsmpi_matrix_bufcopy,&
     & NaturalLinearScalingF12Terms

private

#endif

contains
#ifdef MOD_UNRELEASED 

!> \brief Calculate canonical MP2 energy for full molecular system
!> keeping full AO integrals in memory. Only for testing.
!> \author Thomas Kjaergaard
!> \date 2015
subroutine full_canonical_rimp2_f12(MyMolecule,MyLsitem,Dmat,mp2f12_energy)
   implicit none
   !> Full molecule info
   type(fullmolecule), intent(inout) :: MyMolecule 
   !> Lsitem structure
   type(lsitem), intent(inout) :: mylsitem
   !> HF density matrix
   type(matrix),intent(in) :: Dmat
   !> MP2-F12 correlation energy
   real(realk),intent(inout) :: mp2f12_energy
   !local variables
   integer :: nbasis,nocc,nvirt,ncabsAO,ncabsMO
   !> Canonical MP2 correlation energy
   real(realk) :: mp2_energy
   real(realk) :: E21,Econt(1),E23
   real(realk) :: ExchangeF12V1,CoulombF12V1
   real(realk) :: ExchangeF12X1,CoulombF12X1
   real(realk) :: ExchangeF12B1,CoulombF12B1
   real(realk) :: E_21,E_22,E_23, E_F12, E_21C, E_21Ctmp
   real(realk) :: EV1,EV2,EV3,EV4,EV5,EX1,EX2,EX3,EX4
   real(realk) :: EB1,EB2,EB3,EB4,EB5,EB6,EB7,EB8,EB9
   real(realk) :: EV1tmp,EV2tmp,EV3tmp,EV4tmp,EV5tmp,EX1tmp,EX2tmp,EX3tmp,EX4tmp
   real(realk) :: EB1tmp,EB2tmp,EB3tmp,EB4tmp,EB5tmp,EB6tmp,EB7tmp,EB8tmp,EB9tmp
   real(realk) :: TS,TE,TS2,TE2
   integer :: i,j,a,b,p,q,c,m,mynum,nAtoms,lupri,nbuf1,inode,numnodesstd
   integer(kind=long) :: nsize,nsize2
   integer(kind=ls_mpik) :: node,numnodes
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
   !>   Singles correction
   type(matrix) :: Fic
   type(matrix) :: Fcd
   real(realk)  :: ES2
   !========================================================
   ! RI variables
   !========================================================
   integer :: nAux,NBA,N,K,ncore,NBA2
   real(realk),pointer :: CalphaR(:),CalphaG(:),CalphaF(:),CalphaD(:),CalphaCvirt(:), CalphaT(:)
   real(realk),pointer :: CalphaRcabsMO(:),CalphaGcabsAO(:),CalphaX(:),CalphaCcabs(:), CalphaP(:)
   real(realk),pointer :: CalphaGcabsMO(:),CalphaXcabsAO(:), CalphACcabsT(:), CalphaCocc(:), CalphaCoccT(:)
   real(realk),pointer :: CalphaTvirt(:),UmatTmp(:,:)
   real(realk),pointer :: CalphaTmp(:),EpsOcc(:),EpsVirt(:)
   real(realk),pointer :: Cfull(:,:),ABdecompR(:,:),ABdecompG(:,:),ABdecompC(:,:),ABdecompTvirt(:,:)
   real(realk),pointer :: ABdecompF(:,:),Umat(:,:),Rtilde(:,:),ABdecompX(:,:)
   real(realk),pointer :: CalphaMPI(:),CalphaMPI2(:)
   logical :: master,wakeslaves,ABdecompCreateR,ABdecompCreateG,ABdecompCreateF,ABdecompCreateC,ABdecompCreateTvirt
   logical :: FORCEPRINT,use_bg_buf,ABdecompCreateX
   character :: intspec(5)
   type(matrix) :: CMO_CABS,CMO_RI
   !========================================================
   ! Additional variables
   !========================================================
   integer :: offset, noccfull, nocv,AuxMPIstartMy,iAuxMPIextraMy,AuxMPIstartMPI,iAuxMPIextraMPI
   integer,pointer :: nAuxMPI(:)
   real(realk),pointer :: Taibj(:,:,:,:) !amplitudes not integrals
   real(realk),pointer :: gmo(:,:,:,:)
   real(realk),pointer :: gao(:,:,:,:)
   real(realk) :: eps,factor
   real(realk),pointer :: Fkj(:,:)
   real(realk),pointer :: Co(:,:)
#ifdef VAR_MPI
   real(realk) :: lsmpibufferRIMP2(20)
   lsmpibufferRIMP2=0.0E0_realk
#endif
   if(MyMolecule%mem_distributed)then
      call lsquit("ERROR(full_canonical_rimp2_f12): does not work with PDM type fullmolecule",-1)
   endif
 
   E_21C = 0.0E0_realk
   MP2_energy = mp2f12_energy
   mp2f12_energy = 0.0E0_realk
   use_bg_buf = .FALSE.

   lupri = DECinfo%output
#ifdef VAR_TIME
    FORCEPRINT = .TRUE.
#else
    FORCEPRINT = .FALSE.
#endif    
    call LSTIMER('START ',TS,TE,DECinfo%output,ForcePrint)

#ifdef VAR_MPI
    master= (infpar%mynum == infpar%master)
    mynum = infpar%mynum
    numnodes = infpar%nodtot
    wakeslaves = infpar%nodtot.GT.1
    if(.NOT.master)lupri = 6
#else
    ! If MPI is not used, consider the single node to be "master"
    master=.true.
    mynum = 0
    numnodes = 1
    wakeslaves = .false.
#endif
    numnodesstd = numnodes
   call LSTIMER('START ',TS2,TE2,DECinfo%output,ForcePrint)

   ! Init stuff
   ! **********
   ncore  = MyMolecule%ncore
   nbasis = MyMolecule%nbasis
   nvirt  = MyMolecule%nvirt
   naux   = MyMolecule%nauxbasis
   nAtoms = MyMolecule%nAtoms

   !noccfull = nocc
   !IF(DECinfo%frozencore)call lsquit('RI-MP2-F12 frozen core not implemented',-1)

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
   nocv = nvirt + noccfull

   IF(master)THEN
      call determine_CABS_nbast(ncabsAO,ncabsMO,mylsitem%setting,DECinfo%output)
      call mat_init(CMO_CABS,nCabsAO,ncabsMO)
      call build_CABS_MO(CMO_CABS,nCabsAO,mylsitem%SETTING,lupri)    
      call mat_init(CMO_RI,nCabsAO,nCabsAO)
      call build_RI_MO(CMO_RI,nCabsAO,mylsitem%SETTING,lupri)
   ENDIF

   ! ***********************************************************
   !   Constructing Coefficient matrices 
   ! ***********************************************************   
   ! Cocc
   if(DECinfo%frozencore) then
      call mem_alloc(Co,nbasis,nocc)
      do j=1,nocc
         do i=1,nbasis
            Co(i,j) = MyMolecule%Co%elm2(i,MyMolecule%ncore+j)
         enddo
      enddo
   else
      call mem_alloc(Co,nbasis,nocc)
      call dcopy(nbasis*nocc, MyMolecule%Co%elm2, 1, Co, 1)
   end if

   ! Fkj
   if(DECinfo%frozencore) then
      call mem_alloc(Fkj,nocc,nocc)
      do j=1,nocc
         do i=1,nocc
            Fkj(i,j) = MyMolecule%oofock%elm2(MyMolecule%ncore+i,MyMolecule%ncore+j)
         end do
      end do
   else
      call mem_alloc(Fkj,nocc,nocc)
      call dcopy(nocc*nocc, MyMolecule%oofock%elm2, 1, Fkj, 1)
   endif

   ! Cfull
   call mem_alloc(Cfull,nbasis,nocv)
   do J=1,noccfull
      do I=1,nbasis
         Cfull(I,J) = MyMolecule%Co%elm2(I,J)
      enddo
   enddo
   do P=1,nvirt
      do I=1,nbasis
         Cfull(I,noccfull+P) = MyMolecule%Cv%elm2(I,P)
      enddo
   enddo

   ! ***********************************************************
   !   Printing Input variables 
   ! ***********************************************************
   if(DECinfo%F12debug) then
      IF(master)THEN
         !ncabsAO, ncabsMO not set on slave
         print *, "-------------------------------------------------"
         print *, "     full_rimp2f12.F90                           "
         print *, "-------------------------------------------------"
         print *, "nbasis:   ", nbasis
         print *, "nocc:     ", nocc
         print *, "nvirt:    ", nvirt
         print *, "-------------------------------------------------"
         print *, "noccfull  ", noccfull
         print *, "ncabsAO   ", ncabsAO
         print *, "ncabsMO   ", ncabsMO
         print *, "-------------------------------------------------"
         print *, "offest:   ", offset
         print *, "ncore:    ", ncore
      ENDIF
   end if

   IF(naux.EQ.0)call lsquit('Error no Aux functions in full_canonical_rimp2_f12',-1)

   IF(master)THEN
      call get_F12_mixed_MO_Matrices(MyLsitem,MyMolecule,Dmat,nbasis,ncabsAO,&
           & nocc,noccfull,nvirt,ncabsMO,HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp,Fic,Fcd)
   ENDIF

   call LSTIMER('FULLRIMP2:Init',TS2,TE2,DECinfo%output,ForcePrint)
   !==================================================================
   != Step 1:  Fijkl,Xijkl,Dijkl                                     =
   !=          corresponding to V1,X1,B1                             =
   != These are special since                                        =
   != 1. the have a very simple structure(B1 does require robust DF) =
   != 2. They could all be done in a Linear scaling way outside DEC  =
   !=    without density fitting.                                    =
   != 3. The intermediates are only used once                        =
   != 4. Due to noccEOS,noccEOS,noccEOS,noccEOS very small mem req   =
   != Note                                                           =
   != Fijkl:                                                         =
   != The Gaussian geminal divided by the Coulomb operator g/r12     =
   != Xijkl                                                          =
   != The Gaussian geminal squared g^2                               =
   != Dijkl                                                          =
   != The double commutator [[T,g],g] with g = Gaussian geminal      =
   != Since the integral (alpha|[[T,g],g]|beta) is not positive      =
   != definite this term is done using robust density fitting        =
   != Done according to Eq. 87 of                                    =
   != J Comput Chem 32: 2492–2513, 2011                              =
   !==================================================================
   IF(master)THEN
      write(DECinfo%output,'(/,a)') ' ================================================ '
      write(DECinfo%output,'(a)')   '            FULL-RI-MP2F12 ENERGY TERMS            '
      write(DECinfo%output,'(a,/)') ' ================================================ '
      write(*,'(/,a)') ' ================================================ '
      write(*,'(a)')   '           FULL-RI-MP2F12 ENERGY TERMS             '
      write(*,'(a,/)') ' ================================================ '
   ENDIF
   call LSTIMER('START ',TS2,TE2,DECinfo%output,ForcePrint)

#ifdef VAR_MPI
   StartUpSlaves: if(wakeslaves .and. master) then
      ! Wake up slaves to do the job: slaves awoken up with (RIMP2FULL)
      ! and call full_canonical_rimp2_slave which communicate info 
      ! then calls full_canonical_rimp2.
      call ls_mpibcast(RIMP2F12FULL,infpar%master,infpar%lg_comm)
      ! Communicate fragment information to slaves
      call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call mpicopy_lsitem(MyLsitem,infpar%lg_comm)
      call lsmpi_matrix_bufcopy(Dmat,master)
      call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call mpi_bcast_fullmolecule(MyMolecule)    
   endif StartUpSlaves

   call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
   call lsmpi_matrix_bufcopy(HJir,master)
   call lsmpi_matrix_bufcopy(Krr,master)
   call lsmpi_matrix_bufcopy(Frr,master)
   call lsmpi_matrix_bufcopy(Fac,master)
   call lsmpi_matrix_bufcopy(Fpp,master)
   call lsmpi_matrix_bufcopy(Fii,master)
   call lsmpi_matrix_bufcopy(Fmm,master)
   call lsmpi_matrix_bufcopy(Frm,master)
   call lsmpi_matrix_bufcopy(Fcp,master)
   call lsmpi_matrix_bufcopy(Fic,master)
   call lsmpi_matrix_bufcopy(Fcd,master)
   call lsmpi_matrix_bufcopy(CMO_CABS,master)
   IF(.NOT.Master)THEN
      nCabsAO = CMO_CABS%nrow
      nCabsMO = CMO_CABS%ncol
   ENDIF
   call lsmpi_matrix_bufcopy(CMO_RI,master)
   call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
#endif

   intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
   intspec(2) = 'R' !Regular AO basis function on center 3
   intspec(3) = 'R' !Regular AO basis function on center 4

   IF(DECinfo%NaturalLinearScalingF12TermsV1)THEN
      !This energy contribution have already been calculated so we extract the information 
      EV1 = MyMolecule%EF12NLSV1
#ifdef VAR_MPI 
      IF(master)THEN
         lsmpibufferRIMP2(2)=EV1
      ENDIF
#else
      mp2f12_energy = mp2f12_energy  + EV1
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V1,LS) = ',EV1
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V1,LS) = ',EV1
#endif         
   ELSE
      call mem_alloc(ABdecompF,nAux,nAux)
      ABdecompCreateF = .TRUE.
      ! Calculate the Fitting Coefficients (alpha|F|ij)
      intspec(4) = 'F' !The Gaussian geminal divided by the Coulomb operator g/r12 (GGemCouOperator)
      intspec(5) = 'F' !The Gaussian geminal divided by the Coulomb operator g/r12 (GGemCouOperator)
      !Unique: CalphaF(NBA,nocc,nocc)
      call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,FORCEPRINT, &
           & wakeslaves,Co,nocc,Co,nocc,mynum,numnodesstd,CalphaF,NBA,ABdecompF,&
           & ABdecompCreateF,intspec,use_bg_buf)
      call mem_dealloc(ABdecompF)
      ABdecompCreateF = .FALSE.

      !perform this suborutine on the GPU (async)  - you do not need to wait for the results
      call ContractOne4CenterF12IntegralsRI(NBA,nocc,CalphaF,CoulombF12V1,ExchangeF12V1)

      !The minus is due to the Valeev factor
      EV1 = -1.0E0_realk*((5.0E0_realk*0.25E0_realk)*CoulombF12V1-ExchangeF12V1*0.25E0_realk)
#ifdef VAR_MPI 
      lsmpibufferRIMP2(2)=EV1       !we need to perform a MPI reduction at the end 
#else
      mp2f12_energy = mp2f12_energy  + EV1
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V1,RI) = ',EV1
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V1,RI) = ',EV1
#endif   
      call mem_dealloc(CalphaF)
   ENDIF

   IF(DECinfo%NaturalLinearScalingF12TermsX1)THEN
      !This energy contribution have already been calculated so we extract the information 
      EX1 = MyMolecule%EF12NLSX1
#ifdef VAR_MPI 
      IF(master)THEN
         lsmpibufferRIMP2(3)=EX1
      ENDIF
#else
      mp2f12_energy = mp2f12_energy  + EX1
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X1,LS) = ',EX1
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X1,LS) = ',EX1
#endif         
   ELSE
      call mem_alloc(ABdecompX,nAux,nAux)
      ABdecompCreateX = .TRUE.
      !Calculate the Fitting Coefficients (alpha|g^2|ij) 
      intspec(4) = '2' !The Gaussian geminal operator g^2 (GGemCouOperator)
      intspec(5) = '2' !The Gaussian geminal operator g^2 (GGemCouOperator)
      !Unique: CalphaX(NBA,nocc,nocc)
      call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
           & FORCEPRINT,wakeslaves,Co,nocc,Co,nocc,&
           & mynum,numnodesstd,CalphaX,NBA,ABdecompX,ABdecompCreateX,intspec,use_bg_buf)
      ABdecompCreateX = .FALSE.

      !perform this suborutine on the GPU (async)  - you do not need to wait for the results
      IF(DECinfo%use_canonical)THEN
         call ContractOne4CenterF12IntegralsRI2(NBA,nocc,CalphaX,Fii%elms,CoulombF12X1,ExchangeF12X1)
      ELSE
         nsize = NBA*nocc*nocc
         call mem_alloc(CalphaT,nsize)     ! G_ij = C_ik F_kj 
         M = NBA*nocc         !rows of Output Matrix
         K = nocc             !summation dimension
         N = nocc             !columns of Output Matrix
         !call dgemm('N','N',M,N,K,1.0E0_realk,CalphaX,M,Fmm%elms,K,0.0E0_realk,CalphaT,M)
         call dgemm('N','N',M,N,K,1.0E0_realk,CalphaX,M,Fkj,K,0.0E0_realk,CalphaT,M)
         call ContractOne4CenterF12IntegralsRI2_nc(NBA,nocc,CalphaX,CalphaT,CoulombF12X1,ExchangeF12X1)
         call mem_dealloc(CalphaT)
      ENDIF
      !minus is due to the overall minus from equation (41) and (42) due to
      !contribution from the \bar{B}_{ij}^{ij}
      EX1 = -1.0E0_realk*(0.21875E0_realk*CoulombF12X1 + 0.03125E0_realk*ExchangeF12X1)
#ifdef VAR_MPI 
      lsmpibufferRIMP2(3)=EX1       !we need to perform a MPI reduction at the end 
#else
      mp2f12_energy = mp2f12_energy  + EX1
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X1,RI) = ',EX1
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X1,RI) = ',EX1
#endif   
   ENDIF

   !We need CalphaR both for B1 but also other terms
   call mem_alloc(ABdecompR,nAux,nAux)
   ABdecompCreateR = .TRUE.
   !We need CalphaR(NBA,nocc,nocc) but this is a subset of the CalphaR(NBA,nocc,nbasis) we need later
   !so we calculate the full CalphaR(NBA,nocc,nbasis)
   intspec(4) = 'C' !Regular Coulomb operator 1/r12
   intspec(5) = 'C' !The metric operator = Regular Coulomb operator 1/r12
   !Build the G coefficient of Eq. 90 of J Comput Chem 32: 2492–2513, 2011
   !Unique: CalphaR(NBA,nocc,nocv)
   call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,Co,nocc,Cfull,nocv,&
        & mynum,numnodesstd,CalphaR,NBA,ABdecompR,ABdecompCreateR,intspec,use_bg_buf)
   ABdecompCreateR = .FALSE.

   IF(DECinfo%NaturalLinearScalingF12TermsB1)THEN
      !This energy contribution have already been calculated so we extract the information 
      EB1 = MyMolecule%EF12NLSB1
#ifdef VAR_MPI 
      IF(master)THEN
         lsmpibufferRIMP2(1)=EB1
      ENDIF
      IF(wakeslaves)THEN
         nbuf1=numnodes
         call mem_alloc(nAuxMPI,nbuf1)
         call BuildnAuxMPIUsedRI(nAux,numnodesstd,nAuxMPI)      
      ENDIF
#else
      mp2f12_energy = mp2f12_energy  + EB1
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B1,LS) = ',EB1
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B1,LS) = ',EB1
#endif         
   ELSE
      !Calculate the Fitting Coefficients (alpha|[[T,g],g]|ij) 
      intspec(4) = 'D' !The double commutator [[T,g],g] 
      intspec(5) = 'C' !The metric operator = Regular Coulomb operator 1/r12
      !Build the R coefficient of Eq. 90 of J Comput Chem 32: 2492–2513, 2011
      !Unique: CalphaD(NBA,nocc,nocc)
      call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
           & FORCEPRINT,wakeslaves,Co,nocc,Co,nocc,&
           & mynum,numnodesstd,CalphaD,NBA,ABdecompR,ABdecompCreateR,intspec,use_bg_buf)

      !Build the U matrix in Eq. 88 of J Comput Chem 32: 2492–2513, 2011
      call mem_alloc(Umat,nAux,nAux)
      !perform this suborutine on the GPU (Async)
      call Build_RobustERImatU(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
           & FORCEPRINT,wakeslaves,Co,nocc,Co,nocc,&
           & mynum,numnodesstd,ABdecompR,'D',Umat)      
#ifdef VAR_MPI
      !Build the R tilde coefficient of Eq. 89 of J Comput Chem 32: 2492–2513, 2011
      !We need to do:
      !CalphaD(NBA,nocc,nocc) = -0.5*CalphaD(NBA,nocc,nocc) + Loop_NBA2 Umat(NBA,NBA2)*CalphaR(NBA2,nocc,nocc)
      !Where NBA is the Auxiliary assigned to this node, while NBA2 can be assigned to another node
      IF(wakeslaves)THEN
         nbuf1=numnodes
         call mem_alloc(nAuxMPI,nbuf1)
         call BuildnAuxMPIUsedRI(nAux,numnodesstd,nAuxMPI)      
         call BuildnAuxMPIUsedRIinfo(nAux,numnodesstd,mynum,AuxMPIstartMy,iAuxMPIextraMy)
         DO inode = 1,numnodes
            nbuf1 = nAuxMPI(inode)
            NBA2 = nAuxMPI(inode)
            nsize = nbuf1*nocc*nocv
            nbuf1 = NBA2
            IF(inode.EQ.1)THEN
               factor = -0.5E0_realk
            ELSE
               factor = 1.0E0_realk
            ENDIF
            call BuildnAuxMPIUsedRIinfo(nAux,numnodesstd,inode-1,AuxMPIstartMPI,iAuxMPIextraMPI)
            IF(use_bg_buf)THEN
               call mem_pseudo_alloc(UmatTmp,NBA*i8,NBA2*i8)
            ELSE
               call mem_alloc(UmatTmp,NBA,NBA2)
            ENDIF
            Call BuilDUmatTmpRIF12(Umat,nAux,UmatTmp,NBA,NBA2,AuxMPIstartMy,iAuxMPIextraMy,&
                 & AuxMPIstartMPI,iAuxMPIextraMPI)
            IF(mynum.EQ.inode-1)THEN
               !I Bcast My Own CalphaR (the first NBA*nocc*nocc part )
               node = mynum
               call ls_mpibcast(CalphaR,nsize,node,infpar%lg_comm)
               M = NBA          !rows of Output Matrix
               N = nocc*nocc    !columns of Output Matrix
               K = NBA          !summation dimension
               call dgemm('N','N',M,N,K,1.0E0_realk,UmatTmp,M,CalphaR(1+offset*NBA*nocc),K,factor,CalphaD,M)
            ELSE
               node = inode-1
               !recieve
               IF(use_bg_buf)THEN
                  call mem_pseudo_alloc(CalphaMPI,nsize)
               ELSE
                  call mem_alloc(CalphaMPI,nsize)
               ENDIF
               call ls_mpibcast(CalphaMPI,nsize,node,infpar%lg_comm)
               M = NBA          !rows of Output Matrix
               N = nocc*nocc    !columns of Output Matrix
               K = NBA2          !summation dimension
               call dgemm('N','N',M,N,K,1.0E0_realk,UmatTmp,M,CalphaMPI(1+offset*NBA2*nocc),K,factor,CalphaD,M)
               IF(use_bg_buf)THEN
                  call mem_pseudo_dealloc(CalphaMPI)
               ELSE
                  call mem_dealloc(CalphaMPI)
               ENDIF
            ENDIF
            IF(use_bg_buf)THEN
               call mem_pseudo_dealloc(UmatTmp)
            ELSE
               call mem_dealloc(UmatTmp)
            ENDIF
         ENDDO
      ELSE
         M = NBA          !rows of Output Matrix
         K = NBA          !summation dimension
         N = nocc*nocc    !columns of Output Matrix
         call dgemm('N','N',M,N,K,1.0E0_realk,Umat,M,CalphaR(1+offset*NBA*nocc),K,-0.5E0_realk,CalphaD,M)
      ENDIF
#else
      !Build the R tilde coefficient of Eq. 89 of J Comput Chem 32: 2492–2513, 2011
      !perform this suborutine on the GPU (Async)
      !note CalphaR is actual of dimensions (NBA,nocc,nbasis) but here we only access
      !the first part (NBA,nocc,nocc) 
      M = NBA          !rows of Output Matrix
      K = NBA          !summation dimension
      N = nocc*nocc    !columns of Output Matrix
      call dgemm('N','N',M,N,K,1.0E0_realk,Umat,M,CalphaR(1+offset*NBA*nocc),K,-0.5E0_realk,CalphaD,M)
#endif
      call mem_dealloc(Umat)

      !CalphaD is now the R tilde coefficient of Eq. 89 of J Comput Chem 32: 2492–2513, 2011
      !perform this suborutine on the GPU (Async)
      call ContractOne4CenterF12IntegralsRobustRI(NBA,offset,nocc,nocv,CalphaD,CalphaR,EB1)
#ifdef VAR_MPI 
      lsmpibufferRIMP2(1)=EB1       !we need to perform a MPI reduction at the end 
#else
      mp2f12_energy = mp2f12_energy  + EB1
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B1,RI) = ', EB1
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B1,RI) = ', EB1
#endif   
      call mem_dealloc(CalphaD)
   ENDIF

   call LSTIMER('FULLRIMP2:Step1',TS2,TE2,DECinfo%output,ForcePrint)

   !==============================================================
   !=  B2: sum_c' (ic'|f12^2|jj) hJ_ic' - (jc'|f12^2|ij) hJ_ic'  =
   !=  B3: sum_c' (ii|f12^2|jc') hJ_jc' - (ij|f12^2|ic') hJ_ic'  =
   !==============================================================
   IF(DECinfo%NaturalLinearScalingF12TermsX1)THEN
      !Create CalphaX because this have not been made yet
      call mem_alloc(ABdecompX,nAux,nAux)
      ABdecompCreateX = .TRUE.
      !Calculate the Fitting Coefficients (alpha|g^2|ij) 
      intspec(4) = '2' !The Gaussian geminal operator g^2 (GGemCouOperator)
      intspec(5) = '2' !The Gaussian geminal operator g^2 (GGemCouOperator)
      !Unique: CalphaX(NBA,nocc,nocc) IF(NaturalLinearScalingF12TermsX1)
      call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
           & FORCEPRINT,wakeslaves,Co,nocc,Co,nocc,&
           & mynum,numnodesstd,CalphaX,NBA,ABdecompX,ABdecompCreateX,intspec,use_bg_buf)
      ABdecompCreateX = .FALSE.
   ENDIF

   intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
   intspec(2) = 'R' !Regular AO basis function on center 3
   intspec(3) = 'C' !CABS AO basis function on center 4
   intspec(4) = '2' !The f12 Operator
   intspec(5) = '2' !The f12 Operator
   !Unique: CalphaXcabsAO(NBA,nocc,ncabsAO) 
   call Build_CalphaMO2(mylsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,Co,nocc,CMO_RI%elms,ncabsAO,&
        & mynum,numnodesstd,CalphaXcabsAO,NBA,ABdecompX,ABdecompCreateX,intspec,use_bg_buf)   
   
   call ContractOne4CenterF12IntegralsRIB23(nBA,nocc,ncabsAO,CalphaXcabsAO,CalphaX,&
        & hJir%elms,1.0E0_realk,EB2,EB3)
#ifdef VAR_MPI 
   lsmpibufferRIMP2(4)=EB2       !we need to perform a MPI reduction at the end 
   lsmpibufferRIMP2(5)=EB3       !we need to perform a MPI reduction at the end 
#else
   !1.0E0_realk because that term has an overall pluss in Eqs. 25-26
   mp2f12_energy = mp2f12_energy  + EB2
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B2,RI) = ',EB2
   mp2f12_energy = mp2f12_energy  + EB3
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B3,RI) = ',EB3
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B2,RI) = ',EB2
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B3,RI) = ',EB3
#endif

   call mem_dealloc(CalphaXcabsAO)
   call mem_dealloc(CalphaX)
   call mem_dealloc(ABdecompX)
   !ABdecompR still exists
   call LSTIMER('FULLRIMP2:B2B3',TS2,TE2,DECinfo%output,ForcePrint)

   !=================================================================
   != Step 2: Ripjq*Gipjq+Rimjc*Gimjc+Rjmic*Gjmic                   =
   !=        +Gipjq*Gipjq+Gimjc*Gimjc+Gjmic*Gjmic                   =
   !=         corresponding to V2,V3,V4,X2,X3,X4                    =
   != These are special since                                       =
   != 1. the have (more or less) the same structure                 =
   != 2. They can be built from same intermediates                  =
   !=    CalphaR(alpha,i,p),CalphaR(alpha,i,c)                      =
   !=    CalphaG(alpha,i,p),CalphaG(alpha,i,c)                      =
   != 3. The intermediates are only used once                       =
   != 4. Due to noccEOS,noccEOS,noccEOS,noccEOS very small mem req  =
   !=================================================================

   call mem_alloc(ABdecompG,nAux,nAux)
   ABdecompCreateG = .TRUE.

   !==========================================================
   !=                                                        =
   != V2:                 Ripjq*Gipjq                        =
   != The Coulomb Operator Int multiplied with               =
   != The Gaussian geminal operator g                        =
   != Dim(nocc,nbasis,nocc,nbasis)                           =
   !=                                                        =
   !==========================================================
   intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
   intspec(2) = 'R' !Regular AO basis function on center 3
   intspec(3) = 'R' !Regular AO basis function on center 4
   intspec(4) = 'G' !The Gaussian geminal operator g
   intspec(5) = 'G' !The Gaussian geminal operator g
   !Unique: CalphaG(NBA,nocv,nocc)
   call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,Cfull,nocv,Co,nocc,&
        & mynum,numnodesstd,CalphaG,NBA,ABdecompG,ABdecompCreateG,intspec,use_bg_buf)

   !call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
   !     & FORCEPRINT,wakeslaves,Co,nocc,Cfull,nocv,&
   !     & mynum,numnodesstd,CalphaG,NBA,ABdecompG,ABdecompCreateG,intspec,use_bg_buf)
   ABdecompCreateG = .FALSE.
   
   !Do on GPU (Async) while the CPU starts calculating the next fitting Coef.
#ifdef VAR_MPI 
   nsize = NBA*nocv*nocc
   call mem_alloc(CalphaT,nsize)               ! G_qj = C_qk B_kj 
   M = nocv*NBA         !rows of Output Matrix
   K = nocc             !summation dimension
   N = nocc             !columns of Output Matrix
   call dgemm('N','N',M,N,K,1.0E0_realk,CalphaG,M,Fkj,K,0.0E0_realk,CalphaT,M)

   IF(wakeslaves)THEN
      EV2 = 0.0E0_realk
      EX2 = 0.0E0_realk
      DO inode = 1,numnodes
         nbuf1 = nAuxMPI(inode)
         NBA2 = nAuxMPI(inode)
         nsize = nbuf1*nocv*nocc
         IF(mynum.EQ.inode-1)THEN
            !I Bcast My Own CalphaG
            node = mynum
            IF(size(CalphaG).NE.nsize)call lsquit('MPI Bcast error in Full RIMP2F12 A',-1)
            call ls_mpibcast(CalphaG,nsize,node,infpar%lg_comm)
            call ContractTwo4CenterF12IntegralsRI_pf(nBA,nocc,nocv,CalphaR,CalphaG,NBA2,EV2tmp)
            call ContractTwo4CenterF12IntegralsRIX_nc(nBA,nocc,nocv,CalphaG,CalphaT,EX2tmp)
         ELSE
            node = inode-1
            !recieve
            IF(use_bg_buf)THEN
               call mem_pseudo_alloc(CalphaMPI,nsize)
            ELSE
               call mem_alloc(CalphaMPI,nsize)
            ENDIF
            call ls_mpibcast(CalphaMPI,nsize,node,infpar%lg_comm)
            call ContractTwo4CenterF12IntegralsRI_pf(nBA,nocc,nocv,CalphaR,CalphaMPI,NBA2,EV2tmp)
            call ContractTwo4CenterF12IntegralsRIX_ncMPI(nBA,nocc,nocv,CalphaMPI,NBA2,CalphaG,&
                 & CalphaT,EX2tmp)
            IF(use_bg_buf)THEN
               call mem_pseudo_dealloc(CalphaMPI)
            ELSE
               call mem_dealloc(CalphaMPI)
            ENDIF
         ENDIF
         EV2 = EV2 + EV2tmp
         EX2 = EX2 + EX2tmp
      ENDDO
   ELSE
      call ContractTwo4CenterF12IntegralsRI_pf(nBA,nocc,nocv,CalphaR,CalphaG,NBA,EV2)
      call ContractTwo4CenterF12IntegralsRIX_nc(nBA,nocc,nocv,CalphaG,CalphaT,EX2)  
   ENDIF
   lsmpibufferRIMP2(6)=EV2      !we need to perform a MPI reduction at the end 
   lsmpibufferRIMP2(7)=EX2      !we need to perform a MPI reduction at the end    
   call mem_dealloc(CalphaT)
#else
   call ContractTwo4CenterF12IntegralsRI_pf(nBA,nocc,nocv,CalphaR,CalphaG,NBA,EV2)
   mp2f12_energy = mp2f12_energy + EV2
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V2,RI) = ',EV2       
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V2,RI) = ',EV2
   !==========================================================
   !=                                                        =
   != X2: Gipjq*Gipjq                                        =
   != The Gaussian geminal operator Int multiplied with      =
   != The Gaussian geminal operator g                        =
   != Dim(nocc,nbasis,nocc,nbasis)                           =
   !=                                                        =
   !==========================================================
   !Do on GPU (Async) while the CPU starts calculating the next fitting Coef.
   !call ContractTwo4CenterF12IntegralsRIX(nBA,nocc,nbasis,CalphaG,Fii%elms,EX2)
   !mp2f12_energy = mp2f12_energy  + EX2
   !WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X2,RI) = ',EX2       
   !WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X2,RI) = ',EX2

   !Do on GPU (Async) while the CPU starts calculating the next fitting Coef.
   nsize = NBA*nocv*nocc
   call mem_alloc(CalphaT,nsize)               ! G_qj = C_qk B_kj 
   M = nocv*NBA         !rows of Output Matrix
   K = nocc             !summation dimension
   N = nocc             !columns of Output Matrix
!   IF(DECinfo%USE_CANONICAL)THEN
!      call ContractTwo4CenterF12IntegralsRIX(nBA,nocc,nocv,CalphaG,Fii%elms,EX2)
!   ELSE
      call dgemm('N','N',M,N,K,1.0E0_realk,CalphaG,M,Fkj,K,0.0E0_realk,CalphaT,M)
      call ContractTwo4CenterF12IntegralsRIX_nc(nBA,nocc,nocv,CalphaG,CalphaT,EX2)
!   ENDIF
   call mem_dealloc(CalphaT)
   mp2f12_energy = mp2f12_energy  + EX2
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X2,RI) = ',EX2       
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X2,RI) = ',EX2
#endif

   !==========================================================
   !=                                                        =
   != V3:                 Rimjc*Gimjc                        =
   != V4:                 Rjmic*Gjmic                        =
   != The Coulomb Operator Int multiplied with               =
   != The Gaussian geminal operator g                        =
   != Dim: (nocc,noccfull,nocc,ncabsMO)  need 4 Calphas      =
   !=                                                        =
   !==========================================================
   
   intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
   intspec(2) = 'R' !Regular AO basis function on center 3
   intspec(3) = 'C' !CABS AO basis function on center 4
   intspec(4) = 'C' !The Coulomb Operator
   intspec(5) = 'C' !The Coulomb Operator
   !Unique: CalphaRcabsMO(NBA,nocc,ncabsMO)
   call Build_CalphaMO2(mylsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,Co,nocc,CMO_CABS%elms,ncabsMO,&
        & mynum,numnodesstd,CalphaRcabsMO,NBA,ABdecompR,ABdecompCreateR,intspec,use_bg_buf)

   intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
   intspec(2) = 'R' !Regular AO basis function on center 3
   intspec(3) = 'C' !CABS AO basis function on center 4
   intspec(4) = 'G' !The Gaussian geminal operator g
   intspec(5) = 'G' !The Gaussian geminal operator g
   !Unique: CalphaGcabsMO(NBA,nocc,ncabsMO)
   call Build_CalphaMO2(mylsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,Co,nocc,CMO_CABS%elms,ncabsMO,&
        & mynum,numnodesstd,CalphaGcabsMO,NBA,ABdecompG,ABdecompCreateG,intspec,use_bg_buf)


   IF(DECinfo%F12Ccoupling)THEN 
      call mem_alloc(CalphaTmp,NBA*nocc*nvirt) 
      !CalphaTmp(NBA,nocc,nvirt) = CalphaGcabsMO(NBA,nocc,ncabsMO)*Fac(nvirt,ncabsMO)^T 
      M = NBA*nocc     !Rows of Output Matrix  
      N = nvirt        !Columns of Output Matrix 
      K = ncabsMO      !summation dimension 
      call DGEMM('N','T',M,N,K,1.0E0_realk,CalphaGcabsMO,M,Fac%elms,N,0.0E0_realk,CalphaTmp,M) 
      
      call mem_alloc(EpsOcc,nocc) 
      !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I) & 
      !$OMP SHARED(nocc,MyMolecule,EpsOcc,offset) 
      do I=1,nocc 
         EpsOcc(I) = MyMolecule%oofock%elm2(I+offset,I+offset) 
      enddo
      !$OMP END PARALLEL DO 
      call mem_alloc(EpsVirt,nvirt) 
      !$OMP PARALLEL DO DEFAULT(none) PRIVATE(A) & 
      !$OMP SHARED(nvirt,MyMolecule,EpsVirt) 
      do A=1,nvirt 
         EpsVirt(A) = MyMolecule%vvfock%elm2(A,A) 
      enddo
      !$OMP END PARALLEL DO 

#ifdef VAR_MPI 
      IF(wakeslaves)THEN
         E_21C = 0.0E0_realk
         DO inode = 1,numnodes
            nbuf1 = nAuxMPI(inode)
            NBA2 = nAuxMPI(inode)
            nsize = nbuf1*nocv*nocc
            nsize2 = nbuf1*nvirt*nocc
            IF(mynum.EQ.inode-1)THEN
               !I Bcast My Own CalphaG
               node = mynum            
               IF(size(CalphaG).NE.nsize)call lsquit('MPI Bcast error in Full RIMP2F12 B1',-1)
               call ls_mpibcast(CalphaG,nsize,node,infpar%lg_comm)   !CalphaTmp(NBA,nocv,nocc)
               IF(size(CalphaTmp).NE.nsize2)call lsquit('MPI Bcast error in Full RIMP2F12 B2',-1)
               call ls_mpibcast(CalphaTmp,nsize2,node,infpar%lg_comm) !CalphaTmp(NBA,nocc,nvirt)
               Call FullRIMP2F12_CcouplingEnergyCont(NBA,nocc,nvirt,nbasis,CalphaG,CalphaTmp,E_21Ctmp,EpsOcc,EpsVirt) 
            ELSE
               node = inode-1
               !recieve
               IF(use_bg_buf)THEN
                  call mem_pseudo_alloc(CalphaMPI,nsize)
                  call mem_pseudo_alloc(CalphaMPI2,nsize2)
               ELSE
                  call mem_alloc(CalphaMPI,nsize)
                  call mem_alloc(CalphaMPI2,nsize2)
               ENDIF
               call ls_mpibcast(CalphaMPI,nsize,node,infpar%lg_comm)
               call ls_mpibcast(CalphaMPI2,nsize2,node,infpar%lg_comm)
               call FullRIMP2F12_CcouplingEnergyContMPI(NBA,nocc,nvirt,nbasis,CalphaMPI,&
                    & CalphaMPI2,NBA2,CalphaG,CalphaTmp,E_21Ctmp,EpsOcc,EpsVirt)
               IF(use_bg_buf)THEN
                  call mem_pseudo_dealloc(CalphaMPI2)
                  call mem_pseudo_dealloc(CalphaMPI)
               ELSE
                  call mem_dealloc(CalphaMPI)
                  call mem_dealloc(CalphaMPI2)
               ENDIF
            ENDIF
            E_21C = E_21C + E_21Ctmp
         ENDDO
      ELSE
         Call FullRIMP2F12_CcouplingEnergyCont(NBA,nocc,nvirt,nbasis,CalphaG,CalphaTmp,E_21C,EpsOcc,EpsVirt) 
      ENDIF
      lsmpibufferRIMP2(20)=E_21C      !we need to perform a MPI reduction at the end 
#else
      Call FullRIMP2F12_CcouplingEnergyCont(NBA,nocc,nvirt,nbasis,CalphaG,CalphaTmp,E_21C,EpsOcc,EpsVirt) 
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(CC,RI) = ',E_21C 
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(CC,RI) = ',E_21C 
#endif      
      call mem_dealloc(CalphaTmp) 
      call mem_dealloc(EpsVirt) 
      call mem_dealloc(EpsOcc)      
   ENDIF

#ifdef VAR_MPI 
   IF(wakeslaves)THEN
      EV3 = 0.0E0_realk
      EV4 = 0.0E0_realk
      DO inode = 1,numnodes
         nbuf1 = nAuxMPI(inode)
         NBA2 = nAuxMPI(inode)
         nsize = nbuf1*nocc*ncabsMO  !CalphaRcabsMO(NBA,nocc,ncabsMO)
         nsize2 = nbuf1*nocc*nocv    !CalphaR(NBA,nocc,nocv)
         IF(mynum.EQ.inode-1)THEN
            !I Bcast My Own CalphaG
            node = mynum            
            IF(size(CalphaRcabsMO).NE.nsize)call lsquit('MPI Bcast error in Full RIMP2F12 C1',-1)
            call ls_mpibcast(CalphaRcabsMO,nsize,node,infpar%lg_comm)   
            IF(size(CalphaR).NE.nsize2)call lsquit('MPI Bcast error in Full RIMP2F12 C2',-1)
            call ls_mpibcast(CalphaR,nsize2,node,infpar%lg_comm)        
            call ContractTwo4CenterF12IntegralsRI2V3V4(NBA,NBA,nocc,noccfull,ncabsMO,nocv,&
                 & CalphaRcabsMO,CalphaGcabsMO,CalphaR,CalphaG,EV3tmp,EV4tmp)
         ELSE
            node = inode-1
            !recieve
            IF(use_bg_buf)THEN
               call mem_pseudo_alloc(CalphaMPI,nsize)
               call mem_pseudo_alloc(CalphaMPI2,nsize2)
            ELSE
               call mem_alloc(CalphaMPI,nsize)
               call mem_alloc(CalphaMPI2,nsize2)
            ENDIF
            call ls_mpibcast(CalphaMPI,nsize,node,infpar%lg_comm)
            call ls_mpibcast(CalphaMPI2,nsize2,node,infpar%lg_comm)
            call ContractTwo4CenterF12IntegralsRI2V3V4(NBA,NBA2,nocc,noccfull,ncabsMO,nocv,&
                 & CalphaMPI,CalphaGcabsMO,CalphaMPI2,CalphaG,EV3tmp,EV4tmp)
            IF(use_bg_buf)THEN
               call mem_pseudo_dealloc(CalphaMPI2)
               call mem_pseudo_dealloc(CalphaMPI)
            ELSE
               call mem_dealloc(CalphaMPI)
               call mem_dealloc(CalphaMPI2)
            ENDIF
         ENDIF
         EV3 = EV3 + EV3tmp
         EV4 = EV4 + EV4tmp         
      ENDDO
   ELSE
      call ContractTwo4CenterF12IntegralsRI2V3V4(NBA,NBA,nocc,noccfull,ncabsMO,nocv,&
           & CalphaRcabsMO,CalphaGcabsMO,CalphaR,CalphaG,EV3,EV4)         
   ENDIF
   lsmpibufferRIMP2(8)=EV3      !we need to perform a MPI reduction at the end 
   lsmpibufferRIMP2(9)=EV4      !we need to perform a MPI reduction at the end 
#else
   !Do on GPU (Async)
   call ContractTwo4CenterF12IntegralsRI2V3V4(NBA,NBA,nocc,noccfull,ncabsMO,nocv,&
        & CalphaRcabsMO,CalphaGcabsMO,CalphaR,CalphaG,EV3,EV4)   
   mp2f12_energy = mp2f12_energy  + EV3 + EV4
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V3,RI) = ',EV3       
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V4,RI) = ',EV4       
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V3,RI) = ',EV3       
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V4,RI) = ',EV4       
#endif

   call mem_dealloc(ABdecompR)
   call mem_dealloc(CalphaR)
   call mem_dealloc(CalphaRcabsMO)
   
   !==========================================================
   !=                                                        =
   != V5:     Caibj = (Gcibj*Fac + Gcjai*Fcb)*Taibj          =
   !=                                                        = 
   !========================================================== 
   ABdecompCreateTvirt = .TRUE. 
   call mem_alloc(ABdecompTvirt,naux,naux)
   !FIXME Replace this with subset of CalphaR(NBA,nocc,nbasis)
   intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
   intspec(2) = 'R' !Regular AO basis function on center 3
   intspec(3) = 'R' !CABS AO basis function on center 4
   intspec(4) = 'C' !The Coulomb Operator
   intspec(5) = 'C' !The Coulomb Operator
   !FIXME Replace this with subset of CalphaR(NBA,nocc,nocv)
   call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,MyMolecule%Cv%elm2,nvirt,Co,nocc,&
        & mynum,numnodesstd,CalphaTvirt,NBA,ABdecompTvirt,ABdecompCreateTvirt,intspec,use_bg_buf)
   ABdecompCreateTvirt = .FALSE.
   call mem_dealloc(ABdecompTvirt)
   
   call mem_alloc(ABdecompC,nAux,nAux)
   ABdecompCreateC = .TRUE.

   intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
   intspec(2) = 'R' !Regular AO basis function on center 4 
   intspec(3) = 'C' !Cabs AO basis function on center 3
   intspec(4) = 'G' !The Gaussian geminal operator g
   intspec(5) = 'G' !The Gaussian geminal operator g
   !FIXME Replace this with CalphaGcabsMO(nocc,ncabsMO)
   call Build_CalphaMO2(mylsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,Co,nocc,CMO_CABS%elms,ncabsMO,&
        & mynum,numnodesstd,CalphaCcabs,NBA,ABdecompC,ABdecompCreateC,intspec,use_bg_buf)
   
   m = NBA*nocc       
   k = ncabsMO         ! D_ia = A_ic F_ca
   n = nvirt        
   
   !C(alpha*i,ncabsMO)*F(cabsMO,nvirt)
   nsize = nBA*nocc*ncabsMO
   call mem_alloc(CalphaD, nsize)
   call dgemm('N','T',m,n,k,1.0E0_realk,CalphaCcabs,m,Fac%elms,n,0.0E0_realk,CalphaD,m)
   
   intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
   intspec(2) = 'R' !Regular AO basis function on center 3
   intspec(3) = 'R' !Regular AO basis function on center 4
   intspec(4) = 'G' !The Gaussian geminal operator g
   intspec(5) = 'G' !The Gaussian geminal operator g
   !FIXME Replace this with subset of CalphaG(NBA,nocv,nocc)
   call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,Co,nocc,MyMolecule%Cv%elm2,nvirt,&
        & mynum,numnodesstd,CalphaCvirt,NBA,ABdecompC,ABdecompCreateC,intspec,use_bg_buf)

#ifdef VAR_MPI 
   IF(wakeslaves)THEN
      EV5 = 0.0E0_realk
      DO inode = 1,numnodes
         nbuf1 = nAuxMPI(inode)
         NBA2 = nAuxMPI(inode)
         nsize = nbuf1*nvirt*nocc     !CalphaTvirt(NBA,nvirt,nocc)
         IF(mynum.EQ.inode-1)THEN
            !I Bcast My Own CalphaG
            node = mynum            
            IF(size(CalphaTvirt).NE.nsize)call lsquit('MPI Bcast error in Full RIMP2F12 D',-1)
            call ls_mpibcast(CalphaTvirt,nsize,node,infpar%lg_comm)        
            call ContractTwo4CenterF12IntegralsRIC_pf(MyMolecule,offset,nBA,nocc,nvirt,&
                 & CalphaCvirt,CalphaD,CalphaTvirt,NBA,EV5tmp)
         ELSE
            node = inode-1
            !recieve
            IF(use_bg_buf)THEN
               call mem_pseudo_alloc(CalphaMPI,nsize)
            ELSE
               call mem_alloc(CalphaMPI,nsize)
            ENDIF
            call ls_mpibcast(CalphaMPI,nsize,node,infpar%lg_comm)
            call ContractTwo4CenterF12IntegralsRIC_pf(MyMolecule,offset,nBA,nocc,nvirt,&
                 & CalphaCvirt,CalphaD,CalphaMPI,NBA2,EV5tmp)
            IF(use_bg_buf)THEN
               call mem_pseudo_dealloc(CalphaMPI)
            ELSE
               call mem_dealloc(CalphaMPI)
            ENDIF
         ENDIF
         EV5 = EV5 + EV5tmp
      ENDDO
   ELSE
      call ContractTwo4CenterF12IntegralsRIC_pf(MyMolecule,offset,nBA,nocc,nvirt,CalphaCvirt,&
           & CalphaD,CalphaTvirt,NBA,EV5)
   ENDIF
   lsmpibufferRIMP2(10)=EV5      !we need to perform a MPI reduction at the end 
#else
   call ContractTwo4CenterF12IntegralsRIC_pf(MyMolecule,offset,nBA,nocc,nvirt,CalphaCvirt,&
        & CalphaD,CalphaTvirt,NBA,EV5)
   mp2f12_energy = mp2f12_energy + EV5
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V5,RI) = ', EV5
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V5,RI) = ', EV5
#endif

   call mem_dealloc(CalphaTvirt)                
   !ABdecompCreateG = .FALSE.
   call mem_dealloc(CalphaD)
   call mem_dealloc(ABdecompC)
   call mem_dealloc(CalphaCcabs)
   call mem_dealloc(CalphaCvirt)
   
   !==========================================================
   !=                                                        =
   != X3:         Step 3  Gimjc*Gimjc                        =
   != X4:         Step 4  Gjmic*Gjmic                        =
   != The Coulomb Operator Int multiplied with               =
   != The Gaussian geminal operator g                        =
   != Dim: (nocc,noccfull,nocc,ncabsMO)  need 4 Calphas      =
   !=                                                        =
   !==========================================================
   intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
   intspec(2) = 'R' !Regular AO basis function on center 3
   intspec(3) = 'R' !Regular AO basis function on center 4
   intspec(4) = 'G' !The Gaussian geminal operator g
   intspec(5) = 'G' !The Gaussian geminal operator g
     
   call mem_alloc(ABdecompC,nAux,nAux)
   !FIXME CalphaGocc(NBA,nocc,nocc) is a subset of CalphaG(NBA,nocv,nocc)
   call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,noccfull,Co,nocc,mynum,numnodesstd,CalphaCocc,&
        & NBA,ABdecompC,ABdecompCreateC,intspec,use_bg_buf)
   call mem_dealloc(ABdecompC)

   !Do on GPU (Async) while the CPU starts calculating the next fitting Coef.
   nsize = NBA*noccfull*nocc
   call mem_alloc(CalphaT,nsize)
   M = NBA*noccfull         !rows of Output Matrix
   K = nocc             !summation dimension
   N = nocc             !columns of Output Matrix
   call dgemm('N','N',M,N,K,1.0E0_realk,CalphaCocc,M,Fkj,K,0.0E0_realk,CalphaT,M)
   
   intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
   intspec(2) = 'C' !Regular AO basis function on center 4 
   intspec(3) = 'R' !Cabs AO basis function on center 3
   intspec(4) = 'G' !The Gaussian geminal operator g
   intspec(5) = 'G' !The Gaussian geminal operator g
     
   call mem_alloc(ABdecompC,nAux,nAux)
   !FIXME CalphaCcabsT(NBA,ncabsMO,nocc) is (reordered) CalphaGcabsMO(NBA,nocc,ncabsMO)
   call Build_CalphaMO2(mylsitem,master,ncabsAO,nbasis,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,CMO_CABS%elms,ncabsMO,Co,nocc,&
        & mynum,numnodesstd,CalphaCcabsT,NBA,ABdecompC,ABdecompCreateC,intspec,use_bg_buf)
   call mem_dealloc(ABdecompC)
   
   nsize = NBA*nocc*ncabsMO
   call mem_alloc(CalphaP,nsize)
   M = NBA*ncabsMO         !rows of Output Matrix
   K = nocc                !summation dimension
   N = nocc                !columns of Output Matrix
   call dgemm('N','N',M,N,K,1.0E0_realk,CalphaCcabsT,M,Fkj,K,0.0E0_realk,CalphaP,M)

#ifdef VAR_MPI 
   IF(wakeslaves)THEN
      EX3 = 0.0E0_realk
      EX4 = 0.0E0_realk
      DO inode = 1,numnodes
         nbuf1 = nAuxMPI(inode)
         NBA2 = nAuxMPI(inode)
         nsize = nbuf1*nocc*ncabsMO     !CalphaGcabsMO(NBA,nocc,ncabsMO) 
         nsize2 = nbuf1*noccfull*nocc   !CalphaCocc(NBA,noccfull,nocc) 
         IF(mynum.EQ.inode-1)THEN
            !I Bcast My Own CalphaG
            node = mynum            
            IF(size(CalphaGcabsMO).NE.nsize)call lsquit('MPI Bcast error in Full RIMP2F12 E1',-1)
            call ls_mpibcast(CalphaGcabsMO,nsize,node,infpar%lg_comm)
            IF(size(CalphaCocc).NE.nsize2)call lsquit('MPI Bcast error in Full RIMP2F12 E2',-1)
            call ls_mpibcast(CalphaCocc,nsize2,node,infpar%lg_comm)
            call ContractTwo4CenterF12IntegralsRIX3X4_nc(NBA,nocc,noccfull,ncabsMO,&
                 & CalphaGcabsMO,CalphaCocc,CalphaT,CalphaP,EX3tmp,EX4tmp)
         ELSE
            node = inode-1
            !recieve
            IF(use_bg_buf)THEN
               call mem_pseudo_alloc(CalphaMPI,nsize) 
               call mem_pseudo_alloc(CalphaMPI2,nsize2)
            ELSE
               call mem_alloc(CalphaMPI,nsize)
               call mem_alloc(CalphaMPI2,nsize2)
            ENDIF
            call ls_mpibcast(CalphaMPI,nsize,node,infpar%lg_comm)
            call ls_mpibcast(CalphaMPI2,nsize2,node,infpar%lg_comm)
            call ContractTwo4CenterF12IntegralsRIX3X4_ncMPI(NBA,nocc,noccfull,ncabsMO,&
                 & CalphaMPI,CalphaMPI2,NBA2,&
                 & CalphaGcabsMO,CalphaCocc,CalphaT,CalphaP,EX3tmp,EX4tmp)
            IF(use_bg_buf)THEN
               call mem_pseudo_dealloc(CalphaMPI2)
               call mem_pseudo_dealloc(CalphaMPI)
            ELSE
               call mem_dealloc(CalphaMPI)
               call mem_dealloc(CalphaMPI2)
            ENDIF
         ENDIF
         EX3 = EX3 + EX3tmp
         EX4 = EX4 + EX4tmp         
      ENDDO
   ELSE
      call ContractTwo4CenterF12IntegralsRIX3X4_nc(NBA,nocc,noccfull,ncabsMO,&
           & CalphaGcabsMO,CalphaCocc,CalphaT,CalphaP,EX3,EX4)
   ENDIF
   lsmpibufferRIMP2(11)=EX3      !we need to perform a MPI reduction at the end 
   lsmpibufferRIMP2(12)=EX4      !we need to perform a MPI reduction at the end 
#else
   !Do on GPU (Async) while the CPU starts calculating the next fitting Coef.
   call ContractTwo4CenterF12IntegralsRIX3X4_nc(NBA,nocc,noccfull,ncabsMO,&
        & CalphaGcabsMO,CalphaCocc,CalphaT,CalphaP,EX3,EX4)
   
   mp2f12_energy = mp2f12_energy  + EX3 + EX4
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X3,RI) = ',EX3
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X4,RI) = ',EX4
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X3,RI) = ',EX3
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X4,RI) = ',EX4
#endif

   call mem_dealloc(CalphaCcabsT)
   call mem_dealloc(CalphaCocc)
   call mem_dealloc(CalphaT)
   call mem_dealloc(CalphaP)

   call LSTIMER('FULLRIMP2:Step2',TS2,TE2,DECinfo%output,ForcePrint)

   !=================================================================
   != Step 3: The remaining B terms                                 =
   !=================================================================

   !==============================================================
   !=  B4: (ir|f12|jt)Kst(ir|f12|js)      (r,s,t=CabsAO)         =
   !==============================================================
   intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
   intspec(2) = 'R' !Regular AO basis function on center 3
   intspec(3) = 'C' !Regular AO basis function on center 4
   intspec(4) = 'G' !The f12 Operator
   intspec(5) = 'G' !The f12 Operator   
   !Unique CalphaGcabsAO(NBA,nocc,ncabsAO) 
   call Build_CalphaMO2(mylsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,Co,nocc,CMO_RI%elms,ncabsAO,&
        & mynum,numnodesstd,CalphaGcabsAO,NBA,ABdecompG,ABdecompCreateG,intspec,use_bg_buf)
      
   nsize = nBA*nocc*ncabsAO
   call mem_alloc(CalphaD, nsize)
   m =  nBA*nocc                    ! C_mn = A_mk B_kn
   k =  ncabsAO
   n =  ncabsAO
   
   call dgemm('N','N',m,n,k,1.0E0_realk,CalphaGcabsAO,m,Krr%elms,k,0.0E0_realk,CalphaD,m)

#ifdef VAR_MPI 
   IF(wakeslaves)THEN
      EB4 = 0.0E0_realk
      DO inode = 1,numnodes
         nbuf1 = nAuxMPI(inode)
         NBA2 = nAuxMPI(inode)
         nsize = nbuf1*nocc*ncabsAO     !CalphaGcabsAO(NBA,nocc,ncabsAO) 
         IF(mynum.EQ.inode-1)THEN
            !I Bcast My Own CalphaG
            node = mynum            
            IF(size(CalphaGcabsAO).NE.nsize)call lsquit('MPI Bcast error in Full RIMP2F12 F1',-1)
            call ls_mpibcast(CalphaGcabsAO,nsize,node,infpar%lg_comm)
            call ContractTwo4CenterF12IntegralsRIB4(nBA,nocc,ncabsAO,CalphaGcabsAO,CalphaD,EB4tmp)   
         ELSE
            node = inode-1
            !recieve
            IF(use_bg_buf)THEN
               call mem_pseudo_alloc(CalphaMPI,nsize) 
            ELSE
               call mem_alloc(CalphaMPI,nsize)
            ENDIF
            call ls_mpibcast(CalphaMPI,nsize,node,infpar%lg_comm)
            call ContractTwo4CenterF12IntegralsRIB4MPI(nBA,nocc,ncabsAO,CalphaMPI,NBA2,&
                 & CalphaGcabsAO,CalphaD,EB4tmp)   
            IF(use_bg_buf)THEN
               call mem_pseudo_dealloc(CalphaMPI)
            ELSE
               call mem_dealloc(CalphaMPI)
            ENDIF
         ENDIF
         EB4 = EB4 + EB4tmp         
      ENDDO
   ELSE
      call ContractTwo4CenterF12IntegralsRIB4(nBA,nocc,ncabsAO,CalphaGcabsAO,CalphaD,EB4)   
   ENDIF
   lsmpibufferRIMP2(13)=EB4      !we need to perform a MPI reduction at the end 
#else
   call ContractTwo4CenterF12IntegralsRIB4(nBA,nocc,ncabsAO,CalphaGcabsAO,CalphaD,EB4)   
   mp2f12_energy = mp2f12_energy  + EB4
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B4,RI) = ',EB4
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B4,RI) = ', EB4
#endif
   
   !==============================================================
   !=  B5: (ir|f12|jm)Fsr(si|f12|mj)        (r,s=CabsAO)         =
   !==============================================================   
   !We need CalphaG(NBA,nocc,noccfull) but this is a subset of 
   !CalphaG(NBA,nocc,nbasis) which we already have
   !> Dgemm 
   !nsize = nBA*nocc*ncabsAO
   !IF(size(CalphaD).NE.nsize)call lsquit('dim mismatch CalphaD',-1)
   !m =  nBA*nocc                    ! D_jq = C_jp F_qp
   !k =  ncabsAO
   !n =  ncabsAO
   !Do on GPU (Async)
   !call dgemm('N','N',m,n,k,1.0E0_realk,CalphaGcabsAO,m,Frr%elms,k,0.0E0_realk,CalphaD,m)
   !Do on GPU (Async)
   !call ContractTwo4CenterF12IntegralsRIB5(nBA,nocc,ncabsAO,noccfull,CalphaGcabsAO,CalphaG,CalphaD,EB5)
   !mp2f12_energy = mp2f12_energy  + EB5
   !WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B5,RI) = ',EB5
   !WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B5,RI) = ', EB5
   !call mem_dealloc(CalphaD)


   !We need CalphaG(NBA,nocc,noccfull) but this is a subset of 
   !CalphaG(NBA,nocc,nbasis) which we already have
   !> Dgemm 
   nsize = nBA*nocc*ncabsAO
   IF(size(CalphaD).NE.nsize)call lsquit('dim mismatch CalphaD',-1)
   m =  nBA*nocc                    ! D_jq = C_jp F_qp
   k =  ncabsAO
   n =  ncabsAO

   !Do on GPU (Async)
   call dgemm('N','N',m,n,k,1.0E0_realk,CalphaGcabsAO,m,Frr%elms,k,0.0E0_realk,CalphaD,m)
   !Do on GPU (Async)
#ifdef VAR_MPI 
   IF(wakeslaves)THEN
      EB5 = 0.0E0_realk
      DO inode = 1,numnodes
         nbuf1 = nAuxMPI(inode)
         NBA2 = nAuxMPI(inode)
         nsize = nbuf1*nbasis*nocc   !CalphaG(nBA,nbasis,occ)
         nsize2 = nbuf1*nocc*ncabsAO !CalphaD(nBA,nocc,ncabsAO)
         IF(mynum.EQ.inode-1)THEN
            !I Bcast My Own CalphaG
            node = mynum            
            IF(size(CalphaG).NE.nsize)call lsquit('MPI Bcast error in Full RIMP2F12 G1',-1)
            call ls_mpibcast(CalphaG,nsize,node,infpar%lg_comm)
            IF(size(CalphaD).NE.nsize2)call lsquit('MPI Bcast error in Full RIMP2F12 G2',-1)
            call ls_mpibcast(CalphaD,nsize2,node,infpar%lg_comm)
            call ContractTwo4CenterF12IntegralsRIB5(nBA,nocc,ncabsAO,noccfull,nbasis,&
                 & CalphaGcabsAO,CalphaG,CalphaD,EB5tmp)   
         ELSE
            node = inode-1
            !recieve
            IF(use_bg_buf)THEN
               call mem_pseudo_alloc(CalphaMPI,nsize) 
               call mem_pseudo_alloc(CalphaMPI2,nsize2)
            ELSE
               call mem_alloc(CalphaMPI,nsize)
               call mem_alloc(CalphaMPI2,nsize2)
            ENDIF
            call ls_mpibcast(CalphaMPI,nsize,node,infpar%lg_comm)
            call ls_mpibcast(CalphaMPI2,nsize2,node,infpar%lg_comm)
            call ContractTwo4CenterF12IntegralsRIB5MPI(nBA,nocc,ncabsAO,noccfull,nbasis,&
                 & CalphaMPI,CalphaMPI2,NBA2,&
                 & CalphaGcabsAO,CalphaG,CalphaD,EB5tmp)   
            IF(use_bg_buf)THEN
               call mem_pseudo_dealloc(CalphaMPI2)
               call mem_pseudo_dealloc(CalphaMPI)
            ELSE
               call mem_dealloc(CalphaMPI)
               call mem_dealloc(CalphaMPI2)
            ENDIF
         ENDIF
         EB5 = EB5 + EB5tmp
      ENDDO
   ELSE
      call ContractTwo4CenterF12IntegralsRIB5(nBA,nocc,ncabsAO,noccfull,nbasis,&
           & CalphaGcabsAO,CalphaG,CalphaD,EB5)   
   ENDIF
   lsmpibufferRIMP2(14)=EB5      !we need to perform a MPI reduction at the end 
#else
   call ContractTwo4CenterF12IntegralsRIB5(nBA,nocc,ncabsAO,noccfull,nbasis,CalphaGcabsAO,CalphaG,CalphaD,EB5)   
   mp2f12_energy = mp2f12_energy  + EB5
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B5,RI) = ',EB5
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B5,RI) = ', EB5
#endif
   call mem_dealloc(CalphaD)
   
   !==============================================================
   !=  B6: (ip|f12|ja)Fqp(qi|f12|aj)                             =
   !==============================================================
   !> Dgemm 
   !nsize = nBA*nocc*nbasis
   !call mem_alloc(CalphaD, nsize)
   !m =  nBA*nocc               ! D_jq = C_jp F_pq
   !k =  nbasis
   !n =  nbasis
   !Do on GPU (Async)
   !call dgemm('N','N',m,n,k,1.0E0_realk,CalphaG,m,Fpp%elms,k,0.0E0_realk,CalphaD,m)   
   !Do on GPU (Async)
   !call ContractTwo4CenterF12IntegralsRIB6(nBA,nocc,nvirt,nbasis,CalphaG,CalphaD,EB6)
   !mp2f12_energy = mp2f12_energy  + EB6
   !WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B6,RI) = ',EB6
   !WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B6,RI) = ', EB6
   !call mem_dealloc(CalphaD)

   !> Dgemm 
   intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
   intspec(2) = 'R' !Regular AO basis function on center 1 
   intspec(3) = 'R' !Regular AO basis function on center 2 
   intspec(4) = 'G' !The Gaussian geminal operator g
   intspec(5) = 'G' !The Gaussian geminal operator g

   nsize = NBA*nocc*nocv  
   ! call mem_alloc(ABdecompG,nAux,nAux)
   !FIXME: CalphaP(NBA,nocc,nocv) subst of CalphaG(NBA,nocv,nocc)
   call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
      & FORCEPRINT,wakeslaves,Co,nocc,Cfull,nocv,&
      & mynum,numnodesstd,CalphaP,NBA,ABdecompG,ABdecompCreateG,intspec,use_bg_buf)
   !call mem_dealloc(ABdecompG)

   nsize = nBA*nocc*nocv
   call mem_alloc(CalphaD,nsize)

   m =  nBA*nocc               ! D_jq =P_jp F_pq
   k =  nocv
   n =  nocv
   !Do on GPU (Async)
   call dgemm('N','N',m,n,k,1.0E0_realk,CalphaP,m,Fpp%elms,k,0.0E0_realk,CalphaD,m)   
   call mem_dealloc(CalphaP) 
   !Do on GPU (Async)
#ifdef VAR_MPI 
   IF(wakeslaves)THEN
      EB6 = 0.0E0_realk
      DO inode = 1,numnodes
         nbuf1 = nAuxMPI(inode)
         NBA2 = nAuxMPI(inode)
         nsize = nbuf1*nocv*nocc     !CalphaG(NBA,nocv,nocc) 
         IF(mynum.EQ.inode-1)THEN
            !I Bcast My Own CalphaG
            node = mynum            
            IF(size(CalphaG).NE.nsize)call lsquit('MPI Bcast error in Full RIMP2F12 F1',-1)
            call ls_mpibcast(CalphaG,nsize,node,infpar%lg_comm)
            call ContractTwo4CenterF12IntegralsRIB6(nBA,nocc,nvirt,nocv,noccfull,&
                 & CalphaG,CalphaD,EB6tmp)
         ELSE
            node = inode-1
            !recieve
            IF(use_bg_buf)THEN
               call mem_pseudo_alloc(CalphaMPI,nsize) 
            ELSE
               call mem_alloc(CalphaMPI,nsize)
            ENDIF
            call ls_mpibcast(CalphaMPI,nsize,node,infpar%lg_comm)
            call ContractTwo4CenterF12IntegralsRIB6MPI(nBA,nocc,nvirt,nocv,noccfull,&
                 & CalphaMPI,NBA2,CalphaG,CalphaD,EB6tmp)
            IF(use_bg_buf)THEN
               call mem_pseudo_dealloc(CalphaMPI)
            ELSE
               call mem_dealloc(CalphaMPI)
            ENDIF
         ENDIF
         EB6 = EB6 + EB6tmp         
      ENDDO
   ELSE
      call ContractTwo4CenterF12IntegralsRIB6(nBA,nocc,nvirt,nocv,noccfull,CalphaG,CalphaD,EB6)
   ENDIF
   lsmpibufferRIMP2(15)=EB6      !we need to perform a MPI reduction at the end 
#else
   call ContractTwo4CenterF12IntegralsRIB6(nBA,nocc,nvirt,nocv,noccfull,CalphaG,CalphaD,EB6)
   mp2f12_energy = mp2f12_energy  + EB6
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B6,RI) = ',EB6
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B6,RI) = ', EB6
#endif
   call mem_dealloc(CalphaD)
   
   !==============================================================
   !=  B7: (ic|f12|jm)Fnm(ci|F12|nj)                             =
   !==============================================================
   !   We need CalphaGocc(NBA,nocc,nocc) but this is a subset of CalphaG(NBA,nocc,nbasis)
   intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
   intspec(2) = 'R' !Regular AO basis function on center 3
   intspec(3) = 'R' !Regular AO basis function on center 4
   intspec(4) = 'G' !The Gaussian geminal operator g
   intspec(5) = 'G' !The Gaussian geminal operator g

   call mem_alloc(ABdecompC,nAux,nAux)
   !FIXME: CalphaCoccT(NBA,nocc,noccfull) subst of CalphaG(NBA,nocv,nocc)
   call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
      & FORCEPRINT,wakeslaves,Co,nocc,MyMolecule%Co%elm2,noccfull,mynum,numnodesstd,CalphaCoccT,&
      & NBA,ABdecompC,ABdecompCreateC,intspec,use_bg_buf)
   call mem_dealloc(ABdecompC)

   !> Dgemm 
   nsize = nBA*nocc*noccfull
   call mem_alloc(CalphaD, nsize)
   m =  nBA*nocc                       ! D_jm = Cocc_jn F_nm
   k =  noccfull
   n =  noccfull

   call dgemm('N','N',m,n,k,1.0E0_realk,CalphaCoccT,m,Fmm%elms,k,0.0E0_realk,CalphaD,m)
   call GeneralTwo4CenterF12RICoef1223(nBA,&
        & CalphaD,nocc,noccfull,CalphaGcabsMO,nocc,ncabsMO,&
        & CalphaCoccT,nocc,noccfull,EB7,noccfull,wakeslaves,use_bg_buf,&
        & numnodesstd,nAuxMPI,mynum,F12RIB7,F12RIB7MPI)   
#ifdef VAR_MPI 
   lsmpibufferRIMP2(17)=EB7      !we need to perform a MPI reduction at the end 
#else
   mp2f12_energy = mp2f12_energy  + EB7
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B7,RI) = ',EB7
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B7,RI) = ', EB7
#endif
   call mem_dealloc(CalphaCoccT)
   call mem_dealloc(CalphaD)

   !==============================================================
   !=  B8: (ic|f12|jm)Frm(ci|f12|rj)                             =
   !==============================================================
   !> Dgemm 
   nsize = nBA*nocc*noccfull
   call mem_alloc(CalphaD, nsize)
   m =  nBA*nocc                    ! D_jm = C_jr F_rm
   k =  ncabsAO
   n =  noccfull
   !NB! Changed T to N, dont think it will matter but...
   call dgemm('N','N',m,n,k,1.0E0_realk,CalphaGcabsAO,m,Frm%elms,k,0.0E0_realk,CalphaD,m)
   call GeneralTwo4CenterF12RICoef1223(nBA,CalphaG,nocv,nocc,CalphaGcabsMO,nocc,ncabsMO,&
        & CalphaD,nocc,noccfull,EB8,noccfull,wakeslaves,use_bg_buf,numnodesstd,&
        & nAuxMPI,mynum,F12RIB8,F12RIB8MPI)

#ifdef VAR_MPI 
   lsmpibufferRIMP2(18)=EB8      !we need to perform a MPI reduction at the end 
#else
   mp2f12_energy = mp2f12_energy  + EB8
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B8,RI) = ', EB8
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B8,RI) = ', EB8
#endif
   call mem_dealloc(CalphaGcabsAO)
   call mem_dealloc(CalphaD)

   !==============================================================
   !=  B9: (ip|f12|ja)Fcp(ci|f12|aj)                             =
   !==============================================================
   
   !> Dgemm 
   nsize = nBA*nocc*nocv
   call mem_alloc(CalphaD, nsize)
   m =  nBA*nocc                    ! D_jp = C_jc F_cp
   k =  ncabsMO   
   n =  nocv
   call dgemm('N','N',m,n,k,1.0E0_realk,CalphaGcabsMO,m,Fcp%elms,k,0.0E0_realk,CalphaD,m)
   call mem_dealloc(CalphaGcabsMO)
   call GeneralTwo4CenterF12RICoef1112(nBA,CalphaG,nocv,nocc,CalphaD,nocc,nocv,&
        & EB9,noccfull,wakeslaves,use_bg_buf,numnodesstd,nAuxMPI,mynum,F12RIB9,F12RIB9MPI)
#ifdef VAR_MPI 
   lsmpibufferRIMP2(19)=EB9      !we need to perform a MPI reduction at the end 
#else
   mp2f12_energy = mp2f12_energy  + EB9
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B9,RI) = ', EB9
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B9,RI) = ', EB9
#endif

   call mem_dealloc(CalphaG)
   call mem_dealloc(ABdecompG)
   call mem_dealloc(CalphaD)
   
   call mem_dealloc(Co)
   call mem_dealloc(Fkj)
   
   call mem_dealloc(Cfull)
   call mat_free(CMO_CABS)
   call mat_free(CMO_RI)
   call free_F12_mixed_MO_Matrices(HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp,Fic,Fcd)
   call LSTIMER('FULLRIMP2:Step3',TS2,TE2,DECinfo%output,ForcePrint)
   call LSTIMER('FULLRIMP2F12',TS,TE,DECinfo%output,ForcePrint)

#ifdef VAR_MPI
   nbuf1 = 20
   CALL lsmpi_reduction(lsmpibufferRIMP2,nbuf1,infpar%master,infpar%lg_comm)
   EB1=lsmpibufferRIMP2(1)
   EV1=lsmpibufferRIMP2(2)
   EX1=lsmpibufferRIMP2(3)
   EB2=lsmpibufferRIMP2(4)
   EB3=lsmpibufferRIMP2(5)
   EV2=lsmpibufferRIMP2(6)
   EX2=lsmpibufferRIMP2(7)
   EV3=lsmpibufferRIMP2(8)
   EV4=lsmpibufferRIMP2(9)
   EV5=lsmpibufferRIMP2(10)
   EX3=lsmpibufferRIMP2(11)
   EX4=lsmpibufferRIMP2(12)
   EB4=lsmpibufferRIMP2(13)
   EB5=lsmpibufferRIMP2(14)
   EB6=lsmpibufferRIMP2(15)
   !for some reason missing 16
   EB7=lsmpibufferRIMP2(17)
   EB8=lsmpibufferRIMP2(18)
   EB9=lsmpibufferRIMP2(19)
   E_21C=lsmpibufferRIMP2(20)
   DO I=1,size(lsmpibufferRIMP2)
      mp2f12_energy = mp2f12_energy  + lsmpibufferRIMP2(I)
   ENDDO

   IF(master)THEN
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V1,RI) = ', EV1
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X1,RI) = ', EX1
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B1,RI) = ', EB1
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B2,RI) = ', EB2
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B3,RI) = ', EB3
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V2,RI) = ', EV2
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X2,RI) = ', EX2
      IF(DECinfo%F12Ccoupling)THEN 
         WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(CC,RI) = ', E_21C 
      ENDIF
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V3,RI) = ', EV3
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V4,RI) = ', EV4
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V5,RI) = ', EV5
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X3,RI) = ', EX3
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X4,RI) = ', EX4
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B4,RI) = ', EB4
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B5,RI) = ', EB5
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B6,RI) = ', EB6
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B7,RI) = ', EB7
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B8,RI) = ', EB8
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B9,RI) = ', EB9
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V1,RI) = ', EV1
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X1,RI) = ', EX1
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B1,RI) = ', EB1
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B2,RI) = ', EB2
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B3,RI) = ', EB3
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V2,RI) = ', EV2
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X2,RI) = ', EX2
      IF(DECinfo%F12Ccoupling)THEN 
         WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(CC,RI) = ', E_21C 
      ENDIF
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V3,RI) = ', EV3
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V4,RI) = ', EV4
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V5,RI) = ', EV5
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X3,RI) = ', EX3
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X4,RI) = ', EX4
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B4,RI) = ', EB4
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B5,RI) = ', EB5
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B6,RI) = ', EB6
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B7,RI) = ', EB7
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B8,RI) = ', EB8
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B9,RI) = ', EB9
   ENDIF

   IF(wakeslaves)THEN
      call mem_dealloc(nAuxMPI)
   ENDIF
#endif

    E_21 = 0.0E0_realk
    E_21 = EV1 + EV2 + EV3 + EV4 +EV5 + E_21C

    if(DECinfo%F12debug.AND.master) then
       print *, '----------------------------------------'
       print *, ' E21 V term                             '
       print *, '----------------------------------------'
       write(*,'(1X,a,g25.16)') " E21_CC_term:  ", E_21C
       write(*,'(1X,a,g25.16)') " E21_V_term1:  ", EV1
       write(*,'(1X,a,g25.16)') " E21_V_term2:  ", EV2
       write(*,'(1X,a,g25.16)') " E21_V_term3:  ", EV3
       write(*,'(1X,a,g25.16)') " E21_V_term4:  ", EV4
       write(*,'(1X,a,g25.16)') " E21_V_term5:  ", EV5
       print *, '----------------------------------------'
       write(*,'(1X,a,g25.16)') " E21_Vsum:     ", E_21
       
       write(DECinfo%output,'(1X,a,g25.16)') '----------------------------------------'
       write(DECinfo%output,'(1X,a,g25.16)') ' E21 V term                             '
       write(DECinfo%output,'(1X,a,g25.16)') '----------------------------------------'
       write(DECinfo%output,'(1X,a,g25.16)') " E21_CC_term:  ", E_21C
       write(DECinfo%output,'(1X,a,g25.16)') " E21_V_term1:  ", EV1
       write(DECinfo%output,'(1X,a,g25.16)') " E21_V_term2:  ", EV2
       write(DECinfo%output,'(1X,a,g25.16)') " E21_V_term3:  ", EV3
       write(DECinfo%output,'(1X,a,g25.16)') " E21_V_term4:  ", EV4
       write(DECinfo%output,'(1X,a,g25.16)') " E21_V_term5:  ", EV5
       write(DECinfo%output,'(1X,a,g25.16)') '----------------------------------------'
       write(DECinfo%output,'(1X,a,g25.16)') " E21_Vsum:     ", E_21       
    endif
    E_22 = 0.0E0_realk
    E_22 = EX1 + EX2 + EX3 + EX4 

    if(DECinfo%F12debug.AND.master) then
       print *, '----------------------------------------'
       print *, ' E_22 X term                            '
       print *, '----------------------------------------'
       write(*,'(1X,a,g25.16)') " E22_X_term1: ", EX1
       write(*,'(1X,a,g25.16)') " E22_X_term2: ", EX2
       write(*,'(1X,a,g25.16)') " E22_X_term3: ", EX3
       write(*,'(1X,a,g25.16)') " E22_X_term4: ", EX4
       print *, '----------------------------------------'
       write(*,'(1X,a,g25.16)') " E22_Xsum:    ", E_22

       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) ' E_22 X term                            '
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,'(1X,a,g25.16)') " E22_X_term1: ", EX1
       write(DECinfo%output,'(1X,a,g25.16)') " E22_X_term2: ", EX2
       write(DECinfo%output,'(1X,a,g25.16)') " E22_X_term3: ", EX3
       write(DECinfo%output,'(1X,a,g25.16)') " E22_X_term4: ", EX4
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,'(1X,a,g25.16)') " E22_Xsum: ", E_22
    end if

    ! ***********************************************************
    !   Creating the B matrix 
    ! ***********************************************************

    E_23 = 0.0E0_realk
    E_23 = EB1 + EB2 + EB3 + EB4 + EB5 + EB6 + EB7 + EB8 + EB9

    if(DECinfo%F12debug.AND.master) then
       print *, '----------------------------------------'
       print *, ' E_22 B term                            '
       print *, '----------------------------------------'
       write(*,'(1X,a,g25.16)') " E23_B_term1: ", EB1
       write(*,'(1X,a,g25.16)') " E23_B_term2: ", EB2 
       write(*,'(1X,a,g25.16)') " E23_B_term3: ", EB3   
       write(*,'(1X,a,g25.16)') " E23_B_term4: ", EB4   
       write(*,'(1X,a,g25.16)') " E23_B_term5: ", EB5   
       write(*,'(1X,a,g25.16)') " E23_B_term6: ", EB6  
       write(*,'(1X,a,g25.16)') " E23_B_term7: ", EB7   
       write(*,'(1X,a,g25.16)') " E23_B_term8: ", EB8   
       write(*,'(1X,a,g25.16)') " E23_B_term9: ", EB9  
       print *, '----------------------------------------'
       write(*,'(1X,a,g25.16)') " E23_B_sum:   ", E_23

       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,*) ' E_22 B term                            '
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_term1: ", EB1
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_term2: ", EB2   
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_term3: ", EB3  
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_term4: ", EB4   
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_term5: ", EB5   
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_term6: ", EB6   
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_term7: ", EB7  
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_term8: ", EB8  
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_term9: ", EB9   
       write(DECinfo%output,*) '----------------------------------------'
       write(DECinfo%output,'(1X,a,g25.16)') " E23_B_sum:   ", E_23
    end if

    E_F12 = 0.0E0_realk
    E_F12 = E_21 + E_22 + E_23

    if(DECinfo%F12debug.AND.master) then
       print *,   '----------------------------------------------------------------'
       print *,   '                   DEC-MP2-F12 CALCULATION                      '
       print *,   '----------------------------------------------------------------'
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: MP2 CORRELATION ENERGY (For CC) =  ', MP2_energy
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E21 CORRECTION TO ENERGY =     ', E_21
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E22 CORRECTION TO ENERGY =     ', E_22
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E23 CORRECTION TO ENERGY =     ', E_23
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E22+E23 CORRECTION TO ENERGY = ', E_22+E_23
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 CORRECTION TO ENERGY =         ', E_F12
       print *, '-------------------------------------------------------'
       write(*,'(1X,a,f20.10)') ' WANGY TOYCODE: TOTAL CORRELATION ENERGY (For CC) =', MP2_energy+E_F12
    end if
    if(master)then
       write(DECinfo%output,'(1X,a,f20.10)') '----------------------------------------------------------------'
       write(DECinfo%output,'(1X,a,f20.10)') '                  WANGY DEC-MP2-F12 CALCULATION                 '
       write(DECinfo%output,'(1X,a,f20.10)') '----------------------------------------------------------------'
       write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: MP2 CORRELATION ENERGY (For CC) =  ', MP2_energy
       write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E21 CORRECTION TO ENERGY =     ', E_21
       write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E22 CORRECTION TO ENERGY =     ', E_22
       write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E23 CORRECTION TO ENERGY =     ', E_23
       write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 E22+E23 CORRECTION TO ENERGY = ', E_22+E_23
       write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: F12 CORRECTION TO ENERGY =         ', E_F12
       write(DECinfo%output,'(1X,a,f20.10)') ' WANGY TOYCODE: MP2-F12 CORRELATION ENERGY (CC) =  ', MP2_energy+E_F12
    endif

  end subroutine full_canonical_rimp2_f12

  
  subroutine FullRIMP2F12_CcouplingEnergyCont(NBA,nocc,nvirt,nbasis,Galpha,&
       & Galpha2,E,EpsOcc,EpsVirt)
    implicit none
    real(realk),intent(inout) :: E
    integer,intent(in) :: NBA,nocc,nvirt,nbasis
    real(realk),intent(in) :: Galpha(NBA,nbasis,nocc),Galpha2(NBA,nocc,nvirt)
    real(realk),intent(in) :: EpsOcc(nocc),EpsVirt(nvirt)
    !local variables
    integer :: A,B,J,I,ALPHA
    real(realk) :: TMP,CtmpIAJB,CtmpIBJA,eps,T
    TMP = 0.0E0_realk
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(A,B,J,I,&
    !$OMP ALPHA,CtmpIAJB,CtmpIBJA,eps,T) SHARED(NBA,nocc,nvirt,Galpha,Galpha2,&
    !$OMP EpsOcc,EpsVirt) REDUCTION(+:TMP)
    DO B=1,nvirt
     DO J=1,nocc
      DO A=1,nvirt
       DO I=1,nocc
          eps = EpsOcc(I) + EpsOcc(J) - EpsVirt(A) - EpsVirt(B)
          T = 0.0E0_realk
          DO ALPHA=1,NBA
             T = T  + (Galpha(ALPHA,nocc+A,I)*Galpha2(ALPHA,J,B) + Galpha(ALPHA,nocc+B,J)*Galpha2(ALPHA,I,A))/eps
          ENDDO
          CtmpIAJB = 0.0E0_realk
          DO ALPHA=1,NBA
             CtmpIAJB = CtmpIAJB + Galpha(ALPHA,nocc+A,I)*Galpha2(ALPHA,J,B)+ Galpha(ALPHA,nocc+B,J)*Galpha2(ALPHA,I,A)
          ENDDO
          CtmpIBJA = 0.0E0_realk
          DO ALPHA=1,NBA
             CtmpIBJA = CtmpIBJA + Galpha(ALPHA,nocc+A,J)*Galpha2(ALPHA,I,B) + Galpha(ALPHA,nocc+B,I)*Galpha2(ALPHA,J,A)
          ENDDO
          TMP=TMP+(7.0E0_realk*T*CtmpIAJB + 1.0E0_realk*T*CtmpIBJA)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    E = TMP/32.0E0_realk
  end subroutine FullRIMP2F12_CcouplingEnergyCont

  subroutine FullRIMP2F12_CcouplingEnergyContMPI(NBA,nocc,nvirt,nbasis,GalphaMPI,&
       & Galpha2MPI,NBA2,Galpha,Galpha2,E,EpsOcc,EpsVirt)
    implicit none
    real(realk),intent(inout) :: E
    integer,intent(in) :: NBA,nocc,nvirt,nbasis,NBA2
    real(realk),intent(in) :: GalphaMPI(NBA2,nbasis,nocc),Galpha2MPI(NBA2,nocc,nvirt)
    real(realk),intent(in) :: Galpha(NBA,nbasis,nocc),Galpha2(NBA,nocc,nvirt)
    real(realk),intent(in) :: EpsOcc(nocc),EpsVirt(nvirt)
    !local variables
    integer :: A,B,J,I,ALPHA
    real(realk) :: TMP,CtmpIAJB,CtmpIBJA,eps,T
    TMP = 0.0E0_realk
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(A,B,J,I,&
    !$OMP ALPHA,CtmpIAJB,CtmpIBJA,eps,T) SHARED(NBA,nocc,nvirt,Galpha,Galpha2,&
    !$OMP EpsOcc,EpsVirt,GalphaMPI,Galpha2MPI,NBA2) REDUCTION(+:TMP)
    DO B=1,nvirt
     DO J=1,nocc
      DO A=1,nvirt
       DO I=1,nocc
          eps = EpsOcc(I) + EpsOcc(J) - EpsVirt(A) - EpsVirt(B)
          T = 0.0E0_realk
          DO ALPHA=1,NBA2
             T = T  + (GalphaMPI(ALPHA,nocc+A,I)*Galpha2MPI(ALPHA,J,B) + GalphaMPI(ALPHA,nocc+B,J)*Galpha2MPI(ALPHA,I,A))/eps
          ENDDO
          CtmpIAJB = 0.0E0_realk
          DO ALPHA=1,NBA
             CtmpIAJB = CtmpIAJB + Galpha(ALPHA,nocc+A,I)*Galpha2(ALPHA,J,B)+ Galpha(ALPHA,nocc+B,J)*Galpha2(ALPHA,I,A)
          ENDDO
          CtmpIBJA = 0.0E0_realk
          DO ALPHA=1,NBA
             CtmpIBJA = CtmpIBJA + Galpha(ALPHA,nocc+A,J)*Galpha2(ALPHA,I,B) + Galpha(ALPHA,nocc+B,I)*Galpha2(ALPHA,J,A)
          ENDDO
          TMP=TMP+(7.0E0_realk*T*CtmpIAJB + 1.0E0_realk*T*CtmpIBJA)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    E = TMP/32.0E0_realk
  end subroutine FullRIMP2F12_CcouplingEnergyContMPI

  subroutine lsmpi_matrix_bufcopy(Xmat,master)
    implicit none
    type(matrix) :: Xmat
    logical,intent(in) :: master
#ifdef VAR_MPI
    integer :: nrow,ncol
    real(realk),pointer :: Xfull(:,:)
    IF(master)THEN
       nrow = Xmat%nrow
       ncol = Xmat%ncol
       call mem_alloc(Xfull,nrow,ncol)
       call mat_to_full(Xmat,1.0E0_realk,Xfull)
       CALL ls_mpi_buffer(nrow,infpar%master)
       CALL ls_mpi_buffer(ncol,infpar%master)
       CALL ls_mpi_buffer(Xfull,nrow,ncol,infpar%master)
    ELSE
       CALL ls_mpi_buffer(nrow,infpar%master)
       CALL ls_mpi_buffer(ncol,infpar%master)
       call mat_init(Xmat,nrow,ncol)
       call mem_alloc(Xfull,nrow,ncol)
       CALL ls_mpi_buffer(Xfull,nrow,ncol,infpar%master)           
       call mat_set_from_full(Xfull,1.0E0_realk,Xmat)           
    ENDIF
    call mem_dealloc(Xfull)
#endif
  end subroutine lsmpi_matrix_bufcopy

  !> \brief Calculate F12 terms that are naturally linear scaling
  !> \author Thomas Kjaergaard
  !> \date 2015
  subroutine NaturalLinearScalingF12Terms(MyMolecule,MyLsitem,Dmat)
    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: MyMolecule 
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> HF density matrix
    type(matrix),intent(in) :: Dmat
    !local variables
    integer :: noccfull,nbasis,ncabsAO,nocc,nvirt
    real(realk) :: TS2,TE2,Econt(1),CoulombF12V1,ExchangeF12V1,EV1
    real(realk) :: EB1,CoulombF12B1,ExchangeF12B1,EX1
    real(realk) :: CoulombF12X1,ExchangeF12X1
    type(matrix) :: Ftmp(1),Fcc,Fii
    EV1 = 0.0E0_realk
    EB1 = 0.0E0_realk
    EX1 = 0.0E0_realk
    
    call LSTIMER('START ',TS2,TE2,DECinfo%output)
    IF(DECinfo%NaturalLinearScalingF12TermsV1)THEN
       
       Econt(1) = 0.0E0_realk
       call II_get_CoulombEcont(DECinfo%output,DECinfo%output,mylsitem%setting,&
            & [Dmat],Econt,1,GGemCouOperator)
       CoulombF12V1 = Econt(1) 
       Econt(1) = 0.0E0_realk
       call II_get_exchangeEcont(DECinfo%output,DECinfo%output,mylsitem%setting,&
            & [Dmat],Econt,1,GGemCouOperator)
       ExchangeF12V1 = Econt(1)       
       EV1 = -0.25E0_realk*((5.0E0_realk/2.0E0_realk)*CoulombF12V1+ExchangeF12V1)
       WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'NLSF12 Energy contribution: E(V1,LS) = ',EV1
       WRITE(*,'(A50,F20.13)')'NLSF12 Energy contribution: E(V1,LS) = ', EV1
    ENDIF
    
    !Note X1 can be done in a linear scaling
    !Ftmp(mu,nu) = sum_i( C(mu,i)*F(i,i)*C(nu,i)) 
    !giving Ecoulomb = sum_mu,nu,rho,sigma Ftmp(mu,nu) (mu nu|rho sigma) D(rho,sigma)
    !                + sum_mu,nu,rho,sigma D(mu,nu) (mu nu|rho sigma) Ftmp(rho,sigma)
    IF(DECinfo%NaturalLinearScalingF12TermsX1)THEN
       !The double commutator [[T,g],g] term
       nocc = size(MyMolecule%Co%elm2,2)
       nvirt = size(MyMolecule%Cv%elm2,2)
       noccfull = size(MyMolecule%Co%elm2,2)
       nbasis = Dmat%nrow
       IF(associated(MyMolecule%Fij))THEN
          call mat_init(Ftmp(1),nbasis,nbasis)
          call BuildFtmpRIMP2F12(MyMolecule%Co%elm2,Dmat%ncol,noccfull,&
               & MyMolecule%Fij,Ftmp(1)%elms)
       ELSE
          ! Mixed AO/AO full MO Fock matrix 
          call mat_init(Fcc,nbasis,nbasis)
          call get_AO_Fock(nbasis,nbasis,Fcc,Dmat,MyLsitem,'RRRRC')
          call mat_init(Fii,nocc,nocc)
          call MO_transform_AOMatrix(mylsitem,nbasis,nocc,noccfull,nvirt,&
               & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'ii',Fcc,Fii)
          call mat_free(Fcc)
          call mat_init(Ftmp(1),nbasis,nbasis)
          call BuildFtmpRIMP2F12(MyMolecule%Co%elm2,Dmat%ncol,noccfull,&
               & Fii%elms,Ftmp(1)%elms)
          call mat_free(Fii)
       ENDIF

       Econt(1) = 0.0E0_realk
       call II_get_CoulombEcontF12(DECinfo%output,DECinfo%output,mylsitem%setting,&
            & [Dmat],Ftmp,Econt,1,GGemQuaOperator)
       CoulombF12X1 = Econt(1) 

       Econt(1) = 0.0E0_realk
       call II_get_exchangeEcontF12(DECinfo%output,DECinfo%output,mylsitem%setting,&
            & [Dmat],Ftmp,Econt,1,GGemQuaOperator)
       ExchangeF12X1 = -2.0E0_realk*Econt(1)       

       EX1 = -1.0E0_realk*(0.21875E0_realk*CoulombF12X1 + 0.03125E0_realk*ExchangeF12X1)
       WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'NLSF12 Energy contribution: E(X1,LS) = ',EX1
       WRITE(*,'(A50,F20.13)')'NLSF12 Energy contribution: E(X1,LS) = ', EX1
       call mat_free(Ftmp(1))
    ENDIF
   
    IF(DECinfo%NaturalLinearScalingF12TermsB1)THEN
       Econt(1) = 0.0E0_realk
       call II_get_CoulombEcont(DECinfo%output,DECinfo%output,mylsitem%setting,&
            & [Dmat],Econt,1,GGemGrdOperator)
       CoulombF12B1 = Econt(1) 
       Econt(1) = 0.0E0_realk
       call II_get_exchangeEcont(DECinfo%output,DECinfo%output,mylsitem%setting,&
            & [Dmat],Econt,1,GGemGrdOperator)
       ExchangeF12B1 = Econt(1)       
       EB1 = (1.0E0_realk/32.0E0_realk)*((7.0E0_realk/2.0E0_realk)*CoulombF12B1-ExchangeF12B1)
       WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'NLSF12 Energy contribution: E(B1,LS) = ',EB1
       WRITE(*,'(A50,F20.13)')'NLSF12 Energy contribution: E(B1,LS) = ', EB1
    ENDIF

    MyMolecule%EF12NLSB1 = EB1
    MyMolecule%EF12NLSV1 = EV1
    MyMolecule%EF12NLSX1 = EX1
    call LSTIMER('FULL-LS-F12',TS2,TE2,DECinfo%output)
  end subroutine NaturalLinearScalingF12Terms

  subroutine BuildFtmpRIMP2F12(Co,nbasis,nocc,Fij,Ftmp)
    implicit none
    integer,intent(in)     :: nbasis,nocc
    real(realk),intent(in) :: Fij(nocc,nocc),Co(nbasis,nocc)
    real(realk),intent(inout) :: Ftmp(nbasis,nbasis)
    !local variables
    integer :: I,A,B,J
    real(realk) :: TMP

    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(A,B,I,J,&
    !$OMP TMP) SHARED(Co,nbasis,nocc,Fij,Ftmp)
    DO B=1,nbasis
       DO A=1,nbasis
          TMP = 0.0E0_realk
          DO J=1,nocc
             DO I=1,nocc
                TMP = TMP + Co(A,I)*Fij(I,J)*Co(B,J)
             ENDDO
          ENDDO
          Ftmp(A,B) = TMP 
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  end subroutine BuildFtmpRIMP2F12
#else

  subroutine wangy_dummy_sub12()
      implicit none
  end subroutine wangy_dummy_sub12

#endif

end module fullrimp2f12

#ifdef VAR_MPI
  subroutine full_canonical_rimp2f12_slave
    use fullrimp2f12,only: full_canonical_rimp2_f12, lsmpi_matrix_bufcopy
    use infpar_module !infpar
    use lsmpi_type,only:ls_mpiInitBuffer,ls_mpiFinalizeBuffer,&
         & LSMPIBROADCAST,MPI_COMM_LSDALTON 
    use lsmpi_op,only: mpicopy_lsitem
    use precision
    use typedeftype,only:lsitem
    use lsparameters
    use decmpi_module, only: mpi_bcast_fullmolecule
    use DALTONINFO, only: ls_free
    use matrix_module
    use matrix_operations
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
    real(realk) :: rimp2f12_energy    
    !> The HF density matrix 
    type(matrix) :: Dmat

    ! Init MPI buffer
    ! ***************
    ! Main master:  Prepare for writing to buffer
    ! Local master: Receive buffer
    call ls_mpiInitBuffer(infpar%master,LSMPIBROADCAST,MPI_COMM_LSDALTON)
    ! Integral lsitem
    ! ---------------
    call mpicopy_lsitem(MyLsitem,MPI_COMM_LSDALTON)
    call lsmpi_matrix_bufcopy(Dmat,.false.)    
    call ls_mpiFinalizeBuffer(infpar%master,LSMPIBROADCAST,MPI_COMM_LSDALTON)
    ! Full molecule bcasting
    ! **********************
    call mpi_bcast_fullmolecule(MyMolecule)
    
    ! Finalize MPI buffer
    ! *******************
    ! Main master:  Send stuff to local masters and deallocate temp. buffers
    ! Local master: Deallocate buffer etc.
    call full_canonical_rimp2_f12(MyMolecule,MyLsitem,Dmat,rimp2f12_energy)
    call mat_free(Dmat)
    call ls_free(MyLsitem)
    call molecule_finalize(MyMolecule,.false.)
    
  end subroutine full_canonical_rimp2f12_slave
#endif
