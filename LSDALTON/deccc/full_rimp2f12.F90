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

use dec_fragment_utils
use CABS_operations
use ccintegrals

!WARNING FOR TESTING
!use full_f12contractions

!#ifdef MOD_UNRELEASED
use f12_routines_module
use IntegralInterfaceMOD
use rimp2_module
!#endif 

public :: full_canonical_rimp2_f12

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
   type(fullmolecule), intent(in) :: MyMolecule
   !> Lsitem structure
   type(lsitem), intent(inout) :: mylsitem
   !> HF density matrix
   type(matrix),intent(in) :: Dmat
   !> MP2-F12 correlation energy
   real(realk),intent(inout) :: mp2f12_energy
   !> Canonical MP2 correlation energy
   real(realk) :: mp2_energy
   !local variables
   integer :: nbasis,nocc,nvirt,ncabsAO,ncabsMO
   real(realk) :: E21,Econt(1),E23
   real(realk) :: ExchangeF12V1,CoulombF12V1
   real(realk) :: ExchangeF12X1,CoulombF12X1
   real(realk) :: ExchangeF12B1,CoulombF12B1
   real(realk) :: EV1,EV2,EV3,EV4,EV5,EX1,EX2,EX3,EX4
   real(realk) :: EB1,EB2,EB3,EB4,EB5,EB6,EB7,EB8,EB9
   real(realk) :: TS,TE,TS2,TE2
   integer :: i,j,a,b,p,q,c,m,mynum,numnodes,nAtoms,lupri,nsize
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
   integer :: nAux,NBA,N,K
   real(realk),pointer :: CalphaR(:),CalphaG(:),CalphaF(:),CalphaD(:),CalphaCvirt(:)
   real(realk),pointer :: CalphaRcabsMO(:),CalphaGcabsAO(:),CalphaX(:),CalphaCcabs(:)
   real(realk),pointer :: CalphaGcabsMO(:),CalphaXcabsAO(:)
   real(realk),pointer :: Cfull(:,:),ABdecompR(:,:),ABdecompG(:,:),ABdecompC(:,:)
   real(realk),pointer :: ABdecompF(:,:),Umat(:,:),Rtilde(:,:),ABdecompX(:,:)
   logical :: master,wakeslaves,ABdecompCreateR,ABdecompCreateG,ABdecompCreateF,ABdecompCreateC
   logical :: FORCEPRINT,use_bg_buf,LS,ABdecompCreateX
   character :: intspec(5)
   type(matrix) :: CMO_CABS,CMO_RI
   !========================================================
   ! Additional variables
   !========================================================
   integer :: offset, noccfull
   real(realk),pointer :: Taibj(:,:,:,:) !amplitudes not integrals
   real(realk),pointer :: gmo(:,:,:,:)
   real(realk),pointer :: gao(:,:,:,:)
   real(realk) :: eps
  
   if(MyMolecule%mem_distributed)then
      call lsquit("ERROR(full_canonical_rimp2_f12): does not work with PDM type fullmolecule",-1)
   endif

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
   call LSTIMER('START ',TS2,TE2,DECinfo%output,ForcePrint)

   ! Init stuff
   ! **********
   nbasis = MyMolecule%nbasis
   nvirt  = MyMolecule%nvirt
   naux   = MyMolecule%nauxbasis
   nAtoms = MyMolecule%nAtoms
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

   call determine_CABS_nbast(ncabsAO,ncabsMO,mylsitem%setting,DECinfo%output)

   call mat_init(CMO_CABS,nCabsAO,ncabsMO)
   call build_CABS_MO(CMO_CABS,nCabsAO,mylsitem%SETTING,lupri)    
   call mat_init(CMO_RI,nCabsAO,nCabsAO)
   call build_RI_MO(CMO_RI,nCabsAO,mylsitem%SETTING,lupri)

   !NB Remember to have this!! Else memory leak!
   call free_cabs()

   IF(naux.EQ.0)call lsquit('Error no Aux functions in full_canonical_rimp2_f12',-1)

   noccfull = nocc
   IF(DECinfo%frozencore)call lsquit('RI-MP2-F12 frozen core not implemented',-1)

   call mem_alloc(Cfull,nbasis,nbasis)
   do J=1,nocc
      do I=1,nbasis
         Cfull(I,J) = MyMolecule%Co%elm2(I,J)
      enddo
   enddo
   do P=1,nvirt
      do I=1,nbasis
         Cfull(I,nocc+P) = MyMolecule%Cv%elm2(I,P)
      enddo
   enddo

   call get_F12_mixed_MO_Matrices(MyLsitem,MyMolecule,Dmat,nbasis,ncabsAO,&
      & nocc,noccfull,nvirt,ncabsMO,HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp,Fic,Fcd)

   call LSTIMER('FULLRIMP2:Init',TS2,TE2,DECinfo%output,ForcePrint)
   !=================================================================
   != Step 1:  Fijkl,Xijkl,Dijkl                                    =
   !=          corresponding to V1,X1,B1                            =
   != These are special since                                       =
   != 1. the have a very simple structure(B1 does require robust DF)=
   != 2. They could all be done in a Linear scaling way outside DEC =
   !=    without density fitting. 
   != 3. The intermediates are only used once                       =
   != 4. Due to noccEOS,noccEOS,noccEOS,noccEOS very small mem req  =
   != Note                                                          =
   != Fijkl:                                                        =
   != The Gaussian geminal divided by the Coulomb operator g/r12    =
   != Xijkl                                                         =
   != The Gaussian geminal squared g^2                              =
   != Dijkl                                                         =
   != The double commutator [[T,g],g] with g = Gaussian geminal     =
   != Since the integral (alpha|[[T,g],g]|beta) is not positive     =
   != definite this term is done using robust density fitting       =
   != Done according to Eq. 87 of                                   =
   != J Comput Chem 32: 2492–2513, 2011                             =
   !=================================================================

   write(DECinfo%output,'(/,a)') ' ================================================ '
   write(DECinfo%output,'(a)')   '            FULL-RI-MP2F12 ENERGY TERMS            '
   write(DECinfo%output,'(a,/)') ' ================================================ '
   write(*,'(/,a)') ' ================================================ '
   write(*,'(a)')   '           FULL-RI-MP2F12 ENERGY TERMS             '
   write(*,'(a,/)') ' ================================================ '

   LS = .FALSE.
   if (LS) THEN
      call LSTIMER('START ',TS2,TE2,DECinfo%output,ForcePrint)
      !This is how the V1 Fijkl can be done in a linear scaling manner
      call II_get_CoulombEcont(DECinfo%output,DECinfo%output,mylsitem%setting,&
           & [Dmat],Econt,1,GGemCouOperator)
      CoulombF12V1 = Econt(1) 
      Econt(1) = 0.0E0_realk
      call II_get_exchangeEcont(DECinfo%output,DECinfo%output,mylsitem%setting,&
           & [Dmat],Econt,1,GGemCouOperator)
      ExchangeF12V1 = Econt(1)       
      E21 = -0.25E0_realk*((5.0E0_realk/2.0E0_realk)*CoulombF12V1+ExchangeF12V1)
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V1,LS) = ',E21
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V1,LS) = ', E21

      !Note X1 cannot be done in a linear scaling way due to Fii 

      !This is how the B1 Dijkl can be done in a linear scaling manner
      call II_get_CoulombEcont(DECinfo%output,DECinfo%output,mylsitem%setting,&
           & [Dmat],Econt,1,GGemGrdOperator)
      CoulombF12B1 = Econt(1) 
      Econt(1) = 0.0E0_realk
      call II_get_exchangeEcont(DECinfo%output,DECinfo%output,mylsitem%setting,&
           & [Dmat],Econt,1,GGemGrdOperator)
      ExchangeF12B1 = Econt(1)       
      E21 = -(1.0E0_realk/32.0E0_realk)*((7.0E0_realk/2.0E0_realk)*CoulombF12B1-ExchangeF12B1)
      WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B1,LS) = ',E21
      WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B1,LS) = ', E21
      call LSTIMER('FULLRIMP2:LS',TS2,TE2,DECinfo%output,ForcePrint)
   ENDIF
   call LSTIMER('START ',TS2,TE2,DECinfo%output,ForcePrint)

   !normally I do not like to allocate things at the beginning but 
   !due to an analysis of memory heap performance and the 
   !background buffer features of ordered allocations this is beneficial
   call mem_alloc(ABdecompR,nAux,nAux)
   call mem_alloc(ABdecompF,nAux,nAux)
   call mem_alloc(ABdecompX,nAux,nAux)
   ABdecompCreateF = .TRUE.
   ABdecompCreateX = .TRUE.
   ABdecompCreateR = .TRUE.
   intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
   intspec(2) = 'R' !Regular AO basis function on center 3
   intspec(3) = 'R' !Regular AO basis function on center 4


   ! Calculate the Fitting Coefficients (alpha|F|ij)
   use_bg_buf = .FALSE.
   mp2f12_energy = 0.0E0_realk 
   intspec(4) = 'F' !The Gaussian geminal divided by the Coulomb operator g/r12 (GGemCouOperator)
   intspec(5) = 'F' !The Gaussian geminal divided by the Coulomb operator g/r12 (GGemCouOperator)
   call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,MyMolecule%Co%elm2,nocc,&
        & mynum,numnodes,CalphaF,NBA,ABdecompF,ABdecompCreateF,intspec,use_bg_buf)
   ABdecompCreateF = .FALSE.
   !perform this suborutine on the GPU (async)  - you do not need to wait for the results
   call ContractOne4CenterF12IntegralsRI(NBA,nocc,CalphaF,CoulombF12V1,ExchangeF12V1)

   !Calculate the Fitting Coefficients (alpha|g^2|ij) 
   intspec(4) = '2' !The Gaussian geminal operator g^2 (GGemCouOperator)
   intspec(5) = '2' !The Gaussian geminal operator g^2 (GGemCouOperator)
   call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,MyMolecule%Co%elm2,nocc,&
        & mynum,numnodes,CalphaX,NBA,ABdecompX,ABdecompCreateX,intspec,use_bg_buf)
   ABdecompCreateX = .FALSE.

   !perform this suborutine on the GPU (async)  - you do not need to wait for the results
   call ContractOne4CenterF12IntegralsRI2(NBA,nocc,CalphaX,Fii%elms,CoulombF12X1,ExchangeF12X1)
   
   !Calculate the Fitting Coefficients (alpha|[[T,g],g]|ij) 
   intspec(4) = 'D' !The double commutator [[T,g],g] 
   intspec(5) = 'C' !The metric operator = Regular Coulomb operator 1/r12
   !Build the R coefficient of Eq. 90 of J Comput Chem 32: 2492–2513, 2011
   call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,MyMolecule%Co%elm2,nocc,&
        & mynum,numnodes,CalphaD,NBA,ABdecompR,ABdecompCreateR,intspec,use_bg_buf)
   ABdecompCreateR = .FALSE.
   !We need CalphaR(NBA,nocc,nocc) but this is a subset of the CalphaR(NBA,nocc,nbasis) we need later
   !so we calculate the full CalphaR(NBA,nocc,nbasis)
   intspec(4) = 'C' !Regular Coulomb operator 1/r12
   intspec(5) = 'C' !The metric operator = Regular Coulomb operator 1/r12
   !Build the G coefficient of Eq. 90 of J Comput Chem 32: 2492–2513, 2011
   call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,Cfull,nbasis,&
        & mynum,numnodes,CalphaR,NBA,ABdecompR,ABdecompCreateR,intspec,use_bg_buf)
   !Build the U matrix in Eq. 88 of J Comput Chem 32: 2492–2513, 2011
   call mem_alloc(Umat,nAux,nAux)
   !perform this suborutine on the GPU (Async)
   call Build_RobustERImatU(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,MyMolecule%Co%elm2,nocc,&
        & mynum,numnodes,ABdecompR,'D',Umat)
   !Build the R tilde coefficient of Eq. 89 of J Comput Chem 32: 2492–2513, 2011
   M = NBA          !rows of Output Matrix
   N = nocc*nocc    !columns of Output Matrix
   K = NBA          !summation dimension
   !perform this suborutine on the GPU (Async)
   !note CalphaR is actual of dimensions (NBA,nocc,nbasis) but here we only access
   !the first part (NBA,nocc,nocc) 
   call dgemm('N','N',M,N,K,1.0E0_realk,Umat,M,CalphaR,K,-0.5E0_realk,CalphaD,M)
   call mem_dealloc(Umat)
   !CalphaD is now the R tilde coefficient of Eq. 89 of J Comput Chem 32: 2492–2513, 2011
   !perform this suborutine on the GPU (Async)
   call ContractOne4CenterF12IntegralsRobustRI(nAux,nocc,nbasis,CalphaD,CalphaR,EB1)

   !The minus is due to the Valeev factor
   EV1 = -1.0E0_realk*((5.0E0_realk*0.25E0_realk)*CoulombF12V1-ExchangeF12V1*0.25E0_realk)
   mp2f12_energy = mp2f12_energy  + EV1
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V1,RI) = ',EV1
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V1,RI) = ',EV1
   !minus is due to the overall minus from equation (41) and (42) due to
   !contribution from the \bar{B}_{ij}^{ij}
   EX1 = -1.0E0_realk*(0.21875E0_realk*CoulombF12X1 + 0.03125E0_realk*ExchangeF12X1)
   mp2f12_energy = mp2f12_energy  + EX1
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X1,RI) = ',EX1
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X1,RI) = ',EX1
   mp2f12_energy = mp2f12_energy  + EB1
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B1,RI) = ', EB1
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B1,RI) = ', EB1

   call mem_dealloc(CalphaD)
   call LSTIMER('FULLRIMP2:Step1',TS2,TE2,DECinfo%output,ForcePrint)

   !==============================================================
   !=  B2: sum_c' (ic'|f12^2|jj) hJ_ic' - (jc'|f12^2|ij) hJ_ic'  =
   !=  B3: sum_c' (ii|f12^2|jc') hJ_jc' - (ij|f12^2|ic') hJ_ic'  =
   !==============================================================
   intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
   intspec(2) = 'R' !Regular AO basis function on center 3
   intspec(3) = 'C' !CABS AO basis function on center 4
   intspec(4) = '2' !The f12 Operator
   intspec(5) = '2' !The f12 Operator
   call Build_CalphaMO2(mylsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,CMO_RI%elms,ncabsAO,&
        & mynum,numnodes,CalphaXcabsAO,NBA,ABdecompX,ABdecompCreateX,intspec,use_bg_buf)   
   
   call ContractOne4CenterF12IntegralsRIB23(nBA,nocc,ncabsAO,CalphaXcabsAO,CalphaX,&
        & hJir%elms,1.0E0_realk,EB2,EB3)
   !1.0E0_realk because that term has an overall pluss in Eqs. 25-26
   mp2f12_energy = mp2f12_energy  + EB2
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B2,RI) = ',EB2
   mp2f12_energy = mp2f12_energy  + EB3
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B3,RI) = ',EB3
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B2,RI) = ',EB2
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B3,RI) = ',EB3

   call mem_dealloc(CalphaXcabsAO)
   call mem_dealloc(CalphaX)
   call mem_dealloc(CalphaF)
   call mem_dealloc(ABdecompX)
   call mem_dealloc(ABdecompF)
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
   call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,Cfull,nbasis,&
        & mynum,numnodes,CalphaG,NBA,ABdecompG,ABdecompCreateG,intspec,use_bg_buf)
   ABdecompCreateG = .FALSE.

   !Do on GPU (Async) while the CPU starts calculating the next fitting Coef.
   call ContractTwo4CenterF12IntegralsRI(nBA,nocc,nbasis,CalphaR,CalphaG,EV2)
   mp2f12_energy = mp2f12_energy  + EV2
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V2,RI) = ',EV2       
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V2,RI) = ',EV2

   !==========================================================
   !=                                                        =
   != X2: Gipjq*Gipjq                                        =
   != The Gaussian geminal operator Int multiplied with       =
   != The Gaussian geminal operator g                        =
   != Dim(nocc,nbasis,nocc,nbasis)                           =
   !=                                                        =
   !==========================================================
   !Do on GPU (Async) while the CPU starts calculating the next fitting Coef.
   call ContractTwo4CenterF12IntegralsRIX(nBA,nocc,nbasis,CalphaG,Fii%elms,EX2)
   mp2f12_energy = mp2f12_energy  + EX2
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X2,RI) = ',EX2       
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X2,RI) = ',EX2

   !==========================================================
   !=                                                        =
   != V3:                 Rimjc*Gimjc                        =
   != V4:                 Rjmic*Gjmic                        =
   != The Coulomb Operator Int multiplied with               =
   != The Gaussian geminal operator g                        =
   != Dim: (nocc,noccfull,nocc,ncabsMO)  need 4 Calphas      =
   !=                                                        =
   !==========================================================
   
!   We need CalphaRocc(NBA,nocc,nocc) but this is a subset of CalphaR(NBA,nocc,nbasis)
   intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
   intspec(2) = 'R' !Regular AO basis function on center 3
   intspec(3) = 'C' !CABS AO basis function on center 4
   intspec(4) = 'C' !The Coulomb Operator
   intspec(5) = 'C' !The Coulomb Operator
   call Build_CalphaMO2(mylsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,CMO_CABS%elms,ncabsMO,&
        & mynum,numnodes,CalphaRcabsMO,NBA,ABdecompR,ABdecompCreateR,intspec,use_bg_buf)

!   We need CalphaGocc(NBA,nocc,nocc) but this is a subset of CalphaG(NBA,nocc,nbasis)
   intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
   intspec(2) = 'R' !Regular AO basis function on center 3
   intspec(3) = 'C' !CABS AO basis function on center 4
   intspec(4) = 'G' !The Gaussian geminal operator g
   intspec(5) = 'G' !The Gaussian geminal operator g
   call Build_CalphaMO2(mylsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,CMO_CABS%elms,ncabsMO,&
        & mynum,numnodes,CalphaGcabsMO,NBA,ABdecompG,ABdecompCreateG,intspec,use_bg_buf)
   
     !Do on GPU (Async)
     call ContractTwo4CenterF12IntegralsRI2V3V4(NBA,nocc,noccfull,ncabsMO,nbasis,&
        & CalphaRcabsMO,CalphaGcabsMO,CalphaR,CalphaG,EV3,EV4)

     mp2f12_energy = mp2f12_energy  + EV3 + EV4
     WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V3,RI) = ',EV3       
     WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V4,RI) = ',EV4       
     WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V3,RI) = ',EV3       
     WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V4,RI) = ',EV4       

     call mem_dealloc(ABdecompR)
     call mem_dealloc(CalphaR)
     call mem_dealloc(CalphaRcabsMO)

     !==========================================================
     !=                                                        =
     != V5:     Caibj = (Gcibj*Fac + Gcjai*Fcb)*Taibj          =
     !=                                                        = 
     !=                                                        = 
     !=                                                        =
     !========================================================== 

     ! Get all AO integrals 
     ! ********************
     call mem_alloc(gao,nbasis,nbasis,nbasis,nbasis)
     gao = 0.0E0_realk
     call get_full_AO_integrals(nbasis,ncabsAO,gao,MyLsitem,'RRRRC')
     ! Transform AO integrals to MO integrals (A I | B J)
     call get_4Center_MO_integrals(mylsitem,DECinfo%output,nbasis,nocc,noccfull,nvirt,&
        & MyMolecule%Co%elm2, MyMolecule%Cv%elm2,'aiai',gAO,gMO)
     call mem_dealloc(gao)

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

     call mem_alloc(ABdecompC,nAux,nAux)
     ABdecompCreateC = .TRUE.
     intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
     intspec(2) = 'R' !Regular AO basis function on center 4 
     intspec(3) = 'C' !Cabs AO basis function on center 3
     intspec(4) = 'G' !The Gaussian geminal operator g
     intspec(5) = 'G' !The Gaussian geminal operator g

     call Build_CalphaMO2(mylsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,CMO_CABS%elms,ncabsMO,&
        & mynum,numnodes,CalphaCcabs,NBA,ABdecompC,ABdecompCreateC,intspec,use_bg_buf)

     m = NBA*nocc       
     k = ncabsMO         ! C_mn = A_mk B_kn
     n = nvirt  

     !C(alpha*i,ncabsMO)*F(cabsMO,nvirt)
     nsize = nBA*nocc*ncabsMO
     call mem_alloc(CalphaD, nsize)
     call dgemm('N','T',m,n,k,1.0E0_realk,CalphaCcabs,m,Fac%elms,n,0.0E0_realk,CalphaD,m)

     !m = nvirt      
     !k = ncabsMO         ! C_mn = A_mk B_kn
     !n = nocc*NBA  
     !call dgemm('N','T',m,n,k,1.0E0_realk,Fac%elms,m,CalphaCcabs,k,0.0E0_realk,CalphaD,n)

     intspec(1) = 'D' !Auxuliary DF AO basis function on center 1 (2 empty)
     intspec(2) = 'R' !Regular AO basis function on center 3
     intspec(3) = 'R' !Regular AO basis function on center 4
     intspec(4) = 'G' !The Gaussian geminal operator g
     intspec(5) = 'G' !The Gaussian geminal operator g
     call Build_CalphaMO2(mylsitem,master,nbasis,nbasis,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,MyMolecule%Cv%elm2,nvirt,&
        & mynum,numnodes,CalphaCvirt,NBA,ABdecompC,ABdecompCreateC,intspec,use_bg_buf)

     call ContractTwo4CenterF12IntegralsRIC(nBA,nocc,nvirt,CalphaCvirt,CalphaD,Taibj,EV5)
     mp2f12_energy = mp2f12_energy  + EV5
     WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V5,RI) = ', EV5
     WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(V5,RI) = ', EV5

     ABdecompCreateG = .FALSE.
     call mem_dealloc(CalphaD)
     call mem_dealloc(ABdecompC)
     call mem_dealloc(CalphaCcabs)
     call mem_dealloc(CalphaCvirt)

     !Additional 
     call mem_dealloc(gMO)
     call mem_dealloc(Taibj)

     !==========================================================
     !=                                                        =
     != X3:         Step 3  Gimjc*Gimjc                        =
     != X4:         Step 4  Gjmic*Gjmic                        =
     != The Coulomb Operator Int multiplied with               =
     != The Gaussian geminal operator g                        =
     != Dim: (nocc,noccfull,nocc,ncabsMO)  need 4 Calphas      =
     !=                                                        =
     !==========================================================
     !   We need CalphaGocc(NBA,nocc,nocc) but this is a subset of CalphaG(NBA,nocc,nbasis)

     !Do on GPU (Async) while the CPU starts calculating the next fitting Coef.
     call ContractTwo4CenterF12IntegralsRI2X(NBA,nocc,noccfull,ncabsMO,nbasis,&
        & CalphaGcabsMO,CalphaG,Fii%elms,EX3,EX4)

     mp2f12_energy = mp2f12_energy  + EX3 + EX4
     WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X3,RI) = ',EX3       
     WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X4,RI) = ',EX4       
     WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X3,RI) = ',EX3       
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(X4,RI) = ',EX4       

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
   call Build_CalphaMO2(mylsitem,master,nbasis,ncabsAO,nAux,LUPRI,&
        & FORCEPRINT,wakeslaves,MyMolecule%Co%elm2,nocc,CMO_RI%elms,ncabsAO,&
        & mynum,numnodes,CalphaGcabsAO,NBA,ABdecompG,ABdecompCreateG,intspec,use_bg_buf)
      
   nsize = nBA*nocc*ncabsAO
   call mem_alloc(CalphaD, nsize)
   m =  nBA*nocc                    ! C_mn = A_mk B_kn
   k =  ncabsAO
   n =  ncabsAO
   
   call dgemm('N','N',m,n,k,1.0E0_realk,CalphaGcabsAO,m,Krr%elms,k,0.0E0_realk,CalphaD,m)
   call ContractTwo4CenterF12IntegralsRIB4(nBA,nocc,ncabsAO,CalphaGcabsAO,CalphaD,EB4)
   
   mp2f12_energy = mp2f12_energy  + EB4
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B4,RI) = ',EB4
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B4,RI) = ', EB4
   
   !==============================================================
   !=  B5: (ir|f12|jm)Fsr(si|f12|mj)        (r,s=CabsAO)         =
   !==============================================================   
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
   call ContractTwo4CenterF12IntegralsRIB5(nBA,nocc,ncabsAO,nbasis,CalphaGcabsAO,CalphaG,CalphaD,EB5)
   
   mp2f12_energy = mp2f12_energy  + EB5
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B5,RI) = ',EB5
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B5,RI) = ', EB5
   
   call mem_dealloc(CalphaD)
   
   !==============================================================
   !=  B6: (ip|f12|ja)Fqp(qi|f12|aj)                             =
   !==============================================================
   !> Dgemm 
   nsize = nBA*nocc*nbasis
   call mem_alloc(CalphaD, nsize)
   m =  nBA*nocc               ! D_jq = C_jp F_pq
   k =  nbasis
   n =  nbasis
   
   !Do on GPU (Async)
   call dgemm('N','N',m,n,k,1.0E0_realk,CalphaG,m,Fpp%elms,k,0.0E0_realk,CalphaD,m)   
   !Do on GPU (Async)
   call ContractTwo4CenterF12IntegralsRIB6(nBA,nocc,nvirt,nbasis,CalphaG,CalphaD,EB6)
   
   mp2f12_energy = mp2f12_energy  + EB6
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B6,RI) = ',EB6
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B6,RI) = ', EB6
   
   call mem_dealloc(CalphaD)
   
   !==============================================================
   !=  B7: (ic|f12|jm)Fnm(ci|F12|nj)                             =
   !==============================================================
   !We need CalphaG(NBA,nocc,noccfull) but this is a subset of CalphaG(NBA,nocc,nbasis) 
   !that we already have
   
   !> Dgemm 
   nsize = nBA*nocc*noccfull
   call mem_alloc(CalphaD, nsize)
   m =  nBA*noccfull               ! D_jn = C_jn F_nm
   k =  noccfull
   n =  noccfull
   
   !NB! Changed T to N, dont think it will matter but...
!   call dgemm('N','N',m,n,k,1.0E0_realk,CalphaG,m,Fii%elms,k,0.0E0_realk,CalphaD,m)   
   call ContractOccCalpha(NBA,nocc,noccfull,nbasis,CalphaG,Fii%elms,CalphaD)
   call ContractTwo4CenterF12IntegralsRIB7(nBA,nocc,ncabsMO,nbasis,CalphaGcabsMO,CalphaG,CalphaD,EB7)
   
   mp2f12_energy = mp2f12_energy  + EB7
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B7,RI) = ',EB7
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B7,RI) = ', EB7
   
   call mem_dealloc(CalphaD)
   
   !==============================================================
   !=  B8: (ic|f12|jm)Frm(ci|f12|rj)                             =
   !==============================================================
   
   !> Dgemm 
   nsize = nBA*nocc*nocc
   call mem_alloc(CalphaD, nsize)
   m =  nBA*nocc                    ! D_jm = C_jp F_pm
   n =  nocc   
   k =  ncabsAO
   !NB! Changed T to N, dont think it will matter but...
   call dgemm('N','N',m,n,k,1.0E0_realk,CalphaGcabsAO,m,Frm%elms,k,0.0E0_realk,CalphaD,m)

   !we need CalphaG(NBA,nocc,noccfull) but this is a subset of CalphaG(NBA,nocc,nbasis)
   call ContractTwo4CenterF12IntegralsRIB8(nBA,nocc,ncabsMO,nbasis,CalphaGcabsMO,CalphaG,CalphaD,EB8)
   
   mp2f12_energy = mp2f12_energy  + EB8
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B8,RI) = ', EB8
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B8,RI) = ', EB8
   
   call mem_dealloc(CalphaGcabsAO)
   call mem_dealloc(CalphaD)
   
   !==============================================================
   !=  B9: (ip|f12|ja)Fcp(ci|f12|aj)                             =
   !==============================================================
   !> Dgemm 
   nsize = nBA*nocc*nbasis
   call mem_alloc(CalphaD, nsize)
   m =  nBA*nocc                    ! D_jp = C_jc F_cp
   n =  nbasis
   k =  ncabsMO   
   call dgemm('N','N',m,n,k,1.0E0_realk,CalphaGcabsMO,m,Fcp%elms,k,0.0E0_realk,CalphaD,m)
   call mem_dealloc(CalphaGcabsMO)

   call ContractTwo4CenterF12IntegralsRIB9(nBA,nocc,nvirt,nbasis,CalphaG,CalphaD,EB9)
   
   mp2f12_energy = mp2f12_energy  + EB9
   WRITE(DECINFO%OUTPUT,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B9,RI) = ', EB9
   WRITE(*,'(A50,F20.13)')'RIMP2F12 Energy contribution: E(B9,RI) = ', EB9
   
   call mem_dealloc(CalphaG)
   call mem_dealloc(ABdecompG)
   call mem_dealloc(CalphaD)
   
   call mem_dealloc(Cfull)
   call mat_free(CMO_CABS)
   call mat_free(CMO_RI)
   call free_F12_mixed_MO_Matrices(HJir,Krr,Frr,Fac,Fpp,Fii,Fmm,Frm,Fcp,Fic,Fcd)
   call LSTIMER('FULLRIMP2:Step3',TS2,TE2,DECinfo%output,ForcePrint)
   call LSTIMER('FULLRIMP2F12',TS,TE,DECinfo%output,ForcePrint)

end subroutine full_canonical_rimp2_f12

subroutine ContractOne4CenterF12IntegralsRobustRI(nBA,n,nbasis,Rtilde,CalphaR,EJK)
   implicit none
   integer,intent(in)        :: nBA,n,nbasis
   real(realk),intent(in)    :: Rtilde(nBA,n,n)
   real(realk),intent(in)    :: CalphaR(nBA,n,nbasis)
   real(realk),intent(inout) :: EJK
   !local variables
   integer :: I,ALPHA,J
   real(realk) :: TMP,TMP_IJIJ,TMP_JIIJ,TMPD
   TMPD = 0.0E0_realk
   TMP = 0.0E0_realk
   EJK = 0.0E0_realk
   !  EJK1 = 0.0E0_realk
   !  EJK2 = 0.0E0_realk
   !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I,J,ALPHA,TMP_IJIJ,TMP_JIIJ) SHARED(CalphaR,Rtilde,n,&
   !$OMP nba) REDUCTION(+:TMP,TMPD)
   DO J=1,n
      !Diagonal
      TMP_IJIJ = 0.0E0_realk
      DO ALPHA = 1,NBA
         TMP_IJIJ = TMP_IJIJ + CalphaR(ALPHA,J,J)*Rtilde(ALPHA,J,J) + Rtilde(ALPHA,J,J)*CalphaR(ALPHA,J,J)
      ENDDO
      TMPD = TMPD + TMP_IJIJ
      !Non Diagonal
      DO I=j+1,n
         TMP_IJIJ = 0.0E0_realk
         TMP_JIIJ = 0.0E0_realk
         DO ALPHA = 1,NBA
            TMP_IJIJ = TMP_IJIJ + CalphaR(ALPHA,I,I)*Rtilde(ALPHA,J,J) + Rtilde(ALPHA,I,I)*CalphaR(ALPHA,J,J)
            TMP_JIIJ = TMP_JIIJ + CalphaR(ALPHA,I,J)*Rtilde(ALPHA,J,I) + Rtilde(ALPHA,I,J)*CalphaR(ALPHA,J,I)
         ENDDO
         TMP = TMP + 7.0E0_realk * TMP_IJIJ + TMP_JIIJ
         !        EJK1 = EJK1 + 7.0E0_realk * TMP_IJIJ
         !        EJK2 = EJK2 + TMP_JIIJ
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   EJK = 0.25E0_realk*tmpD + 0.0625E0_realk*TMP
   !  print*,'A',0.25E0_realk*tmpD
   !  print*,'B',0.0625E0_realk*EJK1
   !  print*,'C',0.0625E0_realk*EJK2
end subroutine ContractOne4CenterF12IntegralsRobustRI


!Perform on GPU 
subroutine ContractOne4CenterF12IntegralsRI(nBA,n,Calpha,EJ,EK)
   implicit none
   integer,intent(in)        :: nBA,n
   real(realk),intent(in)    :: Calpha(nBA,n,n)
   real(realk),intent(inout) :: EJ,EK
   !local variables
   integer :: I,ALPHA,J
   real(realk) :: TMP,TMPV(n),TMPI
   !Exchange Fiijj
   EK = 0.0E0_realk
   EJ = 0.0E0_realk
   !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(I,J,&
   !$OMP ALPHA) SHARED(Calpha,n,nba) REDUCTION(+:EK,EJ)
   DO I=1,n
      DO J=1,n
         DO ALPHA = 1,NBA
            EJ = EJ + Calpha(ALPHA,I,I)*Calpha(ALPHA,J,J)
            EK = EK + Calpha(ALPHA,I,J)*Calpha(ALPHA,J,I)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

end subroutine ContractOne4CenterF12IntegralsRI

!Perform on GPU 
subroutine ContractOne4CenterF12IntegralsRI2(nBA,n,Calpha,Fii,EJ,EK)
   implicit none
   integer,intent(in)        :: nBA,n
   real(realk),intent(in)    :: Calpha(nBA,n,n)
   real(realk),intent(inout) :: EJ,EK
   real(realk),intent(IN) :: Fii(n,n)
   !local variables
   integer :: i,alpha,j
   real(realk) :: tmp1,tmp2,TMPV(n),TMPI
   !Exchange Fiijj
   EK = 0.0E0_realk
   EJ = 0.0E0_realk
   DO i=1,n
      DO j=1,n
         tmp1 = 0.0E0_realk
         tmp2 = 0.0E0_realk
         DO alpha = 1, nBA
            tmp1 = tmp1 + Calpha(ALPHA,i,i)*Calpha(ALPHA,j,j)
            tmp2 = tmp2 + Calpha(ALPHA,i,j)*Calpha(ALPHA,j,i)
         ENDDO 
         !print *, "i j i j value: ", i,j,i,j, tmp1
         !print *, "i j j i value: ", i,j,j,i, tmp2

         EJ = EJ + tmp1*(Fii(i,i)+Fii(j,j))
         EK = EK + tmp2*(Fii(i,i)+Fii(j,j))
      ENDDO
   ENDDO

   !print *, "Fij-matrix:"
   !DO i=1,n
   !   DO j=1,n
   !print *, "i j Fii(i,j): ", i,j, Fii(i,j)
   !   ENDDO
   !ENDDO

end subroutine ContractOne4CenterF12IntegralsRI2

!Perform on GPU 
subroutine ContractOne4CenterF12IntegralsRIB23(nBA,n1,n2,CalphaR,CalphaG,HJir,Coeff,EJK2,EJK3)
   implicit none
   integer,intent(in)        :: nBA,n1,n2
   real(realk),intent(in)    :: CalphaR(nBA,n1,n2),CalphaG(nBA,n1,n1)
   real(realk),intent(inout) :: EJK2,EJK3
   real(realk),intent(in)    :: Coeff
   real(realk)               :: ED2,EJ2,EJ3,EK2,EK3
   real(realk),intent(IN)    :: HJir(n1,n2)
   !local variables
   integer :: c,i,j,alpha,beta
   real(realk) :: tmpR2,tmpR21,tmpR22,tmpR31,tmpR32,tmp,tmpR
  
   ED2 = 0.0E0_realk
   EJ2 = 0.0E0_realk
   EK2 = 0.0E0_realk
   EJ3 = 0.0E0_realk
   EK3 = 0.0E0_realk 
   DO c=1,n2
      DO j=1,n1
         !Diagonal
         tmpR2 = 0.0E0_realk
         DO alpha = 1,nBA
            tmpR2 = tmpR2 + CalphaR(alpha,j,c)*CalphaG(alpha,j,j)
         ENDDO
         ED2 = ED2 + tmpR2*hJir(j,c)

         !Non Diagonal
         DO i=j+1,n1
            tmpR21 = 0.0E0_realk
            tmpR22 = 0.0E0_realk

            tmpR31 = 0.0E0_realk
            tmpR32 = 0.0E0_realk
            DO alpha = 1,nBA
               tmpR21 = tmpR21 + CalphaR(alpha,i,c)*CalphaG(alpha,j,j)
               tmpR22 = tmpR22 + CalphaR(alpha,j,c)*CalphaG(alpha,i,j)

               tmpR31 = tmpR31 + CalphaR(alpha,j,c)*CalphaG(alpha,i,i)
               tmpR32 = tmpR32 + CalphaR(alpha,i,c)*CalphaG(alpha,i,j)
            ENDDO
            EJ2 = EJ2 + tmpR21*hJir(i,c)
            EK2 = EK2 + tmpR22*hJir(i,c)

            EJ3 = EJ3 + tmpR31*hJir(j,c)
            EK3 = EK3 + tmpR32*hJir(j,c)
         ENDDO
      ENDDO
   ENDDO

   EJK2 = ED2*1.0/4.0 + (EJ2*(7.0/16.0) + EK2*(1.0/16.0))
   EJK3 = ED2*1.0/4.0 + (EJ3*(7.0/16.0) + EK3*(1.0/16.0)) 

end subroutine ContractOne4CenterF12IntegralsRIB23

subroutine ContractTwo4CenterF12IntegralsRI(nBA,n1,n2,CalphaR,CalphaG,EJK)
   implicit none
   integer,intent(in)        :: nBA,n1,n2
   real(realk),intent(in)    :: CalphaR(nBA,n1,n2),CalphaG(nBA,n1,n2)
   real(realk),intent(inout) :: EJK
   real(realk)               :: ED, EJ,EK
   !local variables
   integer :: q,p,i,j,alpha,beta
   real(realk) :: tmpR,tmpG1,tmpG2,tmp

   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   tmp = 0.0E0_realk
   ED =  0.0E0_realk
   EJ =  0.0E0_realk
   EK =  0.0E0_realk
   EJK = 0.0E0_realk
   !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(i,j,p,q,tmpR,&
   !$OMP tmpG1,tmpG2) SHARED(CalphaR,CalphaG,n2,n1,&
   !$OMP nba) REDUCTION(+:EJ,EK,ED)
   DO q=1,n2
      DO p=1,n2
         DO j=1,n1

            !Diagonal
            tmpR = 0.0E0_realk
            DO alpha = 1,NBA
               tmpR = tmpR + CalphaR(alpha,j,p)*CalphaR(alpha,j,q)
            ENDDO
            tmpG1 = 0.0E0_realk
            DO beta = 1,NBA
               tmpG1 = tmpG1 + CalphaG(beta,j,p)*CalphaG(beta,j,q)
            ENDDO
            ED = ED + tmpR*tmpG1
            !Non Diagonal
            DO i=j+1,n1
               tmpR = 0.0E0_realk
               DO alpha = 1,NBA
                  tmpR = tmpR + CalphaR(alpha,i,p)*CalphaR(alpha,j,q)
               ENDDO
               tmpG1 = 0.0E0_realk
               tmpG2 = 0.0E0_realk
               DO beta = 1,NBA
                  tmpG1 = tmpG1 + CalphaG(beta,i,p)*CalphaG(beta,j,q)
                  tmpG2 = tmpG2 + CalphaG(beta,j,p)*CalphaG(beta,i,q)
               ENDDO
               EJ = EJ + tmpR*tmpG1 
               EK = EK + tmpR*tmpG2
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   EJK = ED + 2.5E0_realk*EJ - 0.50_realk*EK 
end subroutine ContractTwo4CenterF12IntegralsRI

subroutine ContractTwo4CenterF12IntegralsRIC(nBA,n1,n2,CalphaV,CalphaD,Taibj,EJK)
   implicit none 
   integer,intent(in)        :: nBA,n1,n2
   real(realk),intent(in)    :: CalphaV(nBA,n1,n2),CalphaD(nBA,n1,n2)
   real(realk),intent(inout) :: EJK
   real(realk)               :: ED,EJ,EK
   real(realk),pointer       :: Caibj(:,:,:,:)
   real(realk),intent(in)    :: Taibj(:,:,:,:)
   !local variables
   integer :: a,b,i,j,alpha,beta
   real(realk) :: tmp,tmpR,tmpG
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   tmp = 0.0E0_realk
   ED =  0.0E0_realk
   EJ =  0.0E0_realk
   EK =  0.0E0_realk
   EJK = 0.0E0_realk
   DO b=1,n2
      DO a=1,n2
         DO j=1,n1
            !Diagonal
            tmp = 0.0E0_realk
            DO alpha = 1,NBA
               tmp = tmp + CalphaV(alpha,j,a)*CalphaD(alpha,j,b) + CalphaV(alpha,j,b)*CalphaD(alpha,j,a)
            ENDDO
            ED = ED + tmp*Taibj(a,j,b,j)
            !Non Diagonal
            DO i=j+1,n1
               tmpR = 0.0E0_realk
               DO alpha = 1,NBA
                  tmpR = tmpR + CalphaV(alpha,i,a)*CalphaD(alpha,j,b) + CalphaV(alpha,j,b)*CalphaD(alpha,i,a)
               ENDDO
               tmpG = 0.0E0_realk
               DO beta = 1,NBA
                  tmpG = tmpG + CalphaV(beta,j,a)*CalphaD(beta,i,b) + CalphaV(beta,i,b)*CalphaD(beta,j,a)
               ENDDO
               EJ = EJ + tmpR*Taibj(a,i,b,j)
               EK = EK + tmpG*Taibj(a,i,b,j)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   EJK = - ED - 2.5E0_realk*EJ + 0.5E0_realk*EK 
end subroutine ContractTwo4CenterF12IntegralsRIC

subroutine ContractTwo4CenterF12IntegralsRIX(nBA,n1,n2,CalphaG,Fii,EJK)
   implicit none
   integer,intent(in)        :: nBA,n1,n2
   real(realk),intent(in)    :: CalphaG(nBA,n1,n2)
   real(realk),intent(inout) :: EJK
   real(realk)               :: ED, EJ,EK
   real(realk),intent(IN)    :: Fii(n1,n1)
   !local variables
   integer :: q,p,i,j,alpha,beta
   real(realk) :: tmpG1,tmpG2,tmp

   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   tmp = 0.0E0_realk
   ED =  0.0E0_realk
   EJ =  0.0E0_realk
   EK =  0.0E0_realk
   EJK = 0.0E0_realk
   DO q=1,n2
      DO p=1,n2
         DO j=1,n1

            !Diagonal
            tmpG1 = 0.0E0_realk
            DO beta = 1,NBA
               tmpG1 = tmpG1 + CalphaG(beta,j,p)*CalphaG(beta,j,q)
               !We have a factor 2 but it's integrated into the reduction
            ENDDO
            ED = ED + tmpG1*tmpG1*Fii(j,j)
            !Non Diagonal
            DO i=j+1,n1
               tmpG1 = 0.0E0_realk
               tmpG2 = 0.0E0_realk
               DO beta = 1,NBA
                  tmpG1 = tmpG1 + CalphaG(beta,i,p)*CalphaG(beta,j,q)
                  tmpG2 = tmpG2 + CalphaG(beta,j,p)*CalphaG(beta,i,q)
               ENDDO
               EJ = EJ + tmpG1*tmpG1*(Fii(i,i) + Fii(j,j)) 
               EK = EK + tmpG1*tmpG2*(Fii(i,i) + Fii(j,j))
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   EJK = (ED*(0.5E0_realk) + (7.0/16.0*EJ + 1.0/16.0*EK)) 
end subroutine ContractTwo4CenterF12IntegralsRIX

subroutine ContractTwo4CenterF12IntegralsRI2V3V4(nBA,n1,n3,n2,nbasis,&
     & CalphaRcabs,CalphaGcabs,CalphaR,CalphaG,EJK3,EJK4)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,n3,nbasis
   real(realk),intent(in)    :: CalphaRcabs(nBA,n1,n2),CalphaGcabs(nBA,n1,n2)
   real(realk),intent(in)    :: CalphaR(nBA,n1,nbasis),CalphaG(nBA,n1,nbasis)
   real(realk),intent(inout) :: EJK3,EJK4
   real(realk)               :: EJ3, EJ4, EK3, EK4, ED
   !local variables
   integer :: m,c,i,j,alpha,beta
   real(realk) :: tmpR3,tmpG13,tmpG23
   real(realk) :: tmpR4,tmpG14,tmpG24
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   ED =  0.0E0_realk
   EJ3 =  0.0E0_realk
   EK3 =  0.0E0_realk
   EJK3 = 0.0E0_realk
   EJ4 =  0.0E0_realk
   EK4 =  0.0E0_realk
   EJK4 = 0.0E0_realk
   !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(i,j,m,c,tmpR3,tmpR4, &
   !$OMP tmpG13,tmpG23,tmpG14,tmpG24) SHARED(CalphaRcabs,CalphaR,CalphaGcabs,CalphaG,&
   !$OMP n3,n2,n1,nba) REDUCTION(+:EJ3,EK3,EJ4,EK4,ED)
   DO m=1,n3
      DO c=1,n2
         DO j=1,n1
            !Diagonal
            tmpR3 = 0.0E0_realk
            DO alpha = 1,NBA
               tmpR3 = tmpR3 + CalphaR(alpha,j,m)*CalphaRcabs(alpha,j,c)
            ENDDO
            tmpG13 = 0.0E0_realk
            DO beta = 1,NBA
               tmpG13 = tmpG13 + CalphaG(beta,j,m)*CalphaGcabs(beta,j,c)
            ENDDO
            ED = ED + tmpR3*tmpG13
            !Non Diagonal
            DO i=j+1,n1
               tmpR3 = 0.0E0_realk
               tmpR4 = 0.0E0_realk
               DO alpha = 1,NBA
                  tmpR3 = tmpR3 + CalphaR(alpha,i,m)*CalphaRcabs(alpha,j,c)
                  tmpR4 = tmpR4 + CalphaR(alpha,j,m)*CalphaRcabs(alpha,i,c)
               ENDDO              
               tmpG13 = 0.0E0_realk
               tmpG23 = 0.0E0_realk
               tmpG14 = 0.0E0_realk
               tmpG24 = 0.0E0_realk
               DO beta = 1,NBA
                  tmpG13 = tmpG13 + CalphaG(beta,i,m)*CalphaGcabs(beta,j,c)
                  tmpG23 = tmpG23 + CalphaG(beta,j,m)*CalphaGcabs(beta,i,c)

                  ! tmpG14 = tmpG14 + CalphaGocc(beta,j,m)*CalphaG(beta,i,c)
                  ! tmpG24 = tmpG24 + CalphaGocc(beta,i,m)*CalphaG(beta,j,c)
               ENDDO
               EJ3 = EJ3 + tmpR3*tmpG13 
               EK3 = EK3 + tmpR3*tmpG23

               EJ4 = EJ4 + tmpR4*tmpG23
               EK4 = EK4 + tmpR4*tmpG13
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   EJK3 = ED + 2.5E0_realk*EJ3-0.5E0_realk*EK3
   EJK4 = ED + 2.5E0_realk*EJ4-0.5E0_realk*EK4

end subroutine ContractTwo4CenterF12IntegralsRI2V3V4

subroutine ContractTwo4CenterF12IntegralsRI2(nBA,n1,n3,n2,CalphaR,CalphaG,&
      & CalphaRocc,CalphaGocc,EJK3,EJK4)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,n3
   real(realk),intent(in)    :: CalphaR(nBA,n1,n2),CalphaG(nBA,n1,n2)
   real(realk),intent(in)    :: CalphaRocc(nBA,n1,n3),CalphaGocc(nBA,n1,n3)
   real(realk),intent(inout) :: EJK3,EJK4
   real(realk)               :: EJ3, EJ4, EK3, EK4, ED
   !local variables
   integer :: m,c,i,j,alpha,beta
   real(realk) :: tmpR3,tmpG13,tmpG23
   real(realk) :: tmpR4,tmpG14,tmpG24
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   ED =  0.0E0_realk
   EJ3 =  0.0E0_realk
   EK3 =  0.0E0_realk
   EJK3 = 0.0E0_realk
   EJ4 =  0.0E0_realk
   EK4 =  0.0E0_realk
   EJK4 = 0.0E0_realk
   !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(i,j,m,c,tmpR3,tmpR4, &
   !$OMP tmpG13,tmpG23,tmpG14,tmpG24) SHARED(CalphaR,CalphaRocc,CalphaG,CalphaGocc,n3,n2,n1,&
   !$OMP nba) REDUCTION(+:EJ3,EK3,EJ4,EK4,ED)
   DO m=1,n3
      DO c=1,n2
         DO j=1,n1
            !Diagonal
            tmpR3 = 0.0E0_realk
            DO alpha = 1,NBA
               tmpR3 = tmpR3 + CalphaRocc(alpha,j,m)*CalphaR(alpha,j,c)
            ENDDO
            tmpG13 = 0.0E0_realk
            DO beta = 1,NBA
               tmpG13 = tmpG13 + CalphaGocc(beta,j,m)*CalphaG(beta,j,c)
            ENDDO
            ED = ED + tmpR3*tmpG13
            !Non Diagonal
            DO i=j+1,n1
               tmpR3 = 0.0E0_realk
               tmpR4 = 0.0E0_realk
               DO alpha = 1,NBA
                  tmpR3 = tmpR3 + CalphaRocc(alpha,i,m)*CalphaR(alpha,j,c)
                  tmpR4 = tmpR4 + CalphaRocc(alpha,j,m)*CalphaR(alpha,i,c)
               ENDDO              
               tmpG13 = 0.0E0_realk
               tmpG23 = 0.0E0_realk
               tmpG14 = 0.0E0_realk
               tmpG24 = 0.0E0_realk
               DO beta = 1,NBA
                  tmpG13 = tmpG13 + CalphaGocc(beta,i,m)*CalphaG(beta,j,c)
                  tmpG23 = tmpG23 + CalphaGocc(beta,j,m)*CalphaG(beta,i,c)

                  ! tmpG14 = tmpG14 + CalphaGocc(beta,j,m)*CalphaG(beta,i,c)
                  ! tmpG24 = tmpG24 + CalphaGocc(beta,i,m)*CalphaG(beta,j,c)
               ENDDO
               EJ3 = EJ3 + tmpR3*tmpG13 
               EK3 = EK3 + tmpR3*tmpG23

               EJ4 = EJ4 + tmpR4*tmpG23
               EK4 = EK4 + tmpR4*tmpG13
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   EJK3 = ED + 2.5E0_realk*EJ3-0.5E0_realk*EK3
   EJK4 = ED + 2.5E0_realk*EJ4-0.5E0_realk*EK4

end subroutine ContractTwo4CenterF12IntegralsRI2


subroutine ContractTwo4CenterF12IntegralsRIB4(nBA,n1,n2,CalphaG,CalphaD,EJK)
   implicit none
   integer,intent(in)        :: nBA,n1,n2
   real(realk),intent(in)    :: CalphaG(nBA,n1,n2)
   real(realk),intent(in)    :: CalphaD(nBa,n1,n2)
   real(realk),intent(inout) :: EJK
   real(realk)               :: EJ, EK, ED
   !local variables
   integer :: p,q,r,i,j,alpha,beta,alpha1,beta1,alpha2,beta2,alpha3,beta3,alpha4,beta4

   real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
   real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2

   ED = 0.0E0_realk
   EJ = 0.0E0_realk
   EK = 0.0E0_realk
  
   DO q=1,n2 !ncabsAO
      DO r=1,n2 !ncabsAO
         DO j=1,n1 !nocc
            
            !Diagonal
            tmpR = 0.0E0_realk
            tmpG = 0.0E0_realk
            DO beta = 1,nBA
               tmpR = tmpR + CalphaG(beta,j,r)*CalphaD(beta,j,q)
               tmpG = tmpG + CalphaG(beta,j,r)*CalphaG(beta,j,q)
            ENDDO
            
            ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
            
            !Non Diagonal
            DO i=j+1,n1
               tmpRJ1 = 0.0E0_realk
               tmpGJ1 = 0.0E0_realk
               DO alpha1 = 1, nBA
                  tmpRJ1 = tmpRJ1 + CalphaG(alpha1,i,r)*CalphaD(alpha1,j,q) 
                  tmpGJ1 = tmpGJ1 + CalphaG(alpha1,i,r)*CalphaG(alpha1,j,q)
               ENDDO
               tmpRJ2 = 0.0E0_realk
               tmpGJ2 = 0.0E0_realk
               DO alpha2 = 1, nBA
                  tmpRJ2 = tmpRJ2 + CalphaG(alpha2,j,r)*CalphaD(alpha2,i,q)
                  tmpGJ2 = tmpGJ2 + CalphaG(alpha2,j,r)*CalphaG(alpha2,i,q)
               ENDDO
               EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
               EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
            ENDDO
         ENDDO
      ENDDO
   ENDDO   
   EJK = ED*0.5E0_realk + 7.0/16.0*EJ + 1.0/16.0*EK 
end subroutine ContractTwo4CenterF12IntegralsRIB4                

subroutine ContractTwo4CenterF12IntegralsRIB5(nBA,n1,n2,nbasis,CalphaGcabs,CalphaG,CalphaD,EJK)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,nbasis
   real(realk),intent(in)    :: CalphaG(nBA,n1,nbasis)
   real(realk),intent(in)    :: CalphaGcabs(nBA,n1,n2)
   real(realk),intent(in)    :: CalphaD(nBA,n1,n2)
   real(realk),intent(inout) :: EJK
   real(realk)               :: EJ, EK, ED
   !local variables
   integer :: p,q,m,i,j,alpha,beta,alpha1,beta1,alpha2,beta2

   real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
   real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2

   ED = 0.0E0_realk
   EJ = 0.0E0_realk
   EK = 0.0E0_realk

   DO q=1,n2 !ncabsAO
      DO m=1,n1 !nocc
         DO j=1,n1 !nocc
            
            !Diagonal
            tmpR = 0.0E0_realk
            tmpG = 0.0E0_realk
            
            DO alpha = 1,nBA
               tmpR = tmpR + CalphaD(alpha,j,q)*CalphaG(alpha,j,m)
            ENDDO
            DO beta = 1,nBA
               tmpG = tmpG + CalphaGcabs(beta,j,q)*CalphaG(beta,j,m)
            ENDDO
            
            ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
            
            !Non Diagonal
            DO i=j+1,n1
               tmpRJ1 = 0.0E0_realk
               tmpGJ1 = 0.0E0_realk
               DO alpha1 = 1, nBA
                  tmpRJ1 = tmpRJ1 + CalphaD(alpha1,i,q)*CalphaG(alpha1,j,m) 
                  tmpGJ1 = tmpGJ1 + CalphaGcabs(alpha1,i,q)*CalphaG(alpha1,j,m)
               ENDDO
               tmpRJ2 = 0.0E0_realk
               tmpGJ2 = 0.0E0_realk
               DO alpha2 = 1, nBA
                  tmpRJ2 = tmpRJ2 + CalphaD(alpha2,j,q)*CalphaG(alpha2,i,m)
                  tmpGJ2 = tmpGJ2 + CalphaGcabs(alpha2,j,q)*CalphaG(alpha2,i,m)
               ENDDO
               EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
               EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   EJK = -1.0E0_realk*(ED*0.5E0_realk + 7.0/16.0*EJ + 1.0/16.0*EK) 
end subroutine ContractTwo4CenterF12IntegralsRIB5                

subroutine ContractTwo4CenterF12IntegralsRIB6(nBA,n1,n2,n3,CalphaG,CalphaD,EJK)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,n3
   real(realk),intent(in)    :: CalphaG(nBA,n1,n3)
   real(realk),intent(in)    :: CalphaD(nBA,n1,n3)
   real(realk),intent(inout) :: EJK
   real(realk)               :: EJ, EK, ED
   !local variables
   integer :: q,a,i,j,alpha,beta,alpha1,beta1,alpha2,beta2

   real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
   real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2

   ED = 0.0E0_realk
   EK = 0.0E0_realk
   EJ = 0.0E0_realk

   DO q=1,n3 !nbasis
      DO a=n1+1,n3 !nvirt
         DO j=1,n1 !nocc
            
            !Diagonal
            tmpR = 0.0E0_realk
            tmpG = 0.0E0_realk
            
            DO alpha = 1,nBA
               tmpR = tmpR + CalphaD(alpha,j,q)*CalphaG(alpha,j,a)
               tmpG = tmpG + CalphaG(alpha,j,q)*CalphaG(alpha,j,a)
            ENDDO
            ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
            !Non Diagonal
            DO i=j+1,n1
               tmpRJ1 = 0.0E0_realk
               tmpGJ1 = 0.0E0_realk
               DO alpha1 = 1, nBA
                  tmpRJ1 = tmpRJ1 + CalphaD(alpha1,i,q)*CalphaG(alpha1,j,a) 
                  tmpGJ1 = tmpGJ1 + CalphaG(alpha1,i,q)*CalphaG(alpha1,j,a)
               ENDDO
               tmpRJ2 = 0.0E0_realk
               tmpGJ2 = 0.0E0_realk
               DO alpha2 = 1, nBA
                  tmpRJ2 = tmpRJ2 + CalphaD(alpha2,j,q)*CalphaG(alpha2,i,a)
                  tmpGJ2 = tmpGJ2 + CalphaG(alpha2,j,q)*CalphaG(alpha2,i,a)
               ENDDO
               EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
               EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)           
            ENDDO
         ENDDO
      ENDDO
   ENDDO   
   EJK = -1.0E0_realk*(ED*0.5E0_realk + 7.0/16.0*EJ + 1.0/16.0*EK) 
end subroutine ContractTwo4CenterF12IntegralsRIB6
            
!warning I am very unclear about the Frozen core implementation of this 
subroutine ContractOccCalpha(NBA,nocc,noccfull,nbasis,CalphaG,Fii,CalphaD)
  implicit none
  integer,intent(in) :: NBA,nocc,noccfull,nbasis
  real(realk),intent(in) :: CalphaG(NBA,nocc,nbasis),Fii(noccfull,noccfull)
  real(realk),intent(inout) :: CalphaD(NBA,nocc,noccfull)
  !
  integer :: i,j,alpha,m
  real(realk) :: TMP
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,alpha,m,&
  !$OMP TMP) SHARED(NBA,nocc,noccfull,CalphaG,CalphaD,Fii)
  !$OMP DO COLLAPSE(3)
  DO m=1,noccfull
     DO I=1,nocc
        DO ALPHA=1,NBA
           CalphaD(ALPHA,I,m) = CalphaG(ALPHA,I,1)*Fii(1,M)
        ENDDO
     ENDDO
  ENDDO
  !$OMP END DO
  !$OMP DO COLLAPSE(2) 
  DO m=1,noccfull
     DO I=1,nocc
        DO J=2,noccfull
           TMP = Fii(J,M)
           DO ALPHA=1,NBA
              CalphaD(ALPHA,I,m) = CalphaD(ALPHA,I,m) + CalphaG(ALPHA,I,J)*TMP
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
end subroutine ContractOccCalpha

subroutine ContractTwo4CenterF12IntegralsRIB7(nBA,n1,n2,nbasis,CalphaR,CalphaG,CalphaD,EJK)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,nbasis
   real(realk),intent(in)    :: CalphaG(nBA,n1,nbasis)
   real(realk),intent(in)    :: CalphaR(nBA,n1,n2)
   real(realk),intent(in)    :: CalphaD(nBA,n1,n1)
   real(realk),intent(inout) :: EJK
   real(realk)               :: EJ, EK, ED
   !local variables
   integer :: p,n,i,j,alpha,beta,alpha1,beta1,alpha2,beta2

   real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
   real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2

   ED = 0.0E0_realk
   EK = 0.0E0_realk
   EJ = 0.0E0_realk

   DO p=1,n2 !ncabs
      DO n=1,n1 !nocc
         DO j=1,n1 !nocc
            
            !Diagonal
            tmpR = 0.0E0_realk
            tmpG = 0.0E0_realk
            
            DO alpha = 1,nBA
               tmpR = tmpR + CalphaR(alpha,j,p)*CalphaD(alpha,j,n)
            ENDDO
            DO beta = 1,nBA
               tmpG = tmpG + CalphaR(beta,j,p)*CalphaG(beta,j,n)
            ENDDO
            
            ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
            
            !Non Diagonal
            DO i=j+1,n1
               tmpRJ1 = 0.0E0_realk
               DO alpha1 = 1, nBA
                  tmpRJ1 = tmpRJ1 + CalphaR(alpha1,i,p)*CalphaD(alpha1,j,n) 
               ENDDO
               tmpRJ2 = 0.0E0_realk
               DO alpha2 = 1, nBA
                  tmpRJ2 = tmpRJ2 + CalphaR(alpha2,j,p)*CalphaD(alpha2,i,n)
               ENDDO
               tmpGJ1 = 0.0E0_realk
               DO beta1 = 1, nBA
                  tmpGJ1 = tmpGJ1 + CalphaR(beta1,i,p)*CalphaG(beta1,j,n)
               ENDDO
               tmpGJ2 = 0.0E0_realk
               DO beta2 = 1, nBA
                  tmpGJ2 = tmpGJ2 + CalphaR(beta2,j,p)*CalphaG(beta2,i,n)
               ENDDO
               EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
               EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
               
            ENDDO
         ENDDO
      ENDDO
   ENDDO   
   EJK = ED*0.5E0_realk + 7.0_realk/16.0_realk*EJ + 1.0_realk/16.0_realk*EK

end subroutine ContractTwo4CenterF12IntegralsRIB7

subroutine ContractTwo4CenterF12IntegralsRIB8(nBA,n1,n2,nbasis,CalphaR,CalphaG,CalphaD,EJK)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,nbasis
   real(realk),intent(in)    :: CalphaG(nBA,n1,nbasis)
   real(realk),intent(in)    :: CalphaR(nBA,n1,n2) !CalphaGcabsMO
   real(realk),intent(in)    :: CalphaD(nBA,n1,n1)
   real(realk),intent(inout) :: EJK
   real(realk)               :: EJ, EK, ED
   !local variables
   integer :: p,m,i,j,alpha,beta,alpha1,beta1,alpha2,beta2

   real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
   real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2

   ED = 0.0E0_realk
   EK = 0.0E0_realk
   EJ = 0.0E0_realk
   
   !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(p,m,i,j,alpha,beta,alpha1,&
   !$OMP beta1,alpha2,beta2,tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2,tmpG,tmpGJ1,tmpGJ2,&
   !$OMP tmpGK1,tmpGK2) SHARED(CalphaR,CalphaG,CalphaD,n2,n1,nba) REDUCTION(+:ED,EJ,EK)
   DO j=1,n1 !nocc
      DO p=1,n2 !ncabs
         DO m=1,n1 !noccfull
            
            !Diagonal
            tmpR = 0.0E0_realk
            tmpG = 0.0E0_realk
            
            DO alpha = 1,nBA
               tmpR = tmpR + CalphaR(alpha,j,p)*CalphaG(alpha,j,m)
            ENDDO
            DO beta = 1,nBA
               tmpG = tmpG + CalphaR(beta,j,p)*CalphaD(beta,j,m)
            ENDDO
            
            ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
            
            !Non Diagonal
            DO i=j+1,n1
               tmpRJ1 = 0.0E0_realk
               DO alpha1 = 1, nBA
                  tmpRJ1 = tmpRJ1 + CalphaR(alpha1,i,p)*CalphaG(alpha1,j,m) 
               ENDDO
               tmpRJ2 = 0.0E0_realk
               DO alpha2 = 1, nBA
                  tmpRJ2 = tmpRJ2 + CalphaR(alpha2,j,p)*CalphaG(alpha2,i,m)
               ENDDO
               tmpGJ1 = 0.0E0_realk
               DO beta1 = 1, nBA
                  tmpGJ1 = tmpGJ1 + CalphaR(beta1,i,p)*CalphaD(beta1,j,m)
               ENDDO
               tmpGJ2 = 0.0E0_realk
               DO beta2 = 1, nBA
                  tmpGJ2 = tmpGJ2 + CalphaR(beta2,j,p)*CalphaD(beta2,i,m)
               ENDDO
               EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
               EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   
   EJK = -2.0E0_realk*(ED*0.5E0_realk + 7.0_realk/16.0_realk*EJ + 1.0_realk/16.0_realk*EK)

end subroutine ContractTwo4CenterF12IntegralsRIB8

subroutine ContractTwo4CenterF12IntegralsRIB9(nBA,n1,n2,n3,CalphaG,CalphaD,EJK)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,n3
   real(realk),intent(in)    :: CalphaG(nBA,n1,n3)
   real(realk),intent(in)    :: CalphaD(nBA,n1,n3)
   real(realk),intent(inout) :: EJK
   real(realk)               :: EJ, EK, ED
   !local variables
   integer :: p,a,i,j,alpha,beta,alpha1,beta1,alpha2,beta2

   real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
   real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2

   ED = 0.0E0_realk
   EK = 0.0E0_realk
   EJ = 0.0E0_realk

   DO p=1,n3 !ncabs
      DO a=n1+1,n3 !nvirt
         DO j=1,n1 !nocc

               !Diagonal
               tmpR = 0.0E0_realk
               DO alpha = 1,nBA
                  tmpR = tmpR + CalphaG(alpha,j,p)*CalphaG(alpha,j,a)
               ENDDO
               tmpG = 0.0E0_realk
               DO beta = 1,nBA
                  tmpG = tmpG + CalphaD(beta,j,p)*CalphaG(beta,j,a)
               ENDDO

               ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 

               !Non Diagonal
               DO i=j+1,n1
                  tmpRJ1 = 0.0E0_realk
                  DO alpha1 = 1, nBA
                     tmpRJ1 = tmpRJ1 + CalphaG(alpha1,i,p)*CalphaG(alpha1,j,a) 
                  ENDDO      
                  tmpRJ2 = 0.0E0_realk
                  DO alpha2 = 1, nBA
                     tmpRJ2 = tmpRJ2 + CalphaG(alpha2,j,p)*CalphaG(alpha2,i,a)
                  ENDDO
                  tmpGJ1 = 0.0E0_realk
                  DO beta1 = 1, nBA
                     tmpGJ1 = tmpGJ1 + CalphaD(beta1,i,p)*CalphaG(beta1,j,a)
                  ENDDO
                  tmpGJ2 = 0.0E0_realk
                  DO beta2 = 1, nBA
                     tmpGJ2 = tmpGJ2 + CalphaD(beta2,j,p)*CalphaG(beta2,i,a)
                 ENDDO
                  EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                  EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   
   EJK = -2.0E0_realk*(ED*0.5E0_realk + 7.0_realk/16.0_realk*EJ + 1.0_realk/16.0_realk*EK)

end subroutine ContractTwo4CenterF12IntegralsRIB9

subroutine ContractTwo4CenterF12IntegralsRI2X(nBA,n1,n3,n2,nbasis,&
     & CalphaGcabs,CalphaG,Fii,EJK3,EJK4)
  implicit none
  integer,intent(in)        :: nBA,n1,n2,n3,nbasis
  real(realk),intent(in)    :: CalphaGcabs(nBA,n1,n2)
  real(realk),intent(in)    :: CalphaG(nBA,n1,nbasis)
  real(realk),intent(IN)    :: Fii(n1,n1)
  real(realk),intent(inout) :: EJK3,EJK4
  real(realk)               :: EJ3, EJ4, EK3, EK4, ED
  !local variables
  integer :: m,c,i,j,alpha,beta
  real(realk) :: tmpR3,tmpG13,tmpG23
  real(realk) :: tmpR4,tmpG14,tmpG24,tmp
  !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
  ED =  0.0E0_realk
  EJ3 =  0.0E0_realk
  EK3 =  0.0E0_realk
  EJK3 = 0.0E0_realk
  EJ4 =  0.0E0_realk
  EK4 =  0.0E0_realk
  EJK4 = 0.0E0_realk
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(i,j,m,c,tmpR3,tmpR4, &
  !$OMP tmpG13,tmpG23,tmp) SHARED(CalphaG,CalphaGcabs,n3,n2,n1,&
  !$OMP nba, Fii) REDUCTION(+:EJ3,EK3,EJ4,EK4,ED)
  DO m=1,n3
     DO c=1,n2
        DO j=1,n1
           !Diagonal
           tmpG13 = 0.0E0_realk
           DO beta = 1,NBA
              tmpG13 = tmpG13 + CalphaG(beta,j,m)*CalphaGcabs(beta,j,c)
           ENDDO
           ED = ED + tmpG13*tmpG13*Fii(j,j)
           !Non Diagonal
           DO i=j+1,n1
!              tmpR3 = 0.0E0_realk
!              tmpR4 = 0.0E0_realk
!              DO alpha = 1,NBA
!                 tmpR3 = tmpR3 + CalphaG(alpha,i,m)*CalphaGcabs(alpha,j,c)
!                 tmpR4 = tmpR4 + CalphaG(alpha,j,m)*CalphaGcabs(alpha,i,c)
!              ENDDO
              tmpG13 = 0.0E0_realk
              tmpG23 = 0.0E0_realk
              DO beta = 1,NBA
                 tmpG13 = tmpG13 + CalphaG(beta,i,m)*CalphaGcabs(beta,j,c)
                 tmpG23 = tmpG23 + CalphaG(beta,j,m)*CalphaGcabs(beta,i,c)
                 
                 ! tmpG14 = tmpG14 + CalphaGocc(beta,j,m)*CalphaG(beta,i,c)
                 ! tmpG24 = tmpG24 + CalphaGocc(beta,i,m)*CalphaG(beta,j,c)
              ENDDO
              tmp = (Fii(i,i) + Fii(j,j)) 
              EJ3 = EJ3 + tmpG13*tmpG13*tmp
              EK3 = EK3 + tmpG13*tmpG23*tmp
              EJ4 = EJ4 + tmpG23*tmpG23*tmp
              EK4 = EK4 + tmpG23*tmpG13*tmp
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  EJK3 = ED*0.5E0_realk + (7.0_realk/16.0_realk*EJ3+1.0_realk/16.0_realk*EK3)
  EJK4 = ED*0.5E0_realk + (7.0_realk/16.0_realk*EJ4+1.0_realk/16.0_realk*EK4)

 end subroutine ContractTwo4CenterF12IntegralsRI2X

!Exchange Ripjq*Gjpiq
!Ripjq*Gjpiq = (AB|OperR|CD)*(EF|OperR|GH)*Docc(A,G)*Docc(C,E)*Dvirt(B,F)*Dvirt(D,H)
!!$  subroutine ContractTwo4CenterF12IntegralsExchange(OperR,OperG,nocc,nbasis,SETTING,INTSPEC)
!!$    implicit none
!!$    integer,intent(in) :: OperR,OperG,nocc,nbasis
!!$    TYPE(LSSETTING),intent(inout)     :: SETTING
!!$    
!!$    NOFAMILY = ls%setting%SCHEME%NOFAMILY
!!$    ls%setting%SCHEME%NOFAMILY = .TRUE.
!!$    !TODO: 
!!$    nullify(batchsize)
!!$    nullify(batchdim)
!!$    nullify(batchindex)
!!$    nullify(orb2batch)
!!$    nullify(batch2orb)
!!$
!!$    doscreen = ls%setting%SCHEME%CS_SCREEN.OR.ls%setting%SCHEME%PS_SCREEN
!!$
!!$    call build_minimalbatchesofAOS(DECinfo%output,setting,nbasis,&
!!$         & batchsize,batchdim,batchindex,nbatches,orb2Batch,INTSPEC(1))
!!$
!!$    call mem_alloc(batch2orb,nbatchesAB)
!!$    do idx=1,nbatchesAB
!!$       call mem_alloc(batch2orb(idx)%orbindex,batchdim(idx) )
!!$       batch2orb(idx)%orbindex = 0
!!$       batch2orb(idx)%norbindex = 0
!!$    end do
!!$    do iorb=1,nbast
!!$       idx = orb2batch(iorb)
!!$       batch2orb(idx)%norbindex = batch2orb(idx)%norbindex+1
!!$       k = batch2orb(idx)%norbindex
!!$       batch2orb(idx)%orbindex(k) = iorb
!!$    end do
!!$
!!$    call mem_alloc(Dbast,nbasis,nbasis)
!!$    call DGEMM
!!$    call mem_alloc(Docc,nbatches,nbatches)
!!$    call ConvertBASTGabToBatchesGab(nbasis,nbatches,setting,Dbast,Docc,lupri,luerr)
!!$    call mem_dealloc(Dbast)
!!$    MaxDocc = MAXVAL(Docc)
!!$    call mem_alloc(MaxDoccV,nbasis)
!!$    DO J=1,nbasis
!!$       MaxDoccV(J) = MAXVAL(MaxDocc(:,J))
!!$    ENDDO
!!$
!!$    call mem_alloc(Dbast,nbasis,nbasis)
!!$    call DGEMM
!!$    call mem_alloc(Dvirt,nbatches,nbatches)
!!$    call ConvertBASTGabToBatchesGab(nbasis,nbatches,setting,Dbast,Dvirt,lupri,luerr)
!!$    call mem_dealloc(Dbast)
!!$    MaxDvirt = MAXVAL(Dvirt)
!!$    call mem_alloc(MaxDvirtV,nbasis)
!!$    DO J=1,nbasis
!!$       MaxDvirtV(J) = MAXVAL(MaxDvirt(:,J))
!!$    ENDDO
!!$
!!$    call mem_alloc(Rscreen,nbatches,nbatches)
!!$    call II_get_2int_BatchScreenMat(DECinfo%output,DECinfo%output,SETTING,&
!!$         & nbatches,Rscreen,nbasis,OperR)
!!$
!!$    WRITE(lupri,*)'Rscreen'
!!$    call ls_output(Rscreen,1,nbatches,1,nbatches,nbatches,nbatches,1,lupri)
!!$
!!$    MaxRscreen = MAXVAL(Rscreen)
!!$
!!$    call mem_alloc(Gscreen,nbatches,nbatches)
!!$    call II_get_2int_BatchScreenMat(DECinfo%output,DECinfo%output,SETTING,&
!!$         & nbatches,Gscreen,nbasis,OperG)
!!$
!!$    WRITE(lupri,*)'Gscreen'
!!$    call ls_output(Gscreen,1,nbatches,1,nbatches,nbatches,nbatches,1,lupri)
!!$
!!$    MaxGscreen = MAXVAL(Gscreen)
!!$
!!$    Threshold_CS = ls%SETTING%SCHEME%THRESHOLD*ls%SETTING%SCHEME%J_THR
!!$    intThreshold = ls%SETTING%SCHEME%THRESHOLD*ls%SETTING%SCHEME%J_THR
!!$    call II_precalc_DECScreenMat(DecScreenR,lupri,luerr,ls%setting,nbatches,&
!!$         & nbatches,INTSPEC,intThreshold)
!!$    IF(doscreen)then
!!$       call II_getBatchOrbitalScreen(DecScreenR,ls%setting,&
!!$            & nbasis,nbatches,nbatches,batchsize,batchsize,batchindex,batchindex,&
!!$            & batchdim,batchdim,INTSPEC,lupri,luerr)
!!$    endif
!!$
!!$    call II_precalc_DECScreenMat(DecScreenG,lupri,luerr,ls%setting,nbatches,&
!!$         & nbatches,INTSPEC,intThreshold)
!!$    IF(doscreen)then
!!$       call II_getBatchOrbitalScreen(DecScreenG,ls%setting,&
!!$            & nbasis,nbatches,nbatches,batchsize,batchsize,batchindex,batchindex,&
!!$            & batchdim,batchdim,INTSPEC,lupri,luerr)
!!$    endif
!!$    FullRHS = .FALSE.
!!$
!!$    !Ripjq*Gjpiq = (AB|OperR|CD)*(EF|OperR|GH)*Docc(A,G)*Docc(C,E)*Dvirt(B,F)*Dvirt(D,H)
!!$    !starting loop over the most sparse of the 2 operators (create list )
!!$    !MPI PARALLIZE THE D,C Loop
!!$    !E = 0
!!$    !do C
!!$    ! do A
!!$    !   IF(RscreenV(A)*Rscreen(C))THEN
!!$    !    done Nac < nbatchA*nbatchC times 
!!$    !    Calc (A,Bfull|OperR|C,Dfull)
!!$    !    construct Rtensor(A,Bfull,C,Dfull)
!!$    !    construct Rtensor(A,Bfull,C,Hfull)=Dvirt(Dfull,Hfull)*Rtensor(A,Bfull,C,Dfull) !N*N*N*Nac
!!$    !    construct Rtensor(A,Ffull,C,Hfull)=Dvirt(Ffull,Bfull)*Rtensor(A,Bfull,C,Hfull)
!!$    !    do G
!!$    !      construct Rtensor(G,Ffull,C,Hfull)=Docc(G,A)*Rtensor(A,Ffull,C,Hfull)
!!$    !      do E On GPU - While CPU does Integral GPU does DGEMM? 
!!$    !        construct Rtensor(G,Ffull,E,Hfull)=Docc(C,E)*Rtensor(G,Ffull,C,Hfull)
!!$    !        maxR = MAXVAL(Rtensor(G,Ffull,E,Hfull)) 
!!$    !        IF(maxR*GscreenV(E)*GscreenV(G))THEN
!!$    !          Calc (E,Ffull|OperG|G,Hfull) !Modified Screening Threshold with maxR
!!$    !          E = E + (E,Ffull|OperG|G,Hfull)*Rtensor(G,Ffull,E,Hfull)
!!$    !        ENDIF
!!$    !      enddo
!!$    !    enddo
!!$    !   ENDIF
!!$    ! enddo
!!$    !enddo
!!$
!!$
!!$    !TODO VERIFY THAT THE SCREENING IS CORRECT - when you do MaxValRjFiD*MaxGscreen*Gscreen(G,H)*Dvirt(D,H)*MaxCMOV(G)*MaxCMO*maxDvirt
!!$    ! it is not correct as you are summing over the elements - but can do something like what they do in Lapalce 
!!$    !MaxValRjFiD*MaxValGjFiH
!!$    !where MaxValGjFiH = MaxGscreen*SUM_G (MaxCMOV(G)*Gscreen(G,H))
!!$    !Psedo Code 
!!$    !Ripjq*Gjpiq = (AB|OperR|CD)*(EF|OperR|GH)*Docc(A,G)*Docc(C,E)*Dvirt(B,F)*Dvirt(D,H)
!!$    !starting loop over the most sparse of the 2 operators (create list )
!!$    !MPI PARALLIZE THE D,C Loop
!!$    !E = 0
!!$    !do D
!!$    ! do C
!!$    !  IF()
!!$    !   done Ncd < nbatchC*nbatchD times 
!!$    !   (std integral AO to MO speed) expect taking CD screening into account 
!!$    !   construct (AfullBfull|OperR|CD)          
!!$    !   construct Rtensor(Afull,Bfull,C,D)       
!!$    !   construct Rtensor(I,Bfull,C,D)=CMO(I,Afull)*Rtensor(Afull,Bfull,C,D) !DGEMM
!!$    !   Reorder   Rtensor(C,D,I,Bfull)<=Rtensor(I,Bfull,C,D)            
!!$    !   construct Rtensor(J,D,I,Bfull)=Rtensor(J,D,I,Bfull) + CMO(J,C)*Rtensor(C,D,I,Bfull) !DGEMM
!!$    !  ENDIF
!!$    ! enddo
!!$    ! Reorder   Rtensor(J,Bfull,I,D)<=Rtensor(J,D,I,Bfull)
!!$    ! construct Rtensor(J,Ffull,I,D)=Rtensor(J,Bfull,I,D)*Dvirt(Bfull,Ffull)
!!$    ! MaxValRjFiD = MAXVAL(Rtensor(J,Ffull,I,D)) !This is not a small number
!!$    ! do H 
!!$    !  do G
!!$    !   IF(MaxValRjFiD*MaxGscreen*Gscreen(G,H).GT.Threshold_CS)THEN 
!!$    !    done Ngh < nbatchG*nbatchH*nbatchD times 
!!$    !    construct Gtensor(EfullFfull|operG|GH) !USING A MODIFIED SCREENING THRESHOLD *MaxValRjFiD 
!!$    !    MAXVAL = Gtensor(Efull,Ffull,G,H)
!!$    !    IF(MAXVAL*MaxValRjFiD.GT.Threshold_CS)THEN 
!!$    !      construct Gtensor(Efull,Ffull,G,H)
!!$    !      construct Gtensor(J,Ffull,G,H)=CMO(J,Efull)*Gtensor(Efull,Ffull,G,H)                !DGEMM
!!$    !      construct Gtensor(J,Ffull,I,H)=Gtensor(J,Ffull,I,H) + CMO(I,G)*Gtensor(J,Ffull,G,H)
!!$    !    ENDIF
!!$    !   ENDIF
!!$    !  enddo
!!$    !  construct Gtensor(J,Ffull,I,D)=Gtensor(J,Ffull,I,H)*Dvirt(D,H) !DGEMM
!!$    ! enddo
!!$    ! E = E + Rtensor(J,Ffull,I,D)*Gtensor(J,Ffull,I,D)
!!$    !enddo
!!$
!!$    BatchD: do D = 1,nbatches
!!$       dimD = batchdim(D)
!!$       
!!$       BatchC: do C = 1,nbatches
!!$          dimC = batchdim(C)           
!!$          
!!$          !Ripjq*Gjpiq = (AB|OperR|CD)*(EF|OperR|GH)*Docc(A,G)*Docc(C,E)*Dvirt(B,F)*Dvirt(D,H)
!!$          IF(MaxRscreen*Rscreen(C,D)*MaxGscreen*MaxGscreen*MaxDvirt(D)*MaxDocc*MaxDoccV(C)*maxDvirt.GT.Threshold_CS)THEN
!!$
!!$             IF(doscreen)ls%setting%LST_GAB_LHS => DECSCREENR%masterGabLHS
!!$             IF(doscreen)ls%setting%LST_GAB_RHS => DECSCREENR%batchGab(C,D)%p
!!$          
!!$             call II_GET_DECPACKED4CENTER_J_ERI(LUPRI,LUERR,ls%SETTING,&
!!$                  & integralsR,batchindex(C),batchindex(D),batchsize(C),batchsize(D),&
!!$                  & nbast,nbast,dimC,dimD,fullRHS,INTSPEC,intThreshold)
!!$
!!$             !          do batch_iD = 1,dimD
!!$             !           iD = batch2orb(D)%orbindex(batch_iD) !Global index
!!$             !           do batch_iC = 1,dimC
!!$             !            iC = batch2orb(C)%orbindex(batch_iC) !Global index
!!$             !Output integrals(1:nbasis,1:nbasis,batch_iC,batch_iD)
!!$             !Reorder to Radcb <= Rabcd
!!$             MaxValRabcd = MAXVAL(integralsR)
!!$             BatchH: do H = 1,nbatches
!!$                dimH = batchdim(H)
!!$                
!!$                BatchG: do G = 1,nbatches
!!$                   dimG = batchdim(G)           
!!$                   
!!$                   !Ripjq*Gjpiq = (AB|OperR|CD)*(EF|OperR|GH)*Docc(A,G)*Docc(C,E)*Dvirt(B,F)*Dvirt(D,H)
!!$                   IF(MaxValRabcd*MaxGscreen*Gscreen(G,H)*Dvirt(D,H)*MaxDoccV(G)*MaxDoccV(C)*maxDvirt.GT.Threshold_CS)THEN
!!$                      
!!$                      IF(doscreen)ls%setting%LST_GAB_LHS => DECSCREENG%masterGabLHS
!!$                      IF(doscreen)ls%setting%LST_GAB_RHS => DECSCREENG%batchGab(G,H)%p
!!$                      
!!$                      call II_GET_DECPACKED4CENTER_J_ERI(LUPRI,LUERR,ls%SETTING,&
!!$                           & integralsG,batchindex(G),batchindex(H),batchsize(G),batchsize(H),&
!!$                           & nbast,nbast,dimG,dimH,fullRHS,INTSPEC,intThreshold)
!!$                      
!!$                      !          do batch_iH = 1,dimH
!!$                      !           iH = batch2orb(H)%orbindex(batch_iH) !Global index
!!$                      !           do batch_iG = 1,dimG
!!$                      !            iG = batch2orb(G)%orbindex(batch_iG) !Global index
!!$                      !Output integrals(1:nbasis,1:nbasis,batch_iG,batch_iH)
!!$                      !Transform to Gabcd = (EF|OperR|GH)*Docc(A,G)*Docc(C,E)*Dvirt(B,F)*Dvirt(D,H)
!!$
!!$                      !Transform1 Gefgd = (EF|OperR|GH)*Dvirt(H,D)
!!$                      !Transform2 Gcfgd = Docc(C,E)*Gefgd
!!$                      !Reorder Ggdcf <= Gcfgd
!!$                      !Transform1 Ggdcb = Ggdcf*Dvirt(B,F)
!!$                      !Transform2 Gadcb = Docc(A,G)*Ggdcb
!!$                      ! E = Radcb*Gadcb
!!$                      
!!$                   ENDIF
!!$                enddo BatchG
!!$             enddo BatchH
!!$          ENDIF
!!$       enddo BatchC
!!$    enddo BatchD
!!$    
!!$
!!$    !Ripjq*Gjpiq = Rpiqj*Gpjqi
!!$    !starting loop over the most sparse of the 2 operators (create list )
!!$    !do A
!!$    ! do B
!!$    !  do C
!!$    !   do D
!!$    !    construct (AB|OperR|CD)    dim: (nbastOnA,nbastOnB,nbastOnC,nbastOnD)
!!$    !    do E 
!!$    !     do F
!!$    !      do G 
!!$    !       do H
!!$    !        construct (EF|operG|GH) dim: (nbastOnE,nbastOnF,nbastOnG,nbastOnH)  
!!$    !        transform to (AD|operG|CB)  dim: (nbastOnA,nbastOnD,nbastOnC,nbastOnB)
!!$    !        E = (AB|operR|CD)*(AD|operG|CB)   scaling: O(N) 
!!$    !       enddo
!!$    !      enddo
!!$    !     enddo
!!$    !    enddo
!!$    !   enddo
!!$    !  enddo
!!$    ! enddo
!!$    !enddo
!!$    ls%setting%SCHEME%NOFAMILY = NOFAMILY
!!$
!!$  end subroutine ContractTwo4CenterF12IntegralsExchange

#else

  subroutine wangy_dummy_sub12()
      implicit none
  end subroutine wangy_dummy_sub12

#endif

end module fullrimp2f12

