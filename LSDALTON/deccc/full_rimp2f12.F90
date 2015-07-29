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
   real(realk) :: E_21,E_22,E_23, E_F12, E_21C
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

    E_21 = 0.0E0_realk
    E_21 = EV1 + EV2 + EV3 + EV4 +EV5 

    if(DECinfo%F12debug) then
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
       write(DECinfo%output,'(1X,a,f15.16)') " E21_Vsum:     ", E_21
    end if

    E_22 = 0.0E0_realk
    E_22 = EX1 + EX2 + EX3 + EX4 

    if(DECinfo%F12debug) then
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

    if(DECinfo%F12debug) then
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

    if(DECinfo%F12debug) then
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




end subroutine full_canonical_rimp2_f12

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

