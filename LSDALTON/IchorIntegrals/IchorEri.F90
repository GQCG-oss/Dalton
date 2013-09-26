!> @file
!> Contains the main Ichor integral drivers for calculation electron repulsion integrals
!> Ichor is the "Integral Code Hand Optimized for Rapid evaluation" 
!> Ichor is the ethereal golden fluid that is the blood of the greek gods

!> \brief Main Ichor drivers for the calculation of integrals 
!> based on the McMurchie-Davidson(McM)/Obara Saika(OS)/Head-Gordon-Pople(HGP)/Rys 
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorErimodule
  use IchorprecisionModule
  use IchorCommonModule
!  use IchorEriCoulombintegralMod, only: IchorCoulombIntegral_SSSS
!  use IchorEriCoulombintegralGeneralMod, only: IchorCoulombIntegral_McM_general
  use IchorEriCoulombintegralOBSGeneralMod, only: IchorCoulombIntegral_OBS_general
  use IchorCoulombIntegral_seg_seg_SSSS_mod, only: IchorCoulombIntegral_seg_seg_SSSS
  use IchorMemory
  use IchorGammaTabulationModule

!> Spherical Specification
Integer,parameter :: SphericalParam = 1 !spherical harmonic basis set
!Job Specification
integer,parameter :: IcorJobEri = 1     !4 center 2 electronic repulsion integrals
!Input Spec
integer,parameter :: IcorInputNoInput = 1 !no input in inputstorage (no Density matrix)
!Parallelization specification 
integer,parameter :: IchorParNone = 1     !no parallelization
!> Screening specification 
integer,parameter :: IchorScreen = 1      !default screening including Cauchy-Schwarz screening and QQR
integer,parameter :: IchorScreenNone = 2  !no screening
!> Debug info specification 
integer,parameter :: IchorDebugNone = 1 !no debug printout
!> Integral Algorithm specification 
!> IchorAlgoSpec = IchorAlgoOS = 1 means Obara-Saika (Head-Gordon Pople)
Integer,parameter :: IchorAlgoOS = 1

public:: IchorEri
private
CONTAINS
subroutine IchorEri(nTypesA,MaxNatomsA,MaxnPrimA,MaxnContA,&
     & AngmomOfTypeA,nAtomsOfTypeA,nPrimOfTypeA,nContOfTypeA,&
     & startOrbitalOfTypeA,Acenters,exponentsOfTypeA,ContractCoeffOfTypeA,&
     & nTypesB,MaxNatomsB,MaxnPrimB,MaxnContB,&
     & AngmomOfTypeB,nAtomsOfTypeB,nPrimOfTypeB,nContOfTypeB,&
     & startOrbitalOfTypeB,Bcenters,exponentsOfTypeB,ContractCoeffOfTypeB,&
     & nTypesC,MaxNatomsC,MaxnPrimC,MaxnContC,&
     & AngmomOfTypeC,nAtomsOfTypeC,nPrimOfTypeC,nContOfTypeC,&
     & startOrbitalOfTypeC,Ccenters,exponentsOfTypeC,ContractCoeffOfTypeC,&
     & nTypesD,MaxNatomsD,MaxnPrimD,MaxnContD,&
     & AngmomOfTypeD,nAtomsOfTypeD,nPrimOfTypeD,nContOfTypeD,&
     & startOrbitalOfTypeD,Dcenters,exponentsOfTypeD,ContractCoeffOfTypeD,&
     & SphericalSpec,IchorJobSpec,IchorInputSpec,IchorInputDim1,IchorInputDim2,IchorInputDim3,&
     & InputStorage,IchorParSpec,IchorScreenSpec,IchorDebugSpec,&
     & IchorAlgoSpec,filestorageIdentifier,MaxMem,MaxFileStorage,&
     & MaxMemAllocated,MemAllocated,&
     & OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5,&
     & OutputStorage,lupri)
implicit none
!> nTypesA is the number of different types of shells, each type is defined by 
!> an angular momentum, a number of primitives(nPrim), a number of contracted functions
!> (nCont) a set of exponents and a set of contraction coefficients. [For Center A]
integer,intent(in) :: nTypesA
!> MaxnAtomsA is the maximum number of Atoms that have the same type. [For Center A]
integer,intent(in) :: MaxnAtomsA
!> MaxnPrim is the maximum number of Primitives among all the types. [For Center A]
integer,intent(in) :: MaxnPrimA
!> MaxnPrim is the maximum number of Contracted functions among all the types. [For Center A]
integer,intent(in) :: MaxnContA
!> AngmomOfTypeA is the angular momentum for each type. [For Center A]
Integer,intent(in) :: AngmomOfTypeA(ntypesA)
!> nAtomsOfTypeA is the number of atoms that have the given type. [For Center A]
Integer,intent(in) :: nAtomsOfTypeA(ntypesA)
!> nPrimOfTypeA is the number of primitive for the given type. [For Center A]
Integer,intent(in) :: nPrimOfTypeA(ntypesA)
!> nContOfTypeA is the number of contracted function for the given type. [For Center A]
Integer,intent(in) :: nContOfTypeA(ntypesA)
!> startorbital for this atom of this type
Integer,intent(in) :: startOrbitalOfTypeA(MaxNatomsA,ntypesA)
!> Acenters is the centers of all the atoms that have the given type. [For Center A]
Real(realk),intent(in) :: Acenters(3,MaxNatomsA,ntypesA)
!> exponentsOfTypeA is the nPrim exponents for the given type. [For Center A]
Real(realk),intent(in) :: exponentsOfTypeA(MaxnprimA,ntypesA)
!> ContractCoeffOfTypeA is the contrction coefficient matrix (nPrim,nCont) for the given type. [For Center A]
Real(realk),intent(in) :: ContractCoeffOfTypeA(MaxnprimA,MaxnContA,ntypesA)
!
! Same for Center B
!
integer,intent(in) :: nTypesB,MaxnAtomsB,MaxnPrimB,MaxnContB
Integer,intent(in) :: AngmomOfTypeB(ntypesB),nAtomsOfTypeB(ntypesB)
Integer,intent(in) :: nContOfTypeB(ntypesB),nPrimOfTypeB(ntypesB),startOrbitalOfTypeB(MaxNatomsB,ntypesB)
Real(realk),intent(in) :: Bcenters(3,MaxNatomsB,ntypesB),exponentsOfTypeB(MaxnprimB,ntypesB)
Real(realk),intent(in) :: ContractCoeffOfTypeB(MaxnprimB,MaxnContB,ntypesB)
!
! Same for Center C
!
integer,intent(in) :: nTypesC,MaxnAtomsC,MaxnPrimC,MaxnContC
Integer,intent(in) :: AngmomOfTypeC(ntypesC),nAtomsOfTypeC(ntypesC)
Integer,intent(in) :: nContOfTypeC(ntypesC),nPrimOfTypeC(ntypesC),startOrbitalOfTypeC(MaxNatomsC,ntypesC)
Real(realk),intent(in) :: Ccenters(3,MaxNatomsC,ntypesC),exponentsOfTypeC(MaxnprimC,ntypesC)
Real(realk),intent(in) :: ContractCoeffOfTypeC(MaxnprimC,MaxnContC,ntypesC)
!
! Same for Center D
!
integer,intent(in) :: nTypesD,MaxnAtomsD,MaxnPrimD,MaxnContD
Integer,intent(in) :: AngmomOfTypeD(ntypesD),nAtomsOfTypeD(ntypesD)
Integer,intent(in) :: nContOfTypeD(ntypesD),nPrimOfTypeD(ntypesD),startOrbitalOfTypeD(MaxNatomsD,ntypesD)
Real(realk),intent(in) :: Dcenters(3,MaxNatomsD,ntypesD),exponentsOfTypeD(MaxnprimD,ntypesD)
Real(realk),intent(in) :: ContractCoeffOfTypeD(MaxnprimD,MaxnContD,ntypesD)
!
!> Spherical Specification (SphericalSpec = SphericalParam = 1) means to use Spherical Harmonic basis functions
Integer,intent(in) :: SphericalSpec
!> Job Specification (IcorJob = IcorJobEri = 1) means that the 4 center 2 electron repulsion integrals
!> should be calculated. 
Integer,intent(in) :: IchorJobSpec
!> Input Specification (IchorInputSpec = IcorInputNoInput = 1) means no Input have been provided
Integer,intent(in) :: IchorInputSpec
!> Input dimensions assuming InputStorage(IchorInputDim1,IchorInputDim2,IchorInputDim3)
Integer,intent(in) :: IchorInputDim1,IchorInputDim2,IchorInputDim3
!> InputStorage
real(realk),intent(in) :: InputStorage(IchorInputDim1,IchorInputDim2,IchorInputDim3)
!> Parallelization specification (communicator and other stuff should be set up with other call) 
!> IchorParSpec = IchorParNone = 1 means no parallelization - no OpenMP, no MPI, no GPU
Integer,intent(in) :: IchorParSpec
!> Screening specification 
!> IchorScreenSpec = IchorScreen = 1 means default screening including Cauchy-Schwarz screening and QQR
!> IchorScreenSpec = IchorScreenNone = 2 means no screening
Integer,intent(in) :: IchorScreenSpec
!> Debug info specification 
!> IchorDebugSpec = IchorDebugNone = 1 means no debug info
Integer,intent(in) :: IchorDebugSpec
!> Integral Algorithm specification 
!> IchorAlgoSpec = IchorAlgoOS = 1 means Obara-Saika (Head-Gordon Pople)
Integer,intent(in) :: IchorAlgoSpec
!> Identifier to determine which file should be used to save integrals on disk, This 
!> should be logical unit number, if zero no file is open. 
Integer,intent(in) :: filestorageIdentifier
!> Maximum Memory that the integral program is allowed to use. Zero means no restrictions 
Integer(kind=long),intent(in) :: MaxMem
!> Maximum File size, if zero - no file will be written or read. 
Integer,intent(in) :: MaxFileStorage
!> Maximum Memory used in the program. The Ichor program adds to this value.  
Integer(kind=long),intent(inout) :: MaxMemAllocated
!> Memory allocated in the program. If the input value is not equal to the output value the
!> there is a memory leak inside the program. 
Integer(kind=long),intent(inout) :: MemAllocated
!> Output dimensions assuming 
!> OutputStorage(OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5)
Integer,intent(in) :: OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5
!> OutputStorage
real(realk),intent(inout) :: OutputStorage(OutputDim1,OutputDim2,OutputDim3,OutputDim4,OutputDim5)
!> Logical unit number of output file.
Integer,intent(in) :: lupri
! Local variables
integer :: nPrimP,nContP,nPrimQ,nContQ
integer :: nTABFJW1,nTABFJW2,i1,i2,i3,i4,TotalAngmomABC,TotalAngmomAB,TotalAngmom
integer :: i12,nPasses,offset,K,I,iPass,MaxPasses,iPrimQ,iPrimP,icont
integer :: ndimPass,oldmaxangmomABCD
integer :: ItypeA,ItypeB,itypeC,itypeD,AngmomA,AngmomB,AngmomC,AngmomD
integer :: nPrimA,nPrimB,nContA,nAtomsA,nAtomsB,nAtomsC,nAtomsD,nOrbA
integer :: nDimA,nOrbCompA,iPrimA,iContA,iAtomA,iPrimB,iContB,iAtomB
integer :: iPrimC,iContC,iAtomC,iPrimD,iContD,iAtomD,nContB,nOrbB,nDimB
integer :: nPrimC,nContC,nPrimD,nContD,nOrbC,nOrbD,nOrbCompB,nOrbCompC,nDimC
integer :: nDimD,nOrbCompD,INTPRINT,startA,startB,startC,startD,maxangmomABCD
integer,pointer :: Piprim1(:),Piprim2(:),Qiprim1(:),Qiprim2(:)
logical :: Psegmented,Qsegmented,PQorder,Spherical
real(realk),pointer :: expP(:),Pcent(:),PpreExpFac(:),Qdistance12(:,:),CDAB(:)
real(realk),pointer :: expQ(:),Qcent(:,:),QpreExpFac(:,:),inversexpQ(:),QcentC(:)
REAL(realk),pointer :: TABFJW(:,:),reducedExponents(:),integralPrefactor(:)
Real(realk), parameter :: PIFAC = 34.986836655249725E0_realk !Two*PI**TwoHalf
real(realk) :: Pdistance12(3),e1,e2,X,Y,Z,d2,p,q,p_q
real(realk),pointer :: expA(:),ContractCoeffA(:,:),Acenter(:,:)
real(realk),pointer :: expB(:),ContractCoeffB(:,:),Bcenter(:,:)
real(realk),pointer :: expC(:),ContractCoeffC(:,:),Ccenter(:,:)
real(realk),pointer :: expD(:),ContractCoeffD(:,:),Dcenter(:,:)
INTPRINT=0

! GAMMATABULATION 
!     is this needed for (SSSS) ? is it better to build it several times for 
!     different Angmom combis ? 
Spherical = SphericalSpec.EQ.SphericalParam
oldmaxangmomABCD = -25
DO ItypeA=1,nTypesA
 AngmomA = AngmomOfTypeA(ItypeA)
 nPrimA = nPrimOfTypeA(ItypeA)
 nContA = nContOfTypeA(ItypeA)
 IF (spherical) THEN
    nOrbCompA = 2*(AngmomA+1)-1
 ELSE
    nOrbCompA = (AngmomA+1)*(AngmomA+2)/2
 ENDIF
 nAtomsA = nAtomsOfTypeA(ItypeA)
 nOrbA = nContA*nOrbCompA
 nDimA = nContA*nOrbCompA*nAtomsA
 call mem_ichor_alloc(expA,nPrimA)
 call mem_ichor_alloc(ContractCoeffA,nPrimA,nContA)
 call mem_ichor_alloc(Acenter,3,nAtomsA)
 do iPrimA = 1, nPrimA
    expA(iPrimA) = exponentsOfTypeA(iPrimA,iTypeA)
 enddo
 do iContA = 1, nContA
    do iPrimA = 1, nPrimA
       ContractCoeffA(iPrimA,iContA) = ContractCoeffOfTypeA(iPrimA,iContA,iTypeA)
    enddo
 enddo
 do iAtomA = 1, nAtomsA
    Acenter(1,iAtomA) = Acenters(1,iAtomA,itypeA)
    Acenter(2,iAtomA) = Acenters(2,iAtomA,itypeA)
    Acenter(3,iAtomA) = Acenters(3,iAtomA,itypeA)
 enddo
! expA =>  exponentsOfTypeA(ItypeA)%elms
! ContractCoeffA => ContractCoeffOfTypeA(ItypeA)%elms
! Acenter => CentersOfTypeA(ItypeA)%center
 DO ItypeB=1,nTypesB
  AngmomB = AngmomOfTypeB(ItypeB)
  nPrimB = nPrimOfTypeB(ItypeB)
  nContB = nContOfTypeB(ItypeB)
  IF (spherical) THEN
     nOrbCompB = 2*(AngmomB+1)-1
  ELSE
     nOrbCompB = (AngmomB+1)*(AngmomB+2)/2
  ENDIF
  nAtomsB = nAtomsOfTypeB(ItypeB)
  nOrbB = nContB*nOrbCompB
  nDimB = nOrbB*nAtomsB

  call mem_ichor_alloc(expB,nPrimB)
  call mem_ichor_alloc(ContractCoeffB,nPrimB,nContB)
  call mem_ichor_alloc(Bcenter,3,nAtomsB)
  do iPrimB = 1, nPrimB
     expB(iPrimB) = exponentsOfTypeB(iPrimB,iTypeB)
  enddo
  do iContB = 1, nContB
     do iPrimB = 1, nPrimB
        ContractCoeffB(iPrimB,iContB) = ContractCoeffOfTypeB(iPrimB,iContB,iTypeB)
     enddo
  enddo
  do iAtomB = 1, nAtomsB
     Bcenter(1,iAtomB) = Bcenters(1,iAtomB,itypeB)
     Bcenter(2,iAtomB) = Bcenters(2,iAtomB,itypeB)
     Bcenter(3,iAtomB) = Bcenters(3,iAtomB,itypeB)
  enddo
!  expB =>  exponentsOfTypeB(ItypeB)%elms
!  ContractCoeffB => ContractCoeffOfTypeB(ItypeB)%elms
!  Bcenter => CentersOfTypeB(ItypeB)%center

  TotalAngmomAB = AngmomA + AngmomB 
  nPrimP = nPrimA*nPrimB
  nContP = nContA*nContB
  call mem_ichor_alloc(Piprim1,nPrimP) !not used
  call mem_ichor_alloc(Piprim2,nPrimP) !not used
  IF (nContP.EQ. 1)THEN
     Psegmented = .TRUE.
  ELSE
     Psegmented = .FALSE.
  ENDIF
  call mem_ichor_alloc(expP,nPrimP)
  !it may be easy to include primitive screening on the LHS here
  i12 = 0
  DO i2=1,nPrimB
     e2  = expB(i2)       
     DO i1=1,nPrimA
        e1  = expA(i1)
        i12 = i12 + 1
        expP(i12) = e1 + e2
     ENDDO
  ENDDO
  !here you could allocate (ndimA,ndimB,full,full) 
  DO ItypeC=1,nTypesC
   AngmomC = AngmomOfTypeC(ItypeC)
   nPrimC = nPrimOfTypeC(ItypeC)
   nContC = nContOfTypeC(ItypeC)
   IF (spherical) THEN
      nOrbCompC = 2*(AngmomC+1)-1
   ELSE
      nOrbCompC = (AngmomC+1)*(AngmomC+2)/2
   ENDIF
   nAtomsC = nAtomsOfTypeC(ItypeC)
   nOrbC = nContC*nOrbCompC
   nDimC = nOrbC*nAtomsC

   call mem_ichor_alloc(expC,nPrimC)
   call mem_ichor_alloc(ContractCoeffC,nPrimC,nContC)
   call mem_ichor_alloc(Ccenter,3,nAtomsC)
   do iPrimC = 1, nPrimC
      expC(iPrimC) = exponentsOfTypeC(iPrimC,iTypeC)
   enddo
   do iContC = 1, nContC
      do iPrimC = 1, nPrimC
         ContractCoeffC(iPrimC,iContC) = ContractCoeffOfTypeC(iPrimC,iContC,iTypeC)
      enddo
   enddo
   do iAtomC = 1, nAtomsC
      Ccenter(1,iAtomC) = Ccenters(1,iAtomC,itypeC)
      Ccenter(2,iAtomC) = Ccenters(2,iAtomC,itypeC)
      Ccenter(3,iAtomC) = Ccenters(3,iAtomC,itypeC)
   enddo
!   expC =>  exponentsOfTypeC(ItypeC)%elms
!   ContractCoeffC => ContractCoeffOfTypeC(ItypeC)%elms
!   Ccenter => CentersOfTypeC(ItypeC)%center

   TotalAngmomABC = TotalAngmomAB + AngmomC 

   DO ItypeD=1,nTypesD
    AngmomD = AngmomOfTypeD(ItypeD)

    maxangmomABCD = AngmomA + AngmomB + AngmomC + AngmomD
    IF(maxangmomABCD.NE.oldmaxangmomABCD)THEN
       IF(oldmaxangmomABCD.NE.-25)THEN
          call mem_ichor_dealloc(TABFJW)
       ENDIF
       nTABFJW1 = AngmomA + AngmomB + AngmomC + AngmomD + 3 
       !only need + 3 after Branos change in BUILD_RJ000 
       nTABFJW2 = 1200
       !TABFJW(0:nTABFJW1,0:nTABFJW2)
       call mem_ichor_alloc(TABFJW,nTABFJW1,nTABFJW2,.TRUE.,.TRUE.)
       CALL GAMMATABULATION(lupri,maxangmomABCD,nTABFJW1,nTABFJW2,TABFJW)  
       oldmaxangmomABCD = maxangmomABCD
    ENDIF
    nPrimD = nPrimOfTypeD(ItypeD)
    nContD = nContOfTypeD(ItypeD)
    IF (spherical) THEN
       nOrbCompD = 2*(AngmomD+1)-1
    ELSE
       nOrbCompD = (AngmomD+1)*(AngmomD+2)/2
    ENDIF
    nAtomsD = nAtomsOfTypeD(ItypeD)
    nOrbD = nContD*nOrbCompD
    nDimD = nContD*nOrbCompD*nAtomsD
    call mem_ichor_alloc(expD,nPrimD)
    call mem_ichor_alloc(ContractCoeffD,nPrimD,nContD)
    call mem_ichor_alloc(Dcenter,3,nAtomsD)
    do iPrimD = 1, nPrimD
       expD(iPrimD) = exponentsOfTypeD(iPrimD,iTypeD)
    enddo
    do iContD = 1, nContD
       do iPrimD = 1, nPrimD
          ContractCoeffD(iPrimD,iContD) = ContractCoeffOfTypeD(iPrimD,iContD,iTypeD)
       enddo
    enddo
    do iAtomD = 1, nAtomsD
       Dcenter(1,iAtomD) = Dcenters(1,iAtomD,itypeD)
       Dcenter(2,iAtomD) = Dcenters(2,iAtomD,itypeD)
       Dcenter(3,iAtomD) = Dcenters(3,iAtomD,itypeD)
    enddo
!    expD =>  exponentsOfTypeD(ItypeD)%elms
!    ContractCoeffD => ContractCoeffOfTypeD(ItypeD)%elms
!    Dcenter => CentersOfTypeD(ItypeD)%center

    TotalAngmom = TotalAngmomABC + AngmomD     
    nPrimQ = nPrimC*nPrimD
    nContQ = nContC*nContD
    MaxPasses=nAtomsC*nAtomsD
    call mem_ichor_alloc(Qiprim1,nPrimQ) !not used
    call mem_ichor_alloc(Qiprim2,nPrimQ) !not used
    IF (nContQ.EQ. 1)THEN
       Qsegmented = .TRUE.
    ELSE
       Qsegmented = .FALSE.
    ENDIF
    call mem_ichor_alloc(expQ,nPrimQ)
    call mem_ichor_alloc(inversexpQ,nPrimQ)
    !it may be easy to include primitive screening on the LHS here
    i12 = 0
    DO i2=1,nPrimD
       e2  = expD(i2)       
       DO i1=1,nPrimC
          e1  = expC(i1)
          i12 = i12 + 1
          expQ(i12) = e1 + e2
          inversexpQ(i12) = 1.0E0_realk/(e1 + e2)
       ENDDO
    ENDDO
    call mem_ichor_alloc(reducedExponents,nPrimP*nPrimQ)
    call mem_ichor_alloc(integralPrefactor,nPrimP*nPrimQ)
    IF(TotalAngmom.EQ.0)THEN
       IF (Psegmented.AND.Qsegmented) THEN
          PQorder = .TRUE.
       ELSE
          PQorder = .FALSE.
       ENDIF
    ELSE
       PQorder = .FALSE.
    ENDIF

    IF(PQorder)THEN
       DO iPrimQ = 1, nPrimQ
          q  = expQ(iPrimQ)
          offset = (iPrimQ-1)*nPrimP
          DO iPrimP=1, nPrimP
             p  = expP(iPrimP)
             p_q = p + q
             reducedExponents(iPrimP+offset) = p*q/p_q
             integralPrefactor(iPrimP+offset) = PIFAC/(p*q*SQRT(p_q))
          ENDDO
       ENDDO
    ELSE
       DO iPrimP=1, nPrimP
          p  = expP(iPrimP)
          offset = (iPrimP-1)*nPrimQ
          DO iPrimQ = 1, nPrimQ
             q  = expQ(iPrimQ)
             p_q = p + q
             reducedExponents(iPrimQ+offset) = p*q/p_q
             integralPrefactor(iPrimQ+offset) = PIFAC/(p*q*SQRT(p_q))
          ENDDO
       ENDDO
    ENDIF
    call mem_ichor_alloc(pcent,3*nPrimP)
    call mem_ichor_alloc(PpreExpFac,nPrimP)
    call mem_ichor_alloc(QcentC,3*nPrimQ)
    call mem_ichor_alloc(QpreExpFac,nPrimQ,MaxPasses)
    call mem_ichor_alloc(Qcent,3*nPrimQ,MaxPasses)
    call mem_ichor_alloc(Qdistance12,3,MaxPasses)

    IF (INTPRINT .GE. 10) THEN
       WRITE(lupri,*)'ReducedExponents PQorder',PQorder
       do iPrimP=1,nPrimP*nPrimQ
          WRITE(lupri,'(3X,ES18.9)')reducedExponents(iPrimP)
       enddo
       WRITE(lupri,*)'IntegralPrefactor PQorder',PQorder
       do iPrimP=1,nPrimP*nPrimQ
          WRITE(lupri,'(3X,ES18.9)')integralPrefactor(iPrimP)
       enddo
    END IF

    IF (INTPRINT .GT. 0) THEN
       WRITE(lupri,'(2X,A,I1,A,I1,A,I1,A,I1,A)')'Angmom = (',AngmomA,',',AngmomB,',',AngmomC,',',AngmomD,')'
       WRITE(lupri,'(2X,A,I3,A,I3,A,I3,A,I3,A)')'nPrimitives = (',nPrimA,',',nPrimB,',',nPrimC,',',nPrimD,')'
       WRITE(lupri,'(2X,A,I3,A,I3,A,I3,A,I3,A)')'nContracted = (',nContA,',',nContB,',',nContC,',',nContD,')'
       WRITE(lupri,'(2X,A)')'ExpA and ContractCoeffA'
       do iPrimP=1,nPrimA
          WRITE(lupri,'(3X,5ES18.9/,(3X,5ES18.9))')expA(iPrimP),(ContractCoeffA(iPrimP,iCont),iCont=1,nContA)
       enddo
       WRITE(lupri,'(2X,A)')'ExpB and ContractCoeffB'
       do iPrimP=1,nPrimB
          WRITE(lupri,'(3X,5ES18.9/,(3X,5ES18.9))')expB(iPrimP),(ContractCoeffB(iPrimP,iCont),iCont=1,nContB)
       enddo
       WRITE(lupri,'(2X,A)')'ExpC and ContractCoeffC'
       do iPrimP=1,nPrimC
          WRITE(lupri,'(3X,5ES18.9/,(3X,5ES18.9))')expC(iPrimP),(ContractCoeffC(iPrimP,iCont),iCont=1,nContC)
       enddo
       WRITE(lupri,'(2X,A)')'ExpD and ContractCoeffD'
       do iPrimP=1,nPrimD
          WRITE(lupri,'(3X,5ES18.9/,(3X,5ES18.9))')expD(iPrimP),(ContractCoeffD(iPrimP,iCont),iCont=1,nContD)
       enddo
    ENDIF
    IF(TotalAngmom.EQ.0)THEN
       !============================================================================================
       !                                        SSSS integrals
       !At this point we know all info about the integral except for the coordinates. 
       !Can be used for screening? Different screening for different types (SS|SS) 
       !APE type screening? 
       !For SSSS it should be exact to do
       !exp(-mu_AB R_AB^2)exp(-mu_CD R_CD^2) R000 and use
       !R000(alpha*R_pq^2) = F0(alpha*R_pq^2) = sqrt(pi/(4*x))erf(x) ;x=alpha*R_pq^2
       !R000(alpha*R_pq^2) > sqrt(pi/(4*x)) ;x=alpha*R_pq^2
       !look at possible things for SSSP type things. 
       !determine memory requirements
       !determine size of passes : nPasses
       !determine screening noScreen, LINK stuff - but at atomic level not batch level(althoug same)
       !Form active Dmat 4 (using screening) matrices (D_BD,D_AD,D_BC,D_AC)
       !qexp
       !ACC,BCC,CCC,DCC
       !Qsegmented
       !Psegmented
       !ThermiteWorkMem
       !============================================================================================
!     call IchorTimer('START',TSTART,TEND,LUPRI)
     !INFO 
     DO IatomA = 1,nAtomsA
      startA = startOrbitalOfTypeA(iAtomA,ItypeA)
      DO IatomB = 1,nAtomsOfTypeB(ItypeB)
       startB = startOrbitalOfTypeB(iAtomB,ItypeB)
!      IF(noScreenAB(IatomB,IatomA))THEN
         !sort combined atom list with respect to distance ? 
         !only difference is the centers involved. 
       i12 = 0
       DO i2=1,nPrimB
        e2  = expB(i2)       
        DO i1=1,nPrimA
         e1  = expA(i1)
         offset = i12*3
         i12 = i12 + 1
         pcent(1+offset) = (e1*Acenter(1,IatomA) + e2*Bcenter(1,IatomB))/expP(i12)
         Pcent(2+offset) = (e1*Acenter(2,IatomA) + e2*Bcenter(2,IatomB))/expP(i12)
         Pcent(3+offset) = (e1*Acenter(3,IatomA) + e2*Bcenter(3,IatomB))/expP(i12)
         X = Acenter(1,IatomA) - Bcenter(1,IatomB)
         Y = Acenter(2,IatomA) - Bcenter(2,IatomB)
         Z = Acenter(3,IatomA) - Bcenter(3,IatomB)
         Pdistance12(1) = X
         Pdistance12(2) = Y
         Pdistance12(3) = Z
         d2 = X*X
         d2 = d2 + Y*Y
         d2 = d2 + Z*Z
         PpreExpFac(i12) = exp(-e1*e2/(e1+e2)*d2)
         IF (Psegmented) THEN
            PpreExpFac(i12) = PpreExpFac(i12)*ContractCoeffA(i1,1)*ContractCoeffB(i2,1)
         ENDIF
        ENDDO
       ENDDO

       IF (INTPRINT .GE. 10) THEN
          WRITE(lupri,'(A,I7,A)')'PpreExpFac(',nPrimP,')'
          do iPrimP=1,nPrimP
             WRITE(lupri,'(3X,ES18.9)') PpreExpFac(iPrimP)
          enddo
       END IF

       !MAKE LIST 
       nPasses = nAtomsC*nAtomsD
       !
       CALL MEM_ICHOR_ALLOC(CDAB,nContA*nContB*ndimC*ndimD)      
       !make qcenter 
       !make QpreExpFac
        DO IatomC = 1,nAtomsC
         i12 = 0
         DO i2=1,nPrimD
            e2  = expD(i2)       
            DO i1=1,nPrimC
               e1  = expC(i1)
               offset = i12*3
               i12 = i12 + 1
               qcentC(1+offset) = (e1*Ccenter(1,IatomC))*inversexpQ(i12)
               qcentC(2+offset) = (e1*Ccenter(2,IatomC))*inversexpQ(i12)
               qcentC(3+offset) = (e1*Ccenter(3,IatomC))*inversexpQ(i12)
            ENDDO
         ENDDO
         DO IatomD = 1,nAtomsD
            iPass = IatomC + (IatomD-1)*nAtomsC !for now
!          IF(noScreenCD(IatomD,IatomC))THEN
            i12 = 0
            DO i2=1,nPrimD
               e2  = expD(i2)       
               DO i1=1,nPrimC
                  e1  = expC(i1)
                  offset = i12*3
                  i12 = i12 + 1
                  qcent(1+offset,iPass) = qcentC(1+offset) + (e2*Dcenter(1,IatomD))*inversexpQ(i12)
                  qcent(2+offset,iPass) = qcentC(2+offset) + (e2*Dcenter(2,IatomD))*inversexpQ(i12)
                  qcent(3+offset,iPass) = qcentC(3+offset) + (e2*Dcenter(3,IatomD))*inversexpQ(i12)
                  X = Ccenter(1,IatomC) - Dcenter(1,IatomD)
                  Y = Ccenter(2,IatomC) - Dcenter(2,IatomD)
                  Z = Ccenter(3,IatomC) - Dcenter(3,IatomD)
                  Qdistance12(1,iPass) = X
                  Qdistance12(2,iPass) = Y
                  Qdistance12(3,iPass) = Z
                  d2 = X*X
                  d2 = d2 + Y*Y
                  d2 = d2 + Z*Z
                  QpreExpFac(i12,iPass) = exp(-e1*e2/(e1+e2)*d2)
                  IF (Qsegmented) THEN
                   QpreExpFac(i12,iPass) = QpreExpFac(i12,iPass)*ContractCoeffC(i1,1)*ContractCoeffD(i2,1)
                  ENDIF
                  !integralPrefactor(nPrimPQ) could be modified instead of QpreExpFac to avoid nPass cost
               ENDDO
            ENDDO
!          ENDIF
         ENDDO
      ENDDO

      IF (INTPRINT .GE. 10) THEN
         WRITE(lupri,'(A,I7,A1,I4,A1)')'QpreExpFac(',nPrimQ,',',nPasses,')'
         DO IPass = 1,nPasses
            WRITE(lupri,*)'IPass',iPass
            do iPrimQ=1,nPrimQ
               WRITE(lupri,'(3X,ES18.9)') QpreExpFac(iPrimQ,iPass)
            enddo
          ENDDO
       END IF

       IF (Psegmented.AND.Qsegmented) THEN
          call IchorCoulombIntegral_seg_seg_SSSS(nPrimP,nPrimQ,nPasses,&
               & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
               & reducedExponents,integralPrefactor,PQorder,CDAB)
       ELSE


!!$          WRITE(lupri,*)'IchorCoulombIntegral_McM_general'
                call IchorQuit('IchorCoulombIntegral_McM_general',-1)
!!$             call IchorCoulombIntegral_McM_general(nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,&
!!$                  & nPrimQ,nPrimP*nPrimQ,nPasses,MaxPasses,intprint,lupri,&
!!$                  & nContA,nContB,nContC,nContD,nContP,nContQ,expP,expQ,&
!!$                  & ContractCoeffA,ContractCoeffB,ContractCoeffC,ContractCoeffD,&
!!$                  & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
!!$                  & Qiprim1,Qiprim2,Piprim1,Piprim2,expA,expB,expC,expD,&
!!$                  & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
!!$                  & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB)
!          call IchorCoulombIntegral_SSSS(nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,&
!               & nPrimQ,nPrimP*nPrimQ,nPasses,MaxPasses,intprint,lupri,&
!               & nContA,nContB,nContC,nContD,nContP,nContQ,expP,expQ,&
!               & ContractCoeffA,ContractCoeffB,ContractCoeffC,ContractCoeffD,&
!               & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
!               & Qiprim1,Qiprim2,Piprim1,Piprim2,expA,expB,expC,expD,&
!               & Qsegmented,Psegmented,reducedExponents,integralPrefactor,CDAB)
       ENDIF
      !Distribute to Active Kmat 4 matrices (K_AC,K_BC,K_AD,K_BD)
             
        
      !must be written out so to vectorize
      !CDAB(IC,ID) = CDAB(IC,ID,StartA+IA,startB+IB)
      DO IatomC = 1,nAtomsC
       DO IatomD = 1,nAtomsD
        OutputStorage(IatomA,IatomB,IatomC,IatomD,1)=CDAB(IatomC+(IatomD-1)*nAtomsC)
       ENDDO
      ENDDO
      CALL Mem_ichor_dealloc(CDAB)
      ENDDO
     ENDDO
!     call IchorTimer('TYPE',TSTART,TEND,LUPRI)
     !distribute to full Kmat
     !free ThermiteWorkMem
     !

    ELSE
       !============================================================================================
       !                                        General integrals
!       call IchorTimer('START',TSTART,TEND,LUPRI)
       DO IatomA = 1,nAtomsA
          startA = startOrbitalOfTypeA(iAtomA,ItypeA)
          DO IatomB = 1,nAtomsOfTypeB(ItypeB)
             startB = startOrbitalOfTypeB(iAtomB,ItypeB)
             !sort combined atom list with respect to distance ? which distance AB,CD,PQ 
             !only difference is the centers involved. 
             i12 = 0
             DO i2=1,nPrimB
                e2  = expB(i2)       
                DO i1=1,nPrimA
                   e1  = expA(i1)
                   offset = i12*3
                   i12 = i12 + 1
                   pcent(1+offset) = (e1*Acenter(1,IatomA) + e2*Bcenter(1,IatomB))/expP(i12)
                   Pcent(2+offset) = (e1*Acenter(2,IatomA) + e2*Bcenter(2,IatomB))/expP(i12)
                   Pcent(3+offset) = (e1*Acenter(3,IatomA) + e2*Bcenter(3,IatomB))/expP(i12)
                   X = Acenter(1,IatomA) - Bcenter(1,IatomB)
                   Y = Acenter(2,IatomA) - Bcenter(2,IatomB)
                   Z = Acenter(3,IatomA) - Bcenter(3,IatomB)
                   Pdistance12(1) = X
                   Pdistance12(2) = Y
                   Pdistance12(3) = Z
                   d2 = X*X
                   d2 = d2 + Y*Y
                   d2 = d2 + Z*Z
                   PpreExpFac(i12) = exp(-e1*e2/(e1+e2)*d2)
                   IF (Psegmented) THEN
                      PpreExpFac(i12) = PpreExpFac(i12)*ContractCoeffA(i1,1)*ContractCoeffB(i2,1)
                   ENDIF
                ENDDO
             ENDDO

             IF (INTPRINT .GE. 10) THEN
                WRITE(lupri,*)'PpreExpFac'
                do iPrimP=1,nPrimP
                   WRITE(lupri,'(3X,ES18.9)') PpreExpFac(iPrimP)
                enddo
             END IF
             
             !MAKE LIST 
             nPasses = nAtomsC*nAtomsD
             CALL MEM_ICHOR_ALLOC(CDAB,nOrbA*nOrbB*ndimC*ndimD)      

             !make qcenter 
             !make QpreExpFac
             !       IF(noScreenAB(IatomB,IatomA))THEN
             DO IatomC = 1,nAtomsC
                i12 = 0
                DO i2=1,nPrimD
                   e2  = expD(i2)       
                   DO i1=1,nPrimC
                      e1  = expC(i1)
                      offset = i12*3
                      i12 = i12 + 1
                      qcentC(1+offset) = (e1*Ccenter(1,IatomC))*inversexpQ(i12)
                      qcentC(2+offset) = (e1*Ccenter(2,IatomC))*inversexpQ(i12)
                      qcentC(3+offset) = (e1*Ccenter(3,IatomC))*inversexpQ(i12)
                   ENDDO
                ENDDO
                DO IatomD = 1,nAtomsD
                   iPass = IatomC + (IatomD-1)*nAtomsC !for now
                   !          IF(noScreenCD(IatomD,IatomC))THEN
                   i12 = 0
                   DO i2=1,nPrimD
                      e2  = expD(i2)       
                      DO i1=1,nPrimC
                         e1  = expC(i1)
                         offset = i12*3
                         i12 = i12 + 1
                         qcent(1+offset,iPass) = qcentC(1+offset) + (e2*Dcenter(1,IatomD))*inversexpQ(i12)
                         qcent(2+offset,iPass) = qcentC(2+offset) + (e2*Dcenter(2,IatomD))*inversexpQ(i12)
                         qcent(3+offset,iPass) = qcentC(3+offset) + (e2*Dcenter(3,IatomD))*inversexpQ(i12)
                         X = Ccenter(1,IatomC) - Dcenter(1,IatomD)
                         Y = Ccenter(2,IatomC) - Dcenter(2,IatomD)
                         Z = Ccenter(3,IatomC) - Dcenter(3,IatomD)
                         Qdistance12(1,iPass) = X
                         Qdistance12(2,iPass) = Y
                         Qdistance12(3,iPass) = Z
                         d2 = X*X
                         d2 = d2 + Y*Y
                         d2 = d2 + Z*Z
                         QpreExpFac(i12,iPass) = exp(-e1*e2/(e1+e2)*d2)
                         IF (Qsegmented) THEN
                            QpreExpFac(i12,iPass) = QpreExpFac(i12,iPass)*ContractCoeffC(i1,1)*ContractCoeffD(i2,1)
                         ENDIF
                         !integralPrefactor(nPrimPQ) could be modified instead of QpreExpFac to avoid nPass cost
                      ENDDO
                   ENDDO
                   !          ENDIF
                ENDDO
             ENDDO
             
             IF (INTPRINT .GE. 10) THEN
                WRITE(lupri,'(A,I7,A1,I4,A1)')'QpreExpFac(',nPrimQ,',',nPasses,')'
                DO IPass = 1,nPasses
                   WRITE(lupri,*)'IPass',iPass
                   do iPrimQ=1,nPrimQ
                      WRITE(lupri,'(3X,ES18.9)') QpreExpFac(iPrimQ,iPass)
                   enddo
                ENDDO
             END IF
             !which method depend on angular momentum and stuff!
             IF(.TRUE.)THEN !OBS)THEN
                call IchorCoulombIntegral_OBS_general(nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,&
                     & nPrimQ,nPrimP*nPrimQ,nPasses,MaxPasses,intprint,lupri,&
                     & nContA,nContB,nContC,nContD,nContP,nContQ,expP,expQ,&
                     & ContractCoeffA,ContractCoeffB,ContractCoeffC,ContractCoeffD,&
                     & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
                     & Qiprim1,Qiprim2,Piprim1,Piprim2,expA,expB,expC,expD,&
                     & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
                     & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB,&
                     & Acenter(:,1),Bcenter(:,1),Ccenter,Dcenter,nAtomsC,nAtomsD,Spherical)
                !output CDAB(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP,nAtomsC,nAtomsD)
                ndimPass = nOrbCompA*nContA*nOrbCompB*nContB*nOrbCompC*nContC*nOrbCompD*nContD
                !             write(lupri,*)'CDAB'
                !             call output(CDAB,1,ndimPass,1,nPasses,ndimPass,nPasses,1,lupri)
                !This must be optimized as well 
!                1. for segmented all nCont=1
!                2. nOrbComp is known early on (1,3,5,7,9,...) and is the same for all passes
                

                call IchorDistribute(nAtomsC,nAtomsD,startOrbitalOfTypeD(1:nAtomsD,ItypeD),&
                     & startOrbitalOfTypeC(1:nAtomsC,ItypeC),startA,startB,&
                     & AngmomA,AngmomB,AngmomC,AngmomD,nContA,nContB,nContC,nContD,&
                     & OutputStorage,Outputdim1,Outputdim2,Outputdim3,Outputdim4,&
                     & CDAB,nOrbA,nOrbB,nOrbC,nOrbD,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD)

             ELSE
                call IchorQuit('IchorCoulombIntegral_McM_general',-1)
                !             ELSEIF(McM)THEN
!!$                call IchorCoulombIntegral_McM_general(nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,&
!!$                     & nPrimQ,nPrimP*nPrimQ,nPasses,MaxPasses,intprint,lupri,&
!!$                     & nContA,nContB,nContC,nContD,nContP,nContQ,expP,expQ,&
!!$                     & ContractCoeffA,ContractCoeffB,ContractCoeffC,ContractCoeffD,&
!!$                     & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
!!$                     & Qiprim1,Qiprim2,Piprim1,Piprim2,expA,expB,expC,expD,&
!!$                     & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
!!$                     & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB)

                ndimPass = nOrbCompA*nContA*nOrbCompB*nContB*nOrbCompC*nContC*nOrbCompD*nContD
                !             write(lupri,*)'CDAB'
                !             call output(CDAB,1,ndimPass,1,nPasses,ndimPass,nPasses,1,lupri)
                DO IatomD = 1,nAtomsD
                 DO IatomC = 1,nAtomsC
                  DO I1 = 1,nOrbCompA*nContA
                   DO I2 = 1,nOrbCompB*nContB
                    DO I3 = 1,nOrbCompC*nContC
                     DO I4 = 1,nOrbCompD*nContD
                      OutputStorage(startA+I1,startB+I2,I3+(IatomC-1)*nOrbC,I4+(IatomD-1)*nOrbD,1)=&
                        & CDAB(I3+(I4-1)*nOrbC+(I1-1)*nOrbC*nOrbD+(I2-1)*nOrbC*nOrbD*nOrbA+(IatomC-1)*ndimPass+(IatomD-1)*ndimPass*nAtomsC)
                     ENDDO
                    ENDDO
                   ENDDO
                  ENDDO
                 ENDDO
                ENDDO

                !             ELSE !RY
                
                !             ENDIF
             ENDIF
             
             !Distribute to Active Kmat 4 matrices (K_AC,K_BC,K_AD,K_BD)
             
             
             !must be written out so to vectorize
             !CDAB(IC,ID) = CDAB(IC,ID,StartA+IA,startB+IB)
!             print*,'CDAB',CDAB
!             print*,'dim1',dim1
!             print*,'dim2',dim2
!             print*,'dim3',dim3
!             print*,'dim4',dim4
!             print*,'ndimA',ndimA
!             print*,'ndimB',ndimB
!             print*,'ndimC',ndimC
!             print*,'ndimD',ndimD
!             write(lupri,*)'integrals'
!             call output(integrals,1,dim1*dim2,1,dim3*dim4,dim1*dim2,dim3*dim4,1,lupri)

             CALL Mem_ichor_dealloc(CDAB)
          ENDDO
       ENDDO
!       call IchorTimer('TYPE',TSTART,TEND,LUPRI)
       !distribute to full Kmat
       !free ThermiteWorkMem
       !
!       call lsquit('not implemented',-1)
    ENDIF
    call mem_ichor_dealloc(pcent)
    call mem_ichor_dealloc(PpreExpFac)
    call mem_ichor_dealloc(QcentC)
    call mem_ichor_dealloc(QpreExpFac)
    call mem_ichor_dealloc(Qcent)
    call mem_ichor_dealloc(Qdistance12)
    call mem_ichor_dealloc(Qiprim1)
    call mem_ichor_dealloc(Qiprim2)
    call mem_ichor_dealloc(expQ)
    call mem_ichor_dealloc(inversexpQ)
    call mem_ichor_dealloc(expD)
    call mem_ichor_dealloc(ContractCoeffD)
    call mem_ichor_dealloc(Dcenter)
   ENDDO !typeD
   call mem_ichor_dealloc(expC)
   call mem_ichor_dealloc(ContractCoeffC)
   call mem_ichor_dealloc(Ccenter)
  ENDDO !typeC
  !here you could distribute the (ndimA,ndimB,full,full) to (full,full,full,full) ?

  call mem_ichor_dealloc(expP)
  call mem_ichor_dealloc(Piprim1) 
  call mem_ichor_dealloc(Piprim2) 

  call mem_ichor_dealloc(expB)
  call mem_ichor_dealloc(ContractCoeffB)
  call mem_ichor_dealloc(Bcenter)
 ENDDO !typeB
 call mem_ichor_dealloc(expA)
 call mem_ichor_dealloc(ContractCoeffA)
 call mem_ichor_dealloc(Acenter)
ENDDO !typeA
call mem_ichor_dealloc(TABFJW)
end subroutine IchorEri

subroutine IchorDistribute(nAtomsC,nAtomsD,startorbitalD,&
     & startorbitalC,startA,startB,AngmomA,AngmomB,AngmomC,AngmomD,&
     & nContA,nContB,nContC,nContD,integrals,dim1,dim2,dim3,dim4,&
     & CDAB,nOrbA,nOrbB,nOrbC,nOrbD,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD)
  implicit none
  integer,intent(in) :: nAtomsC,nAtomsD,startA,startB,AngmomA,AngmomB,AngmomC,AngmomD
  integer,intent(in) :: nContA,nContB,nContC,nContD,dim1,dim2,dim3,dim4
  integer,intent(in) :: nOrbA,nOrbB,nOrbC,nOrbD,nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD
  integer,intent(in) :: startorbitalD(nAtomsD),startorbitalC(nAtomsC)
  !output CDAB(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContQ,nContP,nAtomsC,nAtomsD)
!  real(realk),intent(in) :: CDAB(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContC*nContD,nContA*nContB,nAtomsC,nAtomsD)
  real(realk),intent(in) :: CDAB(nOrbCompA,nOrbCompB,nOrbCompC,nOrbCompD,nContC*nContD,nContA*nContB,nAtomsC*nAtomsD)
  real(realk),intent(inout) :: integrals(dim1,dim2,dim3,dim4)
  !
  integer :: IatomD,IatomC,I1,I2,I3,I4,startD,startC,iAngA,iAngB,iAngC,iAngD
  integer :: iContA,iContB,iContC,iContD,iPassQ,iContP,iContQ
!  IF(segmented)THEN
 !  IF(AngmomA.EQ.  0)THEN
 !   IF(AngmomB.EQ.  0)THEN
 !    IF(AngmomC.EQ.  0)THEN
 !     IF(AngmomD.EQ.  0)THEN
!     CALL IchorDistributeSeg0000(nAtomsC,nAtomsD,startorbitalD,&
!          & startorbitalC,startA,startB,integrals,dim1,dim2,dim3,dim4,&
!          & CDAB)
 !      CALL IchorDistributeSeg0000(nAtomsC,nAtomsD,startorbitalD,&
 !           & startorbitalC,startA,startB,&
 !           & nContA,nContB,nContC,nContD,integrals,dim1,dim2,dim3,dim4,&
 !           & CDAB)
 !     ENDIF
 !    ENDIF
 !   ENDIF
 !  ENDIF
!  ELSE
   DO IatomD = 1,nAtomsD
    startD = startorbitalD(iAtomD)
    DO IatomC = 1,nAtomsC
     startC = startorbitalC(iAtomC)
     IpassQ = IatomC + (IatomD-1)*nAtomsC
     iContQ = 0
     DO iContD = 1,nContD
      DO iContC = 1,nContC
       iContQ = iContQ + 1
       iContP = 0
       DO iContB = 1,nContB
        DO iContA = 1,nContA
         iContP = iContP + 1

         DO iAngD = 1,nOrbCompD
          I4 = startD + iAngD + (iContD-1)*nOrbCompD
          DO iAngC = 1,nOrbCompC
           I3 = startC + iAngC + (iContC-1)*nOrbCompC
           DO iAngB = 1,nOrbCompB
            I2 = startB + iAngB + (iContB-1)*nOrbCompB
            DO iAngA = 1,nOrbCompA
             I1 = startA + iAngA + (iContA-1)*nOrbCompA
             integrals(I1,I2,I3,I4) = CDAB(iAngA,iAngB,iAngC,iAngD,iContQ,iContP,IpassQ)
!             write(6,'(A,I2,A,I2,A,I2,A,I2,A,F22.10)')'integrals(',I1,',',I2,',',I3,',',I4,')',integrals(I1,I2,I3,I4)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
!  ENDIF
end subroutine IchorDistribute

 !alot of jumping around in mem - no matter what since integral is save (abcd)
subroutine IchorDistributeSeg0000(nAtomsC,nAtomsD,startorbitalD,&
     & startorbitalC,startA,startB,integrals,dim1,dim2,dim3,dim4,&
     & CDAB)
  implicit none
  integer,intent(in) :: nAtomsC,nAtomsD,startA,startB
  integer,intent(in) :: dim1,dim2,dim3,dim4
  integer,intent(in) :: startorbitalD(nAtomsD),startorbitalC(nAtomsC)
  real(realk),intent(in) :: CDAB(nAtomsC,nAtomsD)
  real(realk),intent(inout) :: integrals(dim1,dim2,dim3,dim4)
  !
  integer :: IatomD,IatomC,startD,startC
  DO IatomD = 1,nAtomsD
     startD = startorbitalD(iAtomD)
     DO IatomC = 1,nAtomsC
        startC = startorbitalC(iAtomC)
        integrals(startA+1,startB+1,startC+1,startD+1)=CDAB(IatomC,IAtomD)
     ENDDO
  ENDDO
end subroutine IchorDistributeSeg0000
 
END MODULE IchorErimodule
