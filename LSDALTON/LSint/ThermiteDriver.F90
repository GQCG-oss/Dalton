!> @file
!> Contains the main Thermite integral drivers

!> \brief Main Thermite drivers for the calculation of integrals 
!> based on the McMurchie-Davidson scheme, using (Turbo-)hermite functions. (PCCP 2007 9, 4771) 
!> \author T. Kjaergaard and S. Reine
!> \date 2008 
MODULE integraldriver
  use TYPEDEF
  use READMOLEFILE
  use BuildBasisSet
  use ODbatches
  use OD_type
  use sphcart_matrices
  use precision
  use lstiming
  use thermite_integrals
  use thermite_OD
  use ls_util
  use Thermite_prop
  use ThermiteMem_module
  use LSparameters
  use Fundamental
  use OverlapType
  use thermite_distribute
  use thermite_distribute2
  use thermite_distributeGen
  use thermite_distributeDEC
  use thermite_distributeK
  use thermite_distributeK2
  use math_fun
SAVE
TYPE LINKshell
integer,pointer :: Belms(:)
integer,pointer :: IODelms(:)
integer(kind=short),pointer :: RED_GAB(:)
INTEGER         :: DIM
END TYPE LINKshell
public::MAIN_INTEGRAL_DRIVER, jmatclassical, gradclassical,&
     & jmatclassicalmat, gradclassicalgrad, electronnuclearclassic
private
CONTAINS
!> \brief Main integral driver
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!>  This is the Main integral driver, which sets op the overlap distributions 
!>  (call ODbatches) from the Atomic orbitals given in input INPUT%AO.
!>  The routine then calls either the main McMurchie-Davidson driver or 
!>  specialised routines. For the calculation of exchange integrals 
!>  a specialised LinK routine is used and for the calculation of coulomb 
!>  integrals, using the jengine algorithm, another specialised routine 
!>  have been implemented. 
!> 
!>  \param LUPRI the logical unit number for the output file
!>  \param IPRINT the printlevel, determining how much output should be generated
!>  \param INPUT the integral input specifications, contains all info about what to do
!>  \param INTOUT the integral output specifications, determines how the output should be given
SUBROUTINE MAIN_INTEGRAL_DRIVER(LUPRI,IPRINT,INPUT,INTOUT)
implicit none
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(INTEGRALOUTPUT) :: INTOUT
integer              :: LUPRI,IPRINT
TYPE(ODITEM),target  :: OD_LHSt,OD_RHSt
TYPE(ODITEM),pointer :: OD_LHS,OD_RHS
!Character(len=80)    :: IDENTIFIER
!
Real(realk)         :: TS,TE
INTEGER             :: idmat, natoms
!
IF(IPRINT .GT. 4)call PRINT_INTEGRALINPUT(INPUT,LUPRI)

!Quick Exit if possible
IF (INPUT%sameODs .AND. (INPUT%sameLHSaos.NEQV.INPUT%sameRHSaos))THEN
   CALL LSQUIT('Progamming error, sameODs but sameLHSaos.NE.sameRHSos!',lupri)
ENDIF

IF(IPRINT .GT. 5) WRITE(LUPRI,*)'CALLING CREATE OD'
CALL Create_ODbatches(OD_LHSt,INPUT,'LHS',LUPRI)
nullify(OD_LHS)
OD_LHS => OD_LHSt
IF(OD_LHS%nbatches .NE. 0)THEN
   nullify(OD_RHS) 
   IF(INPUT%sameODs)THEN
      OD_RHS => OD_LHSt
   ELSE
      CALL Create_ODbatches(OD_RHSt,INPUT,'RHS',LUPRI)
      OD_RHS => OD_RHSt
   ENDIF
   IF(OD_RHS%nbatches .NE. 0)THEN      
      CALL PRINT_OD(OD_LHS,LUPRI,IPRINT)
      CALL PRINT_OD(OD_RHS,LUPRI,IPRINT)
      IF (INPUT%DO_JENGINE) THEN
!         CALL LSTIMER('START ',TS,TE,LUPRI)
         CALL Jengine(OD_LHS,OD_RHS,INPUT,INTOUT,LUPRI,IPRINT)
!         CALL LSTIMER('Jengine',TS,TE,LUPRI)
      ELSEIF(INPUT%DO_LINK) THEN
        CALL LINK_DRIVER(OD_LHS,OD_RHS,INPUT,INTOUT,LUPRI,IPRINT)
      ELSE
         CALL McMurchieDavidson(OD_LHS,OD_RHS,INPUT,INTOUT,LUPRI,IPRINT)
      ENDIF
      IF(.NOT.Input%sameODs)THEN
         CALL FREE_ODitem(OD_RHS)
      ENDIF
   ENDIF
   CALL FREE_ODitem(OD_LHS)
ENDIF

END SUBROUTINE MAIN_INTEGRAL_DRIVER

!> \brief McMurchie-Davidson based integral driver
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!>  This is the Main McMurchie-Davidson based integral driver, which 
!>  sets up overlaps P and Q and loops over the left hand side (P) overlaps
!>  and 
!>  The routine then calls either the main McMurchie-Davidson driver or 
!>  specialised routines. For the calculation of exchange integrals 
!>  a specialised LinK routine is used and for the calculation of coulomb 
!>  integrals, using the jengine algorithm, another specialised routine 
!>  have been implemented. 
!> 
!>  \param OD_LHS the ODbatches belonging to the Left hand side
!>  \param OD_RHS the ODbatches belonging to the Reft hand side
!>  \param INPUT the integral input specifications, contains all info about what to do
!>  \param OUTPUT the integral output specifications, determines how the output should be given
!>  \param LUPRI the logical unit number for the output file
!>  \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE McMurchieDavidson(OD_LHS,OD_RHS,INPUT,OUTPUT,LUPRI,IPRINT)
implicit none
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(INTEGRALOUTPUT) :: OUTPUT
integer              :: LUPRI,IPRINT
TYPE(ODITEM)         :: OD_LHS,OD_RHS
!
Integer               :: ILHS,IRHS,Start_RHS,End_RHS
TYPE(TUVitem),target  :: SharedTUV
TYPE(Integralitem)    :: Integral
TYPE(Allocitem)       :: Allocations
TYPE(Overlap)         :: P,Q2
TYPE(Overlap),pointer :: Q(:)
TYPE(Overlap),pointer :: PassQ(:)
TYPE(Integrand)       :: PQ
Integer,pointer       :: maxPassesFortypes(:)
Integer,pointer       :: ODpassesIndex(:),nOverlapOfPassType(:) 
Integer,pointer       :: ODTypeIndex(:) 
logical               :: dopasses,screen 
integer               :: maxpasses,TOTmaxpasses,nPassTypes,iPassType 
integer               :: nLHSODTypes,numpasses,nPrimPL,currentODtype
integer               :: numnodes
Integer,pointer :: Maxpassfortype(:),Myjobs(:)
integer :: iODType,IP,irhstmp,ILHSCOUNT,IODLHS,nthreads,tid
type(integerpointer),pointer :: TypeOverlapIndex(:),PassTypeOverlapIndex(:) 
integer,pointer :: OverlapList(:),ILHSCOUNTINDEX(:)
real(realk) :: maxint,TS,TE
integer(kind=long) :: WORKLENGTH,WORKLENGTH2,WORKLENGTH3,WORK1EST
real(realk)           :: ReductionECONT(input%NDMAT_RHS)

#ifdef VAR_OMP
integer, external :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#endif
IF (IPRINT.GT. 5) THEN
  CALL LSHEADER(LUPRI,'McMurchieDavidson')
  WRITE(LUPRI,'(1X,A)') ''
  WRITE(LUPRI,'(1X,A)') '******************************************************************'
  WRITE(LUPRI,'(1X,A)') '***                     McMurchieDavidson                      ***'
  WRITE(LUPRI,'(1X,A)') '******************************************************************'
  WRITE(LUPRI,'(1X,A)') ''
ENDIF
CALL initTUVitem(sharedTUV,Input,OD_LHS,OD_RHS,LUPRI,IPRINT)
integral%TUV => sharedTUV
CALL Allocitem_init(Allocations,'Both')
CALL Allocitem_zero(Allocations,'Both')

! The LHS batches

call mem_alloc(ODTypeIndex,OD_LHS%nbatches)
call SelectODTypesFromODbatch(ODtypeIndex,Allocations,OD_LHS,OD_LHS%nBatches,&
     &nLHSODTypes,TYPEOVERLAPINDEX,INPUT,IPRINT,LUPRI,'LHS')

IF(IPRINT.GT.5)call Allocitem_print(Allocations,'LHS',lupri)

call mem_alloc(ILHSCOUNTINDEX,OD_LHS%nbatches)
ILHSCOUNT=0
DO iODType=1,nLHSODTypes
   DO IODLHS = 1,TypeOverlapIndex(iODtype)%dim
      ILHS = TypeOverlapIndex(iODtype)%elms(IODLHS)
      ILHSCOUNT=ILHSCOUNT+1
      ILHSCOUNTINDEX(ILHSCOUNT) = ILHS
   ENDDO
ENDDO

!The RHS Batches

call mem_alloc(ODpassesIndex,OD_RHS%nbatches)
call mem_alloc(Q,OD_RHS%nbatches)
call INIT_BUFCOUNTERS(2)
call MEM_OVERLAP(OD_RHS,Input,2)
call ALLOC_ODRHS_BUFFERS
DO IRHS=1,OD_RHS%nbatches
   CALL SET_Overlap(Q(IRHS),Input,SharedTUV,Integral,OD_RHS%BATCH(IRHS),2,LUPRI,IPRINT,.TRUE.)
ENDDO
CALL SelectODPassTypes2(ODpassesIndex,PassTypeOverlapIndex,Q,OD_RHS%nbatches,nPassTypes,Allocations,input,IPRINT,LUPRI)
IF(IPRINT.GT.5)call Allocitem_print(Allocations,'RHS',lupri)
call mem_alloc(nOverlapOfPassType,nPassTypes)
DO Ipasstype = 1,nPassTypes
   nOverlapOfPassType(Ipasstype)=0
ENDDO
DO IRHS=1,OD_RHS%nbatches
   nOverlapOfPassType(ODpassesIndex(IRHS))=nOverlapOfPassType(ODpassesIndex(IRHS))+1
ENDDO
ReductionEcont = 0.0E0_realk
!comment to the below OMP statement: omp is turned off if FMM is to be used in a gradient run. The reason 
!          for this is that the overlap index in the derivative and non-derivative case have to match.
!          A possible solution to this unwanted feature is to calculate both moments and derivative moments at 
!          the same time. 
IF(.NOT.INPUT%noOMP)call mem_TurnONThread_Memory()
!$OMP PARALLEL IF((.NOT.INPUT%noOMP).AND.(.NOT.(INPUT%DO_MULMOM .AND. INPUT%DO_MMGRD))) &
!$OMP DEFAULT(none) PRIVATE(maxint,tid,nthreads,integral,ILHSCOUNT,currentODtype,&
!$OMP iODtype,P,PQ,ILHS,IRHS,iPassType,numpasses,Start_RHS,End_RHS,PassQ,screen,&
!$OMP maxPassesFortypes,TOTmaxpasses,WORKLENGTH,WORK1EST,WORKLENGTH2,WORKLENGTH3,&
!$OMP OverlapList,maxpasses,doPasses,TS,TE,Myjobs,&
!$OMP IRHSTMP) SHARED(sharedTUV,Allocations,Input,OD_LHS,OD_RHS,nOverlapOfPassType,PassTypeOverlapIndex,&
#ifdef VAR_MPI
!$OMP infpar,&
#endif
!$OMP Q,ODpassesIndex,OUTPUT,typeoverlapindex,nPassTypes,nLHSODtypes,IPRINT,LUPRI,ILHSCOUNTINDEX,&
!$OMP ODTypeIndex,ReductionECONT)
IF(.NOT.INPUT%noOMP)call init_threadmemvar()
#ifdef VAR_OMP
nthreads=OMP_GET_NUM_THREADS()
tid = omp_get_thread_num()
#else
nthreads=1
tid=0
#endif
IF(output%FullAlphaCD) call wrapInitThermiteIntThreadID(tid)
call INIT_BUFCOUNTERS(1)
DO iODType=1,nLHSODTypes
   call MEM_INIT_OVERLAP(Allocations,iODtype,Input,1,IPRINT,LUPRI)
ENDDO
call ALLOC_ODLHS_BUFFERS
!reassociate pointers as the may be deassociated inside a OMP PARALLEL REGION
integral%TUV => sharedTUV
IF(INPUT%LHSSameAsRHSDmat)THEN
   INPUT%LST_DLHS => INPUT%LST_DRHS
ENDIF
call mem_alloc(maxPassesFortypes,nPassTypes)

!initial allocation
ILHSCOUNT=MIN(1+tid,OD_LHS%nbatches)
ILHS = ILHSCOUNTINDEX(ILHSCOUNT)
iODType=ODTypeIndex(ILHS)
currentODtype = iODtype
CALL INIT_OVERLAP(P,Allocations,iODtype,Input,1,IPRINT,LUPRI)
CALL allocIntegralsWRAP(PQ,Integral,Input,Allocations,iODtype,nPassTypes,&
     &maxPassesFortypes,1,INPUT%DO_PASSES,nOverlapOfPassType,lupri)
call mem_alloc(Integral%Econt,input%NDMAT_RHS)
IF(INPUT%fullcontraction)Integral%Econt=0.0E0_realk
TOTmaxpasses = 0
call INIT_BUFCOUNTERS(3)
DO iPassType=1,nPassTypes
   TOTmaxpasses = MAX(TOTmaxpasses,maxPassesForTypes(iPassType))
   IRHS = PassTypeOverlapIndex(iPassType)%elms(1)    
   CALL MEM_PASS_FROM_OVERLAP(Q(IRHS),Input,maxPassesForTypes(iPassType))
ENDDO
call ALLOC_ODPASS_BUFFERS
call mem_alloc(PassQ,nPassTypes)

TOTmaxpasses = 0
DO iPassType=1,nPassTypes
   TOTmaxpasses = MAX(TOTmaxpasses,maxPassesForTypes(iPassType))
   IRHS = PassTypeOverlapIndex(iPassType)%elms(1) 
   CALL INIT_PASS_FROM_OVERLAP(PassQ(iPassType),Q(IRHS),Input,maxPassesForTypes(iPassType),IPRINT,LUPRI)
   CALL SetPassOrbitals(PassQ(iPassType),Q(IRHS),maxPassesForTypes(iPassType),LUPRI,IPRINT)
ENDDO

!alloc workarray
!LHS requirements
WORKLENGTH = 0
DO iODType=1,nLHSODtypes
   WORK1EST = Allocations%maxprimLHSA(iODtype)*Allocations%maxContLHSA(iODtype)&
        &+2*Allocations%maxETUVlenLHSA(iODtype)+Allocations%maxijkLHSA(iODtype)*Allocations%maxijkLHSA(iODtype)
   WORKLENGTH = MAX(WORKLENGTH,WORK1EST)
ENDDO
!intermediate
WORKLENGTH = WORKLENGTH+Allocations%maxprimRHS*Allocations%maxprimLHS
!RHS requirements
WORKLENGTH2 = 0
DO iODType=1,nLHSODtypes
   CALL determineMaxPassesForType(Allocations,iODtype,nPassTypes,&
        &maxPassesFortypes,1,INPUT%DO_PASSES,nOverlapOfPassType,lupri)
   DO iPassType=1,nPassTypes
      WORK1EST = Allocations%maxprimRHSA(iPassType)*Allocations%maxContRHSA(iPassType)&
           &+2*Allocations%maxETUVlenRHSA(iPassType)*maxPassesForTypes(iPassType)&
           &+Allocations%maxijkRHSA(iPassType)*Allocations%maxijkRHSA(iPasstype)
      WORKLENGTH2 = MAX(WORKLENGTH2,WORK1EST)
      WORK1EST = Allocations%maxprimRHSA(iPassType)*Allocations%maxprimLHSA(iODtype)&
           &*maxPassesForTypes(iPassType)
      WORKLENGTH2 = MAX(WORKLENGTH2,WORK1EST)
   ENDDO
ENDDO
WORKLENGTH3 = MAX(WORKLENGTH2,WORKLENGTH)
!intermediate
WORK1EST = 5*Allocations%maxPrimRHS
WORKLENGTH3 = MAX(WORKLENGTH3,WORK1EST)
call init_workmem(WORKLENGTH3)

CALL determineMaxPassesForType(Allocations,currentODtype,nPassTypes,&
     & maxPassesFortypes,1,INPUT%DO_PASSES,nOverlapOfPassType,lupri)
iODtype = currentODtype

call mem_alloc(overlaplist,TOTmaxpasses)
#ifdef VAR_MPI
IF(output%FullAlphaCD)THEN
   call mem_alloc(Myjobs,OD_LHS%nbatches)
   call BuildMyjobsList(infpar%lg_nodtot,Myjobs,OD_LHS,Output,infpar%lg_mynum,&
        & OD_LHS%nbatches)
ENDIF
#endif
!$OMP DO SCHEDULE(DYNAMIC,1)
DO ILHSCOUNT=1,OD_LHS%nbatches
   ILHS = ILHSCOUNTINDEX(ILHSCOUNT)
#ifdef VAR_MPI
   IF(output%FullAlphaCD)THEN
      IF(Myjobs(ILHS).NE.infpar%lg_mynum)CYCLE
   ENDIF
#endif
   iODType=ODTypeIndex(ILHS)
   IF(iODType.NE.currentODtype)THEN
      !dealloc
      CALL FREE_OVERLAP(P)
      CALL REINIT_OD_LHS
      CALL deallocIntegrals(PQ,Integral)
      DO iPassType=1,nPassTypes
         CALL Free_overlap(PassQ(iPassType))
      ENDDO
      call mem_dealloc(PassQ)
      call mem_dealloc(overlaplist)
      !realloc
      CALL INIT_OVERLAP(P,Allocations,iODtype,Input,1,IPRINT,LUPRI)
      CALL allocIntegralsWRAP(PQ,Integral,Input,Allocations,iODtype,nPassTypes,&
           &maxPassesFortypes,1,INPUT%DO_PASSES,nOverlapOfPassType,lupri)
      !Pass Memory alloc
      call DEALLOC_ODPASS_BUFFERS
      TOTmaxpasses = 0
      DO iPassType=1,nPassTypes
         TOTmaxpasses = MAX(TOTmaxpasses,maxPassesForTypes(iPassType))
         IRHS = PassTypeOverlapIndex(iPassType)%elms(1) 
         CALL MEM_PASS_FROM_OVERLAP(Q(IRHS),Input,maxPassesForTypes(iPassType))
      ENDDO
      call ALLOC_ODPASS_BUFFERS
      !
      call mem_alloc(PassQ,nPassTypes)
      TOTmaxpasses = 0
      DO iPassType=1,nPassTypes
         TOTmaxpasses = MAX(TOTmaxpasses,maxPassesForTypes(iPassType))
         IRHS = PassTypeOverlapIndex(iPassType)%elms(1) 
         CALL INIT_PASS_FROM_OVERLAP(PassQ(iPassType),Q(IRHS),Input,maxPassesForTypes(iPassType),IPRINT,LUPRI)
         CALL SetPassOrbitals(PassQ(iPassType),Q(IRHS),maxPassesForTypes(iPassType),LUPRI,IPRINT)
      ENDDO
      call mem_alloc(overlaplist,TOTmaxpasses)
      currentODtype = iODtype
   ENDIF
   !maybe set Ecoeff LHS at this stage - while we set the RHS Ecoeff according to input? (mod Worklength) 
   CALL SET_Overlap(P,Input,SharedTUV,Integral,OD_LHS%BATCH(ILHS),1,LUPRI,IPRINT,.FALSE.)
   IF(output%FullAlphaCD) call AOtoMO3CenterDEC(P,output)
   CALL Determine_RHS_loop(INPUT,OD_RHS,ILHS,Start_RHS,End_RHS)
   DO iPassType=1,nPassTypes
      numPasses = 0
      maxPasses = maxPassesFortypes(iPassType)
      doPasses = maxPasses.NE. 1
      !********** RHS loop startes here
      DO IRHS = Start_RHS,End_RHS
         IF (ODpassesIndex(IRHS).EQ.iPassType) THEN
            screen = getScreening(P,Q(IRHS),INPUT,LUPRI,IPRINT)
            IF (.NOT.screen) THEN
               IF (IPRINT.GT. 3) WRITE(LUPRI,'(1X,A,I5,A,I5)') 'Overlap distributions P',ILHS,' and Q',IRHS
               IF (doPasses .AND..NOT. ((ILHS.EQ.IRHS).AND.INPUT%sameODs)) THEN
                  numPasses = numPasses + 1
                  IF(numPasses.EQ. 1)CALL SetPassOrbitals(PassQ(iPassType),Q(IRHS),maxPasses,LUPRI,IPRINT)
                  OverlapList(numPasses) = IRHS
                  IF (numPasses.EQ.maxPasses) THEN
                     DO IP=1,maxPasses
                        IRHSTMP = OverlapList(IP)
                        CALL AddOverlapToPass(PassQ(iPassType),Q(IRHSTMP),&
                             & IP,maxPasses,LUPRI,IPRINT)
                     ENDDO
                     CALL FinalizePass(PassQ(iPassType),Q(IRHS),maxPasses,LUPRI,IPRINT)
                     CALL ExplicitIntegrals(Integral,PQ,P,PassQ(iPassType),INPUT,OUTPUT,ILHS,0,&
                          & LUPRI,IPRINT)
                     CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
                     numPasses = 0
                  ENDIF
               ELSE
                  CALL ExplicitIntegrals(Integral,PQ,P,Q(IRHS),INPUT,OUTPUT,ILHS,IRHS,&
                       &                 LUPRI,IPRINT)
                  CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
               ENDIF
            ENDIF! CS-screening
         ENDIF ! Correct passtype
      ENDDO !RHS
      IF (numPasses.GT. 0) THEN
         DO IP=1,numPasses
            IRHS = OverlapList(IP)
            CALL AddOverlapToPass(PassQ(iPassType),Q(IRHS),&
                 & IP,maxPasses,LUPRI,IPRINT)
         ENDDO
         CALL FinalizePass(PassQ(iPassType),Q(IRHS),numPasses,LUPRI,IPRINT)
         CALL ExplicitIntegrals(Integral,PQ,P,PassQ(iPassType),INPUT,OUTPUT,ILHS,IRHS,LUPRI,IPRINT)
         CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
      ENDIF
   ENDDO !Passtypes
ENDDO !ILHSCOUNT
!$OMP END DO NOWAIT

#ifdef VAR_MPI
IF(output%FullAlphaCD) call mem_dealloc(Myjobs)
#endif
CALL FREE_OVERLAP(P)
call DEALLOC_ODLHS_BUFFERS
IF(INPUT%fullcontraction)THEN
!$OMP CRITICAL
 do ILHS=1,input%NDMAT_RHS
  ReductionECONT(ILHS)=ReductionECONT(ILHS) + Integral%Econt(ILHS)
 enddo
!$OMP END CRITICAL
ENDIF
call mem_dealloc(Integral%Econt)
CALL deallocIntegrals(PQ,Integral)
call free_workmem
DO iPassType=1,nPassTypes 
   CALL Free_Overlap(PassQ(iPassType))
ENDDO
call mem_dealloc(PassQ)
call DEALLOC_ODPASS_BUFFERS
call mem_dealloc(OverlapList)

call mem_dealloc(maxPassesFortypes)

IF(.NOT.INPUT%noOMP)call collect_thread_memory()
!$OMP END PARALLEL
IF(.NOT.INPUT%noOMP)call mem_TurnOffThread_Memory()

IF(INPUT%fullcontraction)THEN
 do ILHS=1,input%NDMAT_RHS
  Output%ResultTensor%LSAO(1)%elms(ILHS) = ReductionECONT(ILHS)
 enddo
ENDIF
call DEALLOC_ODRHS_BUFFERS

call mem_dealloc(ODTypeIndex)
call mem_dealloc(ILHSCOUNTINDEX)

call mem_dealloc(nOverlapOfPassType)
DO iODType=1,nLHSODtypes
   call mem_dealloc(TypeOverlapIndex(iODtype)%elms)
ENDDO
DEALLOCATE(TypeOverlapIndex)
DO iPassType=1,nPasstypes
   call mem_dealloc(PassTypeOverlapIndex(iPassType)%elms)
ENDDO
DEALLOCATE(PassTypeOverlapIndex)
call allocitem_free(Allocations,'Both')
DO IRHS=1,OD_RHS%nbatches
   CALL FREE_OVERLAP(Q(IRHS))
ENDDO
call mem_dealloc(Q)
call mem_dealloc(ODpassesIndex)
CALL freeTUVitem(sharedTUV,Input)

END SUBROUTINE McMurchieDavidson

!> \brief driver routine to calculate the explicit integral for a given P and Q overlap
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!>  This is a routine which calculated the explicit 4 center integral from the 
!>  overlap P and Q given in input. It first build the integrand structure PQ
!>  which contains stuff like the reduced exponent, the integral prefactor, etc. 
!>  It then calculates the integral over zero'th order primitive hermite functions 
!>  Eq. 9.9.14 (and for instance 9.5.41) in the book. and performs 
!>  recurrence relations (see for instance Eq. 9.9.18-9.9.20) in order to obtain the full     
!>  integral over primitive hermite functions. 
!>  Finally it does the 
!>     Ecoefficient contraction  (see Eq. 9.5.1)
!>     Primitive Contraction (from primitve AO integrals to contracted AO integrals)    
!>     Spherical transformation (from cartesian to spherical harmonic AOs)
!>  both for electron 2 (contract_Q) and electron 1 (contract_P)
!>
!>  \param Integral contains arrays to store intermidiates and final integrals
!>  \param PQ contains integrand info like reduced exponents, integral prefactor, etc.
!>  \param P contains overlap distribution for the left hand side, electron 1
!>  \param Q contains overlap distribution for the right hand side, electron 2
!>  \param INPUT the integral input specifications, contains all info about what to do
!>  \param OUTPUT the integral output specifications, determines how the output should be given
!>  \param ILHS the index for which LHS overlap distribution this is, relevant when ILHS = IRHS  
!>  \param IRHS the index for which RHS overlap distribution this is, relevant when ILHS = IRHS  
!>  \param LUPRI the logical unit number for the output file
!>  \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE ExplicitIntegrals(Integral,PQ,P,Q,INPUT,OUTPUT,ILHS,IRHS,&
     &                       LUPRI,IPRINT)
implicit none
Integer              :: ILHS,IRHS,LUPRI,IPRINT
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(INTEGRALOUTPUT) :: OUTPUT
TYPE(Integralitem)   :: Integral
TYPE(Overlap)        :: P,Q
TYPE(Integrand)      :: PQ
REAL(REALK)          :: CPU1,CPU2,WALL1,WALL2

LOGICAL             :: DO_INTEREST,error,P_inter, Q_inter
Integer             :: nPQ,la,lb,lc,ld,iPQ
Real(realk)         :: A(3),B(3),C(3),D(3),an,bn,cn,dn,ea,eb,ec,ed
Real(realk),pointer :: inter(:)
Integer             :: derOrder, startDer, endDer,iDerLHS
#ifdef MOD_UNRELEASED
#if VAR_INTEREST
IF (INPUT%interest) THEN
  P_inter = (.NOT.(P%TYPE_empty.OR.P%single)).AND.(P%nPrimitives.EQ. 1).AND.(P%nAngmom.EQ. 1)&
     &      .AND.(P%endGeoOrder.EQ. 0).AND.(P%orbital1%angmom(1).GE.P%orbital2%angmom(1)).AND.&
     &      (P%nPasses.EQ. 1)
  Q_inter = (.NOT.(Q%TYPE_empty.OR.Q%single)).AND.(Q%nPrimitives.EQ. 1).AND.(Q%nAngmom.EQ. 1)&
     &      .AND.(Q%endGeoOrder.EQ. 0).AND.(Q%orbital1%angmom(1).GE.Q%orbital2%angmom(1)).AND.&
     &      (Q%nPasses.EQ. 1)
  DO_INTEREST = (INPUT%operator.EQ.CoulombOperator).AND.P_inter.AND.Q_inter
ELSE
  DO_INTEREST = .FALSE.
ENDIF
#endif
#endif
IF(INPUT%DO_PROP)THEN
  call ExplicitPropIntegrals(Integral,PQ,P,Q,INPUT,OUTPUT,ILHS,IRHS,&
     &                       LUPRI,IPRINT)
#ifdef MOD_UNRELEASED
#if VAR_INTEREST
ELSEIF (DO_INTEREST) THEN
  CALL Build_Integrand(PQ,P,Q,INPUT,ILHS,IRHS,LUPRI,IPRINT)
  la = P%orbital1%angmom(1) + 1
  lb = P%orbital2%angmom(1) + 1
  lc = Q%orbital1%angmom(1) + 1
  ld = Q%orbital2%angmom(1) + 1
  ea = P%orbital1%exponents(1) 
  eb = P%orbital2%exponents(1) 
  ec = Q%orbital1%exponents(1) 
  ed = Q%orbital2%exponents(1) 
  A  = P%orbital1%center
  B  = P%orbital2%center
  C  = Q%orbital1%center
  D  = Q%orbital2%center
  an = P%orbital1%CC(1)%p%elms(1)
  bn = P%orbital2%CC(1)%p%elms(1)
  cn = Q%orbital1%CC(1)%p%elms(1)
  dn = Q%orbital2%CC(1)%p%elms(1)
  nPQ = 1
  call interest_initialize()
  call interest_eri(INTEGRAL%integralsABCD,nPQ,la,ea,A(1),A(2),A(3),an,lb,eb,B(1),B(2),B(3),bn,&
     &              lc,ec,C(1),C(2),C(3),cn,ld,ed,D(1),D(2),D(3),dn,.false.)
!  call mem_alloc(inter,nPQ)
!  inter = INTEGRAL%integralsABCD(1:nPQ)
!  CALL Build_Integrand(PQ,P,Q,INPUT,ILHS,IRHS,LUPRI,IPRINT)
!  !   Hermite 2-electron integral over general operator w
!  CALL Build_HermiteTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT)
!  CALL Contract_Q(INTEGRAL,PQ,Input,LUPRI,IPRINT)
!  CALL Contract_P(INTEGRAL,PQ,input,LUPRI,IPRINT)
!  error = .FALSE.
!  DO iPQ=1,nPQ
!     IF (ABS(inter(iPQ) - INTEGRAL%integralsABCD(iPQ)).GT. 1E-9_realk) THEN
!        print*,'ABS(inter(iPQ) - INTEGRAL%integralsABCD(iPQ))',ABS(inter(iPQ) - INTEGRAL%integralsABCD(iPQ))
!        error = .TRUE.
!     ENDIF
!  ENDDO
!  IF (error) THEN
!     write(*,*) 'dbg:nPQ ',nPQ
!     write(*,*) 'dbg:    ',P%totOrbitals,Q%totOrbitals
!     write(*,*) 'dbg:a   ',la,ea,A(1),A(2),A(3),an
!     write(*,*) 'dbg:b   ',lb,eb,b(1),b(2),b(3),bn
!     write(*,*) 'dbg:c   ',lc,ec,c(1),c(2),c(3),cn
!     write(*,*) 'dbg:d   ',ld,ed,d(1),d(2),d(3),dn
!     write(*,*) 'dbg:inter'
!     write(*,'(6ES15.6)') inter(1:nPQ)
!     write(*,*) 'dbg:therm'
!     write(*,'(6ES15.6)') INTEGRAL%integralsABCD(1:nPQ)
!     call lsquit('error',-1)
!  ENDIF
!  call mem_dealloc(inter)
#endif
#endif
ELSE
!  Settings for derivative loop structure
   call getInputDerivativeInfo(derOrder,startDer,endDer,input,P%TYPE_Empty,Q%TYPE_Empty)
   Integral%startDerivativeIndex = 0
   Integral%nDerivComp = input%NGEODERIVCOMP
   Integral%dohodi = startDer.NE.endDer
   
   CALL Build_Integrand(PQ,P,Q,INPUT,ILHS,IRHS,LUPRI,IPRINT)
   !   Hermite 2-electron integral over general operator w
   CALL Build_HermiteTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT)
   DO iDerLHS = endDer,startDer,-1
     call setIntegralDerivativOrders(Integral,iDerLHS,derOrder,P%single,Q%single)
     CALL Contract_Q(INTEGRAL,PQ,Input,LUPRI,IPRINT)
     CALL Contract_P(INTEGRAL,PQ,input,LUPRI,IPRINT)
     Integral%startDerivativeIndex = Integral%startDerivativeIndex + Integral%lhsGeoComp*Integral%rhsGeoComp
   ENDDO
ENDIF

END SUBROUTINE ExplicitIntegrals

!> \brief driver routine to calculate the explicit one electron property integrals
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!>  This is a routine which calculated the explicit one electron property integral from the 
!>  overlap P given in input. 
!>  The code is to a large extend based on the abacus/her1int routines originally written 
!>  by Kenneth Ruud and Trygve Helgaker. 
!>
!>  \param Integral contains arrays to store intermidiates and final integrals
!>  \param PQ contains integrand info like reduced exponents, integral prefactor, etc.
!>  \param P contains overlap distribution for the left hand side, electron 1
!>  \param Q contains overlap distribution for the right hand side, electron 2
!>  \param INPUT the integral input specifications, contains all info about what to do
!>  \param OUTPUT the integral output specifications, determines how the output should be given
!>  \param ILHS the index for which LHS overlap distribution this is, relevant when ILHS = IRHS  
!>  \param IRHS the index for which RHS overlap distribution this is, relevant when ILHS = IRHS  
!>  \param LUPRI the logical unit number for the output file
!>  \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE ExplicitPropIntegrals(Integral,PQ,P,Q,INPUT,OUTPUT,ILHS,IRHS,&
     &                       LUPRI,IPRINT)
implicit none
Integer              :: ILHS,IRHS,LUPRI,IPRINT
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(INTEGRALOUTPUT) :: OUTPUT
TYPE(Integralitem)   :: Integral
TYPE(Overlap)        :: P,Q
TYPE(Integrand)      :: PQ
REAL(REALK)          :: CPU1,CPU2,WALL1,WALL2

Integer             :: JmaxA,JmaxB,JmaxP,JMAXAP,JMAXBP,JMAXTP,JMAXD,JMAXM,JMAXT
Integer             :: ITEX,iPrimP,iAngmomP,nPrimP,proptype,ijkP,ijk1,ijk2,lmP
integer             :: nOperatorComp,J,nPrimPQ,jmaxPQ,I,U,V,tuv,T,IX,QnPasses
integer             :: iprimA,iprimB,offsetTUV,offsetTUV2,iPassQ,offset
REAL(REALK)         :: SAAB13,FAC,TEXPB1,FACB,TMPB
REAL(REALK)         :: ORIGIN(3),Apreexpfac
LOGICAL             :: DIFODC,KINODC,ONECEN,DONUC1,DOMOM1
Real(realk),pointer :: ODC(:,:,:,:,:,:,:),EXPA(:),EXPB(:),SHGTF(:)
Real(realk),pointer :: PADIST(:,:),PBDIST(:,:),PODIST(:,:),preexpfacinv3(:),exppi(:)
real(realk),parameter :: inv3=1.0E0_realk/3.0E0_realk,D1=1.0E0_realk
real(realk) :: AODIST(3),BODIST(3),CPX,CPY,CPZ,Z
real(realk),pointer   :: ptemp(:),R2(:),alpha(:),prefactor(:)

PROPTYPE = INPUT%PROPTYPE
DIFODC = INPUT%PropDerivEcoeff
KINODC = INPUT%PropKineticEcoeff
DOMOM1 = INPUT%PropMomentEcoeff
!WRITE(lupri,*)'PROPTYPE  = ',PROPTYPE
IF(PROPTYPE.EQ.17.OR.PROPTYPE.EQ.20.OR.PROPTYPE.EQ.26)THEN
   !Center of nuclei B
   ORIGIN = P%orbital2%center
ELSE
   ORIGIN = INPUT%PROP_ORIGIN
ENDIF
!WRITE(lupri,*)'ORIGIN',ORIGIN
!       CALL PRINT_OVERLAP(P,LUPRI,1000,'LHS')
!       CALL PRINT_OVERLAP(Q,LUPRI,1000,'RHS')
CALL Build_Integrand(PQ,P,Q,INPUT,ILHS,IRHS,LUPRI,IPRINT)
nPrimP = P%nPrimitives
IF(INPUT%PropRequireBoys.GT.-1)THEN
 IF(PQ%Operator.EQ.NucpotOperator)THEN
    NPrimPQ=PQ%nPrimitives
    JMAXPQ=PQ%endAngmom+INPUT%PropRequireBoys
    CALL buildNuclearRJ000(Integral%IN,PQ,nPrimPQ,JMAXPQ,LUPRI,IPRINT,integral,INPUT%HIGH_RJ000_ACCURACY)
    IF (IPRINT .GE. 20) THEN
       CALL LSHEADER(LUPRI,'Output from W000')
       DO I=1,NPrimPQ
          DO J=0,JmaxPQ
             WRITE(LUPRI,'(2X,A6,I4,A1,I4,A2,ES16.8)')'SJ000(',J,',',I,')=',Integral%IN(1+J+(I-1)*(JmaxPQ+1))
          ENDDO
       ENDDO
    ENDIF
    INTEGRAL%nTUV=(JMAXPQ+1)*(JMAXPQ+2)*(JMAXPQ+3)/6
    Integral%nEFG=1 !no cartesian multipole moments
    CALL WTUVrecurrence(Integral%Rtuv,Integral%IN,Integral%TUV,PQ%distance,&
         &              0,jMaxPQ,nPrimPQ,integral%ntuv,lupri,iprint)
    IF (IPRINT .GE. 20) THEN
       CALL LSHEADER(LUPRI,'Output from WTUVrecurrence')
       WRITE (LUPRI,'(2X,A,I10)') 'JMAX  ', JMAXPQ
       WRITE (LUPRI,'(2X,A,I10)') 'NPrim  ', NPrimPQ
       WRITE (LUPRI,'(2X,A,I10)') 'NTUV  ', integral%nTUV
       CALL LSHEADER(LUPRI,'Hermite integrals S(t,u,v)')
       DO J = 0, JMAXPQ
          DO T = J,0,-1
             DO U = J-T,0,-1
                V=J-T-U
                TUV=Integral%TUV%tuvIndex(T,U,V)
                WRITE (LUPRI,'(2X,A2,I3,A1,I3,A1,I3,A1,2X,5ES16.8/,(18X,5ES16.8))')&
                     & 'W(',T,',',U,',',V,')', &
                     &(INTEGRAL%Rtuv(I + (TUV-1)*nPrimPQ),I=1,NPrimPQ)
                WRITE (LUPRI,*) ' '
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    IF(Input%addtointegral.AND.Q%orbital1%TYPE_Nucleus)THEN
       QnPasses = Q%nPasses 
       ! We sum all the charges in Overlap PassQ (all same charge)
       call mem_alloc(ptemp,nPrimPQ*Integral%nTUV )
       Z = -Q%orbital1%CC(1)%p%elms(1) !Charge
       DO TUV = 1,Integral%nTUV 
          offsetTUV=(TUV-1)*nPrimPQ
          DO I=1,NPrimPQ
             ptemp(I+offsetTUV)= Z*INTEGRAL%Rtuv(I+offsetTUV)
          ENDDO
       ENDDO
       DO TUV = 1,Integral%nTUV 
          offsetTUV=(TUV-1)*nPrimP
          DO I=1,NPrimP
             iPassQ=1
             offsetTUV2=(I-1)*Q%nPasses + (TUV-1)*nPrimPQ
             INTEGRAL%Rtuv(I+offsetTUV)= ptemp(iPassQ + offsetTUV2) 
             DO iPassQ=2,QnPasses
                INTEGRAL%Rtuv(I+offsetTUV)= INTEGRAL%Rtuv(I+offsetTUV)+ptemp(iPassQ + offsetTUV2) 
             ENDDO
          ENDDO
       ENDDO
       call mem_dealloc(ptemp)
       Q%nPasses = 1
    ENDIF
 ENDIF
ENDIF
JMAXD = INPUT%PROPMAXD
JMAXM = INPUT%PROPMAXM
DO iAngmomP=1,P%nAngmom
   JMAXA = P%orbital1%angmom(P%indexAng1(iAngmomP))
   IF (PROPTYPE .EQ. 42) JMAXA = JMAXA + 1
   IF (PROPTYPE .EQ. 43) JMAXA = JMAXA + 1
   JMAXB = P%orbital2%angmom(P%indexAng2(iAngmomP))

   JMAXT = JMAXA + JMAXB + JMAXD + JMAXM

!   WRITE(lupri,*)'JMAXA',JMAXA
!   WRITE(lupri,*)'JMAXB',JMAXB
!   WRITE(lupri,*)'JMAXD',JMAXD
!   WRITE(lupri,*)'JMAXM',JMAXM
!   WRITE(lupri,*)'JMAXT',JMAXT

   JMAXAP = JMAXA
   JMAXBP = JMAXB + JMAXD
   IF(DOMOM1)THEN
      JMAXAP = JMAXAP + JMAXM
   ENDIF
   ITEX = MAX(1,JMAXD + JMAXM)
   JMAXTP = JMAXT + ITEX

!   IF ((PROPTYPE .EQ. 49) .OR. (PROPTYPE .EQ. 51)) THEN
!      IF (DOMOM1) THEN
!         JMAXBP = JMAXBP + 1
!      ELSE
!         JMAXAP = JMAXAP + 1
!      END IF
!      JMAXTP = JMAXTP +1
!   ENDIF

!   WRITE(lupri,*)'JMAXAP',JMAXAP
!   WRITE(lupri,*)'JMAXBP',JMAXBP
!   WRITE(lupri,*)'JMAXTP',JMAXTP
   nPrimP = P%nPrimitives
   call mem_alloc(PADIST,nPRimP,3)
   call mem_alloc(PBDIST,nPRimP,3)
   call mem_alloc(preexpfacinv3,nPRimP)
   call mem_alloc(exppi,nPRimP)
   DO IX=1,3
      DO iPrimP=1,nPrimP
         PADIST(iPrimP,IX) = P%center(IX+(iPrimP-1)*3)-P%orbital1%center(IX)
         PBDIST(iPrimP,IX) = P%center(IX+(iPrimP-1)*3)-P%orbital2%center(IX)
      ENDDO
   ENDDO
#ifdef VAR_LSDEBUGINT
   IF(size(P%preexpfac).NE.P%nPrimAlloc*P%nAngAlloc)call lsquit('preexpfac dim error',-1)
   IF(nPrimP.GT.P%nPrimAlloc)call lsquit('preexpfac dim error',-1)
#endif
   offset = (iAngmomP-1)*P%nPrimAlloc
   DO iPrimP=1,nPrimP
      Apreexpfac = ABS(P%preexpfac(iPrimP+offset))
      preexpfacinv3(iPRimP) = SIGN(Apreexpfac**(inv3),P%preexpfac(iPrimP+offset)) 
   ENDDO
!#ifdef VAR_MKL
!   call vdinv(nPrimP,P%exponents,EXPPI)
!#else
   DO iPrimP=1,nPrimP
      EXPPI(iPrimP) = 1/(P%exponents(iPrimP))
   ENDDO
!#endif
   call mem_alloc(expA,nPrimP)
   call mem_alloc(expB,nPrimP)
   iPrimP = 0
   DO iPrimB=1,P%orbital2%nPrimitives
      TMPB = P%orbital2%exponents(iPrimB)
      DO iPrimA=1,P%orbital1%nPrimitives
         iPrimP = iPrimP+1
         EXPA(iPrimP) = P%orbital1%exponents(iPrimA)
         EXPB(iPrimP) = TMPB 
      ENDDO
   ENDDO
   
   !npasses!
   AODIST(1) = P%orbital1%center(1)-ORIGIN(1)
   AODIST(2) = P%orbital1%center(2)-ORIGIN(2)
   AODIST(3) = P%orbital1%center(3)-ORIGIN(3)
   BODIST(1) = P%orbital2%center(1)-ORIGIN(1)
   BODIST(2) = P%orbital2%center(2)-ORIGIN(2)
   BODIST(3) = P%orbital2%center(3)-ORIGIN(3)

   call mem_alloc(ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,3,&
        &         .FALSE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.)
   
!        **********************************************
!        ***** Build of Ecoefficients             *****
!        **********************************************

   call buildpropEcoeff(ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,&
        & JMAXAP,JMAXBP,JMAXTP,DIFODC,KINODC,EXPA,EXPB,&
        & preexpfacinv3,EXPPI,PADIST,PBDIST,DOMOM1,&
        & ORIGIN,PROPTYPE,AODIST,BODIST,IPRINT,LUPRI)

   call mem_dealloc(PADIST)
   call mem_dealloc(PBDIST)
   call mem_dealloc(preexpfacinv3)
   call mem_dealloc(expA)
   call mem_dealloc(expB)


!#ifdef VAR_MKL
!   call mem_alloc(SHGTF,nPrimP)
!   call vdsqrt(nprimP,EXPPI,SHGTF)
!   DO iPrimP=1,nPrimP
!      SHGTF(iPrimP) = SQRTPI*SHGTF(iPrimP)
!   ENDDO
!   call mem_dealloc(exppi)
!#else
   call mem_alloc(SHGTF,nPrimP)
   DO iPrimP=1,nPrimP
      SHGTF(iPrimP) = SQRTPI*SQRT(EXPPI(iPrimP))
   ENDDO
   call mem_dealloc(exppi)
!#endif
!        **********************************************
!        ***** Calculation of Hermitian integrals *****
!        **********************************************

   CALL GET_IJK(JMAXA,JMAXB,ijk1,ijk2,lmP,ijkP,.FALSE.,0,.FALSE.)
   IF(PROPTYPE.EQ.1)THEN  
!     OVERLAP
      nOperatorComp = 1
      IF (P%segmented) THEN      
         CALL OVERLAPINTseg(Integral%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,&
              & JMAXM,SHGTF,JMAXA,JMAXB,ijkP,lupri)
         integral%nPrim = 1
      ELSE
         CALL OVERLAPINTgen(Integral%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,&
              & JMAXM,SHGTF,JMAXA,JMAXB,ijkP,lupri)
         integral%nPrim = nPrimP
      ENDIF
   ELSEIF(PROPTYPE.EQ.2)THEN  
!     DIPLEN
      nOperatorComp = 3

      call mem_alloc(PODIST,nPRimP,3)
      DO iPrimP=1,nPrimP
         PODIST(iPrimP,1) = P%center(1+(iPrimP-1)*3)-ORIGIN(1)
         PODIST(iPrimP,2) = P%center(2+(iPrimP-1)*3)-ORIGIN(2)
         PODIST(iPrimP,3) = P%center(3+(iPrimP-1)*3)-ORIGIN(3)
      ENDDO
      IF (P%segmented) THEN      
         CALL DIPLENINTseg(Integral%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,&
              & JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,POdist,lupri)
         integral%nPrim = 1
      ELSE
         CALL DIPLENINTgen(Integral%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,&
              & JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,POdist,lupri)
         integral%nPrim = nPrimP
      ENDIF
      call mem_dealloc(PODIST)
   ELSEIF(PROPTYPE.EQ.3)THEN  
!     DIPVEL
      nOperatorComp = 3
      IF (P%segmented) THEN      
         CALL DIPVELINTseg(Integral%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,&
              & JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,lupri)
         integral%nPrim = 1
      ELSE
         CALL DIPVELINTgen(Integral%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,&
              & JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,lupri)
         integral%nPrim = nPrimP
      ENDIF
   ELSEIF(PROPTYPE.EQ.7)THEN  
!     THETA
      nOperatorComp = 6
      IF (P%segmented) THEN      
         CALL THETAINTseg(Integral%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,&
              & JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,lupri)
         integral%nPrim = 1
      ELSE
         CALL THETAINTgen(Integral%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,&
              & JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,lupri)
         integral%nPrim = nPrimP
      ENDIF
   ELSEIF(PROPTYPE.EQ.10)THEN
      nOperatorComp = 3
      !NPrimPQ=Q%nPasses*nPrimP
      IF (P%segmented) THEN   
         call PSOINTseg(INTEGRAL%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,&
              & JMAXD,JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,&
              & integral%RTUV,Integral%nTUV,NPrimPQ,Q%nPasses ,integral%TUV,lupri)
         integral%nPrim = 1
      ELSE
         call PSOINTgen(INTEGRAL%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,&
              & JMAXD,JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,&
              & integral%RTUV,Integral%nTUV,NPrimPQ,Q%nPasses ,integral%TUV,lupri)
         integral%nPrim = nPrimP
      ENDIF
      nOperatorComp = 3*Q%nPasses 
   ELSEIF(PROPTYPE.EQ.17.OR.PROPTYPE.EQ.18)THEN  
!     ANGLON(17) OR  ANGMOM(18)
      nOperatorComp = 3
      IF (P%segmented) THEN      
         CALL ANGMOMINTseg(Integral%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,&
              & JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,lupri)
         integral%nPrim = 1
      ELSE
         CALL ANGMOMINTgen(Integral%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,&
              & JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,lupri)
         integral%nPrim = nPrimP
      ENDIF
   ELSEIF(PROPTYPE.EQ.19)THEN
      nOperatorComp = 3
      IF (P%segmented) THEN   
         CALL LONMOM1INTseg(Integral%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,&
              & JMAXD,JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,&
              & P%distance12,lupri)
         integral%nPrim = 1
      ELSE
         CALL LONMOM1INTgen(Integral%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,&
              & JMAXD,JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,&
              & P%distance12,lupri)
         integral%nPrim = nPrimP
      ENDIF
   ELSEIF(PROPTYPE.EQ.26)THEN
      nOperatorComp = 9
      IF (P%segmented) THEN   
         call NSTNOLINTseg(INTEGRAL%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,&
              & JMAXD,JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,&
              & integral%RTUV,Integral%nTUV,NPrimPQ,Q%nPasses ,&
              & integral%TUV,P%distance12,lupri)
         integral%nPrim = 1
      ELSE
         call NSTNOLINTgen(INTEGRAL%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,&
              & JMAXD,JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,&
              & integral%RTUV,Integral%nTUV,NPrimPQ,Q%nPasses ,&
              & integral%TUV,P%distance12,lupri)
         integral%nPrim = nPrimP
      ENDIF
      nOperatorComp = 9*Q%nPasses 
   ELSEIF(PROPTYPE.EQ.27)THEN
      nOperatorComp = 9
      !NPrimPQ=Q%nPasses*nPrimP
      IF (P%segmented) THEN   
         call NSTLONINTseg(INTEGRAL%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,&
              & JMAXD,JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,&
              & integral%RTUV,Integral%nTUV,NPrimPQ,Q%nPasses ,&
              & integral%TUV,P%distance12,lupri)
         integral%nPrim = 1
      ELSE
         call NSTLONINTgen(INTEGRAL%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,&
              & JMAXD,JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,&
              & integral%RTUV,Integral%nTUV,NPrimPQ,Q%nPasses ,&
              & integral%TUV,P%distance12,lupri)
         integral%nPrim = nPrimP
      ENDIF
      nOperatorComp = 9*Q%nPasses 
   ELSEIF(PROPTYPE.EQ.42)THEN
      nOperatorComp = 9
      !D-CM1 X
      IF (P%segmented) THEN
         call DCM1INTseg(INTEGRAL%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,&
              & JMAXD,JMAXM,SHGTF,JMAXA-1,JMAXB,ijkP,nOperatorComp,&
              & P%distance12,AODIST,lupri)
         integral%nPrim = 1
      ELSE
         call DCM1INTgen(INTEGRAL%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,&
              & JMAXD,JMAXM,SHGTF,JMAXA-1,JMAXB,ijkP,nOperatorComp,&
              & P%distance12,AODIST,lupri)
         integral%nPrim = nPrimP
      ENDIF
!ANTI TRUE
   ELSEIF(PROPTYPE.EQ.43)THEN
      nOperatorComp = 18
      !D-CM2 XX,XY,XZ,YY,YZ,ZZ
      IF (P%segmented) THEN
         call DCM2INTseg(INTEGRAL%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,&
              & JMAXD,JMAXM,SHGTF,JMAXA-1,JMAXB,ijkP,nOperatorComp,&
              & P%distance12,AODIST,lupri)
         integral%nPrim = 1
      ELSE
         call DCM2INTgen(INTEGRAL%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,&
              & JMAXD,JMAXM,SHGTF,JMAXA-1,JMAXB,ijkP,nOperatorComp,&
              & P%distance12,AODIST,lupri)
         integral%nPrim = nPrimP
      ENDIF
!ANTI FALSE
   ELSEIF(PROPTYPE.EQ.55)THEN  
!     ROTSTR
      nOperatorComp = 6
      IF (P%segmented) THEN      
         CALL ROTSTRINTseg(Integral%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,&
              & JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,lupri)
         integral%nPrim = 1
      ELSE
         CALL ROTSTRINTgen(Integral%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,JMAXD,&
              & JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,lupri)
         integral%nPrim = nPrimP
      ENDIF
   ELSEIF(PROPTYPE.EQ.190)THEN
      nOperatorComp = 3
      IF (P%segmented) THEN   
         CALL LONMOM2INTseg(Integral%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,&
              & JMAXD,JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,&
              & P%distance12,Integral%Rtuv,integral%nTUV,Integral%TUV,lupri)
         integral%nPrim = 1
      ELSE
         CALL LONMOM2INTgen(Integral%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,&
              & JMAXD,JMAXM,SHGTF,JMAXA,JMAXB,ijkP,nOperatorComp,&
              & P%distance12,Integral%Rtuv,integral%nTUV,Integral%TUV,lupri)
         integral%nPrim = nPrimP
      ENDIF
   ELSEIF(PROPTYPE.EQ.64)THEN
      nOperatorComp = 1
      IF (P%segmented) THEN   
         call POTINTseg(INTEGRAL%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,&
              & JMAXD,JMAXM,SHGTF,JMAXA,JMAXB,ijkP,integral%RTUV,&
              & Integral%nTUV,integral%TUV,lupri)
         integral%nPrim = 1
      ELSE
         call POTINTgen(INTEGRAL%IN,ODC,nPrimP,JMAXA,JMAXB,JMAXT,&
              & JMAXD,JMAXM,SHGTF,JMAXA,JMAXB,ijkP,integral%RTUV,&
              & Integral%nTUV,integral%TUV,lupri)
         integral%nPrim = nPrimP
      ENDIF
   ENDIF
   call mem_dealloc(SHGTF)

   IF(IPRINT.GT.100)THEN
      CALL PrintTensor(INTEGRAL%IN,'PropertyIntegrals.  ',nPrimP,nOperatorComp,&
           & ijkP,lupri,'PrimP ','nOPcom','ijk   ',3)
   ENDIF

   integral%nOrb = nOperatorComp
   integral%nAng = ijkP
   Integral%ngeoDeriv = 1
   Integral%nmagDeriv = 1
   integral%nOperatorComp = nOperatorComp !not used
   
   CALL ContractBasis(Integral,P,iAngmomP,LUPRI,IPRINT)
   !Output(nContP,nPrimQ,ntuvQ,ijkP) = Output(nContP,nOperatorComp,ijkP) 
   CALL SphericalTransform(Integral,P,iAngmomP,LUPRI,IPRINT)
   !Output(nContP,nPrimQ,ntuvQ,lmP) = Output(nContP,nOperatorComp,lmP) 
   CALL AddToP(Integral,P,iAngmomP,nOperatorComp,input%PropAnti,LUPRI,IPRINT)
   call mem_dealloc(ODC)
ENDDO

END SUBROUTINE ExplicitPropIntegrals

SUBROUTINE AddToP(Integral,P,iAngmomP,nOperatorComp,antiAB,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Overlap)      :: P
integer :: iAngmomP,nOperatorComp,LUPRI,IPRINT
logical :: antiAB
!
integer :: nA,nB,nCompA,nCompB,nContA,nContB,startA,startB,nPassP,angA,angB
logical :: permute,Antipermute
angA    = P%indexAng1(iAngmomP)
angB    = P%indexAng2(iAngmomP)
nA      = P%orbital1%totOrbitals
nContA  = P%orbital1%nContracted(angA)
startA  = P%orbital1%startLocOrb(angA)
nCompA  = P%orbital1%nOrbComp(angA)
nB      = P%orbital2%totOrbitals
nContB  = P%orbital2%nContracted(angB)
startB  = P%orbital2%startLocOrb(angB)
nCompB  = P%orbital2%nOrbComp(angB)
permute = P%sameAO.AND.(startA.NE.startB)
Antipermute = antiAB.AND.permute
CALL AddToP1(Integral%integralsABCD,Integral%IN,nA,nB,nCompA,nCompB,nContA,&
   & nContB,startA,startB,nOperatorComp,permute,Antipermute,LUPRI,IPRINT)
END SUBROUTINE ADDTOP

subroutine AddToP1(integralsAB,AddtoP,nA,nB,nCompA,nCompB,nContA,nContB,&
     &startA,startB,nOperatorComp,permute,Antipermute,LUPRI,IPRINT)
implicit none
integer :: nA,nB,nCompA,nCompB,nContA,nContB,startA,startB,nOperatorComp
integer :: lupri,iprint
logical :: permute,Antipermute
Real(realk),intent(IN)   :: AddtoP(nContA,nContB,nOperatorComp,nCompA,nCompB)
Real(realk),intent(INOUT):: integralsAB(nA,nB,nOperatorComp)
!
Integer :: X,iA,iB,iContA,iContB,iAngA,iAngB
 do X=1,nOperatorComp
  iB = startB -1
  DO iContB=1,nContB
   DO iAngB=1,nCompB
    iB = iB + 1
    iA = startA -1
    DO iContA=1,nContA
     DO iAngA=1,nCompA
      iA = iA + 1
      integralsAB(iA,iB,X) = AddtoP(iContA,iContB,X,iAngA,iAngB)
      IF (antipermute) THEN
         integralsAB(iB,iA,X) = -AddtoP(iContA,iContB,X,iAngA,iAngB)
      ELSEIF (permute) THEN
         integralsAB(iB,iA,X) = AddtoP(iContA,iContB,X,iAngA,iAngB)
      ENDIF
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
END SUBROUTINE AddToP1

!> \brief add the overlap P to an overlap, which we call PassP
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param PassP The overlap, to be built
!> \param P the overlap
!> \param numPass the current number of passes in PassP 
!> \param maxpasses the maximum number of passes
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE AddOverlapToPass(PassP,P,numPass,maxPasses,LUPRI,IPRINT)
Implicit none
TYPE(Overlap),intent(inout) :: PassP
TYPE(Overlap),intent(in)    :: P
integer,intent(in)          :: numPass,maxPasses,LUPRI,IPRINT
!
integer :: nprim,nang,na1,na2,iangmom,startPrim,endOrbital,offset
integer :: length,startPassIndex,I,indexETUV,passIndexETUV,ip,ia,offset1,offset2

nprim = P%nprimitives
nAng = P%nAngmom
na1 = P%orbital1%nAngmom
na2 = P%orbital2%nAngmom
offset = (numPass-1)*na1
DO ia = 1,nA1
   PassP%orbital1%startOrbital(ia+offset) = P%orbital1%startOrbital(ia)
ENDDO
offset = (numPass-1)*na2
DO ia = 1,nA2
   PassP%orbital2%startOrbital(ia+offset) = P%orbital2%startOrbital(ia)
ENDDO

PassP%orb1atom(numPass)  = P%orb1atom(1)
PassP%orb2atom(numPass)  = P%orb2atom(1)
PassP%orb1batch(numPass) = P%orb1batch(1)
PassP%orb2batch(numPass) = P%orb2batch(1)
PassP%orb1mol(numPass)   = P%orb1mol(1)
PassP%orb2mol(numPass)   = P%orb2mol(1)
!   startPrim  = 1+(numPass-1)*nprim
!   endOrbital = numPass*nPrim
PassP%distance12(1+(numpass-1)*3) = P%distance12(1)
PassP%distance12(2+(numpass-1)*3) = P%distance12(2)
PassP%distance12(3+(numpass-1)*3) = P%distance12(3)
startPrim  = (numPass-1)*nprim
DO ip=1,nprim
   PassP%center(1+(startPrim+ip-1)*3) = P%center(1+(ip-1)*3)
   PassP%center(2+(startPrim+ip-1)*3) = P%center(2+(ip-1)*3)
   PassP%center(3+(startPrim+ip-1)*3) = P%center(3+(ip-1)*3)
ENDDO
DO ip=1,nprim
   PassP%iPrim1(startPrim+ip) = P%iprim1(ip)
   PassP%iPrim2(startPrim+ip) = P%iprim2(ip)
ENDDO
#ifdef VAR_LSDEBUGINT
   IF(size(P%preexpfac).NE.P%nPrimAlloc*P%nAngAlloc)call lsquit('preexpfac dim errorA',-1)
   IF(size(PassP%preexpfac).NE.PassP%nPrimAlloc*PassP%nAngAlloc)call lsquit('preexpfac dim errorB',-1)
   IF(nAng.GT.P%nAngAlloc)call lsquit('preexpfac dim errorC',-1)
   IF(nPrim.GT.P%nPrimAlloc)call lsquit('preexpfac dim errorD',-1)
   IF(nPrim+startPrim.GT.PassP%nPrimAlloc)call lsquit('preexpfac dim errorF',-1)
#endif
DO ia = 1,nAng
   offset1 = startPrim+(ia-1)*PassP%nPrimAlloc
   offset2 = (ia-1)*P%nPrimAlloc
   DO ip=1,nprim
      PassP%preExpFac(ip+offset1) = P%preExpFac(ip+offset2)
   ENDDO
ENDDO

END SUBROUTINE ADDOVERLAPTOPASS

!> \brief add the overlap P to an overlap, which we call PassP
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param PassP The overlap, to be built
!> \param P the overlap
!> \param numPass the current number of passes in PassP 
!> \param maxpasses the maximum number of passes
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE modifyOverlapCenter(PassQ,nprim,overlaplist,work,numPass,maxpasses,&
     & atom1,atom2,atomIndex1,atomIndex2,batch1,batch2,X1,Y1,Z1,X2,Y2,Z2,&
     & nRHSoverlaps,natoms1,natoms2,LUPRI,IPRINT)
Implicit none
TYPE(Overlap),intent(inout) :: PassQ
!TYPE(ODITEM),intent(in)     :: OD
integer,intent(in)          :: numPass,LUPRI,IPRINT,maxpasses,nprim
integer,intent(in)          :: overlaplist(numpass),nRHSoverlaps,natoms1,natoms2
real(realk),intent(inout)   :: work(5*nprim)
integer,intent(in)     :: atom1(nRHSoverlaps),atom2(nRHSoverlaps)
integer,intent(in)     :: batch1(nRHSoverlaps),batch2(nRHSoverlaps)
integer,intent(in)     :: atomIndex1(nRHSoverlaps),atomIndex2(nRHSoverlaps)
real(realk),intent(in) :: X1(natoms1),Y1(natoms1),Z1(natoms1)
real(realk),intent(in) :: X2(natoms2),Y2(natoms2),Z2(natoms2)
!
integer :: nang,nang1,nang2,startPrim,ipass,i1,i2,ia12
integer :: ip,irhs,work2,work3,work4,work5,nprim1,nprim2,ia1,ia2,startA2
integer :: atomC,atomD,batchC,batchD,atomIndexC,atomIndexD,offset,nprimAlloc
real(realk) :: e1,e2,XC,YC,ZC,XD,YD,ZD,XX,YY,ZZ,d2,tmp,tmp2,tmp1
real(realk),parameter :: D1 = 1E0_realk 
real(realk),pointer :: preexpfac(:),center(:)

work2 = nprim
work3 = 2*nprim
work4 = 3*nprim
work5 = 4*nprim
nAng = PassQ%nAngmom
#ifdef VAR_LSDEBUGINT
IF(size(PassQ%preexpfac).NE.PassQ%nPrimAlloc*PassQ%nAngAlloc)call lsquit('preexpfac dim errorB',-1)
#endif
preexpfac => PassQ%preexpfac
center => PassQ%center
nprim1 = PassQ%orbital1%nprimitives
nprimAlloc = PassQ%nPrimAlloc
nprim2 = PassQ%orbital2%nprimitives
nAng1 = PassQ%orbital1%nAngmom
nAng2 = PassQ%orbital2%nAngmom
DO ip=1,nprim      
   i1 = PassQ%iPrim1(ip)
   i2 = PassQ%iPrim2(ip)
   e1 = PassQ%orbital1%exponents(i1)
   e2 = PassQ%orbital2%exponents(i2)
   work(ip) = e1
   work(ip+work2) = e2
   work(ip+work3) = D1/(e1+e2)
   work(ip+work4) = -e1*e2*work(ip+work3)
ENDDO

DO IPass = 1,numPass
   startPrim  = (IPass-1)*nprim
   IRHS = overlaplist(IPass)
   atomC = atom1(IRHS) 
   atomD = atom2(IRHS) 
   batchC = batch1(IRHS) 
   batchD = batch2(IRHS) 
   atomIndexC = atomIndex1(IRHS)
   atomIndexD = atomIndex2(IRHS)
   PassQ%orb1atom(IPass)  = atomC
   PassQ%orb1batch(IPass) = batchC
   PassQ%orb2atom(IPass)  = atomD
   PassQ%orb2batch(IPass) = batchD
   PassQ%orb1mol(iPass) = atomIndexC
   PassQ%orb2mol(iPass) = atomIndexD
   XC = X1(atomC)
   YC = Y1(atomC)
   ZC = Z1(atomC)
   XD = X2(atomD)
   YD = Y2(atomD)
   ZD = Z2(atomD)
   XX = XC-XD
   YY = YC-YD
   ZZ = ZC-ZD
   PassQ%distance12(1+(IPass-1)*3) = XX
   PassQ%distance12(2+(IPass-1)*3) = YY
   PassQ%distance12(3+(IPass-1)*3) = ZZ
   d2 = XX*XX + YY*YY + ZZ*ZZ
   do ip = 1,nprim
      work(ip+work5) = work(ip+work4)*d2
   enddo
#ifdef VAR_LSDEBUGINT
IF(startPrim+nPrim.GT.PassQ%nPrimAlloc)call lsquit('preexpfac dim errorQ1',-1)
#endif
!call ls_vdexp(nprim,work(1+work5),preExpFac(startPrim+1))
!#if VAR_MKL
!   call vdexp(nprim,work(1+work5),preExpFac(startPrim+1))
!#else
   DO ip = 1,nprim
      preExpFac(startPrim+ip) = exp(work(ip+work5))
   ENDDO
!#endif
   DO ip = 1,nprim
      tmp = work(ip+work3)
      e1 = work(ip)
      e2 = work(ip+work2)
      center(1+(startPrim+ip-1)*3) = (e1*XC+e2*XD)*tmp 
      center(2+(startPrim+ip-1)*3) = (e1*YC+e2*YD)*tmp 
      center(3+(startPrim+ip-1)*3) = (e1*ZC+e2*ZD)*tmp 
   ENDDO
ENDDO
#ifdef VAR_LSDEBUGINT
IF(nAng.GT.PassQ%nAngAlloc)call lsquit('preexpfac dim errorQ1',-1)
IF(nprim*numPass.GT.nprimAlloc)call lsquit('preexpfac dim errorQ1',-1)
#endif
iF(nAng.GT. 1)THEN
   DO ia12 = 1,nAng
      offset =(ia12-1)*nprimAlloc
      DO ip=1,nprim*numPass
         preExpFac(ip+offset) = preExpFac(ip)
      ENDDO
   ENDDO
ENDIF
IF (PassQ%segmented) THEN
#ifdef VAR_LSDEBUGINT
IF(nprim*numPass.GT.nprimAlloc)call lsquit('preexpfac dim errorQ1',-1)
#endif
   iA12 = 0
   DO iA1=1,nAng1
      startA2 = 1
      IF (PassQ%sameAO) startA2 = iA1
      DO iA2=startA2,nAng2
         offset =iA12*nprimAlloc
         iA12  = iA12 + 1
#ifdef VAR_LSDEBUGINT
         IF(iA12.GT.PassQ%nAngAlloc)call lsquit('preexpfac dim errorQ2',-1)
#endif
         DO ip=1,nprim      
            i1 = PassQ%iPrim1(ip)
            i2 = PassQ%iPrim2(ip)
            tmp2 = PassQ%orbital2%CC(iA2)%p%elms(i2) 
            tmp1 = PassQ%orbital1%CC(iA1)%p%elms(i1)
            work(ip) = tmp1*tmp2
         ENDDO
         DO IPass=1,numPass
            startprim = (IPass-1)*nprim
            DO ip = 1,nprim  
               PreExpFac(ip+startPrim+offset) = &
                    PreExpFac(ip+startPrim+offset)*work(ip)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE MODIFYOVERLAPCENTER

!> \brief Jengine based integral driver (CPL 323, 425 and JCP 114, 6572)
!> \author S. Reine
!> \date 2010
!>
!> \param OD_LHS the ODbatches belonging to the Left hand side
!> \param OD_RHS the ODbatches belonging to the Reft hand side
!> \param INPUT the integral input specifications, contains all info about what to do
!> \param OUTPUT the integral output specifications, determines how the output should be given
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE Jengine(OD_LHS,OD_RHS,INPUT,OUTPUT,LUPRI,IPRINT)
implicit none
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(INTEGRALOUTPUT) :: OUTPUT
integer              :: LUPRI,IPRINT
TYPE(ODITEM)         :: OD_LHS,OD_RHS
!
Integer               :: ILHS,IRHS,Start_LHS,End_LHS,NFTUVbatches
TYPE(Integralitem)    :: Integral
TYPE(Allocitem)       :: Allocations
TYPE(TUVitem),target  :: SharedTUV
TYPE(Overlap)         :: P
TYPE(Overlap)         :: PassF
TYPE(Overlap),pointer :: Q(:)
TYPE(Overlap),pointer :: F(:)
TYPE(Integrand)       :: PQ
Integer               :: nPrimP, nPrimQ,nPrimPass,PassType,doLHS
Integer               :: maxPrimPass(20)
Logical               :: screen,dopasses
integer               :: maxPrim,maxTUVdim,maxOrb
integer(kind=long)    :: WORKLENGTH
real(realk)           :: ReductionECONT(input%NDMAT_RHS)
IF (IPRINT.GT. 5) THEN
  CALL LSHEADER(LUPRI,'Jengine')
  WRITE(LUPRI,'(1X,A)') ''
  WRITE(LUPRI,'(1X,A)') '******************************************************************'
  WRITE(LUPRI,'(1X,A)') '***                          Jengine                           ***'
  WRITE(LUPRI,'(1X,A)') '******************************************************************'
  WRITE(LUPRI,'(1X,A)') ''
ENDIF

CALL initTUVitem(sharedTUV,Input,OD_LHS,OD_RHS,LUPRI,IPRINT)
integral%TUV => sharedTUV
doPasses  = INPUT%DO_PASSES

IF(doPasses.AND..NOT.(INPUT%DO_FMM.OR.INPUT%OE_SCREEN)) THEN
  maxPrimPass   = 2
  maxPrimPass(1) = 128
  maxPrimPass(2) =  64
  maxPrimPass(3) =  32
  maxPrimPass(4) =  16
  maxPrimPass(5) =   8
  maxPrimPass(6) =   4
!  WRITE(LUPRI,*)'Jengine integrals are calculated using passes with maxPrimPass=',maxPrimPass
ELSE
  maxPrimPass = 1
ENDIF

CALL AllocItem_init(Allocations,'Both')
CALL AllocItem_zero(Allocations,'Both')
CALL SET_ALLOC(Allocations,Input,OD_LHS,'LHS',IPRINT,LUPRI)
CALL SET_ALLOC(Allocations,Input,OD_RHS,'RHS',IPRINT,LUPRI)

call MEM_OVERLAP(OD_RHS,Input,2)
call ALLOC_ODRHS_BUFFERS
call mem_alloc(Q,OD_RHS%nbatches)

!we need room for a Contraction Matrix
WORKLENGTH = Allocations%maxprimRHS*Allocations%maxContRHS
!we need room for a contraction matrix
WORKLENGTH = Allocations%maxprimRHS*Allocations%maxijkRHS*INPUT%NDMAT_RHS
!When we do a spherical transformation on Etuv we need 2 Etuvs and a spherical transform
WORKLENGTH = WORKLENGTH+2*Allocations%maxETUVlenRHS+Allocations%maxijkRHS*Allocations%maxijkRHS
!for an intermediate
WORKLENGTH = WORKLENGTH+Allocations%maxprimRHS*Allocations%maxprimLHS
call init_workmem(WORKLENGTH) 
CALL SET_FTUVbatches(F,NFTUVbatches,OD_RHS,Q,Input,SharedTUV,Integral,Allocations,&
     &               INPUT%NDMAT_RHS,maxPrimPass,LUPRI,IPRINT)
call free_workmem
call mem_dealloc(Q)
call DEALLOC_ODRHS_BUFFERS
Q=>F
NULLIFY(F)
ReductionEcont = 0.0E0_realk
IF(.NOT.INPUT%noOMP)call mem_TurnONThread_Memory()
!$OMP PARALLEL IF(.NOT.INPUT%noOMP) DEFAULT(none) PRIVATE(integral,P,PQ,ILHS,IRHS,Start_LHS,End_LHS,nPrimP,nPrimQ,&
!$OMP screen,PassF,nPrimPass,passType,doLHS,WORKLENGTH,maxPrim,maxTUVdim,&
!$OMP maxOrb) SHARED(sharedTUV,Allocations,Input,OD_LHS,OD_RHS,&
!$OMP NFTUVbatches,Q,dopasses,maxprimpass,output,IPRINT,LUPRI,ReductionEcont)
IF(.NOT.INPUT%noOMP)call init_threadmemvar()

!reassociate pointers as the may be deassociated inside a OMP PARALLEL REGION
integral%TUV => sharedTUV
IF(INPUT%LHSSameAsRHSDmat)THEN
   INPUT%LST_DLHS => INPUT%LST_DRHS
ENDIF

call INIT_BUFCOUNTERS(1)
call MEM_INIT_OVERLAP(Allocations,0,Input,1,IPRINT,LUPRI)
call ALLOC_ODLHS_BUFFERS
CALL INIT_OVERLAP(P,Allocations,0,Input,1,IPRINT,LUPRI)


!we need room for a LHS or RHS Contraction Matrix
WORKLENGTH = MAX(Allocations%maxprimRHS*Allocations%maxContRHS,Allocations%maxprimLHS*Allocations%maxContLHS) 
!When we do a spherical transformation on Etuv we need 2 Etuvs and a spherical transform
WORKLENGTH = WORKLENGTH+MAX(2*Allocations%maxETUVlenLHS+&
     & Allocations%maxijkLHS*Allocations%maxijkLHS,2*Allocations%maxETUVlenRHS+&
     & Allocations%maxijkRHS*Allocations%maxijkRHS) 
!for an intermediate
WORKLENGTH = WORKLENGTH+Allocations%maxprimRHS*Allocations%maxprimLHS
call init_workmem(WORKLENGTH) 

call getIntegralDimsJengine2(Allocations,maxPrimPass,INPUT%NDMAT_RHS,&
     & INPUT%geoderivorder,INPUT%magderivorder,maxTUVdim,maxPrim)
call init_bufcounters(5)
CALL MemFTUVbatches(maxPrim,maxTUVdim,INPUT%NDMAT_RHS,5)
call ALLOC_ODPASSF_BUFFERS

maxTUVdim = getIntegralDimsJengine(Allocations,maxPrimPass,INPUT%NDMAT_RHS,INPUT%geoderivorder,INPUT%magderivorder)
maxPrim   = Allocations%maxPrimLHS*Allocations%maxPrimRHS
maxOrb    = Allocations%maxTotOrbRHS*Allocations%maxTotOrbLHS*INPUT%NDMAT_RHS
CALL allocIntegrals(PQ,Integral,maxPrim,maxTUVdim,maxorb,lupri)
call mem_alloc(Integral%Econt,input%NDMAT_RHS)
IF(INPUT%fullcontraction)Integral%Econt = 0.0E0_realk   
!$OMP DO SCHEDULE(DYNAMIC,1)
DO ILHS = 1,OD_LHS%nbatches
   CALL SET_Overlap(P,Input,SharedTUV,Integral,OD_LHS%BATCH(ILHS),&
     &  1,LUPRI,IPRINT,.FALSE.)
#ifdef VAR_LSDEBUGINT
    IF(P%nPrimitives*P%nTUV*input%ndmat_rhs.GT.allocIntmaxTUVdim)THEN
       call lsquit('ls_dzero in Jengine alloc error',lupri)
    ENDIF
#endif
  CALL LS_DZERO(Integral%tuvQ,P%nPrimitives*P%nTUV*input%ndmat_rhs)
  nPrimPass = 0
  passType  = 0
  doLHS = 0 ! this is an integer that keeps track of the RHS loop if all the RHS loop is screened away
            ! then doLHS remains 0 and the contract_P should not be called.
  DO IRHS = 1,NFTUVbatches
      screen = getScreening(P,Q(IRHS),INPUT,LUPRI,IPRINT)
      IF (.NOT.screen) THEN
        IF (IPRINT.GT. 3) WRITE(LUPRI,'(1X,A,I3,A,I3)') 'Overlap distributions P',ILHS,' and Q',IRHS
!       ******************************* Passes *****************************
!       Integrals are calculated using passes (collecting similar FTUV-batches 
!       before performing integration to reduce computation overhead)
        IF (doPasses) THEN
!         If new pass-type calculate contribution from previous pass type
          IF ((passType.NE.Q(IRHS)%passType).AND.(nPrimPass.GT. 0)) THEN
              CALL JengineInnerCont(PQ,P,PassF,INTEGRAL,INPUT,ILHS,IRHS, &
     &                              INPUT%NDMAT_RHS,LUPRI,IPRINT)
              nPrimPass = 0
              CALL free_Overlap(PassF)
              CALL REINIT_OD_PASSF
          ENDIF
!         If the OD-batch has more primivites that the maximum do not add to pass
!         but calculate directly
          IF (Q(IRHS)%nPrimitives.GE.maxPrimPass(Q(IRHS)%endAngmom+1)) THEN
              CALL JengineInnerCont(PQ,P,Q(IRHS),INTEGRAL,INPUT,ILHS,IRHS, &
     &                              INPUT%NDMAT_RHS,LUPRI,IPRINT)
          ELSE
            IF (nPrimPass.EQ. 0) THEN
              passType = Q(IRHS)%passType
              CALL InitFTUVbatches(PassF,passType,maxPrimPass(Q(IRHS)%endAngmom+1),&
     &                             Q(IRHS)%startAngmom, Q(IRHS)%endAngmom,Q(IRHS)%nTUV,&
     &                             INPUT%NDMAT_RHS,5)
            ENDIF

!           Calcluate old pass first if number of primitives would exceed
!           the maximal number
            IF ((Q(IRHS)%nPrimitives+nPrimPass).GT.maxPrimPass(Q(IRHS)%endAngmom+1)) THEN
              CALL JengineInnerCont(PQ,P,PassF,INTEGRAL,INPUT,ILHS,IRHS, &
     &                              INPUT%NDMAT_RHS,LUPRI,IPRINT)
              nPrimPass = 0
              PassF%nPrimitives = 0
            ENDIF
            CALL AddODtoFTUV(PassF,Q(IRHS)%nPrimitives,nPrimPass,Q(IRHS),INPUT%NDMAT_RHS,LUPRI,IPRINT)
            nPrimPass = nPrimPass + Q(IRHS)%nPrimitives
          ENDIF
!       ****************************** NO Passes ***************************
!       Integrals are calculated batchwise
        ELSE
            CALL JengineInnerCont(PQ,P,Q(IRHS),INTEGRAL,INPUT,ILHS,IRHS, &
     &                            INPUT%NDMAT_RHS,LUPRI,IPRINT)
        ENDIF
        doLHS = 1
     ENDIF
  ENDDO
  IF(doLHS .EQ. 1)THEN
     IF (doPasses.AND.(nPrimPass.GT. 0)) THEN
        CALL JengineInnerCont(PQ,P,PassF,INTEGRAL,INPUT,ILHS,IRHS, &
             &                     INPUT%NDMAT_RHS,LUPRI,IPRINT)
     ENDIF
     Integral%startDerivativeIndex = 0
     Integral%nDerivComp = input%NGEODERIVCOMP
     CALL setIntegralDerivativOrders(Integral,input%GEODERIVORDER,input%GEODERIVORDER,P%single,PassF%single)
     CALL Contract_P(INTEGRAL,PQ,input,LUPRI,IPRINT)
     CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
     IF (doPasses.AND.(nPrimPass.GT. 0)) THEN
        nPrimPass = 0
        CALL free_Overlap(PassF)
        CALL REINIT_OD_PASSF
     ENDIF
  ENDIF
ENDDO
!$OMP END DO NOWAIT

call DEALLOC_ODPASSF_BUFFERS
IF(INPUT%fullcontraction)THEN
!$OMP CRITICAL
 do iLHS=1,input%NDMAT_RHS
    ReductionEcont(iLHS) = ReductionEcont(iLHS) + Integral%Econt(iLHS)
 enddo   
!$OMP END CRITICAL
ENDIF
call mem_dealloc(Integral%Econt)
CALL deallocIntegrals(PQ,Integral)
CALL FREE_OVERLAP(P)
call DEALLOC_ODLHS_BUFFERS
call free_workmem

IF(.NOT.INPUT%noOMP)call collect_thread_memory()
!$OMP END PARALLEL
IF(.NOT.INPUT%noOMP)call mem_TurnOffThread_Memory()

IF(INPUT%fullcontraction)THEN
   do iLHS=1,input%NDMAT_RHS
      Output%ResultTensor%LSAO(1)%elms(iLHS) = ReductionEcont(iLHS)
   enddo
ENDIF
DO IRHS = 1,NFTUVbatches
  IF (Q(IRHS)%nPrimitives.GT. 0) CALL FREE_OVERLAP(Q(IRHS)) !maybe not needed
ENDDO
call mem_dealloc(Q)
call DEALLOC_ODFTUV_BUFFERS

CALL freeTUVitem(sharedTUV,Input)

END SUBROUTINE Jengine

!> \brief The inner contraction required in a Jengine based integration
!> \author S. Reine
!> \date 2010
!>
!>  \param PQ contains integrand info like reduced exponents, integral prefactor, etc.
!>  \param P contains overlap distribution for the left hand side, electron 1
!>  \param Q contains overlap distribution for the right hand side, electron 2
!>  \param Integral contains arrays to store intermidiates and final integrals
!>  \param INPUT the integral input specifications, contains all info about what to do
!>  \param ILHS the index for which LHS overlap distribution this is, relevant when ILHS = IRHS  
!>  \param IRHS the index for which RHS overlap distribution this is, relevant when ILHS = IRHS  
!>  \param NDMAT the number of density matrices
!>  \param LUPRI the logical unit number for the output file
!>  \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE JengineInnerCont(PQ,P,Q,INTEGRAL,INPUT,ILHS,IRHS,NDMAT,LUPRI,IPRINT)
implicit none
TYPE(INTEGRALINPUT) :: INPUT
integer             :: LUPRI,IPRINT,ILHS,IRHS,NDMAT
TYPE(Integralitem)  :: Integral
TYPE(Overlap)       :: P
TYPE(Overlap)       :: Q
TYPE(Integrand)     :: PQ
!
CALL Build_Integrand(PQ,P,Q,INPUT,ILHS,IRHS,LUPRI,IPRINT)
!   Hermite 2-electron integral over general operator w 
CALL Build_HermiteTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT)
CALL ContractFTUV(INTEGRAL,P,Q,NDMAT,LUPRI,IPRINT)
END SUBROUTINE JengineInnerCont

!> \brief LinK based integral driver (JCP 109,1663) 
!> \author T. Kjaergaard
!> \date 2010
!>
!>  \param OD_LHS the ODbatches belonging to the Left hand side
!>  \param OD_RHS the ODbatches belonging to the Reft hand side
!>  \param INPUT the integral input specifications, contains all info about what to do
!>  \param LSOUTPUT the integral output specifications, determines how the output should be given
!>  \param LUPRI the logical unit number for the output file
!>  \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE LINK_DRIVER(OD_LHS,OD_RHS,INPUT,LSOUTPUT,LUPRI,IPRINT)
  !now just change AB -> BA (with everything else consistent)
implicit none
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(INTEGRALOUTPUT) :: LSOUTPUT
integer              :: LUPRI,IPRINT
TYPE(ODITEM)         :: OD_LHS,OD_RHS
!
Integer,target        :: ILHS,STATICSIZEOFDOINT
Integer               :: IRHS,Start_RHS,End_RHS,i,j
TYPE(TUVitem),target  :: SharedTUV
TYPE(Integralitem)    :: Integral
TYPE(Allocitem)       :: Allocations
TYPE(Overlap)         :: P
TYPE(Overlap),pointer :: Q(:)
TYPE(Overlap),pointer :: PassQ(:)
TYPE(Integrand)       :: PQ
Integer,pointer       :: ODpassesIndex(:)
Integer               :: nPassTypes,nRHSoverlaps,nLHSoverlaps
Integer,pointer       :: Maxpassfortype(:),numpasses(:)
integer(kind=short),pointer :: RED_GAB_LHS(:)
integer(kind=short),pointer :: RED_GAB_RHS(:)
integer(kind=short),pointer :: RED_GAB_TMP(:)
integer(kind=short),pointer :: RED_DMAT_LHS(:,:)
integer(kind=short),pointer :: RED_DMAT_RHS(:,:)
integer(kind=short),pointer :: maxLHSGAB(:),maxRHSGAB(:)
integer(kind=short)   :: maxLHSELM,maxRHSELM,DMATELM1,DMATELM2,MAXDMAT
integer(kind=short),pointer :: SORTING(:,:),SORTING2(:)
integer(kind=short) :: CS_THRLOG,REDGABLHS,REDDMATRHSBD,REDDMATRHSAD,TMP_short
integer(kind=short) :: DALINK_THRLOG
Integer               :: dim1,dim2,dim3,dim4,A,B,C,D,idmat,endC
Integer               :: LISTSIZE!,LISTSIZE1,LISTSIZE2
Integer               :: nB,nC,nD,OLDI,RHSINDEX,startB,startD,I2,JK
Integer,pointer   :: LIST(:)!,LIST2(:),LIST1(:)
LOGICAL,pointer   :: DoINT(:)
Integer(kind=long)  :: WORKLENGTH,WORKLENGTH2,WORKLENGTH3,WORKEST1
real(realk),pointer :: TMPWORK(:)
!real(realk),allocatable :: SORTING(:)
!integer,allocatable     :: bshell(:)
!LOGICAL               :: A_LOOP_DONE,B_LOOP_DONE,C_LOOP_DONE,D_LOOP_DONE,UNIQUE
TYPE(LINKshell),pointer  :: brashell(:),ketshell(:),ML(:)
real(realk)           :: CPUTimeLINK,WALLTIMELINK,WALLTIME3,WALLTIME4
real(realk)           :: CPUTIMESTART,WALLTIMESTART,CPUTIMEEND,WALLTIMEEND
logical               :: set_orbital1
integer               :: maxpasses,iPassType,startc,nLHSODTypes,iODtype,iodLHS
integer,pointer   :: batchindex1(:),batchindex2(:)
integer,pointer   :: SIZEOFDOINT
logical        :: NOELEMENTSADDED,LHS,screen,dalink,doneBuild,DRHS_SYM,DLHS_SYM,sameODs,INPUTDO_PASSES
real(realk)    :: mbieP(2),RAB,CPUTIME3,CPUTIME4
type(lstensor),pointer :: pGAB
type(lstensor) :: localKmat
logical :: MBIE_SCREEN,nograd,nonSR_EXCHANGE,sameLHSaos,PerformCALC
logical,pointer :: dopasses(:)
Integer,pointer  :: ODTypeIndex(:)
type(integerpointer),pointer :: PassTypeOverlapindex(:)
integer,pointer :: nOverlapOfPassType(:),AIndex(:),BIndex(:)
integer :: TOTmaxpasses,natoms3,natoms4
type(integerpointer),pointer :: TypeOverlapIndex(:)
Integer,pointer       :: maxPassesFortypes(:)
integer :: atomA,atomB,batchA,batchB,Gindex,IP,currentODtype,ILHSCOUNT,iODtype2
integer,pointer :: OverlapList(:,:),ILHSCOUNTINDEX(:)

integer,pointer     :: atomC(:),atomD(:),batchC(:),batchD(:),atomIndexC(:),atomIndexD(:)
real(realk),pointer :: X3(:),Y3(:),Z3(:),X4(:),Y4(:),Z4(:)
real(realk) :: TMP,factor,CS_THRESHOLD
integer :: IRHSI(1),nLHSbatches,nA,node,numnodes
integer :: nthreads,tid,nDMAT_RHS,nDMAT_LHS,IOMPLHSCOUNT
integer,pointer :: Belms(:)
integer,pointer :: IODelms(:)
real(realk)           :: ReductionECONT(input%NDMAT_RHS)
#ifdef VAR_OMP
integer, external :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#endif
CALL LS_GETTIM(CPUTIMESTART,WALLTIMESTART)
!IF(.NOT. INPUT%CS_SCREEN)THEN
!   WRITE(LUPRI,*)'Link requires screening, but you have turned off screening'
!   CALL LSQUIT('Link requires screening, but you have turned off screening',lupri)
!ENDIF

dalink = INPUT%DO_DALINK
DRHS_SYM = INPUT%DRHS_SYM
DLHS_SYM = INPUT%DLHS_SYM
nDMAT_RHS = INPUT%nDMAT_RHS
nDMAT_LHS = INPUT%nDMAT_LHS
CS_THRESHOLD = INPUT%CS_THRESHOLD/INPUT%exchangeFactor
sameLHSaos = INPUT%sameLHSaos
call Obtain_CS_THRLOG(CS_THRLOG,INPUT%CS_THRESHOLD/INPUT%exchangeFactor)
DALINK_THRLOG = CS_THRLOG-INPUT%DASCREEN_THRLOG
sameODs = INPUT%sameODs
MBIE_SCREEN = INPUT%MBIE_SCREEN
nonSR_EXCHANGE = .NOT.(Input%Operator.EQ.ErfcOperator)
nLHSoverlaps = OD_LHS%nbatches
nRHSoverlaps = OD_RHS%nbatches
INPUTDO_PASSES = INPUT%DO_PASSES
IF (IPRINT.GT. 5) THEN
  IF(DALINK)THEN
     CALL LSHEADER(LUPRI,'DaLinK')
     WRITE(LUPRI,'(1X,A)') ''
     WRITE(LUPRI,'(1X,A)') '******************************************************************'
     WRITE(LUPRI,'(1X,A)') '***                         DaLINK                             ***'
     WRITE(LUPRI,'(1X,A)') '******************************************************************'
     WRITE(LUPRI,'(1X,A)') ''
  ELSE
     CALL LSHEADER(LUPRI,'LinK')
     WRITE(LUPRI,'(1X,A)') ''
     WRITE(LUPRI,'(1X,A)') '******************************************************************'
     WRITE(LUPRI,'(1X,A)') '***                           LINK                             ***'
     WRITE(LUPRI,'(1X,A)') '******************************************************************'
     WRITE(LUPRI,'(1X,A)') ''
  ENDIF
ENDIF
CALL initTUVitem(sharedTUV,Input,OD_LHS,OD_RHS,LUPRI,IPRINT)
Integral%TUV => sharedTUV
CALL AllocItem_init(Allocations,'Both')
CALL Allocitem_zero(Allocations,'Both')
call mem_alloc(ODTypeIndex,nLHSoverlaps)
call SelectODTypesFromODbatch(ODtypeIndex,Allocations,OD_LHS,nLHSoverlaps,&
     &nLHSODTypes,TYPEOVERLAPINDEX,INPUT,IPRINT,LUPRI,'LHS')

call mem_alloc(ODpassesIndex,nRHSoverlaps)
call SelectODTypesFromODbatch(ODpassesIndex,Allocations,OD_RHS,nRHSoverlaps,&
     &nPassTypes,PassTYPEOVERLAPINDEX,INPUT,IPRINT,LUPRI,'RHS')
IF(nPassTypes.EQ.nRHSoverlaps)INPUTDO_PASSES = .FALSE.

call mem_alloc(nOverlapOfPassType,nPassTypes)
DO Ipasstype = 1,nPassTypes
   nOverlapOfPassType(Ipasstype)=0
ENDDO
DO IRHS=1,nRHSoverlaps
   nOverlapOfPassType(ODpassesIndex(IRHS))=nOverlapOfPassType(ODpassesIndex(IRHS))+1
ENDDO

!================================================================================
! This section is a Link Specific part, building reduced gab matrices
! determining significant brashell pairs and so on 
!================================================================================
dim1=INPUT%AO(1)%p%nbatches
dim2=INPUT%AO(2)%p%nbatches
dim3=INPUT%AO(3)%p%nbatches
dim4=INPUT%AO(4)%p%nbatches

IF(INPUT%CS_SCREEN)THEN
   call mem_alloc(RED_GAB_RHS,nRHSoverlaps)
   pGAB => input%LST_GAB_RHS
   DO IRHS=1,nRHSoverlaps
      RED_GAB_RHS(IRHS) = OD_RHS%BATCH(IRHS)%maxgab
   ENDDO
   call mem_alloc(maxRHSGAB,dim4)
   maxRHSGAB = shortzero
   DO IRHS=1,nRHSoverlaps
      D = OD_RHS%BATCH(IRHS)%IB
      maxRHSGAB(D) = MAX(maxRHSGAB(D),RED_GAB_RHS(IRHS))
   ENDDO
   maxRHSELM = shortzero
   DO D=1,dim4
      maxRHSELM = MAX(maxRHSELM,maxRHSGAB(D))
   ENDDO

   IF(sameODs)THEN
      RED_GAB_LHS => RED_GAB_RHS
      maxLHSELM = maxRHSELM
      maxLHSGAB => maxRHSGAB 
   ELSE
      pGAB => input%LST_GAB_LHS
      call mem_alloc(RED_GAB_LHS,nLHSoverlaps)
      DO ILHS=1,nLHSoverlaps
         RED_GAB_LHS(ILHS) = OD_LHS%BATCH(ILHS)%maxgab
      ENDDO
      call mem_alloc(maxLHSGAB,dim2)
      maxLHSGAB = shortzero
      DO ILHS=1,nLHSoverlaps
         B = OD_LHS%BATCH(ILHS)%IB
         maxLHSGAB(B) = MAX(maxLHSGAB(B),RED_GAB_LHS(ILHS))
      ENDDO
      maxLHSELM = shortzero
      DO B=1,dim2
         maxLHSELM = MAX(maxLHSELM,maxLHSGAB(B))
      ENDDO
   ENDIF
ELSE
   !no screening so we manual set all screening quantities to 1000
   call mem_alloc(RED_GAB_RHS,nRHSoverlaps)
   DO IRHS=1,nRHSoverlaps
      RED_GAB_RHS(IRHS) = 3
   ENDDO
   call mem_alloc(maxRHSGAB,dim4)
   DO D=1,dim4
      maxRHSGAB(D) = 3
   ENDDO
   maxRHSELM = 3

   IF(sameODs)THEN
      RED_GAB_LHS => RED_GAB_RHS
      maxLHSELM = maxRHSELM
      maxLHSGAB => maxRHSGAB 
   ELSE
      call mem_alloc(RED_GAB_LHS,nLHSoverlaps)
      DO ILHS=1,nLHSoverlaps
         RED_GAB_LHS(ILHS) = 3
      ENDDO
      call mem_alloc(maxLHSGAB,dim2)
      DO ILHS=1,dim2
         maxLHSGAB(ILHS) = 3
      ENDDO
      maxLHSELM = 3
   ENDIF
ENDIF

IF(INPUT%CS_SCREEN)THEN
   PerformCALC = maxLHSELM .GT. CS_THRLOG-maxRHSELM
ELSE
   PerformCALC=.TRUE.
ENDIF

IF(PerformCALC)THEN
 ! determine significant ket shell pairs
 NULLIFY(ketshell)
 ALLOCATE(ketshell(dim4))
 call mem_alloc(SORTING,dim3,dim4)
 CALL DETERMINE_SHELL_PAIRS_RHS(dim4,dim3,RED_GAB_RHS,CS_THRLOG,&
         &maxLHSELM,ketshell,OD_RHS,INPUT%sameRHSaos,SORTING,&
         &nRHSoverlaps,lupri)
 DO D=1,dim4
    IF(ketshell(D)%DIM.GT.0)THEN
       call mem_alloc(ketshell(D)%RED_GAB,ketshell(D)%DIM)
       DO nC=1,ketshell(D)%DIM
          C=ketshell(D)%belms(nC)
          IRHS=ketshell(D)%IODelms(nC)
          ketshell(D)%RED_GAB(nC) = RED_GAB_RHS(IRHS)
       ENDDO
    ENDIF
 ENDDO 

 ! determine significant bra shell pairs
 NULLIFY(brashell)
 ALLOCATE(brashell(dim2))
 IF(.NOT.sameODs)THEN
    call mem_dealloc(SORTING)
    call mem_alloc(SORTING,dim1,dim2)
 ENDIF
 CALL DETERMINE_SHELL_PAIRS_LHS(dim1,dim2,RED_GAB_LHS,CS_THRLOG,&
      &maxRHSELM,brashell,OD_LHS,sameLHSaos,SORTING,&
      &nLHSoverlaps,lupri)
 call mem_dealloc(SORTING)
 ! build reduced density matrix
 IF(INPUT%CS_SCREEN)THEN
    call mem_alloc(RED_DMAT_RHS,dim2,dim4)
    call Build_full_shortint2dim_from_lstensor(input%LST_DRHS,RED_DMAT_RHS,dim2,dim4,lupri)
    IF (input%sameODs) call symmetrize_SDMAT(RED_DMAT_RHS,dim2,dim4)
    !WRITE(lupri,*)'The Dmat output'
    !call shortint_output(RED_DMAT_RHS,dim1,dim3,lupri)
 ELSE
    !no screening so we set the screening Dmat quantity to 1000
    call mem_alloc(RED_DMAT_RHS,dim2,dim4)
    DO D=1,dim4
       DO C=1,dim2
          RED_DMAT_RHS(C,D) = 3
       ENDDO
    ENDDO
 ENDIF

 IF(DALINK)THEN
    IF(SameODs)THEN
       RED_DMAT_LHS => RED_DMAT_RHS 
    ELSE
       call mem_alloc(RED_DMAT_LHS,dim1,dim3)
       IF(INPUT%CS_SCREEN)THEN
          call Build_full_shortint2dim_from_lstensor(input%LST_DLHS,RED_DMAT_LHS,dim1,dim3,lupri)
       ELSE
          DO D=1,dim3
             DO C=1,dim1
                RED_DMAT_LHS(C,D) = 3
             ENDDO
          ENDDO
       ENDIF
    ENDIF
    call mem_alloc(batchindex1,nRHSoverlaps)
    call mem_alloc(batchindex2,nRHSoverlaps)
    DO D=1,dim4
       IF(ketshell(D)%DIM.GT.0)THEN
          DO nC=1,ketshell(D)%DIM
             C=ketshell(D)%belms(nC)
             IRHS=ketshell(D)%IODelms(nC)
             batchindex1(IRHS)=C
             batchindex2(IRHS)=D
          ENDDO
       ENDIF
    ENDDO
 ENDIF
 NULLIFY(ML)
 ALLOCATE(ML(dim2))
 call mem_alloc(SORTING2,dim4)
 CALL DETERMINE_BRAKET_PAIRS(dim1,dim2,dim4,RED_DMAT_RHS,&
      &maxLHSGAB,maxRHSGAB,CS_THRLOG,ML,sameODs,SORTING2)
 call mem_dealloc(SORTING2)
 !================================================================================
 ! Done Link Specific part, now we set up the ODtypes and so on
 !================================================================================
 call mem_alloc(AIndex,nLHSoverlaps)
 call mem_alloc(BIndex,nLHSoverlaps)
 call mem_alloc(ILHSCOUNTINDEX,nLHSoverlaps)
 ILHSCOUNT=0
 DO iODType=1,nLHSODtypes
    DO B = 1,dim2
       IODelms => brashell(B)%IODelms
       Belms => brashell(B)%Belms
       DO nA=1,brashell(B)%DIM
          A=Belms(nA)
          IF(sameLHSaos)THEN
             IF(B.GE.A)THEN
                ILHS=IODelms(nA)
                IODType2 = ODtypeIndex(ILHS)
                IF(IODType2.EQ.iODtype)THEN
                   ILHSCOUNT=ILHSCOUNT+1
                   ILHSCOUNTINDEX(ILHSCOUNT) = ILHS
                   Aindex(ILHSCOUNT) = A
                   Bindex(ILHSCOUNT) = B
                ENDIF
             ENDIF
          ELSE
             ILHS=IODelms(nA)
             IODType2 = ODtypeIndex(ILHS)
             IF(IODType2.EQ.iODtype)THEN
                ILHSCOUNT=ILHSCOUNT+1
                ILHSCOUNTINDEX(ILHSCOUNT) = ILHS
                Aindex(ILHSCOUNT) = A
                Bindex(ILHSCOUNT) = B
             ENDIF
          ENDIF
       ENDDO
    ENDDO
 ENDDO
 nLHSbatches = ILHSCOUNT
 CALL LS_GETTIM(CPUTIME3,WALLTIME3)
 ReductionEcont = 0.0E0_realk

 IF(.NOT.INPUT%noOMP) call mem_TurnONThread_Memory()
!$OMP PARALLEL IF(.NOT.INPUT%noOMP) DEFAULT(NONE) PRIVATE(ILHSCOUNT,integral,PQ,IRHS,&
!$OMP RHSINDEX,A,nB,B,ILHS,nC,C,nD,D,LIST,LISTSIZE,TMPWORK,WORKLENGTH,DoINT,&
!$OMP NOELEMENTSADDED,DMATELM1,DMATELM2,MAXDMAT,iPasstype,numpasses,PassQ,Q,&
!$OMP mbieP,rab,SIZEOFDOINT,STATICSIZEOFDOINT,MAXpassesForTypes,dopasses,P,&
!$OMP totmaxpasses,screen,overlaplist,REDGABLHS,REDDMATRHSAD,REDDMATRHSBD,&
!$OMP TMP_short,atomC,atomD,atomIndexC,atomIndexD,batchC,batchD,X3,Y3,Z3,X4,Y4,Z4,IRHSI,localKmat,&
!$OMP iODType,currentODtype,tid,nthreads,WORKLENGTH2,WORKLENGTH3,WORKEST1,IODelms,&
!$OMP Belms,RED_GAB_TMP,node,numnodes,natoms3,natoms4,IOMPLHSCOUNT) SHARED(Input,&
!$OMP nPassTypes,nLHSODtypes,CS_THRESHOLD,nonSR_EXCHANGE,&
!$OMP Allocations,OD_LHS,OD_RHS,iprint,lupri,sharedTUV,nOverlapOfPassType,&
!$OMP PassTypeOverlapIndex,TypeOverlapIndex,ODpassesIndex,dim3,dim4,sameods,&
!$OMP aindex,bindex,mbie_screen,ml,ketshell,brashell,red_dmat_rhs,red_dmat_lhs,&
!$OMP cs_thrlog,lsoutput,drhs_sym,dalink,dalink_thrlog,&
!$OMP batchindex1,batchindex2,red_gab_lhs,red_gab_rhs,inputdo_passes,&
!$OMP nRHSoverlaps,ndmat_rhs,factor,nograd,ILHSCOUNTINDEX,&
!$OMP ODtypeIndex,nLHSbatches,ReductionECONT)
 IF(.NOT.INPUT%noOMP) call init_threadmemvar()
#ifdef VAR_MPI
 node = input%node
 numnodes = input%numnodes
#else
 node = 0
 numnodes = 1
#endif
#ifdef VAR_OMP
 nthreads=OMP_GET_NUM_THREADS()
 tid=OMP_GET_THREAD_NUM()
!$OMP MASTER
 if(tid==0)then
   IF(node.EQ.0.AND.IPRINT.GT.1)WRITE(lupri,'(4X,A,I3,A)')'This is an OpenMP calculation using ',&
     & omp_get_num_threads(),' threads.' 
 endif
!$OMP END MASTER

!$OMP BARRIER
#else
 nthreads=1
 tid=0
#endif

!IF(nograd)then
!   call copy_lstensor_to_lstensor(lsoutput%resultTensor,localKmat)
!ENDIF

!reassociate pointers as the may be deassociated inside a OMP PARALLEL REGION
 integral%TUV => sharedTUV
 IF(INPUT%LHSSameAsRHSDmat)THEN
    INPUT%LST_DLHS => INPUT%LST_DRHS
 ENDIF

 call buildRHS_centerinfo(nAtoms3,nAtoms4,nRHSoverlaps,atomC,atomD,atomIndexC,atomIndexD,&
     & batchC,batchD,X3,Y3,Z3,X4,Y4,Z4,INPUT,OD_RHS,IPRINT,LUPRI)

 call mem_alloc(Q,nPassTypes)
 !RHS OVERLAP MEMORY CALCULATION
 call INIT_BUFCOUNTERS(2)
 DO iPassType=1,nPassTypes
    IRHS = PassTypeOverlapIndex(iPassType)%elms(1) 
    call MEM_SINGLE_OVERLAP(OD_RHS,IRHS,Input,2)
 ENDDO   
 call ALLOC_ODRHS_BUFFERS
 !LHS OVERLAP MEMORY CALCULATION
 call INIT_BUFCOUNTERS(1)
 DO iODType=1,nLHSODTypes
    call MEM_INIT_OVERLAP(Allocations,iODtype,Input,1,IPRINT,LUPRI)
 ENDDO
 !ALLOC LHS AND RHS OVERLAP BUFFERS
 call ALLOC_ODLHS_BUFFERS

 DO iPassType=1,nPassTypes
    IRHS = PassTypeOverlapIndex(iPassType)%elms(1) 
    CALL SET_Overlap(Q(iPassType),Input,SharedTUV,Integral,OD_RHS%BATCH(IRHS),2,LUPRI,IPRINT,.TRUE.)
 ENDDO

 call mem_alloc(numpasses,nPassTypes)
 numpasses = 0
 call mem_alloc(LIST,dim3*dim4+nPassTypes)
 call mem_alloc(DoINT,dim3*dim4)
 NULLIFY(SIZEOFDOINT)
 IF(sameODs)THEN
    SIZEOFDOINT => ILHS
 ELSE
    STATICSIZEOFDOINT = dim3*dim4
    SIZEOFDOINT => STATICSIZEOFDOINT
 ENDIF
 call mem_alloc(maxPassesFortypes,nPassTypes)
 call mem_alloc(doPasses,nPassTypes)
 !initial allocation
 ILHSCOUNT = MIN(1+tid,nLHSbatches)
 ILHS = ILHSCOUNTINDEX(ILHSCOUNT)
 iODType = ODtypeIndex(ILHS)
 currentODtype = iODType
 CALL INIT_OVERLAP(P,Allocations,iODtype,Input,1,IPRINT,LUPRI)
 CALL allocIntegralsWRAP(PQ,Integral,Input,Allocations,iODtype,nPassTypes,&
      &maxPassesFortypes,1,INPUTDO_PASSES,nOverlapOfPassType,lupri)
 call mem_alloc(Integral%Econt,input%NDMAT_RHS)
 IF(INPUT%fullcontraction)Integral%Econt = 0.0E0_realk   
 TOTmaxpasses = 0
 call INIT_BUFCOUNTERS(3)
 DO iPassType=1,nPassTypes
    TOTmaxpasses = MAX(TOTmaxpasses,maxPassesForTypes(iPassType))
    doPasses(iPassType) = maxPassesForTypes(iPassType).NE. 1
    IF(doPasses(iPassType))THEN
       CALL MEM_PASS_FROM_OVERLAP(Q(Ipasstype),Input,maxPassesForTypes(iPassType))
    ENDIF
 ENDDO
 call ALLOC_ODPASS_BUFFERS
 call mem_alloc(PassQ,nPassTypes)
 TOTmaxpasses = 0
 DO iPassType=1,nPassTypes
    TOTmaxpasses = MAX(TOTmaxpasses,maxPassesForTypes(iPassType))
    doPasses(iPassType) = maxPassesForTypes(iPassType).NE. 1
    IF(doPasses(iPassType))THEN
       CALL INIT_PASS_FROM_OVERLAP(PassQ(iPassType),Q(Ipasstype),Input,maxPassesForTypes(iPassType),IPRINT,LUPRI)
       CALL SetPassOrbitals(PassQ(iPassType),Q(Ipasstype),maxPassesForTypes(iPassType),LUPRI,IPRINT)
    ENDIF
 ENDDO
 !alloc workarray
 !LHS requirements
 WORKLENGTH = 0
 DO iODType=1,nLHSODtypes
    WORKEST1 = Allocations%maxprimLHSA(iODtype)*Allocations%maxContLHSA(iODtype)&
         &+2*Allocations%maxETUVlenLHSA(iODtype)+Allocations%maxijkLHSA(iODtype)*Allocations%maxijkLHSA(iODtype)
    WORKLENGTH = MAX(WORKLENGTH,WORKEST1)
 ENDDO
 !intermediate
 WORKLENGTH = WORKLENGTH+Allocations%maxprimRHS*Allocations%maxprimLHS
 
 !RHS requirements
 WORKLENGTH2 = 0
 DO iODType=1,nLHSODtypes
    CALL determineMaxPassesForType(Allocations,iODtype,nPassTypes,&
         &maxPassesFortypes,1,INPUTDO_PASSES,nOverlapOfPassType,lupri)
    DO iPassType=1,nPassTypes
       WORKEST1 = Allocations%maxprimRHSA(iPassType)*Allocations%maxContRHSA(iPassType)&
            &+2*Allocations%maxETUVlenRHSA(iPassType)*maxPassesForTypes(iPassType)&
            &+Allocations%maxijkRHSA(iPassType)*Allocations%maxijkRHSA(iPasstype)
       WORKLENGTH2 = MAX(WORKLENGTH2,WORKEST1)
       WORKEST1 = Allocations%maxprimRHSA(iPassType)*Allocations%maxprimLHSA(iODtype)&
            &*maxPassesForTypes(iPassType)
       WORKLENGTH2 = MAX(WORKLENGTH2,WORKEST1)
    ENDDO
 ENDDO
 WORKLENGTH3 = MAX(WORKLENGTH2,WORKLENGTH)
 !intermediate
 WORKEST1 = 5*Allocations%maxPrimRHS
 WORKLENGTH3 = MAX(WORKLENGTH3,WORKEST1)
 call init_workmem(WORKLENGTH3)
 
 CALL determineMaxPassesForType(Allocations,currentODtype,nPassTypes,&
      & maxPassesFortypes,1,INPUTDO_PASSES,nOverlapOfPassType,lupri)
 iODtype = currentODtype 
 
 call mem_alloc(overlaplist,TOTmaxpasses,npassTypes)
 !!$OMP DO SCHEDULE(DYNAMIC,1)
 IOMPLHSCOUNT = 0
 DO ILHSCOUNT=1+node,nLHSbatches,numnodes
  !each node does different of theses ILHSCOUNT
  IOMPLHSCOUNT = IOMPLHSCOUNT+1
  !Each thread does different of these IOMPLHSCOUNT
  !Without MPI IOMPLHSCOUNT goes (1,2,..,nLHSbatches)
  !and Thread1 takes 1, Thread2 takes 2, ...
  !With 4 MPI procs IOMPLHSCOUNT goes (1,5,9,13,...)
  !and Thread1 takes 1, Thread2 takes 5, ...
  IF(MOD(IOMPLHSCOUNT,nthreads).EQ.tid)THEN
   ILHS = ILHSCOUNTINDEX(ILHSCOUNT)
   iODType = ODtypeIndex(ILHS)
   IF(iODType.NE.currentODtype)THEN
      !dealloc
      CALL FREE_OVERLAP(P)
      CALL REINIT_OD_LHS
      CALL deallocIntegrals(PQ,Integral)
      DO iPassType=1,nPassTypes
         IF(doPasses(iPassType))THEN
            CALL Free_overlap(PassQ(iPassType))
         ENDIF
      ENDDO
      call mem_dealloc(PassQ)
      call mem_dealloc(overlaplist)
      !realloc
      CALL INIT_OVERLAP(P,Allocations,iODtype,Input,1,IPRINT,LUPRI)
      CALL allocIntegralsWRAP(PQ,Integral,Input,Allocations,iODtype,nPassTypes,&
           &maxPassesFortypes,1,INPUTDO_PASSES,nOverlapOfPassType,lupri)

      call DEALLOC_ODPASS_BUFFERS
      TOTmaxpasses = 0
      DO iPassType=1,nPassTypes
         TOTmaxpasses = MAX(TOTmaxpasses,maxPassesForTypes(iPassType))
         doPasses(iPassType) = maxPassesForTypes(iPassType).NE. 1
         IF(doPasses(iPassType))THEN
            CALL MEM_PASS_FROM_OVERLAP(Q(iPassType),Input,maxPassesForTypes(iPassType))
         ENDIF
      ENDDO
      call ALLOC_ODPASS_BUFFERS

      call mem_alloc(PassQ,nPassTypes)
      TOTmaxpasses = 0
      DO iPassType=1,nPassTypes
         TOTmaxpasses = MAX(TOTmaxpasses,maxPassesForTypes(iPassType))
         IRHS = PassTypeOverlapIndex(iPassType)%elms(1) 
         doPasses(iPassType) = maxPassesForTypes(iPassType).NE. 1
         IF(doPasses(iPassType))THEN
            CALL INIT_PASS_FROM_OVERLAP(PassQ(iPassType),Q(Ipasstype),Input,maxPassesForTypes(iPassType),IPRINT,LUPRI)
            CALL SetPassOrbitals(PassQ(iPassType),Q(Ipasstype),maxPassesForTypes(iPassType),LUPRI,IPRINT)
         ENDIF
      ENDDO
      call mem_alloc(overlaplist,TOTmaxpasses,npassTypes)
      currentODtype = iODtype
   ENDIF
   CALL SET_Overlap(P,Input,SharedTUV,Integral,OD_LHS%BATCH(ILHS),1,LUPRI,IPRINT,.FALSE.)
   A = AIndex(ILHSCOUNT)
   B = BIndex(ILHSCOUNT)
   REDGABLHS = RED_GAB_LHS(ILHS)
   !CREATING DoINT LIST from the BD element of the density matrix
   IF(MBIE_SCREEN)THEN
      mbieP(1)=P%MBIE(1)
      mbieP(2)=P%MBIE(2)
      RAB = P%ODextent
   ENDIF
   DO IRHS=1,SIZEOFDOINT
      DoINT(IRHS) = .FALSE.
   ENDDO
   LOOPC: DO nD=1,ML(B)%DIM
      D = ML(B)%Belms(nD) 
      REDDMATRHSBD = RED_DMAT_RHS(B,D)
      TMP_short = REDDMATRHSBD+REDGABLHS
      NOELEMENTSADDED = .TRUE.

      IF(ketshell(D)%DIM.GT.0)THEN
         IODelms => ketshell(D)%IODelms
         Belms => ketshell(D)%Belms
         RED_GAB_TMP => ketshell(D)%RED_GAB
         LOOPD: DO nC=1,ketshell(D)%DIM
            C=Belms(nC)
            IRHS = IODelms(nC)
            !SINCE THINGS ARE NOW IN BATCHES WITH THE SAME VALUES THIS COULD MAYBE BE MODIFIED
            IF(TMP_short+RED_GAB_TMP(nC) .LE. CS_THRLOG )EXIT LOOPD
            NOELEMENTSADDED = .FALSE.
            IF(MBIE_SCREEN)THEN
               screen = .FALSE.
               CALL MBIE_SCREENING(screen,RAB,P,mbieP,OD_RHS%BATCH(IRHS)%ODextent,&
                    & OD_RHS%BATCH(IRHS)%ODCenter,input%lst_GAB_RHS%MBIE(1,C,D),&
                    & input%lst_GAB_RHS%MBIE(2,C,D),CS_THRESHOLD,&
                    & RED_DMAT_RHS(B,D),nonSR_EXCHANGE,Input%AttOmega,lupri)
               DoINT(IRHS) = .NOT.screen
            ELSE
               DoINT(IRHS) = .TRUE.
            ENDIF
         ENDDO LOOPD
      ENDIF
      IF(NOELEMENTSADDED)EXIT LOOPC !NO ELEMENTS ADDED IN D LOOP
   ENDDO LOOPC
   IF(INPUT%sameLHSaos)then
      !Add to DoINT from the AD element of the density matrix
      LOOPC2: DO nD=1,ML(A)%DIM
         D = ML(A)%Belms(nD)
         REDDMATRHSAD = RED_DMAT_RHS(A,D)
         TMP_short = REDDMATRHSAD+REDGABLHS
         NOELEMENTSADDED = .TRUE.
         IF(ketshell(D)%DIM.GT.0)THEN
            IODelms => ketshell(D)%IODelms
            Belms => ketshell(D)%Belms
            RED_GAB_TMP => ketshell(D)%RED_GAB
            LOOPD2: DO nC=1,ketshell(D)%DIM
               C=Belms(nC)
               IRHS = IODelms(nC)
               IF(TMP_short+RED_GAB_TMP(nC) .LE. CS_THRLOG ) EXIT LOOPD2
               NOELEMENTSADDED = .FALSE.
               IF(MBIE_SCREEN)THEN
                  screen = .FALSE.
                  CALL MBIE_SCREENING(screen,RAB,P,mbieP,OD_RHS%BATCH(IRHS)%ODextent,&
                       & OD_RHS%BATCH(IRHS)%ODCenter,input%lst_GAB_RHS%MBIE(1,C,D),&
                       & input%lst_GAB_RHS%MBIE(2,C,D),CS_THRESHOLD,&
                       & RED_DMAT_RHS(A,D),nonSR_EXCHANGE,Input%AttOmega,lupri)
                  DoINT(IRHS) = .NOT.screen
               ELSE
                  DoINT(IRHS) = .TRUE.
               ENDIF
            ENDDO LOOPD2
         ENDIF
         IF(NOELEMENTSADDED)EXIT LOOPC2!NO ELEMENTS ADDED IN D LOOP
      ENDDO LOOPC2
   ENDIF
   !Do diagonal =================================================
   IF(sameODs)THEN
      ipasstype = ODpassesIndex(ILHS)
      IRHSI(1) = ILHS
      call mem_workpointer_alloc(TMPWORK,5*Q(iPassType)%nPrimitives)
      CALL modifyOverlapCenter(Q(iPassType),Q(iPassType)%nPrimitives,&
           & IRHSI,TMPWORK,&
           & 1,1,atomC,atomD,atomIndexC,atomIndexD,batchC,batchD,X3,Y3,Z3,X4,Y4,Z4,&
           & nRHSoverlaps,natoms3,natoms4,LUPRI,IPRINT)
      call mem_workpointer_dealloc(TMPWORK)
      CALL ExplicitIntegrals(Integral,PQ,P,Q(ipasstype),INPUT,LSOUTPUT,&
           & ILHS,ILHS,LUPRI,IPRINT)
      CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,LSOUTPUT,LUPRI,IPRINT)
      DoINT(ILHS) = .FALSE.
   ENDIF
   !Done diagonal ================================================
   CALL BUILD_LIST(LUPRI,LIST,DoINT,LISTSIZE,dim3*dim4,SIZEOFDOINT)
   LOOPRHS: DO RHSINDEX=1,LISTSIZE
      IRHS=LIST(RHSINDEX)
      iPassType = ODpassesIndex(IRHS)
      IF(DALINK)THEN
         !EXTRA DENSITY ACCELERATED SCREENING
         C=batchindex1(IRHS)
         D=batchindex2(IRHS)
         DMATELM1=RED_DMAT_RHS(B,D)+RED_DMAT_LHS(A,C)
         IF(INPUT%sameLHSaos)then
            DMATELM2=RED_DMAT_RHS(A,D)+RED_DMAT_LHS(B,C)
            MAXDMAT=MAX(DMATELM1,DMATELM2)
         ELSE
            MAXDMAT=DMATELM1
         ENDIF
         IF(MAXDMAT.LT.shortzero)MAXDMAT=shortzero
         IF(MAXDMAT+REDGABLHS+RED_GAB_RHS(IRHS) .LT. DALINK_THRLOG ) CYCLE LOOPRHS
      ENDIF
      IF(doPasses(iPassType)) THEN
         numPasses(iPassType) = numPasses(iPassType) + 1
         overlaplist(numPasses(iPassType),ipassType) = IRHS
         IF (numPasses(iPassType).EQ.maxPassesForTypes(iPassType)) THEN

            call mem_workpointer_alloc(TMPWORK,5*Q(iPassType)%nPrimitives)
            CALL modifyOverlapCenter(PassQ(iPassType),PassQ(iPassType)%nPrimitives,&
                 & overlaplist(1:numPasses(iPassType),iPasstype),&
                 & TMPWORK,&
                 & numPasses(iPassType),maxPassesForTypes(iPassType),&
                 & atomC,atomD,atomIndexC,atomIndexD,batchC,batchD,X3,Y3,Z3,X4,Y4,Z4,&
                 & nRHSoverlaps,natoms3,natoms4,LUPRI,IPRINT)
            call mem_workpointer_dealloc(TMPWORK)
            IRHS = overlaplist(1,iPassType)
            CALL FinalizePass(PassQ(iPassType),Q(Ipasstype),maxPassesForTypes(iPassType),LUPRI,IPRINT)
            CALL ExplicitIntegrals(Integral,PQ,P,PassQ(iPassType),INPUT,LSOUTPUT,&
                 & ILHS,IRHS,LUPRI,IPRINT)
            CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,LSOUTPUT,LUPRI,IPRINT)
            numPasses(iPassType) = 0
         ENDIF
      ELSE
         IRHSI(1) = IRHS
         call mem_workpointer_alloc(TMPWORK,5*Q(iPassType)%nPrimitives)
         CALL modifyOverlapCenter(Q(iPassType),Q(iPassType)%nPrimitives,&
              & IRHSI,TMPWORK,&
              & 1,1,atomC,atomD,atomIndexC,atomIndexD,batchC,batchD,X3,Y3,Z3,X4,Y4,Z4,&
              & nRHSoverlaps,natoms3,natoms4,LUPRI,IPRINT)
         call mem_workpointer_dealloc(TMPWORK)
         CALL ExplicitIntegrals(Integral,PQ,P,Q(Ipasstype),INPUT,LSOUTPUT,&
              & ILHS,IRHS,LUPRI,IPRINT)
         CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,LSOUTPUT,LUPRI,IPRINT)
      ENDIF
   ENDDO LOOPRHS
   DO iPassType=1,nPassTypes
      IF(numPasses(iPassType).GT. 0) THEN
         call mem_workpointer_alloc(TMPWORK,5*Q(iPassType)%nPrimitives)
         CALL modifyOverlapCenter(PassQ(iPassType),PassQ(iPassType)%nPrimitives,&
              & overlaplist(1:numPasses(iPassType),iPasstype),&
              & TMPWORK,&
              & numPasses(iPassType),maxPassesForTypes(iPassType),&
              & atomC,atomD,atomIndexC,atomIndexD,batchC,batchD,X3,Y3,Z3,X4,Y4,Z4,&
              & nRHSoverlaps,natoms3,natoms4,LUPRI,IPRINT)
         call mem_workpointer_dealloc(TMPWORK)
         IRHS = overlaplist(1,iPassType)
         CALL FinalizePass(PassQ(iPassType),Q(Ipasstype),numPasses(iPassType),LUPRI,IPRINT)
         IRHS=ILHS+1
         !IRHS has been used to build Q(IRHS) which was then added to 
         !PassQ and for all Q(IRHS), IRHS was different from ILHS.
         !which was treated as a special case.  
         !IRHS do no longer have any meaning, but is not allowed be accidentally be equal ILHS
         !because then a triangular loop will be used which PassQ was not built for
         CALL ExplicitIntegrals(Integral,PQ,P,PassQ(iPassType),INPUT,LSOUTPUT,&
              & ILHS,IRHS,LUPRI,IPRINT)
         CALL DistributeIntegrals(INTEGRAL,PQ,INPUT,LSOUTPUT,LUPRI,IPRINT)
         numPasses(iPassType) = 0
      ENDIF
   ENDDO
   !numpass is now 0 for all ipasstype, and ready for next loop
  ENDIF
 ENDDO
 !!$OMP END DO NOWAIT

 CALL FREE_OVERLAP(P)
 call DEALLOC_ODLHS_BUFFERS
 call DEALLOC_ODRHS_BUFFERS

 IF(INPUT%fullcontraction)THEN
 !$OMP CRITICAL
  do i=1,input%NDMAT_RHS
   ReductionECONT(i) = ReductionECONT(i) + Integral%Econt(i)
  enddo
 !$OMP END CRITICAL
 ENDIF
 call mem_dealloc(Integral%Econt)
 CALL deallocIntegrals(PQ,Integral)
 DO iPassType=1,nPassTypes
    IF(doPasses(iPassType))THEN
       CALL Free_overlap(PassQ(iPassType))
    ENDIF
 ENDDO
 call mem_dealloc(PassQ)
 call DEALLOC_ODPASS_BUFFERS
 call free_workmem
 call mem_dealloc(overlaplist)
 call mem_dealloc(LIST)
 call mem_dealloc(DoINT)
 call mem_dealloc(numpasses)
 call mem_dealloc(maxPassesFortypes)
 call mem_dealloc(doPasses)
 DO iPassType=1,nPassTypes
    CALL FREE_Overlap(Q(iPassType))
 ENDDO
 call mem_dealloc(Q)
 
 call freeRHS_centerinfo(atomC,atomD,atomIndexC,atomIndexD,batchC,batchD,&
      & X3,Y3,Z3,X4,Y4,Z4,INPUT%sameRHSaos)
 !IF(nograd)then
 !!$OMP CRITICAL
 !   call add_lstensor_to_lstensor(localKmat,lsoutput%resultTensor)
 !!$OMP END CRITICAL
 !   call lstensor_free(localKmat)
 !ENDIF

 IF(.NOT.INPUT%noOMP)call collect_thread_memory()
!$OMP END PARALLEL
 IF(.NOT.INPUT%noOMP)call mem_TurnOffThread_Memory()

 IF(INPUT%fullcontraction)THEN
  do ILHS=1,input%NDMAT_RHS
   LsOutput%ResultTensor%LSAO(1)%elms(ILHS) = ReductionECONT(ILHS)
  enddo
 ENDIF

 IF (IPRINT.GT. 5) THEN
   IF (LSTIME_PRINT)THEN
      CALL LS_GETTIM(CPUTIME4,WALLTIME4)
      CPUTIMELINK = CPUTIME4-CPUTIME3
      WALLTIMELINK = WALLTIME4-WALLTIME3
      WRITE(lupri,*)'OVERALL WALL TIMINGS in LinkLoop'
      CALL ls_TIMTXT('>>>  WALL Time used in LinkLoop is         ',WALLTIMELINK,LUPRI)
      WRITE(lupri,*)'OVERALL CPU TIMINGS in LinkLoop'
      CALL ls_TIMTXT('>>>  CPU Time used in LinkLoop is         ',CPUTIMELINK,LUPRI)
   ENDIF
 ENDIF

 call mem_dealloc(ILHSCOUNTINDEX)
 call mem_dealloc(BIndex)
 call mem_dealloc(AIndex)

 DO B=1,dim2
    call linkshell_free(ML(B))
 ENDDO
 DEALLOCATE(ML)
 NULLIFY(ML)

 IF(DALINK)THEN
    IF(.NOT.SameODs)call mem_dealloc(RED_DMAT_LHS)
    call mem_dealloc(batchindex1)
    call mem_dealloc(batchindex2)
 ENDIF
 call mem_dealloc(RED_DMAT_RHS)

 DO B=1,dim2
    IF(brashell(B)%DIM.NE. 0)THEN
       call linkshell_free(brashell(B))
    ENDIF
 ENDDO
 DEALLOCATE(brashell)
 NULLIFY(brashell)

 DO D=1,dim4
    IF(ketshell(D)%DIM.NE. 0)THEN
       call mem_dealloc(ketshell(D)%RED_GAB)
       call linkshell_free(ketshell(D))
    ENDIF
 ENDDO
 DEALLOCATE(ketshell)
 NULLIFY(ketshell)
ENDIF !performCalc

IF(.NOT.sameODs)THEN
   call mem_dealloc(maxLHSGAB)
   call mem_dealloc(RED_GAB_LHS)
ENDIF
call mem_dealloc(maxRHSGAB)
call mem_dealloc(RED_GAB_RHS)
call mem_dealloc(nOverlapofPassType)
call mem_dealloc(ODpassesIndex)
call mem_dealloc(ODTypeIndex)
DO iPassType=1,nPasstypes
   call mem_dealloc(PassTypeOverlapIndex(iPassType)%elms)
ENDDO
DEALLOCATE(PassTypeOverlapIndex)
DO iODType=1,nLHSODtypes
   call mem_dealloc(TypeOverlapIndex(iODType)%elms)
ENDDO
DEALLOCATE(TypeOverlapIndex)
call allocitem_free(Allocations,'Both')
CALL freeTUVitem(sharedTUV,Input)

CALL LS_GETTIM(CPUTIMEEND,WALLTIMEEND)
CPUTIMELINK = CPUTIMEEND-CPUTIMESTART
WALLTIMELINK = WALLTIMEEND-WALLTIMESTART

IF (IPRINT.GT. 5) THEN
   IF (LSTIME_PRINT)THEN
      WRITE(lupri,*)'OVERALL WALL TIMINGS '
      CALL ls_TIMTXT('>>>  WALL Time used in Link is             ',WALLTIMELINK,LUPRI)
      WRITE(lupri,*)'OVERALL CPU TIMINGS '
      CALL ls_TIMTXT('>>>  CPU Time used in Link is             ',CPUTIMELINK,LUPRI)
   ENDIF
ENDIF
END SUBROUTINE LINK_DRIVER

subroutine symmetrize_SDMAT(MAT2,dim1,dim2)
  implicit none
  integer             :: dim1,dim2
  integer(kind=short) :: MAT2(dim1,dim2)
  !
  integer(kind=short) :: TMP(dim1,dim2)
  integer             :: I,J,TMPIJ,TMPJI,MAXTMP
  DO J=1,dim2
     DO I=1,J
        TMPIJ = MAT2(I,J)
        TMPJI = MAT2(J,I)
        MAXTMP = MAX(TMPIJ,TMPJI)
        MAT2(I,J) = MAXTMP
        MAT2(J,I) = MAXTMP
     ENDDO
  ENDDO
end subroutine symmetrize_SDMAT

subroutine buildRHS_centerinfo(nAtoms3,nAtoms4,nRHSoverlaps,atomC,atomD,&
    &atomIndexC,atomIndexD,batchC,batchD,X3,Y3,Z3,X4,Y4,Z4,input,OD_RHS,IPRINT,LUPRI)
implicit none
TYPE(INTEGRALINPUT),intent(in)  :: INPUT
TYPE(ODITEM),intent(in)         :: OD_RHS
integer,intent(inout) :: natoms3,natoms4
integer :: iprint,lupri,nRHSoverlaps,IRHS
logical :: sameRHSaos
integer,pointer :: atomC(:),atomD(:),batchC(:),batchD(:),atomIndexC(:),atomIndexD(:)
real(realk),pointer :: X3(:),Y3(:),Z3(:),X4(:),Y4(:),Z4(:)

nAtoms3 = INPUT%AO(3)%p%natoms
nAtoms4 = INPUT%AO(4)%p%natoms
call mem_alloc(atomC,nRHSoverlaps)
call mem_alloc(atomD,nRHSoverlaps)
call mem_alloc(atomIndexC,nRHSoverlaps)
call mem_alloc(atomIndexD,nRHSoverlaps)
call mem_alloc(batchC,nRHSoverlaps)
call mem_alloc(batchD,nRHSoverlaps)
call mem_alloc(X3,nAtoms3)
call mem_alloc(Y3,nAtoms3)
call mem_alloc(Z3,nAtoms3)
DO IRHS=1,nRHSoverlaps
   atomC(IRHS)     = OD_RHS%BATCH(IRHS)%AO(1)%p%atom
   atomIndexC(IRHS)= OD_RHS%BATCH(IRHS)%AO(1)%p%molecularIndex
   batchC(IRHS)    = OD_RHS%BATCH(IRHS)%AO(1)%p%batch
   atomD(IRHS)     = OD_RHS%BATCH(IRHS)%AO(2)%p%atom
   atomIndexD(IRHS)= OD_RHS%BATCH(IRHS)%AO(2)%p%molecularIndex
   batchD(IRHS)    = OD_RHS%BATCH(IRHS)%AO(2)%p%batch
   X3(atomC(IRHS)) = OD_RHS%BATCH(IRHS)%AO(1)%p%center(1)
   Y3(atomC(IRHS)) = OD_RHS%BATCH(IRHS)%AO(1)%p%center(2)
   Z3(atomC(IRHS)) = OD_RHS%BATCH(IRHS)%AO(1)%p%center(3)
ENDDO
IF(INPUT%sameRHSaos)THEN
   X4 => X3
   Y4 => Y3
   Z4 => Z3
ELSE
   call mem_alloc(X4,nAtoms4)
   call mem_alloc(Y4,nAtoms4)
   call mem_alloc(Z4,nAtoms4)
   DO IRHS=1,nRHSoverlaps
      X4(atomD(IRHS)) = OD_RHS%BATCH(IRHS)%AO(2)%p%center(1)
      Y4(atomD(IRHS)) = OD_RHS%BATCH(IRHS)%AO(2)%p%center(2)
      Z4(atomD(IRHS)) = OD_RHS%BATCH(IRHS)%AO(2)%p%center(3)
   ENDDO
ENDIF

end subroutine buildRHS_centerinfo

subroutine freeRHS_centerinfo(atomC,atomD,atomIndexC,atomIndexD,batchC,batchD,&
     & X3,Y3,Z3,X4,Y4,Z4,sameRHSaos)
implicit none
logical :: sameRHSaos
integer,pointer :: atomC(:),atomD(:),batchC(:),batchD(:),atomIndexC(:),atomIndexD(:)
real(realk),pointer :: X3(:),Y3(:),Z3(:),X4(:),Y4(:),Z4(:)

call mem_dealloc(atomC)
call mem_dealloc(atomIndexC)
call mem_dealloc(batchC)
call mem_dealloc(atomD)
call mem_dealloc(atomIndexD)
call mem_dealloc(batchD)
call mem_dealloc(X3)
call mem_dealloc(Y3)
call mem_dealloc(Z3)
IF(.NOT.sameRHSaos)THEN
   call mem_dealloc(X4)
   call mem_dealloc(Y4)
   call mem_dealloc(Z4)
ENDIF
end subroutine freeRHS_centerinfo

!> \brief build list of RHS ODbatches to calculate, based on logical array
!> \author T. Kjaergaard
!> \date 2010
!>
!> \param LUPRI the logical unit number for the output file
!> \param LIST the list to be built
!> \param DOINT the logical array from which the LIST should be built
!> \param LISTSIZE the size of the list to be built
!> \param SIZE the allocd size of LIST and Doint
!> \param ILHS the ODbatch index for the LHS OD batch
SUBROUTINE BUILD_LIST(LUPRI,LIST,DoINT,LISTSIZE,SIZE,ILHS)
IMPLICIT NONE
Integer  :: LISTSIZE,SIZE,LUPRI,ILHS
Integer  :: LIST(SIZE)
Logical  :: DoINT(SIZE)
!
Integer  :: K,J

K=0
DO J=1,ILHS
   IF(DoINT(J))THEN
      K=K+1
      LIST(K)=J
   ENDIF
ENDDO
LISTSIZE=K

END SUBROUTINE BUILD_LIST

!> \brief determines the maximum gab matrix element, and array of max gab elements for given bacthindex
!> \author T. Kjaergaard
!> \date 2010
!>
!> \param RED_GAB the gab matrix reduced to odbatches
!> \param dim1 the size of the first dimension
!> \param dim2 the size of the second dimension
!> \param maxGab the array of maximum gab matrix elements
!> \param maxElm the maximum gab matrix element
SUBROUTINE MAXGABELM(RED_GAB,dim1,dim2,maxGAB,maxELM)
IMPLICIT NONE
INTEGER      :: dim1,dim2 
REAL(REALK)  :: maxELM,RED_GAB(dim1,dim2),maxGAB(dim1)
!
INTEGER      :: A,B

maxELM=0
DO A=1,dim1
   maxGAB(A)= 0
   DO B=1,dim2
      maxGAB(A)=MAX(maxGAB(A),RED_GAB(A,B)) 
   ENDDO
   maxELM=MAX(MAXELM,maxGAB(A))
ENDDO

END SUBROUTINE MAXGABELM

!> \brief determines the significant shell pairs (used in LinK) see JCP 109,1663
!> \author T. Kjaergaard
!> \date 2010
!>
!> \param dim1 the size of the first dimension
!> \param dim2 the size of the second dimension
!> \param RED_GAB_LHS the gab matrix reduced to odbatches
!> \param CS_THRLOG the Cauchy-Schwarz screening threshold
!> \param maxRHSElm the maximum rhs gab matrix element
!> \param brashell the significant shell pairs 
!> \param OD the ODbatch
!> \param LHS flag to determine if this is a left hand side (LHS) shell or RHS
!> \param sameRHS flag to determine if the RHS atomic orbitals are the same
!> \param nPrimOD number of primitives for each ODbatch
SUBROUTINE DETERMINE_SHELL_PAIRS_LHS(dim1,dim2,RED_GAB_LHS,CS_THRLOG,&
     &maxRHSELM,brashell,OD,sameLHS,SORTING,ODnbatches,lupri)
IMPLICIT NONE
integer                 :: dim1,dim2,ODnbatches,lupri
integer(kind=short)     :: RED_GAB_LHS(ODnbatches),CS_THRLOG,maxRHSELM
TYPE(LINKshell)         :: brashell(:)
TYPE(ODITEM)            :: OD
logical                 :: sameLHS
integer(kind=short)     :: SORTING(dim1,dim2)
!
!real(realk)             :: SORTING(dim2,dim1)
!integer                 :: bshell(dim2,dim1),IODindex(dim1,dim2),number_bshell(dim1)
integer,pointer         :: bshell(:,:),IODindex(:,:),number_bshell(:)
integer                 :: A,startB,I,B,C,IOD,nB,IRHS
integer(kind=short) :: CRIT

call mem_alloc(bshell,dim1,dim2)
call mem_alloc(IODindex,dim2,dim1)
call mem_alloc(number_bshell,dim2)
CRIT = CS_THRLOG-maxRHSELM
IF(CRIT.LT.shortzero)CRIT=shortzero
number_bshell = 0
DO IOD=1,ODnbatches
   A = OD%BATCH(IOD)%IA
   B = OD%BATCH(IOD)%IB
   IF(A.EQ.B)THEN
      IF(RED_GAB_LHS(IOD) .GT. CRIT)THEN
         number_bshell(B) = number_bshell(B)+1
         bshell(number_bshell(B),B) = A
         SORTING(number_bshell(B),B) = RED_GAB_LHS(IOD)
         IODindex(B,A) = IOD
      ENDIF
   ELSE
      IF(RED_GAB_LHS(IOD) .GT. CRIT)THEN
         number_bshell(B) = number_bshell(B)+1
         bshell(number_bshell(B),B) = A
         SORTING(number_bshell(B),B) = RED_GAB_LHS(IOD)
         IODindex(B,A) = IOD
         IF(sameLHS)then
            number_bshell(A) = number_bshell(A)+1
            bshell(number_bshell(A),A) = B
            SORTING(number_bshell(A),A) = RED_GAB_LHS(IOD)
            IODindex(A,B) = IOD   
         ENDIF
      ENDIF
   ENDIF
ENDDO
DO B=1,dim2
   IF(number_bshell(B).GT. 0)THEN
      I = number_bshell(B)
      call linkshell_init(brashell(B),I)         
      CALL linksort(SORTING(1:I,B),bshell(1:I,B),I,CS_THRLOG)
      DO C = 1,I
         brashell(B)%Belms(C)=bSHELL(C,B)
         brashell(B)%IODelms(C)=IODindex(B,bSHELL(C,B))
      ENDDO
   ELSE
      brashell(B)%DIM = 0
   ENDIF
ENDDO
call mem_dealloc(bshell)
call mem_dealloc(IODindex)
call mem_dealloc(number_bshell)

END SUBROUTINE DETERMINE_SHELL_PAIRS_LHS

!> \brief determines the significant shell pairs (used in LinK) see JCP 109,1663
!> \author T. Kjaergaard
!> \date 2010
!>
!> \param dim1 the size of the first dimension
!> \param dim2 the size of the second dimension
!> \param RED_GAB_LHS the gab matrix reduced to odbatches
!> \param CS_THRLOG the Cauchy-Schwarz screening threshold
!> \param maxRHSElm the maximum rhs gab matrix element
!> \param brashell the significant shell pairs 
!> \param OD the ODbatch
!> \param LHS flag to determine if this is a left hand side (LHS) shell or RHS
!> \param sameRHS flag to determine if the RHS atomic orbitals are the same
!> \param nPrimOD number of primitives for each ODbatch
SUBROUTINE DETERMINE_SHELL_PAIRS_RHS(dim4,dim3,RED_GAB_RHS,CS_THRLOG,&
     &maxLHSELM,ketshell,OD,sameRHS,SORTING,ODnbatches,lupri)
IMPLICIT NONE
integer                 :: dim3,dim4,ODnbatches,lupri
integer(kind=short)     :: RED_GAB_RHS(ODnbatches),CS_THRLOG,maxLHSELM
TYPE(LINKshell)         :: ketshell(:)
TYPE(ODITEM)            :: OD
logical                 :: sameRHS
integer(kind=short)     :: SORTING(dim3,dim4)
!
integer,pointer         :: bshell(:,:),IODindex(:,:),number_bshell(:)
integer                 :: C,D,I,A,IOD
integer(kind=short) :: CRIT

call mem_alloc(bshell,dim3,dim4)
call mem_alloc(IODindex,dim4,dim3)
call mem_alloc(number_bshell,dim4)
CRIT = CS_THRLOG-maxLHSELM
IF(CRIT.LT.shortzero)CRIT=shortzero

number_bshell = 0
DO IOD=1,OD%nbatches
   C = OD%BATCH(IOD)%IA
   D = OD%BATCH(IOD)%IB
   IF(C.EQ.D)THEN
      IF(RED_GAB_RHS(IOD) .GT. CRIT)THEN
         number_bshell(D) = number_bshell(D)+1
         bshell(number_bshell(D),D) = C
         SORTING(number_bshell(D),D) = RED_GAB_RHS(IOD)
         IODindex(D,C) = IOD
      ENDIF
   ELSE
      IF(RED_GAB_RHS(IOD) .GT. CRIT)THEN
         number_bshell(D) = number_bshell(D)+1
         bshell(number_bshell(D),D) = C
         SORTING(number_bshell(D),D) = RED_GAB_RHS(IOD)
         IODindex(D,C) = IOD
         IF(sameRHS)then
            number_bshell(C) = number_bshell(C)+1
            bshell(number_bshell(C),C) = D
            SORTING(number_bshell(C),C) = RED_GAB_RHS(IOD)
            IODindex(C,D) = IOD
         endif
      ENDIF
   ENDIF
ENDDO
DO D=1,dim4
   IF(number_bshell(D).GT. 0)THEN
      I = number_bshell(D)
      call linkshell_init(ketshell(D),I)
      CALL linksort(SORTING(1:I,D),bshell(1:I,D),I,CS_THRLOG)
      DO A = 1,I
         ketshell(D)%Belms(A)=bSHELL(A,D)
         ketshell(D)%IODelms(A)=IODindex(D,bSHELL(A,D))
      ENDDO
   ELSE
      ketshell(D)%DIM = 0
   ENDIF
ENDDO

call mem_dealloc(bshell)
call mem_dealloc(IODindex)
call mem_dealloc(number_bshell)

END SUBROUTINE DETERMINE_SHELL_PAIRS_RHS

!> \brief determines the significant braket pairs (used in LinK)see JCP 109,1663
!> \author T. Kjaergaard
!> \date 2010
!>
!> \param dim1 the number of Atomic orbital shells on center 1
!> \param dim2 the number of Atomic orbital shells on center 2
!> \param dim3 the number of Atomic orbital shells on center 3
!> \param RED_DMAT_RHS the density matrix in shell dimensions
!> \param maxLHSGAB the maximum lhs gab matrix elements for given shell index
!> \param maxRHSGAB the maximum rhs gab matrix elements for given shell index
!> \param CS_THRESHOLD the Cauchy-Schwarz screening threshold
!> \param ML the significant braket pairs 
!> \param sameODs flag to determine if the LHS and RHS atomic orbitals are the same
SUBROUTINE DETERMINE_BRAKET_PAIRS(dim1,dim2,dim4,RED_DMAT_RHS,maxLHSGAB,maxRHSGAB,CS_THRLOG,ML,sameODs,SORTING)
IMPLICIT NONE
integer                 :: dim1,dim2,dim4,NDMAT
integer(kind=short)     :: RED_DMAT_RHS(dim2,dim4),CS_THRLOG,maxRHSELM
integer(kind=short)     :: maxLHSGAB(dim2),maxRHSGAB(dim4)
TYPE(LINKshell)         :: ML(dim2)
logical                 :: sameODs
integer(kind=short)     :: SORTING(dim4)
!
integer                 :: bshell(dim4)
integer                 :: A,B,D,I,endD
DO B=1,dim2
   I=0
   endD=dim4
   DO D=1,endD
      IF(RED_DMAT_RHS(B,D)+maxLHSGAB(B)+maxRHSGAB(D).GT. CS_THRLOG)THEN
         I=I+1
         bshell(I)=D
         SORTING(I)=RED_DMAT_RHS(B,D)+maxRHSGAB(D)
      ENDIF
   ENDDO
   IF(I .EQ. 0)THEN
      call linkshell_init(ML(B),I)
   ELSE
      ML(B)%DIM=I
      CALL Linksort(SORTING(1:I),bshell(1:I),I,CS_THRLOG)
      call linkshell_init(ML(B),I)
      DO D=1,I
         ML(B)%Belms(D)=bSHELL(D)
      ENDDO
   ENDIF
ENDDO
END SUBROUTINE DETERMINE_BRAKET_PAIRS

!*************************************************************************
!*
!*                 End of Link
!*
!*************************************************************************

!> \brief wrapper routine which determines which integral distribution routine to call 
!> \author T. Kjaergaard
!> \date 2010
!>
!> \param Integral contains the integral which should be distributed
!> \param PQ contains info about the integrand 
!> \param Input contains info about the requested integral evaluation
!> \param Output contains output lstensor 
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE DistributeIntegrals(Integral,PQ,Input,Output,LUPRI,IPRINT)
implicit none
Type(integralitem)   :: Integral
Type(integrand)      :: PQ
Type(IntegralInput)  :: Input
Type(IntegralOutput) :: Output
Integer              :: LUPRI,IPRINT,dimQ,dimP,nMAT,nSPHMAT,mmorder
Integer              :: nA,nB,nDeriv,nPasses,dimA,dimB
Integer              :: lhsDer,rhsDer

IF(Input%do_prop)THEN
   dimA = PQ%P%p%orbital1%totOrbitals
   dimB = PQ%P%p%orbital2%totOrbitals
   IF(INPUT%fullcontraction)THEN
      call distributePropExpVal(&
           & output%ResultTensor%LSAO(1)%elms,&
           & integral%integralsABCD,input%LST_DLHS,dimA,dimB,&
           & integral%nOperatorComp,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
   ELSE
      CALL DistributePropIntegrals(OUTPUT%resultTensor,Integral%integralsABCD,&
           &dimA,dimB,integral%nOperatorComp,PQ,INPUT,OUTPUT,LUPRI,IPRINT)
   ENDIF
ELSE
   lhsDer = integral%lhsGeoOrder
   rhsDer = integral%rhsGeoOrder
   dimP=PQ%P%p%totOrbitals(lhsDer+1)
   dimQ=PQ%Q%p%totOrbitals(rhsDer+1)

   IF (Input%CS_int) THEN
      IF(OUTPUT%RealGabMatrix)THEN
         !Cauchy-Schwarz screening integrals
         call distributeCS(OUTPUT%resultTensor,PQ,Integral%integralsABCD,&
              & dimQ,dimP,Input,Output,LUPRI,IPRINT)      
      ELSE
         !Cauchy-Schwarz screening integrals
         call distributeCS(OUTPUT%ScreenTensor,PQ,Integral%integralsABCD,&
              & dimQ,dimP,Input,Output,LUPRI,IPRINT)      
      ENDIF
   ELSEIF(Input%PS_int) THEN
      !Cauchy-Schwarz screening integrals
      call distributePS(OUTPUT%ScreenTensor,PQ,Integral%integralsABCD,&
           & dimQ,dimP,Input,Output,LUPRI,IPRINT)
   ELSE
      IF (Input%DO_FOCK) THEN
       IF (Input%DO_Coulomb) THEN
        IF (Input%DO_JENGINE) THEN
         CALL distributeJengine(OUTPUT%resultTensor,PQ,Integral%IntegralsABCD,&
              & dimQ,dimP,Input,output,Integral%Econt,LUPRI,IPRINT)
        ELSE
         IF(INPUT%fullcontraction)THEN
          call distributeJcont(&
               & integral%integralsABCD,input%LST_DLHS,input%LST_DRHS,&
               & input%NDMAT_RHS,PQ,Input,output,Integral%Econt,LUPRI,IPRINT)
         ELSE
          CALL distributeCoulomb(OUTPUT%resultTensor,integral%integralsABCD,&
               & input%LST_DRHS,input%NDMAT_RHS,PQ,Input,output,LUPRI,IPRINT)
         ENDIF
        ENDIF
       ENDIF
       IF(Input%DO_Exchange)THEN
        IF(Input%DO_GRADIENT)THEN
         IF (input%NDMAT_LHS.NE.input%NDMAT_RHS) THEN
          CALL LSQUIT('ndmat mismatch in distributeIntegrals',-1)
         ENDIF
         call distributeKgrad(OUTPUT%resultTensor,integral%integralsABCD,input%LST_DLHS,&
              & input%LST_DRHS,input%NDMAT_RHS,PQ,Input,output,LUPRI,IPRINT)
        ELSEIF(INPUT%fullcontraction)THEN
           call distributeKcont(&
                & integral%integralsABCD,input%LST_DLHS,input%LST_DRHS,&
                & input%NDMAT_RHS,PQ,Input,output,Integral%Econt,LUPRI,IPRINT)
        ELSE
           call distributeExchange(OUTPUT%resultTensor,integral%integralsABCD,input%LST_DRHS,&
                & input%NDMAT_RHS,INPUT%exchangeFactor,PQ,Input,output,LUPRI,IPRINT)
        ENDIF
       ENDIF
      ELSE
       ! No pre-contraction with densities (either full integrals, or contraction
       ! with full integrals)         
       IF(INPUT%DO_MULMOM)THEN
        nMAT    = INPUT%nMultipoleMomentComp
        mmorder = INPUT%MMORDER
        nSPHMAT = (mmorder+1)*(mmorder+1)
        nA      = PQ%P%p%orbital1%totOrbitals
        nB      = PQ%P%p%orbital2%totOrbitals
        nDeriv  = PQ%P%p%ngeoDerivComp * PQ%Q%p%ngeoDerivComp
        nPasses = PQ%P%p%nPasses * PQ%Q%p%nPasses
        !$OMP CRITICAL (distributeMM)
        CALL printMMtoFile(PQ,Integral%integralsABCD,nA,nB,nMAT,nSPHMAT,nDeriv,nPasses,INPUT,OUTPUT,LUPRI,IPRINT)
        !$OMP END CRITICAL (distributeMM)
       ELSE
        IF (output%dograd) THEN
         CALL distributePQgrad(OUTPUT%resultTensor,&
              & PQ,Integral%integralsABCD,dimQ,dimP,Input,output,LUPRI,IPRINT)
        ELSE
         IF (output%decpacked)then
          CALL Explicit4centerDEC(OUTPUT%resultMat,output%ndim(1),&
               & output%ndim(2),output%ndim(3),output%ndim(4),PQ,&
               & Integral%integralsABCD,dimQ,dimP,Input,output,LUPRI,IPRINT)
         ELSEIF (output%decpacked2)then
            CALL Explicit4centerDEC2(OUTPUT%resultMat,output%ndim(1),&
                 & output%ndim(2),output%ndim(3),output%ndim(4),PQ,&
                 & Integral%integralsABCD,dimQ,dimP,Input,output,LUPRI,IPRINT)
         ELSEIF (output%decpackedK)then
            CALL Explicit4centerDECK(OUTPUT%resultMat,output%ndim(1),&
                 & output%ndim(2),output%ndim(3),output%ndim(4),PQ,&
                 & Integral%integralsABCD,dimQ,dimP,Input,output,LUPRI,IPRINT)
         ELSEIF(output%FullAlphaCD)THEN
            call Explicit3CenterDEC(OUTPUT%resultMat,OUTPUT%result3D,&
                 & output%ndim(1),output%ndim(3),output%ndim(4),&
                 & output%ndim3D(1),output%ndim3D(2),output%ndim3D(3),&
                 & PQ,Integral%integralsABCD,dimQ,dimP,&
                 & Input,output,LUPRI,IPRINT)
         ELSE
            !default simple distribute without contraction
            CALL GeneraldistributePQ(OUTPUT%resultTensor,PQ,&
                 & Integral%integralsABCD,dimQ,dimP,Input,output,LUPRI,IPRINT)
         ENDIF
        ENDIF
       ENDIF
      ENDIF
   ENDIF
ENDIF

END SUBROUTINE DistributeIntegrals

!> \brief build the integrand structure PQ from P and Q overlap distributions
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param PQ contains info about the integrand, to be built. 
!> \param P the left hand side overlap distribution
!> \param Q the right hand side overlap distribution
!> \param Input contains info about the requested integral evaluation
!> \param ILHS the left hand side ODbatch index 
!> \param IRHS the right hand side ODbatch index 
!> \param IUNIT the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE Build_Integrand(PQ,P,Q,INPUT,ILHS,IRHS,IUNIT,IPRINT)
TYPE(Integrand)       :: PQ
TYPE(Overlap),target  :: P,Q
TYPE(IntegralInput)   :: INPUT
Integer               :: IUNIT,IPRINT,ILHS,IRHS
!
Integer     :: ipq,ip,iq,idir
Real(realk) :: d2, pref
Real(realk), parameter :: PI=3.14159265358979323846E0_realk, Nill = 0.0E0_realk, OneHalf=1.5E0_realk
Real(realk), parameter :: Two = 2.0E0_realk, TwoHalf=2.5E0_realk 
logical :: coulomb,nucpotLHS,nucpotRHS
!
integer :: i
 PQ%Operator = INPUT%operator
 PQ%kinetic  = PQ%Operator.EQ.KineticOperator
 IF(Q%magderiv.EQ.1)THEN 
    !we always differentiate in contract_P so we interchange the order
    PQ%P%p            => Q
    PQ%Q%p            => P
    PQ%reverseOrder = .TRUE.
 ELSE
    PQ%P%p            => P
    PQ%Q%p            => Q
    PQ%reverseOrder = .FALSE.
 ENDIF
!simen Refine!
 PQ%nPrimitives    = P%nPrimitives*P%nPasses   * Q%nPrimitives*Q%nPasses
 PQ%maxContracted  = P%maxContracted * Q%maxContracted
 PQ%minAngmom      = P%minAngmom     + Q%minAngmom
 PQ%maxAngmom      = P%maxAngmom     + Q%maxAngmom
 PQ%startAngmom    = P%startAngmom   + Q%startAngmom 
 PQ%endAngmom      = PQ%maxAngmom    + INPUT%geoDerivOrder + input%CMorder &
     &               + INPUT%magderOrderP + INPUT%magderOrderQ

 PQ%samePQ         = INPUT%sameODs .AND. (ILHS .EQ. IRHS)
 PQ%nAngmom        = P%nAngmom       * Q%nAngmom 
 IF (PQ%samePQ) THEN
  PQ%nAngmom = PQ%P%p%nAngmom*(PQ%P%p%nAngmom+1)/2
!  PQ%nPrimitives = P%nPrimitives*(P%nPrimitives+1)/2
 ENDIF
!Dimensions according to nPrimitives
coulomb = ((PQ%Operator.EQ.CoulombOperator).OR.(PQ%Operator.EQ.ErfOperator)).OR.&
     & ((PQ%Operator.EQ.GGemOperator).OR.(PQ%Operator.EQ.CAMOperator)).OR.&
     & (PQ%Operator.EQ.ErfcOperator).OR.(PQ%Operator.EQ.GGemCouOperator).OR.&
     & (PQ%Operator.EQ.GGemGrdOperator).OR.(PQ%Operator.EQ.GGemQuaOperator)
nucpotLHS = (PQ%Operator.EQ.NucpotOperator).AND.PQ%P%p%TYPE_Nucleus
nucpotRHS = (PQ%Operator.EQ.NucpotOperator).AND.PQ%Q%p%TYPE_Nucleus
!CALL SetIntegrand(PQ%iprimP,PQ%iprimQ,PQ%distance,PQ%squaredDistance,&
CALL SetIntegrand(PQ%distance,PQ%squaredDistance,&
     &PQ%exponents,PQ%reducedExponents,coulomb,nucpotRHS,nucpotLHS,PQ%integralPrefactor,&
     &PQ%nPrimitives,PQ%P%p%exponents,PQ%Q%p%exponents,PQ%P%p%center,PQ%Q%p%center,&
     &PQ%P%p%nPrimitives*PQ%P%p%nPasses,PQ%Q%p%nPrimitives*PQ%Q%p%nPasses,&
     &PQ%Q%p%nPrimitives,PQ%Q%p%nPasses,IUNIT)

IF(INPUT%ATTFACTOR) CALL DSCAL(PQ%nPrimitives,1E0_realk+INPUT%ATTOMEGA/PI,PQ%integralPrefactor,1)

 IF (IPRINT.GT. 15) THEN
    IF(PQ%reverseOrder)THEN
       WRITE(IUNIT,*)'The LHS and RHS will be used in reverse order'
       WRITE(IUNIT,*)'The here printet LHS overlap, will be used as the RHS overlap'
    ENDIF
    CALL PRINT_OVERLAP(P,IUNIT,IPRINT,'LHS')
    CALL PRINT_OVERLAP(Q,IUNIT,IPRINT,'RHS') 
    CALL PRINT_Integrand(PQ,IUNIT,INPUT)
 ENDIF
END SUBROUTINE Build_Integrand

!> \brief set integrand wrapper routine which branch out depending on the requested order of primitive AOs 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param rPQ x,y and z distance between overlap P and Q, (3,nPrimP*nPrimQ)
!> \param squaredDistance the squared distance between overlap P and Q
!> \param exponents the exponent = p + q (p = exponent of LHS overlap P,etc) 
!> \param reducedExponents the reduced exponent p*q/p_q
!> \param coulomb flag which tells if it is a coulomb repulsion integral
!> \param nucpot flag which tells if it is a electron-nuclear attraction integral
!> \param integralprefactor the integralprefactor which depend on the type of integral 
!> \param npq the number of primitives AOs (nprimitives for P * nprimitives for Q)
!> \param pexp the exponents of overlap distribution P
!> \param qexp the exponents of overlap distribution P
!> \param pcent the centers of overlap distribution P
!> \param qcent the centers of overlap distribution P
!> \param np nprimitives*nPasses for P
!> \param nq nprimitives*nPasses for Q
!> \param nPrimQ nprimitives for Q
!> \param nPassQ number of Passes for Q
!> \param IUNIT the logical unit number for the output file
SUBROUTINE SetIntegrand(rPQ,squaredDistance,exponents,&
     &                  reducedExponents,coulomb,nucpotRHS,nucpotLHS,integralPrefactor,npq,&
     &                  pexps,qexps,pcent,qcent,nPrimPassP,nPrimPassQ,nprimQ,nPassQ,IUNIT)
implicit none
integer,intent(in)     :: nPrimPassP,nPrimPassQ,npq,IUNIT,nprimQ,nPassQ
real(realk),intent(inout) :: rPQ(nPrimPassQ,nPrimPassP,3)
real(realk),intent(inout) :: exponents(nPrimQ,nPassQ,nPrimPassP)
real(realk),intent(inout) :: squaredDistance(nPrimPassQ,nPrimPassP)
real(realk),intent(inout) :: reducedExponents(nPrimQ,nPassQ,nPrimPassP)
real(realk),intent(inout) :: integralPrefactor(nPrimQ,nPassQ,nPrimPassP)
real(realk),intent(in) :: pexps(nPrimPassP),pcent(3,nPrimPassP)
real(realk),intent(in) :: qexps(nPrimPassQ),qcent(3,nPrimPassQ)
logical,intent(in)     :: nucpotRHS,nucpotLHS,coulomb
!
real(realk) :: px,py,pz,qx,qy,qz,pqx,pqy,pqz
integer :: ipq,ip,iq,idir,startq,ipassQ,tmpipq
real(realk) :: pqdist(3),pqdist2,d2,p,q,p_q,pq_inv
Real(realk), parameter :: PI=3.14159265358979323846E0_realk, Nill = 0.0E0_realk, OneHalf=1.5E0_realk
Real(realk), parameter :: Two = 2.0E0_realk, TwoHalf=2.5E0_realk 
Real(realk), parameter :: PIFAC = 34.986836655249725E0_realk !Two*PI**TwoHalf
Real(realk), parameter :: TWOPI = 6.2831853071795862E0_realk 
! Too much stack memory usage
!Real(realk) :: expoTMP(nPrimQ,np),expoTMP2(nPrimQ,np)
!Real(realk) :: invexponents(nPrimQ,np)
!Real(realk) :: sqrtexpoTMP(nPrimQ,np)
!

!Primitives ordered according to (Q,P)
!exponents(nPrimQ,nPassQ,nP)

!use rPQ and squaredDistance as temporary array
IF (coulomb) THEN
   !   exponents not needed for Coulomb
   DO ip=1, nPrimPassP
      p = pexps(ip)
      DO iq=1, nPrimQ
         rPQ(iq,ip,2) = p * qexps(iq)
         rPQ(iq,ip,3) = 1.0E0_realk/(p + qexps(iq))
         squaredDistance(iq,ip) = SQRT(p + qexps(iq))
      ENDDO
   ENDDO
ELSE
   DO ip=1, nPrimPassP
      p = pexps(ip)
      DO iq=1, nPrimQ
         exponents(iq,1,ip) = p + qexps(iq)
         rPQ(iq,ip,2) = p * qexps(iq)
         rPQ(iq,ip,3) = 1.0E0_realk/(p + qexps(iq))
      ENDDO
   ENDDO
ENDIF

!WE BUILD THE PART WHICH IS THE SAME FOR ALL PASSES
DO ip=1, nPrimPassp
   p = pexps(ip)
   DO iq=1, nPrimQ
      reducedExponents(iq,1,ip) = rPQ(iq,ip,2)*rPQ(iq,ip,3)
   ENDDO
ENDDO
IF (coulomb) THEN
   DO ip=1, nPrimPassp
      DO iq=1, nPrimQ
         integralPrefactor(iq,1,ip) = PIFAC/(rPQ(iq,ip,2)*squaredDistance(iq,ip))
      ENDDO
   ENDDO
ELSEIF(nucpotRHS)THEN
   DO ip=1, nPrimPassp
      p  = pexps(ip)
      DO iq=1, nPrimQ
         integralPrefactor(iq,1,ip) = TWOPI/p
      ENDDO
   ENDDO
ELSEIF(nucpotLHS)THEN
   DO iq=1, nPrimQ
      q  = qexps(iq)
      DO ip=1, nPrimPassP
         integralPrefactor(iq,1,ip) = TWOPI/q
      ENDDO
   ENDDO
ELSE
   DO ip=1, nPrimPassp
      DO iq=1, nPrimQ
         integralPrefactor(iq,1,ip) = (PI*rPQ(iq,ip,3))**(OneHalf)
      ENDDO
   ENDDO
ENDIF
!COPY TO ALL PASSES
IF(npassQ.GT.1)then
   IF(.NOT.coulomb)THEN
      !not needed for coulomb
      DO ip=1, nPrimPassp
         do ipassQ=2,npassQ
            DO iq=1, nPrimQ
               exponents(iq,ipassQ,ip) = exponents(iq,1,ip)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
   DO ip=1, nPrimPassp
      do ipassQ=2,npassQ
         DO iq=1, nPrimQ
            reducedExponents(iq,ipassQ,ip) = reducedExponents(iq,1,ip) 
         ENDDO
      ENDDO
   enddo
   DO ip=1, nPrimPassp
      do ipassQ=2,npassQ
         DO iq=1, nPrimQ
            integralPrefactor(iq,ipassQ,ip) = integralPrefactor(iq,1,ip)
         ENDDO
      ENDDO
   ENDDO
ENDIF
!BUILD THE PART WHICH IS UNIQUE FOR EACH OVERLAP
!ipq = 0
DO ip=1, nPrimPassp
   px = pcent(1,ip)
   py = pcent(2,ip)
   pz = pcent(3,ip)
   DO iq=1, nPrimPassq
      pqx = px - qcent(1,iq)
      pqy = py - qcent(2,iq)
      pqz = pz - qcent(3,iq)

      rPQ(iq,ip,1) = pqx
      rPQ(iq,ip,2) = pqy
      rPQ(iq,ip,3) = pqz
       
      squaredDistance(iq,ip) = pqx*pqx+pqy*pqy+pqz*pqz
   ENDDO
!   ipq = ipq+nq
ENDDO

END SUBROUTINE SetIntegrand

!> \brief Print the integrand structure PQ
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param PQ contains info about the integrand, to be built. 
!> \param IUNIT the logical unit number for the output file
!> \param INPUT contains info about the requested integral evaluation
SUBROUTINE PRINT_Integrand(PQ,IUNIT,INPUT)
TYPE(Integrand)       :: PQ
Integer               :: IUNIT
TYPE(IntegralInput)   :: INPUT

!
Integer :: i
!
CALL LSHEADER(IUNIT,'Integrand output')
WRITE(IUNIT,'(1X,A)') ''
WRITE(IUNIT,'(1X,A)') '*************************'
WRITE(IUNIT,'(1X,A)') '***     Integrand     ***'
WRITE(IUNIT,'(1X,A)') '*************************'
WRITE(IUNIT,'(1X,A)') ''
WRITE(IUNIT,'(5X,A,I3)') 'Operator      =', PQ%Operator
WRITE(IUNIT,'(5X,A,I6)') 'nAngmom       =', PQ%nAngmom
WRITE(IUNIT,'(5X,A,I6)') 'nPrimitives   =', PQ%nPrimitives
WRITE(IUNIT,'(5X,A,I6)') 'maxContracted =', PQ%maxContracted
WRITE(IUNIT,'(5X,A,I6)') 'minAngmom     =', PQ%minAngmom
WRITE(IUNIT,'(5X,A,I6)') 'maxAngmom     =', PQ%maxAngmom
WRITE(IUNIT,'(5X,A,I6)') 'startAngmom   =', PQ%startAngmom
WRITE(IUNIT,'(5X,A,I6)') 'endAngmom     =', PQ%endAngmom
WRITE(IUNIT,'(3X,A)')    '----------------------------------------------------------------------------------'
!WRITE(IUNIT,'(3X,A)')    '  prim  iP  iQ   X_PQ    Y_PQ    Z_PQ   R_PQ^2     exp     redExp     intPre     '
WRITE(IUNIT,'(3X,A)')    '  prim        X_PQ        Y_PQ        Z_PQ      R_PQ^2         exp      redExp       intPre'
DO i=1,PQ%nPrimitives
  WRITE(IUNIT,'(5X,I4,7ES12.4)') i,&
       & PQ%distance(i,1),PQ%distance(i,2),PQ%distance(i,3),&
       & PQ%squaredDistance(i),PQ%exponents(i),PQ%reducedExponents(i),PQ%integralPrefactor(i)
ENDDO
WRITE(IUNIT,'(3X,A)')    '----------------------------------------------------------------------------------'

END SUBROUTINE PRINT_Integrand

!> \brief contract the LHS overlap P (primitive to contracted functions, Ecoeff contraction,..)
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param INTEGRAL storage of integrals and intermediates
!> \param PQ contains info about the integrand, to be built. 
!> \param INPUT contains info about the requested integral evaluation
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE Contract_P(INTEGRAL,PQ,input,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: INTEGRAL
TYPE(IntegralInput):: INPUT
TYPE(Integrand)    :: PQ
Integer            :: LUPRI,IPRINT
!
Integer  :: iAngmom,startMomQ,iOrbital,ideriv,nderiv
!
REAL(REALK)          :: CPU1,CPU2,WALL1,WALL2
Real(realk),dimension(:),pointer :: ptemp
!Integer :: ISAVE
!
!Loop over angular contributions sharing the set of primitive functions
DO iAngmom=1,PQ%P%p%nAngmom
   ! Distribute Hermite 2-el integrals WTUV to <tuv1|w|tuv2>
   CALL DistributeHermiteP(Integral,PQ,iAngmom,LUPRI,IPRINT,input)
   !Output(nPrimP,nOrbQ,ntuvP)'
   CALL ContractEcoeffP(Integral,PQ%P%p,iAngmom,LUPRI,IPRINT)
   !Output(nPrimP,nOrbQ,ijkP)
   CALL ContractBasis(Integral,PQ%P%p,iAngmom,LUPRI,IPRINT)
   !Output(nContP,nOrbQ,ijkP)
   CALL extractDifferentiated(Integral,PQ%P%p,iAngmom,integral%lhsGeoOrder,Integral%lhsGeoComp,LUPRI,IPRINT)
   IF (PQ%P%p%magderiv.EQ.1)CALL extractMagDifferentiated(Integral,PQ%P%p,iAngmom,LUPRI,IPRINT)
   CALL SphericalTransform(Integral,PQ%P%p,iAngmom,LUPRI,IPRINT)
   !Output(nContP,nOrbQ,lmP,nderivP)
   CALL AddToPQ(Integral,PQ,iAngmom,input%nMultipoleMomentComp,input%contAng,LUPRI,IPRINT)
   !Output(nContP,nOrbQ,nDerivP,lmP)
ENDDO

END SUBROUTINE Contract_P

!> \brief Extracts the directional derivative components from the Hermite primitive integrals components
!> \author S. Reine
!> \date 2010-03-07
!> \param integral Contains the information about the integrals
!> \param P Contains the information about the product overlap 
!> \param iAngmom The angular component of the overlap
!> \param nComp The number of derivative components
!> \param LUPRI Default output pint unit
!> \param IPRINT Print level (0 no output - high lots of output)
SUBROUTINE extractDifferentiated(Integral,P,iAngmom,iOrder,nComp,LUPRI,IPRINT)
implicit none
TYPE(Integralitem),intent(INOUT) :: INTEGRAL
TYPE(Overlap),intent(IN)         :: P
Integer,intent(IN)               :: iAngmom,iOrder,nComp,LUPRI,IPRINT
!
Integer :: iA1,iA2,l1,l2,ijk,lm,ijkdiff,ijkcart,dim1,ijk1,ijk2
Real(realk),dimension(:),pointer :: ptemp
!
IF (iOrder.EQ.0) RETURN  !No extraction for the undifferentiated case

iA1 = P%indexAng1(iangmom)
iA2 = P%indexAng2(iangmom)
l1 = P%orbital1%angmom(iA1)
l2 = P%orbital2%angmom(iA2)
CALL GET_IJK(l1,l2,ijk1,ijk2,lm,ijkdiff,P%sphericalEcoeff,iOrder,P%single)
ijk = ijk1*ijk2
dim1 = Integral%nOrb*Integral%nPrim

#ifdef VAR_LSDEBUGINT
IF(dim1*ijkdiff.GT.allocIntmaxTUVdim)THEN
   print*,'dim1',dim1
   print*,'ijkdiff',ijkdiff
   print*,'allocIntmaxTUVdim',allocIntmaxTUVdim
   call lsquit('alloc error',lupri)
ENDIF
IF(dim1*P%ngeoDerivComp*P%nCartesianMomentComp*ijk.GT.allocIntmaxTUVdim)THEN
   print*,'dim1',dim1
   print*,'ijk',ijk
   print*,'P%nGeoDerivComp',P%nGeoDerivComp
   print*,'P%nCartesianMomentComp',P%nCartesianMomentComp
   print*,'allocIntmaxTUVdim',allocIntmaxTUVdim
   call lsquit('alloc error : TK1',lupri)
ENDIF
#endif

CALL extractDifferentiated_PA(Integral%IN,INTEGRAL%OUT,l1,l2,ijkdiff,ijk,dim1,&
     & nComp*P%nCartesianMomentComp,iOrder,P%single,P%orbital1%TYPE_empty,P%orbital2%TYPE_empty,&
     & LUPRI,IPRINT)
!
Integral%nAng = ijk
Integral%ngeoDeriv = nComp*P%nCartesianMomentComp

!Swap pointers IN and OUT
ptemp => Integral%IN
Integral%IN  => Integral%OUT
Integral%OUT => ptemp
!
END SUBROUTINE extractDifferentiated

!> \brief Extracts the directional derivative components from the Hermite primitive integrals components
!> \author S. Reine
!> \date 2010-03-07
!> \param HermiteDiff Differentiated integrals in Hermite-primitive basis (including expoentital perfactors)
!> \param DirectionalDiff Extracted Hermite components in the different derivative directions
!> \param l1 Angular momentum of orbital 1
!> \param l2 Angular momentum of orbital 2
!> \param ijkdiff Number of Hermite differentiated ijk-components
!> \param ijk Number of undifferentiated ijk-components
!> \param dim1 The number of contracted functions multiplied by the number of orbitals (of the other electron)
!> \param derivComp The number of directional derivative components
!> \param geoderivorder The order of differentiation (of current electron)
!> \param single True if only one AO
!> \param empty1 True if orbital 1 is type_empty
!> \param empty2 True if orbital 2 is type_empty
!> \param LUPRI Default output pint unit
!> \param IPRINT Print level (0 no output - high lots of output)
SUBROUTINE extractDifferentiated_PA(HermiteDiff,DirectionalDiff,l1,l2,ijkdiff,ijk,dim1,&
     &                              derivComp,geoderivorder,single,empty1,empty2,LUPRI,IPRINT)
implicit none
Integer,intent(IN)      :: l1,l2,ijkdiff,dim1,ijk,derivComp,geoderivorder
Real(realk),intent(IN)  :: HermiteDiff(dim1,ijkdiff)
Real(realk),intent(INOUT) :: DirectionalDiff(dim1,derivComp,ijk)
Logical,intent(IN)      :: single,empty1,empty2
Integer,intent(IN)      :: LUPRI,IPRINT
!
Integer :: increment(2),icent,ncent,iComp
Integer :: iX,iY,iZ,jX,jY,jZ
Integer :: ijkIndex(0:l1+geoderivorder,0:l1+geoderivorder,0:l1+geoderivorder,&
     & 0:l2+geoderivorder,0:l2+geoderivorder,0:l2+geoderivorder)
Integer :: i,iDiff,iDim,P1,iP1,jP1,kP1,P2,iP2,jP2,kP2

ncent = 1 + geoderivorder
IF (single) ncent = 1

iComp = 0
iDiff=0
DO icent=1,ncent
  increment(1) = geoderivorder+1-icent
  increment(2) = icent - 1
  IF (single) THEN
    IF (empty1) THEN
      increment(1) = 0
      increment(2) = geoderivorder
    ELSEIF (empty2) THEN
      increment(1) = geoderivorder
      increment(2) = 0
    ELSE
      CALL LSQUIT('Error in extractDifferentiated_PA. single and not empty1 or empty2',-1)
    ENDIF
  ENDIF

  P2 = l2+increment(2)
  P1 = l1+increment(1)
! Set up differentiated ijk index
  DO iP2=P2,0,-1
    DO jP2=P2-iP2,0,-1
      kP2 = P2-iP2-jP2
      DO iP1=P1,0,-1
        DO jP1=P1-iP1,0,-1
          kP1=P1-iP1-jP1
          iDiff = iDiff+1
          ijkIndex(iP1,jP1,kP1,iP2,jP2,kP2) = iDiff
        ENDDO
      ENDDO
    ENDDO
  ENDDO

! Loop over the three Cartesian directions for both orbitals
      DO iX=increment(1),0,-1
        DO iY=increment(1)-iX,0,-1
          iZ = increment(1)-iX-iY
  DO jX=increment(2),0,-1
    DO jY=increment(2)-jX,0,-1
      jZ=increment(2)-jX-jY
          iComp = iComp+1
!         Regular ijk-loop
          i = 0
          DO iP2=l2,0,-1
            DO jP2=l2-iP2,0,-1
              kP2 = l2-iP2-jP2
              DO iP1=l1,0,-1
                DO jP1=l1-iP1,0,-1
                  kP1=l1-iP1-jP1
                  i=i+1
                  iDiff = ijkIndex(iP1+iX,jP1+iY,kP1+iZ,iP2+jX,jP2+jY,kP2+jZ)
                  DO iDim=1,dim1
                    DirectionalDiff(iDim,iComp,i) = HermiteDiff(iDim,iDiff)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
!         Regular ijk-loop ends here
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO
!
IF (IPRINT.GT. 20) THEN
  CALL PrintTensor(DirectionalDiff,'Directional deriv.  ',dim1,derivComp,ijk,&
     &lupri,'dim1  ','deriv ','ijk   ',3)
ENDIF
!
END SUBROUTINE extractDifferentiated_PA

!> \brief Extracts the directional derivative components from the Hermite primitive integrals components
!> \author S. Reine
!> \date 2010-03-07
!> \param integral Contains the information about the integrals
!> \param P Contains the information about the product overlap 
!> \param iAngmom The angular component of the overlap
!> \param LUPRI Default output pint unit
!> \param IPRINT Print level (0 no output - high lots of output)
SUBROUTINE extractMagDifferentiated(Integral,P,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(Integralitem),intent(INOUT) :: INTEGRAL
TYPE(Overlap),intent(IN)         :: P
Integer,intent(IN)               :: iAngmom,LUPRI,IPRINT
!
Real(realk),dimension(:),pointer :: ptemp
!
call extractMagDifferentiated1(Integral%OUT,INTEGRAL%IN,P%distance12,&
        &Integral%nAng/3,Integral%nOrb,Integral%nPrim,LUPRI,IPRINT) 
Integral%nAng = Integral%nAng/3
!Swap pointers IN and OUT
ptemp => Integral%IN
Integral%IN  => Integral%OUT
Integral%OUT => ptemp

END SUBROUTINE ExtractMagDifferentiated

!> \brief distribute the calculated integrals into magnetic derivative integrals (performs the crossproduct required)
!> \author T. Kjaergaard
!> \date 2010-25-10
!> \param MagDerivInt the actual integrals (ab|x|cd),(ab|y|cd),(ab|z|cd) 
!> \param CARTINT a temporary matrix used to build the cross-product
!> \param distance12 the Xab,Yab,Zab distance between the two atomic orbitals
!> \param nAngP number of angular components for overlap P
!> \param nQ total number of orbitals on overlap Q
!> \param nContP number of primitive*npasses for overlap P
!> \param LUPRI the logical unit number for the output file
subroutine extractMagDifferentiated1(MagDerivInt,CARTINT,distance12,nAngP,nOrbQ,nContP,LUPRI,IPRINT)
implicit none
integer :: nAngP,nOrbQ,nContP,LUPRI,IPRINT
real(realk) :: MagDerivInt(nContP*nOrbQ,3,nAngP)
real(realk) :: CARTINT(nContP*nOrbQ,nAngP,3)
real(realk) :: distance12(3) !(3,npasses)
!
integer,parameter :: betalist(3)=(/2,3,1/),gammalist(3)=(/3,1,2/)
integer :: X,beta,gamma,iAng,iContP
real(realk) :: BETADIST,GAMMADIST

DO X=1,3
   beta = betalist(X)
   gamma = gammalist(X)
   BETADIST=distance12(beta)
   GAMMADIST=distance12(gamma)   
   DO iAng=1,nAngP
      DO iContP=1,nContP*nOrbQ
         MagDerivInt(iContP,X,iAng)= BETADIST*CARTINT(iContP,iAng,GAMMA)&
              &-GAMMADIST*CARTINT(iContP,iAng,BETA)
      ENDDO
   ENDDO
ENDDO
IF(IPRINT.GT.100)THEN
   call printtensor(MagDerivInt,'MAGDERIVCROSSPRODOCT',nContP*nOrbQ,3,nAngP,&
        &lupri,'nP    ','coor  ','nAngP ',2)
ENDIF
end subroutine extractMagDifferentiated1

!> \brief wrapper routine to add the calculated integral into integral%PQ for storage until distribution
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param INTEGRAL storage of integrals and intermediates
!> \param PQ contains info about the integrand, to be built. 
!> \param iAngmom the index of angular momentum 
!> \param iOrbital the current index in the integral%PQ array  
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE AddToPQ(Integral,PQ,iAngmom,nMultipoleMomentComp,contAng,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Integrand)    :: PQ
Integer            :: iAngmom,nMultipoleMomentComp,LUPRI,IPRINT
logical,intent(IN) :: contAng
!
Integer :: nQ,nAng,nCont,ndim5P,startDer
Integer :: angA,angB,nA,nB,nContA,nContB,nCompA,nCompB,startA,startB,nPassP
Integer :: nC,nD,nPassQ,ndim5Q,nDer
TYPE(Overlap),pointer :: Q,P
Logical :: permute

Q => PQ%Q%p
P => PQ%P%p
nQ      = Q%totOrbitals(integral%rhsGeoOrder+1) *nMultipoleMomentComp
nC      = Q%orbital1%totOrbitals  *nMultipoleMomentComp
nD      = Q%orbital2%totOrbitals
nPassQ  = Q%nPasses
IF (PQ%Q%p%TYPE_FTUV) nPassQ = 1

startDer = Integral%startDerivativeIndex
ndim5Q   = integral%rhsGeoComp*Q%nCartesianMomentComp
ndim5P   = integral%lhsGeoComp*P%nCartesianMomentComp

nDer     = Integral%nDerivComp*Q%nCartesianMomentComp*P%nCartesianMomentComp

IF(P%magderiv.EQ.1)ndim5P = ndim5P*3
IF(P%magderiv.EQ.1) nDer = nDer*3

nAng    = P%nOrbComp(iAngmom)
nCont   = P%nContracted(iAngmom)*P%nPasses

angA    = P%indexAng1(iAngmom)
angB    = P%indexAng2(iAngmom)
nA      = P%orbital1%totOrbitals
nContA  = P%orbital1%nContracted(angA)
startA  = P%orbital1%startLocOrb(angA)
nCompA  = P%orbital1%nOrbComp(angA)
nB      = P%orbital2%totOrbitals
nContB  = P%orbital2%nContracted(angB)
startB  = P%orbital2%startLocOrb(angB)
nCompB  = P%orbital2%nOrbComp(angB)
nPassP  = P%nPasses
permute = P%sameAO.AND.(startA.NE.startB)

IF (IPRINT.GT.20) THEN
  CALL print_addPQ(Integral%IN,nC*nD,nPassQ*ndim5Q,ndim5P,nAng,nCont,permute,contAng,lupri,iprint)
ENDIF

IF (.NOT.contAng) THEN !Default AO ordering: angular,contracted
  IF (nQ.EQ. 1) THEN
#ifdef VAR_LSDEBUGINT
  IF(nA*nB*ndim5P*nPassP.GT.allocIntmaxTUVdim)THEN
     print*,'nA',nA
     print*,'nB',nB
     print*,'ndim5P',ndim5P
     print*,'nPassP',nPassP
     print*,'allocIntmaxTUVdim',allocIntmaxTUVdim
     call lsquit('AddToPQ1 alloc error nA*nB*ndim5P*nPassP > allocIntmaxTUVdim',lupri)
  ENDIF
  IF(nCont*ndim5P*nAng.GT.allocIntmaxTUVdim)THEN
     print*,'nCont  ',nCont
     print*,'ndim5P',ndim5P
     print*,'nAng   ',nAng
     print*,'allocIntmaxTUVdim',allocIntmaxTUVdim
     call lsquit('AddToPQ1 alloc error nCont*nAng*ndim5P > allocIntmaxTUVdim',lupri)
  ENDIF
#endif

    CALL AddToPQ1(Integral%integralsABCD,Integral%IN,nA,nB,nCompA,nCompB,&
         &        nContA,nContB,startA,startB,nAng,nCont,nPassP,ndim5P,permute,startDer,nDer,LUPRI,IPRINT)
  ELSE
#ifdef VAR_LSDEBUGINT
  IF(nC*nD*nA*nB*ndim5Q*ndim5P*nPassP*nPassQ.GT.allocIntmaxTUVdim)THEN
     print*,'nQ',nQ
     print*,'nA',nA
     print*,'nB',nB
     print*,'ndim5Q',ndim5Q
     print*,'ndim5P',ndim5P
     print*,'nPassP',nPassP
     print*,'nPassQ',nPassQ
     print*,'nQ*nA*nB*ndim5Q*ndim5P*nPassP*nPassQ',&
          & nQ*nA*nB*ndim5Q*ndim5P*nPassP*nPassQ
     print*,'allocIntmaxTUVdim',allocIntmaxTUVdim
     call lsquit('AddToPQgen alloc error nQ*nA*nB*ndim5Q*ndim5P*nPassP*nPassQ > allocIntmaxTUVdim',lupri)
  ENDIF
  IF(nCont*nC*nD*nPassQ*ndim5Q*ndim5P*nAng.GT.allocIntmaxTUVdim)THEN
     print*,'ndim5Q',ndim5Q
     print*,'nPassQ',nPassQ
     print*,'nQ',nQ
     print*,'nCont  ',nCont
     print*,'ndim5P',ndim5P
     print*,'nAng   ',nAng
     print*,'allocIntmaxTUVdim',allocIntmaxTUVdim
     call lsquit('AddToPQgen alloc error nCont*nQ*nPassQ*ndim5Q*ndim5P*nAng > allocIntmaxTUVdim',lupri)
  ENDIF
#endif
    CALL AddToPQgen(Integral%integralsABCD,Integral%IN,nC*nD,nPassQ,ndim5Q,nA,nB,nCompA,nCompB,&
         &        nContA,nContB,startA,startB,nAng,nCont,nPassP,ndim5P,permute,startDer,nDer,LUPRI,IPRINT)
  ENDIF
ELSE !Optional AO ordering: contracted first then angular components (only applicable to generally
     !contracted basis sets)
  CALL AddToPQgen_ca(Integral%integralsABCD,Integral%IN,nC*nD,nPassQ,ndim5Q,nA,nB,nCompA,nCompB,&
         &           nContA,nContB,startA,startB,nAng,nCont,nPassP,ndim5P,permute,startDer,nDer,LUPRI,IPRINT)
ENDIF

IF (IPRINT.GT. 10) THEN
  CALL PrintPQ(Integral%integralsABCD,nA,nB,nC,nD,ndim5P,ndim5Q,nPassP,nPassQ,nDer,LUPRI,IPRINT)
ENDIF

CONTAINS
  SUBROUTINE print_addPQ(AddPQ,nQ,n5Q,nAng,nCont,n5P,permute,contang,LUPRI,IPRINT)
  implicit none
  INTEGER,intent(IN)     :: nQ,n5Q,nAng,nCont,n5P,LUPRI,IPRINT
  LOGICAL,intent(IN)     :: permute,contang
  Real(realk),intent(IN) :: AddPQ(nCont,nQ,n5Q,n5P,nAng)
  !
  Integer :: iCont,iQ,i5Q,i5P,iAng
  
  write(lupri,'(1X,A)')    'Printing AddPQ entering AddToPQ'
  write(lupri,'(3X,A10,I5)') 'n5Q    =',n5Q
  write(lupri,'(3X,A10,I5)') 'n5P    =',n5P
  write(lupri,'(3X,A10,I5)') 'nQ     =',nQ
  write(lupri,'(3X,A10,I5)') 'nAngP  =',nAng
  write(lupri,'(3X,A10,I5)') 'nContP =',nCont
  write(lupri,'(3X,A10,L5)') 'permute =',permute
  write(lupri,'(3X,A10,L5)') 'contang =',contang
  
  write(lupri,'(5X,A)') '  iQ5  iP5   iQ  iAng   iCont=1,...,nCont'
  DO i5Q = 1,n5Q
  DO i5P = 1,n5P
  DO iQ  = 1,nQ
  DO iAng = 1,nAng
    write(lupri,'(5X,4I5,10F15.9/(25X,10F15.9))') i5Q,i5P,iQ,iAng,(AddPQ(iCont,iQ,i5Q,i5P,iAng),iCont=1,nCont)
  ENDDO
  ENDDO
  ENDDO
  ENDDO
  END SUBROUTINE print_addPQ
END SUBROUTINE AddToPQ

 
!> \brief add the calculated integral into integral%PQ for storage until distribution
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param CDAB  the integrals to be returned
!> \param addPQ the integrals in a form suitable for the inner contraction
!> \param nQ the number of AO functions on the RHS overlap Q
!> \param nA the number of AO functions of the first LHS orbital
!> \param nB the number of AO functions of the second LHS orbital
!> \param nCompA the number of angular components of the first LHS orbital for given iAngmom
!> \param nCompB the number of angular components of the second LHS orbital for given iAngmom
!> \param nContA the number of contracted functions of the first LHS orbital for given iAngmom
!> \param nContB the number of contracted functions of the second LHS orbital for given iAngmom
!> \param startA the local starting orbital index of the first LHS orbital for given iAngmom
!> \param startB the local starting orbital index of the second LHS orbital for given iAngmom
!> \param nAng the number of angular momentum on P
!> \param nCont the number of contracted functions on P
!> \param nPassP the number of passes for LHS overlap P
!> \param ndim5P the number of derivatives for LHS overlap P
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE AddToPQ1(CDAB,AddPQ,nA,nB,nCompA,nCompB,nContA,nContB,&
     &              startA,startB,nAng,nCont,nPassP,ndim5P,permute,startDer,nDer,LUPRI,IPRINT)
implicit none
Integer,intent(IN)        :: nA,nB,nCompA,nCompB,nContA,nContB,startA,startB
Integer,intent(IN)        :: nAng,nCont,nPassP,ndim5P,startDer,nDer,LUPRI,IPRINT
Logical,intent(IN)        :: permute
Real(realk),intent(INOUT) :: CDAB(nA,nB,nDer,nPassP)
Real(realk),intent(IN)    :: AddPQ(nCont,ndim5P,nAng)
!
Integer :: iAng,iContP,iDerivP,iPassP,iContA,iContB,iA,iB,iCompA,iCompB
!
DO iDerivP=1,ndim5P
  iContP=0
  DO iPassP=1,nPassP
    DO iContB=1,nContB
      DO iContA=1,nContA
        iContP = iContP + 1
        iAng=0
        iB = startB + (iContB-1)*nCompB - 1
        DO iCompB=1,nCompB
          iB=iB+1
          iA = startA + (iContA-1)*nCompA - 1
          DO iCompA=1,nCompA
            CDAB(iA+iCompA,iB,startDer+iDerivP,iPassP) = AddPQ(iContP,iDerivP,iAng+iCompA)
          ENDDO !iContA
          IF (permute) THEN
            DO iCompA=1,nCompA
              CDAB(iB,iA+iCompA,startDer+iDerivP,iPassP) = AddPQ(iContP,iDerivP,iAng+iCompA)
            ENDDO !iContA
          ENDIF !permute
          iA   = iA  +  nCompA
          iAng = iAng + nCompA
        ENDDO !iContB
      ENDDO !iPassP
    ENDDO !iDerivP
  ENDDO !iCompA
ENDDO !iCompB

END SUBROUTINE AddToPQ1

!> \brief add the calculated integral into integral%PQ for storage until distribution
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param CDAB  the integrals to be returned
!> \param addPQ the integrals in a form suitable for the inner contraction
!> \param nQ the number of AO functions on the RHS overlap Q
!> \param nA the number of AO functions of the first LHS orbital
!> \param nB the number of AO functions of the second LHS orbital
!> \param nCompA the number of angular components of the first LHS orbital for given iAngmom
!> \param nCompB the number of angular components of the second LHS orbital for given iAngmom
!> \param nContA the number of contracted functions of the first LHS orbital for given iAngmom
!> \param nContB the number of contracted functions of the second LHS orbital for given iAngmom
!> \param startA the local starting orbital index of the first LHS orbital for given iAngmom
!> \param startB the local starting orbital index of the second LHS orbital for given iAngmom
!> \param nAng the number of angular momentum on P
!> \param nCont the number of contracted functions on P
!> \param nPassP the number of passes for LHS overlap P
!> \param ndim5P the number of derivatives for LHS overlap P
!> \param startDer the derivative offset
!> \param nDer the total number of derivative components
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE AddToPQgen(CDAB,AddPQ,nQ,nPassQ,ndim5Q,nA,nB,nCompA,nCompB,nContA,nContB,&
     &              startA,startB,nAng,nCont,nPassP,ndim5P,permute,startDer,nDer,LUPRI,IPRINT)
implicit none
Integer,intent(IN)        :: nQ,nPassQ,ndim5Q,nA,nB,nCompA,nCompB,nContA,nContB,startA,startB
Integer,intent(IN)        :: nAng,nCont,nPassP,ndim5P,startDer,nDer,LUPRI,IPRINT
Logical,intent(IN)        :: permute
Real(realk),intent(INOUT) :: CDAB(nQ,nA,nB,nDer,nPassQ,nPassP)
Real(realk),intent(IN)    :: AddPQ(nCont,nQ,nPassQ*ndim5Q,ndim5P,nAng)
!
Integer :: iQ,iAng,iContP,iDerivP,iPassP,iContA,iContB,iA,iB,iCompA,iCompB
Integer :: iDerivQ,iPassQ,iQpd,iDeriv
!
DO iDerivP=1,ndim5P
  iContP=0
  DO iPassP=1,nPassP
    DO iContB=1,nContB
      DO iContA=1,nContA
        iContP = iContP + 1
        iAng=0
        iB = startB + (iContB-1)*nCompB - 1
        DO iCompB=1,nCompB
          iB=iB+1
          iA = startA + (iContA-1)*nCompA - 1
          DO iCompA=1,nCompA
            iA=iA+1
            iAng = iAng+1
            DO iPassQ=1,nPassQ
             iDeriv = startDer + ndim5Q*(iDerivP-1)
             DO iDerivQ=1,ndim5Q
              iDeriv = iDeriv + 1
              iQpd = iPassQ + (iDerivQ-1)*nPassQ
              DO iQ=1,nQ
               CDAB(iQ,iA,iB,iDeriv,iPassQ,iPassP) = AddPQ(iContP,iQ,iQpd,iDerivP,iAng)
              ENDDO !iQ
              IF (permute) THEN
                DO iQ=1,nQ
                 CDAB(iQ,iB,iA,iDeriv,iPassQ,iPassP) = AddPQ(iContP,iQ,iQpd,iDerivP,iAng)
                ENDDO !iQ
              ENDIF
             ENDDO !iDerivQ
            ENDDO !iPassQ
          ENDDO !iCompA
        ENDDO !iCompB
      ENDDO !iContA
    ENDDO !iContB
  ENDDO !iPassP
ENDDO !iDerivP

END SUBROUTINE AddToPQgen

!> \brief add the calculated integral into integral%PQ for storage until distribution
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param CDAB  the integrals to be returned
!> \param addPQ the integrals in a form suitable for the inner contraction
!> \param nQ the number of AO functions on the RHS overlap Q
!> \param nA the number of AO functions of the first LHS orbital
!> \param nB the number of AO functions of the second LHS orbital
!> \param nCompA the number of angular components of the first LHS orbital for given iAngmom
!> \param nCompB the number of angular components of the second LHS orbital for given iAngmom
!> \param nContA the number of contracted functions of the first LHS orbital for given iAngmom
!> \param nContB the number of contracted functions of the second LHS orbital for given iAngmom
!> \param startA the local starting orbital index of the first LHS orbital for given iAngmom
!> \param startB the local starting orbital index of the second LHS orbital for given iAngmom
!> \param nAng the number of angular momentum on P
!> \param nCont the number of contracted functions on P
!> \param nPassP the number of passes for LHS overlap P
!> \param ndim5P the number of derivatives for LHS overlap P
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE AddToPQgen_ca(CDAB,AddPQ,nQ,nPassQ,ndim5Q,nA,nB,nCompA,nCompB,nContA,nContB,&
     &                   startA,startB,nAng,nCont,nPassP,ndim5P,permute,startDer,nDer,LUPRI,IPRINT)
implicit none
Integer,intent(IN)        :: nQ,nPassQ,ndim5Q,nA,nB,nCompA,nCompB,nContA,nContB,startA,startB
Integer,intent(IN)        :: nAng,nCont,nPassP,ndim5P,startDer,nDer,LUPRI,IPRINT
Logical,intent(IN)        :: permute
Real(realk),intent(INOUT) :: CDAB(nQ,nA,nB,nDer,nPassQ,nPassP)
Real(realk),intent(IN)    :: AddPQ(nCont,nQ,nPassQ*ndim5Q,ndim5P,nAng)
!
Integer :: iQ,iAng,iContP,iDerivP,iPassP,iContA,iContB,iA,iB,iCompA,iCompB
Integer :: iDerivQ,iPassQ,iQpd,iDeriv
!
DO iDerivP=1,ndim5P
  iAng=0
  DO iCompB=1,nCompB
    DO iCompA=1,nCompA
      iAng = iAng+1
      iContP=0
      DO iPassP=1,nPassP
        iB = startB + (iCompB-1)*nContB - 1
        DO iContB=1,nContB
          iB=iB+1
          iA = startA + (iCompA-1)*nContA - 1
          DO iContA=1,nContA
            iA=iA+1
            iContP = iContP + 1
            DO iPassQ=1,nPassQ
             iDeriv = startDer + ndim5Q*(iDerivP-1)
             DO iDerivQ=1,ndim5Q
              iDeriv = iDeriv + 1
              iQpd = iPassQ + (iDerivQ-1)*nPassQ
              DO iQ=1,nQ
               CDAB(iQ,iA,iB,iDeriv,iPassQ,iPassP) = AddPQ(iContP,iQ,iQpd,iDerivP,iAng)
              ENDDO !iQ
              IF (permute) THEN
                DO iQ=1,nQ
                 CDAB(iQ,iB,iA,iDeriv,iPassQ,iPassP) = AddPQ(iContP,iQ,iQpd,iDerivP,iAng)
                ENDDO !iQ
              ENDIF
             ENDDO !iDerivQ
            ENDDO !iPassQ
          ENDDO !iContA
        ENDDO !iContB
      ENDDO !iPassP
    ENDDO !iCompA
  ENDDO !iCompB
ENDDO !iDerivP

END SUBROUTINE AddToPQgen_ca

SUBROUTINE PrintPQ(CDAB,nA,nB,nC,nD,ndim5P,ndim5Q,nPassP,nPassQ,nDer,LUPRI,IPRINT)
implicit none
Integer,intent(IN)     :: nA,nB,nC,nD,ndim5P,ndim5Q,nPassP,nPassQ,nDer,LUPRI,IPRINT
Real(realk),intent(IN) :: CDAB(nC,nD,nA,nB,nDer,nPassQ*nPassP)
!
Integer :: iA,iB,iC,iD,iDerivQ,iPassQ,iDerivP,iPassP,iDeriv,iPass

CALL LSHEADER(LUPRI,'PQ contribution')
iPass = 0
DO iPassP=1,nPassP
 DO iPassQ=1,nPassQ
  iPass=iPass+1
  WRITE(LUPRI,'(3X,A,I2,A,I2,A)') 'Pass numbers (iPassP,iPassQ)=(',iPassP,',',iPassQ,')'
  DO iDeriv=1,nDer
    WRITE(LUPRI,'(3X,A,I2,A)') 'Derivative components (iDeriv)=(',iDeriv,')'
    DO iB=1,nB
     DO iA=1,nA
      DO iD=1,nD
       WRITE(LUPRI,'(5X,A,I2,A,I2,A,I2)') 'iD =',iD,' iA =',iA,' iB =',iB
       WRITE(LUPRI,'(7X,5F28.17)') (CDAB(iC,iD,iA,iB,iDeriv,iPass),iC=1,nC)
      ENDDO !iC
     ENDDO !iA
    ENDDO !iB
  ENDDO !iDeriv
 ENDDO !iPassQ
ENDDO !iPassP
END SUBROUTINE PrintPQ

!> \brief wrapper to distribute primitive hermite integrals to integral%IN 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> distribute primitive hermite integrals to integral%IN 
!> after the contraction of the RHS overlap Q have been performed 
!>
!> \param integral the storage of integrals and intermediates
!> \param PQ contains info about the integrand, to be built
!> \param iAngmom the index of angular momentum
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE DistributeHermiteP(Integral,PQ,iAngmom,LUPRI,IPRINT,input)
implicit none
TYPE(Integralitem)   :: Integral
TYPE(Integrand)      :: PQ
TYPE(integralInput) :: input
Integer             :: LUPRI,IPRINT,iAngmom
!
Integer :: ntuvP,ntuvPfull,nP,nOrbQ,startOrbQ,endOrbQ
Integer :: startP,fullSP,endP,ioffP,fullOP,jP,tP,uP,vP,ituvP,ifullP,iOrbQ,iPrimP
Integer :: ifullp1,ifullp2,ifullp3,end2P,der
!
der = integral%lhsGeoOrder + PQ%P%p%magderiv

startP = 0
IF (PQ%P%p%type_hermite_single) startP = PQ%P%p%angmom(iAngmom) + der

fullSP   = PQ%P%p%startAngmom
endP     = PQ%P%p%angmom(iAngmom) + der
nP       = PQ%P%p%nPrimitives
ioffP    = startP*(startP+1)*(startP+2)/6
fullOP   = fullSP*(fullSP+1)*(fullSP+2)/6
ntuvP    = (endP+1)*(endP+2)*(endP+3)/6 - ioffP
ntuvPfull = PQ%P%p%nTUV

IF (integral%dohodi) THEN !Special case for HODI
  ntuvPfull = ntuvP 
  fullOP    = ioffP
ENDIF

nOrbQ = PQ%Q%p%orbital1%totOrbitals*PQ%Q%p%orbital2%totOrbitals
IF (.NOT.PQ%Q%p%type_ftuv) nOrbQ = nOrbQ*PQ%Q%p%nPasses*integral%rhsGeoComp*PQ%Q%p%nCartesianMomentComp
IF (INPUT%DO_MULMOM) nOrbQ = Input%nMultipoleMomentComp

Integral%nAng  = ntuvP
Integral%nPrim = np
Integral%nOrb  = nOrbQ
Integral%ngeoDeriv  = 1

#ifdef VAR_LSDEBUGINT
IF(nP*nOrbQ*ntuvP.GT.allocIntmaxTUVdim)THEN
   call lsquit('DistributeHermiteP_regular alloc error3',lupri)
ENDIF
IF(nP*PQ%P%p%nTUV*nOrbQ.GT.allocIntmaxTUVdim)THEN
   call lsquit('DistributeHermiteP_regular alloc error4',lupri)
ENDIF
#endif
IF(nP.GT. 2)THEN !General case
   CALL DistributeHermiteP_regularGen(Integral%IN,Integral%tuvQ,Integral%TUV%TUVindex,&
        &                             startP,endP,nP,ioffP,fullOP,ntuvPfull,ntuvP,nOrbQ,LUPRI,IPRINT)
ELSEIF(nP.EQ. 1)THEN !Special case for nPrim=1
   CALL DistributeHermiteP_regular1(Integral%IN,Integral%tuvQ,Integral%TUV%TUVindex,&
        &                             startP,endP,nP,ioffP,fullOP,ntuvPfull,ntuvP,nOrbQ,LUPRI,IPRINT)

ELSE !Special case for nPrim=2
   CALL DistributeHermiteP_regular2(Integral%IN,Integral%tuvQ,Integral%TUV%TUVindex,&
        &                             startP,endP,nP,ioffP,fullOP,ntuvPfull,ntuvP,nOrbQ,LUPRI,IPRINT)

ENDIF

IF (IPRINT.GT. 20) THEN
   CALL PrintTensor(Integral%IN,'DistributeHermiteP  ',&
        &nP,nOrbQ,ntuvP,Lupri,'Prim  ','iOrb  ','iTUV  ',3)
ENDIF

END SUBROUTINE DistributeHermiteP

!> \brief distribute primitive hermite integrals to integral%IN, the default worker routine 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> distribute primitive hermite integrals to integral%IN 
!> after the contraction of the RHS overlap Q have been performed 
!>
!> \param integralIN the storage of the primitive hermite integrals
!> \param TUVQ the intermidiate integrals (at this point the RHS contraction have been performed)
!> \param TUVindex for a given t,u,v values it gives the TUV index (t,u,v)=(000,100,010,001,.. etc) means TUV = (1,2,3,4,..)
!> \param startp start angular momentum  on overlap P
!> \param endp end angular momentum on overlap P 
!> \param nP number of primitve(*nPasses) for the P overlap
!> \param ioffP offset for the overlap P tuv index 
!> \param fullOP offset for the full tuv index 
!> \param ntuvfull number of tuv indexes for the full P and Q overlap
!> \param ntuvP  number of tuv indexes for the P overlap 
!> \param nOrbQ the total number of orbitals
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE DistributeHermiteP_regularGen(IntegralIN,TUVQ,TUVindex,startP,endP,nP,ioffP,&
     &                                fullOP,ntuvFull,ntuvP,nOrbQ,LUPRI,IPRINT)
implicit none
Integer,intent(in) :: startP,endP,nP,ioffP,fullOP,ntuvFull,ntuvP,nOrbQ,LUPRI,IPRINT
Real(realk),intent(inout) :: IntegralIN(nP,nOrbQ,ntuvP)
Real(realk),intent(in)    :: TUVQ(np,ntuvFull,nOrbQ)
Integer,pointer :: TUVindex(:,:,:)
!
Integer :: jP,tP,uP,vP,ituvP,ifullP,iOrbQ,iPrimP
!
DO jP = startP,endP
   DO tP=jP,0,-1
      DO uP=jP-tP,0,-1
         vP=jP-tP-uP
         ituvP  = TUVindex(tP,uP,vP)-ioffP
         ifullP = TUVindex(tP,uP,vP)-fullOP
         DO iOrbQ=1,nOrbQ
            DO iPrimP=1,nP
               IntegralIN(iPrimP,iOrbQ,ituvP) = TUVQ(iPrimP,ifullP,iOrbQ)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE DistributeHermiteP_regularGen

!> \brief distribute primitive hermite integrals to integral%IN, the default worker routine 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> distribute primitive hermite integrals to integral%IN 
!> after the contraction of the RHS overlap Q have been performed 
!>
!> \param integralIN the storage of the primitive hermite integrals
!> \param TUVQ the intermidiate integrals (at this point the RHS contraction have been performed)
!> \param TUVindex for a given t,u,v values it gives the TUV index (t,u,v)=(000,100,010,001,.. etc) means TUV = (1,2,3,4,..)
!> \param startp start angular momentum  on overlap P
!> \param endp end angular momentum on overlap P 
!> \param nP number of primitve(*nPasses) for the P overlap
!> \param ioffP offset for the overlap P tuv index 
!> \param fullOP offset for the full tuv index 
!> \param ntuvfull number of tuv indexes for the full P and Q overlap
!> \param ntuvP  number of tuv indexes for the P overlap 
!> \param nOrbQ the total number of orbitals
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE DistributeHermiteP_regular1(IntegralIN,TUVQ,TUVindex,startP,endP,nP,ioffP,&
     &                                fullOP,ntuvFull,ntuvP,nOrbQ,LUPRI,IPRINT)
implicit none
Integer,intent(in) :: startP,endP,nP,ioffP,fullOP,ntuvFull,ntuvP,nOrbQ,LUPRI,IPRINT
Real(realk),intent(inout) :: IntegralIN(nOrbQ,ntuvP)
Real(realk),intent(in)    :: TUVQ(ntuvFull,nOrbQ)
Integer,pointer :: TUVindex(:,:,:)
!
Integer :: jP,tP,uP,vP,ituvP,ifullP,iOrbQ,iPrimP
!
DO jP = startP,endP
   DO tP=jP,0,-1
      DO uP=jP-tP,0,-1
         vP=jP-tP-uP
         ituvP  = TUVindex(tP,uP,vP)-ioffP
         ifullP = TUVindex(tP,uP,vP)-fullOP
         DO iOrbQ=1,nOrbQ
            IntegralIN(iOrbQ,ituvP) = TUVQ(ifullP,iOrbQ)
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE DistributeHermiteP_regular1

!> \brief distribute primitive hermite integrals to integral%IN, the default worker routine 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> distribute primitive hermite integrals to integral%IN 
!> after the contraction of the RHS overlap Q have been performed 
!>
!> \param integralIN the storage of the primitive hermite integrals
!> \param TUVQ the intermidiate integrals (at this point the RHS contraction have been performed)
!> \param TUVindex for a given t,u,v values it gives the TUV index (t,u,v)=(000,100,010,001,.. etc) means TUV = (1,2,3,4,..)
!> \param startp start angular momentum  on overlap P
!> \param endp end angular momentum on overlap P 
!> \param nP number of primitve(*nPasses) for the P overlap
!> \param ioffP offset for the overlap P tuv index 
!> \param fullOP offset for the full tuv index 
!> \param ntuvfull number of tuv indexes for the full P and Q overlap
!> \param ntuvP  number of tuv indexes for the P overlap 
!> \param nOrbQ the total number of orbitals
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE DistributeHermiteP_regular2(IntegralIN,TUVQ,TUVindex,startP,endP,nP,ioffP,&
     &                                fullOP,ntuvFull,ntuvP,nOrbQ,LUPRI,IPRINT)
implicit none
Integer,intent(in) :: startP,endP,nP,ioffP,fullOP,ntuvFull,ntuvP,nOrbQ,LUPRI,IPRINT
Real(realk),intent(inout) :: IntegralIN(nP,nOrbQ,ntuvP)
Real(realk),intent(in)    :: TUVQ(np,ntuvFull,nOrbQ)
Integer,pointer :: TUVindex(:,:,:)
!
Integer :: jP,tP,uP,vP,ituvP,ifullP,iOrbQ,iPrimP
!
DO jP = startP,endP
   DO tP=jP,0,-1
      DO uP=jP-tP,0,-1
         vP=jP-tP-uP
         ituvP  = TUVindex(tP,uP,vP)-ioffP
         ifullP = TUVindex(tP,uP,vP)-fullOP
         DO iOrbQ=1,nOrbQ
            IntegralIN(1,iOrbQ,ituvP) = TUVQ(1,ifullP,iOrbQ)
            IntegralIN(2,iOrbQ,ituvP) = TUVQ(2,ifullP,iOrbQ)
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE DistributeHermiteP_regular2

!> \brief distribute primitive hermite integrals to integral%IN, special case
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> distribute primitive hermite integrals to integral%IN 
!> after the contraction of the RHS overlap Q have been performed 
!> this is a special case for the kinetic energy integral, where we apply 
!> nabla to the right hand side overlap Q, and we therefor need to sum the x,y and z componets
!>
!> \param integralIN the storage of the primitive hermite integrals
!> \param TUVQ the intermidiate integrals (at this point the RHS contraction have been performed)
!> \param TUVindex for a given t,u,v values it gives the TUV index (t,u,v)=(000,100,010,001,.. etc) means TUV = (1,2,3,4,..)
!> \param startp start angular momentum  on overlap P
!> \param endp end angular momentum on overlap P
!> \param nP number of primitve(*nPasses) for the P overlap
!> \param ioffP offset for the overlap P tuv index 
!> \param fullOP offset for the full tuv index 
!> \param ntuvfull number of tuv indexes for the full P and Q overlap
!> \param ntuvP  number of tuv indexes for the P overlap 
!> \param nOrbQ the total number of orbitals
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE WTUV_extract_kinetic(IntegralIN,TUVQ,TUVindex,startP,endP,nP,ioffP,&
     &                          fullOP,ntuvFull,ntuvP,nOrbQ,deriv,LUPRI,IPRINT)
implicit none
Integer,intent(in) :: startP,endP,nP,ioffP,fullOP,ntuvFull,ntuvP,nOrbQ,LUPRI,IPRINT
Real(realk),intent(inout) :: IntegralIN(nP,nOrbQ,ntuvP)
Real(realk),intent(in)    :: TUVQ(np,ntuvFull,nOrbQ)
Integer,pointer    :: TUVindex(:,:,:)
logical,intent(in) :: deriv
!
Integer :: jP,tP,uP,vP,ituvP,iOrbQ,iPrimP
Integer :: ifullp1,ifullp2,ifullp3
Real(realk), parameter :: Half=0.5E0_realk,One=1E0_realk
Real(realk) :: factor

IF (deriv) THEN
  factor = One
ELSE
  factor = Half
ENDIF

DO jP = startP,endP
   DO tP=jP,0,-1
      DO uP=jP-tP,0,-1
         vP=jP-tP-uP
         ituvP  = TUVindex(tP,uP,vP)-ioffP
         ifullP1 = TUVindex(tP+2,uP,vP)-fullOP
         ifullP2 = TUVindex(tP,uP+2,vP)-fullOP
         ifullP3 = TUVindex(tP,uP,vP+2)-fullOP             
         DO iOrbQ=1,nOrbQ
            DO iPrimP=1,nP
               IntegralIN(iPrimP,iOrbQ,ituvP) =  factor*(&
                    &- TUVQ(iPrimP,ifullP1,iOrbQ) &
                    &- TUVQ(iPrimP,ifullP2,iOrbQ) &
                    &- TUVQ(iPrimP,ifullP3,iOrbQ))
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE WTUV_extract_kinetic

!> \brief contract the RHS overlap Q (primitive to contracted functions, Ecoeff contraction,..)
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param INTEGRAL storage of integrals and intermediates
!> \param PQ contains info about the integrand, to be built. 
!> \param INPUT contains info about the requested integral evaluation
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE Contract_Q(INTEGRAL,PQ,Input,LUPRI,IPRINT)
implicit none
TYPE(Integralitem)  :: INTEGRAL
TYPE(Integrand)     :: PQ
TYPE(IntegralInput) :: Input
Integer             :: LUPRI,IPRINT
!
Integer  :: iAngmom,startMomQ,iorder,ndim5Comp,nDerivQ
!
IF (INPUT%nMultipoleMomentComp.GT. 1.OR.INPUT%DO_MULMOM) THEN
   ndim5Comp = INPUT%nMultipoleMomentComp
   nDerivQ = 1
ELSE
   ndim5Comp=Input%nCartesianMomentComp
   nDerivQ = Integral%rhsGeoComp
ENDIF

!Loop for multipole momemnt - should be revised!
DO iOrder=1,ndim5Comp
   !Loop over angular contributions sharing the set of primitive functions
   DO iAngmom=1,PQ%Q%p%nAngmom
      CALL ContractEcoeffQ(Integral,PQ,iAngmom,iOrder,LUPRI,IPRINT)
      !Output(nPrimQ,nPrimP,ntuvP,ijkQ)
      CALL ContractBasis(Integral,PQ%Q%p,iAngmom,LUPRI,IPRINT)
      !Output(nContQ,nPrimP,ntuvP,nDerivQ,ijkQ)
      CALL extractDifferentiated(Integral,PQ%Q%p,iAngmom,Integral%rhsGeoOrder,Integral%rhsGeoComp,LUPRI,IPRINT)
!     IF (PQ%Q%p%magderiv.EQ.1)CALL extractMagDifferentiated(Integral,PQ%Q%p,iAngmom,LUPRI,IPRINT)
      !Output(nContQ,nPrimP,ntuvP,ijkQ)
      CALL SphericalTransform(Integral,PQ%Q%p,iAngmom,LUPRI,IPRINT)
      !Output(nContQ,nPrimP,ntuvP,nDerivQ,lmQ)
      CALL AddToTUVQ(Integral,PQ,iAngmom,iOrder,ndim5Comp,nDerivQ,Input%contAng,LUPRI,IPRINT)
   ENDDO
ENDDO
    
END SUBROUTINE Contract_Q

!> \brief wrapper contract Ecoefficients for the right hand side overlap Q
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param INTEGRAL storage of integrals and intermediates
!> \param PQ contains info about the integrand, to be built.
!> \param iAngmom the index of angular momentum
!> \param ideriv the derivative index
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE ContractEcoeffQ(Integral,PQ,iAngmom,iDeriv,LUPRI,IPRINT)
implicit none
TYPE(Integralitem)  :: INTEGRAL
TYPE(Integrand)     :: PQ
Integer             :: LUPRI,IPRINT
Integer             :: iAngmom,iDeriv
!
TYPE(Overlap),pointer :: P,Q
Integer               :: jP,jQ,tP,uP,vP,tQ,uQ,vQ,ituvP,ituvQ,ituvPQ,nPrimQ,nPassQ
Integer               :: nP,nQ,ioffP,ioffQ,sPQ,ioffPQ,i1,i2,j1,j2,nprim1,nprim2
Integer               :: iPrimPQ,iP,iQ,iTUV,maxQ,minQ,nEFG,l1,l2,ijk1,ijk2,ijk,ijkcart
Integer               :: ntuvP,ntuvQ,startP,endP,startQ,endQ,endPQ,ntuvPQ,startE
Integer               :: start1,start2,nAng,ipassQ,iangQ,iprimq,derP,derQ
Real(realk),parameter :: signQ = -1E0_realk
Real(realk)           :: DM1 = -1.0E0_realk,sign
Real(realk),pointer   :: Ecoeffs(:)

real(4) :: tarr(2),etime,tnew,told,tstart
Real(realk),dimension(:),pointer :: ptemp
Real(realk),parameter :: D2=2.0E0_realk,D1=1.0E0_realk

P   => PQ%P%p
Q   => PQ%Q%p
derP = Integral%lhsGeoOrder
derQ = Integral%rhsGeoOrder

startP = P%startAngmom
IF (integral%dohodi) THEN
  endP   = P%maxAngmom + derP
  IF (P%type_hermite_single) startP = endP
ELSE
  endP = P%endAngmom
ENDIF

startQ = 0
maxQ   = Q%endAngmom
minQ   = Q%startAngmom
endQ   = Q%angmom(iAngmom) + derQ
IF (Q%type_hermite_single) startQ = endQ
nP     = P%nPrimitives*P%nPasses
nQ     = Q%nPrimitives*Q%nPasses
ioffP  = startP*(startP+1)*(startP+2)/6
ioffQ  = startQ*(startQ+1)*(startQ+2)/6
!sPQ    = startP+Q%startAngmom
sPQ    = PQ%startAngmom
ioffPQ = sPQ*(sPQ+1)*(sPQ+2)/6
endPQ  = maxQ + P%endAngmom
ntuvP  = (endP+1)*(endP+2)*(endP+3)/6-ioffP
ntuvQ  = (endQ+1)*(endQ+2)*(endQ+3)/6-ioffQ
ntuvPQ = (endPQ+1)*(endPQ+2)*(endPQ+3)/6-ioffPQ
Integral%nOrb  = ntuvP*nP
Integral%nAng  = ntuvQ
Integral%nPrim = nQ
Integral%ngeoDeriv = 1
nEFG=Integral%nEFG
!IF ((((maxQ.EQ. 0).OR.((endP.EQ. 0).AND.((maxQ.EQ.endQ).AND.(minQ.EQ.startQ)))).AND.(nEFG.EQ. 1)).AND.Q%type_hermite_Single) THEN
IF (.FALSE.) THEN
   CALL SingleHermiteEcoeff(Integral,Q,signQ,iAngmom,Q%nPasses,ntuvQ,ntuvP*nP,LUPRI,IPRINT)
ELSE
   IF (Q%type_hermite_single) THEN
      CALL DirectSingleHermiteEcoeff(Integral,Q,nP,nQ,startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ,ntuvPQ*nEFG,&
             &                       ntuvPQ,ntuvP,ntuvQ,ideriv,signQ,derQ,iAngmom,lupri)
   ELSE
      l1 = Q%orbital1%angmom(Q%indexAng1(iangmom))
      l2 = Q%orbital2%angmom(Q%indexAng2(iangmom))
      CALL GET_IJK(l1,l2,ijk1,ijk2,ijk,ijkcart,Q%sphericalEcoeff,derQ,Q%single)
!      print*,'debug 53 Ecoeffs',Integral%nPrim*Integral%nAng*ijk
      call mem_workpointer_alloc(Ecoeffs,Integral%nPrim*Integral%nAng*ijk)
      CALL BuildEcoeffTensor(integral%TUV,Q,signQ,Ecoeffs,ijk,ijkcart,&
           & Integral%nAng,Integral%nPrim,iAngmom,derQ,Q%nPasses,LUPRI,IPRINT,1)
      IF (Q%segmented) THEN         
#ifdef VAR_LSDEBUGINT
         IF(Q%nPrimitives*Q%nPasses*nP*ntuvPQ*nEFG.GT.allocIntmaxTUVdim)THEN
            call lsquit('DirectcontractEQseg alloc error1',lupri)
         ENDIF
         IF(Q%nPasses*nP*ntuvP*ijk.GT.allocIntmaxTUVdim)THEN
            print*,'Q%nPasses*nP*ntuvP*ijk',Q%nPasses*nP*ntuvP*ijk
            print*,'allocIntmaxTUVdim',allocIntmaxTUVdim
            print*,'Q%nPasses',Q%nPasses
            print*,'nP',nP
            print*,'ntuvP',ntuvP
            print*,'ijk',ijk
            CALL PRINT_OVERLAP(PQ%P%p,6,1000,'LHS')
            CALL PRINT_OVERLAP(PQ%Q%p,6,1000,'RHS') 
            call lsquit('DirectcontractEQseg alloc error2',lupri)
         ENDIF
#endif
         CALL DirectcontractEQseg(Integral%IN,Integral%RTUV,Integral%TUV%TUVindex,nP*nQ,&
              & nP,Q%nPrimitives,Q%nPasses,startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ,ntuvPQ*nEFG,&
              & ntuvPQ,ntuvP,ntuvQ,ideriv,Ecoeffs,1,ijk,l1,l2,Integral%TUV%SYMindex,Integral%TUV%iPQxyz,lupri)
         Integral%nPrim = Q%nPasses
      ELSE
#ifdef VAR_LSDEBUGINT
         IF(nQ*nP*ntuvPQ*nEFG.GT.allocIntmaxTUVdim)THEN
            call lsquit('DirectcontractEQgen alloc error1',lupri)
         ENDIF
         IF(nQ*nP*ntuvP*ijk.GT.allocIntmaxTUVdim)THEN
            print*,'nQ*nP*ntuvP*ijk',nQ*nP*ntuvP*ijk
            print*,'allocIntmaxTUVdim',allocIntmaxTUVdim
            print*,'nQ',nQ,'Q%nPasses',Q%nPasses,'Q%nPrimitives',Q%nPrimitives
            print*,'nP',nP
            print*,'ntuvP',ntuvP
            print*,'ijk',ijk
            CALL PRINT_OVERLAP(PQ%P%p,6,1000,'LHS')
            CALL PRINT_OVERLAP(PQ%Q%p,6,1000,'RHS') 
            call lsquit('DirectcontractEQgen alloc error2',lupri)
         ENDIF
#endif
         CALL DirectcontractEQgen(Integral%IN,Integral%RTUV,Integral%TUV%TUVindex,nP*nQ,&
              & nP,nQ,startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ,ntuvPQ*nEFG,ntuvPQ,ntuvP,&
              & ntuvQ,ideriv,Ecoeffs,1,ijk,l1,l2,Integral%TUV%SYMindex,Integral%TUV%iPQxyz,lupri)
      ENDIF
      call mem_workpointer_dealloc(Ecoeffs)
      Integral%nAng = ijk
   ENDIF
ENDIF
IF(IPRINT.GE. 100)THEN
  CALL Print_DirectContractOutput(Integral%IN,Integral%nPrim,Integral%nOrb,Integral%nAng,lupri)
ENDIF
END SUBROUTINE ContractEcoeffQ

subroutine Print_DirectContractOutput(WTUVEcoeff,nPrim,nOrb,nAng,lupri)
implicit none
integer :: nPrim,nOrb,nAng,lupri
real(realk) :: WTUVEcoeff(nPrim,nOrb,nAng)
!
integer :: iOrbP,iAngQ,iPrimQ

WRITE(LUPRI,*)'PRINT_DIRECTCONTRACTOUTPUT'
WRITE(LUPRI,'(1X,A6,1X,A6,1X,A11,I5 )') 'iOrbP','iAngQ','iPrimQ = 1,',nPrim
DO iAngQ = 1,nAng
   DO iOrbP = 1,nOrb
      WRITE(LUPRI,'(1X,I4,2X,I4,2X,2ES16.6/,(13X,2ES16.6))') iOrbP,iAngQ,(WTUVEcoeff(iPrimQ,iOrbP,iAngQ),iPrimQ=1,nPrim)
   ENDDO
ENDDO
END subroutine PRINT_DIRECTCONTRACTOUTPUT

!> \brief contraction of single hermite Ecoefficients for the right hand side overlap Q
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param Herint storage of hermite integrals
!> \param P contains overlap distribution for the left or right hand side
!> \param signP the sign (+1 or -1)
!> \param j1 angular momentum of orbital 1
!> \param j2 angular momentum of orbital 2
!> \param nPrim1 number of primitive for orbital 1
!> \param nPrim2 number of primitive for orbital 2
!> \param nPasses number of passes 
!> \param nAng number of angular components 
!> \param nOrb number of orbitals
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SingleHermiteEcoeff1gen(Herint,P,signP,j1,j2,nPrim1,nPrim2,nPasses,nAng,nOrb,preexpfac,&
     & nPrimAlloc,nAngAlloc,iAngmomQ,LUPRI,IPRINT)
implicit none
Integer       :: nAng,nOrb,LUPRI,IPRINT,nPrimAlloc,nAngAlloc,iAngmomQ
Integer       :: nPrim1,nPrim2,j1,j2,nPasses
Real(realk)   :: HerInt(nPrim1,nPrim2,nPasses,nAng*nOrb)
TYPE(Overlap) :: P
Real(realk)   :: signP,preexpfac(nPrimAlloc,nAngAlloc)
!
Real(realk)   :: pref2,pref12
Real(realk),parameter   :: D1=1.0E0_realk,D2=2.0E0_realk
Integer       :: iPrim1,iPrim2,iAngOrb,iPasses,iPrimQ

IF(P%contractbasis)THEN
   !outside or inside?  
   DO iAngOrb=1,nAng*nOrb
      DO iPasses=1,nPasses
         DO iPrim2=1,nPrim2
            pref2=signP/(D2*P%orbital2%exponents(iPrim2))**j2
            DO iPrim1=1,nPrim1
               pref12=D1/(D2*P%orbital1%exponents(iPrim1))**j1*pref2
               HerInt(iPrim1,iPrim2,iPasses,iAngOrb) = HerInt(iPrim1,iPrim2,iPasses,iAngOrb) * pref12
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ELSE
!PS_int
   DO iAngOrb=1,nAng*nOrb
      iPrimQ=0
      DO iPasses=1,nPasses
         DO iPrim2=1,nPrim2
            pref2=signP/(D2*P%orbital2%exponents(iPrim2))**j2
            DO iPrim1=1,nPrim1
               iPrimQ=iPrimQ+1
               pref12=D1/(D2*P%orbital1%exponents(iPrim1))**j1*pref2*preexpfac(iPrimQ,iAngmomQ)
               HerInt(iPrim1,iPrim2,iPasses,iAngOrb) = HerInt(iPrim1,iPrim2,iPasses,iAngOrb) * pref12
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE SingleHermiteEcoeff1gen

!> \brief contraction of single hermite Ecoefficients for the right hand side overlap Q
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param Herint storage of hermite integrals
!> \param P contains overlap distribution for the left or right hand side
!> \param signP the sign (+1 or -1)
!> \param j1 angular momentum of orbital 1
!> \param j2 angular momentum of orbital 2
!> \param nPrim1 number of primitive for orbital 1
!> \param nPrim2 number of primitive for orbital 2
!> \param nPasses number of passes 
!> \param nAng number of angular components 
!> \param nOrb number of orbitals
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SingleHermiteEcoeff1seg(ContInt,Herint,P,signP,j1,j2,nPrim1,nPrim2,nPasses,nAng,nOrb,preexpfac,&
     & nPrimAlloc,nAngAlloc,iAngmom,LUPRI,IPRINT)
implicit none
Integer       :: nAng,nOrb,LUPRI,IPRINT,nPrimAlloc,nAngAlloc,iAngmom
Integer       :: nPrim1,nPrim2,j1,j2,nPasses
Real(realk)   :: ContInt(nPasses,nAng*nOrb)
Real(realk)   :: HerInt(nPrim1,nPrim2,nPasses,nAng*nOrb)
TYPE(Overlap) :: P
Real(realk)   :: signP,preexpfac(nPrimAlloc,nAngAlloc)
!
Real(realk)   :: pref2,pref12,tmp
Real(realk),parameter   :: D1=1.0E0_realk,D2=2.0E0_realk
Integer       :: iPrim1,iPrim2,iAngOrb,iPasses,iPrimQ

DO iAngOrb=1,nAng*nOrb
   iPrimQ=0
   DO iPasses=1,nPasses
      tmp = 0E0_realk
      DO iPrim2=1,nPrim2
         pref2=signP/(D2*P%orbital2%exponents(iPrim2))**j2
         DO iPrim1=1,nPrim1
            iPrimQ=iPrimQ+1
            pref12=D1/(D2*P%orbital1%exponents(iPrim1))**j1*pref2*preexpfac(iPrimQ,iAngmom)
            tmp = tmp + HerInt(iPrim1,iPrim2,iPasses,iAngOrb)*pref12
         ENDDO
      ENDDO
      ContInt(iPasses,iAngOrb) = tmp
   ENDDO
ENDDO

END SUBROUTINE SingleHermiteEcoeff1seg

SUBROUTINE DirectSingleHermiteEcoeff(Integral,Q,nP,nQ,startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ,ntuvPQnEFG,&
             &                       ntuvPQ,ntuvP,ntuvQ,ideriv,signQ,derQ,iAngmom,lupri)
implicit none
Type(integralItem),intent(inout) :: Integral
Type(overlap),intent(in)         :: Q
Integer,intent(IN)               :: nP,nQ,startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ,ntuvPQnEFG
Integer,intent(IN)               :: ntuvPQ,ntuvP,ntuvQ,ideriv,derQ,iAngmom,lupri
Real(realk),intent(IN)           :: signQ

integer     :: l1,l2,nPrim1,nPrim2
Real(realk) :: sign
Real(realk),parameter :: D1 = 1E0_realk

l1 = Q%orbital1%angmom(Q%indexAng1(iangmom))
l2 = Q%orbital2%angmom(Q%indexAng2(iangmom))
IF(l1.EQ. 0 .AND. l2 .EQ. 0 .AND. derQ.EQ.0)THEN
   sign = D1
ELSE
   sign = signQ**(l1+l2+derQ)
ENDIF
nPrim1 = Q%orbital1%nPrimitives
nPrim2 = Q%orbital2%nPrimitives
IF (Q%segmented) THEN
#ifdef VAR_LSDEBUGINT
IF(Q%nPrimitives*nP*Q%nPasses*ntuvPQnEFG.GT.allocIntmaxTUVdim)THEN
   call lsquit('DirectSingleHermiteEcoeffseg alloc error1',lupri)
ENDIF
IF(nP*Q%nPasses*ntuvP*ntuvQ.GT.allocIntmaxTUVdim)THEN
   call lsquit('DirectSingleHermiteEcoeffseg alloc error2',lupri)
ENDIF
IF(Q%nPasses*Q%nPrimitives.GT.Q%nPrimAlloc)THEN
   call lsquit('DirectSingleHermiteEcoeffseg alloc error3',lupri)
ENDIF

#endif
  CALL DirectSingleHermiteEcoeffseg(Integral%IN,Integral%RTUV,Integral%TUV%TUVindex,nP*nQ,&
       &                  nP,Q%nPrimitives,startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ,ntuvPQnEFG,&
       &                  ntuvPQ,ntuvP,ntuvQ,ideriv,Q,sign,l1,l2,nPrim1,nPrim2,&
       &                  Q%nPasses,Q%preexpfac,Q%nPrimAlloc,Q%nAngAlloc,iAngmom,lupri)
  Integral%nPrim = Q%nPasses
ELSE
#ifdef VAR_LSDEBUGINT
IF(nP*nQ*ntuvPQnEFG.GT.allocIntmaxTUVdim)THEN
   call lsquit('DirectSingleHermiteEcoeffgen alloc error1',lupri)
ENDIF
IF(nQ*nP*ntuvP*ntuvQ.GT.allocIntmaxTUVdim)THEN
   call lsquit('DirectSingleHermiteEcoeffgen alloc error2',lupri)
ENDIF
#endif
  CALL DirectSingleHermiteEcoeffgen(Integral%IN,Integral%RTUV,Integral%TUV%TUVindex,nP*nQ,&
       &                  nP,nQ,startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ,ntuvPQnEFG,&
       &                  ntuvPQ,ntuvP,ntuvQ,ideriv,Q,sign,l1,l2,nPrim1,nPrim2,&
       &                  Q%nPasses,Q%preexpfac,Q%nPrimAlloc,Q%nAngAlloc,iAngmom,lupri)
ENDIF
END SUBROUTINE DirectSingleHermiteEcoeff

!> \brief direct contraction of single hermite Ecoefficients for the right hand side overlap Q
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param OUT the output intermidiates
!> \param WTUV2 the primitive hermite integrals 
!> \param nPrim total number of primitive
!> \param nPrimP number of primitive for overlap P
!> \param nPrimQ number of primitive for overlap Q
!> \param startP starting angular momentum for overlap P
!> \param endP end angular momentum for overlap P
!> \param startQ starting angular momentum for overlap Q
!> \param endQ end angular momentum for overlap Q
!> \param ioffP offset for tuvindexing for overlap P
!> \param ioffQ offset for tuvindexing for overlap Q
!> \param ioffPQ offset for tuvindexing for overlap PQ
!> \param nTUVEFGPQ total number of angular components
!> \param ntuvPQ number of angular components for PQ
!> \param ntuvP number of angular components for P
!> \param ntuvQ number of angular components for Q
!> \param ideriv derivative index
!> \param Q the overlap distribution
!> \param signQ the sign for the Q overlap
!> \param l1 angular momentum for orbital 1
!> \param l2 angular momentum for orbital 2
!> \param nPrim1 number of primitives for orbital 1
!> \param nPrim2 number of primitives for orbital 2
!> \param nPasses number of passes
!> \param LUPRI the logical unit number for the output file
SUBROUTINE DirectSingleHermiteEcoeffgen(OUT,WTUV2,TUVindex,nPrim,nPrimP,nPrimQ,startP,endP,startQ,endQ,&
     & ioffP,ioffQ,ioffPQ,nTUVEFGPQ,ntuvPQ,ntuvP,ntuvQ,ideriv,Q,signQ,l1,l2,nPrim1,nPrim2,&
     & nPasses,preexpfac,nPrimAlloc,nAngAlloc,iAngmomQ,lupri)
implicit none
TYPE(Overlap)      :: Q
Integer         :: nPrim,startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ
Integer         :: ntuvEFGPQ,ntuvPQ,ntuvP,ntuvQ,ideriv,nPrimP,nPrimQ,STARTE
Integer         :: j1,j2,nprim1,nprim2,npasses,ipasses,iprim1,iprim2
Integer         :: nPrimAlloc,nAngAlloc,iAngmomQ
!Real(realk)     :: WTUV(ntuvEFGPQ,nPrim),tuvTUV(ntuvQ,ntuvP,nPrim)
Real(realk)     :: WTUV2(nPrim,ntuvEFGPQ)!,tuvTUV(ntuvQ,ntuvP,nPrim)
!Real(realk)     :: OUT(ntuvQ,ntuvP,nPrim),signQ,pref2,pref12
Real(realk)     :: OUT(nPrimQ,nPrimP,ntuvP,ntuvQ),signQ,pref2,pref12
Real(realk)     :: preexpfac(nPrimAlloc,nAngAlloc)
Integer,pointer :: TUVindex(:,:,:)
!Real(realk)     :: Ec(nPrimQ)
Integer             :: l1,l2,ijk
!
Integer     :: TUVQPindex(nTUVP*nTUVQ),ituvQP,ntuvp2,ntuvq2
Integer     :: jP,jQ,tP,uP,vP,tQ,uQ,vQ,ituvP,ituvQ,ituvPQ,iOrbP,ijkQ,M
Integer     :: iPrimPQ,iOFF,lupri,iituvP,iituvQ,inTUVPQ,iE,iPrimQ,iprimP,iprimQP
Real(realk),parameter :: D2=2.0E0_realk,D1=1.0E0_realk
Real(realk),pointer  :: TMP(:)
!print*,'MEM_WORKPOINTER_ALLOC(TMP,nPrim)',nPrim,nPrimQ,nPrimP,nPasses
CALL MEM_WORKPOINTER_ALLOC(TMP,nPrim)
IF(Q%contractbasis)then
   iPrimQ=0
   DO iPasses=1,nPasses
      DO iPrim2=1,nPrim2
         pref2=signQ/(D2*Q%orbital2%exponents(iPrim2))**l2
         DO iPrim1=1,nPrim1
            iPrimQ = iPrimQ + 1
            pref12=D1/(D2*Q%orbital1%exponents(iPrim1))**l1*pref2
            DO iPrimP=1,nPrimP
               iPrimQP = iPrimQ+(iPrimP-1)*nPrimQ
               TMP(iPrimQP)=pref12
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ELSE
!PS_int
   iPrimQ=0
   DO iPasses=1,nPasses
      DO iPrim2=1,nPrim2
         pref2=signQ/(D2*Q%orbital2%exponents(iPrim2))**l2
         DO iPrim1=1,nPrim1
            iPrimQ = iPrimQ + 1
            pref12=D1/(D2*Q%orbital1%exponents(iPrim1))**l1*pref2*preexpfac(iPrimQ,iAngmomQ)
            DO iPrimP=1,nPrimP
               iPrimQP = iPrimQ+(iPrimP-1)*nPrimQ
               TMP(iPrimQP)=pref12
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDIF
!FIRST WE FIND THE TUVQP INDEX BECAUSE FINDING THINGS IN TUVindex IS EXPENSIVE
iOFF = (ideriv-1)*ntuvPQ-ioffPQ
ituvQP = 0
ituvP = 0
DO jP = startP,endP
   DO tP=jP,0,-1
      DO uP=jP-tP,0,-1
         vP=jP-tP-uP
         ituvP=ituvP+1
         ituvQ=0
         DO jQ = startQ,endQ
            DO tQ=jQ,0,-1
               DO uQ=jQ-tQ,0,-1
                  vQ=jQ-tQ-uQ
                  ituvQ=ituvQ+1
                  iTUVQP = iTUVQP+1 
                  TUVQPindex(iTUVQP) = TUVindex(tP+tQ,uP+uQ,vP+vQ)+iOFF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
ntuvp2=ituvp
ntuvq2=ituvq

DO ituvP=1,ntuvP2
   DO ituvQ=1,ntuvQ2
      ituvQP=TUVQPindex(ituvQ+(ituvP-1)*ntuvQ2) 
      DO iPrimP=1,nPrimP
         iPrimQP = (iPrimP-1)*nPrimQ
         DO iPrimQ=1,nPrimQ
            OUT(iPrimQ,iPrimP,ituvP,ituvQ) = WTUV2(iPrimQP+iPrimQ,ituvQP)*TMP(iPrimQP+iPrimQ)
         ENDDO
      ENDDO
   ENDDO
ENDDO
CALL MEM_WORKPOINTER_DEALLOC(TMP)

END SUBROUTINE DirectSingleHermiteEcoeffgen

!> \brief direct contraction of single hermite Ecoefficients for the right hand side overlap Q
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param OUT the output intermidiates
!> \param WTUV2 the primitive hermite integrals 
!> \param nPrim total number of primitive
!> \param nPrimP number of primitive for overlap P
!> \param nPrimQ number of primitive for overlap Q
!> \param startP starting angular momentum for overlap P
!> \param endP end angular momentum for overlap P
!> \param startQ starting angular momentum for overlap Q
!> \param endQ end angular momentum for overlap Q
!> \param ioffP offset for tuvindexing for overlap P
!> \param ioffQ offset for tuvindexing for overlap Q
!> \param ioffPQ offset for tuvindexing for overlap PQ
!> \param nTUVEFGPQ total number of angular components
!> \param ntuvPQ number of angular components for PQ
!> \param ntuvP number of angular components for P
!> \param ntuvQ number of angular components for Q
!> \param ideriv derivative index
!> \param Q the overlap distribution
!> \param signQ the sign for the Q overlap
!> \param l1 angular momentum for orbital 1
!> \param l2 angular momentum for orbital 2
!> \param nPrim1 number of primitives for orbital 1
!> \param nPrim2 number of primitives for orbital 2
!> \param nPasses number of passes
!> \param LUPRI the logical unit number for the output file
SUBROUTINE DirectSingleHermiteEcoeffseg(OUT,WTUV2,TUVindex,nPrim,nPrimP,nPrimQ,startP,endP,startQ,endQ,&
     & ioffP,ioffQ,ioffPQ,nTUVEFGPQ,ntuvPQ,ntuvP,ntuvQ,ideriv,Q,signQ,l1,l2,nPrim1,nPrim2,nPasses,&
     & preexpfac,nPrimAlloc,nAngAlloc,iAngmomQ,lupri)
implicit none
TYPE(Overlap)      :: Q
Integer         :: nPrim,startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ
Integer         :: ntuvEFGPQ,ntuvPQ,ntuvP,ntuvQ,ideriv,nPrimP,nPrimQ,STARTE
Integer         :: j1,j2,nprim1,nprim2,npasses,ipasses,iprim1,iprim2
Integer         :: nPrimAlloc,nAngAlloc,iAngmomQ
Real(realk)     :: WTUV2(nPrimQ,nPrimP,nPasses,ntuvEFGPQ)!,tuvTUV(ntuvQ,ntuvP,nPrim)
Real(realk)     :: OUT(nPrimP,nPasses,ntuvP,ntuvQ),signQ,pref2,pref12
Real(realk)     :: preexpfac(nPrimAlloc,nAngAlloc)
Integer,pointer :: TUVindex(:,:,:)
!Real(realk)     :: Ec(nPrimQ)
Integer             :: l1,l2,ijk
!
Integer     :: TUVQPindex(nTUVP*nTUVQ),ituvQP,ntuvp2,ntuvq2,offset
Integer     :: jP,jQ,tP,uP,vP,tQ,uQ,vQ,ituvP,ituvQ,ituvPQ,iOrbP,ijkQ,M
Integer     :: iPrimPQ,iOFF,lupri,iituvP,iituvQ,inTUVPQ,iE,iPrimQ,iprimP,iprimQP,iQ
Real(realk),parameter :: D2=2.0E0_realk,D1=1.0E0_realk,D0=0E0_realk
Real(realk),pointer  :: prefactor(:)
Real(realk)          :: tmp
!print*,'debug 53 prefactor',nPrimQ*nPasses
CALL MEM_WORKPOINTER_ALLOC(prefactor,nPrimQ*nPasses)
DO iPasses=1,nPasses
   iPrimQ=0
   offset = (iPasses-1)*nPrimQ
   DO iPrim2=1,nPrim2
      pref2=signQ/(D2*Q%orbital2%exponents(iPrim2))**l2
      DO iPrim1=1,nPrim1
         iPrimQ = iPrimQ + 1
         pref12=D1/(D2*Q%orbital1%exponents(iPrim1))**l1*pref2*preexpfac(iPrimQ+offset,iAngmomQ)
         prefactor(iPrimQ+offset)=pref12
      ENDDO
   ENDDO
ENDDO
!FIRST WE FIND THE TUVQP INDEX BECAUSE FINDING THINGS IN TUVindex IS EXPENSIVE
iOFF = (ideriv-1)*ntuvPQ-ioffPQ
ituvQP = 0
ituvP = 0
DO jP = startP,endP
   DO tP=jP,0,-1
      DO uP=jP-tP,0,-1
         vP=jP-tP-uP
         ituvP=ituvP+1
         ituvQ=0
         DO jQ = startQ,endQ
            DO tQ=jQ,0,-1
               DO uQ=jQ-tQ,0,-1
                  vQ=jQ-tQ-uQ
                  ituvQ=ituvQ+1
                  iTUVQP = iTUVQP+1 
                  TUVQPindex(iTUVQP) = TUVindex(tP+tQ,uP+uQ,vP+vQ)+iOFF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
ntuvp2=ituvp
ntuvq2=ituvq

DO ituvP=1,ntuvP2
  DO ituvQ=1,ntuvQ2
    ituvQP=TUVQPindex(ituvQ+(ituvP-1)*ntuvQ2) 
    DO iPasses=1,nPasses
      offset = (iPasses-1)*nPrimQ
      DO iPrimP=1,nPrimP
         tmp = WTUV2(1,iPrimP,iPasses,ituvQP)*prefactor(1+offset)
         DO iPrimQ=2,nPrimQ
            tmp = tmp + WTUV2(iPrimQ,iPrimP,iPasses,ituvQP)*prefactor(iPrimQ+offset)
         ENDDO
         OUT(iPrimP,iPasses,ituvP,ituvQ) = tmp
      ENDDO
    ENDDO
  ENDDO
ENDDO

CALL MEM_WORKPOINTER_DEALLOC(prefactor)

END SUBROUTINE DirectSingleHermiteEcoeffseg

!> \brief direct contraction of hermite Ecoefficients for the right hand side overlap Q
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param OUT the output intermidiates
!> \param WTUV3 the primitive hermite integrals 
!> \param TUVindex for a given t,u,v values it gives the TUV index (t,u,v)=(000,100,010,001,.. etc) means TUV = (1,2,3,4,..)
!> \param nPrim total number of primitive
!> \param nPrimP number of primitive for overlap P
!> \param nPrimQ number of primitive for overlap Q
!> \param startP starting angular momentum for overlap P
!> \param endP end angular momentum for overlap P
!> \param startQ starting angular momentum for overlap Q
!> \param endQ end angular momentum for overlap Q
!> \param ioffP offset for tuvindexing for overlap P
!> \param ioffQ offset for tuvindexing for overlap Q
!> \param ioffPQ offset for tuvindexing for overlap PQ
!> \param nTUVEFGPQ total number of angular components
!> \param ntuvPQ number of angular components for PQ
!> \param ntuvP number of angular components for P
!> \param ntuvQ number of angular components for Q
!> \param ideriv derivative index
!> \param Ecoeffs the ecoefficients 
!> \param startE the starting index for Ecoeffs
!> \param ijk cartesian anugular components
!> \param l1 angular momentum for orbital 1
!> \param l2 angular momentum for orbital 2
!> \param LUPRI the logical unit number for the output file
SUBROUTINE DirectcontractEQgen(OUT,WTUV3,TUVindex,nPrim,nPrimP,nPrimQ,startP,endP,startQ,endQ,ioffP,&
     & ioffQ,ioffPQ,nTUVEFGPQ,ntuvPQ,ntuvP,ntuvQ,ideriv,Ecoeffs,startE,ijk,l1,l2,SYMindex,iPQxyz,lupri)
implicit none
Integer         :: nPrim,startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ,nPrimQ
Integer         :: ntuvEFGPQ,ntuvPQ,ntuvP,ntuvQ,ideriv,nPrimP,STARTE,iPQxyz
Integer         :: l1,l2,ijk
Real(realk)     :: WTUV3(nPrimQ,nPrimP,ntuvEFGPQ)
Real(realk)     :: OUT(nPrimQ,nPrimP,ntuvP,ijk)
Integer,pointer :: TUVindex(:,:,:),SYMINDEX(:)
Real(realk)     :: Ecoeffs(nPrimQ,ntuvQ,ijk)
Real(realk)     :: Ec(nPrimQ)
!
Real(realk) :: maxEcont,THRESHOLD,Etmp
Integer     :: startec,TUVQPindex(ntuvP*ntuvQ),ntuvp2,ntuvq2
Integer     :: jP,jQ,tP,uP,vP,tQ,uQ,vQ,ituvP,ituvQ,iOrbP,ijkQ,M,ituvQP
Integer     :: iPrimPQ,iOFF,lupri,iituvP,iituvQ,inTUVPQ,iE,iPrimQ,iprimP,iprimQP
logical     :: addcontribution
!
!FIRST WE FIND THE TUVQP INDEX BECAUSE FINDING THINGS IN TUVindex IS EXPENSIVE
iOFF = (ideriv-1)*ntuvPQ-ioffPQ
ituvQP = 0
ituvQ = 0
DO jQ = startQ,endQ
   DO tQ=jQ,0,-1
      DO uQ=jQ-tQ,0,-1
         vQ=jQ-tQ-uQ
         ituvQ = ituvQ+1
         ituvP = 0
         DO jP = startP,endP
            DO tP=jP,0,-1
               DO uP=jP-tP,0,-1
                  vP=jP-tP-uP
                  ituvP = ituvP+1
                  iTUVQP = iTUVQP+1 
                  TUVQPindex(iTUVQP) = TUVindex(tP+tQ,uP+uQ,vP+vQ)+iOFF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
ntuvp2=ituvp
ntuvq2=ituvq
   
IF(iPQxyz.EQ. 0)THEN !no symmetry the general case
   THRESHOLD = 1.0E-15_realk
   DO ijkQ = 1, ijk
      maxEcont = 0E0_realk
      DO iPrimQ=1,nPrimQ
         Etmp = Ecoeffs(iPrimQ,1,ijkQ)
         Ec(iPrimQ) = Etmp
         maxEcont = MAX(maxEcont,ABS(Etmp))
      ENDDO
      IF(maxEcont .GE. THRESHOLD) THEN
         DO ituvP = 1,ntuvP2
            ituvQP=TUVQPindex(ituvP)
            DO iPrimQ=1,nPrimQ
               Etmp = Ec(iPrimQ)
               DO iPrimP=1,nPrimP 
                  OUT(iPrimQ,iPrimP,ituvP,ijkQ) = WTUV3(iPrimQ,iPrimP,ituvQP)*Etmp
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO ituvP = 1,ntuvP2
            DO iPrimP=1,nPrimP 
               DO iPrimQ=1,nPrimQ
                  OUT(iPrimQ,iPrimP,ituvP,ijkQ) = 0E0_realk
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      DO ituvQ = 2,ntuvQ2
         maxEcont = 0E0_realk
         DO iPrimQ=1,nPrimQ
            Etmp = Ecoeffs(iPrimQ,ituvQ,ijkQ)
            Ec(iPrimQ) = Etmp
            maxEcont = MAX(maxEcont,ABS(Etmp))
         ENDDO
         IF(maxEcont .LT. THRESHOLD)CYCLE
         IOFF = (ituvQ-1)*ntuvP2
         DO ituvP = 1,ntuvP2
            ituvQP=TUVQPindex(ituvP+IOFF)
            DO iPrimQ=1,nPrimQ
               Etmp = Ec(iPrimQ)
               DO iPrimP=1,nPrimP 
                  OUT(iPrimQ,iPrimP,ituvP,ijkQ) = OUT(iPrimQ,iPrimP,ituvP,ijkQ)+WTUV3(iPrimQ,iPrimP,ituvQP)*Etmp
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ELSE !special case where we exploits local symmetry
   THRESHOLD = 1.0E-15_realk
   DO ijkQ = 1, ijk
      maxEcont = 0E0_realk
      DO iPrimQ=1,nPrimQ
         Etmp = Ecoeffs(iPrimQ,1,ijkQ)
         Ec(iPrimQ) = Etmp
         maxEcont = MAX(maxEcont,ABS(Etmp))
      ENDDO
      IF(maxEcont .GE. THRESHOLD) THEN
         DO ituvP = 1,ntuvP2
            ituvQP=TUVQPindex(ituvP)
            IF(SYMINDEX(ITUVQP).EQ. 0)THEN
               DO iPrimQ=1,nPrimQ
                  Etmp = Ec(iPrimQ)
                  DO iPrimP=1,nPrimP 
                     OUT(iPrimQ,iPrimP,ituvP,ijkQ) = WTUV3(iPrimQ,iPrimP,ituvQP)*Etmp
                  ENDDO
               ENDDO
            ELSE
               DO iPrimP=1,nPrimP 
                  DO iPrimQ=1,nPrimQ
!                     IF(ABS(WTUV3(iPrimQ,iPrimP,ituvQP)).GT. 1E-11_realk)THEN
!                        CALL LSQUIT('local sym error1',-1)
!                     ENDIF
                     OUT(iPrimQ,iPrimP,ituvP,ijkQ) = 0E0_realk
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ELSE
         DO ituvP = 1,ntuvP2
            DO iPrimP=1,nPrimP 
               DO iPrimQ=1,nPrimQ
                  OUT(iPrimQ,iPrimP,ituvP,ijkQ) = 0E0_realk
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      DO ituvQ = 2,ntuvQ2
         maxEcont = 0E0_realk
         DO iPrimQ=1,nPrimQ
            Etmp = Ecoeffs(iPrimQ,ituvQ,ijkQ)
            Ec(iPrimQ) = Etmp
            maxEcont = MAX(maxEcont,ABS(Etmp))
         ENDDO
         IF(maxEcont .LT. THRESHOLD)CYCLE
         IOFF = (ituvQ-1)*ntuvP2
         DO ituvP = 1,ntuvP2
            ituvQP=TUVQPindex(ituvP+IOFF)
            IF(SYMINDEX(ITUVQP).EQ. 0)THEN
               DO iPrimQ=1,nPrimQ
                  Etmp = Ec(iPrimQ)
                  DO iPrimP=1,nPrimP 
                     OUT(iPrimQ,iPrimP,ituvP,ijkQ) = OUT(iPrimQ,iPrimP,ituvP,ijkQ)+WTUV3(iPrimQ,iPrimP,ituvQP)*Etmp
                  ENDDO
               ENDDO
!            ELSE
!               DO iPrimP=1,nPrimP 
!                  DO iPrimQ=1,nPrimQ
!                     IF(ABS(WTUV3(iPrimQ,iPrimP,ituvQP)).GT. 1E-11_realk)THEN
!                        WRITE(lupri,*)'ituvQP',ituvQP
!                        WRITE(lupri,*)'WTUV3(iPrimQ,iPrimP,ituvQP)',WTUV3(iPrimQ,iPrimP,ituvQP)
!                        WRITE(lupri,*)'SYMINDEX(ituvQP)',SYMINDEX(ITUVQP)
!
!                        CALL LSQUIT('local sym error2',-1)
!                     ENDIF
!                  ENDDO
!               ENDDO
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDIF

!!$IF(IPRINT.GT. 10)THEN
!!$   WRITE(LUPRI,*)'DirectcontractEQgen'
!!$   WRITE(LUPRI,'(1X,A6,1X,A6,1X,A6,12X,A6 )') 'iPrimP','ituvP ',' ijkQ ','iPrimQ'
!!$   DO ijkQ = 1, ijk
!!$    DO ituvP = 1,ntuvP2
!!$     DO iPrimP=1,nPrimP 
!!$        WRITE(LUPRI,'(1X,I4,2X,I4,2X,I4,5F9.4/,(17X,5F9.4))') iPrimP,ituvP,ijkQ,(OUT(iPrimQ,iPrimP,ituvP,ijkQ),iPrimQ=1,nPrimQ)
!!$     ENDDO
!!$    ENDDO
!!$   ENDDO
!!$ENDIF

END SUBROUTINE DirectcontractEQgen

!> \brief direct contraction of hermite Ecoefficients for the right hand side overlap Q
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param OUT the output intermidiates
!> \param WTUV3 the primitive hermite integrals 
!> \param TUVindex for a given t,u,v values it gives the TUV index (t,u,v)=(000,100,010,001,.. etc) means TUV = (1,2,3,4,..)
!> \param nPrim total number of primitive
!> \param nPrimP number of primitive for overlap P
!> \param nPrimQ number of primitive for overlap Q
!> \param startP starting angular momentum for overlap P
!> \param endP end angular momentum for overlap P
!> \param startQ starting angular momentum for overlap Q
!> \param endQ end angular momentum for overlap Q
!> \param ioffP offset for tuvindexing for overlap P
!> \param ioffQ offset for tuvindexing for overlap Q
!> \param ioffPQ offset for tuvindexing for overlap PQ
!> \param nTUVEFGPQ total number of angular components
!> \param ntuvPQ number of angular components for PQ
!> \param ntuvP number of angular components for P
!> \param ntuvQ number of angular components for Q
!> \param ideriv derivative index
!> \param Ecoeffs the ecoefficients 
!> \param startE the starting index for Ecoeffs
!> \param ijk cartesian anugular components
!> \param l1 angular momentum for orbital 1
!> \param l2 angular momentum for orbital 2
!> \param LUPRI the logical unit number for the output file
SUBROUTINE DirectcontractEQseg(OUT,WTUV3,TUVindex,nPrim,nPrimP,nPrimQ,nPasses,startP,endP,startQ,endQ,&
     & ioffP,ioffQ,ioffPQ,nTUVEFGPQ,ntuvPQ,ntuvP,ntuvQ,ideriv,Ecoeffs,startE,ijk,l1,l2,SYMindex,iPQxyz,lupri)
implicit none
Integer         :: nPrim,startP,endP,startQ,endQ,ioffP,ioffQ,ioffPQ,nPrimQ,nPasses
Integer         :: ntuvEFGPQ,ntuvPQ,ntuvP,ntuvQ,ideriv,nPrimP,STARTE,iPQxyz
Integer             :: l1,l2,ijk
Real(realk)     :: WTUV3(nPrimQ,nPasses,nPrimP,ntuvEFGPQ)
Real(realk)     :: OUT(nPasses,nPrimP,ntuvP,ijk)
Integer,pointer :: TUVindex(:,:,:),SYMINDEX(:)
Real(realk)     :: Ecoeffs(nPrimQ,nPasses,ntuvQ,ijk)
Real(realk)     :: Ec(nPrimQ,nPasses)
!
Real(realk) :: maxEcont,THRESHOLD,Etmp,tmp
Integer     :: startec,TUVQPindex(ntuvP*ntuvQ),ntuvp2,ntuvq2
Integer     :: jP,jQ,tP,uP,vP,tQ,uQ,vQ,ituvP,ituvQ,iOrbP,ijkQ,M,ituvQP
Integer     :: iPrimPQ,iOFF,lupri,iituvP,iituvQ,inTUVPQ,iE,iPrimQ,iprimP,iprimQP,iPasses
logical     :: addcontribution
real(realk),parameter :: D0 = 0E0_realk
!
!FIRST WE FIND THE TUVQP INDEX BECAUSE FINDING THINGS IN TUVindex IS EXPENSIVE
iOFF = (ideriv-1)*ntuvPQ-ioffPQ
ituvQP = 0
ituvQ = 0
DO jQ = startQ,endQ
   DO tQ=jQ,0,-1
      DO uQ=jQ-tQ,0,-1
         vQ=jQ-tQ-uQ
         ituvQ = ituvQ+1
         ituvP = 0
         DO jP = startP,endP
            DO tP=jP,0,-1
               DO uP=jP-tP,0,-1
                  vP=jP-tP-uP
                  ituvP = ituvP+1
                  iTUVQP = iTUVQP+1 
                  TUVQPindex(iTUVQP) = TUVindex(tP+tQ,uP+uQ,vP+vQ)+iOFF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
ntuvp2=ituvp
ntuvq2=ituvq
   
THRESHOLD = 1.0E-15_realk
IF (nPrimQ.GT. 1) THEN !general case
  CALL DirectcontractEQsegPrimQ(OUT,WTUV3,Ec,Ecoeffs,TUVQPindex,nPrimQ,nPasses,nPrimP,&
     &                          ijk,ntuvEFGPQ,ntuvP2,ntuvQ2,ntuvP,ntuvQ,SYMindex,iPQxyz,THRESHOLD)
ELSE IF (nPasses.GT. 1) THEN
   IF(iPQxyz.EQ. 0)THEN !no symmetry the general case
      CALL DirectcontractEQsegPassQ(OUT,WTUV3,Ec,Ecoeffs,TUVQPindex,nPasses,nPrimP,&
           &                          ijk,ntuvEFGPQ,ntuvP2,ntuvQ2,ntuvP,ntuvQ,THRESHOLD)
   ELSE !use local symmetry special case
      CALL DirectcontractEQsegPassQsym(OUT,WTUV3,Ec,Ecoeffs,TUVQPindex,nPasses,nPrimP,&
           &          ijk,ntuvEFGPQ,ntuvP2,ntuvQ2,ntuvP,ntuvQ,SYMindex,iPQxyz,THRESHOLD)
   ENDIF
ELSE 
  CALL DirectcontractEQsegPrimP(OUT,WTUV3,Ecoeffs,TUVQPindex,nPrimP,&
     &                          ijk,ntuvEFGPQ,ntuvP2,ntuvQ2,ntuvP,ntuvQ,THRESHOLD)
ENDIF

!!$IF(IPRINT.GT. 10)THEN
!!$   WRITE(LUPRI,*)'DirectcontractEQseg'
!!$   WRITE(LUPRI,'(1X,A6,1X,A6,1X,A6,12X,A6 )') 'iPrimP','ituvP ',' ijkQ ','iPrimQ'
!!$   DO ijkQ = 1, ijk
!!$    DO ituvP = 1,ntuvP2
!!$     DO iPrimP=1,nPrimP 
!!$        WRITE(LUPRI,'(1X,I4,2X,I4,2X,I4,5F9.4/,(17X,5F9.4))') iPrimP,ituvP,ijkQ,(OUT(iPrimQ,iPrimP,ituvP,ijkQ),iPrimQ=1,nPrimQ)
!!$     ENDDO
!!$    ENDDO
!!$   ENDDO
!!$ENDIF

END SUBROUTINE DirectcontractEQseg

SUBROUTINE DirectcontractEQsegPrimQ(OUT,WTUV3,Ec,Ecoeffs,TUVQPindex,nPrimQ,nPasses,nPrimP,&
     &                    ijk,ntuvEFGPQ,ntuvP2,ntuvQ2,ntuvP,ntuvQ,SYMindex,iPQxyz,THRESHOLD)
implicit none
Integer,intent(IN)         :: nPrimQ,nPasses,nPrimP,ijk,ntuvEFGPQ,ntuvP2,ntuvQ2,ntuvP,ntuvQ
Integer,intent(IN)         :: TUVQPindex(ntuvP*ntuvQ),iPQxyz
Real(realk), intent(INOUT) :: OUT(nPasses,nPrimP,ntuvP,ijk)
Real(realk), intent(IN)    :: WTUV3(nPrimQ,nPasses,nPrimP,ntuvEFGPQ)
Real(realk), intent(INOUT) :: Ec(nPrimQ,nPasses)
Real(realk), intent(IN)    :: Ecoeffs(nPrimQ,nPasses,ntuvQ,ijk)
Real(realk),intent(IN)     :: THRESHOLD
integer,pointer            :: symindex(:)
!
Integer     :: ijkQ,iPasses,iPrimQ,ituvP,ituvQP,iPrimP,ituvQ,IOFF
Real(realk) :: Etmp,maxEcont,tmp
!
DO ijkQ = 1, ijk
   maxEcont = 0E0_realk
   DO iPasses=1,nPasses
      DO iPrimQ=1,nPrimQ
         Etmp = Ecoeffs(iPrimQ,iPasses,1,ijkQ)
         Ec(iPrimQ,iPasses) = Etmp
         maxEcont = MAX(maxEcont,ABS(Etmp))
      ENDDO
   ENDDO
   IF(maxEcont .GE. THRESHOLD)THEN
     DO ituvP = 1,ntuvP2
        ituvQP=TUVQPindex(ituvP)
        IF(SYMINDEX(ITUVQP).EQ. 0)THEN
           DO iPrimP=1,nPrimP
              DO iPasses=1,nPasses
                 tmp = Ec(1,iPasses)*WTUV3(1,iPasses,iPrimP,ituvQP)
                 DO iPrimQ=2,nPrimQ
                    tmp = tmp + Ec(iPrimQ,iPasses)*WTUV3(iPrimQ,iPasses,iPrimP,ituvQP)
                 ENDDO
                 OUT(iPasses,iPrimP,ituvP,ijkQ) = tmp
              ENDDO
           ENDDO
        ELSE
           DO iPrimP=1,nPrimP
              DO iPasses=1,nPasses
                 OUT(iPasses,iPrimP,ituvP,ijkQ) = 0E0_realk
              ENDDO
           ENDDO
        ENDIF
     ENDDO
   ELSE
     OUT(:,:,:,ijkQ) = 0E0_realk
   ENDIF
   DO ituvQ = 2,ntuvQ2
      maxEcont = 0E0_realk
      DO iPasses=1,nPasses
         DO iPrimQ=1,nPrimQ
            Etmp = Ecoeffs(iPrimQ,iPasses,ituvQ,ijkQ)
            Ec(iPrimQ,iPasses) = Etmp
            maxEcont = MAX(maxEcont,ABS(Etmp))
         ENDDO
      ENDDO
      IF(maxEcont .LT. THRESHOLD)CYCLE
      IOFF = (ituvQ-1)*ntuvP2
      DO ituvP = 1,ntuvP2
         ituvQP=TUVQPindex(ituvP+IOFF)
         IF(SYMINDEX(ITUVQP).EQ. 0)THEN
            DO iPrimP=1,nPrimP
               DO iPasses=1,nPasses
                  tmp = Ec(1,iPasses)*WTUV3(1,iPasses,iPrimP,ituvQP)
                  DO iPrimQ=2,nPrimQ
                     tmp = tmp + Ec(iPrimQ,iPasses)*WTUV3(iPrimQ,iPasses,iPrimP,ituvQP)
                  ENDDO
                  OUT(iPasses,iPrimP,ituvP,ijkQ) = OUT(iPasses,iPrimP,ituvP,ijkQ) + tmp
               ENDDO
            ENDDO
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE DirectcontractEQsegPrimQ

SUBROUTINE DirectcontractEQsegPassQ(OUT,WTUV3,Ec,Ecoeffs,TUVQPindex,nPasses,nPrimP,&
     &                              ijk,ntuvEFGPQ,ntuvP2,ntuvQ2,ntuvP,ntuvQ,THRESHOLD)
implicit none
Integer,intent(IN)         :: nPasses,nPrimP,ijk,ntuvEFGPQ,ntuvP2,ntuvQ2,ntuvP,ntuvQ
Integer,intent(IN)         :: TUVQPindex(ntuvP*ntuvQ)
Real(realk), intent(INOUT) :: OUT(nPasses,nPrimP,ntuvP,ijk)
Real(realk), intent(IN)    :: WTUV3(nPasses,nPrimP,ntuvEFGPQ)
Real(realk), intent(INOUT) :: Ec(nPasses)
Real(realk), intent(IN)    :: Ecoeffs(nPasses,ntuvQ,ijk)
Real(realk),intent(IN)     :: THRESHOLD
!
Integer     :: ijkQ,iPasses,ituvP,ituvQP,iPrimP,ituvQ,IOFF
Real(realk) :: Etmp,maxEcont,tmp
!
DO ijkQ = 1, ijk
!  ituvQ = 1
   maxEcont = 0E0_realk
   DO iPasses=1,nPasses
       Etmp = Ecoeffs(iPasses,1,ijkQ)
       Ec(iPasses) = Etmp
       maxEcont = MAX(maxEcont,ABS(Etmp))
   ENDDO
   IF(maxEcont .GE. THRESHOLD)THEN
     DO ituvP = 1,ntuvP2
        ituvQP=TUVQPindex(ituvP)
        DO iPrimP=1,nPrimP
           DO iPasses=1,nPasses
              OUT(iPasses,iPrimP,ituvP,ijkQ) = Ec(iPasses)*WTUV3(iPasses,iPrimP,ituvQP)
           ENDDO
        ENDDO
     ENDDO
   ELSE
     OUT(:,:,:,ijkQ) = 0E0_realk
   ENDIF
!  ituvQ > 1
   DO ituvQ = 2,ntuvQ2
      maxEcont = 0E0_realk
      DO iPasses=1,nPasses
         Etmp = Ecoeffs(iPasses,ituvQ,ijkQ)
         Ec(iPasses) = Etmp
         maxEcont = MAX(maxEcont,ABS(Etmp))
      ENDDO
      IF(maxEcont .LT. THRESHOLD)CYCLE
      IOFF = (ituvQ-1)*ntuvP2
      DO ituvP = 1,ntuvP2
         ituvQP=TUVQPindex(ituvP+IOFF)
         DO iPrimP=1,nPrimP
            DO iPasses=1,nPasses
               OUT(iPasses,iPrimP,ituvP,ijkQ) = OUT(iPasses,iPrimP,ituvP,ijkQ) &
     &              + Ec(iPasses)*WTUV3(iPasses,iPrimP,ituvQP)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE DirectcontractEQsegPassQ

SUBROUTINE DirectcontractEQsegPassQsym(OUT,WTUV3,Ec,Ecoeffs,TUVQPindex,nPasses,nPrimP,&
     &               ijk,ntuvEFGPQ,ntuvP2,ntuvQ2,ntuvP,ntuvQ,SYMindex,iPQxyz,THRESHOLD)
implicit none
Integer,intent(IN)         :: nPasses,nPrimP,ijk,ntuvEFGPQ,ntuvP2,ntuvQ2,ntuvP,ntuvQ
Real(realk), intent(INOUT) :: OUT(nPasses,nPrimP,ntuvP,ijk)
Real(realk), intent(IN)    :: WTUV3(nPasses,nPrimP,ntuvEFGPQ)
Real(realk), intent(INOUT) :: Ec(nPasses)
Real(realk), intent(IN)    :: Ecoeffs(nPasses,ntuvQ,ijk)
Integer,intent(IN)         :: TUVQPindex(ntuvP*ntuvQ),iPQxyz
Real(realk),intent(IN)     :: THRESHOLD
integer,pointer            :: symindex(:)
!
Integer     :: ijkQ,iPasses,ituvP,ituvQP,iPrimP,ituvQ,IOFF
Real(realk) :: Etmp,maxEcont,tmp
!
DO ijkQ = 1, ijk
!  ituvQ = 1
   maxEcont = 0E0_realk
   DO iPasses=1,nPasses
       Etmp = Ecoeffs(iPasses,1,ijkQ)
       Ec(iPasses) = Etmp
       maxEcont = MAX(maxEcont,ABS(Etmp))
   ENDDO
   IF(maxEcont .GE. THRESHOLD)THEN
     DO ituvP = 1,ntuvP2
        ituvQP=TUVQPindex(ituvP)
        IF(SYMINDEX(ITUVQP).EQ. 0)THEN
           DO iPrimP=1,nPrimP
              DO iPasses=1,nPasses
                 OUT(iPasses,iPrimP,ituvP,ijkQ) = Ec(iPasses)*WTUV3(iPasses,iPrimP,ituvQP)
              ENDDO
           ENDDO
        ELSE
           DO iPrimP=1,nPrimP
              DO iPasses=1,nPasses
                 OUT(iPasses,iPrimP,ituvP,ijkQ) = 0E0_realk
              ENDDO
           ENDDO
        ENDIF
     ENDDO
   ELSE
     DO ituvP = 1,ntuvP2
        DO iPrimP=1,nPrimP
           DO iPasses=1,nPasses
              OUT(iPasses,iPrimP,ituvP,ijkQ) = 0E0_realk
           ENDDO
        ENDDO
     ENDDO
   ENDIF
!  ituvQ > 1
   DO ituvQ = 2,ntuvQ2
      maxEcont = 0E0_realk
      DO iPasses=1,nPasses
         Etmp = Ecoeffs(iPasses,ituvQ,ijkQ)
         Ec(iPasses) = Etmp
         maxEcont = MAX(maxEcont,ABS(Etmp))
      ENDDO
      IF(maxEcont .LT. THRESHOLD)CYCLE
      IOFF = (ituvQ-1)*ntuvP2
      DO ituvP = 1,ntuvP2
         ituvQP=TUVQPindex(ituvP+IOFF)
         IF(SYMINDEX(ITUVQP).EQ. 0)THEN
            DO iPrimP=1,nPrimP
               DO iPasses=1,nPasses
                  OUT(iPasses,iPrimP,ituvP,ijkQ) = OUT(iPasses,iPrimP,ituvP,ijkQ) &
                       &              + Ec(iPasses)*WTUV3(iPasses,iPrimP,ituvQP)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE DirectcontractEQsegPassQsym

SUBROUTINE DirectcontractEQsegPrimP(OUT,WTUV3,Ecoeffs,TUVQPindex,nPrimP,&
     & ijk,ntuvEFGPQ,ntuvP2,ntuvQ2,ntuvP,ntuvQ,THRESHOLD)
implicit none
Integer,intent(IN)         :: nPrimP,ijk,ntuvEFGPQ,ntuvP2,ntuvQ2,ntuvP,ntuvQ
Real(realk), intent(INOUT) :: OUT(nPrimP,ntuvP,ijk)
Real(realk), intent(IN)    :: WTUV3(nPrimP,ntuvEFGPQ)
Real(realk), intent(IN)    :: Ecoeffs(ntuvQ,ijk)
Integer,intent(IN)         :: TUVQPindex(ntuvP*ntuvQ)
Real(realk),intent(IN)     :: THRESHOLD
!
Integer     :: ijkQ,ituvP,ituvQP,iPrimP,ituvQ,IOFF
Real(realk) :: Ec,maxEcont,tmp
!
DO ijkQ = 1, ijk
   Ec = Ecoeffs(1,ijkQ)
   maxEcont = ABS(Ec)
   IF(maxEcont .GE. THRESHOLD)THEN
     DO ituvP = 1,ntuvP2
        ituvQP=TUVQPindex(ituvP)
        DO iPrimP=1,nPrimP
           OUT(iPrimP,ituvP,ijkQ) = Ec*WTUV3(iPrimP,ituvQP)
        ENDDO
     ENDDO
   ELSE
     OUT(:,:,ijkQ) = 0E0_realk
   ENDIF
   DO ituvQ = 2,ntuvQ2
      Ec = Ecoeffs(ituvQ,ijkQ)
      maxEcont = ABS(Ec)
      IF(maxEcont .LT. THRESHOLD)CYCLE
      IOFF = (ituvQ-1)*ntuvP2
      DO ituvP = 1,ntuvP2
         ituvQP=TUVQPindex(ituvP+IOFF)
         DO iPrimP=1,nPrimP
           OUT(iPrimP,ituvP,ijkQ) = OUT(iPrimP,ituvP,ijkQ) + Ec*WTUV3(iPrimP,ituvQP)
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE DirectcontractEQsegPrimP

!> \brief add intermidiate to Integral%integrals  
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param Integral contains arrays to store intermidiates and final integrals
!> \param PQ contains integrand info like reduced exponents, integral prefactor, etc.
!> \param iAngmom the index of angular momentum
!> \param iOrder the current derivative components
!> \param nOrder the number of derivative components to process on by one
!> \param nDerivQ the number of derivative components to process simultaneously
!> \param contAng specifies contrated angular or angular contracted AO ordering
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE AddToTUVQ(Integral,PQ,iAngmom,iOrder,nOrder,nDerivQ,contAng,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Integrand)    :: PQ
Integer            :: iAngmom,iOrder,nOrder,nDerivQ,LUPRI,IPRINT 
logical,intent(IN) :: contAng
!
Integer :: nPrimP,nTUVP,nP,nCompQ,nContQ,nContC,nContD,nCompC,nCompD
!
Integer :: nC,nD,startC,startD,nPassQ,angC,AngD,nDer,nDerTUV,startDer,endDer
Integer :: startP,endP
Integer,target  :: iDer,iOne=1
Integer,pointer :: jDer
TYPE(Overlap),pointer :: Q,P
Logical :: permute

P => PQ%P%p
nPrimP  = P%nPrimitives*PQ%P%p%nPasses
startP  = P%startAngmom
endP    = P%endAngmom
IF (integral%dohodi) THEN
  endP = P%maxAngmom + Integral%lhsGeoOrder
  IF (P%type_hermite_single) startP = endP
ENDIF
  
nTUVP   = (endP+1)*(endP+2)*(endP+3)/6 - startP*(startP+1)*(startP+2)/6

nP      = nPrimP*nTUVP
Q => PQ%Q%p

nCompQ  = Q%nOrbComp(iAngmom)
nContQ  = Q%nContracted(iAngmom)*Q%nPasses
angC    = Q%indexAng1(iAngmom)
angD    = Q%indexAng2(iAngmom)
nC      = Q%orbital1%totOrbitals
nContC  = Q%orbital1%nContracted(angC)
startC  = Q%orbital1%startLocOrb(angC)
nCompC  = Q%orbital1%nOrbComp(angC)
nD      = Q%orbital2%totOrbitals
nContD  = Q%orbital2%nContracted(angD)
startD  = Q%orbital2%startLocOrb(angD)
nCompD  = Q%orbital2%nOrbComp(angD)
nPassQ  = Q%nPasses
permute = Q%sameAO.AND.(startC.NE.startD)

#ifdef VAR_LSDEBUGINT
IF(nP*nC*nD*nPassQ*nOrder.GT.allocIntmaxTUVdim)THEN
   call lsquit('AddToTUVQ1 alloc error1',lupri)
ENDIF
IF(nContQ*nP*nCompQ.GT.allocIntmaxTUVdim)THEN
   writE(*,*) 'nContQ            = ',nContQ
   writE(*,*) 'nP                = ',nP    
   writE(*,*) 'nCompQ            = ',nCompQ
   writE(*,*) 'allocIntmaxTUVdim = ',allocIntmaxTUVdim
   call lsquit('AddToTUVQ1 alloc error2',lupri)
ENDIF
#endif

IF (nDerivQ.GT.1) THEN
  nDer     = nDerivQ
  nDerTUV  = nDerivQ
  startDer = 1
  endDer   = nDerivQ
  jDer     => iDer
ELSE
  nDer     = nOrder
  nDerTUV  = 1
  startDer = iOrder
  endDer   = iOrder
  jDer     => iOne
ENDIF

IF (.NOT.contAng) THEN !Default AO ordering, angular first contracted second
  DO iDer=startDer,endDer
    if(.not.permute)then !default, asssuming permute is false
       IF(PQ%Q%p%segmented)THEN
          CALL AddToTUVQ1seg(Integral%tuvQ,Integral%IN,nP,nCompQ,nContQ,&
               &          nC,nD,startC,startD,nContC,nContD,nCompC,nCompD,nPassQ,&
               &          iDer,nDer,jDer,nDerTUV,LUPRI,IPRINT)
       ELSE
          CALL AddToTUVQ1(Integral%tuvQ,Integral%IN,nP,nCompQ,nContQ,&
               &          nC,nD,startC,startD,nContC,nContD,nCompC,nCompD,nPassQ,&
               &          iDer,nDer,jDer,nDerTUV,LUPRI,IPRINT)
       ENDIF
    else
       !general case 
       CALL AddToTUVQ1gen(Integral%tuvQ,Integral%IN,nP,nCompQ,nContQ,&
            &          nC,nD,startC,startD,nContC,nContD,nCompC,nCompD,nPassQ,&
            &          iDer,nDer,jDer,nDerTUV,permute,LUPRI,IPRINT)
    endif
     IF (IPRINT.GT. 15) THEN
       CALL PrintTuvQ(Integral%tuvQ,nP,nC,nD,nPassQ,iDer,nDer,LUPRI)
     ENDIF
  ENDDO
ELSE !Optional AO ordering, contraced first angular second
  DO iDer=startDer,endDer
     !general case 
     CALL AddToTUVQ1gen_ca(Integral%tuvQ,Integral%IN,nP,nCompQ,nContQ,&
          &          nC,nD,startC,startD,nContC,nContD,nCompC,nCompD,nPassQ,&
          &          iDer,nDer,jDer,nDerTUV,permute,LUPRI,IPRINT)
     IF (IPRINT.GT. 15) THEN
       CALL PrintTuvQ(Integral%tuvQ,nP,nC,nD,nPassQ,iDer,nDer,LUPRI)
     ENDIF
  ENDDO
ENDIF

END SUBROUTINE AddToTUVQ

SUBROUTINE AddToTUVQ1(tuvCD,tuvQ,nP,nCompQ,nContQ,nC,nD,startC,startD,&
     &                nContC,nContD,nCompC,nCompD,nPassQ,iOrder,nOrder,iDer,nDer,&
     &                LUPRI,IPRINT)
implicit none
Integer,intent(IN)        :: nP,nCompQ,nContQ,nC,nD,startC,startD,nContC,nContD
Integer,intent(IN)        :: nCompC,nCompD,nPassQ,iOrder,nOrder,iDer,nDer,LUPRI,IPRINT
Real(realk),intent(INOUT) :: tuvCD(nP,nC,nD,nPassQ,nOrder)
Real(realk),intent(IN)    :: tuvQ(nContQ,nP,nDer,nCompQ)
!
Integer :: iPrimP,iAng,iContQ,iC,iD,iPassQ,iCompC,iCompD,iContC,iContD,iP
!
iContQ = 0
DO iPassQ=1,nPassQ
 DO iContD=1,nContD
  DO iContC=1,nContC
   iContQ = iContQ + 1
   iAng = 0
   iD = startD + (iContD-1)*nCompD - 1
   DO iCompD=1,nCompD
    iD = iD + 1
    iC = startC + (iContC-1)*nCompC - 1
    DO iCompC=1,nCompC
     iC = iC + 1
     iAng = iAng + 1
     DO iP= 1,nP
      TUVCD(iP,iC,iD,iPassQ,iOrder) = TUVQ(iContQ,iP,iDer,iAng)
     ENDDO !iP
    ENDDO !iContC
   ENDDO !iContD
  ENDDO !iPassQ
 ENDDO !iCompC
ENDDO !iCompD
!
END SUBROUTINE AddToTUVQ1

SUBROUTINE AddToTUVQ1seg(tuvCD,tuvQ,nP,nCompQ,nContQ,nC,nD,startC,startD,&
     &                   nContC,nContD,nCompC,nCompD,nPassQ,iOrder,nOrder,iDer,nDer,&
     &                   LUPRI,IPRINT)
implicit none
Integer,intent(IN)        :: nP,nCompQ,nContQ,nC,nD,startC,startD,nContC,nContD
Integer,intent(IN)        :: nCompC,nCompD,nPassQ,iOrder,nOrder,iDer,nDer,LUPRI,IPRINT
Real(realk),intent(INOUT) :: tuvCD(nP,nC,nD,nPassQ,nOrder)
Real(realk),intent(IN)    :: tuvQ(nContQ,nP,nDer,nCompQ)
!
Integer :: iPrimP,iAng,iContQ,iC,iD,iPassQ,iCompC,iCompD,iP
!
DO iPassQ=1,nPassQ
 iAng = 0
 DO iD=startD,startD-1+nCompD
    DO iC=startC,startC-1+nCompC
       iAng = iAng + 1
       DO iP= 1,nP
          TUVCD(iP,iC,iD,iPassQ,iOrder) = TUVQ(iPassQ,iP,iDer,iAng)
       ENDDO !iP
    ENDDO !iCompC
 ENDDO !iCompD
ENDDO !iPassQ

END SUBROUTINE AddToTUVQ1seg

SUBROUTINE AddToTUVQ1gen(tuvCD,tuvQ,nP,nCompQ,nContQ,nC,nD,startC,startD,&
     &                nContC,nContD,nCompC,nCompD,nPassQ,iOrder,nOrder,iDer,nDer,&
     &                permute,LUPRI,IPRINT)
implicit none
Integer,intent(IN)        :: nP,nCompQ,nContQ,nC,nD,startC,startD,nContC,nContD
Integer,intent(IN)        :: nCompC,nCompD,nPassQ,iOrder,nOrder,iDer,nDer,LUPRI,IPRINT
Real(realk),intent(INOUT) :: tuvCD(nP,nC,nD,nPassQ,nOrder)
Real(realk),intent(IN)    :: tuvQ(nContQ,nP,nDer,nCompQ)
Logical,intent(IN)        :: permute
!
Integer :: iPrimP,iAng,iContQ,iC,iD,iPassQ,iCompC,iCompD,iContC,iContD,iP
!
iContQ = 0
DO iPassQ=1,nPassQ
 DO iContD=1,nContD
  DO iContC=1,nContC
   iContQ = iContQ + 1
   iAng = 0
   iD = startD + (iContD-1)*nCompD - 1
   DO iCompD=1,nCompD
    iD = iD + 1
    iC = startC + (iContC-1)*nCompC - 1
    DO iCompC=1,nCompC
     iC = iC + 1
     iAng = iAng + 1
     DO iP= 1,nP
      TUVCD(iP,iC,iD,iPassQ,iOrder) = TUVQ(iContQ,iP,iDer,iAng)
     ENDDO !iP
     IF (permute) THEN
       DO iP= 1,nP
        TUVCD(iP,iD,iC,iPassQ,iOrder) = TUVQ(iContQ,iP,iDer,iAng)
       ENDDO !iP
     ENDIF
    ENDDO !iContC
   ENDDO !iContD
  ENDDO !iPassQ
 ENDDO !iCompC
ENDDO !iCompD
!
END SUBROUTINE AddToTUVQ1gen

SUBROUTINE AddToTUVQ1gen_ca(tuvCD,tuvQ,nP,nCompQ,nContQ,nC,nD,startC,startD,&
     &                      nContC,nContD,nCompC,nCompD,nPassQ,iOrder,nOrder,iDer,nDer,&
     &                      permute,LUPRI,IPRINT)
implicit none
Integer,intent(IN)        :: nP,nCompQ,nContQ,nC,nD,startC,startD,nContC,nContD
Integer,intent(IN)        :: nCompC,nCompD,nPassQ,iOrder,nOrder,iDer,nDer,LUPRI,IPRINT
Real(realk),intent(INOUT) :: tuvCD(nP,nC,nD,nPassQ,nOrder)
Real(realk),intent(IN)    :: tuvQ(nContQ,nP,nDer,nCompQ)
Logical,intent(IN)        :: permute
!
Integer :: iPrimP,iAng,iContQ,iC,iD,iPassQ,iCompC,iCompD,iContC,iContD,iP
!
iAng = 0
DO iCompD=1,nCompD
 DO iCompC=1,nCompC
  iAng = iAng + 1
  iContQ = 0
  DO iPassQ=1,nPassQ
   iD = startD + (iCompD-1)*nContD - 1
   DO iContD=1,nContD
    iD = iD + 1
    iC = startC + (iCompC-1)*nContC - 1
    DO iContC=1,nContC
     iC = iC + 1
     iContQ = iContQ + 1
     DO iP= 1,nP
      TUVCD(iP,iC,iD,iPassQ,iOrder) = TUVQ(iContQ,iP,iDer,iAng)
     ENDDO !iP
     IF (permute) THEN
       DO iP= 1,nP
        TUVCD(iP,iD,iC,iPassQ,iOrder) = TUVQ(iContQ,iP,iDer,iAng)
       ENDDO !iP
     ENDIF
    ENDDO !iPassQ
   ENDDO !iCompC
  ENDDO !iCompD
 ENDDO !iContC
ENDDO !iContD
!
END SUBROUTINE AddToTUVQ1gen_ca

SUBROUTINE PrintTuvQ(tuvCD,nP,nC,nD,nPassQ,iOrder,nOrder,LUPRI)
implicit none
Integer,intent(IN)     :: nP,nC,nD,nPassQ,iOrder,nOrder,LUPRI
Real(realk),intent(IN) :: tuvCD(nP,nC,nD,nPassQ,nOrder)
!
Integer :: iP,iC,iD,iPassQ

CALL LSHEADER(LUPRI,'TUVCD')
DO iPassQ=1,nPassQ
 DO iD=1,nD
  DO iC=1,nC
   WRITE(LUPRI,'(3X,A,I2,A,I2,A,I2)') 'iPassQ =',iPassQ,' iC =',iC,' iD =',iD
   WRITE(LUPRI,'(5F28.17)') (tuvCD(iP,iC,iD,iPassQ,iOrder),iP=1,nP)
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE PrintTuvQ

!> \brief Perform Spherical Transformation of the intermidiate integrals
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param Integral contains arrays to store intermidiates and final integrals
!> \param P the overlap distribution P or Q
!> \param iAngmom the index of angular momentum
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SphericalTransform(Integral,P,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Overlap)      :: P
Integer            :: LUPRI,IPRINT,iAngmom
!
!Real(realk),pointer :: Spherical(:,:)
Integer             :: lm1,lm2,ijk1,ijk2,iA1,iA2,ang1,ang2,lm,ijk
Integer             :: dim1,nOperatorComp
Logical             :: Sph1,Sph2,spherical
Real(realk),parameter :: D1=1.0E0_realk,D0=0.0E0_realk
Real(realk),dimension(:),pointer :: ptemp
Real(realk),pointer :: TMPWORK(:)
!
! Do not perform a spherical transformation if the E-coefficients are spherical
IF (P%sphericalEcoeff) RETURN

iA1  = P%indexAng1(iangmom)
iA2  = P%indexAng2(iangmom)
ang1 = P%orbital1%angmom(iA1)
ang2 = P%orbital2%angmom(iA2)
Sph1 = P%orbital1%spherical.AND.(ang1.GT. 1)
Sph2 = P%orbital2%spherical.AND.(ang2.GT. 1)
spherical = Sph1.OR.Sph2
IF (spherical) THEN
  nOperatorComp=1 
  IF(P%magderiv.EQ.1)nOperatorComp=3
  ijk1 = (ang1+1)*(ang1+2)/2 
  ijk2 = (ang2+1)*(ang2+2)/2
  CALL GET_IJK(ang1,ang2,lm1,lm2,lm,ijk,.TRUE.,0,P%single)
#ifdef VAR_LSDEBUGINT
IF(Integral%nPrim*Integral%nOrb*Integral%nGeoDeriv*ijk1*ijk2.GT.allocIntmaxTUVdim)THEN
   call lsquit('SphericalTransformation alloc error1',lupri)
ENDIF
IF(Integral%nPrim*Integral%nOrb*Integral%nGeoDeriv*lm1*lm2.GT.allocIntmaxTUVdim)THEN
   call lsquit('SphericalTransformation alloc error2',lupri)
ENDIF
#endif
!print*,'mem_workpointer_alloc(TMPWORK,ijk1*ijk2*lm1*lm2)'
call mem_workpointer_alloc(TMPWORK,ijk1*ijk2*lm1*lm2)
dim1 = Integral%nPrim*Integral%nOrb*Integral%ngeoDeriv*nOperatorComp
CALL SphericalTransformation(integral%TUV,TMPWORK,&
     & ijk1,ijk2,lm1,lm2,ang1,ang2,Integral%IN,Integral%OUT,&
     & dim1,lupri,iprint,1)
call mem_workpointer_dealloc(TMPWORK)
! Swap pointers IN and OUT
  ptemp => Integral%IN
  Integral%IN  => Integral%OUT
  Integral%OUT => ptemp
!=============================================

  Integral%nAng = lm  
ENDIF
END SUBROUTINE SphericalTransform

!> \brief Perform contraction with Ecoefficients 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param Integral contains arrays to store intermidiates and final integrals
!> \param P the overlap distribution P or Q
!> \param signP if Q overlap signP=-1 else signP=1
!> \param iAngmom the index of angular momentum
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE ContractEcoeffP(Integral,P,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(Integralitem)    :: Integral
TYPE(Overlap)         :: P
Integer               :: LUPRI,IPRINT,iAngmom
real(realk),parameter :: signP = 1.0_realk
!
Real(realk),pointer :: Ecoeffs(:)
Integer             :: l1,l2,ijk1,ijk2,ijk,ijkcart,iPrimP,nETUV,nPrim1,nPrim2
Integer             :: nOperatorComp,der
Real(realk),parameter  :: D1=1.0E0_realk,D0=0.0E0_realk
Real(realk),dimension(:),pointer :: ptemp
!
der = Integral%lhsGeoOrder

IF (P%type_hermite_single) THEN
   CALL SingleHermiteEcoeff(Integral,P,signP,iAngmom,P%nPasses,Integral%nAng,Integral%nOrb,LUPRI,IPRINT)
ELSE
   l1 = P%orbital1%angmom(P%indexAng1(iangmom))
   l2 = P%orbital2%angmom(P%indexAng2(iangmom))
   CALL GET_IJK(l1,l2,ijk1,ijk2,ijk,ijkcart,P%sphericalEcoeff,Integral%lhsGeoOrder,P%single)
   nOperatorComp = 1
   IF(P%magderiv.EQ.1)nOperatorComp = 3
   nETUV = ijk*Integral%nAng*Integral%nPrim*nOperatorComp
!   print*,'debug 564545465',nETUV
   call mem_workpointer_alloc(Ecoeffs,nETUV)
   CALL BuildEcoeffTensor(integral%TUV,P,signP,Ecoeffs,ijk,ijkcart,Integral%nAng,&
        &Integral%nPrim,iAngmom,der,P%nPasses,LUPRI,IPRINT,nOperatorComp)
#ifdef VAR_LSDEBUGINT
   IF(Integral%nPrim*Integral%nOrb*Integral%nAng.GT.allocIntmaxTUVdim)THEN
      call lsquit('ContractEcoeff1 alloc error1',lupri)
   ENDIF
   IF(Integral%nPrim*Integral%nOrb*ijk.GT.allocIntmaxTUVdim)THEN
      call lsquit('ContractEcoeff1 alloc error2',lupri)
   ENDIF
#endif
   CALL ContractEcoeff1(Integral%IN,Integral%OUT,Ecoeffs,Integral%nOrb,Integral%nAng,ijk*nOperatorComp,&
        &                 Integral%nPrim,P%segmented)      
   IF (P%segmented) Integral%nPrim = P%nPasses
   call mem_workpointer_dealloc(Ecoeffs)

   Integral%nAng = ijk*nOperatorComp
   ptemp => Integral%IN
   Integral%IN  => Integral%OUT
   Integral%OUT => ptemp
ENDIF

   IF (IPRINT.GT. 50) THEN
      CALL PrintTensor(Integral%IN,'ContractEcoeff      ',&
           &Integral%nPrim,Integral%nOrb,Integral%nAng,Lupri,&
           &'Prim  ','iOrb  ','ijk   ',3)
   ENDIF

END SUBROUTINE ContractEcoeffP

SUBROUTINE ContractEcoeff1(intIN,intOUT,Etuv,nOrb,nAng,ijk,nPrim,segmented)
implicit none
Real(realk),pointer :: intIN(:)
Real(realk),pointer :: intOUT(:)
Real(realk),pointer :: Etuv(:)
Integer,intent(IN)  :: nOrb,nAng,ijk,nPrim
Logical,intent(IN)  :: segmented

IF (segmented) THEN
  CALL ContractEcoeffseg(intIN,intOUT,Etuv,nOrb,nAng,ijk,nPrim)
ELSE
  CALL ContractEcoeffgen(intIN,intOUT,Etuv,nOrb,nAng,ijk,nPrim)
ENDIF
END SUBROUTINE ContractEcoeff1

!> \brief Contract Ecoefficients with Integrals \f$ OUT(ijk,nOrb,nPrim) = \sum_{ntuv} Ecoeffs(ijk,ntuv,nPrim)*IN(ntuv,nOrb,nPrim) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param IntegralIN input integrals to be contracted with Ecoeffs Ordering(ntuv,nOrb,nPrim)
!> \param IntegralOUT output of contraction with Ecoeffs Ordering(ntuv,nOrb,nPrim)
!> \param Ecoeffs cartesian or hermite Ecoefficients
!> \param startE start index in the Ecoefficient 
!> \param nOrb either nPrimP*ntuvP (since the P overlap has not been worked on yet) or nContQ*nCompOrbQ (since the Q overlap have been contracted)
!> \param nAng number of tuv components
!> \param ijk number of cartesian components or number of spherical components
!> \param nPrim number of primitive functions for the overlap being worked on
SUBROUTINE contractEcoeffgen(IntegralIN,IntegralOUT,Ecoeffs,nOrb,nTUV,nAng,nPrim)
implicit none           
Integer,intent(in)        :: nOrb,nTUV,nAng,nPrim
Real(realk),intent(in)    :: IntegralIN(nPrim,nOrb,nTUV)
Real(realk),intent(inout) :: IntegralOUT(nPrim,nOrb,nAng)
Real(realk),pointer       :: Ecoeffs(:)
!
Integer      :: endE

IF (nAng.EQ. 1) THEN
  CALL contractEcoeffGenSS(IntegralIN,IntegralOUT,Ecoeffs,nOrb,nPrim)
ELSE
  CALL contractEcoeffGenAng(IntegralIN,IntegralOUT,Ecoeffs,nOrb,nTUV,nAng,nPrim)
ENDIF
END SUBROUTINE contractEcoeffgen

!> \brief Contract Ecoefficients with Integrals \f$ OUT(ijk,nOrb,nPrim) = \sum_{ntuv} Ecoeffs(ijk,ntuv,nPrim)*IN(ntuv,nOrb,nPrim) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param IntegralIN input integrals to be contracted with Ecoeffs Ordering(ntuv,nOrb,nPrim)
!> \param IntegralOUT output of contraction with Ecoeffs Ordering(ntuv,nOrb,nPrim)
!> \param Ecoeffs cartesian or hermite Ecoefficients
!> \param startE start index in the Ecoefficient 
!> \param nOrb either nPrimP*ntuvP (since the P overlap has not been worked on yet) or nContQ*nCompOrbQ (since the Q overlap have been contracted)
!> \param nAng number of tuv components
!> \param ijk number of cartesian components or number of spherical components
!> \param nPrim number of primitive functions for the overlap being worked on
SUBROUTINE contractEcoeffseg(IntegralIN,IntegralOUT,Ecoeffs,nOrb,nTUV,nAng,nPrim)
implicit none           
Integer,intent(in)        :: nOrb,nTUV,nAng,nPrim
Real(realk),intent(in)    :: IntegralIN(nPrim,nOrb,nTUV)
Real(realk),intent(inout) :: IntegralOUT(nOrb,nAng)
Real(realk),pointer       :: Ecoeffs(:)
!
Integer      :: endE

IF (nAng.EQ. 1) THEN
  CALL contractEcoeffSegSS(IntegralIN,IntegralOUT,Ecoeffs,nOrb,nPrim)
ELSE
  IF (nPrim.GT. 1) THEN
    CALL contractEcoeffSegAng(IntegralIN,IntegralOUT,Ecoeffs,nOrb,nTUV,nAng,nPrim)
  ELSE
    CALL contractEcoeffSegAng1(IntegralIN,IntegralOUT,Ecoeffs,nOrb,nTUV,nAng)
  ENDIF
ENDIF
END SUBROUTINE contractEcoeffseg

SUBROUTINE contractEcoeffGenang(IntegralIN,IntegralOUT,Etuv,nOrb,nTUV,nAng,nPrim)
implicit none
Integer,intent(IN) :: nOrb,nTUV,nAng,nPrim
Real(realk),intent(IN)  :: IntegralIN(nPrim,nOrb,nTUV)
Real(realk),intent(IN)  :: Etuv(nPrim,nTUV,nAng)
Real(realk),intent(OUT) :: IntegralOUT(nPrim,nOrb,nAng)
!
Integer      :: iPrimP,iTUV,iAng,iOrb
Real(realk) :: tmp
!
DO iAng= 1,nAng
 DO iOrb = 1,nOrb
  DO iPrimP = 1,nPrim
    tmp = IntegralIN(iPrimP,iOrb,1)*Etuv(iPrimP,1,iAng)
    DO iTUV=2,nTUV
      tmp=tmp+IntegralIN(iPrimP,iOrb,iTUV)*Etuv(iPrimP,iTUV,iAng)
    ENDDO
    IntegralOUT(iPrimP,iOrb,iAng) = tmp
  ENDDO
 ENDDO
ENDDO
!
END SUBROUTINE contractEcoeffGenang


SUBROUTINE contractEcoeffSegang(IntegralIN,IntegralOUT,Etuv,nOrb,nTUV,nAng,nPrim)
implicit none
Integer,intent(IN) :: nOrb,nTUV,nAng,nPrim
Real(realk),intent(IN)  :: IntegralIN(nPrim,nOrb,nTUV)
Real(realk),intent(IN)  :: Etuv(nPrim,nTUV,nAng)
Real(realk),intent(OUT) :: IntegralOUT(nOrb,nAng)
!
Integer      :: iPrimP,iTUV,iAng,iOrb
Real(realk) :: tmp
!
DO iAng= 1,nAng
 DO iOrb = 1,nOrb
  tmp = 0E0_realk
  DO iTUV=1,nTUV
   DO iPrimP = 1,nPrim
    tmp=tmp+IntegralIN(iPrimP,iOrb,iTUV)*Etuv(iPrimP,iTUV,iAng)
   ENDDO
  ENDDO
  IntegralOUT(iOrb,iAng) = tmp
 ENDDO
ENDDO
!
END SUBROUTINE contractEcoeffSegang

SUBROUTINE contractEcoeffSegang1(IntegralIN,IntegralOUT,Etuv,nOrb,nTUV,nAng)
implicit none
Integer,intent(IN) :: nOrb,nTUV,nAng
Real(realk),intent(IN)  :: IntegralIN(nOrb,nTUV)
Real(realk),intent(IN)  :: Etuv(nTUV,nAng)
Real(realk),intent(OUT) :: IntegralOUT(nOrb,nAng)
!
Integer      :: iTUV,iAng,iOrb
Real(realk) :: tmp
real(realk),parameter :: zero=1.0E-16_realk,D0=0E0_realk
!
!DO iAng= 1,nAng
! DO iOrb = 1,nOrb
!  tmp = IntegralIN(iOrb,1)*Etuv(1,iAng)
!  DO iTUV=2,nTUV
!    tmp=tmp+IntegralIN(iOrb,iTUV)*Etuv(iTUV,iAng)
!  ENDDO
!  IntegralOUT(iOrb,iAng) = tmp
! ENDDO
!ENDDO
DO iAng= 1,nAng
   DO iOrb = 1,nOrb
      IntegralOUT(iOrb,iAng) = D0
   ENDDO
   DO iTUV=1,nTUV
      tmp = Etuv(iTUV,iAng)
      IF(ABS(tmp).GT.zero)THEN
         DO iOrb = 1,nOrb
            IntegralOUT(iOrb,iAng) = IntegralOUT(iOrb,iAng) + IntegralIN(iOrb,iTUV)*tmp
         ENDDO
      ENDIF
   ENDDO
ENDDO

END SUBROUTINE contractEcoeffSegang1


!> \brief special case of Contract Ecoefficients with Integrals when ijk=1 \f$ OUT(1,nOrb,nPrim) = \sum_{ntuv} Ecoeffs(1,1,nPrim)*IN(1,nOrb,nPrim) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param IntegralIN input integrals to be contracted with Ecoeffs Ordering(nOrb,nPrim)
!> \param IntegralOUT output of contraction with Ecoeffs Ordering(nOrb,nPrim)
!> \param Ecoeffs cartesian or hermite Ecoefficients
!> \param nOrb either nPrimP*ntuvP (since the P overlap has not been worked on yet) or nContQ*nCompOrbQ (since the Q overlap have been contracted)
!> \param nPrim number of primitive functions for the overlap being worked on
SUBROUTINE contractEcoeffGenss(IntegralIN,IntegralOUT,Ecoeffs,nOrb,nPrim)
 implicit none
 Integer,intent(in)     :: nOrb,nPrim
 Real(realk),intent(in),   dimension(nPrim,nOrb) :: IntegralIN
 Real(realk),intent(inout),dimension(nPrim,nOrb) :: IntegralOUT
 Real(realk),intent(in),dimension(nPrim)         :: Ecoeffs
 !
 Integer     :: iPrim,iOrb

 DO iOrb=1,nOrb
    DO iPrim = 1,nPrim
       IntegralOUT(iPrim,iOrb) = Ecoeffs(iPrim)*IntegralIN(iPrim,iOrb)
    ENDDO
 ENDDO

END SUBROUTINE contractEcoeffGenss

!> \brief special case of Contract Ecoefficients with Integrals when ijk=1 \f$ OUT(1,nOrb,nPrim) = \sum_{ntuv} Ecoeffs(1,1,nPrim)*IN(1,nOrb,nPrim) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param IntegralIN input integrals to be contracted with Ecoeffs Ordering(nOrb,nPrim)
!> \param IntegralOUT output of contraction with Ecoeffs Ordering(nOrb,nPrim)
!> \param Ecoeffs cartesian or hermite Ecoefficients
!> \param nOrb either nPrimP*ntuvP (since the P overlap has not been worked on yet) or nContQ*nCompOrbQ (since the Q overlap have been contracted)
!> \param nPrim number of primitive functions for the overlap being worked on
SUBROUTINE contractEcoeffSegss(IntegralIN,IntegralOUT,Ecoeffs,nOrb,nPrim)
 implicit none
 Integer,intent(in)     :: nOrb,nPrim
 Real(realk),intent(in),dimension(nPrim,nOrb) :: IntegralIN
 Real(realk),intent(inout),dimension(nOrb)    :: IntegralOUT
 Real(realk),intent(in),dimension(nPrim)      :: Ecoeffs
 !
 Integer     :: iPrim,iOrb
 Real(realk) :: tmp

 DO iOrb=1,nOrb
    tmp = Ecoeffs(1)*IntegralIN(1,iOrb)
    DO iPrim = 2,nPrim
       tmp = tmp + Ecoeffs(iPrim)*IntegralIN(iPrim,iOrb)
    ENDDO
    IntegralOUT(iOrb) = tmp
 ENDDO

END SUBROUTINE contractEcoeffSegss

!> \brief wrapper routine for contract with Ecoefficients, special case for single hermite ODbatch
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Herint the hermite integrals
!> \param P the overlap distribution P or Q
!> \param signP the sign (1 for P and -1 for Q)
!> \param iAngmom the angular moment index
!> \param nPasses the number of passes
!> \param nAng the number of angular components
!> \param nOrb the number of orbital components
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SingleHermiteEcoeff(Integral,P,signP,iAngmom,nPasses,nAng,nOrb,LUPRI,IPRINT)
implicit none
Type(IntegralItem) :: Integral
TYPE(Overlap) :: P
Integer       :: iAngmom,nPasses,nAng,nOrb,LUPRI,IPRINT
Real(realk)   :: signP
!
Integer               :: nPrim1,nPrim2
Integer               :: j1,j2,i1,i2
Real(realk)           :: sign
Real(realk),parameter :: D1 = 1.0E0_realk
Real(realk),pointer   :: ptemp(:)
!
i1 = P%indexAng1(iAngmom)
i2 = P%indexAng2(iAngmom)
j1 = P%orbital1%angmom(i1)
j2 = P%orbital2%angmom(i2)
IF(j1.EQ. 0 .AND. j2 .EQ. 0)THEN
   sign = D1
ELSE
   sign = signP**(j1+j2)
ENDIF
nPrim1 = P%orbital1%nPrimitives
nPrim2 = P%orbital2%nPrimitives
IF (P%segmented) THEN
! Multiply in preexponential factor and contract to segmented basis
#ifdef VAR_LSDEBUGINT
IF(nPrim1*nPrim2*nPasses*nAng*nOrb.GT.allocIntmaxTUVdim)THEN
   call lsquit('SingleHermiteEcoeff1seg alloc error1',lupri)
ENDIF
IF(nPasses*nAng*nOrb.GT.allocIntmaxTUVdim)THEN
   call lsquit('SingleHermiteEcoeff1seg alloc error2',lupri)
ENDIF
#endif
  CALL SingleHermiteEcoeff1seg(Integral%OUT,Integral%IN,P,sign,j1,j2,&
       & nPrim1,nPrim2,nPasses,nAng,nOrb,&
       & P%preexpfac,P%nPrimAlloc,P%nAngAlloc,iAngmom,LUPRI,IPRINT)
! Swap pointers
  ptemp => Integral%OUT
  Integral%OUT => Integral%IN
  Integral%IN  => ptemp
! Number of contracted equal to one, therefore
  Integral%nPrim = P%nPasses
ELSE
! Multiply in preexponential factor
#ifdef VAR_LSDEBUGINT
IF(nPrim1*nPrim2*nPasses*nAng*nOrb.GT.allocIntmaxTUVdim)THEN
   call lsquit('SingleHermiteEcoeff1gen alloc error1',lupri)
ENDIF
#endif
  CALL SingleHermiteEcoeff1gen(Integral%in,P,sign,j1,j2,nPrim1,nPrim2,nPasses,nAng,nOrb,&
       & P%preexpfac,P%nPrimAlloc,P%nAngAlloc,iAngmom,LUPRI,IPRINT)
ENDIF
END SUBROUTINE SingleHermiteEcoeff

!!> \brief contract with Ecoefficients, special case for single hermite ODbatch
!!> \author S. Reine and T. Kjaergaard
!!> \date 2009
!!> \param Herint the hermite integrals
!!> \param P the overlap distribution P or Q
!!> \param signP the sign (1 for P and -1 for Q)
!!> \param j1 angular moment for orbital 1
!!> \param j2 angular moment for orbital 2
!!> \param nPrim1 number of primitives for orbital 1
!!> \param nPrim2 number of primitives for orbital 2
!!> \param nPasses the number of passes
!!> \param nAng the number of angular components
!!> \param nOrb the number of orbital components
!!> \param LUPRI the logical unit number for the output file
!!> \param IPRINT the printlevel, determining how much output should be generated
!SUBROUTINE SingleHermiteEcoeff1(HerInt,P,signP,j1,j2,nPrim1,nPrim2,nPasses,nAng,nOrb,PreExpFac,LUPRI,IPRINT)
!implicit none
!Real(realk)   :: HerInt(nPrim1*nPrim2,nPasses*nAng*nOrb)
!TYPE(Overlap) :: P
!Integer       :: nAng,nOrb,LUPRI,IPRINT
!Integer       :: nPrim1,nPrim2,j1,j2,nPasses
!Real(realk)   :: signP,preexpfac(nPrim1*nPrim2*nPasses)
!!
!Real(realk)   :: pref2,pref12,D1=1.0E0_realk,D2=2.0E0_realk
!Integer       :: iPrim1,iPrim2,iAngOrb,iPasses,iPrim12
!
!IF(P%contractBasis)then
!   iPrim12 = 0
!   DO iPrim2=1,nPrim2
!      pref2=signP/(D2*P%orbital2%exponents(iPrim2))**j2
!      DO iPrim1=1,nPrim1
!         iPrim12 = iPrim12 + 1
!         pref12=D1/(D2*P%orbital1%exponents(iPrim1))**j1*pref2
!         DO iAngOrb=1,nAng*nOrb*nPasses
!            HerInt(iPrim12,iAngOrb) = HerInt(iPrim12,iAngOrb) * pref12
!         ENDDO
!      ENDDO
!   ENDDO
!ELSE
!!PS_int
!   iPrim12 = 0
!   DO iPrim2=1,nPrim2
!      pref2=signP/(D2*P%orbital2%exponents(iPrim2))**j2
!      DO iPrim1=1,nPrim1
!         iPrim12 = iPrim12 + 1
!         pref12=D1/(D2*P%orbital1%exponents(iPrim1))**j1*pref2*preexpfac(iPrim12)
!         DO iAngOrb=1,nAng*nOrb*nPasses
!            HerInt(iPrim12,iAngOrb) = HerInt(iPrim12,iAngOrb) * pref12
!         ENDDO
!      ENDDO
!   ENDDO
!ENDIF
!
!END SUBROUTINE SingleHermiteEcoeff1

!> \brief Print after basisset contraction
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param Contracted multiarrray to be printet
!> \param nAng first dimension
!> \param nOrb second dimension
!> \param nPrim third dimension
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE printContracted(Contracted,nAng,nOrb,nPrim,LUPRI,IPRINT)
 Integer,intent(in)     :: nAng,nOrb,nPrim
 Real(realk),intent(in) :: Contracted(nPrim,nOrb,nAng)
 !
 Integer :: iAng,iOrb,iPrim
 IF (IPRINT.GT. 50) THEN
    WRITE(LUPRI,'(5X,A)') '**********************************************************'
    WRITE(LUPRI,'(5X,A)') '***                   Contracted'
    WRITE(LUPRI,'(5X,A)') '**********************************************************'
    DO iAng=1,nAng
       DO iOrb=1,nOrb
          WRITE(LUPRI,'(5X,A,I3,A,I3)') 'iAngP =',iAng,' iOrb =',iOrb
          WRITE(LUPRI,'(2X,5ES13.6/(2X,5ES13.6))') (Contracted(iPrim,iOrb,iAng),iPrim=1,nPrim)
       ENDDO
    ENDDO
 ENDIF
END SUBROUTINE printContracted

!> \brief wrapper routine to perform hermite integrals with FTUVs
!> \author S. Reine
!> \date 2010
!>
!> \param Integral contains arrays to store intermidiates and final integrals
!> \param P the left hand side overlap distribution 
!> \param Q the right hand side overlap distribution 
!> \param ndmat number of density matrices
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE ContractFTUV(Integral,P,Q,ndmat,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Overlap)      :: P,Q
Integer            :: LUPRI,IPRINT,ndmat,tuvOffP,startAngP,tuvOffQ,startAngQ,startAngPQ,tuvOffPQ
!
startAngP  = P%startAngmom
tuvOffP    = startAngP*(startAngP+1)*(startAngP+2)/6
startAngQ  = Q%startAngmom
tuvOffQ    = startAngQ*(startAngQ+1)*(startAngQ+2)/6
startAngPQ = startAngP + startAngQ
tuvOffPQ   = startAngPQ*(startAngPQ+1)*(startAngPQ+2)/6

#ifdef VAR_LSDEBUGINT
IF(P%nPrimitives*P%nTUV*ndmat.GT.allocIntmaxTUVdim)THEN
   call lsquit('ContractFTUVPA alloc error1',lupri)
ENDIF
IF(Integral%ntuv*Integral%nPrim.GT.allocIntmaxTUVdim)THEN
   call lsquit('ContractFTUVPA alloc error2',lupri)
ENDIF
IF(size(Q%FTUV).NE.Q%nPrimAlloc*Q%nTUVAlloc*ndmat)THEN
   call lsquit('ContractFTUVPA alloc error3',lupri)
ENDIF

#endif
CALL ContractFTUVPA(Integral%tuvQ,Integral%RTUV,Q%FTUV,Integral%TUV%TUVindex,&
     &              Integral%TUV%Tindex,Integral%TUV%Uindex,Integral%TUV%Vindex,&
     &              Integral%ntuv,Integral%nPrim,P%nTUV,Q%nTUV,tuvOffP,tuvOffQ,tuvOffPQ,&
     &              P%nPrimitives,Q%nPrimitives,ndmat,Q%nPrimAlloc,Q%nTUVAlloc,lupri)
Integral%nAng  = 1
Integral%nPrim = ndmat
Integral%nOrb  = P%nPrimitives*P%nTUV
!CALL PrintTuvQ_old(Integral%integrals,P%nTUV,P%nPrimitives,1,ndmat,LUPRI,IPRINT)
#ifdef VAR_LSDEBUGINT
IF (IPRINT.GT. 50) THEN
   CALL PrintTensor(Integral%tuvQ,' InnerFTUV          ',P%nPrimitives,P%nTUV,ndmat,Lupri,&
        & 'iPrim ','iTUV  ','idmat ',1)
ENDIF
#endif

END SUBROUTINE ContractFTUV

!> \brief Perform hermite integrals with FTUVs ordering (nPrimPQ,nTUVPQ)
!> \author S. Reine
!> \date 2010
!>
!> \param TUV the output intermediate integral
!> \param WTUV the primitive hermite integrals 
!> \param FTUV the FTUV batches 
!> \param TUVindex for a given t,u,v values it gives the TUV index (t,u,v)=(000,100,010,001,.. etc) means TUV = (1,2,3,4,..)
!> \param Tindex for a given TUV index it gives the T index
!> \param Uindex for a given TUV index it gives the U index
!> \param Vindex for a given TUV index it gives the V index
!> \param ntuvPQ number of angular components for PQ
!> \param nPrimPQ total number of primitive
!> \param ntuvP number of angular components for P
!> \param ntuvQ number of angular components for Q
!> \param tuvoffP offset for tuvindexing for overlap P
!> \param tuvoffQ offset for tuvindexing for overlap Q
!> \param tuvoffPQ offset for tuvindexing for overlap PQ
!> \param nPrimP number of primitive for overlap P
!> \param nPrimQ number of primitive for overlap Q
!> \param nDmat the number of density matrices
!> \param LUPRI the logical unit number for the output file
SUBROUTINE ContractFTUVPA(TUV,WTUV,FTUV,TUVindex,Tindex,Uindex,Vindex,&
     &                    ntuvPQ,nPrimPQ,ntuvP,ntuvQ,tuvOffP,tuvOffQ,tuvOffPQ,&
     &                    nPrimP,nPrimQ,nDmat,nPrimAlloc,nTUVAlloc,lupri)
implicit none
Integer             :: ntuvP,ntuvQ,nPrimP,nPrimQ,ntuvPQ,nPrimPQ,nDmat
Integer             :: tuvOffP,tuvOffQ,tuvOffPQ,nPrimAlloc,nTUVAlloc
Real(realk)         :: TUV(nPrimP,ntuvP,nDmat)
Real(realk)         :: FTUV(nPrimAlloc*nTUVAlloc*ndmat)
Real(realk)         :: WTUV(nPrimPQ,nTUVPQ)
Integer,pointer     :: TUVindex(:,:,:),Tindex(:),Uindex(:),Vindex(:)
!
Integer     :: idmat,tuvP,tP,uP,vP,tuvQ,tQ,uQ,vQ,jQ,tuvPQ
Integer     :: iPrimP,iPrimQ,iPrimPQ,lupri,offset
Real(realk) :: fsum
DO tuvP=1,ntuvP
 tP = Tindex(tuvP+tuvOffP)
 uP = Uindex(tuvP+tuvOffP)
 vP = Vindex(tuvP+tuvOffP)
 DO tuvQ=1,ntuvQ
  tQ = Tindex(tuvQ+tuvOffQ)
  uQ = Uindex(tuvQ+tuvOffQ)
  vQ = Vindex(tuvQ+tuvOffQ)
  jQ = tQ+uQ+vQ
  tuvPQ = TUVindex(tP+tQ,uP+uQ,vP+vQ)-tuvOffPQ
  DO idmat = 1,ndmat
   iPrimPQ = 1
   offset = (tuvQ-1)*nPrimAlloc+(idmat-1)*nPrimAlloc*nTUVAlloc
   DO iPrimP=1,nPrimP
    fsum = WTUV(iPrimPQ,tuvPQ)*FTUV(1+offset)
    iPrimPQ = iPrimPQ+1
    DO iPrimQ=2,nPrimQ
     fsum = fsum + WTUV(iPrimPQ,tuvPQ)*FTUV(iPrimQ+offset)
     iPrimPQ = iPrimPQ+1
    ENDDO
    TUV(iPrimP,tuvP,idmat) = TUV(iPrimP,tuvP,idmat) + fsum
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE ContractFTUVPA

!> \brief determine the start and end index of the RHS loop 
!> \author S. Reine
!> \date 2010
!> \param INPUT the integral input specifications, contains all info about what to do
!> \param OD_RHS the ODbatches belonging to the Reft hand side
!> \param ILHS the left hand side 
!> \param start_RHS the start index of the RHS loop
!> \param end_RHS the end index of the RHS loop
SUBROUTINE Determine_RHS_loop(INPUT,OD_RHS,ILHS,Start_RHS,End_RHS)
implicit none
TYPE(INTEGRALINPUT)  :: INPUT
TYPE(ODITEM)         :: OD_RHS
Integer              :: ILHS,Start_RHS,End_RHS
! Screening-integral loop
  IF (INPUT%CS_int.OR.INPUT%PS_int) THEN
    Start_RHS = ILHS
    End_RHS   = ILHS
! Triangular loop
  ELSE IF (INPUT%sameODs) THEN
    Start_RHS = 1
    End_RHS   = ILHS
! Regular loop
  ELSE
    Start_RHS = 1
    End_RHS   = OD_RHS%nbatches
  ENDIF
END SUBROUTINE Determine_RHS_loop

!!> \brief clear the integrals
!!> \author S. Reine and T. Kjaergaard
!!> \date 2010
!!> \param Integral contains arrays to store intermidiates and final integrals
!SUBROUTINE CLEAR_integral(INTEGRAL)
!TYPE(Integralitem)      :: integral
! 
!DEALLOC(INTEGRAL%Wtuv)
!NULLIFY(INTEGRAL%Wtuv)
! 
!END SUBROUTINE

!> \brief Subroutine that returns A, the transpose of B
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param A the output matrix, to be transpose of B
!> \param B the input matrix
!> \param n1 dimension 1 of the output matrix
!> \param n2 dimension 2 of the output matrix
SUBROUTINE transposition(A,B,n1,n2)
implicit none
integer     :: n1,n2
Real(realk) :: A(n1,n2)
Real(realk) :: B(n2,n1)
A = transpose(B)
END SUBROUTINE transposition

!> \brief determine the maximum distance
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param MAXDIST the maximum distance
!> \param distance the input distance
!> \param nprim the size of distance
SUBROUTINE MAXDISTANCE(MAXDIST,distance,nPrim)
IMPLICIT NONE
INTEGER      :: nPrim,I
REAL(REALK)  :: distance(nPrim),MAXDIST

MAXDIST = sum ( abs ( distance(1:1+(nPrim-1)*1:1) ) )

END SUBROUTINE MAXDISTANCE

!> \brief sort routine, sorts largest element first 
!> \author T. Kjaergaard
!> \date 2010
!> \param A the array according to which the INDEXES should be sorted 
!> \param INDEXES an array of indexes which should be reordered
RECURSIVE SUBROUTINE Linksort(a,INDEXES,ndim,THRLOG)
IMPLICIT NONE
INTEGER :: ndim
integer(kind=short), INTENT(INOUT) :: a(ndim)
INTEGER :: INDEXES(ndim)
integer(kind=short) :: THRLOG
!
integer :: n(THRLOG-2:10),offset(THRLOG-2:10)
integer :: I
integer(kind=short) :: tmpa(ndim)
INTEGER :: tmpINDEXES(ndim),nmin,nmax
nmin=THRLOG-2
nmax=10 !could be set using maxgabelm 
Do I = nmin,nmax
   n(I)=0
ENDDO
Do I = 1, ndim
   n(a(I)) = n(a(I))+1
End Do
offset(nmax)=0
Do I = nmax-1,nmin,-1
   offset(I)=n(I+1)+offset(I+1)
ENDDO
Do I = nmin,nmax
   n(I)=0
ENDDO
Do I = 1, ndim
   n(a(I)) = n(a(I))+1
!   tmpa(offset(a(I))+n(a(I)))=a(I)
   tmpINDEXES(offset(a(I))+n(a(I)))=INDEXES(I)
End Do
Do I = 1, ndim
!   a(I=)tmpa(I)
   INDEXES(I)=tmpINDEXES(I)
End Do

END SUBROUTINE LINKSORT

!> \brief alloc the integral item
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param PQ contains integrand info, to be allocd, using values in alloc
!> \param Integral contains arrays to store intermidiates and final integrals, to be allocd
!> \param input the integral input specifications, contains all info about what to do
!> \param Alloc the allocitem, which contain info about maximum values required for allocation
!> \param maxpasses the maximum number of passes
!> \param ndmat number of density matrix
SUBROUTINE allocIntegrals(PQ,Integral,maxPrim,maxTUVdim,maxOrb,lupri)
  use memory_handling
implicit none
Integer :: lupri
TYPE(Integrand)    :: PQ
TYPE(Integralitem) :: Integral
Integer :: maxTUVdim,maxPrim,maxOrb
logical :: orderAngPrim
integer (kind=long) :: nsize1,nsize2

call mem_alloc(PQ%distance,maxPrim,3)
call mem_alloc(PQ%squaredDistance,maxPrim)
Call Mem_alloc(PQ%exponents,maxPrim)
Call Mem_alloc(PQ%reducedExponents,maxPrim)
Call Mem_alloc(PQ%integralPrefactor,maxPrim)

#ifdef VAR_LSDEBUGINT
call set_allocIntmaxTUVdim(maxprim,maxTUVdim,maxOrb)
#endif

CALL MEM_ALLOC(Integral%IN,maxTUVdim)
CALL MEM_ALLOC(Integral%OUT,maxTUVdim)
CALL MEM_ALLOC(Integral%RTUV,maxTUVdim)
CALL MEM_ALLOC(Integral%tuvQ,maxTUVdim)
CALL MEM_ALLOC(Integral%integralsABCD,maxTUVdim)

nsize1 = mem_realsize*size(PQ%distance)  + mem_realsize*size(PQ%squaredDistance) + &
        & mem_realsize*size(PQ%exponents) + mem_realsize*size(PQ%reducedExponents)+ &
        & mem_realsize*size(PQ%integralPrefactor) 
call mem_allocated_mem_integrand(nsize1)

nsize2 = mem_realsize*size(Integral%IN)  + mem_realsize*size(Integral%OUT) + &
        & mem_realsize*size(Integral%RTUV) + mem_realsize*size(Integral%RTUV) + &
        & mem_realsize*size(Integral%integralsABCD)
call mem_allocated_mem_integralitem(nsize2)

Integral%startDerivativeIndex = 0
Integral%nDerivComp = 1
Integral%dohodi = .FALSE.

END SUBROUTINE allocIntegrals

!> \brief Wrapper routine that first determine the size, and then alloc the integral item
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param PQ contains integrand info, to be allocd, using values in alloc
!> \param Integral contains arrays to store intermidiates and final integrals, to be allocd
!> \param input the integral input specifications, contains all info about what to do
!> \param Alloc the allocitem, which contain info about maximum values required for allocation
!> \param iAlloc the ODtype index
!> \param npassTypes the number of different passtypes
!> \param maxpassesfortype and output which determines the number of passes we collect for each type
!> \param ndmat number of density matrix
!> \param dopasses a logical which specifies for each passtype, if we should collect passes 
!> \param nOverlapofPasstype the number of overlaps of each passtype
!> \param lupri logical unit number
SUBROUTINE allocIntegralsWRAP(PQ,Integral,Input,Alloc,iAlloc,npassTypes,MaxPassesForType,&
     & ndmat,dopasses,nOverlapOfPassType,lupri)
use memory_handling
implicit none
Integer,intent(in)   :: nPassTypes,ndmat,lupri,iAlloc
Integer,intent(in)   :: nOverlapOfPassType(nPassTypes)
integer,intent(inout):: MaxPassesForType(nPassTypes)
TYPE(Integrand)      :: PQ
TYPE(Integralitem)   :: Integral
TYPE(Integralinput)  :: Input
TYPE(Allocitem)      :: Alloc
logical,intent(in) :: dopasses
!
Integer :: maxTUVdim,maxPrim
Integer,parameter :: cachesize=256*1000 !=256000 elements of 8 bytes = 2.048 Mb
!Integer,parameter :: cachesize=64*1000 !=64000 elements of 8 bytes = 512 kb
Integer :: MAXMEM = cachesize 
Integer :: passlist(5),iPassType,maxTUVdimTEMP,maxPrimTEMP,MEMUSAGE,ipass
integer :: trialnpass,maxOrb,maxOrbTEMP,maxTUVdimTEMP2,maxPrimTEMP2,maxOrbTEMP2
integer :: maxTUVdimTEMP3,maxPrimTEMP3,maxOrbTEMP3,MEMUSAGE3,CARTORDER,nEFG

call determineMaxPassesForType(Alloc,iAlloc,npassTypes,MaxPassesForType,&
     & ndmat,dopasses,nOverlapOfPassType,lupri)

maxTUVdim=0
maxPrim=0
maxOrb=0
do iPassType=1,nPassTypes
   maxTUVdimTEMP = Alloc%maxPrimTUVLHSA(iAlloc)*Alloc%maxPrimTUVRHSA(iPassType)*MaxPassesForType(iPassType)
   maxTUVdim = MAX(maxTUVdim,maxTUVdimTEMP)
   maxPrimTEMP = Alloc%maxPrimLHSA(iAlloc)*Alloc%maxPrimRHSA(iPassType)*MaxPassesForType(iPassType)
   maxPrim   = MAX(maxPRim,maxPrimTEMP)
   maxOrbTEMP = Alloc%maxTotOrbLHSA(iAlloc)*ndmat*Alloc%maxTotOrbRHSA(iPassType)*MaxPassesForType(iPassType)
   maxOrb = MAX(maxorb,maxOrbTEMP)
enddo
CARTORDER = MAX(INPUT%GEODERIVORDER,INPUT%MAGDERIVORDER,INPUT%MMORDER,INPUT%CMORDER)
CARTORDER=CARTORDER+Input%PropmaxM+Input%PropmaxD
IF(Input%PropRequireBoys.GT.-1)CARTORDER=CARTORDER+Input%PropRequireBoys
nEFG = (CARTORDER+1)*(CARTORDER+2)*(CARTORDER+3)/6
maxTUVdim = maxTUVdim*nEFG
IF(INPUT%PROPTYPE.EQ.42)THEN
   maxTUVdim = maxTUVdim*9
ENDIF
IF(INPUT%PROPTYPE.EQ.43)THEN
   maxTUVdim = maxTUVdim*18
ENDIF

call allocIntegrals(PQ,Integral,maxPrim,maxTUVdim,maxOrb,lupri)

END SUBROUTINE allocIntegralsWRAP

subroutine determineMaxPassesForType(Alloc,iAlloc,npassTypes,MaxPassesForType,&
     & ndmat,dopasses,nOverlapOfPassType,lupri)
use memory_handling
implicit none
Integer,intent(in)   :: nPassTypes,ndmat,lupri,iAlloc
Integer,intent(in)   :: nOverlapOfPassType(nPassTypes)
integer,intent(inout):: MaxPassesForType(nPassTypes)
TYPE(Allocitem)      :: Alloc
logical,intent(in) :: dopasses
!
Integer :: maxTUVdim,maxPrim
Integer,parameter :: cachesize=256*1000 !=256000 elements of 8 bytes = 2.048 Mb
!Integer,parameter :: cachesize=64*1000 !=64000 elements of 8 bytes = 512 kb
Integer :: MAXMEM = cachesize 
Integer :: passlist(5),iPassType,maxTUVdimTEMP,maxPrimTEMP,MEMUSAGE,ipass
integer :: trialnpass,maxOrb,maxOrbTEMP,maxTUVdimTEMP2,maxPrimTEMP2,maxOrbTEMP2
integer :: maxTUVdimTEMP3,maxPrimTEMP3,maxOrbTEMP3,MEMUSAGE3,CARTORDER,nEFG

passlist(1) = 4 
passlist(2) = 8 
passlist(3) = 16 
passlist(4) = 32 
passlist(5) = 64 
!passlist(6) = 128

maxTUVdim=0
maxPrim=0
maxOrb=0
if(dopasses)then
   do iPassType=1,nPassTypes
      MaxPassesForType(iPassType) = 1
   enddo
   ! we would like to collect overlap together in passes, but we 
   ! do not want to use to much memory. 
   !
   !we find the number of passes which is as large as possible without being larger than MAXMEM for any RHSpassType
   do iPassType=1,nPassTypes
      ! The memory usage for a single overlap of this type
      maxTUVdimTEMP = Alloc%maxPrimTUVLHSA(iAlloc)*Alloc%maxPrimTUVRHSA(iPassType)
      maxPrimTEMP   = Alloc%maxPrimLHSA(iAlloc)*Alloc%maxPrimRHSA(iPassType)
      maxOrbTEMP = ndmat*Alloc%maxTotOrbLHSA(iAlloc)*Alloc%maxTotOrbRHSA(iPassType)
      MEMUSAGE = 7*maxPrimTEMP + 4*maxTUVdimTEMP + maxOrbTEMP
      ! we loop over possible number of collected number of ODtypes 
      ! in order to asses the memory usage for collected overlaps of this type
      do iPass=1,size(passlist)
         trialnpass = passlist(ipass)
         !the collected number of overlaps must not exceed the actual number of overlaps of this type 
         IF(trialnpass.LE.nOverlapOfPassType(iPassType))THEN
            !the memory usage of the collected number of overlaps must not exceed maxmem
            IF(MEMUSAGE*trialnpass .LT. MAXMEM)THEN
               ! as MEMUSAGE is determined from the maximum value of
               ! maxTUVdim,maxPrim and maxOrb
               ! these cannot exceed maxmem either 
               maxTUVdimTEMP2 = maxTUVdimTEMP*trialnpass
               maxPrimTEMP2 = maxPrimTEMP*trialnpass
               maxOrbTEMP2 = maxORbTEMP*trialnpass
               maxTUVdimTEMP3 = MAX(maxTUVdimTEMP2,maxTUVdim)
               maxPrimTEMP3 = MAX(maxPrimTEMP2,maxPrim) 
               maxOrbTEMP3 = MAX(maxOrbTEMP2,maxOrb) 
               MEMUSAGE3 = 7*maxPrimTEMP3 + 4*maxTUVdimTEMP3 + maxOrbTEMP3
               IF(MEMUSAGE3 .LT. MAXMEM)THEN
                  maxTUVdim = maxTUVdimTEMP3 
                  maxPrim = maxPrimTEMP3 
                  maxOrb = maxOrbTEMP3 
                  MaxPassesForType(iPassType) = trialnpass
               ELSE
                  EXIT
               ENDIF
            ELSE
               EXIT
            ENDIF
         else
            EXIT
         endif
      enddo
   enddo
else
   do iPassType=1,nPassTypes
      MaxPassesForType(iPassType) = 1
   enddo
endif   

end subroutine determineMaxPassesForType

!> \brief dealloc the integral item
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param PQ contains integrand info to be deallocd
!> \param Integral contains arrays to store intermidiates and final integrals, to be deallocd
SUBROUTINE deallocIntegrals(PQ,Integral)
  use memory_handling
IMPLICIT NONE
TYPE(Integrand)    :: PQ
TYPE(Integralitem) :: Integral
integer (kind=long) :: nsize

nsize = mem_realsize*size(PQ%distance)  + mem_realsize*size(PQ%squaredDistance) + &
        & mem_realsize*size(PQ%exponents) + mem_realsize*size(PQ%reducedExponents)+ &
        & mem_realsize*size(PQ%integralPrefactor) !- mem_intsize*size(PQ%iprimP) - &
call mem_deallocated_mem_integrand(nsize)

nsize = mem_realsize*size(Integral%IN)  + mem_realsize*size(Integral%OUT) + &
        & mem_realsize*size(Integral%RTUV) + mem_realsize*size(Integral%tuvQ) + &
        & mem_realsize*size(Integral%integralsABCD)
call mem_deallocated_mem_integralitem(nsize)

call mem_dealloc(PQ%distance)
call mem_dealloc(PQ%squaredDistance)
call mem_dealloc(PQ%exponents)
call mem_dealloc(PQ%reducedExponents)
call mem_dealloc(PQ%integralPrefactor)
!call mem_dealloc(PQ%iprimP)
!call mem_dealloc(PQ%iprimQ)
call mem_dealloc(Integral%IN)
call mem_dealloc(Integral%OUT)
call mem_dealloc(Integral%RTUV)
call mem_dealloc(Integral%tuvQ)
call mem_dealloc(Integral%integralsABCD)

END SUBROUTINE deallocIntegrals

!> \brief wrapper routine which branch out and build the primitive hermite integral 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param integral contains arrays to store intermidiates and final integrals
!> \param input the integral input specifications, contains all info about what to do
!> \param PQ contains integrand info like reduced exponents, integral prefactor, etc
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE Build_HermiteTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT)
implicit none     
TYPE(Integrand)         :: PQ
TYPE(Integralitem)      :: integral
TYPE(Integralinput)     :: input
INTEGER                 :: LUPRI,IPRINT

IF(PQ%Operator.NE.CarmomOperator.AND.PQ%Operator.NE.MulmomOperator)THEN
   CALL GET_WTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT,PQ%Operator)
ELSE
   SELECT CASE(PQ%Operator)
   CASE (CarmomOperator)
      !Cartesian multipole order. Center at input%origo
      !Default is 0,0,0 center of coordinat system.
      CALL GETMOMTUV(INTEGRAL,INPUT,PQ,PQ%nPrimitives,LUPRI,IPRINT,INPUT%CMORDER)
   CASE (MulmomOperator)
      !multipolemoments. Same as Carmom but center of the multipole 
      !moment integral is the center of the contracted ODbatch.
      CALL GETMOMTUV(INTEGRAL,INPUT,PQ,PQ%nPrimitives,LUPRI,IPRINT,INPUT%MMORDER)
   CASE DEFAULT
      WRITE(LUPRI,'(1X,A,A16)') 'Programming error! Not a case in Build_HermiteTUV:',PQ%Operator
      CALL LSQUIT('Programming error! Not a case in Build_HermiteTUV',LUPRI)
   END SELECT
ENDIF
END SUBROUTINE Build_HermiteTUV

!> \brief builds the primitive hermite integral  
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param integral contains arrays to store intermidiates and final integrals
!> \param input the integral input specifications, contains all info about what to do
!> \param PQ contains integrand info like reduced exponents, integral prefactor, etc
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
!> \param OPERATORLABEL the label describing the operator 
SUBROUTINE GET_WTUV(INTEGRAL,INPUT,PQ,LUPRI,IPRINT,OPERATORLABEL)
! The final S(T,U,V) integrals are arranged as follows:
! S(000)
! S(100) S(010) S(001)
! S(200) S(110) S(101) S(020) S(011) S(002)
! S(300) S(210) S(201) S(120) S(111) S(102) S(030) S(021) S(012) S(003)
implicit none
integer                 :: Operatorlabel
TYPE(Integrand)         :: PQ
TYPE(Integralitem)      :: integral
TYPE(Integralinput)     :: input
real(realk),pointer     :: SJ0002(:,:)
real(realk),pointer     :: temp(:)
real(realk),pointer     :: wtuv(:)
INTEGER                 :: LUPRI,SUM,J,K,T,U,V,TUV,IOFF
INTEGER                 :: nPrim,IPRINT,ntuv,L,I,offsetTUV,offset,offset2
INTEGER                 :: zeroX,zeroY,zeroZ,Jmax,Jstart,nPrimP,nPassQ,iPass
real(realk)             :: X0,Y0,Z0
Integer                 :: der,ntuvFull,fullOP,ioffP,addDer
Logical                 :: kinetic,facone,efP,efQ

kinetic = INPUT%operator .EQ. KineticOperator
!efP     = PQ%P%p%type_elField
!efQ     = PQ%Q%p%type_elField
!
NPrim=PQ%nPrimitives
JMAX=PQ%endAngmom
Jstart=PQ%startAngmom

addDer = 0
IF (kinetic) THEN
   addDer = 2
!ELSE IF (efP) THEN
!   addDer = 1
!ELSE IF (efQ) THEN
!   addDer = 1
ENDIF

JMAX=PQ%endAngmom + addDer
Jstart=PQ%startAngmom + addDer

#ifdef VAR_LSDEBUGINT
IF(nPrim*(JMAX+1).GT.allocIntmaxTUVdim)THEN
   call lsquit('GET_WTUV alloc error1',lupri)
ENDIF
#endif

SELECT CASE(OPERATORLABEL)
CASE (CoulombOperator)
   CALL buildRJ000(Integral%IN,PQ,nPrim,JMAX,LUPRI,IPRINT,integral,INPUT%HIGH_RJ000_ACCURACY)
CASE (OverlapOperator)
   CALL buildSJ000(Integral%IN,PQ,nPrim,JMAX,LUPRI,IPRINT)
CASE (KineticOperator)
   CALL buildSJ000(Integral%IN,PQ,nPrim,JMAX,LUPRI,IPRINT)
CASE (CAMOperator)
   call mem_alloc(SJ0002,JMAX,nPrim,.TRUE.,.FALSE.)
   CALL buildRJ000(Integral%IN,PQ,nPrim,JMAX,LUPRI,IPRINT,integral,INPUT%HIGH_RJ000_ACCURACY)
   CALL buildErfRJ000(SJ0002,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,INPUT%ATTomega,INPUT%HIGH_RJ000_ACCURACY)
   ! alpha /r12 
   CALL DSCAL(NPrim*(JMAX+1),INPUT%ATTalpha,Integral%IN,1)
   ! (alpha + beta erf(\omega r12) )/r12 
   CALL DAXPY(NPrim*(JMAX+1),INPUT%ATTbeta,SJ0002,1,Integral%IN,1)
   call mem_dealloc(SJ0002)
CASE (ErfcOperator)
   ! ( 1 - erf(\omega r12) )/r12
   call mem_alloc(SJ0002,JMAX,nPrim,.TRUE.,.FALSE.)
   CALL buildRJ000(Integral%IN,PQ,nPrim,JMAX,LUPRI,IPRINT,integral,INPUT%HIGH_RJ000_ACCURACY)
   CALL buildErfRJ000(SJ0002,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,INPUT%ATTomega,INPUT%HIGH_RJ000_ACCURACY)
   CALL DAXPY(NPrim*(JMAX+1),-1.0E0_realk,SJ0002,1,Integral%IN,1)
   call mem_dealloc(SJ0002)
CASE (ErfOperator)
   CALL buildErfRJ000(Integral%IN,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,INPUT%ATTomega,INPUT%HIGH_RJ000_ACCURACY)
CASE (NucpotOperator)
   CALL buildNuclearRJ000(Integral%IN,PQ,nPrim,JMAX,LUPRI,IPRINT,integral,INPUT%HIGH_RJ000_ACCURACY)
CASE (GGemOperator)
   CALL buildGJ000(Integral%IN,nPrim,JMAX,PQ%integralPrefactor,PQ%reducedExponents,PQ%squaredDistance,&
     &             INPUT%GGem%coeff,INPUT%GGem%exponent,input%GGem%n)
CASE (GGemQuaOperator)
   CALL buildGJ000(Integral%IN,nPrim,JMAX,PQ%integralPrefactor,PQ%reducedExponents,PQ%squaredDistance,&
     &             INPUT%GGem%coeff,INPUT%GGem%exponent,input%GGem%n)
CASE (GGemCouOperator)
   CALL buildGJ000Coulomb(Integral%IN,nPrim,JMAX,PQ%integralPrefactor,PQ%reducedExponents,PQ%squaredDistance,&
     &                    INPUT%GGem%coeff,INPUT%GGem%exponent,input%GGem%n,integral,INPUT%HIGH_RJ000_ACCURACY,&
     &                    LUPRI,IPRINT)
CASE (GGemGrdOperator)
   CALL buildGJ000Grad(Integral%IN,nPrim,JMAX,PQ%integralPrefactor,PQ%reducedExponents,PQ%squaredDistance,&
     &             INPUT%GGem%coeff,INPUT%GGem%exponent,INPUT%GGem%expProd,input%GGem%n)
CASE DEFAULT
  print*,'Programming error! Not a case in GET_WTUV:',OPERATORLABEL
  WRITE(LUPRI,'(1X,A,A16)') 'Programming error! Not a case in GET_WTUV:',OPERATORLABEL
  CALL LSQUIT('Programming error! Not a case in GET_WTUV',LUPRI)
END SELECT

IF (IPRINT .GE. 10) THEN
   CALL LSHEADER(LUPRI,'Output from W000')
   DO I=1,NPrim
      DO J=0,Jmax
         WRITE(LUPRI,'(2X,A6,I4,A1,I4,A2,ES16.8)')'SJ000(',J,',',I,')=',Integral%IN(1+J+(I-1)*(Jmax+1))
      ENDDO
   ENDDO
END IF

ntuv=(JMAX+1)*(JMAX+2)*(JMAX+3)/6-Jstart*(Jstart+1)*(Jstart+2)/6
INTEGRAL%nTUV=ntuv
Integral%nEFG=1 !no cartesian multipole moments

#ifdef VAR_LSDEBUGINT
IF(nPrim*nTUV.GT.allocIntmaxTUVdim)THEN
   call lsquit('GET_WTUV alloc error2',lupri)
ENDIF
#endif

IF (kinetic) THEN
  wtuv => Integral%OUT
ELSE
  wtuv => Integral%Rtuv
ENDIF

CALL WTUVrecurrence(wtuv,Integral%IN,Integral%TUV,PQ%distance,&
     &              jStart,jMax,nPrim,ntuv,lupri,iprint)


IF (kinetic) THEN
  der = input%GEODERIVORDER
  facone = (der.EQ.1).AND..NOT.integral%dohodi
  ntuvFull = ntuv
  fullOP = Jstart*(Jstart+1)*(Jstart+2)/6
  JMAX=PQ%endAngmom
  Jstart=PQ%startAngmom
  ntuv=(JMAX+1)*(JMAX+2)*(JMAX+3)/6-Jstart*(Jstart+1)*(Jstart+2)/6
  ioffP = Jstart*(Jstart+1)*(Jstart+2)/6
  CALL WTUV_extract_kinetic(Integral%Rtuv,wtuv,Integral%TUV%TUVindex,&
     &                      Jstart,jmax,nPrim,ioffP,fullOP,ntuvfull,ntuv,1,facone,LUPRI,IPRINT)
  INTEGRAL%nTUV=ntuv
ENDIF

IF (IPRINT .GE. 10) THEN
   CALL LSHEADER(LUPRI,'Output from WTUVrecurrence')
   WRITE (LUPRI,'(2X,A,I10)') 'JSTART', JSTART
   WRITE (LUPRI,'(2X,A,I10)') 'JMAX  ', JMAX
   WRITE (LUPRI,'(2X,A,I10)') 'NPrim  ', NPrim
   WRITE (LUPRI,'(2X,A,I10)') 'NTUV  ', nTUV
   IF (IPRINT .GE. 20) THEN
      CALL LSHEADER(LUPRI,'Hermite integrals S(t,u,v)')
      DO J = Jstart, JMAX
         DO T = J,0,-1
            DO U = J-T,0,-1
               V=J-T-U
               TUV=Integral%TUV%tuvIndex(T,U,V)-Jstart*(Jstart+1)*(Jstart+2)/6
               WRITE (LUPRI,'(2X,A2,I1,A1,I1,A1,I1,A1,2X,5ES16.8/,(12X,5ES16.8))')&
                    & 'W(',T,',',U,',',V,')', &
                    &(INTEGRAL%Rtuv(I + (TUV-1)*nPrim),I=1,NPrim)
               WRITE (LUPRI,*) ' '
            ENDDO
         ENDDO
      ENDDO
   END IF
END IF

Integral%nPrim=nPrim

IF((Input%addtointegral.AND.PQ%Q%p%orbital1%TYPE_Nucleus).AND.&
     &(OPERATORLABEL.EQ.NucpotOperator.AND.PQ%Q%p%nPasses.GT.1))THEN
   IF(INPUT%geoderOrderP.EQ.0)THEN
      ! We sum all the Passes in Overlap PassQ (all same charge)
      ! it is more efficient to do the sum here 
      call mem_alloc(temp,nPrim*Integral%nTUV)
      DO TUV = 1,Integral%nTUV 
         offsetTUV=(TUV-1)*nPrim
         DO I=1,NPrim
            temp(I+offsetTUV)= INTEGRAL%Rtuv(I+offsetTUV)
         ENDDO
      ENDDO
      nPrimP = PQ%P%p%nPrimitives
      nPassQ = PQ%Q%p%nPasses
      !temp(nPassQ,nPrimP,nTUV) => Rtuv(nPrimP,nTUV) 
      DO TUV = 1,Integral%nTUV 
         offset=(TUV-1)*nPrimP
         DO I=1,NPrimP
            iPass=1         
            offset2=(I-1)*NPassQ+(TUV-1)*nPrim
            INTEGRAL%Rtuv(I+offset)= temp(iPass + offset2)
            DO iPass=2,nPassQ
               INTEGRAL%Rtuv(I+offset)= INTEGRAL%Rtuv(I+offset)+temp(iPass + offset2) 
            ENDDO
         ENDDO
      ENDDO
      call mem_dealloc(temp)

      PQ%Q%p%nPasses = 1
      Integral%nPrim = nPrimP
      nPrim = nPrimP
   ENDIF
ENDIF

END SUBROUTINE GET_WTUV

!> \brief build the Classical contribution to the gradient wrt nuclear position through FMM
!> \author Andreas Krapp
!> \date 2010
!> \param INTOUT the integral output specifications, determines how the output should be given
!> \param natoms number of atoms
!> \param input the integral input specifications, contains all info about what to do
SUBROUTINE GradClassical(INTOUT,natoms,Input)
implicit none
TYPE(INTEGRALOUTPUT) :: INTOUT
INTEGER              :: NATOMS
TYPE(INTEGRALINPUT)  :: INPUT
Real(realk),pointer :: gradient(:)
!
integer :: i,j
call mem_alloc(gradient,3*natoms)
gradient = 0.0E0_realk
call GradClassical1(gradient,natoms,Input)
! add to the total gradient
do i = 1, nAtoms
   do j = 1, 3
      INTOUT%resultTensor%LSAO(i)%elms(j) =   &
           &INTOUT%resultTensor%LSAO(i)%elms(j) + gradient((i-1)*3+j)
   end do
end do
call mem_dealloc(gradient)
END SUBROUTINE GradClassical

SUBROUTINE GradClassicalGrad(gradient,natoms,Input)
implicit none
INTEGER              :: NATOMS
TYPE(INTEGRALINPUT)  :: INPUT
Real(realk) :: gradient(3*natoms)
Real(realk),pointer :: gradient2(:)
integer :: i,j
call mem_alloc(gradient2,3*natoms)
gradient2 = 0.0E0_realk
call GradClassical1(gradient2,natoms,Input)
! add to the total gradient
!do i = 1, nAtoms
!   gradient(i) = gradient(i) + gradient2(1+(i-1)*3)
!   gradient(i+nAtoms) = gradient(i+nAtoms) + gradient2(2+(i-1)*3)
!   gradient(i+2*nAtoms) = gradient(i+2*nAtoms) + gradient2(3+(i-1)*3)
!enddo
do i = 1, nAtoms*3
   gradient(i) = gradient(i) + gradient2(i)
enddo
call mem_dealloc(gradient2)
END SUBROUTINE GradClassicalGrad

SUBROUTINE GradClassical1(gradient,natoms,Input)
implicit none
TYPE(INTEGRALINPUT)  :: INPUT
integer,intent(in) :: NATOMS
Real(realk) :: gradient(3*natoms)

Integer,parameter   :: LWRK = 1
Real(realk)         :: WRK(LWRK)
!
INTEGER     :: NLEVEL,LMAX,TLMAX,LEXTRA,ALGORITHM,T_CONTRACTOR,I,J
REAL(REALK) :: GRAIN,RPQMIN,SCREEN,FAC
LOGICAL     :: USEUMAT,BRFREE,GRSRCH,INCNN
LOGICAL     :: DYNLMAX,ALLSQR,TRSRCH,NOBOXP,CONTRACTED
LOGICAL     :: full

! Set default values for the mm driver
CALL mm_init_defs(NLEVEL,GRAIN,TLMAX,ALGORITHM,                            &
     &            USEUMAT,BRFREE,RPQMIN,DYNLMAX,LEXTRA,                    &
     &            ALLSQR,INCNN,TRSRCH,NOBOXP,GRSRCH,T_CONTRACTOR,SCREEN)

!Simen find appropriate variable to determine if the basis set is contracted/uncontracted
contracted = .TRUE.

! Change default settings based on LSint-settings and input
CALL mm_init_scheme(NLEVEL,GRAIN,Input%MM_LMAX,Input%MM_TLMAX,ALGORITHM,  &
     &              USEUMAT,BRFREE,RPQMIN,DYNLMAX,LEXTRA,                  &
     &              ALLSQR,INCNN,TRSRCH,NOBOXP,GRSRCH,T_CONTRACTOR,        &
     &              Input%MM_SCREENTHR,contracted)

full = .TRUE.
IF (Input%MM_NOONE) full = .FALSE.

IF (full) THEN
   ! Both two-electron repulsion and nuclear-electron attraction contributions to the classical Coulomb matrix
   write(*,*) 'Attention: the two electron gradient will include '
   write(*,*) '           the classical part of the nuclear repulsion'
   write(*,*) '           the classical part of the nuclear attraction and'
   write(*,*) '           the total electron repulsion contributions!'
   CALL mm_get_grad('NUC_EL',gradient,natoms,WRK,LWRK)
ELSE
   ! Only two-electron contribution
   CALL mm_get_grad('TWO_EL',gradient,natoms,WRK,LWRK)
ENDIF
END SUBROUTINE GradClassical1

!> \brief build the Classical Coulomb-matrix contribution through FMM
!> \author S. Reine
!> \date 2010
!> \param INTOUT the integral output specifications, determines how the output should be given
!> \param n1 the size of the first dimension
!> \param n2 the size of the second dimension
!> \param input the integral input specifications, contains all info about what to do
SUBROUTINE JmatClassical(INTOUT,n1,n2,Input)
implicit none
Integer     :: n1,n2
TYPE(INTEGRALOUTPUT) :: INTOUT
TYPE(INTEGRALINPUT)  :: INPUT
!
Real(realk),pointer :: Jmat(:,:)
call mem_alloc(Jmat,n1,n2)
Jmat = 0E0_realk
call JmatClassical1(Jmat,n1,n2,Input)
call add_full_2dim_to_lstensor(intout%resultTensor,Jmat,n1,n2,1)
call mem_dealloc(Jmat)
END SUBROUTINE JmatClassical

SUBROUTINE JmatClassicalMAT(MAT,n1,n2,Input)
implicit none
Integer     :: n1,n2
type(matrix) :: MAT
type(matrix) :: TMP
TYPE(INTEGRALINPUT)  :: INPUT
!
Real(realk),pointer :: Jmat(:,:)
call mem_alloc(Jmat,n1,n2)
Jmat = 0E0_realk
call JmatClassical1(Jmat,n1,n2,Input)
call mat_init(TMP,MAT%nrow,MAT%ncol)
call mat_set_from_full(Jmat,1.0E0_realk,TMP)
call mem_dealloc(Jmat)
call mat_daxpy(1.0E0_realk,TMP,MAT)
call mat_free(TMP)
END SUBROUTINE JmatClassicalMAT

subroutine JmatClassical1(Jmat,n1,n2,Input)
  implicit none
  Integer     :: n1,n2
  Real(realk) :: Jmat(n1,n2)
  TYPE(INTEGRALINPUT)  :: INPUT    
  ! Dummy argument for mm_get_J_matrix - should be removed
  Integer,parameter :: LWRK = 1
  Real(realk)       :: WRK(LWRK)
  !
  INTEGER     :: NLEVEL,LMAX,TLMAX,LEXTRA,ALGORITHM,T_CONTRACTOR
  REAL(REALK) :: GRAIN,RPQMIN,SCREEN
  LOGICAL     :: USEUMAT,BRFREE,GRSRCH,INCNN
  LOGICAL     :: DYNLMAX,ALLSQR,TRSRCH,NOBOXP,CONTRACTED
  LOGICAL :: full
  !   Set default values for the mm driver
  CALL mm_init_defs(NLEVEL,GRAIN,TLMAX,ALGORITHM,USEUMAT,BRFREE,RPQMIN,&
       & DYNLMAX,LEXTRA,ALLSQR,INCNN,TRSRCH,NOBOXP,GRSRCH,T_CONTRACTOR,SCREEN)
  
  !Simen find appropriate variable to determine if the basis set is contracted/uncontracted
  contracted = .TRUE.
  
  !Change default settings based on LSint-settings and input
  CALL mm_init_scheme(NLEVEL,GRAIN,Input%MM_LMAX,Input%MM_TLMAX,ALGORITHM,&
       & USEUMAT,BRFREE,RPQMIN,DYNLMAX,LEXTRA,ALLSQR,INCNN,TRSRCH,NOBOXP,&
       & GRSRCH,T_CONTRACTOR,Input%MM_SCREENTHR,contracted)
  ! Calculate moments here?
  
  full = .TRUE.
  !do_increment does not exist, find the appropriate variable
  !   IF (do_increment) full = .FALSE.
  IF (Input%MM_NOONE) full = .FALSE.
  full = .FALSE. !Problem for the energy for closed shell systems
  
  IF (full) THEN
     !Both two-electron repulsion and nuclear-electron attraction contributions to the classical Coulomb matrix
     !Simen Hack to cope with the (strange) order for density-fitting FMM
     CALL mm_get_J_matrix('FULL_J',Jmat,n2,n1,WRK,LWRK)
     !CALL mm_get_J_matrix('FULL_J',Jmat,n1,n2,WRK,LWRK)
  ELSE
     !Only two-electron contribution
     !Simen Hack to cope with the (strange) order for density-fitting FMM
     CALL mm_get_J_matrix('TWO_EL',Jmat,n2,n1,WRK,LWRK)
     !CALL mm_get_J_matrix('TWO_EL',Jmat,n1,n2,WRK,LWRK)
  ENDIF
END SUBROUTINE JmatClassical1

!> \brief build the Classical electron-nuclear attraction contribution through FMM
!> \author S. Reine
!> \date 2010
!> \param INTOUT the integral output specifications, determines how the output should be given
!> \param n1 the size of the first dimension
!> \param n2 the size of the second dimension
!> \param input the integral input specifications, contains all info about what to do
SUBROUTINE electronNuclearClassic(INTOUT,n1,n2,Input)
implicit none
Integer     :: n1,n2
TYPE(INTEGRALOUTPUT) :: INTOUT
TYPE(INTEGRALINPUT)  :: INPUT

! Dummy argument for mm_get_J_matrix - should be removed
Integer,parameter :: LWRK = 1
Real(realk)       :: WRK(LWRK)
Real(realk),pointer :: H1(:,:)
!
INTEGER     :: NLEVEL,LMAX,TLMAX,LEXTRA,ALGORITHM,T_CONTRACTOR
REAL(REALK) :: GRAIN,RPQMIN,SCREEN
LOGICAL     :: USEUMAT,BRFREE,GRSRCH,INCNN
LOGICAL     :: DYNLMAX,ALLSQR,TRSRCH,NOBOXP,CONTRACTED

!
call mem_alloc(H1,n1,n2)
H1 = 0E0_realk
!   Set default values for the mm driver
CALL mm_init_defs(NLEVEL,GRAIN,TLMAX,ALGORITHM,                            &
     &            USEUMAT,BRFREE,RPQMIN,DYNLMAX,LEXTRA,                    &
     &            ALLSQR,INCNN,TRSRCH,NOBOXP,GRSRCH,T_CONTRACTOR,SCREEN)

!Simen find appropriate variable to determine if the basis set is contracted/uncontracted
contracted = .TRUE.

!   Change default settings based on LSint-settings and input
CALL mm_init_scheme(NLEVEL,GRAIN,Input%MM_LMAX,Input%MM_TLMAX,ALGORITHM,&
     USEUMAT,BRFREE,RPQMIN,DYNLMAX,LEXTRA,                  &
     ALLSQR,INCNN,TRSRCH,NOBOXP,GRSRCH,T_CONTRACTOR,        &
     Input%MM_SCREENTHR,contracted)
! Calculate moments here?

!   Classical nuclear-electron attraction contribution
CALL mm_get_J_matrix('ONE_EL',H1,n1,n2,WRK,LWRK)
call add_full_2dim_to_lstensor(intout%resultTensor,H1,n1,n2,1)
call mem_dealloc(H1)

END SUBROUTINE electronNuclearClassic

!> \brief determine if the integral should be calculated or screened away
!> \author S. Reine
!> \date 2010
!> \param P contains overlap distribution for the left hand side, electron 1
!> \param Q contains overlap distribution for the right hand side, electron 2
!> \param input the integral input specifications, contains all info about what to do
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
LOGICAL FUNCTION getScreening(P,Q,INPUT,LUPRI,IPRINT)
implicit none
TYPE(Overlap)        :: P
TYPE(Overlap)        :: Q
TYPE(INTEGRALINPUT)  :: INPUT
integer              :: LUPRI,IPRINT
!
logical     :: screen
Real(realk) :: distance,distance2,extentSum,nomX,denomX,ERFC_Rpq,ERF_Rpq,CAMcorr
Integer     :: dir
!MBIE
real(realk) :: RAB,RCD,X,Y,Z,RPQ,RPRIME,MBIE0,MBIE1,MBIEEST
real(realk),parameter :: D0=0E0_realk,D2=2E0_realk,D3=3E0_realk,D10=10E0_realk
real(realk),parameter :: fourPim = 4E0_realk/3.14159265358979323846E0_realk !4/pi
real(realk),parameter :: twosqPim = 1.1283791670955126E0_realk !2/sqrt(pi)

screen = .FALSE.
getScreening = .FALSE.
!Cauchy-Schwarz screening
IF (INPUT%CS_SCREEN) THEN
  screen = (P%maxGab+Q%maxGab.LE.INPUT%CS_THRLOG)
  IF (screen) THEN
    getScreening = .TRUE.
    RETURN
  ENDIF
ENDIF
!Classical or extent based screening
IF (INPUT%NonClassical_SCREEN.OR.INPUT%OE_SCREEN) THEN
  distance2 = 0.0E0_realk
  DO dir=1,3
    distance = P%ODcenter(dir) - Q%ODcenter(dir)
    distance2 = distance2 + distance*distance
  ENDDO
  extentSum = P%ODextent + Q%ODextent
  screen = distance2.GT.(extentSum*extentSum)
  IF (screen) THEN
    getScreening = .TRUE.
    RETURN
  ENDIF
ENDIF
IF(INPUT%MBIE_SCREEN)THEN
   RAB = P%ODextent
   RCD = Q%ODextent
   X=P%ODcenter(1) - Q%ODCENTER(1) 
   Y=P%ODcenter(2) - Q%ODCENTER(2) 
   Z=P%ODcenter(3) - Q%ODCENTER(3) 
   RPQ = SQRT(X*X+Y*Y+Z*Z)
   IF(RPQ .GT. RAB+RCD+2)THEN 
      ! note that you will get horrible estimates if RPQ is only slightly larger than RAB+RCD+1 we therefore 
      ! use RPQ.GT.RAB+RCD+2 before we use the estimates of MBIE, compared to RPQ.GT.RAB+RCD+1 used in (JCP 123,184101) 
      ! In order to obtain good estimates the Rab and Rcd should be as small as possible while still large enough 
      ! to ensure that the multipole series converges. This is a bit tricky for small distances, but the value of eta
      ! in ODbatches.f90 seem to allow this.   
      Rprime = RPQ-RAB-RCD
      !Note that MBIE-1 can be larger than MBIE-0 if RPQ is not much larger than RAB+RCD
      MBIE0=ABS(P%MBIE(1)*Q%MBIE(1)/(Rprime-1))
      MBIE1=ABS(P%MBIE(1)*Q%MBIE(1)/RPQ)+ABS((P%MBIE(2)*Q%MBIE(1)+P%MBIE(1)*Q%MBIE(2))/(Rprime*Rprime-Rprime))
      MBIEEST = MIN(MBIE0,MBIE1)
      IF(INPUT%operator.EQ.ErfcOperator)THEN
         !   use of short range exchange replacing the 
         !   1/r12 operator with (1-erf(\omega r12))/r12 = erfc(\omega r12))/r12
         !   For x > 0 erfc(x) is bounded by       
         !   2/(sqrt(pi)) \frac{e^{-x^{2}} }{ x + \sqrt{ x^{2} + 4pim }}
         X = Rpq*Input%AttOmega
         nomX = exp(-X*X)
         denomX = X + SQRT(X*X+fourPim)
         ERFC_Rpq   = twosqPim * nomX/denomX
!         print*,'MBIEEST ORIGINAL ',MBIEEST
         MBIEEST = ERFC_Rpq*MBIEEST
!         print*,'P%maxGab*Q%maxGab',10.0E0_realk**(P%maxGab+Q%maxGab)
!         print*,'ERFC_Rpq*MBIEEST ',ERFC_Rpq*MBIEEST,'          ',ERFC_Rpq,' RPQ',RPQ,'RPQ*omega',X
      ENDIF
      IF(INPUT%operator.EQ.CAMOperator)THEN
         call lsquit('THIS COMBI CAM and MBIE do not work, I do not know why. TK',-1)
         !   use the long range exchange replacing the 
         !   1/r12 operator with (alpha+beta*erf(\omega r12))/r12 
         !   For x > 0 erfc(x) is bounded FROM BELOW by       
         !   2/(sqrt(pi)) \frac{e^{-x^{2}} }{ x + \sqrt{ x^{2} + 2 }}
         !   which means that erf(=1-erfc) is bounded FROM ABOVE by
         !   1 - 2/(sqrt(pi)) \frac{e^{-x^{2}} }{ x + \sqrt{ x^{2} + 2 }}
         X = Rpq*Input%AttOmega
         nomX = exp(-X*X)
         denomX = X + SQRT(X*X+D2)
         ERF_Rpq   = 1 - twosqPim * nomX/denomX
         CAMcorr = input%ATTalpha + input%ATTbeta * ERF_Rpq !note that alpha is not  
         !put into the exchangefactor anymore -> TO CHECK for consistency
         MBIEEST = CAMcorr*MBIEEST
      ENDIF
      screen = (MBIEEST.LE.INPUT%CS_THRESHOLD)
      IF (screen) THEN
!         print*,'ORI MBIE EST',MIN(MBIE0,MBIE1), 'Oper   =',INPUT%operator(1:7)
!         IF(INPUT%operator(1:7).EQ.'CAM    ')print*,'MBIE EST',MBIEEST,'CAMcorr',CAMcorr
!         IF (INPUT%CS_SCREEN)print*,'P%maxGab*Q%maxGab',10.0E0_realk**(P%maxGab+Q%maxGab)
         getScreening = .TRUE.
         RETURN
      ENDIF
   ENDIF
ENDIF
END FUNCTION getScreening

!> \brief determine if the integral should be calculated or screened away by MBIE screening
!> \author T. Kjaergaard
!> \date 2010
!> \param screen logical which determine if integral should be screened away
!> \param Rab extent of the P overlap distribution
!> \param P contains overlap distribution for the left hand side, electron 1
!> \param Q contains overlap distribution for the right hand side, electron 2
!> \param mbieP the absolute spherical multipoles
!> \param CS_THRESHOLD the screening threshold
!> \param DRED the reduced density matrix element multiplied on the integral
!> \param LUPRI the logical unit number for the output file
SUBROUTINE MBIE_SCREENING(screen,RAB,P,mbieP,RCD,QODcenter,Qmbie1,Qmbie2,&
     & CS_THRESHOLD,DRED,nonSR_EXCHANGE,omega,lupri)
implicit none
logical,intent(inout) :: screen
real(realk),intent(in)  :: omega
integer,intent(in) :: lupri
type(overlap),intent(in) :: P
logical,intent(in) :: nonSR_EXCHANGE
real(realk),intent(in) :: RAB,mbieP(2),CS_THRESHOLD,QODcenter(3),RCD,Qmbie1,Qmbie2
integer(kind=short),intent(in) :: DRED
!
real(realk) :: RPQ,X,Y,Z,MBIE0,MBIE1,MBIEEST,Dredreal,ERFC_Rpq
real(realk) :: qx,qy,qz,px,py,pz,pqx,pqy,pqz,Rprime,nomX,denomX
integer :: iq,ip
real(realk),parameter :: D0=0E0_realk,D2=2E0_realk,D3=3E0_realk,D10=10E0_realk
real(realk),parameter :: fourPim = 4E0_realk/3.14159265358979323846E0_realk !4/pi
real(realk),parameter :: twosqPim = 1.1283791670955126E0_realk !2/sqrt(pi)

!RCD = Q%ODextent
X=P%ODcenter(1) - QODCENTER(1)
Y=P%ODcenter(2) - QODCENTER(2)
Z=P%ODcenter(3) - QODCENTER(3)
RPQ = SQRT(X*X+Y*Y+Z*Z)!braket distance
IF(RPQ .GT. RAB+RCD+2)THEN 
   DREDreal = D10**DRED
   ! note that you will get horrible estimates if RPQ is only slightly larger than RAB+RCD+1 we therefore 
   ! use RPQ.GT.RAB+RCD+2 before we use the estimates of MBIE, compared to RPQ.GT.RAB+RCD+1 used in (JCP 123,184101) 
   ! In order to obtain good estimates the Rab and Rcd should be as small as possible while still large enough 
   ! to ensure that the multipole series converges. This is a bit tricky for small distances, but the value of eta
   ! in ODbatches.f90 seem to allow this.   
   Rprime = RPQ-RAB-RCD
   !Note that MBIE-1 can be larger than MBIE-0 if RPQ is not much larger than RAB+RCD
   MBIE0=ABS(MBIEP(1)*QMBIE1/(Rprime-1))
   MBIE1=ABS(MBIEP(1)*QMBIE1/RPQ)+ABS((MBIEP(2)*QMBIE1+MBIEP(1)*QMBIE2)/(Rprime*Rprime-Rprime))
   IF(nonSR_EXCHANGE)THEN !default regular MBIE
      MBIEEST = MIN(MBIE0,MBIE1)
   ELSE  
      !   use of short range exchange replacing the 
      !   1/r12 operator with (1-erf(\omega r12))/r12
      !   For x > 0 erfc(x) is bounded by       
      !   2/(sqrt(pi)) \frac{e^{-x^{2}} }{ x + \sqrt{ x^{2} + 4pim }}
      X = Rpq*Omega
      nomX = exp(-X*X)
      denomX = X + SQRT(X*X+fourPim)
      ERFC_Rpq   = twosqPim * nomX/denomX
      MBIEEST = ERFC_Rpq*MIN(MBIE0,MBIE1)
   ENDIF
   screen = (MBIEEST*DREDreal.LE.CS_THRESHOLD)
ENDIF
END SUBROUTINE MBIE_SCREENING

!> \brief determine if the integral should be calculated or screened away
!> \author S. Reine
!> \date 2010
!> \param P contains overlap distribution for the left hand side, electron 1
!> \param Q contains overlap distribution for the right hand side, electron 2
!> \param input the integral input specifications, contains all info about what to do
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
LOGICAL FUNCTION getScreeningOD(P,Q,INPUT,LUPRI,IPRINT)
implicit none
TYPE(ODBATCH)        :: P
TYPE(ODBATCH)        :: Q
TYPE(INTEGRALINPUT)  :: INPUT
integer              :: LUPRI,IPRINT
!
logical     :: screen
Real(realk) :: distance,distance2,extentSum
Integer     :: dir

screen = .FALSE.
getScreeningOD = .FALSE.
!Cauchy-Schwarz screening
IF (INPUT%CS_SCREEN) THEN
  screen = (P%maxGab*Q%maxGab.LE.INPUT%CS_THRESHOLD)
  IF (screen) THEN
    getScreeningOD = .TRUE.
    RETURN
  ENDIF
ENDIF
!Classical screening
IF (INPUT%NonClassical_SCREEN.OR.INPUT%OE_SCREEN) THEN
  distance2 = 0.0E0_realk
  DO dir=1,3
    distance = P%ODcenter(dir) - Q%ODcenter(dir)
    distance2 = distance2 + distance*distance
  ENDDO
  extentSum = P%ODextent + Q%ODextent
  screen = distance2.GT.(extentSum*extentSum)
  IF (screen) THEN
    getScreeningOD = .TRUE.
    RETURN
  ENDIF
ENDIF

END FUNCTION getScreeningOD

!> \brief initialise the Linkshell item used in LinK
!> \author T. Kjaergaard
!> \date 2010
!> \param A the linkshell
!> \param I the size of the linkshell
subroutine LINKshell_init(A,I)
  use memory_handling
implicit none
type(LINKshell), intent(inout) :: A
integer, intent(in)            :: I
integer (kind=long)            :: nsize
   !NULLIFY(A%Belms)   done in mem_alloc
   !NULLIFY(A%IODelms)
   A%DIM = I
   if (I == 0) then
      call mem_alloc(A%Belms,1)
      call mem_alloc(A%IODelms,1)
   else
      call mem_alloc(A%Belms,I)
      call mem_alloc(A%IODelms,I)
   endif

   nsize = size(A%Belms)*mem_intsize + &
        & size(A%IODelms)*mem_intsize + mem_intsize
   call mem_allocated_mem_linkshell(nsize)
end subroutine LINKshell_init

!> \brief free the Linkshell item used in LinK
!> \author T. Kjaergaard
!> \date 2010
!> \param A the linkshell
subroutine LINKshell_free(A)
  use memory_handling
implicit none
type(LINKshell), intent(inout) :: A
integer (kind=long) :: nsize

nsize = size(A%Belms)*mem_intsize + &
        & size(A%IODelms)*mem_intsize + mem_intsize
call mem_deallocated_mem_linkshell(nsize)

call mem_dealloc(A%Belms)
call mem_dealloc(A%IODelms)
end subroutine LINKshell_free

!> \brief Perform contraction with contraction coefficients to obtain contracted integrals
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!>
!> \param Integral contains arrays to store intermidiates and final integrals
!> \param P the overlap distribution P or Q
!> \param iAngmom the index of angular momentum
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE ContractBasis(Integral,P,iAngmom,LUPRI,IPRINT)
implicit none
TYPE(Integralitem) :: Integral
TYPE(Overlap)      :: P
Integer            :: LUPRI,IPRINT,iAngmom
!
Integer               :: iA1,iA2,nCont,nC1,nC2,dim1,nP1,nP2,nPrim
real(realk),pointer   :: ptemp(:)
real(realk),pointer   :: TMPWORK(:)

IF(P%contractbasis)THEN
  nPrim = P%nPrimitives
  dim1 = Integral%nAng*Integral%nOrb
  iA1 = P%indexAng1(iangmom)
  iA2 = P%indexAng2(iangmom)
  nC1 = P%orbital1%nContracted(iA1)
  nC2 = P%orbital2%nContracted(iA2)
  nCont = nC1*nC2
  nP1 = P%orbital1%nPrimitives
  nP2 = P%orbital2%nPrimitives

#ifdef VAR_LSDEBUGINT
IF(nPrim*P%nPasses*Integral%nOrb*Integral%nAng.GT.allocIntmaxTUVdim)THEN
   call lsquit('ContractBasisGen alloc error1',lupri)
ENDIF
IF(nCont*P%nPasses*Integral%nOrb*Integral%nAng.GT.allocIntmaxTUVdim)THEN
   call lsquit('ContractBasisGen alloc error2',lupri)
ENDIF
#endif
!print*,'debug TMPWORK564545465',nCont*nP1*nP2
call mem_workpointer_alloc(TMPWORK,nCont*nP1*nP2)

IF(nCont.LE. 5)THEN !the general case
   CALL ContractBasisGen(Integral%IN,Integral%OUT,P,nPrim,& 
        &    nCont,P%nPasses,Integral%nAng,Integral%nOrb,&
        &    TMPWORK,nC1,nC2,iA1,iA2,LUPRI,IPRINT) 
ELSE
   !here we use dgemm - it is general - no assumptions
   CALL ContractBasisGenDgemm(Integral%IN,Integral%OUT,P,nPrim,& 
        &    nCont,P%nPasses,Integral%nAng,Integral%nOrb,&
        &    TMPWORK,nC1,nC2,iA1,iA2,LUPRI,IPRINT) 
ENDIF
call mem_workpointer_dealloc(TMPWORK)
  ptemp => Integral%IN 
  Integral%IN  => Integral%OUT 
  Integral%OUT => ptemp 
  Integral%nPrim = nCont*P%nPasses
  IF (IPRINT.GT. 50) THEN
     CALL PrintTensor(Integral%IN,'Contracted          ',Integral%nPrim,Integral%nOrb,Integral%nAng,Lupri,&
          & 'iPass ','iOrb  ','iAngP ',1)
  ENDIF
ENDIF

END SUBROUTINE ContractBasis

!> \brief Contract contraction coeffficients PA ordering \f$ OUT(nCont,nPasses,dim) = \sum_{nPrim} CC(nCont,nPrim)*IN(nPrim,nPasses,dim) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param PrimInt to be contracted with contraction coefficients (nPrim,nPasses,dim)
!> \param ContInt output after contraction (nCont,nPasses,dim)
!> \param CC empty input - just allocd memory - to be used for contactionmatrix  
!> \param P overlap which is contracted
!> \param nPrim number of primitive functions
!> \param nCont number of contracted functions
!> \param nPasses number of passes
!> \param nDim dimension of the unaffected dimension
!> \param nC1 number of contracted functions on center A (for P) or C (for Q)
!> \param nC2 number of contracted functions on center B (for P) or D (for Q)
!> \param iA1 angular moment or center A (for P) or C (for Q) (used for family basis set)
!> \param iA2 angular moment or center B (for P) or D (for Q) (used for family basis set)
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE ContractBasisGen(PrimInt,ContInt,P,nPrim,nCont,nPasses,nAng,nOrb,&
    &                       CC,nC1,nC2,iA1,iA2,LUPRI,IPRINT)
implicit none
Integer,intent(in) :: nPrim,nCont,nPasses,nAng,nOrb,nC1,nC2,iA1,iA2,LUPRI,IPRINT
Real(realk),intent(in),dimension(nPrim,nPasses*nOrb*nAng)    :: PrimInt!TEMP
Real(realk),intent(inout),dimension(nCont,nPasses*nOrb*nAng) :: ContInt!TEMP
TYPE(Overlap),intent(in) :: P
!
!Real(realk),intent(inout),dimension(nC1*nC2,nPrim)           :: CC
!Integer     :: iPass,iCont,iPrim
!Real(realk),parameter :: D1=1.0E0_realk,D0=0.0E0_realk
!Real(realk) :: TMP

Real(realk),intent(inout),dimension(nPrim,nC1*nC2)           :: CC
Integer     :: iPass,iCont,iPrim
Real(realk),parameter :: D1=1.0E0_realk,D0=0.0E0_realk
Real(realk) :: TMP

!CALL ConstructContraction_PA(CC,P,nPrim,nC1,nC2,iA1,iA2,LUPRI,IPRINT)
!DO iPass=1,nPasses*nOrb*nAng
!   DO iCont = 1,nCont
!      TMP = CC(iCont,1)*PrimInt(1,iPass)
!      DO iPrim = 2,nPrim
!         TMP = TMP+CC(iCont,iPrim)*PrimInt(iPrim,iPass)
!      ENDDO
!      ContInt(iCont,iPass) = TMP
!   ENDDO
!ENDDO
CALL ConstructContraction_PA2(CC,P,nPrim,nC1,nC2,iA1,iA2,LUPRI,IPRINT)
DO iPass=1,nPasses*nOrb*nAng
   DO iCont = 1,nCont
      TMP = D0
      DO iPrim = 1,nPrim
         TMP = TMP+CC(iPrim,iCont)*PrimInt(iPrim,iPass)
      ENDDO
      ContInt(iCont,iPass) = TMP
   ENDDO
ENDDO
IF (IPRINT.GT. 50) THEN
   CALL PrintTensor(ContInt,'Contracted          ',nCont*nPasses,nOrb,nAng,Lupri,&
        & 'iCont ','iOrb  ','iAngP ',1)
ENDIF
END SUBROUTINE ContractBasisGen

!> \brief Contract contraction coeffficients PA ordering \f$ OUT(nCont,nPasses,dim) = \sum_{nPrim} CC(nCont,nPrim)*IN(nPrim,nPasses,dim) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param PrimInt to be contracted with contraction coefficients (nPrim,nPasses,dim)
!> \param ContInt output after contraction (nCont,nPasses,dim)
!> \param CC empty input - just allocd memory - to be used for contactionmatrix  
!> \param P overlap which is contracted
!> \param nPrim number of primitive functions
!> \param nCont number of contracted functions
!> \param nPasses number of passes
!> \param nDim dimension of the unaffected dimension
!> \param nC1 number of contracted functions on center A (for P) or C (for Q)
!> \param nC2 number of contracted functions on center B (for P) or D (for Q)
!> \param iA1 angular moment or center A (for P) or C (for Q) (used for family basis set)
!> \param iA2 angular moment or center B (for P) or D (for Q) (used for family basis set)
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE ContractBasisGenDgemm(PrimInt,ContInt,P,nPrim,nCont,nPasses,nAng,nOrb,&
    &                       CC,nC1,nC2,iA1,iA2,LUPRI,IPRINT)
implicit none
Integer,intent(in) :: nPrim,nCont,nPasses,nAng,nOrb,nC1,nC2,iA1,iA2,LUPRI,IPRINT
Real(realk),intent(in),dimension(nPrim,nPasses*nOrb*nAng)    :: PrimInt!TEMP
Real(realk),intent(inout),dimension(nCont,nPasses*nOrb*nAng) :: ContInt!TEMP
TYPE(Overlap),intent(in) :: P
!
!Real(realk),intent(inout),dimension(nC1*nC2,nPrim)           :: CC
!Integer     :: iPass,iCont,iPrim
!Real(realk),parameter :: D1=1.0E0_realk,D0=0.0E0_realk
!Real(realk) :: TMP

Real(realk),intent(inout),dimension(nC1*nC2,nPrim)           :: CC
Integer     :: iPass,iCont,iPrim
Real(realk),parameter :: D1=1.0E0_realk,D0=0.0E0_realk
Real(realk) :: TMP

CALL ConstructContraction_PA(CC,P,nPrim,nC1,nC2,iA1,iA2,LUPRI,IPRINT)
CALL DGEMM('N','N',nCont,nOrb*nAng*nPasses,nPrim,D1,CC,nCont,&
        &     PrimInt,nPrim,D0,ContInt,nCont)
IF (IPRINT.GT. 50) THEN
   CALL PrintTensor(ContInt,'Contracted          ',nCont*nPasses,nOrb,nAng,Lupri,&
        & 'iCont ','iOrb  ','iAngP ',1)
ENDIF
END SUBROUTINE ContractBasisGenDgemm

!> \brief construct the contraction matrix order PA
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param CC the contraction matrix
!> \param P the overlap distribution P or Q
!> \param nP the number of primitives
!> \param nC1 the number of contracted functions on orbital 1
!> \param nC2 the number of contracted functions on orbital 2
!> \param iA1 the angular moment index on orbital 1
!> \param iA2 the angular moment index on orbital 2
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE ConstructContraction_PA2(CC,P,nP,nC1,nC2,iA1,iA2,LUPRI,IPRINT)
Implicit none
Type(Overlap) :: P
Integer       :: iA1,iA2,nP,nC1,nC2,LUPRI,IPRINT
Real(realk)   :: CC(nP,nC1,nC2)
!
Integer       :: iP,i1,i2,iC1,iC2,nP1,nP2

nP1 = P%orbital1%nPrimitives
nP2 = P%orbital2%nPrimitives
DO iP=1,nP
   i1 = P%iprim1(iP)
   i2 = P%iprim2(iP)
   DO iC2=1,nC2
      DO iC1=1,nC1
         CC(iP,iC1,iC2)=P%orbital1%CC(iA1)%p%elms(i1+(iC1-1)*nP1)*P%orbital2%CC(iA2)%p%elms(i2+(iC2-1)*nP2) 
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE ConstructContraction_PA2

INTEGER FUNCTION getIntegralDimsJengine(Alloc,maxPrimPass,ndmat,nderiv,magderiv)
implicit none
TYPE(Allocitem),intent(in) :: Alloc
integer,intent(in)         :: maxPrimPass(20),ndmat,nderiv,magderiv
!
Integer :: iLHS,iRHS,L,nTUV,nPrim,getMaxTUVdim

!Maximal number of Wtuv-components 
getMaxTUVdim = 0
DO iLHS=0,Alloc%maxAngmomLHS
  DO iRHS=0,Alloc%maxAngmomRHS
    L = iLHS + iRHS + nderiv + magderiv
    nTUV = (L+1)*(L+2)*(L+3)/6
    IF (iLHS+1.GT. 20) CALL LSQUIT('Error maxPrimAngmomLHS getIntegralDimsJengine',-1)
    IF (iRHS+1.GT. 20) CALL LSQUIT('Error maxPrimAngmomRHS getIntegralDimsJengine',-1)
    nPrim = Alloc%maxPrimAngmomLHS(iLHS+1)*max(maxPrimPass(iRHS+1),Alloc%maxPrimAngmomRHS(iRHS+1))
    getMaxTUVdim = max(getMaxTUVdim,nTUV*nPrim)
  ENDDO
ENDDO

!Maximal number of contracted components
getIntegralDimsJengine = max(getMaxTUVdim,Alloc%maxTotOrbLHS*ndmat,Alloc%maxijkLHS*ndmat)

END FUNCTION getIntegralDimsJengine

subroutine getIntegralDimsJengine2(Alloc,maxPrimPass,ndmat,nderiv,magderiv,maxTUVOUT,maxPrimPassOUT)
implicit none
TYPE(Allocitem),intent(in) :: Alloc
integer,intent(in)         :: maxPrimPass(20),ndmat,nderiv,magderiv
integer,intent(out)        :: maxTUVOUT,maxPrimPassOUT
!
Integer :: iLHS,iRHS,L,nTUV,nPrim

!Maximal number of Wtuv-components 
maxPrimPassOUT = 0
maxTUVOUT = 0 
DO iLHS=0,Alloc%maxAngmomLHS
  DO iRHS=0,Alloc%maxAngmomRHS
    L = iLHS + iRHS + nderiv + magderiv
    nTUV = (L+1)*(L+2)*(L+3)/6
    IF (iLHS+1.GT. 20) CALL LSQUIT('Error maxPrimAngmomLHS getIntegralDimsJengine',-1)
    IF (iRHS+1.GT. 20) CALL LSQUIT('Error maxPrimAngmomRHS getIntegralDimsJengine',-1)
    nPrim = Alloc%maxPrimAngmomLHS(iLHS+1)*max(maxPrimPass(iRHS+1),Alloc%maxPrimAngmomRHS(iRHS+1))
    maxPrimPassOUT = MAX(maxPrimPassOUT,nPrim)
    maxTUVOUT = MAX(maxTUVOUT,nTUV)
  ENDDO
ENDDO

END subroutine getIntegralDimsJengine2

SUBROUTINE setIntegralDerivativOrders(Integral,iDerLHS,derOrder,singleP,singleQ)
implicit none
integer,intent(IN)               :: iDerLHS,derOrder
logical,intent(IN)               :: singleP,singleQ
TYPE(integralitem),intent(INOUT) :: Integral
Integer :: iDerRHS
iDerRHS = derOrder - iDerLHS


Integral%lhsGeoOrder = iDerLHS
Integral%lhsGeoComp  = getODgeoComp(iDerLHS,singleP)
Integral%rhsGeoOrder = iDerRHS
Integral%rhsGeoComp  = getODgeoComp(iDerRHS,singleQ)

END SUBROUTINE setIntegralDerivativOrders

END MODULE integraldriver
