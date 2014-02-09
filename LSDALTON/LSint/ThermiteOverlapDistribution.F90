!> @file
!> Contains the orbital and overlap types, and subroutines to set these up
MODULE Thermite_OD
use typedef
use SphCart_matrices
!use ODbatches
use OD_type
use ls_util
use precision
use Integralparameters
use ThermiteMem_module
use OverlapType
TYPE Allocitem
Integer        :: maxPrimAngmomLHS(20)    !used for Rtuv-tensor
Integer        :: maxPrimAngmomRHS(20)    !used for Rtuv-tensor
!========================================= 
Integer,pointer        :: maxPrimLHSA(:)      !used for integrals+overlap 
Integer,pointer        :: maxPrimTUVLHSA(:)   !used for integrals+overlap 
Integer,pointer        :: maxnAngLHSA(:)      !used for overlap 
Integer,pointer        :: maxContLHSA(:)      !used for orbital 
Integer,pointer        :: maxAngmomOrbLHSA(:) !used for orbital 
Integer,pointer        :: maxTUVLHSA(:)       !used for F-type overlap 
Integer,pointer        :: maxTotOrbLHSA(:)    !used for temporary array (PQ) 
Integer,pointer        :: maxijkLHSA(:)       !used for temporary array (PQ) 
Integer,pointer        :: maxETUVlenLHSA(:)   !used for overlap (Ecoeff)
!========================================= 
Integer,pointer        :: maxPrimRHSA(:) 
Integer,pointer        :: maxPrimTUVRHSA(:) 
Integer,pointer        :: maxnAngRHSA(:) 
Integer,pointer        :: maxContRHSA(:) 
Integer,pointer        :: maxAngmomOrbRHSA(:) 
Integer,pointer        :: maxTUVRHSA(:) 
Integer,pointer        :: maxTotOrbRHSA(:) 
Integer,pointer        :: maxijkRHSA(:) 
Integer,pointer        :: maxETUVlenRHSA(:) 
!========================================== 
Integer                :: maxPrimLHS       !used for integrals+overlap
Integer                :: maxPrimTUVLHS    !used for integrals+overlap
Integer                :: maxnAngLHS       !used for overlap
Integer                :: maxContLHS       !used for orbital
Integer                :: maxAngmomOrbLHS  !used for orbital
Integer                :: maxAngmomLHS     !used for overlap
Integer                :: maxTUVLHS        !used for F-type overlap
Integer                :: maxTotOrbLHS     !used for temporary array (PQ)
Integer                :: maxijkLHS        !used for temporary array (PQ)
Integer                :: maxETUVlenLHS

Integer                :: maxPrimRHS
Integer                :: maxPrimTUVRHS
Integer                :: maxnAngRHS
Integer                :: maxContRHS
Integer                :: maxAngmomOrbRHS
Integer                :: maxAngmomRHS
Integer                :: maxTUVRHS
Integer                :: maxTotOrbRHS
Integer                :: maxijkRHS
Integer                :: maxETUVlenRHS
END TYPE Allocitem

!Simen Should remove all explicit dependence on IntegralItem from this file
!      and the move IntegralItem to ThermiteIntegrals.f90
!****** INTEGRAL
TYPE Integralitem
Real(realk),dimension(:),pointer    :: IN
Real(realk),dimension(:),pointer    :: OUT
Real(realk),dimension(:),pointer    :: RTUV
Real(realk),dimension(:),pointer    :: TUVQ
Real(realk),dimension(:),pointer    :: integralsABCD
Real(realk),dimension(:),pointer    :: ECONT
!
!Real(realk),dimension(:,:),pointer     :: Wtuv       !nPrimitives:ntuv
!Real(realk),dimension(:,:,:),pointer   :: tuvTUV     !ntuvQ(i):ntuvP:nPrimPQ
!Real(realk),dimension(:,:,:),pointer    :: TUVQ(:,:,:)     !ntuvP:nPrimP:nOrbQ
!Generic integral intermadiates used of contraction with Ecoefficients, 
!Spherical transformation, and contraction to contracted basis
!Real(realk),dimension(:,:,:),pointer    :: IntegralIN
!Real(realk),dimension(:,:,:),pointer    :: IntegralOUT
INTEGER                :: nAng
INTEGER                :: nPrim
INTEGER                :: nOrb
INTEGER                :: nTUV
INTEGER                :: nEFG !multipole orders
Logical                :: dohodi
INTEGER                :: nGeoDeriv
INTEGER                :: lhsGeoOrder
INTEGER                :: rhsGeoOrder
INTEGER                :: lhsGeoComp
INTEGER                :: rhsGeoComp
INTEGER                :: nMagDeriv
INTEGER                :: startDerivativeIndex
INTEGER                :: nDerivComp
INTEGER                :: nOperatorComp
LOGICAL                :: Jengine
TYPE(TUVitem),pointer  :: TUV
END TYPE Integralitem

CONTAINS
!> \brief One of two routines to set up the overlap (and the two corresponding orbitals) from the OD-batch
!> \author S. Reine and T. Kjaergaard
!> \date 2009-??-??
!> \param P The overlap
!> \param np The number of (significant) primitives
!> \param Input The integral input
!> \param sharedTUV The TUV indeces
!> \param Integral The integral specifications
!> \param Alloc Information about the maximal sizes for allocation purposes
!> \param ODB The OD-batch information
!> \param iElectron The electron (1 = 'LHS', 2 = 'RHS')
!> \param lupri The default print unit
!> \param iprint The print level (the higher the more information)
!> \param side The side ('LHS' or 'RHS')
!> Simen: either remove side or iElectron!!!
SUBROUTINE SET_OVERLAP(P,Input,sharedTUV,Integral,ODB,IELECTRON,LUPRI,IPRINT,DoAlloc)
use memory_handling
Implicit none
TYPE(IntegralInput),target,intent(in) :: Input
TYPE(IntegralItem),intent(inout)      :: Integral
TYPE(Overlap),intent(inout)           :: P
TYPE(ODBATCH),intent(in)              :: ODB
!Character*(*)                         :: side
Integer,intent(in)                    :: IELECTRON,LUPRI,IPRINT
TYPE(TUVitem),intent(in)              :: SharedTUV
LOGICAL,intent(in) :: DoAlloc
!
Integer             :: nGeoderivComp,idir,i1,i2,i12,start1,start2,end1,end2,orb1,orb2,start,offset
integer             :: l1,l2,ijk1,ijk2,ijk,ijkcart,maxijk,iangmom
Real(realk)         :: e1,e2,d2,signP,DMATmax
!Character(len=80)   :: t1,t2
Type(lstensor),pointer :: DMAT2
Type(lstensor),pointer :: GAB2
LOGICAL             :: LHS, DMATscreen,useFTUV,doscreen,screen,DMATpermute
Integer :: ndmat,idmat,indexETUV,l,nETUV,nTUV,maxangmom,minangmom,offset2
Integer :: iA12,iA1,iA2,startA2,maxnp,norder,nOperatorComp,maxBat,maxAng,offset1
Integer,pointer     :: lenEtuv(:)
Integer :: batchA,batchB,atomA,atomB,sAA,sA,sB,iSA,iSB,Dindex,elms,s1,s2,n1
Integer :: dimA,dimB,Gindex,nmom,nprim1,nprim2,np,ndim5,CMorder,maxgabelm,maxGab
integer :: iDer
integer(kind=short) :: maxPrimGab,maxPrimGabElm,Dmatmax2
integer(kind=short),pointer :: primgabmat(:,:)
!
maxPrimGabElm=shortzero
DMATscreen = .FALSE.
NULLIFY(DMAT2)
NULLIFY(GAB2)
!
NMOM = 1
DMATpermute = .FALSE.
IF(IELECTRON.EQ. 1)THEN
   ndmat = INPUT%NDMAT_LHS
   IF(INPUT%PS_SCREEN.OR.INPUT%CS_SCREEN.OR.Input%MBIE_SCREEN)THEN
      GAB2 => INPUT%LST_GAB_LHS   
   ENDIF
   IF(INPUT%PS_SCREEN)THEN
      maxPrimGabElm=INPUT%PS_MAXELM_RHS
   ENDIF
   IF(INPUT%CS_SCREEN)THEN
!ToDo Make as default for gradient calculations
      IF (Input%DO_JENGINE.AND.Input%DO_DAJENGINE) THEN
         DMATscreen = .TRUE.
         DMAT2 => INPUT%LST_DLHS
      ENDIF
      IF (Input%DO_DACOULOMB) THEN
         DMATscreen = .TRUE.
         DMAT2 => INPUT%LST_DLHS
      ENDIF
   ENDIF
   LHS=.TRUE.
   ndim5=INPUT%nMagDerivCompP
   P%startGeoORder = 0
   P%endGeoOrder = Input%geoDerOrderP
   P%nCartesianMomentComp=1
   CMorder = 0
   signP=1.0E0_realk
   P%DO_MULMOM  = .FALSE.
   P%MMORDER    = 0
   P%magderiv = INPUT%MagderOrderP
ELSE
   ndmat = INPUT%NDMAT_RHS
   IF(INPUT%PS_SCREEN.OR.INPUT%CS_SCREEN.OR.Input%MBIE_SCREEN)THEN
      GAB2 => INPUT%LST_GAB_RHS   
   ENDIF
   IF(INPUT%PS_SCREEN)THEN
      maxPrimGabElm=INPUT%PS_MAXELM_LHS
   ENDIF
   IF(INPUT%CS_SCREEN)THEN
      IF (Input%DO_JENGINE) THEN
        DMATscreen = .TRUE.
        DMAT2 => INPUT%LST_DRHS
        DMATpermute = INPUT%sameRHSaos
      ENDIF
      IF (Input%DO_DACOULOMB) THEN
        DMATscreen = .TRUE.
        DMAT2 => INPUT%LST_DRHS
      ENDIF
   ENDIF
   LHS=.FALSE.
   ndim5=INPUT%nCartesianMomentComp*INPUT%nMagDerivCompQ
   P%startGeoORder = 0
   P%endGeoOrder = input%geoDerOrderQ
   P%nCartesianMomentComp=INPUT%nCartesianMomentComp
   CMorder = INPUT%CMorder
   signP=-1.0E0_realk
   IF(INPUT%OPERATOR .EQ.MulmomOperator) THEN
      P%DO_MULMOM  = .TRUE.
      P%MMORDER    = INPUT%MMORDER
      NMOM         = INPUT%nMultipoleMomentComp
   ELSE
      P%DO_MULMOM  = .FALSE.
      P%MMORDER    = 0
   ENDIF
   P%magderiv = INPUT%MagderOrderQ
ENDIF
IF(DoAlloc)THEN
 CALL ALLOC_ORBITAL(P%orbital1,ODB%AO(1)%p%nprimitives,1,ODB%AO(1)%p%nAngmom,ODB%AO(1)%p%maxAngmom,.FALSE.,IELECTRON)
 CALL ALLOC_ORBITAL(P%orbital2,ODB%AO(2)%p%nprimitives,1,ODB%AO(2)%p%nAngmom,ODB%AO(2)%p%maxAngmom,.FALSE.,IELECTRON)
 call mem_ODpointer_alloc(P%Orb1atom,1,IELECTRON)
 call mem_ODpointer_alloc(P%Orb1batch,1,IELECTRON)
 call mem_ODpointer_alloc(P%Orb1mol,1,IELECTRON)
 call mem_ODpointer_alloc(P%Orb2atom,1,IELECTRON)
 call mem_ODpointer_alloc(P%Orb2batch,1,IELECTRON)
 call mem_ODpointer_alloc(P%Orb2mol,1,IELECTRON)
ENDIF

CALL SET_ORBITAL(P%orbital1,ODB%AO(1)%p,integral,input%HermiteEcoeff,LUPRI)
CALL SET_ORBITAL(P%orbital2,ODB%AO(2)%p,integral,input%HermiteEcoeff,LUPRI)
P%Orb1atom(1) = ODB%AO(1)%p%atom
P%Orb2atom(1) = ODB%AO(2)%p%atom
P%Orb1mol(1) = ODB%AO(1)%p%molecularIndex
P%Orb2mol(1) = ODB%AO(2)%p%molecularIndex
IF(.NOT.(Input%CS_int).OR.Input%RealGabMatrix)THEN !default
   P%Orb1batch(1) = ODB%AO(1)%p%batch
   P%Orb2batch(1) = ODB%AO(2)%p%batch
ELSE !Screening integrals 
   P%Orb1batch(1) = ODB%IA 
   P%Orb2batch(1) = ODB%IB
ENDIF
!==========================================================================
!                     Determine overlap type
!==========================================================================
P%type_Empty = P%orbital1%type_Empty .AND. P%orbital2%type_Empty 
P%maxGab = shortzero
IF(P%type_Empty)THEN
   IF(Input%operator .EQ. MulmomOperator .AND. (.NOT.LHS) )THEN
      ! this is done in order to use screening on multipole moments used for FMM
      P%maxGab = shortzero
   ENDIF
ENDIF
P%type_Hermite_single = (P%orbital1%type_Hermite .AND. P%orbital2%type_Empty).OR.&
     &(P%orbital1%type_Empty .AND. P%orbital2%type_Hermite)
P%type_Hermite = P%orbital1%type_Hermite .AND. P%orbital2%type_Hermite
P%type_Cartesian_single = (P%orbital1%type_Cartesian .AND. P%orbital2%type_Empty).OR.&
     &(P%orbital1%type_Empty .AND. P%orbital2%type_Cartesian)
P%type_Cartesian = P%orbital1%type_Cartesian .AND. P%orbital2%type_Cartesian
P%type_Nucleus = (P%orbital1%type_Nucleus .AND. P%orbital2%type_Empty)&
     &.OR.(P%orbital1%type_Empty .AND. P%orbital2%type_Nucleus)
P%type_FTUV = .FALSE.

P%single = P%type_Hermite_single.OR.P%type_Cartesian_single.OR.P%type_Nucleus

useFTUV = INPUT%DO_JENGINE .AND. IELECTRON.EQ. 2
! Spherical E-coefficients reduces computational scaling wrt. angular momentum, and 
! should be used as much as possible. However, inthe below cases it should be turned off.
! Note here that for J-engine (useFTUV) we have a special case
 P%sphericalEcoeff = (Input%sphericalEcoeff.AND..NOT.P%type_hermite_single).OR.useFTUV
 P%sphericalEcoeff = P%sphericalEcoeff.AND.(P%orbital1%spherical.OR.P%orbital2%spherical)
 P%sphericalEcoeff = P%sphericalEcoeff.AND..NOT.(INPUT%operator .EQ. KineticOperator .AND. LHS) 
!Note that the next line will cause difficulties for jengine derivatives (for RHS differentiation)
 P%sphericalEcoeff = P%sphericalEcoeff.AND.(.NOT.P%endGeoOrder.GE. 1)

!==========================================================================

 CALL getDerivComp(P%nGeoDerivComp,P%endGeoOrder,P%type_empty,P%single)
 

 P%ODcenter = ODB%ODcenter
 P%ODextent = ODB%ODextent
 P%sameAO = ODB%sameAO
 CALL GET_NPRIMITIVES(np,P,Input,ielectron,lupri)
 maxnp = max(np,P%orbital1%nPrimitives,P%orbital2%nPrimitives)
 P%nPrimitives   = np
! Default is to build OD-batches first, and then later collect into passes
 P%nPasses       = 1
#ifdef VAR_LSDEBUGINT
 IF (np.EQ. 0) call lsquit('np=0. this should never happen! TK',lupri)
#endif
!==========================================================================
!   Determine angular momentum to alloc the Overlap structure
!             and some angular momentum related quantities
!==========================================================================

 P%nAngmom       = ODB%nAngmom
 P%maxContracted = ODB%maxContracted
!
 indexETUV = 1
 minAngmom   = 99
 maxangmom = 0
 i12 = 0
 maxijk = 0
 nOperatorComp = 1
 IF(P%Magderiv.EQ.1)nOperatorComp = 3
 CALL MEM_ALLOC(lenETUV,P%orbital1%nAngmom*P%orbital2%nAngmom)
 DO i1=1,P%orbital1%nAngmom
    start2 = 1
    IF (P%sameAO) start2 = i1
    DO i2=start2,P%orbital2%nAngmom
       i12  = i12 + 1
       l1   = P%orbital1%angmom(i1)
       l2   = P%orbital2%angmom(i2)
       CALL GET_IJK(l1,l2,ijk1,ijk2,ijk,ijkcart,P%sphericalEcoeff,P%endGeoOrder,P%single)
       maxijk = max(maxijk,ijk)*nOperatorComp
       l  = l1+l2+P%endGeoOrder+CMorder+P%magderiv
       nTUV = (l+1)*(l+2)*(l+3)/6 !including nEFG if any
       lenETUV(i12) =  nTUV*ijk*P%nPrimitives*nOperatorComp
       indexETUV = indexETUV + nTUV*ijk*P%nPrimitives
       maxangmom = MAX(maxangmom,l)
       minangmom = MIN(minangmom,l)
    ENDDO
 ENDDO

 nETUV = indexETUV - 1
 P%maxAngmom = maxAngmom-(P%endGeoOrder+CMorder+P%magderiv)
 P%minAngmom = minAngmom-(P%endGeoOrder+CMorder+P%magderiv)
 P%startAngmom = 0
 IF (P%type_hermite_single) P%startAngmom = P%minAngmom
 P%endAngmom   = maxAngmom

 P%nTUV = (P%endAngmom+1)*(P%endAngmom+2)*(P%endAngmom+3)/6 - &
    &     P%startAngmom*(P%startAngmom+1)*(P%startAngmom+2)/6
!==========================================================================
!  Alloc the Overlap if specified 
!==========================================================================
IF(DoAlloc)THEN
 CALL ALLOC_OVERLAP(P,np,1,P%nAngmom,useFTUV,P%nTUV,P%endGeoOrder,Input%NDMAT_RHS,ielectron)
ENDIF
!==========================================================================
!  Fill in Angmom quantites in the now allocd overlap 
!==========================================================================

DO iAngmom = 1,P%nAngmom
   P%lenETUV(iAngmom) = lenETUV(iAngmom)  
ENDDO
CALL MEM_DEALLOC(lenETUV)

P%segmented = .TRUE.
P%totOrbitals(:) = 0
i12 = 0
DO i1=1,P%orbital1%nAngmom
   start2 = 1
   IF (P%sameAO) start2 = i1
   DO i2=start2,P%orbital2%nAngmom
      i12  = i12 + 1
      P%angmom(i12)      = P%orbital1%angmom(i1) + P%orbital2%angmom(i2)
      P%indexAng1(i12)   = i1
      P%indexAng2(i12)   = i2
      P%nOrbitals(i12)   = P%orbital1%nOrbitals(i1)*P%orbital2%nOrbitals(i2)
      DO iDer=P%startGeoOrder,P%endGeoOrder
        P%totOrbitals(iDer+1) = P%totOrbitals(iDer+1) + P%nOrbitals(i12)*ndim5*NMOM*getODgeoComp(iDer,P%single)
      ENDDO
      P%nOrbComp(i12)    = P%orbital1%nOrbComp(i1)*P%orbital2%nOrbComp(i2)
      P%nContracted(i12) = P%orbital1%nContracted(i1)*P%orbital2%nContracted(i2)
      IF (P%nContracted(i12).GT. 1.OR.(INPUT%DO_JENGINE.AND.IELECTRON.EQ. 2)) P%segmented = .FALSE.
   ENDDO
ENDDO
!CAREFUL
!IF (P%segmented) THEN
!  IF (IELECTRON.EQ. 1) P%segmented = .FALSE.
!ENDIF
!CAREFUL

!==========================================================================
!   Determine P%distance12 P%squaredDistance
!==========================================================================

 d2 = 0.0E0_realk
 DO idir=1,3
   P%distance12(idir) = P%orbital1%center(idir)-P%orbital2%center(idir)
   d2 = d2 + P%distance12(idir) * P%distance12(idir)
 ENDDO
 P%squaredDistance = d2

!==========================================================================
!   Build primgabmat
!==========================================================================

 doscreen = ((INPUT%PS_SCREEN).AND.(.NOT.P%type_Empty))
 IF(doscreen)then
  !change order! 
  call mem_alloc(primgabmat,P%orbital1%nPrimitives,P%orbital2%nPrimitives)
  AtomA = P%orb1atom(1)
  AtomB = P%orb2atom(1)
  BatchA = P%orb1batch(1)
  BatchB = P%orb2batch(1)
  Gindex = GAB2%INDEX(AtomA,AtomB,1,1)
#ifdef VAR_LSDEBUGINT
  IF(Gindex.EQ.0)THEN
     call lsquit('lstensor allocation error in SET_OVERLAP',lupri)
  ENDIF
#endif
  maxBat = GAB2%SLSAO(Gindex)%maxBat
  s1 = GAB2%SLSAO(Gindex)%startLocalOrb(BatchA)-1
  s2 = GAB2%SLSAO(Gindex)%startLocalOrb(BatchB+maxBat)-1
  n1 = GAB2%SLSAO(Gindex)%nLocal(1)
  DO i2=1,P%orbital2%nPrimitives
   DO i1=1,P%orbital1%nPrimitives
    primgabmat(i1,i2) = &
         & (GAB2%SLSAO(Gindex)%selms(s1+i1 + (s2+i2-1)*n1))
   ENDDO
  ENDDO
 ENDIF
!==========================================================================
!   Determine P%exponents, P%reducedExponents, P%center, P%preExpFac 
!             P%iprim1, P%iprim2, using primgabmat
!==========================================================================
 i12 = 0
 DO i2=1,P%orbital2%nPrimitives
    DO i1=1,P%orbital1%nPrimitives
     IF (doscreen) THEN
       screen = primgabmat(i1,i2) .LE. INPUT%PS_THRLOG-maxPrimGabElm
     ELSE
       screen = .FALSE.
     ENDIF
     IF (.NOT.screen) THEN
#ifdef VAR_LSDEBUGINT
       IF(size(P%preExpFac).NE.P%nPrimAlloc*P%nAngAlloc)call lsquit('PreExpFac dim error',-1)
#endif
       offset = i12*3
       i12 = i12 + 1
       e1  = P%orbital1%exponents(i1)
       e2  = P%orbital2%exponents(i2)       
       P%exponents(i12)        = e1 + e2
       IF (P%exponents(i12) .NE. 0.0E0_realk) THEN
         P%reducedExponents(i12) = e1*e2/(e1+e2)
         P%center(1+offset) = (e1*P%orbital1%center(1) + e2*P%orbital2%center(1))/(e1+e2)
         P%center(2+offset) = (e1*P%orbital1%center(2) + e2*P%orbital2%center(2))/(e1+e2)
         P%center(3+offset) = (e1*P%orbital1%center(3) + e2*P%orbital2%center(3))/(e1+e2)
       ELSEIF(Input%operator .EQ. NucpotOperator)THEN
         P%reducedExponents(i12) = e1
         P%center(1+offset)         = P%orbital1%center(1)
         P%center(2+offset)         = P%orbital1%center(2)
         P%center(3+offset)         = P%orbital1%center(3)
       ELSE
         P%reducedExponents(i12) = 0.0E0_realk
         P%center(1+offset)         = 0.0E0_realk
         P%center(2+offset)         = 0.0E0_realk
         P%center(3+offset)         = 0.0E0_realk
       ENDIF 
       IF (.NOT. P%single) THEN 
         IF (ABS(P%exponents(i12)) .GT. 1.0E-15_realk) THEN
!PGI compiler gives some stange values for P%reducedExponents(i12)*d2 = 728 
            IF (P%reducedExponents(i12)*d2 .GT. 500.0E0_realk) THEN
               P%preExpFac(i12)        = 0.0E0_realk
            ELSE	    
               P%preExpFac(i12)        = exp(-P%reducedExponents(i12)*d2)
            ENDIF	    
         ELSEIF(Input%operator .EQ. NucpotOperator)THEN
            P%preExpFac(i12)        = exp(-e1*d2)
         ELSE
            P%preExpFac(i12)        = 1.0E0_realk
         ENDIF
       ELSE
          P%preExpFac(i12)        = 1.0E0_realk
       ENDIF
       P%iprim1(i12)           = i1
       P%iprim2(i12)           = i2
     ENDIF   
   ENDDO
 ENDDO
 IF (np .NE. i12) THEN
   WRITE(LUPRI,'(1X,A,2I5)') 'Error in set_overlap. Mismatch in number of primitives, (np,i12)=',&
     &  np,i12
   CALL LSQUIT('Error: Mismatch between get_nPrimitives and i12 in set_overlap',lupri)
 ENDIF
 IF(P%nAngmom.GT.1)THEN
    DO iAngmom = 2,P%nAngmom
       offset = (iAngmom-1)*P%nPrimAlloc
       DO i12 = 1,np
          P%preExpFac(i12+offset) = P%preExpFac(i12) 
       ENDDO
    ENDDO
 ENDIF
IF(doscreen) THEN
   call mem_dealloc(primgabmat)
ENDIF
!========================================================================== 
!   Determine if we should perform a basisset contraction, this should be done 
!   before Etensor is calculated (before P%preExpFac(1) is used)  
!========================================================================== 
#ifdef VAR_LSDEBUGINT
IF(size(P%preExpFac).NE.P%nPrimAlloc*P%nAngAlloc)call lsquit('PreExpFac dim error2',-1)
#endif
IF(Input%PS_int)THEN
 P%ContractBasis = .FALSE. 
 nprim1=P%orbital1%nPrimitives
 nprim2=P%orbital2%nPrimitives
 iA12 = 0
 DO iA1=1,P%orbital1%nAngmom
  startA2 = 1
  IF (P%sameAO) startA2 = iA1
  DO iA2=startA2,P%orbital2%nAngmom
   iA12  = iA12 + 1
   offset = (iA12-1)*P%nPrimAlloc
   DO i12=1,P%nPrimitives
     i1 = P%iprim1(i12)
     i2 = P%iprim2(i12)
     P%preExpFac(i12+offset) = P%preExpFac(i12+offset)&
          &*P%orbital1%CC(iA1)%p%elms(i1+(i1-1)*nprim1)&
          &*P%orbital2%CC(iA2)%p%elms(i2+(i2-1)*nprim2) 
   ENDDO
  ENDDO
 ENDDO
ELSEIF (P%segmented) THEN
 P%ContractBasis = .FALSE. 
 nprim1=P%orbital1%nPrimitives
 nprim2=P%orbital2%nPrimitives
 iA12 = 0
 DO iA1=1,P%orbital1%nAngmom
  startA2 = 1
  IF (P%sameAO) startA2 = iA1
  DO iA2=startA2,P%orbital2%nAngmom
   offset = iA12*P%nPrimAlloc
   iA12  = iA12 + 1
   DO i12=1,P%nPrimitives
     i1 = P%iprim1(i12)
     i2 = P%iprim2(i12)
     P%preExpFac(i12+offset) = P%preExpFac(i12+offset)*&
          & P%orbital1%CC(iA1)%p%elms(i1)*P%orbital2%CC(iA2)%p%elms(i2) 
   ENDDO
  ENDDO
 ENDDO
ELSE 
   P%ContractBasis = .TRUE. 
ENDIF

!==========================================================================
!   Determine Screening Quantities : P%maxgab and P%MBIE(1), P%MBIE(2)
!==========================================================================

IF (Input%CS_SCREEN .AND. (.NOT.P%type_Empty)) THEN
   maxGab = shortzero
   IF (DMATscreen)THEN
      maxGab = shortzero
      AtomA = P%orb1atom(1)
      BatchA = P%orb1batch(1)
      AtomB = P%orb2atom(1)
      BatchB = P%orb2batch(1)
      Dindex = DMAT2%INDEX(atomA,atomB,1,1)
      maxBat = DMAT2%LSAO(Dindex)%maxBat
      maxAng = DMAT2%LSAO(Dindex)%maxAng
      offset1 = (batchA-1)*maxAng
      offset2 = (batchB-1)*maxAng+maxAng*maxBat
      IF(Dindex.NE.0)THEN
         DO i1=1,P%orbital1%nAngmom
            start = 1
            IF (P%sameAO) start = i1
            DO i2=start,P%orbital2%nAngmom
               iSA = DMAT2%LSAO(Dindex)%startLocalOrb(i1+offset1)-1
               iSB = DMAT2%LSAO(Dindex)%startLocalOrb(i2+offset2)-1
               dimA = DMAT2%LSAO(Dindex)%nLocal(1)
               dimB = DMAT2%LSAO(Dindex)%nLocal(2)
               DO orb1=1,P%orbital1%nOrbitals(i1)
                  DO orb2=1,P%orbital2%nOrbitals(i2)
                     DMATmax = 0.0E0_realk                  
                     DO idmat=1,ndmat
                        elms = iSA+orb1 + (iSB+orb2-1)*dimA + (idmat-1)*dimA*dimB
                        DMATmax = max(DMATmax,abs(DMAT2%LSAO(Dindex)%elms(elms)))
                     ENDDO
                     !Beware when converting from double precision to short integer
                     !If double precision is less than 10^-33 then you can run into
                     !problems with short integer overflow
                     IF(Dmatmax.GT.ShortIntCrit)THEN
                        DMATMAX2 = CEILING(log10(DMATmax))
                     ELSE
                        DMATMAX2 = shortzero
                     ENDIF
                     maxgabelm = ODB%maxgab + DMATmax2
                     maxGab = MAX(maxGab,maxgabelm)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ELSE
         maxGab = shortzero
      ENDIF
      IF(DMATpermute)THEN
         AtomA = P%orb1atom(1)
         BatchA = P%orb1batch(1)
         AtomB = P%orb2atom(1)
         BatchB = P%orb2batch(1)         
         Dindex = DMAT2%INDEX(atomB,atomA,1,1)
         maxBat = DMAT2%LSAO(Dindex)%maxBat
         maxAng = DMAT2%LSAO(Dindex)%maxAng
         offset1 = (batchB-1)*maxAng                !now offset for B
         offset2 = (batchA-1)*maxAng+maxAng*maxBat  !now offset for A
         IF(Dindex.NE.0)THEN
            DO i1=1,P%orbital1%nAngmom    !angmom for A
               start = 1
               IF (P%sameAO) start = i1
               DO i2=start,P%orbital2%nAngmom 
                  IF(P%orbital1%startOrbital(i1).NE.P%orbital2%startOrbital(i2))THEN
                     iSB = DMAT2%LSAO(Dindex)%startLocalOrb(i2+offset1)-1  !B
                     iSA = DMAT2%LSAO(Dindex)%startLocalOrb(i1+offset2)-1  !A
                     dimB = DMAT2%LSAO(Dindex)%nLocal(1)
                     dimA = DMAT2%LSAO(Dindex)%nLocal(2)
                     DO orb1=1,P%orbital1%nOrbitals(i1)      !A
                        DO orb2=1,P%orbital2%nOrbitals(i2)   !B
                           DO idmat=1,ndmat
                              elms = iSB+orb2 + (iSA+orb1-1)*dimB + (idmat-1)*dimA*dimB
                              DMATmax = max(DMATmax,abs(DMAT2%LSAO(Dindex)%elms(elms)))
                           ENDDO
                           !Beware when converting from double precision to short integer
                           !If double precision is less than 10^-33 then you can run into
                           !problems with short integer overflow
                           IF(Dmatmax.GT.ShortIntCrit)THEN
                              DMATMAX2 = CEILING(log10(DMATmax))
                           ELSE
                              DMATMAX2 = shortzero
                           ENDIF
                           maxgabelm = ODB%maxgab + DMATmax2
                           maxGab = MAX(maxGab,maxgabelm)
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
      ENDIF
   ELSE
      maxGab = ODB%maxgab
   ENDIF
   P%maxGab = maxGab
ENDIF
IF (Input%MBIE_SCREEN .AND. (.NOT.P%type_Empty)) THEN
   P%MBIE(1) = GAB2%MBIE(1,ODB%IA,ODB%IB)
   P%MBIE(2) = GAB2%MBIE(2,ODB%IA,ODB%IB)
ENDIF
NULLIFY(GAB2)

!==========================================================================
!   Jengine specific
!==========================================================================

IF (INPUT%DO_JENGINE .AND. IELECTRON.EQ. 2) THEN
  IF(.NOT.DoAlloc)call lsquit('set_overlap jengine error',lupri)
  CALL SET_FTUV(P,INPUT,Integral,LUPRI,IPRINT)
! Change Overlap P to confirm with FTUV
  P%TYPE_FTUV = .TRUE.
  P%TYPE_Empty = .FALSE.
  P%TYPE_Hermite = .FALSE.
  P%TYPE_Hermite_single = .FALSE.
  P%TYPE_Cartesian = .FALSE.
  P%TYPE_Cartesian_single = .FALSE.
  P%TYPE_Nucleus = .FALSE.
  P%sphericalEcoeff = .FALSE.
  P%nAngmom = 1
  P%angmom(1) = P%maxAngmom
  P%orbital1%nAngmom = 1
  P%orbital1%angmom(1) = 0
  P%orbital1%nContracted(1)    = 1
  P%orbital1%nOrbComp(1)       = 1
  P%orbital1%startOrbital(1) = 1
  P%orbital1%totOrbitals       = ndmat
  P%orbital2%nAngmom = 1
  P%orbital2%angmom(1) = 0
  P%orbital2%nContracted(1)    = 1
  P%orbital2%nOrbComp(1)       = 1
  P%orbital2%startOrbital(1) = 1
  P%orbital2%totOrbitals       = 1
  P%indexAng1(1)   = 1
  P%indexAng2(1)   = 1
  P%nOrbitals(1)   = 1
  P%totOrbitals(:) = ndmat
  P%nGeoDerivComp  = 1
  P%nCartesianMomentComp  = 1
  P%nOrbComp(1)    = 1
  P%nContracted(1) = 1
ENDIF

END SUBROUTINE SET_OVERLAP

SUBROUTINE MEM_OVERLAP(OD,Input,IELECTRON)
Implicit none
TYPE(ODITEM) :: OD
TYPE(IntegralInput),target,intent(in) :: Input
Integer,intent(in)                    :: IELECTRON
!
logical :: useFTUV
integer :: extraAngmom,IRHS,maxAngmom,nTUV,np1,np2,ng
integer(kind=long) :: nint,nrealk,nint2,nrealk2
IF( INPUT%operator .EQ. KineticOperator .AND. IELECTRON.EQ.1)THEN
   extraAngmom=2
ELSE
   extraAngmom=0
ENDIF
useFTUV = INPUT%DO_JENGINE .AND. IELECTRON.EQ. 2
nint2 = 0
nrealk2 = 0
DO IRHS=1,OD%nbatches
   nint = 0
   nrealk = 0
   np1 = OD%BATCH(IRHS)%AO(1)%p%nprimitives
   np2 = OD%BATCH(IRHS)%AO(2)%p%nprimitives
   CALL MEM_ORBITAL(np1,1,OD%BATCH(IRHS)%AO(1)%p%nAngmom,&
        & OD%BATCH(IRHS)%AO(1)%p%maxAngmom,.FALSE.,nrealk,nint)
   CALL MEM_ORBITAL(np2,1,OD%BATCH(IRHS)%AO(2)%p%nAngmom,&
        & OD%BATCH(IRHS)%AO(2)%p%maxAngmom,.FALSE.,nrealk,nint)
   nint = nint + 6
   maxangmom = OD%BATCH(IRHS)%AO(1)%p%maxAngmom + OD%BATCH(IRHS)%AO(2)%p%maxAngmom + extraAngmom
   nTUV = (maxAngmom+1)*(maxAngmom+2)*(maxAngmom+3)/6
   maxAngmom = OD%BATCH(IRHS)%nAngmom
   IF (IELECTRON.EQ.1) THEN
     ng = input%geoDerOrderP
   ELSE
     ng = input%geoDerOrderQ
   ENDIF
   CALL MEM_OVERLAP1(np1*np2,1,maxAngmom,useFTUV,nTUV,ng,Input%NDMAT_RHS,nrealk,nint)
   nint2 = nint2+nint
   nrealk2 = nrealk2 + nrealk   
ENDDO
call SET_BUFCOUNTERS(ielectron,nint2,nrealk2)
END SUBROUTINE MEM_OVERLAP

SUBROUTINE MEM_SINGLE_OVERLAP(OD,IRHS,Input,IELECTRON)
Implicit none
TYPE(ODITEM) :: OD
TYPE(IntegralInput),target,intent(in) :: Input
Integer,intent(in)                    :: IELECTRON,IRHS
!
logical :: useFTUV
integer :: extraAngmom,maxAngmom,nTUV,np1,np2,nAngmom,ng
integer(kind=long) :: nint,nrealk,nint2,nrealk2
IF( INPUT%operator .EQ. KineticOperator .AND. IELECTRON.EQ.1)THEN
   extraAngmom=2
ELSE
   extraAngmom=0
ENDIF
useFTUV = INPUT%DO_JENGINE .AND. IELECTRON.EQ. 2
nint = 0
nrealk = 0

IF (IELECTRON.EQ.1) THEN
  ng = input%geoDerOrderP
ELSE
  ng = input%geoDerOrderQ
ENDIF

np1 = OD%BATCH(IRHS)%AO(1)%p%nprimitives
np2 = OD%BATCH(IRHS)%AO(2)%p%nprimitives
CALL MEM_ORBITAL(np1,1,OD%BATCH(IRHS)%AO(1)%p%nAngmom,&
     & OD%BATCH(IRHS)%AO(1)%p%maxAngmom,.FALSE.,nrealk,nint)
CALL MEM_ORBITAL(np2,1,OD%BATCH(IRHS)%AO(2)%p%nAngmom,&
     & OD%BATCH(IRHS)%AO(2)%p%maxAngmom,.FALSE.,nrealk,nint)
nint = nint + 6
maxangmom = OD%BATCH(IRHS)%AO(1)%p%maxAngmom + OD%BATCH(IRHS)%AO(2)%p%maxAngmom + extraAngmom
nTUV = (maxAngmom+1)*(maxAngmom+2)*(maxAngmom+3)/6
nAngmom = OD%BATCH(IRHS)%nAngmom
CALL MEM_OVERLAP1(np1*np2,1,nAngmom,useFTUV,nTUV,ng,Input%NDMAT_RHS,nrealk,nint)
call ADD_BUFCOUNTERS(ielectron,nint,nrealk)

END SUBROUTINE MEM_SINGLE_OVERLAP

SUBROUTINE MEM_ORBITAL(maxprim,npass,maxangmom,maxA,FTUVorb,nrealk,nint)
IMPLICIT NONE
Integer,intent(in) :: maxprim,maxangmom,maxA,npass
integer(kind=long),intent(inout) :: nrealk,nint
logical,intent(in) :: FTUVorb
nint = nint + 4*maxangmom + maxangmom*npass
IF(.NOT.FTUVorb)THEN
   nrealk = nrealk + maxprim
   nint = nint + maxangmom
ENDIF
END SUBROUTINE MEM_ORBITAL

subroutine MEM_OVERLAP1(nprim12,npass,nAng,allocFTUV,nTUV,geoOrder,ndmat,nrealk,nint)
implicit none
integer,intent(in) :: nprim12,npass,nAng,nTUV,ndmat,geoOrder
logical,intent(in) :: allocFTUV
integer(kind=long),intent(inout) :: nrealk,nint
!
integer(kind=long) :: tmp
nrealk = nrealk  + (5+nAng)*nprim12+3*npass
nint = nint + 2*nprim12 + 7*nAng + (geoOrder+1)
IF (allocFTUV) THEN
   tmp = nTUV*ndmat
   nrealk = nrealk  +  nPrim12*tmp
ENDIF
end subroutine MEM_OVERLAP1

!> \brief Returns the number of two-center derivative components for a given geometric derivative order
!> \author S. Reine
!> \date 2010-03-05
!> \param nDerivComp The number of derivative components
!> \param derOrder The derivative order
!> \param empty True if empty overlap
!> \param single True if only a single AO
SUBROUTINE getDerivComp(nDerivComp,derOrder,empty,single)
implicit none
Integer,intent(OUT) :: nDerivComp
Integer,intent(IN)  :: derOrder
Logical,intent(IN)  :: empty,single
!
Integer :: iDer,jDer
!
!Consistency testing
IF (empty.AND.single) CALL LSQUIT('Error in getDerivComp. Single and empty!',-1)
IF (derOrder.LT. 0) CALL LSQUIT('Error in getDerivComp. derOrder<0',-1)
!
IF (empty) THEN
  nDerivComp = 1
ELSE IF (single) THEN
  nDerivComp = (derOrder+1)*(derOrder+2)/2
ELSE 
  nDerivComp = 0
  DO iDer=0,derOrder
    jDer = derOrder-iDer
    nDerivComp = nDerivComp + (iDer+1)*(iDer+2)*(jDer+1)*(jDer+2)/4
  ENDDO
ENDIF
!
END SUBROUTINE getDerivComp

!> \brief initialise the overlap
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param P The overlap
!> \param Alloc Information about the maximal sizes for allocation purposes
!> \param Input contains info about the integral evaluation requested
!> \param ODbat The OD-batch information
!> \param side The side ('LHS' or 'RHS') 
!> \param iElectron The electron (1 = 'LHS', 2 = 'RHS')
!> \param iprint The print level (the higher the more information)
!> \param lupri The default print unit
SUBROUTINE INIT_OVERLAP(P,Alloc,iAlloc,Input,IELECTRON,IPRINT,LUPRI)
Implicit none
INTEGER             :: lupri,IPRINT,IELECTRON,iAlloc
TYPE(Overlap)       :: P
TYPE(IntegralInput) :: Input
TYPE(AllocItem)     :: Alloc
Integer             :: maxprim,maxtuv,maxnangmom,geoOrder
Integer             :: MAXETUVLEN,maxA
LOGICAL           :: LHS,hermiteSingle,useFTUV 

IF(IELECTRON.EQ.1)THEN
   LHS=.TRUE.
ELSEIF(IELECTRON.EQ.2)THEN
   LHS=.FALSE.
ELSE
   WRITE(LUPRI,'(1X,2A)') 'Wrong case in INIT_OVERLAP =',IELECTRON
   CALL LSQUIT('Wrong case in INIT_OVERLAP',lupri)
END IF

IF(LHS)THEN
   IF(iAlloc.NE. 0)THEN 
      maxprim = Alloc%maxPrimLHSA(iAlloc) 
      maxnAngmom = Alloc%maxnAngLHSA(iAlloc) 
      maxA = Alloc%maxAngmomOrbLHSA(iAlloc) 
      maxTUV = Alloc%maxTUVLHSA(iAlloc) 
      MAXETUVLEN = Alloc%maxetuvlenLHSA(iAlloc) 
   ELSE 
      maxprim = Alloc%maxPrimLHS 
      maxnAngmom = Alloc%maxnAngLHS 
      maxA = Alloc%maxAngmomOrbLHS 
      maxTUV = Alloc%maxTUVLHS 
      MAXETUVLEN = Alloc%maxetuvlenLHS 
   ENDIF
   IF(maxnAngmom .EQ. 0)maxnAngmom=1
   IF(MAXETUVLEN .EQ. 0)MAXETUVLEN=1
   geoOrder = input%geoDerOrderP
ELSE
   IF(iAlloc.NE. 0)THEN 
      maxprim = Alloc%maxPrimRHSA(iAlloc) 
      maxnAngmom = Alloc%maxnAngRHSA(iAlloc) 
      maxA = Alloc%maxAngmomOrbRHSA(iAlloc) 
      maxTUV = Alloc%maxTUVRHSA(iAlloc) 
      MAXETUVLEN = Alloc%maxetuvlenRHSA(iAlloc) 
   ELSE 
      maxprim = Alloc%maxPrimRHS 
      maxnAngmom = Alloc%maxnAngRHS 
      maxA = Alloc%maxAngmomOrbRHS 
      maxTUV = Alloc%maxTUVRHS 
      MAXETUVLEN = Alloc%maxetuvlenRHS 
   ENDIF
   IF(maxnAngmom .EQ. 0)maxnAngmom=1
   IF(MAXETUVLEN .EQ. 0)MAXETUVLEN=1
   geoOrder = input%geoDerOrderQ
ENDIF

CALL ALLOC_ORBITAL(P%orbital1,maxprim,1,maxnangmom,maxA,.FALSE.,IELECTRON)
CALL ALLOC_ORBITAL(P%orbital2,maxprim,1,maxnangmom,maxA,.FALSE.,IELECTRON)
call mem_ODpointer_alloc(P%Orb1atom,1,ielectron)
call mem_ODpointer_alloc(P%Orb1mol,1,ielectron)
call mem_ODpointer_alloc(P%Orb1batch,1,ielectron)
call mem_ODpointer_alloc(P%Orb2atom,1,ielectron)
call mem_ODpointer_alloc(P%Orb2mol,1,ielectron)
call mem_ODpointer_alloc(P%Orb2batch,1,ielectron)

! hermiteSingle = input%hermiteEcoeff.AND.&
!      & ((.NOT.ODbat%BATCH(1)%AO(1)%p%type_Empty.AND.ODbat%BATCH(1)%AO(2)%p%type_Empty).OR.&
!      & (ODbat%BATCH(1)%AO(1)%p%type_Empty.AND..NOT.ODbat%BATCH(1)%AO(2)%p%type_Empty))
 useFTUV = INPUT%DO_JENGINE .AND. IELECTRON.EQ. 2
 CALL ALLOC_OVERLAP(P,maxprim,1,maxnAngmom,useFTUV,maxTUV,geoOrder,Input%NDMAT_RHS,ielectron)

END SUBROUTINE INIT_OVERLAP

!> \brief initialise the overlap
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param P The overlap
!> \param Alloc Information about the maximal sizes for allocation purposes
!> \param Input contains info about the integral evaluation requested
!> \param ODbat The OD-batch information
!> \param side The side ('LHS' or 'RHS') 
!> \param iElectron The electron (1 = 'LHS', 2 = 'RHS')
!> \param iprint The print level (the higher the more information)
!> \param lupri The default print unit
SUBROUTINE MEM_INIT_OVERLAP(Alloc,iAlloc,Input,IELECTRON,IPRINT,LUPRI)
Implicit none
INTEGER             :: lupri,IPRINT,IELECTRON,iAlloc
TYPE(IntegralInput) :: Input
TYPE(AllocItem)     :: Alloc
Integer             :: maxprim,maxtuv,maxnangmom,maxgeoorder
Integer             :: MAXETUVLEN,maxA
LOGICAL           :: LHS,hermiteSingle,useFTUV
integer(kind=long) :: nint,nrealk
nrealk=0
nint=0
IF(IELECTRON.EQ.1)THEN
   LHS=.TRUE.
ELSEIF(IELECTRON.EQ.2)THEN
   LHS=.FALSE.
ELSE
   WRITE(LUPRI,'(1X,2A)') 'Wrong case in MEM_INIT_OVERLAP =',IELECTRON
   CALL LSQUIT('Wrong case in MEM_INIT_OVERLAP',lupri)
END IF

IF(LHS)THEN
   IF(iAlloc.NE. 0)THEN 
      maxprim = Alloc%maxPrimLHSA(iAlloc) 
      maxnAngmom = Alloc%maxnAngLHSA(iAlloc) 
      maxA = Alloc%maxAngmomOrbLHSA(iAlloc) 
      maxTUV = Alloc%maxTUVLHSA(iAlloc) 
      MAXETUVLEN = Alloc%maxetuvlenLHSA(iAlloc) 
   ELSE 
      maxprim = Alloc%maxPrimLHS 
      maxnAngmom = Alloc%maxnAngLHS 
      maxA = Alloc%maxAngmomOrbLHS 
      maxTUV = Alloc%maxTUVLHS 
      MAXETUVLEN = Alloc%maxetuvlenLHS 
   ENDIF
   IF(maxnAngmom .EQ. 0)maxnAngmom=1
   IF(MAXETUVLEN .EQ. 0)MAXETUVLEN=1
   maxgeoorder = input%geoDerOrderP
ELSE
   IF(iAlloc.NE. 0)THEN 
      maxprim = Alloc%maxPrimRHSA(iAlloc) 
      maxnAngmom = Alloc%maxnAngRHSA(iAlloc) 
      maxA = Alloc%maxAngmomOrbRHSA(iAlloc) 
      maxTUV = Alloc%maxTUVRHSA(iAlloc) 
      MAXETUVLEN = Alloc%maxetuvlenRHSA(iAlloc) 
   ELSE 
      maxprim = Alloc%maxPrimRHS 
      maxnAngmom = Alloc%maxnAngRHS 
      maxA = Alloc%maxAngmomOrbRHS 
      maxTUV = Alloc%maxTUVRHS 
      MAXETUVLEN = Alloc%maxetuvlenRHS 
   ENDIF
   IF(maxnAngmom .EQ. 0)maxnAngmom=1
   IF(MAXETUVLEN .EQ. 0)MAXETUVLEN=1
   maxgeoorder = input%geoDerOrderQ
ENDIF

CALL MEM_ORBITAL(maxprim,1,maxnangmom,maxA,.FALSE.,nrealk,nint)
CALL MEM_ORBITAL(maxprim,1,maxnangmom,maxA,.FALSE.,nrealk,nint)
nint = nint + 6
useFTUV = INPUT%DO_JENGINE .AND. IELECTRON.EQ. 2
CALL MEM_OVERLAP1(maxprim,1,maxnAngmom,useFTUV,maxTUV,maxgeoorder,Input%NDMAT_RHS,nrealk,nint)
call MAX_BUFCOUNTERS(ielectron,nint,nrealk)
END SUBROUTINE MEM_INIT_OVERLAP

!> \brief alloc the overlap
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param P The overlap
!> \param nprim1 the number of primitive functions of orbital 1
!> \param nprim2 the number of primitive functions of orbital 2
!> \param npass the number of passes 
!> \param nAng the number of angular moments
!> \param allocFTUV flag to alloc the FTUV
!> \param nTUV the number of angular components 
!> \param geoOrder maximum geometrical differential order
!> \param ndmat the number of density matrices
SUBROUTINE ALLOC_OVERLAP(P,nprim12,npass,nAng,allocFTUV,nTUV,geoOrder,ndmat,ielectron)
use memory_handling
IMPLICIT NONE
TYPE(Overlap)    :: P
INTEGER          :: nprim1,nprim2,npass,nAng,nTUV,geoORder,ndmat,nprim12,ielectron
LOGICAL          :: allocFTUV
integer (kind=long) :: nsize
!!$Call Mem_alloc(P%center,3,nprim12)
!!$Call Mem_alloc(P%exponents,nprim12)        !call mem_alloc(A%Belms,I)
!!$Call Mem_alloc(P%reducedExponents,nprim12)
!!$Call Mem_alloc(P%preExpFac,nprim12,nAng)
!!$Call Mem_alloc(P%iprim1,nprim12)
!!$Call Mem_alloc(P%iprim2,nprim12)
!!$CALL MEM_ALLOC(P%distance12,3,npass)
!!$Call Mem_alloc(P%angmom,nAng)
!!$Call Mem_alloc(P%lenETUV,nAng)
!!$Call Mem_alloc(P%indexAng1,nAng)
!!$Call Mem_alloc(P%indexAng2,nAng)
!!$Call Mem_alloc(P%nOrbitals,nAng)
!!$Call Mem_alloc(P%nOrbComp,nAng)
!!$Call Mem_alloc(P%nContracted,nAng)

Call Mem_ODpointer_alloc(P%center,3*nprim12,ielectron)
Call Mem_ODpointer_alloc(P%exponents,nprim12,ielectron)      
Call Mem_ODpointer_alloc(P%reducedExponents,nprim12,ielectron)
Call Mem_ODpointer_alloc(P%preExpFac,nprim12*nAng,ielectron)
Call Mem_ODpointer_alloc(P%iprim1,nprim12,ielectron)
Call Mem_ODpointer_alloc(P%iprim2,nprim12,ielectron)
CALL MEM_ODPOINTER_ALLOC(P%distance12,3*npass,ielectron)
Call Mem_ODpointer_alloc(P%angmom,nAng,ielectron)
Call Mem_ODpointer_alloc(P%lenETUV,nAng,ielectron)
Call Mem_ODpointer_alloc(P%indexAng1,nAng,ielectron)
Call Mem_ODpointer_alloc(P%indexAng2,nAng,ielectron)
Call Mem_ODpointer_alloc(P%nOrbitals,nAng,ielectron)
Call Mem_ODpointer_alloc(P%nOrbComp,nAng,ielectron)
Call Mem_ODpointer_alloc(P%nContracted,nAng,ielectron)
Call Mem_ODpointer_alloc(P%totOrbitals,geoOrder+1,ielectron)

!nsize = mem_realsize*size(P%center)&
!        & + mem_realsize*size(P%exponents) + mem_realsize*size(P%reducedExponents) &
!        & + mem_realsize*size(P%preExpFac) + mem_intsize*size(P%iprim1)      + mem_intsize*size(P%iprim2)&
!        & + mem_intsize*size(P%lenETUV)    + mem_realsize*size(P%distance12) + mem_intsize*size(P%angmom)&
!        & + mem_intsize*size(P%indexAng1)  + mem_intsize*size(P%indexAng2) &
!        & + mem_intsize*size(P%nOrbitals)  + mem_intsize*size(P%nOrbComp)&
!        & + mem_intsize*size(P%nContracted)
!call mem_allocd_mem_overlap(nsize)

P%nAngAlloc = nAng
P%nTUVAlloc = nTUV
P%nPrimAlloc = nPrim12
P%FTUVisAlloc = .FALSE.
IF (allocFTUV) THEN
!   CALL MEM_ALLOC(P%FTUV,nPrim1,nTUV,ndmat)
   CALL MEM_ODpointer_ALLOC(P%FTUV,nPrim12*nTUV*ndmat,ielectron)
   P%FTUVisAlloc = .TRUE.
!   nsize = mem_realsize*size(P%FTUV)
!   call mem_allocated_mem_overlap(nsize)
ENDIF

END SUBROUTINE ALLOC_OVERLAP

!> \brief alloc the orbital
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Orb The orbital, to alloc
!> \param maxprim the maximum number of primitives
!> \param npass the number of passes 
!> \param maxangmom the maximum number of angular momentum
!> \param maxA the maximum number of spherical transformation matrices
!> \param FTUVorb flag to determine if this is a FTUV orbital
SUBROUTINE ALLOC_ORBITAL(Orb,maxprim,npass,maxangmom,maxA,FTUVorb,ielectron)
use memory_handling
IMPLICIT NONE
TYPE(Orbital) :: Orb
Integer :: maxprim,maxangmom,maxA,npass,ielectron
logical :: FTUVorb
integer (kind=long) :: nsize

CALL MEM_ODpointer_ALLOC(Orb%angmom,maxangmom,ielectron)
CALL MEM_ODpointer_ALLOC(Orb%nContracted,maxangmom,ielectron)
CALL MEM_ODpointer_ALLOC(Orb%nOrbComp,maxangmom,ielectron)
CALL MEM_ODpointer_ALLOC(Orb%startOrbital,maxangmom*npass,ielectron)
CALL MEM_ODpointer_ALLOC(Orb%startLocOrb,maxangmom,ielectron)

!nsize = mem_intsize*size(Orb%angmom)&
!        & + mem_intsize*size(Orb%nContracted) +  mem_intsize*size(Orb%nOrbComp)&
!        & + mem_intsize*size(Orb%startOrbital)&!+ mem_intsize*size(Orb%atom)& 
!        & + mem_intsize*size(Orb%startLocOrb) !+ mem_intsize*size(Orb%batch)
!call mem_allocated_mem_overlap(nsize)

IF(.NOT.FTUVorb)THEN
   CALL MEM_ODpointer_ALLOC(Orb%exponents,maxprim,ielectron)
   CALL MEM_ODpointer_ALLOC(Orb%nOrbitals,maxangmom,ielectron)
   
!   nsize = mem_realsize*size(Orb%exponents) + mem_intsize*size(Orb%startprimOrbital)&
!        & + mem_intsize*size(Orb%nOrbitals)
!   call mem_allocated_mem_overlap(nsize)
   
!   call mem_alloc(Orb%CC,maxangmom)
   Orb%FTUVorb = .FALSE.
ELSE
   Orb%FTUVorb = .TRUE.
!   NULLIFY(Orb%CC)
ENDIF

END SUBROUTINE ALLOC_ORBITAL

!> \brief set the orbitals when passes are used 
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param PassP The overlap, to be built
!> \param P the overlap
!> \param maxpasses the maximum number of passes
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SetPassOrbitals(PassP,P,maxPasses,LUPRI,IPRINT)
Implicit none
TYPE(Overlap)       :: PassP,P
Integer             :: maxPasses,LUPRI,IPRINT
!
Integer :: np,na,mc,ma,mp,ia
!
integer :: i,startPrim,endOrbital,ip

! Orbital 1
np = P%orbital1%nPrimitives
na = P%orbital1%nAngmom
mc = P%orbital1%maxContracted
ma = P%orbital1%maxAngmom

PassP%orbital1%TYPE_Empty = P%orbital1%type_Empty
PassP%orbital1%TYPE_Hermite = P%orbital1%type_Hermite
PassP%orbital1%TYPE_Cartesian = P%orbital1%type_Cartesian
PassP%orbital1%TYPE_Nucleus = P%orbital1%type_Nucleus
PassP%orbital1%TYPE_pCharge = P%orbital1%type_pCharge

PassP%orbital1%spherical     = P%orbital1%spherical     
PassP%orbital1%center        = P%orbital1%center        
PassP%orbital1%nPrimitives   = np
PassP%orbital1%maxContracted = mc
PassP%orbital1%maxAngmom     = ma
PassP%orbital1%nAngmom       = na
PassP%orbital1%CCidentifier  = P%orbital1%CCidentifier  
!replace with dcopy
PassP%orbital1%exponents(1:np) = P%orbital1%exponents(1:np)
PassP%orbital1%angmom(1:na)       = P%orbital1%angmom(1:na)
PassP%orbital1%nContracted(1:na)  = P%orbital1%nContracted(1:na)
PassP%orbital1%nOrbComp(1:na)     = P%orbital1%nOrbComp(1:na)
PassP%orbital1%nOrbitals(1:na)    = P%orbital1%nOrbitals(1:na)
PassP%orbital1%startLocOrb(1:na)  = P%orbital1%startLocOrb(1:na)
PassP%orbital1%totOrbitals        = P%orbital1%totOrbitals
!maybe replace with do loop
!   PassP%orbital1%CC(1:np,1:mc,1:na) = P%orbital1%CC(1:np,1:mc,1:na)
DO ia = 1,na
   PassP%orbital1%CC(ia)%p => P%orbital1%CC(ia)%p
ENDDO

! Orbital 2
np = P%orbital2%nPrimitives
na = P%orbital2%nAngmom
mc = P%orbital2%maxContracted
ma = P%orbital2%maxAngmom

PassP%orbital2%TYPE_Empty = P%orbital2%type_Empty
PassP%orbital2%TYPE_Hermite = P%orbital2%type_Hermite
PassP%orbital2%TYPE_Cartesian = P%orbital2%type_Cartesian
PassP%orbital2%TYPE_Nucleus = P%orbital2%type_Nucleus
PassP%orbital2%TYPE_pCharge = P%orbital2%type_pCharge

PassP%orbital2%spherical     = P%orbital2%spherical     
PassP%orbital2%center        = P%orbital2%center        
PassP%orbital2%nPrimitives   = np
PassP%orbital2%maxContracted = mc
PassP%orbital2%maxAngmom     = ma
PassP%orbital2%nAngmom       = na
PassP%orbital2%CCidentifier  = P%orbital2%CCidentifier  
PassP%orbital2%exponents(1:np) = P%orbital2%exponents(1:np)
PassP%orbital2%angmom(1:na)       = P%orbital2%angmom(1:na)
PassP%orbital2%nContracted(1:na)  = P%orbital2%nContracted(1:na)
PassP%orbital2%nOrbComp(1:na)     = P%orbital2%nOrbComp(1:na)
PassP%orbital2%nOrbitals(1:na)    = P%orbital2%nOrbitals(1:na)
PassP%orbital2%startLocOrb(1:na)  = P%orbital2%startLocOrb(1:na)
PassP%orbital2%totOrbitals        = P%orbital2%totOrbitals
DO ia = 1,na
   PassP%orbital2%CC(ia)%p => P%orbital2%CC(ia)%p
ENDDO

!Common Overlap
mp = P%nPrimitives*maxPasses
np = P%nPrimitives
na = P%nAngmom

PassP%type_Empty            = P%type_Empty
PassP%type_Hermite          = P%type_Hermite
PassP%type_Hermite_single   = P%type_Hermite_single
PassP%type_Cartesian        = P%type_Cartesian
PassP%type_Cartesian_single = P%type_Cartesian_single
PassP%type_Nucleus          = P%type_Nucleus
PassP%type_FTUV             = P%type_FTUV

PassP%single                = P%single
PassP%nGeoDerivComp         = P%nGeoDerivComp
PassP%nCartesianMomentComp  = P%nCartesianMomentComp

PassP%sphericalEcoeff   = P%sphericalEcoeff
PassP%sameAO            = P%sameAO
PassP%nPrimitives       = P%nPrimitives
PassP%maxContracted     = P%maxContracted
PassP%squaredDistance   = P%squaredDistance
!PassP%maxPrimGab        = P%maxPrimGab 
PassP%maxGab            = P%maxGab 
PassP%nAngmom           = P%nAngmom      
PassP%angmom(1:na)      = P%angmom(1:na)
PassP%indexAng1(1:na)   = P%indexAng1(1:na)
PassP%indexAng2(1:na)   = P%indexAng2(1:na)
PassP%nOrbComp(1:na)    = P%nOrbComp(1:na)
PassP%nContracted(1:na) = P%nContracted(1:na)
PassP%minAngmom         = P%minAngmom
PassP%maxAngmom         = P%maxAngmom
PassP%startAngmom       = P%startAngmom
PassP%startGeoOrder     = P%startGeoOrder
PassP%endGeoOrder       = P%endGeoOrder
PassP%endAngmom         = P%endAngmom
PassP%nTUV              = P%nTUV

PassP%mmorder           = P%mmorder
PassP%do_mulmom         = P%do_mulmom
startPrim = 0
endOrbital = nP-1
DO I = 1,maxpasses
   DO IP = 1,np
      PassP%exponents(startPrim+IP) = P%exponents(IP)
   ENDDO
   startPrim = startPrim+np
   endOrbital = endOrbital+np
ENDDO
!tmp
startPrim = 0
endOrbital = nP-1
DO I = 1,maxpasses
   DO IP = 1,np
      PassP%iPrim1(startPrim+IP) = P%iPrim1(IP)
      PassP%iPrim2(startPrim+IP) = P%iPrim2(IP)
   ENDDO
   startPrim = startPrim+np
   endOrbital = endOrbital+np
ENDDO
PassP%contractbasis  = P%contractbasis
PassP%segmented      = P%segmented
PassP%magderiv       = P%magderiv

END SUBROUTINE SetPassOrbitals

!> \brief Set the Final Passp overlap values
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param PassP The overlap, to be built
!> \param P the overlap
!> \param numPass the current number of passes in PassP 
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE FinalizePass(PassP,P,numPass,LUPRI,IPRINT)
Implicit none
TYPE(Overlap)       :: PassP,P
Integer             :: numPass,LUPRI,IPRINT
!
integer :: na,np

PassP%orbital1%nPasses  = numPass 
PassP%orbital2%nPasses  = numPass 
PassP%nPasses = numPass                                                             
PassP%totOrbitals(:) = P%totOrbitals(:)*numPass                                           
na = P%nAngmom
PassP%nOrbitals(1:na)   = P%nOrbitals(1:na)*numPass   

END SUBROUTINE FinalizePass

SUBROUTINE FREE_OVERLAP(P)
Implicit none
TYPE(Overlap) :: P
IF(.NOT. P%orbital1%FTUVorb)THEN
!   call mem_dealloc(P%orbital1%CC)
ENDIF
IF(.NOT. P%orbital2%FTUVorb)THEN
!   call mem_DEALLOC(P%orbital2%CC)
ENDIF
END SUBROUTINE FREE_OVERLAP

!> \brief set up the FTUV batch
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param P the overlap
!> \param INPUT the integral input specifications, contains all info about what to do
!> \param Integral contains arrays to store intermidiates and final integrals
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SET_FTUV(P,Input,Integral,LUPRI,IPRINT)
use memory_handling
Implicit none
TYPE(Overlap)       :: P
TYPE(IntegralInput) :: Input
TYPE(IntegralItem)  :: Integral
Integer             :: LUPRI,IPRINT
!
Integer             :: iAngmom,l1,l2,l,ijk,nTUV,lm1,lm2,lm
Integer             :: iPrim,iA1,iA2,nCont,nC1,nC2,iAng1,iAng2
Integer             :: idmat,start1,start2,i1,i2,iC1,iC2,indlm,iTUV
Integer             :: startAng, tuvOff,offset3,offset4
Real(realk),pointer :: Ecoeffs(:)
Real(realk),pointer :: CC(:),ContractedDmat(:)
Real(realk),parameter :: DM1=-1.0E0_realk
real(realk)          :: dtemp
Integer             :: SAA,atomC,atomD,batchC,batchD,Dcdindex,maxAng
Integer             :: Ddcindex,sC,sD,iSC,iSD,elms,n1,n2,offset,maxBat
Logical             :: spherical

#ifdef VAR_LSDEBUGINT
IF(size(P%FTUV).NE.P%nPrimAlloc*P%nTUVAlloc*Input%NDMAT_RHS)call lsquit('FTUV dim error',-1)
IF(P%nTUV*P%nPrimitives*Input%NDMAT_RHS.GT.size(P%FTUV))call lsquit('FTUV dim error',-1)
#endif
CALL LS_DZERO(P%FTUV,P%nTUV*P%nPrimitives*Input%NDMAT_RHS)

!dimC = P%orbital1%nOrbComp(1)*P%orbital1%nContracted(1)
!DO SAA=2,P%orbital1%nAngmom
!   dimC = dimC+P%orbital1%nOrbComp(SAA)*P%orbital1%nContracted(SAA)
!ENDDO
!dimD = P%orbital2%nOrbComp(1)*P%orbital2%nContracted(1)
!DO SAA=2,P%orbital2%nAngmom
!   dimD = dimD+P%orbital2%nOrbComp(SAA)*P%orbital2%nContracted(SAA)
!ENDDO

spherical = P%orbital1%spherical .OR. P%orbital2%spherical
IF (spherical .AND. .NOT.P%sphericalEcoeff) CALL LSQUIT('Error in SET_FTUV. Assumes sphericalEcoeff!',lupri)

AtomC = P%orb1atom(1)
AtomD = P%orb2atom(1)
BatchC = P%orb1batch(1)
BatchD = P%orb2batch(1)
Dcdindex = INPUT%LST_DRHS%INDEX(atomC,atomD,1,1)
maxBAt = INPUT%LST_DRHS%LSAO(Dcdindex)%maxBat
maxAng = INPUT%LST_DRHS%LSAO(Dcdindex)%maxAng
offset3 = (BatchC-1)*maxAng
offset4 = (BatchD-1)*maxAng+maxAng*maxBat
IF(Dcdindex.NE.0)THEN
   IF(Input%sameRHSaos) Ddcindex = INPUT%LST_DRHS%INDEX(atomD,atomC,1,1)
   n1 = INPUT%LST_DRHS%LSAO(Dcdindex)%nLocal(1)
   n2 = INPUT%LST_DRHS%LSAO(Dcdindex)%nLocal(2)

   DO iAngmom = 1,P%nAngmom
      iA1 = P%indexAng1(iangmom)
      iA2 = P%indexAng2(iangmom)
      l1 = P%orbital1%angmom(iA1)
      l2 = P%orbital2%angmom(iA2)
      l  = l1+l2
      CALL GET_IJK(l1,l2,lm1,lm2,lm,ijk,P%sphericalEcoeff,P%endGeoOrder,P%single)
      startAng = P%startAngmom
      tuvOff   = startAng*(startAng+1)*(startAng+2)/6
      nTUV = (l+1)*(l+2)*(l+3)/6 - tuvOff
!      print*,'debug1 Ecoeffs',P%nPrimitives*nTUV*lm
      call mem_workpointer_alloc(Ecoeffs,P%nPrimitives*nTUV*lm)
      CALL BuildEcoeffTensor(Integral%TUV,P,DM1,Ecoeffs,lm,ijk,nTUV,P%nPrimitives,&
           &                   iAngmom,0,1,LUPRI,IPRINT,1)
      ! Set up contraction matrix
      nC1 = P%orbital1%nContracted(iA1)
      nC2 = P%orbital2%nContracted(iA2)
      nCont = nC1*nC2
!      print*,'debug1 CC',nC1*nC2*P%nPrimitives
      call mem_workpointer_alloc(CC,nC1*nC2*P%nPrimitives)
      CALL ConstructContraction_PA(CC,P,P%nPrimitives,nC1,nC2,iA1,iA2,LUPRI,IPRINT)
      
      ! Make contraction of density-matrix and contraction coefficients
!      print*,'debug1 ContractedDmat',P%nPrimitives*lm*Input%NDMAT_RHS
      call mem_workpointer_alloc(ContractedDmat,P%nPrimitives*lm*Input%NDMAT_RHS)
      CALL LS_DZERO(ContractedDmat,lm*P%nPrimitives*Input%NDMAT_RHS)
      start1 = P%orbital1%startOrbital(iA1)
      start2 = P%orbital2%startOrbital(iA2)
      
      sC = INPUT%LST_DRHS%LSAO(Dcdindex)%startLocalOrb(iA1+offset3)-1
      sD = INPUT%LST_DRHS%LSAO(Dcdindex)%startLocalOrb(iA2+offset4)-1
      IF (.NOT.input%contAng) THEN !Default AO ordering
        call ftuvbasiscontraction(n1,n2,nC2,lm2,sD,nC1,lm1,sC,Input%nDMAT_RHS,INPUT%LST_DRHS,&
             & Input%sameRHSaos,start1,start2,Dcdindex,Ddcindex,Input%CoulombFactor,P%nPrimitives,&
             & lm,CC,ContractedDmat)
      ELSE !Optional AO ordering: contraced first angular second
        call ftuvbasiscontraction_ca(n1,n2,nC2,lm2,sD,nC1,lm1,sC,Input%nDMAT_RHS,INPUT%LST_DRHS,&
             & Input%sameRHSaos,start1,start2,Dcdindex,Ddcindex,Input%CoulombFactor,P%nPrimitives,&
             & lm,CC,ContractedDmat)
      ENDIF

      call mem_workpointer_dealloc(CC)
#ifdef VAR_LSDEBUGINT
      IF(P%nPrimitives*nTUV*Input%NDMAT_RHS.GT.P%nPrimAlloc*P%nTUVAlloc*Input%NDMAT_RHS)call lsquit('FTUV dim error',-1)
#endif
      CALL SetUpFTUVs1(P%FTUV,P%nPrimitives,nTUV,Input%NDMAT_RHS,lm,ContractedDmat,Ecoeffs,P%nTUVAlloc,P%nPrimAlloc)
      call mem_workpointer_dealloc(ContractedDmat)
      call mem_workpointer_dealloc(Ecoeffs)
   ENDDO
!ELSE
   !All Density matrix elements are zero so the FTUV is zero. 
ENDIF

IF (IPRINT.GT. 20) THEN
   CALL LSHEADER(LUPRI,'SET_FTUV')
   WRITE(LUPRI,'(1X,A,I5)') 'Number of primitives', P%nPrimitives
   WRITE(LUPRI,'(1X,A,I5)') "Number of TUV's     ", nTUV
   WRITE(LUPRI,'(1X,A,I5)') "Number of D-mat     ", Input%NDMAT_RHS
   DO idmat=1,Input%NDMAT_RHS
      DO iTUV=1,nTUV
         offset = (iTUV-1)*P%nPrimAlloc+(idmat-1)*P%nPrimAlloc*P%nTUVAlloc
         WRITE(LUPRI,'(1X,A,I5,A,I5)') 'iTUV =', iTUV, ' iDmat=',idmat
         WRITE(LUPRI,'(3X,5ES16.8/,(3X,5ES16.8))') (P%FTUV(iPrim+offset),iPrim=1,P%nPrimitives)
      ENDDO
   ENDDO
ENDIF
END SUBROUTINE SET_FTUV

SUBROUTINE SetUpFTUVs1(FTUV,nPrim,nTUV,NDMAT_RHS,lm,ContractedDmat,Ecoeffs,nTUVAlloc,nPrimAlloc)
  implicit none
  integer,intent(in) :: nPrim,nTUV,NDMAT_RHS,lm,ntuvAlloc,nPrimAlloc
  REAL(REALK),intent(inout) :: FTUV(nPrimAlloc*nTUVAlloc*NDMAT_RHS)
  REAL(REALK),intent(in) :: ContractedDmat(nPrim,lm,ndmat_rhs)
  REAL(REALK),intent(in) :: Ecoeffs(nPrim,nTUV,lm)
  !
  integer :: idmat,indlm,iTUV,iPrim,offset
  ! Set up FTUV's
  DO idmat=1,NDMAT_RHS
   DO indlm=1,lm
    DO iTUV=1,nTUV
     offset = (iTUV-1)*nPrimAlloc+(idmat-1)*nPrimAlloc*nTUVAlloc
     DO iPrim=1,nPrim
      FTUV(iPrim+offset) = FTUV(iPrim+offset) &
           & + ContractedDmat(iPrim,indlm,idmat)*Ecoeffs(iPrim,iTUV,indlm)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
END SUBROUTINE SETUPFTUVS1

subroutine ftuvbasiscontraction(n1,n2,nC2,lm2,sD,nC1,lm1,sC,ndmat,DRHS,&
     & sameRHSaos,start1,start2,Dcdindex,Ddcindex,CoulombFactor,nPrim,&
     & lmAlloc,CC,ContractedDmat)
  implicit none
  integer,intent(in) :: n1,n2,nC2,lm2,sD,nC1,lm1,sC,ndmat,start1,start2
  integer,intent(in) :: Dcdindex,Ddcindex,nPrim,lmAlloc
  type(lstensor),intent(in) :: DRHS
  logical,intent(in) :: sameRHSaos
  real(realk),intent(inout) :: ContractedDmat(nPrim,lmAlloc,ndmat)
  real(realk),intent(in) :: CC(nC1,nC2,nPrim),CoulombFactor
  !
  integer :: i2,iC2,iAng2,iSD,i1,iC1,iAng1,indlm,iSC,elms,iPrim,idmat
  real(realk) :: dtemp
  
  DO idmat=1,nDMAT
     i2 = 0
     DO iC2=1,nC2
        DO iAng2=1,lm2
           i2 = i2 + 1
           iSD = i2+sD
           i1=0
           DO iC1=1,nC1
              DO iAng1=1,lm1
                 indlm = (iAng2-1)*lm1 + iAng1
                 i1 = i1+1
                 iSC = i1+sC
                 elms = iSC + (iSD-1)*n1 + (idmat-1)*n1*n2
                 dtemp = DRHS%LSAO(Dcdindex)%elms(elms)
                 IF ((start1.NE.start2).AND.sameRHSaos) THEN
                    elms = iSD + (iSC-1)*n2 + (idmat-1)*n1*n2
                    dtemp = dtemp + DRHS%LSAO(Ddcindex)%elms(elms)
                 ENDIF
                 dtemp = dtemp*CoulombFactor
                 DO iPrim=1,nPrim
                    ContractedDmat(iPrim,indlm,idmat) = ContractedDmat(iPrim,indlm,idmat) &
                         &                                    + dtemp * CC(iC1,iC2,iPrim)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
end subroutine ftuvbasiscontraction

!This subroutine does the same as ftuvbasiscontraction, except that it assumes
!AO ordering to be contracted first (inner) angular second
subroutine ftuvbasiscontraction_ca(n1,n2,nC2,lm2,sD,nC1,lm1,sC,ndmat,DRHS,&
     & sameRHSaos,start1,start2,Dcdindex,Ddcindex,CoulombFactor,nPrim,&
     & lmAlloc,CC,ContractedDmat)
  implicit none
  integer,intent(in) :: n1,n2,nC2,lm2,sD,nC1,lm1,sC,ndmat,start1,start2
  integer,intent(in) :: Dcdindex,Ddcindex,nPrim,lmAlloc
  type(lstensor),intent(in) :: DRHS
  logical,intent(in) :: sameRHSaos
  real(realk),intent(inout) :: ContractedDmat(nPrim,lmAlloc,ndmat)
  real(realk),intent(in) :: CC(nC1,nC2,nPrim),CoulombFactor
  !
  integer :: i2,iC2,iAng2,iSD,i1,iC1,iAng1,indlm,iSC,elms,iPrim,idmat
  real(realk) :: dtemp
  
  DO idmat=1,nDMAT
     i2 = 0
     DO iAng2=1,lm2
        DO iC2=1,nC2
           i2 = i2 + 1
           iSD = i2+sD
           i1=0
           DO iAng1=1,lm1
              DO iC1=1,nC1
                 indlm = (iAng2-1)*lm1 + iAng1
                 i1 = i1+1
                 iSC = i1+sC
                 elms = iSC + (iSD-1)*n1 + (idmat-1)*n1*n2
                 dtemp = DRHS%LSAO(Dcdindex)%elms(elms)
                 IF ((start1.NE.start2).AND.sameRHSaos) THEN
                    elms = iSD + (iSC-1)*n2 + (idmat-1)*n1*n2
                    dtemp = dtemp + DRHS%LSAO(Ddcindex)%elms(elms)
                 ENDIF
                 dtemp = dtemp*CoulombFactor
                 DO iPrim=1,nPrim
                    ContractedDmat(iPrim,indlm,idmat) = ContractedDmat(iPrim,indlm,idmat) &
                         &                                    + dtemp * CC(iC1,iC2,iPrim)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
end subroutine ftuvbasiscontraction_ca

!> \brief Finds the number of significant product primitives
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param nPrimitives the number of nPrimitives for this overlap 
!> \param P the overlap
!> \param INPUT the integral input specifications, contains all info about what to do
!> \param side LHS or RHS
!> \param LUPRI the logical unit number for the output file
SUBROUTINE GET_NPRIMITIVES(nPrimitives,P,Input,ielectron,lupri)
Implicit none
TYPE(IntegralInput),target :: Input
TYPE(Overlap)       :: P
Integer             :: nPrimitives,lupri,ielectron
!Character*(*)       :: Side
!
Integer              :: i1,i2,start2
integer(kind=short)  :: maxgab,maxelm
type(lstensor), pointer :: GAB
integer :: atomA,atomB,batchA,batchB,Gindex,n1,s1,s2,maxbat
!real(realk),parameter :: TEN=1E1_realk
IF (ielectron.EQ. 1) THEN
  GAB => INPUT%LST_GAB_LHS
  MAXELM=INPUT%PS_MAXELM_RHS
ELSE
  GAB => INPUT%LST_GAB_RHS
  MAXELM=INPUT%PS_MAXELM_LHS
ENDIF
nPrimitives = 0
IF(INPUT%PS_SCREEN.AND.(.NOT.P%type_Empty))THEN
  AtomA = P%orb1atom(1)
  AtomB = P%orb2atom(1)
  BatchA = P%orb1batch(1)
  BatchB = P%orb2batch(1)
  Gindex = GAB%INDEX(atomA,atomB,1,1)
#ifdef VAR_LSDEBUGINT
  IF(Gindex.EQ.0)THEN
     call lsquit('lstensor allocation error in SET_OVERLAP',lupri)
  ENDIF
#endif
  maxbat = GAB%SLSAO(Gindex)%maxbat
  s1 = GAB%SLSAO(Gindex)%startLocalOrb(BatchA)-1
  s2 = GAB%SLSAO(Gindex)%startLocalOrb(BatchB+maxbat)-1
  n1 = GAB%SLSAO(Gindex)%nLocal(1)
  DO i2=1,P%orbital2%nPrimitives
     DO i1=1,P%orbital1%nPrimitives
        maxGab = GAB%SLSAO(Gindex)%selms(s1+i1 + (s2+i2-1)*n1)
        IF( MAXGAB .GT. INPUT%PS_THRLOG-MAXELM)THEN
            nPrimitives = nPrimitives + 1
        ENDIF
      ENDDO
   ENDDO
ELSE
   i1 = P%orbital1%nPrimitives
   i2 = P%orbital2%nPrimitives
   nPrimitives = i1*i2
ENDIF

END SUBROUTINE GET_NPRIMITIVES

!> \brief set the values of the orbital
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Orb the orbital to be set
!> \param AOB the Atomic batch
!> \param INTEGRAL storage of integrals and intermediates
!> \param LUPRI the logical unit number for the output file
SUBROUTINE SET_ORBITAL(Orb,AOB,integral,hermiteEcoeff,LUPRI)
IMPLICIT NONE
TYPE(Orbital) :: Orb
TYPE(AOBATCH) :: AOB
TYPE(IntegralItem)  :: Integral
integer             :: lupri
logical             :: hermiteEcoeff
!
Integer :: ia,na,np,i
!
 np   = AOB%nPrimitives
 na   = AOB%nAngmom

 Orb%TYPE_Empty = AOB%type_Empty
 IF(Orb%TYPE_Empty)then
    Orb%TYPE_Hermite = .FALSE.
    Orb%TYPE_Cartesian = .FALSE.
 ELSE
    IF(HermiteEcoeff)then
       Orb%TYPE_Hermite = .TRUE.
       Orb%TYPE_Cartesian = .FALSE.
    ELSE
       Orb%TYPE_Hermite = .FALSE.
       Orb%TYPE_Cartesian = .TRUE.
    ENDIF
 ENDIF
 Orb%TYPE_Nucleus = AOB%type_Nucleus
 Orb%TYPE_pCharge = AOB%type_pCharge

 Orb%spherical     = AOB%spherical
 Orb%center        = AOB%center
 Orb%nPrimitives   = np 
 Orb%nPasses       = 1 
 Orb%maxContracted = AOB%maxContracted
 Orb%maxAngmom     = AOB%maxAngmom
 Orb%nAngmom       = na
 Orb%CCidentifier  = AOB%CCindex(1)
 Orb%exponents(1:np) = AOB%pExponents%elms(1:np)
 Orb%angmom(1:na) = AOB%angmom(1:na)
 Orb%nContracted(1:na)  = AOB%nContracted(1:na)
 Orb%startOrbital(1:na) = AOB%startOrbital(1:na)
! Orb%atom(1) = AOB%atom
! Orb%batch(1) = AOB%batch
 Orb%nOrbComp(1:na) = AOB%nOrbComp(1:na)
 Orb%nOrbitals(1:na) = AOB%nOrbitals(1:na)
 Orb%totOrbitals   = 0 
 DO ia=1,na
    Orb%startLocOrb(ia) = Orb%totOrbitals + 1
    Orb%totOrbitals = Orb%totOrbitals + Orb%nOrbitals(ia)
 ENDDO
 DO ia=1,na
    Orb%CC(ia)%p => AOB%pCC(ia)%p
 ENDDO
END SUBROUTINE SET_ORBITAL

!> \brief print the overlap
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param P the overlap distribution
!> \param IUNIT the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
!> \param side LHS or RHS
SUBROUTINE PRINT_OVERLAP(P,IUNIT,IPRINT,SIDE)
IMPLICIT NONE
TYPE(Overlap) :: P
Integer       :: IUNIT,IPRINT
Character*(*)      :: Side
!
integer :: i,ipass,ip,ia
CALL LSHEADER(IUNIT,'OVERLAP')
IF (SIDE.EQ.'RHS') THEN
   WRITE(IUNIT,'(1X,A)') ''
   WRITE(IUNIT,'(1X,A)') '*************************'
   WRITE(IUNIT,'(1X,A)') '***   OVERLAP RHS     ***'
   WRITE(IUNIT,'(1X,A)') '*************************'
   WRITE(IUNIT,'(1X,A)') ''
ELSEIF (SIDE.EQ.'LHS') THEN
   WRITE(IUNIT,'(1X,A)') ''
   WRITE(IUNIT,'(1X,A)') '*************************'
   WRITE(IUNIT,'(1X,A)') '***   OVERLAP LHS     ***'
   WRITE(IUNIT,'(1X,A)') '*************************'
   WRITE(IUNIT,'(1X,A)') ''
ELSE
   WRITE(IUNIT,'(1X,A)') ''
   WRITE(IUNIT,'(1X,A)') '*************************'
   WRITE(IUNIT,'(1X,A)') '***   OVERLAP         ***'
   WRITE(IUNIT,'(1X,A)') '*************************'
   WRITE(IUNIT,'(1X,A)') ''
ENDIF
IF(.NOT.P%type_FTUV)THEN
   WRITE(IUNIT,'(3X,A)') '-------------- orbital 1 --------------'
   CALL PRINT_ORBITAL(P%orbital1,P%nPasses,IUNIT,IPRINT)
   WRITE(IUNIT,'(3X,A)') '-------------- orbital 2 --------------'
   CALL PRINT_ORBITAL(P%orbital2,P%nPasses,IUNIT,IPRINT)
ENDIF
 WRITE(IUNIT,'(3X,A)') '--------- overlap specifics -----------'
 WRITE(IUNIT,'(5X,A,L1)')  'type empty            =', P%type_Empty
 WRITE(IUNIT,'(5X,A,L1)')  'type cartesian-signle =', P%type_Cartesian_single
 WRITE(IUNIT,'(5X,A,L1)')  'type cartesian        =', P%type_Cartesian
 WRITE(IUNIT,'(5X,A,L1)')  'type hermite          =', P%type_hermite
 WRITE(IUNIT,'(5X,A,L1)')  'type hermite_single   =', P%type_hermite_single
 WRITE(IUNIT,'(5X,A,L1)')  'type ftuv             =', P%type_FTUV
 WRITE(IUNIT,'(5X,A,L1)')  'type Nucleus          =', P%type_Nucleus
 WRITE(IUNIT,'(5X,A,L1)')  'Single                =', P%single
 WRITE(IUNIT,'(5X,A,L1)')  'Spherical-Ecoeff      =', P%sphericalEcoeff

 WRITE(IUNIT,'(5X,A,I3)') 'nAngmom          =', P%nAngmom
 WRITE(IUNIT,'(5X,A,I3)') 'nPrimitives      =', P%nPrimitives
 WRITE(IUNIT,'(5X,A,I3)') 'nPasses          =', P%nPasses
 WRITE(IUNIT,'(5X,A,10I3)') 'totOrbitals(iDer)=', (P%totOrbitals(i+1),i=P%startGeoORder,P%endGeoOrder)
 WRITE(IUNIT,'(5X,A,I3)') 'nTUV             =', P%nTUV
 WRITE(IUNIT,'(5X,A,I3)') 'maxContracted    =', P%maxContracted
 WRITE(IUNIT,'(5X,A,I3)') 'minAngmom        =', P%minAngmom
 WRITE(IUNIT,'(5X,A,I3)') 'maxAngmom        =', P%maxAngmom
 WRITE(IUNIT,'(5X,A,I3)') 'startAngmom      =', P%startAngmom
 WRITE(IUNIT,'(5X,A,I3)') 'endAngmom        =', P%endAngmom
 WRITE(IUNIT,'(5X,A,I3)') 'startGeoOrder    =', P%startGeoOrder
 WRITE(IUNIT,'(5X,A,I3)') 'endGeoOrder      =', P%endGeoOrder
 WRITE(IUNIT,'(5X,A,I3)') 'nGeoDerivComp    =', P%nGeoDerivComp
 WRITE(IUNIT,'(5X,A,I3)') 'nCartesianMomentComp=', P%nCartesianMomentComp
 WRITE(IUNIT,'(5X,A,L1)') 'do_mulmom        =', P%do_mulmom
 WRITE(IUNIT,'(5X,A,I3)') 'magderiv         =', P%magderiv
 WRITE(IUNIT,'(5X,A,I3)') 'mmorder          =', P%mmorder
IF(.NOT. P%type_FTUV)THEN
   DO ipass=1,P%nPasses
      WRITE(IUNIT,'(5X,A,3F8.4,A,I3)')'distance12(A) =', (P%distance12(i+(ipass-1)*3), i=1,3),' for pass ',ipass
   ENDDO
      WRITE(IUNIT,'(5X,A,F8.4)') 'squaredDistance  =', P%squaredDistance
ENDIF
 WRITE(IUNIT,'(5X,A,8I3)')'angmom           =', (P%angmom(i),i=1,P%nAngmom)
 WRITE(IUNIT,'(5X,A,8I3)')'indexAng1        =', (P%indexAng1(i),i=1,P%nAngmom)
 WRITE(IUNIT,'(5X,A,8I3)')'indexAng2        =', (P%indexAng2(i),i=1,P%nAngmom)
 WRITE(IUNIT,'(5X,A,8I3)')'nOrbitals        =', (P%nOrbitals(i),i=1,P%nAngmom)
 WRITE(IUNIT,'(5X,A,8I3)')'nOrbComp         =', (P%nOrbComp(i),i=1,P%nAngmom)
 WRITE(IUNIT,'(5X,A,8I3)')'nContracted      =', (P%nContracted(i),i=1,P%nAngmom)
 WRITE(IUNIT,'(5X,A,1L2)')  'sameAO           =  ', P%sameAO
 WRITE(IUNIT,'(3X,A)')    '------------------------------------------------------------------------'
IF(.NOT.P%type_FTUV)THEN
   DO ipass=1,P%nPasses
      WRITE(IUNIT,'(3X,A,I5)') '*** Pass number ',ipass
      WRITE(IUNIT,'(3X,A)')    '-------------------------------------------------------------------------'
      WRITE(IUNIT,'(3X,A)')    '  iPrim   Px      Py      Pz        p         mu    iprim1 iprim2  '
      ip = (ipass-1)*P%nPrimitives
      DO i=1+ip,P%nPrimitives+ip
         WRITE(IUNIT,'(5X,I4,3F8.4,2ES10.3,2I7)') i, P%center(1+(i-1)*3), P%center(2+(i-1)*3), P%center(3+(i-1)*3), &
              &  P%exponents(i), P%reducedExponents(i), P%iprim1(i), P%iprim2(i)
      ENDDO
      WRITE(IUNIT,'(3X,A)')    '-------------------------------------------------------------------------'
      WRITE(IUNIT,'(3X,A)')    ' iPrim   iAngmom   PreExpFac '
#ifdef VAR_LSDEBUGINT
      IF(size(P%preExpFac).NE.P%nPrimAlloc*P%nAngAlloc)call lsquit('PreExpFac dim error3',-1)
#endif
      do ia=1,P%nAngmom
       DO i=1+ip,P%nPrimitives+ip
          WRITE(IUNIT,'(5X,2I6,1ES16.8)') i, ia, P%preExpFac(i+(ia-1)*P%nPrimitives*P%nPasses)
       ENDDO
      enddo
   ENDDO
ELSE
   DO ipass=1,P%nPasses
      WRITE(IUNIT,'(3X,A,I5)') '*** Pass number ',ipass
      WRITE(IUNIT,'(3X,A)')    '-------------------------------------------------------------------------'
      WRITE(IUNIT,'(3X,A)')    '  iPrim   Px      Py      Pz        p         mu    '
      ip = (ipass-1)*P%nPrimitives
      DO i=1+ip,P%nPrimitives+ip
         WRITE(IUNIT,'(5X,I4,3F8.4,2ES10.3)') i, P%center(1+(i-1)*3), P%center(2+(i-1)*3), P%center(3+(i-1)*3), &
              &  P%exponents(i), P%reducedExponents(i)
      ENDDO
#ifdef VAR_LSDEBUGINT
      IF(size(P%preExpFac).NE.P%nPrimAlloc*P%nAngAlloc)call lsquit('PreExpFac dim error4',-1)
#endif
      WRITE(IUNIT,'(3X,A)')    '-------------------------------------------------------------------------'
      WRITE(IUNIT,'(3X,A)')    ' iPrim   iAngmom   PreExpFac '
      do ia=1,P%nAngmom
       DO i=1+ip,P%nPrimitives+ip
          WRITE(IUNIT,'(5X,2I6,1ES16.8)') i, ia, P%preExpFac(i+(ia-1)*P%nPrimitives*P%nPasses)
       ENDDO
      enddo
   ENDDO
ENDIF
WRITE(IUNIT,'(3X,A)')    '------------------------------------------------------------------------'
END SUBROUTINE PRINT_OVERLAP

!> \brief print the orbital
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Orb the orbital to be printet
!> \param nPasses the number of passes
!> \param IUNIT the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE PRINT_ORBITAL(Orb,nPasses,IUNIT,IPRINT)
TYPE(Orbital),intent(IN) :: Orb
Integer,intent(IN)       :: IUNIT,IPRINT,nPasses
!
Integer :: I,J,K,firstAtom,iCent
Logical :: printCenter

 WRITE(IUNIT,'(5X,A17,1L2)')'Type Empty     = ', Orb%type_Empty
 WRITE(IUNIT,'(5X,A17,1L2)')'Type Cartesian = ', Orb%type_Cartesian
 WRITE(IUNIT,'(5X,A17,1L2)')'Type hermite   = ', Orb%type_hermite
 WRITE(IUNIT,'(5X,A17,1L2)')'Type nucleus   = ', Orb%type_Nucleus
 WRITE(IUNIT,'(5X,A17,1L2)')'Type pCharge   = ', Orb%type_pCharge

 WRITE(IUNIT,'(5X,A17,1L2)')'Spherical      = ', Orb%spherical
 WRITE(IUNIT,'(5X,A17,I3)') 'maxAngmom      = ', Orb%maxAngmom
 WRITE(IUNIT,'(5X,A17,I3)') 'nAngmom        = ', Orb%nAngmom
 WRITE(IUNIT,'(5X,A17,I3)') 'nPrimitives    = ', Orb%nPrimitives
 WRITE(IUNIT,'(5X,A17,I3)') 'totOrbitals    = ', Orb%totOrbitals
 WRITE(IUNIT,'(5X,A17,I3)') 'maxContracted  = ', Orb%maxContracted
!Print center and atomic information (special case for passes)
 IF (.NOT.Orb%type_Empty) THEN
   IF (nPasses.EQ. 1) THEN
     WRITE(IUNIT,'(5X,A17,3F12.6)') 'center (A)     = ', (Orb%center(i),i=1,3)
   ELSE
     printCenter = .TRUE.
!     firstAtom = Orb%Atom(1)
     DO iCent=2,nPasses
!       IF (firstAtom.NE.Orb%Atom(iCent)) THEN
         printCenter = .FALSE.
!         EXIT
!       ENDIF
     ENDDO
     IF (printCenter) THEN
       WRITE(IUNIT,'(5X,A17,3F12.6)') 'Pass-center (A) =', (Orb%center(i),i=1,3)
     ELSE
       WRITE(IUNIT,'(5X,A)')        'Pass contains orbitals that do not share a common center'
     ENDIF
   ENDIF
!   WRITE(IUNIT,'(5X,A15,10I5,/(20X,10I5))') 'Atomic index = ',(Orb%Atom(iCent),iCent=1,nPasses)
 ENDIF
 WRITE(IUNIT,'(3X,A)')    '-------------------------------------------------------------------'
 WRITE(IUNIT,'(3X,A)')    '-  Block   angmom  #cont.  1.orb.  1.local #orb.   #oc   #pri.oc  1.porb  -'
 DO I=1,Orb%nAngmom
   WRITE(IUNIT,'(1X,7I8)')  I, Orb%angmom(I), Orb%nContracted(I), Orb%startOrbital(I), &
     & Orb%startLocOrb(I), Orb%nOrbitals(I), Orb%nOrbComp(I)
 ENDDO
 WRITE(IUNIT,'(3X,A)')    '-------------------------------------------------------------------'
 WRITE(IUNIT,'(3X,A)')    '----------------- Exponents -----------------'
 WRITE(IUNIT,'(5X,5F12.6)') (Orb%exponents(I),I=1,Orb%nPrimitives)
 WRITE(IUNIT,'(3X,A)')    '----------------- Contraction Coefficients -----------------'
 DO K=1,Orb%nAngmom
   WRITE(IUNIT,'(5X,A,I3,A,I3)') '* Angular block number',K,' with angular momentum', Orb%angmom(K)
   DO J=1,Orb%nContracted(K)
!     WRITE(IUNIT,'(5X,I3,5ES12.6,/(8X,5ES12.6))') J,(Orb%CC(I,J,K),I=1,Orb%nPrimitives)
     WRITE(IUNIT,'(5X,I3,5F12.6,/(8X,5F12.6))') J,(Orb%CC(K)%p%elms(I+(J-1)*Orb%nPrimitives),I=1,Orb%nPrimitives)
   ENDDO
!   IF(Orb%angmom(K) .GE. 2)THEN
!      WRITE(IUNIT,'(7X,A,I3)')'Spherical transformation matrix for angular momentum', Orb%angmom(K)
!      CALL OUTPUT(orb%SPH_MAT(Orb%angmom(K)+1)%p%elms,1,2*Orb%angmom(K)+1,&
!         &1,(Orb%angmom(K)+1)*(Orb%angmom(K)+2)/2,2*Orb%angmom(K)+1,&
!         &(Orb%angmom(K)+1)*(Orb%angmom(K)+2)/2,1,IUNIT)
!   ENDIF
      WRITE(IUNIT,*)' '
 ENDDO

END SUBROUTINE PRINT_ORBITAL

!> \brief wrapper routine that brach out and build Ecoefficient tensor
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param TUV The TUV indeces
!> \param P the overlap distribution
!> \param signP the sign (1 for P and -1 for Q)
!> \param Ecoeffs the ecoefficients to be built
!> \param lm the number of spherical angular components
!> \param Nijk the number of cartesian angular components
!> \param ntuvP the number of hermite angular components
!> \param nP the number of primitives
!> \param iAngmom the angular moment index
!> \param nPasses the number of passes
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE BuildEcoeffTensor(TUV,P,signP,Ecoeffs,lm,Nijk,tuvP,nP,iAngmom,der,nPasses,LUPRI,IPRINT,nOperatorComp)
use memory_handling
implicit none
TYPE(TUVitem)      :: TUV
TYPE(Overlap)      :: P
Integer,intent(IN) :: iAngmom,der,nPasses,LUPRI,IPRINT
Integer            :: lm,Nijk,tuvP,nP
Integer            :: nOperatorComp
Real(realk)        :: signP
Real(realk),target :: Ecoeffs(nP*tuvP*lm*nOperatorComp) 
!
Real(realk),pointer:: EcoeffN(:)
Real(realk),pointer :: SpherMat(:)
Integer                 :: ijk1,ijk2,lm1,lm2,ang1,ang2,iA1,iA2,X
Logical                 :: Sph1,Sph2,Spherical

iA1  = P%indexAng1(iangmom)
iA2  = P%indexAng2(iangmom)
ang1 = P%orbital1%angmom(iA1)
ang2 = P%orbital2%angmom(iA2)
Sph1 = P%orbital1%spherical.AND.(ang1.GT. 1)
Sph2 = P%orbital2%spherical.AND.(ang2.GT. 1)
Spherical = P%sphericalEcoeff.AND.(Sph1.OR.Sph2)
!Test consistency for derivatice case
#ifdef VAR_LSDEBUGINT
IF (Spherical.AND.P%endGeoOrder.GT. 0) CALL LSQUIT('Error in BuildEcoeffTensor. Spherical and deriv>0!',lupri)
#endif
IF (Spherical) THEN
  NULLIFY(EcoeffN)
  CALL mem_workpointer_alloc(EcoeffN,nP*tuvP*Nijk*nOperatorComp)
  ijk1 = (ang1+1)*(ang1+2)/2
  ijk2 = (ang2+1)*(ang2+2)/2
  lm1 = 2*ang1+1
  lm2 = 2*ang2+1
  CALL BuildEcoeffTensor_PA(TUV,P,signP,EcoeffN,Nijk*nOperatorComp,tuvP,nP,&
       & iAngmom,der,nPasses,LUPRI,IPRINT)
  CALL mem_workpointer_alloc(SpherMat,ijk1*ijk2*lm1*lm2)
  CALL SphericalTransformation(TUV,SpherMat,ijk1,ijk2,lm1,lm2,ang1,ang2,&
          & EcoeffN,Ecoeffs,tuvP*nP,LUPRI,IPRINT,nOperatorComp)
  CALL mem_workpointer_dealloc(SpherMat)
  CALL mem_workpointer_dealloc(EcoeffN)
ELSE
  CALL BuildEcoeffTensor_PA(TUV,P,signP,Ecoeffs,Nijk*nOperatorComp,&
       & tuvP,nP,iAngmom,der,nPasses,LUPRI,IPRINT)
ENDIF

IF (IPRINT.GT. 50) THEN
  CALL PrintTensor(Ecoeffs,'E-coefficienttensor ',&
     &             nP,tuvP,Nijk*nOperatorComp,Lupri,&
     &             'prim  ','tuv   ','ijk   ',3)
ENDIF

END SUBROUTINE BuildEcoeffTensor

!> \brief build Ecoefficient tensor (ordering ijk,nAng,nprim)
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param TUV The TUV indeces
!> \param P the overlap distribution P or Q
!> \param signP the sign (1 for P and -1 for Q)
!> \param Ecoeffs the ecoefficients to be built
!> \param Nijk the number of cartesian angular components
!> \param ntuv the number of hermite angular components
!> \param nPrim the number of primitives
!> \param iAngmom the angular moment index
!> \param nPasses the number of passes
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE BuildEcoeffTensor_PA(TUV,P,signP,Ecoeffs,Nijk,ntuv,nprim,iAngmom,der,nPasses,LUPRI,IPRINT)
!different ordering of Ecoeffs(ijk,nAng,nprim)
implicit none
TYPE(TUVitem)      :: TUV
TYPE(Overlap)      :: P
Integer            :: iAngmom,der,LUPRI,IPRINT,nPasses
Integer            :: Nijk,ntuv,nprim
Real(realk)        :: signP
!Allocatable for ECOEFFS - Dimensions are nPrimitives, Zero to Sum of Max Angmom 1 and 2, 
!                          Zero to Max Angmom 1, Zero to Max Angmom 2, 3 (for X,Y,Z)
real(realk),pointer  :: ETIJ(:,:,:,:,:) 
integer             :: P1,P2,tP,uP,vP
integer             :: ituvP,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
integer             :: l1,l2,magder,nX,nmstart,nmend,Xdir,x1,x2,x3,nt,nu,nv
Real(realk)         :: Ecoeffs(nprim,ntuv,Nijk) !ijk1*ijk2,nAng,nprim
Real(realk)         :: signijk
Integer             :: increment(2),icent,ncent,iPrim1,iPrim2,ijk,startAng,ioff,offset
Real(realk),pointer :: pref(:)
Real(realk),parameter :: D1=1E0_realk

l1       = P%orbital1%angmom(P%indexAng1(iangmom))
l2       = P%orbital2%angmom(P%indexAng2(iangmom))
magder   = P%magderiv
startAng = P%startAngmom
ioff     = startAng*(startAng+1)*(startAng+2)/6
nX=3
IF(magder.EQ.1)nX=6
 
IF (der.GT. 0) THEN
  IF (.NOT.P%type_Hermite) CALL LSQUIT('Error in BuildEcoeffTensor_PA. endGeoOrder>0 and not type hermite',lupri)
ENDIF

!For derivative case the ETIJ can be split into different components if time-critical, i.e.
!ETIJ(l1+l2+1,l1+1,l2) and ETIJ(l1+l2+1,l1,l2+1) for gradients etc.
call mem_alloc(ETIJ,nprim,l1+l2+2*der+magder,l1+der+magder,l2+der,nX,.FALSE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.)

CALL GET_ECOEFF(ETIJ,nprim,l1+l2+2*der+magder,l1+der+magder,l2+der,P%nPrimitives,nPasses,P,LUPRI,IPRINT)
IF(magder.GT.0)THEN
   IF(P%magderiv.GT.1)call lsquit('only first order magnetic derivative have been implemented',lupri)
   IF(nPasses.GT.1)CALL LSQUIT('ERROR IN DERLONDON_ETIJ: passes not working in comb. with magderiv and .ETUVIL',-1)
   CALL DERLONDON_ETIJ(ETIJ,nprim,l1+l2+2*der+magder,l1+der+magder,l2+der,nX,l1,l2,P%nPrimitives,p%orbital1%center)
ENDIF

CALL LS_DZERO(Ecoeffs,Nijk*ntuv*nprim)

!Derivative settings
ncent = 1 + der
IF (P%single) ncent = 1

#ifdef VAR_LSDEBUGINT
IF(size(P%preExpFac).NE.P%nPrimAlloc*P%nAngAlloc)call lsquit('PreExpFac dim error5',-1)
#endif
pref => P%preExpFac
IF (der.GT. 0) THEN
  CALL mem_alloc(pref,P%nPrimAlloc*P%nAngmom)
ELSE
  pref => P%preExpFac
ENDIF

offset = (iAngmom-1)*P%nPrimAlloc
#ifdef VAR_LSDEBUGINT
IF(iAngmom.GT.P%nAngAlloc)call lsquit('PreExpFac dim error6',-1)
#endif
ijk=0
!Derivative loop
DO icent=1,ncent
! Nuclear derivative settings (using Hermite gaussians according to PCCP,2007,9,4771-4779)
! Set up increments
  increment(1) = der+1-icent
  increment(2) = icent - 1
  IF (P%single) THEN
    IF (P%orbital1%TYPE_empty) THEN
      increment(1) = 0
      increment(2) = der
    ELSEIF (P%orbital2%TYPE_empty) THEN
      increment(1) = der
      increment(2) = 0
    ELSE
      CALL LSQUIT('Error in BuildEcoeffTensor_PA. overlap%single and not orbital%empty',lupri)
    ENDIF
  ENDIF
  nmstart=0
  nmend=0
  IF(magder.GT.0)THEN
     increment(1) = 0
     increment(2) = 0
     nmstart=1
     nmend=3
  ENDIF
! Set up prefactors that include the (2a)^derA * (2b)^derB
  IF (der.GT. 0) THEN
    DO iPrimP=1,nprim
      iPrim1 = P%iprim1(iPrimP)
      iPrim2 = P%iprim2(iPrimP)
      pref(iPrimP+offset) = P%preExpFac(iPrimP+offset)*(2E0_realk*P%orbital1%exponents(iPrim1))**increment(1)*&
           &               (2E0_realk*P%orbital2%exponents(iPrim2))**increment(2)
    ENDDO
  ENDIF

! Regular loop starts here
  DO Xdir = nmstart,nmend
     X1=1
     X2=2
     X3=3
     IF(Xdir.EQ.1)X1=4
     IF(Xdir.EQ.2)X2=5
     IF(Xdir.EQ.3)X3=6
     NT=0
     NU=0
     NV=0
     IF(Xdir.EQ.1)NT=1
     IF(Xdir.EQ.2)NU=1
     IF(Xdir.EQ.3)NV=1    
     DO P2 = 0,l2+increment(2)
        DO iP2=P2,0,-1
           DO jP2=P2-iP2,0,-1
              kP2=P2-iP2-jP2
              DO P1 = 0,l1+increment(1)
                 DO iP1=P1,0,-1
                    DO jP1=P1-iP1,0,-1
                       kP1=P1-iP1-jP1
                       IF(iP1+jP1+kP1+iP2+jP2+kP2 .EQ. l1+l2+der) then
                          ijk=ijk+1
                          DO tP=0,iP1+iP2+NT
                             DO uP=0,jP1+jP2+NU
                                DO vP=0,kP1+kP2+NV
                                   IF (tP+uP+vP .GE. startAng) THEN
                                      IF(MOD(tP+uP+vP,2).EQ. 0)THEN
                                         signijk = D1
                                      ELSE
                                         signijk = signP
                                      ENDIF
                                      ituvP=TUV%TUVindex(tP,uP,vP)-ioff
                                      DO iPrimP=1,nprim
                                         Ecoeffs(iPrimP,ituvP,ijk) = &
                                              &ETIJ(iPrimP,tP,iP1,iP2,X1)&
                                              &*ETIJ(iPrimP,uP,jP1,jP2,X2)&
                                              &*ETIJ(iPrimP,vP,kP1,kP2,X3)&
                                              &*pref(iPrimP+offset)*signijk
                                      ENDDO
                                   ENDIF
                                ENDDO
                             ENDDO
                          ENDDO
                       ENDIF
                    ENDDo
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
ENDDO
IF (der.GT. 0) THEN
! CALL PrintTensor(Ecoeffs,'E-coefficients-build',&
!      &Nijk,ntuv,nPrim,2,'ijk   ','tuv   ','prim  ',3)
  CALL mem_dealloc(pref)
ENDIF

call mem_dealloc(ETIJ)

END SUBROUTINE BuildEcoeffTensor_PA

!> \brief derivative london ecoefficients placed after std ETIJ
!> \author T. Kjaergaard
!> \date 2010
!> \param ETIJ the ecoefficients
!> \param dim1 first dimension of ETIJ
!> \param dim2 first dimension of ETIJ
!> \param dim3 first dimension of ETIJ
!> \param dim4 first dimension of ETIJ
!> \param dim5 first dimension of ETIJ
!> \param l1 angular moment of orbital 1
!> \param l2 angular moment of orbital 2
!> \param nPrimitives the number of primitive functions
!> \param center of orbital1 (X,Y,Z) coordinates
SUBROUTINE DERLONDON_ETIJ(DERLONETIJ,dim1,dim2,dim3,dim4,dim5,l1,l2,nPrimitives,center)
implicit none
integer :: dim1,dim2,dim3,dim4,dim5,l1,l2,nPasses,nPrimitives
real(realk)  :: DERLONETIJ(dim1,0:dim2,0:dim3,0:dim4,dim5),center(3)
!
integer :: Xdir,iP1,iP2,tP,IprimP,IpassP
DO Xdir=1,3
 DO iP1 = 0,l1
  DO iP2 = 0,l2
   do tP = 0,iP1+iP2
    DO IprimP = 1,nPrimitives
     DERLONETIJ(iPrimP,tP,iP1,iP2,Xdir+3)=DERLONETIJ(iPrimP,tP,iP1+1,iP2,Xdir)+center(Xdir)*DERLONETIJ(iPrimP,tP,iP1,iP2,Xdir) 
    ENDDO
   ENDDO
   tP=iP1+iP2+1
   DO IprimP = 1,nPrimitives
      DERLONETIJ(iPrimP,tP,iP1,iP2,Xdir+3)=DERLONETIJ(iPrimP,tP,iP1+1,iP2,Xdir)
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE DERLONDON_ETIJ

!> \brief build Ecoefficient tensor (ordering nPrim,Nijk,ntuv)
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param TUV The TUV indeces
!> \param P the overlap distribution P or Q
!> \param signP the sign (1 for P and -1 for Q)
!> \param Ecoeffs the ecoefficients to be built
!> \param Nijk the number of cartesian angular components
!> \param ntuv the number of hermite angular components
!> \param nPrim the number of primitives
!> \param iAngmom the angular moment index
!> \param nPasses the number of passes
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE BuildEcoeffTensor_PA_old(TUV,P,signP,Ecoeffs,Nijk,ntuv,nprim,iAngmom,nPasses,LUPRI,IPRINT)
implicit none
TYPE(TUVitem)      :: TUV
TYPE(Overlap)      :: P
Integer            :: iAngmom,LUPRI,IPRINT,nPasses
Integer            :: Nijk,ntuv,nprim
Real(realk)        :: signP
!Allocatable for ECOEFFS - Dimensions are nPrimitives, Zero to Sum of Max Angmom 1 and 2, 
!                          Zero to Max Angmom 1, Zero to Max Angmom 2, 3 (for X,Y,Z)
real(realk),pointer  :: ETIJ(:,:,:,:,:) 
integer             :: P1,P2,tP,uP,vP
integer             :: ituvP,iP1,jP1,kP1,iP2,jP2,kP2,iPrimP
integer             :: l1,l2,ijk,offset
Real(realk)         :: Ecoeffs(nPrim,ntuv,Nijk)
Real(realk)         :: signijk

IF (P%endGeoOrder.GT. 0) CALL LSQUIT('BuildEcoeffTensor_PA_old for endGeoOrder>0',lupri)

l1 = P%orbital1%angmom(P%indexAng1(iangmom))
l2 = P%orbital2%angmom(P%indexAng2(iangmom))

call mem_alloc(ETIJ,nprim,l1+l2,l1,l2,3,.FALSE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.)

CALL GET_ECOEFF(ETIJ,nprim,l1+l2,l1,l2,P%nPrimitives,nPasses,P,LUPRI,IPRINT)

CALL LS_DZERO(Ecoeffs,nPrim*Nijk*ntuv)
offset = (iAngmom-1)*P%nPrimAlloc
#ifdef VAR_LSDEBUGINT
IF(size(P%preExpFac).NE.P%nPrimAlloc*P%nAngAlloc)call lsquit('PreExpFac dim error6',-1)
IF(iAngmom.GT.P%nAngAlloc)call lsquit('PreExpFac dim error6',-1)
#endif
ijk=0
DO P2 = 0,l2
   DO iP2=P2,0,-1
      DO jP2=P2-iP2,0,-1
         kP2=P2-iP2-jP2
         DO P1 = 0,l1
            DO iP1=P1,0,-1
               DO jP1=P1-iP1,0,-1
                  kP1=P1-iP1-jP1
                  IF(iP1+jP1+kP1+iP2+jP2+kP2 .EQ. l1+l2) then
                  ijk=ijk+1
                     DO tP=0,iP1+iP2
                        DO uP=0,jP1+jP2
                           DO vP=0,kP1+kP2
                              signijk = signP**(tP+uP+vP)
                              ituvP=TUV%TUVindex(tP,uP,vP)
                              DO iPrimP=1,nPrim!P%nPrimitives
                                 Ecoeffs(iPrimP,ituvP,ijk) = &
                                      &ETIJ(iPrimP,tP,iP1,iP2,1)&
                                      &*ETIJ(iPrimP,uP,jP1,jP2,2)&
                                      &*ETIJ(iPrimP,vP,kP1,kP2,3)&
                                      &*P%preExpFac(iPrimP+offset)*signijk
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDo
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
call mem_dealloc(ETIJ)

END SUBROUTINE BuildEcoeffTensor_PA_old

!> \brief print Ecoefficients
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param EcoeffCont the Ecoefficient to be printed
!> \param nAng the number of angular components
!> \param nOrb the number of orbital components
!> \param nPrim the number of primitives
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE PrintContractEcoeff(EcoeffCont,nAng,nOrb,nPrim,LUPRI,IPRINT)
implicit none
Integer     :: nAng,nOrb,nPrim,LUPRI,IPRINT
Real(realk) :: EcoeffCont(nAng,nOrb,nPrim)
!
Integer :: iAng,iOrb,iPrim
IF (IPRINT.GT. 50) THEN
  WRITE(LUPRI,'(3X,A)') '***************************************************************'
  WRITE(LUPRI,'(3X,A)') '***                       ContractEcoeff'
  WRITE(LUPRI,'(3X,A)') '***************************************************************'
  DO iPrim=1,nPrim
    DO iOrb=1,nOrb
        WRITE(LUPRI,'(5X,A,I3,A,I3)') 'iPrim =',iPrim,', iOrb =',iOrb
        WRITE(LUPRI,'(7X,5ES10.4/,(7X,5ES10.4))') (EcoeffCont(iAng,iOrb,iPrim),iAng=1,nAng)
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE PrintContractEcoeff

!> \brief wrapper routine for spherical transformation of integrals
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Spherical the Spherical transformation matrix 
!> \param ijk1 the number of cartesian angular components of orbital 1
!> \param ijk2 the number of cartesian angular components of orbital 2
!> \param lm1 the number of spherical angular components of orbital 1
!> \param lm2 the number of spherical angular components of orbital 2
!> \param ang1 the angular momentum for orbital 1
!> \param ang2 the angular momentum for orbital 2
!> \param P the overlap distribution either P or Q
!> \param integralNonSp the input array
!> \param integralSpher the output array
!> \param dim1 the dimension of the array to be transformed
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SphericalTransformation(TUV,Spherical,ijk1,ijk2,lm1,lm2,ang1,ang2,&
    &         integralNonSp,integralSpher,dim1,lupri,iprint,nOperatorComp)
implicit none
TYPE(TUVitem)      :: TUV
Integer            :: LUPRI,IPRINT,dim1
Integer             :: lm1,lm2,ijk1,ijk2,ang1,ang2,nOperatorComp
Real(realk)         :: Spherical(ijk1*ijk2,lm1*lm2)
Real(realk)         :: integralNonSp(*)
Real(realk)         :: integralSpher(*)

! Sets up spherical transformation matrix
CALL Buildsphericaltransformation(Spherical,TUV%SPH_MAT(ang1+1)%elms,&
     &              TUV%SPH_MAT(ang2+1)%elms,lm1,lm2,ijk1,ijk2,LUPRI,IPRINT)

CALL SphericalTransformation_PA(Spherical,ijk1*ijk2,lm1*lm2,integralNonSp,integralSpher,&
     &                            dim1,lupri,iprint,nOperatorComp)
!  
END SUBROUTINE SphericalTransformation

!> \brief spherical transformation of integrals order PA
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Spherical the Spherical transformation matrix 
!> \param ijk the number of cartesian angular components
!> \param lm the number of spherical angular components
!> \param integralNonSp the input array
!> \param integralSpher the output array
!> \param dim1 the dimension of the array to be transformed
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SphericalTransformation_PA(Spherical,ijk,lm,integralNonSp,integralSpher,&
     &                                dim1,lupri,iprint,nOperatorComp)
implicit none
Integer            :: LUPRI,IPRINT
Integer            :: lm,ijk, dim1,nOperatorComp
Real(realk)        :: Spherical(ijk,lm)
Real(realk)        :: integralNonSp(dim1,ijk,nOperatorComp)
Real(realk)        :: integralSpher(dim1,lm,nOperatorComp)
!
real(realk),parameter :: Zero=1.0E-18_realk,D0=0E0_realk,D1=1.0E0_realk
integer :: ilm,ang,i,X
real(realk) :: tmp
do X=1,nOperatorComp
   do ilm = 1,lm
      DO i=1,dim1
         integralSpher(i,ilm,X) = D0
      ENDDO
      do ang = 1,ijk
         tmp = Spherical(ang,ilm)
         IF(ABS(tmp).GT.zero)THEN
            DO i=1,dim1
               integralSpher(i,ilm,X) = integralSpher(i,ilm,X) + tmp*integralNonSp(i,ang,X)
            ENDDO
         ENDIF
      enddo
   enddo
enddo
IF(IPRINT.GT.100)THEN
   do X=1,nOperatorComp
      WRITE(lupri,'(A,I3,A,I3,A,I3)')&
           &'The SphericalTransformed intermidiates lm=',lm,'dim1=',dim1,'X=',X
      call output(integralSpher(:,:,X),1,lm,1,dim1,lm,dim1,1,lupri)
   enddo
ENDIF

END SUBROUTINE SphericalTransformation_PA

!> \brief spherical transformation of orbital 1
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Spher1 the Spherical transformation matrix from orbital1
!> \param ijk1 the number of cartesian angular components of orbital 1
!> \param ijk2 the number of cartesian angular components of orbital 2
!> \param lm1 the number of spherical angular components of orbital 1
!> \param lm2 the number of spherical angular components of orbital 2
!> \param ijk the number of cartesian angular components of overlap
!> \param lm the number of spherical angular components of overlap
!> \param integralNonSp the input array
!> \param integralSpher the output array
!> \param dim1 the unaffected dimension of the array to be transformed
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SphericalTransformation1(Spher1,ijk1,ijk2,lm1,lm2,ijk,lm,integralNonSp,&
     &integralSpher,dim1,lupri,iprint)
implicit none
Integer            :: LUPRI,IPRINT
Integer            :: ijk1,ijk2,lm1,lm2,lm,ijk,dim1
Real(realk)        :: Spher1(lm1,ijk1)
Real(realk)        :: integralNonSp(dim1,ijk1,ijk2)
Real(realk)        :: integralSpher(dim1,lm1,ijk2)
!
Real(realk),parameter        :: D0=0.0E0_realk
integer :: i,indijk1,indlm1,indijk2
Real(realk) :: TMP

DO i = 1,dim1
 DO indijk2=1,ijk2
  DO indlm1=1,lm1
   TMP=D0
   DO indijk1=1,ijk1
      TMP = TMP + Spher1(indlm1,indijk1)*integralNonSp(i,indijk1,indijk2)
   ENDDO
   integralSpher(i,indlm1,indijk2) = TMP
  ENDDO
 ENDDO
ENDDO

END SUBROUTINE SphericalTransformation1

!> \brief spherical transformation of orbital 2
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Spher2 the Spherical transformation matrix from orbital1
!> \param ijk1 the number of cartesian angular components of orbital 1
!> \param ijk2 the number of cartesian angular components of orbital 2
!> \param lm1 the number of spherical angular components of orbital 1
!> \param lm2 the number of spherical angular components of orbital 2
!> \param ijk the number of cartesian angular components of overlap
!> \param lm the number of spherical angular components of overlap
!> \param integralNonSp the input array
!> \param integralSpher the output array
!> \param dim1 the unaffected dimension of the array to be transformed
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE SphericalTransformation2(Spher2,ijk1,ijk2,lm1,lm2,ijk,lm,integralNonSp,&
     &integralSpher,dim1,lupri,iprint)
implicit none
Integer            :: LUPRI,IPRINT
Integer            :: ijk1,ijk2,lm1,lm2,lm,ijk,dim1
Real(realk)        :: Spher2(lm2,ijk2)
Real(realk)        :: integralNonSp(dim1,lm1,ijk2)
Real(realk)        :: integralSpher(dim1,lm1,lm2)
!
Real(realk),parameter        :: D0=0.0E0_realk
integer :: i,indlm2,indlm1,indijk2
real(realk) :: TMP

DO i = 1,dim1
 DO indlm1=1,lm1
  DO indlm2=1,lm2
   TMP=D0
   DO indijk2=1,ijk2
      TMP = TMP + Spher2(indlm2,indijk2)*integralNonSp(i,indlm1,indijk2)
   ENDDO
   integralSpher(i,indlm1,indlm2) = TMP
  ENDDO
 ENDDO
ENDDO

END SUBROUTINE SphericalTransformation2

!> \brief construct the spherical transformation matrix
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Spherical the Spherical transformation matrix 
!> \param Spher1 the Spherical transformation matrix for orbital 1 
!> \param Spher2 the Spherical transformation matrix for orbital 2
!> \param lm1 the number of spherical angular components for orbital 1 
!> \param lm2 the number of spherical angular components for orbital 2 
!> \param ijk1 the number of cartesian angular components for orbital 1 
!> \param ijk2 the number of cartesian angular components for orbital 2 
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE Buildsphericaltransformation(Spherical,Spher1,Spher2,lm1,lm2,&
     &                                      ijk1,ijk2,LUPRI,IPRINT)
Integer,intent(in)     :: ijk1,ijk2,lm1,lm2
Real(realk),intent(inout) :: Spherical(ijk1,ijk2,lm1,lm2)
Real(realk),intent(in) :: Spher1(lm1,ijk1)
Real(realk),intent(in) :: Spher2(lm2,ijk2)
!
Integer     :: indijk1,indijk2,indlm1,indlm2
!
DO indijk2=1,ijk2
   DO indijk1=1,ijk1
      DO indlm2=1,lm2
         DO indlm1=1,lm1
            Spherical(indijk1,indijk2,indlm1,indlm2) = Spher1(indlm1,indijk1)*Spher2(indlm2,indijk2)

         ENDDO
      ENDDO
   ENDDO
ENDDO
CALL PrintSphericalTransformation(Spherical,lm1,lm2,ijk1,ijk2,LUPRI,IPRINT)
END SUBROUTINE Buildsphericaltransformation

!> \brief print the spherical transformation matrix
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param Spherical the Spherical transformation matrix 
!> \param lm1 the number of spherical angular components for orbital 1 
!> \param lm2 the number of spherical angular components for orbital 2 
!> \param ijk1 the number of cartesian angular components for orbital 1 
!> \param ijk2 the number of cartesian angular components for orbital 2 
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE PrintSphericalTransformation(Spherical,lm1,lm2,ijk1,ijk2,LUPRI,IPRINT)
Real(realk) :: Spherical(ijk1,ijk2,lm1,lm2)
Integer     :: ijk1,ijk2,lm1,lm2
!
Integer     :: indijk1,indijk2,indlm1,indlm2
!
IF (IPRINT.GT. 50) THEN
  CALL LSHEADER(LUPRI,'SphericalTransformation')
  DO indlm2=1,lm2
    DO indlm1=1,lm1
      WRITE(LUPRI,'(5X,A,I3,A,I3)') 'lm1 =',indlm1,' lm2 =',indlm2
      DO indijk1=1,ijk1
        WRITE(LUPRI,'(5X,5ES16.4)') &
       &       (Spherical(indijk1,indijk2,indlm1,indlm2),indijk2=1,ijk2)
      ENDDO
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE PrintSphericalTransformation

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
SUBROUTINE ConstructContraction_PA(CC,P,nP,nC1,nC2,iA1,iA2,LUPRI,IPRINT)
Implicit none
Type(Overlap) :: P
Integer       :: iA1,iA2,nP,nC1,nC2,LUPRI,IPRINT
Real(realk)   :: CC(nC1,nC2,nP)
!
Integer       :: iP,i1,i2,iC1,iC2,nP1,nP2

nP1 = P%orbital1%nPrimitives
nP2 = P%orbital2%nPrimitives
DO iP=1,nP
   i1 = P%iprim1(iP)
   i2 = P%iprim2(iP)
   DO iC2=1,nC2
      DO iC1=1,nC1
         CC(iC1,iC2,iP)=P%orbital1%CC(iA1)%p%elms(i1+(iC1-1)*nP1)*P%orbital2%CC(iA2)%p%elms(i2+(iC2-1)*nP2) 
      ENDDO
   ENDDO
ENDDO
!IF (IPRINT.GT. 30) THEN
!   CALL PrintTensor(CC,'Overlap CC          ',nP,nC1,nC2,Lupri,&
!        & 'prim  ','C1    ','C2    ',1)
!ENDIF

END SUBROUTINE ConstructContraction_PA

!> \brief Print the contraction matrix
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param CC the contraction matrix
!> \param nP the number of primitives
!> \param nC1 the number of contracted functions on orbital 1
!> \param nC2 the number of contracted functions on orbital 2
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE PrintOverlapContraction(CC,nP,nC1,nC2,LUPRI,IPRINT)
Implicit none
Integer       :: nP,nC1,nC2,LUPRI,IPRINT
Real(realk)   :: CC(nP,nC1,nC2)
!
Integer       :: iP,iC1,iC2
IF (IPRINT.GT. 30) THEN
  WRITE(LUPRI,'(5X,A)') '***************************************************'
  WRITE(LUPRI,'(5X,A)') '***                  Overlap CC'
  WRITE(LUPRI,'(5X,A)') '***************************************************'
  WRITE(LUPRI,'(7X,A)') 'C1  C2 CC(prim12)'
  DO iC2=1,nC2
    DO iC1=1,nC1
      WRITE(LUPRI,'(5X,2I4,5ES10.4/,(13X,5ES10.4))') iC1,iC2,(CC(iP,iC1,iC2),iP=1,nP)
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE PrintOverlapContraction

!> \brief build the single cartesian Ecoeff (ETIJ)
!> \author A. Teale
!> \date 2009
!> \param ETIJ the single cartesian Ecoeff to be built 
!> \param nprimpass the number of primitives * number of npasses
!> \param maxang12 the combined maximum angular momentum
!> \param maxang1 the maximum angular momentum for orbital 1
!> \param maxang2 the maximum angular momentum for orbital 2
!> \param nprim the number of primitives
!> \param npass the number of passes
!> \param P the overlap distribution P or Q
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
SUBROUTINE GET_ECOEFF(ETIJ,nprimpass,maxang12,maxang1,maxang2,nprim,npass,P,LUPRI,IPRINT)
implicit none
REAL(REALK),PARAMETER    :: D1 = 1.0E0_realk, D2 = 2.0E0_realk, DHALF = 0.5E0_realk
integer                  :: IPRINT, LUPRI, I, J, K, CARTDIR
integer                  :: nprim,maxang12,maxang1,maxang2
TYPE(Overlap)            :: P
integer                  :: nprimpass,npass
real(realk) :: PINV(nPrimPass), APINV(nPrimPass), BPINV(nPrimPass), HPINV(nPrimPass), PA(nPrimPass,3), PB(nPrimPass,3)
real(realk) :: TWOA(nPrimPass), TWOB(nPrimPass)
!Allocatable for ECOEFFS - Dimensions are No. Pairs, Zero to Sum of Max Angmom 1 and 2, 
!                          Zero to Max Angmom 1, Zero to Max Angmom 2, 3 (for X,Y,Z)
!real(realk),allocatable  :: ETIJ(:,:,:,:,:) 
real(realk)  :: ETIJ(nprimpass,0:maxang12,0:maxang1,0:maxang2,3),Xdist,Ydist,Zdist 

!IF (IPRINT.GT. 50) CALL PRN_ECINFO(P,LUPRI)

IF (.NOT.P%type_Empty.AND..NOT.P%type_Nucleus)THEN
   DO K=1,nPrimPass
      I=P%iprim1(K)
      J=P%iprim2(K)
      PINV(K)  = D1/(P%orbital1%exponents(I)+P%orbital2%exponents(J)) ! 1/p
      APINV(K) = P%orbital1%exponents(I) * PINV(K)      ! a/p
      BPINV(K) = P%orbital2%exponents(J) * PINV(K)      ! b/p
      HPINV(K) = DHALF*PINV(K)                     ! 1/2p
      TWOA(K)  = D2*P%orbital1%exponents(I)             ! 2a
      TWOB(K)  = D2*P%orbital2%exponents(J)             ! 2b
   ENDDO
   DO I = 1,nPass
      Xdist=P%distance12(1+(I-1)*3)
      Ydist=P%distance12(2+(I-1)*3)
      Zdist=P%distance12(3+(I-1)*3)
      DO K=1,nPrim
         J = K+(I-1)*nPrim
         !CAMT SETUP REQUIRED DISTANCES
         PA(J,1) = -BPINV(J)*Xdist
         PA(J,2) = -BPINV(J)*Ydist
         PA(J,3) = -BPINV(J)*Zdist
         PB(J,1) =  APINV(J)*Xdist
         PB(J,2) =  APINV(J)*Ydist
         PB(J,3) =  APINV(J)*Zdist
      END DO
   ENDDO
ENDIF

!AMT THEN CALL ROUTINES FOR CARTESIAN OR HERMITE ECOEFFS 
!AMT LOOP OVER CARTESIAN DIRECTONS (1,2,3 -> X,Y,Z)
DO CARTDIR = 1,3
      IF (P%type_hermite .OR. P%type_hermite_single) THEN
           CALL HERM_ECOEFFS(ETIJ(1,0,0,0,CARTDIR),nPrimPass,maxAng1,&
          &             maxAng2,PA(:,CARTDIR),PB(:,CARTDIR),HPINV,TWOA,TWOB,LUPRI)
      ELSE IF (P%type_cartesian .OR. P%type_cartesian_single) THEN
           CALL CART_ECOEFFS(ETIJ(1,0,0,0,CARTDIR),nPrimPass,maxAng1,&
          &             maxAng2,PA(:,CARTDIR),PB(:,CARTDIR),HPINV,LUPRI)
      ELSE IF (P%type_Empty) THEN
           DO i=1,nPrimPass
              ETIJ(i,0,0,0,CARTDIR) = D1
           ENDDO
      ELSE IF (P%type_Nucleus) THEN
           DO i=1,nPrimPass
              ETIJ(i,0,0,0,CARTDIR) = D1
           ENDDO
      ELSE
           WRITE(LUPRI,*)'ERROR : UNRECOGNISED TYPE OF OVERLAP!!!'
      ENDIF
      IF (IPRINT.GT. 20) CALL PRN_ECOEFFS(ETIJ,nPrimPass,maxAng1,&
           &            maxAng2,CARTDIR,LUPRI)
ENDDO
END SUBROUTINE GET_ECOEFF

!> \brief SUBROUTINE TO CALCULATE E's USING THE USUAL CARTESIAN RECURANNCE RELATIONS - SEE P354 OF BOOK
!> \author A. Teale
!> \date 2009
!>
!> SEE ERIECF ROUTINE for F77 EQUIVALENT
!>
!> VARIABLES
!> =========      
!> nPrimPairs NUMBER OF PRIMITIVE PAIRS IN OVERLAP DISTRIBUTION
!> MAXI,MAXJ  MAXIMUM I,J VALUES IN ETIJ's (MAX ANGMOM OF G_i and G_j RESPECTIVELY)
!> I,J,K,IJ   COUNTERS     
!> PA         DISTANCE BETWEEN CENTRE OF CHARGE P AND A (THE CENTRE OF PRIMITIVE GAUSSIAN G_i)
!> PB         DISTANCE BETWEEN CENTRE OF CHARGE P AND B (THE CENTRE OF PRIMITIVE GAUSSIAN G_j)
!> HPINV      1/2p (p = a + b I.E. SUM OF EXPONENTS OF PRIM G_i and G_j)           
!> PINV       1/p
!> APINV      a/p
!> BPINV      b/p
!> ETIJ       E COEFFICIENTS 
!> X          INTEGER SPECIFYING X,Y,Z (1,2,3) COMPONENT OF OVERALL OVERLAP DISTRIBUTION
!> T          COUNTER INDEX USED TO COUNT UPTO SUM OF ANGULAR MOMENTA OF G_i AND G_j 
!>
!> \param ETIJ the single cartesian Ecoeff to be built 
!> \param nprimpass the number of primitives * number of npasses
!> \param MAXI the maximum angular momentum for orbital 1
!> \param MAXJ the maximum angular momentum for orbital 2
!> \param PA the argument -b/p*Rpq (b is exponent on orbital 2, p=a+b)
!> \param PB the argument  a/p*Rpq (a is exponent on orbital 1, p=a+b)
!> \param HPINV the argument 1/2p
!> \param LUPRI the logical unit number for the output file
SUBROUTINE CART_ECOEFFS(ETIJ,nPrimPairs,MAXI,MAXJ,PA,PB,HPINV,LUPRI)
implicit none
REAL(REALK),PARAMETER  ::  D1 = 1.0E0_realk, D2 = 2.0E0_realk
INTEGER                ::  I,J,K
INTEGER                ::  T,nPrimPairs,MAXI,MAXJ,IJ,LUPRI
REAL(REALK)            ::  T1
REAL(REALK)            ::  PA(nPrimPairs), PB(nPrimPairs) 
REAL(REALK)            ::  ETIJ(nPrimPairs,0:MAXI+MAXJ,0:MAXI,0:MAXJ)
REAL(REALK)            ::  HPINV(nPrimPairs)


!C
!C Run over I (J = 0)                                     
!C ==================
!C
   DO I = 0, MAXI                                          
     IF (I .LE. 2) THEN
      CALL ECOEF_ICASES(nPrimPairs,MAXI,MAXJ,ETIJ,PA,HPINV,I)
     ELSE                                                 
!C                                                               
!C            E(I,0)                                            
!C                                                               
      DO K = 1, nPrimPairs                                   
        ETIJ(K,  0,I,0) = PA(K)*ETIJ(K,  0,I-1,0)+ ETIJ(K,  1,I-1,0)            
        ETIJ(K,I-1,I,0) = HPINV(K)*ETIJ(K,I-2,I-1,0)+ PA(K)*ETIJ(K,I-1,I-1,0)        
        ETIJ(K,  I,I,0) = HPINV(K)*ETIJ(K,I-1,I-1,0)           
      END DO                                             
      DO T = 1, I - 2                                    
!       T1 = DFLOAT(T + 1)                           
        T1 = T + 1
        DO K = 1, nPrimPairs                                
          ETIJ(K,T,I,0) = HPINV(K)*ETIJ(K,T-1,I-1,0)+ PA(K)*ETIJ(K,  T,I-1,0)+ T1*ETIJ(K,T+1,I-1,0)    
        END DO                                         
      END DO                                            
     END IF                                               
!C                                                               
!C   Run over J                                          
!C   ==========                                          
!C                                                               
     DO J = 1, MAXJ                                      
       IJ = I + J                                        
       IF (IJ .LE. 2) THEN
         CALL ECOEF_IJCASES(nPrimPairs,MAXI,MAXJ,ETIJ,PA,PB,HPINV,I,IJ)
       ELSE                                              
!C                                                                
!C             E(I,J)                                            
!C                                                               
         DO K = 1, nPrimPairs                                
           ETIJ(K,   0,I,J) = PB(K)*ETIJ(K,   0,I,J-1)+ ETIJ(K,   1,I,J-1)   
           ETIJ(K,IJ-1,I,J) = HPINV(K)*ETIJ(K,IJ-2,I,J-1)+ PB(K)*ETIJ(K,IJ-1,I,J-1)   
           ETIJ(K,  IJ,I,J) = HPINV(K)*ETIJ(K,IJ-1,I,J-1)  
         END DO                                         
         DO T = 1, IJ - 2                               
!          T1 = DFLOAT(T + 1)                       
           T1 = T + 1
           DO K = 1, nPrimPairs                            
             ETIJ(K,T,I,J)=HPINV(K)*ETIJ(K,T-1,I,J-1)+ PB(K)*ETIJ(K,  T,I,J-1)+ T1*ETIJ(K,T+1,I,J-1)   
           END DO                                      
         END DO                                        
       END IF                                            
     END DO                                               
   END DO                                                  
RETURN
END SUBROUTINE CART_ECOEFFS

!> \brief helper routine for building the ETIJ, J=0  and I .LE. 2
!> \author A. Teale
!> \date 2009
!>
!> \param nprimpass the number of primitives * number of npasses
!> \param MAXI the maximum angular momentum for orbital 1
!> \param MAXJ the maximum angular momentum for orbital 2
!> \param ETIJ the single cartesian Ecoeff to be built 
!> \param PA the argument -b/p*Rpq (b is exponent on orbital 2, p=a+b)
!> \param HPINV the argument 1/2p
!> \param LUPRI the logical unit number for the output file
SUBROUTINE ECOEF_ICASES(nPrimPairs,MAXI,MAXJ,ETIJ,PA,HPINV,I)
implicit none
REAL(REALK),PARAMETER::  D1 = 1.0E0_realk, D2 = 2.0E0_realk
INTEGER     :: I,nPrimPairs,MAXI,MAXJ,K
REAL(REALK) :: PA(nPrimPairs),HPINV(nPrimPairs)
REAL(REALK) :: ETIJ(nPrimPairs,0:MAXI+MAXJ,0:MAXI,0:MAXJ)
!C                                                                               
!C           E(0,0)                                              
!C
  IF (I .EQ. 0) THEN                               
    DO K = 1, nPrimPairs                               
      ETIJ(K,0,0,0) = D1                          
    END DO
  ELSE IF (I .EQ. 1) THEN                          
!C                                                             
!C           E(1,0)                                            
!C
    DO K = 1, nPrimPairs                               
      ETIJ(K,0,1,0) = PA(K)                       
      ETIJ(K,1,1,0) = HPINV(K)                      
    END DO                                         
  ELSE IF (I .EQ. 2) THEN                          
!C                                                
!C           E(2,0)                                              
!C                                                              
    DO K = 1, nPrimPairs                                
      ETIJ(K,0,2,0) = PA(K)*PA(K) + HPINV(K)     
      ETIJ(K,1,2,0) = D2*PA(K)*HPINV(K)             
      ETIJ(K,2,2,0) = HPINV(K)*HPINV(K)            
    END DO                                         
  ENDIF
END SUBROUTINE ECOEF_ICASES 

!> \brief helper routine for building the ETIJ
!> \author A. Teale
!> \date 2009
!>
!> \param nprimpass the number of primitives * number of npasses
!> \param MAXI the maximum angular momentum for orbital 1
!> \param MAXJ the maximum angular momentum for orbital 2
!> \param ETIJ the single cartesian Ecoeff to be built 
!> \param PA the argument -b/p*Rpq (b is exponent on orbital 2, p=a+b)
!> \param PB the argument  a/p*Rpq (a is exponent on orbital 1, p=a+b)
!> \param HPINV the argument 1/2p
!> \param I counter
!> \param IJ counter
SUBROUTINE ECOEF_IJCASES(nPrimPairs,MAXI,MAXJ,ETIJ,PA,PB,HPINV,I,IJ)
implicit none
REAL(REALK),PARAMETER::  D1 = 1.0E0_realk, D2 = 2.0E0_realk
INTEGER     :: I,IJ,nPrimPairs,MAXI,MAXJ,K
REAL(REALK) :: PA(nPrimPairs),PB(nPrimPairs),HPINV(nPrimPairs)
REAL(REALK) :: ETIJ(nPrimPairs,0:MAXI+MAXJ,0:MAXI,0:MAXJ)
  IF (IJ .EQ. 1) THEN                               
!C                                                                 
!C              E(0,1)                                           
!C                                                               
    DO K = 1, nPrimPairs                             
      ETIJ(K,0,0,1) = PB(K)                    
      ETIJ(K,1,0,1) = HPINV(K)                   
    END DO                                      
  ELSE IF (IJ .EQ. 2) THEN                         
!C                                                
!C                 E(0,2)                        
!C                                              
    IF (I .EQ. 0) THEN                          
      DO K = 1, nPrimPairs                          
        ETIJ(K,0,0,2) = PB(K)*PB(K) + HPINV(K)    
        ETIJ(K,1,0,2) = D2*PB(K)*HPINV(K)            
        ETIJ(K,2,0,2) = HPINV(K)*HPINV(K)               
      END DO                                   
    ELSE                                        
!C                                                               
!C                    E(1,1)                                    
!C                                                             
      DO K = 1, nPrimPairs                          
        ETIJ(K,0,1,1) = PA(K)*PB(K) + HPINV(K)    
        ETIJ(K,1,1,1) = (PA(K) + PB(K))*HPINV(K)      
        ETIJ(K,2,1,1) = HPINV(K)*HPINV(K)               
      END DO                                  
    END IF                                     
  ENDIF 
END SUBROUTINE ECOEF_IJCASES

!> \brief print the ecoeff tensor
!> \author A. Teale
!> \date 2009
!>
!> \param ETIJ the single cartesian Ecoeff to be built 
!> \param nprimpass the number of primitives * number of npasses
!> \param MAXI the maximum angular momentum for orbital 1
!> \param MAXJ the maximum angular momentum for orbital 2
!> \param X an INTEGER SPECIFYING X,Y,Z (1,2,3) COMPONENT OF OVERALL OVERLAP DISTRIBUTION
!> \param LUPRI the logical unit number for the output file
SUBROUTINE PRN_ECOEFFS(ETIJ,nPrimPairs,MAXI,MAXJ,X,LUPRI)
implicit none
INTEGER          ::  MAXI,MAXJ,X,I,J,T,K,nPrimPairs,LUPRI
CHARACTER(len=4) ::  WORD
REAL(REALK)      :: ETIJ(nPrimPairs,0:MAXI+MAXJ,0:MAXI,0:MAXJ)
!C
!C     *********************************
!C     ***** PRINT E COEFFICIENTS  *****
!C     *********************************
!C
         WRITE(LUPRI,*)'Output from LSint ECFs'
         WRITE(LUPRI,*)'----------------------'
         WRITE (LUPRI,'(2X,A,2I5)') 'MAXI,MAXJ   ', MAXI, MAXJ 
            IF (X .EQ. 1) WORD = 'EX00'
            IF (X .EQ. 2) WORD = 'EY00'
            IF (X .EQ. 3) WORD = 'EZ00'
            DO I = 0, MAXI 
            DO J = 0, MAXJ 
            DO T = 0, I + J
               WRITE (LUPRI,'(/,2X,A4,A1,I1,A1,I1,A1,I1,A1,/)')&
     &              WORD, '(', I, ',', J, ', ',  T, ')'
               WRITE (LUPRI,'(1X,6ES16.8)') (ETIJ(K,T,I,J),K=1,nPrimPairs)
            END DO
            END DO
            END DO
END SUBROUTINE PRN_ECOEFFS

!> \brief print the ecoeff info
!> \author A. Teale
!> \date 2009
!>
!> \param P the overlap distribution P or Q
!> \param LUPRI the logical unit number for the output file
SUBROUTINE PRN_ECINFO(P,LUPRI)
implicit none
integer                  :: LUPRI, I
TYPE(Overlap)            :: P

WRITE(LUPRI,*)'HERE IS THE RELEVANT INFORMATION FROM THE P'
WRITE(LUPRI,*)'-------------------------------------------'
WRITE(LUPRI,*)'nAngmom',P%nAngmom
WRITE(LUPRI,*)'nPrimitives',P%nPrimitives
WRITE(LUPRI,*)'maxContracted',P%maxContracted
WRITE(LUPRI,*)'Distance Between A and B:'
DO I=1,P%npasses
   WRITE(LUPRI,*)'Pass nr:',I
   WRITE(LUPRI,*)'X',P%distance12(1+(I-1)*3)
   WRITE(LUPRI,*)'Y',P%distance12(2+(I-1)*3)
   WRITE(LUPRI,*)'Z',P%distance12(3+(I-1)*3)
enddo
WRITE(LUPRI,*)'Squared Distance',P%squaredDistance

WRITE(LUPRI,*)'HERE IS THE INFORMATION RELEVANT FOR ORBITAL 1'
WRITE(LUPRI,*)'----------------------------------------------'
WRITE(LUPRI,*)'TYPE hermite  ',P%orbital1%TYPE_Hermite
WRITE(LUPRI,*)'TYPE empty    ',P%orbital1%TYPE_Empty
WRITE(LUPRI,*)'TYPE cartesian',P%orbital1%TYPE_Cartesian
WRITE(LUPRI,*)'TYPE nucleus  ',P%orbital1%TYPE_Nucleus
WRITE(LUPRI,*)'TYPE pCharge  ',P%orbital1%TYPE_pCharge

WRITE(LUPRI,*)'SPHERICAL ?',P%orbital1%Spherical
WRITE(LUPRI,*)'Maximum Angular Momentum',P%orbital1%maxAngmom
WRITE(LUPRI,*)'Angular Momentum',P%orbital1%nAngmom
WRITE(LUPRI,*)'Number of Primitives',P%orbital1%nPrimitives
WRITE(LUPRI,*)'Max Contracted',P%orbital1%maxContracted
WRITE(LUPRI,*)'Centre of orbital 1:'
WRITE(LUPRI,*)'X',P%orbital1%center(1)
WRITE(LUPRI,*)'Y',P%orbital1%center(2)
WRITE(LUPRI,*)'Z',P%orbital1%center(3)
WRITE(LUPRI,*)'Primitive Exponents in Orbital1'
WRITE(LUPRI,*)'-------------------------------'
DO I=1,P%orbital1%nPrimitives
   WRITE(LUPRI,*)'Primitive',I,P%orbital1%exponents(I)
ENDDO


WRITE(LUPRI,*)'HERE IS THE INFORMATION RELEVANT FOR ORBITAL 2'
WRITE(LUPRI,*)'----------------------------------------------'
WRITE(LUPRI,*)'TYPE hermite  ',P%orbital2%TYPE_Hermite
WRITE(LUPRI,*)'TYPE empty    ',P%orbital2%TYPE_Empty
WRITE(LUPRI,*)'TYPE cartesian',P%orbital2%TYPE_Cartesian
WRITE(LUPRI,*)'TYPE nucleus  ',P%orbital2%TYPE_Nucleus
WRITE(LUPRI,*)'TYPE pCharge  ',P%orbital2%TYPE_pCharge

WRITE(LUPRI,*)'SPHERICAL ?',P%orbital2%Spherical
WRITE(LUPRI,*)'Maximum Angular Momentum',P%orbital2%maxAngmom
WRITE(LUPRI,*)'Angular Momentum',P%orbital2%nAngmom
WRITE(LUPRI,*)'Number of Primitives',P%orbital2%nPrimitives
WRITE(LUPRI,*)'Max Contracted',P%orbital2%maxContracted
WRITE(LUPRI,*)'Centre of orbital 1:'
WRITE(LUPRI,*)'X',P%orbital2%center(1)
WRITE(LUPRI,*)'Y',P%orbital2%center(2)
WRITE(LUPRI,*)'Z',P%orbital2%center(3)
WRITE(LUPRI,*)'Primitive Exponents in Orbital2'
WRITE(LUPRI,*)'-------------------------------'
DO I=1,P%orbital2%nPrimitives
   WRITE(LUPRI,*)'Primitive',I,P%orbital2%exponents(I)
ENDDO
END SUBROUTINE PRN_ECINFO

!> \brief SUBROUTINE TO CALCULATE E's USING THE HERMITE RECURANNCE RELATIONS
!> \author A. Teale
!> \date 2009
!>
!> VARIABLES
!> =========      
!> nPrimPairs   NUMBER OF PRIMITIVE PAIRS IN OVERLAP DISTRIBUTION
!> MAXI,MAXJ    MAXIMUM I,J VALUES IN ETIJ's (MAX ANGMOM OF G_i and G_j RESPECTIVELY)
!> I,J,K,IJ     COUNTERS     
!> PA           DISTANCE BETWEEN CENTRE OF CHARGE P AND A (THE CENTRE OF PRIMITIVE GAUSSIAN G_i)           
!> PB           DISTANCE BETWEEN CENTRE OF CHARGE P AND B (THE CENTRE OF PRIMITIVE GAUSSIAN G_j)
!> HPINV        1/2p (p = a + b I.E. SUM OF EXPONENTS OF PRIM G_i and G_j)           
!> PINV         1/p
!> APINV, BPINV a/p b/p
!> TWOA, TWOB   2a, 2b
!> ETIJ         E COEFFICIENTS 
!> X            INTEGER SPECIFYING X,Y,Z (1,2,3) COMPONENT OF OVERALL OVERLAP DISTRIBUTION
!> T            COUNTER INDEX USED TO COUNT UPTO SUM OF ANGULAR MOMENTA OF G_i AND G_j 
!>
!> \param ETIJ the single cartesian Ecoeff to be built 
!> \param nprimpass the number of primitives * number of npasses
!> \param MAXI the maximum angular momentum for orbital 1
!> \param MAXJ the maximum angular momentum for orbital 2
!> \param PA the argument -b/p*Rpq (b is exponent on orbital 2, p=a+b)
!> \param PB the argument  a/p*Rpq (a is exponent on orbital 1, p=a+b)
!> \param HPINV the argument 1/2p
!> \param TWOA 2a (a is exponent on orbital 1)
!> \param TWOB 2b (b is exponent on orbital 2)
!> \param LUPRI the logical unit number for the output file
SUBROUTINE HERM_ECOEFFS(ETIJ,nPrimPairs,MAXI,MAXJ,PA,PB,&
                      & HPINV,TWOA,TWOB,LUPRI)
implicit none
REAL(REALK),PARAMETER  ::  D1 = 1.0E0_realk, D2 = 2.0E0_realk
INTEGER                ::  I,J,K
INTEGER                ::  T,nPrimPairs,MAXI,MAXJ,IJ,LUPRI
REAL(REALK)            ::  T1, TIM, TJM
REAL(REALK)            ::  PA(nPrimPairs), PB(nPrimPairs)
REAL(REALK)            ::  ETIJ(nPrimPairs,0:MAXI+MAXJ,0:MAXI,0:MAXJ)
REAL(REALK)            ::  HPINV(nPrimPairs)
REAL(REALK)            ::  TWOA(nPrimPairs),TWOB(nPrimPairs)
!C
!C Run over I (J = 0)                                     
!C ==================
!C
   DO I = 0, MAXI
!    TIM = DFLOAT(I-1)
     TIM = I-1
     IF (I .LE. 2) THEN
      CALL ECOEF_ICASES_HERM(nPrimPairs,MAXI,MAXJ,ETIJ,PA,HPINV,TWOA,I)
     ELSE
!C                                                               
!C            E(I,0)                                            
!C                                                               
      DO K = 1, nPrimPairs
        ETIJ(K,  0,I,0) = PA(K)*ETIJ(K,  0,I-1,0)+ ETIJ(K,  1,I-1,0)&
                        & - TIM*ETIJ(K, 0,I-2,0)/TWOA(K)
        ETIJ(K,I-1,I,0) = HPINV(K)*ETIJ(K,I-2,I-1,0)+ PA(K)*ETIJ(K,I-1,I-1,0)
        ETIJ(K,  I,I,0) = HPINV(K)*ETIJ(K,I-1,I-1,0)
      END DO                                   
      DO T = 1, I - 2                                    
!       T1 = DFLOAT(T + 1)
        T1 = T + 1
        DO K = 1, nPrimPairs                                
          ETIJ(K,T,I,0) = HPINV(K)*ETIJ(K,T-1,I-1,0)+ PA(K)*ETIJ(K,  T,I-1,0)&
                        & + T1*ETIJ(K,T+1,I-1,0) - TIM*ETIJ(K, T,I-2,0)/TWOA(K)
        END DO
      END DO         
     END IF             
!C                                                               
!C   Run over J                                          
!C   ==========                                          
!C                                                               
     DO J = 1, MAXJ
       IJ = I + J
!      TJM = DFLOAT(J-1)
       TJM = J-1

       IF (IJ .LE. 2) THEN
         CALL ECOEF_IJCASES_HERM(nPrimPairs,MAXI,MAXJ,ETIJ,PA,PB,HPINV,TWOB,I,IJ)
       ELSE
!C                                                                
!C             E(I,J)                                            
!C                                                               
          IF(J.LT. 2)THEN
             DO K = 1, nPrimPairs
                ETIJ(K,   0,I,J) = PB(K)*ETIJ(K,   0,I,J-1)+ ETIJ(K,   1,I,J-1)
                ETIJ(K,IJ-1,I,J) = HPINV(K)*ETIJ(K,IJ-2,I,J-1)+ PB(K)*ETIJ(K,IJ-1,I,J-1)
                ETIJ(K,  IJ,I,J) = HPINV(K)*ETIJ(K,IJ-1,I,J-1)
             END DO
             DO T = 1, IJ - 2
                T1 = T + 1
                DO K = 1, nPrimPairs
                   ETIJ(K,T,I,J)= HPINV(K)*ETIJ(K,T-1,I,J-1)+ PB(K)*ETIJ(K,  T,I,J-1)&
                        & + T1*ETIJ(K,T+1,I,J-1)! - TJM*ETIJ(K,  T,I,J-2)/TWOB(K)
                END DO
             END DO
          ELSE
             DO K = 1, nPrimPairs
                ETIJ(K,   0,I,J) = PB(K)*ETIJ(K,   0,I,J-1)+ ETIJ(K,   1,I,J-1)&
                     & -TJM*ETIJ(K,  0,I,J-2)/TWOB(K)
                ETIJ(K,IJ-1,I,J) = HPINV(K)*ETIJ(K,IJ-2,I,J-1)+ PB(K)*ETIJ(K,IJ-1,I,J-1)
                ETIJ(K,  IJ,I,J) = HPINV(K)*ETIJ(K,IJ-1,I,J-1)
             END DO
             DO T = 1, IJ - 2
                T1 = T + 1
                DO K = 1, nPrimPairs
                   ETIJ(K,T,I,J)= HPINV(K)*ETIJ(K,T-1,I,J-1)+ PB(K)*ETIJ(K,  T,I,J-1)&
                        & + T1*ETIJ(K,T+1,I,J-1) - TJM*ETIJ(K,  T,I,J-2)/TWOB(K)
                END DO
             END DO
          ENDIF
       END IF
     END DO
   END DO
RETURN
END SUBROUTINE HERM_ECOEFFS

!> \brief helper routine for building the hermite ETIJ
!> \author A. Teale
!> \date 2009
!>
!> \param nprimpass the number of primitives * number of npasses
!> \param MAXI the maximum angular momentum for orbital 1
!> \param MAXJ the maximum angular momentum for orbital 2
!> \param ETIJ the single cartesian Ecoeff to be built 
!> \param PA the argument -b/p*Rpq (b is exponent on orbital 2, p=a+b)
!> \param HPINV the argument 1/2p
!> \param TWOA 2a (a is exponent on orbital 1)
!> \param I counter
SUBROUTINE ECOEF_ICASES_HERM(nPrimPairs,MAXI,MAXJ,ETIJ,PA,HPINV,TWOA,I)
implicit none
REAL(REALK),PARAMETER::  D1 = 1.0E0_realk, D2 = 2.0E0_realk
INTEGER     :: I,nPrimPairs,MAXI,MAXJ,K 
REAL(REALK) :: PA(nPrimPairs),HPINV(nPrimPairs)
REAL(REALK) :: TWOA(nPrimPairs) 
REAL(REALK) :: ETIJ(nPrimPairs,0:MAXI+MAXJ,0:MAXI,0:MAXJ)
!C                                                                               
!C           E(0,0)                                              
!C
  IF (I .EQ. 0) THEN
    DO K = 1, nPrimPairs
      ETIJ(K,0,0,0) = D1
    END DO
  ELSE IF (I .EQ. 1) THEN
!C                                                             
!C           E(1,0)                                            
!C
    DO K = 1, nPrimPairs
      ETIJ(K,0,1,0) = PA(K)
      ETIJ(K,1,1,0) = HPINV(K)
    END DO
  ELSE IF (I .EQ. 2) THEN
!C                                                
!C           E(2,0)                                              
!C                                                              
    DO K = 1, nPrimPairs
      ETIJ(K,0,2,0) = PA(K)*PA(K) + HPINV(K) - D1/TWOA(K)
      ETIJ(K,1,2,0) = D2*PA(K)*HPINV(K)
      ETIJ(K,2,2,0) = HPINV(K)*HPINV(K)
    END DO
  ENDIF
END SUBROUTINE ECOEF_ICASES_HERM

!> \brief second helper routine for building the hermite ETIJ
!> \author A. Teale
!> \date 2009
!>
!> \param nprimpass the number of primitives * number of npasses
!> \param MAXI the maximum angular momentum for orbital 1
!> \param MAXJ the maximum angular momentum for orbital 2
!> \param ETIJ the single cartesian Ecoeff to be built 
!> \param PA the argument -b/p*Rpq (b is exponent on orbital 2, p=a+b)
!> \param PB the argument  a/p*Rpq (a is exponent on orbital 1, p=a+b)
!> \param HPINV the argument 1/2p
!> \param TWOB 2b (b is exponent on orbital 2)
!> \param I counter
!> \param IJ counter
SUBROUTINE ECOEF_IJCASES_HERM(nPrimPairs,MAXI,MAXJ,ETIJ,PA,PB,HPINV,TWOB,I,IJ)
implicit none
REAL(REALK),PARAMETER::  D1 = 1.0E0_realk, D2 = 2.0E0_realk
INTEGER     :: I,IJ,nPrimPairs,MAXI,MAXJ,K
REAL(REALK) :: PA(nPrimPairs),PB(nPrimPairs),HPINV(nPrimPairs)
REAL(REALK) :: TWOB(nPrimPairs)
REAL(REALK) :: ETIJ(nPrimPairs,0:MAXI+MAXJ,0:MAXI,0:MAXJ)
  IF (IJ .EQ. 1) THEN
!C                                                                 
!C              E(0,1)                                           
!C                                                               
    DO K = 1, nPrimPairs
      ETIJ(K,0,0,1) = PB(K)    
      ETIJ(K,1,0,1) = HPINV(K) 
    END DO
  ELSE IF (IJ .EQ. 2) THEN     
!C                                                
!C                 E(0,2)                        
!C                                              
    IF (I .EQ. 0) THEN
      DO K = 1, nPrimPairs     
        ETIJ(K,0,0,2) = PB(K)*PB(K) + HPINV(K) - D1/TWOB(K)
        ETIJ(K,1,0,2) = D2*PB(K)*HPINV(K)
        ETIJ(K,2,0,2) = HPINV(K)*HPINV(K)
      END DO
    ELSE
!C                                                               
!C                    E(1,1)                                    
!C                                                             
      DO K = 1, nPrimPairs     
        ETIJ(K,0,1,1) = PA(K)*PB(K) + HPINV(K)
        ETIJ(K,1,1,1) = (PA(K) + PB(K))*HPINV(K)
        ETIJ(K,2,1,1) = HPINV(K)*HPINV(K)
      END DO
    END IF
  ENDIF
END SUBROUTINE ECOEF_IJCASES_HERM

!*******************************************************************************************
! ANDY END OF ECOEFFS ?
!*******************************************************************************************

!> \brief set the FTUV batch
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param F The FTUV overlap
!> \param NFTUVbatches the number of FTUV batches
!> \param OD The OD-batch
!> \param Q The right or left hand side overlap distribution
!> \param Input The integral input
!> \param sharedTUV The TUV indeces
!> \param Integral The integral specifications
!> \param Alloc Information about the maximal sizes for allocation purposes
!> \param ndmat The number of density matrix
!> \param maxPrimPass maximum number of primitive*npasses
!> \param lupri The default print unit
!> \param iprint The print level (the higher the more information)
SUBROUTINE SET_FTUVbatches(F,NFTUVbatches,OD,Q,Input,SharedTUV,Integral,Alloc,ndmat,maxPrimPass,LUPRI,IPRINT)
implicit none
TYPE(INTEGRALINPUT)   :: INPUT
TYPE(Overlap),pointer :: F(:)
TYPE(Overlap),pointer :: Q(:)
Integer               :: NFTUVbatches,maxPrimPass(20),LUPRI,IPRINT
TYPE(ODITEM)          :: OD 
TYPE(Integralitem)    :: Integral
TYPE(Allocitem)    :: Alloc
TYPE(TUVitem)         :: SharedTUV
!
Integer :: I,J,K,nUnique,ftuvindex,np,ndmat
Integer  ,pointer :: ODtoFTUVindex(:),nPrimOD(:),FTUVtoPassIndex(:)
Integer(KIND=8),pointer :: Identifier(:), UniqeIdentfiers(:)
Integer  ,pointer :: FTUVprim(:),FTUVPassType(:)
Integer  ,pointer :: FTUVminAng(:),FTUVmaxAng(:),FTUVntuv(:),FTUVprimPass(:)
Integer  ,pointer :: offPrim(:)
Logical :: unique,uniquePassType
Integer :: nTUV,nPassTypes,iPassType
Integer :: CSscreenLOG,tenminusmaxgab

call mem_alloc(ODtoFTUVindex,OD%nbatches)
call mem_alloc(nPrimOD,OD%nbatches)
call mem_alloc(Identifier,OD%nbatches)
call mem_alloc(UniqeIdentfiers,OD%nbatches)
call mem_alloc(FTUVprim,OD%nbatches)
call mem_alloc(FTUVPassType,OD%nbatches)
call mem_alloc(FTUVminAng,OD%nbatches)
call mem_alloc(FTUVmaxAng,OD%nbatches)
call mem_alloc(FTUVntuv,OD%nbatches)
call mem_alloc(FTUVprimPass,OD%nbatches)
call mem_alloc(offPrim,OD%nbatches)

nUnique = 0
nPassTypes=0
DO I = 1,OD%nbatches
 CALL SET_Overlap(Q(I),Input,SharedTUV,Integral,OD%BATCH(I),2,LUPRI,IPRINT,.TRUE.)
 np = Q(I)%nprimitives
 nPrimOD(I) = np
!Simen  Add center information to this identifier - maybe some integer value of 
!Simen  a weighted OD-center divided by BOXSIZE. Note the larger the boxsize, the 
!Simen  less effective CS-screening becomes, whereas the FTUV-contraction step becomes
!Simen  faster. I other words there is a tradeoff between two opposing effects here.
!Simen  Note also that Maxgab should be inherited from OD-batches to FTUV-batches
 IF (Q(I)%maxAngmom.GT. 98) CALL LSQUIT('Non-unique angmom-identifier in SET_FTUVbatches',lupri)
 IF (INPUT%CS_SCREEN) THEN
    tenminusmaxgab = 10-Q(i)%maxGab
    CSscreenLOG = max(0,tenminusmaxgab)
 ELSE
    CSscreenLOG = 0
 ENDIF
 IF (CSscreenLOG.GT. 998) CALL LSQUIT('Non-unique CSscreenLOG-identifier in SET_FTUVbatches',lupri)
 Identifier(I) = (CSscreenLOG+1)*10000+(Q(I)%maxAngmom+1)*100 + (Q(I)%minAngmom+1)
 unique = .true.
 uniquePassType = .true.
 DO J=nUnique,1,-1
    IF (Identifier(I).EQ.UniqeIdentfiers(J)) THEN
      iPassType = FTUVpassType(J)
      uniquePassType = .false.
      IF ((FTUVprim(J)+np).GT.maxPrimPass(Q(I)%endAngmom+1)) EXIT
      ftuvIndex = J
      unique = .false.
      exit
    ENDIF
 ENDDO
 IF (uniquePassType) THEN
    nPassTypes = nPassTypes + 1
    iPassType  = nPassTypes
 ENDIF
 IF (unique) THEN
    nUnique = nUnique + 1
    ftuvIndex = nUnique
    UniqeIdentfiers(nUnique) = Identifier(I)
    FTUVpassType(ftuvIndex) = iPassType
    FTUVprim(ftuvIndex) = Q(I)%nPrimitives
    FTUVminAng(ftuvIndex) = Q(I)%startAngmom
    FTUVmaxAng(ftuvIndex) = Q(I)%endAngmom
    FTUVprimPass(ftuvIndex) = maxPrimPass(Q(I)%endAngmom+1)
    FTUVntuv(ftuvIndex)   = (Q(I)%endAngmom+1)*(Q(I)%endAngmom+2)*(Q(I)%endAngmom+3)/6 &
     &                     - Q(I)%startAngmom*(Q(I)%startAngmom+1)*(Q(I)%startAngmom+2)/6
 ELSE
    FTUVprim(ftuvIndex) = FTUVprim(ftuvIndex) + Q(I)%nPrimitives
 ENDIF
 ODtoFTUVindex(I) = ftuvIndex
ENDDO
NFTUVbatches = nUnique

Alloc%maxPrimRHS    = 0
Alloc%maxPrimTUVRHS = 0
DO I=1,NFTUVbatches
  np =  FTUVprim(I)
  nTUV =  FTUVnTUV(I)
  Alloc%maxPrimRHS    = max(np,Alloc%maxPrimRHS,FTUVprimPass(I))
  Alloc%maxPrimTUVRHS = max(np*nTUV,FTUVprimPass(I)*nTUV,Alloc%maxPrimTUVRHS)
ENDDO

call mem_dealloc(FTUVprimPass)

!Initialize FTUV-batches
call INIT_BUFCOUNTERS(4)
DO iPassType=1,nPassTypes
  DO J=1,NFTUVbatches
    IF (FTUVpassType(J).EQ.iPassType) THEN
      CALL memFTUVbatches(FTUVprim(J),FTUVntuv(J),ndmat,4)
    ENDIF
  ENDDO
ENDDO
CALL ALLOC_ODFTUV_BUFFERS

call mem_alloc(FTUVtoPassIndex,NFTUVbatches)
call mem_alloc(F,NFTUVbatches)
I = 0
DO iPassType=1,nPassTypes
  DO J=1,NFTUVbatches
    IF (FTUVpassType(J).EQ.iPassType) THEN
      I = I + 1
      FTUVtoPassIndex(J) = I
      CALL InitFTUVbatches(F(I),iPassType,FTUVprim(J),FTUVminAng(J),FTUVmaxAng(J),FTUVntuv(J),ndmat,4)
      offPrim(I) = 0
    ENDIF
  ENDDO
ENDDO

!Add OD-batches to corresponding FTUV-batches
DO I=1,OD%nbatches
 J = ODtoFTUVindex(I)
 IF (J.GT. 0) THEN
  K = FTUVtoPassIndex(J)
  np = nPrimOD(I)
  CALL AddODtoFTUV(F(K),np,offPrim(K),Q(I),ndmat,LUPRI,IPRINT)
  CALL FREE_OVERLAP(Q(I))
  offPrim(K) = offPrim(K) + np
 ENDIF
ENDDO

call mem_dealloc(ODtoFTUVindex)
call mem_dealloc(nPrimOD)
call mem_dealloc(FTUVtoPassIndex)
call mem_dealloc(Identifier)
call mem_dealloc(UniqeIdentfiers)
call mem_dealloc(FTUVprim)
call mem_dealloc(FTUVPassType)
call mem_dealloc(FTUVminAng)
call mem_dealloc(FTUVmaxAng)
call mem_dealloc(FTUVntuv)
call mem_dealloc(offPrim)

END SUBROUTINE SET_FTUVbatches

!> \brief add Overlap distribution to FTUV
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param F The FTUV overlap
!> \param np number of primitive functions 
!> \param offp an offset 
!> \param Q The right or left hand side overlap distribution
!> \param ndmat The number of density matrix
!> \param lupri The default print unit
!> \param iprint The print level (the higher the more information)
SUBROUTINE AddODtoFTUV(F,np,offp,Q,ndmat,LUPRI,IPRINT)
implicit none
TYPE(Overlap) :: F,Q
Integer       :: np,offp,ndmat,LUPRI,IPRINT
!
Integer :: ip,idir,idmat,tuv,offsetF,offsetQ
IF (F%nTUV.NE.Q%nTUV) THEN
  CALL LSQUIT('Programming error: nTUV maismatch in AddODtoFTUV',lupri)
ENDIF
F%maxGab = max(INT(F%maxGab),INT(Q%maxGab))
F%ODextent = Q%ODextent
F%ODcenter = Q%ODcenter
F%nPrimitives = F%nPrimitives + np
#ifdef VAR_LSDEBUGINT
IF(size(Q%preExpFac).NE.Q%nPrimAlloc*Q%nAngAlloc)call lsquit('PreExpFac dim error9',-1)
IF(size(F%preExpFac).NE.F%nPrimAlloc*F%nAngAlloc)call lsquit('PreExpFac dim error10',-1)
IF(Q%nPrimAlloc.GT.F%nPrimAlloc)call lsquit('PreExpFac dim error11A',-1)
IF(size(F%FTUV).NE.F%nPrimAlloc*F%nTUVAlloc*ndmat)call lsquit('FTUV dim error9',-1)
IF(size(Q%FTUV).NE.Q%nPrimAlloc*Q%nTUVAlloc*ndmat)call lsquit('FTUV dim error10',-1)
IF(size(Q%FTUV).GT.size(F%FTUV))call lsquit('FTUV dim error11',-1)
#endif
DO ip=1,np
   F%preExpFac(ip+offp) = Q%preExpFac(ip)
   offsetF = (ip+offp-1)*3
   offsetQ = (ip-1)*3
   DO idir=1,3
      F%center(idir+offsetF) = Q%center(idir+offsetQ)
   ENDDO
   F%reducedExponents(ip+offp) = Q%reducedExponents(ip)
   F%exponents(ip+offp) = Q%exponents(ip)
ENDDO
DO idmat = 1,ndmat
   DO tuv=1,F%nTUV
      offsetQ = (TUV-1)*Q%nPrimAlloc+(idmat-1)*Q%nPrimAlloc*Q%nTUV
      offsetF = (TUV-1)*F%nPrimAlloc+(idmat-1)*F%nPrimAlloc*F%nTUV+offp
      DO ip=1,np
         F%FTUV(ip+offsetF) = Q%FTUV(ip+offsetQ)
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE AddODtoFTUV

!> \brief init the FTUVbatches
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param F The FTUV overlap
!> \param iPasstype the type of this overlap
!> \param nPrim number of primitives
!> \param minAng the minimum angular momentum
!> \param maxAng the maximum angular momentum
!> \param ntuv number of angular components
!> \param ndmat The number of density matrix
SUBROUTINE InitFTUVbatches(F,iPassType,nPrim,minAng,maxAng,nTUV,ndmat,ielec)
implicit none
TYPE(Overlap) :: F
Integer       :: nPrim,minAng,maxAng,nTUV,ndmat,iPassType,ielec
CALL ALLOC_ORBITAL(F%orbital1,1,1,1,1,.TRUE.,ielec)
CALL ALLOC_ORBITAL(F%orbital2,1,1,1,1,.TRUE.,ielec)
call mem_ODpointer_alloc(F%Orb1atom,1,ielec)
call mem_ODpointer_alloc(F%Orb1mol,1,ielec)
call mem_ODpointer_alloc(F%Orb1batch,1,ielec)
call mem_ODpointer_alloc(F%Orb2atom,1,ielec)
call mem_ODpointer_alloc(F%Orb2mol,1,ielec)
call mem_ODpointer_alloc(F%Orb2batch,1,ielec)
CALL ALLOC_OVERLAP(F,nprim,1,1,.TRUE.,nTUV,0,ndmat,ielec)

F%TYPE_Empty = .FALSE.
F%TYPE_Hermite = .FALSE.
F%TYPE_Hermite_single = .FALSE.
F%TYPE_Cartesian = .FALSE.
F%TYPE_Cartesian_single = .FALSE.
F%TYPE_Nucleus = .FALSE.
F%TYPE_FTUV = .TRUE.

F%single = .FALSE.

F%passType                 = iPassType
F%nAngmom                  = 1
F%orbital1%nAngmom         = 1
F%orbital2%nAngmom         = 1
F%nPrimitives              = 0
F%nPasses                  = 1
F%minAngmom                = minAng
F%maxAngmom                = maxAng
F%startAngmom              = minAng
F%startGeoOrder            = 0
F%endGeoOrder              = 0
F%ngeoderivcomp            = 1
F%nCartesianMomentComp     = 1
F%endAngmom                = maxAng
F%nTUV                     = nTUV
F%totOrbitals(:)           = ndmat
F%angmom(1)                = maxAng
F%orbital1%angmom(1)       = 0
F%orbital1%nContracted(1)  = 1
F%orbital1%nOrbComp(1)     = 1
F%orbital1%startOrbital(1) = 1
F%orbital1%totOrbitals     = ndmat
F%orbital2%angmom(1)       = 0
F%orbital2%nContracted(1)  = 1
F%orbital2%nOrbComp(1)     = 1
F%orbital2%startOrbital(1) = 1
F%orbital2%totOrbitals     = 1
F%indexAng1(1)             = 1
F%indexAng2(1)             = 1
F%nOrbitals(1)             = 1
F%nOrbComp(1)              = 1
F%nContracted(1)           = 1
F%maxGab                   = shortzero
F%mmorder                  = 0
F%do_mulmom                = .false.
F%magderiv                 = 0

END SUBROUTINE InitFTUVbatches

!> \brief init the FTUVbatches
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param F The FTUV overlap
!> \param nPrim number of primitives
!> \param ntuv number of angular components
!> \param ndmat The number of density matrix
SUBROUTINE MemFTUVbatches(nPrim,nTUV,ndmat,ielec)
implicit none
Integer,intent(in) :: nPrim,nTUV,ndmat,ielec
integer(kind=long) :: nrealk,nint
nrealk=0
nint=0
CALL MEM_ORBITAL(1,1,1,1,.TRUE.,nrealk,nint)
CALL MEM_ORBITAL(1,1,1,1,.TRUE.,nrealk,nint)
nint = nint + 6
CALL MEM_OVERLAP1(nprim,1,1,.TRUE.,nTUV,0,ndmat,nrealk,nint)
call ADD_BUFCOUNTERS(ielec,nint,nrealk)
END SUBROUTINE MemFTUVbatches

!> \brief select overlap distribution pass types from OD Batch
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param ODpassesIndex for a given batch this gives the passtype 
!> \param nPrim the number of primitives
!> \param ODB the overlap batch
!> \param nbatches the number of batches in ODB
!> \param nPassTypes the number of passtypes 
!> \param maxpasses the number of maximum passes 
!> \param the number of passes for each type
!> \param INPUT the integral input specifications, contains all info about what to do
!> \param iprint The print level (the higher the more information)
!> \param lupri The default print unit
!> \param side LHS or RHS specification
SUBROUTINE SelectODTypesFromODbatch(ODtypeIndex,Alloc,ODB,nBatches,&
     &nODTypes,TypeOverlapindex,INPUT,IPRINT,LUPRI,SIDE)
use memory_handling
implicit none
type(allocitem) :: alloc
Integer       :: nBatches,nODTypes,IPRINT,LUPRI
Integer :: ODTYPEIndex(nBatches)
type(integerpointer),pointer :: TypeOverlapindex(:)
TYPE(ODITEM)        :: ODB
TYPE(IntegralInput) :: Input
Character*(*)       :: side
!
Integer       :: I,J,numAng1,numAng2,IODtype,N,dim,maxbat
Character(len=42),pointer  :: UniqeIdentfiers(:)
integer(kind=long),pointer  :: UniqeIdentfiers2(:) 
Integer       :: angmomIdentifier1,angmomIdentifier2,CCidentifier1,CCidentifier2
Logical       :: spher1,spher2,type_hermite_single,spherE,unique,sameAO,single
Logical       :: type_empty,LHS
integer :: NgeoDERIVcomp,NMOM,maxijk,totorbitals,nETUV,i1,i2,start2
integer :: l1,l2,l,ijk1,ijk2,ijk,ijkcart,nTUV,maxangmom,minangmom,ndim5
integer :: startangmom,endangmom,np,nprimp,CMorder,nCartesianMomentComp
Type(lstensor),pointer :: GAB
Integer :: batchA,batchB,atomA,atomB,Gindex,nprim1,nprim2,iprim12,iDprim12,nprim
integer :: geoderivorder,nOrb1,maxnprim,maxncont,maxangid,IntsameAO
integer :: maxCCidentifier,magderiv,nOperatorComp,s1,s2,n1
integer(kind=short) :: maxgab,maxelm
IF (Side.EQ.'LHS') THEN
   !We only devide the LHS overlaps into passes for allocation purposes 
   !so they do not need to have the same primitives and so on. 
   CMorder = 0
   LHS = .TRUE.
   call mem_alloc(UniqeIdentfiers,nBatches)
   nODTypes = 0
   DO I=1,nBatches
      nPrim = ODB%BATCH(I)%AO(1)%p%nPrimitives*ODB%BATCH(I)%AO(2)%p%nPrimitives
      spher1 = .FALSE.
      spher2 = .FALSE.
      IF(.NOT.ODB%BATCH(I)%AO(1)%p%TYPE_Empty)THEN
         spher1        = ODB%BATCH(I)%AO(1)%p%spherical
      ENDIF
      IF(.NOT.ODB%BATCH(I)%AO(2)%p%TYPE_Empty)THEN
         spher2        = ODB%BATCH(I)%AO(2)%p%spherical
      ENDIF
      type_hermite_single = input%hermiteEcoeff.AND.(.NOT.ODB%BATCH(I)%AO(1)%p%type_empty.AND.ODB%BATCH(I)%AO(2)%p%type_empty)&
           & .OR. (ODB%BATCH(I)%AO(1)%p%type_empty.AND..NOT.ODB%BATCH(I)%AO(2)%p%type_empty)
      spherE = Input%sphericalEcoeff      
      sameAO = ODB%BATCH(I)%sameAO
      IF (nprim.GT. 9999) THEN
         CALL LSQUIT('nPrim1>9999 in SelectODTypes',lupri)
      ELSEIF (ODB%BATCH(I)%redTYPE.GT. 99999999) THEN
         CALL LSQUIT('redtype>99999999 in SelectODTypes',lupri)
      ENDIF
!$OMP CRITICAL (ifortwrite)
      WRITE(UniqeIdentfiers(I),'(I11,I8,5L3)')nprim,&
           &ODB%BATCH(I)%redTYPE,spher1,spher2,type_hermite_single,spherE,sameAO
!$OMP END CRITICAL (ifortwrite)
      unique = .TRUE.
      DO J=I-1,1,-1
         IF (UniqeIdentfiers(I).EQ.UniqeIdentfiers(J)) THEN
            ODTYPEIndex(I)=ODTYPEIndex(J)
            unique = .FALSE.
            EXIT
         ENDIF
      ENDDO
      IF (unique) THEN
         nODTypes = nODTypes + 1
         ODTYPEIndex(I) = nODTypes
      ENDIF
   ENDDO
ELSEIF (Side.EQ.'RHS') THEN
   CMorder = INPUT%CMorder
   LHS = .FALSE.
   GAB => INPUT%LST_GAB_RHS
   MAXELM=INPUT%PS_MAXELM_LHS
   call mem_alloc(UniqeIdentfiers,nBatches)
   nODTypes = 0
   DO I=1,nBatches
      nprim1 = ODB%BATCH(I)%AO(1)%p%nPrimitives
      nprim2 = ODB%BATCH(I)%AO(2)%p%nPrimitives
      iprim12 = 0
      iDprim12 = 0
      nPrim = ODB%BATCH(I)%AO(1)%p%nPrimitives*ODB%BATCH(I)%AO(2)%p%nPrimitives
      IF (INPUT%PS_SCREEN.AND.(.NOT.(ODB%BATCH(I)%AO(1)%p%type_empty .AND.ODB%BATCH(I)%AO(2)%p%type_empty))) THEN
         AtomA = ODB%BATCH(I)%AO(1)%p%atom
         AtomB = ODB%BATCH(I)%AO(2)%p%atom
         BatchA = ODB%BATCH(I)%AO(1)%p%batch
         BatchB = ODB%BATCH(I)%AO(2)%p%batch
         Gindex = GAB%INDEX(atomA,atomB,1,1)
         maxBat = GAB%SLSAO(Gindex)%maxBat
         s1 = GAB%SLSAO(Gindex)%startLocalOrb(BatchA)-1
         s2 = GAB%SLSAO(Gindex)%startLocalOrb(BatchB+maxBat)-1
         n1 = GAB%SLSAO(Gindex)%nLocal(1)
         DO i2=1,nprim2
            DO i1=1,nprim1
               maxGab = GAB%SLSAO(Gindex)%selms(s1+i1 + (s2+i2-1)*n1)
               IF( MAXGAB .GT. INPUT%PS_THRLOG-MAXELM)THEN                  
                  iPrim12 = iPrim12 + 1
                  iDprim12 = iDprim12 + i1 + i2*nprim
               ENDIF
            ENDDO
         ENDDO
      ELSE
         iprim12 = nPrim
         DO i2=1,ODB%BATCH(I)%AO(2)%p%nPrimitives
            DO i1=1,ODB%BATCH(I)%AO(1)%p%nPrimitives
               iDprim12 = iDprim12 + i1 + i2*nprim
            ENDDO
         ENDDO
      ENDIF
      nprim = iprim12
      spher1 = .FALSE.
      spher2 = .FALSE.
      IF(.NOT.ODB%BATCH(I)%AO(1)%p%TYPE_Empty)THEN
         spher1        = ODB%BATCH(I)%AO(1)%p%spherical
      ENDIF
      IF(.NOT.ODB%BATCH(I)%AO(2)%p%TYPE_Empty)THEN
         spher2        = ODB%BATCH(I)%AO(2)%p%spherical
      ENDIF
      type_hermite_single = input%hermiteEcoeff.AND.&
           & (.NOT.ODB%BATCH(I)%AO(1)%p%type_empty.AND.ODB%BATCH(I)%AO(2)%p%type_empty)&
           & .OR. (ODB%BATCH(I)%AO(1)%p%type_empty.AND..NOT.ODB%BATCH(I)%AO(2)%p%type_empty)
      spherE = Input%sphericalEcoeff      
      sameAO = ODB%BATCH(I)%sameAO
      IF (nprim.GT. 9999) THEN
         CALL LSQUIT('nPrim1>9999 in SelectODTypes',lupri)
      ELSEIF (iDprim12.GT. 9999999) THEN
         CALL LSQUIT('nPrim2>9999999 in SelectODTypes',lupri)
      ELSEIF (ODB%BATCH(I)%ITYPE.GT. 99999999) THEN
         CALL LSQUIT('itype>99999999 in SelectODTypes',lupri)
      ENDIF
      !   WRITE(UniqeIdentfiers(I),'(I4,I7,4I4,5L3)')nprim,iDprim12,&
      !        & CCidentifier1,CCidentifier2,angmomIdentifier2,angmomIdentifier1,&
      !        & spher1,spher2,type_hermite_single,spherE,sameAO
!$OMP CRITICAL (ifortwrite)
      WRITE(UniqeIdentfiers(I),'(I4,I7,I8,5L3)')nprim,iDprim12,&
           & ODB%BATCH(I)%ITYPE,spher1,spher2,type_hermite_single,spherE,sameAO
!$OMP END CRITICAL (ifortwrite)
      unique = .TRUE.
      DO J=I-1,1,-1
         IF (UniqeIdentfiers(I).EQ.UniqeIdentfiers(J)) THEN
            ODTYPEIndex(I)=ODTYPEIndex(J)
            unique = .FALSE.
            EXIT
         ENDIF
      ENDDO
      IF (unique) THEN
         nODTypes = nODTypes + 1
         ODTYPEIndex(I) = nODTypes
      ENDIF
   ENDDO
ENDIF

if(nodTYPES.EQ. 0)call lsquit('should never happen. TK',lupri)

!struct have the number of overlaps of this type
!plus array of overlap indices with this type
nullify(typeOverlapindex)
allocate(typeOverlapindex(nODTypes))
DO IODTYPE=1,nODTypes
   typeOverlapindex(IODTYPE)%dim = 0
ENDDO
DO I=1,nBatches
   IODTYPE = ODTYPEIndex(I)
   TypeOverlapIndex(IODTYPE)%dim = typeOverlapindex(IODTYPE)%dim + 1
ENDDO
DO IODtype=1,nODTypes
   dim = TypeOverlapIndex(IODTYPE)%dim
   call mem_alloc(TypeOverlapindex(IODtype)%elms,dim)
   typeOverlapindex(IODTYPE)%elms = 0
   N=0
   DO I=1,nBatches
      IF(ODTYPEIndex(I).EQ.IODTYPE)THEN
         N=N+1
         TypeOverlapindex(IODtype)%elms(N)=I
      ENDIF
   ENDDO
   IF(N.NE.TypeOverlapindex(IODtype)%dim)call lsquit('error in selectodtype',lupri)
ENDDO

IF (IPRINT.GT. 0) THEN
   CALL LSHEADER(LUPRI,'Output from SelectODTypesFromODbatch')
   WRITE(LUPRI,'(1X,A,I5)') 'Number of OD-batches', nBatches
   WRITE(LUPRI,'(1X,A,I5)') 'Number of ODtypes   ', nODTypes
   IF (IPRINT.GT. 5) THEN
      WRITE(LUPRI,'(3X,A)') 'Batch TYPE      nPrim  angInd1  angInd2   s1 s2 hS AO'
      DO I=1,nBatches 
         WRITE(LUPRI,'(3X,2I5,2X,1A57)') I,ODTYPEIndex(I),UniqeIdentfiers(I)
      ENDDO
   ENDIF
   IF (IPRINT.GT. 5) THEN
      WRITE(LUPRI,'(3X,A)') 'Batch TYPE      nPrim  angInd1  angInd2   s1 s2 hS AO'
      DO IODTYPE=1,nODTypes
         DO I=1,TypeOverlapindex(IODtype)%dim
            WRITE(LUPRI,'(3X,2I5,2X,1A57)') TypeOverlapindex(IODtype)%elms(I),&
                 &IODtype,UniqeIdentfiers(TypeOverlapindex(IODtype)%elms(I))
         ENDDO
      ENDDO
   ENDIF
ENDIF

NMOM = 1
IF (LHS) THEN
   geoderivorder = INPUT%geoderOrderP
   nCartesianMomentComp = 1
   ndim5 = INPUT%ngeoderivcompP*INPUT%nmagderivcompP
   Magderiv = Input%MagderOrderP
ELSE
   geoderivorder = INPUT%geoderOrderQ
   nCartesianMomentComp = INPUT%nCartesianMomentComp
   ndim5 = INPUT%ngeoderivcompQ*nCartesianMomentComp*INPUT%nmagderivcompQ 
   IF(INPUT%OPERATOR .EQ.MulmomOperator) THEN
      NMOM         = INPUT%nMultipoleMomentComp
   ENDIF
   Magderiv = Input%MagderOrderQ
ENDIF

call Allocitem_alloc(Alloc,nODTypes,SIDE)
ngeoderivcomp = 1
nOperatorComp = 1
IF(Magderiv.EQ.1)nOperatorComp = 3
DO iODtype=1,nODTypes
   I=TypeOverlapindex(IODtype)%elms(1)
   single = ODB%BATCH(I)%AO(1)%p%type_empty.OR.ODB%BATCH(I)%AO(2)%p%type_empty
   type_hermite_single=input%hermiteEcoeff.AND.&
        &((.NOT.ODB%BATCH(I)%AO(1)%p%type_empty.AND.ODB%BATCH(I)%AO(2)%p%type_Empty).OR.&
        &(ODB%BATCH(I)%AO(1)%p%type_empty.AND..NOT.ODB%BATCH(I)%AO(2)%p%type_Empty))
   type_empty = (ODB%BATCH(I)%AO(1)%p%type_empty.AND.ODB%BATCH(I)%AO(2)%p%type_empty)
   IF(geoderivorder.GT. 0)CALL getDerivComp(ngeoderivcomp,geoderivorder,type_empty,single)
   maxijk = 0
   totOrbitals = 0
   nETUV = 0
   maxangmom=0
   minangmom=0
   np = ODB%BATCH(I)%AO(1)%p%nprimitives*ODB%BATCH(I)%AO(2)%p%nprimitives
   DO i1=1,ODB%BATCH(I)%AO(1)%p%nAngmom
      nOrb1 = ODB%BATCH(I)%AO(1)%p%nOrbitals(i1)
      DO i2=1,ODB%BATCH(I)%AO(2)%p%nAngmom
         l1   = ODB%BATCH(I)%AO(1)%p%angmom(i1)
         l2   = ODB%BATCH(I)%AO(2)%p%angmom(i2)
         l  = l1+l2+geoderivorder+CMorder+magderiv
         nTUV = ((l+1)*(l+2)*(l+3)/6) !including nEFG if any
         CALL GET_IJK(l1,l2,ijk1,ijk2,ijk,ijkcart,.FALSE.,geoderivorder,single)
         maxijk = max(maxijk,ijkcart,ijk)*nOperatorComp
         nETUV = nETUV + nTUV*max(ijkcart,ijk)*np*nOperatorComp
         maxangmom = MAX(maxangmom,l)
         minangmom = MIN(minangmom,l)
         totOrbitals = totOrbitals + ndim5*NMOM*(nOrb1*ODB%BATCH(I)%AO(2)%p%nOrbitals(i2))
      ENDDO
   ENDDO
   maxangmom = maxangmom
   minangmom = minangmom
   startAngmom = 0
   IF (type_hermite_single) startAngmom = minAngmom 
   endAngmom = maxAngmom
   nTUV = (endAngmom+1)*(endAngmom+2)*(endAngmom+3)/6
   IF (INPUT%operator .EQ. KineticOperator) nTUV = (endAngmom+3)*(endAngmom+4)*(endAngmom+5)/6

   IF (Side.EQ.'LHS') THEN
      IF (endAngmom+1.GT. 20) CALL LSQUIT('Error maxPrimAngmomLHS SelectODTypesFromODbatch',lupri)
      alloc%maxPrimAngmomLHS(endAngmom+1) = max(alloc%maxPrimAngmomLHS(endAngmom+1),np)
      alloc%maxPrimLHSA(iODtype) = np
      alloc%maxPrimTUVLHSA(iODtype) = MAX(np*nTUV,totOrbitals,np*maxijk)
      Alloc%maxContLHSA(iODtype) = ODB%BATCH(I)%AO(1)%p%maxContracted*ODB%BATCH(I)%AO(2)%p%maxContracted
      Alloc%maxnAngLHSA(iODtype) = ODB%BATCH(I)%nAngmom
      Alloc%maxangmomOrbLHSA(iODtype) = max(ODB%BATCH(I)%AO(1)%p%maxAngmom,ODB%BATCH(I)%AO(2)%p%maxAngmom)
      Alloc%maxTUVLHSA(iODtype) = nTUV
      Alloc%maxTotOrbLHSA(iODtype) = totOrbitals
      Alloc%maxijkLHSA(iODtype) = maxijk
      Alloc%maxETUVlenLHSA(iODtype) = nETUV
   ELSE
      IF (endAngmom+1.GT. 20) CALL LSQUIT('Error maxPrimAngmomRHS SelectODTypesFromODbatch',lupri)
      alloc%maxPrimAngmomRHS(endAngmom+1) = max(alloc%maxPrimAngmomRHS(endAngmom+1),np)
      alloc%maxPrimRHSA(iODtype) = np
      alloc%maxPrimTUVRHSA(iODtype) = MAX(np*nTUV,totOrbitals,np*maxijk)
      Alloc%maxContRHSA(iODtype) = ODB%BATCH(I)%AO(1)%p%maxContracted*ODB%BATCH(I)%AO(2)%p%maxContracted
      Alloc%maxnAngRHSA(iODtype) = ODB%BATCH(I)%nAngmom
      Alloc%maxangmomOrbRHSA(iODtype) = max(ODB%BATCH(I)%AO(1)%p%maxAngmom,ODB%BATCH(I)%AO(2)%p%maxAngmom)
      Alloc%maxTUVRHSA(iODtype) = nTUV
      Alloc%maxTotOrbRHSA(iODtype) = totOrbitals
      Alloc%maxijkRHSA(iODtype) = maxijk
      Alloc%maxETUVlenRHSA(iODtype) = nETUV
   ENDIF
ENDDO
call Allocitem_collect(Alloc,SIDE)
call mem_dealloc(UniqeIdentfiers)
END SUBROUTINE SelectODTypesFromODbatch

!> \brief determine the number of primitives from OD Batch
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param nPrim the number of primitives
!> \param ODB the overlap batch
!> \param INPUT the integral input specifications, contains all info about what to do
!> \param side LHS or RHS specification
!> \param lupri The default print unit
SUBROUTINE GET_NPRIM_FROM_ODBATCH(nPrim,ODB,input,side,lupri)
IMPLICIT NONE
INTEGER :: nPrim,lupri
Character*(*)       :: side
TYPE(IntegralInput),target :: Input
TYPE(ODBATCH)       :: ODB
!
integer             :: i1,start2,i2,maxbat
!Character(len=80)   :: t1,t2
integer(kind=short) :: maxgab,MAXELM
!Real(realk),pointer :: GAB(:,:)
type(lstensor),pointer :: GAB
integer             :: AtomA,atomB,batchA,batchB,Gindex,s1,s2,n1
!real(realk),parameter :: TEN=1E1_realk

IF (Side.EQ.'LHS') THEN
  GAB => INPUT%LST_GAB_LHS
  MAXELM=INPUT%PS_MAXELM_RHS
ELSEIF (Side.EQ.'RHS') THEN
  GAB => INPUT%LST_GAB_RHS
  MAXELM=INPUT%PS_MAXELM_LHS
ELSE
  WRITE(*,'(1X,2A)') 'Error in GET_NPRIM_FROM_ODBATCH. Side =',Side
  CALL LSQUIT('Programming error. Wrong Side in GET_NPRIM_FROM_ODBATCH',-1)
ENDIF
nPrim = 0
IF (.NOT.(ODB%AO(1)%p%type_empty .AND.ODB%AO(2)%p%type_empty)) THEN
   IF(INPUT%PS_SCREEN)THEN   
      AtomA = ODB%AO(1)%p%atom
      AtomB = ODB%AO(2)%p%atom
      BatchA = ODB%AO(1)%p%batch
      BatchB = ODB%AO(2)%p%batch
      Gindex = GAB%INDEX(atomA,atomB,1,1)
      maxbat = GAB%SLSAO(Gindex)%maxbat
      s1 = GAB%SLSAO(Gindex)%startLocalOrb(BatchA)-1
      s2 = GAB%SLSAO(Gindex)%startLocalOrb(BatchB+maxbat)-1
      n1 = GAB%SLSAO(Gindex)%nLocal(1)
      DO i1=1,ODB%AO(1)%p%nPrimitives
         DO i2=1,ODB%AO(2)%p%nPrimitives
            maxGab = GAB%SLSAO(Gindex)%selms(s1+i1 + (s2+i2-1)*n1)
            IF( MAXGAB .GT. INPUT%PS_THRLOG-MAXELM)THEN
               !Add if maxgab greater than threshold
               nPrim = nPrim + 1
            ENDIF
         ENDDO
      ENDDO
   ELSE
      nPrim = ODB%AO(1)%p%nPrimitives*ODB%AO(2)%p%nPrimitives
   ENDIF
ELSE
   nPrim = 1
ENDIF

END SUBROUTINE GET_NPRIM_FROM_ODBATCH

!> \brief initialise the overlap
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param PassP The overlap
!> \param Alloc Information about the maximal sizes for allocation purposes
!> \param Input contains info about the integral evaluation requested
!> \param ODbat The OD-batch information
!> \param side The side ('LHS' or 'RHS') 
!> \param maxpasses the maximum number of passes
!> \param iprint The print level (the higher the more information)
!> \param lupri The default print unit
SUBROUTINE INIT_PASS(PassP,Alloc,iAlloc,Input,ODbat,SIDE,maxpasses,IPRINT,LUPRI)
Implicit none
INTEGER             :: lupri,IPRINT,iAlloc
TYPE(Overlap)       :: PassP
TYPE(IntegralInput) :: Input
TYPE(AllocItem)     :: Alloc
TYPE(ODITEM)        :: ODbat
!TYPE(Integralitem)  :: Integral
Character*(*)       :: side
!TYPE(ODBATCH)       :: ODbat(nbatches)
!
!Character(len=80)   :: t1,t2
!Integer             :: NDERIV
Integer             :: endangmom,IELECTRON
Integer             :: maxpasses
Integer             :: lenETUV,ma,mc,mp,na,np,ng
LOGICAL             :: LHS,type_hermite_single 

SELECT CASE(SIDE)
CASE('LHS')
   LHS=.TRUE.
!ORBITAL
   np = Alloc%maxPrimLHSA(iAlloc) 
   na = Alloc%maxnAngLHSA(iAlloc) 
   mc = Alloc%maxContLHSA(iAlloc) 
   ma = Alloc%maxAngmomOrbLHSA(iAlloc) 
   ng = Input%geoDerOrderP
   IELECTRON = 1
CASE('RHS')
   LHS=.FALSE.
   np = Alloc%maxPrimRHSA(iAlloc)  
   na = Alloc%maxnAngRHSA(iAlloc)  
   mc = Alloc%maxContRHSA(iAlloc) 
   ma = Alloc%maxAngmomOrbRHSA(iAlloc) 
   ng = Input%geoDerOrderQ
   IELECTRON = 2
CASE DEFAULT
   WRITE(LUPRI,'(1X,2A)') 'Wrong case in SET_INITIALIZED_OVERLAP =',SIDE
   CALL LSQUIT('Wrong case in SET_INITIALIZED_OVERLAP',lupri)
END SELECT

CALL ALLOC_ORBITAL(PassP%Orbital1,np,maxpasses,na,ma,.FALSE.,ielectron)
CALL ALLOC_ORBITAL(PassP%Orbital2,np,maxpasses,na,ma,.FALSE.,ielectron)
call mem_ODpointer_alloc(PassP%Orb1atom,maxpasses,ielectron)
call mem_ODpointer_alloc(PassP%Orb1mol,maxpasses,ielectron)
call mem_ODpointer_alloc(PassP%Orb1batch,maxpasses,ielectron)
call mem_ODpointer_alloc(PassP%Orb2atom,maxpasses,ielectron)
call mem_ODpointer_alloc(PassP%Orb2mol,maxpasses,ielectron)
call mem_ODpointer_alloc(PassP%Orb2batch,maxpasses,ielectron)
  
IF(LHS)THEN !OVERLAP
   mp = Alloc%maxPrimLHSA(iAlloc)*maxPasses
   np = Alloc%maxPrimLHSA(iAlloc)  
   na = Alloc%maxnAngLHSA(iAlloc)  
!   lenETUV = Alloc%maxetuvlenLHSA(iAlloc) 
ELSE
   mp = Alloc%maxPrimRHSA(iAlloc)*maxPasses
   np = Alloc%maxPrimRHSA(iAlloc)  
   na = Alloc%maxnAngRHSA(iAlloc) 
!   lenETUV = Alloc%maxetuvlenRHSA(iAlloc)
ENDIF
!type_hermite_single = input%hermiteEcoeff.AND.&
!     & ((.NOT.ODbat%BATCH(1)%AO(1)%p%type_empty .AND. ODbat%BATCH(1)%AO(2)%p%type_empty) &
!     &.OR.(ODbat%BATCH(1)%AO(1)%p%type_empty .AND. .NOT.ODbat%BATCH(1)%AO(2)%p%type_empty))
CALL ALLOC_OVERLAP(PassP,mp,maxpasses,na,.FALSE.,np,ng,1,ielectron)

 END SUBROUTINE INIT_PASS

!> \brief print a general 3 dim tensor
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param tensor the tensor to be printet
!> \param label printing label
!> \param dim1 the size of dimension 1
!> \param dim2 the size of dimension 2
!> \param dim3 the size of dimension 3
!> \param lupri The default print unit
!> \param label1 character string for dimension 1
!> \param label2 character string for dimension 2
!> \param label3 character string for dimension 3
!> \param option determines which dimension it should print in the innermost loop
SUBROUTINE PrintTensor(Tensor,label,dim1,dim2,dim3,Lupri,Label1,Label2,Label3,option)
integer          :: dim1,dim2,dim3,lupri,option
real(realk),dimension(dim1,dim2,dim3) :: Tensor
character(len=20):: Label
character(len=6) :: Label1,Label2,label3

  WRITE(LUPRI,'(5X,A)') '***************************************************'
  WRITE(LUPRI,'(5X,A,A20)') '***               ',LABEL
  WRITE(LUPRI,'(5X,A)') '***************************************************'
IF(option .EQ. 1)THEN
   WRITE(LUPRI,'(1X,A6,1X,A6,12X,A6 )') Label2,label3,label1
   DO j=1,dim2
      DO k=1,dim3
         WRITE(LUPRI,'(1X,I4,2X,I4,5ES18.9/,(11X,5ES18.9))') j,k,(Tensor(i,j,k),i=1,dim1)
      ENDDO
   ENDDO
ELSEIF(option .EQ. 2)THEN
   WRITE(LUPRI,'(1X,A6,1X,A6,12X,A6 )') Label1,label3,label2
   DO i=1,dim1
      DO k=1,dim3
         WRITE(LUPRI,'(1X,I4,2X,I4,5ES18.9/,(11X,5ES18.9))') i,k,(Tensor(i,j,k),j=1,dim2)
      ENDDO
   ENDDO
ELSEIF(option .EQ. 3)THEN
   WRITE(LUPRI,'(1X,A6,1X,A6,12X,A6 )') Label1,label2,label3
   DO i=1,dim1
      DO j=1,dim2
         WRITE(LUPRI,'(1X,I4,2X,I4,5ES18.9/,(11X,5ES18.9))') i,j,(Tensor(i,j,k),k=1,dim3)
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE PrintTensor

!> \brief determine the number of cartesian or spherical angular components
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param l1 the angular momentum for orbital 1
!> \param l2 the angular momentum for orbital 2
!> \param ijk1 the number of angular components for orbital 1
!> \param ijk2 the number of angular components for orbital 2
!> \param ijk the total number of angular components
!> \param ijkcart the total number of cartesian angular components
!> \param spherical describes if the spherical components should be determined
!> \param geoderivorder the order of geometrical derivative
!> \param single if only orbital 1 or 2 is used  
SUBROUTINE GET_IJK(l1,l2,ijk1,ijk2,ijk,ijkcart,spherical,geoderivorder,single)
implicit none
Integer :: l1,l2,ijk1,ijk2,ijk,ijkcart,geoderivorder
Logical :: spherical,single
!
Integer :: iDer,ijk1der,ijk2der,der1,der2
!
!Consistency check
IF ((geoderivorder.GT. 0).AND.spherical) CALL LSQUIT('Error in get_ijk. geoderivorder>0 and spherical',-1)

! Regular case
IF (geoderivorder.EQ. 0) THEN
  ijk1 = (l1 + 1)*(l1 + 2)/2
  ijk2 = (l2 + 1)*(l2 + 2)/2
  ijkcart = ijk1*ijk2
  IF (spherical) THEN
   ijk1 = 2*l1 + 1
   ijk2 = 2*l2 + 1
  ENDIF
  ijk = ijk1*ijk2
! Derivative case. Note here that we assume Hermite Guassians following PCCP,2007,9,4771-4779
ELSEIF (geoderivorder.GT. 0) THEN
  ijk1 = (l1 + 1)*(l1 + 2)/2
  ijk2 = (l2 + 1)*(l2 + 2)/2
  IF (single) THEN
    IF (l1.GE.l2) THEN
      ijk1der = (l1 + 1 + geoderivorder)*(l1 + 2 + geoderivorder)/2
      ijkcart = ijk1der*ijk2
    ELSE
      ijk2der = (l2 + 1 + geoderivorder)*(l2 + 2 + geoderivorder)/2
      ijkcart = ijk1*ijk2der
    ENDIF
  ELSE
    ijk = 0
    ijkcart = 0
    DO iDer=0,geoderivorder
      der1 = iDer
      der2 = geoderivorder - iDer
      ijk1der = (l1 + 1 + der1)*(l1 + 2 + der1)/2
      ijk2der = (l2 + 1 + der2)*(l2 + 2 + der2)/2
      ijkcart = ijkcart + ijk1der*ijk2der
    ENDDO
  ENDIF
  ijk = ijkcart
ELSE
  CALL LSQUIT('Error in GET_IJK. geoderivorder<0',-1)
ENDIF
END SUBROUTINE GET_IJK

 !> \brief select overlap distribution pass types from overlap structure
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param ODpassesIndex for a given batch this gives the passtype 
!> \param P the overlap 
!> \param nBatches the number of batches
!> \param nPassTypes the number of passtypes 
!> \param iprint The print level (the higher the more information)
!> \param lupri The default print unit
SUBROUTINE SelectODPassTypes2(ODpassesIndex,PassTypeOverlapindex,Q,nBatches,nPassTypes,alloc,input,IPRINT,LUPRI)
use memory_handling
implicit none
Integer       :: nBatches,nPassTypes,IPRINT,LUPRI
Integer       :: ODpassesIndex(nBatches)
TYPE(IntegralInput) :: Input
Type(Overlap) :: Q(nBatches)
type(allocitem) :: alloc
type(integerpointer),pointer :: PassTypeOverlapindex(:)

Integer       :: I,J,nPrim,numAng1,numAng2
Character(len=69),pointer  :: UniqeIdentfiers(:)
Integer       :: angmomIdentifier1,angmomIdentifier2,CCidentifier1,CCidentifier2
Logical       :: spher1,spher2,type_hermite_single,unique,sameAO,spherE,succes
!
integer :: maxnp,maxijk,i1,start2,i2,l1,l2,l,ijk1,ijk2,ijk,ijkcart
integer :: ipass,dim,N,nETUV,nTUV,totorb,ip,iprim12,CMorder,nOperatorComp
CMorder = input%CMorder
call mem_alloc(UniqeIdentfiers,nBatches)

nPassTypes = 0
DO I=1,nBatches
  nPrim   = Q(I)%nPrimitives
  !ToDo select according to the same iprim1 and iprim2 set-ups 
  !(used for primitive screening)
  numAng1 = Q(I)%orbital1%nAngmom
  numAng2 = Q(I)%orbital2%nAngmom
  angmomIdentifier1 = 0
  angmomIdentifier2 = 0
  DO J=1,numAng1
     angmomIdentifier1 = angmomIdentifier1 + 2**Q(I)%orbital1%angmom(J)
  ENDDO
  DO J=1,numAng2
     angmomIdentifier2 = angmomIdentifier2 + 2**Q(I)%orbital2%angmom(J)
  ENDDO
  iprim12 = 0
  DO IP=1,nPrim
     iprim12 = iprim12 + Q(I)%iprim1(iP) + Q(I)%iprim2(iP)*nprim
  ENDDO
  
  CCidentifier1 = Q(I)%orbital1%CCidentifier
  CCidentifier2 = Q(I)%orbital2%CCidentifier
  spher1        = Q(I)%orbital1%spherical
  spher2        = Q(I)%orbital2%spherical
  type_hermite_single = Q(I)%type_hermite_single
  spherE        = Q(I)%sphericalEcoeff
  sameAO        = Q(I)%sameAO
  
  IF (nPrim.GT. 999999999) THEN
     CALL LSQUIT('nPrim>999999999 in SelectODPassTypes',lupri)
  ELSE IF (angmomIdentifier1.GT. 999999999) THEN
     CALL LSQUIT('angmomIdentifier1>999999999 in SelectODPassTypes',lupri)
  ELSE IF (angmomIdentifier2.GT. 999999999) THEN
     CALL LSQUIT('angmomIdentifier2>999999999 in SelectODPassTypes',lupri)
  ELSE IF (CCidentifier1.GT. 999999999) THEN
     CALL LSQUIT('CCidentifier1>999999999 in SelectODPassTypes',lupri)
  ELSE IF (CCidentifier2.GT. 999999999) THEN
     CALL LSQUIT('CCidentifier2>999999999 in SelectODPassTypes',lupri)
  ENDIF
  !$OMP CRITICAL (ifortwrite)
  WRITE(UniqeIdentfiers(I),'(6I9,5L3)') nPrim,iprim12,angmomIdentifier2,angmomIdentifier1,&
       & CCidentifier1,CCidentifier2,spher1,spher2,type_hermite_single,spherE,sameAO
  !$OMP END CRITICAL (ifortwrite)
  unique = .TRUE.
  DO J=I-1,1,-1
     IF (UniqeIdentfiers(I).EQ.UniqeIdentfiers(J)) THEN
        ODpassesIndex(I)=ODpassesIndex(J)
        unique = .FALSE.
        EXIT
     ENDIF
  ENDDO
  IF (unique) THEN
     nPassTypes = nPassTypes + 1
     ODpassesIndex(I) = nPassTypes
  ENDIF
ENDDO

nOperatorComp = 1
IF(Q(1)%magderiv.EQ.1)nOperatorComp = 3
call Allocitem_alloc(Alloc,nPassTypes,'RHS')
DO ipass=1,npassTypes
 DO I=1,nBatches
  nPrim   = Q(I)%nPrimitives
  IF(nPrim.GT. 0)THEN
   IF(ipass .EQ. ODpassesIndex(I))THEN
      maxnp = Q(I)%orbital1%nPrimitives*Q(I)%orbital2%nPrimitives
      maxijk = 0
      nETUV = 0
      DO i1=1,Q(I)%orbital1%nAngmom
!         start2 = 1
!         IF (Q(I)%sameAO) start2 = i1
!         DO i2=start2,Q(I)%orbital2%nAngmom
         DO i2=1,Q(I)%orbital2%nAngmom
            l1   = Q(I)%orbital1%angmom(i1)
            l2   = Q(I)%orbital2%angmom(i2)
            l = l1+l2+Q(I)%endGeoOrder+CMorder+Q(I)%magderiv
            nTUV = (l+1)*(l+2)*(l+3)/6 !including nEFG if any
            CALL GET_IJK(l1,l2,ijk1,ijk2,ijk,ijkcart,Q(I)%sphericalEcoeff,Q(I)%endGeoOrder,Q(I)%single)
            maxijk = max(maxijk,ijk)*nOperatorComp
            nETUV = nETUV + nTUV*ijkcart*maxnp*nOperatorComp
         ENDDO
      ENDDO
      alloc%maxPrimRHSA(ipass) = maxnp
      totorb = Q(I)%orbital1%totOrbitals*Q(I)%orbital2%totOrbitals
      alloc%maxPrimTUVRHSA(iPass) = MAX(maxnp*Q(I)%nTUV,totOrb,maxnp*maxijk)
      Alloc%maxContRHSA(ipass) = Q(I)%Orbital1%maxContracted*Q(I)%Orbital2%maxContracted
      Alloc%maxnAngRHSA(ipass) = Q(I)%nAngmom
      Alloc%maxangmomOrbRHSA(ipass) = max(Q(I)%Orbital1%maxAngmom,Q(I)%Orbital2%maxAngmom)
      Alloc%maxTUVRHSA(ipass) = Q(I)%nTUV
      Alloc%maxTotOrbRHSA(ipass) = Q(I)%totOrbitals(Q(I)%endGeoOrder+1)
      Alloc%maxijkRHSA(ipass) = maxijk      
      Alloc%maxETUVlenRHSA(ipass) = nETUV
      EXIT
   ENDIF
  ENDIF
 ENDDO
ENDDO
call Allocitem_collect(Alloc,'RHS')

!struct have the number of overlaps of this type
!plus array of overlap indices with this type
nullify(PasstypeOverlapindex)
allocate(PasstypeOverlapindex(nPassTypes))
DO Ipass=1,nPassTypes
   PasstypeOverlapindex(Ipass)%dim = 0
ENDDO
DO I=1,nBatches
   IF(Q(I)%nPrimitives.GT. 0)THEN
      PassTypeOverlapIndex(ODPASSESIndex(I))%dim = PasstypeOverlapindex(ODPASSESIndex(I))%dim + 1
   ENDIF
ENDDO
DO Ipass=1,nPassTypes
   dim = PassTypeOverlapIndex(Ipass)%dim
   call mem_alloc(PassTypeOverlapindex(Ipass)%elms,dim)
   N=0
   DO I=1,nBatches
      IF((Q(I)%nPrimitives.GT. 0).AND.(ODpassesIndex(I).EQ.Ipass))THEN
         N=N+1
         PassTypeOverlapindex(Ipass)%elms(N)=I
      ENDIF
   ENDDO
   IF(N.NE.PassTypeOverlapindex(Ipass)%dim)call lsquit('error in selectodtype',lupri)
ENDDO

IF (IPRINT.GT. 0) THEN
  CALL LSHEADER(LUPRI,'Output from SelectODPassTypes')
  WRITE(LUPRI,'(1X,A,I5)') 'Number of OD-batches', nBatches
  WRITE(LUPRI,'(1X,A,I5)') 'Number of pass-types', nPassTypes
  IF (IPRINT.GT. 5) THEN
    WRITE(LUPRI,'(3X,A)') 'Batch Pass      nPrim  angInd1  angInd2   CCind1    CCind2 s1 s2 hS AO'
    DO I=1,nBatches 
      WRITE(LUPRI,'(3X,2I5,2X,1A57)') I,ODpassesIndex(I),UniqeIdentfiers(I)
    ENDDO
  ENDIF
ENDIF

call mem_dealloc(UniqeIdentfiers)

END SUBROUTINE SelectODPassTypes2

!=====================================================
! Allocitem routines
!=====================================================
subroutine Allocitem_free(Alloc,Side)
type(allocitem) :: Alloc
Character*(*)      :: Side

IF (SIDE.EQ.'Both'.OR.SIDE.EQ.'LHS') THEN
   call mem_dealloc(alloc%maxPrimLHSA)
   call mem_dealloc(alloc%maxPrimTUVLHSA)
   call mem_dealloc(alloc%maxnAngLHSA)
   call mem_dealloc(alloc%maxContLHSA)
   call mem_dealloc(alloc%maxAngmomOrbLHSA)
   call mem_dealloc(alloc%maxTUVLHSA)
!   call mem_dealloc(alloc%maxPrimTUVijkLHSA)
   call mem_dealloc(alloc%maxTotOrbLHSA)
   call mem_dealloc(alloc%maxijkLHSA)
   call mem_dealloc(alloc%maxETUVlenLHSA)
ENDIF
IF (SIDE.EQ.'Both'.OR.SIDE.EQ.'RHS') THEN
   call mem_dealloc(alloc%maxPrimRHSA)
   call mem_dealloc(alloc%maxPrimTUVRHSA)
   call mem_dealloc(alloc%maxnAngRHSA)
   call mem_dealloc(alloc%maxContRHSA)
   call mem_dealloc(alloc%maxAngmomOrbRHSA)
   call mem_dealloc(alloc%maxTUVRHSA)
!   call mem_dealloc(alloc%maxPrimTUVijkRHSA)
   call mem_dealloc(alloc%maxTotOrbRHSA)
   call mem_dealloc(alloc%maxijkRHSA)
   call mem_dealloc(alloc%maxETUVlenRHSA)
ENDIF

end subroutine Allocitem_free

subroutine Allocitem_alloc(Alloc,nsize,SIDE)
type(allocitem) :: Alloc
integer :: nsize
Character*(*)      :: Side

IF (SIDE.EQ.'Both'.OR.SIDE.EQ.'LHS') THEN
   call mem_alloc(alloc%maxPrimLHSA,nsize)
   call mem_alloc(alloc%maxPrimTUVLHSA,nsize)
   call mem_alloc(alloc%maxnAngLHSA,nsize)
   call mem_alloc(alloc%maxContLHSA,nsize)
   call mem_alloc(alloc%maxAngmomOrbLHSA,nsize)
   call mem_alloc(alloc%maxTUVLHSA,nsize)
!   call mem_alloc(alloc%maxPrimTUVijkLHSA,nsize)
   call mem_alloc(alloc%maxTotOrbLHSA,nsize)
   call mem_alloc(alloc%maxijkLHSA,nsize)
   call mem_alloc(alloc%maxETUVlenLHSA,nsize)
ENDIF
IF (SIDE.EQ.'Both'.OR.SIDE.EQ.'RHS') THEN
   call mem_alloc(alloc%maxPrimRHSA,nsize)
   call mem_alloc(alloc%maxPrimTUVRHSA,nsize)
   call mem_alloc(alloc%maxnAngRHSA,nsize)
   call mem_alloc(alloc%maxContRHSA,nsize)
   call mem_alloc(alloc%maxAngmomOrbRHSA,nsize)
   call mem_alloc(alloc%maxTUVRHSA,nsize)
!   call mem_alloc(alloc%maxPrimTUVijkRHSA,nsize)
   call mem_alloc(alloc%maxTotOrbRHSA,nsize)
   call mem_alloc(alloc%maxijkRHSA,nsize)
   call mem_alloc(alloc%maxETUVlenRHSA,nsize)
ENDIF

end subroutine Allocitem_alloc

!> \brief initialise the allocitem 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param Alloc the allocitem, which contain info about maximum values required for allocation
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
!> \param SIDE either both if both side should be initialised or LHS,RHS if only one side
SUBROUTINE Allocitem_collect(Alloc,SIDE)
implicit none
TYPE(Allocitem)    :: Alloc
Character*(*)      :: Side
!
integer :: nsize,I 

IF (SIDE.EQ.'Both'.OR.SIDE.EQ.'LHS') THEN
!  LHS
   nsize = size(Alloc%maxPrimLHSA)
   call Allocitem_zero(Alloc,'LHS')
   DO I = 1,nsize
      alloc%maxPrimLHS = MAX(alloc%maxPrimLHS,alloc%maxPrimLHSA(I))
      alloc%maxPrimTUVLHS = MAX(alloc%maxPrimTUVLHS,alloc%maxPrimTUVLHSA(I))
      alloc%maxnAngLHS = MAX(alloc%maxnAngLHS,alloc%maxnAngLHSA(I))
      alloc%maxContLHS = MAX(alloc%maxContLHS,alloc%maxContLHSA(I))
      alloc%maxAngmomOrbLHS = MAX(alloc%maxAngmomOrbLHS,alloc%maxAngmomOrbLHSA(I))
      alloc%maxTUVLHS = MAX(alloc%maxTUVLHS,alloc%maxTUVLHSA(I))
!      alloc%maxPrimTUVijkLHS = MAX(alloc%maxPrimTUVijkLHS,alloc%maxPrimTUVijkLHSA(I))
      alloc%maxTotOrbLHS = MAX(alloc%maxTotOrbLHS,alloc%maxTotOrbLHSA(I))
      alloc%maxijkLHS = MAX(alloc%maxijkLHS,alloc%maxijkLHSA(I))
      alloc%maxETUVlenLHS = MAX(alloc%maxETUVlenLHS,alloc%maxETUVlenLHSA(I))
   ENDDO
ENDIF
IF (SIDE.EQ.'Both'.OR.SIDE.EQ.'RHS') THEN
!  LHS
   nsize = size(Alloc%maxPrimRHSA)
   call Allocitem_zero(Alloc,'RHS')
   DO I = 1,nsize
      alloc%maxPrimRHS = MAX(alloc%maxPrimRHS,alloc%maxPrimRHSA(I))
      alloc%maxPrimTUVRHS = MAX(alloc%maxPrimTUVRHS,alloc%maxPrimTUVRHSA(I))
      alloc%maxnAngRHS = MAX(alloc%maxnAngRHS,alloc%maxnAngRHSA(I))
      alloc%maxContRHS = MAX(alloc%maxContRHS,alloc%maxContRHSA(I))
      alloc%maxAngmomOrbRHS = MAX(alloc%maxAngmomOrbRHS,alloc%maxAngmomOrbRHSA(I))
      alloc%maxTUVRHS = MAX(alloc%maxTUVRHS,alloc%maxTUVRHSA(I))
!      alloc%maxPrimTUVijkRHS = MAX(alloc%maxPrimTUVijkRHS,alloc%maxPrimTUVijkRHSA(I))
      alloc%maxTotOrbRHS = MAX(alloc%maxTotOrbRHS,alloc%maxTotOrbRHSA(I))
      alloc%maxijkRHS = MAX(alloc%maxijkRHS,alloc%maxijkRHSA(I))
      alloc%maxETUVlenRHS = MAX(alloc%maxETUVlenRHS,alloc%maxETUVlenRHSA(I))
   ENDDO
ENDIF

END SUBROUTINE ALLOCITEM_COLLECT

!> \brief initialise the allocitem 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param Alloc the allocitem, which contain info about maximum values required for allocation
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
!> \param SIDE either both if both side should be initialised or LHS,RHS if only one side
SUBROUTINE Allocitem_init(Alloc,SIDE)
implicit none
TYPE(Allocitem)    :: Alloc
Integer            :: LUPRI,IPRINT
Character*(*)      :: Side
IF (SIDE.EQ.'Both'.OR.SIDE.EQ.'LHS') THEN
   call Allocitem_zero(Alloc,SIDE)
   nullify(alloc%maxPrimLHSA)
   nullify(alloc%maxPrimTUVLHSA)
   nullify(alloc%maxnAngLHSA)
   nullify(alloc%maxContLHSA)
   nullify(alloc%maxAngmomOrbLHSA)
   nullify(alloc%maxTUVLHSA)
   nullify(alloc%maxTotOrbLHSA)
   nullify(alloc%maxijkLHSA)
   nullify(alloc%maxETUVlenLHSA)
ENDIF
IF (SIDE.EQ.'Both'.OR.SIDE.EQ.'RHS') THEN
   call Allocitem_zero(Alloc,SIDE)
   nullify(alloc%maxPrimRHSA)
   nullify(alloc%maxPrimTUVRHSA)
   nullify(alloc%maxnAngRHSA)
   nullify(alloc%maxContRHSA)
   nullify(alloc%maxAngmomOrbRHSA)
   nullify(alloc%maxTUVRHSA)
   nullify(alloc%maxTotOrbRHSA)
   nullify(alloc%maxijkRHSA)
   nullify(alloc%maxETUVlenRHSA)
ENDIF
END SUBROUTINE Allocitem_init

!> \brief initialise the allocitem 
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param Alloc the allocitem, which contain info about maximum values required for allocation
!> \param LUPRI the logical unit number for the output file
!> \param IPRINT the printlevel, determining how much output should be generated
!> \param SIDE either both if both side should be initialised or LHS,RHS if only one side
SUBROUTINE Allocitem_zero(Alloc,SIDE)
implicit none
TYPE(Allocitem)    :: Alloc
Integer            :: LUPRI,IPRINT
Character*(*)      :: Side
IF (SIDE.EQ.'Both'.OR.SIDE.EQ.'LHS') THEN
!  LHS
   Alloc%maxPrimAngmomLHS = 0
   Alloc%maxPrimLHS    = 0
   Alloc%maxPrimTUVLHS = 0
!   Alloc%maxPrimLHSpass    = 0
!   Alloc%maxPrimTUVLHSpass = 0
   Alloc%maxnAngLHS     = 0
   Alloc%maxContLHS    = 0
   Alloc%maxAngmomLHS = 0
   Alloc%maxAngmomOrbLHS = 0
   Alloc%maxTUVLHS = 0
!   Alloc%maxPrimTUVijkLHS = 0
   Alloc%maxTotOrbLHS = 0
   Alloc%maxijkLHS = 0
   Alloc%maxETUVlenLHS = 0
ENDIF
IF (SIDE.EQ.'Both'.OR.SIDE.EQ.'RHS') THEN
   Alloc%maxPrimAngmomRHS = 0
   Alloc%maxPrimRHS    = 0 
   Alloc%maxPrimTUVRHS = 0 
!   Alloc%maxPrimRHSpass    = 0 
!   Alloc%maxPrimTUVRHSpass = 0 
   Alloc%maxnAngRHS     = 0
   Alloc%maxContRHS    = 0
   Alloc%maxAngmomRHS = 0
   Alloc%maxAngmomOrbRHS = 0
   Alloc%maxTUVRHS = 0
!   Alloc%maxPrimTUVijkRHS = 0
   Alloc%maxTotOrbRHS = 0
   Alloc%maxijkRHS = 0
   Alloc%maxETUVlenRHS = 0
ENDIF

END SUBROUTINE Allocitem_zero

subroutine Allocitem_print(Alloc,SIDE,lupri)
type(allocitem) :: Alloc
integer :: lupri
Character*(*)      :: Side

IF (SIDE.EQ.'Both'.OR.SIDE.EQ.'LHS') THEN
   IF(associated(alloc%maxPrimLHSA))THEN
      WRITE(lupri,*)'size(alloc%maxPrimLHSA)',size(alloc%maxPrimLHSA)
      DO I=1,size(alloc%maxPrimLHSA)
         WRITE(lupri,*)'alloc%maxPrimLHSA(',I,')     ',alloc%maxPrimLHSA(I)
         WRITE(lupri,*)'alloc%maxPrimTUVLHSA(',I,')  ',alloc%maxPrimTUVLHSA(I)
         WRITE(lupri,*)'alloc%maxnAngLHSA(',I,')     ',alloc%maxnAngLHSA(I)
         WRITE(lupri,*)'alloc%maxContLHSA(',I,')     ',alloc%maxContLHSA(I)
         WRITE(lupri,*)'alloc%maxAngmomOrbLHSA(',I,')',alloc%maxAngmomOrbLHSA(I)
         WRITE(lupri,*)'alloc%maxTUVLHSA(',I,')      ',alloc%maxTUVLHSA(I)
         WRITE(lupri,*)'alloc%maxTotOrbLHSA(',I,')   ',alloc%maxTotOrbLHSA(I)
         WRITE(lupri,*)'alloc%maxijkLHSA(',I,')      ',alloc%maxijkLHSA(I)
         WRITE(lupri,*)'alloc%maxETUVlenLHSA(',I,')  ',alloc%maxETUVlenLHSA(I)
      ENDDO
   ENDIF
   write(lupri,*)'alloc%maxPrimLHS',alloc%maxPrimLHS
   write(lupri,*)'alloc%maxPrimTUVLHS',alloc%maxPrimTUVLHS
   write(lupri,*)'alloc%maxnAngLHS',alloc%maxnAngLHS
   write(lupri,*)'alloc%maxContLHS',alloc%maxContLHS
   write(lupri,*)'alloc%maxAngmomOrbLHS',alloc%maxAngmomOrbLHS
   write(lupri,*)'alloc%maxTUVLHS',alloc%maxTUVLHS
   write(lupri,*)'alloc%maxTotOrbLHS',alloc%maxTotOrbLHS
   write(lupri,*)'alloc%maxijkLHS',alloc%maxijkLHS
   write(lupri,*)'alloc%maxETUVlenLHS',alloc%maxETUVlenLHS
ENDIF
IF (SIDE.EQ.'Both'.OR.SIDE.EQ.'RHS') THEN
   IF(associated(alloc%maxPrimRHSA))THEN
      WRITE(lupri,*)'size(alloc%maxPrimRHSA)',size(alloc%maxPrimRHSA)
      DO I=1,size(alloc%maxPrimRHSA)
         WRITE(lupri,*)'alloc%maxPrimRHSA(',I,')     ',alloc%maxPrimRHSA(I)
         WRITE(lupri,*)'alloc%maxPrimTUVRHSA(',I,')  ',alloc%maxPrimTUVRHSA(I)
         WRITE(lupri,*)'alloc%maxnAngRHSA(',I,')     ',alloc%maxnAngRHSA(I)
         WRITE(lupri,*)'alloc%maxContRHSA(',I,')     ',alloc%maxContRHSA(I)
         WRITE(lupri,*)'alloc%maxAngmomOrbRHSA(',I,')',alloc%maxAngmomOrbRHSA(I)
         WRITE(lupri,*)'alloc%maxTUVRHSA(',I,')      ',alloc%maxTUVRHSA(I)
         WRITE(lupri,*)'alloc%maxTotOrbRHSA(',I,')   ',alloc%maxTotOrbRHSA(I)
         WRITE(lupri,*)'alloc%maxijkRHSA(',I,')      ',alloc%maxijkRHSA(I)
         WRITE(lupri,*)'alloc%maxETUVlenRHSA(',I,')  ',alloc%maxETUVlenRHSA(I)
      ENDDO
   ENDIF
   write(lupri,*)'alloc%maxPrimRHS',alloc%maxPrimRHS
   write(lupri,*)'alloc%maxPrimTUVRHS',alloc%maxPrimTUVRHS
   write(lupri,*)'alloc%maxnAngRHS',alloc%maxnAngRHS
   write(lupri,*)'alloc%maxContRHS',alloc%maxContRHS
   write(lupri,*)'alloc%maxAngmomOrbRHS',alloc%maxAngmomOrbRHS
   write(lupri,*)'alloc%maxTUVRHS',alloc%maxTUVRHS
   write(lupri,*)'alloc%maxTotOrbRHS',alloc%maxTotOrbRHS
   write(lupri,*)'alloc%maxijkRHS',alloc%maxijkRHS
   write(lupri,*)'alloc%maxETUVlenRHS',alloc%maxETUVlenRHS
ENDIF

end subroutine Allocitem_print

!> \brief calculate maximum values and store them in the alloc item
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param Alloc the allocitem, which contain info about maximum values required for allocation
!> \param input the integral input specifications, contains all info about what to do
!> \param ODbat the ODbatch
!> \param side LHS or RHS character 
!> \param IPRINT the printlevel, determining how much output should be generated
!> \param LUPRI the logical unit number for the output file
SUBROUTINE SET_ALLOC(Alloc,Input,ODbat,SIDE,IPRINT,LUPRI)
Implicit none
INTEGER             :: nbatches,lupri,IPRINT,maxpasses
TYPE(IntegralInput) :: Input
TYPE(AllocItem)     :: Alloc
TYPE(ODITEM)        :: ODbat
Character*(*)       :: side
!Integer,optional    :: maxpassfortype(nPasstype)
!TYPE(ODBATCH)       :: ODbat(nbatches)
!
Integer             :: i1,i2,l1,l2,ijk1,ijk2,ijk,maxijk,I,nprim
Integer             :: angMom,maxangmom,totOrbitals,maxtotorbitals,NGEODERIVCOMP,NMOM
Integer             :: maxprim,endangmom,maxtuv,maxcont,maxnangmom,nOperatorComp
Integer             :: MAXPRIMTUV,MAXTOTORB,MAXPRIMIJK,MAXnTUV,l,nTUV
Integer             :: MAXPRIMTUVIJK,maxA,maxB,start,maxijk2,maxlenETUV
Integer             :: nDer,iDer,ijk1der,ijk2der,ijkcart,magderiv,NMAGDERIVCOMP
LOGICAL             :: LHS,hermiteSingle,single,type_empty

NMOM = 1
SELECT CASE(SIDE)
CASE('LHS')
   LHS=.TRUE.
   NGEODERIVCOMP = INPUT%NgeoDERIVcomp
   magderiv = INPUT%MagDerOrderP
   NMAGDERIVCOMP = INPUT%NmagDERIVcompP
CASE('RHS')
   LHS=.FALSE.
   NGEODERIVCOMP = INPUT%NgeoDERIVcomp
   IF(INPUT%OPERATOR.EQ.MulmomOperator) THEN
      NMOM = INPUT%nMultipoleMomentComp
   ENDIF
   IF(INPUT%OPERATOR.EQ.CarmomOperator) THEN
      NMOM = INPUT%nCartesianMomentComp
   ENDIF
   magderiv = INPUT%MagDerOrderQ
   NMAGDERIVCOMP = INPUT%NmagDERIVcompQ
CASE DEFAULT
   WRITE(LUPRI,'(1X,2A)') 'Wrong case in SET_ALLOC =',SIDE
   CALL LSQUIT('Wrong case in SET_ALLOC',lupri)
END SELECT

nDer = 0
IF (INPUT%DO_GRADIENT) nDer = INPUT%geoderivorder

maxPrim = 0
maxAngmom = 0
maxnAngmom = 0
maxCont = 0
maxijk = 0
MAXPRIMTUV = 0
MAXPRIMIJK = 0
MAXnTUV = 0
MAXPRIMTUVIJK = 0
maxA = 0
maxB = 0
maxTotorb=0
maxijk2=0
maxlenETUV=0
ngeoderivcomp = 1
nOperatorComp=1
IF(magderiv.EQ.1)nOperatorComp=3
DO I=1,ODbat%nbatches
   CALL GET_NPRIM_FROM_ODBATCH(nPrim,ODBat%BATCH(I),INPUT,SIDE,lupri)
   !In case of fragmentation and other things the number of nprim on the 
   !AO may be larger than nprim on the OD and we do use this maxprim to 
   !alloc the AO as well as the OD. 
   maxPrim = MAX(maxPrim,nPrim,ODBat%BATCH(I)%AO(1)%p%nPrimitives,ODBat%BATCH(I)%AO(2)%p%nPrimitives)
   maxnAngmom = MAX(maxnAngmom,ODbat%BATCH(I)%nAngmom)
   maxCont = MAX(maxCont,ODbat%BATCH(I)%maxContracted)
   maxA = MAX(maxA,ODbat%BATCH(I)%AO(1)%p%maxAngmom)
   maxB = MAX(maxB,ODbat%BATCH(I)%AO(2)%p%maxAngmom)
   single = (ODbat%BATCH(I)%AO(1)%p%TYPE_Empty.OR.ODbat%BATCH(I)%AO(2)%p%TYPE_Empty).AND..NOT.&
        &            (ODbat%BATCH(I)%AO(1)%p%TYPE_Empty.AND.ODbat%BATCH(I)%AO(2)%p%TYPE_Empty)
   IF(Input%GeoderOrderP.GT. 0)then
      type_empty = (ODBat%BATCH(I)%AO(1)%p%type_empty.AND.ODBat%BATCH(I)%AO(2)%p%type_empty)
      CALL getDerivComp(ngeoderivcomp,Input%GeoderOrderP,type_empty,single)
   ENDIF
   totOrbitals = 0
   maxijk = 0
   angMom = 0
   maxlenETUV = 0
   DO I1=1,ODbat%BATCH(I)%AO(1)%p%nAngmom
      start=1
      IF(ODbat%BATCH(I)%sameAO)start=I1
      DO I2=start,ODbat%BATCH(I)%AO(2)%p%nAngmom
         l1 = ODbat%BATCH(I)%AO(1)%p%angmom(I1)
         l2 = ODbat%BATCH(I)%AO(2)%p%angmom(I2)
         l  = l1+l2+nDer+input%CMorder+magderiv
         nTUV = ((l+1)*(l+2)*(l+3)/6) !including nEFG if any
         CALL GET_IJK(l1,l2,ijk1,ijk2,ijk,ijkcart,.FALSE.,nDer,single)
         maxlenETUV = maxlenETUV + nTUV*ijkcart*nPrim*nOperatorComp
         maxijk = max(maxijk,ijk)*nOperatorComp
         angMom = max(angMom,l1+l2)
         totOrbitals = totOrbitals + &
& ODbat%BATCH(I)%AO(1)%p%nOrbitals(i1)*ODbat%BATCH(I)%AO(2)%p%nOrbitals(i2)*NGEODERIVCOMP*NMOM*NMAGDERIVCOMP
      ENDDO
   ENDDO
   MAXPRIMIJK = MAX(MAXPRIMIJK,maxijk*nprim)
   maxTotorb = MAX(maxTotorb,totOrbitals)
   maxijk2 = MAX(maxijk2,maxijk)
   endAngmom  = angMom
   IF(INPUT%operator.EQ.KineticOperator.AND.LHS) endAngmom  = angMom + 2
   endAngmom = endAngmom + nDer + input%CMorder + magderiv
   maxTUV = (endAngmom+1)*(endAngmom+2)*(endAngmom+3)/6
   MAXPRIMTUV = MAX(MAXPRIMTUV,maxTUV*nPrim)
!   MAXnTUV=MAX(MAXnTUV,maxTUV)
   MAXPRIMTUVIJK = MAX(MAXPRIMTUVIJK,maxlenETUV)!maxijk*nPrim*maxTUV*ODbat%BATCH(I)%nAngmom) 
  IF(LHS)THEN
     IF (endAngmom+1.GT. 20) CALL LSQUIT('Error maxPrimAngmomLHS SET_ALLOC',lupri)
     Alloc%maxPrimAngmomLHS(angMom+1) = max(Alloc%maxPrimAngmomLHS(angMom+1),nPrim)
  ELSE 
     IF (endAngmom+1.GT. 20) CALL LSQUIT('Error maxPrimAngmomRHS SET_ALLOC',lupri)
     Alloc%maxPrimAngmomRHS(angMom+1) = max(Alloc%maxPrimAngmomRHS(angMom+1),nPrim)
  ENDIF
  maxAngmom = max(maxAngmom,angMom)
ENDDO

IF(LHS)THEN
!   Alloc%maxJLHS = maxAngmom 
   Alloc%maxPrimLHS    = maxprim
   Alloc%maxPrimTUVLHS = max(MAXPRIMTUV,MAXTOTORB,MAXPRIMIJK)
   Alloc%maxnAngLHS     = maxnAngmom
   Alloc%maxContLHS    = maxCont
   Alloc%maxAngmomOrbLHS  = max(maxA,maxB,Alloc%maxAngmomOrbLHS)
   Alloc%maxAngmomLHS     = max(maxAngmom,Alloc%maxAngmomLHS)
   Alloc%maxTUVLHS        = max(maxTUV,Alloc%maxTUVLHS)
!   Alloc%maxPrimTUVijkLHS = max(MAXPRIMTUVIJK,Alloc%maxPrimTUVijkLHS) 
   Alloc%maxTotOrbLHS = max(maxTotorb,Alloc%maxTotOrbLHS)
   Alloc%maxijkLHS = max(maxijk2,Alloc%maxijkLHS)
   Alloc%maxETUVlenLHS = MAXPRIMTUVIJK!maxlenETUV
ELSE
!   Alloc%maxJRHS = maxAngmom 
   Alloc%maxPrimRHS    = maxprim
   Alloc%maxPrimTUVRHS = max(MAXPRIMTUV,MAXTOTORB,MAXPRIMIJK)
   Alloc%maxnAngRHS     = maxnAngmom
   Alloc%maxContRHS    = maxCont
   Alloc%maxAngmomOrbRHS  = max(maxA,maxB,Alloc%maxAngmomOrbRHS)
   Alloc%maxAngmomRHS     = max(maxAngmom,Alloc%maxAngmomRHS)
   Alloc%maxTUVRHS        = max(maxTUV,Alloc%maxTUVRHS)
!   Alloc%maxPrimTUVijkRHS = max(MAXPRIMTUVIJK,Alloc%maxPrimTUVijkRHS) 
   Alloc%maxTotOrbRHS = max(maxTotorb,Alloc%maxTotOrbRHS)
   Alloc%maxijkRHS = max(maxijk2,Alloc%maxijkRHS)
   Alloc%maxETUVlenRHS = MAXPRIMTUVIJK!maxlenETUV
ENDIF


END SUBROUTINE SET_ALLOC

!> \brief initialise the overlap
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param PassP The overlap
!> \param Alloc Information about the maximal sizes for allocation purposes
!> \param Input contains info about the integral evaluation requested
!> \param ODbat The OD-batch information
!> \param side The side ('LHS' or 'RHS') 
!> \param maxpasses the maximum number of passes
!> \param iprint The print level (the higher the more information)
!> \param lupri The default print unit
SUBROUTINE INIT_PASS_FROM_OVERLAP(PassQ,Q,Input,maxpasses,IPRINT,LUPRI)
Implicit none
INTEGER             :: lupri,IPRINT
TYPE(Overlap)       :: PassQ,Q
TYPE(IntegralInput) :: Input
!
Integer             :: maxpasses
Integer             :: ma,mp,na,np,ng,nETUV
!ielectron = 3 hardcoded => PASS
np = Q%Orbital1%nPrimitives
na = Q%Orbital1%nAngmom
ma = Q%Orbital1%maxAngmom
CALL ALLOC_ORBITAL(PassQ%orbital1,np,maxpasses,na,ma,.FALSE.,3)
np = Q%Orbital2%nPrimitives
na = Q%Orbital2%nAngmom
ma = Q%Orbital2%maxAngmom
CALL ALLOC_ORBITAL(PassQ%orbital2,np,maxpasses,na,ma,.FALSE.,3)
call mem_ODpointer_alloc(PassQ%Orb1atom,maxpasses,3)
call mem_ODpointer_alloc(PassQ%Orb1mol,maxpasses,3)
call mem_ODpointer_alloc(PassQ%Orb1batch,maxpasses,3)
call mem_ODpointer_alloc(PassQ%Orb2atom,maxpasses,3)
call mem_ODpointer_alloc(PassQ%Orb2mol,maxpasses,3)
call mem_ODpointer_alloc(PassQ%Orb2batch,maxpasses,3)
np = Q%nPrimitives
mp = np*maxpasses
na = Q%nAngmom
ng = Q%endGeoOrder
CALL ALLOC_OVERLAP(PassQ,mp,maxpasses,na,.FALSE.,1,ng,1,3)

END SUBROUTINE INIT_PASS_FROM_OVERLAP

!> \brief initialise the overlap
!> \author S. Reine and T. Kjaergaard
!> \date 2009
!> \param PassP The overlap
!> \param Alloc Information about the maximal sizes for allocation purposes
!> \param Input contains info about the integral evaluation requested
!> \param ODbat The OD-batch information
!> \param side The side ('LHS' or 'RHS') 
!> \param maxpasses the maximum number of passes
!> \param iprint The print level (the higher the more information)
!> \param lupri The default print unit
SUBROUTINE MEM_PASS_FROM_OVERLAP(Q,Input,maxpasses)
Implicit none
INTEGER,intent(in)             :: maxpasses
TYPE(Overlap),intent(in)       :: Q
TYPE(IntegralInput),intent(in) :: Input
!
Integer             :: ma,mp,na,np,ng,nETUV
integer(kind=long)  :: nint,nrealk
nint = 0
nrealk = 0
np = Q%Orbital1%nPrimitives
na = Q%Orbital1%nAngmom
ma = Q%Orbital1%maxAngmom
CALL MEM_ORBITAL(np,maxpasses,na,ma,.FALSE.,nrealk,nint)
np = Q%Orbital2%nPrimitives
na = Q%Orbital2%nAngmom
ma = Q%Orbital2%maxAngmom
CALL MEM_ORBITAL(np,maxpasses,na,ma,.FALSE.,nrealk,nint)
nint = nint + 6*maxpasses
np = Q%nPrimitives
mp = np*maxpasses
na = Q%nAngmom
ng = Q%endGeoOrder
CALL MEM_OVERLAP1(mp,maxpasses,na,.FALSE.,1,ng,1,nrealk,nint)
call ADD_BUFCOUNTERS(3,nint,nrealk)

END SUBROUTINE MEM_PASS_FROM_OVERLAP

INTEGER FUNCTION getTotalGeoComp(derOrder,LHS,RHS,singleP,singleQ,emptyP,emptyQ)
implicit none
Integer,intent(IN) :: derOrder
Logical,intent(IN) :: singleP,singleQ,LHS,RHS,emptyP,emptyQ
Integer :: iOrder,geoComp,startLHS,endLHS
geoComp=0
startLHS = 0
endLHS   = derOrder
IF (.NOT.LHS.OR.emptyP) endLHS   = 0
IF (.NOT.RHS.OR.emptyQ) startLHS = derORder

DO iOrder=startLHS,endLHS
  geoComp = geoComp + getODgeoComp(iOrder,singleP)*getODgeoComp(derOrder-iOrder,singleQ)
ENDDO
getTotalGeoComp = geoComp
END FUNCTION getTotalGeoComp

INTEGER FUNCTION getODgeoComp(derOrder,single)
implicit none
Integer,intent(IN) :: derOrder
Logical,intent(IN) :: single
IF (derOrder.EQ.0) THEN
  getODgeoComp = 1
ELSE IF (derOrder.EQ.1) THEN
  IF (single) THEN 
    getODgeoComp = 3       !(100 010 001)(000)
  ELSE
    getODgeoComp = 6       !(100 010 001)(000) (000)(100 010 001)
  ENDIF
ELSE IF (derOrder.EQ.2) THEN
  IF (single) THEN 
    getODgeoComp = 6       !(200 110 101 020 011 002)(000)
  ELSE
    getODgeoComp = 6+9+6   !(200 110 101 020 011 002)(000) (100 010 001)(100 010 001) (000)(200 110 101 020 011 002)
  ENDIF
ELSE IF (derOrder.EQ.3) THEN
  IF (single) THEN 
    getODgeoComp = 10            !(300 210 201 120 111 102 030 021 012 003)(000)
  ELSE
    getODgeoComp = 10+18+18+10   !(x3)(x0) (x2)(x1) (x1)(x2) (x0)(x3)
  ENDIF
ELSE IF (derOrder.EQ.4) THEN
  IF (single) THEN 
    getODgeoComp = 15            !(400 310 301 220 211 202 130 121 112 103 040 031 022 013 004)(000)
  ELSE
    getODgeoComp = 15+30+36+30+15   !(x4)(x0) (x3)(x1) (x2)(x2) (x1)(x3) (x0)(x4)
  ENDIF
ELSE
  CALL LSQUIT('Error in function getODgeoComp',-1)
ENDIF
END FUNCTION getODgeoComp


END MODULE Thermite_OD
