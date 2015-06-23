!> @file 
!> Contains integral distribution routines that places calculated integrals in the proper output
!> Thermite integral distribution module
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
MODULE thermite_distributeGen
  use TYPEDEF
  use READMOLEFILE
  use BuildBasisSet
!  use ODbatches
  use OD_type
  use sphcart_matrices
  use precision
  use lstiming
  use thermite_integrals
  use thermite_OD
  use LSTENSOR_OPERATIONSMOD
  use OverlapType

TYPE derivativeInfo
Integer :: Atom(4)
!For each ngeoDerivcompP this returns the directional component (grad x,y,z=1,2,3, hessian xx,xy,xz,yy,yz,zz=1,2,3,4,5,6, etc.)
Integer :: ngeoderivcomp
Integer :: derivComp
Integer, pointer :: dirComp(:,:) 
!For each ngeoDerivcomp this returns the AO the contribution is added to (1-4 mean A,B,C,D, 1,4 means AD, etc.)
Integer, pointer :: AO(:,:)
Integer :: translate      !0 if no translational symmetry, n to translate other contibutions to n
END TYPE derivativeInfo

CONTAINS
!> \brief Initialize the derivative info
!> \author S. Reine
!> \date 2010-03-23
!> \param derivInfo The derivative info
!> \param PQ The integrand
!> \param Input The integral input
SUBROUTINE initDerivativeOverlapInfo(derivInfo,PQ,Input)
implicit none
Type(integrand),intent(in)          :: PQ
Type(IntegralInput),intent(in)      :: Input
Type(derivativeInfo), intent(INOUT) :: derivInfo
!
Logical :: emptyA,emptyB,emptyC,emptyD,singleP,singleQ,emptyP,emptyQ
Integer :: iP,iQ,nP,nQ,n,derP,derQ,nDer,iDerP,iDer,iDeriv,startDer,endDer

Logical :: dohodi

dohodi = (input%geoderOrderP.GT.0) .AND. (input%geoderOrderQ .GT.0)
!
derivInfo%atom = 0
!
!Set up translation (if no translation translate=0)
IF ((PQ%P%p%type_nucleus.OR.PQ%Q%p%type_nucleus).AND..NOT.dohodi) THEN
!  If one side is equal to nuclei we exploit translational symmetry
  IF (PQ%P%p%type_nucleus) THEN
    IF (PQ%P%p%orbital1%type_nucleus) THEN
      derivInfo%translate = 1
    ELSE
      derivInfo%translate = 2
    ENDIF
  ELSE
    IF (PQ%Q%p%orbital1%type_nucleus) THEN
      derivInfo%translate = 3
    ELSE
      derivInfo%translate = 4
    ENDIF
  ENDIF
ELSE
  derivInfo%translate = 0
ENDIF
!
emptyA=PQ%P%p%orbital1%type_empty
emptyB=PQ%P%p%orbital2%type_empty
emptyC=PQ%Q%p%orbital1%type_empty
emptyD=PQ%Q%p%orbital2%type_empty
emptyP  = emptyA.AND.emptyB
emptyQ  = emptyC.AND.emptyD
singleP = (emptyA.OR.emptyB).AND..NOT.(emptyA.AND.emptyB)
singleQ = (emptyC.OR.emptyD).AND..NOT.(emptyC.AND.emptyD)

derP = Input%geoderorderP
derQ = Input%geoderorderQ
nDer = Input%geoderivorder
n    = getTotalGeoComp(nDer,derP.GT.0,derQ.GT.0,singleP,singleQ,emptyP,emptyQ)
call getInputDerivativeInfo(nDer,startDer,endDer,input,emptyP,emptyQ)
derivInfo%ngeoderivcomp = nDer
derivInfo%derivComp     = n
call mem_alloc(derivInfo%dirComp,nDer,n)
call mem_alloc(derivInfo%AO,nDer,n)

iDeriv = 0
DO iDer=endDer,startDer,-1
  iP = iDer
  iQ = nDer - iDer
  IF (iP.GT.derP) CYCLE
  IF (iQ.GT.derQ) CYCLE

  IF (emptyA.AND.emptyB) THEN
    CALL setDerivativeComponents(derivInfo,nDer,iDeriv,0,0,iQ,emptyC,emptyD)
  ELSEIF (emptyB) THEN
    CALL setDerivativeComponents(derivInfo,nDer,iDeriv,iP,0,iQ,emptyC,emptyD)
  ELSEIF (emptyA) THEN
    CALL setDerivativeComponents(derivInfo,nDer,iDeriv,0,iP,iQ,emptyC,emptyD)
  ELSE
    DO iDerP=0,iP
      CALL setDerivativeComponents(derivInfo,nDer,iDeriv,iP-iDerP,iDerP,iQ,emptyC,emptyD)
    ENDDO
  ENDIF
ENDDO
END SUBROUTINE initDerivativeOverlapInfo

!> \brief Print the derivative info
!> \author S. Reine
!> \date 2010-03-23
!> \param derivInfo The derivative info
!> \param LUPRI Default output unit
SUBROUTINE printDerivativeOverlapInfo(derivInfo,LUPRI)
implicit none
TYPE(derivativeInfo), intent(IN) :: derivInfo
Integer,intent(IN)               :: LUPRI
!
Integer :: iDer,ngeoderivcomp,derComp
character(len=80) :: printformat

CALL LSHEADER(LUPRI,'Derivative Info')
IF (derivInfo%atom(1).EQ. 0) THEN
  write(LUPRI,'(3X,A)') 'Atomic numbers not set'
ELSE
  write(LUPRI,'(3X,A,4I5)') 'Atomic numbers: ',derivInfo%atom
ENDIF

IF (derivInfo%translate.EQ. 0) THEN
  write(LUPRI,'(3X,A)') 'No translational symmetry employed'
ELSE
  write(LUPRI,'(3X,A,I2)') 'Translation symmetry employed. Contributions subtracted to AO index number',derivInfo%translate
ENDIF

ngeoderivcomp = derivInfo%ngeoderivcomp
derComp = derivInfo%derivComp
WRITE(LUPRI,'(3X,A,I1,A,I4)') 'Derivative order ',ngeoderivcomp,', number of derivative components', derComp
WRITE(LUPRI,'(5X,A)') 'Derivative index     Directional componants      AO index'
write(printformat,'(A,I2,A,I1,A,I2,A,I1,A)') "(5X,I10,",25-4*ngeoderivcomp,"X,",ngeoderivcomp,"I4,",&
     &                                         25-4*ngeoderivcomp,"X,",ngeoderivcomp,"I4)"
DO iDer=1,derComp
  WRITE(LUPRI,printformat) iDer,derivInfo%dirComp(:,iDer),derivInfo%AO(:,iDer)
ENDDO

END SUBROUTINE printDerivativeOverlapInfo

!> \brief Sets up the derivative component information
!> \author S. Reine
!> \date 2010-03-23
!> \param derivInfo The derivative info
!> \param nDer The derivative order
!> \param iDeriv The derivative component
!> \param derA The derivative order of orbital A
!> \param derB The derivative order of orbital B
!> \param derQ The derivative order of RHS overlap
!> \param emptyC Specifies if orbital C is empty
!> \param emptyD Specifies if orbital D is empty
SUBROUTINE setDerivativeComponents(derivInfo,nDer,iDeriv,derA,derB,derQ,emptyC,emptyD)
implicit none
TYPE(derivativeInfo), intent(INOUT) :: derivInfo
Integer,intent(INOUT)               :: iDeriv
Integer,intent(IN)                  :: nDer,derA,derB,derQ
Logical,intent(IN)                  :: emptyC,emptyD
!
Integer :: xA,yA,zA,xB,yB,zB
DO xA=derA,0,-1
  DO yA=derA-xA,0,-1
    zA=derA-xA-yA
    DO xB=derB,0,-1
      DO yB=derB-xB,0,-1
        zB=derB-xB-yB
        CALL setDerivativeComp1(derivInfo,nDer,iDeriv,derA,xA,yA,zA,derB,xB,yB,zB,derQ,emptyC,emptyD)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE setDerivativeComponents

!> \brief Sets up the derivative component information - continued 1
!> \author S. Reine
!> \date 2010-03-23
!> \param derivInfo The derivative info
!> \param nDer The derivative order
!> \param iDeriv The derivative component
!> \param derA The derivative order of orbital A
!> \param xA The x-directional cartesian component derivative order of orbital A
!> \param yA The y-directional cartesian component derivative order of orbital A
!> \param zA The z-directional cartesian component derivative order of orbital A
!> \param derB The derivative order of orbital B
!> \param xB The x-directional cartesian component derivative order of orbital B
!> \param yB The y-directional cartesian component derivative order of orbital B
!> \param zB The z-directional cartesian component derivative order of orbital B
!> \param derQ The derivative order of RHS overlap
!> \param emptyC Specifies if orbital C is empty
!> \param emptyD Specifies if orbital D is empty
SUBROUTINE setDerivativeComp1(derivInfo,nDer,iDeriv,derA,xA,yA,zA,derB,xB,yB,zB,derQ,emptyC,emptyD)
implicit none
TYPE(derivativeInfo), intent(INOUT) :: derivInfo
Integer,intent(INOUT)               :: iDeriv
Integer,intent(IN)                  :: nDer,derA,xA,yA,zA,derB,xB,yB,zB,derQ
Logical,intent(IN)                  :: emptyC,emptyD
!
Integer :: iDer,endQ
!
endQ = min(derQ,nDer-derA-derB)
IF (emptyC.AND.emptyD) THEN
  CALL setDerivComp2(derivInfo,nDer,iDeriv,derA,xA,yA,zA,derB,xB,yB,zB,0,0)
ELSEIF (emptyD) THEN
  CALL setDerivComp2(derivInfo,nDer,iDeriv,derA,xA,yA,zA,derB,xB,yB,zB,endQ,0)
ELSEIF (emptyC) THEN
  CALL setDerivComp2(derivInfo,nDer,iDeriv,derA,xA,yA,zA,derB,xB,yB,zB,0,endQ)
ELSE
  DO iDer=endQ,0,-1
    CALL setDerivComp2(derivInfo,nDer,iDeriv,derA,xA,yA,zA,derB,xB,yB,zB,iDer,endQ-iDer)
  ENDDO
ENDIF
END SUBROUTINE setDerivativeComp1

!> \brief Sets up the derivative component information - continued 2
!> \author S. Reine
!> \date 2010-03-23
!> \param derivInfo The derivative info
!> \param nDer The derivative order
!> \param iDeriv The derivative component
!> \param derA The derivative order of orbital A
!> \param xA The x-directional cartesian component derivative order of orbital A
!> \param yA The y-directional cartesian component derivative order of orbital A
!> \param zA The z-directional cartesian component derivative order of orbital A
!> \param derB The derivative order of orbital B
!> \param xB The x-directional cartesian component derivative order of orbital B
!> \param yB The y-directional cartesian component derivative order of orbital B
!> \param zB The z-directional cartesian component derivative order of orbital B
!> \param derC The derivative order of orbital C
!> \param derD The derivative order of orbital D
SUBROUTINE setDerivComp2(derivInfo,nDer,iDeriv,derA,xA,yA,zA,derB,xB,yB,zB,derC,derD)
implicit none
TYPE(derivativeInfo), intent(INOUT) :: derivInfo
Integer,intent(INOUT)               :: iDeriv
Integer,intent(IN)                  :: nDer,derA,xA,yA,zA,derB,xB,yB,zB,derC,derD
!
Integer :: xC,yC,zC,xD,yD,zD,direction(0:nDer,0:nDer,0:nDer),dir,iX,iY,iZ,nX,nY,nZ,iAO
!
dir = 0
DO iX=nDer,0,-1
  DO iY=nDer-iX,0,-1
    iZ=nDer-iX-iY
    dir = dir + 1
    direction(iX,iY,iZ) = dir
  ENDDO
ENDDO
!
DO xC=derC,0,-1
  DO yC=derC-xC,0,-1
    zC=derC-xC-yC
    DO xD=derD,0,-1
      DO yD=derD-xD,0,-1
        zD=derD-xD-yD
        iDeriv = iDeriv + 1
        nX = xA+xB+xC+xD
        nY = yA+yB+yC+yD
        nZ = zA+zB+zC+zD
        iAO = 0
        CALL setDerivAO(derivInfo%AO,iAO,xA,iDeriv,1)
        CALL setDerivAO(derivInfo%AO,iAO,yA,iDeriv,1)
        CALL setDerivAO(derivInfo%AO,iAO,zA,iDeriv,1)
        CALL setDerivAO(derivInfo%AO,iAO,xB,iDeriv,2)
        CALL setDerivAO(derivInfo%AO,iAO,yB,iDeriv,2)
        CALL setDerivAO(derivInfo%AO,iAO,zB,iDeriv,2)
        CALL setDerivAO(derivInfo%AO,iAO,xC,iDeriv,3)
        CALL setDerivAO(derivInfo%AO,iAO,yC,iDeriv,3)
        CALL setDerivAO(derivInfo%AO,iAO,zC,iDeriv,3)
        CALL setDerivAO(derivInfo%AO,iAO,xD,iDeriv,4)
        CALL setDerivAO(derivInfo%AO,iAO,yD,iDeriv,4)
        CALL setDerivAO(derivInfo%AO,iAO,zD,iDeriv,4)
        iAO = 0
        CALL setDirComp(derivInfo%dirComp,iAO,xA,iDeriv,1)
        CALL setDirComp(derivInfo%dirComp,iAO,yA,iDeriv,2)
        CALL setDirComp(derivInfo%dirComp,iAO,zA,iDeriv,3)
        CALL setDirComp(derivInfo%dirComp,iAO,xB,iDeriv,1)
        CALL setDirComp(derivInfo%dirComp,iAO,yB,iDeriv,2)
        CALL setDirComp(derivInfo%dirComp,iAO,zB,iDeriv,3)
        CALL setDirComp(derivInfo%dirComp,iAO,xC,iDeriv,1)
        CALL setDirComp(derivInfo%dirComp,iAO,yC,iDeriv,2)
        CALL setDirComp(derivInfo%dirComp,iAO,zC,iDeriv,3)
        CALL setDirComp(derivInfo%dirComp,iAO,xD,iDeriv,1)
        CALL setDirComp(derivInfo%dirComp,iAO,yD,iDeriv,2)
        CALL setDirComp(derivInfo%dirComp,iAO,zD,iDeriv,3)
      ENDDO
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE setDerivComp2

!> \brief Sets up information about the AO-indices for a given deriviative index
!> \author S. Reine
!> \date 2010-03-23
!> \param AO The AO index (1-4) of each deriviative index
!> \param iAO The directional derivative component (1 for gradients, 2 for Hessians, etc.)
!> \param order The cartesian component derivative order
!> \param iDeriv The (full) derivative index
!> \param AOindex The AO index
!> \param AOindex The directional derivative index (x=1,y=2,z=3)
SUBROUTINE setDerivAO(AO,iAO,order,iDeriv,AOindex)
implicit none
Integer,pointer       :: AO(:,:)
Integer,intent(INOUT) :: iAO
Integer,intent(IN)    :: order,iDeriv,AOindex
!
Integer :: i
DO i=1,order
  iAO=iAO+1
  AO(iAO,iDeriv)      = AOindex
ENDDO
END SUBROUTINE setDerivAO

!> \brief Sets up information about the directional cartesian components for given deriviative index
!> \author S. Reine
!> \date 2010-03-23
!> \param dirComp The dirComp index (1-4) of each directional derivative cartesian component and deriviative index
!> \param iAO The directional derivative component (1 for gradients, 2 for Hessians, etc.)
!> \param order The cartesian component derivative order
!> \param iDeriv The (full) derivative index
!> \param dir The directional derivative index (x=1,y=2,z=3)
SUBROUTINE setDirComp(dirComp,iAO,order,iDeriv,dir)
implicit none
Integer,pointer       :: dirComp(:,:)
Integer,intent(INOUT) :: iAO
Integer,intent(IN)    :: order,iDeriv,dir
!
Integer :: i
DO i=1,order
  iAO=iAO+1
  dirComp(iAO,iDeriv) = dir
ENDDO
END SUBROUTINE setDirComp

!> \brief Frees the derivative info
!> \author S. Reine
!> \date 2010-03-23
!> \param derivInfo The derivative info
SUBROUTINE freeDerivativeOverlapInfo(derivInfo)
implicit none
TYPE(derivativeInfo), intent(INOUT) :: derivInfo
derivInfo%translate = 0
derivInfo%ngeoderivcomp = 0
derivInfo%derivComp = 0
call mem_dealloc(derivInfo%dirComp)
call mem_dealloc(derivInfo%AO)
END SUBROUTINE freeDerivativeOverlapInfo

!> \brief new distributePQ to lstensor
!> \author \latexonly T. Kj{\ae}rgaard  \endlatexonly
!> \date 2009-02-05
!> \param RES contains the result lstensor
!> \param PQ contain info about the overlap distributions
!> \param QPmat matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param Input contain info about the requested integral 
!> \param Lsoutput contain info about the requested output, and actual output
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE GeneraldistributePQ(RES,PQ,QPmat2,dimQ,dimP,Input,Lsoutput,LUPRI,IPRINT)
implicit none 
Type(integrand),intent(in)      :: PQ
Type(IntegralInput),intent(in)  :: Input
Type(IntegralOutput),intent(inout) :: Lsoutput
Integer,intent(in)              :: LUPRI,IPRINT
Integer,intent(in)              :: dimQ,dimP
!QPMAT2(dimQ*dimP) Ordering : nC,nD,nA,nB,ndim5,nPasses
!Warning in some cases like when using family basisset dimQ*dimP 
!is not equal to nC*nD*nA*nB*ndim5*nPasses
!so you should use nC*nD*nA*nB*ndim5*nPasses to be sure. 
REAL(REALK),target,intent(in)   :: QPMAT2(:)
TYPE(lstensor),intent(inout)    :: RES
!
REAL(REALK),pointer             :: QPMAT3(:),QPMAT4(:),QPMAT5(:),QPMAT6(:)
type(overlap),pointer :: P,Q
logical :: SamePQ,SameLHSaos,SameRHSaos,SameODs,nucleiP,nucleiQ,add
Logical :: permuteAB,permuteCD,permuteOD,noContraction,RHScontraction
integer :: nAngmomP,ndim5Q,ngeoderivcompP,ngeoderivcompQ,ndim5,nmat,nPasses,iDer,idim5,ndim5P
integer :: iPass,iPassP,iPassQ,iDeriv,iDerivP,iDerivQ,nderiv,ndim5output,iAtom
integer :: nA,nB,nC,nD,batchA,batchB,batchC,batchD,atomA,atomB,atomC,atomD
integer :: indABCD,indABDC,indBACD,indBADC,indCDAB,indCDBA,indDCAB,indDCBA
integer :: ndimMag
Real(realk),pointer :: ABCD(:),ABDC(:),BACD(:),BADC(:),CDAB(:),CDBA(:),DCAB(:),DCBA(:)
Real(realk)  :: factor,center(3,1)
logical :: antiAB,antiCD,AntipermuteAB,AntipermuteCD,translate,same12,same13,same23
Type(derivativeInfo) :: derivInfo
integer :: n1,n2,n3,n4,sA,sB,sC,sD,CMimat,maxBat,maxAng,itrans,nDerivQ
integer :: iAtom1,iAtom2,iAtom3,i1,i2,i3,iPermute,nPermute,nAtoms,iPack
integer,pointer :: dim5(:)
logical,pointer :: negative(:)
integer :: nDimGeo,nTranslate
integer,pointer :: packIndex(:,:,:)

antiAB=.FALSE.
antiCD=.FALSE.
nAngmomP       = PQ%P%p%nAngmom
SamePQ         = PQ%samePQ
SameLHSaos     = INPUT%SameLHSaos
SameRHSaos     = INPUT%SameRHSaos
SameODs        = INPUT%SameODs
!Hack Simen should fix this:
ngeoderivcompP = Input%ngeoderivcomp
ngeoderivcompQ = 1
ndim5P         = ngeoderivcompP
ndim5Q         = ngeoderivcompQ*Input%nCartesianMomentComp
IF(INPUT%CMimat.EQ.0)THEN
   CMimat = INPUT%CMimat
ELSE
   CMimat = INPUT%CMimat
ENDIF
IF(PQ%reverseOrder)THEN
   IF(INPUT%magderOrderQ.NE.1)THEN
      !so far this is the only option when PQ%reverseOrder
      CALL LSQUIT('PQ%reverseOrder.AND.INPUT%magderOrderQ.NE.1',lupri) 
   ENDIF
   antiCD=.TRUE.
   ! we have magnetic differentiated on the RHS
   ! we have reversed the order when we did the integral evaluation. 
   ! Now we change it back to distribute correctly.
   ! so QPmat2 have the ordering nA,nB,nC,nD,ndim5,nPasses
   ! and we need to transpose it to get the ordering
   ! nC,nD,nA,nB,ndim5,nPasses
   P => PQ%Q%p !this points to the actual LHS overlap 
   Q => PQ%P%p !this points to the actual RHS overlap 
   IF(P%magderiv.EQ.1)ndim5P=ndim5P*3
   IF(Q%magderiv.EQ.1)ndim5Q=ndim5Q*3
   !this also mean we need to transpose.
   nA     = P%orbital1%totOrbitals
   nB     = P%orbital2%totOrbitals
   nC     = Q%orbital1%totOrbitals
   nD     = Q%orbital2%totOrbitals
   nPasses = P%nPasses*Q%nPasses
   call mem_workpointer_alloc(QPmat4,nA*nB*nC*nD*3*nPasses)
   call transposeQPorder(QPmat2,QPmat4,nA*nB,nC*nD,3,nPasses,lupri)
ELSEIF(INPUT%magderOrderP.EQ.1)THEN 
   antiAB=.TRUE.
   ! we have done the integral evaluation in the normal order on the LHS
   P => PQ%P%p !this points to the actual LHS overlap 
   Q => PQ%Q%p !this points to the actual RHS overlap 
   IF(P%magderiv.EQ.1)ndim5P=ndim5P*3
   IF(Q%magderiv.EQ.1)ndim5Q=ndim5Q*3
   QPmat4 => QPmat2
ELSE
   P => PQ%P%p
   Q => PQ%Q%p
   QPmat4 => QPmat2
ENDIF
nA     = P%orbital1%totOrbitals
nB     = P%orbital2%totOrbitals
nC     = Q%orbital1%totOrbitals
nD     = Q%orbital2%totOrbitals
!ndim5          = ndim5P*ndim5Q
nderiv         = Input%ngeoderivcomp
ndim5          = Input%ngeoderivcomp*input%nCartesianMomentComp
IF ((P%magderiv.EQ.1).OR.(Q%magderiv.EQ.1)) ndim5 = ndim5*3
!nderiv         = ngeoderivcompP*ngeoderivcompQ 
nDerivQ        = ngeoderivcompQ
ndim5output    = Lsoutput%ndim(5)
nPasses = P%nPasses*Q%nPasses

nDimGeo = 1
nTranslate = 0
translate = .FALSE.
if(input%geoDerivOrder.GE.1)then
   CALL initDerivativeOverlapInfo(derivInfo,PQ,Input)
   IF (IPRINT.GT. 50) THEN
      CALL printDerivativeOverlapInfo(derivInfo,lupri)
   ENDIF
   IF (input%geoDerivOrder.GE.2) nAtoms = res%nAtom(1)
   if(nderiv.GT.1)then
     translate = derivInfo%translate.GT. 0
   endif
   IF (input%geoDerivOrder.EQ.1) THEN
     nDimGeo = 1
   ELSE IF (input%geoDerivOrder.EQ.2) THEN
     nDimGeo = 2
     call mem_alloc(packIndex,3*nAtoms,3*nAtoms,1)
     iPack=0
     DO i1=1,3*nAtoms
       DO i2=i1,3*nAtoms
           iPack=iPack+1
           packIndex(i1,i2,1) = iPack
           packIndex(i2,i1,1) = iPack
       ENDDO
     ENDDO
   ELSE IF (input%geoDerivOrder.EQ.3) THEN
     nDimGeo = 6
     call mem_alloc(packIndex,3*nAtoms,3*nAtoms,3*nAtoms)
     iPack=0
     DO i1=1,3*nAtoms
       DO i2=i1,3*nAtoms
         DO i3=i2,3*nAtoms
           iPack=iPack+1
           packIndex(i1,i2,i3) = iPack
           packIndex(i1,i3,i2) = iPack
           packIndex(i2,i1,i3) = iPack
           packIndex(i2,i3,i1) = iPack
           packIndex(i3,i1,i2) = iPack
           packIndex(i3,i2,i1) = iPack
         ENDDO
       ENDDO
     ENDDO
   ELSE 
     call lsquit('nDimGeo not yet implemented for geoDerivOrder > 3 !',-1)
   ENDIF
   IF (translate) THEN
     IF (input%geoDerivOrder.EQ.1) THEN
       nTranslate = 1
     ELSE IF (input%geoDerivOrder.EQ.2) THEN
       nTranslate = 6
     ELSE IF (input%geoDerivOrder.EQ.3) THEN
       nTranslate = 7
     ELSE
       call lsquit('Translational invariance not yet implemented for geoDerivOrder > 3 !',-1)
     ENDIF
   ENDIF
endif
call mem_alloc(Dim5,nDimGeo+nTranslate)
call mem_alloc(negative,nDimGeo+nTranslate)

if(Input%LinComCarmomType.GT.0)then 
   !magnetic derivative overlaps and other integrals are a 
   !special case which can be written as a 
   !linear combination of cartesian momentum integrals
   call mem_workpointer_alloc(QPmat3,nA*nB*3*nPasses)
   if(Input%LinComCarmomType.EQ.1)then !magnetic derivative overlap integrals
      call GDPQ_magderiv_crossproduct(QPmat4,nA*nB,Input%nCartesianMomentComp,&
           & nPasses,QPmat3,3,P%distance12,lupri)
   elseif(Input%LinComCarmomType.EQ.2)then !LHS half magnetic derivative overlap integrals
      center(1,1) = P%orbital1%center(1)
      center(2,1) = P%orbital1%center(2)
      center(3,1) = P%orbital1%center(3)
      call GDPQ_magderiv_crossproduct(QPmat4,nA*nB,Input%nCartesianMomentComp,&
           & nPasses,QPmat3,3,center,lupri)
   elseif(Input%LinComCarmomType.EQ.3)then !RHS half magnetic derivative overlap integrals
      center(1,1) = - P%orbital2%center(1)
      center(2,1) = - P%orbital2%center(2)
      center(3,1) = - P%orbital2%center(3)
      call GDPQ_magderiv_crossproduct(QPmat4,nA*nB,Input%nCartesianMomentComp,&
           & nPasses,QPmat3,3,center,lupri)
   else
      call lsquit('unknown case in LinComCarmomType',-1)
   endif
   IF(PQ%reverseOrder)call mem_dealloc(QPmat4)
   ndim5 = 3
   ndim5Q= 3
else
   QPmat3 => QPmat4
endif

nmat           = 1
nucleiP        = (P%orbital1%type_nucleus.AND.P%orbital2%type_empty).OR. &
     &           (P%orbital1%type_empty.AND.P%orbital2%type_nucleus)
nucleiQ        = (Q%orbital1%type_nucleus.AND.Q%orbital2%type_empty).OR. &
     &           (Q%orbital1%type_empty.AND.Q%orbital2%type_nucleus)
add = INPUT%AddToIntegral
iPass = 0
DO iPassP=1,P%nPasses
  atomA  = P%orb1atom(iPassP)
  batchA = P%orb1batch(iPassP)
  atomB  = P%orb2atom(iPassP)
  batchB = P%orb2batch(iPassP)
  IF(input%geoDerivOrder.GE.1)then
    derivInfo%Atom(1)=P%orb1mol(iPassP)
    derivInfo%Atom(2)=P%orb2mol(iPassP)
  ELSE
    iDer  = 1
  ENDIF
  IF (nucleiP) THEN
    atomA = 1
    atomB = 1
  ENDIF
  permuteAB  = SameLHSaos.AND.((batchA.NE.batchB).OR.(atomA.NE.atomB))
  AntipermuteAB  = permuteAB.AND.antiAB
  DO iPassQ=1,Q%nPasses
    iPass = iPass + 1
    atomC  = Q%orb1atom(iPassQ)
    batchC = Q%orb1batch(iPassQ)
    atomD  = Q%orb2atom(iPassQ)
    batchD = Q%orb2batch(iPassQ)
    IF(input%geoDerivOrder.GE.1)then
       derivInfo%Atom(3)=Q%orb1mol(iPassQ)
       derivInfo%Atom(4)=Q%orb2mol(iPassQ)
    ELSE
       iDer  = 1
    ENDIF
    IF (nucleiQ) THEN
      atomC = 1
      atomD = 1
    ENDIF
    permuteCD = SameRHSaos.AND.((batchC.NE.batchD).OR.(atomC.NE.atomD))
    permuteOD = SameODs.AND.( ((batchA.NE.batchC).OR.(atomA.NE.atomC)).OR.&
   &                          ((batchB.NE.batchD).OR.(atomB.NE.atomD)) )

    AntipermuteCD  = permuteCD.AND.antiCD
    indABCD = res%index(atomA,atomB,atomC,atomD)
    IF(indABCD.NE.0)THEN 
     ABCD => res%LSAO(indABCD)%elms
     !    ABCD => res%LSAO(indABCD)%batch(batchA,batchB,batchC,batchD)%elms
     n1 = res%LSAO(indABCD)%nLocal(1)
     n2 = res%LSAO(indABCD)%nLocal(2)
     n3 = res%LSAO(indABCD)%nLocal(3)
     n4 = res%LSAO(indABCD)%nLocal(4)
     maxBat = res%LSAO(indABCD)%maxBat
     maxAng = res%LSAO(indABCD)%maxAng
     sA = res%LSAO(indABCD)%startLocalOrb(1+(batchA-1)*maxAng)-1
     sB = res%LSAO(indABCD)%startLocalOrb(1+(batchB-1)*maxAng+maxAng*maxBat)-1 
     sC = res%LSAO(indABCD)%startLocalOrb(1+(batchC-1)*maxAng+2*maxAng*maxBat)-1 
     sD = res%LSAO(indABCD)%startLocalOrb(1+(batchD-1)*maxAng+3*maxAng*maxBat)-1 

     IF (translate) THEN
       call mem_workpointer_alloc(QPmat5,nA*nB*nC*nD*ndim5*nPasses)
       QPmat5 = 0.0E0_realk
       call daxpy(nA*nB*nC*nD*ndim5*nPasses,-1.0E0_realk,QPmat3,1,QPmat5,1)
     ENDIF

     IF (permuteOD) THEN
        indCDAB = res%index(atomC,atomD,atomA,atomB)
        CDAB => res%LSAO(indCDAB)%elms
     ENDIF
     IF (permuteAB) THEN
        indBACD = res%index(atomB,atomA,atomC,atomD)
        BACD => res%LSAO(indBACD)%elms
     ENDIF
     IF (permuteOD.AND.permuteAB) THEN
        indCDBA = res%index(atomC,atomD,atomB,atomA)
        CDBA => res%LSAO(indCDBA)%elms
     ENDIF
     IF (permuteCD) THEN
        indABDC = res%index(atomA,atomB,atomD,atomC)
        ABDC => res%LSAO(indABDC)%elms
     ENDIF
     IF (permuteAB.AND.permuteCD) THEN
        indBADC = res%index(atomB,atomA,atomD,atomC)
        BADC => res%LSAO(indBADC)%elms
     ENDIF
     IF (permuteOD.AND.permuteCD) THEN
        indDCAB = res%index(atomD,atomC,atomA,atomB)
        DCAB => res%LSAO(indDCAB)%elms
     ENDIF
     IF (permuteOD.AND.permuteAB.AND.permuteCD) THEN
        indDCBA = res%index(atomD,atomC,atomB,atomA)
        DCBA => res%LSAO(indDCBA)%elms
     ENDIF
     DO iDeriv = 1,ndim5
       IF(CMimat.EQ.0)THEN
          iDim5=ideriv
       ELSE
          IF(iDeriv .NE. CMimat)CYCLE
          iDim5=1
       ENDIF
       Dim5(1) = iDim5
       nPermute = 1
       IF (translate) iTrans = derivInfo%Atom(derivInfo%translate)
       negative(:) = .FALSE.
       IF (input%geoDerivOrder.EQ. 1)THEN
          iDer = derivInfo%dirComp(1,iDeriv)
          iAtom = derivInfo%Atom(derivInfo%AO(1,iDeriv))
          Dim5(1) = 3*(iAtom-1)+ider
          IF (translate) THEN
            nPermute = nPermute + 1
            Dim5(2) = 3*(iTrans-1)+ider
            negative(2) = .TRUE.
          ENDIF
       ELSEIF (input%geoDerivOrder.EQ. 2)THEN
          iAtom1 = derivInfo%Atom(derivInfo%AO(1,iDeriv))
          iAtom2 = derivInfo%Atom(derivInfo%AO(2,iDeriv))
          i1 = 3*(iAtom1-1)+derivInfo%dirComp(1,iDeriv)
          i2 = 3*(iAtom2-1)+derivInfo%dirComp(2,iDeriv)
          Dim5(1) = packIndex(i1,i2,1)
          Dim5(2) = packIndex(i1,i2,1)
          same12 = (derivInfo%AO(1,iDeriv).EQ.derivInfo%AO(2,iDeriv))
          nPermute=1
          IF ((i1.EQ.i2).AND..NOT.same12) nPermute=2
          IF (translate) THEN
            IF (iAtom1.EQ.iAtom2) THEN
              IF ((.NOT.same12).AND.(iAtom1.EQ.iTrans)) THEN
                nPermute = nPermute - 1
              ENDIF
              Dim5(nPermute+1) = 3*nAtoms*(3*(iTrans-1)+i2-1) + 3*(iAtom1-1)+i1
              Dim5(nPermute+2) = 3*nAtoms*(3*(iAtom2-1)+i2-1) + 3*(iTrans-1)+i1
              Dim5(nPermute+3) = 3*nAtoms*(3*(iTrans-1)+i1-1) + 3*(iTrans-1)+i2
              Dim5(nPermute+4) = 3*nAtoms*(3*(iTrans-1)+i2-1) + 3*(iTrans-1)+i1
              negative(nPermute+1) = .TRUE.
              negative(nPermute+2) = .TRUE.
              nPermute = nPermute + 4
              IF (same12.OR.(iAtom1.EQ.iTrans)) nPermute = nPermute-1
            ELSE
              Dim5(nPermute+1) = 3*nAtoms*(3*(iTrans-1)+i2-1) + 3*(iAtom1-1)+i1
              Dim5(nPermute+2) = 3*nAtoms*(3*(iTrans-1)+i1-1) + 3*(iAtom2-1)+i2
              negative(nPermute+1) = .TRUE.
              negative(nPermute+2) = .TRUE.
              IF (iAtom1.EQ.iTrans.AND.i1.EQ.i2) THEN
                nPermute = nPermute - 1
              ELSE
                Dim5(nPermute+3) = 3*nAtoms*(3*(iAtom1-1)+i1-1) + 3*(iTrans-1)+i2
                negative(nPermute+3) = .TRUE.
              ENDIF
              IF (iAtom2.EQ.iTrans.AND.i1.EQ.i2) THEN
                nPermute = nPermute - 1
              ELSE
                Dim5(nPermute+4) = 3*nAtoms*(3*(iAtom2-1)+i2-1) + 3*(iTrans-1)+i1
                negative(nPermute+4) = .TRUE.
              ENDIF
              Dim5(nPermute+5) = 3*nAtoms*(3*(iTrans-1)+i2-1) + 3*(iTrans-1)+i1
              IF (i1.EQ.i2.AND.((iAtom1.EQ.iTrans).OR.(iAtom2.EQ.iTrans))) THEN
                nPermute = nPermute - 1
              ELSE
                Dim5(nPermute+6) = 3*nAtoms*(3*(iTrans-1)+i1-1) + 3*(iTrans-1)+i2
              ENDIF
              nPermute = nPermute + 6
            ENDIF
          ENDIF
       ELSEIF (input%geoDerivOrder.EQ. 3)THEN
          iAtom1 = derivInfo%Atom(derivInfo%AO(1,iDeriv))
          iAtom2 = derivInfo%Atom(derivInfo%AO(2,iDeriv))
          iAtom3 = derivInfo%Atom(derivInfo%AO(3,iDeriv))
          i1 = 3*(iAtom1-1)+derivInfo%dirComp(1,iDeriv)
          i2 = 3*(iAtom2-1)+derivInfo%dirComp(2,iDeriv)
          i3 = 3*(iAtom3-1)+derivInfo%dirComp(3,iDeriv)
 
          nPermute = 1
          Dim5(1)  = packIndex(i1,i2,i3)
          Dim5(2)  = packIndex(i1,i2,i3)
          Dim5(3)  = packIndex(i1,i2,i3)
          Dim5(4)  = packIndex(i1,i2,i3)
          Dim5(5)  = packIndex(i1,i2,i3)
          Dim5(6)  = packIndex(i1,i2,i3)
          same12 = (derivInfo%AO(1,iDeriv).EQ.derivInfo%AO(2,iDeriv))
          same13 = (derivInfo%AO(1,iDeriv).EQ.derivInfo%AO(3,iDeriv))
          same23 = (derivInfo%AO(2,iDeriv).EQ.derivInfo%AO(3,iDeriv))
          IF ((i1.EQ.i2).AND.(i1.EQ.i3)) THEN
             IF (same12.AND.same13) THEN !All the same AO index
               nPermute=1
             ELSE IF (same12.OR.same13.OR.same23) THEN !Two are the same AO index
               nPermute=3
             ELSE !All different AO indeces
               nPermute=6
             ENDIF
          ELSE IF (i1.EQ.i2) THEN
             nPermute=2
             IF (same12) nPermute=1
          ELSE IF (i1.EQ.i3) THEN
             nPermute=2
             IF (same13) nPermute=1
          ELSE IF (i2.EQ.i3) THEN
             nPermute=2
             IF (same23) nPermute=1
          ENDIF
          IF (translate) THEN
            call lsquit('Error in generalDistributePQ - geoDerivORder 3 and translate not implemented',-1)
          ENDIF
       ELSE IF (input%geoDerivOrder.GE. 4)THEN
         call lsquit('Error in generalDistributePQ - geoDerivORder > 3 not implemented',-1)
       ENDIF
       DO iPermute=1,nPermute
        iDim5 = Dim5(iPermute)
        IF (negative(iPermute)) THEN
          QPmat6 => QPmat5
        ELSE
          QPmat6 => QPmat3
        ENDIF
        IF (permuteOD.AND.permuteCD.AND.permuteAB) THEN
           CALL GDPQ_intABCD(ABCD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intBACD(BACD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intABDC(ABDC,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intBADC(BADC,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intCDAB(CDAB,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intDCAB(DCAB,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intCDBA(CDBA,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intDCBA(DCBA,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
        ELSE IF (permuteOD.AND.permuteCD) THEN
           CALL GDPQ_intABCD(ABCD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intABDC(ABDC,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intCDAB(CDAB,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intDCAB(DCAB,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
        ELSE IF (permuteOD.AND.permuteAB) THEN
           CALL GDPQ_intABCD(ABCD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intBACD(BACD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intCDAB(CDAB,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intCDBA(CDBA,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
        ELSE IF (antipermuteAB.AND.permuteCD) THEN
           CALL GDPQ_intABCD(ABCD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_AntiIntBACD(BACD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intABDC(ABDC,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_AntiIntBADC(BADC,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
        ELSE IF (permuteAB.AND.antipermuteCD) THEN
           CALL GDPQ_intABCD(ABCD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_IntBACD(BACD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_AntiIntABDC(ABDC,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_AntiIntBADC(BADC,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
        ELSE IF (permuteAB.AND.permuteCD) THEN
           CALL GDPQ_intABCD(ABCD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intBACD(BACD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intABDC(ABDC,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intBADC(BADC,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
        ELSE IF (antipermuteAB) THEN
           CALL GDPQ_intABCD(ABCD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_AntiIntBACD(BACD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
        ELSE IF (permuteAB) THEN
           CALL GDPQ_intABCD(ABCD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intBACD(BACD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
        ELSE IF (antipermuteCD) THEN
           CALL GDPQ_intABCD(ABCD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_AntiIntABDC(ABDC,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
        ELSE IF (permuteCD) THEN
           CALL GDPQ_intABCD(ABCD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intABDC(ABDC,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
        ELSE IF (permuteOD) THEN
           CALL GDPQ_intABCD(ABCD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
           CALL GDPQ_intCDAB(CDAB,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
        ELSE
           CALL GDPQ_intABCD(ABCD,n1,n2,n3,n4,sA,sB,sC,sD,&
                & QPmat6,add,nA,nB,nC,nD,ndim5output,idim5,iPass,ideriv,nPasses,ndim5,lupri)
        ENDIF
       ENDDO !Permute
#ifdef VAR_LSDEBUGINT
       IF (IPRINT.GT.40) THEN
          IF (permuteOD.AND.permuteCD.AND.permuteAB) THEN
             CALL GDPQ_printInt(ABCD,'ABCD',n1,n2,n3,n4,ndim5output,lupri)
             CALL GDPQ_printInt(BACD,'BACD',n2,n1,n3,n4,ndim5output,lupri)
             CALL GDPQ_printInt(ABDC,'ABDC',n1,n2,n4,n3,ndim5output,lupri)
             CALL GDPQ_printInt(BADC,'BADC',n2,n1,n4,n3,ndim5output,lupri)
             CALL GDPQ_printInt(CDAB,'CDAB',n3,n4,n1,n2,ndim5output,lupri)
             CALL GDPQ_printInt(CDBA,'CDBA',n3,n4,n2,n1,ndim5output,lupri)
             CALL GDPQ_printInt(DCAB,'DCAB',n4,n3,n1,n2,ndim5output,lupri)
             CALL GDPQ_printInt(DCBA,'DCBA',n4,n3,n2,n1,ndim5output,lupri)
          ELSE IF (permuteOD.AND.permuteCD) THEN
             CALL GDPQ_printInt(ABCD,'ABCD',n1,n2,n3,n4,ndim5output,lupri)
             CALL GDPQ_printInt(ABDC,'ABDC',n1,n2,n4,n3,ndim5output,lupri)
             CALL GDPQ_printInt(CDAB,'CDAB',n3,n4,n1,n2,ndim5output,lupri)
             CALL GDPQ_printInt(DCAB,'DCAB',n4,n3,n1,n2,ndim5output,lupri)
          ELSE IF (permuteOD.AND.permuteAB) THEN
             CALL GDPQ_printInt(ABCD,'ABCD',n1,n2,n3,n4,ndim5output,lupri)
             CALL GDPQ_printInt(BACD,'BACD',n2,n1,n3,n4,ndim5output,lupri)
             CALL GDPQ_printInt(CDAB,'CDAB',n3,n4,n1,n2,ndim5output,lupri)
             CALL GDPQ_printInt(CDBA,'CDBA',n3,n4,n2,n1,ndim5output,lupri)
          ELSE IF (antipermuteAB.AND.permuteCD) THEN
             CALL GDPQ_printInt(ABCD,'ABCD',n1,n2,n3,n4,ndim5output,lupri)
             CALL GDPQ_printInt(BACD,'BACD',n2,n1,n3,n4,ndim5output,lupri)
             CALL GDPQ_printInt(ABDC,'ABDC',n1,n2,n4,n3,ndim5output,lupri)
             CALL GDPQ_printInt(BADC,'BADC',n2,n1,n4,n3,ndim5output,lupri)
          ELSE IF (permuteAB.AND.antipermuteCD) THEN
             CALL GDPQ_printInt(ABCD,'ABCD',n1,n2,n3,n4,ndim5output,lupri)
             CALL GDPQ_printInt(BACD,'BACD',n2,n1,n3,n4,ndim5output,lupri)
             CALL GDPQ_printInt(ABDC,'ABDC',n1,n2,n4,n3,ndim5output,lupri)
             CALL GDPQ_printInt(BADC,'BADC',n2,n1,n4,n3,ndim5output,lupri)
          ELSE IF (permuteAB.AND.permuteCD) THEN
             CALL GDPQ_printInt(ABCD,'ABCD',n1,n2,n3,n4,ndim5output,lupri)
             CALL GDPQ_printInt(BACD,'BACD',n2,n1,n3,n4,ndim5output,lupri)
             CALL GDPQ_printInt(ABDC,'ABDC',n1,n2,n4,n3,ndim5output,lupri)
             CALL GDPQ_printInt(BADC,'BADC',n2,n1,n4,n3,ndim5output,lupri)
          ELSE IF (antipermuteAB) THEN
             CALL GDPQ_printInt(ABCD,'ABCD',n1,n2,n3,n4,ndim5output,lupri)
             CALL GDPQ_printInt(BACD,'BACD',n2,n1,n3,n4,ndim5output,lupri)
          ELSE IF (permuteAB) THEN
             CALL GDPQ_printInt(ABCD,'ABCD',n1,n2,n3,n4,ndim5output,lupri)
             CALL GDPQ_printInt(BACD,'BACD',n2,n1,n3,n4,ndim5output,lupri)
          ELSE IF (antipermuteCD) THEN
             CALL GDPQ_printInt(ABCD,'ABCD',n1,n2,n3,n4,ndim5output,lupri)
             CALL GDPQ_printInt(ABDC,'ABDC',n1,n2,n4,n3,ndim5output,lupri)
          ELSE IF (permuteCD) THEN
             CALL GDPQ_printInt(ABCD,'ABCD',n1,n2,n3,n4,ndim5output,lupri)
             CALL GDPQ_printInt(ABDC,'ABDC',n1,n2,n4,n3,ndim5output,lupri)
          ELSE IF (permuteOD) THEN
             CALL GDPQ_printInt(ABCD,'ABCD',n1,n2,n3,n4,ndim5output,lupri)
             CALL GDPQ_printInt(CDAB,'CDAB',n3,n4,n1,n2,ndim5output,lupri)
          ELSE
             CALL GDPQ_printInt(ABCD,'ABCD',n1,n2,n3,n4,ndim5output,lupri)
          ENDIF
       ENDIF !IPRINT
#endif
!     ENDDO !iDerivQ
!    ENDDO !iDerivP
     ENDDO !iDeriv
     IF (translate) THEN
       call mem_workpointer_dealloc(QPmat5)
     ENDIF
    ENDIF 
  ENDDO !iPassQ
ENDDO !iPassP
IF (nderiv.GT. 1) CALL freeDerivativeOverlapInfo(derivInfo)
call mem_dealloc(negative)
call mem_dealloc(Dim5)
IF(Input%LinComCarmomType.GT.0)THEN
   call mem_workpointer_dealloc(QPmat3)
ELSEIF(PQ%reverseOrder)THEN
   call mem_workpointer_dealloc(QPmat4)
ENDIF
IF (input%geoDerivOrder.GE.2) call mem_dealloc(packIndex)

end SUBROUTINE GeneraldistributePQ

SUBROUTINE GDPQ_printInt(ABCD,text,n1,n2,n3,n4,n5,lupri)
  implicit none
  Integer,intent(IN)          :: n1,n2,n3,n4,n5,lupri
  Character(len=4),intent(IN) :: text
  Real(realk),intent(IN)      :: ABCD(n1,n2,n3,n4,n5)
  !
  Integer :: i1,i2,i3,i4,i5
  WRITE(lupri,'(A,A,5I3)') 'Integrals from generalDistributePQ ',text,n1,n2,n3,n4,n5
  DO i5=1,n5
  DO i4=1,n4
  DO i3=1,n3
  DO i2=1,n2
    write(lupri,'(A,4I3)') '  --integrals',i2,i3,i4,i5
    write(lupri,'(4X,10F17.12)') (ABCD(i1,i2,i3,i4,i5),i1=1,n1)
  ENDDO
  ENDDO
  ENDDO
  ENDDO
END SUBROUTINE GDPQ_printInt

SUBROUTINE GDPQ_intABCD(ABCD,n1,n2,n3,n4,sA,sB,sC,sD,CDAB,add,nA,nB,nC,nD,ndim5output,idim5output,iPass,iDeriv,nPasses,ndim5,lupri)
  implicit none
  Integer,intent(IN)      :: nA,nB,nC,nD,iPass,iDeriv,nPasses,ndim5,idim5output,ndim5output,lupri
  Integer,intent(IN)      :: n1,n2,n3,n4,sA,sB,sC,sD
  Real(realk),intent(IN)  :: CDAB(nC,nD,nA,nB,ndim5,nPasses)
  Real(realk),intent(INOUT) :: ABCD(n1,n2,n3,n4,ndim5output) !nP,nQ,ndim5output)
  Logical,intent(IN)      :: add
  !
  Integer :: iA,iB,iC,iD
  IF (.not.add) THEN
     DO iD=1,nD
        DO iC=1,nC
           DO iB=1,nB
              DO iA=1,nA
                 ABCD(sA+iA,sB+iB,sC+iC,sD+iD,iDim5output) = CDAB(iC,iD,iA,iB,iDeriv,iPass)                   
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ELSE
     !$OMP CRITICAL (GeneraldistributePQBlock)
     DO iD=1,nD
        DO iC=1,nC
           DO iB=1,nB
              DO iA=1,nA
                 ABCD(sA+iA,sB+iB,sC+iC,sD+iD,iDim5output)=ABCD(sA+iA,sB+iB,sC+iC,sD+iD,iDim5output)+&
                      & CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !$OMP END CRITICAL (GeneraldistributePQBlock)
  ENDIF
END SUBROUTINE GDPQ_intABCD

SUBROUTINE GDPQ_intABDC(ABDC,n1,n2,n3,n4,sA,sB,sC,sD,CDAB,add,nA,nB,nC,nD,ndim5output,idim5output,iPass,iDeriv,nPasses,ndim5,lupri)
  implicit none
  Integer,intent(IN)      :: nA,nB,nC,nD,iPass,iDeriv,nPasses,ndim5,idim5output,ndim5output,lupri
  Integer,intent(IN)      :: n1,n2,n3,n4,sA,sB,sC,sD
  Real(realk),intent(IN)  :: CDAB(nC,nD,nA,nB,ndim5,nPasses)
  Real(realk),intent(INOUT) :: ABDC(n1,n2,n4,n3,ndim5output) !(nP,nD,nC,ndim5output)
  Logical,intent(IN)      :: add
  !
  Integer :: iA,iB,iC,iD
  IF (.NOT.add) THEN
     DO iC=1,nC
        DO iD=1,nD
           DO iB=1,nB
              DO iA=1,nA
                 ABDC(sA+iA,sB+iB,sD+iD,sC+iC,iDim5output) = CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ELSE
     !$OMP CRITICAL (GeneraldistributePQBlock)
     DO iC=1,nC
        DO iD=1,nD
           DO iB=1,nB
              DO iA=1,nA
                 ABDC(sA+iA,sB+iB,sD+iD,sC+iC,iDim5output) = &
                      & ABDC(sA+iA,sB+iB,sD+iD,sC+iC,iDim5output) +CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !$OMP END CRITICAL (GeneraldistributePQBlock)
  ENDIF
END SUBROUTINE GDPQ_intABDC

SUBROUTINE GDPQ_AntiIntABDC(ABDC,n1,n2,n3,n4,sA,sB,sC,sD,CDAB,add,nA,nB,nC,nD,ndim5output,&
     & idim5output,iPass,iDeriv,nPasses,ndim5,lupri)
  implicit none
  Integer,intent(IN)      :: nA,nB,nC,nD,iPass,iDeriv,nPasses,ndim5,idim5output,ndim5output,lupri
  Integer,intent(IN)      :: n1,n2,n3,n4,sA,sB,sC,sD
  Real(realk),intent(IN)  :: CDAB(nC,nD,nA,nB,ndim5,nPasses)
  Real(realk),intent(INOUT) :: ABDC(n1,n2,n4,n3,ndim5output)!(nP,nD,nC,ndim5output)
  Logical,intent(IN)      :: add
  !
  Integer :: iA,iB,iC,iD
  IF (.NOT.add) THEN
     DO iC=1,nC
        DO iD=1,nD
           DO iB=1,nB
              DO iA=1,nA
                 ABDC(sA+iA,sB+iB,sD+iD,sC+iC,iDim5output) = CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ELSE
     !$OMP CRITICAL (GeneraldistributePQBlock)
     DO iC=1,nC
        DO iD=1,nD
           DO iB=1,nB
              DO iA=1,nA
                 ABDC(sA+iA,sB+iB,sD+iD,sC+iC,iDim5output) = &
                      & ABDC(sA+iA,sB+iB,sD+iD,sC+iC,iDim5output) + CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !$OMP END CRITICAL (GeneraldistributePQBlock)
  ENDIF
END SUBROUTINE GDPQ_AntiIntABDC

SUBROUTINE GDPQ_intDCAB(DCAB,n1,n2,n3,n4,sA,sB,sC,sD,CDAB,add,nA,nB,nC,nD,ndim5output,idim5output,iPass,iDeriv,nPasses,ndim5,lupri)
  implicit none
  Integer,intent(IN)      :: nA,nB,nC,nD,iPass,iDeriv,nPasses,ndim5,idim5output,ndim5output,lupri
  Integer,intent(IN)      :: n1,n2,n3,n4,sA,sB,sC,sD
  Real(realk),intent(IN)  :: CDAB(nC,nD,nA,nB,ndim5,nPasses)
  Real(realk),intent(INOUT) :: DCAB(n4,n3,n1,n2,ndim5output)!(nD,nC,nP,ndim5output)
  Logical,intent(IN)      :: add
  !
  Integer :: iA,iB,iC,iD
  IF (.NOT.add) THEN
     DO iC=1,nC
        DO iD=1,nD
           DO iB=1,nB
              DO iA=1,nA
                 DCAB(sD+iD,sC+iC,sA+iA,sB+iB,iDim5output) = CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ELSE
     !$OMP CRITICAL (GeneraldistributePQBlock)
     DO iC=1,nC
        DO iD=1,nD
           DO iB=1,nB
              DO iA=1,nA
                 DCAB(sD+iD,sC+iC,sA+iA,sB+iB,iDim5output) &
                      & = DCAB(sD+iD,sC+iC,sA+iA,sB+iB,iDim5output) + CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !$OMP END CRITICAL (GeneraldistributePQBlock)
  ENDIF
END SUBROUTINE GDPQ_intDCAB

SUBROUTINE GDPQ_intBACD(BACD,n1,n2,n3,n4,sA,sB,sC,sD,CDAB,add,nA,nB,nC,nD,ndim5output,idim5output,iPass,iDeriv,nPasses,ndim5,lupri)
  implicit none
  Integer,intent(IN)      :: nA,nB,nC,nD,iPass,iDeriv,nPasses,ndim5,idim5output,ndim5output,lupri
  Integer,intent(IN)      :: n1,n2,n3,n4,sA,sB,sC,sD
  Real(realk),intent(IN)  :: CDAB(nC,nD,nA,nB,ndim5,nPasses)
  Real(realk),intent(INOUT) :: BACD(n2,n1,n3,n4,ndim5output) !(nB,nA,nQ,ndim5output)
  Logical,intent(IN)      :: add
  !
  Integer :: iA,iB,iC,iD
  IF (.NOT.add) THEN
     DO iB=1,nB
        DO iA=1,nA
           DO iD=1,nD
              DO iC=1,nC
                 BACD(sB+iB,sA+iA,sC+iC,sD+iD,iDim5output) = CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ELSE
     !$OMP CRITICAL (GeneraldistributePQBlock)
     DO iB=1,nB
        DO iA=1,nA
           DO iD=1,nD
              DO iC=1,nC
                 BACD(sB+iB,sA+iA,sC+iC,sD+iD,iDim5output) &
                      &=BACD(sB+iB,sA+iA,sC+iC,sD+iD,iDim5output)+CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !$OMP END CRITICAL (GeneraldistributePQBlock)
  ENDIF
END SUBROUTINE GDPQ_intBACD

SUBROUTINE GDPQ_AntiIntBACD(BACD,n1,n2,n3,n4,sA,sB,sC,sD,CDAB,add,nA,nB,nC,nD,ndim5output,&
     & idim5output,iPass,iDeriv,nPasses,ndim5,lupri)
  implicit none
  Integer,intent(IN)      :: nA,nB,nC,nD,iPass,iDeriv,nPasses,ndim5,idim5output,ndim5output,lupri
  Integer,intent(IN)      :: n1,n2,n3,n4,sA,sB,sC,sD
  Real(realk),intent(IN)  :: CDAB(nC,nD,nA,nB,ndim5,nPasses)
  Real(realk),intent(INOUT) :: BACD(n2,n1,n3,n4,ndim5output)!(nB,nA,nQ,ndim5output)
  Logical,intent(IN)      :: add
  !
  Integer :: iA,iB,iC,iD
  IF (.NOT.add) THEN
     DO iB=1,nB
        DO iA=1,nA
           DO iD=1,nD
              DO iC=1,nC
                 BACD(sB+iB,sA+iA,sC+iC,sD+iD,iDim5output) = CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ELSE
     !$OMP CRITICAL (GeneraldistributePQBlock)
     DO iB=1,nB
        DO iA=1,nA
           DO iD=1,nD
              DO iC=1,nC
                 BACD(sB+iB,sA+iA,sC+iC,sD+iD,iDim5output) &
                      & = BACD(sB+iB,sA+iA,sC+iC,sD+iD,iDim5output) - CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !$OMP END CRITICAL (GeneraldistributePQBlock)
  ENDIF
END SUBROUTINE GDPQ_AntiIntBACD

SUBROUTINE GDPQ_intBADC(BADC,n1,n2,n3,n4,sA,sB,sC,sD,CDAB,add,nA,nB,nC,nD,ndim5output,&
     & idim5output,iPass,iDeriv,nPasses,ndim5,lupri)
  implicit none
  Integer,intent(IN)      :: nA,nB,nC,nD,iPass,iDeriv,nPasses,ndim5,idim5output,ndim5output,lupri
  Integer,intent(IN)      :: n1,n2,n3,n4,sA,sB,sC,sD
  Real(realk),intent(IN)  :: CDAB(nC,nD,nA,nB,ndim5,nPasses)
  Real(realk),intent(INOUT) :: BADC(n2,n1,n4,n3,ndim5output)!(nB,nA,nD,nC,ndim5output)
  Logical,intent(IN)      :: add
  !
  Integer :: iA,iB,iC,iD
  IF (.NOT.add) THEN
     DO iD=1,nD
        DO iC=1,nC
           DO iA=1,nA
              DO iB=1,nB
                 BADC(sB+iB,sA+iA,sD+iD,sC+iC,iDim5output) = CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ELSE
     !$OMP CRITICAL (GeneraldistributePQBlock)
     DO iD=1,nD
        DO iC=1,nC
           DO iA=1,nA
              DO iB=1,nB
                 BADC(sB+iB,sA+iA,sD+iD,sC+iC,iDim5output) &
                      & = BADC(sB+iB,sA+iA,sD+iD,sC+iC,iDim5output) + CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !$OMP END CRITICAL (GeneraldistributePQBlock)
  ENDIF
END SUBROUTINE GDPQ_intBADC

SUBROUTINE GDPQ_AntiIntBADC(BADC,n1,n2,n3,n4,sA,sB,sC,sD,CDAB,add,nA,nB,nC,nD,ndim5output,&
     & idim5output,iPass,iDeriv,nPasses,ndim5,lupri)
  implicit none
  Integer,intent(IN)      :: nA,nB,nC,nD,iPass,iDeriv,nPasses,ndim5,idim5output,ndim5output,lupri
  Integer,intent(IN)      :: n1,n2,n3,n4,sA,sB,sC,sD
  Real(realk),intent(IN)  :: CDAB(nC,nD,nA,nB,ndim5,nPasses)
  Real(realk),intent(INOUT) :: BADC(n2,n1,n4,n3,ndim5output)!(nB,nA,nD,nC,ndim5output)
  Logical,intent(IN)      :: add
  !
  Integer :: iA,iB,iC,iD
  IF (.NOT.add) THEN
     DO iD=1,nD
        DO iC=1,nC
           DO iA=1,nA
              DO iB=1,nB
                 BADC(sB+iB,sA+iA,sD+iD,sC+iC,iDim5output) = CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ELSE
     !$OMP CRITICAL (GeneraldistributePQBlock)
     DO iD=1,nD
        DO iC=1,nC
           DO iA=1,nA
              DO iB=1,nB
                 BADC(sB+iB,sA+iA,sD+iD,sC+iC,iDim5output) &
                      & = BADC(sB+iB,sA+iA,sD+iD,sC+iC,iDim5output) - CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !$OMP END CRITICAL (GeneraldistributePQBlock)
  ENDIF
END SUBROUTINE GDPQ_AntiIntBADC

SUBROUTINE GDPQ_intDCBA(DCBA,n1,n2,n3,n4,sA,sB,sC,sD,CDAB,add,nA,nB,nC,nD,ndim5output,&
     & idim5output,iPass,iDeriv,nPasses,ndim5,lupri)
  implicit none
  Integer,intent(IN)      :: nA,nB,nC,nD,iPass,iDeriv,nPasses,ndim5,idim5output,ndim5output,lupri
  Integer,intent(IN)      :: n1,n2,n3,n4,sA,sB,sC,sD
  Real(realk),intent(IN)  :: CDAB(nC,nD,nA,nB,ndim5,nPasses)
  Real(realk),intent(INOUT) :: DCBA(n4,n3,n2,n1,ndim5output)!(nD,nC,nB,nA,ndim5output)
  Logical,intent(IN)      :: add
  !
  Integer :: iA,iB,iC,iD
  IF (.NOT.add) THEN
     DO iB=1,nB
        DO iA=1,nA
           DO iD=1,nD
              DO iC=1,nC
                 DCBA(sD+iD,sC+iC,sB+iB,sA+iA,iDim5output) = CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ELSE
     !$OMP CRITICAL (GeneraldistributePQBlock)
     DO iB=1,nB
        DO iA=1,nA
           DO iD=1,nD
              DO iC=1,nC
                 DCBA(sD+iD,sC+iC,sB+iB,sA+iA,iDim5output) &
                      & = DCBA(sD+iD,sC+iC,sB+iB,sA+iA,iDim5output) + CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !$OMP END CRITICAL (GeneraldistributePQBlock)
  ENDIF
END SUBROUTINE GDPQ_intDCBA

SUBROUTINE GDPQ_intCDBA(CDBA,n1,n2,n3,n4,sA,sB,sC,sD,CDAB,add,nA,nB,nC,nD,ndim5output,&
     & idim5output,iPass,iDeriv,nPasses,ndim5,lupri)
  implicit none
  Integer,intent(IN)      :: nA,nB,nC,nD,iPass,iDeriv,nPasses,ndim5,idim5output,ndim5output,lupri
  Integer,intent(IN)      :: n1,n2,n3,n4,sA,sB,sC,sD
  Real(realk),intent(IN)  :: CDAB(nC,nD,nA,nB,ndim5,nPasses)
  Real(realk),intent(INOUT) :: CDBA(n3,n4,n2,n1,ndim5output)!(nQ,nB,nA,ndim5output)
  Logical,intent(IN)      :: add
  !
  Integer :: iA,iB,iC,iD
  IF (.NOT.add) THEN
     DO iB=1,nB
        DO iA=1,nA
           DO iD=1,nD
              DO iC=1,nC
                 CDBA(sC+iC,sD+iD,sB+iB,sA+iA,iDim5output) = CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ELSE
     !$OMP CRITICAL (GeneraldistributePQBlock)
     DO iB=1,nB
        DO iA=1,nA
           DO iD=1,nD
              DO iC=1,nC
                 CDBA(sC+iC,sD+iD,sB+iB,sA+iA,iDim5output) &
                      & = CDBA(sC+iC,sD+iD,sB+iB,sA+iA,iDim5output) + CDAB(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !$OMP END CRITICAL (GeneraldistributePQBlock)
  ENDIF
END SUBROUTINE GDPQ_intCDBA

SUBROUTINE GDPQ_intCDAB(CDAB,n1,n2,n3,n4,sA,sB,sC,sD,CDABin,add,nA,nB,nC,nD,ndim5output,&
     & idim5output,iPass,iDeriv,nPasses,ndim5,lupri)
  implicit none
  Integer,intent(IN)      :: nA,nB,nC,nD,iPass,iDeriv,nPasses,ndim5,idim5output,ndim5output,lupri
  Integer,intent(IN)      :: n1,n2,n3,n4,sA,sB,sC,sD
  Real(realk),intent(IN)  :: CDABin(nC,nD,nA,nB,ndim5,nPasses)
  Real(realk),intent(INOUT) :: CDAB(n3,n4,n1,n2,ndim5output) !(nPQ,ndim5output)
  Logical,intent(IN)      :: add
  !
  Integer :: iA,iB,iC,iD

  Integer :: iPQ
  IF (.NOT.add) THEN
     DO iA=1,nA
        DO iB=1,nB
           DO iD=1,nD
              DO iC=1,nC
                 CDAB(sC+iC,sD+iD,sA+iA,sB+iB,iDim5output) = CDABin(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ELSE
     !$OMP CRITICAL (GeneraldistributePQBlock)
     DO iA=1,nA
        DO iB=1,nB
           DO iD=1,nD
              DO iC=1,nC
                 CDAB(sC+iC,sD+iD,sA+iA,sB+iB,iDim5output) &
                      & = CDAB(sC+iC,sD+iD,sA+iA,sB+iB,iDim5output) + CDABin(iC,iD,iA,iB,iDeriv,iPass)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !$OMP END CRITICAL (GeneraldistributePQBlock)
  ENDIF
END SUBROUTINE GDPQ_intCDAB

SUBROUTINE GDPQ_magderiv_crossproduct(QPmatCARMOM,nP,nCarmom,nPasses,QPmatMAG,ndimMag,distanceAB,lupri)
  implicit none
  Integer,intent(IN)      :: nP,nCarmom,nPasses,ndimMag,lupri
  Real(realk),intent(IN)  :: QPmatCARMOM(nP,nCarmom,nPasses),distanceAB(3,nPasses)
  Real(realk),intent(INOUT) :: QPmatMAG(nP,ndimMag,nPasses)
  !
  integer,parameter :: betalist(3)=(/2,3,1/),gammalist(3)=(/3,1,2/)
  integer :: iPassP,X,beta,gamma,beta2,gamma2,IP
  real(realk) :: BETADIST,GAMMADIST
  DO iPassP=1,nPasses
     DO X=1,3
        beta = betalist(X)
        gamma = gammalist(X)
        BETADIST=distanceAB(beta,iPassP)
        GAMMADIST=distanceAB(gamma,iPassP) 
        GAMMA2 = gamma+1 !because the first carmom is normal overlap and then X;Y;Z; moments
        BETA2 = beta+1
        DO IP=1,nP
           QPmatMAG(IP,X,iPassP) =  BETADIST*QPmatCARMOM(iP,GAMMA2,iPassP)&
                &                  -GAMMADIST*QPmatCARMOM(iP,BETA2,iPassP)
        ENDDO
     ENDDO
  ENDDO
end SUBROUTINE GDPQ_magderiv_crossproduct

subroutine transposeQPorder(QPmatIN,QPmatOUT,nQ,nP,ndim5,nPasses,lupri)
  implicit none
  Integer,intent(IN)      :: nP,nQ,ndim5,nPasses,lupri
  Real(realk),intent(IN)  :: QPmatIN(nQ,nP,ndim5*nPasses)
  Real(realk),intent(INOUT) :: QPmatOUT(nP,nQ,ndim5*nPasses)
  !
  Integer :: IP,IQ,I5passes
  DO I5passes=1,npasses*ndim5
     DO IP=1,nP
        DO IQ=1,nQ
           QPmatOUT(iP,iQ,I5Passes)=QPmatIN(IQ,IP,I5Passes)
        ENDDO
     ENDDO
  ENDDO
end subroutine transposeQPorder

!> \brief new distributePQ to lstensor
!> \author \latexonly T. Kj{\ae}rgaard  \endlatexonly
!> \date 2009-02-05
!> \param RES contains the result lstensor
!> \param PQ contain info about the overlap distributions
!> \param QPmat matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param Input contain info about the requested integral 
!> \param Lsoutput contain info about the requested output, and actual output
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE distributePQgrad(RES,PQ,QPmat2,dimQ,dimP,Input,Lsoutput,LUPRI,IPRINT)
implicit none 
Type(integrand),intent(in)      :: PQ
Type(IntegralInput),intent(in)  :: Input
Type(IntegralOutput),intent(inout) :: Lsoutput
Integer,intent(in)              :: LUPRI,IPRINT
Integer,intent(in)              :: dimQ,dimP
REAL(REALK),target,intent(in)   :: QPMAT2(dimQ*dimP)
TYPE(lstensor),intent(inout)    :: RES
!
type(overlap),pointer :: P,Q
logical :: SamePQ,SameLHSaos,SameRHSaos,SameODs,nucleiP,nucleiQ,add,kinetic
Logical :: permuteAB,permuteCD,permuteOD,noContraction,RHScontraction,doGrad,translate
integer :: nAngmomP,ngeoderivcomp,nmat,nPasses
integer :: iPass,iPassP,iPassQ,iDeriv,iDerivP,iDerivQ
integer :: nA,nB,nC,nD,batchA,batchB,batchC,batchD,atomA,atomB,atomC,atomD
!integer :: indABCD,indABDC,indBACD,indBADC,indCDAB,indCDBA,indDCAB,indDCBA
integer :: AC,CA,AB,BA,n1,n2,n3,s1,s2,s3,maxBat,maxAng
Real(realk),pointer :: ABCD(:),ABDC(:),BACD(:),BADC(:),CDAB(:),CDBA(:),DCAB(:),DCBA(:)
Real(realk)  :: factor
Real(realk),pointer :: grad(:),trans(:),DAB(:),DAC(:),DBA(:),DCA(:)
Type(derivativeInfo) :: derivInfo
Real(realk) :: derCont(Input%NDMAT_LHS),center(3,1)
!Real(realk) :: derCont(10)
Integer     :: geoderivorder,iDer,iAtom,iTrans,imat,iLSAO
logical  :: AC_center,antiAB
REAL(REALK),pointer             :: QPMAT3(:)
antiAB=.FALSE.
geoderivorder = Input%geoderivorder
IF(geoderivorder.GT.0)THEN
   CALL initDerivativeOverlapInfo(derivInfo,PQ,Input)
   IF (IPRINT.GT. 50) THEN
      CALL printDerivativeOverlapInfo(derivInfo,LUPRI)
   ENDIF
ENDIF
kinetic        = PQ%kinetic
AC_center = kinetic.OR.INPUT%AC_center
nAngmomP       = PQ%P%p%nAngmom
SamePQ         = PQ%samePQ
SameLHSaos     = INPUT%SameLHSaos
SameRHSaos     = INPUT%SameRHSaos
SameODs        = INPUT%SameODs

ngeoderivcomp  = INPUT%nGeoDerivComp

if(Input%nCartesianMomentComp.GT.1.AND.ngeoderivcomp.GT.1)then   
   print*,'Input%ngeoderivcomp',Input%ngeoderivcomp
   print*,'Input%nCartesianMomentComp',Input%nCartesianMomentComp
   call lsquit('CartesianMoment geo derivatives not implemented yet',lupri)
endif
RHScontraction = .TRUE.
nmat           = 1
if (RHScontraction) nmat = Input%NDMAT_RHS
noContraction  = .NOT. RHScontraction

P => PQ%P%p
Q => PQ%Q%p
IF(INPUT%magderOrderP.EQ.1)THEN 
   antiAB=.TRUE.
ENDIF
nucleiP        = (P%orbital1%type_nucleus.AND.P%orbital2%type_empty).OR. &
     &           (P%orbital1%type_empty.AND.P%orbital2%type_nucleus)
nucleiQ        = (Q%orbital1%type_nucleus.AND.Q%orbital2%type_empty).OR. &
     &           (Q%orbital1%type_empty.AND.Q%orbital2%type_nucleus)
add = INPUT%AddToIntegral

nA     = P%orbital1%totOrbitals
nB     = P%orbital2%totOrbitals
nC     = Q%orbital1%totOrbitals
nD     = Q%orbital2%totOrbitals
nPasses = P%nPasses*Q%nPasses

if(Input%LinComCarmomType.GT.0)then 
   !magnetic derivative overlaps and other integrals are a 
   !special case which can be written as a 
   !linear combination of cartesian momentum integrals
   call mem_alloc(QPmat3,nA*nB*3*nPasses)
   if(Input%LinComCarmomType.EQ.1)then !magnetic derivative overlap integrals
      call magderiv_crossproduct(QPmat2,nA*nB,Input%nCartesianMomentComp,&
           & nPasses,QPmat3,3,P%distance12,lupri)
   elseif(Input%LinComCarmomType.EQ.2)then !LHS half magnetic derivative overlap integrals
      center(1,1) = P%orbital1%center(1)
      center(2,1) = P%orbital1%center(2)
      center(3,1) = P%orbital1%center(3)
      call magderiv_crossproduct(QPmat2,nA*nB,Input%nCartesianMomentComp,&
           & nPasses,QPmat3,3,center,lupri)
   elseif(Input%LinComCarmomType.EQ.3)then !RHS half magnetic derivative overlap integrals
      center(1,1) = - P%orbital2%center(1)
      center(2,1) = - P%orbital2%center(2)
      center(3,1) = - P%orbital2%center(3)
      call magderiv_crossproduct(QPmat2,nA*nB,Input%nCartesianMomentComp,&
           & nPasses,QPmat3,3,center,lupri)
   else
      call lsquit('unknown case in LinComCarmomType',-1)
   endif
else
   QPmat3 => QPmat2
endif

iPass = 0
DO iPassP=1,P%nPasses
  atomA  = P%orb1atom(iPassP)
  batchA = P%orb1batch(iPassP)
  atomB  = P%orb2atom(iPassP)
  batchB = P%orb2batch(iPassP)
  permuteAB  = SameLHSaos.AND.((batchA.NE.batchB).OR.(atomA.NE.atomB))
  DO iPassQ=1,Q%nPasses
    iPass = iPass + 1

    atomC  = Q%orb1atom(iPassQ)
    batchC = Q%orb1batch(iPassQ)
    atomD  = Q%orb2atom(iPassQ)
    batchD = Q%orb2batch(iPassQ)
    permuteCD = SameRHSaos.AND.((batchC.NE.batchD).OR.(atomC.NE.atomD))
    permuteOD = SameODs.AND.( ((batchA.NE.batchC).OR.(atomA.NE.atomC)).OR.&
   &                          ((batchB.NE.batchD).OR.(atomB.NE.atomD)) )

    derivInfo%Atom(1)=atomA
    derivInfo%Atom(2)=atomB
    derivInfo%Atom(3)=atomC
    derivInfo%Atom(4)=atomD
    nmat = Input%NDMAT_LHS

    IF (AC_center) THEN
      AC = Input%LST_DLHS%INDEX(atomA,atomC,1,1)
      DAC => Input%LST_DLHS%LSAO(AC)%elms
      n1 = Input%LST_DLHS%LSAO(AC)%nLocal(1) 
      n3 = Input%LST_DLHS%LSAO(AC)%nLocal(2) 
      maxBat = Input%LST_DLHS%LSAO(AC)%maxBat
      maxAng = Input%LST_DLHS%LSAO(AC)%maxAng
      s1 = Input%LST_DLHS%LSAO(AC)%startLocalOrb(1+(batchA-1)*maxAng) - 1 
      s3 = Input%LST_DLHS%LSAO(AC)%startLocalOrb(1+(batchC-1)*maxAng+maxAng*maxBat)-1 
      IF (permuteOD)THEN
         CA = Input%LST_DLHS%INDEX(atomC,atomA,1,1)
         DCA => Input%LST_DLHS%LSAO(CA)%elms
      ENDIF
    ELSE
       AB = Input%LST_DLHS%INDEX(atomA,atomB,1,1)
       DAB => Input%LST_DLHS%LSAO(AB)%elms
       n1 = Input%LST_DLHS%LSAO(AB)%nLocal(1) 
       n2 = Input%LST_DLHS%LSAO(AB)%nLocal(2) 
       maxBat = Input%LST_DLHS%LSAO(AB)%maxBat
       maxAng = Input%LST_DLHS%LSAO(AB)%maxAng
       s1 = Input%LST_DLHS%LSAO(AB)%startLocalOrb(1+(batchA-1)*maxAng)-1 
       s2 = Input%LST_DLHS%LSAO(AB)%startLocalOrb(1+(batchB-1)*maxAng+maxAng*maxBat)-1 
       IF (permuteAB)THEN
          BA = Input%LST_DLHS%INDEX(atomB,atomA,1,1)
          DBA => Input%LST_DLHS%LSAO(BA)%elms
          IF(size(Input%LST_DLHS%LSAO(BA)%elms).NE.n1*n2)CALL LSQUIT('size eror',-1)
       ENDIF
    ENDIF
    IF((P%magderiv.EQ.1).OR.(Input%LinComCarmomType.GT.0))then 
       grad => RES%LSAO(1)%elms
       DO iDeriv=1,3
          derCont = 0E0_realk
          IF (AC_center)THEN
             ! Kinetic type of integrals <a 0|w|c 0>
             derCont = derCont + sac(QPmat3,DAC,n1,n3,s1,s3,nmat,nA,nC,iPass,iDeriv,nPasses,3,lupri)
             IF (permuteOD)THEN
                derCont = derCont + sca(QPmat3,DCA,n3,n1,s3,s1,nmat,nA,nC,iPass,iDeriv,nPasses,3,lupri)
             ENDIF
          ELSE !AB_center
             ! Nuclear-electronic attraction or overlap type integrals
             derCont = derCont + hab(QPmat3,DAB,n1,n2,s1,s2,nmat,nA,nB,iPass,iDeriv,nPasses,3,lupri)
             IF (permuteAB.AND.AntiAB)THEN
                derCont = derCont + antihba(QPmat3,DBA,n2,n1,s2,s1,nmat,nA,nB,iPass,iDeriv,nPasses,3,lupri)
             ELSEIF(permuteAB)THEN
                derCont = derCont + hba(QPmat3,DBA,n2,n1,s2,s1,nmat,nA,nB,iPass,iDeriv,nPasses,3,lupri)
             ENDIF
          ENDIF
          !$OMP CRITICAL(distributeGrd) 
          DO imat=1,nmat
             grad(iDeriv+(imat-1)*3) = grad(iDeriv+(imat-1)*3) + derCont(imat)
          ENDDO
          !$OMP END CRITICAL(distributeGrd)           
       ENDDO
    ELSE
     DO iDeriv=1,ngeoderivcomp
      
       iAtom = derivInfo%Atom(derivInfo%AO(1,iDeriv))
!The structure of the lstensor when used to store the molecular gradient
!is maybe not the most logical and we refer to discussion in the subroutine
!init_gradientlstensor in lsutil/lstensor_operations.f90
       iLSAO = RES%INDEX(iAtom,1,1,derivInfo%AO(1,iDeriv))
#ifdef VAR_LSDEBUGINT
       IF(iLSAO.EQ. 0)THEN
          WRITE(lupri,*)'iLSAO is zero - for some reason'
          WRITE(lupri,*)'iAtom          ',iAtom
          WRITE(lupri,*)'Center (1-4) = ',derivInfo%AO(1,iDeriv)
          WRITE(lupri,*)'Error in use of gradientlstensor distributePQgrad'
          CALL LSQUIT('Error in use of gradientlstensor distributePQgrad',lupri)
       ENDIF
#endif
       iDer  = derivInfo%dirComp(1,iDeriv)
       grad => RES%LSAO(iLSAO)%elms
       translate = derivInfo%translate.GT. 0
       IF (translate) THEN
         iTrans = derivInfo%Atom(derivInfo%translate)
         !The structure of the lstensor when used to store the molecular gradient
         !is maybe not the most logical and we refer to discussion in the subroutine
         !init_gradientlstensor in lsutil/lstensor_operations.f90
         iLSAO = RES%INDEX(iTrans,1,1,derivInfo%translate)
#ifdef VAR_LSDEBUGINT
         IF(iLSAO.EQ. 0)THEN
          WRITE(lupri,*)'iLSAO is zero - for some reason'
          WRITE(lupri,*)'iTrans         ',iTrans
          WRITE(lupri,*)'Center (1-4) = ',derivInfo%translate
          WRITE(lupri,*)'Error in use of gradientlstensor distributePQgrad'
          CALL LSQUIT('Error in use of gradientlstensor distributePQgrad',lupri)
         ENDIF
#endif
         trans => RES%LSAO(iLSAO)%elms
       ENDIF
       derCont = 0E0_realk
       IF (AC_center)THEN
         ! Kinetic type of integrals <a 0|w|c 0>
         derCont = derCont + sac(QPmat3,DAC,n1,n3,s1,s3,nmat,nA,nC,iPass,iDeriv,nPasses,ngeoderivcomp,lupri)
         IF (permuteOD) derCont = derCont + sca(QPmat3,DCA,n3,n1,s3,s1,nmat,nA,nC,iPass,iDeriv,nPasses,ngeoderivcomp,lupri)
       ELSE
         ! Nuclear-electronic attraction or overlap type integrals
         derCont = derCont + hab(QPmat3,DAB,n1,n2,s1,s2,nmat,nA,nB,iPass,iDeriv,nPasses,ngeoderivcomp,lupri)
         IF (permuteAB) derCont = derCont + hba(QPmat3,DBA,n2,n1,s2,s1,nmat,nA,nB,iPass,iDeriv,nPasses,ngeoderivcomp,lupri)
       ENDIF
      !$OMP CRITICAL(distributeGrd) 
       DO imat=1,nmat
         grad(iDer+(imat-1)*3) = grad(iDer+(imat-1)*3) + derCont(imat)
         IF (translate) trans(iDer+(imat-1)*3) = trans(iDer+(imat-1)*3) - derCont(imat)
       ENDDO
       !$OMP END CRITICAL(distributeGrd) 
     ENDDO !iDeriv
    ENDIF
  ENDDO !iPassQ
ENDDO !iPassP

IF(geoderivorder.GT.0)THEN
   CALL freeDerivativeOverlapInfo(derivInfo)
ENDIF
if(Input%LinComCarmomType.GT.0)then 
   call mem_dealloc(QPmat3)
ENDIF


CONTAINS
FUNCTION hab(ab,dab,n1,n2,s1,s2,nmat,nA,nB,iPass,iDeriv,nPasses,nderivcomp,lupri)
implicit none
Integer,intent(IN)     :: nmat,nA,nB,iPass,iDeriv,nPasses,nderivcomp,lupri,n1,n2,s1,s2
Real(realk),intent(IN) :: ab(nA,nB,nderivcomp,nPasses)
Real(realk),intent(IN) :: dab(n1,n2,nmat)
Real(realk)            :: hab(nmat)
!
Integer :: iA,iB,imat

hab = 0E0_realk
DO imat=1,nmat
  DO iB=1,nB
   DO iA=1,nA
    hab(imat) = hab(imat) + dab(s1+iA,s2+iB,imat)*AB(iA,iB,iDeriv,iPass)
   ENDDO
  ENDDO
ENDDO
END FUNCTION hab

FUNCTION hba(ab,dba,n2,n1,s2,s1,nmat,nA,nB,iPass,iDeriv,nPasses,nderivcomp,lupri)
implicit none
Integer,intent(IN)     :: nmat,nA,nB,iPass,iDeriv,nPasses,nderivcomp,lupri,n1,n2,s1,s2
Real(realk),intent(IN) :: ab(nA,nB,nderivcomp,nPasses)
Real(realk),intent(IN) :: dba(n2,n1,nmat)
Real(realk)            :: hba(nmat)
!
Integer :: iA,iB,imat

hba = 0E0_realk
DO imat=1,nmat
  DO iB=1,nB
    DO iA=1,nA
       hba(imat) = hba(imat) + dba(s2+iB,s1+iA,imat)*AB(iA,iB,iDeriv,iPass)
    ENDDO
  ENDDO
ENDDO
END FUNCTION hba

FUNCTION antihba(ab,dba,n2,n1,s2,s1,nmat,nA,nB,iPass,iDeriv,nPasses,nderivcomp,lupri)
implicit none
Integer,intent(IN)     :: nmat,nA,nB,iPass,iDeriv,nPasses,nderivcomp,lupri,n1,n2,s1,s2
Real(realk),intent(IN) :: ab(nA,nB,nderivcomp,nPasses)
Real(realk),intent(IN) :: dba(n2,n1,nmat)
Real(realk)            :: antihba(nmat)
!
Integer :: iA,iB,imat

antihba = 0E0_realk
DO imat=1,nmat
  DO iB=1,nB
    DO iA=1,nA
     antihba(imat) = antihba(imat) - dba(s2+iB,s1+iA,imat)*AB(iA,iB,iDeriv,iPass)
    ENDDO
  ENDDO
ENDDO
END FUNCTION antihba

FUNCTION sac(ac,dac,n1,n3,s1,s3,nmat,nA,nC,iPass,iDeriv,nPasses,nderivcomp,lupri)
implicit none
Integer,intent(IN)     :: nmat,nA,nC,iPass,iDeriv,nPasses,nderivcomp,lupri,n1,n3,s1,s3
Real(realk),intent(IN) :: ac(nC,nA,nderivcomp,nPasses)
Real(realk),intent(IN) :: dac(n1,n3,nmat)
Real(realk)            :: sac(nmat)
!
Integer :: iA,iC,imat

sac = 0E0_realk
DO imat=1,nmat
  DO iA=1,nA
    DO iC=1,nC
     sac(imat) = sac(imat) + dac(s1+iA,s3+iC,imat)*AC(iC,iA,iDeriv,iPass)
    ENDDO
  ENDDO
ENDDO
END FUNCTION sac

FUNCTION sca(ac,dca,n3,n1,s3,s1,nmat,nA,nC,iPass,iDeriv,nPasses,nderivcomp,lupri)
implicit none
Integer,intent(IN)     :: nmat,nA,nC,iPass,iDeriv,nPasses,nderivcomp,lupri,n3,n1,s3,s1
Real(realk),intent(IN) :: ac(nC,nA,nderivcomp,nPasses)
Real(realk),intent(IN) :: dca(n3,n1,nmat)
Real(realk)            :: sca(nmat)
!
Integer :: iA,iC,imat

sca = 0E0_realk
DO imat=1,nmat
  DO iA=1,nA
   DO iC=1,nC
    sca(imat) = sca(imat) + dca(s3+iC,s1+iA,imat)*AC(iC,iA,iDeriv,iPass)
   ENDDO
  ENDDO
ENDDO
END FUNCTION sca

SUBROUTINE magderiv_crossproduct(QPmatCARMOM,nP,nCarmom,nPasses,QPmatMAG,ndimMag,distanceAB,lupri)
  implicit none
  Integer,intent(IN)      :: nP,nCarmom,nPasses,ndimMag,lupri
  Real(realk),intent(IN)  :: QPmatCARMOM(nP,nCarmom,nPasses),distanceAB(3,nPasses)
  Real(realk),intent(INOUT) :: QPmatMAG(nP,ndimMag,nPasses)
!
  integer,parameter :: betalist(3)=(/2,3,1/),gammalist(3)=(/3,1,2/)
  integer :: iPassP,X,beta,gamma,beta2,gamma2,IP
  real(realk) :: BETADIST,GAMMADIST
  DO iPassP=1,P%nPasses
     DO X=1,3
        beta = betalist(X)
        gamma = gammalist(X)
        BETADIST=distanceAB(beta,iPassP)
        GAMMADIST=distanceAB(gamma,iPassP) 
        GAMMA2 = gamma+1 !because the first carmom is normal overlap and then X;Y;Z; moments
        BETA2 = beta+1
        DO IP=1,nP
           QPmatMAG(IP,X,iPassP) =  BETADIST*QPmatCARMOM(iP,GAMMA2,iPassP)&
                &                  -GAMMADIST*QPmatCARMOM(iP,BETA2,iPassP)
        ENDDO
     ENDDO
  ENDDO
end SUBROUTINE magderiv_crossproduct

END SUBROUTINE distributePQgrad

!> \brief print calculated derivative multipole moments to file
!> \author \latexonly A. Krapp \endlatexonly
!> \date 2010
!> \param PQ contain info about the overlap distributions
!> \param QPmat matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param nMMcomp = dimQ = number of cartesian MM components
!> \param nSPHMAT number of spherical MM components
!> \param Input contain info about the requested integral 
!> \param Output2 contain info about whether or not to use the buffer to store integral before dumping to disk
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> Modified according to new structure - S. Reine 2011-03-17
SUBROUTINE printMMtoFile(PQ,ABmom,nA,nB,nCartMom,nSpherMom,ngeoderivcomp,nPass,INPUT,OUTPUT2,LUPRI,IPRINT)
use memory_handling
implicit none
Type(integrand),intent(in)         :: PQ
Type(IntegralInput),intent(inout)  :: Input
Type(IntegralOutput),intent(inout) :: Output2
Integer,intent(in)                 :: LUPRI,IPRINT
Integer,intent(in)                 :: nA,nB,nCartMom,nSpherMom,ngeoderivcomp,nPass
REAL(REALK),intent(in)             :: ABmom(nCartMom,nA*nB*ngeoderivcomp*nPass)
!
Integer :: nAngmomP,nOrbP,ngeoderivcompP,endOrbP,mmorder
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,I,J
Integer :: iAtom,iDer,startOrb,endOrb,iAngmom,iCart,iSpher
Integer :: batchA,batchB,startA,startB,endA
Integer :: geoderivorder,ideriv,atoma,atomb
Integer :: iPassP,iOrbP,imat,m,ml
Integer :: ICOUNTIND,NBATCH,INPUTSTARTA,INPUTSTARTB
!
REAL(REALK)            :: THRESHOLD
REAL(REALK),pointer    :: SpherMom(:,:)
INTEGER,pointer        :: INDDSTR(:,:,:),IBATCHCont(:,:,:)
Real(realk)            :: tmp,absMom,maxMom,maxMomBatch
!
logical               :: AddBatchCont(nA*nB),AddBatch,dograd,dotriangular
Type(derivativeInfo)  :: derivInfo
Type(overlap),pointer :: P
Integer               :: angA(nA),angB(nB)
Integer               :: orbA(nA,nPass),orbB(nB,nPass),iOrb,iOrbFull,idim2
Integer               :: iAA,iAB,sA,cA,aA,icA,eA,sB,eB,cB,aB,icB

P => PQ%P%p
mmorder  = Input%mmorder
ngeoderivcompP  = Input%ngeoderivcomp

IF (ngeoderivcompP.NE.ngeoderivcomp)   CALL LSQUIT('Error in printMMtoFile. ngeoderivcomp mismatch',lupri)
IF (P%nPasses.ne. 1)      CALL LSQUIT('Error in printMMtoFile. P%nPasses > 1',lupri)
IF (PQ%Q%p%nPasses.GT. 1) CALL LSQUIT('Error in printMMtoFile. Q%nPasses > 1',lupri)


geoderivorder = Input%geoderivorder
IF (geoderivorder.GT. 0) THEN
 CALL initDerivativeOverlapInfo(derivInfo,PQ,Input)
 IF (IPRINT.GT. 50) THEN
  CALL printDerivativeOverlapInfo(derivInfo,LUPRI)
 ENDIF
ENDIF

dograd     = INPUT%DO_GRADIENT
nAngmomP   = P%nAngmom
AddBatch   = .FALSE.
threshold  = INPUT%CS_THRESHOLD*1E-1_realk
IF(INPUT%MM_NOSCREEN) THRESHOLD = 0.0E0_realk
CALL mm_set_screen_threshold(THRESHOLD)


INPUTSTARTA = INPUT%MMstartA
INPUTSTARTB = INPUT%MMstartB

! Make spherical transformation of Cartesian multipole moments
! cSimen Should be made another place - not in the distribution routines!!!
! cSimen I think we should make this transformation in contract_Q
! cSimen by making an orbtial-type 'MULMOM'
call mem_alloc(SpherMom,nSpherMom,nA*nB*ngeoderivcomp*nPass)
SpherMom = 0E0_realk
!DO idim2=1,nA*nB*ngeoderivcomp*nPass
call SPHERICAL_TRANSFORMATION3(ABmom,SpherMom,&
     & nA*nB*ngeoderivcomp*nPass,nCartMom,nSpherMom,mmorder,iprint,lupri)
!ENDDO

!Set up the angular moments for each component
DO iAngmom=1,P%orbital1%nAngmom
  startOrb = P%orbital1%startLocOrb(iAngmom)
  endOrb   = startOrb + P%orbital1%nOrbitals(iAngmom) -1
  angA(startOrb:endOrb) = P%orbital1%angmom(iAngmom)
  DO iPassP=1,P%nPasses
    iOrbfull = INPUTSTARTA + P%orbital1%startOrbital(iAngmom+(iPassP-1)*P%orbital1%nAngmom) - 1
    IF (P%orbital1%type_empty) iOrbfull = -1
    DO iOrb=startOrb,endOrb
      iOrbFull = iOrbFull + 1
      orbA(iOrb,iPassP) = iOrbfull
    ENDDO
  ENDDO
ENDDO

DO iAngmom=1,P%orbital2%nAngmom
  startOrb = P%orbital2%startLocOrb(iAngmom)
  endOrb   = startOrb + P%orbital2%nOrbitals(iAngmom) -1
  angB(startOrb:endOrb) = P%orbital2%angmom(iAngmom)
  DO iPassP=1,P%nPasses
    iOrbfull = INPUTSTARTB + P%orbital2%startOrbital(iAngmom+(iPassP-1)*P%orbital2%nAngmom) - 1
    IF (P%orbital2%type_empty) iOrbfull = -1
    DO iOrb=startOrb,endOrb
      iOrbFull = iOrbFull + 1
      orbB(iOrb,iPassP) = iOrbfull
    ENDDO
  ENDDO
ENDDO

ICOUNTIND = INPUT%MMunique_ID1 !unique identifier ()start with 0 
NBATCH    = INPUT%MMunique_ID2 !BATCH INDEX starts with 0 

call mem_alloc(INDDSTR,nA,nB,P%nPasses)
call mem_alloc(IBATCHCont,nA,nB,P%nPasses)
INDDSTR    = -1
IBATCHCont = -1
!Set up the shell index and ao-pair indexes used for unknown purposes in the fmm code
DO iPassP=1,P%nPasses

  atoma  = P%orb1atom(iPassP)
  batchA = P%orb1batch(iPassP)
  atomb  = P%orb2atom(iPassP)
  batchB = P%orb2batch(iPassP)
  dotriangular  = INPUT%SameLHSaos .AND. ((batchA.EQ.batchB).AND.(atoma.EQ.atomb))

  DO iAngmom=1,P%nAngmom
    iAA = P%indexAng1(iangmom)
    iAB = P%indexAng2(iangmom)
    sA  = P%orbital1%startLocOrb(iAA)
    cA  = P%orbital1%nContracted(iAA)
    aA  = P%orbital1%nOrbComp(iAA)
    cB  = P%orbital2%nContracted(iAB)
    aB  = P%orbital2%nOrbComp(iAB)
    DO icA=1,cA
      eA = sA + aA - 1
      sB = P%orbital2%startLocOrb(iAB)
      DO icB=1,cB
        eB = sB + aB - 1
        maxMomBatch = 0E0_realk
        DO iB=sB,eB
          DO iA=sA,eA
            IF (dotriangular.AND.(iA.GT.iB)) CYCLE
            maxMom      = 0E0_realk
            DO iDerivP=1,ngeoderivcompP
              idim2 = iA+(iB-1)*nA + (iDerivP-1)*nB*nA + (iPassP-1)*ngeoderivcompP*nB*nA
              DO iSpher=1,nSpherMom
                absMom      = abs(spherMom(iSpher,idim2))
                maxMom      = max(maxMom,absMom)
                maxMomBatch = max(maxMomBatch,absMom)
              ENDDO
            ENDDO
            IF (maxMom.GT.THRESHOLD) THEN
              ICOUNTIND = ICOUNTIND + 1
              INDDSTR(iA,iB,iPassP) = ICOUNTIND
              IF (dotriangular) INDDSTR(iB,iA,iPassP) = ICOUNTIND
            ENDIF
          ENDDO
        ENDDO
        IF (maxMomBatch.GT.THRESHOLD) THEN
          NBATCH = NBATCH + 1
          IBATCHCont(sA:eA,sB:eB,iPassP) = NBATCH
          IF (dotriangular) IBATCHCont(sB:eB,sA:eA,iPassP) = NBATCH
        ENDIF
        sB = eB+1
      ENDDO !icB
      sA = eA+1
    ENDDO !icA
  ENDDO !iAngmom
ENDDO !iPassP
INPUT%MMunique_ID1 = ICOUNTIND
INPUT%MMunique_ID2 = NBATCH

DO iPassP = 1, P%nPasses

  atoma  = P%orb1atom(iPassP)
  batchA = P%orb1batch(iPassP)
  atomb  = P%orb2atom(iPassP)
  batchB = P%orb2batch(iPassP)
  dotriangular  = INPUT%SameLHSaos .AND. ((batchA.EQ.batchB).AND.(atoma.EQ.atomb))
    
  IF (geoderivorder.GT. 0) THEN
    derivInfo%Atom(1)=P%orb1mol(iPassP)
    derivInfo%Atom(2)=P%orb2mol(iPassP)
  ELSE
    iDer  = 0
    iAtom = -1
    atomA = -1
    atomB = -1
  ENDIF
  !
  DO iDerivP=1,ngeoderivcompP
   IF (geoderivorder.GT. 0) THEN
     iAtom = derivInfo%Atom(derivInfo%AO(1,iDerivP))
     iDer  = derivInfo%dirComp(1,iDerivP)
   ENDIF
   DO iB=1,nB
    endA = nA
    IF (dotriangular) endA=iB
    DO iA=1,endA

     ICOUNTIND = INDDSTR(iA,iB,iPassP)
     NBATCH   = IBATCHCont(iA,iB,iPassP)

     call printMMtoFile1(SpherMom,nSpherMom,mmorder,nA,nB,iA,iB,angA(iA),angB(iB),&
           &orbA(iA,iPassP),orbB(iB,iPassP),  &
           &OUTPUT2,PQ,INPUT,iDer,iDerivP,iPassP,ngeoderivcomp,nPass,iAtom,atomA,atomB,&
           &NBATCH,ICOUNTIND,THRESHOLD)
    ENDDO
   ENDDO
  ENDDO

ENDDO
call mem_dealloc(INDDSTR)
call mem_dealloc(IBATCHCont)

!
IF(geoderivorder.gt. 0) CALL freeDerivativeOverlapInfo(derivInfo)
call mem_dealloc(SpherMom)

END SUBROUTINE printMMtoFile

!> \brief wrapper routine for the printing of multipole moments
!> \author \latexonly A. Krapp \endlatexonly
!> \date 2010
 SUBROUTINE printMMtoFile1(SpherMom,nSPHMAT,mmorder,nA,nB,iA,iB,angA,angB,orbA,orbB,&
                        & OUTPUT2,PQ,INPUT,iDer,iDerivP,iPassP,ngeoderivcompP,nPassP,iAtom,& 
                        & atomA,atomB,NBATCH,ICOUNTIND,THRESHOLD)
  implicit none
  Type(integrand),intent(in)         :: PQ
  Type(IntegralInput),intent(in)     :: Input
  Type(IntegralOutput),intent(inout) :: Output2
  !
  Integer,intent(in)    :: mmorder,nA,nB,iA,iB,angA,angB,orbA,orbB,nSPHMAT
  integer,intent(in)    :: iDer,iDerivP,iPassP,ngeoderivcompP,nPassP,iAtom,atomA,atomB
  integer,intent(inout) :: NBATCH,ICOUNTIND
  !
  REAL(REALK),intent(in) :: SpherMom(nSPHMAT,nA,nB,ngeoderivcompP,nPassP),THRESHOLD
  !
  ! local variables
  integer  :: ml,LU1,LU2,imat,m,imm

  IF(INPUT%DO_GRADIENT) THEN
     LU1 = INPUT%lu_mmdada
     LU2 = INPUT%lu_mmdadr
  ELSE
     LU1 = INPUT%lu_mmdata
     LU2 = INPUT%lu_mmdatr
  END IF

  IF (NBATCH.GT. 0 .AND. ICOUNTIND.GT. 0) THEN
    imat = 0
    DO m=0,mmorder
     DO imm=0,m
      DO ml = imm,-imm,-max(2*imm,1)
       imat = imat + 1
       ! empty buffer if necessary
       if (OUTPUT2%USEBUFMM .and. (OUTPUT2%IBUFI.GE. OUTPUT2%MMBUFLEN-1)) then
         CALL LS_EMPTYIBUF(OUTPUT2,OUTPUT2%IBUF,LU1)
         CALL LS_EMPTYRBUF(OUTPUT2,OUTPUT2%RBUF,LU2)
       endif
       call printMMtoFile2(-SpherMom(imat,iA,iB,iDerivP,iPassP),&
                  &nA,nB,iA,iB,angA,angB,orbA,orbB,      &
                  &OUTPUT2,m,ml,PQ,INPUT,&
                  &iDerivP,atomA,atomB,NBATCH,ICOUNTIND,THRESHOLD)
      ENDDO
     ENDDO
    ENDDO
  ENDIF
 END SUBROUTINE printMMtoFile1

!> \brief routine for the printing of multipole moments to file
!> \author \latexonly A. Krapp \endlatexonly
!> \date 2010
 SUBROUTINE printMMtoFile2(integral,nA,nB,iA,iB,angA,angB,orbA,orbB,&
                        & OUTPUT2,m,ml,PQ,INPUT,iDerivP,atomA,atomB,NBATCH,& 
                        & ICOUNTIND,THRESHOLD)
  implicit none
  Type(integrand),intent(in)         :: PQ
  Type(IntegralInput),intent(in)     :: Input
  Type(IntegralOutput),intent(inout) :: Output2
  !
  Integer,intent(in)    :: nA,nB,iA,iB,angA,angB,orbA,orbB
  integer,intent(in)    :: m,ml,iDerivP,atomA,atomB
  integer,intent(inout) :: NBATCH,ICOUNTIND
  !
  REAL(REALK),intent(in) :: integral,THRESHOLD
  !
  ! local variables
  integer      :: iscoor, LU1
  logical      :: dograd
  !
  dograd = INPUT%DO_GRADIENT
  !
  IF(DOGRAD) THEN
     iscoor = iDerivP
     LU1 = INPUT%lu_mmdada
  ELSE
     iscoor = -1
     LU1 = INPUT%lu_mmdata
  END IF

  IF (abs(integral).GT.THRESHOLD) THEN
    if (OUTPUT2%USEBUFMM)then
       CALL LS_FILLBUFFER(OUTPUT2,m,ml,orbB,orbA,angB,angA, &
          &NBATCH,ICOUNTIND,&
          &PQ%P%p%ODextent,PQ%P%p%ODcenter(1),PQ%P%p%ODcenter(2),PQ%P%p%ODcenter(3),&
          &integral,atomA,atomB,iscoor)
    else
       WRITE(LU1)m,ml,orbB,orbA,angB,angA, &
          &NBATCH,ICOUNTIND,&
          &PQ%P%p%ODextent,PQ%P%p%ODcenter(1),PQ%P%p%ODcenter(2),PQ%P%p%ODcenter(3),&
          &integral,atomA,atomB,iscoor
    endif
  ENDIF
END SUBROUTINE printMMtoFile2

!> \brief Distribute Coulomb-type integrals using explicit integrals (as opposed to Jengine type integrals)
!> \author \latexonly S. Reine  \endlatexonly
!> \date 2011-03-08
!> \param Jmat contains the result lstensor
!> \param CDAB matrix containing calculated integrals
!> \param Dmat contains the rhs-density matrix/matrices
!> \param ndmat number of rhs-density matries
!> \param PQ information about the calculated integrals (to be distributed)
!> \param Input contain info about the requested integral 
!> \param Lsoutput contain info about the requested output, and actual output
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
! call distributeJcont(&
!               & integral%integralsABCD,input%LST_DLHS,input%LST_DRHS,&
!               & input%NDMAT_RHS,PQ,Input,output,Integral%Econt,LUPRI,IPRINT)
!         ELSE
!
SUBROUTINE distributeJcont(CDAB,Dlhs,Drhs,ndmat,PQ,Input,Lsoutput,Jcont,LUPRI,IPRINT)
implicit none 
Type(integrand),intent(in)         :: PQ
Type(IntegralInput),intent(in)     :: Input
Type(IntegralOutput),intent(inout) :: Lsoutput
Integer,intent(in)                 :: ndmat,LUPRI,IPRINT
REAL(REALK),pointer                :: CDAB(:)
TYPE(lstensor),intent(in)          :: Dlhs,Drhs
REAL(REALK),intent(inout)          :: Jcont(ndmat)
!
type(overlap),pointer :: P,Q
logical :: SamePQ,SameLHSaos,SameRHSaos,SameODs
Logical :: permuteAB,permuteCD,permuteOD
integer :: nAngmomP,nPasses,iPass,iPassP,iPassQ
integer :: nA,nB,nC,nD,batchA,batchB,batchC,batchD,atomA,atomB,atomC,atomD
integer :: indDCD,indDDC,indDAB,indDBA,maxbat,maxang
integer :: nCd,nDd,sCd,sDd,nAd,nBd,sAd,sBd
Real(realk) :: factor
Real(realk),pointer :: dCD(:),dDC(:),dAB(:),dBA(:)
#ifdef VAR_MPI
call lsquit('not implemented',-1)
#else
!Special for same lhs and rhs !!!!
!assumes same LHS and same RHS and symmetric matrices
!THIS DOES NOT HOLD FOR MPI....
SameLHSaos     = INPUT%SameLHSaos
SameRHSaos     = INPUT%SameRHSaos
SameODs        = INPUT%SameODs
IF(PQ%reverseOrder)CALL LSQUIT('PQ%reverseOrder Jcont',lupri) 
P => PQ%P%p
Q => PQ%Q%p
nAngmomP       = P%nAngmom
SamePQ         = PQ%samePQ
nPasses = P%nPasses*Q%nPasses
iPass = 0
nA     = P%orbital1%totOrbitals
nB     = P%orbital2%totOrbitals
nC     = Q%orbital1%totOrbitals
nD     = Q%orbital2%totOrbitals
DO iPassP=1,P%nPasses
  atomA  = P%orb1atom(iPassP)
  batchA = P%orb1batch(iPassP)
  atomB  = P%orb2atom(iPassP)
  batchB = P%orb2batch(iPassP)
  permuteAB  = SameLHSaos.AND.((batchA.NE.batchB).OR.(atomA.NE.atomB))
  indDAB = Dlhs%index(atomA,atomB,1,1)
  dAB => Dlhs%LSAO(indDAB)%elms

  nAd = Dlhs%LSAO(indDAB)%nLocal(1)
  nBd = Dlhs%LSAO(indDAB)%nLocal(2)
  maxBat = Dlhs%LSAO(indDAB)%maxBat
  maxAng = Dlhs%LSAO(indDAB)%maxAng
  sAd = Dlhs%LSAO(indDAB)%startLocalOrb(1+(batchA-1)*maxAng)-1 
  sBd = Dlhs%LSAO(indDAB)%startLocalOrb(1+(batchB-1)*maxAng+maxAng*maxBat)-1

  !can contract the CDAB integral at this point using a dgemm giving
  ! iC,iD,idmat,ipass
  
  DO iPassQ=1,Q%nPasses
    iPass = iPass + 1

    atomC  = Q%orb1atom(iPassQ)
    batchC = Q%orb1batch(iPassQ)
    atomD  = Q%orb2atom(iPassQ)
    batchD = Q%orb2batch(iPassQ)

    permuteCD = SameRHSaos.AND.((batchC.NE.batchD).OR.(atomC.NE.atomD))
    permuteOD = SameODs.AND.( ((batchA.NE.batchC).OR.(atomA.NE.atomC)).OR.&
   &                          ((batchB.NE.batchD).OR.(atomB.NE.atomD)) )

    indDCD = Dlhs%index(atomC,atomD,1,1)
    dCD => Dlhs%LSAO(indDCD)%elms

    nCd = Dlhs%LSAO(indDCD)%nLocal(1)
    nDd = Dlhs%LSAO(indDCD)%nLocal(2)

    maxBat = Dlhs%LSAO(indDCD)%maxBat
    maxAng = Dlhs%LSAO(indDCD)%maxAng
    sCd = Dlhs%LSAO(indDCD)%startLocalOrb(1+(batchC-1)*maxAng)-1 
    sDd = Dlhs%LSAO(indDCD)%startLocalOrb(1+(batchD-1)*maxAng+maxAng*maxBat)-1

    factor = 2.0E0_realk
    IF(permuteAB)factor = 4.0E0_realk
    IF(permuteCD)factor = 2.0E0_realk*factor
    IF(permuteOD)factor = 2.0E0_realk*factor
    CALL JCONT_ABCD(Jcont,dAB,dCD,nAd,nBd,nCd,nDd,sAd,sBd,sCd,sDd,CDAB,nA,nB,nC,nD,&
         & iPass,nPasses,ndmat,factor,lupri)
  ENDDO !iPassQ
ENDDO !iPassP
#endif

CONTAINS
SUBROUTINE JCONT_ABCD(Jcont,dAB,dCD,nAd,nBd,nCd,nDd,sAd,sBd,sCd,sDd,CDAB,&
     &nA,nB,nC,nD,iPass,nPasses,ndmat,factor,lupri)
implicit none
Integer,intent(IN)        :: nA,nB,nC,nD,iPass,nPasses,lupri,ndmat
Integer,intent(IN)        :: nAd,nBd,nCd,nDd,sAd,sBd,sCd,sDd
Real(realk),intent(IN)    :: CDAB(nC,nD,nA,nB,nPasses),factor
Real(realk),intent(IN)    :: dAB(nAd,nBd,ndmat)
Real(realk),intent(IN)    :: dCD(nCd,nDd,ndmat)
Real(realk),intent(INOUT) :: Jcont(ndmat)
Real(realk)    :: dABTMP(nA,nB)
Real(realk)    :: dCDTMP(nC,nD)
!
Integer :: iA,iB,iC,iD,idmat
Real(realk) :: tmp,tmpD
!
!DO iB=1,nB
! DO iA=1,nA
!  DO idmat=1,ndmat
!   TMPD = DAB(sAd+iA,sBd+iB,idmat)
!   TMP = 0.0E0_realk
!   DO iD=1,nD
!    DO iC=1,nC
!     TMP = TMP  + TMPD * CDAB(iC,iD,iA,iB,iPass) * dCD(sCd+iC,sDd+iD,idmat)
!    ENDDO
!   ENDDO
!   Jcont(idmat) = Jcont(idmat) + TMP*factor
!  ENDDO
! ENDDO
!ENDDO
!this can be avoided if sA = 0 and nAd = nA ...
DO idmat=1,ndmat
 DO iB=1,nB
  DO iA=1,nA
    DABTMP(iA,iB) = DAB(sAd+iA,sBd+iB,idmat)
  ENDDO
 ENDDO
 DO iD=1,nD
  DO iC=1,nC
    DCDTMP(iC,iD) = dCD(sCd+iC,sDd+iD,idmat)
  ENDDO
 ENDDO
 CALL JCONT_ABCD2(CDAB,nC*nD,nA*nB,iPass,nPasses,dABTMP,dCDTMP,TMP) 
 Jcont(idmat) = Jcont(idmat) + TMP*factor
ENDDO
END SUBROUTINE JCONT_ABCD

SUBROUTINE JCONT_ABCD2(CDAB,nCnD,nAnB,iPass,nPasses,dABTMP,dCDTMP,TMPJCONT)
implicit none
Integer,intent(IN)        :: nAnB,nCnD,iPass,nPasses
Real(realk),intent(IN)    :: CDAB(nCnD,nAnB,nPasses)
Real(realk),intent(IN)    :: dABTMP(nAnB)
Real(realk),intent(IN)    :: dCDTMP(nCnD)
Real(realk),intent(INOUT) :: TMPJcont
!
Real(realk) :: TMPD
Integer :: iAiB,iCiD
TMPJcont = 0.0d0
DO iAiB=1,nAnB
 TMPD = dABTMP(iAiB)
 DO iCiD=1,nCnD
  TMPJcont = TMPJcont + TMPD * CDAB(iCiD,iAiB,iPass) * dCDTMP(iCiD)
 ENDDO
ENDDO
END SUBROUTINE JCONT_ABCD2

END SUBROUTINE distributeJcont

SUBROUTINE getInputDerivativeInfo(derOrder,startDer,endDer,input,emptyP,emptyQ)
implicit none
Integer,intent(OUT)            :: derOrder,startDer,endDer
TYPE(integralinput),intent(IN) :: input
Logical,intent(IN)             :: emptyP,emptyQ

  derOrder = input%GEODERIVORDER
  startDer = 0
  endDer   = derOrder
  IF ((INPUT%geoderOrderP.EQ.0).OR.emptyP) endDer = 0           !Ie. iDerLHS = 0,        iDerRHS = derOrder
  IF ((INPUT%geoderOrderQ.EQ.0).OR.emptyQ) startDer = derOrder  !Ie. iDerLHS = derOrder, iDerRHS = 0

END SUBROUTINE getInputDerivativeInfo

END MODULE thermite_distributeGen


