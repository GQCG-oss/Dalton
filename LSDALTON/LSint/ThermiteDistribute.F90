!> @file 
!> Contains integral distribution routines that places calculated integrals in the proper output
!> Thermite integral distribution module
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
MODULE thermite_distribute
  use thermite_distributeGen
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

CONTAINS
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
SUBROUTINE distributeCoulomb(Jmat,CDABin,Dmat,ndmat,PQ,Input,Lsoutput,LUPRI,IPRINT)
implicit none 
Type(integrand),intent(in)         :: PQ
Type(IntegralInput),intent(in)     :: Input
Type(IntegralOutput),intent(inout) :: Lsoutput
Integer,intent(in)                 :: ndmat,LUPRI,IPRINT
REAL(REALK),pointer                :: CDABin(:)
!REAL(REALK),intent(inout)          :: Jcont(Input%NDMAT_RHS)
TYPE(lstensor),intent(inout)       :: Jmat
TYPE(lstensor),intent(in)          :: Dmat
!
type(overlap),pointer :: P,Q
REAL(REALK),pointer                :: CDAB(:)
logical :: SamePQ,SameLHSaos,SameRHSaos,SameODs,nucleiP,nucleiQ
Logical :: permuteAB,permuteCD,permuteOD,GederivCoulomb,AntiCD
Logical :: antipermuteCD
integer :: nAngmomP,ngeoderivcompQ,ngeoderivcompP,ngeoderivcomp,nPasses
integer :: iPass,iPassP,iPassQ,iDeriv,iDerivP,iDerivQ
integer :: nA,nB,nC,nD,batchA,batchB,batchC,batchD,atomA,atomB,atomC,atomD
integer :: indJAB,indJBA,indJCD,indJDC,indDCD,indDDC,indDAB,indDBA,maxbat,maxang
integer :: nAj,nBj,nCd,nDd,sAj,sBj,sCd,sDd,nCj,nDj,nAd,nBd,sCj,sDj,sAd,sBd
Real(realk),pointer :: jAB(:),jBA(:),jCD(:),jDC(:),dCD(:),dDC(:),dAB(:),dBA(:)
Real(realk)  :: factor
Type(derivativeInfo) :: derivInfo
integer :: ndim5Q,ndim5P,ndim5,nderiv,ider,iatom,ndim5output,idim5
logical :: SameActualODs
IF(INPUT%fullcontraction)call lsquit('Jcont not implemented without Jengine',-1)
antiCD=.FALSE.
SameLHSaos     = INPUT%SameLHSaos
SameRHSaos     = INPUT%SameRHSaos
SameODs        = INPUT%SameODs
ngeoderivcompQ        = 1
ngeoderivcompP        = Input%ngeoderivcomp
ngeoderivcomp         = ngeoderivcompP*ngeoderivcompQ
ndim5Q = ngeoderivcompQ
ndim5P = ngeoderivcompP
IF(PQ%reverseOrder)THEN
   IF(INPUT%magderOrderQ.NE.1)THEN
      !so far this is the only option when PQ%reverseOrder
      CALL LSQUIT('PQ%reverseOrder.AND.INPUT%magderOrderQ.NE.1',lupri) 
   ENDIF
   antiCD=.TRUE.
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
   call mem_workpointer_alloc(CDAB,nA*nB*nC*nD*3*nPasses)
   call transposeQPorder(CDABin,CDAB,nA*nB,nC*nD,3,nPasses,lupri)
ELSE
   P => PQ%P%p
   Q => PQ%Q%p
   CDAB => CDABin
ENDIF

nAngmomP       = P%nAngmom
SamePQ         = PQ%samePQ
nderiv = ngeoderivcompP*ngeoderivcompQ 
ndim5output = Lsoutput%ndim(5)/ndmat
GederivCoulomb = INPUT%geoderOrderP .EQ. 1
IF(GederivCoulomb)then
   IF(SameODs)call lsquit('Error in distributeCoulomb sameOD derivCoulo',-1)
ENDIF
IF(nderiv.GT.1)then
   CALL initDerivativeOverlapInfo(derivInfo,PQ,Input)
   IF (IPRINT.GT. 50) THEN
      CALL printDerivativeOverlapInfo(derivInfo,6)
   ENDIF
endif
ndim5 = ndim5P*ndim5Q
nucleiP        = (P%orbital1%type_nucleus.AND.P%orbital2%type_empty).OR. &
     &           (P%orbital1%type_empty.AND.P%orbital2%type_nucleus)
nucleiQ        = (Q%orbital1%type_nucleus.AND.Q%orbital2%type_empty).OR. &
     &           (Q%orbital1%type_empty.AND.Q%orbital2%type_nucleus)

factor = 1E0_realk

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
  IF(nderiv.GT.1)then
     derivInfo%Atom(1)=P%orb1mol(iPassP)
     derivInfo%Atom(2)=P%orb2mol(iPassP)
  ENDIF
  IF (nucleiP) THEN
    atomA = 1
    atomB = 1
  ENDIF
  permuteAB  = SameLHSaos.AND.((batchA.NE.batchB).OR.(atomA.NE.atomB))
  DO iPassQ=1,Q%nPasses
    iPass = iPass + 1

    atomC  = Q%orb1atom(iPassQ)
    batchC = Q%orb1batch(iPassQ)
    atomD  = Q%orb2atom(iPassQ)
    batchD = Q%orb2batch(iPassQ)
    IF(nderiv.GT.1)then
       derivInfo%Atom(3)=Q%orb1mol(iPassQ)
       derivInfo%Atom(4)=Q%orb2mol(iPassQ)
    ENDIF
    IF (nucleiQ) THEN
      atomC = 1
      atomD = 1
    ENDIF
    permuteCD = SameRHSaos.AND.((batchC.NE.batchD).OR.(atomC.NE.atomD))
    permuteOD = SameODs.AND.( ((batchA.NE.batchC).OR.(atomA.NE.atomC)).OR.&
   &                          ((batchB.NE.batchD).OR.(atomB.NE.atomD)) )
    AntipermuteCD  = permuteCD.AND.antiCD
    IF(GederivCoulomb)then
       !For the calculation of the geometical derivative of the Coulomb
       !matrix 
       !J^x_{AB} = ((AB)^x|CD)D_{CD} + (AB|(CD)^x)D_{CD} 
       !can be calculated as 
       !J^x_{AB} = ((AB)^x|CD)D_{CD} + D_{CD}((CD)^x|AB)
       !so we set permuteOD to true 
       !WARNING: This only works for non MPI when all 4 AOs are the same. 
#ifdef VAR_MPI
       print*,'calc of J^x_{AB} uses a hack that do not work for MPI'
       print*,'This will be fixed in the next release'
       call lsquit('calc of J^x_{AB} uses a hack that do not work for MPI',-1)
#endif       
       permuteOD = .TRUE.
    ENDIF
    indJAB = Jmat%index(atomA,atomB,1,1)
    indDCD = Dmat%index(atomC,atomD,1,1)
    jAB => Jmat%LSAO(indJAB)%elms
    dCD => Dmat%LSAO(indDCD)%elms

    nAj = Jmat%LSAO(indJAB)%nLocal(1)
    nBj = Jmat%LSAO(indJAB)%nLocal(2)
    nCd = Dmat%LSAO(indDCD)%nLocal(1)
    nDd = Dmat%LSAO(indDCD)%nLocal(2)

    maxBat = Jmat%LSAO(indJAB)%maxBat
    maxAng = Jmat%LSAO(indJAB)%maxAng
    sAj = Jmat%LSAO(indJAB)%startLocalOrb(1+(batchA-1)*maxAng)-1 
    sBj = Jmat%LSAO(indJAB)%startLocalOrb(1+(batchB-1)*maxAng+maxAng*maxBat)-1
    maxBat = Dmat%LSAO(indDCD)%maxBat
    maxAng = Dmat%LSAO(indDCD)%maxAng
    sCd = Dmat%LSAO(indDCD)%startLocalOrb(1+(batchC-1)*maxAng)-1 
    sDd = Dmat%LSAO(indDCD)%startLocalOrb(1+(batchD-1)*maxAng+maxAng*maxBat)-1

    IF (permuteAB) THEN
      indJBA = Jmat%index(atomB,atomA,1,1)
      jBA => Jmat%LSAO(indJBA)%elms
    ENDIF
    IF (permuteCD) THEN
      indDDC = Dmat%index(atomD,atomC,1,1)
      dDC => Dmat%LSAO(indDDC)%elms
    ENDIF
    IF (permuteOD) THEN
      indJCD = Jmat%index(atomC,atomD,1,1)
      indDAB = Dmat%index(atomA,atomB,1,1)
      jCD => Jmat%LSAO(indJCD)%elms
      dAB => Dmat%LSAO(indDAB)%elms
      nCj = Jmat%LSAO(indJCD)%nLocal(1)
      nDj = Jmat%LSAO(indJCD)%nLocal(2)
      nAd = Dmat%LSAO(indDAB)%nLocal(1)
      nBd = Dmat%LSAO(indDAB)%nLocal(2)

      maxBat = Jmat%LSAO(indJCD)%maxBat
      maxAng = Jmat%LSAO(indJCD)%maxAng
      sCj = Jmat%LSAO(indJCD)%startLocalOrb(1+(batchC-1)*maxAng) - 1 
      sDj = Jmat%LSAO(indJCD)%startLocalOrb(1+(batchD-1)*maxAng+maxAng*maxBat)-1
      maxBat = Dmat%LSAO(indDAB)%maxBat
      maxAng = Dmat%LSAO(indDAB)%maxAng
      sAd = Dmat%LSAO(indDAB)%startLocalOrb(1+(batchA-1)*maxAng) - 1 
      sBd = Dmat%LSAO(indDAB)%startLocalOrb(1+(batchB-1)*maxAng+maxAng*maxBat)-1

      IF (permuteAB) THEN
        indDBA = Dmat%index(atomB,atomA,1,1)
        dBA => Dmat%LSAO(indDBA)%elms
      ENDIF
      IF (permuteCD) THEN
        indJDC = Jmat%index(atomD,atomC,1,1)
        jDC => Jmat%LSAO(indJDC)%elms
      ENDIF
    ENDIF
    iDeriv = 0
    DO iDerivP=1,ndim5P!ngeoderivcompP
     DO iDerivQ=1,ndim5Q!ngeoderivcompQ
      iDeriv = iDeriv + 1
      iDim5=ideriv
      IF (nderiv.GT. 1)THEN
         iDer = derivInfo%dirComp(1,iDerivP)
         iAtom = derivInfo%Atom(derivInfo%AO(1,iDerivP))
         iDim5 = 3*(iAtom-1)+ider
      ENDIF
!$OMP CRITICAL (distributeFockBlock)
      IF (permuteOD.AND.permuteCD.AND.permuteAB) THEN
        CALL jAB_CD(jAB,dCD,nAj,nBj,nCd,nDd,sAj,sBj,sCd,sDd,CDAB,factor,nA,nB,nC,nD,&
             & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
        CALL jBA_CD(jBA,dCD,nBj,nAj,nCd,nDd,sBj,sAj,sCd,sDd,CDAB,factor,nA,nB,nC,nD,&
             & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
        CALL jAB_DC(jAB,dDC,nAj,nBj,nDd,nCd,sAj,sBj,sDd,sCd,CDAB,factor,nA,nB,nC,nD,&
             & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
        CALL jBA_DC(jBA,dDC,nBj,nAj,nDd,nCd,sBj,sAj,sDd,sCd,CDAB,factor,nA,nB,nC,nD,&
             & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
        CALL jCD_AB(jCD,dAB,nCj,nDj,nAd,nBd,sCj,sDj,sAd,sBd,CDAB,factor,nA,nB,nC,nD,&
             & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
        CALL jDC_AB(jDC,dAB,nDj,nCj,nAd,nBd,sDj,sCj,sAd,sBd,CDAB,factor,nA,nB,nC,nD,&
             & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
        CALL jCD_BA(jCD,dBA,nCj,nDj,nBd,nAd,sCj,sDj,sBd,sAd,CDAB,factor,nA,nB,nC,nD,&
             & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
        CALL jDC_BA(jDC,dBA,nDj,nCj,nBd,nAd,sDj,sCj,sBd,sAd,CDAB,factor,nA,nB,nC,nD,&
             & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
      ELSE IF (permuteOD.AND.permuteCD) THEN
        CALL jAB_CD(jAB,dCD,nAj,nBj,nCd,nDd,sAj,sBj,sCd,sDd,CDAB,factor,nA,nB,nC,nD,&
             & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
        CALL jAB_DC(jAB,dDC,nAj,nBj,nDd,nCd,sAj,sBj,sDd,sCd,CDAB,factor,nA,nB,nC,nD,&
             & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
        CALL jCD_AB(jCD,dAB,nCj,nDj,nAd,nBd,sCj,sDj,sAd,sBd,CDAB,factor,nA,nB,nC,nD,&
             & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
        CALL jDC_AB(jDC,dAB,nDj,nCj,nAd,nBd,sDj,sCj,sAd,sBd,CDAB,factor,nA,nB,nC,nD,&
             & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
      ELSE IF (permuteOD.AND.permuteAB) THEN
        CALL jAB_CD(jAB,dCD,nAj,nBj,nCd,nDd,sAj,sBj,sCd,sDd,CDAB,factor,nA,nB,nC,nD,&
             & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
        CALL jBA_CD(jBA,dCD,nBj,nAj,nCd,nDd,sBj,sAj,sCd,sDd,CDAB,factor,nA,nB,nC,nD,&
             & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
        CALL jCD_AB(jCD,dAB,nCj,nDj,nAd,nBd,sCj,sDj,sAd,sBd,CDAB,factor,nA,nB,nC,nD,&
             & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
        CALL jCD_BA(jCD,dBA,nCj,nDj,nBd,nAd,sCj,sDj,sBd,sAd,CDAB,factor,nA,nB,nC,nD,&
             & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
      ELSE IF (permuteAB.AND.AntipermuteCD) THEN
         CALL jAB_CD(jAB,dCD,nAj,nBj,nCd,nDd,sAj,sBj,sCd,sDd,CDAB,factor,nA,nB,nC,nD,&
              & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
         CALL jBA_CD(jBA,dCD,nBj,nAj,nCd,nDd,sBj,sAj,sCd,sDd,CDAB,factor,nA,nB,nC,nD,&
              & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
         CALL jAB_DC(jAB,dDC,nAj,nBj,nDd,nCd,sAj,sBj,sDd,sCd,CDAB,-factor,nA,nB,nC,nD,&
              & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
         CALL jBA_DC(jBA,dDC,nBj,nAj,nDd,nCd,sBj,sAj,sDd,sCd,CDAB,-factor,nA,nB,nC,nD,&
              & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
      ELSE IF (permuteAB.AND.permuteCD) THEN
         CALL jAB_CD(jAB,dCD,nAj,nBj,nCd,nDd,sAj,sBj,sCd,sDd,CDAB,factor,nA,nB,nC,nD,&
              & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
         CALL jBA_CD(jBA,dCD,nBj,nAj,nCd,nDd,sBj,sAj,sCd,sDd,CDAB,factor,nA,nB,nC,nD,&
              & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
         CALL jAB_DC(jAB,dDC,nAj,nBj,nDd,nCd,sAj,sBj,sDd,sCd,CDAB,factor,nA,nB,nC,nD,&
              & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
         CALL jBA_DC(jBA,dDC,nBj,nAj,nDd,nCd,sBj,sAj,sDd,sCd,CDAB,factor,nA,nB,nC,nD,&
              & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
      ELSE IF (permuteAB) THEN
        CALL jAB_CD(jAB,dCD,nAj,nBj,nCd,nDd,sAj,sBj,sCd,sDd,CDAB,factor,nA,nB,nC,nD,&
             & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
        CALL jBA_CD(jBA,dCD,nBj,nAj,nCd,nDd,sBj,sAj,sCd,sDd,CDAB,factor,nA,nB,nC,nD,&
             & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
      ELSE IF (AntipermuteCD) THEN
         CALL jAB_CD(jAB,dCD,nAj,nBj,nCd,nDd,sAj,sBj,sCd,sDd,CDAB,factor,nA,nB,nC,nD,&
              & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
         CALL jAB_DC(jAB,dDC,nAj,nBj,nDd,nCd,sAj,sBj,sDd,sCd,CDAB,-factor,nA,nB,nC,nD,&
              & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
      ELSE IF (permuteCD) THEN
         CALL jAB_CD(jAB,dCD,nAj,nBj,nCd,nDd,sAj,sBj,sCd,sDd,CDAB,factor,nA,nB,nC,nD,&
              & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
         CALL jAB_DC(jAB,dDC,nAj,nBj,nDd,nCd,sAj,sBj,sDd,sCd,CDAB,factor,nA,nB,nC,nD,&
              & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
      ELSE IF (permuteOD) THEN
         CALL jAB_CD(jAB,dCD,nAj,nBj,nCd,nDd,sAj,sBj,sCd,sDd,CDAB,factor,nA,nB,nC,nD,&
              & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
         CALL jCD_AB(jCD,dAB,nCj,nDj,nAd,nBd,sCj,sDj,sAd,sBd,CDAB,factor,nA,nB,nC,nD,&
              & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
      ELSE
         CALL jAB_CD(jAB,dCD,nAj,nBj,nCd,nDd,sAj,sBj,sCd,sDd,CDAB,factor,nA,nB,nC,nD,&
              & iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
      ENDIF
!$OMP END CRITICAL (distributeFockBlock)
     ENDDO !iDerivQ
    ENDDO !iDerivP
  ENDDO !iPassQ
ENDDO !iPassP
IF (nderiv.GT. 1) CALL freeDerivativeOverlapInfo(derivInfo)
IF(PQ%reverseOrder)THEN
   call mem_workpointer_dealloc(CDAB)
ENDIF

CONTAINS
SUBROUTINE jAB_CD(jAB,dCD,nAj,nBj,nCd,nDd,sAj,sBj,sCd,sDd,CDAB,factor,&
     &nA,nB,nC,nD,iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
implicit none
Integer,intent(IN)        :: nA,nB,nC,nD,iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,lupri,ndmat
Integer,intent(IN)        :: nAj,nBj,nCd,nDd,sAj,sBj,sCd,sDd
Real(realk),intent(IN)    :: CDAB(nC,nD,nA,nB,ndim5,nPasses)
Real(realk),intent(INOUT) :: jAB(nAj,nBj,ndim5output,ndmat)!(nP,ndim5,ndmat)
Real(realk),intent(IN)    :: dCD(nCd,nDd,ndmat)
Real(realk),intent(IN)    :: factor
!
Integer :: iA,iB,iC,iD,idmat
!
DO idmat=1,ndmat
 DO iD=1,nD
  DO iC=1,nC
   DO iB=1,nB
    DO iA=1,nA
     jAB(sAj+iA,sBj+iB,idim5,idmat) = jAB(sAj+iA,sBj+iB,idim5,idmat) + &
          &factor * CDAB(iC,iD,iA,iB,iDeriv,iPass) * dCD(sCd+iC,sDd+iD,idmat)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE jAB_CD

SUBROUTINE jCD_AB(jCD,dAB,nCj,nDj,nAd,nBd,sCj,sDj,sAd,sBd,CDAB,factor,&
     &nA,nB,nC,nD,iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
implicit none
Integer,intent(IN)        :: nA,nB,nC,nD,iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,lupri,ndmat
Integer,intent(IN)        :: nCj,nDj,nAd,nBd,sCj,sDj,sAd,sBd
Real(realk),intent(IN)    :: CDAB(nC,nD,nA,nB,ndim5,nPasses)
Real(realk),intent(INOUT) :: jCD(nCj,nDj,ndim5output,ndmat)
Real(realk),intent(IN)    :: dAB(nAd,nBd,ndmat)
Real(realk),intent(IN)    :: factor
!
Integer :: iA,iB,iC,iD,idmat
!
DO idmat=1,ndmat
 DO iD=1,nD
  DO iC=1,nC
   DO iB=1,nB
    DO iA=1,nA
     jCD(sCj+iC,sDj+iD,idim5,idmat) = jCD(sCj+iC,sDj+iD,idim5,idmat) + &
          & factor * CDAB(iC,iD,iA,iB,iDeriv,iPass) * dAB(sAd+iA,sBd+iB,idmat)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE jCD_AB

SUBROUTINE jAB_DC(jAB,dDC,nAj,nBj,nDd,nCd,sAj,sBj,sDd,sCd,CDAB,factor,&
     &nA,nB,nC,nD,iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
implicit none
Integer,intent(IN)        :: nA,nB,nC,nD,iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,lupri,ndmat
Integer,intent(IN)        :: nAj,nBj,nDd,nCd,sAj,sBj,sDd,sCd
Real(realk),intent(IN)    :: CDAB(nC,nD,nA,nB,ndim5,nPasses)
Real(realk),intent(INOUT) :: jAB(nAj,nBj,ndim5output,ndmat)
Real(realk),intent(IN)    :: dDC(nDd,nCd,ndmat)
Real(realk),intent(IN)    :: factor
!
Integer :: iA,iB,iC,iD,idmat
!
DO idmat=1,ndmat
 DO iB=1,nB
  DO iA=1,nA
   DO iD=1,nD
    DO iC=1,nC
     jAB(sAj+iA,sBj+iB,idim5,idmat) = jAB(sAj+iA,sBj+iB,idim5,idmat) + &
          & factor * CDAB(iC,iD,iA,iB,iDeriv,iPass) * dDC(sDd+iD,sCd+iC,idmat)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE jAB_DC

SUBROUTINE jDC_AB(jDC,dAB,nDj,nCj,nAd,nBd,sDj,sCj,sAd,sBd,CDAB,factor,&
     & nA,nB,nC,nD,iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
implicit none
Integer,intent(IN)        :: nA,nB,nC,nD,iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,lupri,ndmat
Integer,intent(IN)        :: nDj,nCj,nAd,nBd,sDj,sCj,sAd,sBd
Real(realk),intent(IN)    :: CDAB(nC,nD,nA,nB,ndim5,nPasses)
Real(realk),intent(INOUT) :: jDC(nDj,nCj,ndim5output,ndmat)
Real(realk),intent(IN)    :: dAB(nAd,nBd,ndmat)
Real(realk),intent(IN)    :: factor
!
Integer :: iA,iB,iC,iD,idmat
!
DO idmat=1,ndmat
 DO iD=1,nD
  DO iC=1,nC
   DO iB=1,nB
    DO iA=1,nA
     jDC(sDj+iD,sCj+iC,idim5,idmat) = jDC(sDj+iD,sCj+iC,idim5,idmat) + &
          & factor * CDAB(iC,iD,iA,iB,iDeriv,iPass) * dAB(sAd+iA,sBd+iB,idmat)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE jDC_AB

SUBROUTINE jBA_CD(jBA,dCD,nBj,nAj,nCd,nDd,sBj,sAj,sCd,sDd,CDAB,factor,&
     & nA,nB,nC,nD,iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
implicit none
Integer,intent(IN)        :: nA,nB,nC,nD,iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,lupri,ndmat
Integer,intent(IN)        :: nBj,nAj,nCd,nDd,sBj,sAj,sCd,sDd
Real(realk),intent(IN)    :: CDAB(nC,nD,nA,nB,ndim5,nPasses)
Real(realk),intent(INOUT) :: jBA(nBj,nAj,ndim5output,ndmat)
Real(realk),intent(IN)    :: dCD(nCd,nDd,ndmat)
Real(realk),intent(IN)    :: factor
!
Integer :: iA,iB,iC,iD,idmat
!
DO idmat=1,ndmat
 DO iB=1,nB
  DO iA=1,nA
   DO iD=1,nD
    DO iC=1,nC
     jBA(sBj+iB,sAj+iA,idim5,idmat) = jBA(sBj+iB,sAj+iA,idim5,idmat) &
          & + factor * CDAB(iC,iD,iA,iB,iDeriv,iPass) * dCD(sCd+iC,sDd+iD,idmat)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE jBA_CD

SUBROUTINE jCD_BA(jCD,dBA,nCj,nDj,nBd,nAd,sCj,sDj,sBd,sAd,CDAB,factor,&
     & nA,nB,nC,nD,iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
implicit none
Integer,intent(IN)        :: nA,nB,nC,nD,iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,lupri,ndmat
Integer,intent(IN)        :: nCj,nDj,nBd,nAd,sCj,sDj,sBd,sAd
Real(realk),intent(IN)    :: CDAB(nC,nD,nA,nB,ndim5,nPasses)
Real(realk),intent(INOUT) :: jCD(nCj,nDj,ndim5output,ndmat)
Real(realk),intent(IN)    :: dBA(nBd,nAd,ndmat)
Real(realk),intent(IN)    :: factor
!
Integer :: iA,iB,iC,iD,idmat
!
DO idmat=1,ndmat
 DO iD=1,nD
  DO iC=1,nC
   DO iB=1,nB
    DO iA=1,nA
     jCD(sCj+iC,sDj+iD,idim5,idmat) = jCD(sCj+iC,sDj+iD,idim5,idmat) + &
          & factor * CDAB(iC,iD,iA,iB,iDeriv,iPass) * dBA(sBd+iB,sAd+iA,idmat)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE jCD_BA

SUBROUTINE jBA_DC(jBA,dDC,nBj,nAj,nDd,nCd,sBj,sAj,sDd,sCd,CDAB,factor,&
     & nA,nB,nC,nD,iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
implicit none
Integer,intent(IN)        :: nA,nB,nC,nD,iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,lupri,ndmat
Integer,intent(IN)        :: nBj,nAj,nDd,nCd,sBj,sAj,sDd,sCd
Real(realk),intent(IN)    :: CDAB(nC,nD,nA,nB,ndim5,nPasses)
Real(realk),intent(INOUT) :: jBA(nBj,nAj,ndim5output,ndmat)
Real(realk),intent(IN)    :: dDC(nDd,nCd,ndmat)
Real(realk),intent(IN)    :: factor
!
Integer :: iA,iB,iC,iD,idmat
!
DO idmat=1,ndmat
 DO iB=1,nB
  DO iA=1,nA
   DO iD=1,nD
    DO iC=1,nC
     jBA(sBj+iB,sAj+iA,idim5,idmat) = jBA(sBj+iB,sAj+iA,idim5,idmat) + &
          & factor * CDAB(iC,iD,iA,iB,iDeriv,iPass) * dDC(sDd+iD,sCd+iC,idmat)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE jBA_DC

SUBROUTINE jDC_BA(jDC,dBA,nDj,nCj,nBd,nAd,sDj,sCj,sBd,sAd,CDAB,factor,&
     & nA,nB,nC,nD,iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,ndmat,lupri)
implicit none
Integer,intent(IN)        :: nA,nB,nC,nD,iPass,iDeriv,idim5,nPasses,ndim5,ndim5output,lupri,ndmat
Integer,intent(IN)        :: nDj,nCj,nBd,nAd,sDj,sCj,sBd,sAd
Real(realk),intent(IN)    :: CDAB(nC,nD,nA,nB,ndim5,nPasses)
Real(realk),intent(INOUT) :: jDC(nDj,nCj,ndim5output,ndmat)
Real(realk),intent(IN)    :: dBA(nBd,nAd,ndmat)
Real(realk),intent(IN)    :: factor
!
Integer :: iA,iB,iC,iD,idmat
!
DO idmat=1,ndmat
 DO iD=1,nD
  DO iC=1,nC
   DO iB=1,nB
    DO iA=1,nA
     jDC(sDj+iD,sCj+iC,idim5,idmat) = jDC(sDj+iD,sCj+iC,idim5,idmat) + &
          & factor * CDAB(iC,iD,iA,iB,iDeriv,iPass) * dBA(sBd+iB,sAd+iA,idmat)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE jDC_BA

END SUBROUTINE distributeCoulomb

END MODULE thermite_distribute

