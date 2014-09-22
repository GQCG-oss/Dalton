!> @file 
!> Contains integral distribution routines that places calculated integrals in the proper output
!> Thermite integral distribution module
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
MODULE thermite_distribute2
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
SUBROUTINE distributeJengine(RES,PQ,QPmat2,dimQ,dimP,Input,Lsoutput,Jcont,LUPRI,IPRINT)
  implicit none 
  Type(integrand),intent(in)      :: PQ
  Type(IntegralInput),intent(in)  :: Input
  Type(IntegralOutput),intent(inout) :: Lsoutput
  Integer,intent(in)              :: LUPRI,IPRINT
  Integer,intent(in)              :: dimQ,dimP
  REAL(REALK),intent(in)          :: QPMAT2(dimQ,dimP)
  REAL(REALK),intent(inout)       :: Jcont(Input%NDMAT_RHS)
  TYPE(lstensor),intent(inout)    :: RES
  TYPE(Overlap),pointer           :: P
  !
  Integer :: nAngmomP,nOrbP,ngeoderivcompP,nMAT
  Integer :: iAngmomP,iOrbitalP,iderivP
  Integer :: iA,iB,sA,sB,sC,sD,SAA,SBB
  Integer :: nContA,nContB,nAngA,nAngB,startA,startB
  Integer :: iOrbP,idmat,iPassP,n1,n2,nA,nB,nC,nD
  logical :: SameLHSaos,dograd,dopermutation
  integer :: AtomA,atomB,geoderivorder,AB,BA,dimA,dimB,iLSAO
  real(realk) :: dummy(Input%NDMAT_RHS)
!
  Type(derivativeInfo) :: derivInfo
  Integer     :: iAtom,iDer,passPoffset,batchA,batchB,localA,localB
  Integer     :: Dabind,Dbaind,idim5,ndim5P,ndim5output,maxBat,maxAng
  Real(realk) :: derCont
  Real(realk),pointer :: JAB(:),JBA(:),Jgrad(:),DAB(:),DBA(:)
  LOGICAL :: antiAB,doAntipermutation,doFullContract
  antiAB=.FALSE.

  dograd = INPUT%DO_GRADIENT
  doFullContract = INPUT%fullcontraction
  nMAT=Input%NDMAT_RHS
  SameLHSaos = INPUT%SameLHSaos
  !Replace with lhsGeoOrder?
  ngeoderivcompP=Input%ngeoderivcomp
  geoderivorder = Input%geoderivorder
  ndim5P=NgeoderivcompP
  ndim5output = Lsoutput%ndim(5)/nMAT
  IF (geoderivorder.GT. 0) THEN
     CALL initDerivativeOverlapInfo(derivInfo,PQ,Input)
     IF (IPRINT.GT. 50) THEN
        CALL printDerivativeOverlapInfo(derivInfo,LUPRI)
     ENDIF
  ENDIF
  P => PQ%P%p
  IF(INPUT%magderOrderP.EQ.1)THEN 
     antiAB=.TRUE.
     IF(P%magderiv.EQ.1)ndim5P=ndim5P*3
  ENDIF
  DO iPassP = 1, P%nPasses
    atoma  = P%orb1atom(iPassP)
    batchA = P%orb1batch(iPassP)
    nA     = P%orbital1%totOrbitals
    atomb  = P%orb2atom(iPassP)
    batchB = P%orb2batch(iPassP)
    nB     = P%orbital2%totOrbitals
    dopermutation = SameLHSaos .AND.((batchA.NE.batchB).OR.(atoma.NE.atomb))
    doAntipermutation = dopermutation.AND.antiAB
    IF (dograd) THEN
      derivInfo%Atom(1)=atomA
      derivInfo%Atom(2)=atomB
      AB    = Input%LST_DLHS%INDEX(atomA,atomB,1,1)
      DAB   => Input%LST_DLHS%LSAO(AB)%elms
      n1 = Input%LST_DLHS%LSAO(AB)%nLocal(1)
      n2 = Input%LST_DLHS%LSAO(AB)%nLocal(2)
      maxBat = Input%LST_DLHS%LSAO(AB)%maxBat
      maxAng = Input%LST_DLHS%LSAO(AB)%maxAng
      sA = Input%LST_DLHS%LSAO(AB)%startLocalOrb(1+(batchA-1)*maxAng) - 1
      sB = Input%LST_DLHS%LSAO(AB)%startLocalOrb(1+(batchB-1)*maxAng+maxAng*maxBat) - 1 
      IF (dopermutation) THEN
        BA  = Input%LST_DLHS%INDEX(atomB,atomA,1,1)
        DBA => Input%LST_DLHS%LSAO(BA)%elms
      ENDIF !dopermutation
      DO iDerivP=1,NgeoderivcompP
        iDer  = derivInfo%dirComp(1,iDerivP)
        iAtom = derivInfo%Atom(derivInfo%AO(1,iDerivP))
!The structure of the lstensor when used to store the molecular gradient
!is maybe not the most logical and we refer to discussion in the subroutine
!init_gradientlstensor in lsutil/lstensor_operations.f90
        iLSAO = Lsoutput%resultTensor%INDEX(iAtom,1,1,derivInfo%AO(1,iDerivP))
#ifdef VAR_LSDEBUGINT
        IF(iLSAO.EQ. 0)THEN
           WRITE(lupri,*)'iLSAO is zero - for some reason'
           WRITE(lupri,*)'iAtom          ',iAtom
           WRITE(lupri,*)'Center (1-4) = ',derivInfo%AO(1,iDerivP)
           WRITE(lupri,*)'Error in use of gradientlstensor distributeJengine'
           CALL LSQUIT('Error in use of gradientlstensor',lupri)
        ENDIF
#endif
        Jgrad => Lsoutput%resultTensor%LSAO(iLSAO)%elms

        CALL addJgradN(Jgrad,n1,n2,sA,sB,QPMAT2,DAB,nMAT,nA,nB,P%nPasses,NgeoderivcompP,iPassP,iDerivP,iDer)
        IF (dopermutation) THEN
          CALL addJgradT(Jgrad,n1,n2,sA,sB,QPMAT2,DBA,nMAT,nA,nB,P%nPasses,NgeoderivcompP,iPassP,iDerivP,iDer)
        ENDIF !dopermutation
      ENDDO !iDerivP
    ELSEIF(doFullContract)THEN
       AB    = Input%LST_DLHS%INDEX(atomA,atomB,1,1)
       DAB   => Input%LST_DLHS%LSAO(AB)%elms
       n1 = Input%LST_DLHS%LSAO(AB)%nLocal(1)
       n2 = Input%LST_DLHS%LSAO(AB)%nLocal(2)
       maxBat = Input%LST_DLHS%LSAO(AB)%maxBat
       maxAng = Input%LST_DLHS%LSAO(AB)%maxAng
       sA = Input%LST_DLHS%LSAO(AB)%startLocalOrb(1+(batchA-1)*maxAng) - 1
       sB = Input%LST_DLHS%LSAO(AB)%startLocalOrb(1+(batchB-1)*maxAng+maxAng*maxBat) - 1 
       IF (dopermutation) THEN
          BA  = Input%LST_DLHS%INDEX(atomB,atomA,1,1)
          DBA => Input%LST_DLHS%LSAO(BA)%elms
       ENDIF !dopermutation
       CALL addJcontN(Jcont,n1,n2,sA,sB,QPMAT2,&
            & DAB,nMAT,nA,nB,P%nPasses,iPassP)
       IF (dopermutation) THEN
        CALL addJcontT(Jcont,n1,n2,sA,sB,QPMAT2,&
             & DBA,nMAT,nA,nB,P%nPasses,iPassP)
       ENDIF !dopermutation
    ELSE
      AB = RES%INDEX(atoma,atomb,1,1)
      JAB => RES%LSAO(AB)%elms
      n1 = res%LSAO(AB)%nLocal(1)
      n2 = res%LSAO(AB)%nLocal(2)
      maxBat = res%LSAO(AB)%maxBat
      maxAng = res%LSAO(AB)%maxAng
      sA = res%LSAO(AB)%startLocalOrb(1+(batchA-1)*maxAng) - 1
      sB = res%LSAO(AB)%startLocalOrb(1+(batchB-1)*maxAng+maxAng*maxBat) - 1 
      IF (dopermutation) THEN
        BA = RES%INDEX(atomb,atoma,1,1)
        JBA => RES%LSAO(BA)%elms        
      ENDIF !dopermutation
      DO iDerivP=1,ndim5P
        idim5 = iDerivP        
        CALL addJab(Jab,n1,n2,sA,sB,QPMAT2,nMAT,nA,nB,P%nPasses,ndim5P,ndim5output,idim5,iPassP,iDerivP)
        IF (doAntipermutation) THEN
          CALL addAntiJba(Jba,n1,n2,sA,sB,QPMAT2,nMAT,nA,nB,P%nPasses,ndim5P,ndim5output,idim5,iPassP,iDerivP)
        ELSEIF (dopermutation) THEN
          CALL addJba(Jba,n1,n2,sA,sB,QPMAT2,nMAT,nA,nB,P%nPasses,ndim5P,ndim5output,idim5,iPassP,iDerivP)
        ENDIF !dopermutation
      ENDDO !iDerivP
    ENDIF ! dograd
  ENDDO !iPassP

  IF (geoderivorder.GT. 0) CALL freeDerivativeOverlapInfo(derivInfo)
           
CONTAINS
SUBROUTINE addJab(Jab,n1,n2,sA,sB,Integrals,nMAT,nA,nB,nPasses,ndim5,ndim5output,idim5,iPass,iDeriv)
implicit none
Integer,intent(IN)        :: nMAT,nA,nB,nPasses,ndim5,ndim5output,idim5,iPass,iDeriv
Integer,intent(IN)        :: n1,n2,sA,sB
Real(realk),intent(IN)    :: Integrals(nMat,nA,nB,ndim5,nPasses)
Real(realk),intent(INOUT) :: Jab(n1,n2,nmat,ndim5output) !(nP,nmat,ndim5output)
!
Integer :: iA,iB,iMat

!$OMP CRITICAL (addJabSection)
DO iMat=1,nMat
  DO iB=1,nB
   DO iA=1,nA
    Jab(sA+iA,sB+iB,iMat,idim5) = Jab(sA+iA,sB+iB,iMat,idim5) + Integrals(iMat,iA,iB,iDeriv,iPass)
   ENDDO
  ENDDO
ENDDO
!$OMP END CRITICAL (addJabSection)
END SUBROUTINE addJab

SUBROUTINE addJba(Jba,n1,n2,sA,sB,Integrals,nMAT,nA,nB,nPasses,ndim5,ndim5output,idim5,iPass,iDeriv)
implicit none
Integer,intent(IN)        :: nMAT,nA,nB,nPasses,ndim5,ndim5output,idim5,iPass,iDeriv
Integer,intent(IN)        :: n1,n2,sA,sB
Real(realk),intent(IN)    :: Integrals(nMat,nA,nB,ndim5,nPasses)
Real(realk),intent(INOUT) :: Jba(n2,n1,nmat,ndim5output)!(nB,nA,nMat,ndim5output)
!
Integer :: iA,iB,iMat

!$OMP CRITICAL (addJabSection)
DO iMat=1,nMat
 DO iA=1,nA
  DO iB=1,nB
   Jba(sB+iB,sA+iA,iMat,idim5) = Jba(sB+iB,sA+iA,iMat,idim5) + Integrals(iMat,iA,iB,iDeriv,iPass)
  ENDDO
 ENDDO
ENDDO
!$OMP END CRITICAL (addJabSection)
END SUBROUTINE addJba

SUBROUTINE addAntiJba(Jba,n1,n2,sA,sB,Integrals,nMAT,nA,nB,nPasses,ndim5,ndim5output,idim5,iPass,iDeriv)
implicit none
Integer,intent(IN)        :: nMAT,nA,nB,nPasses,ndim5,ndim5output,idim5,iPass,iDeriv
Integer,intent(IN)        :: n1,n2,sA,sB
Real(realk),intent(IN)    :: Integrals(nMat,nA,nB,ndim5,nPasses)
Real(realk),intent(INOUT) :: Jba(n2,n1,nmat,ndim5output)!(nB,nA,nMat,ndim5output)
!
Integer :: iA,iB,iMat

!$OMP CRITICAL (addJabSection)
DO iMat=1,nMat
 DO iA=1,nA
  DO iB=1,nB
   Jba(sB+iB,sA+iA,iMat,idim5) = Jba(sB+iB,sA+iA,iMat,idim5) - Integrals(iMat,iA,iB,iDeriv,iPass)
  ENDDO
 ENDDO
ENDDO
!$OMP END CRITICAL (addJabSection)
END SUBROUTINE addAntiJba
!
SUBROUTINE addJgradN(Jgrad,n1,n2,sA,sB,Integrals,Dab,nMat,nA,nB,nPasses,ngeoderivcomp,iPass,iDeriv,iDer)
implicit none
Integer,intent(IN)        :: nMat,nA,nB,nPasses,ngeoderivcomp,iPass,iDeriv,iDer
Integer,intent(IN)        :: n1,n2,sA,sB
Real(realk),intent(INOUT) :: Jgrad(3)
Real(realk),intent(IN)    :: Integrals(nMat,nA,nB,ngeoderivcomp,nPasses)
Real(realk),intent(IN)    :: Dab(n1,n2,nMat) !(nP,nMat)
!
Integer :: iA,iB,iMat
Real(realk) :: Jtmp

Jtmp = 0E0_realk
DO iMat=1,nMat
  DO iB=1,nB
   DO iA=1,nA
     Jtmp = Jtmp + Integrals(iMat,iA,iB,iDeriv,iPass)*Dab(sA+iA,sB+iB,iMat)
   ENDDO
  ENDDO
ENDDO
!$OMP CRITICAL (addJgradNSection)
Jgrad(iDer) = Jgrad(iDer) + Jtmp * 0.5E0_realk
!$OMP END CRITICAL (addJgradNSection)

END SUBROUTINE addJgradN

SUBROUTINE addJgradT(Jgrad,n1,n2,sA,sB,Integrals,Dba,nMat,nA,nB,nPasses,ngeoderivcomp,iPass,iDeriv,iDer)
implicit none
Integer,intent(IN)        :: nMat,nA,nB,nPasses,ngeoderivcomp,iPass,iDeriv,iDer
Integer,intent(IN)        :: n1,n2,sA,sB
Real(realk),intent(INOUT) :: Jgrad(3)
Real(realk),intent(IN)    :: Integrals(nMat,nA,nB,ngeoderivcomp,nPasses)
Real(realk),intent(IN)    :: Dba(n2,n1,nMat)!(nB,nA,nMat)
!
Integer     :: iP,iMat
Real(realk) :: Jtmp

Jtmp = 0E0_realk
DO iMat=1,nMat
  DO iB=1,nB
    DO iA=1,nA
      Jtmp = Jtmp + Integrals(iMat,iA,iB,iDeriv,iPass)*Dba(sB+iB,sA+iA,iMat)
    ENDDO
  ENDDO
ENDDO
!$OMP CRITICAL (addJgradNSection)
Jgrad(iDer) = Jgrad(iDer) + Jtmp * 0.5E0_realk
!$OMP END CRITICAL (addJgradNSection)

END SUBROUTINE addJgradT

SUBROUTINE addJcontN(Jcont,n1,n2,sA,sB,Integrals,Dab,nMat,nA,nB,nPasses,iPass)
implicit none
Integer,intent(IN)        :: nMat,nA,nB,nPasses,iPass
Integer,intent(IN)        :: n1,n2,sA,sB
Real(realk),intent(INOUT) :: Jcont(nmat)
Real(realk),intent(IN)    :: Integrals(nMat,nA,nB,nPasses)
Real(realk),intent(IN)    :: Dab(n1,n2,nMat)
!
Integer :: iP,iMat
Real(realk) :: Jtmp

DO iMat=1,nMat
  Jtmp = 0E0_realk
  DO iB=1,nB
   DO iA=1,nA
    Jtmp = Jtmp + Integrals(iMat,iA,iB,iPass)*Dab(sA+iA,sB+iB,iMat)
   ENDDO
  ENDDO
  Jcont(iMat) = Jcont(iMat) + Jtmp 
ENDDO

END SUBROUTINE addJcontN

SUBROUTINE addJcontT(Jcont,n1,n2,sA,sB,Integrals,Dba,nMat,nA,nB,nPasses,iPass)
implicit none
Integer,intent(IN)        :: nMat,nA,nB,nPasses,iPass
Integer,intent(IN)        :: n1,n2,sA,sB
Real(realk),intent(INOUT) :: Jcont(nmat)
Real(realk),intent(IN)    :: Integrals(nMat,nA,nB,nPasses)
Real(realk),intent(IN)    :: Dba(n2,n1,nMat)
!
Integer     :: iP,iMat
Real(realk) :: Jtmp

DO iMat=1,nMat
  Jtmp = 0E0_realk
  DO iB=1,nB
    DO iA=1,nA
      Jtmp = Jtmp + Integrals(iMat,iA,iB,iPass)*Dba(sB+iB,sA+iA,iMat)
    ENDDO
  ENDDO
  Jcont(iMat) = Jcont(iMat) + Jtmp
ENDDO

END SUBROUTINE addJcontT

end SUBROUTINE distributeJengine

!> \brief new distributePQ to lstensor
!> \author T. Kjaergaard
!> \date 2009-02-05
!> \param RES the output lstensor
!> \param PQ contain info about the overlap distributions
!> \param QPmat2 matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param input the integral input contains integral specification 
!> \param output the integral output specification NOT NEEDED 
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE distributeCS(RES,PQ,QPmat2,dimQ,dimP,Input,Output2,LUPRI,IPRINT)
implicit none 
Type(integrand),intent(in)      :: PQ
Type(IntegralInput),intent(in)  :: Input
Type(IntegralOutput),intent(inout) :: Output2
Integer,intent(in)              :: LUPRI,IPRINT
Integer,intent(in)              :: dimQ,dimP
REAL(REALK),intent(in)          :: QPMAT2(dimQ,dimP)
TYPE(lstensor),intent(inout)    :: RES
!
Integer               :: nA,nB, AtomA,atomB
logical               :: SameLHSaos,dopermutation
integer               :: IA,IB,batchA,batchB,n1,n2,maxBat,maxAng,sA,sB,AB,BA
TYPE(Overlap),pointer :: P
real(realk) :: maxgabelm
SameLHSaos = INPUT%SameLHSaos
P => PQ%P%p
IF (P%nPasses.GT. 1) CALL LSQUIT('Error in distributeCS. nPasses > 1',-1)

!Beware when converting from double precision to short integer 
!If double precision is less than 10^-33 then you can run into
!problems with short integer overflow
nA     = P%orbital1%totOrbitals
nB     = P%orbital2%totOrbitals
IF(Output2%RealGabMatrix)THEN
   atoma  = P%orb1atom(1)
   batchA = P%orb1batch(1)
   atomb  = P%orb2atom(1)
   batchB = P%orb2batch(1)
   AB = output2%resultTensor%INDEX(AtomA,AtomB,1,1) 
   n1 = output2%resultTensor%LSAO(AB)%nLocal(1)
   n2 = output2%resultTensor%LSAO(AB)%nLocal(2)
   maxBat = output2%resultTensor%LSAO(AB)%maxBat
   maxAng = output2%resultTensor%LSAO(AB)%maxAng
   sA = output2%resultTensor%LSAO(AB)%startLocalOrb(1+(batchA-1)*maxAng) - 1
   sB = output2%resultTensor%LSAO(AB)%startLocalOrb(1+(batchB-1)*maxAng+maxAng*maxBat) - 1 
   CALL BuildGab(QPMAT2,nA,nB,output2%resultTensor%LSAO(AB)%elms,n1,n2,sA,sB)
   dopermutation = SameLHSaos .AND.((batchA.NE.batchB).OR.(atoma.NE.atomb))
   IF (dopermutation) THEN
      BA = output2%resultTensor%INDEX(AtomB,AtomA,1,1) 
      CALL BuildGabPermute(QPMAT2,nA,nB,output2%resultTensor%LSAO(BA)%elms,n1,n2,sA,sB)
   ENDIF
ELSE
   CALL findMaxAbsGabElm(QPMAT2,nA*nB,maxgabelm,input%CS_THRESHOLD,lupri)

   atomA   = P%orb1atom(1)
   atomB   = P%orb2atom(1)
   IA   = P%orb1batch(1)
   IB   = P%orb2batch(1)
!$OMP CRITICAL (FullGab)
   IF (maxgabelm.GT.shortintCRIT) THEN
      RES%maxgab(IA,IB) = CEILING(LOG10(sqrt(maxgabelm)))
   ELSE
      RES%maxgab(IA,IB) = shortzero !meaning -33=> 10**-33=0
   ENDIF
   dopermutation = SameLHSaos .AND.((IA.NE.IB).OR.(atoma.NE.atomb))
   IF (dopermutation) THEN
      RES%maxgab(IB,IA) =  RES%maxgab(IA,IB)
   ENDIF
!$OMP END CRITICAL (FullGab)
ENDIF
CONTAINS
SUBROUTINE findMaxAbsGabElm(Gabint,nP,maxgabelm,thresh,lupri)
implicit none
Integer,intent(IN)        :: nP,lupri
Real(realk),intent(IN)    :: Gabint(nP,nP),thresh
Real(realk),intent(INOUT) :: maxgabelm
!
real(realk) :: maxgabelm2
Integer :: iP
!write(lupri,*)'Gabint findMaxAbsGabElm'
!call ls_output(Gabint,1,nP,1,nP,nP,nP,1,lupri)
maxgabelm = 0.0E0_realk
DO iP=1,nP
#ifdef VAR_LSDEBUGINT
  IF(Gabint(iP,iP).GT. shortintCrit)THEN
#endif
     IF (Gabint(iP,iP).GT. maxgabelm) THEN
        maxgabelm = Gabint(iP,iP)
     ENDIF
#ifdef VAR_LSDEBUGINT
  ELSE
     IF(abs(Gabint(iP,iP)).GT.thresh) THEN
        write(lupri,*) 'Error in addGab. Negative screening int',Gabint(iP,iP)
        call lsquit('Error in addGab. Negative screening int',lupri)
     ENDIF
  ENDIF
#endif
ENDDO
END SUBROUTINE findMaxAbsGabElm

SUBROUTINE BuildGab(Gabint,nA,nB,OutputMatrix,n1,n2,sA,sB)
implicit none
Integer,intent(IN)        :: n1,n2,sA,sB,nA,nB
Real(realk),intent(IN)    :: Gabint(nA,nB,nA,nB)
Real(realk),intent(INOUT) :: OutputMatrix(n1,n2)
!
Integer :: iA,iB
!write(lupri,*)'Gabint BuildGab'
!call ls_output(Gabint,1,nA*nB,1,nA*nB,nA*nB,nA*nB,1,lupri)
DO iB=1,nB
   DO iA=1,nA
#ifdef VAR_LSDEBUGINT
      IF(Gabint(iA,iB,iA,iB).LT. -1.0E-12_realk)THEN
         print*, 'Error in buildGab. Negative screening int',Gabint(iA,iB,iA,iB)
         call lsquit('Error in buildGab. Negative screening int',-1)
      ENDIF
#endif
      OutputMatrix(sA+iA,sB+iB) = sqrt(abs(Gabint(iA,iB,iA,iB)))      
   ENDDO
ENDDO
END SUBROUTINE BuildGab

SUBROUTINE BuildGabPermute(Gabint,nA,nB,OutputMatrix,n1,n2,sA,sB)
implicit none
Integer,intent(IN)        :: n1,n2,sA,sB,nA,nB
Real(realk),intent(IN)    :: Gabint(nA,nB,nA,nB)
Real(realk),intent(INOUT) :: OutputMatrix(n2,n1)
!
Integer :: iA,iB
DO iB=1,nB
   DO iA=1,nA
#ifdef VAR_LSDEBUGINT
      IF(Gabint(iA,iB,iA,iB).LT. -1.0E-12_realk)THEN
         print*, 'Error in buildGab. Negative screening int',Gabint(iA,iB,iA,iB)
         call lsquit('Error in buildGab. Negative screening int',-1)
      ENDIF
#endif
      OutputMatrix(sB+iB,sA+iA) = sqrt(abs(Gabint(iA,iB,iA,iB)))
   ENDDO
ENDDO
END SUBROUTINE BuildGabPermute

END SUBROUTINE distributeCS

!> \brief new distributePQ to lstensor
!> \author T. Kjaergaard
!> \date 2009-02-05
!> \param RES the output lstensor
!> \param PQ contain info about the overlap distributions
!> \param QPmat2 matrix containing calculated integrals
!> \param dimQ dimension 1 of QPmat
!> \param dimP dimension 2 of QPmat
!> \param input the integral input contains integral specification 
!> \param output the integral output specification 
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE distributePS(RES3,PQ,QPmat2,dimQ,dimP,Input,Output2,LUPRI,IPRINT)
implicit none 
Type(integrand),intent(in)      :: PQ
Type(IntegralInput),intent(in)  :: Input
Type(IntegralOutput),intent(inout) :: Output2
Integer,intent(in)              :: LUPRI,IPRINT
Integer,intent(in)              :: dimQ,dimP
REAL(REALK),intent(in)          :: QPMAT2(dimQ,dimP)
TYPE(lstensor),intent(inout)    :: RES3
!
Integer :: nAngmomP,nOrbP,totOrbQ,endOrbP,totOrbP,nMAT
Integer :: iAngmomP,iOrbitalP,iderivP
Integer :: iA,iB,sA,sB,sC,sD,SAA,SBB,SCC,SDD
Integer :: nContA,nContB,nAngA,nAngB,startA,startB
Integer :: iOrbP,iContP,iderivQ,ideriv
Integer :: iOrbQ,idmat
Integer :: startA1,startB1
real(realk) :: THRESHOLD,TMP
integer :: PASSPOFFSET,PASSQOFFSET
integer :: iPassP
Integer :: n1,n2,dimA,dimB,nA,nB
Integer :: AtomA,atomB!,atomC,atomD
!Integer :: ABCD,BACD,ABDC,BADC,CDAB,CDBA,DCAB,DCBA,ngeoderivcomp
Integer :: ngeoderivcomp
logical :: samePQ,SameLHSaos,dopermutation
integer :: batchA,batchB,AB,BA,cBatchA,cBatchB,maxbat
real(realk),pointer :: dummy(:)
!
integer(kind=short),pointer   :: Gab(:),Gba(:)
TYPE(Overlap),pointer :: P
SameLHSaos = INPUT%SameLHSaos
P => PQ%P%p
atomA   = P%orb1atom(1)
batchA  = P%orb1batch(1)
nA      = P%orbital1%totOrbitals
atomB   = P%orb2atom(1)
batchB  = P%orb2batch(1)
nB      = P%orbital2%totOrbitals

IF(P%nPasses.GT. 1)CALL LSQUIT('Error in distributelstensorCS. nPasses > 1',-1)

dopermutation = SameLHSaos .AND.((batchA.NE.batchB).OR.(atoma.NE.atomb))
AB = RES3%INDEX(atomA,atomB,1,1)
IF(AB.EQ.0)THEN
!   print*,'ERROR atomA,ATomB',AtomA,ATomB
!   call ls_output(QPMAT2,1,dimQ,1,dimP,dimQ,dimP,1,6)
!   call lsquit('PS OD SCREEN something wrong?  input%od_screen',-1)
   RETURN
ENDIF
Gab => RES3%SLSAO(AB)%selms
n1 = RES3%SLSAO(AB)%nLocal(1)
n2 = RES3%SLSAO(AB)%nLocal(2)
maxBat = RES3%SLSAO(AB)%maxBat
sA = RES3%SLSAO(AB)%startLocalOrb(batchA) - 1
sB = RES3%SLSAO(AB)%startLocalOrb(batchB+maxBat) - 1 
nContA = P%orbital1%nContracted(1)
nContB = P%orbital2%nContracted(1)
IF (dopermutation) THEN
   BA = RES3%INDEX(atomB,atomA,1,1)
   Gba => RES3%SLSAO(BA)%selms
ENDIF

!iOrbitalP = 0
DO iAngmomP=1,P%nAngmom
   iA = P%indexAng1(iAngmomP)
   iB = P%indexAng2(iAngmomP)
   nAngA = P%orbital1%nOrbComp(iA)
   nAngB = P%orbital2%nOrbComp(iB)
   startA = P%orbital1%startLocOrb(iA)
   startB = P%orbital2%startLocOrb(iB)
   CALL addPSGab(Gab,n1,n2,sA,sB,QPMAT2,dimQ,dimP,nA,nB,nContA,nAngA,nContB,nAngB,startA,startB,&
        & input%PS_THRESHOLD,input%contang,lupri)
   IF (dopermutation) THEN
      CALL addPSGba(Gba,n1,n2,sA,sB,QPMAT2,dimQ,dimP,nA,nB,nContA,nAngA,nContB,nAngB,startA,startB,&
           & input%PS_THRESHOLD,input%contang,lupri)
   ENDIF
ENDDO

CONTAINS
SUBROUTINE addPSGab(PSGab,n1,n2,sA,sB,PSGabint,dimQ,dimP,nA,nB,nContA,nAngA,&
     &              nContB,nAngB,startA,startB,thresh,contang,lupri)
implicit none
Integer,intent(IN)      :: lupri,nContA,nAngA,nContB,nAngB,dimQ,dimP,nA,nB,startA,startB
Integer,intent(IN)      :: n1,n2,sA,sB
integer(kind=short),intent(INOUT) :: PSGab(n1,n2)!(nContA,nContB)
Real(realk),intent(IN)  :: PSGabint(nA,nB,nA,nB)
Real(realk),intent(IN)  :: thresh
Logical,intent(IN)      :: contang

IF (.NOT.contang) THEN
! Default ordering of angular components first, contracted second
  CALL addPSGab_ac(PSGab,n1,n2,sA,sB,PSGabint,dimQ,dimP,nA,nB,nContA,nAngA,&
     & nContB,nAngB,startA,startB,thresh,lupri)
ELSE
! CONTANG ordering of contracted components first, angular second
  CALL addPSGab_ca(PSGab,n1,n2,sA,sB,PSGabint,dimQ,dimP,nA,nB,nContA,nAngA,&
     & nContB,nAngB,startA,startB,thresh,lupri)
ENDIF
 
END SUBROUTINE addPSGab

SUBROUTINE addPSGba(PSGba,n1,n2,sA,sB,PSGabint,dimQ,dimP,nA,nB,nContA,nAngA,&
     &              nContB,nAngB,startA,startB,thresh,contang,lupri)
implicit none
Integer,intent(IN)      :: lupri,nContA,nAngA,nContB,nAngB,dimQ,dimP,nA,nB,startA,startB
Integer,intent(IN)      :: n1,n2,sA,sB
integer(kind=short),intent(INOUT) :: PSGba(n2,n1)!(nContB,nContA)
Real(realk),intent(IN)  :: PSGabint(nA,nB,nA,nB)
Real(realk),intent(IN)  :: thresh
Logical,intent(IN)      :: contang

IF (.NOT.contang) THEN
! Default ordering of angular components first, contracted second
  CALL addPSGba_ac(PSGba,n1,n2,sA,sB,PSGabint,dimQ,dimP,nA,nB,nContA,nAngA,&
     & nContB,nAngB,startA,startB,thresh,lupri)
ELSE
! CONTANG ordering of contracted components first, angular second
  CALL addPSGba_ca(PSGba,n1,n2,sA,sB,PSGabint,dimQ,dimP,nA,nB,nContA,nAngA,&
     & nContB,nAngB,startA,startB,thresh,lupri)
ENDIF
 
END SUBROUTINE addPSGba

SUBROUTINE addPSGab_ac(PSGab,n1,n2,sA,sB,PSGabint,dimQ,dimP,nA,nB,nContA,nAngA,&
     & nContB,nAngB,startA,startB,thresh,lupri)
implicit none
Integer,intent(IN)      :: lupri,nContA,nAngA,nContB,nAngB,dimQ,dimP,nA,nB,startA,startB
Integer,intent(IN)      :: n1,n2,sA,sB
integer(kind=short),intent(INOUT) :: PSGab(n1,n2)!(nContA,nContB)
Real(realk),intent(IN)  :: PSGabint(nA,nB,nA,nB)
Real(realk),intent(IN)  :: thresh
!
integer(kind=short) :: PSGabIntTmp
Integer :: iContB,iContA,iAngA,iAngB,iOrbP,iA,iB,offset

iB = startB -1
DO iContB=1,nContB
 DO iAngB=1,nAngB
  iB = iB + 1
  iA = startA-1
  DO iContA=1,nContA
   DO iAngA=1,nAngA
    iA = iA + 1
    !Beware when converting from double precision to short integer 
    !If double precision is less than 10^-33 then you can run into
    !problems with short integer overflow
    IF(PSGabint(iA,iB,iA,iB).GT. shortintcrit)THEN
       PSgabIntTmp = CEILING(LOG10(sqrt(PSGabint(iA,iB,iA,iB))))
       PSGab(sA+iContA,sB+iContB)=MAX(PSGab(sA+iContA,sB+iContB),PSgabIntTmp)
    ELSE
       IF (abs(PSGabint(iA,iB,iA,iB)).GT.thresh) THEN
          PSgabIntTmp = CEILING(LOG10(sqrt(ABS(PSGabint(iA,iB,iA,iB)))))
          PSGab(sA+iContA,sB+iContB)=MAX(PSGab(sA+iContA,sB+iContB),PSgabIntTmp)
!          write(lupri,*) 'Error in addPSGab. Negative screening integral',PSGabint(iA,iB,iA,iB)
!          call lsquit('Error in addPSGab. Negative screening integral',lupri)
       ENDIF       
    ENDIF
   ENDDO
  ENDDO
 ENDDO
ENDDO

END SUBROUTINE addPSGab_ac

SUBROUTINE addPSGba_ac(PSGba,n1,n2,sA,sB,PSGabint2,dimQ,dimP,nA,nB,nContA,nAngA,&
     &nContB,nAngB,startA,startB,thresh,lupri)
implicit none
Integer,intent(IN)      :: lupri,nContA,nAngA,nContB,nAngB,dimQ,dimP,nA,nB,startA,startB
Integer,intent(IN)      :: n1,n2,sA,sB
integer(kind=short),intent(INOUT) :: PSGba(n2,n1)!(nContB,nContA)
Real(realk),intent(IN)  :: PSGabint2(nA,nB,nA,nB)
!Real(realk),intent(IN)  :: PSGabint2(dimQ,dimP)
Real(realk),intent(IN)  :: thresh
!
integer(kind=short) :: PSGabIntTmp
Integer :: iA,iB,iContB,iContA,iAngA,iAngB,offset,iOrbP

iB = startB -1
DO iContB=1,nContB
 DO iAngB=1,nAngB
  iB = iB + 1
  iA = startA-1
  DO iContA=1,nContA
   DO iAngA=1,nAngA
    iA = iA + 1
    !Beware when converting from double precision to short integer 
    !If double precision is less than 10^-33 then you can run into
    !problems with short integer overflow
    IF(PSGabint2(iA,iB,iA,iB).GT. shortintcrit)THEN
       PSgabIntTmp = CEILING(LOG10(sqrt(PSGabint2(iA,iB,iA,iB))))
!$OMP CRITICAL (PSGab)
       PSGba(sB+iContB,sA+iContA) = MAX(PSGba(sB+iContB,sA+iContA),PSgabIntTmp)
!$OMP END CRITICAL (PSGab)
    ELSE
       IF(ABS(PSGabint2(iA,iB,iA,iB)).GT.thresh)THEN
          PSgabIntTmp = CEILING(LOG10(sqrt(ABS(PSGabint2(iA,iB,iA,iB)))))
!$OMP CRITICAL (PSGab)
          PSGba(sB+iContB,sA+iContA) = MAX(PSGba(sB+iContB,sA+iContA),PSgabIntTmp)
!$OMP END CRITICAL (PSGab)

!          write(lupri,*)'Error in addPSGba. Negative screening integral',PSGabint2(iA,iB,iA,iB)
!          call lsquit('Error in addPSGba. Negative screening integral',lupri)
       ENDIF
    ENDIF
   ENDDO
  ENDDO
 ENDDO
ENDDO

END SUBROUTINE addPSGba_ac

SUBROUTINE addPSGab_ca(PSGab,n1,n2,sA,sB,PSGabint,dimQ,dimP,nA,nB,nContA,nAngA,&
     & nContB,nAngB,startA,startB,thresh,lupri)
implicit none
Integer,intent(IN)      :: lupri,nContA,nAngA,nContB,nAngB,dimQ,dimP,nA,nB,startA,startB
Integer,intent(IN)      :: n1,n2,sA,sB
integer(kind=short),intent(INOUT) :: PSGab(n1,n2)!(nContA,nContB)
Real(realk),intent(IN)  :: PSGabint(nA,nB,nA,nB)
Real(realk),intent(IN)  :: thresh
!
integer(kind=short) :: PSGabIntTmp
Integer :: iContB,iContA,iAngA,iAngB,iOrbP,iA,iB,offset

iB = startB -1
DO iAngB=1,nAngB
 DO iContB=1,nContB
  iB = iB + 1
  iA = startA-1
  DO iAngA=1,nAngA
   DO iContA=1,nContA
    iA = iA + 1
    !Beware when converting from double precision to short integer 
    !If double precision is less than 10^-33 then you can run into
    !problems with short integer overflow
    IF(PSGabint(iA,iB,iA,iB).GT. shortintcrit)THEN
       PSgabIntTmp = CEILING(LOG10(sqrt(PSGabint(iA,iB,iA,iB))))
       PSGab(sA+iContA,sB+iContB)=MAX(PSGab(sA+iContA,sB+iContB),PSgabIntTmp)
    ELSE
       IF (abs(PSGabint(iA,iB,iA,iB)).GT.thresh) THEN
          PSgabIntTmp = CEILING(LOG10(sqrt(ABS(PSGabint(iA,iB,iA,iB)))))
          PSGab(sA+iContA,sB+iContB)=MAX(PSGab(sA+iContA,sB+iContB),PSgabIntTmp)
!          write(lupri,*) 'Error in addPSGab. Negative screening integral',PSGabint(iA,iB,iA,iB)
!          call lsquit('Error in addPSGab. Negative screening integral',lupri)
       ENDIF       
    ENDIF
   ENDDO
  ENDDO
 ENDDO
ENDDO

END SUBROUTINE addPSGab_ca

SUBROUTINE addPSGba_ca(PSGba,n1,n2,sA,sB,PSGabint2,dimQ,dimP,nA,nB,nContA,nAngA,&
     &nContB,nAngB,startA,startB,thresh,lupri)
implicit none
Integer,intent(IN)      :: lupri,nContA,nAngA,nContB,nAngB,dimQ,dimP,nA,nB,startA,startB
Integer,intent(IN)      :: n1,n2,sA,sB
integer(kind=short),intent(INOUT) :: PSGba(n2,n1)!(nContB,nContA)
Real(realk),intent(IN)  :: PSGabint2(nA,nB,nA,nB)
!Real(realk),intent(IN)  :: PSGabint2(dimQ,dimP)
Real(realk),intent(IN)  :: thresh
!
integer(kind=short) :: PSGabIntTmp
Integer :: iA,iB,iContB,iContA,iAngA,iAngB,offset,iOrbP

iB = startB -1
DO iAngB=1,nAngB
 DO iContB=1,nContB
  iB = iB + 1
  iA = startA-1
  DO iAngA=1,nAngA
   DO iContA=1,nContA
    iA = iA + 1
    !Beware when converting from double precision to short integer 
    !If double precision is less than 10^-33 then you can run into
    !problems with short integer overflow
    IF(PSGabint2(iA,iB,iA,iB).GT. shortintcrit)THEN
       PSgabIntTmp = CEILING(LOG10(sqrt(PSGabint2(iA,iB,iA,iB))))
       PSGba(sB+iContB,sA+iContA) = MAX(PSGba(sB+iContB,sA+iContA),PSgabIntTmp)  
    ELSE
       IF(ABS(PSGabint2(iA,iB,iA,iB)).GT.thresh)THEN
          PSgabIntTmp = CEILING(LOG10(sqrt(ABS(PSGabint2(iA,iB,iA,iB)))))
          PSGba(sB+iContB,sA+iContA) = MAX(PSGba(sB+iContB,sA+iContA),PSgabIntTmp)
!          write(lupri,*)'Error in addPSGba. Negative screening integral',PSGabint2(iA,iB,iA,iB)
!          call lsquit('Error in addPSGba. Negative screening integral',lupri)
       ENDIF
    ENDIF
   ENDDO
  ENDDO
 ENDDO
ENDDO

END SUBROUTINE addPSGba_ca

END SUBROUTINE distributePS

END MODULE thermite_distribute2


