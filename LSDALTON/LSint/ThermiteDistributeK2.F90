!> @file 
!> Contains integral distribution routines that places calculated integrals in the proper output
!> Thermite integral distribution module
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
MODULE thermite_distributeK2
  use thermite_distribute
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
!> \brief Distribute exchange-type integrals
!> \author \latexonly S. Reine  \endlatexonly
!> \date 2011-03-08
!> \param Kmat contains the result lstensor
!> \param CDAB matrix containing calculated integrals
!> \param Dmat contains the rhs-density matrix/matrices
!> \param ndmat number of rhs-density matries
!> \param PQ information about the calculated integrals (to be distributed)
!> \param Input contain info about the requested integral 
!> \param Lsoutput contain info about the requested output, and actual output
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE distributeKgrad(Kgrad,CDAB,Dlhs,Drhs,ndmat,PQ,Input,&
    &                      Lsoutput,LUPRI,IPRINT)
implicit none 
Type(integrand),intent(in)         :: PQ
Type(IntegralInput),intent(in)     :: Input
Type(IntegralOutput),intent(inout) :: Lsoutput
Integer,intent(in)                 :: ndmat,LUPRI,IPRINT
REAL(REALK),pointer                :: CDAB(:)
TYPE(lstensor),intent(inout)       :: Kgrad
TYPE(lstensor),intent(in)          :: Dlhs,Drhs
!
type(overlap),pointer :: P,Q
logical :: SamePQ,SameLHSaos,SameRHSaos,SameODs
Logical :: permuteAB,permuteCD,permuteOD
integer :: nAngmomP,nDeriv,nPasses,iLSAO
integer :: iPass,iPassP,iPassQ,iDeriv,iDerivP,iDerivQ,geoderivorder,iAtom,iDer
integer :: nA,nB,nC,nD,batchA,batchB,batchC,batchD,atomA,atomB,atomC,atomD
!integer :: indLAC,indLBC,indLAD,indLBD,indRBD,indRAD,indRBC,indRAC
!integer :: indLCA,indLCB,indLDA,indLDB,indRDB,indRDA,indRCB,indRCA
integer :: CA,DB,CB,DA,maxAng,maxBat
integer :: AC,BD,niderindex
integer :: nl1,nl2,nl3,nl4,sl1,sl2,sl3,sl4,nr1,nr2,nr3,nr4,sr1,sr2,sr3,sr4
Real(realk),pointer :: lAC(:),lBC(:),lAD(:),lBD(:),rBD(:),rAD(:),rBC(:),rAC(:)
Real(realk),pointer :: lCA(:),lCB(:),lDA(:),lDB(:),rDB(:),rDA(:),rCB(:),rCA(:)
Real(realk),pointer :: kelms(:)
real(realk),pointer :: ltheta(:,:,:,:)
integer :: iLSAOindex(PQ%P%p%nPasses*PQ%Q%p%nPasses*Input%ngeoderivcomp)
integer :: iderindex(PQ%P%p%nPasses*PQ%Q%p%nPasses*Input%ngeoderivcomp)
real(realk) :: idertmp(PQ%P%p%nPasses*PQ%Q%p%nPasses*Input%ngeoderivcomp)
integer,parameter,dimension(4) :: permuteAO = (/1,3,2,4/)

Type(derivativeInfo) :: derivInfo
Real(realk) :: tmp
P => PQ%P%p
Q => PQ%Q%p
nAngmomP   = P%nAngmom
SamePQ     = PQ%samePQ
SameLHSaos = INPUT%SameLHSaos
SameRHSaos = INPUT%SameRHSaos
SameODs    = INPUT%SameODs
nDeriv     = Input%ngeoderivcomp
geoderivorder = Input%geoderivorder
CALL initDerivativeOverlapInfo(derivInfo,PQ,Input)
IF (IPRINT.GT. 50) THEN
  CALL printDerivativeOverlapInfo(derivInfo,LUPRI)
ENDIF

niderindex = 0
nPasses = P%nPasses*Q%nPasses
iPass = 0
DO iPassP=1,P%nPasses
  atomA  = P%orb1atom(iPassP)
  batchA = P%orb1batch(iPassP)
  nA     = P%orbital1%totOrbitals
  atomB  = P%orb2atom(iPassP)
  batchB = P%orb2batch(iPassP)
  nB     = P%orbital2%totOrbitals
  derivInfo%Atom(1)=AtomA
  derivInfo%Atom(2)=AtomB

  permuteAB  = SameLHSaos.AND.((batchA.NE.batchB).OR.(atomA.NE.atomB))
  DO iPassQ=1,Q%nPasses
    iPass = iPass + 1

    atomC  = Q%orb1atom(iPassQ)
    batchC = Q%orb1batch(iPassQ)
    nC     = Q%orbital1%totOrbitals
    atomD  = Q%orb2atom(iPassQ)
    batchD = Q%orb2batch(iPassQ)
    nD     = Q%orbital2%totOrbitals
    derivInfo%Atom(3)=atomC
    derivInfo%Atom(4)=atomD
    permuteCD = SameRHSaos.AND.((batchC.NE.batchD).OR.(atomC.NE.atomD))
    permuteOD = SameODs.AND.( ((batchA.NE.batchC).OR.(atomA.NE.atomC)).OR.&
   &                          ((batchB.NE.batchD).OR.(atomB.NE.atomD)) )
    permuteOD = SameODs.AND..NOT.((batchA.EQ.batchC).AND.(atomA.EQ.atomC).AND.&
   &                          (batchB.EQ.batchD).AND.(atomB.EQ.atomD) )

    call mem_alloc(ltheta,nC,nD,nA,nB)
    ltheta = 0E0_realk
    IF (permuteAB.AND.permuteCD) THEN
       CA = Dlhs%index(atomC,atomA,1,1)
       lCA => Dlhs%lsao(CA)%elms
       DB = Dlhs%index(atomD,atomB,1,1)
       lDB => Dlhs%lsao(DB)%elms
       
       nl3 = Dlhs%lsao(CA)%nLocal(1) 
       nl1 = Dlhs%lsao(CA)%nLocal(2) 
       nl4 = Dlhs%lsao(DB)%nLocal(1) 
       nl2 = Dlhs%lsao(DB)%nLocal(2) 
       maxAng = Dlhs%lsao(CA)%maxAng
       maxBat = Dlhs%lsao(CA)%maxBat
       sl3 = Dlhs%lsao(CA)%startLocalOrb(1+(batchC-1)*maxAng)-1
       sl1 = Dlhs%lsao(CA)%startLocalOrb(1+(batchA-1)*maxAng+maxAng*maxBat)-1
       maxAng = Dlhs%lsao(DB)%maxAng
       maxBat = Dlhs%lsao(DB)%maxBat
       sl4 = Dlhs%lsao(DB)%startLocalOrb(1+(batchD-1)*maxAng)-1
       sl2 = Dlhs%lsao(DB)%startLocalOrb(1+(batchB-1)*maxAng+maxAng*maxBat)-1

       CB = Dlhs%index(atomC,atomB,1,1)
       lCB => Dlhs%lsao(CB)%elms
       DA = Dlhs%index(atomD,atomA,1,1)
       lDA => Dlhs%lsao(DA)%elms

       DB = Drhs%index(atomD,atomB,1,1)
       rDB => Drhs%lsao(DB)%elms
       CA = Drhs%index(atomC,atomA,1,1)
       rCA => Drhs%lsao(CA)%elms
       
       nr4 = Drhs%lsao(DB)%nLocal(1) 
       nr2 = Drhs%lsao(DB)%nLocal(2) 
       nr3 = Drhs%lsao(CA)%nLocal(1) 
       nr1 = Drhs%lsao(CA)%nLocal(2) 
       maxAng = Drhs%lsao(DB)%maxAng
       maxBat = Drhs%lsao(DB)%maxBat
       sr4 = Drhs%lsao(DB)%startLocalOrb(1+(batchD-1)*maxAng)-1
       sr2 = Drhs%lsao(DB)%startLocalOrb(1+(batchB-1)*maxAng+maxAng*maxBat)-1
       maxAng = Drhs%lsao(CA)%maxAng
       maxBat = Drhs%lsao(CA)%maxBat
       sr3 = Drhs%lsao(CA)%startLocalOrb(1+(batchC-1)*maxAng)-1
       sr1 = Drhs%lsao(CA)%startLocalOrb(1+(batchA-1)*maxAng+maxAng*maxBat)-1

       DA = Drhs%index(atomD,atomA,1,1)
       rDA => Drhs%lsao(DA)%elms
       CB = Drhs%index(atomC,atomB,1,1)
       rCB => Drhs%lsao(CB)%elms

       call build_ltheta_TT(ltheta,lCA,rDB,lCB,rDA,lDA,rCB,lDB,rCA,&
            & nl1,nl2,nl3,nl4,sl1,sl2,sl3,sl4,nr1,nr2,nr3,nr4,sr1,sr2,sr3,sr4,nA,nB,nC,nD,ndmat)
    ELSE IF (permuteAB) THEN
       CA = Dlhs%index(atomC,atomA,1,1)
       lCA => Dlhs%lsao(CA)%elms
       CB = Dlhs%index(atomC,atomB,1,1)
       lCB => Dlhs%lsao(CB)%elms
       
       nl3 = Dlhs%lsao(CA)%nLocal(1) 
       nl1 = Dlhs%lsao(CA)%nLocal(2) 
       nl2 = Dlhs%lsao(CB)%nLocal(2) 
       maxAng = Dlhs%lsao(CA)%maxAng
       maxBat = Dlhs%lsao(CA)%maxBat
       sl3 = Dlhs%lsao(CA)%startLocalOrb(1+(batchC-1)*maxAng)-1
       sl1 = Dlhs%lsao(CA)%startLocalOrb(1+(batchA-1)*maxAng+maxAng*maxBat)-1
       maxAng = Dlhs%lsao(CB)%maxAng
       maxBat = Dlhs%lsao(CB)%maxBat
       sl2 = Dlhs%lsao(CB)%startLocalOrb(1+(batchB-1)*maxAng+maxAng*maxBat)-1

       DB = Drhs%index(atomD,atomB,1,1)
       rDB => Drhs%lsao(DB)%elms
       DA = Drhs%index(atomD,atomA,1,1)
       rDA => Drhs%lsao(DA)%elms
       
       nr4 = Drhs%lsao(DB)%nLocal(1) 
       nr2 = Drhs%lsao(DB)%nLocal(2) 
       nr1 = Drhs%lsao(DA)%nLocal(2) 
       maxAng = Drhs%lsao(DB)%maxAng
       maxBat = Drhs%lsao(DB)%maxBat
       sr4 = Drhs%lsao(DB)%startLocalOrb(1+(batchD-1)*maxAng)-1
       sr2 = Drhs%lsao(DB)%startLocalOrb(1+(batchB-1)*maxAng+maxAng*maxBat)-1
       maxAng = Drhs%lsao(DA)%maxAng
       maxBat = Drhs%lsao(DA)%maxBat
       sr1 = Drhs%lsao(DA)%startLocalOrb(1+(batchA-1)*maxAng+maxAng*maxBat)-1

       call build_ltheta_TF(ltheta,lCA,rDB,lCB,rDA,&
            & nl1,nl2,nl3,sl1,sl2,sl3,nr1,nr2,nr4,sr1,sr2,sr4,&
            & nA,nB,nC,nD,ndmat)
    ELSE IF (permuteCD) THEN
       CA = Dlhs%index(atomC,atomA,1,1)
       lCA => Dlhs%lsao(CA)%elms
       DA = Dlhs%index(atomD,atomA,1,1)
       lDA => Dlhs%lsao(DA)%elms
       
       nl3 = Dlhs%lsao(CA)%nLocal(1) 
       nl1 = Dlhs%lsao(CA)%nLocal(2) 
       nl4 = Dlhs%lsao(DA)%nLocal(1) 
       maxAng = Dlhs%lsao(CA)%maxAng
       maxBat = Dlhs%lsao(CA)%maxBat
       sl3 = Dlhs%lsao(CA)%startLocalOrb(1+(batchC-1)*maxAng)-1
       sl1 = Dlhs%lsao(CA)%startLocalOrb(1+(batchA-1)*maxAng+maxAng*maxBat)-1
       maxAng = Dlhs%lsao(DA)%maxAng
       maxBat = Dlhs%lsao(DA)%maxBat
       sl4 = Dlhs%lsao(DA)%startLocalOrb(1+(batchD-1)*maxAng)-1

       DB = Drhs%index(atomD,atomB,1,1)
       rDB => Drhs%lsao(DB)%elms
       CB = Drhs%index(atomC,atomB,1,1)
       rCB => Drhs%lsao(CB)%elms
       
       nr4 = Drhs%lsao(DB)%nLocal(1) 
       nr2 = Drhs%lsao(DB)%nLocal(2) 
       nr3 = Drhs%lsao(CB)%nLocal(1) 
       maxAng = Drhs%lsao(DB)%maxAng
       maxBat = Drhs%lsao(DB)%maxBat
       sr4 = Drhs%lsao(DB)%startLocalOrb(1+(batchD-1)*maxAng)-1
       sr2 = Drhs%lsao(DB)%startLocalOrb(1+(batchB-1)*maxAng+maxAng*maxBat)-1
       maxAng = Drhs%lsao(CB)%maxAng
       maxBat = Drhs%lsao(CB)%maxBat
       sr3 = Drhs%lsao(CB)%startLocalOrb(1+(batchC-1)*maxAng)-1

       call build_ltheta_FT(ltheta,lCA,rDB,lDA,rCB,&
            & nl1,nl3,nl4,sl1,sl3,sl4,nr2,nr3,nr4,sr2,sr3,sr4,&
            & nA,nB,nC,nD,ndmat)
    ELSE
       AC = Dlhs%index(atomA,atomC,1,1)
       lAC => Dlhs%lsao(AC)%elms
       
       nl1 = Dlhs%lsao(AC)%nLocal(1) 
       nl3 = Dlhs%lsao(AC)%nLocal(2) 
       maxAng = Dlhs%lsao(AC)%maxAng
       maxBat = Dlhs%lsao(AC)%maxBat
       sl1 = Dlhs%lsao(AC)%startLocalOrb(1+(batchA-1)*maxAng)-1
       sl3 = Dlhs%lsao(AC)%startLocalOrb(1+(batchC-1)*maxAng+maxAng*maxBat)-1

       BD = Drhs%index(atomB,atomD,1,1)
       rBD => Drhs%lsao(BD)%elms
       
       nr2 = Drhs%lsao(BD)%nLocal(1) 
       nr4 = Drhs%lsao(BD)%nLocal(2) 
       maxAng = Drhs%lsao(BD)%maxAng
       maxBat = Drhs%lsao(BD)%maxBat
       sr2 = Drhs%lsao(BD)%startLocalOrb(1+(batchB-1)*maxAng)-1
       sr4 = Drhs%lsao(BD)%startLocalOrb(1+(batchD-1)*maxAng+maxAng*maxBat)-1
       call build_ltheta_FF(ltheta,lAC,rBD,&
            & nl1,nl3,sl1,sl3,nr2,nr4,sr2,sr4,&
            & nA,nB,nC,nD,ndmat)
    ENDIF
    DO iDeriv=1,nDeriv

      tmp = gtheta(ltheta,CDAB,nA,nB,nC,nD,iPass,iDeriv,nPasses,nDeriv,lupri)

      iDer  = derivInfo%dirComp(1,iDeriv)
      iAtom = derivInfo%Atom(derivInfo%AO(1,iDeriv))
!The structure of the lstensor when used to store the molecular gradient
!is maybe not the most logical and we refer to discussion in the subroutine
!init_gradientlstensor in lsutil/lstensor_operations.f90
      iLSAO = Kgrad%INDEX(iAtom,1,1,permuteAO(derivInfo%AO(1,iDeriv)))
#ifdef VAR_LSDEBUGINT
      IF(iLSAO.EQ. 0)THEN
         WRITE(lupri,*)'iLSAO is zero - for some reason'
         WRITE(lupri,*)'iAtom          ',iAtom
         WRITE(lupri,*)'Center (1-4) = ',derivInfo%AO(1,iDeriv)
         WRITE(lupri,*)'Error in use of gradientlstensor distributeKgrad'
         CALL LSQUIT('Error in use of gradientlstensor distributeKgrad',lupri)
      ENDIF
#endif
      niderindex = niderindex + 1
      iLSAOindex(niderindex) = iLSAO
      iderindex(niderindex) = ider
      idertmp(niderindex) = tmp
!!$OMP ATOMIC
!      Kgrad%LSAO(iLSAO)%elms(iDer) = Kgrad%LSAO(iLSAO)%elms(iDer) - 0.5E0_realk*tmp
    ENDDO !iDeriv
    call mem_dealloc(ltheta)
  ENDDO !iPassQ
ENDDO !iPassP
!!$OMP CRITICAL (distributeKgradBlock)
do iDeriv=1,niderindex
 iLSAO = iLSAOindex(iDeriv)
 ider = iderindex(iDeriv)
 tmp = idertmp(iDeriv)
!$OMP ATOMIC
 Kgrad%LSAO(iLSAO)%elms(iDer) = Kgrad%LSAO(iLSAO)%elms(iDer) - 0.5E0_realk*tmp
enddo
!!$OMP END CRITICAL (distributeKgradBlock)

CALL freeDerivativeOverlapInfo(derivInfo)

CONTAINS
FUNCTION gtheta(ltheta,CDAB,nA,nB,nC,nD,iPass,iDeriv,nPasses,nDeriv,lupri)
implicit none
integer,intent(IN) :: nA,nB,nC,nD,iPass,iDeriv,nPasses,nDeriv,lupri
real(realk),intent(IN) :: ltheta(nC*nD*nA*nB)
Real(realk),intent(IN) :: CDAB(nC*nD*nA*nB,nDeriv,nPasses)
Real(realk)            :: gtheta
!
Integer     :: i
Real(realk) :: tmp
!
gtheta = 0E0_realk
DO i=1,nC*nD*nA*nB
  gtheta = gtheta + ltheta(i)*CDAB(i,iDeriv,iPass)
ENDDO

END FUNCTION gtheta

SUBROUTINE build_ltheta_TT(ltheta,lCA,rDB,lCB,rDA,lDA,rCB,lDB,rCA,&
     & nl1,nl2,nl3,nl4,sl1,sl2,sl3,sl4,nr1,nr2,nr3,nr4,sr1,sr2,sr3,sr4,nA,nB,nC,nD,ndmat)
implicit none
integer,intent(IN) :: nA,nB,nC,nD,ndmat
integer,intent(IN) :: nl1,nl2,nl3,nl4,sl1,sl2,sl3,sl4,nr1,nr2,nr3,nr4,sr1,sr2,sr3,sr4                      
real(realk) :: ltheta(nC,nD,nA,nB)
real(realk) :: lCA(nl3,nl1,ndmat)!(nC,nA,ndmat)
real(realk) :: lDA(nl4,nl1,ndmat)!(nD,nA,ndmat)
real(realk) :: lCB(nl3,nl2,ndmat)!(nC,nB,ndmat)
real(realk) :: lDB(nl4,nl2,ndmat)!(nD,nB,ndmat)
real(realk) :: rDB(nr4,nr2,ndmat)!(nD,nB,ndmat)
real(realk) :: rCB(nr3,nr2,ndmat)!(nC,nB,ndmat)
real(realk) :: rDA(nr4,nr1,ndmat)!(nD,nA,ndmat)
real(realk) :: rCA(nr3,nr1,ndmat)!(nC,nA,ndmat)
!
Integer :: idmat,iA,iB,iC,iD
!
real(realk),pointer :: ltheta_CADB(:,:,:,:)
real(realk),pointer :: ltheta_DACB(:,:,:,:)
real(realk),pointer :: ltheta_CBDA(:,:,:,:)
real(realk),pointer :: ltheta_DBCA(:,:,:,:)

IF (ndmat.GT.4) THEN
 call mem_alloc(ltheta_CADB,nC,nA,nD,nB)
 call mem_alloc(ltheta_DACB,nD,nA,nC,nB)
 call mem_alloc(ltheta_CBDA,nC,nB,nD,nA)
 call mem_alloc(ltheta_DBCA,nD,nB,nC,nA)
 call create_ltheta(ltheta_CADB,lCA,rDB,nl3,nl1,sl3,sl1,nr4,nr2,sr4,sr2,nC,nA,nD,nB,ndmat)
 call create_ltheta(ltheta_CBDA,lCB,rDA,nl3,nl2,sl3,sl2,nr4,nr1,sr4,sr1,nC,nB,nD,nA,ndmat)
 call create_ltheta(ltheta_DACB,lDA,rCB,nl4,nl1,sl4,sl1,nr3,nr2,sr3,sr2,nD,nA,nC,nB,ndmat)
 call create_ltheta(ltheta_DBCA,lDB,rCA,nl4,nl2,sl4,sl2,nr3,nr1,sr3,sr1,nD,nB,nC,nA,ndmat)
 DO iB=1,nB
  DO iD=1,nD
   DO iA=1,nA
    DO iC=1,nC
      ltheta(iC,iD,iA,iB) =   ltheta_CADB(iC,iA,iD,iB) + ltheta_DACB(iD,iA,iC,iB) &
      &                     + ltheta_CBDA(iC,iB,iD,iA) + ltheta_DBCA(iD,iB,iC,iA)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
 call mem_dealloc(ltheta_CADB)
 call mem_dealloc(ltheta_DACB)
 call mem_dealloc(ltheta_CBDA)
 call mem_dealloc(ltheta_DBCA)
ELSE
 DO idmat=1,ndmat
  DO iB=1,nB
   DO iD=1,nD
    DO iA=1,nA
     DO iC=1,nC
      ltheta(iC,iD,iA,iB) = ltheta(iC,iD,iA,iB) + lCA(sl3+iC,sl1+iA,idmat)*rDB(sr4+iD,sr2+iB,idmat) &
      & + lDA(sl4+iD,sl1+iA,idmat)*rCB(sr3+iC,sr2+iB,idmat) + lCB(sl3+iC,sl2+iB,idmat)*rDA(sr4+iD,sr1+iA,idmat) &
      & + lDB(sl4+iD,sl2+iB,idmat)*rCA(sr3+iC,sr1+iA,idmat)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDIF
END SUBROUTINE build_ltheta_TT

SUBROUTINE build_ltheta_TF(ltheta,lCA,rDB,lCB,rDA,&
     & nl1,nl2,nl3,sl1,sl2,sl3,nr1,nr2,nr4,sr1,sr2,sr4,nA,nB,nC,nD,ndmat)
implicit none
integer,intent(IN) :: nA,nB,nC,nD,ndmat
integer,intent(IN) :: nl1,nl2,nl3,sl1,sl2,sl3,nr1,nr2,nr4,sr1,sr2,sr4                      
real(realk) :: ltheta(nC,nD,nA,nB)
real(realk) :: lCA(nl3,nl1,ndmat)!(nC,nA,ndmat)
real(realk) :: lCB(nl3,nl2,ndmat)!(nC,nB,ndmat)
real(realk) :: rDB(nr4,nr2,ndmat)!(nD,nB,ndmat)
real(realk) :: rDA(nr4,nr1,ndmat)!(nD,nA,ndmat)
!
Integer :: idmat,iA,iB,iC,iD
real(realk),pointer :: ltheta_CADB(:,:,:,:)
real(realk),pointer :: ltheta_CBDA(:,:,:,:)

IF (ndmat.GT.4) THEN
 call mem_alloc(ltheta_CADB,nC,nA,nD,nB)
 call mem_alloc(ltheta_CBDA,nC,nB,nD,nA)
 call create_ltheta(ltheta_CADB,lCA,rDB,nl3,nl1,sl3,sl1,nr4,nr2,sr4,sr2,nC,nA,nD,nB,ndmat)
 call create_ltheta(ltheta_CBDA,lCB,rDA,nl3,nl2,sl3,sl2,nr4,nr1,sr4,sr1,nC,nB,nD,nA,ndmat)

 DO iB=1,nB
  DO iD=1,nD
   DO iA=1,nA
    DO iC=1,nC
      ltheta(iC,iD,iA,iB) =   ltheta_CADB(iC,iA,iD,iB) + ltheta_CBDA(iC,iB,iD,iA)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
 call mem_dealloc(ltheta_CADB)
 call mem_dealloc(ltheta_CBDA)
ELSE
 DO idmat=1,ndmat
  DO iB=1,nB
   DO iD=1,nD
    DO iA=1,nA
     DO iC=1,nC
      ltheta(iC,iD,iA,iB) = ltheta(iC,iD,iA,iB) + lCA(sl3+iC,sl1+iA,idmat)*rDB(sr4+iD,sr2+iB,idmat) &
                                              & + lCB(sl3+iC,sl2+iB,idmat)*rDA(sr4+iD,sr1+iA,idmat)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDIF
END SUBROUTINE build_ltheta_TF

SUBROUTINE build_ltheta_FT(ltheta,lCA,rDB,lDA,rCB,&
     & nl1,nl3,nl4,sl1,sl3,sl4,nr2,nr3,nr4,sr2,sr3,sr4,nA,nB,nC,nD,ndmat)
implicit none
integer,intent(IN) :: nA,nB,nC,nD,ndmat
integer,intent(IN) :: nl1,nl3,nl4,sl1,sl3,sl4,nr2,nr3,nr4,sr2,sr3,sr4                      
real(realk) :: ltheta(nC,nD,nA,nB)
real(realk) :: lCA(nl3,nl1,ndmat)!(nC,nA,ndmat)
real(realk) :: lDA(nl4,nl1,ndmat)!(nD,nA,ndmat)
real(realk) :: rDB(nr4,nr2,ndmat)!(nD,nB,ndmat)
real(realk) :: rCB(nr3,nr2,ndmat)!(nC,nB,ndmat)
!
Integer :: idmat,iA,iB,iC,iD
real(realk),pointer :: ltheta_CADB(:,:,:,:)
real(realk),pointer :: ltheta_DACB(:,:,:,:)

IF (ndmat.GT.4) THEN
 call mem_alloc(ltheta_CADB,nC,nA,nD,nB)
 call mem_alloc(ltheta_DACB,nD,nA,nC,nB)
 call create_ltheta(ltheta_CADB,lCA,rDB,nl3,nl1,sl3,sl1,nr4,nr2,sr4,sr2,nC,nA,nD,nB,ndmat)
 call create_ltheta(ltheta_DACB,lDA,rCB,nl4,nl1,sl4,sl1,nr3,nr2,sr3,sr2,nD,nA,nC,nB,ndmat)


 DO iB=1,nB
  DO iD=1,nD
   DO iA=1,nA
    DO iC=1,nC
      ltheta(iC,iD,iA,iB) =   ltheta_CADB(iC,iA,iD,iB) + ltheta_DACB(iD,iA,iC,iB)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
 call mem_dealloc(ltheta_CADB)
 call mem_dealloc(ltheta_DACB)
ELSE
 DO idmat=1,ndmat
  DO iB=1,nB
   DO iD=1,nD
    DO iA=1,nA
     DO iC=1,nC
      ltheta(iC,iD,iA,iB) = ltheta(iC,iD,iA,iB) + lCA(sl3+iC,sl1+iA,idmat)*rDB(sr4+iD,sr2+iB,idmat) &
           & + lDA(sl4+iD,sl1+iA,idmat)*rCB(sr3+iC,sr2+iB,idmat)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDIF
END SUBROUTINE build_ltheta_FT

SUBROUTINE build_ltheta_FF(ltheta,lAC,rBD,&
     & nl1,nl3,sl1,sl3,nr2,nr4,sr2,sr4,nA,nB,nC,nD,ndmat)
implicit none
integer,intent(IN) :: nA,nB,nC,nD,ndmat
integer,intent(IN) :: nl1,nl3,sl1,sl3,nr2,nr4,sr2,sr4                      
real(realk) :: ltheta(nC,nD,nA,nB)
real(realk) :: lAC(nl1,nl3,ndmat)!(nA,nC,ndmat)
real(realk) :: rBD(nr2,nr4,ndmat)!(nB,nD,ndmat)
!
Integer :: idmat,iA,iB,iC,iD
real(realk),pointer :: ltheta_ACBD(:,:,:,:)
real(realk),pointer :: ltheta_DACB(:,:,:,:)

IF (ndmat.GT.4) THEN
 call mem_alloc(ltheta_ACBD,nA,nC,nB,nD)
 call create_ltheta(ltheta_ACBD,lAC,rBD,nl1,nl3,sl1,sl3,nr2,nr4,sr2,sr4,nA,nC,nB,nD,ndmat)

 DO iB=1,nB
  DO iD=1,nD
   DO iA=1,nA
    DO iC=1,nC
      ltheta(iC,iD,iA,iB) =   ltheta_ACBD(iA,iC,iB,iD)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
 call mem_dealloc(ltheta_ACBD)
ELSE
 DO idmat=1,ndmat
  DO iB=1,nB
   DO iD=1,nD
    DO iA=1,nA
     DO iC=1,nC
      ltheta(iC,iD,iA,iB) = ltheta(iC,iD,iA,iB) + lAC(sl1+iA,sl3+iC,idmat)*rBD(sr2+iB,sr4+iD,idmat)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDIF
END SUBROUTINE build_ltheta_FF

SUBROUTINE create_ltheta(ltheta_PQ,P,Q,n1,n2,sA,sB,n3,n4,sC,sD,nA,nB,nC,nD,ndmat)
implicit none
integer,intent(IN)        :: nA,nB,nC,nD,n1,n2,n3,n4,sA,sB,sC,sD,ndmat
real(realk),intent(INOUT) :: ltheta_PQ(nA,nB,nC,nD)
real(realk),intent(IN)    :: P(n1,n2,ndmat)
real(realk),intent(IN)    :: Q(n3,n4,ndmat)
!
integer :: iA,iB,iC,iD,idmat
real(realk) :: TMP
ltheta_PQ = 0E0_realk
DO idmat=1,ndmat
 DO iD=1,nD
  DO iC=1,nC
   TMP = Q(sC+iC,sD+iD,idmat)
   DO iB=1,nB
    DO iA=1,nA
     ltheta_PQ(iA,iB,iC,iD) = ltheta_PQ(iA,iB,iC,iD) + P(sA+iA,sB+iB,idmat)*TMP
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO
END SUBROUTINE create_ltheta

END SUBROUTINE distributeKgrad

!> \brief Distribute exchange-type integrals
!> \author \latexonly S. Reine  \endlatexonly
!> \date 2011-03-08
!> \param Kmat contains the result lstensor
!> \param CDAB matrix containing calculated integrals
!> \param Dmat contains the rhs-density matrix/matrices
!> \param ndmat number of rhs-density matries
!> \param PQ information about the calculated integrals (to be distributed)
!> \param Input contain info about the requested integral 
!> \param Lsoutput contain info about the requested output, and actual output
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE distributeKcont(CDAB,Dlhs,Drhs,ndmat,PQ,Input,&
    &                      Lsoutput,Kcont,LUPRI,IPRINT)
implicit none 
Type(integrand),intent(in)         :: PQ
Type(IntegralInput),intent(in)     :: Input
Type(IntegralOutput),intent(inout) :: Lsoutput
Integer,intent(in)                 :: ndmat,LUPRI,IPRINT
REAL(REALK),pointer                :: CDAB(:)
TYPE(lstensor),intent(in)          :: Dlhs,Drhs
Real(Realk),intent(inout)          :: Kcont(ndmat)
!
type(overlap),pointer :: P,Q
logical :: SamePQ,SameLHSaos,SameRHSaos,SameODs
Logical :: permuteAB,permuteCD,permuteOD
integer :: nAngmomP,nPasses,maxBat,maxAng
integer :: iPass,iPassP,iPassQ,idmat,DAr,CBr
integer :: nA,nB,nC,nD,batchA,batchB,batchC,batchD,atomA,atomB,atomC,atomD
integer :: AC,BD,DA,CB,AD,BC,n1,n2,n3,n4,s1,s2,s3,s4
integer :: CA,DB
Real(realk),pointer :: lAC(:),lBC(:),lAD(:),lBD(:),rBD(:),rAD(:),rBC(:),rAC(:)
Real(realk),pointer :: lCA(:),lCB(:),lDA(:),lDB(:),rDB(:),rDA(:),rCB(:),rCA(:)
Real(realk) :: tmpKcont(ndmat)
IF(.NOT.INPUT%LHSSameAsRHSDmat)THEN
   !We assume sym dmat

   P => PQ%P%p
   Q => PQ%Q%p
   nAngmomP   = P%nAngmom
   SamePQ     = PQ%samePQ
   SameLHSaos = INPUT%SameLHSaos
   SameRHSaos = INPUT%SameRHSaos
   SameODs    = INPUT%SameODs
   nPasses = P%nPasses*Q%nPasses
   iPass = 0
   tmpKcont = 0E0_realk
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
    DO iPassQ=1,Q%nPasses
     iPass = iPass + 1
         
     atomC  = Q%orb1atom(iPassQ)
     batchC = Q%orb1batch(iPassQ)
     atomD  = Q%orb2atom(iPassQ)
     batchD = Q%orb2batch(iPassQ)
     permuteCD = SameRHSaos.AND.((batchC.NE.batchD).OR.(atomC.NE.atomD))
     permuteOD = SameODs.AND.( ((batchA.NE.batchC).OR.(atomA.NE.atomC)).OR.&
          &                          ((batchB.NE.batchD).OR.(atomB.NE.atomD)) )
     permuteOD = SameODs.AND..NOT.((batchA.EQ.batchC).AND.(atomA.EQ.atomC).AND.&
          &                          (batchB.EQ.batchD).AND.(atomB.EQ.atomD) )
     
     
     AC = Dlhs%index(atomA,atomC,1,1)
     lAC => Dlhs%lsao(AC)%elms       
     n1 = Dlhs%lsao(AC)%nLocal(1) 
     n3 = Dlhs%lsao(AC)%nLocal(2) 
     maxAng = Dlhs%lsao(AC)%maxAng
     maxBat = Dlhs%lsao(AC)%maxBat
     s1 = Dlhs%lsao(AC)%startLocalOrb(1+(batchA-1)*maxAng)-1
     s3 = Dlhs%lsao(AC)%startLocalOrb(1+(batchC-1)*maxAng+maxAng*maxBat)-1
     
     BD = Drhs%index(atomB,atomD,1,1)
     rBD => Drhs%lsao(BD)%elms
     
     n2 = Drhs%lsao(BD)%nLocal(1) 
     n4 = Drhs%lsao(BD)%nLocal(2) 
     maxAng = Drhs%lsao(BD)%maxAng
     maxBat = Drhs%lsao(BD)%maxBat
     s2 = Drhs%lsao(BD)%startLocalOrb(1+(batchB-1)*maxAng)-1
     s4 = Drhs%lsao(BD)%startLocalOrb(1+(batchD-1)*maxAng+maxAng*maxBat)-1
     
     IF (permuteCD) THEN
        AD = Dlhs%index(atomA,atomD,1,1)
        lAD => Dlhs%lsao(AD)%elms       
        BC = Drhs%index(atomB,atomC,1,1)
        rBC => Drhs%lsao(BC)%elms
     ENDIF
     IF (permuteAB) THEN
        BD = Dlhs%index(atomB,atomD,1,1)
        lBD => Dlhs%lsao(BD)%elms       
        AC = Drhs%index(atomA,atomC,1,1)
        rAC => Drhs%lsao(AC)%elms
        
        AD = Drhs%index(atomA,atomD,1,1)
        rAD => Drhs%lsao(AD)%elms       
        BC = Dlhs%index(atomB,atomC,1,1)
        lBC => Dlhs%lsao(BC)%elms
     ENDIF
     IF (permuteAB.AND.permuteCD) THEN
        call fullEcont2(tmpKcont,lAC,lBC,lAD,lBD,rBD,rAD,rBC,rAC,n1,n2,n3,n4,s1,s2,s3,s4,&
             & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,permuteOD,lupri)
     ELSE IF (permuteAB) THEN
        call fullEcont3(tmpKcont,lAC,lBC,rBD,rAD,n1,n2,n3,n4,s1,s2,s3,s4,&
             & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,permuteOD,lupri)
     ELSE IF (permuteCD) THEN       
        call fullEcont4(tmpKcont,lAC,lAD,rBD,rBC,n1,n2,n3,n4,s1,s2,s3,s4,&
             & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,permuteOD,lupri)
     ELSE
        call fullEcont5(tmpKcont,lAC,rBD,n1,n2,n3,n4,s1,s2,s3,s4,&
             & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,permuteOD,lupri)
     ENDIF
    ENDDO !iPassQ
   ENDDO !iPassP

   !reduction afterwards
!!$OMP CRITICAL (distributeKgradBlock)
   DO idmat=1,ndmat
      Kcont(idmat) = Kcont(idmat) - tmpKcont(idmat)
   ENDDO
!!$OMP END CRITICAL (distributeKgradBlock)
ELSE

   !Special for same lhs and rhs !!!!
   !assumes symmetric matrices
   
   P => PQ%P%p
   Q => PQ%Q%p
   nAngmomP   = P%nAngmom
   SamePQ     = PQ%samePQ
   SameLHSaos = INPUT%SameLHSaos
   SameRHSaos = INPUT%SameRHSaos
   SameODs    = INPUT%SameODs
   nPasses = P%nPasses*Q%nPasses
   iPass = 0
   tmpKcont = 0E0_realk
   DO iPassP=1,P%nPasses
    atomA  = P%orb1atom(iPassP)
    batchA = P%orb1batch(iPassP)
    nA     = P%orbital1%totOrbitals
    atomB  = P%orb2atom(iPassP)
    batchB = P%orb2batch(iPassP)
    nB     = P%orbital2%totOrbitals    
    permuteAB  = SameLHSaos.AND.((batchA.NE.batchB).OR.(atomA.NE.atomB))
    DO iPassQ=1,Q%nPasses
     iPass = iPass + 1

     atomC  = Q%orb1atom(iPassQ)
     batchC = Q%orb1batch(iPassQ)
     nC     = Q%orbital1%totOrbitals
     atomD  = Q%orb2atom(iPassQ)
     batchD = Q%orb2batch(iPassQ)
     nD     = Q%orbital2%totOrbitals
     permuteCD = SameRHSaos.AND.((batchC.NE.batchD).OR.(atomC.NE.atomD))
     permuteOD = SameODs.AND.( ((batchA.NE.batchC).OR.(atomA.NE.atomC)).OR.&
          &                          ((batchB.NE.batchD).OR.(atomB.NE.atomD)) )
     permuteOD = SameODs.AND..NOT.((batchA.EQ.batchC).AND.(atomA.EQ.atomC).AND.&
          &                          (batchB.EQ.batchD).AND.(atomB.EQ.atomD) )
     
     CA = Dlhs%index(atomC,atomA,1,1)
     lCA => Dlhs%lsao(CA)%elms
     DB = Dlhs%index(atomD,atomB,1,1)
     lDB => Dlhs%lsao(DB)%elms       
     n3 = Dlhs%lsao(CA)%nLocal(1) 
     n1 = Dlhs%lsao(CA)%nLocal(2) 
     n4 = Dlhs%lsao(DB)%nLocal(1) 
     n2 = Dlhs%lsao(DB)%nLocal(2) 
     maxBat = Dlhs%lsao(CA)%maxBat
     maxAng = Dlhs%lsao(CA)%maxAng
     s3 = Dlhs%lsao(CA)%startLocalOrb(1+(batchC-1)*maxAng)-1
     s1 = Dlhs%lsao(CA)%startLocalOrb(1+(batchA-1)*maxAng+maxAng*maxBat)-1
     maxBat = Dlhs%lsao(DB)%maxBat
     maxAng = Dlhs%lsao(DB)%maxAng
     s4 = Dlhs%lsao(DB)%startLocalOrb(1+(batchD-1)*maxAng)-1
     s2 = Dlhs%lsao(DB)%startLocalOrb(1+(batchB-1)*maxAng+maxAng*maxBat)-1
     IF (permuteOD.AND.permuteCD.AND.permuteAB) THEN
        DA = Dlhs%index(atomD,atomA,1,1)
        lDA => Dlhs%lsao(DA)%elms
        CB = Dlhs%index(atomC,atomB,1,1)
        lCB => Dlhs%lsao(CB)%elms
        call full1(tmpKcont,lCA,lDA,lCB,lDB,n1,n2,n3,n4,s1,s2,s3,s4, &
             & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,lupri)
     ELSE IF (permuteOD.AND.permuteCD) THEN
        DA = Dlhs%index(atomD,atomA,1,1)
        lDA => Dlhs%lsao(DA)%elms
        CB = Dlhs%index(atomC,atomB,1,1)
        lCB => Dlhs%lsao(CB)%elms
        call full2(tmpKcont,lCA,lDA,lCB,lDB,n1,n2,n3,n4,s1,s2,s3,s4,&
             & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,lupri)
     ELSE IF (permuteOD.AND.permuteAB) THEN
        DA = Dlhs%index(atomD,atomA,1,1)
        lDA => Dlhs%lsao(DA)%elms
        CB = Dlhs%index(atomC,atomB,1,1)
        lCB => Dlhs%lsao(CB)%elms
        call full2(tmpKcont,lCA,lDA,lCB,lDB,n1,n2,n3,n4,s1,s2,s3,s4,&
             & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,lupri)
     ELSE IF (permuteAB.AND.permuteCD) THEN
        DA = Dlhs%index(atomD,atomA,1,1)
        lDA => Dlhs%lsao(DA)%elms
        CB = Dlhs%index(atomC,atomB,1,1)
        lCB => Dlhs%lsao(CB)%elms
        call full2(tmpKcont,lCA,lDA,lCB,lDB,n1,n2,n3,n4,s1,s2,s3,s4,&
             & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,lupri)
     ELSE IF (permuteAB) THEN
        DA = Dlhs%index(atomD,atomA,1,1)
        lDA => Dlhs%lsao(DA)%elms
        CB = Dlhs%index(atomC,atomB,1,1)
        lCB => Dlhs%lsao(CB)%elms
        call full3(tmpKcont,lCA,lDA,lCB,lDB,n1,n2,n3,n4,s1,s2,s3,s4,&
             & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,lupri)
     ELSE IF (permuteCD) THEN
        DA = Dlhs%index(atomD,atomA,1,1)
        lDA => Dlhs%lsao(DA)%elms
        CB = Dlhs%index(atomC,atomB,1,1)
        lCB => Dlhs%lsao(CB)%elms
        call full3(tmpKcont,lCA,lDA,lCB,lDB,n1,n2,n3,n4,s1,s2,s3,s4,&
             & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,lupri)
     ELSE IF (permuteOD) THEN
        call full4(tmpKcont,lCA,lDB,n1,n2,n3,n4,s1,s2,s3,s4,&
             & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,lupri)
     ELSE
        call full5(tmpKcont,lCA,lDB,n1,n2,n3,n4,s1,s2,s3,s4,&
             & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,lupri)
     ENDIF
    ENDDO !iPassQ
   ENDDO !iPassP
   !not necessary ! 
   !reduction afterwards
!!$OMP CRITICAL (distributeKgradBlock)
   DO idmat=1,ndmat
      Kcont(idmat) = Kcont(idmat) - tmpKcont(idmat)
   ENDDO
!!$OMP END CRITICAL (distributeKgradBlock)
ENDIF
   
CONTAINS
subroutine full1(Kcont,DCA,DDA,DCB,DDB,n1,n2,n3,n4,s1,s2,s3,s4,&
     & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,lupri)
implicit none
Integer,intent(IN)     :: nA,nB,nC,nD,iPass,nPasses,lupri,ndmat
Integer,intent(IN)     :: n1,n2,n3,n4,s1,s2,s3,s4
Real(realk),intent(IN) :: CDAB(nC,nD,nA,nB,nPasses)
Real(realk),intent(IN) :: DCA(n3,n1,ndmat),DDA(n4,n1,ndmat),DCB(n3,n2,ndmat),DDB(n4,n2,ndmat)
Real(realk),intent(INOUT) :: Kcont(ndmat)
!
Integer     :: iA,iB,iC,iD,idmat
Real(realk) :: tmpDB,tmpDA
real(realk),parameter :: D4=4.0E0_realk,D2=2.0E0_realk,D0=0.0E0_realk
!
DO idmat=1,ndmat
 DO iB=1,nB
  DO iA=1,nA
   DO iD=1,nD
    tmpDA = D0
    tmpDB = D0
    DO iC=1,nC
     tmpDB = tmpDB + DCA(s3+iC,s1+iA,idmat) * CDAB(iC,iD,iA,iB,iPass)
     tmpDA = tmpDA + DCB(s3+iC,s2+iB,idmat) * CDAB(iC,iD,iA,iB,iPass)
    ENDDO
    Kcont(idmat) = Kcont(idmat) + D4 * tmpDA * DDA(s4+iD,s1+iA,idmat)
    Kcont(idmat) = Kcont(idmat) + D4 * tmpDB * DDB(s4+iD,s2+iB,idmat)
   ENDDO
  ENDDO
 ENDDO
ENDDO
end subroutine full1

subroutine full2(Kcont,DCA,DDA,DCB,DDB,n1,n2,n3,n4,s1,s2,s3,s4,&
     & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,lupri)
implicit none
Integer,intent(IN)     :: nA,nB,nC,nD,iPass,nPasses,lupri,ndmat
Integer,intent(IN)     :: n1,n2,n3,n4,s1,s2,s3,s4
Real(realk),intent(IN) :: CDAB(nC,nD,nA,nB,nPasses)
Real(realk),intent(IN) :: DCA(n3,n1,ndmat),DDA(n4,n1,ndmat),DCB(n3,n2,ndmat),DDB(n4,n2,ndmat)
Real(realk),intent(INOUT) :: Kcont(ndmat)
!
Integer     :: iA,iB,iC,iD,idmat
Real(realk) :: tmpDB,tmpDA
real(realk),parameter :: D4=4.0E0_realk,D2=2.0E0_realk,D0=0.0E0_realk
!
DO idmat=1,ndmat
 DO iB=1,nB
  DO iA=1,nA
   DO iD=1,nD
    tmpDA = D0
    tmpDB = D0
    DO iC=1,nC
     tmpDB = tmpDB + DCA(s3+iC,s1+iA,idmat) * CDAB(iC,iD,iA,iB,iPass)
     tmpDA = tmpDA + DCB(s3+iC,s2+iB,idmat) * CDAB(iC,iD,iA,iB,iPass)
    ENDDO
    Kcont(idmat) = Kcont(idmat) + D2 * tmpDB * DDB(s4+iD,s2+iB,idmat)
    Kcont(idmat) = Kcont(idmat) + D2 * tmpDA * DDA(s4+iD,s1+iA,idmat)
   ENDDO
  ENDDO
 ENDDO
ENDDO
end subroutine full2

subroutine full3(Kcont,DCA,DDA,DCB,DDB,n1,n2,n3,n4,s1,s2,s3,s4,&
     & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,lupri)
implicit none
Integer,intent(IN)     :: nA,nB,nC,nD,iPass,nPasses,lupri,ndmat
Integer,intent(IN)     :: n1,n2,n3,n4,s1,s2,s3,s4
Real(realk),intent(IN) :: CDAB(nC,nD,nA,nB,nPasses)
Real(realk),intent(IN) :: DCA(n3,n1,ndmat),DDA(n4,n1,ndmat),DCB(n3,n2,ndmat),DDB(n4,n2,ndmat)
Real(realk),intent(INOUT) :: Kcont(ndmat)
!
Integer     :: iA,iB,iC,iD,idmat
Real(realk) :: tmpDB,tmpDA
real(realk),parameter :: D4=4.0E0_realk,D2=2.0E0_realk,D0=0.0E0_realk
!
DO idmat=1,ndmat
 DO iB=1,nB
  DO iA=1,nA
   DO iD=1,nD
    tmpDA = D0
    tmpDB = D0
    DO iC=1,nC
     tmpDB = tmpDB + DCA(s3+iC,s1+iA,idmat) * CDAB(iC,iD,iA,iB,iPass)
     tmpDA = tmpDA + DCB(s3+iC,s2+iB,idmat) * CDAB(iC,iD,iA,iB,iPass)
    ENDDO
    Kcont(idmat) = Kcont(idmat) + tmpDB * DDB(s4+iD,s2+iB,idmat)
    Kcont(idmat) = Kcont(idmat) + tmpDA * DDA(s4+iD,s1+iA,idmat)
   ENDDO
  ENDDO
 ENDDO
ENDDO
end subroutine full3

subroutine full4(Kcont,DCA,DDB,n1,n2,n3,n4,s1,s2,s3,s4,&
     & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,lupri)
implicit none
Integer,intent(IN)     :: nA,nB,nC,nD,iPass,nPasses,lupri,ndmat
Integer,intent(IN)     :: n1,n2,n3,n4,s1,s2,s3,s4
Real(realk),intent(IN) :: CDAB(nC,nD,nA,nB,nPasses)
Real(realk),intent(IN) :: DCA(n3,n1,ndmat),DDB(n4,n2,ndmat)
Real(realk),intent(INOUT) :: Kcont(ndmat)
!
Integer     :: iA,iB,iC,iD,idmat
Real(realk) :: tmpDB,tmpDA
real(realk),parameter :: D2=2.0E0_realk,D0=0.0E0_realk
!
DO idmat=1,ndmat
 DO iB=1,nB
  DO iA=1,nA
   DO iD=1,nD
    tmpDA = D0
    tmpDB = D0
    DO iC=1,nC
     tmpDB = tmpDB + DCA(s3+iC,s1+iA,idmat) * CDAB(iC,iD,iA,iB,iPass)
    ENDDO
    Kcont(idmat) = Kcont(idmat) + D2 * tmpDB * DDB(s4+iD,s2+iB,idmat)
   ENDDO
  ENDDO
 ENDDO
ENDDO
end subroutine full4

subroutine full5(Kcont,DCA,DDB,n1,n2,n3,n4,s1,s2,s3,s4,&
     & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,lupri)
implicit none
Integer,intent(IN)     :: nA,nB,nC,nD,iPass,nPasses,lupri,ndmat
Integer,intent(IN)     :: n1,n2,n3,n4,s1,s2,s3,s4
Real(realk),intent(IN) :: CDAB(nC,nD,nA,nB,nPasses)
Real(realk),intent(IN) :: DCA(n3,n1,ndmat),DDB(n4,n2,ndmat)
Real(realk),intent(INOUT) :: Kcont(ndmat)
!
Integer     :: iA,iB,iC,iD,idmat
Real(realk) :: tmpDB(ndmat)
real(realk),parameter :: D0=0.0E0_realk
!
IF(ndmat.EQ.10)THEN
 DO iB=1,nB
  DO iA=1,nA
   DO iD=1,nD
     tmpDB = D0
     DO iC=1,nC
        tmpDB(1) = tmpDB(1) + DCA(s3+iC,s1+iA,1) * CDAB(iC,iD,iA,iB,iPass)
        tmpDB(2) = tmpDB(2) + DCA(s3+iC,s1+iA,2) * CDAB(iC,iD,iA,iB,iPass)
        tmpDB(3) = tmpDB(3) + DCA(s3+iC,s1+iA,3) * CDAB(iC,iD,iA,iB,iPass)
        tmpDB(4) = tmpDB(4) + DCA(s3+iC,s1+iA,4) * CDAB(iC,iD,iA,iB,iPass)
        tmpDB(5) = tmpDB(5) + DCA(s3+iC,s1+iA,5) * CDAB(iC,iD,iA,iB,iPass)
        tmpDB(6) = tmpDB(6) + DCA(s3+iC,s1+iA,6) * CDAB(iC,iD,iA,iB,iPass)
        tmpDB(7) = tmpDB(7) + DCA(s3+iC,s1+iA,7) * CDAB(iC,iD,iA,iB,iPass)
        tmpDB(8) = tmpDB(8) + DCA(s3+iC,s1+iA,8) * CDAB(iC,iD,iA,iB,iPass)
        tmpDB(9) = tmpDB(9) + DCA(s3+iC,s1+iA,9) * CDAB(iC,iD,iA,iB,iPass)
        tmpDB(10) = tmpDB(10) + DCA(s3+iC,s1+iA,10) * CDAB(iC,iD,iA,iB,iPass)
     ENDDO
     DO idmat=1,10
        Kcont(idmat) = Kcont(idmat) + tmpDB(idmat) * DDB(s4+iD,s2+iB,idmat)
     ENDDO
   ENDDO
  ENDDO
 ENDDO
ELSE
 DO iB=1,nB
  DO iA=1,nA
   DO iD=1,nD
     tmpDB = D0
     DO iC=1,nC
        DO idmat=1,ndmat
           tmpDB(idmat) = tmpDB(idmat) + DCA(s3+iC,s1+iA,idmat) * CDAB(iC,iD,iA,iB,iPass)
        ENDDO
     ENDDO
     DO idmat=1,ndmat
        Kcont(idmat) = Kcont(idmat) + tmpDB(idmat) * DDB(s4+iD,s2+iB,idmat)
     ENDDO    
   ENDDO
  ENDDO
 ENDDO
ENDIF
end subroutine full5

subroutine fullEcont2(Kcont,lAC,lBC,lAD,lBD,rBD,rAD,rBC,rAC,&
     & n1,n2,n3,n4,s1,s2,s3,s4,CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,sym,lupri)
implicit none
Integer,intent(IN)     :: nA,nB,nC,nD,iPass,nPasses,lupri,ndmat
Integer,intent(IN)     :: n1,n2,n3,n4,s1,s2,s3,s4
Real(realk),intent(IN) :: CDAB(nC,nD,nA,nB,nPasses)
Real(realk),intent(IN) :: lAC(n1,n3,ndmat),rBD(n2,n4,ndmat)
Real(realk),intent(IN) :: lBC(n2,n3,ndmat),rAD(n1,n4,ndmat)
Real(realk),intent(IN) :: lAD(n1,n4,ndmat),rBC(n2,n3,ndmat)
Real(realk),intent(IN) :: lBD(n2,n4,ndmat),rAC(n1,n3,ndmat)
Real(realk),intent(INOUT) :: Kcont(ndmat)
logical,intent(IN)     :: sym
!
Integer     :: iA,iB,iC,iD,idmat
Real(realk) :: tmpBD(ndmat),tmpAD(ndmat),tmpBC(ndmat),tmpAC(ndmat)
real(realk),parameter :: D0=0.0E0_realk,D2=2.0E0_realk
!
IF(.NOT.sym)THEN
 DO iB=1,nB
  DO iA=1,nA
   DO iD=1,nD
    tmpBD = D0
    tmpAD = D0
    tmpBC = D0
    tmpAC = D0
    DO iC=1,nC
     DO idmat=1,ndmat
       tmpBD(idmat) = tmpBD(idmat) + lAC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
       tmpAD(idmat) = tmpAD(idmat) + lBC(s2+iB,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
       tmpBC(idmat) = tmpBC(idmat) + rBC(s2+iB,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
       tmpAC(idmat) = tmpAC(idmat) + rAC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
     ENDDO
    ENDDO
    DO idmat=1,ndmat
     Kcont(idmat) = Kcont(idmat) + tmpBD(idmat) * rBD(s2+iB,s4+iD,idmat)
     Kcont(idmat) = Kcont(idmat) + tmpAD(idmat) * rAD(s1+iA,s4+iD,idmat)
     Kcont(idmat) = Kcont(idmat) + tmpBC(idmat) * lAD(s1+iA,s4+iD,idmat)
     Kcont(idmat) = Kcont(idmat) + tmpAC(idmat) * lBD(s2+iB,s4+iD,idmat)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ELSE
 DO iB=1,nB
  DO iA=1,nA
   DO iD=1,nD
    tmpBD = D0
    tmpAD = D0
    tmpBC = D0
    tmpAC = D0
    DO iC=1,nC
     DO idmat=1,ndmat
       tmpBD(idmat) = tmpBD(idmat) + lAC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
       tmpAD(idmat) = tmpAD(idmat) + lBC(s2+iB,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
       tmpBC(idmat) = tmpBC(idmat) + rBC(s2+iB,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
       tmpAC(idmat) = tmpAC(idmat) + rAC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
     ENDDO
    ENDDO
    DO idmat=1,ndmat
     Kcont(idmat) = Kcont(idmat) + D2 * tmpBD(idmat) * rBD(s2+iB,s4+iD,idmat)
     Kcont(idmat) = Kcont(idmat) + D2 * tmpAD(idmat) * rAD(s1+iA,s4+iD,idmat)
     Kcont(idmat) = Kcont(idmat) + D2 * tmpBC(idmat) * lAD(s1+iA,s4+iD,idmat)
     Kcont(idmat) = Kcont(idmat) + D2 * tmpAC(idmat) * lBD(s2+iB,s4+iD,idmat)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDIF
end subroutine fullEcont2

subroutine fullEcont3(Kcont,DAC,DBC,DBD,DAD,n1,n2,n3,n4,s1,s2,s3,s4,&
     & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,sym,lupri)
implicit none
Integer,intent(IN)     :: nA,nB,nC,nD,iPass,nPasses,lupri,ndmat
Integer,intent(IN)     :: n1,n2,n3,n4,s1,s2,s3,s4
Real(realk),intent(IN) :: CDAB(nC,nD,nA,nB,nPasses)
Real(realk),intent(IN) :: DAC(n1,n3,ndmat),DBD(n2,n4,ndmat)
Real(realk),intent(IN) :: DBC(n2,n3,ndmat),DAD(n1,n4,ndmat)
Real(realk),intent(INOUT) :: Kcont(ndmat)
logical,intent(IN)     :: sym
!
Integer     :: iA,iB,iC,iD,idmat
Real(realk) :: tmpBD(ndmat),tmpAD(ndmat)
real(realk),parameter :: D0=0.0E0_realk,D2=2.0E0_realk
!
IF(.NOT.sym)THEN
 DO iB=1,nB
  DO iA=1,nA
   DO iD=1,nD
    tmpBD = D0
    tmpAD = D0
    DO iC=1,nC
     DO idmat=1,ndmat
        tmpBD(idmat) = tmpBD(idmat) + DAC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
        tmpAD(idmat) = tmpAD(idmat) + DBC(s2+iB,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
     ENDDO
    ENDDO
    DO idmat=1,ndmat
     Kcont(idmat) = Kcont(idmat) + tmpBD(idmat) * DBD(s2+iB,s4+iD,idmat)
     Kcont(idmat) = Kcont(idmat) + tmpAD(idmat) * DAD(s1+iA,s4+iD,idmat)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ELSE
 DO iB=1,nB
  DO iA=1,nA
   DO iD=1,nD
    tmpBD = D0
    tmpAD = D0
    DO iC=1,nC
     DO idmat=1,ndmat
        tmpBD(idmat) = tmpBD(idmat) + DAC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
        tmpAD(idmat) = tmpAD(idmat) + DBC(s2+iB,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
     ENDDO
    ENDDO
    DO idmat=1,ndmat
     Kcont(idmat) = Kcont(idmat) + D2 * tmpBD(idmat) * DBD(s2+iB,s4+iD,idmat)
     Kcont(idmat) = Kcont(idmat) + D2 * tmpAD(idmat) * DAD(s1+iA,s4+iD,idmat)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDIF
end subroutine fullEcont3

subroutine fullEcont4(Kcont,DAC,DAD,DBD,DBC,n1,n2,n3,n4,s1,s2,s3,s4,&
     & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,sym,lupri)
implicit none
Integer,intent(IN)     :: nA,nB,nC,nD,iPass,nPasses,lupri,ndmat
Integer,intent(IN)     :: n1,n2,n3,n4,s1,s2,s3,s4
Real(realk),intent(IN) :: CDAB(nC,nD,nA,nB,nPasses)
Real(realk),intent(IN) :: DAC(n1,n3,ndmat),DBD(n2,n4,ndmat)
Real(realk),intent(IN) :: DAD(n1,n4,ndmat),DBC(n2,n3,ndmat)
Real(realk),intent(INOUT) :: Kcont(ndmat)
logical,intent(IN)     :: sym
!
Integer     :: iA,iB,iC,iD,idmat
Real(realk) :: tmpBD(ndmat),tmpBC(ndmat)
real(realk),parameter :: D0=0.0E0_realk,D2=2.0E0_realk
!
IF(.NOT.sym)THEN
 DO iB=1,nB
  DO iA=1,nA
   DO iD=1,nD
    tmpBD = D0
    tmpBC = D0
    DO iC=1,nC
     DO idmat=1,ndmat
       tmpBD(idmat) = tmpBD(idmat) + DAC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
       tmpBC(idmat) = tmpBC(idmat) + DBC(s2+iB,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
     ENDDO
    ENDDO
    DO idmat=1,ndmat
     Kcont(idmat) = Kcont(idmat) + tmpBD(idmat) * DBD(s2+iB,s4+iD,idmat)
     Kcont(idmat) = Kcont(idmat) + tmpBC(idmat) * DAD(s1+iA,s4+iD,idmat)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ELSE
 DO iB=1,nB
  DO iA=1,nA
   DO iD=1,nD
    tmpBD = D0
    tmpBC = D0
    DO iC=1,nC
     DO idmat=1,ndmat
       tmpBD(idmat) = tmpBD(idmat) + DAC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
       tmpBC(idmat) = tmpBC(idmat) + DBC(s2+iB,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
     ENDDO
    ENDDO
    DO idmat=1,ndmat
     Kcont(idmat) = Kcont(idmat) + D2 * tmpBD(idmat) * DBD(s2+iB,s4+iD,idmat)
     Kcont(idmat) = Kcont(idmat) + D2 * tmpBC(idmat) * DAD(s1+iA,s4+iD,idmat)
    ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDIF
end subroutine fullEcont4

subroutine fullEcont5(Kcont,DAC,DBD,n1,n2,n3,n4,s1,s2,s3,s4,&
     & CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,sym,lupri)
implicit none
Integer,intent(IN)     :: nA,nB,nC,nD,iPass,nPasses,lupri,ndmat
Integer,intent(IN)     :: n1,n2,n3,n4,s1,s2,s3,s4
Real(realk),intent(IN) :: CDAB(nC,nD,nA,nB,nPasses)
Real(realk),intent(IN) :: DAC(n1,n3,ndmat),DBD(n2,n4,ndmat)
Real(realk),intent(INOUT) :: Kcont(ndmat)
logical,intent(IN)     :: sym
!
Integer     :: iA,iB,iC,iD,idmat
Real(realk) :: tmpBD(ndmat)
real(realk),parameter :: D0=0.0E0_realk,D2=2.0E0_realk
!
IF(ndmat.EQ.10)THEN
 IF(.NOT.sym)THEN
  DO iB=1,nB
   DO iA=1,nA
    DO iD=1,nD
     tmpBD = D0
     DO iC=1,nC
        DO idmat=1,10
           tmpBD(idmat) = tmpBD(idmat) + DAC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
        ENDDO
     ENDDO
     DO idmat=1,10
        Kcont(idmat) = Kcont(idmat) + tmpBD(idmat) * DBD(s2+iB,s4+iD,idmat)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ELSE
  DO iB=1,nB
   DO iA=1,nA
    DO iD=1,nD
     tmpBD = D0
     DO iC=1,nC
        DO idmat=1,10
           tmpBD(idmat) = tmpBD(idmat) + DAC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
        ENDDO
     ENDDO
     DO idmat=1,10
        Kcont(idmat) = Kcont(idmat) + D2*tmpBD(idmat) * DBD(s2+iB,s4+iD,idmat)
     ENDDO
    ENDDO
   ENDDO
  ENDDO
 ENDIF
ELSE
 IF(.NOT.sym)THEN
  DO iB=1,nB
   DO iA=1,nA
    DO iD=1,nD
     tmpBD = D0
     DO iC=1,nC
        DO idmat=1,ndmat
           tmpBD(idmat) = tmpBD(idmat) + DAC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
        ENDDO
     ENDDO
     DO idmat=1,ndmat
        Kcont(idmat) = Kcont(idmat) + tmpBD(idmat) * DBD(s2+iB,s4+iD,idmat)
     ENDDO    
    ENDDO
   ENDDO
  ENDDO
 ELSE
  DO iB=1,nB
   DO iA=1,nA
    DO iD=1,nD
     tmpBD = D0
     DO iC=1,nC
        DO idmat=1,ndmat
           tmpBD(idmat) = tmpBD(idmat) + DAC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,iPass)
        ENDDO
     ENDDO
     DO idmat=1,ndmat
        Kcont(idmat) = Kcont(idmat) + D2 * tmpBD(idmat) * DBD(s2+iB,s4+iD,idmat)
     ENDDO    
    ENDDO
   ENDDO
  ENDDO
 ENDIF
ENDIF
end subroutine fullEcont5

END SUBROUTINE distributeKcont

!> \brief Distribute exchange-type integrals
!> \author \latexonly S. Reine  \endlatexonly
!> \date 2011-03-08
!> \param Kmat contains the result lstensor
!> \param CDAB matrix containing calculated integrals
!> \param Dmat contains the rhs-density matrix/matrices
!> \param ndmat number of rhs-density matries
!> \param PQ information about the calculated integrals (to be distributed)
!> \param Input contain info about the requested integral 
!> \param Lsoutput contain info about the requested output, and actual output
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE distributemagKcont(CDAB,Dlhs,Drhs,ndmat,PQ,Input,&
    &                      Lsoutput,Kcont,LUPRI,IPRINT)
implicit none 
Type(integrand),intent(in)         :: PQ
Type(IntegralInput),intent(in)     :: Input
Type(IntegralOutput),intent(inout) :: Lsoutput
Integer,intent(in)                 :: ndmat,LUPRI,IPRINT
REAL(REALK),pointer                :: CDAB(:)
TYPE(lstensor),intent(in)          :: Dlhs,Drhs
Real(Realk),intent(inout)          :: Kcont(3,ndmat)
!
type(overlap),pointer :: P,Q
logical :: SamePQ,SameLHSaos,SameRHSaos,SameODs
Logical :: AntipermuteAB,permuteCD
integer :: nAngmomP,nPasses,maxBat,maxAng
integer :: iPass,iPassP,iPassQ,idmat,DAr,CBr
integer :: nA,nB,nC,nD,batchA,batchB,batchC,batchD,atomA,atomB,atomC,atomD
integer :: AC,BD,DA,CB,AD,BC,n1,n2,n3,n4,s1,s2,s3,s4
integer :: CA,DB
Real(realk),pointer :: lAC(:),lBC(:),lAD(:),lBD(:),rBD(:),rAD(:),rBC(:),rAC(:)
Real(realk) :: tmpKcont(3,ndmat)
P => PQ%P%p
Q => PQ%Q%p
nAngmomP   = P%nAngmom
SamePQ     = PQ%samePQ
SameLHSaos = INPUT%SameLHSaos
SameRHSaos = INPUT%SameRHSaos
SameODs    = INPUT%SameODs
nPasses = P%nPasses*Q%nPasses
iPass = 0
tmpKcont = 0E0_realk
DO iPassP=1,P%nPasses
  atomA  = P%orb1atom(iPassP)
  batchA = P%orb1batch(iPassP)
  nA     = P%orbital1%totOrbitals
  atomB  = P%orb2atom(iPassP)
  batchB = P%orb2batch(iPassP)
  nB     = P%orbital2%totOrbitals

  antipermuteAB  = SameLHSaos.AND.((batchA.NE.batchB).OR.(atomA.NE.atomB))

  DO iPassQ=1,Q%nPasses
    iPass = iPass + 1

    atomC  = Q%orb1atom(iPassQ)
    batchC = Q%orb1batch(iPassQ)
    nC     = Q%orbital1%totOrbitals
    atomD  = Q%orb2atom(iPassQ)
    batchD = Q%orb2batch(iPassQ)
    nD     = Q%orbital2%totOrbitals
    permuteCD = SameRHSaos.AND.((batchC.NE.batchD).OR.(atomC.NE.atomD))

    !LHS
    AC = Dlhs%index(atomA,atomC,1,1)
    lAC => Dlhs%lsao(AC)%elms
    n1 = Dlhs%lsao(AC)%nLocal(1) 
    n3 = Dlhs%lsao(AC)%nLocal(2) 
    maxBat = Dlhs%lsao(AC)%maxBat
    maxAng = Dlhs%lsao(AC)%maxAng
    s1 = Dlhs%lsao(AC)%startLocalOrb(1+(batchA-1)*maxAng)-1
    s3 = Dlhs%lsao(AC)%startLocalOrb(1+(batchC-1)*maxAng+maxAng*maxBat)-1
    !RHS 
    BD = Drhs%index(atomB,atomD,1,1)
    rBD => Drhs%lsao(BD)%elms       
    n2 = Drhs%lsao(BD)%nLocal(1) 
    n4 = Drhs%lsao(BD)%nLocal(2) 
    maxBat = Drhs%lsao(BD)%maxBat
    maxAng = Drhs%lsao(BD)%maxAng
    s2 = Drhs%lsao(BD)%startLocalOrb(1+(batchB-1)*maxAng)-1
    s4 = Drhs%lsao(BD)%startLocalOrb(1+(batchD-1)*maxAng+maxAng*maxBat)-1

    IF (AntipermuteAB.AND.permuteCD)THEN
       BC = Dlhs%index(atomB,atomC,1,1)
       lBC => Dlhs%lsao(BC)%elms
       AD = Drhs%index(atomA,atomD,1,1)
       rAD => Drhs%lsao(AD)%elms       

       AD = Dlhs%index(atomA,atomD,1,1)
       lAD => Dlhs%lsao(AD)%elms       
       BC = Drhs%index(atomB,atomC,1,1)
       rBC => Drhs%lsao(BC)%elms

       BD = Dlhs%index(atomB,atomD,1,1)
       lBD => Dlhs%lsao(BD)%elms       
       AC = Drhs%index(atomA,atomC,1,1)
       rAC => Drhs%lsao(AC)%elms
       !Kcont = DLHS(A,C)*Integral(A,B,C,D)^b*DRHS(B,D)
       !Kcont = DLHS(B,C)*Integral(A,B,D,C)^b*DRHS(A,D)
       !Kcont = DLHS(A,D)*Integral(A,B,D,C)^b*DRHS(B,C)
       !Kcont = DLHS(B,D)*Integral(A,B,D,C)^b*DRHS(A,C)
       call magfullABCD(tmpKcont,lAC,rBD,lBC,rAD,lAD,rBC,lBD,rAC,&
            & n1,n2,n3,n4,s1,s2,s3,s4,CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,lupri)
    ELSEIF (AntipermuteAB) THEN
       !
       BC = Dlhs%index(atomB,atomC,1,1)
       lBC => Dlhs%lsao(BC)%elms
       AD = Drhs%index(atomA,atomD,1,1)
       rAD => Drhs%lsao(AD)%elms       
       !Kcont = DLHS(A,C)*Integral(A,B,C,D)^b*DRHS(B,D)
       !Kcont = DLHS(B,C)*Integral(A,B,D,C)^b*DRHS(A,D)
       call magfullAB(tmpKcont,lAC,rBD,lBC,rAD,&
            & n1,n2,n3,n4,s1,s2,s3,s4,CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,lupri)
    ELSEIF (permuteCD) THEN
       AD = Dlhs%index(atomA,atomD,1,1)
       lAD => Dlhs%lsao(AD)%elms       
       BC = Drhs%index(atomB,atomC,1,1)
       rBC => Drhs%lsao(BC)%elms
       !Kcont = DLHS(A,C)*Integral(A,B,C,D)^b*DRHS(B,D)
       !Kcont = DLHS(A,D)*Integral(A,B,D,C)^b*DRHS(B,C)
       call magfullCD(tmpKcont,lAC,rBD,lAD,rBC,&
            & n1,n2,n3,n4,s1,s2,s3,s4,CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,lupri)
    ELSE
       !Kcont = DLHS(A,C)*Integral(A,B,C,D)^b*DRHS(B,D)
       call magfull(tmpKcont,lAC,rBD,&
            & n1,n2,n3,n4,s1,s2,s3,s4,CDAB,nA,nB,nC,nD,iPass,nPasses,ndmat,lupri)
    ENDIF
  ENDDO !iPassQ
ENDDO !iPassP
!reduction afterwards - each thread have its own version of Kcont
DO idmat=1,ndmat
   Kcont(1,idmat) = Kcont(1,idmat) - tmpKcont(1,idmat)
   Kcont(2,idmat) = Kcont(2,idmat) - tmpKcont(2,idmat)
   Kcont(3,idmat) = Kcont(3,idmat) - tmpKcont(3,idmat)
ENDDO

CONTAINS
subroutine magfullABCD(Kcont,DLHS_AC,DRHS_BD,DLHS_BC,DRHS_AD,&
     & DLHS_AD,DRHS_BC,DLHS_BD,DRHS_AC,&
     & n1,n2,n3,n4,s1,s2,s3,s4,CDAB,nA,nB,nC,nD,iPass,nPasses,nLHS,lupri)
implicit none
Integer,intent(IN)     :: nA,nB,nC,nD,iPass,nPasses,lupri,nLHS
Integer,intent(IN)     :: n1,n2,n3,n4,s1,s2,s3,s4
Real(realk),intent(IN) :: CDAB(nC,nD,nA,nB,3,nPasses)
Real(realk),intent(IN) :: DLHS_AC(n1,n3,nLHS),DRHS_BD(n2,n4)
Real(realk),intent(IN) :: DLHS_BC(n2,n3,nLHS),DRHS_AD(n1,n4)
Real(realk),intent(IN) :: DLHS_AD(n1,n4,nLHS),DRHS_BC(n2,n3)
Real(realk),intent(IN) :: DLHS_BD(n2,n4,nLHS),DRHS_AC(n1,n3)
Real(realk),intent(INOUT) :: Kcont(3,nLHS)
!
Integer     :: iA,iB,iC,iD,idmat
Real(realk) :: tmpBD(3),tmpAD(3),tmpBC(3),tmpAC(3)
real(realk),parameter :: D0=0.0E0_realk
DO iB=1,nB
 DO iA=1,nA
  DO iD=1,nD
   DO idmat=1,nLHS
    tmpBD(1) = D0
    tmpBD(2) = D0
    tmpBD(3) = D0
    tmpAD(1) = D0
    tmpAD(2) = D0
    tmpAD(3) = D0
    tmpBC(1) = D0
    tmpBC(2) = D0
    tmpBC(3) = D0
    tmpAC(1) = D0
    tmpAC(2) = D0
    tmpAC(3) = D0
    !Kcont(A,B,C,D) = DLHS(A,C)*Integral((AB)^b|CD)*DRHS(B,D)
    !Kcont(B,A,C,D) = DLHS(B,C)*Integral((BA)^b|CD)*DRHS(A,D) = - DLHS(B,C)*Integral((AB)^b|CD)*DRHS(A,D)
    !Kcont(A,B,D,C) = DLHS(A,D)*Integral((AB)^b|DC)*DRHS(B,C) =   DLHS(A,D)*Integral((AB)^b|CD)*DRHS(B,C)
    !Kcont(B,A,D,C) = DLHS(B,D)*Integral((BA)^b|DC)*DRHS(A,C) = - DLHS(B,D)*Integral((AB)^b|CD)*DRHS(A,C)
    DO iC=1,nC
     tmpBD(1) = tmpBD(1) + DLHS_AC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,1,iPass)
     tmpBD(2) = tmpBD(2) + DLHS_AC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,2,iPass)
     tmpBD(3) = tmpBD(3) + DLHS_AC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,3,iPass)
     !permute AB
     tmpAD(1) = tmpAD(1) + DLHS_BC(s2+iB,s3+iC,idmat) * CDAB(iC,iD,iA,iB,1,iPass)
     tmpAD(2) = tmpAD(2) + DLHS_BC(s2+iB,s3+iC,idmat) * CDAB(iC,iD,iA,iB,2,iPass)
     tmpAD(3) = tmpAD(3) + DLHS_BC(s2+iB,s3+iC,idmat) * CDAB(iC,iD,iA,iB,3,iPass)
     !permute CD
     tmpBC(1) = tmpBC(1) + DRHS_BC(s2+iB,s3+iC) * CDAB(iC,iD,iA,iB,1,iPass)
     tmpBC(2) = tmpBC(2) + DRHS_BC(s2+iB,s3+iC) * CDAB(iC,iD,iA,iB,2,iPass)
     tmpBC(3) = tmpBC(3) + DRHS_BC(s2+iB,s3+iC) * CDAB(iC,iD,iA,iB,3,iPass)
     !permute AB and CD
     tmpAC(1) = tmpAC(1) + DRHS_AC(s1+iA,s3+iC) * CDAB(iC,iD,iA,iB,1,iPass)
     tmpAC(2) = tmpAC(2) + DRHS_AC(s1+iA,s3+iC) * CDAB(iC,iD,iA,iB,2,iPass)
     tmpAC(3) = tmpAC(3) + DRHS_AC(s1+iA,s3+iC) * CDAB(iC,iD,iA,iB,3,iPass)
    ENDDO
    !Kcont = DLHS(A,C)*Integral(A,B,C,D)^b*DRHS(B,D)
    Kcont(1,idmat) = Kcont(1,idmat) + tmpBD(1) * DRHS_BD(s2+iB,s4+iD)
    Kcont(2,idmat) = Kcont(2,idmat) + tmpBD(2) * DRHS_BD(s2+iB,s4+iD)
    Kcont(3,idmat) = Kcont(3,idmat) + tmpBD(3) * DRHS_BD(s2+iB,s4+iD)
    !Kcont = DLHS(B,C)*Integral(B,A,C,D)^b*DRHS(A,D) (permute AB)
    Kcont(1,idmat) = Kcont(1,idmat) - tmpAD(1) * DRHS_AD(s1+iA,s4+iD) 
    Kcont(2,idmat) = Kcont(2,idmat) - tmpAD(2) * DRHS_AD(s1+iA,s4+iD)
    Kcont(3,idmat) = Kcont(3,idmat) - tmpAD(3) * DRHS_AD(s1+iA,s4+iD)
    !Kcont = DLHS(A,D)*Integral(A,B,D,C)^b*DRHS(B,C) (permute CD)
    Kcont(1,idmat) = Kcont(1,idmat) + tmpBC(1) * DLHS_AD(s1+iA,s4+iD,idmat) 
    Kcont(2,idmat) = Kcont(2,idmat) + tmpBC(2) * DLHS_AD(s1+iA,s4+iD,idmat)
    Kcont(3,idmat) = Kcont(3,idmat) + tmpBC(3) * DLHS_AD(s1+iA,s4+iD,idmat)
    !Kcont = DLHS(B,D)*Integral(B,A,C,D)^b*DRHS(A,C) (permute AB and CD)
    Kcont(1,idmat) = Kcont(1,idmat) - tmpAC(1) * DLHS_BD(s2+iB,s4+iD,idmat) 
    Kcont(2,idmat) = Kcont(2,idmat) - tmpAC(2) * DLHS_BD(s2+iB,s4+iD,idmat)
    Kcont(3,idmat) = Kcont(3,idmat) - tmpAC(3) * DLHS_BD(s2+iB,s4+iD,idmat)
   ENDDO    
  ENDDO
 ENDDO
ENDDO
end subroutine magfullABCD

subroutine magfull(Kcont,DLHS_AC,DRHS_BD,&
     & n1,n2,n3,n4,s1,s2,s3,s4,&
     & CDAB,nA,nB,nC,nD,iPass,nPasses,nLHS,lupri)
implicit none
Integer,intent(IN)     :: nA,nB,nC,nD,iPass,nPasses,lupri,nLHS
Integer,intent(IN)     :: n1,n2,n3,n4,s1,s2,s3,s4
Real(realk),intent(IN) :: CDAB(nC,nD,nA,nB,3,nPasses)
Real(realk),intent(IN) :: DLHS_AC(n1,n3,nLHS),DRHS_BD(n2,n4)
Real(realk),intent(INOUT) :: Kcont(3,nLHS)
!
Integer     :: iA,iB,iC,iD,idmat
Real(realk) :: tmpBD(3),tmpDB(3)
real(realk),parameter :: D0=0.0E0_realk
DO iB=1,nB
 DO iA=1,nA
  DO iD=1,nD
   DO idmat=1,nLHS
    tmpBD(1) = D0
    tmpBD(2) = D0
    tmpBD(3) = D0
    !Kcont = DLHS(A,C)*((AB)^b|CD)*DRHS(B,D)
    DO iC=1,nC
     tmpBD(1) = tmpBD(1) + DLHS_AC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,1,iPass)
     tmpBD(2) = tmpBD(2) + DLHS_AC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,2,iPass)
     tmpBD(3) = tmpBD(3) + DLHS_AC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,3,iPass)
    ENDDO
    Kcont(1,idmat) = Kcont(1,idmat) + tmpBD(1) * DRHS_BD(s2+iB,s4+iD)
    Kcont(2,idmat) = Kcont(2,idmat) + tmpBD(2) * DRHS_BD(s2+iB,s4+iD)
    Kcont(3,idmat) = Kcont(3,idmat) + tmpBD(3) * DRHS_BD(s2+iB,s4+iD)
   ENDDO    
  ENDDO
 ENDDO
ENDDO
end subroutine magfull

!include AB permutational symmetry
subroutine magfullAB(Kcont,DLHS_AC,DRHS_BD,DLHS_BC,DRHS_AD,&
     & n1,n2,n3,n4,s1,s2,s3,s4,&
     & CDAB,nA,nB,nC,nD,iPass,nPasses,nLHS,lupri)
implicit none
Integer,intent(IN)     :: nA,nB,nC,nD,iPass,nPasses,lupri,nLHS
Integer,intent(IN)     :: n1,n2,n3,n4,s1,s2,s3,s4
Real(realk),intent(IN) :: CDAB(nC,nD,nA,nB,3,nPasses)
Real(realk),intent(IN) :: DRHS_BD(n2,n4),DLHS_AC(n1,n3,nLHS)
Real(realk),intent(IN) :: DRHS_AD(n1,n4),DLHS_BC(n2,n3,nLHS)
Real(realk),intent(INOUT) :: Kcont(3,nLHS)
!
Integer     :: iA,iB,iC,iD,idmat
Real(realk) :: tmpBD(3),tmpBC(3)
real(realk),parameter :: D0=0.0E0_realk
DO iB=1,nB
 DO iA=1,nA
  DO iD=1,nD
   DO idmat=1,nLHS
    tmpBD(1) = D0
    tmpBD(2) = D0
    tmpBD(3) = D0
    tmpBC(1) = D0
    tmpBC(2) = D0
    tmpBC(3) = D0
    !Kcont(A,B,C,D) =  DLHS(A,C)*((AB)^b|CD)*DRHS(B,D) 
    !Kcont(B,A,C,D) =  DLHS(B,C)*((BA)^b|CD)*DRHS(A,D) = - DLHS(B,C)*((AB)^b|CD)*DRHS(A,D)
    DO iC=1,nC
     tmpBD(1) = tmpBD(1) + DLHS_AC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,1,iPass)
     tmpBD(2) = tmpBD(2) + DLHS_AC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,2,iPass)
     tmpBD(3) = tmpBD(3) + DLHS_AC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,3,iPass)
     tmpBC(1) = tmpBC(1) + DLHS_BC(s2+iB,s3+iC,idmat) * CDAB(iC,iD,iA,iB,1,iPass)
     tmpBC(2) = tmpBC(2) + DLHS_BC(s2+iB,s3+iC,idmat) * CDAB(iC,iD,iA,iB,2,iPass)
     tmpBC(3) = tmpBC(3) + DLHS_BC(s2+iB,s3+iC,idmat) * CDAB(iC,iD,iA,iB,3,iPass)
    ENDDO
    !Kcont = DLHS(A,C)*Integral(A,B,C,D)^b*DRHS(B,D)
    Kcont(1,idmat) = Kcont(1,idmat) + tmpBD(1) * DRHS_BD(s2+iB,s4+iD)
    Kcont(2,idmat) = Kcont(2,idmat) + tmpBD(2) * DRHS_BD(s2+iB,s4+iD)
    Kcont(3,idmat) = Kcont(3,idmat) + tmpBD(3) * DRHS_BD(s2+iB,s4+iD)
    !Kcont =  DLHS(B,C)*Integral(B,A,C,D)^b*DRHS(A,D)
    !Kcont = -DLHS(B,C)*Integral(A,B,C,D)^b*DRHS(A,D)
    Kcont(1,idmat) = Kcont(1,idmat) - tmpBC(1) * DRHS_AD(s1+iA,s4+iD) 
    Kcont(2,idmat) = Kcont(2,idmat) - tmpBC(2) * DRHS_AD(s1+iA,s4+iD)
    Kcont(3,idmat) = Kcont(3,idmat) - tmpBC(3) * DRHS_AD(s1+iA,s4+iD)
   ENDDO    
  ENDDO
 ENDDO
ENDDO
end subroutine magfullAB

!include CD permutational symmetry
subroutine magfullCD(Kcont,DLHS_AC,DRHS_BD,DLHS_AD,DRHS_BC,&
     & n1,n2,n3,n4,s1,s2,s3,s4,&
     & CDAB,nA,nB,nC,nD,iPass,nPasses,nLHS,lupri)
implicit none
Integer,intent(IN)     :: nA,nB,nC,nD,iPass,nPasses,lupri,nLHS
Integer,intent(IN)     :: n1,n2,n3,n4,s1,s2,s3,s4
Real(realk),intent(IN) :: CDAB(nC,nD,nA,nB,3,nPasses)
Real(realk),intent(IN) :: DRHS_BD(n2,n4),DLHS_AC(n1,n3,nLHS)
Real(realk),intent(IN) :: DRHS_BC(n2,n3),DLHS_AD(n1,n4,nLHS)
Real(realk),intent(INOUT) :: Kcont(3,nLHS)
!
Integer     :: iA,iB,iC,iD,idmat
Real(realk) :: tmpBD(3),tmpBC(3)
real(realk),parameter :: D0=0.0E0_realk
DO iB=1,nB
 DO iA=1,nA
  DO iD=1,nD
   DO idmat=1,nLHS
    tmpBD(1) = D0
    tmpBD(2) = D0
    tmpBD(3) = D0
    tmpBC(1) = D0
    tmpBC(2) = D0
    tmpBC(3) = D0
    !Kcont(A,B,C,D) = DLHS(A,C)*Integral((AB)^b|CD)*DRHS(B,D)
    !Kcont(A,B,D,C) = DLHS(A,D)*Integral((AB)^b|DC)*DRHS(B,C) = DLHS(A,D)*Integral((AB)^b|CD)*DRHS(B,C)
    DO iC=1,nC
     tmpBD(1) = tmpBD(1) + DLHS_AC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,1,iPass)
     tmpBD(2) = tmpBD(2) + DLHS_AC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,2,iPass)
     tmpBD(3) = tmpBD(3) + DLHS_AC(s1+iA,s3+iC,idmat) * CDAB(iC,iD,iA,iB,3,iPass)
     tmpBC(1) = tmpBC(1) + DRHS_BC(s2+iB,s3+iC) * CDAB(iC,iD,iA,iB,1,iPass)
     tmpBC(2) = tmpBC(2) + DRHS_BC(s2+iB,s3+iC) * CDAB(iC,iD,iA,iB,2,iPass)
     tmpBC(3) = tmpBC(3) + DRHS_BC(s2+iB,s3+iC) * CDAB(iC,iD,iA,iB,3,iPass)
    ENDDO
    !Kcont(A,B,C,D) = DLHS(A,C)*Integral((AB)^b|CD)*DRHS(B,D)
    Kcont(1,idmat) = Kcont(1,idmat) + tmpBD(1) * DRHS_BD(s2+iB,s4+iD)
    Kcont(2,idmat) = Kcont(2,idmat) + tmpBD(2) * DRHS_BD(s2+iB,s4+iD)
    Kcont(3,idmat) = Kcont(3,idmat) + tmpBD(3) * DRHS_BD(s2+iB,s4+iD)
    !Kcont(A,B,D,C) = DLHS(A,D)*Integral((AB)^b|CD)*DRHS(B,C)
    Kcont(1,idmat) = Kcont(1,idmat) + tmpBC(1) * DLHS_AD(s1+iA,s4+iD,idmat) 
    Kcont(2,idmat) = Kcont(2,idmat) + tmpBC(2) * DLHS_AD(s1+iA,s4+iD,idmat)
    Kcont(3,idmat) = Kcont(3,idmat) + tmpBC(3) * DLHS_AD(s1+iA,s4+iD,idmat)
   ENDDO    
  ENDDO
 ENDDO
ENDDO
end subroutine magfullCD

END SUBROUTINE distributemagKcont

END MODULE thermite_distributeK2


