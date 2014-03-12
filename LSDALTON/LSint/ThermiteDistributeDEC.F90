!> @file 
!> Contains integral distribution routines that places calculated integrals in the proper output
!> Thermite integral distribution module
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
MODULE thermite_distributeDEC
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
!DEC wants the integrals in (nbast,nbast,dim3,dim4) but it is faster 
!to calculate them as (dim3,dim4,nbast,nbast) 
!so we calculate them as (dim3,dim4,nbast,nbast) but store the 
!result as (nbast,nbast,dim3,dim4) 
! so DECRES have the dimension (dim3,dim4,dim1,dim2,1) = (nbast,nbast,dim3,dim4)
SUBROUTINE Explicit4CenterDEC(DECRES,dim1,dim2,dim3,dim4,&
     & PQ,QPmat2,dimQ,dimP,Input,Lsoutput,LUPRI,IPRINT)
implicit none 
Type(integrand),intent(in)      :: PQ
Type(IntegralInput),intent(in)  :: Input
Type(IntegralOutput),intent(inout) :: Lsoutput
Integer,intent(in)              :: LUPRI,IPRINT
Integer,intent(in)              :: dimQ,dimP,dim1,dim2,dim3,dim4
REAL(REALK),intent(in)          :: QPMAT2(dimQ,dimP)
REAL(REALK),intent(inout)       :: DECRES(dim3,dim4,dim1,dim2,1)
!
type(overlap),pointer :: P,Q
logical :: SamePQ,SameLHSaos,SameRHSaos,SameODs,nucleiP,nucleiQ,add
Logical :: permuteAB,permuteCD,permuteOD,noContraction,RHScontraction
integer :: nAngmomP,nmat,nPasses
integer :: iPass,iPassP,iPassQ,iDeriv,iDerivP,iDerivQ,nAngmomQ
integer :: nA,nB,nC,nD,batchA,batchB,batchC,batchD,atomA,atomB,atomC,atomD
integer :: indADBC,indADCB,indDACB,indDABC,indBCAD,indCBAD,indBCDA,indCBDA
Real(realk),pointer :: kADBC(:),kADCB(:),kBCAD(:),kCBAD(:),kDACB(:),kDABC(:),kBCDA(:),kCBDA(:)

integer :: iA,iB,iC,iD,nContA,nContB,nContC,nContD
integer :: nAngA,nAngB,nAngC,nAngD,iAngmomP,iAngmomQ
integer :: startA,startB,startC,startD

nAngmomP       = PQ%P%p%nAngmom
nAngmomQ       = PQ%Q%p%nAngmom
SamePQ         = PQ%samePQ
SameLHSaos     = INPUT%SameLHSaos
SameRHSaos     = INPUT%SameRHSaos
SameODs        = INPUT%SameODs
nmat           = 1

P => PQ%P%p
Q => PQ%Q%p

nPasses = P%nPasses*Q%nPasses
DO iAngmomP=1,nAngmomP
 iA = P%indexAng1(iAngmomP)
 iB = P%indexAng2(iAngmomP)
 nContA = P%orbital1%nContracted(iA)
 nContB = P%orbital2%nContracted(iB)
 nAngA = P%orbital1%nOrbComp(iA)
 nAngB = P%orbital2%nOrbComp(iB)
 DO iPassP=1,P%nPasses
  startA  = P%orbital1%startOrbital(iA+(iPassP-1)*P%orbital1%nAngmom) - 1
  startB  = P%orbital2%startOrbital(iB+(iPassP-1)*P%orbital2%nAngmom) - 1    

  atomA  = P%orb1atom(iPassP)
  batchA = P%orb1batch(iPassP)
  nA     = P%orbital1%totOrbitals
  atomB  = P%orb2atom(iPassP)
  batchB = P%orb2batch(iPassP)
  nB     = P%orbital2%totOrbitals
  permuteAB  = SameLHSaos.AND.((batchA.NE.batchB).OR.(atomA.NE.atomB))
  DO iAngmomQ=1,nAngmomQ
   iC = Q%indexAng1(iAngmomQ)
   iD = Q%indexAng2(iAngmomQ)
   nContC = Q%orbital1%nContracted(iC)
   nContD = Q%orbital2%nContracted(iD)
   nAngC = Q%orbital1%nOrbComp(iC)
   nAngD = Q%orbital2%nOrbComp(iD)
   DO iPassQ=1,Q%nPasses
    iPass = iPassQ + (iPassP-1)*Q%nPasses
    startC  = Q%orbital1%startOrbital(iC+(iPassQ-1)*Q%orbital1%nAngmom) - 1
    startD  = Q%orbital2%startOrbital(iD+(iPassQ-1)*Q%orbital2%nAngmom) - 1    
    atomC  = Q%orb1atom(iPassQ)
    batchC = Q%orb1batch(iPassQ)
    nC     = Q%orbital1%totOrbitals
    atomD  = Q%orb2atom(iPassQ)
    batchD = Q%orb2batch(iPassQ)
    nD     = Q%orbital2%totOrbitals
    permuteCD = SameRHSaos.AND.((batchC.NE.batchD).OR.(atomC.NE.atomD))
    permuteOD = SameODs.AND.( ((batchA.NE.batchC).OR.(atomA.NE.atomC)).OR.&
   &                          ((batchB.NE.batchD).OR.(atomB.NE.atomD)) )
    IF (permuteOD.AND.permuteCD.AND.permuteAB) THEN
       call DECdistributeSUbFULL(DECRES,dim1,dim2,dim3,dim4,QPmat2,nA,nB,nC,nD,iPass,nPasses,&
            & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
    ELSE IF (permuteOD.AND.permuteCD) THEN
       call DECdistributeSUb7(DECRES,dim1,dim2,dim3,dim4,QPmat2,nA,nB,nC,nD,iPass,nPasses,&
            & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
    ELSE IF (permuteOD.AND.permuteAB) THEN
       call DECdistributeSUb6(DECRES,dim1,dim2,dim3,dim4,QPmat2,nA,nB,nC,nD,iPass,nPasses,&
            & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
    ELSE IF (permuteAB.AND.permuteCD) THEN
       call DECdistributeSUb5(DECRES,dim1,dim2,dim3,dim4,QPmat2,nA,nB,nC,nD,iPass,nPasses,&
            & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
    ELSE IF (permuteAB) THEN
       call DECdistributeSUb4(DECRES,dim1,dim2,dim3,dim4,QPmat2,nA,nB,nC,nD,iPass,nPasses,&
            & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
    ELSE IF (permuteCD) THEN
       call DECdistributeSUb3(DECRES,dim1,dim2,dim3,dim4,QPmat2,nA,nB,nC,nD,iPass,nPasses,&
            & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
    ELSE IF (permuteOD) THEN
       call DECdistributeSUb2(DECRES,dim1,dim2,dim3,dim4,QPmat2,nA,nB,nC,nD,iPass,nPasses,&
            & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
    ELSE
       call DECdistributeSUb1(DECRES,dim1,dim2,dim3,dim4,QPmat2,nA,nB,nC,nD,iPass,nPasses,&
            & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
    ENDIF
   ENDDO !iPassQ
  ENDDO !iAngmomQ
 ENDDO !iPassP
ENDDO !iAngmomP

CONTAINS
  subroutine DECdistributeSUbFULL(DECRES,dim1,dim2,dim3,dim4,CDAB,nA,nB,nC,nD,iPass,nPasses,&
       & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
    implicit none
    integer,intent(in) :: dim1,dim2,dim3,dim4,nA,nB,nC,nD,iPass,nPasses,lupri
    integer,intent(in) :: startA,startB,startC,startD,nContA,nContB,nContC,nContD
    integer,intent(in) :: nAngA,nAngB,nAngC,nAngD
    real(realk) :: DECRES(dim3,dim4,dim1,dim2,1)
    real(realk) :: CDAB(nC,nD,nA,nB,nPasses)
    !
    integer :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
    integer :: iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,sC1,sC2
    real(realk) :: int
    
    DO iAngB=1,nAngB
     iB2 = (iAngB-1)*nContB
     iB1 = startB+iB2
      DO iContB=1,nContB
       iB2 = iB2+1
       iB1 = iB1+1
       DO iAngA=1,nAngA
        iA2 = (iAngA-1)*nContA
        iA1 = startA+iA2
        DO iContA=1,nContA
         iA2 = iA2+1
         iA1 = iA1+1
         DO iAngD=1,nAngD
          iD2 = (iAngD-1)*nContD
          iD1 = startD+iD2
          DO iContD=1,nContD
           iD1=iD1+1
           iD2=iD2+1
           DO iAngC=1,nAngC
            sC1 = startC + (iAngC-1)*nContC
            sC2 = (iAngC-1)*nContC
            DO iContC=1,nContC
             int = CDAB(sC2+iContC,iD2,iA2,iB2,iPass)
             DECRES(sC1+iContC,iD1,iA1,iB1,1)=int          !ABCD
             DECRES(sC1+iContC,iD1,iB1,iA1,1)=int          !BACD
             DECRES(iD1,sC1+iContC,iA1,iB1,1)=int          !ABDC
             DECRES(iD1,sC1+iContC,iB1,iA1,1)=int          !BADC
             DECRES(iA1,iB1,sC1+iContC,iD1,1)=int          !CDAB
             DECRES(iA1,iB1,iD1,sC1+iContC,1)=int          !DCAB
             DECRES(iB1,iA1,sC1+iContC,iD1,1)=int          !CDBA
             DECRES(iB1,iA1,iD1,sC1+iContC,1)=int          !DCBA
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
   END subroutine DECDISTRIBUTESUBFULL

  subroutine DECdistributeSUb1(DECRES,dim1,dim2,dim3,dim4,CDAB,nA,nB,nC,nD,iPass,nPasses,&
       & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
    implicit none
    integer,intent(in) :: dim1,dim2,dim3,dim4,nA,nB,nC,nD,iPass,nPasses,lupri
    integer,intent(in) :: startA,startB,startC,startD,nContA,nContB,nContC,nContD
    integer,intent(in) :: nAngA,nAngB,nAngC,nAngD
    real(realk) :: DECRES(dim3,dim4,dim1,dim2,1)
    real(realk) :: CDAB(nC,nD,nA,nB,nPasses)
    !
    integer :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
    integer :: iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,sC1,sC2
    real(realk) :: int
    
    DO iAngB=1,nAngB
     iB2 = (iAngB-1)*nContB
     iB1 = startB+iB2
      DO iContB=1,nContB
       iB2 = iB2+1
       iB1 = iB1+1
       DO iAngA=1,nAngA
        iA2 = (iAngA-1)*nContA
        iA1 = startA+iA2
        DO iContA=1,nContA
         iA2 = iA2+1
         iA1 = iA1+1
         DO iAngD=1,nAngD
          iD2 = (iAngD-1)*nContD
          iD1 = startD+iD2
          DO iContD=1,nContD
           iD1=iD1+1
           iD2=iD2+1
           DO iAngC=1,nAngC
            sC1 = startC + (iAngC-1)*nContC
            sC2 = (iAngC-1)*nContC
            DO iContC=1,nContC
             DECRES(sC1+iContC,iD1,iA1,iB1,1)=CDAB(sC2+iContC,iD2,iA2,iB2,iPass)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
   END subroutine DECDISTRIBUTESUB1

  subroutine DECdistributeSUb7(DECRES,dim1,dim2,dim3,dim4,CDAB,nA,nB,nC,nD,iPass,nPasses,&
       & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
    implicit none
    integer,intent(in) :: dim1,dim2,dim3,dim4,nA,nB,nC,nD,iPass,nPasses,lupri
    integer,intent(in) :: startA,startB,startC,startD,nContA,nContB,nContC,nContD
    integer,intent(in) :: nAngA,nAngB,nAngC,nAngD
    real(realk) :: DECRES(dim3,dim4,dim1,dim2,1)
    real(realk) :: CDAB(nC,nD,nA,nB,nPasses)
    !
    integer :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
    integer :: iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,sC1,sC2
    real(realk) :: int
    
    DO iAngB=1,nAngB
     iB2 = (iAngB-1)*nContB
     iB1 = startB+iB2
      DO iContB=1,nContB
       iB2 = iB2+1
       iB1 = iB1+1
       DO iAngA=1,nAngA
        iA2 = (iAngA-1)*nContA
        iA1 = startA+iA2
        DO iContA=1,nContA
         iA2 = iA2+1
         iA1 = iA1+1
         DO iAngD=1,nAngD
          iD2 = (iAngD-1)*nContD
          iD1 = startD+iD2
          DO iContD=1,nContD
           iD1=iD1+1
           iD2=iD2+1
           DO iAngC=1,nAngC
            sC1 = startC + (iAngC-1)*nContC
            sC2 = (iAngC-1)*nContC
            DO iContC=1,nContC
             int = CDAB(sC2+iContC,iD2,iA2,iB2,iPass)
             DECRES(sC1+iContC,iD1,iA1,iB1,1)=int          !ABCD
             DECRES(iD1,sC1+iContC,iA1,iB1,1)=int          !ABDC
             DECRES(iA1,iB1,sC1+iContC,iD1,1)=int          !CDAB
             DECRES(iA1,iB1,iD1,sC1+iContC,1)=int          !DCAB
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
   END subroutine DECDISTRIBUTESUB7

  subroutine DECdistributeSUb6(DECRES,dim1,dim2,dim3,dim4,CDAB,nA,nB,nC,nD,iPass,nPasses,&
       & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
    implicit none
    integer,intent(in) :: dim1,dim2,dim3,dim4,nA,nB,nC,nD,iPass,nPasses,lupri
    integer,intent(in) :: startA,startB,startC,startD,nContA,nContB,nContC,nContD
    integer,intent(in) :: nAngA,nAngB,nAngC,nAngD
    real(realk) :: DECRES(dim3,dim4,dim1,dim2,1)
    real(realk) :: CDAB(nC,nD,nA,nB,nPasses)
    !
    integer :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
    integer :: iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,sC1,sC2
    real(realk) :: int
    
    DO iAngB=1,nAngB
     iB2 = (iAngB-1)*nContB
     iB1 = startB+iB2
      DO iContB=1,nContB
       iB2 = iB2+1
       iB1 = iB1+1
       DO iAngA=1,nAngA
        iA2 = (iAngA-1)*nContA
        iA1 = startA+iA2
        DO iContA=1,nContA
         iA2 = iA2+1
         iA1 = iA1+1
         DO iAngD=1,nAngD
          iD2 = (iAngD-1)*nContD
          iD1 = startD+iD2
          DO iContD=1,nContD
           iD1=iD1+1
           iD2=iD2+1
           DO iAngC=1,nAngC
            sC1 = startC + (iAngC-1)*nContC
            sC2 = (iAngC-1)*nContC
            DO iContC=1,nContC
             int = CDAB(sC2+iContC,iD2,iA2,iB2,iPass)
             DECRES(sC1+iContC,iD1,iA1,iB1,1)=int          !ABCD
             DECRES(sC1+iContC,iD1,iB1,iA1,1)=int          !BACD
             DECRES(iA1,iB1,sC1+iContC,iD1,1)=int          !CDAB
             DECRES(iB1,iA1,sC1+iContC,iD1,1)=int          !CDBA
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
   END subroutine DECDISTRIBUTESUB6

  subroutine DECdistributeSUb5(DECRES,dim1,dim2,dim3,dim4,CDAB,nA,nB,nC,nD,iPass,nPasses,&
       & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
    implicit none
    integer,intent(in) :: dim1,dim2,dim3,dim4,nA,nB,nC,nD,iPass,nPasses,lupri
    integer,intent(in) :: startA,startB,startC,startD,nContA,nContB,nContC,nContD
    integer,intent(in) :: nAngA,nAngB,nAngC,nAngD
    real(realk) :: DECRES(dim3,dim4,dim1,dim2,1)
    real(realk) :: CDAB(nC,nD,nA,nB,nPasses)
    !
    integer :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
    integer :: iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,sC1,sC2
    real(realk) :: int
    
    DO iAngB=1,nAngB
     iB2 = (iAngB-1)*nContB
     iB1 = startB+iB2
      DO iContB=1,nContB
       iB2 = iB2+1
       iB1 = iB1+1
       DO iAngA=1,nAngA
        iA2 = (iAngA-1)*nContA
        iA1 = startA+iA2
        DO iContA=1,nContA
         iA2 = iA2+1
         iA1 = iA1+1
         DO iAngD=1,nAngD
          iD2 = (iAngD-1)*nContD
          iD1 = startD+iD2
          DO iContD=1,nContD
           iD1=iD1+1
           iD2=iD2+1
           DO iAngC=1,nAngC
            sC1 = startC + (iAngC-1)*nContC
            sC2 = (iAngC-1)*nContC
            DO iContC=1,nContC
             int = CDAB(sC2+iContC,iD2,iA2,iB2,iPass)
             DECRES(sC1+iContC,iD1,iA1,iB1,1)=int          !ABCD
             DECRES(iD1,sC1+iContC,iA1,iB1,1)=int          !ABDC
             DECRES(sC1+iContC,iD1,iB1,iA1,1)=int          !BACD
             DECRES(iD1,sC1+iContC,iB1,iA1,1)=int          !BADC
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
   END subroutine DECDISTRIBUTESUB5

  subroutine DECdistributeSUb4(DECRES,dim1,dim2,dim3,dim4,CDAB,nA,nB,nC,nD,iPass,nPasses,&
       & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
    implicit none
    integer,intent(in) :: dim1,dim2,dim3,dim4,nA,nB,nC,nD,iPass,nPasses,lupri
    integer,intent(in) :: startA,startB,startC,startD,nContA,nContB,nContC,nContD
    integer,intent(in) :: nAngA,nAngB,nAngC,nAngD
    real(realk) :: DECRES(dim3,dim4,dim1,dim2,1)
    real(realk) :: CDAB(nC,nD,nA,nB,nPasses)
    !
    integer :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
    integer :: iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,sC1,sC2
    real(realk) :: int
    
    DO iAngB=1,nAngB
     iB2 = (iAngB-1)*nContB
     iB1 = startB+iB2
      DO iContB=1,nContB
       iB2 = iB2+1
       iB1 = iB1+1
       DO iAngA=1,nAngA
        iA2 = (iAngA-1)*nContA
        iA1 = startA+iA2
        DO iContA=1,nContA
         iA2 = iA2+1
         iA1 = iA1+1
         DO iAngD=1,nAngD
          iD2 = (iAngD-1)*nContD
          iD1 = startD+iD2
          DO iContD=1,nContD
           iD1=iD1+1
           iD2=iD2+1
           DO iAngC=1,nAngC
            sC1 = startC + (iAngC-1)*nContC
            sC2 = (iAngC-1)*nContC
            DO iContC=1,nContC
             int = CDAB(sC2+iContC,iD2,iA2,iB2,iPass)
             DECRES(sC1+iContC,iD1,iA1,iB1,1)=int          !ABCD
             DECRES(sC1+iContC,iD1,iB1,iA1,1)=int          !BACD
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
   END subroutine DECDISTRIBUTESUB4

  subroutine DECdistributeSUb3(DECRES,dim1,dim2,dim3,dim4,CDAB,nA,nB,nC,nD,iPass,nPasses,&
       & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
    implicit none
    integer,intent(in) :: dim1,dim2,dim3,dim4,nA,nB,nC,nD,iPass,nPasses,lupri
    integer,intent(in) :: startA,startB,startC,startD,nContA,nContB,nContC,nContD
    integer,intent(in) :: nAngA,nAngB,nAngC,nAngD
    real(realk) :: DECRES(dim3,dim4,dim1,dim2,1)
    real(realk) :: CDAB(nC,nD,nA,nB,nPasses)
    !
    integer :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
    integer :: iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,sC1,sC2
    real(realk) :: int
    
    DO iAngB=1,nAngB
     iB2 = (iAngB-1)*nContB
     iB1 = startB+iB2
      DO iContB=1,nContB
       iB2 = iB2+1
       iB1 = iB1+1
       DO iAngA=1,nAngA
        iA2 = (iAngA-1)*nContA
        iA1 = startA+iA2
        DO iContA=1,nContA
         iA2 = iA2+1
         iA1 = iA1+1
         DO iAngD=1,nAngD
          iD2 = (iAngD-1)*nContD
          iD1 = startD+iD2
          DO iContD=1,nContD
           iD1=iD1+1
           iD2=iD2+1
           DO iAngC=1,nAngC
            sC1 = startC + (iAngC-1)*nContC
            sC2 = (iAngC-1)*nContC
            DO iContC=1,nContC
             int = CDAB(sC2+iContC,iD2,iA2,iB2,iPass)
             DECRES(sC1+iContC,iD1,iA1,iB1,1)=int          !ABCD
             DECRES(iD1,sC1+iContC,iA1,iB1,1)=int          !ABDC
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
   END subroutine DECDISTRIBUTESUB3

  subroutine DECdistributeSUb2(DECRES,dim1,dim2,dim3,dim4,CDAB,nA,nB,nC,nD,iPass,nPasses,&
       & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
    implicit none
    integer,intent(in) :: dim1,dim2,dim3,dim4,nA,nB,nC,nD,iPass,nPasses,lupri
    integer,intent(in) :: startA,startB,startC,startD,nContA,nContB,nContC,nContD
    integer,intent(in) :: nAngA,nAngB,nAngC,nAngD
    real(realk) :: DECRES(dim3,dim4,dim1,dim2,1)
    real(realk) :: CDAB(nC,nD,nA,nB,nPasses)
    !
    integer :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
    integer :: iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,sC1,sC2
    real(realk) :: int
    DO iAngB=1,nAngB
     iB2 = (iAngB-1)*nContB
     iB1 = startB+iB2
      DO iContB=1,nContB
       iB2 = iB2+1
       iB1 = iB1+1
       DO iAngA=1,nAngA
        iA2 = (iAngA-1)*nContA
        iA1 = startA+iA2
        DO iContA=1,nContA
         iA2 = iA2+1
         iA1 = iA1+1
         DO iAngD=1,nAngD
          iD2 = (iAngD-1)*nContD
          iD1 = startD+iD2
          DO iContD=1,nContD
           iD1=iD1+1
           iD2=iD2+1
           DO iAngC=1,nAngC
            sC1 = startC + (iAngC-1)*nContC
            sC2 = (iAngC-1)*nContC
            DO iContC=1,nContC
             int = CDAB(sC2+iContC,iD2,iA2,iB2,iPass)
             DECRES(sC1+iContC,iD1,iA1,iB1,1)=int          !ABCD
             DECRES(iA1,iB1,sC1+iContC,iD1,1)=int          !CDAB
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
   END subroutine DECDISTRIBUTESUB2

 END SUBROUTINE EXPLICIT4CENTERDEC

! !DEC wants the integrals in (nbast,dim3,dim4,nbast) but it is faster 
! !to calculate them as (dim3,dim4,nbast,nbast) 
! !so we calculate them as (dim3,dim4,nbast,nbast) but store the 
! !result as (nbast,dim3,dim4,nbast) 
! ! so DECRES have the dimension (dim1,dim3,dim4,dim2,1) = (nbast,dim3,dim4,nbast)
! SUBROUTINE Explicit4CenterDEC2(DECRES,dim1,dim2,dim3,dim4,&
!      & PQ,QPmat2,dimQ,dimP,Input,Lsoutput,LUPRI,IPRINT)
! implicit none 
! Type(integrand),intent(in)      :: PQ
! Type(IntegralInput),intent(in)  :: Input
! Type(IntegralOutput),intent(inout) :: Lsoutput
! Integer,intent(in)              :: LUPRI,IPRINT
! Integer,intent(in)              :: dimQ,dimP,dim1,dim2,dim3,dim4
! REAL(REALK),intent(in)          :: QPMAT2(dimQ,dimP)
! REAL(REALK),intent(inout)       :: DECRES(dim1,dim2,dim3,dim4,1)
! !
! type(overlap),pointer :: P,Q
! logical :: SamePQ,SameLHSaos,SameRHSaos,SameODs,nucleiP,nucleiQ,add
! Logical :: permuteAB,permuteCD,permuteOD,noContraction,RHScontraction
! integer :: nAngmomP,nmat,nPasses
! integer :: iPass,iPassP,iPassQ,iDeriv,iDerivP,iDerivQ,nAngmomQ
! integer :: nA,nB,nC,nD,batchA,batchB,batchC,batchD,atomA,atomB,atomC,atomD
! integer :: indADBC,indADCB,indDACB,indDABC,indBCAD,indCBAD,indBCDA,indCBDA
! Real(realk),pointer :: kADBC(:),kADCB(:),kBCAD(:),kCBAD(:),kDACB(:),kDABC(:),kBCDA(:),kCBDA(:)

! integer :: iA,iB,iC,iD,nContA,nContB,nContC,nContD
! integer :: nAngA,nAngB,nAngC,nAngD,iAngmomP,iAngmomQ
! integer :: startA,startB,startC,startD

! nAngmomP       = PQ%P%p%nAngmom
! nAngmomQ       = PQ%Q%p%nAngmom
! SamePQ         = PQ%samePQ
! SameLHSaos     = INPUT%SameLHSaos
! SameRHSaos     = INPUT%SameRHSaos
! SameODs        = INPUT%SameODs
! nmat           = 1

! P => PQ%P%p
! Q => PQ%Q%p

! nPasses = P%nPasses*Q%nPasses
! DO iAngmomP=1,nAngmomP
!  iA = P%indexAng1(iAngmomP)
!  iB = P%indexAng2(iAngmomP)
!  nContA = P%orbital1%nContracted(iA)
!  nContB = P%orbital2%nContracted(iB)
!  nAngA = P%orbital1%nOrbComp(iA)
!  nAngB = P%orbital2%nOrbComp(iB)
!  DO iPassP=1,P%nPasses
!   startA  = P%orbital1%startOrbital(iA+(iPassP-1)*P%orbital1%nAngmom) - 1
!   startB  = P%orbital2%startOrbital(iB+(iPassP-1)*P%orbital2%nAngmom) - 1    

!   atomA  = P%orb1atom(iPassP)
!   batchA = P%orb1batch(iPassP)
!   nA     = P%orbital1%totOrbitals
!   atomB  = P%orb2atom(iPassP)
!   batchB = P%orb2batch(iPassP)
!   nB     = P%orbital2%totOrbitals
!   permuteAB  = SameLHSaos.AND.((batchA.NE.batchB).OR.(atomA.NE.atomB))
!   DO iAngmomQ=1,nAngmomQ
!    iC = Q%indexAng1(iAngmomQ)
!    iD = Q%indexAng2(iAngmomQ)
!    nContC = Q%orbital1%nContracted(iC)
!    nContD = Q%orbital2%nContracted(iD)
!    nAngC = Q%orbital1%nOrbComp(iC)
!    nAngD = Q%orbital2%nOrbComp(iD)
!    DO iPassQ=1,Q%nPasses
!     iPass = iPassQ + (iPassP-1)*Q%nPasses
!     startC  = Q%orbital1%startOrbital(iC+(iPassQ-1)*Q%orbital1%nAngmom) - 1
!     startD  = Q%orbital2%startOrbital(iD+(iPassQ-1)*Q%orbital2%nAngmom) - 1    
!     atomC  = Q%orb1atom(iPassQ)
!     batchC = Q%orb1batch(iPassQ)
!     nC     = Q%orbital1%totOrbitals
!     atomD  = Q%orb2atom(iPassQ)
!     batchD = Q%orb2batch(iPassQ)
!     nD     = Q%orbital2%totOrbitals
!     permuteCD = SameRHSaos.AND.((batchC.NE.batchD).OR.(atomC.NE.atomD))
!     permuteOD = SameODs.AND.( ((batchA.NE.batchC).OR.(atomA.NE.atomC)).OR.&
!    &                          ((batchB.NE.batchD).OR.(atomB.NE.atomD)) )
!     IF (permuteOD.AND.permuteCD.AND.permuteAB) THEN
!        call DECdistributeSUbFULL(DECRES,dim1,dim2,dim3,dim4,QPmat2,nA,nB,nC,nD,iPass,nPasses,&
!             & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
!     ELSE IF (permuteOD.AND.permuteCD) THEN
!        call DECdistributeSUb7(DECRES,dim1,dim2,dim3,dim4,QPmat2,nA,nB,nC,nD,iPass,nPasses,&
!             & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
!     ELSE IF (permuteOD.AND.permuteAB) THEN
!        call DECdistributeSUb6(DECRES,dim1,dim2,dim3,dim4,QPmat2,nA,nB,nC,nD,iPass,nPasses,&
!             & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
!     ELSE IF (permuteAB.AND.permuteCD) THEN
!        call DECdistributeSUb5(DECRES,dim1,dim2,dim3,dim4,QPmat2,nA,nB,nC,nD,iPass,nPasses,&
!             & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
!     ELSE IF (permuteAB) THEN
!        call DECdistributeSUb4(DECRES,dim1,dim2,dim3,dim4,QPmat2,nA,nB,nC,nD,iPass,nPasses,&
!             & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
!     ELSE IF (permuteCD) THEN
!        call DECdistributeSUb3(DECRES,dim1,dim2,dim3,dim4,QPmat2,nA,nB,nC,nD,iPass,nPasses,&
!             & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
!     ELSE IF (permuteOD) THEN
!        call DECdistributeSUb2(DECRES,dim1,dim2,dim3,dim4,QPmat2,nA,nB,nC,nD,iPass,nPasses,&
!             & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
!     ELSE
!        call DECdistributeSUb1(DECRES,dim1,dim2,dim3,dim4,QPmat2,nA,nB,nC,nD,iPass,nPasses,&
!             & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
!     ENDIF
!    ENDDO !iPassQ
!   ENDDO !iAngmomQ
!  ENDDO !iPassP
! ENDDO !iAngmomP

! CONTAINS
!   subroutine DECdistributeSUb1(DECRES,dim1,dim2,dim3,dim4,CDAB,nA,nB,nC,nD,iPass,nPasses,&
!        & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
!     implicit none
!     integer,intent(in) :: dim1,dim2,dim3,dim4,nA,nB,nC,nD,iPass,nPasses,lupri
!     integer,intent(in) :: startA,startB,startC,startD,nContA,nContB,nContC,nContD
!     integer,intent(in) :: nAngA,nAngB,nAngC,nAngD
!     real(realk) :: DECRES(dim1,dim2,dim3,dim4,1)
!     real(realk) :: CDAB(nC,nD,nA,nB,nPasses)
!     !
!     integer :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
!     integer :: iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,sC1,sC2
!     real(realk) :: int
!     print*,'DECRES(dim1,dim2,dim3,dim4,1)',dim1,dim2,dim3,dim4
!     DO iAngB=1,nAngB
!      iB2 = (iAngB-1)*nContB
!      iB1 = startB+iB2
!       DO iContB=1,nContB
!        iB2 = iB2+1
!        iB1 = iB1+1
!        DO iAngA=1,nAngA
!         iA2 = (iAngA-1)*nContA
!         iA1 = startA+iA2
!         DO iContA=1,nContA
!          iA2 = iA2+1
!          iA1 = iA1+1
!          DO iAngD=1,nAngD
!           iD2 = (iAngD-1)*nContD
!           iD1 = startD+iD2
!           DO iContD=1,nContD
!            iD1=iD1+1
!            iD2=iD2+1
!            DO iAngC=1,nAngC
!             sC1 = startC + (iAngC-1)*nContC
!             sC2 = (iAngC-1)*nContC
!             DO iContC=1,nContC
!              DECRES(sC1+iContC,iA1,iB1,iD1,1)=CDAB(sC2+iContC,iD2,iA2,iB2,iPass)
!             ENDDO
!            ENDDO
!           ENDDO
!          ENDDO
!         ENDDO
!        ENDDO
!       ENDDO
!      ENDDO
!    END subroutine DECDISTRIBUTESUB1

!   subroutine DECdistributeSUb3(DECRES,dim1,dim2,dim3,dim4,CDAB,nA,nB,nC,nD,iPass,nPasses,&
!        & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
!     implicit none
!     integer,intent(in) :: dim1,dim2,dim3,dim4,nA,nB,nC,nD,iPass,nPasses,lupri
!     integer,intent(in) :: startA,startB,startC,startD,nContA,nContB,nContC,nContD
!     integer,intent(in) :: nAngA,nAngB,nAngC,nAngD
!     real(realk) :: DECRES(dim1,dim2,dim3,dim4,1)
!     real(realk) :: CDAB(nC,nD,nA,nB,nPasses)
!     !
!     integer :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
!     integer :: iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,sC1,sC2
!     real(realk) :: int
    
!     DO iAngB=1,nAngB
!      iB2 = (iAngB-1)*nContB
!      iB1 = startB+iB2
!       DO iContB=1,nContB
!        iB2 = iB2+1
!        iB1 = iB1+1
!        DO iAngA=1,nAngA
!         iA2 = (iAngA-1)*nContA
!         iA1 = startA+iA2
!         DO iContA=1,nContA
!          iA2 = iA2+1
!          iA1 = iA1+1
!          DO iAngD=1,nAngD
!           iD2 = (iAngD-1)*nContD
!           iD1 = startD+iD2
!           DO iContD=1,nContD
!            iD1=iD1+1
!            iD2=iD2+1
!            DO iAngC=1,nAngC
!             sC1 = startC + (iAngC-1)*nContC
!             sC2 = (iAngC-1)*nContC
!             DO iContC=1,nContC
!              int = CDAB(sC2+iContC,iD2,iA2,iB2,iPass)
!              DECRES(sC1+iContC,iA1,iB1,iD1,1)=int          !ABCD
!              DECRES(iD1,iA1,iB1,sC1+iContC,1)=int          !ABDC
!             ENDDO
!            ENDDO
!           ENDDO
!          ENDDO
!         ENDDO
!        ENDDO
!       ENDDO
!      ENDDO
!    END subroutine DECDISTRIBUTESUB3

!   subroutine DECdistributeSUbFULL(DECRES,dim1,dim2,dim3,dim4,CDAB,nA,nB,nC,nD,iPass,nPasses,&
!        & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
!     implicit none
!     integer,intent(in) :: dim1,dim2,dim3,dim4,nA,nB,nC,nD,iPass,nPasses,lupri
!     integer,intent(in) :: startA,startB,startC,startD,nContA,nContB,nContC,nContD
!     integer,intent(in) :: nAngA,nAngB,nAngC,nAngD
!     real(realk) :: DECRES(dim1,dim2,dim3,dim4,1)
!     real(realk) :: CDAB(nC,nD,nA,nB,nPasses)
!     !
!     integer :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
!     integer :: iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,sC1,sC2
!     real(realk) :: int
    
!     DO iAngB=1,nAngB
!      iB2 = (iAngB-1)*nContB
!      iB1 = startB+iB2
!       DO iContB=1,nContB
!        iB2 = iB2+1
!        iB1 = iB1+1
!        DO iAngA=1,nAngA
!         iA2 = (iAngA-1)*nContA
!         iA1 = startA+iA2
!         DO iContA=1,nContA
!          iA2 = iA2+1
!          iA1 = iA1+1
!          DO iAngD=1,nAngD
!           iD2 = (iAngD-1)*nContD
!           iD1 = startD+iD2
!           DO iContD=1,nContD
!            iD1=iD1+1
!            iD2=iD2+1
!            DO iAngC=1,nAngC
!             sC1 = startC + (iAngC-1)*nContC
!             sC2 = (iAngC-1)*nContC
!             DO iContC=1,nContC
!              int = CDAB(sC2+iContC,iD2,iA2,iB2,iPass)
!              DECRES(sC1+iContC,iA1,iB1,iD1,1)=int          !ABCD
!              DECRES(sC1+iContC,iB1,iA1,iD1,1)=int          !BACD
!              DECRES(iD1,iA1,iB1,sC1+iContC,1)=int          !ABDC
!              DECRES(iD1,iB1,iA1,sC1+iContC,1)=int          !BADC
!              DECRES(iA1,sC1+iContC,iD1,iB1,1)=int          !CDAB
!              DECRES(iA1,iD1,sC1+iContC,iB1,1)=int          !DCAB
!              DECRES(iB1,sC1+iContC,iD1,iA1,1)=int          !CDBA
!              DECRES(iB1,iD1,sC1+iContC,iA1,1)=int          !DCBA
!             ENDDO
!            ENDDO
!           ENDDO
!          ENDDO
!         ENDDO
!        ENDDO
!       ENDDO
!      ENDDO
!    END subroutine DECDISTRIBUTESUBFULL

!   subroutine DECdistributeSUb7(DECRES,dim1,dim2,dim3,dim4,CDAB,nA,nB,nC,nD,iPass,nPasses,&
!        & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
!     implicit none
!     integer,intent(in) :: dim1,dim2,dim3,dim4,nA,nB,nC,nD,iPass,nPasses,lupri
!     integer,intent(in) :: startA,startB,startC,startD,nContA,nContB,nContC,nContD
!     integer,intent(in) :: nAngA,nAngB,nAngC,nAngD
!     real(realk) :: DECRES(dim1,dim2,dim3,dim4,1)
!     real(realk) :: CDAB(nC,nD,nA,nB,nPasses)
!     !
!     integer :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
!     integer :: iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,sC1,sC2
!     real(realk) :: int
    
!     DO iAngB=1,nAngB
!      iB2 = (iAngB-1)*nContB
!      iB1 = startB+iB2
!       DO iContB=1,nContB
!        iB2 = iB2+1
!        iB1 = iB1+1
!        DO iAngA=1,nAngA
!         iA2 = (iAngA-1)*nContA
!         iA1 = startA+iA2
!         DO iContA=1,nContA
!          iA2 = iA2+1
!          iA1 = iA1+1
!          DO iAngD=1,nAngD
!           iD2 = (iAngD-1)*nContD
!           iD1 = startD+iD2
!           DO iContD=1,nContD
!            iD1=iD1+1
!            iD2=iD2+1
!            DO iAngC=1,nAngC
!             sC1 = startC + (iAngC-1)*nContC
!             sC2 = (iAngC-1)*nContC
!             DO iContC=1,nContC
!              int = CDAB(sC2+iContC,iD2,iA2,iB2,iPass)
!              DECRES(sC1+iContC,iA1,iB1,iD1,1)=int          !ABCD
!              DECRES(iD1,iA1,iB1,sC1+iContC,1)=int          !ABDC
!              DECRES(iA1,sC1+iContC,iD1,iB1,1)=int          !CDAB
!              DECRES(iA1,iD1,sC1+iContC,iB1,1)=int          !DCAB
!             ENDDO
!            ENDDO
!           ENDDO
!          ENDDO
!         ENDDO
!        ENDDO
!       ENDDO
!      ENDDO
!    END subroutine DECDISTRIBUTESUB7

!   subroutine DECdistributeSUb6(DECRES,dim1,dim2,dim3,dim4,CDAB,nA,nB,nC,nD,iPass,nPasses,&
!        & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
!     implicit none
!     integer,intent(in) :: dim1,dim2,dim3,dim4,nA,nB,nC,nD,iPass,nPasses,lupri
!     integer,intent(in) :: startA,startB,startC,startD,nContA,nContB,nContC,nContD
!     integer,intent(in) :: nAngA,nAngB,nAngC,nAngD
!     real(realk) :: DECRES(dim1,dim2,dim3,dim4,1)
!     real(realk) :: CDAB(nC,nD,nA,nB,nPasses)
!     !
!     integer :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
!     integer :: iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,sC1,sC2
!     real(realk) :: int
    
!     DO iAngB=1,nAngB
!      iB2 = (iAngB-1)*nContB
!      iB1 = startB+iB2
!       DO iContB=1,nContB
!        iB2 = iB2+1
!        iB1 = iB1+1
!        DO iAngA=1,nAngA
!         iA2 = (iAngA-1)*nContA
!         iA1 = startA+iA2
!         DO iContA=1,nContA
!          iA2 = iA2+1
!          iA1 = iA1+1
!          DO iAngD=1,nAngD
!           iD2 = (iAngD-1)*nContD
!           iD1 = startD+iD2
!           DO iContD=1,nContD
!            iD1=iD1+1
!            iD2=iD2+1
!            DO iAngC=1,nAngC
!             sC1 = startC + (iAngC-1)*nContC
!             sC2 = (iAngC-1)*nContC
!             DO iContC=1,nContC
!              int = CDAB(sC2+iContC,iD2,iA2,iB2,iPass)
!              DECRES(sC1+iContC,iA1,iB1,iD1,1)=int          !ABCD
!              DECRES(sC1+iContC,iB1,iA1,iD1,1)=int          !BACD
!              DECRES(iA1,sC1+iContC,iD1,iB1,1)=int          !CDAB
!              DECRES(iB1,sC1+iContC,iD1,iA1,1)=int          !CDBA
!             ENDDO
!            ENDDO
!           ENDDO
!          ENDDO
!         ENDDO
!        ENDDO
!       ENDDO
!      ENDDO
!    END subroutine DECDISTRIBUTESUB6

!   subroutine DECdistributeSUb5(DECRES,dim1,dim2,dim3,dim4,CDAB,nA,nB,nC,nD,iPass,nPasses,&
!        & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
!     implicit none
!     integer,intent(in) :: dim1,dim2,dim3,dim4,nA,nB,nC,nD,iPass,nPasses,lupri
!     integer,intent(in) :: startA,startB,startC,startD,nContA,nContB,nContC,nContD
!     integer,intent(in) :: nAngA,nAngB,nAngC,nAngD
!     real(realk) :: DECRES(dim1,dim2,dim3,dim4,1)
!     real(realk) :: CDAB(nC,nD,nA,nB,nPasses)
!     !
!     integer :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
!     integer :: iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,sC1,sC2
!     real(realk) :: int
    
!     DO iAngB=1,nAngB
!      iB2 = (iAngB-1)*nContB
!      iB1 = startB+iB2
!       DO iContB=1,nContB
!        iB2 = iB2+1
!        iB1 = iB1+1
!        DO iAngA=1,nAngA
!         iA2 = (iAngA-1)*nContA
!         iA1 = startA+iA2
!         DO iContA=1,nContA
!          iA2 = iA2+1
!          iA1 = iA1+1
!          DO iAngD=1,nAngD
!           iD2 = (iAngD-1)*nContD
!           iD1 = startD+iD2
!           DO iContD=1,nContD
!            iD1=iD1+1
!            iD2=iD2+1
!            DO iAngC=1,nAngC
!             sC1 = startC + (iAngC-1)*nContC
!             sC2 = (iAngC-1)*nContC
!             DO iContC=1,nContC
!              int = CDAB(sC2+iContC,iD2,iA2,iB2,iPass)
!              DECRES(sC1+iContC,iA1,iB1,iD1,1)=int          !ABCD
!              DECRES(iD1,iA1,iB1,sC1+iContC,1)=int          !ABDC
!              DECRES(sC1+iContC,iB1,iA1,iD1,1)=int          !BACD
!              DECRES(iD1,iB1,iA1,sC1+iContC,1)=int          !BADC
!             ENDDO
!            ENDDO
!           ENDDO
!          ENDDO
!         ENDDO
!        ENDDO
!       ENDDO
!      ENDDO
!    END subroutine DECDISTRIBUTESUB5

!   subroutine DECdistributeSUb4(DECRES,dim1,dim2,dim3,dim4,CDAB,nA,nB,nC,nD,iPass,nPasses,&
!        & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
!     implicit none
!     integer,intent(in) :: dim1,dim2,dim3,dim4,nA,nB,nC,nD,iPass,nPasses,lupri
!     integer,intent(in) :: startA,startB,startC,startD,nContA,nContB,nContC,nContD
!     integer,intent(in) :: nAngA,nAngB,nAngC,nAngD
!     real(realk) :: DECRES(dim1,dim2,dim3,dim4,1)
!     real(realk) :: CDAB(nC,nD,nA,nB,nPasses)
!     !
!     integer :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
!     integer :: iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,sC1,sC2
!     real(realk) :: int
    
!     DO iAngB=1,nAngB
!      iB2 = (iAngB-1)*nContB
!      iB1 = startB+iB2
!       DO iContB=1,nContB
!        iB2 = iB2+1
!        iB1 = iB1+1
!        DO iAngA=1,nAngA
!         iA2 = (iAngA-1)*nContA
!         iA1 = startA+iA2
!         DO iContA=1,nContA
!          iA2 = iA2+1
!          iA1 = iA1+1
!          DO iAngD=1,nAngD
!           iD2 = (iAngD-1)*nContD
!           iD1 = startD+iD2
!           DO iContD=1,nContD
!            iD1=iD1+1
!            iD2=iD2+1
!            DO iAngC=1,nAngC
!             sC1 = startC + (iAngC-1)*nContC
!             sC2 = (iAngC-1)*nContC
!             DO iContC=1,nContC
!              int = CDAB(sC2+iContC,iD2,iA2,iB2,iPass)
!              DECRES(sC1+iContC,iA1,iB1,iD1,1)=int          !ABCD
!              DECRES(sC1+iContC,iB1,iA1,iD1,1)=int          !BACD
!             ENDDO
!            ENDDO
!           ENDDO
!          ENDDO
!         ENDDO
!        ENDDO
!       ENDDO
!      ENDDO
!    END subroutine DECDISTRIBUTESUB4

!   subroutine DECdistributeSUb2(DECRES,dim1,dim2,dim3,dim4,CDAB,nA,nB,nC,nD,iPass,nPasses,&
!        & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
!     implicit none
!     integer,intent(in) :: dim1,dim2,dim3,dim4,nA,nB,nC,nD,iPass,nPasses,lupri
!     integer,intent(in) :: startA,startB,startC,startD,nContA,nContB,nContC,nContD
!     integer,intent(in) :: nAngA,nAngB,nAngC,nAngD
!     real(realk) :: DECRES(dim1,dim2,dim3,dim4,1)
!     real(realk) :: CDAB(nC,nD,nA,nB,nPasses)
!     !
!     integer :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
!     integer :: iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,sC1,sC2
!     real(realk) :: int
!     DO iAngB=1,nAngB
!      iB2 = (iAngB-1)*nContB
!      iB1 = startB+iB2
!       DO iContB=1,nContB
!        iB2 = iB2+1
!        iB1 = iB1+1
!        DO iAngA=1,nAngA
!         iA2 = (iAngA-1)*nContA
!         iA1 = startA+iA2
!         DO iContA=1,nContA
!          iA2 = iA2+1
!          iA1 = iA1+1
!          DO iAngD=1,nAngD
!           iD2 = (iAngD-1)*nContD
!           iD1 = startD+iD2
!           DO iContD=1,nContD
!            iD1=iD1+1
!            iD2=iD2+1
!            DO iAngC=1,nAngC
!             sC1 = startC + (iAngC-1)*nContC
!             sC2 = (iAngC-1)*nContC
!             DO iContC=1,nContC
!              int = CDAB(sC2+iContC,iD2,iA2,iB2,iPass)
!              DECRES(sC1+iContC,iA1,iB1,iD1,1)=int          !ABCD
!              DECRES(iA1,sC1+iContC,iD1,iB1,1)=int          !CDAB
!             ENDDO
!            ENDDO
!           ENDDO
!          ENDDO
!         ENDDO
!        ENDDO
!       ENDDO
!      ENDDO
!    END subroutine DECDISTRIBUTESUB2

!  END SUBROUTINE EXPLICIT4CENTERDEC2

 !DEC wants the decpacked K integrals directly in the output (not lstensor)
 SUBROUTINE Explicit4CenterDECK(DECRES,dim1,dim2,dim3,dim4,&
      & PQ,QPmat2,dimQ,dimP,Input,Lsoutput,LUPRI,IPRINT)
   implicit none 
   Type(integrand),intent(in)      :: PQ
   Type(IntegralInput),intent(in)  :: Input
   Type(IntegralOutput),intent(inout) :: Lsoutput
   Integer,intent(in)              :: LUPRI,IPRINT
   Integer,intent(in)              :: dimQ,dimP,dim1,dim2,dim3,dim4
   REAL(REALK),intent(in)          :: QPMAT2(dimQ,dimP)
   REAL(REALK),intent(inout)       :: DECRES(dim1,dim2,dim3,dim4,1)
   !
   type(overlap),pointer :: P,Q
   logical :: SamePQ,SameLHSaos,SameRHSaos,SameODs,nucleiP,nucleiQ,add
   Logical :: permuteAB,permuteCD,permuteOD,noContraction,RHScontraction
   integer :: nAngmomP,nmat,nPasses
   integer :: iPass,iPassP,iPassQ,iDeriv,iDerivP,iDerivQ,nAngmomQ
   integer :: nA,nB,nC,nD,batchA,batchB,batchC,batchD,atomA,atomB,atomC,atomD
   integer :: indADBC,indADCB,indDACB,indDABC,indBCAD,indCBAD,indBCDA,indCBDA
   Real(realk),pointer :: kADBC(:),kADCB(:),kBCAD(:),kCBAD(:),kDACB(:),kDABC(:),kBCDA(:),kCBDA(:)

   integer :: iA,iB,iC,iD,nContA,nContB,nContC,nContD
   integer :: nAngA,nAngB,nAngC,nAngD,iAngmomP,iAngmomQ
   integer :: startA,startB,startC,startD

   nAngmomP       = PQ%P%p%nAngmom
   nAngmomQ       = PQ%Q%p%nAngmom
   SamePQ         = PQ%samePQ
   SameODs        = INPUT%SameODs
   IF(INPUT%SameLHSaos)call lsquit('same LHS aos in DECK',-1)
   IF(INPUT%SameRHSaos)call lsquit('same RHS aos in DECK',-1)
   nmat           = 1

   P => PQ%P%p
   Q => PQ%Q%p

   nPasses = P%nPasses*Q%nPasses
   DO iAngmomP=1,nAngmomP
      iA = P%indexAng1(iAngmomP)
      iB = P%indexAng2(iAngmomP)
      nContA = P%orbital1%nContracted(iA)
      nContB = P%orbital2%nContracted(iB)
      nAngA = P%orbital1%nOrbComp(iA)
      nAngB = P%orbital2%nOrbComp(iB)
      nA     = P%orbital1%totOrbitals
      nB     = P%orbital2%totOrbitals
      nC     = Q%orbital1%totOrbitals
      nD     = Q%orbital2%totOrbitals
      DO iPassP=1,P%nPasses
         startA  = P%orbital1%startOrbital(iA+(iPassP-1)*P%orbital1%nAngmom) - 1
         startB  = P%orbital2%startOrbital(iB+(iPassP-1)*P%orbital2%nAngmom) - 1    

         atomA  = P%orb1atom(iPassP)
         batchA = P%orb1batch(iPassP)
         atomB  = P%orb2atom(iPassP)
         batchB = P%orb2batch(iPassP)
         DO iAngmomQ=1,nAngmomQ
            iC = Q%indexAng1(iAngmomQ)
            iD = Q%indexAng2(iAngmomQ)
            nContC = Q%orbital1%nContracted(iC)
            nContD = Q%orbital2%nContracted(iD)
            nAngC = Q%orbital1%nOrbComp(iC)
            nAngD = Q%orbital2%nOrbComp(iD)
            DO iPassQ=1,Q%nPasses
               iPass = iPassQ + (iPassP-1)*Q%nPasses
               startC  = Q%orbital1%startOrbital(iC+(iPassQ-1)*Q%orbital1%nAngmom) - 1
               startD  = Q%orbital2%startOrbital(iD+(iPassQ-1)*Q%orbital2%nAngmom) - 1    
               atomC  = Q%orb1atom(iPassQ)
               batchC = Q%orb1batch(iPassQ)
               atomD  = Q%orb2atom(iPassQ)
               batchD = Q%orb2batch(iPassQ)
               permuteOD = SameODs.AND.( ((batchA.NE.batchC).OR.(atomA.NE.atomC)).OR.&
                    &                          ((batchB.NE.batchD).OR.(atomB.NE.atomD)) )
               IF (permuteOD) THEN
                  call DECKdistributeSUb2(DECRES,dim1,dim2,dim3,dim4,QPmat2,nA,nB,nC,nD,iPass,nPasses,&
                       & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
               ELSE
                  call DECKdistributeSUb1(DECRES,dim1,dim2,dim3,dim4,QPmat2,nA,nB,nC,nD,iPass,nPasses,&
                       & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
               ENDIF
            ENDDO !iPassQ
         ENDDO !iAngmomQ
      ENDDO !iPassP
   ENDDO !iAngmomP

 CONTAINS
   subroutine DECKdistributeSUb2(DECRES,dim1,dim2,dim3,dim4,CDAB,nA,nB,nC,nD,iPass,nPasses,&
        & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
     implicit none
     integer,intent(in) :: dim1,dim2,dim3,dim4,nA,nB,nC,nD,iPass,nPasses,lupri
     integer,intent(in) :: startA,startB,startC,startD,nContA,nContB,nContC,nContD
     integer,intent(in) :: nAngA,nAngB,nAngC,nAngD
     real(realk) :: DECRES(dim1,dim2,dim3,dim4,1)
     real(realk) :: CDAB(nC,nD,nA,nB,nPasses)
     !
     integer :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
     integer :: iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,sC1,sC2
     real(realk) :: int
     DO iAngB=1,nAngB
        iB2 = (iAngB-1)*nContB
        iB1 = startB+iB2
        DO iContB=1,nContB
           iB2 = iB2+1
           iB1 = iB1+1
           DO iAngA=1,nAngA
              iA2 = (iAngA-1)*nContA
              iA1 = startA+iA2
              DO iContA=1,nContA
                 iA2 = iA2+1
                 iA1 = iA1+1
                 DO iAngD=1,nAngD
                    iD2 = (iAngD-1)*nContD
                    iD1 = startD+iD2
                    DO iContD=1,nContD
                       iD1=iD1+1
                       iD2=iD2+1
                       DO iAngC=1,nAngC
                          sC1 = startC + (iAngC-1)*nContC
                          sC2 = (iAngC-1)*nContC
                          DO iContC=1,nContC
                             int = CDAB(sC2+iContC,iD2,iA2,iB2,iPass)
                             DECRES(iA1,iB1,sC1+iContC,iD1,1)=int          !ABCD
                             DECRES(sC1+iContC,iD1,iA1,iB1,1)=int          !CDAB
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
   END subroutine DECKDISTRIBUTESUB2

   subroutine DECKdistributeSUb1(DECRES,dim1,dim2,dim3,dim4,CDAB,nA,nB,nC,nD,iPass,nPasses,&
        & startA,startB,startC,startD,nContA,nContB,nContC,nContD,nAngA,nAngB,nAngC,nAngD,lupri)
     implicit none
     integer,intent(in) :: dim1,dim2,dim3,dim4,nA,nB,nC,nD,iPass,nPasses,lupri
     integer,intent(in) :: startA,startB,startC,startD,nContA,nContB,nContC,nContD
     integer,intent(in) :: nAngA,nAngB,nAngC,nAngD
     real(realk) :: DECRES(dim1,dim2,dim3,dim4,1)
     real(realk) :: CDAB(nC,nD,nA,nB,nPasses)
     !
     integer :: iContA,iContB,iContC,iContD,iAngA,iAngB,iAngC,iAngD
     integer :: iA1,iA2,iB1,iB2,iC1,iC2,iD1,iD2,sC1,sC2
     real(realk) :: int

     DO iAngB=1,nAngB
        iB2 = (iAngB-1)*nContB
        iB1 = startB+iB2
        DO iContB=1,nContB
           iB2 = iB2+1
           iB1 = iB1+1
           DO iAngA=1,nAngA
              iA2 = (iAngA-1)*nContA
              iA1 = startA+iA2
              DO iContA=1,nContA
                 iA2 = iA2+1
                 iA1 = iA1+1
                 DO iAngD=1,nAngD
                    iD2 = (iAngD-1)*nContD
                    iD1 = startD+iD2
                    DO iContD=1,nContD
                       iD1=iD1+1
                       iD2=iD2+1
                       DO iAngC=1,nAngC
                          sC1 = startC + (iAngC-1)*nContC
                          sC2 = (iAngC-1)*nContC
                          DO iContC=1,nContC
                             DECRES(iA1,iB1,sC1+iContC,iD1,1)=CDAB(sC2+iContC,iD2,iA2,iB2,iPass)
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
   END subroutine DECKDISTRIBUTESUB1

 END SUBROUTINE EXPLICIT4CENTERDECK

END MODULE thermite_distributeDEC


