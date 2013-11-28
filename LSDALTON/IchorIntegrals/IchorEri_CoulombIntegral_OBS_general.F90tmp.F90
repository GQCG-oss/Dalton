!> @file
!> Contains the general McM driver 

!> \brief General Obara Saika Scheme Integral scheme
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorEriCoulombintegralOBSGeneralMod
use precision
use ThermiteMem_module
use memory_handling
private 
public :: IchorCoulombIntegral_OBS_general
CONTAINS
  subroutine IchorCoulombIntegral_OBS_general(nPrimA,nPrimB,nPrimC,nPrimD,&
       & nPrimP,nPrimQ,nPrimPQ,nPasses,MaxPasses,IntPrint,lupri,&
       & nContA,nContB,nContC,nContD,nContP,nContQ,pexp,qexp,ACC,BCC,CCC,DCC,&
       & pcent,qcent,Ppreexpfac,Qpreexpfac,nTABFJW1,nTABFJW2,TABFJW,&
       & Qiprim1,Qiprim2,Piprim1,Piprim2,Aexp,Bexp,Cexp,Dexp,&
       & Qsegmented,Psegmented,reducedExponents,integralPrefactor,&
       & AngmomA,AngmomB,AngmomC,AngmomD,Pdistance12,Qdistance12,PQorder,CDAB,&
       & Acenter,Bcenter,Ccenter,Dcenter)
    implicit none
    integer,intent(in) :: nPrimQ,nPrimP,nPasses,nPrimA,nPrimB,nPrimC,nPrimD
    integer,intent(in) :: nPrimPQ,MaxPasses,IntPrint,lupri
    integer,intent(in) :: nContA,nContB,nContC,nContD,nContP,nContQ,nTABFJW1,nTABFJW2
    integer,intent(in) :: Qiprim1(nPrimQ),Qiprim2(nPrimQ),Piprim1(nPrimP),Piprim2(nPrimP)
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Cexp(nPrimC),Dexp(nPrimD)
    logical,intent(in)     :: Qsegmented,Psegmented
    real(realk),intent(in) :: pexp(nPrimP),qexp(nPrimQ)
    real(realk),intent(in) :: pcent(3*nPrimP)           !qcent(3,nPrimP)
    real(realk),intent(in) :: qcent(3*nPrimQ*MaxPasses) !qcent(3,nPrimQ,MaxPasses)
    real(realk),intent(in) :: QpreExpFac(nPrimQ*MaxPasses),PpreExpFac(nPrimP)
    real(realk),intent(in) :: TABFJW(0:nTABFJW1,0:nTABFJW2)
!    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
!    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(inout) :: CDAB(:)
    real(realk),intent(in) :: integralPrefactor(nPrimPQ)
    !integralPrefactor(nPrimP,nPrimQ)
    real(realk),intent(in) :: reducedExponents(nPrimPQ)
    !reducedExponents(nPrimP,nPrimQ)
    real(realk),intent(in) :: Qdistance12(3*MaxPasses) !Ccenter-Dcenter
    !Qdistance12(3,MaxPasses)
    real(realk),intent(in) :: Pdistance12(3)           !Acenter-Bcenter 
    real(realk),pointer :: squaredDistance(:),Rpq(:),Rpa(:),Rqc(:)
    real(realk) :: Acenter(3),Bcenter(3),Ccenter(3),Dcenter(3)
    !
!    real(realk),target :: TMPWORK((2+2*nPasses)*nPrimP*nPrimQ)
    integer :: nPrimPassQ,AngmomP,AngmomQ,AngmomPQ,nTUV,nTUVQ
    logical :: PQorder
    integer :: ijk1,ijk2,ijkPcart,ijk1s,ijk2s,ijkPsph
    integer :: ijk3,ijk4,ijkQcart,ijk3s,ijk4s,ijkQsph,nPassesP,nTUVP
    logical :: Sph1,Sph2,sphericalGTO,Sph3,Sph4,SphericalTransP,SphericalTransQ
    integer :: AngmomID,currentMaxAngmom
    real(realk),pointer :: RJ000(:),WTUV(:),RE(:),EcoeffN(:),ERE(:),SERE(:),SERES(:),OUTPUTinterest(:),AUXarray(:,:,:)
    real(realk),pointer :: AUXarray2(:),Theta(:),AuxArrayCont(:),FINALOUT(:)
    integer :: TUV,J,T,U,V,I,nsize,LimitAngmom,la,lb,lc,ld,currentMaxAngmomP
    integer :: nContPQ,nTUVA,nTUVB,nTUVAspec,nTUVBspec,ijkA,ijkB
    logical :: RHS

    call interest_initialize()
    !A(1) center
    !an CC for A
    print*,'3*nPrimP',3*nPrimP
    call mem_alloc(Rpa,3*nPrimP)
    call build_Rpa(nPrimP,Pcent,Acenter,Rpa)
    print*,'3*nPrimQ',3*nPrimQ
    call mem_alloc(Rqc,3*nPrimQ)
    call build_Rpa(nPrimQ,Qcent,Ccenter,Rqc)

!    call mem_alloc(Rwp,3*nPrimQ*nPrimC)
!    call build_Rpa(nPrimQ,Qcent,Qcenter,Rwp)
    nsize = size(CDAB)*10
    print*,'nsize',nsize
    call mem_alloc(OUTPUTinterest,nsize)
    nsize = nPrimPQ
    print*,'AngmomA',AngmomA,'ACC',ACC
    print*,'AngmomB',AngmomB,'BCC',BCC
    print*,'AngmomC',AngmomC,'CCC',CCC
    print*,'AngmomD',AngmomD,'DCC',DCC
    la = AngmomA+1
    lb = AngmomB+1
    lc = AngmomC+1
    ld = AngmomD+1
    call interest_eri(OUTPUTinterest,nsize,&
         & la,Aexp,Acenter(1),Acenter(2),Acenter(3),ACC,&
         & lb,Bexp,Bcenter(1),Bcenter(2),Bcenter(3),BCC,&
         & lc,Cexp,Ccenter(1),Ccenter(2),Ccenter(3),CCC,&
         & ld,Dexp,Dcenter(1),Dcenter(2),Dcenter(3),DCC,&
         & lupri)!,&
!         & .false.)
    
    print*,'AFTER AngmomA',AngmomA,'ACC',ACC
    print*,'AFTER AngmomB',AngmomB,'BCC',BCC
    print*,'AFTER AngmomC',AngmomC,'CCC',CCC
    print*,'AFTER AngmomD',AngmomD,'DCC',DCC

    write(lupri,*)'OUTPUTinterest',OUTPUTinterest

!    write(lupri,*)'b*XAB/q',Bexp(1)*Pdistance12(1)/qexp(1)
!    write(lupri,*)'d*XCD/q',Dexp(1)*Qdistance12(1)/qexp(1)
!    write(lupri,*)'b*XAB/q+d*XCD/q',Bexp(1)*Pdistance12(1)/qexp(1)+Dexp(1)*Qdistance12(1)/qexp(1)
!    write(lupri,*)'Xqc + p/q*Xpa',Rqc(1) + Pexp(1)/Qexp(1)*Rpa(1)

    sphericalGTO = .TRUE.
    IF(PQorder)THEN
       call lsquit('PQorder OBS general expect to get QP ordering',-1)
    ENDIF
    nPrimPassQ = nPrimQ*nPasses 
    call mem_alloc(squaredDistance,nPasses*nPrimPQ)
    call mem_alloc(Rpq,3*nPasses*nPrimPQ)
    !builds squaredDistance(nPrimQ,nPrimP,nPassQ)
    !builds Rpq(3,nPrimQ,nPrimP,nPassQ)
    call build_QP_squaredDistance_and_Rpq(nPrimQ,nPasses,nPrimP,Qcent,Pcent,&
         & squaredDistance,Rpq)
    
    AngmomP = AngmomA+AngmomB
    AngmomQ = AngmomC+AngmomD
    AngmomPQ  = AngmomP + AngmomQ
    LimitAngmom = AngmomPQ
    call mem_alloc(RJ000,(LimitAngmom+1)*nPasses*nPrimPQ)
    call buildRJ000_general(nPasses,nPrimPQ,nTABFJW1,nTABFJW2,reducedExponents,&
         & squaredDistance,TABFJW,RJ000,LimitAngmom)

    !builds RJ000(0:AngmomPQ,nPrimQ,nPrimP,nPasses)

    !IF (INTPRINT .GE. 10) THEN
       WRITE(lupri,*)'Output from W000'
       DO I=1,nPrimQ*nPrimP*nPasses
          DO J=0,LimitAngmom
             WRITE(LUPRI,'(2X,A6,I4,A1,I4,A2,ES16.8)')'RJ000(',J,',',I,')=',RJ000(1+J+(I-1)*(AngmomPQ+1))
          ENDDO
       ENDDO
    !END IF

    nTUV = (LimitAngmom+1)*(LimitAngmom+2)*(LimitAngmom+3)/6
    write(lupri,*)'nTUV=',nTUV
    call mem_alloc(AUXarray,nTUV,LimitAngmom+1,nPasses*nPrimPQ)
    !Create ThetaAux(N,1) N=(1:LimitAngmom+1) corresponding to N=(0:LimitAngmom)=(0:AngmomPQ-1)
    !not 1 = (000) 
    call scaleRJ000_general(nPasses,nPrimP,nPrimQ,reducedExponents,&
         & RJ000,LimitAngmom,IntegralPrefactor,AUXarray,nTUV,&
         & PpreExpFac,QpreExpFac,lupri)
    
    WRITE(lupri,*)'AngmomAB ',AngmomP

    !Vertical Recurrence Relations Eq. 9.10.26 of "The Book"

    !Generate i=1 Theta intermediate using Eq. 9.10.26 of "The Book"
    !ThetaAux(i+1,N) = Xpa*ThetaAux(i,N) + (-alpha/p*Xpq)*ThetaAux(i,N+1) 
    !i = 0 last 2 term vanish
    !creates  ThetaAux(2,N),ThetaAux(3,N),ThetaAux(4,N)
    !2 = (100)
    !3 = (010)
    !4 = (001)
    !which all corresponds to i=0 but using different cartesian directions
    currentMaxAngmom=LimitAngmom ! = AngmomPQ
    !CurrentMaxAngmom is the size of N 
    !we build up all ThetaAux(i+1,N) with N=1:currentMaxAngmom+1
    !Once we are done using the Vertical Recurrence Relation we only need N=0 (which corresponds to N=1)
    call VerticalRecurrence1(nPasses,nPrimP,nPrimQ,nTUV,LimitAngmom,currentMaxAngmom,&
         & reducedExponents,Pexp,Rpq,Rpa,AUXarray,lupri)
    WRITE(lupri,*)'DONE VerticalRecurrence1',currentMaxAngmom
    IF(currentMaxAngmom.EQ.1)THEN

!       nContP = nContA*nContB
!       nContQ = nContC*nContD
       nContPQ = nContP*nContQ
       call mem_alloc(AuxArray2,3*nPasses*nPrimPQ)
       call colletAux(AuxArray2,3,nPasses*nPrimPQ,AUXarray,nTUV,LimitAngmom)
       nTUVQ=1
       nTUVP=3
       nTUV=4
!       call mem_alloc(AuxArrayCont,nTUVQ*nTUV*nPasses*nContPQ)

       !Primitive contraction for all A,B,C, and D 
       IF(Qsegmented.AND.Psegmented)THEN
          !already included in preexpfac
          ACC(1,1) = 1.E0_realk
          BCC(1,1) = 1.E0_realk
          CCC(1,1) = 1.E0_realk
          DCC(1,1) = 1.E0_realk
       ENDIF
       call PrimitiveContraction(AUXarray2,CDAB,nTUVQ*nTUVP,nTUVQ*nTUV,nPrimP,nPrimQ,nPasses,&
            & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD,&
            & lupri)
       RETURN
!    ELSE
!       call lsquit('not the angmom combi testet',-1)       
    ENDIF

    !Generate i=2 Theta intermediate using Eq. 9.10.26 of "The Book"
    !ThetaAux(i+1,N) = Xpa*ThetaAux(i,N) + (-alpha/p*Xpq)*ThetaAux(i,N+1) 
    !                    + i/(2p) * (ThetaAux(i-1,N) - alpha/p*ThetaAux(i-1,N+1))
    !4 term recurrence relation
    !creates 
    !5:10 = (200,110,101,020,011,002), all coresponding to i=2
    currentMaxAngmom=currentMaxAngmom-1
    call VerticalRecurrence2(nPasses,nPrimP,nPrimQ,nTUV,LimitAngmom,currentMaxAngmom,&
         & reducedExponents,Pexp,Rpq,Rpa,AUXarray,lupri)
    WRITE(lupri,*)'DONE VerticalRecurrence2',currentMaxAngmom

    IF(currentMaxAngmom.EQ.1)THEN
       !PPSS for instance
       !       nContP = nContA*nContB
       !       nContQ = nContC*nContD
       nContPQ = nContP*nContQ
       call mem_alloc(AuxArray2,9*nPasses*nPrimPQ)
       call colletAux2(AuxArray2,9,nPasses*nPrimPQ,AUXarray,nTUV,LimitAngmom)
       nTUVQ=1
       nTUVP=9
       nTUV=10
!       call mem_alloc(AuxArrayCont,nTUVQ*nTUV*nPasses*nContPQ)

       !Primitive contraction for all A,B,C, and D 
       IF(Qsegmented.AND.Psegmented)THEN
          !already included in preexpfac
          ACC(1,1) = 1.E0_realk
          BCC(1,1) = 1.E0_realk
          CCC(1,1) = 1.E0_realk
          DCC(1,1) = 1.E0_realk
       ENDIF

       call mem_alloc(AuxArrayCont,nTUVQ*nTUVP*nPasses*nContPQ)

       call PrimitiveContraction(AUXarray2,AuxArrayCont,nTUVQ*nTUVP,nTUVQ*nTUV,nPrimP,nPrimQ,nPasses,&
            & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD,&
            & lupri)

       nTUVP = (AngmomP+1)*(AngmomP+2)*(AngmomP+3)/6
       nTUVA = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6
       nTUVB = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6
       nTUVAspec = nTUVA-AngmomA*(AngmomA+1)*(AngmomA+2)/6
       nTUVBspec = nTUVB-AngmomB*(AngmomB+1)*(AngmomB+2)/6
       WRITE(lupri,*)'nTUVP',nTUVP
       WRITE(lupri,*)'nTUVA',nTUVA
       WRITE(lupri,*)'nTUVB',nTUVB
       WRITE(lupri,*)'nTUVAspec',nTUVAspec
       WRITE(lupri,*)'nTUVBspec',nTUVBspec
       call HorizontalRR_LHS(AUXarrayCont,nTUVQ,nTUVP,nTUV,nContPQ,nPasses,Qdistance12,&
            & AngmomP,AngmomA,AngmomB,CDAB,nTUVAspec,nTUVBspec,lupri)

       RETURN
!    ELSE
!       call lsquit('not the angmom combi testet',-1)       
    ENDIF


    !Generate i=3 Theta intermediate using Eq. 9.10.26 of "The Book"
    !11:20 = (300,..,030,...,003), all coresponding to i=3
    currentMaxAngmom=currentMaxAngmom-1
    call VerticalRecurrence3(nPasses,nPrimP,nPrimQ,nTUV,LimitAngmom,currentMaxAngmom,&
         & reducedExponents,Pexp,Rpq,Rpa,AUXarray,lupri)

    !Generate i=4 Theta intermediate using Eq. 9.10.26 of "The Book"
    !21:35 all coresponding to i=4
    currentMaxAngmom=currentMaxAngmom-1
    call VerticalRecurrence4(nPasses,nPrimP,nPrimQ,nTUV,LimitAngmom,currentMaxAngmom,&
         & reducedExponents,Pexp,Rpq,Rpa,AUXarray,lupri)

    !Generate i=5 Theta intermediate using Eq. 9.10.26 of "The Book"
    !36:56 all coresponding to i=5
    currentMaxAngmom=currentMaxAngmom-1
    call VerticalRecurrence5(nPasses,nPrimP,nPrimQ,nTUV,LimitAngmom,currentMaxAngmom,&
         & reducedExponents,Pexp,Rpq,Rpa,AUXarray,lupri)

    !Generate i=6 Theta intermediate using Eq. 9.10.26 of "The Book"
    !57:84 all coresponding to i=5
    currentMaxAngmom=currentMaxAngmom-1
    call VerticalRecurrence6(nPasses,nPrimP,nPrimQ,nTUV,LimitAngmom,currentMaxAngmom,&
         & reducedExponents,Pexp,Rpq,Rpa,AUXarray,lupri)
    
    IF(currentMaxAngmom.EQ.1)THEN
       Write(lupri,*)'currentMaxAngmom',currentMaxAngmom
       Write(lupri,*)'AngmomPQ',AngmomPQ
       currentMaxAngmom = AngmomQ
       nTUVQ = (AngmomQ+1)*(AngmomQ+2)*(AngmomQ+3)/6
       Write(lupri,*)'nTUVQ',nTUVQ
       !obsolete    
       call mem_alloc(AuxArray2,nTUVQ*nTUV*nPasses*nPrimPQ)
       !Builds AuxArray2 which corresponds to i=1:AngmomPQ and k=0:0 nTUVQ=1
       call GenerateThetaFromAuxarray(nTUVQ,nTUV,nPasses*nPrimPQ,AuxArray2,Auxarray,LimitAngmom,lupri)

       !Electron Transfer Relations Eq. 9.10.27 of "The Book"
       !we only Need N=0
       !ThetaAux2(k+1,i) = (Xqc+(p/q)*Xpa)*ThetaAux2(k,i) + i/(2*q)*ThetaAux2(k,i-1) 
       !                    + k/(2q) * ThetaAux2(k-1,i) - (p/q)*ThetaAux2(k,i+1)
       !
       !Note: Xqc+(p/q)*Xpa = (-b*Xab + d*Xcd)/q
       !k=0  means special case 
       !ThetaAux2(1,i) = (Xqc+(p/q)*Xpa)*ThetaAux2(k,i) + i/(2*q)*ThetaAux2(k,i-1) 
       !                    - (p/q)*ThetaAux2(k,i+1)
       !       
       currentMaxAngmomP = AngmomP
       currentMaxAngmom = currentMaxAngmom - 1 
       WRITE(lupri,*)'TransferRecurrence1'
       WRITE(lupri,*)'currentMaxAngmomP',currentMaxAngmomP
       WRITE(lupri,*)'currentMaxAngmom',currentMaxAngmom
       call TransferRecurrence1(nPasses,nPrimP,nPrimQ,nTUVQ,nTUV,LimitAngmom,currentMaxAngmom,&
            & reducedExponents,Pexp,Qexp,Rpq,Rqc,Rpa,AUXarray2,currentMaxAngmomP,lupri)
       !NOW WE HAVE AUXarray2(nTUVQ,nTUV,nPasses*nPrimP*nPrimQ)
       !with nTUVQ and nTUVPQ but do we need all of them? 

!       nContP = nContA*nContB
!       nContQ = nContC*nContD
       nContPQ = nContP*nContQ
       !this should be a ntuvP instead of ntuv?

       nTUVP = (AngmomP+1)*(AngmomP+2)*(AngmomP+3)/6
       call mem_alloc(AuxArrayCont,nTUVQ*nTUVP*nPasses*nContPQ)

       !Primitive contraction for all A,B,C, and D 
       IF(Qsegmented.AND.Psegmented)THEN
          ACC(1,1) = 1.E0_realk
          BCC(1,1) = 1.E0_realk
          CCC(1,1) = 1.E0_realk
          DCC(1,1) = 1.E0_realk
       ENDIF
       call PrimitiveContraction(AUXarray2,AUXarrayCont,nTUVQ*nTUVP,nTUVQ*nTUV,nPrimP,nPrimQ,nPasses,&
            & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD,&
            & lupri)
       
       !Horizontal Recurrence relations Eq. 9.10.28 of "The Book"
       !Horizontal Recurrence relations Eq. 9.10.29 of "The Book"
       !avoid if Xab = 0 / Xcd = 0 
       !avoid if Xab = 0 / Xcd = 0 
       !avoid if Xab = 0 / Xcd = 0 

       !depending on the order!
       !we start by the LHS
!       SELECT CASE(AngmomP)    
       !FD case
!       call HorizontalRR_RHS_PS(AUXarrayCont,nTUV,nContPQ,nPasses,Qdistance12)
       !no need we already have angmomB=0
       !FD case
       nTUVP = (AngmomP+1)*(AngmomP+2)*(AngmomP+3)/6
       nTUVA = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6
       nTUVB = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6
       WRITE(lupri,*)'nTUVP',nTUVP
       WRITE(lupri,*)'nTUVA',nTUVA
       WRITE(lupri,*)'nTUVB',nTUVB
       nTUVAspec = nTUVA-AngmomA*(AngmomA+1)*(AngmomA+2)/6
       nTUVBspec = nTUVB-AngmomB*(AngmomB+1)*(AngmomB+2)/6
       WRITE(lupri,*)'nTUVAspec',nTUVAspec
       WRITE(lupri,*)'nTUVBspec',nTUVBspec
       call mem_alloc(THETA,nTUVQ*nTUVAspec*nTUVBspec*nContPQ*nPasses)
       
!       call HorizontalRR_RHS(AUXarrayCont,nTUVQ,nTUVP,nTUV,nContPQ,nPasses,Qdistance12,&
!       & AngmomP,AngmomA,AngmomB,THETA,nTUVAspec,nTUVBspec,lupri)
!       nTUVQ = 3 
!       call SphericalTransformation_RHS(THETA,nTUVQ,nTUVAspec,nTUVBspec,nContPQ,nPasses,&
!            & AngmomP,AngmomA,AngmomB,lupri)

       call HorizontalRR_LHS(AUXarrayCont,nTUVQ,nTUVP,nTUV,nContPQ,nPasses,Qdistance12,&
            & AngmomP,AngmomA,AngmomB,THETA,nTUVAspec,nTUVBspec,lupri)

!       WRITE(lupri,*)'THETA',THETA
       
       ijkA = 7
       ijkB = 5
!       call mem_alloc(FINALOUT,3*ijkA*ijkB*nContPQ*nPasses)
       call SphericalTransformation_LHS(THETA,nTUVQ,nTUVAspec,nTUVBspec,nContPQ*nPasses,&
!            & AngmomP,AngmomA,AngmomB,FINALOUT,ijkA,ijkB,lupri)             
            & AngmomP,AngmomA,AngmomB,CDAB,ijkA,ijkB,lupri)
       RETURN
    ELSE
       call lsquit('not the angmom combi testet',-1)       
    ENDIF

    !The input should be RJ000(0:AngmomPQ,nPrimQ,nPrimP,nPasses)
    !The Output should be Theta(nTUVPQ,nPrimPQ,nPasses) 
    !then the Theta2(nTUVQ,nTUVP,nPrimPQ,nPasses) is generated
    !then the ContTheta(nTUVQ,nTUVP,nContPQ,nPasses) or ContTheta(nTUVQ,nTUVP,nPasses)
    !then the ThetaCD(nTUVC,nTUVD,nTUVP,nContPQ,nPasses) is generated
    !then the Theta2(lmQ,nTUVP,nContPQ,nPasses) is generated
    !(special code for segmentet nCont=1)


!    call mem_alloc(AUXarray2,nTUV,LimitAngmom+1,LimitAngmom+1,nPasses*nPrimPQ)
   
!    IF( AngmomP-1.eq.6 )THEN
!       if( currentMaxAngmom.eq.1 )then

!       else

!          WRITE(lupri,*)'AngmomP-1',AngmomP-1,'.EQ.6'
!          WRITE(lupri,*)'TransferRecurrence1'
!          currentMaxAngmom=currentMaxAngmom-1
!          call TransferRecurrence1(nPasses,nPrimP,nPrimQ,nTUV,LimitAngmom,currentMaxAngmom,&
!         & reducedExponents,Pexp,Qexp,Rpq,Rqc,Rpa,AUXarray,AUXarray2,lupri)
!       endif
!    ENDIF

!    currentMaxAngmom=currentMaxAngmom-1
!    call VerticalRecurrence6(nPasses,nPrimP,nPrimQ,nTUV,LimitAngmom,currentMaxAngmom,&
!         & reducedExponents,Pexp,Rpq,Rpa,AUXarray,lupri)
    call mem_dealloc(squaredDistance)
    call mem_dealloc(Rpq)
    call mem_dealloc(RJ000)

  end subroutine IchorCoulombIntegral_OBS_general

  subroutine colletAux2(AuxArray2,ntuv9,nP,AUXarray,nTUV,LimitAngmom)
    implicit none
    integer :: ntuv9,nP,nTUV,LimitAngmom
    real(realk) :: AuxArray2(ntuv9,nP),AUXarray(nTUV,LimitAngmom+1,nP)
    !
    integer :: ip
    DO IP=1,NP
       AuxArray2(1,iP) = AUXarray(2,1,iP)
       AuxArray2(2,iP) = AUXarray(3,1,iP)
       AuxArray2(3,iP) = AUXarray(4,1,iP)
       AuxArray2(4,iP) = AUXarray(5,1,iP)
       AuxArray2(5,iP) = AUXarray(6,1,iP)
       AuxArray2(6,iP) = AUXarray(7,1,iP)
       AuxArray2(7,iP) = AUXarray(8,1,iP)
       AuxArray2(8,iP) = AUXarray(9,1,iP)
       AuxArray2(9,iP) = AUXarray(10,1,iP)
    ENDDO
  end subroutine colletAux2

  subroutine colletAux(AuxArray2,ntuv3,nP,AUXarray,nTUV,LimitAngmom)
    implicit none
    integer :: ntuv3,nP,nTUV,LimitAngmom
    real(realk) :: AuxArray2(ntuv3,nP),AUXarray(nTUV,LimitAngmom+1,nP)
         !
    integer :: ip
    DO IP=1,NP
       AuxArray2(1,iP) = AUXarray(2,1,iP)
       AuxArray2(2,iP) = AUXarray(3,1,iP)
       AuxArray2(3,iP) = AUXarray(4,1,iP)
    ENDDO
  end subroutine colletAux
       
  subroutine SphericalTransformation_LHS(THETA,nTUVQ,nTUVAspec,nTUVBspec,nContPasses,&
       & AngmomP,AngmomA,AngmomB,FINALOUT,ijkA,ijkB,lupri)
    implicit none
    integer,intent(in) :: nTUVQ,nTUVAspec,nTUVBspec,nContPasses
    integer,intent(in) :: AngmomP,AngmomA,AngmomB,lupri,ijkA,ijkB
    real(realk),intent(in) :: THETA(nTUVQ,nTUVAspec,nTUVBspec,nContPasses)
    real(realk),intent(inout) :: FINALOUT(3,ijkA,ijkB,nContPasses) 
    !
    real(realk),pointer :: SPHMATA(:,:),SPHMATB(:,:)
    integer :: PA,PNA,nsizea,nsizeb,pb,pnb,ip,iijka,iijkb,ituva,ituvb,ituvq
    PA = (angmomA+1)*(AngmomA+2)/2  !(1,3,6,10,)
    WRITE(lupri,*)'AngmomA',AngmomA
    WRITE(lupri,*)'P=',PA
    PNA = 2*angmomA+1 !(1,3,5,7,9)
    WRITE(lupri,*)'PN=',PNA
    call mem_alloc(SPHMATA,PNA,PA)
    nsizeA=PNA*PA
    SPHMATA = 0.0E0_realk
    call Build_SPHMAT23(angmomA,nSIZEA,SPHMATA) 
    
    PB = (angmomB+1)*(AngmomB+2)/2
    WRITE(lupri,*)'AngmomB',AngmomB
    WRITE(lupri,*)'P=',PB
    PNB = 2*angmomB+1
    WRITE(lupri,*)'PN=',PNB
    call mem_alloc(SPHMATB,PNB,PB)
    nsizeB=PNB*PB
    SPHMATB = 0.0E0_realk
    call Build_SPHMAT23(angmomB,nSIZEB,SPHMATB) 
    
    DO iP=1,nContPasses
     DO iijkA=1,ijkA
      DO iijkB=1,ijkB
         !warning
       DO iTUVQ=1,nTUVQ-1
        FINALOUT(iTUVQ,iijkA,iijkB,iP) = 0.0E0_realk
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    DO iP=1,nContPasses
     DO iijkA=1,ijkA
      DO iTUVA=1,nTUVAspec
       DO iijkB=1,ijkB
        DO iTUVB=1,nTUVBspec
         DO iTUVQ=1,nTUVQ-1
!warning
          FINALOUT(iTUVQ,iijkA,iijkB,iP) = FINALOUT(iTUVQ,iijkA,iijkB,iP) + &
               & THETA(iTUVQ+1,iTUVA,iTUVB,iP)*SPHMATA(iijkA,iTUVA)*SPHMATB(iijkB,iTUVB)          
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    DO iP=1,nContPasses
     DO iijkA=1,ijkA
      DO iijkB=1,ijkB
       DO iTUVQ=1,nTUVQ-1
          WRITE(lupri,'(A,I3,A,I3,A,I3,A,I3,A,ES16.8)')'FINALOUT(',iTUVQ,',',iijkA,',',iijkB,',',iP,')',&
               & FINALOUT(iTUVQ,iijkA,iijkB,iP)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    call mem_dealloc(SPHMATA)
    call mem_dealloc(SPHMATB)
  end subroutine SphericalTransformation_LHS

  !> \brief build spherical transformation matrices
  !> \author T. Kjaergaard
  !> \date 2010
  SUBROUTINE Build_SPHMAT23(L,SIZE,SPHMAT) 
    use math_fun
    IMPLICIT NONE
    INTEGER          :: nAngmom,SIZE
    REAL(realk)      :: SPHMAT(SIZE)
    !
    INTEGER          :: L,I
    Real(realk), parameter :: DM1 = -1.0_4, DO = 0.0_4, D1 = 1.0_4, D2 = 2.0_4
    INTEGER  :: M1,MADR,MABS,V0, NDER, IOFF
    INTEGER  :: EE,FF,GG,BX,BY,BZ,II,JJ,KK,PX,PY,PZ
    REAL(realk)  :: FACNRM,FAC1,FAC2, FAC3,FAC4, FACTOR
    INTEGER  :: T,U,V,A,B,C,P,Q,R,X,Y,Z,TOT,IADR,AX0,AY0,AZ0,NSIZE
    INTEGER  :: M,Ncol,Nrow,INDEX,NDIM,STARTINDEX
    
    IF(L .LE. 1)CALL LSQUIT('ERROR IN Build_PRECALCULATED_SPHMAT',-1)
    NSIZE=0
    NRow = 2*L+1
    NCol = (L+1)*(L+2)/2
    DO M1 = 0, 2*L 
       M = M1 - L
       IF (L.EQ. 1) THEN
          IF (M .EQ. -1) MADR =  0  
          IF (M .EQ.  0) MADR =  1 
          IF (M .EQ.  1) MADR = -1 
       ELSE
          MADR = M
       END IF
       MABS = ABS(M)
       V0 = 0
       IF (M .LT. 0) V0 = 1 
       FACNRM = D1
       IF (M .NE. 0) FACNRM = SQRT(D2*FACULT(6,L+MABS)*&
            &FACULT(6,L-MABS))/(FACULT(6,L)*(D2**MABS))
       FACNRM = FACNRM*DM1**((0-MOD(0,2))/2)*D2**(-0)
       FACNRM = FACNRM/SQRT(FACUL2(6,2*L-1))
       DO T = 0, L - MABS, 2
          DO U = 0, T, 2
             DO V = V0, MABS, 2
                ! almost 6.4.48 in the book
                FAC3 = FACNRM*BINOM(6,L,T/2)*BINOM(6,L-T/2,MABS+T/2)&
                     &                    *BINOM(6,T/2,U/2)*BINOM(6,MABS,V)
                DO A = 0, MIN(0,T+MABS-U-V) 
                   DO B = 0, MIN(0,U+V)
                      DO C = 0, MIN(0,L-T-MABS)
                         ! 6.4.47 in the book
                         DO P = 0, - A, 2
                            DO Q = 0, - B, 2
                               DO R = 0, - C, 2
                                  FACTOR = DM1**(A+B+C+(T+V-V0-P-Q-R)/2)*&
                                       &   D2**(-A-B-C-P-Q-R-T)*FAC3
                                  X = T+MABS-U-V-2*A-P
                                  Y = U+V-2*B-Q
                                  Z = L-T-MABS-2*C-R
                                  TOT = X + Y + Z
                                  IADR = 1 + (2*L+1)*(NCRT(X,Y,Z)-1) + L + MADR+NSIZE
                                  SPHMAT(IADR) = SPHMAT(IADR) + FACTOR 
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
    NSIZE= NSIZE+Nrow*Ncol 
  END SUBROUTINE BUILD_SPHMAT23

  subroutine HorizontalRR_LHS(AUXarrayCont,nTUVQ,nTUVP,nTUV,nContPQ,nPasses,Qdistance12,&
       & AngmomP,AngmomA,AngmomB,THETA,nTUVAspec,nTUVBspec,lupri)
    implicit none
    integer :: nTUV,nTUVP,nTUVQ,nContPQ,nPasses,lupri,AngmomP,AngmomA,AngmomB
    integer :: nTUVAspec,nTUVBspec
    real(realk) :: AUXarrayCont(nTUVQ,nTUV,nContPQ,nPasses),Qdistance12(3,nPasses)
    real(realk) :: THETA(nTUVQ,nTUVAspec,nTUVBspec,nContPQ,nPasses)
    !
    integer,pointer :: TUVindex(:,:,:),Tindex(:),Uindex(:),Vindex(:)
    integer :: NTUVMAX,iTUVplus1y,iTUVplus1x,iTUVplus1z,NTUVAstart
    integer :: J,T,U,V,iTUVA,iTUVB,iTUVAspec
    integer :: iPass,icont,TA,UA,VA,iTUVQ,NTUVAcurrent,currentAngmomA
    real(realk) :: TMP(nTUVQ,nTUVP,nTUVP)
    real(realk) :: TMP2(nTUVQ,10,6),Xab,Yab,Zab
    WRITE(lupri,*)'HorizontalRR_LHS'
    call mem_alloc(TUVindex,AngmomP,AngmomP,AngmomP,.TRUE.,.TRUE.,.TRUE.)
    NTUVMAX = (AngmomP+1)*(AngmomP+2)*(AngmomP+3)/6
    call mem_alloc(Tindex,NTUVMAX)
    call mem_alloc(Uindex,NTUVMAX)
    call mem_alloc(Vindex,NTUVMAX)
    iTUVA=0
    DO J = 0, AngmomP
       DO T = J,0,-1       
          DO U = J-T,0,-1
             iTUVA=iTUVA+1
             V=J-T-U
             TUVindex(T,U,V)=iTUVA
             Tindex(iTUVA)=T
             Uindex(iTUVA)=U
             Vindex(iTUVA)=V
          ENDDO
       ENDDO
    ENDDO
    NTUVAspec = (AngmomA+1)*(AngmomA+2)*(AngmomA+3)/6 - &
         &         (AngmomA)*(AngmomA+1)*(AngmomA+2)/6
    NTUVBspec = (AngmomB+1)*(AngmomB+2)*(AngmomB+3)/6 - &
         &         (AngmomB)*(AngmomB+1)*(AngmomB+2)/6
    Write(lupri,*)'NTUVAspec',NTUVAspec
    Write(lupri,*)'NTUVBspec',NTUVBspec
    NTUVAstart = (AngmomA)*(AngmomA+1)*(AngmomA+2)/6
    Write(lupri,*)'NTUVAstart',NTUVAstart
    
    DO iPass=1,nPasses
       Xab = Qdistance12(1,iPass)
       Yab = Qdistance12(2,iPass)
       Zab = Qdistance12(3,iPass)
       DO icont=1,nContPQ
          !Build Theta(i,j=0)  i=AngmomP => nTUVA = nTUVP, nTUVB=1=(0,0,0)
          !warning 5
          WRITE(lupri,*)'HorizontalRR_LHS  build1'
          DO iTUVA=1,nTUVP
             DO iTUVQ=1,nTUVQ
                TMP(iTUVQ,iTUVA,1) = AUXarrayCont(iTUVQ,iTUVA,iCont,iPass)
                WRITE(lupri,'(A,I3,A,I3,A,ES16.8)')'ATMP(',iTUVQ,',',iTUVA,',1)=',TMP(iTUVQ,iTUVA,1)
             ENDDO
          ENDDO
          !Build Theta(i,j=1)  i=AngmomP-1 => nTUVA < nTUVP, nTUVB=4
          currentAngmomA = AngmomP - 1
          NTUVAcurrent = (currentAngmomA+1)*(currentAngmomA+2)*(currentAngmomA+3)/6          
          !warning 5
          WRITE(lupri,*)'HorizontalRR_LHS  build 2,3,4'
          IF(AngmomB.EQ.1)THEN
             DO iTUVAspec = 1,NTUVAspec
                iTUVA = iTUVAspec+NTUVAstart
                TA = Tindex(iTUVA)
                UA = Uindex(iTUVA)
                VA = Vindex(iTUVA)
                iTUVplus1x = TUVindex(TA+1,UA,VA)
                iTUVplus1y = TUVindex(TA,UA+1,VA)
                iTUVplus1z = TUVindex(TA,UA,VA+1)             
                !(i,100) = (i+1x,000) + X*(i,000)
                !(i,010) = (i+1y,000) + Y*(i,000)
                !(i,001) = (i+1z,000) + Z*(i,000)
                DO iTUVQ=1,nTUVQ
                   TMP(iTUVQ,iTUVAspec,2) = TMP(iTUVQ,iTUVplus1x,1) + Xab*TMP(iTUVQ,iTUVA,1)
                   TMP(iTUVQ,iTUVAspec,3) = TMP(iTUVQ,iTUVplus1y,1) + Yab*TMP(iTUVQ,iTUVA,1)
                   TMP(iTUVQ,iTUVAspec,4) = TMP(iTUVQ,iTUVplus1z,1) + Zab*TMP(iTUVQ,iTUVA,1)
                   WRITE(lupri,'(A,I3,A,I3,A,ES16.8)')'BTMP(',iTUVQ,',',iTUVA,',2)=',TMP(iTUVQ,iTUVA,2)
                   WRITE(lupri,'(A,I3,A,I3,A,ES16.8)')'BTMP(',iTUVQ,',',iTUVA,',3)=',TMP(iTUVQ,iTUVA,3)
                   WRITE(lupri,'(A,I3,A,I3,A,ES16.8)')'BTMP(',iTUVQ,',',iTUVA,',4)=',TMP(iTUVQ,iTUVA,4)
                ENDDO
             ENDDO
          ELSEIF(AngmomB.EQ.2)THEN
             DO iTUVA=1,NTUVAcurrent
                TA = Tindex(iTUVA)
                UA = Uindex(iTUVA)
                VA = Vindex(iTUVA)
                iTUVplus1x = TUVindex(TA+1,UA,VA)
                iTUVplus1y = TUVindex(TA,UA+1,VA)
                iTUVplus1z = TUVindex(TA,UA,VA+1)             
                !(i,100) = (i+1x,000) + X*(i,000)
                !(i,010) = (i+1y,000) + Y*(i,000)
                !(i,001) = (i+1z,000) + Z*(i,000)
                DO iTUVQ=1,nTUVQ
                   TMP(iTUVQ,iTUVA,2) = TMP(iTUVQ,iTUVplus1x,1) + Xab*TMP(iTUVQ,iTUVA,1)
                   TMP(iTUVQ,iTUVA,3) = TMP(iTUVQ,iTUVplus1y,1) + Yab*TMP(iTUVQ,iTUVA,1)
                   TMP(iTUVQ,iTUVA,4) = TMP(iTUVQ,iTUVplus1z,1) + Zab*TMP(iTUVQ,iTUVA,1)
                   WRITE(lupri,'(A,I3,A,I3,A,ES16.8)')'BTMP(',iTUVQ,',',iTUVA,',2)=',TMP(iTUVQ,iTUVA,2)
                   WRITE(lupri,'(A,I3,A,I3,A,ES16.8)')'BTMP(',iTUVQ,',',iTUVA,',3)=',TMP(iTUVQ,iTUVA,3)
                   WRITE(lupri,'(A,I3,A,I3,A,ES16.8)')'BTMP(',iTUVQ,',',iTUVA,',4)=',TMP(iTUVQ,iTUVA,4)
                ENDDO
             ENDDO
             !Build Theta(i,j=1)  i=AngmomP-1 => nTUVA < nTUVP, nTUVB=4
             currentAngmomA = currentAngmomA - 1
             NTUVAcurrent = (currentAngmomA+1)*(currentAngmomA+2)*(currentAngmomA+3)/6          
             !          DO iTUVA=1,NTUVAcurrent
             WRITE(lupri,*)'HorizontalRR_LHS  build A,D'
             DO iTUVAspec = 1,NTUVAspec
                iTUVA = iTUVAspec+NTUVAstart
                TA = Tindex(iTUVA)
                UA = Uindex(iTUVA)
                VA = Vindex(iTUVA)
                iTUVplus1x = TUVindex(TA+1,UA,VA)
                iTUVplus1y = TUVindex(TA,UA+1,VA)
                iTUVplus1z = TUVindex(TA,UA,VA+1)             
                !(i,200) = (i+1x,100) + X*(i,100)
                !(i,110) = (i+1x,010) + X*(i,010)
                !(i,101) = (i+1x,001) + X*(i,001)
                !(i,020) = (i+1y,010) + Y*(i,010)
                !(i,011) = (i+1z,010) + Z*(i,010)
                !(i,002) = (i+1z,001) + Z*(i,001)
                DO iTUVQ=1,nTUVQ
                   TMP2(iTUVQ,iTUVAspec,1) = TMP(iTUVQ,iTUVplus1x,2) + Xab*TMP(iTUVQ,iTUVA,2)
                   TMP2(iTUVQ,iTUVAspec,2) = TMP(iTUVQ,iTUVplus1x,3) + Xab*TMP(iTUVQ,iTUVA,3)
                   TMP2(iTUVQ,iTUVAspec,3) = TMP(iTUVQ,iTUVplus1x,4) + Xab*TMP(iTUVQ,iTUVA,4)
                   TMP2(iTUVQ,iTUVAspec,4) = TMP(iTUVQ,iTUVplus1y,3) + Yab*TMP(iTUVQ,iTUVA,3)
                   TMP2(iTUVQ,iTUVAspec,5) = TMP(iTUVQ,iTUVplus1z,3) + Zab*TMP(iTUVQ,iTUVA,3)
                   TMP2(iTUVQ,iTUVAspec,6) = TMP(iTUVQ,iTUVplus1z,4) + Zab*TMP(iTUVQ,iTUVA,4)
                   WRITE(lupri,'(A,I3,A,I3,A,ES16.8)')'TMP2(',iTUVQ,',',iTUVAspec,',1)=',TMP2(iTUVQ,iTUVAspec,1)
                   WRITE(lupri,'(A,I3,A,I3,A,ES16.8)')'TMP2(',iTUVQ,',',iTUVAspec,',2)=',TMP2(iTUVQ,iTUVAspec,2)
                   WRITE(lupri,'(A,I3,A,I3,A,ES16.8)')'TMP2(',iTUVQ,',',iTUVAspec,',3)=',TMP2(iTUVQ,iTUVAspec,3)
                   WRITE(lupri,'(A,I3,A,I3,A,ES16.8)')'TMP2(',iTUVQ,',',iTUVAspec,',4)=',TMP2(iTUVQ,iTUVAspec,4)
                   WRITE(lupri,'(A,I3,A,I3,A,ES16.8)')'TMP2(',iTUVQ,',',iTUVAspec,',5)=',TMP2(iTUVQ,iTUVAspec,5)
                   WRITE(lupri,'(A,I3,A,I3,A,ES16.8)')'TMP2(',iTUVQ,',',iTUVAspec,',6)=',TMP2(iTUVQ,iTUVAspec,6)
                ENDDO
             ENDDO
          ENDIF


          DO iTUVA = 1,NTUVAspec 
             DO iTUVB = 1,NTUVBspec 
                DO iTUVQ=1,nTUVQ
                   THETA(iTUVQ,iTUVA,iTUVB,icont,iPass) = TMP2(iTUVQ,iTUVA,iTUVB)
                   WRITE(lupri,'(A,I3,A,I3,A,I3,A,ES16.8)')&
                        &'THETA(',iTUVQ,',',iTUVA,',',iTUVB,')',TMP2(iTUVQ,iTUVA,iTUVB)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    call mem_dealloc(TUVindex)
    call mem_dealloc(Tindex)
    call mem_dealloc(Uindex)
    call mem_dealloc(Vindex)
  end subroutine HorizontalRR_LHS

  subroutine PrimitiveContraction(AUXarray2,AUXarrayCont,nTUVfull2,nTUVfull1,nPrimP,nPrimQ,nPasses,&
       & nContP,nContQ,ACC,BCC,CCC,DCC,nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD,&
       & lupri)
    implicit none
    !Warning Primitive screening modifies this!!! 
    integer,intent(in) :: nTUVfull1,nTUVfull2,nPrimP,nPrimQ,nPasses,nContP,nContQ,lupri
    integer,intent(in) :: nPrimA,nContA,nPrimB,nContB,nPrimC,nContC,nPrimD,nContD
    real(realk),intent(in) :: ACC(nPrimA,nContA),BCC(nPrimB,nContB)
    real(realk),intent(in) :: CCC(nPrimC,nContC),DCC(nPrimD,nContD)
    real(realk),intent(in) :: AUXarray2(nTUVfull1,nPrimQ,nPrimP,nPasses)
    real(realk),intent(inout) :: AUXarrayCont(nTUVfull2,nContQ,nContP,nPasses)
    !
    integer :: iPassQ,iContA,iContB,iContC,iContD,iPrimA,iPrimB,iPrimC,iPrimD
    integer :: iTUV,iPrimQ,iPrimP,iContQ,iContP
    real(realk) :: B,ABCDTMP,ABDTMP,ABTMP
    WRITE(lupri,*)'PrimitiveContraction'
    !all passes have same ACCs,BCCs,...
    !maybe construct big CC(nPrimQP,nContQP) matrix and call dgemm nPass times
    !the construction of CC should scale as c**4*p**4 and the 
    !dgemm should scale as c**4*p**4*L**6 but hopefully with efficient FLOP count, although not quadratic matrices....
    !special for nContPQ = 1 
    !special for nContP = 1
    !special for nContQ = 1
    !special for nContA = 1 ...
    !memory should be c**4*p**4 + p**4*L**6 which is fine
    !this would be a simple sum for segmentet! or maybe the sum can be moved into the previous electron transfer reccurence
    do iPassQ = 1,nPasses
     do iContB=1,nContB
      do iContA=1,nContA
       iContP = iContA+(iContB-1)*nContA
       do iContD=1,nContD
        do iContC=1,nContC
         iContQ = iContC+(iContD-1)*nContC
         do iTUV=1,nTUVfull2
          AUXarrayCont(iTUV,iContQ,iContP,iPassQ) = 0.0E0_realk
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
    do iPassQ = 1,nPasses
     do iContB=1,nContB
      do iPrimB=1,nPrimB
       B = BCC(iPrimB,iContB)
       do iContA=1,nContA
        iContP = iContA+(iContB-1)*nContA
        do iPrimA=1,nPrimA
         iPrimP = iPrimA + (iPrimB-1)*nPrimA
         ABTMP = ACC(iPrimA,iContA)*B
         do iContD=1,nContD
          do iPrimD=1,nPrimD
           ABDTMP = DCC(iPrimD,iContD)*ABTMP
           do iContC=1,nContC
            iContQ = iContC+(iContD-1)*nContC
            iPrimQ = (iPrimD-1)*nPrimC
            do iPrimC=1,nPrimC
             ABCDTMP = CCC(iPrimC,iContC)*ABDTMP
             print*,'ABCDTMP',ABCDTMP
             iPrimQ = iPrimQ + 1
             do iTUV=1,nTUVfull2
              AUXarrayCont(iTUV,iContQ,iContP,iPassQ) = AUXarrayCont(iTUV,iContQ,iContP,iPassQ) + &
                   & ABCDTMP*AUXarray2(iTUV,iPrimQ,iPrimP,iPassQ)
              WRITE(lupri,'(A,I6,A,ES16.8)')'cont AUXarrayCont(',iTUV,',iContQ,iContP,iPassQ)',&
                   & ABCDTMP*AUXarray2(iTUV,iPrimQ,iPrimP,iPassQ)
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
    do iPassQ = 1,nPasses
     do iContB=1,nContB
      do iContA=1,nContA
       iContP = iContA+(iContB-1)*nContA
       do iContD=1,nContD
        do iContC=1,nContC
         iContQ = iContC+(iContD-1)*nContC
         do iTUV=1,nTUVfull2
            WRITE(lupri,'(A,I6,A,ES16.8)')'AUXarrayCont(',iTUV,',iContQ,iContP,iPassQ)',&
                 & AUXarrayCont(iTUV,iContQ,iContP,iPassQ)
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo

  end subroutine PrimitiveContraction

  subroutine TransferRecurrence1(nPasses,nPrimP,nPrimQ,nTUVQ,nTUV,LimitAngmom,currentMaxAngmom,&
       & reducedExponents,Pexp,Qexp,Rpq,Rqc,Rpa,AUXarray2,currentMaxAngmomP,lupri)
    implicit none
    integer,intent(in) :: nPasses,nPrimP,nPrimQ,currentMaxAngmom,nTUVQ,nTUV,LimitAngmom,lupri
    integer,intent(in) :: currentMaxAngmomP
    real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Rpq(3,nPrimQ,nPrimP,nPasses)
    real(realk),intent(in) :: Rqc(3,nPrimQ),Rpa(3,nPrimP),Qexp(nPrimQ),Pexp(nPrimP)
    real(realk),intent(inout) :: AUXarray2(nTUVQ,nTUV,nPasses*nPrimP*nPrimQ)
    !
    integer :: iPassQ,iPrimP,iPrimQ,J,I,N,offset,T,U,V,iTUV,ipq,nTUV2
    real(realk) :: expP,invexpQ,alphaQ,alphaXpq,alphaYpq,alphaZpq,Xqc,Yqc,Zqc
    real(realk) :: Xpa,Ypa,Zpa,inv2expQ,pinvq,facX,facY,facZ
    real(realk),parameter :: D1=1.0E0_realk,D2=2.0E0_realk
    integer,pointer :: TUVindex(:,:,:),Tindex(:),Uindex(:),Vindex(:),Iindex(:)
    integer :: NTUVMAX,iTUVplus1y,iTUVplus1x,iTUVplus1z,iTUVminus1x,iTUVminus1y,iTUVminus1z
    WRITE(lupri,*)'TransferRecurrence1'
    
    call mem_alloc(TUVindex,LimitAngmom,LimitAngmom,LimitAngmom,.TRUE.,.TRUE.,.TRUE.)
    NTUVMAX = (LimitAngmom+1)*(LimitAngmom+2)*(LimitAngmom+3)/6
    nTUV2 = (currentMaxAngmomP+1)*(currentMaxAngmomP+2)*(currentMaxAngmomP+3)/6
    write(lupri,*)'currentMaxAngmom',currentMaxAngmom
    write(lupri,*)'nTUV2',nTUV2,'currentMaxAngmomP',currentMaxAngmomP
    write(lupri,*)'nTUV',nTUV,'LimitAngmom',LimitAngmom
    call mem_alloc(Tindex,NTUVMAX)
    call mem_alloc(Uindex,NTUVMAX)
    call mem_alloc(Vindex,NTUVMAX)
    call mem_alloc(Iindex,NTUVMAX)
    iTUV=0
    DO J = 0, LimitAngmom
       DO T = J,0,-1       
          DO U = J-T,0,-1
             iTUV=iTUV+1
             V=J-T-U
             TUVindex(T,U,V)=iTUV
             Iindex(iTUV)=J
             Tindex(iTUV)=T
             Uindex(iTUV)=U
             Vindex(iTUV)=V
          ENDDO
       ENDDO
    ENDDO

    WRITE(lupri,*)'TransferRecurrence1  currentMaxAngmom=',currentMaxAngmom
    !k=0  means special case 
    !ThetaAux2(1,i) = (Xqc+(p/q)*Xpa)*ThetaAux2(0,i) + i/(2*q)*ThetaAux2(0,i-1) 
    !                    - (p/q)*ThetaAux2(0,i+1)
    do iPassQ = 1,nPasses
     offset = (iPassQ-1)*nPrimQ*nPrimP
     do IPrimP = 1,nPrimP         
      expP = Pexp(iPrimP)
      Xpa = Rpa(1,iPrimP)
      Ypa = Rpa(2,iPrimP)
      Zpa = Rpa(3,iPrimP)
      IPQ = (iPrimP-1)*nPrimQ + offset
      !for segmentet you can sum all iprimQ into TMP(2),TMP(3),TMP(4) and plug into AUXarray2(2,iTUV,Ipass)
      do IPrimQ = 1,nPrimQ 
       invexpQ = D1/Qexp(iPrimQ)
       inv2expQ = D1/(D2*Qexp(iPrimQ))
       Xqc = Rqc(1,iPrimQ)
       Yqc = Rqc(2,iPrimQ)
       Zqc = Rqc(3,iPrimQ)
       alphaQ = -reducedExponents(iPrimQ,iPrimP)*invexpQ
       alphaXpq = alphaQ*Rpq(1,iPrimQ,iPrimP,iPassQ)
       alphaYpq = alphaq*Rpq(2,iPrimQ,iPrimP,iPassQ)
       alphaZpq = alphaQ*Rpq(3,iPrimQ,iPrimP,iPassQ)
       pinvq = expP*invexpQ
       IPQ = IPQ + 1
       facX = (Xqc+pinvq*Xpa)
       facY = (Yqc+pinvq*Ypa)
       facZ = (Zqc+pinvq*Zpa)
       
       do ITUV=1,nTUV2
          I = Iindex(iTUV)
          T = Tindex(iTUV) 
          U = Uindex(iTUV) 
          V = Vindex(iTUV)
          iTUVplus1x = TUVindex(T+1,U,V)
          iTUVplus1y = TUVindex(T,U+1,V)
          iTUVplus1z = TUVindex(T,U,V+1)
          
          IF(I.EQ.0)THEN
             AUXarray2(2,iTUV,IPQ) = facX*AUXarray2(1,iTUV,IPQ) - pinvq*AUXarray2(1,iTUVplus1x,IPQ)
             AUXarray2(3,iTUV,IPQ) = facY*AUXarray2(1,iTUV,IPQ) - pinvq*AUXarray2(1,iTUVplus1y,IPQ)
             AUXarray2(4,iTUV,IPQ) = facZ*AUXarray2(1,iTUV,IPQ) - pinvq*AUXarray2(1,iTUVplus1z,IPQ)
          ELSE
             IF(T.EQ.0)THEN
                AUXarray2(2,iTUV,IPQ) = facX*AUXarray2(1,iTUV,IPQ) - pinvq*AUXarray2(1,iTUVplus1x,IPQ)
             ELSE
                iTUVminus1x = TUVindex(T-1,U,V)
                AUXarray2(2,iTUV,IPQ) = facX*AUXarray2(1,iTUV,IPQ) + I*inv2expQ*AUXarray2(1,iTUVminus1x,IPQ) - pinvq*AUXarray2(1,iTUVplus1x,IPQ)
             ENDIF
             IF(U.EQ.0)THEN
                AUXarray2(3,iTUV,IPQ) = facY*AUXarray2(1,iTUV,IPQ) - pinvq*AUXarray2(1,iTUVplus1y,IPQ)
             ELSE
                iTUVminus1y = TUVindex(T,U-1,V)
                AUXarray2(3,iTUV,IPQ) = facY*AUXarray2(1,iTUV,IPQ) + I*inv2expQ*AUXarray2(1,iTUVminus1y,IPQ) - pinvq*AUXarray2(1,iTUVplus1y,IPQ)
             ENDIF
             IF(V.EQ.0)THEN
                AUXarray2(4,iTUV,IPQ) = facZ*AUXarray2(1,iTUV,IPQ) - pinvq*AUXarray2(1,iTUVplus1z,IPQ)             
             ELSE
                iTUVminus1z = TUVindex(T,U,V-1)
                AUXarray2(4,iTUV,IPQ) = facZ*AUXarray2(1,iTUV,IPQ) + I*inv2expQ*AUXarray2(1,iTUVminus1z,IPQ) - pinvq*AUXarray2(1,iTUVplus1z,IPQ)
             ENDIF
          ENDIF
          WRITE(lupri,'(A,I3,A,ES16.8)')'AUXarray2(1,',iTUV,',IPQ) ',AUXarray2(1,iTUV,IPQ)
          WRITE(lupri,'(A,I3,A,ES16.8)')'AUXarray2(2,',iTUV,',IPQ) ',AUXarray2(2,iTUV,IPQ)
          WRITE(lupri,'(A,I3,A,ES16.8)')'AUXarray2(3,',iTUV,',IPQ) ',AUXarray2(3,iTUV,IPQ)
          WRITE(lupri,'(A,I3,A,ES16.8)')'AUXarray2(4,',iTUV,',IPQ) ',AUXarray2(4,iTUV,IPQ)
       enddo
      enddo
     enddo
    enddo
    call mem_alloc(TUVindex,LimitAngmom,LimitAngmom,LimitAngmom)
    call mem_alloc(Tindex,LimitAngmom)
    call mem_alloc(Uindex,LimitAngmom)
    call mem_alloc(Vindex,LimitAngmom)
  end subroutine TransferRecurrence1

   subroutine GenerateThetaFromAuxarray(nTUVQ,nTUV,nPrimPassPQ,AUXarray2,Auxarray,LimitAngmom,lupri)
      implicit none
      integer,intent(in) :: nTUV,nPrimPassPQ,nTUVQ,LimitAngmom,lupri
      real(realk),intent(inout) :: AUXarray2(nTUVQ,nTUV,nPrimPassPQ)
      real(realk),intent(in) :: AUXarray(nTUV,LimitAngmom+1,nPrimPassPQ)
      !
      integer :: IPQ,iTUV
      WRITE(lupri,*)'GenerateThetaFromAuxarray'
      do IPQ=1,nPrimPassPQ
         do iTUV=1,nTUV
            AUXarray2(1,iTUV,IPQ) = AUXarray(iTUV,1,IPQ)
         enddo
      enddo
    end subroutine GenerateThetaFromAuxarray


  subroutine VerticalRecurrence1(nPasses,nPrimP,nPrimQ,nTUV,LimitAngmom,currentMaxAngmom,&
       & reducedExponents,Pexp,Rpq,Rpa,AUXarray,lupri)
    implicit none
    integer,intent(in) :: nPasses,nPrimP,nPrimQ,currentMaxAngmom,nTUV,LimitAngmom,lupri
    real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Rpq(3,nPrimQ,nPrimP,nPasses)
    real(realk),intent(in) :: Rpa(3,nPrimP)

    real(realk),intent(inout) :: AUXarray(nTUV,LimitAngmom+1,nPasses*nPrimP*nPrimQ)
    !
    integer :: iPassQ,iPrimP,iPrimQ,J,I
    real(realk) :: invexpP,alphaP,alphaXpq,alphaYpq,alphaZpq,Xpa,Ypa,Zpa
    real(realk),parameter :: D1=1.0E0_realk
    WRITE(lupri,*)'VerticalRecurrence1  currentMaxAngmom=',currentMaxAngmom
    !ThetaAux(n,1,0,0) = Xpa*ThetaAux(n,0,0,0) + (-alpha/p*Xpq)*ThetaAux(n+1,0,0,0) 
    !i = 0 last 2 term vanish
    do IPrimP = 1,nPrimP
       invexpP = D1/Pexp(iPrimP)
       Xpa = Rpa(1,iPrimP)
       Ypa = Rpa(2,iPrimP)
       Zpa = Rpa(3,iPrimP)
       do iPassQ = 1,nPasses
          I = (iPrimP-1)*nPrimQ + (iPassQ-1)*nPrimQ*nPrimP
          do IPrimQ = 1,nPrimQ
             alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
             alphaXpq = alphaP*Rpq(1,iPrimQ,iPrimP,iPassQ)
             alphaYpq = alphaP*Rpq(2,iPrimQ,iPrimP,iPassQ)
             alphaZpq = alphaP*Rpq(3,iPrimQ,iPrimP,iPassQ)
             write(lupri,*)'alphaXpq',alphaXpq
             write(lupri,*)'alphaXpq',alphaYpq
             write(lupri,*)'alphaXpq',alphaZpq
             I = I + 1
             do J=1,currentMaxAngmom
                AUXarray(2,J,I) = Xpa*AUXarray(1,J,I) + alphaXpq*AUXarray(1,J+1,I)
                AUXarray(3,J,I) = Ypa*AUXarray(1,J,I) + alphaYpq*AUXarray(1,J+1,I)
                AUXarray(4,J,I) = Zpa*AUXarray(1,J,I) + alphaZpq*AUXarray(1,J+1,I)
                write(lupri,*)'AUXarray(2,J,I)',AUXarray(2,J,I)
                write(lupri,*)'AUXarray(3,J,I)',AUXarray(3,J,I)
                write(lupri,*)'AUXarray(4,J,I)',AUXarray(4,J,I)
             enddo
          enddo
       enddo
    enddo
  end subroutine VerticalRecurrence1

  subroutine VerticalRecurrence2(nPasses,nPrimP,nPrimQ,nTUV,LimitAngmom,currentMaxAngmom,&
       & reducedExponents,Pexp,Rpq,Rpa,AUXarray,lupri)
    implicit none
    integer,intent(in) :: nPasses,nPrimP,nPrimQ,currentMaxAngmom,nTUV,LimitAngmom,lupri
    real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Rpq(3,nPrimQ,nPrimP,nPasses)
    real(realk),intent(in) :: Rpa(3,nPrimP)

    real(realk),intent(inout) :: AUXarray(nTUV,LimitAngmom+1,nPasses*nPrimP*nPrimQ)
    !
    integer :: iPassQ,iPrimP,iPrimQ,J,I
    real(realk) :: invexpP,alphaP,alphaXpq,alphaYpq,alphaZpq,Xpa,Ypa,Zpa
    real(realk) :: inv2expP,TwoTerms
    real(realk),parameter :: D1=1.0E0_realk,D2=2.0E0_realk
    WRITE(lupri,*)'VerticalRecurrence2  currentMaxAngmom=',currentMaxAngmom
    !ThetaAux(n,i+1,0,0) = Xpa*ThetaAux(n,i,0,0) + (-alpha/p*Xpq)*ThetaAux(n+1,i,0,0) 
    !                    + i/(2p) * (ThetaAux(n,i-1,0,0) - alpha/p*ThetaAux(n+1,i-1,0,0))
    !i = 1 
    do IPrimP = 1,nPrimP
       invexpP = D1/Pexp(iPrimP)
       inv2expP = D1/(D2*Pexp(iPrimP))   
       Xpa = Rpa(1,iPrimP)
       Ypa = Rpa(2,iPrimP)
       Zpa = Rpa(3,iPrimP)
       do iPassQ = 1,nPasses
          I = (iPrimP-1)*nPrimQ + (iPassQ-1)*nPrimQ*nPrimP
          do IPrimQ = 1,nPrimQ
             alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
             alphaXpq = alphaP*Rpq(1,iPrimQ,iPrimP,iPassQ)
             alphaYpq = alphaP*Rpq(2,iPrimQ,iPrimP,iPassQ)
             alphaZpq = alphaP*Rpq(3,iPrimQ,iPrimP,iPassQ)
             I = I + 1
             do J=1,currentMaxAngmom
                TwoTerms = inv2expP*(AUXarray(1,J,I) + alphaP*AUXarray(1,J+1,I))

                AUXarray(5,J,I) = Xpa*AUXarray(2,J,I) + alphaXpq*AUXarray(2,J+1,I) + TwoTerms
                AUXarray(6,J,I) = Xpa*AUXarray(3,J,I) + alphaXpq*AUXarray(3,J+1,I)
                AUXarray(7,J,I) = Xpa*AUXarray(4,J,I) + alphaXpq*AUXarray(4,J+1,I)
                AUXarray(8,J,I) = Ypa*AUXarray(3,J,I) + alphaYpq*AUXarray(3,J+1,I) + TwoTerms
                AUXarray(9,J,I) = Zpa*AUXarray(3,J,I) + alphaZpq*AUXarray(3,J+1,I)
                AUXarray(10,J,I)= Zpa*AUXarray(4,J,I) + alphaZpq*AUXarray(4,J+1,I) + TwoTerms
                write(lupri,*)'AUXarray(5,J,I) ',AUXarray(5,J,I)
                write(lupri,*)'AUXarray(6,J,I) ',AUXarray(6,J,I)
                write(lupri,*)'AUXarray(7,J,I) ',AUXarray(7,J,I)
                write(lupri,*)'AUXarray(8,J,I) ',AUXarray(8,J,I)
                write(lupri,*)'AUXarray(9,J,I) ',AUXarray(9,J,I)
                write(lupri,*)'AUXarray(10,J,I)',AUXarray(10,J,I)
             enddo
          enddo
       enddo
    enddo
  end subroutine VerticalRecurrence2

  subroutine VerticalRecurrence3(nPasses,nPrimP,nPrimQ,nTUV,LimitAngmom,currentMaxAngmom,&
       & reducedExponents,Pexp,Rpq,Rpa,AUXarray,lupri)
    implicit none
    integer,intent(in) :: nPasses,nPrimP,nPrimQ,currentMaxAngmom,nTUV,LimitAngmom,lupri
    real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Rpq(3,nPrimQ,nPrimP,nPasses)
    real(realk),intent(in) :: Rpa(3,nPrimP)

    real(realk),intent(inout) :: AUXarray(nTUV,LimitAngmom+1,nPasses*nPrimP*nPrimQ)
    !
    integer :: iPassQ,iPrimP,iPrimQ,J,I
    real(realk) :: invexpP,alphaP,alphaXpq,alphaYpq,alphaZpq,Xpa,Ypa,Zpa
    real(realk) :: inv2expP,TwoTerms1,TwoTerms2,TwoTerms3
    real(realk),parameter :: D1=1.0E0_realk,D2=2.0E0_realk
    WRITE(lupri,*)'VerticalRecurrence3  currentMaxAngmom=',currentMaxAngmom
    !ThetaAux(n,i+1,0,0) = Xpa*ThetaAux(n,i,0,0) + (-alpha/p*Xpq)*ThetaAux(n+1,i,0,0) 
    !                    + i/(2p) * (ThetaAux(n,i-1,0,0) - alpha/p*ThetaAux(n+1,i-1,0,0))
    !i = 2
    !i/(2p) = 1/p
    do IPrimP = 1,nPrimP
       invexpP = D1/Pexp(iPrimP)
!       inv2expP = D1/(D2*Pexp(iPrimP))   
       Xpa = Rpa(1,iPrimP)
       Ypa = Rpa(2,iPrimP)
       Zpa = Rpa(3,iPrimP)
       do iPassQ = 1,nPasses
          I = (iPrimP-1)*nPrimQ + (iPassQ-1)*nPrimQ*nPrimP
          do IPrimQ = 1,nPrimQ
             alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
             alphaXpq = alphaP*Rpq(1,iPrimQ,iPrimP,iPassQ)
             alphaYpq = alphaP*Rpq(2,iPrimQ,iPrimP,iPassQ)
             alphaZpq = alphaP*Rpq(3,iPrimQ,iPrimP,iPassQ)
             I = I + 1
             do J=1,currentMaxAngmom
!                TwoTerms1 = inv2expP*(AUXarray(2,J,I) + alphaP*AUXarray(2,J+1,I))
!                TwoTerms2 = inv2expP*(AUXarray(3,J,I) + alphaP*AUXarray(3,J+1,I))
!                TwoTerms3 = inv2expP*(AUXarray(4,J,I) + alphaP*AUXarray(4,J+1,I))
                TwoTerms1 = invexpP*(AUXarray(2,J,I) + alphaP*AUXarray(2,J+1,I))
                TwoTerms2 = invexpP*(AUXarray(3,J,I) + alphaP*AUXarray(3,J+1,I))
                TwoTerms3 = invexpP*(AUXarray(4,J,I) + alphaP*AUXarray(4,J+1,I))

                AUXarray(11,J,I) = Xpa*AUXarray(5,J,I) + alphaXpq*AUXarray(5,J+1,I) + TwoTerms1
                AUXarray(12,J,I) = Ypa*AUXarray(5,J,I) + alphaYpq*AUXarray(5,J+1,I)
                AUXarray(13,J,I) = Zpa*AUXarray(5,J,I) + alphaZpq*AUXarray(5,J+1,I)
                AUXarray(14,J,I) = Xpa*AUXarray(8,J,I) + alphaXpq*AUXarray(8,J+1,I)
                AUXarray(15,J,I) = Xpa*AUXarray(9,J,I) + alphaXpq*AUXarray(9,J+1,I)
                AUXarray(16,J,I) = Xpa*AUXarray(10,J,I)+ alphaXpq*AUXarray(10,J+1,I)
                AUXarray(17,J,I) = Ypa*AUXarray(8,J,I) + alphaYpq*AUXarray(8,J+1,I) + TwoTerms2
                AUXarray(18,J,I) = Zpa*AUXarray(8,J,I) + alphaZpq*AUXarray(8,J+1,I)
                AUXarray(19,J,I) = Ypa*AUXarray(10,J,I)+ alphaYpq*AUXarray(10,J+1,I)
                AUXarray(20,J,I) = Zpa*AUXarray(10,J,I)+ alphaZpq*AUXarray(10,J+1,I) + TwoTerms3
                write(lupri,*)'AUXarray(11,J,I)',AUXarray(11,J,I)
                write(lupri,*)'AUXarray(12,J,I)',AUXarray(12,J,I)
                write(lupri,*)'AUXarray(13,J,I)',AUXarray(13,J,I)
                write(lupri,*)'AUXarray(14,J,I)',AUXarray(14,J,I)
                write(lupri,*)'AUXarray(15,J,I)',AUXarray(15,J,I)
                write(lupri,*)'AUXarray(16,J,I)',AUXarray(16,J,I)
                write(lupri,*)'AUXarray(17,J,I)',AUXarray(17,J,I)
                write(lupri,*)'AUXarray(18,J,I)',AUXarray(18,J,I)
                write(lupri,*)'AUXarray(19,J,I)',AUXarray(19,J,I)
                write(lupri,*)'AUXarray(20,J,I)',AUXarray(20,J,I)
             enddo
          enddo
       enddo
    enddo
  end subroutine VerticalRecurrence3

  subroutine VerticalRecurrence4(nPasses,nPrimP,nPrimQ,nTUV,LimitAngmom,currentMaxAngmom,&
       & reducedExponents,Pexp,Rpq,Rpa,AUXarray,lupri)
    implicit none
    integer,intent(in) :: nPasses,nPrimP,nPrimQ,currentMaxAngmom,nTUV,LimitAngmom,lupri
    real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Rpq(3,nPrimQ,nPrimP,nPasses)
    real(realk),intent(in) :: Rpa(3,nPrimP)

    real(realk),intent(inout) :: AUXarray(nTUV,LimitAngmom+1,nPasses*nPrimP*nPrimQ)
    !
    integer :: iPassQ,iPrimP,iPrimQ,J,I
    real(realk) :: invexpP,alphaP,alphaXpq,alphaYpq,alphaZpq,Xpa,Ypa,Zpa
    real(realk) :: inv2expP,TwoTerms1,TwoTerms2,TwoTerms3
    real(realk),parameter :: D1=1.0E0_realk,D2=2.0E0_realk, D3=3.0E0_realk
    WRITE(lupri,*)'VerticalRecurrence4  currentMaxAngmom=',currentMaxAngmom
    !ThetaAux(n,i+1,0,0) = Xpa*ThetaAux(n,i,0,0) + (-alpha/p*Xpq)*ThetaAux(n+1,i,0,0) 
    !                    + i/(2p) * (ThetaAux(n,i-1,0,0) - alpha/p*ThetaAux(n+1,i-1,0,0))
    !i = 3
    do IPrimP = 1,nPrimP
       invexpP = D1/Pexp(iPrimP)
       inv2expP = D1/(D2*Pexp(iPrimP))   
       Xpa = Rpa(1,iPrimP)
       Ypa = Rpa(2,iPrimP)
       Zpa = Rpa(3,iPrimP)
       do iPassQ = 1,nPasses
          I = (iPrimP-1)*nPrimQ + (iPassQ-1)*nPrimQ*nPrimP
          do IPrimQ = 1,nPrimQ
             alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
             alphaXpq = alphaP*Rpq(1,iPrimQ,iPrimP,iPassQ)
             alphaYpq = alphaP*Rpq(2,iPrimQ,iPrimP,iPassQ)
             alphaZpq = alphaP*Rpq(3,iPrimQ,iPrimP,iPassQ)
             I = I + 1
             do J=1,currentMaxAngmom
                TwoTerms1 = inv2expP*(AUXarray(5,J,I) + alphaP*AUXarray(5,J+1,I))
                TwoTerms2 = inv2expP*(AUXarray(8,J,I) + alphaP*AUXarray(8,J+1,I))
                TwoTerms3 = inv2expP*(AUXarray(10,J,I) + alphaP*AUXarray(10,J+1,I))

                AUXarray(21,J,I) = Xpa*AUXarray(11,J,I) + alphaXpq*AUXarray(11,J+1,I) + D3*TwoTerms1
                AUXarray(22,J,I) = Ypa*AUXarray(11,J,I) + alphaYpq*AUXarray(11,J+1,I)
                AUXarray(23,J,I) = Zpa*AUXarray(11,J,I) + alphaZpq*AUXarray(11,J+1,I)
                AUXarray(24,J,I) = Ypa*AUXarray(12,J,I) + alphaYpq*AUXarray(12,J+1,I) + TwoTerms1
                AUXarray(25,J,I) = Zpa*AUXarray(12,J,I) + alphaZpq*AUXarray(12,J+1,I)
                AUXarray(26,J,I) = Zpa*AUXarray(13,J,I) + alphaZpq*AUXarray(13,J+1,I) + TwoTerms1
                AUXarray(27,J,I) = Xpa*AUXarray(17,J,I) + alphaXpq*AUXarray(17,J+1,I)
                AUXarray(28,J,I) = Xpa*AUXarray(18,J,I) + alphaXpq*AUXarray(18,J+1,I)
                AUXarray(29,J,I) = Xpa*AUXarray(19,J,I) + alphaXpq*AUXarray(19,J+1,I)
                AUXarray(30,J,I) = Xpa*AUXarray(20,J,I) + alphaXpq*AUXarray(20,J+1,I)
                AUXarray(31,J,I) = Ypa*AUXarray(17,J,I) + alphaYpq*AUXarray(17,J+1,I) + D3*TwoTerms2
                AUXarray(32,J,I) = Zpa*AUXarray(17,J,I) + alphaZpq*AUXarray(17,J+1,I)
                AUXarray(33,J,I) = Zpa*AUXarray(18,J,I) + alphaZpq*AUXarray(18,J+1,I) + TwoTerms2
                AUXarray(34,J,I) = Ypa*AUXarray(20,J,I) + alphaYpq*AUXarray(20,J+1,I)
                AUXarray(35,J,I) = Zpa*AUXarray(20,J,I) + alphaZpq*AUXarray(20,J+1,I) + D3*TwoTerms3
                write(lupri,*)'AUXarray(21,J,I)',AUXarray(21,J,I)
                write(lupri,*)'AUXarray(22,J,I)',AUXarray(22,J,I)
                write(lupri,*)'AUXarray(23,J,I)',AUXarray(23,J,I)
                write(lupri,*)'AUXarray(24,J,I)',AUXarray(24,J,I)
                write(lupri,*)'AUXarray(25,J,I)',AUXarray(25,J,I)
                write(lupri,*)'AUXarray(26,J,I)',AUXarray(26,J,I)
                write(lupri,*)'AUXarray(27,J,I)',AUXarray(27,J,I)
                write(lupri,*)'AUXarray(28,J,I)',AUXarray(28,J,I)
                write(lupri,*)'AUXarray(29,J,I)',AUXarray(29,J,I)
                write(lupri,*)'AUXarray(30,J,I)',AUXarray(30,J,I)
                write(lupri,*)'AUXarray(31,J,I)',AUXarray(31,J,I)
                write(lupri,*)'AUXarray(32,J,I)',AUXarray(32,J,I)
                write(lupri,*)'AUXarray(33,J,I)',AUXarray(33,J,I)
                write(lupri,*)'AUXarray(34,J,I)',AUXarray(34,J,I)
                write(lupri,*)'AUXarray(35,J,I)',AUXarray(35,J,I)
             enddo
          enddo
       enddo
    enddo
  end subroutine VerticalRecurrence4

  subroutine VerticalRecurrence5(nPasses,nPrimP,nPrimQ,nTUV,LimitAngmom,currentMaxAngmom,&
       & reducedExponents,Pexp,Rpq,Rpa,AUXarray,lupri)
    implicit none
    integer,intent(in) :: nPasses,nPrimP,nPrimQ,currentMaxAngmom,nTUV,LimitAngmom,lupri
    real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Rpq(3,nPrimQ,nPrimP,nPasses)
    real(realk),intent(in) :: Rpa(3,nPrimP)

    real(realk),intent(inout) :: AUXarray(nTUV,LimitAngmom+1,nPasses*nPrimP*nPrimQ)
    !
    integer :: iPassQ,iPrimP,iPrimQ,J,I,ipr
    real(realk) :: invexpP,alphaP,alphaXpq,alphaYpq,alphaZpq,Xpa,Ypa,Zpa
    real(realk) :: inv2expP,TwoTerms1,TwoTerms2,TwoTerms3
    real(realk),parameter :: D1=1.0E0_realk,D2=2.0E0_realk, D4=4.0E0_realk
    WRITE(lupri,*)'VerticalRecurrence5  currentMaxAngmom=',currentMaxAngmom
    !ThetaAux(n,i+1,0,0) = Xpa*ThetaAux(n,i,0,0) + (-alpha/p*Xpq)*ThetaAux(n+1,i,0,0) 
    !                    + i/(2p) * (ThetaAux(n,i-1,0,0) - alpha/p*ThetaAux(n+1,i-1,0,0))
    !i = 4
    do IPrimP = 1,nPrimP
       invexpP = D1/Pexp(iPrimP)
       inv2expP = D1/(D2*Pexp(iPrimP))   
       Xpa = Rpa(1,iPrimP)
       Ypa = Rpa(2,iPrimP)
       Zpa = Rpa(3,iPrimP)
       do iPassQ = 1,nPasses
          I = (iPrimP-1)*nPrimQ + (iPassQ-1)*nPrimQ*nPrimP
          do IPrimQ = 1,nPrimQ
             alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
             alphaXpq = alphaP*Rpq(1,iPrimQ,iPrimP,iPassQ)
             alphaYpq = alphaP*Rpq(2,iPrimQ,iPrimP,iPassQ)
             alphaZpq = alphaP*Rpq(3,iPrimQ,iPrimP,iPassQ)
             I = I + 1
             do J=1,currentMaxAngmom
                TwoTerms1 = inv2expP*(AUXarray(11,J,I) + alphaP*AUXarray(11,J+1,I))
                TwoTerms2 = inv2expP*(AUXarray(17,J,I) + alphaP*AUXarray(17,J+1,I))
                TwoTerms3 = inv2expP*(AUXarray(20,J,I) + alphaP*AUXarray(20,J+1,I))

                AUXarray(36,J,I) = Xpa*AUXarray(21,J,I) + alphaXpq*AUXarray(21,J+1,I) + D4*TwoTerms1
                AUXarray(37,J,I) = Ypa*AUXarray(21,J,I) + alphaYpq*AUXarray(21,J+1,I)
                AUXarray(38,J,I) = Zpa*AUXarray(21,J,I) + alphaZpq*AUXarray(21,J+1,I)
                AUXarray(39,J,I) = Ypa*AUXarray(22,J,I) + alphaYpq*AUXarray(22,J+1,I) + TwoTerms1
                AUXarray(40,J,I) = Zpa*AUXarray(22,J,I) + alphaZpq*AUXarray(22,J+1,I)
                AUXarray(41,J,I) = Zpa*AUXarray(23,J,I) + alphaZpq*AUXarray(23,J+1,I) + TwoTerms1
                AUXarray(42,J,I) = Xpa*AUXarray(27,J,I) + alphaXpq*AUXarray(27,J+1,I) + TwoTerms2
                AUXarray(43,J,I) = Zpa*AUXarray(24,J,I) + alphaZpq*AUXarray(24,J+1,I)
                AUXarray(44,J,I) = Ypa*AUXarray(26,J,I) + alphaYpq*AUXarray(26,J+1,I)
                AUXarray(45,J,I) = Xpa*AUXarray(30,J,I) + alphaXpq*AUXarray(30,J+1,I) + TwoTerms3
                AUXarray(46,J,I) = Xpa*AUXarray(31,J,I) + alphaXpq*AUXarray(31,J+1,I)
                AUXarray(47,J,I) = Xpa*AUXarray(32,J,I) + alphaXpq*AUXarray(32,J+1,I)
                AUXarray(48,J,I) = Xpa*AUXarray(33,J,I) + alphaXpq*AUXarray(33,J+1,I)
                AUXarray(49,J,I) = Xpa*AUXarray(34,J,I) + alphaXpq*AUXarray(34,J+1,I)
                AUXarray(50,J,I) = Xpa*AUXarray(35,J,I) + alphaXpq*AUXarray(35,J+1,I)
                AUXarray(51,J,I) = Ypa*AUXarray(31,J,I) + alphaYpq*AUXarray(31,J+1,I) + D4*TwoTerms2
                AUXarray(52,J,I) = Zpa*AUXarray(31,J,I) + alphaZpq*AUXarray(31,J+1,I)
                AUXarray(53,J,I) = Zpa*AUXarray(32,J,I) + alphaZpq*AUXarray(32,J+1,I) + TwoTerms2
                AUXarray(54,J,I) = Ypa*AUXarray(34,J,I) + alphaYpq*AUXarray(34,J+1,I) + TwoTerms3
                AUXarray(55,J,I) = Ypa*AUXarray(35,J,I) + alphaYpq*AUXarray(35,J+1,I)
                AUXarray(56,J,I) = Zpa*AUXarray(35,J,I) + alphaZpq*AUXarray(35,J+1,I) + D4*TwoTerms3
                do ipr=36,56
                   write(lupri,*)'AUXarray(',ipr,',J,I)',AUXarray(ipr,J,I)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine VerticalRecurrence5

  subroutine VerticalRecurrence6(nPasses,nPrimP,nPrimQ,nTUV,LimitAngmom,currentMaxAngmom,&
       & reducedExponents,Pexp,Rpq,Rpa,AUXarray,lupri)
    implicit none
    integer,intent(in) :: nPasses,nPrimP,nPrimQ,currentMaxAngmom,nTUV,LimitAngmom,lupri
    real(realk),intent(in) :: reducedExponents(nPrimQ,nPrimP),Pexp(nPrimP),Rpq(3,nPrimQ,nPrimP,nPasses)
    real(realk),intent(in) :: Rpa(3,nPrimP)

    real(realk),intent(inout) :: AUXarray(nTUV,LimitAngmom+1,nPasses*nPrimP*nPrimQ)
    !
    integer :: iPassQ,iPrimP,iPrimQ,J,I,ipr
    real(realk) :: invexpP,alphaP,alphaXpq,alphaYpq,alphaZpq,Xpa,Ypa,Zpa
    real(realk) :: inv2expP,TwoTerms1,TwoTerms2,TwoTerms3,TwoTerms4,TwoTerms5,TwoTerms6
    real(realk) :: TwoTerms7
    real(realk),parameter :: D1=1.0E0_realk,D2=2.0E0_realk, D5=5.0E0_realk
    WRITE(lupri,*)'VerticalRecurrence6  currentMaxAngmom=',currentMaxAngmom
    !ThetaAux(n,i+1,0,0) = Xpa*ThetaAux(n,i,0,0) + (-alpha/p*Xpq)*ThetaAux(n+1,i,0,0) 
    !                    + i/(2p) * (ThetaAux(n,i-1,0,0) - alpha/p*ThetaAux(n+1,i-1,0,0))
    !i = 5

    !NOT TESTET !!!!!
    do IPrimP = 1,nPrimP
       invexpP = D1/Pexp(iPrimP)
       inv2expP = D1/(D2*Pexp(iPrimP))   
       Xpa = Rpa(1,iPrimP)
       Ypa = Rpa(2,iPrimP)
       Zpa = Rpa(3,iPrimP)
       do iPassQ = 1,nPasses
          I = (iPrimP-1)*nPrimQ + (iPassQ-1)*nPrimQ*nPrimP
          do IPrimQ = 1,nPrimQ
             alphaP = -reducedExponents(iPrimQ,iPrimP)*invexpP
             alphaXpq = alphaP*Rpq(1,iPrimQ,iPrimP,iPassQ)
             alphaYpq = alphaP*Rpq(2,iPrimQ,iPrimP,iPassQ)
             alphaZpq = alphaP*Rpq(3,iPrimQ,iPrimP,iPassQ)
             I = I + 1
             do J=1,currentMaxAngmom
                TwoTerms1 = inv2expP*(AUXarray(21,J,I) + alphaP*AUXarray(21,J+1,I))
                TwoTerms4 = D2*inv2expP*(AUXarray(22,J,I) + alphaP*AUXarray(22,J+1,I))
                TwoTerms5 = D2*inv2expP*(AUXarray(23,J,I) + alphaP*AUXarray(23,J+1,I))
                TwoTerms6 = inv2expP*(AUXarray(24,J,I) + alphaP*AUXarray(24,J+1,I))
                TwoTerms2 = inv2expP*(AUXarray(31,J,I) + alphaP*AUXarray(31,J+1,I))
                TwoTerms7 = D2*inv2expP*(AUXarray(32,J,I) + alphaP*AUXarray(32,J+1,I))
                TwoTerms3 = inv2expP*(AUXarray(35,J,I) + alphaP*AUXarray(35,J+1,I))

                AUXarray(57,J,I) = Xpa*AUXarray(36,J,I) + alphaXpq*AUXarray(36,J+1,I) + D5*TwoTerms1
                AUXarray(58,J,I) = Ypa*AUXarray(36,J,I) + alphaYpq*AUXarray(36,J+1,I)
                AUXarray(59,J,I) = Zpa*AUXarray(36,J,I) + alphaZpq*AUXarray(36,J+1,I)
                AUXarray(60,J,I) = Ypa*AUXarray(37,J,I) + alphaYpq*AUXarray(37,J+1,I) + TwoTerms1
                AUXarray(61,J,I) = Zpa*AUXarray(37,J,I) + alphaZpq*AUXarray(37,J+1,I)
                AUXarray(62,J,I) = Zpa*AUXarray(38,J,I) + alphaZpq*AUXarray(38,J+1,I) + TwoTerms1
                AUXarray(63,J,I) = Ypa*AUXarray(39,J,I) + alphaYpq*AUXarray(39,J+1,I) + TwoTerms4
                AUXarray(64,J,I) = Zpa*AUXarray(39,J,I) + alphaZpq*AUXarray(39,J+1,I)
                AUXarray(65,J,I) = Ypa*AUXarray(41,J,I) + alphaYpq*AUXarray(41,J+1,I)
                AUXarray(66,J,I) = Zpa*AUXarray(41,J,I) + alphaZpq*AUXarray(41,J+1,I) + TwoTerms5
                AUXarray(67,J,I) = Xpa*AUXarray(46,J,I) + alphaXpq*AUXarray(46,J+1,I) + TwoTerms2
                AUXarray(68,J,I) = Zpa*AUXarray(42,J,I) + alphaZpq*AUXarray(42,J+1,I)
                AUXarray(69,J,I) = Zpa*AUXarray(43,J,I) + alphaZpq*AUXarray(43,J+1,I) + TwoTerms6
                AUXarray(70,J,I) = Ypa*AUXarray(45,J,I) + alphaYpq*AUXarray(45,J+1,I)
                AUXarray(71,J,I) = Xpa*AUXarray(50,J,I) + alphaXpq*AUXarray(50,J+1,I) + TwoTerms3
                AUXarray(72,J,I) = Xpa*AUXarray(51,J,I) + alphaXpq*AUXarray(51,J+1,I)
                AUXarray(73,J,I) = Xpa*AUXarray(52,J,I) + alphaXpq*AUXarray(52,J+1,I)
                AUXarray(74,J,I) = Xpa*AUXarray(53,J,I) + alphaXpq*AUXarray(53,J+1,I)
                AUXarray(75,J,I) = Xpa*AUXarray(54,J,I) + alphaXpq*AUXarray(54,J+1,I)
                AUXarray(76,J,I) = Xpa*AUXarray(55,J,I) + alphaXpq*AUXarray(55,J+1,I)
                AUXarray(77,J,I) = Xpa*AUXarray(56,J,I) + alphaXpq*AUXarray(56,J+1,I)
                AUXarray(78,J,I) = Ypa*AUXarray(51,J,I) + alphaYpq*AUXarray(51,J+1,I) + D5*TwoTerms2
                AUXarray(79,J,I) = Zpa*AUXarray(51,J,I) + alphaZpq*AUXarray(51,J+1,I)
                AUXarray(80,J,I) = Zpa*AUXarray(52,J,I) + alphaZpq*AUXarray(52,J+1,I) + TwoTerms2
                AUXarray(81,J,I) = Zpa*AUXarray(53,J,I) + alphaZpq*AUXarray(53,J+1,I) + TwoTerms7
                AUXarray(82,J,I) = Ypa*AUXarray(55,J,I) + alphaYpq*AUXarray(55,J+1,I) + TwoTerms3
                AUXarray(83,J,I) = Ypa*AUXarray(56,J,I) + alphaYpq*AUXarray(56,J+1,I)
                AUXarray(84,J,I) = Zpa*AUXarray(56,J,I) + alphaZpq*AUXarray(56,J+1,I) + D5*TwoTerms3
                do ipr=57,84
                   write(lupri,*)'AUXarray(',ipr,',J,I)',AUXarray(ipr,J,I)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine VerticalRecurrence6

  subroutine build_Rpa(nPrimP,Pcent,Acent,Rpa)
    implicit none
    integer,intent(in) :: nPrimP
    real(realk),intent(in) :: Pcent(3,nPrimP),Acent(3)
    real(realk),intent(inout) :: Rpa(3,nPrimP)
    !
    integer :: iPrimP
    real(realk) :: Ax,Ay,Az
    Ax = -Acent(1)
    Ay = -Acent(2)
    Az = -Acent(3)
    DO iPrimP=1, nPrimP
       Rpa(1,iPrimP) = Pcent(1,iPrimP) + Ax
       Rpa(2,iPrimP) = Pcent(2,iPrimP) + Ay
       Rpa(3,iPrimP) = Pcent(3,iPrimP) + Az
    ENDDO
  end subroutine build_Rpa

  subroutine build_QP_squaredDistance_and_Rpq(nPrimQ,nPasses,nPrimP,Qcent,Pcent,squaredDistance,Rpq)
    implicit none
    integer,intent(in) :: nPrimQ,nPasses,nPrimP
    real(realk),intent(in) :: Qcent(3,nPrimQ,nPasses),Pcent(3,nPrimP)
    real(realk),intent(inout) :: squaredDistance(nPrimQ,nPrimP,nPasses)
    real(realk),intent(inout) :: Rpq(3,nPrimQ,nPrimP,nPasses)
    !
    integer :: iPrimQ,iPassQ,iPrimP
    real(realk) :: px,py,pz,pqx,pqy,pqz
    DO iPassQ=1, nPasses
       DO iPrimP=1, nPrimP
          px = Pcent(1,iPrimP)
          py = Pcent(2,iPrimP)
          pz = Pcent(3,iPrimP)
          DO iPrimQ=1, nPrimQ
             pqx = px - Qcent(1,iPrimQ,iPassQ)
             pqy = py - Qcent(2,iPrimQ,iPassQ)
             pqz = pz - Qcent(3,iPrimQ,iPassQ)
             rPQ(1,iPrimQ,iPrimP,iPassQ) = pqx
             rPQ(2,iPrimQ,iPrimP,iPassQ) = pqy
             rPQ(3,iPrimQ,iPrimP,iPassQ) = pqz
             squaredDistance(iPrimQ,iPrimP,iPassQ) = pqx*pqx+pqy*pqy+pqz*pqz
          ENDDO
       ENDDO
    ENDDO
  end subroutine build_QP_squaredDistance_and_Rpq

  SUBROUTINE buildRJ000_general(nPasses,nPrimPQ,nTABFJW1,nTABFJW2,reducedExponents,&
       & squaredDistance,TABFJW,RJ000,JMAX)
    IMPLICIT NONE
    INTEGER,intent(in)         :: nPrimPQ,Jmax,nTABFJW1,nTABFJW2,nPasses
    REAL(REALK),intent(in)     :: reducedExponents(nPrimPQ)
    REAL(REALK),intent(in)     :: squaredDistance(nPrimPQ,nPasses)
    REAL(REALK),intent(in)     :: TABFJW(0:nTABFJW1,0:nTABFJW2)
    REAL(REALK),intent(inout)  :: RJ000(0:Jmax,nPrimPQ,nPasses)
    !
    REAL(REALK)     :: D2JP36,WVAL
    REAL(REALK),PARAMETER :: HALF =0.5E0_realk,D1=1E0_realk,D2 = 2E0_realk, D4 = 4E0_realk, D100=100E0_realk
    Real(realk),parameter :: D12 = 12E0_realk, TENTH = 0.01E0_realk
    REAL(REALK), PARAMETER :: COEF2 = HALF,  COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk
    REAL(REALK), PARAMETER :: COEF5 = - D1/120E0_realk, COEF6 = D1/720E0_realk
    Integer :: IPNT,J
    Real(realk) :: WDIFF,RWVAL,REXPW,GVAL,PREF,D2MALPHA
    REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092E0_realk
    REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686E0_realk
    REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909E0_realk
    REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346E0_realk
    Real(realk), parameter :: PI=3.14159265358979323846E0_realk
    REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk
    REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2
    REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI
    Real(realk) :: W2,W3,R
    REAL(REALK), PARAMETER :: SMALL = 1E-15_realk
    integer :: iPassQ,iPQ
    !make for different values of JMAX => loop unroll  
    !sorting? 
    D2JP36 = 2*JMAX + 36
    DO iPassq=1, nPasses
       DO ipq=1, nPrimPQ
          !(nPrimP,nPrimQ)*(nPrimP,nPrimQ,nPasses)
          WVAL = reducedExponents(ipq)*squaredDistance(ipq,iPassQ)
          !  0 < WVAL < 0.000001
          IF (ABS(WVAL) .LT. SMALL) THEN         
             RJ000(0,ipq,ipassq) = D1
             DO J=1,JMAX
                RJ000(J,ipq,ipassq)= D1/(2*J + 1) !THE BOYS FUNCTION FOR ZERO ARGUMENT
             ENDDO
             !  0 < WVAL < 12 
          ELSE IF (WVAL .LT. D12) THEN
             IPNT = NINT(D100*WVAL)
             WDIFF = WVAL - TENTH*IPNT
             W2    = WDIFF*WDIFF
             W3    = W2*WDIFF
             W2    = W2*COEF2
             W3    = W3*COEF3
             DO J=0,JMAX
                R = TABFJW(J,IPNT)
                R = R -TABFJW(J+1,IPNT)*WDIFF
                R = R + TABFJW(J+2,IPNT)*W2
                R = R + TABFJW(J+3,IPNT)*W3
                RJ000(J,ipq,ipassq) = R
             ENDDO
             !  12 < WVAL <= (2J+36) 
          ELSE IF (WVAL.LE.D2JP36) THEN
             REXPW = HALF*EXP(-WVAL)
             RWVAL = D1/WVAL
             GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
             RJ000(0,ipq,ipassq) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
             DO J=1,JMAX
                RJ000(J,ipq,ipassq) = RWVAL*((J - HALF)*RJ000(J-1,ipq,ipassq)-REXPW)
             ENDDO
             !  (2J+36) < WVAL 
          ELSE
             RWVAL = PID4/WVAL
             RJ000(0,ipq,ipassq) = SQRT(RWVAL)
             RWVAL = RWVAL*PID4I
             DO J = 1, JMAX
                RJ000(J,ipq,ipassq) = RWVAL*(J - HALF)*RJ000(J-1,ipq,ipassq)
             ENDDO
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE buildRJ000_general

  SUBROUTINE scaleRJ000_general(nPasses,nPrimP,nPrimQ,reducedExponents,&
       & RJ000,JMAX,integralPrefactor,AUXarray,nTUV,PpreExpFac,QpreExpFac,lupri)
    IMPLICIT NONE
    INTEGER,intent(in)         :: nPrimP,nPrimQ,Jmax,nPasses,nTUV,lupri
    REAL(REALK),intent(in)     :: QpreExpFac(nPrimQ,nPasses),PpreExpFac(nPrimP)
    REAL(REALK),intent(in)     :: reducedExponents(nPrimQ,nPrimP)
    REAL(REALK),intent(in)     :: integralPrefactor(nprimQ,nPrimP)
    REAL(REALK),intent(in)     :: RJ000(0:Jmax,nPrimQ,nPrimP,nPasses)
    REAL(REALK),intent(inout)  :: AUXarray(nTUV,JMAX+1,nPrimQ,nPrimP,nPasses)
    !
    Real(realk) :: PREF,Pexpfac
    integer :: iPassQ,iP,ipq,iq,pqoffset,j
    ! Scaling
    DO ip=1, nPrimP
       pqoffset = (ip-1)*nPrimQ
       Pexpfac = PpreExpFac(iP)
       DO iPassq=1, nPasses
          DO iq=1, nPrimQ
             PREF = integralPrefactor(iq,ip)*QpreExpFac(iq,iPassQ)*Pexpfac
             DO j=0,jmax
                AUXarray(1,J+1,iq,ip,iPassq) = PREF*RJ000(J,iq,ip,iPassq)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE scaleRJ000_general

end MODULE IchorEriCoulombintegralOBSGeneralMod
