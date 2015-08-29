!> @file
!> 

module f12ri_util_module
  use precision
  use memory_handling
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type,only: ls_mpibcast
#endif
  !  use lstiming

  private
  public :: GeneralTwo4CenterF12RICoef1112,&
       & GeneralTwo4CenterF12RICoef1223,&
       & F12RIB4,F12RIB4MPI,F12RIB6,F12RIB6MPI,F12RIB9,F12RIB9MPI,&
       & F12RIB5,F12RIB5MPI,F12RIB7,F12RIB7MPI,F12RIB8,F12RIB8MPI

contains

  ! (Calpha1*Calpha1)*(Calpha1*Calpha2)
  ! Bcast Calpha1
  subroutine GeneralTwo4CenterF12RICoef1112(nBA,&
       & Calpha1,ndim12,ndim13,Calpha2,ndim22,ndim23,&
       & Energy,offset,SlavesAwake,use_bg_buf,&
       & numnodes,nAuxMPI,mynum,InputSubroutine,InputSubroutineMPI)
    implicit none
    integer :: ndim12,ndim13,ndim22,ndim23,NBA
    integer :: numnodes,mynum,offset
    real(realk) :: Calpha1(NBA*ndim12*ndim13)
    real(realk) :: Calpha2(NBA*ndim22*ndim23)
    real(realk) :: Energy
    logical :: SlavesAwake,use_bg_buf
    integer :: nAuxMPI(numnodes)    
    EXTERNAL InputSubroutine !NAME OF SUBROUTINE TO CALL
    EXTERNAL InputSubroutineMPI !NAME OF SUBROUTINE TO CALL
    !local variables
#ifdef VAR_MPI
    integer :: inode,nbuf1,NBA2
    integer(kind=ls_mpik) :: node
    integer(kind=long) :: nsize1
    real(realk) :: EnergyTmp
    real(realk),pointer :: CalphaMPI1(:)
#endif
    IF(SlavesAwake)THEN
#ifdef VAR_MPI
       Energy  = 0.0E0_realk
       DO inode = 1,numnodes
          nbuf1 = nAuxMPI(inode)
          NBA2 = nAuxMPI(inode)           !dimension of chunk we are working on 
          nsize1 = nbuf1*ndim12*ndim13    !Calpha1(NBA,ndim12,ndim13)
          IF(mynum.EQ.inode-1)THEN
             !I Bcast My Own Calpha1
             node = mynum            
             call ls_mpibcast(Calpha1,nsize1,node,infpar%lg_comm)
             call InputSubroutine(nBA,Calpha1,ndim12,ndim13,Calpha2,ndim22,ndim23,EnergyTmp,offset)
          ELSE
             node = inode-1
             !recieve Calpha1 from a different node
             IF(use_bg_buf)THEN
                call mem_pseudo_alloc(CalphaMPI1,nsize1) 
             ELSE
                call mem_alloc(CalphaMPI1,nsize1)
             ENDIF
             call ls_mpibcast(CalphaMPI1,nsize1,node,infpar%lg_comm)
             call InputSubroutineMPI(nBA,Calpha1,ndim12,ndim13,Calpha2,ndim22,ndim23,&
                  & EnergyTmp,offset,NBA2,CalphaMPI1)
             IF(use_bg_buf)THEN
                call mem_pseudo_dealloc(CalphaMPI1)
             ELSE
                call mem_dealloc(CalphaMPI1)
             ENDIF
          ENDIF
          Energy = Energy + EnergyTmp
       ENDDO
#else
       call lsquit('GeneralTwo4CenterF12RICoef1112 SlavesAwake=.TRUE. require MPI',-1)
#endif
    ELSE
       call InputSubroutine(nBA,Calpha1,ndim12,ndim13,Calpha2,ndim22,ndim23,Energy,offset)       
    ENDIF
  end subroutine GeneralTwo4CenterF12RICoef1112

  ! (Calpha1*Calpha2)*(Calpha2*Calpha3)
  ! Bcast Calpha1,Calpha2
  subroutine GeneralTwo4CenterF12RICoef1223(nBA,&
       & Calpha1,ndim12,ndim13,Calpha2,ndim22,ndim23,&
       & Calpha3,ndim32,ndim33,Energy,offset,SlavesAwake,use_bg_buf,&
       & numnodes,nAuxMPI,mynum,InputSubroutine,InputSubroutineMPI)
    implicit none
    integer :: ndim12,ndim13,ndim22,ndim23,ndim32,ndim33,NBA
    real(realk) :: Calpha1(NBA*ndim12*ndim13)
    real(realk) :: Calpha2(NBA*ndim22*ndim23)
    real(realk) :: Calpha3(NBA*ndim32*ndim33)
    real(realk) :: Energy
    logical :: SlavesAwake,use_bg_buf
    integer :: numnodes,mynum,offset    
    integer :: nAuxMPI(numnodes)
    EXTERNAL InputSubroutine !NAME OF SUBROUTINE TO CALL
    EXTERNAL InputSubroutineMPI !NAME OF SUBROUTINE TO CALL
    !local variables
#ifdef VAR_MPI
    integer :: nbuf1,NBA2,inode
    integer(kind=ls_mpik) :: node
    integer(kind=long) :: nsize1,nsize2
    real(realk) :: EnergyTmp
    real(realk),pointer :: CalphaMPI1(:),CalphaMPI2(:)
#endif
    IF(SlavesAwake)THEN
#ifdef VAR_MPI
       Energy  = 0.0E0_realk
       DO inode = 1,numnodes
          nbuf1 = nAuxMPI(inode)
          NBA2 = nAuxMPI(inode)           !dimension of chunk we are working on 
          nsize1 = nbuf1*ndim12*ndim13    !Calpha1(NBA,ndim12,ndim13)
          nsize2 = nbuf1*ndim22*ndim23    !Calpha2(NBA,ndim22,ndim23)          
          IF(mynum.EQ.inode-1)THEN
             !I Bcast My Own Calpha1,Calpha2
             node = mynum            
             call ls_mpibcast(Calpha1,nsize1,node,infpar%lg_comm)
             call ls_mpibcast(Calpha2,nsize2,node,infpar%lg_comm)
             call InputSubroutine(nBA,Calpha1,ndim12,ndim13,Calpha2,ndim22,ndim23,&
                  & Calpha3,ndim32,ndim33,EnergyTmp,offset)
          ELSE
             node = inode-1
             !recieve Calpha1 and Calpha2 from a different node
             IF(use_bg_buf)THEN
                call mem_pseudo_alloc(CalphaMPI1,nsize1) 
                call mem_pseudo_alloc(CalphaMPI2,nsize2)
             ELSE
                call mem_alloc(CalphaMPI1,nsize1)
                call mem_alloc(CalphaMPI2,nsize2)
             ENDIF
             call ls_mpibcast(CalphaMPI1,nsize1,node,infpar%lg_comm)
             call ls_mpibcast(CalphaMPI2,nsize2,node,infpar%lg_comm)
             call InputSubroutineMPI(nBA,Calpha1,ndim12,ndim13,Calpha2,ndim22,ndim23,&
                  & Calpha3,ndim32,ndim33,EnergyTmp,offset,NBA2,CalphaMPI1,CalphaMPI2)
             IF(use_bg_buf)THEN
                call mem_pseudo_dealloc(CalphaMPI2)
                call mem_pseudo_dealloc(CalphaMPI1)
             ELSE
                call mem_dealloc(CalphaMPI1)
                call mem_dealloc(CalphaMPI2)
             ENDIF
          ENDIF
          Energy = Energy + EnergyTmp
       ENDDO
#else
       call lsquit('GeneralTwo4CenterF12RICoef1112 SlavesAwake=.TRUE. require MPI',-1)
#endif
    ELSE
       call InputSubroutine(nBA,Calpha1,ndim12,ndim13,Calpha2,ndim22,ndim23,&
            & Calpha3,ndim32,ndim33,Energy,offset)
    ENDIF
  end subroutine GeneralTwo4CenterF12RICoef1223

  !============================================================================
  !   Subroutines that can be given as input to GeneralTwo4CenterF12RI3Coef1112
  !============================================================================

  !( CalphaG(beta,i,r)*CalphaG(beta,j,q) ) * ( CalphaG(beta,j,r)*CalphaD(beta,i,q))
  subroutine F12RIB4(nBA,CalphaG,n1,n2,CalphaD,nD1,nD2,EJK,offset)
    implicit none
    integer,intent(in)        :: nBA,n1,n2,offset,nD1,nD2
    real(realk),intent(in)    :: CalphaG(nBA,n1,n2)
    real(realk),intent(in)    :: CalphaD(nBa,n1,n2)
    real(realk),intent(inout) :: EJK
    !local variables
    real(realk)               :: EJ, EK, ED
    integer :: p,q,r,i,j,beta,alpha
    real(realk) :: tmpR,tmpRJ1,tmpRJ2
    real(realk) :: tmpG,tmpGJ1,tmpGJ2
    ED = 0.0E0_realk
    EJ = 0.0E0_realk
    EK = 0.0E0_realk
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(p,q,r,i,j,beta,alpha,&
    !$OMP tmpR,tmpRJ1,tmpRJ2,tmpG,tmpGJ1,tmpGJ2) SHARED(nBA,n1,n2,&
    !$OMP CalphaG,CalphaD) REDUCTION(+:ED,EJ,EK)
    DO q=1,n2 !ncabsAO
       DO r=1,n2 !ncabsAO
          DO j=1,n1 !nocc
             !Diagonal
             tmpR = 0.0E0_realk
             tmpG = 0.0E0_realk
             DO beta = 1,nBA
                tmpG = tmpG + CalphaG(beta,j,r)*CalphaG(beta,j,q)
                tmpR = tmpR + CalphaG(beta,j,r)*CalphaD(beta,j,q)
             ENDDO
             ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
             !Non Diagonal
             DO i=j+1,n1
                tmpRJ1 = 0.0E0_realk
                tmpGJ1 = 0.0E0_realk
                tmpRJ2 = 0.0E0_realk
                tmpGJ2 = 0.0E0_realk
                DO alpha = 1, nBA
                   tmpGJ1 = tmpGJ1 + CalphaG(alpha,i,r)*CalphaG(alpha,j,q)
                   tmpGJ2 = tmpGJ2 + CalphaG(alpha,j,r)*CalphaG(alpha,i,q)
                   tmpRJ1 = tmpRJ1 + CalphaG(alpha,i,r)*CalphaD(alpha,j,q) 
                   tmpRJ2 = tmpRJ2 + CalphaG(alpha,j,r)*CalphaD(alpha,i,q)
                ENDDO
                EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    EJK = ED*0.5E0_realk + 7.0E0_realk/16.0E0_realk*EJ + 1.0E0_realk/16.0E0_realk*EK 
  end subroutine F12RIB4

  !( CalphaG(beta,i,r)*CalphaG(beta,j,q) ) * ( CalphaG(beta,j,r)*CalphaD(beta,i,q))
  subroutine F12RIB4MPI(nBA,CalphaG,n1,n2,CalphaD,nD1,nD2,EJK,offset,NBA2,&
       & CalphaGMPI)
    implicit none
    integer,intent(in)        :: nBA,n1,n2,NBA2,offset,nD1,nD2
    real(realk),intent(in)    :: CalphaGMPI(nBA2,n1,n2)
    real(realk),intent(in)    :: CalphaG(nBA,n1,n2)
    real(realk),intent(in)    :: CalphaD(nBa,n1,n2)
    real(realk),intent(inout) :: EJK
    !local variables
    real(realk)               :: EJ, EK, ED
    integer :: q,r,i,j,beta
    real(realk) :: tmpR,tmpRJ1,tmpRJ2
    real(realk) :: tmpG,tmpGJ1,tmpGJ2
    ED = 0.0E0_realk
    EJ = 0.0E0_realk
    EK = 0.0E0_realk
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(q,r,i,j,beta,&
    !$OMP tmpR,tmpRJ1,tmpRJ2,tmpG,tmpGJ1,tmpGJ2) SHARED(nBA,n1,n2,CalphaGMPI,&
    !$OMP NBA2,CalphaG,CalphaD) REDUCTION(+:ED,EJ,EK)
    DO q=1,n2 !ncabsAO
       DO r=1,n2 !ncabsAO
          DO j=1,n1 !nocc
             !Diagonal
             tmpG = 0.0E0_realk
             DO beta = 1,nBA2
                tmpG = tmpG + CalphaGMPI(beta,j,r)*CalphaGMPI(beta,j,q)
             ENDDO
             tmpR = 0.0E0_realk
             DO beta = 1,nBA
                tmpR = tmpR + CalphaG(beta,j,r)*CalphaD(beta,j,q)
             ENDDO
             ED = ED + tmpR*tmpG !We have a factor 2 which is integrated
             !Non Diagonal
             DO i=j+1,n1
                tmpGJ1 = 0.0E0_realk
                tmpGJ2 = 0.0E0_realk
                DO beta = 1, nBA2
                   tmpGJ1 = tmpGJ1 + CalphaGMPI(beta,i,r)*CalphaGMPI(beta,j,q)
                   tmpGJ2 = tmpGJ2 + CalphaGMPI(beta,j,r)*CalphaGMPI(beta,i,q)
                ENDDO
                tmpRJ1 = 0.0E0_realk
                tmpRJ2 = 0.0E0_realk
                DO beta = 1, nBA
                   tmpRJ1 = tmpRJ1 + CalphaG(beta,i,r)*CalphaD(beta,j,q) 
                   tmpRJ2 = tmpRJ2 + CalphaG(beta,j,r)*CalphaD(beta,i,q)
                ENDDO
                EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)             
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    EJK = ED*0.5E0_realk + 7.0E0_realk/16.0E0_realk*EJ + 1.0E0_realk/16.0E0_realk*EK 
  end subroutine F12RIB4MPI

  !(CalphaG(alpha,q,i)*CalphaG(alpha,a,j)) *(CalphaD(alpha,i,q)*CalphaG(alpha,a,j))
  subroutine F12RIB6(nBA,CalphaG,n3,n1,CalphaD,nD1,nD2,EJK,noccfull)
    implicit none
    integer,intent(in)        :: nBA,n1,n3,noccfull,nD1,nD2
    real(realk),intent(in)    :: CalphaG(nBA,n3,n1)
    real(realk),intent(in)    :: CalphaD(nBA,n1,n3)
    real(realk),intent(inout) :: EJK
    !local variables
    real(realk)               :: ED, EJ, EK
    integer :: q,a,i,j,alpha
    real(realk) :: tmpR,tmpRJ1,tmpRJ2
    real(realk) :: tmpG,tmpGJ1,tmpGJ2
    ED = 0.0E0_realk
    EK = 0.0E0_realk
    EJ = 0.0E0_realk
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(q,a,i,j,alpha,&
    !$OMP tmpR,tmpRJ1,tmpRJ2,tmpG,tmpGJ1,tmpGJ2) SHARED(nBA,n1,n3,&
    !$OMP noccfull,CalphaG,CalphaD,EJK) REDUCTION(+:ED,EJ,EK)
    DO q=1,n3 !nocv 
       DO a=noccfull+1,n3 !nvirt
          DO j=1,n1 !nocc       
             !Diagonal          
             tmpR = 0.0E0_realk   
             tmpG = 0.0E0_realk   
             DO alpha = 1,nBA     
                tmpG = tmpG + CalphaG(alpha,q,j)*CalphaG(alpha,a,j)
                tmpR = tmpR + CalphaD(alpha,j,q)*CalphaG(alpha,a,j)
             ENDDO
             ED = ED + tmpR*tmpG !We have a factor 2 which is integrated
             !Non Diagonal
             DO i=j+1,n1
                tmpRJ1 = 0.0E0_realk 
                tmpGJ1 = 0.0E0_realk 
                tmpRJ2 = 0.0E0_realk                                    
                tmpGJ2 = 0.0E0_realk                                    
                DO alpha = 1, nBA   
                   tmpGJ1 = tmpGJ1 + CalphaG(alpha,q,i)*CalphaG(alpha,a,j)
                   tmpGJ2 = tmpGJ2 + CalphaG(alpha,q,j)*CalphaG(alpha,a,i)
                   tmpRJ1 = tmpRJ1 + CalphaD(alpha,i,q)*CalphaG(alpha,a,j)
                   tmpRJ2 = tmpRJ2 + CalphaD(alpha,j,q)*CalphaG(alpha,a,i)
                ENDDO
                EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)                   
                EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)                   
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    EJK = -1.0E0_realk*(ED*0.5E0_realk + 7.0E0_realk/16.0E0_realk*EJ + 1.0E0_realk/16.0E0_realk*EK) 
  end subroutine F12RIB6

  subroutine F12RIB6MPI(nBA,CalphaG,n3,n1,CalphaD,nD1,nD3,EJK,noccfull,&
       & NBA2,CalphaGMPI)
    implicit none
    integer,intent(in)        :: nBA,n1,n3,noccfull,NBA2,nD1,nD3
    real(realk),intent(in)    :: CalphaGMPI(nBA2,n3,n1)
    real(realk),intent(in)    :: CalphaG(nBA,n3,n1)
    real(realk),intent(in)    :: CalphaD(nBA,n1,n3)
    real(realk),intent(inout) :: EJK
    !local variables
    real(realk)               :: EJ, EK, ED
    integer :: q,a,i,j,alpha
    real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
    real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2
    ED = 0.0E0_realk
    EK = 0.0E0_realk
    EJ = 0.0E0_realk
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(q,a,i,j,alpha,&
    !$OMP tmpR,tmpRJ1,tmpRJ2,tmpG,tmpGJ1,tmpGJ2) SHARED(nBA,n1,n3,&
    !$OMP CalphaGMPI,NBA2,noccfull,CalphaG,CalphaD,EJK) REDUCTION(+:ED,EJ,EK)
    DO q=1,n3 !nocv 
       DO a=noccfull+1,n3 !nvirt
          DO j=1,n1 !nocc       
             !Diagonal          
             tmpG = 0.0E0_realk   
             DO alpha = 1,nBA2
                tmpG = tmpG + CalphaGMPI(alpha,q,j)*CalphaGMPI(alpha,a,j)
             ENDDO
             tmpR = 0.0E0_realk   
             DO alpha = 1,nBA     
                tmpR = tmpR + CalphaD(alpha,j,q)*CalphaG(alpha,a,j)
             ENDDO
             ED = ED + tmpR*tmpG !We have a factor 2 which is integrated          
             !Non Diagonal
             DO i=j+1,n1
                tmpGJ1 = 0.0E0_realk 
                tmpGJ2 = 0.0E0_realk                                    
                DO alpha = 1,nBA2
                   tmpGJ1 = tmpGJ1 + CalphaGMPI(alpha,q,i)*CalphaGMPI(alpha,a,j)
                   tmpGJ2 = tmpGJ2 + CalphaGMPI(alpha,q,j)*CalphaGMPI(alpha,a,i)
                ENDDO
                tmpRJ1 = 0.0E0_realk 
                tmpRJ2 = 0.0E0_realk                                    
                DO alpha = 1, nBA   
                   tmpRJ1 = tmpRJ1 + CalphaD(alpha,i,q)*CalphaG(alpha,a,j)
                   tmpRJ2 = tmpRJ2 + CalphaD(alpha,j,q)*CalphaG(alpha,a,i)
                ENDDO
                EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)                   
                EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)                   
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    EJK = -1.0E0_realk*(ED*0.5E0_realk + 7.0E0_realk/16.0E0_realk*EJ + 1.0E0_realk/16.0E0_realk*EK)

  end subroutine F12RIB6MPI

  !((CalphaG(alpha,p,i)*CalphaG(alpha,a,j)) *  (CalphaD(alpha,i,p)*CalphaG(alpha,a,j)))
  subroutine F12RIB9(nBA,CalphaG,n3,n1,CalphaD,nD1,nD2,EJK,noccfull)
    implicit none
    integer,intent(in)        :: nBA,n1,n3,noccfull,nD1,nD2
    real(realk),intent(in)    :: CalphaG(nBA,n3,n1)
    real(realk),intent(in)    :: CalphaD(nBA,n1,n3)
    real(realk),intent(inout) :: EJK
    !local variables
    real(realk)               :: EJ, EK, ED
    integer :: p,a,i,j,alpha
    real(realk) :: tmpR,tmpRJ1,tmpRJ2
    real(realk) :: tmpG,tmpGJ1,tmpGJ2
    ED = 0.0E0_realk
    EK = 0.0E0_realk
    EJ = 0.0E0_realk
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(p,a,i,j,alpha,&
    !$OMP tmpR,tmpRJ1,tmpRJ2,tmpG,tmpGJ1,tmpGJ2) SHARED(nBA,&
    !$OMP noccfull,n1,n3,CalphaG,CalphaD) REDUCTION(+:ED,EJ,EK)
    DO p=1,n3 !ncabs
       DO a=noccfull+1,n3 !nvirt
          DO j=1,n1 !nocc
             !Diagonal
             tmpR = 0.0E0_realk
             tmpG = 0.0E0_realk
             DO alpha = 1,nBA
                tmpR = tmpR + CalphaG(alpha,p,j)*CalphaG(alpha,a,j)
                tmpG = tmpG + CalphaD(alpha,j,p)*CalphaG(alpha,a,j)
             ENDDO
             ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
             !Non Diagonal
             DO i=j+1,n1
                tmpRJ1 = 0.0E0_realk
                tmpRJ2 = 0.0E0_realk
                tmpGJ1 = 0.0E0_realk
                tmpGJ2 = 0.0E0_realk
                DO alpha = 1, nBA
                   tmpRJ1 = tmpRJ1 + CalphaG(alpha,p,i)*CalphaG(alpha,a,j) 
                   tmpRJ2 = tmpRJ2 + CalphaG(alpha,p,j)*CalphaG(alpha,a,i)
                   tmpGJ1 = tmpGJ1 + CalphaD(alpha,i,p)*CalphaG(alpha,a,j)
                   tmpGJ2 = tmpGJ2 + CalphaD(alpha,j,p)*CalphaG(alpha,a,i)
                ENDDO
                EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)          
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    EJK = -2.0E0_realk*(ED*0.5E0_realk + 7.0_realk/16.0_realk*EJ + 1.0_realk/16.0_realk*EK)
  end subroutine F12RIB9

  !((CalphaG(alpha,p,i)*CalphaG(alpha,a,j)) *  (CalphaD(alpha,i,p)*CalphaG(alpha,a,j)))
  subroutine F12RIB9MPI(nBA,CalphaG,n3,n1,CalphaD,nD1,nD2,EJK,noccfull,&
       & NBA2,CalphaGMPI)
    implicit none
    integer,intent(in)        :: nBA,n1,n3,noccfull,NBA2,nD1,nD2
    real(realk),intent(in)    :: CalphaGMPI(nBA2,n3,n1)
    real(realk),intent(in)    :: CalphaG(nBA,n3,n1)
    real(realk),intent(in)    :: CalphaD(nBA,n1,n3)
    real(realk),intent(inout) :: EJK
    !local variables
    real(realk)               :: EJ, EK, ED
    integer :: p,a,i,j,alpha
    real(realk) :: tmpR,tmpRJ1,tmpRJ2
    real(realk) :: tmpG,tmpGJ1,tmpGJ2
    ED = 0.0E0_realk
    EK = 0.0E0_realk
    EJ = 0.0E0_realk
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(p,a,i,j,alpha,&
    !$OMP tmpR,tmpRJ1,tmpRJ2,tmpG,tmpGJ1,tmpGJ2) SHARED(nBA,&
    !$OMP noccfull,n1,n3,CalphaGMPI,NBA2,CalphaG,CalphaD) REDUCTION(+:ED,EJ,EK)
    DO p=1,n3 !ncabs
       DO a=noccfull+1,n3 !nvirt
          DO j=1,n1 !nocc
             !Diagonal
             tmpR = 0.0E0_realk
             DO alpha = 1,nBA2
                tmpR = tmpR + CalphaGMPI(alpha,p,j)*CalphaGMPI(alpha,a,j)
             ENDDO
             tmpG = 0.0E0_realk
             DO alpha = 1,nBA
                tmpG = tmpG + CalphaD(alpha,j,p)*CalphaG(alpha,a,j)
             ENDDO
             ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
             !Non Diagonal
             DO i=j+1,n1
                tmpRJ1 = 0.0E0_realk
                tmpRJ2 = 0.0E0_realk
                DO alpha = 1, nBA2
                   tmpRJ1 = tmpRJ1 + CalphaGMPI(alpha,p,i)*CalphaGMPI(alpha,a,j) 
                   tmpRJ2 = tmpRJ2 + CalphaGMPI(alpha,p,j)*CalphaGMPI(alpha,a,i)
                ENDDO
                tmpGJ1 = 0.0E0_realk
                tmpGJ2 = 0.0E0_realk
                DO alpha = 1, nBA
                   tmpGJ1 = tmpGJ1 + CalphaD(alpha,i,p)*CalphaG(alpha,a,j)
                   tmpGJ2 = tmpGJ2 + CalphaD(alpha,j,p)*CalphaG(alpha,a,i)
                ENDDO
                EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    EJK = -2.0E0_realk*(ED*0.5E0_realk + 7.0_realk/16.0_realk*EJ + 1.0_realk/16.0_realk*EK)
  end subroutine F12RIB9MPI
  
  !============================================================================
  !   Subroutines that can be given as input to GeneralTwo4CenterF12RI3Coef1223
  !============================================================================
  !CalphaD(alpha1,i,q)*CalphaG(alpha1,m,j) * CalphaG(alpha1,m,j)*CalphaGcabs(alpha1,i,q)
  subroutine F12RIB5(nBA,CalphaD,n1,n2,CalphaG,nbasis,nG1,&
       & CalphaGcabs,nGcabs1,nGcabs2,EJK,n3)
    implicit none
    integer,intent(in)        :: nBA,n1,n2,nbasis,nGcabs1,nGcabs2,nG1,n3
    real(realk),intent(in)    :: CalphaD(nBA,n1,n2)
    real(realk),intent(in)    :: CalphaG(nBA,nbasis,n1)
    real(realk),intent(in)    :: CalphaGcabs(nBA,n1,n2)
    real(realk),intent(inout) :: EJK
    !local variables
    real(realk)               :: EJ, EK, ED
    integer :: q,m,i,j,alpha
    real(realk) :: tmpR,tmpRJ1,tmpRJ2
    real(realk) :: tmpG,tmpGJ1,tmpGJ2
    ED = 0.0E0_realk
    EJ = 0.0E0_realk
    EK = 0.0E0_realk
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(q,m,i,j,alpha,&
    !$OMP tmpR,tmpRJ1,tmpRJ2,tmpG,tmpGJ1,tmpGJ2) SHARED(nBA,&
    !$OMP n1,n2,n3,nbasis,CalphaGcabs,CalphaG,CalphaD) REDUCTION(+:ED,EJ,EK)
    DO q=1,n2 !ncabsAO
       DO m=1,n3 !noccAOS
          DO j=1,n1 !nocc
             !Diagonal
             tmpR = 0.0E0_realk
             tmpG = 0.0E0_realk
             DO alpha = 1,nBA
                tmpR = tmpR + CalphaD(alpha,j,q)*CalphaG(alpha,m,j)
                tmpG = tmpG + CalphaGcabs(alpha,j,q)*CalphaG(alpha,m,j)
             ENDDO
             ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
             !Non Diagonal
             DO i=j+1,n1
                tmpRJ1 = 0.0E0_realk
                tmpGJ1 = 0.0E0_realk
                tmpRJ2 = 0.0E0_realk
                tmpGJ2 = 0.0E0_realk
                DO alpha = 1, nBA
                   tmpRJ1 = tmpRJ1 + CalphaD(alpha,i,q)*CalphaG(alpha,m,j) 
                   tmpRJ2 = tmpRJ2 + CalphaD(alpha,j,q)*CalphaG(alpha,m,i)
                   tmpGJ1 = tmpGJ1 + CalphaGcabs(alpha,i,q)*CalphaG(alpha,m,j)
                   tmpGJ2 = tmpGJ2 + CalphaGcabs(alpha,j,q)*CalphaG(alpha,m,i)
                ENDDO
                EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    EJK = -1.0E0_realk*(ED*0.5E0_realk + 7.0E0_realk/16.0E0_realk*EJ + 1.0E0_realk/16.0E0_realk*EK) 
  end subroutine F12RIB5

  !CalphaD(alpha1,i,q)*CalphaG(alpha1,m,j) * CalphaG(alpha1,m,j)*CalphaGcabs(alpha1,i,q)
  subroutine F12RIB5MPI(nBA,CalphaD,n1,n2,CalphaG,nbasis,nG1,&
       & CalphaGcabs,nGcabs1,nGcabs2,EJK,n3,NBA2,CalphaDMPI,CalphaGMPI)
    implicit none
    integer,intent(in)        :: nBA,n1,n2,n3,nbasis,NBA2,nG1,nGcabs1,nGcabs2
    real(realk),intent(in)    :: CalphaD(nBA,n1,n2)
    real(realk),intent(in)    :: CalphaG(nBA,nbasis,n1)
    real(realk),intent(in)    :: CalphaGcabs(nBA,n1,n2)
    real(realk),intent(in)    :: CalphaDMPI(nBA2,n1,n2)
    real(realk),intent(in)    :: CalphaGMPI(nBA2,nbasis,n1)
    real(realk),intent(inout) :: EJK
    !local variables
    real(realk)               :: EJ, EK, ED
    integer :: p,q,m,i,j,alpha
    real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
    real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2
    ED = 0.0E0_realk
    EJ = 0.0E0_realk
    EK = 0.0E0_realk
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(q,m,i,j,alpha,&
    !$OMP tmpR,tmpRJ1,tmpRJ2,tmpG,tmpGJ1,tmpGJ2) SHARED(nBA,&
    !$OMP n1,n2,n3,nbasis,CalphaGcabs,CalphaG,CalphaD,CalphaGMPI,CalphaDMPI,&
    !$OMP NBA2) REDUCTION(+:ED,EJ,EK)
    DO q=1,n2 !ncabsAO
       DO m=1,n3 !noccAOS
          DO j=1,n1 !nocc
             !Diagonal
             tmpR = 0.0E0_realk
             tmpG = 0.0E0_realk
             DO alpha = 1,nBA2
                tmpR = tmpR + CalphaDMPI(alpha,j,q)*CalphaGMPI(alpha,m,j)
             ENDDO
             DO alpha = 1,nBA
                tmpG = tmpG + CalphaGcabs(alpha,j,q)*CalphaG(alpha,m,j)
             ENDDO
             ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
             !Non Diagonal
             DO i=j+1,n1
                tmpRJ1 = 0.0E0_realk
                tmpRJ2 = 0.0E0_realk
                tmpGJ1 = 0.0E0_realk
                tmpGJ2 = 0.0E0_realk
                DO alpha = 1, nBA2
                   tmpRJ1 = tmpRJ1 + CalphaDMPI(alpha,i,q)*CalphaGMPI(alpha,m,j) 
                   tmpRJ2 = tmpRJ2 + CalphaDMPI(alpha,j,q)*CalphaGMPI(alpha,m,i)
                ENDDO
                DO alpha = 1, nBA
                   tmpGJ1 = tmpGJ1 + CalphaGcabs(alpha,i,q)*CalphaG(alpha,m,j)
                   tmpGJ2 = tmpGJ2 + CalphaGcabs(alpha,j,q)*CalphaG(alpha,m,i)
                ENDDO
                EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    EJK = -1.0E0_realk*(ED*0.5E0_realk + 7.0E0_realk/16.0E0_realk*EJ + 1.0E0_realk/16.0E0_realk*EK) 
  end subroutine F12RIB5MPI

  !CalphaD(alpha1,j,n)*CalphaR(alpha1,i,c) * CalphaR(beta1,i,c)*CalphaG(beta1,j,n)
  subroutine F12RIB7(nBA,CalphaD,n1,n2,CalphaR,nR1,n3,CalphaG,nG1,nG2,&
       & EJK,offset)
    implicit none
    integer,intent(in)        :: nBA,n1,n2,n3,nG1,nG2,nR1,offset
    real(realk),intent(in)    :: CalphaD(nBA,n1,n2)
    real(realk),intent(in)    :: CalphaR(nBA,n1,n3)
    real(realk),intent(in)    :: CalphaG(nBA,n1,n2)
    real(realk),intent(inout) :: EJK
    !local variables
    real(realk)               :: EJ, EK, ED
    integer :: c,n,i,j,alpha,beta,alpha1,alpha2,beta1,beta2
    real(realk) :: tmpR,tmpRJ1,tmpRJ2
    real(realk) :: tmpG,tmpGJ1,tmpGJ2
    ED = 0.0E0_realk
    EK = 0.0E0_realk
    EJ = 0.0E0_realk
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(c,n,i,j,alpha,beta,alpha1,alpha2,&
    !$OMP beta1,beta2,tmpR,tmpRJ1,tmpRJ2,tmpG,tmpGJ1,tmpGJ2) SHARED(nBA,n1,n2,n3,CalphaR,&
    !$OMP CalphaG,CalphaD) REDUCTION(+:ED,EJ,EK)
    DO c=1,n3 !ncabsMO
       DO n=1,n2 !noccfull
          DO j=1,n1 !noccEOS
             !Diagonal
             tmpR = 0.0E0_realk
             tmpG = 0.0E0_realk
             DO alpha = 1,nBA
                tmpR = tmpR + CalphaR(alpha,j,c)*CalphaD(alpha,j,n)
                tmpG = tmpG + CalphaR(alpha,j,c)*CalphaG(alpha,j,n)
             ENDDO
             ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
             !Non Diagonal
             DO i=j+1,n1
                tmpRJ1 = 0.0E0_realk
                tmpRJ2 = 0.0E0_realk
                tmpGJ1 = 0.0E0_realk
                tmpGJ2 = 0.0E0_realk
                DO alpha = 1, nBA
                   tmpRJ1 = tmpRJ1 + CalphaR(alpha,i,c)*CalphaD(alpha,j,n) 
                   tmpRJ2 = tmpRJ2 + CalphaR(alpha,j,c)*CalphaD(alpha,i,n)
                   tmpGJ1 = tmpGJ1 + CalphaR(alpha,i,c)*CalphaG(alpha,j,n)
                   tmpGJ2 = tmpGJ2 + CalphaR(alpha,j,c)*CalphaG(alpha,i,n)
                ENDDO
                EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)             
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    EJK = 1.0E0_realk*(ED*0.5E0_realk + 7.0E0_realk/16.0E0_realk*EJ + 1.0E0_realk/16.0E0_realk*EK) 
    
  end subroutine F12RIB7

  !CalphaD(alpha1,j,n)*CalphaR(alpha1,i,c) * CalphaR(beta1,i,c)*CalphaG(beta1,j,n)
  subroutine F12RIB7MPI(nBA,CalphaD,n1,n2,CalphaR,nR1,n3,CalphaG,nG1,nG2,&
       & EJK,offset,NBA2,CalphaDMPI,CalphaRMPI)
    implicit none
    integer,intent(in)        :: nBA,n1,n2,n3,NBA2,nR1,nG1,nG2,offset
    real(realk),intent(in)    :: CalphaDMPI(nBA2,n1,n2)
    real(realk),intent(in)    :: CalphaRMPI(nBA2,n1,n3)
    real(realk),intent(in)    :: CalphaD(nBA,n1,n2)
    real(realk),intent(in)    :: CalphaR(nBA,n1,n3)
    real(realk),intent(in)    :: CalphaG(nBA,n1,n2)
    real(realk),intent(inout) :: EJK
    !local variables
    real(realk)               :: EJ, EK, ED
    integer :: c,n,i,j,alpha
    real(realk) :: tmpR,tmpRJ1,tmpRJ2
    real(realk) :: tmpG,tmpGJ1,tmpGJ2
    ED = 0.0E0_realk
    EK = 0.0E0_realk
    EJ = 0.0E0_realk
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(c,n,i,j,alpha,&
    !$OMP tmpR,tmpRJ1,tmpRJ2,tmpG,tmpGJ1,tmpGJ2) SHARED(nBA,n1,n2,n3,CalphaR,CalphaG,&
    !$OMP CalphaD,CalphaRMPI,CalphaDMPI,NBA2) REDUCTION(+:ED,EJ,EK)
    DO c=1,n3 !ncabsMO
       DO n=1,n2 !noccfull
          DO j=1,n1 !noccEOS
             !Diagonal
             tmpR = 0.0E0_realk
             tmpG = 0.0E0_realk
             DO alpha = 1,nBA2
                tmpR = tmpR + CalphaRMPI(alpha,j,c)*CalphaDMPI(alpha,j,n)
             ENDDO
             DO alpha = 1,nBA
                tmpG = tmpG + CalphaR(alpha,j,c)*CalphaG(alpha,j,n)
             ENDDO
             ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
             !Non Diagonal
             DO i=j+1,n1
                tmpRJ1 = 0.0E0_realk
                tmpRJ2 = 0.0E0_realk
                DO alpha = 1, nBA2
                   tmpRJ1 = tmpRJ1 + CalphaRMPI(alpha,i,c)*CalphaDMPI(alpha,j,n) 
                   tmpRJ2 = tmpRJ2 + CalphaRMPI(alpha,j,c)*CalphaDMPI(alpha,i,n)
                ENDDO
                tmpGJ1 = 0.0E0_realk
                tmpGJ2 = 0.0E0_realk
                DO alpha = 1, nBA
                   tmpGJ1 = tmpGJ1 + CalphaR(alpha,i,c)*CalphaG(alpha,j,n)
                   tmpGJ2 = tmpGJ2 + CalphaR(alpha,j,c)*CalphaG(alpha,i,n)
                ENDDO
                EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    EJK = 1.0E0_realk*(ED*0.5E0_realk + 7.0E0_realk/16.0E0_realk*EJ + 1.0E0_realk/16.0E0_realk*EK) 

  end subroutine F12RIB7MPI

  !CalphaG(alpha1,m,j)*CalphaR(alpha1,i,p) * CalphaR(beta1,i,p)*CalphaD(beta1,j,m)
  subroutine F12RIB8(nBA,CalphaG,nocv,n1,CalphaR,nR1,n2,CalphaD,nD1,noccfull,EJK,offset)
    implicit none
    integer,intent(in)        :: nBA,n1,n2,nocv,noccfull,nR1,nD1,offset
    real(realk),intent(in)    :: CalphaG(nBA,nocv,n1)
    real(realk),intent(in)    :: CalphaR(nBA,n1,n2) 
    real(realk),intent(in)    :: CalphaD(nBA,n1,noccfull)
    real(realk),intent(inout) :: EJK
    !local variables
    real(realk)               :: EJ, EK, ED
    integer :: p,m,i,j,alpha
    real(realk) :: tmpR,tmpRJ1,tmpRJ2
    real(realk) :: tmpG,tmpGJ1,tmpGJ2
    ED = 0.0E0_realk
    EK = 0.0E0_realk
    EJ = 0.0E0_realk
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(p,m,i,j,alpha,&
    !$OMP tmpR,tmpRJ1,tmpRJ2,tmpG,tmpGJ1,tmpGJ2) SHARED(nBA,n1,n2,&
    !$OMP nocv,noccfull,CalphaR,CalphaG,CalphaD) REDUCTION(+:ED,EJ,EK)
    DO p=1,n2 !ncabs
       DO m=1,noccfull !noccfull
          DO j=1,n1 !nocc
             !Diagonal
             tmpR = 0.0E0_realk
             tmpG = 0.0E0_realk
             DO alpha = 1,nBA
                tmpR = tmpR + CalphaR(alpha,j,p)*CalphaG(alpha,m,j)
                tmpG = tmpG + CalphaR(alpha,j,p)*CalphaD(alpha,j,m)
             ENDDO
             ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
             !Non Diagonal
             DO i=j+1,n1
                tmpRJ1 = 0.0E0_realk
                tmpRJ2 = 0.0E0_realk
                tmpGJ1 = 0.0E0_realk
                tmpGJ2 = 0.0E0_realk
                DO alpha = 1, nBA
                   tmpRJ1 = tmpRJ1 + CalphaR(alpha,i,p)*CalphaG(alpha,m,j) 
                   tmpRJ2 = tmpRJ2 + CalphaR(alpha,j,p)*CalphaG(alpha,m,i)
                   tmpGJ1 = tmpGJ1 + CalphaR(alpha,i,p)*CalphaD(alpha,j,m)
                   tmpGJ2 = tmpGJ2 + CalphaR(alpha,j,p)*CalphaD(alpha,i,m)
                ENDDO
                EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    EJK = -2.0E0_realk*(ED*0.5E0_realk + 7.0_realk/16.0_realk*EJ + 1.0_realk/16.0_realk*EK)
  end subroutine F12RIB8

  !CalphaG(alpha1,m,j)*CalphaR(alpha1,i,p) * CalphaR(beta1,i,p)*CalphaD(beta1,j,m)
  subroutine F12RIB8MPI(nBA,CalphaG,nocv,n1,CalphaR,nR1,n2,CalphaD,nD1,noccfull,EJK,offset,NBA2,CalphaGMPI,CalphaRMPI)
    implicit none
    integer,intent(in)        :: nBA,n1,n2,nocv,noccfull,NBA2,nR1,nD1,offset
    real(realk),intent(in)    :: CalphaG(nBA,nocv,n1)
    real(realk),intent(in)    :: CalphaR(nBA,n1,n2) 
    real(realk),intent(in)    :: CalphaD(nBA,n1,noccfull)
    real(realk),intent(in)    :: CalphaGMPI(nBA2,nocv,n1) 
    real(realk),intent(in)    :: CalphaRMPI(nBA2,n1,n2)
    real(realk),intent(inout) :: EJK
    !local variables
    real(realk)               :: EJ, EK, ED
    integer :: p,m,i,j,alpha,beta
    real(realk) :: tmpR,tmpRJ1,tmpRJ2
    real(realk) :: tmpG,tmpGJ1,tmpGJ2
    ED = 0.0E0_realk
    EK = 0.0E0_realk
    EJ = 0.0E0_realk
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(p,m,i,j,alpha,beta,&
    !$OMP tmpR,tmpRJ1,tmpRJ2,tmpG,tmpGJ1,tmpGJ2) SHARED(nBA,n1,n2,&
    !$OMP nocv,noccfull,CalphaR,CalphaG,CalphaD,CalphaRMPI,CalphaGMPI,&
    !$OMP NBA2) REDUCTION(+:ED,EJ,EK)
    DO p=1,n2 !ncabs
       DO m=1,noccfull !noccfull
          DO j=1,n1 !nocc
             !Diagonal
             tmpR = 0.0E0_realk
             tmpG = 0.0E0_realk
             DO alpha = 1,nBA2
                tmpR = tmpR + CalphaRMPI(alpha,j,p)*CalphaGMPI(alpha,m,j)
             ENDDO
             DO beta = 1,nBA
                tmpG = tmpG + CalphaR(beta,j,p)*CalphaD(beta,j,m)
             ENDDO
             ED = ED + tmpR*tmpG !We have a factor 2 which is integrated 
             !Non Diagonal
             DO i=j+1,n1
                tmpRJ1 = 0.0E0_realk
                tmpRJ2 = 0.0E0_realk
                DO alpha = 1, nBA2
                   tmpRJ1 = tmpRJ1 + CalphaRMPI(alpha,i,p)*CalphaGMPI(alpha,m,j) 
                   tmpRJ2 = tmpRJ2 + CalphaRMPI(alpha,j,p)*CalphaGMPI(alpha,m,i)
                ENDDO
                tmpGJ1 = 0.0E0_realk
                tmpGJ2 = 0.0E0_realk
                DO beta = 1, nBA
                   tmpGJ1 = tmpGJ1 + CalphaR(beta,i,p)*CalphaD(beta,j,m)
                   tmpGJ2 = tmpGJ2 + CalphaR(beta,j,p)*CalphaD(beta,i,m)
                ENDDO
                EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    EJK = -2.0E0_realk*(ED*0.5E0_realk + 7.0_realk/16.0_realk*EJ + 1.0_realk/16.0_realk*EK)
  end subroutine F12RIB8MPI
  
end module f12ri_util_module
