!> @file
!> Full calculation
!> This file is mainly a playground for new developments, not intended to be included in a release.

module decf12_routines_module


  use fundamental
  use precision
  use typedeftype!, only:lsitem
  use typedef
  use dec_typedef_module
  use matrix_module
  use matrix_operations
  use memory_handling
  use IchorErimoduleHost

  ! DEC DEPENDENCIES (within deccc directory)   
  ! *****************************************
  use CABS_operations
  use full_f12contractions
  use ccintegrals  
 
  public :: ContractOne4CenterF12IntegralsRIV1_dec, ContractOne4CenterF12IntegralsRIX1_dec, &
       & ContractOne4CenterF12IntegralsRobustRIB1_dec, & 
       & ContractOne4CenterF12IntegralsRIB23_dec, ContractTwo4CenterF12IntegralsRIV2_dec, &
       & ContractTwo4CenterF12IntegralsRIX2_dec, ContractTwo4CenterF12IntegralsRIV34_dec, &
       & ContractTwo4CenterF12IntegralsRIV5_dec, ContractTwo4CenterF12IntegralsRIX34_dec, &
       & ContractTwo4CenterF12IntegralsRIB4_dec, ContractTwo4CenterF12IntegralsRIB5_dec, &
       & ContractTwo4CenterF12IntegralsRIB6_dec, ContractTwo4CenterF12IntegralsRIB7_dec, &
       & ContractTwo4CenterF12IntegralsRIB8_dec, ContractTwo4CenterF12IntegralsRIB9_dec
  private

contains

!==========================================================
!===           Subroutines for DEC-RIMP2-F12            ===
!==========================================================


!Perform on GPU 
subroutine ContractOne4CenterF12IntegralsRIV1_dec(nBA,n,Calpha,EJ,EK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n
   real(realk),intent(in)    :: Calpha(nBA,n,n)
   real(realk),intent(inout) :: EJ,EK
   !local variables
   integer :: I,ALPHA,J
   real(realk) :: TMP,TMPV(n),TMPI
   !Dopair
   logical,intent(in),optional :: dopair_occ_in(n,n)
   logical :: dopair_occ(n,n)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif 
   
   !Exchange Fiijj
   EK = 0.0E0_realk
   EJ = 0.0E0_realk

   !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(I,J,&
   !$OMP ALPHA) SHARED(Calpha,n,nba,dopair_occ) REDUCTION(+:EK,EJ)
   DO I=1,n
      DO J=1,n
         if(dopair_occ(I,J)) then
            DO ALPHA = 1, NBA
               EJ = EJ + Calpha(ALPHA,I,I)*Calpha(ALPHA,J,J)
               EK = EK + Calpha(ALPHA,I,J)*Calpha(ALPHA,J,I)
            ENDDO
         endif
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   !EJK = 7.0/32.0*EJ + 1.0/32.0*EK 
   !print *,"COULOMBX2:  ",  -1.0E0_realk*(5.0E0_realk*0.25E0_realk)*EJ
   !print *,"EXCHANGEX2: ",   1.0E0_realk*(0.250_realk*EK) 
   !print *,"COULOMBX2+EXCHANGEX2:", -1.0E0_realk*((5.0E0_realk*0.25E0_realk)*EJ-EK*0.25E0_realk)   

end subroutine ContractOne4CenterF12IntegralsRIV1_dec

subroutine ContractOne4CenterF12IntegralsRIX1_dec(nBA,n,Calpha,CalphaT,EKJ,dopair_occ_in)
   implicit none                                                                                                             
   integer,intent(in)        :: nBA,n                                                                                        
   real(realk),intent(in)    :: Calpha(nBA,n,n)                                                                              
   real(realk),intent(in)    :: CalphaT(nBA,n,n)                                                                             
   real(realk),intent(inout) :: EKJ
   real(realk)               :: EJ,EK
   !local variables                                                                                                          
   integer :: i,alpha,beta,j                                                                                                 
   real(realk) :: tmpR1,tmpR2,tmpG1,tmpG2,TMPV(n),TMPI                                                                       
   !Dopair                                                                                                                   
   logical,intent(in),optional :: dopair_occ_in(n,n)                                                                         
   logical :: dopair_occ(n,n)                                                                                                
                                                                                                                             
   if(present(dopair_occ_in)) then                                                                                           
      dopair_occ = dopair_occ_in                                                                                             
   else                                                                                                                      
      dopair_occ = .TRUE.                                                                                                    
   endif                                                                                                                     
                                                                                                                             
   EK = 0.0E0_realk                                                                                                          
   EJ = 0.0E0_realk                                                                                                          
   DO i=1,n                                                                                                                  
      DO j=1,n                                                                                                               
         if(dopair_occ(i,j)) then                                                                                            
            tmpR1 = 0.0E0_realk                                                                                              
            tmpR2 = 0.0E0_realk                                                                                              
            DO alpha = 1, nBA                                                                                                
               tmpR1 = tmpR1 + Calpha(ALPHA,i,i)*CalphaT(ALPHA,j,j)                                                          
               tmpR2 = tmpR2 + Calpha(ALPHA,j,j)*CalphaT(ALPHA,i,i)                                                          
            ENDDO                                                                                                            
            tmpG1 = 0.0E0_realk                                                                                              
            tmpG2 = 0.0E0_realk                                                                                              
            DO beta = 1, nBA                                                                                                 
               tmpG1 = tmpG1 + Calpha(BETA,j,i)*CalphaT(BETA,i,j)                                                            
               tmpG2 = tmpG2 + Calpha(BETA,i,j)*CalphaT(BETA,j,i)                                                            
            ENDDO                                                                                                            
            EJ = EJ + tmpR1+tmpR2                                                                                            
            EK = EK + tmpG1+tmpG2                                                                                            
         endif                                                                                                               
      ENDDO                                                                                                                  
   ENDDO        
   EKJ = -1.0*(7.0/32.0*EJ + 1.0/32.0*EK)

end subroutine ContractOne4CenterF12IntegralsRIX1_dec

subroutine ContractOne4CenterF12IntegralsRobustRIB1_dec(nBA,ncore,n,nocvAOS,Rtilde,CalphaR,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n,nocvAOS,ncore
   real(realk),intent(in)    :: Rtilde(nBA,n,n)
   real(realk),intent(in)    :: CalphaR(nBA,n,nocvAOS)
   real(realk),intent(inout) :: EJK
   !local variables
   integer :: I,ALPHA,J
   real(realk) :: TMP,TMP_IJIJ,TMP_JIIJ
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n,n)
   logical :: dopair_occ(n,n)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif
   !print *, "***************CalphaR"
   !do alpha=1,nba
   !  call ls_output(CalphaR(alpha,:,:),1,n,1,n,n,n,1,6)
   !enddo

   !print *, "*************Rtilde"
   !do alpha=1,nba
   !   call ls_output(Rtilde(alpha,:,:),1,n,1,n,n,n,1,6)
   !enddo

   TMP = 0.0E0_realk
   EJK = 0.0E0_realk
   !$OMP PARALLEL DO DEFAULT(none) PRIVATE(I,J,ALPHA,TMP_IJIJ,TMP_JIIJ) SHARED(CalphaR,Rtilde,n,&
   !$OMP nba,dopair_occ,ncore) REDUCTION(+:TMP)
   !Diagonal
   DO J=1,n !noccEOS
      DO I=1,n !noccEOS
         IF(dopair_occ(I,J)) THEN
            TMP_IJIJ = 0.0E0_realk
            TMP_JIIJ = 0.0E0_realk
            DO ALPHA = 1,NBA
               TMP_IJIJ = TMP_IJIJ + CalphaR(ALPHA,I,I+ncore)*Rtilde(ALPHA,J,J) + Rtilde(ALPHA,I,I)*CalphaR(ALPHA,J,J+ncore)
               TMP_JIIJ = TMP_JIIJ + CalphaR(ALPHA,I,J+ncore)*Rtilde(ALPHA,J,I) + Rtilde(ALPHA,I,J)*CalphaR(ALPHA,J,I+ncore)
            ENDDO
            TMP = TMP + 7.0E0_realk * TMP_IJIJ + TMP_JIIJ
         ENDIF
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   EJK = 1.0/32.0*TMP
end subroutine ContractOne4CenterF12IntegralsRobustRIB1_dec

subroutine ContractOne4CenterF12IntegralsRIB23_dec(nBA,n1,n2,CalphaR,CalphaG,hJir, &
                                                   & Coeff,EJK2,EJK3,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2
   real(realk),intent(in)    :: CalphaR(nBA,n1,n2),CalphaG(nBA,n1,n1)
   real(realk),intent(inout) :: EJK2,EJK3
   real(realk),intent(in)    :: Coeff
   real(realk)               :: EJ2,EJ3,EK2,EK3
   real(realk),intent(IN)    :: HJir(n1,n2)
   !local variables
   integer :: c,i,j,alpha,beta
   real(realk) :: tmpR2,tmpR21,tmpR22,tmpR31,tmpR32,tmp,tmpR
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   EJ2 = 0.0E0_realk
   EK2 = 0.0E0_realk
   EJ3 = 0.0E0_realk
   EK3 = 0.0E0_realk 
   DO c=1,n2
      DO j=1,n1
         DO i=1,n1
            IF(dopair_occ(I,J)) THEN
               tmpR21 = 0.0E0_realk
               tmpR22 = 0.0E0_realk
               tmpR31 = 0.0E0_realk
               tmpR32 = 0.0E0_realk
               DO alpha = 1,nBA
                  tmpR21 = tmpR21 + CalphaR(alpha,i,c)*CalphaG(alpha,j,j)
                  tmpR22 = tmpR22 + CalphaR(alpha,j,c)*CalphaG(alpha,i,j)

                  tmpR31 = tmpR31 + CalphaR(alpha,j,c)*CalphaG(alpha,i,i)
                  tmpR32 = tmpR32 + CalphaR(alpha,i,c)*CalphaG(alpha,i,j)
               ENDDO
               EJ2 = EJ2 + tmpR21*hJir(i,c)
               EK2 = EK2 + tmpR22*hJir(i,c)
               EJ3 = EJ3 + tmpR31*hJir(j,c)
               EK3 = EK3 + tmpR32*hJir(j,c)
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   EJK2 = 7.0/32.0*EJ2 + 1.0/32.0*EK2
   EJK3 = 7.0/32.0*EJ3 + 1.0/32.0*EK3 

end subroutine ContractOne4CenterF12IntegralsRIB23_dec

subroutine ContractTwo4CenterF12IntegralsRIV2_dec(nBA,nBA2,n1,n2,CalphaR,CalphaG,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,nBA2,n1,n2
   real(realk),intent(in)    :: CalphaR(nBA,n1,n2)
   real(realk),intent(in)    :: CalphaG(nBA2,n2,n1)
   real(realk),intent(inout) :: EJK
   real(realk)               :: ED, EJ,EK
   !local variables
   integer :: q,p,i,j,alpha,beta
   real(realk) :: tmpR,tmpG1,tmpG2,tmp
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif  
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   tmp = 0.0E0_realk
   ED =  0.0E0_realk
   EJ =  0.0E0_realk
   EK =  0.0E0_realk
   EJK = 0.0E0_realk
   !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(i,j,p,q,tmpR,&
   !$OMP tmpG1,tmpG2) SHARED(CalphaR,CalphaG,n2,n1,&
   !$OMP nba,nba2,dopair_occ) REDUCTION(+:EJ,EK,ED)
   DO q=1,n2 !nocv
      DO p=1,n2 !nocv
         DO j=1,n1 !nocc
            DO i=1,n1 !nocc
               IF(dopair_occ(I,J)) THEN
                  tmpR = 0.0E0_realk
                  DO alpha = 1,NBA
                     tmpR = tmpR + CalphaR(alpha,i,p)*CalphaR(alpha,j,q)
                  ENDDO
                  tmpG1 = 0.0E0_realk
                  tmpG2 = 0.0E0_realk
                  DO beta = 1,NBA2
                     !tmpG1 = tmpG1 + CalphaG(beta,p,i)*CalphaG(beta,q,j)
                     !tmpG2 = tmpG2 + CalphaG(beta,p,j)*CalphaG(beta,q,i)
                     tmpG1 = tmpG1 + CalphaG(beta,p,i)*CalphaG(beta,q,j)
                     tmpG2 = tmpG2 + CalphaG(beta,p,j)*CalphaG(beta,q,i)
                  ENDDO
                  EJ = EJ + tmpR*tmpG1 
                  EK = EK + tmpR*tmpG2
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   EJK = 5.0/4.0*EJ - 1.0/4.0*EK 
end subroutine ContractTwo4CenterF12IntegralsRIV2_dec

subroutine ContractTwo4CenterF12IntegralsRIX2_dec(nBA,nBA2,n1,n2,CalphaC,CalphaCMPI,CalphaG,EJK,dopair_occ_in)
 implicit none
   integer,intent(in)        :: nBA,nBA2,n1,n2
   real(realk),intent(in)    :: CalphaC(nBA,n2,n1)
   real(realk),intent(in)    :: CalphaCMPI(nBA2,n2,n1)
   real(realk),intent(in)    :: CalphaG(nBA,n2,n1)
   real(realk),intent(inout) :: EJK
   real(realk)               :: ED, EJ,EK
   !local variables
   integer :: q,p,i,j,alpha,beta,gamma
   real(realk) :: tmpG1,tmpG2,tmpG3,tmpG4,tmpR1,tmp
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   tmp = 0.0E0_realk
   ED =  0.0E0_realk
   EJ =  0.0E0_realk
   EK =  0.0E0_realk
   EJK = 0.0E0_realk
   DO q=1,n2 !nocvAOS
      DO p=1,n2 !nocvAOS
         DO i=1,n1 !noccEOS
            DO j=1,n1 !noccEOS
               if(dopair_occ(i,j)) then
                  tmpR1 = 0.0E0_realk
                  DO alpha = 1, nBA2
                     tmpR1 = tmpR1 + CalphaCMPI(ALPHA,p,i)*CalphaCMPI(ALPHA,q,j)
                  ENDDO
                  tmpG1 = 0.0E0_realk
                  tmpG2 = 0.0E0_realk
                  tmpG3 = 0.0E0_realk
                  tmpG4 = 0.0E0_realk
                  DO beta = 1, nBA
                     tmpG1 = tmpG1 + CalphaC(BETA,p,i)*CalphaG(BETA,q,j)
                     tmpG2 = tmpG2 + CalphaG(BETA,p,i)*CalphaC(BETA,q,j)
                     tmpG3 = tmpG3 + CalphaG(BETA,p,j)*CalphaC(BETA,q,i)
                     tmpG4 = tmpG4 + CalphaC(BETA,p,j)*CalphaG(BETA,q,i)
                  ENDDO
                  EJ = EJ + tmpR1*(tmpG1 + tmpG2)
                  EK = EK + tmpR1*(tmpG3 + tmpG4)
               endif
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   EJK = 7.0/32.0*EJ + 1.0/32.0*EK 
end subroutine ContractTwo4CenterF12IntegralsRIX2_dec

!NB!NB!
!EJK3 and EJK4 gives the same energy in the non-diagonal form of the subroutine, same for X3 and X4
subroutine ContractTwo4CenterF12IntegralsRIV34_dec(nBA,n1,n2,n3,nocv,&
      & CalphaRcabs,CalphaGcabs,CalphaR,CalphaG,EJK3,EJK4,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,n3,nocv
   real(realk),intent(in)    :: CalphaRcabs(nBA,n1,n3)
   real(realk),intent(in)    :: CalphaGcabs(nBA,n1,n3)
   real(realk),intent(in)    :: CalphaR(nBA,n1,nocv)
   real(realk),intent(in)    :: CalphaG(nBA,nocv,n1)
   real(realk),intent(inout) :: EJK3,EJK4
   real(realk)               :: EJ3, EK3
   !local variables
   integer :: m,c,i,j,alpha,beta
   real(realk) :: tmpR3,tmpG13,tmpG23
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   EJ3 =  0.0E0_realk
   EK3 =  0.0E0_realk
   EJK3 = 0.0E0_realk
   !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(i,j,m,c,tmpR3, &                                                
   !$OMP tmpG13,tmpG23) SHARED(CalphaRcabs,CalphaR,CalphaGcabs,CalphaG,&                                       
   !$OMP n3,n2,n1,nba,dopair_occ) REDUCTION(+:EJ3,EK3)                                                            
   DO c=1,n3                                                                                                                 
      DO m=1,n2 !noccAOStot                                                                                                   
         DO j=1,n1                                                                                                           
            DO i=1,n1                                                                                                      
               IF(dopair_occ(I,J)) THEN                                                                                      
                  tmpR3 = 0.0E0_realk                                                                                        
                  DO alpha = 1,NBA                                                                                           
                     tmpR3 = tmpR3 + CalphaR(alpha,i,m)*CalphaRcabs(alpha,j,c)                                               
                  ENDDO                                                                                                      
                  tmpG13 = 0.0E0_realk                                                                                       
                  tmpG23 = 0.0E0_realk                                                                                       
                  DO beta = 1,NBA                                                                                            
                     tmpG13 = tmpG13 + CalphaG(beta,m,i)*CalphaGcabs(beta,j,c)                                               
                     tmpG23 = tmpG23 + CalphaG(beta,m,j)*CalphaGcabs(beta,i,c)                                               
                  ENDDO                                                                                                      
                  EJ3 = EJ3 + tmpR3*tmpG13                                                                                   
                  EK3 = EK3 + tmpR3*tmpG23                                                                                   
               ENDIF
            ENDDO                                                                                                            
         ENDDO                                                                                                               
      ENDDO                                                                                                                  
   ENDDO                                                                                                                     
   !$OMP END PARALLEL DO
   EJK3 = 1.25*EJ3-0.25E0_realk*EK3
   EJK4 = EJK3
end subroutine ContractTwo4CenterF12IntegralsRIV34_dec

subroutine ContractTwo4CenterF12IntegralsRIV5_dec(nBA,n1,n2,CalphaV,CalphaD,Taibj,EJK,dopair_occ_in)
   implicit none 
   integer,intent(in)        :: nBA,n1,n2
   real(realk),intent(in)    :: CalphaV(nBA,n1,n2),CalphaD(nBA,n1,n2)
   real(realk),intent(inout) :: EJK
   real(realk)               :: ED,EJ,EK,eps
   real(realk),intent(in)    :: Taibj(:,:,:,:)
   !local variables
   integer :: a,b,i,j,alpha,beta,gamma
   real(realk) :: tmp,tmpR,tmpG,tmpT
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   tmp = 0.0E0_realk
   ED =  0.0E0_realk
   EJ =  0.0E0_realk
   EK =  0.0E0_realk
   EJK = 0.0E0_realk
   DO j=1,n1 !nocc
      DO b=1,n2 !nvirt
         DO i=1,n1 !nocc
            DO a=1,n2 !nvirt 
               IF(dopair_occ(I,J)) THEN
                  tmpR = 0.0E0_realk
                  DO alpha = 1,NBA
                     tmpR = tmpR + CalphaV(alpha,i,a)*CalphaD(alpha,j,b) + CalphaV(alpha,j,b)*CalphaD(alpha,i,a)
                  ENDDO
                  tmpG = 0.0E0_realk
                  DO beta = 1,NBA
                     tmpG = tmpG + CalphaV(beta,j,a)*CalphaD(beta,i,b) + CalphaV(beta,i,b)*CalphaD(beta,j,a)
                  ENDDO
                  EJ = EJ + tmpR*Taibj(a,i,b,j)
                  EK = EK + tmpG*Taibj(a,i,b,j)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   EJK = -5.0/4.0*EJ + 1.0/4.0*EK 
   !print *,"COULOMBE V5:  ", 5.0/4.0*EJ
   !print *,"EXCHANGE V5: ",  1.0/4.0*EK
   !print *,"COULOMBE V5 + EXCHANGE: V5", EJK      
end subroutine ContractTwo4CenterF12IntegralsRIV5_dec

subroutine ContractTwo4CenterF12IntegralsRIX34_dec(nBA,n1,n2,n3,&
      & CalphaGcabs,CalphaC,CalphaG,CalphaP,EJK3,EJK4,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,n3
   real(realk),intent(IN)    :: CalphaGcabs(nBA,n1,n3)
   real(realk),intent(IN)    :: CalphaC(nBA,n2,n1)
   real(realk),intent(IN)    :: CalphaG(nBA,n2,n1)
   real(realk),intent(IN)    :: CalphaP(nBA,n3,n1)
   real(realk),intent(inout) :: EJK3,EJK4
   real(realk)               :: EJ3, EJ4, EK3, EK4
   !local variables
   integer :: m,c,i,j,alpha,beta
   real(realk) :: tmpR3,tmpG31,tmpG32,tmpG33,tmpG34
   real(realk) :: tmpR4,tmpG41,tmpG42,tmpG43,tmpG44
   real(realk) :: tmp
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif
   !Exchange Ripjq*Gjpiq Scaling(N*N*O*O*Naux)
   EJ3 =  0.0E0_realk
   EK3 =  0.0E0_realk
   EJK3 = 0.0E0_realk
   EJ4 =  0.0E0_realk
   EK4 =  0.0E0_realk
   EJK4 = 0.0E0_realk
   !!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(i,j,m,c,tmpR3,tmpR4, &
   !!$OMP tmpG13,tmpG23,tmp) SHARED(CalphaC,CalphaP,CalphaG,CalphaGcabs,n3,n2,n1,dopair_occ) &
   !!$OMP nba, Fii) REDUCTION(+:EJ3,EK3,EJ4,EK4)
   DO c=1,n3 !ncabsMO
      DO m=1,n2 !noccAOStot
         DO j=1,n1 !noccEOS
            DO i=1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpR3 = 0.0E0_realk
                  !tmpR4 = 0.0E0_realk
                  DO alpha = 1,NBA
                     tmpR3 = tmpR3 + CalphaC(alpha,m,i)*CalphaGcabs(alpha,j,c)
                  ENDDO
                  tmpG31 = 0.0E0_realk
                  tmpG32 = 0.0E0_realk
                  tmpG33 = 0.0E0_realk
                  tmpG34 = 0.0E0_realk
                  DO beta = 1,NBA
                     tmpG31 = tmpG31 + CalphaC(beta,m,i)*CalphaP(beta,c,j)
                     tmpG32 = tmpG32 + CalphaG(beta,m,i)*CalphaGcabs(beta,j,c)
                     tmpG33 = tmpG33 + CalphaG(beta,m,j)*CalphaGcabs(beta,i,c)
                     tmpG34 = tmpG34 + CalphaC(beta,m,j)*CalphaP(beta,c,i)
                  ENDDO
                  EJ3 = EJ3 + tmpR3*(tmpG31 + tmpG32)
                  EK3 = EK3 + tmpR3*(tmpG33 + tmpG34)
                  !EJ4 = EJ4 + tmpR4*(tmpG33 + tmpG34)
                  !EK4 = EK4 + tmpR4*(tmpG31 + tmpG32)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !!$OMP END PARALLEL DO
   EJK3 = 7.0/32.0_realk*EJ3+1.0_realk/32.0*EK3
   !EJK4 = 7.0/32.0_realk*EJ4+1.0_realk/32.0*EK4
   EJK4 = EJK3
  !print *,"COULOMBX2:  ", 7.0/32.0*EJ3
  !print *,"EXCHANGEX2: ", 1.0/32.0*EK3
  !print *,"COULOMBX2+EXCHANGEX2:", 7.0/32.0*EJ3 + 1.0/32.0*EK3      
end subroutine ContractTwo4CenterF12IntegralsRIX34_dec

subroutine ContractTwo4CenterF12IntegralsRIB4_dec(nBA,n1,n2,CalphaG,CalphaD,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2
   real(realk),intent(in)    :: CalphaG(nBA,n1,n2)
   real(realk),intent(in)    :: CalphaD(nBa,n1,n2)
   real(realk),intent(inout) :: EJK
   real(realk)               :: EJ, EK, ED
   !local variables
   integer :: p,q,r,i,j,alpha,beta,alpha1,beta1,alpha2,beta2,alpha3,beta3,alpha4,beta4
   real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
   real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   ED = 0.0E0_realk
   EJ = 0.0E0_realk
   EK = 0.0E0_realk
   DO q=1,n2 !ncabsAO
      DO r=1,n2 !ncabsAO
         DO j=1,n1 !nocc
            DO i=1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpRJ1 = 0.0E0_realk
                  tmpGJ1 = 0.0E0_realk
                  DO alpha1 = 1, nBA
                     tmpRJ1 = tmpRJ1 + CalphaG(alpha1,i,r)*CalphaD(alpha1,j,q) 
                     tmpGJ1 = tmpGJ1 + CalphaG(alpha1,i,r)*CalphaG(alpha1,j,q)
                  ENDDO
                  tmpRJ2 = 0.0E0_realk
                  tmpGJ2 = 0.0E0_realk
                  DO alpha2 = 1, nBA
                     tmpRJ2 = tmpRJ2 + CalphaG(alpha2,j,r)*CalphaD(alpha2,i,q)
                     tmpGJ2 = tmpGJ2 + CalphaG(alpha2,j,r)*CalphaG(alpha2,i,q)
                  ENDDO
                  EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                  EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO   
   EJK = 7.0/32.0*EJ + 1.0/32.0*EK 
end subroutine ContractTwo4CenterF12IntegralsRIB4_dec                

subroutine ContractTwo4CenterF12IntegralsRIB5_dec(nBA,n1,n2,n3,nocv,CalphaGcabs,CalphaG,CalphaD,EJK,dopair_occ_in)
   implicit none                                                                                                             
   integer,intent(in)        :: nBA,n1,n2,n3,nocv                                                                          
   !real(realk),intent(in)    :: CalphaG(nBA,nocv,n1)                                                                       
   real(realk),intent(in)    :: CalphaG(nBA,n3,n1)
   real(realk),intent(in)    :: CalphaGcabs(nBA,n1,n2)                                                                       
   real(realk),intent(in)    :: CalphaD(nBA,n1,n2)                                                                           
   real(realk),intent(inout) :: EJK
   real(realk)               :: EJ, EK
   !local variables                                                                                                          
   integer :: r,m,i,j,alpha,beta,alpha1,beta1,alpha2,beta2                                                                 
   real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2                                                                           
   real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2                                                                           
   !Dopair                                                                                                                   
   logical,intent(in),optional :: dopair_occ_in(n1,n1)                                                                       
   logical :: dopair_occ(n1,n1)                                                                                              
   if(present(dopair_occ_in)) then                                                                                           
      dopair_occ = dopair_occ_in                                                                                             
   else                                                                                                                      
      dopair_occ = .TRUE.                                                                                                    
   endif                                                                                                                     
   EJ = 0.0E0_realk                                                                                                          
   EK = 0.0E0_realk                                                                                                          
   DO r=1,n2 !ncabsAO                                                                                                        
      DO m=1,n3 !noccAOStot
         DO j=1,n1 !noccEOS                                                                                                     
            DO i=1,n1 ! !noccEOS
               IF(dopair_occ(I,J)) THEN                                                                                      
                  tmpRJ1 = 0.0E0_realk                                                                                       
                  tmpGJ1 = 0.0E0_realk                                                                                       
                  DO alpha1 = 1, nBA                                                                                         
                     tmpRJ1 = tmpRJ1 + CalphaD(alpha1,i,r)*CalphaG(alpha1,m,j)                                               
                     tmpGJ1 = tmpGJ1 + CalphaGcabs(alpha1,i,r)*CalphaG(alpha1,m,j)                                           
                  ENDDO                                                                                                      
                  tmpRJ2 = 0.0E0_realk                                                                                       
                  tmpGJ2 = 0.0E0_realk                                                                                       
                  DO alpha2 = 1, nBA                                                                                         
                     tmpRJ2 = tmpRJ2 + CalphaD(alpha2,j,r)*CalphaG(alpha2,m,i)                                               
                     tmpGJ2 = tmpGJ2 + CalphaGcabs(alpha2,j,r)*CalphaG(alpha2,m,i)                                           
                  ENDDO                                                                                                      
                  EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                  EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   EJK = -1.0E0_realk*(7.0/32.0*EJ + 1.0/32.0*EK) 
end subroutine ContractTwo4CenterF12IntegralsRIB5_dec                

subroutine ContractTwo4CenterF12IntegralsRIB6_dec(nBA,n1,n2,n3,noccAOS,CalphaG,CalphaD,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,n3,noccAOS
   real(realk),intent(in)    :: CalphaG(nBA,n3,n1)
   real(realk),intent(in)    :: CalphaD(nBA,n1,n3)
   real(realk),intent(inout) :: EJK
   real(realk)               :: EJ, EK
   !local variables
   integer :: q,a,i,j,alpha,beta,alpha1,beta1,alpha2,beta2
   real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
   real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif
   EK = 0.0E0_realk
   EJ = 0.0E0_realk
   DO q=1,n3 !nocvAOS
      DO a=1+noccAOS,n3 !nvirtAOS
         DO j=1,n1 !noccEOS
            DO i=1,n1 !noccEOS
               IF(dopair_occ(I,J)) THEN
                  tmpRJ1 = 0.0E0_realk
                  tmpGJ1 = 0.0E0_realk
                  DO alpha1 = 1, nBA
                     !tmpRJ1 = tmpRJ1 + CalphaD(alpha1,i,q)*CalphaC(alpha1,j,a) 
                     !tmpGJ1 = tmpGJ1 + CalphaG(alpha1,q,i)*CalphaC(alpha1,j,a)
                     
                     tmpRJ1 = tmpRJ1 + CalphaD(alpha1,i,q)*CalphaG(alpha1,a,j) 
                     tmpGJ1 = tmpGJ1 + CalphaG(alpha1,q,i)*CalphaG(alpha1,a,j)
                  ENDDO
                  tmpRJ2 = 0.0E0_realk
                  tmpGJ2 = 0.0E0_realk
                  DO alpha2 = 1, nBA
                     !tmpRJ2 = tmpRJ2 + CalphaD(alpha2,j,q)*CalphaC(alpha2,i,a)
                     !tmpGJ2 = tmpGJ2 + CalphaG(alpha2,q,j)*CalphaC(alpha2,i,a)
                     
                     tmpRJ2 = tmpRJ2 + CalphaD(alpha2,j,q)*CalphaG(alpha2,a,i)
                     tmpGJ2 = tmpGJ2 + CalphaG(alpha2,q,j)*CalphaG(alpha2,a,i)
                  ENDDO
                  EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                  EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)           
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO   
   EJK = -1.0E0_realk*(7.0/32.0*EJ + 1.0/32.0*EK)
   !print *,"COULOMBB6:  ", 7.0/32.0*EJ
   !print *,"EXCHANGEB6: ", 1.0/32.0*EK
   !print *,"COULOMBB6+EXCHANGEB6:", 7.0/32.0*EJ + 1.0/32.0*EK     
end subroutine ContractTwo4CenterF12IntegralsRIB6_dec

subroutine ContractTwo4CenterF12IntegralsRIB7_dec(nBA,n1,n2,n3,CalphaR,CalphaG,CalphaD,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,n3
   real(realk),intent(in)    :: CalphaG(nBA,n1,n2)
   real(realk),intent(in)    :: CalphaR(nBA,n1,n3)
   real(realk),intent(in)    :: CalphaD(nBA,n1,n2)
   real(realk),intent(inout) :: EJK
   real(realk)               :: EJ, EK
   !local variables
   integer :: c,n,i,j,alpha,beta,alpha1,beta1,alpha2,beta2
   real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
   real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif  
   EK = 0.0E0_realk
   EJ = 0.0E0_realk
   DO c=1,n3 !ncabsMO
      DO n=1,n2 !noccfull
         DO j=1,n1 !noccEOS
            DO i=1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpRJ1 = 0.0E0_realk
                  tmpRJ2 = 0.0E0_realk
                  DO alpha = 1, nBA
                     tmpRJ1 = tmpRJ1 + CalphaR(alpha,i,c)*CalphaD(alpha,j,n) 
                     tmpRJ2 = tmpRJ2 + CalphaR(alpha,j,c)*CalphaD(alpha,i,n)
                  ENDDO
                  tmpGJ1 = 0.0E0_realk
                  tmpGJ2 = 0.0E0_realk
                  DO beta = 1, nBA
                     tmpGJ1 = tmpGJ1 + CalphaR(beta,i,c)*CalphaG(beta,j,n)
                     tmpGJ2 = tmpGJ2 + CalphaR(beta,j,c)*CalphaG(beta,i,n)
                  ENDDO
                  EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                  EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO   
   EJK = 7.0/32.0*EJ + 1.0/32.0*EK

end subroutine ContractTwo4CenterF12IntegralsRIB7_dec

subroutine ContractTwo4CenterF12IntegralsRIB8_dec(nBA,n1,n2,nocv,noccAOStot,CalphaR,CalphaG,CalphaD,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,nocv,noccAOStot
   real(realk),intent(in)    :: CalphaG(nBA,nocv,n1)
   real(realk),intent(in)    :: CalphaR(nBA,n1,n2) !CalphaGcabsMO
   real(realk),intent(in)    :: CalphaD(nBA,n1,noccAOStot)
   real(realk),intent(inout) :: EJK
   real(realk)               :: EJ, EK
   !local variables
   integer :: r,m,i,j,alpha,beta,alpha1,beta1,alpha2,beta2
   real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
   real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif   
   EK = 0.0E0_realk
   EJ = 0.0E0_realk
   !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(r,m,i,j,alpha,beta,alpha1,&
   !$OMP beta1,alpha2,beta2,tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2,tmpG,tmpGJ1,tmpGJ2,&
   !$OMP tmpGK1,tmpGK2) SHARED(CalphaR,CalphaG,CalphaD,n2,n1,noccAOStot,nba,dopair_occ) REDUCTION(+:EJ,EK)
   DO r=1,n2 !ncabsAO
      DO m=1,noccAOStot !noccAOS
         DO j=1,n1 !noccEOS
            DO i=1,n1 !noccEOS
               IF(dopair_occ(I,J)) THEN
                  tmpRJ1 = 0.0E0_realk
                  DO alpha1 = 1, nBA
                     tmpRJ1 = tmpRJ1 + CalphaR(alpha1,i,r)*CalphaG(alpha1,m,j) 
                  ENDDO
                  tmpRJ2 = 0.0E0_realk
                  DO alpha2 = 1, nBA
                     tmpRJ2 = tmpRJ2 + CalphaR(alpha2,j,r)*CalphaG(alpha2,m,i)
                  ENDDO
                  tmpGJ1 = 0.0E0_realk
                  DO beta1 = 1, nBA
                     tmpGJ1 = tmpGJ1 + CalphaR(beta1,i,r)*CalphaD(beta1,j,m)
                  ENDDO
                  tmpGJ2 = 0.0E0_realk
                  DO beta2 = 1, nBA
                     tmpGJ2 = tmpGJ2 + CalphaR(beta2,j,r)*CalphaD(beta2,i,m)
                  ENDDO
                  EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                  EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   EJK = -2.0E0_realk*(7.0_realk/32.0_realk*EJ + 1.0_realk/32.0_realk*EK)
end subroutine ContractTwo4CenterF12IntegralsRIB8_dec

subroutine ContractTwo4CenterF12IntegralsRIB9_dec(nBA,n1,n2,n3,noccAOS,CalphaG,CalphaD,EJK,dopair_occ_in)
   implicit none
   integer,intent(in)        :: nBA,n1,n2,n3,noccAOS
   real(realk),intent(in)    :: CalphaG(nBA,n3,n1)
   real(realk),intent(in)    :: CalphaD(nBA,n1,n3)
   real(realk),intent(inout) :: EJK
   real(realk)               :: EJ, EK
   !local variables
   integer :: c,a,i,j,alpha,beta,alpha1,beta1,alpha2,beta2
   real(realk) :: tmpR,tmpRJ1,tmpRJ2,tmpRK1,tmpRK2
   real(realk) :: tmpG,tmpGJ1,tmpGJ2,tmpGK1,tmpGK2
   !Dopair                                                                          
   logical,intent(in),optional :: dopair_occ_in(n1,n1)
   logical :: dopair_occ(n1,n1)
   if(present(dopair_occ_in)) then
      dopair_occ = dopair_occ_in
   else
      dopair_occ = .TRUE.
   endif
   EK = 0.0E0_realk
   EJ = 0.0E0_realk
   DO c=1,n3 !ncabsAO
      DO a=noccAOS+1,n3 !nvirtAOS
         DO j=1,n1 !nocc
            DO i=1,n1
               IF(dopair_occ(I,J)) THEN
                  tmpRJ1 = 0.0E0_realk
                  DO alpha1 = 1, nBA
                     tmpRJ1 = tmpRJ1 + CalphaG(alpha1,c,i)*CalphaG(alpha1,a,j) 
                  ENDDO      
                  tmpRJ2 = 0.0E0_realk
                  DO alpha2 = 1, nBA
                     tmpRJ2 = tmpRJ2 + CalphaG(alpha2,c,j)*CalphaG(alpha2,a,i)
                  ENDDO
                  tmpGJ1 = 0.0E0_realk
                  DO beta1 = 1, nBA
                     tmpGJ1 = tmpGJ1 + CalphaD(beta1,i,c)*CalphaG(beta1,a,j)
                  ENDDO
                  tmpGJ2 = 0.0E0_realk
                  DO beta2 = 1, nBA
                     tmpGJ2 = tmpGJ2 + CalphaD(beta2,j,c)*CalphaG(beta2,a,i)
                  ENDDO
                  EJ = EJ + (tmpRJ1*tmpGJ1 + tmpRJ2*tmpGJ2)
                  EK = EK + (tmpRJ2*tmpGJ1 + tmpRJ1*tmpGJ2)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   EJK = -2.0E0_realk*(7.0_realk/32.0_realk*EJ + 1.0_realk/32.0_realk*EK)
end subroutine ContractTwo4CenterF12IntegralsRIB9_dec

end module decf12_routines_module
