SUBROUTINE LSFORT_WRT(lupri,str, len)
  INTEGER   lupri,len
  CHARACTER str(len)
  WRITE(LUPRI,*) str
END SUBROUTINE LSFORT_WRT


FUNCTION DISERR90(RD,L)
use precision
IMPLICIT NONE                                                                  
!     Provide grid spacing h for given angular momentum L           
!     and discretization error RD                                   
!                                                                   
!     Based on eqs. (17) and (18) of                                
!       R. Lindh, P.-Aa. Malmqvist and L. Gagliardi                 
!       "Molecular integrals by numerical quadrature",              
!       Theor. Chem. Acc. 106 (2001) 178-187                        
!                                                                   
!     The array CF(4,L) contains coefficients of a 3rd order        
!     polynomial fit to provide start values for the                
!     determination of H by a Newton-Raphson search.                
!                                                                   
!     Written by T. Saue July 2002                                  
!     Transformed to fortran 90 by T. kjaergaard                                                              

!      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
REAL(REALK)             :: RD
INTEGER                 :: L

REAL(realk), PARAMETER  :: ACC=1.0E-5_realk
REAL(realk), PARAMETER  :: D0=0.0E0_realk,D1=1.0E0_realk,D2=2.0E0_realk,TI=1.0E1_realk
REAL(realk), PARAMETER  :: PI = 3.141592653589793E00_realk
!REAL, PARAMETER  :: ACC=1.0E-5_realk
!REAL, PARAMETER  :: D0=0.0E0_realk,D1=1.0E0_realk,D2=2.0E0_realk,TI=1.0E1_realk
!REAL, PARAMETER  :: PI = 3.141592653589793E00_realk
INTEGER, PARAMETER      :: MXIT=20
REAL(realk)             :: CF(4,0:4)
REAL(realk)             :: FAC,RDOG,DISERR90,HTLOG,PIH,PIHL
REAL(realk)             :: PIEX,U0,U1,F0,F1,DX,RDLOG,POLVAL90
!REAL             :: CF(4,0:4)
!REAL             :: FAC,RDOG,DISERR90,HTLOG,PIH,PIHL
!REAL             :: PIEX,U0,U1,F0,F1,DX,RD,RDLOG,POLVAL90
INTEGER           :: IFAC,I,LM,IT
CF(1,0) = 0.91570E0_realk
CF(2,0) = 0.78806E-1_realk
CF(3,0) = 0.28056E-2_realk
CF(4,0) = 3.4197E-05_realk
CF(1,1) = 0.74912E0_realk
CF(2,1) = 0.61502E-1_realk
CF(3,1) = 0.21558E-2_realk
CF(4,1) = 2.6100E-05_realk
CF(1,2) = 0.65449E0_realk
CF(2,2) = 0.52322E-1_realk
CF(3,2) = 0.18217E-2_realk
CF(4,2) = 2.2004E-05_realk
CF(1,3) = 0.59321E0_realk
CF(2,3) = 0.46769E-1_realk
CF(3,3) = 0.16261E-2_realk
CF(4,3) = 1.9649E-05_realk
CF(1,4) = 0.55125E0_realk
CF(2,4) = 0.43269E-1_realk
CF(3,4) = 0.15084E-2_realk
CF(4,4) = 1.8270E-05_realk


!DATA CF/0.91570E0_realk,0.78806E-1_realk,0.28056E-2_realk,3.4197E-05_realk,
!     &        0.74912E0_realk,0.61502E-1_realk,0.21558E-2_realk,2.6100E-05_realk,
!     &        0.65449E0_realk,0.52322E-1_realk,0.18217E-2_realk,2.2004E-05_realk,
!     &        0.59321E0_realk,0.46769E-1_realk,0.16261E-2_realk,1.9649E-05_realk,
!     &        0.55125E0_realk,0.43269E-1_realk,0.15084E-2_realk,1.8270E-05_realk/
!
!     Initialization
!
FAC  = SQRT(D2)*D2*D2
IFAC = 1
DO I = 1,L
   FAC   = FAC*D2
   IFAC  = IFAC*(2*I+1)
ENDDO
FAC = FAC/IFAC
LM = MIN(L,4)
RDLOG = LOG(RD)
DISERR90 = POLVAL90(3,CF(1,LM),RDLOG)
HTLOG = LOG(DISERR90)
!     Newton-Raphson search
DO IT = 1,MXIT
   PIH  = PI/DISERR90
   PIHL = PIH
   PIEX = PI*PIH/D2
   DO I = 1,L
      PIHL = PIHL*PIH
   ENDDO
   U0   = FAC*PIHL*EXP(-PIEX)
   U1   = U0*((PIEX/DISERR90)-(L+1)/PIH)
   F0   = LOG(U0)-RDLOG
   F1   = DISERR90*U1/U0
   DX = F0/F1
   HTLOG = HTLOG - DX
   DISERR90 = EXP(HTLOG)
   IF(ABS(DX).LT.ACC) RETURN
ENDDO
STOP 'Error in DISERR90'

END FUNCTION DISERR90

FUNCTION OUTERR90(AL,L,RD)
use precision
IMPLICIT NONE                                                                      
!     Provide outer grid point for given angular momentum L           
!     outer exponent AL and discretization error RD                   
!                                                                     
!     Based on eq. (19) of                                            
!       R. Lindh, P.-Aa. Malmqvist and L. Gagliardi                   
!       "Molecular integrals by numerical quadrature",                
!       Theor. Chem. Acc. 106 (2001) 178-187                          
!                                                                     
!     The variable U = AL*R*R is found by a Newton-Raphson search.    
!                                                                     
!     Written by T. Saue July 2002                                    
!                                                                     
!

!      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
REAL(REALK)            :: OUTERR90,AL,RD
INTEGER                :: L
!
REAL(realk),PARAMETER  :: ACC=1.0E-6_realk
REAL(realk),PARAMETER  :: D0=0.0E0_realk,D1=1.0E0_realk,D2=2.0E0_realk,TI=1.0E1_realk
REAL(realk),PARAMETER  :: PI = 3.141592653589793E00_realk
!REAL,PARAMETER  :: ACC=1.0E-6_realk
!REAL,PARAMETER  :: D0=0.0E0_realk,D1=1.0E0_realk,D2=2.0E0_realk,TI=1.0E1_realk
!REAL,PARAMETER  :: PI = 3.141592653589793E00_realk
INTEGER,PARAMETER      :: MXIT=8
REAL(realk)            :: TOLEN,FAC,EXPL,A,ALN,RLN,U,F0HLN
REAL(realk)            :: F1HLN,DX
!REAL            :: TOLEN,FAC,EXPL,A,ALN,RLN,U,F0HLN
!REAL            :: F1HLN,DX,OUTERR90,AL,RD
INTEGER                :: I,IT
!     Initialization

TOLEN = D2
FAC = D1
DO I = 1,L
   TOLEN = TOLEN*D2
   FAC   = FAC*(2*I+1)
ENDDO
EXPL = (2*L+1)/D2
A = SQRT(PI)*FAC/TOLEN
ALN = LOG(A)
RLN = LOG(RD)
U = 35.0E0_realk
!     Newton-Raphson search
DO IT = 1,MXIT
   F0HLN = ALN+EXPL*LOG(U)-U-RLN
   F1HLN = EXPL/U-D1
   DX = F0HLN/F1HLN
   U = U - DX
   IF(ABS(DX).LT.ACC) THEN
      OUTERR90 = SQRT(U/AL)
      RETURN
   ENDIF
ENDDO
STOP 'Error in OUTERR90'

END FUNCTION OUTERR90

! Trond: polynomial evaluation.
FUNCTION POLVAL90(NORDER,B,XVAL)
use precision
IMPLICIT NONE
REAL(REALK):: B(NORDER+1),POLVAL90,XVAL
!REAL:: B(NORDER+1),POLVAL90,XVAL
INTEGER    :: XBUF,I,NORDER
POLVAL90 = B(1)
XBUF = 1
DO I = 2,(NORDER+1)
   XBUF = XBUF*XVAL
   POLVAL90 = POLVAL90 + XBUF*B(I)
ENDDO
END FUNCTION POLVAL90

