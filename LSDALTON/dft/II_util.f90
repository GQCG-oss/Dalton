SUBROUTINE LSFORT_WRT(lupri,str, len)
  INTEGER   lupri,len
  CHARACTER str(len)
  WRITE(LUPRI,*) str
END SUBROUTINE LSFORT_WRT

SUBROUTINE II_RADLMG(RADNDE,RADWT,NR,RADERR,NRADPT,NUCORB,AA,NHTYP,IPRINT)
!grid-generation subroutine
use precision
IMPLICIT NONE
REAL(REALK), PARAMETER  :: D0 = 0E0_realk, D1 = 1E0_realk,D2 = 2E0_realk, D3 = 3E0_realk
!REAL, PARAMETER  :: D0 = 0E0_realk, D1 = 1E0_realk,D2 = 2E0_realk, D3 = 3E0_realk
INTEGER     :: NHTYP
REAL(REALK) :: AA(2,NHTYP,2),RADERR
!REAL :: AA(2,NHTYP,2),RADERR
REAL(REALK) :: RADWT(NRADPT), RADNDE(NRADPT)
!REAL :: RADWT(NRADPT), RADNDE(NRADPT)
INTEGER     :: NUCORB(NHTYP,2),NR,NRADPT,IPRINT
!
INTEGER     :: LL,L,NBAS,IR,I
REAL(REALK) :: AH,H,EPH,RL,RH,AL,HTMP,RHTMP,GRDC,DUMMY
REAL(realk) :: DISERR90,OUTERR90
!REAL :: AH,H,EPH,RL,RH,AL,HTMP,RHTMP,GRDC,DUMMY,DISERR90,OUTERR90
     
DUMMY=1E+20_realk
!     Grid spacing to H and inner grid point to AH
NR = 0
H  = DUMMY
AH = 0E0_realk
DO LL = 1,NHTYP
   L = LL-1
   NBAS=NUCORB(LL,1)+NUCORB(LL,2)
   IF(NBAS.GT. 0) THEN
      HTMP = DISERR90(RADERR,L) !function
      H = MIN(H,HTMP)
!      IF(IPRINT.GE. 3) WRITE(LUPRI,*) L,'-orbitals --> ',HTMP
   ENDIF
   AH = MAX(AH,AA(1,LL,1))
ENDDO
IF(AH .EQ. 0E0_realk) RETURN
EPH = EXP(H)
!IF(IPRINT.GE. 3)WRITE(LUPRI,'(1X,A,F6.3)') 'Value chosen: ',H
!...  Inner grid point AA->R transformation.
AH = D2*AH
!IF(IPRINT.GE. 3)WRITE(LUPRI,*) 'AH = ',AH
RL = ((1.9E0_realk+LOG(RADERR))/D3)-(LOG(AH)/D2)
RL = EXP(RL)
!IF(IPRINT.GE. 3) WRITE(LUPRI,'(A,1P,E12.5)')'* Inner grid point:',RL
!...  Outer point
!IF(IPRINT.GE. 3) WRITE(LUPRI,'(A)') '* Outer point:'
RH = D0
DO LL = 1,NHTYP
   L = LL-1
   AL=DUMMY
   IF(NUCORB(LL,1).GT. 0) AL=AA(2,LL,1)
   IF(NUCORB(LL,2).GT. 0) AL=MIN(AL,AA(2,LL,2))
   IF(AL.LT.DUMMY) THEN
      AL = AL+AL
      RHTMP = OUTERR90(AL,L,RADERR) !function          
      RH=MAX(RH,RHTMP)
 !     IF(IPRINT.GE. 3) THEN
 !        WRITE(LUPRI,*) L,'-orbitals --> ',RHTMP
 !     ENDIF
   ENDIF
ENDDO
!IF(IPRINT.GE. 3)WRITE(LUPRI,'(1X,A,F6.3)') 'Value chosen: ',RH
GRDC = RL/(EPH-D1)
!IF(IPRINT.GE. 3)WRITE(LUPRI,'(A,E12.5)') 'Constant c:',GRDC
NR = NINT(LOG(D1+(RH/GRDC))/H)
!IF(IPRINT.GE. 3)WRITE(LUPRI,'(A,I5)') 'Number of points:',NR
IF(NR.GT.NRADPT) CALL LSQUIT('Too many radial points.',-1)
RADNDE(NR) = RL
RADWT(NR)  = (RL+GRDC)*RL*RL*H
DO IR = NR-1,1,-1
   RADNDE(IR) = (RADNDE(IR+1)+GRDC)*EPH-GRDC
   RADWT(IR) = (RADNDE(IR)+GRDC)*RADNDE(IR)*RADNDE(IR)*H
ENDDO
END SUBROUTINE II_RADLMG


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

SUBROUTINE II_NUCBAS(NUCORB,AA,NCENT,NHKT,NUCO,NHTYP,NUCIND,KMAX,&
     & MXPRIM,PRIEXP,NSTART,IPRINT)
!grid-generation
use precision
IMPLICIT NONE
!  Extract basis information for all centers
!  Written by T.Saue March 12 2001
INTEGER               :: NHTYP,NUCIND,KMAX,IPRINT,MXPRIM
INTEGER               :: NCENT(KMAX),NHKT(KMAX),NUCO(KMAX),NSTART(KMAX)
REAL(REALK)           :: PRIEXP(MXPRIM)
!REAL           :: PRIEXP(MXPRIM)
!
REAL(REALK),PARAMETER :: D0=0.0E0_realk
REAL(REALK)           :: AA(2,NHTYP,2,NUCIND),A,DUMMY
!REAL,PARAMETER :: D0=0.0E0_realk
!REAL           :: AA(2,NHTYP,2,NUCIND),A,DUMMY
INTEGER               :: NUCORB(NHTYP,2,NUCIND)
INTEGER  :: NDIM,JCENT,JC,JLVAL,ISHELL,ICENT,IC,ILVAL,IPRIM,JPRIM,I,IEXP,LL,L
INTEGER  :: NPRIM,LUPRI
!     Initialize

!      PRINT*,'NHTYP=',NHTYP
!      PRINT*,'KMAX=',KMAX
!      PRINT*,'MXPRIM=',MXPRIM
!      PRINT*,'NUCIND=',NUCIND
!      DO ISHELL = 1,KMAX
!         print*,'NCENT(',ISHELL,')=',NCENT(ISHELL)
!      ENDDO
!      DO ISHELL = 1,KMAX
!         print*,'NHKT(',ISHELL,')=',NHKT(ISHELL)
!      ENDDO
!      DO ISHELL = 1,KMAX
!         print*,'NUCO(',ISHELL,')=',NUCO(ISHELL)
!      ENDDO
!      DO ISHELL = 1,KMAX
!         print*,'NSTART(',ISHELL,')=',NSTART(ISHELL)
!      ENDDO
!      DO ISHELL = 1,MXPRIM
!         print*,'PRIEXP(',ISHELL,')=',PRIEXP(ISHELL)
!      ENDDO

NDIM = 2*NHTYP*NUCIND
CALL LS_IZERO(NUCORB,NDIM)
DUMMY=1.0E+20_realk
JCENT = 0
JPRIM = -1
JC = -1
JLVAL = -1
DO ISHELL = 1,KMAX
   ICENT = NCENT(ISHELL)
 !  print*,'ISHELL',ISHELL,'ICENT',ICENT
   IF(ICENT.NE.JCENT) THEN
      JCENT = ICENT
      JLVAL = 0
   ENDIF
!   IC = LCLASS(ISHELL)
    IC = 1
   IF(IC.NE.JC) THEN
      JC    = IC
      JLVAL = 0
   ENDIF
   ILVAL = NHKT(ISHELL)
   IF(ILVAL.NE.JLVAL) THEN
      JLVAL = ILVAL
      NUCORB(ILVAL,IC,ICENT) = 0
!      print*,'NUCORB(ILVAL,IC,ICENT)',NUCORB(ILVAL,IC,ICENT)
      AA(1,ILVAL,IC,ICENT)=D0
      AA(2,ILVAL,IC,ICENT)=DUMMY
 !     print*,'AA(1,',ILVAL,',',IC,',',ICENT,')',AA(1,ILVAL,IC,ICENT),'AA(2,',ILVAL,',',IC,',',ICENT,')',AA(2,ILVAL,IC,ICENT)
   ENDIF
   NUCORB(ILVAL,IC,ICENT)=NUCORB(ILVAL,IC,ICENT)+1
!   print*,'NUCORB(ILVAL,IC,ICENT)',NUCORB(ILVAL,IC,ICENT)
   IPRIM = NSTART(ISHELL)
   IF(IPRIM.NE.JPRIM) THEN
      JPRIM = IPRIM
      NPRIM = NUCO(ISHELL)
      DO IEXP = 1,NPRIM
         A=PRIEXP(IPRIM+IEXP)            
         AA(1,ILVAL,IC,ICENT)=MAX(AA(1,ILVAL,IC,ICENT),A)
         AA(2,ILVAL,IC,ICENT)=MIN(AA(2,ILVAL,IC,ICENT),A)
!            print*,'AA(1,',ILVAL,',',IC,',',ICENT,')',AA(1,ILVAL,IC,ICENT),'AA(2,',ILVAL,',',IC,',',ICENT,')',AA(2,ILVAL,IC,ICENT)
      ENDDO
   ENDIF
ENDDO
!C
      LUPRI=2
!IF(IPRINT.GE. 2) THEN
!   WRITE(LUPRI,*)'NUCORB:Basis set information:'
!   DO I = 1,NUCIND
!      WRITE(LUPRI,'(/A,I4,A/)') '*** Center: ',I,' ***'
!      WRITE(LUPRI,'(2X,A)') '* Large components:'
!      IC = 1
!      DO LL = 1,NHTYP
!         IF(NUCORB(LL,IC,I).GT. 0) THEN
!            print*,'NUCORB(',LL,',',IC,',',I,') =',NUCORB(LL,IC,I)
!            L=LL-1
!            print*, L,'-orbitals:',NUCORB(LL,IC,I),&
!            &      'Alpha_H(1,',LL,',1,',I,')=',AA(1,LL,IC,I),&
!            &      'Alpha_L(2,',LL,',1,',I,')=',AA(2,LL,IC,I)
!         else
!            write(lupri,*)'NUCORB(',LL,',',IC,',',I,') GT 0',NUCORB(LL,IC,I)
!         ENDIF
!      ENDDO
!   ENDDO
!ENDIF
!C
 END SUBROUTINE II_NUCBAS

!!$ SUBROUTINE II_GET_RADIIOF_SHELLS(RSHEL,KMAX,DFTHRI,JSTRT,PRIEXP,NUCO,NHKT)
!!$! get radii of all shells as defined by specified threshold DFTHRI.
!!$IMPLICIT NONE
!!$#include <aovec.h>
!!$#include <maxorb.h>
!!$#include <primit.h>
!!$#include <shells.h>
!!$#include <priunit.h>
!!$#include <dftcom.h>
!!$!debug below
!!$#include <mxcent.h>
!!$#include <nuclei.h>
!!$REAL(REALK)  :: FACL(10),RSHEL(KMAX),DFTHRI,THLOG,R2
!!$INTEGER      :: ISHELA,KMAX,JSTRT(KMAX),IAOS
!!$INTEGER      :: NUCO(KMAX),NHKT(KMAX),!,NSTART(KMAX),NCENT(KMAX),
!!$!JSTRT is the same as NSTART in II_NUCBAS
!!$FACL(1)= 1E0_realk
!!$FACL(2)= 1.3333E0_realk
!!$FACL(2)= 1.6E0_realk
!!$FACL(2)= 1.83E0_realk
!!$FACL(2)= 2.03E0_realk
!!$FACL(2)= 2.22E0_realk
!!$FACL(2)= 2.39E0_realk
!!$FACL(2)= 2.55E0_realk
!!$FACL(2)= 2.70E0_realk
!!$FACL(2)= 2.84E0_realk
!!$     
!!$THLOG = LOG(DFTHRI)
!!$DO ISHELA=1,KMAX
!!$   JSTA = JSTRT(ISHELA)
!!$   RSHEL(ISHELA) =  0E0_realk
!!$   NUMCFA = NUMCF(ISHELA)!WARNING
!!$   DO IAOS = JSTA+1, JSTA + NUCO(ISHELA)
!!$      IF(ABS(PRICCF(IAOS,NUMCFA)).GT. 0E0_realk) THEN !WARNING
!!$         R2 = (LOG(ABS(PRICCF(IAOS,NUMCFA)))-THLOG)/PRIEXP(IAOS)
!!$         IF(RSHEL(ISHELA).LT.R2) RSHEL(ISHELA) = R2
!!$      END IF
!!$   END DO
!!$END DO
!!$DO ISHELA=1,KMAX
!!$   RSHEL(ISHELA) = SQRT(RSHEL(ISHELA))*FACL(NHKT(ISHELA))
!!$END DO
!!$END SUBROUTINE II_GET_RADIIOF_SHELLS


SUBROUTINE II_GETBLOCKS(CENTER,CELLSZ,RSHEL,KMAX,NCENT,NBAST,NATOM,X,Y,Z,NBLCNT,IBLCKS)
use precision
IMPLICIT NONE
  !     get blocks of active SHELLS in cube of CELLSZ size centered at
  !     CENTER.
  !
  !     RSHEL2 - precomputed shell extents (squared).
  !     NBLCNT (output) - number of active blocks
  !     IBLCKS (output) - pairs of (startindex, stopindex)
  
INTEGER      ::  NBAST,KMAX,NATOM,ISHELA,I,J
INTEGER      ::  IBLCKS(2,NBAST),ISHLEN,NBLCNT,IPREV,ICENT,NCENT(KMAX)
REAL(REALK)  ::  CENTER(3),RSHEL(KMAX),CELLSZ,CELLDG,X(NATOM),Y(NATOM),Z(NATOM)
REAL(REALK)  ::  DST,PX,PY,PZ 
!REAL  ::  CENTER(3),RSHEL(KMAX),CELLSZ,CELLDG,X(NATOM),Y(NATOM),Z(NATOM)
!REAL  ::  DST,PX,PY,PZ 

  NBLCNT = 0
  IPREV = -1111
  CELLDG = CELLSZ*0.5E0_realk*SQRT(3E0_realk)
  DO ISHELA=1,KMAX
     ICENT = NCENT(ISHELA)
     PX = center(1)-X(ICENT)
     PY = center(2)-Y(ICENT)
     PZ = center(3)-Z(ICENT)
     DST = SQRT(PX*PX + PY*PY + PZ*PZ)
     IF(DST.LE.RSHEL(ISHELA)+CELLDG) THEN
!       accepted...
        IF(ISHELA.NE.IPREV+1) THEN
           NBLCNT = NBLCNT + 1
           IBLCKS(1,NBLCNT) = ISHELA
        END IF
        IPREV = ISHELA
        IBLCKS(2,NBLCNT) = ISHELA
!        print "('shell',2I3,' at ',3F6.2,' rad. ',F6.2,' accepted')",&
!     &         ISHELA,ICENT,X(ICENT),Y(ICENT),Z(ICENT),SQRT(RSHEL(ISHELA))

     END IF
!     print "('shell',2I3,' at ',3F6.2,' rad. ',F6.2,' rejected')",&
!     &        ISHELA,ICENT,X(ICENT),Y(ICENT),Z(ICENT),RSHEL(ISHELA)
  END DO !shellloop
 ! print "('cell at:',3F15.5,' blocks:',I3)", CENTER(1),CENTER(2),CENTER(3),NBLCNT
 ! print "(8('[',2I3,']'))", ((IBLCKS(I,J),I=1,2),J=1,NBLCNT)

END SUBROUTINE II_GETBLOCKS

