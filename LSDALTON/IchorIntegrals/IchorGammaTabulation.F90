!> @file
!> Contains common routine for the Ichor Code
!> \brief Contains common routine for the Ichor Code
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorGammaTabulationModule
use IchorprecisionModule
public:: gammaTabulation
private
CONTAINS
!> \brief tabulation of the boys function or incomplete gamma function \f$ F_{n} \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param LUPRI logical unit number for output printing
!> \param JMX tabulate to this order of angular momentum
!> \param sharedTUV stores the tabulated values
!>
!> ***** Tabulation of incomplete gamma function *****
!>   For J = JMX a power series expansion is used, see for
!>   example Eq.(39) given by V. Saunders in "Computational
!>   Techniques in Quantum Chemistry and Molecular Physics",
!>   Reidel 1975.  For J < JMX the values are calculated
!>   using downward recursion in J.
!>
SUBROUTINE gammaTabulation(LUPRI,JMX,nTABFJW1,nTABFJW2,TABFJW)
 IMPLICIT NONE
 INTEGER,intent(in)       :: JMX,LUPRI,nTABFJW1,nTABFJW2
 REAL(realk)              :: TABFJW(0:nTABFJW1,0:nTABFJW2)
!
 INTEGER           :: MAXJ0,IADR,IPOINT,IORDER,JADR,JMAX,J
 REAL(REALK)       :: DENOM,D2MAX1,R2MAX1,TERM,SUM,REXPW,WVAL,D2WAL
 REAL(REALK), PARAMETER :: HALF = 0.5E0_realk,  TEN6 = 1.0E6_realk
 REAL(REALK), PARAMETER :: D1 = 1E0_realk, D10 = 10E0_realk
 REAL(REALK), PARAMETER :: D2 = 2E0_realk, D4 = 4E0_realk, D12 = 12E0_realk, TENTH = 0.01E0_realk
 REAL(REALK), PARAMETER :: COEF2 = HALF,  COEF3 = - D1/6E0_realk, COEF4 = D1/24E0_realk
 REAL(REALK), PARAMETER :: COEF5 = - D1/120E0_realk, COEF6 = D1/720E0_realk

 REAL(REALK), PARAMETER :: PI    = 3.14159265358979323846E00_realk
 REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730E00_realk
 REAL(REALK), PARAMETER :: R2PI52 = 5.91496717279561287782E00_realk
 REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2
 REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI

 REAL(REALK), PARAMETER :: GFAC30 =  .4999489092E0_realk 
 REAL(REALK), PARAMETER :: GFAC31 = -.2473631686E0_realk
 REAL(REALK), PARAMETER :: GFAC32 =  .321180909E0_realk
 REAL(REALK), PARAMETER :: GFAC33 = -.3811559346E0_realk
 REAL(REALK), PARAMETER :: GFAC20 = .4998436875E0_realk
 REAL(REALK), PARAMETER :: GFAC21 = -.24249438E0_realk
 REAL(REALK), PARAMETER :: GFAC22 =  .24642845E0_realk
 REAL(REALK), PARAMETER :: GFAC10 =  .499093162E0_realk
 REAL(REALK), PARAMETER :: GFAC11 = -.2152832E0_realk
 REAL(REALK), PARAMETER :: GFAC00 =  .490E0_realk

 JMAX = JMX + 3!6
 MAXJ0 = JMAX
 !
 !     WVAL = 0.0
 !
 IADR = 1
 DENOM = D1
 DO J = 0,JMAX
    TABFJW(J,0) = D1/DENOM
    IADR = IADR + 1201
    DENOM = DENOM + D2
 ENDDO
 !
 !     WVAL = 0.1, 0.2, 0.3,... 12.0
 !
 IADR = IADR - 1201
 !D2MAX1 = DFLOAT(2*JMAX + 1)
 D2MAX1 = 2*JMAX + 1
 R2MAX1 = D1/D2MAX1
 DO IPOINT = 1,1200
    !  WVAL = TENTH*DFLOAT(IPOINT)
    WVAL = TENTH*IPOINT
    D2WAL = WVAL + WVAL
    IADR = IADR + 1
    TERM = R2MAX1
    SUM = TERM
    DENOM = D2MAX1
    DO IORDER = 2,200
       DENOM = DENOM + D2
       TERM = TERM*D2WAL/DENOM
       SUM = SUM + TERM
       IF (TERM .LE. 1.0E-15_realk) EXIT
    ENDDO
    REXPW = EXP(-WVAL)
    TABFJW(JMAX,IPOINT) = REXPW*SUM
    DENOM = D2MAX1
    JADR = IADR
    DO J = 1,JMAX
       DENOM = DENOM - D2
       TABFJW(JMAX-J,IPOINT) = (TABFJW(JMAX-J+1,IPOINT)*D2WAL + REXPW)/DENOM
       JADR = JADR - 1201
    ENDDO
 ENDDO
END SUBROUTINE gammaTabulation

END MODULE IchorGammaTabulationModule
