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
 INTEGER           :: MAXJ0,IADR,IPOINT,IORDER,JADR,JMAX,J,IADR2
 REAL(REALK)       :: DENOM,D2MAX1,R2MAX1,TERM,SUM,REXPW,WVAL,D2WAL
 JMAX = JMX + 3
 MAXJ0 = JMAX
 !
 !     WVAL = 0.0
 !
 IADR = 1
 DENOM = 1E0_realk
 DO J = 0,JMAX
    TABFJW(J,0) = 1E0_realk/DENOM
    IADR = IADR + 1201
    DENOM = DENOM + 2E0_realk
 ENDDO
 !
 !     WVAL = 0.1, 0.2, 0.3,... 12.0
 !
 IADR = IADR - 1201
 !D2MAX1 = DFLOAT(2*JMAX + 1)
 D2MAX1 = 2*JMAX + 1
 R2MAX1 = 1E0_realk/D2MAX1
!$OMP PARALLEL DO PRIVATE(IPOINT,WVAL,D2WAL,IADR2,TERM,SUM,DENOM,&
!$OMP IORDER,REXPW,JADR,J) FIRSTPRIVATE(IADR,R2MAX1,D2MAX1,&
!$OMP JMAX) SHARED(TABFJW) SCHEDULE(DYNAMIC,5)
 DO IPOINT = 1,1200
    !  WVAL = TENTH*DFLOAT(IPOINT)
    WVAL = 0.01E0_realk*IPOINT
    D2WAL = WVAL + WVAL
    IADR2 = IADR + IPOINT
    TERM = R2MAX1
    SUM = R2MAX1
    DENOM = D2MAX1
    DO IORDER = 2,200
       DENOM = DENOM + 2E0_realk
       TERM = TERM*D2WAL/DENOM
       SUM = SUM + TERM
       IF (TERM .LE. 1.0E-15_realk) EXIT
    ENDDO
    REXPW = EXP(-WVAL)
    TABFJW(JMAX,IPOINT) = REXPW*SUM
    DENOM = D2MAX1
    JADR = IADR2
    DO J = 1,JMAX
       DENOM = DENOM - 2E0_realk
       TABFJW(JMAX-J,IPOINT) = (TABFJW(JMAX-J+1,IPOINT)*D2WAL + REXPW)/DENOM
       JADR = JADR - 1201
    ENDDO
 ENDDO
!$OMP END PARALLEL DO 
END SUBROUTINE gammaTabulation

END MODULE IchorGammaTabulationModule
