C
C  /* Deck so_densai1 */
      SUBROUTINE SO_DENSAI1(DENSAI,LDENSAI,DSRHF,LDSRHF,CMO,LCMO,
     &                     T2M1,LT2M1,ISYMD,ISYDIS,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Start calculating the second-order SOPPA one-electron
C              density matricx Dai by adding up integrals multiplied
C              with T2 amplitudes. 
C              Multiplication with MO-energies and factors is
C              done in SO_DENSAI2
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION DENSAI(LDENSAI), DSRHF(LDSRHF), CMO(LCMO), T2M1(LT2M1)
      DIMENSION WORK(LWORK)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_DENSAI1')
C
C--------------------------------------------------------------------
C     Add the two contribution corresponding to the G and H terms in
C     eqs (29) and (30) in J. Chem. Phys. 104, p. 4157 (1996) (except
C     that the factor 2 in these equations is left out).
C--------------------------------------------------------------------
C
      IONE = 1
C
      IF ( LT2M1 .GT. 0 ) THEN
C
         CALL CCRHS_G(DSRHF,DENSAI,CMO,IONE,CMO,IONE,T2M1,
     &                WORK,LWORK,ISYDIS,ISYMD,IONE)
C
         CALL CCRHS_H(DSRHF,DENSAI,CMO,CMO,T2M1,WORK,LWORK,
     &                ISYDIS,ISYMD,IONE)
C
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_DENSAI1')
C
      RETURN
      END
