C
C  /* Deck so_moenergy */
      SUBROUTINE SO_MOENERGY(FOCKD,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Henrik Koch and Keld Bak, December 1995
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Get MO-energies
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION FOCKD(NORBTS)
      DIMENSION WORK(LWORK)
C     LOGICAL   OPENY
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <inftap.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_MOENERGY')
C
C-------------------------------------
C     Read canonical orbital energies.
C-------------------------------------
C
      IF (LUSIFC .LE. 0) CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ',
     &                               'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC) POTNUC,EMY,EACTIV,EMCSCF,ISTATE,ISPIN,NACTEL,
     &              LSYM
C
      ESCF = EMCSCF
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (FOCKD(I), I=1,NORBTS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C-------------------------------------------------------------------
C     If frozen orbitals delete elements refering to those orbitals.
C-------------------------------------------------------------------
C
      IF (FROIMP .OR. FROEXP) CALL CCSD_DELFRO(FOCKD,WORK,LWORK)
C
C----------------------------------------------------------
C     Change symmetry-ordering of the Fock-matrix diagonal.
C----------------------------------------------------------
C
CSPAS:11/11-2003:
C     CALL FOCK_REORDER(FOCKD,WORK,LWORK,IPRINT)
      CALL FOCK_REORDER(FOCKD,WORK,LWORK)
CKeinSPASmehr
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_MOENERGY')
C
      RETURN
      END
