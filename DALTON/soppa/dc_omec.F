C  /* Deck dc_omes */
      SUBROUTINE DC_OMEC(OMEC,RES2,LRES2,OME,ISYMTR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     PURPOSE: Calculate 2p-2h part of double corrected RPA excitation
C     energy. See eq. (35) in RPA(D) paper.
C
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (ONE    = 1.0D0, TWO = 2.0D0 )
      DIMENSION RES2(LRES2),WORK(LWORK)
C
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "soppinf.h"
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('DC_OMEC')
C
C-----------------------
C     Memory allocation.
C-----------------------
C
      LEDIA1 = NT1AM(ISYMTR)
      LEDIA2 = N2P2HOP(ISYMTR)
C
      KEDIA1  = 1
      KEDIA2  = KEDIA1 + LEDIA1
      KEND    = KEDIA2 + LEDIA2
      LWORK1  = LWORK - KEND
C
      CALL SO_MEMMAX ('DC_OMEC',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('DC_OMEC',' ',KEND,LWORK)
C
C-------------------------------------
C     Read diagonal E[2] elements
C-------------------------------------
C
      CALL GPOPEN  (LUDIAG,'SO_DIAG','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
      REWIND LUDIAG
C
      READ(LUDIAG) ( WORK(KEDIA1+I-1), I = 1,LEDIA1)
      READ(LUDIAG) ( WORK(KEDIA2+I-1), I = 1,LEDIA2)
C
      CALL GPCLOSE (LUDIAG,'KEEP')
C
C----------------------------
C     Calculate contribution.
C----------------------------
C
      OMEC = 0.0D0
C
      DO IABIJ = 1,LRES2
C
         IABIJ2 = IABIJ + LEDIA1
C
         DEN = - ( WORK(IABIJ2) - OME )
C
         XIDEN = ONE / DEN
C
         OMEC = OMEC + RES2(IABIJ) * RES2(IABIJ) * XIDEN
C
C-----------------------------------------------------------
C     Calculate the first-order RPA(D) eigenvector (eq 16)
C-----------------------------------------------------------
C
         RES2(IABIJ) = RES2(IABIJ) * XIDEN
C
      END DO
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('DC_OMEC')
C
      RETURN
      END
