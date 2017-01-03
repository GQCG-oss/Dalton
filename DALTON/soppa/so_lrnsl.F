C
C  /* Deck so_lrnsl */
      SUBROUTINE SO_LRNSL(MODEL,ISYM,T2AM,LT2AM,
     &                    DENSIJ,LDENSIJ,
     &                    DENSAB,LDENSAB,DENSAI,LDENSAI,
     &                    WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Stephan P. A. Sauer: 5.12.2003
C     Rasmus Faber: Nov. 2015 --> Implemented pert. dens. approach.
C
C     PURPOSE: Calculates the frequency dependent linear response
C              properties from the perturbed density matrices and
C              appropriate property integrals with the atomic
C              integral direct SOPPA program.
C
      use so_info, only: fn_rdens, sop_stat_trh, sop_dp
C
C#include "implicit.h"
      implicit none
#include "priunit.h"
C
#include "soppinf.h"
C From infrsp we need
C FREQ <- Array of requecies 
C NFREQ <- Number of frequencies
#include "infrsp.h"
C inflr.h below requires MAXLBL from rspprp.h
#include "rspprp.h"
C We need 
C LBLLR <- Label of operators 
C NGPLR <- Number of operators of each symmetry
#include "inflr.h"
#include "ccsdsym.h"
C REP <- Symbols of the ireducible representations.
#include "pgroup.h"
C

C
      integer, intent(in) :: isym, ! symmetry of this batch of perturbations
     &       LDENSIJ, LDENSAB, LDENSAI, LT2AM, LWORK ! array lengths
      REAL(sop_dp),INTENT(IN) :: DENSIJ(LDENSIJ), DENSAB(LDENSAB),
     &                     DENSAI(LDENSAI), T2AM(LT2AM)
      REAL(sop_dp),INTENT(INOUT) :: WORK(LWORK)
      CHARACTER(len=5), intent(in) :: MODEL
C
      CHARACTER*8 LABEL1, LABEL2
      CHARACTER*8 RTNLBL(2)
      REAL(sop_dp)  ::  PROPVAL, DFACTOR
      INTEGER :: LPRP1, LWORK1, LPDENSIJ, LPDENSAB, LPDENSAI, LPDENSTOT,
     &           LWORK2
      INTEGER :: KPRP1, KEND1, KPDENSIJ, KPDENSAB, KPDENSAI, KEND2
      INTEGER :: IDUMMY, ILAB1, ILAB2, IFREQ, LURDENS
C
      LOGICAL   IMAGPROP
      REAL(sop_dp), PARAMETER :: DP5=0.5D0
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_LRNSL')
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      KEND1 = 1
C
C----------------------------------
C     Space for MO property matrix.
C----------------------------------
C
      LPRP1 = N2BST(ISYM)
      KPRP1 = KEND1
      KEND1 = KPRP1 + LPRP1
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('SO_POLAR.1',LWORK1)
      IF (LWORK1 .LT.0) CALL STOPIT('SO_POLAR.1',' ',KEND1,LWORK)
C
C-----------------------------------------
C     Open files with perturbed densities.
C-----------------------------------------
C
      LURDENS = -1
      CALL GPOPEN(LURDENS, FN_RDENS, 'OLD',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
C
C=============================================
C     Loop over the first property operators.
C=============================================
C
      DO 200 ILAB1 = 1, NGPLR(ISYM)
C
C---------------------------------------------------
C        Find label and symmetry of second operator.
C---------------------------------------------------
C
         LABEL1 = LBLLR(ISYM,ILAB1)
C
C--------------------------------------------------
C           Get the property integrals in MO basis.
C--------------------------------------------------
C
         CALL SO_ONEPMO(WORK(KPRP1),LPRP1,LABEL1,ISYM,
     &                  RTNLBL,WORK(KEND1),LWORK1)
C
         IMAGPROP = RTNLBL(2).EQ.'ANTISYMM'
C
         DFACTOR = -1.0D0
         IF (IMAGPROP) DFACTOR = 1.0D0
C
C  These allocations as such doesn't change...
         LPDENSIJ = NIJDEN(ISYM)
         IF (MODEL.EQ.'AORPA') THEN
            LPDENSAB = 0
         ELSE
            LPDENSAB = NABDEN(ISYM)
         ENDIF
         LPDENSAI = NAIDEN(ISYM)

         KPDENSIJ = KEND1
         KPDENSAB = KPDENSIJ + LPDENSIJ
         KPDENSAI = KPDENSAB + LPDENSAB
         KEND2    = KPDENSAI + LPDENSAI
         LWORK2   = LWORK - KEND2
C
         CALL SO_MEMMAX ('SO_POLAR.2',LWORK2)
         IF (LWORK2 .LT. 0)
     &            CALL STOPIT('SO_POLAR.2',' ',KEND2,LWORK)
C
C===============================================
C           Form second order properties SNDPRP.
C===============================================
C
C
         DO 100 IFREQ = 1, NFREQ ! Change these varibles
            DO ILAB2 = 1, NGPLR(ISYM)
               LABEL2 = LBLLR(ISYM,ILAB2)
C
C-----------------------------------------------------------
C        Get the perturbed density matrix from file
C-----------------------------------------------------------
C
               CALL SO_FNDPOS(LABEL1,FREQ(IFREQ),LURDENS)
               LPDENSTOT = LPDENSIJ + LPDENSAB + LPDENSAI
               IF (MODEL.EQ.'AORPA') THEN
                  CALL SO_REAVE(WORK(KPDENSAI),LPDENSAI,ISYM,
     &                          LABEL1,FREQ(IFREQ),LURDENS)
                  CALL DZERO(WORK(KPDENSIJ),LPDENSIJ)
               ELSE
                  CALL SO_REAVE(WORK(KPDENSIJ),LPDENSTOT,ISYM,
     &                          LABEL1,FREQ(IFREQ),LURDENS)
               ENDIF
C
C---------------------------------------------------------------------
C        Calculate second order properties PROPVAL.
C---------------------------------------------------------------------
C
               CALL SO_PROPMO(ISYM,PROPVAL,
     &                        MODEL.NE.'AORPA',IMAGPROP,
     &                        WORK(KPRP1),LPRP1,
     &                        WORK(KPDENSIJ),LPDENSIJ,
     &                        WORK(KPDENSAB),LPDENSAB,
     &                        WORK(KPDENSAI),LPDENSAI)
C
               WRITE(LUPRI,'(/A,2(/A,A10,I4,3A,F10.6),/,/A,F20.12)')
     &              ' Linear response function in a.u.',
     &              ' A operator, symmetry, frequency: ',
     &              LABEL1, ISYM,'  (',REP(ISYM-1),')',-FREQ(IFREQ),
     &              ' B operator, symmetry, frequency: ',
     &              LABEL2, ISYM,'  (',REP(ISYM-1),')', FREQ(IFREQ),
     &              ' Value of linear response -<<A;B>>(omega): ',
     &              PROPVAL
C
            END DO
C
  100    CONTINUE
C
  200 CONTINUE
C
C--------------------------------------
C           Close files with densities.
C--------------------------------------
C
      CALL GPCLOSE(LURDENS,'KEEP')
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('SO_LRNSL')
      RETURN
      END
