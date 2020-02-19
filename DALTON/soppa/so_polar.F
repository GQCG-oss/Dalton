C
C  /* Deck so_polar */
      SUBROUTINE SO_POLAR(MODEL,ISYMTR,IDIP,LABEL1,NLBTOT,T2AM,LT2AM,
     &                    DENSIJ,LDENSIJ,
     &                    DENSAB,LDENSAB,DENSAI,LDENSAI,POLDD,POLDQ,
     &                    POLDL,POLDA,WORK,LWORK)
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
      use so_info, only: fn_rdens, sop_stat_trh
      use so_data, only: fileinf
C
#include "implicit.h"
#include "priunit.h"
C
#include "soppinf.h"
#include "cbilnr.h"
#include "cbiexc.h"
#include "absorp.h"
#include "ccsdsym.h"
C Get irat
#include "iratdef.h"
C
      DIMENSION POLDD(2,3,3,NFRVAL), POLDQ(2,3,3,3,NFRVAL)
      DIMENSION POLDL(2,3,3,NFRVAL), POLDA(2,3,3,NFRVAL)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), DENSAI(LDENSAI)
      DIMENSION T2AM(LT2AM)
      DIMENSION WORK(LWORK)
C
      CHARACTER*8 LABEL1, LABEL2
      CHARACTER*5 MODEL
      CHARACTER*8 RTNLBL(2)
      DIMENSION   SNDPRP(2)
C
      LOGICAL   IMAGPROP, OLDDX
      PARAMETER (DP5=0.5D0)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_POLAR')
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
      LPRP1 = N2BST(ISYMTR)
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
      LPDENSIJ = NIJDEN(ISYMTR)
      LPDENSAI = NAIDEN(ISYMTR)
      IF (MODEL.EQ.'AORPA') THEN
         LPDENSAB = 0
         LPDENSTOT = LPDENSAI
      ELSE
         LPDENSAB = NABDEN(ISYMTR)
         LPDENSTOT = LPDENSIJ + LPDENSAB + LPDENSAI
      ENDIF
      LURDENS = -1
      CALL GPOPEN(LURDENS, FN_RDENS, 'OLD','DIRECT','UNFORMATTED',
     &            IRAT*LPDENSTOT,OLDDX)
C
C=============================================
C     Loop over the second property operators.
C=============================================
C
      DO 200 IPRLBL = 1, NLBTOT
C      DO 200 IPRLBL = 1, 3
C
C---------------------------------------------------
C        Find label and symmetry of second operator.
C---------------------------------------------------
C
         LABEL2 = LABAPP(IPRLBL)
         KSYM   = LABSYM(IPRLBL)
C
C--------------------------------------------------------
C        If symmetry of first operator equals symmetry of
C        second operator, that is if ISYMTR = KSYM, then.
C--------------------------------------------------------
C
         IF (KSYM .EQ. ISYMTR) THEN
C
C--------------------------------------------------
C           Get the property integrals in MO basis.
C--------------------------------------------------
C
            CALL SO_ONEPMO(WORK(KPRP1),LPRP1,LABEL2,ISYMTR,
     &                     RTNLBL,WORK(KEND1),LWORK1)
C
            IMAGPROP = RTNLBL(2).EQ.'ANTISYMM'
C
            DFACTOR = -1.0D0
            IF (IMAGPROP) DFACTOR = 1.0D0
C
C===============================================
C           Form second order properties SNDPRP.
C===============================================
C
            IF (.NOT. ABSORP) THEN
C
               DO 100 IFRVAL = 1, NFRVAL
C
C  These allocations as such doesn't change...

                  KPDENSIJ = KEND1
                  KPDENSAB = KPDENSIJ + LPDENSIJ
                  KPDENSAI = KPDENSAB + LPDENSAB
                  KEND2    = KPDENSAI + LPDENSAI
                  LWORK2   = LWORK - KEND2
C
                  CALL SO_MEMMAX ('SO_POLAR.2',LWORK2)
                  IF (LWORK2 .LT. 0)
     &                     CALL STOPIT('SO_POLAR.2',' ',KEND2,LWORK)
C
C-----------------------------------------------------------
C           Get the perturbed density matrix from file
C-----------------------------------------------------------
C
                  IF (MODEL.EQ.'AORPA') THEN
                     CALL SO_REAVE(WORK(KPDENSAI),LPDENSAI,ISYMTR,
     &                             LABEL1,FRVAL(IFRVAL),LURDENS)
                     CALL DZERO(WORK(KPDENSIJ),LPDENSIJ)
                  ELSE
                     CALL SO_REAVE(WORK(KPDENSIJ),LPDENSTOT,ISYMTR,
     &                             LABEL1,FRVAL(IFRVAL),LURDENS)
                  ENDIF
C
C---------------------------------------------------------------------
C           Calculate second order properties SNDPRP.
C---------------------------------------------------------------------
C
                  CALL SO_PROPMO(ISYMTR,SNDPRP(1),
     &                           MODEL.NE.'AORPA',IMAGPROP,
     &                           WORK(KPRP1),LPRP1,
     &                           WORK(KPDENSIJ),LPDENSIJ,
     &                           WORK(KPDENSAB),LPDENSAB,
     &                           WORK(KPDENSAI),LPDENSAI)
C
                  IF (IPRSOP .GT. 4) THEN
                     WRITE (LUPRI,'(/,A,F15.8)')
     &                   ' Frequency = ',FRVAL(IFRVAL)
                     WRITE (LUPRI,'(4A,F15.8)')
     &                   ' Second order property for ',
     &                   LABEL2,LABEL1,' = ',SNDPRP(1)
                  ENDIF
C
C---------------------------------------------------------------------
C                 Write properties into the various property matrices.
C---------------------------------------------------------------------
C
CRF  all the label seems to have coordinate first
                  IDX = COORD_IDX(LABEL2(1:1))
                  IF (IDX .LT. 1 .OR. IDX .GT. 3 ) GOTO 888
C
                  IF (LABEL2(2:7).EQ.'DIPLEN') THEN
C
                      POLDD(1,IDIP,IDX,IFRVAL) = SNDPRP(1)
C
                  ELSE IF (LABEL2(3:8).EQ.'THETA ') THEN
C
                     JDX = COORD_IDX(LABEL2(2:2))
                     IF (JDX .LT. 1 .OR. JDX .GT. 3 ) GOTO 888
                     POLDQ(1,IDIP,IDX,JDX,IFRVAL) = SNDPRP(1)
                     IF (IDX .NE. JDX )
     &                   POLDQ(1,IDIP,JDX,IDX,IFRVAL) = SNDPRP(1)
C
                  ELSE IF (LABEL2(2:7).EQ.'LONMAG') THEN
C
C----------------------------------------------------------------
C                 Multiply with minus the Bohr-magneton (-0.5) to
C                 create the magnetic dipole operator from the
C                 angular momentum operator.
C----------------------------------------------------------------
C
                     POLDL(1,IDIP,IDX,IFRVAL) = -DP5*SNDPRP(1)
C
                  ELSE IF (LABEL2(2:7).EQ.'ANGMOM') THEN
C
                     POLDA(1,IDIP,IDX,IFRVAL) = -DP5*SNDPRP(1)
C
                  END IF
C
  100          CONTINUE
C
            END IF
C
         END IF
C
  200 CONTINUE
C
C--------------------------------------
C           Close files with densities.
C--------------------------------------
C
      CALL GPCLOSE(LURDENS,'DELETE')
      call fileinf%empty() 
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('SO_POLAR')
      RETURN
C-------------------
C     Error handling
C-------------------
C
  888 CONTINUE ! FIRST CHARACTER of LABEL2 not X, Y or Z
      WRITE(LUPRI,'(A,A)') 'SO_POLAR:: ERROR in integral label', LABEL2


      CONTAINS
C
C        A function that converts the letters X,Y and Z to
C        the integers 1,2 and3
         PURE FUNCTION COORD_IDX(A)
            CHARACTER(LEN=1), INTENT(IN) :: A
            INTEGER :: COORD_IDX
            COORD_IDX = ICHAR(A) - ICHAR('X') + 1
            RETURN
         END FUNCTION
C
      END
