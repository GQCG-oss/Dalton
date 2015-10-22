C
C  /* Deck so_lnraba */
      SUBROUTINE SO_LNRABA(POLDD,POLDQ,POLDL,POLDA,WORK,LWORK,PASS)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Stephan P. A. Sauer, Februrary 1999
C     Stephan P. A. Sauer: 3.12.2003
C
C     PURPOSE: Main driver routine for the calculation of frequency
C              dependent linear response properties with the atomic 
C              integral direct SOPPA program.
C
#include "implicit.h"
#include "priunit.h"
C
#include "soppinf.h"
#include "cbilnr.h"
#include "cbiexc.h"
#include "absorp.h"
#include "ccsdsym.h"
c#include "infdim.h"
#include "inforb.h"
c#include "iratdef.h"
c#include "maxaqn.h"
c#include "maxorb.h"
c#include "mxcent.h"
c#include "nuclei.h"
c#include "symmet.h"
c#include "wrkrsp.h"
C
c     PARAMETER (HALF = 0.5D0,ESUDIP = 64604.885D0,ESUECD = 471.44360D0)
C
      DIMENSION POLDD(2,3,3,NFRVAL,4), POLDQ(2,3,3,3,NFRVAL,4)
      DIMENSION POLDL(2,3,3,NFRVAL,4), POLDA(2,3,3,NFRVAL,4)
      DIMENSION WORK(LWORK)
C
      CHARACTER*5 MODEL
      CHARACTER*8 LABEL1
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_LNRABA')
C
C-------------------------------------------
C     Start timing the AO-SOPPA calculation.
C-------------------------------------------
C
      CALL TIMER('START ',TIMEIN,TIMOUT)
C
      DTIME   = SECOND()
      TIMTOT  = DTIME
      CALL GETTIM (DUMMY,WTIME)
      TIMWTO  = WTIME
C
C-------------------------------------------
C     Initializes memory statistics routine.
C-------------------------------------------
C
      CALL SO_MEMMAX ('START',LWORK)
C
C---------------------------------------------------------------
C     If only RPA initialize some pointers and sort AO-integrals
C     into distributions.
C---------------------------------------------------------------
C
      IF ( AORPA .AND. (.NOT.DCRPA) .AND. (.NOT.AOSOP) 
     &           .AND. (.NOT.AOSOC) )    CALL RP_INIT(WORK,LWORK)
C
C---------------------------------------------------
C     Initialize a few pointerarrays and dimensions.
C---------------------------------------------------
C
      DTIME     = SECOND()
      CALL SO_INIT
      DTIME     = SECOND()  - DTIME
      SOTIME(2) = SOTIME(2) + DTIME
C
C----------------------------------
C     Set print level for AO-SOPPA.
C----------------------------------
C
      IPRSOP = IPRLNR
C
      IF (IPRLNR .GE. 0) 
     &     CALL TITLER('Solving Linear Response Equations','¤',118)
C
C---------------------------------
C     1. allocation of work space.
C---------------------------------
C
      LFOCKD  = NORBT
C
      KFOCKD  = 1
      KEND1   = KFOCKD  + LFOCKD
      LWORK1  = LWORK   - KEND1
C
      IF (DCRPA .OR. AOSOP .OR. AOSOC) THEN
C
         LT2AM   = NT2AMX
         LDENSIJ = NIJDEN(1)
         LDENSAB = NABDEN(1)
         LDENSAI = NAIDEN(1)
C
         KT2AM   = KFOCKD  + LFOCKD
         KDENSIJ = KT2AM   + LT2AM
         KDENSAB = KDENSIJ + LDENSIJ
         KDENSAI = KDENSAB + LDENSAB
         KEND1   = KDENSAI + LDENSAI
         LWORK1  = LWORK   - KEND1
C
      END IF
C
      CALL SO_MEMMAX ('SO_LNRABA.1',LWORK1)
      IF (LWORK1 .LT .0) CALL STOPIT('SO_LNRABA.1',' ',KEND1,LWORK)
C
C------------------------------------------------
C     Get MO-energies (the fock matrix diagonal).
C------------------------------------------------
C
      DTIME     = SECOND()
      CALL SO_MOENERGY(WORK(KFOCKD),WORK(KEND1),LWORK1)
      DTIME     = SECOND()  - DTIME
      SOTIME(3) = SOTIME(3) + DTIME
C
C------------------------------------------------------
C     Construct property-integrals and write to LUPROP.
C------------------------------------------------------
C
      CALL SO_PRPINT('LINEAR',NLBTOT,WORK(KEND1),LWORK1)
C
C---------------------------------------------------------------------
C     Initialize arrays of transition moments and excitation energies.
C---------------------------------------------------------------------
C
      CALL DZERO(SNDPRP,2)
      CALL DZERO(POLDD,2*9*NFRVAL*4)
      CALL DZERO(POLDQ,2*27*NFRVAL*4)
      CALL DZERO(POLDL,2*9*NFRVAL*4)
      CALL DZERO(POLDA,2*9*NFRVAL*4)
C
C==========================================
C     Determine linear response properties.
C==========================================
C
C----------------------------------------------------------
C     Adjust number of trialvectors for each excitation and
C     write information to output.
C----------------------------------------------------------
C
      IF ( NSAVMX .LT. 2 ) THEN
C
         NSAVMX = 2
C
         WRITE(LUPRI,'(1X,A,/,A,I2,A)') 
     &   'NOTICE: Maximum number of trial vectors for each'//
     &   ' property',' is raised to',NSAVMX, ' as this is'//
     &   ' minimal space allowed.'
C
      END IF
C
      WRITE(LUPRI,'(/,1X,A,I2,/)')
     &'Maximum number of trial vectors for each'//
     &' property is ',NSAVMX
C
C--------------------------------------------------
C     Loop over symmetry of the property operators.
C--------------------------------------------------
C
      DO 200 ISYM = 1, NSYM
C
         DO 100 IOPER = 1, NOPER(ISYM)
C
            LABEL1 = LABOP(IOPER,ISYM)
            IF (LABEL1 .EQ.'XDIPLEN') IDIP = 1
            IF (LABEL1 .EQ.'YDIPLEN') IDIP = 2
            IF (LABEL1 .EQ.'ZDIPLEN') IDIP = 3
C
            ISYMTR = ISYM
C
C====================================
C           For RPA do the following.
C====================================
C
            IF (AORPA) THEN
C
               MODEL = 'AORPA'
C
               IF (IPRSOP .GE. 3) THEN
                  WRITE(LUPRI,9000)
                  WRITE(LUPRI,'(1X,A,2(/A,I8))') 'RPA:',
     &            ' Perturbation symmetry     :',ISYM,
     &            ' p-h + h-p variables.      :',2*NT1AM(ISYM)
                  WRITE(LUPRI,9001)
               END IF
C
C------------------------------------------------------------------
C              Determine RPA solution vectors for the given operator
C              and set of frequencies.
C------------------------------------------------------------------
C
               CALL RP_RSPLEQ(ISYMTR,LABEL1,FRVAL,NFRVAL,WORK(KFOCKD),
     &                        LFOCKD,WORK(KEND1),LWORK1)
C
C----------------------------------------------
C              Calculate second order property.
C----------------------------------------------
C
               CALL SO_POLAR(MODEL,ISYMTR,IDIP,LABEL1,DUMMY,1,DUMMY,1,
     &                       DUMMY,1,
     &                       DUMMY,1,POLDD(1,1,1,1,1),
     &                       POLDQ(1,1,1,1,1,1),POLDL(1,1,1,1,1),
     &                       POLDA(1,1,1,1,1),WORK(KEND1),LWORK1)

C
            END IF
CPFP  2009 Implementation of AOHRP to be done ??
            IF (AOHRP) THEN
C
C
               IF (IPRSOP .GE. 3) THEN
                  WRITE(LUPRI,9000)
                  WRITE(LUPRI,'(1X,A,2(/A,I8))') 'RPA:',
     &            ' Perturbation symmetry     :',ISYM,
     &            ' p-h + h-p variables.      :',2*NT1AM(ISYM)
                  WRITE(LUPRI,9001)
               END IF
C
C------------------------------------------------------------------
C              Determine RPA solution vectors for the given operator
C              and set of frequencies.
C------------------------------------------------------------------
C
c               CALL RP_RSPLEQ(ISYMTR,LABEL1,FRVAL,NFRVAL,WORK(KFOCKD),
c     &                        LFOCKD,WORK(KEND1),LWORK1)
C
C----------------------------------------------
C              Calculate second order property.
C----------------------------------------------
C
c               CALL SO_POLAR(ISYMTR,IDIP,LABEL1,DUMMY,1,DUMMY,1,DUMMY,1,
c     &                       DUMMY,1,POLDD(1,1,1,1,1),
c     &                       POLDQ(1,1,1,1,1,1),POLDL(1,1,1,1,1),
c     &                       POLDA(1,1,1,1,1),WORK(KEND1),LWORK1)
c
C
            END IF
Cend-PFP
C
C======================================
C           For SOPPA do the following.
C======================================
C
            IF (AOSOP) THEN
C
               MODEL = 'AOSOP'
CPFP   2009-Jun
               NVARPT  = 2*(NT1AM(ISYM)+NT2AM(ISYM))
C               NVARPT  = 2*(NT1AM(ISYM)+N2P2HOP(ISYM))
Cend-PFP
               IF (IPRSOP .GE. 3) THEN
                  WRITE(LUPRI,9000)
                  WRITE(LUPRI,'(1X,A)') 'SOPPA:'
                  WRITE(LUPRI,'(3(A,I8,/),A,I8)')
     &            ' Perturbation symmetry     :',ISYM,
     &            ' p-h variables.            :',2*NT1AM(ISYM),
CPFP   2009-Jun
     &            ' 2p-2h variables.          :',2*NT2AM(ISYM),
C     &            ' 2p-2h variables.          :',2*N2P2HOP(ISYM),
Cend-PFP
     &            ' Total number of variables :',NVARPT
                  WRITE(LUPRI,9001)
               END IF
C
C-----------------------------------------------------
C              Get T2 amplitudes and density matrices.
C-----------------------------------------------------
C
               CALL SO_GETT2('MP2 ',WORK(KT2AM),LT2AM,WORK(KDENSAI),
     &                       LDENSAI,WORK(KDENSIJ),LDENSIJ,
     &                       WORK(KDENSAB),LDENSAB,WORK(KEND1),LWORK1)
C
C--------------------------------
C              Initialize DENSAI.
C--------------------------------
C
               CALL DZERO(WORK(KDENSAI),LDENSAI)
C
C------------------------------------------------------------------
C              Determine SOPPA excitation energies and excitation
C              vectors. The excitation vectors are written to file.
C------------------------------------------------------------------
C
               CALL SO_RSPLEQ(MODEL,LABEL1,ISYMTR,FRVAL,NFRVAL,
     &                        WORK(KDENSIJ),LDENSIJ,WORK(KDENSAB),
     &                        LDENSAB,WORK(KDENSAI),LDENSAI,WORK(KT2AM),
     &                        LT2AM,WORK(KFOCKD),LFOCKD,
     &                        WORK(KEND1),LWORK1)
C
C----------------------------------------------
C              Calculate second order property.
C----------------------------------------------
C
               CALL SO_POLAR(MODEL,ISYMTR,IDIP,LABEL1,WORK(KT2AM),LT2AM,
     &                       WORK(KDENSIJ),LDENSIJ,WORK(KDENSAB),
     &                       LDENSAB,WORK(KDENSAI),LDENSAI,
     &                       POLDD(1,1,1,1,2),POLDQ(1,1,1,1,1,2),
     &                       POLDL(1,1,1,1,2),POLDA(1,1,1,1,2),
     &                       WORK(KEND1),LWORK1)
C
            END IF
C     
  100    CONTINUE
C     
  200 CONTINUE
C
C---------------------------------
C     3. allocation of work space.
C---------------------------------
C
      LPARRA = LSOTIM
C
      KPARRA = KEND1
      KEND3  = KPARRA + LPARRA
      LWORK3 = LWORK  - KEND3
C
      CALL SO_MEMMAX ('SO_LNRABA.3     ',LWORK3)
      IF (LWORK3 .LT.0) CALL STOPIT('SO_LNRABA.3',' ',KEND3,LWORK)
C
C---------------------------------------------------
C     Print memory statistics for SOPPA subroutines.
C---------------------------------------------------
C
      CALL SO_MEMMAX('STATISTICS      ',0)
C
C-----------------------------------------
C     Print timings for SOPPA subroutines.
C-----------------------------------------
C
      CALL SO_TIME(TIMTOT,TIMWTO,WORK(KPARRA),LPARRA)
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('SO_LNRABA')
C
      RETURN
C
 9000 FORMAT(/' -----------------------------------')
 9001 FORMAT(' -----------------------------------')
      END
