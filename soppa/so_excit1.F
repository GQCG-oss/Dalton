#include <single.h>
C
C  /* Deck so_excit1 */
CSPAS:23/5-11: second and third moment sum rules
C     SUBROUTINE SO_EXCIT1(TRLEN,TRVEL,TRLON,TRMAG,BSRLON,SECMAT,EXENG,
C    &                     WORK,LWORK)
      SUBROUTINE SO_EXCIT1(TRLEN,TRVEL,TQLEN,TQVEL,TTLEN,TRLON,TRMAG,
     &                     BSRLON,SECMAT,EXENG,WORK,LWORK)
CKeinSPASmehr
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1997
C     SPAS :4/12-1997 : call to ccsd_init1 include. This should allow
C                       to run the RPA part without going through the
C                       CC program.
C
C     PURPOSE: Main driver routine for the calculation of excitation 
C              energies and transition moments with the atomic integral
C              direct SOPPA program.
C
#include <implicit.h>
#include <priunit.h>
C
CSPAS:15/03-2006: merge with Dalton-2.0  
C#include <cbilrs.h>
CKeinSPASmehr
#include <cbiexc.h>
#include <ccsdsym.h>
#include <infdim.h>
#include <inforb.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
CSPAS:10/11-2003: is merged with maxorb.h
C#include <mxorb.h>
CKeinSPASmehr
#include <nuclei.h>
CSPAS:10/11-2003: I cannot find the corresponding include file
C#include <nuctap.h>
CKeinSPASmehr
#include <soppinf.h>
#include <symmet.h>
#include <wrkrsp.h>
C
      PARAMETER (HALF = 0.5D0,ESUDIP = 64604.885D0,ESUECD = 471.44360D0)
C
      DIMENSION TRLEN(3,NSYM,MXNEXI,4),  TRVEL(3,NSYM,MXNEXI,4)
CSPAS:23/5-11: second and third moment sum rules
      DIMENSION TQLEN(3,3,NSYM,MXNEXI,4),  TQVEL(3,3,NSYM,MXNEXI,4)
      DIMENSION TTLEN(10,NSYM,MXNEXI,4)
CKeinSPASmehr
      DIMENSION TRLON(3,NSYM,MXNEXI,4),  TRMAG(3,NSYM,MXNEXI,4)
      DIMENSION BSRLON(3,NSYM,MXNEXI,4), EXENG(NSYM,MXNEXI,4)
      DIMENSION SECMAT(3,MXNEXI,NSYM)
      DIMENSION WORK(LWORK)
C
      CHARACTER MODEL*5
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_EXCIT1')
C
C-------------------------------------------
C     Start timing the AO-SOPPA calculation.
C-------------------------------------------
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
      IPRSOP = IPREXC
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
      CALL SO_MEMMAX ('SO_EXCIT1.1',LWORK1)
      IF (LWORK1 .LT .0) CALL STOPIT('SO_EXCIT1.1',' ',KEND1,LWORK)
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
      CALL SO_PRPINT('EXCITA',NLBTOT,WORK(KEND1),LWORK1)
C
C---------------------------------------------------------------------
C     Initialize arrays of transition moments and excitation energies.
C---------------------------------------------------------------------
C
      CALL DZERO(TRLEN, 3*NSYM*MXNEXI*4)
      CALL DZERO(TRVEL, 3*NSYM*MXNEXI*4)
CSPAS:23/5-11: second and third moment sum rules
      CALL DZERO(TQLEN, 3*3*NSYM*MXNEXI*4)
      CALL DZERO(TQVEL, 3*3*NSYM*MXNEXI*4)
      CALL DZERO(TTLEN, 10*NSYM*MXNEXI*4)
CKeinSPASmehr
      CALL DZERO(TRLON, 3*NSYM*MXNEXI*4)
      CALL DZERO(TRMAG, 3*NSYM*MXNEXI*4)
      CALL DZERO(BSRLON,3*NSYM*MXNEXI*4)
      CALL DZERO(SECMAT,3*NSYM*MXNEXI)
      CALL DZERO(EXENG, NSYM*MXNEXI*4)
C
C========================================================
C     Determine excitation energies, response vectors and 
C     linear transition moments.
C========================================================
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
     &   ' excitation',' is raised to',NSAVMX, ' as this is'//
     &   ' minimal space allowed.'
C
      END IF
C
      WRITE(LUPRI,'(/,1X,A,I2,/)')
     &'Maximum number of trial vectors for each'//
     &' excitation is ',NSAVMX
C
C-------------------------------------------------------------
C     Loop over symmetry of the excitation vector.
C     NOTE: The excited state has the symmtry of the reference
C     state times the symmetry of the excitation vector.
C-------------------------------------------------------------
C
      DO 100 ISYM = 1, NSYM
C
         CALL RSPSET
C
         IF (NEXCIT(ISYM) .GT. 0) THEN
C
C---------------------------------------
C           2. allocation of work space.
C---------------------------------------
C
            KEXVAL = KEND1
            KEND2  = KEXVAL  + NEXCIT(ISYM)
            LWORK2 = LWORK   - KEND2
C
            CALL SO_MEMMAX ('SO_EXCIT1.2',LWORK2)
            IF (LWORK2 .LT.0) CALL STOPIT('SO_EXCIT1.2',' ',KEND2,LWORK)
C
            ISYMTR = ISYM
            LEXVAL = NEXCIT(ISYM)
C
C====================================
C           For RPA do the following.
C====================================
C
            IF (AORPA) THEN
C
               NVARPT = 2*NT1AM(ISYM)
C
               IF (IPRSOP .GE. 3) THEN
                  WRITE(LUPRI,9000)
                  WRITE(LUPRI,'(1X,A,2(/A,I8))') 'RPA:',
     &            ' Perturbation symmetry     :',ISYM,
     &            ' p-h + h-p variables.      :',2*NT1AM(ISYM)
                  WRITE(LUPRI,9001)
               END IF
C
               NEXCI  = MIN(NEXCIT(ISYM),NT1AM(ISYM))
C
               IF ( NEXCI .LT. NEXCIT(ISYM) ) THEN
                  WRITE(LUPRI,'(1X,A,I5,/,1X,A)') 
     &            'NOTICE: The number of excitations are reduced to '
     &            ,NEXCI,'as this is the full dimension of the '//
     &            'excitation space'
               END IF
C
C------------------------------------------------------------------
C              Determine RPA excitation energies and excitation
C              vectors. The excitation vectors are written to file.
C------------------------------------------------------------------
C
               CALL RP_RSPLEX(ISYMTR,NEXCI,WORK(KEXVAL),LEXVAL,
     &                        WORK(KFOCKD),LFOCKD,WORK(KEND2),LWORK2)
C
C----------------------------------------------------
C              Copy RPA excitation energies to EXENG.
C----------------------------------------------------
C
               DO 10 IEXVAL = 1,NEXCI
                  IF (IEXVAL .LE. MXNEXI)
     &               EXENG(ISYM,IEXVAL,1) = WORK(KEXVAL-1+IEXVAL)
                  IF (IEXVAL .GT. MXNEXI)
     &               WRITE(LUPRI,*) 
     &                    'WARNING: IEXVAL greater than MXNEXI.'
   10          CONTINUE
C
C-----------------------------------------------
C              Determine RPA transition moments.
C-----------------------------------------------
C
               CALL RP_TRMOM(ISYMTR,NLBTOT,NEXCI,TRLEN(1,1,1,1),
CSPAS:23/5-11: second and third moment sum rules
C    &                       TRVEL(1,1,1,1),TRLON(1,1,1,1),
C    &                       TRMAG(1,1,1,1),BSRLON(1,1,1,1),
C    &                       EXENG(1,1,1),WORK(KEND2),LWORK2)
     &                       TRVEL(1,1,1,1),TQLEN(1,1,1,1,1),
     &                       TQVEL(1,1,1,1,1),TTLEN(1,1,1,1),
     &                       TRLON(1,1,1,1),TRMAG(1,1,1,1),
     &                       BSRLON(1,1,1,1),EXENG(1,1,1),
     &                       WORK(KEND2),LWORK2)
CKeinSPASmehr
C
C--------------------------------------------------------
C              Determine charge radius of excited states.
C--------------------------------------------------------
C
               CALL RP_CHARGE(ISYMTR,NEXCI,SECMAT,WORK(KEND2),LWORK2)
C
            END IF
C
C=======================================
C           For RPA(D) do the following.
C=======================================
C
            IF (DCRPA) THEN
C
               MODEL = 'DCRPA'
C
               NVARPT = 2*(NT1AM(ISYM)+NT2AM(ISYM))
C
               IF (IPRSOP .GE. 3) THEN
                  WRITE(LUPRI,9000)
                  WRITE(LUPRI,'(1X,A,4(/A,I8))') 'RPA(D):',
     &            ' Perturbation symmetry     :',ISYM,
     &            ' p-h variables.            :',2*NT1AM(ISYM),
     &            ' 2p-2h variables.          :',2*NT2AM(ISYM),
     &            ' Total number of variables :',NVARPT
                  WRITE(LUPRI,9001)
               END IF
C
               NEXCI  = MIN(NEXCIT(ISYM),NT1AM(ISYM))
C
               IF ( NEXCI .LT. NEXCIT(ISYM) ) THEN
                  WRITE(LUPRI,'(1X,A,I5,/,1X,A)') 
     &            'NOTICE: The number of excitations are reduced to '
     &            ,NEXCI,'as this is the full dimension of the '//
     &            'RPA excitation space'
               END IF
C
C-----------------------------------------------------
C              Get T2 amplitudes and density matrices.
C-----------------------------------------------------
C
               CALL SO_GETT2('MP2 ',WORK(KT2AM),LT2AM,WORK(KDENSAI),
     &                       LDENSAI,WORK(KDENSIJ),LDENSIJ,
     &                       WORK(KDENSAB),LDENSAB,WORK(KEND2),LWORK2)
C
C--------------------------------
C              Initialize DENSAI.
C--------------------------------
C
               CALL DZERO(WORK(KDENSAI),LDENSAI)
C
C------------------------------------------------------------------
C              Determine RPA(D) excitation energies and excitation 
C              vectors. The excitation vectors are written to file.
C------------------------------------------------------------------
C
               CALL DC_RSPLEX(ISYMTR,NEXCI,WORK(KEXVAL),LEXVAL,
     &                        WORK(KDENSIJ),LDENSIJ,WORK(KDENSAB),
     &                        LDENSAB,WORK(KDENSAI),LDENSAI,WORK(KT2AM),
     &                        LT2AM,WORK(KFOCKD),LFOCKD,
     &                        WORK(KEND2),LWORK2)
C
C-------------------------------------------------------
C              Copy RPA(D) excitation energies to EXENG.
C-------------------------------------------------------
C
               DO 20 IEXVAL = 1,NEXCI
                  IF (IEXVAL .LE. MXNEXI)
     &               EXENG(ISYM,IEXVAL,2) = WORK(KEXVAL-1+IEXVAL)
                  IF (IEXVAL .GT. MXNEXI)
     &               WRITE(LUPRI,*) 
     &                    'WARNING: IEXVAL greater than MXNEXI.'
   20          CONTINUE
C
C--------------------------------------------------
C              Determine RPA(D) transition moments.
C--------------------------------------------------
C
               CALL SO_TRMOM(MODEL,ISYMTR,NLBTOT,NEXCI,WORK(KT2AM),
     &                       LT2AM,WORK(KDENSIJ),LDENSIJ,WORK(KDENSAB),
     &                       LDENSAB,WORK(KDENSAI),LDENSAI,
     &                       TRLEN(1,1,1,2),TRVEL(1,1,1,2),
CSPAS:23/5-11: second and third moment sum rules
     &                       TQLEN(1,1,1,1,2),TQVEL(1,1,1,1,2),
     &                       TTLEN(1,1,1,2),
CKeinSPASmehr
     &                       TRLON(1,1,1,2),TRMAG(1,1,1,2),
     &                       BSRLON(1,1,1,2),EXENG(1,1,2),
     &                       WORK(KEND2),LWORK2)
C
            END IF
C
C======================================
C           For SOPPA do the following.
C======================================
C
            IF (AOSOP) THEN
C
               MODEL = 'AOSOP'
C
               NVARPT  = 2*(NT1AM(ISYM)+NT2AM(ISYM))
C
               IF (IPRSOP .GE. 3) THEN
                  WRITE(LUPRI,9000)
                  WRITE(LUPRI,'(1X,A)') 'SOPPA:'
                  WRITE(LUPRI,'(3(A,I8,/),A,I8)')
     &            ' Perturbation symmetry     :',ISYM,
     &            ' p-h variables.            :',2*NT1AM(ISYM),
     &            ' 2p-2h variables.          :',2*NT2AM(ISYM),
     &            ' Total number of variables :',NVARPT
                  WRITE(LUPRI,9001)
               END IF
C
               NEXCI = MIN(NEXCIT(ISYM),NVARPT/2)
C
               IF ( NEXCI .LT. NEXCIT(ISYM) ) THEN
                  WRITE(LUPRI,'(1X,A,I5,/,1X,A)') 
     &            'NOTICE: The number of excitations are reduced to '
     &            ,NEXCI,'as this is the full dimension of the '//
     &            'excitation space'
               END IF
C
C-----------------------------------------------------
C              Get T2 amplitudes and density matrices.
C-----------------------------------------------------
C
               CALL SO_GETT2('MP2 ',WORK(KT2AM),LT2AM,WORK(KDENSAI),
     &                       LDENSAI,WORK(KDENSIJ),LDENSIJ,
     &                       WORK(KDENSAB),LDENSAB,WORK(KEND2),LWORK2)
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
               CALL SO_RSPLEX(MODEL,ISYMTR,NEXCI,WORK(KEXVAL),LEXVAL,
     &                        WORK(KDENSIJ),LDENSIJ,WORK(KDENSAB),
     &                        LDENSAB,WORK(KDENSAI),LDENSAI,WORK(KT2AM),
     &                        LT2AM,WORK(KFOCKD),LFOCKD,
     &                        WORK(KEND2),LWORK2)
C
C-----------------------------------------------------
C              Copy SOPPA exciation energies to EXENG.
C-----------------------------------------------------
C
               DO 30 IEXVAL = 1,NEXCI
                  IF (IEXVAL .LE. MXNEXI)
     &               EXENG(ISYM,IEXVAL,3) = WORK(KEXVAL-1+IEXVAL)
                  IF (IEXVAL .GT. MXNEXI)
     &               WRITE(LUPRI,*) 
     &                    'WARNING: IEXVAL greater than MXNEXI.'
   30          CONTINUE
C
C-------------------------------------------------
C              Determine SOPPA transition moments.
C-------------------------------------------------
C
               CALL SO_TRMOM(MODEL,ISYMTR,NLBTOT,NEXCI,WORK(KT2AM),
     &                       LT2AM,WORK(KDENSIJ),LDENSIJ,WORK(KDENSAB),
     &                       LDENSAB,WORK(KDENSAI),LDENSAI,
     &                       TRLEN(1,1,1,3),TRVEL(1,1,1,3),
CSPAS:23/5-11: second and third moment sum rules
     &                       TQLEN(1,1,1,1,3),TQVEL(1,1,1,1,3),
     &                       TTLEN(1,1,1,3),
CKeinSPASmehr
     &                       TRLON(1,1,1,3),TRMAG(1,1,1,3),
     &                       BSRLON(1,1,1,3),EXENG(1,1,3),
     &                       WORK(KEND2),LWORK2)
C
            END IF
C
C============================================
C           For SOPPA(CCSD) do the following.
C============================================
C
            IF (AOSOC) THEN
C
               MODEL = 'AOSOC'
C
               NVARPT  = 2*(NT1AM(ISYM)+NT2AM(ISYM))
C
               IF (IPRSOP .GE. 3) THEN
                  WRITE(LUPRI,9000)
                  WRITE(LUPRI,'(1X,A)') 'SOPPA(CCSD):'
                  WRITE(LUPRI,'(3(A,I8,/),A,I8)')
     &            ' Perturbation symmetry     :',ISYM,
     &            ' p-h variables.            :',2*NT1AM(ISYM),
     &            ' 2p-2h variables.          :',2*NT2AM(ISYM),
     &            ' Total number of variables :',NVARPT
                  WRITE(LUPRI,9001)
               END IF
C
               NEXCI = MIN(NEXCIT(ISYM),NVARPT/2)
C
               IF ( NEXCI .LT. NEXCIT(ISYM) ) THEN
                  WRITE(LUPRI,'(1X,A,I5,/,1X,A)')
     &            'NOTICE: The number of excitations are reduced to '
     &            ,NEXCI,'as this is the full dimension of the '//
     &            'excitation space'
               END IF
C
C-----------------------------------------------------
C              Get T2 amplitudes and density matrices.
C-----------------------------------------------------
C
               CALL SO_GETT2('CCSD',WORK(KT2AM),LT2AM,WORK(KDENSAI),
     &                       LDENSAI,WORK(KDENSIJ),LDENSIJ,
     &                       WORK(KDENSAB),LDENSAB,WORK(KEND2),LWORK2)
C
C----------------------------------------------------------------------
C              Determine SOPPA(CCSD) excitation energies and excitation
C              vectors. The excitation vectors are written to file.
C----------------------------------------------------------------------
C
               CALL SO_RSPLEX(MODEL,ISYMTR,NEXCI,WORK(KEXVAL),LEXVAL,
     &                        WORK(KDENSIJ),LDENSIJ,WORK(KDENSAB),
     &                        LDENSAB,WORK(KDENSAI),LDENSAI,WORK(KT2AM),
     &                        LT2AM,WORK(KFOCKD),LFOCKD,
     &                        WORK(KEND2),LWORK2)
C
C-----------------------------------------------------------
C              Copy SOPPA(CCSD) exciation energies to EXENG.
C-----------------------------------------------------------
C
               DO 40 IEXVAL = 1,NEXCI
                  IF (IEXVAL .LE. MXNEXI)
     &               EXENG(ISYM,IEXVAL,4) = WORK(KEXVAL-1+IEXVAL)
                  IF (IEXVAL .GT. MXNEXI)
     &               WRITE(LUPRI,*) 
     &                    'WARNING: IEXVAL greater than MXNEXI.'
   40          CONTINUE
C
C-------------------------------------------------------
C              Determine SOPPA(CCSD) transition moments.
C-------------------------------------------------------
C
               CALL SO_TRMOM(MODEL,ISYMTR,NLBTOT,NEXCI,WORK(KT2AM),
     &                       LT2AM,WORK(KDENSIJ),LDENSIJ,WORK(KDENSAB),
     &                       LDENSAB,WORK(KDENSAI),LDENSAI,
     &                       TRLEN(1,1,1,4),TRVEL(1,1,1,4),
CSPAS:23/5-11: second and third moment sum rules
     &                       TQLEN(1,1,1,1,4),TQVEL(1,1,1,1,4),
     &                       TTLEN(1,1,1,4),
CKeinSPASmehr
     &                       TRLON(1,1,1,4),TRMAG(1,1,1,4),
     &                       BSRLON(1,1,1,4),EXENG(1,1,4),
     &                       WORK(KEND2),LWORK2)
C
            END IF
C
         END IF
C
  100 CONTINUE
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
      CALL SO_MEMMAX ('SO_EXCIT1.3     ',LWORK3)
      IF (LWORK3 .LT.0) CALL STOPIT('SO_EXCIT1.3',' ',KEND3,LWORK)
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
      CALL QEXIT('SO_EXCIT1')
C
      RETURN
C
 9000 FORMAT(/' -----------------------------------')
 9001 FORMAT(' -----------------------------------')
      END
