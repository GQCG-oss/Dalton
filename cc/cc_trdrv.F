C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
#include <single.h>
C
c*DECK CC_TRDRV
      SUBROUTINE CC_TRDRV(FRHO1,LUFR1,FRHO2,LUFR2,
     *                    FC1AM,LUFC1,FC2AM,LUFC2,TRIPLET,
     *                    CC3DIIS,ITRAN,FREQ,
     *                    ISIDE,IST,NSIM,WORK,LWORK)
C
C-------------------------------------------------------------------
C
C     Written by Ove Christiansen 2-11-1996.
C
C     "transformation drive"
C
C     Directs the transformation of trial vecors multiplyed onto the
C     Jacobian from the right and left.
C
C     Performs NSIM linear transformations.
C
C     The trial vectors is taken from files:
C 
C           (name,unit nr.)
C     c1am: (FC1AM,LUFC1)
C     c2am: (FC2AM,LUFC2)
C     Rho1: (FRHO1,LUFR1)
C     Rho2: (FRHO2,LUFR2)
C
C     Triplet transformation by Christof Haettig.
C
C     It is assumed that on disc are the global intermediates
C     For singlet:
C             CC2: E1, E2, F 
C             CCSD & CC3: E1, E2, F, Gamma, BF, C, D intermedias.
C     For triplet:
C             CCSD : E1,E2,F,Gamma,BF,C,D,CD intermediates
C
C     Local intermediates are put to scratch as well: 
C
C     For right hand transformation: 
C             C-tilde,D-tilde intermediates.
C
C     For left hand side:
C
C             O3V integrals.
C
C     For Triplet : (1)C-tilde, (3)C-tilde, (1)D-tilde
C                   (P)D-tilde, (M)D-tilde.
C
C
C     The variable MINSCR controls where the vectors are to be
C     transformed one by one or all in one integral calculation.
C
C     Variables IVEC is the start vector number to be transformed
C     and ITR is the vector number for the result vector.
C     (Used slightly different in CCS and CC2 relative to CCSD,..
C      since in these cases)
C
C     Biorthonormal basis: 
C
C     Singlet:
C     (L: 1/2*Eia,(2*EiaEjb + EjaEib)/[6*(1 + delta(ab)delta(ij))] )
C     (R: Eai, EaiEbj)
C
C     Triplet: See CC_RHTR3.F
C
C     Version 1.0 2-11-1996 (base on cclr_trr)
C
C-------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <maxash.h>
#include <mxorb.h>
C
#include <ccorb.h>
#include <ccisao.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <cclr.h>
#include <ccsdio.h>
#include <leinf.h>
#include <aovec.h>
      PARAMETER ( TWO = 2.0D00,XHALF=0.5D00 )
      DIMENSION WORK(LWORK), FREQ(*)
C
      LOGICAL ORSAVE,T2TSAV,MSCRS,TRIPLET,CC3DIIS
      CHARACTER*8  FRHO1,FRHO2,FC1AM,FC2AM,FR2SD
      INTEGER ITRAN(NSIM)
C
C-------------------------
C     Save and initialize.
C-------------------------
C
      CALL QENTER('CC_TRDRV')
      MSCRS = MINSCR
      IF ( CCSDT .OR. (ISIDE .EQ. 2)) THEN
         MINSCR = .TRUE.
      ENDIF  
C
      IF (IPRINT .GT. 5) THEN
         CALL AROUND(' START OF CC_TRDRV ')
         IF (DIRECT) WRITE(LUPRI,*) ' Atomic direct calculation'
         WRITE(LUPRI,*) ' Workspace input: ',LWORK
         WRITE(LUPRI,*) ' Triplet flag   : ',TRIPLET
         WRITE(LUPRI,'(/A48,I5)')
     *   ' CC_TRDRV: The total number of trial vectors is: ',NSIM
      ENDIF
C
      T2TSAV = T2TCOR
      IF (CCS .OR. CC2) T2TCOR = .FALSE.
      IF (CC2 .AND.(ISIDE .EQ. 2)) T2TCOR = T2TSAV
      IF (TRIPLET) T2TCOR = .FALSE.
C
      NCAMP = NT1AM(ISYMTR) + NT2AM(ISYMTR)
      IF (TRIPLET) NCAMP = NCAMP + NT2AMA(ISYMTR)
      IF (CCS) NCAMP = NT1AM(ISYMTR)
C
      NC2 = 1
      IF ( CCS ) NC2 = 0
C
      IF (DIRECT .AND. JACEXP ) CALL QUIT(
     * 'CC_TRDRV: JACEXP test flag does not work in DIRECT calc.')
C
C------------------------------------------------------------------
C     Make rho2 file name.
C     For CCSD rho2 has to be stored on different file due to 
C     different length.
C------------------------------------------------------------------
C
      IF (.NOT. (CCS.OR.CC2)) THEN
         LUFSD = 0
         FR2SD = 'CC_TRA2_'
      ELSE
         LUFSD = LUFR2
         FR2SD = FRHO2
      ENDIF
C
C------------------------------------------------------------------
C     Cheat and do CCS left transformation by right hand 
C     transformation.
C------------------------------------------------------------------
C
      NSIDSA = ISIDE 
      IF ((ISIDE .EQ. -1 ) .AND. CCS ) ISIDE = 1
C
C--------------------------------------------------------------------
C     Make total transformation.
C     Orthonormal basis with omegor packed rho in 
C     delta loop in ccsd.
C     NB: It is assumed that also omegor for intermediate BF
C         for ccsd.
C--------------------------------------------------------------------
C
      ONLY21 = .FALSE.
C
      ORSAVE = OMEGOR
      IF ( CCS .OR. CC2) THEN
         OMEGOR = .FALSE.
      ELSE
         OMEGOR = .TRUE.
      ENDIF
 
C
C------------------------------------------------------
C     Transform vectors:  ONE or ALL.
C     If minscr then one at a time, else all vectors in one 
C     integral calculation. 
C 
C     Keep C1 in core and C2 on disk.
C     Keep rho1 in core and rho2 on disk. ?????????????????????
C------------------------------------------------------
C
      IF ( MINSCR ) THEN
         NSIMTR = 1
      ELSE
C        NSIMTR = NSIM
         NSIMTR = MAX(NSIM,1)
      ENDIF
C
      IF (NSIMTR .GT. 100 ) CALL QUIT(
     *   ' Maximum nr of simultaneously trial vectors')
C
      IF (CC3DIIS .AND. (.NOT.MINSCR))
     *  CALL QUIT('CC3DIIS option requires MINSCR=.TRUE.')
C
C--------------------------------------------------------
C     Max is 100 due to limitations in the array IT2DLR
C     set in ccsd.cdk.
C--------------------------------------------------------
C
      NL  = (NSIM -1 )/NSIMTR + 1
C
      DO 100 I = 1, NL
C
         IF ( .NOT.MINSCR) THEN
C
            K1 = IST
C
            IF (IPRINT .GT. 5 ) THEN
               WRITE(LUPRI,'(A24,I3,A37)') ' CC_TRDRV: Transforming ',
     *         NSIMTR,' vectors in one call to CC_RHTR.'
            ENDIF
C
         ELSE
C
            IF (CC3DIIS) THEN
              K1    = ITRAN(I)
              ECURR = FREQ(K1)
            ELSE
              K1 = I + IST - 1
            END IF
            
C
            IF (IPRINT .GT. 5 ) THEN
               WRITE(LUPRI,'(A35,I4)')
     *         ' CC_TRDRV: Transforming vector nr.:',K1
               IF (CCSDT) WRITE(LUPRI,'(A,F8.4)')
     *         ' CC_TRDRV: Frequency in R3 denominators:',ECURR
            ENDIF
C
         ENDIF
C
         CALL FLSHFO(LUPRI)
C
C--------------------
C        Allocations.
C--------------------
C
         IF (.NOT. TRIPLET) THEN
            NRHO2 = MAX(NT2AM(ISYMTR),NT2AM(1),2*NT2ORT(ISYMTR))
            IF (ISIDE .EQ. -1) NRHO2 = MAX(NRHO2,2*NT2ORT(1))
            IF ( CC2 ) NRHO2 = MAX(NT2AM(ISYMTR),NT2AM(1))
            IF ( CCS ) NRHO2 = 2
C
            NC2AM = MAX(NT2SQ(ISYMTR),NT2SQ(1),
     *              NT2AM(ISYMTR)+2*NT2ORT(1))
            IF ( CC2 ) NC2AM = MAX(NT2SQ(ISYMTR),NT2SQ(1))
            IF ( CCS ) NC2AM = 2
C
            NRHO1 = NT1AM(ISYMTR)*NSIMTR
C
            KRHO1 = 1    
            KRHO2 = KRHO1 + NRHO1
            KC1AM = KRHO2 + NRHO2
            KC2AM = KC1AM + NT1AM(ISYMTR)*NSIMTR
            KEND1 = KC2AM + NC2AM
            LWRK1 = LWORK - KEND1             
            IF (LWRK1 .LE. 0 )
     *             CALL QUIT('Too little workspace in cclr_trr')
C
C--------------------------------------------
C        Read the CC trial vectors from disk.
C--------------------------------------------
C
            DO 150 IV = 1, NSIMTR
               KOFF1  = KC1AM + NT1AM(ISYMTR)*(IV - 1)
               CALL CC_RVEC(LUFC1,FC1AM,NT1AM(ISYMTR),NT1AM(ISYMTR),
     *                      IV+K1-1,WORK(KOFF1))
  150       CONTINUE
C
            IF ((.NOT.CCS).AND.( MINSCR)) THEN
               CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                      K1,WORK(KRHO2))
            ENDIF
C
C---------------------------
C        Print C amplitudes.
C---------------------------
C
            IF (IPRINT .GT. 45) THEN
               DO 151 IV = 1, NSIMTR
                  KOFF1  = KC1AM + NT1AM(ISYMTR)*(IV - 1)
                  CALL AROUND('CC_TRDRV: (C1,C2) packed ')
                  WRITE(LUPRI,*) 'Vector nr is = ',IV+K1-1
                  CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,1,NC2)
  151          CONTINUE
            ENDIF
C
C----------------------------------------------------------
C        Test with finited difference CC jacobian if FDJAC.
C----------------------------------------------------------
C
            IF (FDJAC .OR. JACEXP) THEN
C
               CALL DCOPY(NT2AM(ISYMTR),WORK(KRHO2),1,WORK(KC2AM),1)
               CALL CCLR_DUMTRR(WORK(KC1AM),WORK(KRHO1),WORK(KRHO2),
     *                          WORK(KEND1),LWRK1)
C
C-----------------------------------------
C           Print the transformed vectors.
C-----------------------------------------
C
               IF (IPRINT .GT. 50) THEN
                  CALL AROUND('CC_TRDRV: RHO  trans. by DUMTRR .')
                  CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYMTR,1,1)
               ENDIF
C
            ELSE
C
C---------------------------------------------------
C           Prepare the C-amplitudes.
C---------------------------------------------------
C
               IF ( .NOT. CCS ) THEN
                  IF ( ISIDE .GE. 1) THEN
                     CALL CCLR_DIASCL(WORK(KRHO2),TWO,ISYMTR)
                  ENDIF
                  CALL CC_T2SQ(WORK(KRHO2),WORK(KC2AM),ISYMTR)
               ENDIF
C
C-------------------------------------
C           Print the squared vectors.
C-------------------------------------
C
               IF (IPRINT.GT.50) THEN
                  CALL AROUND('CC_TRDRV: (C1,C2) squared ')
                  CALL CC_PRSQ(WORK(KC1AM),WORK(KC2AM),ISYMTR,1,NC2)
               ENDIF
C
C---------------------------
C           Zero rho vector.
C---------------------------
C
               CALL DZERO(WORK(KRHO1),NT1AM(ISYMTR)*NSIMTR)
               CALL DZERO(WORK(KRHO2),NRHO2)
C
C------------------------------
C           File read and save.
C------------------------------
C
               IF (.NOT. (CCS.OR.CC2)) THEN 
                  CALL WOPEN2(LUFSD,FR2SD,64,0)
               ENDIF
C
               NRHO2 = MAX(NT2AM(ISYMTR),2*NT2ORT(ISYMTR))
               IF (CC2 ) NRHO2 = NT2AM(ISYMTR)
               DO 80 IV = 1, NSIMTR
                  NR1 = IV + K1 - 1
                  IF (.NOT. (CCS.OR.CC2)) THEN
                     NR2 = IV 
                  ELSE
                     NR2 = NR1
                  ENDIF
                  CALL CC_WVEC(LUFR1,FRHO1,NT1AM(ISYMTR),
     *                         NT1AM(ISYMTR),NR1,WORK(KRHO1))
                  IF (.NOT.CCS) THEN
                     CALL CC_WVEC(LUFSD,FR2SD,NRHO2,NRHO2,NR2,
     *                            WORK(KRHO2))
                  ENDIF
  80           CONTINUE
C
            ENDIF
C
         ELSE
C
C----------------------------------------
C       For Triplet start at the
C       beginning of the workspace.
C----------------------------------------
C
            NRHO1 = NT1AM(ISYMTR)
            NRHO2 = NT2AM(ISYMTR)+NT2AMA(ISYMTR)
C
            KRHO1 = 1    
            KRHO2 = KRHO1 + NRHO1
            KEND1 = KRHO2 + NRHO2
            LWRK1 = LWORK - KEND1             
            IF (LWRK1 .LE. 0) THEN
                CALL QUIT('Too little workspace in CC_TRDRV ')
            ENDIF
C
C------------------------------------------------
C           File opening for the triplet case
C------------------------------------------------
C
             IF (.NOT. (CCS.OR.CC2)) THEN 
               CALL WOPEN2(LUFSD,FR2SD,64,0)
             ENDIF
C
         ENDIF
C
C----------------------------------------
C           Calculate transformed vectors.
C-----------------------------------------
C
         IF (.NOT. (FDJAC .OR. JACEXP)) THEN
            IVEC = K1
            IF ( .NOT.(CCS.OR.CC2) ) THEN
               ITR  = 1
            ELSE
               ITR  = K1
            ENDIF
C
            LRHO1 = NT1AM(ISYMTR)
            IF (TRIPLET) THEN
               CALL CC_RHTR3(FRHO1,LUFR1,FR2SD,LUFSD,
     *                      FC1AM,LUFC1,FC2AM,LUFC2,
     *                      WORK,LWORK,NSIMTR,
     *                      IVEC,ITR)
C
            ELSE
               IF (ISIDE .EQ. 1) THEN
                  CALL CC_RHTR(FRHO1,LUFR1,FR2SD,LUFSD,
     *                         FC1AM,LUFC1,FC2AM,LUFC2,
     *                         WORK(KRHO1),WORK(KRHO2),
     *                         WORK(KC1AM),WORK(KC2AM),
     *                         WORK(KEND1),LWRK1,NSIMTR,
     *                         IVEC,ITR,LRHO1)
               ELSE IF (ISIDE .EQ. -1) THEN
                  CALL CC_LHTR(FRHO1,LUFR1,FR2SD,LUFSD,
     *                         FC1AM,LUFC1,FC2AM,LUFC2,
     *                         WORK(KRHO1),WORK(KRHO2),
     *                         WORK(KC1AM),WORK(KC2AM),
     *                         WORK(KEND1),LWRK1,NSIMTR,
     *                         IVEC,ITR,LRHO1)
               ELSE
                  CALL QUIT(' ISIDE should be -1 or +1 ')
               ENDIF
C
            ENDIF
C
C-----------------------------------------
C           Print the transformed vectors.
C-----------------------------------------
C
            IF (( IPRINT .GT. 15).OR.(.NOT.(CCS.OR.CC2))) THEN
              IF (TRIPLET) THEN
                 NRHO2 = NT2AM(ISYMTR) + NT2AMA(ISYMTR)
              ELSE
                 NRHO2 = MAX(NT2AM(ISYMTR),2*NT2ORT(ISYMTR))
                 IF (CC2 ) NRHO2 = NT2AM(ISYMTR)
              END IF
              DO 90 IV = 1, NSIMTR
                 NR1 = IV + K1 - 1
                 IF (.NOT. (CCS.OR.CC2)) THEN
                    NR2 = IV 
                 ELSE
                    NR2 = NR1
                 ENDIF
                 CALL CC_RVEC(LUFR1,FRHO1,NT1AM(ISYMTR),NT1AM(ISYMTR),
     *                        NR1,WORK(KRHO1))
                 IF (.NOT.CCS) THEN
                   CALL CC_RVEC(LUFSD,FR2SD,NRHO2,NRHO2,NR2,WORK(KRHO2))
                 END IF
                 IF (IPRINT .GT. 45) THEN
                    CALL AROUND('CC_TRDRV: RHO = trans. Vector ')
                    WRITE (LUPRI,*) 'number of vector on file:',NR2
                    IF (TRIPLET.AND. (.NOT.CCS)) NC2 = 2
                    CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYMTR,1,NC2)
                 ENDIF
!
                 IF (.NOT.(CCS.OR.CC2) .AND. (.NOT. TRIPLET)) THEN
                    CALL CC_WVEC(LUFR2,FRHO2,NT2AM(ISYMTR),
     *                           NT2AM(ISYMTR),NR1,WORK(KRHO2))
                 ENDIF
!
                 IF ((TRIPLET) .AND. (.NOT. (CCS.OR.CC2))) THEN
                    CALL CC_WVEC3(LUFR2,FRHO2,
     *                            NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                            NT2AM(ISYMTR)+NT2AMA(ISYMTR),
     *                            NR1,0,WORK(KRHO2))
                 ENDIF
  90          CONTINUE
            ENDIF
C
C-----------------------
C           Close files.
C-----------------------
C
            IF ( .NOT.(CCS.OR.CC2)) THEN
               CALL WCLOSE2(LUFSD,FR2SD,'DELETE')
            ENDIF
C
         ENDIF
C
 100  CONTINUE
C
C-------------
C     Restore.
C-------------
C
      T2TCOR = T2TSAV
      ISIDE  = NSIDSA
      OMEGOR = ORSAVE
      MINSCR = MSCRS
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND(' END OF CC_TRDRV ')
      ENDIF
C
   1  FORMAT(1x,A35,1X,E20.10)
      CALL QEXIT('CC_TRDRV')
      RETURN
      END
