C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
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
      SUBROUTINE CC_TRDRV(ECURR,FRHO1,LUFR1,FRHO2,LUFR2,FRHO12,LUFR12,
     *                          FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                    TRIPLET,CC3DIIS,ITRAN,FREQ,FS12AM,LUFS12,
     *                    ISIDE,IST,NSIM,WORK,LWORK)
C
C-------------------------------------------------------------------
C
C     "transformation drive"
C
C     Directs the transformation of trial vecors multiplyed onto the
C     Jacobian from the right and left.
C
C     Performs NSIM linear transformations.
C
C     The trial vectors are taken from files:
C 
C             (name,unit nr.)
C     c1am:   (FC1AM,LUFC1)     singles
C     c2am:   (FC2AM,LUFC2)     doubles
C     cr12am: (FC12AM,LUFC12)   R12 doubles
C     Rho1:   (FRHO1,LUFR1)     singles
C     Rho2:   (FRHO2,LUFR2)     doubles
C     RhoR12: (FRHO12,LUFR12)   R12 doubles
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
C     The variable MINSCR controls how the vectors are to be
C     transformed: one by one or all in one integral calculation.
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
C     Written by Ove Christiansen 2-11-1996.
C     Changes for triplet transformation by Christof Haettig, 1999.
C     Changes for CC-R12 by Christof Haettig, Jun 2003.
C
C-------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <maxash.h>
#include <maxorb.h>
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
      CHARACTER*8  FRHO1,FRHO2,FC1AM,FC2AM,FR2SD,FRHO12,FC12AM,FS12AM
      INTEGER ITRAN(NSIM)
C
C-------------------------
C     Save and initialize.
C-------------------------
C
      CALL QENTER('CC_TRDRV')
C
      MSCRS = MINSCR
      IF ( CCSDT .OR. (ISIDE.EQ.2) .OR. FDJAC .OR. JACEXP ) THEN
         MINSCR = .TRUE.
      ENDIF  
C
      IF (IPRINT .GT. 5) THEN
         CALL AROUND(' START OF CC_TRDRV ')
         IF (DIRECT) WRITE(LUPRI,*) ' Atomic direct calculation'
         IF (CCR12)  WRITE(LUPRI,*) ' CC-R12 calculation'
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
      IF (CCR12) NCAMP = NCAMP + NGAMMA(ISYMTR)
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
 
C        -----------------------------------------------
C        Prepare memory depending on which case we have:
C        -----------------------------------------------

         IF (.NOT. TRIPLET) THEN
C           ------------
C           Allocations:
C           ------------
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
C           ---------------------------------
C           Read the CC trial vectors from disk.
C           ---------------------------------
            IF (IPRINT .GT. 45) CALL AROUND('CC_TRDRV: C1')
            DO 150 IV = 1, NSIMTR
               KOFF1  = KC1AM + NT1AM(ISYMTR)*(IV - 1)
               CALL CC_RVEC(LUFC1,FC1AM,NT1AM(ISYMTR),NT1AM(ISYMTR),
     *                      IV+K1-1,WORK(KOFF1))
               IF (IPRINT.GT.45 .AND. (CCS.OR.(.NOT.MINSCR))) THEN
                  WRITE(LUPRI,*) 'Vector nr.',IV+K1-1
                  CALL CC_PRP(WORK(KOFF1),WORK(KRHO2),ISYMTR,1,0)
               ENDIF
  150       CONTINUE
C
            IF ((.NOT.CCS).AND.( MINSCR)) THEN
               CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                      K1,WORK(KRHO2))
               IF (IPRINT .GT. 45) THEN
                  CALL AROUND('CC_TRDRV: (C1,C2) packed ')
                  WRITE(LUPRI,*) 'Vector nr.',K1
                  CALL CC_PRP(WORK(KC1AM),WORK(KRHO2),ISYMTR,1,1)
               ENDIF
            ENDIF
C
C           -----------------------------------------------
C           Test with finited difference CC jacobian if FDJAC.
C           -----------------------------------------------
            IF (FDJAC .OR. JACEXP) THEN
               IF (CCR12) THEN
                 KC12AM = KEND1
                 KS12AM = KC12AM + NGAMMA(ISYMTR)
                 KRHO12 = KS12AM + NGAMMA(ISYMTR)
                 KEND1  = KRHO12 + NGAMMA(ISYMTR)
                 LWRK1  = LWORK  - KEND1             
                 IF (LWRK1 .LE. 0 )
     *              CALL QUIT('Too little workspace in cclr_trr (2)')

                 CALL CC_RVEC(LUFC12,FC12AM,NGAMMA(ISYMTR),
     *                        NGAMMA(ISYMTR),K1,WORK(KC12AM))
                 IF (IPRINT .GT. 45) THEN
                   CALL AROUND('R12 double excitation part of vector')
                   CALL OUTPUT(WORK(KRHO12),1,NGAMMA(ISYMTR),1,1,
     *                         NGAMMA(ISYMTR),1,1,LUPRI)
                 END IF
               END IF

               CALL DCOPY(NT2AM(ISYMTR),WORK(KRHO2),1,WORK(KC2AM),1)
               CALL CCLR_DUMTRR(WORK(KC1AM),WORK(KRHO1),WORK(KRHO2),
     *                          WORK(KRHO12),WORK(KS12AM),
     *                          WORK(KEND1),LWRK1)
C
C              ------------------------------
C              Print the transformed vectors.
C              ------------------------------
               IF (IPRINT .GT. 50) THEN
                  CALL AROUND('CC_TRDRV: RHO  trans. by DUMTRR .')
                  CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYMTR,1,1)
                  IF (CCR12) THEN
                    WRITE(LUPRI,*) 'A * C_12:'
                    CALL OUTPUT(WORK(KRHO12),1,NGAMMA(ISYMTR),1,1,
     *                          NGAMMA(ISYMTR),1,1,LUPRI)
                    WRITE(LUPRI,*) 'S * C_12:'
                    CALL OUTPUT(WORK(KS12AM),1,NGAMMA(ISYMTR),1,1,
     *                          NGAMMA(ISYMTR),1,1,LUPRI)
                  END IF
               ENDIF

               CALL CC_WVEC(LUFR1,FRHO1,NT1AM(ISYMTR),
     *                      NT1AM(ISYMTR),K1,WORK(KRHO1))
               CALL CC_WVEC(LUFR2,FRHO2,NT2AM(ISYMTR),
     *                      NT2AM(ISYMTR),K1,WORK(KRHO2))
               IF (CCR12) THEN
                 CALL CC_WVEC(LUFR12,FRHO12,NGAMMA(ISYMTR),
     *                        NGAMMA(ISYMTR),K1,WORK(KRHO12))
                 CALL CC_WVEC(LUFS12,FS12AM,NGAMMA(ISYMTR),
     *                        NGAMMA(ISYMTR),K1,WORK(KS12AM))
               END IF
            ELSE
C
C              ----------------------------------------
C              Prepare the C-amplitudes.
C              ----------------------------------------
               IF ( .NOT. CCS ) THEN
                  IF ( ISIDE .GE. 1) THEN
                     CALL CCLR_DIASCL(WORK(KRHO2),TWO,ISYMTR)
                  ENDIF
                  CALL CC_T2SQ(WORK(KRHO2),WORK(KC2AM),ISYMTR)
                  IF (IPRINT.GT.50) THEN
                     CALL AROUND('CC_TRDRV: (C1,C2) squared ')
                     CALL CC_PRSQ(WORK(KC1AM),WORK(KC2AM),ISYMTR,1,2)
                  ENDIF
               ENDIF
C
C              ----------------
C              Zero rho vector.
C              ----------------
               CALL DZERO(WORK(KRHO1),NT1AM(ISYMTR)*NSIMTR)
               CALL DZERO(WORK(KRHO2),NRHO2)
C
C              -------------------
C              File read and save.
C               ------------------
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
C           -----------------------------
C           For Triplet start at the
C           beginning of the workspace.
C           -----------------------------
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
C           -------------------------------------
C           File opening for the triplet case
C           -------------------------------------
            IF (.NOT. (CCS.OR.CC2)) THEN 
               CALL WOPEN2(LUFSD,FR2SD,64,0)
            ENDIF
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
               IF (CCR12) CALL QUIT('No triplet yet for CCR12')
               CALL CC_RHTR3(ECURR,FRHO1,LUFR1,FR2SD,LUFSD,
     *                       FC1AM,LUFC1,FC2AM,LUFC2,
     *                       WORK,LWORK,NSIMTR,
     *                       IVEC,ITR)
C
            ELSE
               IF (ISIDE .EQ. 1) THEN
                  CALL CC_RHTR(ECURR,
     *                         FRHO1,LUFR1,FR2SD,LUFSD,FRHO12,LUFR12,
     *                         FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                         WORK(KRHO1),WORK(KRHO2),
     *                         WORK(KC1AM),WORK(KC2AM),
     *                         WORK(KEND1),LWRK1,NSIMTR,
     *                         IVEC,ITR,LRHO1,.FALSE.,DUMMY)
               ELSE IF (ISIDE .EQ. -1) THEN
                  IF (CCR12) CALL QUIT('No left transform. for CCR12')
                  CALL CC_LHTR(ECURR,FRHO1,LUFR1,FR2SD,LUFSD,
     *                         FC1AM,LUFC1,FC2AM,LUFC2,
     *                         WORK(KRHO1),WORK(KRHO2),
     *                         WORK(KC1AM),WORK(KC2AM),
     *                         WORK(KEND1),LWRK1,NSIMTR,
     *                         IVEC,ITR,LRHO1)
               ELSE
                  CALL QUIT(' ISIDE should be -1 or +1 ')
               ENDIF

               IF ( CCR12 ) THEN
                 DO IV = 1, NSIMTR
                  CALL CC_R12METRIC(ISYMTR,WORK(KEND1),LWRK1,
     *                       FC12AM,LUFC12,FS12AM,LUFS12,IVEC-1+IV)
                 END DO
               END IF

C
            ENDIF
C
C           ------------------------------
C           Print the transformed vectors.
C           ------------------------------
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

                 IF (.NOT.(CCS.OR.CC2) .AND. (.NOT. TRIPLET)) THEN
                    CALL CC_WVEC(LUFR2,FRHO2,NT2AM(ISYMTR),
     *                           NT2AM(ISYMTR),NR1,WORK(KRHO2))
                 ENDIF

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
