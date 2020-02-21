C
C  /* Deck so_trial3 */
      SUBROUTINE SO_TRIAL3(MODEL,NNEWTR,ISYMTR,IMAGPROP,
     &                     IFREQ,FRVAL,NFRVAL,NEXCI,
     &                     DENSIJ,LDENSIJ,DENSAB,LDENSAB,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C     Andrea Ligabue December 2003
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Determine the initial trialvectors for the linear
C     response equations using eq. (19) with R = (C +-C).
C
C   INPUT:
C     ISYMTR    Symmetry of trial vectors to generate
C     IMAGPROP  Is property imaginary?
C     NEXCI     Number of trial vectors to create( 1?)
C     IFREQ     Freq counter in the FRVAL(NFRVAL) array
C     FRVAL(NFRVAL) All frequencies to consider
C               Density Matrices:
C     DENSIJ(LDENSIJ), DENSAB(LDENSAB)
C
C   OUTPUT:
C     NNEWTR    Number of trial vectors created
C               (now always == NEXCI)
C
C   SCRATCH:
C     WORK(LWORK)
CW
      use so_info, only: so_has_doubles, sop_stat_trh
C
#include "implicit.h"
#include "priunit.h"
C
#include "soppinf.h"
#include "ccsdsym.h"
C
      PARAMETER (ONE = 1.0D0, STHR = 1.0D-5)
      PARAMETER (HALF = 0.5D0 ,TWO = 2.0D0)
C
      DIMENSION   DENSIJ(LDENSIJ), DENSAB(LDENSAB)
      DIMENSION   FRVAL(NFRVAL)
      DIMENSION   WORK(LWORK)
C
      LOGICAL     IMAGPROP
      LOGICAL     STATIC, DOUBLES
C
      CHARACTER*5 MODEL
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_TRIAL3')
C
      STATIC = ABS(FRVAL(IFREQ)).LT.SOP_STAT_TRH
C
CRF Doubles branches should be usefull to enable HRPA.
CRF However, RPA still has it's own version of the routine,
CRF since in that case we have no explicit
      DOUBLES = SO_HAS_DOUBLES(MODEL)

C
C---------------------------------
C     1. allocation of work space.
C---------------------------------
C
      LEDIA1 = NT1AM(ISYMTR)
      IF (DOUBLES) THEN
         LEDIA2 = N2P2HOP(ISYMTR)
      ELSE
         LEDIA2 = 0
      ENDIF
      LSDIA1 = NT1AM(ISYMTR)
C
      KEDIA1  = 1
      KEDIA2  = KEDIA1 + LEDIA1
      KSDIA1  = KEDIA2 + LEDIA2
      KEND1   = KSDIA1 + LSDIA1
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('SO_TRIAL3.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_TRIAL3.1',' ',KEND1,LWORK)
C
C------------------------------------------
C     Read diagonal E[2] and S[2] elements.
C------------------------------------------
C
      CALL GPOPEN  (LUDIAG,'SO_DIAG','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
      REWIND LUDIAG
C
      READ(LUDIAG) ( WORK(KEDIA1+I-1), I = 1,LEDIA1)
      IF (DOUBLES) THEN
         READ(LUDIAG) ( WORK(KEDIA2+I-1), I = 1,LEDIA2)
      ENDIF
      READ(LUDIAG) ( WORK(KSDIA1+I-1), I = 1,LSDIA1)
C
      CALL GPCLOSE (LUDIAG,'KEEP')
C
C---------------------------------
C     2. allocation of work space.
C---------------------------------
C
cLig  I try to use less memeory as I can so I read the part of the GP I
cLig  need, compute the part of the trial, write on file and then again
cLig  ... so I can use the same space, but I need at least nt2am for
cLig  TV and nt2am * 2 for gp vector. The negative effect of that is
cLig  that I can really write only the TV in this format (1E 1D 2E 2D)
      LTRIAL  = LEDIA2
CRF Generate D part from E part
      LGPVC1  = LEDIA1 !* 2
      IF (DOUBLES) THEN
         LGPVC2  = LEDIA2 !* 2
         LTRIAL  = LEDIA2
      ELSE
         LGPVC2  = 0
         LTRIAL  = LEDIA1
      ENDIF

C
      KTRIAL  = KEND1
      KGPVEC  = KTRIAL + LTRIAL
      KEND2   = KGPVEC + LGPVC2
      LWORK2  = LWORK  - KEND2
C
      CALL SO_MEMMAX ('SO_TRIAL3.2',LWORK2)
      IF (LWORK2 .LT. 0) CALL STOPIT('SO_TRIAL3.2',' ',KEND2,LWORK)
CRF Factor for creating D gradient from E gradient
      DFACTOR  = -ONE
      IF (IMAGPROP) DFACTOR = ONE
C
C--------------------------------------------
C     Open files for storing of trialvectors.
C--------------------------------------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LEDIA1)
      CALL SO_OPEN(LUTR1D,FNTR1D,LEDIA1)
      IF (DOUBLES) THEN
         CALL SO_OPEN(LUTR2E,FNTR2E,LEDIA2)
         CALL SO_OPEN(LUTR2D,FNTR2D,LEDIA2)
      ENDIF
C
C---------------------------------------------------------
C     Create initial new trial vectors and write to files.
C---------------------------------------------------------
C
      NNEWTR = NEXCI
cLig  just another stupid loop, since nexci is always 1 and nnewtr =
cLig  nexci ... remove ?
C
C      DO 100 INEWTR=1,NNEWTR
      INEWTR = 1
C
C----------------------------------------------------------------------
C     For the first frequency: use the GP to create the starrting trial
C                              vector.
C     Otherwise use the last trial vector of the previous frequency
C----------------------------------------------------------------------
C
      IF((IFREQ.EQ.1).OR.STATIC) THEN
C
C------------------------------------------------------------------
C           Read the 1E and 1D part of the  GP vector from the file
C------------------------------------------------------------------
C
            CALL SO_OPEN(LUGPV1,FNGPV1,LGPVC1)
C
            CALL SO_READ(WORK(KGPVEC),LGPVC1,LUGPV1,FNGPV1,1)
C
            CALL SO_CLOSE(LUGPV1,FNGPV1,'KEEP')
C
C-------------------------------------------------------------
C           Calculate the 1E part of the trial vector as
C           (E[2]diag - omega*S[2]diag) * GP and write on file
C-------------------------------------------------------------
C
            DO IELEM=1,LEDIA1
C
              TMP = ONE / ( WORK(KEDIA1+IELEM-1) - FRVAL(IFREQ) *
     &                      WORK(KSDIA1+IELEM-1) )
              WORK(KTRIAL+IELEM-1) = TMP * WORK(KGPVEC+IELEM-1)
C
            END DO
C
            CALL SO_WRITE(WORK(KTRIAL),LEDIA1,LUTR1E,FNTR1E,INEWTR)
C
C-------------------------------------------------------------
C           Calculate the 1D part of the trial vector as
C           (E[2]diag + omega*S[2]diag) * GP and write on file
C           Only if not a static calculation.
C-------------------------------------------------------------
C
            IF(.NOT.STATIC) THEN
C
              DO IELEM=1,LEDIA1
C
                TMP = DFACTOR / ( WORK(KEDIA1+IELEM-1) + FRVAL(IFREQ) *
     &                        WORK(KSDIA1+IELEM-1) )
                WORK(KTRIAL+IELEM-1) = TMP * WORK(KGPVEC+IELEM-1)
C
              END DO
C
            ELSE
               CALL DZERO(WORK(KTRIAL),LTRIAL)
            END IF
C
            CALL SO_WRITE(WORK(KTRIAL),LEDIA1,LUTR1D,FNTR1D,INEWTR)
C
            IF (DOUBLES) THEN
C
C------------------------------------------------------------------
C           Read the 2E and 2D part of the  GP vector from the file
C------------------------------------------------------------------
C
               CALL SO_OPEN(LUGPV2,FNGPV2,LGPVC2)
C
               CALL SO_READ(WORK(KGPVEC),LGPVC2,LUGPV2,FNGPV2,1)
C
               CALL SO_CLOSE(LUGPV2,FNGPV2,'KEEP')
C
C------------------------------------------------------------------
C           Transform the 2E and 2D part of the  GP vector from the
C           right to left side
C------------------------------------------------------------------
C
               IF(IPRSOP.GT.100) THEN
C
                 CALL AROUND('GP before LTR transformation')
                 CALL OUTPUT(WORK(KGPVEC),1,LEDIA2,1,2,LEDIA2,2,1,LUPRI)
C
               ENDIF
C
               IF (.NOT.TRIPLET)
     &              CALL SO_TMLTR(WORK(KGPVEC),HALF,ISYMTR)
C
               IF(IPRSOP.GT.10) THEN
C
                 CALL AROUND('GP after LTR transformation')
                 CALL OUTPUT(WORK(KGPVEC),1,LEDIA2,1,1,LEDIA2,1,1,LUPRI)
C
               ENDIF
C
C-------------------------------------------------------------
C           Calculate the 2E part of the trial vector as
C           (E[2]diag - omega) * GP and write on file
C-------------------------------------------------------------
C
               DO IELEM=1,LEDIA2
C
                  TMP = ONE / ( WORK(KEDIA2+IELEM-1) - FRVAL(IFREQ) )
                  WORK(KTRIAL+IELEM-1) = TMP * WORK(KGPVEC+IELEM-1)
C
               END DO
C
               CALL SO_WRITE(WORK(KTRIAL),LEDIA2,LUTR2E,FNTR2E,INEWTR)
C
C-------------------------------------------------------------
C           Calculate the 2D part of the trial vector as
C           (E[2]diag + omega) * GP and write on file
C-------------------------------------------------------------
C
               IF(.NOT.STATIC) THEN
C
                  DO IELEM=1,LEDIA2
C
                     TMP = DFACTOR /
     &                     ( WORK(KEDIA2+IELEM-1) + FRVAL(IFREQ) )
                     WORK(KTRIAL+IELEM-1) = TMP *
     &                       WORK(KGPVEC+IELEM-1)
C
                  END DO
C
               ELSE
C
                  CALL DZERO(WORK(KTRIAL),LTRIAL)
C
               END IF
C
               CALL SO_WRITE(WORK(KTRIAL),LEDIA2,LUTR2D,FNTR2D,INEWTR)
            ENDIF
C
C---------------------------------------------------------------
C        For other frequencies use the last trial vector on file
C---------------------------------------------------------------
C
         ELSE
C
cspas should NOLDTR not be 1 here ?
cLig  this operation is quite stupid since it goes to read and write
cLig  always in the same space ... since inewtr is also always 1 !
C
CRF --- No need to do anything!!!
            NOLDTR = 1

         ENDIF
C
C  100 CONTINUE
C
      IF( IPRSOP .GE. 15) THEN
C
C--------------------------------------
C     Write new trial vecotrs to output
C--------------------------------------
C
         DO 200 INEWTR=1,NNEWTR

            WRITE(LUPRI,'(I3,A,A)') INEWTR,
     &         '. raw trial vector in SO_TRIAL3 before',
     &         ' orthonormalization'
C
            CALL SO_READ(WORK(KTRIAL),LEDIA1,LUTR1E,FNTR1E,INEWTR)
            WRITE (LUPRI,'(I8,1X,F14.8)')
     &           (I,WORK(KTRIAL+I-1),I=1,LEDIA1)
            CALL SO_READ(WORK(KTRIAL),LEDIA1,LUTR1D,FNTR1D,INEWTR)
            WRITE (LUPRI,'(I8,1X,F14.8)')
     &           (I,WORK(KTRIAL+I-1),I=1,LEDIA1)
            IF (DOUBLES) THEN
               CALL SO_READ(WORK(KTRIAL),LEDIA2,LUTR2E,FNTR2E,INEWTR)
               WRITE (LUPRI,'(I8,1X,F14.8)')
     &              (I,WORK(KTRIAL+I-1),I=1,LEDIA2)
               CALL SO_READ(WORK(KTRIAL),LEDIA2,LUTR2D,FNTR2D,INEWTR)
               WRITE (LUPRI,'(I8,1X,F14.8)')
     &              (I,WORK(KTRIAL+I-1),I=1,LEDIA2)
            END IF

C
  200    CONTINUE
C
      ENDIF
C
C-----------------------------------------------
C     Orthogonalize new trial vectors over S[2].
C-----------------------------------------------
C
Cekd: Same problem as for RP_ORTH_TRN, Cannot send in explicit zeros
C     as they may be assigned inside routine.
Corig  CALL SO_ORTH_TRN(0,NNEWTR,0,ISYMTR,
C     &                 DENSIJ,LDENSIJ,DENSAB,LDENSAB,WORK,LWORK)
C Needed? For normalizing new vector?
      NLINDP = 0
      NOLDTR = 0
      CALL SO_ORTH_TRN(DOUBLES,'LINEAR',NOLDTR,NNEWTR,NLINDP,ISYMTR,
     &                 DENSIJ,LDENSIJ,DENSAB,LDENSAB,WORK(KEND2),LWORK2)
C
      IF( IPRSOP .GE. 10) THEN
C
C------------------------------------------------------
C     Write new trial orthogonalized  vecotrs to output
C------------------------------------------------------
C
         DO 300 INEWTR=1,NNEWTR

            WRITE(LUPRI,'(I3,A,A)') INEWTR,
     &          '. raw trial vector in RP_TRIAL3 after',
     &          ' orthonormalization'
C
            CALL SO_READ(WORK(KTRIAL),LEDIA1,LUTR1E,FNTR1E,INEWTR)
            WRITE (LUPRI,'(I8,1X,F14.8)')
     &            (I,WORK(KTRIAL+I-1),I=1,LEDIA1)
            CALL SO_READ(WORK(KTRIAL),LEDIA1,LUTR1D,FNTR1D,INEWTR)
            WRITE (LUPRI,'(I8,1X,F14.8)')
     &            (I,WORK(KTRIAL+I-1),I=1,LEDIA1)
            IF (DOUBLES) THEN
               CALL SO_READ(WORK(KTRIAL),LEDIA2,LUTR2E,FNTR2E,INEWTR)
               WRITE (LUPRI,'(I8,1X,F14.8)')
     &             (I,WORK(KTRIAL+I-1),I=1,LEDIA2)
               CALL SO_READ(WORK(KTRIAL),LEDIA2,LUTR2D,FNTR2D,INEWTR)
               WRITE (LUPRI,'(I8,1X,F14.8)')
     &             (I,WORK(KTRIAL+I-1),I=1,LEDIA2)
            END IF
C
  300    CONTINUE
C
      ENDIF
C
C------------------------------------
C     Close files with trial vectors.
C------------------------------------
C
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
      IF (DOUBLES) THEN
         CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
         CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
      ENDIF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('SO_TRIAL3')
C
      RETURN
      END
