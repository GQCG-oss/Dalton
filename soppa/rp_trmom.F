C
C  /* Deck rp_trmom */
CSPAS:23/5-11: second and third moment sum rules
C     SUBROUTINE RP_TRMOM(ISYMTR,NLBTOT,NEXCI,TRLEN,TRVEL,TRLON,TRMAG,
C    &                    BSRLON,EXENG,WORK,LWORK)
      SUBROUTINE RP_TRMOM(ISYMTR,NLBTOT,NEXCI,TRLEN,TRVEL,TQLEN,TQVEL,
     &                    TTLEN,TRLON,TRMAG,BSRLON,EXENG,WORK,LWORK)
CKeinSPASmehr
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C     Andrea Ligabue, December 2003: linear response functions 
C                                    implemented
C
C     PURPOSE: Calculate RPA transition moments.
C
#include <implicit.h>
#include <priunit.h>
C
#include <cbiexc.h>
CSPAS:15/03-2006: merge with Dalton-2.0  
C#include <cbilrs.h>
CKeinSPASmehr
#include <ccsdsym.h>
#include <ccorb.h>
#include <soppinf.h>
#include <mxcent.h>
#include <nuclei.h>
C
      DIMENSION TRLEN(3,NSYM,MXNEXI), TRVEL(3,NSYM,MXNEXI)
CSPAS:23/5-11: second and third moment sum rules
      DIMENSION TQLEN(3,3,NSYM,MXNEXI), TQVEL(3,3,NSYM,MXNEXI)
      DIMENSION TTLEN(10,NSYM,MXNEXI)
CKeinSPASmehr
      DIMENSION TRLON(3,NSYM,MXNEXI), TRMAG(3,NSYM,MXNEXI)
      DIMENSION BSRLON(3,NSYM,MXNEXI),EXENG(NSYM,MXNEXI)
      DIMENSION WORK(LWORK)
C
      CHARACTER*8 LABEL
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_TRMOM')
C
C-------------------------------
C     Loop over property labels.
C-------------------------------
C
      DO 200 IPRLBL = 1, NLBTOT
C
         LABEL  = LABAPP(IPRLBL)
         KSYMOP = LABSYM(IPRLBL)
C
C-----------------------------------------------------------------
C        Check if integrals are second moment of charge integrals.
C-----------------------------------------------------------------
C
CSPAS:23/5-11: second and third moment sum rules
C        IF ( LABEL(3:8) .EQ. 'SECMOM' ) GOTO 200
CKeinSPASmehr
C
C----------------------------------------------------------------------
C        If KSYMOP (the operator symmetry) equals ISYMTR (the reference
C        state symmetry times excited state symmmety) calculate the
C        transition moment.
C----------------------------------------------------------------------
C
         IF (KSYMOP .EQ. ISYMTR) THEN
C
C-------------------------------------
C           Allocation of work space.
C-------------------------------------
C
            LGPVEC = 2 * NT1AM(ISYMTR)
            LGPVCH = NT1AM(ISYMTR)
            LTR1E  = NT1AM(ISYMTR)
            LTR1D  = NT1AM(ISYMTR)
            LTRMOM = NEXCI
            LTRME  = NEXCI
            LTRMD  = NEXCI
C
            KGPVEC  = 1
            KTR1E   = KGPVEC + LGPVEC
            KTR1D   = KTR1E  + LTR1E
            KTRMOM  = KTR1D  + LTR1D
            KTRME   = KTRMOM + LTRMOM
            KTRMD   = KTRME  + LTRME 
            KEND1   = KTRMD  + LTRMD
            LWORK1  = LWORK  - KEND1
C
            CALL SO_MEMMAX ('RP_TRMOM',LWORK1)
            IF (LWORK1 .LT. 0) CALL STOPIT('RP_TRMOM',' ',KEND1,LWORK)
C
C----------------------------------------------------------------------
C           Calculate the gradient property vector and write to output.
C----------------------------------------------------------------------
C
            CALL RP_GETGP(WORK(KGPVEC),LGPVEC,LABEL,ISYMTR,WORK(KEND1),
     &                    LWORK1)
C
            IF ( IPRSOP .GE. 5) THEN
               WRITE(LUPRI,'(2A)') 'In RP_TRMOM : GP vector: ',LABEL
               CALL OUTPUT(WORK(KGPVEC),1,LGPVCH,1,2,LGPVCH,2,1,LUPRI)
            ENDIF
C
C--------------------------------------------
C           Open files with solution vectors.
C--------------------------------------------
C
            CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
            CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
C
C-----------------------------------------------------------
C           Loop over excitations of the specified symmetry.
C-----------------------------------------------------------
C
            DO 100 IEXCI = 1,NEXCI
C
C------------------------------------
C              Read solution vectors.
C------------------------------------
C
               CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,IEXCI)
               CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,IEXCI)
C
               IF ( IPRSOP .GE. 5) THEN
                  WRITE(LUPRI,'(A)') 'Solution vector: '
                  WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)')
     &                 (I,WORK(KTR1E+I-1),WORK(KTR1D+I-1),I=1,LTR1E)
               ENDIF
C
C-------------------------------------------------------------
C              Calculate excitation part of transition moment.
C-------------------------------------------------------------
C
               WORK(KTRME+IEXCI-1) =
     &            DDOT(LTR1E,WORK(KGPVEC),1,WORK(KTR1E),1)
C
C----------------------------------------------------------------
C              Calculate de-excitation part of transition moment.
C----------------------------------------------------------------
C
c               IF ( (LABEL(2:7) .EQ. 'DIPLEN') ) THEN
C
c                  WORK(KTRMD+IEXCI-1) = 
c     &               -DDOT(LGPVEC,WORK(KGPVEC),1,WORK(KTR1D),1)
C
c               ELSE IF ( (LABEL(2:7) .EQ. 'DIPVEL') .OR.
c     &                   (LABEL(2:7) .EQ. 'ANGMOM')      ) THEN
C
c                  WORK(KTRMD+IEXCI-1) = 
c     &               DDOT(LGPVEC,WORK(KGPVEC),1,WORK(KTR1D),1)
C
c               ELSE
C
c                  WRITE(LUPRI,*) 'Revision of program needed for',LABEL
c                  GO TO 200
C
c               END IF
               WORK(KTRMD+IEXCI-1) =
     &            DDOT(LTR1D,WORK(KGPVEC+LTR1E),1,WORK(KTR1D),1)
C
C--------------------------------------------------------------------
C              Calculate total transition moment and write to output.
C--------------------------------------------------------------------
C
               WORK(KTRMOM+IEXCI-1)  = WORK(KTRME+IEXCI-1) 
     &                               + WORK(KTRMD+IEXCI-1)
C
               IF ( IPRSOP .GE. 4) THEN
C
                  WRITE(LUPRI,'(/,1X,A,A,I3,A)')
     &            LABEL,' transition moment for state',IEXCI,':'
                  WRITE(LUPRI,'(2(F10.5,A),F10.5)')
     &            WORK(KTRME+IEXCI-1),' + ',WORK(KTRMD+IEXCI-1),
     &            ' = ',WORK(KTRMOM+IEXCI-1)
C
               ENDIF
C
C---------------------------------------------------
C              Collect transition moments in arrays.
C---------------------------------------------------
C
               CALL SO_COLLECT_TM(ISYMTR,LABEL,IEXCI,
     &                            WORK(KTRMOM+IEXCI-1),TRLEN,TRVEL,
CSPAS:23/5-11: second and third moment sum rules
     &                            TQLEN,TQVEL,TTLEN,
CKeinSPASmehr
     &                            TRLON,TRMAG,BSRLON,EXENG)
C
  100       CONTINUE
C
C--------------------------------------------
C           Close files with solution vectors.
C--------------------------------------------
C
            CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
            CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
C
C-------------------------------------------------------
C           Write tabel of transition moments to output.
C-------------------------------------------------------
C
            IF ( IPRSOP .GE. 1) THEN
C
               WRITE(LUPRI,'(/,1X,A,A,I3)')
     &         LABEL,' RPA transition moments for excitation symmetry',
     &         ISYMTR
               WRITE(LUPRI,9001)
               WRITE(LUPRI,'(A)')
     &         ' Exci.    1p-1h        1h-1p        Total'
               WRITE(LUPRI,9002)
C
               DO 150 IEXCI = 1,NEXCI
C
                  WRITE(LUPRI,'(1X,I3,2X,2(F10.5,A),F10.5)') IEXCI,
     &            WORK(KTRME+IEXCI-1),' + ',WORK(KTRMD+IEXCI-1),
     &            ' = ',WORK(KTRMOM+IEXCI-1)
C
  150          CONTINUE
C
               WRITE(LUPRI,9001)
C
            ENDIF
C
         END IF
C
  200 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_TRMOM')
C
      RETURN
C
 9001 FORMAT(1X,'=========================================')
 9002 FORMAT(1X,'-----------------------------------------')
C
      END
