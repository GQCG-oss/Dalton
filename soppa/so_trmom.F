#include <single.h>
C
C  /* Deck so_trmom */
      SUBROUTINE SO_TRMOM(MODEL,ISYMTR,NLBTOT,NEXCI,T2AM,LT2AM,DENSIJ,
     &                    LDENSIJ,DENSAB,LDENSAB,DENSAI,LDENSAI,
CSPAS:23/5-11: second and third moment sum rules
C    &                    TRLEN,TRVEL,TRLON,TRMAG,BSRLON,EXENG,
C    &                    WORK,LWORK)
     &                    TRLEN,TRVEL,TQLEN,TQVEL,TTLEN,TRLON,TRMAG,
     &                    BSRLON,EXENG,WORK,LWORK)
CKeinSPASmehr
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate RPA(D), SOPPA and SOPPA(CCSD) transition 
C              moments.
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
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), DENSAI(LDENSAI)
      DIMENSION TRLEN(3,NSYM,MXNEXI), TRVEL(3,NSYM,MXNEXI)
CSPAS:23/5-11: second and third moment sum rules
      DIMENSION TQLEN(3,3,NSYM,MXNEXI), TQVEL(3,3,NSYM,MXNEXI)
      DIMENSION TTLEN(10,NSYM,MXNEXI)
CKeinSPASmehr
      DIMENSION TRLON(3,NSYM,MXNEXI), TRMAG(3,NSYM,MXNEXI)
      DIMENSION BSRLON(3,NSYM,MXNEXI), EXENG(NSYM,MXNEXI)
      DIMENSION T2AM(LT2AM), WORK(LWORK)
C
      CHARACTER*8 LABEL
      CHARACTER*5 MODEL
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_TRMOM')
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
C---------------------------------------------------------------------
C       If KSYMOP (the operator symmetry) equals ISYMTR (the reference
C       state symmetry times excited state symmmety) calculate the
C       transition moment.
C---------------------------------------------------------------------
C
        IF (KSYMOP .EQ. ISYMTR) THEN
C
C-------------------------------------
C           Allocation of work space.
C-------------------------------------
C
            LGPVC1  = 2 * NT1AM(ISYMTR) 
            LGPVC2  = 2 * NT2AM(ISYMTR) 
            LGPVC1H = NT1AM(ISYMTR) 
            LGPVC2H = NT2AM(ISYMTR) 
            LTR1    = NT1AM(ISYMTR)
            LTR2    = NT2AM(ISYMTR)
            LTRMOM  = NEXCI
            LTRM1E  = NEXCI
            LTRM1D  = NEXCI
            LTRM2E  = NEXCI
            LTRM2D  = NEXCI
C
            KGPVC1  = 1
            KGPVC2  = KGPVC1 + LGPVC1
            KTR1    = KGPVC2 + LGPVC2
            KTR2    = KTR1   + LTR1 
            KTRMOM  = KTR2   + LTR2 
            KTRM1E  = KTRMOM + LTRMOM
            KTRM1D  = KTRM1E + LTRM1E
            KTRM2E  = KTRM1D + LTRM1D
            KTRM2D  = KTRM2E + LTRM2E
            KEND1   = KTRM2D + LTRM2D 
            LWORK1  = LWORK  - KEND1
C
            CALL SO_MEMMAX ('SO_TRMOM',LWORK1)
            IF (LWORK1 .LT. 0) CALL STOPIT('SO_TRMOM',' ',KEND1,LWORK)
C
C----------------------------------------------------------------------
C           Calculate the gradient property vector and write to output.
C----------------------------------------------------------------------
C
            CALL SO_GETGP(WORK(KGPVC1),LGPVC1,WORK(KGPVC2),LGPVC2,
     &                    LABEL,ISYMTR,T2AM,LT2AM,DENSIJ,LDENSIJ,
     &                    DENSAB,LDENSAB,DENSAI,LDENSAI,
     &                    WORK(KEND1),LWORK1)
C
            IF ( IPRSOP .GE. 5) THEN
               WRITE(LUPRI,'(2A)') 'In SO_TRMOM : GP vector: ',LABEL
               CALL OUTPUT(WORK(KGPVC1),1,LGPVC1H,1,2,LGPVC1H,2,1,LUPRI)
               CALL OUTPUT(WORK(KGPVC2),1,LGPVC2H,1,2,LGPVC2H,2,1,LUPRI)
            ENDIF
C
C--------------------------------------------
C           Open files with solution vectors.
C--------------------------------------------
C
            CALL SO_OPEN(LUTR1E,FNTR1E,LTR1)
            CALL SO_OPEN(LUTR1D,FNTR1D,LTR1)
            CALL SO_OPEN(LUTR2E,FNTR2E,LTR2)
            CALL SO_OPEN(LUTR2D,FNTR2D,LTR2)
C
C-----------------------------------------------------------
C           Loop over excitations of the specified symmetry.
C-----------------------------------------------------------
C
            DO 100 IEXCI = 1,NEXCI
C
C-------------------------------------------------------
C              Read excitation part of solution vectors.
C-------------------------------------------------------
C
               CALL SO_READ(WORK(KTR1),LTR1,LUTR1E,FNTR1E,IEXCI)
               CALL SO_READ(WORK(KTR2),LTR2,LUTR2E,FNTR2E,IEXCI)
C
               IF ( IPRSOP .GE. 5) THEN
                  WRITE(LUPRI,'(A)') 'Exci. part of solution vector: '
                  CALL OUTPUT(WORK(KTR1),1,LTR1,1,1,LTR1,1,1,LUPRI)
                  CALL OUTPUT(WORK(KTR2),1,LTR2,1,1,LTR2,1,1,LUPRI)
               ENDIF
C
C-------------------------------------------------------------
C              Calculate excitation part of transition moment.
C-------------------------------------------------------------
C
               WORK(KTRM1E+IEXCI-1) = 
     &            DDOT(LGPVC1H,WORK(KGPVC1),1,WORK(KTR1),1)
               WORK(KTRM2E+IEXCI-1) = 
     &            DDOT(LGPVC2H,WORK(KGPVC2),1,WORK(KTR2),1)
C
C----------------------------------------------------------
C              Read de-excitation part of solution vectors.
C----------------------------------------------------------
C
               CALL SO_READ(WORK(KTR1),LTR1,LUTR1D,FNTR1D,IEXCI)
               CALL SO_READ(WORK(KTR2),LTR2,LUTR2D,FNTR2D,IEXCI)
C
               IF ( IPRSOP .GE. 5) THEN
                  WRITE(LUPRI,'(A)') 'De-ex. part of solution vector: '
                  CALL OUTPUT(WORK(KTR1),1,LTR1,1,1,LTR1,1,1,LUPRI)
                  CALL OUTPUT(WORK(KTR2),1,LTR2,1,1,LTR2,1,1,LUPRI)
               ENDIF
C
C----------------------------------------------------------------
C              Calculate de-excitation part of transition moment.
C----------------------------------------------------------------
C
c               IF ( (LABEL(2:7) .EQ. 'DIPLEN') ) THEN
C
c                  WORK(KTRM1D+IEXCI-1) = 
c     &               -DDOT(LGPVC1,WORK(KGPVC1),1,WORK(KTR1),1)
c                  WORK(KTRM2D+IEXCI-1) = 
c     &               -DDOT(LGPVC2,WORK(KGPVC2),1,WORK(KTR2),1)
C
c               ELSE IF ( (LABEL(2:7) .EQ. 'DIPVEL') .OR.
c     &                   (LABEL(2:7) .EQ. 'ANGMOM')      ) THEN
C
c                  WORK(KTRM1D+IEXCI-1) = 
c     &               DDOT(LGPVC1,WORK(KGPVC1),1,WORK(KTR1),1)
c                  WORK(KTRM2D+IEXCI-1) = 
c     &               DDOT(LGPVC2,WORK(KGPVC2),1,WORK(KTR2),1)
C
c               ELSE
C
c                  WRITE(LUPRI,*) 'Revision of program needed for',LABEL
c                  GO TO 200
C
c               END IF
               WORK(KTRM1D+IEXCI-1) = 
     &             DDOT(LGPVC1H,WORK(KGPVC1+LGPVC1H),1,WORK(KTR1),1)
               WORK(KTRM2D+IEXCI-1) = 
     &             DDOT(LGPVC2H,WORK(KGPVC2+LGPVC2H),1,WORK(KTR2),1)
               
C
C-------------------------------------------------------------------
C              Calculate total transiton moment and write to output.
C-------------------------------------------------------------------
C
               WORK(KTRMOM+IEXCI-1) = WORK(KTRM1E+IEXCI-1) 
     &                              + WORK(KTRM2E+IEXCI-1)
     &                              + WORK(KTRM1D+IEXCI-1) 
     &                              + WORK(KTRM2D+IEXCI-1)
C
               IF ( IPRSOP .GE. 4) THEN
C
                  WRITE(LUPRI,'(/,1X,A,A,I3,A)')
     &            LABEL,' transition moment for state',IEXCI,':'
                  WRITE(LUPRI,'(4(F10.5,A),F10.5)')
     &            WORK(KTRM1E+IEXCI-1),' + ',WORK(KTRM2E+IEXCI-1),' + ',
     &            WORK(KTRM1D+IEXCI-1),' + ',WORK(KTRM2D+IEXCI-1),
     &            ' = ',WORK(KTRMOM)
C
               END IF
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
            CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
            CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
C
C-------------------------------------------------------
C           Write tabel of transition moments to output.
C-------------------------------------------------------
C
            IF ( IPRSOP .GE. 1) THEN
C
               IF (MODEL .EQ. 'DCRPA') THEN
                  WRITE(LUPRI,'(/,1X,A,A,A,I3)') LABEL,
     &                  ' RPA(D) transition moments for excitation',
     &                  ' symmetry',ISYMTR
               ELSE IF (MODEL .EQ. 'AOSOP') THEN
                  WRITE(LUPRI,'(/,1X,A,A,A,I3)') LABEL,
     &                  ' SOPPA transition moments for excitation',
     &                  ' symmetry',ISYMTR
               ELSE IF (MODEL .EQ. 'AOSOC') THEN
                  WRITE(LUPRI,'(/,1X,A,A,A,I3)') LABEL,
     &                 ' SOPPA(CCSD) transition moments for excitation',
     &                 ' symmetry',ISYMTR
               ELSE
                  WRITE(LUPRI,'(A)') 'ERROR: Wrong model in SO_TRMOM'
                  CALL QUIT('ERROR: Wrong model in SO_TRMOM')
               END IF
C
               WRITE(LUPRI,9001) 
               WRITE(LUPRI,'(A)')
     &         ' Exci.    1p-1h        2p-2h        1h-1p'//
     &         '        2h-2p        Total'
               WRITE(LUPRI,9002) 
C
               DO 150 IEXCI = 1,NEXCI
C
                  WRITE(LUPRI,'(1X,I3,2X,4(F10.5,A),F10.5)') IEXCI,
     &            WORK(KTRM1E+IEXCI-1),' + ',WORK(KTRM2E+IEXCI-1),' + ',
     &            WORK(KTRM1D+IEXCI-1),' + ',WORK(KTRM2D+IEXCI-1),
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
      CALL QEXIT('SO_TRMOM')
C
      RETURN
C
 9001 FORMAT(1X,'================================================',
     &       '===================')
 9002 FORMAT(1X,'------------------------------------------------',
     &       '-------------------')
C
      END
