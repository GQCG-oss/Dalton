      SUBROUTINE SO_PROP(DOUBLES,NOLDTR,NNEWTR,ISYMTR,IMAGPROP,
     &                   REDZ,LREDZ,LABEL,FREQ,
     &                   PROP,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     A.L. January 2004
C
C     PURPOSE: Calculates the value of properties befor trial vectors
C     orthonomarlization and write to a file the (1Z 2Z) 
C     solution vector
C
C     NOLDTR    number of trial vectors keep from the preciding iteration
C     NNEWTR    number of new trial vectors added in this iteration
C     ISYMTR    symmetry of the properties      
C     REDZ      REDZ(LREDZ) will be the (1Z^R 2Z) array
C     LREDZ     dimension of the REDZ array 2*(NNEWTR + NOLDTR)
C     LABEL     lablel of the propoerty
C     FREQ      value of the frequency
C     WORK
C     LWORK
C
#include <implicit.h>
#include <soppinf.h>
#include <priunit.h>
#include <ccsdsym.h>
C
      PARAMETER (ZERO = 0.0D+00, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THREE = 3.0D0)
C
C---------------------------------
C     Dimensions of the arguments.
C---------------------------------
C
      DIMENSION WORK(LWORK),REDZ(LREDZ)
      LOGICAL   IMAGPROP, DOUBLES
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_PROP')
C
C---------------------------------
C     1. allocation of work space.
C---------------------------------
C
      NTRIAL   = NNEWTR + NOLDTR
      NTRIALH  = NTRIAL/2
C
      LGPVC1   = 2*NT1AM(ISYMTR)
C      LGPVC1   = NT1AM(ISYMTR)
C      LGPVC2   = NT2AM(ISYMTR)
      LGPVC1H  = NT1AM(ISYMTR)
      IF (DOUBLES) THEN
         LGPVC2   = 2*NT2AM(ISYMTR)
         LGPVC2H  = NT2AM(ISYMTR)
      ELSE ! set doubles lenghts to zero 
         LGPVC2  = 0
         LGPVC2H = 0
      ENDIF
C
      KSOVC1   = 1
      KSOVC2   = KSOVC1 + LGPVC1 
      KGPVC1   = KSOVC2 + LGPVC2
      KGPVC2   = KGPVC1 + LGPVC1H 
      KSOVC1E  = KGPVC2 + LGPVC2H
      KSOVC1D  = KSOVC1E + LGPVC1H
      KSOVC2E  = KSOVC1D + LGPVC1H
      KSOVC2D  = KSOVC2E + LGPVC2H
      KEND1    = KSOVC2D + LGPVC2H
      LWORK1   = LWORK    - KEND1
C
      CALL SO_MEMMAX('SO_PROP',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_PROP',' ',KEND1,LWORK)
C
C-------------------------------------------
C     Open trial and solution vectors files.
C-------------------------------------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LGPVC1H)
      CALL SO_OPEN(LUTR1D,FNTR1D,LGPVC1H)
C
C     Initialize unit numbers before GPOPEN
      LUSV1E = -1      
      CALL GPOPEN(LUSV1E,FNSV1E,'UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      LUSV1D = -1
      CALL GPOPEN(LUSV1D,FNSV1D,'UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      IF(DOUBLES) THEN
         CALL SO_OPEN(LUTR2E,FNTR2E,LGPVC2H)
         CALL SO_OPEN(LUTR2D,FNTR2D,LGPVC2H)
         LUSV2E = -1
         CALL GPOPEN(LUSV2E,FNSV2E,'UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         LUSV2D = -1
         CALL GPOPEN(LUSV2D,FNSV2D,'UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
      ENDIF
C
C      
C      CALL DZERO(WORK(KSOVC1E),LGPVC1H)
C      CALL DZERO(WORK(KSOVC1D),LGPVC1H)
      CALL DZERO(WORK(KSOVC1),LGPVC1)
C
C      CALL DZERO(WORK(KSOVC2E),LGPVC2H)
C      CALL DZERO(WORK(KSOVC2D),LGPVC2H)
      IF(DOUBLES) CALL DZERO(WORK(KSOVC2),LGPVC2)
C
      DO 100 ITRIAL=1,NTRIAL
C
        CALL SO_READ(WORK(KSOVC1E),LGPVC1H,LUTR1E,FNTR1E,ITRIAL)
        CALL SO_READ(WORK(KSOVC1D),LGPVC1H,LUTR1D,FNTR1D,ITRIAL)
        IF (DOUBLES) THEN
           CALL SO_READ(WORK(KSOVC2E),LGPVC2H,LUTR2E,FNTR2E,ITRIAL)
           CALL SO_READ(WORK(KSOVC2D),LGPVC2H,LUTR2D,FNTR2D,ITRIAL)
        ENDIF
C
        IF(IPRSOP.GT.10) THEN
          CALL AROUND('Trial Vectors')
          CALL OUTPUT(WORK(KSOVC1E),1,LGPVC1H,1,2,LGPVC1H,2,1,LUPRI)
!          CALL OUTPUT(WORK(KSOVC1D),1,LGPVC1H,1,2,LGPVC1H,2,1,LUPRI)
          IF(DOUBLES)
     &        CALL OUTPUT(WORK(KSOVC2E),1,LGPVC2H,1,2,LGPVC2H,2,1,LUPRI)
C          CALL OUTPUT(WORK(KSOVC2D),1,LGPVC2H,1,2,LGPVC2H,2,1,LUPRI)
cLig DEBUG <<
cLig          WRITE(LUPRI,'(1X,A,F10.8)') 'C1: ',REDZ(ITRIAL)
cLig          WRITE(LUPRI,'(1X,A,F10.8)') 'C2: ',REDZ(ITRIAL+NTRIAL)
cLig DEBUG >>
        ENDIF
cLig  (1b 2b)i * Zi
CRF     Contribution of the straigt vectors 
        CALL DAXPY(LGPVC1H,REDZ(ITRIAL),WORK(KSOVC1E),1,
     &             WORK(KSOVC1),1)
        CALL DAXPY(LGPVC1H,REDZ(ITRIAL),WORK(KSOVC1D),1,
     &             WORK(KSOVC1+LGPVC1H),1)
        IF(DOUBLES)THEN
           CALL DAXPY(LGPVC2H,REDZ(ITRIAL),WORK(KSOVC2E),1,
     &                WORK(KSOVC2),1)
           CALL DAXPY(LGPVC2H,REDZ(ITRIAL),WORK(KSOVC2D),1,
     &                WORK(KSOVC2+LGPVC2H),1)
        ENDIF
cLig  (2b 1b)i * Zi+1/2 ...
CRF     Contribution of the paired vectors
        CALL DAXPY(LGPVC1H,REDZ(ITRIAL+NTRIAL),WORK(KSOVC1D),1,
     &             WORK(KSOVC1),1)
        CALL DAXPY(LGPVC1H,REDZ(ITRIAL+NTRIAL),WORK(KSOVC1E),1,
     &             WORK(KSOVC1+LGPVC1H),1)
        IF(DOUBLES)THEN
           CALL DAXPY(LGPVC2H,REDZ(ITRIAL+NTRIAL),WORK(KSOVC2D),1,
     &                WORK(KSOVC2),1)
           CALL DAXPY(LGPVC2H,REDZ(ITRIAL+NTRIAL),WORK(KSOVC2E),1,
     &                WORK(KSOVC2+LGPVC2H),1)
        ENDIF
C
100   CONTINUE
      IF(IPRSOP.GT.5) THEN
        CALL AROUND('Solution Vector')
        CALL OUTPUT(WORK(KSOVC1),1,LGPVC1H,1,2,LGPVC1H,2,1,LUPRI)
        IF(DOUBLES) 
     &    CALL OUTPUT(WORK(KSOVC2),1,LGPVC2H,1,2,LGPVC2H,2,1,LUPRI)
      ENDIF
C
C------------------------------------
C     Write Solution Vector on a file
C------------------------------------
C
C     find the right positions on SV files
      CALL SO_FNDPOS(LABEL,FREQ,LUSV1E)
      CALL SO_FNDPOS(LABEL,FREQ,LUSV1D)
C     write on files
      CALL SO_WRTVE(WORK(KSOVC1),LGPVC1H,ISYMTR,LABEL,
     &             FREQ,LUSV1E)
C
      CALL SO_WRTVE(WORK(KSOVC1+LGPVC1H),LGPVC1H,ISYMTR,LABEL,
     &             FREQ,LUSV1D)
C
      IF(DOUBLES)THEN
         CALL SO_FNDPOS(LABEL,FREQ,LUSV2E)
         CALL SO_FNDPOS(LABEL,FREQ,LUSV2D)
         CALL SO_WRTVE(WORK(KSOVC2),LGPVC2H,ISYMTR,LABEL,
     &                 FREQ,LUSV2E)
C
         CALL SO_WRTVE(WORK(KSOVC2+LGPVC2H),LGPVC2H,ISYMTR,LABEL,
     &                 FREQ,LUSV2D)
      ENDIF
C
cLig <> get the GP vec DEBUG ... I don't need it to compute the
cLig    solution vecotr
CRF all this can be skipped
C
      CALL SO_OPEN(LUGPV1,FNGPV1,LGPVC1)
      CALL SO_READ(WORK(KGPVC1),LGPVC1H,LUGPV1,FNGPV1,1)
      CALL SO_CLOSE(LUGPV1,FNGPV1,'KEEP')
      IF(DOUBLES)THEN
         CALL SO_OPEN(LUGPV2,FNGPV2,LGPVC2)
         CALL SO_READ(WORK(KGPVC2),LGPVC2H,LUGPV2,FNGPV2,1)
         CALL SO_CLOSE(LUGPV2,FNGPV2,'KEEP')
      ENDIF
cLig DEBUG <> compute the property
CRF generate for d-ex
      DFACTOR = -1.0D0
      IF (IMAGPROP) DFACTOR = 1.0D0      
      PROP1  =      DDOT(LGPVC1H,WORK(KGPVC1),1,WORK(KSOVC1),1)
     & + DFACTOR * DDOT(LGPVC1H,WORK(KGPVC1),1,WORK(KSOVC1+LGPVC1H),1) 
      PROP = PROP1
      IF(DOUBLES) THEN
         PROP2 =      DDOT(LGPVC2H,WORK(KGPVC2),1,WORK(KSOVC2),1)
     &   + DFACTOR * DDOT(LGPVC2H,WORK(KGPVC2),1,WORK(KSOVC2+LGPVC2H),1)
         PROP = PROP + PROP2
      ENDIF
C
C
C-----------------
C     Close files.
C-----------------
C
      CALL GPCLOSE(LUSV1E,'KEEP')
      CALL GPCLOSE(LUSV1D,'KEEP')
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
      IF(DOUBLES)THEN
         CALL GPCLOSE(LUSV2E,'KEEP')
         CALL GPCLOSE(LUSV2D,'KEEP')
         CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
         CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
      ENDIF
cLig DEBUG <<
C
C-----------------------
C        Print property
C-----------------------
C
      IF( IPRSOP .GT. 5) THEN
C
         CALL AROUND('Second order Property')
C
         WRITE(LUPRI,'(F10.8)') PROP
         IF(DOUBLES) THEN
            WRITE(LUPRI,'(/,1X,A,F10.8)') '1p1h contribution: ',PROP1
            WRITE(LUPRI,'(1X,A,F10.8)') '2p2h contribution: ',PROP2
         ENDIF
C
      ENDIF
C
cLig DEBUG <<
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_PROP')
C
      RETURN
      END
