C
C  /* Deck rp_prop */
      SUBROUTINE RP_PROP(NOLDTR,NNEWTR,ISYMTR,REDZ,LREDZ,LABEL,FREQ,
     &                   PROP,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     A.L. November 2003
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
#include "implicit.h"
#include "soppinf.h"
#include "priunit.h"
#include "ccsdsym.h"
C
      PARAMETER (ZERO = 0.0D+00, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THREE = 3.0D0)
C
C---------------------------------
C     Dimensions of the arguments.
C---------------------------------
C
      DIMENSION REDZ(LREDZ)
      DIMENSION WORK(LWORK)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_PROP')
C
C---------------------------------
C     1. allocation of work space.
C---------------------------------
C
      NTRIAL   = NNEWTR + NOLDTR
C
      LGPVEC   = 2*NT1AM(ISYMTR)
      LGPVCH   = NT1AM(ISYMTR)
C
      KSOVEC   = 1
      KGPVEC   = KSOVEC + LGPVEC
      KSOVEC1E = KGPVEC + LGPVEC
      KSOVEC1D = KSOVEC1E + LGPVCH
      KEND1    = KSOVEC1D + LGPVCH
      LWORK1   = LWORK    - KEND1
C
      CALL SO_MEMMAX('RP_PROP',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('RP_PROP',' ',KEND1,LWORK)
C
C-------------------
C     Open trial vectors files.
C-------------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LGPVCH)
      CALL SO_OPEN(LUTR1D,FNTR1D,LGPVCH)
C
C     Initialize unit numbers before GPOPEN
      LUSV1E = -1      
      CALL GPOPEN(LUSV1E,FNSV1E,'UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      LUSV1D = -1
      CALL GPOPEN(LUSV1D,FNSV1D,'UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
C
cLig <> get the GP vec DEBUG ... I don't need it to compute the
cLig    solution vecotr
C
      CALL SO_OPEN(LUGPVE,FNGPVE,LGPVEC)
      CALL SO_READ(WORK(KGPVEC),LGPVEC,LUGPVE,FNGPVE,1)
      CALL SO_CLOSE(LUGPVE,FNGPVE,'KEEP')
C      
      CALL DZERO(WORK(KSOVEC1E),LGPVCH)
      CALL DZERO(WORK(KSOVEC1D),LGPVCH)
      CALL DZERO(WORK(KSOVEC),LGPVEC)
C
      DO 100 ITRIAL=1,NTRIAL
C
        CALL SO_READ(WORK(KSOVEC1E),LGPVCH,LUTR1E,FNTR1E,ITRIAL)
        CALL SO_READ(WORK(KSOVEC1D),LGPVCH,LUTR1D,FNTR1D,ITRIAL)
C
        IF(IPRSOP.GT.10) THEN
          CALL AROUND('Trial Vectors')
          CALL OUTPUT(WORK(KSOVEC1E),1,LGPVCH,1,2,LGPVCH,2,1,LUPRI)
cLig          CALL OUTPUT(WORK(KSOVEC1D),1,LGPVCH,1,2,LGPVCH,2,1,LUPRI)
        ENDIF
C
cLig  (1b 2b)i * Zi
        CALL DAXPY(LGPVCH,REDZ(ITRIAL),WORK(KSOVEC1E),1,
     &             WORK(KSOVEC),1)
        CALL DAXPY(LGPVCH,REDZ(ITRIAL),WORK(KSOVEC1D),1,
     &             WORK(KSOVEC+LGPVCH),1)
cLig  (2b 1b)i * Zi+1/2 ...
        CALL DAXPY(LGPVCH,REDZ(ITRIAL+NTRIAL),WORK(KSOVEC1D),1,
     &             WORK(KSOVEC),1)
        CALL DAXPY(LGPVCH,REDZ(ITRIAL+NTRIAL),WORK(KSOVEC1E),1,
     &             WORK(KSOVEC+LGPVCH),1)
C
100   CONTINUE
C
      IF(IPRSOP.GT.5) THEN
        CALL AROUND('Solution Vector')
        CALL OUTPUT(WORK(KSOVEC),1,LGPVCH,1,2,LGPVCH,2,1,LUPRI)
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
      CALL SO_WRTVE(WORK(KSOVEC),LGPVCH,ISYMTR,LABEL,
     &             FREQ,LUSV1E)
C
      CALL SO_WRTVE(WORK(KSOVEC+LGPVCH),LGPVCH,ISYMTR,LABEL,
     &             FREQ,LUSV1D)
C
cLig DEBUG <> compute the property
      PROP = DDOT(LGPVEC,WORK(KGPVEC),1,WORK(KSOVEC),1)
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
         WRITE(LUPRI,'(1X,A,A,F12.5,A,I3,A,F10.8)') LABEL,' FREQ ',
     &         FREQ,' SYMM',ISYMTR,' VALUE ', PROP
C
      ENDIF
C
cLig DEBUG <<
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('RP_PROP')
C
      RETURN
      END
