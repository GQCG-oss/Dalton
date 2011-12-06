!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
#include <single.h>
C
      SUBROUTINE CRHYP(VECA,VECB,VECC,VECD,VECBC,VECBD,VECCD,
     *                 RESVEC,CMO,UDV,PV,FOCK,FC,FV,
     *                 XINDX,MJWOP,WRK,LWRK)
C
#include <implicit.h>
#include <dummy.h>
#include <iratdef.h>
C
C PURPOSE:
C CALCULATION OF SECOND HYPERPOLARIZABILITIES
C PROPERTIES
C
      LOGICAL DOHYP, ATEST, DIPLEN, GAMFLG
      DATA ATEST/.FALSE./
C
      CHARACTER*8 ALAB,BLAB,CLAB,DLAB,BLANK
      PARAMETER   (BLANK='        ')
C
      DIMENSION VECA(*),VECB(*),VECC(*),VECD(*)
      DIMENSION VECBC(*),VECCD(*),VECBD(*)
      DIMENSION RESVEC(*), GAMMA(3,3,3,3)
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*)
      DIMENSION XINDX(*),WRK(*)
C
      PARAMETER ( D0 = 0.0D0 , THRFRQ = 1.0D-14 )
C
#include <priunit.h>
#include <infrsp.h>
#include <maxorb.h>
#include <infvar.h>
#include <inforb.h>
#include <infpri.h>
#include <infspi.h>
#include <wrkrsp.h>
#include <tstjep.h>
#include <infhso.h>
#include <qrinf.h>
      DIMENSION MJWOP(2,MAXWOP,8)
#include <rspprp.h>
#include <indcr.h>
#include <infcr.h>
#include <infinp.h>
C
      CALL QENTER('CRHYP')
      WRITE(LUPRI,'(//A,A)')
     *' <<<<< CALCULATING CONTRIBUTIONS TO SECOND',
     *' HYPERPOLARIZABILITY >>>>>'
      CALL FLSHFO(LUPRI)
C
C     We write the results to a complementary output file. This will then
C     both serve as a file for getting a summary of the results, but more
C     importantly, it will serve as a way of avoiding already calculated
C     gamma-components in case of a crashed calculation. Check for calculated
C     components are done in BCDCHK.
C
      LURSPRES = -1
      CALL GPOPEN(LURSPRES,'RESULTS.RSP','UNKNOWN',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
C
      DO 200 IDFR = 1,NDCRFR
      DO 300 ICFR = 1,NCCRFR
      DO 400 IBFR = 1,NBCRFR
C
C
C Set GAMMA to ZERO
C
      DO L=1,3
      DO K=1,3
      DO J=1,3
      DO I=1,3
         GAMMA(I,J,K,L) = D0
      END DO
      END DO
      END DO
      END DO
C
      DO 500 ISYMD = 1,NSYM
      DO 550 ISYMC = 1,NSYM
      DO 600 ISYMB = 1,NSYM
C
      ISYMA = MULD2H(ISYMD,MULD2H(ISYMC,ISYMB))
      IF ( (NDCROP(ISYMD).GT.0) .AND. (NCCROP(ISYMC).GT.0) .AND.
     *     (NBCROP(ISYMB).GT.0) .AND. (NACROP(ISYMA).GT.0) ) THEN
C
C  If a special process has been specified, we only need certain
C  combinations of the specified frequencies
C
      IF (CRSPEC) THEN
         IF (CRKERR) THEN
            IF (ABS(CCRFR(ICFR)).LE.THRFRQ .AND.
     *          ABS(DCRFR(IDFR)).LE.THRFRQ) GO TO 649
         ENDIF
         IF (CRSHG) THEN
            DIFFRQ = BCRFR(IBFR) - CCRFR(ICFR)
            IF (ABS(DIFFRQ).LE.THRFRQ .AND. ABS(DCRFR(IDFR)).LE.THRFRQ)
     *         GO TO 649
         ENDIF
         IF (CRIDRI) THEN
            DIFFRQ = BCRFR(IBFR) - CCRFR(ICFR)
            IF (ABS(DIFFRQ).LE.THRFRQ) THEN
               DIFFRQ = BCRFR(IBFR) + DCRFR(IDFR)
               IF (ABS(DIFFRQ).LE.THRFRQ) GO TO 649
            ENDIF
         ENDIF
         IF (CRTHG) THEN
            DIFFRQ = BCRFR(IBFR) - CCRFR(ICFR)
            IF (ABS(DIFFRQ).LE.THRFRQ) THEN
               DIFFRQ = BCRFR(IBFR) - DCRFR(IDFR)
               IF (ABS(DIFFRQ).LE.THRFRQ) GO TO 649
            ENDIF
         ENDIF
         GO TO 400
      ENDIF
 649  CONTINUE
C
      DO 650 IDOP = 1,NDCROP(ISYMD)
      DO 700 ICOP = 1,NCCROP(ISYMC)
      DO 750 IBOP = 1,NBCROP(ISYMB)
      DO 800 IAOP = 1,NACROP(ISYMA)
C
C     Initialize variables.
C     Check if an equivalent gamma calculation already has been done,
C     DOHYP indicates the result. Read response vectors from disk.
C     Check if some of the response vectors are equal or zero,
C     IBCDEQ indicates the result
C
      CALL BCDCHK(DOHYP,IBCDEQ,LURSPRES,DIPLEN,GAMMA,
     *            ISYMA,ISYMB,ISYMC,ISYMD,ISYMBC,ISYMBD,ISYMCD,
     *            ALAB,BLAB,CLAB,DLAB,IAOP,IBOP,ICOP,IDOP,
     *            IBFR,ICFR,IDFR,FREQA,FREQB,FREQC,FREQD,
     *            KZYVA,KZYVB,KZYVC,KZYVD,KZYVBC,KZYVBD,KZYVCD,
     *            VECA,VECB,VECC,VECD,VECBC,VECBD,VECCD,IDUM)
C
      IF (.NOT.DOHYP) GOTO 800
C
C    Initialize second hyperpolarizability
C
      HYPVAL = 0
C
      IF (IPRRSP.GT.0) WRITE(LUPRI,'(//A15,2A20,/A)')
     *   'Contribution','Term','Accumulated',
     *   ' ------------------------------------------------------'
C
C
C     Calculate Na T[4] Nb Nc Nd
C
C
      IF (IPRRSP.GT.5) CALL TIMER('START ',TIMSTR,TIMEND)
C
      CALL T4DRV(IBCDEQ,ISYMA,ISYMB,ISYMC,ISYMD,VECA,VECB,VECC,VECD,
     *           -FREQB,-FREQC,-FREQD,XINDX,UDV,PV,MJWOP,WRK,LWRK,
     *           CMO,FC)
      VAL = DDOT(KZYVA,WRK,1,VECA,1)
      HYPVAL = HYPVAL + VAL
C
      IF (IPRRSP.GT.0)
     * WRITE(LUPRI,'(A17,F18.8,F20.8)')' Na T[4] Nb Nc Nd',VAL,HYPVAL
C
      IF (IPRRSP.GT.5) CALL TIMER('T4DRV ',TIMSTR,TIMEND)
C
C
C     Calculate Na T[3] Nb Ncd type terms (three permutations)
C
C
      IF (IPRRSP.GT.5) CALL TIMER('START ',TIMSTR,TIMEND)
C
      CALL T3DRV(1,ISYMA,ISYMB,ISYMCD,VECB,VECCD,ATEST,VECA,
     *           -FREQB,-FREQC-FREQD,XINDX,UDV,PV,MJWOP,
     &           WRK,LWRK,CMO,FC,FV)
      VAL = -DDOT(KZYVA,WRK,1,VECA,1)
      TMPVAL = VAL
      HYPVAL = HYPVAL + VAL
C
      IF (IBCDEQ.EQ.24) THEN
         TMPVAL = TMPVAL + 2*VAL
         HYPVAL = HYPVAL + 2*VAL
      ELSE IF (IBCDEQ.EQ.2) THEN
         TMPVAL = TMPVAL + VAL
         HYPVAL = HYPVAL + VAL
         CALL T3DRV(1,ISYMA,ISYMD,ISYMBC,VECD,VECBC,ATEST,VECA,
     *        -FREQD,-FREQB-FREQC,XINDX,UDV,PV,MJWOP,
     &        WRK,LWRK,CMO,FC,FV)
         VAL = -DDOT(KZYVA,WRK,1,VECA,1)
         TMPVAL = TMPVAL + VAL
         HYPVAL = HYPVAL + VAL
      ELSE IF (IBCDEQ.EQ.3) THEN
         TMPVAL = TMPVAL + VAL
         HYPVAL = HYPVAL + VAL
         CALL T3DRV(1,ISYMA,ISYMC,ISYMBD,VECC,VECBD,ATEST,VECA,
     *              -FREQC,-FREQB-FREQD,XINDX,UDV,PV,MJWOP,
     &              WRK,LWRK,CMO,FC,FV)
         VAL = -DDOT(KZYVA,WRK,1,VECA,1)
         TMPVAL = TMPVAL + VAL
         HYPVAL = HYPVAL + VAL
      ELSE IF (IBCDEQ.EQ.4) THEN
         CALL T3DRV(1,ISYMA,ISYMC,ISYMBD,VECC,VECBD,ATEST,VECA,
     *              -FREQC,-FREQB-FREQD,XINDX,UDV,PV,MJWOP,
     &              WRK,LWRK,CMO,FC,FV)
         VAL = -DDOT(KZYVA,WRK,1,VECA,1)
         TMPVAL = TMPVAL + 2*VAL
         HYPVAL = HYPVAL + 2*VAL
      ELSE
         CALL T3DRV(1,ISYMA,ISYMC,ISYMBD,VECC,VECBD,ATEST,VECA,
     *              -FREQC,-FREQB-FREQD,XINDX,UDV,PV,MJWOP,
     &              WRK,LWRK,CMO,FC,FV)
         VAL = -DDOT(KZYVA,WRK,1,VECA,1)
         TMPVAL = TMPVAL + VAL
         HYPVAL = HYPVAL + VAL
         CALL T3DRV(1,ISYMA,ISYMD,ISYMBC,VECD,VECBC,ATEST,VECA,
     *              -FREQD,-FREQB-FREQC,XINDX,UDV,PV,MJWOP,
     &              WRK,LWRK,CMO,FC,FV)
         VAL = -DDOT(KZYVA,WRK,1,VECA,1)
         TMPVAL = TMPVAL + VAL
         HYPVAL = HYPVAL + VAL
      END IF
C
      IF (IPRRSP.GT.0)
     * WRITE(LUPRI,'(A15,2F20.8)')' Na T[3] Nx Nyz',TMPVAL,HYPVAL
C
      IF (IPRRSP.GT.5) CALL TIMER('T3DRV ',TIMSTR,TIMEND)

C
C     Calculate Na B[3] Nc Nd type terms 
C     (two of six permutations in each call)
C
      IF (IPRRSP.GT.5) CALL TIMER('START ',TIMSTR,TIMEND)
C
      CALL X3INIT(KZYVA,KZYVC,KZYVD,ISYMA,ISYMC,ISYMD,BLAB,
     *            ISYMB,VECC,VECD,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVA,RESVEC,1,VECA,1)
      TMPVAL = VAL
      HYPVAL = HYPVAL + VAL
C
      CALL X3INIT(KZYVA,KZYVB,KZYVD,ISYMA,ISYMB,ISYMD,CLAB,
     *            ISYMC,VECB,VECD,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVA,RESVEC,1,VECA,1)
      HYPVAL = HYPVAL + VAL
      TMPVAL = TMPVAL + VAL
C
      CALL X3INIT(KZYVA,KZYVB,KZYVC,ISYMA,ISYMB,ISYMC,DLAB,
     *            ISYMD,VECB,VECC,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVA,RESVEC,1,VECA,1)
      TMPVAL = TMPVAL + VAL
      HYPVAL = HYPVAL + VAL
C
      IF (IPRRSP.GT.0)
     *  WRITE(LUPRI,'(A15,3F20.8)')' Na X[3] Ny Nz ',TMPVAL,HYPVAL
C
C
C     Calculate Nb A[3] Nc Nd type terms
C     (two of six permutations in each call)
C
      CALL A3INIT(KZYVB,KZYVC,KZYVD,ISYMB,ISYMC,ISYMD,ALAB,
     *            ISYMA,VECC,VECD,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVB,RESVEC,1,VECB,1)
      TMPVAL = VAL
      HYPVAL = HYPVAL + VAL
C
      CALL A3INIT(KZYVC,KZYVB,KZYVD,ISYMC,ISYMB,ISYMD,ALAB,
     *            ISYMA,VECB,VECD,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVC,RESVEC,1,VECC,1)
      TMPVAL = TMPVAL + VAL
      HYPVAL = HYPVAL + VAL
C
      CALL A3INIT(KZYVD,KZYVB,KZYVC,ISYMD,ISYMB,ISYMC,ALAB,
     *            ISYMA,VECB,VECC,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVD,RESVEC,1,VECD,1)
      TMPVAL = TMPVAL + VAL
      HYPVAL = HYPVAL + VAL
C
      IF (IPRRSP.GT.0)
     *  WRITE(LUPRI,'(A15,2F20.8)')' Nx A[3] Ny Nz ',TMPVAL,HYPVAL
C
C
C     Calculate Na B[2] Ncd type terms (three permutations)
C
C
      CALL X2INIT(1,KZYVA,KZYVCD,ISYMA,ISPINA,ISYMCD,0,WRK(1),VECCD,
     *            RESVEC,XINDX,UDV,PV,BLAB,ISYMB,0,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVA,RESVEC,1,VECA,1)
      TMPVAL = VAL
      HYPVAL = HYPVAL + VAL
C
      CALL X2INIT(1,KZYVA,KZYVBD,ISYMA,ISPINA,ISYMBD,0,1,VECBD,
     *            RESVEC,XINDX,UDV,PV,CLAB,ISYMC,ISPINC,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVA,RESVEC,1,VECA,1)
      TMPVAL = TMPVAL + VAL
      HYPVAL = HYPVAL + VAL
C
      CALL X2INIT(1,KZYVA,KZYVBC,ISYMA,ISPINA,ISYMBC,0,1,VECBC,
     *            RESVEC,XINDX,UDV,PV,DLAB,ISYMD,ISPIND,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVA,RESVEC,1,VECA,1)
      TMPVAL = TMPVAL + VAL
      HYPVAL = HYPVAL + VAL
C
      IF (IPRRSP.GT.0)
     * WRITE(LUPRI,'(A15,2F20.8)') ' Na X[2] Nyz   ',TMPVAL,HYPVAL
C
C
C     Calculate Nb A[2] Ncd type terms (six permutations)
C
C
C     CALL DCOPY(KZYVB,VECB,1,RESVEC)
      CALL A2INIT(1,KZYVB,KZYVCD,ISYMB,ISPINB,ISYMCD,0,1,VECCD,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVB,RESVEC,1,VECB,1)
      TMPVAL = VAL
      HYPVAL = HYPVAL + VAL
C
C     CALL DCOPY(KZYVCD,VECCD,1,RESVEC)
      CALL A2INIT(1,KZYVCD,KZYVB,ISYMCD,0,ISYMB,ISPINB,1,VECB,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVCD,RESVEC,1,VECCD,1)
      TMPVAL = TMPVAL + VAL
      HYPVAL = HYPVAL + VAL
C
C     CALL DCOPY(KZYVC,VECC,1,RESVEC)
      CALL A2INIT(1,KZYVC,KZYVBD,ISYMC,ISPINC,ISYMBD,0,1,VECBD,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVC,RESVEC,1,VECC,1)
      TMPVAL = TMPVAL + VAL
      HYPVAL = HYPVAL + VAL
C
C     CALL DCOPY(KZYVBD,VECBD,1,RESVEC)
      CALL A2INIT(1,KZYVBD,KZYVC,ISYMBD,0,ISYMC,ISPINC,1,VECC,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVBD,RESVEC,1,VECBD,1)
      TMPVAL = TMPVAL + VAL
      HYPVAL = HYPVAL + VAL
C
C     CALL DCOPY(KZYVD,VECD,1,RESVEC)
      CALL A2INIT(1,KZYVD,KZYVBC,ISYMD,ISPIND,ISYMBC,0,1,VECBC,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVD,RESVEC,1,VECD,1)
      TMPVAL = TMPVAL + VAL
      HYPVAL = HYPVAL + VAL
C
C     CALL DCOPY(KZYVBC,VECBC,1,RESVEC)
      CALL A2INIT(1,KZYVBC,KZYVD,ISYMBC,0,ISYMD,ISPIND,1,VECD,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVBC,RESVEC,1,VECBC,1)
      TMPVAL = TMPVAL + VAL
      HYPVAL = HYPVAL + VAL
C
      IF (IPRRSP.GT.0)
     * WRITE(LUPRI,'(A15,2F20.8)') ' Nx A[2] Nyz   ',TMPVAL,HYPVAL
C
      IF (IPRRSP.GT.5) CALL TIMER('OTHERS',TIMSTR,TIMEND)
C
      WRITE(LUPRI,'(/A,4(/A,A10,I4,F10.6))')
     &    '@ Cubic response function value in a.u. for',
     &    '@ A operator, symmetry, frequency: ',ALAB,ISYMA,-FREQA,
     &    '@ B operator, symmetry, frequency: ',BLAB,ISYMB,FREQB,
     &    '@ C operator, symmetry, frequency: ',CLAB,ISYMC,FREQC,
     &    '@ D operator, symmetry, frequency: ',DLAB,ISYMD,FREQD
      WRITE(LUPRI,'(/A,F20.8/)') '@ << A; B, C, D >>  = ', HYPVAL
      CALL FLSHFO(LUPRI)
C
C     Write to result file
C
      WRITE(LURSPRES,'(/A,4(/A,A10,I4,F10.6))')
     &    ' Cubic response function value in a.u. for',
     &    ' A operator, symmetry, frequency: ',ALAB,ISYMA,-FREQA,
     &    ' B operator, symmetry, frequency: ',BLAB,ISYMB,FREQB,
     &    ' C operator, symmetry, frequency: ',CLAB,ISYMC,FREQC,
     &    ' D operator, symmetry, frequency: ',DLAB,ISYMD,FREQD
      WRITE(LURSPRES,'(/A,F20.8)') ' << A; B, C, D >>  = ', HYPVAL
      CALL FLSHFO(LURSPRES)
C
C If dipole operators only, store result in GAMMA(A,B,C,D)
C     (hjaaj Oct 2001: note that
C      gamma = <<mu; -mu, -mu, -mu>> = - <<r; r, r, r>> )
C
      IF (DIPLEN) THEN
         GAMFLG = .TRUE.
         CALL DIPLAB(ALAB,I1)
         CALL DIPLAB(BLAB,I2)
         CALL DIPLAB(CLAB,I3)
         CALL DIPLAB(DLAB,I4)
         GAMMA(I1,I2,I3,I4) = -HYPVAL
      END IF
C
C Operator loops
C
 800  CONTINUE
 750  CONTINUE
 700  CONTINUE
 650  CONTINUE
C
      END IF
C
C Symmetry loops
C
 600  CONTINUE
 550  CONTINUE
 500  CONTINUE
C
      IF (GAMFLG) CALL GAMPRI(GAMMA,FREQB,FREQC,FREQD)
C
C Frequency loops
C
 400  CONTINUE
 300  CONTINUE
 200  CONTINUE
C
C    End of subroutine CRHYP
C
      CALL GPCLOSE(LURSPRES,'KEEP')
      CALL QEXIT('CRHYP')
      RETURN
      END
      SUBROUTINE READVE(ISYMA,ISYMB,ISYMC,ISYMD,
     *                  ALAB,BLAB,CLAB,DLAB,
     *                  FREQA,FREQB,FREQC,FREQD,
     *                  KZYVA,KZYVB,KZYVC,KZYVD,KZYVBC,KZYVBD,KZYVCD,
     *                  VECA,VECB,VECC,VECD,VECBC,VECBD,VECCD)
C
C PURPOSE: Read in response vectors. Called from BCDCHK.
C
#include <implicit.h>
#include <iratdef.h>
C
      LOGICAL FOUND, CONV
      CHARACTER*8 ALAB,BLAB,CLAB,DLAB,BLANK
      PARAMETER (BLANK='       ', D0=0.0D0)
      DIMENSION VECA(*),VECB(*),VECC(*),VECD(*)
      DIMENSION VECBC(*),VECCD(*),VECBD(*)
C
#include <infrsp.h>
#include <maxorb.h>
#include <priunit.h>
#include <infvar.h>
#include <inforb.h>
#include <infpri.h>
#include <infspi.h>
#include <wrkrsp.h>
#include <tstjep.h>
#include <infhso.h>
#include <inftap.h>
#include <rspprp.h>
#include <qrinf.h>
#include <indcr.h>
#include <infcr.h>
#include <inflr.h>
#include <inftmo.h>
C
      IF (IPRRSP.GT.100) THEN
         WRITE(LUPRI,'(3(/A),4I3,/A,4A10,/A,4F12.8)')
     *         '     Variables in READVE ',
     *         '     =================== ',
     *         'ISYMA,ISYMB,ISYMC,ISYMD: ',ISYMA,ISYMB,ISYMC,ISYMD,
     *         'ALAB,BLAB,CLAB,DLAB    : ',ALAB,BLAB,CLAB,DLAB,
     *         'FREQA,FREQB,FREQC,FREQD: ',FREQA,FREQB,FREQC,FREQD
      END IF
C
C     Read in Na
C
      CALL REARSP(LURSP,KZYVA,VECA,ALAB,BLANK,FREQA,D0,ISYMA,0,
     &            THCLR,FOUND,CONV,ANTSYM)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,F7.5,A,I3,/A)') ' Response label ',ALAB,
     &           ' with frequency ',FREQA, ' and symmetry',
     &           ISYMA,' not found on file LURSP'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,F7.5,/A,I3,A)') ' @WARNING>>>>'//
     &           ' Response label ',ALAB,
     &           ' with frequency ',FREQA, ' and symmetry',
     &           ISYMA,' not converged on file LURSP'
         END IF
      END IF
      IF (FREQA .LT. D0) THEN
         CALL DSWAP(KZYVA/2,VECA,1,VECA(1+KZYVA/2),1)
         IF (ANTSYM .LT. D0) CALL DSCAL(KZYVA,ANTSYM,VECA,1)
      END IF
C
      IF (IPRRSP.GT.100) THEN
         WRITE(LUPRI,'(/A/A)')
     *         '     Na  vector in READVE ',
     *         '     =================== '
         CALL OUTPUT(VECA,1,KZYVA/2,1,2,KZYVA/2,2,1,LUPRI)
      END IF
C
C     Read in Nb
C
      CALL REARSP(LURSP,KZYVB,VECB,BLAB,BLANK,FREQB,D0,ISYMB,0,
     &            THCLR,FOUND,CONV,ANTSYM)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,F7.5,A,I3,/A)') ' Response label ',BLAB,
     &           ' with frequency ',FREQB, ' and symmetry',
     &           ISYMB,' not found on file LURSP'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,F7.5,/A,I3,A)') ' @WARNING>>>>'//
     &           ' Response label ',BLAB,
     &           ' with frequency ',FREQB, ' and symmetry',
     &           ISYMB,' not converged on file LURSP'
         END IF
      END IF
      IF (FREQB .LT. D0) THEN
         CALL DSWAP(KZYVB/2,VECB,1,VECB(1+KZYVB/2),1)
         IF (ANTSYM .LT. D0) CALL DSCAL(KZYVB,ANTSYM,VECB,1)
      END IF
C
      IF (IPRRSP.GT.100) THEN
         WRITE(LUPRI,'(/A)') '     Nb  vector in READVE '
         WRITE(LUPRI,'(A)')  '     =================== '
         CALL OUTPUT(VECB,1,KZYVB/2,1,2,KZYVB/2,2,1,LUPRI)
      END IF
C
C     Read in Nc
C
      CALL REARSP(LURSP,KZYVC,VECC,CLAB,BLANK,FREQC,D0,ISYMC,0,
     &            THCLR,FOUND,CONV,ANTSYM)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,F7.5,A,I3,/A)') ' Response label ',CLAB,
     &           ' with frequency ',FREQC, ' and symmetry',
     &           ISYMC,' not found on file LURSP'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,F7.5,/A,I3,A)') ' @WARNING>>>>'//
     &           ' Response label ',CLAB,
     &           ' with frequency ',FREQC, ' and symmetry',
     &           ISYMC,' not converged on file LURSP'
         END IF
      END IF
      IF (FREQC .LT. D0) THEN
         CALL DSWAP(KZYVC/2,VECC,1,VECC(1+KZYVC/2),1)
         IF (ANTSYM .LT. D0) CALL DSCAL(KZYVC,ANTSYM,VECC,1)
      END IF
C
      IF (IPRRSP.GT.100) THEN
         WRITE(LUPRI,'(/A)') '     Nc  vector in READVE '
         WRITE(LUPRI,'(A)')  '     =================== '
         CALL OUTPUT(VECC,1,KZYVC/2,1,2,KZYVC/2,2,1,LUPRI)
      END IF
C
C     Read in Nd
C
      CALL REARSP(LURSP,KZYVD,VECD,DLAB,BLANK,FREQD,D0,ISYMD,0,
     &            THCLR,FOUND,CONV,ANTSYM)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,F7.5,A,I3,/A)') ' Response label ',DLAB,
     &           ' with frequency ',FREQD, ' and symmetry',
     &           ISYMD,' not found on file LURSP'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,F7.5,/A,I3,A)') ' @WARNING>>>>'//
     &           ' Response label ',DLAB,
     &           ' with frequency ',FREQD, ' and symmetry',
     &           ISYMD,' not converged on file LURSP'
         END IF
      END IF
      IF (FREQD .LT. D0) THEN
         CALL DSWAP(KZYVD/2,VECD,1,VECD(1+KZYVD/2),1)
         IF (ANTSYM .LT. D0) CALL DSCAL(KZYVD,ANTSYM,VECD,1)
      END IF
C
      IF (IPRRSP.GT.100) THEN
         WRITE(LUPRI,'(/A)') '     Nd  vector in READVE '
         WRITE(LUPRI,'(A)')  '     =================== '
         CALL OUTPUT(VECD,1,KZYVD/2,1,2,KZYVD/2,2,1,LUPRI)
      END IF
C
C     Read in Nbc
C
      CALL REARSP(LURSP,KZYVBC,VECBC,BLAB,CLAB,FREQB,FREQC,
     &            ISYMB,ISYMC,THCLR,FOUND,CONV,ANTSYM)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,1X,A,F7.5,1X,F7.5,/A,I3,1X,I3,A)') 
     &           ' Response labels ',BLAB,CLAB,
     &           ' with frequencies ',FREQB, FREQC,
     &           ' and symmetries ',ISYMB,ISYMC,
     &           ' not found on file LURSP'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,1X,A,F7.5,1X,F7.5,/A,I3,1X,I3,A)') 
     &           ' @WARNING>>>> Response labels ',BLAB,CLAB,
     &           ' with frequencies ',FREQB, FREQC,
     &           ' and symmetries ',ISYMB,ISYMC,
     &           ' not converged on file LURSP'
         END IF
      END IF
C
      IF (IPRRSP.GT.100) THEN
         WRITE(LUPRI,'(/A)') '     Nbc vector in READVE '
         WRITE(LUPRI,'(A)')  '     =================== '
         CALL OUTPUT(VECBC,1,KZYVBC/2,1,2,KZYVBC/2,2,1,LUPRI)
      END IF
C
C     Read in Nbd
C
      CALL REARSP(LURSP,KZYVBD,VECBD,BLAB,DLAB,FREQB,FREQD,
     &            ISYMB,ISYMD,THCLR,FOUND,CONV,ANTSYM)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,1X,A,F7.5,1X,F7.5,/A,I3,1X,I3,A)') 
     &           ' Response labels ',BLAB,DLAB,
     &           ' with frequencies ',FREQB, FREQD,
     &           ' and symmetries ',ISYMB,ISYMD,
     &           ' not found on file LURSP'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,1X,A,F7.5,1X,F7.5,/A,I3,1X,I3,A)') 
     &           ' @WARNING>>>> Response labels ',BLAB,DLAB,
     &           ' with frequencies ',FREQB, FREQD,
     &           ' and symmetries ',ISYMB,ISYMD,
     &           ' not converged on file LURSP'
         END IF
      END IF
C
      IF (IPRRSP.GT.100) THEN
         WRITE(LUPRI,'(/A)') '     Nbd vector in READVE '
         WRITE(LUPRI,'(A)')  '     =================== '
         CALL OUTPUT(VECBD,1,KZYVBD/2,1,2,KZYVBD/2,2,1,LUPRI)
      END IF
C
C     Read in Ncd
C
      CALL REARSP(LURSP,KZYVCD,VECCD,CLAB,DLAB,FREQC,FREQD,
     &            ISYMC,ISYMD,THCLR,FOUND,CONV,ANTSYM)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,1X,A,F7.5,1X,F7.5,/A,I3,1X,I3,A)') 
     &           ' Response labels ',CLAB,DLAB,
     &           ' with frequencies ',FREQC, FREQD,
     &           ' and symmetries ',ISYMC,ISYMD,
     &           ' not found on file LURSP'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,1X,A,F7.5,1X,F7.5,/A,I3,1X,I3,A)') 
     &           ' @WARNING>>>> Response labels ',CLAB,DLAB,
     &           ' with frequencies ',FREQC, FREQD,
     &           ' and symmetries ',ISYMC,ISYMD,
     &           ' not converged on file LURSP'
         END IF
      END IF
C
      IF (IPRRSP.GT.100) THEN
         WRITE(LUPRI,'(/A)') '     Ncd vector in READVE '
         WRITE(LUPRI,'(A)')  '     =================== '
         CALL OUTPUT(VECCD,1,KZYVCD/2,1,2,KZYVCD/2,2,1,LUPRI)
      END IF
C
      RETURN
      END
      SUBROUTINE GAMPRI(GAMMA,FREQB,FREQC,FREQD)
C
#include <implicit.h>
C
#include <maxorb.h>
#include <infrsp.h>
#include <priunit.h>
#include <infvar.h>
#include <inforb.h>
#include <infpri.h>
#include <infspi.h>
#include <wrkrsp.h>
#include <tstjep.h>
#include <rspprp.h>
#include <infhso.h>
#include <qrinf.h>
#include <indcr.h>
#include <infcr.h>
#include <inftmo.h>
C
      PARAMETER( D0=0.0D0 )
      CHARACTER*1 DIPOP
      CHARACTER*8 DIQMMM
      DIMENSION GAMMA(3,3,3,3), DIPOP(3), DIQMMM(3)
C
      DIPOP(1) = 'X'
      DIPOP(2) = 'Y'
      DIPOP(3) = 'Z'

      DIQMMM(1) = 'X       '
      DIQMMM(2) = 'Y       '
      DIQMMM(3) = 'Z       '

C
      CALL HEADER(' Summary of gamma values for a set of frequencies',0)
      WRITE(LUPRI,'(3(A,F10.6,/))') '@ B-freq:',FREQB,
     *     '@ C-freq:',FREQC,'@ D-freq:',FREQD
C
      IZERO = 0
C
      GAMAVE = D0
C
      DO I=1,3
      DO J=1,3
         IF (GAMMA(I,I,J,J) .EQ. D0) THEN
            IZERO = IZERO + 1
         ELSE
            WRITE(LUPRI,'(A,8(A1),F16.4)') 
     *           '@ gamma(',DIPOP(I),';',DIPOP(I),',',DIPOP(J),
     *           ',',DIPOP(J),')',GAMMA(I,I,J,J)

            CALL WRIPRO(GAMMA(I,I,J,J),"CR-SCF/DFT",4,
     *                  DIQMMM(I),DIQMMM(I),
     *                  DIQMMM(J),DIQMMM(J),
     *                  FREQB,FREQC,FREQD,
     *                  1,0,0,0)
         END IF
         IF (I.NE.J) THEN
            IF (GAMMA(I,J,J,I) .EQ. D0) THEN
               IZERO = IZERO + 1
            ELSE
               WRITE(LUPRI,'(A,8(A1),F16.4)') 
     *              '@ gamma(',DIPOP(I),';',DIPOP(J),',',DIPOP(J),
     *           ',',DIPOP(I),')',GAMMA(I,J,J,I)
            CALL WRIPRO(GAMMA(I,J,J,I),"CR-SCF/DFT",4,
     *                  DIQMMM(I),DIQMMM(J),
     *                  DIQMMM(J),DIQMMM(I),
     *                  FREQB,FREQC,FREQD,
     *                  1,0,0,0)

            END IF
            IF (GAMMA(I,J,I,J) .EQ. D0) THEN
               IZERO = IZERO + 1
            ELSE
               WRITE(LUPRI,'(A,8(A1),F16.4)') 
     *              '@ gamma(',DIPOP(I),';',DIPOP(J),',',DIPOP(I),
     *           ',',DIPOP(J),')',GAMMA(I,J,I,J)
            CALL WRIPRO(GAMMA(I,J,I,J),"CR-SCF/DFT",4,
     *                  DIQMMM(I),DIQMMM(J),
     *                  DIQMMM(I),DIQMMM(J),
     *                  FREQB,FREQC,FREQD,
     *                  1,0,0,0)
            END IF
         END IF
         GAMAVE=GAMAVE + 
     *        ( GAMMA(I,I,J,J) + 
     *          GAMMA(I,J,J,I) + 
     *          GAMMA(I,J,I,J) )/15
      END DO
      END DO
      IF (IZERO .EQ. 0) THEN
         WRITE(LUPRI,'(/2(A),F16.6)') '@ Averaged gamma parallel ',
     *     'to the applied field is ', GAMAVE
      END IF
C
C End of GAMRPI
C
      RETURN
      END
