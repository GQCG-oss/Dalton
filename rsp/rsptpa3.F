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
      SUBROUTINE CRTPA(VECA,VECB,VECC,VECD,VECBC,VECBD,VECCD,
     *                 RESVEC,CMO,UDV,PV,FOCK,FC,FV,
     *                 XINDX,MJWOP,WRK,LWRK)
C
#include <implicit.h>
#include <iratdef.h>
#include <dummy.h>
C
C PURPOSE:
C CALCULATION OF TWO-PHOTON ABSORPTION BETWEEN EXCITED STATES
C
      LOGICAL DOMOM, DIPLEN, ATEST
C
      CHARACTER*8 ALAB,BLAB,CLAB,DLAB
C
      DIMENSION VECA(*),VECB(*),VECC(*),VECD(*)
      DIMENSION VECBC(*),VECCD(*),VECBD(*)
      DIMENSION RESVEC(*)
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*)
      DIMENSION XINDX(*),WRK(*)
C
      PARAMETER ( D0 = 0.0D0 )
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
#include <inftpa.h>
C
      CALL QENTER('CRTPA')
      ATEST = .FALSE.
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
      WRITE(LUPRI,'(//A,A)')
     *' <<<<< CALCULATING CONTRIBUTIONS TO TWO-PHOTON',
     *' ABSORPTION AMPLITUDE >>>>>'
C
      DO 200 ISYMD = 1,NSYM
      DO 300 ISYMC = 1,NSYM
      DO 400 ISYMB = 1,NSYM
C
      ISYMDX = MULD2H(IREFSY,ISYMD)
      ISYMA = MULD2H(ISYMD,MULD2H(ISYMC,ISYMB))
      IF ( (NTPCN2(ISYMD).GT.0) .AND. (NTPCN1(ISYMC).GT.0) .AND.
     *     (NBTPOP(ISYMB).GT.0) .AND. (NATPOP(ISYMA).GT.0) ) THEN
C
      DO 500 ID   = 1,NTPCN2(ISYMD)
C
      DO 600 IC   = 1,NTPCN1(ISYMC)
C
      DO 700 IBOP = 1,NBTPOP(ISYMB)
      DO 750 IBFR = 1,NBTPFR
C
      DO 800 IAOP = 1,NATPOP(ISYMA)
C
C     Initialize variables.
C     Check if an equivalent moment calculation already has been done,
C     DOMOM indicates the result.
C     Read response vectors and eigen vectors from disk.
C     SWAP VECC
C     Check if some of the response vectors are equal or zero,
C     IBCDEQ indicates the result
C
C If excited state alpha (TPALP) is specified we only compute certain
C components.
C
      IF (TPALP .AND. ISYMD.NE.ISYMC .OR. ID.NE.IC ) GO TO 800
C
      CALL BCDCHK(DOMOM,IBCDEQ,LURSPRES,DIPLEN,DUMMY,
     *            ISYMA,ISYMB,ISYMC,ISYMD,ISYMBC,ISYMBD,ISYMCD,
     *            ALAB,BLAB,CLAB,DLAB,IAOP,IBOP,0,0,
     *            IBFR,IC,ID,FREQA,FREQB,FREQC,FREQD,
     *            KZYVA,KZYVB,KZYVC,KZYVD,KZYVBC,KZYVBD,KZYVCD,
     *            VECA,VECB,VECC,VECD,VECBC,VECBD,VECCD,IDUM)
C
      IF (.NOT.DOMOM) GOTO 800
      IF (TPALP .AND. ALAB.NE.BLAB ) GO TO 800
C
C    Initialize two-photon absorption amplitude
C
      TPARES = 0
C
      IF (IPRRSP.GT.0) WRITE(LUPRI,'(/A15,2A20,/A)')
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
     *           -FREQB,-FREQC,-FREQD,XINDX,UDV,PV,MJWOP,
     &           WRK,LWRK,CMO,FC)
      VAL = -DDOT(KZYVA,WRK,1,VECA,1)
      TPARES = TPARES + VAL
C
      IF (IPRRSP.GT.0)
     * WRITE(LUPRI,'(A17,F18.8,F20.8)')' Na T[4] Nb Nc Nd',VAL,TPARES
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
      VAL = DDOT(KZYVA,WRK,1,VECA,1)
      TMPVAL = VAL
      TPARES = TPARES + VAL
C
      IF (IBCDEQ.EQ.2) THEN
         TMPVAL = TMPVAL + VAL
         TPARES = TPARES + VAL
         CALL T3DRV(1,ISYMA,ISYMD,ISYMBC,VECD,VECBC,ATEST,VECA,
     *              -FREQD,-FREQB-FREQC,XINDX,UDV,PV,MJWOP,
     *              WRK,LWRK,CMO,FC,FV)
         VAL = DDOT(KZYVA,WRK,1,VECA,1)
         TMPVAL = TMPVAL + VAL
         TPARES = TPARES + VAL
      ELSE
         CALL T3DRV(1,ISYMA,ISYMC,ISYMBD,VECC,VECBD,ATEST,VECA,
     *              -FREQC,-FREQB-FREQD,XINDX,UDV,PV,MJWOP,
     *              WRK,LWRK,CMO,FC,FV)
         VAL = DDOT(KZYVA,WRK,1,VECA,1)
         TMPVAL = TMPVAL + VAL
         TPARES = TPARES + VAL
         CALL T3DRV(1,ISYMA,ISYMD,ISYMBC,VECD,VECBC,ATEST,VECA,
     *              -FREQD,-FREQB-FREQC,XINDX,UDV,PV,MJWOP,
     *              WRK,LWRK,CMO,FC,FV)
         VAL = DDOT(KZYVA,WRK,1,VECA,1)
         TMPVAL = TMPVAL + VAL
         TPARES = TPARES + VAL
      END IF
C
      IF (IPRRSP.GT.0)
     * WRITE(LUPRI,'(A15,2F20.8)')' Na T[3] Nx Nyz',TMPVAL,TPARES
C
      IF (IPRRSP.GT.5) CALL TIMER('T3DRV ',TIMSTR,TIMEND)

C
C     Calculate Na B[3] Nc Nd type terms (two permutations)
C
C
      IF (IPRRSP.GT.5) CALL TIMER('START ',TIMSTR,TIMEND)
C
      CALL X3INIT(KZYVA,KZYVC,KZYVD,ISYMA,ISYMC,ISYMD,BLAB,
     *            ISYMB,VECC,VECD,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVA,RESVEC,1,VECA,1)
      TMPVAL = VAL
      TPARES = TPARES + VAL
C
      IF (IPRRSP.GT.0)
     *  WRITE(LUPRI,'(A15,3F20.8)')' Na X[3] Ny Nz ',TMPVAL,TPARES
C
C
C     Calculate Nb A[3] Nc Nd type terms 
C     (two of six permutations in each call)
C
C
      CALL A3INIT(KZYVB,KZYVC,KZYVD,ISYMB,ISYMC,ISYMD,ALAB,
     *            ISYMA,VECC,VECD,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVB,RESVEC,1,VECB,1)
      TMPVAL = VAL
      TPARES = TPARES + VAL
C
      CALL A3INIT(KZYVC,KZYVB,KZYVD,ISYMC,ISYMB,ISYMD,ALAB,
     *            ISYMA,VECB,VECD,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVC,RESVEC,1,VECC,1)
      TMPVAL = TMPVAL + VAL
      TPARES = TPARES + VAL
C
      CALL A3INIT(KZYVD,KZYVB,KZYVC,ISYMD,ISYMB,ISYMC,ALAB,
     *            ISYMA,VECB,VECC,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVD,RESVEC,1,VECD,1)
      TMPVAL = TMPVAL + VAL
      TPARES = TPARES + VAL
C
      IF (IPRRSP.GT.0)
     *  WRITE(LUPRI,'(A15,2F20.8)')' Nx A[3] Ny Nz ',TMPVAL,TPARES
C
C
C     Calculate Na B[2] Ncd type terms (one permutations)
C
C
      CALL X2INIT(1,KZYVA,KZYVCD,ISYMA,ISPINA,ISYMCD,0,WRK(1),VECCD,
     *            RESVEC,XINDX,UDV,PV,BLAB,ISYMB,0,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVA,RESVEC,1,VECA,1)
      TMPVAL = VAL
      TPARES = TPARES + VAL
C
      IF (IPRRSP.GT.0)
     * WRITE(LUPRI,'(A15,2F20.8)') ' Na X[2] Nyz   ',TMPVAL,TPARES
C
C
C     Calculate Nb A[2] Ncd type terms (six permutations)
C
C
      CALL A2INIT(1,KZYVB,KZYVCD,ISYMB,ISPINB,ISYMCD,0,1,VECCD,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVB,RESVEC,1,VECB,1)
      TMPVAL = VAL
      TPARES = TPARES + VAL
C
      CALL A2INIT(1,KZYVCD,KZYVB,ISYMCD,0,ISYMB,ISPINB,1,VECB,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVCD,RESVEC,1,VECCD,1)
      TMPVAL = TMPVAL + VAL
      TPARES = TPARES + VAL
C
      CALL A2INIT(1,KZYVC,KZYVBD,ISYMC,ISPINC,ISYMBD,0,1,VECBD,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVC,RESVEC,1,VECC,1)
      TMPVAL = TMPVAL + VAL
      TPARES = TPARES + VAL
C
      CALL A2INIT(1,KZYVBD,KZYVC,ISYMBD,0,ISYMC,ISPINC,1,VECC,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVBD,RESVEC,1,VECBD,1)
      TMPVAL = TMPVAL + VAL
      TPARES = TPARES + VAL
C
      CALL A2INIT(1,KZYVD,KZYVBC,ISYMD,ISPIND,ISYMBC,0,1,VECBC,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVD,RESVEC,1,VECD,1)
      TMPVAL = TMPVAL + VAL
      TPARES = TPARES + VAL
C
      CALL A2INIT(1,KZYVBC,KZYVD,ISYMBC,0,ISYMD,ISPIND,1,VECD,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVBC,RESVEC,1,VECBC,1)
      TMPVAL = TMPVAL + VAL
      TPARES = TPARES + VAL
C
      IF (IPRRSP.GT.0)
     * WRITE(LUPRI,'(A15,2F20.8)') ' Nx A[2] Nyz   ',TMPVAL,TPARES
C
      IF (IPRRSP.GT.5) CALL TIMER('OTHERS',TIMSTR,TIMEND)
C
      WRITE(LUPRI,'(/A,2(/A,A10,I4,F12.6))')
     *     '@ Third order transition moment in a.u. for',
     *     '@ A operator, symmetry, frequency: ',ALAB,ISYMA,FREQA,
     *     '@ B operator, symmetry, frequency: ',BLAB,ISYMB,-FREQB
      WRITE(LUPRI,'(2(/A,2I4,F12.6))')
     *     '@ State no., symmetry, excitation energy:',IC,ISYMC,-FREQC,
     *     '@ State no., symmetry, excitation energy:',ID,ISYMD,FREQD
      WRITE(LUPRI,'(/A,F20.8/)') '@ < e | AB | f >  = ', TPARES
C
      WRITE(LURSPRES,'(/A,2(/A,A10,I4,F10.6))')
     *     ' Third order transition moment in a.u. for',
     *     ' A operator, symmetry, frequency: ',ALAB,ISYMA,FREQA,
     *     ' B operator, symmetry, frequency: ',BLAB,ISYMB,-FREQB
      WRITE(LURSPRES,'(2(/A,2I4,F10.6))')
     *     ' State no., symmetry, excitation energy:',IC,ISYMC,-FREQC,
     *     ' State no., symmetry, excitation energy:',ID,ISYMD,FREQD
      WRITE(LURSPRES,'(/A,F20.8)') ' < e | AB | f >  = ', TPARES
C
 800  CONTINUE
 750  CONTINUE
 700  CONTINUE
 600  CONTINUE
 500  CONTINUE
C
      END IF
C
 400  CONTINUE
 300  CONTINUE
 200  CONTINUE
      CALL GPCLOSE(LURSPRES,'KEEP')
C
C    End of subroutine CRTPA
C
      CALL QEXIT('CRTPA')
      RETURN
      END
