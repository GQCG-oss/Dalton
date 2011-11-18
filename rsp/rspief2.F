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
C=============================================================================
C    /* Deck E3IEF */
C=============================================================================
      SUBROUTINE E3IEF2(VECA, VEC1, VEC2,ETRS,XINDX,ZYM1,ZYM2,
     *              DEN1,UDV,WORK,LWORK,KZYVR,KZYV1,KZYV2,
     *              IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP)
C
C
C     Purpose:
C     Outer driver routine for IEF-PCM solvent contribution 
C     to E[3] times two vectors. 
C     Completely rewritten from scratch!
C

C now we build the V[3] term of the gradient
#include <implicit.h>
#include <dummy.h>
C
#include <maxorb.h>
#include <mxcent.h>
#include <priunit.h>
#include <inforb.h>
#include <infdim.h>
#include <infinp.h>
#include <infvar.h>
#include <infrsp.h>
#include <infpri.h>
#include "rspprp.h"
#include <infcr.h>
#include <infspi.h>
#include <infden.h>
#include <pcmdef.h>
#include <pcm.h>
      PARAMETER ( D1 = 1.0D0, DM1 = -1.0D0, D2 = 2.0D0, DM2 = -2.0D0,
     $            DP5= 0.5D0, DMP5= -0.5D0)

      DIMENSION ETRS(KZYVR),XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI),DEN1(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),WORK(*),CMO(*)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2),VECA(KZYVR)
      DIMENSION MJWOP(2,MAXWOP,8)
      INTEGER ADDFLG

      LOGICAL DYNCHG,LORB,LCON,LREF

      WRITE(LUPRI,*) 'INITIAL ETRS'
      CALL OUTPUT(ETRS,1,KZYVR,1,1,KZYVR,1,1,LUPRI)

      KFREE = 1
      LFREE = LWORK
      CALL MEMGET('REAL',KCHG,  NTS, WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KQ12,  NTS, WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KQ1,   NTS, WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KQ2,   NTS, WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KQD1,  NTS, WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KQD2,  NTS, WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KQ1D2, NTS, WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KQ2D1, NTS, WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KQD12, NTS, WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCTS,3*NTS, WORK,KFREE,LFREE)
C
      CALL DZERO(WORK(KCHG),  NTS)
      CALL DZERO(WORK(KQ12),  NTS)
      CALL DZERO(WORK(KQ1),   NTS)
      CALL DZERO(WORK(KQ2),   NTS)
      CALL DZERO(WORK(KQD1),  NTS)
      CALL DZERO(WORK(KQD2),  NTS)
      CALL DZERO(WORK(KQ1D2), NTS)
      CALL DZERO(WORK(KQ2D1), NTS)
      CALL DZERO(WORK(KQD12), NTS)
      CALL DZERO(WORK(KCTS),3*NTS)

      DO I=1,NTS
         K = 3 * (I-1)
         WORK(KCTS + K)     = XTSCOR(I)
         WORK(KCTS + K + 1) = YTSCOR(I)
         WORK(KCTS + K + 2) = ZTSCOR(I)
      END DO

C we fetch the transformation vectors
      NSIM = 1
      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,ZYM1,MJWOP)
      CALL GTZYMT(NSIM,VEC2,KZYV2,ISYMV2,ZYM2,MJWOP)

c ISYMDN is the symmmetry of the density
      ISYMDN = 1
      OVLAP = D1
      DYNCHG = .FALSE.
      ISPING = ISPINA
      ISPIN1 = ISPINB
      ISPIN2 = ISPINC
      ISPIN3 = 0

      CALL DAXPY(NTS,D1,QSN,1,WORK(KCHG),1)
      CALL DAXPY(NTS,D1,QSE,1,WORK(KCHG),1)
C
C Expectation values of transformed charges
C
      IF(ISPIN1.EQ.0) THEN
         IKLVL = 1
         CALL RSPASC(IKLVL,NTS,ISPIN1,IDUMMY,IDUMMY,ISYMV1,IDUMMY,
     $        IDUMMY,ISYMDN,DYNCHG,OVLAP,ZYM1,DUMMY,DUMMY,WORK(KQ1),
     $        WORK(KCTS),UDV,CMO,WORK(KFREE),LFREE)
      END IF

      IF(ISPIN2.EQ.0) THEN
         IKLVL = 1
         CALL RSPASC(IKLVL,NTS,ISPIN2,IDUMMY,IDUMMY,ISYMV2,IDUMMY,
     $        IDUMMY,ISYMDN,DYNCHG,OVLAP,ZYM2,DUMMY,DUMMY,WORK(KQ2),
     $        WORK(KCTS),UDV,CMO,WORK(KFREE),LFREE)
      END IF

      IF(ISPIN1.EQ.ISPIN2) THEN
         IKLVL = 2
         CALL RSPASC(IKLVL,NTS,ISPIN1,ISPIN2,IDUMMY,ISYMV1,ISYMV2,
     $        IDUMMY,ISYMDN,DYNCHG,OVLAP,ZYM1,ZYM2,DUMMY,WORK(KQ12),
     $        WORK(KCTS),UDV,CMO,WORK(KFREE),LFREE)
         IKLVL = 2
         CALL RSPASC(IKLVL,NTS,ISPIN2,ISPIN1,IDUMMY,ISYMV2,ISYMV1,
     $        IDUMMY,ISYMDN,DYNCHG,OVLAP,ZYM2,ZYM1,DUMMY,WORK(KQ12),
     $        WORK(KCTS),UDV,CMO,WORK(KFREE),LFREE)
      END IF 

CLF      IF (MZCONF(ISYMV1).GT.0) THEN
      IF(.FALSE.) THEN
         NCASE = 1
         CALL DENGET(NCASE,DEN1)
         IKLVL = 0
         CALL RSPASC(IKLVL,NTS,IDUMMY,IDUMMY,IDUMMY,IDUMMY,IDUMMY,
     $        IDUMMY,ISYMDN,DYNCHG,OVLAP,DUMMY,DUMMY,DUMMY,WORK(KQD1),
     $        WORK(KCTS),DEN1,CMO,WORK(KFREE),LFREE)
      END IF

      IF(.FALSE.) THEN
         WRITE(LUPRI,*) 'TRANSFORMED CHARGES'
         WRITE(LUPRI,*) 'qsn+qse'
         CALL OUTPUT(WORK(KCHG),1,NTS,1,1,NTS,1,1,LUPRI)
         WRITE(LUPRI,*) 'q(1)'
         CALL OUTPUT(WORK(KQ1),1,NTS,1,1,NTS,1,1,LUPRI)
         WRITE(LUPRI,*) 'q(2)'
         CALL OUTPUT(WORK(KQ2),1,NTS,1,1,NTS,1,1,LUPRI)
         WRITE(LUPRI,*) 'q(12)'
         CALL OUTPUT(WORK(KQ12),1,NTS,1,1,NTS,1,1,LUPRI)
      END IF

      IKLVL = 0
      ADDFLG = 1
C
C Transformed potentials contracted with charges
C
      CALL RSPPOT(NSIM,KZYVR,IGRSYM,ISYMDN,IKLVL,NTS,ISPING,
     $     IDUMMY,IDUMMY,
     $     IDUMMY,IDUMMY,IDUMMY,IDUMMY,MJWOP,
     $     ADDFLG,OVLAP,DMP5,DUMMY,DUMMY,DUMMY,UDV,ETRS,WORK(KQ12),CMO,
     $     XINDX,WORK(KFREE),LFREE)

      IKLVL = 1
      ADDFLG = 1
C
C V(k2)*q(k1)
C
      CALL RSPPOT(NSIM,KZYVR,IGRSYM,ISYMDN,IKLVL,NTS,ISPING,
     $     ISPIN1,IDUMMY,
     $     IDUMMY,ISYMV1,ISYMV2,IDUMMY,MJWOP,
     $     ADDFLG,OVLAP,DM1,ZYM1,DUMMY,DUMMY,UDV,ETRS,WORK(KQ2),CMO,
     $     XINDX,WORK(KFREE),LFREE)
C
C V(k1)*q(k2)
C
      CALL RSPPOT(NSIM,KZYVR,IGRSYM,ISYMDN,IKLVL,NTS,ISPING,
     $     ISPIN2,IDUMMY,
     $     IDUMMY,ISYMV2,ISYMV1,IDUMMY,MJWOP,
     $     ADDFLG,OVLAP,DM1,ZYM2,DUMMY,DUMMY,UDV,ETRS,WORK(KQ1),CMO,
     $     XINDX,WORK(KFREE),LFREE)

      IKLVL = 2
      ADDFLG = 1
      CALL RSPPOT(NSIM,KZYVR,IGRSYM,ISYMDN,IKLVL,NTS,ISPING,
     $     ISPIN1,ISPIN2,
     $     IDUMMY,ISYMV1,ISYMV2,IDUMMY,MJWOP,
     $     ADDFLG,OVLAP,DMP5,ZYM1,ZYM2,DUMMY,UDV,ETRS,WORK(KCHG),CMO,
     $     XINDX,WORK(KFREE),LFREE)
      CALL RSPPOT(NSIM,KZYVR,IGRSYM,ISYMDN,IKLVL,NTS,ISPING,
     $     ISPIN2,ISPIN1,
     $     IDUMMY,ISYMV2,ISYMV1,IDUMMY,MJWOP,
     $     ADDFLG,OVLAP,DMP5,ZYM2,ZYM1,DUMMY,UDV,ETRS,WORK(KCHG),CMO,
     $     XINDX,WORK(KFREE),LFREE)

Clf      WRITE(LUPRI,*) 'FINAL ETRS'
Clf      CALL OUTPUT(ETRS,1,KZYVR,1,1,KZYVR,1,1,LUPRI)
      RETURN
      END

C=============================================================================
C    /* Deck RSPPOT */
C=============================================================================
      SUBROUTINE RSPPOT(NSIM,KZYVR,IGRSYM,ISYMDN,IKLVL,NTS,
     $     ISPING,ISPIN1,ISPIN2,ISPIN3,ISYMV1,ISYMV2,ISYMV3,
     $     MJWOP,ADDFLG,OVLAP,DFCTR,
     $     ZYM1,ZYM2,ZYM3,UDV,ETRS,CHG,CMO,XINDX,
     $     WORK,LWORK)
C
#include <implicit.h>
#include <inforb.h>
#include <infvar.h>
#include <inftap.h>
#include <priunit.h>
#include <infrsp.h>
#include <dummy.h>
#include <infspi.h>
C
c     Makes the potential part of the PCM solvent contribution to the
c     quadratic (and cubic) response gradient.
C
C
      INTEGER ADDFLG
      INTEGER MJWOP(2,MAXWOP,8)
      LOGICAL FNDLAB,LORB,LCON,LREF
      DIMENSION UDV(*),XINDX(*)
      DIMENSION CMO(*),CHG(*),WORK(*),ETRS(*)
      DIMENSION ZYM1(*),ZYM2(*),ZYM3(*)
C
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET('REAL',KVAOT, NNBASX,     WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVMOT, NNORBX,     WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVMO,  N2ORBX,     WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVMO2, N2ORBX,     WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVMO3, N2ORBX,     WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KUCMO, NORBT*NBAST,WORK,KFREE,LFREE)
C      CALL MEMGET('REAL',KV3Q0,N2ORBX,WORK,KFREE,LFREE)

      CALL DZERO(WORK(KVAOT),NNBASX)     
      CALL DZERO(WORK(KVMOT),NNORBX)     
      CALL DZERO(WORK(KVMO ),N2ORBX)     
      CALL DZERO(WORK(KVMO2),N2ORBX)     
      CALL DZERO(WORK(KVMO3),N2ORBX)     
      CALL DZERO(WORK(KUCMO),NORBT*NBAST)

      CALL UPKCMO(CMO,WORK(KUCMO))

#ifdef PCM_DEBUG
      print *,'beg of rsppot',iklvl,isymv1,isymv2,isymdn
#endif

Clf this initialization of symmetry needs to be rethought if one is to implement CR!!!!
      ISYM = 1
      IF (IKLVL.EQ.1) ISYM = ISYMV2
      IF (IKLVL.EQ.2) ISYM = 1
c      IF (IKLVL.EQ.3) ISYM = MULD2H(MULD2H(ISYMV1,ISYMV2),ISYMV3)
      ISYMT = MULD2H(ISYM,ISYMDN)
C      ISYMT = 1

c
C Calculate potentials in AO basis, mult. by the charges, transform to MO and unpack
C

Clf: NOTE: probably nosim needs to be changed according to symmetry.
      NOSIM = 1
      

#ifdef PCM_DEBUG
      write(lupri,*) 'chg before pot',iklvl
#endif
      CALL OUTPUT(CHG,1,NTS,1,1,NTS,1,1,LUPRI)
      CALL J1INT(CHG,.FALSE.,WORK(KVAOT),NOSIM,.FALSE.,'NPETES ',
     &     ISYMT,WORK(KFREE),LFREE)
      CALL UTHU(WORK(KVAOT),WORK(KVMOT),WORK(KUCMO),WORK(KFREE),
     $     NBAST,NORBT)
      CALL DSPTSI(NORBT,WORK(KVMOT),WORK(KVMO2))

#ifdef PCM_DEBUG
      write(lupri,*) 'potentialsXcharges', iklvl
#endif
      CALL OUTPUT(WORK(KVMO2),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)

      IGRSPI = 0
C     First transformation of potentials: V^e_{ab} --> V^e_{ab}({}^1\kappa) 
      IF (IKLVL.GE.1) THEN
         IGRSPI = ISPIN1
         CALL DZERO(WORK(KVMO),N2ORBX)
         CALL OITH1(ISYMV1,ZYM1,WORK(KVMO2),WORK(KVMO),ISYM)
         ISYM = ISYMV1
      END IF
C     Second transformation of potentials: V^e_{ab}({}^1\kappa) --> V^e_{ab}({}^1\kappa {}^2\kappa) 
      IF (IKLVL.GE.2) THEN
         IGRSPI = MULSP(IGRSPI,ISPIN2)
         CALL DZERO(WORK(KVMO2),N2ORBX)
         CALL OITH1(ISYMV2,ZYM2,WORK(KVMO),WORK(KVMO2),ISYM)
         ISYM = MULD2H(ISYM,ISYMV2)
      END IF 
C     Third transformation of potentials.....
      IF (IKLVL.GE.3) THEN
         IGRSPI = MULSP(IGRSPI,ISPIN3)
         CALL DZERO(WORK(KVMO),N2ORBX)
         CALL OITH1(ISYMV3,ZYM3,WORK(KVMO2),WORK(KVMO),ISYM)
         ISYM = MULD2H(ISYM,ISYMV3)
      END IF
      IF ((IKLVL.EQ.0).OR.(IKLVL.EQ.2)) THEN
         CALL DCOPY(N2ORBX,WORK(KVMO2),1,WORK(KVMO),1)
      END IF
#ifdef PCM_DEBUG
      write(lupri,*) 'potentialsXcharges after transformations', iklvl
#endif
      CALL OUTPUT(WORK(KVMO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)

C
C ADD A SCALING FACTOR AND SUM THE OPERATOR (ADDFLG.LE.0 ONLY WHEN DEBUGGING)
C
      CALL DZERO(WORK(KVMO2),N2ORBX)
      IF(ADDFLG.GT.0) THEN
         CALL DAXPY(N2ORBX,DFCTR,WORK(KVMO),1,WORK(KVMO2),1)
      END IF

      
      ISYMV = IREFSY
C      ISPIN = 0 (ISPING)
C   THIS MUST BE SET HERE!!!!!!!!!!!!!!!!
C      IGRSPI = 0
      NZYVEC = 0 
      NZCVEC = 0
      LORB = .TRUE.
      LCON = .FALSE.
      LREF = .TRUE.

      IF (.true.) THEN
         WRITE(LUPRI,*) 'Norm of TOPGET ',DNRM2(N2ORBX,work(kvmo2),1)
      END IF

      IF(ADDFLG.GT.0) THEN
         CALL PCM1GR(NSIM,KZYVR,IDUMMY,ISPING,IGRSYM,IGRSPI,ISYMV,ETRS,
     $        DUMMY,NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,WORK(KVMO2),XINDX,
     $        MJWOP,WORK(KFREE),LFREE,LORB,LCON,LREF)
Clf         WRITE(LUPRI,*) 'ETRS after the gradient'
Clf         CALL OUTPUT(ETRS,1,KZYVR,1,1,KZYVR,1,1,LUPRI)
      END IF

      RETURN
      END

C=============================================================================
C    /* Deck RSPASC */
C=============================================================================
      SUBROUTINE RSPASC(IKLVL,NTS,ISPIN1,ISPIN2,ISPIN3,ISYMV1,
     $     ISYMV2,ISYMV3,ISYMDN,DYNCHG,OVLAP,ZYM1,ZYM2,ZYM3,TCHG,
     $     CTS,DEN,CMO,WORK,LWORK)
C
#include <implicit.h>
#include <mxcent.h>
#include <inforb.h>
#include <inftap.h>
#include <infrsp.h>
#include <priunit.h>
#include <orgcom.h>
C
C     Makes the charge part of the PCM solvent contribution of the
C     quadraticn and (possibly) the cubic response. The subroutine
C     takes info about the n. of one-index transformations, spin and
C     symmetry. The required density is constructed inside.
C
C
C
      PARAMETER (D1=1.0D0)
      LOGICAL FNDLAB
      LOGICAL DYNCHG
      LOGICAL TRIMAT,EXP1VL,TOFILE
      DIMENSION CMO(*),DEN(*),TCHG(*),CTS(*),WORK(*)
      DIMENSION ZYM1(*),ZYM2(*),ZYM3(*)
      DIMENSION INTREP(9*MXCENT),INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)

      KFREE = 1
      LFREE = LWORK
Clf dirty fix!!
      NTSIRR = NTS
#ifdef PCM_DEBUG
      print *,'rspasc',lwork,lfree,nnbasx,norbt,nbast,n2orbx,nts,ntsirr
#endif
      CALL MEMGET('REAL',KCHGAO,NNBASX,     WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KUCMO ,NORBT*NBAST,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTELM ,N2ORBX,     WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTLMA ,N2ORBX,     WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTCHG ,NTS,        WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KQCHG ,NTS,        WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KQMAT ,NTS*NTSIRR, WORK,KFREE,LFREE)

      CALL DZERO(WORK(KCHGAO),NNBASX)
      CALL DZERO(WORK(KUCMO) ,NORBT*NBAST)
      CALL DZERO(WORK(KTELM) ,N2ORBX)
      CALL DZERO(WORK(KTLMA) ,N2ORBX)
      CALL DZERO(WORK(KTCHG) ,NTS) 
      CALL DZERO(WORK(KQCHG) ,NTS) 
      CALL DZERO(WORK(KQMAT) ,NTS*NTSIRR) 

C      CALL GPOPEN(LUPROP,'AOPROPER','OLD',' ',
C     *        'UNFORMATTED',IDUMMY,.FALSE.)

C
C     Unpack symmetry blocked CMO
C
      CALL UPKCMO(CMO,WORK(KUCMO))


      IF (IKLVL.LE.0) ISYM = 1
      IF (IKLVL.EQ.1) ISYM = ISYMV1
      IF (IKLVL.EQ.2) ISYM = MULD2H(ISYMV1,ISYMV2)
      IF (IKLVL.EQ.3) ISYM = MULD2H(MULD2H(ISYMV1,ISYMV2),ISYMV3)
      ISYMT = MULD2H(ISYM,ISYMDN)

C      REWIND(LUPROP)
      DO ITS=1,NTS
         ISYM = ISYMT
         CALL DZERO(WORK(KCHGAO),NNBASX)
C         IF (DYNCHG) THEN
C            IF (FNDLAB('J3-PCMIN',LUPROP)) THEN
C               CALL READT(LUPROP,NNBASX,WORK(KCHGAO)) 
C            ELSE
C               WRITE (LUPRI,'(/A)') ' Integral label J3-PCMIN not'//
C     &              'found on file AOPROPER'
C               CALL QUIT('Integral label not found in POPGET')
C            END IF
C         ELSE
C            IF (FNDLAB('J2-PCMIN',LUPROP)) THEN
C               CALL READT(LUPROP,NNBASX,WORK(KCHGAO))
C            ELSE
C               WRITE (LUPRI,'(/A)') ' Integral label J2-PCMIN not'//
C     &              'found on file AOPROPER'
C               CALL QUIT('Integral label not found in POPGET')
C            END IF
C         END IF
         L = 1
         NCOMP = NSYM
         KTS = 3*(ITS-1)
         DIPORG(1) = CTS(KTS+1)
         DIPORG(2) = CTS(KTS+2)
         DIPORG(3) = CTS(KTS+3)
         EXP1VL    = .FALSE.
         TOFILE    = .FALSE.
         KPATOM    = 0
         TRIMAT    = .TRUE.
         CALL GET1IN(WORK(KCHGAO),'NPETES ',NCOMP,WORK(KFREE),LWORK,
     &        LABINT,INTREP,INTADR,L,TOFILE,KPATOM,TRIMAT,
     &        DUMMY,EXP1VL,DUMMY,IPRRSP)
         JCHGAO = KCHGAO
         DO ILOP = 1, NSYM
            ISYM = ILOP
            JTS = (ILOP - 1)*NTSIRR + ITS
            CALL UTHU(WORK(JCHGAO),WORK(KTLMA),WORK(KUCMO),WORK(KFREE),
     $           NBAST,NORBT)
            CALL DSPTSI(NORBT,WORK(KTLMA),WORK(KTELM))
C     First transformation of charges: q^e_{ab} --> q^e_{ab}({}^1\kappa) 
            IF (IKLVL.GE.1) THEN
               CALL DZERO(WORK(KTLMA),N2ORBX)
               CALL OITH1(ISYMV1,ZYM1,WORK(KTELM),WORK(KTLMA),ISYM)
               ISYM = MULD2H(ISYM,ISYMV1)
            END IF
C     Second transformation of charges: q^e_{ab}({}^1\kappa) --> q^e_{ab}({}^1\kappa {}^2\kappa) 
            IF (IKLVL.GE.2) THEN
               CALL DZERO(WORK(KTELM),N2ORBX)
               CALL OITH1(ISYMV2,ZYM2,WORK(KTLMA),WORK(KTELM),ISYM)
               ISYM = MULD2H(ISYM,ISYMV2)
            END IF 
C     Third transformation of charges: hope you can figure out the formula.....
            IF (IKLVL.GE.3) THEN
               CALL DZERO(WORK(KTLMA),N2ORBX)
               CALL OITH1(ISYMV3,ZYM3,WORK(KTELM),WORK(KTLMA),ISYM)
               ISYM = MULD2H(ISYM,ISYMV3)
            END IF
            IF ((IKLVL.EQ.1).OR.(IKLVL.EQ.3)) THEN
               CALL DCOPY(N2ORBX,WORK(KTLMA),1,WORK(KTELM),1)
            END IF
C     Contract transformed charges with the density
            CALL MELONE(WORK(KTELM),ISYM,DEN,OVLAP,FACT,200,'RSPASC')
#ifdef PCM_DEBUG
            print *,its,fact,ilop,jts
#endif
            WORK(KTCHG + JTS - 1) = FACT
            JCHGAO = JCHGAO + NNBASX
         END DO
      END DO
      
      CALL V2Q(WORK(KQMAT),WORK(KTCHG),WORK(KQCHG),QTEXS,.false.)
Clf      CALL GPCLOSE(LUPCMD,'KEEP')
      CALL DAXPY(NTS,D1,WORK(KQCHG),1,TCHG,1)

      RETURN
      END

C=============================================================================
C    /* Deck PCMGR1 */
C=============================================================================
      SUBROUTINE PCM1GR(NSIM,KZYVR,INTSYM,ISPIN,IGRSYM,IGRSPI,
     *                  ISYMV,OTRS,
     *                  VEC,NZYVEC,NZCVEC,OVLAP,ISYMDN,DEN1,OPMAT,
     *                  XINDX,MJWOP,WORK,LWORK,LORB,LCON,LREFST)
C
C     Compute the gradient resulting from multiplying the S matrix
C     with one ore more vectors. Copied from RSP1GR
C
#include <implicit.h>
#include <infdim.h>
#include <inforb.h>
#include <priunit.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include "rspprp.h"
#include <infhyp.h>
#include <infvar.h>
#include <qrinf.h>
#include <infpri.h>
C
      LOGICAL LORB, LCON, LREFST
C
      DIMENSION OTRS (KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION VEC(NZYVEC)
      DIMENSION DEN1(NASHDI,NASHDI)
      DIMENSION OPMAT(NORBT,NORBT)
      DIMENSION XINDX(*)
      DIMENSION WORK(*)
C
#ifdef PCM_DEBUG
      print *,'lorb,lcon,lrefst',lorb,lcon,lrefst
#endif
clf      IF ( IPRRSP .GT. 150 ) THEN
      IF ( .true. ) THEN
         WRITE(LUPRI,'(A)') ' Vector in PCM1GR'
         IF (LREFST) WRITE(LUPRI,'(A)') ' (Reference state)'
         CALL OUTPUT(VEC,1,NZYVEC,1,NSIM,NZYVEC,NSIM,1,LUPRI)
         IF ( LORB ) THEN
            WRITE(LUPRI,'(//A)') ' Density matrix in PCP1GR'
            CALL OUTPUT(DEN1,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         END IF
         WRITE(LUPRI,'(//A)') ' One electron matrix in PCM1GR'
         CALL OUTPUT(OPMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      IF ( LORB ) THEN
         TRPLET = IGRSPI.NE.ISPIN
         CALL ORBSX(NSIM,IGRSYM,KZYVR,OTRS,OPMAT,OVLAP,
     *              ISYMDN,DEN1,MJWOP,WORK,LWORK)
      END IF
C
      IF ( LCON ) THEN
         ISYMJ  = MULD2H( IGRSYM, IREFSY )
         NZCSTJ = MZCONF( IGRSYM )
C
         TRPLET = ISPIN.EQ.1
         CALL CONSX(NSIM,KZYVR,IGRSYM,OPMAT,VEC,NZYVEC,
     *              NZCVEC,ISYMV,NZCSTJ,ISYMJ,LREFST,OTRS,XINDX,
     *              WORK,LWORK)
      END IF
C
      RETURN
      END

C=============================================================================
C    /* Deck DENGET */
C=============================================================================
      SUBROUTINE DENGET(NCASE,DEN1)
#include <implicit.h>
#include <dummy.h>
      LOGICAL LREF
C
C     L. Frediani, Nov 2005. Purpose: initialization of parameters for
C     different densities needed for singlet and triplet HF and MCSCF
C     quadratic response calculations
C

C
C parameters initialization
C
      GO TO (101,102,103,104,105,106,107,108,109,110), NCASE

 101  CONTINUE
c      ILSYM  = MULD2H(IREFSY,ISYMV1)
c      IRSYM  = MULD2H(IREFSY,ISYMV2)
c      NCL    = MZCONF(ISYMV1)
c      NCR    = MZCONF(ISYMV2)
c      KZVARL = MZYVAR(ISYMV1)
c      KZVARR = MZYVAR(ISYMV2)
c      LREF   = .FALSE.
cc      ISYMDN = MULD2H(ILSYM,IRSYM)
cC
cC     We get triplet density in case op. A is triplet
cC
c      JSPIN1 = MULSP(ISPIN1,ISPIN2)
c      JSPIN2 = 0

      GOTO 200
 102  CONTINUE
      GOTO 200
 103  CONTINUE
      GOTO 200
 104  CONTINUE
      GOTO 200
 105  CONTINUE
      GOTO 200
 106  CONTINUE
      GOTO 200
 107  CONTINUE
      GOTO 200
 108  CONTINUE
      GOTO 200
 109  CONTINUE
      GOTO 200
 110  CONTINUE
      GOTO 200
C
C Density calculation
C
 200  CONTINUE

c      CALL DZERO(DEN1,NASHT*NASHT)
c      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
c     *         VEC1,VEC2,OVLAP,DEN1,DUMMY,JSPIN1,JSPIN2,TDM,NORHO2,
c     *         XINDX,WRK,KFREE,LFREE,LREF)
      RETURN
      END
