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
*======================================================================
      SUBROUTINE CC_PQI(CTR2,ISYCTR,T2AM,ISYTAMP,
     &                  FILNAM,LUPQIM,IADRPQ,IADR,IV,WORK,LWORK)
*----------------------------------------------------------------------
*
*     Purpose: Calculate the P and Q intermediates from the 
*              Lagrangian multipiers CTR2
*              and the amplitude vector T2AM
*              and write them to a file direct access file FILNAM
*              
*   P{aik,del} = sum_{dl} (2 t(dl,k;del) - t(dk,l;del)) Zeta(dl;ai)
*   P{aik,del} = sum_{dl} t(dl,k;del) (2Zeta(di;al) + Zeta(dl;ai))
*
*   Christof Haettig & Asger Halkier, August 1998
*
*======================================================================
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <iratdef.h>
#include <second.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER LWORK, LUPQIM
      CHARACTER*(*) FILNAM

#if defined (SYS_CRAY)
      REAL ONE, TWO, THREE, HALF, XNORM, DDOT, ZERO
      REAL WORK(*), CTR2(*), T2AM(*)
      REAL TIMET, TIMIO, TIMZWV, TIMTAM, DTIME, TIMTI
      REAL TIMZET,TIMSCL
#else
      DOUBLE PRECISION ONE, TWO, THREE, HALF, DDOT, ZERO
      DOUBLE PRECISION WORK(*), CTR2(*), T2AM(*)
      DOUBLE PRECISION TIMET,TIMIO,TIMTI,TIMZWV,TIMTAM,DTIME
      DOUBLE PRECISION TIMZET,TIMSCL
#endif
      PARAMETER(ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,HALF=0.5D0,ZERO=0.0D0)

      INTEGER IV, IADR, IADRQ
      INTEGER IADRPQ(MXCORB,IV)
      INTEGER ISYTAMP, ISYCTR, ISYMD, ISYTIN, ISYAIK
      INTEGER KEND1, LWRK1, KEND2, LWRK2
      INTEGER IOPT, IDEL, ILLL
      INTEGER KTINT, KTJNT, KPINT, KQINT, KLAMDH, KLAMDP, KT1AM
     

      CALL QENTER('CC_PQI')
*----------------------------------------------------------------------
* set symmetries and allocate work space:
*----------------------------------------------------------------------
      KLAMDH = 1
      KLAMDP = KLAMDH + NLAMDT
      KT1AM  = KLAMDP + NLAMDT
      KEND1  = KT1AM  + NT1AMX
      LWRK1  = LWORK - KEND1

      IF ( LWRK1 .LT. 0 ) THEN
         WRITE (LUPRI,*) 'Insufficient memory in CC_PQI.'
         CALL QUIT('Insufficient memory in CC_PQI.')
      END IF

      TIMET  = SECOND()
      TIMIO  = ZERO
      TIMTI  = ZERO
      TIMZWV = ZERO
      TIMTAM = ZERO
      TIMZET = ZERO
      TIMSCL = ZERO

      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'Norm of T2AMP:',
     &        DDOT(NT2AM(ISYTAMP),T2AM,1,T2AM,1)
         WRITE (LUPRI,*) 'T2AMP:'
         CALL CC_PRP(WORK,T2AM,ISYTAMP,0,1)
         WRITE (LUPRI,*) 'Norm of CTR2 :',
     &        DDOT(NT2SQ(ISYCTR),CTR2,1,CTR2,1)
         WRITE (LUPRI,*) 'CTR2:'
         CALL CC_PRSQ(WORK,CTR2,ISYCTR,0,1)
      END IF
*----------------------------------------------------------------------
* get XLAMD matrices from zero T1 amplitudes:
*----------------------------------------------------------------------
      CALL DZERO(WORK(KT1AM),NT1AMX)
       
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     &            WORK(KEND1),LWRK1)

*----------------------------------------------------------------------
* calculate the P intermediate in a loop over AO index delta:
*----------------------------------------------------------------------
      DO ISYMD = 1, NSYM
      DO ILLL = 1, NBAS(ISYMD)
        IDEL = IBAS(ISYMD) + ILLL
    
        ISYTIN = MULD2H(ISYMD,ISYTAMP)
        ISYAIK = MULD2H(ISYTIN,ISYCTR)
    
        KTINT = KEND1
        KTJNT = KTINT + NT2BCD(ISYTIN)
        KPINT = KTJNT + NT2BCD(ISYTIN)
        KQINT = KPINT + NT2BCD(ISYAIK)
        KEND2 = KQINT + NT2BCD(ISYAIK)
        LWRK2 = LWORK  - KEND2

        IF  (LWRK2 .LT. 0) THEN
           WRITE (LUPRI,*) 'Insufficient memory in CC_PQIM.'
           CALL QUIT('Insufficient memory in CC_PQIM.')
        END IF
C
C       calculate delta batch of backtransformed amplitudes:
C
        DTIME = SECOND()
        CALL CC_TI(WORK(KTINT),ISYTIN,T2AM,ISYTAMP,WORK(KLAMDH),1,
     &             WORK(KEND2),LWRK2,IDEL,ISYMD)
        TIMTI = TIMTI + SECOND() - DTIME
C
C       calculate 2 x Coul - Exc combination of backtransformed T:
C
        DTIME = SECOND()
        CALL DCOPY(NT2BCD(ISYTIN),WORK(KTINT),1,WORK(KTJNT),1)

        CALL DSCAL(NT2BCD(ISYTIN),-ONE,WORK(KTJNT),1)

        CALL CCLT_P21I(WORK(KTJNT),ISYTIN,WORK(KEND2),LWRK2,
     &                 IT2BCD,NT2BCD,IT1AM,NT1AM,NVIR)
        
        CALL DAXPY(NT2BCD(ISYTIN),TWO,WORK(KTINT),1,WORK(KTJNT),1)
        TIMTAM = TIMTAM + SECOND() - DTIME
C
C       calculate the P intermediate:
C
        IOPT = 1
        DTIME = SECOND()
        CALL CC_ZWVI(WORK(KPINT),CTR2,ISYCTR,WORK(KTJNT),
     &               ISYTIN,WORK(KEND2),LWRK2,IOPT)
        TIMZWV = TIMZWV + SECOND() - DTIME
        
        DTIME = SECOND()
        CALL DSCAL(NT2BCD(ISYAIK),HALF,WORK(KPINT),1)
        TIMSCL = TIMSCL + SECOND() - DTIME
C
C       write the intermediate to file:
C
        IADRPQ(IDEL,IV) = IADR

        DTIME = SECOND()
        CALL  PUTWA2(LUPQIM,FILNAM,WORK(KPINT),IADR,NT2BCD(ISYAIK))
        TIMIO = TIMIO + SECOND() - DTIME

        IADR = IADR + 2*NT2BCD(ISYAIK)

        IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'CC_PQI> P interm. in AO for IDEL = ',IDEL
           WRITE (LUPRI,'(5(F12.8))') 
     &          (WORK(KPINT+I),I=0,NT2BCD(ISYAIK)-1)
        END IF
        
      END DO
      END DO

*----------------------------------------------------------------------
*     calculate (2 x Exc + Coul)/3 combination of ZETA:
*     (we interrupt here the loop over AO to calculate the modified
*      ZETA and accept that we recalculate the backtransformed 
*      amplitude since it the transformation of the amplitudes is
*      less expansive than the transformation and restruction of ZETA 
*      for each delta. both the transformation of TAMP and the 
*      transformation/restruction of ZETA scale with N V^2 O^2. )
*----------------------------------------------------------------------
      DTIME = SECOND()
      CALL CCRHS_T2BT(CTR2,WORK(KEND1),LWRK1,ISYCTR)
      CALL CCSD_T2TP(CTR2,WORK(KEND1),LWRK1,ISYCTR)
      TIMZET = TIMZET + SECOND() - DTIME

*----------------------------------------------------------------------
* calculate the Q intermediate in a loop over AO index delta:
*----------------------------------------------------------------------
      DO ISYMD = 1, NSYM
      DO ILLL = 1, NBAS(ISYMD)
        IDEL = IBAS(ISYMD) + ILLL
    
        ISYTIN = MULD2H(ISYMD,ISYTAMP)
        ISYAIK = MULD2H(ISYTIN,ISYCTR)
    
        KTINT = KEND1
        KTJNT = KTINT + NT2BCD(ISYTIN)
        KPINT = KTJNT + NT2BCD(ISYTIN)
        KQINT = KPINT + NT2BCD(ISYAIK)
        KEND2 = KQINT + NT2BCD(ISYAIK)
        LWRK2 = LWORK  - KEND2

        IF  (LWRK2 .LT. 0) THEN
           WRITE (LUPRI,*) 'Insufficient memory in CC_PQIM.'
           CALL QUIT('Insufficient memory in CC_PQIM.')
        END IF
C
C       calculate delta batch of backtransformed amplitudes:
C
        DTIME = SECOND()
        CALL CC_TI(WORK(KTINT),ISYTIN,T2AM,ISYTAMP,WORK(KLAMDH),1,
     &             WORK(KEND2),LWRK2,IDEL,ISYMD)
        TIMTI = TIMTI + SECOND() - DTIME
C
C       calculate Q intermediate:
C
        IOPT = 2
        DTIME = SECOND()
        CALL CC_ZWVI(WORK(KQINT),CTR2,ISYCTR,WORK(KTINT),
     &               ISYTIN,WORK(KEND2),LWRK2,IOPT)
        TIMZWV = TIMZWV + SECOND() - DTIME

        DTIME = SECOND()
        CALL DSCAL(NT2BCD(ISYAIK),THREE*HALF,WORK(KQINT),1)
        TIMSCL = TIMSCL + SECOND() - DTIME
C
C       write the intermediate to file:
C
        IADRQ = IADRPQ(IDEL,IV) + NT2BCD(ISYAIK)

        DTIME = SECOND()
        CALL  PUTWA2(LUPQIM,FILNAM,WORK(KQINT),IADRQ,NT2BCD(ISYAIK))
        TIMIO = TIMIO + SECOND() - DTIME

        IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'CC_PQI> Q interm. in AO for IDEL = ',IDEL
           WRITE (LUPRI,'(5(F12.8))')
     &          (WORK(KQINT+I),I=0,NT2BCD(ISYAIK)-1)
        END IF
        
      END DO
      END DO

*----------------------------------------------------------------------
* recover Zeta vector:
*----------------------------------------------------------------------
      DTIME = SECOND()
      CALL CCSD_T2TP(CTR2,WORK(KEND1),LWRK1,ISYCTR)
      CALL CCRHS_T2TR(CTR2,WORK(KEND1),LWRK1,ISYCTR)
      TIMZET = TIMZET + SECOND() - DTIME

*----------------------------------------------------------------------
* return:
*----------------------------------------------------------------------
      IF (IPRINT.GE.10) THEN
         TIMET = SECOND() - TIMET
         WRITE(LUPRI,*) '  Timings of CC_PQIM:'
         WRITE(LUPRI,1) 'I/O             ', TIMIO
         WRITE(LUPRI,1) 'CC_TI           ', TIMTI
         WRITE(LUPRI,1) '2C-E of T2      ', TIMTAM
         WRITE(LUPRI,1) '2E+C of L2      ', TIMZET
         WRITE(LUPRI,1) 'CC_ZWV          ', TIMZWV
         WRITE(LUPRI,1) 'scaling etc.    ', TIMSCL
         WRITE(LUPRI,1) 'total CC_PQIM:  ', TIMET  
      END IF

   1  FORMAT(1x,'Time used for',2x,A18,2x,': ',f10.2,' seconds')
 
      CALL QEXIT('CC_PQI')
      RETURN
      END
*======================================================================
