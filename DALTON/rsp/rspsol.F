!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
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
C
#ifdef UNDEF
===========================================================================
/* Comdeck log */
941116-hjaaj: removed tab characters and other formatting
941109-pj+kvm: corrected code for triplet
940708-hjaaj
RSPSLV: stop if soppa
940423-hjaaj:
SLVINP: removed extra blank in '.IPRSLV '
---
RSPSOL.u module written by Poul Joergensen and Kurt Mikkelsen, summer 1993
===========================================================================
#endif
C  /* Deck slvinp */
      SUBROUTINE SLVINP(WORD)
C
#include "implicit.h"
C
#include "priunit.h"
#include "infslv.h"
#include "rspprp.h"
#include "infpri.h"
#include "infrsp.h"
C
      LOGICAL NEWDEF
      PARAMETER ( NTABLE = 2 )
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER*8 LABEL
C
      DATA TABLE /'.SLVANT', '.IPRSLV'/
C
C READ IN  INPUT
C
      NEWDEF = (WORD .EQ. '*RSPSLV ')
      ICHANG = 0
      IF (NEWDEF) THEN
      WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') GO TO 100
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2),I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' KEYWORD "',WORD,
     *            '" NOT RECOGNIZED IN SLVINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT(' ILLEGAL KEYWORD IN LRINP ')
    1          CONTINUE
                  SLVANT = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRSLV
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' PROMPT "',WORD,
     *            '" NOT RECOGNIZED IN SLVINP.'
               CALL QUIT(' ILLEGAL PROMPT IN SLVINP ')
            END IF
         GO TO 100
      END IF
  300 CONTINUE
      WRITE(LUPRI,'(/A)')' ********* SLVINP ********'
      NAQRTO = 0
      NBQRTO = 0
      NCQRTO = 0
      IF (ICHANG .GT. 0) THEN
         WRITE(LUPRI,*)' print level in solvent, IPRSLV =', IPRSLV
         WRITE(LUPRI,*)' test with SLVANT =', SLVANT
      END IF
C
C *** END OF SLVINP
C
      RETURN
      END
C  /* Deck rspslv */
      SUBROUTINE RSPSLV(NCSIM,NOSIM,BCVECS,BOVECS,CMO,INDXCI,
     &                  UDV,EVECS,WRK,LWRK)
C
#if defined (HAS_PCMSOLVER)
      use pcm_config, only: pcm_configuration, pcm_cfg
      use pcm_linrsp, only: pcm_lr_driver
#endif
      use qmcmm_response, only: qmcmm_lr
      use pelib_interface, only: use_pelib, pelib_ifc_lr
C
#include "implicit.h"
#include "maxorb.h"
C
#include "priunit.h"
#include "infpri.h"
#include "infinp.h"
#include "pcmlog.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "inforb.h"
#include "infslv.h"
C---------------------
#include "gnrinf.h"
#include "infpar.h"
#include "mxcent.h"
#include "pcmdef.h"
#include "pcm.h"
#include "qm3.h"
#include "qmmm.h"
C---------------------
C
      PARAMETER ( DM1 = -1.0D0 )
C
      DIMENSION BCVECS(*),BOVECS(*),CMO(*),INDXCI(*)
      DIMENSION UDV(NASHT,*),EVECS(*),WRK(LWRK)
C
C
C
      CALL QENTER('RSPSLV')

C      IF (SOPPA) THEN
C         WRITE (LUPRI,*)
C     &   'RSPSLV-ERROR: solvent response is not implemented for SOPPA'
C         CALL QUIT(
C     &   'RSPSLV-ERROR: solvent response is not implemented for SOPPA')
C      END IF
Clf      IF (IREFSY.EQ.KSYMST .AND. NASHT.GT.1 .AND. KZCONF.NE.NCREF)
Clf     &   THEN
ClfC
ClfChj-May2000: it is not so quick to implement triplet solvent response
ClfC  with CSF's because RSPSLV calls sirius routines, and all the Sirius
ClfC  routines are using CSF's if that was used in Sirius.
ClfC  However, when triplet we need to use determinants to get the triplet
ClfC  configurations included.
ClfC
Clf         WRITE (LUPRI,*) 'RSPSLV-ERROR: '//
Clf     &      'solvent triplet response is not implemented for CSFs'
Clf         WRITE (LUPRI,*) 'KZCONF .ne. NCREF:',KZCONF,NCREF
Clf         WRITE(LUPRI,*)'Use .DETERMINANTS and try again!'
Clf         CALL QUIT('RSPSLV-ERROR: '//
Clf     &      'solvent triplet response is not implemented for CSFs')
Clf      END IF
C
      IPRSAV = IPRRSP
      IPRRSP = IPRSLV
      IF (SLVANT) THEN
         DO 501 I = 1,KZWOPT
            BOVECS(I+KZWOPT) =- BOVECS(I)
 501     CONTINUE
      END IF
      IF (IPRRSP. GT . 101 ) THEN
         IF (NCSIM.GT.0) THEN
            WRITE(LUPRI,*) NCSIM ,' CONFIGURATION TRIAL VECTORS'
            CALL OUTPUT(BCVECS,1,KZCONF,1,NCSIM,KZCONF,NCSIM,-1,LUPRI)
         END IF
         IF (NOSIM.GT.0) THEN
            WRITE(LUPRI,*) NOSIM ,' (Z_i Y_i) ORBITAL TRIAL VECTORS'
            CALL OUTPUT(BOVECS,1,KZWOPT,1,2*NOSIM,
     &           KZWOPT,2*NOSIM,-1,LUPRI)
         END IF
      END IF
C
C WORK ALLOCATION
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWRK
C
      CALL MEMGET('REAL',KCREF ,NCREF,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDVT  ,(N2ASHX*NCSIM),WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDVTTR,(N2ASHX*NCSIM),WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDV   ,NNASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDVTR ,NNASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KUDVTR,N2ASHX,WRK,KFREE,LFREE)

      CALL GETREF(WRK(KCREF),NCREF)
C
C CALCULATE TRANSITION DENSITY MATRIX
C
      IF ( (NCSIM.GT.0) .AND. (.NOT.RSPCI) ) THEN
         JSPIN1 = 0
         JSPIN2 = 0
         CALL RSPTDM(NCSIM,IREFSY,KSYMST,NCREF,KZCONF,WRK(KCREF),
     *                 BCVECS,WRK(KDVT),DUMMY,
     *                 JSPIN1,JSPIN2,.TRUE.,.TRUE.,
     *                 INDXCI,WRK,KFREE,LFREE)
C        CALL RSPTDM(NCSIM,ILRESY,IRSYM,NCLREF,NCRDIM,CLREF,
C    *                 CR, RHO1,RHO2, ISPIN1,ISPIN2,TDM,NORHO2,
C    *                 XNDXCI,WORK,KFREE,LFREE)
C |OR> has minus sign, scale density with minus
         CALL DSCAL(NCSIM*N2ASHX,DM1,WRK(KDVT),1)
         IF (TRPLET) THEN
            JSPIN1 = 1
            JSPIN2 = 0
         CALL RSPTDM(NCSIM,IREFSY,KSYMST,NCREF,KZCONF,WRK(KCREF),
     *                 BCVECS,WRK(KDVTTR),DUMMY,
     *                 JSPIN1,JSPIN2,.TRUE.,.TRUE.,
     *                 INDXCI,WRK,KFREE,LFREE)
C        CALL RSPTDM(NCSIM,ILRESY,IRSYM,NCLREF,NCRDIM,CLREF,
C    *                 CR, RHO1,RHO2, ISPIN1,ISPIN2,TDM,NORHO2,
C    *                 XNDXCI,WORK,KFREE,LFREE)
C |OR> has minus sign, scale density with minus
            CALL DSCAL(NCSIM*N2ASHX,DM1,WRK(KDVTTR),1)
         END IF
      END IF
C
      IF (TRPLET) THEN
C
C CONSTRUCT TRIPLET ONE ELECTRON DENSITY MATRIX
C
         JSPIN1 = 1
         JSPIN2 = 0
         CALL RSPDM(IREFSY,IREFSY,NCREF,NCREF,WRK(KCREF),
     *                 WRK(KCREF),WRK(KUDVTR),DUMMY,
     *                 JSPIN1,JSPIN2,.FALSE.,.TRUE.,
     *                 INDXCI,WRK,KFREE,LFREE)
C
C TRIANGULAR PACKING OF TRIPLET ONE-ELECTRON DENSITY MATRIX
C
         IJ = 0
         DO 103 I = 1,NASHT
            DO 203 J = 1,I
               IJ = IJ + 1
               WRK(KDVTR-1+IJ) = WRK(KUDVTR + (I-1)*NASHT + J -1 )
 203        CONTINUE
 103     CONTINUE
         IF (IPRRSP.GT.50 .AND. NASHT.GT.0) THEN
            WRITE(LUPRI,*)' TRIPLET ONE-ELECTRON DENSITY MATRIX'
            CALL OUTPUT(WRK(KUDVTR),1,NASHT,1,NASHT,NASHT,NASHT,
     *                                            -1,LUPRI)
         END IF
      END IF
C
C TRIANGULAR PACKING OF ONE ELECTRON DENSITY MATRIX
C
      IJ = 0
      DO 100 I = 1,NASHT
         DO 200 J = 1,I
            IJ = IJ + 1
            WRK(KDV-1+IJ) = UDV(I,J)
 200     CONTINUE
 100  CONTINUE
      IF (SLVANT) THEN
         KFREE1 = KFREE + NNASHX*NCSIM
         LFREE1 = LWRK - KFREE1
         DO 300 ICSIM = 1,NCSIM
            KOFDVT = (ICSIM-1)*N2ASHX
            KOFD2T = (ICSIM-1)*NNASHX
            IJ = 0
            DO 400 I = 1,NASHT
               DO 500 J = 1,I
                  IJ = IJ + 1
                  WRK(KFREE+KOFD2T+IJ-1) =
     *                        WRK(KDVT+KOFDVT+(I-1)*NASHT+J-1) +
     *                        WRK(KDVT+KOFDVT+(J-1)*NASHT+I-1)
 500           CONTINUE
 400        CONTINUE
 300     CONTINUE
         CALL DZERO(EVECS,(NOSIM+NCSIM)*KZYVAR)
         WRITE(LUPRI,*)' TEST CALCULATION WITH SLVANT =',SLVANT
         IF (IPRRSP .GT. 2) WRITE(LUPRI,*)
     *   ' CONTRIBUTIONS FROM SOLLIN TO LINEAR TRANSFORMATION'
         CALL SOLLIN(NCSIM,NOSIM,BCVECS,BOVECS,WRK(KCREF),CMO,INDXCI,
     &               WRK(KDV),WRK(KFREE),EVECS,EVECS(NCSIM*KZYVAR+1),
     &               .FALSE.,WRK(KFREE1),LFREE1)
         IF ((NCSIM.GT.0).AND.(IPRRSP.GT.101)) THEN
            WRITE(LUPRI,*)' Linear transformed CONFIGURATION VECTORS'
            CALL OUTPUT(EVECS,1,KZVAR,1,2*NCSIM,KZVAR,2*NCSIM,-1,LUPRI)
         END IF
         IF ((NOSIM.GT.0).AND.(IPRRSP.GT.101)) THEN
            WRITE(LUPRI,*)' Linear transformed ORBITAL VECTORS'
            CALL OUTPUT(EVECS(1+KZYVAR*NCSIM),1,KZVAR,1,2*NOSIM,
     *                                          KZVAR,2*NOSIM,-1,LUPRI)
         END IF
         CALL DZERO(EVECS,(NOSIM+NCSIM)*KZYVAR)
      END IF
      IF (QM3 .AND. (.NOT. MMPCM)) THEN
         IF (IPRRSP .GT. 2) WRITE(LUPRI,*)
     *   ' CONTRIBUTIONS FROM QM3LNO TO LINEAR TRANSFORMATION'
#if defined(VAR_MPI)
         IF (NODTOT .GE. 1) THEN
            CALL QM3LNO_P(NOSIM,BOVECS,WRK(KCREF),CMO,INDXCI,UDV,
     &               WRK(KDV),WRK(KUDVTR),EVECS(NCSIM*KZYVAR+1),
     &               WRK(KFREE),LFREE)
         ELSE
#endif
            CALL QM3LNO(NOSIM,BOVECS,WRK(KCREF),CMO,INDXCI,UDV,
     &               WRK(KDV),WRK(KUDVTR),WRK(KDVTR),
     &               EVECS(NCSIM*KZYVAR+1),
     &               WRK(KFREE),LFREE)
#if defined(VAR_MPI)
         END IF
#endif
      ENDIF
      IF (PCM .AND. (.NOT. MMPCM)) THEN
         LADDMM = .TRUE.
         IF (IPRRSP .GT. 2) WRITE(LUPRI,*)
     *   ' CONTRIBUTIONS FROM PCMLTR TO LINEAR TRANSFORMATION'
         CALL PCMLTR(NCSIM,NOSIM,BCVECS,BOVECS,WRK(KCREF),CMO,INDXCI,
     &               UDV,WRK(KDV),WRK(KUDVTR),WRK(KDVTR),
     &               WRK(KDVT),WRK(KDVTTR),
     &               EVECS,EVECS(NCSIM*KZYVAR+1),WRK(KFREE),LFREE)
      ENDIF

      IF (MMPCM) THEN
         IF (NOSIM .GT. MXEPCM) THEN
           WRITE(LUPRI,*) 'Error Increase MXEPCM to at least: ',NOSIM
           CALL QUIT('Increase MXEPCM')
         ENDIF
         LADDMM = .FALSE.
         DMMSAVE = 0.0D0
         FIRST1 = .TRUE.
         DO WHILE (.NOT. LADDMM)
            CALL QM3LNO(NOSIM,BOVECS,WRK(KCREF),CMO,INDXCI,UDV,
     &                 WRK(KDV),WRK(KUDVTR),WRK(KDVTR),
     &                 EVECS(NCSIM*KZYVAR+1),
     &                 WRK(KFREE),LFREE)
            IF (FIRST1) THEN
               CALL PCM_INIT_QVFLAGS('EL1')
               CALL PCM_COMP_POT(CMO, UDV, WRK(KUDVTR), BOVECS, DUMMY,
     &                          DUMMY, NOSIM, IDUMMY, WRK(KFREE), LFREE)
               CALL PCM_COMP_CHG(NOSIM, NEQRSP, WRK(KFREE), LFREE)
            END IF
            CALL PCM_INIT_QVFLAGS('MM1')
            CALL PCM_COMP_POT(DUMMY, DUMMY, DUMMY, DUMMY, DUMMY,
     &                       DUMMY,NOSIM,IDUMMY, WRK(KFREE), LFREE)
            CALL PCM_COMP_CHG(NOSIM, NEQRSP, WRK(KFREE), LFREE)
            CALL PCM_WRITE_PCM2MM(NOSIM,NTS,XTSCOR,YTSCOR,ZTSCOR,
     &                           QSE1,QSMM1)
            FIRST1 = .FALSE.
         ENDDO
         CALL PCMLTR(NCSIM,NOSIM,BCVECS,BOVECS,WRK(KCREF),CMO,
     &        INDXCI,
     &        UDV,WRK(KDV),WRK(KUDVTR),WRK(KDVTR),
     &        WRK(KDVT),WRK(KDVTTR),
     &        EVECS,EVECS(NCSIM*KZYVAR+1),WRK(KFREE),LFREE)
         FIRST1 = .FALSE.
      ENDIF
      IF (USE_PELIB()) THEN
        CALL PELIB_IFC_LR(NCSIM, NOSIM, BCVECS, BOVECS, WRK(KCREF), CMO,
     &                    INDXCI, UDV, WRK(KDV), WRK(KUDVTR),
     &                    WRK(KDVTR), WRK(KDVT), WRK(KDVTTR), EVECS,
     &                    EVECS(NCSIM*KZYVAR+1), WRK(KFREE), LFREE)
      END IF
#if defined (HAS_PCMSOLVER)
      if (pcm_cfg%do_pcm) then
         call pcm_lr_driver(nosim,bovecs,cmo,indxci,
     &               udv,wrk(kdv),wrk(kudvtr),wrk(kdvtr),
     &               wrk(kdvt),wrk(kdvttr),
     &               evecs(ncsim*kzyvar+1),wrk(kfree),lfree)
      end if
#endif
      IF (QMMM) THEN
         IF (IPQMMM .GT. 1) CALL TIMER('START ',TIMSTR,TIMEND)
         IF (IPRRSP .GT. 2) WRITE(LUPRI,*)
     *   ' CONTRIBUTIONS FROM QMMMLTR TO LINEAR TRANSFORMATION'
         CALL QMMMLTR(NCSIM,NOSIM,BCVECS,BOVECS,WRK(KCREF),CMO,INDXCI,
     &               UDV,WRK(KDV),WRK(KUDVTR),WRK(KDVTR),
     &               WRK(KDVT),WRK(KDVTTR),
     &               EVECS,EVECS(NCSIM*KZYVAR+1),WRK(KFREE),LFREE)
         IF (IPQMMM .GT. 1) CALL TIMER('QMMMLNO',TIMSTR,TIMEND)
      ENDIF
      IF (FLAG(16)) THEN
         IF (IPRRSP .GT. 2) WRITE(LUPRI,*)
     *   ' CONTRIBUTIONS FROM SLVLTR TO LINEAR TRANSFORMATION'
         CALL SLVLTR(NCSIM,NOSIM,BCVECS,BOVECS,WRK(KCREF),CMO,INDXCI,
     &               UDV,WRK(KDV),WRK(KUDVTR),WRK(KDVTR),
     &               WRK(KDVT),WRK(KDVTTR),
     &               EVECS,EVECS(NCSIM*KZYVAR+1),WRK(KFREE),LFREE)
      END IF
      IF (QMNPMM) THEN
         IF (IPRRSP .GT. 2) WRITE(LUPRI,*)
     *   ' CONTRIBUTIONS FROM QMCMM_LR TO LINEAR TRANSFORMATION'
         CALL qmcmm_lr(NOSIM,BOVECS,WRK(KCREF),CMO,
     &                 UDV,WRK(KDV),WRK(KUDVTR),WRK(KDVTR),
     &                 EVECS(NCSIM*KZYVAR+1),WRK(KFREE),LFREE)
      ENDIF
C
      IF ((NCSIM.GT.0).AND.(IPRRSP.GT.101)) THEN
         WRITE(LUPRI,*)' LINEAR TRANSFORMED CONFIGURATION VECTORS'
         CALL OUTPUT(EVECS,1,KZVAR,1,2*NCSIM,KZVAR,2*NCSIM,-1,LUPRI)
      END IF
      IF ((NOSIM.GT.0).AND.(IPRRSP.GT.101)) THEN
         WRITE(LUPRI,*)' LINEAR TRANSFIRMED ORBITAL VECTOR'
         CALL OUTPUT(EVECS(1+NCSIM*KZYVAR),1,KZVAR,1,2*NOSIM,
     *                                       KZVAR,2*NOSIM,-1,LUPRI)
      END IF
      IF (SLVANT) THEN
         WRITE(LUPRI,*)' TEST CALCULATION : SLVANT = ', SLVANT
         CALL QUIT(' RSPSLV; TEST calculation with SLVANT finished ')
      END IF
      CALL MEMREL('RSPSLV',WRK,1,KFRSAV,KFREE,LFREE)
      IPRRSP = IPRSAV
      CALL QEXIT('RSPSLV')
      RETURN
C END RSPSLV
      END
C  /* Deck slvltr */
      SUBROUTINE SLVLTR(NCSIM,NOSIM,BCVECS,BOVECS,CREF,CMO,INDXCI,
     &                  UDV,DV,UDVTR,DVTR,DTV,DTVTR,SCVECS,SOVECS,
     &                  WRK,LWRK)
C
C KM + PJ 6dec92 based on SOLLIN
C
C Common driver for SLVLNC and SLVLNO
C
#include "implicit.h"
#include "dummy.h"
      DIMENSION BCVECS(*),BOVECS(*),CREF(*),CMO(*),INDXCI(*)
      DIMENSION UDV(*),DV(*),DTV(*),SCVECS(*),SOVECS(*),WRK(LWRK)
      DIMENSION UDVTR(*),DVTR(*),DTVTR(*)
C
C Used from common blocks:
C   INFINP : LSOLMX,NLMSOL, EPSOL, INERSF
C
#include "maxorb.h"
#include "priunit.h"
#include "infpri.h"
#include "infinp.h"
#include "infrsp.h"
#include "inftap.h"
#include "wrkrsp.h"
C     logical inersf
      inersf = .false.

C
      CALL QENTER('SLVLTR')
      KFRSAV = 1
      LFRSAV = LWRK
      KFREE  = KFRSAV
      LFREE  = LWRK
C
C     Calculate f(l) factors;
C     Rsol and EPsol dependence is located here.
C
      CALL MEMGET('REAL',KFLVEC,NLMSOL,WRK,KFREE,LFREE)
      CALL SOLFL(WRK(KFLVEC),EPSOL,RSOL,LSOLMX)
      IF (INERSI) THEN
         CALL MEMGET('REAL',KFLINR,NLMSOL,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KTLMSI,NLMSOL,WRK,KFREE,LFREE)
         CALL SOLINR(WRK(KFLVEC),WRK(KFLINR),WRK(KTLMSI))
      ELSE
         CALL MEMGET('REAL',KFLINR,0,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KTLMSI,0,WRK,KFREE,LFREE)
      END IF
      CALL MEMGET('INTE',KSYRLM,(2*LSOLMX+1),WRK,KFREE,LFREE)
C
      IF (LUSOL .LE. 0)
     &  CALL GPOPEN(LUSOL,FNSOL,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
      IF (NCSIM .GT. 0) THEN
         IF (IPRRSP.GT.101) THEN
            WRITE(LUPRI,*)' LINEAR TRANSFORMED CONFIGURATION VECTOR'
            WRITE(LUPRI,*)' **** BEFORE SLVLNC **** '
            CALL OUTPUT(SCVECS,1,KZYVAR,1,NCSIM,KZYVAR,NCSIM,1,LUPRI)
         END IF
         CALL SLVLNC(NCSIM,BCVECS,CREF,CMO,INDXCI,
     *               WRK(KFLVEC),WRK(KFLINR),WRK(KTLMSI),
     *               UDV,DV,UDVTR,DVTR,DTV,DTVTR,SCVECS,WRK(KSYRLM),
     *               WRK(KFREE),LFREE)
         IF (IPRRSP.GT.101) THEN
            WRITE(LUPRI,*)' LINEAR TRANSFORMED CONFIGURATION VECTOR'
            WRITE(LUPRI,*)' **** AFTER SLVLNC **** '
            CALL OUTPUT(SCVECS,1,KZYVAR,1,NCSIM,KZYVAR,NCSIM,1,LUPRI)
         END IF
      END IF
      IF ( NOSIM .GT.0 ) THEN
         IF (IPRRSP.GT.101) THEN
            WRITE(LUPRI,*)' LINEAR TRANSFORMED ORBITAL VECTOR'
            WRITE(LUPRI,*)' **** AFTER SLVLNO **** '
            CALL OUTPUT(SOVECS,1,KZYVAR,1,NOSIM,KZYVAR,NOSIM,1,LUPRI)
         END IF
         CALL SLVLNO(NCSIM,NOSIM,BOVECS,CREF,CMO,INDXCI,
     *               WRK(KFLVEC),WRK(KFLINR),WRK(KTLMSI),
     *               UDV,DV,UDVTR,DVTR,SOVECS,WRK(KSYRLM),
     *               WRK(KFREE),LFREE)
         IF (IPRRSP.GT.101) THEN
            WRITE(LUPRI,*)' LINEAR TRANSFORMED ORBITAL VECTOR'
            WRITE(LUPRI,*)' **** AFTER SLVLNO **** '
            CALL OUTPUT(SOVECS,1,KZYVAR,1,NOSIM,KZYVAR,NOSIM,1,LUPRI)
         END IF
      END IF
      CALL GPCLOSE(LUSOL,'KEEP')
      CALL MEMREL('SLVLTR',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('SLVLTR')
      RETURN
      END
C  /* Deck slvlnc */
      SUBROUTINE SLVLNC(NCSIM,BCVEC,CREF,CMO,INDXCI,
     *                  FLVEC,FLINR,TLMSI,
     *                  UDV,DV,UDVTR,DVTR,DTV,DTVTR,SVEC,ISYRLM,
     *                  WRK,LFREE)
C
C 6dec92
C
C  Purpose:  Calculate MCSCF E2  contribution from a
C            surrounding medium to a csf trial vector.
C            Cavity radius = Rsol and dielectric constant = EPsol.
C
#include "implicit.h"
      DIMENSION BCVEC(*),  CREF(*), CMO(*),  ISYRLM(*)
      DIMENSION FLVEC(NLMSOL), FLINR(NLMSOL), TLMSI(NLMSOL)
      DIMENSION INDXCI(*), UDV(*), DV(*),   DTV(N2ASHX,*)
      DIMENSION SVEC(KZYVAR,*),     WRK(*)
      DIMENSION UDVTR(*),DVTR(*),DTVTR(N2ASHX,*)
C
#include "iratdef.h"
C
      PARAMETER ( D0 = 0.0D0 , D2 = 2.0D0 , THRZER = 1.0D-14 )
C
C
C  Used from common blocks:
C    INFINP : NLMSOL, LSOLMX, INERSI
C    INFORB : NNASHX, NNORBX, NNBASX, etc.
C    INFTAP : LUSOL, LBSYMB
C
#include "maxash.h"
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "inforb.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "inftap.h"
#include "infpri.h"

C
C     Statement functions
C
      TN(   LM) = WRK((KTNLM -1)+LM)
C
      CALL QENTER('SLVLNC')
C
C     Core allocation
C
      KRXC   = 1
      KRYC   = KRXC   + NCSIM*NNORBX
      KTXCAC = KRYC   +       NNORBX
      KW10   = KTXCAC + NCSIM
C
      KTNLM  = KW10
      KUCMO  = KTNLM  + NLMSOL
      KRLMAC = KUCMO  + NORBT*NBAST
      KRLM   = KRLMAC + NNASHX
      KW20   = KRLM   + NNORBX
C     2.1 read rlmao in ao basis and transform to rlm in mo basis
      KRLMAO = KW20
      KW21   = KRLMAO + NNBASX
C     3.0 SOLSC
      KRXCAC = KW10
      KRYCAC = KRXCAC + NCSIM*NNASHX
      KW30   = KRYCAC +       NNASHX
      LW30   = LFREE  + 1 - KW30
C
C     4.0 rspsor
      KURXC  = KRYC   + NNORBX
      KURYC  = KURXC  + N2ORBX
      KW40   = KURYC  + N2ORBX
      KNEED  = MAX(KW21,KW30)
      KNEED  = MAX(KNEED,KW40)
      IF (KNEED .GT. LFREE) CALL ERRWRK('SLVLNC',KNEED,LFREE)
C
C
C     Read nuclear contributions to ERLM (TN(l,m))
C     and position for R(l,m).
C
      REWIND LUSOL
      CALL MOLLAB('SOLVRLM ',LUSOL,LUERR)
      READ (LUSOL)
      CALL READT(LUSOL,NLMSOL,WRK(KTNLM))
C
C
C     Construct effective operators Rxc and Ryc
C     =========================================
C
C     Zero Rxc and Ryc storage
C
      CALL DZERO(WRK(KRXC), NCSIM * NNORBX )
      CALL DZERO(WRK(KRYC),         NNORBX )
      TYCAC = D0
C
C     Unpack symmetry blocked CMO
C
      CALL UPKCMO(CMO,WRK(KUCMO))
C
C     Loop over l,m expansion.
C
      LM = 0
      DO 520 L = 0,LSOLMX
         READ (LUSOL) L1,(ISYRLM(M),M=1,2*L+1)
      DO 500 M = -L,L
         LM = LM + 1
         IF (ISYRLM(L+M+1) .NE. KSYMOP .AND. ISYRLM(L+M+1) .NE. 1) THEN
            READ (LUSOL)
            GO TO 500
         END IF
C
C        Read R(l,m) in ao basis and transform to mo basis.
C        Extract active-active block in RLMAC = WRK(KRLMAC).
C        Electronic contribution TE(l,m).
C
         CALL READT(LUSOL,NNBASX,WRK(KRLMAO))
         CALL UTHU(WRK(KRLMAO),WRK(KRLM),WRK(KUCMO),WRK(KW21),
     *             NBAST,NORBT)
         IF (NASHT .GT. 0) THEN
            CALL GETAC2(WRK(KRLM),WRK(KRLMAC))
         END IF
         TLM   = TN(LM) - SOLELM(DV,WRK(KRLMAC),WRK(KRLM),TELMAC)
C
         JRXC  = KRXC
         DO 200 ICSIM = 1,NCSIM
            TTT = SLVELM(DTV(1,ICSIM),WRK(KRLMAC),WRK(KRLM),TAC)
            FACX  = D2 * FLVEC(LM) * TAC
            IF (ABS(FACX) .GT. THRZER)
     *         CALL DAXPY(NNORBX,FACX,WRK(KRLM),1,WRK(JRXC),1)
            JRXC  = JRXC + NNORBX
  200    CONTINUE
C
         FACY  = - D2 * FLVEC(LM) * TLM
         IF (INERSI) THEN
            FACY = FACY - D2*FLINR(LM) * TLMSI(LM)
         END IF
         TYCAC = TYCAC + FACY*TELMAC
         IF (ABS(FACY) .GT. THRZER)
     *      CALL DAXPY(NNORBX,FACY,WRK(KRLM),1,WRK(KRYC),1)
C
  500 CONTINUE
  520 CONTINUE
C
C
C     Calculate Rxc and Ryc contributions to SCVECS(NVAR,NCSIM)
C     =========================================================
C
C     Rxc contribution is gradient-like,
C     Ryc contribution is "lintrn"-like:
C
C     scvecs(i) = scvecs(i) + [ d(Rxc) / di ]
C                           + sum(j) [ d2(Ryc) / di dj ] * bcvecs(j)
C
C
C     ... CSF part of sigma vectors
C
      JRXC   = KRXC
      JRXCAC = KRXCAC
      DO 600 ICSIM = 1,NCSIM
         CALL GETAC2(WRK(JRXC),WRK(JRXCAC))
         IF (TRPLET) THEN
            TXC    = SOLELM(DVTR,WRK(JRXCAC),WRK(JRXC),TXCAC)
            TXC    = TXCAC
         ELSE
            TXC    = SOLELM(DV,WRK(JRXCAC),WRK(JRXC),TXCAC)
         END IF
         WRK(KTXCAC-1+ICSIM) = TXCAC
         JRXC   = JRXC   + NNORBX
         JRXCAC = JRXCAC + NNASHX
  600 CONTINUE
C
      CALL GETAC2(WRK(KRYC),WRK(KRYCAC))
C
      CALL SLVSC(NCSIM,0,NNASHX,BCVEC,CREF,SVEC,WRK(KRXCAC),WRK(KRYCAC),
     *           WRK(KTXCAC),TYCAC,INDXCI,WRK(KW30),LW30)
C     CALL SLVSC(NCSIM,NOSIM,NNASHX,BCVECS,CREF,SVECS,
C    *           RXAC,RYAC,TRXAC,TRYAC,INDXCI,WRK,LWRK)
C
      IF (IPRRSP.GT.101) THEN
         WRITE(LUPRI,*)' LINEAR TRANSFIRMED CONFIGURATION VECTOR'
         WRITE(LUPRI,*)' **** AFTER SLVSC in SLVLNC **** '
         CALL OUTPUT(SVEC,1,KZYVAR,1,NCSIM,KZYVAR,NCSIM,1,LUPRI)
      END IF
C
C     ... orbital part of sigma vector(s)
C
      IF (KZWOPT .GT. 0) THEN
         JRXC   = KRXC
         CALL DSPTSI(NORBT,WRK(KRYC),WRK(KURYC))
         DO 800 ICSIM = 1,NCSIM
            CALL DSPTSI(NORBT,WRK(JRXC),WRK(KURXC))
            IF (TRPLET) THEN
               CALL SLVSOR(.TRUE.,.FALSE.,1,UDVTR,
     *                     SVEC(1,ICSIM),WRK(KURXC))
            ELSE
               CALL SLVSOR(.TRUE.,.TRUE.,1,UDV,
     *                     SVEC(1,ICSIM),WRK(KURXC))
            END IF
            IF (IPRRSP.GT.101) THEN
               WRITE(LUPRI,*)' **** AFTER SLVSOR  in SLVLNC **** '
               WRITE(LUPRI,*)
     *         ' orbital part of LINEAR TRANSFORMED CONF VEC No',ICSIM
               WRITE(LUPRI,*)' TXC CONTRIBUTION'
               CALL OUTPUT(SVEC(1,ICSIM),1,KZYVAR,1,1,KZYVAR,1,1,LUPRI)
            END IF
            IF (TRPLET) THEN
               CALL SLVSOR(.FALSE.,.FALSE.,1,DTVTR(1,ICSIM),
     *                      SVEC(1,ICSIM),WRK(KURYC))
            ELSE
               CALL SLVSOR(.FALSE.,.FALSE.,1,DTV(1,ICSIM),
     *                     SVEC(1,ICSIM),WRK(KURYC))
            END IF
            IF (IPRRSP.GT.101) THEN
               WRITE(LUPRI,*)
     *         ' orbital part of LINEAR TRANSFORMED CONF VEC No',ICSIM
               WRITE(LUPRI,*)' TG CONTRIBUTION'
               CALL OUTPUT(SVEC(1,ICSIM),1,KZYVAR,1,1,KZYVAR,1,1,LUPRI)
            END IF

            JRXC   = JRXC   + NNORBX
  800    CONTINUE
         IF (IPRRSP.GT.101) THEN
            WRITE(LUPRI,*)' LINEAR TRANSFORMED CONFIGURATION VECTOR'
            WRITE(LUPRI,*)' **** AFTER SLVSOR  in SLVLNC **** '
            CALL OUTPUT(SVEC,1,KZYVAR,1,NCSIM,KZYVAR,NCSIM,1,LUPRI)
         END IF
      END IF
C
      CALL QEXIT('SLVLNC')
      RETURN
C     end of SLVLNC.
      END
C  /* Deck slvlno */
      SUBROUTINE SLVLNO(NCSIM,NOSIM,BOVECS,CREF,CMO,INDXCI,
     *                  FLVEC,FLINR,TLMSI,
     *                  UDV,DV,UDVTR,DVTR,SVEC,ISYRLM,
     *                  WRK,LFREE)
C
C
C  Purpose:  Calculate MCSCF E2 contribution from a
C            surrounding medium to an orbital trial vector.
C            Cavity radius = Rsol and dielectric constant = EPsol.
C
C
#include "implicit.h"
      DIMENSION BOVECS(*), CREF(*), CMO(*)
      DIMENSION INDXCI(*),   UDV(*),   DV(*),   ISYRLM(*)
      DIMENSION FLVEC(NLMSOL), FLINR(NLMSOL), TLMSI(NLMSOL)
      DIMENSION SVEC(KZYVAR,*),    WRK(*)
      DIMENSION UDVTR(*),DVTR(*)
C
#include "iratdef.h"
C
      PARAMETER ( D0 = 0.0D0 , D2 = 2.0D0 , THRZER = 1.0D-14 )
#include "dummy.h"
C
C
C  Used from common blocks:
C    INFINP : NLMSOL, LSOLMX, INERSI
C    INFORB : NNASHX, NNORBX, NNBASX, etc.
C    INFVAR : JWOP
C    INFRSP :
C    WRKRSP :
C    INFTAP : LUSOL,  LBSYMB
C    INFPRI : IPRSOL
C
#include "maxash.h"
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "inftap.h"
#include "infpri.h"
C
C
C     Statement functions
C
      TN(   LM) = WRK((KTNLM -1)+LM)
C
      CALL QENTER('SLVLNO')
C
C     Determine if full Hessian or only orbital Hessian
C
C
      IF (IPRRSP .GE. 40) THEN
         WRITE (LUPRI,'(//A)') ' --- TEST OUTPUT FROM SLVLNO ---'
      END IF
      IF (IPRRSP .GE. 140) THEN
         WRITE (LUPRI,'(/A)') ' --- SLVLNO - svec(1,nosim) on entry'
         CALL OUTPUT(SVEC,1,KZYVAR,1,NOSIM,KZYVAR,NOSIM,1,LUPRI)
      END IF
C
C
C     Core allocation
C
      KRXYO  = 1
      IF (TRPLET) THEN
         KRXYOT = KRXYO  + NOSIM*N2ORBX
         KTEXY  = KRXYOT + NOSIM*N2ORBX
      ELSE
         KTEXY  = KRXYO  + NOSIM*N2ORBX
      ENDIF
      KUBO   = KTEXY  + NOSIM
      KW10   = KUBO   + NOSIM*N2ORBX
C
      KTNLM  = KW10
      KUCMO  = KTNLM  + NLMSOL
      KRLMAC = KUCMO  + NORBT*NBAST
      KRLM   = KRLMAC + NNASHX
      KW20   = KRLM   + NNORBX
C     2.1 read rlmao in ao basis and transform to rlm in mo basis
      KRLMAO = KW20
      KURLM  = KRLMAO + NNBASX
      KURXLM = KURLM  + N2ORBX
      KURXAC = KURXLM + N2ORBX
      KW21   = KURXAC + N2ASHX
C
      KRXLM  = KW20
      KW22   = KRXLM  + NNORBX
C     3.0 SOLSC
      KRXYOA = KW10
      IF (TRPLET) THEN
         KRXYAT = KRXYOA + NOSIM*N2ASHX
         KTXYOA = KRXYAT + NOSIM*N2ASHX
      ELSE
         KTXYOA = KRXYOA + NOSIM*N2ASHX
      END IF
      IF (TRPLET) THEN
         KTXYAT = KTXYOA + NOSIM
         KOVLP  = KTXYAT + NOSIM
      ELSE
         KOVLP  = KTXYOA + NOSIM
      END IF
      KW30   = KOVLP  + NOSIM
      LW30   = LFREE  + 1 - KW30
C
      KNEED = MAX(KW21,KW22,KW30)
      IF (KNEED .GT. LFREE) CALL ERRWRK('SOLLNO',KNEED,LFREE)
C
C     Read nuclear contributions to ERLM (TN(l,m))
C     and position for R(l,m).
C
      REWIND LUSOL
      CALL MOLLAB('SOLVRLM ',LUSOL,LUERR)
      READ (LUSOL)
      CALL READT(LUSOL,NLMSOL,WRK(KTNLM))
C
C
C
C     Construct effective operators Rxo and Ryo
C     =========================================
C
C     Zero Rxo and Ryo storage
C
      CALL DZERO(WRK(KRXYO), NOSIM * N2ORBX )
      IF (TRPLET) THEN
         CALL DZERO(WRK(KRXYOT), NOSIM * N2ORBX )
      END IF
      CALL DZERO(WRK(KTEXY), NOSIM )
C
C     Unpack symmetry blocked CMO
C
      CALL UPKCMO(CMO,WRK(KUCMO))
C
C     Calculate unpacked orbital trial vectors in UBO
C
C
      IF (NOSIM.GT.0) THEN
         CALL RSPZYM(NOSIM,BOVECS,WRK(KUBO))
         CALL DSCAL(NOSIM*N2ORBX,-1.0D0,WRK(KUBO),1)
         IF (IPRRSP .GE. 55) THEN
            DO 210 IOSIM = 1,NOSIM
               JUBO = KUBO + (IOSIM-1)*N2ORBX
               WRITE (LUPRI,2110) IOSIM,NOSIM
               CALL OUTPUT(WRK(JUBO),1,NORBT,1,NORBT,NORBT,NORBT,1,
     &                     LUPRI)
  210       CONTINUE
         END IF
      END IF
 2110 FORMAT (/,' Orbital trial vector unpacked to matrix form (no.',
     *        I3,' of',I3,')')
C
C
C     Loop over l,m expansion.
C
      LM = 0
      DO 520 L = 0,LSOLMX
         READ (LUSOL) L1,(ISYRLM(M),M=1,2*L+1)
      DO 500 M = -L,L
         LM = LM + 1
         IF (ISYRLM(L+M+1) .NE. KSYMOP .AND. ISYRLM(L+M+1) .NE. 1) THEN
            READ (LUSOL)
            GO TO 500
         END IF
C
C        Read R(l,m) in ao basis and transform to mo basis.
C        Extract active-active block in RLMAC = WRK(KRLMAC).
C        Electronic contribution TE(l,m).
C
C
         CALL READT(LUSOL,NNBASX,WRK(KRLMAO))
         CALL UTHU(WRK(KRLMAO),WRK(KRLM),WRK(KUCMO),WRK(KW21),
     *             NBAST,NORBT)
         IF (NASHT .GT. 0) THEN
            CALL GETAC2(WRK(KRLM),WRK(KRLMAC))
         END IF
         TELM  = SOLELM(DV,WRK(KRLMAC),WRK(KRLM),TELMAC)
         TLM   = TN(LM) - TELM
C
         IF (IPRRSP .GE. 40) THEN
            WRITE (LUPRI,'(/A,I5,4(/A,F15.10))')
     *         ' LM       ',LM,
     *         ' TELM     ',TELM,
     *         ' TELMAC   ',TELMAC,
     *         ' TLM      ',TLM,
     *         ' FLVEC(LM)',FLVEC(LM)
         END IF
         IF (IPRRSP .GE. 70) THEN
            WRITE (LUPRI,'(/A)') ' Rlm_mo matrix:'
            CALL OUTPAK(WRK(KRLM),  NORBT,1,LUPRI)
            IF (NASHT .GT. 0) THEN
               WRITE (LUPRI,'(/A)') ' Rlm_ac matrix:'
               CALL OUTPAK(WRK(KRLMAC),NASHT,1,LUPRI)
            END IF
         END IF
C
         CALL DSPTSI(NORBT,WRK(KRLM),WRK(KURLM))
         FACY  = - D2 * FLVEC(LM) * TLM
         IF (INERSI) THEN
            FACY  = FACY - D2*FLINR(LM) * TLMSI(LM)
         END IF
         DO 200 IOSIM = 1,NOSIM
            JRXYO = KRXYO + (IOSIM-1)*N2ORBX
            IF (TRPLET) THEN
               JRXYOT = KRXYOT + (IOSIM-1)*N2ORBX
            END IF
            JUBO = KUBO + (IOSIM-1)*N2ORBX
C
C           Calculate one-index transformed integrals
C
            CALL DZERO(WRK(KURXLM),N2ORBX)
            CALL ONEXH1(WRK(JUBO),WRK(KURLM),WRK(KURXLM))
C           CALL DGETSP(NORBT,WRK(KURXLM),WRK(KRXLM))
C           TODO: what if LSYMPT .ne. 1 ???
C           900720: I think OK, i.e. that only RLM of sym 1
C           such that RXLM is of sym LSYMPT is needed.
C           CALL TR1UH1(UBOVEC,H1,H1X,IH1SYM)
C
            IF (NASHT .GT. 0) CALL GETACQ(WRK(KURXLM),WRK(KURXAC))
            IF (IPRRSP .GE. 60) THEN
               WRITE (LUPRI,'(/A,I5)') ' Rlm_X_mo matrix, IOSIM =',IOSIM
               CALL OUTPUT(WRK(KURXLM),1,NORBT,1,NORBT,NORBT,NORBT,
     *                                 1,LUPRI)
               IF (NASHT .GT. 0) THEN
                  WRITE (LUPRI,'(/A)') ' Rlm_X_ac matrix:'
                  CALL OUTPUT(WRK(KURXAC),1,NASHT,1,NASHT,NASHT,NASHT,
     *                                 1,LUPRI)
               END IF
            END IF
            IF (TRPLET) THEN
               FACX  = D2 * FLVEC(LM)
     *                 * SLVTLM(UDVTR,WRK(KURXAC),WRK(KURXLM),TAC)
            ELSE
               FACX  = D2 * FLVEC(LM)
     *                 * SLVQLM(UDV,WRK(KURXAC),WRK(KURXLM),TAC)
            END IF
            WRK(KTEXY-1+IOSIM) = WRK(KTEXY-1+IOSIM) + FACX * TELM
            IF (IPRRSP .GE. 40)
     *         WRITE (LUPRI,'(/A,2I5,2F15.10)')
     *         ' LM, IOSIM, FACX, FACY :',LM,IOSIM,FACX,FACY
            IF (ABS(FACX) .GT. THRZER) THEN
               CALL DAXPY(N2ORBX,FACX,WRK(KURLM), 1,WRK(JRXYO),1)
            END IF
            IF (ABS(FACY) .GT. THRZER) THEN
               IF (TRPLET) THEN
                  CALL DAXPY(N2ORBX,FACY,WRK(KURXLM),1,WRK(JRXYOT),1)
               ELSE
                  CALL DAXPY(N2ORBX,FACY,WRK(KURXLM),1,WRK(JRXYO),1)
               END IF
            END IF
  200    CONTINUE
C
  500 CONTINUE
  520 CONTINUE
C
C
      IF (IPRRSP .GE. 50) THEN
         IF (TRPLET) THEN
            DO 603 IOSIM = 1,NOSIM
               JRXYO = KRXYO + (IOSIM-1)*N2ORBX
               WRITE (LUPRI,'(/2A,I3,A,I3)')
     *         ' --- SLVLNO - TRIPLET ',
     *         '(Rxo) matrix no.',IOSIM,' of',NOSIM
               CALL OUTPUT(WRK(JRXYO),1,NORBT,1,NORBT,NORBT,NORBT,
     *                              1,LUPRI)
               JRXYOT = KRXYOT + (IOSIM-1)*N2ORBX
               WRITE (LUPRI,'(/A,I3,A,I3)')
     *         ' --- SLVLNO - (Ryo) matrix no.',IOSIM,' of',NOSIM
               CALL OUTPUT(WRK(JRXYOT),1,NORBT,1,NORBT,NORBT,NORBT,
     *                              1,LUPRI)
  603       CONTINUE
         ELSE
            DO 600 IOSIM = 1,NOSIM
               JRXYO = KRXYO + (IOSIM-1)*N2ORBX
               WRITE (LUPRI,'(/A,I3,A,I3)')
     *         ' --- SLVLNO - (Rxo + Ryo) matrix no.',IOSIM,' of',NOSIM
               CALL OUTPUT(WRK(JRXYO),1,NORBT,1,NORBT,NORBT,NORBT,
     *                              1,LUPRI)
  600       CONTINUE
         END IF
      END IF
      IF (.NOT.TDHF) THEN
C
C
C
C        Calculate Rxo and Ryo contributions to SVEC(NSVEC,NOSIM)
C        =========================================================
C
C        Rxo and Ryo contributions are gradient-like:
C
C        svec(i) = svec(i) + [ d(Rxo) / di ] + [ d(Ryo) / di ]
C             = svec(i) + [ d(Rxo + Ryo) / di ]
C
C        Note: the 1/2 sum(t) [ gsol(ts)b(rt) - gsol(tr)b(st) ]
C           correction to the one-index transformation in Ryo
C           is added in lintrn, where g(rs) = gstd(rs) + gsol(rs).
C
C        ... CSF part of sigma vectors
C
         DO 700 IOSIM = 1,NOSIM
            JRXYO  = KRXYO  + (IOSIM-1)*N2ORBX
            JRXYOA = KRXYOA + (IOSIM-1)*N2ASHX
            IF (TRPLET) THEN
               JRXYOT = KRXYOT + (IOSIM-1)*N2ORBX
               JRXYAT = KRXYAT + (IOSIM-1)*N2ASHX
            END IF
            IF (IREFSY .EQ. KSYMST) THEN
               WRK(KOVLP-1+IOSIM) = DDOT(KZCONF,CREF,1,SVEC(1,IOSIM),1)
            ELSE
               WRK(KOVLP-1+IOSIM) = D0
            END IF
            CALL GETACQ(WRK(JRXYO),WRK(JRXYOA))
            TXYO   = SLVQLM(UDV,WRK(JRXYOA),WRK(JRXYO),TXYOAC)
            WRK(KTXYOA-1+IOSIM) = TXYOAC
            IF (TRPLET) THEN
               CALL GETACQ(WRK(JRXYOT),WRK(JRXYAT))
               TXYOT  = SLVQLM(UDVTR,WRK(JRXYAT),WRK(JRXYOT),TXYACT)
               WRK(KTXYAT-1+IOSIM) = TXYACT
               IF (IPRRSP .GE. 40) THEN
                  WRITE (LUPRI,*) ' TRIPLET CALCULATION'
                  WRITE (LUPRI,'(/A,I5,3F15.10)')
     *            ' IOSIM, C_OVLP, TXYOAC, TXYO :',
     *            IOSIM,WRK(KOVLP-1+IOSIM),TXYOAC,TXYO
                  WRITE (LUPRI,'(/A,I5,3F15.10)')
     *            ' IOSIM, C_OVLP, TXYACT, TXYOT :',
     *            IOSIM,WRK(KOVLP-1+IOSIM),TXYACT,TXYOT
               END IF
            ELSE
               IF (IPRRSP .GE. 40) THEN
                  WRITE (LUPRI,'(/A,I5,3F15.10)')
     *            ' IOSIM, C_OVLP, TXYOAC, TXYO :',
     *            IOSIM,WRK(KOVLP-1+IOSIM),TXYOAC,TXYO
               END IF
            END IF
  700    CONTINUE
         CALL SLVSC(0,NOSIM,N2ASHX,DUMMY,CREF,SVEC,WRK(KRXYOA),
     *              WRK(KRXYAT),
     *              WRK(KTXYOA),DUMMY,INDXCI,WRK(KW30),LW30)
C        CALL SLVSC(NCSIM,NOSIM,N2ASHX,BCVECS,CREF,SVECS,
C    *              RXAC,RYAC,TRXAC,TRYAC,INDXCI,WRK,LWRK)
C
         IF (IPRRSP.GT.101) THEN
            WRITE(LUPRI,*)' LINEAR TRANSFORMED ORBITAL VECTOR'
            WRITE(LUPRI,*)' **** AFTER SLVSC in SLVLNO **** '
            CALL OUTPUT(SVEC,1,KZYVAR,1,NOSIM,KZYVAR,NOSIM,1,LUPRI)
         END IF
      END IF
      IF (TRPLET) THEN
         CALL SLVSOR(.TRUE.,.FALSE.,NOSIM,UDVTR,SVEC(1,1),WRK(KRXYO))
         CALL SLVSOR(.TRUE.,.TRUE.,NOSIM,UDV,SVEC(1,1),WRK(KRXYOT))
      ELSE
         CALL SLVSOR(.TRUE.,.TRUE.,NOSIM,UDV,SVEC(1,1),WRK(KRXYO))
      ENDIF
      IF (IPRRSP.GT.101) THEN
         WRITE(LUPRI,*)' LINEAR TRANSFORMED ORBITAL VECTOR'
         WRITE(LUPRI,*)' **** AFTER SLVSOR in SLVLNO **** '
         CALL OUTPUT(SVEC,1,KZYVAR,1,NOSIM,KZYVAR,NOSIM,1,LUPRI)
      END IF
      IF ( KZCONF.GT.0 .AND. IREFSY.EQ.KSYMST) THEN
         IF (NCREF .NE. KZCONF) CALL QUIT('SLVLNO: NCREF .ne. KZCONF')
C
C        ... test orthogonality
C
         DO 900 IOSIM = 1,NOSIM
            TEST = DDOT(KZCONF,CREF,1,SVEC(1,IOSIM),1)
     *           - WRK(KOVLP-1+IOSIM)
            IF (ABS(TEST) .GT. 1.D-8) THEN
               NWARN = NWARN + 1
               WRITE (LUPRI,'(/A,I5,/A,1P,D12.4)')
     *            ' --- SLVLNO WARNING, for IOSIM =',IOSIM,
     *            ' <CREF | SVEC_solvent(iosim) > =',TEST
            END IF
  900    CONTINUE
      END IF
C
C        ... test print
C
      IF (IPRRSP .GE. 140) THEN
         WRITE (LUPRI,'(/A)') ' --- SLVLNO - svec(ci,1) on exit'
         DO 930 I = 1,KZCONF
            IF (SVEC(I,1) .NE. D0) WRITE (LUPRI,'(A,I10,F15.10)')
     *         ' conf #',I,SVEC(I,1)
  930    CONTINUE
      END IF
      IF (IPRRSP .GE. 140) THEN
         WRITE (LUPRI,'(/A)') ' --- SLVLNO - svec(orb,1) on exit'
         WRITE (LUPRI,'(/A)') ' Z - PART OR VECTOR '
         CALL OUTPUT(SVEC(KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,1,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' Y - PART OR VECTOR '
         CALL OUTPUT(SVEC(KZVAR+KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,
     *               1,1,LUPRI)
      END IF
C
      CALL QEXIT('SLVLNO')
      RETURN
C     ... end of SLVLNO.
      END
C  /* Deck slvsc */
      SUBROUTINE SLVSC(NCSIM,NOSIM,NXASHX,BCVECS,CREF,SVECS,
     *                 RXAC,RYAC,TRXAC,TRYAC,XNDXCI,WRK,LWRK)
C
C   5-May-1987 Hans Joergen Aa. Jensen
C   Rewritten for determinant CI 19-Jul-1990 hjaaj
C
C   Purpose:
C      Perform CI transformation of a one-electron operator,
C      as the solvent operator for "cavity in dielectric" model.
C
C
C  Input:   ordinary and 1-index
C           transformed inactive solvent-matrices (RXAC and RYAC),
C           trial and reference CI-vectors (BCVECS and CREF)
C
C  Output:  New CI-part of sigma vector (SVECS)
C
C Scratch:  WRK
C
#include "implicit.h"
      DIMENSION BCVECS(KZCONF,*),  SVECS(KZYVAR,*)
      DIMENSION CREF(*),  RXAC(NXASHX,*),RYAC(NNASHX)
      DIMENSION TRXAC(*), XNDXCI(*), WRK(LWRK)
C
      PARAMETER ( D2 = 2.0D0 , DM1 = -1.0D0 )
C     ... note: factors of 2 (D2) below takes care of the "2" in
C         2 ( <j / R / B> - T bj )
C
C Used from common blocks:
C   INFORB : NNASHX,N2ASHX,NASHT
C   INFRSP : IREFSY,KSYMST,KSYMOP,NCREF,KZCONF
C
#include "priunit.h"
#include "inforb.h"
#include "infpri.h"
#include "infrsp.h"
#include "wrkrsp.h"
C
      CALL QENTER('SLVSC ')
C
C     ... input check
C
      IF (NCSIM .GT. 0 .AND. NOSIM .GT. 0) THEN
         WRITE (LUERR,'(//A/A,2I5)')
     *   ' --- SLVSC ERROR, both NCSIM .gt. 0 and NOSIM .gt. 0',
     *   '     NCSIM, NOSIM =',NCSIM,NOSIM
         CALL QTRACE(LUERR)
         CALL QUIT('SLVSC ERROR, both NCSIM .gt. 0 and NOSIM .gt. 0')
      END IF
C
      IF (NCSIM .GT. 0) THEN
         KURYAC = 1
         KURXAC = KURYAC + N2ASHX
         KW1    = KURXAC + N2ASHX
         LW1    = LWRK + 1 - KW1
         CALL DSPTSI(NASHT,RYAC,WRK(KURYAC))
C        CALL DSCAL(N2ASHX,DM1,WRK(KURYAC),1)
         DO 18 ICSIM = 1,NCSIM
            CALL DSPTSI(NASHT,RXAC(1,ICSIM),WRK(KURXAC))
C           CALL DSCAL(N2ASHX,DM1,WRK(KURXAC),1)
C rxac defined with opposite sign relative to SOLSC
C
            IF (IPRRSP.GT.100) THEN
               WRITE(LUPRI,*)' Z part,RXAC MATRIX, ICSIM=', ICSIM
               CALL OUTPUT(WRK(KURXAC),1,NASHT,1,NASHT,
     *                                 NASHT,NASHT,1,LUPRI)
            END IF
C           First SVECS(NA) = SVECS(NA) + < NA | 2*RX | CREF >
C
C Z PART OF LINEAR TRANSFORMATION
C
            JSPIN1 = 0
            JSPIN2 = 0
            CALL CISIGD(IREFSY,KSYMST,NCREF,KZCONF,CREF,SVECS(1,ICSIM),
     &                  WRK(KURXAC),DUMMY,
     &                  .TRUE.,.FALSE.,XNDXCI,JSPIN1,JSPIN2,
     &                  WRK(KW1),LW1)
            IF (IPRRSP.GT.101) THEN
               WRITE(LUPRI,*)
     *       ' TXC Z CONF  PART OF LINEAR TRANSFORMED CONF TRIAL VECTOR'
               CALL OUTPUT(SVECS(1,ICSIM),1,KZCONF,1,1,KZCONF,1,1,LUPRI)
            END IF
C           then  SVECS(NA) = SVECS(NA) + < NA | 2*RY | BCVECS >
            IF (IPRRSP.GT.100) THEN
               WRITE(LUPRI,*)' Z part,-1.0*RYAC MATRIX, ICSIM=', ICSIM
               CALL OUTPUT(WRK(KURYAC),1,NASHT,1,NASHT,
     *                                 NASHT,NASHT,1,LUPRI)
            END IF
            JSPIN1 = 0
            JSPIN2 = 0
            CALL CISIGD(KSYMST,KSYMST,KZCONF,KZCONF,
     &                  BCVECS(1,ICSIM),SVECS(1,ICSIM),
     &                  WRK(KURYAC),DUMMY,
     &                  .TRUE.,.FALSE.,XNDXCI,JSPIN1,JSPIN2
     &                  ,WRK(KW1),LW1)
            IF (IPRRSP.GT.101) THEN
               WRITE(LUPRI,*)
     *       ' TG Z CONF  PART OF LINEAR TRANSFORMED CONF TRIAL VECTOR'
               CALL OUTPUT(SVECS(1,ICSIM),1,KZCONF,1,1,KZCONF,1,1,LUPRI)
            END IF
            DO 200 NA = 1,KZCONF
               SVECS(NA,ICSIM) = SVECS(NA,ICSIM)
     *                         -   TRYAC        * BCVECS(NA,ICSIM)
C    *                         - D2 * TRXAC(ICSIM) * CREF(NA)
  200          CONTINUE
            IF (IPRRSP.GT.101) THEN
               WRITE(LUPRI,*)
     *       ' TG Z CONF  PART OF LINEAR TRANSFORMED CONF TRIAL VECTOR'
               WRITE(LUPRI,*)' <0/Tg/0> S(j) subtracted '
               CALL OUTPUT(SVECS(1,ICSIM),1,KZCONF,1,1,KZCONF,1,1,LUPRI)
            END IF
            IF (IREFSY .EQ. KSYMST) THEN
               FAC = DDOT(KZCONF,SVECS(1,ICSIM),1,CREF,1)
               CALL DAXPY(KZCONF,(-FAC),CREF,1,SVECS(1,ICSIM),1)
               IF (IPRRSP.GT.101) THEN
                  WRITE(LUPRI,*)
     *            ' TG+TXC Z CONF  PART OF LINEAR TRANSFORMED',
     *            'CONF TRIAL VECTOR'
                  WRITE(LUPRI,*)' REFERENCE COMPONENT PROJECTED OUT'
                  CALL OUTPUT(SVECS(1,ICSIM),1,KZCONF,1,1,KZCONF,1,
     *                                                      1,LUPRI)
               END IF
            END IF
C
C Y PART OF LINEAR TRANSFORMATION
C
            CALL DGETRN(WRK(KURXAC),NASHT,NASHT)
            CALL DSCAL(N2ASHX,DM1,WRK(KURXAC),1)
            IF (IPRRSP.GT.100) THEN
               WRITE(LUPRI,*)' Y PART, -RXAC MATRIX, '//
     *         ' ICSIM=', ICSIM
               CALL OUTPUT(WRK(KURXAC),1,NASHT,1,NASHT,
     *                                 NASHT,NASHT,-1,LUPRI)
            END IF
            JSPIN1 = 0
            JSPIN2 = 0
            CALL CISIGD(IREFSY,KSYMST,NCREF,KZCONF,CREF,
     &                  SVECS(1+KZVAR,ICSIM),WRK(KURXAC),DUMMY,
     &                  .TRUE.,.FALSE.,XNDXCI,JSPIN1,JSPIN2,
     &                  WRK(KW1),LW1)
            IF (IPRRSP.GT.101) THEN
               WRITE(LUPRI,*)
     *       ' TXC Y CONF  PART OF LINEAR TRANSFORMED CONF TRIAL VECTOR'
               CALL OUTPUT(SVECS(1+KZVAR,ICSIM),1,KZCONF,1,1,
     *                                            KZCONF,1,1,LUPRI)
            END IF
C           then  SVECS(NA) = SVECS(NA) + < NA | 2*RY | BCVECS >
C
C <0L| CONTRIBUTION VANISH SINCE CONF TRIAL VECTORS HAVE ZERO Y
C COMPONENT
C
C           CALL CISIGD(IREFSY,KSYMST,NCREF,KZCONF,
C    &                  BCVECS(1,ICSIM),SVECS(1,ICSIM),
C    &                  WRK(KURYAC),DUMMY,
C    &                  .TRUE.,.FALSE.,XNDXCI,0,0,WRK(KW1),LW1)
            IF (IREFSY .EQ. KSYMST) THEN
               FAC = DDOT(KZCONF,SVECS(1+KZVAR,ICSIM),1,CREF,1)
               CALL DAXPY(KZCONF,(-FAC),CREF,1,SVECS(1+KZVAR,ICSIM),1)
C              DO 201 NA = 1,KZCONF
C                 SVECS(NA,ICSIM) = SVECS(NA,ICSIM)
C    *                            - D2 * TRXAC(ICSIM) * CREF(NA)
C NO Y COMPONENT IN CONF TRIAL VECTORS
C    *                            - D2 * TRYAC        * BCVECS(NA,ICSIM)
C 201          CONTINUE
               IF (IPRRSP.GT.101) THEN
                  WRITE(LUPRI,*)
     *            ' TXC Y CONF  PART OF LINEAR TRANSFORMED',
     *            'CONF TRIAL VECTOR'
                  WRITE(LUPRI,*)' REFERENCE COMPONENT PROJECTED OUT'
                  CALL OUTPUT(SVECS(1+KZVAR,ICSIM),1,KZCONF,1,1,
     *                                             KZCONF,1,1,LUPRI)
               END IF
            ELSE
C NO Y COMPONENT IN CONF TRIAL VECTORS
C
C              FAC = - D2 * TRYAC
C              CALL DAXPY(KZCONF,FAC,BCVECS(1,ICSIM),1,SVECS(1,ICSIM),1)
            END IF
   18    CONTINUE
      END IF
      IF (NOSIM .GT. 0) THEN
         KURXAC = 1
         KW1    = KURXAC + N2ASHX
         LW1    = LWRK + 1 - KW1
         DO 22 IOSIM = 1,NOSIM
C           CALL DSCAL(N2ASHX,D2,RXAC(1,IOSIM),1)
            IF (IPRRSP.GT.100) THEN
               WRITE(LUPRI,*)' Z PART,RXAC MATRIX,
     *         IOSIM=', IOSIM
               CALL OUTPUT(RXAC(1,IOSIM),1,NASHT,1,NASHT,
     *                                 NASHT,NASHT,1,LUPRI)
            END IF
C
C Z PART OF LINEAR TRANSFORMATION
C
            CALL DGETRN(RXAC(1,IOSIM),NASHT,NASHT)
            JSPIN1 = 0
            JSPIN2 = 0
            CALL CISIGD(IREFSY,KSYMST,NCREF,KZCONF,CREF,SVECS(1,IOSIM),
     &               RXAC(1,IOSIM),DUMMY,
     &               .TRUE.,.FALSE.,XNDXCI,JSPIN1,JSPIN2,WRK(KW1),LW1)
            IF (IPRRSP.GT.101) THEN
               WRITE(LUPRI,*)
     *       ' TXO+TYO Z CONF  PART OF LINEAR TRANSF ORB TRIAL VECTOR'
               CALL OUTPUT(SVECS(1,IOSIM),1,KZCONF,1,1,KZCONF,1,1,LUPRI)
            END IF
C           Remove CREF component of SVECS (inactive contribution
C           not needed because of this projection)
            IF (IREFSY .EQ. KSYMST) THEN
               FAC = DDOT(KZCONF,SVECS(1,IOSIM),1,CREF,1)
               CALL DAXPY(KZCONF,(-FAC),CREF,1,SVECS(1,IOSIM),1)
               IF (IPRRSP.GT.101) THEN
                  WRITE(LUPRI,*)
     *            ' TXO+TYO Z CONF  PART OF LINEAR TRANSFORMED',
     *            'ORB TRIAL VECTOR'
                  WRITE(LUPRI,*)' REFERENCE COMPONENT PROJECTED OUT'
                  CALL OUTPUT(SVECS(1,IOSIM),1,KZCONF,1,1,KZCONF,1,
     *                                                      1,LUPRI)
               END IF
            END IF
            IF (TRPLET) THEN
C              CALL DSCAL(N2ASHX,D2,RXAC(1,IOSIM+NOSIM),1)
               IF (IPRRSP.GT.100) THEN
                  WRITE(LUPRI,*)' Z PART triplet, RXAC MATRIX, '//
     *            'IOSIM=', IOSIM
                  CALL OUTPUT(RXAC(1,IOSIM+NOSIM),1,NASHT,1,NASHT,
     *                                 NASHT,NASHT,-1,LUPRI)
               END IF
C
C Z PART OF LINEAR TRANSFORMATION
C
               CALL DGETRN(RXAC(1,IOSIM+NOSIM),NASHT,NASHT)
               JSPIN1 = 1
               JSPIN2 = 0
               CALL CISIGD(IREFSY,KSYMST,NCREF,KZCONF,CREF,
     &               SVECS(1,IOSIM),
     &               RXAC(1,IOSIM+NOSIM),DUMMY,
     &               .TRUE.,.FALSE.,XNDXCI,JSPIN1,JSPIN2,WRK(KW1),LW1)
               IF (IPRRSP.GT.101) THEN
                  WRITE(LUPRI,*)
     *          ' TXO+TYO Z CONF  TRIPLET PART OF LINEAR TRANSF',
     *          '  ORB TRIAL VECTOR'
                  CALL OUTPUT(SVECS(1,IOSIM),1,KZCONF,1,1,KZCONF,1,
     *                                                    1,LUPRI)
               END IF
            END IF
C
C Y PART OF LINEAR TRANSFORMATION
C
            CALL DGETRN(RXAC(1,IOSIM),NASHT,NASHT)
            CALL DSCAL(N2ASHX,DM1,RXAC(1,IOSIM),1)
            IF (IPRRSP.GT.100) THEN
               WRITE(LUPRI,*)' Y PART,RXAC MATRIX, SCALED WITH -1.0
     *         IOSIM=', IOSIM
               CALL OUTPUT(RXAC(1,IOSIM),1,NASHT,1,NASHT,
     *                                 NASHT,NASHT,1,LUPRI)
            END IF
C           now SVECS(NA) = SVECS(NA) + <NA | 2*RXY | CREF >
            JSPIN1 = 0
            JSPIN2 = 0
            CALL CISIGD(IREFSY,KSYMST,NCREF,KZCONF,CREF,
     &               SVECS(KZVAR+1,IOSIM),
     &               RXAC(1,IOSIM),DUMMY,
     &               .TRUE.,.FALSE.,XNDXCI,JSPIN1,JSPIN2,
     &                WRK(KW1),LW1)
            IF (IPRRSP.GT.101) THEN
               WRITE(LUPRI,*)
     *      ' TXO+TYO Y CONF  PART OF LINEAR TRANSF ORB TRIAL VECTOR'
               CALL OUTPUT(SVECS(1+KZVAR,IOSIM),1,KZCONF,1,1,
     *                                          KZCONF,1,1,LUPRI)
            END IF
C           Remove CREF component of SVECS (inactive contribution
C           not needed because of this projection)
            IF (IREFSY .EQ. KSYMST) THEN
               FAC = DDOT(KZCONF,SVECS(1+KZVAR,IOSIM),1,CREF,1)
               CALL DAXPY(KZCONF,(-FAC),CREF,1,SVECS(1+KZVAR,IOSIM),1)
               IF (IPRRSP.GT.101) THEN
                  WRITE(LUPRI,*)
     *            ' TXO+TYO Y CONF  PART OF LINEAR TRANSFORMED',
     *            ' ORB TRIAL VECTOR'
                  WRITE(LUPRI,*)' REFERENCE COMPONENT PROJECTED OUT'
                  CALL OUTPUT(SVECS(1+KZVAR,IOSIM),
     *                        1,KZCONF,1,1,KZCONF,1,1,LUPRI)
               END IF
            END IF
            IF (TRPLET) THEN
               CALL DGETRN(RXAC(1,IOSIM+NOSIM),NASHT,NASHT)
               CALL DSCAL(N2ASHX,DM1,RXAC(1,IOSIM+NOSIM),1)
               IF (IPRRSP.GT.100) THEN
                  WRITE(LUPRI,*)' Y PART triplet ,RXAC MATRIX,',
     *            'SCALED WITH -1.0
     *            IOSIM=', IOSIM
                  CALL OUTPUT(RXAC(1,IOSIM+NOSIM),1,NASHT,1,NASHT,
     *                                 NASHT,NASHT,1,LUPRI)
               END IF
C              now SVECS(NA) = SVECS(NA) + <NA | 2*RXY | CREF >
               JSPIN1 = 1
               JSPIN2 = 0
               CALL CISIGD(IREFSY,KSYMST,NCREF,KZCONF,CREF,
     &               SVECS(KZVAR+1,IOSIM),
     &               RXAC(1,IOSIM+NOSIM),DUMMY,
     &               .TRUE.,.FALSE.,XNDXCI,JSPIN1,JSPIN2,
     &               WRK(KW1),LW1)
               IF (IPRRSP.GT.101) THEN
               WRITE(LUPRI,*)
     *         ' TXO+TYO Y CONF TRIPLET  PART OF LINEAR TRANSF',
     *         ' ORB TRIAL VECTOR'
               CALL OUTPUT(SVECS(1+KZVAR,IOSIM),1,KZCONF,1,1,
     *                                          KZCONF,1,1,LUPRI)
               END IF
            END IF
   22    CONTINUE
      END IF
C
C     CALL CISIGD(ICSYM,IHCSYM,NCDET,NHCDET, C,HC, FCAC,H2AC,
C    &            NOH2,IH8SM,XNDXCI,ISPIN1,ISPIN2,WRK,LFREE)
C
      CALL QEXIT('SLVSC ')
      RETURN
      END
C  /* Deck slvqlm */
      FUNCTION SLVQLM(DV,RLMAC,RLM,TELMAC)
C
C   30-Nov-1986 Hans Jorgen Aa. Jensen
C   l.r. 6-May-1990 hjaaj
C
C QUADRATIC PACKED DENSITY MATRIX AND ACTIVE INTEGRALS
C
#include "implicit.h"
      DIMENSION DV(*), RLMAC(*), RLM(NORBT,*)
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0 )
C
C Used from common blocks:
C  INFORB: NISHT,NASHT
C  INFIND: IROW(*), ISX(*)
C
#include "maxash.h"
#include "maxorb.h"
#include "inforb.h"
#include "infind.h"
C
      IF (NASHT .GT. 0) THEN
         TELM = DDOT(N2ASHX,DV,1,RLMAC,1)
      ELSE
         TELM   = D0
      END IF
      TELMAC = TELM
C
      DO 300 IW = 1,NISHT
         I = ISX(IW)
         TELM = TELM + D2*RLM(I,I)
  300 CONTINUE
C
      SLVQLM = TELM
      RETURN
      END
C  /* Deck slvelm */
      FUNCTION SLVELM(DTV,RLMAC,RLM,TELMAC)
C
C   30-Nov-1986 Hans Jorgen Aa. Jensen
C   l.r. 6-May-1990 hjaaj
C   6DEC92 KM + PJ SOLELM MODIFIED TO QUADRATIC DTV
C
#include "implicit.h"
      DIMENSION DTV(NASHT,*), RLMAC(*), RLM(*)
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0 )
C
C Used from common blocks:
C  INFORB: NISHT,NASHT
C  INFIND: IROW(*), ISX(*)
C
#include "maxash.h"
#include "maxorb.h"
#include "inforb.h"
#include "infind.h"
C
      TELM =  D0
      IF (NASHT .GT. 0) THEN
         DO 100 I = 1,NASHT
            TELM = TELM + DTV(I,I)*RLMAC(I*(I+1)/2)
            DO 200 J = 1,I-1
               IJ = I*(I-1)/2+J
               TELM = TELM +(DTV(I,J)+DTV(J,I))*RLMAC(IJ)
  200       CONTINUE
  100    CONTINUE
      END IF
      TELMAC = TELM
C
      DO 300 IW = 1,NISHT
         I = ISX(IW)
         TELM = TELM + D2*RLM(IROW(I+1))
  300 CONTINUE
C
      SLVELM = TELM
      RETURN
      END
C  /* Deck slvtlm */
      FUNCTION SLVTLM(DV,RLMAC,RLM,TELMAC)
C
C   30-Nov-1986 Hans Jorgen Aa. Jensen
C   l.r. 6-May-1990 hjaaj
C
C QUADRATIC PACKED DENSITY MATRIX AND ACTIVE INTEGRALS
C
#include "implicit.h"
      DIMENSION DV(*), RLMAC(*), RLM(NORBT,*)
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0 )
C
C Used from common blocks:
C  INFORB: NISHT,NASHT
C  INFIND: IROW(*), ISX(*)
C
#include "maxash.h"
#include "maxorb.h"
#include "inforb.h"
#include "infind.h"
C
      IF (NASHT .GT. 0) THEN
         TELM = DDOT(N2ASHX,DV,1,RLMAC,1)
      ELSE
         TELM   = D0
      END IF
      TELMAC = TELM
C
C
      SLVTLM = TELM
      RETURN
      END
C  /* Deck getacq */
      SUBROUTINE GETACQ(HH,HHAC)
C
C  6-May-1990 Hans Joergen Aa. Jensen
C
C Module : SIRGP
C
C Purpose: Extract NNASHX block with active-active orbital indices out
C          of NNORBX matrix.
C
#include "implicit.h"
      DIMENSION HH(NORBT,NORBT), HHAC(NASHT,NASHT)
C
C INFORB : NORBT, NASHT, NNORBX,NNASHX
C INFIND : IOBTYP(*), ICH(*)
C
#include "maxash.h"
#include "maxorb.h"
#include "inforb.h"
#include "infind.h"
C
      CALL QENTER('GETACQ')
      IF (NASHT .EQ. 0) GO TO 9999
C
      DO 200 J = 1,NORBT
      IF (IOBTYP(J) .EQ. JTACT) THEN
         NJ  = ICH(J)
         DO 100 I = 1,NORBT
         IF (IOBTYP(I) .EQ. JTACT) THEN
            NI = ICH(I)
            HHAC(NI,NJ) = HH(I,J)
         END IF
  100    CONTINUE
      END IF
  200 CONTINUE
C
 9999 CALL QEXIT('GETACQ')
      RETURN
C
C     End of GETACQ.
C
      END
C  /* Deck onexh1 */
      SUBROUTINE ONEXH1(OITMAT,H1,H1X)
C
C
C One-index transform H1(norbt,norbt)
C to H1X(norbt,norbt) using OITMAT(norbt,norbt)
C
C
C Input : H1(a,b) of symmetry IH1SYM
C Output: H1X(a,b) = H1X(a,b)
C                  + H1(a,b) one-index transformed with OITMAT
C                  = H1X(a,b)
C                  + sum(c) [ OITMAT(a,c) H1(c,b)
C                           - H1(a,c) OITMAT(c,b) ]
C
#include "implicit.h"
      DIMENSION OITMAT(NORBT,NORBT), H1(NORBT,NORBT), H1X(NORBT,NORBT)
C
C Used from common blocks:
C  INFORB : NSYM, MULD2H, NORBT, NORB(*), IORB(*)
C
#include "inforb.h"
C
C
C One-index transform both indices and add to previous content.
C
C  (a~ b~) =  SUM(c) [ OITMAT(a,c)*(c b)
C                    - (a c) OITMAT(c,b) ]
C
      CALL DGEMM('N','N',NORBT,NORBT,NORBT,1.D0,
     &           OITMAT,NORBT,
     &           H1,NORBT,1.D0,
     &           H1X,NORBT)
      CALL DGEMM('N','N',NORBT,NORBT,NORBT,-1.D0,
     &           H1,NORBT,
     &           OITMAT,NORBT,1.D0,
     &           H1X,NORBT)
C
C     End of ONEXH1
C
      RETURN
      END
C  /* Deck slvsor */
      SUBROUTINE SLVSOR(ONEIND,SNGDN,NSIM,UDV,EVECS,ZYMAT)
C
C WRITTEN 14-april 1993 (modified RSPSOR)
C
C PURPOSE:
C
C   CREATE ORBITAL PART OF LINEAR TRANSFORMATION WITH
C   solvent contribution to E(2)
C
C    ONEIND = .TRUE. FOR NORMAL ONE ELECTRON DENSITY MATRIX
C
C    ONEIND = .FALSE. FOR TRANSITION DENSITY MATRICES WHERE THE
C                     INACTIVE INACTIVE PART OF THE ONE ELECTRON
C                     DENSITY MATRIX IS ZERO
C                ******************************
C
C    SNGDN  = .TRUE.  FOR SINGLET DENSITIES
C
C             .FALSE. FOR TRIPLET DENSITIES WHERE INACTIVE
C                     CONTRIBUTIONS ARE ZERO
C
C    [ E(P,Q) , ZYM ] =
C
C    ACTIVE-INACTIVE (P,Q) = (T,I)
C                             L,K
C    ZYM(I,X)*DV(T,X)-2*ZYM(I,T)
C        K       L          K,L
C
C    INACTIVE-ACTIVE         (I,U)
C                             K,L
C    -ZYM(X,I)*DV(X,U)+2*ZYM(U,I)
C           K       L        L,K
C
C    ACTIVE-ACTIVE           (T,U)
C                             K,L
C                             L,K
C    ZYM(U,X)*DV(T,X)-ZYM(X,T)*DV(X,U)
C        L       K          K       L
C        K       L          L       K
C
C    ACTIVE-SECONDARY        (T,A)
C                             K,L
C    ZYM(A,X)*DV(T,X)
C        L       K
C
C    SECONDARY-ACTIVE        (A,U)
C                             L,K
C    -ZYM(X,A)*DV(X,U)
C           L       K
C
C    INACTIVE-SECONDARY      (I,A)
C                             K,L
C    2*ZYM(A,I)
C          L,K
C
C    SECONDARY-INACTIVE      (A,J)
C                             L,K
C    -2*ZYM(J,A)
C           K,L
C
C ZYM IS ONE ELECTRON OPERATOR WITH ZYMAT AS ELEMENTS
C
#include "implicit.h"
C
      LOGICAL ONEIND,SNGDN
C
      DIMENSION UDV(NASHDI,NASHDI,*), EVECS(KZYVAR,*)
      DIMENSION ZYMAT(NORBT,NORBT,*)
C
C  INFDIM : NASHDI
C
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "infvar.h"
#include "inforb.h"
#include "infind.h"
#include "infdim.h"
#include "infpri.h"
#include "infrsp.h"
#include "wrkrsp.h"
C
C -- local constants
C
      PARAMETER (  D2 =2.0D0  )
C
C
      KYCONF = KZCONF + KZVAR
C
C DISTRIBUTE FOCK MATRICES IN EVECS
C
         KSYM1 = 0
         DO 1300 IG = 1,KZWOPT
            K     = JWOP(1,IG)
            L     = JWOP(2,IG)
            KSYM  = ISMO(K)
            LSYM  = ISMO(L)
            IF( KSYM.NE.KSYM1 ) THEN
               KSYM1 = KSYM
               IORBK = IORB(KSYM)
               NASHK = NASH(KSYM)
               NISHK = NISH(KSYM)
               IASHK = IASH(KSYM)
               IIORBK= IIORB(KSYM)
               IORBL = IORB(LSYM)
               NASHL = NASH(LSYM)
               NISHL = NISH(LSYM)
               IASHL = IASH(LSYM)
            END IF
            ITYPK = IOBTYP(K)
            ITYPL = IOBTYP(L)
            IF ( ITYPK.EQ.JTINAC ) THEN ! inactive-active or inactive-secondary
               DO 2000 ISIM = 1 ,NSIM
                  IF (ONEIND) THEN
                   IF ( SNGDN ) THEN
                     EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                       + D2 * ZYMAT(L,K,ISIM)
                     IF ((.NOT. TDA) .AND. (.NOT. CISRPA)) THEN
                       EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                         - D2 * ZYMAT(K,L,ISIM)
                     END IF
                   END IF
                  END IF
                  IF ( ITYPL.EQ.JTACT ) THEN ! inactive-active
                   IF (ONEIND) THEN
                     NWL = ISW(L) - NISHT
                     EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                  - DDOT(NASHL,ZYMAT(IORBL+NISHL+1,K,ISIM),1,
     *                                        UDV(IASHL+1,NWL,1),1)
                     IF ((.NOT. TDA) .AND. (.NOT. CISRPA)) THEN
                        DO 735 IX = 1,NASHL
                           EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                        + ZYMAT(K,IORBL+NISHL+IX,ISIM)*
     *                          UDV(NWL,IASHL+IX,1)
 735                    CONTINUE
                     ENDIF
                   ELSE
                     NWL = ISW(L) - NISHT
                     EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                  - DDOT(NASHK,ZYMAT(IORBK+NISHK+1,K,ISIM),1,
     *                                        UDV(IASHK+1,NWL,1),1)
                     IF ((.NOT. TDA) .AND. (.NOT. CISRPA)) THEN
                        DO 736 IX = 1,NASHK
                           EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                        + ZYMAT(K,IORBK+NISHK+IX,ISIM)*
     *                          UDV(NWL,IASHK+IX,1)
 736                    CONTINUE
                     END IF
                   END IF
                  END IF
 2000          CONTINUE
            ELSE
              IF (ITYPL.EQ.JTACT) THEN ! active-active
                IF (ONEIND) THEN
                  NWL = ISW(L) - NISHT
                  NWK = ISW(K) - NISHT
                  DO 2020 ISIM=1,NSIM
                     EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                     -DDOT(NASHL,ZYMAT(IORBL+NISHL+1,K,ISIM),1,
     *                                        UDV(IASHL+1,NWL,1),1)
                     IF ((.NOT. TDA) .AND. (.NOT. CISRPA)) THEN
                        DO 745 IX = 1,NASHL
                           EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                        +ZYMAT(K,IORBL+NISHL+IX,ISIM)*
     *                                        UDV(NWL,IASHL+IX,1)
 745                    CONTINUE
                        EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                     -DDOT(NASHK,ZYMAT(IORBK+NISHK+1,L,ISIM),1,
     *                                      UDV(IASHK+1,NWK,1),1)
                     END IF
                     DO 755 IX = 1,NASHK
                        EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                        +ZYMAT(L,IORBK+NISHK+IX,ISIM)*
     *                                      UDV(NWK,IASHK+IX,1)
 755                 CONTINUE
 2020             CONTINUE
                ELSE
                  NWL = ISW(L) - NISHT
                  NWK = ISW(K) - NISHT
                  DO 2021 ISIM=1,NSIM
                     EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                     -DDOT(NASHK,ZYMAT(IORBK+NISHK+1,K,ISIM),1,
     *                                        UDV(IASHK+1,NWL,1),1)
                     IF ((.NOT. TDA) .AND. (.NOT. CISRPA)) THEN
                        DO 746 IX = 1,NASHK
                           EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                        +ZYMAT(K,IORBK+NISHK+IX,ISIM)*
     *                                        UDV(NWL,IASHK+IX,1)
 746                    CONTINUE
                        EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                     -DDOT(NASHL,ZYMAT(IORBL+NISHL+1,L,ISIM),1,
     *                                      UDV(IASHL+1,NWK,1),1)
                     END IF
                     DO 756 IX = 1,NASHL
                        EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                        +ZYMAT(L,IORBL+NISHL+IX,ISIM)*
     *                                      UDV(NWK,IASHL+IX,1)
 756                 CONTINUE
 2021             CONTINUE
                END IF
             ELSE ! active-secondary
               IF (ONEIND) THEN
                  NWK = ISW(K) - NISHT
                  DO 2030 ISIM=1,NSIM
                     IF ((.NOT. TDA) .AND. (.NOT. CISRPA)) THEN
                        EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                     -DDOT(NASHK,ZYMAT(IORBK+NISHK+1,L,ISIM),1,
     *                                        UDV(IASHK+1,NWK,1),1)
                     END IF
                     DO 765 IX = 1,NASHK
                        EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                        +ZYMAT(L,IORBK+NISHK+IX,ISIM)*
     *                                        UDV(NWK,IASHK+IX,1)
 765                 CONTINUE
 2030             CONTINUE
              ELSE
                  NWK = ISW(K) - NISHT
                  DO 2031 ISIM=1,NSIM
                     IF ((.NOT. TDA) .AND. (.NOT. CISRPA)) THEN
                        EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                     -DDOT(NASHL,ZYMAT(IORBL+NISHL+1,L,ISIM),1,
     *                                        UDV(IASHL+1,NWK,1),1)
                     END IF
                     DO 766 IX = 1,NASHL
                        EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                        +ZYMAT(L,IORBL+NISHL+IX,ISIM)*
     *                                        UDV(NWK,IASHL+IX,1)
 766                 CONTINUE
 2031             CONTINUE
              END IF
             ENDIF
            ENDIF
 1300    CONTINUE
C
      IF (IPRRSP.GT.110) THEN
         WRITE(LUPRI,*)' ONEIND =',ONEIND
         WRITE(LUPRI,*)
     &   ' After LINEAR TRANSFORMATION with solvent part of E[2]'
         CALL OUTPUT(EVECS,1,KZYVAR,1,NSIM,KZYVAR,NSIM,1,LUPRI)
      END IF
C
C END OF SLVSOR
C
      RETURN
      END
