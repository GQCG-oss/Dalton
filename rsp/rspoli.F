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
#ifdef UNDEF
===========================================================================
/* Comdeck rev_log */
Revision 1.2  2000/05/24 19:04:56  hjj
new GETREF call; LUINTD .gt. 0 check for GPCLOSE call

990427-hjaaj: use if DXAO is symmetric or antisymmetric in IFCTYP
950513-kr Added comdeck INFVAR and and changed to ISYMDM()=JWOPSY
961228-kr Transfered SOPPA particle-hole matrices in paramater list to RSPOLI
9412230hjaaj: revised SIRFCK call
940708-hjaaj: SOPPA changes
940602-hjaaj: do not call QONEDI and QTD for nasht.eq.1
931124-martin packer
RSPOLI: added call to HRPA which adds A(2) and B(2) matrices to FCX
920721-hinne hettema
RSPOLI: implemented RSPSUP (super symmetry averaging)
===========================================================================
#endif
C  /* Deck rspoli */
      SUBROUTINE RSPOLI(NCSIM,NOSIM,UDV,ZYCVEC,LZYCVEC,FC,FV,PVX,ZYMAT,
     *                  FCX,FVX,QAX,QBX, FVTD,QATD,QBTD, EVECS,
     *                  XINDX,CMO,WRK,LWRK)
C
C Copyright 13 FEB 1986 Poul Joergensen and Hans Joergen Aa. Jensen
C
C 21-1-1992: Changed for Direct RPA (KEYWORD: DIRFCK), HA.
C 21-7-1992: Supersymmetry averaging added. Hinne Hettema
C 28-2-1995: Changed for Direct RPA with AO integrals read
C            from disk, NOITRA=.TRUE. P&D /is not included?/hjaaj
C Oct-2007-hjaaj: added LZYCVEC because in lagrang.F RSPOLI is called
C            with ZYCVEC=CREF, and for some RHF cases KZCONF=0 and NCREF=1
C
C DFT modifications T. Helgaker 
C
C PURPOSE:
C    CALCULATE ORBITAL PART OF THE LINEAR TRANSFORMATIONS
C      E[2]*X
C
C FLOW:
C    1) READ IN MULLIKEN INTEGRAL DISTRIBUTION AND ADD CONTRIBUTIONS
C       TO FOCK CORE , FOCK VALENCE , QA AND QB MATRICES
C    2) READ IN DIRAC INTEGRAL DISTRIBUTION AND ADD CONTRIBUTIONS
C       TO QAX AND QBX MATRICES
C    3) DISTRIBUTE FOCK CORE , FOCK VALENCE , QA AND QB MATRICES
C       INTO ORBITAL PART OF LINEAR TRANSFORMED VECTORS
C
C Naming conventions:
C    FCX, FVX, ... are FC, FV, ... with one-index transformed integrals
C    FVTD, ...      are FV calculated with Transition dens.mat.
C
#include <implicit.h>
C
      DIMENSION UDV(NASHDI,*),ZYCVEC(*),FC(*),FV(*),PVX(*)
      DIMENSION ZYMAT(N2ORBX,*),FCX(NORBT,NORBT,*)
      DIMENSION FVX(NORBT,NORBT,*),QAX(NORBT,NASHDI,*)
      DIMENSION QBX(NORBT,NASHDI,*),FVTD(NORBT,NORBT,*)
      DIMENSION QATD(NORBT,NASHDI,*),QBTD(NORBT,NASHDI,*)
      DIMENSION EVECS(*)
      DIMENSION XINDX(*),WRK(*),CMO(*)
C
#include <thrzer.h>
#include <dummy.h>
C
C  INFDIM : NASHDI
C  INFINP : DIRFCK
C  INFRSP : ?
C  INFSOP : A2EXIST,KABSAD,KABTAD,...
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infinp.h>
#include <infdim.h>
#include <infpri.h>
#include <infrsp.h>
#include <infvar.h>
#include <inftra.h>
#include <orbtypdef.h>
#include <wrkrsp.h>
#include <infsop.h>
#include <dftcom.h>
CSONIA SONIA SONIA
CSONIA SONIA SONIA control direct for kappabar of CCSD(T)
CSONIA SONIA SONIA
#include <grdccpt.h>
CSONIA SONIA SONIA
CSONIA SONIA SONIA end control direct for kappabar of CCSD(T)
CSONIA SONIA SONIA
C
C -- local constants
C
      PARAMETER ( D1 = 1D0, DP5 = 0.5D0, D2 = 2.0D0, DM1 = -1.0D0 )
      PARAMETER ( MXFMAT = 100 )
      LOGICAL   DFTADX, ROKS
      LOGICAL   DOFXC, DOFXV, DOFTC, DOFTV
      DIMENSION ISYMDM(MXFMAT), IFCTYP(MXFMAT)
      DIMENSION NEEDMU(6),NEEDDI(6)
      DATA NEEDMU /1,1,1,0,0,0/
      DATA NEEDDI /1,1,1,0,0,0/
      CALL QENTER('RSPOLI')
C
      ROKS = NASHT.GT.0 .AND. DODFT
C
      IF (SOPPA) THEN
         IF (DIRFCK) THEN
           CALL QUIT('RSPOLI: SOPPA not implemented for direct Fock!')
         ENDIF
         NEEDMU(4) = 1
C        exclude virt-virt contributions to B(2) matrix
C        for .OPTORB iterations from ORPCTL
         IF (NOSIM .GT. 0 .AND. KZCONF .GT. 0) NEEDMU(6) = 1
      ENDIF

      DOFXC  = NISHT .GT. 0 .OR. DODFT .OR. DOHFSRDFT .OR. DOMCSRDFT
      DOFXV  = NASHT .GT. 0
      DOFTC  = DOMCSRDFT
      DOFTV  = NASHT .GT. 0

      IF (DOMCSRDFT) THEN
         NEEDMU(1) = 0
         NEEDMU(2) = 0
         NEEDMU(3) = 1
         NEEDDI(1) = 0
         NEEDDI(2) = 0
         NEEDDI(3) = 1
      END IF

C
C ALLOCATE WORK SPACE FOR DTV and PTVD and SOPPA
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWRK
      CALL MEMGET('REAL',KDTV ,NCSIM*N2ASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KPTVD,NCSIM*N2ASHX*N2ASHX,WRK,KFREE,LFREE)
      IF (KZCONF.GT.0) THEN
         LH2XAC = NOSIM*N2ASHX*NNASHX
      ELSE
         LH2XAC = 0
      END IF
      CALL MEMGET('REAL',KH2XAC,LH2XAC,WRK,KFREE,LFREE)
      IF (SOPPA) THEN
         CALL MEMGET('REAL',KH2XP ,N2ORBX,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KCOEUN,N2ORBX,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KTZYMT,NOSIM*N2ORBX,WRK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KH2XP ,0,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KCOEUN,0,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KTZYMT,0,WRK,KFREE,LFREE)
      ENDIF
      KWRK0 = KFREE
C
C CALCULATE TRANSITION DENSITY MATRIX
C
      IF ( (NCSIM.GT.0) .AND. (.NOT.RSPCI) .AND. (.NOT.SOPPA) ) THEN
         CALL MEMGET('REAL',KCREF,NCREF,WRK,KFREE,LFREE)
         IF (DOMCSRDFT) KWRK0 = KFREE ! must keep CREF for later
         CALL GETREF(WRK(KCREF),NCREF)
         IF (TRPLET) THEN
            ISPIN1 = 1
            ISPIN2 = 0
         ELSE
            ISPIN1 = 0
            ISPIN2 = 0
         END IF
         CALL RSPTDM(NCSIM,IREFSY,KSYMST,NCREF,LZYCVEC,WRK(KCREF),
     *                 ZYCVEC,WRK(KDTV),WRK(KPTVD),
     *                 ISPIN1,ISPIN2,.TRUE.,.FALSE.,
     *                 XINDX,WRK,KFREE,LFREE)
C        CALL RSPTDM(NCSIM,ILRESY,IRSYM,NCLREF,NCRDIM,CLREF,
C    *                 CR, RHO1,RHO2, ISPIN1,ISPIN2,TDM,NORHO2,
C    *                 XNDXCI,WORK,KFREE,LFREE)
         CALL MEMREL('RSPOLI.tdm',WRK,KFRSAV,KWRK0,KFREE,LFREE)
      END IF

C
C INITIALIZE ONE INDEX TRANSFORMED (X) AND TRANSITION DENSITY (TD)
C Q, FC AND FV MATRICES
C
      IF (NOSIM.GT.0) THEN
         IF (SOPPA) THEN
C        ... A(2) matrix is saved permanently in XINDX
C        Initialisation of XINDX is done in SET2SOPPA now.
C            IF (.NOT.A2EXIST) CALL DZERO(XINDX,N2ORBX)
            DO I = 1,NOSIM
              JTZYMT = KTZYMT + (I-1)*N2ORBX
              CALL MTRSP(NORBT,NORBT,ZYMAT(1,I),NORBT,WRK(JTZYMT),NORBT)
C             CALL MTRSP(NROWA,NCOLA,A,NRDIMA,B,NRDIMB)
            END DO
         END IF
         CALL DZERO(FCX,N2ORBX*NOSIM)
         IF (NASHT.GT.0) THEN
            CALL DZERO(FVX,N2ORBX*NOSIM)
            NTOT = NORBT*NASHT*NOSIM
            CALL DZERO(QAX,NTOT)
            CALL DZERO(QBX,NTOT)
            IF (KZCONF.GT.0 ) THEN
               NTOT = NOSIM*N2ASHX*NNASHX
               CALL DZERO(WRK(KH2XAC),NTOT)
            END IF
         END IF
      END IF
      IF (NCSIM.GT.0) THEN
         NTOT =  N2ORBX * NCSIM
         IF (DOMCSRDFT) NTOT = 2*NTOT
C        ... we also need FCTD for MC-srDFT
         CALL DZERO(FVTD,NTOT)
         NTOT = NORBT*NASHT*NCSIM
         IF (NTOT .GT. 0) THEN
C        ... is zero for SOPPA
            CALL DZERO(QATD,NTOT)
            CALL DZERO(QBTD,NTOT)
         END IF
      END IF
C
C If DIRFCK compute Fock matrices in AO basis
C If also single determinant, then nothing to do in MO part
C
      IF (DIRFCK .AND. (TDHF .OR. DODFT)) GO TO 1000
      IF (DODFT) GO TO 1000
C     ... KS-DFT only works with direct code, even if DIRFCK is false /hjaaj

CSONIA SONIA SONIA
CSONIA SONIA SONIA control direct for kappabar of CCSD(T)
CSONIA SONIA SONIA
      IF (LGRDCCPT) THEN
         write(lupri,*)'Warning RSPOLI: LGRDCCPT = ', LGRDCCPT
         GO TO 1000
      END IF
CSONIA SONIA SONIA
CSONIA SONIA SONIA end control direct for kappabar of CCSD(T)
CSONIA SONIA SONIA
C
C If not DIRFCK, i.e.
C calculate Fock matrices etc. in MO basis, then
C ALLOCATE WORK SPACE FOR MULLIKEN and DIRAC DISTRIBUTIONs
C to Fock matrices
C
      IF (.NOT.DIRFCK) THEN
         CALL MEMGET('REAL',KDENA ,NOSIM*NORBT*NASHT,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KDENB ,NOSIM*NORBT*NASHT,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KFCOCO,NOSIM*N2ORBX,WRK,KFREE,LFREE)
         CALL DZERO(WRK(KFCOCO),N2ORBX*NOSIM)
         IF (NASHT.GT.0) THEN
            CALL MEMGET('REAL',KFVOCO,NOSIM*N2ORBX,WRK,KFREE,LFREE)
            CALL DZERO(WRK(KFVOCO),N2ORBX*NOSIM)
         ELSE
            CALL MEMGET('REAL',KFVOCO,0,WRK,KFREE,LFREE)
         END IF
         IF ((NASHT .GT. 0) .AND. (NOSIM .GT.0)) THEN
            CALL RSPTR1(NOSIM,UDV,ZYMAT,WRK(KDENA),WRK(KDENB))
         END IF
      ELSE
         CALL MEMGET('REAL',KDENA ,0,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KDENB ,0,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KFCOCO,0,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KFVOCO,0,WRK,KFREE,LFREE)
      END IF
      CALL MEMGET('REAL',KH2   ,      N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2X  ,      N2ORBX,WRK,KFREE,LFREE)
      KWRK1 = KFREE
C
C setup for reading Mulliken MO integrals ...
C NEEDED DISTRIBUTIONS DEFINED IN NEEDMU(6)
C
      CALL MEMGET('REAL',KPVCD ,      N2ASHX,WRK,KFREE,LFREE)
      IDIST = 0
 90   CALL NXTH2M(IC,ID,WRK(KH2),NEEDMU,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 95
C     ... IF IDIST.LT.0 NO MORE DISTRIBUTIONS
         KWRK2 = KFREE
         LWRK2 = LFREE
C        (KFREE,LFREE are updated inside NXTH2M)
      ICDSYM= MULD2H(ISMO(IC),ISMO(ID))
      ICW   = ISW(IC)
      IDW   = ISW(ID)
      IDI   = ID
      ICI   = IC
      ICIW  = ICW
      IDIW  = IDW
C
C     find distribution type ITYPCD =
C     1:inactive-inactive  2:active-inactive  3:active-active
C     4:secondary-inactive 5:secondary-active 6:secondary-secondary
C     We do not need type 6 (except for SOPPA).
C
      ITYPC  = IOBTYP(IC)
      ITYPD  = IOBTYP(ID)
      ITYPCD = IDBTYP(ITYPC,ITYPD)
C
C ORDER INDICES OF C AND D SUCH THAT C < D
C
      IF ( ICW.GT.IDW ) THEN
         ICIW   = IDW
         IDIW   = ICW
         ICI    = ID
         IDI    = IC
         ISWAP  = ITYPC
         ITYPC  = ITYPD
         ITYPD  = ISWAP
      ENDIF
      NCIW  = ICIW - NISHT
      NDIW  = IDIW - NISHT
C

      IF (NOSIM.GT.0) THEN
C
C CONTRIBUTIONS TO FCX AND FVX
C
         IF (.NOT.DIRFCK .AND. ITYPCD .LE. 3)
     &      CALL FONEMU(NOSIM,ICI,IDI,WRK(KH2),
     &               FCX,FVX,ZYMAT,WRK(KDENA),WRK(KDENB),
     &               WRK(KWRK2),LWRK2)
C           CALL FONEMU(NSIM,ICI1,IDI1,H2,
C    *                  FCOEX,FVOEX,ZYMAT,DENA,DENB,WRK,LWRK)
C
Cmjp CONTRIBTUIONS FORM THE HRPA TERMS
         IF (SOPPA) THEN
            CALL HRPA(FCX,UDV,PVX,WRK(KH2),ZYMAT,WRK(KTZYMT),
     *                WRK(KH2X),WRK(KH2XP),WRK(KCOEUN),ICI,IDI,ICDSYM,
     *                ITYPCD,NOSIM,XINDX(KIADR1),WRK(KWRK2),LWRK2)
            IF (ITYPCD .EQ. 4) THEN
               IF (.NOT.HIRPA .AND. KZCONF .GT. 0) THEN
C                 ... skip SOPH2X for .OPTORB iterations
                  CALL SOPH2X(EVECS(1+NCSIM*KZYVAR),ZYMAT,WRK(KTZYMT),
     *                      WRK(KH2),ICI,IDI,ICDSYM,NOSIM,XINDX(KABSAD),
     &                      XINDX(KABTAD),XINDX(KIJSAD),XINDX(KIJTAD),
     &                      XINDX(KIJ1AD),XINDX(KIJ2AD),XINDX(KIJ3AD),
     *                      WRK(KWRK2),LWRK2)
               END IF
C  Construct the A(2) matrix explicitly (but only once)
C  for the symmetrization of A(2)b
               IF (.NOT.A2EXIST) THEN
                  CALL HRPAA2(PVX,WRK(KH2),XINDX(KAB2),WRK(KCOEUN),
     &                        XINDX(KIADR1),ICI,IDI,ICDSYM)
               END IF
            END IF
         END IF
         IF ((NASHT.GT.1).AND.(ITYPCD.EQ.3)) THEN
C
C CONTRIBUTIONS TO QAX AND QBX
C
            CALL QONEMU(NOSIM,NCIW,NDIW,ICDSYM,
     *                  QAX,QBX,ZYMAT,WRK(KH2),WRK(KH2X),WRK(KH2XAC),
     *                  PVX,WRK(KPVCD),WRK(KWRK2),LWRK2)
C           CALL QONEMU(NOSIM,NCIW,NDIW,ICDSYM,QAX,QBX,
C    *                  ZYMAT,H2,H2X,H2XAC,PVX,PVCD,WRK,LWRK)
         END IF
      END IF
      IF (NCSIM.GT.0) THEN
         IF (SOPPA) THEN
            IF (.NOT.HIRPA .AND. ITYPCD.EQ.4) THEN
               CALL SOPPAF(FVTD,WRK(KH2),ZYCVEC,ICI,IDI,ICDSYM,NCSIM,
     *                     XINDX,WRK(KWRK2),LWRK2)
            END IF
         ELSE
            IF (.NOT. DIRFCK .AND.
     &         (ITYPCD.EQ.2 .OR. ITYPCD.EQ.3 .OR. ITYPCD.EQ.5) ) THEN
C
C CONTRIBUTIONS TO FVTD
C
               CALL FTDMU(NCSIM,ICI,IDI,FVTD,WRK(KDTV),WRK(KH2))
C              CALL FTDMU(NSIM,ICI1,IDI1,FVTD,DTV,WRK(KH2))
            END IF
            IF (ITYPCD.EQ.3) THEN
C
C CONTRIBUTIONS TO QATD AND QBTD
C
              CALL QTD(NCSIM,NCIW,NDIW,ICDSYM,QATD,QBTD,WRK(KH2),
     *                 WRK(KPTVD),WRK(KPVCD),WRK(KWRK2),LWRK2)
C             CALL QTD(NCSIM,NCIW,NDIW,ICDSYM,QATD,QBTD,H2,PTVD,PVDEN,
C    *                  WRK,LWRK)
            END IF
         END IF
      END IF
      GO TO 90
C
C ALL MULLIKEN DISTRIBUTIONS HAVE BEEN READ IN
C
 95   CONTINUE
      IF ((IPRRSP.GT.50).AND.(NOSIM.GT.0).AND.(NASHT.GT.1)) THEN
         DO 1001 IOSIM = 1,NOSIM
         WRITE(LUPRI,'(/A)') ' MULLIKEN DISTRIBUTION CONTRIBUTION'
         WRITE(LUPRI,'(/A,I5,A)')
     *      ' QAX FOR',IOSIM,' ORBITAL TRIAL VECTOR'
         CALL OUTPUT(QAX(1,1,IOSIM),
     *                     1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
         WRITE(LUPRI,'(/A,I5,A)')
     *      ' QBX FOR',IOSIM,' ORBITAL TRIAL VECTOR'
         CALL OUTPUT(QBX(1,1,IOSIM),
     *                     1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
 1001    CONTINUE
      END IF
      IF ((IPRRSP.GT.50).AND.(NOSIM.GT.0).AND. .NOT.DIRFCK) THEN
         DO 1002 IOSIM = 1,NOSIM
         WRITE(LUPRI,'(/2A,I5)')' FCOEX MULLIKEN DISTRIBUTION CONTRB'
     *      ,' FOR VECTOR ',IOSIM
         CALL OUTPUT(FCX(1,1,IOSIM),
     *                     1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(/2A,I5)')' FVOEX MULLIKEN DISTRIBUTION CONTRB'
     *      ,' FOR VECTOR ',IOSIM
         CALL OUTPUT(FVX(1,1,IOSIM),
     *                     1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
 1002    CONTINUE
      END IF
      IF ((IPRRSP.GT.50).AND.(NCSIM.GT.0).AND.(NASHT.GT.1)) THEN
         DO 1003 ICSIM = 1,NCSIM
         WRITE(LUPRI,'(/A)') ' MULLIKEN DISTRIBUTION CONTRIBUTION'
         WRITE(LUPRI,'(/A,I5,A)')
     *      ' QATD FOR',ICSIM,' CONFIGURATION TRIAL VECTOR'
         CALL OUTPUT(QATD(1,1,ICSIM),
     *                     1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
         WRITE(LUPRI,'(/A,I5,A)')
     *      ' QBTD FOR',ICSIM,' CONFIGURATION TRIAL VECTOR'
         CALL OUTPUT(QBTD(1,1,ICSIM),
     *                     1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
 1003    CONTINUE
      END IF
      IF ((IPRRSP.GT.50).AND.(NCSIM.GT.0).AND. .NOT.DIRFCK) THEN
         DO 1004 ICSIM = 1,NCSIM
         WRITE(LUPRI,'(/2A,I5)')' FVTD MULLIKEN DISTRIBUTION CONTRB'
     *      ,' FOR VECTOR ',ICSIM
         CALL OUTPUT(FVTD(1,1,ICSIM),
     *                     1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
 1004    CONTINUE
      END IF
      CALL MEMREL('RSPOLI.h2m',WRK,KFRSAV,KWRK1,KFREE,LFREE)
C
C   **********************************************************
C
C ADD CONTRIBUTIONS  FROM DIRAC INTEGRAL DISTRIBUTIONS
C
C   *********************************************************
C
C ALLOCATE WORK SPACE FOR DIRAC DISTRIBUTION
C
      CALL MEMGET('REAL',KPVCD2,N2ASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KPVCD3,N2ASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KPVDC2,N2ASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KPVDC3,N2ASHX,WRK,KFREE,LFREE)
C
      IDIST  = 0
      LUINTD = -1
      IF (.NOT. NEWTRA) CALL GPOPEN(LUINTD,'DRCINT','OLD',' ',
     &                              'UNFORMATTED',IDUMMY,.FALSE.)
 1090 CALL NXTH2D(IC,ID,WRK(KH2),NEEDDI,LUINTD,WRK,KFREE,LFREE,IDIST)
         IF (IDIST .LT. 0) GO TO 1095
C  IF IDIST.LT.0 NO MORE DISTRIBUTIONS
         ICDSYM= MULD2H(ISMO(IC),ISMO(ID))
         ICW   = ISW(IC)
         IDW   = ISW(ID)
         IDI   = ID
         ICI   = IC
         ICIW  = ICW
         IDIW  = IDW
C
C ORDER INDICES OF C AND D SUCH THAT C =>D
C
         IF ( IDW.GT.ICW ) THEN
            ICIW   = IDW
            IDIW   = ICW
            ICI    = ID
            IDI    = IC
            CALL DGETRN(WRK(KH2),NORBT,NORBT)
         ENDIF
         NCIW  = ICIW - NISHT
         NDIW  = IDIW - NISHT
C
         IF ( IPRRSP.GT.150 ) THEN
            WRITE(LUPRI,'(/A,2I8)')' DIRAC DISTRIBUTION : IC,ID ',IC,ID
            CALL OUTPUT(WRK(KH2),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            write(lupri,*) 'NOSIM=', NOSIM
         ENDIF
         IF (NOSIM.GT.0) THEN
            IF (.NOT.DIRFCK)
     &         CALL FONEDR(NOSIM,ICI,IDI,WRK(KH2),
     &                  WRK(KFCOCO),WRK(KFVOCO),FCX,FVX,ZYMAT,
     &                  WRK(KDENA),WRK(KDENB),WRK(KFREE),LFREE)
C           CALL FONEDR(NSIM,ICI1,IDI1,H2D,FCOCO,FVOCO,
C    *                     FCOEX,FVOEX,ZYMAT,DENA,DENB,WRK,LWRK)
C
C CONTRIBUTIONS TO QAX AND QBX
C
            IF ( (IOBTYP(ICI).EQ.JTACT) .AND. (IOBTYP(IDI).EQ.JTACT))
     &         CALL QONEDI(NOSIM,NCIW,NDIW,ICDSYM,QAX,QBX,ZYMAT,
     &            WRK(KH2),WRK(KH2X),PVX,WRK(KPVCD2),WRK(KPVCD3),
     &            WRK(KPVDC2),WRK(KPVDC3),WRK(KFREE),LFREE)
C           CALL QONEDI(NOSIM,NCIW,NDIW,ICDSYM,
C    *               QAX,QBX,ZYMAT,H2,H2X,
C    *               PVX,PVCD2,PVCD3,PVDC2,PVDC3,WRK,LWRK)
         END IF
#ifdef DEBUG_RSPOLI

         IF (NOSIM.GT.0) THEN
            DO IOSIM = 1,NOSIM
            WRITE(LUPRI,'(/A)') ' DIRAC DISTRIBUTION CONTRIBUTION added'
            WRITE(LUPRI,'(/A,I5,A)')
     *      ' QAX FOR',IOSIM,' ORBITAL TRIAL VECTOR'
            CALL OUTPUT(QAX(1,1,IOSIM),
     *                  1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
            WRITE(LUPRI,'(/A,I5,A)')
     *      ' QBX FOR',IOSIM,' ORBITAL TRIAL VECTOR'
            CALL OUTPUT(QBX(1,1,IOSIM),
     *                     1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
            END DO
         END IF

#endif
         IF (NCSIM.GT.0 .AND. .NOT.SOPPA .AND. .NOT.DIRFCK) THEN
            CALL FTDDR(NCSIM,ICI,IDI,FVTD,WRK(KDTV),WRK(KH2))
C           CALL FTDDR(NSIM,ICI1,IDI1,FVTD,DTV,H2D)
         END IF
C
         GO TO 1090
C
C ALL DIRAC DISTRIBUTIONS HAVE BEEN READ IN
C
 1095 CONTINUE
      IF (LUINTD .GT. 0) CALL GPCLOSE(LUINTD,'KEEP')
C
      IF (DIRFCK) THEN
         CALL MEMREL('RSPOLI.h2m+d',WRK,KFRSAV,KWRK0,KFREE,LFREE)
         GO TO 1000
C we now go to the direct part of rspoli in order to build
C the Fock operators
      ENDIF
C
C ADD CONTRIBUTION TO FCX AND FVX FROM ONE INDEX TRANSFORMED
C TOTAL SYMMETRIC FOCK MATRICES
C
      IF ( NOSIM .GT. 0 )  THEN
C
C SUM UP COULOMN CONTRIBUTIONS WHICH ARE OBTAINED USING THE
C TERM IS SYMMETRIC
C
         IF (IPRRSP.GT.50) THEN
            DO 1010 IOSIM = 1,NOSIM
               WRITE(LUPRI,'(/2A,I5)')' FCOEX  CONTRB'
     *                              ,' FOR VECTOR ',IOSIM
               CALL OUTPUT(FCX(1,1,IOSIM),
     *                     1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
               WRITE(LUPRI,'(/2A,I5)')' FCOCO CONTRB'
     *                              ,' FOR VECTOR ',IOSIM
               CALL OUTPUT(WRK(KFCOCO+ (IOSIM-1)*NORBT*NORBT),
     *                     1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
               IF (NASHT.GT.0) THEN
                  WRITE(LUPRI,'(/2A,I5)')' FVOEX  CONTRB'
     *                                 ,' FOR VECTOR ',IOSIM
                  CALL OUTPUT(FVX(1,1,IOSIM),
     *                        1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                  WRITE(LUPRI,'(/2A,I5)')' FVOCO CONTRB'
     *                                 ,' FOR VECTOR ',IOSIM
                  CALL OUTPUT(WRK(KFVOCO+ (IOSIM-1)*NORBT*NORBT),
     *                        1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
               END IF
 1010       CONTINUE
         END IF
         CALL FONEDN(NOSIM,WRK(KFCOCO),WRK(KFVOCO),FCX,FVX)
      ENDIF
C
      CALL MEMREL('RSPOLI.h2m+d',WRK,KFRSAV,KWRK0,KFREE,LFREE)
      GO TO 2000
 1000 CONTINUE
C
C end up here if DIRFCK
C       
C       Allocate work memory:
C       FXCAO,FXVAO,FTCAO,FTVAO must be stored consecutively
C       DXCAO,DXVAO,DTCAO,DTVAO must be stored consecutively
C       (requirement for SIRFCK call)
C
      NFMAT  = 0
      IF (DOFXC) NFMAT = NFMAT + NOSIM
C     ... for FXCAO (+ VXxcAO when DFT or srDFT)
      IF (DOFXV) NFMAT = NFMAT + NOSIM
C     ... for FXVAO
      IF (DOFTC) NFMAT = NFMAT + NCSIM
C     ... for FTCAO (with V[2c]xcAO when srDFT)
      IF (DOFTV) NFMAT = NFMAT + NCSIM
C     ... for FTVAO
C
      NDMAT = NFMAT
      IF (DOHFSRDFT .OR. DOMCSRDFT) NDMAT = NDMAT + 2
C     ... for DtotAO, DXtotAO for SRDFT:SRDFTLTR
C
      CALL MEMGET('REAL',KFXCAO,NFMAT*N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDXCAO,NDMAT*N2BASX,WRK,KFREE,LFREE)
      KFXVAO = KFXCAO + NOSIM*N2BASX
      KDXVAO = KDXCAO + NOSIM*N2BASX
      IF (DOFXV) THEN
         KFTCAO  = KFXVAO + NOSIM*N2BASX
         KDTCAO  = KDXVAO + NOSIM*N2BASX
      ELSE
         KFTCAO  = KFXVAO
         KDTCAO  = KDXVAO
      END IF
      IF (DOFTC) THEN
         KFTVAO  = KFTCAO + NCSIM*N2BASX
         KDTVAO  = KDTCAO + NCSIM*N2BASX
      ELSE
         KFTVAO  = KFTCAO
         KDTVAO  = KDTCAO
      END IF
      KDTOTAO = KDTVAO + NCSIM*N2BASX

#ifdef DEBUG_RSPOLI
      write(lupri,*) 'rspoli: ncsim,nosim',ncsim,nosim
      write(lupri,*) 'rspoli: nfmat,ndmat',nfmat,ndmat
      write(lupri,*) 'rspoli: KDXCAO,KDXVAO,KDTCAO,KDTVAO,
     &                        KDTOTAO,kfree',
     &                KDXCAO,KDXVAO,KDTCAO,KDTVAO,KDTOTAO,KFREE
#endif

#ifdef MOD_SRDFT
      IF (DOHFSRDFT .OR. DOMCSRDFT) THEN
         CALL SRDFT_DENS(WRK(KDTOTAO),UDV,PVX,CMO,WRK,KFREE,LFREE)
C        Get DtotAO, allocation for DXtotAO in SIRFCK.SRDFT can here be
C        used temporarily for DVAO if NASHT.gt.0
C
         IF (IPRRSP .GT. 50) THEN
            WRITE(LUPRI,*) 'SRDFT: DTCAO matrix'
            CALL OUTPUT(WRK(KDTCAO),1,NBAST,1,NBAST,
     &         NBAST,NBAST,-1,LUPRI)
            WRITE(LUPRI,*) 'SRDFT: DTVAO matrix'
            CALL OUTPUT(WRK(KDTVAO),1,NBAST,1,NBAST,
     &         NBAST,NBAST,-1,LUPRI)
            WRITE(LUPRI,*) 'SRDFT: total density matrix AO basis DTOTAO'
            CALL OUTPUT(WRK(KDTOTAO),1,NBAST,1,NBAST,
     &         NBAST,NBAST,-1,LUPRI)
         END IF
      END IF
#endif
C
C
C Manu about MC-srDFT: Note that for ncsim >0 and doftc 
C       the core transition DM (DTCAO) is set 
C       to zero --> the corresponding IFCTYP will therefore be zero.
C       When calculating the srH contribution to the Hessian in sirfck
C       DTVAO is added to DTCAO and the srH term calculated with IFCTYP=0, in
C       other words nothing is calculated ... IFCTYP must be modified before
C       calling sirfck. see the next "Manu about srdft" ...   
C
C
C     Transform any transition DMs to the AO basis 
      IF (NCSIM.GT.0) THEN
         IF (DOFTC) CALL DZERO(WRK(KDTCAO),NCSIM*N2BASX)
         JAO = KDTVAO
         JMO = KDTV
         DO ICSIM = 1, NCSIM
            CALL FCKDEN2(.FALSE.,.TRUE.,DUMMY,WRK(JAO),CMO,WRK(JMO),
     &          WRK(KFREE),LFREE)
            JAO = JAO + N2BASX
            JMO = JMO + N2ASHX
         END DO
      END IF

#ifdef DEBUG_RSPOLI
      IF (IPRRSP .GT. 50) THEN
            write(lupri,*) 'after calling fckden2'
            WRITE(LUPRI,*) 'SRDFT: DTCAO matrix'
            CALL OUTPUT(WRK(KDTCAO),1,NBAST,1,NBAST,
     &         NBAST,NBAST,-1,LUPRI)
            WRITE(LUPRI,*) 'SRDFT: DTVAO matrix'
            CALL OUTPUT(WRK(KDTVAO),1,NBAST,1,NBAST,
     &         NBAST,NBAST,-1,LUPRI)
            WRITE(LUPRI,*) 'SRDFT: total density matrix AO basis DTOTAO'
            CALL OUTPUT(WRK(KDTOTAO),1,NBAST,1,NBAST,
     &         NBAST,NBAST,-1,LUPRI)
      END IF
#endif
C
C Construct NOSIM Inactive and Active density matrices according to eq.27
C in Chem Phys 119, 227 (1988).
C
C They should be consecutive in core. We transpose the inactive
C and active D-matrices and multiply the inactive D-matrix by 2.0
C to get things to work for RPA !!
C
      JOFFAO = 0
      DO 1400 IOSIM = 1,NOSIM
         CALL DEQ27(CMO,ZYMAT(1,IOSIM),UDV,WRK(KDXCAO+JOFFAO),
     *              WRK(KDXVAO+JOFFAO),WRK(KFREE),LFREE)
        JOFFAO = JOFFAO + N2BASX     
 1400 CONTINUE
      IF (NISHT .GT. 0) THEN
         CALL DSCAL(NOSIM*N2BASX,D2,WRK(KDXCAO),1)
      ELSE IF (DOFXC) THEN
         CALL DZERO(WRK(KDXCAO),NOSIM*N2BASX)
      END IF
C
C     Set ISYMDM and IFCTYP information
C
      CALL DZERO(WRK(KFXCAO),NFMAT*N2BASX)
      IF (NFMAT .GT. MXFMAT) THEN
         CALL QUIT( 'RSPOLI ERROR: NFMAT .gt. MXFMAT')
      END IF
      JDXCAO = KDXCAO
      IF (ROKS .OR. HSROHF) THEN
         CALL DAXPY(NOSIM*N2BASX,D1,WRK(KDXVAO),1,WRK(KDXCAO),1)
         CALL DSCAL(NOSIM*N2BASX,-D1,WRK(KDXVAO),1)
      END IF
      DO I = 1,NFMAT
C
C     Changed ISYMDM() from 1 to JWOPSY by hint from hjj, kr-may-95
C
         ISYMDM(I) = JWOPSY
C
C     IFCTYP = XY
C       X indicates symmetry about diagonal
C         X = 0 No symmetry
C         X = 1 Symmetric
C         X = 2 Anti-symmetric
C       Y indicates contributions
C         Y = 0 no contribution !
C         Y = 1 Coulomb
C         Y = 2 Exchange
C         Y = 3 Coulomb + Exchange
C
C     Check if density matrix is unsymmetric (IX=0),
C     symmetric (IX=10), antisymmetric (IX=20), or zero matrix (IX=30)
C     to threshold THRZER
C
         IX = 10 * MATSYM(NBAST,NBAST,WRK(JDXCAO),THRZER)
C        INTEGER FUNCTION MATSYM(N,NDIM,AMAT,THRZER)
C
         IF (IPRRSP .GT. 60) THEN
         write(lupri,*) 'I, IX = ', I, IX 
         END IF 
C
         IF (IX .EQ. 30) THEN
C           zero density matrix, do nothing !
            IFCTYP(I) = 0
         ELSE IF (IX .EQ. 20) THEN
C           Only exchange if antisymmetric density matrix
            IFCTYP(I) = IX + 2
         ELSE IF (TRPLET) THEN
C           only exchange
            IFCTYP(I) = IX + 2
         ELSE
C           Coulomb+exchange
            IFCTYP(I) = IX + 3
         END IF
         JDXCAO = JDXCAO + N2BASX
      END DO ! I = 1,NFMAT

#ifdef MOD_SRDFT
C
C Manu about MC-srDFT: in order to calculate the srH Kernel we must reset all values of
C      IFCTYP for FTCAO to the one obtained for FTVAO (because when calculating
C      IFCTYP, DTVAO should be considered, and not DTCAO (which is set to zero !) 
C
      IF ((NCSIM.GT.0).AND.DOFTC.AND.DOMCSRDFT) THEN     
#ifdef DEBUG_RSPOLI
         write(lupri,*) '****** IFCTYP modified for srH Kernel !!!!' 
         write(lupri,*) '****** with DOFTC true                !!!!' 
#endif
         DO I=1,NCSIM
            IFCTYP(NFMAT-2*NCSIM+I)=IFCTYP(NFMAT-NCSIM+I) 
         END DO
      END IF
C
C
C Manu about MC-srDFT: used same trick in order to calculate the srH contribution to
C                      the sigma vector when there is an orbital trial vector
C
      IF ((NOSIM.GT.0).AND.DOFXC.AND.DOMCSRDFT) THEN     
#ifdef DEBUG_RSPOLI
         write(lupri,*) '****** IFCTYP modified for srH Kernel !!!!' 
         write(lupri,*) '****** with DOFXC true                !!!!' 
#endif
         DO I=1,NOSIM
            IFCTYP(I)=IFCTYP(NOSIM+I) 
         END DO
      END IF
#endif
C
C
      IF (ROKS.OR.HSROHF) THEN
C     Change IFCTYP for all Fa matrices 
        DO I=1,NOSIM  
          IF (TRPLET) THEN
             IFCTYP(I)=10*(IFCTYP(I)/10) + 2
             IFCTYP(NOSIM+I)=10*(IFCTYP(NOSIM+I)/10) + 3
          ELSE
             IFCTYP(I)=10*(IFCTYP(I)/10) + 3
             IFCTYP(NOSIM+I)=10*(IFCTYP(NOSIM+I)/10) + 2
          END IF
        END DO
      END IF  

      IF (IPRRSP .GT. 60) THEN
         WRITE(LUPRI,*) 'RSPOLI:SIRFCK NFMAT DIRFCK',NFMAT,DIRFCK
         WRITE(LUPRI,*) 'ISYMDM ',(ISYMDM(I),I=1,NFMAT)
         WRITE(LUPRI,*) 'IFCTYP ',(IFCTYP(I),I=1,NFMAT)
         DO I=1,NOSIM
            JDXCAO = KDXCAO + (I-1)*N2BASX
            WRITE(LUPRI,*) 'RSPOLI:SIRFCK DXCAO no.',I
            CALL OUTPUT(WRK(JDXCAO),1,NBAST,1,NBAST,
     &                  NBAST,NBAST,-1,LUPRI)
         END DO
         write(lupri,*) 'calling sirfck in rspoli'
      END IF
      DFTADX = DFTADD
      DFTADD = .FALSE.
#ifdef DEBUG_RSPOLI
      write(lupri,*) 'RSPOLI: before calling sirfck'
      write(lupri,*) 'RSPOLI: '
#endif
! Manu test !!!!!!!!!!!
!     IFCTYP(1)=1  
!     IFCTYP(3)=1  
!     ISYMDM(1)=5 
#ifdef DEBUG_RSPOLI
      print *, 'just called SIRFCK'
      write(lupri,*) 'before call sirfck in rspoli'
      IF (DOFXC.AND.NOSIM.gt.0) THEN
       write(lupri,*) 'WRK(KFXCAO)='
       CALL OUTPUT(WRK(KFXCAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      IF (DOFXV.and.NOSIM.gt.0) THEN
       write(lupri,*) 'WRK(KFXVAO)='
       CALL OUTPUT(WRK(KFXVAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      IF (DOFTC.and.NCSIM.gt.0) THEN
       write(lupri,*) 'WRK(KFTCAO)='
       CALL OUTPUT(WRK(KFTCAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      IF (NCSIM.gt.0) THEN
       write(lupri,*) 'WRK(KFTVAO)='
       CALL OUTPUT(WRK(KFTVAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
!!!   
!!!   print also density matrices
!!!
      IF (DOFXC.AND.NOSIM.gt.0) THEN
       write(lupri,*) 'WRK(KDXCAO)='
       CALL OUTPUT(WRK(KDXCAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      IF (DOFXV.and.NOSIM.gt.0) THEN
       write(lupri,*) 'WRK(KDXVAO)='
       CALL OUTPUT(WRK(KDXVAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      IF (DOFTC.and.NCSIM.gt.0) THEN
       write(lupri,*) 'WRK(KDTCAO)='
       CALL OUTPUT(WRK(KDTCAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      IF (NCSIM.gt.0) THEN
       write(lupri,*) 'WRK(KDTVAO)='
       CALL OUTPUT(WRK(KDTVAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
#endif
! Manu test !!!!!!!!!!!
      CALL SIRFCK(WRK(KFXCAO),WRK(KDXCAO),NFMAT,ISYMDM,IFCTYP,
     &            DIRFCK,WRK(KFREE),LFREE)
#ifdef DEBUG_RSPOLI
      print *, 'just called SIRFCK'
      write(lupri,*) 'after call sirfck in rspoli'
      IF (DOFXC.AND.NOSIM.gt.0) THEN
       write(lupri,*) 'WRK(KFXCAO)='
       CALL OUTPUT(WRK(KFXCAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      IF (DOFXV.and.NOSIM.gt.0) THEN
       write(lupri,*) 'WRK(KFXVAO)='
       CALL OUTPUT(WRK(KFXVAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      IF (DOFTC.and.NCSIM.gt.0) THEN
       write(lupri,*) 'WRK(KFTCAO)='
       CALL OUTPUT(WRK(KFTCAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      IF (NCSIM.gt.0) THEN
       write(lupri,*) 'WRK(KFTVAO)='
       CALL OUTPUT(WRK(KFTVAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
#endif
      DFTADD = DFTADX
      IF (ROKS.OR.HSROHF) THEN
C
C Calculated matrices are (FC+FV,-FV(exch)=FV-Q)
C Input fock matrices (form SIRIFC) are of the form
C (FC+Q,FV-Q), where Q=FV(coul) + 2*FV(exch)
C adapt to this form (which works in RSPORB)
C 
C FC+Q 
         CALL DAXPY(NOSIM*N2BASX,-D1,WRK(KFXVAO),1,WRK(KFXCAO),1)
      END IF
C ***** Transform FXCAO and FXVAO to MO basis and add to FCX
C ***** and FVX, resp. (the one-index transformations of FC
C ***** and FV are added after this routine).
C           CALL AUTPV(ISYM,JSYM,U,V,
C    &                 PRPAO,NBAS,NBAST,PRPMO,NORB,NORBT,WRK,LWRK)
C
C
      IF (NOSIM .GT. 0) THEN
      JOFFAO = 0
      DO 5900 IOSIM = 1,NOSIM
         DO 300 ISYM=1,NSYM
            JSYM   = MULD2H(ISYM,KSYMOP)
            NORBI  = NORB(ISYM)
            NORBJ  = NORB(JSYM)
         IF (NORBI.EQ.0 .OR. NORBJ.EQ.0) GO TO 300
            IF (DOFXC) THEN
               CALL AUTPV(ISYM,JSYM,CMO(ICMO(ISYM)+1),CMO(ICMO(JSYM)+1),
     &              WRK(KFXCAO+JOFFAO),NBAS,NBAST,FCX(1,1,IOSIM),NORB,
     &                  NORBT,WRK(KFREE),LFREE)
            END IF
            IF (DOFXV) THEN
               CALL AUTPV(ISYM,JSYM,CMO(ICMO(ISYM)+1),CMO(ICMO(JSYM)+1),
     &              WRK(KFXVAO+JOFFAO),NBAS,NBAST,FVX(1,1,IOSIM),NORB,
     &                  NORBT,WRK(KFREE),LFREE)
            END IF
  300    CONTINUE
C
         JOFFAO = JOFFAO + N2BASX
 5900 CONTINUE
      END IF
C
C     Manu+hjaaj Aug 09:
C     transform FTCAO and FTVAO to FTC(in FTV) and FTV
C     (for DOMCSRDFT and maybe soon also DOMC (standard MCSCF))
C
      IF (NCSIM .GT. 0) THEN
      JOFFAO = 0
      DO ICSIM = 1,NCSIM
         DO 320 ISYM=1,NSYM
            JSYM   = MULD2H(ISYM,KSYMOP)
            NORBI  = NORB(ISYM)
            NORBJ  = NORB(JSYM)
         IF (NORBI.EQ.0 .OR. NORBJ.EQ.0) GO TO 320
            IF (DOFTC) THEN
C              FCT exists with srDFT contribution for srDFT!!!
C              NB! saved in FVTD(*,*,NCSIM+ICSIM) after FVT
               CALL AUTPV(ISYM,JSYM,CMO(ICMO(ISYM)+1),CMO(ICMO(JSYM)+1),
     &              WRK(KFTCAO+JOFFAO),NBAS,NBAST,FVTD(1,1,NCSIM+ICSIM),
     &              NORB,NORBT,WRK(KFREE),LFREE)
            END IF
            IF (DOFTV) THEN
               CALL AUTPV(ISYM,JSYM,CMO(ICMO(ISYM)+1),CMO(ICMO(JSYM)+1),
     &              WRK(KFTVAO+JOFFAO),NBAS,NBAST,FVTD(1,1,ICSIM),
     &              NORB,NORBT,WRK(KFREE),LFREE)
            END IF
  320    CONTINUE
C
         JOFFAO = JOFFAO + N2BASX
      END DO ! ICSIM = 1,NCSIM
      END IF ! NCSIM .gt. 0
C
      CALL MEMREL('RSPOLI.sirfck',WRK,KFRSAV,KWRK0,KFREE,LFREE)
C
C Manu debug1
         IF (IPRRSP.GT.30) THEN
            DO 1015 ICSIM = 1,NCSIM
            WRITE(LUPRI,'(/2A,I5)')' FVTD before aotomo'
     *         ,' FOR VECTOR ',ICSIM
            CALL OUTPUT(FVTD(1,1,ICSIM),
     *                     1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            WRITE(LUPRI,'(/2A,I5)')'V^[2c]sr_Hxc before aotomo'
     *         ,' FOR VECTOR ',ICSIM
            CALL OUTPUT(FVTD(1,1,NCSIM+ICSIM),
     *                     1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
 1015       CONTINUE
         END IF
C Manu debug2
C
C The following would comply for a call from LAGRAN in HF
C
      IF (NCSIM.EQ.1 .AND. TDHF  .AND. TRPLET) THEN
C
C For high spin - dft we just unpack FC and FV
C
         IF (HSROHF.OR.DFTADD) THEN
            CALL MEMGET('REAL',KFC,N2ORBX,WRK,KFREE,LFREE)
            CALL MEMGET('REAL',KTMP,NNORBX,WRK,KFREE,LFREE)
            CALL PKSYM1(WRK(KTMP),FV,NORB,NSYM,-1)
            CALL DSPTSI(NORBT,WRK(KTMP),FVTD)
            CALL PKSYM1(WRK(KTMP),FC,NORB,NSYM,-1)
            CALL DSPTSI(NORBT,WRK(KTMP),WRK(KFC))
         ELSE
            CALL MEMGET('REAL',KDA,N2BASX,WRK,KFREE,LFREE)
            CALL MEMGET('REAL',KFA,N2BASX,WRK,KFREE,LFREE)
            CALL FCKDEN(
     &         .FALSE.,.TRUE.,DUMMY,WRK(KDA),CMO,WRK(KDTV),
     &          WRK(KFREE),LFREE
     &         )
            CALL DZERO(WRK(KFA),N2BASX)
            ISYMDM(1)=1
            IFCTYP(1)=12
            CALL SIRFCK(
     &        WRK(KFA),WRK(KDA),1,ISYMDM,IFCTYP,DIRFCK,WRK(KFREE),LFREE)
            CALL AOTOMO(WRK(KFA),FVTD,CMO,1,WRK(KFREE),LFREE)
C           CALL AOTOMO(XAO,XMO,CMO,XSYM,WRK,LWRK)
         END IF
      END IF
 2000 CONTINUE
C
C From here the code is as before DIRFCK
C
      IF (NOSIM .GT. 0) THEN
C
C     1) add DFT contributions
C
         IF (DODFT) THEN
#ifdef DEBUG_RSPOLI
           write(lupri,*) 'DODFT: before call dft_lin_resp*'
           IF (DOFXC) THEN
             write(lupri,*) 'FCX before dft_lin_resp*='
             CALL OUTPUT(FCX,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
           END IF
           IF (DOFXV) THEN
             write(lupri,*) 'FCV before dft_lin_resp*='
             CALL OUTPUT(FCV,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
           END IF
#endif
C          Special case of alpha density only (old code)  
           IF((NISHT.EQ.0) .AND. (NASHT.GT.0)) THEN
             DO IOSIM = 1, NOSIM
               call dft_lin_respab(fcx(1,1,iosim),fvx(1,1,iosim),
     &              cmo,zymat(1,iosim),trplet,ksymop,wrk(kfree),
     &              lfree)
             END DO
           ENDIF 
C          Proceed normal way otherwise 
           IF ((NASHT.GT.0) .AND. (NISHT.GT.0)) THEN 
             call dft_lin_respab_b(NOSIM,fcx,fvx,cmo,zymat,
     &                           trplet, ksymop,wrk(kfree),lfree) 
           ELSE
             call dft_lin_respf(NOSIM,fcx,cmo,zymat,
     &                         trplet, ksymop,wrk(kfree),lfree)
           ENDIF
#ifdef DEBUG_RSPOLI
           write(lupri,*) 'DODFT: after call dft_lin_resp*'
           IF (DOFXC.AND.NOSIM.gt.0) THEN
             write(lupri,*) 'FCX after dft_lin_resp*='
             CALL OUTPUT(FCX,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
           END IF
           IF (DOFXV.and.NOSIM.gt.0) THEN
             write(lupri,*) 'FCV after dft_lin_resp*='
             CALL OUTPUT(FCV,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
           END IF
#endif
         END IF  ! IF (DODFT) THEN
C
C        ADD CONTRIBUTION FROM THE A(0)S(2) TERM IF HRPA
C
         IF (SOPPA) THEN
           CALL A0S2(FCX,FC,UDV,ZYMAT,WRK(KTZYMT),WRK(KH2X),WRK(KH2XP),
     *             WRK(KCOEUN),NOSIM, WRK,LWRK)
C
C  Cancel out the nonsymmetric parts of A(2) matrix. Use the explicitly
C  constructed A(2) matrix in XINDX.
C
           IF (.NOT. A2EXIST) THEN
             CALL DSCAL(N2ORBX,DP5,XINDX(KAB2),1)
             IF (IPRRSP .GT. 40) THEN
                  WRITE(LUPRI,'(/A)')
     &      ' SOPPA A(2)(ai,bj) packed as .5*A(2)(i,j) and .5*A(2)(a,b)'
                  CALL OUTPUT(XINDX,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
             END IF
             A2EXIST = .TRUE.
           END IF
           CALL HRPAHM(FCX,XINDX(KAB2),ZYMAT,NOSIM)
         END IF
C
C ADD CONTRIBUTION TO FCX AND FVX FROM ONE-INDEX TRANSFORMED
C TOTAL SYMMETRIC FOCK MATRICES
C
         CALL FCKOIN(NOSIM,FC,FV,ZYMAT,FCX,FVX)
C
C
         IF ( (KZCONF .GT. 0) .AND. .NOT.SOPPA) THEN
C
C HALF TRANSFORMED INTEGRALS ARE NOW IN H2XAC AND CONFIGURATION
C PART OF LINEAR TRANSFORMATION CAN BE CARRIED OUT FOR EACH ORBITAL
C TRIAL VECTOR
C
            CALL H2XSIG(NOSIM,FCX,ZYMAT,WRK(KH2XAC),
     *                  EVECS(1+NCSIM*KZYVAR),
     *                  XINDX,WRK(KFREE),LFREE)
C
         END IF
C
C
C DISTRIBUTE FOCK AND Q MATRICES IN LINEAR TRANSFORMED VECTORS
C
         IF ((IPRRSP.GT.30).AND.(NASHT.GT.1)) THEN
            DO 700 IOSIM = 1,NOSIM
               WRITE(LUPRI,'(/A,I5,A)') ' QAX FOR ',IOSIM,
     *             ' ORBITAL TRIAL VECTOR'
               CALL OUTPUT(QAX(1,1,IOSIM),1,NORBT,1,NASHT,
     *                     NORBT,NASHT,1,LUPRI)
               WRITE(LUPRI,'(/A,I5,A)') ' QBX FOR',IOSIM,
     *             ' ORBITAL TRIAL VECTOR'
               CALL OUTPUT(QBX(1,1,IOSIM),1,NORBT,1,NASHT,
     *                     NORBT,NASHT,1,LUPRI)
 700        CONTINUE
         END IF
         IF (IPRRSP.GT.50) THEN
            DO 1012 IOSIM = 1,NOSIM
               WRITE(LUPRI,'(/2A,I5)')' FCX total contribution'
     *                              ,' FOR VECTOR ',IOSIM
               CALL OUTPUT(FCX(1,1,IOSIM),
     *                     1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
               IF (DOFXV) THEN
                  WRITE(LUPRI,'(/2A,I5)')' FVX total contribution'
     *                                 ,' FOR VECTOR ',IOSIM
                  CALL OUTPUT(FVX(1,1,IOSIM),
     *                           1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
               END IF
 1012       CONTINUE
         END IF
#ifdef DEBUG_RSPOLI
         write(lupri,*) 'I call RSPORB with ONEIND true'  
#endif
         CALL RSPORB(.TRUE.,NOSIM,FC,FCX,FVX,QAX,QBX,UDV,
     *               EVECS(NCSIM*KZYVAR+1),ZYMAT)
C        CALL RSPORB(ONEIND,NSIM,FC,FCX,FVX,QAX,QBX,UDV,EVECS,ZYMAT)

      END IF   ! IF (NOSIM.GT.0) THEN
C ------------------------------------------------------------------
      IF (NCSIM.GT.0) THEN
         IF (DOMCSRDFT) THEN
C           NOTE, the TDM passed to SIRFCK is really
C               minus the TDM, thus we need to multiply
C               the FTC matrices in FVTD(1,1,NCSIM+ICSIM) by minus 1
C               to get the right contribution from the DFT
C               module in srDFT.
C
             CALL DSCAL(NCSIM*N2ORBX,DM1,FVTD(1,1,NCSIM+1),1)
         END IF
C
         IF (IPRRSP.GT.30 .AND. .NOT.SOPPA) THEN
            DO 800 ICSIM = 1,NCSIM
               WRITE(LUPRI,'(/A,I5,A)') ' QATD FOR ',ICSIM,
     *             ' CONFIGURATION TRIAL VECTOR'
               CALL OUTPUT(QATD(1,1,ICSIM),1,NORBT,1,NASHT,
     *                     NORBT,NASHT,1,LUPRI)
               WRITE(LUPRI,'(/A,I5,A)') ' QBTD FOR',ICSIM,
     *             ' CONFIGURATION TRIAL VECTOR'
               CALL OUTPUT(QBTD(1,1,ICSIM),1,NORBT,1,NASHT,
     *                     NORBT,NASHT,1,LUPRI)
 800        CONTINUE
         END IF
         IF (IPRRSP.GT.30) THEN
            DO 1014 ICSIM = 1,NCSIM
            WRITE(LUPRI,'(/2A,I5)')' FVTD total contribution'
     *         ,' FOR VECTOR ',ICSIM
            CALL OUTPUT(FVTD(1,1,ICSIM),
     *                     1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            WRITE(LUPRI,'(/2A,I5)')'V^[2c]sr_Hxc total contribution'
     *         ,' FOR VECTOR ',ICSIM
            CALL OUTPUT(FVTD(1,1,NCSIM+ICSIM),
     *                     1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
 1014       CONTINUE
         END IF
         IF (HSROHF.OR.DODFT) THEN
            CALL RSPORB(.TRUE.,NCSIM,FC,WRK(KFC),FVTD,QATD,QBTD,UDV,
     *               EVECS,DUMMY)
C
C Note: This is to cancel the change of sign in RSPORB
C       In ordinary ROHF/MCSCF, RSPORB is called with negative
C       matrices (due to negative active densities produced by
C       the call to RSPTDM above when the CI vector is the reference state)
C      
C
            CALL DSCAL(KZYVAR,DM1,EVECS,1)
         ELSE
C           (FVTD(1,1,NCSIM+1) contains the V[2c]xcsr for srDFT,
C            it is NOT used when normal MCSCF)
            CALL RSPORB(.FALSE.,NCSIM,FC,FVTD(1,1,NCSIM+1),FVTD,
     *                  QATD,QBTD,WRK(KDTV),EVECS,DUMMY)
         END IF
#ifdef MOD_SRDFT
         IF (DOMCSRDFT) THEN
C           ... special MCSRDFT contribution to csf sigma-vectors,
C               stored in FTV(1,NCSIM+ICSIM) in SIRTR1 after the 
C               standard FTV matrices.
C               In paper: V^([2c]xc-SR) matrix.
C
C           Aug. 09: structure inspired from rspsol:SLVSC routine /hjaaj
C
            CALL MEMGET('REAL',KSRCVEC ,KZCONF,WRK,KFREE,LFREE)
            CALL MEMGET('REAL',KSRW    ,N2ASHX,WRK,KFREE,LFREE)
            JEVECS = 1
            DO ICSIM = 1,NCSIM
               CALL GETAC1(FVTD(1,1,NCSIM+ICSIM),WRK(KSRW))
               IF (IPRRSP .GT. 40) THEN
                  WRITE(LUPRI,'(/A,I3)')
     &            'MCSRDFT "FTC"="V^([2c],xc-SR)" no.',ICSIM
                  CALL OUTPUT(FVTD(1,1,NCSIM+ICSIM),
     &               1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                  WRITE(LUPRI,'(/A,I3)') 'MCSRDFT SRLTRAC no.',ICSIM
                  CALL OUTPUT(WRK(KSRW),
     &               1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
               END IF
               JSPIN1 = 0
               JSPIN2 = 0
C
#ifdef DEBUG_RSPOLI
               IF (IPRRSP.GT.101) THEN
                  WRITE(LUPRI,*)
     *            ' Z CONF  PART OF LINEAR TRANSFORMED',
     *            ' before Vsr[2c] contr.',
     *            'CONF TRIAL VECTOR'
                  CALL OUTPUT(EVECS(JEVECS),1,KZVAR,1,2,
     *                                        KZVAR,2,1,LUPRI)
                  write(lupri,*) 'kzvar = ', kzvar  
                  write(lupri,*) 'kzconf = ', kzconf  
                END IF
#endif
C
C     Z part of linear transformation (operator: V[2c])
C     (The Y part configuration contribution is zero) <<<---- Manu: no, it is
C                                                                   NOT !
C
               CALL CISIGD(IREFSY,KSYMST,NCREF,KZCONF,WRK(KCREF),
     &                  EVECS(JEVECS),WRK(KSRW),DUMMY,
     &                  .TRUE.,.FALSE.,XINDX,JSPIN1,JSPIN2,
     &                  WRK(KFREE),LFREE)
#ifdef DEBUG_RSPOLI
               write(lupri,*) 'Z part ... V[2c]'
               write(lupri,*) 'IREFSY = ', IREFSY
               write(lupri,*) 'KSYMST = ', KSYMST
#endif

               IF (IREFSY .EQ. KSYMST) THEN
                  FAC = DDOT(KZCONF,EVECS(JEVECS),1,WRK(KCREF),1)
                  CALL DAXPY(KZCONF,(-FAC),WRK(KCREF),1,EVECS(JEVECS),1)
               END IF
C
               IF (IPRRSP.GT.101) THEN
                    WRITE(LUPRI,*)
     *              ' srDFT Z CONF  PART OF LINEAR TRANSFORMED ',
     *              'CONF TRIAL VECTOR'
                    WRITE(LUPRI,*)' REFERENCE COMPONENT PROJECTED OUT',
     &              ', factor = ',-FAC
                    CALL OUTPUT(EVECS(JEVECS),1,KZVAR,1,2,
     *                                        KZVAR,2,1,LUPRI)
               END IF
C
! Manu: build the Y part of the sigma vector (due to the Vsr[2c] operator)
C
C
               JYEVECS=JEVECS+KZVAR  
!                                         
C
C        
               CALL CISIGD(IREFSY,KSYMST,NCREF,KZCONF,WRK(KCREF),
     &                  EVECS(JYEVECS),WRK(KSRW),DUMMY,
     &                  .TRUE.,.FALSE.,XINDX,JSPIN1,JSPIN2,
     &                  WRK(KFREE),LFREE)
               IF (IREFSY .EQ. KSYMST) THEN
                  FAC = DDOT(KZCONF,EVECS(JYEVECS),1,WRK(KCREF),1)
                  CALL DAXPY(KZCONF,(-FAC),WRK(KCREF),1,EVECS(JYEVECS),
     &                       1)
!
!

               END IF
!
! Manu:        For the Y conf. part, the transition density matrix changes sign.
!              Therefore we have to multiply the Vsr[2c] contribution to 
!              the Y part of the sigma vector by -1 
!
               CALL DSCAL(KZCONF,DM1,EVECS(JYEVECS),1)  
!
               IF (IPRRSP.GT.101) THEN
                    WRITE(LUPRI,*)
     *              ' srDFT Z and Y CONF  PART OF LINEAR TRANSFORMED ',
     *              'CONF TRIAL VECTOR'
                    WRITE(LUPRI,*)' REFERENCE COMPONENT PROJECTED OUT',
     &              ', factor = ',-FAC
                    CALL OUTPUT(EVECS(JEVECS),1,KZVAR,1,2,
     *                                        KZVAR,2,1,LUPRI)
               END IF
C
C
C     Special MCSRDFT correction, corresponding to the effective operator
C        V^([2c],xc-SR)
C     is saved in FVTD(*,*,NCSIM+ICSIM) in RSPOLI.
C
               CALL RSP_SRDFTSO( .FALSE.,.TRUE.,UDV,
     &             FVTD(1,1,NCSIM+ICSIM),EVECS(JEVECS))
C
               JEVECS = JEVECS + KZYVAR
            END DO ! ICSIM = 1,NCSIM
         END IF ! DOMCSRDFT
#endif
      END IF ! NCSIM .gt. 0
      CALL MEMREL('RSPOLI.bottom',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
C
C  *** Perform supersymmetry averaging
C
      IF(RSPSUP .AND. (KZWOPT.GT. 0) ) THEN
         CALL RFANTI(NOSIM,EVECS(NCSIM*KZYVAR+1),ZYMAT,WRK,LWRK)
      END IF
      IF ( RSPSUP .AND. (KSYMOP. EQ. 1) .AND. (KZWOPT .GT. 0)) THEN
         CALL RSPAVE(EVECS(KZCONF+1),KZVAR,2*(NCSIM+NOSIM))
      END IF
C
C END OF RSPOLI
C
      CALL QEXIT('RSPOLI')
      RETURN
      END
