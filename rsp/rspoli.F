C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
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
RSPOLI: added call to HRPA which adds A(2) and B(2) matrices to FCONE
920721-hinne hettema
RSPOLI: implemented RSPSUP (super symmetry averaging)
===========================================================================
#endif
C  /* Deck rspoli */
      SUBROUTINE RSPOLI(NCSIM,NOSIM,UDV,ZYCVEC,FC,FV,PVX,ZYMAT,FCONE,
     *                  FVONE,QAONE,QBONE,FVTD,QATD,QBTD,EVECS,
     *                  H2,H2X,XINDX,CMO,WRK,LWRK)
C
C Copyright 13 FEB 1986 Poul Joergensen and Hans Joergen Aa. Jensen
C
C 21-1-1992: Changed for Direct RPA (KEYWORD: DIRFCK), HA.
C 21-7-1992: Supersymmetry averaging added. Hinne Hettema
C 28-2-1995: Changed for Direct RPA with AO integrals read
C            from disk, NOITRA=.TRUE. P&D /is not included?/hjaaj
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
C       TO QAONE AND QBONE MATRICES
C    3) DISTRIBUTE FOCK CORE , FOCK VALENCE , QA AND QB MATRICES
C       INTO ORBITAL PART OF LINEAR TRANSFORMED VECTORS
C
#include <implicit.h>
C
      DIMENSION UDV(NASHDI,*),ZYCVEC(*),FC(*),FV(*),PVX(*)
      DIMENSION ZYMAT(N2ORBX,*),FCONE(NORBT,NORBT,*)
      DIMENSION FVONE(NORBT,NORBT,*),QAONE(NORBT,NASHDI,*)
      DIMENSION QBONE(NORBT,NASHDI,*),FVTD(NORBT,NORBT,*)
      DIMENSION QATD(NORBT,NASHDI,*),QBTD(NORBT,NASHDI,*)
      DIMENSION EVECS(*),H2(NORBT,*),H2X(NORBT,*)
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
      LOGICAL DFTADX, ROKS
      DIMENSION ISYMDM(MXFMAT), IFCTYP(MXFMAT)
      DIMENSION NEEDMU(6),NEEDDI(6)
      DATA NEEDMU /1,1,1,0,0,0/
      DATA NEEDDI /1,1,1,0,0,0/
      ROKS = NASHT.GT.0 .AND. DFTADD
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
      IF(ROKS.AND..NOT.TDHF) THEN
         CALL QUIT('RSPOLI: ROKS not implemented for '
     &        //   'MO-integral method')
      END IF

C
C ALLOCATE WORK SPACE FOR DVT and PVTD and SOPPA
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWRK
      CALL MEMGET('REAL',KDVT ,NCSIM*N2ASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KPVTD,NCSIM*N2ASHX*N2ASHX,WRK,KFREE,LFREE)
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
         CALL GETREF(WRK(KCREF),NCREF)
         IF (TRPLET) THEN
            ISPIN1 = 1
            ISPIN2 = 0
         ELSE
            ISPIN1 = 0
            ISPIN2 = 0
         END IF
         CALL RSPTDM(NCSIM,IREFSY,KSYMST,NCREF,KZCONF,WRK(KCREF),
     *                 ZYCVEC,WRK(KDVT),WRK(KPVTD),
     *                 ISPIN1,ISPIN2,.TRUE.,.FALSE.,
     *                 XINDX,WRK,KFREE,LFREE)
C        CALL RSPTDM(NCSIM,ILRESY,IRSYM,NCLREF,NCRDIM,CLREF,
C    *                 CR, RHO1,RHO2, ISPIN1,ISPIN2,TDM,NORHO2,
C    *                 XNDXCI,WORK,KFREE,LFREE)
         CALL MEMREL('RSPOLI.tdm',WRK,KFRSAV,KWRK0,KFREE,LFREE)
      END IF

C
C INITIALIZE ONE INDEX TRANSFORMED AND TRANSITION DENSITY Q,FC
C AND FV MATRICES
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
         CALL DZERO(FCONE,N2ORBX*NOSIM)
         IF (NASHT.GT.0) THEN
            CALL DZERO(FVONE,N2ORBX*NOSIM)
            NTOT = NORBT*NASHT*NOSIM
            CALL DZERO(QAONE,NTOT)
            CALL DZERO(QBONE,NTOT)
            IF (KZCONF.GT.0 ) THEN
               NTOT = NOSIM*N2ASHX*NNASHX
               CALL DZERO(WRK(KH2XAC),NTOT)
            END IF
         END IF
      END IF
      IF (NCSIM.GT.0) THEN
         NTOT = NORBT * NORBT * NCSIM
         CALL DZERO(FVTD,NTOT)
         NTOT = NORBT*NASHT*NCSIM
         IF (NTOT .GT. 0) THEN
C        ... is zero for SOPPA
            CALL DZERO(QATD,NTOT)
            CALL DZERO(QBTD,NTOT)
         END IF
      END IF
C
C IF DIRFCK (direct RPA) compute new densities
C and Fock matrices in AO basis
C
      IF (DIRFCK .OR. HSROHF .OR. DODFT) GO TO 1000

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
C not DIRFCK, i.e.
C calculate Fock matrices etc. in MO basis
C
C ALLOCATE WORK SPACE FOR MULLIKEN and DIRAC DISTRIBUTIONs
C
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
      KWRK1 = KFREE
      IF ((NASHT .GT. 0) .AND. (NOSIM .GT.0)) THEN
         CALL RSPTR1(NOSIM,UDV,ZYMAT,WRK(KDENA),WRK(KDENB))
      END IF
C
C setup for reading Mulliken MO integrals ...
C NEEDED DISTRIBUTIONS DEFINED IN NEEDMU(6)
C
      CALL MEMGET('REAL',KPVCD ,      N2ASHX,WRK,KFREE,LFREE)
      IDIST = 0
 90   CALL NXTH2M(IC,ID,H2,NEEDMU,WRK,KFREE,LFREE,IDIST)
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
C CONTRIBUTIONS TO FCONE AND FVONE
C
         IF (ITYPCD .LE. 3) CALL FONEMU(NOSIM,ICI,IDI,H2,
     *               FCONE,FVONE,ZYMAT,WRK(KDENA),WRK(KDENB),
     *               WRK(KWRK2),LWRK2)
C        CALL FONEMU(NSIM,ICI1,IDI1,H2,
C    *                  FCOEX,FVOEX,ZYMAT,DENA,DENB,WRK,LWRK)
C
Cmjp CONTRIBTUIONS FORM THE HRPA TERMS
         IF (SOPPA) THEN
            CALL HRPA(FCONE,UDV,PVX,H2,ZYMAT,WRK(KTZYMT),H2X,WRK(KH2XP),
     *                WRK(KCOEUN),ICI,IDI,ICDSYM,ITYPCD,NOSIM,
     *                XINDX(KIADR1),WRK(KWRK2),LWRK2)
            IF (ITYPCD .EQ. 4) THEN
               IF (.NOT.HIRPA .AND. KZCONF .GT. 0) THEN
C                 ... skip SOPH2X for .OPTORB iterations
                  CALL SOPH2X(EVECS(1+NCSIM*KZYVAR),ZYMAT,WRK(KTZYMT),
     *                        H2,ICI,IDI,ICDSYM,NOSIM,XINDX(KABSAD),
     &                        XINDX(KABTAD),XINDX(KIJSAD),XINDX(KIJTAD),
     &                        XINDX(KIJ1AD),XINDX(KIJ2AD),XINDX(KIJ3AD),
     *                        WRK(KWRK2),LWRK2)
               END IF
C  Construct the A(2) matrix explicitly (but only once)
C  for the symmetrization of A(2)b
               IF (.NOT.A2EXIST) THEN
                  CALL HRPAA2(PVX,H2,XINDX(KAB2),WRK(KCOEUN),
     &                        XINDX(KIADR1),ICI,IDI,ICDSYM)
               END IF
            END IF
         END IF
         IF ((NASHT.GT.1).AND.(ITYPCD.EQ.3)) THEN
C
C
C CONTRIBUTIONS TO QAONE AND QBONE
C
            CALL QONEMU(NOSIM,NCIW,NDIW,ICDSYM,
     *                  QAONE,QBONE,ZYMAT,H2,H2X,WRK(KH2XAC),
     *                  PVX,WRK(KPVCD),WRK(KWRK2),LWRK2)
C           CALL QONEMU(NOSIM,NCIW,NDIW,ICDSYM,QAONE,QBONE,
C    *                  ZYMAT,H2,H2X,H2XAC,PVX,PVCD,WRK,LWRK)
         END IF
      END IF
      IF (NCSIM.GT.0) THEN
         IF (SOPPA) THEN
            IF (.NOT.HIRPA .AND. ITYPCD.EQ.4) THEN
               CALL SOPPAF(FVTD,H2,ZYCVEC,ICI,IDI,ICDSYM,NCSIM,
     *                     XINDX,WRK(KWRK2),LWRK2)
            END IF
         ELSE
            IF (ITYPCD.EQ.2 .OR. ITYPCD.EQ.3 .OR. ITYPCD.EQ.5) THEN
C
C CONTRIBUTIONS TO FVTD
C
               CALL FTDMU(NCSIM,ICI,IDI,FVTD,WRK(KDVT),H2)
C              CALL FTDMU(NSIM,ICI1,IDI1,FVTD,DVT,H2)
            END IF
            IF (ITYPCD.EQ.3) THEN
C
C CONTRIBUTIONS TO QATD AND QBTD
C
              CALL QTD(NCSIM,NCIW,NDIW,ICDSYM,QATD,QBTD,H2,WRK(KPVTD),
     *                 WRK(KPVCD),WRK(KWRK2),LWRK2)
C             CALL QTD(NCSIM,NCIW,NDIW,ICDSYM,QATD,QBTD,H2,PVTD,PVDEN,
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
     *      ' QAONE FOR',IOSIM,' ORBITAL TRIAL VECTOR'
         CALL OUTPUT(QAONE(1,1,IOSIM),
     *                     1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
         WRITE(LUPRI,'(/A,I5,A)')
     *      ' QBONE FOR',IOSIM,' ORBITAL TRIAL VECTOR'
         CALL OUTPUT(QBONE(1,1,IOSIM),
     *                     1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
 1001    CONTINUE
      END IF
      IF ((IPRRSP.GT.50).AND.(NOSIM.GT.0)) THEN
         DO 1002 IOSIM = 1,NOSIM
         WRITE(LUPRI,'(/2A,I5)')' FCOEX MULLIKEN DISTRIBUTION CONTRB'
     *      ,' FOR VECTOR ',IOSIM
         CALL OUTPUT(FCONE(1,1,IOSIM),
     *                     1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(/2A,I5)')' FVOEX MULLIKEN DISTRIBUTION CONTRB'
     *      ,' FOR VECTOR ',IOSIM
         CALL OUTPUT(FVONE(1,1,IOSIM),
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
      IF ((IPRRSP.GT.50).AND.(NCSIM.GT.0)) THEN
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
      LUINTD = 0
      IF (.NOT. NEWTRA) CALL GPOPEN(LUINTD,'DRCINT','OLD',' ',
     &                              'UNFORMATTED',IDUMMY,.FALSE.)
 1090 CALL NXTH2D(IC,ID,H2,NEEDDI,LUINTD,WRK,KFREE,LFREE,IDIST)
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
            CALL DGETRN(H2,NORBT,NORBT)
         ENDIF
         NCIW  = ICIW - NISHT
         NDIW  = IDIW - NISHT
C
         IF ( IPRRSP.GT.150 ) THEN
            WRITE(LUPRI,'(/A,2I8)')' DIRAC DISTRIBUTION : IC,ID ',IC,ID
            CALL OUTPUT(H2(1,1),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         ENDIF
         IF (NOSIM.GT.0) THEN
            CALL FONEDR(NOSIM,ICI,IDI,H2,WRK(KFCOCO),WRK(KFVOCO),
     *                  FCONE,FVONE,ZYMAT,WRK(KDENA),WRK(KDENB),
     *                  WRK(KFREE),LFREE)
C           CALL FONEDR(NSIM,ICI1,IDI1,H2D,FCOCO,FVOCO,
C    *                     FCOEX,FVOEX,ZYMAT,DENA,DENB,WRK,LWRK)
C
C CONTRIBUTIONS TO QAONE AND QBONE
C
            IF( (IOBTYP(ICI).EQ.JTACT) .AND. (IOBTYP(IDI).EQ.JTACT))
     *      CALL QONEDI(NOSIM,NCIW,NDIW,ICDSYM,QAONE,QBONE,ZYMAT,
     *            H2,H2X,PVX,WRK(KPVCD2),WRK(KPVCD3),WRK(KPVDC2),
     *            WRK(KPVDC3),WRK(KFREE),LFREE)
C           CALL QONEDI(NOSIM,NCIW,NDIW,ICDSYM,
C    *               QAONE,QBONE,ZYMAT,H2,H2X,
C    *               PVX,PVCD2,PVCD3,PVDC2,PVDC3,WRK,LWRK)
         END IF
         IF (NCSIM.GT.0 .AND. .NOT.SOPPA) THEN
            CALL FTDDR(NCSIM,ICI,IDI,FVTD,WRK(KDVT),H2)
C           CALL FTDDR(NSIM,ICI1,IDI1,FVTD,DVT,H2D)
         END IF
C
         GO TO 1090
C
C ALL DIRAC DISTRIBUTIONS HAVE BEEN READ IN
C
 1095 CONTINUE
      IF (LUINTD .GT. 0) CALL GPCLOSE(LUINTD,'KEEP')
C
C ADD CONTRIBUTION TO FCONE AND FVONE FROM ONE INDEX TRANSFORMED
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
               CALL OUTPUT(FCONE(1,1,IOSIM),
     *                     1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
               WRITE(LUPRI,'(/2A,I5)')' FCOCO CONTRB'
     *                              ,' FOR VECTOR ',IOSIM
               CALL OUTPUT(WRK(KFCOCO+ (IOSIM-1)*NORBT*NORBT),
     *                     1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
               IF (NASHT.GT.0) THEN
                  WRITE(LUPRI,'(/2A,I5)')' FVOEX  CONTRB'
     *                                 ,' FOR VECTOR ',IOSIM
                  CALL OUTPUT(FVONE(1,1,IOSIM),
     *                        1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                  WRITE(LUPRI,'(/2A,I5)')' FVOCO CONTRB'
     *                                 ,' FOR VECTOR ',IOSIM
                  CALL OUTPUT(WRK(KFVOCO+ (IOSIM-1)*NORBT*NORBT),
     *                        1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
               END IF
 1010       CONTINUE
         END IF
         CALL FONEDN(NOSIM,WRK(KFCOCO),WRK(KFVOCO),FCONE,FVONE)
      ENDIF
C
      CALL MEMREL('RSPOLI.h2m+d',WRK,KFRSAV,KWRK0,KFREE,LFREE)
      GO TO 2000
 1000 CONTINUE
C
C end up here if DIRFCK
C       
C       Allocate work memory:
C       FXCAO,FXVAO must be stored consecutively
C       DXCAO,DXVAO must be stored consecutively
C       (requirement for SIRFCK call)
C
      NFMAT  = 0
      IF (NISHT .GT. 0) NFMAT = NFMAT + NOSIM
C     ... for FXCAO
      IF (NASHT .GT. 0) NFMAT = NFMAT + NOSIM
C     ... for FXVAO
C
      CALL MEMGET('REAL',KFXCAO,NFMAT*N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDXCAO,NFMAT*N2BASX,WRK,KFREE,LFREE)
      IF (NASHT .GT. 0) THEN
         KFXVAO = KFXCAO + NOSIM*N2BASX
         KDXVAO = KDXCAO + NOSIM*N2BASX
      ELSE
         KFXVAO = KFXCAO
         KDXVAO = KDXCAO
      END IF

C
C Construct NOSIM Inactive and Active density matrices according to eq.27
C CP, 119, 227 (1988).
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
      IF (NISHT .GT. 0)
     &   CALL DSCAL(NOSIM*N2BASX,D2,WRK(KDXCAO),1)
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
      END DO
C     Change IFCTYP for all Fa matrices 
      IF (ROKS.OR.HSROHF) THEN
      DO 1443 I=1,NOSIM  
          IF (TRPLET) THEN
             IFCTYP(I)=10*(IFCTYP(I)/10) + 2
             IFCTYP(NOSIM+I)=10*(IFCTYP(NOSIM+I)/10) + 3
          ELSE
             IFCTYP(I)=10*(IFCTYP(I)/10) + 3
             IFCTYP(NOSIM+I)=10*(IFCTYP(NOSIM+I)/10) + 2
          END IF
 1443 CONTINUE   
      END IF  
      DFTADX = DFTADD
      DFTADD = .FALSE.
      CALL SIRFCK(WRK(KFXCAO),WRK(KDXCAO),NFMAT,ISYMDM,IFCTYP,
     &            DIRFCK,WRK(KFREE),LFREE)
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
C ***** Transform FXCAO and FXVAO to MO basis and add to FCONE
C ***** and FVONE, resp. (the one-index transformations of FCONE
C ***** and FVONE are added after this routine).
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
            IF (NISHT .GT. 0) THEN
               CALL AUTPV(ISYM,JSYM,CMO(ICMO(ISYM)+1),CMO(ICMO(JSYM)+1),
     &              WRK(KFXCAO+JOFFAO),NBAS,NBAST,FCONE(1,1,IOSIM),NORB,
     &                  NORBT,WRK(KFREE),LFREE)
            END IF
            IF (NASHT .GT. 0) THEN
               CALL AUTPV(ISYM,JSYM,CMO(ICMO(ISYM)+1),CMO(ICMO(JSYM)+1),
     &              WRK(KFXVAO+JOFFAO),NBAS,NBAST,FVONE(1,1,IOSIM),NORB,
     &                  NORBT,WRK(KFREE),LFREE)
            END IF
  300    CONTINUE
C
         JOFFAO = JOFFAO + N2BASX
 5900 CONTINUE
      END IF
C
      CALL MEMREL('RSPOLI.sirfck',WRK,KFRSAV,KWRK0,KFREE,LFREE)
C
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
     &         .FALSE.,.TRUE.,DUMMY,WRK(KDA),CMO,WRK(KDVT),
     &          WRK(KFREE),LFREE
     &         )
            CALL DZERO(WRK(KFA),N2BASX)
            ISYMDM(1)=1
            IFCTYP(1)=12
            CALL SIRFCK(
     &         WRK(KFA),WRK(KDA),1,ISYMDM,IFCTYP,DIRFCK,WRK(KFREE),LFREE
     &         )
            CALL AOTOMO(WRK(KFA),FVTD,CMO,1)
         END IF
      END IF
 2000 CONTINUE
C
C From here the code is as before DIRFCK
C
C add DFT contribution 
C
      IF (DODFT) THEN
         DO IOSIM = 1, NOSIM
            if(NASHT.GT.0) THEN
               call dft_lin_respab(fcone(1,1,iosim),fvone(1,1,iosim),
     &                   cmo,zymat(1,iosim),trplet,ksymop,wrk(kwrk0),
     &                   lfree)
            ELSE
               call dft_lin_respf(fcone(1,1,iosim),cmo,zymat(1,iosim),
     &                           trplet, ksymop,wrk(kwrk0),lfree)
c              CALL DFTATR(FCONE(1,1,IOSIM),CMO,ZYMAT(1,IOSIM),TRPLET,
c     &                    KSYMOP,WRK(KWRK0),LFREE)
            ENDIF
         END DO
      END IF
C
C ADD CONTRIBUTION FROM THE A(0)S(2) TERM IF HRPA
C
      IF (SOPPA .AND. NOSIM .GT. 0) THEN
         CALL A0S2(FCONE,FC,UDV,ZYMAT,WRK(KTZYMT),H2X,WRK(KH2XP),
     *             WRK(KCOEUN),NOSIM, WRK,LWRK)
C
C  Cancel out the nonsymmetric parts of A(2) matrix. Use the explicitly
C  contructed A(2) matrix in XINDX.
C
         IF (.NOT. A2EXIST) THEN
            CALL DSCAL(N2ORBX,DP5,XINDX(KAB2),1)
            IF (IPRRSP .GT. 40) THEN
               WRITE(LUPRI,'(/A)')
     &      ' SOPPA A(2)(ai,bj) packed as .5*A(2)(i,j) and .5*A(2)(a,b)'
               CALL OUTPUT(XINDX, 1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            END IF
            A2EXIST = .TRUE.
         END IF
         CALL HRPAHM(FCONE,XINDX(KAB2),ZYMAT,NOSIM)
      ENDIF
C
C ADD CONTRIBUTION TO FCONE AND FVONE FROM ONE INDEX TRANSFORMED
C TOTAL SYMMETRIC FOCK MATRICES
C
      IF (NOSIM. GT. 0) THEN
         CALL FCKOIN(NOSIM,FC,FV,ZYMAT,FCONE,FVONE)
C        CALL FCKOIN(NSIM, FC,FV,ZYMAT,FCONE,FVONE)
C
      END IF
C
      IF ( ( NOSIM.GT.0) .AND. (KZCONF .GT. 0) .AND. .NOT.SOPPA) THEN
C
C HALF TRANSFORMED INTEGRALS ARE NOW IN H2XAC AND CONFIGURATION
C PART OF LINEAR TRANSFORMATION CAN BE CARRIED OUT FOR EACH ORBITAL
C TRIAL VECTOR
C
         CALL H2XSIG(NOSIM,FCONE,ZYMAT,WRK(KH2XAC),
     *               EVECS(1+NCSIM*KZYVAR),
     *               XINDX,WRK(KFREE),LFREE)
C
      END IF
C
C
C DISTRIBUTE FOCK AND Q MATRICES IN LINEAR TRANSFORMED VECTORS
C
      IF (NOSIM.GT.0) THEN
         IF ((IPRRSP.GT.30).AND.(NASHT.GT.1)) THEN
            DO 700 IOSIM = 1,NOSIM
               WRITE(LUPRI,'(/A,I5,A)') ' QAONE FOR ',IOSIM,
     *             ' ORBITAL TRIAL VECTOR'
               CALL OUTPUT(QAONE(1,1,IOSIM),1,NORBT,1,NASHT,
     *                     NORBT,NASHT,1,LUPRI)
               WRITE(LUPRI,'(/A,I5,A)') ' QBONE FOR',IOSIM,
     *             ' ORBITAL TRIAL VECTOR'
               CALL OUTPUT(QBONE(1,1,IOSIM),1,NORBT,1,NASHT,
     *                     NORBT,NASHT,1,LUPRI)
 700        CONTINUE
         END IF
         IF (IPRRSP.GT.50) THEN
            DO 1012 IOSIM = 1,NOSIM
               WRITE(LUPRI,'(/2A,I5)')' FCONE total contribution'
     *                              ,' FOR VECTOR ',IOSIM
               CALL OUTPUT(FCONE(1,1,IOSIM),
     *                     1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
               IF (NASHT .GT. 0) THEN
                  WRITE(LUPRI,'(/2A,I5)')' FVONE total contribution'
     *                                 ,' FOR VECTOR ',IOSIM
                  CALL OUTPUT(FVONE(1,1,IOSIM),
     *                           1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
               END IF
 1012       CONTINUE
         END IF
         CALL RSPORB(.TRUE.,NOSIM,FC,FCONE,FVONE,QAONE,QBONE,UDV,
     *               EVECS(NCSIM*KZYVAR+1),ZYMAT)
C        CALL RSPORB(ONEIND,NSIM,FC,FCX,FVX,QAX,QBX,UDV,
C    *                  EVECS,ZYMAT)
      END IF
      IF (NCSIM.GT.0) THEN
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
 1014       CONTINUE
         END IF
         IF (HSROHF.OR.DFTADD) THEN
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
            CALL RSPORB(.FALSE.,NCSIM,FC,DUMMY,FVTD,QATD,QBTD,WRK(KDVT),
     *               EVECS,DUMMY)
         END IF
      END IF
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
      RETURN
      END
