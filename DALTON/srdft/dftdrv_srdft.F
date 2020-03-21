!
!...   Copyright (c) 2014 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2014 (2014), see http://daltonprogram.org"
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


#ifdef INCLUDE_DFTDRV_SRDFT
C*****************************************************************************
      SUBROUTINE DFTEXC(DMAT,NDMAT,EXCMAT,NXCMAT,MAKDEN,DOERG,DOGRD,
     &                  DOLND,DOATR,TRPLET,HES,DOHES,WORK,LWORK,KPRINT)
C*****************************************************************************
C
C     T. Helgaker sep 99 / feb 01
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C 
      PARAMETER (FACTHR = 1.0D-8) 
      LOGICAL MAKDEN, DOERG, DOGRD, DOLND, DOATR, DOHES, TRPLET, DOGGA,
     &        DODRC, DOVWN, DOBCK, DOLYP
      DIMENSION DMAT(NBAST,NBAST,NDMAT), EXCMAT(NBAST,NBAST,NXCMAT),
     &          HES(*), WORK(LWORK)
C
#include "inforb.h"
#include "nuclei.h"
#include "dftcom.h"
#include "dftinf.h"
#include "dfterg.h"
C
      IPRINT = 6
      NBUF   = 1000000
C
C     Square DMAT if necessary
C 
      IF (MAKDEN) CALL DFTDNS(DMAT,WORK,LWORK,IPRINT)
C
      NC0  = 0
      IF (DOHES) NC0  = NORBT
C
      FACDRC = WDFTX
      FACVWN = WDFTC
      FACBCK = WDFTB
      FACLYP = WDFTL
C
      DODRC = abs(FACDRC) .GT. FACTHR
      DOVWN = abs(FACVWN) .GT. FACTHR 
      DOBCK = abs(FACBCK) .GT. FACTHR 
      DOLYP = abs(FACLYP) .GT. FACTHR 
      DOGGA = DOBCK .OR. DOLYP
C
C     Number of AOs and their addresses
C
      NDER = 0
      IF (DOGRD) NDER = 1
      IF (DOGGA) THEN
         NDER = NDER + 1 
         IF (DFTPOT) NDER = NDER + 1 
      END IF
C
      IF (NDER.EQ.0) NTYPSO =  1
      IF (NDER.EQ.1) NTYPSO =  4
      IF (NDER.EQ.2) NTYPSO = 10
      NSO0 = 1
      NSO1 = 2
      NSO2 = 5
      IF (DOLND) THEN
         NTYPSO = NTYPSO + 3
         NSOB   = NTYPSO - 2 
         IF (DOGGA) THEN
            NTYPSO = NTYPSO + 9
            NSOB1  = NTYPSO - 8 
         END IF
      END IF
      KSO0 = (NSO0-1)*NBAST + 1
      KSO1 = (NSO1-1)*NBAST + 1
      KSO2 = (NSO2-1)*NBAST + 1
      KSOB = (NSOB-1)*NBAST + 1
      KSOB1 = (NSOB1-1)*NBAST + 1
C
C     Allocations 
C
      KX    = 1
      KY    = KX    + NBUF
      KZ    = KY    + NBUF
      KW    = KZ    + NBUF
      KGSO  = KW    + NBUF
      KCNT  = KGSO  + NTYPSO*NBAST 
      KDST  = KCNT  + NBAST
      KC0   = KDST  + NATOMS 
      KC1   = KC0   + NC0
      KC2   = KC1   + 3*NC0
      KDGA  = KC2   + NC0
      KLST  = KDGA  + NBAST
      LWRK  = LWORK - KLST + 1
      IF (KLST.GT.LWORK) CALL STOPIT('DFTEXC','DFTDRV',KLST,LWORK)
#if defined (VAR_MPI)
      CALL DFTPAR(DMAT(1,1,1),DMAT(1,1,2),EXCMAT,NXCMAT,DOERG,DOGRD,
     &     DOLND,DOATR,WORK(KX),WORK(KY),WORK(KZ),WORK(KW),NBUF,
     &     WORK(KGSO),WORK(KCNT),WORK(KDST),WORK(KC0),WORK(KC1),
     &     WORK(KC2),WORK(KDGA),HES,DOHES,TRPLET,FACDRC,FACVWN,
     &     FACBCK,FACLYP,DODRC,DOVWN,DOBCK,DOLYP,DOGGA,
     &     WORK(KLST),LWRK,IPRINT)
#else 
      CALL DFTDRV(DMAT(1,1,1),DMAT(1,1,2),EXCMAT,NXCMAT,DOERG,DOGRD,
     &     DOLND,DOATR,WORK(KX),WORK(KY),WORK(KZ),WORK(KW),NBUF,
     &     WORK(KGSO),WORK(KCNT),WORK(KDST),WORK(KC0),WORK(KC1),
     &     WORK(KC2),WORK(KDGA),HES,DOHES,TRPLET,FACDRC,FACVWN,
     &     FACBCK,FACLYP,DODRC,DOVWN,DOBCK,DOLYP,DOGGA,
     &     WORK(KLST),LWRK,IPRINT)
#endif
      RETURN
      END


C*****************************************************************************
      SUBROUTINE DFTDRV(DMAT,DTRMAT,EXCMAT,NXCMAT,DOERG,DOGRD,DOLND,
     &                  DOATR,CORX,CORY,CORZ,WEIGHT,NBUF,GSO,NCNT,DST,
     &                  C0,C1,C2,DMAGAO,HES,DOHES,TRPLET,FACDRC,FACVWN,
     &                  FACBCK,FACLYP,DODRC,DOVWN,DOBCK,DOLYP,DOGGA,
     &                  WORK,LWORK,IPRINT)
C*****************************************************************************
C     T. Helgaker sep 99 / feb 01
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "pi.h"
#include "dummy.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, DP5 = 0.5D0)
C     PARAMETER (DFTHR0 = 1.0D-8, DFTHRL = 1.0D-10, DFTHRI = 1.0D-13) 
C
#include "inforb.h"
#include "infvar.h"
#include "nuclei.h"
#include "dfterg.h"
#include "energy.h"
#include "dftcom.h"
#include "dftinf.h"
#include "orgcom.h"
#include "maxorb.h"
#include "infpar.h"
C
      INTEGER A, B
      LOGICAL DOERG, DOGRD, DODRC, DOVWN, DOBCK, DOLYP, 
     &        DOLND, DOATR, FROMVX, DOGGA, TRPLET, DOHES 
      DIMENSION CORX(NBUF), CORY(NBUF), CORZ(NBUF), WEIGHT(NBUF), 
     &          DMAT(NBAST,NBAST), DTRMAT(NBAST,NBAST),
     &          GSO(NBAST*NTYPSO), EXCMAT(NBAST,NBAST,NXCMAT), 
     &          NCNT(NBAST), DST(NATOMS), C0(NORBT), C1(NORBT,3), 
     &          C2(NORBT),DMAGAO(NBAST), HES(NOCCT,NVIRT,NOCCT,NVIRT), 
     &          RHG(3), WORK(LWORK), COR(3)
      CHARACTER*16 QUADNAME
C
#include "chrnos.h"
C
      CALL TIMER('START ',TIMSTR,TIMEND)
      QUADNAME='                '
C
C
C     Calculate grid
C
      IF (.NOT.DFTGRID_DONE_OLD) THEN
         CALL MAKE_DFTGRID(WORK,LWORK,NTOT_DFTGRID,1,.FALSE.)
         CALL CONDFT 
         DFTGRID_DONE_OLD = .TRUE.
      END IF
C
      IF (IPRINT .GT. 100) THEN
         CALL HEADER('DMAT in DFTDRV ',-1)
         CALL OUTPUT(DMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         IF (DOATR) THEN
            CALL HEADER('DTRMAT in DFTDRV ',-1)
            CALL OUTPUT(DTRMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
      END IF
C
      IF (DOERG) EXCTRO = DP5*DDOT(NBAST*NBAST,DMAT,1,EXCMAT,1)
      IF (DOGRD) CALL DZERO(GRADFT,MXCOOR)
C
      FROMVX = DOGGA .AND. DFTPOT
C
      ELCTRN = D0
C
      ERGDRC = D0
      ERGVWN = D0
      ERGBCK = D0
      ERGLYP = D0
C
      ED = D0
      EV = D0
      EB = D0
      EL = D0
C
      VD  = D0
      VV  = D0
      VB1 = D0 
      VB2 = D0 
      VL1 = D0 
      VL2 = D0 
      VXB = D0
C
      LUQUAD = -1
C            
C     Make quadname: Can take 9999 procs 
C     
!     QUADNAME = 'DALTON.QUAD.'//chrnos(mynum/1000)
!    &     //chrnos((mynum-(mynum/1000)*1000)/100)
!    &     //chrnos((mynum-(mynum/100)*100)/10)
!    &     //chrnos(mynum-(mynum/10)*10)
      QUADNAME = 'DALTON.QUAD' ! GPOPEN adds slave no.
C     
      CALL GPOPEN(LUQUAD,QUADNAME,'OLD','SEQUENTIAL',
     &     'UNFORMATTED',IDUMMY,LDUMMY)      
C
      R2TOT = D0
      NPNTS = 0
  100 CONTINUE
      READ(LUQUAD) NPOINT
      IF (NPOINT.GT.0) THEN
         NPNTS = NPNTS + NPOINT
         CALL REAQUA_OLD(CORX,CORY,CORZ,WEIGHT,LUQUAD,NPOINT)
#if defined (VAR_MPI)
         DO 300 IPNT = 1+MYNUM, NPOINT, NODTOT+1 
#else 
         DO 300 IPNT = 1, NPOINT
#endif
            IF (IPRINT .GT. 100) THEN
               WRITE (LUPRI,'(2X,I6,4F15.6)') 
     &            IPNT,CORX(IPNT),CORY(IPNT),CORZ(IPNT),WEIGHT(IPNT)
            END IF
C
            WGHT = WEIGHT(IPNT)
            WDRC = FACDRC*WGHT
            WVWN = FACVWN*WGHT
            WBCK = FACBCK*WGHT
            WLYP = FACLYP*WGHT
C
C           AOs
C           ===
C
            THRINT = DFTHRI/WGHT
            COR(1) = CORX(IPNT)
            COR(2) = CORY(IPNT)
            COR(3) = CORZ(IPNT)
            CALL GETSOS(GSO,NCNT,COR,WORK,LWORK,
     &                  NBAST,DOLND,DOGGA,THRINT,IPRINT)
C
C           Density
C           =======
C
            CALL GETRHO_OLD(DMAT,GSO,RHO,RHO13,DMAGAO,THRINT,IPRINT)
C
            IF (RHO.GT.DFTHR0) THEN 
C
C              Gradient of density
C              ===================
C
               IF (DOGGA) THEN
                  CALL DGEMV('T',NBAST,3,D2,GSO(KSO1),NBAST,DMAGAO,1,
     &                       D0,RHG,1)
                  RHOGRD = DSQRT(RHG(1)**2 + RHG(2)**2 + RHG(3)**2)
               END IF
C
C              Numerical test
C              ==============
C
               IF (.FALSE.) THEN
                  CALL DFTNUM(DODRC,DOVWN,DOLYP,DOBCK,RHO,RHO13,RHOGRD)
               END IF
C
C              Hessian of density
C              ===================
C
               IF (FROMVX) THEN
                  CALL DFTRHH(DMAT,DMAGAO,GSO(KSO0),GSO(KSO1),GSO(KSO2),
     &                        RHG,RHOLAP,RHOGHG)
               END IF
C
C              Number of electrons
C              ===================
C
               ELCTRN = ELCTRN + WGHT*RHO
C
C              Energy
C              ======
C
               IF (DOERG) THEN
                  IF (DODRC) CALL EDRC(ED,RHO,RHO13)
                  IF (DOVWN) CALL EVWN(EV,RHO,RHO13)
                  IF (DOLYP) CALL ELYP(EL,RHO,RHO13,RHOGRD)
                  IF (DOBCK) CALL EBCK(EB,RHO,RHO13,RHOGRD)
                  ERGDRC = ERGDRC + WDRC*ED
                  ERGVWN = ERGVWN + WVWN*EV
                  ERGBCK = ERGBCK + WBCK*EB 
                  ERGLYP = ERGLYP + WLYP*EL
               END IF
C
C              Exchange-correlation potential 
C              ==============================
C
               IF (DOERG.OR.DOGRD.OR.DOLND) THEN
                  IF (DODRC) CALL VDRC(VD,RHO13)
                  IF (DOVWN) CALL VVWN(VV,RHO,RHO13)
                  IF (FROMVX) THEN
                     IF(DOBCK) CALL VBCK(VB1,RHO,RHO13,RHOGRD,RHOLAP,
     &                                   RHOGHG)
                     IF(DOLYP) CALL VLYP(VL1,RHO,RHO13,RHOGRD,RHOLAP)
                  ELSE
                     IF (DOBCK) CALL GBCK(VB1,VB2,RHO,RHO13,RHOGRD)
                     IF (DOLYP) CALL GLYP(VL1,VL2,RHO,RHO13,RHOGRD)
                     VXB = WBCK*VB2 + WLYP*VL2
                  END IF
                  VXC = WDRC*VD + WVWN*VV + WBCK*VB1 + WLYP*VL1
               END IF
C
C              Exchange-correlation contribution to Kohn-Sham matrix
C              =====================================================
C
               IF (DOERG) THEN
                  CALL DFTKSM(EXCMAT,GSO(KSO0),GSO(KSO1),RHG,VXC,VXB,
     &                        DOGGA,FROMVX,DFTHRL)
               END IF
C
C              Exchange-correlation contribution to molecular gradient
C              =======================================================
C
               IF (DOGRD) THEN
                   CALL DFTFRC(DMAT,GSO(KSO0),GSO(KSO1),GSO(KSO2),
     &                         VXC,VXB,RHG(1),RHG(2),RHG(3),DOGGA)
               END IF
C
C              London contributions to Kohn-Sham matrix 
C              ========================================
C
               IF (DOLND) THEN
                  CALL DFTMAG(EXCMAT,
     &                        CORX(IPNT),CORY(IPNT),CORZ(IPNT),
     &                        GSO(KSO0),GSO(KSO1),GSO(KSOB),GSO(KSOB1),
     &                        VXC,VXB,RHG,DOGGA,FROMVX)
               END IF
C
C              Hessian transformation
C              ======================
C
               IF (DOATR) THEN
                  CALL DFTLTR(JWOPSY,DTRMAT,EXCMAT,GSO(KSO0),GSO(KSO1),
     &                        C0,C1,C2,HES,RHO,RHO13,RHOGRD,RHG,
     &                        WDRC,WVWN,WBCK,WLYP,DODRC,DOVWN,DOBCK,
     &                        DOLYP,DOGGA,TRPLET,DOHES,DMAGAO)
               END IF
            END IF
C
  300    CONTINUE
C
         GO TO 100
      ELSE IF (NPOINT .EQ.0 ) THEN
         GO TO 100
      END IF
  200 CONTINUE
C
      CALL GPCLOSE(LUQUAD,'KEEP')
C
#if defined  (VAR_MPI)       
      CALL UPDTENE(ELCTRN,ERGDRC,ERGVWN,ERGBCK,ERGLYP,EXCMAT,
     &     NBAST,NXCMAT,WORK,LWORK)
C     
      IF (DOGRD) CALL UPDTGRD(WORK,LWORK)
      IF (MYNUM.NE.0) RETURN
#endif
      IF (DOERG) THEN
         EDFTY = ERGDRC + ERGVWN + EXCTRO + ERGBCK + ERGLYP
     &         - DP5*DDOT(NBAST*NBAST,DMAT,1,EXCMAT,1)
      END IF
      IF (DOERG) THEN
         DO I = 1, NBAST
         DO J = 1, I - 1 
            AVERAG = DP5*(EXCMAT(I,J,1) + EXCMAT(J,I,1))
            EXCMAT(I,J,1) = AVERAG 
            EXCMAT(J,I,1) = AVERAG 
         END DO
         END DO
      END IF 
      IF (DOATR) THEN
         IF (DOHES) THEN
            DO B = 1, NVIRT
            DO J = 1, NOCCT 
            DO A = B + 1, NVIRT
            DO I = 1, NOCCT 
               HES(I,A,J,B) = HES(J,B,I,A)
            END DO
            END DO
            END DO
            END DO
         ELSE
            DO I = 1, NBAST
            DO J = 1, I - 1 
               AVERAG = EXCMAT(I,J,1) + EXCMAT(J,I,1)
               EXCMAT(I,J,1) = AVERAG 
               EXCMAT(J,I,1) = AVERAG 
            END DO
            END DO
         END IF
      END IF 
      IF (DOLND) THEN
         DO I = 1, NBAST
         DO J = 1, I - 1 
            AVERAG = DP5*(EXCMAT(I,J,1) - EXCMAT(J,I,1))
            EXCMAT(I,J,1) =   AVERAG 
            EXCMAT(J,I,1) = - AVERAG 
         END DO
         END DO
      END IF 
C
C     Test on the number of electrons
C
      ELCTRX = FLOAT(2*NRHFT)
      ERROR  = ELCTRN - ELCTRX
      IF (ABS(ERROR) .GT. DFTELS) THEN
         WRITE (LUPRI,'(4(/2X,A,F14.6),/2X,A)')
     &   ' Number of electrons from numerical integration:',ELCTRN,
     &   ' Number of electrons from orbial occupations:   ',ELCTRX,
     &   ' Error in the number of electrons:              ',ERROR,
     &   ' Error larger than DFTELS (set input):          ',DFTELS,
     &   ' Calculation aborted.'
         CALL QUIT
     &    ('Wrong number of electrons in DFTDRV. Calculation aborted.')
      END IF
C
C     Print section
C
      IF (IPRINT .GT. 2) THEN 
         WRITE (LUPRI,'(/2X,A,F14.7,1P,D14.2,I14)')
     &      ' Number of electrons/abscissas:  ', 
     &        ELCTRN,ELCTRN-ELCTRX,NPNTS
         IF (DOERG) THEN
            WRITE (LUPRI,'(2X,A,3F14.6)') 
     &         ' Dirac and Becke exchange energy:',       
     &           ERGDRC, ERGBCK, ERGDRC+ERGBCK
             WRITE (LUPRI,'(2X,A,3F14.6)') 
     &          ' VWN and LYP correlation energy :', 
     &          ERGVWN, ERGLYP, ERGVWN + ERGLYP
            WRITE (LUPRI,'(2X,A,28X,F14.6)') 
     &          ' DFT exchange-correlation energy:', 
     &          ERGDRC + ERGBCK + ERGVWN + ERGLYP
         END IF
      END IF
      IF (IPRINT .GT. 100) THEN
         IF (DOERG) THEN
            CALL HEADER('EXCMAT in DFTDRV ',-1)
            CALL OUTPUT(EXCMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
         IF (DOLND) THEN
            CALL HEADER('EXCMAT in DFTDRV (x direction)',-1)
            CALL OUTPUT(EXCMAT(1,1,1),
     &                  1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            CALL HEADER('EXCMAT in DFTDRV (y direction)',-1)
            CALL OUTPUT(EXCMAT(1,1,2),
     &                  1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            CALL HEADER('EXCMAT in DFTDRV (z direction)',-1)
            CALL OUTPUT(EXCMAT(1,1,3),
     &                  1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
         IF (DOGRD) THEN
            CALL HEADER('GRADFT in DFTDRV ',-1)
            CALL OUTPUT(GRADFT,1,3,1,NATOMS,3,NATOMS,1,LUPRI)
         END IF
         IF (DOATR) THEN
            CALL HEADER('EXCMAT in DFTDRV ',-1)
            CALL OUTPUT(EXCMAT(1,1,1),
     &                  1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
      END IF
      CALL TIMER('DFTDRV',TIMSTR,TIMEND)
      RETURN 
      END


C*****************************************************************************
      SUBROUTINE DFTATR(FMAT,CMO,DV,ZYMAT,TRPLET,KSYMOP,WORK,LWORK)
C*****************************************************************************
C     T. Helgaker Sep 99 / Feb 01
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "dummy.h"
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0)
C
      DIMENSION CMO(*),ZYMAT(N2ORBX),DV(*),FMAT(N2ORBX),WORK(LWORK)
      LOGICAL DOPRNT, TRPLET
C
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
#include "dftcom.h"
C
      DOPRNT = .FALSE.
      IPRINT = 6
C
      KDAO = 1
      KDA1 = KDAO + N2BASX
      KFAO = KDA1 + N2BASX
      KLST = KFAO + N2BASX
      LWRK = LWORK - KLST + 1
      IF (KLST .GT. LWORK) CALL STOPIT('DFTATR','LWORK',KLST,LWORK)
C
C     AO density matrix
C
      JKEEP = JWOPSY
      JWOPSY = 1
      CALL FCKDEN(.TRUE.,(NASHT.GT.0),WORK(KDAO),WORK(KDA1),CMO,DV,
     &            WORK(KLST),LWRK)
      IF (NASHT .GT. 0) THEN
C        make total rho DTAO in WORK(KDAO)
         CALL DAXPY(N2BASX,D1,WORK(KDA1),1,WORK(KDAO),1)
      END IF
      JWOPSY = JKEEP
      IF (DOPRNT) THEN
         CALL HEADER('AO density matrix in DFTATR ',-1)
         CALL OUTPUT(WORK(KDAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C     AO transformation matrix
C
      IF (DOPRNT) THEN
         CALL HEADER('ZYMAT in DFTATR ',-1)
         CALL OUTPUT(ZYMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      if (nasht.gt.0) call quit('implement nasht.gt.0 in deq27 call')
      CALL DEQ27(CMO,ZYMAT,DUMMY,WORK(KDA1),DUMMY,WORK(KLST),LWRK)
      IF (DOPRNT) THEN
         CALL HEADER('AO transformation matrix in DFTATR ',-1)
         CALL OUTPUT(WORK(KDA1),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C     AO Fock matrix contribution
C
      CALL DZERO(WORK(KFAO),N2BASX)
      IF (DODFT) THEN
         CALL DFTEXC(WORK(KDAO),2,WORK(KFAO),1,
     &            .FALSE.,.FALSE.,.FALSE.,.FALSE.,.TRUE.,
     &            TRPLET,DUMMY,.FALSE.,WORK(KLST),LWRK,IPRINT)
      ELSE
C
C HFSRDFT --------------------------------------------------------
C        SRDFT ...
         IF (.NOT. DOHFSRDFT) CALL QUIT('DOHFSRDFT false in DFTATR')
         IF (NASHT.GT.0) CALL QUIT('NASHT .gt. 0 for SRDFT in DFTATR')
C        CALL SRDFT(ND_SIM,EXCMAT,DMAT,EDFT(1:3),DOERG,DOGRD,DOATR,WORK,LWORK,IPRINT)
!         WRITE(LUPRI,*) 'DFTATR: calling SRDFT'  !JT
         CALL SRDFT(1,WORK(KFAO),WORK(KDAO),VDUMMY,
     &              .FALSE.,.FALSE.,.TRUE.,.FALSE.,
     &              WORK(KLST),LWRK,IPRINT)
C-------------------------------------------------------------------
      END IF

      IF (DOPRNT) THEN
         CALL HEADER('AO Fock matrix contribution in DFTATR ',-1)
         CALL OUTPUT(WORK(KFAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C     MO Fock matrix contribution 
C
      CALL DZERO(WORK(KDA1),N2ORBX)
      DO ISYM = 1, NSYM
         JSYM  = MULD2H(ISYM,KSYMOP)
         NORBI = NORB(ISYM)
         NORBJ = NORB(JSYM)
         IF (NORBI.GT.0 .AND. NORBJ.GT.0) THEN 
            CALL AUTPV(ISYM,JSYM,CMO(ICMO(ISYM)+1),CMO(ICMO(JSYM)+1),
     &                 WORK(KFAO),NBAS,NBAST,WORK(KDA1),NORB,
     &                 NORBT,WORK(KLST),LWRK)
         END IF
      END DO
      IF (DOPRNT) THEN
         CALL HEADER('MO Fock matrix contribution in DFTATR ',-1)
         CALL OUTPUT(WORK(KDA1),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     Add to Fock matrix
C
      IF (DOPRNT) THEN
         CALL HEADER('Original MO Fock matrix in DFTATR ',-1)
         CALL OUTPUT(FMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      CALL DAXPY(NORBT*NORBT,D2,WORK(KDA1),1,FMAT,1) 
      IF (DOPRNT) THEN
         CALL HEADER('MO Fock matrix in DFTATR ',-1)
         CALL OUTPUT(FMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      RETURN
      END


C*****************************************************************************
      SUBROUTINE DFTATX(FMAT,CMO,ZYMAT,TRPLET,WORK,LWORK)
C*****************************************************************************
C
C     T. Helgaker Sep 99
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "dummy.h"
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0)
      PARAMETER (NOC = 7, NVT = 1)
C
      INTEGER A, B
      DIMENSION CMO(*),ZYMAT(NORBT,NORBT),FMAT(NORBT,NORBT),WORK(LWORK)
      DIMENSION HES(NOC*NVT*NOC*NVT) 
      LOGICAL DOPRNT,FIRST,TRPLET
      SAVE FIRST, HES
      DATA FIRST/.TRUE./
C
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
C
      DOPRNT = .FALSE.
      IPRINT = 2
C
      NDIM1 = (NOCCT*NVIRT)**2
      NDIM2 = (NOC*NVT)**2
      IF (NDIM1 .GT. NDIM2) CALL STOPIT('DFTATX','HES',NDIM1,NDIM2)
C
      KDAO = 1
      KDA1 = KDAO + N2BASX
      KFAO = KDA1 + N2BASX
      KLST = KFAO + N2BASX
      IF (KLST .GT. LWORK) CALL STOPIT('DFTATX','LWORK',KLST,LWORK)
      LWRK = LWORK - KLST + 1
C
C     AO density matrix
C
      IF (FIRST) THEN
         FIRST = .FALSE.
         CALL FCKDEN(.TRUE.,.FALSE.,WORK(KDAO),DUMMY,CMO,DUMMY,
     &               WORK(KLST),LWRK)
         CALL DZERO(HES,NOCCT*NVIRT*NOCCT*NVIRT)
         CALL DCOPY(NBAST*NORBT,CMO,1,WORK(KDA1),1)
         IF (DOPRNT) THEN
            CALL HEADER('AO density matrix in DFTATX ',-1)
            CALL OUTPUT(WORK(KDAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            CALL HEADER('MO coefficient matrix in DFTATX ',-1)
            CALL OUTPUT(WORK(KDA1),1,NBAST,1,NORBT,NBAST,NORBT,1,LUPRI)
         END IF
         CALL DFTEXC(WORK(KDAO),2,WORK(KFAO),1,
     &               .FALSE.,.FALSE.,.FALSE.,.FALSE.,.TRUE.,
     &               TRPLET,HES,.TRUE.,WORK(KLST),LWRK,IPRINT)
         IF (DOPRNT) THEN
            NPAIR = NOCCT*NVIRT
            CALL HEADER('Hessian matrix in DFTATX ',-1)
            CALL OUTPUT(HES,1,NPAIR,1,NPAIR,NPAIR,NPAIR,1,LUPRI)
         END IF
      END IF
      IF (DOPRNT) THEN
         CALL HEADER('ZYMAT in DFTATX ',-1)
         CALL OUTPUT(ZYMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      IHES = 0
      DO A = 1, NVIRT
      DO I = 1, NOCCT
         SUM = D0 
         DO B = 1, NVIRT
         DO J = 1, NOCCT
            IHES = IHES + 1
            SUM = SUM + HES(IHES)*(ZYMAT(J,B+NOCCT) - ZYMAT(B+NOCCT,J))
         END DO
         END DO
         WORK(KDA1 + NORBT*(NOCCT + A - 1) + I - 1) = SUM
         WORK(KDA1 + NORBT*(I - 1) + NOCCT + A - 1) = SUM 
      END DO
      END DO
      IF (DOPRNT) THEN
         CALL HEADER('MO Fock matrix contribution in DFTATX ',-1)
         CALL OUTPUT(WORK(KDA1),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     Add to Fock matrix
C
      IF (DOPRNT) THEN
         CALL HEADER('Original MO Fock matrix in DFTATX ',-1)
         CALL OUTPUT(FMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      CALL DAXPY(NORBT*NORBT,D2,WORK(KDA1),1,FMAT,1) 
      IF (DOPRNT) THEN
         CALL HEADER('MO Fock matrix in DFTATX ',-1)
         CALL OUTPUT(FMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      RETURN
      END


C*****************************************************************************
      SUBROUTINE DFTLND(FX,FY,FZ,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     T. Helgaker oct 99
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      DIMENSION FX(N2BASX), FY(N2BASX), FZ(N2BASX), WORK(LWORK)
#include "inforb.h"
C
      KXCMAT = 1
      KDAO   = KXCMAT + 3*N2BASX
      KLAST  = KDAO   +   N2BASX
      IF (KLAST .GT. LWORK) CALL STOPIT('DFTLND',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL DZERO(WORK(KXCMAT),3*N2BASX)
      CALL DFTEXC(WORK(KDAO),1,WORK(KXCMAT),3,.TRUE.,.FALSE.,.FALSE.,
     &            .TRUE.,.FALSE.,.FALSE.,
     &            DUMMY,.FALSE.,WORK(KLAST),LWRK,IPRINT)
C
      KXC = KXCMAT - 1
      KYC = KXCMAT - 1 +   N2BASX
      KZC = KXCMAT - 1 + 2*N2BASX
      DO I = 1, N2BASX
         FX(I) = FX(I) + WORK(KXC + I)
         FY(I) = FY(I) + WORK(KYC + I)
         FZ(I) = FZ(I) + WORK(KZC + I)
      END DO
C
      RETURN
      END


C*****************************************************************************
      SUBROUTINE DFTSOL(WORK,LWORK)
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      DIMENSION WORK(*)
#include "infvar.h"
#include "inflin.h"
#include "inforb.h"
C
      KVEC = 1
      KFCK = KVEC + NWOPPT
      KLST = KFCK + NNORBT - 1
      IF (KLST .GT. LWORK) CALL STOPIT('DFTSOL',' ',KLST,LWORK)
      CALL DFTSO1(WORK(KVEC),WORK(KFCK))
      RETURN
      END


C*****************************************************************************
      SUBROUTINE DFTSO1(VEC,FCK)
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxorb.h"
      PARAMETER (FAC = -1.0D0/4.0D0)
      DIMENSION VEC(*), FCK(*)
#include "inftap.h"
#include "infvar.h"
#include "inflin.h"
#include "inforb.h"
      IA(I) = JWOP(1,I)
      IB(I) = JWOP(2,I)
      IAA(I) = IA(I)*(IA(I)+1)/2
      IBB(I) = IB(I)*(IB(I)+1)/2
C
      REWIND (LUSIFC)
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC) 
      READ (LUSIFC)
      READ (LUSIFC) 
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC)
      CALL READI(LUSIFC,IRAT*NNORBT,FCK)
      DO 100 I = 1,NWOPPT
         VEC(I) = FAC*VEC(I)/(FCK(IAA(I)) - FCK(IBB(I)))
  100 CONTINUE
      IF (.FALSE.) THEN
         WRITE(LUPRI,'(//,A)')' Orbital part of solution'
         PRFAC = 0.1D0
         CALL PRKAP(NWOPPT,VEC,PRFAC,LUPRI)
      END IF
      RETURN
      END
#else
      subroutine dummy_dftdrv_srdft ! to avoid warnings about no symbols
      end
#endif /* INCLUDE_DFTDRV_SRDFT */
