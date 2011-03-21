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
#ifdef OLD_LOG
===========================================================================
Revision 1.2  2001/01/17 10:41:29  vebjornb
Calls to *MPA*B* in arhpack.F have been replaced with DGEMM calls

Revision 1.3  2000/05/24 19:00:24  hjj
new GETREF and CREF allocation (fix of triplet problem for CSF:s)
Do not read S(2) CI diagonal but use ilicit that itit is D1

961011-hjaaj
merged with SOPPA changes by ekd
LRST  : modification needed for proper handling of cubic response
950705-hjaaj
RSPPP : New code for GOS oscillator strengths.
940523-hjaaj
RSPCTL: removed CALL FLUSH(LURSP1) (system dependent); instead REWIND LURSP1
941117-hjaaj
QONEMU: changed a SMPATB to SMPAB to change DDOT to DAXPY
FCKOIN: optimized
940708-hjaaj
c: SOPPA changes: test KOFFTY.eq.1, not KSYMOP.eq.1 for CREF
LRST  : only orb (ph) trial vectors for SOPPA (no conf., i.e. no 2p2h)
RSPDIA: skip S[2] reads if SHIFT .eq. D0
931119-hjaaj
RSPCTL: improved output when insufficient memory stop
931007-hjaaj
removed RSPMC and GETREF (now in rspmai.u)
extracted **FILENAME=ORPCTL.UPD module as rsporp.u
931004-pj+hjaaj
merged previous updates (including some of the comments below and more),
merged/inserted solvent and RSPSUP changes
920613-ov: In rspe2c.u:
           A lot of explicit DOUBLE PRECISION changed to REAL for CRAY
           Statement functions DABS(X)=ABS(X) etc.
920624-ov: New conversion factors in RSPPP from CODATA
921125-hj: RSPMC  (rspe2c.u) consistency check with NCONRF, LSYMRF
           instead of NCONF,LSYM (because of problems with calls from
           ABACUS).  The NWOPT test has been disabled.
920722-hjaaj
extracted EISPACK RGG routines to file "rggeis.u"
RSPNEX: new THRLDV definition for lin.dep.
LRST: stop if RSPSUP true and restart (need some additional coding)
920721-Hinne Hettema
RSPCTL: "NSIM=0" for RSPRED before GO TO 1000 (2 places)
ORBDIA: inserted averaging code (if RSPSUP)
===========================================================================
#endif
C  /* Deck e2sym */
      SUBROUTINE E2SYM(IBTYP,CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,
     &                 WRK,LWRK)
C
#include <implicit.h>
C
C  THIS IS A TEST ROUTINE THAT
C  1)  READ LINEAR TRANSFORMED VECTORS FOR (Z,Y) VECTORS FROM LURSP3 AND
C      LURSP5 AND CALCULATE LINEAR TRANSFORMATION FOR (Y,Z) VECTORS
C      USING RSPLIN
C
C  2)  SET UP REDUCED E(2) AND S(2) AS SQUARE MATRICES
C      AND CHECK THE STRUCTURE
C
#include <priunit.h>
#include <wrkrsp.h>
#include <inftap.h>
#include <infrsp.h>
#include <infpri.h>
C
      PARAMETER ( D0 = 0.0D0 , DM1=-1.0D0 , D1 = 1.0D0 , DTOL = 1.0D-6 )
      PARAMETER ( DUMMY = 1.0D20 )
      DIMENSION IBTYP(*)
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
#include <ibndxdef.h>
C
      CALL HEADER('Test output requested with .ABSYM',-1)
C
C DETERMINE NUMBER OF ORBITAL AND CSF TRIAL VECTORS
C
      KZDIM = KZRED + KOFFTY
      KZYDIM= 2*KZDIM
      NCTOT = 0
      NOTOT = 0
      DO 1000 J = 1,KZDIM
         IF (IBTYP(J).EQ.JBCNDX) THEN
            NCTOT = NCTOT + 1
         ELSE
            NOTOT = NOTOT + 1
         END IF
 1000 CONTINUE
C
C READ IN TRIAL (Z,Y) TRIAL VECTORS AND LINEAR TRANSFORMED VECTORS
C
      KBCVEC= 1
      KBOVEC= KBCVEC + 2*NCTOT*KZYCON
      KE2LIN= KBOVEC + 2*NOTOT*KZYWOP
      KS2LIN= KE2LIN + KZYVAR*KZYDIM
      KE2MAT= KS2LIN + KZYVAR*KZYDIM
      KS2MAT= KE2MAT + KZYDIM*KZYDIM
      KDIAE = KS2MAT + KZYDIM*KZYDIM
      IF (SOPPA) THEN
         KWRK1 = KDIAE + KZCONF
      ELSE
         KWRK1 = KDIAE
      ENDIF
      LWRK1 = LWRK   - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('E2SYM',KWRK1-1,LWRK)
C
C
      CALL REWSPL(LURSP3)
      CALL REWSPL(LURSP5)
#if defined (VAR_SPLITFILES)
      JRSP3 = 0
      JRSP5 = 0
#endif
      JBOOFF = KBOVEC
      JBCOFF = KBCVEC
      JE2OFF = KE2LIN
      JS2OFF = KS2LIN
      JOTOT  = NOTOT*KZYWOP
      JCTOT  = NCTOT*KZYCON
      IF (SOPPA) THEN
         JDIAE = KDIAE
         CALL REWSPL(LURSP4)
#if defined (VAR_SPLITFILES)
         JRSP4 = 0
         CALL READST(LURSP4,LRSP4,WRK(JDIAE),JRSP4,.TRUE.)
         IF (LRSP4 .NE. KZCONF) CALL QUIT('Inconsistent record length'//
     &        ' on unit LURSP4 in E2SYM')
#else
         CALL READT(LURSP4,KZCONF,WRK(JDIAE))
#endif
      ENDIF
      CALL DZERO(WRK(KBCVEC),2*NCTOT*KZYCON)
      DO 15 J = 1,KZDIM
         IF (IBTYP(J).EQ.JBONDX) THEN
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP3,LRSP3,WRK(JBOOFF),JRSP3,.TRUE.)
            IF (LRSP3 .NE. KZYWOP) CALL QUIT('Inconsistent record '//
     &           'length on unit LURSP3 in E2SYM')
#else
            CALL READT(LURSP3,KZYWOP,WRK(JBOOFF))
#endif
            CALL DCOPY(KZWOPT,WRK(JBOOFF),1,WRK(JBOOFF+JOTOT+KZWOPT),1)
            CALL DCOPY(KZWOPT,WRK(JBOOFF+KZWOPT),1,WRK(JBOOFF+JOTOT),1)
            CALL DZERO(WRK(JS2OFF),KZYVAR)
            IF (KZWOPT.GT.0)
     *      CALL RSPSLI(0,1,WRK(JBOOFF),WRK(JBOOFF),UDV,WRK(JS2OFF),
     *                  XINDX,WRK(KWRK1),LWRK1)
            JBOOFF = JBOOFF + KZYWOP
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP5,LRSP5,WRK(JE2OFF),JRSP5,.TRUE.)
            IF (LRSP5 .NE. KZYVAR) CALL QUIT('Inconsistent record '//
     &           'length on unit LURSP5 in E2SYM')
#else
            CALL READT(LURSP5,KZYVAR,WRK(JE2OFF))
#endif
         ELSE
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP3,LRSP3,WRK(JBCOFF),JRSP3,.TRUE.)
            IF (LRSP3 .NE. KZCONF) CALL QUIT('Inconsistent record '//
     &           'length on unit LURSP3 in E2SYM (2)')
#else
            CALL READT(LURSP3,KZCONF,WRK(JBCOFF))
#endif
            CALL DCOPY(KZCONF,WRK(JBCOFF),1,WRK(JBCOFF+JCTOT+KZCONF),1)
            CALL DZERO(WRK(JS2OFF),KZYVAR)
            IF (KZCONF.GT.0)
     *      CALL RSPSLI(1,0,WRK(JBCOFF),WRK(JBCOFF),UDV,WRK(JS2OFF),
     *                  XINDX,WRK(KWRK1),LWRK1)
            IF (SOPPA) THEN
#if defined (VAR_SPLITFILES)
               CALL READST(LURSP5,LRSP5,WRK(JE2OFF+KZCONF),JRSP5,.TRUE.)
               IF (LRSP5 .NE. KZYWOP) CALL QUIT('Inconsistent record '//
     &           'length on unit LURSP5 in E2SYM (2)')
#else
               CALL READT(LURSP5,KZYWOP,WRK(JE2OFF+KZCONF))
#endif
               CALL DCOPY(KZWOPT,WRK(JE2OFF+KZVAR),1,
     *                    WRK(JE2OFF+KZVAR+KZCONF),1)
               DO I=0,KZCONF-1
                  WRK(JE2OFF+I) = WRK(JDIAE+I) * WRK(JBCOFF+I)
               ENDDO
               CALL DZERO(WRK(JE2OFF+KZVAR),KZCONF)
            ELSE
#if defined (VAR_SPLITFILES)
               CALL READST(LURSP5,LRSP5,WRK(JE2OFF),JRSP5,.TRUE.)
               IF (LRSP5 .NE. KZYVAR) CALL QUIT('Inconsistent record '//
     &           'length on unit LURSP5 in E2SYM (2)')
#else
               CALL READT(LURSP5,KZYVAR,WRK(JE2OFF))
#endif
            ENDIF
            CALL DCOPY(KZVAR,WRK(JE2OFF+KZVAR),1,
     *                 WRK(JE2OFF+KZDIM*KZYVAR),1)
            CALL DCOPY(KZVAR,WRK(JE2OFF),1,
     *                 WRK(JE2OFF+KZDIM*KZYVAR+KZVAR),1)
            CALL DCOPY(KZVAR,WRK(JS2OFF+KZVAR),1,
     *                 WRK(JS2OFF+KZDIM*KZYVAR),1)
            CALL DCOPY(KZVAR,WRK(JS2OFF),1,
     *                 WRK(JS2OFF+KZDIM*KZYVAR+KZVAR),1)
            CALL DSCAL(KZYVAR,DM1,WRK(JS2OFF+KZDIM*KZYVAR),1)
            JBCOFF = JBCOFF + KZYCON
         ENDIF
         JE2OFF = JE2OFF + KZYVAR
         JS2OFF = JS2OFF + KZYVAR
 15   CONTINUE
      IF (IPRRSP.GT.110) THEN
         WRITE(LUPRI,'(2A,I8)')' E(2) LINEAR TRANSFORMED VECTORS',
     *   ' : DIMENSION ',KZDIM
         CALL OUTPUT(WRK(KE2LIN),1,KZYVAR,1,KZDIM,KZYVAR,KZDIM,1,LUPRI)
         WRITE(LUPRI,'(2A,I8)')' S(2) LINEAR TRANSFORMED VECTORS',
     *   ' : DIMENSION ',KZDIM
         CALL OUTPUT(WRK(KS2LIN),1,KZYVAR,1,KZDIM,KZYVAR,KZDIM,1,LUPRI)
         IF ( NOTOT.GT.0 ) THEN
            WRITE(LUPRI,'(2A,I8)')' ORBITAL TRIAL VECTORS ',
     *      ' : DIMENSION ',2*NOTOT
            CALL OUTPUT(WRK(KBOVEC),1,KZYWOP,1,2*NOTOT,KZYWOP,
     *                                              2*NOTOT,1,LUPRI)
         END IF
         IF ( NCTOT.GT.0 ) THEN
            WRITE(LUPRI,'(2A,I8)')' CONFIGURATION TRIAL VECTORS ',
     *      ' : DIMENSION ',2*NCTOT
            CALL OUTPUT(WRK(KBCVEC),1,KZYCON,1,2*NCTOT,KZYCON,
     *                                              2*NCTOT,1,LUPRI)
         END IF
      END IF
      JBOOFF = KBOVEC + JOTOT
      JBCOFF = KBCVEC + JCTOT
C
C CARRY OUT LINEAR TRANSFORMATION ON (Y,Z) VECTORS
C
      DO 20 J = 1,KZDIM
         IF (IBTYP(J).EQ.JBONDX) THEN
            NCSIM = 0
            NOSIM = 1
            KBVEC = JBOOFF
            JBOOFF= JBOOFF + KZYWOP
            CALL RSPLIN(NCSIM,NOSIM,WRK(KBVEC),WRK(KBVEC),
     *                  CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *                  XINDX,WRK(KWRK1),LWRK1)
C
C           CALL RSPLIN(NCSIM,NOSIM,ZYCVEC,ZYOVEC,
C    *                  CMO,UDV,PV,FC,FV,FCAC,H2AC,
C    *                  XINDX,WRK,LWRK)
C
C
            CALL DCOPY(KZYVAR,WRK(KWRK1),1,
     *              WRK(KE2LIN+(J-1+KZDIM)*KZYVAR),1)
            CALL DCOPY(KZYVAR,WRK(KWRK1+KZYVAR),1,
     *              WRK(KS2LIN+(J-1+KZDIM)*KZYVAR),1)
         END IF
 20   CONTINUE
      IF (IPRRSP.GT.110) THEN
         WRITE(LUPRI,'(2A,I8)')' E(2) LINEAR TRANSFORMED VECTORS',
     *   ' : DIMENSION ',KZYDIM
         CALL OUTPUT(WRK(KE2LIN),1,KZYVAR,1,KZYDIM,KZYVAR,KZYDIM,1,
     &               LUPRI)
         WRITE(LUPRI,'(2A,I8)')' S(2) LINEAR TRANSFORMED VECTORS',
     *   ' : DIMENSION ',KZYDIM
         CALL OUTPUT(WRK(KS2LIN),1,KZYVAR,1,KZYDIM,KZYVAR,KZYDIM,1,
     &               LUPRI)
         IF ( NOTOT.GT.0 ) THEN
            WRITE(LUPRI,'(2A,I8)')' ORBITAL TRIAL VECTORS ',
     *      ' : DIMENSION ',2*NOTOT
            CALL OUTPUT(WRK(KBOVEC),1,KZYWOP,1,2*NOTOT,KZYWOP,
     *                                              2*NOTOT,1,LUPRI)
         END IF
         IF ( NCTOT.GT.0 ) THEN
            WRITE(LUPRI,'(2A,I8)')' CONFIGURATION TRIAL VECTORS ',
     *      ' : DIMENSION ',2*NCTOT
            CALL OUTPUT(WRK(KBCVEC),1,KZYCON,1,2*NCTOT,KZYCON,
     *                                              2*NCTOT,1,LUPRI)
         END IF
      END IF
C
C CHECK IF TRIAL VECTORS ARE ORTHONORMAL
C
      IF ( NOTOT.GT.0) THEN
         IJ = 0
         XMAX = D0
         IMAX = 0
         JMAX = 0
         DO 30 I = 1,2*NOTOT
            IOFF = KBOVEC+(I-1)*KZYWOP
            DO 40 J = 1,I
               IJ = IJ + 1
               JOFF = KBOVEC+(J-1)*KZYWOP
               XIJ = DDOT(KZYWOP,WRK(IOFF),1,WRK(JOFF),1)
               WRK(KWRK1-1+IJ) = XIJ
               IF ((ABS(XIJ).GT.XMAX).AND.(I.NE.J)) THEN
                  IMAX = I
                  JMAX = J
                  XMAX = XIJ
               ENDIF
 40         CONTINUE
 30      CONTINUE
         WRITE(LUPRI,'(/A)')' OVERLAP FOR ORBITAL TRIAL VECTORS'
         CALL OUTPAK(WRK(KWRK1),2*NOTOT,1,LUPRI)
         IF ( XMAX.GT.DTOL )
     *   WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *   ' LARGEST ELEMENT (I,J), I=',IMAX,' J=',JMAX,' :',XMAX
      END IF
      IF ( NCTOT.GT.0 ) THEN
         IJ = 0
         XMAX = D0
         IMAX = 0
         JMAX = 0
         DO 50 I = 1,2*NCTOT
            IOFF = KBCVEC+(I-1)*KZYCON
            DO 60 J = 1,I
               IJ = IJ + 1
               JOFF = KBCVEC+(J-1)*KZYCON
               XIJ = DDOT(KZYCON,WRK(IOFF),1,WRK(JOFF),1)
               WRK(KWRK1-1+IJ) = XIJ
               IF ((ABS(XIJ).GT.XMAX).AND.(I.NE.J)) THEN
                  IMAX = I
                  JMAX = J
                  XMAX = XIJ
               ENDIF
 60         CONTINUE
 50      CONTINUE
         WRITE(LUPRI,'(/A)')' OVERLAP FOR CONFIGURATION TRIAL VECTORS'
         CALL OUTPAK(WRK(KWRK1),2*NCTOT,1,LUPRI)
         IF ( XMAX.GT.DTOL )
     *   WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *   ' LARGEST ELEMENT (I,J), I=',IMAX,' J=',JMAX,' :',XMAX
      END IF
C
C SET UP REDUCED E2 AND S2 MATRICES
C
      JBOVEC = KBOVEC
      JBCVEC = KBCVEC
      DO 70 J = 1,KZYDIM
         JCHK = J
         IF ( J.GT.KZDIM ) THEN
            JCHK = JCHK - KZDIM
         END IF
         IF ( IBTYP(JCHK).EQ.JBCNDX) THEN
            JZYVAR = KZYCON
            JZVAR  = KZCONF
            JZOFF  = 0
            JYOFF  = KZVAR
            JBVEC  = JBCVEC
            JBCVEC = JBCVEC + KZYCON
         ELSE
            JZYVAR = KZYWOP
            JZVAR  = KZWOPT
            JZOFF  = KZCONF
            JYOFF  = KZVAR  + KZCONF
            JBVEC  = JBOVEC
            JBOVEC = JBOVEC + KZYWOP
         END IF
         DO 80 I = 1,KZYDIM
            X1 = DDOT(JZVAR,WRK(JBVEC),1,
     *                WRK(KE2LIN+(I-1)*KZYVAR+JZOFF),1)
            X2 = DDOT(JZVAR,WRK(JBVEC+JZVAR),1,
     *                WRK(KE2LIN+(I-1)*KZYVAR+JYOFF),1)
            WRK(KE2MAT-1+(I-1)*KZYDIM+J) = X1 + X2
            X1 = DDOT(JZVAR,WRK(JBVEC),1,
     *                WRK(KS2LIN+(I-1)*KZYVAR+JZOFF),1)
            X2 = DDOT(JZVAR,WRK(JBVEC+JZVAR),1,
     *                WRK(KS2LIN+(I-1)*KZYVAR+JYOFF),1)
            WRK(KS2MAT-1+(I-1)*KZYDIM+J) = X1 + X2
 80      CONTINUE
 70   CONTINUE
      WRITE(LUPRI,'(A,I8)')' REDUCED E(2) MATRIX : DIMENSION ',KZYDIM
      CALL OUTPUT(WRK(KE2MAT),1,KZYDIM,1,KZYDIM,KZYDIM,KZYDIM,1,LUPRI)
      WRITE(LUPRI,'(A,I8)')' REDUCED S(2) MATRIX : DIMENSION ',KZYDIM
      CALL OUTPUT(WRK(KS2MAT),1,KZYDIM,1,KZYDIM,KZYDIM,KZYDIM,1,LUPRI)
C
C CHECK BLOCK STRUCTURE OF E2 AND S2 AND WRITE OUT LARGEST DEVIATIONS
C
      ZAMAX = D0
      IA    = 0
      JA    = 0
      ZBMAX = D0
      IB    = 0
      JB    = 0
      ZSMAX = D0
      IS    = 0
      JS    = 0
      ZDMAX = D0
      ID    = 0
      JD    = 0
      DO 150 IZ = 1,KZDIM
         IY = IZ +KZDIM
         DO 160 JZ = 1,IZ
            JY = JZ + KZDIM
            ZAIJ  = WRK(KE2MAT-1+(IZ-1)*KZYDIM+JZ)
            ZAJI  = WRK(KE2MAT-1+(JZ-1)*KZYDIM+IZ)
            YAIJ  = WRK(KE2MAT-1+(IY-1)*KZYDIM+JY)
            YAJI  = WRK(KE2MAT-1+(JY-1)*KZYDIM+IY)
            ZDEV  = MAX(ABS(ZAIJ-ZAJI),ABS(ZAIJ-YAIJ))
            ZDEV  = MAX(ZDEV,ABS(ZAIJ-YAJI))
            IF (ZDEV.GT.ZAMAX) THEN
               ZAMAX = ZDEV
               IA    = IZ
               JA    = JZ
            END IF
            ZBIJ  = WRK(KE2MAT-1+(IZ-1)*KZYDIM+JY)
            ZBJI  = WRK(KE2MAT-1+(JZ-1)*KZYDIM+IY)
            YBIJ  = WRK(KE2MAT-1+(IY-1)*KZYDIM+JZ)
            YBJI  = WRK(KE2MAT-1+(JY-1)*KZYDIM+IZ)
            ZDEV  = MAX(ABS(ZBIJ-ZBJI),ABS(ZBIJ-YBIJ))
            ZDEV  = MAX(ZDEV,ABS(ZBIJ-YBJI))
            IF (ZDEV.GT.ZBMAX) THEN
               ZBMAX = ZDEV
               IB    = IZ
               JB    = JZ
            END IF
            ZSIJ  = WRK(KS2MAT-1+(IZ-1)*KZYDIM+JZ)
            ZSJI  = WRK(KS2MAT-1+(JZ-1)*KZYDIM+IZ)
            YSIJ  =-WRK(KS2MAT-1+(IY-1)*KZYDIM+JY)
            YSJI  =-WRK(KS2MAT-1+(JY-1)*KZYDIM+IY)
            ZDEV  = MAX(ABS(ZSIJ-ZSJI),ABS(ZSIJ-YSIJ))
            ZDEV  = MAX(ZDEV,ABS(ZSIJ-YSJI))
            IF (ZDEV.GT.ZSMAX) THEN
               ZSMAX = ZDEV
               IS    = IZ
               JS    = JZ
            END IF
            ZDIJ  = WRK(KS2MAT-1+(IZ-1)*KZYDIM+JY)
            ZDJI  =-WRK(KS2MAT-1+(JZ-1)*KZYDIM+IY)
            YDIJ  =-WRK(KS2MAT-1+(IY-1)*KZYDIM+JZ)
            YDJI  = WRK(KS2MAT-1+(JY-1)*KZYDIM+IZ)
            ZDEV  = MAX(ABS(ZDIJ-ZDJI),ABS(ZDIJ-YDIJ))
            ZDEV  = MAX(ZDEV,ABS(ZDIJ-YDJI))
            IF (ZDEV.GT.ZDMAX) THEN
               ZDMAX = ZDEV
               ID    = IZ
               JD    = JZ
            END IF
 160     CONTINUE
 150  CONTINUE
      IZ=IA
      JZ=JA
      IY=IZ+KZDIM
      JY=JZ+KZDIM
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' A(I,J)-MATRIX : MAXIMUM DEVIATION, I=',IA,' J=',JA,' DEV=',ZAMAX
      IF ( ZAMAX.GT.DTOL )
     *WRITE(LUPRI,'(/4(/A,I5,A,I5,A,1P,G16.8))')
     *' IZ=',IZ,' JZ=',JZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(IZ-1)*KZYDIM+JZ),
     *' IZ=',JZ,' JZ=',IZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(JZ-1)*KZYDIM+IZ),
     *' IZ=',IY,' JZ=',JY,' E2(IZ,JZ)',WRK(KE2MAT-1+(IY-1)*KZYDIM+JY),
     *' IZ=',JY,' JZ=',IY,' E2(IZ,JZ)',WRK(KE2MAT-1+(JY-1)*KZYDIM+IY)
      IZ=IB
      JZ=JB
      IY=IZ+KZDIM
      JY=JZ+KZDIM
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' B(I,J)-MATRIX : MAXIMUM DEVIATION, I=',IB,' J=',JB,' DEV=',ZBMAX
      IF ( ZBMAX.GT.DTOL )
     *WRITE(LUPRI,'(/4(/A,I5,A,I5,A,1P,G16.8))')
     *' IZ=',IZ,' JZ=',JY,' E2(IZ,JZ)',WRK(KE2MAT-1+(IZ-1)*KZYDIM+JY),
     *' IZ=',JZ,' JZ=',IY,' E2(IZ,JZ)',WRK(KE2MAT-1+(JZ-1)*KZYDIM+IY),
     *' IZ=',IY,' JZ=',JZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(IY-1)*KZYDIM+JZ),
     *' IZ=',JY,' JZ=',IZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(JY-1)*KZYDIM+IZ)
      IZ=IS
      JZ=JS
      IY=IZ+KZDIM
      JY=JZ+KZDIM
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' S(I,J)-MATRIX : MAXIMUM DEVIATION, I=',IS,' J=',JS,' DEV=',ZSMAX
      IF ( ZSMAX.GT.DTOL )
     *WRITE(LUPRI,'(/4(/A,I5,A,I5,A,1P,G16.8))')
     *' IZ=',IZ,' JZ=',JY,' E2(IZ,JZ)',WRK(KE2MAT-1+(IZ-1)*KZYDIM+JY),
     *' IZ=',JZ,' JZ=',IY,' E2(IZ,JZ)',WRK(KE2MAT-1+(JZ-1)*KZYDIM+IY),
     *' IZ=',IY,' JZ=',JZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(IY-1)*KZYDIM+JZ),
     *' IZ=',JY,' JZ=',IZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(JY-1)*KZYDIM+IZ)
      IZ=ID
      JZ=JD
      IY=IZ+KZDIM
      JY=JZ+KZDIM
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' D(I,J)-MATRIX : MAXIMUM DEVIATION, I=',IA,' J=',JA,' DEV=',ZDMAX
      IF ( ZDMAX.GT.DTOL )
     *WRITE(LUPRI,'(/4(/A,I5,A,I5,A,1P,G16.8))')
     *' IZ=',IZ,' JZ=',JY,' E2(IZ,JZ)',WRK(KE2MAT-1+(IZ-1)*KZYDIM+JY),
     *' IZ=',JZ,' JZ=',IY,' E2(IZ,JZ)',WRK(KE2MAT-1+(JZ-1)*KZYDIM+IY),
     *' IZ=',IY,' JZ=',JZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(IY-1)*KZYDIM+JZ),
     *' IZ=',JY,' JZ=',IZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(JY-1)*KZYDIM+IZ)
C
C CHECK EIGENVALUES OF E2
C
C ALLOCATE WORK SPACE
C
      KE2    = KWRK1
      KEIVAL = KE2   + (KZYDIM*KZYDIM+KZYDIM)/2
      KEIVEC = KEIVAL+ KZYDIM
      KBIG   = KEIVEC+ KZYDIM*KZYDIM
      KJBIG  = KBIG  + KZYDIM
      KWRK1  = KJBIG + KZYDIM
C
      WRITE(LUPRI,'(/A,/A)')' ********************************',
     *' CHECK EIGENVALUES OF E2 WITH JACO ROUTINE '
      IF ( KOFFTY.EQ.1 ) WRITE(LUPRI,'(/A,/A)')
     *' *** E2 *** HAVE TWO ZERO EIGENVALUES BECAUSE REFERENCE STATE',
     *'            IS INCLUDED IN REDUCED SPACE'
      IJ = 0
      DO 500 J = 1,KZYDIM
         DO 600 I = 1,J
            IJ = IJ + 1
            WRK(KE2-1+IJ) = WRK(KE2MAT-1+(I-1)*KZYDIM+J)
 600     CONTINUE
 500  CONTINUE
      CALL JACO(WRK(KE2),WRK(KEIVEC),KZYDIM,KZYDIM,0,
     *          WRK(KBIG),WRK(KJBIG))
C     CALL JACO(F,V,NB,NMAX,NROWV,BIG,JBIG)
      DO 610 I = 1,KZYDIM
 610     WRK(KEIVAL-1+I) = WRK(KE2-1 + (I*I+I)/2)
      CALL ORDER(WRK(KEIVEC),WRK(KEIVAL),KZYDIM,0)
C     CALL ORDER(EVEC,EVAL,N,NEVEC)
      WRITE(LUPRI,'(/A)') '   The eigenvalues of E2 :'
      DO 620 I = 1,KZYDIM,2
         WRITE(LUPRI,'(1P,2(I10,D20.10))')
     *      I,WRK(KEIVAL-1+I), I+1, WRK(KEIVAL+I)
 620  CONTINUE
C
C CHECK EIGENVALUES OF S2
C
      WRITE(LUPRI,'(/A,/A)')' ********************************',
     *' CHECK EIGENVALUES OF S2 WITH JACO ROUTINE '
      IJ = 0
      DO 700 J = 1,KZYDIM
         DO 800 I = 1,J
            IJ = IJ + 1
            WRK(KE2-1+IJ) = WRK(KS2MAT-1+(I-1)*KZYDIM+J)
 800     CONTINUE
 700  CONTINUE
      CALL JACO(WRK(KE2),WRK(KEIVEC),KZYDIM,KZYDIM,0,
     *          WRK(KBIG),WRK(KJBIG))
C     CALL JACO(F,V,NB,NMAX,NROWV,BIG,JBIG)
      DO 810 I = 1,KZYDIM
 810     WRK(KEIVAL-1+I) = WRK(KE2-1 + (I*I+I)/2)
      CALL ORDER(WRK(KEIVEC),WRK(KEIVAL),KZYDIM,0)
C     CALL ORDER(EVEC,EVAL,N,NEVEC)
      WRITE(LUPRI,'(/A)') '   The eigenvalues of S2 :'
      DO 820 I = 1,KZDIM
         WRITE(LUPRI,'(1P,2(I10,D20.10))')
     *      I,WRK(KEIVAL-1+I), KZYDIM+1-I, WRK(KEIVAL+KZYDIM-I)
 820  CONTINUE
C
C END OF E2SYM
C
      RETURN
      END
C  /* Deck fckode */
      SUBROUTINE FCKODE(NSIM,ICI,IDI,NCIW,NDIW,ICDTYP,H2,ZYMAT,
     *                  FCONE,FVONE,DENA,DENB)
C
C WRITTEN 13 FEB 1986
C
C PURPOSE:
C    CALCULATE ONE-INDEX TRANSFORMED DENSITY CONTRIBUTION TO
C    ONE-INDEX TRANSFORMED INACTIVE AND ACTIVE FOCK MATRICES
C    FOR A GENERAL TRANSFORMATION MATRIX
C
C    FCONE(P,Q) = SUM(R,J) ( 2*(PQ,RJ)*(ZYMAT(J,R) - ZYMAT(R,J)) +
C                 (PR,JQ)*ZYMAT(R,J) - (PJ,RQ)*ZYMAT(J,R) )
C
C    FVONE(P,Q) = SUM(R,Y) (   (PQ,YR)*(DENB(R,Y) - DENA(R,Y)) +
C                 (PR,YQ)*DENA(R,Y)*0.5 - (PY,RQ)*DENB(R,Y)*0.5 )
C
C THE ROUTINE IS CALLED FOR EACH INTEGRAL DENSITY (CD,**) (C<D)
C
#include <implicit.h>
      DIMENSION H2(NORBT,*),         ZYMAT(NORBT,NORBT,*)
      DIMENSION FCONE(NORBT,NORBT,*),FVONE(NORBT,NORBT,*)
      DIMENSION DENA(NORBT,NASHDI,*),DENB(NORBT,NASHDI,*)
C
C  INFDIM : NASHDI
C
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infpri.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
C -- local constants
C
      PARAMETER ( DP5 = 0.5D0 , D2 = 2.0D0  , DTOL = 1.0D-9 )
C
      ICSYM = ISMO(ICI)
      IDSYM = ISMO(IDI)
      ICDSYM= MULD2H(ICSYM,IDSYM)
C
      DO 200 ISIM = 1,NSIM
         IF (ICDTYP.LE.3) THEN
C
C ADD CONTRIBUTIONS TO FCONE
C  ( INDEX C IS INACTIVE )
C
            IF ((ICDTYP.NE.1) .AND. (.NOT. TRPLET)) THEN
C
C   FCONE(P,Q) = (JR,PQ)*(ZYMAT(J,R)-ZYMAT(R,J))*2
C
               ZYCD = ( ZYMAT(ICI,IDI,ISIM) -
     *                           ZYMAT(IDI,ICI,ISIM))*D2
               IF ( ABS(ZYCD).GT.DTOL ) THEN
                  CALL DAXPY(N2ORBX,ZYCD,H2,1,FCONE(1,1,ISIM),1)
               END IF
            END IF
            IRSYM = MULD2H(KSYMOP,ICSYM)
            NORBR = NORB(IRSYM)
            IROFF = IORB(IRSYM) + 1
            IREND = IORB(IRSYM) + NORBR
            IPSYM = MULD2H(IRSYM,ICDSYM)
            IPOFF = IORB(IPSYM) + 1
            IPEND = IORB(IPSYM) + NORB(IPSYM)
            DO 300 IP =IPOFF,IPEND
C
C   FCONE(P,Q) = SUM(R) (JQ,RP)*ZYMAT(R,J)
C
               FCONE(IP,IDI,ISIM) = FCONE(IP,IDI,ISIM) +
     *            DDOT(NORBR,H2(IROFF,IP),1,ZYMAT(IROFF,ICI,ISIM),1)
C
C   FCONE(Q,P) = SUM(R) (JQ,RP)*ZYMAT(J,R)
C
               DO 400 IR = IROFF,IREND
                  FCONE(IDI,IP,ISIM) = FCONE(IDI,IP,ISIM) -
     *                  H2(IR,IP)*ZYMAT(ICI,IR,ISIM)
 400           CONTINUE
 300        CONTINUE
            IF ( (ICDTYP.EQ.1) .AND. (ICI.NE.IDI) ) THEN
C
C ADD CONTRIBUTIONS FROM DISTRIBUTION (DC,**) WHEN BOTH C AND D ARE
C INACTIVE
               IRSYM = MULD2H(KSYMOP,IDSYM)
               NORBR = NORB(IRSYM)
               IROFF = IORB(IRSYM) + 1
               IREND = IORB(IRSYM) + NORBR
               IPSYM = MULD2H(IRSYM,ICDSYM)
               IPOFF = IORB(IPSYM) + 1
               IPEND = IORB(IPSYM) + NORB(IPSYM)
               DO 500 IP =IPOFF,IPEND
C
C   FCONE(P,Q) = SUM(R) (JQ,RP)*ZYMAT(R,J)
C
                  FCONE(IP,ICI,ISIM) = FCONE(IP,ICI,ISIM) +
     *               DDOT(NORBR,H2(IROFF,IP),1,ZYMAT(IROFF,IDI,ISIM),1)
C
C   FCONE(Q,P) = SUM(R) (JQ,RP)*ZYMAT(J,R)
C
                  DO 600 IR = IROFF,IREND
                     FCONE(ICI,IP,ISIM) = FCONE(ICI,IP,ISIM) -
     *                  H2(IR,IP)*ZYMAT(IDI,IR,ISIM)
 600              CONTINUE
 500           CONTINUE
            END IF
         END IF
C
C ADD CONTRIBUTIONS TO FVONE
C   ( INDEX C IS ACTIVE )
C
         IF (NASHT.GT.0) THEN
            IF (ICDTYP.GE.3) THEN
               IF (ICDTYP.EQ.3) THEN
                  ISWAP = ICI
                  ICI   = IDI
                  IDI   = ISWAP
                  ISWAP = NCIW
                  NCIW  = NDIW
                  NDIW  = ISWAP
                  ISWAP = ICSYM
                  ICSYM = IDSYM
                  IDSYM = ISWAP
               END IF
               IF (.NOT. TRPLET) THEN
C
C   FVONE(P,Q) = (YR,PQ)*(DENB(R,Y)-DENA(R,Y))
C
                  DENBA = (DENB(IDI,NCIW,ISIM) -
     *                              DENA(IDI,NCIW,ISIM))
                  IF (ABS(DENBA).GT.DTOL) THEN
                     CALL DAXPY(N2ORBX,DENBA,H2,1,FVONE(1,1,ISIM),1)
                  END IF
               END IF
               IRSYM = MULD2H(KSYMOP,ICSYM)
               NORBR = NORB(IRSYM)
               IROFF = IORB(IRSYM) + 1
               IREND = IORB(IRSYM) + NORBR
               IPSYM = MULD2H(IRSYM,ICDSYM)
               IPOFF = IORB(IPSYM) + 1
               IPEND = IORB(IPSYM) + NORB(IPSYM)
               DO 700 IP =IPOFF,IPEND
C
C   FVONE(P,Q) = SUM(R) (XQ,RP)*DENA(R,X)*0.5
C
                  FVONE(IP,IDI,ISIM) = FVONE(IP,IDI,ISIM) +
     *            DP5*DDOT(NORBR,H2(IROFF,IP),1,DENA(IROFF,NCIW,ISIM),1)
C
C   FVONE(Q,P) = -SUM(R) (XQ,RP)*DENB(R,X)*0.5
C
                  FVONE(IDI,IP,ISIM) = FVONE(IDI,IP,ISIM) -
     *            DP5*DDOT(NORBR,H2(IROFF,IP),1,DENB(IROFF,NCIW,ISIM),1)
 700           CONTINUE
               IF (ICDTYP.EQ.3) THEN
                  ISWAP = ICI
                  ICI   = IDI
                  IDI   = ISWAP
                  ISWAP = NCIW
                  NCIW  = NDIW
                  NDIW  = ISWAP
                  ISWAP = ICSYM
                  ICSYM = IDSYM
                  IDSYM = ISWAP
               END IF
C
               IF ((ICDTYP.EQ.4) .AND. (ICI.NE.IDI)) THEN
C
C ADD CONTRIBUTIONS FROM (DC,**) DISTRIBUTION WHEN C AND D ARE BOTH
C ACTIVE
C
C
C   FVONE(P,Q) = (YR,PQ)*(DENB(R,Y)-DENA(R,Y))
C
                  IF ( .NOT. TRPLET ) THEN
                     DENBA = (DENB(ICI,NDIW,ISIM) -
     *                              DENA(ICI,NDIW,ISIM))
                     IF (ABS(DENBA).GT.DTOL) THEN
                        CALL DAXPY(N2ORBX,DENBA,H2,1,FVONE(1,1,ISIM),1)
                     END IF
                  END IF
                  IRSYM = MULD2H(KSYMOP,IDSYM)
                  NORBR = NORB(IRSYM)
                  IROFF = IORB(IRSYM) + 1
                  IREND = IORB(IRSYM) + NORBR
                  IPSYM = MULD2H(IRSYM,ICDSYM)
                  IPOFF = IORB(IPSYM) + 1
                  IPEND = IORB(IPSYM) + NORB(IPSYM)
                  DO 800 IP =IPOFF,IPEND
C
C   FVONE(P,Q) = SUM(R) (XQ,RP)*DENA(R,X)*0.5
C
                     FVONE(IP,ICI,ISIM) = FVONE(IP,ICI,ISIM) +
     *               DP5*DDOT(NORBR,H2(IROFF,IP),1,
     *                              DENA(IROFF,NDIW,ISIM),1)
C
C   FVONE(Q,P) = -SUM(R) (XQ,RP)*DENB(R,X)*0.5
C
                     FVONE(ICI,IP,ISIM) = FVONE(ICI,IP,ISIM) -
     *               DP5*DDOT(NORBR,H2(IROFF,IP),1,
     *                              DENB(IROFF,NDIW,ISIM),1)
 800              CONTINUE
               END IF
            END IF
         END IF
 200  CONTINUE
C
C END OF FCKODE
C
      RETURN
      END
C  /* Deck fckoin */
      SUBROUTINE FCKOIN(NSIM,FC,FV,ZYMAT,FCONE,FVONE)
C
C WRITTEN 13 FEB 1986
C revised 941117-hjaaj
C
C PURPOSE:
C    ADD CONTRIBUTION TO ONE-INDEX TRANSFORMED FOCK MATRICES WHICH
C    ORIGINATE FROM TRANSFORMING THE FOCK MATRICES
C
#include <implicit.h>
C
      DIMENSION FC(*),FV(*),ZYMAT(NORBT,NORBT,*)
      DIMENSION FCONE(NORBT,NORBT,*),FVONE(NORBT,NORBT,*)
C
#include <thrzer.h>
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
#include <infrsp.h>
#include <wrkrsp.h>
      IF (IPRRSP.GT.90) THEN
         DO 3000 ISIM = 1,NSIM
            WRITE(LUPRI,*) ISIM,
     &           ' ONE-INDEX TRANSFORMED FOCK CORE MATRIX',
     *           ' BEFORE DIRECT [FC,ZYMAT] CONTRIBUTION'
            CALL OUTPUT(FCONE(1,1,ISIM),1,NORBT,1,NORBT,NORBT,NORBT,
     *                  1,LUPRI)
            IF (NASHT.GT.0) THEN
               WRITE(LUPRI,*) ISIM,
     *           ' ONE-INDEX TRANSFORMED FOCK VALENCE MATRIX',
     *           ' BEFORE DIRECT [FC,ZYMAT] CONTRIBUTION'
               CALL OUTPUT(FVONE(1,1,ISIM),1,NORBT,1,NORBT,NORBT,NORBT
     *                     ,1,LUPRI)
            END IF
 3000    CONTINUE
      END IF
C
C ADD ONE INDEX TRANSFORMED TOTAL SYMMETRIC FOCK MATRIX
C
      DO 1000 ISIM=1,NSIM
         DO 800 IQSYM=1,NSYM
            IORBIQ = IORB(IQSYM)
            NORBIQ = NORB(IQSYM)
            IIORBQ = IIORB(IQSYM)
            IPSYM  = MULD2H(IQSYM,KSYMOP)
            IORBIP = IORB(IPSYM)
            NORBIP = NORB(IPSYM)
            IIORBP = IIORB(IPSYM)
            DO 700 IQ=1,NORBIQ
C
C   FCONE(P,Q) = FCONE(P,Q)   + SUM(R)  ZYMAT(P,R) * FC(R,Q)
C   FCONE(Q,P) = FCONE(Q,P)   - SUM(R)  FC(Q,R) * ZYMAT(R,P)
C
               IQROFF = IIORBQ + IROW(IQ)
               IRQOFF = IIORBQ + IQ
               DO 600 IR=1,NORBIQ
                  IF (IR .GT. IQ) THEN
                     IRQ = IRQOFF + IROW(IR)
                  ELSE
                     IRQ = IQROFF + IR
                  END IF
                  IF ( ABS(FC(IRQ)) .GT. THRZER) THEN
                  DO 510 IP=1,NORBIP
                     FCONE(IORBIP+IP,IORBIQ+IQ,ISIM) =
     *                  FCONE(IORBIP+IP,IORBIQ+IQ,ISIM) +
     *                  ZYMAT(IORBIP+IP,IORBIQ+IR,ISIM) * FC(IRQ)
                     FCONE(IORBIQ+IQ,IORBIP+IP,ISIM) =
     *                  FCONE(IORBIQ+IQ,IORBIP+IP,ISIM) -
     *                  ZYMAT(IORBIQ+IR,IORBIP+IP,ISIM) * FC(IRQ)
 510              CONTINUE
                  END IF
                  IF ( NASHT.GT.0 ) THEN
                     IF (ABS(FV(IRQ)) .GT. THRZER) THEN
                        DO 520 IP=1,NORBIP
                           FVONE(IORBIP+IP,IORBIQ+IQ,ISIM) =
     *                          FVONE(IORBIP+IP,IORBIQ+IQ,ISIM) +
     *                          ZYMAT(IORBIP+IP,IORBIQ+IR,ISIM) *FV(IRQ)
                           FVONE(IORBIQ+IQ,IORBIP+IP,ISIM) =
     &                          FVONE(IORBIQ+IQ,IORBIP+IP,ISIM) -
     &                          ZYMAT(IORBIQ+IR,IORBIP+IP,ISIM) *FV(IRQ)
 520                    CONTINUE
                     END IF
                  END IF
 600           CONTINUE
C
C
 700        CONTINUE
 800     CONTINUE
 1000 CONTINUE
C
      IF (IPRRSP.GT.90) THEN
         DO 2000 ISIM = 1,NSIM
            WRITE(LUPRI,*) ISIM,' ONE-INDEX TRANSFORMED FOCK CORE '//
     &                          'MATRIX'
            CALL OUTPUT(FCONE(1,1,ISIM),1,NORBT,1,NORBT,NORBT,NORBT
     *                  ,1,LUPRI)
            IF (NASHT.GT.0) THEN
               WRITE(LUPRI,*) ISIM,
     *                   ' ONE-INDEX TRANSFORMED FOCK VALENCE MATRIX'
               CALL OUTPUT(FVONE(1,1,ISIM),1,NORBT,1,NORBT,NORBT,NORBT
     *                     ,1,LUPRI)
            END IF
 2000    CONTINUE
         WRITE(LUPRI,*)' CORE FOCK MATRIX '
         CALL OUTPKB(FC(1),NORB,NSYM,1,LUPRI)
         IF (NASHT.GT.0) THEN
            WRITE (LUPRI,*)' VALENCE FOCK MATRIX '
            CALL OUTPKB(FV(1),NORB,NSYM,1,LUPRI)
         END IF
      END IF
C
C END OF FCKOIN
C
      RETURN
      END
C  /* Deck fvtsd */
      SUBROUTINE FVTSD(NSIM,ICI,IDI,NCIW,NDIW,ICDTYP,FVTD,H2,DVT)
C
C WRITTEN 13 FEB 1986
C
C PURPOSE:
C    CALCULATE ACTIVE FOCK MATRICES WITH TRANSITION DENSITY MATRIX
C
C    FVTD(P,Q) = SUM(X,Y) [(PQ,XY)-0.5*(PY,XQ)] * DVT(X,Y)
C
C THE ROUTINE IS CALLED FOR EACH INTEGRAL DENSITY (CD,AB) (C<D)
C
#include <implicit.h>
      DIMENSION FVTD(NORBT,NORBT,*),H2(NORBT,*),DVT(NASHDI,NASHDI,*)
C
C  INFDIM : NASHDI
C
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infpri.h>
#include <wrkrsp.h>
#include <infrsp.h>
C
C -- local constants
C
      PARAMETER ( DP5 = 0.5D0 , D0 = 0.0D0 , D2 = 2.0D0  )
      PARAMETER ( DTOL = 1.0D-9 )
C
      IF (ICDTYP.EQ.3) THEN
         ISWAP = ICI
         ICI   = IDI
         IDI   = ISWAP
         ISWAP = NCIW
         NCIW  = NDIW
         NDIW  = ISWAP
      END IF
      ICSYM = ISMO(ICI)
      IDSYM = ISMO(IDI)
      ICDSYM= MULD2H(ICSYM,IDSYM)
      IF ((ICDTYP.EQ.4) .AND. (.NOT. TRPLET)) THEN
C
C   FVONE(P,Q) = (XY,PQ)*DVT(X,Y)
C
         DO 100 ISIM = 1,NSIM
            IF (ABS(DVT(NCIW,NDIW,ISIM)).GT.DTOL) THEN
               CALL DAXPY(N2ORBX,DVT(NCIW,NDIW,ISIM),H2,1,
     *                                      FVTD(1,1,ISIM),1)
            END IF
            IF (ICI.NE.IDI) THEN
               IF (ABS(DVT(NDIW,NCIW,ISIM)).GT.DTOL) THEN
                  CALL DAXPY(N2ORBX,DVT(NDIW,NCIW,ISIM),H2,1,
     *                                      FVTD(1,1,ISIM),1)
               END IF
            END IF
 100     CONTINUE
      END IF
C
C   FVTD(P,Q) = SUM(X) (YP,XQ)*DVT(X,Y)*DP5
C
      DO 150 ISIM = 1,NSIM
         IXSYM = MULD2H(KSYMOP,ICSYM)
         IXOFF = IORB(IXSYM)+1+NISH(IXSYM)
         IXWOFF= ISW(IXOFF) - NISHT
         NASHX = NASH(IXSYM)
         IQSYM = MULD2H(IXSYM,ICDSYM)
         IQOFF = IORB(IQSYM) + 1
         IQEND = IORB(IQSYM) + NORB(IQSYM)
         DO 200 IQ =IQOFF,IQEND
            FVTD(IDI,IQ,ISIM) = FVTD(IDI,IQ,ISIM) -
     *         DP5* DDOT(NASHX,H2(IXOFF,IQ),1,DVT(IXWOFF,NCIW,ISIM),1)
 200     CONTINUE
         IF ((ICDTYP.EQ.4).AND.(ICI.NE.IDI)) THEN
C
C ADD CONTRIBUTION FOR (DC,**) DISTRIBUTION WHEN D AND C BOTH ARE
C ACTIVE ORBITALS
C
            IXSYM = MULD2H(KSYMOP,IDSYM)
            IXOFF = IORB(IXSYM)+1+NISH(IXSYM)
            IXWOFF= ISW(IXOFF) - NISHT
            NASHX = NASH(IXSYM)
            IQSYM = MULD2H(IXSYM,ICDSYM)
            IQOFF = IORB(IQSYM) + 1
            IQEND = IORB(IQSYM) + NORB(IQSYM)
C
C   FVTD(P,Q) = SUM(X) (YP,XQ)*DVT(X,Y)*DP5
C
            DO 300 IQ =IQOFF,IQEND
               FVTD(ICI,IQ,ISIM) = FVTD(ICI,IQ,ISIM) -
     *         DP5* DDOT(NASHX,H2(IXOFF,IQ),1,DVT(IXWOFF,NDIW,ISIM),1)
 300        CONTINUE
         END IF
 150  CONTINUE
      IF (ICDTYP.EQ.3) THEN
         ISWAP = ICI
         ICI   = IDI
         IDI   = ISWAP
         ISWAP = NCIW
         NCIW  = NDIW
         NDIW  = ISWAP
      END IF
C
C END OF FVTSD
C
      RETURN
      END
C  /* Deck h2xsig */
      SUBROUTINE H2XSIG(NOSIM,FCONE,ZYMAT,H2XAC,
     *                  EVECS,XINDX,WRK,LWRK)
C
C CALCULATE FOR AN ORBITAL TRIAL VECTOR CONFIGURATION PART OF
C LINEAR TRANSFORMATION
C
C                ( <J,H2X,0> )
C     E[2]*N = - (           )
C                (-<0,H2X,J> )
C
C H2X ARE ONE INDEX TRANSFORMED INTEGRALS
C
C H2XAC CONTAIN HALF TRANSFORMED INTEGRALS. THE ONE INDEX TRANSFORMED
C INTEGRALS ARE OBTAINED FROM H2XAC IN GETINT.
C
#include <implicit.h>
C
      DIMENSION FCONE(NORBT,NORBT,*),ZYMAT(NORBT,NORBT,*)
      DIMENSION H2XAC(N2ASHX,NNASHX,*)
      DIMENSION EVECS(KZYVAR,*),XINDX(*),WRK(*)
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infopt.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
#include <infdim.h>
#include <cbgetdis.h>
C
      LOGICAL NOH2
C
      PARAMETER ( DM1 = -1.0D0 )
C
C ALLOCATE WORK SPACE
C
      IF (IREFSY .EQ. KSYMST) THEN
         NDREF = KZCONF
C        ... if KZCONF .ne. NCREF, we need CREF in determinants
      ELSE
         NDREF = NCREF
      END IF
      KUFCAC = 1
      KREFCO = KUFCAC + N2ASHX
      KFREE  = KREFCO + NDREF
      LFREE  = LWRK   - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('H2XSIG',KFREE-1,LWRK)
      IADH2X = -1
C     ... IADH2X .lt. 0: H2XAC in core
C
      IF ( TRPLET ) THEN
         ISPIN1 = 1
         ISPIN2 = 0
      ELSE
         ISPIN1 = 0
         ISPIN2 = 0
      END IF
C
      CALL GETREF(WRK(KREFCO),NDREF)
      DO 100 ISIM = 1,NOSIM
C
C CREATE Z PART  -<J,H2X,0> OF LINEAR TRANSFORMATION WITH E(2)
C
C        Get transposed FCXAC matrix for Z sigma vector
C        (note: CISIGD requires UFCAC(I,J) = FCXAC(J,I))
C
         DO 220 IW = 1,NASHT
            IX = ISX(NISHT+IW)
            DO 230 JW = 1,NASHT
               JX = ISX(NISHT+JW)
               IJW = (IW-1) * NASHT + JW
               WRK(KUFCAC-1+IJW) = FCONE(IX,JX,ISIM)
 230        CONTINUE
 220     CONTINUE
C
         DISTYP = 2
C
         CALL CISIGD(IREFSY,KSYMST,NDREF,KZCONF,WRK(KREFCO),
     *               EVECS(1,ISIM),
     *               WRK(KUFCAC),H2XAC(1,1,ISIM),.FALSE.,.FALSE.,
     *               XINDX,ISPIN1,ISPIN2,WRK(KFREE),LFREE)
C        CALL RSPSIG(ICSYM,IHCSYM,NCDET,NHCDET,C,HC,UFCAC,H2AC,IFLAG,
C    *               NOH2,WORK,KFREE,LFREE)
C
C        IF ((NDREF.NE.NCDET).OR.(KZCONF.NE.NHCDET)) THEN
C           WRITE(LUPRI,'(/2(A,I5),/3(A,I5))')
C    *      ' NUMBER OF REFERENCE DETERMINANTS ,NDREF:',NDREF,
C    *      ' CALCULATED NUMBER ,NCDET:',NCDET,
C    *      ' NUMBER OF DETERMINANTS FOR SYMMETRY',KSYMOP,
C    *      '  IS:',KZCONF,'  CALCULATED NUMBER,NHCDET:',NHCDET
C           CALL QUIT(' H2XSIG,INCORRECT CALCULATION OF DETERMINANTS')
C        END IF
         CALL DSCAL(KZCONF,DM1,EVECS(1,ISIM),1)
         IF (IREFSY .EQ. KSYMST) THEN
C           ... remove CREF component of E[2] vector
            T1 = DDOT(KZCONF,WRK(KREFCO),1,EVECS(1,ISIM),1)
            CALL DAXPY(KZCONF,(-T1),WRK(KREFCO),1,EVECS(1,ISIM),1)
         END IF
         IF (IPRRSP.GT.110) THEN
            WRITE(LUPRI,'(/A)')
     *        ' CSF PART OF E(2) TRANSFORMED ORBITAL VECTOR: Z PART'
            CALL OUTPUT(EVECS(1,ISIM),1,KZCONF,1,1,
     *                  KZCONF,1,1,LUPRI)
         END IF
C
C CREATE Y PART  <0,H2X,J> OF LINEAR TRANSFORMATION WITH E(2)
C
C        Get FCXAC matrix for Y sigma vector
C
         DO 120 IW = 1,NASHT
            IX = ISX(NISHT+IW)
            DO 130 JW = 1,NASHT
               JX = ISX(NISHT+JW)
               IJW = (IW-1) * NASHT + JW
               WRK(KUFCAC-1+IJW) = FCONE(JX,IX,ISIM)
 130        CONTINUE
 120     CONTINUE
         IF (IPRRSP.GT.110) THEN
            WRITE(LUPRI,'(/A)')' ONE INDEX TRANSFORMED FCAC '
            CALL OUTPUT(WRK(KUFCAC),1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         END IF
C
         DISTYP = 3
         CALL CISIGD(IREFSY,KSYMST,NDREF,KZCONF,WRK(KREFCO),
     *               EVECS(1+KZVAR,ISIM),
     *               WRK(KUFCAC),H2XAC(1,1,ISIM),.FALSE.,.FALSE.,
     *               XINDX,ISPIN1,ISPIN2,WRK(KFREE),LFREE)
         IF (IREFSY .EQ. KSYMST) THEN
C           ... remove CREF component of E[2] vector
            T1 = DDOT(KZCONF,WRK(KREFCO),1,EVECS(1+KZVAR,ISIM),1)
            CALL DAXPY(KZCONF,(-T1),WRK(KREFCO),1,EVECS(1+KZVAR,ISIM),1)
         END IF
         IF (IPRRSP.GT.110) THEN
            WRITE(LUPRI,'(/A)')
     *        ' CSF PART OF E(2) TRANSFORMED ORBITAL VECTOR: Y PART'
            CALL OUTPUT(EVECS(1+KZVAR,ISIM),1,KZCONF,1,1,
     *                  KZCONF,1,1,LUPRI)
         END IF
C
 100  CONTINUE
      RETURN
C
C     ... END OF H2XSIG.
C
      END
C  /* Deck lrst */
      SUBROUTINE LRST(BINMEM,MAXSIM,IBTYP,EIVAL,EIVEC,GD,
     *                REDGD,REDS,REDE,UDV,FCAC,H2AC,XINDX,
     *                LAB1,LAB2,WRK,LWRK)
C
C PURPOSE: CREATE START VECTOR(S) FOR SOLUTION OF A LINEAR
C          SET OF EQUATIONS
C
C
C          USE GRADIENT VECTOR MULTIPLIED WITH INVERSE
C          DIAGONAL HESSIAN MATRIX ELEMENTS
C
C          IF NTYPE.EQ.1 USE OLD TRIAL VECTORS
C
#include <implicit.h>
#include <dummy.h>
      CHARACTER*8 LAB1, LAB2
      DIMENSION IBTYP(*),EIVAL(*),EIVEC(*),GD(*),REDGD(*)
      DIMENSION REDS(*),REDE(*),UDV(*),FCAC(*),H2AC(*),XINDX(*),WRK(*)
      PARAMETER ( DTEST = 1.0D-4, D0=0.0D0, D1=1.0D0, D1TEST = 1.0D-6)
      PARAMETER ( DM1 = -1.0D0 )
#include <thrldp.h>
#include <ibndxdef.h>
C
C Used from common blocks:
C   INFRSP : ???,SOPPA
C
#include <priunit.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <inftap.h>
#include <infpri.h>
#include <infdim.h>
#include <infopt.h>
#include <phpinf.h>
C
      LOGICAL BINMEM
      LOGICAL FNDLAB, FOUND, CONV
C
      KZRED = 0
      IF ( RESTLR ) THEN
C
         IF (RSPSUP .AND. (KSYMOP.EQ.1)) THEN
            WRITE (LUERR,*)
     *      'LRST: RSPAVE necessary for RESTLR but not yet implemented'
            CALL QUIT('LRST error: Implement RSPAVE for RESTLR')
         END IF
C
         LURSP1 = -1
         CALL GPOPEN(LURSP1,'RSPRST.E2C','UNKNOWN',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND (LURSP1)
         IF ( .NOT.FNDLAB('RESTART ',LURSP1) ) THEN
            WRITE (LUPRI,'(/2A,I4)')
     *      ' LRST: Restart not possible, RESTART label not found on',
     *      ' LURSP1 =',LURSP1
            CALL QUIT('LRST: Restart impossible, '//
     &                'RESTART label not found.')
         END IF
         READ (LURSP1)JSYMOP, KZDIM,(IBTYP(I),I=1,KZDIM)
         IF ( IPRRSP.GT. 4 ) THEN
            WRITE (LUPRI,'(/A,2(/A,I5),/A,(T10,10I4))')
     *      ' >>> Restart in LRST >>>',
     *      ' OPERATOR SYMMETRY          :',JSYMOP,
     *      ' DIMENSION OF REDUCED SPACE :',KZDIM,
     *      ' IBTYP :',(IBTYP(I),I=1,KZDIM)
         END IF
         IF (JSYMOP.NE.KSYMOP) THEN
            WRITE (LUPRI,'(/,(A,I4))')
     *      ' LRST: Restart not possible, symmetry on restart file',
     *      JSYMOP,' symmetry of current operator',KSYMOP
            CALL QUIT('LRST: Restart impossible, no match of symmetry ')
         END IF
         CALL GPCLOSE(LURSP1,'KEEP')
C
C Check how many sigma vectors have been written to LURSP5
C
         CALL REWSPL(LURSP5)
         JRSP5 = 0
         JZDIM = 0
         DO 50 I = 1,KZDIM
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP5,LRSP5,DUMMY,JRSP5,.FALSE.)
            IF (LRSP5 .EQ. -1) GOTO 51
#else
            READ (LURSP5,END=51,ERR=51) JUNK
#endif
            JZDIM = I
  50     CONTINUE
  51     KZDIM = JZDIM
C
C
C READ IN OLD TRIAL VECTORS AND SET UP REDUCED MATRICES
C
C REPEAT UNTILL
C
         NLOAD  = 0
         JLRSTV = KZDIM
         CALL REWSPL(LURSP3)
         CALL REWSPL(LURSP5)
#if defined (VAR_SPLITFILES)
         JRSP3 = 0
         JRSP5 = 0
#endif
 200     CONTINUE
            NLOAD  = NLOAD +1
            IF ((NLOAD.EQ.1).AND.(KOFFTY.EQ.1)) THEN
#if defined (VAR_SPLITFILES)
               CALL READST(LURSP3,LRSP3,DUMMY,JRSP3,.FALSE.)
               CALL READST(LURSP5,LRSP5,DUMMY,JRSP5,.FALSE.)
#else
               READ (LURSP3)
               READ (LURSP5)
#endif
               JLRSTV = JLRSTV - 1
            ENDIF
            IF ( NLOAD .GT. 1) THEN
               CALL REWSPL(LURSP3)
#if !defined (VAR_SPLITFILES)
               IF ( KOFFTY.EQ.1)READ (LURSP3)
#else
               JRSP3 = 0
               IF (KOFFTY.EQ.1)
     &              CALL READST(LURSP3,LRSP3,DUMMY,JRSP3,.FALSE.)
#endif
               DO 61 IVEC = 1,KZRED
#if defined (VAR_SPLITFILES)
                  CALL READST(LURSP3,LRSP3,DUMMY,JRSP3,.FALSE.)
#else
                  READ (LURSP3)
#endif
 61            CONTINUE
            END IF
            NSIM   = MIN(MAXSIM,JLRSTV)
            NCSIM  = 0
            NOSIM  = 0
            DO 410 ISIM=1,NSIM
               IF (IBTYP(KOFFTY+KZRED+ISIM).EQ.JBCNDX) THEN
                  NCSIM=NCSIM+1
               ELSE
                  NOSIM=NOSIM+1
               ENDIF
 410        CONTINUE
            KBCVEC = 1
            KBOVEC = KBCVEC + NCSIM*KZCONF
            KECVEC = KBOVEC + NOSIM*KZYWOP
            KEOVEC = KECVEC + NCSIM*KZYVAR
            KDIAE  = KEOVEC + NOSIM*KZYVAR
            IF (SOPPA) THEN
               KWRK = KDIAE + KZCONF
            ELSE
               KWRK = KDIAE
            ENDIF
            LWRKRE = LWRK   - KWRK
            IF (LWRKRE.LT.0) CALL ERRWRK('LRST 1',KWRK-1,LWRK)
            ISTBO  = KBOVEC
            ISTBC  = KBCVEC
            ISTEC  = KECVEC
            ISTEO  = KEOVEC
            IF (SOPPA) THEN
               CALL REWSPL(LURSP4)
#if defined (VAR_SPLITFILES)
               JRSP4 = 0
               CALL READST(LURSP4,LRSP4,WRK(KDIAE),JRSP4,.TRUE.)
               IF (LRSP4 .NE. KZCONF) CALL QUIT('Inconsistent record '//
     &              'length on unit LURSP4 in LRST')
#else
               CALL READT(LURSP4,KZCONF,WRK(KDIAE))
#endif
            ENDIF
            DO 400 ISIM=1,NSIM
               IF (IBTYP(KOFFTY+KZRED+ISIM).EQ.JBCNDX) THEN
#if defined (VAR_SPLITFILES)
                  CALL READST(LURSP3,LRSP3,WRK(ISTBC),JRSP3,.TRUE.)
                  IF (LRSP3 .NE. KZCONF) CALL QUIT('Inconsistent '//
     &                 'record length on unit LURSP3 in LRST')
#else
                  CALL READT(LURSP3,KZCONF,WRK(ISTBC))
#endif
C
C   If SOPPA read in the p-h part of the transformed vector
C   and construct the 2p-2h part and put everything as if read
C
                  IF (SOPPA) THEN
#if defined (VAR_SPLITFILES)
                     CALL READST(LURSP5,LRSP5,
     &                           WRK(ISTEC+KZCONF),JRSP5,.TRUE.)
                     IF (LRSP5 .NE. KZYWOP) CALL QUIT('Inconsistent '//
     &                    'record length on unit LURSP5 in LRST')
#else
                     CALL READT(LURSP5,KZYWOP,WRK(ISTEC+KZCONF))
#endif
                     CALL DCOPY(KZWOPT,WRK(ISTEC+KZVAR),1,
     *                          WRK(ISTEC+KZVAR+KZCONF),1)
                     DO I=0,KZCONF-1
                        WRK(ISTEC+I) = WRK(KDIAE+I) * WRK(ISTBC+I)
                     ENDDO
                     CALL DZERO(WRK(ISTEC+KZVAR),KZCONF)
                  ELSE
#if defined (VAR_SPLITFILES)
                     CALL READST(LURSP5,LRSP5,WRK(ISTEC),JRSP5,.TRUE.)
                     IF (LRSP5 .NE. KZYVAR) CALL QUIT('Inconsistent '//
     &                    'record length on unit LURSP5 in LRST (2)')
#else
                     CALL READT(LURSP5,KZYVAR,WRK(ISTEC))
#endif
                  ENDIF
                  ISTBC = ISTBC + KZCONF
                  ISTEC = ISTEC + KZYVAR
               ELSE
#if defined (VAR_SPLITFILES)
                  CALL READST(LURSP3,LRSP3,WRK(ISTBO),JRSP3,.TRUE.)
                  IF (LRSP3 .NE. KZYWOP) CALL QUIT('Inconsistent '//
     &                 'record length on unit LURSP3 in LRST')
                  CALL READST(LURSP5,LRSP5,WRK(ISTEO),JRSP5,.TRUE.)
                     IF (LRSP5 .NE. KZYVAR) CALL QUIT('Inconsistent '//
     &                    'record length on unit LURSP5 in LRST (3)')
#else
                  CALL READT(LURSP3,KZYWOP,WRK(ISTBO))
                  CALL READT(LURSP5,KZYVAR,WRK(ISTEO))
#endif
                  ISTBO = ISTBO + KZYWOP
                  ISTEO = ISTEO + KZYVAR
               ENDIF
 400        CONTINUE
            KZRED=KZRED+NSIM
            KZYRED=2*KZRED
            JLRSTV=JLRSTV-NSIM
C
C           CALL RSPRED(1,..) INCREASE DIMENSION OF REDUCED RSP EQUATION
C
            CALL RSPRED (1,.TRUE.,NSIM,IBTYP,GD,REDGD,REDE,REDS,
     *                 EIVAL,EIVEC,WRK(KBCVEC),WRK(KBOVEC),
     *                 UDV,WRK(KECVEC),XINDX,WRK(KWRK),LWRKRE)
C
C           CALL RSPRED (ICTL,LINEQ,N,IBTYP,GD,REDGD,REDE,REDS,
C    *                 EIVAL,EIVEC,BCVEC,BOVEC,SVECS,EVECS,XINDX,WRK)
C
            IF (JLRSTV.LE.0) THEN
                IF (IPRRSP.GE.3)WRITE (LUPRI,'(//A,I5,A/)')
     *             ' ***',KZRED,' OLD TRIAL VECTORS FOUND ON LURSP3'
            ELSE
               GO TO 200
C        ^--------------
            END IF
      END IF
      IF ((.NOT. RESTLR). OR. (KZRED .EQ. 0)) THEN
         IF (KOFFTY.EQ.1) THEN
            CALL REWSPL(LURSP3)
            CALL REWSPL(LURSP5)
#if defined (VAR_SPLITFILES)
            JRSP3 = 0
            JRSP5 = 0
#endif
            CALL GETREF(WRK,KZCONF)
#if defined (VAR_SPLITFILES)
            CALL WRITST(LURSP3,KZCONF,WRK,JRSP3)
            CALL WRITST(LURSP5,1,DUMMY,JRSP5)
#else
            CALL WRITT(LURSP3,KZCONF,WRK)
            WRITE (LURSP5) D0, D0, D0, D0
#endif
Chj         ... rec. 1 on LURSP5 should never be read when KOFFTY .eq. 1
            IBTYP(1) = JBCNDX
            IF (IPRRSP.GT.75) THEN
               WRITE(LUPRI,'(/A)')
     *         ' ** LRST ** REFERENCE VECTOR WRITTEN ON LURSP3 rec. # 1'
            END IF
         ENDIF
      ENDIF
C
C CREATE NEW ORBITAL AND CSF TRIAL VECTORS BY DIVIDING GRADIENTS
C WITH DIAGONAL HESSIAN ELEMENTS
C
C WORK SPACE ALLOCATION
C
C     CALL PHPINI(LPHPMX,KZCONF,KZWOPT,MAXPHP)
      KDIAE  = 1
      KDIASO = KDIAE  + LPHPMX
      KTOT   = KDIASO + KZWOPT
      LTOT   = LWRK   - KTOT
      IF (LTOT.LT.0) CALL ERRWRK('LRST 2',KTOT-1,LWRK)
      IF (KZCONF.GT.0) THEN
         CALL PHPDSK(WRK(KDIAE),LPHPMX,.FALSE.)
      ELSE
         CALL RSPEDG(WRK(KDIAE))
      END IF
      IF (KZWOPT.GT.0) CALL RSPSOD(WRK(KDIASO))
      NLOAD = 0
      NTSIM = 0
      IF (IPRRSP.GT.110) THEN
         WRITE(LUPRI,*)'  GRADIENT VECTOR (Z-part and Y-part)'
         CALL OUTPUT(GD,1,KZVAR,1,2,KZVAR,2,1,LUPRI)
      ENDIF
      IF (IPRRSP.GT.150) THEN
         WRITE(LUPRI,*)' DIAGONAL HESSIAN MATRIX WVAL=0.0'
         CALL OUTPUT(WRK(KDIAE),1,KZVAR,1,1,KZVAR,1,1,LUPRI)
      ENDIF
      DO 150 ISIM = 1,KEXSIM,MAXSIM
         NBX = MIN(MAXSIM,(KEXSIM+1-ISIM))
         KBCVEC = KTOT
         IF (SOPPA) THEN
            KBOVEC = KBCVEC
         ELSE
            KBOVEC = KBCVEC + 2*NBX*KZCONF
         END IF
         KWRK1  = KBOVEC + NBX*KZYWOP
         LWRK1  = LWRK   - KWRK1 - KZYVAR
C
C     We might read in a response vector from file, and do this in
C     WRK(KWRK1).
C
         IF (LWRK1.LT.0) CALL ERRWRK('LRST 3',KWRK1-1,LWRK)
         ISTBC  = KBCVEC
         ISTBO  = KBOVEC
         NLOAD = NLOAD + 1
         NCSIM = 0
         NOSIM = 0
         DO 160 IR = 1,NBX
            IROFF = ISIM-1+IR
            WVAL  = EIVAL(IROFF)
            CALL REARSP(LURSP,KLEN,WRK(KWRK1),LAB1,LAB2,WVAL,D0,
     &                  KSYMOP,0,THCRSP,FOUND,CONV,ANTSYM)
            IF (FOUND .AND. .NOT. CONV) THEN
               IF (KZCONF .GT. 0) THEN
                  CALL DCOPY(KZCONF,WRK(KWRK1),1,WRK(ISTBC),1)
                  NCSIM = NCSIM + 1
                  ISTBC = ISTBC + KZCONF
                  IF (WVAL .NE. D0) THEN
                     CALL DCOPY(KZCONF,WRK(KWRK1+KZVAR),1,WRK(ISTBC),1)
                     NCSIM = NCSIM + 1
                     ISTBC = ISTBC + KZCONF
                  END IF
               END IF
               IF (KZWOPT .GT. 0) THEN
                  CALL DCOPY(KZWOPT,WRK(KWRK1+KZCONF),1,WRK(ISTBO),1)
                  CALL DCOPY(KZWOPT,WRK(KWRK1+KZVAR+KZCONF),1,
     &                 WRK(ISTBO+KZWOPT),1)
                  NOSIM = NOSIM + 1
                  ISTBO = ISTBO + KZYWOP
               END IF
               IF (WVAL .LT. 0) CALL DSWAP(KZWOPT,WRK(KWRK1),1,
     &                                     WRK(KWRK1 + KZWOPT),1)
               GO TO 160
            END IF
            IF (KZCONF.GT.0 .AND. .NOT.SOPPA) THEN
               TXNRM = DNORM2(KZCONF,GD,1)
               TYNRM = DNORM2(KZCONF,GD(KZVAR+1),1)
               IF (TXNRM .GT. TYNRM) THEN
                  KGDC = 1
               ELSE
                  KGDC = KZVAR + 1
               END IF
               CALL NEXCI(.FALSE.,WVAL,KZCONF,WRK(ISTBC),DUMMY,
     *                    GD(KGDC),WRK(KDIAE),IPRRSP,WRK(KWRK1),LWRK1)
               NCSIM = NCSIM + 1
               ISTBC = ISTBC + KZCONF
            ENDIF
            IF (KZWOPT.GE.1) THEN
               DO 105 I=1,KZWOPT
                  IZ = I + KZCONF
                  ZDIA = WRK(KDIAE-1+IZ) - WVAL * WRK(KDIASO-1+I)
                  YDIA = WRK(KDIAE-1+IZ) + WVAL * WRK(KDIASO-1+I)
                  IF (ABS(ZDIA).LT.DTEST ) THEN
                     WRK(ISTBO-1+I) = GD(KZCONF+I)
     *                           /SIGN(DTEST,ZDIA)
                  ELSE
                     WRK(ISTBO-1+I) = GD(KZCONF+I) / ZDIA
                  ENDIF
                  IF (ABS(YDIA).LT.DTEST ) THEN
                     WRK(ISTBO-1+KZWOPT+I) = GD(KZVAR+KZCONF+I)
     *                          /SIGN(DTEST,YDIA)
                  ELSE
                     WRK(ISTBO-1+KZWOPT+I) = GD(I+KZVAR+KZCONF) / YDIA
                  ENDIF
 105           CONTINUE
               NOSIM  = NOSIM  + 1
               ISTBO = ISTBO + KZYWOP
            ENDIF
 160     CONTINUE
         NLSIM = NOSIM + NCSIM
         DO 170 IR = 1 , NLSIM
            IF (IR.LE.NCSIM) IBTYP(KOFFTY+KZRED+NTSIM+IR)=JBCNDX
            IF (IR.GT.NCSIM) IBTYP(KOFFTY+KZRED+NTSIM+IR)=JBONDX
 170     CONTINUE
         IF ((IPRRSP.GT.110).AND.(NCSIM.GT.0)) THEN
            WRITE(LUPRI,'(/A,I5)')' NCSIM CSF TRIAL VECTORS',NCSIM
            CALL OUTPUT(WRK(KBCVEC),1,KZCONF,1,NCSIM,
     &                  KZCONF,NCSIM,1,LUPRI)
         ENDIF
         IF ((IPRRSP.GT.100).AND.(NOSIM.GT.0)) THEN
            WRITE(LUPRI,'(/A,I5)')' NOSIM (Z,Y) ORBITAL TRIAL VECTORS',
     &           NOSIM
            CALL OUTPUT(WRK(KBOVEC),1,KZWOPT,1,2*NOSIM,
     *                  KZWOPT,2*NOSIM,1,LUPRI)
         ENDIF
         IF (IPRRSP.GT.80) THEN
            WRITE(LUPRI,'(/A,I5)')
     *         '    I       IBTYP(I)          I=1,NLSIM NLSIM',NLSIM
            DO 3007 I=1+KOFFTY,NLSIM+KOFFTY
               WRITE(LUPRI,'(I10,I10)')I,IBTYP(I)
 3007       CONTINUE
         END IF
         THRLDV = KZYVAR*THRLDP
C
C        hj-aug2000: we must make sure that configuration and orbital
C        trial vectors are consecutive in memory before CALL RSPORT
C
         ISTBO = KBCVEC+NCSIM*KZCONF
         IF (NOSIM .GT. 0 .AND. KBOVEC .NE. ISTBO) THEN
            DO I = 1, NOSIM*KZYWOP
               WRK(ISTBO-1+I) = WRK(KBOVEC-1+I)
            END DO
         END IF
         NBPREV_here = KOFFTY + KZRED + NTSIM
         CALL RSPORT(WRK(KBCVEC),NLSIM,NBPREV_here,IBTYP,
     *               THRLDV,WRK(KWRK1),LURSP3)
C        CALL RSPORT (BVECS,NBX,NBPREV,IBTYP,THRLDP,OLDVEC,LU3)
         IF (IPRRSP.GT.80) THEN
            WRITE(LUPRI,'(/A//A,I5)')' AFTER RSPORT IN LRST',
     *         '    I       IBTYP(I)          I=1,NLSIM NLSIM=',NLSIM
            DO 3008 I=1+KOFFTY,NLSIM+KOFFTY
               WRITE(LUPRI,'(I5,I10)')I,IBTYP(I)
 3008       CONTINUE
         END IF
         NCSIM = 0
         NOSIM = 0
         DO 3010 IX = 1,NLSIM
            IF (IBTYP(KZRED+KOFFTY+NTSIM+IX).EQ.JBCNDX) THEN
               NCSIM = NCSIM + 1
            ELSE
               NOSIM = NOSIM + 1
            ENDIF
 3010    CONTINUE
         IF ((IPRRSP.GT.110).AND.(NCSIM.GT.0)) THEN
            WRITE(LUPRI,*)' NCSIM CSF TRIAL VECTORS',NCSIM
            CALL OUTPUT(WRK(KBCVEC),1,KZCONF,1,NCSIM,KZCONF,NCSIM,1,
     &                  LUPRI)
         ENDIF
         IF ((IPRRSP.GT.110).AND.(NOSIM.GT.0)) THEN
            WRITE(LUPRI,*)' NOSIM ORBITAL TRIAL VECTORS',NOSIM
            CALL OUTPUT(WRK(KBCVEC+NCSIM*KZCONF),1,KZYWOP,
     *                1,NOSIM,KZYWOP,NOSIM,1,LUPRI)
         ENDIF
         NTSIM = NTSIM + NLSIM
 150  CONTINUE
      KEXSTV = NTSIM
      IF ((NLOAD.EQ.1).AND.(KEXSTV.LE.MAXSIM)) THEN
         BINMEM = .TRUE.
         NTOT = NCSIM*KZCONF + NOSIM*KZYWOP
         IF (KBCVEC .NE. 1) THEN
            DO 156 II = 1,NTOT
               WRK(II) = WRK(KBCVEC-1+II)
 156        CONTINUE
C           CALL DCOPY(NTOT,WRK(KBCVEC),1,WRK(1),1)
         END IF
      ELSE
         BINMEM = .FALSE.
      ENDIF
C
      IF (KEXSTV .LE. 0) THEN
         IF (IPRRSP.GE.3) WRITE (LUPRI,'(//A,I5)')
     *     ' LRST, START VECTOR IS NOT LINEAR INDEPENDENT.',KEXSTV
         IF (KZRED.LE.0) THEN
            CALL QUIT(' LRST, START VECTORS NOT LINEAR INDEPENDENT.')
         END IF
      END IF
C
C     END OF LRST.
C
      RETURN
      END
C  /* Deck lrsopr */
      SUBROUTINE LRSOPR(NOP,SOPR,FC,CMO,UDV,PV,XINDX,LURSP2,WRK,LWRK)
C
C CALCULATE SECOND ORDER PROPERTIES.
C QUADRATIC ACCURACY IS OBTAINED WITH SELLERS FORMULA
C
C     SOPR(I,J) = GP(I)*SOL(J) - SOL(I)*RES(J)
C
C GP(I)     : PROPERTY GRADIENT
C SOL(I)    : SOLUTION VECTOR
C RES(I)    : RESIDUE
C SOPR(I,J) : SECOND ORDER PROPERTY FOR PROPERTY I,J
C
#include <implicit.h>
#include <iratdef.h>
C
#include <codata.h>
      PARAMETER ( DTOL   = 1.0D-10 )
C
#include <priunit.h>
#include <rspprp.h>
#include <inflr.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <infpri.h>
#include <mxcent.h>
#include <gtensor.h>
#include <elweak.h>
#include <dummy.h>
C
      DIMENSION FC(*),CMO(*),UDV(*),PV(*),XINDX(*),WRK(*)
      DIMENSION SOPR(NFREQ,NOP,*)
C
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KISOL,KZYVAR,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KJSOL,KZYVAR,WRK,KFREE,LFREE)
      KJRES = KJSOL
      KGP   = KFREE
      LWRKGP= LWRK  - KGP
      KWRK1 = KGP   + KZYVAR
      IF (KWRK1 .GT. LWRK) CALL ERRWRK('LRSOPR',-(KWRK1-1),LWRK)
C
      CALL HEADER('Final output of second order properties from'//
     &            ' linear response',-1)
      IF (TRPLET) THEN
         WRITE (LUPRI,'(/A)') '@ Spin symmetry of operators: triplet'
      ELSE
         WRITE (LUPRI,'(/A)') '@ Spin symmetry of operators: singlet'
      END IF
      WRITE (LUPRI,'(/A/A)')
     &   ' Note that minus the linear response function:'//
     &   ' - << A; B >>(omega) is printed.',
     &   ' The results are of quadratic accuracy using Sellers formula.'
C
      DO 100 IOP = 1,NOP
         IF (SOPRSY) THEN
            JST = 1
         ELSE
            JST = IOP
         END IF
         CALL GETGPV(LBLLR(KSYMOP,IOP),FC,DUMMY,CMO,UDV,PV,XINDX,ANTSYM,
     *               WRK(KGP),LWRKGP)
         DO 200 IFREQ = 1,NFREQ
            IFRSOL = (IOP-1)*NFREQ*2 + (IFREQ-1)*2 + 1
            CALL READDX(LURSP2,IFRSOL,IRAT*KZYVAR,WRK(KISOL))
            DO 300 JOP = JST,NOP
               JFRSOL = (JOP-1)*NFREQ*2 + (IFREQ-1)*2 + 1
               CALL READDX(LURSP2,JFRSOL,IRAT*KZYVAR,WRK(KJSOL))
               SOPR(IFREQ,IOP,JOP) =
     *            DDOT(KZYVAR,WRK(KGP),1,WRK(KJSOL),1)
               IF (IOP.NE.JOP) THEN
                  IF (SOPRSY) THEN
                     WRITE(LUPRI,'(/A/A,1P,D15.7,A/,A,3X,A,D15.7)')
     *               ' TEST OUTPUT, linear accuracy',
     *               ' FREQUENCY = ',FREQ(IFREQ),' au',
     *               LBLLR(KSYMOP,IOP),
     *               LBLLR(KSYMOP,JOP),SOPR(IFREQ,IOP,JOP)
                  END IF
                  JFRRES = JFRSOL + 1
                  CALL READDX(LURSP2,JFRRES,IRAT*KZYVAR,WRK(KJRES))
                  SOPR(IFREQ,IOP,JOP) = SOPR(IFREQ,IOP,JOP)
     *                 - DDOT(KZYVAR,WRK(KISOL),1,WRK(KJRES),1)
               END IF
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
      DO 500 IFREQ = 1,NFREQ
         IF ( FREQ(IFREQ).LT.DTOL ) THEN
            WRITE(LUPRI,'(/A/)')
     *      '@ FREQUENCY INDEPENDENT SECOND ORDER PROPERTIES'
         ELSE
            WRITE(LUPRI,'(/A/,5(/A,1P,D15.7),/)')
     * '@ FREQUENCY DEPENDENT SECOND ORDER PROPERTIES WITH FREQUENCY :',
     *      '@    a.u.:',FREQ(IFREQ),
     *      '     cm-1:',XTKAYS*FREQ(IFREQ),
     *      '     e.V.:',XTEV*FREQ(IFREQ),
     *      '  kJ/mol.:',XKJMOL*FREQ(IFREQ),
     *      '     nm  :',XTNM/FREQ(IFREQ)
         END IF

         DO 400 IOP = 1,NOP
            IF (SOPRSY) THEN
               JST = 1
            ELSE
               JST = IOP
            END IF
            IF (SOPW4) THEN
               CALL GETGPV(LBLLR(KSYMOP,IOP),FC,DUMMY,CMO,UDV,PV,XINDX,
     &            ANTSYM,WRK(KGP),LWRKGP)
               CALL SOPDW4(WRK(KGP),WRK(KISOL),WRK(KJSOL),FREQ(IFREQ))
            END IF
c ach
            DO 600 JOP = JST,NOP
               XPROP = SOPR(IFREQ,IOP,JOP)
               WRITE(LUPRI,'(5A,1P,D20.12)')
     * '@ -<< ',LBLLR(KSYMOP,IOP),' ; ',LBLLR(KSYMOP,JOP),' >> =',XPROP
C JK
                CALL WRIPRO(XPROP,"LR-SCF/DFT   ",2,
     *                      LBLLR(KSYMOP,IOP),LBLLR(KSYMOP,JOP),
     *                      LBLLR(KSYMOP,IOP),LBLLR(KSYMOP,JOP),
     *                      FREQ(IFREQ),FREQ(IFREQ),FREQ(IFREQ),
     *                      1,0,0,0)
C JK
               IF (SOPW4) THEN
                  CALL GETGPV(LBLLR(KSYMOP,JOP),FC,DUMMY,CMO,UDV,PV,
     &                        XINDX,ANTSYM,WRK(KGP),LWRKGP)
                  W4TERM = DDOT(KZYVAR,WRK(KGP),1,WRK(KISOL),1)
                  WRITE(LUPRI,'(A,1P,D20.12)')
     &            '@ SOPPA W4(E) term  ',W4TERM
               END IF
               IF (GCALC) THEN
                  CALL GOZHSO(IOP,JOP,XPROP)
               END IF
 600        CONTINUE
 400     CONTINUE
         if (elweak) call printelweak(nop,sopr,ksymop)
 500  CONTINUE
      RETURN
      END
C  /* Deck orbdia */
      SUBROUTINE ORBDIA(FOCK,FC,FV,UDV,EODIA,SODIA)
C
C WRITTEN 13-MAY 1987
C
C PURPOSE:
C    CALCULATE FOCK CONTRIBUTIONS TO DIAGONAL E(2) AND S(2) MATRICES
C
C       EODIA(L,K) = < [E(K,L),H,E(L,K)] > L>K
C
C       SODIA(L,K) = < [E(K,L),E(L,K)] >   L>K
C
C    SECONDARY - INACTIVE
C     (L,K) = (A,I)
C    EODIA(A,I) = -FOCK(I,I) + 2*FC(A,A) +2FV(A,A)
C    SODIA(A,I) = 2
C
C    SECONDARY - ACTIVE
C             (A,M)
C    EODIA(A,M) = UDV(M,M)*FC(A,A) - FOCK(M,M)
C    SODIA(A,M) = UDV(M,M)
C
C    ACTIVE - INACTIVE
C             (M,I)
C    EODIA(M,I) = 2*FC(M,M) + UDV(M,M)*FC(I,I) - FOCK(I,I) -FOCK(M,M) +
C                2*FV(M,M)
C    SODIA(M,I) = 2 - UDV(M,M)
C
C    ACTIVE - ACTIVE
C             (M,N)
C    EODIA(M,N) = UDV(N,N)*FC(M,M) + UDV(M,M)*FC(N,N) - FOCK(M,M) -FOCK(N,N)
C    SODIA(M,N) = UDV(N,N) - UDV(M,M)
C
C    AVDIA = .TRUE. , ADD FV CONTRIBUTIONS TO EODIA
C            WHICH ORIGINATE FROM FOCK TYPE DECOUPLING OF THE TWO
C            ELECTRON DENSITY MATRIX.
C            ALL UDV(*,*)*FC(*,*) CONTRIBUTIONS THEN BECOME
C                UDV(*,*)*(FC(*,*)+FV(*,*))
C
C                ******************************
C
#include <implicit.h>
C
C
      DIMENSION FC(*),FV(*),FOCK(*),UDV(NASHDI,NASHDI)
      DIMENSION EODIA(*),SODIA(*)
C
      PARAMETER ( D0 = 0.0D0 )
C  INFDIM : NASHDI
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <infvar.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infpri.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
C -- local constants
C
      PARAMETER ( D2 = 2.0D0 )
C
      CALL DZERO(EODIA,KZWOPT)
      CALL DZERO(SODIA,KZWOPT)
      KSYM1 = 0
      DO 100 IG = 1,KZWOPT
         K     = JWOP(1,IG)
         L     = JWOP(2,IG)
         KSYM  = ISMO(K)
         LSYM  = ISMO(L)
         IF ( KSYM.NE.KSYM1 ) THEN
            KSYM1 = KSYM
            NORBK = NORB(KSYM)
            IORBK = IORB(KSYM)
            IIORBK= IIORB(KSYM)
            I2ORBK= I2ORB(KSYM)
            IORBL = IORB(LSYM)
            NORBL = NORB(LSYM)
            IIORBL= IIORB(LSYM)
            I2ORBL= I2ORB(LSYM)
         END IF
         NK    = K - IORBK
         NL    = L - IORBL
         ITYPK = IOBTYP(K)
         ITYPL = IOBTYP(L)
         IIKK  = IIORBK + NK*(NK+1)/2
         IILL  = IIORBL + NL*(NL+1)/2
         I2LL  = I2ORBL+NORBL*(NL-1)+NL
         I2KK  = I2ORBK+NORBK*(NK-1)+NK
         IF ( ITYPK.EQ.JTINAC )THEN
            EODIA(IG) = EODIA(IG) - FOCK(I2KK) + D2*FC(IILL) +
     *                  D2*FV(IILL)
            SODIA(IG) = SODIA(IG) + D2
            IF ( ITYPL.EQ.JTACT ) THEN
               NWL = ISW(L) - NISHT
               EODIA(IG) = EODIA(IG) - FOCK(I2LL) +
     *                     UDV(NWL,NWL)*FC(IIKK)
               IF (AVDIA) THEN
                  EODIA(IG) = EODIA(IG)+ UDV(NWL,NWL)*FV(IIKK)
               END IF
               SODIA(IG) = SODIA(IG) - UDV(NWL,NWL)
            END IF
         ELSE IF ( ITYPK .EQ. JTACT ) THEN
            IF (ITYPL.EQ.JTACT) THEN
               NWL = ISW(L) - NISHT
               NWK = ISW(K) - NISHT
               EODIA(IG) = EODIA(IG) + UDV(NWK,NWK)*FC(IILL)+
     *                   UDV(NWL,NWL)*FC(IIKK) - FOCK(I2LL) - FOCK(I2KK)
               IF (AVDIA) THEN
                  EODIA(IG) = EODIA(IG) + UDV(NWK,NWK)*FV(IILL) +
     *                               UDV(NWL,NWL)*FV(IIKK)
               END IF
               SODIA(IG) = SODIA(IG) + UDV(NWK,NWK) - UDV(NWL,NWL)
            ELSE
               NWK = ISW(K) - NISHT
               EODIA(IG) = EODIA(IG) + UDV(NWK,NWK)*FC(IILL) -
     *                     FOCK(I2KK)
               IF (AVDIA) THEN
                  EODIA(IG) = EODIA(IG) + UDV(NWK,NWK)*FV(IILL)
               END IF
               SODIA(IG) = SODIA(IG) + UDV(NWK,NWK)
            ENDIF
         ELSE
            WRITE (LUPRI,*) 'ORBDIA FATAL ERROR: k not occ. for '//
     &                      'k->l rot.'
            WRITE (LUPRI,*) 'IG, K, L, ITYPK, ITYPL',IG,K,L,ITYPK,ITYPL
            WRITE (LUPRI,*) 'Possible cause: inconsistent common '//
     *         'blocks, check /INFIND/ and MXCORB'
            CALL QUIT('ORBDIA FATAL ERROR: k not occ. for k->l rot.')
         ENDIF
 100  CONTINUE
C
C *** Perform averaging of EODIA and SODIA if necessary
C
      IF (RSPSUP .AND. (KSYMOP .EQ. 1)) THEN
         CALL RSPAVE(EODIA,KZWOPT,1)
         CALL RSPAVE(SODIA,KZWOPT,1)
      END IF
C
C CHECK IF DIAGONAL ORBITAL HESSIAN HAS NEGATIVE ELEMENT
C ADD ORBSFT TO REMOVE SPURIOUS NEGATIVE ELEMENTS
C
      NNEG = 0
      DO 300 I = 1,KZWOPT
         EODIA(I) = EODIA(I) +  ORBSFT
         IF ( EODIA(I).LT.D0 ) NNEG = NNEG + 1
 300  CONTINUE
      IF (NNEG.GT.0) THEN
         NWARN = NWARN + 1
         WRITE(LUPRI,'(/A,/A,I5,A)')
     *   '  ********* WARNING ORBDIA *********** ',
     *   '  DIAGONAL ORBITAL HESSIAN HAS ',NNEG,'  NEGATIVE ELEMENTS'
      END IF
      IF (IPRRSP.GT.101) THEN
         WRITE(LUPRI,*)' DIAGONAL ORBITAL PART OF E(2) '
         CALL OUTPUT(EODIA,1,KZWOPT,1,1,KZWOPT,1,1,LUPRI)
         WRITE(LUPRI,*)' DIAGONAL ORBITAL PART OF S(2)'
         CALL OUTPUT(SODIA,1,KZWOPT,1,1,KZWOPT,1,1,LUPRI)
      END IF
C
      RETURN
      END
C  /* Deck qonedi */
      SUBROUTINE QONEDI(NOSIM,NCIW,NDIW,ICDSYM,QAONE,QBONE,ZYMAT,H2,
     *                  H2X,PVX,PVCD2,PVCD3,PVDC2,PVDC3,WRK,LWRK)
C
C PURPOSE:
C  CALCULATE THE CONTRIBUTIONS TO QA AND QB MATRICES WHICH
C  ARE DETERMINED FROM ACTIVE-ACTIVE DIRAC INTEGRAL DISTRIBUTION
C  <AB,CD> C>D (DIRAC NOTATION).
C  CONTRIBUTIONS FROM MULLIKEN DISTRIBUTION ARE CALCULATED IN QONEMU
C
#include <implicit.h>
C
      DIMENSION QAONE(NORBT,NASHDI,*),QBONE(NORBT,NASHDI,*)
      DIMENSION ZYMAT(NORBT,NORBT,*),H2(NORBT,*),H2X(NORBT,*)
      DIMENSION PVX(*),PVCD2(NASHDI,*),PVCD3(NASHDI,*),PVDC2(NASHDI,*)
      DIMENSION PVDC3(NASHDI,*),WRK(*)
C
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infdim.h>
#include <infind.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
      PARAMETER ( DM1 = -1.0D0 )
C
      IF (TRPLET) THEN
         IPOIDI = 1 + N2ASHX*N2ASHX
      ELSE
         IPOIDI = 1
      END IF
      CALL PVXDIS(NCIW,NDIW,PVCD2,PVX(IPOIDI),2)
      CALL PVXDIS(NCIW,NDIW,PVCD3,PVX(IPOIDI),3)
      IF (NCIW.NE.NDIW) THEN
         CALL PVXDIS(NDIW,NCIW,PVDC2,PVX(IPOIDI),2)
         CALL PVXDIS(NDIW,NCIW,PVDC3,PVX(IPOIDI),3)
      END IF
C
C     CALL PVXDIS(K,L,PVDEN,PVX,IPVDIS)
C
      DO 1000 IOSIM = 1,NOSIM
C
         CALL DZERO(H2X,NORBT*NORBT)
         DO 1100 IASYM = 1,NSYM
            IBSYM = MULD2H(IASYM,ICDSYM)
            IPSYM = MULD2H(IASYM,KSYMOP)
            NORBA = NORB(IASYM)
            NORBB = NORB(IBSYM)
            NORBP = NORB(IPSYM)
            IOFFA = IORB(IASYM) + 1
            IOFFB = IORB(IBSYM) + 1
            IOFFP = IORB(IPSYM) + 1
C
C ONE INDEX TRANSFORM FIRST INTEGRAL INDEX
C
C  <P~B,CD> =  SUM(A) ZYMAT(P,A)*<AB,CD>
C
            IF ( (NORBA.GT.0) .AND. (NORBB.GT.0) .AND. (NORBP.GT.0) )
     *         CALL DGEMM('N','N',NORBP,NORBB,NORBA,1.D0,
     &                    ZYMAT(IOFFP,IOFFA,IOSIM),NORBT,
     &                    H2(IOFFA,IOFFB),NORBT,1.D0,
     &                    H2X(IOFFP,IOFFB),NORBT)
 1100    CONTINUE
C
C
         DO 1200 IQSYM = 1,NSYM
            IPSYM = MULD2H(IQSYM,KSYMOP)
            IYSYM = MULD2H(IQSYM,ICDSYM)
            NORBP = NORB(IPSYM)
            NASHQ = NASH(IQSYM)
            NASHY = NASH(IYSYM)
            IF ((NORBP.GT.0).AND.(NASHQ.GT.0).AND.(NASHY.GT.0)) THEN
               IOFFP = IORB(IPSYM) + 1
               IOFFY = IORB(IYSYM) + 1
               NISHQ = NISH(IQSYM)
               NISHY = NISH(IYSYM)
               ISTY  = ISW( IORB(IYSYM) + NISH(IYSYM) +1 ) - NISHT
               ISTQ  = ISW( IORB(IQSYM) + NISH(IQSYM) +1 ) - NISHT
C
C ADD SUM(Y) <Y~P,CD> * [YC,QD] TO QBONE(P,Q)
C
               CALL DGEMM('T','N',NORBP,NASHQ,NASHY,1.D0,
     &                    H2X(IOFFY+NISHY,IOFFP),NORBT,
     &                    PVCD2(ISTY,ISTQ),NASHT,1.D0,
     &                    QBONE(IOFFP,ISTQ,IOSIM),NORBT)
C
C ADD SUM(X) <Y~P,CD> * [YC,DQ] TO QAONE(P,Q)
C
               CALL DGEMM('T','N',NORBP,NASHQ,NASHY,1.D0,
     &                    H2X(IOFFY+NISHY,IOFFP),NORBT,
     &                    PVCD3(ISTY,ISTQ),NASHT,1.D0,
     &                    QAONE(IOFFP,ISTQ,IOSIM),NORBT)
            END IF
 1200    CONTINUE
         IF (NCIW.NE.NDIW) THEN
            CALL DZERO(H2X,NORBT*NORBT)
            DO 1400 IASYM = 1,NSYM
               IBSYM = MULD2H(IASYM,ICDSYM)
               IPSYM = MULD2H(IASYM,KSYMOP)
               NORBA = NORB(IASYM)
               NORBB = NORB(IBSYM)
               NORBP = NORB(IPSYM)
               IOFFA = IORB(IASYM) + 1
               IOFFB = IORB(IBSYM) + 1
               IOFFP = IORB(IPSYM) + 1
C
C ONE INDEX TRANSFORM FIRST INTEGRAL INDEX
C
C  <P~B,DC> =  SUM(A) ZYMAT(P,A)*<BA,CD>
C
               IF ( (NORBA.GT.0) .AND. (NORBB.GT.0) .AND. (NORBP.GT.0) )
     *             CALL DGEMM('N','T',NORBP,NORBB,NORBA,1.D0,
     &                        ZYMAT(IOFFP,IOFFA,IOSIM),NORBT,
     &                        H2(IOFFB,IOFFA),NORBT,1.D0,
     &                        H2X(IOFFP,IOFFB),NORBT)
 1400       CONTINUE
C
C
            DO 1500 IQSYM = 1,NSYM
               IPSYM = MULD2H(IQSYM,KSYMOP)
               IYSYM = MULD2H(IQSYM,ICDSYM)
               NORBP = NORB(IPSYM)
               NASHQ = NASH(IQSYM)
               NASHY = NASH(IYSYM)
               IF ( (NORBP.GT.0) .AND. (NASHQ.GT.0) .AND. (NASHY.GT.0) )
     *                                                             THEN
                  IOFFP = IORB(IPSYM) + 1
                  IOFFY = IORB(IYSYM) + 1
                  NISHQ = NISH(IQSYM)
                  NISHY = NISH(IYSYM)
                  ISTY  = ISW( IORB(IYSYM) + NISH(IYSYM) +1 ) - NISHT
                  ISTQ  = ISW( IORB(IQSYM) + NISH(IQSYM) +1 ) - NISHT
C
C ADD SUM(Y) <Y~P,DC> * [YD,QC] TO QBONE(P,Q)
C
                  CALL DGEMM('T','N',NORBP,NASHQ,NASHY,1.D0,
     &                       H2X(IOFFY+NISHY,IOFFP),NORBT,
     &                       PVDC2(ISTY,ISTQ),NASHT,1.D0,
     &                       QBONE(IOFFP,ISTQ,IOSIM),NORBT)
C
C ADD SUM(X) <Y~P,DC> * [YD,CQ] TO QAONE(P,Q)
C
                  CALL DGEMM('T','N',NORBP,NASHQ,NASHY,1.D0,
     &                       H2X(IOFFY+NISHY,IOFFP),NORBT,
     &                       PVDC3(ISTY,ISTQ),NASHT,1.D0,
     &                       QAONE(IOFFP,ISTQ,IOSIM),NORBT)
               END IF
 1500       CONTINUE
         END IF
C
C ADD ONE INDEX TRANSFORMATION ON SECOND INTEGRAL INDEX
C
         CALL DZERO(H2X,NORBT*NORBT)
         DO 1600 IASYM = 1,NSYM
            IBSYM = MULD2H(IASYM,ICDSYM)
            IPSYM = MULD2H(IASYM,KSYMOP)
            IOFFA = IORB(IASYM) + 1
            IOFFB = IORB(IBSYM) + 1
            IOFFP = IORB(IPSYM) + 1
            NORBA = NORB(IASYM)
            NORBB = NORB(IBSYM)
            NORBP = NORB(IPSYM)
C
C  <CD,P~B> = -SUM(A) ZYMAT(A,P)*<CD,AB>
C
            IF ( (NORBA.GT.0) .AND. (NORBB.GT.0) .AND. (NORBP.GT.0) )
     *         CALL DGEMM('T','N',NORBP,NORBB,NORBA,1.D0,
     &                    ZYMAT(IOFFA,IOFFP,IOSIM),NORBT,
     &                    H2(IOFFA,IOFFB),NORBT,1.D0,
     &                    H2X(IOFFP,IOFFB),NORBT)
 1600    CONTINUE
         CALL DSCAL(N2ORBX,DM1,H2X,1)
         DO 1700 IQSYM = 1,NSYM
            IPSYM = MULD2H(IQSYM,KSYMOP)
            IYSYM = MULD2H(IQSYM,ICDSYM)
            NORBP = NORB(IPSYM)
            NASHQ = NASH(IQSYM)
            NASHY = NASH(IYSYM)
            IF ( (NORBP.GT.0) .AND. (NASHQ.GT.0) .AND. (NASHY.GT.0) )
     *                                                            THEN
               IOFFP = IORB(IPSYM) + 1
               IOFFY = IORB(IYSYM) + 1
               NISHQ = NISH(IQSYM)
               NISHY = NISH(IYSYM)
               ISTY  = ISW( IORB(IYSYM) + NISH(IYSYM) +1 ) - NISHT
               ISTQ  = ISW( IORB(IQSYM) + NISH(IQSYM) +1 ) - NISHT
C
C ADD SUM(Y) <CD,Y~P> * [YC,DQ] TO QBONE(P,Q)
C
               CALL DGEMM('T','N',NORBP,NASHQ,NASHY,1.D0,
     &                    H2X(IOFFY+NISHY,IOFFP),NORBT,
     &                    PVCD3(ISTY,ISTQ),NASHT,1.D0,
     &                    QBONE(IOFFP,ISTQ,IOSIM),NORBT)
C
C ADD SUM(Y) <CD,Y~P> * [YC,QD] TO QAONE(P,Q)
C
               CALL DGEMM('T','N',NORBP,NASHQ,NASHY,1.D0,
     &                    H2X(IOFFY+NISHY,IOFFP),NORBT,
     &                    PVCD2(ISTY,ISTQ),NASHT,1.D0,
     &                    QAONE(IOFFP,ISTQ,IOSIM),NORBT)
            END IF
 1700    CONTINUE
         IF (NCIW.NE.NDIW) THEN
         CALL DZERO(H2X,NORBT*NORBT)
            DO 1800 IASYM = 1,NSYM
               IBSYM = MULD2H(IASYM,ICDSYM)
               IPSYM = MULD2H(IASYM,KSYMOP)
               IOFFA = IORB(IASYM) + 1
               IOFFB = IORB(IBSYM) + 1
               IOFFP = IORB(IPSYM) + 1
               NORBA = NORB(IASYM)
               NORBB = NORB(IBSYM)
               NORBP = NORB(IPSYM)
C
C  <DC,P~B> = -SUM(A) ZYMAT(A,P)*<CD,BA>
C
               IF ( (NORBA.GT.0) .AND. (NORBB.GT.0) .AND. (NORBP.GT.0) )
     *            CALL DGEMM('T','T',NORBP,NORBB,NORBA,1.D0,
     &                       ZYMAT(IOFFA,IOFFP,IOSIM),NORBT,
     &                       H2(IOFFB,IOFFA),NORBT,0.D0,
     &                       H2X(IOFFP,IOFFB),NORBT)
 1800       CONTINUE
            CALL DSCAL(N2ORBX,DM1,H2X,1)
            DO 1900 IQSYM = 1,NSYM
               IPSYM = MULD2H(IQSYM,KSYMOP)
               IYSYM = MULD2H(IQSYM,ICDSYM)
               NORBP = NORB(IPSYM)
               NASHQ = NASH(IQSYM)
               NASHY = NASH(IYSYM)
               IF ( (NORBP.GT.0) .AND. (NASHQ.GT.0) .AND. (NASHY.GT.0) )
     *                                                            THEN
                  IOFFP = IORB(IPSYM) + 1
                  IOFFY = IORB(IYSYM) + 1
                  NISHQ = NISH(IQSYM)
                  NISHY = NISH(IYSYM)
                  ISTY  = ISW( IORB(IYSYM) + NISH(IYSYM) +1 ) - NISHT
                  ISTQ  = ISW( IORB(IQSYM) + NISH(IQSYM) +1 ) - NISHT
C
C ADD SUM(Y) <DC,Y~P> * [YD,CQ] TO QBONE(P,Q)
C
                  CALL DGEMM('T','N',NORBP,NASHQ,NASHY,1.D0,
     &                       H2X(IOFFY+NISHY,IOFFP),NORBT,
     &                       PVDC3(ISTY,ISTQ),NASHT,1.D0,
     &                       QBONE(IOFFP,ISTQ,IOSIM),NORBT)
C
C ADD SUM(Y) <DC,Y~P> * [YD,QC] TO QAONE(P,Q)
C
                  CALL DGEMM('T','N',NORBP,NASHQ,NASHY,1.D0,
     &                       H2X(IOFFY+NISHY,IOFFP),NORBT,
     &                       PVDC2(ISTY,ISTQ),NASHT,1.D0,
     &                       QAONE(IOFFP,ISTQ,IOSIM),NORBT)
               END IF
 1900       CONTINUE
         END IF
 1000 CONTINUE
C
C END OF QONEDI
C
      RETURN
      END
C  /* Deck qonemu */
      SUBROUTINE QONEMU(NOSIM,NCIW,NDIW,ICDSYM,QAONE,QBONE,
     *                  ZYMAT,H2,H2X,H2XAC,PVX,PVCD,WRK,LWRK)
C
C PURPOSE:
C  CALCULATE THE CONTRIBUTIONS TO QA AND QB MATRICES WHICH
C  ARE DETERMINED FROM ACTIVE-ACTIVE MULLIKEN INTEGRAL DISTRIBUTION
C  (**,CD) C>D .
C  CONTRIBUTIONS FROM DIRAC DISTRIBUTION XY CALCULATED IN QONEDI
C
#include <implicit.h>
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infpri.h>
C
      DIMENSION QAONE(NORBT,NASHT,*),QBONE(NORBT,NASHT,*)
      DIMENSION ZYMAT(NORBT,NORBT,*),PVX(*),PVCD(NASHT,*),H2(NORBT,*)
      DIMENSION H2X(NORBT,*),H2XAC(N2ASHX,NNASHX,*),WRK(*)
      IF (SOPPA) THEN
         WRITE(LUPRI,'(/A)') 'ERROR: QONEMU called in SOPPA'
         CALL QTRACE(LUERR)
         CALL QUIT('ERROR: QONEMU called in SOPPA')
      END IF
      DO 1000 IOSIM = 1,NOSIM
         CALL DZERO(H2X,NORBT*NORBT)
         DO 1100 IASYM = 1,NSYM
            IBSYM = MULD2H(IASYM,ICDSYM)
            IPSYM = MULD2H(IASYM,KSYMOP)
            IQSYM = MULD2H(IBSYM,KSYMOP)
            IOFFA = IORB(IASYM) + 1
            IOFFB = IORB(IBSYM) + 1
            IOFFP = IORB(IPSYM) + 1
            IOFFQ = IORB(IQSYM) + 1
            NORBA = NORB(IASYM)
            NORBB = NORB(IBSYM)
            NORBP = NORB(IPSYM)
            NORBQ = NORB(IQSYM)
C
C ONE INDEX TRANSFORM FIRST INTEGRAL INDEX
C
C  (P~B,CD) =  SUM(A) ZYMAT(P,A)*(AB,CD)
C
            IF ((NORBP.NE.0) .AND. (NORBA.NE.0) .AND. (NORBB.NE.0))
     *      CALL DGEMM('N','N',NORBP,NORBB,NORBA,1.D0,
     &                 ZYMAT(IOFFP,IOFFA,IOSIM),NORBT,
     &                 H2(IOFFA,IOFFB),NORBT,1.D0,
     &                 H2X(IOFFP,IOFFB),NORBT)
C
C ADD ONE INDEX TRANSFORMATION ON SECOND INTEGRAL INDEX
C
C  (AQ~,CD) = -SUM(B) (AB,CD)*ZYMAT(B,Q)
C
        IF ((NORBB.NE.0) .AND. (NORBA.NE.0) .AND. (NORBQ.NE.0))
     *      CALL DGEMM('N','N',NORBA,NORBQ,NORBB,-1.D0,
     &                 H2(IOFFA,IOFFB),NORBT,
     &                 ZYMAT(IOFFB,IOFFQ,IOSIM),NORBT,1.D0,
     &                 H2X(IOFFA,IOFFQ),NORBT)
 1100    CONTINUE
         IF ( IPRRSP.GT.150 ) THEN
            WRITE(LUPRI,'(/A,2I5)')
     *       ' NON-TRANSFORMED INTEGRALS: DISTRBR: NCIW,NDIW',NCIW,NDIW
             CALL OUTPUT(H2,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            WRITE(LUPRI,'(/A,2I5)')
     *       ' HALF TRANSFORMED INTEGRALS: DISTRBR: NCIW,NDIW',NCIW,NDIW
             CALL OUTPUT(H2X,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
C
C DISTRIBUTE ACTIVE HALF ONE INDEX TRANSFORMED INTEGRALS INTO H2XAC
C
         IF (KZCONF.GT.0) THEN
            IF (NCIW.GE.NDIW) THEN
               NCDW = NCIW*(NCIW-1)/2 + NDIW
            ELSE
               NCDW = NDIW*(NDIW-1)/2 + NCIW
            END IF
            DO 1200 IBW = 1,NASHT
               IB = ISX(NISHT+IBW)
               DO 1300 IAW = 1,NASHT
                  IA = ISX(NISHT+IAW)
                  IABW = (IBW-1) * NASHT + IAW
                  H2XAC(IABW,NCDW,IOSIM) = H2XAC(IABW,NCDW,IOSIM) +
     *                                     H2X(IA,IB)
 1300          CONTINUE
 1200       CONTINUE
         IF ( IPRRSP.GT.150 ) THEN
            WRITE(LUPRI,'(/A,I5)')
     *       ' HALF TRANSFORMED ACTIVE INTEGRALS: DISTRB: NCDW',NCDW
             CALL OUTPUT(H2XAC(1,NCDW,IOSIM),1,NASHT,1,NASHT,
     *                   NASHT,NASHT,1,LUPRI)
         END IF
         END IF
C
C GET DENSITY DISTRIBUTION
C
      CALL PVXDIS(NCIW,NDIW,PVCD,PVX,1)
C     CALL PVXDIS(K,L,PVDEN,PVX,IPVDIS)
C
C ADD CONTRIBUTIONS TO QAONE AND QBONE
C
         DO 3000 IQSYM = 1,NSYM
            IPSYM = MULD2H(IQSYM,KSYMOP)
            IYSYM = MULD2H(IQSYM,ICDSYM)
            NORBP = NORB(IPSYM)
            NASHQ = NASH(IQSYM)
            NASHY = NASH(IYSYM)
            IF ( (NORBP.GT.0) .AND. (NASHQ.GT.0) .AND. (NASHY.GT.0) )
     *                                                           THEN
               IOFFP = IORB(IPSYM) + 1
               IOFFY = IORB(IYSYM) + 1
               NISHY = NISH(IYSYM)
               ISTY  = ISW( IORB(IYSYM) + NISH(IYSYM) +1 ) - NISHT
               ISTQ  = ISW( IORB(IQSYM) + NISH(IQSYM) +1 ) - NISHT
C
C ADD SUM(Y) (P~Y~,CD)*{[QY,CD]+(1-DELTA(C,D)*[YQ,DC]} TO QBONE(P,Q)
C
               CALL DGEMM('N','T',NORBP,NASHQ,NASHY,1.D0,
     &                    H2X(IOFFP,IOFFY+NISHY),NORBT,
     &                    PVCD(ISTQ,ISTY),NASHT,1.D0,
     &                    QBONE(IOFFP,ISTQ,IOSIM),NORBT)
C
C ADD SUM(Y) (Y~P~,CD)*{[YQ,CD]+(1-DELTA(D,C)*[YQ,DC]} TO QAONE(P,Q)
C
               CALL DGEMM('T','N',NORBP,NASHQ,NASHY,1.D0,
     &                    H2X(IOFFY+NISHY,IOFFP),NORBT,
     &                    PVCD(ISTY,ISTQ),NASHT,1.D0,
     &                    QAONE(IOFFP,ISTQ,IOSIM),NORBT)
            END IF
 3000    CONTINUE
C
 1000 CONTINUE
C
C END OF QONEMU
C
      RETURN
      END
C  /* Deck qtd */
      SUBROUTINE QTD(NCSIM,NCIW,NDIW,ICDSYM,QATD,QBTD,H2,PVTD,PVDEN,
     *               WRK,LWRK)
C
C PURPOSE:
C  CALCULATE THE CONTRIBUTIONS TO QATD AND QBTD MATRICES WHICH
C  ARE DETERMINED FROM INTEGRAL DISTRIBUTION  (AB,CD) ,D>=C
C  (MULLIKEN NOTATION)
C
#include <implicit.h>
C
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infdim.h>
#include <infind.h>
#include <wrkrsp.h>
C
      DIMENSION QATD(NORBT,NASHDI,*),QBTD(NORBT,NASHDI,*)
      DIMENSION PVTD(N2ASHX,N2ASHX,*),PVDEN(NASHT,*)
      DIMENSION H2(NORBT,NORBT),WRK(*)
      DO 1000 ICSIM = 1,NCSIM
         CALL PVXDIS(NCIW,NDIW,PVDEN,PVTD(1,1,ICSIM),1)
C        CALL PVXDIS(K,L,PVDEN,PVX,IPVDIS)
         DO 3000 IQSYM = 1,NSYM
            IPSYM = MULD2H(IQSYM,KSYMOP)
            IYSYM = MULD2H(ICDSYM,IPSYM)
            NORBP = NORB(IPSYM)
            NASHQ = NASH(IQSYM)
            NASHY = NASH(IYSYM)
            IF ( (NORBP.GT.0) .AND. (NASHQ.GT.0) .AND. (NASHY.GT.0) )
     *                                                          THEN
               IOFFP = IORB(IPSYM)
               IOFFY = IORB(IYSYM)
               NISHY = NISH(IYSYM)
               ISTQ  = ISW( IORB(IQSYM) + NISH(IQSYM) +1 ) - NISHT
               ISTY  = ISW( IORB(IYSYM) + NISH(IYSYM) +1 ) - NISHT
C
C ADD SUM(Y) (PY,CD)*{[QY,CD]+(1-DELTA(D,C)*[QY,DC]} TO QBTD(P,Q)
C
               CALL DGEMM('N','T',NORBP,NASHQ,NASHY,1.D0,
     &                    H2(IOFFP+1,IOFFY+NISHY+1),NORBT,
     &                    PVDEN(ISTQ,ISTY),NASHT,1.D0,
     &                    QBTD(IOFFP+1,ISTQ,ICSIM),NORBT)
C
C ADD SUM(Y) (PY,CD)*{[YQ,CD]+(1-DELTA(D,C)*[YQ,DC]} TO QATD(P,Q)
C
               CALL DGEMM('N','N',NORBP,NASHQ,NASHY,1.D0,
     &                    H2(IOFFP+1,IOFFY+NISHY+1),NORBT,
     &                    PVDEN(ISTY,ISTQ),NASHT,1.D0,
     &                    QATD(IOFFP+1,ISTQ,ICSIM),NORBT)
            END IF
 3000    CONTINUE
 1000 CONTINUE
C
C END OF QTD
C
      RETURN
      END
C  /* Deck rspcli */
      SUBROUTINE RSPCLI(NCSIM,ZYCVEC,FCAC,H2AC,
     *                  EVECS,XINDX,WRK,LWRK)
C
C CALCULATE CONFIGURATION PART OF LINEAR TRANSFORMATION
C
C                ( <J,H,0(R)> )          ( S(J) )
C     E[2]*N = - (            ) + EMCSCF*(      )
C                (-<0(L),H,J> )          ( S'(J))
C
C     |0(R)> = -SUM(N) S(N) *|N>
C     |0(L)> =  SUM(N) S'(N)*|N>
C
C MODIFIED 5-5-1988 : ASSUME |0(L)> = 0 i.e. S' = 0 AND
C                     NOT STORED'
#include <implicit.h>
C
      DIMENSION ZYCVEC(KZCONF,*),FCAC(*),H2AC(*)
      DIMENSION EVECS(KZYVAR,*),XINDX(*),WRK(*)
C
#include <priunit.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infopt.h>
#include <inforb.h>
#include <infpri.h>
#include <infdim.h>
#include <cbgetdis.h>
C
      LOGICAL NOH2, IH8SM
C
      PARAMETER ( D1 = 1.0D0 , DM1 = -1.0D0 )
C
C ALLOCATE WORK SPACE
C
      KUFCAC = 1
      KFREE  = KUFCAC + N2ASHX
      IF (IREFSY .EQ. KSYMST) THEN
         KNEED  = KFREE + KZCONF
      ELSE
         KNEED  = KFREE
      END IF
      LFREE  = LWRK  - KFREE
      IF (KNEED .GT. LWRK) CALL ERRWRK('RSPCLI',KNEED,LWRK)
C
C DEFINE PACKING TO GET INTEGRAL DISTRIBUTION FROM H2AC
C AND OTHER PARAMETERS TO BE USED IN CISIGD
C
      DISTYP = 1
      IADH2  = -1
C     ... IADH2 .lt. 0: H2AC in core
      NOH2   = .FALSE.
      IH8SM  = .TRUE.
      ISPIN1 = 0
      ISPIN2 = 0
      DO 100 ISIM = 1,NCSIM
         CALL DSPTSI(NASHT,FCAC,WRK(KUFCAC))
C     ... unpack FCAC using CALL DSPTSI(N,ASP,ASI)
         CALL CISIGD(KSYMST,KSYMST,KZCONF,KZCONF,ZYCVEC(1,ISIM),
     *               EVECS(1,ISIM),WRK(KUFCAC),H2AC,NOH2,IH8SM,
     *               XINDX,ISPIN1,ISPIN2,WRK(KFREE),LFREE)
C        CALL RSPSIG(ICSYM,IHCSYM,NCDET,NHCDET,C,HC,FCAC,H2AC,IFLAG,
C    *                  NOH2,WORK,KFREE,LFREE)
         IF ( IPRRSP.GT.110 ) THEN
            WRITE(LUPRI,'(/A,A)')
     *      ' LINEAR TRANSFORMED VECTOR WITHOUT MCSCF ENERGY' ,
     *      ' PROJECTED OUT : Z COMPONENT'
            CALL OUTPUT(EVECS(1,ISIM),1,KZCONF,1,1,KZCONF,1,1,LUPRI)
         END IF
         CALL DAXPY(KZCONF,-EACTIV,ZYCVEC(1,ISIM),1,EVECS(1,ISIM),1)
         IF ( IPRRSP.GT.110 ) THEN
            WRITE(LUPRI,'(/A,A)')
     *      ' LINEAR TRANSFORMED VECTOR AFTER MCSCF ENERGY' ,
     *      ' PROJECTED OUT : ZCOMPONENT '
            CALL OUTPUT(EVECS(1,ISIM),1,KZCONF,1,1,KZCONF,1,1,LUPRI)
         END IF
         IF (IREFSY .EQ. KSYMST) THEN
C           ... remove CREF component of E[2] vector
            KDREF = KFREE
            CALL GETREF(WRK(KDREF),KZCONF)
Chj         ... KZCONF as we may need CREF in determinants
            T1 = DDOT(KZCONF,WRK(KDREF),1,EVECS(1,ISIM),1)
            CALL DAXPY(KZCONF,(-T1),WRK(KDREF),1,EVECS(1,ISIM),1)
         END IF
         IF (IPRRSP.GT.110) THEN
            WRITE(LUPRI,'(/A)')
     *        ' CSF PART OF E(2) LINEAR TRANSFORMED CSF VECTOR '
            CALL OUTPUT(EVECS(1,ISIM),1,KZCONF,1,1,
     *                  KZCONF,1,1,LUPRI)
         END IF
C        CALL DSPTSI(NASHT,FCAC,WRK(KUFCAC))
C        CALL RSPSIG(KSYMST,KSYMST,KZCONF,KZCONF,ZYCVEC(1+KZCONF,ISIM),
C    *               EVECS(1+KZVAR,ISIM),
C    *               WRK(KUFCAC),H2AC,IFLAG,
C    *               0,WRK(KWRK1),KFREE,LFREE)
C        IF ( IPRRSP.GT.110 ) THEN
C           WRITE(LUPRI,'(/A,A)')
C    *      ' LINEAR TRANSFORMED VECTOR WITHOUT MCSCF ENERGY' ,
C    *      ' PROJECTED OUT : Y COMPONENT'
C           CALL OUTPUT(EVECS(1+KZVAR,ISIM),
C    *                  1,KZCONF,1,1,KZCONF,1,1,LUPRI)
C        END IF
C         CALL DAXPY(KZCONF,-EACTIV,ZYCVEC(1+KZCONF,ISIM),1,
C    *              EVECS(1+KZVAR,ISIM),1)
C        IF ( IPRRSP.GT.110 ) THEN
C           WRITE(LUPRI,'(/A,A)')
C    *      ' LINEAR TRANSFORMED VECTOR AFTER MCSCF ENERGY' ,
C    *      ' PROJECTED OUT : Y COMPONENT '
C           CALL OUTPUT(EVECS(1+KZVAR,ISIM),
C    *                  1,KZCONF,1,1,KZCONF,1,1,LUPRI)
C        END IF
 100  CONTINUE
C
C END OF RSPCLI
C
      RETURN
      END
C---------------------------------------------------------
C  /* Deck rspctl */
      SUBROUTINE RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *                  LINEQ,LAB1,LAB2,GD,REDGD,REDE,REDS,
     *                  IBTYP,EIVAL,RESID,EIVEC,XINDX,WRK,LWRK)
C
C Written 6-Nov-1984 by Poul Joergensen
C Revisions:
C  7-Nov-1984 hjaaj
C
C PURPOSE:
C  LINEQ.EQ.FALSE  DIRECT THE SOLUTION THE GENERALIZED RSP
C                  EIGENVALUE PROBLEM
C
C                  ( E(2)-W(I)S(2))X(I) = 0
C
C  LINEQ.EQ.TRUE  DIRECT THE SOLUTION THE GENERALIZED RSP
C                  LINEAR EQUATIONS
C
C                  ( E(2)-W(I)S(2))X(I) - GD = 0
C
C THE PAIRED STRUCTURE OF THE EIGENVECTORS
C (Z(I),Y(I)) AND (Y(I),Z(I)) IS USED
C TO SET UP THE REDUCED GENERALIZED RSP EIGENVALUE
C AND LINEAR EQUATIONS
C
C TRIAL VECTORS ARE SPLIT SUCH THAT THEY ONLY CONTAIN EITHER
C NON-VANISHING ORBITAL OR CSF COMPONENTS
C
C KEXSTV: NUMBER OF START VECTORS
C KEXSIM: DESIRED NUMBER OF SIMULTANEOUS ROOTS IN A
C         MICROITERATION
C JEXSIM: THE NUMBER OF TRIAL VECTORS IN A MACROITERATION
C  ( TRIAL VECTORS ARE ONLY ADDED FOR THE NON CONVERGED VECTORS AND
C    LINEAR DEPENDENCE IS REMOVED)
C
C MAXRIT: MAXIMUM NUMBER OF MICROITERATIONS
C
#include <implicit.h>
#include <dummy.h>
C
      PARAMETER ( MAXVEC = 12, DM1 = -1.0D0 )
      DIMENSION CMO(*),UDV(NASHDI,*),PV(*),FC(*),FV(*)
      DIMENSION FCAC(*), H2AC(*), GD(*),REDGD(*)
      DIMENSION REDE(*),REDS(*),IBTYP(*),EIVAL(*),RESID(*),EIVEC(*)
      DIMENSION XINDX(*),WRK(*)
C
#include <ibndxdef.h>
C
C Used from common blocks:
C  /WRKRSP/: ??
C  /INFPRI/: LUERR
C
#include <priunit.h>
#include <inforb.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infdim.h>
#include <infpri.h>
#include <inftap.h>
#include <infsop.h>
C
      LOGICAL BINMEM, LINEQ, OPTSAV, LMAXIT
      CHARACTER*8 LAB123(3),RESTAR, LAB1, LAB2
      DATA LAB123/'********','********','RESPONPP'/
      DATA RESTAR/'RESTART '/
C
      CALL QENTER('RSPCTL')
      IF (IPRRSP .GE. 2) THEN
      IF (LINEQ) THEN
         WRITE (LUPRI,'(///2A//A,I3,A,L3/)')
     *   ' <<<  SOLVING SETS OF LINEAR EQUATIONS ',
     *   'FOR LINEAR RESPONSE PROPERTIES >>>',
     &   ' Operator symmetry =',KSYMOP,'; triplet = ',TRPLET
      ELSE
         WRITE (LUPRI,'(///2A//A,I3,A,L3/)')
     *   ' <<< EXCITATION ENERGIES',
     *   ' AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>',
     &   ' Operator symmetry =',KSYMOP,'; triplet = ',TRPLET
      END IF
      END IF
C
C  DEFINE MAXIMUM NUMBER OF TRIAL VECTORS IN RSPLIN AND RSPNEX
C
C  WORK SPACE REQUIREMENT IN RSPLIN
C
C  SPECIAL REQUIREMENT FOR AN ORBITAL TRIAL VECTOR
C
C  FROM RSPELI + TRIAL VECTORS
C
      KORB   =  KZYWOP + KZYVAR + 3*NORBT*NORBT + 2*NASHT*NORBT
C
C  FROM RSPOLI
C
      KORB   = KORB + 2*NORBT*NASHT + N2ASHX*NNASHX + 2*NORBT*NORBT
      IF (RSPCI) KORB = 0
C
C  SPECIAL REQUIREMENT FOR A CONFIGURATION TRIAL VECTOR
C
C  FROM RSPELI + TRIAL VECTORS
C
      KCSF   =  KZCONF + KZYVAR + NORBT*NORBT + 2*NORBT*NASHT
C
C  FROM RSPOLI
C
      IF (.NOT.RSPCI) KCSF = KCSF + N2ASHX + N2ASHX*N2ASHX
C
C  MAXIMUM SPACE FOR ONE TRIAL VECTOR
C
      KMAXVE = MAX(KORB,KCSF)
C
C  WORK SPACE REQUIREMENT FOR EACH RSPLIN CALL
C
      KUSE   = 2*NORBT*NORBT + 2*LBINTM
C
C  FROM RSPOLI
C
      KUSE   = KUSE + 2*N2ASHX
C                     2 :TO ASSURE SPACE FOR ALL PVCD MATRICES
C
C  WORK SPACE USED EACH TIME A CSF LINEAR TRANSFORMATION
C  OR A DENSITY MATRIX IS CONSTRUCTED
C
      KLIDEN = N2ASHX  + LACIMX
C  REFERENCE VECTOR ONLY NEEDED FOR MCSCF CALCULATION
      IF (.NOT.RSPCI) KLIDEN = KLIDEN + MAX(KZCONF,NCREF)
C
      KTOUSE = LWRK - KUSE - KLIDEN
C     SOPCLI reads 2p-2h diagonal (D(0) matrix)
      IF (SOPPA) KTOUSE = KTOUSE - KZCONF
      MAXSIM = MIN(MAXVEC,(KTOUSE-100)/KMAXVE)
      IF (IPRRSP.GT.37 .OR. MAXSIM.LE.0) THEN
         WRITE(LUPRI,*)
     *     ' KLIDEN,KUSE,KMAXVE,LBINTM,NORBT,NASHT,KZYVAR,LWRK '
         WRITE(LUPRI,*)
     *       KLIDEN,KUSE,KMAXVE,LBINTM,NORBT,NASHT,KZYVAR,LWRK
         WRITE(LUPRI,*)' LACIMX',LACIMX
      END IF
      IF (IPRRSP.GE.10) WRITE(LUPRI,'(/A,I8)')
     *     ' MAXIMUM NUMBER OF SIMULTANEOUS TRIAL VECTORS:',MAXSIM
C
      IF (MAXSIM.LE.0) THEN
         WRITE (LUPRI ,8000) LWRK,(100+KMAXVE-KTOUSE)
         WRITE (LUERR,8000) LWRK,(100+KMAXVE-KTOUSE)
         CALL QTRACE(LUERR)
         CALL QUIT('RSPCTL: INSUFFICIENT MEMORY ALLOCATION')
      ENDIF
 8000 FORMAT(/' RSPCTL, not enough memory for one trial vector',
     *       /'         LWRK =',I10,
     *       /'         need >',I10,' additional work memory space.')
C
C Space for RSPNEX:
C
      MAXNEX = 30
      MAXORP = 100
C
C REQUIREMENT FOR EACH TRIAL VECTOR
C
      KNEX   = KZYVAR + KZYWOP
C
C REQUIREMENT FOR EACH CALL TO RSPNEX
C
      KUSE   = KZYVAR + NORBT*NORBT + 2*N2ASHX
      IF (.NOT.RSPCI) KUSE = KUSE + MAX(NCREF,KZCONF) + LACIMX
      IF (OPTORB) THEN
         KOPT = 2*MAXORP*(MAXORP+1)/2 + 2*MAXORP + 3*KZYWOP*MAXORP
         KOPT = KOPT + KZYWOP
C
C  MAXIMUM FOR ORPLIN AND ORPSVE
C
         KTOT = KORB - KZYCON  -N2ASHX*NNASHX
         KOPT = KOPT + MAX(KZYVAR,KTOT)
      ELSE
         KOPT = 0
      END IF
C REQUIREMENT FOR EACH CALL TO RSPNEX
      IF ( KOPT .GT. KUSE ) KUSE = KOPT
      KLEFT  = LWRK - KUSE
C RSPNEX reads SOPPA 2p-2h diagonal of length KZCONF
      IF (SOPPA) KLEFT = KLEFT - KZCONF
      NSIM3 = MIN(MAXNEX,(KLEFT-100)/KNEX)
      IF(IPRRSP.GE.10 .OR. NSIM3 .LE. 0)  WRITE(LUPRI,'(/A,I8)')
     *' MAXIMUM NUMBER OF SIMULTANEOUS TRIAL VECTORS IN RSPNEX:',NSIM3
      IF (NSIM3.LE.0) THEN
         WRITE (LUPRI ,8000) LWRK,(100+KLEFT+KNEX)
         WRITE (LUERR,8000) LWRK,(100+KLEFT+KNEX)
         CALL QTRACE(LUERR)
         CALL QUIT('RSPCTL: INSUFFICIENT SPACE')
      ENDIF
      KVECS = 1
C
C CONSTRUCT E2 AND S2 EXPLICITLY USING UNIT TRIAL VECTORS
C
      IF (ABCHK ) THEN
         CALL RSPES2(GD,REDGD,REDE,REDS,EIVAL,EIVEC,
     *               CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *               XINDX,WRK,LWRK)
C        CALL RSPES2(GD,REDGD,REDE,REDS,EIVAL,EIVEC,
C    *               CMO,UDV,PV,FC,FV,FCAC,H2AC,
C    *               XINDX,WRK,LWRK)
         WRITE(LUPRI,'(/A)')' .ABCHK E2 AND S2 CONSTUCTED : STOP'
         CALL QUIT('RSPCTL: ABCHK FINISHED')
      END IF
C
C     Initialize micro iteration counter
C
      ITMIC = 0
C
C     SET UP INITIAL TRIAL VECTORS
C
      IF (IPRRSP.GT.15) THEN
         WRITE(LUPRI,*)' LINEQ,RESTPP,ORBSPC',LINEQ,RESTPP,ORBSPC
         IF (LINEQ) WRITE(LUPRI,*) ' FREQ =',(EIVAL(I),I=1,MAXSIM)
      END IF
C
      IF (LINEQ) THEN
         CALL LRST(BINMEM,MAXSIM,IBTYP,EIVAL,EIVEC,GD,
     *             REDGD,REDS,REDE,UDV,FCAC,H2AC,XINDX,
     *             LAB1,LAB2,WRK,LWRK)
         IF ((KZRED.GT.0) .AND. (KEXSTV.LE.0)) THEN
            ICTL = 2
            KBCVEC = 1
            KBOVEC = 1
            KEVECS = 1
            KWRK2  = 1
            LWRK2  = LWRK
            NSIM   = 0
            GO TO 1000
         END IF
         IF (ORBSPC) THEN
            DO 313 II = 1,KZWOPT
               GD(II) = GD(KZCONF+II)
 313        CONTINUE
C           CALL DCOPY(KZWOPT,GD(1+KZCONF),1,GD,1)
            DO 314 II = 1,KZWOPT
               GD(KZWOPT+II) = GD(KZVAR+KZCONF+II)
 314        CONTINUE
C           CALL DCOPY(KZWOPT,GD(1+KZVAR+KZCONF),1,
C    *                 GD(1+KZWOPT),1)
            CALL OROPTI(EIVAL,IBTYP,GD,CMO,UDV,PV,FC,FV,FCAC,XINDX,
     *                  WRK,LWRK)
            WRITE(LUPRI,'(/A)')
     * ' STOP: CHECK OF ORBITAL BLOCK OF E(2) AND S(2) FINISHED'
            CALL QUIT(' END OF ORBSPC CHECK')
         END IF
      ELSE
         IF (RESTPP) THEN
C           Restart in PP is only possible in first call
C           (in any later call of RSPCTL the restart vectors will be of
C            another symmetry).
            RESTPP = .FALSE.
            CALL  PPRST(IBTYP,REDS,REDE,UDV,XINDX,WRK,LWRK)
C           CALL  PPRST(IBTYP,REDS,REDE,UDV,XINDX,WRK,LWRK)
            IF (KZRED .LT. KEXSIM) THEN
               WRITE (LUPRI,'(//A,/A,I5,/A,I5)')
     *         ' RSPCTL: ".RESTPP" impossible, too few trial vectors.',
     *         '     number of trial vectors on file :',KZRED,
     *         '     number of roots (".ROOTS")     :',KEXSIM
               CALL QUIT('RESTART PP impossible,too few trial vectors.')
            END IF
            ICTL = 2
            KBCVEC = 1
            KBOVEC = 1
            KEVECS = 1
            KWRK2  = 1
            LWRK2  = LWRK
            NSIM   = 0
            GO TO 1000
         ELSE
            CALL PPST(IBTYP,XINDX,FCAC,H2AC,WRK,MAXSIM,BINMEM,LWRK)
         END IF
      ENDIF
C
      IF (IPRRSP.GT.5) WRITE(LUPRI,*)'AFTER START ROUTINE BINMEM',BINMEM
      JEXSIM = KEXSTV
C     ... JEXSIM is redefined in RSPNEX
C
C     --- START OF MICROITERATION LOOP ---
C
 100  CONTINUE
         CALL FLSHFO(LUPRI)
         ITMIC = ITMIC + 1
         IF (IPRRSP.GE.3)
     *   WRITE(LUPRI,'(/A,I5)') ' ** RSPCTL MICROITERATION NUMBER',ITMIC
C
C        CHECK IF ALL B VECTORS IN MEMORY
C        IF NOT POSITION LURSP3
C
 200     CONTINUE
         NSIM = MIN(MAXSIM,JEXSIM)
         IF (IPRRSP.GT.15) WRITE(LUPRI,'(A,3I5)')' NSIM,MAXSIM,JEXSIM',
     *      NSIM,MAXSIM,JEXSIM
         NCSIM = 0
         NOSIM = 0
         DO 105 I=1,NSIM
            IF (IBTYP(KZRED+I+KOFFTY).EQ.JBCNDX) THEN
               NCSIM = NCSIM + 1
            ELSE
               NOSIM = NOSIM +1
            ENDIF
 105     CONTINUE
         IF (IPRRSP.GE.15) WRITE(LUPRI,*)' NCSIM,NOSIM',NCSIM,NOSIM
         KBCVEC= KVECS
         KBOVEC= KBCVEC + KZCONF*NCSIM
         KWRK1 = KBOVEC + KZYWOP*NOSIM
         LWRK1 = LWRK   - KWRK1
         IF (LWRK1.LT.0) CALL ERRWRK('RSPCTL 1',KWRK1-1,LWRK)
         IF (.NOT.BINMEM) THEN
            CALL REWSPL(LURSP3)
#if !defined (VAR_SPLITFILES)
            IF (KOFFTY.EQ.1) READ (LURSP3)
#else
            JRSP3 = 0
            IF (KOFFTY.EQ.1)
     &           CALL READST(LURSP3,LRSP3,DUMMY,JRSP3,.FALSE.)
#endif
            DO 210 I=1,KZRED
#if defined (VAR_SPLITFILES)
 210           CALL READST(LURSP3,LRSP3,DUMMY,JRSP3,.FALSE.)
#else
 210           READ (LURSP3)
#endif
            ISTBC = KBCVEC
            ISTBO = KBOVEC
            DO 400 ISIM=1,NSIM
               IF (IBTYP(KZRED+ISIM+KOFFTY).EQ.JBCNDX) THEN
#if defined (VAR_SPLITFILES)
                  CALL READST(LURSP3,LRSP3,WRK(ISTBC),JRSP3,.TRUE.)
                  IF (LRSP3 .NE. KZCONF) CALL QUIT('Inconsistent '//
     &                 'record length on unit LURSP3 in RSPCTL')
#else
                  CALL READT(LURSP3,KZCONF,WRK(ISTBC))
#endif
                  ISTBC = ISTBC + KZCONF
               ELSE
#if defined (VAR_SPLITFILES)
                  CALL READST(LURSP3,LRSP3,WRK(ISTBO),JRSP3,.TRUE.)
                  IF (LRSP3 .NE. KZYWOP) CALL QUIT('Inconsistent '//
     &                 'record length on unit LURSP3 in RSPCTL')
#else
                  CALL READT(LURSP3,KZYWOP,WRK(ISTBO))
#endif
                  ISTBO = ISTBO + KZYWOP
               ENDIF
 400        CONTINUE
         END IF
C
         IF (IPRRSP.GT.110) THEN
            WRITE(LUPRI,'(/A)')' BEFORE RSPLIN'
            IF (NOSIM .GT. 0) THEN
               WRITE(LUPRI,'(/I5,A)')NOSIM,' ORBITAL TRIAL VECTORS'
               CALL OUTPUT(WRK(KBOVEC),1,KZYWOP,1,NOSIM,
     *                     KZYWOP,NOSIM,1,LUPRI)
            ENDIF
            IF (NCSIM .GT. 0) THEN
               WRITE(LUPRI,'(/I5,A)')NCSIM,' CONFIG. TRIAL VECTORS'
               CALL OUTPUT(WRK(KBCVEC),1,KZCONF,1,NCSIM,
     *                     KZCONF,NCSIM,1,LUPRI)
            ENDIF
         END IF
         CALL FLSHFO(LUPRI)

         CALL RSPELI(NCSIM,NOSIM,WRK(KBCVEC),WRK(KBOVEC),
     *               CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *               XINDX,WRK(KWRK1),LWRK1)
         CALL FLSHFO(LUPRI)
C
C        CALL RSPELI(NCSIM,NOSIM,ZYCVEC,ZYOVEC,
C    *               CMO,UDV,PV,FC,FV,FCAC,H2AC,
C    *               XINDX,WRK,LWRK)
C
          KEVECS = KWRK1
          KSOP1  = KEVECS + (NCSIM+NOSIM)*KZYVAR
          IF (SOPPA) THEN
             KWRK2 = KSOP1 + KZYWOP
          ELSE
             KWRK2 = KSOP1
          ENDIF
          LWRK2  = LWRK   - KWRK2
          IF (LWRK2.LT.0) CALL ERRWRK('RSPCTL 2',KWRK2-1,LWRK)
C
C        WRITE E(2)*X ON LURSP5
C
         CALL REWSPL(LURSP5)
         JRSP5 = 0
         IF (KOFFTY.EQ.1) THEN
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP5,LRSP5,DUMMY,JRSP5,.FALSE.)
#else
            READ (LURSP5)
#endif
         ENDIF
         DO 215 I=1,KZRED
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP5,LRSP5,DUMMY,JRSP5,.FALSE.)
#else
            READ (LURSP5)
#endif
 215     CONTINUE
         ISTEC = KEVECS
         ISTEO = ISTEC + NCSIM*KZYVAR
         DO 216 ISIM=1,NSIM
            IF (IBTYP(KZRED+ISIM+KOFFTY).EQ.JBCNDX) THEN
C
C   If SOPPA only write p-h part of transformed vector to disk
C   and construct 2p-2h part when necessary
C
               IF (SOPPA) THEN
                  CALL DCOPY(KZWOPT,WRK(ISTEC+KZCONF),1,
     *                       WRK(KSOP1),1)
                  CALL DCOPY(KZWOPT,WRK(ISTEC+KZVAR+KZCONF),1,
     *                       WRK(KSOP1+KZWOPT),1)
#if defined (VAR_SPLITFILES)
                  CALL WRITST(LURSP5,KZYWOP,WRK(KSOP1),JRSP5)
               ELSE
                  CALL WRITST(LURSP5,KZYVAR,WRK(ISTEC),JRSP5)
#else
                  CALL WRITT(LURSP5,KZYWOP,WRK(KSOP1))
               ELSE
                  CALL WRITT(LURSP5,KZYVAR,WRK(ISTEC))
#endif
               ENDIF
               ISTEC = ISTEC + KZYVAR
            ELSE
#if defined (VAR_SPLITFILES)
               CALL WRITST(LURSP5,KZYVAR,WRK(ISTEO),JRSP5)
#else
               CALL WRITT(LURSP5,KZYVAR,WRK(ISTEO))
#endif
               ISTEO = ISTEO + KZYVAR
            ENDIF
 216     CONTINUE
         KZRED =KZRED+NSIM
         KZYRED=2*KZRED
C
C PREPARE FOR RESTART
C
         IF (IPRRSP.GT.30) THEN
            WRITE(LUPRI,*)' BEFORE WRITING IBTYP ON LURSP1 LINEQ:',LINEQ
         END IF
         KTOT = KZRED + KOFFTY
         IF (IPRRSP.GT.7)
     *      WRITE(LUPRI,'(/A,I5)')' TRIAL VECTORS FOR RESTART:',KTOT
         LURSP1 = -1
         CALL GPOPEN(LURSP1,'RSPRST.E2C','UNKNOWN',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND (LURSP1)
         WRITE(LURSP1)LAB123,RESTAR
         WRITE(LURSP1)KSYMOP,KTOT,(IBTYP(I),I=1,KTOT)
C        flush any output buffers
         REWIND (LURSP1)
         CALL GPCLOSE(LURSP1,'KEEP')
C
         IF (IPRRSP.GT.110) THEN
            IF (NOSIM .GT. 0) THEN
               WRITE(LUPRI,'(/I5,A)')NOSIM,' ORBITAL TRIAL VECTORS'
               CALL OUTPUT(WRK(KBOVEC),1,KZYWOP,1,NOSIM,
     *                     KZYWOP,NOSIM,1,LUPRI)
            ENDIF
            IF (NCSIM .GT. 0) THEN
               WRITE(LUPRI,'(/I5,A)')NCSIM,' CONFIG. TRIAL VECTORS'
               CALL OUTPUT(WRK(KBCVEC),1,KZCONF,1,NCSIM,
     *                     KZCONF,NCSIM,1,LUPRI)
            ENDIF
            WRITE(LUPRI,'(/I5,A)')NSIM, ' E(2) LINEAR TRANSF. VECTORS'
            CALL OUTPUT(WRK(KEVECS),1,KZYVAR,1,NSIM,
     *                     KZYVAR,NSIM,1,LUPRI)
         ENDIF
         JEXSIM =JEXSIM-NSIM
         IF (JEXSIM.GT.0) THEN
C
C           CALL RSPRED(1,..) INCREASE DIMENSION OF REDUCED RSP EQUATION
C
            CALL RSPRED(1,LINEQ,NSIM,IBTYP,GD,REDGD,REDE,REDS,
     *                  EIVAL,EIVEC,WRK(KBCVEC),WRK(KBOVEC),
     *                  UDV,WRK(KEVECS),XINDX,WRK(KWRK2),LWRK2)
C           CALL RSPRED (ICTL,LINEQ,N,IBTYP,GD,REDGD,REDE,REDS,
C    *                 EIVAL,EIVEC,BCVEC,BOVEC,UDV,EVECS,
C    *                 XINDX,WRK,LWRK)
            GO TO 200
C       ^-----------
         END IF
         ICTL = 3
 1000    CONTINUE
C
C        CALL RSPRED(,3,..) INCREASE DIMENSION OF REDUCED RSP EQUATION
C        AND SOLVE FOR EIGENVALUES AND EIGENVECTORS
C
         CALL RSPRED(ICTL,LINEQ,NSIM,IBTYP,GD,REDGD,REDE,REDS,
     *               EIVAL,EIVEC,WRK(KBCVEC),WRK(KBOVEC),
     *               UDV,WRK(KEVECS),XINDX,WRK(KWRK2),LWRK2)
C
C        CALCULATE TRANSITION MOMENT
C
C        CALL RSPMOM(IBTYP,EIVAL,EIVEC,BVECS,SVECS,PRVEC,NSIM)
C
C SET UP REDUCED MATRIX EXPLICLY AND CHECK BLOCK STRUCTURE IN E2 AND S2
C
         IF (ABSYM) THEN
            CALL E2SYM(IBTYP,CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,
     &                 WRK,LWRK)
         END IF
C        CREATE KEXSIM NEW LINEAR INDEPENDENT TRIAL VECTORS IN NEXT()
C        FROM THE EIGENVECTORS OF THE REDUCED  RSP EQUATION
C
         IF (KEXSIM.LE.MAXSIM) THEN
            BINMEM=.TRUE.
         ELSE
            BINMEM=.FALSE.
         ENDIF
         OPTSAV = OPTORB
         ITCIN = 0
         DO 510 IIB = 1,(KOFFTY+KZRED)
            IF (IBTYP(IIB).EQ.JBCNDX) ITCIN = ITCIN + 1
 510     CONTINUE
         ITCIN = ITCIN - KOFFTY
         IF (.NOT.LINEQ .AND. ITCIN.LE.2) OPTORB = .FALSE.
C        880226, PJ+HJAAJ: big problems with optorb for
C        RSPPP if no previous orbital trial vectors.
         LMAXIT = (ITMIC .GE. MAXIT)
         CALL RSPNEX(LINEQ,LMAXIT,NSIM3,IBTYP,EIVAL,RESID,EIVEC,GD,
     *               CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
C        CALL RSPNEX(LINEQ,LMAXIT,NSIM,IBTYP,EIVAL,EIVEC,GD,
C    *               CMO,UDV,PVX,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
         OPTORB = OPTSAV
         IF (KCONV.EQ.-2) THEN
C           (MAXIMUM DIMENSION OF REDUCED SPACE EXCEEDED)
            NWARN = NWARN + 1
            WRITE (LUERR,'(/A/A)')
     *          ' *** RSPCTL WARNING-MICROITERATIONS STOPPED BECAUSE',
     *          '     MAXIMUM DIMENSION OF REDUCED SPACE EXCEEDED.'
            WRITE (LUPRI,'(/A/A)')
     *          ' *** RSPCTL WARNING-MICROITERATIONS STOPPED BECAUSE',
     *          '     MAXIMUM DIMENSION OF REDUCED SPACE EXCEEDED.'
         ELSE IF (KCONV.LT.0) THEN
C           (LINEAR DEPENDENCE BETWEEN NEW TRIAL VECTOR )
            NWARN = NWARN + 1
            WRITE (LUERR,'(/A/A)')
     *          ' *** RSPCTL WARNING-MICROITERATIONS STOPPED BECAUSE',
     *          '     OF LINEAR DEPENDENCE BETWEEN NEW TRIAL VECTORS'
            WRITE (LUPRI,'(/A/A)')
     *          ' *** RSPCTL WARNING-MICROITERATIONS STOPPED BECAUSE',
     *          '     OF LINEAR DEPENDENCE BETWEEN NEW TRIAL VECTORS'
         ELSE IF(KCONV.GT.0)THEN
C           (CONVERGED)
            IF (IPRRSP .GE. 0)
     &      WRITE(LUPRI,'(/A)')' *** RSPCTL MICROITERATIONS CONVERGED'
         ELSE
C           (NOT CONVERGED)
            IF (ITMIC.GE.MAXIT) THEN
C              (MAX NO OF MICROITERATIONS REACHED)
               NWARN = NWARN + 1
               WRITE(LUERR,'(/A,I4,A)')
     *         ' *** RSPCTL WARNING-MAXIMUM NUMBER OF MICROITERATIONS,',
     *         ITMIC,', REACHED'
               WRITE(LUPRI,'(/A,I4,A)')
     *         ' *** RSPCTL WARNING-MAXIMUM NUMBER OF MICROITERATIONS,',
     *         ITMIC,', REACHED'
            ELSE
               GO TO 100
C     ^-----------------
            END IF
         END IF
C
C     --- END OF MICROITERATION LOOP ---
C
C
C     END OF RSPCTL
C
      CALL QEXIT('RSPCTL')
      RETURN
      END
C  /* Deck rspdia */
      SUBROUTINE RSPDIA(A,SHIFT)
C
C CONSTRUCT DIAGONAL OF E(2) + SHIFT* S(2)
C
#include <implicit.h>
      DIMENSION A(*)
      PARAMETER ( DMP5 = -0.5D0, D1 = 1.0D0 , DM2=-2.0D0 )
      PARAMETER ( D0 = 0.0D0 )
C
#include <infrsp.h>
#include <wrkrsp.h>
#include <inftap.h>
#include <infpri.h>
C
      CALL REWSPL(LURSP4)
#if defined (VAR_SPLITFILES)
      JRSP4 = 0
#endif
C
C *** now READ diagonals of E[2]-matrix needed for RSPLIN
C
      IF (KZCONF .GT. 0) THEN
#if defined (VAR_SPLITFILES)
         CALL READST(LURSP4,LRSP4,A,JRSP4,.TRUE.)
         IF (LRSP4 .NE. KZCONF) CALL QUIT('Inconsistent record length'//
     &        ' on unit LURSP4 in RSPDIA')
#else
         CALL READT (LURSP4,KZCONF,A)
#endif
      END IF
      IF (KZWOPT .GT. 0) THEN
#if defined (VAR_SPLITFILES)
         CALL READST(LURSP4,LRSP4,A(1+KZCONF),JRSP4,.TRUE.)
         IF (LRSP4 .NE. KZWOPT) CALL QUIT('Inconsistent record length'//
     &        ' on unit LURSP4 in RSPDIA (2)')
#else
         CALL READT (LURSP4,KZWOPT,A(1+KZCONF))
#endif
      END IF
C
C A(1:KZVAR):         (E+SHIFT*S)
C A(KZVAR+1:KZYVAR):  (E-SHIFT*S)
C
C
      IF (SHIFT .NE. D0) THEN
         IF (KZCONF .GT. 0) THEN
            DO 100 I = 1,KZCONF
Chj            A(I+KZVAR) = A(I+KZVAR)*SHIFT
               A(I+KZVAR) = SHIFT
 100        CONTINUE
         END IF
         IF (KZWOPT .GT. 0) THEN
C *** now READ orbital diagonal of S[2]-matrix
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP4,LRSP4,A(1+KZVAR+KZCONF),JRSP4,.TRUE.)
            IF (LRSP4 .NE. KZWOPT) CALL QUIT('Inconsistent record'//
     &           ' length on unit LURSP4 in RSPDIA (3)')
#else
            CALL READT (LURSP4,KZWOPT,A(1+KZVAR+KZCONF))
#endif
            CALL DSCAL(KZWOPT,SHIFT,A(1+KZVAR+KZCONF),1)
         END IF
         CALL DAXPY(KZVAR,D1,A(KZVAR+1),1,A,1)
         CALL DAXPY(KZVAR,DMP5,A,1,A(KZVAR+1),1)
         CALL DSCAL(KZVAR,DM2,A(KZVAR+1),1)
      ENDIF
C
C *** END OF RSPDIA
C
      RETURN
      END
C  /* Deck rspedg */
      SUBROUTINE RSPEDG(A)
C
C CONSTRUCT DIAGONAL OF E(2) (z part)
C
#include <implicit.h>
      DIMENSION A(*)
C
#include <infrsp.h>
#include <wrkrsp.h>
#include <inftap.h>
#include <infpri.h>
C
      CALL REWSPL(LURSP4)
#if defined (VAR_SPLITFILES)
      JRSP4 = 0
#endif
C
C *** now READ diagonals of E[2]-matrix needed for RSPLIN
C
      IF (KZCONF .GT. 0) THEN
#if defined (VAR_SPLITFILES)
         CALL READST(LURSP4,LRSP4,A,JRSP4,.TRUE.)
         IF (LRSP4 .NE. KZCONF) THEN
            WRITE(LUPRI,*) 'ERROR, LRSP4 .ne. KZCONF',LRSP4,KZCONF
            CALL QUIT('Inconsistent record length'//
     &        ' on unit LURSP4 in RSPEDG - conf')
         END IF
#else
         CALL READT (LURSP4,KZCONF,A)
#endif
      END IF
      IF (KZWOPT .GT. 0) THEN
#if defined (VAR_SPLITFILES)
         CALL READST(LURSP4,LRSP4,A(1+KZCONF),JRSP4,.TRUE.)
         IF (LRSP4 .NE. KZWOPT) THEN
            WRITE(LUPRI,*) 'ERROR, LRSP4 .ne. KZWOPT',LRSP4,KZWOPT
            CALL QUIT('Inconsistent record length'//
     &        ' on unit LURSP4 in RSPEDG - orb')
         END IF
#else
         CALL READT (LURSP4,KZWOPT,A(1+KZCONF))
#endif
      END IF
C
C *** END OF RSPEDG
C
      RETURN
      END
C  /* Deck rspsod */
      SUBROUTINE RSPSOD(A)
C
C READ DIAGONAL OF ORBITAL PART OF S(2)
C
#include <implicit.h>
#include <dummy.h>
      DIMENSION A(*)
C
#include <infrsp.h>
#include <wrkrsp.h>
#include <inftap.h>
#include <infpri.h>
C
      CALL REWSPL(LURSP4)
      JRSP4 = 0
C
C *** now SKIP diagonals of E[2]-matrix needed for RSPLIN
C
      IF (KZCONF .GT. 0) THEN
#if defined (VAR_SPLITFILES)
         CALL READST(LURSP4,LRSP4,DUMMY,JRSP4,.FALSE.)
#else
         READ(LURSP4)
#endif
      END IF
      IF (KZWOPT .GT. 0) THEN
#if defined (VAR_SPLITFILES)
         CALL READST(LURSP4,LRSP4,DUMMY,JRSP4,.FALSE.)
#else
         READ(LURSP4)
#endif
      END IF
C
C *** now READ orbital diagonal of S[2]-matrix needed for RSPTRN
C
      IF (KZWOPT .GT. 0) THEN
#if defined (VAR_SPLITFILES)
         CALL READST(LURSP4,LRSP4,A(1),JRSP4,.TRUE.)
         IF (LRSP4 .NE. KZWOPT) CALL QUIT('Inconsistent record length'//
     &        ' on unit LURSP4 in RSPSOD')
#else
         CALL READT (LURSP4,KZWOPT,A(1))
#endif
      END IF
C
C *** END OF RSPSOD
C
      RETURN
      END
C  /* Deck rsppp */
      SUBROUTINE RSPPP(CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
C
C  Purpose:
C     CONTROL CALCULATION OF EXCITATION ENERGIES AND TRANSITION
C     MOMENTS
C     (PP for Polarization Propagator)
C
#include <implicit.h>
#include "maxorb.h"
#include <dummy.h>
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
      CHARACTER BLANK*8, PRI_LINE*80
C
      PARAMETER ( MAXSIM = 15, LMXFR = 1000, BLANK = '        ' )
      PARAMETER ( TOLRTT = 1.0D-12)
      DIMENSION OSCSTR(MAXSIM,12),TRSMOM(MAXSIM,60)
      DIMENSION RMLEN(3,3,LMXFR), RMVEL(3,3,LMXFR)
      DIMENSION RQLEN(3,3,LMXFR), RQVEL(3,3,LMXFR)
      DIMENSION RTLEN(3,3,LMXFR), RTVEL(3,3,LMXFR)
      DIMENSION XLEVICI(3,3,3)
      DIMENSION XNACS2(MAXSIM+1) ,XNACGD(MAXSIM+1)
      LOGICAL REFST
      LOGICAL DOMO2
      PARAMETER ( D2 = 2.0D0, D100 = 100.0D0, D0 = 0.0D0)
      PARAMETER ( D2R3 = (2.0D0/3.0D0))
C
C
#include <codata.h>
      PARAMETER (RAU2CGS = ECHARGE*ECHARGE*XTANG*CCM*1D36*HBAR/EMASS)
      PARAMETER (ESUECD = RAU2CGS)
C
C Used from common blocks:
C  /INFOPT/ : EMCSCF
C  /INFRSP/ : most items (/INFRSP/ gives control information for
C                         the response calculation(s) )
C  /WRKRSP/ :
C
#include <priunit.h>
#include <infopt.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <infpp.h>
#include <inflr.h>
#include <inforb.h>
#include <infdim.h>
#include "infinp.h"
#include <infpri.h>
#include <inftap.h>
#include <esg.h>
C allocate work space for (oriented) electronic CD and
C for reduced E(2) and reduced S(2)
!     DOMO2 = DODFT .OR. DOHFSRDFT .OR. DOMCSRDFT
      DOMO2 = DODFT
      IF (RSPOCD) THEN
         IF (LMXFR .LT. NPPCNV(KSYMOP)) THEN
            NINFO = NINFO + 1
            WRITE(LUPRI,*)
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) '***** NOTICE:'
            WRITE(LUPRI,*) 'RSPPP: Rotatory strength tensors',
     &                     ' will only be calculated for the'
            WRITE(LUPRI,*) '       lowest ',LMXFR,' states.'
            WRITE(LUPRI,*)
            WRITE(LUPRI,*)
         END IF
         LFREX  = NPPCNV(KSYMOP)
         LTRVEC = 3*NPPCNV(KSYMOP)
         LTRTEN = 3*3*NPPCNV(KSYMOP)
         LROTL  = NPPCNV(KSYMOP)
         LROTV  = NPPCNV(KSYMOP)
         CALL DZERO(XLEVICI,3*3*3)
         XLEVICI(2,3,1) =  1.0D0
         XLEVICI(3,2,1) = -1.0D0
         XLEVICI(3,1,2) =  1.0D0
         XLEVICI(1,3,2) = -1.0D0
         XLEVICI(1,2,3) =  1.0D0
         XLEVICI(2,1,3) = -1.0D0
      ELSE IF (RSPECD) THEN
         LFREX  = NPPCNV(KSYMOP)
         LTRVEC = 3*NPPCNV(KSYMOP)
         LTRTEN = 0
         LROTL  = 0
         LROTV  = 0
      ELSE
         LFREX  = 0
         LTRVEC = 0
         LTRTEN = 0
         LROTL  = 0
         LROTV  = 0
      END IF
      IF (DOMO2) THEN
         KCMO2 = 1
         KLAST = KCMO2 + 2*N2ORBX
         LWRK1 = LWRK + 1 - KLAST
         CALL DFTMOMO(WRK(KCMO2),WRK(KLAST),LWRK1,IPRPP)
      ELSE
         KLAST = 1
      END IF
      KFREX  = KLAST
      KTRLEN = KFREX  + LFREX
      KTRVEL = KTRLEN + LTRVEC
      KTRMAG = KTRVEL + LTRVEC
      KTQLEN = KTRMAG + LTRVEC
      KTQVEL = KTQLEN + LTRTEN
      KROTL  = KTQVEL + LTRTEN
      KROTV  = KROTL  + LROTL
      KREDE  = KROTV  + LROTV
      KREDS  = KREDE  + MAXRM*MAXRM
      KIBTYP = KREDS  + MAXRM*MAXRM
      KEIVAL = KIBTYP + MAXRM
      KRESID = KEIVAL + MAXRM
      KEIVEC = KRESID + MAXRM
      KWRK1  = KEIVEC + MAXRM*MAXRM
      LWRK1  = LWRK + 1 - KWRK1
      IF (IPRPP .GT. 2 .OR. LWRK1 .LT. 2*KZYVAR) THEN
         WRITE(LUPRI,'(/A)') ' >>> IN RSPPP:'
         WRITE(LUPRI,*)' THCPP,MAXRM ',THCPP,MAXRM
         WRITE(LUPRI,*)' KSYMOP,NGPPP(KSYMOP) ',KSYMOP,NGPPP(KSYMOP)
         WRITE(LUPRI,*)' LWRK ,LWRK1 ',LWRK,LWRK1
      END IF

      CALL DZERO(WRK(KFREX),LFREX)
      CALL DZERO(WRK(KTRLEN),LTRVEC)
      CALL DZERO(WRK(KTRVEL),LTRVEC)
      CALL DZERO(WRK(KTRMAG),LTRVEC)
      CALL DZERO(WRK(KTQLEN),LTRTEN)
      CALL DZERO(WRK(KTQVEL),LTRTEN)
      CALL DZERO(WRK(KROTL),LROTL)
      CALL DZERO(WRK(KROTV),LROTV)
C
C ALLOCATE WORK SPACE FOR EIGENVECTORS AND TRANSITION MOMENTS
C
C 941003-hjaaj
C Requirement for 'HDO' and REFST true
C (Was previously always allocated, also when not needed)
C
      REFST = .FALSE.
      IF (KZCONF .GT. 0 .AND. KOFFTY.EQ.1) THEN
         DO 100 IOP = 1,NGPPP(KSYMOP)
  100       IF (LBLPP(KSYMOP,IOP)(1:3) .EQ. 'HDO') REFST = .TRUE.
      END IF
      LNEED = 100
      IF (REFST) LNEED = LNEED + KZCONF
C
C REQUIREMENT in GETGPV or RSPEVE (RSPEVE needs less than GETGPV)
C
      IF (NGPPP(KSYMOP) .GT. 0) THEN
         LNEED = LNEED + KZYVAR + N2ORBX + N2BASX
         IF (SOPPA) THEN
            LNEED = LNEED + 2*N2ORBX + 2*NISHT*NSSHT
         ELSE
            LNEED = LNEED + N2ASHX + LACIMX
         END IF
      ELSE
         LNEED = LNEED + MAX(KZCONF,KZYWOP)
      END IF
C
C MAXIMUM NUMBER OF SIMULTANEOUS SOLUTION VECTORS
C
      NSIM = MIN(KEXCNV, MAXSIM, (LWRK1-LNEED)/KZYVAR )
      IF (IPRPP .GT. 2 .OR. NSIM .LE. 0) THEN
         LWRK2 = KWRK1 + LNEED + KZYVAR
C        ... need at least space for one KZYVAR (NSIM = 1)
         WRITE (LUPRI,*) ' KEXCNV,NSIM,LWRK2 ',KEXCNV,NSIM,LWRK2
         IF (REFST) WRITE (LUPRI,*) ' HDO and REFST ',REFST
         IF (NSIM.LE.0) CALL ERRWRK('RSPPP work space',-LWRK2,LWRK)
      END IF
C
      KBVECS = KWRK1
      KCREF  = KBVECS + NSIM*KZYVAR
      IF (REFST) THEN
         KWRK2  = KCREF  + KZCONF
         CALL GETREF(WRK(KCREF),KZCONF)
      ELSE
         KWRK2  = KCREF
      END IF
      LWRK2  = LWRK   - KWRK2
C
C
C     940926-hjaaj: old lwrk1 check removed
C                   (3*kzyvar is much more than needed),
C                   we now let RSPCTL check if sufficient memory.
C
      KZRED  = 0
      KZYRED = 0
      THCRSP = THCPP
      IPRRSP = IPRPP
      MAXIT  = MAXITP
C
C     Call RSPCTL to solve propagator eigen problem
C
      CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *            .FALSE.,BLANK,BLANK,DUMMY,DUMMY,WRK(KREDE),WRK(KREDS),
     *            WRK(KIBTYP),WRK(KEIVAL),WRK(KRESID),WRK(KEIVEC),
     *            XINDX,WRK(KWRK1),LWRK1)
C     CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
C    *            LINEQ,GD,REDGD,REDE,REDS,
C    *            IBTYP,EIVAL,EIVEC,XINDX,WRK,LWRK)
C
C   Save the excitation energies for excited state gradient 
C   (singlet only for the moment)
C
      ISPN = 1
      DO I = 1,KEXCNV
         EXCITA(KSYMOP,I,ISPN) = WRK(KEIVAL-1+I)
      END DO
C
C CALCULATE EIGENVECTORS AND TRANSITION MOMENTS
C
      DO 500 ISIM = 1,KEXCNV,NSIM
         NBX = MIN( NSIM,(KEXCNV+1-ISIM) )
         CALL RSPEVE(WRK(KIBTYP),WRK(KEIVAL),WRK(KEIVEC),WRK(KBVECS),
     *               WRK(KWRK2),NBX,(ISIM-1))
C        CALL RSPEVE(IBTYP,EIVAL,EIVEC,BVECS,WRK,NBX,IBOFF)
         IF ( NGPPP(KSYMOP) .GT. 0 ) THEN
            JOP = 1
            DO 600 IOP = 1,NGPPP(KSYMOP)
               CALL GETGPV(LBLPP(KSYMOP,IOP),FC,FV,CMO,UDV,PV,
     *                     XINDX,ANTSYM,WRK(KWRK2),LWRK2 )
               IF (LBLPP(KSYMOP,IOP)(1:3) .EQ. 'HDO') THEN
                  IF (REFST) THEN
                     XNACGD(NBX+1) = DDOT(KZCONF,WRK(KWRK2),1,
     *                               WRK(KCREF),1)
                  ENDIF
                  DO 450 INUM = 1,NBX
                     XNACGD(INUM) = DDOT(KZYVAR,WRK(KWRK2),1,
     *                              WRK(KBVECS+(INUM-1)*KZYVAR),1)
 450              CONTINUE
                  CALL RDS2(LBLPP(KSYMOP,IOP),UDV,XINDX,
     *                      WRK(KWRK2),LWRK2)
                  DO 455 INUM = 1,NBX
                     XNACS2(INUM) = DDOT(KZYVAR,WRK(KWRK2),1,
     *                              WRK(KBVECS+(INUM-1)*KZYVAR),1)
 455              CONTINUE
                  IF (REFST) THEN
                     XNACS2(NBX+1) = DDOT(KZCONF,WRK(KWRK2),1,
     *                               WRK(KCREF),1)
                  ENDIF
                  WRITE(LUPRI,'(/A,A,2(/A))')
     &               ' OPERATOR TYPE:    ',LBLPP(KSYMOP,IOP),
     &               ' NON ADIABATIC COUPLING ELEMENTS',
     &             ' STATE  Excit.en.    HDO term   RDS2 term     total'
                  IF (REFST) THEN
                     WRITE(LUPRI,'(/I5,F12.8,3F12.6)')
     *               0,0.D0,XNACGD(NBX+1),XNACS2(NBX+1),XNACGD(NBX+1)
     *               +XNACS2(NBX+1)
                  ENDIF
                  DO 456 INUM = 1,NBX
                     WRITE(LUPRI,'(/I5,F12.8,3F12.6)')
     *               (ISIM-1+INUM),WRK(KEIVAL-1+(ISIM-1+INUM)),
     *               XNACGD(INUM),XNACS2(INUM),XNACGD(INUM)+XNACS2(INUM)
 456              CONTINUE
                  GO TO 600
               END IF
               IF (IPRRSP.GT.0) WRITE(LUPRI,'(/A,A)')
     *            '@ Transition operator type:    ',LBLPP(KSYMOP,IOP)
               IF (IPRRSP.GT.11) WRITE(LUPRI,'(A,I4,A,L4)')
     *              ' Operator symmetry:',KSYMOP,'; triplet = ',TRPLET
               DO 550 INUM = 1,NBX
                 LOCST = ISIM-1+INUM
                  IF ((RSPECD.OR.RSPOCD) .AND. (LOCST.GT.NPPCNV(KSYMOP))
     &               ) THEN
                     CALL QUIT('State count error in RSPPP')
                  END IF
                  XMOM = DDOT(KZYVAR,WRK(KWRK2),1,
     *                        WRK(KBVECS+(INUM-1)*KZYVAR),1)
                  IF (IPRRSP.GT.0) WRITE(LUPRI,751)
     &               (ISIM-1+INUM),XMOM, WRK(KEIVAL-1+ISIM-1+INUM)*XTEV
 751              FORMAT('@ STATE NO:',I5,' *TRANSITION MOMENT:'
     *                 ,1P,G16.8,' *ENERGY(eV):',1P,G16.8)
                  IF (RSPOCD .OR. RSPECD) THEN
                     WRK(KFREX+LOCST-1) = WRK(KEIVAL-1+ISIM-1+INUM)
                  END IF
                  IF ( LBLPP(KSYMOP,IOP) .EQ. 'XDIPLEN ' .OR.
     *                 LBLPP(KSYMOP,IOP) .EQ. 'YDIPLEN ' .OR.
     *                 LBLPP(KSYMOP,IOP) .EQ. 'ZDIPLEN ' ) THEN
                      OSCIL = D2R3*XMOM*XMOM*WRK(KEIVAL-1+ISIM-1+INUM)
                      IF (IPRRSP.GT.11) WRITE(LUPRI,'(/A,1P,G16.8)')
     *                   ' OSCILLATOR STRENGTH (LENGTH)   :',OSCIL
                      TRSMOM(INUM,JOP) = XMOM
                      OSCSTR(INUM,JOP) = OSCIL
                      IF (RSPECD .OR. RSPOCD) THEN
                         IF (LBLPP(KSYMOP,IOP).EQ.'XDIPLEN ') THEN
                            KOFF = KTRLEN + 3*(LOCST - 1)
                            WRK(KOFF) = XMOM
                         ELSE IF (LBLPP(KSYMOP,IOP).EQ.'YDIPLEN ') THEN
                            KOFF = KTRLEN + 3*(LOCST - 1) + 1
                            WRK(KOFF) = XMOM
                         ELSE IF (LBLPP(KSYMOP,IOP).EQ.'ZDIPLEN ') THEN
                            KOFF = KTRLEN + 3*(LOCST - 1) + 2
                            WRK(KOFF) = XMOM
                         END IF
                      END IF
                  ELSE IF ( LBLPP(KSYMOP,IOP) .EQ. 'XDIPVEL ' .OR.
     *                 LBLPP(KSYMOP,IOP) .EQ. 'YDIPVEL ' .OR.
     *                 LBLPP(KSYMOP,IOP) .EQ. 'ZDIPVEL ' ) THEN
                      OSCIL = D2R3*XMOM*XMOM/WRK(KEIVAL-1+ISIM-1+INUM)
                      IF (IPRRSP.GT.11) WRITE(LUPRI,'(/A,1P,G16.8)')
     *                   ' OSCILLATOR STRENGTH (VELOCITY) :',OSCIL
                      TRSMOM(INUM,JOP) = XMOM
                      OSCSTR(INUM,JOP) = OSCIL
                      IF (RSPECD .OR. RSPOCD) THEN
                         IF (LBLPP(KSYMOP,IOP).EQ.'XDIPVEL ') THEN
                            KOFF = KTRVEL + 3*(LOCST - 1)
                            WRK(KOFF) = XMOM
                         ELSE IF (LBLPP(KSYMOP,IOP).EQ.'YDIPVEL ') THEN
                            KOFF = KTRVEL + 3*(LOCST - 1) + 1
                            WRK(KOFF) = XMOM
                         ELSE IF (LBLPP(KSYMOP,IOP).EQ.'ZDIPVEL ') THEN
                            KOFF = KTRVEL + 3*(LOCST - 1) + 2
                            WRK(KOFF) = XMOM
                         END IF
                      END IF
                  ELSE IF ( LBLPP(KSYMOP,IOP) .EQ. 'XANGMOM ' .OR.
     *                 LBLPP(KSYMOP,IOP) .EQ. 'YANGMOM ' .OR.
     *                 LBLPP(KSYMOP,IOP) .EQ. 'ZANGMOM ' ) THEN
                      TRSMOM(INUM,JOP) = XMOM
                      IF (RSPECD .OR. RSPOCD) THEN
                         IF (LBLPP(KSYMOP,IOP).EQ.'XANGMOM ') THEN
                            KOFF = KTRMAG + 3*(LOCST - 1)
                            WRK(KOFF) = XMOM
                         ELSE IF (LBLPP(KSYMOP,IOP).EQ.'YANGMOM ') THEN
                            KOFF = KTRMAG + 3*(LOCST - 1) + 1
                            WRK(KOFF) = XMOM
                         ELSE IF (LBLPP(KSYMOP,IOP).EQ.'ZANGMOM ') THEN
                            KOFF = KTRMAG + 3*(LOCST - 1) + 2
                            WRK(KOFF) = XMOM
                         END IF
                      END IF
                  ELSE IF ( LBLPP(KSYMOP,IOP) .EQ. 'XXSECMOM' .OR.
     *                 LBLPP(KSYMOP,IOP) .EQ. 'XYSECMOM' .OR.
     *                 LBLPP(KSYMOP,IOP) .EQ. 'XZSECMOM' .OR.
     *                 LBLPP(KSYMOP,IOP) .EQ. 'YYSECMOM' .OR.
     *                 LBLPP(KSYMOP,IOP) .EQ. 'YZSECMOM' .OR.
     *                 LBLPP(KSYMOP,IOP) .EQ. 'ZZSECMOM' ) THEN
                      TRSMOM(INUM,JOP) = XMOM
                      IF (RSPECD .OR. RSPOCD) THEN
                         IF (LBLPP(KSYMOP,IOP).EQ.'XXSECMOM') THEN
                            I1 = 1
                            I2 = 1
                            KOFF = KTQLEN + 3*3*(LOCST - 1) + 3*(I2 - 1)
     &                           + I1 - 1
                            WRK(KOFF) = XMOM
                         ELSE IF (LBLPP(KSYMOP,IOP).EQ.'XYSECMOM') THEN
                            I1 = 1
                            I2 = 2
                            KOFF1 = KTQLEN + 3*3*(LOCST - 1)
     &                            + 3*(I2 - 1) + I1 - 1
                            KOFF2 = KTQLEN + 3*3*(LOCST - 1)
     &                            + 3*(I1 - 1) + I2 - 1
                            WRK(KOFF1) = XMOM
                            WRK(KOFF2) = XMOM
                         ELSE IF (LBLPP(KSYMOP,IOP).EQ.'XZSECMOM') THEN
                            I1 = 1
                            I2 = 3
                            KOFF1 = KTQLEN + 3*3*(LOCST - 1)
     &                            + 3*(I2 - 1) + I1 - 1
                            KOFF2 = KTQLEN + 3*3*(LOCST - 1)
     &                            + 3*(I1 - 1) + I2 - 1
                            WRK(KOFF1) = XMOM
                            WRK(KOFF2) = XMOM
                         ELSE IF (LBLPP(KSYMOP,IOP).EQ.'YYSECMOM') THEN
                            I1 = 2
                            I2 = 2
                            KOFF = KTQLEN + 3*3*(LOCST - 1) + 3*(I2 - 1)
     &                           + I1 - 1
                            WRK(KOFF) = XMOM
                         ELSE IF (LBLPP(KSYMOP,IOP).EQ.'YZSECMOM') THEN
                            I1 = 2
                            I2 = 3
                            KOFF1 = KTQLEN + 3*3*(LOCST - 1)
     &                            + 3*(I2 - 1) + I1 - 1
                            KOFF2 = KTQLEN + 3*3*(LOCST - 1)
     &                            + 3*(I1 - 1) + I2 - 1
                            WRK(KOFF1) = XMOM
                            WRK(KOFF2) = XMOM
                         ELSE IF (LBLPP(KSYMOP,IOP).EQ.'ZZSECMOM') THEN
                            I1 = 3
                            I2 = 3
                            KOFF = KTQLEN + 3*3*(LOCST - 1) + 3*(I2 - 1)
     &                           + I1 - 1
                            WRK(KOFF) = XMOM
                         END IF
                      END IF
                  ELSE IF ( LBLPP(KSYMOP,IOP) .EQ. 'XXROTSTR' .OR.
     *                 LBLPP(KSYMOP,IOP) .EQ. 'XYROTSTR' .OR.
     *                 LBLPP(KSYMOP,IOP) .EQ. 'XZROTSTR' .OR.
     *                 LBLPP(KSYMOP,IOP) .EQ. 'YYROTSTR' .OR.
     *                 LBLPP(KSYMOP,IOP) .EQ. 'YZROTSTR' .OR.
     *                 LBLPP(KSYMOP,IOP) .EQ. 'ZZROTSTR' ) THEN
                      TRSMOM(INUM,JOP) = XMOM
                      IF (RSPECD .OR. RSPOCD) THEN
                         IF (LBLPP(KSYMOP,IOP).EQ.'XXROTSTR') THEN
                            I1 = 1
                            I2 = 1
                            KOFF = KTQVEL + 3*3*(LOCST - 1) + 3*(I2 - 1)
     &                           + I1 - 1
                            WRK(KOFF) = XMOM
                         ELSE IF (LBLPP(KSYMOP,IOP).EQ.'XYROTSTR') THEN
                            I1 = 1
                            I2 = 2
                            KOFF1 = KTQVEL + 3*3*(LOCST - 1)
     &                            + 3*(I2 - 1) + I1 - 1
                            KOFF2 = KTQVEL + 3*3*(LOCST - 1)
     &                            + 3*(I1 - 1) + I2 - 1
                            WRK(KOFF1) = XMOM
                            WRK(KOFF2) = XMOM
                         ELSE IF (LBLPP(KSYMOP,IOP).EQ.'XZROTSTR') THEN
                            I1 = 1
                            I2 = 3
                            KOFF1 = KTQVEL + 3*3*(LOCST - 1)
     &                            + 3*(I2 - 1) + I1 - 1
                            KOFF2 = KTQVEL + 3*3*(LOCST - 1)
     &                            + 3*(I1 - 1) + I2 - 1
                            WRK(KOFF1) = XMOM
                            WRK(KOFF2) = XMOM
                         ELSE IF (LBLPP(KSYMOP,IOP).EQ.'YYROTSTR') THEN
                            I1 = 2
                            I2 = 2
                            KOFF = KTQVEL + 3*3*(LOCST - 1) + 3*(I2 - 1)
     &                           + I1 - 1
                            WRK(KOFF) = XMOM
                         ELSE IF (LBLPP(KSYMOP,IOP).EQ.'YZROTSTR') THEN
                            I1 = 2
                            I2 = 3
                            KOFF1 = KTQVEL + 3*3*(LOCST - 1)
     &                            + 3*(I2 - 1) + I1 - 1
                            KOFF2 = KTQVEL + 3*3*(LOCST - 1)
     &                            + 3*(I1 - 1) + I2 - 1
                            WRK(KOFF1) = XMOM
                            WRK(KOFF2) = XMOM
                         ELSE IF (LBLPP(KSYMOP,IOP).EQ.'ZZROTSTR') THEN
                            I1 = 3
                            I2 = 3
                            KOFF = KTQVEL + 3*3*(LOCST - 1) + 3*(I2 - 1)
     &                           + I1 - 1
                            WRK(KOFF) = XMOM
                         END IF
                      END IF
                  ELSE IF ( LBLPP(KSYMOP,IOP)(1:5) .EQ. 'COS K' .OR.
     *                 LBLPP(KSYMOP,IOP)(1:5) .EQ. 'SIN K' ) THEN
                      OSCIL = D2*XMOM*XMOM*WRK(KEIVAL-1+ISIM-1+INUM)
                      IF (IPRRSP.GT.11) WRITE(LUPRI,'(/A,1P,G16.8)')
     *                   ' GENERALIZED OSCILLATOR STRENGTH :',OSCIL
                      TRSMOM(INUM,JOP) = XMOM
                      OSCSTR(INUM,JOP) = OSCIL
                  ELSE IF (LBLPP(KSYMOP,IOP)(3:8).EQ.'SPNORB' .OR.
     &                     LBLPP(KSYMOP,IOP)(3:8).EQ.'MNF-SO') THEN
                     XMOMAU = XMOM
                     TRSMOM(INUM,JOP) = XMOMAU
                     IF (IPRRSP.GT.11) THEN
                        XMOMCM = XMOMAU*XTKAYS
                        WRITE(LUPRI,'(/3A,G16.8,A/41X,G16.8,A)')
     *                  ' SPIN-ORBIT COUPLING CONSTANT (',
     *                    LBLPP(KSYMOP,IOP), '):',
     *                    XMOMAU,' (AU)', XMOMCM,' (CM-1)'
                     END IF
                  ELSE IF (LBLPP(KSYMOP,IOP)(1:3) .NE. 'HDO') THEN
                      TRSMOM(INUM,JOP) = XMOM
                  END IF
 550           CONTINUE
               IF (LBLPP(KSYMOP,IOP)(1:3) .NE. 'HDO' ) THEN
                  JOP = JOP + 1
               END IF
 600        CONTINUE
         END IF
         PRI_LINE = ' '
         WRITE (PRI_LINE,'(3(A,I2))')
     &   '@ Excit. operator sym',KSYMOP,' & ref. state sym',IREFSY,
     &   ' => excited state symmetry',KSYMST
         CALL TITLER(PRI_LINE,'*',103)
         DO 555 INUM = 1,NBX
            PRI_LINE = ' '
            WRITE (PRI_LINE,'(A,I5,A,I3)')
     *      '@ Excited state no:',(ISIM-1+INUM),' in symmetry',KSYMST
            CALL HEADER(PRI_LINE,0)
            WRITE(LUPRI,'( A,1P,G16.8,A,3(/A1,19X,G16.8,A),
     *                   //A,T21,G16.8,A)')
     *      '@ Excitation energy :',WRK(KEIVAL-1+ISIM-1+INUM),'au',
     *      '@',WRK(KEIVAL-1+ISIM-1+INUM)*XTEV,  ' eV',
     *      '@',WRK(KEIVAL-1+ISIM-1+INUM)*XTKAYS,' cm-1',
     *      ' ',WRK(KEIVAL-1+ISIM-1+INUM)*XKJMOL,' kJ / mol',
     *      '@ Total energy :',WRK(KEIVAL-1+ISIM-1+INUM)+EMCSCF,' au'
            IF (SOPPA) THEN
              TXNRM = DNORM2(KZCONF,WRK(KBVECS+(INUM-1)*KZYVAR),1)
              TYNRM = DNORM2(KZCONF,WRK(KBVECS+(INUM-1)*KZYVAR+KZVAR),1)
              T2P2HN = D100*(TXNRM*TXNRM + TYNRM*TYNRM)
              TPHNRM = D100 - T2P2HN
              WRITE(LUPRI,'(2(/A,F8.2,A))')
     &       '@ SOPPA  p-h  weight in excitation operator:',TPHNRM,' %',
     &       '@ SOPPA 2p-2h weight in excitation operator:',T2P2HN,' %'
            END IF
            IF ( NGPPP(KSYMOP) .GT. 0 ) THEN
               JOP = 0
               DO 545 IOP = 1,NGPPP(KSYMOP)
               IF (LBLPP(KSYMOP,IOP)(1:3) .EQ. 'HDO' ) GO TO 545
C                 ... skip HDO
C
                  JOP = JOP + 1
                  WRITE(LUPRI,'(/A,A)')
     *                  '@ Operator type:    ',LBLPP(KSYMOP,IOP)
C
                  IF ( LBLPP(KSYMOP,IOP)(2:8) .EQ. 'DIPLEN ' ) THEN
                     WRITE(LUPRI,'(A,1P,G16.8,A,G16.8,A)')
     *                  '@ Oscillator strength (LENGTH)   :',
     *                  OSCSTR(INUM,JOP),
     *                  '  (Transition moment :',TRSMOM(INUM,JOP),' )'
                  ELSE IF ( LBLPP(KSYMOP,IOP)(2:8) .EQ. 'DIPVEL ' ) THEN
                     WRITE(LUPRI,'(A,1P,G16.8,A,G16.8,A)')
     *                  '@ Oscillator strength (VELOCITY) :',
     *                  OSCSTR(INUM,JOP),
     *                  '  (Transition moment :',TRSMOM(INUM,JOP),' )'
                  ELSE IF ( LBLPP(KSYMOP,IOP)(1:5) .EQ. 'COS K' .OR.
     &                      LBLPP(KSYMOP,IOP)(1:5) .EQ. 'SIN K' ) THEN
                     WRITE(LUPRI,'(A,1P,G16.8)')
     *                  '@ Generalized oscillator strength:',
     *                  OSCSTR(INUM,JOP)
                     IF ( LBLPP(KSYMOP,IOP)(1:5) .EQ. 'COS K' ) THEN
                        WRITE(LUPRI,'(A,1P,G16.8)')
     *                  '@ Transition moment  (cos(kr)/k) :',
     *                  TRSMOM(INUM,JOP)
                     ELSE
                        WRITE(LUPRI,'(A,1P,G16.8)')
     *                  '@ Transition moment  (sin(kr)/k) :',
     *                  TRSMOM(INUM,JOP)
                     END IF
                  ELSE IF (LBLPP(KSYMOP,IOP)(3:8).EQ.'SPNORB' .OR.
     &                     LBLPP(KSYMOP,IOP)(3:8).EQ.'MNF-SO' .OR.
     &                     LBLPP(KSYMOP,IOP)(3:8).EQ.'SPNSCA') THEN
                     WRITE(LUPRI,'(3A,1P,G16.8,A)')
     *                  '@ Spin-orbit coupling constant (',
     *                    LBLPP(KSYMOP,IOP),') :',
     *               TRSMOM(INUM,JOP)*XTKAYS,' cm-1'
                     WRITE(LUPRI,'(A,1P,G16.8,A)')
     *                  '                                          ',
     *               TRSMOM(INUM,JOP),' au'
                  ELSE
                     WRITE(LUPRI,'(A,1P,G16.8)')
     *                  '@ Transition moment              :',
     *                  TRSMOM(INUM,JOP)
                  END IF
  545          CONTINUE
            END IF
            IF (IPRRSP .GT. 3) THEN
               WRITE (LUPRI,'(/A,I3)')
     &              ' Eigenvector for state no.',(ISIM-1+INUM)
               CALL RSPPRO(WRK(KBVECS+(INUM-1)*KZYVAR+KZCONF),KZVAR,
     *                     LUPRI)
C              CALL RSPPRC(WRK(KBVECS+(INUM-1)*KZYVAR),KZCONF,KZVAR,LUPRI)
               CALL RSPANC(WRK(KBVECS+(INUM-1)*KZYVAR),KZCONF,KZVAR,
     *                     MULD2H(KSYMOP,IREFSY),XINDX,MULD2H,LUPRI)
            END IF
            IF (DOMO2) CALL RSPPMO(WRK(KBVECS+(INUM-1)*KZYVAR+KZCONF),
     &                             KZVAR,WRK(KCMO2),LUPRI)
            CALL WRTRSP(LURSP,KZYVAR,WRK(KBVECS+(INUM-1)*KZYVAR),
     &                  'EXCITLAB',BLANK,WRK(KEIVAL-1+ISIM-1+INUM),D0,
     &                  KSYMOP,0,WRK(KRESID-1+ISIM-1+INUM),D0)
 555     CONTINUE
 500  CONTINUE
C
C     Print oscillator/rotatory strengths
C
      IF (RSPOCD .OR. RSPECD) THEN

         CALL AROUND('Absorption and Circular Dichroism Spectra')

         WRITE (LUPRI,'(/,3(/,5X,A,I3))')
     &   ' Reference  state    symmetry:',IREFSY,
     &   ' Excitation operator symmetry:',KSYMOP,
     &   ' Excited    state    symmetry:',KSYMST

         WRITE (LUPRI,'(/,5X,A,/5X,A,//)')
     &      ' Units: 10**(-40) (esu**2)*(cm**2) (rotatory   strength)',
     &      '        dimensionless              (oscillator strength)'

         CALL HEADER('Oscillator and Scalar Rotatory Strength',18)
         WRITE(LUPRI,'(/,8X,A,/,8X,A,/,8X,A)')
     &   ' State Frequency   Oscillator Strength    Rotatory  Strength',
     &   '         (eV)      Velocity     Length    Velocity    Length',
     &   ' -----------------------------------------------------------'

         DO LOCST = 1,NPPCNV(KSYMOP)

            IOFF0  = LOCST - 1
            IOFF1  = 3*IOFF0
            EXENG  = WRK(KFREX+IOFF0)
            EXEV   = EXENG*XTEV

            OSCLEN = D2R3*DDOT(3,WRK(KTRLEN+IOFF1),1,
     &                           WRK(KTRLEN+IOFF1),1)*EXENG
            ROTLEN = 0.5D0*DDOT(3,WRK(KTRLEN+IOFF1),1,
     &                            WRK(KTRMAG+IOFF1),1)

            IF (ABS(EXENG) .GT. 0.0D0) THEN

               OSCVEL = D2R3*DDOT(3,WRK(KTRVEL+IOFF1),1,
     &                              WRK(KTRVEL+IOFF1),1)/EXENG
               ROTVEL = -0.5D0*DDOT(3,WRK(KTRVEL+IOFF1),1,
     &                                WRK(KTRMAG+IOFF1),1)/EXENG


               WRITE(LUPRI,'(8X,I5,5(1X,F10.4))')
     &         LOCST,EXEV,OSCVEL,OSCLEN,ROTVEL*ESUECD,ROTLEN*ESUECD

            ELSE

             WRITE(LUPRI,'(8X,I5,1X,F10.4,1X,A,1X,F10.4,1X,A,1X,F10.4)')
     &         LOCST,EXEV,' undefined',OSCLEN,' undefined',ROTLEN*ESUECD

               ROTVEL = 0.0D0

            END IF

            IF (RSPOCD) THEN
               WRK(KROTV+IOFF0) = ROTVEL*ESUECD
               WRK(KROTL+IOFF0) = ROTLEN*ESUECD
            END IF

         END DO

      END IF

      IF (RSPOCD) THEN

         LOCWRN = 0
         NUMST  = MIN(LMXFR,NPPCNV(KSYMOP))

         CALL DZERO(RMVEL,3*3*NUMST)
         CALL DZERO(RMLEN,3*3*NUMST)
         CALL DZERO(RQVEL,3*3*NUMST)
         CALL DZERO(RQLEN,3*3*NUMST)
         CALL DZERO(RTVEL,3*3*NUMST)
         CALL DZERO(RTLEN,3*3*NUMST)

         CALL HEADER('Magnetic Dipole Rotatory Strength Tensor',18)
         WRITE(LUPRI,'(/,A,/,A,A,/,A,A)')
     &   ' State Gauge       Magnetic Dipole Rotatory Strength Tensor ',
     &   '                  xx        xy        xz        yy        yz',
     &   '        zz',
     &   ' -----------------------------------------------------------',
     &   '------------'
         DO LOCST = 1,NUMST

            IOFF0  = LOCST - 1
            IOFF1  = 3*IOFF0
            EXENG  = WRK(KFREX+IOFF0)

            PDOTL  = DDOT(3,WRK(KTRVEL+IOFF1),1,
     &                      WRK(KTRMAG+IOFF1),1)
            RDOTL  = DDOT(3,WRK(KTRLEN+IOFF1),1,
     &                      WRK(KTRMAG+IOFF1),1)

            DO J = 1,3
               DO I = 1,3
                  RMVEL(I,J,LOCST) =
     &                       WRK(KTRVEL+IOFF1+J-1)*WRK(KTRMAG+IOFF1+I-1)
                  RMLEN(I,J,LOCST) =
     &                      -WRK(KTRLEN+IOFF1+J-1)*WRK(KTRMAG+IOFF1+I-1)
                  IF (I .EQ. J) THEN
                     RMVEL(I,I,LOCST) = RMVEL(I,I,LOCST) - PDOTL
                     RMLEN(I,I,LOCST) = RMLEN(I,I,LOCST) + RDOTL
                  END IF
               END DO
            END DO

            DO J = 1,3
               DO I = 1,J
                  RMVEL(I,J,LOCST) = RMVEL(I,J,LOCST) + RMVEL(J,I,LOCST)
                  RMVEL(J,I,LOCST) = RMVEL(I,J,LOCST)
                  RMLEN(I,J,LOCST) = RMLEN(I,J,LOCST) + RMLEN(J,I,LOCST)
                  RMLEN(J,I,LOCST) = RMLEN(I,J,LOCST)
               END DO
            END DO

            FACL = 0.375D0*ESUECD
            CALL DSCAL(3*3,FACL,RMLEN(1,1,LOCST),1)

            IF (ABS(EXENG) .GT. 0.0D0) THEN

               FACV = FACL/EXENG
               CALL DSCAL(3*3,FACV,RMVEL(1,1,LOCST),1)

               WRITE(LUPRI,'(I5,3X,A,6(1X,F9.4))')
     &         LOCST,'Vel ',RMVEL(1,1,LOCST),RMVEL(1,2,LOCST),
     &                      RMVEL(1,3,LOCST),
     &                      RMVEL(2,2,LOCST),RMVEL(2,3,LOCST),
     &                      RMVEL(3,3,LOCST)

            ELSE

               CALL DZERO(RMVEL(1,1,LOCST),3*3)
               WRITE(LUPRI,'(I5,3X,A,1X,A)')
     &         LOCST,'Vel ','    ------ undefined ------'

            END IF

            TSTVEL = RMVEL(1,1,LOCST)
            TSTLEN = RMLEN(1,1,LOCST)
            DO I = 2,3
               TSTVEL = TSTVEL + RMVEL(I,I,LOCST)
               TSTLEN = TSTLEN + RMLEN(I,I,LOCST)
            END DO
            TSTVEL = TSTVEL/3.0D0
            TSTLEN = TSTLEN/3.0D0

            IF (ABS(TSTVEL-WRK(KROTV+IOFF0)) .GT. TOLRTT) THEN
               WRITE(LUPRI,'(1X,A,1P,D16.8,/,1X,A,1P,D16.8)')
     &         'WARNING: wrong average: ',TSTVEL,
     &         '         expected     : ',WRK(KROTV+IOFF0)
               WRITE(LUPRI,'(1X,A,1P,D16.8)')
     &         '         difference   : ',TSTVEL-WRK(KROTV+IOFF0)
               LOCWRN = LOCWRN + 1
            END IF

            WRITE(LUPRI,'(8X,A,6(1X,F9.4))')
     &            'Len ',RMLEN(1,1,LOCST),RMLEN(1,2,LOCST),
     &                   RMLEN(1,3,LOCST),
     &                   RMLEN(2,2,LOCST),RMLEN(2,3,LOCST),
     &                   RMLEN(3,3,LOCST)

            IF (ABS(TSTLEN-WRK(KROTL+IOFF0)) .GT. TOLRTT) THEN
               WRITE(LUPRI,'(1X,A,1P,D16.8,/,1X,A,1P,D16.8)')
     &         'WARNING: wrong average: ',TSTLEN,
     &         '         expected     : ',WRK(KROTL+IOFF0)
               WRITE(LUPRI,'(1X,A,1P,D16.8)')
     &         '         difference   : ',TSTLEN-WRK(KROTL+IOFF0)
               LOCWRN = LOCWRN + 1
            END IF

         END DO

         CALL HEADER('Electric Quadrupole Rotatory Strength Tensor',18)
         WRITE(LUPRI,'(/,A,/,A,A,/,A,A)')
     &   ' State Gauge   Electric Quadrupole Rotatory Strength Tensor ',
     &   '                  xx        xy        xz        yy        yz',
     &   '        zz',
     &   ' -----------------------------------------------------------',
     &   '------------'
         DO LOCST = 1,NUMST

            IOFF0  = LOCST - 1
            IOFF1  = 3*IOFF0
            IOFF2  = 3*IOFF1
            EXENG  = WRK(KFREX+IOFF0)

            DO K = 1,3
               DO J = 1,3
                  DO L = 1,3
                     DO M = 1,3
                        RQVEL(J,K,LOCST) =
     &                  RQVEL(J,K,LOCST) - XLEVICI(M,L,J)*
     &                  WRK(KTRVEL+IOFF1+L-1)*
     &                  WRK(KTQVEL+IOFF2+3*(K-1)+M-1)
                        RQLEN(J,K,LOCST) =
     &                  RQLEN(J,K,LOCST) - XLEVICI(M,L,J)*
     &                  WRK(KTRLEN+IOFF1+L-1)*
     &                  WRK(KTQLEN+IOFF2+3*(K-1)+M-1)
                     END DO
                  END DO
               END DO
            END DO

            DO J = 1,3
               DO I = 1,J
                  RQVEL(I,J,LOCST) = RQVEL(I,J,LOCST)
     &                             + RQVEL(J,I,LOCST)
                  RQVEL(J,I,LOCST) = RQVEL(I,J,LOCST)
                  RQLEN(I,J,LOCST) = RQLEN(I,J,LOCST)
     &                             + RQLEN(J,I,LOCST)
                  RQLEN(J,I,LOCST) = RQLEN(I,J,LOCST)
               END DO
            END DO

            FACL = -0.375D0*ESUECD*EXENG
            CALL DSCAL(3*3,FACL,RQLEN(1,1,LOCST),1)

            IF (ABS(EXENG) .GT. 0.0D0) THEN

               FACV = -0.375D0*ESUECD/EXENG
               CALL DSCAL(3*3,FACV,RQVEL(1,1,LOCST),1)

               WRITE(LUPRI,'(I5,3X,A,6(1X,F9.4))')
     &         LOCST,'Vel ',RQVEL(1,1,LOCST),RQVEL(1,2,LOCST),
     &                      RQVEL(1,3,LOCST),
     &                      RQVEL(2,2,LOCST),RQVEL(2,3,LOCST),
     &                      RQVEL(3,3,LOCST)

            ELSE

               CALL DZERO(RQVEL(1,1,LOCST),3*3)
               WRITE(LUPRI,'(I5,3X,A,1X,A)')
     &         LOCST,'Vel ','    ------ undefined ------'

            END IF

            TSTVEL = RQVEL(1,1,LOCST)
            TSTLEN = RQLEN(1,1,LOCST)
            DO I = 2,3
               TSTVEL = TSTVEL + RQVEL(I,I,LOCST)
               TSTLEN = TSTLEN + RQLEN(I,I,LOCST)
            END DO
            TSTVEL = TSTVEL/3.0D0
            TSTLEN = TSTLEN/3.0D0

            IF (ABS(TSTVEL) .GT. TOLRTT) THEN
               WRITE(LUPRI,'(1X,A,1P,D16.8)')
     &         'WARNING: wrong average: ',TSTVEL
               LOCWRN = LOCWRN + 1
            END IF

            WRITE(LUPRI,'(8X,A,6(1X,F9.4))')
     &            'Len ',RQLEN(1,1,LOCST),RQLEN(1,2,LOCST),
     &                   RQLEN(1,3,LOCST),
     &                   RQLEN(2,2,LOCST),RQLEN(2,3,LOCST),
     &                   RQLEN(3,3,LOCST)

            IF (ABS(TSTLEN) .GT. TOLRTT) THEN
               WRITE(LUPRI,'(1X,A,1P,D16.8)')
     &         'WARNING: wrong average: ',TSTLEN
               LOCWRN = LOCWRN + 1
            END IF

         END DO

         CALL HEADER('Total Rotatory Strength Tensor',18)
         WRITE(LUPRI,'(/,A,/,A,A,/,A,A)')
     &   ' State Gauge                 Total Rotatory Strength Tensor ',
     &   '                  xx        xy        xz        yy        yz',
     &   '        zz',
     &   ' -----------------------------------------------------------',
     &   '------------'
         LENTOT = 3*3*NPPCNV(KSYMOP)
         CALL DCOPY(LENTOT,RMVEL,1,RTVEL,1)
         CALL DAXPY(LENTOT,1.0D0,RQVEL,1,RTVEL,1)
         CALL DCOPY(LENTOT,RMLEN,1,RTLEN,1)
         CALL DAXPY(LENTOT,1.0D0,RQLEN,1,RTLEN,1)
         DO LOCST = 1,NUMST

            IOFF0 = LOCST - 1

            IF (ABS(EXENG) .GT. 0.0D0) THEN

               WRITE(LUPRI,'(I5,3X,A,6(1X,F9.4))')
     &         LOCST,'Vel ',RTVEL(1,1,LOCST),RTVEL(1,2,LOCST),
     &                      RTVEL(1,3,LOCST),
     &                      RTVEL(2,2,LOCST),RTVEL(2,3,LOCST),
     &                      RTVEL(3,3,LOCST)

            ELSE

               CALL DZERO(RTVEL(1,1,LOCST),3*3)
               WRITE(LUPRI,'(I5,3X,A,1X,A)')
     &         LOCST,'Vel ','    ------ undefined ------'

            END IF

            TSTVEL = RTVEL(1,1,LOCST)
            TSTLEN = RTLEN(1,1,LOCST)
            DO I = 2,3
               TSTVEL = TSTVEL + RTVEL(I,I,LOCST)
               TSTLEN = TSTLEN + RTLEN(I,I,LOCST)
            END DO
            TSTVEL = TSTVEL/3.0D0
            TSTLEN = TSTLEN/3.0D0

            IF (ABS(TSTVEL-WRK(KROTV+IOFF0)) .GT. TOLRTT) THEN
               WRITE(LUPRI,'(1X,A,1P,D16.8,/,1X,A,1P,D16.8)')
     &         'WARNING: wrong average: ',TSTVEL,
     &         '         expected     : ',WRK(KROTV+IOFF0)
               WRITE(LUPRI,'(1X,A,1P,D16.8)')
     &         '         difference   : ',TSTVEL-WRK(KROTV+IOFF0)
               LOCWRN = LOCWRN + 1
            END IF

            WRITE(LUPRI,'(8X,A,6(1X,F9.4))')
     &            'Len ',RTLEN(1,1,LOCST),RTLEN(1,2,LOCST),
     &                   RTLEN(1,3,LOCST),
     &                   RTLEN(2,2,LOCST),RTLEN(2,3,LOCST),
     &                   RTLEN(3,3,LOCST)

            IF (ABS(TSTLEN-WRK(KROTL+IOFF0)) .GT. TOLRTT) THEN
               WRITE(LUPRI,'(1X,A,1P,D16.8,/,1X,A,1P,D16.8)')
     &         'WARNING: wrong average: ',TSTLEN,
     &         '         expected     : ',WRK(KROTL+IOFF0)
               WRITE(LUPRI,'(1X,A,1P,D16.8)')
     &         '         difference   : ',TSTLEN-WRK(KROTL+IOFF0)
               LOCWRN = LOCWRN + 1
            END IF

         END DO

         IF (LOCWRN .NE. 0) THEN
            NWARN = NWARN + 1
            WRITE(LUPRI,'(//,1X,A,I6,A,A,/,1X,A,//)')
     &      'WARNING:',LOCWRN,' warnings were issued for rotatory ',
     &      ' strength tensors.',
     &      '         Averages are incorrect!!!!'
         END IF

      END IF
C
C *** END OF RSPPP
C
      RETURN
      END
C  /* Deck rsplr */
      SUBROUTINE RSPLR(CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
C
C     Driver Routine for Linear Response:
C     DETERMINE SECOND ORDER MOLECULAR PROPERTIES
C
#include <implicit.h>
#include <iratdef.h>
#include <thrzer.h>
C
      CHARACTER*8 BLANK
      PARAMETER (D0 = 0.0D0, BLANK = '        ')
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
#include <priunit.h>
#include <infpri.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <infvar.h>
#include <inftap.h>
#include <rspprp.h>
#include <inflr.h>
#include <infdim.h>
#include <inforb.h>
C
C Local variables
C
      LOGICAL OLDDX, FOUND, CONV, NOSEL
C
C
      NOSEL  = .FALSE.
C
C     NOSEL is activated if we have a converged response vector on file,
C     but then obviously no residuals. We calculate in this case the
C     property without the use of Sellers formula.
C
      KREDE  = 1
      KREDS  = KREDE  + MAXRM*MAXRM
      KIBTYP = KREDS  + MAXRM*MAXRM
      KEIVAL = KIBTYP + MAXRM
      KRESID = KEIVAL + MAXRM
      KEIVEC = KRESID + MAXRM
      KREDGD = KEIVEC + MAXRM*MAXRM
      KGD    = KREDGD + MAXRM*NFREQ
      KWRK1  = KGD    + KZYVAR
      LWRK1  = LWRK   - KWRK1
C
C WORK SPACE FOR RSPEVE AND RSPRES
C
C WORK SPACE FOR EACH RESIDUAL VECTOR
C
      MRES   = MAX(NORBT*NORBT,2*N2ASHX) + KZYVAR
C
C WORK SPACE EACH TIME RSPRES IS CALLED
C
      MWRK   = KZYVAR
      IF (.NOT. RSPCI ) MWRK = MWRK + MAX(KZCONF,NCREF) +  LACIMX
C NUMBER OF SIMULTANEOUS VECTORS
      MXFREQ = (LWRK1- MWRK)/MRES
      MXFREQ = MIN(NFREQ,MXFREQ)
      KWRKE  = KWRK1
      KBVECS = KWRKE + KZYVAR
C
      IF (MXFREQ .LE. 0) THEN
        WRITE (LUERR,9100) LWRK1,KWRK1+MWRK+MRES
        WRITE (LUPRI,9100) LWRK1,KWRK1+MWRK+MRES
        CALL QTRACE(LUERR)
        CALL QUIT('RSPLR: INSUFFICIENT SPACE TO SOLVE LINEAR EQUATIONS')
      ENDIF
 9100 FORMAT(/' RSPLR, work space too small for RSPRES',
     *       /'        had:',I10,', need more than:',I10)
C
C     Open LURSP2 for NSOLVC solution vectors and NSOLVC
C        residual vectors.
C
      NSOLVC = NFREQ * NGPLR(KSYMOP)
C
      LURSP2 = -1
      IF (NSOLVC .GE. 1) CALL GPOPEN(LURSP2,'RSPSOLVEC.DA','NEW',
     &     'DIRECT',' ',IRAT*KZYVAR,OLDDX)
C
C
      KZRED  = 0
      KZYRED = 0
      THCRSP = THCLR
      IPRRSP = IPRLR
      MAXIT  = MAXITL
C
C     Call RSPCTL to solve linear set of response equations
C
      DO 600 IOP = 1,NGPLR(KSYMOP)
         WRITE (LUPRI,'(//A,I3,/2A,/A,(T25,5F10.6))')
     &   ' RSPLR -- linear response calculation for symmetry',KSYMOP,
     &   ' RSPLR -- operator label : ',LBLLR(KSYMOP,IOP),
     &   ' RSPLR -- frequencies :',(FREQ(I),I=1,NFREQ)
         CALL GETGPV(LBLLR(KSYMOP,IOP),FC,FV,CMO,UDV,PV,XINDX,ANTSYM,
     *               WRK(KGD),LWRK1)
         DNORM = DNRM2(KZYVAR,WRK(KGD),1)
         write(lupri,*) 'FREQ 2 ',FREQ(1:NFREQ)
         IF (IPRRSP .GT. 30) THEN
            CALL HEADER('Gradient vectors in RSPLR for operator '//
     &                  LBLLR(KSYMOP,IOP)//':',-1)
            WRITE (LUPRI,'(/A,1P,D12.4)') ' Norm : ', DNORM
            CALL OUTPUT(WRK(KGD),1,KZVAR,1,2,KZVAR,2,1,LUPRI)
         END IF
         JFREQ  = 0
         CALL DZERO(WRK(KEIVAL),NFREQ)
         DO 410 I = 1, NFREQ
            AFRQ = ABS(FREQ(I))
            CALL REARSP(LURSP,KLEN,WRK(KBVECS),LBLLR(KSYMOP,IOP),
     &                  BLANK,FREQ(I),D0,KSYMOP,0,THCLR,
     &                  FOUND,CONV,ANTSYM)
            IF (FOUND .AND. CONV) THEN
               WRITE(LUPRI,'(/A,/A12,A10,F22.6,/A42)')
     *              ' Converged solution vector already on file',
     *               LBLLR(KSYMOP,IOP),' freq',FREQ(I),
     *              ' -----------------------------------------'
            ELSE
C              ! Remove duplicates, they will cause linear dependency
               DO IJ = 0, JFREQ-1
                  IF (ABS(WRK(KEIVAL+IJ)-AFRQ) .LT. THRZER) GOTO 410
               END DO
               WRK(KEIVAL+JFREQ) = FREQ(I)
               JFREQ = JFREQ + 1
            END IF
 410     CONTINUE
         IF (JFREQ .EQ. 0) THEN
            ! Found converged solution vectors, Sellers formula not possible
            NOSEL = .TRUE.
            GOTO 600
         END IF
C
         IF (DNORM .LT. THRNRM) THEN
            WRITE (LUPRI,*) ' >>> RSPCTL skipped because norm of'
            WRITE (LUPRI,*) '     property vector '//LBLLR(KSYMOP,IOP)//
     &                      ' is only',DNORM
            CALL DZERO(WRK(KBVECS),KZYVAR)
            DO 450 IFREQ = 1,JFREQ
               CALL WRTRSP(LURSP,KZYVAR,WRK(KBVECS),LBLLR(KSYMOP,IOP),
     *                     BLANK,WRK(KEIVAL-1+IFREQ),D0,KSYMOP,0,D0,
     *                     ANTSYM)
               IFRSOL = (IOP-1)*NFREQ*2 + (IFREQ-1)*2 + 1
               IFRRES = IFRSOL + 1
               CALL WRITDX(LURSP2,IFRSOL,IRAT*KZYVAR,WRK(KBVECS))
               CALL WRITDX(LURSP2,IFRRES,IRAT*KZYVAR,WRK(KBVECS))
  450       CONTINUE
            GO TO 600
         END IF
C
         KEXSIM = JFREQ
         KEXCNV = JFREQ
C
         CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *               .TRUE.,LBLLR(KSYMOP,IOP),BLANK,WRK(KGD),
     *               WRK(KREDGD),WRK(KREDE),WRK(KREDS),
     *               WRK(KIBTYP),WRK(KEIVAL),WRK(KRESID),WRK(KEIVEC),
     *               XINDX,WRK(KWRK1),LWRK1)
C        CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
C                    LINEQ,GD,REDGD,REDE,REDS,
C    *               IBTYP,EIVAL,EIVEC,XINDX,WRK,LWRK)
C
         DO 700 IFREQ = 1,JFREQ,MXFREQ
            NBX = MIN(MXFREQ,(JFREQ+1-IFREQ))
            IBOFF = IFREQ - 1
            CALL RSPEVE(WRK(KIBTYP),WRK(KEIVAL),WRK(KEIVEC),
     *                  WRK(KBVECS),WRK(KWRKE),NBX,IBOFF)
C           CALL RSPEVE(IBTYP,EIVAL,EIVEC,BVECS,WRK,NBX,IBOFF)
            DO 750 IVEC = 1,NBX
               IBV    = (IVEC-1)*KZYVAR + KBVECS
               NVEC   = (IOP-1)*2*NFREQ + (IBOFF+IVEC-1)*2 + 1
               IF (IPRRSP .GT. 3) THEN
                  WRITE (LUPRI,'(/A,I5)')' EIGENVECTOR NUMBER ',
     &                 IFREQ-1+IVEC
C                 CALL RSPPRC(WRK(KBVECS+(IVEC-1)*KZYVAR),KZCONF,KZVAR,LUPRI)
                  CALL RSPANC(WRK(KBVECS+(IVEC-1)*KZYVAR),KZCONF,KZVAR,
     *                        MULD2H(KSYMOP,IREFSY),XINDX,MULD2H,LUPRI)
                  CALL RSPPRO(WRK(KBVECS+(IVEC-1)*KZYVAR+KZCONF),KZVAR,
     *                        LUPRI)
               END IF
               CALL WRITDX(LURSP2,NVEC,IRAT*KZYVAR,WRK(IBV))
               ANTSYM = 1.0d0
               CALL WRTRSP(LURSP,KZYVAR,WRK(IBV),LBLLR(KSYMOP,IOP),
     &                     BLANK,WRK(KEIVAL-1+IVEC),D0,KSYMOP,0,
     &                     WRK(KRESID-1+IVEC),ANTSYM)
               IF (IPRRSP .GT. 30) THEN
                  CALL HEADER('Solution vectors in RSPLR for operator '
     &                        //LBLLR(KSYMOP,IOP)//':',-1)
                  CALL OUTPUT(WRK(IBV),1,KZVAR,1,2,KZVAR,2,1,LUPRI)
                  DNORM = DNRM2(KZYVAR,WRK(IBV),1)
                  WRITE (LUPRI,'(A,D20.14)') ' Norm:', DNORM
               END IF
  750       CONTINUE
            CALL RSPRES(WRK(KIBTYP),WRK(KEIVAL),WRK(KEIVEC),
     *                  WRK(KGD),WRK(KWRK1),
     *                  XINDX,LWRK1,UDV,NBX,IBOFF)
C           CALL RSPRES(IBTYP,EIVAL,EIVEC,GD,RD,
C    *                  XINDX,LWRK,UDV,NBX,IBOFF)
            DO 760 IVEC = 1,NBX
               IBV    = (IVEC-1)*KZYVAR + KWRK1
               NVEC   = (IOP-1)*2*NFREQ + (IBOFF+IVEC-1)*2 + 2
               CALL WRITDX(LURSP2,NVEC,IRAT*KZYVAR,WRK(IBV))
               IF (IPRRSP .GT. 30) THEN
                  CALL HEADER('Residuum vectors in RSPLR for operator '
     &                        //LBLLR(KSYMOP,IOP)//':',-1)
                  CALL OUTPUT(WRK(IBV),1,KZVAR,1,2,KZVAR,2,1,LUPRI)
                  DNORM = DNRM2(KZYVAR,WRK(IBV),1)
                  WRITE (LUPRI,'(A,D20.14)') ' Norm:', DNORM
               END IF
 760        CONTINUE
  700    CONTINUE
  600 CONTINUE
      IF (LWRK1.LT.0) CALL ERRWRK('RSPLR 2',KWRK1-1,LWRK)
      IF (NOSEL) THEN
         KVECB = 1
         KWRK1 = KVECB + KZYVAR
         LWRK1 = LWRK - KWRK1
         NINFO = NINFO + 1
         WRITE (LUPRI,'(/A)')'@ INFO >>>>> Linear response functions'//
     &        ' calculated using presolved response vectors',
     &        '@ INFO >>>>> These results will not have quadratic '//
     &        'accuracy'
         CALL LRNSL(CMO,UDV,PV,FC,XINDX,WRK(KVECB),WRK(KWRK1),LWRK1)
      ELSE
         NOP   = NGPLR(KSYMOP)
         KSOPR = 1
         KWRK1 = KSOPR + NFREQ*NOP*NOP
         LWRK1 = LWRK  - KWRK1
         IF (LWRK1.LT.0) CALL ERRWRK('RSPLR 2',KWRK1-1,LWRK)
         CALL LRSOPR(NOP,WRK(KSOPR),FC,CMO,UDV,PV,XINDX,LURSP2,
     &               WRK(KWRK1),LWRK1)
      END IF
C
C     Close and delete LURSP2
C
      CALL GPCLOSE(LURSP2,'DELETE')
C
C *** end of RSPLR --
C
      RETURN
      END
C  /* Deck rds2 */
      SUBROUTINE RDS2(LABEL,UDV,XINDX,WRK,LWRK)
C 28-Jun-1991 PJ+HJJ
#include <implicit.h>
#include <iratdef.h>
      CHARACTER*8 LABEL
      DIMENSION XINDX(*),UDV(*),WRK(*)
      PARAMETER ( D1 = 1.0D0 )
C
C Used from common blocks:
C   INFDIM : NVARMA
C
#include <priunit.h>
#include <inftap.h>
#include <infdim.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
      LOGICAL OLDDX
C
      KS2LIN = 1
      KRD    = KS2LIN + KZYVAR
      KWRK1  = KRD    + KZVAR + KZWOPT
      LWRK1  = LWRK   - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('RDS2',KWRK1-1,LWRK)
      READ(LABEL(4:6),'(I3)') NLAB
      NREC = 2*NLAB -1
C
C     Open LURD file from ABACUS
C
      LURD = -1
      CALL GPOPEN(LURD,ABARDR,'OLD','DIRECT',' ',IRAT*NVARMA,OLDDX)
C     Read RD_Z from ABACUS
      CALL READDX(LURD,NREC,IRAT*KZVAR,WRK(KRD))
      CALL GPCLOSE(LURD,'KEEP')
      IF (IPRRSP .GT. 120) THEN
         WRITE (LUPRI,*) 'RDS2 test print of ABACUS RD vector'
         WRITE (LUPRI,*) 'label : ',LABEL
         WRITE (LUPRI,*) 'NREC  : ',NREC
         WRITE (LUPRI,*) 'NVARMA: ',NVARMA
         WRITE (LUPRI,*) 'KZVAR : ',KZVAR
         CALL OUTPUT(WRK(KRD),1,KZVAR,1,1,KZVAR,1,1,LUPRI)
      END IF
C     Calculate S[2] * RD where we use RD_Y = -RD_Z
      CALL DZERO(WRK(KS2LIN),KZYVAR)
      IF (KZCONF.GT.0) THEN
         CALL RSPSLI(1,0,WRK(KRD),DUMMY,UDV,WRK(KS2LIN),
     *   XINDX,WRK(KWRK1),LWRK1)
      END IF
      IF (KZWOPT.GT.0) THEN
         CALL DZERO(WRK(KRD+KZVAR),KZWOPT)
         CALL RSPSLI(0,1,DUMMY,WRK(KRD+KZCONF),UDV,WRK(KS2LIN),
     *   XINDX,WRK(KWRK1),LWRK1)
      END IF
      CALL DAXPY(KZVAR,D1,WRK(KS2LIN),1,WRK(KS2LIN+KZVAR),1)
      CALL DCOPY(KZVAR,WRK(KS2LIN+KZVAR),1,WRK(KS2LIN),1)
      RETURN
      END
C  /* Deck rds2x */
      SUBROUTINE RDS2X(NLAB,UDV,XINDX,WRK,LWRK)
C Made from RDS2 by klb feb. 1997.
#include <implicit.h>
#include <iratdef.h>
      DIMENSION XINDX(*),UDV(*),WRK(*)
      PARAMETER ( D1 = 1.0D0 )
C
C Used from common blocks:
C   INFDIM : NVARMA
C
#include <priunit.h>
#include <inftap.h>
#include <infdim.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
      LOGICAL OLDDX
C
      KS2LIN = 1
      KRD    = KS2LIN + KZYVAR
      KWRK1  = KRD    + KZVAR + KZWOPT
      LWRK1  = LWRK   - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('RDS2X',KWRK1-1,LWRK)
      NREC = 2*NLAB -1
C
C     Open LURD file from ABACUS
C   
      LURD = -1
      CALL GPOPEN(LURD,ABARDR,'OLD','DIRECT',' ',IRAT*NVARMA,OLDDX)
C     Read RD_Z from ABACUS
      CALL READDX(LURD,NREC,IRAT*KZVAR,WRK(KRD))
      CALL GPCLOSE(LURD,'KEEP')
      IF (IPRRSP .GT. 120) THEN
         WRITE (LUPRI,*) 'RDS2X test print of ABACUS RD vector'
         WRITE (LUPRI,*) 'NLAB  : ',NLAB
         WRITE (LUPRI,*) 'NREC  : ',NREC
         WRITE (LUPRI,*) 'NVARMA: ',NVARMA
         WRITE (LUPRI,*) 'KZVAR : ',KZVAR
         CALL OUTPUT(WRK(KRD),1,KZVAR,1,1,KZVAR,1,1,LUPRI)
      END IF
C     Calculate S[2] * RD where we use RD_Y = -RD_Z
      CALL DZERO(WRK(KS2LIN),KZYVAR)
      IF (KZCONF.GT.0) THEN
         CALL RSPSLI(1,0,WRK(KRD),DUMMY,UDV,WRK(KS2LIN),
     *   XINDX,WRK(KWRK1),LWRK1)
      END IF
      IF (KZWOPT.GT.0) THEN
         CALL DZERO(WRK(KRD+KZVAR),KZWOPT)
         CALL RSPSLI(0,1,DUMMY,WRK(KRD+KZCONF),UDV,WRK(KS2LIN),
     *   XINDX,WRK(KWRK1),LWRK1)
      END IF
      CALL DAXPY(KZVAR,D1,WRK(KS2LIN),1,WRK(KS2LIN+KZVAR),1)
      CALL DCOPY(KZVAR,WRK(KS2LIN+KZVAR),1,WRK(KS2LIN),1)
      RETURN
      END
C  /* Deck rspres */
      SUBROUTINE RSPRES(IBTYP,EIVAL,EIVEC,GD,RD,
     *                  XINDX,LWRK,UDV,NBX,IBOFF)
C
C PURPOSE: CONSTRUCT RESIDUAL RD(I) = (E(2)-W(I)*S(2))*X(I) - GD
C
C On input, vector X(I) is placed in RD(I+1), I = 1,NBX.
C On output, vector RD(I) is placed in RD(I), I = 1,NBX.
C LWRK is available work space in RD.
C
#include <implicit.h>
#include <dummy.h>
      DIMENSION IBTYP(*),EIVAL(*),EIVEC(KZYRED,*)
      DIMENSION GD(*),UDV(*),RD(*),XINDX(*)
      PARAMETER (  D0 = 0.0D0 , DM1=-1.0D0)
C
C Used from common blocks:
C  INFRSP: ?
C  WRKRSP: KZVAR,KZYVAR
C          LURSP5
C
#include <priunit.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <inftap.h>
C
#include <ibndxdef.h>
C
C CONSTRUCT -W(I) * S[2]*X(I) IN FIRST NBX ELEMENTS OF RD
C
      KDIAE = (NBX+1)*KZYVAR + 1
      IF (SOPPA) THEN
         KSOPO = KDIAE + KZCONF
         KWRK1 = KSOPO + KZYWOP
      ELSE
         KSOPO = KDIAE
         KWRK1 = KSOPO
      ENDIF
      LWRK1 = LWRK - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('RSPRES',KWRK1-1,LWRK)
      IF (SOPPA) THEN
         CALL REWSPL(LURSP4)
#if defined (VAR_SPLITFILES)
         JRSP4 = 0
         CALL READST(LURSP4,LRSP4,RD(KDIAE),JRSP4,.TRUE.)
         IF (LRSP4 .NE. KZCONF) CALL QUIT('Inconsistent record length'//
     &        ' on unit LURSP4 in RSPRES')
#else
         CALL READT(LURSP4,KZCONF,RD(KDIAE))
#endif
      ENDIF
      DO 50 IBX = 1,NBX
         IRDOFF = (IBX-1)*KZYVAR + 1
         WIBX   = -EIVAL(IBOFF+IBX)
         CALL DZERO(RD(IRDOFF),KZYVAR)
         IF (WIBX .EQ. D0) GOTO 49
         IF (KZCONF.GT.0) THEN
            CALL RSPSLI(1,0,RD(IRDOFF+KZYVAR+KZVAR),DUMMY,
     *                  UDV,RD(IRDOFF),XINDX,RD(KWRK1),LWRK1)
            CALL DSWAP(KZVAR,RD(IRDOFF),1,RD(IRDOFF+KZVAR),1)
            CALL DSCAL(KZYVAR,DM1,RD(IRDOFF),1)
            CALL RSPSLI(1,0,RD(IRDOFF+KZYVAR),DUMMY,
     *                  UDV,RD(IRDOFF),XINDX,RD(KWRK1),LWRK1)
         END IF
         IF (KZWOPT.GT.0) THEN
            IF (SOPPA) THEN
               DO 312 II = 0,KZWOPT-1
                  RD(KSOPO+II) =
     *                 RD(IRDOFF+KZYVAR+KZCONF+II)
                  RD(KSOPO+KZWOPT+II) =
     *                 RD(IRDOFF+KZYVAR+KZVAR+KZCONF+II)
 312           CONTINUE
            ELSE
               KSOPO = IRDOFF + KZYVAR + KZCONF
               DO 313 II = 0,KZWOPT-1
                  RD(IRDOFF+KZYVAR+KZVAR+II) =
     *                 RD(IRDOFF+KZYVAR+KZVAR+KZCONF+II)
 313           CONTINUE
            ENDIF
            CALL RSPSLI(0,1,DUMMY,RD(KSOPO),
     *                  UDV,RD(IRDOFF),XINDX,RD(KWRK1),LWRK1)
         END IF
         CALL DSCAL(KZYVAR,WIBX,RD(IRDOFF),1)
 49      CONTINUE
         IF (SOPPA) THEN
            DO I=0,KZCONF-1
               RD(IRDOFF+I) = RD(IRDOFF+I) +
     *                        RD(KDIAE+I)*RD(IRDOFF+KZYVAR+I)
               RD(IRDOFF+KZVAR+I) = RD(IRDOFF+KZVAR+I) +
     *                        RD(KDIAE+I)*RD(IRDOFF+KZYVAR+KZVAR+I)
            ENDDO
         ENDIF
 50   CONTINUE
C
C ADD  E[2]*X(I)
C
      KRD   = 1
      KWRK1 = KRD + NBX*KZYVAR
      CALL REWSPL(LURSP5)
#if !defined (VAR_SPLITFILES)
      IF (KOFFTY.EQ.1) READ(LURSP5)
#else
      JRSP5 = 0
      IF (KOFFTY .EQ. 1)
     &     CALL READST(LURSP5,LRSP5,RD(KWRK1),JRSP5,.FALSE.)
#endif
      DO 600 K = 1,KZRED
         IF (SOPPA) THEN
            IF (IBTYP(KOFFTY+K) .EQ. JBCNDX) THEN
#if defined (VAR_SPLITFILES)
               CALL READST(LURSP5,LRSP5,RD(KWRK1),JRSP5,.TRUE.)
               IF (LRSP5 .NE. KZYWOP) CALL QUIT('Inconsistent record '//
     &              'length on unit LURSP5 in RSPRES')
            ELSE
               CALL READST(LURSP5,LRSP5,RD(KWRK1),JRSP5,.TRUE.)
               IF (LRSP5 .NE. KZYVAR) CALL QUIT('Inconsistent record '//
     &              'length on unit LURSP5 in RSPRES (2)')
#else
               CALL READT(LURSP5,KZYWOP,RD(KWRK1))
            ELSE
               CALL READT(LURSP5,KZYVAR,RD(KWRK1))
#endif
            ENDIF
         ELSE
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP5,LRSP5,RD(KWRK1),JRSP5,.TRUE.)
            IF (LRSP5 .NE. KZYVAR) CALL QUIT('Inconsistent record '//
     &              'length on unit LURSP5 in RSPRES (3)')
#else
            CALL READT(LURSP5,KZYVAR,RD(KWRK1))
#endif
         ENDIF
         DO 700 JR = 1,NBX
            JRDOFF  = (JR-1)*KZYVAR + 1
            JROOTJ = IBOFF+JR
            FAC1   = EIVEC(2*K-1,JROOTJ)
            FAC2   = EIVEC(2*K,JROOTJ)
            IF (SOPPA .AND. IBTYP(K+KOFFTY) .EQ. JBCNDX) THEN
                  CALL DAXPY(KZWOPT,FAC1,RD(KWRK1),1,
     *                       RD(JRDOFF+KZCONF),1)
                  CALL DAXPY(KZWOPT,FAC1,RD(KWRK1+KZWOPT),1,
     *                       RD(JRDOFF+KZVAR+KZCONF),1)
                  CALL DAXPY(KZWOPT,FAC2,RD(KWRK1),1,
     *                       RD(JRDOFF+KZVAR+KZCONF),1)
                  CALL DAXPY(KZWOPT,FAC2,RD(KWRK1+KZWOPT),1,
     *                       RD(JRDOFF+KZCONF),1)
            ELSE
               CALL DAXPY(KZYVAR,FAC1,RD(KWRK1),1,RD(JRDOFF),1)
               CALL DAXPY(KZVAR,FAC2,RD(KWRK1+KZVAR),1,
     *                    RD(JRDOFF),1)
               CALL DAXPY(KZVAR,FAC2,RD(KWRK1),1,
     *                    RD(KZVAR+JRDOFF),1)
            ENDIF
 700     CONTINUE
 600  CONTINUE
C
C      ADD -GD
C
      DO 1050 JR=1,NBX
         JRDOFF = (JR-1)*KZYVAR + 1
         CALL DAXPY(KZYVAR,DM1,GD(1),1,RD(JRDOFF),1)
 1050 CONTINUE
      IF (IPRRSP.GT.110) THEN
         WRITE (LUPRI,*) ' KEXSIM RESIDUAL VECTORS',KEXSIM
         CALL OUTPUT(RD,1,KZYVAR,1,NBX,KZYVAR,NBX,1,LUPRI)
      END IF
C
C        THE RESIDUAL IS NOW CONSTRUCTED
C
      RETURN
C
C     END OF RSPRES
C
      END
C  /* Deck rspeve */
      SUBROUTINE RSPEVE(IBTYP,EIVAL,EIVEC,BVECS,WRK,NBX,IBOFF)
C
C PURPOSE:
C   CALCULATION OF EIGENVECTORS FROM REDUCED SPACE
C
#include <implicit.h>
#include <dummy.h>
      DIMENSION BVECS(KZYVAR,*),WRK(*)
      DIMENSION IBTYP(*),EIVAL(*),EIVEC(KZYRED,*)
C
#include <ibndxdef.h>
C
#include <priunit.h>
#include <infrsp.h>
#include <inftap.h>
#include <wrkrsp.h>
#include <infpri.h>
C
      NTOT= NBX*KZYVAR
      CALL DZERO(BVECS,NTOT)
      CALL REWSPL(LURSP3)
#if !defined (VAR_SPLITFILES)
      IF (KOFFTY.EQ.1) READ(LURSP3)
#else
      JRSP3 = 0
      IF (KOFFTY.EQ.1) CALL READST(LURSP3,LRSP3,DUMMY,JRSP3,.FALSE.)
#endif
      DO 200 K = 1,KZRED
         IBTYPK = IBTYP(K+KOFFTY)
         IF (IBTYPK.EQ.JBCNDX) THEN
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP3,LRSP3,WRK,JRSP3,.TRUE.)
            IF (LRSP3 .NE. KZCONF) CALL QUIT('Inconsistent record '//
     &           'length on unit LURSP3 in RSPEVE')
#else
            CALL READT(LURSP3,KZCONF,WRK)
#endif
            DO 300 JR = 1,NBX
               CALL DAXPY(KZCONF,EIVEC((2*K-1),IBOFF+JR),
     *                    WRK(1),1,BVECS(1,JR),1)
               JYOFF = 1 + KZVAR
               CALL DAXPY(KZCONF,EIVEC(2*K,IBOFF+JR),
     *                 WRK(1),1,BVECS(JYOFF,JR),1)
 300        CONTINUE
         ELSE
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP3,LRSP3,WRK,JRSP3,.TRUE.)
            IF (LRSP3 .NE. KZYWOP) CALL QUIT('Inconsistent record '//
     &           'length on unit LURSP3 in RSPEVE')
#else
            CALL READT(LURSP3,KZYWOP,WRK)
#endif
            JZOFF = 1 + KZCONF
            JYOFF = 1 + KZVAR +  KZCONF
            DO 301 JR = 1,NBX
               CALL DAXPY(KZWOPT,EIVEC((2*K-1),IBOFF+JR),
     *                 WRK(1),1,BVECS(JZOFF,JR),1)
               CALL DAXPY(KZWOPT,EIVEC((2*K-1),IBOFF+JR),
     *                 WRK(1+KZWOPT),1,BVECS(JYOFF,JR),1)
               CALL DAXPY(KZWOPT,EIVEC(2*K,IBOFF+JR),
     *                 WRK(1+KZWOPT),1,BVECS(JZOFF,JR),1)
               CALL DAXPY(KZWOPT,EIVEC(2*K,IBOFF+JR),
     *                 WRK(1),1,BVECS(JYOFF,JR),1)
 301        CONTINUE
         ENDIF
 200  CONTINUE
C
      IF (IPRRSP .GE. 35) THEN
         DO 400 JR=1,NBX
            WRITE(LUPRI,'(//A,I5/T11,A,T24,A,T44,A/)')
     *         ' EIGENVECTOR NO.',(IBOFF+JR),
     *         'TRIAL VEC.','coeff. for ZY','coeff. for YZ'
            WRITE(LUPRI,'(I15,1P,2D20.6)')
     *         (K,EIVEC(2*K-1,IBOFF+JR),EIVEC(2*K,IBOFF+JR), K=1,KZRED)
 400     CONTINUE
      END IF
C
C     END OF RSPEVE
C
      RETURN
      END
C  /* Deck rspes2 */
      SUBROUTINE RSPES2(GD,REDGD,REDE,REDS,EIVAL,EIVEC,
     *                  CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,
     *                  WRK,LWRK)
C
#include <implicit.h>
C
C  THIS IS A TEST ROUTINE THAT
C  1)  CALCULATE THE E(2) AND S(2) MATRICES EXPLICITLY BY CARRYING
C      OUT LINEAR TRANSFORMATIONS ON UNIT VECTORS
C  2)  DIAGONALIZE THE E(2) AND S(2)
C
#include <priunit.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infpri.h>
C
      PARAMETER ( D0 = 0.0D0 , D1 = 1.0D0 , DTOL = 1.0D-8 )
      DIMENSION GD(*),REDGD(*),REDE(*),REDS(*),EIVAL(*),EIVEC(*)
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
C  ALLOCATE WORK SPACE
C
      KE2   = 1
      KS2   = KE2    + KZYVAR*KZYVAR
      KBVEC = KS2    + KZYVAR*KZYVAR
      KWRK1 = KBVEC  + KZYVAR
      IF (KZCONF.GT.0 .AND. KOFFTY.EQ.1) THEN
         KCREF = KWRK1
         KWRK1 = KCREF  + KZCONF
      ELSE
         KCREF = -999 999 999
      END IF
      LWRK1 = LWRK   - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('RSPES2',KWRK1-1,LWRK)
C
      IF (KZCONF.GT.0 .AND. KOFFTY.EQ.1) THEN
         CALL GETREF(WRK(KCREF),KZCONF)
      END IF
      CALL DZERO(WRK(KBVEC),KZYVAR)
      DO 100 I = 1,KZYVAR
         IF (I.LE.KZCONF) THEN
            NCSIM = 1
            NOSIM = 0
            IOFF  = I
         ELSE IF (I.LE.KZVAR) THEN
            NCSIM = 0
            NOSIM = 1
            IOFF  = I - KZCONF
         ELSE IF (I.LE.KZVAR+KZCONF) THEN
            GO TO 100
         ELSE
            NCSIM = 0
            NOSIM = 1
            IOFF  = I - KZVAR - KZCONF + KZWOPT
         ENDIF
         WRK(KBVEC-1+IOFF) = D1
         IF (( NCSIM.GT.0 ).AND.( KOFFTY.EQ.1)) THEN
            IBOFF = 0
            ICREF  = IOFF
            IF (IOFF.GT.KZCONF) THEN
               IBOFF = KZCONF
               ICREF = ICREF - KZCONF
            END IF
            CALL DAXPY(KZCONF,-WRK(KCREF+ICREF-1),WRK(KCREF),1,
     *                 WRK(KBVEC+IBOFF),1)
         END IF
         IF (IPRRSP.GT.110) THEN
            IF(NOSIM.GT.0)  THEN
               KZYDIM = KZYWOP
               WRITE(LUPRI,'(/A)')' ORBITAL TRIAL VECTOR'
            END IF
            IF(NCSIM.GT.0) THEN
               KZYDIM = KZCONF
               WRITE(LUPRI,'(/A)')' CONFIGURATION TRIAL VECTOR'
            END IF
            CALL OUTPUT(WRK(KBVEC),1,KZYDIM,1,1,KZYDIM,1,1,LUPRI)
         END IF
         CALL RSPLIN(NCSIM,NOSIM,WRK(KBVEC),WRK(KBVEC),
     *               CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *               XINDX,WRK(KWRK1),LWRK1)
C
C        CALL RSPLIN(NCSIM,NOSIM,ZYCVEC,ZYOVEC,
C    *               CMO,UDV,PV,FC,FV,FCAC,H2AC,
C    *               XINDX,WRK,LWRK)
C
C PROJECT OUT RERERENCE STATE COMPONENTS FROM LINEAR TRANSFORMED
C E2 AND S2 VECTORS
C
         IF ((.NOT.TDHF ).AND.( KSYMOP.EQ.1)) THEN
            E2OVL = DDOT(KZCONF,WRK(KCREF),1,WRK(KWRK1),1)
            CALL DAXPY(KZCONF,-E2OVL,WRK(KCREF),1,WRK(KWRK1),1)
            S2OVL = DDOT(KZCONF,WRK(KCREF),1,WRK(KWRK1+KZYVAR),1)
            CALL DAXPY(KZCONF,-S2OVL,WRK(KCREF),1,WRK(KWRK1+KZYVAR),1)
            E2OVL = DDOT(KZCONF,WRK(KCREF),1,WRK(KWRK1+KZVAR),1)
            CALL DAXPY(KZCONF,-E2OVL,WRK(KCREF),1,WRK(KWRK1+KZVAR),1)
            S2OVL = DDOT(KZCONF,WRK(KCREF),1,
     *                   WRK(KWRK1+KZYVAR+KZVAR),1)
            CALL DAXPY(KZCONF,-S2OVL,WRK(KCREF),1,
     *                 WRK(KWRK1+KZYVAR+KZVAR),1)
         END IF
         CALL DCOPY(KZYVAR,WRK(KWRK1),1,WRK(KE2+(I-1)*KZYVAR),1)
         CALL DCOPY(KZYVAR,WRK(KWRK1+KZYVAR),1,
     *              WRK(KS2+(I-1)*KZYVAR),1)
         IF (( NCSIM.GT.0 ).AND.( KOFFTY.EQ.1)) THEN
            CALL DZERO(WRK(KBVEC),KZCONF)
         ELSE
            WRK(KBVEC-1+IOFF) = D0
         END IF
 100  CONTINUE
      KTOT = KZVAR * KZYVAR
      DO 105 J = 1,KZCONF
         JTOT = (J-1)*KZYVAR
         CALL DCOPY(KZVAR,WRK(KE2+JTOT),1,WRK(KE2+JTOT+KTOT+KZVAR),1)
         CALL DCOPY(KZVAR,WRK(KE2+JTOT+KZVAR),1,WRK(KE2+JTOT+KTOT),1)
         CALL DCOPY(KZVAR,WRK(KS2+JTOT),1,WRK(KS2+JTOT+KTOT+KZVAR),1)
         CALL DCOPY(KZVAR,WRK(KS2+JTOT+KZVAR),1,WRK(KS2+JTOT+KTOT),1)
         CALL DSCAL(KZYVAR,-1.0D0,WRK(KS2+JTOT+KTOT),1)
 105  CONTINUE
      WRITE(LUPRI,'(A,I8)')' E(2) MATRIX : DIMENSION ',KZYVAR
      CALL OUTPUT(WRK(KE2),1,KZYVAR,1,KZYVAR,KZYVAR,KZYVAR,1,LUPRI)
      WRITE(LUPRI,'(A,I8)')' S(2) MATRIX : DIMENSION ',KZYVAR
      CALL OUTPUT(WRK(KS2),1,KZYVAR,1,KZYVAR,KZYVAR,KZYVAR,1,LUPRI)
      ZAMAX = D0
      IA    = 0
      JA    = 0
      ZBMAX = D0
      IB    = 0
      JB    = 0
      ZSMAX = D0
      IS    = 0
      JS    = 0
      ZDMAX = D0
      ID    = 0
      JD    = 0
      DO 150 IZ = 1,KZVAR
         IY = IZ +KZVAR
         DO 160 JZ = 1,IZ
            JY = JZ + KZVAR
            ZAIJ  = WRK(KE2-1+(IZ-1)*KZYVAR+JZ)
            ZAJI  = WRK(KE2-1+(JZ-1)*KZYVAR+IZ)
            YAIJ  = WRK(KE2-1+(IY-1)*KZYVAR+JY)
            YAJI  = WRK(KE2-1+(JY-1)*KZYVAR+IY)
            ZDEV  = MAX(ABS(ZAIJ-ZAJI),ABS(ZAIJ-YAIJ))
            ZDEV  = MAX(ZDEV,ABS(ZAIJ-YAJI))
            IF (ZDEV.GT.ZAMAX) THEN
               ZAMAX = ZDEV
               IA    = IZ
               JA    = JZ
            END IF
            ZBIJ  = WRK(KE2-1+(IZ-1)*KZYVAR+JY)
            ZBJI  = WRK(KE2-1+(JZ-1)*KZYVAR+IY)
            YBIJ  = WRK(KE2-1+(IY-1)*KZYVAR+JZ)
            YBJI  = WRK(KE2-1+(JY-1)*KZYVAR+IZ)
            ZDEV  = MAX(ABS(ZBIJ-ZBJI),ABS(ZBIJ-YBIJ))
            ZDEV  = MAX(ZDEV,ABS(ZBIJ-YBJI))
            IF (ZDEV.GT.ZBMAX) THEN
               ZBMAX = ZDEV
               IB    = IZ
               JB    = JZ
            END IF
            ZSIJ  = WRK(KS2-1+(IZ-1)*KZYVAR+JZ)
            ZSJI  = WRK(KS2-1+(JZ-1)*KZYVAR+IZ)
            YSIJ  =-WRK(KS2-1+(IY-1)*KZYVAR+JY)
            YSJI  =-WRK(KS2-1+(JY-1)*KZYVAR+IY)
            ZDEV  = MAX(ABS(ZSIJ-ZSJI),ABS(ZSIJ-YSIJ))
            ZDEV  = MAX(ZDEV,ABS(ZSIJ-YSJI))
            IF (ZDEV.GT.ZSMAX) THEN
               ZSMAX = ZDEV
               IS    = IZ
               JS    = JZ
            END IF
            ZDIJ  = WRK(KS2-1+(IZ-1)*KZYVAR+JY)
            ZDJI  =-WRK(KS2-1+(JZ-1)*KZYVAR+IY)
            YDIJ  =-WRK(KS2-1+(IY-1)*KZYVAR+JZ)
            YDJI  = WRK(KS2-1+(JY-1)*KZYVAR+IZ)
            ZDEV  = MAX(ABS(ZDIJ-ZDJI),ABS(ZDIJ-YDIJ))
            ZDEV  = MAX(ZDEV,ABS(ZDIJ-YDJI))
            IF (ZDEV.GT.ZDMAX) THEN
               ZDMAX = ZDEV
               ID    = IZ
               JD    = JZ
            END IF
 160     CONTINUE
 150  CONTINUE
      IZ=IA
      JZ=JA
      IY=IZ+KZVAR
      JY=JZ+KZVAR
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' A(I,J)-MATRIX : MAXIMUM DEVIATION, I=',IA,' J=',JA,' DEV=',ZAMAX
      IF ( ZAMAX.GT.DTOL)
     *WRITE(LUPRI,'(/4(/A,I5,A,I5,A,1P,G16.8))')
     *' IZ=',IZ,' JZ=',JZ,' E2(IZ,JZ)',WRK(KE2-1+(IZ-1)*KZYVAR+JZ),
     *' IZ=',JZ,' JZ=',IZ,' E2(IZ,JZ)',WRK(KE2-1+(JZ-1)*KZYVAR+IZ),
     *' IZ=',IY,' JZ=',JY,' E2(IZ,JZ)',WRK(KE2-1+(IY-1)*KZYVAR+JY),
     *' IZ=',JY,' JZ=',IY,' E2(IZ,JZ)',WRK(KE2-1+(JY-1)*KZYVAR+IY)
      IZ=IB
      JZ=JB
      IY=IZ+KZVAR
      JY=JZ+KZVAR
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' B(I,J)-MATRIX : MAXIMUM DEVIATION, I=',IB,' J=',JB,' DEV=',ZBMAX
      IF (ZBMAX.GT.DTOL)
     *WRITE(LUPRI,'(/4(/A,I5,A,I5,A,1P,G16.8))')
     *' IZ=',IZ,' JZ=',JY,' E2(IZ,JZ)',WRK(KE2-1+(IZ-1)*KZYVAR+JY),
     *' IZ=',JZ,' JZ=',IY,' E2(IZ,JZ)',WRK(KE2-1+(JZ-1)*KZYVAR+IY),
     *' IZ=',IY,' JZ=',JZ,' E2(IZ,JZ)',WRK(KE2-1+(IY-1)*KZYVAR+JZ),
     *' IZ=',JY,' JZ=',IZ,' E2(IZ,JZ)',WRK(KE2-1+(JY-1)*KZYVAR+IZ)
      IZ=IS
      JZ=JS
      IY=IZ+KZVAR
      JY=JZ+KZVAR
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' S(I,J)-MATRIX : MAXIMUM DEVIATION, I=',JS,' J=',IS,' DEV=',ZSMAX
      IF(ZSMAX.GT.DTOL)
     *WRITE(LUPRI,'(/4(/A,I5,A,I5,A,1P,G16.8))')
     *' IZ=',IZ,' JZ=',JZ,' E2(IZ,JZ)',WRK(KE2-1+(IZ-1)*KZYVAR+JZ),
     *' IZ=',JZ,' JZ=',IZ,' E2(IZ,JZ)',WRK(KE2-1+(JZ-1)*KZYVAR+IZ),
     *' IZ=',IY,' JZ=',JY,' E2(IZ,JZ)',WRK(KE2-1+(IY-1)*KZYVAR+JY),
     *' IZ=',JY,' JZ=',IY,' E2(IZ,JZ)',WRK(KE2-1+(JY-1)*KZYVAR+IY)
      IZ=ID
      JZ=JD
      IY=IZ+KZVAR
      JY=JZ+KZVAR
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' D(I,J)-MATRIX : MAXIMUM DEVIATION, I=',ID,' J=',JD,' DEV=',ZDMAX
      IF(ZDMAX.GT.DTOL)
     *WRITE(LUPRI,'(/4(/A,I5,A,I5,A,1P,G16.8))')
     *' IZ=',IZ,' JZ=',JY,' E2(IZ,JZ)',WRK(KE2-1+(IZ-1)*KZYVAR+JY),
     *' IZ=',JZ,' JZ=',IY,' E2(IZ,JZ)',WRK(KE2-1+(JZ-1)*KZYVAR+IY),
     *' IZ=',IY,' JZ=',JZ,' E2(IZ,JZ)',WRK(KE2-1+(IY-1)*KZYVAR+JZ),
     *' IZ=',JY,' JZ=',IZ,' E2(IZ,JZ)',WRK(KE2-1+(JY-1)*KZYVAR+IZ)
      IF ( KZYVAR.GT.MAXRM ) RETURN
      KZYRED = KZYVAR
      KZRED  = KZVAR
      IJ = 0
      DO 200 I = 1,KZYVAR
         DO 300 J = 1,I
            IJ = IJ + 1
            REDE(IJ) = WRK(KE2-1+(I-1)*KZYVAR+J)
            REDS(IJ) = WRK(KS2-1+(I-1)*KZYVAR+J)
 300     CONTINUE
 200  CONTINUE
C
C     CALL RSPRED(,2,..) SOLVE FOR EIGENVALUES AND EIGENVECTORS
C
      ICTL = 2
      NSIM = 0
      CALL RSPRED(ICTL,.FALSE.,NSIM,IBTYP,GD,REDGD,REDE,REDS,
     *            EIVAL,EIVEC,WRK(KWRK1),WRK(KWRK1),
     *            UDV,WRK(KWRK1),XINDX,WRK(KWRK1),LWRK1)
C     CALL RSPRED (ICTL,LINEQ,N,IBTYP,GD,REDGD,REDE,REDS,
C    *            EIVAL,EIVEC,BCVEC,BOVEC,UDV,EVECS,XINDX,
C    *            WRK,LWRK)
C
C CHECK EIGENVALUES OF E2
C
      WRITE(LUPRI,'(/A,/A)')' ********************************',
     *' CHECK EIGENVALUES OF E2 WITH RSPRED ROUTINE '
      IJ = 0
      DO 500 J = 1,KZYVAR
         DO 600 I = 1,J
            IJ = IJ + 1
            IF (I.EQ.J) THEN
               REDS(IJ) = D1
            ELSE
               REDS(IJ) = D0
            END IF
 600     CONTINUE
 500  CONTINUE
      CALL RSPRED(ICTL,.FALSE.,NSIM,IBTYP,GD,REDGD,REDE,REDS,
     *            EIVAL,EIVEC,WRK(KWRK1),WRK(KWRK1),
     *            UDV,WRK(KWRK1),XINDX,WRK(KWRK1),LWRK1)
C
C CHECK EIGENVALUES OF S2
C
      WRITE(LUPRI,'(/A,/A)')' ********************************',
     *' CHECK EIGENVALUES OF S2 WITH RSPRED ROUTINE '
      IJ = 0
      DO 700 I = 1,KZYVAR
         DO 800 J = 1,I
            IJ = IJ + 1
            REDE(IJ) = WRK(KS2-1+(I-1)*KZYVAR+J)
            IF (I.EQ.J) THEN
               REDS(IJ) = D1
            ELSE
               REDS(IJ) = D0
            END IF
 800     CONTINUE
 700  CONTINUE
      CALL RSPRED(ICTL,.FALSE.,NSIM,IBTYP,GD,REDGD,REDE,REDS,
     *            EIVAL,EIVEC,WRK(KWRK1),WRK(KWRK1),
     *            UDV,WRK(KWRK1),XINDX,WRK(KWRK1),LWRK1)
      RETURN
      END
C  /* Deck oropti */
      SUBROUTINE OROPTI(EIVAL,IBTYP,A1,CMO,UDV,PV,FC,FV,FCAC,XINDX,
     *                  WRK,LWRK)
C
C PURPOSE SOLVE ORBITAL PART OF LINEAR SET OF EQUATIONS
C
#include <implicit.h>
C
      DIMENSION EIVAL(*),IBTYP(*),A1(*)
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FCAC(*),XINDX(*),WRK(*)
C
      PARAMETER ( D0 = 0.0D0 , D1 = 1.0D0 , DM1 = -1.0D0 )
C
#include <infrsp.h>
#include <wrkrsp.h>
C
C PUT DIAGONAL ORBITAL PART OF E(2) AND S(2) IN WRK(1) AND
C AND WRK(1+KZYCON)
C
      CALL RSPDIA(WRK,D0)
      CALL RSPDIA(WRK(1+KZYVAR),D1)
      DO 313 II = 1,KZWOPT
         WRK(II) = WRK(KZCONF+II)
 313  CONTINUE
C     CALL DCOPY(KZWOPT,WRK(1+KZCONF),1,WRK,1)
      DO 314 II = 1,KZWOPT
         WRK(KZWOPT+II) = WRK(KZVAR+KZCONF+II)
 314  CONTINUE
C     CALL DCOPY(KZWOPT,WRK(1+KZVAR+KZCONF),1,WRK(1+KZWOPT),1)
      DO 315 II = 1,KZWOPT
         WRK(KZYWOP+II) = WRK(KZYVAR+KZCONF+II)
 315  CONTINUE
C     CALL DCOPY(KZWOPT,WRK(1+KZYVAR+KZCONF),1,WRK(1+KZYWOP),1)
      DO 316 II = 1,KZWOPT
         WRK(KZYWOP+KZWOPT+II) = WRK(KZYVAR+KZVAR+KZCONF+II)
 316  CONTINUE
C     CALL DCOPY(KZWOPT,WRK(1+KZYVAR+KZVAR+KZCONF),1,
C    *                                  WRK(1+KZYWOP+KZWOPT),1)
      CALL DAXPY(KZYWOP,DM1,WRK,1,WRK(1+KZYWOP),1)
      LWRKPA = LWRK - 2*KZYWOP - 1
      IF (LWRKPA.LT.0) CALL ERRWRK('OROPTI',2*KZYWOP,LWRK)
      NOSIM  = 1
      CALL ORPPAR(NOSIM,THCRSP,EIVAL,IBTYP,A1,WRK,
     *            WRK(1+KZYWOP),CMO,UDV,PV,FC,FV,FCAC,XINDX,
     *            WRK(1+2*KZYWOP),LWRKPA)
C
C     CALL ORPPAR(NOSIM,THCORP,EOVAL,IBTYP,A1,ORBDIE,
C    *            ORBDIS,CMO,UDV,PVX,FC,FV,FCAC,XINDX,
C    *            WRK,LWRK)
C
      RETURN
      END
C  /* Deck e2osym */
      SUBROUTINE E2OSYM(KOZRED,OBVEC,OEVEC,OSVEC,
     *           CMO,UDV,PV,FC,FV,FCAC,XINDX,WRK,LWRK)
C
#include <implicit.h>
C
C  THIS IS A TEST ROUTINE TO CHECK THE ORBITAL PART OF THE
C  RESPONSE MATRICES
C
C  1)  SET UP ORBITAL PART OF REDUCED E(2) AND S(2)
C      AS SQUARE MATRICES AND CHECK THE STRUCTURE
C
C  2)  CHECK ORTHOGONALITY OF TRIAL VECTORS
C
#include <priunit.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infpri.h>
C
      PARAMETER ( D0 = 0.0D0 , D1 = 1.0D0 , DTOL = 1.0D-6 )
      DIMENSION OBVEC(KZYWOP,*),OEVEC(KZYWOP,*),OSVEC(KZYWOP,*)
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*)
      DIMENSION XINDX(*),WRK(*)
C
C
      KZYDIM = 2 * KOZRED
      KZDIM  = KOZRED
      DO 15 J = 1,KZDIM
         CALL DCOPY(KZWOPT,OBVEC(1,J),1,OBVEC(1+KZWOPT,J+KZDIM),1)
         CALL DCOPY(KZWOPT,OBVEC(1+KZWOPT,J),1,OBVEC(1,J+KZDIM),1)
 15   CONTINUE
      IF (IPRRSP.GT.110) THEN
         WRITE(LUPRI,'(2A,I8)')
     *   ' ORBITAL E(2) LINEAR TRANSFORMED VECTORS',
     *   ' : DIMENSION ',KZDIM
         CALL OUTPUT(OEVEC,1,KZYWOP,1,KZDIM,KZYWOP,KZDIM,1,LUPRI)
         WRITE(LUPRI,'(2A,I8)')
     *   ' ORBITAL S(2) LINEAR TRANSFORMED VECTORS',
     *   ' : DIMENSION ',KZDIM
         CALL OUTPUT(OSVEC,1,KZYWOP,1,KZDIM,KZYWOP,KZDIM,1,LUPRI)
            WRITE(LUPRI,'(2A,I8)')' ORBITAL TRIAL VECTORS ',
     *      ' : DIMENSION ',2*KZDIM
         CALL OUTPUT(OBVEC,1,KZYWOP,1,KZYDIM,KZYWOP,KZYDIM,1,LUPRI)
      END IF
C
C CARRY OUT LINEAR TRANSFORMATION ON (Y,Z) VECTORS
C
      NSIM = 1
      DO 20 J = 1,KZDIM
         CALL ORPLIN(NSIM,OBVEC(1,KZDIM+J),CMO,UDV,PV,FC,FV,FCAC,
     *               XINDX,WRK,LWRK)
C        CALL ORPLIN(NOSIM,ZYOVEC,CMO,UDV,PVX,FC,FV,FCAC,
C    *               XINDX,WRK,LWRK)
         CALL DCOPY(KZYWOP,WRK(1),1,OEVEC(1,KZDIM+J),1)
         CALL DCOPY(KZYWOP,WRK(1+KZYWOP),1,OSVEC(1,KZDIM+J),1)
C
C        CALL ORPLIN(KZDIM,OBVEC,OEVEC,OSVEC,WRK,LWRK1)
C
 20   CONTINUE
      SMAX = D0
      IS = 0
      JS = 0
      KS = 0
      EMAX = D0
      IE = 0
      JE = 0
      KE = 0
      DO 25 I = 1,KZDIM
         DO 26 J = 1,KZYWOP
            IF (J.GT.KZWOPT) THEN
               KADD = -KZWOPT
            ELSE
               KADD = KZWOPT
            END IF
            XE =ABS(ABS(OEVEC(J,I))-ABS(OEVEC(J+KADD,KZDIM+I)))
            IF (XE.GT.EMAX) THEN
               IE = I
               JE = J
               KE = J+KADD
               EMAX = XE
            END IF
            XS =ABS(ABS(OSVEC(J,I))-ABS(OSVEC(J+KADD,KZDIM+I)))
            IF (XS.GT.SMAX) THEN
               IS = I
               JS = J
               KS = J + KADD
               SMAX = XS
            END IF
 26      CONTINUE
 25   CONTINUE
      WRITE(LUPRI,'(/A)')
     *' LARGEST DEVIATIONS IN LINEAR TRANSFORMED VECTORS'
         WRITE(LUPRI,'(/A)')' E(2) ORBITAL PART'
         WRITE(LUPRI,'(/A,1P,G16.8)')
     *   ' LARGEST DEVIATION  :',EMAX
      IF ( EMAX.GT.DTOL ) THEN
         WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *   ' OEVEC(J,I) J:',JE,'   I:',IE,' VALUE:',OEVEC(JE,IE)
         WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *   ' OEVEC(J,I) J:',KE,'   I:',IE+KZDIM,' VALUE:',
     *     OEVEC(KE,IE+KZDIM)
      END IF
         WRITE(LUPRI,'(/A)')' S(2) ORBITAL PART'
         WRITE(LUPRI,'(/A,1P,G16.8)')
     *   ' LARGEST DEVIATION  :',SMAX
      IF ( SMAX.GT.DTOL ) THEN
         WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *   ' OEVEC(J,I) J:',JS,'   I:',IS,' VALUE:',OSVEC(JS,IS)
         WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *   ' OEVEC(J,I) J:',KS,'   I:',IS+KZDIM,' VALUE:',
     *     OEVEC(KS,IS+KZDIM)
      END IF
      IF (IPRRSP.GT.130) THEN
         WRITE(LUPRI,'(2A,I8)')
     *   ' ORBITAL E(2) LINEAR TRANSFORMED VECTORS',
     *   ' : DIMENSION ',KZYDIM
         CALL OUTPUT(OEVEC,1,KZYWOP,1,KZYDIM,KZYWOP,KZYDIM,1,LUPRI)
         WRITE(LUPRI,'(2A,I8)')
     *   ' ORBITAL S(2) LINEAR TRANSFORMED VECTORS',
     *   ' : DIMENSION ',KZYDIM
         CALL OUTPUT(OSVEC,1,KZYWOP,1,KZYDIM,KZYWOP,KZYDIM,1,LUPRI)
            WRITE(LUPRI,'(2A,I8)')' ORBITAL TRIAL VECTORS ',
     *      ' : DIMENSION ',2*KZDIM
         CALL OUTPUT(OBVEC,1,KZYWOP,1,KZYDIM,KZYWOP,KZYDIM,1,LUPRI)
      END IF
C
C CHECK IF TRIAL VECTORS ARE ORTHONORMAL
C
      KE2MAT = 1
      KS2MAT = KE2MAT + KZYDIM*KZYDIM
      KWRK1  = KS2MAT + KZYDIM*KZYDIM
      XMAX = 0.0D0
      IMAX = 0
      JMAX = 0
      IJ   = 0
      DO 30 I = 1,2*KZDIM
         DO 40 J = 1,I
            IJ = IJ + 1
            XIJ = DDOT(KZYWOP,OBVEC(1,I),1,OBVEC(1,J),1)
            IF((I.NE.J).AND.(ABS(XIJ).GT.XMAX)) THEN
               IMAX = I
               JMAX = J
               XMAX = XIJ
            END IF
            WRK(KWRK1-1+IJ) = XIJ
 40      CONTINUE
 30   CONTINUE
            WRITE(LUPRI,'(/A)')
     *      ' OVERLAP FOR ORBITAL TRIAL VECTORS'
         IF (IPRRSP.GT.90) THEN
            CALL OUTPAK(WRK(KWRK1),2*KZDIM,1,LUPRI)
         END IF
         IF ( XMAX.GT.DTOL )
     *   WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *   ' LARGEST ELEMENT (I,J), I=',IMAX,' J=',JMAX,' :',XMAX
C
C SET UP REDUCED E2 AND S2 MATRICES
C
      DO 70 J = 1,KZYDIM
         DO 80 I = 1,KZYDIM
            X1 = DDOT(KZYWOP,OBVEC(1,J),1,OEVEC(1,I),1)
            WRK(KE2MAT-1+(I-1)*KZYDIM+J) = X1
            X1 = DDOT(KZYWOP,OBVEC(1,J),1,OSVEC(1,I),1)
            WRK(KS2MAT-1+(I-1)*KZYDIM+J) = X1
 80      CONTINUE
 70   CONTINUE
      WRITE(LUPRI,'(A,I8)')' REDUCED E(2) MATRIX : DIMENSION ',KZYDIM
      CALL OUTPUT(WRK(KE2MAT),1,KZYDIM,1,KZYDIM,KZYDIM,KZYDIM,1,LUPRI)
      WRITE(LUPRI,'(A,I8)')' REDUCED S(2) MATRIX : DIMENSION ',KZYDIM
      CALL OUTPUT(WRK(KS2MAT),1,KZYDIM,1,KZYDIM,KZYDIM,KZYDIM,1,LUPRI)
C
C CHECK BLOCK STRUCTURE OF E2 AND S2 AND WRITE OUT LARGEST DEVIATIONS
C
      ZAMAX = D0
      IA    = 0
      JA    = 0
      ZBMAX = D0
      IB    = 0
      JB    = 0
      ZSMAX = D0
      IS    = 0
      JS    = 0
      ZDMAX = D0
      ID    = 0
      JD    = 0
      DO 150 IZ = 1,KZDIM
         IY = IZ +KZDIM
         DO 160 JZ = 1,IZ
            JY = JZ + KZDIM
            ZAIJ  = WRK(KE2MAT-1+(IZ-1)*KZYDIM+JZ)
            ZAJI  = WRK(KE2MAT-1+(JZ-1)*KZYDIM+IZ)
            YAIJ  = WRK(KE2MAT-1+(IY-1)*KZYDIM+JY)
            YAJI  = WRK(KE2MAT-1+(JY-1)*KZYDIM+IY)
            ZDEV  = MAX(ABS(ZAIJ-ZAJI),ABS(ZAIJ-YAIJ))
            ZDEV  = MAX(ZDEV,ABS(ZAIJ-YAJI))
            IF (ZDEV.GT.ZAMAX) THEN
               ZAMAX = ZDEV
               IA    = IZ
               JA    = JZ
            END IF
            ZBIJ  = WRK(KE2MAT-1+(IZ-1)*KZYDIM+JY)
            ZBJI  = WRK(KE2MAT-1+(JZ-1)*KZYDIM+IY)
            YBIJ  = WRK(KE2MAT-1+(IY-1)*KZYDIM+JZ)
            YBJI  = WRK(KE2MAT-1+(JY-1)*KZYDIM+IZ)
            ZDEV  = MAX(ABS(ZBIJ-ZBJI),ABS(ZBIJ-YBIJ))
            ZDEV  = MAX(ZDEV,ABS(ZBIJ-YBJI))
            IF (ZDEV.GT.ZBMAX) THEN
               ZBMAX = ZDEV
               IB    = IZ
               JB    = JZ
            END IF
            ZSIJ  = WRK(KS2MAT-1+(IZ-1)*KZYDIM+JZ)
            ZSJI  = WRK(KS2MAT-1+(JZ-1)*KZYDIM+IZ)
            YSIJ  =-WRK(KS2MAT-1+(IY-1)*KZYDIM+JY)
            YSJI  =-WRK(KS2MAT-1+(JY-1)*KZYDIM+IY)
            ZDEV  = MAX(ABS(ZSIJ-ZSJI),ABS(ZSIJ-YSIJ))
            ZDEV  = MAX(ZDEV,ABS(ZSIJ-YSJI))
            IF (ZDEV.GT.ZSMAX) THEN
               ZSMAX = ZDEV
               IS    = IZ
               JS    = JZ
            END IF
            ZDIJ  = WRK(KS2MAT-1+(IZ-1)*KZYDIM+JY)
            ZDJI  =-WRK(KS2MAT-1+(JZ-1)*KZYDIM+IY)
            YDIJ  =-WRK(KS2MAT-1+(IY-1)*KZYDIM+JZ)
            YDJI  = WRK(KS2MAT-1+(JY-1)*KZYDIM+IZ)
            ZDEV  = MAX(ABS(ZDIJ-ZDJI),ABS(ZDIJ-YDIJ))
            ZDEV  = MAX(ZDEV,ABS(ZDIJ-YDJI))
            IF (ZDEV.GT.ZDMAX) THEN
               ZDMAX = ZDEV
               ID    = IZ
               JD    = JZ
            END IF
 160     CONTINUE
 150  CONTINUE
      IZ=IA
      JZ=JA
      IY=IZ+KZDIM
      JY=JZ+KZDIM
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' A(I,J)-MATRIX : MAXIMUM DEVIATION, I=',IA,' J=',JA,' DEV=',ZAMAX
      IF ( ZAMAX.GT.DTOL )
     *WRITE(LUPRI,'(/4(/A,I5,A,I5,A,1P,G16.8))')
     *' IZ=',IZ,' JZ=',JZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(IZ-1)*KZYDIM+JZ),
     *' IZ=',JZ,' JZ=',IZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(JZ-1)*KZYDIM+IZ),
     *' IZ=',IY,' JZ=',JY,' E2(IZ,JZ)',WRK(KE2MAT-1+(IY-1)*KZYDIM+JY),
     *' IZ=',JY,' JZ=',IY,' E2(IZ,JZ)',WRK(KE2MAT-1+(JY-1)*KZYDIM+IY)
      IZ=IB
      JZ=JB
      IY=IZ+KZDIM
      JY=JZ+KZDIM
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' B(I,J)-MATRIX : MAXIMUM DEVIATION, I=',IB,' J=',JB,' DEV=',ZBMAX
      IF ( ZBMAX.GT.DTOL )
     *WRITE(LUPRI,'(/4(/A,I5,A,I5,A,1P,G16.8))')
     *' IZ=',IZ,' JZ=',JY,' E2(IZ,JZ)',WRK(KE2MAT-1+(IZ-1)*KZYDIM+JY),
     *' IZ=',JZ,' JZ=',IY,' E2(IZ,JZ)',WRK(KE2MAT-1+(JZ-1)*KZYDIM+IY),
     *' IZ=',IY,' JZ=',JZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(IY-1)*KZYDIM+JZ),
     *' IZ=',JY,' JZ=',IZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(JY-1)*KZYDIM+IZ)
      IZ=IS
      JZ=JS
      IY=IZ+KZDIM
      JY=JZ+KZDIM
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' S(I,J)-MATRIX : MAXIMUM DEVIATION, I=',IS,' J=',JS,' DEV=',ZSMAX
      IF ( ZSMAX.GT.DTOL )
     *WRITE(LUPRI,'(/4(/A,I5,A,I5,A,1P,G16.8))')
     *' IZ=',IZ,' JZ=',JY,' E2(IZ,JZ)',WRK(KE2MAT-1+(IZ-1)*KZYDIM+JY),
     *' IZ=',JZ,' JZ=',IY,' E2(IZ,JZ)',WRK(KE2MAT-1+(JZ-1)*KZYDIM+IY),
     *' IZ=',IY,' JZ=',JZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(IY-1)*KZYDIM+JZ),
     *' IZ=',JY,' JZ=',IZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(JY-1)*KZYDIM+IZ)
      IZ=ID
      JZ=JD
      IY=IZ+KZDIM
      JY=JZ+KZDIM
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' D(I,J)-MATRIX : MAXIMUM DEVIATION, I=',IA,' J=',JA,' DEV=',ZDMAX
      IF ( ZDMAX.GT.DTOL )
     *WRITE(LUPRI,'(/4(/A,I5,A,I5,A,1P,G16.8))')
     *' IZ=',IZ,' JZ=',JY,' E2(IZ,JZ)',WRK(KE2MAT-1+(IZ-1)*KZYDIM+JY),
     *' IZ=',JZ,' JZ=',IY,' E2(IZ,JZ)',WRK(KE2MAT-1+(JZ-1)*KZYDIM+IY),
     *' IZ=',IY,' JZ=',JZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(IY-1)*KZYDIM+JZ),
     *' IZ=',JY,' JZ=',IZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(JY-1)*KZYDIM+IZ)
C
C CHECK EIGENVALUES OF E2
C
C ALLOCATE WORK SPACE
C
      KE2    = KWRK1
      KS2    = KE2   + KZYDIM*KZYDIM
      KEIVAL = KS2   + KZYDIM*KZYDIM
      KEIVEC = KEIVAL+ KZYDIM
      KWRK1  = KEIVEC+ KZYDIM*KZYDIM
      LWRK1  = LWRK  - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('E2OSYM',KWRK1-1,LWRK)
C
      WRITE(LUPRI,'(/A,/A)')' ********************************',
     *' CHECK EIGENVALUES OF E2 WITH RSPRED ROUTINE '
      IJ = 0
      DO 500 J = 1,KZYDIM
         DO 600 I = 1,J
            IJ = IJ + 1
            WRK(KE2-1+IJ) = WRK(KE2MAT-1+(I-1)*KZYDIM+J)
            IF (I.EQ.J) THEN
               WRK(KS2-1+IJ) = D1
            ELSE
               WRK(KS2-1+IJ) = D0
            END IF
 600     CONTINUE
 500  CONTINUE
      ICTL = 2
      NSIM = 0
      KZSAV = KZRED
      KZYSAV= KZYRED
      KZRED = KZDIM
      KZYRED= KZYDIM
      CALL RSPRED(ICTL,.FALSE.,NSIM,IBTYP,GD,REDGD,WRK(KE2),WRK(KS2),
     *            WRK(KEIVAL),WRK(KEIVEC),WRK(KWRK1),WRK(KWRK1),
     *            UDV,WRK(KWRK1),XINDX,WRK(KWRK1),LWRK1)
      KZRED = KZSAV
      KZYRED= KZYSAV
C
C CHECK EIGENVALUES OF S2
C
      WRITE(LUPRI,'(/A,/A)')' ********************************',
     *' CHECK EIGENVALUES OF S2 WITH RSPRED ROUTINE '
      IJ = 0
      DO 700 J = 1,KZYDIM
         DO 800 I = 1,J
            IJ = IJ + 1
            WRK(KE2-1+IJ) = WRK(KS2MAT-1+(I-1)*KZYDIM+J)
            IF (I.EQ.J) THEN
               WRK(KS2-1+IJ) = D1
            ELSE
               WRK(KS2-1+IJ) = D0
            END IF
 800     CONTINUE
 700  CONTINUE
      ICTL = 2
      KZSAV = KZRED
      KZYSAV= KZYRED
      KZRED = KZDIM
      KZYRED= KZYDIM
      CALL RSPRED(ICTL,.FALSE.,NSIM,IBTYP,GD,REDGD,WRK(KE2),WRK(KS2),
     *            WRK(KEIVAL),WRK(KEIVEC),WRK(KWRK1),WRK(KWRK1),
     *            UDV,WRK(KWRK1),XINDX,WRK(KWRK1),LWRK1)
      KZRED = KZSAV
      KZYRED= KZYSAV
C
C END OF E2OSYM
C
      RETURN
      END
C  /* Deck e2ochk */
      SUBROUTINE E2OCHK(OBVEC,OEVEC,OSVEC,CMO,UDV,PV,FC,FV,FCAC,XINDX,
     *                  WRK,LWRK)
C
#include <implicit.h>
C
C  THIS IS A TEST ROUTINE TO CHECK THE ORBITAL PART OF THE
C  RESPONSE MATRICES
C
C  1)  SET UP BLOCK OF ORBITAL E(2) AND S(2)
C      AS SQUARE MATRICES AND CHECK THE STRUCTURE.
C      THE BLOCK IS FORMED USING UNIT VECTORS AS TRIAL
C      VECTORS. TOTALLY MAXCHK VECTORS ARE USED WITH ONE IN
C      ELEMENTS ISTCHK TO ISTCHK+MAXCHK-1
C
#include <priunit.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infpri.h>
C
      PARAMETER ( D0 = 0.0D0 , D1 = 1.0D0 , DTOL = 1.0D-6 )
      DIMENSION OBVEC(KZYWOP,*),OEVEC(KZYWOP,*),OSVEC(KZYWOP,*)
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*)
      DIMENSION XINDX(*),WRK(*)
C
C
      KZYDIM = 2 * MAXOCK
      KZDIM  = MAXOCK
      NTOT   = KZDIM * KZYWOP
      CALL DZERO(OBVEC,NTOT)
      DO 10 J = 1,KZDIM
         OBVEC(ISTOCK-1+J,J) = D1
 10   CONTINUE
      DO 15 J = 1,KZDIM
         CALL DCOPY(KZWOPT,OBVEC(1,J),1,OBVEC(1+KZWOPT,J+KZDIM),1)
         CALL DCOPY(KZWOPT,OBVEC(1+KZWOPT,J),1,OBVEC(1,J+KZDIM),1)
 15   CONTINUE
      IF (IPRRSP.GT.130) THEN
         WRITE(LUPRI,'(2A,I8)')
     *   ' ORBITAL E(2) LINEAR TRANSFORMED VECTORS',
     *   ' : DIMENSION ',KZDIM
         CALL OUTPUT(OEVEC,1,KZYWOP,1,KZDIM,KZYWOP,KZDIM,1,LUPRI)
         WRITE(LUPRI,'(2A,I8)')
     *   ' ORBITAL S(2) LINEAR TRANSFORMED VECTORS',
     *   ' : DIMENSION ',KZDIM
         CALL OUTPUT(OSVEC,1,KZYWOP,1,KZDIM,KZYWOP,KZDIM,1,LUPRI)
            WRITE(LUPRI,'(2A,I8)')' ORBITAL TRIAL VECTORS ',
     *      ' : DIMENSION ',2*KZDIM
         CALL OUTPUT(OBVEC,1,KZYWOP,1,KZYDIM,KZYWOP,KZYDIM,1,LUPRI)
      END IF
C
C CARRY OUT LINEAR TRANSFORMATION ON (Y,Z) VECTORS
C
      NSIM = 1
      DO 20 J = 1,KZYDIM
         CALL ORPLIN(NSIM,OBVEC(1,J),CMO,UDV,PV,FC,FV,FCAC,XINDX,
     *               WRK,LWRK)
C        CALL ORPLIN(NOSIM,ZYOVEC,CMO,UDV,PVX,FC,FV,FCAC,XINDX,
C    *               WRK,LWRK)
         CALL DCOPY(KZYWOP,WRK(1),1,OEVEC(1,J),1)
         CALL DCOPY(KZYWOP,WRK(1+KZYWOP),1,OSVEC(1,J),1)
C
C        CALL ORPLIN(KZDIM,OBVEC,OEVEC,OSVEC,WRK,LWRK1)
C
 20   CONTINUE
      SMAX = D0
      IS = 0
      JS = 0
      KS = 0
      EMAX = D0
      IE = 0
      JE = 0
      KE = 0
      DO 25 I = 1,KZDIM
         DO 26 J = 1,KZYWOP
            IF (J.GT.KZWOPT) THEN
               KADD = -KZWOPT
            ELSE
               KADD = KZWOPT
            END IF
            XE =ABS(ABS(OEVEC(J,I))-ABS(OEVEC(J+KADD,KZDIM+I)))
            IF (XE.GT.EMAX) THEN
               IE = I
               JE = J
               KE = J+KADD
               EMAX = XE
            END IF
            XS =ABS(ABS(OSVEC(J,I))-ABS(OSVEC(J+KADD,KZDIM+I)))
            IF (XS.GT.SMAX) THEN
               IS = I
               JS = J
               KS = J + KADD
               SMAX = XS
            END IF
 26      CONTINUE
 25   CONTINUE
      WRITE(LUPRI,'(/A)')
     *' LARGEST DEVIATIONS IN LINEAR TRANSFORMED VECTORS'
         WRITE(LUPRI,'(/A)')' E(2) ORBITAL PART'
         WRITE(LUPRI,'(/A,1P,G16.8)')
     *   ' LARGEST DEVIATION  :',EMAX
      IF ( EMAX.GT.DTOL ) THEN
         WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *   ' OEVEC(J,I) J:',JE,'   I:',IE,' VALUE:',OEVEC(JE,IE)
         WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *   ' OEVEC(J,I) J:',KE,'   I:',IE+KZDIM,' VALUE:',
     *     OEVEC(KE,IE+KZDIM)
      END IF
         WRITE(LUPRI,'(/A)')' S(2) ORBITAL PART'
         WRITE(LUPRI,'(/A,1P,G16.8)')
     *   ' LARGEST DEVIATION  :',SMAX
      IF ( SMAX.GT.DTOL ) THEN
         WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *   ' OEVEC(J,I) J:',JS,'   I:',IS,' VALUE:',OSVEC(JS,IS)
         WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *   ' OEVEC(J,I) J:',KS,'   I:',IS+KZDIM,' VALUE:',
     *     OEVEC(KS,IS+KZDIM)
      END IF
      IF (IPRRSP.GT.130) THEN
         WRITE(LUPRI,'(2A,I8)')
     *   ' ORBITAL E(2) LINEAR TRANSFORMED VECTORS',
     *   ' : DIMENSION ',KZYDIM
         CALL OUTPUT(OEVEC,1,KZYWOP,1,KZYDIM,KZYWOP,KZYDIM,1,LUPRI)
         WRITE(LUPRI,'(2A,I8)')
     *   ' ORBITAL S(2) LINEAR TRANSFORMED VECTORS',
     *   ' : DIMENSION ',KZYDIM
         CALL OUTPUT(OSVEC,1,KZYWOP,1,KZYDIM,KZYWOP,KZYDIM,1,LUPRI)
            WRITE(LUPRI,'(2A,I8)')' ORBITAL TRIAL VECTORS ',
     *      ' : DIMENSION ',2*KZDIM
         CALL OUTPUT(OBVEC,1,KZYWOP,1,KZYDIM,KZYWOP,KZYDIM,1,LUPRI)
      END IF
C
C CHECK IF TRIAL VECTORS ARE ORTHONORMAL
C
      KE2MAT = 1
      KS2MAT = KE2MAT + KZYDIM*KZYDIM
      KWRK1  = KS2MAT + KZYDIM*KZYDIM
      XMAX = 0.0D0
      IMAX = 0
      JMAX = 0
      IJ   = 0
      DO 30 I = 1,2*KZDIM
         DO 40 J = 1,I
            IJ = IJ + 1
            XIJ = DDOT(KZYWOP,OBVEC(1,I),1,OBVEC(1,J),1)
            IF((I.NE.J).AND.(ABS(XIJ).GT.XMAX)) THEN
               IMAX = I
               JMAX = J
               XMAX = XIJ
            END IF
            WRK(KWRK1-1+IJ) = XIJ
 40      CONTINUE
 30   CONTINUE
            WRITE(LUPRI,'(/A)')
     *      ' OVERLAP FOR ORBITAL TRIAL VECTORS'
         IF (IPRRSP.GT.130) THEN
            CALL OUTPAK(WRK(KWRK1),2*KZDIM,1,LUPRI)
         END IF
         IF ( XMAX.GT.DTOL )
     *   WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *   ' LARGEST ELEMENT (I,J), I=',IMAX,' J=',JMAX,' :',XMAX
C
C SET UP REDUCED E2 AND S2 MATRICES
C
      DO 70 J = 1,KZYDIM
         DO 80 I = 1,KZYDIM
            X1 = DDOT(KZYWOP,OBVEC(1,J),1,OEVEC(1,I),1)
            WRK(KE2MAT-1+(I-1)*KZYDIM+J) = X1
            X1 = DDOT(KZYWOP,OBVEC(1,J),1,OSVEC(1,I),1)
            WRK(KS2MAT-1+(I-1)*KZYDIM+J) = X1
 80      CONTINUE
 70   CONTINUE
      WRITE(LUPRI,'(A,I8)')' REDUCED E(2) MATRIX : DIMENSION ',KZYDIM
      CALL OUTPUT(WRK(KE2MAT),1,KZYDIM,1,KZYDIM,KZYDIM,KZYDIM,1,LUPRI)
      WRITE(LUPRI,'(A,I8)')' REDUCED S(2) MATRIX : DIMENSION ',KZYDIM
      CALL OUTPUT(WRK(KS2MAT),1,KZYDIM,1,KZYDIM,KZYDIM,KZYDIM,1,LUPRI)
C
C CHECK BLOCK STRUCTURE OF E2 AND S2 AND WRITE OUT LARGEST DEVIATIONS
C
      ZAMAX = D0
      IA    = 0
      JA    = 0
      ZBMAX = D0
      IB    = 0
      JB    = 0
      ZSMAX = D0
      IS    = 0
      JS    = 0
      ZDMAX = D0
      ID    = 0
      JD    = 0
      DO 150 IZ = 1,KZDIM
         IY = IZ +KZDIM
         DO 160 JZ = 1,IZ
            JY = JZ + KZDIM
            ZAIJ  = WRK(KE2MAT-1+(IZ-1)*KZYDIM+JZ)
            ZAJI  = WRK(KE2MAT-1+(JZ-1)*KZYDIM+IZ)
            YAIJ  = WRK(KE2MAT-1+(IY-1)*KZYDIM+JY)
            YAJI  = WRK(KE2MAT-1+(JY-1)*KZYDIM+IY)
            ZDEV  = MAX(ABS(ZAIJ-ZAJI),ABS(ZAIJ-YAIJ))
            ZDEV  = MAX(ZDEV,ABS(ZAIJ-YAJI))
            IF (ZDEV.GT.ZAMAX) THEN
               ZAMAX = ZDEV
               IA    = IZ
               JA    = JZ
            END IF
            ZBIJ  = WRK(KE2MAT-1+(IZ-1)*KZYDIM+JY)
            ZBJI  = WRK(KE2MAT-1+(JZ-1)*KZYDIM+IY)
            YBIJ  = WRK(KE2MAT-1+(IY-1)*KZYDIM+JZ)
            YBJI  = WRK(KE2MAT-1+(JY-1)*KZYDIM+IZ)
            ZDEV  = MAX(ABS(ZBIJ-ZBJI),ABS(ZBIJ-YBIJ))
            ZDEV  = MAX(ZDEV,ABS(ZBIJ-YBJI))
            IF (ZDEV.GT.ZBMAX) THEN
               ZBMAX = ZDEV
               IB    = IZ
               JB    = JZ
            END IF
            ZSIJ  = WRK(KS2MAT-1+(IZ-1)*KZYDIM+JZ)
            ZSJI  = WRK(KS2MAT-1+(JZ-1)*KZYDIM+IZ)
            YSIJ  =-WRK(KS2MAT-1+(IY-1)*KZYDIM+JY)
            YSJI  =-WRK(KS2MAT-1+(JY-1)*KZYDIM+IY)
            ZDEV  = MAX(ABS(ZSIJ-ZSJI),ABS(ZSIJ-YSIJ))
            ZDEV  = MAX(ZDEV,ABS(ZSIJ-YSJI))
            IF (ZDEV.GT.ZSMAX) THEN
               ZSMAX = ZDEV
               IS    = IZ
               JS    = JZ
            END IF
            ZDIJ  = WRK(KS2MAT-1+(IZ-1)*KZYDIM+JY)
            ZDJI  =-WRK(KS2MAT-1+(JZ-1)*KZYDIM+IY)
            YDIJ  =-WRK(KS2MAT-1+(IY-1)*KZYDIM+JZ)
            YDJI  = WRK(KS2MAT-1+(JY-1)*KZYDIM+IZ)
            ZDEV  = MAX(ABS(ZDIJ-ZDJI),ABS(ZDIJ-YDIJ))
            ZDEV  = MAX(ZDEV,ABS(ZDIJ-YDJI))
            IF (ZDEV.GT.ZDMAX) THEN
               ZDMAX = ZDEV
               ID    = IZ
               JD    = JZ
            END IF
 160     CONTINUE
 150  CONTINUE
      IZ=IA
      JZ=JA
      IY=IZ+KZDIM
      JY=JZ+KZDIM
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' A(I,J)-MATRIX : MAXIMUM DEVIATION, I=',IA,' J=',JA,' DEV=',ZAMAX
      IF ( ZAMAX.GT.DTOL )
     *WRITE(LUPRI,'(/4(/A,I5,A,I5,A,1P,G16.8))')
     *' IZ=',IZ,' JZ=',JZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(IZ-1)*KZYDIM+JZ),
     *' IZ=',JZ,' JZ=',IZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(JZ-1)*KZYDIM+IZ),
     *' IZ=',IY,' JZ=',JY,' E2(IZ,JZ)',WRK(KE2MAT-1+(IY-1)*KZYDIM+JY),
     *' IZ=',JY,' JZ=',IY,' E2(IZ,JZ)',WRK(KE2MAT-1+(JY-1)*KZYDIM+IY)
      IZ=IB
      JZ=JB
      IY=IZ+KZDIM
      JY=JZ+KZDIM
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' B(I,J)-MATRIX : MAXIMUM DEVIATION, I=',IB,' J=',JB,' DEV=',ZBMAX
      IF ( ZBMAX.GT.DTOL )
     *WRITE(LUPRI,'(/4(/A,I5,A,I5,A,1P,G16.8))')
     *' IZ=',IZ,' JZ=',JY,' E2(IZ,JZ)',WRK(KE2MAT-1+(IZ-1)*KZYDIM+JY),
     *' IZ=',JZ,' JZ=',IY,' E2(IZ,JZ)',WRK(KE2MAT-1+(JZ-1)*KZYDIM+IY),
     *' IZ=',IY,' JZ=',JZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(IY-1)*KZYDIM+JZ),
     *' IZ=',JY,' JZ=',IZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(JY-1)*KZYDIM+IZ)
      IZ=IS
      JZ=JS
      IY=IZ+KZDIM
      JY=JZ+KZDIM
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' S(I,J)-MATRIX : MAXIMUM DEVIATION, I=',IS,' J=',JS,' DEV=',ZSMAX
      IF ( ZSMAX.GT.DTOL )
     *WRITE(LUPRI,'(/4(/A,I5,A,I5,A,1P,G16.8))')
     *' IZ=',IZ,' JZ=',JY,' E2(IZ,JZ)',WRK(KE2MAT-1+(IZ-1)*KZYDIM+JY),
     *' IZ=',JZ,' JZ=',IY,' E2(IZ,JZ)',WRK(KE2MAT-1+(JZ-1)*KZYDIM+IY),
     *' IZ=',IY,' JZ=',JZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(IY-1)*KZYDIM+JZ),
     *' IZ=',JY,' JZ=',IZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(JY-1)*KZYDIM+IZ)
      IZ=ID
      JZ=JD
      IY=IZ+KZDIM
      JY=JZ+KZDIM
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' D(I,J)-MATRIX : MAXIMUM DEVIATION, I=',IA,' J=',JA,' DEV=',ZDMAX
      IF ( ZDMAX.GT.DTOL )
     *WRITE(LUPRI,'(/4(/A,I5,A,I5,A,1P,G16.8))')
     *' IZ=',IZ,' JZ=',JY,' E2(IZ,JZ)',WRK(KE2MAT-1+(IZ-1)*KZYDIM+JY),
     *' IZ=',JZ,' JZ=',IY,' E2(IZ,JZ)',WRK(KE2MAT-1+(JZ-1)*KZYDIM+IY),
     *' IZ=',IY,' JZ=',JZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(IY-1)*KZYDIM+JZ),
     *' IZ=',JY,' JZ=',IZ,' E2(IZ,JZ)',WRK(KE2MAT-1+(JY-1)*KZYDIM+IZ)
C
C CHECK EIGENVALUES OF E2
C
C ALLOCATE WORK SPACE
C
      KE2    = KWRK1
      KS2    = KE2   + KZYDIM*KZYDIM
      KEIVAL = KS2   + KZYDIM*KZYDIM
      KEIVEC = KEIVAL+ KZYDIM
      KWRK1  = KEIVEC+ KZYDIM*KZYDIM
      LWRK1  = LWRK  - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('E2OCHK',KWRK1-1,LWRK)
C
      WRITE(LUPRI,'(/A,/A)')' ********************************',
     *' CHECK EIGENVALUES OF E2 WITH RSPRED ROUTINE '
      IJ = 0
      DO 500 J = 1,KZYDIM
         DO 600 I = 1,J
            IJ = IJ + 1
            WRK(KE2-1+IJ) = WRK(KE2MAT-1+(I-1)*KZYDIM+J)
            IF (I.EQ.J) THEN
               WRK(KS2-1+IJ) = D1
            ELSE
               WRK(KS2-1+IJ) = D0
            END IF
 600     CONTINUE
 500  CONTINUE
      ICTL = 2
      NSIM = 0
      KZSAV = KZRED
      KZYSAV= KZYRED
      KZRED = KZDIM
      KZYRED= KZYDIM
      IBTYP = -1000
      CALL RSPRED(ICTL,.FALSE.,NSIM,IBTYP,GD,REDGD,WRK(KE2),WRK(KS2),
     *            WRK(KEIVAL),WRK(KEIVEC),WRK(KWRK1),WRK(KWRK1),
     *            UDV,WRK(KWRK1),XINDX,WRK(KWRK1),LWRK1)
      KZRED = KZSAV
      KZYRED= KZYSAV
C
C CHECK EIGENVALUES OF S2
C
      WRITE(LUPRI,'(/A,/A)')' ********************************',
     *' CHECK EIGENVALUES OF S2 WITH RSPRED ROUTINE '
      IJ = 0
      DO 700 J = 1,KZYDIM
         DO 800 I = 1,J
            IJ = IJ + 1
            WRK(KE2-1+IJ) = WRK(KS2MAT-1+(I-1)*KZYDIM+J)
            IF (I.EQ.J) THEN
               WRK(KS2-1+IJ) = D1
            ELSE
               WRK(KS2-1+IJ) = D0
            END IF
 800     CONTINUE
 700  CONTINUE
      ICTL = 2
      KZSAV = KZRED
      KZYSAV= KZYRED
      KZRED = KZDIM
      KZYRED= KZYDIM
      IBTYP = -1000
      CALL RSPRED(ICTL,.FALSE.,NSIM,IBTYP,GD,REDGD,WRK(KE2),WRK(KS2),
     *            WRK(KEIVAL),WRK(KEIVEC),WRK(KWRK1),WRK(KWRK1),
     *            UDV,WRK(KWRK1),XINDX,WRK(KWRK1),LWRK1)
      KZRED = KZSAV
      KZYRED= KZYSAV
C
C END OF EO2CHK
C
      RETURN
      END
C  /* Deck lrnsl */
      SUBROUTINE LRNSL(CMO,UDV,PV,FC,XINDX,VECB,WRK,LWRK)
C
#include <implicit.h>
#include <dummy.h>
C
C PURPOSE:
C CALCULATION OF LINEAR POLARIZABILITIES WITHOUT USING
C SELLERS FORMULA.
C
C Modified version of LRHYP, done by K.Ruud on the plane from SLC to SAN
C on July 1, 2000
C
      LOGICAL FOUND, CONV
      CHARACTER*8 BLANK
      PARAMETER (BLANK = '        ', D0 = 0.0D0)
      DIMENSION CMO(*),UDV(*),PV(*),FC(*)
      DIMENSION XINDX(*),VECB(*), WRK(LWRK)
C
#include <infrsp.h>
#include <maxorb.h>
#include <infvar.h>
#include <priunit.h>
#include <inftap.h>
#include <rspprp.h>
#include <inflr.h>
#include <inforb.h>
#include <infpri.h>
#include <wrkrsp.h>
C
C
      WRITE(LUPRI,'(A)') ' '
      CALL PRSYMB(LUPRI,'=',70,1)
      WRITE(LUPRI,'(2A)')
     &     ' >>>>>>>>    L I N E A R   R E S P O N S E   ',
     &     'F U N C T I O N   <<<<<<<<'
      CALL PRSYMB(LUPRI,'=',70,1)
C
      KVECA = 1
C
      IF (NGPLR(KSYMOP).LE.0) RETURN
C
C     Define variables that depend on symmetry
C
      DO IOPA = 1, NGPLR(KSYMOP)
         CALL GETGPV(LBLLR(KSYMOP,IOPA),FC,DUMMY,CMO,UDV,PV,XINDX,
     *               ANTSYM,WRK(KVECA),LWRK)
C
C     Loop over solution vectors
C
         DO IFREQ = 1, NFREQ
            DO IOPB = 1, NGPLR(KSYMOP)
               CALL REARSP(LURSP,KLEN,VECB,LBLLR(KSYMOP,IOPB),
     &                     BLANK,FREQ(IFREQ),D0,KSYMOP,0,THCLR,
     &                     FOUND,CONV,ANTSYM)
               IF (.NOT. (FOUND .AND. CONV))
     &              CALL RIOERR(FOUND, LBLLR(KSYMOP,IOPB),
     &                          FREQ(IFREQ),KSYMOP)
C
C     Compute linear response value
C
               VAL = DDOT(KZYVAR,WRK(KVECA),1,VECB,1)
C              
               WRITE(LUPRI,'(/A,2(/A,A10,I4,F10.6),/,/A,F20.12)')
     &              ' Linear response function in a.u.',
     &              ' A operator, symmetry, frequency: ',
     &              LBLLR(KSYMOP,IOPA),KSYMOP,-FREQ(IFREQ),
     &              ' B operator, symmetry, frequency: ',
     &                LBLLR(KSYMOP,IOPB),KSYMOP, FREQ(IFREQ),
     &              ' Value of linear response -<<A;B>(omega): ',VAL
            END DO
         END DO
      END DO
C
      RETURN
      END

      subroutine printelweak(nop,sopr,ksymop)
#include <implicit.h>
#include <priunit.h>
#include <codata.h>
#include <rspprp.h>
#include <inflr.h>
#include <mxcent.h>
#include <elweak.h>
#include <nuclei.h>
#include <chrnos.h>

      dimension sopr(1,nop,nop),pvodd(mxcent)
      character*8 lab1, lab2
      character*3 prelab,num
      write(lupri,*) 'pvso,pvso2,pvpso,natoms'
      write(lupri,*) pvso,pvso2,pvpso,natoms
      if (PVSO) THEN
         tot=0.0d0
         do I=1,NOP
            do j=1,NOP
               lab1=lbllr(ksymop,i)
               lab2=lbllr(ksymop,j)
               if ((lab2.eq.'X1SPNORB' .and. lab1 .eq.'PVIOLA X')) then

C     &             (lab1.eq.'X1SPNORB' .and. lab2 .eq.'PVIOLA X')
                  tot=tot+sopr(ksymop,i,j)*2867.081D0
               endif
               if ((lab2.eq.'Y1SPNORB' .and. lab1 .eq.'PVIOLA Y')) then
C     .or.
C     &             (lab1.eq.'Y1SPNORB' .and. lab2 .eq.'PVIOLA Y')) then
                  tot=tot+sopr(ksymop,i,j)*2867.081D0
               endif
               if ((lab2.eq.'Z1SPNORB' .and. lab1 .eq.'PVIOLA Z')) then
c     .or.
C     &             (lab1.eq.'Z1SPNORB' .and. lab2 .eq.'PVIOLA Z')) then
                  tot=tot+sopr(ksymop,i,j)*2867.081D0
               endif
            enddo
         enddo

         tot1el=tot
         call around('One-el pv term based on one-el so-operator')
         write(lupri,'(F20.8,A10/)') tot,'  10^-20 Eh'

      endif

      IF (PVPSO) THEN
         tot=0.0d0
         call dzero(pvodd,natoms)
         do I=1,NOP
            do j=1,NOP
               lab1=lbllr(ksymop,i)
               lab2=lbllr(ksymop,j)
               read(lab2(1:3),'(A3)') PRELAB
C
               if (lab1.eq.'PVIOLA X' .and.PRELAB.eq.'PSO') then
                  read(lab2(5:8),'(I4)') INUM
                  iatom=((inum-1)/3)+1
                  if (mod(inum,3).eq.1) then
C                     write(0,*) lab1, lab2,inum,iatom,mod(inum,3)
                     pvodd(iatom)= pvodd(iatom)+sopr(ksymop,i,j)
                  endif
               endif
               if (lab1.eq.'PVIOLA Y' .and.PRELAB.eq.'PSO') then
                  read(lab2(5:8),'(I4)') INUM
                  iatom=((inum-1)/3)+1
                  if (mod(inum,3).eq.2) then
                     pvodd(iatom)= pvodd(iatom)+sopr(ksymop,i,j)
C                     write(0,*) lab1,lab2,inum,iatom,mod(inum,3)
                  endif
               endif
               if (lab1.eq.'PVIOLA Z' .and.PRELAB.eq.'PSO') then
                  read(lab2(5:8),'(I4)') INUM
                  iatom=((inum-1)/3)+1
                  if (mod(inum,3).eq.0) then
C                     write(0,*)lab1,lab2,inum,iatom,mod(inum,3)
                     pvodd(iatom)= pvodd(iatom)+sopr(ksymop,i,j)
                  endif
               endif
            enddo
         enddo
         write(lupri,*)
         call around('One-el term based on pso-integrals')
         write(lupri,*)
         write(lupri,*) 'Name   Charge             PV '
         do iatom=1,natoms
            pvodd(iatom) = pvodd(iatom)*charge(iatom)/4.0d0
            pvodd(iatom) = pvodd(iatom)*2867.081D0*alpha2
            write(lupri,'(A7,F8.0,F20.8,A10)') namn(iatom),charge(iatom)
     &           ,pvodd(iatom),' 10^-20 Eh'
            tot = tot +  pvodd(iatom)
         enddo
         call around('Total One-el term based on the pso integrals')
         write(lupri,'(F20.8,A10/)') tot,'  10^-20 Eh'

         if (.not. pvso) tot1el = tot
      endif

      if (PVSO2) THEN
         tot = 0.0d0
         do I=1,NOP
            do j=1,NOP
               lab1=lbllr(ksymop,i)
               lab2=lbllr(ksymop,j)
               if (lab2.eq.'X2SPNORB' .and. lab1 .eq.'PVIOLA X') then
                  tot=tot+sopr(ksymop,i,j)*2867.081D0
               endif
               if (lab2.eq.'Y2SPNORB' .and. lab1 .eq.'PVIOLA Y') then
                  tot=tot+sopr(ksymop,i,j)*2867.081D0
               endif
               if (lab2.eq.'Z2SPNORB' .and. lab1 .eq.'PVIOLA Z') then
                  tot=tot+sopr(ksymop,i,j)*2867.081D0
               endif
            enddo
         enddo

         tot2el=tot
         call around('Two-el term based on two-el spin orbit operator')
         write(lupri,'(F20.8,A10)') tot,'  10^-20 Eh'
         write(lupri,*)
      endif

      call around('Total PV-INTERACTION')
      write(lupri,'(F20.8,A10)') (tot1el+tot2el),' 10^-20 Eh'
      write(lupri,*)

      return
      end
C  /* Deck RSPPMO */
      SUBROUTINE RSPPMO (WOP,IOFFY,CMO2,LUPRI)
C
C     June 2008; PBHT MO overlap diagnostic
C     Reference: 
C       Excitation energies in density functional theory: An evaluation
C       and a diagnostic test
C       Michael J. G. Peach, Peter Benfield, Trygve Helgaker, and David J. Tozer,
C       J.  Chem. Phys. 128, 044118 (2008), DOI:10.1063/1.2831900
C
C     (subroutine based on RSPPRO tuh)
C
#include <implicit.h>
      PARAMETER(PTHRS = 1.D-08)
      DIMENSION WOP(*)
      DIMENSION CMO2(NORBT,NORBT,2)
C
#include <maxorb.h>
#include <maxash.h>
#include <infind.h>
#include <inforb.h>
#include <infvar.h>
 
      CALL HEADER('PBHT MO Overlap Diagnostic',-1)
      WRITE(LUPRI,'(A/A/)')
     &  '  Reference: MJG Peach, P Benfield, T Helgaker, and DJ Tozer.',
     &  '             J Chem Phys 128, 044118 (2008)'
      IF (NWOPT.GT.0) THEN 
C        WRITE (LUPRI,1010) JWOPSY
C        WRITE (LUPRI, '(2X,2A,/)')  '            I    A       ',
C    &       'K_IA           K_AI         <|I|*|A|>      <I^2*A^2>'
         CMOTOT=0.d0
         WOPTOT=0.d0
         DO IG = 1,NWOPT
            K  = JWOP(1,IG)
            L  = JWOP(2,IG)
            KS = ISMO(K)
            LS = ISMO(L)
            WOP1 = WOP(IG)
            WOP2 = WOP(IOFFY+IG)
            IF (ABS(WOP1).GT.PTHRS .OR. ABS(WOP2).GT.PTHRS) THEN
C              WRITE (LUPRI,'(2X,A,3I5,4F15.10)') 
C    &            '@@@', IG,K,L,WOP1,WOP2,CMO2(K,L,1),CMO2(K,L,2)
                  WMJGP=(WOP1-WOP2)**2.0d0
                  CMOMJGP=WMJGP*CMO2(K,L,1)
                  CMOTOT=CMOTOT+CMOMJGP
                  WOPTOT=WOPTOT+WMJGP
            END IF
         END DO
         ALAMBDA=CMOTOT/WOPTOT
         WRITE(LUPRI,'(25X,"Lambda:",2X,F15.10)')ALAMBDA
      END IF
      RETURN
 1010 FORMAT(5X,'Response orbital operator symmetry =',I2,/)
      END
