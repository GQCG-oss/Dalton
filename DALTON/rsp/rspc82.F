!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
C
#ifdef UNDEF
========================================================================
/* Comdeck rev_log */
Revision 1.2  2001/01/17 10:41:29  vebjornb
Calls to *MPA*B* in arhpack.F have been replaced with DGEMM calls

Revision 1.3  2000/05/24 18:48:53  hjj
new GETREF calls with NDREF (fixing triplet with CSF)
do not save S(2) CI diagonal on LURSP4 ( all ones)

revision 1.2 2000/04/14 14:33:01 hjj
partial bug fix for triplet and CSFs (more need to be done)
========================================================================
950208-hjaaj
RSPORT: s/ZEQLY/OVLMIN/ (OVLMIN = 1.0D-4; ZEQLY was 1.0D-6 but they
  tested the same thing; Hans Agren had problems for cubic response with
  ZEQLY=1.0D-6 and OVLMIN=1.0D-4);  improved code for Z=Y or Z=-Y tests.
940708-hjaaj: SOPPA changes
931004-pj/hjaaj PPST  : fixed DCOPY problem for RHF
930910-hjaaj    RSPNEX: always use relative convergence
930712-hjaaj
RSPORT: ZEQLY = 1.0D-6 instead of 1.0D-10; lowered print level;
        discard trial vector for linear dep. if OVLPI close to 1;
        do symmetric orthonormalization twice
920929-ov      RSPSLI: fixed work dimension check
920904-hjaaj
RSPNEX: 1) fixed error in print statement
        -- corrected print of solution and residual vectors
        2) fixed error for multiple orbital trial vectors, when
        one or more of the lower root(s) converged by
        setting KEX(i) = -1 for converged roots.
        3) skip S[2]*X(i) calculation if EIVAL(i) .eq. D0
        4) only DZERO(BVECS(KINFIX)) if OPTORB
920721-hinne hettema
  RSPORT,RSPSOR: implemented RSPSUP (super symmetry averaging)
========================================================================
#endif
C  /* Deck rsplin */
      SUBROUTINE RSPLIN(NCSIM,NOSIM,ZYCVEC,ZYOVEC,
     *                  CMO,UDV,PVX,FC,FV,FCAC,H2AC,
     *                  XINDX,WRK,LWRK)
C
C PURPOSE:
C  CARRY OUT LINEAR TRANSFORMATIONS
C
C     EVECS(I) = E[2]*N(I)
C     SVECS(I) = S[2]*N(I)
C
C  THE TRIAL VECTORS HAVE EITHER CONFIGURATION (ZYCVEC) OR
C  (ZYOVEC) ORBITAL COMPONENTS
C
C OUTPUT :
C
C     EVECS(I) AND SVECS(I)
C
C  THE FIRST (NCSIM+NOSIM)*KZYVAR VARIABLES IN WRK CONTAIN EVECS
C  THE NEXT  (NCSIM+NOSIM)*KZYVAR VARIABLES IN WRK CONTAIN SVECS
C
#include "implicit.h"
      INTEGER ABSADR, ABTADR
C
      DIMENSION ZYCVEC(*),ZYOVEC(*)
      DIMENSION CMO(*),UDV(*),PVX(*),FC(*),FV(*),FCAC(*)
      DIMENSION H2AC(*),XINDX(*),WRK(*)
C
C  INFDIM : NASHDI
C
#include "inforb.h"
#include "infdim.h"
#include "inftap.h"
#include "wrkrsp.h"
#include "infrsp.h"
C
C ALLOCATE WORK SPACE FOR RSPLIN
C
      CALL QENTER('RSPLIN')
      KELI   = 1
      KSLI   = KELI + (NCSIM+NOSIM) * KZYVAR
      KWRKS  = KSLI + (NCSIM+NOSIM) * KZYVAR
      LWRKS  = LWRK   - KWRKS
      IF (LWRKS.LT.0) CALL ERRWRK('RSPLIN',KWRKS-1,LWRK)
      CALL RSPELI(NCSIM,NOSIM,ZYCVEC,ZYOVEC,
     *            CMO,UDV,PVX,FC,FV,FCAC,H2AC,
     *            XINDX,WRK,LWRK)
C     CALL RSPELI(NCSIM,NOSIM,ZYCVEC,ZYOVEC,
C    *            CMO,UDV,PVX,FC,FV,FCAC,H2AC,
C    *            XINDX,WRK,LWRK)
      NTOT = (NCSIM+NOSIM) * KZYVAR
      CALL DZERO(WRK(KSLI),NTOT)
      CALL RSPSLI(NCSIM,NOSIM,ZYCVEC,ZYOVEC,
     *            UDV,WRK(KSLI),XINDX,WRK(KWRKS),LWRKS)
      CALL QEXIT('RSPLIN')
      RETURN
      END
C  /* Deck rspsli */
      SUBROUTINE RSPSLI(NCSIM,NOSIM,ZYCVEC,ZYOVEC,
     *                  UDV,SVEC,XINDX,WRK,LWRK)
C
C PURPOSE:
C  CARRY OUT LINEAR TRANSFORMATIONS
C
C     SVEC(I) = S[2]*N(I)
C
C  THE TRIAL VECTORS HAVE EITHER CONFIGURATION (ZYCVEC) OR
C  (ZYOVEC) ORBITAL COMPONENTS
C
C OUTPUT :
C
C      SVEC(I)
C
C  THE FIRST (NCSIM+NOSIM)*KZYVAR VARIABLES IN WRK CONTAIN SVECS
C
#include "implicit.h"
C
      DIMENSION ZYCVEC(*),ZYOVEC(*)
      DIMENSION UDV(*),SVEC(*),XINDX(*),WRK(*)
C
C  INFDIM : NASHDI
C
#include "dummy.h"
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "inforb.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "infind.h"
#include "infpri.h"
#include "infdim.h"
C
      LOGICAL NOH2
C
      PARAMETER ( D1 = 1.0D0 , DM1 = -1.0D0 )
C
      CALL QENTER('RSPSLI')
C
C
C ALLOCATE WORK SPACE
C
      IF (TDHF .OR. KZCONF.EQ.0) THEN
C        ... no CREF needed for TDHF and for "OPTORB" /HJAaJ
         NDREF = 0
      ELSE IF (IREFSY .EQ. KSYMST) THEN
         NDREF = KZCONF
C        ... either KZCONF .eq. NCREF or
C            KZCONF .gt. NCREF and then we need CREF in KZCONF determinants
C            instead of in NCREF CSFs.
      ELSE
         NDREF = NCREF
      END IF
      KZYMAT = 1
      KDVT   = KZYMAT + NOSIM * N2ORBX
      KCREF  = KDVT   + NCSIM * N2ASHX
      KUFCAC = KCREF  + NDREF
      KFREE  = KUFCAC + N2ASHX
      LFREE  = LWRK   - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('RSPSLI',KFREE-1,LWRK)
      IF (NDREF.GT.0) CALL GETREF(WRK(KCREF),NDREF)
      KSTRT  = 1
C
C CONSTRUCT CONFIGURATION PART OF LINEAR TRANSFORMED
C CONFIGURATION TRIAL VECTORS
C S[2] configuration block is a unit matrix both for
C MCSCF and SOPPA
C
      IF ( NCSIM .GT. 0 ) THEN
         DO 100 ISIM = 1,NCSIM
            KSZOFF = 1 + (ISIM-1)*KZYVAR
            CALL DAXPY(KZCONF,D1,ZYCVEC(1+(ISIM-1)*KZCONF),
     &                 1,SVEC(KSZOFF),1)
 100     CONTINUE
         IF (IPRRSP.GT.110) THEN
            WRITE(LUPRI,'(/A)') ' S(2) linear transformed'//
     &         ' configuration part (X1 Y1 X2 Y2 ...)'
            CALL OUTPUT(SVEC,1,KZVAR,1,2*NCSIM,
     &                  KZVAR,2*NCSIM,-1,LUPRI)
         END IF
      END IF
C
C UNPACK ORBITAL VECTORS
C
      IF ( NOSIM .GT. 0 ) THEN
         CALL RSPZYM(NOSIM,ZYOVEC,WRK(KZYMAT))
C        CALL RSPZYM(NSIM,ZYVEC,ZYMAT)
C
      END IF
C
C CONSTRUCT CONFIGURATION PART OF LINEAR TRANSFORMED
C ORBITAL TRIAL VECTORS
C The ph-2p2h coupling block of S[2] is zero for SOPPA
C
      IF ((NOSIM.GT.0).AND.(KZCONF.GT.0).AND.(.NOT.SOPPA)) THEN
         KCVETO = NCSIM*KZYVAR
         DO 200 ISIM = 1,NOSIM
C
C CREATE Z PART OF LINEAR TRANSFORMATION WITH S(2)
C
            DO 420 IW = 1,NASHT
               IX = ISX(NISHT+IW)
               DO 430 JW = 1,NASHT
                  JX = ISX(NISHT+JW)
                  IJW = (IW-1) * NASHT + JW
                  IXJX = (ISIM-1)*NORBT*NORBT + (JX-1)*NORBT +IX
                  WRK(KUFCAC-1+IJW) = WRK(KZYMAT-1+IXJX)
 430           CONTINUE
 420        CONTINUE
C
            IF (IPRRSP.GT.100) THEN
               WRITE(LUPRI,'(/A)')' ACTIVE PART OF  ZYMAT '
               CALL OUTPUT(WRK(KUFCAC),1,NASHT,1,NASHT,
     *                     NASHT,NASHT,-1,LUPRI)
            END IF
            IF (TRPLET) THEN
               ISPIN1 = 0
               ISPIN2 = 1
            ELSE
               ISPIN1 = 0
               ISPIN2 = 0
            END IF
            CALL CISIGD(IREFSY,KSYMST,NDREF,KZCONF,WRK(KCREF),
     *                  SVEC(KCVETO+1+(ISIM-1)*KZYVAR),
     *                  WRK(KUFCAC),DUMMY,.TRUE.,.TRUE.,
     *                  XINDX,ISPIN1,ISPIN2,WRK(KFREE),LFREE)
C
            IF (IREFSY .EQ. KSYMST .AND. .NOT.TRPLET) THEN
C           ... remove CREF component of S[2] vector
               T1 = DDOT(KZCONF,WRK(KCREF),1,
     *                   SVEC(KCVETO+1+(ISIM-1)*KZYVAR),1)
               CALL DAXPY(KZCONF,(-T1),WRK(KCREF),1,
     *                    SVEC(KCVETO+1+(ISIM-1)*KZYVAR),1)
            END IF
            IF (IPRRSP.GT.110) THEN
               WRITE(LUPRI,'(/A)')
     *        ' CSF PART OF S(2) TRANSFORMED ORBITAL VECTOR: Z PART'
               CALL OUTPUT(SVEC(KCVETO+1+(ISIM-1)*KZYVAR),1,1,1,KZCONF,
     *                     1,KZCONF,-1,LUPRI)
            END IF
C
C CREATE Y PART OF S(2) LINEAR TRANSFORMED TRIAL VECTOR
C
C        Get transposed ZYMAC matrix for Y sigma vector
C        (note: CISIGD requires UFCAC(I,J) = ZYMAC(J,I))
C
            DO 320 IW = 1,NASHT
               IX = ISX(NISHT+IW)
               DO 330 JW = 1,NASHT
                  JX = ISX(NISHT+JW)
                  IJW = (IW-1) * NASHT + JW
                  JXIX = (ISIM-1)*NORBT*NORBT + (IX-1)*NORBT +JX
                  WRK(KUFCAC-1+IJW) = -WRK(KZYMAT-1+JXIX)
C THE MINUS SIGN IN LINEAR TRANSFORMATION IS PUT INTO THE ONE
C ELECTRON OPERATOR
 330           CONTINUE
 320        CONTINUE
            IF (IPRRSP.GT.100) THEN
               WRITE(LUPRI,'(/A)')' ACTIVE PART OF TRANSPOSED ZYMAT '
               CALL OUTPUT(WRK(KUFCAC),1,NASHT,1,NASHT,
     *                     NASHT,NASHT,-1,LUPRI)
            END IF
            IF (TRPLET) THEN
               ISPIN1 = 0
               ISPIN2 = 1
            ELSE
               ISPIN1 = 0
               ISPIN2 = 0
            END IF
            CALL CISIGD(IREFSY,KSYMST,NDREF,KZCONF,WRK(KCREF),
     *                  SVEC(KCVETO+1+(ISIM-1)*KZYVAR+KZVAR),
     *                  WRK(KUFCAC),DUMMY,.TRUE.,.TRUE.,
     *                  XINDX,ISPIN1,ISPIN2,WRK(KFREE),LFREE)
            IF (IREFSY .EQ. KSYMST .AND. .NOT.TRPLET) THEN
C
               T1 = DDOT(KZCONF,WRK(KCREF),1,
     *                   SVEC(KCVETO+1+(ISIM-1)*KZYVAR+KZVAR),1)
               CALL DAXPY(KZCONF,(-T1),WRK(KCREF),1,
     *                    SVEC(KCVETO+1+(ISIM-1)*KZYVAR+KZVAR),1)
            END IF
            IF (IPRRSP.GT.110) THEN
               WRITE(LUPRI,'(/A)')
     *        ' CSF PART OF S(2) TRANSFORMED ORBITAL VECTOR: Y PART'
               CALL OUTPUT(SVEC(KCVETO+1+(ISIM-1)*KZYVAR+KZVAR),
     *                     1,1,1,KZCONF,1,KZCONF,-1,LUPRI)
            END IF
 200     CONTINUE
      END IF
C
C CONSTRUCT ORBITAL PART OF LINEAR TRANSFORMED VECTORS
C
C CALCULATE TRANSITION DENSITY MATRIX
C The ph-2p2h coupling block of S[2] is zero for SOPPA
C
      IF ( (NCSIM.GT.0) .AND. (.NOT.RSPCI) .AND. (.NOT.SOPPA)) THEN
         IF (TRPLET) THEN
            ISPIN1 = 1
            ISPIN2 = 0
         ELSE
            ISPIN1 = 0
            ISPIN2 = 0
         END IF
         CALL RSPTDM(NCSIM,IREFSY,KSYMST,NDREF,KZCONF,WRK(KCREF),
     *                 ZYCVEC,WRK(KDVT),DUMMY,
     *                 ISPIN1,ISPIN2,.TRUE.,.TRUE.,
     *                 XINDX,WRK(KFREE),KSTRT,LFREE)
C        CALL RSPTDM(NCSIM,ILRESY,IRSYM,NCLREF,NCRDIM,CLREF,
C    *                 CR, RHO1,RHO2, ISPIN1,ISPIN2,TDM,NORHO2,
C    *                 XNDXCI,WORK,KFREE,LFREE)
         CALL DSCAL(NCSIM*N2ASHX,DM1,WRK(KDVT),1)
C THE MINUS SIGN IN THE LINEAR TRANSFORMATION IS PUT
C INTO THE ONE ELECTRON TRANSITION DENSITY MATRIX
         CALL RSPSOR(.FALSE.,NCSIM,WRK(KDVT),SVEC,VDUMMY)
      END IF
      IF (NOSIM.GT.0) THEN
         IF (SOPPA) THEN
            CALL RSPS2M(NOSIM,UDV,SVEC(1+NCSIM*KZYVAR),
     &                  WRK(KZYMAT),WRK(KFREE))
         ELSE
            CALL RSPSOR(.TRUE.,NOSIM,UDV,SVEC(1+NCSIM*KZYVAR),
     *                  WRK(KZYMAT))
C           CALL RSPSOR(ONEIND,NSIM,UDV,SVECS,ZYMAT)
         ENDIF
      END IF
C      IF ( (NCSIM.GT.0) .AND. (.NOT.RSPCI) ) THEN
C         CALL RSPSOR(.FALSE.,NCSIM,WRK(KDVT),SVEC,VDUMMY)
C      END IF
C
C END OF RSPSLI
C
      CALL QEXIT('RSPSLI')
      RETURN
      END
C  /* Deck rspsor */
      SUBROUTINE RSPSOR(ONEIND,NSIM,UDV,SVECS,ZYMAT)
C
C WRITTEN 14-FEB 1986
C List of updates
C 21-Jul-1992 Hinne Hettema (super) symmetry averaging added
C
C PURPOSE:
C
C   CREATE ORBITAL PART OF LINEAR TRANSFORMATION WITH S(2)
C
C    ONEIND = .TRUE. FOR A ORBITAL TRIAL VECTOR
C
C                ******************************
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
      LOGICAL ONEIND
C
      DIMENSION UDV(NASHDI,NASHDI,*)
      DIMENSION SVECS(KZYVAR,*)
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
      CALL QENTER('RSPSOR')
      KYCONF = KZCONF + KZVAR
C
C DISTRIBUTE FOCK MATRICES IN EVECS
C
         K_SYM1 = 0
         DO 1300 IG = 1,KZWOPT
            K     = JWOP(1,IG)
            L     = JWOP(2,IG)
            K_SYM  = ISMO(K)
            L_SYM  = ISMO(L)
            IF( K_SYM.NE.K_SYM1 ) THEN
               K_SYM1 = K_SYM
               IORBK = IORB(K_SYM)
               NASHK = NASH(K_SYM)
               NISHK = NISH(K_SYM)
               IASHK = IASH(K_SYM)
               IIORBK= IIORB(K_SYM)
               IORBL = IORB(L_SYM)
               NASHL = NASH(L_SYM)
               NISHL = NISH(L_SYM)
               IASHL = IASH(L_SYM)
            END IF
            ITYPK = IOBTYP(K)
            ITYPL = IOBTYP(L)
            IF ( ITYPK.EQ.JTINAC )THEN
               DO 2000 ISIM = 1 ,NSIM
                  IF (ONEIND) THEN
                     SVECS(KZCONF+IG,ISIM) = SVECS(KZCONF+IG,ISIM)
     *                       + D2 * ZYMAT(L,K,ISIM)
                     IF ((.NOT. TDA) .AND. (.NOT. CISRPA))
     *               SVECS(KYCONF+IG,ISIM) = SVECS(KYCONF+IG,ISIM)
     *                       - D2 * ZYMAT(K,L,ISIM)
                  END IF
                  IF ( ITYPL.EQ.JTACT ) THEN
                     NWL = ISW(L) - NISHT
                     IF (ONEIND) THEN
                        SVECS(KZCONF+IG,ISIM) = SVECS(KZCONF+IG,ISIM)
     *                  - DDOT(NASHL,ZYMAT(IORBL+NISHL+1,K,ISIM),1,
     *                               UDV(IASHL+1,NWL,1),1)
                        IF ((.NOT. TDA) .AND. (.NOT. CISRPA)) THEN
                           DO IX = 1,NASHL
                              SVECS(KYCONF+IG,ISIM)
     *                        = SVECS(KYCONF+IG,ISIM)
     *                        + ZYMAT(K,IORBL+NISHL+IX,ISIM)*
     *                          UDV(NWL,IASHL+IX,1)
                           END DO
                        END IF
                     END IF
                  END IF
 2000          CONTINUE
            ELSE
              IF (ITYPL.EQ.JTACT) THEN
                  NWL = ISW(L) - NISHT
                  NWK = ISW(K) - NISHT
                  DO 2020 ISIM=1,NSIM
                     IF (ONEIND) THEN
                        SVECS(KZCONF+IG,ISIM) = SVECS(KZCONF+IG,ISIM)
     *                     -DDOT(NASHL,ZYMAT(IORBL+NISHL+1,K,ISIM),1,
     *                                        UDV(IASHL+1,NWL,1),1)
                        IF ((.NOT. TDA) .AND. (.NOT. CISRPA)) THEN
                           DO IX = 1,NASHL
                              SVECS(KYCONF+IG,ISIM)
     *                        = SVECS(KYCONF+IG,ISIM)
     *                        + ZYMAT(K,IORBL+NISHL+IX,ISIM)*
     *                          UDV(NWL,IASHL+IX,1)
                           END DO
                           SVECS(KYCONF+IG,ISIM) = SVECS(KYCONF+IG,ISIM)
     *                        -DDOT(NASHK,ZYMAT(IORBK+NISHK+1,L,ISIM),1,
     *                                      UDV(IASHK+1,NWK,1),1)
                        END IF
                        DO IX = 1,NASHK
                           SVECS(KZCONF+IG,ISIM) = SVECS(KZCONF+IG,ISIM)
     *                        + ZYMAT(L,IORBK+NISHK+IX,ISIM)*
     *                          UDV(NWK,IASHK+IX,1)
                        END DO
                     ELSE
                        SVECS(KZCONF+IG,ISIM) =
     *                     SVECS(KZCONF+IG,ISIM) + UDV(NWK,NWL,ISIM)
                        IF ((.NOT. TDA) .AND. (.NOT. CISRPA))
     *                  SVECS(KYCONF+IG,ISIM) =
     *                     SVECS(KYCONF+IG,ISIM) + UDV(NWL,NWK,ISIM)
                     END IF
 2020             CONTINUE
               ELSE
                  NWK = ISW(K) - NISHT
                  DO 2030 ISIM=1,NSIM
                     IF (ONEIND) THEN
                        IF ((.NOT. TDA) .AND. (.NOT. CISRPA))
     *                     SVECS(KYCONF+IG,ISIM) = SVECS(KYCONF+IG,ISIM)
     *                     -DDOT(NASHK,ZYMAT(IORBK+NISHK+1,L,ISIM),1,
     *                                 UDV(IASHK+1,NWK,1),1)
                        DO IX = 1,NASHK
                           SVECS(KZCONF+IG,ISIM) = SVECS(KZCONF+IG,ISIM)
     *                        + ZYMAT(L,IORBK+NISHK+IX,ISIM)*
     *                          UDV(NWK,IASHK+IX,1)
                        END DO
                     END IF
 2030             CONTINUE
               ENDIF
            ENDIF
 1300    CONTINUE
C
C *** Perform supersymmetry averaging
C
      IF (RSPSUP .AND. (KSYMOP .EQ. 1)) THEN
         IF (IPRRSP.GT.210) THEN
            WRITE(LUPRI,*)' TRANSFORMATION WITH S(2) BEFORE RSPAVE'//
     &      ' (Z1 Y1 Z2 Y2 ...)'
            CALL OUTPUT(SVECS,1,KZVAR,1,2*NSIM,KZVAR,2*NSIM,-1,LUPRI)
         END IF
         CALL RSPAVE(SVECS(KZCONF+1,1),KZVAR,2*NSIM)
      END IF
C
      IF (IPRRSP.GT.110) THEN
         WRITE(LUPRI,*)' LINEAR TRANSFORMATION WITH  S(2)'//
     &      ' (Z1 Y1 Z2 Y2 ...)'
         CALL OUTPUT(SVECS,1,KZVAR,1,2*NSIM,KZVAR,2*NSIM,-1,LUPRI)
      END IF
C
C END OF RSPSOR
C
      CALL QEXIT('RSPSOR')
      RETURN
      END
C  /* Deck rspeli */
      SUBROUTINE RSPELI(NCSIM,NOSIM,ZYCVEC,ZYOVEC,CMO,UDV,PVX,FC,FV,
     *                  FCAC,H2AC,XINDX,WRK,LWRK)
C
C PURPOSE:
C  CARRY OUT LINEAR TRANSFORMATIONS
C
C     EVECS(I) = E[2]*N(I)
C
C  THE TRIAL VECTORS HAVE EITHER CONFIGURATION (ZYCVEC) OR
C  (ZYOVEC) ORBITAL COMPONENTS
C
C OUTPUT :
C
C     EVECS(I)
C
C  THE FIRST (NCSIM+NOSIM)*KZYVAR VARIABLES IN WRK CONTAIN EVECS
C
#ifdef HAS_PCMSOLVER
      use pcm_config, only: pcm_configuration, pcm_cfg
#endif
      use pelib_interface, only: use_pelib
#include "implicit.h"
C
      DIMENSION ZYCVEC(*),ZYOVEC(*)
      DIMENSION CMO(*),UDV(NASHDI,*),PVX(*),FC(*),FV(*),FCAC(*)
      DIMENSION H2AC(*),XINDX(*),WRK(*)
C
C  INFDIM : NASHDI
C
#include "priunit.h"
#include "maxorb.h"
#include "infinp.h"
#include "pcmlog.h"
#include "inforb.h"
#include "infdim.h"
#include "inftap.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "gnrinf.h"
#include "mxcent.h"
C
      logical :: is_pcmsolver = .false.
#ifdef HAS_PCMSOLVER
      if (pcm_cfg%do_pcm) then
        is_pcmsolver = .true.
      end if
#endif
C
      CALL QENTER('RSPELI')
C
C
C ALLOCATE WORK SPACE
C
      KEVECS = 1
      KRSPCL = KEVECS + (NCSIM+NOSIM) * KZYVAR
      LWRKRS = LWRK   - KRSPCL
      IF (LWRKRS.LT.0) CALL ERRWRK('RSPELI',KRSPCL-1,LWRK)
C
C INITIALIZE EVECS
C
      CALL DZERO(WRK,(NCSIM+NOSIM)*KZYVAR)
C
C ALLOCATE WORK SPACE FOR RSPCLI
C
C
C CONSTRUCT CONFIGURATION PART OF LINEAR TRANSFORMED VECTORS
C FOR A CONFIGURATION TRIAL VECTOR
C
      IF ( NCSIM .GT. 0 ) THEN
         IF (SOPPA) THEN
            CALL SOPCLI(NCSIM,ZYCVEC,WRK(KEVECS),WRK(KRSPCL),LWRKRS)
         ELSE
            CALL RSPCLI(NCSIM,ZYCVEC,FCAC,H2AC,
     *                  WRK(KEVECS),XINDX,WRK(KRSPCL),LWRKRS)
C           CALL RSPCLI(NCSIM,ZYCVEC,FCAC,H2AC,EVECS,XINDX,WRK,LWRK)
         END IF
      END IF
C
C ALLOCATE WORK SPACE FOR RSPOLI
C
      KZYMAT = KRSPCL
      KFCX   = KZYMAT + NOSIM * NORBT * NORBT
      KFVX   = KFCX   + NOSIM * NORBT * NORBT
      KQAX   = KFVX   + NOSIM * NORBT * NORBT
      KQBX   = KQAX   + NOSIM * NORBT * NASHT
C
      LFVTD  = NCSIM * N2ORBX
      IF (DOMCSRDFT) LFVTD = 2*LFVTD
C     ... Extra allocation for "FCTD" in MCSCF-SRDFT,
C         is needed for the VxcTD matrices.
      KFVTD  = KQBX   + NOSIM * NORBT * NASHT
      KQATD  = KFVTD  + LFVTD
C
      KQBTD  = KQATD  + NCSIM * NORBT * NASHT
      KWRK1  = KQBTD  + NCSIM * NORBT * NASHT
      LWRK1  = LWRK   - KWRK1
C
C UNPACK ORBITAL VECTORS
C
      IF ( NOSIM .GT. 0 ) THEN
         CALL RSPZYM(NOSIM,ZYOVEC,WRK(KZYMAT))
C        CALL RSPZYM(NSIM,ZYVEC,ZYMAT)
C
      END IF
C
C CONSTRUCT ORBITAL PART OF LINEAR TRANSFORMED VECTORS
C
      IF (.NOT.RSPCI) THEN
        CALL RSPOLI(NCSIM,NOSIM,UDV,ZYCVEC,KZCONF,FC,FV,PVX,WRK(KZYMAT),
     *              WRK(KFCX),WRK(KFVX), WRK(KQAX),WRK(KQBX),
     *              WRK(KFVTD),WRK(KQATD),WRK(KQBTD),WRK(KEVECS),
     *              XINDX,CMO,WRK(KWRK1),LWRK1)
C       CALL RSPOLI(NCSIM,NOSIM,UDV,ZYCVEC,LZYCVEC,FC,FV,PVX,ZYMAT,
C    *              FCX,FVX,QAX,QBX, FVTD,QATD,QBTD,EVECS,
C    *              XINDX,CMO,WRK,LWRK)
      ENDIF

      CALL FLSHFO(LUPRI)
C
C FLAG(16) MEANS SOLVENT CALCULATION with spherical cavity
C
C EMBEDDING = (FLAG(16) .OR. PCM .OR. QM3 .OR. QMMM .OR. QMNPMM)
C
      IF (EMBEDDING .OR. USE_PELIB() .or. is_pcmsolver) THEN
         CALL RSPSLV(NCSIM,NOSIM,ZYCVEC,ZYOVEC,CMO,XINDX,UDV,
     *               WRK(KEVECS),WRK(KWRK1),LWRK1)
      END IF
C
C *** END OF RSPELI
C
      CALL QEXIT('RSPELI')
      RETURN
      END
C  /* Deck rsporb */
      SUBROUTINE RSPORB(ONEIND,NSIM,FC,FCX,FVX,QAX,QBX,UDV, EVECS)
C
C WRITTEN 14-FEB 1986
C
C PURPOSE:
C  1)DISTRIBUTE INACTIVE FCX AND ACTIVE FVX FOCK MATRICES AND QAX AND
C    QBX MATRICES INTO ORBITAL PART OF LINEAR TRANSFORMED VECTOR
C
C                  ( [ E(K,L) , H ] )    K<L
C       E[2]*N = - (                )
C                  ( [ E(L,K) , H ] )    K>L
C    X MAY REFER TO EITHER ONE INDEX TRANSFORMED MATRICES (N IS A
C    ORBITAL TRIAL VECTOR ) OR TRANSITION DENSITY MATRICES (N IS A
C    CONFIGURATION TRIAL VECTOR). FOR CONFIGURATION TRIAL VECTORS
C    FCX IS IDENTICALLY ZERO. FURTHER OVERLAP IS ZERO BECAUSE TRIAL
C    VECTORS ARE CHOSEN ORTHOGONAL TO REFERENCE STATE.
C
C  2)CREATE LINEAR TRANSFORMED VECTOR S[2]*N FOR N EQUAL TO EITHER A
C    ORBITAL OR A CONFIGURATION TRIAL VECTOR
C
C    ONEIND = .TRUE. FOR AN ORBITAL TRIAL VECTOR
C
C                ******************************
C
C    [ E(P,Q) , H ] =
C
C    ACTIVE-INACTIVE (P,Q) = (T,I)
C                             L,K
C    FCX(I,X)*DV(T,X)-2*FCX(I,T)-2FVX(I,T)+QBX(I,T)
C        K       L          K,L       K,L      K,L
C
C    INACTIVE-ACTIVE         (I,U)
C                             K,L
C    -FCX(X,I)*DV(X,U)+2*FCX(U,I)+2*FVX(U,I)-QAX(I,U)
C           K       L        L,K        L,K      K,L
C
C    ACTIVE-ACTIVE           (T,U)
C                             K,L
C                             L,K
C    FCX(U,X)*DV(T,X)-FCX(X,T)*DV(X,U)+QBX(U,T)-QAX(T,U)
C        L       K          K       L      L,K      K,L
C        K       L          L       K      K,L      L,K
C
C    ACTIVE-SECONDARY        (T,A)
C                             K,L
C    FCX(A,X)*DV(T,X)+QBX(A,T)
C        L       K        L,K
C
C    SECONDARY-ACTIVE        (A,U)
C                             L,K
C    -FCX(X,A)*DV(X,U)-QAX(A,U)
C           L       K      L,K
C
C    INACTIVE-SECONDARY      (I,A)
C                             K,L
C    2*FCX(A,I)+2FVX(A,I)
C          L,K       L,K
C
C    SECONDARY-INACTIVE      (A,J)
C                             L,K
C    -2*FCX(J,A)-2*FVX(J,A)
C           K,L        K,L
C
C
#include "implicit.h"
#include "priunit.h"
C
      LOGICAL ONEIND
C
      DIMENSION FC(*),FCX(NORBT,NORBT,*),FVX(NORBT,NORBT,*)
      DIMENSION QAX(NORBT,NASHDI,*),QBX(NORBT,NASHDI,*)
      DIMENSION UDV(NASHDI,NASHDI,*), EVECS(KZYVAR,*)
C
C  INFINP : DOMCSRDFT
C  INFDIM : NASHDI
C
#include "maxorb.h"
#include "maxash.h"
#include "infinp.h"
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
      PARAMETER ( D0 = 0.0D0 , D2 =2.0D0 , DM1 = -1.0D0 )
C
      CALL QENTER('RSPORB')
C
      KYCONF = KZCONF + KZVAR
      if (KYCONF + kzwopt .gt. kzyvar) then
         write (lupri,*) 'hjaaj kzconf, kzvar, kzwopt, kzyvar',
     &                    kzconf, kzvar, kzwopt, kzyvar
         call quit('hjaaj problem in rsporb')
      end if
C
C DISTRIBUTE FOCK MATRICES IN EVECS
C
         K_SYM1 = 0
         DO 1300 IG = 1,KZWOPT
            K     = JWOP(1,IG)
            L     = JWOP(2,IG)
            K_SYM  = ISMO(K)
            L_SYM  = ISMO(L)
            IF( K_SYM.NE.K_SYM1 ) THEN
               K_SYM1 = K_SYM
               NORBK = NORB(K_SYM)
               IORBK = IORB(K_SYM)
               NASHK = NASH(K_SYM)
               NISHK = NISH(K_SYM)
               IASHK = IASH(K_SYM)
               IIORBK= IIORB(K_SYM)
               IORBL = IORB(L_SYM)
               NASHL = NASH(L_SYM)
               NISHL = NISH(L_SYM)
               NORBL = NORB(L_SYM)
               IASHL = IASH(L_SYM)
               IIORBL= IIORB(L_SYM)
            END IF
            NK    = K - IORBK
            NL    = L - IORBL
            ITYPK = IOBTYP(K)
            ITYPL = IOBTYP(L)
            IF ( ITYPK.EQ.JTINAC )THEN
               DO 2000 ISIM = 1 ,NSIM
                  IF (ONEIND) THEN
                     EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                       + D2 * FCX(L,K,ISIM)
                     IF ((.NOT. TDA) .AND. (.NOT. CISRPA))
     *                  EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                       - D2 * FCX(K,L,ISIM)
                     IF ( NASHT . GT . 0 ) THEN
                        EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                       + D2 * FVX(L,K,ISIM)
                     IF ((.NOT. TDA) .AND. (.NOT. CISRPA))
     *                  EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                       - D2 * FVX(K,L,ISIM)
                     END IF
                  ELSE
                     EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                       + D2 * FVX(L,K,ISIM)
                     IF ((.NOT. TDA) .AND. (.NOT. CISRPA))
     *                  EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                       - D2 * FVX(K,L,ISIM)
                  END IF
                  IF ( ITYPL.EQ.JTACT ) THEN
                     NWL = ISW(L) - NISHT
                     IF (ONEIND) THEN
                        EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                  - DDOT(NASHL,FCX(IORBL+NISHL+1,K,ISIM),1,
     *                                        UDV(IASHL+1,NWL,1),1)
                        IF ((.NOT. TDA) .AND. (.NOT. CISRPA)) THEN
                           DO 730 IX = 1,NASHL
                              EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,
     *                                                   ISIM)
     *                        + FCX(K,IORBL+NISHL+IX,ISIM)*
     *                                    UDV(NWL,IASHL+IX,1)
 730                       CONTINUE
                        END IF
                     ELSE
                        TEMP1 = D0
                        TEMP2 = D0
                        NX  = NISHK
                        DO 825 NWX = IASHK+1,IASHK+NASHK
                           NX = NX + 1
                           IF (NX.LE.NK) THEN
                              FCXK = FC(IIORBK+IROW(NK)+NX)
                           ELSE
                              FCXK = FC(IIORBK+IROW(NX)+NK)
                           END IF
                           TEMP1 = TEMP1 + FCXK * UDV(NWX,NWL,ISIM)
                           TEMP2 = TEMP2 + FCXK * UDV(NWL,NWX,ISIM)
 825                    CONTINUE
                        EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                                          - TEMP1
                        IF ((.NOT. TDA) .AND. (.NOT. CISRPA))
     *                     EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                                             + TEMP2
                     END IF
                     EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                     - QAX(K,NWL,ISIM)
                     IF ((.NOT. TDA) .AND. (.NOT. CISRPA))
     *                  EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                        + QBX(K,NWL,ISIM)
                  END IF
 2000          CONTINUE
            ELSE
              IF (ITYPL.EQ.JTACT) THEN
                  NWL = ISW(L) - NISHT
                  NWK = ISW(K) - NISHT
                  DO 2020 ISIM=1,NSIM
                     IF (ONEIND) THEN
                        EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                  -DDOT(NASHL,FCX(IORBL+NISHL+1,K,ISIM),1,
     *                                        UDV(IASHL+1,NWL,1),1)
                        IF ((.NOT. TDA) .AND. (.NOT. CISRPA)) THEN
                           DO IX = 1,NASHL
                              EVECS(KYCONF+IG,ISIM) =
     *                          EVECS(KYCONF+IG,ISIM)
     *                          +FCX(K,IORBL+NISHL+IX,ISIM)*
     *                                          UDV(NWL,IASHL+IX,1)
                           END DO
                           EVECS(KYCONF+IG,ISIM) =
     *                          EVECS(KYCONF+IG,ISIM)
     *                          -DDOT(NASHK,FCX(IORBK+NISHK+1,L,ISIM),
     *                                        1,UDV(IASHK+1,NWK,1),1)
                        END IF
                        DO IX = 1,NASHK
                           EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                        +FCX(L,IORBK+NISHK+IX,ISIM)*
     *                                        UDV(NWK,IASHK+IX,1)
                        END DO
                     ELSE
                        TEMP1 = D0
                        TEMP2 = D0
                        NX  = NISHL
                        DO 835 NWX = IASHL+1,IASHL+NASHL
                           NX = NX + 1
                           IF (NX.LE.NL) THEN
                              FCXL = FC(IIORBL+IROW(NL)+NX)
                           ELSE
                              FCXL = FC(IIORBL+IROW(NX)+NL)
                           END IF
                           TEMP1 = TEMP1 + FCXL * UDV(NWX,NWK,ISIM)
                           TEMP2 = TEMP2 + FCXL * UDV(NWK,NWX,ISIM)
 835                    CONTINUE
                        EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                                          + TEMP2
                        IF ((.NOT. TDA) .AND. (.NOT. CISRPA))
     *                      EVECS(KYCONF+IG,ISIM) =
     *                             EVECS(KYCONF+IG,ISIM) - TEMP1
                        TEMP1 = D0
                        TEMP2 = D0
                        NX  = NISHK
                        DO 845 NWX = IASHK+1,IASHK+NASHK
                           NX = NX + 1
                           IF (NX.LE.NK) THEN
                              FCXK = FC(IIORBK+IROW(NK)+NX)
                           ELSE
                              FCXK = FC(IIORBK+IROW(NX)+NK)
                           END IF
                           TEMP1 = TEMP1 + FCXK * UDV(NWX,NWL,ISIM)
                           TEMP2 = TEMP2 + FCXK * UDV(NWL,NWX,ISIM)
 845                    CONTINUE
                        EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                                          - TEMP1
                        IF ((.NOT. TDA) .AND. (.NOT. CISRPA))
     *                      EVECS(KYCONF+IG,ISIM) =
     *                             EVECS(KYCONF+IG,ISIM) + TEMP2
                     END IF
                     EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                  +QBX(L,NWK,ISIM) -QAX(K,NWL,ISIM)
                     IF ((.NOT. TDA) .AND. (.NOT. CISRPA))
     *                  EVECS(KYCONF+IG,ISIM) =
     *                          EVECS(KYCONF+IG,ISIM)
     *                          +QBX(K,NWL,ISIM) - QAX(L,NWK,ISIM)
 2020             CONTINUE
               ELSE
                  NWK = ISW(K) - NISHT
                  DO 2030 ISIM=1,NSIM
                     IF (ONEIND) THEN
                        IF ((.NOT. TDA) .AND. (.NOT. CISRPA))
     *                      EVECS(KYCONF+IG,ISIM) =
     *                         EVECS(KYCONF+IG,ISIM)
     *                         -DDOT(NASHK,FCX(IORBK+NISHK+1,L,ISIM),1,
     *                                        UDV(IASHK+1,NWK,1),1)
                        DO 760 IX = 1,NASHK
                           EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                        +FCX(L,IORBK+NISHK+IX,ISIM)*
     *                                        UDV(NWK,IASHK+IX,1)
 760                    CONTINUE
                     ELSE
                        TEMP1 = D0
                        TEMP2 = D0
                        NX  = NISHL
                        DO 860 NWX = IASHL+1,IASHL+NASHL
                           NX = NX + 1
                           IF (NX.LE.NL) THEN
                              FCXL = FC(IIORBL+IROW(NL)+NX)
                           ELSE
                              FCXL = FC(IIORBL+IROW(NX)+NL)
                           END IF
                           TEMP2 = TEMP2 + FCXL * UDV(NWX,NWK,ISIM)
                           TEMP1 = TEMP1 + FCXL * UDV(NWK,NWX,ISIM)
 860                    CONTINUE
                        EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                                          + TEMP1
                        IF ((.NOT. TDA) .AND. (.NOT. CISRPA))
     *                      EVECS(KYCONF+IG,ISIM) =
     *                              EVECS(KYCONF+IG,ISIM) - TEMP2
                     END IF
                     IF ((.NOT. TDA) .AND. (.NOT. CISRPA))
     *                  EVECS(KYCONF+IG,ISIM) =
     *                           EVECS(KYCONF+IG,ISIM)
     *                           -QAX(L,NWK,ISIM)
                     EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                     +QBX(L,NWK,ISIM)
 2030             CONTINUE
               ENDIF
            ENDIF
 1300    CONTINUE
C
C CHANGE SIGN ON ORBITAL PART OF LINEAR TRANSFORMATION E[2]*N
C
      DO ISIM = 1,NSIM
         CALL DSCAL(KZWOPT,DM1,EVECS(KZCONF+1,ISIM),1)
         CALL DSCAL(KZWOPT,DM1,EVECS(KYCONF+1,ISIM),1)
      END DO
C
      IF (IPRRSP.GT.110) THEN
         WRITE(LUPRI,*)' LINEAR TRANSFORMATION WITH  E(2)'//
     &      ' (Z1 Y1 Z2 Y2 ...)'
         CALL OUTPUT(EVECS,1,KZVAR,1,2*NSIM,KZVAR,2*NSIM,-1,LUPRI)
      END IF
C
C END OF RSPORB
C
      CALL QEXIT('RSPORB')
      RETURN
      END
C  /* Deck rspnex */
      SUBROUTINE RSPNEX(LINEQ,LMAXIT,NSIM,IBTYP,EIVAL,RESID,EIVEC,GP,
     *                  CMO,UDV,PVX,FC,FV,FCAC,H2AC,XINDX,
     *                  BVECS,LWRK)
C
C PURPOSE: 1)CONSTRUCT RESIDUAL (E(2)-W(I)*S(2))*X(I)
C            AND MODIFIED GRADIENT  (E(2)-W(I)*S(2))*X(I)(csf part)
C            FOR KEXSIM EIGENVECTORS X(I) OF REDUCED RSP PROBLEM
C          2)TEST FOR CONVERGENCE OF KEXCNV EIGENVECTORS,
C            CONVERGENCE CRITERIUM:
C            //(E(2)-W(I)*S(2))*X(I)// .LE. THCRSP * //X(I)//
C          3)GENERATE NEW TRIAL VECTORS
C             A) CSF TRIAL VECTORS, USE GENERALIZED DAVIDSON ALGORITHM
C                CSF TRIAL VECTORS CONTAIN KZCONF ELEMENTS
C             B) ORBITAL TRIAL VECTORS
C                IF OPTORB.EQ.TRUE, USE OPTIMAL ORBITAL ALGORITHM ELSE
C                USE GENERALIZED DAVIDSON ALGORITHM
C
C Parameters: LMAXIT true: max iterations reached, do not
C                          calculate new trial vectors
C
C PJ NOV-1984
C LAST REV 12 NOV 1984, 22 APR 1988, 28 Oct 93 hjaaj
C
#include "implicit.h"
#include "mxcent.h"
#include "dummy.h"
      DIMENSION IBTYP(*),EIVAL(*),RESID(*),EIVEC(KZYRED,*)
      DIMENSION GP(*),BVECS(*)
      DIMENSION CMO(*),UDV(*),PVX(*),FC(*),FV(*),FCAC(*)
      DIMENSION XINDX(*),H2AC(*)
#include "ibndxdef.h"
      PARAMETER ( D0=0.0D0 , D1=1.0D0 , DM1=-1.0D0 )
      PARAMETER ( REDFAC = 0.3D0 , DTEST=1.0D-4  )
C
C Used from common blocks:
C  INFRSP: ?
C  WRKRSP: KZVAR,KZYVAR,KEXSIM,KEXCNV,
C          LURSP5
C  INFOPT: EACTIV
C
#include "priunit.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "inftap.h"
#include "infopt.h"
#include "infdim.h"
#include "abares.h"
#include "thrldp.h"
#include "qm3.h"
C
      LOGICAL LINEQ,OLSEN1,LMAXIT
      real(8), allocatable :: THCORP(:),EOVAL(:)
      integer, allocatable :: KEX(:)
C
      CALL QENTER('RSPNEX')
C
      IF (KEXSIM .LE. 0) THEN
         WRITE (LUPRI,*) 'ERROR in RSPNEX: KEXSIM .le. 0'
         WRITE (LUPRI,*) 'KEXSIM =',KEXSIM
         CALL QUIT('RSPNEX error: KEXSIM .le. 0')
      END IF

      allocate(KEX(KEXSIM))
      if (OPTORB) then
         allocate(EOVAL(NSIM))
         allocate(THCORP(NSIM))
      end if

C
C     KEX(*) is used for defining which trial vectors are
C     configuration type (=1) and which are orbital type (=0).
C     920904-hjaaj: KEX(i) = -1 means this root is converged.
C     RESID(*) is used for saving residuals for final print.
C
      NTEST = 0
      IBOFF = 0
      NTCONV = 0
      DO 4000 ISIMC = 1,KEXSIM,NSIM
         NBX  = MIN(NSIM,(KEXSIM+1-ISIMC))
         NFIN = 0
C
C RSPEVE constructs solution vectors from reduced space
C CONSTRUCT -W(I) * S[2]*X(I) IN FIRST NBX ELEMENTS OF
C BVECS. THE PART WHICH IS USED IN THE OPTIMAL ORBITAL
C TRIAL VECTOR ALGORITHM WHICH IS STORED IN THE CONSECUTIVE
C NBX*KZYWOP ELEMENTS
C
C First check if all frequencies zero (W(i)=0, all i) because then
C we may skip the -W(I) * S[2] term completely
C If SOPPA we will construct the 2p-2h diagonal times the 2p-2h trial
C vector explicitly so we need the solution vectors
C
         IF (SOPPA) GOTO 25
         DO 24 I = ISIMC, ISIMC-1+NBX
            IF ( EIVAL(I) .NE. D0 ) GO TO 25
   24    CONTINUE
         LEN = NBX*KZYVAR
         IF (OPTORB) LEN = LEN + NBX*KZYWOP
         CALL DZERO(BVECS,LEN)
         GO TO 200
C
   25    CALL RSPEVE(IBTYP,EIVAL,EIVEC,BVECS(1+KZYVAR),
     *               BVECS,NBX,(ISIMC-1))
         IF (IPRRSP.GT.101) THEN
            WRITE (LUPRI,*) ' RSPNEX: Solution vectors no.',
     *         ISIMC,' to',(ISIMC-1+NBX),' (Z1 Y1 ...)'
            CALL OUTPUT(BVECS(1+KZYVAR),1,KZVAR,1,2*NBX,
     *                  KZVAR,2*NBX,-1,LUPRI)
         ENDIF
         KRES   = 1
         KINFIX = KRES + (NBX+1)*KZYVAR
         KBORB  = KINFIX + NBX*KZYWOP
         KDIA   = KBORB + KZYWOP
C
C Make room for the D-matrix in SOPPA
C
         IF (SOPPA) THEN
            KWRK1 = KDIA + KZCONF
         ELSE
            KWRK1 = KDIA
         ENDIF
         LWRK1 = LWRK   - KWRK1
         IF (LWRK1.LT.0) CALL ERRWRK('RSPNEX 1',KWRK1-1,LWRK)
         IF (OPTORB) CALL DZERO(BVECS(KINFIX),NBX*KZYWOP)
C
C  Read in the D-matrix in SOPPA from LURSP4
C
         IF (SOPPA) THEN
            REWIND (LURSP4)
            CALL READT(LURSP4,KZCONF,BVECS(KDIA))
         ENDIF
         DO 50 IBX = 1,NBX
            IRESO = (IBX-1)*KZYVAR + 1
            IFIXOF= KINFIX + (IBX-1)*KZYWOP
            WIBX  = -EIVAL(ISIMC-1+IBX)
            CALL DZERO(BVECS(IRESO),KZYVAR)
            IF (WIBX .EQ. D0) GO TO 49
C           ... skip S[2] term if zero frequency W(i) = 0
C           We still need to do the D-matrix term in SOPPA
            IF ( KZCONF .GT. 0 ) THEN
               CALL RSPSLI(1,0,BVECS(IRESO+KZYVAR+KZVAR),VDUMMY,
     *                  UDV,BVECS(IRESO),XINDX,BVECS(KWRK1),LWRK1)
               CALL DSWAP(KZVAR,BVECS(IRESO),1,BVECS(IRESO+KZVAR),1)
               CALL DSCAL(KZYVAR,DM1,BVECS(IRESO),1)
               CALL RSPSLI(1,0,BVECS(IRESO+KZYVAR),VDUMMY,
     *                     UDV,BVECS(IRESO),XINDX,BVECS(KWRK1),LWRK1)
               IF ( OPTORB ) THEN
                  CALL DAXPY(KZWOPT,WIBX,BVECS(IRESO+KZCONF),1,
     *                       BVECS(IFIXOF),1)
                  CALL DZERO(BVECS(IRESO+KZCONF),KZWOPT)
                  CALL DAXPY(KZWOPT,WIBX,BVECS(IRESO+KZVAR+KZCONF),1,
     *                       BVECS(IFIXOF+KZWOPT),1)
                  CALL DZERO(BVECS(IRESO+KZVAR+KZCONF),KZWOPT)
               END IF
            END IF
            IF ( KZWOPT.GT.0 ) THEN
               CALL DCOPY(KZWOPT,BVECS(IRESO+KZYVAR+KZCONF),1,
     *                    BVECS(KBORB),1)
               CALL DCOPY(KZWOPT,BVECS(IRESO+KZYVAR+KZVAR+KZCONF),1,
     *                    BVECS(KBORB+KZWOPT),1)
               CALL RSPSLI(0,1,VDUMMY,BVECS(KBORB),
     *                     UDV,BVECS(IRESO),XINDX,BVECS(KWRK1),LWRK1)
            ENDIF
            CALL DSCAL(KZYVAR,WIBX,BVECS(IRESO),1)
C
 49         IF (SOPPA) THEN
               DO I=0,KZCONF-1
                  BVECS(IRESO+I) = BVECS(IRESO+I) +
     *                             BVECS(KDIA+I) * BVECS(IRESO+KZYVAR+I)
                  BVECS(IRESO+KZVAR+I) = BVECS(IRESO+KZVAR+I) +
     *                       BVECS(KDIA+I) * BVECS(IRESO+KZVAR+KZYVAR+I)
               ENDDO
            ENDIF
 50      CONTINUE
         IF (OPTORB) THEN
            NTT = NBX*KZYWOP
            DO 52 II = 1,NTT
               BVECS(NBX*KZYVAR+II) = BVECS(KINFIX-1+II)
 52         CONTINUE
C           CALL DCOPY(NBX*KZYWOP,BVECS(KINFIX),1,
C    *                            BVECS(1+NBX*KZYVAR),1)
         END IF
C
C ALLOCATE WORK SPACE
C
  200    KBVECS = 1
         KFIXOR = KBVECS + NBX*KZYVAR
         KWRK1  = KFIXOR + NBX*KZYWOP
         LWRK1  = LWRK   - KWRK1
         IF (LWRK1.LT.0) CALL ERRWRK('RSPNEX 2',KWRK1-1,LWRK)
C
C CONSTRUCT  E[2]*X(I) WHERE X(I) IS THE I:th
C EIGENVECTOR AND EIGENVALUE RESPECTIVELY OF THE REDUCED
C RSP PROBLEM
C
         REWIND (LURSP5)
         IF (KOFFTY.EQ.1) READ(LURSP5)
         DO 600 K = 1,KZRED
C
C   If SOPPA and the b-vector is 2p2h then only read the p-h part
C   from LURSP5 since the 2p-2h is constructed explicitly above
C
            IF (SOPPA) THEN
               IF (IBTYP(K+KOFFTY) .EQ. JBCNDX) THEN
                  CALL READT(LURSP5,KZYWOP,BVECS(KWRK1))
               ELSE
                  CALL READT(LURSP5,KZYVAR,BVECS(KWRK1))
               ENDIF
            ELSE
               CALL READT(LURSP5,KZYVAR,BVECS(KWRK1))
            ENDIF
            DO 700 JR = 1,NBX
               JROFF  = (JR-1)*KZYVAR
               JROOTJ = IBOFF+JR
               FAC1   = EIVEC(2*K-1,JROOTJ)
               FAC2   = EIVEC(2*K,JROOTJ)
               IF (SOPPA .AND. IBTYP(K+KOFFTY).EQ.JBCNDX
     *                   .AND. .NOT.OPTORB) THEN
                  CALL DAXPY(KZWOPT,FAC1,BVECS(KWRK1),1,
     *                       BVECS(1+JROFF+KZCONF),1)
                  CALL DAXPY(KZWOPT,FAC1,BVECS(KWRK1+KZWOPT),1,
     *                       BVECS(1+JROFF+KZVAR+KZCONF),1)
                  CALL DAXPY(KZWOPT,FAC2,BVECS(KWRK1+KZWOPT),1,
     *                       BVECS(1+JROFF+KZCONF),1)
                  CALL DAXPY(KZWOPT,FAC2,BVECS(KWRK1),1,
     *                       BVECS(1+JROFF+KZVAR+KZCONF),1)
               ELSE IF ( IBTYP(K+KOFFTY) .EQ. JBONDX
     *                              .OR. .NOT.OPTORB) THEN
                  CALL DAXPY(KZYVAR,FAC1,BVECS(KWRK1),1,
     *                       BVECS(1+JROFF),1)
                  CALL DAXPY(KZVAR,FAC2,BVECS(KWRK1+KZVAR),1,
     *                                       BVECS(1+JROFF),1)
                  CALL DAXPY(KZVAR,FAC2,BVECS(KWRK1),1,
     *                                       BVECS(1+KZVAR+JROFF),1)
               ELSE
C
C CSF PART OF LINEAR TRANSFORMED VECTOR
C
                  IF (.NOT. SOPPA) THEN
                     CALL DAXPY(KZCONF,FAC1,BVECS(KWRK1),1,
     *                          BVECS(1+JROFF),1)
                     CALL DAXPY(KZCONF,FAC1,BVECS(KWRK1+KZVAR),1,
     *                                      BVECS(1+KZVAR+JROFF),1)
                     CALL DAXPY(KZCONF,FAC2,BVECS(KWRK1),1,
     *                                      BVECS(1+KZVAR+JROFF),1)
                     CALL DAXPY(KZCONF,FAC2,BVECS(KWRK1+KZVAR),1,
     *                                      BVECS(1+JROFF),1)
                  ENDIF
C
C                 ORBITAL PART OF LINEAR TRANSFORMED VECTOR
C
                  IF (KZWOPT.GT.0) THEN
                     KKAPZY = KFIXOR + (JR-1)*KZYWOP
                     KKAPYZ = KKAPZY + KZWOPT
                     CALL DAXPY(KZWOPT,FAC1,BVECS(KWRK1+KZCONF),1,
     *                                      BVECS(KKAPZY),1)
                     CALL DAXPY(KZWOPT,FAC1,BVECS(KWRK1+KZVAR+KZCONF),
     *                                      1,BVECS(KKAPYZ),1)
                     CALL DAXPY(KZWOPT,FAC2,BVECS(KWRK1+KZCONF),1,
     *                                      BVECS(KKAPYZ),1)
                     CALL DAXPY(KZWOPT,FAC2,BVECS(KWRK1+KZVAR+KZCONF),
     *                                      1,BVECS(KKAPZY),1)
                  ENDIF
               ENDIF
 700        CONTINUE
 600     CONTINUE
C
         IF (OPTORB) THEN
            KFIXOF = KFIXOR
            DO 1000 JR = 1,NBX
               JRESOF = (JR-1)*KZYVAR + 1
               CALL DAXPY(KZWOPT,D1,BVECS(KFIXOF),1,
     *                              BVECS(JRESOF+KZCONF),1)
               CALL DAXPY(KZWOPT,D1,BVECS(KFIXOF+KZWOPT),1,
     *                              BVECS(JRESOF+KZVAR+KZCONF),1)
               KFIXOF = KFIXOF + KZYWOP
 1000       CONTINUE
         END IF
         IF (LINEQ) THEN
            KFIXOF = KFIXOR
            DO 1050 JR=1,NBX
               JROFF = (JR-1)*KZYVAR
               CALL DAXPY(KZYVAR,DM1,GP,1,
     *                        BVECS(1+JROFF),1)
               IF (OPTORB) THEN
                  CALL DAXPY(KZWOPT,DM1,GP(1+KZCONF),1,
     *                        BVECS(KFIXOF),1)
                  CALL DAXPY(KZWOPT,DM1,GP(1+KZVAR+KZCONF),1,
     *                        BVECS(KFIXOF+KZWOPT),1)
               ENDIF
            KFIXOF = KFIXOF + KZYWOP
 1050       CONTINUE
         ENDIF
         IF (IPRRSP.GT.101) THEN
            WRITE (LUPRI,*) ' RSPNEX: (Z,Y) residual vectors no.',
     *         ISIMC,' to',(ISIMC-1+NBX)
            CALL OUTPUT(BVECS,1,KZVAR,1,2*NBX,
     *                  KZVAR,2*NBX,-1,LUPRI)
         ENDIF
C
C        THE RESIDUAL IS NOW CONSTRUCTED
C
C        TEST FOR CONVERGENCE OF THE EIGENVECTORS
C        AND FORM NEW TRIAL VECTORS
C
         CALL PHPINI(LPHPMX,KZCONF,KZWOPT,MAXPHP)
         KDIAE   = KWRK1
         KDIASO  = KDIAE  + LPHPMX
         IF (OLSEN) THEN
            KCVEC = KDIASO + KZWOPT
            KTOT  = KCVEC  + KZCONF
         ELSE
            KCVEC = KDIASO
            KTOT  = KDIASO + KZWOPT
         END IF
         LTOT  = LWRK  - KTOT
         IF (LTOT.LT.0) CALL ERRWRK('RSPNEX 3',KTOT-1,LWRK)
         ECORE  = -EACTIV
         IF (PHPRES) THEN
            IPWAY = 2
         ELSE IF (KZCONF.GT.0) THEN
            CALL PHPDSK(BVECS(KDIAE),LPHPMX,.FALSE.)
         ELSE
            CALL RSPEDG(BVECS(KDIAE))
         END IF
         CALL RSPSOD(BVECS(KDIASO))
         NCSIM = 0
         NOSIM = 0
         NOTCNV = 0
         KCOFF = 0
         IF (IPRRSP.GE.1 .AND. IPRRSP.LE.5) THEN
            IF (LINEQ) THEN
               WRITE(LUPRI,'(/6X,A)')
     &' No.  Residual tot.,    conf., and orb.    Bnorm      Frequency '
            ELSE
               WRITE(LUPRI,'(/6X,A)')
     &'Root  Residual tot.,    conf., and orb.    Bnorm      Eigenvalue'
            END IF
            WRITE(LUPRI,'(6X,A)')
     &'----------------------------------------------------------------'
         END IF
         DO 2000 JR = 1,NBX
            JROFF = (JR-1)*KZYVAR
            JROOTJ = IBOFF+JR
            QONORM = DDOT(KZWOPT,BVECS(1+KZCONF+JROFF),1,
     *                               BVECS(1+KZCONF+JROFF),1)
     *             + DDOT(KZWOPT,BVECS(1+KZVAR+KZCONF+JROFF),1,
     *                               BVECS(1+KZVAR+KZCONF+JROFF),1)
            QCZNOR = DDOT(KZCONF,BVECS(1+JROFF),1,
     *                               BVECS(1+JROFF),1)
            QCYNOR = DDOT(KZCONF,BVECS(1+KZVAR+JROFF),1,
     *                               BVECS(1+KZVAR+JROFF),1)
            QCNORM = QCZNOR + QCYNOR
            QNORM  = SQRT(QONORM+QCNORM)
            QCNORM = SQRT(QCNORM)
            QONORM = SQRT(QONORM)
            QCZNOR = SQRT(QCZNOR)
            QCYNOR = SQRT(QCYNOR)
            BNORM  = DDOT(KZYRED,EIVEC(1,JROOTJ),1,EIVEC(1,JROOTJ),1)
            BNORM  = SQRT(BNORM)
            IF (IPRRSP.GE.6) THEN
               WRITE(LUPRI,'(/A,I3,1P,3(A,D12.5),/T10,2(A,D12.5))')
     *           ' ROOT:',JROOTJ,' RESIDUAL TOT:',QNORM,
     *           ' CONF:',QCNORM,' ORB:',QONORM,
     *           ' BNORM:',BNORM,' EIVAL:',EIVAL(JROOTJ)
               IF (QCNORM .NE. D0) WRITE(LUPRI,'(1P,T10,2(A,D12.5))')
     *           ' RESIDUAL CONF Z-COMP',QCZNOR,'    Y-COMP',QCYNOR
            ELSE IF (IPRRSP.GE.1) THEN
               WRITE(LUPRI,'(I10,1P,D15.5,3D10.2,D15.5)')
     *           JROOTJ,QNORM,QCNORM,QONORM,BNORM,EIVAL(JROOTJ)
            END IF
CHJAAJ      BNORM = MIN(BNORM,D1)
C           ... 880309: ask for absolute convergence if BNORM
C               is greater than 1, otherwise relative convergence.
C               930910-hjaaj: NO, gives problems with number of digits
C               when e.g. BNORM = 1.0D6 for a FC term (remember E[2]
C               is asymmetric of the order of the mcscf gradient).
            QBTEST = QNORM/BNORM
            RESID(JROOTJ) = QBTEST
            IF (QBTEST .LE. THCRSP) THEN
C           ... this root is converged
               KEX(JR) = -1
            ELSE
C           ... this root is not converged
               IF (JROOTJ .LE. KEXCNV) NOTCNV = NOTCNV + 1
               EXCITA  =  EIVAL(JROOTJ)
              IF (.NOT.LMAXIT) THEN
               IF (QCNORM.GT.QONORM) THEN
                  NCSIM = NCSIM + 1
                  KEX(IBOFF+JR) = 1
C
C              USE DAVIDSON ALGORITHM TO FORM NEW TRIAL VECTORS
C
                  IF ( QCZNOR .GT. QCYNOR ) THEN
                     IOFF  = 0
                     EXITA =   EXCITA
                  ELSE
                     IOFF  = KZVAR
                     EXITA =  -EXCITA
                  END IF
                  IF ((.NOT.LINEQ).AND.OLSEN) THEN
                     OLSEN1 = .TRUE.
                  ELSE
                     OLSEN1 = .FALSE.
                  ENDIF
                  IF (OLSEN1) THEN
                     CALL RSPCVE(IOFF,JROOTJ,IBTYP,EIVAL,EIVEC,
     *                  BVECS(KCVEC),BVECS(KTOT),LTOT)
                  END IF
                  IF (PHPRES) THEN
                     DO 1003 J = 1,KZCONF
                        BVECS(KDIAE-1+J) = BVECS(JROFF+IOFF+J) /
     *                                     BVECS(KDIAE-1+J)
 1003                CONTINUE
                     CALL PHPGET(KSYMST,KZCONF,XINDX,FCAC,
     *                           H2AC,BVECS(KDIAE),
     *                           ECORE,IPWAY,IPRRSP,BVECS(KTOT),LTOT)
                     CALL RSPEDG(BVECS(KDIAE))
                  END IF
                  CALL NEXCI(OLSEN1,EXITA,KZCONF,BVECS(KCOFF+1),
     *               BVECS(KCVEC),BVECS(JROFF+IOFF+1),
     *               BVECS(KDIAE),IPRRSP,BVECS(KTOT),LTOT)
C                 CALL NEXCI(OLSEN,ENER,NCVAR,D,  XVEC,RES,
C    &                 DIAG,IPRPHP,WRK,LWRK)
                  KCOFF = KCOFF + KZCONF
               ELSE
                  NOSIM = NOSIM + 1
                  KEX(IBOFF+JR) = 0
                  IF (PHPRES) CALL RSPEDG(BVECS(KDIAE))
                  IF (OPTORB) THEN
C
C        USE OPTIMAL ORBITAL ALGORITHM TO FORM NEW TRIAL VECTORS
C        WILL BE DONE IN LATER SECTION
C
                     IF (QCNORM .GT. THCRSP) THEN
                        THCORP(NOSIM) = REDFAC * QCNORM/BNORM
                     ELSE
C                       ... CI part is converged, attempt
C                           quadratic convergence in orbital part.
                        THCORP(NOSIM) = REDFAC *
     &                     MAX(THCRSP,QBTEST*QBTEST)
                     END IF
                     EOVAL(NOSIM)  = EIVAL(JROOTJ)
                  ELSE
C
C        USE DAVIDSON ALGORITHM TO FORM NEW TRIAL VECTORS
C
                     KSJR = KFIXOR - 1 + (JR-1)*KZYWOP
                     DO I=1,KZWOPT
                        IZ = I+KZCONF

                        ZDIA = BVECS(KDIAE-1+IZ)
     &                       - EXCITA*BVECS(KDIASO-1+I)
                        IF (ABS(ZDIA).LT.DTEST) ZDIA = SIGN(DTEST,ZDIA)
                        BVECS(KSJR+I) = BVECS(IZ+JROFF) / ZDIA

                        YDIA = BVECS(KDIAE-1+IZ)
     &                       + EXCITA*BVECS(KDIASO-1+I)
                        IF (ABS(YDIA).LT.DTEST) YDIA = SIGN(DTEST,YDIA)
                        BVECS(KSJR+KZWOPT+I) =
     &                         BVECS(IZ+KZVAR+JROFF) / YDIA
                     END DO

                     IF (TDA .OR. CISRPA) THEN
                        ! Only Z component allowed, otherwise the removal of Y component in RSPOR will not work correctly
                        IF (LINEQ) THEN
                           ! Copy Y component to Z component if Y component has largest element,
                           ! this is necessary to converge frequency dependent properties
                           ! (for excitation energies, just zero Y component).
                           IB_MAX = IDAMAX(KZYWOP,BVECS(KSJR),1)
                           IF ( IB_MAX .GT. KZWOPT ) THEN
                              CALL DCOPY(KZWOPT,BVECS(KSJR+KZWOPT),1,
     &                                          BVECS(KSJR),1)
                           END IF
                        END IF
                        ! zero Y component
                        CALL DZERO(BVECS(KSJR+KZWOPT+1),KZWOPT)
                     END IF ! (TDA .OR. CISRPA)

                  ENDIF ! (OPTORB) ... ELSE ...
               ENDIF ! (QCNORM.GT.QONORM) ... ELSE ...
              ENDIF ! .NOT. LMAXIT
            ENDIF ! (QBTEST .LE. THCRSP) ... ELSE ...
 2000    CONTINUE
         NTCONV = NTCONV + NOTCNV

         IF (MMPCM) THEN
             IF (NSTATES.LE.0) THEN
               CALL QUIT('FATAL ERROR: NSTATES.LE.0 for MMPCM')
            ENDIF
            DO KL=1,NBX
               ICQM3(KL) = KEX(KL)
            END DO
         ENDIF

         IF (NTCONV.EQ.0 .OR. NCSIM+NOSIM.EQ.0) GO TO 3999
C        ... NTCONV is number of non-converged response vectors
C
      IF ((IPRRSP.GT.110).AND.(NCSIM.GT.0)) THEN
         WRITE(LUPRI,*)' NCSIM CSF TRIAL VECTORS',NCSIM
         CALL OUTPUT(BVECS,1,KZCONF,1,NCSIM,KZCONF,NCSIM,-1,LUPRI)
      ENDIF
      IF ((IPRRSP.GT.110).AND.(NOSIM.GT.0)) THEN
         WRITE(LUPRI,*)' NBX ORBITAL VECTORS WRITTEN OUT',NBX
         WRITE(LUPRI,*)' NOSIM ORBITAL TRIAL VECTORS    ',NOSIM
         WRITE(LUPRI,*)' KZWOPT                         ',KZWOPT
         WRITE(LUPRI,*)' (Z1 Y1 Z2 Y2 ...)'
         CALL OUTPUT(BVECS(KFIXOR),1,KZWOPT,1,2*NBX,
     &               KZWOPT,2*NBX,-1,LUPRI)
      ENDIF
         KCVEC = 1
         KOVEC = KCVEC + NCSIM * KZCONF
         KWRK1 = KOVEC + NOSIM * KZYWOP
         LWRK1 = LWRK  - KWRK1
         IF (LWRK1.LT.0) CALL ERRWRK('RSPNEX 4',KWRK1-1,LWRK)
         KOOFF = KOVEC
         DO 2010 JR = 1,NBX
            IF (KEX(IBOFF+JR) .EQ. 0) THEN
               DO 2011 II = 1,KZYWOP
                  BVECS(KOOFF-1+II) = BVECS(KFIXOR+(JR-1)*KZYWOP-1+II)
 2011          CONTINUE
C              CALL DCOPY(KZYWOP,BVECS(KFIXOR+(JR-1)*KZYWOP),1,
C    *                           BVECS(KOOFF),1)
               KOOFF = KOOFF + KZYWOP
            ENDIF
            IF (JR.LE.NCSIM) THEN
               IBTYP(KOFFTY+KZRED+NTEST+JR) = JBCNDX
            ELSE
               IBTYP(KOFFTY+KZRED+NTEST+JR) = JBONDX
            ENDIF
 2010    CONTINUE
         IF (OPTORB.AND.(NOSIM.GT.0)) THEN
C
C GET DIAGONAL ORBITAL PART OF E(2) AND S(2)
C
            CALL DCOPY(KZWOPT,BVECS(KDIAE+KZCONF),1,BVECS(KTOT),1)
            CALL DCOPY(KZWOPT,BVECS(KDIAE+KZCONF),1,
     *                        BVECS(KTOT+KZWOPT),1)
            CALL DCOPY(KZWOPT,BVECS(KDIASO),1,BVECS(KTOT+KZYWOP),1)
            CALL DCOPY(KZWOPT,BVECS(KDIASO),1,
     *                        BVECS(KTOT+KZYWOP+KZWOPT),1)
            CALL DSCAL(KZWOPT,DM1,BVECS(KTOT+KZYWOP+KZWOPT),1)
            LWRKPA = LTOT - 2*KZYWOP
            CALL ORPPAR(NOSIM,THCORP,EOVAL,IBTYP,BVECS(KOVEC),
     *                  BVECS(KTOT), BVECS(KTOT+KZYWOP),
     *                  CMO,UDV,PVX,FC,FV,FCAC,
     *                  XINDX,BVECS(KTOT+2*KZYWOP),LWRKPA)
C
C           CALL ORPPAR(NOSIM,THCORP,EOVAL,IBTYP,A1,ORBDIE,ORBDIS,
C    *                  CMO,UDV,PVX,FC,FV,FCAC,XINDX,WRK,LWRK)
C
         ENDIF
C
C        ORTHOGONALIZE TRIAL VECTORS AND EXAMINE FOR LINEAR DEPENDENCE
C
         NFIN = NCSIM + NOSIM
         NBPREV_here = KZRED + NTEST + KOFFTY
         CALL RSPORT(BVECS,NFIN,NBPREV_here,IBTYP,
     *               THRLDP,BVECS(KTOT),LURSP3)
C        CALL RSPORT(BVECS,NBX,NBPREV,IBTYP,THRLDP,OLDVEC,LU3)
C
C        NUMBER OF NEW TRIAL VECTORS (NFIN is updated by RSPORT)
C
      IF (IPRRSP.GT.105) THEN
         NCSIM = 0
         NOSIM = 0
         NTSIM = KZRED + KOFFTY + NTEST
         DO 3010 IX = 1,NFIN
            IF (IBTYP(NTSIM+IX).EQ.JBCNDX) THEN
               NCSIM = NCSIM + 1
            ELSE
               NOSIM = NOSIM + 1
            ENDIF
 3010    CONTINUE
      ENDIF
      IF ((IPRRSP.GT.105).AND.(NCSIM.GT.0)) THEN
         WRITE(LUPRI,*)' NCSIM CSF trial vectors after RSPORT',NCSIM
         CALL OUTPUT(BVECS,1,KZCONF,1,NCSIM,KZCONF,NCSIM,-1,LUPRI)
      ENDIF
      IF ((IPRRSP.GT.105).AND.(NOSIM.GT.0)) THEN
         WRITE(LUPRI,'(/A,I5/A)')
     +      ' NOSIM orbital trial vectors after RSPORT',NOSIM,
     +      ' ( columns: Z_1 Y_1 Z_2 Y_2 ... Z_NOSIM Y_NOSIM)'
         CALL OUTPUT(BVECS(1+NCSIM*KZCONF),1,KZWOPT,
     *             1,NOSIM*2,KZWOPT,NOSIM*2,-1,LUPRI)
      ENDIF
C
 3999 CONTINUE
         IF (IPRRSP.GT.5) THEN
            WRITE(LUPRI,'(/A,I5,A,I5)')
     *      ' NUMBER OF TRIAL VECTORS IN THIS LOAD',NBX,
     *      '   OFFSET FOR THIS LOAD',IBOFF
            WRITE(LUPRI,'(I5,A,I5)') NOTCNV,
     *      ' NON-CONVERGED SOLUTION VECTORS THIS LOAD, TOTAL:',NTCONV
            WRITE(LUPRI,'(I5,A,I5)') NFIN,
     *      ' LINEAR INDEPENDENT TRIAL VECTORS ADDED TO PREVIOUS',NTEST
         END IF
         NTEST  = NTEST + NFIN
         IBOFF  = IBOFF + NBX
         IF (KZYRED + 2*NTEST .GT. MAXRM) GO TO 4010
C        ... test if maximum dimension of reduced space exceeded
 4000 CONTINUE
 4010 CONTINUE
      JEXSIM = NTEST
C
C JEXSIM: THE ACTUAL NUMBER OF NEW TRIAL VECTORS IN THE
C        NEXT MICROITERATION
C
C ---
C
C Output:
C
      IF (NTCONV.EQ.0) THEN
C
C        EIGENVECTORS CONVERGED
C
         IF (IPRRSP .GE. 0) WRITE (LUPRI,5030) KEXCNV
 5030    FORMAT(/' *** THE REQUESTED',I5,' SOLUTION VECTORS CONVERGED')
         KCONV=1
      ELSE IF (NTEST.EQ.0) THEN
C
C        LINEAR DEPENDENCE BETWEEN TRIAL VECTORS
C
         IF (LMAXIT) THEN
            WRITE(LUPRI,5020)
            KCONV=0
         ELSE
            WRITE(LUPRI,5010)
            KCONV=-1
         END IF
 5010    FORMAT(/' *** MICROITERATIONS STOPPED DUE TO LINEAR',
     *           ' DEPENDENCE BETWEEN NEW TRIAL VECTORS')
 5020    FORMAT(/' *** MICROITERATIONS STOPPED DUE TO MAX',
     *           ' ITERATIONS REACHED.')
         WRITE (LUPRI,5031) KEXCNV
      ELSE IF (KZYRED + 2*NTEST .GT. MAXRM) THEN
C
C        MAXIMUM DIMENSION OF REDUCED SPACE EXCEEDED
C
         WRITE (LUPRI,'(/A/A,I5,A/)')
     *   ' *** MICROITERATIONS STOPPED BEFORE CONVERGENCE BECAUSE',
     *   '     MAXIMUM DIMENSION OF REDUCED SPACE',MAXRM,' EXCEEDED.'
         WRITE (LUPRI,5031) KEXCNV
         KCONV = -2
      ELSE
C
C        MICROITERATIONS NOT CONVERGED
C
         KCONV=0
      ENDIF
 5031 FORMAT(/' *** WARNING : REQUESTED',I5,
     &        ' SOLUTION VECTORS NOT CONVERGED')
      IF (LMAXIT .OR. KCONV .LT. 0 .OR. IPRRSP .GE. 11
     &    .OR. (IPRRSP .GE. 2 .AND. KCONV .EQ. 1) ) THEN
         WRITE (LUPRI,96) THCRSP,KZYRED,(I,RESID(I),I=1,KEXCNV)
   96    FORMAT(/' Convergence of RSP solution vectors,',
     &           ' threshold =',1P,D9.2,
     *          /' ------------------------------------',
     &           '---------------------------',
     *          /' (dimension of paired reduced space:',I5,')',
     *          /,(' RSP solution vector no.',I5,'; norm of residual',
     *          1P,D11.2))
      END IF
C
C     Transfer information about the convergence back to ABACUS
C
      IF (KCONV .NE. 1) THEN
         FINRES = RESID(1)
         NOCONV = .TRUE.
         KZRED_ABA = KZRED
      ELSE
         FINRES = RESID(1)
         NOCONV = .FALSE.
         KZRED_ABA = KZRED
      END IF
C
C     END OF RSPNEX.
C
      deallocate(KEX)
      if (OPTORB) then
         deallocate(EOVAL)
         deallocate(THCORP)
      end if
      CALL QEXIT('RSPNEX')
      RETURN
      END
C  /* Deck rspord */
      SUBROUTINE RSPORD(SRED,EIVEC,ALFAR,ALFAI,BETA,WRK1,WRK2,ISNDX)
C
C Analyze and order eigenvectors from RSPRED
C
C      JO 1984-11-8
C
#include "implicit.h"
      DIMENSION SRED(KZYRED,KZYRED), EIVEC(KZYRED,KZYRED)
      DIMENSION WRK1(KZYRED,KZYRED), WRK2(KZYRED,KZYRED)
      DIMENSION ALFAR(*),ALFAI(*),BETA(*), ISNDX(3)
C
C eigenvalue k is (ALFAR(k)+i*ALFAI(k))/BETA(k)
C
C Used from common blocks:
C  /WRKRSP/: KZYRED
C
#include "priunit.h"
#include "infrsp.h"
#include "wrkrsp.h"
C
      PARAMETER (ZEROT=1.0D-9)
      PARAMETER (COMPLX = 1.0D7 )
      PARAMETER (D0=0.0D0, D1=1.0D0)
C
      CALL QENTER('RSPORD')
C
      IF (IPRRSP .GT. 20) THEN
         WRITE (LUPRI,'(//A/A)')
     *      '     (ALFAR + i ALFAI) / BETA are eigenvalues;',
     *      '   I      ALFAR           ALFAI          BETA'
         WRITE (LUPRI,'(1P,I5,3D15.6)')
     *      (I,ALFAR(I),ALFAI(I),BETA(I),I=1,KZYRED)
      END IF
      DO 10 I=1,KZYRED
         IF (ABS(BETA(I)).GE.ZEROT) THEN
            ALFAR(I)=ALFAR(I)/BETA(I)
            ALFAI(I)=ALFAI(I)/BETA(I)
         ELSE
C           singularities
            WRITE(LUPRI,1010)KSYMOP,KZYRED,I,ALFAR(I),ALFAI(I),BETA(I)
 1010       FORMAT(/' *** WARNING Singularity in reduced',
     +          ' response equation of symmetry',I3,
     +         /' *** WARNING KZYRED,I,ALFAR(I),ALFAI(I),BETA(I):',
     +          2I4,1P,3D15.6)
         END IF
C        IF(ABS(ALFAR(I)).GT.ZEROT) THEN
C    +   .OR.(ABS(ALFAR(I)).LE.ZEROT.AND.ABS(ALFAI(I)).GT.ZCRIT)) THEN
C           complex eigenvalues
C        RATIO= ABS(ALFAI(I)/ALFAR(I))
C        IF(RATIO.GT.ZCRIT) WRITE(LUPRI,1020) ALFAR(I),ALFAI(I)
         IF(ABS(ALFAI(I)).GT.D0) THEN
            WRITE(LUPRI,1020) KSYMOP,KZYRED,I,ALFAR(I),ALFAI(I)
 1020       FORMAT(/' *** WARNING Complex eigenvalue in reduced',
     +          ' response eq. of symmetry',I3,
     +         /' *** WARNING real and imaginary part :',2I4,1P,2D15.6)
C
C SET EIGENVALUE EQUAL TO COMPLX IN ORDER TO BE ABLE TO SKIP
C CONTRIBUTIONS FROM THIS ROOT WHEN SUMMING UP TERMS
C IN THE CALCULATION OF THE EFFECTIVE SPECTRUM IN C6 CALCULATIONS
C
            ALFAR(I) = COMPLX
         END IF
  10  CONTINUE
      IF (IPRRSP .GT. 15) THEN
         WRITE (LUPRI,'(/A)')
     &      ' Unsorted eigenvalues of [ E(2) , S(2) ] :'
         WRITE (LUPRI,'(I10,1P,D15.6)') (I,ALFAR(I),I=1,KZYRED)
      END IF
C
C     reduced S(2) in eigenvector basis
C
      CALL DGEMM('N','N',KZYRED,KZYRED,KZYRED,1.D0,
     &           SRED,KZYRED,
     &           EIVEC,KZYRED,0.D0,
     &           WRK1,KZYRED)
      CALL DGEMM('T','N',KZYRED,KZYRED,KZYRED,1.D0,
     &           EIVEC,KZYRED,
     &           WRK1,KZYRED,0.D0,
     &           WRK2,KZYRED)
      IF (IPRRSP .GT. 20) THEN
         WRITE (LUPRI,'(/A)') ' Reduced S(2) in eigenvector basis :'
         IF (IPRRSP .GE. 25) THEN
            CALL OUTPUT(WRK2,1,KZYRED,1,KZYRED,KZYRED,KZYRED,-1,LUPRI)
         ELSE
            WRITE (LUPRI,'(I10,1P,D20.4)') (I,WRK2(I,I),I=1,KZYRED)
         END IF
      END IF
C
C     Select eigenvectors with positive normalization and store them
C        as the first eigenvectors.
C
      IPOS = 0
      IZER = 0
      INEG = 0
      DO 20 I=1,KZYRED
         IF (BETA(I) .LE. ZEROT) THEN
            IZER  = IZER + 1
         ELSE IF (WRK2(I,I) .GT. D0) THEN
            IPOS  = IPOS + 1
            SCALE = D1/SQRT(WRK2(I,I))
            CALL DSCAL(KZYRED,SCALE,EIVEC(1,I),1)
            IF (IPOS.NE.I) THEN
               CALL DSWAP(KZYRED,EIVEC(1,I),1,EIVEC(1,IPOS),1)
               XSAVE           = ALFAR(IPOS)
               ALFAR(IPOS)     = ALFAR(I)
               ALFAR(I)        = XSAVE
               XSAVE           = WRK2(IPOS,IPOS)
               WRK2(IPOS,IPOS) = WRK2(I,I)
               WRK2(I,I)       = XSAVE
               XSAVE           = BETA(IPOS)
               BETA(IPOS)      = BETA(I)
               BETA(I)         = XSAVE
            END IF
         ELSE IF (WRK2(I,I) .LT. D0) THEN
            INEG  = INEG + 1
            SCALE = D1/SQRT(ABS(WRK2(I,I)))
            CALL DSCAL(KZYRED,SCALE,EIVEC(1,I),1)
            IF (ALFAR(I) .EQ. COMPLX) ALFAR(I) = -COMPLX
         END IF
   20 CONTINUE
      NNEG = IPOS
      DO 22 I=IPOS+1,KZYRED
         IF (ABS(BETA(I)) .GT. ZEROT) THEN
            NNEG = NNEG + 1
            IF (NNEG.NE.I) THEN
               CALL DSWAP(KZYRED,EIVEC(1,I),1,EIVEC(1,NNEG),1)
               XSAVE           = ALFAR(NNEG)
               ALFAR(NNEG)     = ALFAR(I)
               ALFAR(I)        = XSAVE
               XSAVE           = WRK2(NNEG,NNEG)
               WRK2(NNEG,NNEG) = WRK2(I,I)
               WRK2(I,I)       = XSAVE
            END IF
         END IF
   22 CONTINUE
      ISNDX(1) = IPOS
      ISNDX(2) = IZER
      ISNDX(3) = INEG
      IF (IPOS.NE.(KZYRED/2) .OR. IZER.NE.0)
     &   WRITE (LUPRI,2020) IPOS,IZER,INEG
 2020 FORMAT(/' *** EIGENVECTORS WITH POSITIVE METRIC:',I6,
     *       /' *** EIGENVECTORS WITH ZERO METRIC:    ',I6,
     *       /' *** EIGENVECTORS WITH NEGATIVE METRIC:',I6)
C
C Order eigensolutions in ascending order of eigenvalues.
C
      DO 100 I=1,IPOS
         JMIN = I
         AMIN = ALFAR(I)
         DO, J=I+1,IPOS
            IF(ALFAR(J).LT.AMIN) THEN
               AMIN = ALFAR(J)
               JMIN = J
            ENDIF
         END DO
         IF (JMIN.NE.I) THEN
            ALFAR(JMIN)=ALFAR(I)
            ALFAR(I)=AMIN
            CALL DSWAP(KZYRED,EIVEC(1,I),1,EIVEC(1,JMIN),1)
         ENDIF
  100 CONTINUE
      DO 110 I=IPOS+1,IPOS+INEG
         JMIN = I
         AMIN = ALFAR(I)
         DO, J=I+1,IPOS+INEG
            IF(ALFAR(J).GT.AMIN) THEN
               AMIN = ALFAR(J)
               JMIN = J
            ENDIF
         END DO
         IF (JMIN.NE.I) THEN
            ALFAR(JMIN)=ALFAR(I)
            ALFAR(I)=AMIN
            CALL DSWAP(KZYRED,EIVEC(1,I),1,EIVEC(1,JMIN),1)
         ENDIF
  110 CONTINUE
      IF (INEG.NE.IPOS) THEN
         WRITE(LUPRI,'(/A/A,I5)')' **WARNING**'//
     &    ' number of eigenvalues with negative metric differ from'//
     &    ' number with positive metric.',
     &    ' Dimension of reduced matrix =', KZYRED
         IF (IPRRSP.GT.20) THEN
            WRITE(LUPRI,'(/A)') '   NUMBER    EIGENVALUE '
            WRITE(LUPRI,'(I10,1P,D20.8)') (I,ALFAR(I),I=1,KZYRED)
         END IF
      ELSE
         IF (IPRRSP.GT.20) WRITE(LUPRI,'(/A)')
     &      '      NUMBER      EIGENVALUE      PAIRED EIGENVALUE'
         DO, I=1,IPOS
            RELDIF = ABS(ABS(ALFAR(I))-ABS(ALFAR(IPOS+I)))
            IF (ABS(ALFAR(I)) .GT. ZEROT)
     &         RELDIF = RELDIF / ABS(ALFAR(I))
            IF (RELDIF .GT. ZEROT) THEN
               J = 1
               WRITE(LUPRI,'(/A,I5)')
     &         ' **WARNING**  EIGENVALUES NOT PAIRED.'//
     &         ' Dimension of reduced matrix =', KZYRED
            ELSE
               J = 0
            END IF
            IF (IPRRSP.GT.20 .OR. J.NE.0) WRITE(LUPRI,'(I10,1P,2D20.8)')
     &         I,ALFAR(I),ALFAR(IPOS+I)
         END DO
         IF (IZER .GT. 0) THEN
            WRITE(LUPRI,'(/A/A,I5/A)')
     &         ' WARNING: ZERO METRIC EIGENVALUE AND ZERO EIGENVALUE.',
     &         ' Dimension of reduced matrix =', KZYRED,
     &         '     NUMBER       EIGENVALUE'
            WRITE(LUPRI,'(I10,1P,D20.8)')(I,ALFAR(IPOS+INEG+I),I=1,IZER)
         END IF
         CALL FLSHFO(LUPRI)
      END IF
C
      CALL QEXIT('RSPORD')
      RETURN
      END
C  /* Deck rsport */
      SUBROUTINE RSPORT (BVECS,NBNEW,NBPREV,IBTYP,THRLDP,OLDVEC,LU3)
C
C Written 25-Oct-1984 by Hans Jorgen Aa. Jensen
C (based on ORTVEC from MSCSF section)
C REVISED: 22-Jan-1986 hjaaj (Z = Y or Z = -Y check).
C           8-Jun-1986 hjaaj (pure ci or orbital trial vectors, IBTYP)
C          17-Jun-1987 hjaaj (fixed orthogonalization when NBNEW .gt. 1)
C           5-MAY-1988 PJ (CONFIGURATION TRIAL VECTORS HAVE ONLY
C                          KZCONF VARIABLE )
C          21-Jul-1992 HH Do supersymmetry averaging if necessary.
C          12-Jul-1993 hjaaj (ZEQLY = 1.0D-6 instead of 1.0D-10;
C              discard trial vector for linear dep. if OVLPI close to 1;
C              do symmetric orthonormalization twice)
C           8-Feb-1995 hjaaj (changed ZEQLY to OVLMIN; improved
C              code for Z = Y or Z = -Y tests).
C
C Purpose:
C    Orthogonalize new b-vectors against all previous b-vectors
C  and among themselves, and renormalize.
C    Each b-vector is in (Z, Y) form, the (Y, Z) vector to be obtained
C  by switching first and second half. Each of the new b-vectors are
C  orthogonalized both against (Z, Y) and (Y, Z) for each previous
C  b-vector.
C  (Orthogonalization is performed twice if round-off is large,
C   if larger than THRRND).
C
C Input:
C  BVECS,  non-orthogonal b-vectors
C  NBNEW,  number of b-vectors in BVECS
C  NBPREV, number of previous b-vectors on LU3
C  IBTYP,  b-vector type (CI or obital)
C  THRLDP, threshold for linear dependence
C
C Output:
C  BVECS,  orthogonal b-vectors (also written to LU3)
C  NBNEW,  number of b-vectors written to LU3
C
C Scratch:
C  OLDVEC(KZYVAR), scratch array for old b-vectors on LU3
C
#include "implicit.h"
#include "dummy.h"
      DIMENSION BVECS(*), IBTYP(*), OLDVEC(*)
C
#include "ibndxdef.h"
C
C    Change value of ZEQLY for problem with linear equation
C    95-02-08
C     PARAMETER (ZEQLY = 1.D-6, T1MIN = 1.D-8, THRRND = 1.D-4)
C
      PARAMETER (OVLMIN = 1.D-4, T1MIN = 1.D-8, THRRND = 1.D-4)
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, D2 = 2.0D0)
C
C Used from common blocks:
C  INFRSP: RSPSUP
C  WRKRSP: KZCONF,KZWOPT,KSYMOP
C
#include "priunit.h"
#include "infrsp.h"
#include "wrkrsp.h"
C
C     Local arrays:
C
      integer, allocatable :: JLDP(:)
C
      IF (NBNEW .LE. 0) RETURN

      allocate( JLDP(NBNEW) )
C
      CALL QENTER('RSPORT')
C
C *** Average new trial vectors for supersymmetry
C
      IF ( RSPSUP .AND. (KSYMOP .EQ. 1)) THEN
         IBVEC = 1
         DO 2600 IRX = 1, NBNEW
            IF (IBTYP (NBPREV+IRX) .EQ. JBCNDX) THEN
C           -- Configuration trial vector; shift the starting point
               IBVEC = IBVEC + KZCONF
            ELSE IF (IBTYP(NBPREV+IRX) .EQ. JBONDX) THEN
C           -- Orbital trial vector; do the averaging
               CALL RSPAVE(BVECS(IBVEC),KZWOPT,2)
               IBVEC = IBVEC + KZYWOP
            END IF
 2600    CONTINUE
      END IF
C
C
C     For zero frequency linear response, the result will
C     be Z = Y for real perturbations and Z = -Y for imaginary
C     perturbations; trial vectors will have the same structure.
C
C
C     However, if Z = Y or Z = -Y then (Z Y) and (Y Z) are linear
C     dependent and we want instead (Z 0) and (0 Z) as trial vectors.
C
      IBVEC = 1
      NCNEW = 0
      NONEW = 0
      DO 1400 IRX = 1,NBNEW
         JLDP(IRX) = 1 ! initialize to not linear dependent
         IF (IBTYP(NBPREV+IRX) .EQ. JBCNDX) THEN
            NCNEW = NCNEW + 1
            IBVEC = IBVEC + KZCONF
         ELSE IF (IBTYP(NBPREV+IRX) .EQ. JBONDX) THEN
            NONEW = NONEW + 1
            T1    = DDOT(KZYWOP,BVECS(IBVEC),1,BVECS(IBVEC),1)
            JBVEC = IBVEC + KZWOPT
            OVLPI = D2*DDOT(KZWOPT,BVECS(IBVEC),1,BVECS(JBVEC),1)
C           test if Z + Y is zero
            TZPY  = T1 + OVLPI
            IF (IPRRSP .GE. 10) WRITE (LUPRI,'(/A/,(A,1P,D15.8))')
     *         ' ** TEST OF Z = Y OR Z = -Y:',
     *         ' (Z Y) norm squared =',T1,
     *         ' Z + Y norm squared =',TZPY
            IF (TZPY .LE. OVLMIN*T1) THEN
               IF (IPRRSP .GE. 10) WRITE (LUPRI,'(/A,I4/A)')
     *            ' Z = -Y in trial vector no.',IRX,
     *            ' Y component removed.'
               CALL DZERO(BVECS(JBVEC),KZWOPT)
            ELSE IF (TZPY .GT. T1) THEN
C              test if Z - Y is zero;
C              only if TZPY = 2 * T1 can Z - Y be zero,
C              so TZPY .gt. T1 is a safe test.
               TZMY = T1 - OVLPI
               IF (IPRRSP .GE. 10) WRITE (LUPRI,'(A,1P,D15.8)')
     *            ' Z - Y norm squared =',TZMY
               IF (TZMY .LE. OVLMIN*T1) THEN
                  IF (IPRRSP .GE. 10) WRITE (LUPRI,'(/A,I4/A)')
     *            ' Z = Y in trial vector no.',IRX,
     *            ' Y component removed.'
                  CALL DZERO(BVECS(JBVEC),KZWOPT)
               END IF
            END IF
            IBVEC = IBVEC + KZYWOP
         ELSE
            GO TO 9100
         END IF
 1400 CONTINUE
C
      IROUND = 0
      ITURN  = 0
 1500 ITURN  = ITURN + 1
C
C        Orthogonalize new b-vectors agains previous b-vectors
C        (both (Z, Y) and (Y, Z))
C
         REWIND (LU3)
         DO 2000 K = 1,NBPREV
            IBTYPK = IBTYP(K)
            IF (IBTYPK .EQ. JBCNDX) THEN
               IF (NCNEW .EQ. 0) THEN
                  READ (LU3)
                  GO TO 2000
C        v------------------
               END IF
               LBZ = KZCONF
               LBZY= KZCONF
               LBSKIP = 2*KZWOPT
            ELSE IF (IBTYPK .EQ. JBONDX) THEN
               IF (NONEW .EQ. 0) THEN
                  READ (LU3)
                  GO TO 2000
C        v------------------
               END IF
               LBZ = KZWOPT
               LBZY= 2*KZWOPT
               LBSKIP = KZCONF
            ELSE
               GO TO 9100
            END IF
            CALL READT(LU3,LBZY,OLDVEC)
            IBVEC = 1
            DO 1900 IRX = 1,NBNEW
               IF (IBTYP(NBPREV+IRX) .NE. IBTYPK) THEN
                  IBVEC = IBVEC + LBSKIP
               ELSE IF (JLDP(IRX).EQ.0) THEN
                  IBVEC = IBVEC + LBZY
               ELSE
C                 Orthogonalize to (Z Y)
                  TT = -DDOT(LBZY,   OLDVEC,1,BVECS(IBVEC),1)
                  CALL DAXPY(LBZY,TT,OLDVEC,1,BVECS(IBVEC),1)
C                 Orthogonalize to (Y Z)
                  IF ( IBTYPK .EQ. JBONDX ) THEN
                     JBVEC = IBVEC + LBZ
                     TT = -DDOT(LBZ,   OLDVEC,       1,BVECS(JBVEC),1)
     *                    -DDOT(LBZ,   OLDVEC(1+LBZ),1,BVECS(IBVEC),1)
                     CALL DAXPY(LBZ,TT,OLDVEC(1+LBZ),1,BVECS(IBVEC),1)
                     CALL DAXPY(LBZ,TT,OLDVEC,       1,BVECS(JBVEC),1)
                  END IF
                  IBVEC = IBVEC + LBZY
               END IF
 1900       CONTINUE
 2000    CONTINUE
C
C        Orthogonalize new vectors against each other
C
         IBVEC  = 1
         DO 2400 IRX = 1,NBNEW
            IBTYPI = IBTYP(NBPREV+IRX)
            IF (IBTYPI .EQ. JBCNDX) THEN
               LBZ = KZCONF
               LBZY= KZCONF
               LBSKIP = 2*KZWOPT
            ELSE
               LBZ = KZWOPT
               LBZY= 2*KZWOPT
               LBSKIP = KZCONF
            END IF
         IF (JLDP(IRX) .NE. 0) THEN
            JBVEC = 1
            DO 2300 JRX = 1,(IRX-1)
            IF (IBTYP(NBPREV+JRX) .NE. IBTYPI) THEN
               JBVEC = JBVEC + LBSKIP
            ELSE IF (JLDP(JRX).EQ.0) THEN
               JBVEC = JBVEC + LBZY
            ELSE
               T1 =  DDOT(LBZY,BVECS(JBVEC),1,BVECS(JBVEC),1)
               TT = -DDOT(LBZY,BVECS(JBVEC),1,BVECS(IBVEC),1)/T1
               CALL DAXPY(LBZY,TT,BVECS(JBVEC),1,BVECS(IBVEC),1)
               IF (IBTYPI.EQ.JBONDX) THEN
                  TT=(-DDOT(LBZ,BVECS(JBVEC+LBZ),1,BVECS(IBVEC),1)
     *                -DDOT(LBZ,BVECS(JBVEC),1,BVECS(IBVEC+LBZ),1))/T1
                  CALL DAXPY(LBZ,TT,BVECS(JBVEC+LBZ),1,BVECS(IBVEC),1)
                  CALL DAXPY(LBZ,TT,BVECS(JBVEC),1,BVECS(IBVEC+LBZ),1)
               END IF
               JBVEC = JBVEC + LBZY
            END IF
 2300       CONTINUE
         END IF
C
C NORMALIZE VECTOR NUMBER IRX
C
         IF (JLDP(IRX) .NE. 0) THEN
            TT = DDOT(LBZY,BVECS(IBVEC),1,BVECS(IBVEC),1)
            IF (TT .LE. THRLDP) THEN
               IF (IPRRSP.GT.2) WRITE (LUPRI,3250) IRX,TT
               JLDP(IRX) = 0
            ELSE IF (TT .LT. THRRND) THEN
               IF (ITURN .EQ. 1) THEN
                  IROUND = IROUND + 1
               ELSE
                  IF (IPRRSP.GT.5) WRITE (LUPRI,3255) IRX,TT
                  JLDP(IRX) = 0
               END IF
            END IF
            IF (JLDP(IRX) .NE. 0) THEN
               IF (IPRRSP .GE. 8) WRITE (LUPRI,'(/A,I6,1P,D12.5)')
     &          ' RSPORT: b-vector, norm**2:',IRX,TT
               IF (TT .LT. T1MIN) THEN
                  TT = D1 / SQRT(TT)
                  CALL DSCAL(LBZY,TT,BVECS(IBVEC),1)
                  TT = DDOT(LBZY,BVECS(IBVEC),1,BVECS(IBVEC),1)
               END IF
               TT = D1 / SQRT(TT)
               CALL DSCAL(LBZY,TT,BVECS(IBVEC),1)
            END IF
         END IF
C
C Perform symmetric orthonormalization of (Z Y) and (Y Z) pair
C for vector number IRX
C
C   -1/2       ( C1   C2 )               (  1     OVLPI )
C  S      =   (           ) where S  =  (                )
C              ( C2   C1 )               ( OVLPI     1  )
C
         IF ((JLDP(IRX) .NE. 0).AND.(IBTYPI.EQ.JBONDX)) THEN
            JBVEC = IBVEC + LBZ
            JTURN = 0
 2350       JTURN = JTURN + 1
            OVLPI = D2*DDOT(LBZ,BVECS(IBVEC),1,BVECS(JBVEC),1)
            IF (IPRRSP .GE. 22) THEN
               X1 = DNRM2(LBZY,BVECS(IBVEC),1)
               WRITE (LUPRI,'(A,1P,D14.6,D19.10)')
     *         ' RSPORT, norm and <ZY/YZ> overlap before S-1/2',X1,OVLPI
            END IF
            X2 = ABS(OVLPI)
            IF (ABS(D1-X2) .LT. OVLMIN) THEN
               IF (IPRRSP.GT.5) WRITE (LUPRI,3260) IRX,OVLPI
               JLDP(IRX) = 0
               GO TO 2390
            END IF
            X1 = D1+OVLPI
            X2 = D1-OVLPI
            IF (ABS(X1 - X2) .GT. THRLDP) THEN
               X1 = DP5 / SQRT(X1)
               X2 = DP5 / SQRT(X2)
               C1 = X1 + X2
               C2 = X1 - X2
               ! write(lupri,*) 'X1, X2, C1, C2',X1,X2,C1,C1
               CALL DCOPY(LBZY,BVECS(IBVEC),1,OLDVEC,1)
               CALL DSCAL(LBZY,C1,BVECS(IBVEC),1)
               CALL DAXPY(LBZ,C2,OLDVEC,1,BVECS(JBVEC),1)
               CALL DAXPY(LBZ,C2,OLDVEC(1+LBZ),1,BVECS(IBVEC),1)
            IF (JTURN .EQ. 1) GO TO 2350
C           ...930712-hjaaj: an output with LBZY only = 480 showed
C              that OVLPI ca. 1.0D-14 after first S-1/2;
C              but ca. 1.0D-16 after second S-1/2
               IF (IPRRSP .GE. 22) THEN
                  X1 = DNRM2(LBZY,BVECS(IBVEC),1)
                  OVLPI= D2*DDOT(LBZ,BVECS(IBVEC),1,BVECS(JBVEC),1)
                  WRITE (LUPRI,'(A,1P,D14.6,D19.10)')
     *            ' RSPORT, norm and <ZY/YZ> overlap after  S-1/2',
     *            X1,OVLPI
               END IF
            END IF
         END IF
C
C
 2390    IBVEC = IBVEC + LBZY
 2400    CONTINUE
         IF (IROUND.GT.0) THEN
            IROUND = 0
            IF (ITURN.EQ.1) GO TO 1500
C     ^-------------------------------
            CALL QUIT('RSPORT sw-error; IROUND .gt. 0 and ITURN .gt. 1')
         END IF
 3250 FORMAT(/' (RSPORT) b-vector no.',I3,
     &        ' is removed because of linear dependence;',1P,D12.5)
 3255 FORMAT(/' (RSPORT) b-vector no.',I3,' is removed because of small'
     &        ' norm**2 in 2. Gram-Schmidt orthonormalization:',
     &        1P,D12.5)
 3260 FORMAT(/' (RSPORT) b-vector no.',I3,' is removed because of'
     &        ' linear dependence, <Z Y|Y Z> overlap =',
     &        1P,D14.6)
C
C Save new b-vectors after the old ones on LU3
C
      IBVECX = 1
      IBVEC  = 1
      IBPREV = NBPREV
      DO 3300 IRX = 1,NBNEW
         IBTYPI = IBTYP(NBPREV+IRX)
         IBTYP(NBPREV+IRX) = 0
         IF (IBTYPI .EQ. JBCNDX) THEN
            LBZY= KZCONF
         ELSE
            LBZY= 2*KZWOPT
         END IF
         IF (JLDP(IRX) .NE. 0) THEN
            IBPREV = IBPREV + 1
            IBTYP(IBPREV) = IBTYPI
C
            CALL WRITT(LU3,LBZY,BVECS(IBVEC))
            IF (IBVECX .NE. IBVEC) THEN
               DO 5000 I = 1,LBZY
                  BVECS(IBVECX-1+I) = BVECS(IBVEC-1+I)
 5000          CONTINUE
C    *         CALL DCOPY(LBZY,BVECS(IBVEC),1,BVECS(IBVECX),1)
            END IF
            IBVECX = IBVECX + LBZY
         END IF
         IBVEC = IBVEC + LBZY
 3300 CONTINUE
C
C     Set NBNEW to actual number of new trial vectors
C
      IF (IBPREV-NBPREV .NE. NBNEW) THEN
         WRITE(LUPRI,'(A,I5,A,I5,A)')
     &      'RSPORT:',NBNEW - IBPREV + NBPREV,' out of',
     &       NBNEW,' new trial vectors linear dependent'
      END IF
      NBNEW = IBPREV - NBPREV
C
      deallocate( JLDP )

      CALL QEXIT('RSPORT')
      RETURN
C
C     ERROR : Illegal IBTYP found:
C
 9100 CONTINUE
         WRITE (LUPRI,'(//A,/,(5X,20I3))')
     *      ' RSPORT, illegal IBTYP found, dump of IBTYP:',
     *      (IBTYP(I),I=1,(NBPREV+NBNEW))
         CALL QTRACE(LUERR)
         CALL QUIT('RSPORT ERROR, illegal IBTYP found.')
C
C *** End of subroutine RSPORT
C
      END
C  /* Deck rspprc */
      SUBROUTINE RSPPRC (RSPVEC,NZCONF,IOFFY,IUNITX)
C
C Written 27-Feb-1988 hjaaj
C
C Purpose:
C  Write Z and Y parts of response configuration operator
C  array RSPVEC to unit IUNIT.
C
C  IUNIT = IABS(IUNITX)
C
C  If IUNITX .lt. 0, write all elements;
C  if IUNITX .gt. 0, write only elements with absolut value
C                    .ge. PRFAC * (max. absolut value)
C
#include "implicit.h"
      DIMENSION RSPVEC(*)
      PARAMETER (PRFAC = 1.D-1)
C
C Used from common blocks:
C   INFVAR : JWOPSY
! infrsp.h : TDHF
C
#include "maxorb.h"
#include "infvar.h"
#include "infrsp.h"
C
      IUNIT = IABS(IUNITX)
C
C Any elements?
C
      IF (NZCONF.LE.0) GO TO 300
C
C 1. find maximum absolute value in RSPVEC
C
      IF (IUNITX.GT.0) THEN
         MAXX   = IDAMAX(NZCONF,RSPVEC,1)
         MAXY   = IDAMAX(NZCONF,RSPVEC(IOFFY+1),1)
         THPLIM = PRFAC * MAX( ABS(RSPVEC(MAXX)), ABS(RSPVEC(MAXY)) )
C        WRITE (IUNIT,1000) JWOPSY,PRFAC*100.
      ELSE
C        WRITE (IUNIT,1010) JWOPSY
      ENDIF
      WRITE (IUNIT,1020)
      INPR = 0
      DO 200 IG = 1,NZCONF
         IF (IUNITX.LT.0 .OR. ABS(RSPVEC(IG))       .GT. THPLIM
     *        .OR. ABS(RSPVEC(IOFFY+IG)) .GT. THPLIM) THEN
            WRITE (IUNIT,1100) IG,RSPVEC(IG),RSPVEC(IOFFY+IG)
         ELSE
            INPR = INPR + 1
         ENDIF
  200 CONTINUE
      IF (INPR.GT.0) WRITE (IUNIT,1200) INPR,THPLIM
C
      RETURN
C
  300 CONTINUE
      IF (.NOT.TDHF) WRITE (IUNIT,3000) JWOPSY
      RETURN
C
C1000 FORMAT(/5X,'Configuration operator symmetry =',I2,/,
C    *  5X,'(only elements abs greater than',F8.2,
C    *     ' % of max abs value)')
C1010 FORMAT(/5X,'Response configuration operator symmetry =',I2)
 1020 FORMAT(/,5X,' Index(conf)                         Z operator',
     *            '          Y operator',
     *       /,5X,' ----------                      --------------',
     *            '      --------------')
 1100 FORMAT(I12,20X,2F20.10)
 1200 FORMAT(/I9,' elements with absolute value less than',
     *       1P,D9.2,' not printed.',
     *       //' The numbers in parenthesis',
     *         ' give the configuration symmetry.')
 3000 FORMAT(/5X,'Configuration operator symmetry =',I2,
     *       /5X,'>> NO ELEMENTS <<')
C
C -- end of subroutine RSPPRC
C
      END
C  /* Deck rspanc */
      SUBROUTINE RSPANC (RSPVEC,NZCONF,IOFFY,ISYMV,XNDXCI,SYMPRO,IUNITX)
C
C Based on RSPPRC
C Analyse CI-vector with ANACIN  - ov 920320
C
C
C Purpose:
C  Write Z and Y parts of response configuration operator
C  array RSPVEC to unit IUNIT.
C
C  IUNIT = IABS(IUNITX)
C
C  If IUNITX .lt. 0, write all elements;
C  if IUNITX .gt. 0, write only elements with absolut value
C                    .ge. PRFAC * (max. absolut value)
C
#include "implicit.h"
      DIMENSION RSPVEC(*),XNDXCI(*),SYMPRO(*)
      PARAMETER (PRFAC = 1.D-1, PRSOP1 = 0.3D0, PRSOP2 = 0.25D0)
C
C Used from common blocks:
C   INFINP : FLAG(27), LSYM
C   INFVAR : JWOPSY
C   INFRSP : SOPPA, TDHF
C   INFSOP : KABSAD,KABTAD,......
C
#include "maxorb.h"
#include "infinp.h"
#include "infvar.h"
#include "infrsp.h"
#include "infsop.h"
C
      CHARACTER*1 STTYP
C
      IUNIT = IABS(IUNITX)
C
C Any elements?
C
      IF (NZCONF.LE.0) GO TO 300
C
C 1. find maximum absolute value in RSPVEC
C
      IF (IUNITX.GT.0) THEN
         MAXX   = IDAMAX(NZCONF,RSPVEC,1)
         MAXY   = IDAMAX(NZCONF,RSPVEC(IOFFY+1),1)
         IF (SOPPA) THEN
            THPLIM = PRSOP1 * MAX( ABS(RSPVEC(MAXX)), ABS(RSPVEC(MAXY)))
            THPSOP = PRSOP2 * MAX( ABS(RSPVEC(MAXX)), ABS(RSPVEC(MAXY)))
         ELSE
            THPLIM = PRFAC * MAX( ABS(RSPVEC(MAXX)), ABS(RSPVEC(MAXY)) )
         END IF
C        WRITE (IUNIT,1000) JWOPSY,PRFAC*100.
      ELSE
         THPLIM = 0.0D0
C        WRITE (IUNIT,1010) JWOPSY
      ENDIF
      IF (SOPPA) THEN
         WRITE(IUNIT,1011)
         CALL SOPIJ(ISYMV,IOFFY,IUNIT,RSPVEC,XNDXCI(KABSAD),
     *              XNDXCI(KABTAD),XNDXCI(KIJSAD),XNDXCI(KIJTAD),
     *              XNDXCI(KIJ1AD),XNDXCI(KIJ2AD),XNDXCI(KIJ3AD),THPSOP)
         IF (IPRRSP.GT.20) WRITE(IUNIT,1021)
      ELSE
         WRITE(IUNIT,1020)
      ENDIF
      INPR = 0
      DO 200 IG = 1,NZCONF
         IF (ABS(RSPVEC(IG)) .GT. THPLIM .OR.
     &       ABS(RSPVEC(IOFFY+IG)) .GT. THPLIM) THEN
            IF (SOPPA) THEN
               IF (IPRRSP .GT. 20) THEN
                  CALL SOPIJAB(IG,ISYMV,II,IISYM,IJ,
     &                         IJSYM,IA,IASYM,IB,IBSYM,STTYP,
     &                         XNDXCI(KABSAD),XNDXCI(KABTAD),
     &                         XNDXCI(KIJSAD),XNDXCI(KIJTAD),
     &                         XNDXCI(KIJ1AD),XNDXCI(KIJ2AD),
     &                         XNDXCI(KIJ3AD))
                  WRITE (IUNIT,1110) IG,II,IISYM,IJ,IJSYM,IA,IASYM,IB,
     &                IBSYM,STTYP,RSPVEC(IG),RSPVEC(IOFFY+IG)
               ENDIF
            ELSE
               WRITE (IUNIT,1100) IG,RSPVEC(IG),RSPVEC(IOFFY+IG)
            END IF
         ELSE
            INPR = INPR + 1
         ENDIF
  200 CONTINUE
      IF (INPR.GT.0 .AND. .NOT.SOPPA) WRITE (IUNIT,1200) INPR,THPLIM
      IF (INPR.GT.0 .AND. SOPPA .AND. IPRRSP.GT.20)
     &                                WRITE (IUNIT,1210) INPR,THPLIM
      IF (.NOT. SOPPA) THEN
C        980820-hjaaj: Find out if vector is in CSFs or in det.s
         ICSF = 0
         IF (.NOT. FLAG(27) .AND. .NOT.TRPFLG .AND. ISYMV .EQ. LSYM)
     &      ICSF = 1
         CALL ANACIN(RSPVEC,ISYMV,THPLIM,100,XNDXCI,SYMPRO,IUNITX,ICSF)
      END IF
C
      RETURN
C
  300 CONTINUE
      IF (.NOT.TDHF) WRITE (IUNIT,3000) JWOPSY
      RETURN
C
C1000 FORMAT(/5X,'Configuration operator symmetry =',I2,/,
C    *  5X,'(only elements abs greater than',F8.2,
C    *     ' % of max abs value)')
C1010 FORMAT(/5X,'Response configuration operator symmetry =',I2)
 1011 FORMAT(/5X,'The total double excitation out of the pair i,j')
 1020 FORMAT(/,5X,' Index(conf)                         Z operator',
     *            '          Y operator',
     *       /,5X,' ----------                      --------------',
     *            '      --------------')
 1021 FORMAT(/' Index(2p2h)    i      j      a      b   Type',
     *            '    Z operator       Y operator',
     *       /' -----------   ---    ---    ---    ---  ----',
     *            '  --------------   --------------')
 1100 FORMAT(I12,20X,2F20.10)
 1110 FORMAT(I10,I6,'(',I1,')',I4,'(',I1,')',I4,'(',I1,')',I4,'(',
     *       I1,')',A3,2F17.10)
 1200 FORMAT(/I9,' elements with absolute value ≤',
     *       1P,D9.2,' are not printed.',
     *       //' The numbers in parenthesis',
     *         ' give the configuration symmetry.')
 1210 FORMAT(/I9,' elements with absolute value ≤',
     *       1P,D9.2,' are not printed.',
     *       //' The numbers in parenthesis give the symmetry of the ',
     *         '2p-2h elements.')
 3000 FORMAT(/5X,'Configuration operator symmetry =',I2,
     *       /5X,'>> NO ELEMENTS <<')
C
C -- end of subroutine RSPANC
C
      END
C  /* Deck rsppro */
      SUBROUTINE RSPPRO (WOP,IOFFY,UDV,IUNITX)
C
C Written 31-May-1987 hjaaj
C
C Purpose:
C  Write (rs) and (sr) parts (Z and Y parts) of
C  response orbital operator array WOP to unit IUNIT.
C
C  IUNIT = IABS(IUNITX)
C
C  If IUNITX .lt. 0, write all elements;
C  if IUNITX .gt. 0, write only elements with absolut value
C                    .ge. PRFAC * (max. absolut value)
C
#include "implicit.h"
      PARAMETER (PRFAC = 1.D-1)
      DIMENSION WOP(*), UDV(NASHDI,*)
      REAL*8  :: WOP_SCALED(NWOPT,2)   ! allocatable array
C
C Used from common blocks:
C   INFIND : ISMO(*)
C   INFVAR : NWOPT,JWOPSY,JWOP(2,*)
C   infrsp.h : RSPCI
C
#include "maxorb.h"
#include "maxash.h"
#include "infdim.h"
#include "infind.h"
#include "infvar.h"
#include "infrsp.h"
C
      IUNIT = IABS(IUNITX)
C
C Any elements?
C
      IF (NWOPT.LE.0) THEN
         IF (.NOT.RSPCI) WRITE (IUNIT,3000) JWOPSY
         GO TO 300
      END IF
C
C Scale vector with square root of metric diagonal elements ( S[2] )
C
      WOP_SC_AMAX = 0.0D0
      DO IG = 1, NWOPT
         K  = JWOP(1,IG)
         L  = JWOP(2,IG)
         KA = ICH(K)
         IF (KA .LT. 0) THEN
            OCC_K = 2.0D0
         ELSE IF (KA .GT. 0) THEN
            OCC_K = UDV(KA,KA)
         ELSE
            OCC_K = 0.0D0
         END IF
         LA = ICH(L)
         IF (LA .LT. 0) THEN
            OCC_L = 2.0D0
         ELSE IF (LA .GT. 0) THEN
            OCC_L = UDV(LA,LA)
         ELSE
            OCC_L = 0.0D0
         END IF
         FACTOR = ABS(OCC_K - OCC_L)
         FACTOR = SQRT( FACTOR )
         WOP_SCALED(IG,1) = FACTOR*WOP(IG)
         WOP_SCALED(IG,2) = FACTOR*WOP(IOFFY+IG)
         WOP_SC_AMAX = MAX(WOP_SC_AMAX, ABS(WOP_SCALED(IG,1)),
     &                                  ABS(WOP_SCALED(IG,2)) )
      END DO
C
C 1. find maximum absolute value in scaled WOP
C
      IF (IUNITX.GT.0) THEN
         THPWOP = PRFAC * WOP_SC_AMAX
         WRITE (IUNIT,1000) JWOPSY,PRFAC*100.D0
      ELSE
         THPWOP = 0.0D0
         WRITE (IUNIT,1010) JWOPSY
      ENDIF
      WRITE (IUNIT,1020)
      INPR = 0
      DO 200 IG = 1,NWOPT
         IF (ABS(WOP_SCALED(IG,1)) .GT. THPWOP .OR.
     &       ABS(WOP_SCALED(IG,2)) .GT. THPWOP) THEN
            K  = JWOP(1,IG)
            L  = JWOP(2,IG)
            KS = ISMO(K)
            LS = ISMO(L)
            WRITE (IUNIT,1100) IG,K,KS,L,LS,WOP(IG),WOP(IOFFY+IG),
     &         WOP_SCALED(IG,1),WOP_SCALED(IG,2)
         ELSE
            INPR = INPR + 1
         ENDIF
  200 CONTINUE
      IF (INPR.GT.0) WRITE (IUNIT,1200) INPR,THPWOP
C
  300 CONTINUE
      RETURN
C
 1000 FORMAT(/5X,'Response orbital operator symmetry =',I2,
     &       /5X,'(only scaled elements abs greater than',F8.2,
     &           ' % of max abs value)')
 1010 FORMAT(/5X,'Response orbital operator symmetry =',I2,
     &       /5X,'(all non-zero values)')
 1020 FORMAT(/5X,' Index(r,s)      r      s  ',
     &           '      (r s) operator      (s r) operator',
     &           '      (r s) scaled        (s r) scaled',
     &       /5X,' ----------    -----  -----',
     &           '      --------------      --------------',
     &           '      --------------      --------------')
 1100 FORMAT(I12,I10,'(',I1,')',I4,'(',I1,')',4F20.10)
 1200 FORMAT(/I9,' elements with absolute value ≤',
     *       1P,D9.2,' not printed.',
     *       //' The numbers in parenthesis give the orbital symmetry.')
 3000 FORMAT(/5X,'Orbital operator symmetry =',I2,
     *       /5X,'>> NO ELEMENTS <<')
C
C -- end of subroutine RSPPRO
C
      END
C  /* Deck rspprw */
      SUBROUTINE RSPPRW (WOP,MJWOP,NZWOPT,NWOPSY,IOFFY,IUNITX)
C
C Written Jan 92, pj and ov, based on RSPPRO
C Modified Mar 920314 ov, symmetry in parameter list (more RSPPRO-like)
C
C Purpose:
C  Write (rs) and (sr) parts (Z and Y parts) of
C  response orbital operator array WOP to unit IUNIT.
C
C  IUNIT = IABS(IUNITX)
C
C  If IUNITX .lt. 0, write all elements;
C  if IUNITX .gt. 0, write only elements with absolut value
C                    .ge. PRFAC * (max. absolut value)
C
! infrsp.h : RSPCI
#include "implicit.h"
#include "maxorb.h"
#include "maxash.h"
#include "infdim.h"
#include "infind.h"
#include "infvar.h"
#include "infrsp.h"
      DIMENSION WOP(*), MJWOP(2,MAXWOP,8)
      PARAMETER (PRFAC = 1.D-1)
C
C Used from common blocks:
C   INFIND : ISMO(*)
C
      IUNIT = IABS(IUNITX)
C
C Any elements?
C
      IF (NZWOPT.LE.0) GO TO 300
C
C 1. find maximum absolute value in WOP
C
      IF (IUNITX.GT.0) THEN
         MAXX   = IDAMAX(NZWOPT,WOP,1)
         MAXY   = IDAMAX(NZWOPT,WOP(IOFFY+1),1)
         THPWOP = PRFAC * MAX( ABS(WOP(MAXX)) , ABS(WOP(MAXY)) )
         WRITE (IUNIT,1000) NWOPSY,PRFAC*100.D0
      ELSE
         WRITE (IUNIT,1010) NWOPSY
      ENDIF
      WRITE (IUNIT,1020)
      INPR = 0
      DO 200 IG = 1,NZWOPT
         IF (IUNITX.LT.0 .OR. ABS(WOP(IG))       .GT. THPWOP
     *                   .OR. ABS(WOP(IOFFY+IG)) .GT. THPWOP) THEN
            K  = MJWOP(1,IG,NWOPSY)
            L  = MJWOP(2,IG,NWOPSY)
            KS = ISMO(K)
            LS = ISMO(L)
            WRITE (IUNIT,1100) IG,K,KS,L,LS,WOP(IG),WOP(IOFFY+IG)
         ELSE
            INPR = INPR + 1
         ENDIF
  200 CONTINUE
      IF (INPR.GT.0) WRITE (IUNIT,1200) INPR,THPWOP
C
      RETURN
C
  300 CONTINUE
      IF (.NOT.RSPCI) WRITE (IUNIT,3000) NWOPSY
      RETURN
C
 1000 FORMAT(/5X,'Orbital operator symmetry =',I2,/,
     *  5X,'(only elements abs greater than',F8.2,
     *     ' % of max abs value)')
 1010 FORMAT(/5X,'Response orbital operator symmetry =',I2)
 1020 FORMAT(/,5X,' Index(r,s)      r      s        (r s) operator',
     *            '      (s r) operator',
     *       /,5X,' ----------    -----  -----      --------------',
     *            '      --------------')
 1100 FORMAT(I12,I10,'(',I1,')',I4,'(',I1,')',2F20.10)
 1200 FORMAT(/I9,' elements with absolute value less than',
     *       1P,D9.2,' not printed.',
     *       //' The numbers in parenthesis give the orbital symmetry.')
 3000 FORMAT(/5X,'Orbital operator symmetry =',I2,
     *       /5X,'>> NO ELEMENTS <<')
C
C -- end of subroutine RSPPRW
C
      END
C  /* Deck rspred */
      SUBROUTINE RSPRED (ICTL,LINEQ,N,IBTYP,GP,REDGP,REDE,REDS,
     *                   EIVAL,EIVEC,BCVEC,BOVEC,UDV,EVECS,
     *                   XINDX,WRK,LWRK)
C
C Written 30-oct-1984 by Jeppe Olsen.
C Modified 6-june-1986 to handle split orbital and csf
C trial vectors
C
C PURPOSE:
C Solve mcrsp/rpa problem in reduced subspace.
C The structue of E[2] and S[2] is used to obtain
C a reduced problem with pairing of eigensolutions.
C Although only KZRED vectors are known, a reduced
C problem with KZRED*2 vectors is used.
C
C Input:
C  ICTL, flow control
C         =1, extend the reduced E2 and S2 matrices
C         =2, find the new eigenvalues and -vectors
C         =3, both
C  LINEQ, true for linear equations, false for eigenvalues
C  N, number of new b-vectors
C  IBTYP, types of the N new b-vectors (in BCVEC or BOVEC)
C  REDE,REDS, the old reduced matrices (dimension: KZYRED-2*N)
C  BCVEC,BOVEC, the N new b-vectors
C  EVECS, E[2] times the N b-vectors
C  S[2] times the N b-vectors are created in RSPSLI and
C    stored in EVECS after EVECS is used
C  GP, the property gradient (when LINEQ true)
C  REDGP, the old reduced GP vector (dimension: KZYRED-2*N) when LINEQ
C  UDV, square (unpacked) active density matrix
C  XINDX, information for CI module
C
C
C Output:
C  REDE, the new, extended reduced E[2]-matrix (dimension: KZYRED)
C  REDS, rhe new, extended reduced S[2]-matrix (dimension: KZYRED)
C  REDGP, rhe new, extended reduced GP vector (dimension: KZYRED) when LINEQ true
C  EIVEC, eigenvectors of the new reduced mcrsp/rpa eq.
C  EVALR, real part of eigenvalues of the new reduced  mcrsp/rpa eq.
C
C Scratch:
C  WRK, real scratch array dimension kzyvar
C  BCVEC ( in second part) dimension 3*kzyred**2+2*kzyred
C
#include "implicit.h"
#include "dummy.h"
      DIMENSION IBTYP(*),GP(*),REDGP(*),REDE(*),REDS(*),
     &          UDV(*),BOVEC(*),BCVEC(*),EVECS(KZYVAR,*)
      DIMENSION EIVEC(KZYRED,*),EIVAL(*),XINDX(*),WRK(*)
C
#include "ibndxdef.h"
C
      PARAMETER (D0 = 0.0D0 , D1 = 1.0D0 , DRTEST = 1.0D-5)
C
C Used from common blocks:
C  /WRKRSP/: KZVAR,KZYRED,KZRED,KZYVAR,?
C  /INFIND/: IROW(*). ?
C  /INFPRI/: ??
C
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "inftap.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "infind.h"
#include "infpri.h"
#include "absorp.h"
      INTEGER, ALLOCATABLE ::  IVENU(:)
C
      LOGICAL LINEQ
C
      DIMENSION ISNDX(3),DET(2)
C
      CALL QENTER('RSPRED')
      IF (KZYRED.GT.MAXRM) THEN
         WRITE(LUPRI,'(//A/A,I5,/A,I5)')
     *   ' ERROR IN RSPRED ---',
     *   ' DIMENSION OF REDUCED SPACE IS  ',KZYRED,
     *   ' WHICH EXCEEDS ALLOWED DIMENSION',MAXRM
         CALL QUIT('RSPRED: TOO LARGE DIMENSION OF REDUCED SPACE')
      END IF
CSPAS:01.05.2009 bug fix for large number of  excitation energies
      IF (KZYRED.GT.LIROW) THEN
         WRITE(LUPRI,'(//A/A,A,I5,A,I5,/A)')
     &   ' ERROR IN RSPRED ---',
     &   ' DIMENSION OF INDEX VECTOR IROW IS TO SMALL; ',
     &   ' NEEDED ',KZYRED,', GIVEN :',LIROW,
     &   ' CHANGE DIMENSION IN infind.h AND RECOMPILE'
         CALL QUIT('RSPRED: TOO SMALL DIMENSION OF IROW: RECOMPILE')
      END IF
CKeinSPASmehr
C
C Section 1: extend reduced E[2]-AND S[2]-matrices
C            with N new b-vectors
C
      IF (ICTL.LT.1 .OR. ICTL.GT.3) GO TO 2000
      IF (ICTL.EQ.2) GO TO 1000
      IF (N.LT.1) GO TO 1000
      allocate( IVENU(1:N) )
      IZRED = KZRED - N
C
C Memory check
C
      LNEED = MAX(KZYWOP,KZCONF,NCREF)
C     ... CREF needed in RSPSLI (+ more)
      IF (LNEED .GT. LWRK) CALL ERRWRK('RSPRED',-LNEED,LWRK)
C
C CREATE POINTER VECTOR FOR E(2) LINEAR TRANSFORMED VECTORS
C IN CORE
C
      NCTOT = 0
      NOTOT = 0
      DO INUM = 1,N
         IF ( IBTYP(KOFFTY+IZRED+INUM) .EQ. JBONDX ) THEN
            NOTOT = NOTOT + 1
         ELSE
            NCTOT = NCTOT + 1
         END IF
      END DO
      ICTOT = 0
      IOTOT = 0
      DO INUM = 1,N
         IF ( IBTYP(KOFFTY+IZRED+INUM) .EQ. JBONDX ) THEN
            IOTOT = IOTOT + 1
            IVENU(INUM) = NCTOT + IOTOT
         ELSE
            ICTOT = ICTOT + 1
            IVENU(INUM) = ICTOT
         END IF
      END DO
      IF ( IPRRSP .GT. 20 ) THEN
         WRITE(LUPRI,' (/A)')
     &   ' INDEX VECTOR FOR E(2) TRANSFORMED TRIAL VECTORS'
         WRITE(LUPRI,*) ( IVENU(INUM),INUM=1,N )
         WRITE(LUPRI,' (/A)')
     &   ' IBTYP VECTOR FOR E(2) TRANSFORMED TRIAL VECTORS'
         WRITE(LUPRI,*)
     &   ( IBTYP(INUM),INUM=KOFFTY+IZRED+1,KOFFTY+IZRED+N )
      END IF
      IF (LINEQ) THEN
C
C EXTEND REDUCED GRADIENTS
C
         IZYRED= 2*IZRED
         IF (IPRRSP.GT.100) THEN
            WRITE(LUPRI,*)'RSPRED: GRADIENT VECTOR (Z-part and Y-part)'
            CALL OUTPUT(GP,1,KZVAR,1,2,KZVAR,2,1,LUPRI)
            WRITE(LUPRI,*)KEXSIM,' REDUCED GRADIENT VECTORS DIM IZYRED',
     &         IZYRED
         ENDIF
         IF (IZYRED.GT.0) THEN
            DO I=1,KEXSIM
               IF (IPRRSP.GT.100) CALL OUTPUT(REDGP((I-1)*IZYRED+1),
     &            1,2,1,IZRED,2,IZRED,-1,LUPRI)
               CALL DCOPY(IZYRED,REDGP((I-1)*IZYRED+1),1,EIVEC(1,I),1)
            END DO
         END IF
         DO 12 J = 1,KEXSIM
            ISTBC = 1
            ISTBO = 1
            DO 10 I=1,N
               IF (IBTYP(IZRED+KOFFTY+I).EQ.JBCNDX) THEN
                  PZY = DDOT(KZCONF,BCVEC(ISTBC),1,GP,1)
                  PYZ = DDOT(KZCONF,BCVEC(ISTBC),1,GP(1+KZVAR),1)
                  ISTBC = ISTBC + KZCONF
               ELSE
                  PZY = DDOT(KZWOPT,BOVEC(ISTBO),1,GP(1+KZCONF),1)
     &                + DDOT(KZWOPT,BOVEC(ISTBO+KZWOPT),1,
     &                       GP(1+KZVAR+KZCONF),1)
                  PYZ = DDOT(KZWOPT,BOVEC(ISTBO),1,
     &                       GP(1+KZVAR+KZCONF),1)
     &             + DDOT(KZWOPT,BOVEC(ISTBO+KZWOPT),1,GP(1+KZCONF),1)
                  ISTBO = ISTBO + KZYWOP
               ENDIF
               EIVEC(IZYRED+2*I-1,J) = PZY
               EIVEC(IZYRED+2*I,J)   = PYZ
 10         CONTINUE
 12      CONTINUE
         NTOT = KEXSIM*KZYRED
         CALL DCOPY(NTOT,EIVEC,1,REDGP,1)
         IF (IPRRSP.GE.15) THEN
            DO I = 1,KEXSIM
               WRITE(LUPRI,'(/I3,A,I5)') I,
     &            '. REDUCED GRADIENT VECTOR (GZ, GY) DIM KZRED', KZRED
               CALL OUTPUT(REDGP(1+(I-1)*KZYRED),1,2,1,KZRED,
     &            2,KZRED,-1,LUPRI)
            END DO
         ENDIF
      ENDIF
C
C  EXTEND REDUCED S- AND E-MATRICES
C
C  b-vectors times new S-and E-vectors
C
      REWIND (LURSP3)
      IF (KOFFTY.EQ.1) READ (LURSP3)
      KBOVEC = 1
      KBCVEC = 1
      DO 15 J = 1,KZRED
         IBTYPJ = IBTYP(J+KOFFTY)
         IF (IBTYPJ.EQ.JBONDX) THEN
            JZYVAR = KZYWOP
            JZVAR = KZWOPT
            JZOFF = 1+KZCONF
            JYOFF = KZVAR+JZOFF
         ELSE
            JZYVAR = KZCONF
            JZVAR = KZCONF
            JZOFF = 1
            JYOFF = KZVAR+1
         ENDIF
         IF(J.LE.IZRED) THEN
            CALL READT(LURSP3,JZYVAR,WRK)
            KMIN = 1
         ELSE
            IF (IBTYP(J+KOFFTY).EQ.JBONDX) THEN
               CALL DCOPY(JZYVAR,BOVEC(KBOVEC),1,WRK,1)
               KBOVEC = KBOVEC + JZYVAR
            ELSE
               CALL DCOPY(JZYVAR,BCVEC(KBCVEC),1,WRK,1)
               KBCVEC = KBCVEC + JZYVAR
            ENDIF
            KMIN = J-IZRED
         END IF
         DO 13 K=KMIN,N
            KEFF = IVENU(K)
            X1 = DDOT(JZVAR,WRK,         1,EVECS(JZOFF,KEFF),1)
            IF ( IBTYPJ .EQ. JBONDX ) X1 = X1
     *         + DDOT(JZVAR,WRK(1+JZVAR),1,EVECS(JYOFF,KEFF),1)
            X2 = DDOT(JZVAR,WRK,         1,EVECS(JYOFF,KEFF),1)
            IF ( IBTYPJ .EQ. JBONDX ) X2 = X2
     *         + DDOT(JZVAR,WRK(1+JZVAR),1,EVECS(JZOFF,KEFF),1)
            REDE(2*J  +IROW(2*(IZRED+K)))    = X1
            REDE(2*J-1+IROW(2*(IZRED+K)-1) ) = X1
            IF(J.NE.(K+IZRED)) REDE(2*J+IROW(2*(IZRED+K)-1)) = X2
            REDE(2*J-1+IROW(2*(IZRED+K)))    = X2
   13    CONTINUE
   15 CONTINUE
C
C *********************************************************
C
      NCSIM = 0
      NOSIM = 0
      DO 200 I = 1,N
         IF ( IBTYP(KOFFTY+IZRED+I) .EQ. JBONDX ) THEN
            NOSIM = NOSIM + 1
         ELSE
            NCSIM = NCSIM + 1
         END IF
 200  CONTINUE
      NTOT = (NOSIM+NCSIM) * KZYVAR
      CALL DZERO(EVECS,NTOT)
      CALL RSPSLI(NCSIM,NOSIM,BCVEC,BOVEC,
     *            UDV,EVECS,XINDX,WRK,LWRK)
C     CALL RSPSLI(NCSIM,NOSIM,ZYCVEC,ZYOVEC,
C    *                  UDV,SVEC,XINDX,WRK,LWRK)
      REWIND (LURSP3)
      IF (KOFFTY.EQ.1) READ (LURSP3)
      KBOVEC = 1
      KBCVEC = 1
      DO 25 J = 1,KZRED
         IBTYPJ = IBTYP(J+KOFFTY)
         IF (IBTYPJ.EQ.JBONDX) THEN
            JZYVAR = KZYWOP
            JZVAR = KZWOPT
            JZOFF = 1+KZCONF
            JYOFF = KZVAR+JZOFF
         ELSE
            JZYVAR = KZCONF
            JZVAR = KZCONF
            JZOFF = 1
            JYOFF = KZVAR+1
         ENDIF
         IF(J.LE.IZRED) THEN
            CALL READT(LURSP3,JZYVAR,WRK)
            KMIN = 1
         ELSE
            IF (IBTYP(J+KOFFTY).EQ.JBONDX) THEN
               CALL DCOPY(JZYVAR,BOVEC(KBOVEC),1,WRK,1)
               KBOVEC = KBOVEC + JZYVAR
            ELSE
               CALL DCOPY(JZYVAR,BCVEC(KBCVEC),1,WRK,1)
               KBCVEC = KBCVEC + JZYVAR
            ENDIF
            KMIN = J-IZRED
         END IF
         DO 23 K=KMIN,N
            KEFF = IVENU(K)
            X1 = DDOT(JZVAR,WRK,         1,EVECS(JZOFF,KEFF),1)
            IF (IBTYPJ .EQ. JBONDX) X1 = X1
     *         + DDOT(JZVAR,WRK(1+JZVAR),1,EVECS(JYOFF,KEFF),1)
            X2 = DDOT(JZVAR,WRK,         1,EVECS(JYOFF,KEFF),1)
            IF (IBTYPJ .EQ. JBONDX) X2 = X2
     *         + DDOT(JZVAR,WRK(1+JZVAR),1,EVECS(JZOFF,KEFF),1)
            REDS(2*J-1+IROW(2*(IZRED+K)-1)) =  X1
            REDS(2*J  +IROW(2*(IZRED+K)))   = -X1
            IF(J.NE.(IZRED+K)) REDS(2*J+IROW(2*(IZRED+K)-1)) = X2
            REDS(2*J-1+IROW(2*(IZRED+K)))   = -X2
   23    CONTINUE
   25 CONTINUE

      deallocate( IVENU )
C
C *********************************************************
C Section 2: find eigenvalues and -vectors of reduced L-matrix
C
 1000 CONTINUE
      IF (IPRRSP .GE. 10) THEN
         WRITE (LUPRI,'(/A)') ' RSPRED: REDUCED HESSIAN MATRIX:'
         CALL OUTPAK(REDE,KZYRED,-1,LUPRI)
         WRITE (LUPRI,'(/A)') ' RSPRED: REDUCED METRIC MATRIX:'
         CALL OUTPAK(REDS,KZYRED,-1,LUPRI)
      END IF
      IF (ICTL .EQ. 1) GO TO 2000
C
      IF (LINEQ) THEN
C
C     Solve reduced linear response problem in subspace
C
         DO 65 IGP = 1,KEXSIM
            AFREQ = EIVAL(IGP)
            IJ=0
            DO 70 I=1,KZYRED
               DO 75 J=1,I
                  IJ = IJ + 1
                  BCVEC(IJ) = REDE(IJ) - AFREQ*REDS(IJ)
   75          CONTINUE
   70       CONTINUE
            CALL DCOPY(KZYRED,REDGP((IGP-1)*KZYRED+1),1,
     *                  EIVEC(1,IGP),1)
         IF (IPRRSP.GE.9) THEN
            WRITE(LUPRI,*)' AFREQ,KZYRED,IGP',AFREQ,KZYRED,IGP
         ENDIF
         IF (IPRRSP.GT.20) THEN
            WRITE(LUPRI,'(/A)') ' REDUCED E(2)-W*S(2) MATRIX:'
            CALL OUTPAK(BCVEC,KZYRED,-1,LUPRI)
            WRITE(LUPRI,*)'REDUCED GRADIENT VECTORS DIM (2,KZRED)',KZRED
            CALL OUTPUT(EIVEC(1,IGP),1,2,1,KZRED,2,KZRED,-1,LUPRI)
         ENDIF
C
C           USE DSPSLI ( WHICH CALLS LINPACK ROUTINES ) TO SOLVE
C           SYSTEM OF LINEAR EQUATIONS
C
            NSIM = 1
            KK1  = 1 + KZYRED*(KZYRED+1)/2
            CALL DSPSLI(KZYRED,NSIM,BCVEC,EIVEC(1,IGP),BCVEC(KK1),INFO,
     *                  DET,ISNDX)
C
C           CALL DSPSLI(N,NSIM,AP,B,KPIVOT,INFO,DET,INERT)
C
            IF (INFO .NE. 0) THEN
               WRITE (LUPRI,'(/A,I5)')
     *              ' *** ERROR IN RSPRED.DSPSLI, INFO =',INFO
               CALL QUIT('*** ERROR IN RSPRED.DSPSLI')
            END IF
C
C
            INEG     = ISNDX(2)
            ISNDX(2) = ISNDX(3)
            ISNDX(3) = INEG
            IF (INEG.GT.0 .AND. .NOT.ABSORP)
     *           WRITE(LUPRI,'(/2A,I3,A,L8,A,I6)')
     *         ' *** INFO, negative eigenvalues in reduced matrix.',
     &         ' Symmetry',KSYMOP,',  triplet',TRPLET,
     &         ' Dimension of reduced matrix =',KZYRED
            IF (IPRRSP.GE.4 .OR. (.NOT.ABSORP .AND. INEG.GT.0)) THEN
               POL=DDOT(KZYRED,EIVEC(1,IGP),1,REDGP((IGP-1)*KZYRED+1),1)
               WRITE(LUPRI,1800) AFREQ,KZRED,POL,ISNDX,DET
            END IF
 1800       FORMAT(/' GP * SOLUTION vector at frequency',F12.5,' a.u.',
     *             /' after',I5,' linear transformations is',1P,D20.8,
     *             /' INERTIA (POS,ZER,NEG) of reduced matrix is',3I5,
     *             /' Determinant of reduced matrix is',0P,F6.2,
     *              ' * 10 **',F5.1)
C
            IF (IPRRSP .GE. 9) THEN
               WRITE (LUPRI,'(/A)')
     *            ' REDUCED PROPERTY VECTOR AND REDUCED SOLUTION:'
               WRITE (LUPRI,'(1P,/I5,2D15.6/I5,2D15.6)')
     *           (K,REDGP((IGP-1)*KZYRED+K), EIVEC(K,IGP), K=1,KZYRED)
            END IF
 65      CONTINUE
      ELSE
C
C
C Solve reduced general problem in subspace
C
C Expand reduced matrices to square matrices
C
         N2ZYRED = KZYRED*KZYRED
         CALL DSPTGE(KZYRED,REDE,BCVEC)
         CALL DSPTGE(KZYRED,REDS,BCVEC(1+N2ZYRED))
C
         KK1 = 1   + 2*N2ZYRED
         KK2 = KK1 + KZYRED
         KK3 = KK2 + KZYRED
         KK4 = KK3 + N2ZYRED
C
C        Use EISPACK routine for real general matrices in
C        generalized eigenvalue problem
C
         CALL RGG(KZYRED,KZYRED,BCVEC,BCVEC(1+N2ZYRED),EIVAL,
     *            BCVEC(KK1),BCVEC(KK2),1,EIVEC,IERR)
C        CALL RGG(NM,N,A,B,ALFR,ALFI,BETA,MATZ,Z,IERR)
C
         IF ( IERR.NE.0 ) THEN
            WRITE(LUPRI,'(/A)') ' FATAL ERROR **RESPONS:'//
     &      ' Eigenvalue routine RGG failed'
            CALL QUIT(' RSPRED: Eigenvalue routine failed')
         END IF
C
         IF (IPRRSP .GE. 25) THEN
           WRITE (LUPRI,'(/A)') ' REDUCED EIGENVECTORS BEFORE ORDERING:'
           CALL OUTPUT(EIVEC,1,KZYRED,1,KZYRED,KZYRED,KZYRED,-1,LUPRI)
           WRITE (LUPRI,'(/A)') ' REDUCED EIGENVALUES (2,KZRED):'
           CALL OUTPUT(EIVAL,1,2,1,KZRED,2,KZRED,-1,LUPRI)
         END IF
C
         CALL DSPTGE(KZYRED,REDS,BCVEC)
         CALL RSPORD(BCVEC,EIVEC,EIVAL,BCVEC(KK1),BCVEC(KK2),
     &               BCVEC(KK3),BCVEC(1+N2ZYRED),ISNDX)
C        CALL RSPORD(SRED,EIVEC,ALFAR,ALFI,BETA,BVECS1,BVECS2,ISNDX)
C
         IF (ISNDX(2) .GT. 0 .OR. ISNDX(1) .NE. ISNDX(3) .OR.
     &       IPRRSP .GE. 30) THEN
            NNZYRED = (KZYRED*(KZYRED+1))/2
            CALL DCOPY(NNZYRED,REDS,1,BCVEC,1)
            WRITE (lupri,*) 'S[2] reduced matrix'
            CALL outpak(bcvec,kzyred,-1,lupri)
            CALL JACO(BCVEC,BCVEC(KK1),KZYRED,KZYRED,0,
     &          BCVEC(KK2),BCVEC(KK3))
C           CALL JACO(F,V,NB,NMAX,NROWV,BIG,JBIG)
            DO, I = 1,KZYRED
               BCVEC(I) = BCVEC( (I*I+I)/2 )
            END DO
            ! CALL ORDER(WRK,BCVEC,KZYRED,0)
            WRITE(LUPRI,*) 'The eigenvalues of S[2] reduced matrix'
            CALL OUTPUT(BCVEC,1,2,1,KZRED,2,KZRED,-1,LUPRI)

            CALL DCOPY(NNZYRED,REDE,1,BCVEC,1)
            WRITE (lupri,*) 'E[2] reduced matrix'
            CALL outpak(bcvec,kzyred,-1,lupri)
            CALL JACO(BCVEC,BCVEC(KK1),KZYRED,KZYRED,0,
     &          BCVEC(KK2),BCVEC(KK3))
C           CALL JACO(F,V,NB,NMAX,NROWV,BIG,JBIG)
            DO, I = 1,KZYRED
               BCVEC(I) = BCVEC( (I*I+I)/2 )
            END DO
            ! CALL ORDER(WRK,BCVEC,KZYRED,0)
            WRITE(LUPRI,*) 'The eigenvalues of E[2] reduced matrix'
            CALL OUTPUT(BCVEC,1,2,1,KZRED,2,KZRED,-1,LUPRI)

         END IF
         IF (IPRRSP .GE. 10) THEN
            IF (IPRRSP .GE. 25) THEN
               WRITE (LUPRI,'(/A)')
     &               ' REDUCED EIGENVECTORS AFTER ORDERING:'
               CALL OUTPUT(EIVEC,1,KZYRED,1,KZYRED,KZYRED,KZYRED,
     &                     -1,LUPRI)
            END IF
            WRITE(LUPRI,'(//2A,3I5//A/A)')
     *         ' NUMBER OF EIGENVALUES WITH POSITIVE,',
     *         ' ZERO, AND NEGATIVE METRIC:',ISNDX,
     *         ' THE EIGENVALUES WITH POSITIVE METRIC:',
     *         '      NUMBER      EIGENVALUE'
            IPOS = ISNDX(1)
            DO, I=1,IPOS
               IF (ABS(EIVAL(I)) .LT. DRTEST) THEN
                 WRITE(LUPRI,'(I10,1P,D20.8,A/)') I,EIVAL(I),
     *            ' *** RSPRED  WARNING **** ZERO EIGENVALUE.'
               ELSE
                  WRITE(LUPRI,'(I10,1P,D20.8)') I,EIVAL(I)
               END IF
            END DO
         END IF
      ENDIF
C
C *** End of subroutine RSPRED
C
 2000 CONTINUE
      CALL QEXIT('RSPRED')
      RETURN
      END
C  /* Deck ppst */
      SUBROUTINE PPST(IBTYP,XINDX,FCAC,H2AC,WRK,NSIM,BINMEM,LWRK)
C
C PURPOSE: CREATE START VECTORS FOR SOLUTION OF GENERALIZED
C          RSP EIGENVALUE PROBLEM
C
#include "implicit.h"
#include "dummy.h"
C
      DIMENSION IBTYP(*),WRK(*)
      LOGICAL BINMEM
      DIMENSION XINDX(*),FCAC(*),H2AC(*)
C
      PARAMETER (MAXST=400)
      DIMENSION STVAL(MAXST),ISTART(MAXST),IEQUAL(MAXST)
C
#include "ibndxdef.h"
C
      PARAMETER ( TOLEQL=1.0D-6 )
      PARAMETER ( D0=0.0D0 , D1=1.0D0 , DM1 = -1.0D0 )
C
C wrkrsp.h : DETERM, ?
C infinp.h : ISPIN
C
#include "priunit.h"
#include "maxorb.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "infinp.h"
#include "inftap.h"
#include "infdim.h"
#include "infpri.h"
#include "infopt.h"
#include "phpinf.h"
#include "thrldp.h"
      LOGICAL DO_PHPVEC
C
C WORK SPACE ALLOCATION
C
      CALL QENTER('PPST')
C
      CALL PHPINI(LPHPMX,KZCONF,KZWOPT,MAXPHP)
      KDIAE  = 1
      KDIASO = KDIAE  + LPHPMX
      KTOT   = KDIASO + KZWOPT
      LTOT   = LWRK   - KTOT
      IF (LTOT.LT.0) CALL ERRWRK('PPST 1',KTOT-1,LWRK)
      IF (KZCONF.GT.0) THEN
         CALL PHPDSK(WRK(KDIAE),LPHPMX,.FALSE.)
      ELSE
         CALL RSPEDG(WRK(KDIAE))
      END IF
      IF (IPRRSP.GT.10) THEN
         WRITE(LUPRI,*)' DIMENSION SPECIFIED OF PHP, MAXPHP=',MAXPHP
         WRITE(LUPRI,*)' DIMENSION USED , NPCVAR=',NPCVAR
      END IF
      CALL RSPSOD(WRK(KDIASO))
      NCSIM = 0
      KEXSTV = MIN(KEXSTV+KOFFTY,KZVAR)
      IF (KEXSTV.GT.MAXST ) THEN
         WRITE (LUPRI,90) KEXSTV,MAXST
 90      FORMAT(/' ** PPST ** NUMBER OF START VECTORS',I5,
     *          /'   EXCEED THE MAXIMUM ALLOWED NUMBER',I5)
         CALL QUIT('PPST: TOO MANY START VECTORS')
      END IF
      KEXST2 = MIN(5*KEXSTV+8,KZVAR)
      KEXST2 = MIN(KEXST2,  MAXST)
      IF (IPRRSP .GT.110 ) THEN
         WRITE(LUPRI,'(/A)') ' DIAGONAL of E(2)'
         CALL OUTPUT(WRK(KDIAE),1,KZVAR,1,1,KZVAR,1,1,LUPRI)
         WRITE(LUPRI,'(/A)') ' DIAGONAL of orbital part of S(2)'
         CALL OUTPUT(WRK(KDIASO),1,KZWOPT,1,1,KZWOPT,1,1,LUPRI)
      END IF
      DO 53 I = 1,KZWOPT
         IF (ABS(WRK(I+KDIASO-1)) .GT. TOLEQL) THEN
            WRK(I+KDIAE-1+KZCONF) = WRK(I+KDIAE-1+KZCONF) /
     *                              ABS(WRK(I+KDIASO-1))
         END IF
 53   CONTINUE
      IF (IPRRSP .GT.150 ) THEN
         WRITE(LUPRI,'(/A)') ' DIAGONAL E(2)/ABS(S(2))  '
         CALL OUTPUT(WRK(KDIAE),1,KZVAR,1,1,KZVAR,1,1,LUPRI)
      END IF
      DO 100 K=1,KEXST2
         ISTART(K) = K
         STVAL(K)  = WRK(K)
 100  CONTINUE
      KK=KEXST2
C
C     Find KEXST2 lowest diagonal elements
C     REPEAT UNTIL ...
 105  CONTINUE
         DO 110 I=1,KEXST2-1
            DO 120 J=I+1,KEXST2
               IF ((STVAL(J)-STVAL(I)) .LT. D0) THEN
                  X =STVAL(I)
                  II=ISTART(I)
                  STVAL(I) =STVAL(J)
                  ISTART(I)=ISTART(J)
                  STVAL(J) =X
                  ISTART(J)=II
               ENDIF
 120        CONTINUE
 110     CONTINUE
         STMAX=STVAL(KEXST2)
C
 125     CONTINUE
            KK=KK+1
            IF (KK.LE.KZVAR) THEN
               IF ((STMAX-WRK(KK)).GT.D0) THEN
                  ISTART(KEXST2) = KK
                  STVAL(KEXST2)  = WRK(KK)
                  GO TO 105
C     ^--------------------
               ENDIF
               GO TO 125
C        ^--------------
            END IF

C 16. Dec. 2014-hjaaj: No plus combinations for CSFs nor triplet excitations
C     nor non-singlet reference state
C     (for determinants the plus combinations of degenerate elements are used
C     to avoid triplet CI states in calculations of singlet excitations,
C     note that determinants are always used when KSYMOP .ne. 1)
C
      IF ((.NOT.DETERM .AND. KSYMOP.EQ.1)  ! no plus combinations for CSFs in KSYMOP.eq.1, ground state symmetry
     &      .OR. TRPLET                    ! never plus combinations for triplet excitations!
     &      .OR. ISPIN.NE.1 ) THEN         ! not a singlet reference state

         DO JEXST2 = 1, KEXST2
            IEQUAL(JEXST2) = 1 ! no plus combinations!
         END DO

C        Do include ALL (near)-degenerate diagonal elements as start vectors
C        (Otherwise we will probably break spatial and/or spin symmetry) /hjaaj Feb. 2019
         JEXSTV = KEXSTV
         DO IEXST2 = KEXSTV+1, KEXST2
            IF ( (STVAL(IEXST2)-STVAL(IEXST2-1)).GT.TOLEQL ) EXIT
            JEXSTV = IEXST2
         END DO
         IF (KEXSTV .NE. JEXSTV) THEN
            WRITE(LUPRI,'(/A,I0,A,I0,A,I0,A/A)')
     &     ' INFO Symmetry ',KSYMOP,
     &      ': Number of roots increased from ',KEXSTV,' to ',
     &      JEXSTV,' because of (near)-degenerate diagonal elements.',
     &     ' (Otherwise the calculation would probably break symmetry.)'
            KEXSTV = JEXSTV
         END IF

         GO TO 2000 ! skip next block where "plus combinations" are made

      END IF

C --- block for constructing plus combinations of (near)-degenerate
C     diagonal elements among determinants when singlet reference and singlet excitations

      JEXSTV = 0
      JEXST2 = 0

C 19. Feb. 1988: Find number of (near)-degenerate diagonal elements
C 19. Jun. 1995-hjaaj: not for orbital trial vectors
  160 CONTINUE
         IEXST2 = JEXST2
  165    IEXST2 = IEXST2 + 1
         IF (IEXST2 .LT. KEXST2 .AND.
     &      ISTART(IEXST2) .LE. KZCONF .AND. ! no plus combinations of orbital trial vectors
     &      (STVAL(IEXST2+1)-STVAL(IEXST2)).LE.TOLEQL) GO TO 165
         JEXSTV = JEXSTV + 1
         IEQUAL(JEXSTV) = IEXST2 - JEXST2
         JEXST2 = IEXST2
      IF (JEXSTV .LT. KEXSTV .AND. JEXST2 .LT. KEXST2) GO TO 160
C
      IF (JEXSTV .LT. KEXSTV .OR. JEXST2 .EQ. KEXST2) THEN
         IF (KEXST2 .LT. KZVAR) THEN
            WRITE(LUPRI,'(//A/)')
     &         ' ERROR in PPST: too few diagonal elements sorted.'
            IPRRSP = 7
         END IF
      END IF
      IF (IPRRSP.GE.7) THEN
         WRITE(LUPRI,*) 'KEXST2,JEXST2 :',KEXST2,JEXST2
         WRITE(LUPRI,*) '(I,(ISTART(I),STVAL(I)),I=1,KEXST2)'
         DO I = 1,KEXST2
            WRITE(LUPRI,*) I,ISTART(I),STVAL(I)
         END DO
         WRITE(LUPRI,*) 'KEXSTV,JEXSTV :',KEXSTV,JEXSTV
         WRITE(LUPRI,*) '(I, IEQUAL(I), I = 1,JEXSTV)'
         DO I = 1,JEXSTV
            WRITE(LUPRI,*) I,IEQUAL(I)
         END DO
      END IF
      IF (JEXSTV .LT. KEXSTV .OR. JEXST2 .EQ. KEXST2) THEN
         IF (KEXST2 .LT. KZVAR) THEN
            CALL QUIT('ERROR in PPST, too few diagonal elements sorted')
         ELSE IF ( JEXSTV.LT.KEXSTV ) THEN
            WRITE(LUPRI,*) 'KEXSTV reduced from',KEXSTV,' to',JEXSTV
            WRITE(LUPRI,*) 'because of degenerate diagonal elements.'
            KEXSTV = JEXSTV
            KEXCNV = JEXSTV
            KEXSIM = JEXSTV
         END IF
      END IF
C
C     Can we use PHP vectors for CI start vectors ?
C
 2000 CONTINUE
      IF ( KEXSTV .GT. NPCVAR .OR. KEXSTV .GE. KZCONF) THEN
         DO_PHPVEC = .FALSE.
      ELSE
         DO_PHPVEC = .TRUE.
      END IF

      NOVLMX = 0 ! to become index of CI start vector with max overlap with reference vector
      REWIND (LURSP3)
      REWIND (LURSP5)
      IF (KOFFTY.EQ.1) THEN
C
C PUT REFERENCE CONFIGURATION AS FIRST RECORD ON LURSP3
C Put zero vector as first record on LURSP5
C
         KCREF  = KDIASO
         KCRFRE = KCREF  + KZCONF
Chj      KTOT   = MAX(KCRFRE + MAXPHP, KCREF + KZYVAR)
         KTOT   = KCRFRE + MAXPHP
         LTOT   = LWRK   - KTOT
         IF (LTOT.LT.0) CALL ERRWRK('PPST 2',KTOT-1,LWRK)
         CALL GETREF(WRK(KCREF),KZCONF)
         IF (IPRRSP.GT.115) THEN
            WRITE(LUPRI,'(/A,2I12)')
     *         ' REFERENCE CONFIGURATION. KZCONF,NCREF =',KZCONF,NCREF
            CALL OUTPUT(WRK(KCREF),1,KZCONF,1,1,KZCONF,1,1,LUPRI)
         END IF
         CALL WRITT(LURSP3,KZCONF,WRK(KCREF))
         NCNUM = 0
         IRX   = 0
         ISIMX = 1
         OVLPMX = D0
         IF (DO_PHPVEC) THEN
            CALL GATVEC(WRK(KCRFRE),WRK(KCREF), WRK(KPNTR),NPCVAR)
Chj         ... PHP pointer in KPNTR is also based on KZCONF
            IF (IPRRSP.GT.80) THEN
               WRITE(LUPRI,*)' GATHERED REFERENCE VECTOR; length',NPCVAR
               CALL OUTPUT(WRK(KCRFRE),1,NPCVAR,1,1,NPCVAR,1,1,LUPRI)
            END IF
         END IF
         DO  80 ISIM = 1,KEXSTV
            IF (ISTART(ISIMX) .LE. KZCONF) THEN
               IF (DO_PHPVEC) THEN
                  NCNUM  = NCNUM + 1
                  KRDVEC = KPEVEC + (NCNUM-1)*NPCVAR
                  OVLPRF = DDOT(NPCVAR,WRK(KCRFRE),1,WRK(KRDVEC),1)
               ELSE
                  NCNUM = NCNUM + 1
                  OVLPRF = D0
                  DO 82 IEQL = 1,IEQUAL(ISIM)
                     INUM = ISTART(ISIMX+IEQL-1)
                     OVLPRF = OVLPRF +  WRK(KCREF-1+INUM)
   82             CONTINUE
                  FEQUAL = IEQUAL(ISIM)
                  OVLPRF = OVLPRF/SQRT(FEQUAL)
               END IF
               IF (ABS(OVLPRF).GT.OVLPMX) THEN
                  OVLPMX = ABS(OVLPRF)
                  NOVLMX = NCNUM
               END IF
            ENDIF
            ISIMX = ISIMX + IEQUAL(ISIM)
 80      CONTINUE
         IF (IPRRSP.GT.6) THEN
            WRITE(LUPRI,*)' REFERENCE VECTOR HAS LARGEST OVERLAP WITH'
            IF (DO_PHPVEC) THEN
               WRITE(LUPRI,*)'EIGENVECTOR NO',NOVLMX,' OF PHP'
            ELSE
               WRITE(LUPRI,*)'SORTED VECTOR NO',NOVLMX
            END IF
            WRITE(LUPRI,*)' OVERLAP =',OVLPMX
         END IF
         WRITE (LURSP5) D0,D0,D0,D0
Chj      ... record 1 on LURSP5 should never be read when KOFFTY.eq.1 !
         IBTYP(1) = JBCNDX
      ENDIF
      NLOAD = 0
      NTSIM = KOFFTY
      ISIMX = 1
      NCNUM = 0
      NCCOUN = 0

      NCTOT = 0
      NOTOT = 0

      DO 150 ISIM = 1,KEXSTV,NSIM
         NBX    = MIN(NSIM,(KEXSTV+1-ISIM))
         NLOAD  = NLOAD + 1
         NCSIM  = 0
         NOSIM  = 0
         IRX    = 0
         IRJUMP = 0
         DO 200 IR = 1,NBX
            IF (ISTART(ISIMX+IRX) .LE. KZCONF) THEN
               NCNUM = NCNUM + 1
               IF(  DO_PHPVEC .AND. (NCNUM.EQ.NOVLMX) ) THEN
                  IRJUMP = IR
               ELSE
                  NCSIM = NCSIM + 1
                  IBTYP(NTSIM+IR) = JBCNDX
               END IF
            ELSE
               NOSIM = NOSIM + 1
               IBTYP(NTSIM+IR) = JBONDX
            ENDIF
            IRX = IRX + IEQUAL(ISIM-1+IR)
 200     CONTINUE
         IF (IPRRSP.GT.40) THEN
            IF (IRJUMP.NE.0) THEN
               WRITE(LUPRI,*)' IRJUMP =',IRJUMP
            END IF
         END IF
         KCOFF = KDIASO
         KOOFF = KCOFF + NCSIM*KZCONF
         KWRK  = KOOFF + NOSIM*KZYWOP
         NTOT  = NCSIM*KZCONF + NOSIM*KZYWOP
         CALL DZERO(WRK(KCOFF),NTOT)
         IRX   = 0
         DO 250 IR = 1,NBX
            IF ( IR .EQ. IRJUMP ) THEN
               NCCOUN = NCCOUN + 1
            ELSE
               IF (IBTYP(NTSIM+IR) .EQ. JBCNDX) THEN
                  IF (DO_PHPVEC) THEN
                     NCCOUN = NCCOUN + 1
                     KEVSC = KPEVEC + (NCCOUN-1)*NPCVAR
                     IF (IPRRSP.GT.85) THEN
                        WRITE(LUPRI,*)' START VECTOR BEFORE SCATTERING'
                        CALL OUTPUT(WRK(KEVSC),1,NPCVAR,1,1,
     *                              NPCVAR,1,1,LUPRI)
                     END IF
                     CALL SCAVEC(WRK(KCOFF),WRK(KEVSC),
     *                           WRK(KPNTR),NPCVAR)
                     IF (IPRRSP.GT.109) THEN
                        WRITE(LUPRI,*)' START VECTOR AFTER SCATTERING'
                        CALL OUTPUT(WRK(KCOFF),1,KZCONF,1,1,
     *                              KZCONF,1,1,LUPRI)
                     END IF
                  ELSE
                     DO 242 IEQL = 1,IEQUAL(ISIM-1+IR)
                       INUM = ISTART(ISIMX+IRX+IEQL-1)
                       WRK(KCOFF-1+INUM) = D1
                       IF ( TRPLET .AND. (IEQUAL(ISIM-1+IR).EQ.2) .AND.
     *                   (IEQL.EQ.2)) WRK(KCOFF-1+INUM) = -D1
  242                CONTINUE
                  END IF
                  KCOFF = KCOFF + KZCONF
               ELSE
                  DO 244 IEQL = 1,IEQUAL(ISIM-1+IR)
                     INUM = ISTART(ISIMX+IRX+IEQL-1) - KZCONF
                     WRK(KOOFF-1+INUM) = D1
  244             CONTINUE
                  KOOFF = KOOFF + KZYWOP
               END IF
            END IF
            IRX = IRX + IEQUAL(ISIM-1+IR)
 250     CONTINUE
         NLSIM = NOSIM + NCSIM
         DO 270 IR = 1,NLSIM
            IF (IR.LE.NCSIM) THEN
               IBTYP(NTSIM+IR) = JBCNDX
            ELSE
               IBTYP(NTSIM+IR) = JBONDX
            ENDIF
 270     CONTINUE
         IF ((IPRRSP.GT.110).AND.(NCSIM.GT.0)) THEN
            WRITE(LUPRI,*)' NCSIM CSF TRIAL VECTORS',NCSIM
            CALL OUTPUT(WRK(KDIASO),1,KZCONF,1,NCSIM,
     *                  KZCONF,NCSIM,1,LUPRI)
         ENDIF
         IF ((IPRRSP.GT.110).AND.(NOSIM.GT.0)) THEN
            WRITE(LUPRI,*)' NOSIM ORBITAL TRIAL VECTORS',NOSIM
            CALL OUTPUT(WRK(KDIASO+NCSIM*KZCONF),1,KZYWOP,
     *             1,NOSIM,KZYWOP,NOSIM,1,LUPRI)
         ENDIF
         CALL RSPORT(WRK(KDIASO),NLSIM,NTSIM,IBTYP,THRLDP,
     *               WRK(KWRK),LURSP3)
C        CALL RSPORT(BVECS,NBX,NBPREV,IBTYP,THRLDP,OLDVEC,LU3)
         IF (IPRRSP.GT.110) THEN
            NCSIM = 0
            NOSIM = 0
            DO 3010 IX = 1,NLSIM
               IF (IBTYP(NTSIM+IX).EQ.JBCNDX) THEN
                  NCSIM = NCSIM + 1
               ELSE
                  NOSIM = NOSIM + 1
               ENDIF
 3010       CONTINUE
            IF (NCSIM.GT.0) THEN
               WRITE(LUPRI,*)' NCSIM CSF TRIAL VECTORS',NCSIM
               CALL OUTPUT(WRK(KDIASO),1,KZCONF,1,NCSIM,
     *                     KZCONF,NCSIM,-1,LUPRI)
            ENDIF
            IF (NOSIM.GT.0) THEN
               WRITE(LUPRI,*)' NOSIM ORBITAL TRIAL VECTORS',NOSIM
               CALL OUTPUT(WRK(KDIASO+NCSIM*KZCONF),1,KZYWOP,
     *                   1,NOSIM,KZYWOP,NOSIM,-1,LUPRI)
            END IF
         END IF
         NTSIM = NTSIM + NLSIM
         ISIMX = ISIMX + IRX
         NCTOT = NCTOT + NCSIM
         NOTOT = NOTOT + NOSIM
 150  CONTINUE
      KEXSTV = NTSIM - KOFFTY
      IF (NLOAD.EQ.1) THEN
         BINMEM = .TRUE.
         IF (KDIASO .NE. 1) THEN
            NTOT = NCSIM*KZCONF + NOSIM*KZYWOP
            DO 300 I = 1,NTOT
               WRK(I) = WRK((KDIASO-1)+I)
  300       CONTINUE
C           CALL DCOPY(NTOT,WRK(KDIASO),1,WRK,1)
         END IF
      ELSE
         BINMEM = .FALSE.
      ENDIF
      IF (IPRRSP.GT.2) THEN
         IF (DO_PHPVEC) THEN
            WRITE(LUPRI,*)NCTOT,' START configuration VECTORS USING'//
     &         ' LOWEST EIGENVALUES OF PHP'
         ELSE
            WRITE(LUPRI,*)NCTOT,' START configuration VECTORS USING'//
     &          ' LOWEST DIAGONAL HESSIAN ELEMENTS'
         END IF
         WRITE(LUPRI,*)NOTOT,' START orbital VECTORS'
      END IF
C
C     END OF PPST
C
      CALL QEXIT('PPST')
      RETURN
      END
C  /* Deck rsptr1 */
      SUBROUTINE RSPTR1(NSIM,UDV,ZYMAT,DENA,DENB)
C
C PURPOSE:
C   CALCULATE MODIFIED DENSITY MATRICES
C
C   DENB(R,Y) = SUM(X)  UDV(X,Y) * ZYMAT(X,R)
C
C   DENA(R,Y) = SUM(X)  UDV(Y,X) * ZYMAT(R,X)
C
C   THE MATRIX DENB IS STORED WITHOUT USING SYMMETRY
C
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION UDV(NASHT,*),ZYMAT(NORBT,NORBT,*)
      DIMENSION DENA(NORBT,NASHT,*),DENB(NORBT,NASHT,*)
C
#include "maxorb.h"
#include "maxash.h"
#include "inforb.h"
#include "infind.h"
#include "infrsp.h"
#include "wrkrsp.h"
C
      CALL QENTER('RSPTR1')
C
      CALL DZERO(DENB,NORBT*NASHT*NSIM)
      CALL DZERO(DENA,NORBT*NASHT*NSIM)
      DO 100 ISIM=1,NSIM
         IASHOF =0
         DO 1000 IYSYM = 1,NSYM
            IASHIY = IASH(IYSYM)
            IORBIY = IORB(IYSYM)
            NISHIY = NISH(IYSYM)
            NASHIY = NASH(IYSYM)
            IRSYM  = MULD2H(IYSYM,KSYMOP)
            IORBIR = IORB(IRSYM)
            NORBIR = NORB(IRSYM)
            DO 1010 IY = 1,NASHIY
               IYOFF = ISX( NISHT + IASHOF + IY )
               DO 2020 IX = 1,IY
                  IXOFF = ISX( NISHT + IASHOF + IX )
                  CALL DAXPY(NORBIR,UDV(IASHOF+IY,IASHOF+IX),
     *                  ZYMAT(1+IORBIR,IXOFF,ISIM),1,
     *                  DENA(1+IORBIR,IASHOF+IY,ISIM),1)
                  IF ( IX .NE. IY ) THEN
                     CALL DAXPY(NORBIR,UDV(IASHOF+IX,IASHOF+IY),
     *                     ZYMAT(1+IORBIR,IYOFF,ISIM),1,
     *                     DENA(1+IORBIR,IASHOF+IX,ISIM),1)
                  ENDIF
 2020          CONTINUE
               DO 1020 IR = 1,NORBIR
                  DENB(IORBIR+IR,IASHOF+IY,ISIM) =
     *             DENB(IORBIR+IR,IASHOF+IY,ISIM) +
     *             DDOT(NASHIY,ZYMAT(IORBIY+NISHIY+1,IORBIR+IR,ISIM),1,
     *               UDV(IASHOF+1,IASHOF+IY),1)
 1020          CONTINUE
 1010       CONTINUE
            IASHOF = IASHOF + NASHIY
 1000    CONTINUE
 100  CONTINUE
      IF ( IPRRSP.GT.90 ) THEN
         WRITE(LUPRI,*)NSIM,' DENA MATRICES'
         DO ISIM = 1, NSIM
            CALL OUTPUT(DENA(1,1,ISIM),1,NORBT,1,NASHT,NORBT,NASHT,
     &                  -1,LUPRI)
         END DO
         WRITE(LUPRI,*)NSIM,' DENB MATRICES'
         DO ISIM = 1, NSIM
            CALL OUTPUT(DENB(1,1,ISIM),1,NORBT,1,NASHT,NORBT,NASHT,
     &                  -1,LUPRI)
         END DO
      END IF
      CALL QEXIT('RSPTR1')
      RETURN
      END
C  /* Deck rspzym */
      SUBROUTINE RSPZYM(NSIM,ZYVEC,ZYMAT)
C
C WRITTEN 9-2-1986
C
C This subroutine packs the ( NON TOTAL SYMMETRIC )
C orbital part of a trial vector into matrix form.
C The matrix is stored without use of symmetry.
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION ZYVEC(KZYWOP,*), ZYMAT(NORBT,NORBT,*)
C
#include "maxorb.h"
#include "maxash.h"
#include "infvar.h"
#include "inforb.h"
#include "infind.h"
#include "infrsp.h"
#include "wrkrsp.h"
C
      PARAMETER ( D1 = 1.0D0 )
C
      CALL QENTER('RSPZYM')
C
      CALL DZERO(ZYMAT,NORBT*NORBT*NSIM)
      DO 100 ISIM=1,NSIM
         DO 200 IG=1,KZWOPT
            I=JWOP(1,IG)
            J=JWOP(2,IG)
            ZYMAT(J,I,ISIM)=ZYVEC(IG,ISIM)
            IF (ANTTES) THEN
               ZYMAT(I,J,ISIM)=-ZYVEC(IG,ISIM)
            ELSE
               ZYMAT(I,J,ISIM)=ZYVEC(IG+KZWOPT,ISIM)
            END IF
  200    CONTINUE
  100 CONTINUE
C      CALL DZERO(ZYMAT,NORBT*NORBT*NSIM)
C      DO 300 ISIM = 1,NSIM
C         DO 310 I = 1,NORBT
C           ZYMAT(I,I,ISIM) = D1
C 310     CONTINUE
C 300  CONTINUE
      IF (IPRRSP.GT.100) THEN
         WRITE(LUPRI,*)NSIM,' ZYMAT MATRICES'
         DO 400 ISIM = 1,NSIM
            CALL OUTPUT(ZYMAT(1,1,ISIM),1,NORBT,1,NORBT,NORBT,NORBT,
     &                  -1,LUPRI)
 400     CONTINUE
      END IF
C
      CALL QEXIT('RSPZYM')
      RETURN
      END
C  /* Deck rspcve */
      SUBROUTINE RSPCVE(IOFFC,IBNUM,IBTYP,EIVAL,EIVEC,BVECS,WRK,LWRK)
C
C PURPOSE:
C   CALCULATION OF Z OR Y COMPONENT OF CONFIGUTATION PART OF AN
C   EIGENVECTOR FROM REDUCED SPACE
C
#include "implicit.h"
#include "dummy.h"
      DIMENSION BVECS(*),WRK(*)
      DIMENSION IBTYP(*),EIVAL(*),EIVEC(KZYRED,*)
C
#include "ibndxdef.h"
C
#include "priunit.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "inftap.h"
#include "infpri.h"
C
      CALL QENTER('RSPCVE')
C
      CALL DZERO(BVECS,KZCONF)
      REWIND (LURSP3)
      IF (KOFFTY.EQ.1) READ (LURSP3)
      DO 200 K = 1,KZRED
         IBTYPK = IBTYP(K+KOFFTY)
         IF (IBTYPK.EQ.JBCNDX) THEN
            CALL READT(LURSP3,KZCONF,WRK)
            IF (IOFFC.EQ.0) THEN
               CALL DAXPY(KZCONF,EIVEC((2*K-1),IBNUM),
     *                    WRK,1,BVECS,1)
            ELSE
               CALL DAXPY(KZCONF,EIVEC(2*K,IBNUM),
     *                 WRK,1,BVECS,1)
            END IF
         ELSE
            READ(LURSP3)
         ENDIF
 200  CONTINUE
C
      IF (IPRRSP .GE. 35) THEN
         WRITE(LUPRI,*)'  SOLUTION VECTOR FOR OLSEN ALGORITHM'
         CALL OUTPUT(BVECS,1,KZCONF,1,1,KZCONF,1,1,LUPRI)
      END IF
C
C     END OF RSPCVE
C
      CALL QEXIT('RSPCVE')
      RETURN
      END
C  /* Deck phpdsk */
      SUBROUTINE PHPDSK(A,LPHPMX,PHPWRT)
C
C E(2) DIAGONAL + PHP INFORMATION AS LAST RECORD ON LURSP4
C  PHPWRT = .TRUE.  write information
C         = .false. read information
C
#include "implicit.h"
#include "dummy.h"
      DIMENSION A(*)
C
#include "infrsp.h"
#include "wrkrsp.h"
#include "infpri.h"
#include "inftap.h"
      LOGICAL PHPWRT
C
      CALL QENTER('PHPDSK')
C
      REWIND (LURSP4)
      JRSP4 = 0
C
C *** now SKIP diagonals of E[2]-matrix needed for RSPLIN
C
      IF (KZCONF .GT. 0) THEN
         READ (LURSP4)
      END IF
      IF (KZWOPT .GT. 0) THEN
         READ (LURSP4)
      END IF
C
C *** now SKIP diagonals of S[2]-matrix needed for RSPTRN
      IF (KZWOPT .GT. 0) THEN
         READ (LURSP4)
      END IF
C
C    read or write diagonal phpinformation
C
      IF (PHPWRT) THEN
         CALL WRITT(LURSP4,LPHPMX,A)
      ELSE
          CALL READT(LURSP4,LPHPMX,A)
      END IF
C
C *** END OF PHPDSK
C
      CALL QEXIT('PHPDSK')
      RETURN
      END
!  -- end of rsp/rspc82.F --
