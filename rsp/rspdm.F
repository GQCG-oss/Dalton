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
C  /* Deck rspdm */
      SUBROUTINE RSPDM(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR, RHO1,RHO2,
     *                 ISPIN1,ISPIN2,TDM,NORHO2,XINDX,WORK,
     *                 KFREE,LFREE)
C
C     CONSTRUCT ONE (RHO1) AND  TWO (RHO2) PARTICLE DENSITY
C     MATRICES TO BE USED IN RESPONSE CALCULATION.
C
C  OUTPUT:
C
C  RHO1(KL)   = (L/ E(KL) /R)
C
C  RHO2       :  TWO-BODY TRANSITION DENSITY MATRIX FOR STATE
C               (L/ AND /R) PACKED WITH SQUARE DISTRIBUTIONS
C
C  RHO2(IJKL) =    RHO2[NASHT,NASHT,NASHT,NASHT]
C
C  RHO2(IJKL) =    (L/ E(IJ,KL) - DELTA(JK) E(IL) /R )
C
#include <implicit.h>
C
      DIMENSION CL(*), CR(*),RHO1(NASHT,*)
      DIMENSION RHO2(NASHT,NASHT,NASHT,*)
      DIMENSION XINDX(*),WORK(*)
C
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0 , SMALL = 1.0D-8)
C
C Used from common blocks:
C   INFINP : LSYM,FLAG(*)
C   INFORB : NASHT,N2ASHX,NNASHX
C   INFVAR : ?? not used ?? /920920-hjaaj
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
C
C
      INTEGER RHOTYP
      LOGICAL USEPSM, CSFEXP, TDM, NORHO2
      CALL QENTER('RSPDM')
C
C     Treat NASHT .eq. 1 as a special case
C
      IF (NASHT .EQ. 0) GO TO 9999
      IF (NASHT .EQ. 1) THEN
         IF ( (ILSYM.EQ.IRSYM) .AND. (ABS(CL(1)-D1).LT.SMALL)
     *                         .AND. (ABS(CR(1)-D1).LT.SMALL) ) THEN
            RHO1(1,1) = D1
C           RHO1 = 1.0 both for ISPIN1 = 0 and ISPIN1 = 1
         ELSE
            RHO1(1,1) = D0
         END IF
         IF (.NOT. NORHO2)
     *   RHO2(1,1,1,1) = D0
         GO TO 9999
      ELSE IF (HSROHF) THEN
         IF ( (ILSYM.EQ.IRSYM) .AND. (ABS(CL(1)-D1).LT.SMALL)
     *                         .AND. (ABS(CR(1)-D1).LT.SMALL) ) THEN
            CALL DUNIT(RHO1,NASHT)
         ELSE
            CALL DZERO(RHO1,N2ASHX)
         END IF
         IF (.NOT. NORHO2) THEN
            DO I=1,NASHT
               DO J=1,NASHT
                  RHO2(I,I,J,J)=D1
                  RHO2(I,J,J,I)=RHO2(I,J,J,I)-D1
               END DO
            END DO
         END IF
         GO TO 9999
      END IF
C
C     Set RHOTYP and USEPSM
C     RHOTYP = 2 for non-symmetrized density matrix
C     (PV(ij,kl) usually .ne. PV(ij,lk))
C     USEPSM tells densid if it is allowed to use permutation symmetry
C     for Ms = 0.
C
      RHOTYP = 2
      IF (FLAG(66)) THEN
         USEPSM = .FALSE.
      ELSE
         USEPSM = .TRUE.
      END IF
C
C
      CSFEXP = .NOT.FLAG(27)
C
      IF ( IPRRSP.GT.65 ) THEN
         WRITE(LUPRI,'(/A)')' ***RSPDM BEFORE CALLING DENSID'
         WRITE(LUPRI,'(/A,/2I6,2I4,4I8)')
     *     ' ILSYM IRSYM ISPIN1/2  NCLDIM  NCRDIM   KFREE   LFREE:',
     *       ILSYM,IRSYM,ISPIN1,ISPIN2,NCLDIM,NCRDIM,KFREE,LFREE
      END IF
      IF ( IPRRSP.GT.120 ) THEN
         WRITE(LUPRI,'(/A)')' *RSPDM* LEFT REFERENCE VECTOR'
         CALL OUTPUT(CL,1,NCLDIM,1,1,NCLDIM,1,1,LUPRI)
         WRITE(LUPRI,'(/A)')' *RSPDM* RIGHT REFERENCE VECTOR'
         CALL OUTPUT(CR,1,NCRDIM,1,1,NCRDIM,1,1,LUPRI)
      END IF
      IF (.NOT. NORHO2)
     *CALL SETVEC(RHO2,D0,N2ASHX*N2ASHX)
      CALL SETVEC(RHO1,D0,N2ASHX)
      CALL DENSID(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR,RHO1,RHO2,
     *            RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,TDM,
     *            XINDX,WORK,KFREE,LFREE)
C     CALL DENSID(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR,RHO1,RHO2,
C    *            RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,TDM,
C    *            XNDXCI,WORK,KFREE,LFREE)
C
      IF ( IPRRSP.GT.90 ) THEN
         WRITE(LUPRI,'(/A)')' *** LEAVING RSPDM'
         WRITE(LUPRI,'(/A,/2I6,2I4,4I8)')
     *     ' ILSYM IRSYM ISPIN1/2  NCLDIM  NCRDIM   KFREE   LFREE:',
     *       ILSYM,IRSYM,ISPIN1,ISPIN2,NCLDIM,NCRDIM,KFREE,LFREE
      END IF
      IF ( IPRRSP.GT.90 ) THEN
         WRITE (LUPRI,1100)
         CALL OUTPUT(RHO1,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
      ENDIF
 1100 FORMAT(/' RSPDM: One-el. density matrix, active part, MO-basis')
C
      IF ( (.NOT. NORHO2 ) .AND. IPRRSP .GE. 120 ) THEN
         WRITE(LUPRI,'(/A)')
     *      ' RSPDM: Two-body density matrix RHO2 (non-zero elements):'
         DO 240 L = 1, NASHT
            DO 240 K = 1, NASHT
               WRITE(LUPRI,'(/A,2I5)') ' RHO2 distribution K,L:',K,L
               CALL OUTPUT(RHO2(1,1,K,L),1,NASHT,1,NASHT,
     *                     NASHT,NASHT,1,LUPRI)
  240    CONTINUE
      END IF
C
 9999 CONTINUE
      CALL QEXIT('RSPDM')
C     ... end of RSPDM.
      END
C  /* Deck rsptdm */
      SUBROUTINE RSPTDM(NCSIM,ILRESY,IRSYM,NCLREF,NCRDIM,CLREF,
     *                 CR, RHO1,RHO2, ISPIN1,ISPIN2,TDM,NORHO2,
     *                 XINDX,WORK,KFREE,LFREE)
C
C     CONSTRUCT ONE (RHO1) AND  TWO (RHO2) PARTICLE TRANSITION DENSITY
C     MATRICES TO BE USED IN RESPONSE CALCULATION.
C
C  OUTPUT:
C
C  RHO1(KL)   = (CLREF/ E(KL) /-CR)
C
C  RHO2       :  TWO-BODY TRANSITION DENSITY MATRIX FOR STATE
C               (CLREF/ AND /CR) PACKED WITH SQUARE DISTRIBUTIONS
C
C  RHO2(IJKL) =    RHO2[NASHT,NASHT,NASHT,NASHT]
C
C  RHO2(IJKL) =    (CLREF/ E(IJ,KL) - DELTA(JK) E(IL) /-CR )
C
#include <implicit.h>
C
      DIMENSION CLREF(*), CR(KZCONF,*),RHO1(NASHDI,NASHDI,*)
      DIMENSION RHO2(NASHDI*NASHDI,NASHDI*NASHDI,*)
      DIMENSION XINDX(*),WORK(*)
C
#include <priunit.h>
#include <wrkrsp.h>
#include <infdim.h>
#include <infpri.h>
#include <infrsp.h>
#include <inforb.h>
C
      PARAMETER ( DM1 = -1.0D0 )
      LOGICAL  TDM, NORHO2
C
      CALL QENTER('RSPTDM')
      IF (SOPPA) THEN
         WRITE (LUPRI,*) 'RSPTDM FATAL ERROR: called in SOPPA calc.'
         CALL QUIT('RSPTDM FATAL ERROR: called in SOPPA calc.')
      END IF
      DO 100 ICSIM = 1,NCSIM
         IF ( IPRRSP.GT.65 ) THEN
            WRITE(LUPRI,*) '*** RSPTDM: calling RSPDM for ICSIM=',ICSIM
         END IF
         CALL RSPDM(ILRESY,IRSYM,NCLREF,NCRDIM,CLREF,CR(1,ICSIM),
     *              RHO1(1,1,ICSIM),RHO2(1,1,ICSIM),
     *              ISPIN1,ISPIN2,TDM,NORHO2,XINDX,WORK,
     *              KFREE,LFREE)
C        CALL RSPDM(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR, RHO1,RHO2,
C    *                 ISPIN1,ISPIN2,TDM,NORHO2,XINDX,WORK,
C    *                 KFREE,LFREE)
C
C TAKE CARE OF MINUS SIGN IN CR
C
         IF( .NOT.NORHO2 )
     *   CALL DSCAL(N2ASHX*N2ASHX,DM1,RHO2(1,1,ICSIM),1)
         CALL DSCAL(N2ASHX,       DM1,RHO1(1,1,ICSIM), 1)
 100  CONTINUE
C
      CALL QEXIT('RSPTDM')
      RETURN
C     ... end of RSPTDM.
      END
C  /* Deck rspgdm */
      SUBROUTINE RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *                 CL,CR,OVLAP,RHO1,RHO2,ISPIN1,ISPIN2,TDM,NORHO2,
     *                 XINDX,WRK,KFREE,LFREE,LREF)
C
C     CONSTRUCT ONE (RHO1) AND  TWO (RHO2) PARTICLE TRANSITION DENSITY
C     MATRICES TO BE USED IN RESPONSE CALCULATION.
C
C     Adapted from RSPTDM
C
C     Instead of RSPTDM, this routine can be called with a whole
C     vector, whereas RSPTDM needs to be called with the configura-
C     tion part only. Any state can be put left and right, the
C     logical LREFR determining whether we have the reference state
C     on the right hand side, LREFL determines whether we have it on the
C     right hand side.
C
C  OUTPUT: <0(L) /../0R> + <0L/../0(L)>
C  L,R between brackets may be reference state
C
C  RHO1(KL)   = (CL/ E(KL) /-CR) + (CR / E(KL) / -CL)
C
C  RHO2       :  TWO-BODY TRANSITION DENSITY MATRIX FOR STATE
C               (CL/ AND /CR) PACKED WITH SQUARE DISTRIBUTIONS
C
C  RHO2(IJKL) =    RHO2[NASHT,NASHT,NASHT,NASHT]
C
C  RHO2(IJKL) =    (CL/ E(IJ,KL) - DELTA(JK) E(IL) /-CR )
C
#include <implicit.h>
C
      DIMENSION CL(KZVARL,*), CR(KZVARR,*)
      DIMENSION RHO1(NASHDI,NASHDI,*)
      DIMENSION RHO2(NASHDI*NASHDI,NASHDI*NASHDI,*)
      DIMENSION XINDX(*),WRK(*)
C
#include <maxorb.h>
#include <priunit.h>
#include <infvar.h>
#include <inforb.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <infdim.h>
#include <qrinf.h>
C
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0 )
      LOGICAL  TDM, NORHO2, LREF
C
      IF (SOPPA) THEN
         WRITE (LUPRI,*) 'RSPGDM FATAL ERROR: called in SOPPA calc.'
         CALL QUIT('RSPGDM FATAL ERROR: called in SOPPA calc.')
      END IF
      IF ( IPRRSP .GT. 200) THEN
         WRITE(LUPRI,'(//A)') ' Left hand side vector in RSPGDM'
         IF ( LREF ) WRITE(LUPRI,'(A)') ' (Reference state)'
         CALL OUTPUT(CL,1,KZVARL,1,NSIM,KZVARR,NSIM,1,LUPRI)
         WRITE(LUPRI,'(//A)') ' Right hand side vector in RSPGDM'
         CALL OUTPUT(CR,1,KZVARR,1,NSIM,KZVARR,NSIM,1,LUPRI)
      ENDIF
C
      IF (.NOT.TDM ) THEN
         WRITE(LUPRI,'(A)') ' FATAL ERROR: Illegal call of RSPGDM;'
         WRITE(LUPRI,'(A)') ' TDM is FALSE, must be TRUE'
         CALL QUIT(' Illegal call of RSPGDM; TDM = FALSE')
      END IF
C
      OVLAP = D0
C
C     Treat the case NASHT <= 1
C
      IF ( NASHT .EQ. 0 ) RETURN
C
      IF (NCL .EQ. 0 .OR. NCR .EQ. 0 .OR. NASHT .EQ. 1) THEN
         CALL DZERO(RHO1,NSIM*N2ASHX)
         IF (.NOT. NORHO2) CALL DZERO(RHO2,NSIM*N2ASHX*N2ASHX)
         RETURN
      END IF
C
C     Do <0(L)  ..  0R>
C
         IF( LREF ) THEN
            IOFF = 1
            NDREF = MZCONF(1)
C           ... may be NCDETS instead of NCONF inf triplet
C               therefore not NCREF which always is NCONF
            IF (KZVARL .NE. NDREF ) THEN
               CALL QUIT(' Illegal call of RSPGDM: KZVARL ne NDREF')
            ENDIF
            IF (NCL .NE. NDREF ) THEN
               CALL QUIT(' Illegal call of RSPGDM: NCL ne NDREF')
            ENDIF
         ELSE
            IOFF = MZVAR(MULD2H(IREFSY,ILSYM)) + 1
         ENDIF
C
      DO 100 ISIM = 1,NSIM
         IF ( IPRRSP.GT.65 ) THEN
            WRITE(LUPRI,*) '*** RSPGDM: 1. call of RSPDM for ISIM=',ISIM
         END IF
         CALL RSPDM(ILSYM,IRSYM,NCL,NCR,CL(IOFF,ISIM),CR(1,ISIM),
     *              RHO1(1,1,ISIM),RHO2(1,1,ISIM),
     *              ISPIN1,ISPIN2,TDM,NORHO2,XINDX,WRK,
     *              KFREE,LFREE)
C        CALL RSPDM(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR, RHO1,RHO2,
C    *                 ISPIN1,ISPIN2,TDM,NORHO2,XINDX,WRK,
C    *                 KFREE,LFREE)
C
         FACT = DM1
         CALL DSCAL(N2ASHX,FACT,RHO1(1,1,ISIM),1)
         IF (.NOT.NORHO2)  THEN
            CALL DSCAL(N2ASHX*N2ASHX,FACT,RHO2(1,1,ISIM),1)
         ENDIF
         IF (ILSYM .EQ. IRSYM) THEN
            OVLAP = OVLAP - DDOT(NCL,CL(IOFF,ISIM),1,CR(1,ISIM),1)
         END IF
C
C      Do <0L  ..  0(R)>
C
         IOFF = MZVAR(MULD2H(IREFSY,IRSYM)) + 1
         KRHO1 = KFREE
         KRHO2 = KRHO1 + N2ASHX
         IF (NORHO2) THEN
            KFREE1 = KRHO2
         ELSE
            KFREE1 = KRHO2 + N2ASHX*N2ASHX
         END IF
         LFREE1 = LFREE - KFREE1
         IF (LFREE1.LT.0) CALL ERRWRK('RSPGDM',KFREE1-1,LFREE)
C
         IF ( IPRRSP.GT.65 ) THEN
            WRITE(LUPRI,*) '*** RSPGDM: 2. call of RSPDM for ISIM=',ISIM
         END IF
         CALL RSPDM(IRSYM,ILSYM,NCR,NCL,CR(IOFF,ISIM),CL(1,ISIM),
     *              WRK(KRHO1),WRK(KRHO2),
     *              ISPIN1,ISPIN2,TDM,NORHO2,XINDX,WRK,
     *              KFREE1,LFREE1)
C
C Take care of minus sign in CR
C ( in case of no LREF both terms have a CR, scale the result)
C
         IF ( LREF ) THEN
            FACT = D1
         ELSE
            FACT = DM1
         END IF
         CALL DAXPY(N2ASHX,FACT,WRK(KRHO1),1,RHO1(1,1,ISIM),1)
         IF (.NOT.NORHO2)  THEN
            CALL DAXPY(N2ASHX*N2ASHX,FACT,WRK(KRHO2),1,RHO2(1,1,ISIM),1)
         ENDIF
         IF (ILSYM .EQ. IRSYM) THEN
            OVLAP = OVLAP + FACT*DDOT(NCR,CR(IOFF,ISIM),1,CL(1,ISIM),1)
         END IF
 100  CONTINUE
C
      IF (IPRRSP .GT. 60) THEN
         WRITE (LUPRI,*) 'RSPGDM: Overlap factor =',OVLAP
      END IF
C
      RETURN
C     ... end of RSPGDM.
      END
C  /* Deck pvxdis */
      SUBROUTINE PVXDIS(K,L,PVDEN,PVX,IPVDIS)
C
C GET 2-ELECTRON DENSITY DISTRIBUTIONS OF VARIOUS TYPE FROM PVX
C
C     IPVDIS = 1  [**,KL] + { 1 - DELTA(K,L)  }  [**,LK]
C                 IN PVDEN[*,*]
C                 SQARE PACKED DISTRIBUTIONS IN PVX
C
C     IPVDIS = 2  [*K,*L] IN PVDEN[*,*]
C                 SQARE PACKED DISTRIBUTIONS IN PVX
C
C     IPVDIS = 3  [*K,L*] IN PVDEN[*,*]
C                 SQARE PACKED DISTRIBUTIONS IN PVX
C
C     IPVDIS = 4  [**,KL] IN PVDEN[*,*]
C                 SQARE PACKED DISTRIBUTIONS IN PVX
C
#include <implicit.h>
C
      DIMENSION PVDEN(NASHDI,*),PVX(NASHDI,NASHDI,NASHDI,*)
C
#include <priunit.h>
#include <infdim.h>
#include <inforb.h>
#include <infrsp.h>
C
      PARAMETER ( D1 = 1.0D0 )
      IF ( IPRRSP.GT.2000 ) THEN
         WRITE(LUPRI,'(/A)')' ********* PVXDIS **********'
         DO 2000 I=1,NASHT
            DO 2100 J=1,NASHT
            WRITE(LUPRI,'(/A,2I6)')
     *        ' TWO-BODY DENSITY DISTRIBUTION: I,J ',I,J
              CALL OUTPUT(PVX(1,1,J,I),1,NASHT,1,NASHT,
     *                    NASHT,NASHT,1,LUPRI)
 2100       CONTINUE
 2000    CONTINUE
      END IF
C
      GO TO (1,2,3,4),IPVDIS
         WRITE(LUERR,'(/A,I5)')
     *      ' PVXDIS :INCORRECT SPECIFICATION OF IPVDIS ,IPVDIS:',IPVDIS
         CALL QUIT('PVXDIS INCORRECT SPECIFICATION OF IPVDIS')
 1    CONTINUE
         CALL DCOPY(N2ASHX,PVX(1,1,K,L),1,PVDEN(1,1),1)
         IF (K.NE.L) CALL DAXPY(N2ASHX,D1,PVX(1,1,L,K),1,PVDEN(1,1),1)
      GO TO 100
 2    CONTINUE
         DO 200 I=1,NASHT
           CALL DCOPY(NASHT,PVX(1,K,I,L),1,PVDEN(1,I),1)
 200     CONTINUE
      GO TO 100
 3    CONTINUE
         DO 300 I=1,NASHT
           CALL DCOPY(NASHT,PVX(1,K,L,I),1,PVDEN(1,I),1)
 300     CONTINUE
      GO TO 100
 4    CONTINUE
         CALL DCOPY(N2ASHX,PVX(1,1,K,L),1,PVDEN(1,1),1)
      GO TO 100
 100  CONTINUE
      IF (IPRRSP.GT.150) THEN
         WRITE(LUPRI,'(/A,I5,A,2I5)') ' PVXDIS DISTRIBUTION TYPE',
     *      IPVDIS,'     DISTRIBUTION: K,L',K,L
         CALL OUTPUT(PVDEN,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
      END IF
      RETURN
      END
