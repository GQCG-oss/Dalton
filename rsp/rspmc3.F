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
#ifdef REV_LOG
===========================================================================
Revision 1.2  2000/05/24 19:04:06  hjj
new getref calls with appropriate NDREF instead of NCREF
(fixing error for triplet with CSF)
===========================================================================
#endif
      SUBROUTINE MELONE(PRPMO,IOPSYM,DEN1,OVLAP,ONETOT,IPRONE,PRPLBL)
C
C CALCULATE AVERAGE VALUE OF ONE ELECTRON OPERATOR
C
#include <implicit.h>
C
      CHARACTER*(*)PRPLBL
C
      DIMENSION PRPMO(NORBT,*),DEN1(NASHDI,*)
C
#include <priunit.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
C
      PARAMETER ( D2 = 2.0D0 , D0 = 0.0D0 )
      PARAMETER ( BIGLIM = 100000.0D0, SMLLIM = 0.01D0 )
C
         ONEACT = D0
         ONEINA = D0
         DO 50 ISYM = 1,NSYM
            JSYM = MULD2H(IOPSYM,ISYM)
            NASHI = NASH(ISYM)
            NISHI = NISH(ISYM)
            IORBI = IORB(ISYM)
            IASHI = IASH(ISYM)
            NASHJ = NASH(JSYM)
            NISHJ = NISH(JSYM)
            IORBJ = IORB(JSYM)
            IASHJ = IASH(JSYM)
            DO 80 IINAC = 1,NISHI
               ONEINA = ONEINA + PRPMO(IORBI+IINAC,IORBI+IINAC)
 80         CONTINUE
            DO 60 IA = 1,NASHI
               DO 70 JA = 1,NASHJ
                  ONEACT = ONEACT + DEN1(IASHI+IA,IASHJ+JA) *
     *                  PRPMO(IORBI+NISHI+IA,IORBJ+NISHJ+JA)
 70            CONTINUE
 60         CONTINUE
 50      CONTINUE
C
      ONEINA = ONEINA * D2 * OVLAP
      ONETOT = ONEINA + ONEACT
      IF (IPRRSP.GE.IPRONE) THEN
         IF (ABS(ONETOT) .GT. SMLLIM .AND. ABS(ONETOT) .LT. BIGLIM)
     *                                                       THEN
            WRITE(LUPRI,'(3(/5X,2A,F15.8))')
     *      PRPLBL,' INACTIVE PART:',ONEINA,
     *      PRPLBL,' ACTIVE PART  :',ONEACT,
     *      PRPLBL,' TOTAL        :',ONETOT
         ELSE
            WRITE(LUPRI,'(3(/5X,2A,1P,D15.8))')
     *      PRPLBL,' INACTIVE PART:',ONEINA,
     *      PRPLBL,' ACTIVE PART  :',ONEACT,
     *      PRPLBL,' TOTAL        :',ONETOT
         END IF
      END IF
C
C END OF MELONE
C
      RETURN
      END
      SUBROUTINE MELTWO(H1,DEN1,DEN2,OVLAP,ISYMDN,IDAX,IOPSYM,
     *                  ISPIN1,ISPIN2,ISPIN3,IKLVL,
     *                  ZYMAT1,ISYM1,ZYMAT2,ISYM2,ZYMAT3,ISYM3,
     *                  RES,CMO,WRK,LWRK)
C
C Calculation of the two-electron matrix element <L|K|R>.
C The generalized hamiltonian K can be one-index transformed 0, 1 or 2
C times according to IKLVL. (0 gives the ordinary untransformed hamiltonian.)
C
#include <implicit.h>
C
#include <priunit.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
#include <maxash.h>
#include <maxorb.h>
#include <infind.h>
C
      PARAMETER (D1=1.0D0, D2=2.0D0)
      DIMENSION WRK(*)
      DIMENSION H1(NORBT,NORBT)
      DIMENSION DEN1(NASHDI,NASHDI),DEN2(N2ASHX,N2ASHX)
      DIMENSION ZYMAT1(NORBT,NORBT),ZYMAT2(NORBT,NORBT)
      DIMENSION ZYMAT3(NORBT,NORBT)
      DIMENSION CMO(*)
C
      LOGICAL LCON,LORB
C
      CALL QENTER('MELTWO')
C
      IGRSPI = 0
C
      IF (IPRRSP.GT.200) THEN
         WRITE(LUPRI,'(/A)') 'ENTER MELTWO'
         WRITE(LUPRI,'(A,I5)') 'IDAX =',IDAX
         WRITE(LUPRI,'(A,I5)') 'IOPSYM =',IOPSYM
         WRITE(LUPRI,'(A,I5)') 'ISYMDN =',ISYMDN
         WRITE(LUPRI,*) 'OVLAP  =',OVLAP
         WRITE(LUPRI,'(A,I5)') 'IGRSPI =',IGRSPI
         WRITE(LUPRI,'(A,I5)') 'ISPIN1 =',ISPIN1
         WRITE(LUPRI,'(A,I5)') 'ISPIN2 =',ISPIN2
         WRITE(LUPRI,'(A,I5)') 'ISPIN3 =',ISPIN3
         WRITE(LUPRI,'(A,I5)') 'IKLVL =',IKLVL
         WRITE(LUPRI,'(A,I5)') 'ISYM1 =',ISYM1
         WRITE(LUPRI,'(A,I5)') 'ISYM2 =',ISYM2
         WRITE(LUPRI,'(A,I5)') 'ISYM3 =',ISYM3
         WRITE(LUPRI,'(A)') 'One-electron H1'
         CALL OUTPUT(H1,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(A)') 'One-electron density'
         CALL OUTPUT(DEN1,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         WRITE(LUPRI,'(A)') 'Two-electron density'
         CALL PRIAC2(DEN2,NASHT,LUPRI)
      END IF
C
C  KONE = 1 is needed in the call to RSPFX. Putting the constant 1
C  in that call gives a MEMGET error in NXTH2X.
C  28/3-1996-ekd
C
      KONE   = 1
C
      KH2AX  = 1
      LH2AX  = N2ASHX * N2ASHX
      IF (DIROIT) THEN
         IF (IKLVL.EQ.2) LH2AX = LH2AX * 2
         IF (IKLVL.EQ.3) LH2AX = LH2AX * 4
      END IF
      KFA    = KH2AX  + LH2AX
      KFI    = KFA    + N2ORBX
      KQA    = KFI    + N2ORBX
      KQB    = KQA    + NORBT  * NASHDI
      KPVD   = KQB    + NORBT  * NASHDI
      KH2    = KPVD   + N2ASHX * N2ASHX
      KH2X   = KH2    + N2ORBX
      KFREE  = KH2X   + N2ORBX
      LFREE  = LWRK   - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('MELTWO',KFREE-1,LWRK)
C
      CALL DCOPY(N2ORBX,H1,1,WRK(KFI),1)
      CALL DZERO(WRK(KFA),N2ORBX)
      IF (NASHT .GT. 0) THEN
         CALL DZERO(WRK(KQA),NASHT*NORBT)
         CALL DZERO(WRK(KQB),NASHT*NORBT)
         CALL DZERO(WRK(KH2AX),LH2AX)
      END IF
C
      IF (DIROIT) THEN
         INTSYM = 1
      ELSE
         INTSYM = IOPSYM
      END IF
C
C     We only need the transformed inactive Fock matrix and the
C     tranformed active two-electron integrals.
C     However, FA, QA, and QB must also be initialized because
C     they are referenced in RSPFXD and this has given
C     floating point exception on DEC alpha. /July 2001 HJAaJ
C
      LCON = .TRUE.
      LORB = .FALSE.
      IF (DIROIT) THEN
         CALL RSPFXD(WRK(KFI),WRK(KFA),WRK(KQA),WRK(KQB),WRK(KH2AX),
     *                  IDAX,INTSYM,ISYMDN,DEN1,DEN2,WRK(KPVD),
     *                  WRK(KH2),WRK(KH2X),WRK(KFREE),KONE,LFREE,
     *                  LCON,LORB,IPRRSP,IGRSPI,ISPIN1,ISPIN2,ISPIN3,
     *                  IKLVL,ZYMAT1,ISYM1,ZYMAT2,ISYM2,ZYMAT3,ISYM3)
      ELSE
         CALL RSPFX(WRK(KFI),WRK(KFA),WRK(KQA),WRK(KQB),WRK(KH2AX),IDAX,
     *              INTSYM,ISYMDN,DEN1,DEN2,WRK(KPVD),WRK(KH2X),
     *              WRK(KFREE),KONE,LFREE,LCON,LORB,IPRRSP) 
      END IF
C     
C     The transformed two-electron integrals do not have the factor 1/2
C     that MEL2 assumes. If IKLVL = 2 we have two densities. For the
C     spin free case they are equal and we just add them.
C
      IF (DIROIT.AND.(IKLVL.GT.0)) THEN
        IF (IKLVL.EQ.2) THEN
           CALL DAXPY(N2ASHX*N2ASHX,D1,WRK(KH2AX+N2ASHX*N2ASHX),1,
     *                WRK(KH2AX),1)
        END IF
        CALL DSCAL(N2ASHX*N2ASHX,D2,WRK(KH2AX),1)
      END IF
C
      CALL MEL2(H1,WRK(KFI),WRK(KH2AX),IOPSYM,DEN1,DEN2,OVLAP,RES)
C
      CALL QEXIT('MELTWO')
C
      RETURN
      END
      SUBROUTINE MEL2(H1,FI,H2A,IOPSYM,DEN1,DEN2,OVLAP,RES)
C
C Calculate the two-electron matrix element <L|H|R>.
C
C RES = OVLAP * sum(i) [FI(i,i) + H1(i,i)] + sum(x,y) FI(x,y) DEN1(x,y)
C         + 1/2 * sum(x,y,u,v) H2A(x,y,u,v) DEN2(x,y,u,v)
C
#include <implicit.h>
C
#include <priunit.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
#include <maxash.h>
#include <maxorb.h>
#include <infind.h>
C
      PARAMETER ( D0 = 0.0D0, DH = 0.5D0 )
C
      DIMENSION H1(NORBT,NORBT),FI(NORBT,NORBT)
      DIMENSION H2A(N2ASHX*N2ASHX)
      DIMENSION DEN1(NASHDI,NASHDI),DEN2(N2ASHX*N2ASHX)
C
      IF (IPRRSP.GT.200) THEN
         WRITE(LUPRI,'(/A)') 'ENTER MEL2'
         WRITE(LUPRI,'(A,I5)') 'IOPSYM =',IOPSYM
         WRITE(LUPRI,*) 'OVLAP  =',OVLAP
         WRITE(LUPRI,'(A)') 'One-electron H1'
         CALL OUTPUT(H1,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(A)') 'Inactive Fock matrix FI'
         CALL OUTPUT(FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(A)') 'One-electron density'
         CALL OUTPUT(DEN1,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         WRITE(LUPRI,'(A)') 'Two-electron density'
         CALL PRIAC2(DEN2,NASHT,LUPRI)
         WRITE(LUPRI,'(A)') 'Two-electron H2A'
         CALL PRIAC2(H2A,NASHT,LUPRI)
         IF (DIROIT) THEN
            CALL PRIAC2(H2A(N2ASHX*N2ASHX+1),NASHT,LUPRI)
         END IF
      END IF
C
      RES = D0
      DO 10 ISYM = 1,NSYM
         NISHI = NISH(ISYM)
         IORBI = IORB(ISYM)
         DO 20 IINAC = 1,NISHI
            RES = RES + FI(IORBI+IINAC,IORBI+IINAC)
            RES = RES + H1(IORBI+IINAC,IORBI+IINAC)
   20    CONTINUE
   10 CONTINUE
      RES = RES * OVLAP
C
      DO 30 ISYM = 1,NSYM
         JSYM  = MULD2H(IOPSYM,ISYM)
         NASHI = NASH(ISYM)
         NISHI = NISH(ISYM)
         IORBI = IORB(ISYM)
         IASHI = IASH(ISYM)
         NASHJ = NASH(JSYM)
         NISHJ = NISH(JSYM)
         IORBJ = IORB(JSYM)
         IASHJ = IASH(JSYM)
         DO 40 IAC = 1,NASHI
         DO 50 JAC = 1,NASHJ
            RES = RES + FI(IORBI+NISHI+IAC,IORBJ+NISHJ+JAC)*
     *                DEN1(IASHI+IAC,IASHJ+JAC)
   50    CONTINUE
   40    CONTINUE
   30 CONTINUE
C
      RES = RES + DH*DDOT(N2ASHX*N2ASHX,H2A,1,DEN2,1)
C
      IF (IPRRSP.GT.10) THEN
         WRITE(LUPRI,'(/A,F15.8)') ' Result in MEL2:',RES
      END IF
C
      RETURN
      END
      SUBROUTINE S4DRV(KZYVR,KZYV1,KZYV2,KZYV3,
     *                 IGRSYM,ISYMV1,ISYMV2,ISYMV3,
     *                 S4TRS,VEC1,VEC2,VEC3,
     *                 UDV,ZYMAT,DEN1,XINDX,MJWOP,WRK,LWRK)
C
C     Drive the computation of S[4] times three vectors
C
#include <implicit.h>
#include <infdim.h>
#include <inforb.h>
#include <maxorb.h>
#include <maxash.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <infhyp.h>
#include <infvar.h>
#include <infind.h>
#include <infpri.h>
#include <qrinf.h>
C
      PARAMETER( D1= 1.0D0, D6=6.0D0, DH=0.5D0 )
C
      DIMENSION WRK(*)
      DIMENSION S4TRS(KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION VEC1(KZYV1), VEC2(KZYV2), VEC3(KZYV3)
      DIMENSION ZYMAT(NORBT,NORBT)
      DIMENSION XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION DEN1(NASHDI,NASHDI)
C
      LOGICAL LORB, LCON, LREF, TDM, NORHO2
C
C     Layout some workspace
C
      KCREF  = 1
      KRES   = KCREF + MZCONF(1)
      KFREE  = KRES  + KZYVR
      LFREE  = LWRK  - KFREE + 1
      IF (LFREE.LT.0) CALL ERRWRK('S4DRV 1',KFREE-1,LWRK)
C
C
C     Initialize variables
C
      TDM    = .TRUE.
      NORHO2 = .TRUE.
      NSIM = 1
      ISPIN = 0
C
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
      CALL SCASE1(KZYVR,KZYV1,KZYV2,KZYV3,
     *                 IGRSYM,ISYMV1,ISYMV2,ISYMV3,
     *                 S4TRS,VEC1,VEC2,VEC3,WRK(KCREF),
     *                 UDV,DEN1,XINDX,MJWOP,WRK(KFREE),LFREE)
C      
      CALL SCASE1(KZYVR,KZYV2,KZYV1,KZYV3,
     *                 IGRSYM,ISYMV2,ISYMV1,ISYMV3,
     *                 S4TRS,VEC2,VEC1,VEC3,WRK(KCREF),
     *                 UDV,DEN1,XINDX,MJWOP,WRK(KFREE),LFREE)
C      
      CALL SCASE1(KZYVR,KZYV3,KZYV1,KZYV2,
     *                 IGRSYM,ISYMV3,ISYMV1,ISYMV2,
     *                 S4TRS,VEC3,VEC1,VEC2,WRK(KCREF),
     *                 UDV,DEN1,XINDX,MJWOP,WRK(KFREE),LFREE)
C      
      CALL SCASE2(KZYVR,KZYV1,KZYV2,KZYV3,
     *                 IGRSYM,ISYMV1,ISYMV2,ISYMV3,
     *                 S4TRS,VEC1,VEC2,VEC3,WRK(KCREF),
     *                 UDV,ZYMAT,DEN1,XINDX,MJWOP,WRK(KFREE),LFREE)
C
      CALL SCASE2(KZYVR,KZYV1,KZYV3,KZYV2,
     *                 IGRSYM,ISYMV1,ISYMV3,ISYMV2,
     *                 S4TRS,VEC1,VEC3,VEC2,WRK(KCREF),
     *                 UDV,ZYMAT,DEN1,XINDX,MJWOP,WRK(KFREE),LFREE)
C
      IF (MZWOPT(ISYMV1).EQ.0 .OR.
     *    MZCONF(ISYMV2).EQ.0 .OR.
     *    MZCONF(ISYMV3).EQ.0) GOTO 4000
C
C     /   <02L| [qj,K(k1)] |03R>  + <03L| [qj,K(k1)] |02R>  \
C     |                       0                             |
C     |   <02L| [qj+,K(k1)] |03R> + <03L| [qj+,K(k1)] |02R> |
C     \                       0                             /
C
C     Construct <03L|..|02R> + <02L|..|03R> density
C
      ILSYM  = MULD2H(IREFSY,ISYMV2)
      IRSYM  = MULD2H(IREFSY,ISYMV3)
      NCL    = MZCONF(ISYMV2)
      NCR    = MZCONF(ISYMV3)
      KZVARL = MZYVAR(ISYMV2)
      KZVARR = MZYVAR(ISYMV3)
      LREF   = .FALSE.
      ISYMDN = MULD2H(ILSYM,IRSYM)
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *         VEC2,VEC3,OVLAP,DEN1,DUMMY,ISPIN,ISPIN,TDM,NORHO2,
     *         XINDX,WRK,KFREE,LFREE,LREF)
C
C     Make the gradient
C
      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,ZYMAT,MJWOP)
C
      IF ( MZWOPT(IGRSYM) .GT. 0 ) THEN
         CALL ORBSX(NSIM,IGRSYM,KZYVR,S4TRS,ZYMAT,OVLAP,ISYMDN,
     *              DEN1,MJWOP,WRK(KFREE),LFREE)
      END IF
      IF (ISYMV2.NE.ISYMV3) GOTO 4000
C
C     /   <0| [qj ,K(k1)] |0> \
C     | 1/2<j| K(k1) |0>      | * ( S(2)S(3)' + S(2)'S(3) )
C     |   <0| [qj+,K(k1)] |0> |
C     \ -1/2<0| K(k1) |j>     /
C
      NZCONF = MZCONF(ISYMV2)
      NZVAR  = MZVAR(ISYMV2)
      F1 = DDOT(NZCONF,VEC2,1,VEC3(NZVAR+1),1) +
     *     DDOT(NZCONF,VEC3,1,VEC2(NZVAR+1),1)
C
      ISYMDN = 1
      OVLAP  = D1
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF(ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      LREF = .TRUE.
      IF (LCON) THEN
         CALL DZERO(WRK(KRES),KZYVR)
         CALL CONSX(NSIM,KZYVR,IGRSYM,ZYMAT,WRK(KCREF),
     *              MZCONF(1),MZCONF(1),IREFSY,MZCONF(IGRSYM),ISYMST,
     *              LREF,WRK(KRES),XINDX,WRK(KFREE),LFREE)
         CALL DSCAL(KZYVR,DH,WRK(KRES),1)
         CALL DAXPY(KZYVR,F1,WRK(KRES),1,S4TRS,1)
      END IF
      IF (LORB) THEN
         CALL DZERO(WRK(KRES),KZYVR)
         CALL ORBSX(NSIM,IGRSYM,KZYVR,WRK(KRES),ZYMAT,OVLAP,ISYMDN,
     *              UDV,MJWOP,WRK(KFREE),LFREE)
         CALL DAXPY(KZYVR,F1,WRK(KRES),1,S4TRS,1)
      END IF
C
 4000 CONTINUE
C
      CALL SCASE3(KZYVR,KZYV1,KZYV2,KZYV3,
     *                 IGRSYM,ISYMV1,ISYMV2,ISYMV3,
     *                 S4TRS,VEC1,VEC2,VEC3,WRK(KCREF),
     *                 UDV,ZYMAT,DEN1,XINDX,MJWOP,WRK(KFREE),LFREE)
C
      CALL SCASE3(KZYVR,KZYV1,KZYV3,KZYV2,
     *                 IGRSYM,ISYMV1,ISYMV3,ISYMV2,
     *                 S4TRS,VEC1,VEC3,VEC2,WRK(KCREF),
     *                 UDV,ZYMAT,DEN1,XINDX,MJWOP,WRK(KFREE),LFREE)
C
      IF (MZWOPT(ISYMV1).EQ.0 .OR.
     *    MZWOPT(ISYMV2).EQ.0 .OR.
     *    MZWOPT(ISYMV3).EQ.0) RETURN
C
C     / <0| [qj ,K(k1,k2,k3)+K(k1,k3,k2)] |0> \
C     | <j| K(k1,k2,k3)+K(k1,k3,k2) |0>       | * 1/6
C     | <0| [qj+,K(k1,k2,k3)+K(k1,k3,k2)] |0> |
C     \ -<0| K(k1,k2,k3)+K(k1,k3,k2) |j>      /
C
      ISYMDN = 1
      OVLAP  = D1
C
C     Make the gradient
C
      ISYMV  = IREFSY
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      LREF = .TRUE.
      NZYVEC = MZCONF(1)
      NZCVEC = MZCONF(1)
      CALL TRZYM2(VEC1,VEC2,VEC3,KZYV1,KZYV2,KZYV3,
     *           ISYMV1,ISYMV2,ISYMV3,ZYMAT,MJWOP,
     *           WRK(KFREE),LFREE)
      CALL DSCAL(NORBT*NORBT,1/D6,ZYMAT,1)
      CALL RSP1GR(NSIM,KZYVR,IDUM,ISPIN,IGRSYM,ISPIN,ISYMV,S4TRS,
     *            WRK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,ZYMAT,
     *            XINDX,MJWOP,WRK(KFREE),LFREE,LORB,LCON,LREF)
C
      CALL TRZYM2(VEC1,VEC3,VEC2,KZYV1,KZYV3,KZYV2,
     *           ISYMV1,ISYMV3,ISYMV2,ZYMAT,MJWOP,
     *           WRK(KFREE),LFREE)
      CALL DSCAL(NORBT*NORBT,1/D6,ZYMAT,1)
      CALL RSP1GR(NSIM,KZYVR,IDUM,ISPIN,IGRSYM,ISPIN,ISYMV,S4TRS,
     *            WRK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,ZYMAT,
     *            XINDX,MJWOP,WRK(KFREE),LFREE,LORB,LCON,LREF)
C
      RETURN
      END
      SUBROUTINE SCASE1(KZYVR,KZYV1,KZYV2,KZYV3,
     *                 IGRSYM,ISYMV1,ISYMV2,ISYMV3,
     *                 S4TRS,VEC1,VEC2,VEC3,CREF,
     *                 UDV,DEN1,XINDX,MJWOP,WRK,LWRK)
C      
C
C     /  <01L| qj |0> + <0| qj |01R>  \
C     |               0               | * -2/3*(S(2)S(3)' + S(2)'S(3))
C     |  <01L| qj+|0> + <0| qj+|01R>  |
C     \               0               /
C      
C      
#include <implicit.h>
#include <infdim.h>
#include <inforb.h>
#include <maxorb.h>
#include <maxash.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <infhyp.h>
#include <infvar.h>
#include <infind.h>
#include <infpri.h>
#include <qrinf.h>
C
      PARAMETER( D2=2.0D0, D3=3.0D0 )
C
      DIMENSION WRK(*)
      DIMENSION S4TRS(KZYVR)
      DIMENSION VEC1(KZYV1), VEC2(KZYV2), VEC3(KZYV3)
      DIMENSION XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION DEN1(NASHDI,NASHDI)
      DIMENSION CREF(*), MJWOP(2,MAXWOP,8)
C
      LOGICAL LORB, LCON, LREF, TDM, NORHO2
C
      IF (MZCONF(ISYMV1).EQ.0 .OR.
     *    MZCONF(ISYMV2).EQ.0 .OR.
     *    MZCONF(ISYMV3).EQ.0) RETURN
C
      IF (ISYMV2.NE.ISYMV3) RETURN
C
C     Initialize variables
C
      TDM    = .TRUE.
      NORHO2 = .TRUE.
      NSIM = 1
      ISPIN = 0
      KONE = 1
C
C     Construct the density matrix <01L|..|0> + <0|..|01R>
C
      ILSYM  = IREFSY
      IRSYM  = MULD2H(IREFSY,ISYMV1)
      NCL    = MZCONF(1)
      NCR    = MZCONF(ISYMV1)
      KZVARL = MZCONF(1)
      KZVARR = MZYVAR(ISYMV1)
      LREF   = .TRUE.
      ISYMDN = MULD2H(ILSYM,IRSYM)
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *         CREF,VEC1,OVLAP,DEN1,DUMMY,ISPIN,ISPIN,TDM,
     *         NORHO2,XINDX,WRK,KONE,LWRK,LREF)
C
C     Put the factor into the density matrix
C
      NZCONF = MZCONF(ISYMV2)
      NZVAR  = MZVAR(ISYMV2)
      F1 = DDOT(NZCONF,VEC2,1,VEC3(NZVAR+1),1) +
     *     DDOT(NZCONF,VEC3,1,VEC2(NZVAR+1),1)
      F1 = -F1*D2/D3
      CALL DSCAL(NASHT*NASHT,F1,DEN1,1)
C
C     Make the gradient on the fly:
C
      NZCONF = MZCONF(IGRSYM)
      NYCONF = MZCONF(IGRSYM) + MZVAR(IGRSYM)
C
      DO 500 IC = 1, MZWOPT(IGRSYM)
         K = MJWOP(1,IC,IGRSYM)
         L = MJWOP(2,IC,IGRSYM)
         ITYPK = IOBTYP(K)
         ITYPL = IOBTYP(L)
         IF(ITYPK.EQ.JTACT .AND. ITYPL. EQ. JTACT) THEN
            NWK  = ISW(K) - NISHT
            NWL  = ISW(L) - NISHT
               S4TRS(NZCONF+IC) = S4TRS(NZCONF+IC)
     *                                      + DEN1(NWK,NWL)
               S4TRS(NYCONF+IC) = S4TRS(NYCONF+IC)
     *                                      + DEN1(NWL,NWK)
         END IF
500   CONTINUE
C
      RETURN
      END
      SUBROUTINE SCASE2(KZYVR,KZYV1,KZYV2,KZYV3,
     *                 IGRSYM,ISYMV1,ISYMV2,ISYMV3,
     *                 S4TRS,VEC1,VEC2,VEC3,CREF,
     *                 UDV,ZYMAT,DEN1,XINDX,MJWOP,WRK,LWRK)
C
C
#include <implicit.h>
#include <infdim.h>
#include <inforb.h>
#include <maxorb.h>
#include <maxash.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <infhyp.h>
#include <infvar.h>
#include <infind.h>
#include <infpri.h>
#include <qrinf.h>
C
      PARAMETER( DH=0.5D0, D0=0.0D0, D1=1.0D0, D6=6.0D0 )
C
      DIMENSION WRK(*)
      DIMENSION S4TRS(KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION VEC1(KZYV1), VEC2(KZYV2), VEC3(KZYV3)
      DIMENSION XINDX(*)
      DIMENSION ZYMAT(NORBT,NORBT)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION DEN1(NASHDI,NASHDI)
      DIMENSION CREF(*)
C
      LOGICAL LORB, LCON, LREF, TDM, NORHO2
C
C     Initialize variables
C
      TDM    = .TRUE.
      NORHO2 = .TRUE.
      NSIM = 1
      ISPIN = 0
      IPRONE = 75
      F = D0
      F1 = D0
      F2 = D0
      KONE = 1
C
      IF (MZCONF(ISYMV1).EQ.0 .OR.
     *    MZCONF(ISYMV2).EQ.0 .OR.
     *    MZCONF(ISYMV3).EQ.0) GOTO 1000
C
C     /   0     \
C     |  Sj(1)  | * 1/6*S(2)S(3)'
C     |   0     |
C     \ -Sj(1)' /
C
      IF (ISYMV2.EQ.ISYMV3) THEN
         NZCONF = MZCONF(ISYMV2)
         NZVAR  = MZVAR(ISYMV2)
         FACT = 1/D6*DDOT(NZCONF,VEC2,1,VEC3(NZVAR+1),1)
         NZCONF = MZCONF(ISYMV1)
         NZVAR  = MZVAR(ISYMV1)
         CALL DAXPY(NZCONF,FACT,VEC1,1,S4TRS,1)
         CALL DAXPY(NZCONF,-FACT,VEC1(NZVAR+1),1,S4TRS(NZVAR+1),1)
      END IF
C
      IF (ISYMV1.NE.ISYMV3) RETURN
C
C     /  0            \
C     | (F+F1)*Sj(2)  |    6*F1 = S(3)'S(1) - 2*S(1)'S(3)
C     |  0            |
C     \ (F+F2)*Sj(2)' /    6*F2 = 2*S(3)'S(1) - S(1)'S(3)
C
C     F = 1/2*<03L|K(k1)|0> + <0|K(k1)|03R> + 1/2*<0|K(k1,k3)|0>
C
      NZCONF = MZCONF(ISYMV1)
      NZVAR  = MZVAR(ISYMV1)
      F1 = DDOT(NZCONF,VEC1,1,VEC3(NZVAR+1),1) -
     *     2*DDOT(NZCONF,VEC3,1,VEC1(NZVAR+1),1)
      F2 = 2*DDOT(NZCONF,VEC1,1,VEC3(NZVAR+1),1) -
     *     DDOT(NZCONF,VEC3,1,VEC1(NZVAR+1),1)
      F1 = F1/D6
      F2 = F2/D6
C
 1000 CONTINUE
      IF (MZCONF(ISYMV3).EQ.0 .OR. IGRSYM.NE.ISYMV2) RETURN
C
      IF (MZWOPT(ISYMV1).GT.0 .AND. MZCONF(ISYMV3).GT.0) THEN
C
         CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,ZYMAT,MJWOP)
C
C        F3R = <0|K(k1)|-03R>
C
         ILSYM  = IREFSY
         IRSYM  = MULD2H(IREFSY,ISYMV3)
         NCL    = MZCONF(1)
         NCR    = MZCONF(ISYMV3)
         IOFF   = 1
         CALL DZERO(DEN1,NASHT*NASHT)
         CALL RSPDM(ILSYM,IRSYM,NCL,NCR,CREF,VEC3(IOFF),
     *              DEN1,DUMMY,ISPIN,ISPIN,TDM,NORHO2,XINDX,WRK,
     *              KONE,LWRK)
         OVLAP = D0
         IF (ILSYM.EQ.IRSYM)
     *         OVLAP = DDOT(NCL,CREF,1,VEC3(IOFF),1)
C   
         CALL MELONE(ZYMAT,ISYMV1,DEN1,OVLAP,F3R,IPRONE,'F3R in SCASE2')
C
C        F3L = <03L|K(k1)|0>
C
         ILSYM  = MULD2H(IREFSY,ISYMV3)
         IRSYM  = IREFSY
         NCL    = MZCONF(ISYMV3)
         NCR    = MZCONF(1)
         IOFF   = MZVAR(ISYMV3) + 1
         CALL DZERO(DEN1,NASHT*NASHT)
         CALL RSPDM(ILSYM,IRSYM,NCL,NCR,VEC3(IOFF),CREF,
     *              DEN1,DUMMY,ISPIN,ISPIN,TDM,NORHO2,XINDX,WRK,
     *              KONE,LWRK)
         OVLAP = D0
         IF (ILSYM.EQ.IRSYM)
     *      OVLAP = DDOT(NCL,CREF,1,VEC1(IOFF),1)
C
         CALL MELONE(ZYMAT,ISYMV1,DEN1,OVLAP,F3L,IPRONE,'F3L in SCASE2')
         F = DH*F3L - F3R
      END IF
C
      IF (MZWOPT(ISYMV1).GT.0 .AND. MZWOPT(ISYMV3).GT.0) THEN
C
C        FACT = <0|K(k1,k3)|0>
C
         CALL TRZYMT(NSIM,VEC3,VEC1,KZYV3,KZYV1,ISYMV3,ISYMV1,ZYMAT,
     *            MJWOP,WRK,LWRK)
         OVLAP = D1
         CALL MELONE(ZYMAT,1,UDV,OVLAP,FACT,IPRONE,'FACT in SCASE2')
         F = F + DH*FACT
      END IF
C
      NZCONF = MZCONF(IGRSYM)
      NZVAR  = MZVAR(IGRSYM)
      CALL DAXPY(NZCONF,F+F1,VEC2,1,S4TRS,1)
      CALL DAXPY(NZCONF,F+F2,VEC2(NZVAR+1),1,S4TRS(NZVAR+1),1)
C
      RETURN
      END
      SUBROUTINE SCASE3(KZYVR,KZYV1,KZYV2,KZYV3,
     *                 IGRSYM,ISYMV1,ISYMV2,ISYMV3,
     *                 S4TRS,VEC1,VEC2,VEC3,CREF,
     *                 UDV,ZYMAT,DEN1,XINDX,MJWOP,WRK,LWRK)
C
C     /   <0| [qj,K(k1,k3)] |02R>  + <02L| [qj,K(k1,k3)] |0>  \
C     |   <j| K(k1,k3) |02R>                                  |*1/2
C     |   <0| [qj+,K(k1,k3)] |02R> + <02L| [qj+,K(k1,k3)] |0> |
C     \  -<02L| K(k1,k3) |j>                                  /
C
#include <implicit.h>
#include <infdim.h>
#include <inforb.h>
#include <maxorb.h>
#include <maxash.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <infhyp.h>
#include <infvar.h>
#include <infind.h>
#include <infpri.h>
#include <qrinf.h>
C
      PARAMETER( DH=0.5D0 )
C
      DIMENSION WRK(*)
      DIMENSION S4TRS(KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION VEC1(KZYV1), VEC2(KZYV2), VEC3(KZYV3)
      DIMENSION XINDX(*)
      DIMENSION ZYMAT(NORBT,NORBT)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION DEN1(NASHDI,NASHDI)
      DIMENSION CREF(*)
C
      LOGICAL LORB, LCON, LREF, TDM, NORHO2
C
      IF (MZWOPT(ISYMV1).EQ.0 .OR.
     *    MZCONF(ISYMV2).EQ.0 .OR.
     *    MZWOPT(ISYMV3).EQ.0) RETURN
C
C     Initialize variables
C
      TDM    = .TRUE.
      NORHO2 = .TRUE.
      NSIM = 1
      ISPIN = 0
      KONE = 1
C
C     Construct the density matrix <02L|..|0> + <0|..|02R>
C
      ILSYM  = IREFSY
      IRSYM  = MULD2H(IREFSY,ISYMV2)
      NCL    = MZCONF(1)
      NCR    = MZCONF(ISYMV2)
      KZVARL = MZCONF(1)
      KZVARR = MZYVAR(ISYMV2)
      LREF   = .TRUE.
      ISYMDN = MULD2H(ILSYM,IRSYM)
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *         CREF,VEC2,OVLAP,DEN1,DUMMY,ISPIN,ISPIN,TDM,
     *         NORHO2,XINDX,WRK,KONE,LWRK,LREF)
C
C     Put the factor into the operator
C
      CALL TRZYMT(NSIM,VEC3,VEC1,KZYV3,KZYV1,ISYMV3,ISYMV1,ZYMAT,MJWOP,
     *            WRK,LWRK)
      CALL DSCAL(NORBT*NORBT,DH,ZYMAT,1)
C
C     Make the gradient
C
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      LREF = .FALSE.
      NZYVEC = MZYVAR(ISYMV2)
      NZCVEC = MZCONF(ISYMV2)
      CALL RSP1GR(NSIM,KZYVR,IDUM,ISPIN,IGRSYM,ISPIN,ISYMV2,S4TRS,
     *            VEC2,NZYVEC,NZCVEC,OVLAP,ISYMDN,DEN1,ZYMAT,
     *            XINDX,MJWOP,WRK,LWRK,LORB,LCON,LREF)
C
      RETURN
      END
