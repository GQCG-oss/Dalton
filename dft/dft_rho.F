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
C  /* Deck getrho */
      SUBROUTINE GETRHO(DMAT,GSO,RHO,DMAGAO,DFTHRI)
C
C     T. Helgaker feb 01
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0,
     &           DP5 = 0.5D0)
#include <inforb.h>
C
      DIMENSION DMAT(NBAST,NBAST), GSO(NBAST,*), DMAGAO(NBAST)
C
#include <symmet.h>
C
C
      IF (NSYM.EQ.1) THEN
         CALL DSYMV('U',NBAST,D1,DMAT,NBAST,GSO,1,D0,DMAGAO,1)
      ELSE
         CALL DZERO(DMAGAO,NBAST)
         DO ISYM = 1, NSYM
            ISTR = IBAS(ISYM) + 1
            IEND = IBAS(ISYM) + NBAS(ISYM)
            DO J = ISTR, IEND 
               GSOJ = GSO(J,1)
               IF (DABS(GSOJ).GT.DFTHRI) THEN
                  DO I = ISTR, IEND 
                     DMAGAO(I) = DMAGAO(I) + GSOJ*DMAT(I,J)
                  END DO
               END IF
            END DO
         END DO
      END IF
      RHO   = DDOT(NBAST,GSO,1,DMAGAO,1)
      RETURN
      END      
C /* Deck dftrhh */
      SUBROUTINE DFTRHH(DMAT,DMAGAO,GAO,GAO1,GAO2,RH,RHOLAP,RHOGHG)
C
C     T. Helgaker oct 2000
C
#include <implicit.h>
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0) 
C
#include <inforb.h>
C
      DIMENSION DMAT(NBAST,NBAST), DMAGAO(NBAST),RH(3),
     &          GAO(NBAST), GAO1(NBAST,3), GAO2(NBAST,6)
C    GAO2 projected onto DMGAO
      DIMENSION GAO2PR(6)
C
      RXX = D0
      RXY = D0
      RXZ = D0
      RYY = D0
      RYZ = D0
      RZZ = D0
      DO J = 1, NBAST
         GJX  = GAO1(J,1)
         GJY  = GAO1(J,2)
         GJZ  = GAO1(J,3)
         DGIX = DDOT(NBAST,DMAT(1,J),1,GAO1(1,1),1)
         DGIY = DDOT(NBAST,DMAT(1,J),1,GAO1(1,2),1)
         DGIZ = DDOT(NBAST,DMAT(1,J),1,GAO1(1,3),1)
         RXX = RXX + D2*GJX*DGIX
         RXY = RXY + DGIX*GJY + DGIY*GJX
         RXZ = RXZ + DGIX*GJZ + DGIZ*GJX
         RYY = RYY + D2*GJY*DGIY
         RYZ = RYZ + DGIY*GJZ + DGIZ*GJY
         RZZ = RZZ + D2*GJZ*DGIZ
      END DO
      CALL DGEMV('T',NBAST,6,D2,GAO2,NBAST,DMAGAO,1,0D0,GAO2PR,1)
      RXX = RXX + GAO2PR(1)
      RXY = RXY + GAO2PR(2)
      RXZ = RXZ + GAO2PR(3)
      RYY = RYY + GAO2PR(4)
      RYZ = RYZ + GAO2PR(5)
      RZZ = RZZ + GAO2PR(6)
      RHOLAP = RXX + RYY + RZZ
      RHOGHG = RH(1)*RXX*RH(1) + RH(2)*RYY*RH(2) + RH(3)*RZZ*RH(3)
     &   + D2*(RH(1)*RXY*RH(2) + RH(1)*RXZ*RH(3) + RH(2)*RYZ*RH(3))
      END
c
c ===================================================================
c BLOCKED VERSION OF PROPERTY EVALUATION ROUTINES
c Written by Pawel Salek, closely based on the above.
c ===================================================================
      SUBROUTINE getrho_blocked_lda(DMAT,GAO,NBLOCKS,IBLOCKS,LDAIB,
     &                              TMP,NVCLEN,RHO)
c     computes  <o|dmat|o'>
c     i.e rho_a where dmat is a density matrix.
#include <implicit.h>
#include <inforb.h>
#include <maxorb.h>
#include <shells.h>
      DIMENSION DMAT(NBAST,NBAST), GAO(NVCLEN,NBAST,*)
      DIMENSION NBLOCKS(NSYM), IBLOCKS(2,LDAIB,NSYM), RHO(NVCLEN)
      DIMENSION TMP(NVCLEN,NBAST)
#include <priunit.h>
c
      CALL zeroorbs(TMP,NBLOCKS,IBLOCKS,LDAIB,NVCLEN)
c      CALL DZERO(TMP,NVCLEN*NBAST)
      DO ISYM=1,NSYM
         DO IBL=1, NBLOCKS(ISYM)
         DO IDX=IBLOCKS(1,IBL,ISYM), IBLOCKS(2,IBL,ISYM)
            DO JBL=1, NBLOCKS(ISYM)
            JTOP = MIN(IBLOCKS(2,JBL,ISYM),IDX-1)
            DO JDX=IBLOCKS(1,JBL,ISYM), JTOP
               DO K = 1, NVCLEN
                  TMP(K,JDX) = TMP(K,JDX) + GAO(K,IDX,1)*DMAT(IDX,JDX)
               END DO
            END DO
            END DO
            DO K = 1, NVCLEN
               TMP(K,IDX) = TMP(K,IDX) + GAO(K,IDX,1)*DMAT(IDX,IDX)*0.5
            END DO
         END DO
         END DO
      END DO
      call dzero(RHO, NVCLEN)
      DO ISYM=1,NSYM
         DO IBL=1, NBLOCKS(ISYM)
         DO IDX=IBLOCKS(1,IBL,ISYM), IBLOCKS(2,IBL,ISYM)
            DO K = 1, NVCLEN
               RHO(K) = RHO(K) + GAO(K,IDX,1)*TMP(K,IDX)*2D0
            END DO
         END DO
         END DO
      END DO
      END
c      getexp_blocked_lda does almost the same thing as getrho
c     but it does not use the simplification that the DMAT is symmetric.
      SUBROUTINE getexp_blocked_lda(IDSYM,DMAT,GAO,NBLOCKS,IBLOCKS,
     &                              LDAIB,TMP,NVCLEN,RHO)
#include <implicit.h>
#include <inforb.h>
      DIMENSION DMAT(NBAST,NBAST), GAO(NVCLEN,NBAST,*)
      DIMENSION NBLOCKS(NSYM), IBLOCKS(2,LDAIB,NSYM), RHO(NVCLEN)
      DIMENSION TMP(NVCLEN,NBAST)
#include <priunit.h>
c
      CALL zeroorbs(TMP,NBLOCKS,IBLOCKS,LDAIB,NVCLEN)
c     CALL DZERO(TMP,NVCLEN*NBAST)
      DO ISYM = 1, NSYM
      DO IBL=1, NBLOCKS(ISYM)
      DO IDX=IBLOCKS(1,IBL,ISYM), IBLOCKS(2,IBL,ISYM)
         JSYM = MULD2H(ISYM,IDSYM)
         DO JBL=1, NBLOCKS(JSYM)
         DO JDX=IBLOCKS(1,JBL,JSYM), IBLOCKS(2,JBL,JSYM)
            DO K = 1, NVCLEN
               TMP(K,JDX) = TMP(K,JDX) + GAO(K,IDX,1)*DMAT(IDX,JDX)
            END DO
         END DO
         END DO
      END DO
      END DO
      END DO
      DO K = 1, NVCLEN
         RHO(K) = 0D0
      END DO
      DO ISYM = 1, NSYM
      DO IBL=1, NBLOCKS(ISYM)
      DO IDX=IBLOCKS(1,IBL,ISYM), IBLOCKS(2,IBL,ISYM)
         DO K = 1, NVCLEN
            RHO(K) = RHO(K) + GAO(K,IDX,1)*TMP(K,IDX)
         END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE getrho_blocked_gga(DMAT,GAO,NBLOCKS,IBLOCKS,LDAIB,
     &                              TMP,NVCLEN,RHO,GRAD)
c     computes  <o|dmat|o'>
c     i.e rho_a where dmat is a density matrix (it can be a total
c     density end then one will get total density, or it can be an
c     alpha/beta density.
c     assert(NTYPSO>=NRHO)
#include <implicit.h>
#include <inforb.h>
      DIMENSION DMAT(NBAST,NBAST), GAO(NVCLEN,NBAST,*)
      DIMENSION NBLOCKS(NSYM),IBLOCKS(2,LDAIB,NSYM)
      DIMENSION RHO(NVCLEN), GRAD(3,NVCLEN),TMP(NVCLEN,NBAST)
#include <priunit.h>
c
      CALL zeroorbs(TMP,NBLOCKS,IBLOCKS,LDAIB,NVCLEN)
c     CALL DZERO(TMP,NVCLEN*NBAST)
      DO ISYM = 1, NSYM
         DO IBL=1, NBLOCKS(ISYM)
            ISTART=IBLOCKS(1,IBL,ISYM)
            ILEN=IBLOCKS(2,IBL,ISYM)-ISTART+1
            DO JBL=1, NBLOCKS(ISYM)
               JSTART=IBLOCKS(1,JBL,ISYM)
               JLEN=IBLOCKS(2,JBL,ISYM)-JSTART+1
               call dgemm('N','N',NVCLEN,JLEN,ILEN,1.0d0,
     &                    GAO(1,ISTART,1),NVCLEN,
     &                    DMAT(ISTART,JSTART),NBAST,1.0d0,
     &                    TMP(1,JSTART),NVCLEN)
            END DO
         END DO
      END DO
      call dzero(RHO, NVCLEN)
      call dzero(GRAD, 3*NVCLEN)
      DO ISYM = 1, NSYM
      DO IBL=1, NBLOCKS(ISYM)
      DO IDX=IBLOCKS(1,IBL,ISYM), IBLOCKS(2,IBL,ISYM)
         DO K = 1, NVCLEN
            RHO(K)    = RHO(K)    + GAO(K,IDX,1)*TMP(K,IDX)
            GRAD(1,K) = GRAD(1,K) + 2*GAO(K,IDX,2)*TMP(K,IDX)
            GRAD(2,K) = GRAD(2,K) + 2*GAO(K,IDX,3)*TMP(K,IDX)
            GRAD(3,K) = GRAD(3,K) + 2*GAO(K,IDX,4)*TMP(K,IDX)
         END DO
      END DO
      END DO
      END DO
      END
c
      SUBROUTINE getexp_blocked_gga(IDSYM,DMAT,GAO,NBLOCKS,IBLOCKS,
     &                              LDAIB,TMP,NVCLEN,GRAD)
c     very similar to getrho_blocked_gga but the computed values
c     are symmetrized i.e it computes <o|dmat|o'> + <o'|dmat|o>
c     GRAD(1,:) contains rho.
c     GRAD(2:4,:) contains gradient components.
c     assert(NTYPSO>=NRHO)
#include <implicit.h>
#include <inforb.h>
      DIMENSION DMAT(NBAST,NBAST), GAO(NVCLEN,NBAST,*)
      DIMENSION NBLOCKS(NSYM),IBLOCKS(2,LDAIB,NSYM), GRAD(4,NVCLEN)
      DIMENSION TMP(NVCLEN,NBAST)
#include <priunit.h>
c
      CALL zeroorbs(TMP,NBLOCKS,IBLOCKS,LDAIB,NVCLEN)
c      CALL DZERO(TMP,NVCLEN*NBAST)
      DO ISYM = 1, NSYM
      DO IBL=1, NBLOCKS(ISYM)
         JSYM = MULD2H(ISYM,IDSYM)
         DO JBL=1, NBLOCKS(JSYM)
            DO IDX=IBLOCKS(1,IBL,ISYM), IBLOCKS(2,IBL,ISYM)
               DO JDX=IBLOCKS(1,JBL,JSYM), IBLOCKS(2,JBL,JSYM)
                  DIJSYM = DMAT(IDX,JDX) + DMAT(JDX,IDX)
                  DO K = 1, NVCLEN
                     TMP(K,JDX) = TMP(K,JDX) + GAO(K,IDX,1)*DIJSYM
                  END DO
               END DO
            END DO
         END DO
      END DO
      END DO
      call dzero(GRAD, 4*NVCLEN)
      DO ISYM = 1, NSYM
      DO IBL=1, NBLOCKS(ISYM)
      DO IDX=IBLOCKS(1,IBL,ISYM), IBLOCKS(2,IBL,ISYM)
         DO K = 1, NVCLEN
            GRAD(1,K) = GRAD(1,K) + GAO(K,IDX,1)*TMP(K,IDX)*0.5
            GRAD(2,K) = GRAD(2,K) + GAO(K,IDX,2)*TMP(K,IDX)
            GRAD(3,K) = GRAD(3,K) + GAO(K,IDX,3)*TMP(K,IDX)
            GRAD(4,K) = GRAD(4,K) + GAO(K,IDX,4)*TMP(K,IDX)
         END DO
      END DO
      END DO
      END DO
      END
c
      SUBROUTINE GETBLOCKS(CENTER,CELLSZ,RSHEL,NBLCNT,IBLCKS)
c     get blocks of active SHELLS in cube of CELLSZ size centered at
C     CENTER.
c
c     RSHEL2 - precomputed shell extents (squared).
C     NBLCNT (output) - number of active blocks
c     IBLCKS (output) - pairs of (startindex, stopindex)
c     
c     symmetry handling takes some effort.
c
#include <implicit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <aovec.h>
#include <maxorb.h>
#include <primit.h>
#include <shells.h>
#include <dftcom.h>
#include <inforb.h>
#include <maxaqn.h>
#include <symmet.h>
      DIMENSION CENTER(3), RSHEL(KMAX), IBLCKS(2,NBAST)
      DIMENSION  ISSS(KMAX)
#include <ibtfun.h>

!#define NO_SHELL_SCREENING
#ifdef NO_SHELL_SCREENING
      NBLCNT = 1
      IBLCKS(1,1) = 1
      IBLCKS(2,1) = KMAX
      RETURN
#endif

      ISHLEN = 0
      NBLCNT = 0
      IPREV = -1111
      CELLDG = CELLSZ*0.5D0*SQRT(3D0)
      DO ISHELA=1,KMAX
         ICENT = NCENT(ISHELA)
         MULCNT = ISTBNU(ICENT)
c        
c        try different symmetry equivalent atoms:
c        
         DO ISYMOP = 0, MAXOPR         
            IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
               CORX = PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
               CORY = PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
               CORZ = PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
               PX = center(1)-CORX
               PY = center(2)-CORY
               PZ = center(3)-CORZ
               DST = SQRT(PX*PX + PY*PY + PZ*PZ)
c              RSHEL2(ISHELA)
               IF(DST.LE.RSHEL(ISHELA)+CELLDG) THEN
c                 accepted...
                  IF(ISHELA.NE.IPREV+1) THEN
                     NBLCNT = NBLCNT + 1
                     IBLCKS(1,NBLCNT) = ISHELA
                  END IF
                  IPREV = ISHELA
                  IBLCKS(2,NBLCNT) = ISHELA
                  ISSS(ISHELA) = ISYMOP
c
c                 since this shell has been accepted, there is
c                 no reason to try other symmetry-dependent
c                 atoms.
c         print "('shell',2I3,' at ',3F6.2,' rad. ',F6.2,' accepted')",
c     &        ISHELA,ICENT,CORD(1,ICENT),CORD(2,ICENT),CORD(3,ICENT),
c     &        SQRT(RSHEL2(ISHELA))
                  GO TO 10
               END IF
            END IF
         END DO
c         print "('shell',2I3,' at ',3F6.2,' rad. ',F6.2,' rejected')",
c     &        ISHELA,ICENT,CORD(1,ICENT),CORD(2,ICENT),CORD(3,ICENT),
c     &        RSHEL(ISHELA)
 10   CONTINUE
      END DO   ! ISHELA=1,KMAX
c      print "('cell at:',3F15.5,' blocks:',I3)", CENTER,NBLCNT
c      print "(8('[',2I3,']'))", ((IBLCKS(I,J),I=1,2),J=1,NBLCNT)
      END
c
      SUBROUTINE GTEXTS(RSHEL)
c     get radii of all shells as defined by specified threshold DFTHRI.
#include <implicit.h>
#include <aovec.h>
#include <maxorb.h>
#include <primit.h>
#include <shells.h>
#include <priunit.h>
#include <dftcom.h>
cdebug below
#include <mxcent.h>
#include <nuclei.h>
      DIMENSION FACL(10)
      DATA FACL /1D0,1.3333D0,1.6D0,1.83D0,2.03D0,
     &     2.22D0,2.39D0,2.55D0,2.70D0,2.84D0/
      DIMENSION RSHEL(KMAX)
      THLOG = LOG(DFTHRI)
      DO ISHELA=1,KMAX
         JSTA = JSTRT(ISHELA)
         RSHEL(ISHELA) =  0D0
         NUMCFA = NUMCF(ISHELA)
         DO IAOS = JSTA+1, JSTA + NUCO(ISHELA)
            IF(ABS(PRICCF(IAOS,NUMCFA)).GT.0D0) THEN
               R2 = (LOG(ABS(PRICCF(IAOS,NUMCFA)))-THLOG)/PRIEXP(IAOS)
               IF(RSHEL(ISHELA).LT.R2) RSHEL(ISHELA) = R2
            END IF
         END DO
      END DO
      DO ISHELA=1,KMAX
         RSHEL(ISHELA) = SQRT(RSHEL(ISHELA))*FACL(NHKT(ISHELA))
      END DO
      END
c
c     transform shell block indices to orbital block indices.
c     IORIDX contains preprocessed information about
c     where given shell begins and ends in given symmetry.

      SUBROUTINE SHLTOORB(NSHLBL,shlblock,NORBBL,orbblock,IORIDX)
#include <implicit.h>
#include <inforb.h>
      INTEGER SHLBLOCK(2,NSHLBL), ORBBLOCK(2,NSHLBL,NSYM)
      INTEGER IORIDX(2,KMAX,NSYM), NORBBL(NSYM)
#include <maxorb.h>
#include <shells.h>
#include <mxcent.h>
#include <nuclei.h>
      ILO = 0
      IHI = 0
      DO ISYM = 1, NSYM
         NORBBL(ISYM) = 0
         DO I = 1, NSHLBL
            DO ISHELL = SHLBLOCK(1,I), SHLBLOCK(2,I)
               ILO = IORIDX(1,ISHELL,ISYM)
               IF(ILO.NE.0) GO TO 10
            END DO
 10         continue
            DO ISHELL = SHLBLOCK(2,I), SHLBLOCK(1,I),-1
               IHI = IORIDX(2,ISHELL,ISYM)
               IF(IHI.GT.0) GO TO 20
            END DO
 20         continue
#if 0
            print '(I2,":shell ",2I4," translated to orbital ",2I4)',
     &           ISYM, SHLBLOCK(1,I), SHLBLOCK(2,I),ILO,IHI
#endif
            IF(ILO.LE.IHI) THEN
               NORBBL(ISYM) = NORBBL(ISYM) + 1
               INI = NORBBL(ISYM)
               ORBBLOCK(1,INI,ISYM) = ILO
               ORBBLOCK(2,INI,ISYM) = IHI
            END IF
         END DO
      END DO
      END
