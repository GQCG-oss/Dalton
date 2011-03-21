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
C  FILE: abacus/abaltra.F
C  PRINCIPAL AUTHOR: Hans Joergen Aa. Jensen 1992, 1995
C  PURPOSE: MO transformation of London two-electron integrals
C           for H2LAC(u,v,x,y) and FQL(p,u)
C
#ifdef UNDEF
===========================================================================
/* Comdeck log */
950102-hjaaj TRLTDS: improved vectorization (now not NASHA in inner loop)
941215-hjaaj removed TRFDIS (routine was duplicated in abadtra.F)
921117:to do: skip pack/unpack if NSYM .eq. 1
              option for only FQL if open-shell Hartree-Fock ?
              probably not worth the effort
===========================================================================
#endif
C  /* Deck lontra */
      SUBROUTINE LONTRA(IBXYZ,IBOPSY,H2LAC,FQL,CMO,PVSQ,WRK,LWRK,IPRINT)
C
C Copyright (c) 1992, 1995 Hans Joergen Aa. Jensen
C
C Purpose: MO transformation of London two-electron integrals
C
C Input:
C  IBXYZ = 1,2,3 for Bx, By, Bz component, respectively.
C  IBOPSY : symmetry of B operator component
C  CMO(NCMOT) : MO coefficients
C  PVSQ(NASHT,NASHT,NASHT,NASHT) : full 2-el. density matrix
C  IPRINT : print level
C Output:
C  H2LAC(NASHT,NASHT,NASHT,NASHT) : active 2-el. London int.s
C  FQL(NORBT,NASHT) : FQ matrix for London integrals
C Scratch:
C  WRK(LWRK)
C
#include <implicit.h>
#include <iratdef.h>
#include <priunit.h>
      DIMENSION H2LAC(*), FQL(*), CMO(*), PVSQ(N2ASHX,NASHT,NASHT)
      DIMENSION WRK(LWRK)
C
C Used from common blocks:
C  INFORB : NASHT,NBAST,N2BASX,N2ASHX
C  INFIND : NSM()
C  INFTRL : NDAC(),MXNDAO
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include "inftrl.h"
C
      CALL QENTER('LONTRA')
      if (iprint .gt. 2) then
         write (lupri,*) 'OUTPUT from LONTRA, iprint =',iprint
         write (lupri,*) 'ibxyz,ibopsy',ibxyz,ibopsy
         write (lupri,*) 'lwrk =',lwrk
         CALL FLSHFO(lupri)
      end if
      KFRSAV = 1
      LFRSAV = LWRK
      KFREE  = KFRSAV
      LFREE  = LFRSAV
      CALL MEMGET('INTE',KINDAO,N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDAC,N2ASHX,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDRS,64,    WRK,KFREE,LFREE)
C
      CALL TRLSET(WRK(KINDAO),WRK(KINDAC),WRK(KINDRS),IPRINT)
      if (iprint .gt. 2) then
         write (lupri,*) 'finished TRLSET'
         CALL FLSHFO(lupri)
      end if
C
C     Pack PVSQ
C
      CALL MEMGET('REAL',KPVFPK,N2ASHX,WRK,KFREE,LFREE)
      if (iprint .gt. 25) then
         write (lupri,*) 'Packing PV matrix'
      end if
      DO 280 KY = 1,NASHT
         KYSYM = NSM(KY)
         DO 260 KX = 1,NASHT
            KXSYM = NSM(KX)
            KXYSYM = MULD2H(KXSYM,KYSYM)
C           ... KUVSYM = KXYSYM for PVSQ(KU,KV,KX,KY) .ne. 0
            CALL TRLPAK(PVSQ(1,KX,KY),WRK(KPVFPK),NASH,IASH,NASHT,
     *                  KXYSYM,1)
            if (iprint .gt. 30) then
               write (lupri,*)
     *         'PV(u,v) square for KX,KY,KXYSYM',KX,KY,KXYSYM
               CALL OUTPUT(PVSQ(1,KX,KY),1,NASHT,1,NASHT,
     &            NASHT,NASHT,1,LUPRI)
               write (lupri,*)
     *         'PV(u,v) packed NDAC(KXYSYM)=',NDAC(KXYSYM)
               CALL wrtmat(WRK(KPVFPK),1,NDAC(KXYSYM),1,NDAC(KXYSYM),0)
            end if
            CALL DCOPY(NDAC(KXYSYM),WRK(KPVFPK),1,PVSQ(1,KX,KY),1)
  260    CONTINUE
  280 CONTINUE
      CALL MEMREL('LONTRA.PVPAK',WRK,KFRSAV,KPVFPK,KFREE,LFREE)
      if (iprint .gt. 2) then
         write (lupri,*) 'finished pack PVSQ'
         CALL FLSHFO(LUPRI)
      end if
C
      CALL MEMGET('REAL',KFQAOA,NBAST*NASHT,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDDI,N2BASX,WRK,KFREE,LFREE)
C     need LW1 in TRLDRV:
      LBINT = 600
C     record: BUF(600),IBUF(600),LENGTH
      LBUF = LBINT*(5*IRAT+1) + 1
      LW1  = MAX(IRAT*NBAST*NASHT,IRAT*2*MXNDAC,LBUF)
      LW1  = (LW1-1)/IRAT + 1
C
      MXDIS = (LFREE-LW1-10) / MXNDAO
      MXDIS = MIN(MXDIS,NNBASX)
      CALL MEMGET('REAL',KH2LBU,MXDIS*MXNDAO,WRK,KFREE,LFREE)
      CALL DZERO(WRK(KFQAOA),NBAST*NASHT)
      CALL TRLDRV(IBXYZ,IBOPSY,H2LAC,WRK(KFQAOA),CMO,PVSQ,
     *            MXDIS,WRK(KH2LBU),
     *            WRK(KINDAO),WRK(KINDAC),WRK(KINDRS),WRK(KINDDI),
     *            WRK,KFREE,LFREE,IPRINT)
C
      if (iprint .gt. 4) then
         write (lupri,*) 'finished TRLDRV'
         write (lupri,*) 'FQLAO(NBAST,NASHT) matrix:'
         CALL OUTPUT(WRK(KFQAOA),1,NBAST,1,NASHT,NBAST,NASHT,1,LUPRI)
         CALL FLSHFO(LUPRI)
      end if
      CALL MEMREL('LONTRA.TRLDRV',WRK,KFRSAV,KINDDI,KFREE,LFREE)
C
C     Final transformation of first index in FQAOAC from
C     AO basis to MO basis
C
      CALL TRLFQ(IBOPSY,FQL,WRK(KFQAOA),CMO)
      if (iprint .gt. 4) then
         write (lupri,*) 'finished TRLFQ'
         write (lupri,*) 'FQL(NORBT,NASHT) matrix:'
         CALL OUTPUT(FQL,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
         CALL FLSHFO(LUPRI)
      end if
C
C     921005-hjaaj: Because CISIGD is not modified for
C        symmetry packed H2AC yet, we must unpack H2LAC
C        and complete the matrix using (uv/xy) = -(vu/yx)
C
      CALL TRLH2U(IBOPSY,H2LAC,WRK(KINDAC),WRK(KFREE))
      if (iprint .gt. 4) then
         write (lupri,*) 'finished TRLH2U'
         if (iprint .gt. 5) then
            write (lupri,*) 'H2LAC matrix'
            call output(h2lac,1,n2ashx,1,n2ashx,n2ashx,n2ashx,1,LUPRI)
         end if
         CALL FLSHFO(LUPRI)
      end if
C
C     Unpack PVSQ again
C
      CALL MEMGET('REAL',KPVFPK,N2ASHX,WRK,KFREE,LFREE)
      DO 880 KY = 1,NASHT
         KYSYM = NSM(KY)
         DO 860 KX = 1,NASHT
            KXSYM = NSM(KX)
            KXYSYM = MULD2H(KXSYM,KYSYM)
C           ... KUVSYM = KXYSYM for PVSQ(KU,KV,KX,KY) .ne. 0
            CALL DCOPY(NDAC(KXYSYM),PVSQ(1,KX,KY),1,WRK(KPVFPK),1)
            CALL TRLPAK(PVSQ(1,KX,KY),WRK(KPVFPK),NASH,IASH,NASHT,
     *                  KXYSYM,-1)
  860    CONTINUE
  880 CONTINUE
      if (iprint .gt. 2) then
         write (lupri,*) 'finished PVSQ unpack'
         CALL FLSHFO(LUPRI)
      end if
C     CALL MEMREL('LONTRA.PVUPK',WRK,KFRSAV,KPVFPK,KFREE,LFREE)
C
      CALL MEMREL('LONTRA',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('LONTRA')
      RETURN
      END
C  /* Deck trldrv */
      SUBROUTINE TRLDRV(IBXYZ,IBOPSY,H2LAC,FQAOAC,CMO,PVF, MXDIS,
     *                  H2LBUF,INDAO,INDAC,INDRS,INDDIS,
     *                  WRK,KFRSAV,LFRSAV,IPRINT)
C
C Written 2.Oct.1992 hjaaj
C
C Driver for transformation of two-electron London integrals (H2L)
C
#include <implicit.h>
#include <iratdef.h>
#include <priunit.h>
      DIMENSION H2LAC(*),FQAOAC(NBAST,NASHT),CMO(NCMOT),PVF(*)
      DIMENSION H2LBUF(MXNDAO,*),INDDIS(NBAST,NBAST),WRK(*)
      DIMENSION INDAO(NBAST,NBAST), INDAC(NASHT,NASHT), INDRS(8,8)
C
C Used from common blocks:
C   INFORB : NCMOT,NBAST,NASHT,...
C   INFIND : ISAO(*)
C   INFTRL : MXNDAO
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <inftap.h>
#include <infind.h>
#include "inftrl.h"
C
      PARAMETER (LBINT = 600)
C
      KFREE = KFRSAV
      LFREE = LFRSAV
C     record: BUF(600),IBUF(600),LENGTH
      LBUF = LBINT*(5*IRAT+1) + 1
      LW1  = MAX(IRAT*NBAST*NASHT,IRAT*2*MXNDAC,LBUF)
      LW1    = (LW1-1)/IRAT + 1
      CALL MEMGET('REAL',KW1,LW1,WRK,KFREE,LFREE)
      KBUFW1 = KW1
      KBUFW1 = KW1
      KIBUF  = KBUFW1 + LBINT*(IRAT + 1) + 1
      KUVRY  = KW1
      KVUYR  = KUVRY + MXNDAC
      KPVXR  = KW1
      KPVXRT = KVUYR
C
      if (iprint .gt. 20) then
         DO 10 ISYM = 1,NSYM
            JCMO1 = ICMO(ISYM) + 1
            NBASI = NBAS(ISYM)
            NORBI = NORB(ISYM)
            write (lupri,*) 'TRLDRV: MO coefficients symmetry',ISYM
            CALL OUTPUT(CMO(JCMO1),1,NBASI,1,NORBI,NBASI,NORBI,1,LUPRI)
   10    CONTINUE
         CALL FLSHFO(LUPRI)
      end if
C
C     Initialize for TRFDIS
C
      ICASE = -1
      LRL   = 0
      LSL   = 0
C
C     Begin passes over AO file
C
      NPASS = 0
  100 CONTINUE
C
C        Read as many /rs) AO distributions as fit in memory.
C        INDDIS(r,s) is index for these distributions.
C        H2LBUF(pq,inddis(r,s)) all (pq/ for given /rs)
C
         CALL TRFDIS(ICASE,LRL,LSL,INDRS,INDDIS,MXDIS,NDIS)
         IF (IPRINT .GE. 15) THEN
            WRITE (lupri,*) 'TRLDRV: ICASE,LRL,LSL =',ICASE,LRL,LSL
         END IF
         IF (ICASE .LT. 0) GO TO 9999
C        ... finished when ICASE .lt. 0 !
         NPASS = NPASS + 1
         if (iprint .gt. 4) then
            write (lupri,*) 'TRLDRV: LU2DER pass no.',NPASS
            write (lupri,*) 'number of distributions:',NDIS
            IF (IPRINT .GE. 15) THEN
               WRITE (lupri,*) 'INDRS  array:'
               CALL IWRTMA(INDRS,NSYM,NSYM,8,8)
               IF (IPRINT .GE. 20) THEN
                  WRITE (lupri,*) 'INDDIS array:'
                  CALL IWRTMA(INDDIS,NBAST,NBAST,NBAST,NBAST)
            END IF
            END IF
            CALL FLSHFO(LUPRI)
         end if
C
         CALL DZERO(H2LBUF,NDIS*MXNDAO)
         CALL TRLRAO(LU2DER,IBXYZ,WRK(KIBUF),H2LBUF,
     *               INDAO,INDDIS,WRK(KBUFW1))
         if (iprint .gt. 20) then
            write (lupri,*) 'TRLDRV: finished TRLRAO in pass no.',NPASS
            CALL FLSHFO(LUPRI)
         end if
C
C        Transform (pq/rs) to (uv/rs) for the /rs) distributions
C        in memory.  H2LBUF(*,rs) is overwritten with (uv/rs)
C        integrals.
C
         CALL TRLPQ(IBOPSY,H2LBUF,CMO,INDDIS,WRK(KW1),IPRINT)
C        CALL TRLPQ(IBOPSY,H2LBUF,CMO,INDDIS,H2LWRK,IPRINT)
         if (iprint .gt. 20) then
            write (LUPRI,*) 'TRLDRV: finished TRLPQ in pass no.',NPASS
            CALL FLSHFO(LUPRI)
         end if
C
C        Use (uv/rs) for FQL and H2LAC
C
         KRSYM = 0
         DO 800 KR = 1,NBAST
            IF (ISAO(KR) .NE. KRSYM) THEN
               KRSYM = ISAO(KR)
               NASHR = NASH(KRSYM)
               NBASR = NBAS(KRSYM)
               IASHR = IASH(KRSYM)
               IBASR = IBAS(KRSYM)
               ICMORA = ICMO(KRSYM) + NISH(KRSYM)*NBASR + 1
            END IF
            IF (NASHR .EQ. 0) GO TO 800
            DO 700 KSSYM = 1,NSYM
               NASHS = NASH(KSSYM)
               KXSYM = MULD2H(KSSYM,IBOPSY)
               NASHX = NASH(KXSYM)
               IF (NASHX .EQ. 0 .AND. NASHS .EQ. 0) GO TO 700
C              Neither contribution for FQL or for H2LAC
               NBASS = NBAS(KSSYM)
               IBASS = IBAS(KSSYM)
               KS1 = IBASS + 1
               KSL = IBASS + NBASS
               IDOFQ = 0
               IDOH2 = 0
               DO 310 KS = KS1,KSL
                  IF (INDDIS(KR,KS) .GT. 0) THEN
                     IDOH2 = 1
                     IDOFQ = 1
                  ELSE IF (INDDIS(KR,KS) .NE. 0) THEN
                     IDOFQ = 1
                  END IF
  310          CONTINUE
               IF (IDOH2 .EQ. 0 .AND. IDOFQ .EQ. 0) GO TO 700
C
               KRSSYM = MULD2H(KRSYM,KSSYM)
               KUVSYM = MULD2H(KRSSYM,IBOPSY)
C   also :     KUVSYM = MULD2H(KXSYM,KRSYM)
               NDACUV = NDAC(KUVSYM)
               IF (NDACUV .EQ. 0) GO TO 700
               if (iprint .gt. 10) then
                  write (lupri,*) 'KR,KRSYM,KSSYM,KRSSYM',
     &               KR,KRSYM,KSSYM,KRSSYM
                  write (lupri,*) 'KUVSYM,NDACUV',KUVSYM,NDACUV
                  write (lupri,*) 'IDOH2, IDOFQ ',IDOH2, IDOFQ
               end if
C              ... No act-act elements in this symmetry block
C
C              First FQL contributions
C              NOTE: r and s are switched in this section,
C                    distributions are numbered (uv / sr)
C
C  INDSR .gt. 0:
C    FQAOAC(s,x) += sum(uv) (uv/sr) * PV(uv,rx)
C                           H2LBUF(uv,indsr) * PVXR(uv)
C  INDSR .lt. 0
C    FQAOAC(s,x) += sum(uv)  (uv/sr) * PV(uv,xr)
C                   sum(uv) -(vu/rs) * PV(vu,rx)
C                   sum(uv)  (vu/rs) * (-PV(vu,rx))
C                           H2LBUF(vu,-indsr) * PVXRT(vu)
C
C
C
               IF (NASHX .EQ. 0 .OR. IDOFQ .EQ. 0) GO TO 600
               IASHX = IASH(KXSYM)
               KS1 = IBASS + 1
               KSL = IBASS + NBASS
               DO 480 KX = IASHX+1,IASHX+NASHX
                  CALL TRLPV(KX,KR,KRSYM,CMO(ICMORA),NBASR,NASHR,
     *                       KUVSYM,PVF,WRK(KPVXR))
                  DO 460 KS = KS1,KSL
                     INDSR = INDDIS(KS,KR)
                     IF (INDSR .GT. 0) THEN
                        FQAOAC(KS,KX) = FQAOAC(KS,KX) +
     *                     DDOT(NDACUV,H2LBUF(1, INDSR),1,WRK(KPVXR),1)
                     ELSE IF (INDSR .LT. 0) THEN
                        FQAOAC(KS,KX) = FQAOAC(KS,KX) +
     *                     DDOT(NDACUV,H2LBUF(1,-INDSR),1,WRK(KPVXRT),1)
                     END IF
  460             CONTINUE
                  if (iprint .gt. 25) then
                     write (lupri,*) 'PVXR matrices for KR,KX=',KR,KX
                     CALL OUTPUT(WRK(KPVXR),1,NDACUV,1,2,MXNDAC,2,1,
     &                           LUPRI)
                     write (lupri,*) 'FQAOAC for KSSYM,KX=',KSSYM,KX
                     CALL OUTPUT(FQAOAC(1,1),KS1,KSL,KX,KX,
     &                           NBAST,NASHT,1,LUPRI)
                  end if
  480          CONTINUE
C
C
C              Second and last : H2LAC contributions
C
  600          IF (NASHS .EQ. 0 .OR. IDOH2 .EQ. 0) GO TO 700
               NBASS  = NBAS(KSSYM)
               ICMOSA = ICMO(KSSYM) + NISH(KSSYM)*NBASS + 1
               CALL TRLAH2(KR,KRSYM,KSSYM,KUVSYM,
     *              CMO(ICMORA),CMO(ICMOSA),NBASR,NASHR,NBASS,NASHS,
     *              INDDIS,H2LBUF,H2LAC,WRK(KUVRY),WRK(KVUYR),IPRINT)
  700       CONTINUE
  800    CONTINUE
         if (iprint .gt. 20) then
            write (lupri,*)
     &       'TRLDRV: finished transformation in pass no.',NPASS
            CALL FLSHFO(LUPRI)
         end if
C
      GO TO 100
C
 9999 CONTINUE
      IF (IPRINT .GE. 2) THEN
         WRITE (LUPRI,'(/A,I5,A)')
     &    'TRLDRV: finished transformation in',NPASS,
     &    ' passes over London AO integral file.'
         CALL FLSHFO(LUPRI)
      END IF
      RETURN
      END
C  /* Deck trlset */
      SUBROUTINE TRLSET(INDAO,INDAC,INDRS,IPRINT)
C
C Copyright 25-Sep-1992 Hans Joergen Aa. Jensen
C
C Set index inforamtion INDAO and INDAC and
C set index information for matrices in /INFTRL/.
C Set INDRS(rsym,ssym) = 0 for sym.dist. desired
C                      = 654321 for sym.dist. not needed
C
C IH1XX(IBSYM,IABSYM) is off-set for block (IASYM,IBSYM) in H1XX
C                     of symmetry IABSYM (IASYM = MULD2H(IBSYM,IABSYM))
C
#include <implicit.h>
      DIMENSION INDAO(NBAST,NBAST), INDAC(NASHT,NASHT), INDRS(8,8)
C
C Used from common blocks
C   INFORB: NSYM, NASH(8), NBAS(8), MULD2H(8,8)
C   INFTRL: IH1AO(8,8),IH1AC(8,8),NDAO(8), NDAC(8)
C
#include <priunit.h>
#include <inforb.h>
#include "inftrl.h"
C
      MXNDAO = 0
      MXNDAC = 0
      DO 500 IABSYM = 1,NSYM
         IH1AO(1,IABSYM) = 0
         IH1AC(1,IABSYM) = 0
         JNDAO = 0
         JNDAC = 0
         DO 400 IBSYM = 1,NSYM
            IASYM = MULD2H(IABSYM,IBSYM)
C           ... 1. IH1AO and IH1AC
            IF (IBSYM .LT. NSYM) THEN
               IH1AO(IBSYM+1,IABSYM) = IH1AO(IBSYM,IABSYM) +
     *            NBAS(IASYM)*NBAS(IBSYM)
               IH1AC(IBSYM+1,IABSYM) = IH1AC(IBSYM,IABSYM) +
     *            NASH(IASYM)*NASH(IBSYM)
            END IF
C           ... 2. INDAO
            DO 140 JB = IBAS(IBSYM)+1,IBAS(IBSYM)+NBAS(IBSYM)
               DO 120 JA = IBAS(IASYM)+1,IBAS(IASYM)+NBAS(IASYM)
                  JNDAO = JNDAO + 1
                  INDAO(JA,JB) = JNDAO
  120          CONTINUE
  140       CONTINUE
            NDAO(IABSYM) = JNDAO
            MXNDAO = MAX(MXNDAO,JNDAO)
C           ... 3. INDAC
            DO 240 JB = IASH(IBSYM)+1,IASH(IBSYM)+NASH(IBSYM)
               DO 220 JA = IASH(IASYM)+1,IASH(IASYM)+NASH(IASYM)
                  JNDAC = JNDAC + 1
                  INDAC(JA,JB) = JNDAC
  220          CONTINUE
  240       CONTINUE
            NDAC(IABSYM) = JNDAC
            MXNDAC = MAX(MXNDAC,JNDAC)
  400    CONTINUE
  500 CONTINUE
C
      DO 840 IBSYM = 1,NSYM
         DO 820 IASYM = 1,NSYM
            IF (NASH(IASYM) .GT. 0 .OR. NASH(IBSYM) .GT. 0) THEN
               INDRS(IASYM,IBSYM) = 0
            ELSE
               INDRS(IASYM,IBSYM) = 654321
            END IF
  820    CONTINUE
  840 CONTINUE
C
      IF (IPRINT .GT. 4) THEN
         write (lupri,*) 'output from TRLSET'
         WRITE (lupri,*) 'NASH(*)',(NASH(I),I=1,NSYM)
         WRITE (lupri,*) 'NBAS(*)',(NBAS(I),I=1,NSYM)
         write (lupri,*) 'MXNDAC,MXNDAO',MXNDAC,MXNDAO
         write (lupri,*) 'NDAC(*)',(NDAC(I),I=1,NSYM)
         write (lupri,*) 'NDAO(*)',(NDAO(I),I=1,NSYM)
         IF (IPRINT .GT. 10) THEN
            write (lupri,*) 'INDAC matrix'
            CALL IWRTMA(INDAC,NASHT,NASHT,NASHT,NASHT)
            write (lupri,*) 'INDAO matrix'
            CALL IWRTMA(INDAO,NBAST,NBAST,NBAST,NBAST)
         END IF
         write (lupri,*) 'IH1AC matrix'
         CALL IWRTMA(IH1AC,NSYM,NSYM,8,8)
         write (lupri,*) 'IH1AO matrix'
         CALL IWRTMA(IH1AO,NSYM,NSYM,8,8)
      END IF
      RETURN
      END
C  /* Deck trlpq */
      SUBROUTINE TRLPQ(IBOPSY,H2LBUF,CMO,INDDIS,H2LWRK,IPRINT)
C
C 6-Oct-1992 Hans Joergen Aa. Jensen
C
C Transform (pq/rs) to (uv/rs) for the available
C /rs) distributions according to INDDIS(r,s).
C
C Input : H2LBUF(pq,indrs) = (pq / rs)
C Output: H2LBUF(uv,indrs) = (uv / rs)
C
#include <implicit.h>
C
      DIMENSION H2LBUF(MXNDAO,*), CMO(NCMOT)
      DIMENSION INDDIS(NBAST,NBAST), H2LWRK(*)
C
C Used from common blocks:
C   INFORB : NBAST, NCMOT
C   INFIND : ISAO(*)
C   INFTRL : MXNDAO
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include "inftrl.h"
C
      DO 800 KR = 1,NBAST
         KRSYM = ISAO(KR)
         DO 700 KS = 1,NBAST
            INDRS = INDDIS(KR,KS)
            IF (INDRS .GT. 0) THEN
               KSSYM = ISAO(KS)
               IRSSYM = MULD2H(KRSYM,KSSYM)
               IPQSYM = MULD2H(IRSSYM,IBOPSY)
               IF (IPRINT .GT. 25) THEN
                  NDAOPQ = NDAO(IPQSYM)
                  WRITE (LUPRI,*) 'H2CDAO(pq,rs) for r,s=',KR,KS
                  WRITE (LUPRI,*) 'indrs,sym(pq), ndao(pqsym)',
     *               INDRS,IPQSYM,NDAOPQ
                  CALL WRTMAT(H2LBUF(1,INDRS),1,NDAOPQ,1,NDAOPQ,0)
               END IF
               CALL TRLTDS(IPQSYM,H2LBUF(1,INDRS),H2LBUF(1,INDRS),
     *                     CMO,H2LWRK)
C              Note: equivalence(H2CDAO,H2CDAC)
C              CALL TRLTDS(IABSYM,H2CDAO,H2CDAC,CMO,WRK)
               IF (IPRINT .GT. 25) THEN
                  NDACUV = NDAC(IPQSYM)
                  WRITE (LUPRI,*) 'H2CDAC(uv,rs) for r,s=',KR,KS
                  WRITE (LUPRI,*) 'sym(uv), ndac(uvsym)',IPQSYM,NDACUV
                  CALL WRTMAT(H2LBUF(1,INDRS),1,NDACUV,1,NDACUV,0)
               END IF
            END IF
  700    CONTINUE
  800 CONTINUE
C
      RETURN
      END
C  /* Deck trltds */
      SUBROUTINE TRLTDS(IABSYM,H2CDAO,H2CDAC,CMO,WRK)
C
C Copyright (c) 24-Sep-1992/2-Jan-1995 Hans Joergen Aa. Jensen
C
C Transform H2CDAO(nbast,nbast) of symmetry IABSYM
C        to H2CDAC(nasht,nasht) using CMO.
C H2CDAO(a,b) .ne. H2CDAO(b,a) is assumed.
C H2CDAO and H2CDAC are both symmetry packed.
C
C H2CDAO and H2CDAC may share same memory location.
C
C Input : H2CDAO(a,b) = (a b | c d)
C         CMO = m.o. coefficients
C         IABSYM = symmetry of H2CDAO, H2CDAC matrices
C Output: H2CDAC(u,v) = (u v | c d) =
C         sum(b) [ sum(a) [ CMO(a,u) H2CDAO(a,b) ] CMO(b,v) ]
C
#include <implicit.h>
      DIMENSION CMO(NCMOT), H2CDAO(*), H2CDAC(*), WRK(*)
C
C Used from common blocks:
C  INFORB : NSYM, MULD2H, NASH(*), NBAS(*), NISH(*), ICMO(*)
C  INFTRL : IH1AO(8,8), IH1AC(8,8)
C
#include <inforb.h>
#include "inftrl.h"
C
      DO 100 IBSYM = 1,NSYM
         IASYM = MULD2H(IBSYM,IABSYM)
         NASHB = NASH(IBSYM)
         NASHA = NASH(IASYM)
         IF ((NASHA.NE.0) .AND. (NASHB.NE.0)) THEN
            NBASB = NBAS(IBSYM)
            NBASA = NBAS(IASYM)
            JCMOB = ICMO(IBSYM) + NISH(IBSYM)*NBASB + 1
            JCMOA = ICMO(IASYM) + NISH(IASYM)*NBASA + 1
            JH2AO = IH1AO(IBSYM,IABSYM) + 1
            JH2AC = IH1AC(IBSYM,IABSYM) + 1
            CALL DGEMM('N','N',NBASA,NASHB,NBASB,1.D0,
     &                 H2CDAO(JH2AO),NBASA,
     &                 CMO(JCMOB),NBASB,0.D0,
     &                 WRK,NBASA)
            CALL DGEMM('T','N',NASHA,NASHB,NBASA,1.D0,
     &                 CMO(JCMOA),NBASA,
     &                 WRK,NBASA,0.D0,
     &                 H2CDAC(JH2AC),NASHA)
         END IF
  100 CONTINUE
C
C     End of TRLTDS
C
      RETURN
      END
C  /* Deck trlrao */
      SUBROUTINE TRLRAO(LULONA,IBXYZ,IBUF,H2LAO,INDAO,INDDIS,WRK)
C
C 25-Sep-1992 hjaaj / revised 6-Oct-1992
C
C Reads two-electron London integrals from file to H2LAO buffer.
C LULONA: File to read
C IBXYZ : =1,2,3 for Bx,By,Bz respectively.
C IBUF  : Points to start element in label buffer.
C H2LAO : H2LAO(r,s,pq), pq=(mpqoff+1:mpqoff+npq); dim (nbasr,nbass,npq)
C
#include <implicit.h>
#include <iratdef.h>
#include <priunit.h>
#include <mxcent.h>
      DIMENSION H2LAO(MXNDAO,*), IBUF(4,600), WRK(600*(IRAT + 1) + 1)
      DIMENSION INDAO(NBAST,NBAST), INDDIS(NBAST,NBAST)
C
C Used from common blocsk:
C   INFORB : IBAS(8)
C   INFIND : IROW(*),ISAO(*)
C   INFTRL : MXNDAO
C
#include <nuclei.h>
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include "inftrl.h"
C
C     Local variables:
C
      CHARACTER*8 KEY
      LOGICAL DOCOOR, FIRST
C
#include <ibtdef.h>
C
      CALL REWSPL(LULONA)
      KEY = 'AO2MGINT'
      CALL MOLLAB(KEY,LULONA,LUPRI)
#if defined (VAR_INT64)
      NIBUF = 1
      NBITS = 16
#else
      IF (NBASIS .LE. 255) THEN
         NIBUF = 1
         NBITS = 8
      ELSE
         call quit('TRLRAO only implemented for < 256 basis functions')
         NIBUF = 2
         NBITS = 16
      END IF
#endif
      LBUF = 600
      JBUF = 256
C
      LENINT = IRAT*LBUF + LBUF + 1
      KINT = 1
      KIINT = KINT + LBUF
C
      DOCOOR = .FALSE.
 150  CONTINUE
#if defined (VAR_SPLITFILES)         
         CALL READSI(LULONA,LENINT,WRK(KINT),JBUF)
#else
         CALL READI(LULONA,LENINT,WRK(KINT))
#endif
         CALL AOLAB4(WRK(KIINT),LBUF,NIBUF,NBITS,IBUF,LENGTH)
         IF (LENGTH .GT. 0) THEN
            DO 100 I = 1, LENGTH
               JS    = IBUF(4,I)
               IF (JS .EQ. 0) THEN
                  ICOOR = IBUF(3,I)
                  DOCOOR = ICOOR.EQ.IBXYZ
               ELSE IF (DOCOOR) THEN
                  JR = IBUF(3,I)
                  JP = IBUF(1,I)
                  JQ = IBUF(2,I)
C                 BUF(I) =  {P Q | R S} =  {R S | P Q}
C                        = -{Q P | S R} = -{S R | Q P}
                  JRS = INDDIS(JR,JS)
                  IF (JRS .GT. 0) THEN
                     H2LAO(INDAO(JP,JQ),JRS) = WRK(KINT + I - 1)
                     IF (JR .EQ. JS) H2LAO(INDAO(JQ,JP),JRS) = 
     &                    -WRK(KINT + I - 1)
                  ELSE IF (JRS .LT. 0) THEN
                     H2LAO(INDAO(JQ,JP),-JRS) = -WRK(KINT + I - 1)
C                    {Q P | S R} = -{P Q | R S}
                  END IF
                  JPQ = INDDIS(JP,JQ)
                  IF (JPQ .GT. 0) THEN
                     H2LAO(INDAO(JR,JS),JPQ) = WRK(KINT + I - 1)
                     IF (JP .EQ. JQ) H2LAO(INDAO(JS,JR),JPQ) = 
     &                    -WRK(KINT + I -1)
                  ELSE IF (JPQ .LT. 0) THEN
                     H2LAO(INDAO(JS,JR),-JPQ) = -WRK(KINT + I - 1)
C                    {S R | Q P} = -{P Q | R S}
                  END IF
               ENDIF
 100        CONTINUE
         ELSE IF (LENGTH .LT. 0 ) THEN
            GO TO 300
         END IF
         GO TO 150
C
 300  CONTINUE
      RETURN
      END
C  /* Deck trlah2 */
      SUBROUTINE TRLAH2(KR,KRSYM,KSSYM,KUVSYM,
     *                  CMOR,CMOS,NBASR,NASHR,NBASS,NASHS,
     *                  INDDIS,H2UVRS,H2LAC,H2UVRY,H2VUYR,IPRINT)
C
C 5-Oct-1992 Hans Joergen Aagaard Jensen
C 'AH2' for Add to H2LAC
C Add (uv / rs) contribution to (uv / xy) (x .ge. y)
C
#include <implicit.h>
      DIMENSION CMOR(NBASR,NASHR), CMOS(NBASS,NASHS)
      DIMENSION INDDIS(NBAST,NBAST)
      DIMENSION H2UVRS(MXNDAO,*), H2LAC(N2ASHX,NASHT,NASHT)
      DIMENSION H2UVRY(MXNDAC), H2VUYR(MXNDAC)
      PARAMETER ( DP5 = 0.5D0, D0 = 0.0D0 )
C
C Used from common blocks:
C   INFORB : ??
C   INFTRL : NDAC(8),MXNDAO,NXNDAC
C
#include <priunit.h>
#include <inforb.h>
#include "inftrl.h"
C
      IF (IPRINT .GT. 25) THEN
         write (lupri,*) 'output from trlah2, kr,kssym=',kr,kssym
      END IF
      JR     = KR - IBAS(KRSYM)
      NDACUV = NDAC(KUVSYM)
      IBASS  = IBAS(KSSYM)
      DO 400 JY = 1,NASHS
         KY = IASH(KSSYM) + JY
         CALL DZERO(H2UVRY,NDACUV)
         DO 200 JS = 1,NBASS
            KS = IBASS + JS
            INDRS = INDDIS(KR,KS)
            IF (INDRS .GT. 0) THEN
               CMOSY = CMOS(JS,JY)
               IF (CMOSY .NE. D0) THEN
                  IF (KS .EQ. KR) CMOSY = DP5 * CMOSY
C                 ... KR .eq. KS contribution is included both
C                     in 300 loop and 350 loop, therefore we must
C                     multiply by 0.5 (921019-hjaaj)
                  CALL DAXPY(NDACUV,CMOSY,H2UVRS(1,INDRS),1,H2UVRY,1)
               END IF
            END IF
  200    CONTINUE
         IF (IPRINT .GT. 25) THEN
            write (lupri,*) 'TRLAH2 H2UVRY for r,y=',KR,KY
            CALL OUTPUT(H2UVRY,1,NDACUV,1,1,NDACUV,1,1,LUPRI)
         END IF
C
         IF (KRSYM .EQ. KSSYM) THEN
            JX1 = JY
         ELSE
            JX1 = 1
         END IF
         DO 300 JX = JX1,NASHR
            KX = IASH(KRSYM) + JX
            CMORX = CMOR(JR,JX)
            IF (CMORX .NE. D0) THEN
               CALL DAXPY(NDACUV,CMORX,H2UVRY,1,H2LAC(1,KX,KY),1)
            END IF
  300    CONTINUE
C     Special case for KRSYM .eq. KSSYM
C     we must include JX .le. JY to get all contributions
C     in diagonal blocks.
         IF (KRSYM .EQ. KSSYM) THEN
            CALL TRLTRS(H2UVRY,H2VUYR,KUVSYM,NASH,IH1AC)
            DO 350 JX = 1,JY
               KX = IASH(KRSYM) + JX
               CMORX = -CMOR(JR,JX)
C              ... minus because H2VUYR = -{vu|yr}
               IF (CMORX .NE. D0) THEN
                  CALL DAXPY(NDACUV,CMORX,H2VUYR,1,H2LAC(1,KY,KX),1)
               END IF
  350       CONTINUE
         END IF
  400 CONTINUE
      RETURN
C     end of TRLAH2
      END
C  /* Deck trltrs */
      SUBROUTINE TRLTRS(H2UVRY,H2VUYR,IABSYM,NBLK,INDBLK)
C
C 8-Oct-1992 hjaaj
C   NBLK(8) will typically be NASH(8)
C   INDBLK(8,8) will typically be IH1AC(8,8)
C
#include <implicit.h>
      DIMENSION H2UVRY(*), H2VUYR(*), NBLK(8), INDBLK(8,8)
C
C Used from common blocks:
C   INFORB : NSYM, MULD2H()
C
#include <inforb.h>
C
      DO 700 IBSYM = 1,NSYM
         IASYM = MULD2H(IBSYM,IABSYM)
         JUVRY = INDBLK(IBSYM,IABSYM) + 1
         JVUYR = INDBLK(IASYM,IABSYM) + 1
         NBLKA = NBLK(IASYM)
         NBLKB = NBLK(IBSYM)
         CALL MTRSP(NBLKA,NBLKB,H2UVRY(JUVRY),NBLKA,
     *              H2VUYR(JVUYR),NBLKB)
  700 CONTINUE
      RETURN
      END
C  /* Deck trlh2u */
      SUBROUTINE TRLH2U(IBOPSY,H2LAC,INDAC,WRK)
C
C 6-Oct-1992 hjaaj
C
#include <implicit.h>
      DIMENSION H2LAC(N2ASHX,NASHT,NASHT)
      DIMENSION INDAC(NASHT,NASHT), WRK(*)
C
C Used from common blocks:
C  INFORB : NASHT,NASH(),IASH()
C  INFIND : NSM()
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
C
      PARAMETER (IWAY = -1, DM1 = -1.0D0)
C     ... IWAY=-1 to unpack
      DO 800 KY = 1,NASHT
         KYSYM = NSM(KY)
C        Unpack KX .ge. KY
         DO 600 KX = KY,NASHT
            KXSYM = NSM(KX)
            KXYSYM = MULD2H(KXSYM,KYSYM)
            KUVSYM = MULD2H(KXYSYM,IBOPSY)
            CALL TRLPAK(WRK,H2LAC(1,KX,KY),
     *                  NASH,IASH,NASHT,KUVSYM,IWAY)
C           CALL TRLPAK(H1SQ,H1PK,NBLK,IBLK,NBLKT,IH1SYM,IWAY)
            CALL DCOPY(N2ASHX,WRK,1,H2LAC(1,KX,KY),1)
C
C           Fill in KX .lt. KY using {uv/xy} = -{vu/yx}
C           HJMAERK: implement new option in GETIN2 for
C                    H2AC(N2ASHX,NNASHX) (921007)
C
            IF (KX .GT. KY) THEN
               CALL DSCAL(N2ASHX,DM1,WRK,1)
               CALL MTRSP(NASHT,NASHT,WRK,NASHT,H2LAC(1,KY,KX),NASHT)
            END IF
  600    CONTINUE
  800 CONTINUE
      RETURN
      END
C  /* Deck trlpv */
      SUBROUTINE TRLPV(KX,KS,KSSYM,CMOS,NBASS,NASHS,
     *                 KUVSYM,PVF,PVXS)
C
C 6-Oct-1992 Hans Joergen Aa. Jensen
C
C PVXS(uv,1) =  sum(y) PVF(uv,x,y) * CMO(s,y) =  [uv|xs]
C PVXS(uv,2) = -sum(y) PVF(uv,y,x) * CMO(s,y) = -[uv|sx] = -[vu|xs]
C
#include <implicit.h>
      DIMENSION CMOS(NBASS,NASHS), PVF(N2ASHX,NASHT,NASHT)
      DIMENSION PVXS(MXNDAC,2)
C
      PARAMETER (DM1 = -1.0D0)
C
C Used from common blocks:
C  INFORB : NASHT
C  INFTRL : NDAC(8), IH1AC(8,8)
C
#include <inforb.h>
#include "inftrl.h"
C
      NDACUV = NDAC(KUVSYM)
      CALL DZERO(PVXS,2*MXNDAC)
      IASHS  = IASH(KSSYM)
      JS     = KS - IORB(KSSYM)
C     Note: if PVF packed KX .ge. KY we need to check
C           and use PVF(1,KXY) transposed if KX .LT. KY
      DO 800 JY = 1,NASHS
         KY = JY + IASHS
         CALL DAXPY(NDACUV,CMOS(JS,JY),PVF(1,KX,KY),1,PVXS(1,1),1)
  800 CONTINUE
      CALL TRLTRS(PVXS(1,1),PVXS(1,2),KUVSYM,NASH,IH1AC)
      CALL DSCAL(NDACUV,DM1,PVXS(1,2),1)
      RETURN
      END
C  /* Deck trlfq */
      SUBROUTINE TRLFQ(IBOPSY,FQL,FQAOAC,CMO)
C
C 7-Oct-1992 hjaaj
C
C FQL(a,x) = sum(r) FQAOAC(r,x) * CMO(r,a)
C
#include <implicit.h>
      DIMENSION FQL(NORBT,NASHT), FQAOAC(NBAST,NASHT), CMO(NCMOT)
C
C Used from common blocks:
C  INFORB : NORBT,NASHT,NBAST,NCMOT,NSYM
C
#include <inforb.h>
C
      DO 800 KSYM = 1,NSYM
         NASHK  = NASH(KSYM)
      IF (NASHK.EQ.0) GO TO 800
         IASYM  = MULD2H(KSYM,IBOPSY)
         NORBA  = NORB(IASYM)
         NBASA  = NBAS(IASYM)
         IASHK1 = IASH(KSYM) + 1
         IORBA1 = IORB(IASYM) + 1
         IBASA1 = IBAS(IASYM) + 1
         JCMOA1 = ICMO(IASYM) + 1
         CALL DGEMM('T','N',NORBA,NASHK,NBASA,1.D0,
     &              CMO(JCMOA1),NBASA,
     &              FQAOAC(IBASA1,IASHK1),NBAST,0.D0,
     &              FQL(IORBA1,IASHK1),NORBT)
  800 CONTINUE
      RETURN
      END
C  /* Deck trlpak */
      SUBROUTINE TRLPAK(H1SQ,H1PK,NBLK,IBLK,NBLKT,IH1SYM,IWAY)
C
C Copyright 6-Oct-1992 Hans Joergen Aa. Jensen
C (based on OITPAK)
C
C NBLK,IBLK,NBLKT may e.g. be NORB(),IORB(),NORBT or NASH(),IASH(),NASHT
C
C IWAY .ge. 0 : Pack H1SQ in H1PK
C      else   : Unpack H1PK in H1SQ
C
#include <implicit.h>
      DIMENSION H1SQ(NBLKT,NBLKT), H1PK(*), NBLK(8), IBLK(8)
C
C Used from common blocks:
C   INFORB : NSYM,MULD2H
C
#include <inforb.h>
C
      N2BLKX = NBLKT*NBLKT
      IF (IWAY .LT. 0) CALL DZERO(H1SQ,N2BLKX)
      ISTH1P = 1
      DO 300 IBSYM = 1,NSYM
         IASYM  = MULD2H(IBSYM,IH1SYM)
         NBLKA  = NBLK(IASYM)
         NBLKB  = NBLK(IBSYM)
      IF (NBLKA .EQ. 0 .OR. NBLKB .EQ. 0) GO TO 300
         IAST   = IBLK(IASYM) + 1
         IBST   = IBLK(IBSYM) + 1
         IF (IWAY .GE. 0) THEN
            CALL MCOPY(NBLKA,NBLKB,H1SQ(IAST,IBST),NBLKT,
     *                 H1PK(ISTH1P),NBLKA)
         ELSE
            CALL MCOPY(NBLKA,NBLKB,H1PK(ISTH1P),NBLKA,
     *                 H1SQ(IAST,IBST),NBLKT)
         END IF
         ISTH1P = ISTH1P + NBLKA*NBLKB
  300 CONTINUE
C
C         MCOPY(nrowa,ncola,A,nrdima,B,nrdimb)
C
      RETURN
      END
