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
/* Comdeck log */
NOTE: H2DAO and H2DAC submatrices are transposed compared to abaltra
  because of possible triangular NNASHX packing.

TO DO: s/INFTRL/CBDTRA/ ????

from abaltra.u:
921117:to do: skip pack/unpack if NSYM .eq. 1
              option for only FQD if open-shell Hartree-Fock ?
              probably not worth the effort
===========================================================================
#endif
C  /* Deck dertra */
      SUBROUTINE DERTRA(IATOM,CMO,PV,H2DAC,FQD,WORK,KWORK,MWORK)
C
C     ***************************************************************
C     *****   Transformation of symmetry integral derivatives   *****
C     ***************************************************************
C
C Written 1992 by Hans Joergen Aa. Jensen.
C
C IPRTRD in common /CBDTRA/ controls printout,
C IPRTRD .le. 0, no printout
C IPRTRD .gt. 0, JTER is used to prevent repetition of timing, number
C                of integrals and other statistics.
C
C Input:
C  WORK, work array
C  KWORK, actual length of work array
C  MWORK, desired length of work array (NOT USED IN THIS VERSION)
C         (dependent on size of physical memory)
C
C
C     Purpose: transform nuclear derivative integrals to mo basis,
C              returned in H2DAC and FQD.
C
C     Called from nuclear Hessian calculation.
C
C
#include <implicit.h>
#include <maxash.h>
#include <maxorb.h>
#include <mxcent.h>
      DIMENSION CMO(*),PV(*),FQD(*),H2DAC(*),WORK(KWORK)
      LOGICAL OPNSTA
#include <iratdef.h>
C
#include <priunit.h>
#include <nodint.h>
#include <cbdtra.h>
#include <nuclei.h>
#include <inforb.h>
#include <infind.h>
C
      CHARACTER*1 OPRTYP(3)
      PARAMETER (LBINT = 600)
      SAVE JATOM
      DATA JATOM /0/
      DATA OPRTYP /'x', 'y', 'z'/
C
C     This version runs from a file which contains only one
C     symmetry-adapted integral derivative per label
C
      NUMSIM = 1
C
      IF (IPRTRD.GT.0) THEN
        IF (JATOM.NE.IATOM) THEN
          JTER = 1
        ELSE
          JTER = JTER + 1
        ENDIF
        IF (IPRTRD .GT. 10) JTER = 1
      ELSE
        JTER = 0
      ENDIF
C
C     Decide direction and symmetry type of differentiation.  This
C     is held in the LOGICAL arrays SRTCOO and SRTSYM in /CBDTRA/
C
      JCOOR = 0
      DO 10 KCOOR = 1,3
         IF (SRTCOO(KCOOR)) JCOOR = KCOOR
10    CONTINUE
      JCORSY = 0
      DO 20 KOPSYM = 1,NSYM
         IF (SRTSYM(KOPSYM)) JCORSY = KOPSYM
20    CONTINUE
      IF (JCOOR .EQ. 0 .OR. JCORSY .EQ. 0) THEN
         WRITE(LUPRI,'(/A,//,2(A,I4,/))')
     1      ' DERTRA: error in coordinate or symmetry specification',
     2      ' Direction: ', JCOOR,
     3      ' Symmetry:  ', JCORSY
         CALL QUIT(' Invalid specification in DERTRA')
      ENDIF
      IF (NASHT.EQ.0) THEN
       WRITE(LUPRI,'(/A)')' DERTRA: NASHT = 0, transformation skipped.'
       GO TO 9999
      ENDIF
C
C
      IF (IPRTRD .GT. 1) THEN
         WRITE(LUPRI,'(//,A,//,2A,/,A,I1,/)')
     1      ' Transformation of SO integral derivatives to MO basis',
     2      ' Coordinate direction:   ',OPRTYP(JCOOR),
     3      ' Symmetry of derivative: ',JCORSY
      ENDIF
      LWORK = KWORK
      IF((IPRTRD.GT.2) .AND. (JTER.EQ.1)) WRITE(LUPRI,903) LWORK
  903 FORMAT(/' (DERTRA) Working area in transformation section',I8)
C
C     ***** TYPE OF TRANSFORMATION:                           *****
C     *****         FIRST INDEX,ALL ORBITALS                  *****
C     *****         OTHER INDICES,ACTIVE ORBITALS ONLY        *****
C          (always first order)
C
C
C     The total number of integral derivatives for this case
C     is passed in /NODINT/
C
      LPQRS = NOINTS
      IF((IPRTRD.GT.2) .AND. (JTER.EQ.1)) WRITE(LUPRI,902) LPQRS
  902 FORMAT(' Maximum number of integrals =',I10)
C
C     max MXCORB orbitals (and max 255 because of index packing)
C
      MXBAST = MIN(255,MXCORB)
      IF (NORBT.GT.MXBAST .OR. NBAST.GT.MXBAST) THEN
         WRITE (LUPRI,9020) MXBAST, NORBT, NBAST
         CALL QUIT('DERTRA error: too many orbitals for this version')
      ENDIF
 9020 FORMAT(' (DERTRA) ERROR, DERTRA is limited to',I5,' orbitals,'
     *      /'          NORBT is',I5,' and NBAST is',I5)
C
      JATOM = IATOM
      CALL TRDNEW(JCOOR,JCORSY,H2DAC,FQD,CMO,PV,WORK,LWORK,IPRTRD)
C     CALL TRDNEW(JCOOR,JCORSY,H2DAC,FQD,CMO,PV,WRK,LWRK,IPRINT)
C
 9999 CALL FLSHFO(LUPRI)
      RETURN
C
      END
C  /* Deck trdnew */
      SUBROUTINE TRDNEW(JCOOR,JCORSY,H2DAC,FQD,CMO,PV,WRK,LWRK,IPRINT)
C
C 23-Sep-1992 hjaaj
C
C Input:
C  JCOOR = 1,2,3 for d/dx, d/dy, d/dz component, respectively.
C  JCORSY symmetry of JCOOR component
C  CMO(NCMOT) MO coefficients
C  PV(NNASHX,NNASHX) 2-el. density matrix
C Output:
#if defined (VAR_H2DPACK)
C  H2DAC(NNASHX,NNASHX) active 2-el. derivative int.s
#else
C  H2DAC(N2ASHX,N2ASHX) active 2-el. derivative int.s
#endif
C  FQD(NORBT,NASHT) FQ matrix for derivative integrals
C Scratch:
C  WRK(LWRK)
C
#include <implicit.h>
#include <iratdef.h>
C Used from common blocks:
C  INFORB : NASHT,NBAST,N2BASX,N2ASHX,NNASHX
C  INFIND : NSM(),IROW()
C  INFTRL : NDAC(),MXNDAO
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include "inftrl.h"
C
      DIMENSION H2DAC(*), FQD(*), CMO(*), PV(NNASHX,NNASHX)
      DIMENSION WRK(LWRK)
C
      CALL QENTER('TRDNEW')
      if (iprint .gt. 6) then
         write (lupri,*) 'OUTPUT from TRDNEW, iprint = ',iprint
         write (lupri,*) 'jcoor,jcorsy ',jcoor,jcorsy
         write (lupri,*) 'lwrk = ',lwrk
         CALL FLSHFO(LUPRI)
      end if
      KFRSAV = 1
      LFRSAV = LWRK
      KFREE  = KFRSAV
      LFREE  = LFRSAV
      CALL MEMGET('INTE',KINDAO,N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDAC,N2ASHX,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDRS,64,    WRK,KFREE,LFREE)
C
C NECyuji
      CALL IZERO(WRK(KINDRS),64)
      CALL TRDSET(JCORSY,WRK(KINDAO),WRK(KINDAC),WRK(KINDRS),IPRINT)
      if (iprint .gt. 6) then
         write (lupri,*) 'finished TRDSET'
         CALL FLSHFO(LUPRI)
      end if
C
C     Pack PV
C
      CALL MEMGET('REAL',KPVFPK,N2ASHX,WRK,KFREE,LFREE)
      if (iprint .gt. 25) then
         write (lupri,*) 'Packing PV matrix'
      end if
      DO 280 KY = 1,NASHT
         KYSYM = NSM(KY)
         DO 260 KX = 1,KY
            KXY = IROW(KY) + KX
            KXSYM = NSM(KX)
            KXYSYM = MULD2H(KXSYM,KYSYM)
C           ... KUVSYM = KXYSYM for PV(KU.KV,KX.KY) .ne. 0
C           Fold PV for FQ constructions
            CALL DSPTSI(NASHT,PV(1,KXY),WRK(KPVFPK))
            CALL DGEFSP(NASHT,WRK(KPVFPK),PV(1,KXY))
C
            CALL TRDPAK(PV(1,KXY),WRK(KPVFPK),NASH,IASH,NASHT,
     *                  KXYSYM,1)
            if (iprint .gt. 30) then
               write (lupri,*)
     *         'PV(u,v) for KX,KY,KXYSYM ',KX,KY,KXYSYM
               CALL OUTPAK(PV(1,KXY),NASHT,1,LUPRI)
               write (lupri,*)
     *         'PV(u,v) packed NDAC(KXYSYM) = ',NDAC(KXYSYM)
               CALL wrtmat(WRK(KPVFPK),1,NDAC(KXYSYM),1,NDAC(KXYSYM),0)
            end if
            CALL DCOPY(NDAC(KXYSYM),WRK(KPVFPK),1,PV(1,KXY),1)
  260    CONTINUE
  280 CONTINUE
      CALL MEMREL('TRDNEW.PVPAK',WRK,KFRSAV,KPVFPK,KFREE,LFREE)
      if (iprint .gt. 6) then
         write (lupri,*) 'finished pack PV'
         CALL FLSHFO(LUPRI)
      end if
C
      CALL MEMGET('REAL',KFQAOA,NBAST*NASHT,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDDI,N2BASX,WRK,KFREE,LFREE)
C     need LW1 in TRDDRV:
      LBINT = 600
C     record: BUF(LBINT),IBUF(LBINT),LENGTH
C      LBUF = LBINT*(IRAT+1) + 1
      LBUF = LBINT*(5*IRAT+1) + 1
      LW1  = MAX(IRAT*NBAST*NASHT,IRAT*2*MXNDAC,LBUF)
      LW1  = (LW1-1)/IRAT + 1
C
      MXDIS = (LFREE-LW1-10) / MXNDAO
      MXDIS = MIN(MXDIS,NNBASX)
      CALL MEMGET('REAL',KH2DBU,MXDIS*MXNDAO,WRK,KFREE,LFREE)
      CALL DZERO(WRK(KFQAOA),NBAST*NASHT)
      CALL TRDDRV(JCOOR,JCORSY,H2DAC,WRK(KFQAOA),CMO,PV,
     *            MXDIS,WRK(KH2DBU),
     *            WRK(KINDAO),WRK(KINDAC),WRK(KINDRS),WRK(KINDDI),
     *            WRK,KFREE,LFREE,IPRINT)
C
      if (iprint .gt. 4) then
         write (lupri,*) 'finished TRDDRV'
         if (iprint .gt. 5) then
            write (lupri,*) 'FQDAO(NBAST,NASHT) matrix:'
            CALL OUTPUT(WRK(KFQAOA),1,NBAST,1,NASHT,NBAST,NASHT,1,LUPRI)
         end if
         CALL FLSHFO(LUPRI)
      end if
      CALL MEMREL('TRDNEW.TRDDRV',WRK,KFRSAV,KINDDI,KFREE,LFREE)
C
C     Final transformation of first index in FQAOAC from
C     AO basis to MO basis
C
      CALL TRDFQ(JCORSY,FQD,WRK(KFQAOA),CMO)
      if (iprint .gt. 4) then
         write (lupri,*) 'finished TRDFQ'
         write (lupri,*) 'FQD(NORBT,NASHT) matrix:'
         CALL OUTPUT(FQD,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
         CALL FLSHFO(LUPRI)
      end if
C
C     921005-hjaaj: Because CISIGD is not modified for
C        symmetry packed H2AC yet, we must unpack H2DAC
C        and complete the matrix
C
      CALL TRDH2U(JCORSY,H2DAC,WRK(KINDAC),WRK(KFREE))
      if (iprint .gt. 4) then
         write (lupri,*) 'finished TRDH2U'
         if (iprint .ge. 5) then
            write (lupri,*) 'H2DAC matrix'
#if defined (VAR_H2DPACK)
            call output(h2dac,1,nnashx,1,nnashx,nnashx,nnashx,1,lupri)
#else
            call output(h2dac,1,n2ashx,1,n2ashx,n2ashx,n2ashx,1,lupri)
#endif
         end if
         if (iprint .ge. 10) then
            write (lupri,*) 'Antisymmetric component of H2DAC matrix'
#if defined (VAR_H2DPACK)
            nnashy = (nnashx*(nnashx+1))/2
            CALL MEMGET('REAL',KH2DAP,NNASHY,WRK,KFREE,LFREE)
            call dgetap(NNASHX,H2DAC,WRK(KH2DAP))
            call outpak(WRK(KH2DAP),NNASHX,1,LUPRI)
#else
            CALL MEMGET('REAL',KH2DAP,N2ASHX*N2ASHX,WRK,KFREE,LFREE)
            call dgetap(N2ASHX,H2DAC,WRK(KH2DAP))
            call outpak(WRK(KH2DAP),N2ASHX,1,LUPRI)
#endif
            CALL MEMREL('TRDNEW.H2DAP',WRK,KFRSAV,KH2DAP,KFREE,LFREE)
         end if
         CALL FLSHFO(LUPRI)
      end if
C
C     Unpack PV again
C
      CALL MEMGET('REAL',KPVFPK,N2ASHX,WRK,KFREE,LFREE)
      DO 880 KY = 1,NASHT
         KYSYM = NSM(KY)
         DO 860 KX = 1,KY
            KXY = IROW(KY) + KX
            KXSYM = NSM(KX)
            KXYSYM = MULD2H(KXSYM,KYSYM)
C           ... KUVSYM = KXYSYM for PV(KU.KV,KX.KY) .ne. 0
            CALL DCOPY(NDAC(KXYSYM),PV(1,KXY),1,WRK(KPVFPK),1)
            CALL TRDPAK(PV(1,KXY),WRK(KPVFPK),NASH,IASH,NASHT,
     *                  KXYSYM,-1)
C           Unfold PV again
            CALL DUNFLD(NASHT,PV(1,KXY),WRK(KPVFPK))
            CALL DSITSP(NASHT,WRK(KPVFPK),PV(1,KXY))
  860    CONTINUE
  880 CONTINUE
      if (iprint .gt. 4) then
         write (lupri,*) 'finished PV unpack'
         if (iprint .gt. 10) then
            write (lupri,*) 'PV matrix after unpack:'
            CALL OUTPUT(PV,1,NNASHX,1,NNASHX,NNASHX,NNASHX,1,LUPRI)
         end if
         CALL FLSHFO(LUPRI)
      end if
C     CALL MEMREL('TRDNEW.PVUPK',WRK,KFRSAV,KPVFPK,KFREE,LFREE)
C
      CALL MEMREL('TRDNEW',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TRDNEW')
      RETURN
      END
C  /* Deck trddrv */
      SUBROUTINE TRDDRV(JCOOR,JCORSY,H2DAC,FQAOAC,CMO,PVF, MXDIS,
     *                  H2DBUF,INDAO,INDAC,INDRS,INDDIS,
     *                  WRK,KFRSAV,LFRSAV,IPRINT)
C
C 2.Oct.1992 hjaaj
C
#include <implicit.h>
#include <iratdef.h>
C
C Used from common blocks:
C   INFORB : NCMOT,NBAST,NASHT,...
C   INFIND : ISAO(*)
C   INFTRL : MXNDAO
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <inftap.h>
#include <infind.h>
#include "inftrl.h"
C
      DIMENSION H2DAC(*),FQAOAC(NBAST,NASHT),CMO(NCMOT),PVF(*)
      DIMENSION H2DBUF(MXNDAO,*),INDDIS(NBAST,NBAST),WRK(*)
      DIMENSION INDAO(NBAST,NBAST), INDAC(NASHT,NASHT), INDRS(8,8)
C-B-HJTEMP-MAERKE
      PARAMETER (LBINT = 600)
C-E-HJTEMP-MAERKE
C
      KFREE = KFRSAV
      LFREE = LFRSAV
C
C     record: BUF(LBINT),IBUF(LBINT),LENGTH
      LBUF = LBINT*(5*IRAT+1) + 1
C     LW1  = maxdim(H2DWRK,H2UVRY and PVXR,intbuffer)
C            H2DWRK: NBAST*NASHT  if IABSYM .ne. 1
C                    MXNDAC+NBAST if IABSYM .eq. 1
C            H2UVRY,PVXR: MXNDAC
      LW1  = MAX(IRAT*NBAST*NASHT,IRAT*(MXNDAC+NBAST),LBUF)
      LW1  = (LW1-1)/IRAT + 1
      CALL MEMGET('REAL',KW1,LW1,WRK,KFREE,LFREE)
      KBUFW1 = KW1
      KIBUF  = KBUFW1 + LBINT*(IRAT + 1) + 1
      KUVRY  = KW1
      KPVXR  = KW1
C
      if (iprint .gt. 20) then
         DO 10 ISYM = 1,NSYM
            JCMO1 = ICMO(ISYM) + 1
            NBASI = NBAS(ISYM)
            NORBI = NORB(ISYM)
            write (lupri,*) 'TRDDRV: MO coefficients symmetry ',ISYM
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
C        H2DBUF(pq,inddis(r,s)) all (pq/ for given /rs)
C
         CALL TRFDIS(ICASE,LRL,LSL,INDRS,INDDIS,MXDIS,NDIS)
         IF (IPRINT .GE. 15) THEN
            WRITE (LUPRI,*) 'TRDDRV: ICASE,LRL,LSL =',ICASE,LRL,LSL
         END IF
         IF (ICASE .LT. 0) GO TO 9999
C        ... finished when ICASE .lt. 0 !
         NPASS = NPASS + 1
         if (iprint .gt. 4) then
            write (lupri,*) 'TRDDRV: LU2DER pass no.  ',NPASS
            write (lupri,*) 'number of distributions: ',NDIS
            IF (IPRINT .GE. 15) THEN
               WRITE (LUPRI,*) 'INDRS  array:'
               CALL IWRTMA(INDRS,NSYM,NSYM,8,8)
               IF (IPRINT .GE. 20) THEN
                  WRITE (LUPRI,*) 'INDDIS array:'
                  CALL IWRTMA(INDDIS,NBAST,NBAST,NBAST,NBAST)
            END IF
            END IF
            CALL FLSHFO(LUPRI)
         end if
C
         CALL DZERO(H2DBUF,NDIS*MXNDAO)
         CALL TRDRAO(LU2DER,JCOOR,JCORSY,WRK(KIBUF),H2DBUF,
     *               INDAO,INDDIS,WRK(KBUFW1))
         if (iprint .gt. 20) then
            write (lupri,*) 'TRDDRV: finished TRDRAO in pass no. ',NPASS
            CALL FLSHFO(LUPRI)
         end if
C
C        Transform (pq/rs) to (uv/rs) for the /rs) distributions
C        in memory.  H2DBUF(*,rs) is overwritten with (uv/rs)
C        integrals.
C
         CALL TRDPQ(JCORSY,H2DBUF,CMO,INDDIS,WRK(KW1),IPRINT)
C        CALL TRDPQ(JCORSY,H2DBUF,CMO,INDDIS,H2DWRK,IPRINT)
         if (iprint .gt. 20) then
            write (lupri,*) 'TRDDRV: finished TRDPQ in pass no. ',NPASS
            CALL FLSHFO(LUPRI)
         end if
C
C        Use (uv/rs) for FQD and H2DAC
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
               KXSYM = MULD2H(KSSYM,JCORSY)
               NASHX = NASH(KXSYM)
               IF (NASHX .EQ. 0 .AND. NASHS .EQ. 0) GO TO 700
C              Neither contribution for FQD or for H2DAC
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
               KUVSYM = MULD2H(KRSSYM,JCORSY)
C   also :     KUVSYM = MULD2H(KXSYM,KRSYM)
               NDACUV = NDAC(KUVSYM)
               IF (NDACUV .EQ. 0) GO TO 700
C              ... No act-act elements in this symmetry block
               if (iprint .gt. 10) then
                  write (lupri,*) 'KR,KRSYM,KSSYM,KRSSYM ',
     &               KR,KRSYM,KSSYM,KRSSYM
                  write (lupri,*) 'KUVSYM,NDACUV ',KUVSYM,NDACUV
                  write (lupri,*) 'IDOH2, IDOFQ  ',IDOH2, IDOFQ
               end if
C
C              First FQD contributions
C              NOTE: r and s are switched in this section,
C                    distributions are numbered (uv / sr)
C
C  INDSR .ne. 0:
C    FQAOAC(s,x) += sum(uv) (uv/sr) * PV(uv,rx)
C                           H2DBUF(uv,abs(indsr)) * PVXR(uv)
C
C
C
               IF (NASHX .EQ. 0 .OR. IDOFQ .EQ. 0) GO TO 600
               IASHX = IASH(KXSYM)
               KS1 = IBASS + 1
               KSL = IBASS + NBASS
               DO 480 KX = IASHX+1,IASHX+NASHX
                  CALL TRDPV(KX,KR,KRSYM,CMO(ICMORA),NBASR,NASHR,
     *                       KUVSYM,PVF,WRK(KPVXR))
                  DO 460 KS = KS1,KSL
                     INDSR = ABS( INDDIS(KS,KR) )
                     IF (INDSR .GT. 0) THEN
                        FQAOAC(KS,KX) = FQAOAC(KS,KX) +
     *                     DDOT(NDACUV,H2DBUF(1, INDSR),1,WRK(KPVXR),1)
                     END IF
  460             CONTINUE
                  if (iprint .gt. 25) then
                     write (lupri,*) 'PVXR matrix for KR,KX = ',KR,KX
                     CALL OUTPUT(WRK(KPVXR),1,NDACUV,1,1,MXNDAC,1,1,
     &                           LUPRI)
                     write (lupri,*) 'FQAOAC for KSSYM,KX = ',KSSYM,KX
                     CALL OUTPUT(FQAOAC(1,1),KS1,KSL,KX,KX,
     &                           NBAST,NASHT,1,LUPRI)
                  end if
  480          CONTINUE
C
C
C              Second and last : H2DAC contributions
C
  600       IF (NASHS .EQ. 0 .OR. IDOH2 .EQ. 0) GO TO 700
               NBASS  = NBAS(KSSYM)
               ICMOSA = ICMO(KSSYM) + NISH(KSSYM)*NBASS + 1
               CALL TRDAH2(KR,KRSYM,KSSYM,KUVSYM,
     *              CMO(ICMORA),CMO(ICMOSA),NBASR,NASHR,NBASS,NASHS,
     *              INDDIS,H2DBUF,H2DAC,WRK(KUVRY),IPRINT)
  700       CONTINUE
  800    CONTINUE
         if (iprint .gt. 20) then
            write (lupri,*)
     &       'TRDDRV: finished transformation in pass no. ',NPASS
            CALL FLSHFO(LUPRI)
         end if
C
      GO TO 100
C
 9999 CONTINUE
      IF (IPRINT .GE. 2) THEN
         WRITE (LUPRI,'(/A,I5,A)')
     &    'TRDDRV: finished transformation in',NPASS,
     &    ' passes over derivative AO integral file.'
         CALL FLSHFO(LUPRI)
      END IF
      RETURN
      END
C  /* Deck trdset */
      SUBROUTINE TRDSET(JH2DSY,INDAO,INDAC,INDRS,IPRINT)
C
C Copyright 12-Jun-1994 Hans Joergen Aa. Jensen
C
C Set index information INDAO and INDAC and
C set index information for matrices in /INFTRL/.
C Set INDRS(rsym,ssym) = 0 for sym.dist. desired
C                      = 654321 for sym.dist. not needed
C
C IH1XX(IBSYM,IABSYM) is off-set for block (IASYM,IBSYM) in H1XX
C                     of symmetry IABSYM (IASYM = MULD2H(IBSYM,IABSYM))
C
#include <implicit.h>
C
C Used from common blocks
C   INFORB: NSYM, NASH(8), NBAS(8), MULD2H(8,8)
C   INFTRL: IH1AO(8,8),IH1AC(8,8),NDAO(8), NDAC(8)
C
#include <priunit.h>
#include <inforb.h>
#include "inftrl.h"
C
      DIMENSION INDAO(NBAST,NBAST), INDAC(NASHT,NASHT), INDRS(8,8)
C
      PARAMETER (IMEMERR = -999999999)
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
            IF (IBSYM .LT. IASYM) THEN
               IH1AO(IBSYM+1,IABSYM) = IH1AO(IBSYM,IABSYM) +
     *            NBAS(IASYM)*NBAS(IBSYM)
               IH1AC(IBSYM+1,IABSYM) = IH1AC(IBSYM,IABSYM) +
     *            NASH(IASYM)*NASH(IBSYM)
            ELSE IF (IBSYM .EQ. IASYM) THEN
               IH1AO(IBSYM+1,IABSYM) = IH1AO(IBSYM,IABSYM) +
     *            (NBAS(IBSYM)*(NBAS(IBSYM)+1))/2
               IH1AC(IBSYM+1,IABSYM) = IH1AC(IBSYM,IABSYM) +
     *            (NASH(IBSYM)*(NASH(IBSYM)+1))/2
            ELSE
C           ... IBSYM .GT. IASYM temporary assignment,
C               changed in 420 loop
               IH1AO(IBSYM+1,IABSYM) = IH1AO(IBSYM,IABSYM)
               IH1AC(IBSYM+1,IABSYM) = IH1AC(IBSYM,IABSYM)
            END IF
            END IF
C           ... 2. INDAO
            DO 140 JA = IBAS(IASYM)+1,IBAS(IASYM)+NBAS(IASYM)
               JBEND = MIN(JA,IBAS(IBSYM)+NBAS(IBSYM))
               DO 120 JB = IBAS(IBSYM)+1,JBEND
                  JNDAO = JNDAO + 1
                  INDAO(JA,JB) = JNDAO
                  INDAO(JB,JA) = JNDAO
  120          CONTINUE
  140       CONTINUE
            NDAO(IABSYM) = JNDAO
            MXNDAO = MAX(MXNDAO,JNDAO)
C           ... 3. INDAC
            DO 240 JA = IASH(IASYM)+1,IASH(IASYM)+NASH(IASYM)
               JBEND = MIN(JA,IASH(IBSYM)+NASH(IBSYM))
               DO 220 JB = IASH(IBSYM)+1,JBEND
                  JNDAC = JNDAC + 1
                  INDAC(JA,JB) = JNDAC
                  INDAC(JB,JA) = JNDAC
  220          CONTINUE
  240       CONTINUE
            NDAC(IABSYM) = JNDAC
            MXNDAC = MAX(MXNDAC,JNDAC)
  400    CONTINUE
C
C        reset values for upper triangle to give
C        memory access error (bus errer, segmentation violation)
C        if used.
C
         DO 420 IBSYM = 1,NSYM
            IASYM = MULD2H(IABSYM,IBSYM)
            IF (IBSYM .GT. IASYM) THEN
               IH1AO(IBSYM,IABSYM) = IMEMERR
               IH1AC(IBSYM,IABSYM) = IMEMERR
            END IF
  420    CONTINUE
  500 CONTINUE
C
      DO 840 IBSYM = 1,NSYM
         DO 820 IASYM = 1,IBSYM
            IABSYM = MULD2H(IASYM,IBSYM)
            IUVSYM = MULD2H(IABSYM,JH2DSY)
            IF (NDAC(IUVSYM) .GT. 0 .AND.
     &          (NASH(IASYM) .GT. 0 .OR. NASH(IBSYM) .GT. 0) ) THEN
               INDRS(IASYM,IBSYM) = 0
               INDRS(IBSYM,IASYM) = 0
            ELSE
               INDRS(IASYM,IBSYM) = 654321
               INDRS(IBSYM,IASYM) = 654321
            END IF
  820    CONTINUE
  840 CONTINUE
C
      IF (IPRINT .GT. 4) THEN
         write (lupri,*) 'output from TRDSET'
         WRITE (lupri,*) 'NASH(*) ',(NASH(I),I=1,NSYM)
         WRITE (lupri,*) 'NBAS(*) ',(NBAS(I),I=1,NSYM)
         write (lupri,*) 'MXNDAC,MXNDAO ',MXNDAC,MXNDAO
         write (lupri,*) 'NDAC(*) ',(NDAC(I),I=1,NSYM)
         write (lupri,*) 'NDAO(*) ',(NDAO(I),I=1,NSYM)
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
C  /* Deck trdpq */
      SUBROUTINE TRDPQ(JCORSY,H2DBUF,CMO,INDDIS,H2DWRK,IPRINT)
C
C 6-Oct-1992 Hans Joergen Aa. Jensen
C
C Transform (pq/rs) to (uv/rs) for the available
C /rs) distributions according to INDDIS(r,s).
C
C Input : H2DBUF(pq,indrs) = (pq / rs)
C Output: H2DBUF(uv,indrs) = (uv / rs)
C
#include <implicit.h>
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
      DIMENSION H2DBUF(MXNDAO,*), CMO(NCMOT)
      DIMENSION INDDIS(NBAST,NBAST), H2DWRK(*)
C
      DO 800 KR = 1,NBAST
         KRSYM = ISAO(KR)
         DO 700 KS = 1,NBAST
            INDRS = INDDIS(KR,KS)
            IF (INDRS .GT. 0) THEN
               KSSYM = ISAO(KS)
               IRSSYM = MULD2H(KRSYM,KSSYM)
               IPQSYM = MULD2H(IRSSYM,JCORSY)
               IF (IPRINT .GT. 25) THEN
                  NDAOPQ = NDAO(IPQSYM)
                  WRITE (LUPRI,*) 'H2CDAO(pq,rs) for r,s = ',KR,KS
                  WRITE (LUPRI,*) 'indrs,sym(pq), ndao(pqsym) ',
     *               INDRS,IPQSYM,NDAOPQ
                  CALL WRTMAT(H2DBUF(1,INDRS),1,NDAOPQ,1,NDAOPQ,0)
               END IF
               CALL TRDTDS(IPQSYM,H2DBUF(1,INDRS),H2DBUF(1,INDRS),
     *                     CMO,H2DWRK)
C              Note: equivalence(H2CDAO,H2CDAC)
C              CALL TRDTDS(IABSYM,H2CDAO,H2CDAC,CMO,WRK)
               IF (IPRINT .GT. 25) THEN
                  NDACUV = NDAC(IPQSYM)
                  WRITE (LUPRI,*) 'H2CDAC(uv,rs) for r,s = ',KR,KS
                  WRITE (LUPRI,*) 'sym(uv), ndac(uvsym)',IPQSYM,NDACUV
                  CALL WRTMAT(H2DBUF(1,INDRS),1,NDACUV,1,NDACUV,0)
               END IF
            END IF
  700    CONTINUE
  800 CONTINUE
C
      RETURN
      END
C  /* Deck trdtds */
      SUBROUTINE TRDTDS(IABSYM,H2CDAO,H2CDAC,CMO,WRK)
C
C Copyright (c) 4-Jul-1994/2-Jan-1995 Hans Joergen Aa. Jensen
C 950102 revision : improved vectorization (now not NASHA in inner loop)
C
C Transform H2CDAO(nbast,nbast) of symmetry IABSYM
C        to H2CDAC(nasht,nasht) using CMO.
C H2CDAO(a,b) .eq. H2CDAO(b,a) is assumed.
C H2CDAO and H2CDAC are both symmetry packed.
C
C H2CDAO and H2CDAC may share same memory location.
C
C Input : H2CDAO(a.b) = (a b | c d)
C         CMO = m.o. coefficients
C         IABSYM = symmetry of H2CDAO, H2CDAC matrices
C Output: H2CDAC(u.v) = (u v | c d) =
C         sum(b) [ sum(a) [ CMO(a,u) H2CDAO(a.b) ] CMO(b,v) ]
C
#include <implicit.h>
C
C Used from common blocks:
C  INFORB : NSYM, MULD2H, NASH(*), NBAS(*), NISH(*), ICMO(*)
C  INFTRL : IH1AO(8,8), IH1AC(8,8)
C
#include <inforb.h>
#include "inftrl.h"
C
      DIMENSION CMO(NCMOT), H2CDAO(*), H2CDAC(*), WRK(*)
C
C     Allocation for IABSYM .eq. 1
      KAC = 1
      KWRK = KAC + NDAC(1)
C
      DO 100 IBSYM = 1,NSYM
         IASYM = MULD2H(IBSYM,IABSYM)
      IF (IBSYM .GT. IASYM) GO TO 100
C        ... only saving lower triangle
         NASHB = NASH(IBSYM)
         NASHA = NASH(IASYM)
         IF ((NASHA.NE.0) .AND. (NASHB.NE.0)) THEN
            NBASB = NBAS(IBSYM)
            NBASA = NBAS(IASYM)
            JCMOB = ICMO(IBSYM) + NISH(IBSYM)*NBASB + 1
            JCMOA = ICMO(IASYM) + NISH(IASYM)*NBASA + 1
            JH2AO = IH1AO(IBSYM,IABSYM) + 1
            JH2AC = IH1AC(IBSYM,IABSYM) + 1
         IF (IABSYM .EQ. 1) THEN
C           ... temporary matrix WRK(KAC) is used in order that
C               H2CDAO and H2CDAC may share same memory location.
            CALL UTHU(H2CDAO(JH2AO),WRK(KAC),CMO(JCMOA),
     *                WRK(KWRK),NBASA,NASHA)
C           CALL UTHU(H,HT,U,WRK,NAO,NMO)
            NH2AC = (NASHA*(NASHA+1))/2
            CALL DCOPY(NH2AC,WRK(KAC),1,H2CDAC(JH2AC),1)
         ELSE
            CALL DGEMM('N','N',NBASB,NASHA,NBASA,1.D0,
     &                 H2CDAO(JH2AO),NBASB,
     &                 CMO(JCMOA),NBASA,0.D0,
     &                 WRK,NBASB)
            CALL DGEMM('T','N',NASHB,NASHA,NBASB,1.D0,
     &                 CMO(JCMOB),NBASB,
     &                 WRK,NBASB,0.D0,
     &                 H2CDAC(JH2AC),NASHB)
         END IF
         END IF
  100 CONTINUE
C
C     End of TRDTDS
C
      RETURN
      END
C  /* Deck trdrao */
      SUBROUTINE TRDRAO(LUDERA,JCOOR,JCORSY,IBUF,H2DAO,INDAO,INDDIS,WRK)
C
C 25-Sep-1992 hjaaj / revised 6-Oct-1992
C
C Reads two-electron derivative integrals from file to H2DAO buffer.
C LUDERA: File to read
C JCOOR : =1,2,3 for d/dx,d/dy,d/dz respectively.
C BUF   : Points to start element in integral buffer.
C IBUF  : Points to start element in label buffer.
C H2DAO : H2DAO(rs,pq), pq=(mpqoff+1:mpqoff+npq); dim (mxndao,npq)
C
#include <implicit.h>
#include <iratdef.h>
#include <priunit.h>
#include <mxcent.h>
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
      DIMENSION H2DAO(MXNDAO,*), IBUF(4,600), WRK(600*(IRAT + 1) + 1)
      DIMENSION INDAO(NBAST,NBAST), INDDIS(NBAST,NBAST)
C
C     Local variables:
C
      CHARACTER*8 KEY
      LOGICAL DOCOOR
C
#include <ibtdef.h>
C
      CALL REWSPL(LUDERA)
C     KEY = 'AO2DRINT'
C     CALL MOLLAB(KEY,LUDERA,LUPRI)
#if defined (VAR_INT64)
      NIBUF = 1
      NBITS = 16
#else
      IF (NBASIS .LE. 255) THEN
         NIBUF = 1
         NBITS = 8
      ELSE
         NIBUF = 2
         NBITS = 16
         call quit('TRDRAO only implement for .le. 255 basis functions')
      END IF
#endif
      LBUF = 600
      JBUF = 256
C
      LENINT = IRAT*LBUF + LBUF + 1
      KINT  = 1
      KIINT = KINT + LBUF
C
      DOCOOR = .FALSE.
 150  CONTINUE
#if defined (VAR_SPLITFILES)         
         CALL READSI(LUDERA,LENINT,WRK(KINT),JBUF)
#else
         CALL READI(LUDERA,LENINT,WRK(KINT))
#endif
         CALL AOLAB4(WRK(KIINT),LBUF,NIBUF,NBITS,IBUF,LENGTH)
         IF (LENGTH .GT. 0) THEN
            DO 100 I = 1, LENGTH
               JS    = IBUF(4,I)
               IF (JS .EQ. 0) THEN
                  ICOOR  = IBUF(3,I)
                  ICORSY = IBUF(2,I) + 1
                  JP = IBUF(1,I)
                  DOCOOR = ICOOR.EQ.JCOOR .AND. JCORSY.EQ.ICORSY
               ELSE IF (DOCOOR) THEN
                  JR = IBUF(3,I)
                  JP = IBUF(1,I)
                  JQ = IBUF(2,I)
C                 BUF(I) =  (P Q | R S)
                  JRS = ABS( INDDIS(JR,JS) )
                  IF (JRS .GT. 0) THEN
                     H2DAO(INDAO(JP,JQ),JRS) = WRK(KINT + I - 1)
                  END IF
                  JPQ = ABS( INDDIS(JP,JQ) )
                  IF (JPQ .GT. 0) THEN
                     H2DAO(INDAO(JR,JS),JPQ) = WRK(KINT + I - 1)
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
C  /* Deck trfdis */
      SUBROUTINE TRFDIS(ICASE,LRL,LSL,INDRS,INDDIS,MXDIS,NDISX)
C
C Copyright 28-Sep-1992 Hans Joergen Aa. Jensen
C
C General routine for finding next set of distributions to read into
C memory.
C
C Initializations before first call:
C   INDRS = 0      for desired symmetry distributions
C         = 654321 for sym.dist. not needed
C   ICASE .le. 0: new round
C
#include <implicit.h>
C
C Used from common blocks:
C   INFORB : NBAST, N2BASX
C
#include <inforb.h>
C
      DIMENSION INDRS(8,8), INDDIS(NBAST,NBAST)
C
      IF (ICASE .LE. 0) THEN
C     ... if (first)
         ICASE = 0
      END IF
      ICASE = ICASE + 1
      CALL IZERO(INDDIS,N2BASX)
C
      NDIS = 0
      DO 800 KRSYM = 1, NSYM
         DO 700 KSSYM = 1, KRSYM
C           information about symmetry distributions:
C           INDRS(KRSYM,KSSYM) .GT. 0: finished previously or not needed
C                              .LT. 0: partly done previously
C                               EQ. 0: not treated yet
C
C           NOTE: only one element of INDRS(KRSYM,KSSYM) may be
C                 .lt. 0, or the LRL,LSL tests will not work!
C
         IF (INDRS(KRSYM,KSSYM) .GT. 0) GO TO 700
            KR1 = IBAS(KRSYM) + 1
            KRL = IBAS(KRSYM) + NBAS(KRSYM)
            KS1 = IBAS(KSSYM) + 1
            KS1X= KS1
            KSL = IBAS(KSSYM) + NBAS(KSSYM)
            IF (INDRS(KRSYM,KSSYM) .LT. 0) THEN
               IF (LRL .EQ. 0 .AND. LSL .EQ. 0) THEN
C              ... this sym.dist. was finished last time.
                  INDRS(KRSYM,KSSYM) = -INDRS(KRSYM,KSSYM)
                  GO TO 700
               END IF
               IF (LSL .EQ. KSL .OR. LSL .EQ. LRL) THEN
                  KR1 = LRL + 1
               ELSE
                  KR1 = LRL
                  KS1X= LSL + 1
               END IF
               INDRS(KRSYM,KSSYM) = -ICASE
            ELSE
C           ... INDRS(KRSYM,KSSYM) .eq. 0
               IF (KSSYM .EQ. KRSYM) THEN
                  NRSDIS = ( NBAS(KRSYM)*NBAS(KRSYM) + NBAS(KRSYM) ) / 2
               ELSE
                  NRSDIS = NBAS(KRSYM)*NBAS(KSSYM)
               END IF
               IF (NDIS+NRSDIS .GT. MXDIS) THEN
                  IF (NDIS .GT. 0) GO TO 700
C                 only complete sym.dist. if we already have something
                  INDRS(KRSYM,KSSYM) = -ICASE
               ELSE
                  INDRS(KRSYM,KSSYM) = ICASE
               END IF
            END IF
C
            DO 300 KR = KR1,KRL
               IF (KR .GT. KR1) KS1X = KS1
               KSLX = MIN(KR,KSL)
               DO 200 KS = KS1X,KSLX
                  NDIS = NDIS + 1
                  INDDIS(KR,KS) = NDIS
                  IF (KR .NE. KS) INDDIS(KS,KR) = -NDIS
                  IF (NDIS .EQ. MXDIS) THEN
                     LRL = KR
                     LSL = KS
                     GO TO 9999
                  END IF
  200          CONTINUE
               IF (KR .LT. KRL .AND. KRSYM .NE. KSSYM) THEN
                  IF ( NDIS .GT. 0 .AND.
     &                (NDIS + NBAS(KSSYM)) .GT. MXDIS ) THEN
                     LRL = KR
                     LSL = KSLX
                     GO TO 9999
                  END IF
               END IF
  300       CONTINUE
C           mark that this sym.dist. has been completed
            LRL = 0
            LSL = 0
  700    CONTINUE
  800 CONTINUE
C
C
      IF (NDIS .EQ. 0) ICASE = -1
C     ... finished !
C
 9999 CONTINUE
      NDISX = NDIS
      RETURN
      END
C  /* Deck trdah2 */
      SUBROUTINE TRDAH2(KR,KRSYM,KSSYM,KUVSYM,
     *                  CMOR,CMOS,NBASR,NASHR,NBASS,NASHS,
     *                  INDDIS,H2UVRS,H2DAC,H2UVRY,IPRINT)
C
C 5-Oct-1992 Hans Joergen Aagaard Jensen
C 'AH2' for Add to H2DAC
C Add (uv / rs) contribution to (uv / xy) (x .ge. y)
C
#include <implicit.h>
C
C Used from common blocks:
C   INFORB : ??
C   INFTRL : NDAC(8),MXNDAO,NXNDAC
C
#include <priunit.h>
#include <inforb.h>
#include "inftrl.h"
C
      DIMENSION CMOR(NBASR,NASHR), CMOS(NBASS,NASHS)
      DIMENSION INDDIS(NBAST,NBAST)
#if defined (VAR_H2DPACK)
      DIMENSION H2UVRS(MXNDAO,*), H2DAC(NNASHX,NNASHX)
#else
      DIMENSION H2UVRS(MXNDAO,*), H2DAC(N2ASHX,NASHT,NASHT)
#endif
      DIMENSION H2UVRY(MXNDAC)
      PARAMETER ( DP5 = 0.5D0, D0 = 0.0D0 )
C
#if defined (VAR_H2DPACK)
      IROW(I) = (I*(I-1))/2
#endif
C
      IF (IPRINT .GT. 25) THEN
         write (lupri,*) 'output from trdah2, kr,kssym = ',kr,kssym
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
            write (lupri,*) 'TRDAH2 H2UVRY for r,y = ',KR,KY
            CALL OUTPUT(H2UVRY,1,NDACUV,1,1,NDACUV,1,1,LUPRI)
         END IF
C
C        If XY_sym = 1 (and thus RS_sym = 1) then
C        we fold the x .lt. y into x .gt. y.
C        The x .eq. y term is missing a factor of 2 in this way,
C        this is fixed in TRDH2U
C        C_x
C
         DO 300 JX = 1,NASHR
            CMORX = CMOR(JR,JX)
            IF (CMORX .NE. D0) THEN
               KX = IASH(KRSYM) + JX
#if defined (VAR_H2DPACK)
               IF (KX .GT. KY) THEN
                  KXY = IROW(KX) + KY
               ELSE
                  KXY = IROW(KY) + KX
               END IF
               CALL DAXPY(NDACUV,CMORX,H2UVRY,1,H2DAC(1,KXY),1)
#else
               IF (KX .GT. KY) THEN
                  CALL DAXPY(NDACUV,CMORX,H2UVRY,1,H2DAC(1,KX,KY),1)
               ELSE
                  CALL DAXPY(NDACUV,CMORX,H2UVRY,1,H2DAC(1,KY,KX),1)
               END IF
#endif
            END IF
  300    CONTINUE
  400 CONTINUE
      RETURN
C     end of TRDAH2
      END
C  /* Deck trdh2u */
      SUBROUTINE TRDH2U(JCORSY,H2DAC,INDAC,WRK)
C
C 6-Oct-1992 hjaaj
C
#include <implicit.h>
C
C Used from common blocks:
C  INFORB : NASHT,NASH(),IASH()
C  INFIND : NSM()
C  INFTRL : NDAC()
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include "inftrl.h"
C
#if defined (VAR_H2DPACK)
      DIMENSION H2DAC(NNASHX,NNASHX)
#else
      DIMENSION H2DAC(N2ASHX,NASHT,NASHT)
#endif
      DIMENSION INDAC(NASHT,NASHT), WRK(*)
C
      PARAMETER (IWAY = -1, D2 = 2.0D0)
C     ... IWAY=-1 to unpack
      DO 800 KY = 1,NASHT
         KYSYM = NSM(KY)
         DO 600 KX = 1,KY
            KXSYM = NSM(KX)
            KXYSYM = MULD2H(KXSYM,KYSYM)
            KUVSYM = MULD2H(KXYSYM,JCORSY)
C           finish H2DAX(1,X.X) terms with missing factor 2
C           (stemming from triangular code in TRDAH2)
#if defined (VAR_H2DPACK)
            KXY = IROW(KY) + KX
            IF (KX .EQ. KY) THEN
               NDACUV = NDAC(KUVSYM)
               CALL DSCAL(NDACUV,D2,H2DAC(1,KXY),1)
            END IF
            CALL TRDPAK(WRK,H2DAC(1,KXY),
     *                  NASH,IASH,NASHT,KUVSYM,IWAY)
C           CALL TRDPAK(H1SP,H1PK,NBLK,IBLK,NBLKT,IH1SYM,IWAY)
            CALL DCOPY(NNASHX,WRK,1,H2DAC(1,KXY),1)
#else
            IF (KX .EQ. KY) THEN
               NDACUV = NDAC(KUVSYM)
               CALL DSCAL(NDACUV,D2,H2DAC(1,KX,KX),1)
            END IF
            CALL TRDPAK(WRK,H2DAC(1,KY,KX),
     *                  NASH,IASH,NASHT,KUVSYM,IWAY)
C           CALL TRDPAK(H1SP,H1PK,NBLK,IBLK,NBLKT,IH1SYM,IWAY)
C           square distribution into H2DAC
            CALL DSPTSI(NASHT,WRK,H2DAC(1,KY,KX))
C           CALL DSPTSI(NASHT,ASP,ASI)
            IF (KX .NE. KY) THEN
               CALL DCOPY(N2ASHX,H2DAC(1,KY,KX),1,H2DAC(1,KX,KY),1)
            END IF
#endif
C
  600    CONTINUE
  800 CONTINUE
      RETURN
      END
C  /* Deck trdpv */
      SUBROUTINE TRDPV(KX,KS,KSSYM,CMOS,NBASS,NASHS,
     *                 KUVSYM,PVF,PVXS)
C
C Copyright 6-Oct-1992 Hans Joergen Aa. Jensen
C
C PVXS(uv) =  sum(y) PVF(uv,x,y) * CMO(s,y) =  [uv|xs]
C
#include <implicit.h>
C
C Used from common blocks:
C  INFORB : NASHT
C  INFTRL : NDAC(8)
C
#include <inforb.h>
#include "inftrl.h"
C
      DIMENSION CMOS(NBASS,NASHS), PVF(NNASHX,NNASHX)
      DIMENSION PVXS(MXNDAC)
C
      NDACUV = NDAC(KUVSYM)
      CALL DZERO(PVXS,MXNDAC)
      IASHS  = IASH(KSSYM)
      JS     = KS - IORB(KSSYM)
C     Note: if PVF packed KX .ge. KY we need to check
C           and use PVF(1,KXY) transposed if KX .LT. KY
      DO 800 JY = 1,NASHS
         KY = JY + IASHS
         IF (KX .GT. KY) THEN
            KXY = KX*(KX-1)/2 + KY
         ELSE
            KXY = KY*(KY-1)/2 + KX
         END IF
         CALL DAXPY(NDACUV,CMOS(JS,JY),PVF(1,KXY),1,PVXS,1)
  800 CONTINUE
      RETURN
      END
C  /* Deck trdfq */
      SUBROUTINE TRDFQ(JCORSY,FQD,FQAOAC,CMO)
C
C 7-Oct-1992 hjaaj
C
C FQD(a,x) = sum(r) FQAOAC(r,x) * CMO(r,a)
C
#include <implicit.h>
C
C Used from common blocks:
C  INFORB : NORBT,NASHT,NBAST,NCMOT,NSYM
C
#include <inforb.h>
C
      DIMENSION FQD(NORBT,NASHT), FQAOAC(NBAST,NASHT), CMO(NCMOT)
C
      DO 800 KSYM = 1,NSYM
         NASHK  = NASH(KSYM)
      IF (NASHK.EQ.0) GO TO 800
         IASYM  = MULD2H(KSYM,JCORSY)
         NORBA  = NORB(IASYM)
         NBASA  = NBAS(IASYM)
         IASHK1 = IASH(KSYM) + 1
         IORBA1 = IORB(IASYM) + 1
         IBASA1 = IBAS(IASYM) + 1
         JCMOA1 = ICMO(IASYM) + 1
         IF (NORBA*NBASA .NE. 0) 
     &        CALL DGEMM('T','N',NORBA,NASHK,NBASA,1.D0,
     &                   CMO(JCMOA1),NBASA,
     &                   FQAOAC(IBASA1,IASHK1),NBAST,0.D0,
     &                   FQD(IORBA1,IASHK1),NORBT)
  800 CONTINUE
      RETURN
      END
C  /* Deck trdpak */
      SUBROUTINE TRDPAK(H1SP,H1PK,NBLK,IBLK,NBLKT,IH1SYM,IWAY)
C
C Copyright 4-Jul-1994 Hans Joergen Aa. Jensen
C
C NBLK,IBLK,NBLKT may e.g. be NORB(),IORB(),NORBT or NASH(),IASH(),NASHT
C
C IWAY .ge. 0 : Pack H1SP in H1PK
C      else   : Unpack H1PK in H1SP
C
#include <implicit.h>
C
C Used from common blocks:
C   INFORB : NSYM,MULD2H
C
#include <inforb.h>
C
      DIMENSION H1SP(*), H1PK(*), NBLK(8), IBLK(8)
C
      IF (IH1SYM .EQ. 1) THEN
         CALL PKSYM1(H1SP,H1PK,NBLK,NSYM,IWAY)
C        CALL PKSYM1(ASP,APK,NDIM,NBLK,IWAY)
C        dimension NDIM(NBLK)
         GO TO 9999
      END IF
      NNBLKX = (NBLKT*(NBLKT+1))/2
      IF (IWAY .LT. 0) CALL DZERO(H1SP,NNBLKX)
      IABPK = 1
      DO 300 IBSYM = 1,NSYM
         IASYM  = MULD2H(IBSYM,IH1SYM)
      IF (IBSYM .GT. IASYM) GO TO 300
         NBLKA  = NBLK(IASYM)
         NBLKB  = NBLK(IBSYM)
      IF (NBLKA .EQ. 0 .OR. NBLKB .EQ. 0) GO TO 300
         IBLKA  = IBLK(IASYM)
         IB1    = IBLK(IBSYM) + 1
         DO 200 IA = IBLKA+1,IBLKA+NBLKA
            IABSP  = (IA*(IA-1))/2 + IB1
            IF (IWAY .GE. 0) THEN
               CALL DCOPY(NBLKB,H1SP(IABSP),1,H1PK(IABPK),1)
            ELSE
               CALL DCOPY(NBLKB,H1PK(IABPK),1,H1SP(IABSP),1)
            END IF
            IABPK = IABPK + NBLKB
  200    CONTINUE
  300 CONTINUE
C
C
 9999 RETURN
      END
