C
C...   Copyright (c) 1997 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dalton, an ab initio electronic structure program, Release 1.0
C...   (1997), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, T. Andersen, K. L. Bak, V. Bakken,
C...   O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
C...   H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi, H. Koch,
C...   K. V. Mikkelsen, P. Norman, M. J. Packer, T.Saue,
C...   P. R. Taylor, and O. Vahtras"
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
C
C$Id: eri2out.F,v 1.2 2001-02-12 18:17:54 vebjornb Exp $
C
#include <single.h>
C
C  /* Deck eriout */
      SUBROUTINE ERIOUT(SO,ITYPE,IPNTCR,IODDCC,IPNTUV,BUF,IBUF,INDXBT,
     &                  NCOUNT,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <eridst.h>
      DIMENSION SO(NAOINT,*), 
     &          IPNTCR(MAXBCH,4), IPNTUV(KC2MAX,0:NRDER,2),
     &          WORK(LWORK), IODDCC(NRTOP), BUF(*), IBUF(*),
     &          NINTS(MXDIST), NCOUNT(*), INDXBT(*)
#include <cbieri.h>
#include <ericom.h>
#include <eribuf.h>
#include <symmet.h>
#include <hertop.h>
#include <inftap.h>
C
C     Allocation for ERIOUT
C
      IF (LAST) THEN
         LBIN   = 0
         KBIN   = 1
         KIBIN  = 1
         KSOPRT = 1
         KNDORB = 1
         KLABEL = 1
         KPNCHN = 1
         KNDALL = 1
         KNDXTR = 1
         KJCMPA = 1
         KJCMPC = 1
         KJADDA = 1
         KJADDC = 1
         KIOFCM = 1
         KLAST  = 1
      ELSE IF (GDER .OR. BDER) THEN
         LBIN   = NCCT*KHKTA*KHKTB*KHKTC*KHKTD
C
         KBIN   = 1
         KIBIN  = KBIN   +  LBIN
         KSOPRT = KIBIN  + (NIBUF*LBIN - 1)/IRAT + 1
         KNDORB = KSOPRT + NAOINT 
         KLABEL = KNDORB + (NIBUF*LBIN - 1)/IRAT + 1
         KPNCHN = KLABEL + (4*NAOINT - 1)/IRAT + 1
         KNDALL = KPNCHN + (LBIN - 1)/IRAT + 1
         KNDXTR = KNDALL + (2*NCCS - 1)/IRAT + 1
         KJCMPA = KNDXTR + (4*NCCS*(MAXREP + 1) - 1)/IRAT + 1
         KJCMPC = KJCMPA + 64*MXAQN*MXAQN
         KJADDA = KJCMPC + 64*MXAQN*MXAQN
         KJADDC = KJADDA + 64*MXAQN*MXAQN
         KIOFCM = KJADDC + 64*MXAQN*MXAQN
         KLAST  = KIOFCM + MXAQN*MXAQN*MXAQN*MXAQN
      ELSE IF (DODIST) THEN
         LBIN   = NCCT*KHKTA*KHKTB*KHKTC*KHKTD
         KBIN   = 1
         KIBIN  = KBIN   +  LBIN
         KSOPRT = KIBIN  + (NIBUF*LBIN - 1)/IRAT + 1
         KNDORB = KSOPRT 
         KLABEL = KNDORB + (NIBUF*LBIN - 1)/IRAT + 1
         KPNCHN = KLABEL 
         KNDALL = KPNCHN + (LBIN - 1)/IRAT + 1
         KNDXTR = KNDALL + (2*NCCS - 1)/IRAT + 1
         KJCMPA = KNDXTR + (4*NCCS*(MAXREP + 1) - 1)/IRAT + 1
         KJCMPC = KJCMPA + 64*MXAQN*MXAQN
         KJADDA = KJCMPC + 64*MXAQN*MXAQN
         KJADDC = KJADDA + 64*MXAQN*MXAQN
         KIOFCM = KJADDC + 64*MXAQN*MXAQN
         KLAST  = KIOFCM + MXAQN*MXAQN*MXAQN*MXAQN
      ELSE
         LBIN   = NCCT*KHKTA*KHKTB*KHKTC*KHKTD
         IF (BDER) LBIN = 2*LBIN + 1
         KBIN   = 1
         KIBIN  = KBIN   +  LBIN
         KSOPRT = KIBIN  + (NIBUF*LBIN - 1)/IRAT + 1
         KNDORB = KSOPRT
         KLABEL = KNDORB + (NIBUF*LBIN - 1)/IRAT + 1
         KPNCHN = KLABEL
         KNDALL = KPNCHN
         KNDXTR = KNDALL
         KJCMPA = KNDXTR
         KJCMPC = KJCMPA
         KJADDA = KJCMPC
         KJADDC = KJADDA
         KIOFCM = KJADDC
         KLAST  = KIOFCM
      END IF
      IF (KLAST .GT. LWORK) CALL STOPIT('ERIOUT',' ',KLAST,LWORK)
      IF ((GDER .OR. BDER) .AND. UNDIFF) THEN
         IF (GDER) N0 = 13 
         IF (BDER) N0 =  7 
         CALL ERIOU1(SO,SO(1,N0),ITYPE,WORK(KNDORB),WORK(KLABEL),
     &               WORK(KPNCHN),WORK(KNDALL),WORK(KNDXTR),IPNTCR,
     &               IODDCC,IPNTUV,NINTS,NCOUNT,WORK(KBIN),WORK(KIBIN),
     &               LBIN,BUF,IBUF,WORK(KJCMPA),WORK(KJCMPC),
     &               WORK(KJADDA),WORK(KJADDC),WORK(KIOFCM),INDXBT,
     &               IPRINT)
      ELSE
         CALL ERIOU1(SO,WORK(KSOPRT),ITYPE,WORK(KNDORB),WORK(KLABEL),
     &               WORK(KPNCHN),WORK(KNDALL),WORK(KNDXTR),IPNTCR,
     &               IODDCC,IPNTUV,NINTS,NCOUNT,WORK(KBIN),WORK(KIBIN),
     &               LBIN,BUF,IBUF,WORK(KJCMPA),WORK(KJCMPC),
     &               WORK(KJADDA),WORK(KJADDC),WORK(KIOFCM),INDXBT,
     &               IPRINT)
      END IF
C
      RETURN
      END
C  /* Deck eriou1 */
      SUBROUTINE ERIOU1(SO,SOPRT,ITYPE,INDORB,LABEL,IPNCHN,INDALL,
     &                  INDXTR,IPNTCR,IODDCC,IPNTUV,NINTS,NCOUNT,
     &                  BIN,IBIN,LBIN,BUF,IBUF,JCMPAB,JCMPCD,JADDAB,
     &                  JADDCD,IOFCMP,INDXBT,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <dummy.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
      DIMENSION SO(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD,*), SOPRT(*),
     &          INDORB(LBIN,NIBUF), IPNCHN(*), INDALL(*), INDXTR(*),
     &          LABEL(*), IPNTCR(MAXBCH,4),
     &          IODDCC(NRTOP), IPNTUV(KC2MAX,0:NRDER,2),
     &          BUF(LBUF,NBUFS,0:NCORS), IBUF(LBUF,NIBUF,NBUFS,0:NCORS),
     &          BIN(LBIN), IBIN(LBIN,NIBUF),
     &          NINTS(NDISTR), NCOUNT(NDISTR,0:NCORS), 
     &          JCMPAB(*), JCMPCD(*),
     &          JADDAB(*), JADDCD(*), IOFCMP(*), INDXBT(*)
#include <cbieri.h>
#include <ericom.h>
#include <erithr.h>
#include <nuclei.h>
#include <eribuf.h>
#include <aobtch.h>
#include <hertop.h>
#include <eridst.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERIOU1',-1)
C
C     Last call to dump buffers
C     =========================
C
      IF (LAST) THEN
         IF (GDER .OR. BDER) THEN
            ISTART = 1
            IF (UNDIFF) ISTART = 0
            DO ISCOOR = ISTART, NPERTS 
               CALL ERIWRT(ISCOOR,NINTS,NCOUNT,BIN,IBIN,LBIN,
     &                     BUF(1,1,ISCOOR),IBUF(1,1,1,ISCOOR),IPRINT)
            END DO
         ELSE
             CALL ERIWRT(0,NINTS,NCOUNT,BIN,IBIN,LBIN,BUF(1,1,0),
     &                   IBUF(1,1,1,0),IPRINT)
         END IF
C
C     Geometrical derivatives in distributions
C     ========================================
C
      ELSE IF (GDER .OR. BDER) THEN
C
C        orbital labels and centers
C
         CALL ERILAB(DUMMY,LABEL,INDORB,IPNTCR,IODDCC,IPNTUV,INDXBT,
     &               INDXTR,IPRINT)
C
C        loop over Cartesian directions
C
         NBCHS = MLTPX*KHKTB*KHKTCD
C
         ISTART = 1
         IF (UNDIFF) ISTART = 0
         DO ISCOOR = ISTART, NPERTS 
            IATOM = (ISCOOR + 2)/3
            ICART = MOD(ISCOOR - 1,3) + 1
            IF (IPRINT.GT.10) THEN
               WRITE (LUPRI,'(2X,A,I5)') ' ISCOOR in ERIOU1 ', ISCOOR
            END IF
C
            IRPDST = 0 
C
C           construct derivatives
C
            IF (ISCOOR .GT. 0) THEN
               IF (GDER) THEN
                  CALL ERIPRT(SO,ISCOOR,IATOM,ICART,SOPRT,LABEL,INDORB,
     &                        IPRINT)
               ELSE
                  CALL ERIBRT(SO,ISCOOR,SOPRT,IPRINT)
               END IF
            END IF
C
C           sort distributions
C
            CALL ERIGTH(BIN,IBIN,LBIN,NIBUF,INDORB,SOPRT,IPNTCR,IODDCC,
     &                  IPNTUV,IPNCHN,NBITS,NINTS,NBCHS,IOFCMP,INDXBT,
     &                  ICART,IPRINT)
C
C           write integrals
C
            CALL ERIWRT(ISCOOR,NINTS,NCOUNT,BIN,IBIN,LBIN,
     &                  BUF(1,1,ISCOOR),IBUF(1,1,1,ISCOOR),IPRINT)
         END DO
C
C     Undifferentiated integrals in distributions 
C     ===========================================
C
      ELSE IF (DODIST) THEN
C
C        oribtal labels and integrals
C
         CALL ERINDD(SO,INDORB,IPNCHN,INDALL,INDXTR,IPNTCR,IODDCC,
     &               IPNTUV,BIN,IBIN,NINTS,LBIN,NIBUF,NBITS,JCMPAB,
     &               JCMPCD,JADDAB,JADDCD,IOFCMP,INDXBT,IPRINT)
C
C        write out
C
         CALL ERIWRT(0,NINTS,NCOUNT,BIN,IBIN,LBIN,BUF(1,1,0),
     &               IBUF(1,1,1,0),IPRINT)
C
C     Undifferentiated integrals (all)
C     ================================
C
      ELSE
         CALL ERINDX(SO,ITYPE,INDORB,IPNTCR,IODDCC,IPNTUV,BIN,IBIN,
     &               NINTS,LBIN,NIBUF,INDXBT,NBITS,IPRINT)
         CALL ERIWRT(0,NINTS,NCOUNT,BIN,IBIN,LBIN,BUF(1,1,0),
     &               IBUF(1,1,1,0),IPRINT)
      END IF
C
      RETURN
      END
C  /* Deck eriwrt */
      SUBROUTINE ERIWRT(ISCOOR,NINTS,NCOUNT,BIN,IBIN,LBIN,BUF,IBUF,
     &                  IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <dummy.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
C
#include <eridst.h>
      DIMENSION BUF(LBUF,NBUFS), IBUF(LBUF,NIBUF,NBUFS),
     &          BIN(LBIN), IBIN(LBIN,NIBUF),
     &          NINTS(NDISTR), NCOUNT(NDISTR,0:NCORS)
      CHARACTER*8 FAODER
#include <cbieri.h>
#include <chrnos.h>
#include <nuclei.h>
#include <ericom.h>
#include <erithr.h>
#include <eribuf.h>
#include <aobtch.h>
#include <hertop.h>
#include <inftap.h>
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERIWRT',-1)
C
      IBIT1 = 2**NBITS     - 1
      IBIT2 = 2**(2*NBITS) - 1
C
C     Initialization of LUINTX
C     ========================
C
      IF (FIRST) THEN
         IF (DODIST) THEN
            IF (LUINTR .LE. 0) THEN
               CALL GPOPEN(LUINTR,'AOTWODIS','UNKNOWN',' ',' ',
     *                     IDUMMY,.FALSE.)
            ENDIF
            REWIND LUINTR
         ELSE IF (GDER .OR. BDER) THEN
            DO I = 1, NPERTS 
               LUINTD(I) = 0
               FAODER = 'AO2DER'//CHRNOS(I/10)
     &                         //CHRNOS(MOD(I,10))
               CALL GPOPEN(LUINTD(I),FAODER,'UNKNOWN',' ',
     &                     'UNFORMATTED',IDUMMY,.FALSE.)
               CALL REWSPL(LUINTD(I))
               CALL NEWLAB(FAODER,LUINTD(I),LUPRI)
            END DO
C        ELSE IF (BDER) THEN
C           LUINTX = LU2DER 
C           CALL GPOPEN(LUINTX,'AO2MGINT','UNKNOWN',' ',
C    &                  'UNFORMATTED',IDUMMY,.FALSE.)
C           CALL REWSPL(LUINTX)
C           CALL NEWLAB('AO2MGINT',LUINTX,LUPRI)
         ELSE
            LUINTX = LUINTA 
            CALL REWSPL(LUINTX)
            CALL NEWLAB('BASTWOEL',LUINTX,LUPRI)
         END IF
         CALL IZERO(NCOUNT(1,0),NDISTR*NCORS)
         CALL IZERO(NBUFX(0),NCORS)
         FIRST = .FALSE.
      END IF
C
      IF (GDER .OR. BDER) THEN
         LUINTX = LUINTD(ISCOOR)
      END IF
C
C     Run over distributions
C     ======================
C
      IF (.NOT.LAST) THEN
         IOFF = 0
         DO 400 IDIST = 1, NDISTR
            ICOUNT = NCOUNT(IDIST,ISCOOR)
            NTOT   = NINTS(IDIST)
C
C           Empty whole buffers directly from BIN and IBIN
C           ==============================================
C
            KBUFS = NTOT/LBUF
            DO 500 I = 1, KBUFS
               NBUFX(ISCOOR) = NBUFX(ISCOOR) + 1
               NBUF = NBUFX(ISCOOR)
               KSTR = IOFF + (I - 1)*LBUF + 1
               KEND = IOFF + I*LBUF
               IF (.NOT.NOWRIT) THEN
                  IF (DODIST) THEN
                     WRITE (LUINTR) INDDST(IDIST), ISCOOR 
                     WRITE (LUAORC(ISCOOR),REC=NBUF)
     &                     (  BIN(K  ),K=KSTR,KEND),
     &                     ((IBIN(K,L),K=KSTR,KEND),L=1,NIBUF),
     &                     LBUF
                  ELSE
                     WRITE (LUINTX)
     &                     (  BIN(K  ),K=KSTR,KEND),
     &                     ((IBIN(K,L),K=KSTR,KEND),L=1,NIBUF),
     &                     LBUF
                  END IF
               END IF
               IF (INTPRI) THEN
                  CALL ERIPRI(BIN(KSTR),IBIN(KSTR,1),LBIN,NIBUF,
     &                        LBUF,NBITS,NDISTR,IDIST,ISCOOR)
               END IF
  500       CONTINUE
C
C           Case A: Transfer remainder to buffer
C           ====================================
C
            ISTR = KBUFS*LBUF + 1
C           IF (BDER.AND.KBUFS.GT.0) THEN
C              ISTR = ISTR - 1
C              BIN (ISTR) =  BIN(1)
C              IBIN(ISTR,1) = IBIN(1,1)
C           END IF 
            IF (ICOUNT + (NTOT - ISTR + 1) .LT. LBUF) THEN
               IF (NIBUF .EQ. 1) THEN
                  DO 611 I = IOFF + ISTR, IOFF + NTOT
                     ICOUNT = ICOUNT + 1
                     BUF (ICOUNT,  IDIST) =  BIN(I)
                     IBUF(ICOUNT,1,IDIST) = IBIN(I,1)
  611             CONTINUE
               ELSE
                  DO 612 I = IOFF + ISTR, IOFF + NTOT
                     ICOUNT = ICOUNT + 1
                     BUF (ICOUNT,  IDIST) =  BIN(I)
                     IBUF(ICOUNT,1,IDIST) = IBIN(I,1)
                     IBUF(ICOUNT,2,IDIST) = IBIN(I,2)
  612             CONTINUE
               END IF
C
C           Case B: Transfer remainder to buffer and empty when full
C           ========================================================
C
            ELSE
               NLEFT = LBUF - ICOUNT
               ILAST = ISTR + NLEFT - 1
               IF (NIBUF .EQ. 1) THEN
                  DO 621 I = IOFF + ISTR, IOFF + ILAST
                     ICOUNT = ICOUNT + 1
                     BUF (ICOUNT,  IDIST) =  BIN(I)
                     IBUF(ICOUNT,1,IDIST) = IBIN(I,1)
  621             CONTINUE
               ELSE
                  DO 622 I = IOFF + ISTR, IOFF + ILAST
                     ICOUNT = ICOUNT + 1
                     BUF (ICOUNT,  IDIST) =  BIN(I)
                     IBUF(ICOUNT,1,IDIST) = IBIN(I,1)
                     IBUF(ICOUNT,2,IDIST) = IBIN(I,2)
  622             CONTINUE
               END IF
               NBUFX(ISCOOR) = NBUFX(ISCOOR) + 1
               NBUF = NBUFX(ISCOOR)
               IF (.NOT.NOWRIT) THEN
                  IF (DODIST) THEN
                     WRITE (LUINTR) INDDST(IDIST), ISCOOR
                     WRITE (LUAORC(ISCOOR), REC=NBUF)
     &                     (  BUF(K,  IDIST),K=1,LBUF),
     &                     ((IBUF(K,L,IDIST),K=1,LBUF),L=1,NIBUF),
     &                     LBUF
                  ELSE
                     WRITE (LUINTX)
     &                     (  BUF(K,  IDIST),K=1,LBUF),
     &                     ((IBUF(K,L,IDIST),K=1,LBUF),L=1,NIBUF),
     &                     LBUF
                  END IF
               END IF
               IF (INTPRI) THEN
                  CALL ERIPRI(BUF(1,IDIST),IBUF(1,1,IDIST),LBUF,NIBUF,
     &                        LBUF,NBITS,NDISTR,IDIST,ISCOOR)
               END IF
               ISTR = ILAST + 1
C              IF (BDER) THEN
C                 ISTR = ISTR - 1
C                 BIN (ISTR) =  BIN(1)
C                 IBIN(ISTR,1) = IBIN(1,1)
C              END IF 
               ICOUNT = 0
               IF (NIBUF .EQ. 1) THEN
                  DO 631 I = IOFF + ISTR, IOFF + NTOT
                     ICOUNT = ICOUNT + 1
                     BUF (ICOUNT,  IDIST) =  BIN(I)
                     IBUF(ICOUNT,1,IDIST) = IBIN(I,1)
  631             CONTINUE
               ELSE
                  DO 632 I = IOFF + ISTR, IOFF + NTOT
                     ICOUNT = ICOUNT + 1
                     BUF (ICOUNT,  IDIST) =  BIN(I)
                     IBUF(ICOUNT,1,IDIST) = IBIN(I,1)
                     IBUF(ICOUNT,2,IDIST) = IBIN(I,2)
  632             CONTINUE
               END IF
            END IF
            NCOUNT(IDIST,ISCOOR) = ICOUNT
            IOFF = IOFF + NTOT
  400    CONTINUE
      END IF
C
C     Empty last buffer
C     =================
C
      IF (LAST) THEN
C
         DO 700 IDIST = 1, NDISTR
            ICOUNT = NCOUNT(IDIST,ISCOOR)
            IF (ICOUNT .GT. 0) THEN
               NBUFX(ISCOOR) = NBUFX(ISCOOR) + 1
               NBUF = NBUFX(ISCOOR)
               IF (.NOT.NOWRIT) THEN
                  IF (DODIST) THEN
                     WRITE (LUINTR) INDDST(IDIST), ISCOOR
                     WRITE (LUAORC(ISCOOR),REC=NBUF)
     &                     (  BUF(K,  IDIST),K=1,LBUF),
     &                     ((IBUF(K,L,IDIST),K=1,LBUF),L=1,NIBUF),
     &                     ICOUNT
                  ELSE
                     WRITE (LUINTX)
     &                     (  BUF(K,  IDIST),K=1,LBUF),
     &                     ((IBUF(K,L,IDIST),K=1,LBUF),L=1,NIBUF),
     &                     ICOUNT
                  END IF
               END IF
               IF (INTPRI) THEN
                  CALL ERIPRI(BUF(1,IDIST),IBUF(1,1,IDIST),LBUF,NIBUF,
     &                        ICOUNT,NBITS,NDISTR,IDIST,ISCOOR)
               END IF
            END IF
  700    CONTINUE
         IF (.NOT.DODIST) THEN
            WRITE (LUINTX) (  BUF(K,1),  K=1,LBUF),
     &                     ((IBUF(K,L,1),K=1,LBUF),L=1,NIBUF),-1
            REWIND LUINTX
            IF (.NOT.BDER) CALL GPCLOSE(LUINTX,'KEEP')
         END IF
C
      END IF
C
      RETURN
      END
C  /* Deck erican */
      SUBROUTINE ERICAN(IBIN,INTS,LBIN,NIBUF,NBITS)
#include <implicit.h>
#include <priunit.h>
#include <ibtpar.h>
      DIMENSION IBIN(LBIN,NIBUF)
#include <ibtfun.h>
C
      IBIT1 = 2**   NBITS  - 1
      IBIT2 = 2**(2*NBITS) - 1
C
      IF (NIBUF .EQ. 1) THEN
         DO 100 INT = 1, INTS
            I = IBTAND(IBTSHR(IBIN(INT,1),3*NBITS),IBIT1)
            J = IBTAND(IBTSHR(IBIN(INT,1),2*NBITS),IBIT1)
            K = IBTAND(IBTSHR(IBIN(INT,1),  NBITS),IBIT1)
            L = IBTAND(       IBIN(INT,1),         IBIT1)
            IJ = MAX(I,J)*IBIT1 + I + J
            KL = MAX(K,L)*IBIT1 + K + L
            IBIN(INT,1) = MAX(IJ,KL)*IBIT2 + IJ + KL
  100    CONTINUE
      ELSE
         DO 200 INT = 1, INTS
            I = IBTAND(IBTSHR(IBIN(INT,1),NBITS),IBIT1)
            J = IBTAND(       IBIN(INT,1),       IBIT1)
            K = IBTAND(IBTSHR(IBIN(INT,2),NBITS),IBIT1)
            L = IBTAND(       IBIN(INT,2),       IBIT1)
            IJ = MAX(I,J)*IBIT1 + I + J
            KL = MAX(K,L)*IBIT1 + K + L
            IBIN(INT,1) = MAX(IJ,KL)
            IBIN(INT,2) = MIN(IJ,KL)
  200    CONTINUE
      END IF
      RETURN
      END
C  /* Deck eripri */
      SUBROUTINE ERIPRI(BUF,IBUF,NDIM,NIBUF,ICOUNT,NBITS,NDISTR,IDIST,
     &                  ISCOOR)
#include <implicit.h>
#include <priunit.h>
#include <ibtpar.h>
      DIMENSION BUF(NDIM), IBUF(NDIM,NIBUF)
#include <ibtfun.h>
      IBIT1 = 2**NBITS - 1
      DO 100 INT = 1, ICOUNT
         IF (NIBUF .EQ. 1) THEN
            I = IBTAND(IBTSHR(IBUF(INT,1),3*NBITS),IBIT1)
            J = IBTAND(IBTSHR(IBUF(INT,1),2*NBITS),IBIT1)
            K = IBTAND(IBTSHR(IBUF(INT,1),  NBITS),IBIT1)
            L = IBTAND(       IBUF(INT,1),         IBIT1)
         ELSE
            I = IBTAND(IBTSHR(IBUF(INT,1),NBITS),IBIT1)
            J = IBTAND(       IBUF(INT,1),       IBIT1)
            K = IBTAND(IBTSHR(IBUF(INT,2),NBITS),IBIT1)
            L = IBTAND(       IBUF(INT,2),       IBIT1)
         END IF
         IF (NDISTR .GT. 1) THEN
            WRITE (LUPRI,'(10X,A,2X,I5,5X,4I4,5X,1P,D16.8,5X,I5)')
     &                   ' @@ ', ISCOOR, I, J, K, L, BUF(INT), IDIST
         ELSE
            WRITE (LUPRI,'(10X,A,2X,I5,5X,4I4,5X,1P,D16.8)')
     &                   ' @@ ', ISCOOR, I, J, K, L, BUF(INT)
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck erindx */
      SUBROUTINE ERINDX(SO,ITYPE,INDORB,IPNTCR,IODDCC,IPNTUV,
     &                  BIN,IBIN,INT,LBIN,NIBUF,INDXBT,NBITS,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
      PARAMETER (D0 = 0.0D0, TEN14 = 1.1D14)
      INTEGER A, B, C, D, I, A1, B1, C1, D1, AB, CD, R, S, T, E
      LOGICAL DOREP(0:7,4), CTRIAB, CTRICD, CTRIAC, CTRIBD, CTRIPQ
      DIMENSION SO(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD,*),
     &          INDORB(NCCS,NIBUF),
     &          IPNTCR(MAXBCH,4),
     &          IODDCC(NRTOP), IPNTUV(KC2MAX,0:NRDER,2),
     &          BIN(LBIN), IBIN(LBIN,NIBUF),
     &          IPNRST(0:7,3), INDXBT(MXSHEL*MXCONT,0:7),
     &          IADCMP(MXAQN,MXAQN,2)
#include <cbieri.h>
#include <ericom.h>
#include <erithr.h>
#include <aobtch.h>
#include <hertop.h>
#include <symmet.h>
C
#include <ibtfun.h>
      IBTEST(I,J,K,L) = IBTAND(I,IBTXOR(J,ISYMAO(K,L)))
      IBTREP(I,J,K,L) = IBTXOR(J,IBTAND(I,ISYMAO(K,L)))
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERINDX',-1)
C
      IF (IPRINT .GT. 15) THEN
         WRITE (LUPRI,'(2X,A,2I5)') ' NPQBCX,NPPBCX ',NPQBCX,NPPBCX
         WRITE (LUPRI,'(2X,A,2I5)') ' NPQBCS,NPPBCS ',NPQBCS,NPPBCS
      END IF
C
      CALL PRPREP(DOREP(0,1),NHKTA,KHKTA,ISTBLA)
      CALL PRPREP(DOREP(0,2),NHKTB,KHKTB,ISTBLB)
      CALL PRPREP(DOREP(0,3),NHKTC,KHKTC,ISTBLC)
      CALL PRPREP(DOREP(0,4),NHKTD,KHKTD,ISTBLD)
C
      CALL CMPADR(IADCMP(1,1,1),KHKTA,KHKTB,TKMPAB)
      CALL CMPADR(IADCMP(1,1,2),KHKTC,KHKTD,TKMPCD)
C
      CALL GETRST(IPNRST(0,1),ISTBLR)
      CALL GETRST(IPNRST(0,2),ISTBLS)
      CALL GETRST(IPNRST(0,3),ISTBLT)
C
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(/,2X,A,8L2)')'DOREP A  ',(DOREP(I,1),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP B  ',(DOREP(I,2),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP C  ',(DOREP(I,3),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP D  ',(DOREP(I,4),I=0,MAXREP)
      END IF
C
      NCCS1 = (NPQBCS - NPPBCS)*NCTFAB*NCTFCD
      NCCS2 = NPPBCS*NCTFAB*NCTFCD
C
      IF (NIBUF .EQ. 1) THEN
         IBITA = 2**(3*NBITS)
         IBITB = 2**(2*NBITS)
         IBITC = 2**(  NBITS)
         IBITD = 1
      ELSE
         IBITA = 2**NBITS
         IBITB = 1
         IBITC = 2**NBITS
         IBITD = 1
      END IF
C
      INT = 0
      IF (BDER) THEN
         IATOM = 0
         IREPE = 0
         BIN(1)  = TEN14
         IF (NIBUF.EQ.1) THEN
            IBIN(1,1) = IATOM*2**24+IREPE*2**16+ITYPE*2**8
         ELSE
            IBIN(1,1) = IATOM*2**24+IREPE*2**16+ITYPE*2**8
            IBIN(1,2) = IATOM*2**24+IREPE*2**16+ITYPE*2**8
         END IF
         IF (NIBUF.NE.1) THEN
            WRITE (LUPRI,'(1X,A,/,1X,A)') 
     &         ' ERINDX has only been implemented for NIBUF.EQ.1.',
     &         ' Program aborted.'
            CALL QUIT('ERINDX implemented only for NIBUF.EQ.1')
         END IF 
         INT = 1
         ICRB = ITYPE
         E = ISYMAX(ICRB,2)
      ELSE
         ICRB = 0
         E = 0
      END IF 
C
      IF (NCCS1 .GT. 0) THEN
C
         DO 100 A = 0, MAXREP
         IF (DOREP(A,1)) THEN
         DO 110 B = 0, MAXREP
         IF (DOREP(B,2)) THEN
         DO 120 C = 0, MAXREP
         IF (DOREP(C,3) .AND. 
     &       DOREP(IBTXOR(IBTXOR(IBTXOR(A,B),C),E),4)) THEN
            D = IBTXOR(IBTXOR(IBTXOR(A,B),C),E)
            CD = IBTXOR(C,D)
C
            IF (DIAGAB .AND. B.GT.A) GO TO 190
            IF (DIAGCD .AND. D.GT.C) GO TO 190
C
            IF (NIBUF .EQ. 1) THEN
               DO 201 I = 1, NCCS1
                  INDORB(I,1) = (INDXBT(KNDXBT(IPNTCR(I,1)),A)-1)*IBITA
     &                        + (INDXBT(KNDXBT(IPNTCR(I,2)),B)-1)*IBITB
     &                        + (INDXBT(KNDXBT(IPNTCR(I,3)),C)-1)*IBITC
     &                        + (INDXBT(KNDXBT(IPNTCR(I,4)),D)-1)
  201          CONTINUE
            ELSE
               DO 202 I = 1, NCCS1
                  INDORB(I,1) = (INDXBT(KNDXBT(IPNTCR(I,1)),A)-1)*IBITA
     &                        + (INDXBT(KNDXBT(IPNTCR(I,2)),B)-1)
                  INDORB(I,2) = (INDXBT(KNDXBT(IPNTCR(I,3)),C)-1)*IBITC
     &                        + (INDXBT(KNDXBT(IPNTCR(I,4)),D)-1)
  202          CONTINUE
            END IF
C
            R = IPNRST(B,1)
            S = IPNRST(D,2)
            T = IPNRST(CD,3)
C
            IA   = 0
            MAXB = KHKTB
            MAXD = KHKTD
            CTRIAB = DIAGAB .AND. A.EQ.B
            CTRICD = DIAGCD .AND. C.EQ.D
            DO 300 ICMPA = 1, KHKTA
            IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
            IF (IVARA.EQ.0) THEN
               IA  = IA + IBITA
               IAB = IA
               IF (CTRIAB) MAXB = ICMPA
               DO 310 ICMPB = 1, MAXB
               IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
               IF (IVARB.EQ.0) THEN
                  IAB  = IAB + IBITB
                  IF (NIBUF.EQ.1) THEN
                     IABC = IAB
                  ELSE
                     IABC = 0
                  END IF
                  ICMPAB = IADCMP(ICMPA,ICMPB,1)
                  IODDAB = IODDCC(IPNTUV(ICMPAB,ICRB,1))
                  IF (BDER) THEN
                     JODDAB = IODDCC(IPNTUV(ICMPAB,0,1))
                  ELSE
                     JODDAB = IODDAB 
                  END IF
                  DO 320 ICMPC = 1, KHKTC
                  IVARC = IBTEST(ISTBLC,C,NHKTC,ICMPC)
                  IF (IVARC.EQ.0) THEN
                     IABC  = IABC + IBITC
                     IABCD = IABC
                     IF (CTRICD) MAXD = ICMPC
                     DO 330 ICMPD = 1, MAXD
                     IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                     IF (IVARD.EQ.0) THEN
                        IABCD  = IABCD + 1
                        ICMPCD = IADCMP(ICMPC,ICMPD,2)
                        IODDCD = IODDCC(IPNTUV(ICMPCD,0,2))
                        IF (IODDAB .EQ. IODDCD .OR.
     &                      JODDAB .EQ. IODDCD) THEN
C
C                       Standard integrals
C
                        IF (.NOT.BDER) THEN
                           IF (WRTSCR) THEN
                              IF (NIBUF .EQ. 1) THEN
                                 DO 401 I = 1, NCCS1
                                    SOABCD = SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    IF (ABS(SOABCD) .GT. THRSH) THEN
                                       INT = INT + 1
                                       BIN (INT  ) = SOABCD
                                       IBIN(INT,1) = INDORB(I,1) + IABCD
                                    END IF
  401                            CONTINUE
                              ELSE
                                 DO 402 I = 1, NCCS1
                                    SOABCD = SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    IF (ABS(SOABCD) .GT. THRSH) THEN
                                       INT = INT + 1
                                       BIN (INT  ) = SOABCD
                                       IBIN(INT,1) = INDORB(I,1) + IAB
                                       IBIN(INT,2) = INDORB(I,2) + IABCD
                                    END IF
  402                            CONTINUE
                              END IF
                           ELSE
                              IF (NIBUF .EQ. 1) THEN
                                 DO 451 I = 1, NCCS1
                                    INT = INT + 1
                                    BIN(INT)=SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    IBIN(INT,1) = INDORB(I,1) + IABCD
  451                            CONTINUE
                              ELSE
                                 DO 452 I = 1, NCCS1
                                    INT = INT + 1
                                    BIN(INT)=SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    IBIN(INT,1) = INDORB(I,1) + IAB
                                    IBIN(INT,2) = INDORB(I,2) + IABCD
  452                            CONTINUE
                              END IF
                           END IF
C
C                      B derivative integrals
C
                        ELSE
                           IF (WRTSCR) THEN
                              IF (NIBUF .EQ. 1) THEN
                                 DO 403 I = 1, NCCS1
                                    SO1  = SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    SO2  = SO(I,R,S,T,ICMPAB,ICMPCD,4)
                                    SOAB = SO1 + SO2
                                    SOBA = SO1 - SO2
                                    IJKL = INDORB(I,1) + IABCD
                                    IJLK = INDBEX(IJKL) 
                                    IF (ABS(SOAB) .GT. THRSH) THEN
                                       INT = INT + 1
                                       BIN (INT  ) = SOAB 
                                       IBIN(INT,1) = IJKL 
                                    END IF
                                    IF (IJLK.GT.0 .AND. 
     &                                  ABS(SOBA) .GT. THRSH) THEN
                                       INT = INT + 1
                                       BIN (INT  ) = SOBA 
                                       IBIN(INT,1) = IJLK
                                    END IF
  403                            CONTINUE
                              ELSE
                                 DO 404 I = 1, NCCS1
                                    SO1 = SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    SO2 = SO(I,R,S,T,ICMPAB,ICMPCD,4)
                                    SOAB = SO1 + SO2
                                    SOBA = SO1 - SO2
                                    IF (ABS(SOAB) .GT. THRSH) THEN
                                       INT = INT + 1
                                       BIN (INT  ) = SOAB
                                       IBIN(INT,1) = INDORB(I,1) + IAB
                                       IBIN(INT,2) = INDORB(I,2) + IABCD
                                    END IF
                                    IF (ABS(SOBA) .GT. THRSH) THEN
                                       INT = INT + 1
                                       BIN (INT  ) = SOBA
                                       IBIN(INT,1) = INDORB(I,1) + IAB
                                       IBIN(INT,2) = INDORB(I,2) + IABCD
                                    END IF
  404                            CONTINUE
                              END IF
                           ELSE
                              IF (NIBUF .EQ. 1) THEN
                                 DO 453 I = 1, NCCS1
                                    SO1 = SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    SO2 = SO(I,R,S,T,ICMPAB,ICMPCD,4)
                                    IJKL = INDORB(I,1) + IABCD
                                    IJLK = INDBEX(IJKL) 
                                    INT = INT + 1
                                    BIN(INT) = SO1 + SO2 
                                    IBIN(INT,1) = IJKL 
                                    IF (IJLK.GT.0) THEN
                                       INT = INT + 1
                                       BIN(INT) = SO1 - SO2 
                                       IBIN(INT,1) = IJLK
                                    END IF
  453                            CONTINUE
                              ELSE
                                 DO 454 I = 1, NCCS1
                                    SO1 = SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    SO2 = SO(I,R,S,T,ICMPAB,ICMPCD,4)
                                    INT = INT + 1
                                    BIN(INT) = SO1 + SO2 
                                    IBIN(INT,1) = INDORB(I,1) + IAB
                                    IBIN(INT,2) = INDORB(I,2) + IABCD
                                    INT = INT + 1
                                    BIN(INT) = SO1 - SO2 
                                    IBIN(INT,1) = INDORB(I,1) + IAB
                                    IBIN(INT,2) = INDORB(I,2) + IABCD
  454                            CONTINUE
                              END IF
                           END IF
                        END IF 
                        END IF
                     END IF
  330                CONTINUE
                  END IF
  320             CONTINUE
               END IF
  310          CONTINUE
            END IF
  300       CONTINUE
C
  190       CONTINUE
C
         END IF
  120    CONTINUE
         END IF
  110    CONTINUE
         END IF
  100    CONTINUE
      END IF
C
      IF (NCCS2 .GT. 0) THEN
C
         DO 105 A = 0, MAXREP
         IF (DOREP(A,1)) THEN
         DO 115 B = 0, MAXREP
         IF (DOREP(B,2)) THEN
         DO 125 C = 0, MAXREP
         IF (DOREP(C,3) .AND. 
     &       DOREP(IBTXOR(IBTXOR(IBTXOR(A,B),C),E),4)) THEN
            D = IBTXOR(IBTXOR(IBTXOR(A,B),C),E)
            CD = IBTXOR(C,D)
C
            IF (DIAGAB .AND. B.GT.A) GO TO 195
            IF (DIAGCD .AND. D.GT.C) GO TO 195
            IF (C.GT.A .OR. (C.EQ.A .AND. D.GT.B)) GO TO 195
C
            IF (NIBUF .EQ. 1) THEN
               DO 206 I = NCCS1 + 1, NCCS
                  INDORB(I,1) = (INDXBT(KNDXBT(IPNTCR(I,1)),A)-1)*IBITA
     &                        + (INDXBT(KNDXBT(IPNTCR(I,2)),B)-1)*IBITB
     &                        + (INDXBT(KNDXBT(IPNTCR(I,3)),C)-1)*IBITC
     &                        + (INDXBT(KNDXBT(IPNTCR(I,4)),D)-1)
  206          CONTINUE
            ELSE
               DO 207 I = NCCS1 + 1, NCCS
                  INDORB(I,1) = (INDXBT(KNDXBT(IPNTCR(I,1)),A)-1)*IBITA
     &                        + (INDXBT(KNDXBT(IPNTCR(I,2)),B)-1)
                  INDORB(I,2) = (INDXBT(KNDXBT(IPNTCR(I,3)),C)-1)*IBITC
     &                        + (INDXBT(KNDXBT(IPNTCR(I,4)),D)-1)
  207          CONTINUE
            END IF
C
            R = IPNRST(B,1)
            S = IPNRST(D,2)
            T = IPNRST(CD,3)
C
            MAXB = KHKTB
            MAXC = KHKTC
            MAXD = KHKTD
            CTRIAB = DIAGAB .AND. A.EQ.B
            CTRICD = DIAGCD .AND. C.EQ.D
            CTRIAC = A.EQ.C
            CTRIBD = B.EQ.D
            CTRIPQ = A.EQ.C .AND. B.EQ.D
C
            IA   = 0
            DO 305 ICMPA = 1, KHKTA
            IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
            IF (IVARA.EQ.0) THEN
               IA  = IA + IBITA
               IAB = IA
               IF (CTRIAB) MAXB = ICMPA
               DO 315 ICMPB = 1, MAXB
               IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
               IF (IVARB.EQ.0) THEN
                  IAB    = IAB + IBITB
                  IF (NIBUF.EQ.1) THEN
                     IABC = IAB
                  ELSE
                     IABC = 0
                  END IF
                  ICMPAB = IADCMP(ICMPA,ICMPB,1)
                  IODDAB = IODDCC(IPNTUV(ICMPAB,ICRB,1))
                  IF (BDER) THEN
                     JODDAB = IODDCC(IPNTUV(ICMPAB,0,1))
                  ELSE
                     JODDAB = 0
                  END IF 
                  IF (CTRIAC) MAXC = ICMPA
                  DO 325 ICMPC = 1, MAXC
                  IVARC = IBTEST(ISTBLC,C,NHKTC,ICMPC)
                  IF (IVARC.EQ.0) THEN
                     IABC  = IABC + IBITC
                     IABCD = IABC
                     IF (CTRIPQ .AND. ICMPA.EQ.ICMPC) THEN
                        MAXD = ICMPB
                     ELSE
                        MAXD = KHKTD
                        IF (CTRICD) MAXD = ICMPC
                     END IF
                     DO 335 ICMPD = 1, MAXD
                     IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                     IF (IVARD.EQ.0) THEN
                        IABCD  = IABCD + 1
                        ICMPCD = IADCMP(ICMPC,ICMPD,2)
                        IODDCD = IODDCC(IPNTUV(ICMPCD,0,2))
                        IF (IODDAB .EQ. IODDCD .OR.
     &                      JODDAB .EQ. IODDCD) THEN
C
C                       standard integrals
C
                        IF (.NOT.BDER) THEN
                           IF (WRTSCR) THEN
                              IF (NIBUF.EQ.1) THEN
                                 DO 406 I = NCCS1 + 1, NCCS
                                    SOABCD = SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    IF (ABS(SOABCD) .GT. THRSH) THEN
                                       INT = INT + 1
                                       BIN (INT  ) = SOABCD
                                       IBIN(INT,1) = INDORB(I,1) + IABCD
                                    END IF
  406                            CONTINUE
                              ELSE
                                 DO 407 I = NCCS1 + 1, NCCS
                                    SOABCD = SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    IF (ABS(SOABCD) .GT. THRSH) THEN
                                       INT = INT + 1
                                       BIN (INT  ) = SOABCD
                                       IBIN(INT,1) = INDORB(I,1) + IAB
                                       IBIN(INT,2) = INDORB(I,2) + IABCD
                                    END IF
  407                            CONTINUE
                              END IF
                           ELSE
                              IF (NIBUF.EQ.1) THEN
                                 DO 456 I = NCCS1 + 1, NCCS
                                    INT = INT + 1
                                    BIN(INT)=SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    IBIN(INT,1) = INDORB(I,1) + IABCD
  456                            CONTINUE
                              ELSE
                                 DO 457 I = NCCS1 + 1, NCCS
                                    INT = INT + 1
                                    BIN(INT)=SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    IBIN(INT,1) = INDORB(I,1) + IAB
                                    IBIN(INT,2) = INDORB(I,2) + IABCD
  457                            CONTINUE
                              END IF
                           END IF
C
C                       B-derivative integrals
C
                        ELSE
                           IF (WRTSCR) THEN
                              IF (NIBUF.EQ.1) THEN
                                 DO 408 I = NCCS1 + 1, NCCS
                                    SO1 = SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    SO2 = SO(I,R,S,T,ICMPAB,ICMPCD,4)
                                    SOAB = SO1 + SO2
                                    SOBA = SO1 - SO2
                                    IJKL = INDORB(I,1) + IABCD
                                    IJLK = INDBEX(IJKL) 
                                    IF (ABS(SOAB) .GT. THRSH) THEN
                                       INT = INT + 1
                                       BIN (INT  ) = SOAB
                                       IBIN(INT,1) = IJKL 
                                    END IF
                                    IF (IJLK.GT.0 .AND. 
     &                                  ABS(SOBA).GT.THRSH) THEN
                                       INT = INT + 1
                                       BIN (INT  ) = SOBA
                                       IBIN(INT,1) = IJLK
                                    END IF
  408                            CONTINUE
                              ELSE
                                 DO 409 I = NCCS1 + 1, NCCS
                                    SO1 = SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    SO2 = SO(I,R,S,T,ICMPAB,ICMPCD,4)
                                    SOAB = SO1 + SO2
                                    SOBA = SO1 - SO2
                                    IF (ABS(SOAB) .GT. THRSH) THEN
                                       INT = INT + 1
                                       BIN (INT  ) = SOABCD
                                       IBIN(INT,1) = INDORB(I,1) + IAB
                                       IBIN(INT,2) = INDORB(I,2) + IABCD
                                    END IF
                                    IF (ABS(SOBA) .GT. THRSH) THEN
                                       INT = INT + 1
                                       BIN (INT  ) = SOBA
                                       IBIN(INT,1) = INDORB(I,1) + IBA
                                       IBIN(INT,2) = INDORB(I,2) + IABCD
                                    END IF
  409                            CONTINUE
                              END IF
                           ELSE
                              IF (NIBUF.EQ.1) THEN
                                 DO 458 I = NCCS1 + 1, NCCS
                                    SO1 = SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    SO2 = SO(I,R,S,T,ICMPAB,ICMPCD,4)
                                    IJKL = INDORB(I,1) + IABCD
                                    IJLK = INDBEX(IJKL) 
                                    INT = INT + 1
                                    BIN(INT) = SO1 + SO2
                                    IBIN(INT,1) = IJKL 
                                    IF (IJLK.GT.0) THEN
                                       INT = INT + 1
                                       BIN(INT) = SO1 - SO2
                                       IBIN(INT,1) = IJLK
                                    END IF
  458                            CONTINUE
                              ELSE
                                 DO 459 I = NCCS1 + 1, NCCS
                                    SO1 = SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    SO2 = SO(I,R,S,T,ICMPAB,ICMPCD,4)
                                    INT = INT + 1
                                    BIN(INT) = SO1 + SO2 
                                    IBIN(INT,1) = INDORB(I,1) + IAB
                                    IBIN(INT,2) = INDORB(I,2) + IABCD
                                    INT = INT + 1
                                    BIN(INT) = SO1 - SO2 
                                    IBIN(INT,1) = INDORB(I,1) + IBA
                                    IBIN(INT,2) = INDORB(I,2) + IABCD
  459                            CONTINUE
                              END IF
                           END IF
                        END IF
                        END IF
                     END IF
  335                CONTINUE
                  END IF
  325             CONTINUE
               END IF
  315          CONTINUE
            END IF
  305       CONTINUE
C
  195       CONTINUE
C
         END IF
  125    CONTINUE
         END IF
  115    CONTINUE
         END IF
  105    CONTINUE
      END IF
C
      IF (BDER) THEN 
         IF (INT.EQ.1) INT = 0
      END IF
C
C     Canonical ordering of indices if requested
C     ==========================================
C
      IF (.NOT.BDER.AND.CANIND) CALL ERICAN(IBIN,INT,LBIN,NIBUF,NBITS)
C
      RETURN
      END
C  /* Deck prprep */
      SUBROUTINE PRPREP(DOREP,NHKTA,KHKTA,ISTBLA)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
      INTEGER A
      LOGICAL DOREP(0:7)
#include <symmet.h>
#include <ibtfun.h>
C
      DO 100 A = 0, MAXREP
         DOREP(A) = .FALSE.
  100 CONTINUE
C
      DO 200 ICMP = 1, KHKTA
         DO 300 A = 0, MAXREP
         IF (IBTAND(ISTBLA,IBTXOR(A,ISYMAO(NHKTA,ICMP))).EQ.0) THEN
            DOREP(A) = .TRUE.
         END IF
  300    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck cmpadr */
      SUBROUTINE CMPADR(IADCMP,KHKTA,KHKTB,TKMPAB)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      LOGICAL TKMPAB
      DIMENSION IADCMP(MXAQN,MXAQN)
C
      DO 100 ICMPA = 1, KHKTA
      DO 100 ICMPB = 1, KHKTB
         IF (TKMPAB) THEN
            MAXAB = MAX(ICMPA,ICMPB)
            IADCMP(ICMPA,ICMPB) = MAXAB*(MAXAB-1)/2 + MIN(ICMPA,ICMPB)
         ELSE
            IADCMP(ICMPA,ICMPB) = (ICMPA - 1)*KHKTB + ICMPB
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck getrst */
      SUBROUTINE GETRST(IPNRST,ISTBLR)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
      INTEGER R, R1, X
      DIMENSION IPNRST(0:7)
#include <symmet.h>
#include <ibtfun.h>
      DO 100 X = 0, MAXREP
         R1 = 0
         DO 200 R = 0, MAXREP
         IF (IBTAND(R,ISTBLR) .EQ. 0) THEN
            R1 = R1 + 1
            IF (IBTOR(X,ISTBLR) .EQ. IBTOR(R,ISTBLR)) THEN
               IPNRST(X) = R1
               GO TO 1000
            END IF
         END IF
  200    CONTINUE
 1000    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck erindd */
      SUBROUTINE ERINDD(SO,INDORB,IPNCHN,INDALL,INDXTR,IPNTCR,IODDCC,
     &                  IPNTUV,BIN,IBIN,NINTS,LBIN,NIBUF,NBITS,JCMPAB,
     &                  JCMPCD,JADDAB,JADDCD,IOFCMP,INDXBT,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
C
#include <eridst.h>
      DIMENSION INDORB(NCCT,KHKTAB,KHKTCD,NIBUF),
     &          SO    (NCCT,KHKTAB,KHKTCD), INDXBT(*),
     &          INDALL(NCCS,2), INDXTR(*),  IPNCHN(*),
     &          IPNTCR(MAXBCH,4), IODDCC(NRTOP), 
     &          IPNTUV(KC2MAX,0:NRDER,2),
     &          BIN(LBIN), IBIN(LBIN,NIBUF), NINTS(NDISTR),
     &          JADDAB(*), JADDCD(*), JCMPAB(*), JCMPCD(*), IOFCMP(*)
#include <cbieri.h>
#include <ericom.h>
#include <erithr.h>
#include <aobtch.h>
#include <hertop.h>
#include <symmet.h>
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERINDD',-1)
C
C     Construct labels
C
      CALL ERILBL(INDORB,INDXTR,INDALL(1,1),INDALL(1,2),
     &            IPNTCR,IODDCC,IPNTUV,NBITS,NIBUF,JCMPAB,JCMPCD,
     &            JADDAB,JADDCD,INDXBT,IPRINT)
C
      NBCHS = MLTPX*KHKTB*KHKTCD
      ICART = 0
      CALL ERIGTH(BIN,IBIN,LBIN,NIBUF,INDORB,SO,
     &            IPNTCR,IODDCC,IPNTUV,IPNCHN,NBITS,NINTS,NBCHS,
     &            IOFCMP,INDXBT,ICART,IPRINT)
      RETURN
      END
C  /* Deck erilbl */
      SUBROUTINE ERILBL(INDORB,INDXTR,INDALL,INDCD,IPNTCR,
     &                  IODDCC,IPNTUV,NBITS,NIBUF,JCMPAB,JCMPCD,
     &                  JADDAB,JADDCD,INDXBT,IPRINT)
C
C     constructs INDORB
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
C
#include <eridst.h>
      INTEGER A, B, C, D, AB, CD, R, S, T
      DIMENSION INDORB(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD,NIBUF),
     &          INDALL(NCCS), INDCD(NCCS), INDXTR(NCCS,MAXREP+1,4),
     &          IPNTCR(MAXBCH,4), IODDCC(NRTOP), 
     &          IPNTUV(KC2MAX,0:NRDER,2), IPNRST(0:7,3), 
     &          IBITX(4), NREP(4), IPNREP(8,4),
     &          JCMPAB(MXAQN**2,0:7,0:7), JCMPCD(MXAQN**2,0:7,0:7),
     &          JADDAB(MXAQN**2,0:7,0:7), JADDCD(MXAQN**2,0:7,0:7),
     &          NCMPAB(0:7,0:7), NCMPCD(0:7,0:7), 
     &          INDXBT(MXSHEL*MXCONT,0:7)
#include <cbieri.h>
#include <ericom.h>
#include <erithr.h>
#include <aobtch.h>
#include <hertop.h>
#include <symmet.h>
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERILBL',-1)
C
      IF (NIBUF .EQ. 1) THEN
         IBITX(1) = 2**(3*NBITS)
         IBITX(2) = 2**(2*NBITS)
         IBITX(3) = 2**(  NBITS)
         IBITX(4) = 1
      ELSE
         IBITX(1) = 2**NBITS
         IBITX(2) = 1
         IBITX(3) = 2**NBITS
         IBITX(4) = 1
      END IF
C
      CALL GETREP(NREP(1),IPNREP(1,1),NHKTA,KHKTA,ISTBLA)
      CALL GETREP(NREP(2),IPNREP(1,2),NHKTB,KHKTB,ISTBLB)
      CALL GETREP(NREP(3),IPNREP(1,3),NHKTC,KHKTC,ISTBLC)
      CALL GETREP(NREP(4),IPNREP(1,4),NHKTD,KHKTD,ISTBLD)
C
      CALL GETRST(IPNRST(0,1),ISTBLR)
      CALL GETRST(IPNRST(0,2),ISTBLS)
      CALL GETRST(IPNRST(0,3),ISTBLT)
C
C     Component loops
C     ===============
C
      CALL GETCMP(JCMPAB,JADDAB,NCMPAB,IBITX(1),IBITX(2),
     &            NHKTA,NHKTB,KHKTA,KHKTB,ISTBLA,ISTBLB,TKMPAB,
     &            NREP(1),IPNREP(1,1),IPRINT)
C
      CALL GETCMP(JCMPCD,JADDCD,NCMPCD,IBITX(3),IBITX(4),
     &            NHKTC,NHKTD,KHKTC,KHKTD,ISTBLC,ISTBLD,TKMPCD,
     &            NREP(3),IPNREP(1,3),IPRINT)
C
C     Extract indices
C     ===============
C
      IRPDST = IPNREP(1,1)
C
      DO 100 J = 1, 4
         NSHFT = IBITX(J)
         DO 110 IA = 1, NREP(J)
            A = IPNREP(IA,J)
            DO 120 I = 1, NCCS
               INDXTR(I,IA,J)=(INDXBT(KNDXBT(IPNTCR(I,J)),A)-1)*NSHFT
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      DO 200 ID = 1, NREP(4)
         D = IPNREP(ID,4)
         S = IPNRST( D,2)
         DO 210 IC = 1, NREP(3)
            C  = IPNREP(IC,3)
            CD = IBTXOR(C,D)
            T  = IPNRST(CD,3)
            DO 215 I = 1, NCCS
               INDCD(I) = INDXTR(I,IC,3) + INDXTR(I,ID,4)
  215       CONTINUE
C
            DO 220 IB = 1, NREP(2)
               B = IPNREP(IB,2)
               R = IPNRST( B,1)
               A = IBTXOR(B,CD)
               DO 230 IA = 1, NREP(1)
               IF (IPNREP(IA,1).EQ.A) THEN
                  IF (NIBUF.EQ.1) THEN
                     DO 301 I = 1, NCCS
                        INDALL(I) = INDXTR(I,IA,1)
     &                            + INDXTR(I,IB,2)+INDCD(I)
  301                CONTINUE
                     DO 401 ICD = 1, NCMPCD(C,D)
                        ICMPCD = JCMPCD(ICD,C,D)
                        IADDCD = JADDCD(ICD,C,D)
                        DO 411 IAB = 1, NCMPAB(A,B)
                           ICMPAB = JCMPAB(IAB,A,B)
                           IABCD  = JADDAB(IAB,A,B) + IADDCD
                           DO 501 I = 1, NCCS
                              INDORB(I,R,S,T,ICMPAB,ICMPCD,1) =
     &                                        INDALL(I) + IABCD
  501                      CONTINUE
  411                   CONTINUE
  401                CONTINUE
                  ELSE
                     DO 302 I = 1, NCCS
                        INDALL(I) = INDXTR(I,IA,1)
     &                            + INDXTR(I,IB,2)
  302                CONTINUE
                     DO 402 ICD = 1, NCMPCD(C,D)
                        ICMPCD = JCMPCD(ICD,C,D)
                        IADDCD = JADDCD(ICD,C,D)
                        DO 412 IAB = 1, NCMPAB(A,B)
                           ICMPAB = JCMPAB(IAB,A,B)
                           IADDAB = JADDAB(IAB,A,B)
                           DO 502 I = 1, NCCS
                              INDORB(I,R,S,T,ICMPAB,ICMPCD,1) =
     &                                       INDALL(I) + IADDAB
                              INDORB(I,R,S,T,ICMPAB,ICMPCD,2) =
     &                                       INDCD(I)  + IADDCD
  502                      CONTINUE
  412                   CONTINUE
  402                CONTINUE
                  END IF
               END IF
  230          CONTINUE
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Output from ERILBL',-1)
         DO 900 ICMPCD = 1, KHKTCD
         DO 900 ICMPAB = 1, KHKTAB
         IF (IODDCC(IPNTUV(ICMPAB,0,1)).EQ.
     &       IODDCC(IPNTUV(ICMPCD,0,2)))THEN
            DO 910 T = 1, MLTPT
            DO 910 S = 1, MLTPS
            DO 910 R = 1, MLTPR
            DO 910 I = 1, NCCS
               IF (NIBUF.EQ.1) THEN
                  IJKL = INDORB(I,R,S,T,ICMPAB,ICMPCD,1)
                  I1 = IBTAND(IBTSHR(IJKL,24),2**NBITS - 1)
                  I2 = IBTAND(IBTSHR(IJKL,16),2**NBITS - 1)
                  I3 = IBTAND(IBTSHR(IJKL,8 ),2**NBITS - 1)
                  I4 = IBTAND(       IJKL,    2**NBITS - 1)
               ELSE
                  IJ = INDORB(I,R,S,T,ICMPAB,ICMPCD,1)
                  KL = INDORB(I,R,S,T,ICMPAB,ICMPCD,2)
                  I1 = IBTAND(IBTSHR(IJ,16),2**NBITS - 1)
                  I2 = IBTAND(       IJ,    2**NBITS - 1)
                  I3 = IBTAND(IBTSHR(KL,16),2**NBITS - 1)
                  I4 = IBTAND(       KL,    2**NBITS - 1)
               END IF
               WRITE (LUPRI,'(1X,A,5X,4I5)') ' INDORB ', I1,I2,I3,I4
  910       CONTINUE
         END IF
  900    CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck getcmp */
      SUBROUTINE GETCMP(JCMPAB,JADDAB,NCMPAB,IBITA,IBITB,
     &                  NHKTA,NHKTB,KHKTA,KHKTB,ISTBLA,ISTBLB,TKMPAB,
     &                  NREP,IPNREP,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
      INTEGER A, B
      LOGICAL TKMPAB
      DIMENSION JCMPAB(MXAQN**2,0:7,0:7),
     &          JADDAB(MXAQN**2,0:7,0:7),
     &          NCMPAB(0:7,0:7), NCMP(8,2), JCMP(MXAQN,8,2),
     &          IADCMP(MXAQN,MXAQN),
     &          NREP(2), IPNREP(8,2)
#include <symmet.h>
#include <hertop.h>
#include <ibtfun.h>
C
      CALL CMPADR(IADCMP,KHKTA,KHKTB,TKMPAB)
C
      DO 10 I = 1, NREP(1)
         A = IPNREP(I,1)
         IKMPA = 0
         DO 20 ICMPA = 1, KHKTA
         IF (IBTAND(ISTBLA,IBTXOR(A,ISYMAO(NHKTA,ICMPA))).EQ.0) THEN
            IKMPA = IKMPA + 1
            JCMP(IKMPA,I,1) = ICMPA
         END IF
   20    CONTINUE
         NCMP(I,1) = IKMPA
   10 CONTINUE
C
      DO 30 I = 1, NREP(2)
         B = IPNREP(I,2)
         IKMPB = 0
         DO 40 ICMPB = 1, KHKTB
         IF (IBTAND(ISTBLB,IBTXOR(B,ISYMAO(NHKTB,ICMPB))).EQ.0) THEN
            IKMPB = IKMPB + 1
            JCMP(IKMPB,I,2) = ICMPB
         END IF
   40    CONTINUE
         NCMP(I,2) = IKMPB
   30 CONTINUE
C
      DO 100 I = 1, NREP(1)
         A = IPNREP(I,1)
         DO 200 J = 1, NREP(2)
            B = IPNREP(J,2)
            IA = 0
            IKMPAB =  0
            DO 300 IKMPA = 1, NCMP(I,1)
               ICMPA = JCMP(IKMPA,I,1)
               IA  = IA + IBITA
               IAB = IA
               DO 400 IKMPB = 1, NCMP(J,2)
                  IAB    = IAB + IBITB
                  IKMPAB = IKMPAB + 1
                  JCMPAB(IKMPAB,A,B) = IADCMP(ICMPA,JCMP(IKMPB,J,2))
                  JADDAB(IKMPAB,A,B) = IAB
  400          CONTINUE
  300       CONTINUE
            NCMPAB(A,B) = IKMPAB
  200    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck getrep */
      SUBROUTINE GETREP(NREP,IPNREP,NHKTA,KHKTA,ISTBLA)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
      INTEGER A
      LOGICAL DOREP(0:7)
      DIMENSION NREP(1), IPNREP(8)
#include <symmet.h>
#include <ibtfun.h>
C
      DO 100 A = 0, MAXREP
         DOREP(A) = .FALSE.
  100 CONTINUE
C
      DO 200 ICMP = 1, KHKTA
         DO 300 A = 0, MAXREP
         IF (IBTAND(ISTBLA,IBTXOR(A,ISYMAO(NHKTA,ICMP))).EQ.0) THEN
            DOREP(A) = .TRUE.
         END IF
  300    CONTINUE
  200 CONTINUE
C
      IA = 0
      DO 400 A = 0, MAXREP
      IF (DOREP(A)) THEN
         IA = IA + 1
         IPNREP(IA) = A
      END IF
  400 CONTINUE
      NREP(1) = IA
C
      RETURN
      END
C  /* Deck erigth */
      SUBROUTINE ERIGTH(BIN,IBIN,LBIN,NIBUF,INDORB,SO,IPNTCR,IODDCC,
     &                  IPNTUV,IPNBCH,NBITS,NINTS,NBCHS,IOFCMP,
     &                  INDXBT,ICART,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
C
#include <eridst.h>
      DIMENSION BIN(LBIN), IBIN(LBIN,NIBUF),
     &          INDORB(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD,NIBUF),
     &          SO    (NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD),
     &          IPNTCR(MAXBCH,4), IODDCC(NRTOP), 
     &          IPNTUV(KC2MAX,0:NRDER,2), IOFCMP(*)
      DIMENSION NCCDST(MXDIST), NCCFST(MXDIST),
     &          NBTCHX(MXDIST), LBTCHX(MXDIST),
     &          IPNBCH(NBCHS,*), NINTS(NDISTR)
      DIMENSION INDXBT(MXSHEL*MXCONT,0:7)
#include <cbieri.h>
#include <ericom.h>
#include <aobtch.h>
#include <hertop.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERIGTH',-1)
C
      IRPDST = 0
C
      IANEW = 0
      IDIST = 0
      DO 100 I = 1, NCCS
         IAOLD = IANEW
         IANEW = INDXBT(KNDXBT(IPNTCR(I,1)),IRPDST)
         IF (IAOLD. NE. IANEW) THEN
            IDIST = IDIST + 1
            NCCDST(IDIST) = 1
            NCCFST(IDIST) = I
         ELSE
            NCCDST(IDIST) = NCCDST(IDIST) + 1
         END IF
 100  CONTINUE
      NDST = IDIST
C
      CALL ERICHN(INDORB,NBTCHX,LBTCHX,IPNBCH,NBCHS,NDST,NCCFST,NCCDST,
     &            IODDCC,IPNTUV,NBITS,NIBUF,IOFCMP,ICART,IPRINT)
C
      ISTR = 1
      DO 500 I = 1, NDISTR
         NBTCH = NBTCHX(I)
         LBTCH = LBTCHX(I)
         CALL ERITRF(BIN(ISTR),IBIN(ISTR,1),NTOT,SO,INDORB,NBTCH,LBTCH,
     &               IPNBCH(1,I),NDISTR,LBIN,NIBUF)
         NINTS(I) = NTOT 
         ISTR = ISTR + NTOT 
  500 CONTINUE
C
      RETURN
      END
C  /* Deck eritrf */
      SUBROUTINE ERITRF(BIN,IBIN,NTOT,SO,INDORB,NBTCH,LBTCH,IPNBCH,
     &                  NDISTR,LBIN,NIBUF)
#include <implicit.h>
#include <priunit.h>
      DIMENSION BIN(LBIN), IBIN(LBIN,NIBUF), IPNBCH(*),
     &          SO(*), INDORB(LBIN,NIBUF)
#include <cbieri.h>
#include <erithr.h>
C
      IF (COMPRS) THEN
         IF (NIBUF.EQ.1) THEN
            INT = 0
            DO I = 1, NBTCH
               IADR = IPNBCH(I)
               DO J = 1, LBTCH
                  SOINT = SO(IADR + J)
                  IF (ABS(SOINT).GT.THRSH) THEN
                     INT = INT + 1
                     BIN (INT)   = SOINT 
                     IBIN(INT,1) = INDORB(IADR + J,1)
                  END IF
               END DO  
            END DO   
         ELSE
            INT = 0
            DO I = 1, NBTCH
               IADR = IPNBCH(I)
               DO J = 1, LBTCH
                  SOINT = SO(IADR + J)
                  IF (ABS(SOINT).GT.THRSH) THEN
                     INT = INT + 1
                     BIN (INT)   = SOINT 
                     IBIN(INT,1) = INDORB(IADR + J,1)
                     IBIN(INT,2) = INDORB(IADR + J,2)
                  END IF
               END DO  
            END DO  
         END IF
         NTOT = INT
      ELSE
         IF (NIBUF.EQ.1) THEN
            INT = 0
            DO I = 1, NBTCH
               IADR = IPNBCH(I)
               DO J = 1, LBTCH
                  BIN (INT + J  ) = SO    (IADR + J)
                  IBIN(INT + J,1) = INDORB(IADR + J,1)
               END DO  
               INT = INT + LBTCH
            END DO   
         ELSE
            INT = 0
            DO I = 1, NBTCH
               IADR = IPNBCH(I)
               DO J = 1, LBTCH
                  BIN (INT + J  ) = SO    (IADR + J)
                  IBIN(INT + J,1) = INDORB(IADR + J,1)
                  IBIN(INT + J,2) = INDORB(IADR + J,2)
               END DO  
               INT = INT + LBTCH
            END DO  
         END IF
         NTOT = NBTCH*LBTCH
      END IF
C
      RETURN
      END
C  /* Deck erichn */
      SUBROUTINE ERICHN(INDORB,NBTCHX,LBTCHX,IPNBCH,NBCHS,NDST,NCCFST,
     &                  NCCDST,IODDCC,IPNTUV,NBITS,NIBUF,IOFCMP,ICART,
     &                  IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <eridst.h>
      DIMENSION NBTCHX(*), LBTCHX(*), IPNBCH(NBCHS,*),
     &          NCCFST(*), NCCDST(*), IODDCC(*), 
     &          IPNTUV(KC2MAX,0:NRDER,2),
     &          IOFCMP(MXAQN**4), INDORB(NCCS,*)
#include <ericom.h>
#include <ibtfun.h>
C
      NCMPX = 0
      MLTAB = MLTPX
      MLTCD = MLTPX*KHKTAB
      DO 100 ICMPCD = 1, KHKTCD
         IODDCD = IODDCC(IPNTUV(ICMPCD,ICART,2))
         DO 200 ICMPAB = 1, KHKTAB
            IODDAB = IODDCC(IPNTUV(ICMPAB,0,1))
            IF (IODDAB .EQ. IODDCD) THEN
               NCMPX = NCMPX + 1
               IOFCMP(NCMPX) = (ICMPCD - 1)*MLTCD + (ICMPAB - 1)*MLTAB
            END IF
  200    CONTINUE
  100 CONTINUE
C
      IBT = 2**NBITS - 1
      IF (NIBUF.EQ.1) THEN
         NB3 = 3*NBITS
      ELSE
         NB3 =   NBITS
      END IF
      CALL IZERO(NBTCHX,NDISTR)
      DO 300 I = 1, NDST
         ISTR   = NCCFST(I)
         LENGTH = NCCDST(I)
         DO 400 ICMPX = 1, NCMPX
            IOFF = IOFCMP(ICMPX)
            IADR = IOFF*NCCS + ISTR - 1
            DO 500 K = IOFF + 1, IOFF + MLTPX
               IDIST = INDXDS(IBTAND(IBTSHR(INDORB(ISTR,K),NB3),IBT))
               IBTCH = NBTCHX(IDIST) + 1
               NBTCHX(IDIST) = IBTCH
               LBTCHX(IDIST) = LENGTH
               IPNBCH(IBTCH,IDIST) = IADR
               IADR = IADR + NCCS
  500       CONTINUE
  400    CONTINUE
  300 CONTINUE
      RETURN
      END
C  /* Deck indswp */
      FUNCTION INDSWP(IJKL)
#include <implicit.h>
#include <ibtpar.h>
#include <ibtfun.h>
      KL = IBTAND(IJKL,65535)
      INDSWP = IJKL - KL + 256*IBTAND(KL,255) + IBTAND(IBTSHR(KL,8),255)
      RETURN
      END 
C  /* Deck indbex */
      FUNCTION INDBEX(INDEX)
#include <implicit.h>
#include <priunit.h>
#include <ibtpar.h>
      PARAMETER (NBIT1 = 8,        NBIT2 = 2*NBIT1,  NBIT3 = 3*NBIT1)
      PARAMETER (IBIT1 = 2**NBIT1, IBIT2 = 2**NBIT2, IBIT3 = 2**NBIT3)
      PARAMETER (IBITS = IBIT1 - 1)
#include <ibtfun.h>
C
      I = IBTAND(IBTSHR(INDEX,NBIT3),IBITS)
      J = IBTAND(IBTSHR(INDEX,NBIT2),IBITS)
      K = IBTAND(IBTSHR(INDEX,NBIT1),IBITS)
      L = IBTAND(       INDEX,       IBITS)
      IF (L.EQ.0) THEN
         INDBEX = 0
      ELSE IF (I.NE.J .AND. K.NE.L) THEN
         INDBEX = I*IBIT3 + J*IBIT2 + L*IBIT1 + K
      ELSE
         INDBEX = 0
      END IF
      RETURN
      END
C  /* Deck eriprt */
      SUBROUTINE ERIPRT(SO,ISCOOR,IATOM,ICART,SOPRT,LABEL,INDORB,IPRINT)
C
C     T. Helgaker
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
      PARAMETER (D0 = 0.0D0)
#include <aobtch.h>
#include <symmet.h>
      DIMENSION SOPRT(NAOINT), SO(NAOINT,3,4), LABEL(NAOINT,4),
     &          INDORB(*)
#include <ericom.h>
#include <eribuf.h>
#include <erithr.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERIPRT',-1)
C
      DO I = 1, NAOINT
         GINT = D0 
         IF (LABEL(I,1) .EQ. IATOM) GINT = GINT + SO(I,ICART,1) 
         IF (LABEL(I,2) .EQ. IATOM) GINT = GINT + SO(I,ICART,2) 
         IF (LABEL(I,3) .EQ. IATOM) GINT = GINT + SO(I,ICART,3) 
         IF (LABEL(I,4) .EQ. IATOM) GINT = GINT + SO(I,ICART,4) 
         SOPRT(I) = GINT
      END DO
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Output from ERIPRT',-1)
         CALL ERIPRI(SOPRT,INDORB,NAOINT,NIBUF,NAOINT,NBITS,0,IDUMMY,
     &               ISCOOR)
      END IF
C
      RETURN
      END
      SUBROUTINE ERIBRT(SO,ISCOOR,SOPRT,IPRINT)
C
C     T. Helgaker
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aobtch.h>
#include <symmet.h>
      DIMENSION SOPRT(NAOINT), SO(NAOINT,6)
#include <ericom.h>
#include <eribuf.h>
#include <erithr.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERIPRT',-1)
C
      IF (ISCOOR .LE. 3) THEN
         DO I = 1, NAOINT
            SOPRT(I) = SO(I,ISCOOR) + SO(I,ISCOOR+3) 
         END DO
      ELSE
         DO I = 1, NAOINT
            SOPRT(I) = SO(I,ISCOOR-3) - SO(I,ISCOOR) 
         END DO
      END IF
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Output from ERIPRT',-1)
         CALL ERIPRI(SOPRT,INDORB,NAOINT,NIBUF,NAOINT,NBITS,0,IDUMMY,
     &               ISCOOR)
      END IF
C
      RETURN
      END
C  /* Deck eriave */
      SUBROUTINE ERIAVE(SO,PMAT,DMAT,D2MAT,ID2MAT,IPNTCR,IODDCC,IPNTUV,
     &                  INDXBT,WORK,LWORK,IPRINT)
C
C     tuh march 99
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
C
      DIMENSION SO(*), PMAT(*), DMAT(*), D2MAT(*), ID2MAT(*),
     &          IPNTCR(*), IPNTUV(*), IODDCC(*),
     &          INDXBT(*), WORK(*)
#include <ericom.h>
C
C     Allocations
C
      LBIN   = NCCT*KHKTA*KHKTB*KHKTC*KHKTD
      KWGHT  = 1
      KLABEL = KWGHT + NAOINT
      KPOINT = KLABEL + (4*NAOINT - 1)/IRAT + 1
      KLAST  = KPOINT + (4*LBIN - 1)/IRAT + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('ERIAVE',' ',KLAST,LWORK)
      LWRK = LWORK - KLAST + 1
C
C     Generate weights and labels
C     ===========================
C 
      CALL ERILAB(WORK(KWGHT),WORK(KLABEL),IDUMMY,IPNTCR,IODDCC,IPNTUV,
     &            INDXBT,WORK(KPOINT),IPRINT)
C
C     Generate density matrix 
C     =======================
C
      CALL ERIDMT(PMAT,WORK(KWGHT),WORK(KLABEL),IPNTUV,IODDCC,
     &            DMAT,D2MAT,ID2MAT,IPRINT)
C
C     Multiply integrals and densities 
C     ================================
C
      CALL ERIGRD(SO,PMAT,WORK(KLABEL),IPNTUV,IODDCC,IPRINT)
C
      RETURN
      END
C  /* Deck erilab */
      SUBROUTINE ERILAB(WEIGHT,LABEL,INDORB,IPNTCR,IODDCC,IPNTUV,INDXBT,
     &                  IPOINT,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
      PARAMETER (D8=8.0D0, D4=4.0D0, D2=2.0D0, D1=1.0D0, DP5=0.5D0)
      INTEGER A, B, C, D, AB, CD, R, S, T
      LOGICAL DOREP(0:7,4), CTRIAB, CTRICD
      DIMENSION IPOINT(NCCS,4), 
     &          IPNTCR(MAXBCH,4), INDXBT(MXSHEL*MXCONT,0:7),
     &          IODDCC(NRTOP), IPNTUV(KC2MAX,0:NRDER,2),
     &          IPNRST(0:7,3),
     &          IADCMP(MXAQN,MXAQN,2),
     &          WEIGHT(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD),
     &          LABEL(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD,4),
     &          INDORB(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD,NIBUF)
#include <cbieri.h>
#include <ericom.h>
#include <erithr.h>
#include <eribuf.h>
#include <aobtch.h>
#include <hertop.h>
#include <symmet.h>
C
#include <ibtfun.h>
      IBTEST(I,J,K,L) = IBTAND(I,IBTXOR(J,ISYMAO(K,L)))
      IBTREP(I,J,K,L) = IBTXOR(J,IBTAND(I,ISYMAO(K,L)))
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERILA1',-1)
C
      IF (IPRINT .GT. 15) THEN
         WRITE (LUPRI,'(2X,A,2I5)') ' NPQBCS,NPPBCS ',NPQBCS,NPPBCS
      END IF
C
      CALL PRPREP(DOREP(0,1),NHKTA,KHKTA,ISTBLA)
      CALL PRPREP(DOREP(0,2),NHKTB,KHKTB,ISTBLB)
      CALL PRPREP(DOREP(0,3),NHKTC,KHKTC,ISTBLC)
      CALL PRPREP(DOREP(0,4),NHKTD,KHKTD,ISTBLD)
C
      CALL CMPADR(IADCMP(1,1,1),KHKTA,KHKTB,TKMPAB)
      CALL CMPADR(IADCMP(1,1,2),KHKTC,KHKTD,TKMPCD)
C
      CALL GETRST(IPNRST(0,1),ISTBLR)
      CALL GETRST(IPNRST(0,2),ISTBLS)
      CALL GETRST(IPNRST(0,3),ISTBLT)
C
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(/,2X,A,8L2)')'DOREP A  ',(DOREP(I,1),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP B  ',(DOREP(I,2),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP C  ',(DOREP(I,3),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP D  ',(DOREP(I,4),I=0,MAXREP)
      END IF
C
      NCCS1 = (NPQBCS - NPPBCS)*NCTFAB*NCTFCD
      NCCS2 = NPPBCS*NCTFAB*NCTFCD
C
      IF (EXPERI) CALL DZERO(WEIGHT,NAOINT)
C
      IBIT1 = 2**NBITS     - 1
      IBIT2 = 2**(2*NBITS) - 1
C
      IF (NCCS1 .GT. 0) THEN
C
         DO 100 A = 0, MAXREP
         IF (DOREP(A,1)) THEN
            DO 101 I = 1, NCCS1
               IPOINT(I,1) = INDXBT(KNDXBT(IPNTCR(I,1)),A)
  101       CONTINUE
         DO 110 B = 0, MAXREP
         IF (DOREP(B,2)) THEN
            DO 111 I = 1, NCCS1
               IPOINT(I,2) = INDXBT(KNDXBT(IPNTCR(I,2)),B)
  111       CONTINUE
         DO 120 C = 0, MAXREP
         IF (DOREP(C,3) .AND. DOREP(IBTXOR(IBTXOR(A,B),C),4)) THEN
            D = IBTXOR(IBTXOR(A,B),C)
            CD = IBTXOR(C,D)
C
            IF (DIAGAB .AND. B.GT.A) GO TO 190
            IF (DIAGCD .AND. D.GT.C) GO TO 190
C
            DO 121 I = 1, NCCS1
               IPOINT(I,3) = INDXBT(KNDXBT(IPNTCR(I,3)),C)
               IPOINT(I,4) = INDXBT(KNDXBT(IPNTCR(I,4)),D)
  121       CONTINUE
C
            R = IPNRST(B,1)
            S = IPNRST(D,2)
            T = IPNRST(CD,3)
C
            IA   = -1
            MAXB = KHKTB
            MAXD = KHKTD
            CTRIAB = DIAGAB .AND. A.EQ.B .AND. EXPERI 
            CTRICD = DIAGCD .AND. C.EQ.D .AND. EXPERI 
            DO 300 ICMPA = 1, KHKTA
            IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
            IF (IVARA.EQ.0) THEN
               IA = IA + 1
               IB = -1
               IF (CTRIAB) MAXB = ICMPA
               DO 310 ICMPB = 1, MAXB
               IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
               IF (IVARB.EQ.0) THEN
                  IB = IB + 1
                  IC = -1
                  ICMPAB = IADCMP(ICMPA,ICMPB,1)
                  IODDAB = IODDCC(IPNTUV(ICMPAB,0,1))
                  DO 320 ICMPC = 1, KHKTC
                  IVARC = IBTEST(ISTBLC,C,NHKTC,ICMPC)
                  IF (IVARC.EQ.0) THEN
                     IC = IC + 1
                     ID = -1
                     IF (CTRICD) MAXD = ICMPC
                     DO 330 ICMPD = 1, MAXD
                     IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                     IF (IVARD.EQ.0) THEN
                        ID = ID + 1
                        ICMPCD = IADCMP(ICMPC,ICMPD,2)
                        IF (EXPERI) THEN
                           IODDCX = IODDCC(IPNTUV(ICMPCD,1,2))
                           IODDCY = IODDCC(IPNTUV(ICMPCD,2,2))
                           IODDCZ = IODDCC(IPNTUV(ICMPCD,3,2))
                           IF (IODDAB.EQ.IODDCX .OR.  
     &                         IODDAB.EQ.IODDCY .OR.  
     &                         IODDAB.EQ.IODDCZ) THEN
                            DO 400 I = 1, NCCS1
                               KA = IPOINT(I,1)+IA
                               KB = IPOINT(I,2)+IB
                               KC = IPOINT(I,3)+IC
                               KD = IPOINT(I,4)+ID
                               IF (CCRUN) THEN
                                  FAC = D2
                                  IF (KC.EQ.KD) FAC = DP5*FAC
                               ELSE
                                  FAC = D8
                                  IF (KA.EQ.KB) FAC = DP5*FAC
                                  IF (KC.EQ.KD) FAC = DP5*FAC
                               END IF
                               WEIGHT(I,R,S,T,ICMPAB,ICMPCD) = FAC
                               LABEL(I,R,S,T,ICMPAB,ICMPCD,1) = KA
                               LABEL(I,R,S,T,ICMPAB,ICMPCD,2) = KB
                               LABEL(I,R,S,T,ICMPAB,ICMPCD,3) = KC
                               LABEL(I,R,S,T,ICMPAB,ICMPCD,4) = KD
  400                        CONTINUE
                          END IF
                        ELSE
                           IODDC0 = IODDCC(IPNTUV(ICMPCD,0,2))
                           IODDCX = IODDCC(IPNTUV(ICMPCD,1,2))
                           IODDCY = IODDCC(IPNTUV(ICMPCD,2,2))
                           IODDCZ = IODDCC(IPNTUV(ICMPCD,3,2))
                           IF (IODDAB.EQ.IODDC0 .OR.  
     &                         IODDAB.EQ.IODDCX .OR.  
     &                         IODDAB.EQ.IODDCY .OR.  
     &                         IODDAB.EQ.IODDCZ) THEN
                             DO 401 I = 1, NCCS1
                               KA = IPOINT(I,1)+IA
                               KB = IPOINT(I,2)+IB
                               KC = IPOINT(I,3)+IC
                               KD = IPOINT(I,4)+ID
                               KAB = KA*IBIT1 + KA + KB
                               KCD = KC*IBIT1 + KC + KD
                               IF (NIBUF.EQ.1) THEN
                                  INDORB(I,R,S,T,ICMPAB,ICMPCD,1) =
     &                                   KAB*IBIT2+KAB+KCD
                               ELSE
                                  INDORB(I,R,S,T,ICMPAB,ICMPCD,1) =
     &                                                KAB
                                  INDORB(I,R,S,T,ICMPAB,ICMPCD,2) =
     &                                                  KCD
                               END IF
                               LABEL(I,R,S,T,ICMPAB,ICMPCD,1)=ICNTAO(KA)
                               LABEL(I,R,S,T,ICMPAB,ICMPCD,2)=ICNTAO(KB)
                               LABEL(I,R,S,T,ICMPAB,ICMPCD,3)=ICNTAO(KC)
                               LABEL(I,R,S,T,ICMPAB,ICMPCD,4)=ICNTAO(KD)
  401                        CONTINUE
                           END IF
                        END IF
                     END IF
  330                CONTINUE
                  END IF
  320             CONTINUE
               END IF
  310          CONTINUE
            END IF
  300       CONTINUE
C
  190       CONTINUE
C
         END IF
  120    CONTINUE
         END IF
  110    CONTINUE
         END IF
  100    CONTINUE
      END IF
C
      IF (NCCS2 .GT. 0) THEN
C
         DO 105 A = 0, MAXREP
         IF (DOREP(A,1)) THEN
            DO 106 I = NCCS1 + 1, NCCS
               IPOINT(I,1) = INDXBT(KNDXBT(IPNTCR(I,1)),A)
  106       CONTINUE
         DO 115 B = 0, MAXREP
         IF (DOREP(B,2)) THEN
            DO 116 I = NCCS1 + 1, NCCS
               IPOINT(I,2) = INDXBT(KNDXBT(IPNTCR(I,2)),B)
  116       CONTINUE
         DO 125 C = 0, MAXREP
         IF (DOREP(C,3) .AND. DOREP(IBTXOR(IBTXOR(A,B),C),4)) THEN
            D = IBTXOR(IBTXOR(A,B),C)
            CD = IBTXOR(C,D)
C
            IF (DIAGAB .AND. B.GT.A) GO TO 195
            IF (DIAGCD .AND. D.GT.C) GO TO 195
C           IF (C.GT.A .OR. (C.EQ.A .AND. D.GT.B)) GO TO 195
C
            DO 126 I = NCCS1 + 1, NCCS
               IPOINT(I,3) = INDXBT(KNDXBT(IPNTCR(I,3)),C)
               IPOINT(I,4) = INDXBT(KNDXBT(IPNTCR(I,4)),D)
  126       CONTINUE
C
            R = IPNRST(B,1)
            S = IPNRST(D,2)
            T = IPNRST(CD,3)
C
            MAXB = KHKTB
            MAXD = KHKTD
            CTRIAB = DIAGAB .AND. A.EQ.B .AND. EXPERI 
            CTRICD = DIAGCD .AND. C.EQ.D .AND. EXPERI 
C
            IA = -1
            DO 305 ICMPA = 1, KHKTA
            IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
            IF (IVARA.EQ.0) THEN
               IA = IA + 1
               IB = -1
C              IF (CTRIAB) MAXB = ICMPA
               DO 315 ICMPB = 1, MAXB
               IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
               IF (IVARB.EQ.0) THEN
                  IB = IB + 1
                  IC = -1
                  ICMPAB = IADCMP(ICMPA,ICMPB,1)
                  IODDAB = IODDCC(IPNTUV(ICMPAB,0,1))
                  DO 325 ICMPC = 1, KHKTC
                  IVARC = IBTEST(ISTBLC,C,NHKTC,ICMPC)
                  IF (IVARC.EQ.0) THEN
                     IC = IC + 1
                     ID = -1
C                    IF (CTRICD) MAXD = ICMPC
                     DO 335 ICMPD = 1, MAXD
                     IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                     IF (IVARD.EQ.0) THEN
                        ID = ID + 1
                        ICMPCD = IADCMP(ICMPC,ICMPD,2)
                        IF (EXPERI) THEN
                           IODDCX = IODDCC(IPNTUV(ICMPCD,1,2))
                           IODDCY = IODDCC(IPNTUV(ICMPCD,2,2))
                           IODDCZ = IODDCC(IPNTUV(ICMPCD,3,2))
                           IF (IODDAB.EQ.IODDCX .OR.  
     &                         IODDAB.EQ.IODDCY .OR.  
     &                         IODDAB.EQ.IODDCZ) THEN
                             DO 405 I = NCCS1 + 1, NCCS
                               KA = IPOINT(I,1)+IA
                               KB = IPOINT(I,2)+IB
                               KC = IPOINT(I,3)+IC
                               KD = IPOINT(I,4)+ID
                               IF (CCRUN) THEN
                                  FAC = D1
                                  IF (KC.EQ.KD) FAC = DP5*FAC
                               ELSE
                                  FAC = D4
                                  IF (KA.EQ.KB) FAC = DP5*FAC
                                  IF (KC.EQ.KD) FAC = DP5*FAC
                               END IF
                               WEIGHT(I,R,S,T,ICMPAB,ICMPCD) = FAC
                               LABEL(I,R,S,T,ICMPAB,ICMPCD,1) = KA
                               LABEL(I,R,S,T,ICMPAB,ICMPCD,2) = KB
                               LABEL(I,R,S,T,ICMPAB,ICMPCD,3) = KC
                               LABEL(I,R,S,T,ICMPAB,ICMPCD,4) = KD
  405                        CONTINUE
                           END IF
                        ELSE
                           IODDC0 = IODDCC(IPNTUV(ICMPCD,0,2))
                           IODDCX = IODDCC(IPNTUV(ICMPCD,1,2))
                           IODDCY = IODDCC(IPNTUV(ICMPCD,2,2))
                           IODDCZ = IODDCC(IPNTUV(ICMPCD,3,2))
                           IF (IODDAB.EQ.IODDC0 .OR.  
     &                         IODDAB.EQ.IODDCX .OR.  
     &                         IODDAB.EQ.IODDCY .OR.  
     &                         IODDAB.EQ.IODDCZ) THEN
                             DO 406 I = NCCS1 + 1, NCCS
                               KA = IPOINT(I,1)+IA
                               KB = IPOINT(I,2)+IB
                               KC = IPOINT(I,3)+IC
                               KD = IPOINT(I,4)+ID
                               KAB = KA*IBIT1 + KA + KB
                               KCD = KC*IBIT1 + KC + KD
                               IF (NIBUF.EQ.1) THEN
                                  INDORB(I,R,S,T,ICMPAB,ICMPCD,1) =
     &                                   KAB*IBIT2+KAB+KCD
                               ELSE
                                  INDORB(I,R,S,T,ICMPAB,ICMPCD,1) =
     &                                                  KAB
                                  INDORB(I,R,S,T,ICMPAB,ICMPCD,2) =
     &                                                  KCD
                               END IF
                               LABEL(I,R,S,T,ICMPAB,ICMPCD,1)=ICNTAO(KA)
                               LABEL(I,R,S,T,ICMPAB,ICMPCD,2)=ICNTAO(KB)
                               LABEL(I,R,S,T,ICMPAB,ICMPCD,3)=ICNTAO(KC)
                               LABEL(I,R,S,T,ICMPAB,ICMPCD,4)=ICNTAO(KD)
  406                        CONTINUE
                           END IF
                        END IF
                     END IF
  335                CONTINUE
                  END IF
  325             CONTINUE
               END IF
  315          CONTINUE
            END IF
  305       CONTINUE
C
  195       CONTINUE
C
         END IF
  125    CONTINUE
         END IF
  115    CONTINUE
         END IF
  105    CONTINUE
      END IF
      RETURN
      END
C  /* Deck eridmt */
      SUBROUTINE ERIDMT(PMAT,WEIGHT,LABEL,IPNTUV,IODDCC,DMAT,
     &                  D2MAT,ID2MAT,IPRINT)
C
C     T. Helgaker
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
      PARAMETER (DP5 = 0.5D0, DP25 = 0.25D0)
      LOGICAL FROMD1
      INTEGER A, B, C, D
#include <aobtch.h>
#include <symmet.h>
      DIMENSION PMAT(NAOINT), WEIGHT(NAOINT), LABEL(NAOINT,4),
     &          DMAT(NBASE,NBASE), 
     &          D2MAT(NBASE,NBASE,NPDIMB,NPDIMA),
     &          ID2MAT(MXCORB,2),
     &          IPNTUV(KC2MAX,0:NRDER,2), IODDCC(NRTOP)
#include <ericom.h>
#include <hertop.h>
#include <ibtfun.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERIDMT',-1)
C
      FROMD1 = NPDIMA.EQ.0 .AND. NPDIMB.EQ.0
      IF (.false.) THEN
         DO 100 I = 1, NAOINT
            A = LABEL(I,1) 
            B = LABEL(I,2) 
            C = LABEL(I,3) 
            D = LABEL(I,4) 
            LABEL(I,1) = ICNTAO(A)
            LABEL(I,2) = ICNTAO(B)
            LABEL(I,3) = ICNTAO(C)
            LABEL(I,4) = ICNTAO(D)
            PMAT(I) = DP5*WEIGHT(I)*(DMAT(A,B)*DMAT(C,D) 
     &                        - DP25*DMAT(A,C)*DMAT(B,D)
     &                        - DP25*DMAT(A,D)*DMAT(B,C))
  100    CONTINUE
      ELSE
         IOFF = 0
         DO 400 K = 1, KHKTCD
            IODDCX = IODDCC(IPNTUV(K,1,2))
            IODDCY = IODDCC(IPNTUV(K,2,2))
            IODDCZ = IODDCC(IPNTUV(K,3,2))
            DO 410 J = 1, KHKTAB
               IODDAB = IODDCC(IPNTUV(J,0,1))
               IF (IODDAB.EQ.IODDCX .OR. IODDAB.EQ.IODDCY .OR.  
     &                                   IODDAB.EQ.IODDCZ) THEN
                  IF (FROMD1) THEN
                     DO 500 I = IOFF + 1, IOFF + NCCT
                        A = LABEL(I,1) 
                        B = LABEL(I,2) 
                        C = LABEL(I,3) 
                        D = LABEL(I,4) 
                        LABEL(I,1) = ICNTAO(A)
                        LABEL(I,2) = ICNTAO(B)
                        LABEL(I,3) = ICNTAO(C)
                        LABEL(I,4) = ICNTAO(D)
                        PMAT(I) = DP5*WEIGHT(I)*(DMAT(A,B)*DMAT(C,D) 
     &                                    - DP25*DMAT(A,C)*DMAT(B,D)
     &                                    - DP25*DMAT(A,D)*DMAT(B,C))
  500                CONTINUE
                   ELSE
                     DO 510 I = IOFF + 1, IOFF + NCCT
                        A = LABEL(I,1) 
                        B = LABEL(I,2) 
                        C = LABEL(I,3) 
                        D = LABEL(I,4) 
                        LABEL(I,1) = ICNTAO(A)
                        LABEL(I,2) = ICNTAO(B)
                        LABEL(I,3) = ICNTAO(C)
                        LABEL(I,4) = ICNTAO(D)
                        PMAT(I) = DP5*WEIGHT(I)*
     &                            D2MAT(D,C,ID2MAT(B,2),ID2MAT(A,1))
  510                CONTINUE
                   END IF
               END IF
               IOFF = IOFF + NCCT
  410       CONTINUE
  400    CONTINUE
      END IF
      IF (IPRINT .GT. 30) THEN
         CALL HEADER('DMAT in ERIDMT',-1)
         CALL OUTPUT(DMAT,1,NBASE,1,NBASE,NBASE,NBASE,1,LUPRI)
         CALL HEADER('PMAT in ERIDMT',-1)
         CALL OUTPUT(PMAT,1,1,1,NAOINT,1,NAOINT,1,LUPRI)
      END IF
      RETURN
      END 
C  /* Deck erigrd */
      SUBROUTINE ERIGRD(SO,PMAT,LABEL,IPNTUV,IODDCC,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <mxorb.h>
      INTEGER A, X
      DIMENSION SO(NAOINT,3,4), PMAT(NAOINT), LABEL(NAOINT,4),
     &          IPNTUV(KC2MAX,0:NRDER,2), IODDCC(NRTOP)
#include <nuclei.h>
#include <energy.h>
#include <cbieri.h>
#include <ericom.h>
#include <hertop.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERIGRD',-1)
C
      IF (GRDZER) CALL DZERO(GRADEE,MXCOOR)
C
      IF (GDER) THEN
         IF (.false.) THEN
            DO 100 A = 1, 4
               DO 200 I = 1, NAOINT
                  IAX = 3*LABEL(I,A) - 2
                  IAY = IAX + 1
                  IAZ = IAX + 2
                  GRADEE(IAX) = GRADEE(IAX) + PMAT(I)*SO(I,1,A)
                  GRADEE(IAY) = GRADEE(IAY) + PMAT(I)*SO(I,2,A)
                  GRADEE(IAZ) = GRADEE(IAZ) + PMAT(I)*SO(I,3,A)
  200          CONTINUE
  100       CONTINUE
         ELSE
            DO 300 A = 1, 4
            DO 300 X = 1, 3
               IOFF = 0
               DO 400 K = 1, KHKTCD
                  IODDCD = IODDCC(IPNTUV(K,X,2))
                  DO 410 J = 1, KHKTAB
                     IF (IODDCC(IPNTUV(J,0,1)).EQ.IODDCD) THEN
                        DO 500 I = IOFF + 1, IOFF + NCCT
                           IA = 3*LABEL(I,A) - 3 + X
                           GRADEE(IA) = GRADEE(IA) + PMAT(I)*SO(I,X,A)
  500                   CONTINUE
                     END IF
                     IOFF = IOFF + NCCT
  410             CONTINUE
  400          CONTINUE
  300       CONTINUE
         END IF
      ELSE
         DO 600 J = 1, 3
         DO 600 I = 1, NAOINT
            GRADEE(J) = GRADEE(J) + PMAT(I)*SO(I,J,1)
  600    CONTINUE
      END IF
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('GRADEE in ERIGRD',-1) 
         CALL OUTPUT(GRADEE,1,3,1,NUCDEP,3,NUCDEP,1,LUPRI)
      END IF
      RETURN
      END
