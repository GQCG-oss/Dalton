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
C$Id: eri2out.F,v 1.7 2001-04-18 07:09:09 helgaker Exp $
C
#include <single.h>
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
      SUBROUTINE ERIOU1(SO,SOPRT,ITYPE,INDORB,LABEL,IPNBCH,INDALL,
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
     &          INDORB(LBIN,NIBUF), IPNBCH(*), INDALL(NCCS,2), 
     &          INDXTR(*), LABEL(*), IPNTCR(MAXBCH,4),
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
#include <symmet.h>
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
         CALL ERILAB(DUMMY,LABEL,INDORB,IPNTCR,IODDCC,IPNTUV,INDXBT,
     &               INDXTR,0,0,IPRINT)
C
C        loop over Cartesian directions
C
         ISTART = 1
         IF (UNDIFF) ISTART = 0
         DO ISCOOR = ISTART, NPERTS 
            IF (IPRINT.GT.10) THEN
               WRITE (LUPRI,'(2X,A,I5)') ' ISCOOR in ERIOU1 ', ISCOOR
            END IF
            IATOM = (ISCOOR + 2)/3
            ICART = MOD(ISCOOR - 1,3) + 1
C
C           construct derivatives
C
            IF (ISCOOR .GT. 0) THEN
               IF (GDER) THEN
                  CALL ERIPRT(SO,ISCOOR,IATOM,ICART,SOPRT,LABEL,INDORB,
     &                        IPRINT)
               ELSE
                  CALL ERILAB(DUMMY,LABEL,INDORB,IPNTCR,IODDCC,IPNTUV,
     &                        INDXBT,INDXTR,ISYMAX(ICART,2),1,IPRINT)
                  CALL ERIBRT(SO,ISCOOR,SOPRT,INDORB,IPRINT)
               END IF
            END IF
C
C           sort distributions
C
            CALL ERIGTH(BIN,IBIN,LBIN,NIBUF,INDORB,SOPRT,IPNTCR,IODDCC,
     &                  IPNTUV,IPNBCH,NBITS,NINTS,IOFCMP,INDXBT,
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
C        construct labels
C
         CALL ERILBL(INDORB,INDXTR,INDALL(1,1),INDALL(1,2),
     &               IPNTCR,IODDCC,IPNTUV,NBITS,NIBUF,JCMPAB,JCMPCD,
     &               JADDAB,JADDCD,INDXBT,IPRINT)
C
C        sort distributions
C
         CALL ERIGTH(BIN,IBIN,LBIN,NIBUF,INDORB,SO,
     &               IPNTCR,IODDCC,IPNTUV,IPNBCH,NBITS,NINTS,
     &               IOFCMP,INDXBT,0,IPRINT)
C
C        write integrals
C
         CALL ERIWRT(0,NINTS,NCOUNT,BIN,IBIN,LBIN,BUF(1,1,0),
     &               IBUF(1,1,1,0),IPRINT)
C
C     Undifferentiated integrals (all)
C     ================================
C
      ELSE
C
C        collect integrals and labels
C
         CALL ERINDX(SO,ITYPE,INDORB,IPNTCR,IODDCC,IPNTUV,BIN,IBIN,
     &               NINTS,LBIN,NIBUF,INDXBT,NBITS,IPRINT)
C
C        write integrals
C
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
      CHARACTER*8 FAODER
      DIMENSION BUF(LBUF,NBUFS), IBUF(LBUF,NIBUF,NBUFS),
     &          BIN(LBIN), IBIN(LBIN,NIBUF),
     &          NINTS(NDISTR), NCOUNT(NDISTR,0:NCORS)
#include <cbieri.h>
#include <chrnos.h>
#include <symmet.h>
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
         ELSE
            LUINTX = LUINTA
            CALL REWSPL(LUINTX)
            CALL NEWLAB('BASINFO ',LUINTA,LUPRI)
            WRITE (LUINTA) MAXREP+1,(NAOS(I),I=1,8),LBUF,NIBUF,NBITS
            CALL NEWLAB('BASTWOEL',LUINTX,LUPRI)
         END IF
         CALL IZERO(NCOUNT(1,0),NDISTR*NCORS)
         CALL IZERO(NBUFX(0),NCORS)
         FIRST = .FALSE.
      END IF
C
      IF (GDER .OR. BDER) THEN
         LUINTX = LUINTD(ISCOOR)
      ELSE
         LUINTX = LUINTA
      END IF
C
C     Run over distributions
C     ======================
C
      IF (.NOT.LAST) THEN
         IOFF = 0
         DO IDIST = 1, NDISTR
            ICOUNT = NCOUNT(IDIST,ISCOOR)
            NTOT   = NINTS(IDIST)
C
C           Empty whole buffers directly from BIN and IBIN
C           ==============================================
C
            KBUFS = NTOT/LBUF
            DO I = 1, KBUFS
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
            END DO
C
C           Case A: Transfer remainder to buffer
C           ====================================
C
            ISTR = KBUFS*LBUF + 1
            IF (ICOUNT + (NTOT - ISTR + 1) .LT. LBUF) THEN
               IF (NIBUF .EQ. 1) THEN
                  DO I = IOFF + ISTR, IOFF + NTOT
                     ICOUNT = ICOUNT + 1
                     BUF (ICOUNT,  IDIST) =  BIN(I)
                     IBUF(ICOUNT,1,IDIST) = IBIN(I,1)
                  END DO
               ELSE
                  DO I = IOFF + ISTR, IOFF + NTOT
                     ICOUNT = ICOUNT + 1
                     BUF (ICOUNT,  IDIST) =  BIN(I)
                     IBUF(ICOUNT,1,IDIST) = IBIN(I,1)
                     IBUF(ICOUNT,2,IDIST) = IBIN(I,2)
                  END DO
               END IF
C
C           Case B: Transfer remainder to buffer and empty when full
C           ========================================================
C
            ELSE
               NLEFT = LBUF - ICOUNT
               ILAST = ISTR + NLEFT - 1
               IF (NIBUF .EQ. 1) THEN
                  DO I = IOFF + ISTR, IOFF + ILAST
                     ICOUNT = ICOUNT + 1
                     BUF (ICOUNT,  IDIST) =  BIN(I)
                     IBUF(ICOUNT,1,IDIST) = IBIN(I,1)
                  END DO
               ELSE
                  DO I = IOFF + ISTR, IOFF + ILAST
                     ICOUNT = ICOUNT + 1
                     BUF (ICOUNT,  IDIST) =  BIN(I)
                     IBUF(ICOUNT,1,IDIST) = IBIN(I,1)
                     IBUF(ICOUNT,2,IDIST) = IBIN(I,2)
                  END DO
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
               ICOUNT = 0
               IF (NIBUF .EQ. 1) THEN
                  DO I = IOFF + ISTR, IOFF + NTOT
                     ICOUNT = ICOUNT + 1
                     BUF (ICOUNT,  IDIST) =  BIN(I)
                     IBUF(ICOUNT,1,IDIST) = IBIN(I,1)
                  END DO
               ELSE
                  DO I = IOFF + ISTR, IOFF + NTOT
                     ICOUNT = ICOUNT + 1
                     BUF (ICOUNT,  IDIST) =  BIN(I)
                     IBUF(ICOUNT,1,IDIST) = IBIN(I,1)
                     IBUF(ICOUNT,2,IDIST) = IBIN(I,2)
                  END DO
               END IF
            END IF
            NCOUNT(IDIST,ISCOOR) = ICOUNT
            IOFF = IOFF + NTOT
         END DO
      END IF
C
C     Empty last buffer
C     =================
C
      IF (LAST) THEN
C
         DO IDIST = 1, NDISTR
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
         END DO
         IF (.NOT.DODIST) THEN
            WRITE (LUINTX) (  BUF(K,1),  K=1,LBUF),
     &                     ((IBUF(K,L,1),K=1,LBUF),L=1,NIBUF),-1
            REWIND LUINTX
            IF (.NOT.BDER) CALL GPCLOSE (LUINTX,'KEEP')
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
         DO INT = 1, INTS
            I = IBTAND(IBTSHR(IBIN(INT,1),3*NBITS),IBIT1)
            J = IBTAND(IBTSHR(IBIN(INT,1),2*NBITS),IBIT1)
            K = IBTAND(IBTSHR(IBIN(INT,1),  NBITS),IBIT1)
            L = IBTAND(       IBIN(INT,1),         IBIT1)
            IJ = MAX(I,J)*IBIT1 + I + J
            KL = MAX(K,L)*IBIT1 + K + L
            IBIN(INT,1) = MAX(IJ,KL)*IBIT2 + IJ + KL
         END DO
      ELSE
         DO INT = 1, INTS
            I = IBTAND(IBTSHR(IBIN(INT,1),NBITS),IBIT1)
            J = IBTAND(       IBIN(INT,1),       IBIT1)
            K = IBTAND(IBTSHR(IBIN(INT,2),NBITS),IBIT1)
            L = IBTAND(       IBIN(INT,2),       IBIT1)
            IJ = MAX(I,J)*IBIT1 + I + J
            KL = MAX(K,L)*IBIT1 + K + L
            IBIN(INT,1) = MAX(IJ,KL)
            IBIN(INT,2) = MIN(IJ,KL)
         END DO
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
      DO INT = 1, ICOUNT
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
      END DO
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
      INTEGER A, B, C, D, I, A1, B1, C1, D1, AB, CD, R, S, T
      LOGICAL DOREP(0:7,4), CTRIAB, CTRICD, CTRIAC, CTRIBD, CTRIPQ
      LOGICAL DAB, DCD, DPQ, NRDINT
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
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERINDX',-1)
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
      IF (IPRINT .GT. 15) THEN
         WRITE (LUPRI,'(2X,A,2I5)') ' NPQBCX,NPPBCX ',NPQBCX,NPPBCX
         WRITE (LUPRI,'(2X,A,2I5)') ' NPQBCS,NPPBCS ',NPQBCS,NPPBCS
         WRITE (LUPRI,'(2X,A,4I5)') 
     &                    ' NCTFX         ',NCTFA,NCTFB,NCTFC,NCTFD 
         WRITE (LUPRI,'(2X,A,2I5)') ' NCCS1 ,NCCS2  ',NCCS1, NCCS2
         WRITE (LUPRI,'(2X,A,2I5)') ' NIBUF, NBITS  ',NIBUF, NBITS 
      END IF
C
      INT = 0
C
      IF (NCCS1 .GT. 0) THEN
C
         DO A = 0, MAXREP
         IF (DOREP(A,1)) THEN
         DO B = 0, MAXREP
         IF (DOREP(B,2)) THEN
         DO C = 0, MAXREP
         D = IBTXOR(IBTXOR(A,B),C)
         IF (DOREP(C,3) .AND. DOREP(D,4)) THEN
            CD = IBTXOR(C,D)
C
            IF (DIAGAB .AND. B.GT.A) GO TO 100
            IF (DIAGCD .AND. D.GT.C) GO TO 100
C
            IF (NIBUF .EQ. 1) THEN
               IADR = 0
               DO N = 1, NPQBCS - NPPBCS
                  NA = KNDXBT(IPNTCR(N,1)) - 1
                  NB = KNDXBT(IPNTCR(N,2)) - 1
                  NC = KNDXBT(IPNTCR(N,3)) - 1
                  ND = KNDXBT(IPNTCR(N,4)) - 1
                  DO I = 1, NCTFA 
                  DO J = 1, NCTFB 
                  DO K = 1, NCTFC 
                  DO L = 1, NCTFD 
                     IADR = IADR + 1
                     INDORB(IADR,1) = (INDXBT(NA + I,A)-1)*IBITA
     &                              + (INDXBT(NB + J,B)-1)*IBITB
     &                              + (INDXBT(NC + K,C)-1)*IBITC
     &                              + (INDXBT(ND + L,D)-1)
                  END DO 
                  END DO 
                  END DO 
                  END DO
               END DO
            ELSE
               IADR = 0
               DO N = 1, NPQBCS - NPPBCS
                  NA = KNDXBT(IPNTCR(N,1)) - 1
                  NB = KNDXBT(IPNTCR(N,2)) - 1
                  NC = KNDXBT(IPNTCR(N,3)) - 1
                  ND = KNDXBT(IPNTCR(N,4)) - 1
                  DO I = 1, NCTFA 
                  DO J = 1, NCTFB 
                  DO K = 1, NCTFC 
                  DO L = 1, NCTFD 
                     IADR = IADR + 1
                     INDORB(IADR,1) = (INDXBT(NA + I,A)-1)*IBITA
     &                              + (INDXBT(NB + J,B)-1)
                     INDORB(IADR,2) = (INDXBT(NC + K,C)-1)*IBITC
     &                              + (INDXBT(ND + L,D)-1)
                  END DO
                  END DO
                  END DO
                  END DO
               END DO
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
            DO ICMPA = 1, KHKTA
            IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
            IF (IVARA.EQ.0) THEN
               IA  = IA + IBITA
               IAB = IA
               IF (CTRIAB) MAXB = ICMPA
               DO ICMPB = 1, MAXB
               IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
               IF (IVARB.EQ.0) THEN
                  IAB  = IAB + IBITB
                  IF (NIBUF.EQ.1) THEN
                     IABC = IAB
                  ELSE
                     IABC = 0
                  END IF
                  ICMPAB = IADCMP(ICMPA,ICMPB,1)
                  IODDAB = IODDCC(IPNTUV(ICMPAB,0,1))
                  DO ICMPC = 1, KHKTC
                  IVARC = IBTEST(ISTBLC,C,NHKTC,ICMPC)
                  IF (IVARC.EQ.0) THEN
                     IABC  = IABC + IBITC
                     IABCD = IABC
                     IF (CTRICD) MAXD = ICMPC
                     DO ICMPD = 1, MAXD
                     IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                     IF (IVARD.EQ.0) THEN
                        IABCD  = IABCD + 1
                        ICMPCD = IADCMP(ICMPC,ICMPD,2)
                        IODDCD = IODDCC(IPNTUV(ICMPCD,0,2))
                        IF (IODDAB .EQ. IODDCD) THEN
C
                           DAB = GCONAB .AND. CTRIAB.AND.ICMPA.EQ.ICMPB 
     &                                  .AND .NGTOAB.EQ.1 
                           DCD = GCONCD .AND. CTRICD.AND.ICMPC.EQ.ICMPD 
     &                                  .AND .NGTOCD.EQ.1
                           DPQ = .FALSE. 
C
                           IF (WRTSCR) THEN
                              IF (NIBUF .EQ. 1) THEN
                                 DO I = 1, NCCS1
                                    SOABCD = SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    IF (ABS(SOABCD) .GT. THRSH) THEN
                                    KABCD = INDORB(I,1) + IABCD
                                    IF(NRDINT(KABCD,DAB,DCD,DPQ)) THEN
                                       INT = INT + 1
                                       BIN (INT  ) = SOABCD
                                       IBIN(INT,1) = INDORB(I,1) + IABCD
C                                      CALL ERIPRI(BIN(INT),IBIN(INT,1),
C    &                                         LBIN,NIBUF,1,NBITS,0,0,0)
                                    END IF
                                    END IF
                                 END DO
                              ELSE
                                 DO I = 1, NCCS1
                                    SOABCD = SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    IF (ABS(SOABCD) .GT. THRSH) THEN
                                       INT = INT + 1
                                       BIN (INT  ) = SOABCD
                                       IBIN(INT,1) = INDORB(I,1) + IAB
                                       IBIN(INT,2) = INDORB(I,2) + IABCD
                                    END IF
                                 END DO
                              END IF
                           ELSE
                              IF (NIBUF .EQ. 1) THEN
                                 DO I = 1, NCCS1
                                 KABCD = INDORB(I,1) + IABCD
                                 IF(NRDINT(KABCD,DAB,DCD,DPQ)) THEN
                                    INT = INT + 1
                                    BIN(INT)=SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    IBIN(INT,1) = INDORB(I,1) + IABCD
C                                   CALL ERIPRI(BIN(INT),IBIN(INT,1),
C    &                                      LBIN,NIBUF,1,NBITS,0,0,0)
                                 END IF
                                 END DO
                              ELSE
                                 DO I = 1, NCCS1
                                    INT = INT + 1
                                    BIN(INT)=SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    IBIN(INT,1) = INDORB(I,1) + IAB
                                    IBIN(INT,2) = INDORB(I,2) + IABCD
                                 END DO
                              END IF
                           END IF
                        END IF
                     END IF
                     END DO
                  END IF
                  END DO
               END IF
               END DO
            END IF
            END DO
C
  100       CONTINUE
C
         END IF
         END DO
         END IF
         END DO
         END IF
         END DO
      END IF
C
      IF (NCCS2 .GT. 0) THEN
C
         DO A = 0, MAXREP
         IF (DOREP(A,1)) THEN
         DO B = 0, MAXREP
         IF (DOREP(B,2)) THEN
         DO C = 0, MAXREP
         D = IBTXOR(IBTXOR(A,B),C)
         IF (DOREP(C,3) .AND. DOREP(D,4)) THEN
            CD = IBTXOR(C,D)
C
            IF (DIAGAB .AND. B.GT.A) GO TO 200
            IF (DIAGCD .AND. D.GT.C) GO TO 200
            IF(C.GT.A .OR. (C.EQ.A .AND. D.GT.B)) GO TO 200
C
            IF (NIBUF .EQ. 1) THEN
               IADR = NCCS1
               DO N = NPQBCS - NPPBCS + 1, NPQBCS
                  NA = KNDXBT(IPNTCR(N,1)) - 1
                  NB = KNDXBT(IPNTCR(N,2)) - 1
                  NC = KNDXBT(IPNTCR(N,3)) - 1
                  ND = KNDXBT(IPNTCR(N,4)) - 1
                  DO I = 1, NCTFA
                  DO J = 1, NCTFB 
                  DO K = 1, NCTFC
                  DO L = 1, NCTFD 
                     IADR = IADR + 1
                     INDORB(IADR,1) = (INDXBT(NA + I,A)-1)*IBITA
     &                              + (INDXBT(NB + J,B)-1)*IBITB
     &                              + (INDXBT(NC + K,C)-1)*IBITC
     &                              + (INDXBT(ND + L,D)-1)
                  END DO
                  END DO
                  END DO
                  END DO
               END DO
            ELSE
               IADR = NCCS1
               DO N = NPQBCS - NPPBCS + 1, NPQBCS
                  NA = KNDXBT(IPNTCR(N,1)) - 1
                  NB = KNDXBT(IPNTCR(N,2)) - 1
                  NC = KNDXBT(IPNTCR(N,3)) - 1
                  ND = KNDXBT(IPNTCR(N,4)) - 1
                  DO I = 1, NCTFA 
                  DO J = 1, NCTFB 
                  DO K = 1, NCTFC 
                  DO L = 1, NCTFD 
                     IADR = IADR + 1
                     INDORB(IADR,1) = (INDXBT(NA + I,A)-1)*IBITA
     &                              + (INDXBT(NB + J,B)-1)
                     INDORB(IADR,2) = (INDXBT(NC + K,C)-1)*IBITC
     &                              + (INDXBT(ND + L,D)-1)
                  END DO 
                  END DO 
                  END DO 
                  END DO
               END DO
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
            DO ICMPA = 1, KHKTA
            IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
            IF (IVARA.EQ.0) THEN
               IA  = IA + IBITA
               IAB = IA
               IF (CTRIAB) MAXB = ICMPA
               DO ICMPB = 1, MAXB
               IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
               IF (IVARB.EQ.0) THEN
                  IAB    = IAB + IBITB
                  IF (NIBUF.EQ.1) THEN
                     IABC = IAB
                  ELSE
                     IABC = 0
                  END IF
                  ICMPAB = IADCMP(ICMPA,ICMPB,1)
                  IODDAB = IODDCC(IPNTUV(ICMPAB,0,1))
                  IF (CTRIAC) MAXC = ICMPA
                  DO ICMPC = 1, MAXC
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
                     DO ICMPD = 1, MAXD
                     IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                     IF (IVARD.EQ.0) THEN
                        IABCD  = IABCD + 1
                        ICMPCD = IADCMP(ICMPC,ICMPD,2)
                        IODDCD = IODDCC(IPNTUV(ICMPCD,0,2))
                        IF (IODDAB .EQ. IODDCD) THEN
C
                           DAB = CTRIAB .AND. ICMPA.EQ.ICMPB
                           DCD = CTRICD .AND. ICMPC.EQ.ICMPD
                           DPQ = CTRIPQ .AND. ICMPA.EQ.ICMPC
     &                                  .AND. ICMPB.EQ.ICMPD
C
                           IF (WRTSCR) THEN
                              IF (NIBUF.EQ.1) THEN
                                 DO I = NCCS1 + 1, NCCS
                                    SOABCD = SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    IF (ABS(SOABCD) .GT. THRSH) THEN
                                       KABCD = INDORB(I,1) + IABCD
                                    IF (NRDINT(KABCD,DAB,DCD,DPQ)) THEN
                                       INT = INT + 1
                                       BIN (INT  ) = SOABCD
                                       IBIN(INT,1) = INDORB(I,1) + IABCD
                                    END IF
                                    END IF
                                 END DO
                              ELSE
                                 DO I = NCCS1 + 1, NCCS
                                    SOABCD = SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    IF (ABS(SOABCD) .GT. THRSH) THEN
                                       INT = INT + 1
                                       BIN (INT  ) = SOABCD
                                       IBIN(INT,1) = INDORB(I,1) + IAB
                                       IBIN(INT,2) = INDORB(I,2) + IABCD
                                    END IF
                                 END DO
                              END IF
                           ELSE
                              IF (NIBUF.EQ.1) THEN
                                 DO I = NCCS1 + 1, NCCS
                                   KABCD = INDORB(I,1) + IABCD
                                   IF(NRDINT(KABCD,DAB,DCD,DPQ)) THEN
                                    INT = INT + 1
                                    BIN(INT)=SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    IBIN(INT,1) = INDORB(I,1) + IABCD
C                                   CALL ERIPRI(BIN(INT),IBIN(INT,1),
C    &                                         LBIN,NIBUF,1,NBITS,0,0,0)
                                   END IF 
                                 END DO
                              ELSE
                                 DO I = NCCS1 + 1, NCCS
                                    INT = INT + 1
                                    BIN(INT)=SO(I,R,S,T,ICMPAB,ICMPCD,1)
                                    IBIN(INT,1) = INDORB(I,1) + IAB
                                    IBIN(INT,2) = INDORB(I,2) + IABCD
C                                   CALL ERIPRI(BIN(INT),IBIN(INT,1),
C    &                                       LBIN,NIBUF,1,NBITS,0,0,0)
                                 END DO
                              END IF
                           END IF
                        END IF
                     END IF
                     END DO
                  END IF
                  END DO
               END IF
               END DO
            END IF
            END DO
C
  200       CONTINUE
C
         END IF
         END DO
         END IF
         END DO
         END IF
         END DO
      END IF
C
C     Canonical ordering of indices if requested
C     ==========================================
C
      IF (CANIND) CALL ERICAN(IBIN,INT,LBIN,NIBUF,NBITS)
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
      DO A = 0, MAXREP
         DOREP(A) = .FALSE.
      END DO
C
      DO ICMP = 1, KHKTA
         DO A = 0, MAXREP
         IF (IBTAND(ISTBLA,IBTXOR(A,ISYMAO(NHKTA,ICMP))).EQ.0) THEN
            DOREP(A) = .TRUE.
         END IF
         END DO
      END DO
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
      DO ICMPA = 1, KHKTA
      DO ICMPB = 1, KHKTB
         IF (TKMPAB) THEN
            MAXAB = MAX(ICMPA,ICMPB)
            IADCMP(ICMPA,ICMPB) = MAXAB*(MAXAB-1)/2 + MIN(ICMPA,ICMPB)
         ELSE
            IADCMP(ICMPA,ICMPB) = (ICMPA - 1)*KHKTB + ICMPB
         END IF
      END DO
      END DO
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
      DO X = 0, MAXREP
         R1 = 0
         DO R = 0, MAXREP
         IF (IBTAND(R,ISTBLR) .EQ. 0) THEN
            R1 = R1 + 1
            IF (IBTOR(X,ISTBLR) .EQ. IBTOR(R,ISTBLR)) THEN
               IPNRST(X) = R1
               GO TO 100
            END IF
         END IF
         END DO
  100    CONTINUE
      END DO
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
      DO M = 1, 4
         NSHFT = IBITX(M)
         DO IA = 1, NREP(M)
            A = IPNREP(IA,M)
            IADR = 0
            DO N = 1, NPQBCS
               NA = KNDXBT(IPNTCR(N,M)) - 1
               DO I = 1, NCTFA
               DO J = 1, NCTFB
               DO K = 1, NCTFC
               DO L = 1, NCTFD
                  IADR = IADR + 1
                  IF (M.EQ.1) NAI = NA + I
                  IF (M.EQ.2) NAI = NA + J
                  IF (M.EQ.3) NAI = NA + K
                  IF (M.EQ.4) NAI = NA + L
                  INDXTR(IADR,IA,M) = (INDXBT(NAI,A) - 1)*NSHFT
               END DO
               END DO
               END DO
               END DO
            END DO
         END DO
      END DO
C
      DO ID = 1, NREP(4)
         D = IPNREP(ID,4)
         S = IPNRST( D,2)
         DO IC = 1, NREP(3)
            C  = IPNREP(IC,3)
            CD = IBTXOR(C,D)
            T  = IPNRST(CD,3)
            DO I = 1, NCCS
               INDCD(I) = INDXTR(I,IC,3) + INDXTR(I,ID,4)
            END DO
C
            DO IB = 1, NREP(2)
               B = IPNREP(IB,2)
               R = IPNRST( B,1)
               A = IBTXOR(B,CD)
               DO IA = 1, NREP(1)
               IF (IPNREP(IA,1).EQ.A) THEN
                  IF (NIBUF.EQ.1) THEN
                     DO I = 1, NCCS
                        INDALL(I) = INDXTR(I,IA,1)
     &                            + INDXTR(I,IB,2)+INDCD(I)
                     END DO
                     DO ICD = 1, NCMPCD(C,D)
                        ICMPCD = JCMPCD(ICD,C,D)
                        IADDCD = JADDCD(ICD,C,D)
                        DO IAB = 1, NCMPAB(A,B)
                           ICMPAB = JCMPAB(IAB,A,B)
                           IABCD  = JADDAB(IAB,A,B) + IADDCD
                           DO I = 1, NCCS
                              INDORB(I,R,S,T,ICMPAB,ICMPCD,1) =
     &                                        INDALL(I) + IABCD
                           END DO
                        END DO
                     END DO
                  ELSE
                     DO I = 1, NCCS
                        INDALL(I) = INDXTR(I,IA,1)
     &                            + INDXTR(I,IB,2)
                     END DO
                     DO ICD = 1, NCMPCD(C,D)
                        ICMPCD = JCMPCD(ICD,C,D)
                        IADDCD = JADDCD(ICD,C,D)
                        DO IAB = 1, NCMPAB(A,B)
                           ICMPAB = JCMPAB(IAB,A,B)
                           IADDAB = JADDAB(IAB,A,B)
                           DO I = 1, NCCS
                              INDORB(I,R,S,T,ICMPAB,ICMPCD,1) =
     &                                       INDALL(I) + IADDAB
                              INDORB(I,R,S,T,ICMPAB,ICMPCD,2) =
     &                                       INDCD(I)  + IADDCD
                           END DO
                        END DO
                     END DO
                  END IF
               END IF
               END DO
            END DO
         END DO
      END DO
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Output from ERILBL',-1)
         DO ICMPCD = 1, KHKTCD
         DO ICMPAB = 1, KHKTAB
         IF (IODDCC(IPNTUV(ICMPAB,0,1)).EQ.
     &       IODDCC(IPNTUV(ICMPCD,0,2)))THEN
            DO T = 1, MLTPT
            DO S = 1, MLTPS
            DO R = 1, MLTPR
               DO I = 1, NCCS
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
               END DO
            END DO
            END DO
            END DO
         END IF
         END DO
         END DO
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
      DO I = 1, NREP(1)
         A = IPNREP(I,1)
         IKMPA = 0
         DO ICMPA = 1, KHKTA
         IF (IBTAND(ISTBLA,IBTXOR(A,ISYMAO(NHKTA,ICMPA))).EQ.0) THEN
            IKMPA = IKMPA + 1
            JCMP(IKMPA,I,1) = ICMPA
         END IF
         END DO
         NCMP(I,1) = IKMPA
      END DO
C
      DO I = 1, NREP(2)
         B = IPNREP(I,2)
         IKMPB = 0
         DO ICMPB = 1, KHKTB
         IF (IBTAND(ISTBLB,IBTXOR(B,ISYMAO(NHKTB,ICMPB))).EQ.0) THEN
            IKMPB = IKMPB + 1
            JCMP(IKMPB,I,2) = ICMPB
         END IF
         END DO
         NCMP(I,2) = IKMPB
      END DO
C
      DO I = 1, NREP(1)
         A = IPNREP(I,1)
         DO J = 1, NREP(2)
            B = IPNREP(J,2)
            IA = 0
            IKMPAB =  0
            DO IKMPA = 1, NCMP(I,1)
               ICMPA = JCMP(IKMPA,I,1)
               IA  = IA + IBITA
               IAB = IA
               DO IKMPB = 1, NCMP(J,2)
                  IAB    = IAB + IBITB
                  IKMPAB = IKMPAB + 1
                  JCMPAB(IKMPAB,A,B) = IADCMP(ICMPA,JCMP(IKMPB,J,2))
                  JADDAB(IKMPAB,A,B) = IAB
               END DO
            END DO
            NCMPAB(A,B) = IKMPAB
         END DO
      END DO
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
      DO A = 0, MAXREP
         DOREP(A) = .FALSE.
      END DO
C
      DO ICMP = 1, KHKTA
         DO A = 0, MAXREP
         IF (IBTAND(ISTBLA,IBTXOR(A,ISYMAO(NHKTA,ICMP))).EQ.0) THEN
            DOREP(A) = .TRUE.
         END IF
         END DO
      END DO
C
      IA = 0
      DO A = 0, MAXREP
         IF (DOREP(A)) THEN
            IA = IA + 1
            IPNREP(IA) = A
         END IF
      END DO
      NREP(1) = IA
C
      RETURN
      END
C  /* Deck erigth */
      SUBROUTINE ERIGTH(BIN,IBIN,LBIN,NIBUF,INDORB,SO,IPNTCR,IODDCC,
     &                  IPNTUV,IPNBCH,NBITS,NINTS,IOFCMP,INDXBT,ICART,
     &                  IPRINT)
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
     &          IPNBCH(NCTFA*MLTPX*KHKTB*KHKTCD,*), NINTS(NDISTR)
      DIMENSION INDXBT(MXSHEL*MXCONT,0:7)
#include <cbieri.h>
#include <ericom.h>
#include <aobtch.h>
#include <hertop.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERIGTH',-1)
C
      NBCHS = NCTFA*MLTPX*KHKTB*KHKTCD
C
      IANEW = 0
      IDIST = 0
      NIJKL = 0
      DO N = 1, NPQBCS
         NA = KNDXBT(IPNTCR(N,1)) - 1
         DO I = 1, NCTFA
         DO J = 1, NCTFB
         DO K = 1, NCTFC
         DO L = 1, NCTFD
            NIJKL = NIJKL + 1
            IAOLD = IANEW
            IANEW = INDXBT(NA + I,IRPDST)
            IF (IAOLD. NE. IANEW) THEN
               IDIST = IDIST + 1
               NCCDST(IDIST) = 1
               NCCFST(IDIST) = NIJKL 
            ELSE
               NCCDST(IDIST) = NCCDST(IDIST) + 1
            END IF
         END DO
         END DO
         END DO
         END DO
      END DO
      NDST = IDIST
C
      CALL ERICHN(INDORB,NBTCHX,LBTCHX,IPNBCH,NBCHS,NDST,NCCFST,NCCDST,
     &            IODDCC,IPNTUV,NBITS,NIBUF,IOFCMP,ICART,IPRINT)
C
      ISTR = 1
      DO I = 1, NDISTR
         NBTCH = NBTCHX(I)
         LBTCH = LBTCHX(I)
         CALL ERITRF(BIN(ISTR),IBIN(ISTR,1),NTOT,SO,INDORB,NBTCH,LBTCH,
     &               IPNBCH(1,I),NDISTR,LBIN,NIBUF)
         NINTS(I) = NTOT 
         ISTR = ISTR + NTOT 
      END DO
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
      DIMENSION NBTCHX(*), LBTCHX(*), 
     &          IPNBCH(NCTFA*MLTPX*KHKTB*KHKTCD,*),
     &          NCCFST(*), NCCDST(*), IODDCC(*), 
     &          IPNTUV(KC2MAX,0:NRDER,2),
     &          IOFCMP(MXAQN**4), INDORB(NCCS,*)
#include <ericom.h>
#include <ibtfun.h>
C
      NCMPX = 0
      MLTAB = MLTPX
      MLTCD = MLTPX*KHKTAB
      DO ICMPCD = 1, KHKTCD
         IODDCD = IODDCC(IPNTUV(ICMPCD,ICART,2))
         DO ICMPAB = 1, KHKTAB
            IODDAB = IODDCC(IPNTUV(ICMPAB,0,1))
            IF (IODDAB .EQ. IODDCD) THEN
               NCMPX = NCMPX + 1
               IOFCMP(NCMPX) = (ICMPCD - 1)*MLTCD + (ICMPAB - 1)*MLTAB
            END IF
         END DO
      END DO
C
C     Sets up NBTCHX, LBTCHX, and IPNBCH, 
C     to be used in the following way:
C
C     DO I = 1, NDISTR
C     DO J = 1, NBTCHX(I)
C        IOFF = IPNBCH(J,I)
C        DO K = 1, LBTCHX(I)
C           IADR = IOFF + K
C        END DO
C     END DO
C     END DO
C
      IBT = 2**NBITS - 1
      IF (NIBUF.EQ.1) THEN
         NB3 = 3*NBITS
      ELSE
         NB3 =   NBITS
      END IF
      CALL IZERO(NBTCHX,NDISTR)
      DO I = 1, NDST
         ISTR   = NCCFST(I)
         LENGTH = NCCDST(I)
         DO ICMPX = 1, NCMPX
            IOFF = IOFCMP(ICMPX)
            IADR = IOFF*NCCS + ISTR - 1
            DO K = IOFF + 1, IOFF + MLTPX
               IDIST = INDXDS(IBTAND(IBTSHR(INDORB(ISTR,K),NB3),IBT))
               IBTCH = NBTCHX(IDIST) + 1
               NBTCHX(IDIST) = IBTCH
               LBTCHX(IDIST) = LENGTH
               IF (IDIST.GT.0) IPNBCH(IBTCH,IDIST) = IADR
               IADR = IADR + NCCS
            END DO
         END DO
      END DO
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
      SUBROUTINE ERIBRT(SO,ISCOOR,SOPRT,INDORB,IPRINT)
C
C     T. Helgaker
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D0)
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
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERIBRT',-1)
C
      IF (NCTFAB.EQ.1 .OR. NCTFCD.EQ.1) THEN
         IF (ISCOOR .LE. 3) THEN
            DO I = 1, NAOINT
               SOPRT(I) = SO(I,ISCOOR) + SO(I,ISCOOR+3) 
            END DO
         ELSE
            DO I = 1, NAOINT
               SOPRT(I) = SO(I,ISCOOR-3) - SO(I,ISCOOR) 
            END DO
         END IF
      ELSE
         IF (ISCOOR .LE. 3) THEN
            ICOR1 = ISCOOR
            ICOR2 = ISCOOR + 3
            INT = 0
            IOFF = 0
            DO IAB = 1, KHKTAB 
            DO ICD = 1, KHKTCD 
               DO N = 1, MLTPX*NPQBCS 
                  DO I = 1, NCTFA
                  DO J = 1, NCTFB
                  DO K = 1, NCTFC
                  DO L = 1, NCTFD
                     INT = INT + 1
                     JNT = (I - 1)*NCTFB + J
     &                   + (K - 1)*NCTFAB*NCTFD
     &                   + (L - 1)*NCTFAB + IOFF
                     SOPRT(INT) = SO(JNT,ICOR1) + SO(INT,ICOR2) 
                  END DO
                  END DO
                  END DO
                  END DO
                  IOFF = IOFF + NCTFPQ
               END DO
            END DO
            END DO
         ELSE
            ICOR1 = ISCOOR - 3
            ICOR2 = ISCOOR
            INT = 0
            IOFF = 0
            DO IAB = 1, KHKTAB 
            DO ICD = 1, KHKTCD 
               DO N = 1, MLTPX*NPQBCS 
                  DO I = 1, NCTFA
                  DO J = 1, NCTFB
                  DO K = 1, NCTFC
                  DO L = 1, NCTFD
                     INT = INT + 1
                     JNT = (I - 1)*NCTFB + J
     &                   + (K - 1)*NCTFAB*NCTFD
     &                   + (L - 1)*NCTFAB + IOFF
                     SOPRT(INT) = SO(JNT,ICOR1) - SO(INT,ICOR2) 
                  END DO
                  END DO
                  END DO
                  END DO
                  IOFF = IOFF + NCTFPQ
               END DO
            END DO
            END DO
         END IF
      END IF
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Output from ERIBRT',-1)
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
      KPOINT = KLABEL + (5*NAOINT - 1)/IRAT + 1
      KLAST  = KPOINT + (4*LBIN - 1)/IRAT + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('ERIAVE',' ',KLAST,LWORK)
      LWRK = LWORK - KLAST + 1
      CALL ERIAV1(SO,DMAT,D2MAT,ID2MAT,IPNTCR,IODDCC,IPNTUV,
     &            INDXBT,WORK(KPOINT),IPRINT)
      RETURN
      END
C  /* Deck erilab */
      SUBROUTINE ERILAB(WEIGHT,LABEL,INDORB,IPNTCR,IODDCC,IPNTUV,INDXBT,
     &                  IPOINT,IREPE,ITYPE,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
      PARAMETER (D2=2.0D0, D1=1.0D0)
      INTEGER A, B, C, D, AB, CD, R, S, T
      LOGICAL DOREP(0:7,4), CTRIAB, CTRICD, DCD
      DIMENSION IPOINT(NCCS,4), 
     &          IPNTCR(MAXBCH,4), INDXBT(MXSHEL*MXCONT,0:7),
     &          IODDCC(NRTOP), IPNTUV(KC2MAX,0:NRDER,2),
     &          IPNRST(0:7,3), IODXYZ(0:3),
     &          IADCMP(MXAQN,MXAQN,2),
     &          WEIGHT(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD),
     &          LABEL(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD,4),
     &          INDORB(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD,NIBUF)
#include <cbieri.h>
#include <ericom.h>
#include <erithr.h>
#include <eridst.h>
#include <eribuf.h>
#include <aobtch.h>
#include <hertop.h>
#include <symmet.h>
C
#include <ibtfun.h>
C
      IBTEST(I,J,K,L) = IBTAND(I,IBTXOR(J,ISYMAO(K,L)))
C
      IF (IPRINT .GT. 6) THEN
         CALL HEADER('Subroutine ERILA1',-1)
         IF (IPRINT .GT. 15) THEN
            WRITE (LUPRI,'(2X,A,2I5)') ' NPQBCS,NPPBCS ',NPQBCS,NPPBCS
         END IF
      END IF
C
      CALL PRPREP(DOREP(0,1),NHKTA,KHKTA,ISTBLA)
      CALL PRPREP(DOREP(0,2),NHKTB,KHKTB,ISTBLB)
      CALL PRPREP(DOREP(0,3),NHKTC,KHKTC,ISTBLC)
      CALL PRPREP(DOREP(0,4),NHKTD,KHKTD,ISTBLD)
C
      IRPDST = 0
      DO A = 0, MAXREP
         IF (DOREP(A,1)) THEN
            IRPDST = A
            GO TO 10
         END IF
      END DO
   10 CONTINUE
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
      KCCS1 = 1
      KCCS2 = NCCS1 + 1
C
      IF (EXPERI) CALL DZERO(WEIGHT,NAOINT)
C
      IBIT1 = 2**NBITS     - 1
      IBIT2 = 2**(2*NBITS) - 1
C
      DO A = 0, MAXREP
      IF (DOREP(A,1)) THEN
      DO B = 0, MAXREP
      IF (DOREP(B,2)) THEN
      DO C = 0, MAXREP
      D = IBTXOR(IBTXOR(IBTXOR(A,B),C),IREPE)
      IF (DOREP(C,3) .AND. DOREP(D,4)) THEN
         CD = IBTXOR(C,D)
C
         CTRIAB = DIAGAB .AND. A.EQ.B .AND. EXPERI 
         CTRICD = DIAGCD .AND. C.EQ.D .AND. EXPERI 
C
         R = IPNRST(B,1)
         S = IPNRST(D,2)
         T = IPNRST(CD,3)
         CALL ERIPNT(IPOINT,A,B,C,D,IPNTCR,INDXBT,ITYPE)
         IF (NCCS1 .GT. 0) THEN
C
C           IF (DIAGAB .AND. B.GT.A) GO TO 100
C           IF (DIAGCD .AND. D.GT.C) GO TO 100
C
            IA   = -1
            MAXB = KHKTB
            MAXD = KHKTD
            DO ICMPA = 1, KHKTA
            IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
            IF (IVARA.EQ.0) THEN
               IA = IA + 1
               IB = -1
               IF (CTRIAB) MAXB = ICMPA
               DO ICMPB = 1, MAXB
               IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
               IF (IVARB.EQ.0) THEN
                  IB = IB + 1
                  IC = -1
                  ICMPAB = IADCMP(ICMPA,ICMPB,1)
                  IODXYZ(0) = IODDCC(IPNTUV(ICMPAB,0,1))
                  IODXYZ(1) = IODDCC(IPNTUV(ICMPAB,1,1))
                  IODXYZ(2) = IODDCC(IPNTUV(ICMPAB,2,1))
                  IODXYZ(3) = IODDCC(IPNTUV(ICMPAB,3,1))
                  DO ICMPC = 1, KHKTC
                  IVARC = IBTEST(ISTBLC,C,NHKTC,ICMPC)
                  IF (IVARC.EQ.0) THEN
                     IC = IC + 1
                     ID = -1
                     IF (CTRICD) MAXD = ICMPC
                     DO ICMPD = 1, MAXD
                     IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                     IF (IVARD.EQ.0) THEN
                        ID = ID + 1
                        ICMPCD = IADCMP(ICMPC,ICMPD,2)
                        IODDCD = IODDCC(IPNTUV(ICMPCD,0,2))
                        DCD = GCONCD .AND. CTRICD.AND.ICMPC.EQ.ICMPD
                        CALL ERIBRC(IA,IB,IC,ID,
     &                              WEIGHT(KCCS1,R,S,T,ICMPAB,ICMPCD),
     &                              LABEL (KCCS1,R,S,T,ICMPAB,ICMPCD,1),
     &                              INDORB(KCCS1,R,S,T,ICMPAB,ICMPCD,1),
     &                              IPOINT(KCCS1,1),NCCS1,D2,IODDCD,
     &                              IODXYZ,DCD,IBIT1,IBIT2,IPRINT)
                     END IF
                     END DO
                  END IF
                  END DO
               END IF
               END DO
            END IF
            END DO
  100       CONTINUE
         END IF
         IF (NCCS2 .GT. 0) THEN
C
C           IF (DIAGAB .AND. B.GT.A) GO TO 200
C           IF (DIAGCD .AND. D.GT.C) GO TO 200
C           IF (C.GT.A .OR. (C.EQ.A .AND. D.GT.B)) GO TO 200
C
            IA = -1
            MAXB = KHKTB
            MAXD = KHKTD
            DO ICMPA = 1, KHKTA
            IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
            IF (IVARA.EQ.0) THEN
               IA = IA + 1
               IB = -1
C              IF (CTRIAB) MAXB = ICMPA
               DO ICMPB = 1, MAXB
               IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
               IF (IVARB.EQ.0) THEN
                  IB = IB + 1
                  IC = -1
                  ICMPAB = IADCMP(ICMPA,ICMPB,1)
                  IODXYZ(0) = IODDCC(IPNTUV(ICMPAB,0,1))
                  IODXYZ(1) = IODDCC(IPNTUV(ICMPAB,1,1))
                  IODXYZ(2) = IODDCC(IPNTUV(ICMPAB,2,1))
                  IODXYZ(3) = IODDCC(IPNTUV(ICMPAB,3,1))
                  DO ICMPC = 1, KHKTC
                  IVARC = IBTEST(ISTBLC,C,NHKTC,ICMPC)
                  IF (IVARC.EQ.0) THEN
                     IC = IC + 1
                     ID = -1
C                    IF (CTRICD) MAXD = ICMPC
                     DO ICMPD = 1, MAXD
                     IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                     IF (IVARD.EQ.0) THEN
                        ID = ID + 1
                        ICMPCD = IADCMP(ICMPC,ICMPD,2)
                        IODDCD = IODDCC(IPNTUV(ICMPCD,0,2))
                        DCD = GCONCD .AND. CTRICD.AND.ICMPC.EQ.ICMPD
                        CALL ERIBRC(IA,IB,IC,ID,
     &                              WEIGHT(KCCS2,R,S,T,ICMPAB,ICMPCD),
     &                              LABEL (KCCS2,R,S,T,ICMPAB,ICMPCD,1),
     &                              INDORB(KCCS2,R,S,T,ICMPAB,ICMPCD,1),
     &                              IPOINT(KCCS2,1),NCCS2,D1,IODDCD,
     &                              IODXYZ,DCD,IBIT1,IBIT2,IPRINT)
                     END IF
                     END DO
                  END IF
                  END DO
               END IF
               END DO
            END IF
            END DO
  200       CONTINUE
         END IF
      END IF
      END DO
      END IF
      END DO
      END IF
      END DO
C
      RETURN
      END
C  /* Deck eribrc */
      SUBROUTINE ERIBRC(IA,IB,IC,ID,WEIGHT,LABEL,INDORB,IPOINT,NINTS,
     &                  DFAC,IODDAB,IODXYZ,DCD,IBIT1,IBIT2,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
      PARAMETER (DP5 = 0.5D0, D4 = 4.0D0)
      LOGICAL DCD, DOXYZ
      DIMENSION WEIGHT(NAOINT), LABEL(NAOINT,4), INDORB(NAOINT,2),
     &          IPOINT(NCCS,4), IODXYZ(0:3)
#include <cbieri.h>
#include <ericom.h>
#include <eribuf.h>
#include <aobtch.h>
#include <hertop.h>
#include <nuclei.h>
#include <symmet.h>
C    
      DOXYZ = IODDAB.EQ.IODXYZ(1) .OR. IODDAB.EQ.IODXYZ(2) 
     &                            .OR. IODDAB.EQ.IODXYZ(3)
C
      IF (EXPERI .AND. DOXYZ) THEN
         DO I = 1, NINTS
            KA = IPOINT(I,1) + IA
            KB = IPOINT(I,2) + IB
            KC = IPOINT(I,3) + IC
            KD = IPOINT(I,4) + ID
            IF (CCRUN) THEN
               FAC = DFAC
               IF (DCD .OR. KC.EQ.KD) FAC = DP5*FAC
            ELSE
               FAC = DFAC*D4
               IF (KA.EQ.KB) FAC = DP5*FAC
               IF (DCD .OR. KC.EQ.KD) FAC = DP5*FAC
            END IF
            WEIGHT(I)  = FAC
            LABEL(I,1) = KA
            LABEL(I,2) = KB
            LABEL(I,3) = KC
            LABEL(I,4) = KD
         END DO
      ELSE IF (DOXYZ .OR. IODDAB.EQ.IODXYZ(0)) THEN 
         DO I = 1, NINTS
            KA  = IPOINT(I,1) + IA
            KB  = IPOINT(I,2) + IB
            KC  = IPOINT(I,3) + IC
            KD  = IPOINT(I,4) + ID
            KAB = KA*IBIT1 + KA + KB
            KCD = KC*IBIT1 + KC + KD
            IF (NIBUF.EQ.1) THEN
               INDORB(I,1) = KAB*IBIT2 + KAB + KCD
            ELSE
               INDORB(I,1) = KAB
               INDORB(I,2) = KCD
            END IF
            LABEL(I,1) = ICNTAO(KA)
            LABEL(I,2) = ICNTAO(KB)
            LABEL(I,3) = ICNTAO(KC)
            LABEL(I,4) = ICNTAO(KD)
         END DO
      END IF
      RETURN
      END
C  /* Deck erilbd */
      SUBROUTINE ERIAV1(SO,DMAT,D2MAT,ID2MAT,IPNTCR,IODDCC,IPNTUV,
     &                  INDXBT,IPOINT,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
      PARAMETER (D2=2.0D0, D1=1.0D0)
      INTEGER A, B, C, D, AB, CD, R, S, T
      LOGICAL DOREP(0:7,4), CTRIAB, CTRICD, DCD
      DIMENSION IPOINT(NCCS,5), 
     &          IPNTCR(MAXBCH,4), INDXBT(MXSHEL*MXCONT,0:7),
     &          IODDCC(NRTOP), IPNTUV(KC2MAX,0:NRDER,2),
     &          IPNRST(0:7,3), IODXYZ(3),
     &          IADCMP(MXAQN,MXAQN,2),
     &          SO(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD,12),
     &          DMAT(NBASE,NBASE),
     &          D2MAT(NBASE,NBASE,NPDIMB,NPDIMA),
     &          ID2MAT(MXCORB,2)
#include <cbieri.h>
#include <ericom.h>
#include <erithr.h>
#include <eribuf.h>
#include <aobtch.h>
#include <hertop.h>
#include <symmet.h>
#include <nuclei.h>
#include <energy.h>
C
#include <ibtfun.h>
      IBTEST(I,J,K,L) = IBTAND(I,IBTXOR(J,ISYMAO(K,L)))
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERIAV1',-1)
C
      IF (IPRINT .GT. 15) THEN
         WRITE (LUPRI,'(2X,A,2I5)') ' NPQBCS,NPPBCS ',NPQBCS,NPPBCS
      END IF
C
      IF (GRDZER) CALL DZERO(GRADEE,MXCOOR)
C
      CALL PRPREP(DOREP(0,1),NHKTA,KHKTA,ISTBLA)
      CALL PRPREP(DOREP(0,2),NHKTB,KHKTB,ISTBLB)
      CALL PRPREP(DOREP(0,3),NHKTC,KHKTC,ISTBLC)
      CALL PRPREP(DOREP(0,4),NHKTD,KHKTD,ISTBLD)
C
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(/,2X,A,8L2)')'DOREP A  ',(DOREP(I,1),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP B  ',(DOREP(I,2),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP C  ',(DOREP(I,3),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP D  ',(DOREP(I,4),I=0,MAXREP)
      END IF
C
      CALL CMPADR(IADCMP(1,1,1),KHKTA,KHKTB,TKMPAB)
      CALL CMPADR(IADCMP(1,1,2),KHKTC,KHKTD,TKMPCD)
C
      CALL GETRST(IPNRST(0,1),ISTBLR)
      CALL GETRST(IPNRST(0,2),ISTBLS)
      CALL GETRST(IPNRST(0,3),ISTBLT)
C
      NCCS1 = (NPQBCS - NPPBCS)*NCTFAB*NCTFCD
      NCCS2 = NPPBCS*NCTFAB*NCTFCD
      KCCS1 = 1
      KCCS2 = NCCS1 + 1
C
C     loop over irreps
C
      DO A = 0, MAXREP
      IF (DOREP(A,1)) THEN
      DO B = 0, MAXREP
      IF (DOREP(B,2)) THEN
      DO C = 0, MAXREP
      IF (DOREP(C,3) .AND. DOREP(IBTXOR(IBTXOR(A,B),C),4)) THEN
         D = IBTXOR(IBTXOR(A,B),C)
         CD = IBTXOR(C,D)
C
         CTRIAB = DIAGAB .AND. A.EQ.B 
         CTRICD = DIAGCD .AND. C.EQ.D 
C
         R = IPNRST(B,1)
         S = IPNRST(D,2)
         T = IPNRST(CD,3)
C
         CALL ERIPNT(IPOINT,A,B,C,D,IPNTCR,INDXBT,0)
C
         IF (NCCS1 .GT. 0) THEN
C
C           IF (DIAGAB .AND. B.GT.A) GO TO 100
            IF (DIAGCD .AND. D.GT.C) GO TO 100
C
            IA   = -1
            MAXB = KHKTB
            MAXD = KHKTD
            DO ICMPA = 1, KHKTA
            IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
            IF (IVARA.EQ.0) THEN
               IA = IA + 1
               IB = -1
               IF (CTRIAB) MAXB = ICMPA
               DO ICMPB = 1, MAXB
               IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
               IF (IVARB.EQ.0) THEN
                  IB = IB + 1
                  IC = -1
                  ICMPAB = IADCMP(ICMPA,ICMPB,1)
                  IODXYZ(1) = IODDCC(IPNTUV(ICMPAB,1,1))
                  IODXYZ(2) = IODDCC(IPNTUV(ICMPAB,2,1))
                  IODXYZ(3) = IODDCC(IPNTUV(ICMPAB,3,1))
                  DO ICMPC = 1, KHKTC
                  IVARC = IBTEST(ISTBLC,C,NHKTC,ICMPC)
                  IF (IVARC.EQ.0) THEN
                     IC = IC + 1
                     ID = -1
                     IF (CTRICD) MAXD = ICMPC
                     DO ICMPD = 1, MAXD
                     IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                     IF (IVARD.EQ.0) THEN
                        ID = ID + 1
                        ICMPCD = IADCMP(ICMPC,ICMPD,2)
                        IODDCD = IODDCC(IPNTUV(ICMPCD,0,2))
                        DCD = GCONCD .AND. CTRICD.AND.ICMPC.EQ.ICMPD
                        CALL ERIFRC(IA,IB,IC,ID,
     &                              SO(KCCS1,R,S,T,ICMPAB,ICMPCD,1),
     &                              IPOINT(KCCS1,1),NCCS1,
     &                              DMAT,D2MAT,ID2MAT,D1,
     &                              IODDCD,IODXYZ,DCD,IPRINT)
                     END IF
                     END DO
                  END IF
                  END DO
               END IF
               END DO
            END IF
            END DO
C
  100       CONTINUE
C
         END IF
         IF (NCCS2 .GT. 0) THEN
C
C           IF (DIAGAB .AND. B.GT.A) GO TO 200
C           IF (DIAGCD .AND. D.GT.C) GO TO 200
C           IF (C.GT.A .OR. (C.EQ.A .AND. D.GT.B)) GO TO 200
C
            IA = -1
            MAXB = KHKTB
            MAXD = KHKTD
            DO ICMPA = 1, KHKTA
            IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
            IF (IVARA.EQ.0) THEN
               IA = IA + 1
               IB = -1
C              IF (CTRIAB) MAXB = ICMPA
               DO ICMPB = 1, MAXB
               IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
               IF (IVARB.EQ.0) THEN
                  IB = IB + 1
                  IC = -1
                  ICMPAB = IADCMP(ICMPA,ICMPB,1)
                  IODXYZ(1) = IODDCC(IPNTUV(ICMPAB,1,1))
                  IODXYZ(2) = IODDCC(IPNTUV(ICMPAB,2,1))
                  IODXYZ(3) = IODDCC(IPNTUV(ICMPAB,3,1))
                  DO ICMPC = 1, KHKTC
                  IVARC = IBTEST(ISTBLC,C,NHKTC,ICMPC)
                  IF (IVARC.EQ.0) THEN
                     IC = IC + 1
                     ID = -1
C                    IF (CTRICD) MAXD = ICMPC
                     DO ICMPD = 1, MAXD
                     IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                     IF (IVARD.EQ.0) THEN
                        ID = ID + 1
                        ICMPCD = IADCMP(ICMPC,ICMPD,2)
                        IODDCD = IODDCC(IPNTUV(ICMPCD,0,2))
                        DCD = GCONCD .AND. CTRICD.AND.ICMPC.EQ.ICMPD
                        CALL ERIFRC(IA,IB,IC,ID,
     &                              SO(KCCS2,R,S,T,ICMPAB,ICMPCD,1),
     &                              IPOINT(KCCS2,1),NCCS2,
     &                              DMAT,D2MAT,ID2MAT,D2,
     &                              IODDCD,IODXYZ,DCD,IPRINT)
                     END IF
                     END DO
                  END IF
                  END DO
               END IF
               END DO
            END IF
            END DO
C
  200       CONTINUE
         END IF
C
      END IF
      END DO !C
      END IF
      END DO !B
      END IF
      END DO !A
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('GRADEE in ERIAV1',-1) 
         CALL OUTPUT(GRADEE,1,3,1,NUCDEP,3,NUCDEP,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck eripnt */
      SUBROUTINE ERIPNT(IPOINT,A,B,C,D,IPNTCR,INDXBT,ITYPE)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
      INTEGER A,B,C,D
      DIMENSION IPOINT(NCCS,5), 
     &          IPNTCR(MAXBCH,4), INDXBT(MXSHEL*MXCONT,0:7)
#include <cbieri.h>
#include <ericom.h>
#include <aobtch.h>
C
      IADR = 0
      DO N = 1, NPQBCS
         NA = KNDXBT(IPNTCR(N,1)) - 1
         NB = KNDXBT(IPNTCR(N,2)) - 1
         NC = KNDXBT(IPNTCR(N,3)) - 1
         ND = KNDXBT(IPNTCR(N,4)) - 1
         IF (ITYPE.NE.2) THEN
            DO I = 1, NCTFA
            DO J = 1, NCTFB
            DO K = 1, NCTFC
            DO L = 1, NCTFD
               IADR = IADR + 1
               IPOINT(IADR,1) = INDXBT(NA + I,A)
               IPOINT(IADR,2) = INDXBT(NB + J,B)
               IPOINT(IADR,3) = INDXBT(NC + K,C)
               IPOINT(IADR,4) = INDXBT(ND + L,D)
            END DO
            END DO
            END DO
            END DO
         ELSE
            DO K = 1, NCTFC
            DO L = 1, NCTFD
            DO I = 1, NCTFA
            DO J = 1, NCTFB
               IADR = IADR + 1
               IPOINT(IADR,1) = INDXBT(NA + I,A)
               IPOINT(IADR,2) = INDXBT(NB + J,B)
               IPOINT(IADR,3) = INDXBT(NC + K,C)
               IPOINT(IADR,4) = INDXBT(ND + L,D)
            END DO
            END DO
            END DO
            END DO
         END IF
      END DO
C
      IF (ITYPE.EQ.0) THEN
         IADR = 0
         DO N = 1, NPQBCS
            DO I = 1, NCTFA
            DO J = 1, NCTFB
            DO K = 1, NCTFC
            DO L = 1, NCTFD
               IADR = IADR + 1
               IPOINT(IADR,5) = (I - 1)*NCTFB + J
     &                        + (K - 1)*NCTFAB*NCTFD
     &                        + (L - 1)*NCTFAB 
     &                        + (N - 1)*NCTFPQ
            END DO
            END DO
            END DO
            END DO
         END DO
      END IF
      RETURN
      END
C  /* Deck erifrc */
      SUBROUTINE ERIFRC(IA,IB,IC,ID,SO,IPOINT,NINTS,DMAT,D2MAT,ID2MAT,
     &                  DFAC,IODDAB,IODXYZ,DCD,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
      PARAMETER (D4=4.0D0, D2=2.0D0, D1=1.0D0, DP5=0.5D0, DP25=0.25D0)
      LOGICAL DCD,  NOTRAS
      DIMENSION SO(NAOINT,3,4), IPOINT(NCCS,5), DMAT(NBASE,NBASE), 
     &          D2MAT(NBASE,NBASE,NPDIMB,NPDIMA),
     &          ID2MAT(MXCORB,2),
     &          IODXYZ(3)
#include <cbieri.h>
#include <ericom.h>
#include <eribuf.h>
#include <aobtch.h>
#include <hertop.h>
#include <nuclei.h>
#include <symmet.h>
#include <energy.h>
C    
      IF (IODDAB.EQ.IODXYZ(1) .OR. IODDAB.EQ.IODXYZ(2) 
     &                        .OR. IODDAB.EQ.IODXYZ(3)) THEN
C
         NOTRAS = .TRUE.
C
         DO I = 1, NINTS
C
            KA = IPOINT(I,1) + IA
            KB = IPOINT(I,2) + IB
            KC = IPOINT(I,3) + IC
            KD = IPOINT(I,4) + ID
            IF (CCRUN) THEN
               PMAT = DFAC*D2MAT(KD,KC,ID2MAT(KB,2),ID2MAT(KA,1))
               IF (DCD .OR. KC.EQ.KD) PMAT = DP5*PMAT
            ELSE
               PMAT = DFAC*D4*(DMAT(KA,KB)*DMAT(KC,KD) 
     &                  - DP25*DMAT(KA,KC)*DMAT(KB,KD)
     &                  - DP25*DMAT(KA,KD)*DMAT(KB,KC))
               IF (KA .EQ. KB)        PMAT = DP5*PMAT
               IF (DCD .OR. KC.EQ.KD) PMAT = DP5*PMAT
            END IF
C
            JA = 3*ICNTAO(KA) - 3
            JB = 3*ICNTAO(KB) - 3
            JC = 3*ICNTAO(KC) - 3
            JD = 3*ICNTAO(KD) - 3
C
            IF (NOTRAS) THEN
               J = IPOINT(I,5)
               DO ICOOR = 1, 3
               IF (IODDAB .EQ. IODXYZ(ICOOR)) THEN
                  ISCORA = IPTCNT(JA+ICOOR,0,1)
                  ISCORB = IPTCNT(JB+ICOOR,0,1)
                  ISCORC = IPTCNT(JC+ICOOR,0,1)
                  ISCORD = IPTCNT(JD+ICOOR,0,1)
                  IF (ISCORA.GT.0) GRADEE(ISCORA) 
     &                           = GRADEE(ISCORA) + PMAT*SO(J,ICOOR,1)
                  IF (ISCORB.GT.0) GRADEE(ISCORB) 
     &                           = GRADEE(ISCORB) + PMAT*SO(J,ICOOR,2)
                  IF (ISCORC.GT.0) GRADEE(ISCORC) 
     &                           = GRADEE(ISCORC) + PMAT*SO(I,ICOOR,3)
                  IF (ISCORD.GT.0) GRADEE(ISCORD) 
     &                           = GRADEE(ISCORD) + PMAT*SO(I,ICOOR,4)
               END IF
               END DO
            ELSE
               DO ICOOR = 1, 3
               IF (IODDAB .EQ. IODXYZ(ICOOR)) THEN
                  FA = PMAT*SO(J,ICOOR,1)
                  FB = PMAT*SO(I,ICOOR,2)
                  FC = PMAT*SO(I,ICOOR,3)
                  FD = - FA - FB - FC
                  ISCORA = IPTCNT(JA+ICOOR,0,1)
                  ISCORB = IPTCNT(JB+ICOOR,0,1)
                  ISCORC = IPTCNT(JC+ICOOR,0,1)
                  ISCORD = IPTCNT(JD+ICOOR,0,1)
                  IF (ISCORA.GT.0) GRADEE(ISCORA) = GRADEE(ISCORA) + FA 
                  IF (ISCORB.GT.0) GRADEE(ISCORB) = GRADEE(ISCORB) + FB 
                  IF (ISCORC.GT.0) GRADEE(ISCORC) = GRADEE(ISCORC) + FC 
                  IF (ISCORD.GT.0) GRADEE(ISCORD) = GRADEE(ISCORD) + FD 
               END IF
               END DO
            END IF
         END DO
      END IF
      RETURN
      END
C  /* Deck nrdint */
      FUNCTION NRDINT(INDEX,DAB,DCD,DPQ)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      LOGICAL NRDINT, DAB, DCD, DPQ
#include <eribuf.h>
#include <ibtpar.h>
#include <ibtfun.h>
      IBIT1 = 2**NBITS - 1
      I = IBTAND(IBTSHR(INDEX,3*NBITS),IBIT1)
      J = IBTAND(IBTSHR(INDEX,2*NBITS),IBIT1)
      K = IBTAND(IBTSHR(INDEX,  NBITS),IBIT1)
      L = IBTAND(       INDEX,         IBIT1)
      IJ = MAX(I,J)*IBIT1 + I + J
      KL = MAX(K,L)*IBIT1 + K + L
      NRDINT = .NOT.((DAB .AND. I.GT.J) .OR.
     &               (DCD .AND. K.GT.L) .OR.
     &               (DPQ .AND. IJ.GT.KL))
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
         DO I = 1, NAOINT
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
         END DO
      ELSE
         IOFF = 0
         DO K = 1, KHKTCD
            IODDCX = IODDCC(IPNTUV(K,1,2))
            IODDCY = IODDCC(IPNTUV(K,2,2))
            IODDCZ = IODDCC(IPNTUV(K,3,2))
            DO J = 1, KHKTAB
               IODDAB = IODDCC(IPNTUV(J,0,1))
               IF (IODDAB.EQ.IODDCX .OR. IODDAB.EQ.IODDCY .OR.  
     &                                   IODDAB.EQ.IODDCZ) THEN
                  IF (FROMD1) THEN
                     DO I = IOFF + 1, IOFF + NCCT
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
                     END DO
                   ELSE
                     DO I = IOFF + 1, IOFF + NCCT
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
                     END DO
                   END IF
               END IF
               IOFF = IOFF + NCCT
            END DO
         END DO
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
            DO A = 1, 4
               DO I = 1, NAOINT
                  IAX = 3*LABEL(I,A) - 2
                  IAY = IAX + 1
                  IAZ = IAX + 2
                  GRADEE(IAX) = GRADEE(IAX) + PMAT(I)*SO(I,1,A)
                  GRADEE(IAY) = GRADEE(IAY) + PMAT(I)*SO(I,2,A)
                  GRADEE(IAZ) = GRADEE(IAZ) + PMAT(I)*SO(I,3,A)
               END DO
            END DO
         ELSE
            DO A = 1, 4
            DO X = 1, 3
               IOFF = 0
               DO K = 1, KHKTCD
                  IODDCD = IODDCC(IPNTUV(K,X,2))
                  DO J = 1, KHKTAB
                     IF (IODDCC(IPNTUV(J,0,1)).EQ.IODDCD) THEN
                        DO I = IOFF + 1, IOFF + NCCT
                           IA = 3*LABEL(I,A) - 3 + X
                           GRADEE(IA) = GRADEE(IA) + PMAT(I)*SO(I,X,A)
                        END DO
                     END IF
                     IOFF = IOFF + NCCT
                  END DO
               END DO
            END DO
            END DO
         END IF
      ELSE
         DO J = 1, 3
         DO I = 1, NAOINT
            GRADEE(J) = GRADEE(J) + PMAT(I)*SO(I,J,1)
         END DO
         END DO
      END IF
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('GRADEE in ERIGRD',-1) 
         CALL OUTPUT(GRADEE,1,3,1,NUCDEP,3,NUCDEP,1,LUPRI)
      END IF
      RETURN
      END
