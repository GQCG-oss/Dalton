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
C$Id: eri2dis.F,v 1.1.1.1 2001-02-08 13:33:13 hjj Exp $
C
#include <single.h>
C  /* Deck eridi1 */
      SUBROUTINE ERIDI1(KODCL1,KODCL2,
     &                  KODBC1,KODBC2,KRDBC1,KRDBC2,
     &                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     &                  KFREE,LFREE,KEND,CCFBT,INDXBT,
     &                  WORK,LWORK,IPRINT)
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxaqn.h>
#include <mxorb.h>
      DIMENSION CCFBT(*), INDXBT(*), WORK(LWORK)
#include <ccom.h>
#include <cbieri.h>
#include <eribuf.h>
#include <ericom.h>
#include <erithr.h>
#include <erimem.h>
#include <aobtch.h>
#include <odbtch.h>
#include <symmet.h>
#include <inftap.h>
#include <iratdef.h>
C
      IPRINT = IPRERI
      THRSH  = MAX(THRS,1.00D-15)
C
      DODIST = .TRUE.
      PMSAB  = .FALSE.
      PMS12  = .FALSE.
      CANIND = .FALSE.
C
      IF (LBFINP .EQ. 600) LBFINP = 20 000
C
C     Memory
C
      MEMOK  = .TRUE.
      MEMADD = 0
      MODAB  = 0
      MODCD  = 0
C
C     Symmetry transformation
C
      SYMTRA = DOSTRA .OR. MAXREP .GT. 0
C
      WRTINT = .TRUE.
      ERITWO = .TRUE.
      FCKINT = .FALSE.
C
C     Pi
C
      CALL GETPI
C
C     AO batches
C     ==========
C
      CALL SETAOB(CCFBT,INDXBT,WORK,LWORK,IPRINT)
C
C     Open files
C
C     Note: NBASIS is unknown until after SETAOB has been called
C
      CALL GPOPEN(LUINTR,'AOTWODIS','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      CALL GPOPEN(LUAORC,' ','UNKNOWN','DIRECT','UNFORMATTED',
     &            LBFINP,OLDDX)
C
C     OD batches
C     ==========
C
C     This subroutine returns several arrays for each electron
C     starting at addresses K????1 and K????2. These are to be
C     transferred to ODCDRV.
C
      CALL ODBCHS(KODCL1,KODCL2,
     &            KODBC1,KODBC2,KRDBC1,KRDBC2,
     &            KODPP1,KODPP2,KRDPP1,KRDPP2,
     &            KFREE,LFREE,CCFBT,WORK,LWORK,IPRINT)
C
      KODCL1 = KODCL1 + KEND - 1
      KODCL2 = KODCL2 + KEND - 1
      KODBC1 = KODBC1 + KEND - 1
      KODBC2 = KODBC2 + KEND - 1
      KRDBC1 = KRDBC1 + KEND - 1
      KRDBC2 = KRDBC2 + KEND - 1
      KODPP1 = KODPP1 + KEND - 1
      KODPP2 = KODPP2 + KEND - 1
      KRDPP1 = KRDPP1 + KEND - 1
      KRDPP2 = KRDPP2 + KEND - 1
      KFREE  = KFREE  + KEND - 1
C
      IF (IPRINT .GT. 2) THEN
         WRITE (LUPRI,'(2(/,2X,A,I10))')
     &      ' Memory requirements for ODBCHS:',LWORK - LFREE,
     &      ' Memory left for ODCDRV:        ',LFREE
      END IF
C
      ICALL = 0
      CALL GETDST(ICALL,IPRINT)
C
      RETURN
      END
C  /* Deck eridi2 */
      SUBROUTINE ERIDI2(ICALL,INDEXA,NUMDIS,
     &                  IODCL1,IODCL2,
     &                  IODBC1,IODBC2,RODBC1,RODBC2,
     &                  IODPP1,IODPP2,RODPP1,RODPP2,
     &                  CCFBT,INDXBT,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <eridst.h>
#include <dummy.h>
      DIMENSION IODCL1(NODCL1,NITCL), IODCL2(NODCL2,NITCL),
     &          IODBC1(NODBC1,NITBC), IODBC2(NODBC2,NITBC),
     &          RODBC1(NODBC1,NRTBC), RODBC2(NODBC2,NRTBC),
     &          IODPP1(NODPP1,NITPP), IODPP2(NODPP2,NITPP),
     &          RODPP1(NODPP1,NRTPP), RODPP2(NODPP2,NRTPP),
     &          INDEXA(*), CCFBT(*), INDXBT(*), WORK(LWORK)
#include <ccom.h>
#include <cbieri.h>
#include <ericom.h>
#include <erithr.h>
#include <erimem.h>
#include <aobtch.h>
#include <odbtch.h>
#include <odclss.h>
#include <symmet.h>
#include <eribuf.h>
#include <iratdef.h>
C
      IPRINT = IPRERI
C
      NEWDIS = .TRUE.
C
C     Memory
C
      MEMOK  = .TRUE.
      MEMADD = 0
      MODAB  = 0
      MODCD  = 0
C
C     Symmetry transformation
C
      SYMTRA = DOSTRA .OR. MAXREP .GT. 0
C
      CALL GETDST(ICALL,IPRINT)
C
C     Select integrals to be calculated
C     =================================
C
      CALL PICKAO(IPRINT)
C
C     Information about distributions
C     ===============================
C
      CALL ERIDSI(INDXBT,IPRINT)
C
C     Transfer information
C     ====================
C
      NUMDIS = NDISTR
      CALL ICOPY(NUMDIS,INDDST,1,INDEXA,1)
      IF ( IPRINT .GT. 2 ) THEN
         CALL HEADER('Output from ERIDI2',2)
         WRITE (LUPRI,'(2X,I3,A,8I5/,(40X,8I5))')
     &         NDISTR,' distributions in this HERMIT call:',
     &         (INDEXA(I),I=1,NDISTR)
      ENDIF
C
C     Calculate integrals
C     ===================
C
      IF (.NOT.INTSKP) THEN
         CALL TIMER('START ',TIMSTR,TIMEND)
         CALL ODCDRV(IODCL1,IODCL2,
     &               IODBC1,IODBC2,RODBC1,RODBC2,
     &               IODPP1,IODPP2,RODPP1,RODPP2,
     &               DUMMY,DUMMY,IDUMMY,IDUMMY,
     &               CCFBT,INDXBT,WORK,LWORK,IPRINT)
         IF ( IPRINT .GT. 2) CALL TIMER('HERMIT',TIMSTR,TIMEND)
         CALL FLSHFO(LUPRI)
C
C        Error message in case of insufficient memory
C
         IF (.NOT.MEMOK) THEN
            WRITE (LUPRI,'(//,1X,A,3(/,1X,A,I10))')
     &         ' Not enough memory for this run of ERIDI2.',
     &         ' Available memory in ERIDI2:',LWORK,
     &         ' Required memory for ERIDI2:',LWORK + MEMADD,
     &         ' Increase memory (LWORK) by:',MEMADD
            WRITE (LUPRI,'(/,1X,A,2I5)')
     &         ' Memory requirements largest for OD classes :',
     &           MODAB,MODCD
            CALL QUIT('Insufficient memory in ERIDI2.')
         END IF
      END IF
C
      RETURN
      END
C  /* Deck eriact */
      SUBROUTINE ERIACT(NACTOD,SCREEN,RODBCH,IODBCH,NODBCH,IODCLS,
     &                  NODCLS,ODTRI,IELCTR,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
      PARAMETER (D0 = 0.0D0)
      LOGICAL ODTRI
      DIMENSION NACTOD(NODCLS), IODCLS(NODCLS),
     &          SCREEN(NODBCH), RODBCH(NODBCH), IODBCH(NODBCH,*)
#include <aobtch.h>
C
      IF (IELCTR .EQ. 1) THEN
         IA = 1
         IB = 2
      ELSE
         IA = 3
         IB = 4
      END IF
C
      ISTART = 1
      NACTIV = 0
      DO 100 I = 1, NODCLS
         NODSAB = 0
         NODTAB = IODCLS(I)
         IF (.NOT.ODTRI) THEN
            DO 200 J = ISTART, ISTART + NODTAB - 1
               IF(ACTVBT(IODBCH(J,2),IA).AND.ACTVBT(IODBCH(J,3),IB))THEN
                  SCREEN(J) = RODBCH(J)
                  NODSAB = NODSAB + 1
               ELSE
                  SCREEN(J) = D0
               END IF
  200       CONTINUE
         ELSE
            DO 300 J = ISTART, ISTART + NODTAB - 1
               IF(
     &          (ACTVBT(IODBCH(J,2),IA).AND.ACTVBT(IODBCH(J,3),IB)) .OR.
     &          (ACTVBT(IODBCH(J,3),IA).AND.ACTVBT(IODBCH(J,2),IB)))THEN
                  SCREEN(J) = RODBCH(J)
                  NODSAB = NODSAB + 1
               ELSE
                  SCREEN(J) = D0
               END IF
  300       CONTINUE
         END IF
         NACTOD(I) = NODSAB
         NACTIV = NACTIV + NODSAB
         ISTART = ISTART + NODTAB
  100 CONTINUE
      NTOTAL = ISTART - 1
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Output from ERIACT',-1)
         WRITE (LUPRI,'(1X,A, I5)') ' Electron: ',IELCTR
         WRITE (LUPRI,'(1X,A,3I5)') ' NODCLS, NTOTAL, NACTIV ',
     &                                NODCLS, NTOTAL, NACTIV
         WRITE (LUPRI,'(/,5X,A)') '     class     total    active'
         WRITE (LUPRI,'(5X,A,/)') '------------------------------'
         DO 400 I = 1, NODCLS
            WRITE (LUPRI,'(5X,3I10)') I, IODCLS(I), NACTOD(I)
  400    CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck pickao */
      SUBROUTINE PICKAO(IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <aovec.h>
#include <mxcent.h>
#include <erisel.h>
#include <aobtch.h>
C
      DO 100 J = 1, 4
      DO 100 I = 1, NAOBCH
         ACTVBT(I,J) = .TRUE.
  100 CONTINUE
C
C     First orbital
C
      DO 200 J = 1, 4
         NSET = NSELCT(J)
         IF (NSET .GT. 0) THEN
            DO 300 I = 1, NAOBCH
               ACTVBT(I,J) = .FALSE.
  300       CONTINUE
            DO 400 I = 1, NSET
               ACTVBT(NACTAO(I,J),J) = .TRUE.
  400       CONTINUE
         END IF
  200 CONTINUE
C
      IF (IPRINT .GT. 3) THEN
         CALL HEADER('Output from PICKAO',-1)
         IF (NSELCT(1)+NSELCT(2)+NSELCT(3)+NSELCT(4).GT.0) THEN
            WRITE (LUPRI,'(1X,A,4I5,/)')
     &         ' Number of active AO batches for each electron:',
     &         (NSELCT(I),I=1,4)
            WRITE (LUPRI,'(1X,A,/)') ' ACTVBT in PICKAO:'
            DO 500 I = 1, NAOBCH
               WRITE (LUPRI,'(10X,I5,5X,4L5))') I, (ACTVBT(I,J),J=1,4)
  500       CONTINUE
         ELSE
            WRITE (LUPRI,'(1X,A)') ' All AO batches active.'
         END IF
      END IF
C
      RETURN
      END
C  /* Deck eridsi */
      SUBROUTINE ERIDSI(INDXBT,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <aovec.h>
#include <mxcent.h>
C
      LOGICAL FIRST
      DIMENSION INDXBT(MXSHEL*MXCONT,0:7)
C
#include <cbieri.h>
#include <eridst.h>
#include <erisel.h>
#include <aobtch.h>
#include <symmet.h>
C
#include <ibtfun.h>
C
      IF (DODIST) THEN
         NDISTR = 0
         NACDST = 0
         FIRST  = .TRUE.
         DO 100 I = 1, NAOBCH
            IF (ACTVBT(I,1)) THEN
               NDISTR = NDISTR + NORBBT(I)
               IF (NDISTR .GT. MXDIST) THEN
                  WRITE (LUPRI,'(/,1X,A,/,A,2I5,/,A)')
     &               ' Too many distributions in ERIDSI:',
     &               ' NDISTR, MXDIST =', NDISTR,MXDIST,
     &               ' Calculation stopped.'
                  CALL QUIT('Error in ERIDSI')
               END IF
               NACDST = NACDST + 1
               IACDST(NACDST) = I

               IF (FIRST) THEN
                  NHKDST = NHKTBT(I)
                  KHKDST = KHKTBT(I)
                  MLTDST = MULTBT(I)
C                 NPRDST = NPRFBT(I)
                  ISTABL = ISTBBT(I)
                  FIRST = .FALSE.
               ELSE
                  IF (
     &                   KHKTBT(I).NE.KHKDST .OR.
     &                   MULTBT(I).NE.MLTDST .OR.
C    &                   NPRFBT(I).NE.NPRDST .OR.
     &                   ISTBBT(I).NE.ISTABL
     &                ) THEN
                     WRITE (LUPRI,'(/,1X,A,/,A,2I5,/A,2I5,/,A)')
     &                  ' Inconsistent specification of distributions:',
     &                  ' KHKDST, KHKTBT(I) :', KHKDST, KHKTBT(I),
     &                  ' MLTDST, MULTBT(I) :', MLTDST, MULTBT(I),
     &                  ' Calculation stopped.'
                     CALL QUIT('Error in ERIDSI')
                  END IF
               END IF
            END IF
  100    CONTINUE
C
         IDIST = 0
         DO 200 IREP = 0, MAXREP
            IA = 0
            DO 300 ICMP = 1, KHKDST
            IF(IBTAND(ISTABL,IBTXOR(IREP,ISYMAO(NHKDST,ICMP))).EQ.0)THEN
               IA = IA + 1
               DO 400 J = 1, NACDST
                  IDIST = IDIST + 1
                  INDA = INDXBT(KNDXBT(IACDST(J)),IREP)+IA-1
                  INDDST(IDIST) = INDA
                  INDXDS(INDA) = IDIST
  400          CONTINUE
            END IF
  300       CONTINUE
  200    CONTINUE
C
      ELSE
         NDISTR = 1
         NACDST = 0
      END IF
C
      IF (IPRINT .GT. 1) THEN
         CALL HEADER('Output from ERIDSI',-1)
         WRITE (LUPRI,'(/,1X,A,I5)')
     &      ' Number of active AO batches: ', NACDST
         WRITE (LUPRI,'(1X,A,8I5,/,(30X,8I5))')
     &      ' Active AO batches:           ',(IACDST(I),I=1,NACDST)
         WRITE (LUPRI,'(/,1X,A,I5)')
     &      ' Number of distributions:     ', NDISTR
         WRITE (LUPRI,'(1X,A,8I5,/,(30X,8I5))')
     &      ' Distributions:              ',(INDDST(I),I=1,NDISTR)
         WRITE (LUPRI,'(/,1X,A,2I5)')
     &      ' Number of ang. components and multiplicity of batches:',
     &      KHKDST, MLTDST
      END IF
C
      RETURN
      END
C  /* Deck getdst */
      SUBROUTINE GETDST(ICLX,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <eridst.h>
C
      LOGICAL AOBTGT
C
#include <cbieri.h>
#include <erisel.h>
#include <distcl.h>
#include <aobtch.h>
C
      IF (ICLX.EQ.0) THEN
C
C        MAXCML
C
         MAXCML = 0
         DO 100 I = 1, NAOBCH
            MAXCML = MAX(MAXCML,NORBBT(KAOSRT(I)))
  100    CONTINUE
C
         IF (MAXCML .GT. MXDIST) THEN
            WRITE (LUPRI,'(2X,A,/,2X,A,2I5,/,2X,A)')
     &         'MAXCML exceeds MXDIST in GETDST',
     &         'MAXCML, MXDIST:', MAXCML, MXDIST,
     &         'Calculation cannot proceed.'
            CALL QUIT('MXDIST too small in GETDST')
         END IF
C
         IF (MAXCML .GT. MAXDST) THEN
            WRITE (LUPRI,'(2X,A,/,2X,A,2I5,/,2X,A)')
     &         'MAXCML exceeds MAXDST in GETDST',
     &         'MAXCML, MAXDST:', MAXCML, MAXDST,
     &         'Calculation cannot proceed.'
            CALL QUIT('MAXDST too small in GETDST')
         END IF
C
         MXBTCH = MAXDST
C
C        Classes of AO batches that can be calculated in
C        the same distribution call
C
         ICLASS = 1
         NCLASS(ICLASS) = 1
         IOLD = KAOSRT(1)
         MCLASS(ICLASS) = NORBBT(IOLD)
         KCLASS(ICLASS) = 1
         KLAOBT(1) = 1
         DO 200 I = 2, NAOBCH
            INEW = KAOSRT(I)
            IF (AOBTGT(INEW,IOLD)) THEN
               ICLASS = ICLASS + 1
               NCLASS(ICLASS) = 1
               MCLASS(ICLASS) = NORBBT(INEW)
               KCLASS(ICLASS) = I
               IOLD = INEW
            ELSE
               NCLASS(ICLASS) = NCLASS(ICLASS) + 1
            END IF
            KLAOBT(I) = ICLASS
  200    CONTINUE
         NTCLAS = ICLASS
C
         IF (IPRINT .GT. 2) THEN
            WRITE (LUPRI,'(/,2X,A,I5)')
     &         'Number of classes in GETDST:  ',NTCLAS
            WRITE (LUPRI,'(2X,A,9I5,/,(30X,9I5))')
     &         'Number of batches in classes: ',
     &         (NCLASS(I),I=1,NTCLAS)
            WRITE (LUPRI,'(2X,A,9I5,/,(30X,9I5),/)')
     &         'Multiplicity of batches:      ',
     &         (MCLASS(I),I=1,NTCLAS)
            WRITE (LUPRI,'(2X,A,9I5,/,(30X,9I5),/)')
     &         'Start adress of classes:      ',
     &         (KCLASS(I),I=1,NTCLAS)
            WRITE (LUPRI,'(2X,A,9I5,/,(30X,9I5),/)')
     &         'Classes of AO batches:        ',
     &         (KLAOBT(I),I=1,NAOBCH)
         END IF
C
C        Determine number of calls to integral program
C
         IPRVDS = 0
         ICALL = 0
         DO 300 I = 1, NAOBCH
            IFIRST = IPRVDS + 1
            ICLASS = KLAOBT(IFIRST)
            MULTPL = MCLASS(ICLASS)
            IBTMAX = MXBTCH/MULTPL
C
            IF (ICLASS .EQ. NTCLAS) THEN
               ILAST = MIN(IFIRST + IBTMAX - 1, NAOBCH)
            ELSE
               JCLASS = KLAOBT(IFIRST + IBTMAX - 1)
               IF (ICLASS .EQ. JCLASS) THEN
                  ILAST = IFIRST + IBTMAX - 1
               ELSE
                  ILAST = KCLASS(ICLASS + 1) - 1
               END IF
            END IF
            NACDST = ILAST - IFIRST + 1
C
            ICALL = ICALL + 1
            ICLFRS(ICALL) = IFIRST
            ICLLST(ICALL) = ILAST
            ICLBCH(ICALL) = NACDST
            ICLDST(ICALL) = MULTPL*NACDST
C
            IPRVDS = ILAST
C
            IF (ILAST .EQ. NAOBCH) GO TO 399
C
  300    CONTINUE
  399    CONTINUE
C
         MXCALL = ICALL
         IF (IPRINT .GT. 1) THEN
            WRITE (LUPRI,'(/,2X,A,I5)')
     &         'Number of calls to integral program ',MXCALL
            IF (IPRINT .GT. 2) THEN
               DO 400 I = 1, MXCALL
                  WRITE (LUPRI,'(5X,5I5)')
     &                I, ICLFRS(I), ICLLST(I), ICLBCH(I), ICLDST(I)
  400          CONTINUE
            END IF
         END IF
      ELSE
         NSELCT(1) = ICLBCH(ICLX)
         NSELCT(2) = 0
         NSELCT(3) = 0
         NSELCT(4) = 0
         IFIRST = ICLFRS(ICLX)
         DO 600 I = 1, ICLBCH(ICLX)
            NACTAO(I,1) = KAOSRT(IFIRST + I - 1)
  600    CONTINUE
C
      END IF
C
      RETURN
      END
C  /* Deck er2dis */
      SUBROUTINE ER2DIS(WORK,LWORK)
C
#include <implicit.h>
#include <iratdef.h>
#include <priunit.h>
C
#include <maxorb.h>
#include <maxash.h>
#include <mxorb.h>
#include <aovec.h>
C
      DIMENSION INDEXA(MXCORB), WORK(LWORK)
C
#include <cbieri.h>
#include <distcl.h>
C
      KCCFB1 = 1
      KINDXB = KCCFB1 + MXPRIM*MXCONT
      KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
      LWRK1  = LWORK  - KEND1
      IF (KEND1 .GT. LWORK) CALL STOPIT('ERIDIS',' ',KEND1,LWORK)
      CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     &            KODPP1,KODPP2,KRDPP1,KRDPP2,
     &            KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     &            WORK(KEND1),LWRK1,IPRERI)
      KEND1 = KFREE
      LWRK1 = LFREE
C
      DO 100 IDIST = 1, MXCALL
         WRITE (LUPRI,'(A,I5)')
     &      ' Calling ERIDI2 for distribution #',IDIST
         CALL TIMER('START ',TIMSTR,TIMEND)
         CALL ERIDI2(IDIST,INDEXA,NUMDIS,
     &               WORK(KODCL1),WORK(KODCL2),
     &               WORK(KODBC1),WORK(KODBC2),
     &               WORK(KRDBC1),WORK(KRDBC2),
     &               WORK(KODPP1),WORK(KODPP2),
     &               WORK(KRDPP1),WORK(KRDPP2),
     &               WORK(KCCFB1),WORK(KINDXB),
     &               WORK(KEND1), LWRK1,IPRERI)
         WRITE (LUPRI,'(2X,A,I5)') ' NUMDIS from this call:',NUMDIS
         WRITE (LUPRI,'(1X,A,8I5,/,(30X,8I5))')
     &      ' Distributions:              ',(INDEXA(I),I=1,NUMDIS)
         CALL TIMER('ERIDI2',TIMSTR,TIMEND)
  100 CONTINUE
C
      RETURN
      END
