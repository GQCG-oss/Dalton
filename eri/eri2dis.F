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
C  /* Deck eridi1 */
      SUBROUTINE ERIDI1(KODCL1,KODCL2,
     &                  KODBC1,KODBC2,KRDBC1,KRDBC2,
     &                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     &                  KFREE,LFREE,KEND,CCFBT,INDXBT,
     &                  WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxaqn.h>
#include <maxorb.h>
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
#include <iratdef.h>
#include <inftap.h>
C
      IPRINT = IPRERI
      THRSH  = MAX(THRS,1.00D-15)
C
C
      DODIST = .TRUE.
      PMSAB  = .FALSE.
      PMS12  = .FALSE.
C
C
C     Memory
C
      MEMOK  = .TRUE.
      MEMADD = 0
      MODAB  = 0
      MODCD  = 0
C
      WRTINT = .TRUE.
      FCKINT = .FALSE.
C
C     AO batches
C     ==========
C
      CALL SETAOB(CCFBT,INDXBT,WORK,LWORK,IPRINT)
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
      CALL GETDST(ICALL,ICALL,IPRINT)
C
      RETURN
      END
C  /* Deck eridi2 */
      SUBROUTINE ERIDI2(ICALL,INDEXA,NUMDIS,NGDER,NBDER,
     &                  IODCL1,IODCL2,
     &                  IODBC1,IODBC2,RODBC1,RODBC2,
     &                  IODPP1,IODPP2,RODPP1,RODPP2,
     &                  CCFBT,INDXBT,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <eridst.h>
#include <dummy.h>
      CHARACTER*8 FAODER
      DIMENSION IODCL1(NODCL1,NITCL), IODCL2(NODCL2,NITCL),
     &          IODBC1(NODBC1,NITBC), IODBC2(NODBC2,NITBC),
     &          RODBC1(NODBC1,NRTBC), RODBC2(NODBC2,NRTBC),
     &          IODPP1(NODPP1,NITPP), IODPP2(NODPP2,NITPP),
     &          RODPP1(NODPP1,NRTPP), RODPP2(NODPP2,NRTPP),
     &          INDEXA(*), CCFBT(*), INDXBT(*), WORK(LWORK)
      LOGICAL   OLDDX
#include <ccom.h>
#include <cbieri.h>
#include <ericom.h>
#include <erithr.h>
#include <erimem.h>
#include <eribuf.h>
#include <aobtch.h>
#include <odbtch.h>
#include <nuclei.h>
#include <chrnos.h>
#include <odclss.h>
#include <symmet.h>
#include <iratdef.h>
#include <inftap.h>
C
      IPRINT = IPRERI
C
      GDER   = NGDER .GT. 0 
      BDER   = NBDER .GT. 0 
      IF (GDER .OR. BDER) THEN
         UNDIFF = .TRUE.
      ELSE
         UNDIFF = .FALSE.
      END IF
      WRTINT = .TRUE.
      FCKINT = .FALSE.
      EXPERI = .FALSE.
      CCRUN  = .TRUE.
      NPDIMA = 0
      NPDIMB = 0
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
      IF (LBFINP .EQ. 600) LBFINP = 20 000
C
C     Open files
C
      IF (LUINTR.LT.0) CALL GPOPEN(LUINTR,'AOTWODIS','UNKNOWN',' ',
     &     'UNFORMATTED',IDUMMY,.FALSE.)
#ifdef VAR_INT64
      LRECL = LBFINP + LBFINP + 1
#elif defined (SYS_NEC)
      LRECL = LBFINP + LBFINP/IRAT + 1
#else
      IF (NBASIS .LE. 255) THEN
         LRECL = 2*LBFINP + LBFINP + 1
      ELSE
         LRECL = 2*LBFINP + 2*LBFINP + 1
      END IF
#endif
      NSCOOR = 0
      IF (GDER) NSCOOR = 3*NUCDEP
      IF (BDER) NSCOOR = 6
      DO ISCOOR = 0, NSCOOR
         FAODER = 'AO2DIS'//CHRNOS(ISCOOR/10)
     &                    //CHRNOS(MOD(ISCOOR,10))
         CALL GPOPEN(LUAORC(ISCOOR),FAODER,'UNKNOWN','DIRECT',
     &        'UNFORMATTED',LRECL,OLDDX)
      END DO
C
      CALL GETDST(ICALL,0,IPRINT)
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
!fox
!fox  25-04-2003
!fox  CALL ICOPY(NUMDIS,INDDST,1,INDEXA,1)
      IF (NUMDIS .GT. 0) CALL ICOPY(NUMDIS,INDDST,1,INDEXA,1)
!fox
      IF ( IPRINT .GT. 2 ) THEN
         CALL HEADER('Output from ERIDI2',2)
         WRITE (LUPRI,'(2X,I3,A,8I5/,(40X,8I5))')
     &         NDISTR,' distributions in this ERI call:',
     &         (INDEXA(I),I=1,NDISTR)
      ENDIF
C
C     Calculate integrals
C     ===================
C
      IF (.NOT. INTSKP) THEN
         CALL TIMER('START ',TIMSTR,TIMEND)
         CALL ODCDRV(IODCL1,IODCL2,
     &               IODBC1,IODBC2,RODBC1,RODBC2,
     &               IODPP1,IODPP2,RODPP1,RODPP2,
     &               DUMMY,DUMMY,IDUMMY,IDUMMY,DUMMY,IDUMMY,CCFBT,
     &               INDXBT,WORK,LWORK,IPRINT)
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
      DO ISCOOR = 0, NSCOOR
         CALL GPCLOSE(LUAORC(ISCOOR),'KEEP')
      END DO
C
      RETURN
      END
C  /* Deck eridid */
      SUBROUTINE ERIDID(ICALL1,ICALL2,D2MAT,ID2MAT,NDIMA,NDIMB,
     &                  IODCL1,IODCL2,IODBC1,IODBC2,RODBC1,RODBC2,
     &                  IODPP1,IODPP2,RODPP1,RODPP2,
     &                  CCFBT,INDXBT,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <eridst.h>
#include <dummy.h>
      DIMENSION IODCL1(NODCL1,NITCL), IODCL2(NODCL2,NITCL),
     &          IODBC1(NODBC1,NITBC), IODBC2(NODBC2,NITBC),
     &          RODBC1(NODBC1,NRTBC), RODBC2(NODBC2,NRTBC),
     &          IODPP1(NODPP1,NITPP), IODPP2(NODPP2,NITPP),
     &          RODPP1(NODPP1,NRTPP), RODPP2(NODPP2,NRTPP),
     &          CCFBT(*), INDXBT(*), WORK(LWORK)
      DIMENSION D2MAT(*), ID2MAT(*)
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
      GDER   = .TRUE.
      BDER   = .FALSE.
      WRTINT = .FALSE.
      FCKINT = .FALSE.
      EXPERI = .TRUE.
c
c     Cristian, change it
      UNDIFF = .TRUE.
c
      CCRUN  = .TRUE.
      NPDIMA = NDIMA
      NPDIMB = NDIMB
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
      CALL GETDST(ICALL1,ICALL2,IPRINT)
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
C     Calculate integrals
C     ===================
C
      IF (.NOT.INTSKP) THEN
         CALL TIMER('START ',TIMSTR,TIMEND)
         CALL ODCDRV(IODCL1,IODCL2,
     &               IODBC1,IODBC2,RODBC1,RODBC2,
     &               IODPP1,IODPP2,RODPP1,RODPP2,
     &               DUMMY,DUMMY,IDUMMY,IDUMMY,D2MAT,ID2MAT,CCFBT,
     &               INDXBT,WORK,LWORK,IPRINT)
         IF ( IPRINT .GT. 2) CALL TIMER('HERMIT',TIMSTR,TIMEND)
         CALL FLSHFO(LUPRI)
C
C        Error message in case of insufficient memory
C
         IF (.NOT.MEMOK) THEN
            WRITE (LUPRI,'(//,1X,A,3(/,1X,A,I10))')
     &         ' Not enough memory for this run of ERIDID.',
     &         ' Available memory in ERIDID:',LWORK,
     &         ' Required memory for ERIDID:',LWORK + MEMADD,
     &         ' Increase memory (LWORK) by:',MEMADD
            WRITE (LUPRI,'(/,1X,A,2I5)')
     &         ' Memory requirements largest for OD classes :',
     &           MODAB,MODCD
            CALL QUIT('Insufficient memory in ERIDID.')
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
#include <maxorb.h>
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
#include <maxorb.h>
#include <aovec.h>
#include <mxcent.h>
#include <cbirea.h>
#include <r12int.h>
#include <shells.h>
#include <erisel.h>
#include <aobtch.h>
C
      DO J = 1, 4
      DO I = 1, NAOBCH
         IF (KCLSBT(I).EQ.1) THEN
             ACTVBT(I,J) = .TRUE.
          ELSE
             ACTVBT(I,J) = .FALSE.
          ENDIF
      END DO
      END DO
C
C     First orbital
C
      DO J = 1, 4
         NSET = NSELCT(J)
         IF (NSET .GT. 0) THEN
            DO I = 1, NAOBCH
               ACTVBT(I,J) = .FALSE.
            END DO
            DO I = 1, NSET
               ACTVBT(NACTAO(I,J),J) = .TRUE.
            END DO
         END IF
      END DO
C
      IF (LMULBS .AND. LOOPDP) THEN
         DO I = 1, NAOBCH
            IF (MBIDBT(I) .GT. 1) THEN
               DO J = 2, 4
                  ACTVBT(I,J) = .FALSE.
               END DO
            END IF
            IF (MBIDBT(I) .EQ. 1) THEN
               ACTVBT(I,1) = .FALSE.
            END IF
         END DO
      ELSE IF (LMULBS .AND. .NOT. R12TRA) THEN
C        Do not compute integrals that involve auxiliary
C        basis functions (WK/UniKA/04-11-2002).
         DO I = 1, NAOBCH
            IF (MBIDBT(I) .GT. 1) THEN
               DO J = 1, 4
                  ACTVBT(I,J) = .FALSE.
               END DO
            END IF
         END DO
      END IF
C
      IF (IPRINT .GT. 3) THEN
         CALL HEADER('Output from PICKAO',-1)
         IF (NSELCT(1)+NSELCT(2)+NSELCT(3)+NSELCT(4).GT.0) THEN
            WRITE (LUPRI,'(A,4I5/A/)')
     &         '  Number of active AO batches for each electron:',
     &         (NSELCT(I),I=1,4),
     &         '  ACTVBT(1:NAOBCH,1:4) in PICKAO:'
            DO I = 1, NAOBCH
               WRITE (LUPRI,'(10X,I5,5X,4L5)') I, (ACTVBT(I,J),J=1,4)
            END DO
         ELSE
            WRITE (LUPRI,'(A)') '  All AO batches active.'
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
#include <maxorb.h>
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
         NCTDST = 0
         NDISTR = 0
         NACDST = 0
         FIRST  = .TRUE.
         DO I = 1, NAOBCH
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
C                 NCTDST = NCTFBT(I) ! tew: number of contractns not always constant
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
            else
c             print * ,"skipping i = ", i
            END IF
         END DO
C
C tew: allow for different numbers of contractions for each J
C tew: NB// The order of IDIST is not the same as before.
         IDIST = 0
         IF (NACDST .GT. 0) THEN
           DO IREP = 0, MAXREP
             IA = 0
             DO ICMP = 1, KHKDST
               IF(IBTAND(ISTABL,IBTXOR(IREP,ISYMAO(NHKDST,ICMP)))
     &          .EQ. 0) THEN
                 IA = IA + 1
                 DO J = 1, NACDST
                   NCTDST=NCTFBT(IACDST(J))
                   DO ICNT = 1, NCTDST
                     IDIST = IDIST + 1
                     INDA = INDXBT(KNDXBT(IACDST(J))-1+ICNT,IREP)+IA-1
                     INDDST(IDIST) = INDA
                     INDXDS(INDA) = IDIST
                   END DO
                 END DO
               END IF
             END DO
           END DO
         END IF
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
      SUBROUTINE GETDST(ICLX1,ICLX2,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
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
C
C     ICLX .EQ. 0: initialization call
C
      ICLX = ICLX1
      IF (ICLX.EQ.0) THEN
C
C        MAXCML
C
         MAXCML = 0
         DO I = 1, NAOBCH
            MAXCML = MAX(MAXCML,NORBBT(KAOSRT(I)))
         END DO
C
C        MXDIST: value used in parameter statement for 
C                common blocks etc. 
C        MAXDST: input parameter with default
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
         DO I = 2, NAOBCH
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
         END DO
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
         DO I = 1, NAOBCH
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
            IF (ILAST .EQ. NAOBCH) GO TO 100
C
         END DO
  100    CONTINUE
C
         MXCALL = ICALL
         IF (IPRINT .GT. 1) THEN
            WRITE (LUPRI,'(/,2X,A,I5)')
     &         'Number of calls to integral program ',MXCALL
            IF (IPRINT .GT. 2) THEN
               DO I = 1, MXCALL
                  WRITE (LUPRI,'(5X,5I5)')
     &                I, ICLFRS(I), ICLLST(I), ICLBCH(I), ICLDST(I)
               END DO
            END IF
         END IF
      ELSE IF (ICLX2.EQ.0) THEN
         NSELCT(1) = ICLBCH(ICLX)
         NSELCT(2) = 0
         NSELCT(3) = 0
         NSELCT(4) = 0
         IFIRST = ICLFRS(ICLX)
         DO I = 1, ICLBCH(ICLX)
            NACTAO(I,1) = KAOSRT(IFIRST + I - 1)
         END DO
      ELSE 
         NSELCT(1) = ICLBCH(ICLX1)
         NSELCT(2) = ICLBCH(ICLX2)
         NSELCT(3) = 0
         NSELCT(4) = 0
         IFIRST = ICLFRS(ICLX1)
         DO I = 1, ICLBCH(ICLX1)
            NACTAO(I,1) = KAOSRT(IFIRST + I - 1)
         END DO
         IFIRST = ICLFRS(ICLX2)
         DO I = 1, ICLBCH(ICLX2)
            NACTAO(I,2) = KAOSRT(IFIRST + I - 1)
         END DO
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
      DO IDIST = 1, MXCALL
         WRITE (LUPRI,'(A,I5)')
     &      ' Calling ERIDI2 for distribution #',IDIST
         CALL TIMER('START ',TIMSTR,TIMEND)
         CALL ERIDI2(IDIST,INDEXA,NUMDIS,0,0,
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
      END DO
C
      RETURN
      END
C  /* Deck eriidx */
      SUBROUTINE ERIIDX(ICALL,INDEXA,NUMDIS,INDXBT,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <eridst.h>
#include <dummy.h>
      DIMENSION INDEXA(*), INDXBT(*)
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
      CALL GETDST(ICALL,0,IPRINT)
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
      IF (NUMDIS .GT. 0) CALL ICOPY(NUMDIS,INDDST,1,INDEXA,1)
      IF ( IPRINT .GT. 2 ) THEN
         CALL HEADER('Output from ERIIDX',2)
         WRITE (LUPRI,'(2X,I3,A,8I5/,(40X,8I5))')
     &         NDISTR,' distributions in this ERI call:',
     &         (INDEXA(I),I=1,NDISTR)
      ENDIF
C
      RETURN
      END
