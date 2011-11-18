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
C  /* Deck er2inp */
      SUBROUTINE ER2INP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (NTABLE = 35)
C
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
#include <cbiher.h>
#include <cbieri.h>
#include <erisel.h>
#include <eribuf.h>
#include <gnrinf.h>
      DATA TABLE /'.SKIP  ', '.PRINT ', '.TIME  ', '.RETURN', '.AOBTCH',
     &            '.NOPSAB', '.NOPSCD', '.NOPS12', '.OFFCNT', '.NO12GS',
     &            '.NEWCR1', '.NONCAN', '.WRITEA', '.NOWRIT', '.GABAB ',
     &            '.EXTPRI', '.SELCT1', '.SELCT2', '.SELCT3', '.SELCT4',
     &            '.INTPRI', '.DISTRI', '.INTSKI', '.BUFFER', '.MAXDIS',
     &            '.DISTST', '.MXBCH ', '.NSPMAX', '.NOLOCS', '.GRDZER',
     &            '.XXXXXX', '.DOERIP', '.NOSCRE', '.GENCON', '.NCLERI'/
C
      NEWDEF = WORD .EQ. '*ER2INT'
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
     &                      18,19,20,21,22,23,24,25,26,27,28,29,30,31,
     &                      32,33,34,35),I
                  END IF
               END DO
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     &            '" not recognized in ER2INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in ER2INP.')
    1          CONTINUE
                  RUNERI = .FALSE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRERI, IPRNT1, IPRNT2
                  IPRSUM = IPRNT1 + IPRNT2
                  IF (IPRERI .EQ. IPRDEF .AND. IPRSUM .EQ. 0) THEN
                     ICHANG = ICHANG - 1
                  END IF
               GO TO 100
    3          CONTINUE
                  TIMERI = .TRUE.
               GO TO 100
    4          CONTINUE
                  RTNERI = .TRUE.
               GO TO 100
    5          CONTINUE
                  READ (LUCMD,*) IAOBCH
               GO TO 100
    6          CONTINUE
                  PMSAB = .FALSE.
               GO TO 100
    7          CONTINUE
                  PMSCD = .FALSE.
               GO TO 100
    8          CONTINUE
                  PMS12 = .FALSE.
               GO TO 100
    9          CONTINUE
                  OFFCNT = .TRUE.
               GO TO 100
   10          CONTINUE
                  DIASRT = .FALSE.
               GO TO 100
   11          CONTINUE
                  OLDCR1 = .FALSE.
               GO TO 100
   12          CONTINUE
                  CANIND = .FALSE.
               GO TO 100
   13          CONTINUE
                  WRTSCR = .FALSE.
               GO TO 100
   14          CONTINUE
                  NOWRIT = .TRUE.
               GO TO 100
   15          CONTINUE
                  DGABAB = .TRUE.
               GO TO 100
   16          CONTINUE
                  EXTPRI = .TRUE.
                  READ (LUCMD,*) IPROD1, IPROD2
               GO TO 100
   17          CONTINUE
                  READ (LUCMD,*) NSELCT(1)
                  READ (LUCMD,*) (NACTAO(I,1),I=1,NSELCT(1))
               GO TO 100
   18          CONTINUE
                  READ (LUCMD,*) NSELCT(2)
                  READ (LUCMD,*) (NACTAO(I,2),I=1,NSELCT(2))
               GO TO 100
   19          CONTINUE
                  READ (LUCMD,*) NSELCT(3)
                  READ (LUCMD,*) (NACTAO(I,3),I=1,NSELCT(3))
               GO TO 100
   20          CONTINUE
                  READ (LUCMD,*) NSELCT(4)
                  READ (LUCMD,*) (NACTAO(I,4),I=1,NSELCT(4))
               GO TO 100
   21          CONTINUE
                  INTPRI = .TRUE.
               GO TO 100
   22          CONTINUE
                  DODIST = .TRUE.
               GO TO 100
   23          CONTINUE
                  INTSKP = .TRUE.
               GO TO 100
   24          CONTINUE
                  READ (LUCMD,*) LBFINP
                  LBFINP = MAX(1,LBFINP)
               GO TO 100
   25          CONTINUE
                  READ (LUCMD,*) MAXDST
               GO TO 100
   26          CONTINUE
                  DISTST = .TRUE.
               GO TO 100
   27          CONTINUE
                  READ (LUCMD,*) MXBCH
               GO TO 100
   28          CONTINUE
                  READ (LUCMD,*) NSPMAX
               GO TO 100
   29          CONTINUE
                  NOLOCS = .TRUE. 
               GO TO 100
   30          CONTINUE
                  GRDZER = .TRUE. 
               GO TO 100
   31          CONTINUE
               GO TO 100
   32          CONTINUE
                  DOERIP = .TRUE. 
               GO TO 100
   33          CONTINUE
                  COMPRS = .FALSE. 
               GO TO 100
   34          CONTINUE
                  GENCON_ERI = .TRUE. 
               GO TO 100
   35          CONTINUE
                  NCLERI = .TRUE. 
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in ER2INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in ER2INP.')
            END IF
      END IF
  300 CONTINUE
      IF (PMSAB .NEQV. PMSCD) PMS12 = .FALSE.
      IF (IPRERI .GT. 10) INTPRI = .TRUE.
      IF (ICHANG .EQ. 0) RETURN
      IF (NEWDEF) THEN
         CALL HEADER('Changes of defaults for ER2INP:',1)
         IF (.NOT.RUNERI .AND. .NOT. (DIRCAL .OR. PARCAL)) THEN
            WRITE (LUPRI,'(A)') ' No two-electron integrals calculated.'
         ELSE
            IF (IPRERI .NE. IPRDEF) WRITE (LUPRI,'(A,I5)')
     &         ' Print level in ERIINT:',IPRERI
            IF (IPRSUM .GT. 0) THEN
                 WRITE (LUPRI,'(A,2I3)')
     &             ' Extra output for the following ODs:',IPRNT1, IPRNT2
                IF (RTNERI) WRITE (LUPRI,'(A)')
     &               ' Program will exit ERIINT after these shells.'
            END IF
            IF (INTPRI) WRITE (LUPRI,'(A)')' Integrals will be printed.'
            IF (TIMERI) WRITE (LUPRI,'(A)')
     &         ' Detailed timing for two-electron integrals.'
            IF (INTSKP) WRITE (LUPRI,'(A)')
     &         ' Integrals (call to ODCDRV) skipped.'
            IF (DISTST) WRITE (LUPRI,'(A)')
     &         ' Distribution test.'
            IF (OFFCNT) WRITE (LUPRI,'(A)')
     &         ' No sorting of one- and two-center ODs.'
            IF (.NOT.DIASRT) WRITE (LUPRI,'(A)')
     &         ' No sorting of one- and two-GTO ODs.'
            IF (.NOT.PMSAB) WRITE (LUPRI,'(A)')
     &         ' No permutational symmetry for first electron.'
            IF (.NOT.PMSCD) WRITE (LUPRI,'(A)')
     &         ' No permutational symmetry for second electron.'
            IF (.NOT.PMSAB) WRITE (LUPRI,'(A)')
     &         ' No permutational symmetry between electron 1 and 2.'
            IF (IAOBCH.NE.0) WRITE (LUPRI,'(A,I4,A)')
     &         ' Only integrals with first index belonging to'//
     &         ' AO batch',IAOBCH,' calculated.'
            IF (.NOT.OLDCR1) WRITE (LUPRI,'(A)')
     &         ' New algorithm for Cartesians on first electron.'
            IF (.NOT.CANIND) WRITE (LUPRI,'(A)')
     &         ' Indices not sorted in canonical order.'
            IF (.NOT.WRTSCR) WRITE (LUPRI,'(A)')
     &         ' Integrals not screened before write.'
            IF (NOWRIT) WRITE (LUPRI,'(A)')
     &         ' Integrals will not be written to file.'
            IF (DGABAB) WRITE (LUPRI,'(A)')
     &         ' Only gabab integrals are calculated.' 
            IF (NOLOCS) WRITE (LUPRI,'(A)')
     &         ' No use of local symmetry.'
            IF (GRDZER) WRITE (LUPRI,'(A)')
     &         ' GRADEE is set to zero at each call to ERIGRD.'
            IF (DOERIP) WRITE (LUPRI,'(A)')
     &         ' ERIPRO used instead of TWOINT.'
            IF (EXTPRI) WRITE (LUPRI,'(A,2I5)')
     &         ' Extra high print for the following OD pair:',
     &           IPROD1, IPROD2
            IF (.NOT.COMPRS) WRITE (LUPRI,'(A)')
     &         ' No integral screening before write.'
            IF (GENCON_ERI) WRITE (LUPRI,'(A)')
     &         ' All contractions treated as general.' 
            IF (NCLERI) WRITE (LUPRI,'(A)')
     &         ' Only nonclassical integrals are calculated.' 
            IF (DODIST) THEN
               WRITE (LUPRI, '(A)')
     &         ' Distributions for electron 1.'
               IF (NSELCT(1).EQ.0) THEN
                  WRITE (LUPRI,'(/,A,/A)')
     &            ' Distributions specified but no active '//
     &            ' orbitals for electron. Use .SELCT1.',
     &            ' Calculation cannot proceed.'
                  CALL QUIT('Inconsistent input in ER2INP.')
               END IF
            END IF
            DO IELCTR = 1, 4
               IF (NSELCT(IELCTR).GT.0) THEN
                  WRITE (LUPRI,'(/A,I3,A,I2,A,5I5,/,(43X,5I5))')
     &            ' Only ',NSELCT(IELCTR),
     &            ' AO batches active for electron',IELCTR,':',
     &            (NACTAO(I,IELCTR),I=1,NSELCT(IELCTR))
               END IF
            END DO
            IF (LBFINP .NE. 600) THEN
               WRITE (LUPRI,'(A,I5)')
     &            ' Non-default buffer length:',LBFINP
            END IF
            IF (MAXDST .NE. MAXDSD) THEN
               WRITE (LUPRI,'(A,I5)')
     &            ' Maximum number of distribution in each call:',
     &              MAXDST
            END IF
            IF (MXBCH .NE. MXBCH0) THEN
               WRITE (LUPRI,'(A,I5)')
     &            ' MXBCH for this run:', MXBCH
            END IF
            IF (NSPMAX .NE. 8**3) THEN
               NSPMAX = 2**INT(LOG(REAL(NSPMAX))/LOG(REAL(2)))
               NSPMAX = MIN(MAX(1,NSPMAX),8**3)
               WRITE (LUPRI,'(A,I5)')
     &            ' NSPMAX for this run:', NSPMAX
            END IF
         END IF
      END IF
      RETURN
      END
C  /* Deck er2ini */
      SUBROUTINE ER2INI
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <cbiher.h>
#include <cbieri.h>
#include <erisel.h>
#include <eribuf.h>
#include <ericom.h>
#include <gnrinf.h>
C
C     Initialize /CBIERI/
C
      IF (PARCAL) RUNERI = .FALSE.
      IPRERI = IPRDEF
      IPRNT1 = 0
      IPRNT2 = 0
      RTNERI = .FALSE.
      TIMERI = .FALSE.
      IAOBCH = 0
      PMSAB  = .TRUE.
      PMSCD  = .TRUE.
      PMS12  = .TRUE.
      OFFCNT = .FALSE.
      DIASRT = .TRUE.
      OLDCR1 = .TRUE.
      CANIND = .TRUE.
      WRTSCR = .TRUE.
      DGABAB = .FALSE.
      NOWRIT = .FALSE.
      EXTPRI = .FALSE.
      INTPRI = .FALSE.
      DODIST = .FALSE.
      INTSKP = .FALSE.
      NOLOCS = .FALSE.
      GRDZER = .FALSE.
      COMPRS = .TRUE.
      GDER   = .FALSE.
      BDER   = .FALSE.
      EXPERI = .FALSE.
      LBFINP = 600
      MAXDSD = 40
      MAXDST = MAXDSD
      DISTST = .FALSE.
      MXBCH0 = 1000000000
      MXBCH  = MXBCH0
      NSPMAX = 8**3 
      DOERIP = .FALSE.
      DO I = 1, 4
         NSELCT(I) = 0
      END DO
      GENCON_ERI = .FALSE.
      NCLERI = .FALSE.

      RETURN
      END
C  /* Deck er2int */
      SUBROUTINE ER2INT(CCFBT,INDXBT,WORK,LWORK)
#include <implicit.h>
#include <iratdef.h>
#include <priunit.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <dummy.h>
      DIMENSION CCFBT(*), INDXBT(*), WORK(LWORK), CSMAT(NBAST,NBAST)
      LOGICAL WRTINT_SAVE
#include <ccom.h>
#include <cbieri.h>
#include <ericom.h>
#include <erithr.h>
#include <erimem.h>
#include <aobtch.h>
#include <odbtch.h>
#include <symmet.h>
#include <inforb.h>
C
      CALL QENTER('ER2INT')
C
      IPRINT = IPRERI
      THRSH  = MAX(THRS,1.00D-15)
C
C     CALL ER2INI
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
     &            KFREE,LFREE,CCFBT,WORK,
     &            LWORK,IPRINT)
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
      KLAST = KFREE
      LWRK  = LFREE
C
C     Cauchy-Schwarz factors
C
      IF (.TRUE.) THEN
         DGABAB = .TRUE.
         WRTINT_SAVE = WRTINT
         WRTINT = .FALSE.
         CALL ODCDRV(WORK(KODCL1),WORK(KODCL2),
     &               WORK(KODBC1),WORK(KODBC2),
     &               WORK(KRDBC1),WORK(KRDBC2),
     &               WORK(KODPP1),WORK(KODPP2),
     &               WORK(KRDPP1),WORK(KRDPP2),
     &               CSMAT,DUMMY,1,IDUMMY,DUMMY,IDUMMY,
     &               CCFBT,INDXBT,WORK(KLAST),LWRK,IPRINT)
C
         IF (IPRINT.GT.5) THEN
            CALL HEADER('Cauchy-Schwarz matrix in ER2INT',-1)
            CALL OUTPUT(CSMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
         CALL CSMEXT(CSMAT)
         DGABAB = .FALSE.
         WRTINT = WRTINT_SAVE
      END IF
C
C     Calculate integrals
C     ===================
C
      IF (.NOT.INTSKP) THEN
         CALL ODCDRV(WORK(KODCL1),WORK(KODCL2),
     &               WORK(KODBC1),WORK(KODBC2),
     &               WORK(KRDBC1),WORK(KRDBC2),
     &               WORK(KODPP1),WORK(KODPP2),
     &               WORK(KRDPP1),WORK(KRDPP2),
     &               DUMMY,DUMMY,IDUMMY,IDUMMY,DUMMY,IDUMMY,
     &               CCFBT,INDXBT,WORK(KLAST),LWRK,IPRINT)
C
C        Error message in case of insufficient memory
C
         IF (.NOT.MEMOK) THEN
            WRITE (LUPRI,'(//,1X,A,3(/,1X,A,I10))')
     &         ' Not enough memory for this run of ER2INT.',
     &         ' Available memory in ER2INT:',LWORK,
     &         ' Required memory for ER2INT:',LWORK + MEMADD,
     &         ' Increase memory (LWORK) by:',MEMADD
            WRITE (LUPRI,'(/,1X,A,2I5)')
     &         ' Memory requirements largest for OD classes :',
     &           MODAB,MODCD
            CALL QUIT('Insufficient memory in ER2INT.')
         END IF
      END IF
C
      CALL QEXIT('ER2INT')
      RETURN
      END
C  /* Deck odcdrv */
      SUBROUTINE ODCDRV(IODCL1,IODCL2,
     &                  IODBC1,IODBC2,RODBC1,RODBC2,
     &                  IODPP1,IODPP2,RODPP1,RODPP2,
     &                  FMAT,DMAT,NDMT,IFCTYP,D2MAT,ID2MAT,
     &                  CCFBT,INDXBT,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
      DIMENSION IODCL1(NODCL1,NITCL), IODCL2(NODCL2,NITCL),
     &          IODBC1(NODBC1,NITBC), IODBC2(NODBC2,NITBC),
     &          RODBC1(NODBC1,NRTBC), RODBC2(NODBC2,NRTBC),
     &          IODPP1(NODPP1,NITPP), IODPP2(NODPP2,NITPP),
     &          RODPP1(NODPP1,NRTPP), RODPP2(NODPP2,NRTPP),
     &          FMAT(*), DMAT(*), IFCTYP(*), D2MAT(*), ID2MAT(*),
     &          CCFBT(*), INDXBT(*), WORK(LWORK)
#include <cbieri.h>
#include <odclss.h>
#include <ccom.h>
#include <ericom.h>
#include <eridst.h>
#include <r12int.h>
#include <hertop.h>
#include <eribuf.h>
#include <nuclei.h>
#include <aobtch.h>
#include <iratdef.h>
#include <infpar.h>
C
      MAXDER = 0
      IF (GDER) MAXDER = 1
      IF (BDER) MAXDER = 1
C
      JTOP  = 4*(NHTYP - 1) + MAXDER
C    
C     JTOP must be increased for [T1,r12] integrals (WK/UniKA/04-11-2002).
      IF (U12INT) JTOP = JTOP + 1
C     JTOP must be increased for integrals #5 (WK/UniKA/04-11-2002).
      IF (INTGAC .EQ. 5) JTOP = JTOP + 1
C     JTOP must be at least 1 for r12 integrals (WK/UniKA/04-11-2002).
      IF (R12INT) JTOP = MAX(JTOP,1)
C
      IF (BDER) JTOP = JTOP + 1
      JTOP3 = (JTOP + 1)**3
      NRTOP = (JTOP + 1)*(JTOP + 2)*(JTOP + 3)/6
C
      NHKMAX = NHTYP + MAXDER
      KCKMAX = NHKMAX*(NHKMAX + 1)/2
      KC2MAX = KCKMAX**2
      NRDER  = 3*MAXDER
C
C     Buffers and bits
C     ================
C
      LBUF  = LBFINP
      NBUFS = NDISTR
      IF (GDER .AND. WRTINT) THEN
         NPERTS = 3*NUCDEP
         NCORS = NPERTS 
         IF (UNDIFF) NCORS = NCORS + 1
      ELSE IF (BDER .AND. WRTINT) THEN
         NPERTS = 6
         NCORS = NPERTS
         IF (UNDIFF) NCORS = NCORS + 1
      ELSE
         NPERTS = 1
         NCORS = 1
      END IF
#if defined (VAR_INT64)
      NIBUF = 1
      NBITS = 16
#else
      IF (NBASE .LE. 255) THEN
         NIBUF = 1
         NBITS = 8
      ELSE
         NIBUF = 2
         NBITS = 16
      END IF
#endif
      IF (DGABAB) THEN
         NIBUF = 1
         NBITS = 16
      END IF
C
C     Dimensions for CSQ
C     ==================
C
      CALL SPHDIM
C
C     Allocate
C     ========
C
      KACTAB = 1
      KACTCD = KACTAB + (NODCL1 - 1)/IRAT + 1
      KSCRAB = KACTCD + (NODCL2 - 1)/IRAT + 1
      KSCRCD = KSCRAB + NODBC1
      KINDHR = KSCRCD + NODBC2
      KINDSQ = KINDHR + (JTOP3 - 1)/IRAT + 1
      KIODHR = KINDSQ + (NRTOP - 1)/IRAT + 1
      KLMNVL = KIODHR + (8*NRTOP - 1)/IRAT + 1
      KLMNSM = KLMNVL + (3*NHKMAX*KCKMAX - 1)/IRAT + 1
      KPNTUV = KLMNSM + (16*NHKMAX*KCKMAX - 1)/IRAT + 1
      KRBUF  = KPNTUV + (4*KC2MAX*(NRDER+1) - 1)/IRAT + 1
      KIBUF  = KRBUF  + LBUF*NBUFS*NCORS
      KCSQ   = KIBUF  + (LBUF*NBUFS*NIBUF*NCORS - 1)/IRAT + 1
      KCOUNT = KCSQ   + NCSQ1*NCSQ2
      KLAST  = KCOUNT + (NDISTR*NCORS - 1)/IRAT + 1
C
      IF (KLAST .GT. LWORK) CALL STOPIT('ODCDRV','HERPRP',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
C
      IF (SLAVE) THEN
         MXCLS2 = NODCL1*(NODCL2 + 1)/2
C
         NTASK  = INT(NDEGDI*MXCLS2/(100*NODTOT))
         IF (NTASK .EQ. 0) NTASK = 1
C
         KCLCPU = KLAST
         KWHICH = KCLCPU +  MXCLS2
         KIJS   = KWHICH + (MXCLS2 + 1)/IRAT
         KLAST  = KIJS   + (NTASK  + 1)/IRAT
C
         IF (KLAST .GT. LWORK)
     &       CALL STOPIT('ODCDRV','ERI_NODDRV',KLAST,LWORK)
         LWRK = LWORK - KLAST + 1
C
         CALL ERI_NODDRV(IODCL1,IODCL2,
     &                   IODBC1,IODBC2,RODBC1,RODBC2,
     &                   IODPP1,IODPP2,RODPP1,RODPP2,
     &                   FMAT, DMAT,NDMT,IFCTYP,D2MAT,ID2MAT,
     &                   WORK(KACTAB),WORK(KACTCD),
     &                   WORK(KSCRAB),WORK(KSCRCD),
     &                   WORK(KINDHR),WORK(KINDSQ),WORK(KIODHR),
     &                   WORK(KLMNVL),WORK(KLMNSM),WORK(KPNTUV),
     &                   WORK(KRBUF), WORK(KIBUF),
     &                   WORK(KCLCPU),WORK(KWHICH),WORK(KIJS),MXCLS2,
     &                   CCFBT,INDXBT,WORK(KCOUNT),WORK(KCSQ),
     &                   WORK(KLAST),LWRK,IPRINT)
      ELSE
         CALL ODCDR1(IODCL1,IODCL2,
     &               IODBC1,IODBC2,RODBC1,RODBC2,
     &               IODPP1,IODPP2,RODPP1,RODPP2,
     &               FMAT, DMAT,NDMT,IFCTYP,D2MAT,ID2MAT,
     &               WORK(KACTAB),WORK(KACTCD),
     &               WORK(KSCRAB),WORK(KSCRCD),
     &               WORK(KINDHR),WORK(KINDSQ),WORK(KIODHR),
     &               WORK(KLMNVL),WORK(KLMNSM),WORK(KPNTUV),
     &               WORK(KRBUF), WORK(KIBUF),CCFBT,INDXBT,WORK(KCOUNT),
     &               WORK(KCSQ),  WORK(KLAST),LWRK,IPRINT)
      END IF
C
      RETURN
      END
C  /* Deck odcdr1 */
      SUBROUTINE ODCDR1(IODCL1,IODCL2,
     &                  IODBC1,IODBC2,RODBC1,RODBC2,
     &                  IODPP1,IODPP2,RODPP1,RODPP2,
     &                  FMAT,DMAT,NDMT,IFCTYP,D2MAT,ID2MAT,
     &                  NACTAB,NACTCD,SCRNAB,SCRNCD,
     &                  INDHER,INDHSQ,IODDHR,
     &                  LMNPWR,LMNSYM,IPNTUV,
     &                  BUF,IBUF,CCFBT,INDXBT,NCOUNT,
     &                  CSQ,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <eribuf.h>
#include <dummy.h>
#include <cbieri.h>
#include <erithr.h>
#include <odclss.h>
#include <ccom.h>
#include <aobtch.h>
#include <doxyz.h>
#include <ericom.h>
#include <erimem.h>
#include <erista.h>
#include <hertop.h>
#include <symmet.h>
#include <r12int.h>
#include <denfit.h> 
#include <iratdef.h> 
      DIMENSION IODCL1(NODCL1,NITCL), IODCL2(NODCL2,NITCL),
     &          IODBC1(NODBC1,NITBC), IODBC2(NODBC2,NITBC),
     &          RODBC1(NODBC1,NRTBC), RODBC2(NODBC2,NRTBC),
     &          IODPP1(NODPP1,NITPP), IODPP2(NODPP2,NITPP),
     &          RODPP1(NODPP1,NRTPP), RODPP2(NODPP2,NRTPP),
     &          NACTAB(NODCL1), NACTCD(NODCL2),
     &          SCRNAB(NODBC1), SCRNCD(NODBC2),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          INDHSQ(NRTOP), IODDHR(NRTOP,0:7),
     &          LMNPWR(KCKMAX,NHKMAX,3), LMNSYM(KCKMAX,NHKMAX,2,0:7),
     &          IPNTUV(KC2MAX,0:NRDER,2,2),
     &          FMAT(*), DMAT(*), IFCTYP(NDMT), D2MAT(*), ID2MAT(*),
     &          BUF(*), IBUF(*),
     &          CCFBT(*), INDXBT(*), NCOUNT(*),
     &          CSQ(NCSQ1,NCSQ2), WORK(LWORK)
#include <ibtfun.h>
C
C     Information in arrays:
C
C     IODCL1 and IODCL2 give information specific to each class
C     =========================================================
C
C     NITCL = 22 (WK/UniKA/04-11-2002).
C
C     There are NODCL1 and NODCL2 different classes.
C
C     IODCL1(*,1) : NODTAB - total number of OD batches in this class
C     IODCL1(*,2) : KODSAB - address for first OD batch of this class
C                            (used for addressing IODBC1 and RODBC1)
C     IODCL1(*,3) : NHKTA  - quantum number for orbital a
C     IODCL1(*,4) : NHKTB  - quantum nubmer for orbital b
C     IODCL1(*,5) : KHKTA  - number of spherical components for orb. a
C     IODCL1(*,6) : KHKTB  - number of spherical components for orb. b
C     IODCL1(*,7) : KCKTA  - number of Cartesian components for orb. a
C     IODCL1(*,8) : KCKTB  - number of Cartesian components for orb. b
C     IODCL1(*,9) : ISTBLA - stabilizer for orbital a
C     IODCL1(*,10): ISTBLB - stabilizer for orbital b
C
C     IODCL1(*,11): NPRFAB - number of primitives in each OD batch
C     IODCL1(*,12): NCTFAB - number of contracted functions in each OD batch
C     IODCL1(*,13): NCNTAB - one-center distributions: 1
C                          - two-center distributions: 2
C     IODCL1(*,14): NGTOAB - identical GTOs: 1
C                          - different GTOs: 2
C     IODCL1(*,15): NPRFA  - number of primitives for orb. a 
C     IODCL1(*,16): NPRFB  - number of primitives for orb. b 
C     IODCL1(*,17): NCTFA  - number of contracted for orb. a 
C     IODCL1(*,18): NCTFB  - number of contracted for orb. b 
C     IODCL1(*,19): ICMATA - contraction coefficient matrix for orb. a 
C     IODCL1(*,20): ICMATB - contraction coefficient matrix for orb. b 
C
C     IODCL1(*,21): MBSIDA - multiple basis set identifier for orb. a (WK/UniKA/04-11-2002).
C     IODCL1(*,22): MBSIDB - multiple basis set identifier for orb. b (WK/UniKA/04-11-2002).
C
C
C     IODBC1 and IODBC2 give "integer" information for each OD batch
C     ==============================================================
C
C     NITBC = 4
C
C     IODBC1(*,1): OD class of this OD batch
C     IODBC1(*,2): AO batch for orbital a in this OD batch
C     IODBC1(*,3): AO batch for orbital b in this OD batch
C     IODBC1(*,4): start address for primitive OD distributions ab
C                  (used for accessing IODPP1 and RODPP1)
C
C     RODBC1 and RODBC2 give "real" information for each OD batch
C     ===========================================================
C
C     NRTBC = 1
C
C     RODBC1(*,1): Schwarz's inequality factor for this OD batch
C
C     IODPP1 and IODPP2 give "integer" information for each OD
C     ========================================================
C
C     NITPP = 0 (presently unused)
C
C     RODPP1 and RODPP2 give "real" information for each OD
C     =====================================================
C
C     NRTPP = 3
C
C     RODPP1(*,1): exponent for primitive a
C     RODPP1(*,2): exponent for primitive b
C     RODPP1(*,3): preexponential factors for ab
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Output from subroutine ODCDR1',-1)
         WRITE (LUPRI,'(2X,A,2I4)') ' NODCL1, NODCL2 ', NODCL1, NODCL2
         WRITE (LUPRI,'(2X,A,2I4)') ' NODPP1, NODPP2 ', NODPP1, NODPP2
         WRITE (LUPRI,'(2X,A,2L4)') ' ODTRI1, ODTRI2 ', ODTRI1, ODTRI2
         IF (IPRINT .GT. 8) THEN
            CALL HEADER('IODCL1 in ODCDR1',1)
            DO I = 1, NODCL1
               WRITE (LUPRI,'(5X,I4,5X,16I4)') I,(IODCL1(I,J),J=1,NITCL)
            END DO
            CALL HEADER('IODCL2 in ODCDR1',1)
            DO I = 1, NODCL2
               WRITE (LUPRI,'(5X,I4,5X,16I4)') I,(IODCL2(I,J),J=1,NITCL)
            END DO
            CALL HEADER('IODBC1 in ODCDR1',1)
            DO I = 1, NODBC1
               WRITE (LUPRI,'(5X,I4,5X,7I4)') I,(IODBC1(I,J),J=1,NITBC)
            END DO
            CALL HEADER('IODBC2 in ODCDR1',1)
            DO I = 1, NODBC2
               WRITE (LUPRI,'(5X,I4,5X,7I4)') I,(IODBC2(I,J),J=1,NITBC)
            END DO
            CALL HEADER('RODBC1 in ODCDR1',-1)
            CALL OUTPUT(RODBC1,1,NODBC1,1,NRTBC,NODBC1,NRTBC,1,LUPRI)
            CALL HEADER('RODBC2 in ODCDR1',-1)
            CALL OUTPUT(RODBC2,1,NODBC2,1,NRTBC,NODBC2,NRTBC,1,LUPRI)
            CALL HEADER('RODPP1 in ODCDR1',-1)
            CALL OUTPUT(RODPP1,1,NODPP1,1,NRTPP,NODPP1,NRTPP,1,LUPRI)
            CALL HEADER('RODPP2 in ODCDR1',-1)
            CALL OUTPUT(RODPP2,1,NODPP2,1,NRTPP,NODPP2,NRTPP,1,LUPRI)
         END IF
      END IF
C
C     COMMON /DOXYZ/
C
      DOX = .TRUE.
      DOY = .TRUE.
      DOZ = .TRUE.
C
      FIRST = .TRUE.
      LAST  = .FALSE.
C
      NPPA  = 0
      NCPA  = 0
      NPCA  = 0
      NCCA  = 0
      NODPQ = 0
C
      NTPAS = 0
      NPPXA = 0
      NPCXA = 0
      NCPXA = 0
      NCCXA = 0
      NPPSA = 0
      NPCSA = 0
      NCPSA = 0
      NCCSA = 0
      MPRFPQ = 0
      MPQBCS = 0
C
      NPPX1 = 1000000
      NPCX1 = 1000000
      NCPX1 = 1000000
      NCCX1 = 1000000
      NPPS1 = 1000000
      NPCS1 = 1000000
      NCPS1 = 1000000
      NCCS1 = 1000000
C
      NPPX2 = 0
      NPCX2 = 0
      NCPX2 = 0
      NCCX2 = 0
      NPPS2 = 0
      NPCS2 = 0
      NCPS2 = 0
      NCCS2 = 0
C
      IPASSA = 0
      IPASS1 = 1000
      IPASS2 = 0
      IPASSN = 0
C
C     Spherical-harmonic transformation matrices
C     ==========================================
C
      CALL SPHMAT(CSQ,NCSQ1,NCSQ2,MAXDER,GDER,WORK,LWORK,IPRINT)
C
C     Arrays for Hermite integrals
C     ============================
C
      CALL HERPRP(INDHER,INDHSQ,IODDHR,IPRINT)
      CALL LMNPRP(LMNPWR,LMNSYM,IPRINT)
C
C     Active/inactive batches
C     =======================
C
      IELCTR = 1
      CALL ERIACT(NACTAB,SCRNAB,RODBC1,IODBC1,NODBC1,IODCL1,
     &            NODCL1,ODTRI1,IELCTR,IPRINT)
C
      IELCTR = 2
      CALL ERIACT(NACTCD,SCRNCD,RODBC2,IODBC2,NODBC2,IODCL2,
     &            NODCL2,ODTRI2,IELCTR,IPRINT)
C
      IPRKEP = IPRINT
C
C     Loop over pairs of OD classes
C
      DO IODAB = 1, NODCL1
      IF (NACTAB(IODAB) .GT. 0) THEN
         NODSAB = NACTAB(IODAB)
C
         NODTAB = IODCL1(IODAB, 1)
         KODSAB = IODCL1(IODAB, 2)
         SCRMAB = SCRNAB(KODSAB+IDAMAX(NODTAB,SCRNAB(KODSAB),1)-1)
C
         NHKTA  = IODCL1(IODAB, 3)
         NHKTB  = IODCL1(IODAB, 4)
         KHKTA  = IODCL1(IODAB, 5)
         KHKTB  = IODCL1(IODAB, 6)
         KCKTA  = IODCL1(IODAB, 7)
         KCKTB  = IODCL1(IODAB, 8)
         ISTBLA = IODCL1(IODAB, 9)
         ISTBLB = IODCL1(IODAB,10)
         NPRFAB = IODCL1(IODAB,11)
         NCTFAB = IODCL1(IODAB,12)
         NCNTAB = IODCL1(IODAB,13)
         NGTOAB = IODCL1(IODAB,14)
         NPRFA  = IODCL1(IODAB,15)
         NPRFB  = IODCL1(IODAB,16)
         NCTFA  = IODCL1(IODAB,17)
         NCTFB  = IODCL1(IODAB,18)
         ICMATA = IODCL1(IODAB,19)
         ICMATB = IODCL1(IODAB,20)
C 
         IF (ONEAUX .AND. MBSMAX .EQ. 5) THEN
             IF (.NOT. (R12HYB .AND. COMBSS)) THEN
                IF (IODCL1(IODAB,22) .EQ. 2) THEN
                   MAXCD = 0
                   GOTO 1919
                END IF
             ELSE
                IF (IODCL1(IODAB,21) .EQ. 2) THEN
                   MAXCD = 0
                   GOTO 1919
                END IF
            END IF
         END IF
C
C        Sum of basis-set identifiers for a and b (WK/UniKA/04-11-2002).
         MBSMAB = IODCL1(IODAB,21) + IODCL1(IODAB,22)
         MINCD = 1
         IF (DGABAB) MINCD = IODAB
         MAXCD = NODCL2
         IF (ODTR12.OR.DGABAB) MAXCD = IODAB
C
         DO 1918 IODCD = MINCD, MAXCD
         IF (ONEAUX .AND. MBSMAX .EQ. 5) THEN
            IF (.NOT. (R12HYB .AND. COMBSS)) THEN
               IF (IODCL2(IODCD,21) + IODCL2(IODCD,22) .NE. 2) GOTO 1918
            END IF
         END IF
         IF (NACTCD(IODCD) .GT. 0) THEN
            NODSCD = NACTCD(IODCD)
C
            NODTCD = IODCL2(IODCD, 1)
            KODSCD = IODCL2(IODCD, 2)
            SCRMCD = SCRNCD(KODSCD+IDAMAX(NODTCD,SCRNCD(KODSCD),1)-1)
C           Sum of basis-set identifiers for a, b, c and d (WK/UniKA/04-11-2002).
            MBSMCD = MBSMAB + IODCL2(IODCD,21) + IODCL2(IODCD,22)
C
         IF (SCRMAB*SCRMCD .GT. THRSH .AND. MBSMCD .LE. MBSMAX) THEN
C           Only compute integrals for which the sum of basis-set identifiers
C           is less than its prescribed maximum value MBSMAX (WK/UniKA/04-11-2002).
C
            IF (EXTPRI .AND. IODAB.EQ.IPROD1 .AND. IODCD.EQ.IPROD2) THEN
               IPRINT = 1000
            ELSE
               IPRINT = IPRKEP
            END IF
C
            NHKTC  = IODCL2(IODCD, 3)
            NHKTD  = IODCL2(IODCD, 4)
            KHKTC  = IODCL2(IODCD, 5)
            KHKTD  = IODCL2(IODCD, 6)
            KCKTC  = IODCL2(IODCD, 7)
            KCKTD  = IODCL2(IODCD, 8)
            ISTBLC = IODCL2(IODCD, 9)
            ISTBLD = IODCL2(IODCD,10)
            NPRFCD = IODCL2(IODCD,11)
            NCTFCD = IODCL2(IODCD,12)
            NCNTCD = IODCL2(IODCD,13)
Cfp+tuh
C To avoid artificialy many contributions from the same atom in DensFit
C when running diatom scheme
         IF (DIATOM .AND. IFITDM.EQ.-2 .AND. NCNTCD.EQ.2) GO TO 1917
         IF (DIATOM .AND. IFITDM.EQ. 2 .AND. NCNTCD.EQ.1) GO TO 1917
C
            NGTOCD = IODCL2(IODCD,14)
            NPRFC  = IODCL2(IODCD,15)
            NPRFD  = IODCL2(IODCD,16)
            NCTFC  = IODCL2(IODCD,17)
            NCTFD  = IODCL2(IODCD,18)
            ICMATC = IODCL2(IODCD,19)
            ICMATD = IODCL2(IODCD,20)
C
            NODSPQ = NODSAB*NODSCD
            NPRFPQ = NPRFAB*NPRFCD
            NCTFPQ = NCTFAB*NCTFCD
C
C           Statistics
C
            NODPQ  = NODPQ + 1
            NPPI   = NODSPQ*NPRFPQ
            NCPI   = NODSPQ*NCTFAB*NPRFCD
            NPCI   = NODSPQ*NPRFAB*NCTFCD
            NCCI   = NODSPQ*NCTFAB*NCTFCD
            NPPA   = NPPA  + NPPI
            NCPA   = NCPA  + NCPI
            NPCA   = NPCA  + NPCI
            NCCA   = NCCA  + NCCI
C
            JMAXA  = NHKTA - 1
            JMAXB  = NHKTB - 1
            JMAXC  = NHKTC - 1
            JMAXD  = NHKTD - 1
            JMAXAB = JMAXA + JMAXB + MAXDER
            JMAXCD = JMAXC + JMAXD + MAXDER
            JMAX   = JMAXA + JMAXB + JMAXC + JMAXD + MAXDER
C
C           For [r12,T1+T2] integrals, JMAX must be increased by 1. For
C           r12 integrals, JMAX must be 1 for (ss|ss) integrals, but is
C           not increased otherwise. JMAX is also increased for #5 integrals 
C           (WK/UniKA/04-11-2002).
            IF (U12INT) JMAX = JMAX + 1
            IF (INTGAC .EQ. 5) JMAX = JMAX + 1
            IF (R12INT)        JMAX = MAX(1,JMAX)
C
            NTUVAB = (JMAXAB + 1)*(JMAXAB + 2)*(JMAXAB + 3)/6
            NTUVCD = (JMAXCD + 1)*(JMAXCD + 2)*(JMAXCD + 3)/6
            NTUV   = (JMAX   + 1)*(JMAX   + 2)*(JMAX   + 3)/6
C
            SPHRA  = KHKTA .NE. KCKTA .OR. MAXDER.GT.0
            SPHRB  = KHKTB .NE. KCKTB .OR. MAXDER.GT.0
            SPHRC  = KHKTC .NE. KCKTC .OR. MAXDER.GT.0
            SPHRD  = KHKTD .NE. KCKTD .OR. MAXDER.GT.0
C
            SPHRAB = SPHRA .OR. SPHRB
            SPHRCD = SPHRC .OR. SPHRD
C
            GCONA  = GENCON_ERI .OR. NCTFA .GT. 1 
            GCONB  = GENCON_ERI .OR. NCTFB .GT. 1 
            GCONC  = GENCON_ERI .OR. NCTFC .GT. 1 
            GCOND  = GENCON_ERI .OR. NCTFD .GT. 1 
C
            GCONAB = GCONA .OR. GCONB
            GCONCD = GCONC .OR. GCOND
C
            ISTBLR = IBTOR(ISTBLA,ISTBLB)
            ISTBLS = IBTOR(ISTBLC,ISTBLD)
            ISTBLT = IBTOR(IBTAND(ISTBLA,ISTBLB),IBTAND(ISTBLC,ISTBLD))
C
            MLTPA = MULT(ISTBLA)
            MLTPB = MULT(ISTBLB)
            MLTPC = MULT(ISTBLC)
            MLTPD = MULT(ISTBLD)
C
C
            MLTPR = MULT(ISTBLR)
            MLTPS = MULT(ISTBLS)
            MLTPT = MULT(ISTBLT)
            MLTPX = MLTPR*MLTPS*MLTPT
C
            DIAGAB = NGTOAB .EQ. 1 .AND. ODTRI1
            DIAGCD = NGTOCD .EQ. 1 .AND. ODTRI2
C
            IF (MAXDER .EQ. 0) THEN
              TKMPAB = NCNTAB.EQ.1 .AND. KHKTA.EQ.KHKTB .AND..NOT.DODIST
              TKMPCD = NCNTCD.EQ.1 .AND. KHKTC.EQ.KHKTD .AND..NOT.DODIST
            ELSE IF (MAXDER .EQ. 1) THEN
              TKMPAB = .FALSE. 
              TKMPCD = .FALSE. 
            END IF
            TCMPAB = TKMPAB .AND. .NOT.SPHRAB
            TCMPCD = TKMPCD .AND. .NOT.SPHRCD
C
            KHKTAB = KHKTA*KHKTB
            KHKTCD = KHKTC*KHKTD
            KCKTAB = KCKTA*KCKTB
            KCKTCD = KCKTC*KCKTD
            IF (TKMPAB) KHKTAB = KHKTA*(KHKTA + 1)/2
            IF (TKMPCD) KHKTCD = KHKTC*(KHKTC + 1)/2
            IF (TCMPAB) KCKTAB = KCKTA*(KCKTA + 1)/2
            IF (TCMPCD) KCKTCD = KCKTC*(KCKTC + 1)/2
C
            DOPTH1 = .TRUE.
            DOPTH2 = .FALSE.
C
            MLTPZ = MIN(MLTPX,NSPMAX)
            NREDZ = MLTPX/MLTPZ
            IF (MLTPX .NE. NREDZ*MLTPZ) THEN
               WRITE (LUPRI,'(/,1X,A,3(/1X,A,I5))')
     &            ' Error in ODCDR1: MLTPZ incorrectly determined:',
     &            ' NSPMAX = ', NSPMAX, 
     &            ' MLTPZ  = ', MLTPZ, 
     &            ' NREDZ  = ', NREDZ,
     &            ' MLTPX  = ', MLTPX
               CALL QUIT('Error in ODCDR1')
            END IF
C           N2POW = INT(LOG(REAL(MLTPX))/LOG(REAL(2))+0.10D0)
C
C           Memory requirements
C
            CALL MEMPQB(MEMBCH,NALLOI,IPRINT)
            IF (MEMBCH .GT. LWORK) THEN
               WRITE (LUPRI,'(/,1X,A,2I5,2(/1X,A,I10))')
     &              ' Not enough memory for one PQ batch for'//
     &              ' OD classes: ',IODAB, IODCD,
     &              ' Available memory in ODCPAR:',LWORK,
     &              ' Required memory for ODCPAR:',MEMBCH
               MEMOK = .FALSE.
               NEEDMR = MEMBCH - LWORK
               IF (NEEDMR .GT. MEMADD) THEN
                  MEMADD = NEEDMR
                  MODAB  = IODAB
                  MODCD  = IODCD
               END IF
            ELSE IF (MEMOK) THEN
C              MEMBC1 = (IRAT*LWORK-(IRAT-1)*NALLOI)/MEMBCH
C              MXBCHD = MXBCH/(MLTPZ*NPRFAB*NPRFCD)
C              MIXBCH = MIN(IJ,MEMBC1)
C              MINPRZ = MIN(MXBCHD,MIXBCH)
C              MAXBCH = MAX(MINPRZ,1)
               IJ = 0
               DO I = KODSAB, KODSAB + NODTAB - 1
                  MAXJ = KODSCD + NODTCD - 1
                  IF (ODTR12 .AND. IODAB.EQ.IODCD) MAXJ = I
                  DO J = KODSCD, MAXJ
                     IF (SCRNAB(I)*SCRNCD(J).GT.THRSH) IJ = IJ + 1
                  END DO
               END DO
C              MXBCH1 = MIN(IJ,(IRAT*LWORK-(IRAT-1)*NALLOI)/MEMBCH)
C              MAXBCH = MIN(MXBCH/(MLTPZ*NPRFAB*NPRFCD),MXBCH1)
C              IF (MAXBCH .LT. 1) MAXBCH = 1
               MIXBCH = MIN(IJ,(IRAT*LWORK-(IRAT-1)*NALLOI)/MEMBCH)
               MIN1 =   MIN(MXBCH/(MLTPZ*NPRFAB*NPRFCD),MIXBCH)
               MAXBCH = MAX(1,MIN1,MIXBCH)
            END IF
C
            IF (IPRINT .GT. 10) THEN
               CALL TITLER('New pair of OD classes','*',103)
            END IF
            IF (IPRINT .GT. 5) THEN
               WRITE (LUPRI,'(1X,A,2I5)')
     &            ' New OD classes:        ',IODAB,IODCD
            END IF
            IF (IPRINT .GT. 10) THEN
               WRITE (LUPRI,'(1X,A,4I5)')
     &            ' Angular momenta (NHKT):',NHKTA,NHKTB,NHKTC,NHKTD
               WRITE (LUPRI,'(1X,A,3I5)')
     &            ' Number of active OD batches:  ',NODSAB,NODSCD,NODSPQ
               WRITE (LUPRI,'(1X,A,2I5)')
     &            ' NODTAB, NODTCD:        ',NODTAB,NODTCD
               WRITE (LUPRI,'(1X,A,2I5)')
     &            ' KODSAB, KODSCD:        ',KODSAB,KODSCD
               WRITE (LUPRI,'(1X,A,4I5)')
     &            ' Number of primitives:  ',NPRFA,NPRFB,NPRFC,NPRFD
               WRITE (LUPRI,'(1X,A,3I5)')
     &            ' NPRFAB, NPRFCD:        ',NPRFAB,NPRFCD,NPRFPQ
               WRITE (LUPRI,'(1X,A,4I5)')
     &            ' Number of contracted:  ',NCTFA,NCTFB,NCTFC,NCTFD
               WRITE (LUPRI,'(1X,A,3I5)')
     &            ' NCTFAB, NCTFCD:        ',NCTFAB,NCTFCD,NCTFPQ
               WRITE (LUPRI,'(1X,A,4I5)')
     &            ' ICMATA...:             ',ICMATA,ICMATB,ICMATC,ICMATD
               WRITE (LUPRI,'(1X,A,4L5)')
     &            ' GCON...:               ',GCONA,GCONB,GCONC,GCOND
               WRITE (LUPRI,'(1X,A,2L5)')
     &            ' GCONAB,GCONCD:         ',GCONAB,GCONCD
               WRITE (LUPRI,'(1X,A,4I5)')
     &            ' NPP,NCPI,NPCI,NCCI:    ',NPPI,NCPI,NPCI,NCCI
               WRITE (LUPRI,'(1X,A,2I5)')
     &            ' NCNTAB, NCNTCD:        ',NCNTAB,NCNTCD
               WRITE (LUPRI,'(1X,A,2I5)')
     &            ' NGTOAB, NGTOCD:        ',NGTOAB,NGTOCD
               WRITE (LUPRI,'(1X,A,2L4)')
     &            ' TKMPAB, TKMPCD:        ',TKMPAB,TKMPCD
               WRITE (LUPRI,'(1X,A,2I5)')
     &            ' MEMBCH, MAXBCH:        ',MEMBCH,MAXBCH
            END IF
C
            IF (MEMOK) THEN
C
C              Allocate ODCPAR
C
               NITPQ  = 6
               KPNTCR = 1
               KPNTPP = KPNTCR + (4*MAXBCH - 1)/IRAT + 1
               KLAST  = KPNTPP + (2*MAXBCH - 1)/IRAT + 1
               IF (KLAST.GT.LWORK) CALL STOPIT('ODCDR1',' ',KLAST,LWORK)
               LWRK  = LWORK - KLAST + 1
               CALL ODCPAR(IODBC1,IODBC2,
     &                     IODPP1,IODPP2,RODPP1,RODPP2,
     &                     SCRNAB,SCRNCD,
     &                     INDHER,INDHSQ,IODDHR,
     &                     LMNPWR,LMNSYM,IPNTUV,
     &                     WORK(KPNTCR),WORK(KPNTPP),
     &                     FMAT,DMAT,NDMT,IFCTYP,D2MAT,ID2MAT,BUF,IBUF,
     &                     CCFBT,INDXBT,NCOUNT,CSQ,WORK(KLAST),
     &                     LWRK,IPRINT)
            END IF
 1917       CONTINUE
         END IF
         END IF
 1918    CONTINUE
 1919    CONTINUE
      END IF
      END DO
C
C     Empty last buffer
C     =================
C
      IF (WRTINT) THEN
         LAST = .TRUE.
         CALL ERIOUT(DUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,BUF,IBUF,INDXBT,
     &               NCOUNT,0,WORK,LWORK,IPRINT)
      END IF
C
      IF (IPRINT .GT.  1) THEN
         CALL HEADER('Statistics from ODCDR1',1)
         WRITE (LUPRI,'(1X,A,I10)') ' Number of OD pairs:',NODPQ
         WRITE (LUPRI,'(/,1X,A,I10)') ' Number of passes:  ',NTPAS
         IF (NTPAS .GT. 0) THEN
            WRITE (LUPRI,'(/,1X,A,4I10)')
     &          ' NPPXA, NCPXA, NPCXA, NCCXA ',
     &          NINT(FLOAT(NPPXA)/FLOAT(NTPAS)),
     &          NINT(FLOAT(NCPXA)/FLOAT(NTPAS)),
     &          NINT(FLOAT(NPCXA)/FLOAT(NTPAS)),
     &          NINT(FLOAT(NCCXA)/FLOAT(NTPAS))
            WRITE (LUPRI,'(1X,A,4I10)')
     &          ' min values:                ',
     &          NPPX1, NCPX1, NPCX1, NCCX1
            WRITE (LUPRI,'(1X,A,4I10)') ' max values:                ',
     &          NPPX2, NCPX2, NPCX2, NCCX2
            IF (MAXREP .GT. 0) THEN
               WRITE (LUPRI,'(/,1X,A,4I10)')
     &             ' NPPSA, NCPSA, NPCSA, NCCSA ',
     &             NINT(FLOAT(NPPSA)/FLOAT(NTPAS)),
     &             NINT(FLOAT(NCPSA)/FLOAT(NTPAS)),
     &             NINT(FLOAT(NPCSA)/FLOAT(NTPAS)),
     &             NINT(FLOAT(NCCSA)/FLOAT(NTPAS))
               WRITE (LUPRI,'(1X,A,4I10)')
     &             ' min values:                ',
     &             NPPS1, NCPS1, NPCS1, NCCS1
               WRITE (LUPRI,'(1X,A,4I10)')
     &             ' max values:                ',
     &             NPPS2, NCPS2, NPCS2, NCCS2
               WRITE (LUPRI,'(/,1X,A,3I5)')
     &             ' IPASSA, (min and max):',
     &             NINT(FLOAT(IPASSA)/FLOAT(IPASSN)),IPASS1, IPASS2
            END IF
            WRITE (LUPRI,'(/,1X,A,2I10)')
     &          ' NPQBCS, NPRFPQ ',
     &          NINT(FLOAT(MPQBCS)/FLOAT(NTPAS)),
     &          NINT(FLOAT(MPRFPQ)/FLOAT(NTPAS))
         END IF
      END IF
C
      IPRINT = IPRKEP
C
      RETURN
      END
C  /* Deck odcpar */
      SUBROUTINE ODCPAR(IODBC1,IODBC2,
     &                  IODPP1,IODPP2,RODPP1,RODPP2,
     &                  SCRNAB,SCRNCD,
     &                  INDHER,INDHSQ,IODDHR,LMNPWR,LMNSYM,IPNTUV,
     &                  IPNTCR,IPNTPP,
     &                  FMAT,DMAT,NDMT,IFCTYP,D2MAT,ID2MAT,BUF,IBUF,
     &                  CCFBT,INDXBT,NCOUNT,CSQ,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxaqn.h>
#include <maxorb.h>
      LOGICAL LSTBCH
      DIMENSION IODBC1(NODBC1,NITBC), IODBC2(NODBC2,NITBC),
     &          IODPP1(NODPP1,NITPP), IODPP2(NODPP2,NITPP),
     &          RODPP1(NODPP1,NRTPP), RODPP2(NODPP2,NRTPP),
     &          SCRNAB(NODBC1), SCRNCD(NODBC2),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          INDHSQ(NRTOP), IODDHR(NRTOP,0:7),
     &          LMNPWR(KCKMAX,NHKMAX,3), LMNSYM(KCKMAX,NHKMAX,2,0:7),
     &          IPNTUV(KC2MAX,0:NRDER,2,2),
     &          IPNTCR(MAXBCH,4), IPNTPP(MAXBCH,2),
     &          FMAT(*), DMAT(*), IFCTYP(*), D2MAT(*), ID2MAT(*),
     &          BUF(*), IBUF(*),
     &          CCFBT(*), INDXBT(*), NCOUNT(*), CSQ(*), WORK(LWORK)
#include <aobtch.h>
#include <cbieri.h>
#include <ericom.h>
#include <erimem.h>
#include <erithr.h>
#include <eribuf.h>
#include <odclss.h>
#include <erista.h>
#include <hertop.h>
C
C
C     Collect OD batches in passes
C
      IJ = 0
      II = 0
C
      IPASS  = 0
      NBCHES = 0
C
      MINI = KODSAB
      MAXI = KODSAB + NODTAB - 1
      MINJ = KODSCD
      MAXJ = KODSCD + NODTCD - 1
C
      DIAGPQ = (ODTR12.OR.DGABAB) .AND. IODAB.EQ.IODCD
      DIACLS = DIAGPQ .AND. .NOT.EXPERI
C
      KTYPE = 1
      IF (DGABAB) KTYPE = 2
      NTYPE = 1
      IF (DIACLS) NTYPE = 2
      DO ITYPE = KTYPE, NTYPE
C
         DO I = MINI, MAXI
         SCRAB = SCRNAB(I)
         IF (SCRAB*SCRMCD .GT. THRSH) THEN
C
            IF (DIACLS) THEN
               IF (ITYPE .EQ. 1) THEN
                  MAXJ = I - 1
               ELSE
                  MINJ = I
                  MAXJ = I
               END IF
            END IF
C
            DO J = MINJ, MAXJ
C
               IF (SCRAB*SCRNCD(J).GT.THRSH) THEN
                  IJ = IJ + 1
                  IF (ITYPE .EQ. 2) II = II + 1
                  NBCHES = NBCHES + 1
C
C                 Addresses
C
                  IPNTCR(IJ,1) = IODBC1(I,2)
                  IPNTCR(IJ,2) = IODBC1(I,3)
                  IPNTCR(IJ,3) = IODBC2(J,2)
                  IPNTCR(IJ,4) = IODBC2(J,3)
                  IPNTPP(IJ,1) = IODBC1(I,4)
                  IPNTPP(IJ,2) = IODBC2(J,4)
C
C                 Calculate integrals
C
                  IF (IJ.EQ.MAXBCH) THEN
                     CALL ODCPA2(IJ,II,IODPP1,IODPP2,RODPP1,RODPP2,
     &                           IPNTCR,IPNTPP,INDHER,INDHSQ,IODDHR,
     &                           LMNPWR,LMNSYM,IPNTUV,WORK,LWORK,FMAT,
     &                           DMAT,NDMT,IFCTYP,D2MAT,ID2MAT,BUF,
     &                           IBUF,CSQ,CCFBT,
     &                           INDXBT,NCOUNT,IPASS,IPRINT)
                  END IF
               END IF
            END DO
         END IF
         END DO
      END DO
C
      IF (IJ.GT.0) THEN
         CALL ODCPA2(IJ,II,IODPP1,IODPP2,RODPP1,RODPP2,
     &               IPNTCR,IPNTPP,INDHER,INDHSQ,IODDHR,LMNPWR,LMNSYM,
     &               IPNTUV,WORK,LWORK,FMAT,DMAT,NDMT,IFCTYP,D2MAT,
     &               ID2MAT,BUF,IBUF,
     &               CSQ,CCFBT,INDXBT,NCOUNT,IPASS,IPRINT)
      END IF
C
      IPASSA = IPASSA + IPASS
      IPASS1 = MIN(IPASS1,IPASS)
      IPASS2 = MAX(IPASS2,IPASS)
      IPASSN = IPASSN + 1
C
      IF (IPRINT .GT. 10) THEN
         CALL TITLER('Output from ODCPAR','*',103)
         WRITE (LUPRI,'(A,2I5)') '  IPASS,  NTPAS  ', IPASS,  NTPAS
         WRITE (LUPRI,'(A,2I5)') '  NODSAB, NODSCD ', NODSAB, NODSCD
         WRITE (LUPRI,'(A,2I5)') '  NODTAB, NODTCD ', NODTAB, NODTCD
         WRITE (LUPRI,'(A,2I5)') '  KODSAB, KODSCD ', KODSAB, KODSCD
         WRITE (LUPRI,'(A,2I5)') '  NBCHES, NODSPQ ', NBCHES, NODSPQ
         WRITE (LUPRI,'(A,I10)') ' LWORK in ODCPAR: ',LWORK
      END IF
C
      RETURN
      END
C  /* Deck odcpa2 */
      SUBROUTINE ODCPA2(IJ,II,IODPP1,IODPP2,RODPP1,RODPP2,
     &                  IPNTCR,IPNTPP,INDHER,INDHSQ,IODDHR,
     &                  LMNPWR,LMNSYM,
     &                  IPNTUV,WORK,LWORK,FMAT,DMAT,NDMT,IFCTYP,
     &                  D2MAT,ID2MAT,
     &                  BUF,IBUF,CSQ,CCFBT,INDXBT,NCOUNT,IPASS,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxaqn.h>
#include <maxorb.h>
      LOGICAL LSTBCH
      DIMENSION IODPP1(NODPP1,NITPP), IODPP2(NODPP2,NITPP),
     &          RODPP1(NODPP1,NRTPP), RODPP2(NODPP2,NRTPP),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          INDHSQ(NRTOP), IODDHR(NRTOP,0:7),
     &          LMNPWR(KCKMAX,NHKMAX,3), LMNSYM(KCKMAX,NHKMAX,2,0:7),
     &          IPNTUV(KC2MAX,0:NRDER,2,2),
     &          IPNTCR(MAXBCH,4), IPNTPP(MAXBCH,2),
     &          FMAT(*), DMAT(*), IFCTYP(*), D2MAT(*), ID2MAT(*),
     &          BUF(*), IBUF(*),
     &          CCFBT(*), CSQ(*), INDXBT(*), NCOUNT(*), WORK(LWORK)
#include <aobtch.h>
#include <ericom.h>
#include <erimem.h>
#include <erithr.h>
#include <eribuf.h>
#include <odclss.h>
#include <erista.h>
#include <hertop.h>
#include <iratdef.h>
C
      IPASS = IPASS + 1
      NTPAS = NTPAS + 1
C
      NPQBCS = IJ
      NPPBCS = II
C
      IJ = 0
      II = 0
C
      NPPS = NPQBCS*NPRFAB*NPRFCD
      NPCS = NPQBCS*NPRFAB*NCTFCD
      NCPS = NPQBCS*NCTFAB*NPRFCD
      NCCS = NPQBCS*NCTFAB*NCTFCD
C
      NPQBCX = MLTPZ*NPQBCS
      NPPBCX = MLTPZ*NPPBCS
C
      NPPAB = NPQBCX*NPRFAB
      NPPCD = NPQBCX*NPRFCD
C
      NPPX = MLTPZ*NPPS
      NPCX = MLTPZ*NPCS
      NCPX = MLTPZ*NCPS
      NCCX = MLTPZ*NCCS
      NCCT = MLTPX*NCCS
C
      NPPXA = NPPXA + NPPX
      NPCXA = NPCXA + NPCX
      NCPXA = NCPXA + NCPX
      NCCXA = NCCXA + NCCX
C
      MPRFPQ = MPRFPQ + NPRFPQ 
      MPQBCS = MPQBCS + NPQBCS 
C
      NPPX1 = MIN(NPPX1,NPPX)
      NPCX1 = MIN(NPCX1,NPCX)
      NCPX1 = MIN(NCPX1,NCPX)
      NCCX1 = MIN(NCCX1,NCCX)
C
      NPPX2 = MAX(NPPX2,NPPX)
      NPCX2 = MAX(NPCX2,NPCX)
      NCPX2 = MAX(NCPX2,NCPX)
      NCCX2 = MAX(NCCX2,NCCX)
C
      NPPSA = NPPSA + NPPS
      NPCSA = NPCSA + NPCS
      NCPSA = NCPSA + NCPS
      NCCSA = NCCSA + NCCS
C
      NPPS1 = MIN(NPPS1,NPPS)
      NPCS1 = MIN(NPCS1,NPCS)
      NCPS1 = MIN(NCPS1,NCPS)
      NCCS1 = MIN(NCCS1,NCCS)
C
      NPPS2 = MAX(NPPS2,NPPS)
      NPCS2 = MAX(NPCS2,NPCS)
      NCPS2 = MAX(NCPS2,NCPS)
      NCCS2 = MAX(NCCS2,NCCS)
C
      IF (BDER) THEN
         LCENTR = (4*NPQBCX - 1)/IRAT + 1
      ELSE
         LCENTR = 0
      END IF
C
C     Allocate ODCPAS
C
      KCORTR = 1 
      KFCINT = KCORTR + 12*NPQBCS
      KEXPAB = KFCINT +    NPPX
      KEXPCD = KEXPAB +  3*NPPAB
      KALPHA = KEXPCD +  3*NPPCD
      KFACPQ = KALPHA +    NPPX
      KHEXPP = KFACPQ +    NPPS
      KHEXPQ = KHEXPP +    NPPX
      KCENTR = KHEXPQ +    NPPX
      KCORAO = KCENTR +    LCENTR 
      KCORAB = KCORAO + 12*NPQBCX
      KCORCD = KCORAB +  6*NPPAB
      KCORPQ = KCORCD +  6*NPPCD
#if defined (VAR_CRY)
      KROD1 = KCORPQ  +  3*NPPX
      KROD2 = KROD1   +  3*NPQBCS*NPRFAB
      KLAST = KROD2   +  3*NPQBCS*NPRFCD
#else
      KROD1 = KCORPQ  + 3*NPPX
      KROD2 = KROD1
      KLAST = KROD2
#endif
      LWRK  = LWORK - KLAST + 1
      IF (KLAST .GT. LWORK) THEN
         CALL STOPIT('ODCPAR','ODCPAS',KLAST,LWORK)
      END IF
      CALL ODCPAS(IODPP1,IODPP2,RODPP1,RODPP2,IPNTCR,IPNTPP,INDHER,
     &            INDHSQ,IODDHR,LMNPWR,LMNSYM,IPNTUV,WORK(KFCINT),
     &            WORK(KCORTR),WORK(KEXPAB),WORK(KEXPCD),WORK(KALPHA),
     &            WORK(KFACPQ),WORK(KHEXPP),WORK(KHEXPQ),
     &            WORK(KCORAO),WORK(KCORAB),WORK(KCORCD),
     &            WORK(KCORPQ),WORK(KROD1),WORK(KROD2),WORK(KCENTR),
     &            FMAT,DMAT,NDMT,IFCTYP,D2MAT,ID2MAT,BUF,IBUF,CSQ,CCFBT,
     &            INDXBT,NCOUNT,WORK(KLAST),LWRK,IPASS,IPRINT)
C
      RETURN
      END
C  /* Deck odcpas */
      SUBROUTINE ODCPAS(IODPP1,IODPP2,RODPP1,RODPP2,
     &                  IPNTCR,IPNTPP,
     &                  INDHER,INDHSQ,IODDHR,LMNPWR,LMNSYM,IPNTUV,
     &                  FACINT,COORTR,
     &                  EXPAB,EXPCD,ALPHA,FACPQ,HEXPP,HEXPQ,
     &                  COORAO,COORAB,COORCD,COORPQ,
     &                  ROD1,ROD2,NCENTR,FMAT,DMAT,NDMT,IFCTYP,D2MAT,
     &                  ID2MAT,BUF,IBUF,CSQ,CCFBT,INDXBT,NCOUNT,
     &                  WORK,LWORK,IPASS,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
      DIMENSION IODPP1(NODPP1,NITPP), IODPP2(NODPP2,NITPP),
     &          RODPP1(NODPP1,NRTPP), RODPP2(NODPP2,NRTPP),
     &          ROD1(NPQBCS,NPRFAB,3),
     &          ROD2(NPQBCS,NPRFCD,3),
     &          NCENTR(NPQBCX,4),
     &          FMAT(*), DMAT(*), IFCTYP(*), D2MAT(*), ID2MAT(*),
     &          BUF(*), IBUF(*),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          INDHSQ(NRTOP), IODDHR(NRTOP,0:7),
     &          LMNPWR(KCKMAX,NHKMAX,3), LMNSYM(KCKMAX,NHKMAX,2,0:7),
     &          IPNTUV(KC2MAX,0:NRDER,2,2),
     &          IPNTCR(MAXBCH,4), IPNTPP(MAXBCH,2),
     &          FACINT(NPPX), 
     &          COORAO(NPQBCX,3,4),
     &          COORTR(NPQBCS,3,4), 
     &          COORAB(NPPAB,3,2), 
     &          COORCD(NPPCD,3,2), 
     &          EXPAB(NPPAB,3), EXPCD(NPPCD,3), 
     &          ALPHA(NPPX), FACPQ(NPPS), HEXPP(NPPX), HEXPQ(NPPX),
     &          COORPQ(NPPX,3), CSQ(*), CCFBT(*), INDXBT(*), NCOUNT(*),
     &          WORK(LWORK)
#include <cbieri.h>
#include <odclss.h>
#include <ericom.h>
#include <eriao.h>
#include <eribuf.h>
#include <hertop.h>
#include <aobtch.h>
C
      IF (IPRINT .GT. 10) THEN
         CALL TITLER('Output from ODCPAS','*',103)
         WRITE (LUPRI,'(22X,A,I3/)') ' Pass number ',IPASS
         WRITE (LUPRI,'(1X,A,4I5  )')' NHKT? : ',NHKTA,NHKTB,NHKTC,NHKTD
         WRITE (LUPRI,'(1X,A,3I5  )')' NPRF??: ',NPRFAB,NPRFCD,NPRFPQ
         WRITE (LUPRI,'(1X,A,3I5  )')' NCTF??: ',NCTFAB,NCTFCD,NCTFPQ
         WRITE (LUPRI,'(1X,A,2I5  )')' NCNT??: ',NCNTAB,NCNTCD
         WRITE (LUPRI,'(1X,A, I5  )')' NPQBCS: ',NPQBCS
         WRITE (LUPRI,'(1X,A,I10  )')' LWORK in ODCPAS: ',LWORK
C
         WRITE (LUPRI,'(/,1X,A,I5,10X)')' NPPX:   ',NPPX
         WRITE (LUPRI,'(  1X,A,I5,10X)')' NPCX:   ',NPCX
         WRITE (LUPRI,'(  1X,A,I5,10X)')' NPCX:   ',NCPX
         WRITE (LUPRI,'(  1X,A,I5,10X)')' NCCX:   ',NCCX
C
         IF (IPRINT .GT. 15) THEN
            WRITE (LUPRI,'(/,2X,A,13I5/,(15X,13I5))')
     &            'IPNTCR(*,1)  ',(IPNTCR(I,1),I=1,NPQBCS)
            WRITE (LUPRI,'(/,2X,A,13I5/,(15X,13I5))')
     &            'IPNTCR(*,2)  ',(IPNTCR(I,2),I=1,NPQBCS)
            WRITE (LUPRI,'(/,2X,A,13I5/,(15X,13I5))')
     &            'IPNTCR(*,3)  ',(IPNTCR(I,3),I=1,NPQBCS)
            WRITE (LUPRI,'(/,2X,A,13I5/,(15X,13I5))')
     &            'IPNTCR(*,4)  ',(IPNTCR(I,4),I=1,NPQBCS)
            WRITE (LUPRI,'(/,2X,A,13I5/,(15X,13I5))')
     &            'IPNTPP(*,1)  ',(IPNTPP(I,1),I=1,NPQBCS)
            WRITE (LUPRI,'(/,2X,A,13I5/,(15X,13I5))')
     &            'IPNTPP(*,2)  ',(IPNTPP(I,2),I=1,NPQBCS)
         END IF
      END IF
C
C     Collect coordinates and exponents
C     =================================
C
      CALL ERICOR(COORTR,IPNTCR,ROD1,IPRINT)
C     HEXPP and HEXPQ have been added (WK/UniKA/12-11-2002).
      CALL ERIEXP(EXPAB,EXPCD,ALPHA,FACPQ,HEXPP,HEXPQ,RODPP1,RODPP2,
     &            IPNTPP,IPRINT)
      CALL ERISSO(IPRINT)
      IF (BDER) CALL ERICNT(NCENTR,IPNTCR,IPRINT)
C
C     Set IPNTUV
C     ==========
C
      CALL PNTUV(IPNTUV,LMNPWR,INDHER,IPRINT)
C
C     Integrals
C     =========
C
      CALL ERINT0(IPNTCR,IPNTPP,INDHER,INDHSQ,IODDHR,LMNPWR,LMNSYM,
     &            IPNTUV,RODPP1,RODPP2,ROD1,ROD2,COORTR,
     &            FACINT,EXPAB,EXPCD,ALPHA,FACPQ,HEXPP,HEXPQ,
     &            COORAO,COORAB,COORCD,COORPQ,NCENTR,
     &            FMAT,DMAT,NDMT,IFCTYP,D2MAT,ID2MAT,CSQ,BUF,IBUF,
     &            CCFBT,INDXBT,NCOUNT,WORK,LWORK,IPRINT)
      RETURN
      END
C  /* Deck erint0 */
      SUBROUTINE ERINT0(IPNTCR,IPNTPP,INDHER,INDHSQ,IODDHR,LMNPWR,
     &                  LMNSYM,IPNTUV,RODPP1,RODPP2,ROD1,ROD2,
     &                  COORTR,FACINT,EXPAB,EXPCD,ALPHA,FACPQ,
     &                  HEXPP,HEXPQ,
     &                  COORAO,COORAB,COORCD,COORPQ,NCENTR,
     &                  FMAT,DMAT,NDMT,IFCTYP,D2MAT,ID2MAT,CSQ,BUF,IBUF,
     &                  CCFBT,INDXBT,NCOUNT,
     &                  WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <dummy.h>
      DIMENSION FMAT(*), DMAT(*), IFCTYP(*), D2MAT(*), ID2MAT(*), 
     &          CSQ(*), BUF(*), IBUF(*),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          INDHSQ(NRTOP), IODDHR(NRTOP,0:7),
     &          LMNPWR(KCKMAX,NHKMAX,3), LMNSYM(KCKMAX,NHKMAX,2,0:7),
     &          IPNTUV(KC2MAX,0:NRDER,2),
     &          IPNTCR(MAXBCH,4), IPNTPP(MAXBCH,2),
     &          RODPP1(NODPP1,NRTPP), RODPP2(NODPP2,NRTPP),
     &          ROD1(NPQBCS,NPRFAB,3),
     &          ROD2(NPQBCS,NPRFCD,3),
     &          NCENTR(NPQBCX,4),
     &          FACINT(NPPX), COORAO(NPQBCX,3,4),
     &          COORTR(NPQBCS,3,4), 
     &          COORAB(NPPAB,3,2), 
     &          COORCD(NPPCD,3,2), 
     &          EXPAB(NPPAB,3), EXPCD(NPPCD,3), 
     &          ALPHA(NPPX), FACPQ(NPPS), HEXPP(NPPX), HEXPQ(NPPX),
     &          COORPQ(NPPX,3), CCFBT(*), INDXBT(*), 
     &          NCOUNT(*), WORK(LWORK)
#include <cbieri.h>
#include <odclss.h>
#include <ericom.h>
#include <eriao.h>
#include <eribuf.h>
#include <hertop.h>
#include <nuclei.h>
#include <aobtch.h>
#include <symmet.h>
#include <r12int.h>
C
      IF (IPRINT .GT. 10) THEN
         CALL TITLER('Output from ERINT0','*',103)
         WRITE (LUPRI,'(1X,A,I10  )')' LWORK in ERINT0: ',LWORK
      END IF
C
C     Allocate Memory 
C     ---------------
C
C     CCINT    ... | CPINT
C     HERINT | ... |
C            | RJ000
C
      LRJ000 = NPPX*(JMAX + 1)
      LHRIN1 = NPPX*NTUV
      LHRINT = LHRIN1
C     Path 2 not implemented for R12 method (WK/UniKA/04-11-2002).
      IF ((R12INT .OR. R12EIN .OR. U12INT) .AND. DOPTH2) 
     &CALL QUIT('ERINT0: DOPTH2 not implemented for r12 integrals')
C     Extra space (WK/UniKA/04-11-2002).
      IF (R12INT .OR. U12INT .OR. DOPTH2 .OR. MAXDER.GT.0) 
     &LHRINT = 2*LHRIN1
C
      LCPIN1 = 0 
      LCPIN2 = 0 
      IF (DOPTH1 .OR. MAXDER.GT.0) LCPIN1 = NCPX*NTUVCD*KHKTAB
      IF (DOPTH2 .OR. MAXDER.GT.0) LCPIN2 = NPCX*NTUVAB*KHKTCD
      LCPINT = MAX(LCPIN1,LCPIN2) 
      IF (OLDCR1) THEN
         LCPMAX = LCPINT 
      ELSE
         LCPMAX = NCPX*NTUVCD*KCREC1
      END IF
C
      NINTYP = 1
      IF (MAXDER .EQ. 1) THEN
         IF (GDER) THEN
            IF (UNDIFF) THEN
               NINTYP = 13
            ELSE
               NINTYP = 12
            END IF
         END IF
         IF (BDER) THEN
            IF (UNDIFF) THEN
               NINTYP = 7
            ELSE
               NINTYP = 6
            END IF
         END IF
      END IF
      NAOINT = NCCT*KHKTAB*KHKTCD
      LAOINT = NINTYP*NAOINT 
      LSOINT = 0
      LCENTR = 0
      IOFFAO = 0
      IF (MLTPX .GT. 1) LSOINT = LAOINT
      IF (NREDZ .GT. 1) IOFFAO = LAOINT
      IF (MAXDER .EQ. 1) IOFFAO = LAOINT
      IF (MAXDER .EQ. 0) IOFFAO = LAOINT
      LPMAT = 0
      IF (MAXDER .GT. 0) LPMAT = NAOINT 
C
      KAOINT = 1 
      KSOINT = KAOINT + LSOINT
      KHRINT = KAOINT + IOFFAO
      KHRIN2 = KHRINT + LHRIN1
      IF (R12EIN) THEN
C        Space for Gaussian-damped r12 integrals (WK/UniKA/04-11-2002).
         KRJ000 = KAOINT + IOFFAO + 2*LHRIN1
         KCPINT = KAOINT + IOFFAO + MAX(2*LHRIN1,LAOINT)
         KHARGE = KRJ000 + LRJ000
         KHARGF = KHARGE + NPPX
         KHALPH = KHARGF + NPPX
         KHBETA = KHALPH + NPPX
         KHPRE1 = KHBETA + NPPX
         KHPRE2 = KHPRE1 + NPPX
         KHPRE3 = KHPRE2 + NPPX
         KRO000 = KHPRE3 + NPPX
         KHRIN3 = KRO000 + LRJ000
         KHRIN4 = KHRIN3 + LHRIN1
         KHRIN5 = KHRIN4 + LHRIN1
         LSTEIN = KHRIN5 + LHRIN1
      ELSE
         KRJ000 = KAOINT + IOFFAO + LHRINT
         KCPINT = KAOINT + IOFFAO + MAX(LHRINT,LAOINT)
         KHARGE = KRJ000
         KHARGF = KHARGE
         KHALPH = KHARGF
         KHBETA = KHALPH
         KHPRE1 = KHBETA
         KHPRE2 = KHPRE1
         KHPRE3 = KHPRE2
         KRO000 = KHPRE3
         KHRIN3 = KRO000
         KHRIN4 = KHRIN3
         KHRIN5 = KHRIN4
         LSTEIN = KHRIN5
      END IF
C
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Allocations HERCAR in ERINT0',1)
         WRITE (LUPRI,'(1X,A, I10)') ' LRJ000 ', LRJ000
         WRITE (LUPRI,'(1X,A,2I10)') ' LHRINT ', LHRINT,LHRIN1
         WRITE (LUPRI,'(1X,A, I10)') ' LCPINT ', LCPINT
         WRITE (LUPRI,'(1X,A, I10)') ' LCPMAX ', LCPMAX
         WRITE (LUPRI,'(1X,A, I10)') ' LAOINT ', LAOINT
         WRITE (LUPRI,'(1X,A, I10)') ' LSOINT ', LSOINT
         WRITE (LUPRI,'(1X,A, I10)') ' KAOINT ', KAOINT
         WRITE (LUPRI,'(1X,A, I10)') ' KSOINT ', KSOINT
         WRITE (LUPRI,'(1X,A, I10)') ' KHRINT ', KHRINT
         WRITE (LUPRI,'(1X,A, I10)') ' KRJ000 ', KRJ000
         WRITE (LUPRI,'(1X,A, I10)') ' KCPINT ', KCPINT
C        Output for Gaussian-damped r12 integrals (WK/UniKA/04-11-2002).
         IF (R12EIN) THEN
            WRITE (LUPRI,'(1X,A, I10)') ' KHARGE ', KHARGE
            WRITE (LUPRI,'(1X,A, I10)') ' KHARGF ', KHARGF
            WRITE (LUPRI,'(1X,A, I10)') ' KHALPH ', KHALPH
            WRITE (LUPRI,'(1X,A, I10)') ' KHBETA ', KHBETA
            WRITE (LUPRI,'(1X,A, I10)') ' KHPRE1 ', KHPRE1
            WRITE (LUPRI,'(1X,A, I10)') ' KHPRE2 ', KHPRE2
            WRITE (LUPRI,'(1X,A, I10)') ' KHPRE3 ', KHPRE3
            WRITE (LUPRI,'(1X,A, I10)') ' KRO000 ', KRO000
            WRITE (LUPRI,'(1X,A, I10)') ' KHRIN3 ', KHRIN3
            WRITE (LUPRI,'(1X,A, I10)') ' KHRIN4 ', KHRIN4
            WRITE (LUPRI,'(1X,A, I10)') ' KHRIN5 ', KHRIN5
         END IF
      END IF
C
      IF (R12EIN) THEN
C        Space for Gaussian-damped r12 integrals (WK/UniKA/04-11-2002).
         LSTHER = KRJ000 + 7*NPPX + 2*LRJ000 + 3*LHRIN1
         LSTGAM = KRJ000
      ELSE
         LSTHER = KRJ000 + LRJ000
         LSTGAM = KRJ000 + LRJ000
      END IF
      LSTCR1 = KCPINT + LCPMAX
      LSTCR2 = KCPINT + LCPINT
      LSTOUT = KSOINT + LAOINT
C
      LWKGAM = LWORK - LSTGAM + 1
      LWKHER = LWORK - LSTHER + 1
      LWKCR1 = LWORK - LSTCR1 + 1
      LWKCR2 = LWORK - LSTCR2 + 1
      LWKOUT = LWORK - LSTOUT + 1
C
      IF (LSTGAM.GT.LWORK) CALL STOPIT('ERINT0','LSTGAM',LSTGAM,LWORK)
      IF (LSTHER.GT.LWORK) CALL STOPIT('ERINT0','LSTHER',LSTHER,LWORK)
      IF (LSTCR1.GT.LWORK) CALL STOPIT('ERINT0','LSTCR1',LSTCR1,LWORK)
      IF (LSTCR2.GT.LWORK) CALL STOPIT('ERINT0','LSTCR2',LSTCR2,LWORK)
      IF (LSTOUT.GT.LWORK) CALL STOPIT('ERINT0','LSTOUT',LSTOUT,LWORK)
      IF (LSTEIN.GT.LWORK) CALL STOPIT('ERINT0','LSTEIN',LSTEIN,LWORK)
C
C     **********************************
C     *** Undifferentiated integrals ***
C     **********************************
C
      IF (MAXDER .EQ. 0) THEN
C
C        Calculate integrals
C        ===================
C
C        Extra work space for Gaussian-damped R12 integrals 
C        has been added (WK/UniKA/04-11-2002).
C
         NDIMD  = 1
         CALL ERINTS(WORK(KAOINT),WORK(KHRINT),WORK(KHRIN2),
     &               WORK(KRJ000),WORK(KCPINT),
     &               IPNTPP,INDHER,INDHSQ,IODDHR,LMNPWR,LMNSYM,IPNTUV,
     &               RODPP1,RODPP2,ROD1,ROD2,
     &               COORTR,FACINT,CCFBT,
     &               EXPAB,EXPCD,ALPHA,FACPQ,
     &               COORAO,COORAB,COORCD,COORPQ,CSQ,NCENTR,
     &               WORK(LSTGAM),WORK(LSTHER),WORK(LSTCR1),
     &               WORK(LSTCR2),
     &               WORK(KHARGE),WORK(KHARGF),WORK(KHALPH),
     &               WORK(KHBETA),WORK(KHPRE1),WORK(KHPRE2),
     &               WORK(KHPRE3),WORK(KRO000),
     &               WORK(KHRIN3),WORK(KHRIN4),WORK(KHRIN5),
     &               HEXPP,HEXPQ,LWKGAM,LWKHER,LWKCR1,LWKCR2,
     &               NDIMD,IPRINT)
C
C        Symmetrize
C        ==========
C
         IF (MLTPX.GT.1) THEN
            CALL ERISYM(WORK(KAOINT),WORK(KSOINT),0,0,0,
     &                  IODDHR(1,ICCXYZ),IPNTUV,IPRINT)
         END IF
C
C        Write out integrals
C        ===================
C
         IF (WRTINT) THEN
            CALL ERIOUT(WORK(KSOINT),IDUMMY,IPNTCR,IODDHR(1,ICCXYZ),
     &                  IPNTUV,BUF,IBUF,INDXBT,NCOUNT,0,WORK(LSTOUT),
     &                  LWKOUT,IPRINT)
         END IF
C
C        Construct Fock matrix
C        =====================
C
         IF (FCKINT) THEN
            CALL ERIFOK(WORK(KSOINT),IPNTCR,IODDHR(1,ICCXYZ),IPNTUV,
     &                  FMAT,DMAT,NDMT,IFCTYP,CCFBT,INDXBT,
     &                  WORK(LSTOUT),LWKOUT,IPRINT)
         END IF
C
C        Add to Cauchy-Schwarz matrix
C
         IF (DGABAB) THEN
            CALL ERICSM(WORK(KSOINT),FMAT,IPNTCR,INDXBT,WORK(LSTOUT),
     &                  LWKOUT,IPRINT)
         END IF
C
C     ********************************
C     *** Differentiated integrals ***
C     ********************************
C
      ELSE
         IF (GDER) THEN
            NDIMD  = 12
            IF (UNDIFF) NDIMD = 13
         ELSE IF (BDER) THEN
            NDIMD = 6
            IF (UNDIFF) NDIMD = 7
         END IF
C
C        Extra work space for Gaussian-damped R12 integrals 
C        has been added (WK/UniKA/04-11-2002).
C
         CALL ERINTS(WORK(KAOINT),WORK(KHRINT),WORK(KHRIN2),
     &               WORK(KRJ000),WORK(KCPINT),
     &               IPNTPP,INDHER,INDHSQ,IODDHR,LMNPWR,LMNSYM,IPNTUV,
     &               RODPP1,RODPP2,ROD1,ROD2,
     &               COORTR,FACINT,CCFBT,
     &               EXPAB,EXPCD,ALPHA,FACPQ,
     &               COORAO,COORAB,COORCD,COORPQ,CSQ,NCENTR,
     &               WORK(LSTGAM),WORK(LSTHER),WORK(LSTCR1),
     &               WORK(LSTCR2),
     &               WORK(KHARGE),WORK(KHARGF),WORK(KHALPH),
     &               WORK(KHBETA),WORK(KHPRE1),WORK(KHPRE2),
     &               WORK(KHPRE3),WORK(KRO000),
     &               WORK(KHRIN3),WORK(KHRIN4),WORK(KHRIN5),
     &               HEXPP,HEXPQ,LWKGAM,LWKHER,LWKCR1,LWKCR2,
     &               NDIMD,IPRINT)
C
C        Symmetrization
C        ==============
C
         IF (WRTINT .AND. GDER) THEN
            IF (MLTPX.GT.1) THEN
               DO IREPE = 0, MAXREP
                  KADD = 0
                  DO I = 1, NINTYP
                     ICRB = MOD(I-1,3)+1
                     IATOM = (I + 2)/3 
                     IF (UNDIFF .AND. I.EQ.NINTYP) ICRB = 0
                     CALL ERISYM(WORK(KAOINT+KADD),WORK(KSOINT+KADD),
     &                           IATOM,ICRB,IREPE,IODDHR(1,ICCXYZ),
     &                           IPNTUV,IPRINT)
                     KADD = KADD + NAOINT
                  END DO
                  KINT = KSOINT
                  CALL ERIOUT(WORK(KINT),I,IPNTCR,IODDHR(1,ICCXYZ),
     &                        IPNTUV,BUF,IBUF,INDXBT,NCOUNT,IREPE,
     &                        WORK(LSTOUT),LWKOUT,IPRINT)
              END DO
            ELSE
               KINT = KAOINT
               DO IREPE = 0, MAXREP
                  CALL ERIOUT(WORK(KINT),I,IPNTCR,IODDHR(1,ICCXYZ),
     &                        IPNTUV,BUF,IBUF,INDXBT,NCOUNT,IREPE,
     &                        WORK(LSTOUT),LWKOUT,IPRINT)
               END DO
            END IF
         ELSE
C
C           Symmetrization
C           ==============
C
            IF (MLTPX.GT.1) THEN
               KADD = 0
               DO I = 1, NINTYP
                  ICRB = MOD(I-1,3)+1
                  IATOM = (I + 2)/3 
                  IF (UNDIFF .AND. I.EQ.NINTYP) ICRB = 0
                  CALL ERISYM(WORK(KAOINT+KADD),WORK(KSOINT+KADD),
     &                        IATOM,ICRB,0,IODDHR(1,ICCXYZ),IPNTUV,
     &                        IPRINT)
                  KADD = KADD + NAOINT
               END DO
               KINT = KSOINT
            ELSE
               KINT = KAOINT
            END IF
C
C           Write integrals
C           ===============
C
            IF (WRTINT) THEN
               IF (GDER .OR. BDER) THEN
                  CALL ERIOUT(WORK(KINT),I,IPNTCR,IODDHR(1,ICCXYZ),
     &                        IPNTUV,BUF,IBUF,INDXBT,NCOUNT,0,
     &                        WORK(LSTOUT),LWKOUT,IPRINT)
               ELSE
                  CALL QUIT('Fix code for this case before ERIOUT')
               END IF
            END IF
C
C           Expectation value
C           =================
C
            IF (EXPERI) THEN
               KPMAT = KINT  + LAOINT
               KLAST = KPMAT + NAOINT
               LWRK  = LWORK - KLAST + 1
               IF (KLAST.GT.LWORK) CALL STOPIT('ERINT0','ERIAVE',
     &                                          LWRK,LWORK)
               CALL ERIAVE(WORK(KINT),WORK(KPMAT),DMAT,D2MAT,ID2MAT,
     &                     IPNTCR,IODDHR(1,ICCXYZ),IPNTUV,INDXBT,
     &                     WORK(KLAST),LWRK,IPRINT)
            END IF
         END IF
      END IF
C
      RETURN
      END
C  /* Deck erints */
      SUBROUTINE ERINTS(AOINT,HRINT1,HRINT2,RJ000,CPINT,
     &                  IPNTPP,INDHER,INDHSQ,IODDHR,LMNPWR,LMNSYM,
     &                  IPNTUV,RODPP1,RODPP2,ROD1,ROD2,
     &                  COORTR,FACINT,CCFBT,EXPAB,EXPCD,ALPHA,FACPQ,
     &                  COORAO,COORAB,COORCD,COORPQ,CSQ,NCENTR,
     &                  WRKGAM,WRKHER,WRKCR1,WRKCR2,
     &                  HARGE,HARGF,HALPH,HBETA,HPRE1,HPRE2,
     &                  HPRE3,RO000,HRIN1,HRIN2,HRIN3,
     &                  HEXPP,HEXPQ,LWKGAM,LWKHER,LWKCR1,LWKCR2,
     &                  NDIMD,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
C
C     Arrays for Gaussian-damped R12 integrals (WK/UniKA/04-11-2002).
      DIMENSION HARGE(*), HARGF(*), HALPH(*), HBETA(*), HPRE1(*), 
     &          HPRE2(*), HPRE3(*), RO000(*), HRIN1(*),
     &          HRIN2(*), HRIN3(*), HEXPP(*), HEXPQ(*)
      DIMENSION AOINT(NCCS,MLTPX,KHKTAB*KHKTCD,NDIMD),
     &          HRINT1(*),HRINT2(*),RJ000(*),CPINT(*),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          INDHSQ(NRTOP), IODDHR(NRTOP,0:7),
     &          LMNPWR(KCKMAX,NHKMAX,3), LMNSYM(KCKMAX,NHKMAX,2,0:7),
     &          IPNTUV(KC2MAX,0:NRDER,2), IPNTPP(MAXBCH,2),
     &          FACINT(NPPX), COORAO(NPQBCX,3,4), CCFBT(*),
     &          RODPP1(NODPP1,NRTPP), RODPP2(NODPP2,NRTPP),
     &          ROD1(NPQBCS,NPRFAB,3),
     &          ROD2(NPQBCS,NPRFCD,3),
     &          COORTR(NPQBCS,3,4),
     &          COORAB(NPPAB,3,2), 
     &          COORCD(NPPCD,3,2), 
     &          EXPAB(NPPAB,3), EXPCD(NPPCD,3), 
     &          ALPHA(NPPX), FACPQ(NPPS),
     &          COORPQ(NPPX,3),
     &          CSQ(*), NCENTR(NPQBCS,MLTPZ,4)
c     C. Villani Uni-KA, Feb 2005
      LOGICAL   MIOTST
#include <cbieri.h>
#include <odclss.h>
#include <ericom.h>
#include <eriao.h>
#include <erisop.h>
#include <eribuf.h>
#include <hertop.h>
#include <r12int.h>
#include <aobtch.h>
#include <clsfmm.h>
C
      IF (IPRINT .GT. 10) THEN
         CALL TITLER('Output from ERINTS','*',103)
      END IF
C
C     IF (NREDZ .GT. 1) CALL DZERO(AOINT,NCCT*KHKTAB*KHKTCD*NDIMD)
      CALL DZERO(AOINT,NCCT*KHKTAB*KHKTCD*NDIMD)
C     Try to avoid zeroing - check sph2 rouintes.
c
c     Variable for integral printing (debug)
      miotst=.false.
c     Initializing variable for [T1,r12] derivative integrals
c     C. Villani Uni-KA, Feb 2005
      miau12=.false.
C
      IOFRST = 0
      DO IOPS = 1, NREDZ 
C
         IF (IPRINT .GT. 10) THEN
            WRITE (LUPRI,'(//,1X,A,I5,A,I5,A,/)')
     &         ' Symmetry pass ',IOPS, ' out of ',NREDZ,' in ERINTS.'
            WRITE (LUPRI,'(1X,A,I5)') ' IOFRST ', IOFRST 
         END IF
C
C        Coordinates and exponents
C        =========================
C
         CALL ERICRS(COORAO,COORTR,IPRINT)
C
C        Preexponential factors
C        =======================
C
         CALL ERIPFS(FACINT,FACPQ,COORAB,COORCD,
     &               RODPP1(1,3),RODPP2(1,3),
     &               IPNTPP,ROD1,ROD2,IPRINT)
C
C        Set up coordinate vectors
C        =========================
C
         CALL ERIVCS(EXPAB,EXPCD,COORAO,COORAB,COORCD,COORPQ,IPRINT)
C
C        Test for classical and nonclassical integrals
C        =============================================
C
         IF (NCLERI) CALL ERICLS(FACINT,COORPQ,EXPAB,EXPCD,IPRINT)
C
C        Incomplete gamma function
C        =========================
C
         IF (.NOT. R12EIN)
     &   CALL ERIGAM(RJ000,FACINT,ALPHA,COORPQ,WRKGAM,LWKGAM,IPRINT)
C        ERIGAM is not called for Gaussian-damped r12 integrals (WK/UniKA/04-11-2002).
C
C        Hermite integrals
C        =================
C
         CALL ERIHER(HRINT1,HRINT2,RJ000,COORPQ,INDHER,
     &               IODDHR(1,IPQXYZ),FACINT,ALPHA,
     &               HARGE,HARGF,HALPH,HBETA,HPRE1,HPRE2,
     &               HPRE3,RO000,HRIN1,HRIN2,HRIN3,
     &               HEXPP,HEXPQ,WRKHER,LWKHER,IPRINT)
C        ERIHER replaces call to ERITUV (WK/UniKA/04-11-2002).
C 
C        A: Undifferentiated Cartesian integrals
C        =======================================
C
         IF (MAXDER .EQ. 0) THEN
            NDER3 = 0 
            NDER4 = 0 
            IF (DOPTH1) THEN
               IPATH = 1
               CALL SETCAO(IPRINT)
               CALL CR1DRV(HRINT1,HRINT2,CPINT,INDHER,IODDHR(1,IPQXYZ),
     &                     IODDHR(1,IHCBIT),INDHSQ,LMNPWR,IPNTUV,
     &                     COORAB,EXPAB,CSQ,CCFBT,WRKCR1,LWKCR1,IPRINT)
               CALL CR2DRV(AOINT(1,IOFRST+1,1,1),CPINT,INDHER,
     &                     IODDHR(1,IHCBIT),IODDHR(1,ICDXYZ),
     &                     IODDHR(1,ICCXYZ),LMNPWR,IPNTUV,
     &                     COORCD,EXPCD,CSQ,CCFBT,NCENTR(1,IOFRST+1,1),
     &                     NDER3,NDER4,NDIMD,
     &                     WRKCR2,LWKCR2,IPRINT)
            ELSE
               IPATH = 2
               CALL SETCAO(IPRINT)
               CALL ERISWP(HRINT1,HRINT2,INDHER,IODDHR(1,IPQXYZ),IPRINT)
               CALL CR1DRV(HRINT2,HRINT1,CPINT,INDHER,IODDHR(1,IPQXYZ),
     &                     IODDHR(1,IHCBIT),INDHSQ,LMNPWR,IPNTUV,
     &                     COORCD,EXPCD,CSQ,CCFBT,WRKCR1,LWKCR1,IPRINT)
               CALL CR2DRV(AOINT(1,IOFRST+1,1,1),CPINT,INDHER,
     &                     IODDHR(1,IHCBIT),IODDHR(1,IABXYZ),
     &                     IODDHR(1,ICCXYZ),LMNPWR,IPNTUV,
     &                     COORAB,EXPAB,CSQ,CCFBT,NCENTR(1,IOFRST+1,1),
     &                     NDER3,NDER4,NDIMD,WRKCR2,LWKCR2,IPRINT)
            END IF
C
C        B: Magnetic field derivatives 
C        =============================
C
         ELSE IF (BDER) THEN
C
C           Derivatives for second electron 
C           -------------------------------
C
            IPATH  = 1
            CALL SETCAO(IPRINT)
            CALL CR1DRV(HRINT1,HRINT2,CPINT,INDHER,IODDHR(1,IPQXYZ),
     &                  IODDHR(1,IHCBIT),INDHSQ,LMNPWR,IPNTUV,
     &                  COORAB,EXPAB,CSQ,CCFBT,WRKCR1,LWKCR1,IPRINT)
C
            IF (IPRINT.GT.5) CALL HEADER('Undifferentiated integrals',1)
C
            NDER3 = 0
            NDER4 = 0
            NDIMD1 = 1
            CALL CR2DRV(AOINT(1,IOFRST+1,1,7),CPINT,INDHER,
     &                  IODDHR(1,IHCBIT),IODDHR(1,ICDXYZ),
     &                  IODDHR(1,ICCXYZ),LMNPWR,IPNTUV,
     &                  COORCD,EXPCD,CSQ,CCFBT,NCENTR(1,IOFRST+1,3),
     &                  NDER3,NDER4,NDIMD1,WRKCR2,LWKCR2,IPRINT)
            CALL ERIBUN(AOINT(1,IOFRST+1,1,1),AOINT(1,IOFRST+1,1,7),
     &                  COORAO,IODDHR(1,ICCXYZ),IPNTUV,IPRINT)
C
            IF (IPRINT.GT.5) CALL HEADER('C derivatives',1)
C
            NDER3 = 1
            NDER4 = 0
            NDIMD1 = 3
            CALL CR2DRV(AOINT(1,IOFRST+1,1,4),CPINT,INDHER,
     &                  IODDHR(1,IHCBIT),IODDHR(1,ICDXYZ),
     &                  IODDHR(1,ICCXYZ),LMNPWR,IPNTUV,
     &                  COORCD,EXPCD,CSQ,CCFBT,NCENTR(1,IOFRST+1,3),
     &                  NDER3,NDER4,NDIMD1,WRKCR2,LWKCR2,IPRINT)
C
C           Derivatives for first electron 
C           ------------------------------
C
            IF (IPRINT.GT.5) CALL HEADER('A derivatives',1)
C
            IPATH  = 2
            NDIMD2 = 3
            CALL SETCAO(IPRINT)
            CALL ERISWP(HRINT1,HRINT2,INDHER,IODDHR(1,IPQXYZ),IPRINT)
            CALL CR1DRV(HRINT2,HRINT1,CPINT,INDHER,IODDHR(1,IPQXYZ),
     &                  IODDHR(1,IHCBIT),INDHSQ,LMNPWR,IPNTUV,
     &                  COORCD,EXPCD,CSQ,CCFBT,WRKCR1,LWKCR1,IPRINT)
            NDER3 = 1
            NDER4 = 0
            CALL CR2DRV(AOINT(1,IOFRST+1,1,1),CPINT,INDHER,
     &                  IODDHR(1,IHCBIT),IODDHR(1,IABXYZ),
     &                  IODDHR(1,ICCXYZ),LMNPWR,IPNTUV,
     &                  COORAB,EXPAB,CSQ,CCFBT,NCENTR(1,IOFRST+1,1),
     &                  NDER3,NDER4,NDIMD2,WRKCR2,LWKCR2,IPRINT)
C
C        B: Differentiated Cartesian integrals
C        =====================================
C
         ELSE
           IF(IPRINT.GT.5) THEN
c            Header for integral derivatives over different operators
c            C. Villani Uni-KA, Feb 2005
             WRITE(LUPRI,*) "-----------------------------------------"
             IF(U12INT) THEN
c              if(npqbcs.ne.1.or.mltpx.ne.1.or.mltpz.ne.1)
c    .         write(lupri,*) "warning: may not work"
               WRITE(LUPRI,*) "   Derivatives of [T1,r12] integrals "
               XMIADR=0
             ELSE IF(R12INT) THEN
               WRITE(LUPRI,*) "   Derivatives of r12 integrals      "
             ELSE
               WRITE(LUPRI,*) "   Derivatives of 1/r12 integrals    "
             END IF
             WRITE(LUPRI,*) "-----------------------------------------"
           END IF
C
C           CD derivatives
C           --------------
C
            IPATH  = 1
            NDIMD1 = 3
            CALL SETCAO(IPRINT)
            IF(U12INT) THEN
c             symmetry control (ioddhr) is here always switched off,
c             only that for the hermite integrals (ipqxyz) is kept.
c             C. Villani Uni-KA, Feb 2005
              IHCBIT=0
              ICCXYZ=0
            END IF
            CALL CR1DRV(HRINT1,HRINT2,CPINT,INDHER,IODDHR(1,IPQXYZ),
     &                  IODDHR(1,IHCBIT),INDHSQ,LMNPWR,IPNTUV,
     &                  COORAB,EXPAB,CSQ,CCFBT,WRKCR1,LWKCR1,IPRINT)
            IF (UNDIFF) THEN
               IF (IPRINT.GT.5) CALL HEADER('Undiff. integrals',1)
               NDER3 = 0
               NDER4 = 0
               IAOCR = 13 
c              NDIMD1 = 3
               NDIMD1=1
c              Here it is just needed ndimd1=1, oder????????
c              C. Villani Uni-KA, Feb 2005
               CALL CR2DRV(AOINT(1,IOFRST+1,1,IAOCR),CPINT,INDHER,
     &                     IODDHR(1,IHCBIT),IODDHR(1,ICDXYZ),
     &                     IODDHR(1,ICCXYZ),LMNPWR,IPNTUV,
     &                     COORCD,EXPCD,CSQ,CCFBT,NCENTR(1,IOFRST+1,3),
     &                     NDER3,NDER4,NDIMD1,WRKCR2,LWKCR2,IPRINT)
               NDIMD1=3
            END IF
            IF (IPRINT.GT.5) CALL HEADER('Calculating C derivatives',1)
            NDER3 = 1
            NDER4 = 0
            IAOCR = 7
            CALL CR2DRV(AOINT(1,IOFRST+1,1,IAOCR),CPINT,INDHER,
     &                  IODDHR(1,IHCBIT),IODDHR(1,ICDXYZ),
     &                  IODDHR(1,ICCXYZ),LMNPWR,IPNTUV,
     &                  COORCD,EXPCD,CSQ,CCFBT,NCENTR(1,IOFRST+1,3),
     &                  NDER3,NDER4,NDIMD1,WRKCR2,LWKCR2,IPRINT)
            IF (IPRINT.GT.5) CALL HEADER('Calculating D derivatives',1)
            NDER3 = 0
            NDER4 = 1
            IAOCR = 10
            CALL CR2DRV(AOINT(1,IOFRST+1,1,IAOCR),CPINT,INDHER,
     &                  IODDHR(1,IHCBIT),IODDHR(1,ICDXYZ),
     &                  IODDHR(1,ICCXYZ),LMNPWR,IPNTUV,
     &                  COORCD,EXPCD,CSQ,CCFBT,NCENTR(1,IOFRST+1,3),
     &                  NDER3,NDER4,NDIMD1,WRKCR2,LWKCR2,IPRINT)
C
C           AB derivatives
C           --------------
C
            IF(U12INT) THEN
c             half-computed A and B derivatives [ (d/dP)(d/dR)(ab|r12|cd) ]
c             C. Villani (Uni Karlsruhe Feb 2005)
              IPATH  = 1
              NDIMD2 = 1
c             decreasing ntuv34, just for this derivative part
              nsave34=ntuv34
              ntuv34= max(1,(JMAXCD + 1)*(JMAXCD + 2)*(JMAXCD)/6)
c             calculates half-transformed integrals before eriswp call
              DO XMIADR=1,3
c               loop over x,y,z coordinates
c 
c               A half-derivatives
                AENONB=.TRUE.
                if(iprint.gt.5)
     .            write(lupri,*) 'half-computed A derivatives ', xmiadr
                CALL CR1MIA(HRINT1,HRINT2,CPINT,INDHER,IODDHR(1,IPQXYZ),
     &                     IODDHR(1,0),INDHSQ,LMNPWR,IPNTUV,COORAB,
     &                     EXPAB,CSQ,CCFBT,WRKCR1,LWKCR1,IPRINT)
                NDER3 = 0
                NDER4 = 0
                IAOCR = XMIADR
                CALL CR2DRV(AOINT(1,IOFRST+1,1,IAOCR),CPINT,INDHER,
     &                      IODDHR(1,0),IODDHR(1,0),
     &                      IODDHR(1,0),LMNPWR,IPNTUV,
     &                      COORCD,EXPCD,CSQ,CCFBT,NCENTR(1,IOFRST+1,3),
     &                      NDER3,NDER4,NDIMD2,WRKCR2,LWKCR2,IPRINT)
c               B half-derivatives
                AENONB=.FALSE.
                if(iprint.gt.5)
     .            write(lupri,*) 'half-computed B derivatives ', xmiadr
c 
                CALL CR1MIA(HRINT1,HRINT2,CPINT,INDHER,IODDHR(1,IPQXYZ),
     &                 IODDHR(1,0),INDHSQ,LMNPWR,IPNTUV,COORAB,
     &                 EXPAB,CSQ,CCFBT,WRKCR1,LWKCR1,IPRINT)
                NDER3 = 0
                NDER4 = 0
                IAOCR = 3 + XMIADR
                CALL CR2DRV(AOINT(1,IOFRST+1,1,IAOCR),CPINT,INDHER,
     &                      IODDHR(1,0),IODDHR(1,0),
     &                      IODDHR(1,0),LMNPWR,IPNTUV,
     &                      COORCD,EXPCD,CSQ,CCFBT,NCENTR(1,IOFRST+1,3),
     &                      NDER3,NDER4,NDIMD2,WRKCR2,LWKCR2,IPRINT)
c               just changing sign for B half-derivatives
                call miach(aoint(1,iofrst+1,1,iaocr),ncct*khktab*khktcd)
              END DO
c             giving ntuv34 its value back
              ntuv34= nsave34
              XMIADR=0
c             switching on the option for the calculation
c             of (a-b)/(a+b) d/dA (ab|1/r12|cd)
c             and (a-b)/(a+b) d/dB (ab|1/r12|cd) integral derivatives
c             contribution
              MIAU12=.TRUE.
c             U12INT need to be false in the next call to CR1DRV
c             R12INT=.FALSE.
              U12INT=.FALSE.
            END IF
c
            IF(MIAU12.AND.NGTOAB.EQ.1.AND.NPRFAB.EQ.1) THEN
c             case a=b, uncontracted,      (a-b)/(a+b) = 0
c             For derivatives of [T1,r12] integrals,
c             Just skip the computation of the second contribution
c
c             C. Villani (Uni Karlsruhe Feb 2005)
              if(iprint.gt.5)
     .        write(lupri,*) 'skipping A and B second-half derivatives'
            ELSE
              IPATH  = 2
              NDIMD2 = 3
              CALL SETCAO(IPRINT)
              CALL ERISWP(HRINT1,HRINT2,INDHER,IODDHR(1,IPQXYZ),IPRINT)
              CALL CR1DRV(HRINT2,HRINT1,CPINT,INDHER,IODDHR(1,IPQXYZ),
     &                    IODDHR(1,IHCBIT),INDHSQ,LMNPWR,IPNTUV,
     &                    COORCD,EXPCD,CSQ,CCFBT,WRKCR1,LWKCR1,IPRINT)
C
              IF(IPRINT.GT.5) CALL HEADER('Calculating A derivatives',1)
C
              NDER3 = 1
              NDER4 = 0
              IAOCR = 1
              CALL CR2DRV(AOINT(1,IOFRST+1,1,IAOCR),CPINT,INDHER,
     &                    IODDHR(1,IHCBIT),IODDHR(1,IABXYZ),
     &                    IODDHR(1,ICCXYZ),LMNPWR,IPNTUV,
     &                    COORAB,EXPAB,CSQ,CCFBT,NCENTR(1,IOFRST+1,1),
     &                    NDER3,NDER4,NDIMD2,WRKCR2,LWKCR2,IPRINT)
 
              IF(IPRINT.GT.5) CALL HEADER('Calculating B derivatives',1)
C
              NDER3 = 0
              NDER4 = 1
              IAOCR = 4
              CALL CR2DRV(AOINT(1,IOFRST+1,1,IAOCR),CPINT,INDHER,
     &                    IODDHR(1,IHCBIT),IODDHR(1,IABXYZ),
     &                    IODDHR(1,ICCXYZ),LMNPWR,IPNTUV,
     &                    COORAB,EXPAB,CSQ,CCFBT,NCENTR(1,IOFRST+1,1),
     &                    NDER3,NDER4,NDIMD2,WRKCR2,LWKCR2,IPRINT)
            END IF
            IF(MIAU12) THEN
c             After [T1,r12] integrals derivatives are computed,
c             the logical variable U12INT is given its true value back
c             C. Villani (Uni Karlsruhe Feb 2005)
              U12INT=.TRUE.
              MIAU12=.FALSE.
            END IF
            if(miotst) then
c             C. Villani (Uni Karlsruhe Feb 2005)
              IF(UNDIFF) call mia2pr(aoint(1,iofrst+1,1,13),
     .                               ncct,khktab,khktcd,1,lupri)
              call miaprt(aoint(1,iofrst+1,1,1),iofrst,iprint,lupri)
           end if
         END IF
         IOFRST = IOFRST + MLTPZ
      END DO
C
      RETURN
      END
C  /* Deck pntuv */
      SUBROUTINE PNTUV(IPNTUV,LMNPWR,INDHER,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      DIMENSION IPNTUV(KC2MAX,0:NRDER,2,2), LMNPWR(KCKMAX,NHKMAX,3),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP)
#include <ccom.h>
#include <ericom.h>
#include <sphtrm.h>
#include <hertop.h>
#include <ibtfun.h>
C
      ICMPAB = 0
      MAXB = KCKTB
      DO ICMPA = 1, KCKTA
         LA = LMNPWR(ICMPA,NHKTA,1)
         MA = LMNPWR(ICMPA,NHKTA,2)
         NA = LMNPWR(ICMPA,NHKTA,3)
         IF (TCMPAB) MAXB = ICMPA
         DO ICMPB = 1, MAXB
            ICMPAB = ICMPAB + 1
            LB = LMNPWR(ICMPB,NHKTB,1)
            MB = LMNPWR(ICMPB,NHKTB,2)
            NB = LMNPWR(ICMPB,NHKTB,3)
            IPNTUV(ICMPAB,0,1,2) = INDHER(LA+LB,MA+MB,NA+NB)
            IF (GDER) THEN
               IPNTUV(ICMPAB,1,1,2) = INDHER(LA+LB+1,MA+MB,NA+NB)
               IPNTUV(ICMPAB,2,1,2) = INDHER(LA+LB,MA+MB+1,NA+NB)
               IPNTUV(ICMPAB,3,1,2) = INDHER(LA+LB,MA+MB,NA+NB+1)
            ELSE IF (BDER) THEN
               IPNTUV(ICMPAB,1,1,2) = INDHER(LA+LB,MA+MB+1,NA+NB+1)
               IPNTUV(ICMPAB,2,1,2) = INDHER(LA+LB+1,MA+MB,NA+NB+1)
               IPNTUV(ICMPAB,3,1,2) = INDHER(LA+LB+1,MA+MB+1,NA+NB)
            END IF 
         END DO
      END DO
C
      ICMPCD = 0
      MAXD = KCKTD
      DO ICMPC = 1, KCKTC
         LC = LMNPWR(ICMPC,NHKTC,1)
         MC = LMNPWR(ICMPC,NHKTC,2)
         NC = LMNPWR(ICMPC,NHKTC,3)
         IF (TCMPCD) MAXD = ICMPC
         DO ICMPD = 1, MAXD
            ICMPCD = ICMPCD + 1
            LD = LMNPWR(ICMPD,NHKTD,1)
            MD = LMNPWR(ICMPD,NHKTD,2)
            ND = LMNPWR(ICMPD,NHKTD,3)
            IPNTUV(ICMPCD,0,2,2) = INDHER(LC+LD,MC+MD,NC+ND)
            IF (GDER) THEN
               IPNTUV(ICMPCD,1,2,2) = INDHER(LC+LD+1,MC+MD,NC+ND)
               IPNTUV(ICMPCD,2,2,2) = INDHER(LC+LD,MC+MD+1,NC+ND)
               IPNTUV(ICMPCD,3,2,2) = INDHER(LC+LD,MC+MD,NC+ND+1)
            ELSE IF (BDER) THEN
               IPNTUV(ICMPCD,1,2,2) = INDHER(LC+LD,MC+MD+1,NC+ND+1)
               IPNTUV(ICMPCD,2,2,2) = INDHER(LC+LD+1,MC+MD,NC+ND+1)
               IPNTUV(ICMPCD,3,2,2) = INDHER(LC+LD+1,MC+MD+1,NC+ND)
            END IF
         END DO
      END DO
C
      IF (DOCART) THEN
         DO I = 0, NRDER
            CALL ICOPY(KHKTAB,IPNTUV(1,I,1,2),1,IPNTUV(1,I,1,1),1)
            CALL ICOPY(KHKTCD,IPNTUV(1,I,2,2),1,IPNTUV(1,I,2,1),1)
         END DO
      ELSE
         ICMPAB = 0
         MAXB = KHKTB
         DO ICMPA = 1, KHKTA
            IA = IDMAX(KCKTA,CSP(ISPADR(NHKTA)+ICMPA-1),KHKTA)
            IF (TKMPAB) MAXB = ICMPA
            DO ICMPB = 1, MAXB
               ICMPAB = ICMPAB + 1
               IB = IDMAX(KCKTB,CSP(ISPADR(NHKTB)+ICMPB-1),KHKTB)
               IF (TCMPAB) THEN
                  IAB = MAX(IA,IB)*(MAX(IA,IB) - 1)/2 + MIN(IA,IB)
               ELSE
                  IAB = (IA - 1)*KCKTB + IB
               END IF
               DO I = 0, NRDER 
                  IPNTUV(ICMPAB,I,1,1) = IPNTUV(IAB,I,1,2)
               END DO
            END DO
         END DO
C
         ICMPCD = 0
         MAXD = KHKTD
         DO ICMPC = 1, KHKTC
            IC = IDMAX(KCKTC,CSP(ISPADR(NHKTC)+ICMPC-1),KHKTC)
            IF (TKMPCD) MAXD = ICMPC
            DO ICMPD = 1, MAXD
               ICMPCD = ICMPCD + 1
               ID = IDMAX(KCKTD,CSP(ISPADR(NHKTD)+ICMPD-1),KHKTD)
               IF (TCMPCD) THEN
                  ICD = MAX(IC,ID)*(MAX(IC,ID) - 1)/2 + MIN(IC,ID)
               ELSE
                  ICD = (IC - 1)*KCKTD + ID
               END IF
               DO I = 0, NRDER
                  IPNTUV(ICMPCD,I,2,1) = IPNTUV(ICD,I,2,2)
               END DO
            END DO
         END DO
      END IF
C
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('IPNTUV in PNTUV',-1)
         WRITE (LUPRI,'(2X,A,I2,/)') ' NRDER ',NRDER
         DO J = 0,NRDER
            WRITE (LUPRI,'(15I5  )') (IPNTUV(I,J,1,2),I=1,KCKTAB)
            WRITE (LUPRI,'(15I5  )') (IPNTUV(I,J,1,1),I=1,KHKTAB)   
            WRITE (LUPRI,'(15I5  )') (IPNTUV(I,J,2,2),I=1,KCKTCD)
            WRITE (LUPRI,'(15I5,/)') (IPNTUV(I,J,2,1),I=1,KHKTCD)
         END DO
      END IF
C
      RETURN
      END
C  /* Deck setcao */
      SUBROUTINE SETCAO(IPRINT)
#include <implicit.h>
#include <ericom.h>
#include <eriao.h>
C
      IF (IPATH .EQ. 1) THEN
C
         IELCT1 = 1
         IELCT2 = 2
C
         NHKT1  = NHKTA
         NHKT2  = NHKTB
         NHKT3  = NHKTC
         NHKT4  = NHKTD
C
         JMAX1  = JMAXA
         JMAX2  = JMAXB
         JMAX3  = JMAXC
         JMAX4  = JMAXD
         JMAX12 = JMAXAB
         JMAX34 = JMAXCD
C
         NTUV12 = NTUVAB
         NTUV34 = NTUVCD
C
         KHKT1  = KHKTA
         KHKT2  = KHKTB
         KHKT3  = KHKTC
         KHKT4  = KHKTD
         KHKT12 = KHKTAB
         KHKT34 = KHKTCD
C
         KCKT1  = KCKTA
         KCKT2  = KCKTB
         KCKT3  = KCKTC
         KCKT4  = KCKTD
         KCKT12 = KCKTAB
         KCKT34 = KCKTCD
C
         SPHR1  = SPHRA
         SPHR2  = SPHRB
         SPHR3  = SPHRC
         SPHR4  = SPHRD
         SPHR12 = SPHRAB
         SPHR34 = SPHRCD
C
         TKMP12 = TKMPAB
         TKMP34 = TKMPCD
         TCMP12 = TCMPAB
         TCMP34 = TCMPCD
C
         NPRF1  = NPRFA
         NPRF2  = NPRFB
         NPRF3  = NPRFC
         NPRF4  = NPRFD
         NPRF12 = NPRFAB
         NPRF34 = NPRFCD
C
         NCTF1  = NCTFA
         NCTF2  = NCTFB
         NCTF3  = NCTFC
         NCTF4  = NCTFD
         NCTF12 = NCTFAB
         NCTF34 = NCTFCD
C
         ICMAT1  = ICMATA
         ICMAT2  = ICMATB
         ICMAT3  = ICMATC
         ICMAT4  = ICMATD
C
         GCON1  = GCONA
         GCON2  = GCONB
         GCON3  = GCONC
         GCON4  = GCOND
         GCON12 = GCONAB
         GCON34 = GCONCD
C
         NCNT12 = NCNTAB
         NCNT34 = NCNTCD
C
         NODS12 = NODSAB
         NODS34 = NODSCD
C
         KODS12 = KODSAB
         KODS34 = KODSCD
C
         NPPPP  = NPPX
         NPPCC  = NPCX
         NCCPP  = NCPX
         NCCCC  = NCCX
C
         NPP12 = NPPAB
         NPP34 = NPPCD
C
         I12XYZ  = IABXYZ
         I120(1) = IAB0(1)
         I120(2) = IAB0(2)
         I120(3) = IAB0(3)
C
         I34XYZ  = ICDXYZ
         I340(1) = ICD0(1)
         I340(2) = ICD0(2)
         I340(3) = ICD0(3)
C
         IHHBIT = IHHXYZ
         IHCBIT = IHCXYZ
         ICHBIT = ICHXYZ
         ICCBIT = ICCXYZ
      ELSE
C
         IELCT1 = 2
         IELCT2 = 1
C
         NHKT1  = NHKTC
         NHKT2  = NHKTD
         NHKT3  = NHKTA
         NHKT4  = NHKTB
C
         JMAX1  = JMAXC
         JMAX2  = JMAXD
         JMAX3  = JMAXA
         JMAX4  = JMAXB
         JMAX12 = JMAXCD
         JMAX34 = JMAXAB
C
         NTUV12 = NTUVCD
         NTUV34 = NTUVAB
C
         KHKT1  = KHKTC
         KHKT2  = KHKTD
         KHKT3  = KHKTA
         KHKT4  = KHKTB
         KHKT12 = KHKTCD
         KHKT34 = KHKTAB
C
         KCKT1  = KCKTC
         KCKT2  = KCKTD
         KCKT3  = KCKTA
         KCKT4  = KCKTB
         KCKT12 = KCKTCD
         KCKT34 = KCKTAB
C
         SPHR1  = SPHRC
         SPHR2  = SPHRD
         SPHR3  = SPHRA
         SPHR4  = SPHRB
         SPHR12 = SPHRCD
         SPHR34 = SPHRAB
C
         TKMP12 = TKMPCD
         TKMP34 = TKMPAB
         TCMP12 = TCMPCD
         TCMP34 = TCMPAB
C
         NPRF1  = NPRFC
         NPRF2  = NPRFD
         NPRF3  = NPRFA
         NPRF4  = NPRFB
         NPRF12 = NPRFCD
         NPRF34 = NPRFAB
C
         NCTF1  = NCTFC
         NCTF2  = NCTFD
         NCTF3  = NCTFA
         NCTF4  = NCTFB
         NCTF12 = NCTFCD
         NCTF34 = NCTFAB
C
         ICMAT1  = ICMATC
         ICMAT2  = ICMATD
         ICMAT3  = ICMATA
         ICMAT4  = ICMATB
C
         GCON1  = GCONC
         GCON2  = GCOND
         GCON3  = GCONA
         GCON4  = GCONB
         GCON12 = GCONCD
         GCON34 = GCONAB
C
         NCNT12 = NCNTCD
         NCNT34 = NCNTAB
C
         NODS12 = NODSCD
         NODS34 = NODSAB
C
         KODS12 = KODSCD
         KODS34 = KODSAB
C
         NPPPP  = NPPX
         NPPCC  = NCPX
         NCCPP  = NPCX
         NCCCC  = NCCX
C
         NPP12 = NPPCD
         NPP34 = NPPAB
C
         I12XYZ  = ICDXYZ
         I120(1) = ICD0(1)
         I120(2) = ICD0(2)
         I120(3) = ICD0(3)
C
         I34XYZ  = IABXYZ
         I340(1) = IAB0(1)
         I340(2) = IAB0(2)
         I340(3) = IAB0(3)
C
         IHHBIT = IHHXYZ
         IHCBIT = ICHXYZ
         ICHBIT = IHCXYZ
         ICCBIT = ICCXYZ
      END IF
      RETURN
      END
C  /* Deck lmnprp */
      SUBROUTINE LMNPRP(LMNPWR,LMNSYM,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      DIMENSION LMNPWR(KCKMAX,NHKMAX,3), LMNSYM(KCKMAX,NHKMAX,2,0:7)
#include <ericom.h>
#include <eriao.h>
#include <ccom.h>
#include <sphtrm.h>
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from LMNPRP','*',103)
C
      DO I = 1, NHKMAX
         KCKTX = I*(I + 1)/2
         KHKTX = 2*I - 1
         IF (DOCART) KHKTX = KCKTX
C
C        Angular powers
C
         CALL LMNVAL(I,KCKTX,LMNPWR(1,I,1),LMNPWR(1,I,2),LMNPWR(1,I,3))
C
C        Local symmetries for Cartesian functions
C
         DO ISYM = 0, 7
         DO J = 1, KCKTX
            IBITS =   MOD(LMNPWR(J,I,1),2)
     &            + 2*MOD(LMNPWR(J,I,2),2)
     &            + 4*MOD(LMNPWR(J,I,3),2)
            LMNSYM(J,I,1,ISYM) = IBTAND(ISYM,IBITS)
         END DO
         END DO
C
C        Local symmetries for spherical functions
C
         DO ISYM = 0, 7
         DO J = 1, KHKTX
            IF (DOCART) THEN
               INDMAX = J
            ELSE
               INDMAX = IDMAX(KCKTX,CSP(ISPADR(I)+J-1),KHKTX)
            END IF
            LMNSYM(J,I,2,ISYM) = LMNSYM(INDMAX,I,1,ISYM)
         END DO
         END DO
      END DO
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('LMNPWR from LMNPRP',-1)
         DO I = 1, NHKMAX
         DO J = 1, I*(I + 1)/2
            WRITE (LUPRI,'(5X,2I5,5X,3I5)') I,J,(LMNPWR(J,I,K),K=1,3)
         END DO
         END DO
C
         CALL HEADER('LMNSYM (Cartesian) from LMNPRP',-1)
         DO I = 1, NHKMAX
         DO J = 1, I*(I + 1)/2
            WRITE (LUPRI,'(5X,2I5,5X,8I5)') I,J,(LMNSYM(J,I,1,K),K=0,7)
         END DO
         END DO
C
         CALL HEADER('LMNSYM (Spherical) from LMNPRP',-1)
         DO I = 1, NHKMAX
            KHKMX = 2*I - 1
            IF (DOCART) KHKMX = I*(I + 1)/2
            DO J = 1, KHKMX
               WRITE(LUPRI,'(5X,2I5,5X,8I5)')I,J,(LMNSYM(J,I,2,K),K=0,7)
            END DO
         END DO
      END IF
      RETURN
      END
C  /* Deck npairs */
      INTEGER FUNCTION NPAIRS(J)
      INTEGER J
      NPAIRS = J*(J+1)/2
      RETURN
      END
C  /* Deck kcdim */
      INTEGER FUNCTION KCDIM(JMAX1,JMAX2)
#include <implicit.h>
      ITRI1(I) = (I + 1)*(I + 2)/2
      ITRI2(I) = (I + 1)*(I + 2)*(I + 3)/6
C
      NMAX = 0
      NOFF = ITRI2(MAX(JMAX1,JMAX2) - 1)
      DO J = 0, MIN(JMAX1,JMAX2)
         NDIM1 = ITRI1(J)
         NDIM2 = ITRI2(JMAX1 + JMAX2 - J) - NOFF
         NMAX = MAX(NMAX,NDIM1*NDIM2)
      END DO
      KCDIM = NMAX
C
      RETURN
      END
C  /* Deck mempqb */
      SUBROUTINE MEMPQB(MEMORY,NALLOI,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <mxcent.h>
      LOGICAL DOREC, FRSTAO
#include <cbieri.h>
#include <ericom.h>
#include <eridst.h>
#include <eribuf.h>
#include <symmet.h>
#include <r12int.h>
C
      MEMI(I) = I
      MEMR(I) = IRAT*I
C------------------------------------
C
C     ODCPAR
C     ======
C
      MEMPAR = MEMI(NITPQ*NPRFPQ)
      NAIPAR = 2
C
C     ODCPAS
C     ======
C
#if defined (VAR_CRY)
      MEMPAS = MEMR(27*NPRFPQ + 6*NPRFPQ + 3*NPRFAB + 3*NPRFCD)
#else
      MEMPAS = MEMR(27*NPRFPQ)
#endif
      NAIPAS = 0
C
C     INTDRV
C     ======
C
      LRJ000 = NPRFPQ*(JMAX + 1)
      LHRINT = NPRFPQ*NTUV
      LCPIN1 = 0
      LCPIN2 = 0
      IF (DOPTH1) LCPIN1 = NCTFAB*NPRFCD*NTUVCD*KHKTAB
      IF (DOPTH2) LCPIN2 = NCTFCD*NPRFAB*NTUVAB*KHKTCD
      LCPINT = MAX(LCPIN1,LCPIN2) 
      IF (OLDCR1) THEN
         LCPMAX = LCPINT 
      ELSE
         KCREC1 = KCDIM(JMAXA,JMAXB)
         LCPMAX = NCTFAB*NPRFCD*NTUVCD*KCREC1
      END IF
      LAOINT = NCTFPQ*KHKTAB*KHKTCD
C
      NGAM = MEMR(LHRINT + LRJ000)
      NHER = MEMR(LHRINT + LRJ000)
      NCR1 = MEMR(MAX(LAOINT,LHRINT) + LCPMAX)
      NCR2 = MEMR(MAX(LAOINT,LHRINT) + LCPINT)
      IF (R12INT .OR. U12INT) THEN
C        Space for R12 method (WK/UniKA/04-11-2002).
         NQ52 = MEMR(2*LHRINT + NPRFPQ)
         IF (U12INT) THEN
            NCR1 = MEMR(MAX(LAOINT,2*LHRINT) + LCPMAX)
            NCR2 = MEMR(MAX(LAOINT,2*LHRINT) + LCPINT)
         END IF
      ELSE
         NQ52 = 0
      END IF
      NOUT = MEMR(LAOINT)
      IF (MLTPX .GT. 1) NOUT = NOUT + NOUT
C
C     **************
C     *** ERIGAM ***
C     **************
C
C     ... to be added to NGAM
C
C     ERIGAM
C     ======
C
      NREAL  = 2
      NINTG  = 1
      MEMGAM = MEMR(NREAL*NPRFPQ) + MEMI(NINTG*NPRFPQ)
      NAIGAM = 1
C
C     GETGAM
C     ======
C
      NREAL  = 4 + (JMAX + 1)
      NINTG  = 3
      MEMGGM = MEMR(NREAL*NPRFPQ) + MEMI(NINTG*NPRFPQ)
      NAIGGM = 1
C
C     **************
C     *** ERITUV ***
C     **************
C
C     ... to be added to NHER
C
C     ERITUV
C     ======
C
      NREAL  = 0
      IF (JMAX .GT. 0) NREAL  = NTUV
      MEMTUV = MEMR(NREAL*NPRFPQ)
C
C     **************
C     *** CR1DRV ***
C     **************
C
C     ... to be added to NCR1
C
C     CR1TWO
C     ======
C
      IF (OLDCR1) THEN
         IF (DOPTH1) THEN
            NREAL1 = MEMCR1(KHKTAB,KCKTAB,KHKTA,KHKTB,
     &                      JMAXAB,JMAXA,JMAXB,
     &                      SPHRAB,SPHRA,SPHRB,
     &                      NPRFCD,NTUVCD,NPRFPQ)
         ELSE
            NREAL1 = 0
         END IF
         IF (DOPTH2) THEN
            NREAL2 = MEMCR1(KHKTCD,KCKTCD,KHKTC,KHKTD,
     &                      JMAXCD,JMAXC,JMAXD,
     &                      SPHRCD,SPHRC,SPHRD,
     &                      NPRFAB,NTUVAB,NPRFPQ)
         ELSE
            NREAL2 = 0
         END IF
         MEMC1T = MEMR(MAX(NREAL1,NREAL2))
      ELSE
         MEMC1T = 0
      END IF
C
C     CR1ONE
C     ======
C
      IF (.NOT.OLDCR1) THEN
         LCCONT = 0
         LXDIF  = 0
         LCSINT = 0
         LECOEF = 0
         LHPI   = 0
         LPA    = 0
         LPB    = 0
         LEUV   = 0
         LETUV  = 0
         LHCPRM = 0
         LFIRST = 0
C
         DOREC  = MIN(NHKTA,NHKTB) .NE. 1
         FRSTAO = NHKTA .LE. NHKTB
C
         IF (KHKTAB .GT. 1 .OR. U12INT) THEN
C           No symmetry for [T1,r12] integrals (WK/UniKA/04-11-2002).
            IF (DOREC) THEN
               LCCONT = NPRFCD*NTUVCD*KCREC1
               LXDIF  = 3*NPRFCD
            END IF
C           IF (SPHRA .AND. SPHRB) THEN
C              LCSINT = NPRFCD*NTUVCD*KHKTB
C           END IF
C
            LECOEF = 3*NPRFPQ*(JMAXAB + 1)*(JMAXAB + 1)
            IF (U12INT) LECOEF = 2*LECOEF
            LHPI   =   NPRFPQ
            IF (FRSTAO) THEN
               LPB = 3*NPRFPQ
            ELSE
               LPA = 3*NPRFPQ
            END IF
C
            LEUV   = NPRFPQ
            LETUV  = NPRFPQ
C
            LHCPRM = NPRFPQ*NTUVCD
            LFIRST = LHCPRM
         END IF
C
         NREAL = LCCONT
     &          + MAX( LXDIF,
     &                 LCSINT,
     &                 LECOEF + 2*NPRFPQ + MAX(
     &                              LHPI+LPA+LPB + 2*NPRFPQ,
     &                              LEUV+LETUV+LHCPRM+LFIRST
     &                             )
     &               )
         MEMC1O = MEMR(NREAL)
      ELSE
         MEMC1O = 0
      END IF
C
C     **************
C     *** CR2DRV ***
C     **************
C
C     ... to be added to NCR2
C
C     CR2TWO
C     ======
C
      IF (.TRUE.) THEN
         IF (DOPTH1) THEN
            NREAL1 = MEMCR2(KHKTAB,KHKTCD,KHKTD,KCKTCD,
     &                      JMAXCD,JMAXC,JMAXD,
     &                      SPHRCD,SPHRC,SPHRD,
     &                      NCTFAB,NPRFCD,NCTFPQ)
         ELSE
            NREAL1 = 0
         END IF
         IF (DOPTH2) THEN
            NREAL2 = MEMCR2(KHKTCD,KHKTAB,KHKTB,KCKTAB,
     &                      JMAXAB,JMAXA,JMAXB,
     &                      SPHRAB,SPHRA,SPHRB,
     &                      NCTFCD,NPRFAB,NCTFPQ)
         ELSE
            NREAL2 = 0
         END IF
         MEMC2T = MEMR(MAX(NREAL1,NREAL2))
      ELSE
         MEMC2T = 0
      END IF
C
C     **************
C     *** ERIOUT ***
C     **************
C
C     ... to be added to NOUT
C
C     ERIOUT
C     ======
C
C
      LBIN = NCTFPQ*KHKTA*KHKTB*KHKTC*KHKTD
      IF (WRTINT) THEN
         IF (DODIST) THEN
            NREAL  =             LBIN
            NINTG  = (2*NIBUF+1)*LBIN + 2*NCTFPQ
     &             + 4*(MAXREP+1)*NCTFPQ
            MEMOUT = MEMR(NREAL) + MEMI(NINTG)
            NAIOUT = 5
         ELSE
            NREAL  =   LBIN
            NINTG  = 2*LBIN*NIBUF
            MEMOUT = MEMR(NREAL) + MEMI(NINTG)
            NAIOUT = 2
         END IF
      ELSE
         NREAL  =   LBIN
         NINTG  = 8*LBIN
         MEMOUT = MEMR(NREAL) + MEMI(NINTG)
         NAIOUT = 2
      END IF
C
C     Total memory
C     ============
C
      MEM_CR1  = NCR1 + MAX(MEMC1T,MEMC1O)
      MEMORY = MEMPAR + MEMPAS + MAX(
     &                               NQ52,
     &                               NGAM + MEMGAM + MEMGGM,
     &                               NHER + MEMTUV,
     &                               MEM_CR1,
     &                               NCR2 + MEMC2T,
     &                               NOUT + MEMOUT
     &                              )
      MEMORY = MLTPZ*MEMORY +  MEMR(12*NPRFPQ)
      NALLOI = NAIPAR + NAIPAS + NAIGAM + NAIGGM + NAIOUT
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Output from MEMPQB',-1)
         WRITE (LUPRI,'(2X,A,I12)') 'MEMPAR',MEMPAR
         WRITE (LUPRI,'(2X,A,I12)') 'MEMPAS',MEMPAS
         WRITE (LUPRI,'(2X,A,I12)') 'NQ52  ',NQ52
         WRITE (LUPRI,'(2X,A,I12)') 'NGAM  ',NGAM
         WRITE (LUPRI,'(2X,A,I12)') 'MEMGAM',MEMGAM
         WRITE (LUPRI,'(2X,A,I12)') 'MEMGGM',MEMGGM
         WRITE (LUPRI,'(2X,A,I12)') 'NHER  ',NHER
         WRITE (LUPRI,'(2X,A,I12)') 'MEMTUV',MEMTUV
         WRITE (LUPRI,'(2X,A,I12)') 'NCR1  ',NCR1
         WRITE (LUPRI,'(2X,A,I12)') 'MEMC1T',MEMC1T
         WRITE (LUPRI,'(2X,A,I12)') 'MEMC1O',MEMC1O
         WRITE (LUPRI,'(2X,A,I12)') 'NCR2  ',NCR2
         WRITE (LUPRI,'(2X,A,I12)') 'MEMC2T',MEMC2T
         WRITE (LUPRI,'(2X,A,I12)') 'NOUT  ',NOUT
         WRITE (LUPRI,'(2X,A,I12)') 'MEMOUT',MEMOUT
         WRITE (LUPRI,'(2X,A,I12)') 'MLTPZ ',MLTPZ
         WRITE (LUPRI,'(2X,A,I12)') 'NPRFPQ',NPRFPQ
         WRITE (LUPRI,'(A,I12)')  '**MEMORY',MEMORY
         WRITE (LUPRI,'(A,I12)')  '**NALLOI',NALLOI
      END IF
C
      RETURN
      END
C  /* Deck memcr1 */
      INTEGER FUNCTION MEMCR1(KHKT12,KCKT12,KHKT1,KHKT2,
     &        JMAX12,JMAX1,JMAX2,
     &        SPHR12,SPHR1,SPHR2, NPRF34,NTUV34,NPRFPQ)
C
#include <implicit.h>
#include <r12int.h>
C
      LOGICAL SPHR12, SPHR1, SPHR2
C
      IF (KHKT12 .GT. 1 .OR. U12INT) THEN
C        No symmetry for [T1,r12] integrals (WK/UniKA/04-11-2002).
         IF (SPHR12) THEN
            LCCONT = NPRF34*NTUV34*KCKT12
            LFCINT = NTUV34*KHKT12
         ELSE
            LCCONT = 0
            LFCINT = 0
         END IF
         IF (SPHR1 .AND. SPHR2) THEN
            LCSINT = NPRF34*NTUV34*KHKT2
            LFCSNT = NTUV34*KHKT2
         ELSE
            LCSINT = 0
            LFCSNT = 0
         END IF
C
         LECOEF = 3*NPRFPQ*(JMAX12 + 1)*(JMAX1 + 1)*(JMAX2+1)
         LHPI   =   NPRFPQ
         LPA    = 3*NPRFPQ
         LPB    = 3*NPRFPQ
C
         LEUV   = NPRFPQ
         LETUV  = NPRFPQ
C
         LHCPRM = NPRFPQ*NTUV34
         MEM_A  = MAX(LHPI+LPA+LPB,LEUV+LETUV+LHCPRM+NTUV34)
         MEMCR1 = LCCONT
     &       + MAX(LCSINT + LFCSNT + LFCINT, LECOEF + MEM_A)
      ELSE
         MEMCR1 = NTUV34
      END IF
C
      RETURN
      END
C  /* Deck memcr2 */
      INTEGER FUNCTION MEMCR2(KHKT12,KHKT34,KHKT4,KCKT34,
     &                 JMAX34,JMAX3,JMAX4,
     &                 SPHR34,SPHR3,SPHR4, NCTF12,NPRF34,NCTFPQ)
C
#include <implicit.h>
      LOGICAL SPHR34,SPHR3,SPHR4
C
      IF (KHKT34 .GT. 1) THEN
         IF (SPHR34) THEN
            LCCONT = NCTFPQ*KHKT12*KCKT34
            LFAONT = KHKT12*KHKT34
         ELSE
            LCCONT = 0
            LFAONT = 0
         END IF
         IF (SPHR3 .AND. SPHR4) THEN
            LCSINT = NCTFPQ*KHKT12*KHKT4
            LFCSNT = KHKT12*KHKT4
         ELSE
            LCSINT = 0
            LFCSNT = 0
         END IF
C
         LECOEF = 3*NCTF12*NPRF34*(JMAX34+1)*(JMAX3+1)*(JMAX4+1)
         LHPI   =   NCTF12*NPRF34
         LPC    = 3*NCTF12*NPRF34
         LPD    = 3*NCTF12*NPRF34
C
         LEUV   = NCTF12*NPRF34
         LETUV  = NCTF12*NPRF34
C
         LCCPRM = NCTF12*NPRF34*KHKT12
         MEM_A  = MAX(LHPI+LPC+LPD,LEUV+LETUV+LCCPRM+KHKT12)
         MEMCR2 = LCCONT + MAX(LCSINT + LFCSNT + LFAONT,LECOEF + MEM_A)
      ELSE
         MEMCR2 = KHKT12
      END IF
C
      RETURN
      END

      subroutine miaprt(aoint,iof,iprint,iu)
#include <implicit.h>
#include <ericom.h>
#include <eriao.h>
      integer iprint
      dimension aoint(ncct,khktab,khktcd,*)

c     Not sure, here
      if(ncct.ne.ncccc) stop "ncct.ne.ncccc, unexpected"

c     if(iprint.ge.99) then
        call header('Final spherical integrals - MIAPRT',-1)
        write(iu,*) " "
        write(iu,*) "==========================================="
        write(iu,*) "C derivatives: "
        call mia2pr(aoint(1,iof+1,1,7),ncccc,khktab,khktcd,3,iu)
        write(iu,*) " "
        write(iu,*) "==========================================="
        write(iu,*) "D derivatives: "
        call mia2pr(aoint(1,iof+1,1,10),ncccc,khktab,khktcd,3,iu)
        write(iu,*) " "
        write(iu,*) "==========================================="
        write(iu,*) "A derivatives: "
        call mia2pr(aoint(1,iof+1,1,1),ncccc,khktab,khktcd,3,iu)
        write(iu,*) " "
        write(iu,*) "==========================================="
        write(iu,*) "B derivatives: "
        call mia2pr(aoint(1,iof+1,1,4),ncccc,khktab,khktcd,3,iu)
        write(iu,*) "==========================================="
        write(iu,*) "Checking sum of derivatives... "
c     end if
      call mia2cps(aoint(1,iof+1,1,1),ncccc,khktab,khktcd,3,iu)

      return
      end
c
      subroutine mia2pr(aoint,nccmio,khktab,khktcd,ndimd,u)
      implicit none
      integer nccmio, khktab, khktcd, ndimd, u, i, j, k, icc
      double precision aoint(nccmio,khktab,khktcd,ndimd)
      character*1 xyz(3)

      xyz(1)='x'
      xyz(2)='y'
      xyz(3)='z'

      if(ndimd.gt.1) then
        do k = 1, ndimd
          write(u, '(/1x,a1,a)') xyz(k), ' derivatives '
          do i = 1, khktab
            do j = 1, khktcd
              if(nccmio.eq.1) then
                write(u,505) i, j, aoint(1,i,j,k)
              else
                stop 'may not work...'
                do icc=1,nccmio
                  write(u,505) i, j, aoint(1+icc-1,i,j,k)
                end do
              end if
            end do
          end do
        end do
      else
        do i = 1, khktab
          do j = 1, khktcd
            if(nccmio.eq.1) then
              write(u,505) i, j, aoint(1,i,j,1)
            else
              stop 'may not work...'
              do icc=1,nccmio
                write(u,505) i, j, aoint(1+icc-1,i,j,1)
              end do
            end if
          end do
        end do
      end if

      return
 505  format(2i4,f22.14)
      end

      subroutine miach(v,n)
      implicit none
      integer i, n
      double precision v(n)
c
c     The half-computed B derivatives integrals have the wrong sign
c     it would be more efficient to change sign inside the eri2cft,
c     however this additional call should not affect efficiency much.
c
c     C. Villani (Uni Karlsruhe Feb 2005)
c
      do i=1,n
        v(i)=-v(i)
      end do
      return
      end

      subroutine mia2cps(aoint,nccmio,khktab,khktcd,ndimd,u)
      implicit none
      logical prt
      integer nccmio, khktab, khktcd, ndimd, u, i, j, k, icc
      double precision aoint(nccmio,khktab,khktcd,ndimd), s, thr
      character*1 xyz(3)
c
c     checks translational invariancy of the integrals derivatives
c
c     C. Villani, Uni Karlsruhe, Feb 2005

      xyz(1)='x'
      xyz(2)='y'
      xyz(3)='z'

      prt=.false.
      thr=0.00000001d0

      if(ndimd.gt.1) then
        do k = 1, ndimd
          if(prt) write(u, '(/1x,a1,a)') xyz(k), ' derivatives '
          do i = 1, khktab
            do j = 1, khktcd
              s=aoint(1,i,j,k)+aoint(1,i,j,k+3)+
     .          aoint(1,i,j,k+6)+aoint(1,i,j,k+9)
              if(nccmio.eq.1) then
                if(prt) write(u,505) i, j, s
                if(dabs(s).gt.thr) then
                  print *, i, j, s, 'non inv'
                  write(u,*) i, j, s, 'non inv'
                  call quit('Integral deriv. are trasl. non invariant')
                end if
              else
c               not sure, here
                stop 'Integral deriv. ncc>1, may not work...'
              end if
            end do
          end do
        end do
      else
        do i = 1, khktab
          do j = 1, khktcd
            s=aoint(1,i,j,k)+aoint(1,i,j,k+3)+
     .        aoint(1,i,j,k+6)+aoint(1,i,j,k+9)
            if(nccmio.eq.1) then
              if(prt) write(u,505) i, j, s
              if(dabs(s).gt.thr) then
                print *, i, j, s, 'non inv'
                write(u,*) i, j, s, 'non inv'
                call quit('Integral deriv. are trasl. non invariant')
              end if
            else
c               not sure, here
              stop 'Integral deriv. ncc>1, may not work...'
            end if
          end do
        end do
      end if

      return
 505  format(2i4,f22.14)
      end
