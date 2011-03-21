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
C  /* Deck hsoctl */
      SUBROUTINE HSOCTL (WORD,GP,CMO,UDV,PV,XINDX,ANTSYM,WRK,LWRK)
C
C Olav Vahtras
C Jun 20, 1990
C
C Driver routine for the construction of spin-orbit property vector
C
C Input: Spin-orbit component (as read from LUPROP), WORD
C        MO coefficients,CMO
C        First order reduced density matrix, UDV
C        Second order reduced density matrix, PV
C        Index vector, XINDX
C
C Output: Spin-orbit property vector returned in first elements of WRK
C
#include <implicit.h>
#include <dummy.h>
      DIMENSION GP(*),CMO(*),UDV(*),PV(*),XINDX(*),WRK(LWRK),NEEDSO(6)
      CHARACTER*8 WORD
      DATA NEEDSO/1,1,1,1,1,0/
C
      PARAMETER(D0=0D0, D1=1D0, DM1=-1D0)
      LOGICAL SO2TRA, MOEXIS, FOPEN, DO1, DO2, LORB, LCON
      CHARACTER*8 LABEL
      LOGICAL FNDLAB
      EXTERNAL FNDLAB
C
C Used from common blocks:
C  MULD2H,ISMO,ISW,IOBTYP,IDBTYP,NISHT,NORBT,NASHT,NNASHX
C  INFORB: MULD2H,NISHT,NASHT,NORBT,NNASHX
C  INFIND: ISW,IOBTYP
C  TRHSO : ILXYZ, KSYMSO, OLDTRA
C  INFHSO: IPRHSO, TESTZY, DOSO1, DOSO2
C  INFHYP: HYPCAL
C  INFSMO: SOMOM
C  INFPP: EXMOM
C  INFRSP: ???,SOPPA
C
C-- common blocks:
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <inftap.h>
#include <infind.h>
#include <infvar.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <infpri.h>
#include <trhso.h>
#include <infhso.h>
#include <infhyp.h>
#include <infsmo.h>
#include <infpp.h>
#include <iratdef.h>
#include <codata.h>
C
      CALL QENTER('HSOCTL')
C
      IF (SOPPA) CALL QUIT('HSOCTL: SOPPA not implemented yet!')
C
      CALL DZERO(GP,KZYVAR)
      HSOFAC = ALPHAC**2 / 4
      IF (WORD(2:2).EQ.'1') THEN
         CALL PRPCTL(WORD,GP,CMO,UDV,PV,XINDX,ANTSYM,WRK,LWRK)
         CALL DSCAL(KZYVAR,HSOFAC,GP,1)
         GO TO 9999
      ELSE IF (WORD(2:2).EQ.'2') THEN
         DO1 = .FALSE.
         DO2 = .TRUE.
      ELSE IF (WORD(2:2).EQ.' ') THEN
         DO1 = DOSO1
         DO2 = DOSO2
      ELSE
         WRITE(LUPRI,*) 'Programming error: WORD is not valid: ',WORD
         CALL QUIT('Programming error: WORD is not valid')
      END IF
      ANTSYM = 1.0D0
C
      IF (IPRHSO.GT.0) THEN
         CALL TIMER('START ',HSOSTA,HSOTIM)
         CALL HEADER('Output from HSOCTL',-1)
         WRITE(LUPRI,'(/2A,3X,A)')
     *   ' Spin-orbit property vector calculation',
     *   ' component = ',WORD
         WRITE(LUPRI,'(/A,I5)')' Print level in HSOCTL: ',IPRHSO
         IF (TESTZY) WRITE(LUPRI,'(/A)')
     *' Z and Y parts of configurational property vector explicitly'
         IF (.NOT.DO1) WRITE(LUPRI,'(/A)')
     *' Skip one-electron spin-orbit contributions'
         IF (.NOT.DO2) WRITE(LUPRI,'(/A)')
     *' Skip two-electron spin-orbit contributions'
      END IF
      LORB = KZWOPT.GT.0
      IF (KSYMOP.EQ.1) THEN
         LCON = KZCONF.GT.1
      ELSE
         LCON = KZCONF.GT.0
      END IF
C
C Check if gradient is on file
C
      LUHSO = -1
      CALL GPOPEN(LUHSO,'HSOGRAD','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUHSO
      IF (FNDLAB(WORD,LUHSO)) THEN
         CALL READT(LUHSO,KZYVAR,GP)
         CALL GPCLOSE(LUHSO,'KEEP')
         GO TO 9999
      END IF
C
C ALLOCATE WORK SPACE do it smarter here
C
      KTEST  = 1
      WORK(KTEST) = -999.9D0
      KFC    = KTEST + 1
      IF (LORB) THEN
         KFV    = KFC    + NORBT*NORBT
         KQA    = KFV    + NORBT*NORBT
         KQB    = KQA    + NORBT*NASHT
         KPV12  = KQB    + NORBT*NASHT
         KPV21  = KPV12  + N2ASHX
         KH1    = KPV21  + N2ASHX
      ELSE
         KH1    = KFC    + NORBT*NORBT
         KFV    = KTEST
         KQA    = KTEST
         KQB    = KTEST
         KPV12  = KTEST
         KPV21  = KTEST
      END IF
      KH2    = KH1 + NORBT*NORBT
      IF (LCON) THEN
         KH2AC  = KH2   + NORBT*NORBT
         KWRK1  = KH2AC + NNASHX*NNASHX
      ELSE
         KH2AC  = KTEST
         KWRK1  = KH2   + NORBT*NORBT
      END IF
      LWRK1 = LWRK - KWRK1 + 1
      IF (LWRK1.LT.0) CALL ERRWRK('HSOCTL',KWRK1-1,LWRK)
      CALL DZERO(WRK,KWRK1-1)
C
      IF (X2GRAD) THEN
         KGP2   = 1
         KMJWOP = KGP2 + KZYVAR
         KWRK2  = KMJWOP + (16*MAXWOP + 1)/IRAT
         LWRK2  = LWRK - KWRK2 + 1
         IF (LWRK2.LT.0) CALL ERRWRK('HSOCTL',KWRK2-1,LWRK)
      END IF
C
C Set spin-orbit component ILXYZ for TRHSO
C
      IF (WORD(1:1).EQ.'X') THEN
         ILXYZ = 1
      ELSE IF (WORD(1:1).EQ.'Y') THEN
         ILXYZ = 2
      ELSE IF (WORD(1:1).EQ.'Z') THEN
         ILXYZ = 3
      ELSE
         WRITE(LUPRI,'(/2A)') 'Wrong property in HSOCTL, WORD = ',WORD
         CALL QUIT(' ')
      END IF
C
      IF (IPRHSO.GT.0) THEN
         IF (IPRHSO.GT.20 .AND. LORB .AND. NASHT.GT.0) THEN
            WRITE(LUPRI,'(/A)')
     *   ' First order density matrix:'
            CALL OUTPUT(UDV,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
            WRITE(LUPRI,'(/A)')
     *   ' Second order density matrix:'
            CALL PRIAC2(PV,NASHT,LUPRI)
         END IF
      END IF
C
C
C Read AO one-electron property integrals and transform to MO basis.
C
      IF (DO1) THEN
         LABEL = WORD
         LABEL(2:2) = '1'
         CALL PRPGET (LABEL,CMO,WRK(KH1),ANTSYM,WRK(KWRK1),LWRK1)
      END IF
C
C Print atomic and molecular property integrals, if desired
C
      IF (IPRHSO.GT.20 .AND. DO1) THEN
         WRITE (LUPRI,'(/2A)')' Atomic property integrals:', LABEL
         CALL OUTPUT(WRK(KWRK1),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE (LUPRI,'(/2A)')' Molecular property integrals:', LABEL
         CALL OUTPUT(WRK(KH1),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C Transform spin-orbit two-electron integrals, we need sec-occ in linear
C response and sec-sec in quadratic response
C
      IF (HYPCAL.OR.SOMOM.OR.EXMOM) THEN
         ITRLSO = 2
      ELSE
         ITRLSO = 1
      END IF
      KSYMSO = KSYMOP
C defined in RSPSET (920922-ov)
      IF (.NOT.DO2) GOTO 95
      SO2TRA = .TRUE.
C
C If .OLDTRA has been specified check that MO2SOINT exists
C If it does not exist transform AO2SOINT as usual
C
      IF (OLDTRA) THEN
         CALL GPINQ('MOHSOINT','EXIST',MOEXIS)
         IF (MOEXIS) THEN
            SO2TRA = .FALSE.
         ELSE
            WRITE(LUPRI,'(/3A/A)') ' WARNING: Expected transformed ',
     *                       ' two-electron spin-orbit integrals not',
     *                       ' found.' ,
     *                       ' - file MOHSOINT does not exist'
            NWARN = NWARN + 1
         END IF
      END IF
      IF (SO2TRA) THEN
         IF (IPRHSO.GT.0) THEN
            CALL TIMER('START ',TRASTA,TRATIM)
         END IF
         CALL GPOPEN(LUAHSO,'AO2SOINT','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         CALL TRAHSO(ITRLSO,CMO,WRK(KWRK1),LWRK1)
         CALL GPCLOSE(LUAHSO,'KEEP')
         IF (IPRHSO.GT.0) CALL TIMER('TRAHSO',TRASTA,TRATIM)
      ELSE
         CALL GPINQ('MOHSOINT','OPENED',FOPEN)
         IF (.NOT.FOPEN) THEN
            CALL DAINTF(LUMHSO,'MOHSOINT')
         END IF
      END IF
C
C Read distributions
C
      IDIST = 0
      KFREE = 1
      LFREE = LWRK1
 90   CALL NXTHSO(IC,ID,WRK(KH2),NEEDSO,WRK(KWRK1),KFREE,LFREE,IDIST)
      IF (IDIST.LT.0) GOTO 95
      IF (IC.EQ.ID) GOTO 90
      IF (IPRHSO.GT.20) THEN
         WRITE(LUPRI,'(//A,2I5)')' Integral distribution ',IC,ID
         CALL OUTPUT(WRK(KH2),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C   Construct inactive and active spin-orbit Fock matrices
C
      CALL FSOMU(IC,ID,WRK(KH2),WRK(KFC),WRK(KFV),UDV,LORB,
     * WRK(KWRK1-1+KFREE),LFREE)
C     CALL FSOMU(ICI1,IDI1,H2,FCSO,FVSO,UDV,LORB,WRK,LWRK)
C
C
C   Construct spin-orbit Q matrices
C
      IF (NASHT.GT.0 .AND. LORB)
     *   CALL QSOMU(IC,ID,WRK(KQA),WRK(KQB),WRK(KH2),
     *           PV,WRK(KPV12),WRK(KPV21),
     *           WRK(KWRK1-1+KFREE),LFREE)
C
C Add active-active distributions to H2AC(uv,xy)
C     CALL ADH2AC(H2ACXY,H2XY,IUVSYM)
C
      IF (IOBTYP(IC).EQ.JTACT .AND. IOBTYP(ID).EQ.JTACT .AND.
     *      LCON) THEN
         ICSYM = ISMO(IC)
         IDSYM = ISMO(ID)
         ICDSYM = MULD2H(ICSYM,IDSYM)
         KCDSYM = MULD2H(KSYMOP,ICDSYM)
         NCW = ICH(IC)
         NDW = ICH(ID)
         IF (NCW.GT.NDW) THEN
            NCDW = IROW(NCW) + NDW
            KH2XY = KH2AC + (NCDW-1)*NNASHX
            CALL ADH2AC(WRK(KH2XY),WRK(KH2),KCDSYM)
         ELSE
            NCDW = IROW(NDW) + NCW
            KH2XY = KH2AC + (NCDW-1)*NNASHX
            CALL DAXPY(N2ORBX,-D1,WRK(KH2),1,WRK(KWRK1-1+KFREE),1)
            CALL ADH2AC(WRK(KH2XY),WRK(KWRK1-1+KFREE),KCDSYM)
         END IF
      END IF
C
      IF (IPRHSO.GT.20) THEN
         WRITE(LUPRI,'(//A)') ' Partial inactive Fock matrix'
         CALL OUTPUT(WRK(KFC),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      GOTO 90
 95   CONTINUE
C
      IF (IPRHSO.GT.10) THEN
         WRITE(LUPRI,'(//A)') ' Inactive Fock matrix, one-electron part'
         CALL OUTPUT(WRK(KH1),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(//A)') ' Inactive Fock matrix, two-electron part'
         CALL OUTPUT(WRK(KFC),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C Add the one-electron and two electron parts of the inactive Fock matrix
C
      IF (DO1) CALL DAXPY(N2ORBX,D1,WRK(KH1),1,WRK(KFC),1)
C
      IF (IPRHSO.GT.10) THEN
         WRITE(LUPRI,'(//A)') ' Total inactive Fock matrix'
         CALL OUTPUT(WRK(KFC),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
C
         IF (LORB .AND. NASHT.GT.0) THEN
            WRITE(LUPRI,'(//A)') ' Active Fock matrix'
            CALL OUTPUT(WRK(KFV),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
C
            WRITE(LUPRI,'(//A)') ' Spin-orbit QA matrix'
            CALL OUTPUT(WRK(KQA),1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
            WRITE(LUPRI,'(//A)') ' Spin-orbit QB matrix'
            CALL OUTPUT(WRK(KQB),1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
         END IF
C
         IF (LORB .AND. NASHT.GT.0) THEN
            DO 71 K=1,NASHT
               DO 72 L=1,K
                  KL = IROW(K)+L
                  KLOFF = (KL-1)*NNASHX
                  WRITE(LUPRI,'(/A,2I3,A)')
     *         ' Active integral distribution (**', K,L,')'
                  CALL OUTPAK(WRK(KH2AC+KLOFF),NASHT,1,LUPRI)
 72            CONTINUE
 71         CONTINUE
          END IF
      END IF
C
C
C
C Construct orbital property vector
C
      IF (LORB) THEN
         IF (IPRHSO.GT.0) CALL TIMER('START ',ORBSTA,ORBTIM)
         CALL RSPORB(.TRUE.,1,DUMMY1,
     *    WRK(KFC),WRK(KFV),WRK(KQA),WRK(KQB),
     *    UDV,GP,DUMMY2)
C        CALL RSPORB(ONEIND,NSIM,FC,FCX,FVX,QAX,QBX,UDV,
C    *                  EVECS,ZYMAT)
C
C
         IF (IPRHSO.GT.0) CALL TIMER('RSPORB',ORBSTA,ORBTIM)
      END IF
C
C  Compensate for the sign in RSPORB
C
      CALL DSCAL(KZWOPT,DM1,GP(1+KZCONF),1)
      CALL DSCAL(KZWOPT,DM1,GP(1+KZCONF+KZVAR),1)
C
C Construct configuration property vector
C
      IF (LCON) THEN
         IF (IPRHSO.GT.0) CALL TIMER('START ',ORBSTA,ORBTIM)
         CALL HSOSIG(WRK(KFC),WRK(KH2AC),
     *      GP,XINDX,WRK(KWRK1),LWRK1)
C
C        CALL HSOSIG(FC,H2AC, GP,XINDX,WRK,LWRK)
C
C
         IF (IPRHSO.GT.0) CALL TIMER('HSOSIG',ORBSTA,ORBTIM)
      END IF
      IF (IPRHSO.GT.10) THEN
         WRITE(LUPRI,'(//A//A)')
     *      ' Configuration property vector: ',
     *      '               Z part         Y part'
         CALL OUTPUT(GP,1,KZCONF,1,2,KZVAR,2,1,LUPRI)
         WRITE(LUPRI,'(//A//A)')
     *      ' Orbital property vector: ',
     *      '               Z part         Y part'
         CALL OUTPUT(GP(1+KZCONF),1,KZWOPT,1,2,KZVAR,2,1,LUPRI)
      ELSE IF (IPRHSO.GT.0) THEN
         IF (LORB) CALL RSPPRO (GP(1+KZCONF),KZVAR,LUPRI)
         IF (LCON) CALL RSPPRC (GP,KZCONF,KZVAR,LUPRI)
      END IF
C
      IF (IPRHSO.GT.0) CALL TIMER('HSOCTL',HSOSTA,HSOTIM)
C     Test if work(ktest) has been touched:
      IF (WORK(KTEST) .ne. -999.9D0) THEN
         CALL QUIT('Programming error, WORK(KTEST) has been touched!')
      END IF
C
      IF (X2GRAD) THEN
         CALL HEADER('X2GRAD TEST FOR SPIN-ORBIT GRADIENT ELEMENTS',3)
         CALL DZERO(WRK,KWRK2-1)
         CALL SETZY(WRK(KMJWOP))
         CALL HSOAL2 (WORD,CMO,UDV,PV,XINDX,WRK(KMJWOP),
     &                WRK(KGP2),LWRK2)
         WRITE(LUPRI,'(/A)')'               GP1           GP2'
         WRITE(LUPRI,'(/A)')'      ---------------------------------'
         WRITE(LUPRI,'(/A)')' Orbital Z'
         CALL OUTPUT(GP(1+KZCONF),1,KZWOPT,1,2,KZYVAR,2,1,LUPRI)
         WRITE(LUPRI,'(/A)')' Orbital Y'
         CALL OUTPUT(GP(1+KZCONF+KZVAR),1,KZWOPT,1,2,KZYVAR,2,1,LUPRI)
         WRITE(LUPRI,'(/A)')' Configuration Z'
         CALL OUTPUT(GP(1),1,KZCONF,1,2,KZYVAR,2,1,LUPRI)
         WRITE(LUPRI,'(/A)')' Configuration Y'
         CALL OUTPUT(GP(1+KZVAR),1,KZCONF,1,2,KZYVAR,2,1,LUPRI)
         DMAXGP = D0
         DO 77, K = 0,KZYVAR-1
            DMAXGP = MAX(DMAXGP,ABS(GP(1+K)-WRK(KGP2+K)))
77       CONTINUE
         WRITE(LUPRI,'(/A,E20.8)')
     *      'LARGEST DIFFERENCE OF SPIN-ORBIT GRADIENT ELEMENTS' ,
     *      DMAXGP
         CALL QUIT('END OF X2GRAD TEST')
      END IF
      CALL DSCAL(KZYVAR,HSOFAC,GP,1)
C
C Save on file
C
      CALL NEWLAB(WORD,LUHSO,LUERR)
      CALL WRITT(LUHSO,KZYVAR,GP)
      CALL GPCLOSE(LUHSO,'KEEP')
C
 9999 CALL QEXIT('HSOCTL')
      RETURN
      END
C  /* Deck hsoinp */
      SUBROUTINE HSOINP(WORD)
C
#include <implicit.h>
C
#include <priunit.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <infpri.h>
#include <infdim.h>
#include <trhso.h>
#include <infhso.h>
C
      LOGICAL NEWDEF
      PARAMETER ( NTABLE =  8 )
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
      DATA TABLE /'.TESTZY', '.SO1ONL', '.SO2ONL', '.PRINT', 
     *            '.OLDTRA', '.X2MAT ', '.A2MAT ', '.X2GRAD'/
C
      NEWDEF = (WORD .EQ. '*QUADRA')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') GO TO 100
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO I=1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8), I
                  END IF
               END DO
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' KEYWORD "',WORD,
     *            '" NOT RECOGNIZED IN HSOINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT(' ILLEGAL KEYWORD IN HSOINP ')
 1             CONTINUE
                  TESTZY = .TRUE.
               GO TO 100
 2             CONTINUE
                  DOSO2 = .FALSE.
               GO TO 100
 3             CONTINUE
                  DOSO1 = .FALSE.
               GO TO 100
 4             CONTINUE
                  READ(LUCMD,*)IPRHSO
               GO TO 100
 5             CONTINUE
                  OLDTRA = .TRUE.
               GO TO 100
 6             CONTINUE
                  X2MAT = .TRUE.
               GO TO 100
 7             CONTINUE
                  A2MAT = .TRUE.
               GO TO 100
 8             CONTINUE
                  X2GRAD = .TRUE.
               GO TO 100
            ELSE IF (PROMPT.EQ.'*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' PROMPT "',WORD,
     *            '" NOT RECOGNIZED IN RSPINP.'
               CALL QUIT(' ILLEGAL PROMPT IN HSOINP ')
            END IF
         GO TO 100
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('CHANGES OF DEFAULTS FOR HSOINP:',0)
         IF (TESTZY) WRITE(LUPRI,'(/A/A,L1)')
     *      ' Both parts of configuration property vector explicitely',
     *      ' TESTZY = ',TESTZY
         IF (IPRHSO.NE.0) WRITE(LUPRI,'(/A/A,I5)')
     *      ' Print level in spin-orbit property vector calculation',
     *      ' IPRHSO = ',IPRHSO
         IF (OLDTRA)  WRITE(LUPRI,'(/A/A,L1)')
     *      ' Use existing two-electron spin-orbit integral file',
     *      ' OLDTRA = ',OLDTRA
         IF (.NOT.DOSO1) WRITE(LUPRI,'(/A)')
     *      'Skip one-particle part in HSOCTL'
         IF (.NOT.DOSO2) WRITE(LUPRI,'(/A)')
     *      'Skip two-particle part in HSOCTL'
         IF (X2MAT) WRITE(LUPRI,'(/A)')
     *      'X2MAT: Calculate full X2 matrix (quadratic response)'
         IF (A2MAT) WRITE(LUPRI,'(/A)')
     *      'A2MAT: not implemented'
      END IF
C
C *** END OF HSOINP
C
      RETURN
      END
C  /* Deck fsomu */
      SUBROUTINE FSOMU(ICI1,IDI1,H2,FCSO,FVSO,UDV,LORB,WRK,LWRK)
C
C Olav Vahtras
C Apr 11, 1990
C
C CALCULATE ALL CONTRIBUTIONS TO INACTIVE AND ACTIVE SPIN-ORBIT
C FOCK MATRICES FROM MULLIKEN DISTRIBUTIONS
C
C  FCSO(P,Q) = SUM (K) 2*(KK|P^Q) - SUM (K) 3*(PK|K^Q)
C                                 - SUM (K) 3*(KQ|P^K)
C
C  FVSO(P,Q) = SUM (X,Y) (XY|P^Q) D(XY) - SUM (K) 3/2*(PX|Y^Q) D(XY)
C                                     - SUM (K) 3/2*(XQ|P^Y) D(XY)
C
C
#include <implicit.h>
C
C Used from common blocks:
C
C   INFORB : NORBT,NISHT,NASHT,MULD2H
C   INFIND : ISMO,IOBTYP,ICH,ISMO,IASH,IORB,NISH,NASH,NOCC,NORB
C   INFHSO :
C   WRKRSP : KSYMOP
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infpri.h>
#include <orbtypdef.h>
#include <trhso.h>
#include <infhso.h>
C
      DIMENSION H2(NORBT,NORBT)
      DIMENSION FCSO(NORBT,NORBT),FVSO(NORBT,NORBT)
      DIMENSION UDV(NASHT,NASHT),WRK(*)
      LOGICAL LORB
C
      PARAMETER (D1=1.0D0, D1P5=1.5D0)
C
      CALL QENTER('FSOMU')
C
C
C  Local print level
C
      IPRFSO = 20
C
C     Order (C,D) index such that C .ge. D
C     in inactive-active-secondary order (using ISW)
C
      IF (ISW(ICI1) .GE. ISW(IDI1)) THEN
         ICI = ICI1
         IDI = IDI1
         DISFAC=D1
      ELSE
         ICI = IDI1
         IDI = ICI1
         DISFAC=-D1
      END IF
C
C     Find distribution type ITYPCD
C
      ITYPC  = IOBTYP(ICI)
      ITYPD  = IOBTYP(IDI)
      ITYPCD = IDBTYP(ITYPC,ITYPD)
C
C     We only need secondary-occupied distributions
C
      IF (ITYPCD .EQ. JTSESE) GO TO 9999
C
      IF (ITYPC .EQ. JTACT) NCIW = ICH(ICI)
      IF (ITYPD .EQ. JTACT) NDIW = ICH(IDI)
C
      ICSYM = ISMO(ICI)
      IDSYM = ISMO(IDI)
      ICDSYM = MULD2H(ICSYM,IDSYM)
      KCDSYM = MULD2H(KSYMOP,ICDSYM)
C
      IF (IPRHSO.GT.IPRFSO) THEN
         WRITE(LUPRI,'(/A)')' <<<<<< Output from FSOMU >>>>>>'
         WRITE(LUPRI,'(/A,2I5,5X,A,2X,A)')' Distribution CD',ICI1,IDI1,
     *                          COBTYP(ITYPC),COBTYP(ITYPD)
         WRITE(LUPRI,'(A,2I5)')' Reordered      ',ICI,IDI
         WRITE(LUPRI,'(A,2I5)')' Symmetry       ',ICSYM,IDSYM
      ENDIF
C
C
C Inactive Fock matrix
C
C Direct terms
C
C  FCSO(P,Q) = FCSO(P,Q) + SUM(K) 2*(KK|P^Q)
C
C here:
C        SUM(K) 2*(KK|C^D) -> FCSO(C,D)
C      - SUM(K) 2*(KK|C^D) -> FCSO(D,C)
C
         IF ( KSYMOP.EQ.ICDSYM ) THEN
            DO 10 I=1,NISHT
               IX=ISX(I)
               WRK(I)=H2(IX,IX)
10          CONTINUE
               IF (IPRHSO.GT.IPRFSO) THEN
                  WRITE(LUPRI,'(/A)')' Inactive direct terms'
                  WRITE(LUPRI,'(A)')' (KK|CD) diagonal'
                  CALL OUTPUT(WRK(1),1,NISHT,1,1,NISHT,1,1,LUPRI)
               END IF
               FAC=2*DISFAC*DSUM(NISHT,WRK(1),1)
               FCSO(ICI,IDI)=FCSO(ICI,IDI) + FAC
               FCSO(IDI,ICI)=FCSO(IDI,ICI) - FAC
         ENDIF
C
C Exchange terms: rearranged with D position inactive
C
C  FCSO(P,Q) = FCSO(P,Q) + SUM(K) 3*(PK|Q^K) - 3*(QK|P^K)
C
C
         IF (IPRHSO.GT.IPRFSO) THEN
            WRITE(LUPRI,'(/A)')' Loop over symmetries'
         END IF
         DO 200 ISYM = 1,NSYM
            IPSYM = ISYM
            IOFFP = IORB(IPSYM)
            NASHP = NASH(IPSYM)
            NOCCP = NOCC(IPSYM)
            NORBP = NORB(IPSYM)
            IASHP=IASH(IPSYM)
            IOFFPA=IOFFP+NISH(IPSYM)
            IF (IPRHSO.GT.IPRFSO) THEN
               WRITE(LUPRI,'(/A,I5)')' IPSYM',IPSYM
            ENDIF
            IF ( (NORBP.EQ.0) ) GO TO 200
C
            ICPSYM = MULD2H(ICSYM,IPSYM)
            IDPSYM = MULD2H(IDSYM,IPSYM)
C
C For the case D inactive
C
            IF ((ITYPCD.EQ.JTININ).OR.
     *          (ITYPCD.EQ.JTACIN).OR.
     *          (ITYPCD.EQ.JTSEIN)) THEN
C
C here:
C         + 3*(PD|C^D) -> FCSO(P,C)
C         - 3*(PD|C^D) -> FCSO(C,P)
C
C such that (PC) is at most secondary-active
C
               IF (IPRHSO.GT.IPRFSO) THEN
                  IF (ICPSYM.EQ.KSYMOP .AND.
     *                (ITYPCD.EQ.JTSEIN .AND. NOCCP.GT.0
     *                        .OR.
     *                      ITYPCD.NE.JTSEIN )
     *                        .OR.
     *                IDPSYM.EQ.KSYMOP .AND. ITYPCD.EQ.JTININ) THEN
                     WRITE(LUPRI,'(/A)')' Inactive exchange terms'
                  END IF
               END IF
               FAC = 3*DISFAC
               IF ( ICPSYM.EQ.KSYMOP ) THEN
                  IF (ITYPCD.EQ.JTSEIN) THEN
                     NDIMP = NOCCP
                  ELSE
                     NDIMP = NORBP
                  ENDIF
                  IF (NDIMP.GT.0) THEN
                     CALL DAXPY(NDIMP,FAC,H2(IOFFP+1,IDI),1,
     *                                 FCSO(IOFFP+1,ICI),1)
                     CALL DAXPY(NDIMP,-FAC,H2(IOFFP+1,IDI),1,
     *                                 FCSO(ICI,IOFFP+1),NORBT)
                     IF (IPRHSO.GT.IPRFSO) THEN
                       WRITE(LUPRI,'(A,I3,A,I3,I5)')' PC contribution ',
     *                    IOFFP+1,':',IOFFP+NDIMP,ICI
                       CALL OUTPUT(H2(IOFFP+1,IDI),1,NDIMP,1,1,
     *                                             NORBT,NORBT,1,LUPRI)
                     ENDIF
                  END IF
               END IF
C
C if both C and D are inactive we also have
C
C        - 3*(PC|C^D) -> FCSO(P,D)
C        + 3*(PC|C^D) -> FCSO(D,P)
C
               IF ( IDPSYM.EQ.KSYMOP .AND. ITYPCD.EQ.JTININ ) THEN
                  IF (IPRHSO.GT.IPRFSO) THEN
                     WRITE(LUPRI,'(A,I3,A,I3,I5)')' PD contribution ',
     *                  IOFFP+1,':',IOFFP+NDIMP,IDI
                     CALL OUTPUT(H2(IOFFP+1,ICI),1,NORBP,1,1,
     *                                         NORBT,NORBT,1,LUPRI)
                  ENDIF
                  CALL DAXPY(NORBP,-FAC,H2(IOFFP+1,ICI),1,
     *                                 FCSO(IOFFP+1,IDI),1)
                  CALL DAXPY(NORBP,FAC,H2(IOFFP+1,ICI),1,
     *                                 FCSO(IDI,IOFFP+1),NORBT)
               ENDIF
            ENDIF
C
         IF (LORB) THEN
C
C Active Fock matrix
C
C Direct terms:
C
C  FVSO(P,Q) = FVSO(P,Q)+ SUM(X,Y) (XY|P^Q)*D(XY)
C
C
C here:
C        + SUM(X,Y) (XY|C^D)*D(XY) -> FVSO(C,D)
C        - SUM(X,Y) (XY|C^D)*D(XY) -> FVSO(D,C)
C
C where the sum is taken over diagonal symmetry blocks (X,Y)
C
            FAC = DISFAC
            IXSYM = ISYM
            IASHX = IASH(IXSYM)
            NASHX = NASH(IXSYM)
            IOFFXA = IORB(IXSYM) + NISH(IXSYM)
            IF (KSYMOP.EQ.ICDSYM) THEN
               IF (IPRHSO.GT.IPRFSO) THEN
                  WRITE(LUPRI,'(/A)')' Active direct terms'
               END IF
               DO 20 IX=1,NASHX
                  WRK(IX) = DDOT(NASHX,H2(IOFFXA+1,IOFFXA+IX),1,
     *                           UDV(IASHX+1,IASHX+IX),1)
 20            CONTINUE
               XYSUM = DSUM(NASHX,WRK(1),1)
               FVSO(ICI,IDI) = FVSO(ICI,IDI) + FAC*XYSUM
               FVSO(IDI,ICI) = FVSO(IDI,ICI) - FAC*XYSUM
            END IF
C
C Exchange terms: rearranged with C position active
C
C  FVSO(P,Q) = FVSO(P,Q) - SUM(X,Y) 3/2*(PX|Y^Q)*D(X,Y)
C                        + SUM(X,Y) 3/2*(QX|Y^P)*D(X,Y)
C
C
C here:
C       - SUM(X) 3/2*(PX|C^D)*D(X,C) -> FVSO(P,D)
C       + SUM(X) 3/2*(PX|C^D)*D(X,C) -> FVSO(D,P)
C
C for active-active and active-inactive distributions
C
            FAC=D1P5*DISFAC
            IXSYM = MULD2H(IPSYM,KCDSYM)
            IF (IPRHSO.GT.IPRFSO) THEN
               IF (IXSYM.EQ.ICSYM .AND.
     *                   (ITYPCD.EQ.JTACIN .OR. ITYPCD.EQ.JTACAC)
     *                          .OR.
     *             IXSYM.EQ.IDSYM .AND.
     *                   (ITYPCD.EQ.JTACAC .OR.
     *                          NOCCP.GT.0 .AND. ITYPCD.EQ.JTSEAC)) THEN
                  WRITE(LUPRI,'(/A)')' Active exchange terms'
               END IF
            END IF
            IF (IXSYM.EQ.ICSYM) THEN
               IOFFXA = IORB(IXSYM) + NISH(IXSYM)
               NASHX = NASH(IXSYM)
               IASHX = IASH(IXSYM)
               IF (ITYPCD.EQ.JTACIN .OR. ITYPCD.EQ.JTACAC) THEN
                  CALL DGEMM('N','N',NORBP,1,NASHX,1.D0,
     &                       H2(IOFFP+1,IOFFXA+1),NORBT,
     &                       UDV(IASHX+1,NCIW),NASHT,0.D0,
     &                       WRK(1),NORBP)
                  CALL DAXPY(NORBP,-FAC,WRK(1),1,
     *                              FVSO(IOFFP+1,IDI),1)
                  CALL DAXPY(NORBP,FAC,WRK(1),1,
     *                             FVSO(IDI,IOFFP+1),NORBT)
                  IF (IPRHSO.GT.IPRFSO) THEN
                     WRITE(LUPRI,'(A,I3,A,I3,I5)')' PD contribution ',
     *                  IOFFP+1,':',IOFFP+NORBP,IDI
                     CALL OUTPUT(WRK(1),1,NORBP,1,1,NORBP,1,1,LUPRI)
                  ENDIF
               END IF
            END IF
C
C Exchange terms: rearranged with D position active
C
C  FVSO(P,Q) = FVSO(P,Q) + SUM(X,Y) 3/2*(PX|Q^Y)*D(X,Y)
C                        - SUM(X,Y) 3/2*(QX|P^Y)*D(X,Y)
C
C
C here:
C        SUM(X) 3/2*(PX|C^D)*D(X,D) -> FVSO(P,C)
C      - SUM(X) 3/2*(PX|C^D)*D(X,D) -> FVSO(C,P)
C
C such that (CP) is at most secondary-active
C
C for active-active and secondary-active distributions
C
            IF (IXSYM.EQ.IDSYM) THEN
               IOFFXA = IORB(IXSYM) + NISH(IXSYM)
               NASHX = NASH(IXSYM)
               IASHX = IASH(IXSYM)
               IF (ITYPCD.EQ.JTACAC .OR. ITYPCD.EQ.JTSEAC) THEN
                  IF (ITYPCD.EQ.JTSEAC) THEN
                     NDIMP = NOCCP
                  ELSE
                     NDIMP = NORBP
                  END IF
                  IF (NDIMP .GT. 0) THEN
                     CALL DGEMM('N','N',NDIMP,1,NASHX,1.D0,
     &                    H2(IOFFP+1,IOFFXA+1),NORBT,
     &                    UDV(IASHX+1,NDIW),NASHT,0.D0,
     &                    WRK(1),NDIMP)
                     CALL DAXPY(NDIMP,FAC,WRK(1),1,
     *                    FVSO(IOFFP+1,ICI),1)
                     CALL DAXPY(NDIMP,-FAC,WRK(1),1,
     *                    FVSO(ICI,IOFFP+1),NORBT)
                     IF (IPRHSO.GT.IPRFSO) THEN
                        WRITE(LUPRI,'(A,I3,A,I3,I5)')
     *                       ' PC contribution ',
     *                       IOFFP+1,':',IOFFP+NDIMP,ICI
                        CALL OUTPUT(WRK(1),NDIMP,1,1,1,NDIMP,1,1,LUPRI)
                     ENDIF
                  END IF
               END IF
            END IF
C
         END IF
C        (RSPCI)
200      CONTINUE
9999  CALL QEXIT('FSOMU')
      RETURN
      END
C  /* Deck qsomu */
      SUBROUTINE QSOMU(ICI,IDI,QASO,QBSO,
     *                  H2,PVX,PV12,PV21,
     *                  WRK,LWRK)
C
C Olav Vahtras
C Apr 11, 1990
C
C Purpose:
C  Calculate all contributions to QA and QB spin-orbit matrices
C  from Mulliken (**|C^D) integral distributions
C
C  In general:
C
C  QBSO(P,Q) = SUM(X,Y,W) (PW|X^Y)*( 2*PV(++)(Q,W,X,Y) + PV(--)(Q,W,X,Y) )
C            + SUM(X,Y,W) (XY|P^W)*( 2*PV(--)(X,Y,Q,W) + PV(++)(X,Y,Q,W) )
C
C  QASO(P,Q) = SUM(X,Y,W) (WP|X^Y)*( 2*PV(++)(W,Q,X,Y) + PV(--)(W,Q,X,Y) )
C            + SUM(X,Y,W) (XY|W^P)*( 2*PV(--)(X,Y,W,Q) + PV(++)(X,Y,W,Q) )
C
C   where PV(++)(P,Q,R,S) = <0| e(+,+)(pqrs) |0>
C    and  PV(--)(P,Q,R,S) = <0| e(-,-)(pqrs) |0>
C
#include <implicit.h>
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infpri.h>
#include <orbtypdef.h>
#include <trhso.h>
#include <infhso.h>
C
C
      DIMENSION QASO(NORBT,NASHT),QBSO(NORBT,NASHT)
      DIMENSION PV12(NASHT,NASHT),PV21(NASHT,NASHT)
      DIMENSION H2(NORBT,NORBT),PVX(*),WRK(*)
C
      PARAMETER(D2=2.0D0)
C
      CALL QENTER('QSOMU')
C
C  Local print level
C
      IPRQSO = 20
C
      IF (LWRK.LT.N2ASHX) CALL ERRWRK('QSOMU',N2ASHX,LWRK)
C
C
C Symmetry to type order
C
      ICIW = ISW(ICI)
      IDIW = ISW(IDI)
C
C Orbital type, at least one has to be active to contribute
C
      ITYPC = IOBTYP(ICI)
      ITYPD = IOBTYP(IDI)
      IF (ITYPC.NE.JTACT .AND. ITYPD.NE.JTACT) GO TO 9999
C
C Distribution type and symmetry
C
      ICSYM = ISMO(ICI)
      IDSYM = ISMO(IDI)
      ICDSYM = MULD2H(ICSYM,IDSYM)
      KCDSYM = MULD2H(KSYMOP,ICDSYM)
      ITYPCD=IDBTYP(ITYPC,ITYPD)
C
C     Order within actives for the case that C or D are active
C
      IF (ITYPC .EQ. JTACT) NCIW = ICIW - NISHT
      IF (ITYPD .EQ. JTACT) NDIW = IDIW - NISHT
C
      IF (IPRHSO.GT.IPRQSO) THEN
         WRITE(LUPRI,'(/A)') ' <<<<<< Output from QSOMU >>>>>>'
         WRITE(LUPRI,'(/A,2I5,5X,2A)')' Distribution CD',ICI,IDI,
     *                          COBTYP(ITYPC),COBTYP(ITYPD)
         WRITE(LUPRI,'(A,2I5)')       ' Symmetry       ',ICSYM,IDSYM
      ENDIF
C
C
            IPP = 1
            IMM = N2ASHX*N2ASHX+1
C
C Both C and D are active:
C
      IF (ITYPCD.EQ.JTACAC) THEN
C
C Get (C,D) density distributions in the form  2PV(++)+PV(--)
C
            NCDOFF = (NCIW-1 + (NDIW-1)*NASHT)*N2ASHX
            CALL DCOPY(N2ASHX,PVX(NCDOFF+IMM),1,PV12,1)
            CALL DAXPY(N2ASHX,D2,PVX(NCDOFF+IPP),1,PV12,1)
            IF (IPRHSO.GT.IPRQSO) THEN
               WRITE(LUPRI,'(/A)')' Active diagonal terms'
               WRITE(LUPRI,'(/A,2I5)')
     *        ' Total CD density distribution' ,NCIW,NDIW
               CALL OUTPUT(PV12,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
            END IF
C
C Get (D,C) density distributions
C
            NDCOFF = (NDIW-1 + (NCIW-1)*NASHT)*N2ASHX
            CALL DCOPY(N2ASHX,PVX(NDCOFF+IMM),1,PV21,1)
            CALL DAXPY(N2ASHX,D2,PVX(NDCOFF+IPP),1,PV21,1)
            IF (IPRHSO.GT.IPRQSO) THEN
               WRITE(LUPRI,'(/A,2I5)')
     *       ' Total DC density distribution' ,NDIW,NCIW
               CALL OUTPUT(PV21,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
            END IF
C
C Add contibutions to QASO and QBSO from (C,D) and (D,C) density
C distributions
C
C
               IF (IPRHSO.GT.IPRQSO) THEN
                  WRITE(LUPRI,'(/A)') ' Loop over symmetry blocks '
               END IF
               DO 100 IPSYM = 1,NSYM
                  IWSYM = MULD2H(IPSYM,KCDSYM)
                  IQSYM = MULD2H(IWSYM,ICDSYM)
                  NORBP = NORB(IPSYM)
                  NASHP = NASH(IPSYM)
                  NASHQ = NASH(IQSYM)
                  NASHW = NASH(IWSYM)
                  IF (IPRHSO.GT.IPRQSO) THEN
                     WRITE(LUPRI,'(/A,I5)')' IPSYM',IPSYM
                  END IF
                  IF ( NORBP.GT.0 .AND. NASHW.GT.0 .AND. NASHQ.GT.0 )
     *                                                           THEN
                     IOFFP = IORB(IPSYM)
                     IOFFW = IORB(IWSYM)
                     IOFFPA = IOFFP + NISH(IPSYM)
                     IOFFWA = IOFFW + NISH(IWSYM)
                     IASHP = IASH(IPSYM)
                     IASHQ = IASH(IQSYM)
                     IASHW = IASH(IWSYM)
                     IF (IPRHSO.GT.IPRQSO) THEN
                        WRITE(LUPRI,'(/A,2I5)')
     *                  ' Integral symmetry block PW',IPSYM,IWSYM
                        CALL OUTPUT(H2,IOFFP+1,IOFFP+NORBP,
     *                                 IOFFWA+1,IOFFWA+NASHW,
     *                                 NORBT,NORBT,1,LUPRI)
                        WRITE(LUPRI,'(/A,2I5)')
     *                  ' Density symmetry block QW',IQSYM,IWSYM
                        CALL OUTPUT(PV12,IASHQ+1,IASHQ+NASHQ,
     *                                   IASHW+1,IASHW+NASHW,
     *                                   NASHT,NASHT,1,LUPRI)
                     END IF
C
C QBSO(P,Q) += SUM(W) (PW|C^D)*(2*PV(++)(Q,W,C,D) + PV(--)(Q,W,C,D))
C
C
               CALL DGEMM('N','T',NORBP,NASHQ,NASHW,1.D0,
     &                    H2(IOFFP+1,IOFFWA+1),NORBT,
     &                    PV12(IASHQ+1,IASHW+1),NASHT,1.D0,
     &                    QBSO(IOFFP+1,IASHQ+1),NORBT)
C
C QBSO(P,Q) += SUM(W) (PW|D^C)*(2*PV(++)(Q,W,D,C) + PV(--)(Q,W,D,C))
C        (  -     ...     C^D   )
C
               CALL DGEMM('N','T',NORBP,NASHQ,NASHW,-1.D0,
     &                    H2(IOFFP+1,IOFFWA+1),NORBT,
     &                    PV21(IASHQ+1,IASHW+1),NASHT,1.D0,
     &                    QBSO(IOFFP+1,IASHQ+1),NORBT)
C
C QASO(P,Q) += SUM(W) (WP|C^D)*(2*PV(++)(W,Q,C,D) + PV(--)(W,Q,C,D))
C
               CALL DGEMM('T','N',NORBP,NASHQ,NASHW,1.D0,
     &                    H2(IOFFWA+1,IOFFP+1),NORBT,
     &                    PV12(IASHW+1,IASHQ+1),NASHT,1.D0,
     &                    QASO(IOFFP+1,IASHQ+1),NORBT)
C
C QASO(P,Q) += SUM(W) (WP|D^C)*(2*PV(++)(W,Q,D,C) + PV(--)(W,Q,D,C))
C        (  -     ...     C^D   )
C
               CALL DGEMM('T','N',NORBP,NASHQ,NASHW,-1.D0,
     &                    H2(IOFFWA+1,IOFFP+1),NORBT,
     &                    PV21(IASHW+1,IASHQ+1),NASHT,1.D0,
     &                    QASO(IOFFP+1,IASHQ+1),NORBT)
C
                  END IF
 100           CONTINUE
      END IF
C
C Following integral distributions contribute if either
C C or D is active.
C
C
C QBSO(C,Q) += SUM(X,Y) (XY|C^D)*(2*PV(--)(X,Y,Q,D)+PV(++)(X,Y,Q,D))
C
      IF (ITYPD.EQ.JTACT) THEN
         IQSYM = MULD2H(KSYMOP,ICSYM)
         NASHQ = NASH(IQSYM)
         IASHQ = IASH(IQSYM)
         IF (NASHQ.GT.0) THEN
            IF (IPRHSO.GT.IPRQSO) THEN
               WRITE(LUPRI,'(/A,I5)')
     *         ' Loop over actives in symmetry',IQSYM
            END IF
            DO 200 IQ=1,NASHQ
C
C For each Q read a new (Q,D) density distribution in the form
C 2PV(--) + PV(++)
C
               IF (IPRHSO.GT.IPRQSO) THEN
                  WRITE(LUPRI,'(/A,2I5)')' IQ',IASHQ+IQ
               END IF
               NQDOFF = (IASHQ+IQ-1 + (NDIW-1)*NASHT)*N2ASHX
               CALL DCOPY(N2ASHX,PVX(NQDOFF+IPP),1,PV12,1)
               CALL DAXPY(N2ASHX,D2,PVX(NQDOFF+IMM),1,PV12,1)
               IF (IPRHSO.GT.IPRQSO) THEN
                  WRITE(LUPRI,'(/A,2I5)')
     *           ' Total QD density distribution', IASHQ+IQ,NDIW
                  CALL OUTPUT(PV12,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
               END IF
C
C
                  DO 20 IYSYM = 1,NSYM
                     NASHY = NASH(IYSYM)
                     IASHY = IASH(IYSYM)
                     IOFFYA = IORB(IYSYM) + NISH(IYSYM)
                     IXSYM = MULD2H(IYSYM,KCDSYM)
                     NASHX = NASH(IXSYM)
                     IASHX = IASH(IXSYM)
                     IOFFXA = IORB(IXSYM) + NISH(IXSYM)
                     IF (NASHX.GT.0 .AND .NASHY .GT.0) THEN
C
C QBSO(C,Q) += SUM(X,Y) (XY|C^D)*(2*PV(--)(X,Y,Q,D)+PV(++)(X,Y,Q,D))
C
C
C QASO(C,Q) += SUM(X,Y) (XY|D^C)*(2*PV(--)(X,Y,D,Q)+PV(++)(X,Y,D,Q))
C        (  -       ...     C^D   )           (Q,D)           (Q,D)
C
                     IF (IPRHSO.GT.IPRQSO) THEN
                        WRITE(LUPRI,'(/A,2I5)')
     *                  ' Density symmetry block XY',IXSYM,IYSYM
                        CALL OUTPUT(PV12,IASHX+1,IASHX+NASHX,
     *                                   IASHY+1,IASHY+NASHY,
     *                                   NASHT,NASHT,1,LUPRI)
                     END IF
                        DO 21 IY=1,NASHY
                          WRK(IY) = DDOT(NASHX,H2(IOFFXA+1,IOFFYA+IY),1,
     *                                         PV12(IASHX+1,IASHY+IY),1)
 21                     CONTINUE
                        XYSUM = DSUM(NASHY,WRK,1)
                        QBSO(ICI,IASHQ+IQ)
     *                 = QBSO(ICI,IASHQ+IQ) + XYSUM
                        QASO(ICI,IASHQ+IQ)
     *                 = QASO(ICI,IASHQ+IQ) - XYSUM
                     END IF
 20               CONTINUE
 200        CONTINUE
         END IF
      ENDIF
C
C QBSO(D,Q) += SUM(X,Y) (XY|D^C)*(2*PV(--)(X,Y,Q,C)+PV(++)(X,Y,Q,C))
C
      IF (ITYPC.EQ.JTACT) THEN
         IQSYM = MULD2H(KSYMOP,IDSYM)
         NASHQ = NASH(IQSYM)
         IASHQ = IASH(IQSYM)
         IF (NASHQ.GT.0) THEN
            IF (IPRHSO.GT.IPRQSO) THEN
               WRITE(LUPRI,'(/A,I5)')
     *         ' Loop over actives in symmetry',IQSYM
            END IF
            DO 300 IQ = 1,NASHQ
               IF (IPRHSO.GT.IPRQSO) THEN
                  WRITE(LUPRI,'(/A,2I5)')' IQ',IASHQ+IQ
               END IF
C
C For each Q read a new (Q,C) density distribution in the form
C 2PV(--) + PV(++)
C
               NQCOFF = (IASHQ+IQ-1 + (NCIW-1)*NASHT)*N2ASHX
               CALL DCOPY(N2ASHX,PVX(NQCOFF+IPP),1,PV21,1)
               CALL DAXPY(N2ASHX,D2,PVX(NQCOFF+IMM),1,PV21,1)
               IF (IPRHSO.GT.IPRQSO) THEN
                  WRITE(LUPRI,'(/A,2I5)')
     *           ' Total QC density distribution', IASHQ+IQ,NCIW
                  CALL OUTPUT(PV21,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
               END IF
C
C
                  DO 30 IYSYM = 1,NSYM
                     NASHY = NASH(IYSYM)
                     IASHY = IASH(IYSYM)
                     IOFFYA = IORB(IYSYM) + NISH(IYSYM)
                     IXSYM = MULD2H(IYSYM,KCDSYM)
                     NASHX = NASH(IXSYM)
                     IASHX = IASH(IXSYM)
                     IOFFXA = IORB(IXSYM) + NISH(IXSYM)
                     IF (NASHX.GT.0 .AND. NASHY.GT.0) THEN
C
C QBSO(D,Q) += SUM(X,Y) (XY|D^C)*(2*PV(--)(X,Y,Q,C)+PV(++)(X,Y,Q,C))
C        (  -       ...     C^D   )
C
C
C QASO(D,Q) += SUM(X,Y) (XY|C^D)*(2*PV(--)(X,Y,C,Q)+PV(++)(X,Y,C,Q))
C                                             (Q,C)           (Q,C)
                     IF (IPRHSO.GT.IPRQSO) THEN
                        WRITE(LUPRI,'(/A,2I5)')
     *                  ' Density symmetry block XY',IXSYM,IYSYM
                        CALL OUTPUT(PV21,IASHX+1,IASHX+NASHX,
     *                                   IASHY+1,IASHY+NASHY,
     *                                   NASHT,NASHT,1,LUPRI)
                     END IF
                        DO 31 IY=1,NASHY
                          WRK(IY) = DDOT(NASHX,H2(IOFFXA+1,IOFFYA+IY),1,
     *                                         PV21(IASHX+1,IASHY+IY),1)
 31                     CONTINUE
                        XYSUM = DSUM(NASHY,WRK,1)
                        QBSO(IDI,IASHQ+IQ)
     *                 = QBSO(IDI,IASHQ+IQ) - XYSUM
                        QASO(IDI,IASHQ+IQ)
     *                 = QASO(IDI,IASHQ+IQ) + XYSUM
                     END IF
 30               CONTINUE
C
 300        CONTINUE
         END IF
      END IF
C
 9999 CALL QEXIT('QSOMU')
      RETURN
      END
C  /* Deck hsosig */
      SUBROUTINE HSOSIG(FC,H2AC,
     *                  GP,XINDX,WRK,LWRK)
C
C  Calculate the configuration part of the spin-orbit property vector
C
C                ( <J,HSO,0> )
C      GP(J) =   (           )
C                (-<0,HSO,J> )
C
C
C H2AC contain active two-electron spin-orbit integrals in doubly
C triangularly packed form
C
C Olav Vahtras
C Sep 14, 1990
#include <implicit.h>
C
      DIMENSION FC(NORBT,NORBT)
      DIMENSION H2AC(NNASHX,NNASHX)
      DIMENSION GP(KZYVAR),XINDX(*),WRK(*)
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infopt.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
#include <infdim.h>
#include <cbgetdis.h>
#include <trhso.h>
#include <infhso.h>
C
C
      PARAMETER ( DM1 = -1.0D0 )
C
      CALL QENTER('HSOSIG')
C
C  Local print level
C
      IPRLOC = 20
      IF (IPRHSO.GT.IPRLOC) THEN
         WRITE(LUPRI,'(/A)') ' <<<<<< Output from HSOSIG >>>>>>'
      END IF
C
C ALLOCATE WORK SPACE
C
      IF (IREFSY .EQ. KSYMST) THEN
         NDREF = KZCONF
C        ... if KZCONF .ne. NCREF, we need CREF in determinants
      ELSE
         NDREF = NCREF
      END IF
      KFCAC = 1
      KREFCO = KFCAC + N2ASHX
      KHSQSQ  = KREFCO + NDREF
      KFREE  = KHSQSQ + N2ASHX*N2ASHX
      IF (TESTZY) THEN
         KHSQTR = KFREE
         KFREE  = KHSQTR + N2ASHX*N2ASHX
      END IF
      LFREE  = LWRK   - KFREE+1
      IADINT = -1
C
      IF (LWRK.LT.0) CALL ERRWRK('HSOSIG',KFREE-1,LWRK)
      CALL DZERO(WRK,KFREE-1)
C
C Obtain square packed combined two-electron integrals
C
      CALL H2ACSO(H2AC,WRK(KHSQSQ))
C
      ISPIN1 = 0
      ISPIN2 = 1
C
      CALL GETREF(WRK(KREFCO),NDREF)
C
C CREATE Z PART  <0,HSO,J> OF SPIN-ORBIT GP VECTOR
C
C        Get FCAC matrix for Z sigma vector
C        (note: CISIGD requires UFCAC(I,J) = FCXAC(J,I))
C
         DO 220 IW = 1,NASHT
            IX = ISX(NISHT+IW)
            DO 230 JW = 1,NASHT
               JX = ISX(NISHT+JW)
               IJW = (IW-1) * NASHT + JW
               WRK(KFCAC-1+IJW) = FC(IX,JX)
 230        CONTINUE
 220     CONTINUE
C
         IF (IPRHSO.GT.IPRLOC) THEN
            WRITE(LUPRI,'(/A)')' Inactive Fock matrix'
            CALL OUTPUT(FC,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            WRITE(LUPRI,'(/A)')' Active part of inactive Fock matrix'
            CALL OUTPUT(WRK(KFCAC),1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         END IF
         DISTYP = 6
C
         KSYMST = MULD2H(IREFSY,KSYMSO)
         CALL CISIGD(IREFSY, KSYMST, NDREF, KZCONF, WRK(KREFCO),GP,
     *               WRK(KFCAC), WRK(KHSQSQ), .FALSE., .FALSE.,
     *               XINDX, ISPIN1, ISPIN2, WRK(KFREE), LFREE)
C        CALL RSPSIG(ICSYM,IHCSYM,NCDET,NHCDET,C,HC,UFCAC,H2AC,IFLAG,
C    *               NOH2,WORK,KFREE,LFREE)
C
C        IF ((NDREF.NE.NCDET).OR.(KZCONF.NE.NHCDET)) THEN
C           WRITE(LUPRI,'(/2(A,I5),/3(A,I5))')
C    *      ' NUMBER OF REFERENCE DETERMINANTS ,NDREF:',NDREF,
C    *      ' CALCULATED NUMBER ,NCDET:',NCDET,
C    *      ' NUMBER OF DETERMINANTS FOR SYMMETRY',KSYMOP,
C    *      '  IS:',KZCONF,'  CALCULATED NUMBER,NHCDET:',NHCDET
C        CALL QUIT(' H2XSIG,INCORRECT CALCULATION OF DETERMINANTS')
C        END IF
         IF (IREFSY .EQ. KSYMST) THEN
C           ... remove CREF component of E[2] vector
            T1 = DDOT(KZCONF,WRK(KREFCO),1,GP(1),1)
            CALL DAXPY(KZCONF,(-T1),WRK(KREFCO),1,GP(1),1)
         END IF
         IF (IPRHSO.GT.IPRLOC) THEN
            WRITE(LUPRI,'(/A)')
     *        ' CSF part of spin-orbit property vector: Z part'
            CALL OUTPUT(GP,1,KZCONF,1,1,KZCONF,1,1,LUPRI)
         END IF
C
C CREATE Y PART  <0,HSO,J> OF SPIN-ORBIT GP VECTOR
C
C Normal procedure: copy the Z part
C Test procedure: do it explicitely
C
C        Get transposed FCXAC matrix for Y sigma vector
C        (note: CISIGD requires UFCAC(I,J) = FCXAC(J,I))
C
         IF (TESTZY) THEN
C
C           ISPIN1 = 0
C           ISPIN2 = 1
C
            DO 120 IW = 1,NASHT
               IX = ISX(NISHT+IW)
               DO 130 JW = 1,NASHT
                  JX = ISX(NISHT+JW)
                  IJW = (IW-1) * NASHT + JW
                  WRK(KFCAC-1+IJW) = FC(JX,IX)
 130           CONTINUE
 120        CONTINUE
            IF (IPRHSO.GT.IPRLOC) THEN
               WRITE(LUPRI,'(/A)')' Inactive Fock matrix'
               CALL OUTPUT(FC,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
               WRITE(LUPRI,'(/A)')' Active part of inactive Fock matrix'
               CALL OUTPUT(WRK(KFCAC),1,NASHT,1,NASHT,NASHT,NASHT,1,
     &                    LUPRI)
            END IF
C
C Get transposed integrals
C
C           CALL MTRSP(N2ASHX,N2ASHX,WRK(KHSQSQ),N2ASHX,
C    *                               WRK(KHSQTR),N2ASHX)
            CALL DAXPY(N2ASHX*N2ASHX,DM1,WRK(KHSQSQ),1,WRK(KHSQTR),1)
            DISTYP = 6
            CALL CISIGD(IREFSY, KSYMST, NDREF, KZCONF, WRK(KREFCO),
     *               GP(KZVAR+1),
     *               WRK(KFCAC), WRK(KHSQTR), .FALSE., .FALSE.,
     *               XINDX, ISPIN1, ISPIN2, WRK(KFREE), LFREE)
            IF (IREFSY .EQ. KSYMST) THEN
C              ... remove CREF component of E[2] vector
               T1 = DDOT(KZCONF,WRK(KREFCO),1,GP(1+KZVAR),1)
               CALL DAXPY(KZCONF,(-T1),WRK(KREFCO),1,GP(1+KZVAR),1)
            END IF
            CALL DSCAL(KZCONF,DM1,GP(1+KZVAR),1)
            IF (IPRHSO.GT.IPRLOC) THEN
               WRITE(LUPRI,'(/A)')
     &           ' CSF part of spin-orbit property vector: Y part'
               CALL OUTPUT(GP(1+KZVAR),1,KZCONF,1,1,KZCONF,1,1,LUPRI)
            END IF
C
         ELSE
            CALL DCOPY(KZCONF,GP,1,GP(1+KZVAR),1)
         END IF
 100  CONTINUE
      CALL QEXIT('HSOSIG')
      RETURN
      END
C  /* Deck h2acso */
      SUBROUTINE H2ACSO(H2ACPK,H2SOAC)
C
C 20-Jun-1990 hjaaj
C
#include <implicit.h>
      DIMENSION H2ACPK(NNASHX,NNASHX), H2SOAC(N2ASHX,NASHT,NASHT)
C
      PARAMETER ( DM1 = -1.0D0 )
C
C Used from common blocks:
C  INFORB : NASHT,NNASHX,N2ASHX
C  INFIND : IROW()
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
C
      CALL QENTER('H2ACSO')
C
C     1. Unpack H2ACPK into H2SOAC
C
      DO 140 K = 1,NASHT
         DO 130 L = 1,K-1
            KL = IROW(K) + L
            CALL DSPTSI(NASHT,H2ACPK(1,KL),H2SOAC(1,K,L))
            CALL DAXPY(N2ASHX,DM1,H2SOAC(1,K,L),1,H2SOAC(1,L,K),1)
  130    CONTINUE
         CALL DZERO(H2SOAC(1,K,K),N2ASHX)
  140 CONTINUE
      CALL QEXIT('H2ACSO')
      RETURN
      END
