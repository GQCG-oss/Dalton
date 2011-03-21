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
      DIMENSION GP(*),CMO(*),UDV(*),PV(*),XINDX(*),WRK(LWRK),NEEDSO(6)
      CHARACTER*8 WORD
      DATA NEEDSO/1,1,1,1,1,0/
C
      PARAMETER(D0=0D0, D1=1D0, DM1=-1D0, D2=2D0)
      LOGICAL SO2TRA, FILE_EXISTS, FOPEN, DO1, DO2, LORB, LCON, TRPDEN
      CHARACTER*8 LABEL
      LOGICAL  FNDLAB
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
#include <aovec.h>
#include <priunit.h>
#include <dummy.h>
#include <inforb.h>
#include <inftap.h>
#include <infind.h>
#include <infvar.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <infpri.h>
#include <trhso.h>
#include "rspprp.h"
#include <infhso.h>
#include <infhyp.h>
#include <infsmo.h>
#include <infpp.h>
#include <iratdef.h>
#include <codata.h>
#include <inflr.h>
#include <qrinf.h>
#include <cbihr2.h>
#include "infesr.h"
#include <dftcom.h>
C
      CALL QENTER('HSOCTL')
C
      IF (SOPPA) CALL QUIT('HSOCTL: SOPPA not implemented yet!')
      IPRHSO=MAX(IPRRSP,IPRPP,IPRLR,IPRESR)
      KSYMSO=KSYMOP
      CALL DZERO(GP,KZYVAR)
      HSOFAC=ALPHAC**2/4
      IF (WORD(2:2).EQ.'1') THEN
         IF (LQRGP) THEN
            CALL QRGP(WORD,GP,CMO,XINDX,ANTSYM,WRK,LWRK)
         ELSE
            CALL PRPCTL(WORD,GP,CMO,UDV,PV,XINDX,ANTSYM,WRK,LWRK)
         END IF
         CALL DSCAL(KZYVAR,HSOFAC,GP,1)
         GO TO 9999
      END IF
      ANTSYM = 1.0D0
      IF (WORD(2:2).EQ.'2') THEN
         DO1 = .FALSE.
         DO2 = .TRUE.
      END IF
      IF (WORD(2:2).EQ.' ') THEN
         DO1 = DOSO1
         DO2 = DOSO2
      END IF
C
      IF (IPRHSO.GT.2) THEN
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
      INQUIRE(FILE='HSOGRAD',EXIST=FILE_EXISTS)
      LUHSO = -1
      CALL GPOPEN(LUHSO,'HSOGRAD','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUHSO
      IF (FILE_EXISTS) THEN
C     ... aug07-hjaaj: gfortran doesn't like BACKSPACE(LUHSO) in FNDLAB on empty file!
      IF (FNDLAB(WORD,LUHSO)) THEN
         CALL READT(LUHSO,KZYVAR,GP)
         CALL GPCLOSE(LUHSO,'KEEP')
         GO TO 9999
      END IF
      END IF
C
C ALLOCATE WORK SPACE
C
      KWRK1 = 1
      LWRK1 = LWRK
      CALL MEMGET('REAL', KFC    , NORBT*NORBT,WRK,KWRK1,LWRK1)
      IF (LORB) THEN
         CALL MEMGET('REAL', KFV    , NORBT*NORBT,WRK,KWRK1,LWRK1)
         CALL MEMGET('REAL', KQA    , NORBT*NASHT,WRK,KWRK1,LWRK1)
         CALL MEMGET('REAL', KQB    , NORBT*NASHT,WRK,KWRK1,LWRK1)
      END IF
      CALL MEMGET('REAL', KH1 , NORBT*NORBT,WRK,KWRK1,LWRK1)
      CALL MEMGET('REAL', KH2 , NORBT*NORBT,WRK,KWRK1,LWRK1)
      IF (LCON) THEN
         CALL MEMGET('REAL',  KH2AC , NNASHX*NNASHX,WRK,KWRK1,LWRK1)
      ELSE
         CALL MEMGET('REAL',  KH2AC , 0            ,WRK,KWRK1,LWRK1)
      END IF
C
      IF (IPRHSO.GT.20 .AND. LORB .AND. NASHT.GT.0) THEN
            WRITE(LUPRI,'(/A)')
     *   ' First order density matrix:'
            CALL OUTPUT(UDV,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
            WRITE(LUPRI,'(/A)')
     *   ' Second order density matrix:'
            CALL PRIAC2(PV,NASHT,LUPRI)
      END IF
C
C
C Read AO one-electron property integrals and transform to MO basis.
C
      IF (DO1) THEN
         LABEL = WORD
         LABEL(2:2) = '1'
         KSYMP = -1
         CALL PRPGET (LABEL,CMO,WRK(KH1),KSYMP,ANTSYM,WRK(KWRK1),LWRK1,
     &                IPRHSO)
         IF (KSYMP.NE.KSYMOP) THEN
            WRITE (LUPRI,'(/A/2A/A,2I5/A,F10.2)')
     &           'FATAL ERROR: KSYMOP .ne. KSYMP from PRPGET',
     &           '   Property label  : ',WORD,
     &           '   KSYMOP and KSYMP:',KSYMOP,KSYMP,
     &           '   ANTSYM          :',ANTSYM
            CALL QUIT('KSYMOP .ne. KSYMP from PRPGET')
         END IF
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
C     IF (.NOT.DO2) GOTO 95
      IF (.NOT.TDHF) CALL TRANSFORM_HSO(WORD,CMO,WRK,KWRK1,LWRK1)
C
      IF (TRPLET) THEN
         TRPDEN=.FALSE.
         CALL HSOFCK(WORD,WRK(KFC),WRK(KFV),WRK(KQA),WRK(KQB),
     &               WRK(KH2AC),UDV,PV,CMO,LORB,LCON,WRK,KWRK1,LWRK1)
      ELSE
         TRPDEN=.TRUE.
         CALL MEMGET('REAL',KD,N2ASHX,WRK,KWRK1,LWRK1)
         CALL MEMGET('REAL',KP,2*N2ASHX*N2ASHX,WRK,KWRK1,LWRK1)
         CALL DZERO(WRK(KD),N2ASHX)
         CALL DZERO(WRK(KP),2*N2ASHX*N2ASHX)
         IF (TDHF) THEN
            IF (NASHT.GE.1) THEN
               CALL DUNIT(WRK(KD),NASHT)
            END IF
         ELSE
C
C Need new density matrices
C
            CALL MEMGET('REAL',KCREF,MZCONF(1),WRK,KWRK1,LWRK1)
            CALL GETREF(WRK(KCREF),MZCONF(1))
            CALL RSPDM(IREFSY,IREFSY,MZCONF(1),MZCONF(1),
     &         WRK(KCREF),WRK(KCREF),
     &         WRK(KD),WRK(KP),1,0,.FALSE.,.FALSE.,XINDX,
     &         WRK,KWRK1,LWRK1)
            KP1=KP
            KP2=KP+N2ASHX*N2ASHX
            CALL MTRSP(N2ASHX,N2ASHX,WRK(KP1),N2ASHX,WRK(KP2),N2ASHX)
            CALL MEMREL('HSOCTL<-RSPDM',WRK,KD,KCREF,KWRK1,LWRK1)
         END IF
         CALL HSOFCK(WORD,WRK(KFC),WRK(KFV),WRK(KQA),WRK(KQB),
     &   WRK(KH2AC),WRK(KD),WRK(KP),CMO,LORB,LCON,WRK,KWRK1,LWRK1)
         CALL MEMREL('HSOCTL<-HSOFCK',WRK,KD,KP,KWRK1,LWRK1)
      END IF
 95   CONTINUE
C
C
C Add the one-electron and two electron parts of the inactive Fock matrix
C
      IF (DO1) THEN
         CALL DAXPY(N2ORBX,D1,WRK(KH1),1,WRK(KFC),1)
         IF (TDHF.AND.NASHT.GT.0) THEN
            CALL DAXPY(N2ORBX,-D1,WRK(KH1),1,WRK(KFV),1)
         END IF
      END IF
C
C
C
C Construct orbital property vector
C
      IF (LORB) THEN
         IF (IPRHSO.GT.2) CALL TIMER('START ',ORBSTA,ORBTIM)
         IF (TRPLET) THEN
            CALL HSOORB(.TRUE.,1,DUMMY1,
     *         WRK(KFC),WRK(KFV),WRK(KQA),WRK(KQB),
     *         UDV    ,GP,DUMMY2,TRPDEN)
         ELSE
C
C Quick fix for high spin matrices which do not fit neatly into HSOORB
C
            IF (TDHF.AND.NASHT.GT.0) TRPDEN=.FALSE.
C
            CALL HSOORB(.TRUE.,1,DUMMY1,
     *         WRK(KFC),WRK(KFV),WRK(KQA),WRK(KQB),
     *         WRK(KD),GP,DUMMY2,TRPDEN)
         END IF
C        CALL HSOORB(ONEIND,NSIM,FC,FCX,FVX,QAX,QBX,UDV,
C    *                  EVECS,ZYMAT,TRPDEN)
C
C
         IF (IPRHSO.GT.2) CALL TIMER('RSPORB',ORBSTA,ORBTIM)
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
         IF (IPRHSO.GT.2) CALL TIMER('START ',ORBSTA,ORBTIM)
         CALL HSOSIG(WRK(KFC),WRK(KH2AC),
     *      GP,XINDX,WRK(KWRK1),LWRK1)
C
C        CALL HSOSIG(FC,H2AC, GP,XINDX,WRK,LWRK)
C
C
         IF (IPRHSO.GT.2) CALL TIMER('HSOSIG',ORBSTA,ORBTIM)
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
      ELSE IF (IPRHSO.GT.2) THEN
         IF (LORB) CALL RSPPRO (GP(1+KZCONF),KZVAR,LUPRI)
         IF (LCON) CALL RSPPRC (GP,KZCONF,KZVAR,LUPRI)
      END IF
C
      IF (IPRHSO.GT.2) CALL TIMER('HSOCTL',HSOSTA,HSOTIM)
      IF (X2GRAD) THEN
         CALL MEMREL('X2GRAD',WRK,KFC,KFC,KWRK1,LWRK1)
         CALL MEMGET('INTE',KMJWOP,16*MAXWOP,WRK,KWRK1,LWRK1)
         CALL HEADER('X2GRAD TEST FOR SPIN-ORBIT GRADIENT ELEMENTS',3)
         CALL MEMGET('REAL',KGP2,KZYVAR,WRK,KWRK1,LWRK1)
         CALL DZERO(WRK(KGP2),KZYVAR)
         CALL SETZY(WRK(KMJWOP))
         CALL HSOAL2 (WORD,WRK(KGP2),CMO,UDV,PV,XINDX,WRK(KMJWOP),
     &                WRK(KWRK1),LWRK1)
         WRITE(LUPRI,'(/A)')'               GP1           GP2'
         WRITE(LUPRI,'(/A)')'      ---------------------------------'
         DMAXGP = D0
         DO 77, K = 0,KZYVAR-1
            WRITE(LUPRI,'(10X,2F14.8)') GP(1+K), WRK(KGP2+K)
            DMAXGP = MAX(DMAXGP,ABS(GP(1+K)-WRK(KGP2+K)))
77       CONTINUE
         WRITE(LUPRI,'(/A)')'      ---------------------------------'
         WRITE(LUPRI,'(/A,E20.8)')
     *      'LARGEST DIFFERENCE OF SPIN-ORBIT GRADIENT ELEMENTS' ,
     *      DMAXGP
         CALL MEMREL('X2GRAD',WRK,KGP2,KGP2,KWRK1,LWRK1)
      END IF
      CALL DSCAL(KZYVAR,HSOFAC,GP,1)
C
C Save on file
C
      REWIND LUHSO
      CALL NEWLAB(WORD,LUHSO,LUERR)
      CALL WRITT(LUHSO,KZYVAR,GP)
      CALL MEMREL('HSOCTL',WRK,1,1,KWRK1,LWRK1)
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
#include "infs0.h"
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
      NEWDEF = (WORD .EQ. '*SPIN-O')
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
         IF (IPRHSO.NE.2) WRITE(LUPRI,'(/A/A,I5)')
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
      PARAMETER (D1=1.0D0, D1P5=1.5D0, D2=2.0D0)
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
      IF (TRPLET) THEN
         COUFAC=D1
      ELSE
         COUFAC=D2
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
            FAC = DISFAC*COUFAC
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
                  CALL DGEMM('N','N',NORBP,1,NASHX,1.0D0,
     &                       H2(IOFFP+1,IOFFXA+1),NORBT,
     &                       UDV(IASHX+1,NCIW),NASHT,0.0D0,WRK(1),NORBP)
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
                     CALL DGEMM('N','N',NDIMP,1,NASHX,1.0D0,
     &                          H2(IOFFP+1,IOFFXA+1),NORBT,
     &                          UDV(IASHX+1,NDIW),NASHT,0.0D0,
     &                          WRK(1),NDIMP)
                     CALL DAXPY(NDIMP,FAC,WRK(1),1,
     *                          FVSO(IOFFP+1,ICI),1)
                     CALL DAXPY(NDIMP,-FAC,WRK(1),1,
     *                          FVSO(ICI,IOFFP+1),NORBT)
                     IF (IPRHSO.GT.IPRFSO) THEN
                       WRITE(LUPRI,'(A,I3,A,I3,I5)')' PC contribution ',
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
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL', KFCAC , N2ASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL', KREFCO , NDREF,WRK,KFREE,LFREE)
      CALL MEMGET('REAL', KHSQSQ , N2ASHX*N2ASHX,WRK,KFREE,LFREE)
      CALL DZERO(WRK(KFCAC),N2ASHX)
      CALL DZERO(WRK(KREFCO),NDREF)
      CALL DZERO(WRK(KHSQSQ),N2ASHX*N2ASHX)
      IF (TESTZY) THEN
         CALL MEMGET('REAL', KHSQTR , N2ASHX*N2ASHX,WRK,KFREE,LFREE)
         CALL DZERO(WRK(KHSQTR),N2ASHX*N2ASHX)
      END IF
      IADINT = -1
C
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
      CALL MEMREL('HSOSIG',WRK,1,1,KFREE,LFREE)
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
      SUBROUTINE TRANSFORM_HSO(WORD,CMO,WRK,KWRK1,LWRK1)
#include <implicit.h>
      CHARACTER*8 WORD
      DIMENSION  WRK(*), CMO(*)
#include <dummy.h>
#include "rspprp.h"
#include <infhso.h>
#include <infhyp.h>
#include <infpp.h>
#include <infsmo.h>
#include <inftap.h>
#include <priunit.h>
#include <trhso.h>
#include <wrkrsp.h>
      LOGICAL FOPEN
      LOGICAL MOEXIS, SO2TRA
C
C Set spin-orbit component ILXYZ for TRHSO
C
      CALL QENTER('TRANSFORM_HSO')
      IF (WORD(1:1).EQ.'X') THEN
         ILXYZ = 1
      ELSE IF (WORD(1:1).EQ.'Y') THEN
         ILXYZ = 2
      ELSE IF (WORD(1:1).EQ.'Z') THEN
         ILXYZ = 3
      ELSE
         CALL QUIT('Wrong property in TRANSFORM_HSO, WORD = '//WORD)
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
      SO2TRA = .TRUE.
C
C If .OLDTRA has been specified check that MO2SOINT exists
C If it does not exist transform AO2SOINT as usual
C
      IF (OLDTRA) THEN
         INQUIRE(FILE='MOHSOINT',EXIST=MOEXIS)
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
         INQUIRE(FILE='MOHSOINT',OPENED=FOPEN)
         IF (.NOT.FOPEN) THEN
            CALL DAOPEN(LUMHSO,'MOHSOINT')
         END IF
      END IF
      CALL QEXIT('TRANSFORM_HSO')
      END
      SUBROUTINE HSOFCK(WORD,FC,FV,QA,QB,H2AC,UDV,PV,CMO,LORB,LCON,
     &   WRK,KWRK1,LWRK1)
#include <implicit.h>
      CHARACTER*8 WORD
      DIMENSION FC(*), FV(*), QA(*), QB(*), H2AC(*), UDV(*), PV(*),
     &   CMO(*), WRK(*)
      LOGICAL LORB, LCON
#include <infrsp.h>
#include <inforb.h>
      CALL QENTER('HSOFCK')
      IF (TDHF) THEN
         IF (NASHT.GT.0) THEN
            CALL DZERO(QA,NORBT*NASHT)
            CALL DZERO(QB,NORBT*NASHT)
         END IF
         CALL HSOFCKAO(WORD,CMO,UDV,FC,FV,WRK,KWRK1,LWRK1)
      ELSE
         CALL HSOFCKMO(FC,FV,QA,QB,H2AC,UDV,PV,LORB,LCON,
     &   WRK,KWRK1,LWRK1)
      END IF
      CALL QEXIT('HSOFCK')
      END 
C
C /* Deck hsofckao */
C
      SUBROUTINE HSOFCKAO(WORD,CMO,UDV,FC,FV,WRK,KFREE,LFREE)
      IMPLICIT NONE
      CHARACTER*(*) WORD
      REAL*8    CMO(*), UDV(*), FC(*), FV(*), WRK(*)
      INTEGER KFREE,LFREE
C
C Build HSO inactive fockmatrix from AO integrals
C If DIRFCK=.TRUE. call TWOINT
C else read from AO2SOINT
C
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <dummy.h>
#include <infinp.h>
#include <wrkrsp.h>
#include <inforb.h>
#include <aovec.h>
#include <infrsp.h>
#include <dftcom.h>
#include <symmet.h>
#ifdef VAR_DEBUG
#define DEBUG .TRUE.
#else
#define DEBUG .FALSE.
#endif
      REAL*8 D1
      PARAMETER (D1=1.0D0)
      INTEGER NDMAT,KDMAO,KFMAO,KDMSO,I
      INTEGER KDC,KDV,KFC,KFV,KFMO,KFAO
      INTEGER KJSTRS , KNPRIM , KNCONT , KIORBS , KJORBS , NPAO , KKORBS
      INTEGER I2TYP
      INTEGER IPRTWO,IPRNTA,IPRNTB,IPRNTC,IPRNTD
      INTEGER IRNTYP , NUMDIS, IFCTYP(2), ISYMDM(2)
      REAL*8  TIMEND,TIMSTR,SO1WAL,SO1CPU
      REAL*8  SO2WAL,SO2CPU,SOCPU,SOWAL
      REAL*8  HFXSAV
      LOGICAL TKTIME,RTNTWO
      LOGICAL TRIPLET(2)
      CALL QENTER('HSOFCKAO')

      NDMAT=1
      IF (NASHT.GT.0) NDMAT = 2
      CALL MEMGET('REAL',KFC,NDMAT*N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDC,NDMAT*N2BASX,WRK,KFREE,LFREE)
      KFV = KFC + N2BASX
      KDV = KDC + N2BASX
C
C Full density in 1, spin (=active) density in 2, swap for singlet
C
      CALL GTDMSO(UDV,CMO,WRK(KDC),WRK(KDV),WRK(KFREE))
      IF (NASHT.GT.0) THEN
         CALL DAXPY(N2BASX,D1,WRK(KDV),1,WRK(KDC),1)
         IF (.NOT.TRPLET) THEN
            CALL DSWAP(N2BASX,WRK(KDC),1,WRK(KDV),1)
         END IF
      END IF
      IF (DEBUG) THEN
         CALL MAOPRI(WRK(KDC),'HSOFCKAO:Input  D1')
         IF (NASHT.GT.0) THEN
            CALL MAOPRI(WRK(KDC+N2BASX),'HSOFCKAO:Input  D2')
         END IF
      END IF
      IF (DIRFCK) THEN
         IF (WORD(1:1).EQ.'X') IFCTYP(1)=1
         IF (WORD(1:1).EQ.'Y') IFCTYP(1)=2
         IF (WORD(1:1).EQ.'Z') IFCTYP(1)=3
         IF (TRPLET) IFCTYP(1) = -IFCTYP(1)
         IFCTYP(2) = - IFCTYP(1)
C
C        Transform density to AO basis 
C
         CALL DSOTAO(WRK(KDC),WRK(KFC),NBAST,0,IPRRSP)
         IF (NASHT.GT.0) THEN
            CALL DSOTAO(WRK(KDV),WRK(KFV),NBAST,0,IPRRSP)
         END IF
         CALL DCOPY(NDMAT*N2BASX,WRK(KFC),1,WRK(KDC),1)
C
C Setup for TWOINT
C
         NPAO = MXSHEL*MXAOVC
         CALL MEMGET('INTE',KJSTRS,NPAO*2,WRK,KFREE,LFREE)     
         CALL MEMGET('INTE',KNPRIM,NPAO*2,WRK,KFREE,LFREE)     
         CALL MEMGET('INTE',KNCONT,NPAO*2,WRK,KFREE,LFREE)     
         CALL MEMGET('INTE',KIORBS,NPAO  ,WRK,KFREE,LFREE)     
         CALL MEMGET('INTE',KJORBS,NPAO  ,WRK,KFREE,LFREE)     
         CALL MEMGET('INTE',KKORBS,NPAO  ,WRK,KFREE,LFREE)     
         CALL PAOVEC(WRK(KJSTRS),WRK(KNPRIM),WRK(KNCONT),
     &               WRK(KIORBS),WRK(KJORBS),WRK(KKORBS),0,
     &               .FALSE.,IPRRSP)
         CALL MEMREL('HERINT.PAOVEC',WRK,1,KJORBS,KFREE,LFREE)
         CALL TIMER('START ',TIMSTR,TIMEND)
         CALL GETTIM(SO1CPU,SO1WAL)
         I2TYP = 0
         IRNTYP = -20
         IPRTWO = 0
         TKTIME = .FALSE.
         CALL DZERO(WRK(KFC),NDMAT*N2BASX)
C        Always full exchange for spin-orbit integrals:
         HFXSAV=HFXFAC
         HFXFAC=D1
         CALL TWOINT(WRK(KFREE),LFREE,WRK(KFREE),
     &               WRK(KFC),WRK(KDC),NDMAT,
     &               IDUMMY, IFCTYP,
     &               DUMMY,IDUMMY,NUMDIS,1,IRNTYP,0,0,.TRUE.,.TRUE.,
     &               .FALSE.,TKTIME,IPRTWO,IPRNTA,IPRNTB,IPRNTC,
     &               IPRNTD,RTNTWO,IDUMMY,I2TYP,WRK(KJSTRS),
     &               WRK(KNPRIM),WRK(KNCONT),WRK(KIORBS),
     &               IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &               .FALSE.)
         HFXFAC=HFXSAV
         CALL MEMREL('HSOFCKAO.TWOINT',WRK,KFC,KJSTRS,KFREE,LFREE)
         IF (DEBUG) THEN
            CALL MAOPRI(WRK(KFC),':TWOINT:Calculated F1')
            IF (NASHT.GT.0) THEN
               CALL MAOPRI(WRK(KFC+N2BASX),':TWOINT:Calculated F2')
            END IF
         END IF
C
C  Transform to SO
C
         ISYMDM(1)=0
         ISYMDM(2)=0
         CALL SKLFCK(WRK(KFC),DUMMY,WRK(KFREE),LFREE,
     &               IPRTWO,.FALSE.,
     &               .FALSE.,.FALSE.,.FALSE.,.TRUE.,IDUMMY,.TRUE.,NDMAT,
     &               ISYMDM,IFCTYP,IDUMMY,.TRUE.)
C
         IF (DEBUG) THEN
            CALL MAOPRI(WRK(KFC),':SKLFCK:Calculated F1')
            IF (NASHT.GT.0) THEN
               CALL MAOPRI(WRK(KFC+N2BASX),':SKLFCK:Calculated F2')
            END IF
         END IF
         CALL MEMCHK('HSOFCKAO.SKLFCK',WRK,KFC)
         CALL GETTIM(SO2CPU,SO2WAL)
         SOCPU=SO2CPU-SO1CPU
         SOWAL=SO2WAL-SO1WAL
         CALL TIMER('TWOINT',TIMSTR,TIMEND)
         CALL FLSHFO(LUPRI)
         WRITE(LUPRI,'(2(/A),2(/A,F5.0,A))')
     &      '   Two-electron spin-orbit integrals',
     &   '   =================================',
     &   '   Spin-orbit 2-electron CPU  time ',SOCPU,' seconds',
     &   '   Spin-orbit 2-electron wall time ',SOWAL,' seconds'
      ELSE
C
C Read AO integrals from disk (AO2SOINT)
C
         TRIPLET(1) = TRPLET
         TRIPLET(2) = .NOT.TRPLET
         CALL DZERO(WRK(KFC),NDMAT*N2BASX)
         CALL GETFMATX(WORD,TRIPLET,WRK(KFC),WRK(KDC),NDMAT)
      END IF
      IF (DEBUG) THEN
         CALL MAOPRI(WRK(KFC),'HSOFCKAO:Calculated F1')
         IF (NASHT.GT.0) THEN
            CALL MAOPRI(WRK(KFC+N2BASX),'HSOFCKAO:Calculated F2')
         END IF
      END IF
C
C Reorder matrices to match gradient routine
C (fa+fb,fa-fb)/2 = (fc,fo-fc) -> (fo,fc-fo)
C
      IF (NASHT.GT.0) THEN
         CALL DAXPY(N2BASX,D1,WRK(KFV),1,WRK(KFC),1) !now (fo,fo-fc)
C        CALL DSWAP(N2BASX,WRK(KFC),1,WRK(KFV),1)
         CALL DSCAL(N2BASX,-D1,WRK(KFV),1) !now (fo,fc-fo)
      END IF
C
      CALL DZERO(FC,N2ORBX)
      CALL AOTOMO(WRK(KFC),FC,CMO,KSYMOP,WRK(KFREE),LFREE)
      IF (NASHT.GT.0) THEN
         CALL DZERO(FV,N2ORBX)
         CALL AOTOMO(WRK(KFV),FV,CMO,KSYMOP,WRK(KFREE),LFREE)
      END IF
      IF (DEBUG) THEN
         CALL MAOPRI(WRK(KFC),'HSOFCKAO:FC(AO)')
         IF (NASHT.GT.0) CALL MAOPRI(WRK(KFV),'HSOFCKAO:FC-FO(AO)')
         CALL MMOPRI(FC,'HSOFCKAO:FO(MO)')
         IF (NASHT.GT.0) CALL MMOPRI(FV,'HSOFCKAO:FC-FO(MO)')
      END IF
      CALL MEMREL('HSOFCKAO',WRK,KFC,KFC,KFREE,LFREE)
      CALL QEXIT('HSOFCKAO')
      END
C  /* Deck hsofckmo  */
      SUBROUTINE HSOFCKMO(FC,FV,QA,QB,H2AC,UDV,PV,LORB,LCON,
     &   WRK,KWRK1,LWRK1)
#include <implicit.h>
      DIMENSION FC(*), FV(*), QA(*), QB(*), H2AC(*), UDV(*), PV(*), 
     &   WRK(*)
      LOGICAL LORB, LCON
#include <infhso.h>
#include <maxorb.h>
#include <maxash.h>
#include <infind.h>
#include <inforb.h>
#include <priunit.h>
#include <wrkrsp.h>
      PARAMETER (D1=1.0D0)
      INTEGER NEEDSO(6)
      DATA NEEDSO/1,1,1,1,1,0/
      CALL QENTER('HSOFCK')
      CALL DZERO(FC,N2ORBX)
      IF (LORB) THEN
         CALL DZERO(FV,N2ORBX)
         CALL DZERO(QA,NORBT*NASHT)
         CALL DZERO(QB,NORBT*NASHT)
      END IF
      IF (LCON) THEN
         CALL DZERO(H2AC,NNASHX*NNASHX)
      END IF
      IF (.NOT.DOSO2) THEN
         CALL QEXIT('HSOFCK')
         RETURN
      END IF
      CALL MEMGET('REAL',KH2,N2ORBX,WRK,KWRK1,LWRK1)
      CALL MEMGET('REAL',KPV12,N2ASHX,WRK,KWRK1,LWRK1)
      CALL MEMGET('REAL',KPV21,N2ASHX,WRK,KWRK1,LWRK1)
      IDIST = 0
      KFREE = 1
      LFREE = LWRK1
C
C Read distributions
C
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
      CALL FSOMU(IC,ID,WRK(KH2),FC,FV,UDV,LORB,
     * WRK(KWRK1-1+KFREE),LFREE)
C     CALL FSOMU(ICI1,IDI1,H2,FCSO,FVSO,UDV,LORB,WRK,LWRK)
C
C
C   Construct spin-orbit Q matrices
C
      IF (NASHT.GT.0 .AND. LORB)
     *   CALL QSOMU(IC,ID,QA,QB,WRK(KH2),
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
            KH2XY = 1 + (NCDW-1)*NNASHX
            CALL ADH2AC(H2AC(KH2XY),WRK(KH2),KCDSYM)
         ELSE
            NCDW = IROW(NDW) + NCW
            KH2XY = 1 + (NCDW-1)*NNASHX
            CALL DAXPY(N2ORBX,-D1,WRK(KH2),1,WRK(KWRK1-1+KFREE),1)
            CALL ADH2AC(H2AC(KH2XY),WRK(KWRK1-1+KFREE),KCDSYM)
         END IF
      END IF
      GOTO 90
 95   CONTINUE
C
      IF (IPRHSO.GT.10) THEN
         WRITE(LUPRI,'(//A)') ' Final Inactive Fock matrix'
         CALL OUTPUT(FC,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
C
         IF (LORB .AND. NASHT.GT.0) THEN
            WRITE(LUPRI,'(//A)') ' Final Active Fock matrix'
            CALL OUTPUT(FV,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
C
            WRITE(LUPRI,'(//A)') ' Spin-orbit QA matrix'
            CALL OUTPUT(QA,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
            WRITE(LUPRI,'(//A)') ' Spin-orbit QB matrix'
            CALL OUTPUT(QB,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
         END IF
C
         IF (LCON) THEN
            DO 71 K=1,NASHT
               DO 72 L=1,K
                  KL = IROW(K)+L
                  KLOFF = (KL-1)*NNASHX
                  WRITE(LUPRI,'(/A,2I3,A)')
     *         ' Active integral distribution (**', K,L,')'
                  CALL OUTPAK(H2AC(1+KLOFF),NASHT,1,LUPRI)
 72            CONTINUE
 71         CONTINUE
         END IF
      END IF
      CALL MEMREL('HSOFCK',WRK,KH2,KH2,KWRK1,LWRK1)
      CALL QEXIT('HSOFCK')
      END
C  /* Deck hsoorb */
      SUBROUTINE HSOORB(ONEIND,NSIM,FC,FCX,FVX,QAX,QBX,UDV,
     *                  EVECS,ZYMAT,TRPDEN)
C
C WRITTEN 14-FEB 1986
C adapted 5-Apr 2001
C
C PURPOSE:
C  1)DISTRIBUTE INACTIVE FCX AND ACTIVE FVX FOCK MATRICES AND QAX AND
C    QBX MATRICES INTO ORBITAL PART OF LINEAR TRANSFORMED VECTOR
C
C                  ( [ E(K,L) , H ] )    K<L
C       E[2]*N = - (                )
C                  ( [ E(L,K) , H ] )    K>L
C    X MAY REFER TO EITHER ONE INDEX TRANSFORMED MATRICES (N IS A
C    ORBITAL TRIAL VECTOR ) OR TRANSITION DENSITY MATRICES (N IS A
C    CONFIGURATION TRIAL VECTOR). FOR CONFIGURATION TRIAL VECTORS
C    FCX IS IDENTICALLY ZERO. FURTHER OVERLAP IS ZERO BECAUSE TRIAL
C    VECTORS ARE CHOSEN ORTHOGONAL TO REFERENCE STATE.
C
C  2)CREATE LINEAR TRANSFORMED VECTOR S[2]*N FOR N EQUAL TO EITHER A
C    ORBITAL OR A CONFIGURATION TRIAL VECTOR
C
C    ONEIND = .TRUE. FOR A ORBITAL TRIAL VECTOR
C
C                ******************************
C
C    [ E(P,Q) , H ] =
C
C    ACTIVE-INACTIVE (P,Q) = (T,I)
C                             L,K
C    FCX(I,X)*DV(T,X)-2*FCX(I,T)-2FVX(I,T)+QBX(I,T)
C        K       L          K,L       K,L      K,L
C
C    INACTIVE-ACTIVE         (I,U)
C                             K,L
C    -FCX(X,I)*DV(X,U)+2*FCX(U,I)+2*FVX(U,I)-QAX(I,U)
C           K       L        L,K        L,K      K,L
C
C    ACTIVE-ACTIVE           (T,U)
C                             K,L
C                             L,K
C    FCX(U,X)*DV(T,X)-FCX(X,T)*DV(X,U)+QBX(U,T)-QAX(T,U)
C        L       K          K       L      L,K      K,L
C        K       L          L       K      K,L      L,K
C
C    ACTIVE-SECUNDARY        (T,A)
C                             K,L
C    FCX(A,X)*DV(T,X)+QBX(A,T)
C        L       K        L,K
C
C    SECUNDARY-ACTIVE        (A,U)
C                             L,K
C    -FCX(X,A)*DV(X,U)-QAX(A,U)
C           L       K      L,K
C
C    INACTIVE-SECUNDARY      (I,A)
C                             K,L
C    2*FCX(A,I)+2FVX(A,I)
C          L,K       L,K
C
C    SECUNDARY-INACTIVE      (A,J)
C                             L,K
C    -2*FCX(J,A)-2*FVX(J,A)
C           K,L        K,L
C
C
#include <implicit.h>
#include <priunit.h>
C
      LOGICAL ONEIND, TRPDEN
C
      DIMENSION FC(*),FCX(NORBT,NORBT,*),FVX(NORBT,NORBT,*)
      DIMENSION QAX(NORBT,NASHDI,*),QBX(NORBT,NASHDI,*)
      DIMENSION UDV(NASHDI,NASHDI,*), EVECS(KZYVAR,*)
      DIMENSION ZYMAT(NORBT,NORBT,*)
C
C  INFDIM : NASHDI
C
#include <maxorb.h>
#include <maxash.h>
#include <infvar.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infpri.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
C -- local constants
C
      PARAMETER ( D0 = 0.0D0 , D2 =2.0D0 , DM1 = -1.0D0 )
C
C
      KYCONF = KZCONF + KZVAR
      IF (TRPDEN) THEN
         DC = D0
      ELSE
         DC = D2
      END IF
C
C DISTRIBUTE FOCK MATRICES IN EVECS
C
         KSYM1 = 0
         DO 1300 IG = 1,KZWOPT
            K     = JWOP(1,IG)
            L     = JWOP(2,IG)
            KSYM  = ISMO(K)
            LSYM  = ISMO(L)
            IF( KSYM.NE.KSYM1 ) THEN
               KSYM1 = KSYM
               NORBK = NORB(KSYM)
               IORBK = IORB(KSYM)
               NASHK = NASH(KSYM)
               NISHK = NISH(KSYM)
               IASHK = IASH(KSYM)
               IIORBK= IIORB(KSYM)
               IORBL = IORB(LSYM)
               NASHL = NASH(LSYM)
               NISHL = NISH(LSYM)
               NORBL = NORB(LSYM)
               IASHL = IASH(LSYM)
               IIORBL= IIORB(LSYM)
            END IF
            NK    = K - IORBK
            NL    = L - IORBL
            ITYPK = IOBTYP(K)
            ITYPL = IOBTYP(L)
            IF ( ITYPK.EQ.JTINAC )THEN
               DO 2000 ISIM = 1 ,NSIM
                  IF (ONEIND) THEN
                     EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                       + DC * FCX(L,K,ISIM)
                     EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                       - DC * FCX(K,L,ISIM)
                     IF ( NASHT . GT . 0 ) THEN
                        EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                       + D2 * FVX(L,K,ISIM)
                        EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                       - D2 * FVX(K,L,ISIM)
                     END IF
                  ELSE
                     EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                       + D2 * FVX(L,K,ISIM)
                     EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                       - D2 * FVX(K,L,ISIM)
                  END IF
                  IF ( ITYPL.EQ.JTACT ) THEN
                     NWL = ISW(L) - NISHT
                     IF (ONEIND) THEN
                        EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                  - DDOT(NASHL,FCX(IORBL+NISHL+1,K,ISIM),1,
     *                                        UDV(IASHL+1,NWL,1),1)
                        DO 730 IX = 1,NASHL
                           EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                     + FCX(K,IORBL+NISHL+IX,ISIM)*
     *                                        UDV(NWL,IASHL+IX,1)
 730                    CONTINUE
                     ELSE
                        TEMP1 = D0
                        TEMP2 = D0
                        NX  = NISHK
                        DO 825 NWX = IASHK+1,IASHK+NASHK
                           NX = NX + 1
                           IF (NX.LE.NK) THEN
                              FCXK = FC(IIORBK+IROW(NK)+NX)
                           ELSE
                              FCXK = FC(IIORBK+IROW(NX)+NK)
                           END IF
                           TEMP1 = TEMP1 + FCXK * UDV(NWX,NWL,ISIM)
                           TEMP2 = TEMP2 + FCXK * UDV(NWL,NWX,ISIM)
 825                    CONTINUE
                        EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                                          - TEMP1
                        EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                                          + TEMP2
                     END IF
                     EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                     - QAX(K,NWL,ISIM)
                     EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                     + QBX(K,NWL,ISIM)
                  END IF
 2000          CONTINUE
            ELSE
              IF (ITYPL.EQ.JTACT) THEN
                  NWL = ISW(L) - NISHT
                  NWK = ISW(K) - NISHT
                  DO 2020 ISIM=1,NSIM
                     IF (ONEIND) THEN
                        EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                  -DDOT(NASHL,FCX(IORBL+NISHL+1,K,ISIM),1,
     *                                        UDV(IASHL+1,NWL,1),1)
                        DO 740 IX = 1,NASHL
                           EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                        +FCX(K,IORBL+NISHL+IX,ISIM)*
     *                                        UDV(NWL,IASHL+IX,1)
 740                    CONTINUE
                        EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                     -DDOT(NASHK,FCX(IORBK+NISHK+1,L,ISIM),1,
     *                                        UDV(IASHK+1,NWK,1),1)
                        DO 750 IX = 1,NASHK
                           EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                        +FCX(L,IORBK+NISHK+IX,ISIM)*
     *                                        UDV(NWK,IASHK+IX,1)
 750                    CONTINUE
                     ELSE
                        TEMP1 = D0
                        TEMP2 = D0
                        NX  = NISHL
                        DO 835 NWX = IASHL+1,IASHL+NASHL
                           NX = NX + 1
                           IF (NX.LE.NL) THEN
                              FCXL = FC(IIORBL+IROW(NL)+NX)
                           ELSE
                              FCXL = FC(IIORBL+IROW(NX)+NL)
                           END IF
                           TEMP1 = TEMP1 + FCXL * UDV(NWX,NWK,ISIM)
                           TEMP2 = TEMP2 + FCXL * UDV(NWK,NWX,ISIM)
 835                    CONTINUE
                        EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                                          + TEMP2
                        EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                                          - TEMP1
                        TEMP1 = D0
                        TEMP2 = D0
                        NX  = NISHK
                        DO 845 NWX = IASHK+1,IASHK+NASHK
                           NX = NX + 1
                           IF (NX.LE.NK) THEN
                              FCXK = FC(IIORBK+IROW(NK)+NX)
                           ELSE
                              FCXK = FC(IIORBK+IROW(NX)+NK)
                           END IF
                           TEMP1 = TEMP1 + FCXK * UDV(NWX,NWL,ISIM)
                           TEMP2 = TEMP2 + FCXK * UDV(NWL,NWX,ISIM)
 845                    CONTINUE
                        EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                                          - TEMP1
                        EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                                          + TEMP2
                     END IF
                     EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                  +QBX(L,NWK,ISIM) -QAX(K,NWL,ISIM)
                     EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                  +QBX(K,NWL,ISIM) - QAX(L,NWK,ISIM)
 2020             CONTINUE
               ELSE
                  NWK = ISW(K) - NISHT
                  DO 2030 ISIM=1,NSIM
                     IF (ONEIND) THEN
                        EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                     -DDOT(NASHK,FCX(IORBK+NISHK+1,L,ISIM),1,
     *                                        UDV(IASHK+1,NWK,1),1)
                        DO 760 IX = 1,NASHK
                           EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                        +FCX(L,IORBK+NISHK+IX,ISIM)*
     *                                        UDV(NWK,IASHK+IX,1)
 760                    CONTINUE
                     ELSE
                        TEMP1 = D0
                        TEMP2 = D0
                        NX  = NISHL
                        DO 860 NWX = IASHL+1,IASHL+NASHL
                           NX = NX + 1
                           IF (NX.LE.NL) THEN
                              FCXL = FC(IIORBL+IROW(NL)+NX)
                           ELSE
                              FCXL = FC(IIORBL+IROW(NX)+NL)
                           END IF
                           TEMP2 = TEMP2 + FCXL * UDV(NWX,NWK,ISIM)
                           TEMP1 = TEMP1 + FCXL * UDV(NWK,NWX,ISIM)
 860                    CONTINUE
                        EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                                          + TEMP1
                        EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                                          - TEMP2
                     END IF
                     EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *                     -QAX(L,NWK,ISIM)
                     EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *                     +QBX(L,NWK,ISIM)
 2030             CONTINUE
               ENDIF
            ENDIF
 1300    CONTINUE
C
C CHANGE SIGN ON ORBITAL PART OF LINEAR TRANSFORMATION E[2]*N
C
      DO 3000 ISIM = 1,NSIM
         CALL DSCAL(KZWOPT,DM1,EVECS(KZCONF+1,ISIM),1)
         CALL DSCAL(KZWOPT,DM1,EVECS(KYCONF+1,ISIM),1)
 3000 CONTINUE
C
      IF (IPRRSP.GT.110) THEN
         WRITE(LUPRI,*)' LINEAR TRANSFORMATION WITH  E(2)'
         CALL OUTPUT(EVECS,1,KZYVAR,1,NSIM,KZYVAR,NSIM,1,LUPRI)
      END IF
C
C END OF RSPORB
C
      RETURN
      END
C  /* Deck hsoal2 */
      SUBROUTINE HSOAL2 (WORD,GP,CMO,UDV,PV,XINDX,MJWOP,WRK,LWRK)
C
C Alternative routine for constructing spin-orbit gradient
C by way of X2HSO Case 2
C
#include <implicit.h>
      DIMENSION GP(*), CMO(*),UDV(*),PV(*),XINDX(*),WRK(*)
      CHARACTER*8 WORD
C
      PARAMETER (D1=1.0D0, DM1=-1.0D0, DH=0.5D0)
      LOGICAL DO1,DO2,LORB,LCON,LREFST,TRPSAV
      CHARACTER*8 LABEL
C
C Used from common blocks
C
C INFINP: FLAG()
C INFORB: NORBT...
C INFVAR: NCONF
C WRKRSP: KZYVAR...,KSYMOP
C INFRSP: NCREF,IREFSY,DIROIT,SOPPA
C INFHSO: DOSO1
C INFTAP: LUINDX, LUMHSO
C
#include <priunit.h>
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infvar.h>
#include <qrinf.h>
      DIMENSION MJWOP(2,MAXWOP,8)
#include <infhso.h>
#include <trhso.h>
#include <inftap.h>
C
      CALL QENTER('HSOAL2')
      CALL HEADER('HSOAL2: TEST OF SPIN-ORBIT GRADIENT',1)
      IF (SOPPA) CALL QUIT('HSOAL2: SOPPA not implemented yet!')
      IF (.NOT.TDHF .AND. .NOT.FLAG(27))
     &   CALL QUIT('HSOAL2 is only implemented for .DETERMINANTS')
      IF (.NOT.TDHF .AND. NCONF.NE.NCREF)
     &   CALL QUIT('HSOAL2 inconsistency : NCONF.ne.NCREF')
C
C Initialise MZYVAR... for symmetries KSYMOP and 1 (from INFVAR)
C
C Symmetry KSYMOP:
C
      CALL SETZY(MJWOP)
C
C Symmetry 1
C
      LUINDX = -1
      CALL GPOPEN(LUINDX,'LUINDF','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUINDX
      CALL MOLLAB('EXOPSYM1',LUINDX,LUERR)
      READ (LUINDX) IWOPT,IWOP,IWOPH
      CALL GPCLOSE(LUINDX,'KEEP')
      IVAR = IWOPT + NCONF
      MZVAR(1)  = IVAR
      MZCONF(1) = NCONF
      MZWOPT(1) = IWOPT
      MZYVAR(1) = 2*IVAR
      MZYCON(1) = 2*NCONF
      MZYWOP(1) = 2*IWOPT
      WRITE(LUPRI,'(/A)')    ' Number of variables symmetry 1'
      WRITE(LUPRI,'( A,I5)') '  Rotations:     ',MZWOPT(1)
      WRITE(LUPRI,'( A,I5)') '  Configurations:',MZCONF(1)
      WRITE(LUPRI,'(/A,I2)') ' Number of variables symmetry',KSYMOP
      WRITE(LUPRI,'( A,I5)') '  Rotations:     ',MZWOPT(KSYMOP)
      WRITE(LUPRI,'( A,I5)') '  Configurations:',MZCONF(KSYMOP)
      WRITE(LUPRI,'(/A)')    '  Rotational indices     K    L'
      DO 30 I=1,MZWOPT(KSYMOP)
30       WRITE(LUPRI,'(3I5)') I,MJWOP(1,I,KSYMOP),MJWOP(2,I,KSYMOP)
C
C Allocate work space
C
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KVEC,2*IVAR,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KFI,N2ORBX,WRK,KFREE,LFREE)
C     KFI  = KVEC + 2*IVAR
C     KFREE = KFI  + N2ORBX
C     LFREE = LWRK - KFREE + 1
C     IF (LFREE.LT.0) CALL ERRWRK('HSOAL2',KFREE-1,LWRK)
C
C Convention: X1SPNORB: one-electron parts only
C             X2SPNORB: two-electron parts only
C             X SPNORB: both (default)
C                       overridden by SO1ONLY or SO2ONLY in HSOINP
C
      IF (WORD(2:2).EQ.'1' .OR. WORD(2:2).EQ.' '.AND.DOSO1) THEN
         LABEL = WORD
         LABEL(2:2) = '1'
         KSYMP = -1
         CALL PRPGET(LABEL,CMO,WRK(KFI),KSYMP,ANTSYM,
     &               WRK(KFREE),LFREE,IPRHSO)
         IF (KSYMP.NE.KSYMOP) THEN
            WRITE (LUPRI,'(/A/2A/A,2I5/A,F10.2)')
     &           'FATAL ERROR: KSYMOP .ne. KSYMP from PRPGET',
     &           '   Property label  : ',LABEL,
     &           '   KSYMOP and KSYMP:',KSYMOP,KSYMP,
     &           '   ANTSYM          :',ANTSYM
            CALL QUIT('KSYMOP .ne. KSYMP from PRPGET')
         END IF
      ELSE
         CALL DZERO(WRK(KFI),N2ORBX)
      END IF
C
C
C Set up the vector to be  N_o = ( 0 ) , N_c = ( - Ref )
C                                ( 0 )         (   Ref )
C
      CALL DZERO(WRK(KVEC),2*IVAR)
      CALL GETREF(WRK(KVEC + IVAR),NCONF)
      CALL DAXPY(NCONF,DM1,WRK(KVEC + IVAR),1,WRK(KVEC),1)
      WRITE(LUPRI,'(/A)') 'REFERENCE CI  VECTOR'
      CALL OUTPUT(WRK(KVEC + IVAR),1,NCONF,1,1,NCONF,1,1,LUPRI)
      WRITE(LUPRI,'(/A)') 'CALLING HSO2CR WITH VECTOR'
      CALL RSPPRC(WRK(KVEC),NCONF,IVAR,LUPRI)
      CALL RSPPRW(WRK(KVEC+NCONF),MJWOP,IWOPT,1,IVAR,LUPRI)
      ISYMV = 1
      ISPINV = 0
      IKLVL = 0
      DIROIT = .TRUE.
C
C Gradient setup
C
      IGRSYM = KSYMOP
      IF (TRPLET) THEN
         IGRSPI = 1
      ELSE
         IGRSPI = 0
      END IF
C
C Operator setup
C
      IOPSYM = KSYMOP
      IOPSPI = 1
C
C Orbital gradient
C
      LORB = KZWOPT.GT.0
      LPV = 2*N2ASHX*N2ASHX
      CALL MEMGET('REAL',KD,N2ASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KP,LPV,WRK,KFREE,LFREE)
      KP1=KP
      KP2=KP+N2ASHX*N2ASHX
      IF (TRPLET) THEN
         CALL DCOPY(N2ASHX,UDV,1,WRK(KD),1)
         CALL DCOPY(LPV,PV,1,WRK(KP),1)
         CALL IPSET(0,0,1,1)
      ELSE
         IF (TDHF) THEN
            IF (NASHT.EQ.1) THEN
               WRK(KD)=1.0D0
            END IF
         ELSE
            CALL MEMGET('REAL',KCREF,MZCONF(1),WRK,KFREE,LFREE)
            CALL GETREF(WRK(KCREF),MZCONF(1))
            CALL RSPDM(IREFSY,IREFSY,MZCONF(1),MZCONF(1),
     &      WRK(KCREF),WRK(KCREF),
     &      WRK(KD),WRK(KP),0,1,.FALSE.,.FALSE.,XINDX,
     &      WRK,KFREE,LFREE)
            CALL MEMREL('HSOCTL<-RSPDM',WRK,KD,KCREF,KFREE,LFREE)
            CALL MTRSP(N2ASHX,N2ASHX,WRK(KP1),N2ASHX,WRK(KP2),N2ASHX)
         END IF
         CALL IPSET(0,1,1,0)
      END IF
      OVLAP = D1
      ISYMDN = 1
C
C Configurational gradient
C
      IF ( IGRSYM.EQ.1) THEN
         LCON = KZCONF.GT.1
      ELSE
         LCON = KZCONF.GT.0
      END IF
      LREFST = .FALSE.
C
C Create the gradient
C
      TRPSAV=TRPLET
      CALL HSO2CR(IGRSYM,IGRSPI,GP,WRK(KVEC),
     *            2*IVAR,NCONF,ISYMV,ISYMDN,
     *            XINDX,OVLAP,WRK(KD),WRK(KP),WRK(KFI),
     *            WRK(KFREE),LFREE,
     *            KZYVAR,LCON,LORB,LREFST,LUMHSO,KSYMOP,
     *            IOPSPI,ISPINV,IKLVL,DUM,IDUM,DUM,IDUM,MJWOP)
      CALL HEADER('SPIN-ORBIT GRADIENT FROM HSOAL2',3)
      TRPLET=TRPSAV
C
      IF (IPRRSP.GT.100) THEN
         WRITE(LUPRI,'(/A)') ' Orbital part'
         CALL OUTPUT(GP(1+KZCONF),1,KZWOPT,1,2,KZVAR,2,1,LUPRI)
         WRITE(LUPRI,'(/A)') ' Configuration part'
         CALL OUTPUT(GP,1,KZCONF,1,2,KZVAR,2,1,LUPRI)
      ELSE
         CALL RSPPRW(GP(1+KZCONF),MJWOP,KZWOPT,KSYMOP,KZVAR,LUPRI)
         CALL RSPPRC(GP,KZCONF,KZVAR,LUPRI)
      END IF
      IF (KSYMOP.EQ.1 .AND. LCON) THEN
         WRITE(LUPRI,'(/A)') 'Project out reference component'
         KCREF = KVEC + IVAR
         OVLAP = DDOT(KZCONF,WRK(KCREF),1,GP,1)
         CALL DAXPY(KZCONF,-OVLAP,WRK(KCREF),1,GP,1)
         OVLAP = DDOT(KZCONF,WRK(KCREF),1,GP(1+KZVAR),1)
         CALL DAXPY(KZCONF,-OVLAP,WRK(KCREF),1,GP(1+KZVAR),1)
         CALL RSPPRW(GP(1+KZCONF),MJWOP,KZWOPT,KSYMOP,KZVAR,LUPRI)
         CALL RSPPRC(GP,KZCONF,KZVAR,LUPRI)
      END IF
      CALL MEMREL('HSOAL2',WRK,1,1,KFREE,LFREE)
      CALL QEXIT('HSOAL2')
      RETURN
      END
C
C /* Deck getfmatx */
C
      SUBROUTINE GETFMATX(WORD,TRIPLET,F,D,ND)
#include <implicit.h>
      CHARACTER*8 WORD
      INTEGER ND
      LOGICAL TRIPLET(ND)
#include <inforb.h>
      DIMENSION F(NBAST,NBAST,ND), D(NBAST,NBAST,ND)
C
C Generalize GETFMAT to handle non-symmetric densities (quadratic
C response)
C
#include <priunit.h>
#include <iratdef.h>
C
C LOCAL
C
#if defined (VAR_SPLITFILES)
      PARAMETER (LBUF=600*(IRAT+1)+1)
      PARAMETER (KBUF=1, KIBUF=1 + 600*IRAT)
      DIMENSION IWRK(LBUF)
#endif
      DIMENSION BUF(600)
      INTEGER IBUF(600)
      INTEGER IB, LENGTH, I,J,K,L,ICOOR, LUHSOAO, COMPONENT,LABEL,ID
      DIMENSION LEFT(2), RIGHT(2)
cF90   REAL(DOUBLE) GIJKL,XIJKL,DIJ,DKL,
cF90   REAL(DOUBLE) DP5, D1P5, D1, D2, 
      PARAMETER(DP5=0.5D0, D1P5=1.5D0, D1=1.0D0, D2=2.0D0)
C
      CALL QENTER('GETFMATX')
      IF (ND .GT. 2) THEN
         CALL QUIT('ND .gt. fixed dimension of 2')
      END IF
C
C READ AND DISTRIBUTE SPIN-ORBIT AO INTEGRALS
C
      DO ID=1,ND
         IF (TRIPLET(ID)) THEN
            LEFT(ID)=D1
            RIGHT(ID)=D2
         ELSE
            LEFT(ID)=D2
            RIGHT(ID)=D1
         END IF
      END DO
      LUHSOAO=-1
      CALL GPOPEN(
     &   LUHSOAO,'AO2SOINT','OLD','SEQUENTIAL','UNFORMATTED',-1,.FALSE.
     &   )
      CALL REWSPL(LUHSOAO)
#if defined (VAR_SPLITFILES)
      JRSP = 0
#endif
      CALL MOLLAB('AO2SOINT',LUHSOAO,LUPRI)
      IF (WORD(1:1).EQ.'X') COMPONENT = 1
      IF (WORD(1:1).EQ.'Y') COMPONENT = 2
      IF (WORD(1:1).EQ.'Z') COMPONENT = 3
cF90  READFILE: DO 
 100  CONTINUE
#if defined (VAR_SPLITFILES)
         CALL READSI(LUHSOAO,LBUF,IWRK,JRSP,.TRUE.)
         CALL DCOPY(600,IWRK(KBUF),1,BUF,1)
         CALL ICOPY(600,IWRK(KIBUF),1,IBUF,1)
         LENGTH=IWRK(LBUF)
#else
         READ(LUHSOAO) BUF, IBUF, LENGTH
#endif
cF90     IF (LENGTH.LE.0) EXIT
         IF (LENGTH.LE.0) GO TO 300
cF90     BUFFER: DO IB=1,LENGTH
         DO IB=1,LENGTH
            LABEL=IBUF(IB)
#ifdef VAR_INT64
            I = IBITS(LABEL,48,16)
            J = IBITS(LABEL,32,16)
            K = IBITS(LABEL,16,16)
            L = IBITS(LABEL, 0,16)
#else
            I = IBITS(LABEL,24,8)
            J = IBITS(LABEL,16,8)
            K = IBITS(LABEL, 8,8)
            L = IBITS(LABEL, 0,8)
#endif
            IF (L.EQ.0) THEN
               ICOOR = K
cF90           CYCLE BUFFER
               GO TO 200
            END IF
cF90        IF (ICOOR.NE.COMPONENT) CYCLE BUFFER
            IF (ICOOR.NE.COMPONENT) GO TO 200
            GIJKL=BUF(IB)
            IF (I.EQ.J) GIJKL = DP5*GIJKL
            XIJKL=D1P5*GIJKL
            DO ID=1,ND
               DIJ = LEFT(ID)*(D(I,J,ID) + D(J,I,ID))
               DKL = RIGHT(ID)*(D(L,K,ID) - D(K,L,ID))
               F(I,J,ID) = F(I,J,ID) + GIJKL*DKL
               F(J,I,ID) = F(J,I,ID) + GIJKL*DKL
               F(K,L,ID) = F(K,L,ID) + GIJKL*DIJ
               F(L,K,ID) = F(L,K,ID) - GIJKL*DIJ
               F(I,L,ID) = F(I,L,ID) - XIJKL*D(J,K,ID)
               F(L,I,ID) = F(L,I,ID) + XIJKL*D(K,J,ID)
               F(I,K,ID) = F(I,K,ID) + XIJKL*D(J,L,ID)
               F(K,I,ID) = F(K,I,ID) - XIJKL*D(L,J,ID)
               F(J,L,ID) = F(J,L,ID) - XIJKL*D(I,K,ID)
               F(L,J,ID) = F(L,J,ID) + XIJKL*D(K,I,ID)
               F(J,K,ID) = F(J,K,ID) + XIJKL*D(I,L,ID)
               F(K,J,ID) = F(K,J,ID) - XIJKL*D(L,I,ID)
            END DO
cF90     END DO BUFFER
  200    CONTINUE
         END DO
cF90  END DO READFILE
      GO TO 100
  300 CONTINUE
C
      CALL GPCLOSE(LUHSOAO,'KEEP')
      CALL QEXIT('GETFMATX')
cF90  END SUBROUTINE GETFMATX
      END
