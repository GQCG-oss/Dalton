!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
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
C
C  /* Deck hfcinp */
      SUBROUTINE HFCINP(WORD)
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      LOGICAL NEWDEF
      PARAMETER (NTABLE = 6)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include "esrhfc.h"
      DATA TABLE /'.HFC-FC', '.HFC-SD', '.HFC-SO', '.BRT-SO',
     &            '.EFF-SO', '.PRINT '/
      LOGICAL HFCOLD
C
      NEWDEF = WORD .EQ. '*HFC   '
      ICHANG = 0
C
      IF (NEWDEF) THEN
         WORD1 = WORD
 100          CONTINUE
         READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
            GOTO 100
         ELSE IF (PROMPT .EQ. '.') THEN
            ICHANG = ICHANG + 1
            DO 200 I = 1, NTABLE
               IF (TABLE(I) .EQ. WORD) THEN
                  GOTO (1,2,3,4,5,6), I
               END IF
 200                   CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               GOTO 100
            END IF
            WRITE (LUPRI,'(/,3A,/)') 'Keyword "', WORD,
     &           '" not recognized in HFC.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in HFC')
 1          CONTINUE
               HFCFC  = .TRUE.
               GOTO 100
 2          CONTINUE
               HFCSD  = .TRUE.
               GOTO 100
 3          CONTINUE
               HFCSO  = .TRUE.
               GOTO 100
 4          CONTINUE
               HSOBRT = .TRUE.
               GOTO 100
 5          CONTINUE
               HSOEFF = .TRUE.
               GOTO 100
 6          CONTINUE
               READ (LUCMD,*) IPRHFC
               GOTO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            GOTO 300
         ELSE
            WRITE (LUPRI,'(/,3A,/)') 'Prompt "', WORD,
     &           '" not recognized in HFC.'
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords', LUPRI)
         CALL QUIT('Illegal prompt in HFC')
      END IF
 300    CONTINUE
      IF (HSOBRT .AND. HSOEFF) THEN
         CALL QUIT('Inconsistent input in HFC only one type of '//
     &              'the spin-orbit operator can be selected.')
      END IF
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults in HFC:',0)
         IF (IPRHFC .NE. 0) WRITE(LUPRI,'(A,I5)')
     &       ' Print level in HFC:', IPRHFC
         IF (HFCFC) WRITE (LUPRI,'(A)')
     &                 ' Isotropic Fermi contact contribution '//
     &                 'will be computed'
         IF (HFCSD) WRITE (LUPRI,'(A)')
     &                 ' Traceless spin-dipolar contribution '//
     &                 'will be computed'
         IF (HFCSO) THEN
            WRITE (LUPRI,'(A)') ' Spin-orbit contribution '//
     &             'will be computed:'
            IF (HSOBRT) WRITE (LUPRI,'(A)')
     &         '  using full Breit-Pauli spin-orbit operator'
            IF (HSOEFF) WRITE (LUPRI,'(A)')
     &         '  using effective charge spin-orbit operator'
         END IF
      END IF
C
      ESRHFC = .TRUE.
C
      RETURN
      END
C  /* Deck hfcdrv */
      SUBROUTINE HFCDRV(CMO,UDV,PVX,FOCK,FC,FV,FCAC,H2AC,XINDX,
     *                WORK,LWORK)
C
C     Driver for evaluation of hyperfine coupling constants
C
#include "implicit.h"
C
      DIMENSION CMO(*), UDV(*), PVX(*), FOCK(*),FC(*), FV(*),
     *          FCAC(*), H2AC(*), XINDX(*), WORK(*)
#include "iratdef.h"
#include "maxorb.h"
#include "priunit.h"
#include "inforb.h"
#include "infpri.h"
#include "infmp2.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "infdim.h"
#include "inftra.h"
#include "inftap.h"
#include "rspprp.h"
#include "gnrinf.h"
#include "infinp.h"
#include "inflr.h"
#include "mxcent.h"
#include "esrhfc.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "symmet.h"
#include "chrnos.h"
#include "chrxyz.h"

#include "codata.h"
#include "gfac.h"
#include "cbiher.h"
C
      LOGICAL TRPSAVE, OLDDX, FNDLB2
      EXTERNAL FNDLB2
      CHARACTER*7 RTNLBL(2)
      CHARACTER LABINT*8
C
      CALL QENTER('HFCDRV')
C
      CALL DZERO(AFCEXP,MXCENT)
      CALL DZERO(AFCPOL,MXCENT)
      CALL DZERO(ASDEXP,9*MXCENT)
      CALL DZERO(ASDPOL,9*MXCENT)
      CALL DZERO(ASOVAL,9*MXCENT)
C
      THCRSP = 1.0D-7
      THCESR = 1.0D-7
C
      KFREE = 1
      LFREE = LWORK

C    Compute FC,SD expectation values

      IF (HFCFC) THEN
         CALL AFCAVE(CMO,XINDX,WORK(KFREE),LFREE)
      END IF
      IF (HFCSD) THEN
         CALL ASDAVE(CMO,XINDX,WORK(KFREE),LFREE)
      END IF
C
C    Density relaxation contribution to FC, SD
C
C
      IF (HFCFC .OR. HFCSD) THEN
         TRPSAVE = TRPLET
         TRPLET = .TRUE.
         KSYMOP = 1
         CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WORK(KFREE),
     &               LFREE)
         CALL RSPLAN(CMO,UDV,PVX,FC,FV,FCAC,H2AC,XINDX,
     &               WORK(KFREE),LFREE)
         KLAGR = KFREE
         KGP   = KLAGR + KZYVAR
         KFREE = KGP   + KZYVAR
         LFREE = LFREE - 2*KZYVAR
C        Spin polarization correction to Fermi contact contribution
         IF (HFCFC) THEN
            DO IDX = 1,NUCIND
               LABINT = 'FC '//NAMDEP(NUCPRE(IDX)+1)(1:3)
     &                   //CHRNOS(IPTNUC(IDX,0)/10)//
     &                   CHRNOS(MOD(IPTNUC(IDX,0),10))
               CALL GETGPV(LABINT,FC,FV,CMO,UDV,PVX,XINDX,
     &                     ANTSYM,WORK(KGP),LFREE)
               AFCPOL(IDX) = -DDOT(KZYVAR,WORK(KLAGR),1,WORK(KGP),1)
            END DO
         END IF
C        Spin polarization correction to Spin-dipolar contribution
         IF (HFCSD) THEN
            DO IDX = 1,NUCIND
               DO IREP = 0,MAXREP
                  DO ICOOR1 = 1, 3
                     ISCOR1 = IPTCNT(3*(IDX-1)+ICOOR1,IREP,2)
                     IF (ISCOR1 .GT. 0) THEN
                        DO ICOOR2 = 1, 3
                            IDOCOMP = IEOR(IREP,ISYMAX(ICOOR2,2))
                            IF (IDOCOMP .EQ. 0) THEN
                               IFIRST = ISCOR1/100
                               ISECND = MOD(ISCOR1,100)/10
                               ITHIRD = MOD(MOD(ISCOR1,100),10)
                               LABINT = 'SD '//CHRNOS(IFIRST)
     &                                 //CHRNOS(ISECND)
     &                                 //CHRNOS(ITHIRD)//' '
     &                                 //CHRXYZ(-ICOOR2)
                               CALL GETGPV(LABINT,FC,FV,CMO,UDV,PVX,
     &                                     XINDX,ANTSYM,WORK(KGP),LFREE)
                               ASDPOL(IDX,ICOOR1,ICOOR2) = -DDOT(KZYVAR,
     &                                     WORK(KLAGR),1,WORK(KGP),1)
                            END IF
                        END DO
                     END IF
                  END DO
               END DO
            END DO
         END IF
         TRPLET = TRPSAVE
      END IF
C
C     Spin-orbit contribution to hyperfine coupling tensor
C
      IF (HFCSO) THEN
C     get PSO and AMFI SO integrals
         IF (LUPROP .LE. 0) THEN
            CALL GPOPEN(LUPROP,'AOPROPER','OLD','SEQUENTIAL',
     &                  'UNFORMATTED',IDUMMY,OLDDX)
         END IF
         REWIND LUPROP
         IF (.NOT.FNDLB2('PSO    ',RTNLBL,LUPROP)) THEN
            CALL PR1INT('PSO    ',WORK(KFREE),LFREE,IDUMMY,
     &                  IDUMMY,.FALSE.,.FALSE.,0)
         END IF
         CALL GPCLOSE(LUPROP,'KEEP')
         IF (.NOT. (HSOBRT .OR. HSOEFF)) MNF_SO = .TRUE.
C        solve response equations and compute values
         TRPSAVE = TRPLET
         TRPLET = .FALSE.
         TRPFLG = .FALSE.
         IF (HSOBRT) THEN
            LLROP(INDPRP('X SPNORB'))  = .TRUE.
            LLROP(INDPRP('Y SPNORB'))  = .TRUE.
            LLROP(INDPRP('Z SPNORB'))  = .TRUE.
         ELSE IF (HSOEFF) THEN
            LLROP(INDPRP('X1SPNSCA'))  = .TRUE.
            LLROP(INDPRP('Y1SPNSCA'))  = .TRUE.
            LLROP(INDPRP('Z1SPNSCA'))  = .TRUE.
         ELSE
            LLROP(INDPRP('X1MNF-SO'))  = .TRUE.
            LLROP(INDPRP('Y1MNF-SO'))  = .TRUE.
            LLROP(INDPRP('Z1MNF-SO'))  = .TRUE.
         END IF
         DO 100 ISYM=1,NSYM
            KSYMOP = ISYM
            CALL RSPSYM(WORK(KFREE),LFREE)
            CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,
     &                  WORK(KFREE),LFREE)
            IF (KZVAR.EQ.0) GO TO 100
            IF (NGPLR(ISYM).GT.0) THEN
               CALL RSPLR(CMO,UDV,PVX,FC,FV,FCAC,H2AC,XINDX,
     &                    WORK(KFREE),LFREE)
               CALL ASOCOMP(CMO,UDV,PVX,XINDX,WORK(KFREE),LFREE)
            END IF
 100             CONTINUE
         CALL IZERO(NGPLR,8)
         TRPFLG = .TRUE.
         TRPLET = TRPSAVE
      END IF
C     Convert to MHz
      CALL HFCCON(WORK(KFREE),LFREE)
C     Print results
      CALL HFCRES(WORK(KFREE),LFREE)
C
      CALL QEXIT('HFCDRV')
      RETURN
      END
C  /* Deck afcave */
      SUBROUTINE AFCAVE(CMO, XINDX, WORK, LWORK)
#include "implicit.h"
      DIMENSION CMO(*), XINDX(*), WORK(*)
#include "dummy.h"
#include "inftap.h"
#include "inforb.h"
#include "priunit.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "maxorb.h"
#include "infinp.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "mxcent.h"
#include "nuclei.h"
#include "symmet.h"
#include "chrnos.h"
#include "esrhfc.h"
      LOGICAL OLDDX, FNDLB2, TRPSAVE
      EXTERNAL FNDLB2
      CHARACTER*7 RTNLBL(2)
      CHARACTER LABINT*8
C
      CALL QENTER('AFCAVE')
C     Get Fermi contact integrals (AO basis)
      IF (LUPROP .LE. 0) THEN
         CALL GPOPEN(LUPROP,'AOPROPER','OLD','SEQUENTIAL','UNFORMATTED',
     &               IDUMMY,OLDDX)
      END IF
      REWIND LUPROP
      IF (.NOT.FNDLB2('FERMI C',RTNLBL,LUPROP)) THEN
         CALL PR1INT('FERMI C',WORK,LWORK,IDUMMY,
     &               IDUMMY,.FALSE.,.FALSE.,0)
      END IF
      CALL GPCLOSE(LUPROP,'KEEP')
C
      KFREE=1
      LFREE=LWORK
      CALL MEMGET('REAL',KPRPMO,N2ORBX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDV,N2ASHX,WORK,KFREE,LFREE)
C     Get triplet density
      IF (TDHF) THEN
         CALL DUNIT(WORK(KDV),NASHT)
      ELSE
         CALL MEMGET('REAL',KCREF,NCREF,WORK,KFREE,LFREE)
         CALL GETREF(WORK(KCREF),NCREF)
         CALL RSPDM(IREFSY,IREFSY,NCREF,NCREF,WORK(KCREF),WORK(KCREF),
     &              WORK(KDV),DUMMY,1,0,.FALSE.,.TRUE.,XINDX,WORK,
     &              KFREE,LFREE)
      END IF
C     Compute expectation values
      KSYMOP=1
      TRPSAVE = TRPLET
      TRPLET = .TRUE.
      DO IK = 1, NUCIND
         LABINT = 'FC '//NAMDEP(NUCPRE(IK)+1)(1:3)
     &                //CHRNOS(IPTNUC(IK,0)/10)//
     &                  CHRNOS(MOD(IPTNUC(IK,0),10))
         CALL PRPGET(LABINT,CMO,WORK(KPRPMO),KSYMOP,ANTSYM,WORK(KFREE),
     &               LFREE,0)
         CALL PRPONE(WORK(KPRPMO),WORK(KDV),AFCEXP(IK),1,LABINT)
      END DO
      TRPLET = TRPSAVE
C
      CALL MEMREL('AFCAVE',WORK,1,1,KFREE,LFREE)
C
      CALL QEXIT('AFCAVE')
      RETURN
      END
C  /* Deck asdave */
      SUBROUTINE ASDAVE(CMO, XINDX, WORK, LWORK)
#include "implicit.h"
      DIMENSION CMO(*), XINDX(*), WORK(*)
#include "dummy.h"
#include "inftap.h"
#include "inforb.h"
#include "priunit.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "maxorb.h"
#include "infinp.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "mxcent.h"
#include "nuclei.h"
#include "symmet.h"
#include "chrnos.h"
#include "chrxyz.h"
#include "esrhfc.h"

      LOGICAL OLDDX, FNDLB2, TRPSAVE
      EXTERNAL FNDLB2
      CHARACTER*7 RTNLBL(2)
      CHARACTER LABINT*8
C
      CALL QENTER('ASDAVE')
C     Get Spin-dipolar integrals (AO basis)
      IF (LUPROP .LE. 0) THEN
         CALL GPOPEN(LUPROP,'AOPROPER','OLD','SEQUENTIAL','UNFORMATTED',
     &               IDUMMY,OLDDX)
      END IF
      REWIND LUPROP
      IF (.NOT.FNDLB2('SPIN-DI',RTNLBL,LUPROP)) THEN
         CALL PR1INT('SPIN-DI',WORK,LWORK,IDUMMY,
     &               IDUMMY,.FALSE.,.FALSE.,0)
      END IF
      CALL GPCLOSE(LUPROP,'KEEP')
C
      KFREE=1
      LFREE=LWORK
      CALL MEMGET('REAL',KPRPMO,N2ORBX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDV,N2ASHX,WORK,KFREE,LFREE)
C     Get triplet density
      IF (TDHF) THEN
         CALL DUNIT(WORK(KDV),NASHT)
      ELSE
         CALL MEMGET('REAL',KCREF,NCREF,WORK,KFREE,LFREE)
         CALL GETREF(WORK(KCREF),NCREF)
         CALL RSPDM(IREFSY,IREFSY,NCREF,NCREF,WORK(KCREF),WORK(KCREF),
     &              WORK(KDV),DUMMY,1,0,.FALSE.,.TRUE.,XINDX,WORK,
     &              KFREE,LFREE)
      END IF
C     Compute expectation values
      KSYMOP=1
      TRPSAVE = TRPLET
      TRPLET = .TRUE.
      DO IDX = 1,NUCIND
         DO IREP = 0,MAXREP
            DO ICOOR1 = 1, 3
               ISCOR1 = IPTCNT(3*(IDX-1)+ICOOR1,IREP,2)
               IF (ISCOR1 .GT. 0) THEN
                  DO ICOOR2 = 1, 3
                     IF (IEOR(IREP,ISYMAX(ICOOR2,2)) .EQ. 0) THEN
                        IFIRST = ISCOR1/100
                        ISECND = MOD(ISCOR1,100)/10
                        ITHIRD = MOD(MOD(ISCOR1,100),10)
                        LABINT = 'SD '//CHRNOS(IFIRST)
     &                           //CHRNOS(ISECND)
     &                           //CHRNOS(ITHIRD)//' '
     &                           //CHRXYZ(-ICOOR2)
                        CALL PRPGET(LABINT,CMO,WORK(KPRPMO),KSYMOP,
     &                              ANTSYM,WORK(KFREE),LFREE,0)
                        CALL PRPONE(WORK(KPRPMO),WORK(KDV),
     &                              ASDEXP(IDX,ICOOR1,ICOOR2),1,LABINT)
                     END IF
                  END DO
               END IF
            END DO
         END DO
      END DO
      TRPLET = TRPSAVE
C
      CALL MEMREL('ASDAVE',WORK,1,1,KFREE,LFREE)
C
      CALL QEXIT('ASDAVE')
      RETURN
      END
C Deck /* asocomp */
      SUBROUTINE ASOCOMP(CMO,UDV,PV,XINDX,WORK,LWORK)
#include"implicit.h"
      DOUBLE PRECISION CMO(*), UDV(*), PV(*), XINDX(*), WORK(*)
      INTEGER LWORK
#include "qrinf.h"
#include "inforb.h"
#include "rspprp.h"
#include "inflr.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "symmet.h"
#include "wrkrsp.h"
#include "nuclei.h"
#include "chrnos.h"
#include "esrhfc.h"
#include "priunit.h"
C
      INTEGER KFREE, LFREE
      LOGICAL LDUMMY, FOUND, CONV, OPENED
      INTEGER I, J, KLEN, LU, KSO, ISYM, IPRP, ISO
      DOUBLE PRECISION RESID, ANTSYM, D0, DNORM
      PARAMETER (D0=0.0D0)
      CHARACTER*8 BLANK, LABINT, SOLAB(3)
      DATA BLANK/'        '/
C
      CALL QENTER('ASOCOMP')
C
      KFREE=1
      LFREE=LWORK
C
      INQUIRE(FILE='RSPVEC',OPENED=OPENED,NUMBER=LU)
      IF (.NOT. OPENED) THEN
         LU=-1
         CALL GPOPEN(LU,'RSPVEC','OLD','SEQUENTIAL','UNFORMATTED',0,
     &               LDUMMY)
      END IF
      REWIND(LU)
C
      IF (HSOBRT) THEN
         SOLAB(1) = 'X SPNORB'
         SOLAB(2) = 'Y SPNORB'
         SOLAB(3) = 'Z SPNORB'
      ELSE IF (HSOEFF) THEN
         SOLAB(1) = 'X1SPNSCA'
         SOLAB(2) = 'Y1SPNSCA'
         SOLAB(3) = 'Z1SPNSCA'
      ELSE
         SOLAB(1) = 'X1MNF-SO'
         SOLAB(2) = 'Y1MNF-SO'
         SOLAB(3) = 'Z1MNF-SO'
      END IF
C
      CALL MEMGET('REAL',KSO,KZYVAR,WORK,KFREE,LFREE)
      DO ISN = 1,3
         DO IPRP=1,NGPLR(KSYMOP)
            IF (LBLLR(KSYMOP,IPRP) .EQ. SOLAB(ISN)) THEN
               CALL REARSP(LU,KLEN,WORK(KSO),SOLAB(ISN),BLANK,D0,DUMMY,
     &                     KSYMOP, IDUMMY,THCLR,FOUND,CONV,ANTSYM)
               IF (KLEN.NE.KZYVAR) THEN
                  CALL QUIT('ASOCOMP:error while reading RSPVEC')
               END IF
               DNORM = DNRM2(KLEN,WORK(KSO),1)
               IF (FOUND .AND. CONV .AND. DNORM.GT.THRNRM) THEN
                  DO IDX = 1,NUCIND
                     DO ICOORX = 1,3
                        ISCORX = IPTCNT(3*(IDX-1)+ICOORX,KSYMOP-1,2)
                        IF (ISCORX .GT. 0) THEN
                           IFIRST = ISCORX/100
                           ISECND = MOD(ISCORX,100)/10
                           ITHIRD = MOD(MOD(ISCORX,100),10)
                           LABINT = 'PSO '//CHRNOS(IFIRST)//
     &                               CHRNOS(ISECND)//
     &                               CHRNOS(ITHIRD)//' '
                           CALL GETGPV(LABINT,DUMMY,DUMMY,CMO,UDV,PV,
     &                                 XINDX,ANTSYM,WORK(KFREE),LFREE)
                           ASOVAL(IDX,ICOORX,ISN) = DDOT(KLEN,WORK(KSO),
     &                                                 1,WORK(KFREE),1)
                        END IF
                     END DO
                  END DO
               END IF
            END IF
         END DO
      END DO

      IF (.NOT.OPENED) CALL GPCLOSE(LU,'KEEP')
      CALL MEMREL('ASOCOMP',WORK,KSO,KSO,KFREE,LFREE)
      CALL QEXIT('ASOCOMP')
      RETURN
      END
C  /* Deck hfccon */
      SUBROUTINE HFCCON(WORK,LWORK)
#include "implicit.h"
      DIMENSION WORK(LWORK)
#include "iratdef.h"
#include "maxorb.h"
#include "infinp.h"
#include "mxcent.h"
#include "esrhfc.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "codata.h"
#include "gfac.h"
      PARAMETER ( D1M = -1.0D0, HFCTHR = 1.0D-8)
C
      REFSPIN = DBLE(ISPIN-1)
      RFACTOR = (XTHZ*1.0D-6*ALPHA2)/(2.0D0*XPRTMAS*REFSPIN)
C     Conversion to Gauss
      IF (UNGAUSS) THEN
         RFACTOR = RFACTOR/2.802494
      END IF
C

      DO IDX = 1,NUCDEP
         GNVAL = DISOTP(NINT(CHARGE(IDX)),ISOTOP(IDX),'GVAL')
         RMULT = GNVAL*RFACTOR/NUCDEG(IDX)
         IF (HFCFC) THEN
             AFCEXP(IDX) = RMULT*AFCEXP(IDX)
             AFCPOL(IDX) = RMULT*AFCPOL(IDX)
             IF (DABS(AFCEXP(IDX)) .LE. HFCTHR)
     &          AFCEXP(IDX) = 0.0D0
             IF (DABS(AFCPOL(IDX)) .LE. HFCTHR)
     &          AFCPOL(IDX) = 0.0D0
         END IF
         IF (HFCSD) THEN
            DO I = 1,3
               DO  J = 1,3
                   ASDEXP(IDX,I,J) = RMULT*ASDEXP(IDX,I,J)
                   ASDPOL(IDX,I,J) = RMULT*ASDPOL(IDX,I,J)
                   IF (DABS(ASDEXP(IDX,I,J)) .LE. HFCTHR)
     &                ASDEXP(IDX,I,J) = 0.0D0
                   IF (DABS(ASDPOL(IDX,I,J)) .LE. HFCTHR)
     &                ASDPOL(IDX,I,J) = 0.0D0
               END DO
            END DO
         END IF
         IF (HFCSO) THEN
            RMULT = D1M*GFAC*RMULT
            DO I = 1,3
               DO  J = 1,3
                   ASOVAL(IDX,I,J) = RMULT*ASOVAL(IDX,I,J)
                   IF (DABS(ASOVAL(IDX,I,J)) .LE. HFCTHR)
     &                ASOVAL(IDX,I,J) = 0.0D0
               END DO
            END DO
         END IF
      END DO
C
      RETURN
      END
C  /* Deck hfcres */
      SUBROUTINE HFCRES(WORK,LWORK)
#include "implicit.h"
      DIMENSION WORK(LWORK)
#include "priunit.h"
#include "mxcent.h"
#include "cbiwlk.h"
#include "esrhfc.h"
C
      CALL TITLER('RESPONSE - HYPERFINE COUPLING CONSTANTS','*',116)
C     contributions  printing
      IPRHFC = 11
C     for testing only
      IF (IPRHFC .GE. 10) THEN
         IF (HFCFC) THEN
            IF (UNGAUSS) THEN
               CALL AROUND('Fermi contact contribution to HFC (in G)')
            ELSE
               CALL AROUND('Fermi contact contribution to HFC (in MHz)')
            END IF
            CALL AFCPRT(AFCEXP,AFCPOL)
         END IF
         IF (HFCSD) THEN
            IF (UNGAUSS) THEN
               CALL AROUND('Spin-dipolar contribution to HFC (in G )')
            ELSE
               CALL AROUND('Spin-dipolar contribution to HFC (in MHz)')
            END IF
            CALL ASDPRT(ASDEXP,ASDPOL)
         END IF
         IF (HFCSO) THEN
            IF (UNGAUSS) THEN
                CALL AROUND('Spin-orbit contribution to HFC (in G)')
            ELSE
                CALL AROUND('Spin-orbit contribution to HFC (in MHz)')
            END IF
            CALL ASOPRT(ASOVAL)
         END IF
      END IF
C     total hyperfine coupling constants
      IF (UNGAUSS) THEN
         CALL AROUND('Symmetrized  hyperfine coupling tensor (in G)')
      ELSE
         CALL AROUND('Symmetrized  hyperfine coupling tensor (in MHz)')
      END IF
      CALL ATOTPRT(AFCEXP,AFCPOL,ASDEXP,ASDPOL,ASOVAL,WORK,LWORK)
C
      RETURN
      END
C  /* Deck afcprt */
      SUBROUTINE AFCPRT(XEXP,XPOL)
#include "implicit.h"
      PARAMETER(GNTHR = 1.0D-5)
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
C
      DIMENSION XEXP(*), XPOL(*)
C
      WRITE (LUPRI,'(2X,A,/,2X,A,/,2X,A)')
     & '-------------------------------------------------',
     & '| Atom | Isot. |  <FC>    | <<FC;H>> |  Total   |',
     & '-------------------------------------------------'
      DO IDX = 1, NUCIND
          GNVAL=DABS(DISOTP(NINT(CHARGE(IDX)),ISOTOP(IDX),'GVAL'))
          IF (GNVAL .GT. GNTHR) THEN
             INUCMAS = NINT(DISOTP(NINT(CHARGE(IDX)),ISOTOP(IDX),
     &                      'MASS'))
             WRITE(LUPRI,'(2X,A,I3,A,I3,A,F8.3,A,F8.3,A,F8.3,A)')
     &              '|',IDX,NAMDEP(NUCPRE(IDX)+1)(1:3)//'|  ', INUCMAS,
     &               '  | ', XEXP(IDX),' | ', XPOL(IDX), ' | ',
     &               XEXP(IDX)+XPOL(IDX), ' | '
          END IF
      END DO
      WRITE (LUPRI,'(2X,A,/)')
     & '-------------------------------------------------'
C
      RETURN
      END
      SUBROUTINE ASDPRT(XEXP,XPOL)
#include "implicit.h"
      PARAMETER(GNTHR = 1.0D-5)
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
      DIMENSION XEXP(MXCENT,3,3), XPOL(MXCENT,3,3)
C
        WRITE (LUPRI,'(/,2X,A,/,2X,A,/,2X,A)')
     & '------------------------------------------------------------'//
     & '--------------------------------',
     & '| Atom | Isot. | Contrib. |  A_xx   |  A_yy    |  A_zz    '//
     & '|  A_xy    |   A_xz   |  A_yz    |',
     & '------------------------------------------------------------'//
     & '--------------------------------'
      DO IDX = 1, NUCIND
          GNVAL=DABS(DISOTP(NINT(CHARGE(IDX)),ISOTOP(IDX),'GVAL'))
          IF (GNVAL .GT. GNTHR) THEN
             INUCMAS = NINT(DISOTP(NINT(CHARGE(IDX)),ISOTOP(IDX),
     &                      'MASS'))
             WRITE(LUPRI,'(2X,A,I3,A,I3,A,F8.2,A,F8.2,A,F8.2,A,F8.2,'//
     &              'A,F8.2,A,F8.2,A)')
     &              '|',IDX,NAMDEP(NUCPRE(IDX)+1)(1:3)//'|  ', INUCMAS,
     &               '  |   <SD>   |', XEXP(IDX,1,1),' | ',
     &               XEXP(IDX,2,2), ' | ', XEXP(IDX,3,3), ' | ',
     &               XEXP(IDX,1,2), ' | ', XEXP(IDX,1,3), ' | ',
     &               XEXP(IDX,2,3), ' |'
              WRITE(LUPRI,'(2X,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,'//
     &              'F8.2,A)') '|--------------| <<SD;H>> |',
     &               XPOL(IDX,1,1),' | ',XPOL(IDX,2,2), ' | ',
     &               XPOL(IDX,3,3), ' | ',XPOL(IDX,1,2), ' | ',
     &               XPOL(IDX,1,3), ' | ', XPOL(IDX,2,3), ' |'
              WRITE(LUPRI,'(2X,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,'//
     &              'F8.2,A)') '|              |   Total  |',
     &               XPOL(IDX,1,1)+XEXP(IDX,1,1),' | ',
     &               XPOL(IDX,2,2)+XEXP(IDX,2,2), ' | ',
     &               XPOL(IDX,3,3)+XEXP(IDX,3,3), ' | ',
     &               XPOL(IDX,1,2)+XEXP(IDX,1,2), ' | ',
     &               XPOL(IDX,1,3)+XEXP(IDX,1,3), ' | ',
     &               XPOL(IDX,2,3)+XEXP(IDX,2,3), ' |'
              WRITE (LUPRI,'(2X,A)')
     &        '-----------------------------------------------------'//
     &        '---------------------------------------'
          END IF
      END DO
C
      RETURN
      END
C  /* Deck asoprt */
      SUBROUTINE ASOPRT(XVAL)
#include "implicit.h"
      PARAMETER(GNTHR = 1.0D-5)
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
      DIMENSION  XVAL(MXCENT,3,3)
C
      DO IDX = 1,NUCIND
         GNVAL=DABS(DISOTP(NINT(CHARGE(IDX)),ISOTOP(IDX),'GVAL'))
         IF (GNVAL .GT. GNTHR) THEN
            INUCMAS = NINT(DISOTP(NINT(CHARGE(IDX)),ISOTOP(IDX),
     &                     'MASS'))
            WRITE(LUPRI,'(2X,A,I3,A,I3,A,/,2X,A)') 'Atom: ', IDX,
     &            NAMDEP(NUCPRE(IDX)+1)(1:3)//'Isot.:', INUCMAS,
     &            ' <<PSO;SO>>:',
     &            '================================='
            WRITE(LUPRI,'(2X,A)') '-----------------------------'//
     &      '---------'
            WRITE(LUPRI,'(2X,A,3X,A,3X,A,3X,A,3X,A,3X,A,3X,A)')
     &            '|    |','Sx',' | ', 'Sy', ' | ','Sz',' |'
            WRITE(LUPRI,'(2X,A)') '-----------------------------'//
     &      '---------'
            WRITE(LUPRI,'(2X,A,F8.2,A,F8.2,A,F8.2,A)')
     &            '| Ix |',XVAL(IDX,1,1),' | ', XVAL(IDX,1,2),' | ',
     &             XVAL(IDX,1,3),' |'
            WRITE(LUPRI,'(2X,A,F8.2,A,F8.2,A,F8.2,A)')
     &            '| Iy |', XVAL(IDX,2,1),' | ',XVAL(IDX,2,2),' | ',
     &            XVAL(IDX,2,3),' |'
            WRITE(LUPRI,'(2X,A,F8.2,A,F8.2,A,F8.2,A)')
     &            '| Iz |',XVAL(IDX,3,1),' | ',XVAL(IDX,3,2),
     &            ' | ',XVAL(IDX,3,3),' |'
            WRITE(LUPRI,'(2X,A,/)') '-----------------------------'//
     &      '---------'
          END IF
      END DO
C
      RETURN
      END
C  /* Deck atotprt */
      SUBROUTINE ATOTPRT(XFCEXP,XFCPOL,XSDEXP,XSDPOL,XSOVAL,WORK,LWORK)
#include "implicit.h"
      DIMENSION  WORK(*)
      PARAMETER(GNTHR = 1.0D-5, DP5 = 0.5D0, D3 = 3.0D0)
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
      DIMENSION XFCEXP(MXCENT), XFCPOL(MXCENT), XSDPEXP(MXCENT,3,3)
      DIMENSION XSDPOL(MXCENT,3,3), XSOVAL(MXCENT,3,3)
C
      DOUBLE PRECISION ATOTAL, ATOTVAL, PVALVEC
      DIMENSION ATOTAL(MXCENT,3,3), ATOTVAL(MXCENT,3)
      DIMENSION PVALVEC(MXCENT,3,3)
C
      CALL AHFCSUM(ATOTAL,XFCEXP,XFCPOL,XSDEXP,XSDPOL,XSOVAL)
      CALL HFCPVAL(ATOTAL,ATOTVAL,PVALVEC,WORK,LWORK)
      DO IDX = 1,NUCIND
          GNVAL=DABS(DISOTP(NINT(CHARGE(IDX)),ISOTOP(IDX),'GVAL'))
          IF (GNVAL .GT. GNTHR) THEN
             INUCMAS = NINT(DISOTP(NINT(CHARGE(IDX)),ISOTOP(IDX),
     &                     'MASS'))
             WRITE(LUPRI,'(2X,A,I3,A,I3,A,/,2X,A)') 'Atom: ', IDX,
     &            NAMDEP(NUCPRE(IDX)+1)(1:3)//' Isot.: ', INUCMAS,
     &            ' Sym. HFC tensor:',
     &            '========================================'
             WRITE(LUPRI,'(2X,A)') '-----------------------------'//
     &       '---------------------------------------------------'//
     &       '--------'
             WRITE(LUPRI,'(2X,A,3X,A,3X,A,3X,A,3X,A,3X,A,3X,A,A)')
     &            '|    |','Sx',' | ', 'Sy', ' | ','Sz',' |',
     &            '  Principal values |       Principal axes        |'
             WRITE(LUPRI,'(2X,A)') '-----------------------------'//
     &       '---------------------------------------------------'//
     &       '--------'
             WRITE(LUPRI,'(2X,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,F5.2,A,'//
     &             'F5.2,A,F5.2,A)')
     &            '| Ix |',ATOTAL(IDX,1,1),' | ', ATOTAL(IDX,1,2),' | ',
     &             ATOTAL(IDX,1,3),' | A_(11) | ', ATOTVAL(IDX,1),
     &             ' | (11) | ', PVALVEC(IDX,1,1),' x ',
     &             PVALVEC(IDX,2,1), ' y ', PVALVEC(IDX,3,1), ' z |'
             WRITE(LUPRI,'(2X,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,F5.2,A,'//
     &             'F5.2,A,F5.2,A)')
     &            '| Ix |',ATOTAL(IDX,2,1),' | ', ATOTAL(IDX,2,2),' | ',
     &             ATOTAL(IDX,2,3),' | A_(22) | ', ATOTVAL(IDX,2),
     &             ' | (22) | ', PVALVEC(IDX,1,2),' x ',
     &             PVALVEC(IDX,2,2), ' y ',PVALVEC(IDX,3,2), ' z |'
             WRITE(LUPRI,'(2X,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,F5.2,A,'//
     &             'F5.2,A,F5.2,A)')
     &            '| Ix |',ATOTAL(IDX,3,1),' | ', ATOTAL(IDX,3,2),' | ',
     &             ATOTAL(IDX,3,3),' | A_(11) | ', ATOTVAL(IDX,3),
     &             ' | (11) | ', PVALVEC(IDX,1,3),' x ',
     &             PVALVEC(IDX,2,3), ' y ',PVALVEC(IDX,3,3), ' z |'
             WRITE(LUPRI,'(2X,A,/)') '-----------------------------'//
     &       '---------------------------------------------------'//
     &       '--------'
          END IF
      END DO
C     Print Summary
      CALL AROUND('Summary of hyperfine coupling constants')
      WRITE (LUPRI,'(2X,A,/,2X,A,/,2X,A)')
     & '-------------------------------------------------------------',
     & '| Atom  | Isot. |  A_(11)  |  A_(22)  |  A_(33)  |  A_iso   |',
     & '-------------------------------------------------------------'
      DO IDX =1, NUCIND
          GNVAL=DABS(DISOTP(NINT(CHARGE(IDX)),ISOTOP(IDX),'GVAL'))
          IF (GNVAL .GT. GNTHR) THEN
             INUCMAS = NINT(DISOTP(NINT(CHARGE(IDX)),ISOTOP(IDX),
     &                      'MASS'))
             WRITE(LUPRI,'(2X,A,I3,A,I3,A,F8.2,A,F8.2,A,F8.2,A,'//
     &             'F8.2,A)') '|',IDX,NAMDEP(NUCPRE(IDX)+1)(1:3)//
     &             ' |  ', INUCMAS,'  | ', ATOTVAL(IDX,1),' | ',
     &             ATOTVAL(IDX,2), ' | ', ATOTVAL(IDX,3), ' | ',
     &             (ATOTVAL(IDX,1)+ATOTVAL(IDX,2)+ATOTVAL(IDX,3))/D3,
     &             ' |'
          END IF
      END DO
      WRITE (LUPRI,'(2X,A,/)')
     & '------------------------------------------------------------'
C
      RETURN
      END
C  /* Deck ahfcsum */
      SUBROUTINE AHFCSUM(ATOT,AEXP,APOL,BEXP,BPOL,CVAL)
#include "implicit.h"
#include "inforb.h"
#include "maxorb.h"
#include "infinp.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "esrhfc.h"
      DIMENSION ATOT(MXCENT,3,3), AEXP(MXCENT), APOL(MXCENT)
      DIMENSION BEXP(MXCENT,3,3), BPOL(MXCENT,3,3), CVAL(MXCENT,3,3)
C
      CALL DZERO(ATOT,9*MXCENT)
      DO IDX = 1,NUCIND
         IF (HFCFC) THEN
            DO I = 1,3
               ATOT(IDX,I,I) = AEXP(IDX) + APOL(IDX)
            END DO
         END IF
         IF (HFCSD) THEN
            DO I = 1,3
               DO J = 1,3
                  ATOT(IDX,I,J) = ATOT(IDX,I,J) + BEXP(IDX,I,J)
     &                          + BPOL(IDX,I,J)
               END DO
            END DO
          END IF
          IF (HFCSO) THEN
             DO I = 1,3
                DO J= 1,3
                   ATOT(IDX,I,J) =  ATOT(IDX,I,J) + CVAL(IDX,I,J)
                END DO
             END DO
          END IF
       END DO
C
       RETURN
       END
C  /* Deck hfcpval */
      SUBROUTINE HFCPVAL(ATEN,APVAL,PAXIS,WORK,LWORK)
#include "implicit.h"
       PARAMETER(DP5 = 0.5D0, GNTHR = 1.0D-5)
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
C
       DOUBLE PRECISION ATEN, APVAL, PAXIS
       DIMENSION ATEN(MXCENT,3,3), APVAL(MXCENT,3)
       DIMENSION PAXIS(MXCENT,3,3), WORK(*)
       DIMENSION XVAL(6), XBASIS(3,3)
C
       DO IDX = 1,NUCIND
          GNVAL=DABS(DISOTP(NINT(CHARGE(IDX)),ISOTOP(IDX),'GVAL'))
          IF (GNVAL .GT. GNTHR) THEN
             IJ = 1
             DO I = 1,3
                DO J = 1,I
                   IF (.NOT.(I.EQ.J)) THEN
                      TMPVAL = DP5*(ATEN(IDX,I,J)+ATEN(IDX,J,I))
                      ATEN(IDX,I,J) = TMPVAL
                      ATEN(IDX,J,I) = TMPVAL
                   END IF
                   XVAL(IJ) = ATEN(IDX,I,J)
                   IJ = IJ + 1
                END DO
             END DO
             CALL DUNIT(XBASIS,3)
             CALL JACO(XVAL,XBASIS,3,3,3,WORK,WORK(10))
C            Principal values
             APVAL(IDX,1) = XVAL(1)
             APVAL(IDX,2) = XVAL(3)
             APVAL(IDX,3) = XVAL(6)
C            Principal axes
             DO I = 1,3
                DO J = 1,3
                    PAXIS(IDX,I,J) = XBASIS(I,J)
                END DO
             END DO
          END IF
       END DO
C
       RETURN
       END
! -- end of rsp/rsphfc.F --
