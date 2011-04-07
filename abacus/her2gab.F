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
#ifdef REV_LOG
===========================================================================
Revision 1.2  2001/01/31 17:51:20  vebjornb
Parallel modifications from pdahle, MPI and PVM

Revision 1.3  2000/05/04 12:26:07  hjj
(1) merged with Dirac version:
    -- code for 1. derivative screening by Joern Thyssen
    -- bug fix for SUSCEP GAB matrix
    -- bug fix for HKABCD value

Revision 1.2  2000/05/04 10:19:53  hjj
(1) LUPROP .eq. 0 tests changed to LUPROP .le. 0
(2) open LUPROP as unknown in GABGEN
===========================================================================
#endif
C/* Deck gabdrv */
      SUBROUTINE GABDRV(GABM,WORK,LWORK,IJOB,ITYPE,IGTYP,
     &                  JSTRSH,NPRIMS,NCONTS,IORBSH,
     &                  MAXDIF,JATOM,NOCONT,IPRINT)
C*****************************************************************************
C
C     Generate (ij|ij) integrals for screening purposes.
C     IJOB = 0 GAB-integrals in AO-basis
C     IJOB = 1 GAB-integrals in SO-basis
C     This routine based on TWOINT etc.
C
C     Written by T.Saue Oct 12 1996
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <aovec.h>
#include <mxcent.h>
#include <maxorb.h>
#include <dummy.h>
      PARAMETER (D0 = 0.00 D00)
#ifdef PRG_DIRAC
#include <dcbgen.h>
#else
#include <gnrinf.h>
#endif
#include <ccom.h>
#include <symmet.h>
#include <nuclei.h>
#include <blocks.h>
#include <shells.h>
#include <twosta.h>
#include <twocom.h>
#include <hertop.h>
#include <infpar.h>
C
      LOGICAL NOPV, NODV, PERTUR, EXPECT, UNDIFF, DIRFCK, RETUR, TTIME,
     &        NOCONT, SPNORB, DISTRI, SQ12EL, LONDON, SUSCEP, DDFOCK,
     &        ADISTR,DOINT(2,2),SOFOCK, DIA2SO, ZFS2EL
      DIMENSION GABM(NBASIS,NBASIS,*), IORBSH(MXSHEL,MXAOVC),
     &          JSTRSH(MXSHEL,MXAOVC,2), NPRIMS(MXSHEL,MXAOVC,2),
     &          NCONTS(MXSHEL,MXAOVC,2), WORK(LWORK)
C
      CALL QENTER('GABDRV')
#include <memint.h>
      IF (IPRINT .GT. 5) CALL TITLER('Output from GABDRV','*',103)
C
      TKTIME = .FALSE.
C
C     Threshold for integrals
C
C     Because we take the square root of the integrals
C     we need to have a higher threshold.
C     Because elements of value 1.0D10 are not unlikely
C     we set threshold to 1.0D-40, since
C         SQRT (1.0D10 * 1.0D-40) = 1.0D-15
C     which is the threshold we want.
C     The change of threshold have no influence on the time
C     taken to calculate the GAB matrix; we are just going
C     to avoid some numerical instabilities. /Jan 9 1998 - jth
C
C     THRESH = 1.0D-15
      THRESH = 1.0D-40
C
C     Determine run type
C
      IF(ITYPE.EQ.3  .OR.
     &   ITYPE.EQ.4  .OR.
     &   ITYPE.EQ.6  .OR.  ITYPE.EQ.-6  .OR.
     &   ITYPE.EQ.7  .OR.
     &   ITYPE.EQ.8 ) THEN
        WRITE(LUPRI,'(A,I3)') 
     &    'GABDRV does not support  integral type ',ITYPE
        CALL QUIT('GABDRV: Wrong integral type !')
      ENDIF
      CALL TWORUN(ITYPE,MAXDIF,JATOM,PERTUR,EXPECT,UNDIFF,
     &            DIRFCK,SOFOCK,SPNORB,DIA2SO,ZFS2EL,DISTRI,SQ12EL,
     &            LONDON,SUSCEP,DDFOCK,ADISTR,MAXDER,IATOM,MULE,MULTE,
     &            IPRINT)
      CALL WHTREP(ITYPE,MULE,IPRINT)
      IF(IGTYP.EQ.0) THEN
        DOINT(1,1) = .TRUE.
        DOINT(2,1) = .FALSE.
        DOINT(1,2) = .FALSE.
        DOINT(2,2) = .FALSE.
      ELSEIF(IGTYP.EQ.1) THEN
        DOINT(1,1) = .TRUE.
        DOINT(2,1) = .FALSE.
        DOINT(1,2) = .FALSE.
        DOINT(2,2) = .TRUE.  
      ELSE
        WRITE(LUPRI,'(A,I5)') 'GABDRV:Unknown IGTYP ',IGTYP
        CALL QUIT('GABDRV:Unknown IGTYP !')
      ENDIF        
C
C
C     *******************************
C     ***** Calculate integrals *****
C     *******************************
C     
      NGMAT = 1
      IF (MAXDER.GE.1) NGMAT = 3
C     ... MAXDER is set in TWORUN
      CALL DZERO(GABM,NBASIS*NBASIS*NGMAT)
      JTOP  = 4*(NHTYP - 1) + MAXDER
      JTOP3 = (JTOP + 1)**3
      NRTOP = (JTOP + 1)*(JTOP + 2)*(JTOP + 3)/6
      CALL MEMGET('INTE',KINDHR,  JTOP3,WORK,KFREE,LFREE)      
      CALL MEMGET('INTE',KINDSQ,  NRTOP,WORK,KFREE,LFREE)      
      CALL MEMGET('INTE',KIODHR,8*NRTOP,WORK,KFREE,LFREE)      
      CALL HERPRP(WORK(KINDHR),WORK(KINDSQ),WORK(KIODHR))
C
C     *****************************
C     ***** First Shell Index *****
C     *****************************
C
      KBUF = KFREE
      DO 100 ISHELA = 1,MAXSHL
C
        ICA    = LCLASH(ISHELA)
        NHKTA  = NHKTSH(ISHELA)
        KHKTA  = KHKTSH(ISHELA)
        KCKTA  = KCKTSH(ISHELA)
        SPHRA  = SPHRSH(ISHELA)
        NCENTA = NCNTSH(ISHELA)
        MULA   = ISTBSH(ISHELA)
        MULTA  = MULT(MULA)
        NSTRA  = IORBSB(IORBSH(ISHELA,1))
        NUCA   = NUCOSH(ISHELA)
        NORBA  = NORBSH(ISHELA)
        IF (.NOT.BIGVEC) THEN
           CORAX0 = CENTSH(ISHELA,1)
           CORAY0 = CENTSH(ISHELA,2)
           CORAZ0 = CENTSH(ISHELA,3)
        END IF
C
C        Copy everything into third index
C
        ISHELC = ISHELA
        NHKTC  = NHKTA
        KHKTC  = KHKTA
        KCKTC  = KCKTA
        SPHRC  = SPHRA
        NCENTC = NCENTA
        MULC   = MULA
        MULTC  = MULTA
        NSTRC  = NSTRA
        NUCC   = NUCA
        NORBC  = NORBA
        IF (.NOT.BIGVEC) THEN
           CORCX0 = CORAX0
           CORCY0 = CORAY0
           CORCZ0 = CORAZ0
        END IF
C
C       ******************************
C       ***** Second Shell Index *****
C       ******************************
C
        DO 200 ISHELB = 1,ISHELA
C
          ICB    = LCLASH(ISHELB)
          IF(.NOT.DOINT(ICA,ICB)) GOTO 200
          NHKTB  = NHKTSH(ISHELB)
          KHKTB  = KHKTSH(ISHELB)
          KCKTB  = KCKTSH(ISHELB)
          SPHRB  = SPHRSH(ISHELB)
          NCENTB = NCNTSH(ISHELB)
          MULB   = ISTBSH(ISHELB)
          MULTB  = MULT(MULB)
          NSTRB  = IORBSB(IORBSH(ISHELB,1))
          NUCB   = NUCOSH(ISHELB)
          NORBB  = NORBSH(ISHELB)
          IF (.NOT.BIGVEC) THEN
            CORBX0 = CENTSH(ISHELB,1)
            CORBY0 = CENTSH(ISHELB,2)
            CORBZ0 = CENTSH(ISHELB,3)
          END IF
          GENAB  = .NOT.(SEGMSH(ISHELA) .AND. SEGMSH(ISHELB))
          IGENAB = 1
          IF (.NOT.GENAB) IGENAB = 2
          NSETA  = NSETSH(ISHELA,IGENAB)
          NSETB  = NSETSH(ISHELB,IGENAB)
C
C         Copy everything into fourth index
C
          ISHELD = ISHELB
          NHKTD  = NHKTB
          KHKTD  = KHKTB
          KCKTD  = KCKTB
          SPHRD  = SPHRB
          NCENTD = NCENTB
          MULD   = MULB
          MULTD  = MULTB
          NSTRD  = NSTRB
          NUCD   = NUCB
          NORBD  = NORBB
          IF (.NOT.BIGVEC) THEN
            CORDX0 = CORBX0
            CORDY0 = CORBY0
            CORDZ0 = CORBZ0
          END IF
          GENCD  = GENAB
          IGENCD = IGENAB
          NSETC  = NSETA
          NSETD  = NSETB
C
          SHAEQB = ISHELA .EQ. ISHELB
          SHCEQD = SHAEQB
          SHABAB = .TRUE.
          JMAXA  = NHKTA - 1 + MAXDER
          JMAXB  = NHKTB - 1 + MAXDER
          IF (LONDON) JMAXB = NHKTB - 1
          IF (LONDON) JMAXD = NHKTD - 1
          JMAXC  = JMAXA
          JMAXD  = JMAXB
          MAXAB  = NHKTA + NHKTB - 2
          MAXCD  = MAXAB
          TCONAB = SHAEQB .AND. MAXAB .EQ. 0
          TCONCD = TCONAB
          DIAGAB = SHAEQB .AND. .NOT.BIGVEC
          DIAGCD = DIAGAB
          SPHRAB = SPHRA .OR. SPHRB
          SPHRCD = SPHRAB
          DIACAB = DIAGAB .AND. .NOT.SPHRAB
          DIACCD = DIACAB
C
          NODCAB = NODSYM(MAXOPR,MULA,MULB)
          NODCCD = NODSYM(MAXOPR,MULC,MULD)
C
          KCKMAX = MAX(KCKTA,KCKTB)
          CALL MEMGET('INTE',KLMNVL,20*KCKMAX        ,
     &                              WORK,KFREE,LFREE)      
C
C         First electron
C
          LNPCOA = 2*NSETA*(NODCAB + 1)
          LNPCOB = 2*NSETB*(NODCAB + 1)
          CALL MEMGET('INTE',KNPCOA,LNPCOA,WORK,KFREE,LFREE)      
          CALL MEMGET('INTE',KNPCOB,LNPCOB,WORK,KFREE,LFREE)     
          CALL IZERO(WORK(KNPCOA),LNPCOA)
          CALL IZERO(WORK(KNPCOB),LNPCOB)
          CALL MEMGET('INTE',KJSTRA,NSETA               ,
     &                              WORK,KFREE,LFREE)      
          CALL MEMGET('INTE',KJSTRB,NSETB               ,
     &                              WORK,KFREE,LFREE)      
C
          CALL MEMGET('REAL',KCORAB,9*NUCA*NUCB*NODCAB,
     &                              WORK,KFREE,LFREE)      
          CALL MEMGET('REAL',KEXPAB,3*NUCA*NUCB*NODCAB,
     &                              WORK,KFREE,LFREE)      
          CALL MEMGET('REAL',KFACAB,NUCA*NUCB*NODCAB,
     &                              WORK,KFREE,LFREE)      
          IF (GENAB) THEN
            CALL MEMGET('REAL',KCONTA,2*NORBA*NUCA*NODCAB,
     &                                WORK,KFREE,LFREE)      
            CALL MEMGET('REAL',KCONTB,2*NORBB*NUCB*NODCAB,
     &                                WORK,KFREE,LFREE)      
            CALL MEMGET('INTE',KPNTAB,2*NUCA*NUCB*NODCAB ,
     &                                WORK,KFREE,LFREE)      
            CALL MEMGET('INTE',KREDAB,NORBA*NORBB        ,
     &                                WORK,KFREE,LFREE)      
            KNCSAB = KFREE
          ELSE
            KCONTA = KFREE
            KCONTB = KFREE
            KPNTAB = KFREE
            KREDAB = KFREE
            CALL MEMGET('INTE',KNCSAB,NORBA*NORBB*NODCAB,
     &                                WORK,KFREE,LFREE)      
          ENDIF
C
C         Second electron
C
          LNPCOC = 2*NSETC*(NODCCD + 1)
          LNPCOD = 2*NSETD*(NODCCD + 1)
          CALL MEMGET('INTE',KNPCOC,LNPCOC,WORK,KFREE,LFREE)      
          CALL MEMGET('INTE',KNPCOD,LNPCOD,WORK,KFREE,LFREE)      
          CALL IZERO(WORK(KNPCOC),LNPCOC)
          CALL IZERO(WORK(KNPCOD),LNPCOD)
          CALL MEMGET('INTE',KJSTRC,NSETC               ,
     &                              WORK,KFREE,LFREE)      
          CALL MEMGET('INTE',KJSTRD,NSETD               ,
     &                              WORK,KFREE,LFREE)      
C
          CALL MEMGET('REAL',KCORCD,9*NUCC*NUCD*NODCCD,
     &                              WORK,KFREE,LFREE)      
          CALL MEMGET('REAL',KEXPCD,3*NUCC*NUCD*NODCCD,
     &                              WORK,KFREE,LFREE)      
          CALL MEMGET('REAL',KFACCD,NUCC*NUCD*NODCCD,
     &                              WORK,KFREE,LFREE)      
          IF (GENCD) THEN
            CALL MEMGET('REAL',KCONTC,2*NORBC*NUCC*NODCCD,
     &                                WORK,KFREE,LFREE)      
            CALL MEMGET('REAL',KCONTD,2*NORBD*NUCD*NODCCD,
     &                                WORK,KFREE,LFREE)      
            CALL MEMGET('INTE',KPNTCD,2*NUCC*NUCD*NODCCD ,
     &                                WORK,KFREE,LFREE)      
            CALL MEMGET('INTE',KREDCD,NORBC*NORBD        ,
     &                                WORK,KFREE,LFREE)      
            KNCSCD = KFREE
          ELSE
            KCONTC = KFREE
            KCONTD = KFREE
            KPNTCD = KFREE
            KREDCD = KFREE
            CALL MEMGET('INTE',KNCSCD,NORBC*NORBD*NODCCD,
     &                                WORK,KFREE,LFREE)      
          ENDIF
          CALL MEMGET('INTE',KINDAB,2*NORBA*NORBB,
     &                              WORK,KFREE,LFREE)      
          CALL MEMGET('INTE',KINDCD,2*NORBC*NORBD,
     &                              WORK,KFREE,LFREE)      
          CALL GABDR1(GABM,IJOB,WORK(KFREE),LFREE,
     &          UNDIFF,PERTUR,LONDON,SPNORB,DIA2SO,ZFS2EL,EXPECT,
     &          SUSCEP,DDFOCK,DIRFCK,SOFOCK,DISTRI,IATOM,MULE,MULTE,
     &          MAXDER,NOCONT,NODV,NOPV,THRESH,IPRINT,SQ12EL,
     &          JSTRSH,NPRIMS,NCONTS,IORBSH,
     &          WORK(KLMNVL),WORK(KINDSQ),WORK(KIODHR),WORK(KINDHR),
C               Electron 1:
     &          WORK(KCORAB),WORK(KEXPAB),WORK(KFACAB),WORK(KCONTA),
     &          WORK(KCONTB),WORK(KNPCOA),WORK(KNPCOB),WORK(KNCSAB),
     &          WORK(KJSTRA),WORK(KJSTRB),WORK(KINDAB),WORK(KPNTAB),
     &          WORK(KREDAB),
C               Electron 2:
     &          WORK(KCORCD),WORK(KEXPCD),WORK(KFACCD),WORK(KCONTC),
     &          WORK(KCONTD),WORK(KNPCOC),WORK(KNPCOD),WORK(KNCSCD),
     &          WORK(KJSTRC),WORK(KJSTRD),WORK(KINDCD),WORK(KPNTCD),
     &          WORK(KREDCD))
            CALL MEMREL('GABDRV.200',WORK,KWORK,KLMNVL,KFREE,LFREE)
  200    CONTINUE
  100 CONTINUE
C
C     Print section
C
      IF(IPRINT.GE.5) THEN
        IF(IJOB.EQ.0) THEN
          CALL HEADER('Full GAB-matrix in AO basis',-1)
        ELSE
          CALL HEADER('Full GAB-matrix in SO basis',-1)
        ENDIF
        CALL OUTPUT(GABM,1,NBASIS,1,NBASIS,NBASIS,NBASIS,
     &              -4,LUPRI)
        DO IMAT = 2,NGMAT
          WRITE(LUPRI,'(A,I2)')
     &    ' Full derivative GAB-matrix in direction',IMAT-1
          CALL OUTPUT(GABM(1,1,IMAT),1,NBASIS,1,NBASIS,NBASIS,NBASIS,
     &              -4,LUPRI)
         END DO
      ENDIF
      CALL MEMREL('GABDRV',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('GABDRV')
      RETURN
      END
C  /* Deck gabdr1 */
      SUBROUTINE GABDR1(GABM,IJOB,WORK,LWORK,
     &               UNDIFF,PERTUR,LONDON,SPNORB,DIA2SO,ZFS2EL,EXPECT,
     &               SUSCEP,DDFOCK,DIRFCK,SOFOCK,DISTRI,IATOM,MULE,
     &               MULTE,MAXDER,NOCONT,NODV,NOPV,THRESH,IPRINT,SQ12EL,
     &               JSTRSH,NPRIMS,NCONTS,IORBSH,
     &               LMNVLS,INDHSQ,IODDHR,INDHER,
     &               COORAB,EXPAB,FACAB,CONTA,CONTB,NPCOA,NPCOB,
     &               NUCSAB,JSTRA,JSTRB,NINDAB,NPNTAB,NREDAB,
     &               COORCD,EXPCD,FACCD,CONTC,CONTD,NPCOC,NPCOD,
     &               NUCSCD,JSTRC,JSTRD,NINDCD,NPNTCD,NREDCD)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
#include <maxorb.h>
      LOGICAL NOPV, NODV, PERTUR, EXPECT, UNDIFF, DDFOCK, DIRFCK,
     &        DISTRI, NOCONT, LONDON, SPNORB, DIA2SO, ZFS2EL, ONECEN,
     &        SOFOCK, SQ12EL, SUSCEP, TPRAB, TPRCD
      DIMENSION GABM(*),WORK(LWORK),
     &          LMNVLS(KCKMAX,5,4),INDHSQ(*), IODDHR(*), INDHER(*),
     &          IORBSH(MXSHEL,MXAOVC),JSTRSH(MXSHEL,MXAOVC,2), 
     &          NPRIMS(MXSHEL,MXAOVC,2), NCONTS(MXSHEL,MXAOVC,2),
C               Electron 1
     &          COORAB(NUCA*NUCB,3,3,NODCAB),
     &          EXPAB(NUCA*NUCB,3,NODCAB), FACAB(NUCA*NUCB,NODCAB),
     &          CONTA(*), CONTB(*),
     &          NPCOA(NSETA,2,0:NODCAB),NPCOB(NSETB,2,0:NODCAB), 
     &          NUCSAB(NORBA*NORBB,NODCAB), NUCTAB(8),
     &          JSTRA(NSETA), JSTRB(NSETB),
     &          NINDAB(NORBA*NORBB,2), NORTAB(8),
     &          NPNTAB(*),NREDAB(*),
C               Electron 2
     &          COORCD(NUCC*NUCD,3,3,NODCCD),
     &          EXPCD(NUCC*NUCD,3,NODCCD), FACCD(NUCC*NUCD,NODCCD),
     &          CONTC(*), CONTD(*),
     &          NPCOC(NSETC,2,0:NODCCD),NPCOD(NSETD,2,0:NODCCD), 
     &          NUCSCD(NORBC*NORBD,NODCCD), NUCTCD(8),
     &          JSTRC(NSETC), JSTRD(NSETD),
     &          NINDCD(NORBC*NORBD,2), NORTCD(8),
     &          NPNTCD(*),NREDCD(*)
#include <twocom.h>
#include <blocks.h>
#include <symmet.h>
#include <dorps.h>
#include <nuclei.h>
#include <twosta.h>
#include <ibtfun.h>
C
      CALL QENTER('GABDR1')
#include <memint.h>
      IF (IPRINT .GE. 5) CALL TITLER('Output from GABDR1','*',103)
      CALL MEMGET('INTE',KCORBA,NORBA,WORK,KFREE,LFREE)      
      CALL MEMGET('INTE',KCORBB,NORBB,WORK,KFREE,LFREE)      
      CALL MEMGET('INTE',KCORBC,NORBC,WORK,KFREE,LFREE)      
      CALL MEMGET('INTE',KCORBD,NORBD,WORK,KFREE,LFREE)      
C
C     Electron 1
C
      CALL ICOPY(NSETA,NPRIMS(ISHELA,1,IGENAB),MXSHEL,NPCOA(1,1,0),1)
      CALL ICOPY(NSETA,NCONTS(ISHELA,1,IGENAB),MXSHEL,NPCOA(1,2,0),1)
      CALL ICOPY(NSETB,NPRIMS(ISHELB,1,IGENAB),MXSHEL,NPCOB(1,1,0),1)
      CALL ICOPY(NSETB,NCONTS(ISHELB,1,IGENAB),MXSHEL,NPCOB(1,2,0),1)
      CALL ICOPY(NSETA,JSTRSH(ISHELA,1,IGENAB),MXSHEL,JSTRA,1)
      CALL ICOPY(NSETB,JSTRSH(ISHELB,1,IGENAB),MXSHEL,JSTRB,1)
C
C     Overlap distributions for first electron
C     ========================================
C
      TPRAB = .FALSE.
      CALL ODCVEC(COORAB,EXPAB,FACAB,CONTA,CONTB,JMAXA,JMAXB,NSETA,
     &            NSETB,NUCA,NUCB,NUCTAB,NORBA,NORBB,NPCOA,NPCOB,NUCSAB,
     &            JSTRA,JSTRB,TCONAB,TPRAB,GENAB,12,THRESH,MAXDER,
     &            MULA,MULB,NODCAB,NORTAB,NINDAB,NPNTAB,NREDAB,
     &            KHKTA,KHKTB,EXPECT,DIRFCK,WORK(KFREE),LFREE,
     &            RPRIAB,RCNTAB,IDUMMY,IPRINT)
      IF (ISUM(NODCAB,NUCTAB,1) .EQ. 0) GOTO 100
C
C     Electron 2
C
      CALL ICOPY(NSETC,NPRIMS(ISHELC,1,IGENCD),MXSHEL,NPCOC(1,1,0),1)
      CALL ICOPY(NSETC,NCONTS(ISHELC,1,IGENCD),MXSHEL,NPCOC(1,2,0),1)
      CALL ICOPY(NSETD,NPRIMS(ISHELD,1,IGENCD),MXSHEL,NPCOD(1,1,0),1)
      CALL ICOPY(NSETD,NCONTS(ISHELD,1,IGENCD),MXSHEL,NPCOD(1,2,0),1)
      CALL ICOPY(NSETC,JSTRSH(ISHELC,1,IGENCD),MXSHEL,JSTRC,1)
      CALL ICOPY(NSETD,JSTRSH(ISHELD,1,IGENCD),MXSHEL,JSTRD,1)
C
C     Overlap distributions for second electron
C     =========================================
C
      TPRCD = .FALSE.
      CALL ODCVEC(COORCD,EXPCD,FACCD,CONTC,CONTD,JMAXC,JMAXD,NSETC,
     &            NSETD,NUCC,NUCD,NUCTCD,NORBC,NORBD,NPCOC,NPCOD,NUCSCD,
     &            JSTRC,JSTRD,TCONCD,TPRCD,GENCD,34,THRESH,MAXDER,
     &            MULC,MULD,NODCCD,NORTCD,NINDCD,NPNTCD,NREDCD,
     &            KHKTC,KHKTD,EXPECT,DIRFCK,WORK(KFREE),LFREE,
     &            RPRICD,RCNTCD,IDUMMY,IPRINT)
      IF (ISUM(NODCCD,NUCTCD,1) .EQ. 0) GOTO 100
C
C
C
      CALL GETLMN(LMNVLS,IPRINT)
      KHKTAB = KHKTA*KHKTB
      KCKTAB = KCKTA*KCKTB
      IF (DIAGAB) KHKTAB = KHKTA*(KHKTA + 1)/2
      IF (DIACAB) KCKTAB = KCKTA*(KCKTA + 1)/2
      KHKTCD = KHKTAB
      KCKTCD = KCKTAB
      NORBAB = IMXVEC(NORTAB,NODCAB)
      NORBCD = IMXVEC(NORTCD,NODCCD)
      NOABCD = NORBAB*NORBCD
      IF (MAXDER .EQ. 0) THEN
         NCFTYP = 1
      ELSE IF (MAXDER .EQ. 1) THEN
         NCFTYP = 3
      ELSE IF (MAXDER .EQ. 2) THEN
         NCFTYP = 6
      END IF
C
C     Allocate work space
C
C     Electron 1
      MXUCAB = IMXVEC(NUCTAB,NODCAB)
      LCOFAB = MXUCAB*(JMAXA+JMAXB+1)*(JMAXA+1)*(JMAXB+1)*3*NCFTYP
      CALL MEMGET('REAL',KCOFAB,LCOFAB,WORK,KFREE,LFREE)      
C
C     Electron 2
      MXUCCD = IMXVEC(NUCTCD,NODCCD)
      LCOFCD = MXUCCD*(JMAXC+JMAXD+1)*(JMAXC+1)*(JMAXD+1)*3*NCFTYP
      CALL MEMGET('REAL',KCOFCD,LCOFCD,WORK,KFREE,LFREE)      
C
C     Number of SO integrals
C
      IF (EXPECT .AND. MAXDER .GE. 1) THEN
         NINTS     = NMBOSO(SQ12EL,0)
         NINTMX    = NINTS
         NINTSR(1) = NINTS
C         CALL NINTSO(MULE,LONDON,SPNORB,UNDIFF,SOFOCK,DISTRI,SQ12EL,
C     &               IPRINT)
      ELSE IF (EXPECT .OR. DIRFCK .OR. DDFOCK) THEN
         NINTS  = 0
         NINTMX = 0
      ELSE
         CALL NINTSO(MULE,LONDON,SPNORB,UNDIFF,SOFOCK,DISTRI,SQ12EL,
     &               IPRINT)
      END IF
      NSOINT = NOABCD*NINTS
      IF (NSOINT .GT. 0) THEN
C
        CALL MEMGET('REAL',KSOINT,NSOINT         ,WORK,KFREE,LFREE)      
        CALL MEMGET('INTE',KPNTAO,  NINTMX*NOPREP,WORK,KFREE,LFREE)      
        CALL MEMGET('INTE',KPNTOP,3*NINTMX*NOPREP,WORK,KFREE,LFREE)      
        CALL MEMGET('INTE',KPNTNO,4*NINTMX*NOPREP,WORK,KFREE,LFREE)      
        CALL MEMGET('INTE',KPNTRP,3*NINTMX*NOPREP,WORK,KFREE,LFREE)      
        CALL MEMGET('INTE',KPNTLG,3*NINTMX*NOPREP,WORK,KFREE,LFREE)      
        CALL GABDR2(GABM,IJOB,WORK(KSOINT),NSOINT,WORK(KFREE),LFREE,
     &                  UNDIFF,PERTUR,LONDON,SPNORB,DIA2SO,ZFS2EL,
     &                  EXPECT,SUSCEP,DDFOCK,DIRFCK,SOFOCK,DISTRI,IATOM,
     &                  MULE,MULTE,MAXDER,NOCONT,NODV,NOPV,THRESH,
     &                  IPRINT,SQ12EL,COORAB,EXPAB,FACAB,CONTA,CONTB,
     &                  NPCOA,NPCOB,NUCSAB,JSTRA,JSTRB,NINDAB,NPNTAB,
     &                  NREDAB,WORK(KCOFAB),WORK(KCORBA),WORK(KCORBB),
     &                  NUCTAB,COORCD,EXPCD,FACCD,CONTC,CONTD,NPCOC,
     &                  NPCOD,NUCSCD,JSTRC,JSTRD,NINDCD,NPNTCD,NREDCD,
     &                  WORK(KCOFCD),WORK(KCORBC),WORK(KCORBD),NUCTCD,
     &                  WORK(KPNTAO),WORK(KPNTOP),WORK(KPNTNO),
     &                  WORK(KPNTRP),WORK(KPNTLG),
     &                  JSTRSH,NPRIMS,NCONTS,IORBSH,
     &                  INDHSQ,IODDHR,INDHER,LMNVLS)
      END IF
 100  CONTINUE
      CALL MEMREL('GABDR1',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('GABDR1')
      RETURN
      END
C  /* Deck Gabdr2 */
      SUBROUTINE GABDR2(GABM,IJOB,SOINT,NSOINT,WORK,LWORK,
     &                  UNDIFF,PERTUR,LONDON,SPNORB,DIA2SO,ZFS2EL,
     &                  EXPECT,SUSCEP,DDFOCK,DIRFCK,SOFOCK,DISTRI,IATOM,
     &                  MULE,MULTE,MAXDER,
     &                  NOCONT,NODV,NOPV,THRESH,IPRINT,SQ12EL,
     &                  COORAB,EXPAB,FACAB,CONTA,CONTB,NPCOA,NPCOB,
     &                  NUCSAB,JSTRA,JSTRB,NINDAB,NPNTAB,NREDAB,
     &                  COEFAB,ICORBA,ICORBB,NUCTAB,
     &                  COORCD,EXPCD,FACCD,CONTC,CONTD,NPCOC,NPCOD,
     &                  NUCSCD,JSTRC,JSTRD,NINDCD,NPNTCD,NREDCD,
     &                  COEFCD,ICORBC,ICORBD,NUCTCD,
     &                  IPNTAO,IPNTOP,IPNTNO,IPNTRP,IPNTLG,
     &                  JSTRSH,NPRIMS,NCONTS,IORBSH,
     &                  INDHSQ,IODDHR,INDHER,LMNVLS)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxorb.h>
#include <dummy.h>
      PARAMETER (D1 = 1.0D0, D0 = 0.0D0)
C
#include <twocom.h>
#include <blocks.h>
#include <symmet.h>
#include <dorps.h>
#include <nuclei.h>
#include <twosta.h>
#include <expopt.h>
#ifdef PRG_DIRAC
#include <dcbgrd.h>
#else
#include <energy.h>
#endif
#include <doxyz.h>
C
      CHARACTER COMP(2)*1,SPDCAR*1
      LOGICAL AACDX, AACDY, AACDZ, ABCCX, ABCCY, ABCCZ, NOPV, NODV,
     &        PERTUR, EXPECT, INTS, UNDIFF, CAEQCB, CCEQCD, ABADX,
     &        ABADY,  ABADZ, DDFOCK, DIRFCK, DISTRI, NOCONT, SPNORB,
     &        DIA2SO, ZFS2EL, ONECEN, SQ12EL, SOP000, IPNTLG(3,*),
     &        SOFOCK, LONDON, SUSCEP, ADISTR, DOGAB
      DIMENSION GABM(*),SOINT(NSOINT),WORK(LWORK),
     &          INDHSQ(*),IODDHR(*),INDHER(*),LMNVLS(*),
     &          IPNTAO(*), IPNTOP(3,*), IPNTNO(4,*), IPNTRP(3,*),
     &          COORAB(NUCA*NUCB,3,3,NODCAB),
     &          EXPAB(NUCA*NUCB,3,NODCAB), FACAB(NUCA*NUCB,NODCAB),
     &          CONTA(NORBA*NUCA,2,NODCAB), CONTB(NORBB*NUCB,2,NODCAB),
     &          NPCOA(NSETA,2,0:NODCAB), NPCOB(NSETB,2,0:NODCAB),
     &          NUCSAB(NORBA*NORBB,NODCAB),JSTRA(NSETA),JSTRB(NSETB),
     &          NINDAB(*),NPNTAB(NUCA*NUCB,2,NODCAB),NREDAB(*),
     &          COEFAB(*),ICORBA(NORBA),ICORBB(NORBB),NUCTAB(8),
     &          COORCD(NUCC*NUCD,3,3,NODCCD),
     &          EXPCD(NUCC*NUCD,3,NODCCD), FACCD(NUCC*NUCD,NODCCD),
     &          CONTC(NORBC*NUCC,2,NODCCD), CONTD(NORBD*NUCD,2,NODCCD),
     &          NPCOC(NSETC,2,0:NODCCD), NPCOD(NSETD,2,0:NODCCD),
     &          NUCSCD(NORBC*NORBD,NODCCD),JSTRC(NSETC),JSTRD(NSETD),
     &          NINDCD(*),NPNTCD(NUCC*NUCD,2,NODCCD),NREDCD(*),
     &          COEFCD(*),ICORBC(NORBC),ICORBD(NORBD),NUCTCD(8),
     &          IEFFA(0:7),IEFFB(0:7), IEFFC(0:7), IEFFD(0:7),SIGNT(3),
     &          IORBSH(MXSHEL,MXAOVC),JSTRSH(MXSHEL,MXAOVC,2), 
     &          NPRIMS(MXSHEL,MXAOVC,2), NCONTS(MXSHEL,MXAOVC,2)
C
#include <ibtfun.h>
      IBTAXO(I,J) = IBTAND(I,IBTXOR(I,J))
      XAND(I)     = PT(IBTAND(ISYMAX(1,1),I))
      YAND(I)     = PT(IBTAND(ISYMAX(2,1),I))
      ZAND(I)     = PT(IBTAND(ISYMAX(3,1),I))
C
      CALL QENTER('GABDR2')
#include <memint.h>
      IF (IPRINT .GE. 5) CALL TITLER('Output from GABDR2','*',103)
C
      IF (.NOT.DDFOCK.AND. NSOINT.EQ.0) GOTO 400
C
C     This is a GAB-run
      DOGAB = .TRUE.
C
      IOFF = 0
      DO I = 0, MAXOPR
      IF(IBTAND(I,MULA).EQ.0) THEN
        IEFFA(I) = IOFF
        IOFF = IOFF + KHKTA
        ELSE
          IEFFA(I) = IEFFA(IBTXOR(I,IBTAND(I,MULA)))
        END IF
      ENDDO
      CALL SETEFF(IEFFB,IEFFC,IEFFD)
C
      MULAB  = IBTAND(MULA,MULB)
      MULCD  = IBTAND(MULC,MULD)
C
      IF(DIRAC) THEN
        COMP(1) = 'L'
        COMP(2) = 'S'
      ELSE
        COMP(1) = ' '
        COMP(2) = ' '
      ENDIF
C
C     ***** Both distributions *****
C
      JMAX0  = MAXAB + MAXCD
c     HKABCD is a symmetry factor, which determines how many different integrals
c     have the same value. This symmetry factor is not factorizable 
c     between AB and CD, and must therefore not be included in GAB/YAB matrix. 
c     HKABCD = FMULT(IBTAND(MULAB,MULCD)). /Feb 19 1998 - jth
      HKABCD = D1
      IF (SPNORB) HKABCD = - HKABCD
      IF (UNDIFF .OR. LONDON .OR. SPNORB .OR. EXPGRA) THEN
         DOX = .TRUE.
         DOY = .TRUE.
         DOZ = .TRUE.
      ELSE
         DOX = DOCOOR(1,NCENTA) .OR. DOCOOR(1,NCENTB)
         DOY = DOCOOR(2,NCENTA) .OR. DOCOOR(2,NCENTB)
         DOZ = DOCOOR(3,NCENTA) .OR. DOCOOR(3,NCENTB)
      END IF
      DOXYZ(1) = DOX
      DOXYZ(2) = DOY
      DOXYZ(3) = DOZ
C
      INTS = .FALSE.
C
      ICENTA = -1
      ICENTB = -2
      ICENTC = -3
      ICENTD = -4
      SIGNAX = D1
      SIGNAY = D1
      SIGNAZ = D1
      IF (NCENTA.GT.0) ICENTA = NUCNUM(NCNTSH(ISHELA),1)
C
C     Symmetrization loop
C     ===================
C
      IF (.NOT.(EXPECT.OR.DIRFCK.OR.DDFOCK)) THEN
         CALL CPRLOP(IPNTAO,IPNTOP,IPNTNO,IPNTRP,IPNTLG,SQ12EL,IPRINT)
      END IF
C
C     **********************************
C     ***** First Symmetry Index R *****
C     **********************************
C
C     Generates distinct overlap distributions A*R(B)
C
      IF (UNDIFF) THEN
         SOP000 = .TRUE.
      ELSE
         SOP000 = .FALSE.
         IF (PERTUR .OR. LONDON .OR. SPNORB)
     &       CALL DZERO(SOINT,NSOINT)
      END IF
      NSYMR = 0
      DO 100 ISYMR = 0,MAXOPR
      IF (IBTAND(ISYMR,IBTOR(MULA,MULB)) .EQ. 0) THEN
         NSYMR  = NSYMR + 1
         IF(NCENTB.GT.0) ICENTB = NUCNUM(NCNTSH(ISHELB),
     &        IBTAXO(ISYMR,MULB)+1)
         NUCAB  = NUCTAB(NSYMR)
         SIGNBX = XAND(ISYMR)
         SIGNBY = YAND(ISYMR)
         SIGNBZ = ZAND(ISYMR)
C
C        ***************************************************
C        ***** Charge Distributions for First Electron *****
C        ***************************************************
C
         TPRIAB = .FALSE.
c not CAEQCB for derivative GAB matrix, because we
c need one-center integrals to be calculated by
c two-center code because no one-center derivative code
c is programmed (normally one-center is skipped because
c the total contribution is zero).  /jan.98-hjaaj,jth
         CAEQCB = .NOT.BIGVEC .AND. ICENTA .EQ. ICENTB
     &                        .AND. IATOM  .NE. NCENTA
     &                        .AND. MAXDER .EQ. 0
         CALL ODCOEF(COEFAB,COORAB(1,1,1,NSYMR),EXPAB(1,1,NSYMR),
     &             WORK(KFREE),LFREE,
     &             JMAXA,JMAXB,NHKTA,NHKTB,NSETA,NSETB,NUCA,
     &             NUCB,NUCAB,MXUCAB,NORBA,NORBB,NPCOA,NPCOB,
     &             NUCSAB(1,NSYMR),JSTRA,JSTRB,D1,D1,D1,SIGNBX,SIGNBY,
     &             SIGNBZ,CORAX0,CORAY0,CORAZ0,CORBX0,CORBY0,CORBZ0,
     &             AACDX,AACDY,AACDZ,IAB0X,IAB0Y,IAB0Z,CAEQCB,.TRUE.,
     &             .TRUE.,BIGVEC,UNDIFF,LONDON,SPNORB,DIA2SO,ZFS2EL,
     &             12,THRESH,MAXDER,IPRINT)
         NSYMS = 0
         DO 200 ISYMS = 0, MAXOPR
         IF (IBTAND(ISYMS,IBTOR(MULC,MULD)) .EQ. 0) THEN
           IF(IJOB.EQ.0) THEN
             IF(ISYMS.NE.ISYMR) THEN
               GOTO 200
             ELSE
               NSYMS = NSYMR
             ENDIF
           ELSE
             NSYMS = NSYMS + 1
           ENDIF
           TPRICD = .FALSE.
C
C        ****************************************************
C        ***** Charge Distributions for Second Electron *****
C        ****************************************************
C
           SIGNDX = XAND(ISYMS)
           SIGNDY = YAND(ISYMS)
           SIGNDZ = ZAND(ISYMS)
           NUCCD = NUCTCD(NSYMS)
           CALL ODCOEF(COEFCD,COORCD(1,1,1,NSYMS),EXPCD(1,1,NSYMS),
     &                 WORK,LWORK,JMAXC,JMAXD,NHKTC,NHKTD,NSETC,NSETD,
     &                 NUCC,NUCD,NUCCD,MXUCCD,NORBC,NORBD,NPCOC,NPCOD,
     &                 NUCSCD(1,NSYMS),JSTRC,JSTRD,D1,D1,D1,SIGNDX,
     &                 SIGNDY,SIGNDZ,CORCX0,CORCY0,CORCZ0,CORDX0,
     &                 CORDY0,CORDZ0,ABCCX,ABCCY,ABCCZ,ICD0X,ICD0Y,
     &                 ICD0Z,.FALSE.,.TRUE.,.TRUE.,BIGVEC,UNDIFF,
     &                 LONDON,SPNORB,DIA2SO,ZFS2EL,34,THRESH,
     &                 MAXDER,IPRINT)
            IF (NOCONT) THEN
              NCCINT = NUCAB *NUCCD *KHKTAB*KHKTCD
            ELSE
              NCCINT = NORBAB*NORBCD*KHKTAB*KHKTCD
            END IF
            DO 300 ISYMT = 0, MAXOPR
            IF (IBTAND(ISYMT,IBTOR(MULAB,MULCD)) .EQ. 0) THEN
              IF(IJOB.EQ.0) THEN
                IF(ISYMT.NE.0) GOTO 300
              ENDIF
              ISYMTS = IBTXOR(ISYMT,ISYMS)
              IF(NCENTC.GT.0) ICENTC = NUCNUM(NCNTSH(ISHELC),
     &            IBTAXO(ISYMT,MULC)+1)
              IF(NCENTD.GT.0) ICENTD = NUCNUM(NCNTSH(ISHELD),
     &            IBTAXO(ISYMTS,MULD)+1)
C
C             ****************************************************
C             ***** Charge Distributions for Second Electron *****
C             ****************************************************
C
              SIGNT(1) = XAND(ISYMT)
              SIGNT(2) = YAND(ISYMT)
              SIGNT(3) = ZAND(ISYMT)
              SIGNCX = XAND(ISYMT)
              SIGNCY = YAND(ISYMT)
              SIGNCZ = ZAND(ISYMT)
              SIGNDX = XAND(ISYMTS)
              SIGNDY = YAND(ISYMTS)
              SIGNDZ = ZAND(ISYMTS)
              CCEQCD = .NOT.BIGVEC .AND. (ICENTC .EQ. ICENTD)
     &                             .AND. (NCENTC .NE. IATOM)
C
C             If necessary change sign of expansion coefficients
C             ==================================================
C
              CALL EXCSGN(COEFCD,JMAXC,JMAXD,NHKTC,NHKTD,NUCCD,
     &                    MXUCCD,ICD0X,ICD0Y,ICD0Z,MAXDER,LONDON,
     &                    SPNORB,DIA2SO,EXPGRA,ISYMT,IPRINT)
C
C             Check whether this integral gives zero contribution
C             ===================================================
C
C
C             b) atom to be differentiated does not enter integral
C
              IF (PERTUR .AND. MAXDER .GT. 0) THEN
                IF (IATOM .NE. NCENTA .AND. IATOM .NE. NCENTB)
     &                       GO TO 100
              END IF
C
C             c) one-center integrals
C                We actually need one-center integrals for
C                the GAB-matrix for screening of the gradient
C
              ONECEN = (ICENTA .EQ. ICENTB) .AND.
     &                 (ICENTA .EQ. ICENTC) .AND.
     &                 (ICENTA .EQ. ICENTD) .AND.
     &                 (ICENTA .NE. 0)
              ONECEN = ONECEN .AND. .NOT. DOGAB
              IF (ONECEN) THEN
                IF (MAXDER .EQ. 0 .OR. EXPGRA) THEN
                  IF (IBTAND(JMAX0,1).EQ.1)              GO TO 100
                ELSE
                  IF (PERTUR .OR. EXPECT .OR. LONDON)    GO TO 100
                END IF
             END IF
             IF (LONDON) THEN
               IF ((ICENTA .EQ. ICENTB) .AND.
     &             (ICENTC .EQ. ICENTD)) GO TO 100
             END IF
C
C            d) none of the directions are requested
C
             IF (.NOT. (DOX .OR. DOY .OR. DOZ)) GO TO 100
C
C            Integral contributes
C            ====================
C
C            Local symmetries
C            ================
C
             IF (BIGVEC) THEN
               CALL CORDIF(NSETA,NSETC,THRESH,ABADX,ABADY,ABADZ,
     &                 IPRINT,NPCOA(1,1,NSYMR),NPCOC(1,1,NSYMS),
     &                 JSTRA,JSTRC)
             ELSE
               ABADX  = ABS(CORAX0 - XAND(ISYMT)*CORCX0) .LT. THRESH
               ABADY  = ABS(CORAY0 - YAND(ISYMT)*CORCY0) .LT. THRESH
               ABADZ  = ABS(CORAZ0 - ZAND(ISYMT)*CORCZ0) .LT. THRESH
             END IF
             ISAMEX = 0
             ISAMEY = 0
             ISAMEZ = 0
             IF (AACDX .AND. ABADX .AND. ABCCX) ISAMEX = 1
             IF (AACDY .AND. ABADY .AND. ABCCY) ISAMEY = 1
             IF (AACDZ .AND. ABADZ .AND. ABCCZ) ISAMEZ = 1
              ISMXYZ = ISAMEX + 2*ISAMEY + 4*ISAMEZ
C
CTROND           IF (EXPECT .OR. SUSCEP .OR. DIRFCK .OR. DDFOCK) THEN
             DO 55 I = 1, NORBA
               ICORBA(I) = IORBSH(ISHELA,I)
   55        CONTINUE
             DO 65 I = 1, NORBB
               ICORBB(I) = IORBSH(ISHELB,I) + IEFFB(ISYMR)
   65        CONTINUE
             DO 75 I = 1, NORBC
               ICORBC(I) = IORBSH(ISHELC,I) + IEFFC(ISYMT)
   75        CONTINUE
             DO 85 I = 1, NORBD
               ICORBD(I) = IORBSH(ISHELD,I) + IEFFD(ISYMTS)
   85        CONTINUE
CTROND           END IF
C
C            *******************************
C            ***** Integral Directives *****
C            *******************************
C
             CALL DIRECT(BIGVEC,ICENTA,ICENTB,ICENTC,ICENTD,
     &                   NCENTA,NCENTB,NCENTC,NCENTD,
     &                   0,     ISYMR, ISYMT, ISYMTS,
     &                   SIGNAX,SIGNAY,SIGNAZ,SIGNBX,SIGNBY,SIGNBZ,
     &                   SIGNCX,SIGNCY,SIGNCZ,SIGNDX,SIGNDY,SIGNDZ,
     &                   NCCINT,MAXDER,EXPECT,LONDON,SPNORB,DIA2SO,
     &                   ZFS2EL,EXPGRA,IATOM,MULTE,NINTYP,IPRINT)
C
C            ************************
C            ***** AO integrals *****
C            ************************
C
             LAOINT = NCCINT*NINTYP
             CALL MEMGET('REAL',KAOINT,LAOINT,WORK,KFREE,LFREE)      
             CALL CAOINT(SOINT,DUMMY,1,DUMMY,DUMMY,DUMMY,WORK(KAOINT),
     &                   WORK(KFREE),LFREE,COEFAB,COEFCD,
     &                   COORAB(1,1,1,NSYMR),COORCD(1,1,1,NSYMS),
     &                   EXPAB(1,1,NSYMR),EXPCD(1,1,NSYMS),
     &                   FACAB(1,NSYMR),FACCD(1,NSYMS),
     &                   CONTA(1,1,NSYMR),CONTB(1,1,NSYMR),
     &                   CONTC(1,1,NSYMS),CONTD(1,1,NSYMS),
     &                   NCCINT,NINTYP,UNDIFF,PERTUR,LONDON,SPNORB,
     &                   DIA2SO,ZFS2EL,EXPECT,SUSCEP,DDFOCK,DIRFCK,
     &                   SOFOCK,DISTRI,IATOM,
     &                   MULE,MULTE,MAXDER,BIGVEC,NOCONT,NODV,NOPV,
     &                   ISMXYZ,THRESH,ONECEN,IPRINT,ICORBA,ICORBB,
     &                   ICORBC,ICORBD,INTS,HKABCD,JMAX0,ISYMR,ISYMS,
     &                   ISYMT,ISYMTS,SQ12EL,SOP000,IPNTAO,IPNTOP,
     &                   NPCOA(1,1,NSYMR),NPCOB(1,1,NSYMR),
     &                   NPCOC(1,1,NSYMS),NPCOD(1,1,NSYMS),
     &                   NUCSAB(1,NSYMR), NUCSCD(1,NSYMS),
     &                   SIGNT,INDHSQ,IODDHR,INDHER,LMNVLS,NINDAB,
     &                   NINDCD,NPNTAB(1,1,NSYMR),NPNTCD(1,1,NSYMS),
     &                   NREDAB,NREDCD,IDUMMY,DUMMY,DUMMY,DUMMY,
     &                   NSOINT,IDUMMY,IDUMMY)
C
C            ********************************
C            ***** Process AO integrals *****
C            ********************************
C
             IF(IJOB.EQ.0) THEN
                IF (MAXDER.EQ.0) THEN
                   CALL GAOGAT(GABM,WORK(KAOINT),DUMMY,
     &                         NCCINT,NINTYP,IORBSH,HKABCD,
     &                         ICORBA,ICORBB,ICORBC,ICORBD,IEFFA,IEFFB,
     &                         ISYMR,IPRINT,NINDAB,NINDCD,SUSCEP,MAXDER)
                ELSE
                   DO 111 I = 1,3
                      IF (.NOT. DOXYZ(I)) GOTO 111
                      JGABM = 1 + (I-1)*NBASIS*NBASIS
                      JAOINT= KAOINT + (I-1)*NCCINT
                      JAOIN2= KAOINT + (I+2)*NCCINT
                      CALL GAOGAT(GABM(JGABM),WORK(JAOINT),WORK(JAOIN2),
     &                         NCCINT,NINTYP,IORBSH,HKABCD,
     &                         ICORBA,ICORBB,ICORBC,ICORBD,IEFFA,IEFFB,
     &                         ISYMR,IPRINT,NINDAB,NINDCD,SUSCEP,MAXDER)
 111               CONTINUE
                END IF
             ENDIF
             CALL MEMREL('GABDR2',WORK,KWORK,KWORK,KFREE,LFREE)
           END IF
  300      CONTINUE
         END IF
  200    CONTINUE
      END IF
  100 CONTINUE
C
C     ********************************
C     ***** Process SO integrals *****
C     ********************************
C
  400 CONTINUE
      IF (IJOB.EQ.1.AND.INTS) THEN
CTROND Legg inn utpaking av integraler her !!!!!
C
C        A) Undifferentiated integrals
C        =============================
C
         IF (UNDIFF) THEN
           CALL GABGAT(GABM,SOINT,IPNTNO,IPNTRP,NINDAB,
     &                 NINDCD,IPRINT)
C
C        B) Differentiated integrals
C        ===========================
C
         ELSE IF (PERTUR .OR. SPNORB .OR. DIA2SO .OR. ZFS2EL) THEN
            WRITE(LUPRI,'(A)') 'Not yet implemented !'
            CALL QUIT('Not yet implemented !')
C
C        C) London orbitals
C        ==================
C
         ELSE IF (LONDON) THEN
            WRITE(LUPRI,'(A)') 'Not yet implemented !'
            CALL QUIT('Not yet implemented !')
         END IF
      ENDIF
      CALL QEXIT('GABDR2')
      RETURN
C
      END
C  /* Deck gabgat */
      SUBROUTINE GABGAT(GABM,SO,IPNTNO,IPNTRP,NINDAB,
     &                  NINDCD,IPRINT)
C
C     Collect elements of GAB-matrix
C     Written by T.Saue Oct 7 1996
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <aovec.h>
#include <nuclei.h>
      DIMENSION GABM(NBASIS,NBASIS),
     &          SO(*),IPNTNO(4,*),
     &          IPNTRP(3,*),
     &          NINDAB(NORBA*NORBB,2),NINDCD(NORBC*NORBD,2)
#include <twocom.h>
#include <symmet.h>
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine GABPCK',-1)
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(2X,A,4I5)') 'NHKT? ', NHKTA, NHKTB, NHKTC, NHKTD
         WRITE (LUPRI,'(2X,A,4I5)') 'MUL?  ', MULA,  MULB,  MULC,  MULD
         WRITE (LUPRI,'(2X,A,4I5)') 'NORB? ', NORBA, NORBB, NORBC, NORBD
         WRITE (LUPRI,'(2X,A,4I5)') 'NSTR? ', NSTRA, NSTRB, NSTRC, NSTRD
         WRITE (LUPRI,'(2X,A,2I5)') 'NORBCD', NORBCD
         WRITE (LUPRI,'(2X,A,2I5)') 'NOABCD', NOABCD
         WRITE (LUPRI,'(2X,A,2L5)') 'DIAGAB/CD', DIAGAB, DIAGCD
         WRITE (LUPRI,'(2X,A,2L5)') 'TCONAB/CD', TCONAB, TCONCD
         WRITE (LUPRI,'(2X,A,2L5)') 'SHAEQB/CD', SHAEQB, SHCEQD
         WRITE (LUPRI,'(2X,A, L5)') 'SHABAB', SHABAB
      END IF
C
C
      ISOFF = 0
      DO 100 I = 1, NINTS
         NSTRNA = IPNTNO(1,I)
         NSTRNB = IPNTNO(2,I)
         NSTRNC = IPNTNO(3,I)
         NSTRND = IPNTNO(4,I)
         IREPA  = IPNTRP(1,I)
         IREPB  = IPNTRP(2,I)
         IREPC  = IPNTRP(3,I)
         IREPD  = IBTXOR(IBTXOR(IREPA,IREPB),IREPC)
         INT = 0
         DO 200 IAB = 1, NORBAB
           IA = KHKTA*(NINDAB(IAB,1) - 1)
           IB = KHKTB*(NINDAB(IAB,2) - 1)
           INDA  = IPTSYM(NSTRNA + IA,IREPA)
           INDB  = IPTSYM(NSTRNB + IB,IREPB)
           IMAX  = MAX(INDA,INDB)
           INDB  = MIN(INDA,INDB)
           INDA  = IMAX
           DO 300 ICD = 1,NORBCD
             INT = INT + 1
             IC = KHKTC*(NINDCD(ICD,1) - 1)
             ID = KHKTD*(NINDCD(ICD,2) - 1)
             INDC = IPTSYM(NSTRNC + IC,IREPC)
             INDD = IPTSYM(NSTRND + ID,IREPD)
             IMAX = MAX(INDC,INDD)
             INDD = MIN(INDC,INDD)
             INDC = IMAX
             IF(INDA.EQ.INDC.AND.INDB.EQ.INDD) THEN
               GABM(INDA,INDB) = SO(ISOFF+INT)
               GABM(INDB,INDA) = SO(ISOFF+INT)
             ENDIF
  300      CONTINUE
  200    CONTINUE
         ISOFF = ISOFF + NOABCD
  100 CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Gaogat */
      SUBROUTINE GAOGAT(GABM,AOINT,AOINT2,NCCINT,NINTYP,IORBSH,SYMFAC,
     &                  ICORBA,ICORBB,ICORBC,ICORBD,IEFFA,IEFFB,ISYMR,
     &                  IPRINT,NINDAB,NINDCD,SUSCEP,MAXDER)
C*****************************************************************************
C
C  Gather AO-Cauchy-Schwarz-integrals
C
C  Written by T.Saue Oct 10 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <aovec.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER(D1 = 1.0D0)
      INTEGER A,B,C,D
      LOGICAL SUSCEP
#include <nuclei.h>
      DIMENSION ICORBA(NORBA),ICORBB(NORBB),
     &          ICORBC(NORBC),ICORBD(NORBD),
     &          AOINT(NCCINT), GABM(NBASIS,NBASIS),
     &          AOINT2(NCCINT),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2),
     &          IEFFA(0:7),IEFFB(0:7),IORBSH(MXSHEL,MXAOVC)
#include <twocom.h>
#include <expcom.h>
#include <symmet.h>
#include <blocks.h>
C
#include <ibtfun.h>
C
Chj/980617: changed ISYMA test from MULA to MULAB, we must
C           include all different operations where at least one of
C           A and B move (because we need the GAB integrals
C           for sqrt( T(C S(D)) | T(C S(D)) ) ).
      IF(SYMFAC.NE.D1) CALL DSCAL(NCCINT,SYMFAC,AOINT,1)
      IOFF = 0
      MULAB = IBTAND(MULA,MULB)
      IF(SUSCEP) THEN
        DO ICOMPA = 1,KHKTA
                      KHKTBB = KHKTB
          IF (DIAGAB) KHKTBB = ICOMPA
          DO ICOMPB = 1,KHKTBB
            DO ICOMPC = 1,KHKTC
                          KHKTDD = KHKTD
              IF (DIAGCD) KHKTDD = ICOMPC
              DO ICOMPD = 1,KHKTDD
                DO IORBAB = 1,NORBAB
                  IORBA = 1 + NINDAB(IORBAB,1)/KHKTA
                  IORBB = 1 + NINDAB(IORBAB,2)/KHKTB
                  A = ICORBA(IORBA) + ICOMPA
                  B = ICORBB(IORBB) + ICOMPB
                  DO IORBCD = 1, NORBCD
                    IORBC = 1 + NINDCD(IORBCD,1)/KHKTC
                    IORBD = 1 + NINDCD(IORBCD,2)/KHKTD
                    C = ICORBC(IORBC) + ICOMPC
                    D = ICORBD(IORBD) + ICOMPD
                    IF(A.EQ.C.AND.B.EQ.D) THEN
                      AOBUF     = AOINT(IOFF+IORBCD)
                      DO ISYMA = 0, MAXOPR
                      IF(IBTAND(ISYMA,MULAB).EQ.0) THEN
                         ISYMB = IBTXOR(ISYMR,ISYMA)
                         IA = IORBSH(ISHELA,IORBA) 
     &                       + IEFFA(ISYMA) + ICOMPA
                         IB = IORBSH(ISHELB,IORBB) 
     &                       + IEFFB(ISYMB) + ICOMPB
                         GABM(IA,IB) = AOBUF
                         GABM(IB,IA) = AOBUF
                      END IF
                      ENDDO
                    ENDIF
                  ENDDO
                  IOFF = IOFF + NORBCD
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO ICOMPA = 1,KHKTA
                      KHKTBB = KHKTB
          IF (DIAGAB) KHKTBB = ICOMPA
          DO ICOMPB = 1,KHKTBB
            DO ICOMPC = 1,KHKTC
                          KHKTDD = KHKTD
              IF (DIAGCD) KHKTDD = ICOMPC
              DO ICOMPD = 1,KHKTDD
                DO IORBAB = 1,NORBAB
                  IORBA = NINDAB(IORBAB,1)
                  IORBB = NINDAB(IORBAB,2)
                  A = ICORBA(IORBA) + ICOMPA
                  B = ICORBB(IORBB) + ICOMPB
                  DO IORBCD = 1, NORBCD
                    IORBC = NINDCD(IORBCD,1)
                    IORBD = NINDCD(IORBCD,2)
                    C = ICORBC(IORBC) + ICOMPC
                    D = ICORBD(IORBD) + ICOMPD
                    IF(A.EQ.C.AND.B.EQ.D) THEN
                      AOBUF     = AOINT(IOFF+IORBCD)
                      IF (NINTYP .EQ. 6 .AND. MAXDER .EQ. 2) THEN
                         AOBUF2    = AOINT2(IOFF+IORBCD)
                      ELSE
                         AOBUF2    = AOBUF
                      END IF
                      DO ISYMA = 0, MAXOPR
                      IF(IBTAND(ISYMA,MULAB).EQ.0) THEN
                         ISYMB = IBTXOR(ISYMR,ISYMA)
                         IA = IORBSH(ISHELA,IORBA) 
     &                     + IEFFA(ISYMA) + ICOMPA
                         IB = IORBSH(ISHELB,IORBB) 
     &                     + IEFFB(ISYMB) + ICOMPB
                         GABM(IA,IB) = AOBUF
                         GABM(IB,IA) = AOBUF2
                      END IF
                      ENDDO
                    ENDIF
                  ENDDO
                  IOFF = IOFF + NORBCD
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mkdrao */
      SUBROUTINE MKDRAO(DMAT,DRAO,NDMAT,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Compress a matrix DMAT in AO basis to a reduced matrix DRAO
C     over symmetry unique block indices as defined in 
C     COMMON block BLOCKS (dimension NSYMBL)
C
C     Written by T.Saue Oct 12 1995
C
C     Last revision: Oct 3 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0=0.0D0)
C
#include <aovec.h>
#include <maxorb.h>
#include <mxcent.h>
C
#include <shells.h>
#include <blocks.h>
#include <nuclei.h>
      DIMENSION DMAT(NBASIS*NBASIS,NDMAT),
     &          DRAO(NSYMBL,NSYMBL,NDMAT),WORK(LWORK)
C
      CALL QENTER('MKDRAO')
#include <memint.h>
C
C     Initialization
C
      NDRAO   = NSYMBL*NSYMBL*NDMAT
      CALL DZERO(DRAO,NDRAO)
C
C     Make index array from AOs to blocks
C
      CALL MEMGET('INTE',KINRAO,NBASIS,WORK,KFREE,LFREE)
      CALL PTRDAO(WORK(KINRAO))
C
C     Make reduced AO matrix
C
      DO I = 1,NDMAT
        CALL GATMAT(0,NBASIS,NBASIS,DMAT(1,I),
     &              DRAO(1,1,I),WORK(KINRAO),NSYMBL)
      ENDDO
C
C     Print section
C     =============
C
      IF(IPRINT.GE.4) CALL PRDRAO(DRAO,NDMAT,WORK(KINRAO))
C
      CALL QEXIT('MKDRAO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Prdrao */
      SUBROUTINE PRDRAO(DRAO,NDMAT,INDEX)
C*****************************************************************************
C
C     Print reduced matrix DRAO
C
C     Written by T.Saue Oct 11 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <aovec.h>
#include <maxorb.h>
#include <mxcent.h>
C
#include <shells.h>
#include <blocks.h>
#include <nuclei.h>
      DIMENSION DRAO(NSYMBL,NSYMBL,NDMAT),INDEX(NSYMBL)
      CHARACTER SPDCAR*1,COMP(2)*1
      COMP(1) = 'L'
      COMP(2) = 'S'
C
C     First make pointer from symmetry dependent blocks to 
C     symmetry independent blocks
C     
      NBL = 0
      DO I = 1,MAXSHL
        NCENTA = NCNTSH(I)
        NDEG = NUCDEG(NCENTA)
        DO IDEG = 1,NDEG
          NBL = NBL + 1
          INDEX(NBL) = I
        ENDDO
      ENDDO
C
C     Print section
C
      CALL HEADER('Output from MKDRAO',-1)
      DO K = 1,NDMAT
        WRITE(LUPRI,'(A,I5)') '* Reduced AO-matrix no.',K
        JEND = MIN(NSYMBL,8)
        DO J = 1,NSYMBL,8
          WRITE(LUPRI,'(12X,8(I4,A1,A3,1X,3A1))')
     &         (JJ,'(',NAMN(NCNTSH(INDEX(JJ)))(1:3),
     &          COMP(LCLASH(INDEX(JJ))),
     &          SPDCAR(NHKTSH(INDEX(JJ))-1),')',JJ = J,JEND)
          DO I = 1,NSYMBL
             WRITE(LUPRI,'(I4,A1,A3,1X,3A1,1P,8D12.3)')
     &         I,'(',NAMN(NCNTSH(INDEX(I)))(1:3),
     &         COMP(LCLASH(INDEX(I))),
     &         SPDCAR(NHKTSH(INDEX(I))-1),')',
     &        (DRAO(I,II,K),II = J,JEND)
          ENDDO
          JEND = MIN(JEND+8,NSYMBL)
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Indrao */
      SUBROUTINE PTRDAO(INDRAO)
C*****************************************************************************
C
C     Make index array pointing from AOs to all blocks
C     (that is, not only the symmetry unique blocks of MAXSHL,
C      but the complete list in NSYMBL)
C
C     Written by T.Saue Oct 4 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <aovec.h>
#include <maxorb.h>
#include <mxcent.h>
C
#include <shells.h>
#include <blocks.h>
#include <nuclei.h>
      DIMENSION INDRAO(NBASIS)
C
      IOFF = 0
      DO IA = 1,KMAX
        ISHELA = IPTSHL(IA)
        KHKTA  = KHKTSH(ISHELA)
        ICENTA = NCNTSH(ISHELA)
        NDEG   = NUCDEG(ICENTA)
        DO IDEG = 1,NDEG
          ISYBLA = ISYMBL(ISHELA,IDEG)
          DO II = 1,KHKTA
            INDRAO(IOFF+II) = ISYBLA
          ENDDO
          IOFF = IOFF + KHKTA
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Indrso */
      SUBROUTINE PTRDSO(INDRSO)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Make index array pointing from SOs to symmetry unique blocks
C     (dimension MAXSHL) defined in COMMON block BLOCKS
C
C     Written by T.Saue Oct 11 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <aovec.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
C
#include <symmet.h>
#include <pincom.h>
#include <blocks.h>
#include <shells.h>
#include <nuclei.h>
      DIMENSION INDRSO(NBASIS)
C
#include <ibtfun.h>
C
      ISTRA = 1
      DO IREPA = 0, MAXREP
        NORBA = NAOS(IREPA+1)
        DO I = ISTRA,ISTRA + NORBA - 1
          ISHELA   = IPTSHL(IBTAND(IBTSHR(IPIND(I),16),65535))
          INDRSO(I) = ISHELA
        ENDDO
        ISTRA = ISTRA + NORBA
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mkdrso */
      SUBROUTINE MKDRSO(DMAT,DRSO,NDMAT,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Compress a matrix DMAT in SO basis to a reduced matrix DRSO
C     over symmetry unique block indices as defined in 
C     COMMON block BLOCKS (dimension MAXSHL)
C
C     Written by T.Saue Oct 11 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0=0.0D0)
C
#include <aovec.h>
#include <maxorb.h>
#include <mxcent.h>
C
#include <shells.h>
#include <blocks.h>
#include <nuclei.h>
      DIMENSION DMAT(NBASIS*NBASIS,NDMAT),
     &          DRSO(MAXSHL,MAXSHL,NDMAT),WORK(LWORK)
      CHARACTER SPDCAR*1,COMP(2)*1
C
      CALL QENTER('MKDRSO')
#include <memint.h>
C
C     Initialization
C
      COMP(1) = 'L'
      COMP(2) = 'S'
      NDRSO   = MAXSHL*MAXSHL*NDMAT
      CALL DZERO(DRSO,NDRSO)
C
C     Make index array from SOs to blocks
C
      CALL MEMGET('INTE',KINRSO,NBASIS,WORK,KFREE,LFREE)
      CALL PTRDSO(WORK(KINRSO))
C
C     Make reduced SO matrix
C
      DO I = 1,NDMAT
        CALL GATMAT(0,NBASIS,NBASIS,DMAT(1,I),
     &              DRSO(1,1,I),WORK(KINRSO),MAXSHL)
      ENDDO
C
C     Print section
C     =============
C
      IF(IPRINT.GE.4) THEN
        CALL HEADER('Output from MKDRSO',-1)
        DO K = 1,NDMAT
          WRITE(LUPRI,'(A,I5)') '* Reduced density matrix no.',K
          JEND = MIN(MAXSHL,8)
          DO J = 1,MAXSHL,8
             WRITE(LUPRI,'(12X,8(I4,A1,A3,1X,3A1))')
     &           (JJ,'(',NAMN(NCNTSH(JJ))(1:3),COMP(LCLASH(JJ)),
     &            SPDCAR(NHKTSH(JJ)-1),')',JJ = J,JEND)
            DO I = 1,MAXSHL
               WRITE(LUPRI,'(I4,A1,A3,1X,3A1,1P,8D12.3)')
     &           I,'(',NAMN(NCNTSH(I))(1:3),COMP(LCLASH(I)),
     &            SPDCAR(NHKTSH(I)-1),')',
     &           (DRSO(II,I,K),II = J,JEND)
            ENDDO
            JEND = MIN(JEND+8,MAXSHL)
          ENDDO
        ENDDO
      ENDIF
C
      CALL QEXIT('MKDRSO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck scrsta */
      SUBROUTINE SCRSTA(TEXT,DINTSKP,CPUSEC)
C*****************************************************************************
C
C     Provide statistics about screening
C
C       Screening proceeds in three steps as documented by DINTSKP:
C         Step 1: Screening on integral batches
C           DINTSKP(1,1) - total number of integrals
C           DINTSKP(2,1) - number of integrals skipped (batchwise)
C         Step 2: Screening on individual integrals 
C                 while unpacking indices
C           DINTSKP(1,2) - number of integrals remaining after step 1
C           DINTSKP(2,2) - number of integrals skipped
C         Step 3a: Screening on Coulomb contributions
C           NCM is the number of Fock matrices with Coulomb contribution
C           DINTSKP(1,3) - NCM times number of integrals remaining 
C                         after step 2
C           DINTSKP(2,3) - number of integrals skipped,
C                         the number is obtained by counting the number of
C                         integrals skipped for each NCM separately
C         Step 3b: Screening on exchange contributions
C           DINTSKP(1,4) - NEM times number of integrals remaining 
C                         after step 2
C           DINTSKP(2,4) - number of integrals skipped
C                         the number is obtained by counting the number of
C                         integrals skipped for each NEM separately
C
C     Written by T.Saue Oct 11 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DC = 100.0D0,D1 = 1.0D0,D0 = 0.0D0)
      DIMENSION DINTSKP(2,4),DSTEP(4)
      CHARACTER TEXT*(*)
      CHARACTER SECTID*12, INTTID*12
C
      IF(TEXT(1:6).EQ.'Header') THEN
        WRITE(LUPRI,'(A3,6X,4(5X,A8),3X,A8)') 
     &    'SCR','Step1   ','Step2   ',
     &    'Coulomb ','Exchange','CPU-time'
        RETURN
      ENDIF
      CALL DZERO(DSTEP,4)
      IF(DINTSKP(1,1).GT.D0) DSTEP(1) = DINTSKP(2,1)/DINTSKP(1,1)
      DA       = D1 - DSTEP(1)
      DSTEP(1) = DC*DSTEP(1)
      IF(DINTSKP(1,2).GT.D0) DSTEP(2) = 
     &         (DINTSKP(2,2)/DINTSKP(1,2))*DA
      DA       = DC*DA*(D1 - DSTEP(2))
      DSTEP(2) = DC*DSTEP(2)
      IF(DINTSKP(1,3).GT.D0) DSTEP(3) 
     &         = (DINTSKP(2,3)/DINTSKP(1,3))*DA
      IF(DINTSKP(1,4).GT.D0) 
     &           DSTEP(4) = (DINTSKP(2,4)/DINTSKP(1,4))*DA
C
      INTTID = SECTID(CPUSEC)
#ifndef PRG_DIRAC
C
C
         WRITE(LUPRI,'(2A)') '* Screening statistics ',TEXT
         WRITE(LUPRI,'(4(5X,A8),3X,A9)') 'Step1   ','Step2   ',
     &        'Coulomb ','Exchange','WALL-time'
         WRITE(LUPRI,'(4(6X,F6.2,A1),3X,A12)') 
     &        (DSTEP(I),'%',I=1,4),INTTID
#else /* PRG_DIRAC */
      WRITE(LUPRI,'(A9,4(6X,F6.2,A1),3X,A12)') 
     &     TEXT,(DSTEP(I),'%',I=1,4),INTTID
#endif
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck scrgds */
      SUBROUTINE SCRGDS(TEXT,DINTSKP,IPRINT)
C*****************************************************************************
C
C     Provide statistics about screening of gradient
C
C          DINTSKP(1,I,1) - total number of integrals, x direction
C          DINTSKP(1,I,2) - total number of integrals, y direction
C          DINTSKP(1,I,3) - total number of integrals, z direction
C
C          Step n (n = 1,2)
C          DINTSKP(n+1,I,1) - number of integrals skipped, x direction
C          DINTSKP(n+1,I,2) - number of integrals skipped, y direction
C          DINTSKP(n+1,I,3) - number of integrals skipped, z direction
C
C     Written by J. Thyssen, Jan 16, 1998
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DC = 100.0D0,D1 = 1.0D0,D0 = 0.0D0)
      DIMENSION DINTSKP(3,14,3)
      CHARACTER TEXT*6
C
      CHARACTER*7 CSTEP(7)
      CHARACTER*7 CTEMP
      CHARACTER*6 CCENT(4)
      CHARACTER*4 CTYPE(7,3)
      DIMENSION NTYPE(3)
      DIMENSION DTOTAL(2,0:3,3)
      DATA CCENT /'ONECEN','TWOCEN','THRCEN','FOUCEN'/
      DATA CTYPE /'1122','1212','1221','1222','1211','1112','1121',
     &            '1213','1223','1231','1232','1123','1233','xxxx',
     &            '1234','xxxx','xxxx','xxxx','xxxx','xxxx','xxxx'/
      DATA NTYPE /7,6,1/
C
#ifdef VAR_DEBUG
      DO I =1,3 
      DO J =1,14
      DO K =1,3
         WRITE(LUPRI,*) DINTSKP(I,J,K)
      END DO
      END DO
      END DO
#endif
      ITURN = 0
 10   CONTINUE
      IF (ITURN .EQ. 0) THEN
         WRITE(LUPRI,'(A)') 'Step 1 Screening'
      ELSE
         WRITE(LUPRI,'(A)') 'Step 1 + Step 2 Screening'
      END IF
      WRITE(LUPRI,'(A)') '* Gradient screening statistics:'
      CALL DZERO(DTOTAL,2*4*3)
      IOFF = 0
      DO ICEN = 1,3
         DO J = 1,NTYPE(ICEN)
            IOFF = IOFF + 1
            DO I = 1,3
               DO K = 1,2 
                  DTOTAL(K,0,I) = DTOTAL(K,0,I) + DINTSKP(K,IOFF,I)
                  DTOTAL(K,ICEN,I) = DTOTAL(K,ICEN,I) 
     &                             + DINTSKP(K,IOFF,I)
               END DO
            END DO
         END DO
      END DO
C
C     Total screening in X,Y and Z direction
C
      DO I = 1,3
         IF(DTOTAL(1,0,I).GT.D0) THEN
            DTEMP = DTOTAL(2,0,I) / DTOTAL(1,0,I) * DC
            WRITE(CSTEP(I),'(F6.2,A1)') DTEMP,'%'
         ELSE
            WRITE(CSTEP(I),'(A7)') '  N/A  '
         ENDIF
      END DO
      WRITE(LUPRI,'(6X,3(5X,A8))') '    X   ','    Y   ','    Z   '
      WRITE(LUPRI,'(A6,3(6X,A7))') TEXT,(CSTEP(I),I=1,3)
C
C     Screening on center type ICEN
C
c      IF (IPRINT .GE. 3) THEN
      IF (IPRINT .GE. 0) THEN
         WRITE(LUPRI,'(/A/)') 'Detailed screening statistics'
         WRITE(LUPRI,'(8X,A7,2X,A4,2X,A7,3(2X,A7))')
     &      '  %int  ','Type','  %int  ','   X   ','   Y   ','   Z   '
         IOFF = 0
         DO ICEN = 1,3
C
C           DTEMP1 is number of integrals of type ICEN
C
            DTEMP1 = DTOTAL(1,ICEN,1) +
     &               DTOTAL(1,ICEN,2) +
     &               DTOTAL(1,ICEN,3)
C
C           DTEMP2 is total number of integrals
C
            DTEMP2 = DTOTAL(1,0,1) +
     &               DTOTAL(1,0,2) +
     &               DTOTAL(1,0,3)
C
            IF (DTEMP2 .GT. D0) THEN
               WRITE(CTEMP,'(F6.2,A1)') DTEMP1/DTEMP2*DC,'%'
            ELSE
               WRITE(CTEMP,'(A7)') '  N/A  '
            END IF
            DO I = 1,3
               IF (DTOTAL(1,ICEN,I) .GT. 0) THEN
                  WRITE(CSTEP(I),'(F6.2,A1)') 
     &               DTOTAL(2,ICEN,I)/DTOTAL(1,ICEN,I)*DC,'%'
               ELSE
                  WRITE(CSTEP(I),'(A7)') '  N/A  '
               END IF
            END DO
            WRITE(LUPRI,'(A6,2X,A7,15X,3(2X,A7))') 
     &         CCENT(ICEN+1),CTEMP,(CSTEP(I),I=1,3)
C
C
C
            WRITE(LUPRI,'(/)')
            DO I = 1,NTYPE(ICEN)
               IOFF = IOFF + 1
               DTEMP2 = DINTSKP(1,IOFF,1) +
     &                  DINTSKP(1,IOFF,2) +
     &                  DINTSKP(1,IOFF,3)
               IF (DTEMP1 .GT. D0) THEN
                  WRITE(CTEMP,'(F6.2,A1)') DTEMP2/DTEMP1*DC,'%'
               ELSE
                  WRITE(CSTEP(I),'(A7)') '  N/A  '
               END IF
               DO K = 1,3
                  IF (DINTSKP(1,IOFF,K) .GT. 0) THEN
                     WRITE(CSTEP(K),'(F6.2,A1)') 
     &                  DINTSKP(2,IOFF,K)/DINTSKP(1,IOFF,K)*DC,'%'
                  ELSE
                     WRITE(CSTEP(K),'(A7)') '  N/A  '
                  END IF
               END DO
               WRITE(LUPRI,'(17X,A4,2X,A7,3(2X,A7))') 
     &            CTYPE(I,ICEN),CTEMP,(CSTEP(K),K=1,3)
            END DO
            WRITE(LUPRI,'(/)')
         END DO
      ENDIF
      ITURN = ITURN + 1
      IF (ITURN .EQ. 1) THEN
         DO I = 1,14
            DO J = 1,3
               DINTSKP(2,I,J) = DINTSKP(2,I,J) + DINTSKP(3,I,J)
            END DO
         END DO
         GOTO 10
      END IF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Gabge1 */
      SUBROUTINE GABGE1(IJOB,GABM,ITYPE,IGTYP,
     &                  MAXDIF,JATOM,NOCONT,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Generate (ij|ij) - integrals for screening purposes
C     IJOB = 0 GAB-matrix in AO-basis
C     IJOB = 1 GAB-matrix in SO-basis
C
C     Written by T.Saue Oct 11 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <aovec.h>
C
#include <blocks.h>
#include <nuclei.h>
#include <cbihr2.h>
      LOGICAL NOCONT
      DIMENSION GABM(NBASIS*NBASIS),WORK(LWORK)
#include <memint.h>
C
      NPAO = MXSHEL*MXAOVC
      CALL MEMGET('INTE',KJSTRS,NPAO*2,WORK,KFREE,LFREE)     
      CALL MEMGET('INTE',KNPRIM,NPAO*2,WORK,KFREE,LFREE)     
      CALL MEMGET('INTE',KNCONT,NPAO*2,WORK,KFREE,LFREE)     
      CALL MEMGET('INTE',KIORBS,NPAO  ,WORK,KFREE,LFREE)     
      CALL MEMGET('INTE',KJORBS,NPAO  ,WORK,KFREE,LFREE)     
      CALL MEMGET('INTE',KKORBS,NPAO  ,WORK,KFREE,LFREE)     
C
      CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            WORK(KJORBS),WORK(KKORBS),0,.FALSE.,IPRINT)
      CALL MEMREL('GABGE1',WORK,KWORK,KJORBS,KFREE,LFREE)
      CALL GABDRV(GABM,WORK(KFREE),LFREE,IJOB,ITYPE,IGTYP,
     &            WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),
     &            WORK(KIORBS),MAXDIF,JATOM,NOCONT,IPRTWO)
      CALL MEMREL('GABGE1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Gabdim */
      SUBROUTINE GABDIM(NGAB)
C
C     Return MAXSHL ...make mini-version of PAOVEC
C     T.Saue March 1997
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <shells.h>
      NGAB = KMAX
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Gabgen */
      SUBROUTINE GABGEN(IJOB,ITYPE,IGTYP,MAXDIF,JATOM,NOCONT,
     &                  IPRINT,WORK,LWORK)
C***********************************************************************
C
C     This routine will generate the reduced matrix of 
C     Cauchy-Schwarz integrals and write it to file AOPROPER
C
C     Written by T.Saue Mar 19 1997
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
C Used from
C  gnrinf.h: SRINTS
C
#include <maxorb.h>
#include <mxcent.h>
#include <gnrinf.h>
#include <nuclei.h>
#include <blocks.h>
#include <inftap.h>
      LOGICAL FNDLAB,NOCONT
      CHARACTER GABLAB(3)*8,RTNLBL(2)*8
      DIMENSION WORK(*)
#include <chrxyz.h>
      CALL QENTER('GABGEN')
#include <memint.h>
      IF(ITYPE.EQ.0) THEN
         NGAB = 1
         GABLAB(1) = 'GABXXXXX'
      ELSE IF (ITYPE .EQ.10) THEN
         NGAB = 3
         DO I = 1,NGAB
            GABLAB(I) = 'GABXXX'//CHAR(MAXDIF+48)//CHRXYZ(-I)
         END DO
      ELSE
         WRITE(LUPRI,'(A,I12)') 'GABGEN: Invalid ITYPE =',ITYPE
         CALL QUIT('GABGEN: Invalid ITYPE')
      ENDIF
      DO I = 1,NGAB
         IF    (IJOB.EQ.0) THEN
            IF (SRINTS) THEN
C           ... short range integrals
               GABLAB(I)(4:5) = 'SR'
            ELSE
               GABLAB(I)(4:5) = 'AO'
            END IF
         ELSEIF(IJOB.EQ.1) THEN
            GABLAB(I)(4:5) = 'SO'
         ELSE
            WRITE(LUPRI,'(A,I12)') 'GABGEN: Invalid IJOB =',IJOB
            CALL QUIT('GABGEN: Invalid IJOB')
         ENDIF
      END DO
      CALL GETDAT(RTNLBL(1),RTNLBL(2))     
C
C     Open file AOPROPER if not opened
C
      IF (LUPROP .LE. 0) CALL GPOPEN(LUPROP,'AOPROPER',
     &        'UNKNOWN',' ',' ',IDUMMY,.FALSE.)
C
C       Search for GAB-integrals; return if they are there
C
      DO I = 1,NGAB
         REWIND LUPROP
         IF (.NOT. FNDLAB(GABLAB(I),LUPROP)) GOTO 10
      END DO
      GOTO 999
C     ... GAB matrix(ces) was already on file
 10   CONTINUE
      REWIND LUPROP
      IF(FNDLAB('EOFLABEL',LUPROP)) BACKSPACE LUPROP
C
C     Generate GAB-integrals
C
      IF (MAXDIF.EQ.0) THEN
         NRGAB = 1
      ELSE
         NRGAB = 3
      END IF
      IF (NRGAB .NE. NGAB) CALL QUIT('GABGEN ERROR: NGAB .ne. NRGAB')
      NBUF = NRGAB*NBASIS*NBASIS
      CALL MEMGET('REAL',KBUF,NBUF,WORK,KFREE,LFREE)     
      CALL GABGE1(IJOB,WORK(KBUF),ITYPE,IGTYP,
     &            MAXDIF,JATOM,NOCONT,WORK(KFREE),LFREE,IPRINT)
C
C
C     GAB-matrix in AO-basis
C
      IF(IJOB.EQ.0) THEN
         N2GAB = NSYMBL*NSYMBL
         CALL MEMGET('REAL',KGAB,NGAB*N2GAB,WORK,KFREE,LFREE)     
         DO 123 IGAB = 1,NGAB
            IOFF1 = (IGAB-1)*NBASIS*NBASIS
            IOFF2 = (IGAB-1)*N2GAB
            CALL MKDRAO(WORK(KBUF+IOFF1),WORK(KGAB+IOFF2),1,
     &                  WORK(KFREE),LFREE,IPRINT)
            DO I = 0,N2GAB-1
               WORK(KGAB+I+IOFF2) = SQRT(ABS(WORK(KGAB+I+IOFF2)))
            ENDDO
            IF(IPRINT.GE.4) THEN
               CALL HEADER('GABGEN:Reduced GAB-matrix in AO basis',-1)
               WRITE (LUPRI,'(2A)') ' Label : ',GABLAB(IGAB)
               CALL OUTPUT(WORK(KGAB+IOFF2),1,NSYMBL,1,NSYMBL,
     &                   NSYMBL,NSYMBL,-4,LUPRI)
            ENDIF
 123     CONTINUE
      ELSEIF(IJOB.EQ.1) THEN
C
C     GAB-matrix in SO-basis
C
         IF (NGAB .GT. 1) CALL QUIT('GABGEN: only NGAB .eq. 1 for SO')
         N2GAB = MAXSHL*MAXSHL
         CALL MEMGET('REAL',KGAB,N2GAB,WORK,KFREE,LFREE)
         CALL MKDRSO(WORK(KGAB),WORK(KBUF),1,WORK(KFREE),LFREE,IPRINT)
         DO I = 0,N2GAB-1
            WORK(KGAB+I) = SQRT(ABS(WORK(KGAB+I)))
         ENDDO
         IF(IPRINT.GE.4) THEN
            CALL HEADER('GABGEN:Reduced GAB-matrix in SO basis',-1)
            WRITE (LUPRI,'(2A)') ' Label : ',GABLAB(1)
            CALL OUTPUT(WORK(KGAB),1,MAXSHL,1,MAXSHL,MAXSHL,MAXSHL,
     &                  -4,LUPRI)
         ENDIF
      ELSE
         WRITE(LUPRI,'(A,I5)') 'GABGEN: Unknown IJOB = ',IJOB
         CALL QUIT('GABGEN: Unknown IJOB !')
      ENDIF
C
C
C     Write label
C
      LEN = MAX(4,N2GAB)
      IOFF = 0
      DO I = 1,NGAB
         CALL NEWLB2(GABLAB(I),RTNLBL,LUPROP,LUPRI)
         CALL WRITT(LUPROP,LEN,WORK(KGAB+IOFF))
         IOFF = IOFF + N2GAB
      END DO
      CALL NEWLB2('EOFLABEL',RTNLBL,LUPROP,LUPRI)
C
      CALL GPCLOSE(LUPROP,'KEEP')
      CALL MEMREL('GABGEN',WORK,KWORK,KWORK,KFREE,LFREE)
C
 999  CONTINUE
      CALL QEXIT('GABGEN')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Getgab */
      SUBROUTINE GETGAB(GABM,GABLAB,N2GAB,WORK,LWORK)
C***********************************************************************
C
C     Fetch compressed matrix of Cauchy-Scwartz-integrals from file
C
C     Written by T.Saue Mar 20 1997
C
C***********************************************************************
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
      LOGICAL EXTST,FNDLAB
      CHARACTER GABLAB*8
      DIMENSION GABM(N2GAB), WORK(*)
#include <inftap.h>
C
      CALL QENTER('GETGAB')
C
C     Check if file AOPROPER is present
C
      CALL GPINQ('AOPROPER','EXIST',EXTST)
      IF(EXTST) THEN
C     File AOPROPER present; see if it is opened
        IF (LUPROP .LE. 0) CALL GPOPEN(LUPROP,
     &     'AOPROPER','OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
      ELSE
        WRITE(LUPRI,'(/A)') 'GETGAB: AOPROPER not found !'
        CALL QUIT('GETGAB: AOPROPER not found !')
      ENDIF
C
C     Get GAB-matrix
C
      REWIND LUPROP
      IF(FNDLAB(GABLAB,LUPROP)) THEN
         READ(LUPROP,ERR=10,END=20) GABM
      ELSE IF (GABLAB(6:8) .EQ. 'XXX') THEN
         WRITE(LUPRI,'(/A,A8,A)') '*** WARNING  GETGAB: ',
     &        GABLAB,' not found on AOPROPER. Regenerating.'
         CALL GABGEN(0,0,0,0,0,.FALSE.,0,WORK,LWORK)
         IF (LUPROP .LE. 0) CALL GPOPEN(LUPROP,'AOPROPER',
     &      'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
         REWIND LUPROP
         IF(FNDLAB(GABLAB,LUPROP)) THEN
            READ(LUPROP,ERR=10,END=20) GABM
         ELSE
            WRITE(LUPRI,'(/A,A8,A)')
     &           'GETGAB: Label ',GABLAB,' not found on AOPROPER !'
            CALL QUIT('GETGAB: No GAB on AOPROPER !')
         END IF
      ELSE
         WRITE(LUPRI,'(/A,A8,A)')
     &   'GETGAB: Label ',GABLAB,' not found on AOPROPER !'
         CALL QUIT('GETGAB: No derivative GAB on AOPROPER !')
      ENDIF
C
      CALL GPCLOSE(LUPROP,'KEEP')
      CALL QEXIT('GETGAB')
      RETURN
 10   CONTINUE
      WRITE(LUPRI,'(A,A8)') '* GETGAB: ERROR reading ',GABLAB
        CALL QUIT('GETGAB: Error reading GAB from AOPROPER !')
 20   CONTINUE
      WRITE(LUPRI,'(A,A8)') '* GETGAB: EOF   reading ',GABLAB
        CALL QUIT('GETGAB: EOF reading GAB from AOPROPER !')
      END
C -- end of her2gab.F --
