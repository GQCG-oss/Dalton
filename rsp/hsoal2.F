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
C  /* Deck hsoal2 */
      SUBROUTINE HSOAL2 (WORD,CMO,UDV,PV,XINDX,MJWOP,WRK,LWRK)
C
C Alternative routine for constructing spin-orbit gradient
C by way of X2HSO Case 2
C
#include <implicit.h>
      DIMENSION CMO(*),UDV(*),PV(*),XINDX(*),WRK(*)
      CHARACTER*8 WORD
C
      PARAMETER (D1=1.0D0, DM1=-1.0D0, DH=0.5D0)
      LOGICAL DO1,DO2,LORB,LCON,LREFST
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
      IF (.NOT.FLAG(27))
     &   CALL QUIT('HSOAL2 is only implemented for .DETERMINANTS')
      IF (NCONF.NE.NCREF)
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
      KGP   = 1
      KVEC  = KGP  + KZYVAR
      KFI   = KVEC + 2*IVAR
      KFREE = KFI  + N2ORBX
      LFREE = LWRK - KFREE + 1
      IF (LFREE.LT.0) CALL ERRWRK('HSOAL2',KFREE-1,LWRK)
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
     &               WRK(KFREE),LFREE,IPRRSP)
         IF (KSYMP .NE. KSYMOP) CALL QUIT(
     &      'ERROR: Unexpected symmetry of property matrix')
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
      IGRSPI = 1
C
C Operator setup
C
      IOPSYM = KSYMOP
      IOPSPI = 1
C
C Orbital gradient
C
      LORB = KZWOPT.GT.0
      CALL IPSET(0,0,1,1)
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
      CALL HSO2CR(IGRSYM,IGRSPI,WRK(KGP),WRK(KVEC),
     *            2*IVAR,NCONF,ISYMV,ISYMDN,
     *            XINDX,OVLAP,UDV,PV,WRK(KFI),WRK(KFREE),LFREE,
     *            KZYVAR,LCON,LORB,LREFST,LUMHSO,KSYMOP,
     *            IOPSPI,ISPINV,IKLVL,DUM,IDUM,DUM,IDUM,MJWOP)
      CALL HEADER('SPIN-ORBIT GRADIENT FROM HSOAL2',3)
C
      IF (IPRRSP.GT.100) THEN
         WRITE(LUPRI,'(/A)') ' Orbital part'
         CALL OUTPUT(WRK(KGP+KZCONF),1,KZWOPT,1,2,KZVAR,2,1,LUPRI)
         WRITE(LUPRI,'(/A)') ' Configuration part'
         CALL OUTPUT(WRK(KGP),1,KZCONF,1,2,KZVAR,2,1,LUPRI)
      ELSE
         CALL RSPPRW(WRK(KGP+KZCONF),MJWOP,KZWOPT,KSYMOP,KZVAR,LUPRI)
         CALL RSPPRC(WRK(KGP),KZCONF,KZVAR,LUPRI)
      END IF
      IF (KSYMOP.EQ.1 .AND. LCON) THEN
         WRITE(LUPRI,'(/A)') 'Project out reference component'
         KCREF = KVEC + IVAR
         OVLAP = DDOT(KZCONF,WRK(KCREF),1,WRK(KGP),1)
         CALL DAXPY(KZCONF,-OVLAP,WRK(KCREF),1,WRK(KGP),1)
         OVLAP = DDOT(KZCONF,WRK(KCREF),1,WRK(KGP+KZVAR),1)
         CALL DAXPY(KZCONF,-OVLAP,WRK(KCREF),1,WRK(KGP+KZVAR),1)
         CALL RSPPRW(WRK(KGP+KZCONF),MJWOP,KZWOPT,KSYMOP,KZVAR,LUPRI)
         CALL RSPPRC(WRK(KGP),KZCONF,KZVAR,LUPRI)
      END IF
      CALL QEXIT('HSOAL2')
      RETURN
      END

