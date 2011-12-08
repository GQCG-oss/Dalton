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
      SUBROUTINE CRLRV2(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,XINDX,MJWOP,
     &                  WRK,LWRK)
C
C
C   Purpose: use RSPCTL and RSPEVE to solve the eigenvalue problem
C            to get the two-index response vectors
C
C
#include <implicit.h>
#include <iratdef.h>
C
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*), FOCK(*)
C
      LOGICAL FOUND, CONV
      PARAMETER (D0 = 0.0D0)
C
#include <priunit.h>
#include <infpri.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <inflr.h>
#include <indqr.h>
#include <maxorb.h>
#include <infvar.h>
#include <qrinf.h>
      DIMENSION MJWOP(2,MAXWOP,8)
#include <inftap.h>
#include <inforb.h>
#include <infcr.h>
#include <indcr.h>
#include <infdim.h>
C
      CALL QENTER('CRLRV2')
C
C   Initialize variables
C
      THCRSP = THCLR
      MAXIT  = MAXITL
C
C   Check number of two index linear vector equations
C
      IF (NLRLB2 .LE. 0) GO TO 9999
C
C   Loop over number of two-index vectors
C   Define symmetry dependent variables
C   Allocate core memory
C   Read XY vector from LUXYVE and check norm of vector
C   Solve ( E[2] - w*S[2] ) * Nxy = XY
C
      DO 100 ILRLB2 = 1,NLRLB2
C
C   Define symmetry dependent variables
C
         KSYMOP=ISMCR2(ILRLB2,3)
         CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
         CALL SETZY(MJWOP)
C
         WRITE (LUPRI,'(//A,I3/A,A/A,A/1P,A,D15.6,/A,D15.6)')
     *   ' CRLRV2 -- linear response calc for sym:',KSYMOP,
     *   ' CRLRV2 -- operator label1: ',CRLB2(ILRLB2,1),
     *   ' CRLRV2 -- operator label2: ',CRLB2(ILRLB2,2),
     *   ' CRLRV2 -- freqr1 :',CRFRQ2(ILRLB2,1),
     *   ' CRLRV2 -- freqr2 :',CRFRQ2(ILRLB2,2)
C
C   Allocate core memory for
C   KGD two-index response vector
C   KEIVEC frequency
C
         KREDE  = 1
         KREDS  = KREDE  + MAXRM*MAXRM
         KIBTYP = KREDS  + MAXRM*MAXRM
         KEIVAL = KIBTYP + MAXRM
         KRESID = KEIVAL + 1
         KEIVEC = KRESID + 1
         KREDGD = KEIVEC + MAXRM
         KGD    = KREDGD + MAXRM
         KWRK1  = KGD    + KZYVAR
         LWRK1  = LWRK   - KWRK1
         IF (LWRK1.LT.0) CALL ERRWRK('CRLRV1 1',KWRK1-1,LWRK)
C
C   Reuse space for solution vectors.
C
         KBVECS = KGD
         KWRKE  = KWRK1
C
         WRK(KEIVAL) = CRFRQ2(ILRLB2,1) + CRFRQ2(ILRLB2,2)
C
         CALL REARSP(LURSP,KLEN,WRK(KGD),
     &               CRLB2(ILRLB2,1),CRLB2(ILRLB2,2),
     *               CRFRQ2(ILRLB2,1),CRFRQ2(ILRLB2,2),
     *               ISMCR2(ILRLB2,1),ISMCR2(ILRLB2,2),
     *               THCLR,FOUND,CONV,ANTSYM)
         IF (FOUND .AND. CONV) THEN
            WRITE(LUPRI,'(/A)')
     *           ' Solution vector already on file'
            GO TO 200
         END IF 
C
C    Read XY vector from disk
C
         CALL READDX(LUXYVE,ILRLB2,IRAT*KZYVAR,WRK(KGD))
C
         IF (IPRRSP.GT.10) THEN
            WRITE(LUPRI,'(/A)') '     XY vector in CRLRV2 '
            WRITE(LUPRI,'(A)')  '     =================== '
            WRITE (LUPRI,*) 'Column 1 = Z, Column 2 = Y'
            CALL OUTPUT(WRK(KGD),1,KZVAR,1,2,KZVAR,2,1,LUPRI)
         END IF
C
C    Check norm of XY vector
C
         GDNORM = DNRM2(KZYVAR,WRK(KGD),1)
         IF (GDNORM .LT. THRNRM) THEN
            WRITE (LUPRI,*) ' >>> RSPCTL skipped because norm of'
            WRITE (LUPRI,*) '     property vector is only',GDNORM
            CALL DZERO(WRK(KBVECS),KZYVAR)
            CALL WRTRSP(LURSP,KZYVAR,WRK(KBVECS),
     *                  CRLB2(ILRLB2,1),CRLB2(ILRLB2,2),
     *                  CRFRQ2(ILRLB2,1),CRFRQ2(ILRLB2,2),
     *                  ISMCR2(ILRLB2,1),ISMCR2(ILRLB2,2),D0,D0)
            GO TO 100
         END IF
C
C    Solve ( E[2] - w*S[2] ) * Nxy = XY
C    If INVEXP the linear equation is solved by matrix inversion
C
         KZRED  = 0
         KZYRED = 0
         KEXSIM = 1
         KEXCNV = KEXSIM
         CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *               .TRUE.,CRLB2(ILRLB2,1),CRLB2(ILRLB2,2),WRK(KGD),
     *               WRK(KREDGD),WRK(KREDE),WRK(KREDS),WRK(KIBTYP),
     *               WRK(KEIVAL),WRK(KRESID),WRK(KEIVEC),XINDX,
     *               WRK(KWRK1),LWRK1)
         CALL RSPEVE(WRK(KIBTYP),WRK(KEIVAL),WRK(KEIVEC),
     *                     WRK(KBVECS),WRK(KWRKE),1,0)
C
         IF (INVEXP) THEN
            WRITE(LUPRI,'(/A)') '*** GETNXY, INVEXP = TRUE'
            CALL GETNXY(WRK(KEIVAL),WRK(KGD),
     *                  CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
         END IF 
C
         CALL WRTRSP(LURSP,KZYVAR,WRK(KBVECS),
     *               CRLB2(ILRLB2,1),CRLB2(ILRLB2,2),
     *               CRFRQ2(ILRLB2,1),CRFRQ2(ILRLB2,2),
     *               ISMCR2(ILRLB2,1),ISMCR2(ILRLB2,2),WRK(KRESID),D0)
         VAL = DNRM2(KZYVAR,WRK(KBVECS),1)
C
  200    WRITE(LUPRI,'(2A12,2A10,A20,/A42,A22)')
     *   CRLB2(ILRLB2,1),CRLB2(ILRLB2,2),'freq1','freq2','Norm',
     *   ' -----------------------------------------',
     *   '----------------------'
         WRITE(LUPRI,'(F34.6,F10.6,F20.8)') 
     *          CRFRQ2(ILRLB2,1),CRFRQ2(ILRLB2,2),VAL
         IF (IPRRSP.GT.10 .AND. (.NOT.FOUND)) THEN
            WRITE(LUPRI,'(/A)') '     Final result in CRVEC2 '
            WRITE(LUPRI,'(A)')  '     ====================== '
            WRITE (LUPRI,*) 'Column 1 = Z, Column 2 = Y'
            CALL OUTPUT(WRK(KBVECS),1,KZVAR,1,2,KZVAR,2,1,LUPRI)
         END IF
  100 CONTINUE
C
C *** end of CRLRV2 --
C
 9999 CALL QEXIT('CRLRV2')
      RETURN
      END
