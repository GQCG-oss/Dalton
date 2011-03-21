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
#include "single.h"
C /* Deck gdrv */
      SUBROUTINE GDRV(WORK,LWORK)
#include <implicit.h>
      INTEGER LWORK
      REAL*8 WORK(LWORK)
#include <infdim.h>
#include <inforb.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include "rspprp.h"
#include <inflr.h> 
#include <inflin.h>
#include <dummy.h>
#include <inftap.h>
#include <maxorb.h>
#include <infinp.h>
#include <gtensor.h>
      INTEGER KFREE, LFREE, LWRK1, KWRK1
      INTEGER KINDX, KCMO, KUDV, KPVX, KFOCK, KFC, KFV, KFCAC, KH2AC
      INTEGER ISYM
C
      CALL QENTER('GDRV')
C
C Define variables for linear response
C
      CALL GINIT
C
      KFREE=1
      LFREE=LWORK
      CALL MEMGET('REAL',KINDX,LCINDX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCMO,NCMOT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCREF,NCREF,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP,NNASHX,WORK,KFREE,LFREE)   
      CALL MEMGET('REAL',KUDV,N2ASHX,WORK,KFREE,LFREE)
      IF (RSPCI) THEN
         CALL MEMGET('REAL',KPVX,0,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KFOCK,0,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KFC,0,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KFV,0,WORK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KPVX,2*LPVMAT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KFOCK,N2ORBT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KFC,NNORBT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KFV,NNORBT,WORK,KFREE,LFREE)
      END IF
      CALL MEMGET('REAL',KFCAC,NNASHX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2AC,NNASHX*NNASHX,WORK,KFREE,LFREE)
      CALL MEMGET('WORK',KWRK1,LWRK1,WORK,KFREE,LFREE)
C
C Read interface information
C
      CALL RSPMC(WORK(KCMO),WORK(KUDV),WORK(KPVX),WORK(KFOCK),WORK(KFC),
     *           WORK(KFV),WORK(KFCAC),WORK(KH2AC),WORK(KINDX),
     *           WORK(KWRK1),LWRK1)
C
      IF (.NOT.TDHF .AND. .NOT.RSPCI) THEN
         CALL GETREF(WORK(KCREF),NCREF)
         CALL RSPDM(IREFSY,IREFSY,NCREF,NCREF,WORK(KCREF),WORK(KCREF),
     &      WORK(KUDV),WORK(KPVX),0,0,.FALSE.,.FALSE.,WORK(KINDX),
     &      WORK(KWRK1),1,LWRK1)
         CALL RSPDM(IREFSY,IREFSY,NCREF,NCREF,WORK(KCREF),WORK(KCREF),
     &     WORK(KUDV),WORK(KPVX+LPVMAT),1,1,.FALSE.,.FALSE.,WORK(KINDX),
     &      WORK(KWRK1),1,LWRK1)
      END IF
C
C Calculate first-order properties
C
      CALL GSET(WORK(KCMO),WORK(KUDV),WORK(KINDX),WORK(KWRK1),LWRK1)
      CALL MEMCHK('GDRV<-GSET',WORK,1)
C
C Calculate second-order properties
C
      DO ISYM = 1,NSYM
         KSYMOP = ISYM
C
C Define variables that depend on symmetry
C
         CALL RSPVAR(WORK(KUDV),WORK(KFOCK),WORK(KFC),WORK(KFV),
     &      WORK(KFCAC),WORK(KH2AC),WORK(KINDX),WORK(KWRK1),LWRK1)
         CALL MEMCHK('GDRV<-RSPVAR',WORK,1)
         IF (KZVAR.EQ.0) GO TO 100
C
C Solve linear response equations
C
         IF (NGPLR(ISYM).GT.0) THEN
            CALL RSPLR(
     &         WORK(KCMO),WORK(KUDV),WORK(KPVX),WORK(KFC),WORK(KFV),
     &         WORK(KFCAC),WORK(KH2AC),WORK(KINDX),WORK(KWRK1),LWRK1)
            CALL MEMCHK('GDRV<-RSPLR',WORK,1)
            CALL GSET2(WORK(KCMO),WORK(KUDV),WORK(KPVX),WORK(KINDX),
     &         WORK(KWRK1),LWRK1)
         END IF
  100 CONTINUE
      END DO
C
      CALL GPRINT(WORK(KWRK1),LWRK1)
C
C Clear LR arrays to avoid duplicate calculations
C
      CALL IZERO(NGPLR,8)
C
      CALL MEMREL('GDRV',WORK,1,1,KFREE,LFREE)
      CALL QEXIT('GDRV')
      RETURN
      END
C
C Deck /* gset2 */
C
      SUBROUTINE GSET2(CMO,UDV,PV,XINDX,WORK,LWORK)
      IMPLICIT NONE
      REAL*8 CMO(*), UDV(*), PV(*), XINDX(*), WORK(*)
      INTEGER LWORK, MAXLBL, NPRLBL
C
      REAL*8 DDOT, DNRM2
      EXTERNAL DDOT, DNRM2
      INTEGER INDPRP
      EXTERNAL INDPRP
#include <gtensor.h>
#include <qrinf.h>
#include <inforb.h>
#include "rspprp.h"
#include <inflr.h>
#include <dummy.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <symmet.h>
#include <wrkrsp.h>
      INTEGER KFREE, LFREE
      LOGICAL LDUMMY, FOUND, CONV, OPENED
      INTEGER I, J, KLEN, LU, KANG, ISYM, IPRP, ISO
      REAL*8 RESID, ANTSYM, D0, DNORM
      PARAMETER (D0=0.0D0)
      CHARACTER*8 BLANK
      DATA BLANK/'        '/
      CALL QENTER('GSET2')
      INQUIRE(FILE='RSPVEC',OPENED=OPENED,NUMBER=LU)
      write(2,*)'opened with lu',opened,lu
      IF (.NOT. OPENED) THEN
         LU=-1
         CALL GPOPEN(
     &      LU,'RSPVEC','OLD','SEQUENTIAL','UNFORMATTED',0,LDUMMY
     &      )
      END IF
      REWIND(LU)
      KFREE=1
      LFREE=LWORK
      CALL MEMGET('REAL',KANG,KZYVAR,WORK,KFREE,LFREE)
      DO IPRP=1,NGPLR(KSYMOP)
         DO I=1,3
            IF (ANGLBL(I).EQ.LBLLR(KSYMOP,IPRP)
     &         .AND. LLROP(INDPRP(ANGLBL(I)))) THEN
               CALL REARSP(LU,KLEN,WORK(KANG),ANGLBL(I),BLANK,
     &            D0,DUMMY,KSYMOP, IDUMMY, THCLR,FOUND,CONV,ANTSYM)
               IF (KLEN.NE.KZYVAR) THEN
                  CALL QUIT('GSET2:Serious error reading RSPVEC')
               END IF
               DNORM = DNRM2(KLEN,WORK(KANG),1)
               IF (FOUND .AND. CONV .AND. DNORM.GT.THRNRM) THEN
                  DO ISO=1,NSO
                     DO J=1,3
                        IF (KSYMOP.EQ.(ISYMAX(J,2)+1)) THEN
                           CALL GETGPV(
     &                     SOLBL(J,ISO),DUMMY,DUMMY,CMO,UDV,PV,
     &                     XINDX,ANTSYM,WORK(KFREE),LFREE
     &                     )
                           GOZSO(I,J,ISO) = DDOT(
     &                     KLEN,WORK(KANG),1,WORK(KFREE),1)
                           GOZSO(I,J,ISO)=GOZSO(I,J,ISO)/REFSPIN
                        END IF
                     END DO
                  END DO
               END IF
            END IF
         END DO
      END DO
C
C Restore original state of file
C
      IF (.NOT.OPENED) CALL GPCLOSE(LU,'KEEP')
      CALL MEMREL('GSET2',WORK,KANG,KANG,KFREE,LFREE)
      CALL QEXIT('GSET2')
      RETURN
      END

C /* Deck ginit */
      SUBROUTINE GINIT
#include <implicit.h>
      DIMENSION AVER(3)
#include <gtensor.h>
#include <maxorb.h>
#include <infinp.h>
#include "rspprp.h"
#include <inflr.h>
#include "infave.h"
#include <priunit.h>
      LOGICAL FNDLB2
      EXTERNAL FNDLB2
      LQRGP = .TRUE.
      REFSPIN=DBLE(ISPIN-1)/2
      IF (.NOT.(ADDSO.OR.DOOZSO1.OR.DOOZSO2)) RETURN
C
C Options
C 1. (Default) One and two-electron spin-orbit separate
C 2. If either is specified in input only that one is calculated
C 3. If  'MEAN FIELD' or 'SCALED' is specified the two-electron part is replaced
C
      IF (ECC) THEN
         ANGLBL(1) = 'XANGECC '
         ANGLBL(2) = 'YANGECC '
         ANGLBL(3) = 'ZANGECC '
C
C Make RSPSYM calculate dipole integrals if mising
C
         I=INDPRP('XDIPLEN ')
      ELSE
         ANGLBL(1) = 'XANGMOM '
         ANGLBL(2) = 'YANGMOM '
         ANGLBL(3) = 'ZANGMOM '
      END IF
      IF (ADDSO) THEN
         IF (MNFSO) THEN
            SOLBL(1,1) = 'X1MNF-SO' 
            SOLBL(2,1) = 'Y1MNF-SO' 
            SOLBL(3,1) = 'Z1MNF-SO'
         ELSE IF (SCALED_CHARGES) THEN
            SOLBL(1,1) = 'X1SPNSCA'
            SOLBL(2,1) = 'Y1SPNSCA'
            SOLBL(3,1) = 'Z1SPNSCA'
         ELSE
            SOLBL(1,1) = 'X SPNORB'
            SOLBL(2,1) = 'Y SPNORB'
            SOLBL(3,1) = 'Z SPNORB'
         END IF
      ELSE
         IF (DOOZSO1) THEN
            SOLBL(1,1) = 'X1SPNORB'
            SOLBL(2,1) = 'Y1SPNORB'
            SOLBL(3,1) = 'Z1SPNORB'
         END IF
         IF (DOOZSO2) THEN
            IF (MNFSO) THEN
               SOLBL(1,2) = 'X1MNF-SO' ! contains x1spnorb - subtracted
               SOLBL(2,2) = 'Y1MNF-SO' ! in hso
               SOLBL(3,2) = 'Z1MNF-SO'
            ELSE IF (SCALED_CHARGES) THEN
               SOLBL(1,2) = 'X1SPNSCA'
               SOLBL(2,2) = 'Y1SPNSCA'
               SOLBL(3,2) = 'Z1SPNSCA'
            ELSE IF (SPLITSO2) THEN
               SOLBL(1,2) = 'X2SPNOWN'
               SOLBL(2,2) = 'Y2SPNOWN'
               SOLBL(3,2) = 'Z2SPNOWN'
               SOLBL(1,3) = 'X2SPNOTH'
               SOLBL(2,3) = 'Y2SPNOTH'
               SOLBL(3,3) = 'Z2SPNOTH'
            ELSE
               SOLBL(1,2) = 'X2SPNORB'
               SOLBL(2,2) = 'Y2SPNORB'
               SOLBL(3,2) = 'Z2SPNORB'
            END IF
         END IF
      END IF
      CALL DZERO(GTENSOR,9)
      CALL DZERO(GC1,9)
      CALL DZERO(GC2,9)
      CALL DZERO(GOZSO,9*NSO)
C
C Orbital-Zeeman + Spin-orbit contribution
C
      DO IX=1,3
         LLROP( INDPRP(ANGLBL(IX))) = .TRUE.
      END DO
      NSO=0
      IF (DOOZSO1) NSO=NSO+1
      IF (DOOZSO2) NSO=NSO+1
      IF (SPLITSO2) NSO=NSO+1
      IF (ADDSO) NSO=1
      DO IX=1,3
         DO J=1,NSO
            LLROP( INDPRP(SOLBL(IX,J))) = .FALSE.
            LRRANK( INDPRP(SOLBL(IX,J))) = 1
         END DO
      END DO
      IF (.NOT.G_ALL) THEN
         IF (INDEX(G_LINE,"X").GT.0) THEN
            LLROP(INDPRP(ANGLBL(1))) = .FALSE.
            LRRANK(INDPRP(ANGLBL(1))) = 0
            DO J=1,NSO
C              LLROP( INDPRP(SOLBL(1,J))) = .FALSE.
C              LRRANK( INDPRP(SOLBL(1,J))) = 0
            END DO
         END IF
         IF (INDEX(G_LINE,"Y").GT.0) THEN
            LLROP(INDPRP(ANGLBL(2))) = .FALSE.
            LRRANK(INDPRP(ANGLBL(2))) = 0
            DO J=1,NSO
C              LLROP( INDPRP(SOLBL(2,J))) = .FALSE.
C              LRRANK( INDPRP(SOLBL(2,J))) = 0
            END DO
         END IF
         IF (INDEX(G_LINE,"Z").GT.0) THEN
            LLROP(INDPRP(ANGLBL(3))) = .FALSE.
            LRRANK(INDPRP(ANGLBL(3))) = 0
            DO J=1,NSO
C              LLROP( INDPRP(SOLBL(3,J))) = .FALSE.
C              LRRANK( INDPRP(SOLBL(3,J))) = 0
            END DO
         END IF
      END IF
      RETURN
      END
C /* Deck gset */
      SUBROUTINE GSET(CMO,UDV,XINDX,WRK,LWRK)
#include <implicit.h>
      DIMENSION CMO(*), UDV(*), XINDX(*), WRK(LWRK)
#include <gtensor.h>
#include <inftap.h>
#include <mxcent.h>
#include <cbihr1.h>
#include <cbiher.h>
#include <orgcom.h>
      LOGICAL FNDLAB,FNDLB2,PRINTF
      EXTERNAL FNDLAB,FNDLB2
      CHARACTER*7 RTNLBL(2)
      DOUBLE PRECISION GAGSAV(3)
      IF (.NOT.GCALC) RETURN
      CALL QENTER('GSET')
C
C  Recalculate angular momentum for electronic charge centroid (ECC)
C
      IF (ECC) THEN
         CALL DCOPY(3,GAGORG,1,GAGSAV,1)
         CALL DCOPY(3,ECCORG,1,GAGORG,1)
      END IF
C
C  Relativistic mass correction
C
      IF (DORMC) THEN
         CALL GETRMC(CMO,XINDX,RMC,WRK,LWRK)
      END IF 
C
C  One-electron diamagnetic term ("gauge-correction")
C
      IF (DOGC1) THEN
          CALL GETGC1(CMO,XINDX,WRK,LWRK)
      END IF
C
C  Two-electron diamagnetic term ("gauge-correction")
C
      IF (DOGC2) THEN
          IF (SCALED_CHARGES) THEN 
              CALL GETGC2SCA(CMO,XINDX,WRK,LWRK)
          ELSE    
              CALL AVEDSO(CMO,XINDX,WRK,LWRK)
          END IF 
      END IF
      IF (ECC) CALL DCOPY(3,GAGSAV,1,GAGORG,1)

      CALL QEXIT('GSET')
      RETURN
      END
C
      SUBROUTINE SETECC(WORK,KFREE,LFREE)
#include <implicit.h>
#include <inforb.h>
#include <inftap.h>
#include <gtensor.h>
#include <priunit.h>
      DIMENSION WORK(*)
      INTEGER KFREE, LFREE  
      CHARACTER*8 RTNLBL(2)
      INTEGER N, KD, KP
C       
      KFRSAV = KFREE  
      LFRSAV = LFREE 
      CALL MEMGET('REAL',KD,N2BASX,WORK,KFREE,LFREE)      
      CALL MEMGET('REAL',KP,N2BASX,WORK,KFREE,LFREE)  
C
C Get total density
C
      CALL GET_TD(WORK(KD),WORK,KFREE,LFREE)
C
C Get dipole
C
      CALL GET_P('XDIPLEN ',WORK(KP))
      ECCORG(1)=DDOT(N2BASX,WORK(KP),1,WORK(KD),1)
      CALL GET_P('YDIPLEN ',WORK(KP))
      ECCORG(2)=DDOT(N2BASX,WORK(KP),1,WORK(KD),1)
      CALL GET_P('ZDIPLEN ',WORK(KP))
      ECCORG(3)=DDOT(N2BASX,WORK(KP),1,WORK(KD),1)
      N=2*NISHT+NASHT
      CALL HEADER('OUTPUT FROM SETECC',2)
      WRITE(LUPRI,'(/A,I5    )') '   Number of electrons:',N
      WRITE(LUPRI,'(A,3F10.6)') '   Dipole moment  :',ECCORG
      ECCORG(1)=ECCORG(1)/N
      ECCORG(2)=ECCORG(2)/N
      ECCORG(3)=ECCORG(3)/N
      WRITE(LUPRI,'(A,3F10.6)') '   Charge centroid:',ECCORG
      CALL MEMREL('SETECC',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END

      SUBROUTINE GET_TD(D,WORK,KFREE,LFREE)
#include <implicit.h>   
      DIMENSION WORK(*), D(*)
      INTEGER KFREE, LFREE   
#include <inforb.h>
#include <inftap.h>
#include <dummy.h>
#include <priunit.h>
      PARAMETER ( D1 = 1.0D0)
      LOGICAL CLOSED
C      
      CALL QENTER('GET_D')
C
      KFRSAV = KFREE  
      LFRSAV = LFREE 
      CALL MEMGET('REAL',KCMO,NCMOT,WORK,KFREE,LFREE)      
      CALL MEMGET('REAL',KTMP,NNASHX,WORK,KFREE,LFREE) 
      CALL MEMGET('REAL',KUDV,N2ASHX,WORK,KFREE,LFREE)  
      CALL MEMGET('REAL',KDV,N2BASX,WORK,KFREE,LFREE) 
C
      CLOSED=LUSIFC.LT.0
      IF (CLOSED) CALL GPOPEN(
     &   LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,.FALSE. 
     &   )
      REWIND(LUSIFC)
      CALL MOLLAB(LBSIFC,LUSIFC,LUPRI)
      READ(LUSIFC)
      READ(LUSIFC)
      CALL READT(LUSIFC,NCMOT,WORK(KCMO))
      READ(LUSIFC)
      CALL READT(LUSIFC,NNASHX,WORK(KTMP))
      CALL DSPTSI(NASHT,WORK(KTMP),WORK(KUDV))
      CALL GTDMSO(WORK(KUDV),WORK(KCMO),D,WORK(KDV),WORK(KFREE))
      CALL DAXPY(N2BASX,D1,WORK(KDV),1,D,1)

      IF (CLOSED) CALL GPCLOSE(LUSIFC,'KEEP')
C
      CALL MEMREL('GET_TD',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
C       
      CALL QEXIT('GET_D')
      RETURN
      END

      SUBROUTINE ECCSET(CMO,UDV,WRK,LWRK)
#include <implicit.h>
      DIMENSION CMO(*), UDV(*), WRK(LWRK)
#include <maxorb.h>
#include <priunit.h>
#include <gtensor.h>
#include <inforb.h>
#include <infinp.h>
#include <orgcom.h>
#include <inftap.h>
#include <infpri.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <symmet.h>
#include <wrkrsp.h>
#include <iratdef.h>
#include <infrsp.h>
      LOGICAL FNDLB2
      CHARACTER*7 RTNLBL(2)
      DIMENSION AVER(3)
      CALL QENTER('ECCSET')
C
C  Determine electronic charge centroid
C
C  Find dipole integrals
C
      !IPRRSP=IPRLR
      IPRG=IPRRSP
      CALL GPOPEN(LUPROP,'AOPROPER','OLD','SEQUENTIAL','UNFORMATTED',
     &   0,.FALSE.)
      IF (.NOT.FNDLB2('XDIPLEN ',RTNLBL,LUPROP)) THEN
         CALL PR1INT('DIPLEN ',WRK,LWRK,IDUMMY,
     &               IDUMMY,.FALSE.,.FALSE.,0)
      END IF
      CALL GPCLOSE(LUPROP,'KEEP')
C
C  Calculate average
C
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KPRPMO,N2ORBX,WRK,KFREE,LFREE)
      CALL DZERO(AVER,3)
C
C x
C
      KSYMP = ISYMAX(1,1)+1
      IF (KSYMP.EQ.1) THEN
         CALL PRPGET(
     &      'XDIPLEN ',CMO,WRK(KPRPMO),KSYMP,ANTSYM,
     &      WRK(KFREE),LFREE,IPRG
     &      )
         CALL PRPONE(WRK(KPRPMO),UDV,AVER(1),0,'XDIPLEN')
      END IF
C
C y
C
      KSYMP = ISYMAX(2,1)+1
      IF (KSYMP.EQ.1) THEN
         CALL PRPGET(
     &      'YDIPLEN ',CMO,WRK(KPRPMO),KSYMP,ANTSYM,
     &      WRK(KFREE),LFREE,IPRG
     &      )
         CALL PRPONE(WRK(KPRPMO),UDV,AVER(2),0,'YDIPLEN')
      END IF
C
C z
C
      KSYMP = ISYMAX(3,1)+1
      IF (KSYMP.EQ.1) THEN
         CALL PRPGET(
     &      'ZDIPLEN ',CMO,WRK(KPRPMO),KSYMP,ANTSYM,
     &      WRK(KFREE),LFREE,IPRG
     &      )
         CALL PRPONE(WRK(KPRPMO),UDV,AVER(3),0,'ZDIPLEN')
      END IF
C
C  Electronic charge cenroid
C
      N=2*NISHT+NACTEL
      GAGORG(1)=AVER(1)/N
      GAGORG(2)=AVER(2)/N
      GAGORG(3)=AVER(3)/N
      CALL HEADER('OUTPUT FROM ECCSET',2)
         WRITE(LUPRI,'(/A,I5    )') '   Number of electrons:',N
         WRITE(LUPRI,'(A,3F10.6)') '   Dipole moment  :',AVER
         WRITE(LUPRI,'(A,3F10.6)') '   Charge centroid:',GAGORG
C
C  Recalculate gauge-dependent integrals.
C
      LUPROP=-1
      CALL GPOPEN(LUPROP,'AOPROPER','OLD','SEQUENTIAL','UNFORMATTED',
     &   0,.FALSE.)
      CALL PR1INT('ANGMOM ',WRK,LWRK,IDUMMY,IDUMMY,.FALSE.,.FALSE.,0)
      CALL PR1INT('ELGDIAN',WRK,LWRK,IDUMMY,IDUMMY,.FALSE.,.FALSE.,0)
      CALL GPCLOSE(LUPROP,'KEEP')
C
      CALL QEXIT('ECCSET')
      RETURN
      END
C
      SUBROUTINE GINP(WORD)
C
#include <implicit.h>
C
#include <infrsp.h>
#include <wrkrsp.h>
#include "infs0.h"
#include <infpri.h>
#include <infdim.h>
#include <trhso.h>
#include <infhso.h>
#include <gtensor.h>
#include <priunit.h>
C
      LOGICAL NEWDEF
      PARAMETER ( NTABLE =  12 )
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      LOGICAL NOGC2
      DATA NOGC2/.FALSE./
C
      DATA TABLE /'.ZERO  ', '.ECC   ', '.RMC   ', '.GC1   ', '.GC2   ',
     &            '.OZSO1 ', '.OZSO2 ', '.MEAN-F', '.TEST A', '.ADD-SO',
     &            '.SCALED', '.OWN-OT'/
      NEWDEF = (WORD .EQ. '.G-TENS')
C
C Defaults
C
      GCALC = .TRUE.
      G_ALL = .TRUE.
C
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, *) WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') GO TO 100
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO I=1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12) , I
                  END IF
               END DO
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
C
C When called as sub input from other input routine (now ESRINP)
C
               BACKSPACE LUCMD
               RETURN
C              WRITE (LUPRI,'(/,3A,/)') ' KEYWORD "',WORD,
C    *            '" NOT RECOGNIZED IN GINP.'
C              CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
C              CALL QUIT(' ILLEGAL KEYWORD IN GINP ')
 1             CONTINUE
                  G_ALL = .FALSE.
                  READ(LUCMD,'(A80)')G_LINE
               GO TO 100
 2             CONTINUE
                  ECC = .TRUE.
               GO TO 100
 3             CONTINUE
                  DORMC = .TRUE.
                  DOALL = .FALSE.
               GO TO 100
 4             CONTINUE
                  DOGC1 = .TRUE.
                  DOALL = .FALSE.
               GO TO 100
 5             CONTINUE
                  DOGC2 = .TRUE.
                  DOALL = .FALSE.
               GO TO 100
 6             CONTINUE
                  DOOZSO1 = .TRUE.
                  DOALL = .FALSE.
               GO TO 100
 7             CONTINUE
                  DOOZSO2 = .TRUE.
                  DOALL = .FALSE.
               GO TO 100
 8             CONTINUE
                  MNFSO = .TRUE.
               GO TO 100
 9             CONTINUE
                  MNFSO = .TRUE.
                  TEST_AMFI = .TRUE.
               GO TO 100
10             CONTINUE
                  ADDSO = .TRUE.
               GO TO 100
11             CONTINUE
                  SCALED_CHARGES = .TRUE.
               GO TO 100
12             CONTINUE
                  SPLITSO2 = .TRUE.
               GO TO 100
            ELSE IF (PROMPT.EQ.'*') THEN
               BACKSPACE LUCMD
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' PROMPT "',WORD,
     *            '" NOT RECOGNIZED IN GINP  .'
               CALL QUIT(' ILLEGAL PROMPT IN GINP ')
            END IF
         GO TO 100
      END IF
  300 CONTINUE
      IF (DOALL) THEN
         DORMC = .TRUE.
         DOGC1 = .TRUE.
         DOGC2 = .TRUE.
         DOOZSO1 = .TRUE.
         DOOZSO2 = .TRUE.
      END IF
C     DOOZSO=DOOZSO1.OR.DOOZSO2
      CALL GINIT
      RETURN
      END
C
C *** END OF GINP
C
      SUBROUTINE GETRMC(CMO,XINDX,AVE,WRK,LWRK)
#include <implicit.h>
      DIMENSION CMO(*),XINDX(*),WRK(*)
#include <gtensor.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <codata.h>
#include <gfac.h>
#include <dummy.h>
#include <maxorb.h>
#include <infinp.h>
C
      LOGICAL OLDDX
      LOGICAL FNDLB2
      EXTERNAL FNDLB2
      CHARACTER*7 RTNLBL(2)
      CALL QENTER('GETRMC')
C
C  Check that kinetic energy integrals exist
C
      LUPROP = -1
      CALL GPOPEN(LUPROP,'AOPROPER','OLD','SEQUENTIAL','UNFORMATTED',
     &   IDUMMY,OLDDX)
      REWIND LUPROP
      IF (.NOT.FNDLB2('KINENERG',RTNLBL,LUPROP)) THEN
         CALL PR1INT('KINENER',WRK,LWRK,IDUMMY,
     &               IDUMMY,.FALSE.,.FALSE.,0)
      END IF
      CALL GPCLOSE(LUPROP,'KEEP')
C
C Get kinetic energy integrals
C
      KFREE=1
      LFREE=LWRK
      CALL MEMGET('REAL',KPRPMO,N2ORBX,WRK,KFREE,LFREE)
      KSYMOP=1
      IPRG=IPRRSP
      CALL PRPGET(
     &   'KINENERG',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &   WRK(KFREE),LFREE,IPRG
     &   ) 
C
C
C Get triplet density
C
      CALL MEMGET('REAL',KDEN1,N2ASHX,WRK,KFREE,LFREE)
      IF (HSROHF) THEN  
         CALL DUNIT(WRK(KDEN1),NASHT)
      ELSE 
         CALL MEMGET('REAL',KCREF,NCREF,WRK,KFREE,LFREE)
         CALL GETREF(WRK(KCREF),NCREF)
         CALL RSPDM(IREFSY,IREFSY,NCREF,NCREF,WRK(KCREF),WRK(KCREF),
     &           WRK(KDEN1),DUMMY,
     &           1,0,.FALSE.,.TRUE.,XINDX,WRK,KFREE,LFREE)
      END IF 
        
C     CALL RSPDM(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR, RHO1,RHO2,
C    *                 ISPIN1,ISPIN2,TDM,NORHO2,XINDX,WORK,
C    *                 KFREE,LFREE)
C
C Expectation value
C
      TRPLET = .TRUE.
      CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GRELMC,10,'KINENERG')
      TRPLET = .FALSE. 
C     CALL PRPONE(PRPMO,UDV,ONETOT,IPRONE,PRPLBL)
C
C Scale
C
      GRELMC=-GRELMC*GFAC*ALPHA2/(2*REFSPIN)
C
C Finish
C
      CALL MEMREL('GETRMC',WRK,1,1,KFREE,LFREE)
      CALL QEXIT('GETRMC')
      RETURN
      END
      SUBROUTINE GETGC1(CMO,XINDX,WRK,LWRK)
#include <implicit.h>
      DIMENSION CMO(*),XINDX(*),WRK(*)
#include <gtensor.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <codata.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <infinp.h>
#include <mxcent.h>
#include <symmet.h>
#include <dummy.h>
      LOGICAL FNDLB2
      LOGICAL OLDDX
      CHARACTER*7 RTNLBL(2)
#include <ibtfun.h>
C
      CALL QENTER('GETGC1')
C
C  Check that diamagnetic spin-orbit integrals exist
C
      LUPROP = -1
      CALL GPOPEN(LUPROP,'AOPROPER','OLD','SEQUENTIAL','UNFORMATTED',
     &   IDUMMY,OLDDX)
      REWIND LUPROP
      IF (.NOT.FNDLB2('D1-SO XX',RTNLBL,LUPROP)) THEN
         CALL PR1INT('ELGDIAN',WRK,LWRK,IDUMMY,
     &               IDUMMY,.FALSE.,.FALSE.,0)
      END IF
      CALL GPCLOSE(LUPROP,'KEEP')
C
C Get triplet density
C
      KFREE=1
      LFREE=LWRK
      CALL MEMGET('REAL',KDEN1,N2ASHX,WRK,KFREE,LFREE)
       IF (HSROHF) THEN  
         CALL DUNIT(WRK(KDEN1),NASHT)
      ELSE   
         CALL MEMGET('REAL',KCREF,NCREF,WRK,KFREE,LFREE)
         CALL GETREF(WRK(KCREF),NCREF)
         CALL RSPDM(IREFSY,IREFSY,NCREF,NCREF,WRK(KCREF),WRK(KCREF),
     &           WRK(KDEN1),DUMMY,
     &           1,0,.FALSE.,.TRUE.,XINDX,WRK,KFREE,LFREE)
      END IF 
C     CALL RSPDM(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR, RHO1,RHO2,
C    *                 ISPIN1,ISPIN2,TDM,NORHO2,XINDX,WORK,
C    *                 KFREE,LFREE)
C
C Get diamagnetic spin-orbit integrals
C
      CALL MEMGET('REAL',KPRPMO,N2ORBX,WRK,KFREE,LFREE)
C
C Expectation value
C
      TRPLET = .TRUE.
      KSYMOP=1
      IPRG=IPRRSP
      CALL PRPGET(
     &   'D1-SO XX',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &   WRK(KFREE),LFREE,IPRG
     &   )
      CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GC1(1,1),10,'D1-SO XX')
      CALL PRPGET(
     &   'D1-SO YY',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &   WRK(KFREE),LFREE,IPRG
     &   )
      CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GC1(2,2),10,'D1-SO YY')
      CALL PRPGET(
     &   'D1-SO ZZ',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &   WRK(KFREE),LFREE,IPRG
     &   )
      CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GC1(3,3),10,'D1-SO ZZ')
      IF (IBTXOR(ISYMAX(1,1),ISYMAX(2,1)).EQ.0) THEN
         CALL PRPGET(
     &      'D1-SO XY',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &      WRK(KFREE),LFREE,IPRG
     &      )
         CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GC1(1,2),10,'D1-SO XY')
         CALL PRPGET(
     &      'D1-SO YX',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &      WRK(KFREE),LFREE,IPRG
     &      )
         CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GC1(2,1),10,'D1-SO YX')
      END IF
      IF (IBTXOR(ISYMAX(1,1),ISYMAX(3,1)).EQ.0) THEN
         CALL PRPGET(
     &      'D1-SO XZ',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &      WRK(KFREE),LFREE,IPRG
     &      )
         CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GC1(1,3),10,'D1-SO XZ')
         CALL PRPGET(
     &      'D1-SO ZX',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &      WRK(KFREE),LFREE,IPRG
     &      )
         CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GC1(3,1),10,'D1-SO ZX')
      END IF
      IF (IBTXOR(ISYMAX(2,1),ISYMAX(3,1)).EQ.0) THEN
         CALL PRPGET(
     &      'D1-SO YZ',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &      WRK(KFREE),LFREE,IPRG
     &      )
         CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GC1(2,3),10,'D1-SO YZ')
         CALL PRPGET(
     &      'D1-SO ZY',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &      WRK(KFREE),LFREE,IPRG
     &      )
         CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GC1(3,2),10,'D1-SO ZY')
      END IF
      TRPLET = .FALSE.
C
C Scale
C
      CALL DSCAL(9,ALPHA2/(2*REFSPIN),GC1,1)
C
C Finish
C
      CALL MEMREL('GETGC1',WRK,1,1,KFREE,LFREE)
      CALL QEXIT('GETGC1')
      RETURN
      END
C Deck Scaled effective gauge correction 
      SUBROUTINE GETGC2SCA(CMO,XINDX,WRK,LWRK)
#include <implicit.h>
      DIMENSION CMO(*),XINDX(*),WRK(*)
#include <gtensor.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <codata.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <symmet.h>
#include <dummy.h>
#include <infinp.h>
      LOGICAL FNDLB2
      LOGICAL OLDDX
      CHARACTER*7 RTNLBL(2)
#include <ibtfun.h>
C
      CALL QENTER('GETGC2SCA')
C
C  Check that scalled diamagnetic spin-orbit integrals exist
C
      LUPROP = -1
      CALL GPOPEN(LUPROP,'AOPROPER','OLD','SEQUENTIAL','UNFORMATTED',
     &   IDUMMY,OLDDX)
      REWIND LUPROP
      IF (.NOT.FNDLB2('D-SSO XX',RTNLBL,LUPROP)) THEN
         CALL PR1INT('ELGDSCA',WRK,LWRK,IDUMMY,
     &               IDUMMY,.FALSE.,.FALSE.,0)
      END IF
      CALL GPCLOSE(LUPROP,'KEEP')
C
C Get triplet density
C
      KFREE=1
      LFREE=LWRK
      CALL MEMGET('REAL',KDEN1,N2ASHX,WRK,KFREE,LFREE)   
      IF (HSROHF) THEN  
         CALL DUNIT(WRK(KDEN1),NASHT)
      ELSE 
          CALL MEMGET('REAL',KCREF,NCREF,WRK,KFREE,LFREE)
          CALL GETREF(WRK(KCREF),NCREF)
          CALL RSPDM(IREFSY,IREFSY,NCREF,NCREF,WRK(KCREF),WRK(KCREF),
     &           WRK(KDEN1),DUMMY,
     &           1,0,.FALSE.,.TRUE.,XINDX,WRK,KFREE,LFREE)
      END IF
C     CALL RSPDM(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR, RHO1,RHO2,
C    *                 ISPIN1,ISPIN2,TDM,NORHO2,XINDX,WORK,
C    *                 KFREE,LFREE)
C
C Get diamagnetic spin-orbit integrals
C
      CALL MEMGET('REAL',KPRPMO,N2ORBX,WRK,KFREE,LFREE)
C
C Expectation value
C
      TRPLET = .TRUE.
      KSYMOP=1
      IPRG=IPRRSP
      CALL PRPGET(
     &   'D-SSO XX',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &   WRK(KFREE),LFREE,IPRG
     &   )
      CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GC2(1,1),10,'D-SSO XX')
      CALL PRPGET(
     &   'D-SSO YY',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &   WRK(KFREE),LFREE,IPRG
     &   )
      CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GC2(2,2),10,'D-SSO YY')
      CALL PRPGET(
     &   'D-SSO ZZ',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &   WRK(KFREE),LFREE,IPRG
     &   )
      CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GC2(3,3),10,'D-SSO ZZ')
      IF (IBTXOR(ISYMAX(1,1),ISYMAX(2,1)).EQ.0) THEN
         CALL PRPGET(
     &      'D-SSO XY',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &      WRK(KFREE),LFREE,IPRG
     &      )
         CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GC2(1,2),10,'D-SSO XY')
         CALL PRPGET(
     &      'D-SSO YX',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &      WRK(KFREE),LFREE,IPRG
     &      )
         CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GC2(2,1),10,'D-SSO YX')
      END IF
      IF (IBTXOR(ISYMAX(1,1),ISYMAX(3,1)).EQ.0) THEN
         CALL PRPGET(
     &      'D-SSO XZ',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &      WRK(KFREE),LFREE,IPRG
     &      )
         CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GC2(1,3),10,'D-SSO XZ')
         CALL PRPGET(
     &      'D-SSO ZX',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &      WRK(KFREE),LFREE,IPRG
     &      )
         CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GC2(3,1),10,'D-SSO ZX')
      END IF
      IF (IBTXOR(ISYMAX(2,1),ISYMAX(3,1)).EQ.0) THEN
         CALL PRPGET(
     &      'D-SSO YZ',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &      WRK(KFREE),LFREE,IPRG
     &      )
         CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GC2(2,3),10,'D-SSO YZ')
         CALL PRPGET(
     &      'D-SSO ZY',CMO,WRK(KPRPMO),KSYMOP,ANTSYM,
     &      WRK(KFREE),LFREE,IPRG
     &      )
         CALL PRPONE(WRK(KPRPMO),WRK(KDEN1),GC2(3,2),10,'D-SSO ZY')
      END IF
      TRPLET = .FALSE.
C
C Scale
C
      CALL DSCAL(9,ALPHA2/(2*REFSPIN),GC2,1)
C
C Substract one-electron gauge correction
C 
      DO I=1,3
        DO J=1,3
          GC2(I,J)=GC2(I,J)-GC1(I,J)
        END DO
      END DO
C
C Finish
C
      CALL MEMREL('GETGC2SCA',WRK,1,1,KFREE,LFREE)
      CALL QEXIT('GETGC2SCA')
      RETURN
      END
C /* Deck gset  */
      SUBROUTINE GOZHSO(I,J,XPROP)
#include <implicit.h>
#include "rspprp.h"
#include <inflr.h>
#include <wrkrsp.h>
#include <gtensor.h>
      CHARACTER*8 ALAB,BLAB
      ALAB=LBLLR(KSYMOP,I)
      BLAB=LBLLR(KSYMOP,J)
      IF ((ALAB(2:7).EQ.'ANGMOM'.OR.ALAB(2:7).EQ.'ANGECC')
     &   .AND. 
     &   (BLAB(3:8).EQ.'SPNORB'
     &    .OR.BLAB(3:8).EQ.'MNF-SO'
     &    .OR.BLAB(3:8).EQ.'SPNSCA'
     &    ) ) THEN
         IF (ALAB(1:1).EQ.'X') IA = 1
         IF (ALAB(1:1).EQ.'Y') IA = 2
         IF (ALAB(1:1).EQ.'Z') IA = 3
         IF (BLAB(1:1).EQ.'X') IB = 1
         IF (BLAB(1:1).EQ.'Y') IB = 2
         IF (BLAB(1:1).EQ.'Z') IB = 3
         IF (ADDSO) THEN
            GOZSO(IA,IB,1) = XPROP/REFSPIN
         ELSE
            IF (DOOZSO1.AND.BLAB(2:8).EQ.'1SPNORB') 
     &         GOZSO(IA,IB,1) = XPROP/REFSPIN
            IF (DOOZSO2.AND.(BLAB(2:8).EQ.'2SPNORB'
     &               .OR. BLAB(2:8).EQ.'1MNF-SO'
     &               .OR. BLAB(2:8).EQ.'1SPNSCA'))
     &         GOZSO(IA,IB,2) = XPROP/REFSPIN
         END IF
      END IF
      RETURN
      END
C  /* Deck gprint  */
      SUBROUTINE GPRINT(WORK,LWORK)
#include <implicit.h>
      DIMENSION WORK(*)
#include <infpri.h>
#include <gtensor.h>
#include <codata.h>
#include <gfac.h>
#include <orgcom.h>
#include <priunit.h>
#include <inforb.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <symmet.h>
      INTEGER SWRITE
      EXTERNAL SWRITE
C
      PARAMETER (D1=1.0D0, PPM=1.0D6, DP5=0.5D0, D0=0.0D0)
      CHARACTER*10 FMT
      CHARACTER*80 GHEAD
      LOGICAL LINEAR, DIAG, PR(5)
      DIMENSION GSHIFT(3,3), GEIG(3), PVAL(6)
C
      IF (.NOT.GCALC) RETURN
      CALL QENTER('GPRINT')
      IF (NINT(REFSPIN).NE.0) THEN
         CALL TITLER(' Electronic g-tensor (ppm)','=',3)
C
         IF (ECC) THEN
            WRITE(LUPRI,'(A,3F10.6/)')
     &      '    Gauge origin (electronic charge centroid)',
     &      (ECCORG(I),I=1,3)
         END IF
C
         IF (DORMC) THEN
            CALL HEADER('Relativistic mass contribution',3)
            WRITE(LUPRI,'(F10.0)') GRELMC*PPM
         END IF
C
         IF (DOGC1) THEN
            CALL HEADER('One-electron gauge correction ',3)
            WRITE(LUPRI,'(3F10.0)')((PPM*GC1(I,J),J=1,3),I=1,3)
         END IF
C
         IF (DOGC2) THEN
            CALL HEADER('Two-electron gauge correction ',3)
            WRITE(LUPRI,'(3F10.0)')((PPM*GC2(I,J),J=1,3),I=1,3)
         END IF
C
         IF (ADDSO) THEN
            CALL HEADER('Spin-orbit+orbital-Zeeman contribution',3)
            WRITE(LUPRI,'(3F10.0)')((PPM*GOZSO(I,J,1),J=1,3),I=1,3)
         ELSE
C
            IF (DOOZSO1) THEN
               CALL HEADER('One-electron spin-orbit' //
     &         '+orbital-Zeeman contribution',3)
               WRITE(LUPRI,'(3F10.0)')((PPM*GOZSO(I,J,1),J=1,3),I=1,3)
            END IF
C
            IF (DOOZSO2) THEN
               IF (MNFSO.OR.SCALED_CHARGES) THEN
                  DO I=1,3
                     DO J=1,3
                        GOZSO(I,J,2)=GOZSO(I,J,2)-GOZSO(I,J,1)
                     END DO
                  END DO
               END IF
               IF (SPLITSO2) THEN
                  CALL HEADER('Two-electron spin-own-orbit' //
     &            '+orbital-Zeeman contribution',3)
                  WRITE(LUPRI,'(3F10.0)')
     &               ((PPM*GOZSO(I,J,2),J=1,3),I=1,3)
                  CALL HEADER('Two-electron spin-other-orbit' //
     &            '+orbital-Zeeman contribution',3)
                  WRITE(LUPRI,'(3F10.0)')
     &               ((PPM*GOZSO(I,J,3),J=1,3),I=1,3)
               ELSE
                  CALL HEADER('Two-electron spin-orbit' //
     &            '+orbital-Zeeman contribution',3)
                  WRITE(LUPRI,'(3F10.0)')
     &               ((PPM*GOZSO(I,J,2),J=1,3),I=1,3)
               END IF
            END IF
         END IF
C
C        IF (DOALL) THEN
         IF (.TRUE.) THEN
            CALL HEADER('Total g-tensor shift',3)
            CALL DZERO(GTENSOR,9)
            DO I=1,3
               GTENSOR(I,I)=GTENSOR(I,I)+GRELMC
            END DO
            DO I=1,3
               DO J=1,3
                  GTENSOR(I,J)=GTENSOR(I,J)
     &            +GC1(I,J)+GC2(I,J)
                  DO ISO=1,NSO
                     GTENSOR(I,J)=GTENSOR(I,J)+GOZSO(I,J,ISO)
                  END DO
               END DO
            END DO
            WRITE(LUPRI,'(3F10.0)')((PPM*GTENSOR(I,J),J=1,3),I=1,3)
C
            CALL HEADER('Total g-tensor',3)
            DO I=1,3
               GTENSOR(I,I)=GTENSOR(I,I)+GFAC
            END DO
            WRITE(LUPRI,'(3F10.6)')((GTENSOR(I,J),J=1,3),I=1,3)
            GTRACE=GTENSOR(1,1)+GTENSOR(2,2)+GTENSOR(3,3) 
C
C All shifts in one table
C
            LFMT=7
            FMT='(FX.0)'
            WRITE(FMT(3:3),'(I1)')LFMT
            CALL HEADER('G-shift components (ppm)',0)
            DIAG=NSYM.GE.4
            LINEAR=DIAG.AND.ABS(GTENSOR(1,1)-GTENSOR(2,2)).LT.1D-8
            GHEAD='@G'
            LOFF=14
            PR(1)=LINEAR
            PR(2)=DIAG
            PR(3)=(ISYMAX(1,2).EQ.ISYMAX(2,2)) 
            PR(4)=(ISYMAX(1,2).EQ.ISYMAX(3,2)) 
            PR(5)=(ISYMAX(2,2).EQ.ISYMAX(3,2)) 
            IF (LINEAR) THEN
               LOFF=SWRITE(GHEAD,LOFF,LFMT,'_|_')
               LOFF=SWRITE(GHEAD,LOFF,LFMT,'||')
            ELSE IF (DIAG) THEN
               LOFF=SWRITE(GHEAD,LOFF,LFMT,'xx')
               LOFF=SWRITE(GHEAD,LOFF,LFMT,'yy')
               LOFF=SWRITE(GHEAD,LOFF,LFMT,'zz')
            ELSE
               LOFF=SWRITE(GHEAD,LOFF,LFMT,'xx')
               LOFF=SWRITE(GHEAD,LOFF,LFMT,'yy')
               LOFF=SWRITE(GHEAD,LOFF,LFMT,'zz')
               IF (PR(3)) THEN
                  LOFF=SWRITE(GHEAD,LOFF,LFMT,'xy')
                  LOFF=SWRITE(GHEAD,LOFF,LFMT,'yx')
               END IF
               IF (PR(4)) THEN
                  LOFF=SWRITE(GHEAD,LOFF,LFMT,'xz')
                  LOFF=SWRITE(GHEAD,LOFF,LFMT,'zx')
               END IF
               IF (PR(5)) THEN
                  LOFF=SWRITE(GHEAD,LOFF,LFMT,'yz')
                  LOFF=SWRITE(GHEAD,LOFF,LFMT,'zy')
               END IF
            END IF
            WRITE(LUPRI,'(A)') GHEAD
            CALL FLUSH(LUPRI)
            GRMC=GRELMC*PPM
            CALL DZERO(GSHIFT,9)
            DO I=1,3
               GSHIFT(I,I)=GSHIFT(I,I)+GRMC
            END DO
            IF (DORMC)CALL TENPRI(LUPRI,FMT,'@G RMC',PR,GSHIFT)
            CALL DSCAL(9,PPM,GC1,1)
            CALL DSCAL(9,PPM,GC2,1)
            CALL DSCAL(9*NSO,PPM,GOZSO,1)
            DO I=1,3
               DO J=1,3
                  GSHIFT(I,J) = GSHIFT(I,J) + GC1(I,J)+GC2(I,J)
     &               + GOZSO(I,J,1) + GOZSO(I,J,2) + GOZSO(I,J,3)
               END DO
            END DO
            IF (DOGC1) 
     &         CALL TENPRI(LUPRI,FMT,'@G GC1    ',PR,GC1)
            IF (DOGC2) 
     &         CALL TENPRI(LUPRI,FMT,'@G GC2    ',PR,GC2)
            IF (ADDSO) THEN
               CALL TENPRI(
     &            LUPRI,FMT,'@G OZ-SO  ',PR,GOZSO(1,1,1)
     &            )
            ELSE
               IF (DOOZSO1) THEN
                  CALL TENPRI(
     &               LUPRI,FMT,'@G OZ-SO1 ',PR,GOZSO(1,1,1)
     &               )
               END IF
               IF (DOOZSO2) THEN
                  IF (SPLITSO2) THEN
                     CALL TENPRI(
     &                  LUPRI,FMT,'@G OZ-SOWN',PR,GOZSO(1,1,2)
     &                  )
                     CALL TENPRI(
     &                  LUPRI,FMT,'@G OZ-SOTH',PR,GOZSO(1,1,3)
     &                  )
                  ELSE
                     CALL TENPRI(
     &                  LUPRI,FMT,'@G OZ-SO2 ',PR,GOZSO(1,1,2)
     &                  )
                  END IF
               END IF
            END IF
            CALL TENPRI(LUPRI,FMT,'@G Total  ',PR,GSHIFT)
C
C Get principal axes by diagonalizing g
C Since the tensor may in princple be non-symmetric, we first symmetrize the
C tensor. We diagonalize g(T)*g which is the relevant quantity for the
C energy levels in the ESR spectrum (equivalently we find the right
C eigenvalues/vectors
C
C           GTENSOR(2,1)=(GTENSOR(1,2)+GTENSOR(2,1))*DP5
C           GTENSOR(3,1)=(GTENSOR(1,3)+GTENSOR(3,1))*DP5
C           GTENSOR(3,2)=(GTENSOR(2,3)+GTENSOR(3,2))*DP5
            CALL DCOPY(9,GTENSOR,1,WORK,1)
            CALL DGEMM('T','N',3,3,3,
     &         D1,WORK,3,
     &            WORK,3,
     &         D0,GTENSOR,3
     &         )
            CALL DUNIT(GSHIFT,3)
            IJ = 1
            DO I = 1, 3
               DO J = 1, I
                  PVAL(IJ) = GTENSOR(I,J)
                  IJ = IJ + 1
               END DO
            END DO
            CALL JACO(PVAL,GSHIFT,3,3,3,WORK,WORK(10))
            GEIG(1)=SQRT(PVAL(1))
            GEIG(2)=SQRT(PVAL(3))
            GEIG(3)=SQRT(PVAL(6))
            CALL HEADER('G-shift/tensor eigenvalues and cosines',0)
            WRITE(LUPRI,'(8X,2A14,A24)')'g-shift','g-tensor',
     &         'direction cosines'
            DO K=1,3
               WRITE(LUPRI,'(3X,I5,2F14.6,3F8.4)') 
     &          K,GEIG(K)-GFAC,GEIG(K),
     &         (GSHIFT(J,K),J=1,3) 
            END DO
            GEIGSUM=GEIG(1)+GEIG(2)+GEIG(3)
            IF (ABS(GTRACE-GEIGSUM).GT.1E-6) THEN
               WRITE(LUERR,*)'GPRINT:GTRACE=',GTRACE
               WRITE(LUERR,*)'GPRINT:GAVE  =',GEIGSUM
               CALL QUIT('GPRINT:TRACE NOT EQUAL TO AVERAGE')
            END IF
            WRITE(LUPRI,'(A8,2F14.6)') 'ave' , 
     &         GEIGSUM/3-GFAC,GEIGSUM/3
         END IF
      END IF
      CALL QEXIT('GPRINT')
      RETURN
      END
C /* Deck bdg */
      BLOCK DATA BDG
#include <implicit.h>
#include <gtensor.h>
      DATA GCALC,G_ALL /.FALSE.,.FALSE./
      DATA G_LINE/' '/
      DATA DOALL,DORMC,DOGC1,DOGC2,DOOZSO,DOOZSO1,DOOZSO2,SCALED_CHARGES
     &   /.TRUE.,7*.FALSE./
      END
C  /* Deck avetwo */
      SUBROUTINE AVEDSO(CMO,XINDX,WORK,LWORK)
C
C CALCULATE AVERAGE VALUE OF PROPERTIES
C
C
#include <implicit.h>
      DIMENSION CMO(*), XINDX(*), WORK(LWORK)
#include <dummy.h>
C
C
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infrsp.h>
#include <infind.h>
#include <wrkrsp.h>
#include "rspprp.h"
#include "infave.h"
#include <infpri.h>
#include <inftap.h>
#include <iratdef.h>
#include <aovec.h>
#include <gtensor.h>
#include <codata.h>
#include <priunit.h>
#include <blocks.h>
#include <infinp.h>
C
      PARAMETER ( D1 = 1.0D0 , D2 = 2.0D0 , D8 = 8.0D0 )
      PARAMETER ( D1INF = 0.99999D0 , CKMXPR = 1.0D-6 )

      DIMENSION ISYMDM(2), IFCTYP(2)
      CHARACTER SPD(7)
      DATA SPD/'S','P','D','F','G','H','I'/
      LOGICAL ANTI, PANTI, NODPTR, DIA2SO, NODV, NOPV, NOCONT, TTIME
      LOGICAL FNDLAB, NORHO2, TDM, OLDDX, RETUR, NOBLK, DEBUG
      CHARACTER*5 STRING

      CALL QENTER('AVEDSO')
C
      ANTI = .FALSE.
      NODPTR = .FALSE.
      DIA2SO = .TRUE.
      NOBLK = .FALSE.
      STRING="     "
      CALL GETENV("PRINT",STRING)
      IF (STRING.NE."     " ) THEN
         READ(STRING,'(I5)') IPRINT
      ELSE
         IPRINT = 0
      END IF
      DEBUG = IPRINT.GT.10
C
C     Generate one-electron density matrices
C
C     Inactive
C
      
C     Active
      ISPIN1 = 1
      ISPIN2 = 0
      ILSYM  = IREFSY
      IRSYM  = IREFSY
      NCLDIM = NCREF
      NCRDIM = NCREF
      NORHO2 = .TRUE.
      TDM    = .FALSE.
C
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET('REAL',KUDV,N2ASHX,WORK,KFREE,LFREE)
C     
      NOPV = NASHT.EQ.1 .OR. HSROHF
      IF (NOPV) THEN
         CALL DUNIT(WORK(KUDV),NASHT)
      ELSE     
         CALL MEMGET('REAL',KPV, 0 ,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KCL,NCREF,WORK,KFREE,LFREE)
         CALL GETREF(WORK(KCL),NCREF)
         CALL RSPDM(ILSYM,IRSYM,NCLDIM,NCRDIM,WORK(KCL),WORK(KCL),
     &           WORK(KUDV),WORK(KPV),
     &           ISPIN1,ISPIN2,TDM,NORHO2,XINDX,WORK,
     &           KFREE,LFREE)
         CALL MEMREL('AVEDSO:RSPDM',WORK,KCL,KCL,KFREE,LFREE)
      END IF  
      IF (DEBUG) THEN
         CALL HEADER('AVEDSO|One electron spin density',0)
         CALL OUTPUT(WORK(KUDV),1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
      END IF
C
C     Allocate for active and inactive density matrix (AO)
C
      NDMAT = 2 
      CALL MEMGET('REAL',KDMAO,NDMAT*N2BASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMSO,NDMAT*N2BASX,WORK,KFREE,LFREE)
C
C     Transform to SO basis 
C
      CALL GTDMSO(WORK(KUDV),CMO,WORK(KDMSO),WORK(KDMSO+N2BASX),
     &   WORK(KFREE))
      KLAST = KFREE
      IF  (DEBUG) THEN
         CALL HEADER('Inactive density SO basis',0)
         CALL OUTPUT(WORK(KDMSO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('Active density SO basis',0)
         CALL OUTPUT(WORK(KDMSO+N2BASX),1,NBAST,1,NBAST,NBAST,NBAST,
     &      1,LUPRI)
      END IF
      CALL DSOTAO(WORK(KDMSO),WORK(KDMAO),NBAST,0,IPRINT)
      CALL DSOTAO(WORK(KDMSO+N2BASX),WORK(KDMAO+N2BASX),NBAST,0,IPRINT)
      CALL MEMREL('AVEDSO<-DSOTAO',WORK,KDMAO,KDMSO,KFREE,LFREE)
      IF  (DEBUG) THEN
         CALL HEADER('Inactive density AO basis',0)
         CALL OUTPUT(WORK(KDMAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('Active density AO basis',0)
         CALL OUTPUT(WORK(KDMAO+N2BASX),1,NBAST,1,NBAST,NBAST,NBAST,
     &      1,LUPRI)
      END IF
C
C Write interface file for Lushingtons 2ge2c
C
      IF (NASHT.EQ.1) THEN
         CALL LUSHIF(CMO,WORK(KDMAO),WORK(KDMAO+N2BASX),
     &      WORK,KFREE,LFREE)
      END IF


C
C     *******************************************************
C     ***** Set up COMMON /BLOCKS/ for PSORT and TWOINT *****
C     *******************************************************
C
C     KJSTRS = KLAST
C     KNPRIM = KJSTRS + (MXSHEL*MXAOVC*2 + 1)/IRAT
C     KNCONT = KNPRIM + (MXSHEL*MXAOVC*2 + 1)/IRAT
C     KIORBS = KNCONT + (MXSHEL*MXAOVC*2 + 1)/IRAT
C     KJORBS = KIORBS + (MXSHEL*MXAOVC + 1)/IRAT
C     KKORBS = KJORBS + (MXSHEL*MXAOVC + 1)/IRAT
C     KLAST  = KKORBS + (MXSHEL*MXAOVC + 1)/IRAT
C     IF (KLAST .GT. LWORK) CALL STOPIT('AVEDSO','PAOVEC',KLAST,LWORK)
C     LWORK   = LWORK - KLAST + 1
C
      CALL MEMGET('INTE',KJSTRS,MXSHEL*MXAOVC*2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNPRIM,MXSHEL*MXAOVC*2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCONT,MXSHEL*MXAOVC*2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIORBS,MXSHEL*MXAOVC  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJORBS,MXSHEL*MXAOVC  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KKORBS,MXSHEL*MXAOVC  ,WORK,KFREE,LFREE)
C
      IPRPAO=MAX(0,IPRINT)
      CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            WORK(KJORBS),WORK(KKORBS),0,NOBLK,IPRPAO)
C
C     KLAST = KJORBS
C     LWORK = LWORK - KLAST + 1
C
      CALL MEMREL('AVEDIA:PAOVEC',WORK,KJORBS,KJORBS,KFREE,LFREE)
      
      IF (.NOT.NOPV) THEN
C
C     Transform two-electron density matrix to SO basis
C
         PANTI = .FALSE.
         LUPAO=-1
         CALL GPOPEN(LUPAO,'GC2PAO','NEW',' ',' ',IDUMMY,.FALSE.)
         CALL PTRAN(NODPTR,WORK(KFREE),LFREE,IPRINT,ANTI,PANTI,DIA2SO,
     &      ZFS2EL)
         CALL PSORG(WORK(KFREE),WORK(KFREE),LFREE,WORK(KNCONT),IPRINT,
     &   ANTI,PANTI)
         PANTI = .TRUE.
         LUPAS=-1
         CALL GPOPEN(LUPAS,'GC2PAS','NEW',' ',' ',IDUMMY,.FALSE.)
         CALL PTRAN(NODPTR,WORK(KFREE),LFREE,IPRINT,ANTI,PANTI,DIA2SO,
     &      ZFS2EL)
         CALL PSORG(WORK(KFREE),WORK(KFREE),LFREE,WORK(KNCONT),IPRINT,
     &   ANTI,PANTI)
C
C Now we have two files with symmetrized and anti-symmetrized densities
C with respect to particle permutation. These must be combined to match
C square loop in TWOINT
C
         NNSHL=MAXSHL*(MAXSHL+1)/2
         NNSHL2=NNSHL*NNSHL
         CALL MEMGET('INTE',ISHIND,NNSHL2,WORK,KFREE,LFREE)
         LUPS=-1
         LUPA=-1
         CALL GPOPEN(LUPS,'GC2PS','NEW',' ',' ',IDUMMY,.FALSE.)
         CALL GPOPEN(LUPA,'GC2PA','NEW',' ',' ',IDUMMY,.FALSE.)
         CALL P2ORG(WORK(KNCONT),WORK(ISHIND),NNSHL,LUPAO,LUPAS,LUPS,
     &      WORK,KFREE,LFREE,.FALSE.)
         CALL P2ORG(WORK(KNCONT),WORK(ISHIND),NNSHL,LUPAO,LUPAS,LUPA,
     &      WORK,KFREE,LFREE,.TRUE.)
         CALL GPCLOSE(LUPAO,'DELETE')
         CALL GPCLOSE(LUPAS,'DELETE')
         CALL GPCLOSE(LUPS,'KEEP')
         CALL GPCLOSE(LUPA,'KEEP')
         CALL GPOPEN(LUPAO,'GC2PS','OLD',' ','UNFORMATTED',IDUMMY,
     &      .FALSE.)
         CALL GPOPEN(LUPAS,'GC2PA','OLD',' ','UNFORMATTED',IDUMMY,
     &      .FALSE.)
      END IF
C
C     Call HERMIT to evaluate expectation value. A lot of these variables 
C     may be of interest to control through a input routine.
C
      ISYMDM(1) = 1
      ISYMDM(2) = 1
      IFCTYP(1) = 13
      IFCTYP(2) = 13
      ITYPE     = 10
      MAXDIF    = 2
      JATOM     = 0
      NODV      = .FALSE.
      NOCONT    = .FALSE.
      TTIME     = .FALSE.
      RETUR     = .FALSE.
      IPRNTA    = 0
      IPRNTB    = 0
      IPRNTC    = 0
      IPRNTD    = 0
      I2TYP     = 0
      CALL TWOINT(WORK(KFREE),LFREE,WORK(KFREE),WORK(KFREE),WORK(KDMAO),
     &            NDMAT,IREPDM,IFCTYP,DUMMY,IDUMMY,IDUMMY,1,ITYPE,
     &            MAXDIF,JATOM,NODV,NOPV,NOCONT,TTIME,IPRINT,IPRNTA,
     &            IPRNTB,IPRNTC,IPRNTD,RETUR,IDUMMY,I2TYP,WORK(KJSTRS),
     &            WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,
     &            DUMMY,.FALSE.)
      FAC=-ALPHA2/(4*REFSPIN)
      CALL DSCAL(9,FAC,GC2,1)
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Two-electron diamagnetic spin-orbit part to '//
     &               'g tensor',-1)
         CALL OUTPUT(GC2,1,3,1,3,3,3,1,LUPRI)
      END IF
C
      CALL MEMREL('AVEDSO',WORK,1,1,KFREE,LFREE)
      CALL QEXIT('AVEDSO')
      RETURN
      END
      SUBROUTINE LUSHIF(CMO,DENSNG,DENTRP,WORK,KFREE,LFREE)
C
C Write an interface file for Gerry Lushingtons 2esg2 program
C for testing ROHF evaluation of 2-electron gauge correction
C
#include <implicit.h>
      DIMENSION CMO(*), DENSNG(*), DENTRP(*), WORK(*)
#include <inforb.h> /* NBAST, NISHT, NORBT, N2ORBX */
#include <mxcent.h> /* MXCENT */
#include <nuclei.h> /* NUCIND */
#include <maxorb.h>
#include <infinp.h> /* ISPIN  */
#include <blocks.h> /* MAXSHL, NHKTSH, NUCOSH, NRCSH, NCNTSH */
#include <molinp.h> /* MLINE  */
#include <aovec.h>
#include <primit.h> /* PRIEXP, PRICCF */
#include <inftap.h> /* LUSIFC, FNSIFC, LBSIFC */
#include <infpri.h> /* LUERR */
#include <iratdef.h>/* IRAT */
#include <maxash.h>
#include <infind.h> /* ISMO */
#include <pgroup.h> /* REP */
      PARAMETER (D0=0.0D0, DH=0.5D0, D1=1.0D0, D2=2.0D0)
      CHARACTER SPD(8)
      CHARACTER*10 FMT
      CHARACTER*80 LINE
      DATA SPD/'S','P','D','F','G','H','I','J'/
      DIMENSION IVEC(MXCORB)
      DIMENSION IBF(MXCENT)
      LOGICAL DONE(MXSHEL)
      DO I=1,MAXSHL
         DONE(I)=.FALSE.
      END DO
C
C Logical units used by 2ESG2
C
      LUSHB=81
      LUSHPA=87
      LUSHPB=88
      LUSHP =89
      LUSHPS=90
C
C Basis set info
C
      CALL NEWOPN(LUSHB,'LUSH11')
      REWIND LUSHB
C Number of basis functions
      WRITE(LUSHB,'(A)') '$NBF  '
      WRITE(LUSHB,'(I4)') NBAST
C Number of electrons
      WRITE(LUSHB,'(A)') '$NEL  '
      WRITE(LUSHB,'(I4)') 2*NISHT+1
C Number of centers
      WRITE(LUSHB,'(A)') '$NOC  '
      WRITE(LUSHB,'(I4)') NUCIND ! or NUCDEP?
C Flag ?
      WRITE(LUSHB,'(A)') '$ICAL '
      WRITE(LUSHB,'(I4)') 0
C Multiplicity
      WRITE(LUSHB,'(A)') '$MUL  '
      WRITE(LUSHB,'(I4)') ISPIN
C Number of molecular oribtals
      WRITE(LUSHB,'(A)') '$NMO  '
      WRITE(LUSHB,'(I4)') NORBT
C IBF ? - something with centers - shells per center
      DO ICEN=1,NCNTSH(MAXSHL)
         IBF(ICEN)=0
      END DO
      DO ISHL=1,MAXSHL
         IBF(NCNTSH(ISHL)) = IBF(NCNTSH(ISHL)) + 1
      END DO
      WRITE(LUSHB,'(A)') '$IBF  '
      WRITE(LUSHB,'(2I4)') (IBF(I),I=1,NUCIND)
C
C Basis info  - append LUSH.BAS written by MOLOUT
C
      OPEN(71,FILE='LUSH.BAS',STATUS='OLD')
 10   CONTINUE
         READ(71,'(A80)',END=20) LINE
         WRITE(LUSHB,'(A80)') LINE
         GOTO 10
 20   CONTINUE
      CLOSE(71)

#ifdef extralush
      /* Not needed
C
C Following from interface file (it should exist in response)
C
C     print *,'Sirius inteface'
C     print *,'LUSIFC',LUSIFC
C     print *,'FNSIFC',FNSIFC
C     print *,'LBSIFC',LBSIFC
      OPEN (LUSIFC, FILE=FNSIFC, FORM='UNFORMATTED')
      REWIND LUSIFC
      CALL MOLLAB(LBSIFC,LUSIFC,LUERR)
      READ(LUSIFC)DUM,DUM,DUM,ENERGY
      READ(LUSIFC)
      READ(LUSIFC)
      READ(LUSIFC)
      READ(LUSIFC)
C
C Extract diagonal of Fock matrix for eigenvalues
C
C  F -> EIG -> EIGSRT (sorted eigenvalues)
C
      CALL MEMGET('REAL',KF,N2ORBT,WORK,KFREE,LFREE)     
      CALL MEMGET('REAL',KEIG,NORBT,WORK,KFREE,LFREE)   
      CALL READSQ(LUSIFC,IRAT*N2ORBT,WORK(KF))
      DO ISYM=1,NSYM
         NORBI=NORB(ISYM)
         IF=I2ORB(ISYM)+KF
         IEIG=IORB(ISYM)+KEIG
         CALL DCOPY(NORBI,WORK(IF),NORBI+1,WORK(IEIG),1)
      END DO
C
C EIG contains 2*e(i) for inactive and e(i) for active
C
C Scale the open shell orbital energy so we can sort them
C
      IOP = ISW(NISHT+1)
      WORK(KEIG+IOP-1)=2*WORK(KEIG+IOP-1)
C
C New order to index vector
C
      CALL INDSRT(WORK(KEIG),IVEC,NORBT)
C
C Eigenvalues
      CALL DSCAL(NORBT,DH,WORK(KEIG),1)
      WRITE(LUSHB,'(A)')'$EIG  '
      DO I=1,NORBT
         IRREP=ISMO(IVEC(I))-1
         WRITE(LUSHB,'(I11,F20.12,A5)') I,WORK(KEIG+I-1),REP(IRREP)
      END DO
      CALL MEMREL('LUSHIF',WORK,KF,KF,KFREE,LFREE)
C Total energy
      WRITE(LUSHB,'(A)')'$ENGT '
      WRITE(LUSHB,'(F17.9)') ENERGY
C Orbitals
      WRITE(LUSHB,'(A)')'$VEC  '
C Lines per orbital
      NLINES=NBAST/5
C Entries in last line
      NREM=MOD(NBAST,5)
      IA=0
      DO I=1,NORBT
         DO NN=1,NLINES
            WRITE(LUSHB,'(I2,I3,5D15.8)')I,NN,(CMO(IA+J),J=1,5)
            IA=IA+5
         END DO
         IF (NREM.GT.0) THEN
            WRITE(LUSHB,'(I2,I3,5D15.8)')I,NN,(CMO(IA+J),J=1,NREM)
            IA=IA+NREM
         END IF
      END DO
C Order?  
      WRITE(LUSHB,'(A)')'$ORD  '
      WRITE(LUSHB,'(20I4)')(IVEC(I),I=1,NORBT)
C Core?  
      WRITE(LUSHB,'(A)')'$CORE '
      DO I=1,NORBT
         WRITE(LUSHB,'(I4,F20.10)')I,0.0D0
      END DO
C ?  
      WRITE(LUSHB,'(A)')'$COLEX'
      DO I=1,NORBT
         DO J=1,I
            WRITE(LUSHB,'(2I4,2F15.8)') I,J,0.,0.
         END DO
      END DO
       
      Not needed */ 
#endif
      WRITE(LUSHB,'(/A)')'$END  '
C
C Write densities
C
C  input 1.inactive singlet DI = Da + Db (doubly occupied)
C        2.active   triplet DA = Da - Db (open shells)
C  write
C
C  1. alpha Pa
C  2. beta  Pb
C  3. total P =DI+DA (for one open shell singlet DA = triplet DA)
C  4. spin  PS=DA
C
C     P = Pa + Pb => Pa = (P + Ps)/2 = DI/2 + DA
C     Ps = Pa - Pb =>  Pb = (P - Ps)/2 = DI/2
C
      CALL NEWOPN(LUSHPA,'LUSH17')
      CALL NEWOPN(LUSHPB,'LUSH18')
      CALL NEWOPN(LUSHP ,'LUSH19')
      CALL NEWOPN(LUSHPS,'LUSH20')
      FMT='(5D15.8)'

      CALL MEMGET('REAL',KP,N2BASX,WORK,KFREE,LFREE)
      CALL DCOPY(N2BASX,DENTRP,1,WORK(KP),1)
      CALL DAXPY(N2BASX,DH,DENSNG,1,WORK(KP),1)
      REWIND LUSHPA
      CALL WRFMTD(LUSHPA,FMT,WORK(KP),NBAST)
      CLOSE(LUSHPA)

      CALL DCOPY(N2BASX,DENSNG,1,WORK(KP),1)
      CALL DSCAL(N2BASX,DH,WORK(KP),1)
      REWIND LUSHPB
      CALL WRFMTD(LUSHPB,FMT,WORK(KP),NBAST)
      CLOSE(LUSHPB)

      REWIND LUSHP
      CALL DCOPY(N2BASX,DENSNG,1,WORK(KP),1)
      CALL DAXPY(N2BASX,D1,DENTRP,1,WORK(KP),1)
      CALL WRFMTD(LUSHP,FMT,WORK(KP),NBAST)
      CLOSE(LUSHP)

      REWIND LUSHPS
      CALL WRFMTD(LUSHPS,FMT,DENTRP,NBAST)
      CLOSE(LUSHPS)
      RETURN
      END

      SUBROUTINE WRFMTI(LU,FMT,IVEC,N)
      INTEGER LU, N, IVEC(N)
      CHARACTER*(*) FMT
      WRITE(LU,FMT)(IVEC(I),I=1,N)
      RETURN
      END

      SUBROUTINE WRFMTD(LU,FMT,RMAT,N)
#include <implicit.h>
      INTEGER LU, N
      DIMENSION RMAT(N,N)
      CHARACTER*(*) FMT
      DO I=1,N
         WRITE(LU,FMT)(RMAT(I,J),J=1,N)
      END DO
      RETURN
      END

      SUBROUTINE INDSRT(EIG,IVEC,N)
#include <implicit.h>
C
C Sort vector EIG(N)
C IVEC(N) will give the index in the original vector
C for the sorted values
C
      DIMENSION EIG(N), IVEC(N)
      DO I=1,N
         IVEC(I)=I
      END DO
      DO I=1,N-1
         EIGMIN=EIG(I)
         JMIN=I
         DO J=I+1,N
            IF (EIG(J).LT.EIGMIN) THEN
               JMIN=J
               EIGMIN=EIG(J)
            END IF
         END DO
         TMP=EIG(I)
         EIG(I)=EIG(JMIN)
         EIG(JMIN)=TMP
         ITMP=IVEC(I)
         IVEC(I)=IVEC(JMIN)
         IVEC(JMIN)=ITMP
      END DO
      RETURN
      END

      SUBROUTINE NEWOPN(LU,NAME)
      INTEGER LU
      CHARACTER*(*) NAME
      OPEN (LU, FILE=NAME, STATUS='UNKNOWN')
      CLOSE(LU, STATUS='DELETE')
      OPEN (LU, FILE=NAME, STATUS='NEW')
      END
      SUBROUTINE P2ORG(NCONTS,ISHIND,NNSHL,LUPAO,LUPAS,LUP,
     &   WORK,KFREE,LFREE,PANTI)
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <infind.h>
#include <priunit.h>
#include <aovec.h>
#include <blocks.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <symmet.h>
      PARAMETER ( D1 = 1.0D0 , D2 = 2.0D0 )
      DIMENSION NCONTS(MXSHEL,MXAOVC,2)
      DIMENSION ISHIND(NNSHL,NNSHL)
      DIMENSION WORK(*)
      LOGICAL PANTI
      LOGICAL DEBUG 
      DEBUG=.FALSE.
C
C Establish relationship between shell-block quartet and record number NREC
C on file genereated by PSORg
C Save maximum record length in MXNP
C
      CALL QENTER('P2ORG ')
      REWIND LUPAO
      MXNP=0
      NREC=0
      DO ISHLA=1,MAXSHL
         DO ISHLB=1,ISHLA
            IAB=IROW(ISHLA)+ISHLB
            DO ISHLC=1,ISHLA
               IF (ISHLC.EQ.ISHLA) THEN
                  MXSHLD=ISHLB
               ELSE
                  MXSHLD=ISHLC
               END IF
               DO ISHLD=1,MXSHLD
                  ICD=IROW(ISHLC)+ISHLD
                  NREC=NREC+1
                  ISHIND(IAB,ICD)=NREC
                  ISHIND(ICD,IAB)=NREC
                  READ(LUPAO) NPMAT
                  READ(LUPAO)
                  MXNP=MAX(NPMAT,MXNP)
                  IF (DEBUG) THEN
                     WRITE(LUPRI,'(A,8I5)')
     &               'ISHL[A..D],IAB,ICD,NREC,NPMAT',
     &               ISHLA,ISHLB,ISHLC,ISHLD,
     &               IAB,ICD,NREC,NPMAT
                  END IF
               END DO
            END DO
         END DO
      END DO   
      IF (DEBUG) THEN
         WRITE(LUPRI,*)'  ISHIND'
         DO I=1,NNSHL
            WRITE(LUPRI,*)(ISHIND(I,J),J=1,NNSHL)
         END DO
      END IF
C
C Square loop
C
      CALL MEMGET('REAL',KPAO,MXNP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPAS,MXNP,WORK,KFREE,LFREE)
      REWIND LUP
      DO ISHLA=1,MAXSHL
         KHKTA=KHKTSH(ISHLA)
         MULA=MULT(ISTBSH(ISHLA))
         NSETA=NSETSH(ISHLA,1)
         NRCA = 0
         DO IA=1,NSETA
            NRCAI = NCONTS(ISHLA,IA,1)
            IF (NRCAI.LT.0) NRCAI = 1
            NRCA = NRCA + NRCAI
         END DO
         DO ISHLB=1,ISHLA
            KHKTB=KHKTSH(ISHLB)
            MULB=MULT(ISTBSH(ISHLB))
            NSETB=NSETSH(ISHLB,1)
            NRCB = 0
            DO IB=1,NSETB
               NRCBI = NCONTS(ISHLB,IB,1)
               IF (NRCBI.LT.0) NRCBI = 1
               NRCB = NRCB + NRCBI
            END DO
            IAB=IROW(ISHLA)+ISHLB
            DO ISHLC=1,MAXSHL
               KHKTC=KHKTSH(ISHLC)
               MULC=MULT(ISTBSH(ISHLC))
               NSETC=NSETSH(ISHLC,1)
               NRCC = 0
               DO IC=1,NSETC
                  NRCCI = NCONTS(ISHLC,IC,1)
                  IF (NRCCI.LT.0) NRCCI = 1
                  NRCC = NRCC + NRCCI
               END DO
               DO ISHLD=1,ISHLC
                  KHKTD=KHKTSH(ISHLD)
                  MULD=MULT(ISTBSH(ISHLD))
                  NSETD=NSETSH(ISHLD,1)
                  NRCD = 0
                  DO ID=1,NSETD
                     NRCDI = NCONTS(ISHLD,ID,1)
                     IF (NRCDI.LT.0) NRCDI = 1
                     NRCD = NRCD + NRCDI
                  END DO
                  ICD=IROW(ISHLC)+ISHLD
                  NREC=ISHIND(IAB,ICD)
                  REWIND LUPAO
                  REWIND LUPAS
                  DO IREC=1,NREC-1
                     READ(LUPAO)
                     READ(LUPAS)
                     READ(LUPAO)
                     READ(LUPAS)
                  END DO
                  READ(LUPAO) NPMAT1
                  READ(LUPAS) NPMAT2
                  IF (DEBUG) THEN
                     WRITE(LUPRI,'(A,9I5)')
     &               'ISHL[A..D],IAB,ICD,NREC,NPMAT1/2',
     &               ISHLA,ISHLB,ISHLC,ISHLD,
     &               IAB,ICD,NREC,NPMAT1,NPMAT2
                  END IF
                  IF (NPMAT1.NE.NPMAT2) THEN
                     CALL QUIT('AVEDSO: error in density files')
                  END IF
                  CALL READT(LUPAO,NPMAT1,WORK(KPAO))
                  CALL READT(LUPAS,NPMAT2,WORK(KPAS))
                  NAB=KHKTA*MULA*NRCA*KHKTB*MULB*NRCB
                  NCD=KHKTC*MULC*NRCC*KHKTD*MULD*NRCD
                  IF (DEBUG) THEN
                     WRITE (LUPRI,'(A,2I5)')'NAB,NCD',NAB,NCD
                  END IF
                  IF (IAB.GE.ICD) THEN
                     IF (DEBUG) THEN
                        WRITE (LUPRI,'(A)')'Pao'
                        CALL OUTPUT(WORK(KPAO),1,NAB,1,NCD,NAB,NCD,
     &                     1,LUPRI)
                        WRITE (LUPRI,'(A)')'Pas'
                        CALL OUTPUT(WORK(KPAS),1,NAB,1,NCD,NAB,NCD,
     &                     1,LUPRI)
                     END IF
                     CALL DAXPY(NPMAT1,D1,WORK(KPAS),1,WORK(KPAO),1)
                     IF (DEBUG) THEN
                        WRITE (LUPRI,'(A)')'Pao+Pas'
                        CALL OUTPUT(WORK(KPAO),1,NAB,1,NCD,NAB,NCD,
     &                     1,LUPRI)
                     END IF
                  ELSE
                     IF (DEBUG) THEN
                        WRITE (LUPRI,'(A)')'Pao'
                        CALL OUTPUT(WORK(KPAO),1,NCD,1,NAB,NCD,NAB,
     &                     1,LUPRI)
                        WRITE (LUPRI,'(A)')'Pas'
                        CALL OUTPUT(WORK(KPAS),1,NCD,1,NAB,NCD,NAB,
     &                     1,LUPRI)
                     END IF
                     CALL DAXPY(NPMAT1,-D1,WORK(KPAO),1,WORK(KPAS),1)
                     CALL DSCAL(NPMAT1,-D1,WORK(KPAS),1)
                     IF (DEBUG) THEN
                        WRITE (LUPRI,'(A)')'Pao-Pas'
                        CALL OUTPUT(WORK(KPAS),1,NCD,1,NAB,NCD,NAB,
     &                     1,LUPRI)
                     END IF
                     CALL MTRSP(NCD,NAB,WORK(KPAS),NCD,WORK(KPAO),NAB)
                     IF (DEBUG) THEN
                        WRITE (LUPRI,'(A)')'(Pao-Pas)(T)'
                        CALL OUTPUT(WORK(KPAO),1,NAB,1,NCD,NAB,NCD,
     &                     1,LUPRI)
                     END IF
                  END IF
C
C Scale density with 2, yet unclear why
C
                  CALL DSCAL(NPMAT1,D2,WORK(KPAO),1)
                  WRITE (LUP) NPMAT1
                  CALL WRITT(LUP,NPMAT1,WORK(KPAO))
                  IF (DEBUG) THEN
                     WRITE (LUPRI,'(A)')'Psum'
                     CALL OUTPUT(WORK(KPAO),1,NAB,1,NCD,NAB,NCD,1,LUPRI)
                  END IF
               END DO
            END DO
         END DO
      END DO
      CALL MEMREL('AVEDSO:NPLOOP',WORK,KPAO,KPAO,KFREE,LFREE)
      CALL QEXIT('P2ORG ')
      RETURN
      END

      SUBROUTINE TENPRI(LU,FMT,LABEL,PR,G)
      IMPLICIT NONE
      INTEGER LU
      CHARACTER*(*) FMT,LABEL
      LOGICAL LINEAR,DIAGONAL,PR(5)
      REAL*8 G(3,3)
      CHARACTER*80 LINE
      INTEGER OFFSET, LFMT, SWRITE, RWRITE
      EXTERNAL SWRITE, RWRITE

      LFMT=7
      LINE=LABEL
      IF (PR(1)) THEN
         OFFSET=10
         OFFSET=RWRITE(LINE,OFFSET,LFMT,FMT,G(1,1))
         OFFSET=RWRITE(LINE,OFFSET,LFMT,FMT,G(3,3))
      ELSEIF (PR(2)) THEN
         OFFSET=10
         OFFSET=RWRITE(LINE,OFFSET,LFMT,FMT,G(1,1))
         OFFSET=RWRITE(LINE,OFFSET,LFMT,FMT,G(2,2))
         OFFSET=RWRITE(LINE,OFFSET,LFMT,FMT,G(3,3))
      ELSE
         OFFSET=10
         OFFSET=RWRITE(LINE,OFFSET,LFMT,FMT,G(1,1))
         OFFSET=RWRITE(LINE,OFFSET,LFMT,FMT,G(2,2))
         OFFSET=RWRITE(LINE,OFFSET,LFMT,FMT,G(3,3))
         IF (PR(3)) THEN
            OFFSET=RWRITE(LINE,OFFSET,LFMT,FMT,G(1,2))
            OFFSET=RWRITE(LINE,OFFSET,LFMT,FMT,G(2,1))
         END IF
         IF (PR(4)) THEN
            OFFSET=RWRITE(LINE,OFFSET,LFMT,FMT,G(1,3))
            OFFSET=RWRITE(LINE,OFFSET,LFMT,FMT,G(3,1))
         END IF
         IF (PR(5)) THEN
            OFFSET=RWRITE(LINE,OFFSET,LFMT,FMT,G(2,3))
            OFFSET=RWRITE(LINE,OFFSET,LFMT,FMT,G(3,2))
         END IF
      END IF
      WRITE(LU,'(A)')LINE
      RETURN
      END

      INTEGER FUNCTION SWRITE(STRING,OFFSET,LENGTH,INSTRING)
      IMPLICIT NONE
      CHARACTER*(*) STRING, INSTRING
      INTEGER OFFSET, LENGTH
      STRING(OFFSET:OFFSET+LENGTH)=INSTRING 
      SWRITE=OFFSET+LENGTH
      RETURN
      END

      INTEGER FUNCTION RWRITE(STRING,OFFSET,LENGTH,FMT,INREAL)
      IMPLICIT NONE
      CHARACTER*(*) STRING, FMT
      REAL*8 INREAL
      INTEGER OFFSET, LENGTH
      CHARACTER*80 BUF
      BUF=' '
      WRITE(BUF,FMT)INREAL
      STRING(OFFSET:OFFSET+LENGTH)=BUF(1:LENGTH)
      RWRITE=OFFSET+LENGTH
      RETURN
      END
