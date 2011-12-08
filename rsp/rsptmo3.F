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
      SUBROUTINE CREXSR
C
C SORT THE EXCITATION ENERGIES  AFTER SYMMETRY
C THE NUMBER OF EXCITATION ENERGIES OF A GIVEN SYMMETRY ARE
C STORED IN NEXCR WITH OFFSET IN IEXCR
C
#include <implicit.h>
C
#include <priunit.h>
#include <infrsp.h>
#include <inforb.h>
#include <rspprp.h>
#include <indqr.h>
#include <indcr.h>
#include <infpri.h>
C
      IF (IPRRSP.GT.40) THEN
         WRITE(LUPRI,'(/A,I5,2X,A)')
     *   ' LIST OF NEXLBL= ',NEXLBL,'EXCITATIONS FOR CUBIC RESPONSE'
         WRITE(LUPRI,'(A)')' BEFORE SORTING'
         DO 90 I = 1,NEXLBL
            WRITE(LUPRI,'(/A)')'     I   ISEXCR(I)  JEXCR(I)'
            WRITE(LUPRI,'(/3I10)')I,ISEXCR(I),JEXCR(I)
 90      CONTINUE
      END IF
      ISMMAX = 1
      DO 100 ISYM = 1,NSYM
         I = ISMMAX
         DO 200 ISRT = I,NEXLBL
            IF ( ISEXCR(ISRT).EQ.ISYM) THEN
                ICR   = ISEXCR(ISRT)
                JCR   = JEXCR(ISRT)
                ISEXCR(ISRT) = ISEXCR(ISMMAX)
                JEXCR(ISRT)  = JEXCR(ISMMAX)
                ISEXCR(ISMMAX) = ICR
                JEXCR(ISMMAX)  = JCR
                ISMMAX = ISMMAX + 1
            END IF
 200     CONTINUE
 100  CONTINUE
      ITOTOP = 0
      DO 300 ISYM = 1,NSYM
         INUMOP = 0
         DO 400 I = 1,NEXLBL
            IF(ISEXCR(I).EQ.ISYM) THEN
               INUMOP = INUMOP + 1
            END IF
 400     CONTINUE
         IEXCR(ISYM) = ITOTOP
         ITOTOP = ITOTOP + INUMOP
         NEXCR(ISYM) = INUMOP
 300  CONTINUE
      IF (IPRRSP.GT.15) THEN
         WRITE(LUPRI,'(/A,I5,2X,A)')
     *   ' LIST OF NEXLBL= ',NEXLBL,'EXCITATIONS FOR QUADRATIC RESPONSE'
         WRITE(LUPRI,'(A)')' AFTER SORTING'
         DO 310 I = 1,NEXLBL
            WRITE(LUPRI,'(/A)')'     I   ISEXCR(I)  JEXCR(I)'
            WRITE(LUPRI,'(/3I10)')I,ISEXCR(I),JEXCR(I)
 310     CONTINUE
      END IF
      IF (IPRRSP.GT.10) THEN
         WRITE(LUPRI,'(/A)')
     *   ' NUMBER OF EXCITATIONS IN VARIOUS SYMMETRIES'
         WRITE(LUPRI,'(/A)')' NEXCR(ISYM),ISYM=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NEXCR(ISYM),ISYM=1,NSYM )
      END IF
      RETURN
      END
      SUBROUTINE CRPPVE(CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,
     *                  WRK,LWRK)
C
C  Purpose:
C     CALCULATION OF EXCITATION ENERGIES AND EIGENVECTORS
C     FOR CUBIC RESPONSE
C
#include <implicit.h>
#include <dummy.h>
#include <iratdef.h>
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
      CHARACTER*8 BLANK
      PARAMETER (D0 = 0.0D0)
      PARAMETER   (BLANK='        ')
      PARAMETER ( MAXSIM = 15 )
C
C Used from common blocks:
C  /INFRSP/ : most items (/INFRSP/ gives control information for
C                         the response calculation(s) )
C  /WRKRSP/ :
C
C  /INFORB/ : MULD2H
#include <priunit.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <infpp.h>
#include <maxorb.h>
#include <infvar.h>
#include <qrinf.h>
#include <indcr.h>
#include <indqr.h>
#include <infpri.h>
#include <inforb.h>
#include <inftap.h>
C
C     space allocation for reduced E(2) and reduced S(2)
      KREDE  = 1
      KREDS  = KREDE  + MAXRM*MAXRM
      KIBTYP = KREDS  + MAXRM*MAXRM
      KEIVAL = KIBTYP + MAXRM
      KRESID = KEIVAL + MAXRM
      KEIVEC = KRESID + MAXRM
      KWRK1  = KEIVEC + MAXRM*MAXRM
      LWRK1  = LWRK + 1 - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('CRPPVE 1',KWRK1-1,LWRK)
      IF (IPRRSP .GT. 2) THEN
         WRITE(LUPRI,*)' IN CRPPVE: MAXRM      ',MAXRM
         WRITE(LUPRI,*)' IN CRPPVE: LWRK ,LWRK1',LWRK,LWRK1
         WRITE(LUPRI,*)' IN CRPPVE: THCPP      ',THCPP
      END IF
C
      IF (LWRK1 .LT. 3*KZYVAR) THEN
         WRITE (LUERR,9000) LWRK1,3*KZYVAR
         CALL QTRACE(LUERR)
         CALL QUIT('ERROR, INSUFFICIENT SPACE FOR CRPPVE')
      ENDIF
 9000    FORMAT(/' CRPPVE, work space too small for 3 (z,y)-vectors',
     *          /'         had:',I10,', need more than:',I10)
C
      KZRED  = 0
      KZYRED = 0
      THCRSP = THCPP
      MAXIT  = MAXITP
C
C     Call RSPCTL to solve propagator eigen problem
C
      CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *            .FALSE.,BLANK,BLANK,DUMMY,DUMMY,WRK(KREDE),WRK(KREDS),
     *            WRK(KIBTYP),WRK(KEIVAL),WRK(KRESID),WRK(KEIVEC),
     *            XINDX,WRK(KWRK1),LWRK1)
C     CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
C    *            LINEQ,GD,REDGD,REDE,REDS,
C    *            IBTYP,EIVAL,EIVEC,XINDX,WRK,LWRK)
C
C CALCULATE EIGENVECTORS
C
C MAXIMUM NUMBER OF TRIAL VECTORS
C
      NSIMMA = MIN(MAXSIM, INT((LWRK1-KZVAR)/KZYVAR))
      IF (NSIMMA.GT.KEXCNV) THEN
         NSIM = KEXCNV
      ELSE
         NSIM = NSIMMA
      END IF
      IF ( NSIM.EQ.0) THEN
         WRITE(LUPRI,'(/A,I5)')
     *   ' CRPPVE: TOO LITTLE WORK SPACE, NSIM= ',NSIM
         CALL QUIT(' CRPPVE: TOO LITTLE WORK SPACE')
      END IF
      KBVECS = KWRK1
      KWRK2  = KBVECS + NSIM*KZYVAR
      LWRK2  = LWRK   - KWRK2
      IF (LWRK2.LT.0) CALL ERRWRK('CRPPVE 2',KWRK2-1,LWRK)
C
C SAVE EXCITATION ENERGIES
C
      DO 400 I = 1,KEXCNV
         EXCIT2(KSYMOP,I) = WRK(KEIVAL-1+I)
 400  CONTINUE
C
      DO 500 ISIM = 1,KEXCNV,NSIM
         NBX = MIN( NSIM,(KEXCNV+1-ISIM) )
         IBOFF = ISIM - 1
         CALL RSPEVE(WRK(KIBTYP),WRK(KEIVAL),WRK(KEIVEC),WRK(KBVECS),
     *               WRK(KWRK2),NBX,IBOFF)
C        CALL RSPEVE(IBTYP,EIVAL,EIVEC,BVECS,WRK,NBX,IBOFF)
         DO 750 IVEC = 1,NBX
            IBV  = (IVEC-1)*KZYVAR + KBVECS
            NVEC = INCRLR('EXCITLAB',EXCIT2(KSYMOP,ISIM-1+IVEC),KSYMOP)
            IF (IPRRSP.GT.10) THEN
               WRITE(LUPRI,'(/A)')' SINGLET EXCITATION '
               WRITE(LUPRI,'(/A,I5,3X,A,I5,/A,2X,D12.6)')
     *         ' CRPPVE: EIGENVECTOR ',ISIM-1+IVEC,
     *         ' SYMMETRY',KSYMOP,' EXCITATION ENERGY',
     *         EXCIT2(KSYMOP,ISIM-1+IVEC)
            END IF
            IF (IPRRSP.GT.200)
     *         CALL OUTPUT(WRK(IBV),1,KZVAR,1,2,KZVAR,2,1,LUPRI)
            CALL RSPPRO(WRK(IBV+KZCONF),KZVAR,LUPRI)
            CALL RSPANC(WRK(IBV),KZCONF,KZVAR,
     *                  MULD2H(KSYMOP,IREFSY),XINDX,MULD2H,LUPRI)
            CALL WRTRSP(LURSP,KZYVAR,WRK(IBV),'EXCITLAB',BLANK,
     *                  EXCIT2(KSYMOP,ISIM-1+IVEC),D0,KSYMOP,
     &                  ISIM-1+IVEC,WRK(KRESID - 1 + IVEC),D0)
 750     CONTINUE
 500  CONTINUE
C
C *** END OF CRPPVE
C
      RETURN
      END
