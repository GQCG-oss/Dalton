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
Revision 1.2  2001/01/17 16:09:33  vebjornb
Added dimension checks before calls to DGEMM in response module

Revision 1.2  2000/05/24 18:52:10  hjj
reuse Y2X for KX2Y
change 2 to * in DIMENSION
new GETREF calls
===========================================================================
#endif
      SUBROUTINE CRLRV1(CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
C
#include <implicit.h>
#include <iratdef.h>
#include <thrzer.h>
C
C     Purpose:
C     Solve linear equations for one-index response vectors
C
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
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
#include <inftap.h>
#include <qrinf.h>
#include <inforb.h>
#include <infcr.h>
#include <indcr.h>
#include <infdim.h>
C
C Local variables
C
      LOGICAL FOUND, CONV
      CHARACTER*8 LRLAB(MXLRCR),BLANK
C
      PARAMETER ( BLANK='        ', D0=0.0D0 )
C
      CALL QENTER('CRLRV1')
      IF (NLRCR(KSYMOP) .LE. 0) GO TO 9999
C
C     Count number of different property labels
C     and max number of frequencies, maximum number of frequencies
C     is used for memory allocation.
C
      NLRLAB = 0
      MXFRQ  = 0
      IF (IPRLR .GT. 20) THEN
         WRITE (LUPRI,*)
     *   'Test output from CRLRV1 of all labels and frequencies'
      END IF
C
      DO 220 IOP = 1,NLRCR(KSYMOP)
         IOPVEC = ILRCR(KSYMOP) + IOP
         IF (IPRLR .GT. 20) THEN
            WRITE (LUPRI,*) IOP, IOPVEC, CRLBL(IOPVEC),CRFREQ(IOPVEC)
         END IF
         DO 200 ILRLAB = 1,NLRLAB
            IF (CRLBL(IOPVEC) .EQ. LRLAB(ILRLAB)) GO TO 220
  200    CONTINUE
         NLRLAB = NLRLAB + 1
         LRLAB(NLRLAB) = CRLBL(IOPVEC)
C
         NFRQ = 1
         DO 210 JOP = IOP+1,NLRCR(KSYMOP)
            JOPVEC = ILRCR(KSYMOP) + JOP
            IF (CRLBL(JOPVEC) .EQ. LRLAB(ILRLAB)) NFRQ = NFRQ + 1
  210    CONTINUE
         IF (IPRLR .GT. 20) WRITE (LUPRI,*) 'New label, NFREQ =',NFRQ
         MXFRQ = MAX(MXFRQ,NFRQ)
  220 CONTINUE
C
C     Allocate memory
C     3*KZYVAR is an estimate of space needed in RSPCTL
C
      LNEEDA = 2*MAXRM*MAXRM + MAXRM + 4*KZYVAR
      LNEEDB = 1 + 2*MAXRM
      NFRQMX = (LWRK - LNEEDA) / LNEEDB
      IF (NFRQMX .LT. MXFRQ) THEN
         WRITE (LUERR,9100) LWRK,LNEEDA+LNEEDB*MXFRQ
         CALL QTRACE(LUERR)
         CALL QUIT('CRLRV1: INSUFFICIENT SPACE TO SOLVE '//
     &             'LINEAR EQUATIONS')
      ENDIF
 9100 FORMAT(/' CRLRV1, work space too small for 3 (z,y)-vectors',
     *       /'         had:',I10,', need more than:',I10)
C
      KREDE  = 1
      KREDS  = KREDE  + MAXRM*MAXRM
      KIBTYP = KREDS  + MAXRM*MAXRM
      KEIVAL = KIBTYP + MAXRM
      KRESID = KEIVAL + MXFRQ
      KEIVEC = KRESID + MXFRQ
      KREDGD = KEIVEC + MAXRM*MXFRQ
      KGD    = KREDGD + MAXRM*MXFRQ
      KBVECS  = KGD    + KZYVAR
      NFRQMX = (LWRK - KBVECS - KZYVAR) / KZYVAR
      NFRQMX = MIN(NFRQMX,MXFRQ)
      KWRK1  = KBVECS + NFRQMX*KZYVAR
      LWRK1  = LWRK   - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('CRLRV1 1',KWRK1-1,LWRK)
C
      THCRSP = THCLR
      MAXIT  = MAXITL
C
C     Call RSPCTL to solve linear set of response equations
C
      DO 680 ILRLAB = 1,NLRLAB
         CALL GETGPV(LRLAB(ILRLAB),FC,FV,CMO,UDV,PV,XINDX,ANTSYM,
     &               WRK(KGD),LWRK1)
         GDNORM = DNRM2(KZYVAR,WRK(KGD),1)
         NFREQ  = 0
         CALL DZERO(WRK(KEIVAL),NLRCR(KSYMOP))
         DO 410 IOP = 1,NLRCR(KSYMOP)
            IOPVEC = ILRCR(KSYMOP) + IOP
            AFRQ = ABS(CRFREQ(IOPVEC))
            IF (CRLBL(IOPVEC) .EQ. LRLAB(ILRLAB)) THEN
               CALL REARSP(LURSP,KLEN,WRK(KBVECS),CRLBL(IOPVEC),
     &                     BLANK,AFRQ,D0,KSYMOP,0,THCLR,
     &                     FOUND,CONV,ATNSYM)
               IF (FOUND .AND. CONV) THEN
                  WRITE(LUPRI,'(/A,/A12,A10,/A42)')
     *                 ' Converged solution vector already on file',
     *                 CRLBL(IOPVEC),'freq',
     *                 ' -----------------------------------------'
                  WRITE(LUPRI,'(F22.6)') AFRQ
               ELSE
                  DO IJ = 0, NFREQ
                     IF (ABS(WRK(KEIVAL+IJ)-AFRQ) .LT. THRZER .AND.
     &                    .NOT. AFRQ .LT. THRZER) GOTO 410
                  END DO
                  WRK(KEIVAL+NFREQ) = AFRQ
                  NFREQ = NFREQ + 1
               END IF
            END IF
  410    CONTINUE
         IF (NFREQ.EQ.0) GO TO 680
C
         WRITE (LUPRI,'(//A,I3,/2A,/A,(T25,5F10.6))')
     &   ' CRLRV1 -- linear response calculation for symmetry',KSYMOP,
     &   ' CRLRV1 -- operator label : ',LRLAB(ILRLAB),
     &   ' CRLRV1 -- frequencies :',(WRK(KEIVAL+I),I=0,NFREQ-1)
         IF (IPRRSP.GT.10) THEN
            WRITE (LUPRI,*) 'Column 1 = Z, Column 2 = Y'
            CALL OUTPUT(WRK(KGD),1,KZVAR,1,2,KZVAR,2,1,LUPRI)
         END IF
         IF (GDNORM .LT. THRNRM) THEN
            WRITE (LUPRI,*) ' >>> RSPCTL skipped because norm of'
            WRITE (LUPRI,*) '     property vector is only',GDNORM
            CALL DZERO(WRK(KBVECS),KZYVAR)
            DO 450 IFREQ = 1,NFREQ
               CALL WRTRSP(LURSP,KZYVAR,WRK(KBVECS),LRLAB(ILRLAB),BLANK,
     *                     WRK(KEIVAL-1+IFREQ),D0,KSYMOP,0,D0,ANTSYM)
  450       CONTINUE
            GO TO 680
         END IF
C
         KZRED  = 0
         KZYRED = 0
         KEXSIM = NFREQ
         KEXCNV = KEXSIM
         CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,.TRUE.,LRLAB(ILRLAB),
     *               BLANK,WRK(KGD),WRK(KREDGD),WRK(KREDE),WRK(KREDS),
     *               WRK(KIBTYP),WRK(KEIVAL),WRK(KRESID),WRK(KEIVEC),
     *               XINDX,WRK(KWRK1),LWRK1)
C
         DO 580 IFREQ = 1,NFREQ,NFRQMX
            NBX   = MIN(NFRQMX,(NFREQ+1-IFREQ))
            IBOFF = IFREQ - 1
            CALL RSPEVE(WRK(KIBTYP),WRK(KEIVAL),WRK(KEIVEC),
     *                  WRK(KBVECS),WRK(KWRK1),NBX,IBOFF)
C
            JBVECS = KBVECS
            WRITE(LUPRI,'(/A12,A10,A20,/A42)')
     *      LRLAB(ILRLAB),'freq','LR value',
     *      ' -----------------------------------------'
            DO 560 JFREQ = IFREQ,IFREQ-1+NBX
               JEIVEC = KEIVEC + (JFREQ-1)*KZYRED
               VAL = DDOT(KZYRED,WRK(KREDGD),1,WRK(JEIVEC),1)
               WRITE(LUPRI,'(F22.6,F20.8)') WRK(KEIVAL-1+JFREQ),VAL
               IF (IPRRSP.GT.10) THEN
                  WRITE (LUPRI,*) 'Column 1 = Z, Column 2 = Y'
                  CALL OUTPUT(WRK(JBVECS),1,KZVAR,1,2,KZVAR,2,1,LUPRI)
               END IF
               CALL WRTRSP(LURSP,KZYVAR,WRK(JBVECS),LRLAB(ILRLAB),BLANK,
     *                     WRK(KEIVAL-1+JFREQ),D0,KSYMOP,0,
     *                     WRK(KRESID - 1 + JFREQ),ANTSYM)
               JBVECS = JBVECS + KZYVAR
  560       CONTINUE
  580    CONTINUE
  680 CONTINUE
C
C     End of CRLRV1 --
C
 9999 CALL QEXIT('CRLRV1')
      RETURN
      END
      SUBROUTINE CRVEC2(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *                 XINDX,MJWOP,WRK,LWRK)
C
C PURPOSE:
C  CREATE THE SOLUTION VECTORS THAT ARE REQUIRED TO DO
C  CUBIC RESPONSE
C  SOLUTION OF LINEAR EQUATIONS ARE DETERMINED IN CRLRV1.
C
#include <implicit.h>
#include <iratdef.h>
C
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
#include <maxorb.h>
#include <inforb.h>
#include <infvar.h>
#include <infdim.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <indqr.h>
#include <qrinf.h>
      DIMENSION MJWOP(2,MAXWOP,8)
#include <inftap.h>      
#include <infsmo.h>
#include <infhyp.h>
#include <infpp.h>
#include <infpri.h>
#include <infspi.h>
#include <infcr.h>
#include <indcr.h>
#include <inftmo.h>
#include <inftpa.h>
C
      CALL QENTER('CRVEC2')
      IF (TOMOM) THEN
         NLRLB2 = 0
         CALL TM2IND
      ELSE IF (TPAMP) THEN
         NLRLB2 = 0
         CALL TP2IND
      END IF
C
      CALL GPOPEN(LUXYVE,'CRXYVE','UNKNOWN','DIRECT',' ',IRAT*2*NVARMA,
     &            OLDDX)
C
      CALL XYVEC(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *           XINDX,MJWOP,WRK,LWRK)
C
      CALL CRLRV2(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *            XINDX,MJWOP,WRK,LWRK)
C
      CALL QEXIT('CRVEC2')
      RETURN
      END
      SUBROUTINE XYVEC(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *                 XINDX,MJWOP,WRK,LWRK)
C
#include <implicit.h>
#include <iratdef.h>
C
      LOGICAL FOUND, CONV
      CHARACTER*8 BLANK
      PARAMETER   (BLANK='        ')
      PARAMETER (DM1 = -1.0D0, D0 = 0.0D0)
C
C PURPOSE:
C CALCULATION OF THE XY-VECTORS THAT GIVE THE RIGHT-HAND SIDE OF THE
C EIGENVALUE PROBLEM FOR THE LAST N^XY-RESPONSE VECTORS
C
C
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
      CHARACTER*8 XOPLBL,YOPLBL
C
C
#include <priunit.h>
#include <infrsp.h>
#include <maxorb.h>
#include <infvar.h>
#include <qrinf.h>
      DIMENSION MJWOP(2,MAXWOP,8)
#include <rspprp.h>
#include <indqr.h>
#include <infhyp.h>
#include <inforb.h>
#include <infpri.h>
#include <infspi.h>
#include <inftap.h>
#include <wrkrsp.h>
#include <indcr.h>
#include <infcr.h>
#include <infdim.h>
#include <inftpa.h>
#include <inflr.h>
C
C
      CALL QENTER('XYVEC')
      DO 200 IOP = 1,NLRLB2
         XOPLBL = CRLB2(IOP,1)
         YOPLBL = CRLB2(IOP,2)
         FREQX  = CRFRQ2(IOP,1)
         FREQY  = CRFRQ2(IOP,2)
         ISYMX  = ISMCR2(IOP,1)
         ISYMY  = ISMCR2(IOP,2)
         ISYMXY = ISMCR2(IOP,3)
         IXYFNU = INCR2(XOPLBL,YOPLBL,FREQX,FREQY,ISYMX,ISYMY)
         ISPINX = 0
         ISPINY = 0
         ISPNXY = 0
C
         KVECX  = MZYVAR(ISYMX)
         KVECY  = MZYVAR(ISYMY)
         KVECXY = MZYVAR(ISYMXY)
         KZCONF = MZCONF(ISYMXY)
         KZVAR  = MZVAR(ISYMXY)
C
C     Allocate memory for vectors, and check if sufficient memory
C
         KX    = 1
         KY    = KX    + KVECX
         KX3BC = KY    + KVECY
         KY2X  = KX3BC + KVECXY
         KX2Y  = KY2X
         KFREE = KY2X  + KVECXY
         IF (ISYMXY.EQ.1) THEN
            KCREF = KY2X
            KFREE = MAX(KFREE,KCREF + KZCONF)
         ELSE
            KCREF = -999 999 999
         END IF
         LFREE = LWRK - KFREE
         IF (LFREE.LT.0) CALL ERRWRK('XYVEC Y2X',KFREE,LWRK)
C
         CALL REARSP(LURSP,KLEN,WRK(KX),XOPLBL,YOPLBL,FREQX,FREQY,
     &               ISYMX,ISYMY,THCLR,FOUND,CONV,ANTSYM)
         IF (FOUND .AND. CONV) THEN
            WRITE(LUPRI,'(/A,/2A12,2F12.6)')
     *           ' Skip XYVEC since converged solution vector'//
     *           ' already on file',XOPLBL,YOPLBL,FREQX,FREQY
            GO TO 200
         END IF
C
C     Read in Nx and Ny
C
         CALL REARSP(LURSP,KVECX,WRK(KX),XOPLBL,BLANK,FREQX,D0,
     &               ISYMX,0,THCLR,FOUND,CONV,ANTSYM)
         IF (.NOT. (FOUND .AND. CONV)) THEN
            IF (.NOT. FOUND) THEN
               WRITE (LUPRI,'(/3A,F7.5,A,I3,/A)') ' XYVEC: Response '//
     &              'label ',XOPLBL,' with frequency ',FREQX,
     &              ' and symmetry',ISYMX,' not found on file LURSP'
               CALL QUIT('Response vector not found on file')
            ELSE
               WRITE (LUPRI,'(/3A,F7.5,A,I3,/A)') ' XYVEC: Response '//
     &              'label ',XOPLBL,' with frequency ',FREQX,
     &              ' and symmetry',ISYMX,' not converged on file LURSP'
               CALL QUIT('Non-converged response vector on file')
            END IF
         END IF
         IF (FREQX .LT. D0) THEN
            CALL DSWAP(KVECX/2,WRK(KX),1,WRK(KX+KVECX/2),1)
            IF (ANTSYM .LT. D0) CALL DSCAL(KVECX,ANTSYM,WRK(KX),1)
         END IF
C
         CALL REARSP(LURSP,KVECY,WRK(KY),YOPLBL,BLANK,FREQY,D0,
     &               ISYMY,0,THCLR,FOUND,CONV,ANTSYM)
         IF (.NOT. (FOUND .AND. CONV)) THEN
            IF (.NOT. FOUND) THEN
               WRITE (LUPRI,'(/3A,F7.5,A,I3,/A)') ' XYVEC: Response '//
     &              'label ',YOPLBL,' with frequency ',FREQY,
     &              ' and symmetry',ISYMY,' not found on file LURSP'
               CALL QUIT('Response vector not found on file')
            ELSE
               WRITE (LUPRI,'(/3A,F7.5,A,I3,/A)') ' XYVEC: Response '//
     &              'label ',YOPLBL,' with frequency ',FREQY,
     &              ' and symmetry',ISYMY,' not converged on file LURSP'
               CALL QUIT('Non-converged response vector on file')
            END IF
         END IF
         IF (FREQY .LT. D0) THEN
            CALL DSWAP(KVECY/2,WRK(KY),1,WRK(KY+KVECY/2),1)
            IF (ANTSYM .LT. D0) CALL DSCAL(KVECY,ANTSYM,WRK(KY),1)
         END IF
C
C     Calculate T[3] Nx Ny
C
         CALL T3DRV(1,ISYMXY,ISYMX,ISYMY,WRK(KX),WRK(KY),.FALSE.,
     *               WRK(1),-FREQX,-FREQY,XINDX,UDV,PV,MJWOP,
     *               WRK(KX3BC),LWRK-KX3BC,CMO,FC,FV)
C
C     Calculate Y[2] Nx term and store in WRK(KY2X)

C     Add to X2BC
C
         IF (YOPLBL.NE.'EXCITLAB') THEN
            CALL X2INIT(1,KVECXY,KVECX,ISYMXY,ISPNXY,ISYMX,ISPINX,1,
     *             WRK(KX),WRK(KY2X),XINDX,UDV,PV,YOPLBL,ISYMY,ISPINY,
     *                  CMO,MJWOP,WRK(KFREE),LFREE)
            CALL DAXPY(KVECXY,DM1,WRK(KY2X),1,WRK(KX3BC),1)
         END IF
C
C     Calculate X[2] Ny term and store in WRK(KX2Y)
C     Add to X2BC
C
         IF (XOPLBL.NE.'EXCITLAB') THEN
            CALL X2INIT(1,KVECXY,KVECY,ISYMXY,ISPNXY,ISYMY,ISPINY,1,
     *             WRK(KY),WRK(KX2Y),XINDX,UDV,PV,XOPLBL,ISYMX,ISPINX,
     *                  CMO,MJWOP,WRK(KFREE),LFREE)
            CALL DAXPY(KVECXY,DM1,WRK(KX2Y),1,WRK(KX3BC),1)
         END IF
C
C     Project out reference state if sym(|k>) = IREFSY, XY = <0|op|k>
C     Write vector to file
C
         IF (ISYMXY.EQ.1 .AND. KZCONF.GT.0) THEN
            CALL GETREF(WRK(KCREF),KZCONF)
            T1 = DDOT(KZCONF,WRK(KCREF),1,WRK(KX3BC),1)
            CALL DAXPY(KZCONF,(-T1),WRK(KCREF),1,WRK(KX3BC),1)
            T1 = DDOT(KZCONF,WRK(KCREF),1,WRK(KX3BC+KZVAR),1)
            CALL DAXPY(KZCONF,(-T1),WRK(KCREF),1,WRK(KX3BC+KZVAR),1)
         END IF
C
         CALL WRITDX(LUXYVE,IXYFNU,IRAT*KVECXY,WRK(KX3BC))
C
C     Printout of result
C
         IF (IPRRSP.GT.10) THEN
            WRITE(LUPRI,'(//A,2(/A,A12),6(/A,I8),2(/A,D12.6),//A,/A)')
     *               'Characteristics in XYVEC routine',
     *               'XOPLBL =', CRLB2(IOP,1),
     *               'YOPLBL =', CRLB2(IOP,2),
     *               'ISYMX  =', ISYMX,
     *               'ISYMY  =', ISYMY,
     *               'ISYMXY =', ISYMXY,
     *               'KVECX  =', KVECX,
     *               'KVECY  =', KVECY,
     *               'KVECXY =', KVECXY,
     *               'FREQX  =', CRFRQ2(IOP,1),
     *               'FREQY  =', CRFRQ2(IOP,2),
     * '    Final result in XYVEC        Column 1 = Z, Column 2 = Y', 
     * '    ===================== '
            CALL OUTPUT(WRK(KX3BC),1,KZVAR,1,2,KZVAR,2,1,LUPRI)
         END IF
200   CONTINUE
C
      CALL QEXIT('XYVEC')
      RETURN
      END
