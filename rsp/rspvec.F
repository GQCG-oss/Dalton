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
#ifdef REVLOG
===========================================================================
Revision 1.3  2000/05/24 19:13:28  hjj
Inserted @ in output of final results
Polish

960301-hjaaj: modify MULSP for ISPINA=ISPINB=ISPINC=1
930705-hjaaj: lowered IPRRSP from 5 to 2 in QRVEC for print of lr sym.
910411-hjaaj: RSPPRC and RSPPRO use wrong symmetries in QRSMO,
QRHYP, and QREXM (it uses the last KSYMOP from QRLRVE). TODO
===========================================================================
#endif
C  /* Deck qrcalc */
      SUBROUTINE QRCALC(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *                  XINDX,WRK,LWRK)
C
#include <implicit.h>
#include <iratdef.h>
C
C PURPOSE:
C DRIVER ROUTINE FOR QUADRATIC RESPONSE
C
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
#include <mxcent.h>
#include <maxorb.h>
#include <priunit.h>
#include <infvar.h>
#include <infrsp.h>
#include <qrinf.h>
#include <rspprp.h>
#include <infpp.h>
#include <infhyp.h>
#include <infpri.h>
#include <inforb.h>
#include <infsmo.h>
#include <infspi.h>
#include <inftap.h>
C
      CALL QENTER('QRCALC')
      IF (SOPPA) THEN
         WRITE (LUPRI,*)
     &   'QRCALC ERROR: quadratic response not implemented for SOPPA'
         CALL QUIT(
     &   'QRCALC ERROR: quadratic response not implemented for SOPPA')
      END IF
C
C DETERMINE SOLUTION VECTORS FOR EIGENVALUE AND LINEAR EQUATIONS
C WHICH ARE USED IN QUADRATIC RESPONSE
C
      IF (ISPINA*ISPINB*ISPINC .EQ. 1) THEN
         MULSP(1,1) = 1
         WRITE (LUPRI,'(//A/A/A)')
     &   'QRCALC: MULSP(1,1) reset to triplet instead of singlet',
     &   ' because all three operators are of triplet symmetry.',
     &   ' (MULSP(1,1) tells spin symmetry of the product of two'//
     &   ' triplet operators.)'
      END IF
      CALL FLSHFO(LUPRI)
      KMJWOP = 1
      KWRK1  = KMJWOP + (16*MAXWOP + 1)/IRAT
      LWRK1  = LWRK   - KWRK1
      CALL QRVEC(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *           XINDX,WRK(KMJWOP),WRK(KWRK1),LWRK1)
      CALL FLSHFO(LUPRI)
      KVECB  = KWRK1
      KVECC  = KVECB  + MZYVMX
      KWRK1  = KVECC  + MZYVMX
      LWRK1  = LWRK   - KWRK1
C
      IF (LWRK1.LT.0) CALL ERRWRK('RSPVEC 1',KWRK1-1,LWRK)
C
      IF (HYPCAL) THEN
         CALL LRHYP(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *                 XINDX,WRK(KVECB),WRK(KVECC),WRK(KMJWOP),
     *                 WRK(KWRK1),LWRK1)
         CALL FLSHFO(LUPRI)
         IF (SSCOLL) THEN
            KFC1 = KWRK1
            KFC2 = KFC1 + MXCOOR**2
            KFC  = KFC2 + MXCOOR**2
            KSD1 = KFC  + MXCOOR**2
            KSD2 = KSD1 + MXCOOR**2
            KSD  = KSD2 + MXCOOR**2
            KFS1 = KSD  + MXCOOR**2
            KFS2 = KFS1 + MXCOOR**2
            KFS  = KFS2 + MXCOOR**2
            KLAST= KFS  + MXCOOR**2
            LLEFT = LWRK1 - KLAST
            IF (LLEFT.LT.0) CALL ERRWRK('RSPVEC 2',LLEFT-1,LWRK)
         ELSE
            KFC1 = KWRK1
            KFC2 = KWRK1
            KFC  = KWRK1
            KSD1 = KWRK1
            KSD2 = KWRK1
            KSD  = KWRK1
            KFS1 = KWRK1
            KFS2 = KWRK1
            KFS  = KWRK1
            KLAST= KWRK1
            LLEFT=LWRK1
         END IF
         CALL QRHYP(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *                 XINDX,WRK(KVECB),WRK(KVECC),WRK(KMJWOP),
     *                 WRK(KFC1),WRK(KFC2),WRK(KFC),WRK(KSD1),WRK(KSD2),
     *                 WRK(KSD),WRK(KFS1),WRK(KFS2),WRK(KFS),
     *                 WRK(KLAST),LLEFT)
      END IF
      IF (SOMOM) CALL QRSMO(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *                 XINDX,WRK(KVECB),WRK(KVECC),WRK(KMJWOP),
     *                 WRK(KWRK1),LWRK1)
      IF (EXMOM) CALL QREXM(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *                 XINDX,WRK(KVECB),WRK(KVECC),WRK(KMJWOP),
     *                 WRK(KWRK1),LWRK1)
      CALL FLSHFO(LUPRI)
      CALL QEXIT('QRCALC')
      RETURN
      END
C  /* Deck qrvec */
      SUBROUTINE QRVEC(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *                 XINDX,MJWOP,WRK,LWRK)
C
C PURPOSE:
C  CREATE THE SOLUTION VECTORS THAT ARE REQUIRED TO DO
C  QUADRATIC RESPONSE
C  SOLUTION OF LINEAR EQUATIONS ARE DETERMINED IN QRLRVE.
C  EIGENVECTORS ARE DETERMINED IN QRPPVE
C
#include <implicit.h>
#include <iratdef.h>
C
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
C VECTOR TO STORE NPPCNV(8)
      DIMENSION NKEEP(8)
C
C Used from common blocks:
C   INFRSP : LURSP*
C
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infvar.h>
#include <infdim.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <indqr.h>
      DIMENSION MJWOP(2,MAXWOP,8)
#include <qrinf.h>
#include <infsmo.h>
#include <infhyp.h>
#include <infpp.h>
#include <infpri.h>
#include <infspi.h>
#include <inftap.h>
C
C OPEN FILE FOR SOLUTION VECTORS
C
      CALL QENTER('QRVEC')
      IF (EXMOM .AND. (ISPINC.NE.ISPINB)) THEN
         NEXLB2 = 2*NEXLBL
      ELSE
         NEXLB2 = NEXLBL
      END IF
      NSOLVC = NEXLB2 + NTRLBL + NLRLBL
      NVARMA = NCONMA + NWOPMA
      MZYVMX = 2*NVARMA
      IF (NSOLVC .LE. 0) THEN
         WRITE(LUPRI,'(/A)')
     *     ' rspvec: number of solution vectors incorrect'
         WRITE(LUPRI,'(/A)') 'NEXLBL NEXLB2 NLRLBL NTRLBL:'
         WRITE(LUPRI,'(/A,4I5)') NEXLBL,NEXLB2,NLRLBL,NTRLBL
         CALL QUIT(' rspvec: number of solution vectors incorrect')
      END IF
C
C ORDER THE POINTERS TO THE EIGENVECTORS ACCORDING TO SYMMETRY
C
      CALL QREXSR
C
C SET UP VECTORS DEFINING THE NUMBER OF EIGENVECTORS THAT NEED TO
C BE SOLVED
C
      DO 90 ISYM = 1,NSYM
         NKEEP(ISYM)  = NPPCNV(ISYM)
         NPPCNV(ISYM) = NEXQR(ISYM)
         NPPSIM(ISYM) = NEXQR(ISYM)
         NPPSTV(ISYM) = NEXQR(ISYM)
 90   CONTINUE
C
C
C SOLVE THE EIGENVALUE PROBLEM
C
      IOFFSP = 0
      IF ( SOMOM .OR. ( EXMOM .AND. (ISPINC .EQ. ISPINB) ) ) THEN
         TRPLET = ISPINC.NE.0
      ELSE
         TRPLET = .FALSE.
      END IF
C
 80   CONTINUE
C
C REPEAT FOR BOTH SINGLET AND TRIPLET EXCITATION ENERGIES
C IF CALCULATION OF TRANSITION MOMENT CALCULATION BETWEEN
C EXCITED STATES OF SINGLET AND TRIPLET SYMMETRY
C
      DO 100 KSYMOP = 1,NSYM
      IF ( NEXQR(KSYMOP).EQ.0 ) GO TO 100
C     ... skip this symmetry if no operators
         IF (IPRRSP.GE.2) THEN
            WRITE(LUPRI,'(//A/A,I5/A,L5)')
     *         ' Linear response excitations for quadratic response',
     *         ' - symmetry of excitation operator',KSYMOP,
     *         ' - is the operator a triplet operator ? ',TRPLET
         END IF
C
C        DEFINE VARIABLES THAT DEPEND ON SYMMETRY
C
         CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
C
C SET VARIABLES IN SYMMETRY KSYMOP
C
         CALL SETZY(MJWOP)
         IF ( KZVAR.EQ.0) THEN
            NWARN = NWARN + 1
            WRITE(LUPRI,'(/2A,I3)')' ****WARNING******',
     *         ' QRVEC EXCITATIONS: NO VARIABLES IN SYMMETRY',KSYMOP
            GO TO 100
         END IF
         IF (QRREST) GO TO 100
C
C        ... FIND EXCITATION ENERGIES AND TRANSITION MOMENTS
C
         CALL QRPPVE(CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,IOFFSP,
     *               WRK,LWRK)
         CALL FLSHFO(LUPRI)
 100  CONTINUE
C
C RESTORE INFORMATION ABOUT NUMBER OF EXCITATIONS BETWEEN
C WHICH WE WANT TO CALCULATE TRANSITION MOMENTS
C
      DO 110 ISYM = 1,NSYM
         NPPCNV(ISYM) = NKEEP(ISYM)
 110  CONTINUE
      IF (( EXMOM .AND. (ISPINC.NE.ISPINB)) .AND. (.NOT.TRPLET))
     *                                                      THEN
         IOFFSP = IOFFSP + NEXLBL
         TRPLET = .TRUE.
         GO TO 80
      END IF
C
C SET UP POINTER FOR THE LINEAR RESPONSE EQUATIONS THAT MUST
C BE SOLVED FOR QUADRATIC RESPONSE. THE PREVIOUS SET UP GAVE
C THE RIGHT NUMBER OF VECTORS BUT THE FREQUENCIES WERE NOT
C CORRECT FOR FREQUENCIES CONTAINING EXCITATION ENERGIES.
C THE CORRECT EXCITATION ENERGIES ARE STORED IN EXCITA
C
      NTRLBL = 0
      NLRLBL = 0
      IF (HYPCAL) CALL HYPIND
      IF (SOMOM)  CALL SOMIND
      IF (EXMOM)  CALL EXMIND
C
      IF ( NLRLBL.GT.0 ) THEN
C
C SORT SINGLET LINEAR RESPONSE EQUATIONS ACCORDING TO SYMMETRY
C
         CALL QRLRSR
C
C SOLVE SINGLET LINEAR RESPONSE EQUATIONS
C
         TRPLET = .FALSE.
         KEXSIM = 1
         KEXCNV = 1
         DO 200 KSYMOP = 1,NSYM
         IF ( NLRQR(KSYMOP).EQ.0 ) GO TO 200
C     ... skip this symmetry if no operators
            IF (IPRRSP.GE.2) THEN
               WRITE(LUPRI,'(//A/A,I5)')
     *         ' Linear response calculations for quadratic response',
     *         ' - singlet property operator of symmetry',KSYMOP
            END IF
C
C        DEFINE VARIABLES THAT DEPEND ON SYMMETRY
C
            CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,
     *                  WRK,LWRK)
C           CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
C
C SET VARIABLES IN SYMMETRY KSYMOP
C
            CALL SETZY(MJWOP)
            IF ( KZVAR.EQ.0) THEN
               NWARN = NWARN + 1
               WRITE(LUPRI,'(/2A,I3)')' ****WARNING******',
     *            ' QRVEC SINGLET LR: NO VARIABLES IN SYMMETRY',KSYMOP
               GO TO 200
            END IF
            IF (QRREST) GO TO 200
            CALL QRLRVE(CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,
     *                  WRK,LWRK)
            CALL FLSHFO(LUPRI)
  200    CONTINUE
      END IF
C
      IF ( NTRLBL.GT.0 ) THEN
C
C SORT TRIPLET LINEAR RESPONSE EQUATIONS ACCORDING TO SYMMETRY
C
         CALL QRTRSR
C
C SOLVE LINEAR RESPONSE EQUATIONS
C
         TRPLET = .TRUE.
         KEXSIM = 1
         KEXCNV = 1
         DO 210 KSYMOP = 1,NSYM
         IF ( NTRQR(KSYMOP).EQ.0 ) GO TO 210
C     ... skip this symmetry if no operators
            IF (IPRRSP.GE.2) THEN
               WRITE(LUPRI,'(//A/A,I5)')
     *         ' Linear response calculations for quadratic response',
     *         ' - triplet property operator(s) of symmetry',KSYMOP
            END IF
C
C DEFINE VARIABLES THAT DEPEND ON SYMMETRY
C
            CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
C
C SET VARIABLES IN SYMMETRY KSYMOP
C
            CALL SETZY(MJWOP)
            IF ( KZVAR.EQ.0) THEN
               NWARN = NWARN + 1
               WRITE(LUPRI,'(/2A,I3)')' ****WARNING******',
     *            ' QRVEC TRIPLET LR: NO VARIABLES IN SYMMETRY',KSYMOP
               GO TO 210
            END IF
            IF (QRREST) GO TO 210
               CALL QRTRVE(CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,
     *               WRK,LWRK)
            CALL FLSHFO(LUPRI)
  210    CONTINUE
      END IF
C
C     Delete LR files to save disk space /910109-hjaaj
C
      LURSP1 = -1
      CALL GPOPEN(LURSP1,'RSPRST.E2C','UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      CALL GPCLOSE(LURSP1,'DELETE')
      CALL GPCLOSE(LURSP3,'DELETE')
      CALL GPCLOSE(LURSP5,'DELETE')
C
C     *******************************************
C
      CALL QEXIT('QRVEC')
      RETURN
      END
C  /* Deck qrppve */
      SUBROUTINE QRPPVE(CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,
     *                  IOFFSP,WRK,LWRK)
C
C  Purpose:
C     CALCULATION OF EXCITATION ENERGIES AND EIGENVECTORS
C     FOR QUADRATIC RESPONSE
C
#include <implicit.h>
#include <dummy.h>
#include <iratdef.h>
      CHARACTER*8 BLANK
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
      PARAMETER ( MAXSIM = 15, D0 = 0.0D0, BLANK='        ')
C
C Used from common blocks:
C  /INFRSP/ : most items (/INFRSP/ gives control information for
C                         the response calculation(s) )
C  /WRKRSP/ :
C
C  /INFORB/ : MULD2H
#include <priunit.h>
#include <inftap.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <infpp.h>
#include <maxorb.h>
#include <infvar.h>
#include <qrinf.h>
#include <indqr.h>
#include <infpri.h>
#include <inforb.h>
C
      CALL QENTER('QRPPVE')
C     space allocation for reduced E(2) and reduced S(2)
      KREDE  = 1
      KREDS  = KREDE  + MAXRM*MAXRM
      KIBTYP = KREDS  + MAXRM*MAXRM
      KEIVAL = KIBTYP + MAXRM
      KRESID = KEIVAL + MAXRM
      KEIVEC = KRESID + MAXRM
      KWRK1  = KEIVEC + MAXRM*MAXRM
      LWRK1  = LWRK + 1 - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('QRPPVE 1',KWRK1-1,LWRK)
      IF (IPRRSP .GT. 2) THEN
         WRITE(LUPRI,*)' IN QRPPVE: MAXRM      ',MAXRM
         WRITE(LUPRI,*)' IN QRPPVE: LWRK ,LWRK1',LWRK,LWRK1
         WRITE(LUPRI,*)' IN QRPPVE: THCPP      ',THCPP
      END IF
C
      IF (LWRK1 .LT. 3*KZYVAR) THEN
         WRITE (LUPRI,9000) LWRK1,3*KZYVAR
         CALL QTRACE(LUPRI)
         CALL QUIT('ERROR, INSUFFICIENT WORK MEMORY FOR QRPPVE')
      ENDIF
 9000    FORMAT(/' QRPPVE, work space too small for 3 (z,y)-vectors',
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
     *            .FALSE.,BLANK,BLANK,DUMMY,DUMMY,WRK(KREDE),
     *            WRK(KREDS),WRK(KIBTYP),WRK(KEIVAL),WRK(KRESID),
     *            WRK(KEIVEC),XINDX,WRK(KWRK1),LWRK1)
C     CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
C    *            LINEQ,GP,REDGP,REDE,REDS,
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
     *   ' QRPPVE: TOO LITTLE WORK MEMORY, NSIM= ',NSIM
         CALL QUIT(' QRPPVE: TOO LITTLE WORK MEMORY')
      END IF
      KBVECS = KWRK1
      KWRK2  = KBVECS + NSIM*KZYVAR
      LWRK2  = LWRK   - KWRK2
      IF (LWRK2.LT.0) CALL ERRWRK('QRPPVE 2',KWRK2-1,LWRK)
C
C SAVE EXCITATION ENERGIES
C
      IF ( TRPLET ) THEN
         ISPN = 2
      ELSE
         ISPN = 1
      END IF
      DO 400 I = 1,KEXCNV
         EXCITA(KSYMOP,I,ISPN) = WRK(KEIVAL-1+I)
 400  CONTINUE
C
      DO 500 ISIM = 1,KEXCNV,NSIM
         NBX = MIN( NSIM,(KEXCNV+1-ISIM) )
         IBOFF = ISIM - 1
         CALL RSPEVE(WRK(KIBTYP),WRK(KEIVAL),WRK(KEIVEC),WRK(KBVECS),
     *               WRK(KWRK2),NBX,IBOFF)
C        CALL RSPEVE(IBTYP,EIVAL,EIVEC,BVECS,WRK,NBX,IBOFF)
         DO 750 IVEC = 1,NBX
            IBV    = (IVEC-1)*KZYVAR + KBVECS
            NVEC   = IOFFSP + INQREX(KSYMOP,ISIM-1+IVEC)
            IF (IPRRSP.GT.10) THEN
               IF ( TRPLET ) THEN
                  WRITE(LUPRI,'(/A,I5,3X,A,I5,/A,1P,D14.6)')
     *            '@ QRPPVE: TRIPLET EXCITATION VECTOR',ISIM-1+IVEC,
     *            ' SYMMETRY',KSYMOP,'@ Triplet excitation energy',
     *            EXCITA(KSYMOP,ISIM-1+IVEC,ISPN)
               ELSE
                  WRITE(LUPRI,'(/A,I5,3X,A,I5,/A,1P,D14.6)')
     *            '@ QRPPVE: SINGLET EXCITATION VECTOR',ISIM-1+IVEC,
     *            ' SYMMETRY',KSYMOP,'@ Singlet excitation energy',
     *            EXCITA(KSYMOP,ISIM-1+IVEC,ISPN)
               END IF
            END IF
            IF (IPRRSP.GE.3) THEN
               CALL RSPPRO(WRK(IBV+KZCONF),KZVAR,LUPRI)
               CALL RSPANC(WRK(IBV),KZCONF,KZVAR,
     *                  MULD2H(KSYMOP,IREFSY),XINDX,MULD2H,LUPRI)
               IF (IPRRSP.GT.200)
     *         CALL OUTPUT(WRK(IBV),1,KZVAR,1,2,KZVAR,1,2,LUPRI)
            END IF
            ANTSYM = 1.0D0
            CALL WRTRSP(LURSP,KZYVAR,WRK(IBV),'EXCITLAB',BLANK,
     &                  EXCITA(KSYMOP,ISIM-1+IVEC,ISPN),D0,KSYMOP,
     &                  ISIM-1+IVEC,WRK(KRESID - 1 + IVEC),ANTSYM)
 750     CONTINUE
 500  CONTINUE
C
C *** END OF QRPPVE
C
      CALL QEXIT('QRPPVE')
      RETURN
      END
C  /* Deck qrlrve */
      SUBROUTINE QRLRVE(CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,
     &                  WRK,LWRK)
C
C Revised 28-Feb-1990 HJAaJ
C
#include <implicit.h>
#include <iratdef.h>
#include <dummy.h>
C
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
      LOGICAL FOUND, CONV
      CHARACTER*8 BLANK
      PARAMETER ( THRSML = 1.0D-9, D0 = 0.0D0, BLANK='        ')
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
#include <inforb.h>
#include <inftap.h>
!SONIA: MCDPRJ
#include <infsmo.h>
C
C Local variables
C
      CHARACTER*8 LRLAB(MXLRQR)
C
C DETERMINE SECOND ORDER MOLECULAR PROPERTIES
C
      CALL QENTER('QRLRVE')
      IF (NLRQR(KSYMOP) .LE. 0) GO TO 9999
C
C     count number of different property labels
C     and max number of frequencies
      NLRLAB = 0
      MXFRQ  = 0
      IF (IPRLR .GT. 20) THEN
         WRITE (LUPRI,*)
     *   'Test output from QRLRVE of all labels and frequencies'
      END IF
      DO 220 IOP = 1,NLRQR(KSYMOP)
         IOPVEC = ILRQR(KSYMOP) + IOP
         IF (IPRLR .GT. 20) THEN
            WRITE (LUPRI,*) IOP, IOPVEC, QRLBL(IOPVEC),QRFREQ(IOPVEC)
         END IF
         DO 200 ILRLAB = 1,NLRLAB
            IF (QRLBL(IOPVEC) .EQ. LRLAB(ILRLAB)) GO TO 220
  200    CONTINUE
         NLRLAB = NLRLAB + 1
         LRLAB(NLRLAB) = QRLBL(IOPVEC)
C
         NFRQ = 1
         DO 210 JOP = IOP+1,NLRQR(KSYMOP)
            JOPVEC = ILRQR(KSYMOP) + JOP
            IF (QRLBL(JOPVEC) .EQ. LRLAB(ILRLAB)) NFRQ = NFRQ + 1
  210    CONTINUE
         IF (IPRLR .GT. 20) WRITE (LUPRI,*) 'New label, NFRQ =',NFRQ
         MXFRQ = MAX(MXFRQ,NFRQ)
  220 CONTINUE
C     allocate core
      LNEEDA = 2*MAXRM*MAXRM + MAXRM + 4*KZYVAR
C     ... 3*KZYVAR is an estimate of space needed in RSPCTL
      LNEEDB = 1 + 2*MAXRM
      NFRQMX = (LWRK - LNEEDA) / LNEEDB
      IF (NFRQMX .LT. MXFRQ) THEN
         WRITE (LUPRI,9100) LWRK,LNEEDA+LNEEDB*MXFRQ
         CALL QTRACE(LUPRI)
         CALL QUIT('QRLRVE: INSUFFICIENT WORK MEMORY TO '//
     &             'SOLVE LINEAR EQUATIONS ')
      ENDIF
 9100 FORMAT(/' QRLRVE, work space too small for 3 (z,y)-vectors',
     *       /'         had:',I10,', need more than:',I10)
C
      KREDE  = 1
      KREDS  = KREDE  + MAXRM*MAXRM
      KIBTYP = KREDS  + MAXRM*MAXRM
      KEIVAL = KIBTYP + MAXRM
      KRESID = KEIVAL + MXFRQ
      KEIVEC = KRESID + MXFRQ
      KREDGP = KEIVEC + MAXRM*MXFRQ
      KGP    = KREDGP + MAXRM*MXFRQ
      KWRK1  = KGP    + KZYVAR
      LWRK1  = LWRK   - KWRK1
      !Sonia: extra allocation for projecting out redundant excitation
      !it does not work for MCSCF wfs.
      IF (KZCONF .GT. 0) then
         IF (MCDPRJ) write(lupri,*)'SONIA: RESET MCDPRJ=false for MCSCF'
         MCDPRJ = .false.
      END IF
      if (MCDPRJ) then
        write(lupri,*)'SONIA: allocating for projection'
        KXf   = KWRK1
        KS2Xf = KXf   + KZYVAR
        KWRK1 = KS2Xf + KZYVAR
        LWRK1 = LWRK  - KWRK1
      end if
      !----
C     reuse KGP space for solution vectors.
      KBVECS = KGP
      LGP2   = KZYVAR + 6*N2BASX
      NFRQMX = (LWRK - LGP2 - KBVECS) / KZYVAR
C     ... allocated LGP2 for KGP2 (includes KZYVAR for KWRKE)
      NFRQMX = MIN(NFRQMX,MXFRQ)
      IF (NFRQMX .LT. 1) THEN
         KWRKE  = KBVECS + KZYVAR + LGP2
         CALL ERRWRK('QRLRVE',-KWRKE,LWRK)
      END IF
      KWRKE  = KBVECS + NFRQMX*KZYVAR
      KGP2   = KWRKE
      LGP2   = LWRK   - KGP2
C
      THCRSP = THCLR
      IPRRSP = IPRLR
      MAXIT  = MAXITL
C
C     Call RSPCTL to solve linear set of response equations
C
      DO 680 ILRLAB = 1,NLRLAB
         NFREQ  = 0
         CALL DZERO(WRK(KEIVAL),NLRQR(KSYMOP))
         DO 410 IOP = 1,NLRQR(KSYMOP)
            IOPVEC = ILRQR(KSYMOP) + IOP
            AFRQ = ABS(QRFREQ(IOPVEC))
            IF (QRLBL(IOPVEC) .EQ. LRLAB(ILRLAB)) THEN
               CALL REARSP(LURSP,KLEN,WRK(KBVECS),QRLBL(IOPVEC),
     &                     BLANK,QRFREQ(IOPVEC),D0,KSYMOP,0,THCLR,
     &                     FOUND,CONV,ANTSYM)
               IF (FOUND .AND. CONV) THEN
                  WRITE(LUPRI,'(/A,/A12,A10,/A42)')
     *                 ' Converged solution vector already on file',
     *                 QRLBL(IOPVEC),'freq',
     *                 ' -----------------------------------------'
                  WRITE(LUPRI,'(F22.6)') QRFREQ(IOPVEC)
               ELSE
                  DO IJ = 0, NFREQ
                     IF (ABS(WRK(KEIVAL+IJ)-AFRQ) .LT. THRSML .AND.
     &                    .NOT. AFRQ .LT. THRSML) GOTO 410
                  END DO
                  WRK(KEIVAL+NFREQ) = QRFREQ(IOPVEC)
                  NFREQ = NFREQ + 1
               END IF
            END IF
  410    CONTINUE
         if (mcdprj)
     &   write(lupri,*)'SONIA: we requested NFREQ frequencies=',NFREQ

         IF (NFREQ .EQ. 0) GOTO 680
C
         !SONIA: usual code
         if (.not.(mcdprj.and.LRLAB(ILRLAB)(2:7).eq.'DIPLEN')) then


         CALL GETGPV(LRLAB(ILRLAB),FC,FV,CMO,UDV,PV,XINDX,
     &               ANTSYM,WRK(KGP),LWRK1)
         GPNORM = DNRM2(KZYVAR,WRK(KGP),1)
C
         WRITE (LUPRI,'(//A,I3,/2A,/A,(T25,5F10.6))')
     &   ' QRLRVE -- linear response calculation for symmetry',KSYMOP,
     &   ' QRLRVE -- operator label : ',LRLAB(ILRLAB),
     &   ' QRLRVE -- frequencies :',(WRK(KEIVAL+I),I=0,NFREQ-1)
         IF (GPNORM .LT. THRSML) THEN
            WRITE (LUPRI,*) ' >>> RSPCTL skipped because norm of'
            WRITE (LUPRI,*) '     property vector is only',GPNORM
            CALL DZERO(WRK(KBVECS),KZYVAR)
            DO 450 IFREQ = 1,NFREQ
               CALL WRTRSP(LURSP,KZYVAR,WRK(KBVECS),LRLAB(ILRLAB),BLANK,
     &                     WRK(KEIVAL+IFREQ-1),D0,KSYMOP,0,D0,ANTSYM)
 450        CONTINUE
            GO TO 680
         END IF
C
         KZRED  = 0
         KZYRED = 0
         KEXSIM = NFREQ
         KEXCNV = KEXSIM
         CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *               .TRUE.,LRLAB(ILRLAB),BLANK,WRK(KGP),
     *               WRK(KREDGP),WRK(KREDE),WRK(KREDS),
     *               WRK(KIBTYP),WRK(KEIVAL),WRK(KRESID),WRK(KEIVEC),
     *               XINDX,WRK(KWRK1),LWRK1)
!        CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
!    *               LINEQ,GP,REDGP,REDE,REDS,
!    *               IBTYP,EIVAL,EIVEC,XINDX,WRK,LWRK)
!
         DO 580 IFREQ = 1,NFREQ,NFRQMX
         if (mcdprj)
     &      write(lupri,*)'SONIA: NFREQ=', NFREQ, ' NBX = ', NBX
            NBX   = MIN(NFRQMX,(NFREQ+1-IFREQ))
            IBOFF = IFREQ - 1
            CALL RSPEVE(WRK(KIBTYP),WRK(KEIVAL),WRK(KEIVEC),
     &                  WRK(KBVECS),WRK(KWRKE),NBX,IBOFF)
!           CALL RSPEVE(IBTYP,EIVAL,EIVEC,BVECS,WRK,NBX,IBOFF)
!
            JBVECS = KBVECS
            DO 560 JFREQ = IFREQ,IFREQ-1+NBX
CTODO        IF (IPRRSP.GE.5) THEN
                  WRITE(LUPRI,'(/A,3X,A,3X,A,2X,D12.6,/A,I5)')
     *            ' QRLRVE: SINGLET SOLUTION   LABEL',LRLAB(ILRLAB),
     *            ' FREQUENCY ',WRK(KEIVAL-1+JFREQ),
     *            ' SYMMETRY',KSYMOP
                  JEIVEC = KEIVEC + (JFREQ-1)*KZYRED
               DO JLRLAB = 1, NLRLAB
               CALL GETGPV(LRLAB(JLRLAB),FC,FV,CMO,UDV,PV,XINDX,
     &                     ANTSM2,WRK(KGP2),LGP2)
               VAL = DDOT(KZYVAR,WRK(KGP2),1,WRK(JBVECS),1)
               WRITE (LUPRI,'(/5A,F10.5,A,1P,G22.12)')
     &              '@QRLRVE:  << ', LRLAB(JLRLAB),' ; '
     &              ,LRLAB(ILRLAB),' >> (',WRK(KEIVAL-1+JFREQ),'):',VAL
               END DO
CTODO        END IF
               IF (IPRRSP.GT.100) THEN
                  WRITE (LUPRI,*) 'Column 1 = Z, Column 2 = Y'
                  CALL OUTPUT(WRK(JBVECS),1,KZVAR,1,2,KZVAR,2,1,LUPRI)
               END IF
               CALL WRTRSP(LURSP,KZYVAR,WRK(JBVECS),LRLAB(ILRLAB),
     &                     BLANK,WRK(KEIVAL-1+JFREQ),D0,KSYMOP,0,
     &                     WRK(KRESID - 1 + JFREQ),ANTSYM)
               JBVECS = JBVECS + KZYVAR
  560       CONTINUE
  580    CONTINUE
         else !SONIA: MCDPRJ special code

       !SONIA: MCD, we want to ensure we solve for 1 freq at at time
       !in the resonant case (diplen)
       do ifreq = 1, NFREQ
         WRITE(LUPRI,*)'SONIA: CALL SOLVER 1 freq at a time', ifreq
         !get property gradient vector
         CALL GETGPV(LRLAB(ILRLAB),FC,FV,CMO,UDV,PV,XINDX,
     &                  ANTSYM,WRK(KGP),LWRK1)
         GPNORM = DNRM2(KZYVAR,WRK(KGP),1)
         WRITE (LUPRI,'(//A,I3,/2A,/A,(T25,5F10.6))')
     &   ' QRLRVE PRJ -- linear response calculation for sym',KSYMOP,
     &   ' QRLRVE PRJ-- operator label : ',LRLAB(ILRLAB),
     &   ' QRLRVE PRJ-- frequency :', WRK(KEIVAL+ifreq-1)

         IF (GPNORM .LT. THRSML) THEN
           WRITE (LUPRI,*) ' >>> RSPCTL skipped because norm of'
           WRITE (LUPRI,*) '     property vector is only',GPNORM
           CALL DZERO(WRK(KBVECS),KZYVAR)
           CALL WRTRSP(LURSP,KZYVAR,WRK(KBVECS),LRLAB(ILRLAB),BLANK,
     &                     WRK(KEIVAL+IFREQ-1),D0,KSYMOP,0,D0,ANTSYM)
           GO TO 680
         END IF
         !read in the excitation vector for the resonant excitation
         afrq=WRK(KEIVAL+ifreq-1)
         CALL REARSP(LURSP,KLEN,WRK(KXf),'EXCITLAB',BLANK,
     &               afrq,D0,KSYMOP,0,THCLR,FOUND,CONV,ANTSYM)

         if ((found).and.(conv)) then
           write(lupri,*)'SONIA: found redundant exci X_f on file'
           !CALL RSPPRO(WRK(KXf),KZVAR,-LUPRI)
           call dzero(WRK(KS2Xf),KZYVAR)
           CALL RSPSLI(0,1,DUMMY,WRK(KXf),
     *                 UDV,WRK(KS2Xf),XINDX,WRK(KWRK1),LWRK1)
           write(lupri,*)'SONIA: S2X_f of excitation vector'
           CALL RSPPRO(WRK(KS2Xf),KZVAR,-LUPRI)
        write(lupri,*)'XfS2Xf=',ddot(KZYVAR,WRK(KXf),1,WRK(KS2Xf),1)
           tt = ddot(KZYVAR,WRK(KXf),1,WRK(KGP),1)
           write(lupri,*)'SONIA: ddot{Xf,GD}=', tt
           !write(lupri,*)'SONIA: Property gradient before projection'
           !CALL RSPPRO(WRK(KGP),KZVAR,-LUPRI)
           !project out X_f as: GDp= GD-(ddot{Xf,GD})S[2]X_f
           CALL DAXPY(KZYVAR,-TT,WRK(KS2Xf),1,WRK(KGP),1)
           write(lupri,*)'SONIA: property gradient after projection'
           !CALL RSPPRO(WRK(KGP),KZVAR,-LUPRI)
           tt2 = ddot(KZYVAR,WRK(KXf),1,WRK(KGP),1)
           write(lupri,*)'SONIA: hopefully zero, ddot{Xf|GDp}=', tt2
         end if
         KZRED  = 0
         KZYRED = 0
         KEXSIM = 1
         KEXCNV = 1
         kfreque = KEIVAL+ifreq-1
         !sonia: call response control to start the LR solution...
         !if (.true.) then
         if (abs(tt).le.1.0d-16) then
            !if gradient does not have components along Xf no need
            !to go thru this misery? Maybe not...
            write(lupri,*)'Calling regular RSPCTL,tt=', tt
            CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     &               .TRUE.,LRLAB(ILRLAB),BLANK,WRK(KGP),
     &               WRK(KREDGP),WRK(KREDE),WRK(KREDS),
     &               WRK(KIBTYP),WRK(Kfreque),WRK(KRESID),WRK(KEIVEC),
     &               XINDX,WRK(KWRK1),LWRK1)

         else
            !use the projecting algorithm...
            write(lupri,*)'Calling RSPCTLMCD, tt=',tt
            CALL RSPCTLMCD(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     &               .TRUE.,LRLAB(ILRLAB),BLANK,WRK(KGP),
     &               WRK(KREDGP),WRK(KREDE),WRK(KREDS),
     &               WRK(KIBTYP),WRK(kfreque),WRK(KRESID),WRK(KEIVEC),
     &               XINDX,WRK(KXf),WRK(KS2Xf),WRK(KWRK1),LWRK1)
            write(lupri,*)'Exit projecting response control'
            call flshfo(lupri)
         end if

!         DO 580 IFREQ = 1,NFREQ,NFRQMX
!            NBX   = MIN(NFRQMX,(NFREQ+1-IFREQ))
!            IBOFF = IFREQ - 1
          NBX   = 1 !this is the number of freqs (sets of red space elements)
          IBOFF = 0
           !IBOFF = IFREQ - 1
!          CALL RSPEVE(WRK(KIBTYP),WRK(KEIVAL),WRK(KEIVEC),
!     &                  WRK(KBVECS),WRK(KWRKE),NBX,IBOFF)
          CALL RSPEVE(WRK(KIBTYP),WRK(kfreque),WRK(KEIVEC),
     &                  WRK(KBVECS),WRK(KWRKE),NBX,IBOFF)
          write(lupri,*)'Exit RSPEVE PRJ iboff,ifreq', iboff, ifreq
          call flshfo(lupri)

          JBVECS = KBVECS
          !compute the LR function of the current solution vector
          !DO 560 JFREQ = IFREQ,IFREQ-1+NBX
!TODO     IF (IPRRSP.GE.5) THEN
          JFREQ = IFREQ
          WRITE(LUPRI,'(/A,3X,A,3X,A,2X,D12.6,/A,I5)')
     &    ' QRLRVE: SINGLET SOLUTION   LABEL',LRLAB(ILRLAB),
     &    ' FREQUENCY ',WRK(KEIVAL-1+IFREQ),
     &    ' SYMMETRY',KSYMOP
          CALL FLSHFO(LUPRI)
          JEIVEC = KEIVEC + (JFREQ-1)*KZYRED
          DO JLRLAB = 1, NLRLAB
             CALL GETGPV(LRLAB(JLRLAB),FC,FV,CMO,UDV,PV,XINDX,
     &                     ANTSM2,WRK(KGP2),LGP2)
             VAL = DDOT(KZYVAR,WRK(KGP2),1,WRK(JBVECS),1)
             WRITE (LUPRI,'(/5A,F10.5,A,1P,G22.12)')
     &              '@QRLRVE:  << ', LRLAB(JLRLAB),' ; '
     &              ,LRLAB(ILRLAB),' >> (',WRK(KEIVAL-1+JFREQ),'):',VAL
          END DO
!TODO        END IF
          IF (IPRRSP.GT.100) THEN
            WRITE (LUPRI,*) 'Column 1 = Z, Column 2 = Y'
            CALL OUTPUT(WRK(JBVECS),1,KZVAR,1,2,KZVAR,2,1,LUPRI)
          END IF
          !sonia: I dump the solution vector on file LURSP
          CALL WRTRSP(LURSP,KZYVAR,WRK(JBVECS),LRLAB(ILRLAB),
     &                  BLANK,WRK(KEIVAL-1+JFREQ),D0,KSYMOP,0,
     &                  WRK(KRESID - 1 + JFREQ),ANTSYM)
          JBVECS = JBVECS + KZYVAR
!  560      CONTINUE
!  580    CONTINUE
       end do !ifreq
      end if  !SONIA extra MCDPRJ control...
  680 CONTINUE
C
C *** end of QRLRVE --
C
 9999 CALL QEXIT('QRLRVE')
      RETURN
      END
      SUBROUTINE RIOERR(FOUND,LABEL,FREQ,MSYM)
#include <implicit.h>
#include <priunit.h>
      LOGICAL FOUND, CONV
      CHARACTER*8 LABEL
      IF (.NOT. FOUND) THEN
         WRITE (LUPRI,'(/3A,F12.8,A,I3/A)') 
     &        ' Response label ',LABEL,' with frequency ',
     &        FREQ, ' and symmetry',MSYM,' not found on file LURSP'
         CALL QUIT('Response vector not found on file')
      ELSE
         WRITE (LUPRI,'(/3A,F12.8/A,I3,A)') 
     &        ' @WARNING>>>>'//
     &        ' Response label ',LABEL,' with frequency ', FREQ,
     &        ' and symmetry',MSYM,' not converged on file LURSP'
      END IF
      END
C  /* Deck qrhyp */
      SUBROUTINE QRHYP(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *                 XINDX,VECB,VECC,MJWOP,SPNFC1,SPNFC2,SPNFC,SPNSD1,
     *                 SPNSD2,SPNSD,SPNSU1,SPNSU2,SPNSU,WRK,LWRK)
C
#include <implicit.h>
#include <dummy.h>
#include <iratdef.h>
#include <dftcom.h>
#include <mxcent.h>
C
C PURPOSE:
C CALCULATION OF FIRST HYPERPOLARIZABILITIES AND RELATED PROPERTIES
C
      CHARACTER*1 ALABP, BLABP, CLABP
      CHARACTER*8 BLANK, LABEL
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),VECB(*),VECC(*),WRK(*)
      DIMENSION SPNFC1(*), SPNFC2(*), SPNFC(*), SPNSD1(*), SPNSD2(*),
     &          SPNSD(*), SPNSU1(*), SPNSU2(*), SPNSU(*)
C
      PARAMETER ( D0 = 0.0D0, THRFRQ = 1.0D-14, BLANK = '        ')
C
#include <priunit.h>
#include <infrsp.h>
#include <maxorb.h>
#include <infvar.h>
#include <qrinf.h>
#include <rspprp.h>
#include <indqr.h>
#include <infhyp.h>
#include <inforb.h>
#include <infpri.h>
#include <infspi.h>
#include <wrkrsp.h>
#include <tstjep.h>
#include <infhso.h>
#include <inftap.h>
#include <inflr.h>
C
C     Local variables:
C
      DIMENSION HYPVAL(MXQROP), MJWOP(2,MAXWOP,8)
      LOGICAL   DIPLEN, DOCAL, FOUND, CONV, DFTADX
C
      CALL QENTER('QRHYP')
      IF (SOCOLL) CALL SOINIT
      IF (SSCOLL) CALL SSINIT(SPNFC1,SPNFC2,SPNFC,SPNSD1,SPNSD2,SPNSD,
     &                        SPNSU1,SPNSU2,SPNSU)
      CALL HEADER(' Results from quadratic response calculation',0)
C
      LURSPRES = -1
      CALL GPOPEN(LURSPRES,'RESULTS.RSP','UNKNOWN',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      IPRRSP = IPRHYP
      DO 300 MSYMC = 1,NSYM
         DO 400 MSYMB = 1,NSYM
            MSYMA = MULD2H(MSYMC,MSYMB)
            IF ( (NCQROP(MSYMC) .GT. 0) .AND. (NBQROP(MSYMB) .GT. 0)
     *          .AND. (NAQROP(MSYMA).GT.0) ) THEN
               DO 500 ICOP = 1,NCQROP(MSYMC)
                  DO 600 ICFR = 1,NCQRFR
                     IF ( ISPINC.EQ.1 ) THEN
C                       ICFRNU = NEXLB2 +
C    *                           INQRTR(CQRLB(MSYMC,ICOP),
C    *                                           CQRFR(ICFR),MSYMC)
C                       FREQC = TRFREQ(ICFRNU-NEXLB2-NLRLBL)
C                       LABEL = TRLBL(ICFRNU-NEXLB2)
                        FREQC=CQRFR(ICFR)
                        LABEL=CQRLB(MSYMC,ICOP)
                     ELSE
                        ICFRNU = NEXLB2 + INQRLR(CQRLB(MSYMC,ICOP),
     *                                           CQRFR(ICFR),MSYMC)
                        FREQC = QRFREQ(ICFRNU-NEXLB2)
                        LABEL = QRLBL(ICFRNU-NEXLB2)
                     END IF
                     CALL REARSP(LURSP,KLEN,VECC,LABEL,BLANK,FREQC,D0,
     &                           MSYMC,0,THCLR,FOUND,CONV,ANTSYM)
                     IF (.NOT. (FOUND .AND. CONV))
     &                    CALL RIOERR(FOUND,LABEL, FREQC, MSYMC)
                     IF (FREQC .LT. D0) THEN
                        CALL DSWAP(KLEN/2,VECC,1,VECC(1+KLEN/2),1)
                        IF (ANTSYM .LT. 0.0D0)
     &                       CALL DSCAL(KLEN,ANTSYM,VECC,1)
                     END IF
                     IF (IPRRSP.GT.10) THEN
                        WRITE(LUPRI,'(/A)')'--- C OPERATOR ---'
                        IF ( ISPINC.EQ.1 ) THEN
                           WRITE(LUPRI,'(/A)') ' TRIPLET VECTOR '
                           WRITE(LUPRI,'(/A,3X,A,3X,A,2X,
     *                                    D12.6,/A,I5)')
     *                      ' QRHYP: SOLUTION VECTOR  LABEL',
C    *                      TRLBL(ICFRNU-NEXLB2-NLRLBL),
     *                      LABEL,
C    *                      ' FREQUENCY ',TRFREQ(ICFRNU-NEXLB2-NLRLBL),
     *                      ' FREQUENCY ',FREQC,
     *                      ' SYMMETRY',MSYMC
                        ELSE
                           WRITE(LUPRI,'(/A)') ' SINGLET VECTOR '
                           WRITE(LUPRI,'(/A,3X,A,3X,A,2X,
     *                                    D12.6,/A,I5)')
     *                      ' QRHYP: SOLUTION VECTOR  LABEL',
     *                      QRLBL(ICFRNU-NEXLB2),
     *                      ' FREQUENCY ',QRFREQ(ICFRNU-NEXLB2),
     *                      ' SYMMETRY',MSYMC
                        END IF
                        WRITE (LUPRI,'(/A)')
     *                   ' SOLUTION VECTOR OF LINEAR EQUATIONS'
                        CALL RSPPRC(VECC,MZCONF(MSYMC),
     *                               MZVAR(MSYMC),LUPRI)
                        CALL RSPPRW(VECC(1+MZCONF(MSYMC)),MJWOP,
     *                                MZWOPT(MSYMC),MSYMC,MZVAR(MSYMC),
     *                                LUPRI)
                     END IF
                     DO 700 IBOP = 1,NBQROP(MSYMB)
                        DO 800 IBFR = 1,NBQRFR
C
C  Check for special processes.
C
                           IF ( QRSPEC ) THEN
                              IF (QRPOCK) THEN
                                 IF (ABS(CQRFR(ICFR)).LE.THRFRQ)
     &                                 GO TO 805
                              ENDIF
                              IF (QRSHG) THEN
                                 DIFFRQ = BQRFR(IBFR) - CQRFR(ICFR)
                                 IF (ABS(DIFFRQ).LE.THRFRQ)
     &                                 GO TO 805
                              ENDIF
                              IF (QROPRF) THEN
                                 DIFFRQ = BQRFR(IBFR) + CQRFR(ICFR)
                                 IF (ABS(DIFFRQ).LE.THRFRQ)
     &                                 GO TO 805
                              ENDIF
                              GO TO 800
                           ENDIF
C
 805                       CONTINUE
                           IF ( ISPINB.EQ.1 ) THEN
                              IBFRNU = NEXLB2 +
     *                                 INQRTR(BQRLB(MSYMB,IBOP),
     *                                              BQRFR(IBFR),MSYMB)
!                             FREQB = TRFREQ(IBFRNU-NEXLB2-NLRLBL)
!                             LABEL = TRLBL(IBFRNU-NEXLB2)
                              FREQB=BQRFR(IBFR)
                              LABEL=BQRLB(MSYMB,IBOP)
                           ELSE
                              IBFRNU = NEXLB2 +
     *                                 INQRLR(BQRLB(MSYMB,IBOP),
     *                                              BQRFR(IBFR),MSYMB)
                              LABEL = QRLBL(IBFRNU-NEXLB2)
                              FREQB = QRFREQ(IBFRNU-NEXLB2)
                           END IF
                           BPCFR = CQRFR(ICFR) + BQRFR(IBFR)
                           CALL REARSP(LURSP,KLEN,VECB,LABEL,BLANK,
     &                                 FREQB,D0,MSYMB,0,THCLR,
     &                                 FOUND,CONV,ANTSYM)
                           IF (.NOT. (FOUND .AND. CONV))
     &                          CALL RIOERR(FOUND,LABEL, FREQB, MSYMB)
                           IF (FREQB .LT. D0) THEN
                              CALL DSWAP(KLEN/2,VECB,1,VECB(1+KLEN/2),1)
                              IF (ANTSYM .LT. D0)
     &                             CALL DSCAL(KLEN,ANTSYM,VECB,1)
                           END IF
                           IF (IPRRSP.GT.10) THEN
                              WRITE(LUPRI,'(/A)') ' --- B OPERATOR --- '
                              IF ( ISPINB.EQ.1 ) THEN
                                 WRITE(LUPRI,'(/A)') ' TRIPLET VECTOR '
                                 WRITE(LUPRI,'(/A,3X,A,3X,A,2X,
     *                                        D12.6,/A,I5)')
     *                           ' QRHYP: SOLUTION VECTOR  LABEL',
C    *                             TRLBL(IBFRNU-NEXLB2-NLRLBL),
     *                           LABEL,
     *                           ' FREQUENCY ',
C    *                             TRFREQ(IBFRNU-NEXLB2-NLRLBL),
     *                           FREQB,
     *                           ' SYMMETRY',MSYMB
                              ELSE
                                 WRITE(LUPRI,'(/A)') ' SINGLET VECTOR '
                                 WRITE(LUPRI,'(/A,3X,A,3X,A,2X,
     *                                        D12.6,/A,I5)')
     *                           ' QRHYP: SOLUTION VECTOR  LABEL',
     *                             QRLBL(IBFRNU-NEXLB2),
     *                           ' FREQUENCY ',QRFREQ(IBFRNU-NEXLB2),
     *                           ' SYMMETRY',MSYMB
                              END IF
                              WRITE (LUPRI,'(/A)')
     *                        ' SOLUTION VECTOR OF LINEAR EQUATIONS'
                              CALL RSPPRC(VECB,MZCONF(MSYMB),
     *                                     MZVAR(MSYMB),LUPRI)
                              CALL RSPPRW(VECB(1+MZCONF(MSYMB)),MJWOP,
     *                                    MZWOPT(MSYMB),MSYMB,
     *                                    MZVAR(MSYMB),LUPRI)
                           END IF
                           CALL FLSHFO(LUPRI)
C
C     Find test AVEC for E3TEST and X2TEST
C
                           IAOP   = 1
                           IF ( ISPINA.EQ.0 ) THEN
                              IAFTST = NEXLB2 +
     *                           INQRLR(AQRLB(MSYMA,IAOP),BPCFR,MSYMA)
                           ELSE
                              IAFTST = NEXLB2 + NLRLBL +
     *                           INQRTR(AQRLB(MSYMA,IAOP),BPCFR,MSYMA)
                           END IF
C
                           IF ( RSPCI ) THEN
C                          ... E[3] term is zero for CI
                              CALL DZERO(HYPVAL,NAQROP(MSYMA))
                           ELSE
                              IF (E3TEST) THEN
                                 KAVEC  = 1
                                 KX3WRK = KAVEC + MZYVAR(MSYMA)
                                 LX3WRK = LWRK  - KX3WRK + 1
                                 IF (LX3WRK.LT.0)
     *                             CALL ERRWRK('QRHYP E3',KX3WRK-1,LWRK)
C                                FREQA  = QRFREQ(IAFTST-NEXLB2)
                                 FREQA = BPCFR
                                 CALL REARSP(LURSP,KLEN,WRK(KAVEC),
     &                                AQRLB(MSYMA,IAOP),BLANK,FREQA,D0,
     &                                MSYMA,0,THCLR,FOUND,CONV,ANTSYM)
                                 IF (.NOT. (FOUND .AND. CONV))
     &                              CALL RIOERR(FOUND,AQRLB(MSYMA,IAOP),
     &                                          FREQA,MSYMA)
                                 IF (FREQA .LT. D0) THEN
                                    CALL DSWAP(KLEN/2,WRK(KAVEC),1,
     &                                         WRK(KAVEC+KLEN/2),1)
                                    IF (ANTSYM .LT. D0)
     &                                   CALL DSCAL(KLEN,ANTSYM,
     &                                              WRK(KAVEC),1)
                                 END IF
                              ELSE
                                 KAVEC  = 1
                                 KX3WRK = KAVEC
                                 LX3WRK = LWRK
                              END IF
                              M1 = 1
                              IF (REFCHK) THEN
                                 IF (MSYMB.EQ.IREFSY .AND.
     *                               MSYMC.EQ.IREFSY) THEN
                                    KSYMOP = 1
                                    CALL RSPVAR(UDV,FC,FC,FV,
     *                              FCAC,H2AC,XINDX,WRK(KX3WRK),LX3WRK)
                                    CALL DCOPY(KZWOPT,VECB(1+KZCONF),1,
     *                                         WRK(KX3WRK),1)
                                    CALL DCOPY(KZWOPT,
     *                                         VECB(1+KZCONF+KZVAR),1,
     *                                         WRK(KX3WRK+KZWOPT),1)
                                  IF (ISPINA.EQ.ISPINB) THEN
                                    TRPLET = ISPINA.EQ.1.AND.ISPINB.EQ.1
                                    WRITE(LUPRI,'(/3A)')
     *                                 ' LINEAR TRANSFORMATION WITH',
     *                                 ' RSPELI, ',
     *                                 ' ORBITAL PART OF VECTOR B'
                                    CALL RSPELI(0,1,DUM,WRK(KX3WRK),
     *                                          CMO,UDV,PV,FC,FV,FCAC,
     *                                          H2AC,XINDX,
     *                                          WRK(KX3WRK+KZYWOP),
     *                                          LX3WRK-KZYWOP)
                                  END IF
                                    CALL DCOPY(KZWOPT,VECC(1+KZCONF),1,
     *                                         WRK(KX3WRK),1)
                                    CALL DCOPY(KZWOPT,
     *                                         VECC(1+KZCONF+KZVAR),1,
     *                                         WRK(KX3WRK+KZWOPT),1)
                                  IF (ISPINA.EQ.ISPINC) THEN
                                    TRPLET = ISPINA.EQ.1.AND.ISPINC.EQ.1
                                    WRITE(LUPRI,'(/3A)')
     *                                 ' LINEAR TRANSFORMATION WITH',
     *                                 ' RSPELI, ',
     *                                 ' ORBITAL PART OF VECTOR C'
                                    CALL RSPELI(0,1,DUM,WRK(KX3WRK),
     *                                          CMO,UDV,PV,FC,FV,FCAC,
     *                                          H2AC,XINDX,
     *                                          WRK(KX3WRK+KZYWOP),
     *                                          LX3WRK-KZYWOP)
                                  END IF
                                    CALL GETREF(VECB,KZCONF)
                                    CALL DCOPY(KZCONF,VECB,1,
     &                                         VECB(KZVAR+1),1)
                                    CALL DSCAL(KZCONF,-1.0D0,
     *                                         VECB(KZVAR+1),1)
                                    CALL DCOPY(KZCONF,VECB,1, VECC,1)
                                    CALL DSCAL(KZCONF,-1.0D0, VECC,1)
                                    CALL DCOPY(KZCONF,VECB,1,
     *                                         VECC(KZVAR+1),1)
                                    WRITE(LUPRI,'(/3A)')
     *                               ' CONFIGURATION PART OF VECTOR B',
     *                               ' AND VECTOR C OVERWRITTEN WITH',
     *                               ' REFERENCE VECTOR'
                                 ELSE
                                    CALL QUIT('QRHYP:REFCHK '//
     &                                        'SYMMETRY MISMATCH')
                                 END IF
                              END IF
                              IF (IAABB.NE.0) THEN
                                 IF (IAABB.EQ.1) THEN
                                    NWARN = NWARN + 1
                                    WRITE(LUPRI,'(/A//3A)')
     *                                 ' ****WARNING******',
     *                              'WARNING: ONLY ALFA-ALFA AND',
     *                              ' BETA-BETA COMPONENTS IN',
     *                              ' DENSITY AND CI ROUTINES'
                                 ELSE IF (IAABB.EQ.2) THEN
                                    NWARN = NWARN + 1
                                    WRITE(LUPRI,'(/A//3A)')
     *                                 ' ****WARNING******',
     *                              'WARNING: ONLY ALFA-BETA',
     *                              ' COMPONENTS IN',
     *                              ' DENSITY AND CI ROUTINES'
                                 ELSE
                                    WRITE(LUPRI,'(/A,I12)')
     *                                 ' ****ERROR****** IAABB=',IAABB
                                    CALL QUIT('QRHYP: IAABB ERROR')
                                 END IF
                                 ITSTLP = IAABB
                              END IF
C
C Check if equivalent calculation has been done before,
C DOCAL indicates the result.
C
                              CALL BCCHK(DOCAL,LURSPRES,MSYMA,MSYMB,
     *                             MSYMC,FREQB,FREQC,BQRLB(MSYMB,IBOP),
     *                             CQRLB(MSYMC,ICOP),SPNFC1,SPNFC2,
     *                             SPNFC,SPNSD1,SPNSD2,SPNSD,SPNSU1,
     *                             SPNSU2,SPNSU)
                              IF (.NOT.DOCAL) GO TO 800
C
                              CALL T3DRV(M1,MSYMA,MSYMB,MSYMC,
     *                            VECB,VECC,E3TEST,WRK(KAVEC),
     *                            -FREQB,-FREQC,XINDX,
     *                           UDV,PV,MJWOP,WRK(KX3WRK),LX3WRK,
     &                           CMO,FC,FV)
                              IF (E3TEST) THEN
                                 CALL DCOPY(MZYVAR(MSYMA),WRK(KX3WRK),1,
     *                                      WRK(1),1)
                              END IF
                              KX3BC = 1
                              KAVEC = KX3BC + MZYVAR(MSYMA)
                              CALL FLSHFO(LUPRI)
                              IF (ITSTLP.NE.0)
     *                           CALL QUIT('QRHYP:TSTJEP COMPLETED')
                              DO 910 IAOP = 1,NAQROP(MSYMA)
                                 IF ( ISPINA.EQ.0 ) THEN
                                    IAFRNU = NEXLB2 + INQRLR(
     *                                    AQRLB(MSYMA,IAOP),BPCFR,MSYMA)
                                    LABEL = QRLBL(IAFRNU-NEXLB2)
                                 ELSE
                                    IAFRNU = NEXLB2 + INQRTR(
     *                                    AQRLB(MSYMA,IAOP),BPCFR,MSYMA)
                                    LABEL = TRLBL(IAFRNU-NEXLB2)
                                 ENDIF
C                                FREQA  = QRFREQ(IAFRNU-NEXLB2)
                                 FREQA = BPCFR
                                 CALL REARSP(LURSP,KLEN,WRK(KAVEC),
     &                                LABEL,BLANK,FREQA,D0,MSYMA,0,
     &                                THCLR,FOUND,CONV,ANTSYM)
                                 IF (.NOT. (FOUND .AND. CONV))
     &                             CALL RIOERR(FOUND,LABEL,FREQA,MSYMA)
                                 IF (FREQA .LT. D0) THEN
                                    CALL DSWAP(KLEN/2,WRK(KAVEC),1,
     &                                         WRK(KAVEC+KLEN/2),1)
                                    IF (ANTSYM .LT. D0)
     &                                   CALL DSCAL(KLEN,ANTSYM,
     &                                              WRK(KAVEC),1)
                                 END IF
                                 IF (IPRRSP.GT.10) THEN
                                    WRITE(LUPRI,'(/A)')
     *                                 ' --- A OPERATOR ---'
                                    IF (ISPINA.EQ.0) THEN
                                       WRITE(LUPRI,'(/A)')
     *                                  ' SINGLET VECTOR'
                                    ELSE
                                       WRITE(LUPRI,'(/A)')
     *                                  ' TRIPLET VECTOR'
                                    END IF
                                    WRITE(LUPRI,'(/A,3X,A,3X,A,2X,
     *                                        D12.6,/A,I5)')
     *                              ' QRHYP: SOLUTION VECTOR  LABEL',
     *                              LABEL,' FREQUENCY ',FREQA,
     *                              ' SYMMETRY',MSYMA
                                    WRITE (LUPRI,'(/A)')
     *                              ' LINEAR RESPONSE SOLUTION VECTOR'
                                    CALL RSPPRC(WRK(KAVEC),
     *                                       MZCONF(MSYMA),
     *                                       MZVAR(MSYMA),LUPRI)
                                    CALL RSPPRW(WRK(KAVEC+MZCONF(MSYMA))
     *                                       ,MJWOP,MZWOPT(MSYMA),MSYMA,
     *                                       MZVAR(MSYMA),LUPRI)
                                 END IF
                                 HYPVAL(IAOP) = -DDOT(MZYVAR(MSYMA),
     *                              WRK(KX3BC),1,WRK(KAVEC),1)
                                 IF( IPRRSP .GE. 5) THEN
                                    WRITE(LUPRI,'(/A,2F20.8)')
     *                              '   E3  CONTRIBUTION TO HYPVAL',
     *                              HYPVAL(IAOP),HYPVAL(IAOP)
chjaaj Mar 2005: DFTADX and QDRFT not set
c                                   IF (DFTADX) THEN
c                                      IF (E3TEST.AND.IAOP.EQ.1) THEN
c                                         HYPVAL(1)=HYPVAL(1)+QRDFT
c                                      WRITE(LUPRI,'(/A,2F20.8)')
c    *                              '   DFT CONTRIBUTION TO HYPVAL',
c    *                                  QRDFT,HYPVAL(1)
c                                      END IF
c                                   END IF
                                 END IF
                                 CALL FLSHFO(LUPRI)
  910                         CONTINUE
                           END IF
C
C  B[2] Nc term
C  Calculate B[2] Nc and store in WRK(KB2C)
C  Loop over A operators and calculate Na B[2] Nc term
C
                           KB2C  = 1
                           KAVEC = KB2C + MZYVAR(MSYMA)
                           IF (X2TEST) THEN
                              LABEL = QRLBL(IAFTST-NEXLB2)
                              FREQA  = QRFREQ(IAFTST-NEXLB2)
                              CALL REARSP(LURSP,KLEN,WRK(KAVEC),
     &                             LABEL,BLANK,FREQA,D0,MSYMA,0,
     &                             THCLR,FOUND,CONV,ANTSYM)
                              IF (.NOT. (FOUND .AND. CONV))
     &                          CALL RIOERR(FOUND,LABEL,FREQA,MSYMA)
                              IF (FREQA .LT. D0) THEN
                                 CALL DSWAP(KLEN/2,WRK(KAVEC),1,
     &                                      WRK(KAVEC+KLEN/2),1)
                                 IF (ANTSYM .LT. D0)
     &                                CALL DSCAL(KLEN,ANTSYM,
     &                                WRK(KAVEC),1)
                              END IF
                              KFREE = KAVEC + MZYVAR(MSYMA)
                           ELSE
                              KFREE = KAVEC
                           END IF
                           LFREE = LWRK - KFREE + 1
                           IF (LFREE.LT.0)
     *                        CALL ERRWRK('QRHYP B2C',KFREE-1,LWRK)
                           IF (X2MAT) THEN
                              IPRRSP = IPRHYP
                              CALL
     *                        X2EXPL(1,MZYVAR(MSYMA),MZYVAR(MSYMC),
     *                              MSYMA,ISPINA,MSYMC,ISPINC,
     *                              WRK(KAVEC),VECC,WRK(KB2C),XINDX,
     *                              UDV,PV,
     *                              BQRLB(MSYMB,IBOP),MSYMB,ISPINB,
     *                              CMO,MJWOP,WRK(KFREE),LFREE)
                              CALL QUIT('X2EXPL TEST')
                           END IF
                           CALL X2INIT(1,MZYVAR(MSYMA),MZYVAR(MSYMC),
     *                              MSYMA,ISPINA,MSYMC,ISPINC,
     *                              WRK(KAVEC),VECC,WRK(KB2C),XINDX,
     *                              UDV,PV,
     *                              BQRLB(MSYMB,IBOP),MSYMB,ISPINB,
     *                              CMO,MJWOP,WRK(KFREE),LFREE)
                           DO 920 IAOP = 1,NAQROP(MSYMA)
                              IF ( ISPINA.EQ.0 ) THEN
                                 IAFRNU = NEXLB2 + INQRLR(
     *                                 AQRLB(MSYMA,IAOP),BPCFR,MSYMA)
                                 LABEL = QRLBL(IAFRNU-NEXLB2)
                              ELSE
                                 IAFRNU = NEXLB2 + INQRTR(
     *                                 AQRLB(MSYMA,IAOP),BPCFR,MSYMA)
                                 LABEL = TRLBL(IAFRNU-NEXLB2)
                              ENDIF
C                             FREQA  = QRFREQ(IAFRNU-NEXLB2)
                              FREQA = BPCFR
                              CALL REARSP(LURSP,KLEN,WRK(KAVEC),
     &                             LABEL,BLANK,FREQA,D0,MSYMA,0,
     &                             THCLR,FOUND,CONV,ANTSYM)
                              IF (.NOT. (FOUND .OR. CONV))
     &                             CALL RIOERR(FOUND,LABEL,FREQA,MSYMA)
                              IF (FREQA .LT. D0) THEN
                                 CALL DSWAP(KLEN/2,WRK(KAVEC),1,
     &                                      WRK(KAVEC+KLEN/2),1)
                                 IF (ANTSYM .LT. D0)
     &                              CALL DSCAL(KLEN,ANTSYM,
     &                                         WRK(KAVEC),1)
                              END IF
                              IF (IPRRSP.GT.10) THEN
                                 WRITE(LUPRI,'(/A)')
     *                              ' --- A OPERATOR ---'
                                 IF (ISPINA.EQ.0) THEN
                                    WRITE(LUPRI,'(/A)')
     *                               ' SINGLET VECTOR'
                                 ELSE
                                    WRITE(LUPRI,'(/A)')
     *                               ' TRIPLET VECTOR'
                                 END IF
                                 WRITE(LUPRI,'(/A,3X,A,3X,A,2X,
     *                                     D12.6,/A,I5)')
     *                           ' QRHYP: SOLUTION VECTOR  LABEL',
     *                           LABEL,' FREQUENCY ',FREQA,
     *                           ' SYMMETRY',MSYMA
                                 IF (RSPCI) THEN
                                    WRITE (LUPRI,'(/A)')
     *                              ' LINEAR RESPONSE SOLUTION VECTOR'
                                    CALL RSPPRC(WRK(KAVEC),
     *                                       MZCONF(MSYMA),
     *                                       MZVAR(MSYMA),LUPRI)
                                 ELSE
                                    WRITE (LUPRI,'(/A)')
     *                   ' See linear response solution vector above'
                                 END IF
                              END IF
                              VAL = DDOT(MZYVAR(MSYMA),
     *                                   WRK(KB2C),1,WRK(KAVEC),1)
                              HYPVAL(IAOP) = HYPVAL(IAOP) + VAL
                              IF( IPRRSP .GE. 5) THEN
                                 WRITE(LUPRI,'(/A,2F20.8)')
     *                           ' + B2C CONTRIBUTION TO HYPVAL',
     *                           VAL,HYPVAL(IAOP)
                              END IF
                              CALL FLSHFO(LUPRI)
  920                      CONTINUE
C
C  C[2] Nb term
C  Calculate C[2] Nb and store in WRK(KC2B)
C  Loop over A operators and calculate Na C[2] Nb term
C
                           KC2B  = 1
                           KAVEC = KC2B + MZYVAR(MSYMA)
                           IF (X2TEST) THEN
                              LABEL = QRLBL(IAFTST-NEXLB2)
                              FREQA  = QRFREQ(IAFTST-NEXLB2)
                              CALL REARSP(LURSP,KLEN,WRK(KAVEC),
     &                             LABEL,BLANK,FREQA,D0,MSYMA,0,
     &                             THCLR,FOUND,CONV,ANTSYM)
                              IF (.NOT. (FOUND .AND. CONV))
     &                             CALL RIOERR(FOUND,LABEL,FREQA,MSYMA)
                              IF (FREQA .LT. D0) THEN
                                 CALL DSWAP(KLEN/2,WRK(KAVEC),1,
     &                                      WRK(KAVEC+KLEN/2),1)
                                 IF (ANTSYM .LT. D0)
     &                                CALL DSCAL(KLEN,ANTSYM,
     &                                           WRK(KAVEC),1)
                              END IF
                              KFREE = KAVEC + MZYVAR(MSYMA)
                           ELSE
                              KFREE = KAVEC
                           END IF
                           LFREE = LWRK - KFREE + 1
                           IF (LFREE.LT.0)
     *                        CALL ERRWRK('QRHYP C2B',KFREE-1,LWRK)
                           IF (X2MAT) THEN
                              IPRRSP = IPRHYP
                              CALL X2EXPL(1,MZYVAR(MSYMA),MZYVAR(MSYMB),
     *                              MSYMA,ISPINA,MSYMB,ISPINB,
     *                              WRK(KAVEC),VECB,WRK(KC2B),XINDX,
     *                              UDV,PV,
     *                              CQRLB(MSYMC,ICOP),MSYMC,ISPINC,
     *                              CMO,MJWOP,WRK(KFREE),LFREE)
                              CALL QUIT('X2EXPL TEST')
                           END IF
                           CALL X2INIT(1,MZYVAR(MSYMA),MZYVAR(MSYMB),
     *                           MSYMA,ISPINA,MSYMB,ISPINB,
     *                           WRK(KAVEC),VECB,WRK(KC2B),XINDX,
     *                           UDV,PV,
     *                           CQRLB(MSYMC,ICOP),MSYMC,ISPINC,
     *                           CMO,MJWOP,WRK(KFREE),LFREE)
                           DO 930 IAOP = 1,NAQROP(MSYMA)
                              IF ( ISPINA.EQ.0 ) THEN
                                 IAFRNU = NEXLB2 + INQRLR(
     *                                 AQRLB(MSYMA,IAOP),BPCFR,MSYMA)
                                 LABEL = QRLBL(IAFRNU-NEXLB2)
                              ELSE
                                 IAFRNU = NEXLB2 + INQRTR(
     *                                 AQRLB(MSYMA,IAOP),BPCFR,MSYMA)
                                 LABEL = TRLBL(IAFRNU-NEXLB2)
                              ENDIF
C                             FREQA  = QRFREQ(IAFRNU-NEXLB2)
                              FREQA  = BPCFR
                              CALL REARSP(LURSP,KLEN,WRK(KAVEC),
     &                             LABEL,BLANK,FREQA,D0,MSYMA,0,
     &                             THCLR,FOUND,CONV,ANTSYM)
                              IF (.NOT. (FOUND .AND. CONV))
     &                             CALL RIOERR(FOUND,LABEL,FREQA,MSYMA)
                              IF (FREQA .LT. D0) THEN
                                 CALL DSWAP(KLEN/2,WRK(KAVEC),1,
     &                                      WRK(KAVEC+KLEN/2),1)
                                 IF (ANTSYM .LT. D0) CALL DSCAL(KLEN,
     &                                ANTSYM,WRK(KAVEC),1)
                              END IF
                              IF (IPRRSP.GE.10) THEN
                                 WRITE(LUPRI,'(/A)')
     *                              ' --- A OPERATOR ---'
                                 IF (ISPINA.EQ.0) THEN
                                    WRITE(LUPRI,'(/A)')
     *                               ' SINGLET VECTOR'
                                 ELSE
                                    WRITE(LUPRI,'(/A)')
     *                               ' TRIPLET VECTOR'
                                 END IF
                                 WRITE(LUPRI,'(/A,3X,A,3X,A,2X,
     *                                     D12.6,/A,I5)')
     *                           ' QRHYP: SOLUTION VECTOR  LABEL',
     *                           LABEL,' FREQUENCY ',FREQA,
     *                           ' SYMMETRY',MSYMA
                                 WRITE (LUPRI,'(/A)')
     *                   ' See linear response solution vector above'
                              END IF
                              VAL = DDOT(MZYVAR(MSYMA),
     *                                   WRK(KC2B),1,WRK(KAVEC),1)
                              HYPVAL(IAOP) = HYPVAL(IAOP) + VAL
                              IF( IPRRSP .GE. 5) THEN
                                 WRITE(LUPRI,'(/A,2F20.8)')
     *                           ' + C2B CONTRIBUTION TO HYPVAL',
     *                           VAL,HYPVAL(IAOP)
                              END IF
                              CALL FLSHFO(LUPRI)
  930                      CONTINUE
C
C A2TEST
C
                           DO 940 IAOP = 1,NAQROP(MSYMA)
                              KA2B = 1
                              KFREE = KA2B + MZYVAR(MSYMC)
                              LFREE = LWRK - KFREE + 1
                              IF (LFREE.LT.0)
     *                           CALL ERRWRK('QRHYP A2B',KFREE-1,LWRK)
                              CALL A2INIT(1,MZYVAR(MSYMC),
     *                              MZYVAR(MSYMB),
     *                              MSYMC,ISPINC,MSYMB,ISPINB,
     *                              VECC,VECB,WRK(KA2B),XINDX,
     *                              UDV,PV,
     *                              AQRLB(MSYMA,IAOP),MSYMA,ISPINA,
     *                              CMO,MJWOP,WRK(KFREE),LFREE)
                              VAL = DDOT(MZYVAR(MSYMC),
     *                                   WRK(KA2B),1,VECC,1)
                              HYPVAL(IAOP) = HYPVAL(IAOP) + VAL
                              IF( IPRRSP .GE. 5) THEN
                                 WRITE(LUPRI,'(/A,2F20.8)')
     *                           ' + A2B CONTRIBUTION TO HYPVAL',
     *                           VAL,HYPVAL(IAOP)
                              END IF
                              CALL FLSHFO(LUPRI)
C
                              KA2C = 1
                              KFREE = KA2C + MZYVAR(MSYMB)
                              LFREE = LWRK - KFREE
                              IF (LFREE.LT.0)
     *                           CALL ERRWRK('QRHYP A2C',KFREE-1,LWRK)
                              CALL A2INIT(1,MZYVAR(MSYMB),MZYVAR(MSYMC),
     *                              MSYMB,ISPINB,MSYMC,ISPINC,
     *                              VECB,VECC,WRK(KA2C),XINDX,
     *                              UDV,PV,
     *                              AQRLB(MSYMA,IAOP),MSYMA,ISPINA,
     *                              CMO,MJWOP,WRK(KFREE),LFREE)
                              VAL = DDOT(MZYVAR(MSYMB),
     *                                   WRK(KA2C),1,VECB,1)
                              HYPVAL(IAOP) = HYPVAL(IAOP) + VAL
                              IF( IPRRSP .GE. 5) THEN
                                 WRITE(LUPRI,'(/A,2F20.8)')
     *                           ' + A2C CONTRIBUTION TO HYPVAL',
     *                           VAL,HYPVAL(IAOP)
                              END IF
                              CALL FLSHFO(LUPRI)
 940                       CONTINUE
                           DO 950 IAOP = 1,NAQROP(MSYMA)
C
                              DIPLEN = .FALSE.
                              IF (AQRLB(MSYMA,IAOP)(2:7).EQ.'DIPLEN'
     &                      .AND. BQRLB(MSYMB,IBOP)(2:7).EQ.'DIPLEN'
     &                      .AND. CQRLB(MSYMC,ICOP)(2:7).EQ.'DIPLEN')
     &                           THEN
                                 DIPLEN = .TRUE.
                                 ALABP = AQRLB(MSYMA,IAOP)(1:1)
                                 BLABP = BQRLB(MSYMB,IBOP)(1:1)
                                 CLABP = CQRLB(MSYMC,ICOP)(1:1)
                              WRITE(LUPRI,'(2(A,F9.6),A,6(A1),A,F16.8)')
     &              '@ B-freq =',BQRFR(IBFR),'  C-freq =',CQRFR(ICFR),
     &                   '     beta(',ALABP,';',BLABP,',',CLABP,')',
     &                        ' =',-HYPVAL(IAOP)
                           WRITE(LURSPRES,'(2(A,F9.6),A,6(A1),A,F16.8)')
     &              '@ B-freq =',BQRFR(IBFR),'  C-freq =',CQRFR(ICFR),
     &                   '     beta(',ALABP,';',BLABP,',',CLABP,')',
     &                        ' =',-HYPVAL(IAOP)
C      hjaaj Oct 2001: note that
C            beta = <<mu; -mu, -mu>> = -<<r; r, r>>
                              END IF
C
                              IF (IPRRSP.GE.3 .OR. .NOT.DIPLEN) THEN
                                 WRITE(LUPRI,'(//A,3(/2A,2I5))')
     * '@ Quadratic response function value in a.u. for',
     * '@ A operator, symmetry, spin: ',AQRLB(MSYMA,IAOP),MSYMA,ISPINA,
     * '@ B operator, symmetry, spin: ',BQRLB(MSYMB,IBOP),MSYMB,ISPINB,
     * '@ C operator, symmetry, spin: ',CQRLB(MSYMC,ICOP),MSYMC,ISPINC
                                 WRITE(LUPRI,'(/A,3F15.8)')
     *                           '@omega B, omega C, QR value :',
     *                           BQRFR(IBFR),CQRFR(ICFR),HYPVAL(IAOP)
                                 WRITE(LURSPRES,'(//A,3(/2A,2I5))')
     * '@ Quadratic response function value in a.u. for',
     * '@ A operator, symmetry, spin: ',AQRLB(MSYMA,IAOP),MSYMA,ISPINA,
     * '@ B operator, symmetry, spin: ',BQRLB(MSYMB,IBOP),MSYMB,ISPINB,
     * '@ C operator, symmetry, spin: ',CQRLB(MSYMC,ICOP),MSYMC,ISPINC
                                 WRITE(LURSPRES,'(/A,3F15.8)')
     *                           '@omega B, omega C, QR value :',
     *                           BQRFR(IBFR),CQRFR(ICFR),HYPVAL(IAOP)
                              END IF
                             IF (SOCOLL) CALL SODIST(AQRLB(MSYMA,IAOP),
     *                             MSYMA,BQRLB(MSYMB,IBOP),MSYMB,
     *                             CQRLB(MSYMC,ICOP),MSYMC,HYPVAL(IAOP))
                             IF (SSCOLL) CALL SSDIST(AQRLB(MSYMA,IAOP),
     *                             MSYMA,BQRLB(MSYMB,IBOP),MSYMB,
     *                             CQRLB(MSYMC,ICOP),MSYMC,HYPVAL(IAOP),
     *                             SPNFC1,SPNFC2,SPNFC,SPNSD1,SPNSD2,
     *                             SPNSD,SPNSU1,SPNSU2,SPNSU)
 950                       CONTINUE
 800                    CONTINUE
 700                 CONTINUE
 600              CONTINUE
 500           CONTINUE
            END IF
 400     CONTINUE
 300  CONTINUE
      IF (SOCOLL) THEN 
         KMAT1 = 1
         KMAT2 = KMAT1 + 3*MXCOOR
         KCSTRA = KMAT2 + MXCOOR*MXCOOR
         KSCTRA = KCSTRA + MXCOOR*MXCOOR
         CALL SOSOUT(WRK(KMAT1),WRK(KMAT2),WRK(KCSTRA),WRK(KSCTRA))
      END IF
      IF (SSCOLL) THEN 
         KCSTRA = 1
         KSCTRA = KCSTRA + MXCOOR*MXCOOR
         CALL SPNSO(SPNFC1,SPNFC2,SPNFC,SPNSD1,SPNSD2,SPNSD,SPNSU1,
     &               SPNSU2,SPNSU,WRK(KCSTRA),WRK(KSCTRA))
      END IF
      CALL GPCLOSE(LURSPRES,'KEEP')
      CALL QEXIT('QRHYP')
      RETURN
      END
C  /* Deck qrsmo */
      SUBROUTINE QRSMO(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *                  XINDX,VECB,VECC,MJWOP,WRK,LWRK)
C
C PURPOSE:
C CALCULATION OF SECOND ORDER TRANSITION MOMENTS
C
C Revision 910405-hjaaj: FACB and FACMOM for mixed and velocity repr.
C
C Modified 980722-sonia: MCD. Included common block WRKRSP and
C call to GETGPV.
C
#include <implicit.h>
C
      LOGICAL FOUND, CONV
      CHARACTER*8 BLANK, LABEL
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),VECB(*),VECC(*),WRK(*)
C
      PARAMETER (BLANK = '        ')
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, D4D3 = 4.0D0/3.0D0,
     &            IMXEXC=20, D30 = 30.0D0)
      PARAMETER ( DP5 = 0.5D0 , D3 = 3.0D0)
      PARAMETER (AU_TO_GM=80*3.1415926D0**3 / 137.036**2 * 
     &     0.591177**4 *2.418884)
      DIMENSION RESSOM(3,3,IMXEXC,8), RESFRE(IMXEXC,8)
C
#include <iratdef.h>
C
#include <rspprp.h>
#include <infrsp.h>
#include <infsmo.h>
#include <indqr.h>
#include <inforb.h>
#include <maxorb.h>
#include <priunit.h>
#include <infvar.h>
#include <qrinf.h>
      DIMENSION MJWOP(2,MAXWOP,8)
#include <infpp.h>
#include <infpri.h>
#include <infspi.h>
#include <codata.h>
#include <infhso.h>
#include <inftap.h>
#include <inflr.h>
C Sonia: MCD variables
#include <wrkrsp.h>
C
      DIMENSION JWOPSV(2,MAXWOP)
      CHARACTER A1*1, B1*1, CLABEL*8
      LOGICAL LYESGP
C
      IPRRSP = IPRSMO
      LYESGP = .FALSE.
C
C end lvariables for MCD
C
      CALL QENTER('QRSMO')
C
C
C Counters NCEXC and ICEXC are defined in order to identify the correct 
C response calculations in a two-photon process. RESSOM and RESFRE are 
C temporary storages for the result in order to print out nicely after the 
C calculation.
C
      IF (TWOPHO) THEN
         CALL DZERO(RESSOM,3*3*IMXEXC*8)
         CALL DZERO(RESFRE,IMXEXC*8)
      END IF
      DO 300 MSYMC = 1,NSYM
         MSYMCX = MULD2H(IREFSY,MSYMC)
         IF (PHOSPH) CALL DZERO(PHOSMAT,9*NSMCNV(MSYMC))
         IF (MNFPHO) CALL DZERO(PHOSMN,9*NSMCNV(MSYMC))
         DO 400 MSYMB = 1,NSYM
            MSYMA = MULD2H(MSYMC,MSYMB)
            IF ( (NSMCNV(MSYMC) .GT. 0) .AND. (NBSMOP(MSYMB) .GT. 0)
     *          .AND. (NASMOP(MSYMA).GT.0) ) THEN
               NCEXC = 0
               DO I=1,MSYMC-1
                  NCEXC = NCEXC + NSMCNV(I)
               END DO
               DO 500 IC = 1,NSMCNV(MSYMC)
                  ICEXC = NCEXC + IC
                  ICEXNU = INQREX(MSYMC,IC)
                  FREQC = EXCITA(MSYMC,IC,ISPINC+1)
                  CALL REARSP(LURSP,KLEN,VECC,'EXCITLAB',BLANK,
     &                        FREQC,D0,MSYMC,IC,THCPP,FOUND,CONV,ANTSYM)
                  IF (.NOT. (FOUND .AND. CONV))
     &                 CALL RIOERR(FOUND,'EXCITLAB',FREQC,MSYMC)
                  IF (FREQC .LT. D0) THEN
                     CALL DSWAP(KLEN/2,VECC,1,VECC(1+KLEN/2),1)
                     IF (ANTSYM .LT. D0) CALL DSCAL(KLEN,ANTSYM,VECC,1)
                  END IF
                  IF (IPRRSP.GT.10) THEN
                     WRITE(LUPRI,'(/A)')
     *                  ' --- C OPERATOR , EXCITATION ENERGY'
                     IF ( ISPINC.EQ.1 ) THEN
                        WRITE(LUPRI,'(/A)') ' TRIPLET VECTOR '
                     ELSE
                        WRITE(LUPRI,'(/A)') ' SINGLET VECTOR '
                     END IF
                     WRITE(LUPRI,'(/A,I5,3X,A,I5,/A,2X,D12.6)')
     *               ' QRSMO: EIGENVECTOR ',ICEXNU,
     *               ' SYMMETRY',MSYMC,' EXCITATION ENERGY',FREQC
                     CALL RSPPRC(VECC,MZCONF(MSYMC),
     *                                  MZVAR(MSYMC),LUPRI)
                     CALL RSPPRW(VECC(1+MZCONF(MSYMC)),MJWOP,
     *                                  MZWOPT(MSYMC),MSYMC,
     *                                  MZVAR(MSYMC),LUPRI)
                  END IF
                  DO 700 IBOP = 1,NBSMOP(MSYMB)
                     DO 800 IBFR = 1,NBSMFR
C
C If two-photon process then the B-frequency should equal half the
C excitation energy.
C
                        IF (TWOPHO .AND. IBFR.NE.ICEXC) GO TO 800
                        IF (ISPINB.EQ.1) THEN
                           IBFRNU = NEXLB2 + NLRLBL +
     *                              INQRTR(BSMLB(MSYMB,IBOP),
     *                                            -BSMFR(IBFR),MSYMB)
C                          LABEL  = TRLBL(IBFRNU-NEXLB2)
                           LABEL=BSMLB(MSYMB,IBOP)
C                          FREQB  = TRFREQ(IBFRNU-NEXLB2-NLRLBL)
                           FREQB=BSMFR(IBFR)
                        ELSE
                           IBFRNU = NEXLB2 + INQRLR(BSMLB(MSYMB,IBOP),
     *                                            -BSMFR(IBFR),MSYMB)
                           LABEL  = QRLBL(IBFRNU-NEXLB2)
                           FREQB  = QRFREQ(IBFRNU-NEXLB2)
                        END IF
                        IF (INDEX(BSMLB(MSYMB,IBOP),'DIPVEL') .NE. 0)
     *                     THEN
                           FACB = D1 / BSMFR(IBFR)
                        ELSE
                           FACB = D1
                        END IF
                        CALL REARSP(LURSP,KLEN,VECB,LABEL,BLANK,
     &                              FREQB,D0,MSYMB,0,THCLR,FOUND,
     &                              CONV,ANTSYM)
                        IF (.NOT. (FOUND .AND. CONV))
     &                       CALL RIOERR(FOUND,LABEL,FREQB,MSYMB)
                        IF (FREQB .LT. D0) THEN
                           CALL DSWAP(KLEN/2,VECB,1,VECB(1+KLEN/2),1)
                           IF (ANTSYM .LT. D0) CALL DSCAL(KLEN,ANTSYM,
     &                                                    VECB,1)
                        END IF
                        IF (IPRRSP.GT.10) THEN
                           WRITE(LUPRI,'(/A)')' --- B OPERATOR ---'
                           IF (ISPINB.EQ.1) THEN
                              WRITE(LUPRI,'(/A)') ' TRIPLET VECTOR '
                              WRITE(LUPRI,'(/A,3X,A,3X,A,2X,
     *                                     D12.6,/A,I5)')
     *                        ' QRSMO: SOLUTION VECTOR  LABEL',
     *                          TRLBL(IBFRNU-NEXLB2-NLRLBL),
C    *                      ' FREQUENCY ',TRFREQ(IBFRNU-NEXLB2-NLRLBL),
     *                      ' FREQUENCY ',FREQB,
     *                        ' SYMMETRY',MSYMB
                           ELSE
                              WRITE(LUPRI,'(/A)') ' SINGLET VECTOR '
                              WRITE(LUPRI,'(/A,3X,A,3X,A,2X,
     *                                     D12.6,/A,I5)')
     *                        ' QRSMO: SOLUTION VECTOR  LABEL',
     *                          QRLBL(IBFRNU-NEXLB2),
     *                        ' FREQUENCY ',QRFREQ(IBFRNU-NEXLB2),
     *                        ' SYMMETRY',MSYMB
                           END IF
                           WRITE (LUPRI,'(/A)')
     *                     ' SOLUTION VECTOR OF LINEAR EQUATIONS'
                           CALL RSPPRC(VECB,MZCONF(MSYMB),
     *                                  MZVAR(MSYMB),LUPRI)
                           CALL RSPPRW(VECB(1+MZCONF(MSYMB)),MJWOP,
     *                                  MZWOPT(MSYMB),MSYMB,
     *                                  MZVAR(MSYMB),LUPRI)
                        END IF
                        CALL FLSHFO(LUPRI)
C
C If two-photon process then the A-frequency should equal half the
C excitation energy.
C
                        IF (TWOPHO) THEN 
                           CEXMBF=BSMFR(IBFR)
                        ELSE
                           CEXMBF=EXCITA(MSYMC,IC,ISPINC+1)-BSMFR(IBFR)
                        END IF 
                        IF(.NOT.RSPCI) THEN
                           IF (E3TEST) THEN
                              IAOP   = 1
                              KAVEC  = 1
                              KX3WRK = KAVEC + MZYVAR(MSYMA)
                              LX3WRK = LWRK  - KX3WRK + 1
                              IF (LX3WRK.LT.0)
     *                           CALL ERRWRK('QRSMO 1',KX3WRK-1,LWRK)
                              IF (ISPINA.EQ.0) THEN
                                 IAFRNU = NEXLB2 +
     *                           INQRLR(ASMLB(MSYMA,IAOP),CEXMBF,MSYMA)
                              ELSE
                                 IAFRNU = NEXLB2 + NLRLBL +
     *                           INQRTR(ASMLB(MSYMA,IAOP),CEXMBF,MSYMA)
                              END IF
                              LABEL = QRLBL(IAFRNU-NEXLB2)
                              FREQA = QRFREQ(IAFRNU-NEXLB2)
                              CALL REARSP(LURSP,KLEN,WRK(KAVEC),
     &                             LABEL,BLANK,FREQA,D0,MSYMA,0,
     &                             THCLR,FOUND,CONV,ANTSYM)
                              IF (.NOT. (FOUND .AND. CONV))
     &                             CALL RIOERR(FOUND,LABEL,FREQA,MSYMA)
                              IF (FREQA .LT. D0) THEN
                                 CALL DSWAP(KLEN/2,WRK(KAVEC),1,
     &                                      WRK(KAVEC+KLEN/2),1)
                                 IF (ANTSYM .LT. D0) CALL DSCAL(KLEN,
     &                                ANTSYM,WRK(KAVEC),1)
                              END IF
                           ELSE
                              KAVEC  = 1
                              KX3WRK = KAVEC
                              LX3WRK = LWRK
                           END IF
                           M1 = 1
C
                           CALL T3DRV(M1,MSYMA,MSYMB,MSYMC,VECB,VECC,
     *                         E3TEST,WRK(KAVEC),
     *                         BSMFR(IBFR),-EXCITA(MSYMC,IC,ISPINC+1),
     *                         XINDX,UDV,PV,MJWOP,
     *                         WRK(KX3WRK),LX3WRK,CMO,FC,FV)
                           IF (E3TEST) THEN
                              CALL DCOPY(MZYVAR(MSYMA),WRK(KX3WRK),1,
     *                                   WRK(1),1)
                           END IF
                           CALL FLSHFO(LUPRI)
                        END IF
                        DO 900 IAOP = 1,NASMOP(MSYMA)
                           IF ( ISPINA.EQ.0 ) THEN
                              IAFRNU = NEXLB2 +
     *                           INQRLR(ASMLB(MSYMA,IAOP),CEXMBF,MSYMA)
                           ELSE
                              IAFRNU = NEXLB2 + NLRLBL +
     *                           INQRTR(ASMLB(MSYMA,IAOP),CEXMBF,MSYMA)
                           END IF
                           LABEL = QRLBL(IAFRNU-NEXLB2)
                           FREQA = QRFREQ(IAFRNU-NEXLB2)
                           IF (INDEX(ASMLB(MSYMA,IAOP),'DIPVEL') .NE. 0)
     *                     THEN
                              FACMOM = FACB / CEXMBF
                           ELSE
                              FACMOM = FACB
                           END IF
                           IF ( RSPCI) THEN
                              SMOM = D0
                           ELSE
                              CALL REARSP(LURSP,KLEN,
     &                             WRK(1+MZYVAR(MSYMA)),LABEL,BLANK,
     &                             FREQA,D0,MSYMA,0,THCLR,FOUND,
     &                             CONV,ANTSYM)
                              IF (.NOT. (FOUND .AND. CONV))
     &                             CALL RIOERR(FOUND,LABEL,FREQA,MSYMA)
                              IF (FREQA .LT. D0) THEN
                                 CALL DSWAP(KLEN/2,WRK(1+MZYVAR(MSYMA)),
     &                                1,WRK(1+MZYVAR(MSYMA)+KLEN/2),1)
                                 IF (ANTSYM .LT. D0) CALL DSCAL(KLEN,
     &                                ANTSYM,WRK(1+MZYVAR(MSYMA)),1)
                              END IF
                              IF (IPRRSP.GT.10) THEN
                                 IF ( ISPINA.EQ.0) THEN
                                    WRITE(LUPRI,'(/A//A)')
     &                              ' --- A OPERATOR ---',
     *                              ' SINGLET SOLUTION VECTOR'
                                 ELSE
                                    WRITE(LUPRI,'(/A//A)')
     &                              ' --- A OPERATOR ---',
     *                              ' TRIPLET SOLUTION VECTOR'
                                 END IF
                                 WRITE(LUPRI,'(/A,3X,A,3X,A,2X,
     *                                        D12.6,/A,I5)')
     *                           ' QRSMO: SOLUTION VECTOR  LABEL',
     *                           LABEL, ' FREQUENCY ',FREQA,
     *                           ' SYMMETRY',MSYMA
                                 WRITE (LUPRI,'(/A)')
     *                           ' SOLUTION VECTOR OF LINEAR EQUATIONS'
                                CALL RSPPRC(WRK(MZYVAR(MSYMA)+1),
     *                                   MZCONF(MSYMA),
     *                                   MZVAR(MSYMA),LUPRI)
                                IWRKOF = MZYVAR(MSYMA)+MZCONF(MSYMA)+1
                                CALL RSPPRW(WRK(IWRKOF),MJWOP,
     *                                   MZWOPT(MSYMA),MSYMA,
     *                                   MZVAR(MSYMA),LUPRI)
                              END IF
                              SMOM = FACMOM*DDOT(MZYVAR(MSYMA),WRK(1),1,
     *                                 WRK(1+MZYVAR(MSYMA)),1)
                              IF( IPRRSP .GE. 5) THEN
                                 WRITE(LUPRI,'(/A,2F20.8)')
     *                           ' E3             CONTRIBUTION TO SMOM'
     *                           ,SMOM,SMOM
                              END IF
                              CALL FLSHFO(LUPRI)
                           END IF
C
C X2TEST
C
                           IF (RSPCI) THEN
                              KAVEC = 1
                           ELSE
                              KAVEC = 1 + MZYVAR(MSYMA)
                           END IF
                           IF (X2TEST) THEN
                              KB2C = KAVEC + MZYVAR(MSYMA)
                              CALL REARSP(LURSP,KLEN,WRK(KAVEC),
     &                             LABEL,BLANK,FREQA,D0,MSYMA,0,
     &                             THCLR,FOUND,CONV,ANTSYM)
                              IF (.NOT. (FOUND .AND. CONV))
     &                             CALL RIOERR(FOUND,LABEL,FREQA,MSYMA)
                              IF (FREQA .LT. D0) THEN
                                 CALL DSWAP(KLEN/2,WRK(KAVEC),1,
     &                                      WRK(KAVEC+KLEN/2),1)
                                 IF (ANTSYM .LT. D0) CALL DSCAL(KLEN,
     &                                ANTSYM,WRK(KAVEC),1)
                              END IF
                           ELSE
                              KB2C = KAVEC
                           END IF
                           KFREE = KB2C + MZYVAR(MSYMA)
                           LFREE = LWRK - KFREE + 1
                           IF (LFREE.LT.0)
     *                        CALL ERRWRK('QRSMO 2',KFREE-1,LWRK)
                           CALL X2INIT(1,MZYVAR(MSYMA),MZYVAR(MSYMC),
     *                              MSYMA,ISPINA,MSYMC,ISPINC,
     *                              WRK(KAVEC),VECC,WRK(KB2C),XINDX,
     *                              UDV,PV,
     *                              BSMLB(MSYMB,IBOP),MSYMB,ISPINB,
     *                              CMO,MJWOP,WRK(KFREE),LFREE)
                           IF (.NOT.X2TEST) THEN
                              KAVEC = KFREE
                              CALL REARSP(LURSP,KLEN,WRK(KAVEC),
     &                             LABEL,BLANK,FREQA,D0,MSYMA,0,
     &                             THCLR,FOUND,CONV,ANTSYM)
                              IF (.NOT. (FOUND .AND. CONV))
     &                             CALL RIOERR(FOUND,LABEL,FREQA,MSYMA)
                              IF (FREQA .LT. D0) THEN
                                 CALL DSWAP(KLEN/2,WRK(KAVEC),1,
     &                                      WRK(KAVEC+KLEN/2),1)
                                 IF (ANTSYM .LT. D0) CALL DSCAL(KLEN,
     &                                ANTSYM,WRK(KAVEC),1)
                              END IF
                           END IF
                           VAL = -FACMOM*DDOT(MZYVAR(MSYMA),
     *                                        WRK(KB2C),1,WRK(KAVEC),1)
                           SMOM = SMOM + VAL
                           IF( IPRRSP .GE. 5) THEN
                              WRITE(LUPRI,'(/A,2F20.8)')
     *                        ' E3+B2C         CONTRIBUTION TO SMOM'
     *                        ,VAL,SMOM
                           END IF
                           CALL FLSHFO(LUPRI)
C
C A2TEST
C
                           IF (RSPCI) THEN
                              KA2B = 1
                           ELSE
                              KA2B = 1 + MZYVAR(MSYMA)
                           END IF
                           KFREE = KA2B + MZYVAR(MSYMC)
                           LFREE = LWRK - KFREE + 1
                           IF (LFREE.LT.0)
     *                        CALL ERRWRK('QRSMO 3',KFREE-1,LWRK)
                           CALL A2INIT(1,MZYVAR(MSYMC),MZYVAR(MSYMB),
     *                              MSYMC,ISPINC,MSYMB,ISPINB,
     *                              VECC,VECB,WRK(KA2B),XINDX,
     *                              UDV,PV,
     *                              ASMLB(MSYMA,IAOP),MSYMA,ISPINA,
     *                              CMO,MJWOP,WRK(KFREE),LFREE)
                           VAL = - FACMOM *
     *                             DDOT(MZYVAR(MSYMC),WRK(KA2B),
     *                                  1,VECC,1)
                           SMOM = SMOM + VAL
                           IF( IPRRSP .GE. 5) THEN
                              WRITE(LUPRI,'(/A,2F20.8)')
     *                        ' E3+B2C+A2B     CONTRIBUTION TO SMOM',
     *                        VAL,SMOM
                           END IF
                           CALL FLSHFO(LUPRI)
                           IF (RSPCI) THEN
                              KA2C = 1
                           ELSE
                              KA2C = 1 + MZYVAR(MSYMA)
                           END IF
                           KFREE = KA2C + MZYVAR(MSYMB)
                           LFREE = LWRK - KFREE
                           IF (LFREE.LT.0)
     *                        CALL ERRWRK('QRSMO 4',KFREE-1,LWRK)
                           CALL A2INIT(1,MZYVAR(MSYMB),MZYVAR(MSYMC),
     *                              MSYMB,ISPINB,MSYMC,ISPINC,
     *                              VECB,VECC,WRK(KA2C),XINDX,
     *                              UDV,PV,
     *                              ASMLB(MSYMA,IAOP),MSYMA,ISPINA,
     *                              CMO,MJWOP,WRK(KFREE),LFREE)
                           VAL = - FACMOM*DDOT(MZYVAR(MSYMB),
     *                                WRK(KA2C),1,VECB,1)
                           SMOM = SMOM + VAL
                           IF( IPRRSP .GE. 5) THEN
                              WRITE(LUPRI,'(/A,2F20.8)')
     *                        ' E3+B2C+A2B+A2C CONTRIBUTION TO SMOM'
     *                        ,VAL,SMOM
                           END IF
                        WRITE(LUPRI,'(//A,2(/A,5X,A,2I5),/A,5X,I8,2I5)')
     *   ' Second order moment in a.u. for',
     *   ' A operator label,  symmetry, spin: ',
     *      ASMLB(MSYMA,IAOP),MSYMA,ISPINA,
     *   ' B operator label,  symmetry, spin: ',
     *      BSMLB(MSYMB,IBOP),MSYMB,ISPINB,
     *   ' Excited state no., symmetry, spin: ',IC,MSYMCX,ISPINC
                        WRITE(LUPRI,'(/A,3F12.6)')
     *   ' omega B, excitation energy, moment :',
     *      BSMFR(IBFR),EXCITA(MSYMC,IC,ISPINC+1),SMOM
                              IF (PHOSPH .AND.
     &                        ASMLB(MSYMA,IAOP)(2:8).EQ.'DIPLEN ' .AND.
     &                        BSMLB(MSYMB,IBOP)(2:8).EQ.' SPNORB' ) THEN
                                IDIP=INDEX('XYZ',ASMLB(MSYMA,IAOP)(1:1))
                                IHSO=INDEX('XYZ',BSMLB(MSYMB,IBOP)(1:1))
                                PHOSMAT(IDIP,IHSO,IC)=SMOM
                              END IF
                              IF (MNFPHO .AND.
     &                        ASMLB(MSYMA,IAOP)(2:8).EQ.'DIPLEN ' .AND.
     &                        BSMLB(MSYMB,IBOP)(2:8).EQ.'1MNF-SO' ) THEN
                                IDIP=INDEX('XYZ',ASMLB(MSYMA,IAOP)(1:1))
                                IHSO=INDEX('XYZ',BSMLB(MSYMB,IBOP)(1:1))
                                PHOSMN(IDIP,IHSO,IC)=SMOM
                               END IF
                           IF (TWOPHO) THEN
                              IF (IC.LE.IMXEXC) THEN
                                 IA=INDEX('XYZ',ASMLB(MSYMA,IAOP)(1:1))
                                 IB=INDEX('XYZ',BSMLB(MSYMB,IBOP)(1:1))
                                 RESFRE(IC,MSYMC)=
     &                                EXCITA(MSYMC,IC,ISPINC+1)
                                 RESSOM(IA,IB,IC,MSYMC)=SMOM
                                 RESSOM(IB,IA,IC,MSYMC)=SMOM
                              ELSE
                                 WRITE(LUPRI,'(/A)') 
     &             'WARNING!! not enough space in RESSOM result matrix'
                              END IF
                           END IF
                           CALL FLSHFO(LUPRI)

C------------------------------------------------------------------------*
C Calculation of MCD B term: the full residue SMOM*TM1 (Sonia)
C
                           IF (MCDCAL) THEN
                              LYESGP  = .FALSE.
                              A1 = ASMLB(MSYMA,IAOP)(1:1)
                              B1 = BSMLB(MSYMB,IBOP)(1:1)
                              IF (((A1.EQ.'X') .AND. (B1.EQ.'Y'))
     *                        .OR.((A1.EQ.'Y').AND.(B1.EQ.'X'))) THEN
                                 CLABEL = 'ZDIPLEN '
                                 LYESGP  = .TRUE.
                              END IF
                              IF (((A1.EQ.'X') .AND. (B1.EQ.'Z'))
     *                        .OR.((A1.EQ.'Z').AND.(B1.EQ.'X'))) THEN
                                 CLABEL = 'YDIPLEN '
                                 LYESGP  = .TRUE.
                              END IF
                              IF (((A1.EQ.'Y') .AND. (B1.EQ.'Z'))
     *                        .OR.((A1.EQ.'Z').AND.(B1.EQ.'Y'))) THEN
                                 CLABEL = 'XDIPLEN '
                                 LYESGP  = .TRUE.
                              END IF
                              IF (LYESGP) THEN
C
C Calculate propriety gradient. Returned as first element in WRK
C
                                 !backup variables in wrkrsp.h
                                 KSYMOPSV = KSYMOP
                                 KZYVARSV = KZYVAR
                                 KZVARSV  = KZVAR
                                 KZCONFSV = KZCONF
                                 KZYCONSV = KZYCON
                                 KZYWOPSV = KZYWOP
                                 KZWOPTSV = KZWOPT
                                 KSYMSTSV = KSYMST
                                 DO IOP = 1,KZWOPT
                                   JWOPSV(1,IOP) = JWOP(1,IOP)
                                   JWOPSV(2,IOP) = JWOP(2,IOP)
                                 ENDDO
C
                                 !set variables in wrkrsp.h to current value
                                 KSYMOP = MSYMCX
                                 KZYVAR = MZYVAR(MSYMC)
                                 KZVAR  = MZVAR(MSYMC)
                                 KZCONF = MZCONF(MSYMC)
                                 KZYCON = MZYCON(MSYMC)
                                 KZYWOP = MZYWOP(MSYMC)
                                 KZWOPT = MZWOPT(MSYMC)
                                 KSYMST = MSYMCX
                                 DO IOP = 1,KZWOPT
                                   JWOP(1,IOP) = MJWOP(1,IOP,KSYMOP)
                                   JWOP(2,IOP) = MJWOP(2,IOP,KSYMOP)
                                 ENDDO
c
                                 !GETGPV uses variables in wrkrsp.h
                                 CALL GETGPV(CLABEL,FC,FV,CMO,UDV,PV,
     &                                       XINDX,ANTSYM,
     &                                       WRK(KFREE),LFREE)

                                 TM1 = DDOT(KZYVAR,WRK(KFREE),1,VECC,1)
                                 BCONTR = SMOM*TM1*DP5
c
                      WRITE(LUPRI,'(/A,(/A,5X,A,2I5),/A,5X,I8,2I5)')
     *   ' First order moment in a.u. for',
     *   ' C operator label,  symmetry, spin: ',
     *                       CLABEL,MSYMCX,ISPINC,
     *   ' Excited state no., symmetry, spin: ',IC,MSYMCX,ISPINC
                           WRITE(LUPRI,'(A,2F12.6)')
     *   ' Excitation energy in au,    moment :            ',
     *                     EXCITA(MSYMC,IC,ISPINC+1),TM1
                           WRITE(LUPRI,'(/A,F12.6,/A)')
     *   ' B term contribution:          ', BCONTR,
     *   ' ------------------------------------------- '
                                 !restore variables in wrkrsp.h as originally
                                 KSYMOP = KSYMOPSV
                                 KZYVAR = KZYVARSV
                                 KZVAR  = KZVARSV
                                 KZCONF = KZCONFSV
                                 KZYCON = KZYCONSV
                                 KZYWOP = KZYWOPSV
                                 KZWOPT = KZWOPTSV
                                 KSYMST = KSYMSTSV
                                 DO IOP = 1,KZWOPT
                                   JWOP(1,IOP) = JWOPSV(1,IOP)
                                   JWOP(2,IOP) = JWOPSV(2,IOP)
                                 ENDDO
                              ELSE
                                 WRITE(LUPRI,*)'No MCD for this SOMOM'
                              ENDIF
                           END IF
                           CALL FLSHFO(LUPRI)
C ---------------------------------------------------------------endofmcd-------
 900                    CONTINUE
 800                 CONTINUE
 700              CONTINUE
 500           CONTINUE
            END IF
 400     CONTINUE
         IF (PHOSPH .OR. MNFPHO) THEN
            WRITE (LUPRI,'(/A/A,I3/A,I3/A/)')
     &         ' Phosphorescence electric dipole transition rates',
     &         ' from excited states of symmetry',MULD2H(MSYMC,IREFSY),
     &         ' to the ground state of symmetry',IREFSY,
     &         ' ----------------------------------'
            IF (PHOSPH) WRITE(LUPRI,'(A)')
     &        ' (H_SO: Full spin-orbit integrals used.)'
            IF (MNFPHO) WRITE(LUPRI,'(A)')
     &        ' (AMFI: Atomic Mean field spin-orbit Integrals used.)'
            IF (NSMCNV(MSYMC) .LE. 0) WRITE(LUPRI,'(/A)')
     &        ' No excited states calculated of this symmetry.'
            DO IC = 1,NSMCNV(MSYMC)
               WRITE(LUPRI,'(/A,I3/A/)')
     &  ' Phosphorescence transition rate from excited state no.',IC,
     &  ' (Triplet->singlet transition, high-temperature limit)'
               W=EXCITA(MSYMC,IC,ISPINC+1)
               RATE =D0
               RATEA=D0
               DO I=1,3
                  XDIP=D0
                  XDIPA=D0
                  DO J=1,3
                     IF (PHOSPH) XDIP =XDIP +PHOSMAT(I,J,IC)**2
                     IF (MNFPHO) XDIPA=XDIPA+ PHOSMN(I,J,IC)**2
                  END DO
                  IF (PHOSPH) THEN
                     XRATE=D4D3*(ALPHAC*W)**3*XDIP/XFSEC
                     WRITE(LUPRI,'(2A,1P,G17.5)')
     &               ' Partial rates (H_SO): ',
     &               CHAR(ICHAR('X')+I-1)//'-polarization', XRATE
                     RATE=RATE+XRATE
                   END IF
                   IF (MNFPHO) THEN
                     XRATEA=D4D3*(ALPHAC*W)**3*XDIPA/XFSEC
                     WRITE(LUPRI,'(2A,1P,G17.5)')
     &               ' Partial rates (AMFI): ',
     &               CHAR(ICHAR('X')+I-1)//'-polarization', XRATEA
                     RATEA=RATEA+XRATEA
                   END IF
               END DO
               IF (PHOSPH) WRITE(LUPRI,'(1P,2(/A,G12.5,A))')
     &        ' Total transition       rate    (H_SO)', RATE/D3,' s-1',
     &        ' Total phosphorescence lifetime (H_SO)', D3/RATE,' s'
               IF (MNFPHO) WRITE(LUPRI,'(1P,2(/A,G12.5,A))')
     &        ' Total transition       rate    (AMFI)', RATEA/D3,' s-1',
     &        ' Total phosphorescence lifetime (AMFI)', D3/RATEA,' s'
            END DO
         END IF
 300  CONTINUE
C
      IF (TWOPHO) THEN
         CALL TITLER('FINAL RESULTS FROM TWO-PHOTON CALCULATION',
     &        '*',112)
         WRITE(LUPRI,'(A64,5(/A64))') 
     & ' The two-photon absorption strength for an average molecular  ',
     & ' orientation is computed according to formulas given by       ',
     & ' P.R. Monson and W.M. McClain in J. Chem. Phys. 53:29, 1970   ',
     & ' and W.M. McClain in J. Chem. Phys. 55:2789, 1971.            ',
     & ' The absorption depends on the light polarization.            ',
     & ' A monochromatic light source is assumed.                     '
         WRITE(LUPRI,'(3(/A64))')
     & ' All results are presented in atomic units, except the        ',
     & ' excitation energy which is given in eV and two-photon cross  ',
     & ' section which is given in GM.                                '
         WRITE(LUPRI,'(3(/A))')
     &        '   Conversion factors:',
     &        '      1 a.u. = 1.896788 10^{-50} cm^4 s/photon',
     &        '      1 GM = 10^{-50} cm^4 s/photon'
         WRITE(LUPRI,'(/,3(/A52))') 
     &        '+--------------------------------+',
     &        '| Two-photon transition tensor S |',
     &        '+--------------------------------+'
         WRITE(LUPRI,'(A68)')
     &'---------------------------------------------------------------'     
         WRITE(LUPRI,'(A8,A4,A8,6A8,/A68)') 'Sym','No','Energy',
     &        'Sxx','Syy','Szz','Sxy','Sxz','Syz',
     &'---------------------------------------------------------------'
         DO ISYMC=1,NSYM
         DO IC = 1,NSMCNV(ISYMC)
            IF (IC.LE.IMXEXC) THEN
               WRITE(LUPRI,'(A4,2I4,F8.2,6F8.1)') ' ',ISYMC,IC,
     &              RESFRE(IC,ISYMC)*XTEV,
     &              (RESSOM(I,I,IC,ISYMC),I=1,3),
     &              RESSOM(1,2,IC,ISYMC),RESSOM(1,3,IC,ISYMC),
     &              RESSOM(2,3,IC,ISYMC)
            ELSE
               WRITE(LUPRI,'(/A)') 
     &            'WARNING!! NOT ENOUGH WORK MEMORY TO IN RESULT MATRIX'
            END IF
         END DO
         END DO
         WRITE(LUPRI,'(A68)')
     &'---------------------------------------------------------------'     
         WRITE(LUPRI,'(/,6(/A60))') 
     &        ' Transition probabilities (a.u.)         ',
     &        '-----------------------------------------',
     &        ' D  =  2*Df + 4*Dg, Linear   polarization',
     &        ' D  = -2*Df + 6*Dg, Circular polarization',
     &        ' Df = sum(i,j){ S_ii * S_jj }/30         ',
     &        ' Dg = sum(i,j){ S_ij * S_ij }/30         '
         WRITE(LUPRI,'(3(/A60))') 
     &        '          Two-photon cross sections                ',
     &        '---------------------------------------------------',
     &        ' sigma  =  8*pi^3*alpha^2*hbar/e^4 * E^2*D   (a.u.)'
         WRITE(LUPRI,'(3(/A53))') 
     &        '       Polarization ratio      ',
     &        '-------------------------------',
     &        '    R  = (-Df+3*Dg)/(Df+2*Dg)  '  
         WRITE(LUPRI,'(/,3(/A56))') 
     &        '+-----------------------------------+',
     &        '|   Two-photon absorption summary   |',
     &        '+-----------------------------------+'
         WRITE(LUPRI,'(A62,A)')
     &'-----------------------------------------------------------',
     &'----------------------'
         WRITE(LUPRI,'(A6,A4,A8,A14,4A11,A8/A62,A)') 'Sym','No',
     &        'Energy','Polarization','Df','Dg','D','sigma','R',
     &'-----------------------------------------------------------',
     &'----------------------'
         DO ISYMC=1,NSYM
         DO IC = 1,NSMCNV(ISYMC)
            IF (IC.LE.IMXEXC) THEN
               DF=0.0D0
               DG=0.0D0
               DO I=1,3
               DO J=1,3
                  DF=DF+RESSOM(I,I,IC,ISYMC)*RESSOM(J,J,IC,ISYMC)
                  DG=DG+RESSOM(I,J,IC,ISYMC)**2
               END DO
               END DO
               DF = DF/D30
               DG = DG/D30
               R = (-DF+3*DG)/(DF+2*DG)
               D = 2*DF+4*DG
               s = D * RESFRE(IC,ISYMC)**2 * AU_TO_GM
               WRITE(LUPRI,'(A2,2I4,F8.2,A14,4E11.3,F8.2)')' ',ISYMC,IC,
     &              RESFRE(IC,ISYMC)*XTEV,'Linear     ',DF,DG,D,s,R
               D=-2*DF+6*DG
               s = D * RESFRE(IC,ISYMC)**2 * AU_TO_GM
               WRITE(LUPRI,'(A2,2I4,F8.2,A14,4E11.3,F8.2)')' ',ISYMC,IC,
     &              RESFRE(IC,ISYMC)*XTEV,'Circular   ',DF,DG,D,s,R
               IF (IPRRSP.GT.2) THEN
                  WRITE(LUPRI,'(/,6x,a)')
     &    ' Rotationally averaged two-photon transition strengths',
     &    '                 and rate constants'
                  DELTA_TPlpa = 0.0D0
                  DELTA_TPlpa = ((2*DF)+(4*DG))
                  DELTA_TPlpe = 0.0D0
                  DELTA_TPlpe = ((-1*DF)+(5*DG))
                  DELTA_TPc = 0.0D0
                  DELTA_TPc = ((-2*DF)+(6*DG))
                  AUTIME = HBAR/XTJ
                  dK_1 = 0.0D0
                  dK_1 = (XTANG*XTANG*XTANG*XTANG)*AUTIME
                  dK_2 = 0.0D0
                  dK_2 = (8*PI*PI*PI)*(ALPHAC*ALPHAC)
                  dK_3 = 0.0D0
                  dK_3 = RESFRE(IC,ISYMC)*RESFRE(IC,ISYMC)
                  dK_lpa = 0.0D0
                  dK_lpa = dK_1*dK_2*dK_3*DELTA_TPlpa
                  dK_lpe = 0.0D0
                  dK_lpe = dK_1*dK_2*dK_3*DELTA_TPlpe
                  dK_c = 0.0D0
                  dK_c = dK_1*dK_2*dK_3*DELTA_TPc
                  WRITE(LUPRI,'(/,6x,a)')
     &    ' +-------------------+----------------+----------------+'
                  WRITE(LUPRI,'(6x,a)')
     &    ' |   Polarization    |    DELTA_TP    |   K (0 -> f)   |'
                  WRITE(LUPRI,'(6x,a)')
     &    ' +-------------------+----------------+----------------+'
                  WRITE(LUPRI,'(6x,a,f13.8,a,E14.6,a)')
     &    ' |  linear (para)    |', DELTA_TPlpa, '   |', dK_lpa,'  |'
                  WRITE(LUPRI,'(6x,a)')
     &    ' +-------------------+----------------+----------------+'
                  WRITE(LUPRI,'(6x,a,f13.8,a,E14.6,a)')
     &    ' |  linear (perp)    |', DELTA_TPlpe, '   |', dK_lpe,'  |'
                  WRITE(LUPRI,'(6x,a)')
     &    ' +-------------------+----------------+----------------+'
                  WRITE(LUPRI,'(6x,a,f13.8,a,E14.6,a)')
     &    ' |    circular       |', DELTA_TPc,   '   |', dK_c,'  |'
                  WRITE(LUPRI,'(6x,a)')
     &    ' +-------------------+----------------+----------------+'
               END IF
*********************************************************************
            ELSE
               WRITE(LUPRI,'(/A)') 
     &            'WARNING!! NOT ENOUGH WORK MEMORY TO IN RESULT MATRIX'
            END IF
         END DO
         END DO
         WRITE(LUPRI,'(A62,A)')
     &'-----------------------------------------------------------',
     &'----------------------'
      END IF
C
C
      CALL QEXIT('QRSMO')
      RETURN
      END
C  /* Deck qrexm */
      SUBROUTINE QREXM(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *                  XINDX,VECB,VECC,MJWOP,WRK,LWRK)
C
C PURPOSE:
C CALCULATION OF TRANSITION MOMENTS BETWEEN EXCITED STATES
C < B | A | C>
C
#include <implicit.h>
C
      LOGICAL FOUND, CONV
      CHARACTER*8 BLANK, LABEL
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),VECB(*),VECC(*),WRK(*)
C
      PARAMETER ( D0 = 0.0D0, BLANK = '        ' )
#include <iratdef.h>
C
#include <priunit.h>
#include <rspprp.h>
#include <infrsp.h>
#include <infpp.h>
#include <infsmo.h>
#include <indqr.h>
#include <inforb.h>
#include <maxorb.h>
#include <infvar.h>
#include <qrinf.h>
      DIMENSION MJWOP(2,MAXWOP,8)
#include <infpri.h>
#include <infspi.h>
#include <inftap.h>
#include <inflr.h>
C
CTODO? HJ-910405: exchange 400 and 900 loop
C In 900 loop, if NGPPP(MSYMC).GT.0 calculate <0/C/i> moments
C Calculate oscillator strengths,
C    if DIPVEL then divide by frequency
C
C panor:
C     DOESA=.TRUE.
C     For excited state absorption we require the B-state to be the
C     first state of symmetry ESASYM and that the C-state is separate
C     from the B-state.
C
      CALL QENTER('QREXM')
      IPRRSP = IPREXM
      DO 300 MSYMC = 1,NSYM
         MSYMCX = MULD2H(IREFSY,MSYMC)
         IF (DOESA .OR. EXMTES) THEN
            MSYMBE = NSYM
         ELSE
            MSYMBE = MSYMC
         END IF
         DO 400 MSYMB = 1,MSYMBE
            IF (DOESA .AND. MSYMB.NE.ESASYM) GOTO 400
            MSYMBX = MULD2H(IREFSY,MSYMB)
            MSYMA = MULD2H(MSYMC,MSYMB)
            IF ( (NPPCNV(MSYMC) .GT. 0) .AND. (NPPCNV(MSYMB) .GT. 0)
     *          .AND. (NGPPP(MSYMA).GT.0) ) THEN
               DO 900 IC = 1,NPPCNV(MSYMC)
                  IF (DOESA .AND. MSYMC.EQ.MSYMB .AND. IC.EQ.1) GOTO 900
                  IF(.NOT.LCMEXC(IC,MSYMC)) GO TO 900
                  IF (ISPINC.EQ.ISPINB) THEN
                     ICEXNU = INQREX(MSYMC,IC)
                  ELSE
                     IF ( ISPINC.EQ.0 ) THEN
                        ICEXNU = INQREX(MSYMC,IC)
                     ELSE
                        ICEXNU = NEXLBL + INQREX(MSYMC,IC)
                     END IF
                  END IF
                  FREQC = EXCITA(MSYMC,IC,ISPINC+1)
                  CALL REARSP(LURSP,KLEN,VECC,
     &                 'EXCITLAB',BLANK,FREQC,D0,MSYMC,IC,THCPP,
     &                 FOUND,CONV,ANTSYM)
                  IF (.NOT. (FOUND .AND. CONV)) THEN
                     IF (.NOT. FOUND) THEN
                        WRITE (LUPRI,'(/3A,F7.5,A,I3,/A)') 
     &                       ' Response label ','EXCITLAB',
     &                       ' with frequency ',FREQC,
     &                       ' and symmetry',MSYMC,
     &                       ' not found on file RSPVEC'
                        CALL QUIT('Response vector not '//
     &                       'found on file')
                     ELSE
                        WRITE (LUPRI,'(/3A,F7.5,/A,I3,A)') 
     &                       ' @WARNING>>>>'//
     &                       ' Response label ','EXCITLAB',
     &                       ' with frequency ',FREQC,
     &                       ' and symmetry',MSYMC,
     &                       ' not converged on file RSPVEC'
                     END IF
                  END IF
                  IF (FREQC .LT. D0) THEN
                     CALL DSWAP(KLEN/2,VECC,1,VECC(1+KLEN/2),1)
                     IF (ANTSYM .LT. D0) CALL DSCAL(KLEN,ANTSYM,VECC,1)
                  END IF
                  IF (IPRRSP.GT.10) THEN
                     WRITE(LUPRI,'(/A)')' --- C OPERATOR --- '
                     IF ( ISPINC.EQ.0 ) THEN
                        WRITE(LUPRI,'(/A)')' SINGLET VECTOR '
                     ELSE
                        WRITE(LUPRI,'(/A)')' TRIPLET VECTOR '
                     END IF
                     WRITE(LUPRI,'(/A,I5,3X,A,I5,/A,2X,D12.6)')
     *               ' QREXM: EIGENVECTOR ',ICEXNU,
     *               ' SYMMETRY',MSYMC,' EXCITATION ENERGY'
     *               ,FREQC
                     CALL RSPPRC(VECC,MZCONF(MSYMC),
     *                                  MZVAR(MSYMC),LUPRI)
                     CALL RSPPRW(VECC(1+MZCONF(MSYMC)),MJWOP,
     *                                  MZWOPT(MSYMC),MSYMC,
     *                                  MZVAR(MSYMC),LUPRI)
                  END IF
                  IF (DOESA) THEN
                     IBE = 1
                  ELSE IF (EXMTES .OR. MSYMB.NE.MSYMC) THEN
                     IBE = NPPCNV(MSYMB)
                  ELSE
                     IBE = IC
                  END IF
                  DO 1000 IB = 1,IBE
                     IF(.NOT.LCMEXC(IB,MSYMB)) GO TO 1000
                     IF (ISPINC.EQ.ISPINB) THEN
                        IBEXNU = INQREX(MSYMB,IB)
                     ELSE
                        IF ( ISPINB.EQ.0 ) THEN
                           IBEXNU = INQREX(MSYMB,IB)
                        ELSE
                           IBEXNU = NEXLBL + INQREX(MSYMB,IB)
                        END IF
                     END IF
                     FREQB = EXCITA(MSYMB,IB,ISPINB+1)
                     CALL REARSP(LURSP,KLEN,VECB,
     &                    'EXCITLAB',BLANK,FREQB,D0,MSYMB,IB,THCPP,
     &                    FOUND,CONV,ANTSYM)
                     IF (.NOT. (FOUND .AND. CONV)) THEN
                        IF (.NOT. FOUND) THEN
                           WRITE (LUPRI,'(/3A,F7.5,A,I3,/A)') 
     &                          ' Response label ','EXCITLAB',
     &                          ' with frequency ',FREQB,
     &                          ' and symmetry',MSYMB,
     &                          ' not found on file RSPVEC'
                           CALL QUIT('Response vector not '//
     &                          'found on file')
                        ELSE
                           WRITE (LUPRI,'(/3A,F7.5,/A,I3,A)') 
     &                          ' @WARNING>>>>'//
     &                          ' Response label ','EXCITLAB',
     &                          ' with frequency ',FREQB,
     &                          ' and symmetry',MSYMB,
     &                          ' not converged on file RSPVEC'
                        END IF
                     END IF
                     IF (FREQB .LT. D0) THEN
                        CALL DSWAP(KLEN/2,VECB,1,VECB(1+KLEN/2),1)
                        IF (ANTSYM .LT. D0) CALL DSCAL(KLEN,ANTSYM,
     &                                                 VECB,1)
                     END IF
                     IF (IPRRSP.GT.10) THEN
                        WRITE(LUPRI,'(/A)')' --- B OPERATOR --- '
                        IF ( ISPINB.EQ.0 ) THEN
                           WRITE(LUPRI,'(/A)')' SINGLET VECTOR'
                        ELSE
                           WRITE(LUPRI,'(/A)')' TRIPLET VECTOR'
                        END IF
                        WRITE(LUPRI,'(/A,I5,3X,A,I5,
     *                  /A,D14.6)')
     *                  ' QREXM: EIGENVECTOR ',IBEXNU,
     *                  ' SYMMETRY',MSYMB,' EXCITATION ENERGY'
     *                  ,FREQA
                        CALL RSPPRC(VECB,MZCONF(MSYMB),
     *                                  MZVAR(MSYMB),LUPRI)
                        CALL RSPPRW(VECB(1+MZCONF(MSYMB)),MJWOP,
     *                                  MZWOPT(MSYMB),MSYMB,
     *                                  MZVAR(MSYMB),LUPRI)
                     END IF
                     CALL DSWAP(MZVAR(MSYMB),VECB,1,
     *                          VECB(1+MZVAR(MSYMB)),1)
                     CALL FLSHFO(LUPRI)
                     IF(.NOT.RSPCI) THEN
                        CEXMBE = EXCITA(MSYMC,IC,ISPINC+1)
     *                         - EXCITA(MSYMB,IB,ISPINB+1)
                        IF (E3TEST) THEN
                           IAOP   = 1
                           KAVEC  = 1
                           KX3WRK = KAVEC + MZYVAR(MSYMA)
                           LX3WRK = LWRK  - KX3WRK + 1
                           IF (LX3WRK.LT.0)
     *                        CALL ERRWRK('QREXM 1',KX3WRK-1,LWRK)
                           IF ( ISPINA.EQ.0 ) THEN
                              IAFRNU = NEXLB2
     *                          + INQRLR(LBLPP(MSYMA,IAOP),CEXMBE,MSYMA)
                           ELSE
                              IAFRNU = NEXLB2 + NLRLBL
     *                          + INQRTR(LBLPP(MSYMA,IAOP),CEXMBE,MSYMA)
                           END IF
                        ELSE
                           KAVEC  = 1
                           KX3WRK = KAVEC
                           LX3WRK = LWRK
                        END IF
                        M1 = 1
                        CALL T3DRV(M1,MSYMA,MSYMB,MSYMC,VECB,VECC,
     *                      E3TEST,WRK(KAVEC),
     *                      EXCITA(MSYMB,IB,ISPINB+1),
     *                     -EXCITA(MSYMC,IC,ISPINC+1),
     *                      XINDX,UDV,PV,MJWOP,
     *                      WRK(KX3WRK),LX3WRK,CMO,FC,FV)
                        IF (E3TEST) THEN
                           CALL DCOPY(MZYVAR(MSYMA),WRK(KX3WRK),1,WRK,1)
                        END IF
                        CALL FLSHFO(LUPRI)
                     END IF
                     DO 1100 IAOP = 1,NGPPP(MSYMA)
                        IF (IPRRSP .GE. 5)
     *                  WRITE(LUPRI,'(/A,I5,A,I5,A,I5,/A,I5,A,I5,A,I5,
     *                     /A,I5,A,I5,/3A)')
     *                   ' SYMMETRY OF OPERATOR C ',MSYMC,
     *                   ' SPIN OF OPERATOR C ',ISPINC,
     *                   ' STATE NUMBER',IC,
     *                   ' SYMMETRY OF OPERATOR B ',MSYMB,
     *                   ' SPIN OF OPERATOR B ',ISPINB,
     *                   ' STATE NUMBER',IB,
     *                   ' SYMMETRY OF OPERATOR A ',MSYMA,
     *                   ' SPIN OF OPERATOR A ',ISPINA,
     *                   ' LABEL OF OPERATOR A: ',LBLPP(MSYMA,IAOP),
     *                   ' FOR WHICH TRANSITION MOMENT IS CALCULATED'
                        IF (RSPCI) THEN
                           XMOM   = D0
                        ELSE
                           IF ( ISPINA.EQ.0 ) THEN
                              IAFRNU = NEXLB2
     *                          + INQRLR(LBLPP(MSYMA,IAOP),CEXMBE,MSYMA)
                           LABEL = QRLBL(IAFRNU-NEXLB2)
                           FREQA  = QRFREQ(IAFRNU-NEXLB2)
                           ELSE
                              IAFRNU = NEXLB2 + NLRLBL
     *                          + INQRTR(LBLPP(MSYMA,IAOP),CEXMBE,MSYMA)
                           LABEL = TRLBL(IAFRNU-NEXLB2)
                           FREQA  = TRFREQ(IAFRNU-NEXLB2)
                           END IF
                           CALL REARSP(LURSP,KLEN,WRK(1+MZYVAR(MSYMA)),
     &                          LABEL,BLANK,FREQA,D0,MSYMA,0,THCLR,
     &                          FOUND,CONV,ANTSYM)
                           IF (.NOT. (FOUND .AND. CONV))
     &                          CALL RIOERR(FOUND,LABEL,FREQA,MSYMA)
                           IF (FREQA .LT. D0) THEN
                              CALL DSWAP(KLEN/2,WRK(1+MZYVAR(MSYMA)),1,
     &                             WRK(1+MZYVAR(MSYMA)+KLEN/2),1)
                              IF (ANTSYM .LT. D0) CALL DSCAL(KLEN,
     &                             ANTSYM,WRK(1+MZYVAR(MSYMA)),1)
                           END IF
                           IF (IPRRSP.GT.5) THEN
                             WRITE(LUPRI,'(/A)')' --- A OPERATOR ---'
                             IF ( ISPINA.EQ.0) THEN
                                WRITE(LUPRI,'(/A)')
     *                          ' SINGLET SOLUTION VECTOR'
                             ELSE
                                WRITE(LUPRI,'(/A)')
     *                          ' TRIPLET SOLUTION VECTOR'
                             END IF
                             WRITE(LUPRI,'(/A,3X,A,3X,A,2X,
     *                                 D12.6,/A,I5)')
     *                       ' QREXM: SOLUTION VECTOR  LABEL',LABEL,
     *                       ' FREQUENCY ',FREQA,' SYMMETRY',MSYMA
                             WRITE (LUPRI,'(/A)') ' EIGENVECTOR'
                             CALL RSPPRC(WRK(MZYVAR(MSYMA)+1),
     *                                MZCONF(MSYMA),
     *                                MZVAR(MSYMA),LUPRI)
                             IWRKOF = MZYVAR(MSYMA)+MZCONF(MSYMA)+1
                             CALL RSPPRW(WRK(IWRKOF),MJWOP,
     *                                MZWOPT(MSYMA),MSYMA,MZVAR(MSYMA),
     *                                LUPRI)
                           END IF
                           XMOM = DDOT(MZYVAR(MSYMA),WRK(1),1,
     *                                 WRK(1+MZYVAR(MSYMA)),1)
                           IF( IPRRSP .GE. 5) THEN
                              WRITE(LUPRI,'(/A,2F20.8)')
     *                     ' E3         CONTRIBUTION TO XMOM',XMOM,XMOM
                           END IF
                        END IF
                        CALL FLSHFO(LUPRI)
C
C A2TEST
C
                        IF (RSPCI) THEN
                           KA2B = 1
                        ELSE
                           KA2B = 1 + MZYVAR(MSYMA)
                        END IF
                        KFREE = KA2B + MZYVAR(MSYMC)
                        LFREE = LWRK - KFREE + 1
                        IF (LFREE.LT.0)
     *                     CALL ERRWRK('QREXM 2',KFREE-1,LWRK)
                        CALL A2INIT(1,MZYVAR(MSYMC),MZYVAR(MSYMB),
     *                              MSYMC,ISPINC,MSYMB,ISPINB,
     *                              VECC,VECB,WRK(KA2B),XINDX,
     *                              UDV,PV,
     *                              LBLPP(MSYMA,IAOP),MSYMA,ISPINA,
     *                              CMO,MJWOP,WRK(KFREE),LFREE)
                        VAL = - DDOT(MZYVAR(MSYMC),WRK(KA2B),1,VECC,1)
                        XMOM = XMOM + VAL
                        IF( IPRRSP .GE. 5) THEN
                           WRITE(LUPRI,'(/A,2F20.8)')
     *                     ' E3+A2B     CONTRIBUTION TO XMOM', VAL,XMOM
                        END IF
                        CALL FLSHFO(LUPRI)
                        IF (RSPCI) THEN
                           KA2C = 1
                        ELSE
                           KA2C = 1 + MZYVAR(MSYMA)
                        END IF
                        KFREE = KA2C + MZYVAR(MSYMB)
                        LFREE = LWRK - KFREE + 1
                        IF (LFREE.LT.0)
     *                     CALL ERRWRK('QREXM 3',KFREE-1,LWRK)
                        CALL A2INIT(1,MZYVAR(MSYMB),MZYVAR(MSYMC),
     *                              MSYMB,ISPINB,MSYMC,ISPINC,
     *                              VECB,VECC,WRK(KA2C),XINDX,
     *                              UDV,PV,
     *                              LBLPP(MSYMA,IAOP),MSYMA,ISPINA,
     *                              CMO,MJWOP,WRK(KFREE),LFREE)
                        VAL =  - DDOT(MZYVAR(MSYMB),WRK(KA2C),
     *                               1,VECB,1)
                        XMOM = XMOM + VAL
                        IF( IPRRSP .GE. 5) THEN
                           WRITE(LUPRI,'(/A,2F20.8)')
     *                     ' E3+A2B+A2C CONTRIBUTION TO XMOM', VAL,XMOM
                        END IF
                        WRITE(LUPRI,'(//A/A,5X,A,2I5,2(/A,5X,I8,2I5))')
     *   '@Transition moment <B | A | C> in a.u. for',
     *   '@A operator label,    symmetry, spin: ',
     *      LBLPP(MSYMA,IAOP),MSYMA,ISPINA,
     *   '@B excited state no., symmetry, spin: ',IB,MSYMBX,ISPINB,
     *   '@C excited state no., symmetry, spin: ',IC,MSYMCX,ISPINC
                           WRITE(LUPRI,'(/A,3F13.8)')
     *   '@B and C excitation energies, moment:',
     *   EXCITA(MSYMB,IB,ISPINB+1),EXCITA(MSYMC,IC,ISPINC+1),XMOM
                        CALL FLSHFO(LUPRI)
 1100                CONTINUE
 1000             CONTINUE
 900           CONTINUE
            END IF
 400     CONTINUE
 300  CONTINUE
      CALL QEXIT('QREXM')
      RETURN
      END
C  /* Deck setzy */
      SUBROUTINE SETZY(MJWOP)
C
#include <implicit.h>
C
#include <wrkrsp.h>
#include <maxorb.h>
#include <infvar.h>
#include <qrinf.h>
      DIMENSION MJWOP(2,MAXWOP,8)
C
C SET VARIABLES FOR SYMMETRY KSYMOP
C
      MZVAR(KSYMOP) = KZVAR
      MZYVAR(KSYMOP)= KZYVAR
      MZCONF(KSYMOP)= KZCONF
      MZYCON(KSYMOP)= KZYCON
      MZWOPT(KSYMOP)= KZWOPT
      MZYWOP(KSYMOP)= KZYWOP
      DO 100 IOP = 1,KZWOPT
         MJWOP(1,IOP,KSYMOP) = JWOP(1,IOP)
         MJWOP(2,IOP,KSYMOP) = JWOP(2,IOP)
 100  CONTINUE
      RETURN
      END
C  /* Deck qrtrve */
      SUBROUTINE QRTRVE(CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,
     &                  WRK,LWRK)
C
#include <implicit.h>
#include <iratdef.h>
#include <dummy.h>
C
      CHARACTER*8 BLANK
      LOGICAL FOUND, CONV
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
      PARAMETER ( THRSML = 1.0D-9, D0 = 0.0D0, BLANK = '        ' )
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
#include <inforb.h>
#include <inftap.h>
C
C Local variables
C
      CHARACTER*8 LRLAB(MXLRQR)
C
C DETERMINE SECOND ORDER MOLECULAR PROPERTIES
C
      CALL QENTER('QRTRVE')
      IF (NTRQR(KSYMOP) .LE. 0) GO TO 9999
C
C     count number of different property labels
C     and max number of frequencies
      NLRLAB = 0
      MXFRQ  = 0
      IF (IPRLR .GT. 20) THEN
         WRITE (LUPRI,*)
     *   'Test output from QRTRVE of all labels and frequencies'
      END IF
      DO 220 IOP = 1,NTRQR(KSYMOP)
         IOPVEC = ITRQR(KSYMOP) + IOP
         IF (IPRLR .GT. 20) THEN
            WRITE (LUPRI,*) IOP, IOPVEC, TRLBL(IOPVEC),TRFREQ(IOPVEC)
         END IF
         DO 200 ILRLAB = 1,NLRLAB
            IF (TRLBL(IOPVEC) .EQ. LRLAB(ILRLAB)) GO TO 220
  200    CONTINUE
         NLRLAB = NLRLAB + 1
         LRLAB(NLRLAB) = TRLBL(IOPVEC)
C
         NFRQ = 1
         DO 210 JOP = IOP+1,NTRQR(KSYMOP)
            JOPVEC = ITRQR(KSYMOP) + JOP
            IF (TRLBL(JOPVEC) .EQ. LRLAB(ILRLAB)) NFRQ = NFRQ + 1
  210    CONTINUE
         IF (IPRLR .GT. 20) WRITE (LUPRI,*) 'New label, NFREQ =',NFRQ
         MXFRQ = MAX(MXFRQ,NFRQ)
  220 CONTINUE
C     allocate core
      LNEEDA = 2*MAXRM*MAXRM + MAXRM + 4*KZYVAR
C     ... 3*KZYVAR is an estimate of space needed in RSPCTL
      LNEEDB = 1 + 2*MAXRM
      NFRQMX = (LWRK - LNEEDA) / LNEEDB
      IF (NFRQMX .LT. MXFRQ) THEN
         WRITE (LUPRI,9100) LWRK,LNEEDA+LNEEDB*MXFRQ
         CALL QTRACE(LUPRI)
         CALL QUIT('QRTRVE: INSUFFICIENT WORK MEMORY TO '//
     &             'SOLVE LINEAR EQUATIONS ')
      ENDIF
 9100 FORMAT(/' QRTRVE, work space too small for 3 (z,y)-vectors',
     *       /'         had:',I10,', need more than:',I10)
C
      KREDE  = 1
      KREDS  = KREDE  + MAXRM*MAXRM
      KIBTYP = KREDS  + MAXRM*MAXRM
      KEIVAL = KIBTYP + MAXRM
      KRESID = KEIVAL + MXFRQ
      KEIVEC = KRESID + MXFRQ
      KREDGP = KEIVEC + MAXRM*MXFRQ
      KGP    = KREDGP + MAXRM*MXFRQ
      KWRK1  = KGP    + KZYVAR
      LWRK1  = LWRK   - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('QRTRVE',KWRK1-1,LWRK)
C     reuse KGP space for solution vectors.
      KBVECS = KGP + KZYVAR
      NFRQMX = (LWRK - KBVECS - KZYVAR) / KZYVAR
C     ... allocated KZYVAR for KWRKE
      NFRQMX = MIN(NFRQMX,MXFRQ)
      KWRKE  = KBVECS + NFRQMX*KZYVAR
C
      THCRSP = THCLR
      IPRRSP = IPRLR
      MAXIT  = MAXITL
C
C     Call RSPCTL to solve linear set of response equations
C
      DO 680 ILRLAB = 1,NLRLAB
         NFREQ  = 0
         DO 410 IOP = 1,NTRQR(KSYMOP)
            IOPVEC = ITRQR(KSYMOP) + IOP
            IF (TRLBL(IOPVEC) .EQ. LRLAB(ILRLAB)) THEN
               CALL REARSP(LURSP,KLEN,WRK(KBVECS),TRLBL(IOPVEC),BLANK,
     &                     TRFREQ(IOPVEC),D0,KSYMOP,0,THCLR,FOUND,
     &                     CONV,ANTSYM)
               IF (FOUND .AND. CONV) THEN
                  WRITE(LUPRI,'(/A,/A12,A10,/A42)')
     *                 ' Converged solution vector already on file',
     *                 TRLBL(IOPVEC),'freq',
     *                 ' -----------------------------------------'
                  WRITE(LUPRI,'(F22.6)') TRFREQ(IOPVEC)
               ELSE
                  WRK(KEIVAL+NFREQ) = TRFREQ(IOPVEC)
                  NFREQ = NFREQ + 1
               END IF
            END IF
  410    CONTINUE
         IF (NFREQ .EQ. 0) GOTO 680
C
         CALL GETGPV(LRLAB(ILRLAB),FC,FV,CMO,UDV,PV,XINDX,
     &               ANTSYM,WRK(KGP),LWRK1)
         GPNORM = DNRM2(KZYVAR,WRK(KGP),1)
C
         WRITE (LUPRI,'(//A,I3,/2A,/A,(T25,5F10.6))')
     &   ' QRTRVE -- linear response calculation for symmetry',KSYMOP,
     &   ' QRTRVE -- operator label : ',LRLAB(ILRLAB),
     &   ' QRTRVE -- frequencies :',(WRK(KEIVAL+I),I=0,NFREQ-1)
         IF (GPNORM .LT. THRSML) THEN
            WRITE (LUPRI,*) ' >>> RSPCTL skipped because norm of'
            WRITE (LUPRI,*) '     property vector is only',GPNORM
            CALL DZERO(WRK(KBVECS),KZYVAR)
            DO 450 IFREQ = 1,NFREQ
               CALL WRTRSP(LURSP,KZYVAR,WRK(KBVECS),LRLAB(ILRLAB),BLANK,
     &                     WRK(KEIVAL + IFREQ - 1),D0,KSYMOP,0,D0,
     &                     ANTSYM)
 450        CONTINUE
            GO TO 680
         END IF
C
         KZRED  = 0
         KZYRED = 0
         KEXSIM = NFREQ
         KEXCNV = KEXSIM
         CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *               .TRUE.,LRLAB(ILRLAB),BLANK,WRK(KGP),
     *               WRK(KREDGP),WRK(KREDE),WRK(KREDS),
     *               WRK(KIBTYP),WRK(KEIVAL),WRK(KRESID),WRK(KEIVEC),
     *               XINDX,WRK(KWRK1),LWRK1)
C        CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
C    *               LINEQ,GP,REDGP,REDE,REDS,
C    *               IBTYP,EIVAL,EIVEC,XINDX,WRK,LWRK)
C
         DO 580 IFREQ = 1,NFREQ,NFRQMX
            NBX   = MIN(NFRQMX,(NFREQ+1-IFREQ))
            IBOFF = IFREQ - 1
            CALL RSPEVE(WRK(KIBTYP),WRK(KEIVAL),WRK(KEIVEC),
     *                  WRK(KBVECS),WRK(KWRKE),NBX,IBOFF)
C           CALL RSPEVE(IBTYP,EIVAL,EIVEC,BVECS,WRK,NBX,IBOFF)
C
            JBVECS = KBVECS
            DO 560 JFREQ = IFREQ,IFREQ-1+NBX
               IF (IPRRSP.GE.3) THEN
                  WRITE(LUPRI,'(/A,3X,A,3X,A,1P,D14.6,/A,I5)')
     *            '@QRTRVE: TRIPLET SOLUTION   LABEL',LRLAB(ILRLAB),
     *            ' FREQUENCY ',WRK(KEIVAL-1+JFREQ),
     *            ' SYMMETRY',KSYMOP
                  JEIVEC = KEIVEC + (JFREQ-1)*KZYRED
                  VAL = DDOT(KZYRED,WRK(KREDGP),1,WRK(JEIVEC),1)
                  WRITE(LUPRI,'(/A,1P,G20.12)')
     *            '@QRTRVE: Value of linear response property:',VAL
               END IF
               IF (IPRRSP.GT.100) THEN
                  WRITE (LUPRI,*) 'Column 1 = Z, Column 2 = Y'
                  CALL OUTPUT(WRK(JBVECS),1,KZVAR,1,2,KZVAR,2,1,LUPRI)
               END IF
               CALL WRTRSP(LURSP,KZYVAR,WRK(JBVECS),LRLAB(ILRLAB),BLANK,
     &                     WRK(KEIVAL - 1 + JFREQ),D0,KSYMOP,0,
     &                     WRK(KRESID - 1 + JFREQ),ANTSYM)
               JBVECS = JBVECS + KZYVAR
  560       CONTINUE
  580    CONTINUE
  680 CONTINUE
C
C *** end of QRTRVE --
C
 9999 CALL QEXIT('QRTRVE')
      RETURN
      END
      SUBROUTINE BCCHK(DOCAL,LURSPRES,ISYMA,ISYMB,ISYMC,FREQB,
     &                 FREQC,BLAB,CLAB,SPNFC1,SPNFC2,SPNFC,SPNSD1,
     *                 SPNSD2,SPNSD,SPNSU1,SPNSU2,SPNSU)
C
#include <implicit.h>
#include <iratdef.h>
C
C PURPOSE:
C Check if this is a unique beta calculation or equal to an already
C calculated due to permutation symmetry
C
C
      PARAMETER (THD = 1.0D-8)
      PARAMETER (MXCALC = 1000)
C 
      LOGICAL DOCAL, MEMFLG, DIPLEN, ONFIL
      CHARACTER*8 LAB(MXCALC,3), ALAB, BLAB, CLAB
      CHARACTER DIR(3)*1, STRING*75, LABEL(3)*8
      DIMENSION FRQ(MXCALC,3), ISYM(MXCALC,3)
      DIMENSION SPNFC1(*), SPNFC2(*), SPNFC(*), SPNSD1(*),
     *          SPNSD2(*), SPNSD(*), SPNSU1(*), SPNSU2(*), SPNSU(*)
#include <maxorb.h>
#include <priunit.h>
#include <infvar.h>
#include <infrsp.h>
#include <qrinf.h>
#include <rspprp.h>
#include <infpp.h>
#include <infhyp.h>
#include <infpri.h>
#include <inforb.h>
#include <infsmo.h>
#include <infspi.h>
C
      SAVE LAB, FRQ, N, ISYM
      DATA N/0/
C
      FREQA = -FREQB-FREQC
C
C MEMFLG is false only if all A-operators are done before,
C DOCAL controls a certain A-operator in the  IAOP loop
C
      ONFIL  = .FALSE.
      MEMFLG = .FALSE.
      DO IAOP = 1,NAQROP(ISYMA)
         DOCAL = .TRUE.
         ALAB = AQRLB(ISYMA,IAOP)
C
         IF (ALAB(2:7).EQ.'DIPLEN'
     *        .AND. BLAB(2:7).EQ.'DIPLEN'
     *        .AND. CLAB(2:7).EQ.'DIPLEN') THEN
            DIPLEN = .TRUE.
         ELSE
            DIPLEN = .FALSE.
         END IF
C
Clf         goto 1234
         DO I = 1,N
         DO J = 1,3
         DO K = 1,3
         IF (K.NE.J) THEN
            DO L = 1,3
            IF (L.NE.K .AND. L.NE.J) THEN
C
               IF ( ALAB.EQ.LAB(I,J) .AND.
     *              BLAB.EQ.LAB(I,K) .AND.
     *              CLAB.EQ.LAB(I,L) .AND.
     *              ISYMA.EQ.ISYM(I,J) .AND.
     *              ISYMB.EQ.ISYM(I,K) .AND.
     *              ISYMC.EQ.ISYM(I,L) .AND.
     *              ABS( FREQA-FRQ(I,J)).LT.THD .AND.
     *              ABS( FREQB-FRQ(I,K)).LT.THD .AND.
     *              ABS( FREQC-FRQ(I,L)).LT.THD ) THEN
                  IF (DOCAL) THEN
                     DOCAL = .FALSE.
                     IF (DIPLEN) THEN
                        WRITE(LUPRI,'(2(A,F9.6),2(A,6(A1)))')
     &              '@ B-freq =',FREQB,'  C-freq =',FREQC,
     &                   '     beta(',ALAB(1:1),';',BLAB(1:1),',',
     &                        CLAB(1:1),')',' = beta(',
     &                       LAB(I,1)(1:1),',',LAB(I,2)(1:1),',',
     &                       LAB(I,3)(1:1),')'
                     ELSE
                        WRITE(LUPRI,'(//A,3(/A,A10,I4,F10.6))')
     *     ' Quadratic response function value in a.u. for',
     *     ' A operator, symmetry, frequency: ',ALAB,ISYMA,FREQA,
     *     ' B operator, symmetry, frequency: ',BLAB,ISYMB,FREQB,
     *     ' C operator, symmetry, frequency: ',CLAB,ISYMC,FREQC
            WRITE(LUPRI,'(/A)') ' is equal to the already calculated'
            WRITE(LUPRI,'(/A,3(/A,A10,I4,F10.6))')
     *' Quadratic response function value in a.u. for',
     *' A operator, symmetry, frequency: ',LAB(I,1),ISYM(I,1),FRQ(I,1),
     *' B operator, symmetry, frequency: ',LAB(I,2),ISYM(I,2),FRQ(I,2),
     *' C operator, symmetry, frequency: ',LAB(I,3),ISYM(I,3),FRQ(I,3)
                     END IF
                  END IF
               END IF
            END IF
            END DO
         END IF
         END DO
         END DO
         END DO
C
C     Finally we check if the result is available from a previous
C     calculation
C     
         IF (DOCAL) THEN
            REWIND (LURSPRES)
 346        READ (LURSPRES,'(A75)',END=347) STRING
            IF (STRING(1:8) .EQ. '@ B-freq') THEN
               READ (STRING,'(10X,F9.6,10X,F9.6,10X,A1,1X,A1,1X,
     &              A1,3X,F16.9)') FRQB, FRQC, DIR(1), DIR(2), DIR(3),
     &              QUAVAL
               DO I = 1, 3
                  IF (DIR(I) .EQ. 'X') THEN
                     LABEL(I) = 'XDIPLEN '
                  ELSE IF (DIR(I) .EQ. 'Y') THEN
                     LABEL(I) = 'YDIPLEN '
                  ELSE
                     LABEL(I) = 'ZDIPLEN '
                  END IF
               END DO
            ELSE IF (STRING(1:8) .EQ. '@ Quadra') THEN
               READ (LURSPRES,'(30X,A8,2I5)') LABEL(1),KSYMA,KSPINA
               READ (LURSPRES,'(30X,A8,2I5)') LABEL(2),KSYMB,KSPINB
               READ (LURSPRES,'(30X,A8,2I5)') LABEL(3),KSYMC,KSPINC
               READ (LURSPRES,'(29X,3F15.8)') FRQB, FRQC, QUAVAL
            ELSE
                GOTO 346
             END IF
             IF ((LABEL(1) .EQ. ALAB) .AND. (LABEL(2) .EQ. BLAB) .AND.
     &            (LABEL(3) .EQ. CLAB) .AND.
     &            ABS(FREQB-FRQB).LT.THD
     &            .AND. ABS(FREQC-FRQC).LT.THD) THEN
                IF (DIPLEN) THEN
                   WRITE(LUPRI,'(/A)') ' The following first '//
     &                 'hyperpolarizability has already been calculated'
                   WRITE(LUPRI,'(2(A,F9.6),A,6(A1),A,F16.8)')
     &                  '@ B-freq =',FRQB,'  C-freq =',FRQC,
     &               '     beta(',DIR(1),';',DIR(2),',',DIR(3),')',
     &               ' =',QUAVAL
                ELSE
                   WRITE(LUPRI,'(/A)') ' The following quadratic '//
     &                 'response function has already been calculated'
                   WRITE(LUPRI,'(//A,3(/2A,2I5))')
     * '@ Quadratic response function value in a.u. for',
     * '@ A operator, symmetry, spin: ',LABEL(1),KSYMA,KSPINA,
     * '@ B operator, symmetry, spin: ',LABEL(2),KSYMB,KSPINB,
     * '@ C operator, symmetry, spin: ',LABEL(3),KSYMC,KSPINC
                   WRITE(LUPRI,'(/A,3F15.8)')
     *                  '@omega B, omega C, QR value :',
     *                  FRQB,FRQC,QUAVAL
                   IF (SOCOLL) CALL SODIST(ALAB,KSYMA,BLAB,KSYMB,
     *                  CLAB,KSYMC,QUAVAL)
                   IF (SSCOLL) CALL SSDIST(ALAB,KSYMA,BLAB,KSYMB,CLAB,
     *                  KSYMC,QUAVAL,SPNFC1,SPNFC2,SPNFC,SPNSD1,SPNSD2,
     *                  SPNSD,SPNSU1,SPNSU2,SPNSU)
                END IF
                CALL FLSHFO(LUPRI)
                DOCAL = .FALSE.
                ONFIL = .TRUE.
             END IF
             GOTO 346
 347         CONTINUE
#if defined (VAR_MFDS)
             BACKSPACE (LURSPRES)
#endif
          END IF
C
C If not done before, put in list
C

C     lf 1234 continue there is a minor bug here: for some components
C     both the number is given and afterwards beta is set equal to a
C     previously calculated component It affects just the printout, I
C     hope...
C     -- hjaaj 17oct07: discovered ONFIL was not always defined, maybe
C     that was the reason for this error ? ONFIL is now always
C     initialized to .false..

         IF (DOCAL .OR. ONFIL) THEN
            N = N + 1
            IF (N .GT. MXCALC) THEN
               WRITE (LUPRI,'(/A,I3)')
     &              'BCCHK: # unique calculations .gt. MXCALC =',MXCALC
               CALL QUIT('BCCHK: # unique calculations .gt. MXCALC')
            END IF
            LAB(N,1) = ALAB
            LAB(N,2) = BLAB
            LAB(N,3) = CLAB
            ISYM(N,1) = ISYMA
            ISYM(N,2) = ISYMB
            ISYM(N,3) = ISYMC
            FRQ(N,1) = FREQA
            FRQ(N,2) = FREQB
            FRQ(N,3) = FREQC
         END IF
C
         MEMFLG=(MEMFLG .OR. DOCAL)
C
      END DO
C
      DOCAL = MEMFLG
C
      RETURN
      END
      SUBROUTINE LRHYP(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *                 XINDX,VECA,VECB,MJWOP,WRK,LWRK)
C
#include <implicit.h>
#include <dummy.h>
#include <iratdef.h>
C
C PURPOSE:
C CALCULATION OF LINEAR POLARIZABILITIES IN A QR CALCULATION
C
      LOGICAL FOUND, CONV
      CHARACTER*8 BLANK
      PARAMETER ( D1 = 1.0D0, BLANK = '        ', D0 = 0.0D0 )
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),VECA(*),VECB(*),MJWOP(*),WRK(*)
C
#include <infrsp.h>
#include <maxorb.h>
#include <infvar.h>
#include <qrinf.h>
#include <priunit.h>
#include <rspprp.h>
#include <indqr.h>
#include <infhyp.h>
#include <inforb.h>
#include <infpri.h>
#include <infspi.h>
#include <wrkrsp.h>
#include <infhso.h>
#include <inftap.h>
#include <inflr.h>
C
      IPRRSP = IPRHYP
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KRES,MZYVMX,WRK,KFREE,LFREE)
C
      WRITE(LUPRI,'(A)') ' '
      CALL PRSYMB(LUPRI,'=',70,1)
      WRITE(LUPRI,'(A)') 
     &     ' >>>>>>>>    L I N E A R   R E S P O N S E   '//
     &     'F U N C T I O N   <<<<<<<<'
      CALL PRSYMB(LUPRI,'=',70,1)
C
      WRITE (LUPRI,'(/A/A//A/A,1P,D10.2)')
     &   ' The -<<A;B>>(omega_b) functions from vectors generated',
     &   ' in a *QUADRA calculation of <<A;B,C>>(omega_b,omega_c)',
     &   ' Note: the accuracy of off-diagonal elements will be linear',
     &   ' in the convergence threshold THCLR =',THCLR
      IF (ISPINA .NE. ISPINB) THEN
         WRITE (LUPRI,'(/A/A,2I4)')
     &   ' All zero because spin symmetry of A and B is different.',
     &   ' ISPINA and ISPINB : ',ISPINA,ISPINB
         GO TO 9999
      END IF
C
      DO ISYM = 1,NSYM
C
         IF ((NAQROP(ISYM).GT.0) .AND. (NBQROP(ISYM).GT.0)) THEN
C
C     Define variables that depend on symmetry
C
            KZYVAR = MZYVAR(ISYM)
            KSYMOP = ISYM
            CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
            CALL SETZY(MJWOP)
C
            DO IBOP = 1,NBQROP(ISYM)
            DO IBFR = 1,NBQRFR
C
C     Read in the linear response vector from file
C
               FRQ = BQRFR(IBFR)
               CALL REARSP(LURSP,KLEN,VECB,BQRLB(KSYMOP,IBOP),BLANK,
     &                     FRQ,D0,KSYMOP,0,THCLR,FOUND,CONV,ANTSYM)
               IF (.NOT. (FOUND .AND. CONV))
     &              CALL RIOERR(FOUND,BQRLB(KSYMOP,IBOP),FRQ,KSYMOP)
C
C     Print out response vector
C
               IF (IPRRSP.GT.10) THEN
                  WRITE(LUPRI,'(/A,3X,A,3X,A,2X,D12.6,/A,I5)')
     *                 ' LRHYP: B solution vector;  LABEL',
     *                 BQRLB(ISYM,IBOP),' FREQUENCY ',FRQ,
     *                 ' SYMMETRY',ISYM
                  CALL RSPPRC(VECB,MZCONF(ISYM),MZVAR(ISYM),LUPRI)
                  CALL RSPPRW(VECB(1+MZCONF(ISYM)),MJWOP,
     *                 MZWOPT(ISYM),ISYM,MZVAR(ISYM),LUPRI)
               END IF
C
               DO IAOP = 1,NAQROP(ISYM)
C
C     Get gradient
C
                  CALL GETGPV(AQRLB(ISYM,IAOP),FC,FV,CMO,UDV,PV,XINDX,
     &                        ANTSYM,WRK(KRES),LFREE)
                  CALL DCOPY(KZYVAR,WRK(KRES),1,VECA,1)
C
C     Print out gradient
C
                  IF (IPRRSP.GT.10) THEN
                     WRITE(LUPRI,'(/A,3X,A,3X,A,2X,D12.6,/A,I5)')
     *                    ' LRHYP: A gradient vector;  LABEL',
     *                    AQRLB(ISYM,IAOP),' FREQUENCY ',-FRQ,
     *                    ' SYMMETRY',ISYM
                     CALL RSPPRC(VECA,MZCONF(ISYM),MZVAR(ISYM),LUPRI)
                     CALL RSPPRW(VECA(1+MZCONF(ISYM)),MJWOP,
     *                    MZWOPT(ISYM),ISYM,MZVAR(ISYM),LUPRI)
                  END IF
C
C     Compute linear response value 
C
                  VAL = DDOT(KZYVAR,VECA,1,VECB,1)
C
                  IF (ISPINA .EQ. 0) THEN
                     WRITE(LUPRI,'(/A)')
     &               '@ Singlet linear response function in a.u.'
                  ELSE
                     WRITE(LUPRI,'(/A)')
     &               '@ Triplet linear response function in a.u.'
                  END IF
                  WRITE(LUPRI,'(2(/A,A10,I4,F10.6),/,/A,F20.12)')
     &               '@ A operator, symmetry, frequency: ',
     &               AQRLB(ISYM,IAOP),ISYM,-FRQ,
     &               '@ B operator, symmetry, frequency: ',
     &               BQRLB(ISYM,IBOP),ISYM, FRQ,
     &               '@ Value of linear response -<<A;B>(omega): ',VAL
               END DO
            END DO
            END DO
         END IF
      END DO
C
 9999 CONTINUE
      RETURN
      END
C  /* Deck sodist */
      SUBROUTINE SODIST(AOPLAB,IASYM,BOPLAB,IBSYM,COPLAB,ICSYM,VALUE)
C
C     This routine collects contributions to the SO-corrected shielding 
C     tensor in suitable arrays. K.Ruud, July-03
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <codata.h>
#include <gfac.h>
      LOGICAL ONESO, TWOSO, SUMSO, FRMONL, SDONL, FCSD
      CHARACTER*8 AOPLAB, BOPLAB, COPLAB
      CHARACTER DIRB*1, DIRSO*1, NUMNUC*2, NUCDIR*3, DIR2*1
#include <nuclei.h>
#include <symmet.h>
#include <sigmaso.h>
C
      ONESO = .FALSE.
      TWOSO = .FALSE.
      SUMSO = .FALSE.
      FRMONL = .FALSE.
      SDONL  = .FALSE.
      FCSD   = .FALSE.
      VALUE = VALUE*GFAC/8.0D0
C
      IF (AOPLAB(2:7).EQ.'ANGMOM') DIRB = AOPLAB(1:1)
      IF (BOPLAB(2:7).EQ.'ANGMOM') DIRB = BOPLAB(1:1)
      IF (COPLAB(2:7).EQ.'ANGMOM') DIRB = COPLAB(1:1)
C
      IF (AOPLAB(3:8).EQ.'SPNORB' .OR. AOPLAB(3:8) .EQ. 'MNF-SO') THEN
         DIRSO = AOPLAB(1:1)
         ONESO = AOPLAB(2:2).EQ.'1'
         TWOSO = AOPLAB(2:2).EQ.'2'
         SUMSO = AOPLAB(2:2).EQ.' '
         IREPSO = IASYM
      ELSE IF (BOPLAB(3:8).EQ.'SPNORB'.OR.BOPLAB(3:8).EQ.'MNF-SO') THEN
         DIRSO = BOPLAB(1:1)
         ONESO = BOPLAB(2:2).EQ.'1'
         TWOSO = BOPLAB(2:2).EQ.'2'
         SUMSO = BOPLAB(2:2).EQ.' '
         IREPSO = IBSYM
      ELSE IF (COPLAB(3:8).EQ.'SPNORB'.OR.BOPLAB(3:8).EQ.'MNF-SO') THEN
         DIRSO = COPLAB(1:1)
         ONESO = COPLAB(2:2).EQ.'1'
         TWOSO = COPLAB(2:2).EQ.'2'
         SUMSO = COPLAB(2:2).EQ.' '
         IREPSO = ICSYM
      END IF
C
      IF (AOPLAB(1:3) .EQ. 'FC ') THEN 
         NUMNUC = AOPLAB(7:8)
         FRMONL = .TRUE.
      ELSE IF (BOPLAB(1:3) .EQ. 'FC ') THEN 
         NUMNUC = BOPLAB(7:8)
         FRMONL = .TRUE.
      ELSE IF (COPLAB(1:3) .EQ. 'FC ') THEN 
         NUMNUC = COPLAB(7:8)
         FRMONL = .TRUE.
      END IF
      IF (AOPLAB(1:3) .EQ. 'SD ') THEN 
         NUCDIR = AOPLAB(4:6)
         DIR2 = AOPLAB(8:8)
         SDONL = .TRUE.
      ELSE IF (BOPLAB(1:3) .EQ. 'SD ') THEN 
         NUCDIR = BOPLAB(4:6)
         DIR2 = BOPLAB(8:8)
         SDONL = .TRUE.
      ELSE IF (COPLAB(1:3) .EQ. 'SD ') THEN 
         NUCDIR = COPLAB(4:6)
         DIR2 = COPLAB(8:8)
         SDONL = .TRUE.
      END IF
      IF (AOPLAB(1:3) .EQ. 'SDC') THEN 
         NUCDIR = AOPLAB(4:6)
         DIR2 = AOPLAB(8:8)
         FCSD = .TRUE.
      ELSE IF (BOPLAB(1:3) .EQ. 'SDC') THEN 
         NUCDIR = BOPLAB(4:6)
         DIR2 = BOPLAB(8:8)
         FCSD = .TRUE.
      ELSE IF (COPLAB(1:3) .EQ. 'SDC') THEN 
         NUCDIR = COPLAB(4:6)
         DIR2 = COPLAB(8:8)
         FCSD = .TRUE.
      END IF
C
C     Add contributions to the proper place of the shielding tensor
C     
      IF (DIRB .EQ. 'X') IDIRB = 1
      IF (DIRB .EQ. 'Y') IDIRB = 2
      IF (DIRB .EQ. 'Z') IDIRB = 3
      IF (DIRSO .EQ. 'X') IDIRSO = 1
      IF (DIRSO .EQ. 'Y') IDIRSO = 2
      IF (DIRSO .EQ. 'Z') IDIRSO = 3
      IF (FRMONL) THEN
         READ (NUMNUC,'(I2)') INUCNM
         ISCOOR = IPTCNT(3*(INUCNM - 1) + IDIRSO,IREPSO - 1,2)
         IF (ONESO) SIGFC1(IPTAX(IDIRB,2),ISCOOR) = 
     &                        SIGFC1(IPTAX(IDIRB,2),ISCOOR) - VALUE
         IF (TWOSO) SIGFC2(IPTAX(IDIRB,2),ISCOOR) = 
     &                        SIGFC2(IPTAX(IDIRB,2),ISCOOR) - VALUE
         IF (SUMSO) SIGFC(IPTAX(IDIRB,2),ISCOOR) = 
     &                         SIGFC(IPTAX(IDIRB,2),ISCOOR) - VALUE
      ELSE
         IF (DIR2 .EQ. 'x') IDIR2 = 1
         IF (DIR2 .EQ. 'y') IDIR2 = 2
         IF (DIR2 .EQ. 'z') IDIR2 = 3
         READ (NUCDIR, '(I3)') INUCDR 
         IF (IDIR2 .EQ. IDIRSO) THEN
            IF (SDONL) THEN
               IF (ONESO) SIGSD1(IPTAX(IDIRB,2),INUCDR) =
     &                         SIGSD1(IPTAX(IDIRB,2),INUCDR) - VALUE
               IF (TWOSO) SIGSD2(IPTAX(IDIRB,2),INUCDR) = 
     &                         SIGSD2(IPTAX(IDIRB,2),INUCDR) - VALUE
               IF (SUMSO) SIGSD(IPTAX(IDIRB,2),INUCDR) = 
     &                         SIGSD(IPTAX(IDIRB,2),INUCDR) - VALUE
            ELSE IF (FCSD) THEN
               IF (ONESO) SIGFS1(IPTAX(IDIRB,2),INUCDR) =
     &                         SIGFS1(IPTAX(IDIRB,2),INUCDR) - VALUE
               IF (TWOSO) SIGFS2(IPTAX(IDIRB,2),INUCDR) =
     &                         SIGFS2(IPTAX(IDIRB,2),INUCDR) - VALUE
               IF (SUMSO) SIGFS(IPTAX(IDIRB,2),INUCDR) =
     &                         SIGFS(IPTAX(IDIRB,2),INUCDR) - VALUE
            END IF
         END IF
      END IF
C
      RETURN
      END
C  /* Deck soinit */
      SUBROUTINE SOINIT
#include <implicit.h>
#include <mxcent.h>
#include <sigmaso.h>
C
      CALL DZERO(SIGFC1,3*MXCOOR)
      CALL DZERO(SIGFC2,3*MXCOOR)
      CALL DZERO(SIGFC,3*MXCOOR)
      CALL DZERO(SIGSD1,3*MXCOOR)
      CALL DZERO(SIGSD2,3*MXCOOR)
      CALL DZERO(SIGSD,3*MXCOOR)
      CALL DZERO(SIGFS1,3*MXCOOR)
      CALL DZERO(SIGFS2,3*MXCOOR)
      CALL DZERO(SIGFS,3*MXCOOR)
      RETURN
      END
C  /* Deck sosout */
      SUBROUTINE SOSOUT(SMAT1,SMAT2,CSTRA,SCTRA)
C
C     Collects and prints the final results of the SO-corrections to the
C     shielding constants. K.Ruud, July-03
C
#include <implicit.h>
#include <priunit.h>
#include <thrzer.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <codata.h>
      DIMENSION SMAT1(3,3,MXCENT), SMAT2(3,3,MXCENT), CSTRA(*), SCTRA(*)
#include <abainf.h>
#include <nuclei.h>
#include <symmet.h>
#include <sigmaso.h>
#include <ibtfun.h>
C
      FC1SUM = DNORM2(3*MXCOOR,SIGFC1,1)
      FC2SUM = DNORM2(3*MXCOOR,SIGFC2,1)
      FCSUM = DNORM2(3*MXCOOR,SIGFC,1)
      SD1SUM = DNORM2(3*MXCOOR,SIGSD1,1)
      SD2SUM = DNORM2(3*MXCOOR,SIGSD2,1)
      SDSUM = DNORM2(3*MXCOOR,SIGSD,1)
      FS1SUM = DNORM2(3*MXCOOR,SIGFS1,1)
      FS2SUM = DNORM2(3*MXCOOR,SIGFS2,1)
      FSSUM = DNORM2(3*MXCOOR,SIGFS,1)
      DO I = 1, 8
         DOSYM(I) = .TRUE.
      END DO
C
      IF (FCSUM .LT. THRZER .AND. (FC1SUM .GT. THRZER .AND. 
     &                             FC2SUM .GT. THRZER)) THEN
         DO I = 1, 3*MXCOOR
            SIGFC(I,1) = SIGFC1(I,1) + SIGFC2(I,1)
         END DO
         FCSUM = DNORM2(3*MXCOOR,SIGFC,1)
      END IF
      IF (SDSUM .LT. THRZER .AND. (SD1SUM .GT. THRZER .AND. 
     &                             SD2SUM .GT. THRZER)) THEN
         DO I = 1, 3*MXCOOR
            SIGSD(I,1) = SIGSD1(I,1) + SIGSD2(I,1)
         END DO
         SDSUM = DNORM2(3*MXCOOR,SIGSD,1)
      END IF
      IF (FSSUM .LT. THRZER .AND. (FS1SUM .GT. THRZER .AND. 
     &                             FS2SUM .GT. THRZER)) THEN
         DO I = 1, 3*MXCOOR
            SIGFS(I,1) = SIGFS1(I,1) + SIGFS2(I,1)
         END DO
         FSSUM = DNORM2(3*MXCOOR,SIGFS,1)
      END IF
C
      CALL TITLER('ABACUS - SO CORRECTIONS TO CHEMICAL SHIELDINGS',
     &     '*',124)
      IF (FC1SUM .GT. THRZER) THEN
         CALL HEADER('FC(1) contribution to shielding tensors '//
     &        'in symmetry coordinates (ppm)',-1)
         CALL FCPRI(SIGFC1,'SIGMA ',CSTRA,SCTRA)
      END IF
      IF (FC2SUM .GT. THRZER) THEN
         CALL HEADER('FC(2) contribution to shielding tensors '//
     &        'in symmetry coordinates (ppm)',-1)
         CALL FCPRI(SIGFC2,'SIGMA ',CSTRA,SCTRA)
      END IF
      IF (FCSUM .GT. THRZER) THEN
         CALL HEADER('FC contribution to shielding tensors '//
     &        'in symmetry coordinates (ppm)',-1) 
         CALL FCPRI(SIGFC,'SIGMA ',CSTRA,SCTRA)
      END IF
      IF (SD1SUM .GT. THRZER) THEN
         CALL HEADER('SD(1) contribution to shielding tensors '//
     &        'in symmetry coordinates (ppm)',-1)
         CALL FCPRI(SIGSD1,'SIGMA ',CSTRA,SCTRA)
      END IF
      IF (SD2SUM .GT. THRZER) THEN
         CALL HEADER('SD(2) contribution to shielding tensors '//
     &        'in symmetry coordinates (ppm)',-1)
         CALL FCPRI(SIGSD2,'SIGMA ',CSTRA,SCTRA)
      END IF
      IF (SDSUM .GT. THRZER) THEN
         CALL HEADER('SD contribution to shielding tensors '//
     &        'in symmetry coordinates (ppm)',-1)
         CALL FCPRI(SIGSD,'SIGMA ',CSTRA,SCTRA)
      END IF
      IF (FS1SUM .GT. THRZER) THEN
         CALL HEADER('SD(1)+FC(1) contribution to shielding tensors '//
     &        'in symmetry coordinates (ppm)',-1)
         CALL FCPRI(SIGFS1,'SIGMA ',CSTRA,SCTRA)
      END IF
      IF (FS2SUM .GT. THRZER) THEN
         CALL HEADER('SD(2)+FC(2) contribution to shielding tensors '//
     &        'in symmetry coordinates (ppm)',-1)
         CALL FCPRI(SIGFS2,'SIGMA ',CSTRA,SCTRA)
      END IF
      IF (FSSUM .GT. THRZER) THEN
         CALL HEADER('SD+FC contribution to shielding tensors '//
     &        'in symmetry coordinates (ppm)',-1)
         CALL FCPRI(SIGFS,'SIGMA ',CSTRA,SCTRA)
      END IF
C
      CALL TRADIP(SIGFC1,SMAT1,CSTRA,SCTRA,3*NUCDEP,2,2)
      CALL TRADIP(SIGFC2,SMAT2,CSTRA,SCTRA,3*NUCDEP,2,2)
      CALL AROUND('Summary of isotropic SO-FC contribution to '//
     &            'chemical shieldings')
      WRITE (LUPRI,'(1X,A,/,1X,A)')
      WRITE (LUPRI,'(1X,A,/,1X,A)')
     & '@1atom   1-electron FC    2-electron FC    Total FC',
     & '@1----------------------------------------------'/
     &  /'-------------'
      IATOM = 0
      FAC= ALPHA2*1.D6/3.0D0
      DO I = 1, NUCIND
         DO ISYMOP = 0, MAXOPR
            IF (IBTAND(ISTBNU(I),ISYMOP).EQ.0) THEN
               IATOM = IATOM + 1
               SIGMA1 = (SMAT1(1,1,IATOM) + SMAT1(2,2,IATOM)
     &                + SMAT1(3,3,IATOM))*FAC
               SIGMA2 = (SMAT2(1,1,IATOM) + SMAT2(2,2,IATOM)
     &                + SMAT2(3,3,IATOM))*FAC
               SIGMAT = SIGMA1 + SIGMA2
               WRITE (LUPRI,'(1X,2A,3F15.9)') '@1',
     &              NAMDEP(IATOM), SIGMA1, SIGMA2, SIGMAT
            END IF
         END DO
      END DO
C
      CALL TRADIP(SIGSD1,SMAT1,CSTRA,SCTRA,3*NUCDEP,2,2)
      CALL TRADIP(SIGSD2,SMAT2,CSTRA,SCTRA,3*NUCDEP,2,2)
      CALL AROUND('Summary of isotropic SO-SD contribution to '//
     &            'chemical shieldings')
      WRITE (LUPRI,'(1X,A,/,1X,A)')
     & '@1atom   1-electron SD    2-electron SD    Total SD',
     & '@1----------------------------------------------'/
     &  /'-------------'
      IATOM = 0
      DO I = 1, NUCIND
         DO ISYMOP = 0, MAXOPR
            IF (IBTAND(ISTBNU(I),ISYMOP).EQ.0) THEN
               IATOM = IATOM + 1
               SIGMA1 = (SMAT1(1,1,IATOM) + SMAT1(2,2,IATOM)
     &                + SMAT1(3,3,IATOM))*FAC
               SIGMA2 = (SMAT2(1,1,IATOM) + SMAT2(2,2,IATOM)
     &                + SMAT2(3,3,IATOM))*FAC
               SIGMAT = SIGMA1 + SIGMA2
               WRITE (LUPRI,'(1X,2A,3F15.9)') '@1',
     &              NAMDEP(IATOM), SIGMA1, SIGMA2, SIGMAT
            END IF
         END DO
      END DO
C
      RETURN
      END
C  /* Deck ssdist */
      SUBROUTINE SSDIST(AOPLAB,IASYM,BOPLAB,IBSYM,COPLAB,ICSYM,VALUE,
     &                  SPNFC1,SPNFC2,SPNFC,SPNSD1,SPNSD2,SPNSD,
     &                  SPNSU1,SPNSU2,SPNSU)
C
C     This routine collects contributions to the SO-corrected shielding 
C     tensor in suitable arrays. K.Ruud, July-03
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <codata.h>
#include <gfac.h>
      PARAMETER (AUTOHZ = ALPHA2*ALPHA2*XTHZ/(XFAMU*XFAMU))
      DIMENSION SPNFC1(MXCOOR,MXCOOR), SPNFC2(MXCOOR,MXCOOR),
     &          SPNFC (MXCOOR,MXCOOR), SPNSD1(MXCOOR,MXCOOR),
     &          SPNSD2(MXCOOR,MXCOOR), SPNSD (MXCOOR,MXCOOR),
     &          SPNSU1(MXCOOR,MXCOOR), SPNSU2(MXCOOR,MXCOOR),
     &          SPNSU(MXCOOR,MXCOOR)
      LOGICAL ONESO, TWOSO, SUMSO, FRMONL, SDONL, FCSD
      CHARACTER*8 AOPLAB, BOPLAB, COPLAB
      CHARACTER DIRSO*1, NUMNUC*2, NUCDIR*3, DIR2*1
#include <nuclei.h>
#include <symmet.h>
#include <ibtfun.h>
C
      ONESO = .FALSE.
      TWOSO = .FALSE.
      SUMSO = .FALSE.
      FRMONL = .FALSE.
      SDONL  = .FALSE.
      FCSD   = .FALSE.
      VALUE = -VALUE*AUTOHZ*GFAC/(2.0D0*PI)
C
      IF (AOPLAB(1:3).EQ.'PSO') READ (AOPLAB(5:7),'(I3)') IPSO
      IF (BOPLAB(1:3).EQ.'PSO') READ (BOPLAB(5:7),'(I3)') IPSO
      IF (COPLAB(1:3).EQ.'PSO') READ (COPLAB(5:7),'(I3)') IPSO
C
      IF (AOPLAB(3:8).EQ.'SPNORB' .OR. AOPLAB(3:8) .EQ. 'MNF-SO') THEN
         DIRSO = AOPLAB(1:1)
         ONESO = AOPLAB(2:2).EQ.'1'
         TWOSO = AOPLAB(2:2).EQ.'2'
         SUMSO = AOPLAB(2:2).EQ.' '
         IREPSO = IASYM
      ELSE IF (BOPLAB(3:8).EQ.'SPNORB'.OR.BOPLAB(3:8).EQ.'MNF-SO') THEN
         DIRSO = BOPLAB(1:1)
         ONESO = BOPLAB(2:2).EQ.'1'
         TWOSO = BOPLAB(2:2).EQ.'2'
         SUMSO = BOPLAB(2:2).EQ.' '
         IREPSO = IBSYM
      ELSE IF (COPLAB(3:8).EQ.'SPNORB'.OR.BOPLAB(3:8).EQ.'MNF-SO') THEN
         DIRSO = COPLAB(1:1)
         ONESO = COPLAB(2:2).EQ.'1'
         TWOSO = COPLAB(2:2).EQ.'2'
         SUMSO = COPLAB(2:2).EQ.' '
         IREPSO = ICSYM
      END IF
C
      IF (AOPLAB(1:3) .EQ. 'FC ') THEN 
         NUMNUC = AOPLAB(7:8)
         FRMONL = .TRUE.
      ELSE IF (BOPLAB(1:3) .EQ. 'FC ') THEN 
         NUMNUC = BOPLAB(7:8)
         FRMONL = .TRUE.
      ELSE IF (COPLAB(1:3) .EQ. 'FC ') THEN 
         NUMNUC = COPLAB(7:8)
         FRMONL = .TRUE.
      END IF
      IF (AOPLAB(1:3) .EQ. 'SD ') THEN 
         NUCDIR = AOPLAB(4:6)
         DIR2 = AOPLAB(8:8)
         SDONL = .TRUE.
      ELSE IF (BOPLAB(1:3) .EQ. 'SD ') THEN 
         NUCDIR = BOPLAB(4:6)
         DIR2 = BOPLAB(8:8)
         SDONL = .TRUE.
      ELSE IF (COPLAB(1:3) .EQ. 'SD ') THEN 
         NUCDIR = COPLAB(4:6)
         DIR2 = COPLAB(8:8)
         SDONL = .TRUE.
      END IF
      IF (AOPLAB(1:3) .EQ. 'SDC') THEN 
         NUCDIR = AOPLAB(4:6)
         DIR2 = AOPLAB(8:8)
         FCSD = .TRUE.
      ELSE IF (BOPLAB(1:3) .EQ. 'SDC') THEN 
         NUCDIR = BOPLAB(4:6)
         DIR2 = BOPLAB(8:8)
         FCSD = .TRUE.
      ELSE IF (COPLAB(1:3) .EQ. 'SDC') THEN 
         NUCDIR = COPLAB(4:6)
         DIR2 = COPLAB(8:8)
         FCSD = .TRUE.
      END IF
C
C     Add contributions to the proper place of the shielding tensor
C     
      IF (DIRSO .EQ. 'X') IDIRSO = 1
      IF (DIRSO .EQ. 'Y') IDIRSO = 2
      IF (DIRSO .EQ. 'Z') IDIRSO = 3
      IF (FRMONL) THEN
         READ (NUMNUC,'(I2)') INUCNM
         DO IREP = 0, MAXREP
            DO IATOM = 1, NUCIND
               IF (IPTNUC(IATOM,IREP) .EQ. INUCNM) THEN
                  IREPF = IBTXOR(IREPSO - 1,IREP)
                  ISCOOR = IPTCNT(3*(IATOM - 1) + IDIRSO,IREPF,2)
               END IF
            END DO
         END DO
         IF (ONESO) SPNFC1(IPSO,ISCOOR) = SPNFC1(IPSO,ISCOOR) + VALUE
         IF (TWOSO) SPNFC2(IPSO,ISCOOR) = SPNFC2(IPSO,ISCOOR) + VALUE
         IF (SUMSO) SPNFC(IPSO,ISCOOR) = SPNFC(IPSO,ISCOOR) + VALUE
         IF (ONESO) SPNFC1(ISCOOR,IPSO) = SPNFC1(ISCOOR,IPSO) + VALUE
         IF (TWOSO) SPNFC2(ISCOOR,IPSO) = SPNFC2(ISCOOR,IPSO) + VALUE
         IF (SUMSO) SPNFC(ISCOOR,IPSO) = SPNFC(ISCOOR,IPSO) + VALUE
      ELSE
         IF (DIR2 .EQ. 'x') IDIR2 = 1
         IF (DIR2 .EQ. 'y') IDIR2 = 2
         IF (DIR2 .EQ. 'z') IDIR2 = 3
         READ (NUCDIR, '(I3)') INUCDR 
         IF (IDIR2 .EQ. IDIRSO) THEN
            IF (SDONL) THEN
               IF (ONESO) SPNSD1(IPSO,INUCDR)=SPNSD1(IPSO,INUCDR)+VALUE
               IF (TWOSO) SPNSD2(IPSO,INUCDR)=SPNSD2(IPSO,INUCDR)+VALUE
               IF (SUMSO) SPNSD(IPSO,INUCDR)=SPNSD(IPSO,INUCDR) + VALUE
               IF (ONESO) SPNSD1(INUCDR,IPSO)=SPNSD1(INUCDR,IPSO)+VALUE
               IF (TWOSO) SPNSD2(INUCDR,IPSO)=SPNSD2(INUCDR,IPSO)+VALUE
               IF (SUMSO) SPNSD(INUCDR,IPSO)=SPNSD(INUCDR,IPSO) + VALUE
            ELSE IF (FCSD) THEN
               IF (ONESO) SPNSU1(IPSO,INUCDR)=SPNSU1(IPSO,INUCDR)+VALUE
               IF (TWOSO) SPNSU2(IPSO,INUCDR)=SPNSU2(IPSO,INUCDR)+VALUE
               IF (SUMSO) SPNSU(IPSO,INUCDR)=SPNSU(IPSO,INUCDR) + VALUE
               IF (ONESO) SPNSU1(INUCDR,IPSO)=SPNSU1(INUCDR,IPSO)+VALUE
               IF (TWOSO) SPNSU2(INUCDR,IPSO)=SPNSU2(INUCDR,IPSO)+VALUE
               IF (SUMSO) SPNSU(INUCDR,IPSO)=SPNSU(INUCDR,IPSO) + VALUE
            END IF
         END IF
      END IF
C
      RETURN
      END
C  /* Deck ssinit */
      SUBROUTINE SSINIT(SPNFC1,SPNFC2,SPNFC,SPNSD1,SPNSD2,SPNSD,SPNSU1,
     &                  SPNSU2,SPNSU)
#include <implicit.h>
#include <mxcent.h>
      DIMENSION SPNFC1(*), SPNFC2(*), SPNFC(*), SPNSD1(*), SPNSD2(*),
     &          SPNSD(*), SPNSU1(*), SPNSU2(*), SPNSU(*)
C
      CALL DZERO(SPNFC1,MXCOOR**2)
      CALL DZERO(SPNFC2,MXCOOR**2)
      CALL DZERO(SPNFC,MXCOOR**2)
      CALL DZERO(SPNSD1,MXCOOR**2)
      CALL DZERO(SPNSD2,MXCOOR**2)
      CALL DZERO(SPNSD,MXCOOR**2)
      CALL DZERO(SPNSU1,MXCOOR**2)
      CALL DZERO(SPNSU2,MXCOOR**2)
      CALL DZERO(SPNSU,MXCOOR**2)
C
      RETURN
      END
C  /* Deck spnso */
      SUBROUTINE SPNSO(SPNFC1,SPNFC2,SPNFC,SPNSD1,SPNSD2,SPNSD,SPNSU1,
     &                 SPNSU2,SPNSU,CSTRA,SCTRA)
C
C     Collects and prints the final results of the SO-corrections to the
C     shielding constants. K.Ruud, July-03
C
#include <implicit.h>
#include <priunit.h>
#include <thrzer.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <codata.h>
      DIMENSION CSTRA(*), SCTRA(*)
      DIMENSION SPNFC1(MXCOOR,MXCOOR), SPNFC2(MXCOOR,MXCOOR),
     &          SPNFC(MXCOOR,MXCOOR), SPNSD1(MXCOOR,MXCOOR),
     &          SPNSD2(MXCOOR,MXCOOR), SPNSD(MXCOOR,MXCOOR),
     &          SPNSU1(MXCOOR,MXCOOR), SPNSU2(MXCOOR,MXCOOR),
     &          SPNSU(MXCOOR,MXCOOR)
#include <abainf.h>
#include <nuclei.h>
#include <symmet.h>
#include <sigmaso.h>
#include <ibtfun.h>
C
      FC1SUM = DNORM2(MXCOOR**2,SPNFC1,1)
      FC2SUM = DNORM2(MXCOOR**2,SPNFC2,1)
      FCSUM = DNORM2(MXCOOR**2,SPNFC,1)
      SD1SUM = DNORM2(MXCOOR**2,SPNSD1,1)
      SD2SUM = DNORM2(MXCOOR**2,SPNSD2,1)
      SDSUM = DNORM2(MXCOOR**2,SPNSD,1)
      FS1SUM = DNORM2(MXCOOR**2,SPNSU1,1)
      FS2SUM = DNORM2(MXCOOR**2,SPNSU2,1)
      FSSUM = DNORM2(MXCOOR**2,SPNSU,1)
      DO I = 1, 8
         DOSYM(I) = .TRUE.
      END DO
C
      IF (FCSUM .LT. THRZER .AND. (FC1SUM .GT. THRZER .AND. 
     &                             FC2SUM .GT. THRZER)) THEN
         DO I = 1, MXCOOR*MXCOOR
            SPNFC(I,1) = SPNFC1(I,1) + SPNFC2(I,1)
         END DO
         FCSUM = DNORM2(MXCOOR**2,SPNFC,1)
      END IF
      IF (SDSUM .LT. THRZER .AND. (SD1SUM .GT. THRZER .AND. 
     &                             SD2SUM .GT. THRZER)) THEN
         DO I = 1, MXCOOR*MXCOOR
            SPNSD(I,1) = SPNSD1(I,1) + SPNSD2(I,1)
         END DO
         SDSUM = DNORM2(MXCOOR**2,SPNSD,1)
      END IF
      IF (FSSUM .LT. THRZER .AND. (FS1SUM .GT. THRZER .AND. 
     &                             FS2SUM .GT. THRZER)) THEN
         DO I = 1, MXCOOR*MXCOOR
            SPNSU(I,1) = SPNSU1(I,1) + SPNSU2(I,1)
         END DO
         FSSUM = DNORM2(MXCOOR**2,SPNSU,1)
      END IF
C
      CALL TITLER('ABACUS - SO CORRECTIONS TO REDUCED J COUPLINGS',
     &     '*',124)
      IF (FC1SUM .GT. THRZER) THEN
         CALL HEADER('FC(1) contribution to reduced J couplings ',-1)
         CALL PRIHES (SPNFC1,'SPNSPN',CSTRA,SCTRA)
      END IF
      IF (FC2SUM .GT. THRZER) THEN
         CALL HEADER('FC(2) contribution to reduced J couplings ',-1)
         CALL PRIHES (SPNFC2,'SPNSPN',CSTRA,SCTRA)
      END IF
      IF (FCSUM .GT. THRZER) THEN
         CALL HEADER('FC contribution to reduced J couplings ',-1)
         CALL PRIHES (SPNFC,'SPNSPN',CSTRA,SCTRA)
      END IF
      IF (SD1SUM .GT. THRZER) THEN
         CALL HEADER('SD(1) contribution to reduced J couplings ',-1)
         CALL PRIHES (SPNSD1,'SPNSPN',CSTRA,SCTRA)
      END IF
      IF (SD2SUM .GT. THRZER) THEN
         CALL HEADER('SD(2) contribution to reduced J couplings ',-1)
         CALL PRIHES (SPNSD2,'SPNSPN',CSTRA,SCTRA)
      END IF
      IF (SDSUM .GT. THRZER) THEN
         CALL HEADER('SD contribution to reduced J couplings ',-1)
         CALL PRIHES (SPNSD,'SPNSPN',CSTRA,SCTRA)
      END IF
      IF (FS1SUM .GT. THRZER) THEN
         CALL HEADER('FC+SD(1) contribution to reduced J couplings ',-1)
         CALL PRIHES (SPNSU1,'SPNSPN',CSTRA,SCTRA)
      END IF
      IF (FS2SUM .GT. THRZER) THEN
         CALL HEADER('FC+SD(2) contribution to reduced J couplings ',-1)
         CALL PRIHES (SPNSU2,'SPNSPN',CSTRA,SCTRA)
      END IF
      IF (FSSUM .GT. THRZER) THEN
         CALL HEADER('FC+SD contribution to reduced J couplings ',-1)
         CALL PRIHES (SPNSU,'SPNSPN',CSTRA,SCTRA)
      END IF
C
      RETURN
      END
