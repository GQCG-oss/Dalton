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
C
#ifdef REV_LOG
C===========================================================================
CRevision 1.4  2000/05/24 18:45:09  hjj
Cnew getref calls with appropriate NDREF (fixing CSF and triplet)
Cstop if solvent and CSF
Cbugfix: too little was allocated for KGP2
C===========================================================================
#endif
C  /* Deck lagran */
      SUBROUTINE LAGRAN(WORD,FC,FV,CMO,UDV,PV,XINDX,WRK,LWRK)
C
C 18-AUG-1991
C
C Called from GETGPV when label(5:8) = 'LAGR'
C
#include "implicit.h"
#include "dummy.h"
      DIMENSION CMO(*),UDV(NASHT,*),FC(*),FV(*),PV(*),XINDX(*),WRK(*)
      CHARACTER*8 WORD
C-- common blocks:
#include "maxorb.h"
#include "priunit.h"
#include "infdim.h"
#include "infvar.h"
#include "inforb.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "infpri.h"
#include "infinp.h"
#include <dftcom.h>
#include <pcmlog.h>
C
      PARAMETER (DHALF=0.5D0,D1=1.0D0,DM1=-1.0D0)
C
      CALL QENTER('LAGRAN')
C
      IF (IPRRSP .GE. 10) WRITE (LUPRI,7010) WORD
 7010 FORMAT(//' Output from LAGRAN:'/' ==================='
     *       //' Property label = ',A/)
      IF (SOPPA) THEN
         WRITE (LUPRI,*) 'SOPPA not implemented in LAGRAN yet, sorry!'
         CALL QUIT('LAGRAN-ERROR: SOPPA not implemented')
      END IF
clf      IF ((FLAG(16).OR.PCM).AND.(NASHT.GT.1).AND.(NCREF.NE.KZCONF)) THEN
Chj:     see comments in RSPSLV in rspsol.F
clf         WRITE(LUPRI,*)'Solvent not implemented in LAGRAN with CSFs'
clf         WRITE(LUPRI,*)'Use .DETERMINANTS and try again!'
clf         CALL QUIT('Solvent not implemented in LAGRAN with CSFs')
clf      END IF
C
C ALLOCATE WORK SPACE
C
      KGP   = 1
      KCREF = KGP + KZYVAR
      NDREF = MAX(KZCONF,NCREF)
Chj   ... RSPOLI requires determinants for ZYCVEC when triplet,
Chj   ... thus KZCONF. However, for ROHF NCREF = 1, KZCONF = 0
Chj       and we need CREF(1) for RSPDM. Thus MAX(KZCONF,NCREF)
      KTOT  = KCREF  + NDREF
C
      CALL DZERO(WRK(KGP),KZYVAR)
      CALL GETREF(WRK(KCREF),NDREF)
C
C CONSTRUCT the Lagrangian vector, the GP vector,
C using routine for the ORBITAL PART OF LINEAR TRANSFORMED VECTORS
C
      IF (.NOT.RSPCI) THEN
C
C     ALLOCATE WORK SPACE FOR RSPOLI
C        NCSIM  = 1, NOSIM = 0
C
         KFVTD  = KTOT
         LFVTD  = N2ORBX
         IF (DOMCSRDFT) LFVTD = 2*LFVTD
C     ... Extra allocation for "FCTD" in MCSCF-SRDFT,
C         is needed for the VxcTD matrix.
         KQATD  = KFVTD  + LFVTD
         KQBTD  = KQATD  + NORBT * NASHT
         KWRK1  = KQBTD  + NORBT * NASHT
         LWRK1  = LWRK   - KWRK1
         IF (LWRK1.LT.0) CALL ERRWRK('LAGRAN',KWRK1-1,LWRK)
         IF (DOHFSRDFT .OR. DOMCSRDFT) THEN
            CALL QUIT('LAGRANGE is not implemented yet for srDFT')
         END IF
         IF ((DFTADD.OR.HSROHF).AND.TRPLET) THEN
            CALL DAXPY(NNORBT,D1,FC,1,FV,1)
            CALL DSCAL(NNORBT,DM1,FV,1)
         END IF
         CALL RSPOLI(1,0,UDV,WRK(KCREF),NDREF,FC,FV,PV,DUMMY,
     *            DUMMY,DUMMY, DUMMY,DUMMY,WRK(KFVTD),
     *            WRK(KQATD),WRK(KQBTD),WRK(KGP),
     *            XINDX,CMO,WRK(KWRK1),LWRK1)
C        CALL RSPOLI(NCSIM,NOSIM,UDV,ZYCVEC,LZYCVEC,FC,FV,PVX,ZYMAT,
C    *            FCONE,FVONE,QAONE,QBONE,FVTD,QATD,QBTD,EVECS,
C    *            XINDX,CMO,WRK,LWRK)
         IF ((DFTADD.OR.HSROHF).AND.TRPLET) THEN
            CALL DAXPY(NNORBT,D1,FC,1,FV,1)
            CALL DSCAL(NNORBT,DM1,FV,1)
         END IF
C
      END IF
      IF (IPRRSP.GT.120) THEN
         WRITE(LUPRI,'(/A)') ' LAGRAN: LAGRANGIAN VECTOR '
         CALL OUTPUT(WRK(KGP),1,KZVAR,1,2,KZVAR,2,1,LUPRI)
      END IF
C
C     section for calculating esr-properties with solvent
C     contributions.
C
      IF (FLAG(16).OR.PCM) THEN
C
C        workspace allocation
C
         KGP2   = KTOT
         KTDV   = KGP2  + NVARH
Chj      KTDV   = KGP2  + KZVAR
Chj      ... SOLGDT uses NVARH (.gt. KZVAR !) for GP2
         KDV    = KTDV  + N2ASHX
         KWRK1  = KDV   + NNASHX
         IF (FLAG(16)) THEN
            KWRK2  = KWRK1 + NLMSOL * 2
         ELSE IF (PCM) THEN
            KWRK2  = KWRK1
         ELSE
            KWRK2  = KWRK1
         END IF
         LWRK1  = LWRK  - KWRK2
C
C        initialize solvent gradient.
C
         CALL DZERO(WRK(KGP2),NVARH)
C
C
         ISPIN1 = 1
         ISPIN2 = 0
         CALL RSPDM(IREFSY,IREFSY,NDREF,NDREF,WRK(KCREF),WRK(KCREF),
     *              WRK(KTDV),DUMMY, ISPIN1,ISPIN2,.FALSE.,.TRUE.,
     *              XINDX,WRK,KWRK2,LWRK1)
C        CALL RSPDM(ISYM,ISYM,NCDIM.NCDIM,REF,REF,
C                   TUDV,DUMMY, ISPIN1,ISPIN2,TDM,NORHO2,
C                   XINDX,WRK,KFREE,LFREE)
         IF (IPRRSP.GT.110) THEN
            WRITE(LUPRI,'(/A)')'LAGRANG: ONE ELECTRON SPIN DENSITY'
            CALL OUTPUT(WRK(KTDV),1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         END IF
C
C TRIANGULAR PACKING OF ONE ELECTRON DENSITY MATRIX
C
         CALL DGETSP(NASHT,WRK(KTDV),WRK(KDV))
C
         IF (IPRRSP.GT.50 .AND. NASHT.GT.0) THEN
           WRITE(LUPRI,'(/A)')
     &        'LAGRANG solvent: PACKED ONE ELECTRON SPIN DENSITY MATRIX'
           CALL OUTPAK(WRK(KDV),NASHT,1,LUPRI)
         END IF
C
C        use wavefunction routine to calculate solvent gradient.
C
         IF (FLAG(16)) THEN
            CALL SOLGDT(WRK(KCREF),CMO,XINDX,WRK(KDV),WRK(KGP2),
     *           ESOLT,WRK(KWRK1),WRK(KWRK2),LWRK1,NDREF,IREFSY)
         ELSE IF (PCM) THEN
            CALL PCMGDT(WRK(KCREF),CMO,XINDX,WRK(KDV),WRK(KGP2),
     *           ESOLT,WRK(KWRK2),LWRK1,NDREF,IREFSY)
         ELSE
            CALL QUIT("Lagran: Either FLAG(16) or PCM must be TRUE")
         END IF
               
C     
C        scale by a half for consistency with response and add to GP.
C
         CALL DAXPY(KZWOPT,DHALF,WRK(KGP2+KZCONF),1,WRK(KGP+KZCONF),1)
         CALL DAXPY(KZWOPT,(-DHALF),WRK(KGP2+KZCONF),1,
     *              WRK(KGP+KZVAR+KZCONF),1)
C
         IF (IPRRSP.GT.120) THEN
            WRITE(LUPRI,'(/A)') ' LAGRAN: 2 * SOLVENT LAGRANGIAN VECTOR'
            CALL OUTPUT(WRK(KGP2),1,NVARH,1,1,NVARH,1,1,LUPRI)
            WRITE(LUPRI,'(/A)')
     &      ' LAGRAN: Electronic + solvent LAGRANGIAN VECTOR'
            CALL OUTPUT(WRK(KGP),1,KZVAR,1,2,KZVAR,2,1,LUPRI)
         END IF
C
      ENDIF
C
C *** END OF LAGRAN
C
      CALL QEXIT('LAGRAN')
      RETURN
      END
C  /* Deck rsplan */
      SUBROUTINE RSPLAN(CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
C
C     Lagrangian solution vector is returned in WRK(1)
C
#include "implicit.h"
#include "dummy.h"
#include "iratdef.h"
C
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
#include "priunit.h"
#include "infpri.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "rspprp.h"
#include "inflr.h"
#include "infdim.h"
#include "inforb.h"
#include "infesr.h"
C
C Local variables
C
      CHARACTER*8 LABEL, BLANK
      PARAMETER ( D0 = 0.0D0 , DM8 = 1.0D-8, BLANK = '        ' )
C
      CALL QENTER('RSPLAN')
C
C DETERMINE SECOND ORDER MOLECULAR PROPERTIES
C
      KREDE  = 1
      KREDS  = KREDE  + MAXRM*MAXRM
      KIBTYP = KREDS  + MAXRM*MAXRM
      KEIVAL = KIBTYP + MAXRM
      KRESID = KEIVAL + MAXRM
      KEIVEC = KRESID + MAXRM
      KREDGP = KEIVEC + MAXRM*MAXRM
      KGP    = KREDGP + MAXRM
      KWRK1  = KGP    + KZYVAR
      LWRK1  = LWRK   - KWRK1
C
      IF (LWRK1.LT.3*KZYVAR) THEN
         WRITE (LUERR,9100) LWRK1,3*KZYVAR
         CALL QTRACE(LUERR)
         CALL QUIT('RSPLAN ERROR, INSUFFICIENT MEMORY')
      ENDIF
C
C WORK SPACE FOR RSPEVE
C
      KWRKE  = KWRK1
      KBVECS = KWRKE + KZYVAR
 9100 FORMAT(/' RSPLAN, work space too small for 3 (Z,Y)-vectors',
     *       /'         had:',I10,', need more than:',I10)
C
      KZRED  = 0
      KZYRED = 0
      THCRSP = THCESR
      IPRRSP = IPRESR
      MAXIT  = MAXESR
C
C     Call RSPCTL to solve linear set of response equations
C
      IF (TRPLET) THEN
         LABEL = 'TRIPLAGR'
      ELSE
         LABEL = 'SINGLAGR'
      END IF
      WRITE (LUPRI,'(//A,I3,/2A)')
     & ' RSPLAN -- linear response calculation for symmetry',KSYMOP,
     & ' RSPLAN -- operator label : ',LABEL
      CALL GETGPV(LABEL,FC,FV,CMO,UDV,PV,XINDX,ANTSYM,WRK(KGP),LWRK1)
      IF (IPRRSP.GT.120) THEN
         WRITE(LUPRI,'(/A)') ' RSPLAN: LAGRANGIAN VECTOR '
         CALL OUTPUT(WRK(KGP),1,KZVAR,1,2,KZVAR,2,1,LUPRI)
      END IF
      XNOR = DDOT(KZYVAR,WRK(KGP),1,WRK(KGP),1)
      IF ( XNOR.LT.DM8 ) THEN
         WRITE(LUPRI,'(/A)')
     *   ' LAGRANGIAN VECTOR IS ZER0, SOLUTION IS SET TO ZERO'
         CALL DZERO(WRK(1),KZYVAR)
         GO TO 9999
      END IF
      WRK(KEIVAL) = D0
      KEXSIM = 1
      KEXCNV = 1
      CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *            .TRUE.,LABEL,BLANK,WRK(KGP),WRK(KREDGP),
     *            WRK(KREDE),WRK(KREDS),
     *            WRK(KIBTYP),WRK(KEIVAL),WRK(KRESID),WRK(KEIVEC),
     *            XINDX,WRK(KWRK1),LWRK1)
C     CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
C    *            LINEQ,GP,REDGP,REDE,REDS,
C    *            IBTYP,EIVAL,EIVEC,XINDX,WRK,LWRK)
C
      NBX = 1
      IBOFF = 0
      CALL RSPEVE(WRK(KIBTYP),WRK(KEIVAL),WRK(KEIVEC),
     *            WRK(KBVECS),WRK(KWRKE),NBX,IBOFF)
C     CALL RSPEVE(IBTYP,EIVAL,EIVEC,BVECS,WRK,NBX,IBOFF)
      IF ( IPRRSP.GE.10 ) THEN
         WRITE (LUPRI,'(/A)') ' SOLUTION VECTOR FOR LAGRANGIAN  '
         CALL RSPPRC(WRK(KBVECS),KZCONF,KZVAR,LUPRI)
         CALL RSPPRO(WRK(KBVECS+KZCONF),KZVAR,LUPRI)
      END IF
      CALL DCOPY(KZYVAR,WRK(KBVECS),1,WRK(1),1)
C
C *** end of RSPLAN --
C
 9999 CALL QEXIT('RSPLAN')
      RETURN
      END
C  /* Deck esrinp */
      SUBROUTINE ESRINP(WORD)
C
C     Module for *ESR    input under **RESPONSE
C
#include "implicit.h"
C
#include "priunit.h"
#include <gnrinf.h>
#include "infrsp.h"
#include "rspprp.h"
#include "inflr.h"
#include "infesr.h"
#include "infpri.h"
#include "zfs.h"
#include "gtensor.h"
#include "maxorb.h"
#include "infinp.h"
#include "parnmr.h"
C
      LOGICAL NEWDEF
      PARAMETER ( NTABLE = 10 )
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER*8 LABEL
      CHARACTER LABDAT(9*ATMNUM)*8
C
      DATA TABLE /'.TRPPRP', '.SNGPRP', '.MAX IT', '.THCESR',
     *            '.PRINT ', '.G-TENS', '.ZFS   ', '.FCCALC',
     *            '.ATOMS ','.SDCALC' /
C
C Initialize new HFC input
      FCFLG  = .FALSE.
      SDFLG  = .FALSE.
      ESRNUC = 0
      DO 50 I = 1, ATMNUM
         DO 55 IH = 1, 7
                ISODAT(IH,I) = 0
 55      CONTINUE
 50   CONTINUE
C      
      NEWDEF = (WORD .EQ. '*ESR   ')
      ICHANG = 0
      ESRCAL = .FALSE.
      GCALC  = .FALSE.
      ZFSCAL = .FALSE.
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') GO TO 100
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' KEYWORD "',WORD,
     *            '" NOT RECOGNIZED IN ESRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT(' ILLEGAL KEYWORD IN LRINP ')
    1          CONTINUE
                  READ(LUCMD,'( A )')LABEL
                  ESROPT( INDPRP(LABEL)) = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ(LUCMD,'( A )')LABEL
                  ESROPS( INDPRP(LABEL)) = .TRUE.
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,*) MAXESR
               GO TO 100
    4          CONTINUE
                  READ (LUCMD,*) THCESR
               GO TO 100
    5          CONTINUE
                  READ (LUCMD,*) IPRESR
               GO TO 100
    6          CONTINUE
                  CALL GINP(WORD) !?
               GO TO 100
    7          CONTINUE
                  ZFSCAL = .TRUE.
               GO TO 100 
    8          CONTINUE
                  FCFLG = .TRUE.
               GO TO 100
    9          CONTINUE
                  READ (LUCMD,*) ESRNUC
                  IF (ESRNUC .GT. ATMNUM) THEN
                     WRITE(LUPRI,'(/A,I4,A,I4/A)')
     &         ' Sorry, but .ATOMS under *ESR is',ESRNUC,
     &           ' which is bigger than ATMNUM in parnmr.h',ATMNUM,
     &         ' Either reduce .ATOMS or increase ATMNUM and recompile.'
                     CALL QUIT('.ATOMS too big, see output')
                  END IF
                  READ (LUCMD,*) (NUCINF(IG), IG = 1, ESRNUC)
               GO TO 100
   10          CONTINUE
                  SDFLG = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/3A/)') ' PROMPT "',WORD,
     *            '" NOT RECOGNIZED IN ESRINP FOR *ESR.'
               CALL QUIT(' ILLEGAL PROMPT under *ESR')
            END IF
         GO TO 100
      END IF
  300 CONTINUE 
      IF (THR_REDFAC .GT. 0.0D0) THEN
         ICHANG = ICHANG + 1
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         THRESR = THRESR*THR_REDFAC
      END IF
C
C FC labels generation
C
      IF (FCFLG) THEN
      DO 400 I = 1, ESRNUC
         CALL FCOPER(NUCINF(I),LABEL)
         ESROPT( INDPRP(LABEL)) = .TRUE.
 400  CONTINUE
      IF (.NOT. SDFLG) CALL SETDEG
      END IF
C
C SD labels generation
C
      IF (SDFLG) THEN
      DO 410 I = 1, ESRNUC
         CALL SDOPER(NUCINF(I),LABDAT,INUM)
         DO 420 IG = 1, INUM
             ESROPT( INDPRP(LABDAT(IG))) = .TRUE.
 420     CONTINUE
 410  CONTINUE
      END IF
C
C  Check degen array
C
      DO 430 I = 1, ESRNUC
         IF (DEGEN(I) .LT. 1.0) DEGEN(I)=1.0
 430  CONTINUE
C
C Setting isotopes data for HFC/pNMR printing
C
      CALL SETISO
C
      WRITE(LUPRI,'(/A)')' ********* ESRINP ********'
C
C
      NTOESR = 0
      DO 500 I = 1,NPRLBL
         IF (ESROPS(I)) NTOESR = NTOESR + 1
         IF (ESROPT(I)) NTOESR = NTOESR + 1
 500  CONTINUE
      ESRCAL = NTOESR .GT. 0 .OR. FCFLG .OR. SDFLG
      IF (ESRCAL .OR. GCALC .OR. ZFSCAL) THEN
         CALL HEADER('Changes of defaults under *ESR   :',0)
         IF ( ESRCAL ) THEN
            WRITE (LUPRI,'(/A)')
     *         ' ESR - hyperfine calculation carried out'
            IF (MAXESR.NE.60) WRITE(LUPRI,'(/A,I6)')
     *         ' MAXIMUM NUMBER OF ITERATIONS. MAXESR =',MAXESR
            IF (THCESR.NE.1.0D-5) WRITE(LUPRI,'(/A,D12.6)')
     *         ' THRESHOLD FOR CONVERGENCE. THCESR =',THCESR
            IF (IPRESR.NE.2) WRITE(LUPRI,'(/A,I5)')
     *         ' PRINT LEVEL IN ESR ROUTINES. IPRESR =',IPRESR
         END IF
         IF (ZFSCAL) THEN
            WRITE (LUPRI,'(/A)')
     *         ' ESR - zero-field splitting  calculation carried out'
         END IF
         IF (GCALC) THEN 
            WRITE (LUPRI,'(/A)')
     *         ' ESR - g-tensor calculation carried out'
            IF (ECC) WRITE(LUPRI,'(A)')
     *         ' Electronic charge centroid used as gauge origin'
            IF (MNFSO) WRITE(LUPRI,'(A)')
     *         ' Mean-field spin-orbit approximation is used'
            IF (SCALED_CHARGES) WRITE(LUPRI,'(A)')
     *         ' Two-electron spin-orbit by scaled nuclear charges'
            IF (.NOT.DOALL) THEN
               WRITE(LUPRI,'(/A)')
     *         ' Selected contributions to electronic g_tensor:'
               IF (DORMC) WRITE(LUPRI,'(/A)') 
     *            '   Relativistic mass correction'
               IF (DOGC1) WRITE(LUPRI,'(/A)') 
     *            '   One-electron gauge correction'
               IF (DOGC2) WRITE(LUPRI,'(/A)') 
     *            '   Two-electron gauge correction'
               IF (DOOZSO1) WRITE(LUPRI,'(/A)') 
     *          '   One-electron spin-orbit + orbital Zeeman correction'
               IF (DOOZSO2) WRITE(LUPRI,'(/A)') 
     *          '   Two-electron spin-orbit + orbital Zeeman correction'
            END IF
         END IF
      ELSE
         WRITE (LUPRI,'(/A)')
     *      ' *ESR   input ignored because no operators requested.'
      END IF
      IF (ISPIN .EQ. 1) THEN
        NWARN = NWARN + 1
         WRITE (LUPRI,'(/A/A)')
     *      ' WARNING: *ESR   input ignored because',
     &      ' singlet reference state has no ESR properties !!!'
        ESRCAL = .FALSE.
        FCFLG  = .FALSE.
        SDFLG  = .FALSE.
        GCALC  = .FALSE.
        ZFSCAL = .FALSE.
      END IF
C
C *** END OF ESRINP
C
      RETURN
      END
C  /* Deck rspesr */
      SUBROUTINE RSPESR(CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
C
C Calculate triplet expectation values with Lagrangian correction
C and singlet expectation values for open shell systems with MCSCF or CI
C (non-zero spin densities).
C
C Revised March 2003 hjaaj
C
#include "implicit.h"
#include "iratdef.h"
C
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
      PARAMETER ( D0 = 0.0D0, DUMMY = 1.0D+20 )
      LOGICAL TRPSAVE
C
#include <maxorb.h>
#include <infinp.h>
C
#include "priunit.h"
#include "infpri.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "rspprp.h"
#include "inflr.h"
#include "infdim.h"
#include "inforb.h"
#include "infesr.h"
C
C
C Common blocks for HFC/pNMR printing
C
#include <codata.h>
#include <gfac.h>
#include <mxcent.h>
#include <nuclei.h>
#include <chrxyz.h>
#include <parnmr.h>
C
C Local variables for HFC values handling
C
      CHARACTER*8 LABEL
      CHARACTER*2 CTMP
C
      CALL QENTER('RSPESR')
      CALL HEADER('Output from RSPESR module',0)
C HFC printing counters
      IFCIND = 0
      REFSPIN=DBLE(ISPIN-1)
C
C      Zero SDVAL 
C
      DO 50 I = 1, ATMNUM
         DO 55 IH = 1, 3
             DO 60 IG = 1, 3
                SDVAL(IG,IH,I) = 0.0D0
 60          CONTINUE
 55      CONTINUE
 50   CONTINUE
C
C DETERMINE SECOND ORDER MOLECULAR PROPERTIES
C
      IF ( NESRT(KSYMOP).GT. 0 ) THEN
         NDREF = KZCONF
C        ... we use determinants when triplet, thus not NCREF
         NDREF = MAX(KZCONF,NCREF)
Chj   ... we use determinants triplet, thus KZCONF,
Chj   ... However, for ROHF NCREF = 1, KZCONF = 0 and
Chj       we need CREF(1) for RSPDM. Thus MAX(KZCONF,NCREF)
C
C        Note that GETGPV uses WRK(KGP) as scratch, thus KGP last allocation !
C
         IF (RSPCI) THEN
            KTUDV = 1
            KCREF = KTUDV + N2ASHX
            KWRK1 = KCREF + NDREF
            KLAGR = KCREF
            KGP   = KCREF
         ELSE
            KLAGR = 1
            CALL RSPLAN(CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
            IF (IPRRSP.GT.120) THEN
               WRITE(LUPRI,'(/A)')
     &            'RSPESR: solution vector for lagrangian :'
               CALL OUTPUT(WRK,1,KZVAR,1,2,KZVAR,2,1,LUPRI)
            END IF
            KTUDV = KLAGR + KZYVAR
            KCREF = KTUDV + N2ASHX
            KGP   = KCREF
            KWRK1 = KGP   + MAX(KZYVAR,NDREF)
         END IF
         LWRK1 = LWRK - KWRK1
         IF (LWRK1.LT.0) CALL ERRWRK('RSPESR',KWRK1-1,LWRK)
C
C        Get one electron spin density (triplet MS=0 density)
C
         CALL GETREF(WRK(KCREF),NDREF)
         KFREE  = 1
         LFREE  = LWRK1
         ISPIN1 = 1
         ISPIN2 = 0
         CALL RSPDM(IREFSY,IREFSY,NDREF,NDREF,WRK(KCREF),WRK(KCREF),
     *              WRK(KTUDV),DUMMY, ISPIN1,ISPIN2,.FALSE.,.TRUE.,
     *              XINDX,WRK(KWRK1),KFREE,LFREE)
C        CALL RSPDM(ISYM,ISYM,NCDIM.NCDIM,REF,REF,
C                   TUDV,DUMMY, ISPIN1,ISPIN2,TDM,NORHO2,
C                   XINDX,WRK,KFREE,LFREE)
         IF (IPRRSP.GE.5) THEN
            WRITE(LUPRI,'(/A)')'RSPESR: one electron spin density '
            CALL OUTPUT(WRK(KTUDV),1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         END IF
         TRPSAVE= TRPLET
         TRPLET = .TRUE.
         DO 100 IOP = 1,NESRT(KSYMOP)
            IF (.NOT. RSPCI ) THEN
               IF (IPRRSP.GT.150) THEN
                  WRITE(LUPRI,'(/A)')
     *              'RSPESR: Lagrangian vector before product'
                  CALL OUTPUT(WRK(KLAGR),1,KZVAR,1,2,KZVAR,2,1,LUPRI)
               END IF
               CALL GETGPV(LBESRT(KSYMOP,IOP),FC,FV,CMO,UDV,PV,XINDX,
     *                     ANTSYM,WRK(KGP),LWRK1 )
C              CALL GETGPV(LABELOP,FC,CMO,UDV,PV,XINDX,ANTSYM,WRK,LWRK)
               IF (IPRRSP.GT.120) THEN
                  WRITE(LUPRI,'(/2A)')
     *              'RSPESR: GP vector with label: ',LBESRT(KSYMOP,IOP)
                  CALL OUTPUT(WRK(KGP),1,KZVAR,1,2,KZVAR,2,1,LUPRI)
               END IF
               FOPLG = -DDOT(KZYVAR,WRK(1),1,WRK(KGP),1)
               IF (IPRRSP.GT.5) WRITE(LUPRI,'(/A,1P,D18.10)')
     &            ' GP * LAGRANG SOLUTION:',FOPLG
C
            ELSE
               FOPLG = D0
            END IF
            CALL PRP1AVE(LBESRT(KSYMOP,IOP),AVEVAL,CMO,
     *                  WRK(KTUDV),WRK(KWRK1),LWRK1,IPRRSP)
            HYPFIN = AVEVAL + FOPLG
            WRITE(LUPRI,'(/2A,3(A,F10.6))')
     *        'TRIPLET OPERATOR: "',LBESRT(KSYMOP,IOP),
     *        '" LAGRANGIAN:',FOPLG,' AVERAGE:',AVEVAL,' TOTAL:',HYPFIN
C
C FC and SD values arrays formation
C
            IF (LBESRT(KSYMOP,IOP)(1:3) .EQ. 'FC ') THEN
                IFCIND=IFCIND+1
                HYPVAL(IFCIND,1)=AVEVAL
                HYPVAL(IFCIND,2)=FOPLG
            END IF
            IF (LBESRT(KSYMOP,IOP)(1:3) .EQ. 'SD ') THEN
                LABEL=LBESRT(KSYMOP,IOP)
                READ (LABEL,'(3X,I3)') INDSD1
                CTMP=LABEL(7:8)
                indsd2 = -9898
                IF (CTMP .EQ. ' x') INDSD2=1
                IF (CTMP .EQ. ' y') INDSD2=2
                IF (CTMP .EQ. ' z') INDSD2=3
                if (indsd2 .eq. -9898) then
                   write(lupri,*) 'hjaaj label=',label, indsd1
                   write(lupri,*) 'hjaaj no x or y or z in label(7:8) !'
                   call quit('hjaaj: problem in lagrang.F')
                end if
                SDVAL(SDIND(2,INDSD1),INDSD2,SDIND(1,INDSD1))=HYPFIN
            END IF
C
100      CONTINUE
         TRPLET = TRPSAVE
      END IF
      DO 300 ISYM = 2,NSYM
         DO 350 IOP = 1,NESRT(ISYM)
            WRITE(LUPRI,'(/3A/A,I3)') 'TRIPLET OPERATOR: "',
     &         LBESRT(ISYM,IOP),'" contribution = 0.0 by symmetry.',
     &         '- Operator is of symmetry no.',ISYM
350      CONTINUE
300   CONTINUE
C
C     Singlet operators:
C
      TRPSAVE= TRPLET
      TRPLET = .FALSE.
C     ... so PRP1AVE calculates singlet expectation values /hjaaj march 2003 ...
C         (if TRPLET true, then inactive density matrix is omitted!)
      DO 400 IOP = 1,NESRS(KSYMOP)
         CALL PRP1AVE(LBESRS(KSYMOP,IOP),AVEVAL,CMO,
     *               UDV,WRK(KWRK1),LWRK1,IPRRSP)
C        CALL PRP1AVE(LABELOP,AVEVAL,CMO, UDV,WRK,LWRK,IPRINT)
         WRITE(LUPRI,'(/3A,F10.6)')
     *      'SINGLET OPERATOR: "',LBESRS(KSYMOP,IOP),'" AVERAGE:',AVEVAL
400   CONTINUE
C
      DO 500 ISYM = 2,NSYM
         DO 550 IOP = 1,NESRS(ISYM)
            WRITE(LUPRI,'(/3A/A,I3)') 'SINGLET OPERATOR: "',
     &         LBESRS(ISYM,IOP),'" contribution = 0.0 by symmetry.',
     &         '- Operator is of symmetry no.',ISYM
550      CONTINUE
500   CONTINUE
      TRPLET = TRPSAVE
C
C Isotropic hyperfine coupling printing
C
      IF (FCFLG) THEN
          CALL AROUND('Isotropic  Hyperfine Coupling')
          WRITE(LUPRI,'(A)') '    Atom      Mass      G-val   '
     *                    //'   A, Mhz        A, G    '
          WRITE(LUPRI,'(A)') ' -----------------------------'
     *                    //'--------------------------'
          DO 600 ITMP = 1, ESRNUC
              DO 610 ISONM = 1, ISODAT(1,NUCINF(ITMP))
                  IATIS=ISODAT(ISONM+1,NUCINF(ITMP))
                  XATGV=DISOTP(INT(CHARGE(NUCINF(ITMP))),IATIS,'GVAL')
                  XISMAS=DISOTP(INT(CHARGE(NUCINF(ITMP))),IATIS,'MASS')
                  HYPFIN=((HYPVAL(ITMP,1)+HYPVAL(ITMP,2))
     *                   *XATGV*XTHZ*1.0D-6*ALPHA2)
     *                   /(2.0D0*XPRTMAS*DEGEN(ITMP)*REFSPIN)
                  WRITE(LUPRI,'(A,F6.2,A,F7.5,A,F9.4,A,F9.4,A)')
     *               ' *   '//NAMN(NUCINF(ITMP))(1:4)//'  * ', XISMAS,
     *               ' * ', XATGV,' * ',HYPFIN,'  * ',HYPFIN/2.8025D0,
     *               ' *'
 610          CONTINUE
 600      CONTINUE
          WRITE(LUPRI,'(A)') ' -----------------------------'
     *                      //'--------------------------'
          WRITE(LUPRI,'(/A)') ' *** NOTE: Results printed only for'
     *                         //' symmetry unique atoms!'
          CALL AROUND('R-U contributions to Isotropic hyperfine'
     *                //' coupling')
          WRITE(LUPRI,'(A)') '    Atom    G-val    Average, MHz '
     *                    //' Response, MHz   Total, MHz '
          WRITE(LUPRI,'(A)') ' ------------------------------'
     *                    //'---------------------------------'
          DO 601 ITMP = 1, ESRNUC
              DO 611 ISONM = 1, ISODAT(1,NUCINF(ITMP))
                  IATIS=ISODAT(ISONM+1,NUCINF(ITMP))
                  XATGV=DISOTP(INT(CHARGE(NUCINF(ITMP))),IATIS,'GVAL')
                  HYPFIN=(XATGV*XTHZ*1.0D-6*ALPHA2)
     *                   /(2.0D0*XPRTMAS*DEGEN(ITMP)*REFSPIN)
                  WRITE(LUPRI,'(A,F7.5,A,F11.4,A,F11.4,A,F11.4,A)')
     *               ' *  '//NAMN(NUCINF(ITMP))(1:4)//' * ', XATGV,
     *               ' * ',HYPFIN*HYPVAL(ITMP,1),' * ',
     *               HYPFIN*HYPVAL(ITMP,2),'  * ',
     *               HYPFIN*(HYPVAL(ITMP,1)+HYPVAL(ITMP,2)),'  * '
 611          CONTINUE
 601      CONTINUE
          WRITE(LUPRI,'(A)') ' ------------------------------'
     *                    //'---------------------------------'
          WRITE(LUPRI,'(/A)') ' *** NOTE: Results printed only for'
     *                         //' symmetry unique atoms!'
      ENDIF
C
C Anisotropic hyperfine coupling printing
C
      IF (SDFLG) THEN
            CALL AROUND('Anisotropic Hyperfine Coupling')
      DO 700 ITMP=1,ESRNUC
         DO 710 ISONM=1,ISODAT(1,NUCINF(ITMP))
            IATIS=ISODAT(ISONM+1,NUCINF(ITMP))
            XATGV=DISOTP(INT(CHARGE(NUCINF(ITMP))),IATIS,
     *      'GVAL')
            XISMAS=DISOTP(INT(CHARGE(NUCINF(ITMP))),IATIS,
     *      'MASS')
            HYPFIN=(XATGV*XTHZ*1.0D-6*ALPHA2)/
     *             (2.0D0*XPRTMAS*DEGEN(ITMP)*REFSPIN)
            WRITE(LUPRI,'(/A,F6.2,A,F7.5,A)')
     *                  ' A Tensor components for '//
     *                  NAMN(NUCINF(ITMP))(1:4)//' ( Mass =',
     *                  XISMAS, ' G-val =', XATGV,') in MHz: '
            WRITE(LUPRI,'(A)') ' ================================='
     *                       //'=================================='
            WRITE(LUPRI,'(/A)') '  '
            WRITE(LUPRI,'(A)') '           SX          SY          SZ '
     *                        //'     '
            DO 720 ICRD=1,3
               WRITE(LUPRI,'(A,F11.4,A,F11.4,A,F11.4,A)')
     *         ' I'//CHRXYZ(ICRD)//' *  ',
     *         HYPFIN*SDVAL(ICRD,1,NUCINF(ITMP)), ' * ',
     *         HYPFIN*SDVAL(ICRD,2,NUCINF(ITMP)), ' * ',
     *         HYPFIN*SDVAL(ICRD,3,NUCINF(ITMP)), ' * '
 720        CONTINUE
 710     CONTINUE
 700  CONTINUE
      WRITE(LUPRI,'(/A)') ' *** NOTE: Results printed only for'
     *                         //' symmetry unique atoms!'
      ENDIF
C
C *** end of RSPESR --
      CALL QEXIT('RSPESR')
      RETURN
      END
C  /* Deck solgdt */
      SUBROUTINE SOLGDT(CREF,CMO,INDXCI,DV,G,ESOLT,ERLM,WRK,
     &                  LFREE,NHCREF,KREFSY)
C
C   Based on SOLGRD:
C   Copyright 29-Nov-1986 Hans Joergen Aa. Jensen
C
C   Purpose:  calculate MCSCF energy and gradient contribution
C             from a surrounding medium, cavity radius = Rsol
C             and dielectric constant = EPsol.
C
C   Output:
C    G          MCSCF gradient with solvation contribution added
C    ESOLT      total solvation energy
C    ERLM(lm,1) contains Esol(l,m) contribution to ESOLT
C    ERLM(lm,2) contains Tsol(l,m)
C
#include "implicit.h"
#include "dummy.h"
C
      DIMENSION CREF(*), CMO(*), INDXCI(*)
      DIMENSION DV(*),   G(*),   ERLM(NLMSOL,2),  WRK(*)
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0 )
#include "thrzer.h"
#include "iratdef.h"
#include "priunit.h"
#include "infrsp.h"
C
C Used from common blocks:
C   INFINP: NLMSOL, LSOLMX, EPSOL, RSOL(3)
C   INFVAR: NCONF,  NWOPT,  NVAR,   NVARH
C   INFORB: NNASHX, NNBASX, NNORBX, etc.
C   INFIND: IROW(*)
C   INFTAP: LUSOL,  LUIT2
C   INFPRI: IPRSOL
C
#include "maxash.h"
#include "maxorb.h"
#include "infinp.h"
#include "infvar.h"
#include "inforb.h"
#include "infind.h"
#include "inftap.h"
#include "infpri.h"
C
      LOGICAL     FIRST
      PARAMETER  (MXLMAX = 50)
      DIMENSION   ISYRLM(2*MXLMAX+1)
      CHARACTER*8 STAR8, SOLVDI, EODATA
      SAVE        FIRST
      DATA        FIRST/.TRUE./, STAR8/'********'/
      DATA        SOLVDI/'SOLVDIAG'/, EODATA/'EODATA  '/
C
C     Statement functions;
C     define automatic arrays (dynamic core allocation)
C
      FLVEC(LM) = WRK(LM)
      FLINR(LM) = WRK(KFLINR-1+LM)
      TLMSI(LM) = WRK(KTLMSI-1+LM)
C
      CALL QENTER('SOLGDT')
C
      IF (LSOLMX .GT. MXLMAX) THEN
         WRITE (LUERR,*) 'ERROR SOLGDT, increase MXLMAX parameter'
         WRITE (LUERR,*) ' LSOLMX =',LSOLMX
         WRITE (LUERR,*) ' MXLMAX =',MXLMAX
         CALL QUIT('ERROR SOLGDT, increase MXLMAX parameter')
      END IF
C
C     Core allocation
C        FLVEC  f(l) factors in solvent energy expression
C        DIASH  diagonal of solvent contribution to Hessian
C        GRDLM  TELM gradient for current l,m value in the l,m loop
C        UCMO   CMO unpacked (i.e. no symmetry blocking)
C        RLMAC  active-active subblock of RLM
C        RLM    R(l,m) integrals for current l,m value in l,m loop
C
C     If (INERSF)
C       (i.e. If (inertial polarization contribution to final state))
C        FLINR  f(l) factors for inertial pol. contrib.
C        TLMSI  T(lm) values for initial state
C     end if
C
      KFLVEC = 1
C     ... NOTE: KFLVEC = 1 assumed in FLVEC(LM) definition above.
      IF (INERSF) THEN
         KFLINR = KFLVEC + NLMSOL
         KTLMSI = KFLINR + NLMSOL
         KDIASH = KTLMSI + NLMSOL
      ELSE
         KFLINR = 1
         KTLMSI = 1
         KDIASH = KFLVEC + NLMSOL
      END IF
      KGRDLM = KDIASH + NVAR
      KUCMO  = KGRDLM + NVARH
      KRLMAC = KUCMO  + NORBT*NBAST
      KRLM   = KRLMAC + NNASHX
      KW10   = KRLM   + NNORBX
C     1.1 read rlmao in ao basis and transform to rlm in mo basis
      KRLMAO = KW10
      KW20   = KRLMAO + NNBASX
C     1.2 diagonal contribution for current l,m value in the l,m loop
      KDIALM = KW10
      KW21   = KDIALM + NVAR
      LW21   = LFREE  - KW21
C     1.3 rest of CSF contribution
      KW22   = KW10
C
      KTDV  = MAX(KW20,KW21,KW22)
      KWRK1  = KTDV + NASHT * NASHT
      LWRK1  = LFREE  - KWRK1
      IF (LWRK1 .LT. 0) CALL ERRWRK('SOLGDT',-KWRK1,LFREE)
C
      IF (IPRSOL .GE. 130) THEN
         WRITE (LUPRI,'(/A/A,2I10)')
     *        ' >>> SOLGDT - gtot (input) - non-zero elements',
     *        '     NCONF, NWOPT =',NCONF,NWOPT
         DO 40 I = 1,NCONF
            IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *           ' conf #',I,G(I)
 40      CONTINUE
         DO 50 I = NCONF+1,NVAR
            IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *           ' orb  #',I,G(I)
 50      CONTINUE
      END IF
C
C     Calculate f(l) factors
C     If (INERSF) FLVEC factors describe the optical polarization
C             and FLINR factors describe the inertial polarization
C     else        FLVEC may describe optical or static polarization
C
      CALL SOLFL(WRK(KFLVEC),EPSOL,RSOL,LSOLMX)
      IF (INERSF) THEN
         CALL SOLINR(WRK(KFLVEC),WRK(KFLINR),WRK(KTLMSI))
      END IF
      IF ((IPRSOL .GE. 5 .AND. FIRST) .OR. IPRSOL .GE. 15) THEN
         IF (.NOT.INERSF) THEN
            WRITE (LUPRI,'(//A/A)')
     *      ' >>> SOLGDT:  l     f(l) factor',
     *      '             === ================='
         ELSE
            WRITE (LUPRI,'(//A/A)')
     *      ' >>> SOLGDT:  l  optical f(l) factor inertial f(l) factor',
     *      '             === =================== ===================='
         END IF
         DO 140 L = 0,LSOLMX
            LL = (L+1)*(L+1)
            FL = FLVEC(LL)
            IF (INERSF) THEN
               FLI = FLINR(LL)
               WRITE (LUPRI,'(I15,F17.10,F21.10)') L, FL, FLI
            ELSE
               WRITE (LUPRI,'(I15,F16.10)') L, FL
            END IF
  140    CONTINUE
      END IF
C
C     Read and check dimension information (if first read) and
C     nuclear contributions to ERLM (always).
C
      CALL GPOPEN(LUSOL,FNSOL,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUSOL
      CALL MOLLAB('SOLVRLM ',LUSOL,LUERR)
      IF (FIRST) THEN
         READ (LUSOL) LMAXSS, LMTOT, NNNBAS
         NERR = 0
         IF (LMAXSS .LT. LSOLMX) THEN
            NERR = NERR + 1
            WRITE (LUPRI,'(//2A,2(/A,I5))') ' >>> SOLGDT ERROR,',
     *      ' insufficient number of integrals on LUSOL',
     *      ' l max from SIRIUS input :',LSOLMX,
     *      ' l max from LUSOL  file  :',LMAXSS
         END IF
         IF ((LMAXSS+1)**2 .NE. LMTOT) THEN
            NERR = NERR + 1
            WRITE (LUPRI,'(//2A,3(/A,I5))') ' >>> SOLGDT ERROR,',
     *      ' LUSOL file info inconsistent',
     *      ' l_max               :',LMAXSS,
     *      ' (l_max + 1) ** 2    :',(LMAXSS+1)**2,
     *      ' LMTOT               :',LMTOT
         END IF
         IF (NNNBAS .NE. NBAST) THEN
            NERR = NERR + 1
            WRITE (LUPRI,'(//2A,3(/A,I5))') ' >>> SOLGDT ERROR,',
     *      ' LUSOL file info inconsistent with SIRIUS input',
     *      ' NBAST - LUSOL       :',NNNBAS,
     *      ' NBAST - SIRIUS      :',NBAST
         END IF
         IF (NERR .GT. 0) THEN
            CALL QUIT('SOLGDT ERROR: LUSOL file not OK for this calc.')
         END IF
      ELSE
         READ (LUSOL)
      END IF
      CALL READT(LUSOL,NLMSOL,ERLM(1,2))
C
      IF (IPRSOL .GE. 20 .AND. NASHT .GT. 0) THEN
         WRITE (LUPRI,'(/A)') ' >>> SOLGDT - DV matrix :'
         CALL OUTPAK(DV,NASHT,1,LUPRI)
      END IF
      IF (IPRSOL .GE. 7) THEN
         WRITE (LUPRI,'(/A/)')
     *      ' l, m, Tn(lm) - the nuclear contributions :'
         LM = 0
         DO 220 L = 0,LSOLMX
            DO 210 M = -L,L
               LM = LM + 1
               WRITE (LUPRI,'(2I5,F15.10)') L,M,ERLM(LM,2)
  210       CONTINUE
            WRITE (LUPRI,'()')
  220    CONTINUE
      END IF
C
C     Unpack symmetry blocked CMO
C     Loop over l,m expansion
C
      CALL UPKCMO(CMO,WRK(KUCMO))
      IF (IPRSOL .GE. 6)
     *   WRITE (LUPRI, '(//A/)') ' >>> SOLGDT: START LOOP OVER LM'
      CALL DZERO(WRK(KDIASH),NVAR)
      LM = 0
      DO 520 L = 0,LSOLMX
         READ (LUSOL) L1,(ISYRLM(M),M=1,2*L+1)
         IF (L1 .NE. L) THEN
            WRITE (LUERR,*) 'ERROR SOLGDT: L from LUSOL not as expected'
            WRITE (LUERR,*) 'L from 520 loop:',L
            WRITE (LUERR,*) 'L from LUSOL   :',L1
            CALL QUIT('ERROR SOLGDT: L from LUSOL not as expected')
         END IF
      DO 500 M = -L,L
         LM = LM + 1
         IF (IPRSOL .GE. 15) THEN
            WRITE (LUPRI,'(/A,2I5/A)') ' >>> l, m :',L,M,
     *                                ' ===================='
            WRITE (LUPRI,'(A,I2)') ' Symmetry :',ABS(ISYRLM(L+M+1))
         END IF
         IF (ISYRLM(L+M+1) .NE. 1) THEN
            IF (ABS(ERLM(LM,2)) .GT. 1000.D0*THRZER) THEN
               WRITE (LUPRI,*) 'ERROR SOLGDT for l,m',L,M
               WRITE (LUPRI,*) 'Symmetry :',ISYRLM(L+M+1)
               WRITE (LUPRI,*) 'Tn(l,m) .ne. 0, but =',ERLM(LM,2)
               CALL QUIT('ERROR SOLGDT: Tn(l,m) not 0 as expected')
            END IF
            ERLM(LM,2) = D0
C           ... to fix round-off errors in Tn(l,m) calculation
            IF (ISYRLM(L+M+1) .GT. 1) READ (LUSOL)
            GO TO 500
         END IF
C
C        Read R(l,m) in ao basis and transform to mo basis.
C        Extract active-active block in RLMAC(1) = WRK(KRLMAC).
C
         CALL READT(LUSOL,NNBASX,WRK(KRLMAO))
         CALL UTHU(WRK(KRLMAO),WRK(KRLM),WRK(KUCMO),WRK(KWRK1),
     &             NBAST,NORBT)
         IF (NASHT .GT. 0) THEN
            CALL GETAC2(WRK(KRLM),WRK(KRLMAC))
         END IF
         IF (IPRSOL .GE. 15) THEN
            WRITE (LUPRI,'(/A)') ' Rlm_ao matrix:'
            CALL OUTPAK(WRK(KRLMAO),NBAST,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' Rlm_mo matrix:'
            CALL OUTPAK(WRK(KRLM),  NORBT,1,LUPRI)
            IF (NASHT .GT. 0) THEN
               WRITE (LUPRI,'(/A)') ' Rlm_ac matrix:'
               CALL OUTPAK(WRK(KRLMAC),NASHT,1,LUPRI)
            END IF
         END IF
C
C        Add electronic contribution TE(l,m) to T(l,m)
C
         KFREE=1
         ISPIN1=0
         ISPIN2=0
         CALL RSPDM(KREFSY,KREFSY,NHCREF,NHCREF,CREF,CREF,
     *      WRK(KTDV),DUMMY, ISPIN1,ISPIN2,.FALSE.,.TRUE.,
     *      INDXCI,WRK(KWRK1),KFREE,LWRK1)
C
C  TRIANGULAR PACKING OF ONE ELECTRON DENSITY MATRIX
C
         CALL DSITSP(NASHT,WRK(KTDV),DV)
C
         TELM     = SOLELM(DV,WRK(KRLMAC),WRK(KRLM),TELMAC)
C
C     construct again triplet density matrix
         KFREE =1
         ISPIN1=1
         ISPIN2=0
         CALL RSPDM(KREFSY,KREFSY,NHCREF,NHCREF,CREF,CREF,
     *      WRK(KTDV),DUMMY, ISPIN1,ISPIN2,.FALSE.,.TRUE.,
     *      INDXCI,WRK(KWRK1),KFREE,LWRK1)
C
C  TRIANGULAR PACKING OF ONE ELECTRON DENSITY MATRIX
C
         CALL DSITSP(NASHT,WRK(KTDV),DV)
C
         IF (IPRSOL .GE. 6) THEN
            WRITE (LUPRI,'(A,2I5,/A,3F17.8)')
     *      ' >>> l, m :',L,M,
     *      '     Te(lm), Tn(lm), T(lm) :',
     *         TELM,ERLM(LM,2),ERLM(LM,2)-TELM
            IF (IPRSOL .GE. 10) WRITE (LUPRI,'(A,F17.8)')
     *      ' --- active part of Te(lm) :',TELMAC
            IF (INERSF) WRITE (LUPRI,'(A,F17.8)')
     *      ' --- inertial T(lm) value  :',TLMSI(LM)
         END IF
         ERLM(LM,2) = ERLM(LM,2) - TELM
      IF (ABS(ERLM(LM,2)) .LE. THRZER) THEN
         ERLM(LM,2) = D0
         GO TO 500
      END IF
C     ... test introduced 880109 hjaaj
C         (the only possible problem is the DO 420 loop,
C          but I think (w.o. having checked) that this
C          contribution to the Hessian diagonal also will be
C          zero if ERLM(LM,2) zero).
C
C        Calculate orbital TE(l,m) gradient contribution
C        and part of csf contribution.
C
         CALL DZERO(WRK(KGRDLM),NVARH)
         IF (NCONF .GT. 1) THEN
            CALL SOLGC(CREF,WRK(KRLMAC),TELMAC,WRK(KGRDLM),INDXCI,
     &                 WRK(KWRK1),LWRK1)
C           CALL SOLGC(CREF,RLMAC,TELMAC,GLMCI,INDXCI,WRK,LWRK)
         END IF
         IF (NWOPT .GT. 0) THEN
            CALL SOLGO(D0,DV,WRK(KRLM),WRK(KGRDLM+NCONF))
         END IF
C
C
C        Obtain DIALM = diagonal TE(l,m) Hessian
C                     = 2 ( <i|R(l,m)|i> - TE(l,m) )
C        Add the DIALM contribution and the GRDLM contribution
C        to solvent Hessian diagonal.
C
Clf the diagonal hessian is not needed for triplet gradients
Clf         CALL SOLDIA(TELMAC,WRK(KRLMAC),INDXCI,
Clf     *               WRK(KRLM),DV,WRK(KDIALM),WRK(KW21),LW21)
C        CALL SOLDIA(TELM,RLMAC,INDXCI,RLM,DV,DIAG,WRK,LFREE)
C
         FAC1 = - D2 * FLVEC(LM) * ERLM(LM,2)
         IF (INERSF) THEN
            FAC1 = FAC1 - FLINR(LM) * D2 * TLMSI(LM)
         END IF
         FAC2 =   D2 * FLVEC(LM)
         DO 420 I = 0,(NVAR-1)
            WRK(KDIASH+I) = WRK(KDIASH+I)
     *                    + FAC1 * WRK(KDIALM+I)
     *                    + FAC2 * WRK(KGRDLM+I) * WRK(KGRDLM+I)
  420    CONTINUE
C
C        test orthogonality
C
         IF (IPRSOL .GE. 120) THEN
           WRITE (LUPRI,'(/A)')' >>> SOLGDT - grdlm, dialm, diash, cref'
           DO 430 I = 1,NCONF
              WRITE (LUPRI,'(A,I10,4F12.6)') ' conf #',I,
     *        WRK(KGRDLM-1+I),WRK(KDIALM-1+I),WRK(KDIASH-1+I),CREF(I)
  430      CONTINUE
         END IF
         TEST = DDOT(NCONF,CREF,1,WRK(KGRDLM),1)
         IF (ABS(TEST) .GT. 1.D-8) THEN
            NWARN = NWARN + 1
            WRITE (LUPRI,'(/A,I5,/A,1P,D12.4)')
     *      ' >>> SOLGDT WARNING, for LM =',LM,
     *      ' <CREF | GRADlm > =',TEST
         END IF
C
C        Add TE(l,m) gradient contribution to MCSCF gradient
C        g  =  g  -  2 f(l) * T(l,m) * (dTE(l,m)/d(lambda))
C
         FAC = - D2 * FLVEC(LM) * ERLM(LM,2)
         IF (INERSF) THEN
            FAC = FAC - FLINR(LM) * D2 * TLMSI(LM)
         END IF
         CALL DAXPY(NVARH,FAC,WRK(KGRDLM),1,G,1)
         IF (IPRSOL .GE. 140) THEN
            WRITE (LUPRI,'(/A/A,2I10)')
     *         ' >>> SOLGDT - grdlm, gtot (accum) - non-zero grdlm',
     *         '     NCONF, NWOPT =',NCONF,NWOPT
            DO 440 I = 1,NCONF
               IF (WRK(KGRDLM-1+I) .NE. D0)
     *            WRITE (LUPRI,'(A,I10,3F15.10)')
     *            ' conf #',I,FAC*WRK(KGRDLM-1+I),G(I)
  440       CONTINUE
            DO 450 I = NCONF+1,NVAR
               IF (WRK(KGRDLM-1+I) .NE. D0)
     *            WRITE (LUPRI,'(A,I10,3F15.10)')
     *            ' orb  #',I,FAC*WRK(KGRDLM-1+I),G(I)
  450       CONTINUE
         END IF
C
  500 CONTINUE
  520 CONTINUE
C
      CALL GPCLOSE(LUSOL,'KEEP')
C
C     500 is end of (l,m) loop.
C
C
         IF (IPRSOL .GE. 130) THEN
            WRITE (LUPRI,'(/A/A,2I10)')
     *         ' >>> SOLGDT - gtot (output) - non-zero elements',
     *         '     NCONF, NWOPT =',NCONF,NWOPT
            DO 840 I = 1,NCONF
               IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *         ' conf #',I,G(I)
  840       CONTINUE
            DO 850 I = NCONF+1,NVAR
               IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *         ' orb  #',I,G(I)
  850       CONTINUE
         END IF
C
C
C     Calculate ER(l,m) energy contributions and add them up
C
      ESOLT = D0
      DO 900 LM = 1,NLMSOL
         ERLM(LM,1) = FLVEC(LM) * ERLM(LM,2) * ERLM(LM,2)
         IF (INERSF) THEN
            ERLM(LM,1) = ERLM(LM,1)
     *                 + FLINR(LM) * ERLM(LM,2) * D2 * TLMSI(LM)
     *                 - FLINR(LM) * TLMSI(LM) * TLMSI(LM)
         END IF
         ESOLT    = ESOLT     + ERLM(LM,1)
  900 CONTINUE
C
      FIRST = .FALSE.
      CALL QEXIT('SOLGDT')
      RETURN
C     end of solgdt.
      END
C  /* Deck FCOPER */
      SUBROUTINE FCOPER(ATMIND,LABINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      CHARACTER*8 LABINT
      INTEGER ATMIND
#include <nuclei.h>
#include <chrnos.h>
C
       LABINT = 'FC '//NAMN(ATMIND)(1:3)//CHRNOS(ATMIND/10)
     *          //CHRNOS(MOD(ATMIND,10))
C
      RETURN
      END
C  /* Deck SETDEG */
      SUBROUTINE SETDEG
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <maxorb.h>
#include <parnmr.h>
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <chrxyz.h>
#include <ibtfun.h>
C
      DO 50 I = 1, ESRNUC
          DEGEN(I) = 0.0D0
          DO 100 IREP=0,MAXREP
                ISCOR1 = IPTCNT(3*(NUCINF(I) - 1) + 1,IREP,2)
                IF (ISCOR1 .GT. 0)  DEGEN(I)=DEGEN(I)+1.0D0
 100      CONTINUE
 50   CONTINUE
C
      RETURN
      END
C  /* Deck SETISO */
      SUBROUTINE SETISO
#include <implicit.h>
#include <priunit.h>
#include <parnmr.h>
#include <mxcent.h>
#include <nuclei.h>
C
      DO 100 I=1,ESRNUC
         IF (ISODAT(1,NUCINF(I)) .EQ. 0) THEN
             ISODAT(1,NUCINF(I))=1
             IATIS=1
 200         XATGV=DISOTP(INT(CHARGE(NUCINF(I))),IATIS,'GVAL')
             IF (DABS(XATGV) .LT. 1.0D-5) THEN
                 IATIS=IATIS+1
             GO TO 200
             ENDIF
             ISODAT(2,NUCINF(I))=IATIS
         ENDIF
 100  CONTINUE
C
      RETURN
      END
C  /* Deck SDOPER */
      SUBROUTINE SDOPER(ATMIND,LABINT,NCOOR)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <maxorb.h>
#include <parnmr.h>
      INTEGER ATMIND, NCOOR
      CHARACTER LABINT(9*ATMNUM)*8
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <chrxyz.h>
#include <ibtfun.h>
C
      NCOOR=0
      DEGEN(ATMIND) = 0.0D0
      DO 100 IREP=0,MAXREP
         DO 200 ICOOR1=1,3
            ISCOR1 = IPTCNT(3*(NUCINF(ATMIND ) - 1) + ICOOR1,IREP,2)
            IF (ISCOR1 .GT. 0) THEN
              IF (ICOOR1 .EQ. 1) DEGEN(ATMIND)=DEGEN(ATMIND)+1.0D0
              DO 300 ICOOR2 = 1, 3
                 NCOOR=NCOOR+1
                 LABINT(NCOOR) = 'SD '//CHRNOS(ISCOR1/100)
     *                            //CHRNOS(ISCOR1/10)
     *                            //CHRNOS(MOD(ISCOR1,10))
     *                            //' '//CHRXYZ(-ICOOR2)
 300          CONTINUE
              SDIND(1,ISCOR1)=ATMIND
              SDIND(2,ISCOR1)=ICOOR1
            ENDIF
 200     CONTINUE
 100  CONTINUE
C
      RETURN
      END


C  /* Deck pcmgdt */
      SUBROUTINE PCMGDT(CREF,CMO,INDXCI,DV,G,ESOLT,WRK,LFREE,
     $                  NHCREF,KREFSY)
C
C   Based on SOLGDT:
C   13-02-2006 Luca Frediani
C
C   Purpose:  calculate MCSCF energy and gradient contribution
C             from a surrounding medium with PCM
C
C   Output:
C    G          MCSCF gradient with solvation contribution added
C
#include "implicit.h"
#include "dummy.h"
C
      DIMENSION CREF(*), CMO(*), INDXCI(*)
      DIMENSION DV(*),   G(*),   WRK(*)
      PARAMETER ( D1=1.0d0, DM1=-1.0d0, D0 = 0.0D0, D2 = 2.0D0 )
#include "thrzer.h"
#include "iratdef.h"
#include "priunit.h"
#include "infrsp.h"
#include "mxcent.h"
#include "orgcom.h"
#include "pcmdef.h"
#include "pcmlog.h"
#include "pcm.h"
C
C Used from common blocks:
C   INFINP: NLMSOL, LSOLMX, EPSOL, RSOL(3)
C   INFVAR: NCONF,  NWOPT,  NVAR,   NVARH
C   INFORB: NNASHX, NNBASX, NNORBX, etc.
C   INFIND: IROW(*)
C   INFTAP: LUSOL,  LUIT2
C   INFPRI: IPRSOL
C
#include "maxash.h"
#include "maxorb.h"
#include "infinp.h"
#include "infvar.h"
#include "inforb.h"
#include "infind.h"
#include "inftap.h"
#include "infpri.h"
C
      LOGICAL     FNDLAB, EXP1VL, TOFILE, TRIMAT
      PARAMETER  (MXLMAX = 50)
      CHARACTER*8 STAR8, SOLVDI, EODATA, LABINT(9*MXCENT)
      DATA        STAR8/'********'/
      DATA        SOLVDI/'SOLVDIAG'/, EODATA/'EODATA  '/
C
C     Statement functions;
C     define automatic arrays (dynamic core allocation)
C
      CALL QENTER('PCMGDT')
C
C
C     Core allocation
C        DIASH  diagonal of solvent contribution to Hessian
C        GRDLM  TELM gradient for current l,m value in the l,m loop
C        UCMO   CMO unpacked (i.e. no symmetry blocking)
C
      KJENAO = 1
      KJEN   = KJENAO + NNBASX
      KJ1AO  = KJEN   + NNORBX
      KJ1    = KJ1AO  + NNBASX*NSYM
      KJENAC = KJ1    + NNORBX
      KJ1AC  = KJENAC + NNASHX
      KDIASH  = KJ1AC  + NNASHX
      KGRDLM = KDIASH + NVAR
      KUCMO  = KGRDLM + NVARH
      KJ2GRD = KUCMO  + NORBT*NBAST
      KPOT   = KJ2GRD + NVAR
      KDENC  = KPOT   + NTS
      KDENV  = KDENC  + N2BASX
      KW10   = KDENV  + N2BASX
C     1.3 rest of CSF contribution
      KDIALM = KW10
      KTDV   = KDIALM + NVAR
      KW20   = KTDV + NASHT * NASHT
C
C     Allocations for non-equlibrium energy solvation
      IF (NONEQ) THEN
         KMPOT  = KW20
         KQSEGR = KMPOT  + NTS * NTS 
         KPOTGR = KQSEGR + NTS
         KW21   = KPOTGR + NTS
      ELSE
         KW21   = KW20
      END IF
C
      LW21   = LFREE - KW21
C
      KWRK1  = KW21
      LWRK1  = LFREE  - KWRK1
      IF (LWRK1 .LT. 0) CALL ERRWRK('PCMGDT',-KWRK1,LFREE)
C

      KFREE=1
      ISPIN1=0
      ISPIN2=0
      CALL RSPDM(KREFSY,KREFSY,NHCREF,NHCREF,CREF,CREF,
     *     WRK(KTDV),DUMMY, ISPIN1,ISPIN2,.FALSE.,.TRUE.,
     *     INDXCI,WRK(KWRK1),KFREE,LWRK1)
C     
C     TRIANGULAR PACKING OF ONE ELECTRON DENSITY MATRIX
C     
      CALL DSITSP(NASHT,WRK(KTDV),DV)
C     
clf      IF (IPRPCM .GE. 130) THEN
      IF (.true.) THEN
         WRITE (LUPRI,'(/A/A,2I10)')
     *        ' >>> PCMGDT - gtot (input) - non-zero elements',
     *        '     NCONF, NWOPT =',NCONF,NWOPT
         DO 40 I = 1,NCONF
            IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *           ' conf #',I,G(I)
 40      CONTINUE
         DO 50 I = NCONF+1,NVAR
            IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *           ' orb  #',I,G(I)
 50      CONTINUE
      END IF

      IF (IPRSOL .GE. 20 .AND. NASHT .GT. 0) THEN
         WRITE (LUPRI,'(/A)') ' >>> SOLGDT - DV matrix :'
         CALL OUTPAK(DV,NASHT,1,LUPRI)
      END IF
C
C     Unpack symmetry blocked CMO
C     Loop over l,m expansion
C
      CALL UPKCMO(CMO,WRK(KUCMO))
      CALL DZERO(WRK(KDIASH),NVAR)
      CALL DZERO(WRK(KDIALM),NVAR)
      CALL DZERO(WRK(KJEN),NNORBX)
C
C     Read JEN = J(en) + J(ne) in ao basis and transform to mo basis.
C     Extract active-active block in WRK(KJENAC).
C     
      LUPBKP = LUPROP
      IF (LUPROP .LT. 0) CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',
     &     'SEQUENTIAL','UNFORMATTED',IDUMMY,.FALSE.) 
CLF      IF (.NOT. (FNDLAB('NE-PCMIN',LUPROP))) THEN
      IF (.TRUE.) THEN
         CALL PCMJMAT(WRK(KJENAO),NNBASX,WRK(KWRK1),LWRK1)
      END IF
      REWIND (LUPROP)
      CALL REAPCM('NE-PCMIN','PCMGRD  ',LUPROP,WRK(KJENAO),NNBASX)
      CALL UTHU(WRK(KJENAO),WRK(KJEN),WRK(KUCMO),WRK(KWRK1),
     &          NBAST,NORBT)
      IF (NASHT .GT. 0) THEN
         CALL GETAC2(WRK(KJEN),WRK(KJENAC))
      END IF
      IF (IPRPCM .GE. 15) THEN
         WRITE (LUPRI,'(/A)') ' JEN_ao matrix:'
         CALL OUTPAK(WRK(KJENAO),NBAST,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' JEN_mo matrix:'
         CALL OUTPAK(WRK(KJEN),  NORBT,1,LUPRI)
         IF (NASHT .GT. 0) THEN
            WRITE (LUPRI,'(/A)') ' JEN_ac matrix:'
            CALL OUTPAK(WRK(KJENAC),NASHT,1,LUPRI)
         END IF
      END IF
C     
C     Expextation value of JEN (=PB)
C     
      TJEN     = SOLELM(DV,WRK(KJENAC),WRK(KJEN),TJENAC)
      IF (IPRPCM .GE. 6) THEN
         WRITE (LUPRI,'(A,F17.8)')
     *   ' --- JEN expectation value(=PB) :',TJEN
         WRITE (LUPRI,'(A,F17.8)')
     *   ' --- active part of JEN(=PB)    :',TJENAC
      END IF
Cbm   PB=-TJEN
Cbm   WRITE(LUPRI,*)'PB =',PB
C     
C     Read J2 in ao basis and transform to mo basis.
C     Extract active-active block in WRK(KJ2AC).
C     
      CALL DZERO(WRK(KDENC),N2BASX)
      CALL DZERO(WRK(KDENV),N2BASX)
      CALL FCKDEN((NISHT.GT.0),(NASHT.GT.0),WRK(KDENC),WRK(KDENV),
     &            CMO,DV,WRK(KWRK1),LW21)
      CALL DAXPY(N2BASX,1.0D0,WRK(KDENV),1,WRK(KDENC),1)
      CALL DZERO(WRK(KDENV),N2BASX)
      CALL DGEFSP(NBAST,WRK(KDENC),WRK(KDENV))
ckr      CALL DZERO(WRK(KDENC),N2BASX)
ckr      CALL PKSYM1(WRK(KDENV),WRK(KDENC),NBAS,NSYM,1)

      EXP1VL = .TRUE.
      TOFILE = .FALSE.
      CALL J1INT(WRK(KPOT),EXP1VL,WRK(KDENV),1,TOFILE,'NPETES ',
     &           1,WRK(KWRK1),LW21)
      CALL DZERO(QSE,MXTS)
      CALL V2Q(WRK(KWRK1),WRK(KPOT),QSE,QET,.FALSE.)
      CALL GPCLOSE(LUPCMD,'KEEP')
C      CALL DSCAL(NTS,-1.0D0,QSE,1)
C     
C     Read J1 (=ELECTROSTATIC POTENTIAL) in ao basis and transform to mo
C     basis.  Extract active-active block in WRK(KJ1AC).
C     
      XI = DIPORG(1)
      YI = DIPORG(2)
      ZI = DIPORG(3)
Ckr
Ckr should be parallelized, but a lot of information to send
Ckr However, in general not used for SCF and DFT optimizations
Ckr
      DO I = 1 , NTSIRR
         L = 1
         NCOMP = NSYM
         DIPORG(1) = XTSCOR(I)
         DIPORG(2) = YTSCOR(I)
         DIPORG(3) = ZTSCOR(I)
         EXP1VL    = .FALSE.
         TOFILE    = .FALSE.
         KPATOM    = 0
         TRIMAT    = .TRUE.
         CALL GET1IN(WRK(KJ1AO),'NPETES ',NCOMP,WRK(KWRK1),LW21,LABINT,
     &               INTREP,INTADR,L,TOFILE,KPATOM,TRIMAT,DUMMY,EXP1VL,
     &               DUMMY,IPRPCM)
         CALL UTHU(WRK(KJ1AO),WRK(KJ1),WRK(KUCMO),WRK(KWRK1),
     &        NBAST,NORBT)
         IF (NASHT .GT. 0) THEN
            CALL GETAC2(WRK(KJ1),WRK(KJ1AC))
         END IF
         IF (NONEQ) THEN
Clf here we will need the singlet density, I beleive.....
            WRK(KPOT + I - 1) = SOLELM(DV,WRK(KJ1AC),WRK(KJ1),TJ1AC)
         END IF

         IF (IPRPCM .GE. 15) THEN
            WRITE (LUPRI,'(/A)') ' J1_ao matrix:'
            CALL OUTPAK(WRK(KJ1AO),NBAST,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' J1_mo matrix:'
            CALL OUTPAK(WRK(KJ1),  NORBT,1,LUPRI)
            IF (NASHT .GT. 0) THEN
               WRITE (LUPRI,'(/A)') ' J1_ac matrix:'
               CALL OUTPAK(WRK(KJ1AC),NASHT,1,LUPRI)
            END IF
         END IF
         CALL DAXPY(NNORBX,-QSE(I),WRK(KJ1),1,WRK(KJEN),1)
C
C     Due to the computational cost of building the CI gradient, and the lack
C     of importance on convergence rates, we skip the construction of the 
C     solvent contribution to the diagonal Hessian. K.Ruud, Oct.-01
C
C         CALL DZERO(WRK(KGRDLM),NVARH)
C         IF (NCONF .GT. 1) THEN
C            CALL SOLGC(CREF,WRK(KJ1AC),TJ1AC,WRK(KGRDLM),INDXCI,
C     &                 WRK(KWRK1),LWRK1)
CC           CALL SOLGC(CREF,RLMAC,TELMAC,GLMCI,INDXCI,WRK,LWRK)
CC         END IF
C         IF (NWOPT .GT. 0) THEN
C            CALL SOLGO(DCVAL,DV,WRK(KJ1),WRK(KGRDLM+NCONF))
C         END IF
C         READ(LUGRDQ)WRK(KJ2GRD)
C         DO J =  NCONF, NVAR - 1 
C            WRK(KDIASH+J) = WRK(KDIASH+J) - WRK(KGRDLM+J)*WRK(KJ2GRD+J)
C         ENDDO
      ENDDO
Ckr
Ckr   End of parallelization loop
Ckr
C
C     Expextation value of J + X(0) = PB + PX
C     
      IF (NASHT .GT. 0) THEN
         CALL GETAC2(WRK(KJEN),WRK(KJENAC))
      END IF
      TJEN     = SOLELM(DV,WRK(KJENAC),WRK(KJEN),TJENAC)

      IF (IPRPCM .GE. 6) THEN
         CALL OUTPAK(WRK(KJENAC),NASHT,1,LUPRI)
         WRITE (LUPRI,'(A,F17.8)')
     *        ' --- TJEN expectation value :',TJEN
         WRITE (LUPRI,'(A,F17.8)')
     *        ' --- active part of TJEN    :',TJENAC
      END IF
C     
      KFREE=1
      ISPIN1=1
      ISPIN2=0
      CALL RSPDM(KREFSY,KREFSY,NHCREF,NHCREF,CREF,CREF,
     *     WRK(KTDV),DUMMY, ISPIN1,ISPIN2,.FALSE.,.TRUE.,
     *     INDXCI,WRK(KWRK1),KFREE,LWRK1)
Clf
      IF (.true.) THEN
         WRITE (LUPRI,'(/A)') ' JEN_ao matrix:'
         CALL OUTPAK(WRK(KJENAO),NBAST,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' JEN_mo matrix:'
         CALL OUTPAK(WRK(KJEN),  NORBT,1,LUPRI)
         IF (NASHT .GT. 0) THEN
            WRITE (LUPRI,'(/A)') ' JEN_ac matrix:'
            CALL OUTPAK(WRK(KJENAC),NASHT,1,LUPRI)
         END IF
      END IF


C     
C     TRIANGULAR PACKING OF ONE ELECTRON DENSITY MATRIX
C     
      CALL DSITSP(NASHT,WRK(KTDV),DV)
      CALL DZERO(WRK(KGRDLM),NVARH)
      IF (NCONF .GT. 1) THEN
         CALL SOLGC(CREF,WRK(KJENAC),TJENAC,WRK(KGRDLM),INDXCI,
     &              WRK(KWRK1),LWRK1)
C        CALL SOLGC(CREF,RLMAC,TELMAC,GLMCI,INDXCI,WRK,LWRK)
      END IF
      IF (NWOPT .GT. 0) THEN
         CALL SOLGO(D0,DV,WRK(KJEN),WRK(KGRDLM+NCONF))
      END IF
C     
C     
C     Obtain DIALM = diagonal TE(l,m) Hessian
C                  = 2 ( <i|R(l,m)|i> - TE(l,m) )
C     Add the DIALM contribution and the GRDLM contribution
C     to solvent Hessian diagonal.
C     
C      CALL SOLDIA(TJENAC,WRK(KJENAC),INDXCI,
C     *            WRK(KJEN),DV,WRK(KDIALM),WRK(KW21),LW21)
C     CALL SOLDIA(TELM,RLMAC,INDXCI,RLM,DV,DIAG,WRK,LFREE)
C     
      DO 420 I = 0,(NVAR-1)
         WRK(KDIASH+I) = WRK(KDIASH+I)
     *                 - WRK(KDIALM+I)
  420 CONTINUE
C     
C     test orthogonality
C     
      IF (IPRPCM .GE. 120) THEN
         WRITE (LUPRI,'(/A)')' >>> PCMGRD - grdj1, grdj2, dialm, '//
     &                      'diash, cref'
         DO 430 I = 1,NCONF
            WRITE (LUPRI,'(A,I10,3F10.6)') ' conf #',I,
     *            WRK(KDIALM-1+I),
     *            WRK(KDIASH-1+I),CREF(I)
  430    CONTINUE
      END IF
C     
      TEST = DDOT(NCONF,CREF,1,WRK(KGRDLM),1)
      IF (ABS(TEST) .GT. 1.D-8) THEN
         NWARN = NWARN + 1
         WRITE (LUPRI,'(/A,/A,1P,D12.4)')
     *   ' >>> PCMGRD WARNING, for B',
     *   ' <CREF | GRADB > =',TEST
      END IF
C     
C      TEST = DDOT(NCONF,CREF,1,WRK(KGRDJ1),1)
      IF (ABS(TEST) .GT. 1.D-8) THEN
         NWARN = NWARN + 1
         WRITE (LUPRI,'(/A,/A,1P,D12.4)')
     *   ' >>> PCMGRD WARNING, for J1 ',
     *   ' <CREF | GRADJ1 > =',TEST
      END IF
C      TEST = DDOT(NCONF,CREF,1,WRK(KGRDJ2),1)
      IF (ABS(TEST) .GT. 1.D-8) THEN
         NWARN = NWARN + 1
         WRITE (LUPRI,'(/A,/A,1P,D12.4)')
     *   ' >>> PCMGRD WARNING, for J2',
     *   ' <CREF | GRADJ2 > =',TEST
      END IF
C     
C     Add PCM gradient contribution to MCSCF gradient
C     
      CALL DAXPY(NVARH,DM1,WRK(KGRDLM),1,G,1)
clf      IF (IPRPCM .GE. 140) THEN
      IF (.true.) THEN
         WRITE (LUPRI,'(/A/A,2I10)')
     *      ' >>> PCMGRD - grdB, gtot (accum) - non-zero grdlm',
     *      '     NCONF, NWOPT =',NCONF,NWOPT
         DO 440 I = 1,NCONF
            IF (WRK(KGRDLM-1+I) .NE. D0)
     *         WRITE (LUPRI,'(A,I10,3F15.10)')
     *         ' conf #',I,WRK(KGRDLM-1+I),G(I)
  440    CONTINUE
         DO 450 I = NCONF+1,NVAR
            IF (WRK(KGRDLM-1+I) .NE. D0)
     *         WRITE (LUPRI,'(A,I10,3F15.10)')
     *         ' orb  #',I,WRK(KGRDLM-1+I),G(I)
  450    CONTINUE
      END IF
C     
C     
C     
      IF (IPRPCM .GE. 130) THEN
         WRITE (LUPRI,'(/A/A,2I10)')
     *      ' >>> PCMGRD - gtot (output) - non-zero elements',
     *      '     NCONF, NWOPT =',NCONF,NWOPT
         DO 840 I = 1,NCONF
            IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *      ' conf #',I,G(I)
  840    CONTINUE
         DO 850 I = NCONF+1,NVAR
            IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *      ' orb  #',I,G(I)
  850    CONTINUE
      END IF
      IF (LUPBKP .LT. 0) CALL GPCLOSE(LUPROP,'KEEP')
      CALL QEXIT('PCMGDT')
      RETURN
C     end of pcmgdt.
      END
