C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
#include "single.h"
C
#ifdef REV_LOG
===========================================================================
Revision 1.4  2000/05/24 18:45:09  hjj
new getref calls with appropriate NDREF (fixing CSF and triplet)
stop if solvent and CSF
bugfix: too little was allocated for KGP2
other minor changes
===========================================================================
#endif
C  /* Deck lagran */
      SUBROUTINE LAGRAN(WORD,FC,CMO,UDV,PV,XINDX,WRK,LWRK)
C
C 18-AUG-1991
C
C Called from GETGPV when label(5:8) = 'LAGR'
C
#include "implicit.h"
#include "dummy.h"
      DIMENSION CMO(*),UDV(NASHT,*),FC(*),PV(*),XINDX(*),WRK(*)
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
C
      PARAMETER (DHALF=0.5D0,ONE=1.0D0,DMONE=-1.0D0)
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
      IF (FLAG(16) .AND. NASHT.GT.1 .AND. NCREF .NE. KZCONF) THEN
Chj:     see comments in RSPSLV in rspsol.F
         WRITE(LUPRI,*)'Solvent not implemented in LAGRAN with CSFs'
         WRITE(LUPRI,*)'Use .DETERMINANTS and try again!'
         CALL QUIT('Solvent not implemented in LAGRAN with CSFs')
      END IF
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
         KQATD  = KFVTD  + NORBT * NORBT
         KQBTD  = KQATD  + NORBT * NASHT
         KH2    = KQBTD  + NORBT * NASHT
         KH2X   = KH2    + NORBT * NORBT
         KWRK1  = KH2X   + NORBT * NORBT
         LWRK1  = LWRK   - KWRK1
         IF (LWRK1.LT.0) CALL ERRWRK('LAGRAN',KWRK1-1,LWRK)
         CALL RSPOLI(1,0,UDV,WRK(KCREF),FC,DUMMY,PV,DUMMY,
     *            DUMMY,DUMMY, DUMMY,DUMMY,WRK(KFVTD),
     *            WRK(KQATD),WRK(KQBTD),WRK(KGP),
     *            WRK(KH2),WRK(KH2X),XINDX,CMO,WRK(KWRK1),LWRK1)
C        CALL RSPOLI(NCSIM,NOSIM,UDV,ZYCVEC,FC,FV,PVX,ZYMAT,FCONE,
C    *                  FVONE,QAONE,QBONE,FVTD,QATD,QBTD,EVECS,
C    *                  H2,H2X,XINDX,CMO,WRK,LWRK)
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
      IF (FLAG(16)) THEN
C
C        workspace allocation
C
         KGP2   = KTOT
         KTDV   = KGP2  + NVARH
Chj      KTDV   = KGP2  + KZVAR
Chj      ... SOLGDT uses NVARH (.gt. KZVAR !) for GP2
         KDV    = KTDV  + N2ASHX
         KWRK1  = KDV   + NNASHX
         KWRK2  = KWRK1 + NLMSOL * 2
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
         CALL SOLGDT(WRK(KCREF),CMO,XINDX,WRK(KDV),WRK(KGP2),
     *               ESOLT,WRK(KWRK1),WRK(KWRK2),LWRK1,NDREF,IREFSY)
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
      CALL GETGPV(LABEL,FC,CMO,UDV,PV,XINDX,ANTSYM,WRK(KGP),LWRK1)
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
     *            WRK(KIBTYP),WRK(KRESID),WRK(KEIVAL),WRK(KEIVEC),
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
#include "infrsp.h"
#include "inflr.h"
#include "infesr.h"
#include "rspprp.h"
#include "infpri.h"
C
      LOGICAL NEWDEF
      PARAMETER ( NTABLE = 6 )
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER*8 LABEL
C
      DATA TABLE /'.TRPPRP', '.SNGPRP', '.MAXIT ', '.THCESR',
     *            '.PRINT ', 'XXXXXXX'/
C
      NEWDEF = (WORD .EQ. '*ESR   ')
      ICHANG = 0
      IF (NEWDEF) THEN
      ESRCAL = .TRUE.
      WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') GO TO 100
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6), I
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
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' PROMPT "',WORD,
     *            '" NOT RECOGNIZED IN ESRINP.'
               CALL QUIT(' ILLEGAL PROMPT IN ESRINP ')
            END IF
         GO TO 100
      END IF
  300 CONTINUE
      WRITE(LUPRI,'(/A)')' ********* ESRINP ********'
      NTOESR = 0
      DO 500 I = 1,NPRLBL
         IF (ESROPS(I)) NTOESR = NTOESR + 1
         IF (ESROPT(I)) NTOESR = NTOESR + 1
 500  CONTINUE
      IF ( NTOESR.EQ. 0) THEN
         WRITE (LUPRI,'(/A)')
     *      ' *ESR   input ignored because no operators requested.'
         ESRCAL = .FALSE.
      ELSE
         CALL HEADER('Changes of defaults under *ESR   :',0)
         IF (ESRCAL) WRITE (LUPRI,'(/A,L1)')
     *      ' ESR - calculation carried out : ESRCAL='
     *      ,ESRCAL
         IF (MAXESR.NE.20) WRITE(LUPRI,'(/A,I6)')
     *      ' MAXIMUM NUMBER OF ITERATIONS. MAXESR =',MAXESR
         IF (THCESR.NE.1.0D-5) WRITE(LUPRI,'(/A,D12.6)')
     *      ' THRESHOLD FOR CONVERGENCE. THCESR =',THCESR
         IF (IPRESR.NE.2) WRITE(LUPRI,'(/A,I5)')
     *      ' PRINT LEVEL IN ESR ROUTINES. IPRESR =',IPRESR
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
#include "priunit.h"
#include "infpri.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "inflr.h"
#include "infdim.h"
#include "inforb.h"
#include "infesr.h"
C
      CALL QENTER('RSPESR')
      CALL HEADER('Output from RSPESR module',0)
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
               CALL GETGPV(LBESRT(KSYMOP,IOP),FC,CMO,UDV,PV,XINDX,
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
            CALL PRPAVE(LBESRT(KSYMOP,IOP),AVEVAL,CMO,
     *                  WRK(KTUDV),WRK(KWRK1),LWRK1)
            HYPFIN = AVEVAL + FOPLG
            WRITE(LUPRI,'(/2A,3(A,F10.6))')
     *        'TRIPLET OPERATOR: "',LBESRT(KSYMOP,IOP),
     *        '" LAGRANGIAN:',FOPLG,' AVERAGE:',AVEVAL,' TOTAL:',HYPFIN
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
C     ... so PRPAVE calculates singlet expectation values /hjaaj march 2003 ...
C         (if TRPLET true, then inactive density matrix is omitted!)
      DO 400 IOP = 1,NESRS(KSYMOP)
         CALL PRPAVE(LBESRS(KSYMOP,IOP),AVEVAL,CMO,
     *               UDV,WRK(KWRK1),LWRK1)
C        CALL PRPAVE(LABELOP,AVEVAL,CMO, UDV,WRK,LWRK)
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
C *** end of RSPESR --
      CALL QEXIT('RSPESR')
      RETURN
      END
C  /* Deck solgdt */
      SUBROUTINE SOLGDT(CREF,CMO,INDXCI,DV,G,ESOLT,ERLM,WRK,
     &                  LFREE,NHCREF,KREFSY)
C
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
     *         ' >>> SOLGDT - gtot (input) - non-zero elements',
     *         '     NCONF, NWOPT =',NCONF,NWOPT
            DO 40 I = 1,NCONF
               IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *         ' conf #',I,G(I)
   40       CONTINUE
            DO 50 I = NCONF+1,NVAR
               IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *         ' orb  #',I,G(I)
   50       CONTINUE
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
         CALL SOLDIA(TELMAC,WRK(KRLMAC),INDXCI,
     *               WRK(KRLM),DV,WRK(KDIALM),WRK(KW21),LW21)
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
C     end of solgrd.
      END
