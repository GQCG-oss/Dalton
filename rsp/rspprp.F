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
C===========================================================================
CRevision 1.2  2000/05/24 19:05:49  hjj
Cnew getref calls with appropriate NDREF instead of NCREF
C(fixing error for triplet with CSF)
C
C931130-hjaaj
C  merged VCD changes from Oslo with rspprp.u from here
C  renamed GPLON to PRPABA
C  DECK PRPABA moved to abacus 'abarspn.u' because it uses aba.cdk ...
C===========================================================================
#endif
C  /* Deck pprst */
      SUBROUTINE PPRST(IBTYP,REDS,REDE,UDV,XINDX,WRK,LWRK)
C
C PURPOSE: CREATE START VECTOR(S) FOR SOLUTION OF A LINEAR
C          SET OF EQUATIONS BASED ON OLD TRIAL VECTORS
C
C List of updates
C 21-Jul-1992 Hinne Hettema error if RSPAVE necessary (SUPSYM in MCSCF)
C
C   Used from common blocks:
C     INFORB : ...,N2ORBX,N2ASHX,....
C
#include <implicit.h>
#include <dummy.h>
C
      DIMENSION IBTYP(*),REDS(*),REDE(*),UDV(*)
      DIMENSION XINDX(*),WRK(*)
C
#include <ibndxdef.h>
C
#include <priunit.h>
#include <wrkrsp.h>
#include <inftap.h>
#include <infrsp.h>
#include <infpri.h>
#include <inforb.h>
C
      LOGICAL FNDLAB
C
      IF (RSPSUP .AND. (KSYMOP.EQ.1)) THEN
         WRITE(LUERR,*)
     *   'PPRST: RSPAVE necessary for restart but not yet implemented'
         CALL QUIT('PPRST: Implement RSPAVE for restart')
      END IF
C
      LMAXB  = MAX(KZYWOP,KZCONF)
C     ... max. length of a BOVEC and a BCVEC
      LWRKA  = LMAXB + NCREF + N2ASHX + 100
C     ... needed in RSPRED and RSPRED.RSPSLI irrespective of NSIM
      LWRKB  = LMAXB + KZYVAR + N2ORBX
C     ... needed for each trial vector here and in RSPRED.RSPSLI
      MAXSIM = (LWRK - LWRKA)/LWRKB
C
      KZRED  = 0
      LURSP1 = -1
      CALL GPOPEN(LURSP1,'RSPRST.E2C','UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND (LURSP1)
      IF ( .NOT.FNDLAB('RESTART ',LURSP1) ) THEN
         WRITE (LUPRI,'(/2A,I4)')
     *   ' PPRST: Restart not possible, RESTART label not found on',
     *   ' LURSP1 =',LURSP1
         CALL QUIT('PPRST: Restart impossible, RESTART label not found')
      END IF
      READ (LURSP1)JSYMOP, KZDIM,(IBTYP(I),I=1,KZDIM)
      CALL GPCLOSE(LURSP1,'KEEP')
      IF (JSYMOP.NE.KSYMOP) THEN
         WRITE (LUPRI,'(/,(A,I4))')
     *   ' PPRST: Restart not possible, symmetry on restart file',
     *   JSYMOP,' symmetry of current operator',KSYMOP
         CALL QUIT('PPRST: Restart impossible, no match of symmetry ')
      END IF
C
C Check how many sigma vectors have been written to LURSP5
C
      CALL REWSPL(LURSP5)
      JZDIM = 0
      JRSP5 = 0
      DO 100 I = 1,KZDIM
#if defined (VAR_SPLITFILES)
         CALL READST(LURSP5,LRSP5,DUMMY,JRSP5,.FALSE.)
         IF (LRSP5 .EQ. -1) GOTO 110
#else
         READ (LURSP5,END=110,ERR=110) JUNK
#endif
         JZDIM = I
  100 CONTINUE
  110 KZDIM = JZDIM
C
C READ IN OLD TRIAL VECTORS AND SET UP REDUCED MATRICES
C
C REPEAT UNTILL
C
         NLOAD  = 0
         JZDIM  = KZDIM
         CALL REWSPL(LURSP3)
         CALL REWSPL(LURSP5)
#if defined (VAR_SPLITFILES)
         JRSP3 = 0
         JRSP5 = 0
#endif
 200     CONTINUE
            NLOAD  = NLOAD +1
            IF ((NLOAD.EQ.1).AND.(KOFFTY.EQ.1)) THEN
#if defined (VAR_SPLITFILES)
               CALL READST(LURSP3,LRSP3,DUMMY,JRSP3,.FALSE.)
               IF (LRSP3 .EQ. -1) GOTO 910
               CALL READST(LURSP5,LRSP5,DUMMY,JRSP5,.FALSE.)
               IF (LRSP5 .EQ. -1) GOTO 920
#else
               READ (LURSP3,ERR=910)
               READ (LURSP5,ERR=920)
#endif
               JZDIM = JZDIM - 1
            ENDIF
            IF ( NLOAD .GT. 1) THEN
               CALL REWSPL(LURSP3)
#if defined (VAR_SPLITFILES)
               JRSP3 = 0
               IF (KOFFTY.EQ.1) CALL READST(LURSP3,KZCONF,DUMMY,JRSP3,
     &                                      .FALSE.)
#else
               IF ( KOFFTY.EQ.1)READ (LURSP3)
#endif
               DO 61 IVEC = 1,KZRED
#if defined (VAR_SPLITFILES)
                  CALL READST(LURSP3,LRSP3,DUMMY,JRSP3,.FALSE.)
#else
                  READ (LURSP3)
#endif
 61            CONTINUE
            END IF
            NSIM   = MIN(MAXSIM,JZDIM)
            NCSIM  = 0
            NOSIM  = 0
            DO 410 ISIM=1,NSIM
               IF (IBTYP(KOFFTY+KZRED+ISIM).EQ.JBCNDX) THEN
                  NCSIM=NCSIM+1
               ELSE
                  NOSIM=NOSIM+1
               ENDIF
 410        CONTINUE
            KBCVEC = 1
            KBOVEC = KBCVEC + NCSIM*KZCONF
            KECVEC = KBOVEC + NOSIM*KZYWOP
            KEOVEC = KECVEC + NCSIM*KZYVAR
            KDIAE  = KEOVEC + NOSIM*KZYVAR
            IF (SOPPA) THEN
               KWRK = KDIAE + KZCONF
            ELSE
               KWRK = KDIAE
            ENDIF
            LWRKRE = LWRK   - KWRK
            IF (LWRKRE.LT.0) CALL ERRWRK('PPRST',KWRK-1,LWRK)
            ISTBO  = KBOVEC
            ISTBC  = KBCVEC
            ISTEC  = KECVEC
            ISTEO  = KEOVEC
            IF (SOPPA) THEN
               CALL REWSPL(LURSP4)
               JRSP4 = 0
#if defined (VAR_SPLITFILES)
               CALL READST(LURSP4,LRSP4,WRK(KDIAE),JRSP4,.TRUE.)
               IF (LRSP4 .NE. KZCONF) CALL QUIT('Inconsistent record '//
     &              'length on unit LURSP4 in PPRST')
#else
               CALL READT(LURSP4,KZCONF,WRK(KDIAE))
#endif
            ENDIF
            DO 400 ISIM=1,NSIM
               IF (IBTYP(KOFFTY+KZRED+ISIM).EQ.JBCNDX) THEN
#if defined (VAR_SPLITFILES)
                  CALL READST(LURSP3,LRSP3,WRK(ISTBC),JRSP3,.TRUE.)
                  IF (LRSP3 .NE. KZCONF) CALL QUIT('Inconsistent '//
     &                 'record length on unit LURSP3 in PPRST')
#else
                  CALL READT(LURSP3,KZCONF,WRK(ISTBC))
#endif
C
C   If SOPPA read in the p-h part of the transformed vector
C   and construct the 2p-2h part and put everything as if read
C
                  IF (SOPPA) THEN
#if defined (VAR_SPLITFILES)
                     CALL READST(LURSP5,LRSP5,
     &                           WRK(ISTEC+KZCONF),JRSP5,.TRUE.)
                     IF (LRSP5 .NE. KZYWOP) CALL QUIT('Inconsistent '//
     &                    'record length on unit LURSP5 in PPRST')
#else
                     CALL READT(LURSP5,KZYWOP,WRK(ISTEC+KZCONF))
#endif
                     CALL DCOPY(KZWOPT,WRK(ISTEC+KZVAR),1,
     *                          WRK(ISTEC+KZVAR+KZCONF),1)
                     DO I=0,KZCONF-1
                        WRK(ISTEC+I) = WRK(KDIAE+I) * WRK(ISTBC+I)
                     ENDDO
                     CALL DZERO(WRK(ISTEC+KZVAR),KZCONF)
                  ELSE
#if defined (VAR_SPLITFILES)
                     CALL READST(LURSP5,LRSP5,WRK(ISTEC),JRSP5,.TRUE.)
                     IF (LRSP5 .NE. KZYVAR) CALL QUIT('Inconsistent '//
     &                    'record length on unit LURSP5 in PPRST (2)')
#else
                     CALL READT(LURSP5,KZYVAR,WRK(ISTEC))
#endif
                  ENDIF
                  ISTBC = ISTBC + KZCONF
                  ISTEC = ISTEC + KZYVAR
               ELSE
#if defined (VAR_SPLITFILES)
                  CALL READST(LURSP3,LRSP3,WRK(ISTBO),JRSP3,.TRUE.)
                  IF (LRSP3 .NE. KZYWOP) CALL QUIT('Inconsistent '//
     &                 'record length on unit LURSP3 in PPRST (2)')
                  CALL READST(LURSP5,LRSP5,WRK(ISTEO),JRSP5,.TRUE.)
                  IF (LRSP5 .NE. KZYVAR) CALL QUIT('Inconsistent '//
     &                 'record length on unit LURSP5 in PPRST (3)')
#else
                  CALL READT(LURSP3,KZYWOP,WRK(ISTBO))
                  CALL READT(LURSP5,KZYVAR,WRK(ISTEO))
#endif
                  ISTBO = ISTBO + KZYWOP
                  ISTEO = ISTEO + KZYVAR
               ENDIF
 400        CONTINUE
            KZRED  = KZRED+NSIM
            KZYRED = 2*KZRED
            JZDIM  = JZDIM-NSIM
C
C           CALL RSPRED(1,..) INCREASE DIMENSION OF REDUCED RSP EQUATION
C
            CALL RSPRED (1,.FALSE.,NSIM,IBTYP,DUMMY,DUMMY,REDE,REDS,
     *                 DUMMY,DUMMY,WRK(KBCVEC),WRK(KBOVEC),
     *                 UDV,WRK(KECVEC),XINDX,WRK(KWRK),LWRKRE)
C
C           CALL RSPRED (ICTL,LINEQ,N,IBTYP,GD,REDGD,REDE,REDS,
C    *                 EIVAL,EIVEC,BCVEC,BOVEC,UDV,EVECS,XINDX,
C    *                 WRK,LWRK)
C
            IF (JZDIM.LE.0) THEN
                IF (IPRRSP.GE.10) WRITE (LUPRI,'(//A,I5,A/A/)')
     *             ' ***',KZDIM,' OLD TRIAL VECTORS FOUND',
     *             ' CALLING RSPNEX TO GET NEXT TRIAL VECTOR.'
            ELSE
               GO TO 200
C        ^--------------
            END IF
C
      RETURN
  910 CONTINUE
      WRITE (LUPRI,'(/A)')
     *     ' PPRST: premature end-of-file on LURSP3 during restart.'
      CALL QUIT('PPRST: premature end-of-file on LURSP3 during restart')
  920 CONTINUE
      WRITE (LUPRI,'(/A)')
     *     ' PPRST: premature end-of-file on LURSP5 during restart.'
      CALL QUIT('PPRST: premature end-of-file on LURSP5 during restart')
C
C     END OF PPRST.
C
      END
C  /* Deck prpcon */
      SUBROUTINE PRPCON(NSIM,PRPMO,GP,XINDX,WRK,LWRK)
C
C CALCULATE CONFIGURATION PART OF GRADIENT VECTOR
C
C                (  <J,M,0> )
C                (         )
C                ( -<0,M,J> )
C
C  |0> REFERENCE STATE OF SYMMETRY IREFSY
C  |J> DETERMINANT OF SYMMETRY KSYMST
C
C  M IS A ONE ELECTRON OPERATOR OF SYMMETRY KSYMOP
C  WHICH IS STORED IN PRPMO
C
C  OUTPUT: CONFIGURATION PART OF GP VECTOR
C
C
C
#include <implicit.h>
C
      DIMENSION PRPMO(NORBT,*),GP(KZYVAR,*),XINDX(*),WRK(*)
C
#include <priunit.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infopt.h>
#include <inforb.h>
#include <infpri.h>
#include <infdim.h>
C
      LOGICAL NOH2, IH8SM
C
      PARAMETER ( D1 = 1.0D0 , DM1 = -1.0D0 , DUMMY = 1.0D+25 )
C
C ALLOCATE WORK SPACE
C
      KUPRP  = 1
      KCREF  = KUPRP + N2ASHX
      IF (IREFSY .EQ. KSYMST) THEN
         NDREF = KZCONF
Chj      ... if triplet, we need CREF in determinants
      ELSE
         NDREF = NCREF
      END IF
      KFREE  = KCREF + NDREF
      LFREE  = LWRK  - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('PRPCON',KFREE-1,LWRK)
      CALL GETREF(WRK(KCREF),NDREF)
C
C DEFINE PARAMETERS FOR LINEAR TRANSFORMATION
C
      NOH2   = .TRUE.
      IH8SM  = .TRUE.
      DO 100 ISIM = 1,NSIM
C
C EXTRACT ACTIVE PART OF PROPERTY MATRIX
C        (note: CISIGD requires UPRP(I,J) = PRPMO(J,I))
C
         CALL DZERO(WRK(KUPRP),N2ASHX)
         DO 50 ISYM = 1,NSYM
            JSYM = MULD2H(KSYMOP,ISYM)
            NASHI = NASH(ISYM)
            NISHI = NISH(ISYM)
            IORBI = IORB(ISYM)
            NASHJ = NASH(JSYM)
            NISHJ = NISH(JSYM)
            IORBJ = IORB(JSYM)
            DO 70 JA = 1,NASHJ
               DO 60 IA = 1,NASHI
                  WRK(KUPRP-1+(IA+IASH(ISYM)-1)*NASHT+JA+IASH(JSYM)) =
     &                PRPMO(IORBI+NISHI+IA,IORBJ+NISHJ+JA)
 60            CONTINUE
 70         CONTINUE
 50      CONTINUE
         IF (IPRRSP.GT.110) THEN
            WRITE(LUPRI,'(/A)')' ACTIVE PART OF PROPERTY MATRIX'
            CALL OUTPUT(WRK(KUPRP),1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         END IF
         ISPIN1 = 0
         IF (TRPLET) THEN
            ISPIN2 = 1
         ELSE
            ISPIN2 = 0
         END IF
         CALL CISIGD(IREFSY,KSYMST,NDREF,KZCONF,WRK(KCREF),
     *               GP(1,ISIM),WRK(KUPRP),DUMMY,
     *               NOH2,IH8SM,XINDX,ISPIN1,ISPIN2,WRK(KFREE),LFREE)
C        CALL RSPSIG(ICSYM,IHCSYM,NCDET,NHCDET,C,HC,FCAC,H2AC,IFLAG,
C    *                  NOH2,WORK,KFREE,LFREE)
         IF (IREFSY .EQ. KSYMST) THEN
C           ... remove CREF component of GP vector
            IF ( IPRRSP.GT.110 ) THEN
               WRITE(LUPRI,'(/A)')
     *         ' LINEAR TRANSFORMED PROPERTY VECTOR BEFORE PROJECT.'
               CALL OUTPUT(GP(1,ISIM),1,KZCONF,1,1,KZCONF,1,1,LUPRI)
            END IF
            T1 = DDOT(KZCONF,WRK(KCREF),1,GP(1,ISIM),1)
            CALL DAXPY(KZCONF,(-T1),WRK(KCREF),1,GP(1,ISIM),1)
            IF (IPRRSP.GT.110) THEN
               WRITE(LUPRI,'(/A)')
     *        ' FROM Z COMPONENT OF CONFIGURATION PROPERTY VECTOR'
               WRITE(LUPRI,'(A,1P,D15.8)')
     *        ' AVERAGE VALUE OF ONE ELECTRON OPERATOR, ACTIVE PART:',T1
            END IF
         END IF
         IF ( IPRRSP.GT.110 ) THEN
            WRITE(LUPRI,'(/A)')
     *      ' LINEAR TRANSFORMED PROPERTY VECTOR '
            CALL OUTPUT(GP(1,ISIM),1,KZCONF,1,1,KZCONF,1,1,LUPRI)
         END IF
C
C EXTRACT ACTIVE PART OF TRANSPOSED PROPERTY MATRIX
C        (note: CISIGD requires UPRPT(I,J) = PRPMOT(J,I) = PRPMO(I,J))
C
         CALL DZERO(WRK(KUPRP),N2ASHX)
         DO 150 ISYM = 1,NSYM
            JSYM = MULD2H(KSYMOP,ISYM)
            NASHI = NASH(ISYM)
            NISHI = NISH(ISYM)
            IORBI = IORB(ISYM)
            NASHJ = NASH(JSYM)
            NISHJ = NISH(JSYM)
            IORBJ = IORB(JSYM)
            IASHJ = IASH(JSYM)
            DO 170 JA = 1,NASHJ
               DO 160 IA = 1,NASHI
                  WRK(KUPRP-1+(JA+IASH(JSYM)-1)*NASHT+IA+IASH(ISYM)) =
     *                PRPMO(IORBI+NISHI+IA,IORBJ+NISHJ+JA)
 160           CONTINUE
 170        CONTINUE
 150     CONTINUE
         IF (IPRRSP.GT.110) THEN
            WRITE(LUPRI,'(/A)')
     *            ' ACTIVE PART OF TRANSPOSED PROPERTY MATRIX'
            CALL OUTPUT(WRK(KUPRP),1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         END IF
         ISPIN1 = 0
         IF (TRPLET) THEN
            ISPIN2 = 1
         ELSE
            ISPIN2 = 0
         END IF
         CALL CISIGD(IREFSY,KSYMST,NDREF,KZCONF,WRK(KCREF),
     *               GP(1+KZVAR,ISIM),WRK(KUPRP),DUMMY,
     *               NOH2,IH8SM,XINDX,ISPIN1,ISPIN2,WRK(KFREE),LFREE)
         CALL DSCAL(KZCONF,DM1,GP(1+KZVAR,ISIM),1)
         IF (IREFSY .EQ. KSYMST) THEN
C           ... remove CREF component of GP vector
            IF ( IPRRSP.GT.110 ) THEN
               WRITE(LUPRI,'(/A)')
     *        ' LINEAR TRANSFORMED PROPERTY VECTOR BEFORE PROJECT.'
             CALL OUTPUT(GP(1+KZVAR,ISIM),1,KZCONF,1,1,KZCONF,1,1,LUPRI)
            END IF
            T1 = DDOT(KZCONF,WRK(KCREF),1,GP(1+KZVAR,ISIM),1)
            CALL DAXPY(KZCONF,(-T1),WRK(KCREF),1,GP(1+KZVAR,ISIM),1)
            IF (IPRRSP.GT.110) THEN
               WRITE(LUPRI,'(/A)')
     *        ' FROM Y COMPONENT OF CONFIGURATION PROPERTY VECTOR'
               WRITE(LUPRI,'(A,1P,D15.8)')
     *       ' AVERAGE VALUE OF ONE ELECTRON OPERATOR, ACTIVE PART:',T1
            END IF
         END IF
         IF ( IPRRSP.GT.110 ) THEN
            WRITE(LUPRI,'(/2A)')
     *      ' LINEAR TRANSFORMED GRADIENT VECTOR'
     *      ,' WITH TRANSPOSED PROPERTY MATRIX'
            CALL OUTPUT(GP(1+KZVAR,ISIM),
     *                  1,KZCONF,1,1,KZCONF,1,1,LUPRI)
         END IF
 100  CONTINUE
C
C END OF PRPCON
C
      RETURN
      END
C  /* Deck getgpv */
      SUBROUTINE GETGPV(WORD,FC,FV,CMO,UDV,PV,XINDX,ANTSYM,WRK,LWRK)
C
C PURPOSE
C  CALCULATE GP VECTOR AND RETURN AS FIRST ELEMENTS IN WRK
C
#include <implicit.h>
#include <maxorb.h>
C
      DIMENSION FC(*),CMO(*),UDV(*),PV(*),XINDX(*),WRK(LWRK)
      CHARACTER*8 WORD
C
#include <infvar.h>
#include <wrkrsp.h>
C
      CALL QENTER('GETGPV')
      IF (WORD(3:8).EQ.'SPNORB' .OR. WORD(3:8) .EQ. 'MNF-SO' .OR. 
     &   WORD(3:8).EQ.'SPNSCA') THEN
         KGP   = 1
         KWRK1 = KGP + KZYVAR
         LWRK1 = LWRK + 1 - KWRK1
         CALL HSOCTL(WORD,WRK(KGP),CMO,UDV,PV,XINDX,
     &               ANTSYM,WRK(KWRK1),LWRK1)
      ELSE IF (WORD(5:8).EQ.'LAGR') THEN
         CALL LAGRAN(WORD,FC,FV,CMO,UDV,PV,XINDX,WRK,LWRK)
         ANTSYM = -1.0D0
      ELSE IF (WORD(2:7).EQ.'LONMAG') THEN
         CALL PRPABA(WORD,WRK,LWRK)
      ELSE
         KGP   = 1
         KWRK1 = KGP + KZYVAR
         LWRK1 = LWRK + 1 - KWRK1
         CALL PRPCTL(WORD,WRK(KGP),CMO,UDV,PV,XINDX,ANTSYM,
     &               WRK(KWRK1),LWRK1)
      ENDIF
      CALL QEXIT('GETGPV')
      RETURN
      END
C  /* Deck prpctl */
      SUBROUTINE PRPCTL (WORD,GP,CMO,UDV,PV,XINDX,ANTSYM,WRK,LWRK)
C
C 18-Feb-1987
C
#include <implicit.h>
      DIMENSION GP(*),CMO(*),UDV(*),PV(*),XINDX(*),WRK(*)
      CHARACTER*8 WORD
C-- common blocks:
#include <priunit.h>
#include <infdim.h>
#include <inforb.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <infpri.h>
#include <infsop.h>
C
      CALL QENTER('PRPCTL')
      IF (IPRRSP .GT. 10) WRITE (LUPRI,7010) WORD
 7010 FORMAT(//,' Output from PRPCTL:',
     *        /,' ===================',//,' Property label = ',A,/)
C
C ALLOCATE WORK SPACE
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWRK
      CALL MEMGET('REAL',KPRPMO,N2ORBX,WRK,KFREE,LFREE)
C If SOPPA then make workspace for the second order ph transition moments.
      IF (SOPPA) THEN
         CALL MEMGET('REAL',KPRPSE,N2ORBX,WRK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KPRPSE,     0,WRK,KFREE,LFREE)
      ENDIF
C
      CALL DZERO(GP,KZYVAR)
C
C Read AO property integrals and transform to MO basis.
C
      KSYMP = -1
      CALL PRPGET (WORD,CMO,WRK(KPRPMO),KSYMP,ANTSYM,
     &             WRK(KFREE),LFREE,IPRRSP)
C     CALL PRPGET (WORD,CMO,PRPMO,KSYMP,ANTSYM,WRK,LWRK,IPRINT)
      IF (KSYMP.NE.KSYMOP) THEN
         WRITE (LUPRI,'(/A/2A/A,2I5/A,F10.2)')
     &      'FATAL ERROR: KSYMOP .ne. KSYMP from PRPGET',
     &      '   Property label  : ',WORD,
     &      '   KSYMOP and KSYMP:',KSYMOP,KSYMP,
     &      '   ANTSYM          :',ANTSYM
         CALL QUIT('KSYMOP .ne. KSYMP from PRPGET')
      END IF
C
C Print atomic and molecular property integrals, if desired
C
      IF (IPRRSP.GT.125) THEN
         WRITE (LUPRI,'(/2A)')
     &      'PRPCTL test output: MOLECULAR PROPERTY INTEGRALS: ', WORD
         CALL OUTPUT(WRK(KPRPMO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C CALCULATE AVERAGE VALUE OF ONE ELECTRON OPERATOR
C
      IF ((KSYMOP.EQ.1).AND.(IREFSY.EQ.KSYMST))THEN
         IPRONE = 10
         IF (IPRRSP .GE. IPRONE) THEN
            WRITE (LUPRI,'(/2A)') 'PRPCTL test output: '//
     *      ' Reference state expectation value for ',WORD
            CALL PRPONE(WRK(KPRPMO),UDV,ONETOT,IPRONE,WORD)
         END IF
      END IF
C
C Construct orbital property vector
C
      IF (KZWOPT.GT.0) THEN
         IF (SOPPA .AND. .NOT. HIRPA) THEN
C        ... make second order TMs in WRK(KPRPSE).
            CALL PRPOMP(WRK(KPRPMO),WRK(KPRPSE),UDV)
C        ... CONSTRUCT TMS TO SECOND ORDER
            CALL PRPORS(WRK(KPRPMO),WRK(KPRPSE),GP)
         ELSE
            CALL PRPORB(WRK(KPRPMO),UDV,GP)
C           CALL PRPORB(PRPMO,UDV,GP)
         ENDIF
      ENDIF
C
C Construct configuration property vector
C
      IF (KZCONF.GT.0) THEN
         IF (SOPPA) THEN
           IF (.NOT.HIRPA) THEN
C        ... CONSTRUCT 2P-2H TRANSITION MOMENTS
            IF (TRPLET) THEN
               CALL PRPTCMP(WRK(KPRPMO),PV,GP,XINDX(KABSAD),
     &                      XINDX(KABTAD),XINDX(KIJ1AD),XINDX(KIJ2AD),
     &                      XINDX(KIJ3AD),XINDX(KIADR1),
     &                      WRK,KFREE,LFREE)
            ELSE
               CALL PRPCMP(WRK(KPRPMO),PV,GP,XINDX(KABSAD),
     &                     XINDX(KABTAD),XINDX(KIJSAD),XINDX(KIJTAD),
     &                     XINDX(KIADR1),WRK,KFREE,LFREE)
            ENDIF
           ENDIF
         ELSE
            CALL PRPCON(1,WRK(KPRPMO),GP,XINDX,WRK(KFREE),LFREE)
C           CALL PRPCON(NSIM,PRPMO,GP,XINDX,WRK,LWRK)
         ENDIF
      ENDIF
C ---------------------------------------------------------------
C
      CALL MEMREL('PRPCTL',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('PRPCTL')
      RETURN
      END
C  /* Deck prpget */
      SUBROUTINE PRPGET (WORD,CMO,PRPMO,KSYMP,ANTSYM,WRK,LWRK,IPRINT)
C
C Written 18-Feb-1987 PJ
C
C 1. READ IN SYMMETRY TRANSFORMED AO PROPERTY INTEGRALS FROM LUPROP
C 2. TRANSFORM PROPERTY INTEGRALS FROM AO BASIS TO MO BASIS
C
C Revised Sep. 2003 hjaaj:  IPRINT in param.list,
C  if KSYMP .le. 0 then return symmetry in KSYMP, otherwise use
C  user-specified KSYMP
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
      CHARACTER*8 WORD, RTNLBL(2)
      DIMENSION PRPMO(NORBT,*), CMO(*), WRK(*)
C
C -- common blocks:
C inforb.h : N2BASX, ...
C infind.h : ISAO(*)
C inftap.h : LUPROP
#include <inforb.h>
#include <maxorb.h>
#include <maxash.h>
#include <infind.h>
#include <inftap.h>
C -- local constants:
      PARAMETER (DM1 = -1.0D0, D0 = 0.0D0, D1 =1.0D0 )
C
      LOGICAL FNDLB2
C
      CALL QENTER('PRPGET')
C
C ALLOCATE WORK SPACE
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWRK
      CALL MEMGET('REAL',KPRPAO,N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KWRK1,NNBASX,WRK,KFREE,LFREE)
C
C 1. READ IN AO PROPERTY INTEGRALS
C
      IF (LUPROP .LE. 0) CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &               'UNFORMATTED',IDUMMY,.FALSE.)
C
      REWIND LUPROP
      IF ( FNDLB2(WORD,RTNLBL,LUPROP)) THEN
         IF (RTNLBL(2).EQ.'SQUARE  ') THEN
            ANTSYM = D0
            CALL READT(LUPROP,N2BASX,WRK(KPRPAO))
         ELSE IF (RTNLBL(2).EQ.'SYMMETRI') THEN
            ANTSYM = D1
            CALL READT(LUPROP,NNBASX,WRK(KWRK1))
            CALL DSPTSI(NBAST,WRK(KWRK1),WRK(KPRPAO))
         ELSE IF (RTNLBL(2).EQ.'ANTISYMM') THEN
            ANTSYM = DM1
            CALL READT(LUPROP,NNBASX,WRK(KWRK1))
            CALL DAPTGE(NBAST,WRK(KWRK1),WRK(KPRPAO))
         ELSE
            CALL QUIT('Error: No antisymmetry label on LUPROP')
         END IF
      ELSE
         WRITE (LUPRI,'(//3A)') ' >>> PRPGET: PROPERTY "',WORD,
     *        '" NOT FOUND ON AOPROPER.'
         CALL QUIT('PRPGET: PROPERTY NOT FOUND ON AOPROPER.')
      END IF
      CALL GPCLOSE(LUPROP,'KEEP')
C
      IF (IPRINT.GT.120) THEN
         WRITE(LUPRI,'(/3A)') ' PRPGET: PROPERTY MATRIX ',WORD,
     *         ' IN AO BASIS AS READ FROM LUPROP'
         IF (ANTSYM .NE. D0) THEN
            CALL OUTPAK(WRK(KWRK1),NBAST,1,LUPRI)
            WRITE(LUPRI,'(/3A)') ' PROPERTY MATRIX ',WORD,
     *         ' IN AO BASIS AFTER (ANTI)SYMMETRIZATION '
         END IF
         CALL OUTPUT(WRK(KPRPAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C     Determine spatial symmetry of property matrix
C     (as identified by abs largest element)
C
C     FIXME TODO  /hjaaj sep 2003:
C     What if matrix is a zero matrix ? should we return KSYMP = -1 ??
C
      IF (KSYMP .LE. 0) THEN
C        if KSYMP .gt. 0 use user-specified symmetry
         INDPMX = IDAMAX(N2BASX,WRK(KPRPAO),1)
         IIND = 1 + (INDPMX-1)/NBAST
         JIND = INDPMX - (IIND-1)*NBAST
         KSYMP = MULD2H( ISAO(IIND) , ISAO(JIND) )
       END IF
C
C 2. TRANSFORM PROPERTY INTEGRALS FROM AO SYMMETRY TO MO BASIS
C
      CALL DZERO(PRPMO,N2ORBX)
      DO 600 ISYM=1,NSYM
         JSYM=MULD2H(ISYM,KSYMP)
      IF((NORB(ISYM).LE.0).OR.(NORB(JSYM).LE.0))GO TO 600
         CALL UTHV(CMO(ICMO(ISYM)+1),WRK(KPRPAO),CMO(ICMO(JSYM)+1),
     &             ISYM,JSYM,NBAS(ISYM),NBAS(JSYM),
     &             PRPMO,WRK(KFREE))
C        CALL UTHV(U,PRPAO,V,ISYM,JSYM,NBASI,NBASJ,PRPMO,WRK)
         IF (IPRINT.GT.110) THEN
            WRITE(LUPRI,'(/3A/A,I5,A,I5)')
     &         ' PROPERTY: ',WORD,' IN MO. BASIS after',
     &         ' ISYM=',ISYM,' JSYM=',JSYM
            CALL OUTPUT(PRPMO,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
         IF (IPRINT.GT.140) THEN
            WRITE(LUPRI,'(/A,I5,A)')' MO. COEFFICIENTS FOR SYMMETRY',
     &                                ISYM
            CALL OUTPUT(CMO(ICMO(ISYM)+1),1,NBAS(ISYM),1,NORB(ISYM),
     $                  NBAS(ISYM),NORB(ISYM),1,LUPRI)
            IF (ISYM .NE. JSYM) THEN
            WRITE(LUPRI,'(/A,I5,A)')' MO. COEFFICIENTS FOR SYMMETRY',
     &                                JSYM
            CALL OUTPUT(CMO(ICMO(JSYM)+1),1,NBAS(JSYM),1,NORB(JSYM),
     $                  NBAS(JSYM),NORB(JSYM),1,LUPRI)
            END IF
         END IF
 600  CONTINUE
C
C     If "Half Differentiated Overlap"
C     then extract antisymmetric component of matrix
C
      IF (WORD(1:3) .EQ. 'HDO') THEN
         CALL DGETAP(NORBT,PRPMO,WRK(KPRPAO))
         CALL DAPTGE(NORBT,WRK(KPRPAO),PRPMO)
      END IF
      IF (IPRINT.GE.90) THEN
         WRITE(LUPRI,'(/3A)')' PROPERTY: ',WORD,' IN MO. BASIS'
         CALL OUTPUT(PRPMO,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     Whereas ANTSYM above has referred matrix symmetry, it will from here on
C     refer to the antisymmetry of the response vector, thus we change sign
C     on it. K.Ruud, July 10 2000.
C
      ANTSYM = -ANTSYM
      CALL MEMREL('PRPGET',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('PRPGET')
      RETURN
      END
C  /* Deck prporb */
      SUBROUTINE PRPORB(PRPMO,UDV,GP)
C
C WRITTEN 14-FEB 1986
C
C PURPOSE:
C    DISTRIBUTE PROPERTY MO INTEGRALS INTO GP VECTORS
C
C List of updates
C 21-Jul-1992 Hinne Hettema Average orbital part if necessary.
C
#include <implicit.h>
C
      DIMENSION PRPMO(NORBT,NORBT)
      DIMENSION UDV(NASHDI,NASHDI),GP(KZYVAR)
C
C  INFDIM : NASHDI
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <infvar.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infpri.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
C -- local constants
C
      PARAMETER ( D2 = 2.0D0 )
C
C DISTRIBUTE PROPERTY MATRIX IN GP
C
      CALL DZERO(GP,KZYVAR)
C
      KSYM1 = 0
      DO 1400 IG = 1,KZWOPT
         K = JWOP(1,IG)
         L = JWOP(2,IG)
         KSYM  = ISMO(K)
         LSYM  = ISMO(L)
         IF (KSYM.NE.KSYM1) THEN
            KSYM1 = KSYM
            NORBK = NORB(KSYM)
            IORBK = IORB(KSYM)
            NASHK = NASH(KSYM)
            NISHK = NISH(KSYM)
            IASHK = IASH(KSYM)
            IORBL = IORB(LSYM)
            NASHL = NASH(LSYM)
            NORBL = NORB(LSYM)
            NISHL = NISH(LSYM)
            IASHL = IASH(LSYM)
         END IF
         ITYPK = IOBTYP(K)
         ITYPL = IOBTYP(L)
         IF ( ITYPK.EQ.JTINAC )THEN
            GP(KZCONF+IG) = GP(KZCONF+IG) + D2 * PRPMO(L,K)
            GP(KZVAR+KZCONF+IG) = GP(KZVAR+KZCONF+IG)
     *             - D2 * PRPMO(K,L)
            IF ( ITYPL.EQ.JTACT ) THEN
               NWL = ISW(L) - NISHT
               GP(KZCONF+IG) = GP(KZCONF+IG)
     *         -DDOT(NASHL,PRPMO(IORBL+NISHL+1,K),1,
     *                     UDV(IASHL+1,NWL),1)
               DO 730 IX = 1,NASHL
                  GP(KZVAR+KZCONF+IG) = GP(KZVAR+KZCONF +IG)
     *            +PRPMO(K,IORBL+NISHL+IX)*UDV(NWL,IASHL+IX)
 730           CONTINUE
            ENDIF
         ELSE
           IF (ITYPL.EQ.JTACT) THEN
               NWL = ISW(L) - NISHT
               NWK = ISW(K) - NISHT
               GP(KZCONF+IG) = GP(KZCONF+IG)
     *            -DDOT(NASHL,PRPMO(IORBL+NISHL+1,K),1,
     *                        UDV(IASHL+1,NWL),1)
               DO 740 IX = 1,NASHL
                  GP(KZVAR+KZCONF+IG) = GP(KZVAR+KZCONF +IG)
     *               +PRPMO(K,IORBL+NISHL+IX)*UDV(NWL,IASHL+IX)
 740           CONTINUE
            ELSE
               NWK = ISW(K) - NISHT
            ENDIF
            GP(KZVAR+KZCONF+IG) = GP(KZVAR+KZCONF+ IG)
     *           -DDOT(NASHK,PRPMO(IORBK+NISHK+1,L),1,
     *                 UDV(IASHK+1,NWK),1)
            DO 750 IX = 1,NASHK
               GP(KZCONF+IG) = GP(KZCONF+IG)
     *              +PRPMO(L,IORBK+NISHK+IX)*UDV(NWK,IASHK+IX)
 750        CONTINUE
         ENDIF
 1400 CONTINUE
C
C *** Perform supersymmetry averaging ?
C
      IF (RSPSUP .AND. (KSYMOP .EQ. 1)) THEN
         CALL RSPAVE(GP(KZCONF+1),KZVAR,2)
      END IF
C
      IF (IPRRSP.GT.120) THEN
         WRITE(LUPRI,*)' (PRPORB) Z and Y ORBITAL PROPERTY GP VECTOR'
         CALL OUTPUT(GP,1,KZVAR,1,2,KZVAR,2,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck uthv */
      SUBROUTINE UTHV(U,PRPAO,V,ISYM,JSYM,NBASI,NBASJ,PRPMO,WRK)
C
C WRITTEN 28 JAN 1987 PJ
C
C TRANSFORM (ISYM,JSYM) SYMMETRY BLOCK OF
C THE MATRIX PRPAO FROM AO SYMMETRY ORBITALS
C TO MO BASIS
C
C Input
C  PRPAO ONE-ELECTRON property MATRIX OVER SYMMETRY AO'S
C  U     MO COEFFICIENTS FOR SYMMETRY ISYM
C  V     MO COEFFICIENTS FOR SYMMETRY JSYM
C
C Output
C  PRPMO  (ISYM,JSYM) BLOCK OF THE ONE-ELECTRON PRPAO
C         MATRIX TRANSFORMED TO MO BASIS.
C         PRPMO HAS DIMENSION NORBT*NORBT
C
C Scratch
C  WRK  dimension: NBAS(ISYM)
C
#include <implicit.h>
      DIMENSION U(NBASI,*), V(NBASJ,*)
      DIMENSION PRPAO(NBAST,*),PRPMO(NORBT,*),WRK(*)
C
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
C
      PARAMETER (D0 = 0.0D0)
C
      NORBI=NORB(ISYM)
      IBASI=IBAS(ISYM)
      IORBI=IORB(ISYM)
C
C  TRANSFORM TO MO BASIS
C
      NORBJ=NORB(JSYM)
      IBASJ=IBAS(JSYM)
      IORBJ=IORB(JSYM)
      DO 100 I=1,NORBI
         DO 200 JB=1,NBASJ
            WRK(JB)=DDOT(NBASI,PRPAO(IBASI+1,IBASJ+JB),1,U(1,I),1)
 200     CONTINUE
         DO 300 J=1,NORBJ
            PRPMO(IORBI+I,IORBJ+J)=DDOT(NBASJ,WRK(1),1,V(1,J),1)
 300     CONTINUE
 100  CONTINUE
      RETURN
      END
