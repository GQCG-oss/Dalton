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
c===========================================================================
c/* old abaedrv_log */
c930209-hjaaj: RESMC: defined IPRSOL=IPRINT
c930204-hjaaj: READMC: subtract solvent contribution from FC
c920324-hjaaj: RESINP: increased MAXNR from 20 to 40
c920302-hjaaj: corrected NCONST checks (was .gt. 1, should be .gt. 0)
c900216-hjaaj: reimplementing INGD (not finished)
c-- RESST: check for zero norm of GD vectors
c900219-hjaaj:
c-- RESST: corrected error (removed unnec., erron. ref. to IBNDX)
c      END
c===========================================================================
C  /* Deck resinp */
      SUBROUTINE RESINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxash.h>
      PARAMETER (NTABLE = 16)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <cbires.h>
C
C     Used from common blocks:
C       GNRINF : THR_REDFAC
C       ANRINF : THRNR, MAXNR
C       NUCLEI : NDCORD
C
#include <gnrinf.h>
#include <abainf.h>
#include <anrinf.h>
#include <dorps.h>
#include <infind.h>
#include <nuclei.h>
      DATA TABLE /'.SKIP  ', '.PRINT ','.MAX IT','.THRESH',
     *            '.MAXSIM', '.NEWRD ','.NOTRIA','.RDVECS',
     *            '.NONEXT', '.NRREST','.NOAVER','.D1DIAG',
     *            '.MCHESS', '.DONEXT','.STOP  ','.MAXRED'/
C
      NEWDEF = (WORD .EQ. '*RESPON')
C
      MAXRDM = MAXRED
      MAXSMD = MAXSIM
      MAXNRD = MAXNR
      THRNRD = THRNR
C
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in RESINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in RESINP.')
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,*) MAXNR
                  IF (MAXNR .EQ. MAXNRD) ICHANG = ICHANG - 1
               GO TO 100
    4          CONTINUE
                  READ (LUCMD,*) THRNR
                  IF (THRNR .EQ. THRNRD) ICHANG = ICHANG - 1
               GO TO 100
    5          CONTINUE
                  READ (LUCMD,*) MAXSIM
                  IF (MAXSIM .EQ. MAXSMD)ICHANG = ICHANG - 1
               GO TO 100
    6          CONTINUE
                  NEWRD  = .TRUE.
               GO TO 100
    7             NOTRIA = .TRUE.
               GO TO 100
    8          CONTINUE
                  READ (LUCMD,*) NRDT
                  READ (LUCMD,*) (NRDCO(I), I = 1,NRDT)
               GO TO 100
    9             KAPTST = -1
               GO TO 100
   10             RSTNR = .TRUE.
               GO TO 100
   11             NOAVER = .TRUE.
               GO TO 100
   12             D1DIAG = .TRUE.
               GO TO 100
   13             PRIL2  = .TRUE.
               GO TO 100
   14             KAPTST = 1
               GO TO 100
   15             CUT    = .TRUE.
               GO TO 100
   16             READ (LUCMD,*) MAXRED
                  IF (MAXRED .EQ. MAXRDM)ICHANG = ICHANG - 1
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in RESINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in RESINP.')
            END IF
      END IF
  300 CONTINUE
      IF (THR_REDFAC .GT. 0.0D0) THEN
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         THRNR = THRNR*THR_REDFAC
      END IF
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for RESPON in ABACUS:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' RESPON skipped in this run.'
         ELSE
            IF (IPRINT .NE. IPRDEF) WRITE (LUPRI,'(A,I5)')
     *         ' Print level in RESPON        :',IPRINT
            IF (THRNR .NE. THRNRD) WRITE (LUPRI,'(A,1P,E9.2)')
     *            ' Threshold in RESPON          :',THRNR
            IF (MAXNRD .NE. MAXNR) THEN
               WRITE(LUPRI,'(A,I5)')' Maximum iterations in RESPON :',
     *                              MAXNR
            END IF
            IF (MAXSIM .NE. MAXSMD) THEN
               WRITE(LUPRI,'(A,I5)')' Maximum vectors in ANRCTL    :',
     *                              MAXSIM
            END IF
            IF (MAXRED .GT. LIROW) THEN
               WRITE (LUPRI,'(/A/A,I5,A,I5)')
     &            ' Too large maximum reduced space requested',
     &            ' The requested value of',MAXRED,' is larger '//
     &              'than LIROW parameter in file infind.h:',LIROW
               CALL QUIT('Too large reduced space requested')
            ELSE IF (MAXRED .NE. MAXRDM) THEN
               WRITE(LUPRI,'(A,I5)')
     &           ' Maximum dimension of reduced space (MAXRED):', MAXRED
            END IF
            IF (NEWRD) THEN
               WRITE (LUPRI,'(/A)') ' New RD file will be used.'
            ELSE IF (NOTRIA) THEN
               WRITE (LUPRI,'(/2A)') ' Old solution vectors not',
     *            ' used as trial vectors in RESPON.'
            END IF
            IF (NRDT .GT. 0) THEN
               WRITE (LUPRI,'(/A,(T45,10I3))')
     *            ' Only solve specified gradient vectors:',
     *            (NRDCO(I), I = 1,NRDT)
            END IF
            IF (KAPTST .GT. 0) THEN
               WRITE (LUPRI,'(/A)')
     *         ' Forced optimal orbital trial vectors in this run!'
            ELSE IF (KAPTST .LT. 0) THEN
               WRITE (LUPRI,'(/A)')
     *         ' No optimal orbital trial vectors in this run!'
            END IF
            IF (RSTNR) WRITE (LUPRI,'(/A)')
     *         ' Restart of RESPON from saved trial vectors.'
            IF (D1DIAG) WRITE (LUPRI,'(/2A)')
     *         ' Diagonal Hessian elements not used when',
     *         ' generating trial vectors.'
            IF (NOAVER) WRITE (LUPRI,'(/2A)')
     *         ' Diagonal orbital Hessian for trial vectors',
     *         ' approximated by Fock contributions.'
            IF (PRIL2) WRITE (LUPRI,'(/2A)')
     *         ' The electronic Hessian will be calculated explicitly',
     *         ' and tested for symmetry.'
            IF (CUT) THEN
               WRITE (LUPRI,'(/,A)') ' Program is stopped after RESPON.'
            END IF
         END IF
      END IF
      IF (NRDT .LT. 0 .OR. NRDT .GT. NDCORD(1)+6) THEN
         WRITE (LUPRI,'(A,I10)') ' *** RESINP ERROR: Illegal NRDT:',NRDT
         CALL QUIT(' *** ERROR (RESINP) Illegal NRDT in input')
      END IF
      IPRNR = IPRINT
      RETURN
      END
C  /* Deck resini */
      SUBROUTINE RESINI
C
C     Initialize /CBIRES/
C
#include <implicit.h>
#include <mxcent.h>
#include <anrinf.h>
#include <cbires.h>
#include <abainf.h>
#include <nuclei.h>
C
      IPRINT = IPRDEF
      MAXSIM = 15
      NRDT   = 0
      KAPTST = 0
      SKIP   = .FALSE.
      CUT    = .FALSE.
      NOTRIA = .FALSE.
      NEWRD  = .FALSE.
      RSTNR  = .FALSE.
      NOAVER = .FALSE.
      D1DIAG = .FALSE.
      PRIL2  = .FALSE.
      MAXRED = MAX(400,25*NUCIND)
C     hjaaj Jan 2004: try to allocate sufficient dimension of reduced space
C     for calculation of response contributions to molecular hessian,
C     also for large molecules in C1 symmetry (for NSYM .gt. 1
C     the maximum number of nuclear coordinates in any symmetry must
C     be less than 3*NUCIND).
C
C     Initialize /ANRINF/
C
      THRNR  = 1.D-03
      MAXNR  = 60
C
      THRNRD = THRNR
      MAXNRD = MAXNR
      MAXSMD = MAXSIM
      MAXRDM = MAXRED
      IF (.NOT. (MOLHES .OR. DIPDER .OR. POLAR .OR. QPGRAD))
     &     SKIP = .TRUE.
      RETURN
      END
C  /* Deck respon */
      SUBROUTINE RESPON(WORK,LWORK,PASS)
C
C     Written 23-jan-1985 Hans Joergen Aa. Jensen
C     Modified 14-jun-1985 TUH
C
C     Purpose: Solve response equations
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <iratdef.h>
#include <priunit.h>
      LOGICAL PASS
      DIMENSION WORK(LWORK)
C
C      Used from common blocks:
C        INFORB : NNORBT
C        INFDIM : NVARMA
C        NUCLEI : NUCDEP, DCORD(MXCENT,3)
C        INFTAP : LUGDR,LURDR, LRGDR,LRRDR, NBGDR,NBRDR
C
#include <cbires.h>
#include <abainf.h>
#include <inflin.h>
#include <infvar.h>
#include <infdim.h>
#include <inforb.h>
#include <dorps.h>
#include <nuclei.h>
#include <inftap.h>
C
      LOGICAL OLDDX
C
      IF (SKIP) RETURN
      CALL QENTER('RESPON')
      IF (IPRINT .GT. 0) THEN
         CALL TIMER('START ',TIMEIN,TIMOUT)
         WRITE (LUPRI,'(//A,/)')
     *    '  <<<<<<<<<< Output from RESPON >>>>>>>>>> '
      END IF
C
      CALL GPOPEN(LUGDR,ABAGDR,'OLD','DIRECT',' ',IRAT*NVARMA,OLDDX)
      IF (NEWRD) THEN
         CALL GPOPEN(LURDR,ABARDR,'NEW','DIRECT',' ',IRAT*NVARMA,OLDDX)
      ELSE
         CALL GPOPEN(LURDR,ABARDR,'UNKNOWN','DIRECT',' ',IRAT*NVARMA,
     &               OLDDX)
      END IF
      NEWRD = .NOT. OLDDX
C
      IF (NEWRD) NOTRIA = .TRUE.
      IF (NRDT .EQ. 0) THEN
         NRDT = 0
         IF (MOLHES) THEN
            DO 20 IATOM = 1, NUCIND
               DO 21 ICOOR = 1, 3
                  IF (DCORD(IATOM,ICOOR,1)) THEN
                     NRDT = NRDT + 1
                     NRDCO(NRDT) = 3*(IATOM - 1) + ICOOR
                  END IF
   21          CONTINUE
   20       CONTINUE
         END IF
         IF (POLAR) THEN
            NRDCO(NRDT + 1) = 3*NUCDEP + 1
            NRDCO(NRDT + 2) = 3*NUCDEP + 2
            NRDCO(NRDT + 3) = 3*NUCDEP + 3
            NRDT = NRDT + 3
         END IF
         IF (QPGRAD .AND. .NOT. MOLHES) THEN
            NRDCO(NRDT + 1) = 3*NUCDEP + 4
            NRDCO(NRDT + 2) = 3*NUCDEP + 5
            NRDCO(NRDT + 3) = 3*NUCDEP + 6
            NRDCO(NRDT + 4) = 3*NUCDEP + 7
            NRDCO(NRDT + 5) = 3*NUCDEP + 8
            NRDCO(NRDT + 6) = 3*NUCDEP + 9
            NRDT = NRDT + 6
         END IF
      END IF
      CALL RESMC(.FALSE.,WORK,LWORK)
      CALL GPCLOSE(LURDR,'KEEP')
      CALL GPCLOSE(LUGDR,'KEEP')
      IF (IPRINT .GT. 0) CALL TIMER ('RESPON',TIMEIN,TIMOUT)
      PASS = .TRUE.
      IF (CUT) THEN
         WRITE (LUPRI,'(/,A)')
     &          ' Program stopped after RESPON as required.'
         WRITE (LUPRI,'(A)') ' No restart file has been written.'
         CALL QUIT(' ***** End of ABACUS (in RESPON) *****')
      END IF
      CALL QEXIT('RESPON')
      RETURN
C
      END
C  /* Deck resmc */
      SUBROUTINE RESMC(INGD,WRK,LWRK)
C
C     Purpose:
C
C      Calculate solutions (RD) to the NRDT set of Newton-
C      Raphson equations (A-B)*(RD)-(GD) = 0.
C
C      If INGD is TRUE one set of linear equations are solved.
C      The routine then assumes GD is the first NVAR elements in WRK
C      and returns RD as the first NVAR elements in WRK.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxash.h>
#include <maxorb.h>
#include <ibndxdef.h>
#include <iratdef.h>
      PARAMETER (D1 = 1.0D0)
C
      DIMENSION WRK(*)
      LOGICAL INGD
C
C Used from common blocks:
C   ABAINF : POLAR, DOSYM
C   INFLIN : NCONRF,?
C   INFDIM : LACIMX, LBCIMX, LCINDX, NVARMA,...
C
#include <abainf.h>
#include <inflin.h>
#include <linaba.h>
#include <cbires.h>
#include <cbisol.h>
#include <gdvec.h>
#include <energy.h>
#include <infinp.h>
#include <infvar.h>
#include <inforb.h>
#include <infdim.h>
#include <infind.h>
#include <infpri.h>
#include <inftap.h>
#include <gnrinf.h>
C
      CALL QENTER('RESMC')
C
C     Define for INFLIN and for solvent:
C
      IPRSOL = IPRINT
      IPRLIN = IPRINT
C
      IF (INGD) THEN
         LMAXSI = 1
      ELSE
         LMAXSI = MAXSIM
      END IF
C
C     Define for LINABA:
C
      L2PRI  = PRIL2
      AVEODG = .NOT.NOAVER
C
      IF (IPRLIN.GT.10) WRITE(LUPRI,'(/A,L8)')' INGD in RESMC :',INGD
C
C     Make sure perturbation symmetry is 1 before work space
C     allocations and call readmc (with IPRINT = -1) :
C
      CALL ABAVAR(1,.FALSE.,-1,WRK,LWRK)
C
      IF (SOLVNT) CALL GPOPEN(LUSOL,FNSOL,'OLD',' ','UNFORMATTED',
     &                        IDUMMY,.FALSE.)
C
C     Work space allocations:
C
      KGD    = 1
      IF (INGD) THEN
         KREDH  = KGD    + NVARPT
      ELSE
         KREDH  = 1
      END IF
      KREDGD = KREDH  + MAXRED*(MAXRED + 1)/2
      KSOLEQ = KREDGD + MAXRED*LMAXSI
      KCMO   = KSOLEQ + MAXRED*LMAXSI
      KCREF  = KCMO   + NCMOT
      KGORB  = KCREF  + NCONRF
      KDV    = KGORB  + NWOPH
      KUDV   = KDV    + NNASHX
      KPV    = KUDV   + N2ASHX
      KFOCK  = KPV    + NNASHX*NNASHX
      KFC    = KFOCK  + N2ORBT
      KFV    = KFC    + NNORBT
      KFCAC  = KFV    + NNORBT
      KH2AC  = KFCAC  + NNASHX
      KCINDX = KH2AC  + NNASHX*NNASHX
      KIBNDX = KCINDX + LCINDX
      KWRK1  = KIBNDX + MAXRED
      LWRK1 = LWRK - KWRK1 + 1
      IF (LWRK1 .LE. 0) CALL STOPIT('RESMC','need more than',KWRK1,LWRK)
C
C     If new LURDR, initialize file with zero vectors
C
      IF (.NOT.INGD .AND. NEWRD) THEN
         IF (NVARMA .GT. LWRK) CALL STOPIT('RESMC','LURDR',NVARMA,LWRK)
         CALL DZERO(WRK,NVARMA)
         DO 80 IREC = 1, 2*NGDTOT(1)
            CALL WRITDX (LURDR,IREC,IRAT*NVARMA,WRK)
   80    CONTINUE
         NEWRD = .FALSE.
      END IF
C
C     Read MC information:
C
      CALL READMC(WRK(KCMO),WRK(KCREF),WRK(KDV),WRK(KPV),WRK(KFOCK),
     *            WRK(KFC),WRK(KFV),WRK(KFCAC),WRK(KH2AC),
     *            WRK(KGORB),WRK(KWRK1),LWRK1)
      IF (NASHT.GT.0) CALL DSPTSI(NASHT,WRK(KDV),WRK(KUDV))
C
      IF (IPRLIN.GT.100) THEN
            WRITE(LUPRI,'(/A)')' <<<<< Output from RESMC >>>>> '
         IF (NASHT.GT.0) THEN
            CALL HEADER('One-electron active density matrix packed',-1)
            CALL OUTPAK(WRK(KDV),NASHT,1,LUPRI)
            CALL HEADER('One-electron active density matrix unpacked',
     *                   -1)
            CALL OUTPUT(WRK(KUDV),1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
            CALL HEADER('Fock valence matrix',-1)
            CALL OUTPKB(WRK(KFV),NORB,NSYM,1,LUPRI)
         END IF
         CALL HEADER('Fock core matrix',-1)
         CALL OUTPKB(WRK(KFC),  NORB,NSYM,1,LUPRI)
         CALL HEADER('Total Fock matrix',-1)
         CALL OUTPTB(WRK(KFOCK),NORB,NSYM,1,LUPRI)
      END IF
C
C     Get CI indexing information
C     920217-hjaaj: ICSYM=IHCSYM=LSYMRF means that .TOTSYM
C        always will work for CSF's.  (Both ICSYM and IHCSYM have
C        no effect for determinants: they are overwritten by
C        CALL ABAVAR(ISYM) below.  However, only one CSF <-> det
C        transformation is defined in the current version.)
C
      IF (.NOT. ABAHF) THEN
         CALL GETCIX(WRK(KCINDX),LSYMRF,LSYMRF,WRK(KWRK1),LWRK1,0)
      END IF
C     CALL GETCIX(XNDXCI,ICSYM,IHCSYM,WRK,LWRK,NOSYM)
C
C     ************************************
C     ***** Solve response equations *****
C     ************************************
C
      DO 100 ISYM = 1,NSYM
      IF (DOSYM(ISYM) .AND. (NGDVEC(ISYM,1) .GT. 0)) THEN
         CALL ABAVAR(ISYM,.FALSE.,IPRLIN,WRK(KWRK1),LWRK1)
         IF (IPRLIN .GE. 3) WRITE (LUPRI,'(/A,I5)')
     &      ' Maximum number of response vectors for this symmetry :',
     &       NGDVEC(ISYM,1)
         IF (D1DIAG) THEN
            WRITE (LUPRI,'(A/)') ' Diagonal Hessian set to unity.'
            DO 200 I = 1, NVARPT
               WRK(KWRK1 - 1 + I) = D1
  200       CONTINUE
         ELSE
            CALL ABADIA(WRK(KUDV),WRK(KFOCK),WRK(KFC),WRK(KFV),
     *                  WRK(KFCAC),WRK(KH2AC),WRK(KCINDX),
     *                  WRK(KWRK1),LWRK1)
C           CALL ABADIA(UDV,FOCK,FC,FV,FCAC,H2AC,XNDXCI,WRK,LWRK)
         END IF
         KDIACI = KWRK1
         KDIAOR = KDIACI + NCONST
         KWRK1  = KDIAOR + NWOPPT
         LWRK1  = LWRK   - KWRK1
C
C WORK SPACE REQUIREMENT FOR ABALIN
C  EACH CALL:                        LEACH
C  EACH ORBITAL TRIAL VECTOR:        LORB
C  EACH CONFIGURATION TRIAL VECTOR : LCONF
C
C        CALL LINMEM(NCSIM,NOSIM,KNEED)
         CALL LINMEM(1,0,LC1)
         CALL LINMEM(2,0,LC2)
         CALL LINMEM(0,1,LO1)
         CALL LINMEM(0,2,LO2)
         LCONF  = LC2 - LC1
         LORB   = LO2 - LO1
         LEACH  = MAX(LC1-LCONF,LO1-LORB)
C MAXIMUM FOR EACH TRIAL VECTOR
         LMAX   = MAX(LORB,LCONF)
C GD VECTORS AND MORE FROM ANRCTL
         LMAX   = LMAX + NVARPT + N2ASHX + LMAXSI
C
         LLEFT  = LWRK1 - LEACH - LMAX
         LMAXVE = MIN((LLEFT/NVARPT)+1,LMAXSI)
C        ... LMAXVE is the maximum number of GD vectors we can
C            have in core, with sufficient space for linear
C            transformation of one trial vector
         IF ((LMAXVE .LE. 0) .OR. (LWRK1 .LT. LEACH + LMAX)) THEN
            WRITE(LUPRI,'(/A/3(A,I10))')
     *           ' RESMC: TOO LITTLE WORK SPACE',
     *           ' LMAX+LEACH: ',(LMAX+LEACH),
     *           ' LWRK1: ',LWRK1,' LMAXVE: ',LMAXVE
            CALL QUIT('Work space exceeded in RESMC ')
         END IF
C
         CALL ABALR(LMAXVE,INGD,WRK(KREDH),WRK(KREDGD),WRK(KSOLEQ),
     *              WRK(KIBNDX),WRK(KCMO),WRK(KUDV),WRK(KPV),WRK(KFC),
     *              WRK(KFV),WRK(KFCAC),WRK(KH2AC),WRK(KCREF),
     *              WRK(KGORB),WRK(KDIAOR),WRK(KDIACI),
     *              WRK(KCINDX),WRK(KGD),WRK(KWRK1),LWRK1)
C        CALL ABALR(LMAXSI,INGD,REDH,REDGD,SOLEQ,IBNDX,CMO,UDV,PV,
C    *              FC,FV,FCAC,H2AC,CREF,GORB,DIAOR,DIACI,INDXCI,
C    *              GD,WRK,LWRK)
C
      END IF
 100  CONTINUE
C
      IF (SOLVNT) CALL GPCLOSE(LUSOL,'KEEP')
      CALL QEXIT('RESMC')
      RETURN
      END
C  /* Deck readmc */
      SUBROUTINE READMC(CMO,CREF,DV,PV,FOCK,FC,FV,FCAC,H2AC,GORB,
     *                  WRK,LWRK)
C
C 14-Feb-1985 hjaaj
C l.r. 4-Feb-1993 hjaaj
C
C Purpose:
C   Read that MC information written by WRSIFC in SIRIUS which
C   is needed for RESMC.
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0)
      DIMENSION CMO(*),CREF(*),DV(*),PV(*),FOCK(*),
     *          FC(*),FV(*),FCAC(*),H2AC(*),GORB(*),WRK(*)
C
C Used from common blocks:
C  INFINP : LSYM,FLAG(1:NFLAG)
C
#include <infinp.h>
#include <inflin.h>
#include <inftap.h>
#include <infvar.h>
#include <infdim.h>
#include <infpri.h>
#include <infopt.h>
#include <inforb.h>
C
      LOGICAL FNDLAB
C
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC) EPOT,EMY,EACTIV,EMCSCF
      POTNUC = EPOT
      READ (LUSIFC) NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,
     *            NCDETS,NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT
      CALL READT (LUSIFC,NCMOT,CMO)
      IF (NASHT .GT. 0) THEN
         CALL READT (LUSIFC,NCONF,CREF)
         CALL READT (LUSIFC,NNASHX,DV)
         CALL READT (LUSIFC,N2ORBT,FOCK)
         CALL READT (LUSIFC,NNASHX*NNASHX,PV)
         CALL READT (LUSIFC,NNORBT,FC)
         CALL READT (LUSIFC,NNORBT,FV)
         CALL READT (LUSIFC,NNASHX,FCAC)
         CALL READT (LUSIFC,(NNASHX*NNASHX),H2AC)
      ELSE
         READ (LUSIFC)
         READ (LUSIFC)
         CALL READT (LUSIFC,N2ORBT,FOCK)
         READ (LUSIFC)
         CALL READT (LUSIFC,NNORBT,FC)
         READ (LUSIFC)
         READ (LUSIFC)
         READ (LUSIFC)
         CALL DZERO(FV,NNORBT)
      END IF
C
C     Read GORB (including redundant active-active rotations)
C     for ABALIN.
C
      IF ( NWOPH.GT.0) CALL READT(LUSIFC,NWOPH,GORB)
C
C     If (solvent) read solvent t matrix and subtract from FC
C     (because LINTRN in Sirius expects FC without solvent
C      contributions).
C
      IF (FLAG(16)) THEN
         IF (LWRK .LT. NNORBT) CALL STOPIT('READMC','need',NNORBT,LWRK)
         CALL MOLLAB('SOLVINFO',LUSIFC,LUPRI)
         DO 300 I = 1,6
  300       READ (LUSIFC)
         CALL READT (LUSIFC,NNORBT,WRK)
         DO 310 I = 1,NNORBT
            FC(I) = FC(I) - WRK(I)
  310    CONTINUE
      END IF
      IF (IPRLIN.GT.100) THEN
         WRITE(LUPRI,'(/A)')' <<<<< Output from READMC >>>>> '
         IF (NASHT.GT.0) THEN
            CALL AROUND('One-electron active density matrix packed')
            CALL OUTPAK(DV,NASHT,1,LUPRI)
            CALL AROUND('Fock valence matrix in READMC')
            CALL OUTPKB(FV,NORB,NSYM,1,LUPRI)
         END IF
         CALL AROUND('Fock core matrix in READMC')
         CALL OUTPKB(FC,  NORB,NSYM,1,LUPRI)
         CALL AROUND('Total Fock matrix in READMC')
         CALL OUTPTB(FOCK,NORB,NSYM,1,LUPRI)
      END IF
C
C *** End of READMC
C
      RETURN
      END
C  /* Deck resst */
      SUBROUTINE RESST(INGD,NRSTV,CREF,DIACI,DIAOR,IBNDX,
     *                 WRK,LWRK)
C
C Purpose: Create NRSTV start vectors for solution of
C          linear set of equations.
C
C  Input : NRSTV is equal the total number of gradients.
C          The gradients are split into an orbital and an
C          CSF part which each are used as trial vectors.
C
C          IF INGD the GD vector to use for start guess is at WRK(1)
C          ELSE when NOTRIA is FALSE solution vectors from a previous
C          is used as initial trial vectors (if non-zero norm)
C          together with GD vectors from LUGDR.
C
C  Output: NRSTV is the total number of linearly independent
C          trial vectors created
C
C  Modified Jan 16 1988 tuh -  ENDRD introduced since Alliant
C             sometimes aborts when trying to read DA record
C             after end of file
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxash.h>
#include <maxorb.h>
#include <iratdef.h>
#include <ibndxdef.h>
      PARAMETER ( VECTOL=1.0D-3 , THRLDP=1.0D-20 , DTEST = 1.D-3 )
      PARAMETER ( D0 = 0.0D0 )
      DIMENSION IBNDX(*),CREF(*),DIACI(*),DIAOR(*),WRK(*)
      LOGICAL INGD, FNDTRI, FINDDX, ENDRD
#include <cbires.h>
C
#include <inftap.h>
#include <gdvec.h>
#include <infvar.h>
#include <infind.h>
#include <infpri.h>
#include <inflin.h>
#include <linaba.h>
C
C     Define maximum number of simultaneous vectors
C
      IF (INGD) THEN
         KW1 = 1 + NRSTV*NVARPT
         IF (NRSTV .GT. 1 .OR. NRSTV .LE. 0) THEN
            WRITE (LUPRI,'(/A,I5)')
     *         ' *** ERROR (RESST), INGD = .TRUE. and NRSTV =',NRSTV
            CALL QUIT('ERROR in RESST.')
         END IF
      ELSE
         KW1 = 1
      END IF
      NMAX = MAX(NCONST,NWOPPT)
      NSIM = (LWRK-2*NRSTV-NVARPT-KW1) / NVARPT
      IF (NSIM.LE.0) THEN
         WRITE(LUPRI,'(//A,I5,A,I5)')
     $   ' ***RESST*** not enough work space NSIM =',NSIM,'LWRK=',LWRK
         CALL QUIT('***ERROR, Insufficient space in RESST')
      ENDIF
C
      NCONFX=MAX(4,NCONST)
      NWOPTX=MAX(4,NWOPPT)
C
C     CSF coefficients at expansion point is first record
C
      REWIND LUNR3
      REWIND LUNR5
      IF (LOFFRD.EQ.1) THEN
         NCONRX = MAX(4,NCONRF)
         CALL WRITT(LUNR3,NCONRX,CREF)
         CALL DZERO(WRK(KW1),NMAX)
         IF (NCONST.GT.0) CALL WRITT(LUNR5,NCONFX,WRK(KW1))
         IF (NWOPPT.GT.0) CALL WRITT(LUNR5,NWOPTX,WRK(KW1))
         IBNDX(1) = JBCNDX
         IF (IPRLIN .GE. 4) WRITE(LUPRI,'(/A)')
     &      ' Reference vector written on LUNR3 (LOFFRD = 1)'
      END IF
C
C     Use gradient vector as trial vectors
C
      FNDTRI = .FALSE.
      ENDRD  = .FALSE.
      NCSIM  = 0
      NOSIM  = 0
      DO 4000 ISIM = 1,NRSTV,NSIM
        NBX = MIN(NSIM,(NRSTV+1-ISIM))
        NBC = 0
        NBO = 0
        KBCVEC = 1
        KBOVEC = KBCVEC+NBX*NCONST
        IBCVEC = KBCVEC
        IBOVEC = IBCVEC+NBX*NCONST
        IWRK   = KBCVEC+NBX*NVARPT
        DO 2000 I = 1,NBX
          IREC   = IGDREC(ISIM + I - 1,LSYMPT,1)
          ICOOR  = IGDCOR(ISIM + I - 1,LSYMPT,1)
          IF (.NOT.INGD .AND. .NOT. NOTRIA) THEN
             IF (ENDRD) THEN
                FNDTRI = .FALSE.
             ELSE
                FNDTRI = FINDDX(LURDR,2*IREC-1,NVARPT*IRAT,WRK(IWRK))
             END IF
             IF (.NOT. FNDTRI) THEN
                IF (IPRLIN.GT.0) THEN
                   WRITE (LUPRI,'(2(A,I3),A)')
     &               ' No vector found on LURDR for coordinate ',ICOOR,
     &               ' (record ',IREC,').'
                   WRITE (LUPRI,'(A)') ' End of LURDR encountered. '
                END IF
                ENDRD = .TRUE.
             END IF
          END IF
          IF (FNDTRI) THEN
             VECNR2=DDOT(NVARPT,WRK(IWRK),1,WRK(IWRK),1)
             IF (IPRLIN.GT.5) THEN
                WRITE(LUPRI,'(/A,D15.7)')' VECNR2:',VECNR2
             END IF
             IF(ABS(VECNR2).GT.VECTOL)THEN
                IF (NCONST .GT. 0) THEN
                   CALL DCOPY(NCONST,WRK(IWRK),1,WRK(IBCVEC),1)
                   IBCVEC = IBCVEC + NCONST
                   NBC    = NBC    + 1
                END IF
                IF (NWOPPT .GT. 0) THEN
                   CALL DCOPY(NWOPPT,WRK(IWRK+NCONST),1,WRK(IBOVEC),1)
                   IBOVEC = IBOVEC + NWOPPT
                   NBO    = NBO    + 1
                END IF
                GO TO 2000
             END IF
          END IF
          IF (INGD) THEN
             CALL DCOPY(NVARPT,WRK,1,WRK(IWRK),1)
          ELSE
             CALL READDX(LUGDR,IREC,NVARPT*IRAT,WRK(IWRK))
          END IF
          IF ( IPRLIN.GT.40 ) THEN
             IF (INGD) THEN
                CALL AROUND('Input GD vector in RESST (INGD .true.)')
             ELSE
                CALL AROUND('GD vector in RESST')
                WRITE (LUPRI,'(/A,I5)') ' Coordinate ICOOR',ICOOR
                WRITE (LUPRI,'( A,I5)') ' Record IREC     ',IREC
             END IF
             CALL OUTPUT(WRK(IWRK),1,NVARPT,1,1,NVARPT,1,1,LUPRI)
          END IF
          GDNRM = DNRM2(NVARPT,WRK(IWRK),1)
          IF (GDNRM .EQ. D0) THEN
             IF (INGD) THEN
                WRITE (LUPRI,'(/A)')
     &             ' Input GD vector in RESST zero vector (INGD .true.)'
             ELSE
                WRITE (LUPRI,'(/A,I5,A)')
     &             ' GD vector ICOOR =',ICOOR,' is zero vector in RESST'
                CALL WRITDX(LURDR,2*IREC-1,NVARPT*IRAT,WRK(IWRK))
                CALL WRITDX(LURDR,2*IREC  ,NVARPT*IRAT,WRK(IWRK))
             END IF
             GO TO 2000
          END IF
C
C         CALL DCOPY(NCONF,WRK(IWRK),1,WRK(IBCVEC),1)
C         CALL DCOPY(NWOPT,WRK(IWRK+NCONF),1,WRK(IBOVEC),1)
          IWADD = IWRK - 1
          IF (NCONST.GT.0) THEN
             IBCADD = IBCVEC - 1
             DO 1200 K = 1,NCONST
                D = DIACI(K)
                IF (ABS(D) .LT. DTEST) THEN
                   WRK(IBCADD+K) = WRK(IWADD+K) / SIGN(DTEST,D)
                ELSE
                   WRK(IBCADD+K) = WRK(IWADD+K) / D
                END IF
 1200        CONTINUE
             IBCVEC=IBCVEC+NCONST
             NBC=NBC+1
          END IF
          IWADD = IWADD + NCONST
          IF (NWOPPT .GT. 0) THEN
             IBOADD = IBOVEC - 1
             DO 1400 K = 1,NWOPPT
                D = DIAOR(K)
                IF (ABS(D) .LT. DTEST) THEN
                   WRK(IBOADD+K) = WRK(IWADD+K) / SIGN(DTEST,D)
                ELSE
                   WRK(IBOADD+K) = WRK(IWADD+K) / D
                END IF
 1400        CONTINUE
             IBOVEC=IBOVEC+NWOPPT
             NBO=NBO+1
          END IF
 2000   CONTINUE
        NBPREV = NOSIM + NCSIM + LOFFRD
        IWRKOC = IWRK + NBO + NBC
        IF ( IPRLIN.GT.40 ) THEN
          IF (NBC.GT.0) THEN
            CALL AROUND('CI trial vectors in RESST before ORTBVC')
            WRITE(LUPRI,'(/A,I5)')' NBC ',NBC
            CALL OUTPUT(WRK(KBCVEC),1,NCONST,1,NBC,NCONST,NBC,1,LUPRI)
          END IF
          IF (NBO.GT.0) THEN
            CALL AROUND('Orbital trial vectors in RESST before ORTBVC')
            WRITE(LUPRI,'(/A,I5)')' NBO ',NBO
            CALL OUTPUT(WRK(KBOVEC),1,NWOPPT,1,NBO,NWOPPT,NBO,1,LUPRI)
          END IF
        END IF
        NWOPPM = MAX(NWOPPT,1)
        NCONSM = MAX(NCONST,1)
        CALL ORTBVC(NBC,NBO,NCONSM,NWOPPM,NBPREV,IBNDX,LUNR3,
     &              WRK(KBCVEC),WRK(KBOVEC),THRLDP,WRK(IWRK),
     &              WRK(IWRKOC))
C       CALL ORTBVC(NBC,NBO,NDMBC,NDMBO,NBPREV,IBNDX,LUBVC,
C    $              BCVECS,BOVECS,THRLDP,TT0,OLDVEC)
        IF ( IPRLIN.GT.40 ) THEN
          IF (NBC.GT.0) THEN
            CALL AROUND('CI trial vectors in RESST after ORTBVC')
            WRITE(LUPRI,'(/A,I5)')' NBC ',NBC
            CALL OUTPUT(WRK(KBCVEC),1,NCONST,1,NBC,NCONST,NBC,1,LUPRI)
          END IF
          IF (NBO.GT.0) THEN
            CALL AROUND('Orbital trial vectors in RESST after ORTBVC')
            WRITE(LUPRI,'(/A,I5)')' NBO ',NBO
            CALL OUTPUT(WRK(KBOVEC),1,NWOPPT,1,NBO,NWOPPT,NBO,1,LUPRI)
          END IF
        END IF
        NOSIM=NOSIM+NBO
        NCSIM=NCSIM+NBC
 4000 CONTINUE
      NRSTV=NOSIM+NCSIM
      RETURN
C
C     End of RESST
C
      END
