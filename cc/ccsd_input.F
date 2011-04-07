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
C  /* Deck cc_input */
      SUBROUTINE CC_INPUT(WORD,IREST,MSYM)
C 
C-------------------------------------------------
C
C     30-May 1994 Written by Henrik Koch
C
C     Input sections for the coupled cluster code.
C
C-------------------------------------------------------------
C
C     1994-96 input flags by Ove Christiansen
C     for excitation energies, linear response 
C     and various CC models; 
C     Keywords have been added for R12 method (WK/UniKA/04-11-2002).
C
C     Overview over all keywords:
C
C     Model Keywords:
C     ===============
C     
C     CCS, CC2, CCSD, CC3,
C     CCR(3), CCR(A), CCR(B), CCR(T) (CCSDR()variants)
C     CC(2)(gives CIS(D) excitation energies)
C     CC(3), CC(T)(gives CCSD(T) energy)
C     CC1A, CC1B (for CCSDT-1a and CCSDT-1b models)
C     CCD, MP2
C
C
C     Frozen core and finite diff. Keywords:
C     ======================================
C
C     FROIMP, FROEXP
C     FCORE, FSECON (obsolete)
C     FIELD 
C
C     Control Keywords for energy:
C     ============================
C
C     SKIP, PRINT, DIRECT (the three original)
C     RESTART, NOCCIT, NOT2TC
C     THRENR,  THRLEQ, NSIMLE
C     MAXITE, MXDIIS, MAXRED, MXLRV 
C     MINSCR, MINMEM
C
C-------------------------------------------------------------
C            
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <r12int.h>
#include <maxorb.h>
#include <mxcent.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <inftap.h>
#include <ccfield.h>
#include <cclr.h>
#include <ccfop.h>
#include <leinf.h>
#include <gnrinf.h>
#include <ccrspprp.h>
#include <ccpack.h>
#include <eribuf.h>
#include <cbieri.h>
#include <ccroper.h>
#include <cch2d.h>
Cholesky
#include <cc_cho.h>
#include <ccdeco.h>
#include <chodbg.h>
#include <chomp2.h>
#include <chocc2.h>
C
#include <center.h>
Cholesky  
C
      PARAMETER (NTABLE = 120)
      LOGICAL   SET, NEWDEF, SIRFF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER*(80) LINE
C
      SAVE SET
CSPAS 22.10.2003 keyword .SOPPA(CCSD) added at position 52
      DATA TABLE /'.SKIP  ','.PRINT ','.DIRECT','.RESTAR','.CC3   ',
     &            '*CCEXCI','*CCLRSD','.CCSTST','.NSYM  ','.MAXRED',
     &            '*CCEXGR','.R1SKIP','.L1SKIP','.RESKIP','.LESKIP',
     &            '.F1SKIP','.MAX IT','*CCXOPA','.E0SKIP','.L0SKIP',
     &            '.LISKIP','.CC2   ','.MP2   ','.CC(2) ','.CC1B  ',
     &            '.CC(T) ','.CC(3) ','.CCS   ','.FCORE ','.FSECON',
     &            '.CCD   ','.CC1A  ','.CIS   ','.THRENR','.NOCCIT',
     &            '.IMSKIP','.M1SKIP','.FRSKIP','.MINSCR','.MINMEM',
     &            '.BESKIP','.NEWCAU','.NOT2TC','*CCGR  ','.FROEXP',
     &            '.FROIMP','.MXDIIS','.CCSD  ','.CCR(A)','.CCR(B)',
     &            '*CCFOP ','.SOPPA(','.CCR(3)','.CCR(T)','.FIELD ',
     &            '.DEBUG ','*CCQR2R','.HERDIR','.BUFLEN','*CCLR  ',
     &            '*CCEXLR','.NSIMLE','.THRLEQ','.MXLRV ','*CCTM  ',
     &            '.XXXXXX','*CC5R  ','*CC4R  ','*CCQR  ','*CCCR  ',
     &            '.O2SKIP','.R2SKIP','.X2SKIP','.F2SKIP','.L2SKIP',
     &            '*CCMCD ','.ANAAOD','.PACK  ','.CONNEC','.THRLDP',
     &            '.RCSKIP','.FCSKIP','.LCSKIP','.CO2SKI','.CX2SKI',
     &            '.CR2SKI','.CF2SKI','.CL2SKI','*DERIVA','.N2SKIP',
     &            '.BRSKIP','.FREEZE','*CCSLV ','*R12   ','*R12 IN',
     &            '.PAIRS ','.ETAPTI','.DKABAR','*CCOPA ','*NODDY ',
     &            '.NOEONL','.DIRDER','*CCTPA ','.INT4V ','.ONLYMO',
     &            '.THRVEC','.MTRIP ','.SOPPA2','.AO-SOP','.NOSORT',
     &            '.KEPAOI','*CHO(T)','*CHOCC2','*CHOMP2','*CHODBG',
     &            '.D01SKI','*ACTIVE','.CHO(T)','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
C
      IF (SET) RETURN
      SET = .TRUE.
C
C     Initialize /CCSDINP/ ,/CCLR / and /CCSDSYM/
C
      LGLO = .FALSE.
      SKIP   = .FALSE.
      CCRSTR = .FALSE.
      CCSDT  = .FALSE.
      CC2    = .FALSE.
      MP2    = .FALSE.
      CIS    = .FALSE.
      CCS    = .FALSE.
      CCD    = .FALSE.
      CCP2   = .FALSE.
      CC1A   = .FALSE.
      CC1B   = .FALSE.
      CCPT   = .FALSE.
      CCP3   = .FALSE.
      CCRT   = .FALSE.
      CCR3   = .FALSE.
      CCR1A  = .FALSE.
      CCR1B  = .FALSE.
      RSPIM  = .FALSE.
      TRIPIM = .FALSE.
      LSEC   = .FALSE.
      LCOR   = .FALSE.
      NOCCIT = .FALSE.
      FROIMP = .FALSE.
      FROEXP = .FALSE.
      CCFOP  = .FALSE.
      CCSTST = .FALSE.
      T2TCOR = .TRUE.
      DEBUG  = .FALSE.
      CCLR   = .FALSE.
      CCQR   = .FALSE.
      CCCR   = .FALSE.
      CC4R   = .FALSE.
      CC5R   = .FALSE.
      CCTM   = .FALSE.
      CCMCD  = .FALSE.
      CCLRSD = .FALSE.
      CCQR2R = .FALSE.
      CCOPA  = .FALSE.
      CCXOPA = .FALSE.
      CCTPA  = .FALSE.
      CCEXLR = .FALSE.
      CCDERI = .FALSE.
C     
Cholesky
C 
      CHOPT  = .FALSE.
      CALL CC_CHOPTINIT
      CHODBG = .FALSE.
C     
Cholesky
C 
      LVVVV  = .FALSE.

      L0SKIP = .FALSE.

      F1SKIP = .FALSE.
      R1SKIP = .FALSE.
      L1SKIP = .FALSE.

      RCSKIP = .FALSE.
      FCSKIP = .FALSE.
      LCSKIP = .FALSE.

      O2SKIP = .FALSE.
      R2SKIP = .FALSE.
      X2SKIP = .FALSE.
      F2SKIP = .FALSE.
      L2SKIP = .FALSE.

      CO2SKIP= .FALSE.
      CR2SKIP= .FALSE.
      CX2SKIP= .FALSE.
      CF2SKIP= .FALSE.
      CL2SKIP= .FALSE.

      RESKIP = .FALSE.
      LESKIP = .FALSE.
      E0SKIP = .FALSE.
      IMSKIP = .FALSE.
      LISKIP = .FALSE.
      M1SKIP = .FALSE.
      FRSKIP = .FALSE.
      BESKIP = .FALSE.
      N2SKIP = .FALSE.
      BRSKIP = .FALSE.
      ETADSC = .FALSE.

      FREEZE = .FALSE.
      ONLYMO = .FALSE.
      CCSLV  = .FALSE.
      NCCSLV = 0

      HERDIR = .FALSE.

      NOSORT = .FALSE.

      ANAAOD   = .FALSE.
      MAXRED   = 200
      MXLRV    = MAXRED - 1
      MXDIIS   = 8
      MAXITE   = 40
      NOEONL   = .FALSE.
      THRLEQ   = 1.0D-05
      THRENR   = 1.0D-08
      THRVEC   = 1.0D9
      MTRIP    = .FALSE.
      THRLDPHF = 1.0D-13
CWMK  NFIELD = 0
      NSIMLE = 0
      NEWCAU = .FALSE.
      LBUF   = 250000
C
      KEEPAOTWO = 0
      KEEPAOIN  = .FALSE.
      LPACKINT  = .FALSE.
      THRPCKINT = 1.0D-15
C
CSPAS 22.10.2003 implementing SOPPA(CCSD)
      SIRSOP    = .FALSE.
CKeinSPASmehr
      AOSOPPA   = .FALSE.
C
      CONNECTION = 'SYMMETR'
C
      CALL CC_NODINP(WORD,.TRUE.)
C
      MSYMS = MSYM 
C----------------------------------------------------------------
C Initialize additional labels for CCSD(T)-gradient optimization
C and FOP
C----------------------------------------------------------------

      ETACCPT = .FALSE.
      DIRKAPB = .FALSE.
C
C-----------------------------
C     Cholesky initializations
C-----------------------------
C
      CALL CC_CHOMP2INIT
      CALL CC_CHOCC2INIT
C
C---------------------------------------------------------------------
C     CCSAVE is true for all case, but can be turned of in *ccexci
C     Only used for linear response - mandatory for oscstr and 
C     polarizabilities
C---------------------------------------------------------------------
C
      CCSAVE = .TRUE. 
C
C     Initialize the length of the list of response operator labels
C
      NPRLBL = 0
C
C     Default for evaluation of first derivative integrals
C
      DIRGRD = DIRCAL
C
C     Default values based on gnrinf common block
C
      DIRECT = DIRCAL
      IPRINT = IPRUSR
C
      CALL IZERO(NRHFFR,8)
      CALL IZERO(NVIRFR,8)
      CALL IZERO(KFRRHF,8*MAXFRO)
      CALL IZERO(KFRVIR,8*MAXFRO)
      CALL IZERO(NOMINP,3*8)
C
      ITEST = 0
      DO I=1,8
        NRHFFR(I) = LOCFRO(I)
        ITEST = ITEST + NRHFFR(I)
      END DO
      FROIMP = ITEST .NE. 0
C
      CALL DZERO(EFIELD,MXFELT)
C
      ITEST  = 0
C
      MINSCR = .TRUE.
      MINMEM = .FALSE.
      CCPAIR = .FALSE.
ccn      MKVAJKL = .FALSE.
ckr      R12CAL = .TRUE.
ckr      R12NOA = .FALSE.
ckr      R12NOP = .FALSE.
ckr      R12NOB = .FALSE.
ckr      R12HYB = .TRUE.
ckr      NORXR  = .FALSE.
cwk      VCLTHR =  0D0
cwk      SVDTHR =  1D-15
ckr      R12XXL = .FALSE.

      BRASCL = 1.0D0
      KETSCL = 1.0D0
C
C     If this is a restart run, we read MSYM from SIRIFC
C
      IF (IREST .EQ. 1) THEN
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND LUSIFC
C
         CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
         READ (LUSIFC) MSYM
         CALL GPCLOSE(LUSIFC,'KEEP')
C
      END IF
C-----------------------------------------------
C     SIRIUS values for field is transferred if 
C     there is any.
C-----------------------------------------------
C
      SIRFF = .FALSE.
      NONHF = .FALSE.
      CALL CC_FSIR(MXFELT,NFIELD,LFIELD,EFIELD,NHFFIELD)
      IF (NFIELD .GT. 0) NONHF = .FALSE.
      IF (NFIELD .GT. 0) SIRFF = .TRUE. 
C
C-----------------------------------------------
C  of default section.
C default set after input for minscr and minmem.
C-----------------------------------------------
C
      ICHANG = 0
C
      NEWDEF = (WORD .EQ. '*CC INP' .OR. WORD .EQ. '**CC   '.OR.
     *          WORD .EQ. '*CC    ')
      IF (NEWDEF) THEN
         WORD1 = WORD
 1000    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)

            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 1000
            ELSE IF (PROMPT .EQ. '.' .OR. PROMPT .EQ. '*') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
     *                17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,
     *                32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,
     *                47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,
     *                62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,
     *                77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,
     *                92,93,94,95,96,97,98,99,100,101,102,103,104,
     *                105,106,107,108,109,110,111,112,113,114,115,
     *                116,117,118,119,120), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                GO TO 1000
               END IF
               IF (WORD(1:1) .EQ. '*') THEN
   
                 ! either '*END OF' for '**CC   ' section or a sirius
                 ! keyword. in the former case read next input line
                 IF (WORD.EQ.'*END OF' .AND. WORD1.EQ.'**CC   ') THEN
                   READ (LUCMD, '(A7)') WORD
                   CALL UPCASE(WORD)
                 END IF

                 GO TO 300

               ELSE
                WRITE (LUPRI,'(/3A,/)') ' Keyword "',WORD,
     *             '" not recognized in CCSD_INPUT.'
                CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                CALL QUIT('Illegal keyword in CCSD_INPUT.')
               END IF

    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 1000
    2          CONTINUE
                  READ (LUCMD,*) IPRINT
               GO TO 1000
    3          CONTINUE
                  DIRECT = .TRUE.
               GO TO 1000
    4          CONTINUE
                  CCRSTR = .TRUE.
               GO TO 1000
    5          CONTINUE
                  CC3    = .TRUE.
                  CCSDT  = .TRUE.
               GO TO 1000
    6          CONTINUE
                  CALL CC_EXCINP(WORD,MSYM)
               GO TO 1000
    7          CONTINUE
                  CALL CC_LRSINP(WORD,MSYM)
               GO TO 1000
    8          CONTINUE
                  CCSTST = .TRUE.
               GO TO 1000
    9          CONTINUE
                  READ (LUCMD,*) MSYM2 
                  IF (((MSYM.LT.8).AND.(MSYM.GT.0))
     *               .AND.(MSYM2.NE.MSYM)) THEN
                    CALL QUIT(' Symmetry mismatch in input')
                  ELSE
                     MSYM  = MSYM2
                  ENDIF 
               GO TO 1000
   10          CONTINUE
                  READ (LUCMD, *) MAXRED
               GO TO 1000
   11          CONTINUE
                  CALL CC_EXGRIN(WORD,MSYM)
               GO TO 1000
   12          CONTINUE
                  R1SKIP = .TRUE. 
               GO TO 1000
   13          CONTINUE
                  L1SKIP = .TRUE. 
               GO TO 1000
   14          CONTINUE
                  RESKIP = .TRUE. 
               GO TO 1000
   15          CONTINUE
                  LESKIP = .TRUE. 
               GO TO 1000
   16          CONTINUE
                  F1SKIP = .TRUE.
               GO TO 1000
   17          CONTINUE
                  READ (LUCMD,*) MAXITE
               GO TO 1000
   18          CONTINUE
C                '*CCXOPA'
                  CALL CC_OPAINP(WORD)
               GO TO 1000
   19          CONTINUE
                  E0SKIP = .TRUE.
               GO TO 1000
   20          CONTINUE
                  L0SKIP = .TRUE.
               GO TO 1000
   21          CONTINUE
                  LISKIP = .TRUE.
               GO TO 1000
   22          CONTINUE
                  CC2   = .TRUE.
               GO TO 1000
   23          CONTINUE
                  MP2   = .TRUE.
               GO TO 1000
   24          CONTINUE
                  CCP2  = .TRUE.
               GO TO 1000
   25          CONTINUE
                  CCSDT = .TRUE.
                  CC1B  = .TRUE.
               GO TO 1000
   26          CONTINUE
                  CCPT  = .TRUE.
               GO TO 1000
   27          CONTINUE
                  CCP3  = .TRUE.
               GO TO 1000
   28          CONTINUE
                  CCS   = .TRUE.
               GO TO 1000
   29          CONTINUE
                  LCOR  = .TRUE.
                  READ (LUCMD,*) (ICOR(ISYM),ISYM=1,MSYM)
               GO TO 1000
   30          CONTINUE
                  LSEC  = .TRUE.
                  READ (LUCMD,*) (ISEC(ISYM),ISYM=1,MSYM)
               GO TO 1000
   31          CONTINUE
                  CCD = .TRUE.
               GO TO 1000
   32          CONTINUE
                  CCSDT = .TRUE.
                  CC1A  = .TRUE.
               GO TO 1000
   33          CONTINUE
                  CIS = .TRUE.
               GO TO 1000
   34          CONTINUE
                  READ (LUCMD, *) THRENR
               GO TO 1000
   35          CONTINUE
                  NOCCIT = .TRUE.
               GO TO 1000
   36          CONTINUE
                  IMSKIP = .TRUE.
               GO TO 1000
   37          CONTINUE
                  M1SKIP = .TRUE.
               GO TO 1000
   38          CONTINUE
                  FRSKIP = .TRUE.
               GO TO 1000
   39          CONTINUE
                  READ (LUCMD, *) MINSCR
                  ITEST = ITEST + 1
               GO TO 1000
   40          CONTINUE
                  READ (LUCMD, *) MINMEM
                  ITEST = ITEST + 1
               GO TO 1000
   41          CONTINUE
                  BESKIP = .TRUE.
               GO TO 1000
   42          CONTINUE
                  NEWCAU = .TRUE.
               GO TO 1000
   43          CONTINUE
                  T2TCOR = .FALSE.
               GO TO 1000
   44          CONTINUE
                  CALL CC_GRIN(WORD,MSYM)
               GO TO 1000
   45          CONTINUE
                  FROEXP = .TRUE.
                  IF (FROIMP) FROIMP = .FALSE.
                  IF (FREEZE) CALL QUIT(' Only one of FREEZE - FROEXP')
                  READ(LUCMD,*) (NRHFFR(I),I=1,MSYM)
                  DO 451 ISYM = 1,MSYM
                     IF (NRHFFR(ISYM) .NE. 0) THEN
                        IF (NRHFFR(ISYM) .GT. MAXFRO) THEN
                           WRITE(LUPRI,'(1X,2A,I4)')
     *                          'ERROR: Maximum number of frozen ',
     *                          'orbitals per symmetry is:',MAXFRO
                           CALL QUIT('Too many frozen orbitals')
                        END IF
                        READ(LUCMD,*) (KFRRHF(J,ISYM),J=1,NRHFFR(ISYM))
                     END IF
  451             CONTINUE
                  READ(LUCMD,*) (NVIRFR(I),I=1,MSYM)
                  DO 452 ISYM = 1,MSYM
                     IF (NVIRFR(ISYM) .NE. 0) THEN
                        IF (NVIRFR(ISYM) .GT. MAXFRO) THEN
                           WRITE(LUPRI,'(1X,2A,I4)')
     *                          'ERROR: Maximum number of frozen ',
     *                          'orbitals per symmetry is:',MAXFRO
                           CALL QUIT('Too many frozen orbitals')
                        END IF
                        READ(LUCMD,*) (KFRVIR(J,ISYM),J=1,NVIRFR(ISYM))
                     END IF
  452             CONTINUE
               GO TO 1000
   46          CONTINUE
                  FROIMP = .TRUE.
                  IF (FROEXP) FROEXP = .FALSE.
                  IF (FREEZE) CALL QUIT(' Only one of FREEZE - FROEXP')
                  READ(LUCMD,*) (NRHFFR(I),I=1,MSYM)
                  READ(LUCMD,*) (NVIRFR(I),I=1,MSYM)
               GO TO 1000
   47          CONTINUE
                  READ(LUCMD,*) MXDIIS
               GO TO 1000
   48          CONTINUE
                  CCSD   = .TRUE.
               GO TO 1000
   49          CONTINUE
                  CCR1A  = .TRUE.
               GO TO 1000
   50          CONTINUE
                  CCR1B  = .TRUE.
               GO TO 1000
   51          CONTINUE
                  CALL CC_FOPINP(WORD)
               GO TO 1000
   52          CONTINUE
CSPAS 22.10.2003 implementing SOPPA(CCSD)
                  SIRSOP = .TRUE.
                  CCSD   = .TRUE.
                  KEEPAOTWO = MAX(KEEPAOTWO,2)
CKeinSPASmehr
               GO TO 1000
   53          CONTINUE
                  CCR3   = .TRUE.
                  CCSD   = .TRUE.
               GO TO 1000
   54          CONTINUE
                  CCRT   = .TRUE.
               GO TO 1000
   55          CONTINUE
                  NFIELD = NFIELD + 1
                  IF (NFIELD .LE. MXFELT) THEN
                     READ(LUCMD,*) EFIELD(NFIELD)
                     READ(LUCMD,*) LFIELD(NFIELD)
                     NHFFIELD(NFIELD) = .TRUE.
                  ELSE
                     WRITE(LUPRI,*) 'Too many fields in cc input'
                     CALL QUIT('Too many fields !')
                  ENDIF
                  NONHF = .TRUE.
                  IF (SIRFF ) THEN
                    WRITE(LUPRI,*) ' FF not allowed in '
     *                         //'both Hartree Fock and CC input'
                    CALL QUIT('FF not allowed in both HF and CC input' )
                  ENDIF
               GO TO 1000
   56          CONTINUE
                  DEBUG = .TRUE.
               GO TO 1000
   57          CONTINUE
                  CALL CC_QR2RINP(WORD)
               GO TO 1000
   58          CONTINUE
                  HERDIR = .TRUE.
               GO TO 1000
   59          CONTINUE
                  READ (LUCMD, *) LBUF
               GO TO 1000
   60          CONTINUE
                  CALL CC_LRINP(WORD)
               GO TO 1000
   61          CONTINUE
                  CALL CC_EXLRINP(WORD)
               GO TO 1000
   62          CONTINUE
                  READ (LUCMD, *) NSIMLE
               GO TO 1000
   63          CONTINUE
                  READ (LUCMD, *) THRLEQ
               GO TO 1000
   64          CONTINUE
                  READ (LUCMD, *) MXLRV
               GO TO 1000
   65          CONTINUE
                  CALL CC_TMINP(WORD)
               GO TO 1000
   66          CONTINUE
               GO TO 1000
   67          CONTINUE
                  CALL CC_5RINP(WORD)
               GO TO 1000
   68          CONTINUE
                  CALL CC_4RINP(WORD)
               GO TO 1000
   69          CONTINUE
                  CALL CC_QRINP(WORD)
               GO TO 1000
   70          CONTINUE
                  CALL CC_CRINP(WORD)
               GO TO 1000
   71          CONTINUE
                  O2SKIP = .TRUE.
               GO TO 1000
   72          CONTINUE
                  R2SKIP = .TRUE.
               GO TO 1000
   73          CONTINUE
                  X2SKIP = .TRUE.
               GO TO 1000
   74          CONTINUE
                  F2SKIP = .TRUE.
               GO TO 1000
   75          CONTINUE
                  L2SKIP = .TRUE.
               GO TO 1000
   76          CONTINUE
                 CALL CC_MCDINP(WORD)
               GO TO 1000
   77          CONTINUE
                  ANAAOD = .TRUE.
               GO TO 1000
   78          CONTINUE
C                '.PACK  ' 
                  LPACKINT = .TRUE.
                  READ (LUCMD, *) THRPCKINT
               GO TO 1000
   79          CONTINUE
C                '.CONNEC' 
c                 CONNECTION = 'SYMMETR' / 'NATURAL' 
                  READ (LUCMD, '(A7)') CONNECTION
               GO TO 1000
   80          CONTINUE
C                '.THRLDP' 
                  READ (LUCMD, *) THRLDPHF
               GO TO 1000
   81          CONTINUE
C              '.RCSKIP' 
                 RCSKIP = .TRUE.
               GO TO 1000
   82          CONTINUE
C              '.FCSKIP' 
                 FCSKIP = .TRUE.
               GO TO 1000
   83          CONTINUE
C              '.LCSKIP' 
                 LCSKIP = .TRUE.
               GO TO 1000
   84          CONTINUE
C              '.CO2SKI' 
                 CO2SKIP = .TRUE.
               GO TO 1000
   85          CONTINUE
C              '.CX2SKI' 
                 CX2SKIP = .TRUE.
               GO TO 1000
   86          CONTINUE
C              '.CR2SKI' 
                 CR2SKIP = .TRUE.
               GO TO 1000
   87          CONTINUE
C              '.CF2SKI' 
                 CF2SKIP = .TRUE.
               GO TO 1000
   88          CONTINUE
C              '.CL2SKI' 
                 CL2SKIP = .TRUE.
               GO TO 1000
   89          CONTINUE
C                '*DERIVA' 
                 CCDERI = .TRUE.
                 RELORB = .TRUE.
               GO TO 1000
   90          CONTINUE
C                '.N2SKIP'
                  N2SKIP = .TRUE.
               GO TO 1000
   91          CONTINUE
C                '.BRSKIP'
                  BRSKIP = .TRUE.
               GO TO 1000
   92          CONTINUE
C                '.FREEZE'
                  FREEZE = .TRUE.
                  IF (FROIMP.OR.FROEXP)
     *            CALL QUIT(' Only one of FREEZE - FROEXP - FROIMP')
                  READ(LUCMD,*) NFC,NFV
                  FROIMP = .TRUE.
               GO TO 1000
   93          CONTINUE
C                '*CCSLV '
                  CALL CC_SLVINP(WORD)
               GO TO 1000
   94          CONTINUE
   95          CONTINUE
C                '*R12   ' OR '*R12 IN'
                 CALL CC_R12IN(WORD)
               GO TO 1000
   96          CONTINUE
C                '.PAIRS '
                 CCPAIR = .TRUE.
               GO TO 1000
   97          CONTINUE
C                '.ETAPTI', extra integrals for CCSD(T) geopt (redundant!)
                  ETACCPT = .TRUE.
               GO TO 1000
   98          CONTINUE
C                '.DKABAR', direct KappaBar calculation in nondir CC
                  DIRKAPB = .TRUE.
               GO TO 1000
   99          CONTINUE
C                '*CCOPA ' one-photon absorption strengths
                  CALL CC_OPAINP(WORD)
               GO TO 1000
  100          CONTINUE
C                '*NODDY '
                  CALL CC_NODINP(WORD,.FALSE.)
               GO TO 1000
 101           CONTINUE
                  NOEONL = .TRUE.
               GO TO 1000
 102           CONTINUE
C                 '.DIRDER' direct calculation of derivative integrals
                  DIRGRD = .TRUE. 
               GO TO 1000
 103           CONTINUE
C                 '.CCTPA ' two-photon absorption strengths
                  CALL CC_OPAINP(WORD)
               GO TO 1000
 104           CONTINUE
C                 '.INT4V ' use VVVV integrals in CC3 left transformation
                  LVVVV = .TRUE.
               GO TO 1000
 105           CONTINUE
                  ONLYMO = .TRUE.
               GO TO 1000
 106           CONTINUE
C                 '.THRVEC' convergence threshold for norm of vector function
                  READ(LUCMD,*) RDTHVC
                  IF (RDTHVC .GT. 0.0D0) THRVEC = RDTHVC
               GO TO 1000
 107           CONTINUE
                  MTRIP  = .TRUE.
               GO TO 1000
 108           CONTINUE
                  SIRSOP = .TRUE.
                  CC2   = .TRUE.
                  KEEPAOTWO = MAX(KEEPAOTWO,2)
               GO TO 1000
 109           CONTINUE
C                '.AO-SOPPA'
                  AOSOPPA  = .TRUE.
                  KEEPAOIN = .TRUE.
               GO TO 1000
 110           CONTINUE 
C                 'NOSORT'
                  NOSORT = .TRUE.
               GO TO 1000
 111           CONTINUE 
C                 'KEPAOI'
                  KEEPAOIN = .TRUE.
               GO TO 1000
 112           CONTINUE 
C                '*CHO(T)'
                  CHOPT = .TRUE.
                  CCPT  = .TRUE.
                  CALL CC_CHOPTINP(WORD)
               GO TO 1000
 113           CONTINUE 
C                '*CHOCC2'
                 CALL CC_CHOCC2INP(WORD)
               GO TO 1000
 114           CONTINUE 
C                '*CHOMP2'
                 CALL CC_CHOMP2INP(WORD)
               GO TO 1000
 115           CONTINUE 
C                '*CHODBG'
                 CALL CC_CHODBINP(WORD)
               GO TO 1000
 116           CONTINUE 
C                '.D01SKI'
                 D01SKIP = .TRUE.
               GO TO 1000
 117           CONTINUE 
C                '*ACTIVE'
                 ACTIVE = .TRUE.
                 CALL CC_ACTINP(WORD,MSYM)
               GO TO 1000
 118           CONTINUE 
C                '.CHO(T)'
                  CHOPT = .TRUE.
                  CCPT  = .TRUE.
               GO TO 1000
 119           CONTINUE 
 120           CONTINUE 
               GO TO 1000
            ELSE
               WRITE (LUPRI,'(/3A,/)') ' Prompt "',WORD,
     *            '" not recognized in CC2INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in CC2INP.')
            END IF
      END IF
  300 CONTINUE
C
C---------------------------------------------------
C     set some defaults...
C---------------------------------------------------
C
      MSYM = MSYMS
      IF (ITEST .EQ. 0 ) THEN
         MINSCR = .TRUE.
         IF (DIRECT) MINSCR = .FALSE.
         MINMEM = .FALSE.
         IF (DIRECT) MINMEM = .FALSE.
      ENDIF
      IF ( .NOT. MINSCR ) MINMEM = .TRUE.

      DIRGRD = ( DIRGRD .OR. DIRECT )
       
      !Sonia: replace ETAPTI keyword....
      ETACCPT = ( (OPTNEW.OR.CCDERI).AND.(CCPT) )

Cho
      IF (CHOINT) IPRINT = MAX(IPRINT,1)
      IF (CHOINT) THRVEC = THRENR*1.0D2
Cho
      RETURN
      END
C  /* Deck CC_PRTINP */
      SUBROUTINE CC_PRTINP(IWUNIT)
C
C     K.Ruud, Jan.-96. Split from CC_INPUT in order to place CC output more 
C     adequatly in SIRIUS Print of input processing
C
#include <implicit.h>
#include <priunit.h>
#include <cclrinf.h>
#include <cclr.h>
#include <ccfop.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <cclres.h>
#include <ccqr2r.h>
#include <ccfield.h>
#include <ccsdsym.h>
C
      IF (ICHANG .GT. 0) THEN
         IF (IWUNIT .eq. LUPRI) THEN
            CALL HEADER('Changes of defaults for CC:',0)
         ELSE
            WRITE(IWUNIT,'(//10X,A/10X,A/)')
     &         'Changes of defaults for CC:',
     &         '---------------------------'
         END IF
         IF (SKIP) THEN
            WRITE (IWUNIT,'(A)') ' -CCSD skipped in this run.'
         ELSE
C
            IF (DIRECT) WRITE (IWUNIT,'(/A/A)')
     *         ' -This is a direct atomic orbital integral based',
     *         '  calculation for coupled cluster wave function'
C
            IF (CCSDT) WRITE (IWUNIT,'(/A)')
     *         ' -Iterative triple excitations included '
C
            IF (FROIMP) WRITE(IWUNIT,'(/A)')
     *         ' -Implicit frozen core calculation'
C
            IF (FROEXP) WRITE(IWUNIT,'(/A)')
     *         ' -Explicit frozen core calculation'
C
            IF (.NOT. T2TCOR) WRITE(IWUNIT,'(/A)')
     *         ' -Transposed t2-amplitudes not hold in core'
C
            IF (NFIELD.GT.0) THEN
               WRITE(IWUNIT,'(A)') ' -Calculation with finite fields:'
               DO IF = 1, NFIELD
                  WRITE(IWUNIT,'(A,F10.6,A,A8)')
     *               '    Field strength: ',EFIELD(IF),
     *               'Field Label:  ',LFIELD(IF)
               END DO
               CALL CC_FIELD_PRTINP(IWUNIT)
               IF (CCSDT.AND.NONHF) THEN
                  IF (.NOT. CC3) THEN
                    WRITE(IWUNIT,*)'No triples unrelaxed FF possibility'
                    CALL QUIT('No triples unrelaxed FF possibility')
                  ENDIF
               ENDIF
               IF (CCS.AND.NONHF) THEN
                  WRITE(IWUNIT,*) 'No CCS unrelaxed FF possibility '
                  WRITE(IWUNIT,*) 'Use instead CC2 with CCSTST option.'
                  CALL QUIT('No CCS unrelaxed FF possibility')
               ENDIF                                             
               IF (NONHF .AND. RELORB) THEN
                  WRITE(IWUNIT,*) 'Inconsistency: Non HF reference and '
     *                        //'relaxed derivative requested'
                  CALL QUIT('Inconsistency: in FF '//
     &                 'and relaxed derivative')
               ENDIF
               ! put operators for "unrelaxed" fields on common CCRSPOP
               ! (needed for CCR12 with unrelaxed finite fields)
               IF (NONHF) THEN
                  DO IFIELD = 1, NFIELD
                     IDX = INDPRPCC(LFIELD(IFIELD))
                  END DO
               END IF
            END IF 
C
            IF (LHTR.AND.((CCLRSD).OR.
     *         (CCR3.OR.CCRT.OR.CCR1A.OR.CC1B.OR.CC1A.OR.CC1B)))
     *      THEN
               WRITE(IWUNIT,*) 'Input inconsistent due to LHTR '
               CALL QUIT('Do not use LHTR for this '//
     &              'type of calculation ')
            ENDIF
C
C
            IF ( DEBUG ) WRITE(IWUNIT,'(A)')
     *         ' -Debug printout activated '
            IF ( CCEXCI ) WRITE(IWUNIT,'(A)')
     *         ' -Excitation energies calculated'
            IF ( CCLRSD .OR. CCOPA) WRITE(IWUNIT,'(A)')
     *         ' -One-photon absorption strengths will be calculated'
            IF ( CCTPA ) WRITE(IWUNIT,'(A)')
     *         ' -Two-photon absorption strengths will be calculated'
            IF ( JACTST ) WRITE(IWUNIT,'(A)')
     *         ' -Jacobian tested agains finite difference Jacobian'
            IF ( JACEXP ) WRITE(IWUNIT,'(A)')
     *         ' -Jacobian constructed explicit'
            IF ( FDEXCI ) WRITE(IWUNIT,'(A)')
     *         ' -Excitation energies of finite diff. Jacobian calc.'
            IF ( CCLR  ) WRITE(IWUNIT,'(A)')
     *         ' -Linear response properties calculated'
            IF ( CAUCHY) WRITE(IWUNIT,'(A)')
     *         ' -Dispersion coefficients for linear response calc.'
            IF ( CCQR  ) WRITE(IWUNIT,'(A)')
     *         ' -Quadratic response properties calculated'
            IF ( CCCR  ) WRITE(IWUNIT,'(A)')
     *         ' -Cubic response properties calculated'
            IF ( OSCSTR) WRITE(IWUNIT,'(A)')
     *         ' -Oscillator strengths calculated'
            IF ( CCQR2R .OR. CCXOPA) WRITE(IWUNIT,'(A)')
     *         ' -Transition strengths between two excited states '//
     *         'calculated.'
            IF (CCEXGR) WRITE(IWUNIT,'(A)')
     *         ' -Excited state properties calculated'
            IF ( CCMCD  ) WRITE(IWUNIT,'(A)')
     *         ' -Magnetic circular dichroism B calculated'
            IF (DIPMOM) WRITE(IWUNIT,'(A)')
     *         ' -Dipole moment calculated'
            IF (QUADRU) WRITE(IWUNIT,'(A)')
     *         ' -Traceless quadrupole moment calculated'
            IF (NQCC) WRITE(IWUNIT,'(A)')
     *         ' -Electric field gradient calculated'
            IF (RELCOR) WRITE(IWUNIT,'(A)')
     *         ' -Relativistic corrections to energy calculated'
            IF (SECMOM) WRITE(IWUNIT,'(A)')
     *         ' -Electronic second moment of charge calculated'
            IF (DAR2EL) WRITE(IWUNIT,'(A)')
     *         ' -Relativistic two-electron Darwin term calculated'
            IF (DPTECO) WRITE(IWUNIT,'(A)')
     *         ' -First-order DPT energy corrections calculated'
            IF (SIRSOP .AND. CCSD) WRITE (IWUNIT,'(A)')
     *         ' -CCSD Amplitudes appended to Sirius interface'//
     *         ' for SOPPA(CCSD)'
            IF (SIRSOP .AND. CC2) WRITE (IWUNIT,'(A)')
     *         ' -CC2 Amplitudes appended to Sirius interface'//
     *         ' for SOPPA(CC2)'
            IF (AOSOPPA) WRITE (IWUNIT,'(/A,A)')
     *         ' MP2 Amplitudes written for atom integral direct',
     *         ' SOPPA calculations'
C
         END IF
         WRITE (IWUNIT,'(A)') '  '
      END IF
C
      RETURN
      END
c /* deck CC_FIELD_PRTINP */
      SUBROUTINE CC_FIELD_PRTINP(IWUNIT)
C
C     Calculate nuclear contribution to energy in electric field
C
C     The dipole moment origin is the center of charge.
C     It is assumed that the molecule is properly oriented.
C
C     ASM & JCh  February 1996
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <nuclei.h>
#include <symmet.h>
#include <ccfield.h>
#include <ccorb.h>
C
      PARAMETER (ZERO = 0.0D0)
C
      DIMENSION GEOM(3,MXCENT), QCHAR(MXCENT), ELEFLD(3)
C
      CHARACTER*6 FLDTYP
      CHARACTER*1 FLDDIR
C
#include <ibtfun.h>
C
C----------------------------------------------
C     Calculate total electric field.
C----------------------------------------------
C
      CALL DZERO(ELEFLD,3)
C
      DO 200 I = 1,NFIELD
C
         FLDTYP = LFIELD(I)(2:7)
         FLDDIR = LFIELD(I)(1:1)
C
         IF (FLDTYP .EQ. 'DIPLEN') THEN
C
            IF (FLDDIR .EQ. 'X') THEN
               JDIR = 1
            ELSE IF (FLDDIR .EQ. 'Y') THEN
               JDIR = 2
            ELSE
               JDIR = 3
            END IF
C
            ELEFLD(JDIR) = ELEFLD(JDIR) + EFIELD(I)
C
         END IF
C
  200 CONTINUE
C
      ELFLNR = DSQRT(DDOT(3,ELEFLD,1,ELEFLD,1))
C
      IF (ELFLNR .NE. 0.0D0) THEN
C
         WRITE(IWUNIT,'(/A,3F14.8/A,F14.8)')
     &      ' Electric field:  ',(ELEFLD(I),I=1,3),
     &      ' Total norm:      ',ELFLNR
C
      END IF
C
C----------------------------------------------
C     Cartesian coordinates of dependent atoms.
C----------------------------------------------
C
      JATOM = 0
      DO 300 ICENT = 1, NUCIND
C
         MULCNT = ISTBNU(ICENT)
C
         IF (MULT(MULCNT) .EQ. 1) THEN
C
            JATOM = JATOM + 1
C
            QCHAR(JATOM) = CHARGE(ICENT)
C
            DO 310 I = 1,3
               GEOM(I,JATOM) = CORD(I,ICENT)
  310       CONTINUE
C
         ELSE
C
            DO 320 ISYOPR = 0,MAXOPR
               IF (IBTAND(ISYOPR,MULCNT) .EQ. 0) THEN
C
                  JATOM = JATOM + 1
C
                  QCHAR(JATOM) = CHARGE(ICENT)
C
                  DO 330 I = 1,3
C
                     PTAT          = PT(IBTAND(ISYMAX(I,1),ISYOPR))
                     GEOM(I,JATOM) = PTAT*CORD(I,ICENT)
C
  330             CONTINUE
C
               END IF
  320       CONTINUE
C
         END IF
C
  300 CONTINUE
C
C-----------------------------------------
C     Coordinates of the center of charge.
C-----------------------------------------
C
      XCQ  = ZERO
      YCQ  = ZERO
      ZCQ  = ZERO
      SUMQ = ZERO
C
      DO 400 I = 1,NUCDEP
C
         XCQ  = XCQ  + GEOM(1,I)*QCHAR(I)
         YCQ  = YCQ  + GEOM(2,I)*QCHAR(I)
         ZCQ  = ZCQ  + GEOM(3,I)*QCHAR(I)
         SUMQ = SUMQ + QCHAR(I)
C
  400 CONTINUE
C
      CORR = -(XCQ*ELEFLD(1) + YCQ*ELEFLD(2) + ZCQ*ELEFLD(3))  
C
      XCQ = XCQ/SUMQ
      YCQ = YCQ/SUMQ
      ZCQ = ZCQ/SUMQ
C
C-----------------------------------------------
C     Contribution relative to center of charge.
C-----------------------------------------------
C
C     QTOT = DFLOAT(NRHFTS*2)
C     CORR = -QTOT*(XCQ*ELEFLD(1) + YCQ*ELEFLD(2) + ZCQ*ELEFLD(3))
C
      WRITE(IWUNIT,'(A,3F14.8)') 'Center of charge:',XCQ,YCQ,ZCQ
      WRITE(IWUNIT,'(2A,F14.8)') 'Charge correction to interaction ',
     *                      'with electric field:', CORR
      WRITE(IWUNIT,*)
C
      RETURN
      END
      SUBROUTINE CC_FSIR(MXFEL,NFIEL,LFIEL,EFIEL,NHFFIELD)
C
C     If field was set in Hartree-Fock transfer to CC.
C
C     Ove Christiansen 11-6-1996
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <infinp.h>
C
      DIMENSION EFIEL(MXFEL)
      CHARACTER*8 LFIEL(MXFEL)
      LOGICAL NHFFIELD(MXFEL)
C
      IF (NFIELD .GT. 0 ) THEN
        IF (NFIELD.GT.MXFEL) THEN
           WRITE (LUPRI,*) 
     *           'CC_FSIR: Too many fields added in Hartree Fock.'
           CALL QUIT('CC_FSIR: Too many fields added in Hartree Fock.')
        END IF
        NFIEL = NFIELD
        DO IF = 1, NFIELD
           LFIEL(IF)    = LFIELD(IF)
           EFIEL(IF)    = EFIELD(IF)
           NHFFIELD(IF) = .FALSE.
        END DO
      ENDIF
C
      RETURN
      END
c /* deck cc_excinp */
C=====================================================================*
       SUBROUTINE CC_EXCINP(WORD,MSYM)
C---------------------------------------------------------------------*
C
C    Purpose: Read input for CC excited state calculations.
C
C    if (WORD .eq '*CCEXCI ') read & process input and set defaults, 
C    else set only defaults 
C 
C    Ove Christiansen 24-10 1996
C    Kasper Hald & Christof Haettig 12-08-99, changes for triplet
C
C=====================================================================*
#include <implicit.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccsdsym.h>
#include <cclr.h>
#include <cclres.h>
#include <leinf.h>
#include <cclrinf.h>
#include <ccrspprp.h>
#include <ccexci.h>

* local parameters:
      CHARACTER SECNAM*(9)
      PARAMETER (SECNAM='CC_EXCINP')

      LOGICAL LSTVEC
      INTEGER NTABLE
      PARAMETER (NTABLE = 25)

      DIMENSION NSTAR(8)

* variables:
      LOGICAL SET
      SAVE SET

      CHARACTER WORD*(7)
      CHARACTER TABLE(NTABLE)*(8)

      INTEGER IJUMP, NTRIP 

* data:
      DATA SET /.FALSE./
      DATA TABLE /'.NCCEXC','.R3DIIS','.FDJAC ','.FDEXCI','.JACEXP',
     *            '.JACTST','.LHTR  ','.NOSCOM','.STSD  ','.TOLSC ',
     *            '.OMEINP','.STVEC ','.STOLD ','.CCTREN','.THREXC',
     *            '.CCSPIC','.CC2PIC','.CCSDPI','.MARGIN','.SQROVL',
     *            '.ANALYS','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (SET) RETURN
      SET = .TRUE.

*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*
C
      CCSDT_DIIS = .FALSE.
      FDJAC  = .FALSE.
      FDEXCI = .FALSE.
      JACEXP = .FALSE.
      JACTST = .FALSE.
      LHTR   = .FALSE.
      OMESC  = .TRUE.
      STSD   = .FALSE.
      TOLSC  = 1.0D-04 
      STVEC  = .FALSE. 
      STOLD  = .FALSE.
      THREXC = 1.0D-04
      CCSPIC = .FALSE.
      CC2PIC = .FALSE.
      CCSDPI = .FALSE.
      OMPCCS = 0.0D0
      OMPCC2 = 0.0D0
      OMPCCSD= 0.0D0
      MARGIN = .FALSE.
      EXCI_CONT = .FALSE.
      XMARGIN = 1.0
C
      SQROVLP    = .FALSE.
      CCSDTRENRM = .FALSE.
C
      CALL IZERO(NCCEXCI,3*8)
      CALL IZERO(NOMINP,3*8)
C
C     Other initializations
C
      NSIDE  = 1
      STCCS  = .FALSE.
C
*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCEXCI') THEN

100   CONTINUE
       
* get new input line:
        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
          CALL UPCASE(WORD)
        END DO

        IF (WORD(1:1) .EQ. '.') THEN

c         table look up:
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO
     
c         jump to the appropiate input section:
          IF (IJUMP .LE. NTABLE) THEN
            ICHANG = ICHANG + 1
            GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     *            21,22,23,24,25),IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_EXCINP.')
      
C           -------------------------------------------
C           .NCCEXC # excitation energies to solve for:
C           -------------------------------------------
1           CONTINUE
               ! read # singlett states
               READ (LUCMD,*) (NCCEXCI(ISYM,1),ISYM=1,MSYM)

               WRITE (LUPRI,*) 'NCCEXCI for singlet:',
     &              (NCCEXCI(ISYM,1),ISYM=1,MSYM)

               ! check for further excitation energy input:
               READ (LUCMD,'(A7)') WORD
               CALL UPCASE(WORD)
               BACKSPACE(LUCMD)
               IF (WORD(1:1).NE.'.' .AND. WORD(1:1).NE.'*'
     &                              .AND. WORD(1:1).NE.'!' ) THEN

                  ! read # triplett states
                  READ (LUCMD,*) (NCCEXCI(ISYM,3),ISYM=1,MSYM)
                  WRITE (LUPRI,*) 'NCCEXCI for triplet:',
     &                 (NCCEXCI(ISYM,3),ISYM=1,MSYM)
               END IF                                 
            GO TO 100

C
C---------------------------------------------------------------
C           .R3DIIS  use CCDIIS_SOL for iterative triples models
C---------------------------------------------------------------
C
2           CONTINUE
            CCSDT_DIIS = .TRUE.
            GO TO 100

C
C--------------------------------------------------------
C           .FDJAC Calculate Finited difference jacobian:
C--------------------------------------------------------
C
3           CONTINUE
               FDJAC  = .TRUE.
            GO TO 100

C
C----------------------------------------------------------
C           .FDEXCI Diagonalize finite difference jacobian:
C----------------------------------------------------------
C
4           CONTINUE
              FDEXCI = .TRUE.
            GO TO 100

C
C-------------------------------------------------
C           .JACEXP : Construct jacobian explicit:
C-------------------------------------------------
C
5           CONTINUE
               JACEXP = .TRUE.
            GO TO 100
C
C-------------------------
C           .JACTST : Jacobian test.
C-------------------------
C
6           CONTINUE
               JACTST = .TRUE. 
            GO TO 100
C
C           ---------------------------------------------------
C           .LHTR : Use left hand transformation in calculation
C                   of excitation energies.
C           ---------------------------------------------------     
C
7           CONTINUE
               LHTR   = .TRUE. 
            GO TO 100
C
C           ---------------------------------------------
C           .NOSCOM : Do not solve self-consistently for
C                     triples excitation energies.
C           ---------------------------------------------       
8           CONTINUE
               OMESC  = .FALSE.
            GO TO 100
C
C           ----------------------------------------------------
C           .STSD  Start with calculation of singles and doubles
C                  excitation energies with triples amplitudes.
C           ----------------------------------------------------        
9           CONTINUE
               STSD   = .TRUE. 
            GO TO 100
C
C--------------------------------------------------------------
C           .TOLSC : Set threshold for solving selfconsitently.
C--------------------------------------------------------------
C
10          CONTINUE
               READ (LUCMD, *) TOLSC
            GO TO 100
C
C----------------------------------------------------------
C           .OMEINP : Readin omega for triples calculation.  
C----------------------------------------------------------
C
11          CONTINUE
C              Read the singlet states
               READ (LUCMD,*) (NOMINP(ISYM,1),ISYM=1,MSYM)
               OMEINP = .TRUE.
               DO 131 ISYM = 1, MSYM
                  DO 132 IOM = 1, NOMINP(ISYM,1)
                    READ (LUCMD,*) IOMINP(IOM,ISYM,1),
     *                             EOMINP(IOM,ISYM,1)
  132             CONTINUE
  131          CONTINUE
C
C              Check for further excitation energy input:
               READ(LUCMD,'(A7)') WORD
               CALL UPCASE(WORD)
               BACKSPACE(LUCMD)
               IF (WORD(1:1).NE.'.' .AND. WORD(1:1).NE.'*'
     *            .AND. WORD(1:1).NE.'#' .AND. WORD(1:1).NE.'!' ) THEN
C
C                 Readin for the triplet states
                  READ (LUCMD,*) (NOMINP(ISYM,3),ISYM=1,MSYM)
                  DO 133 ISYM = 1, MSYM
                     DO 134 IOM = 1, NOMINP(ISYM,3)
                       READ (LUCMD,*) IOMINP(IOM,ISYM,3),
     *                                EOMINP(IOM,ISYM,3)
  134                CONTINUE
  133             CONTINUE
               END IF
C
            GO TO 100
C
C--------------------------------
C           .STVEC : Choose start vectors.  
C--------------------------------
C
12          CONTINUE
               STVEC = .TRUE.
               READ (LUCMD,*) (NSTAR(ISYM),ISYM=1,MSYM)
               DO 331 ISYM = 1, MSYM
                  READ (LUCMD,*) (ISTVEC(K,ISYM),K=1,NSTAR(ISYM))
  331          CONTINUE
            GO TO 100
C
C
C--------------------------------------------------
C           .STOLD : Start from old vectors on file.
C--------------------------------------------------
C
13          CONTINUE
               STOLD = .TRUE.
            GO TO 100
C-----------------------------------------------------------------------
C           .CCTREN : normalize right eigenvectors for triples methods
C                     such that ( RE S+D+T | RE S+D+T ) = 1, default is
C                     to normalize as ( RE S+D | RE S+D ) = 1
C                     (see routine CCEXNORM)
C-----------------------------------------------------------------------
C
14          CONTINUE
               CCSDTRENRM = .TRUE.
            GO TO 100
C
C------------------------------------------------------------------------
C           .THREXC Set threshold for calculation of excitation energies.
C------------------------------------------------------------------------
C
15          CONTINUE
               READ (LUCMD, *) THREXC
            GO TO 100
C
C---------------------------------------------------------------
C           .CCSPIC Pick istate with right CCS excitation energy
C---------------------------------------------------------------
C

16          CONTINUE
               CCSPIC = .TRUE.
               READ(LUCMD,*) OMPCCS
            GO TO 100
C
C---------------------------------------------------------------
C           .CC2PIC Pick istate with right CC2 excitation energy
C---------------------------------------------------------------
C
17          CONTINUE
               CC2PIC = .TRUE.
               READ(LUCMD,*) OMPCC2
            GO TO 100
C
C-----------------------------------------------------------------
C           .CCSDPIC Pick istate with right CCSD excitation energy
C-----------------------------------------------------------------
C
18          CONTINUE
               CCSDPI = .TRUE.
               READ(LUCMD,*) OMPCCSD
            GO TO 100
C
C------------------------------------------------------------
C           .MARGIN;  Give margin in the 'picking' of states.
C------------------------------------------------------------
C
19          CONTINUE
               MARGIN = .TRUE.
               READ(LUCMD,*) XMARGIN
            GO TO 100
C
C----------------------------------------------------------------
C           .SQROVL  Compute full overlap matrix for eigenvectors
C                    (test option, see subroutine CCEXNORM)
C----------------------------------------------------------------
C
20          CONTINUE
               SQROVLP = .TRUE.
            GO TO 100
C
C----------------------------------------------------------------
C           .ANALYS  unused
C----------------------------------------------------------------
C
21          CONTINUE
               EXCI_CONT = .TRUE.
            GO TO 100
C
C----------------------------------------------------------------
C           .XXXXXX  unused
C----------------------------------------------------------------
C
22          CONTINUE
            GO TO 100
C
C----------------------------------------------------------------
C           .XXXXXX  unused
C----------------------------------------------------------------
C
23          CONTINUE
            GO TO 100
C
C----------------------------------------------------------------
C           .XXXXXX  unused
C----------------------------------------------------------------
C
24          CONTINUE
            GO TO 100
C
C----------------------------------------------------------------
C           .XXXXXX  unused
C----------------------------------------------------------------
C
25          CONTINUE
            GO TO 100
C
          ELSE
           WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF

      END IF

200   CONTINUE

*---------------------------------------------------------------------*
* post processing: consistency check, symmetry set up, etc:
*---------------------------------------------------------------------*

C     ------------------------------------------------------------
C     number of start vectors must equal to the number of required
C     excitation energies... (why?! why does a mix not work?)
C     ------------------------------------------------------------
      IF (STVEC ) THEN
         LSTVEC = .TRUE.
         DO ISYM = 1, MSYM
            IF (NSTAR(ISYM).NE.(NCCEXCI(ISYM,1)+NCCEXCI(ISYM,3))) THEN
               LSTVEC=.FALSE.
            END IF
         END DO

         IF (.NOT. LSTVEC) THEN
            CALL QUIT('Inconsistent input in *CCEXCI : '//
     &           'NSTAR .ne. NCCEXCI ')
         END IF
      ENDIF
                         
C     ------------------------------------------------------------
C     omega for triples calculation must be specified for all
C     states (singlet or triplet at the moment)
C     ------------------------------------------------------------
      DO IMULT = 1, 3, 2
C
         NOME  = 0
         DO ISYM = 1, MSYM
            NOME = NOME + NOMINP(ISYM,IMULT)
            IF (NOMINP(ISYM,IMULT) .GT. NCCEXCI(ISYM,IMULT)) THEN
               WRITE(LUPRI,*) ' NOMINP .GT. NCCEXCI for symmetry ',ISYM
               WRITE(LUPRI,*) ' and multiplicity ',IMULT
               CALL QUIT(' NOMINP .GT. NCCEXCI')
            ENDIF
         END DO                              
C
         IF ((OMESC.OR.(CCRT.OR.CCR3.OR.CCR1A.OR.CCR1B))
     *            .AND.(NOME .EQ. 0)) THEN
            MXTOMN = .TRUE.
            DO ISYM = 1, MSYM
               NOMINP(ISYM,IMULT) = NCCEXCI(ISYM,IMULT)
               DO IOM = 1, NOMINP(ISYM,IMULT)
                  IOMINP(IOM,ISYM,IMULT) = NOMINP(ISYM,IMULT) + 1 - IOM
                  EOMINP(IOM,ISYM,IMULT) = 0.0
               END DO
            END DO
C
         ENDIF                                       
C
      ENDDO
C
C     ----------------------------
C     set up symmetry information:
C     ----------------------------
      NEXCI  = 0
      NTRIP  = 0
      DO ISYM = 1,MSYM
         ISYOFE(ISYM) = NEXCI
         ITROFE(ISYM) = ISYOFE(ISYM) + NCCEXCI(ISYM,1)
         NEXCI        = ITROFE(ISYM) + NCCEXCI(ISYM,3)
         NTRIP        = NTRIP        + NCCEXCI(ISYM,3)
         DO IEX = ISYOFE(ISYM)+1, NEXCI
            ISYEXC(IEX) = ISYM
         END DO
         DO IEX = ISYOFE(ISYM)+1, ITROFE(ISYM)
            IMULTE(IEX) = 1
         END DO
         DO IEX = ITROFE(ISYM)+1, NEXCI
            IMULTE(IEX) = 3
         END DO
      END DO
C
      IF (IPRINT.GT.15) THEN
         WRITE(LUPRI,*) 'IN CC_EXCINP: '
         WRITE(LUPRI,*) 'NEXCI: ',NEXCI
         WRITE(LUPRI,*) 'Singlet: ',(NCCEXCI(J,1),J=1,MSYM)
         WRITE(LUPRI,*) 'Triplet: ',(NCCEXCI(J,3),J=1,MSYM)
         WRITE(LUPRI,*) 'ISYOFE:',(ISYOFE(J), J=1,MSYM)
         WRITE(LUPRI,*) 'ITROFE:',(ISYOFE(J), J=1,MSYM)
         WRITE(LUPRI,*) 'ISYEXC:',(ISYEXC(J), J=1,NEXCI)
         WRITE(LUPRI,*) 'IMULTE:',(IMULTE(J), J=1,NEXCI)
         WRITE(LUPRI,*) 'EIGVAL:',(EIGVAL(J), J=1,NEXCI)
      ENDIF
C
C     ---------------------------------------------------------------
C     if we are going for triplett states set flag for intermediates:
C     ---------------------------------------------------------------
C
      IF (NTRIP.GT.0) TRIPIM = .TRUE.
C
C     ----------------------------------------------------------
C     initialize eigenvalues with (non-degenerate) dummy values:
C     ----------------------------------------------------------
      DO IEXCI = 1, NEXCI
         EIGVAL(IEXCI) = 1.0D6 + NEXCI
      END DO
C
C---------------------------------------------------------------------
C     Finally if we are to calculate anything at all, put CCEXCI true.
C---------------------------------------------------------------------
C
      CCEXCI  = ((NEXCI.GT.0).OR.JACTST.OR.JACEXP.OR.FDJAC.OR.FDEXCI)
      IF (CCEXCI) RSPIM = .TRUE.
      IF (NEXCI .EQ. 0) THEN
         OSCSTR = .FALSE.
         NINFO = NINFO + 1
         WRITE(LUPRI,'(/A)') '@ INFO: No excitation energy requested'//
     &   ' even though CCEXCI is set - right?'
      END IF
C                                                          
      RETURN
      END
C=====================================================================*
c/* deck cc_lrsinp */
       SUBROUTINE CC_LRSINP(WORD,MSYM)
C---------------------------------------------------------------------*
C
C    Purpose: Read input for CC excited state calculations.
C
C    if (WORD .eq '*CCLRSD ') read & process input and set defaults,
C    else set only defaults
C
C    Ove Christiansen 24-10 1996
C
C=====================================================================*
#include <implicit.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccsdsym.h>
#include <cclr.h>
#include <leinf.h>
#include <cclrinf.h>
#include <ccrspprp.h>
#include <ccexci.h>
#include <cclres.h>

* local parameters:
      CHARACTER SECNAM*(9)
      PARAMETER (SECNAM='CC_LRSINP')

      INTEGER NTABLE
      PARAMETER (NTABLE = 20)

* variables:
      LOGICAL SET
      SAVE SET

      CHARACTER WORD*(7)
      CHARACTER LABELA*(8),LABELB*(8),LABHELP*70
      CHARACTER TABLE(NTABLE)*(8)

      INTEGER IJUMP,IDIP(3),IANG(3),IQUA(6)
* data:
      DATA SET /.FALSE./
      DATA TABLE /'.DIPOLE','.ECDLEN','.DIPLEN','.NO2N+1','.OPERAT',
     *            '.SELEXC','.DIPVEL','.DIPMIX','.ECDVEL','.OECDLE',
     *            '.OECDVE','.OLD_LR','.BOTHLR','.NEW_LR','.ECD   ',
     *            '.OECD  ','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/

*--------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (SET) RETURN
      SET = .TRUE.

*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*
C
      OSCSTR = .FALSE.
      VELSTR = .FALSE.
      MIXSTR = .FALSE.
      ROTLEN = .FALSE.
      ROTVEL = .FALSE.
      RTNLEN = .FALSE.
      RTNVEL = .FALSE.
      LRS2N1 = .TRUE. 
      SELLRS = .FALSE.
      OLDLRS = .FALSE.
      BOTHLRS = .FALSE.
C
      NSELRS = 0
      NLRSOP = 0
C
C     Other initializations
C

      ICHANG = 0

*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCLRSD') THEN

100   CONTINUE

* get new input line:
        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
          CALL UPCASE(WORD)
        END DO

        IF (WORD(1:1) .EQ. '.') THEN

c         table look up:
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO

c         jump to the appropiate input section:
          IF (IJUMP .LE. NTABLE) THEN
            ICHANG = ICHANG + 1
            GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
     &           IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_LRRINP.')

C
C-------------------------------------------------
C           Calculate dipole oscillator strengths.
C-------------------------------------------------
C
1           CONTINUE
              IF (OSCSTR) GO TO 100
              OSCSTR = .TRUE.
              IF (NLRSOP+9 .GT. MXLRSO) THEN
                WRITE(LUPRI,'(2(/A,I5))')
     &          ' NO. OF OPERATOR DOUBLES SPECIFIED  : ',NLRSOP+9,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSO
                CALL QUIT('TOO MANY OPERATOR DOUBLES IN CCLRSD.')
              END IF
              IDIP(1) = INDPRPCC('XDIPLEN ')
              IDIP(2) = INDPRPCC('YDIPLEN ')
              IDIP(3) = INDPRPCC('ZDIPLEN ')
              DO IDXA=1,3
                 DO IDXB=1,3
                   IDX = NLRSOP + (IDXA-1)*3+IDXB
                   IALRSOP(IDX) = IDIP(IDXA)
                   IBLRSOP(IDX) = IDIP(IDXB)
                 END DO
              END DO
              NLRSOP = NLRSOP + 9
            GO TO 100
C
C--------------------------------------------------------------
C           .ECDLEN: calculate length gauge rotatory strengths.
C--------------------------------------------------------------
C
2           CONTINUE
              IF (ROTLEN) GO TO 100
              ROTLEN = .TRUE.
              IF (NLRSOP+3 .GT. MXLRSO) THEN
                WRITE(LUPRI,'(2(/A,I5))')
     &          ' NO. OF OPERATOR DOUBLES SPECIFIED  : ',NLRSOP+3,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSO
                CALL QUIT('TOO MANY OPERATOR DOUBLES IN CCLRSD.')
              END IF
              IDIP(1) = INDPRPCC('XDIPLEN ')
              IDIP(2) = INDPRPCC('YDIPLEN ')
              IDIP(3) = INDPRPCC('ZDIPLEN ')
              IANG(1) = INDPRPCC('XANGMOM ')
              IANG(2) = INDPRPCC('YANGMOM ')
              IANG(3) = INDPRPCC('ZANGMOM ')
              DO IDXAB=1,3
                IDX = NLRSOP + IDXAB
                IALRSOP(IDX) = IDIP(IDXAB)
                IBLRSOP(IDX) = IANG(IDXAB)
              END DO
              NLRSOP = NLRSOP + 3
            GO TO 100
C
C-------------------------------------
C           .DIPLEN: alias for .DIPOLE
C-------------------------------------
C
3           CONTINUE
              IF (OSCSTR) GO TO 100
              GO TO 1
c           GO TO 100
C
C--------------------------------------------------------------------------
C           Do NOT Use 2n+1 rule expression for transition matrix elements.
C--------------------------------------------------------------------------
C
4           CONTINUE
              LRS2N1 = .FALSE.
            GO TO 100
C
C---------------------------
C           Input OPERATors.
C---------------------------
C
5           CONTINUE
              READ (LUCMD,'(2A)') LABELA, LABELB
              DO WHILE (LABELA(1:1).NE.'.' .AND. LABELA(1:1).NE.'*')
                IF (LABELA(1:1).NE.'!') THEN
                  IF (NLRSOP.LT.MXLRSO) THEN
                    NLRSOP = NLRSOP + 1
                    IALRSOP(NLRSOP) = INDPRPCC(LABELA)
                    IBLRSOP(NLRSOP) = INDPRPCC(LABELB)
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)')
     &               ' NO. OF OPERATOR DOUBLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSO
                    CALL QUIT('TOO MANY OPERATOR DOUBLES IN CCLRSD.')
                  END IF
                END IF
                READ (LUCMD,'(2A)') LABELA, LABELB
              END DO
              BACKSPACE(LUCMD)
            GO TO 100
C
C-------------------------
C           Select states.
C-------------------------
C
6           CONTINUE
              SELLRS =.TRUE.
              READ (LUCMD,'(A70)') LABHELP
              DO WHILE (LABHELP(1:1).NE.'.' .AND. LABHELP(1:1).NE.'*')
                IF (LABHELP(1:1).NE.'!') THEN
                  READ(LABHELP,*) IXSYM,IXST
                  IF (NSELRS.LT.MXLRSST) THEN
                    NSELRS = NSELRS + 1
                    ISELRS(NSELRS,1) = IXSYM
                    ISELRS(NSELRS,2) = IXST
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)')
     &               ' NO. OF STATES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSST
                    CALL QUIT('TOO MANY STATES IN CCLRS.')
                  END IF
                END IF
                READ (LUCMD,'(A70)') LABHELP
              END DO
              BACKSPACE(LUCMD)
            GO TO 100
C
C---------------------------------------------------------
C           Calculate velocity gauge oscillator strengths.
C---------------------------------------------------------
C
7           CONTINUE
              IF (VELSTR) GO TO 100
              VELSTR = .TRUE.
              IF (NLRSOP+9 .GT. MXLRSO) THEN
                WRITE(LUPRI,'(2(/A,I5))')
     &          ' NO. OF OPERATOR DOUBLES SPECIFIED  : ',NLRSOP+9,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSO
                CALL QUIT('TOO MANY OPERATOR DOUBLES IN CCLRSD.')
              END IF
              IDIP(1) = INDPRPCC('XDIPVEL ')
              IDIP(2) = INDPRPCC('YDIPVEL ')
              IDIP(3) = INDPRPCC('ZDIPVEL ')
              DO IDXA=1,3
                 DO IDXB=1,3
                   IDX = NLRSOP + (IDXA-1)*3+IDXB
                   IALRSOP(IDX) = IDIP(IDXA)
                   IBLRSOP(IDX) = IDIP(IDXB)
                 END DO
              END DO
              NLRSOP = NLRSOP + 9
            GO TO 100
C
C---------------------------------------------------------------
C           .DIPMIX: calculate mixed gauge oscillator strengths.
C---------------------------------------------------------------
C
8           CONTINUE
              IF (MIXSTR) GO TO 100
              MIXSTR = .TRUE.
              IF (NLRSOP+9 .GT. MXLRSO) THEN
                WRITE(LUPRI,'(2(/A,I5))')
     &          ' NO. OF OPERATOR DOUBLES SPECIFIED  : ',NLRSOP+9,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSO
                CALL QUIT('TOO MANY OPERATOR DOUBLES IN CCLRSD.')
              END IF
              IDIP(1) = INDPRPCC('XDIPLEN ')
              IDIP(2) = INDPRPCC('YDIPLEN ')
              IDIP(3) = INDPRPCC('ZDIPLEN ')
              IANG(1) = INDPRPCC('XDIPVEL ')
              IANG(2) = INDPRPCC('YDIPVEL ')
              IANG(3) = INDPRPCC('ZDIPVEL ')
              DO IDXA=1,3
                 DO IDXB=1,3
                   IDX = NLRSOP + (IDXA-1)*3+IDXB
                   IALRSOP(IDX) = IDIP(IDXA)
                   IBLRSOP(IDX) = IANG(IDXB)
                 END DO
              END DO
              NLRSOP = NLRSOP + 9
            GO TO 100
C
C----------------------------------------------------------------
C           .ECDVEL: calculate velocity gauge rotatory strengths.
C----------------------------------------------------------------
C
9           CONTINUE
              IF (ROTVEL) GO TO 100
              ROTVEL = .TRUE.
              IF (NLRSOP+3 .GT. MXLRSO) THEN
                WRITE(LUPRI,'(2(/A,I5))')
     &          ' NO. OF OPERATOR DOUBLES SPECIFIED  : ',NLRSOP+3,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSO
                CALL QUIT('TOO MANY OPERATOR DOUBLES IN CCLRSD.')
              END IF
              IDIP(1) = INDPRPCC('XDIPVEL ')
              IDIP(2) = INDPRPCC('YDIPVEL ')
              IDIP(3) = INDPRPCC('ZDIPVEL ')
              IANG(1) = INDPRPCC('XANGMOM ')
              IANG(2) = INDPRPCC('YANGMOM ')
              IANG(3) = INDPRPCC('ZANGMOM ')
              DO IDXAB=1,3
                IDX = NLRSOP + IDXAB
                IALRSOP(IDX) = IDIP(IDXAB)
                IBLRSOP(IDX) = IANG(IDXAB)
              END DO
              NLRSOP = NLRSOP + 3
            GO TO 100
C
C---------------------------------------------------------------------
C           .OECDLE: calculate length gauge rotatory strength tensors.
C---------------------------------------------------------------------
C
10          CONTINUE
              IF (RTNLEN) GO TO 100
              RTNLEN = .TRUE.
              ROTLEN = .TRUE.
              IF (NLRSOP+27 .GT. MXLRSO) THEN
                WRITE(LUPRI,'(2(/A,I5))')
     &          ' NO. OF OPERATOR DOUBLES SPECIFIED  : ',NLRSOP+27,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSO
                CALL QUIT('TOO MANY OPERATOR DOUBLES IN CCLRSD.')
              END IF
              IDIP(1) = INDPRPCC('XDIPLEN ')
              IDIP(2) = INDPRPCC('YDIPLEN ')
              IDIP(3) = INDPRPCC('ZDIPLEN ')
              IQUA(1) = INDPRPCC('XXSECMOM')
              IQUA(2) = INDPRPCC('XYSECMOM')
              IQUA(3) = INDPRPCC('XZSECMOM')
              IQUA(4) = INDPRPCC('YYSECMOM')
              IQUA(5) = INDPRPCC('YZSECMOM')
              IQUA(6) = INDPRPCC('ZZSECMOM')
              IANG(1) = INDPRPCC('XANGMOM ')
              IANG(2) = INDPRPCC('YANGMOM ')
              IANG(3) = INDPRPCC('ZANGMOM ')
              DO IDXA=1,3
                 DO IDXB=1,6
                   IDX = NLRSOP + (IDXA-1)*6+IDXB
                   IALRSOP(IDX) = IDIP(IDXA)
                   IBLRSOP(IDX) = IQUA(IDXB)
                 END DO
              END DO
              NLRSOP = NLRSOP + 18
              DO IDXA=1,3
                 DO IDXB=1,3
                   IDX = NLRSOP + (IDXA-1)*3+IDXB
                   IALRSOP(IDX) = IDIP(IDXA)
                   IBLRSOP(IDX) = IANG(IDXB)
                 END DO
              END DO
              NLRSOP = NLRSOP + 9
            GO TO 100
C
C-----------------------------------------------------------------------
C           .OECDVE: calculate velocity gauge rotatory strength tensors.
C-----------------------------------------------------------------------
C
11          CONTINUE
              IF (RTNVEL) GO TO 100
              RTNVEL = .TRUE.
              ROTVEL = .TRUE.
              IF (NLRSOP+27 .GT. MXLRSO) THEN
                WRITE(LUPRI,'(2(/A,I5))')
     &          ' NO. OF OPERATOR DOUBLES SPECIFIED  : ',NLRSOP+27,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSO
                CALL QUIT('TOO MANY OPERATOR DOUBLES IN CCLRSD.')
              END IF
              IDIP(1) = INDPRPCC('XDIPVEL ')
              IDIP(2) = INDPRPCC('YDIPVEL ')
              IDIP(3) = INDPRPCC('ZDIPVEL ')
              IQUA(1) = INDPRPCC('XXROTSTR')
              IQUA(2) = INDPRPCC('XYROTSTR')
              IQUA(3) = INDPRPCC('XZROTSTR')
              IQUA(4) = INDPRPCC('YYROTSTR')
              IQUA(5) = INDPRPCC('YZROTSTR')
              IQUA(6) = INDPRPCC('ZZROTSTR')
              IANG(1) = INDPRPCC('XANGMOM ')
              IANG(2) = INDPRPCC('YANGMOM ')
              IANG(3) = INDPRPCC('ZANGMOM ')
              DO IDXA=1,3
                 DO IDXB=1,6
                   IDX = NLRSOP + (IDXA-1)*6+IDXB
                   IALRSOP(IDX) = IDIP(IDXA)
                   IBLRSOP(IDX) = IQUA(IDXB)
                 END DO
              END DO
              NLRSOP = NLRSOP + 18
              DO IDXA=1,3
                 DO IDXB=1,3
                   IDX = NLRSOP + (IDXA-1)*3+IDXB
                   IALRSOP(IDX) = IDIP(IDXA)
                   IBLRSOP(IDX) = IANG(IDXB)
                 END DO
              END DO
              NLRSOP = NLRSOP + 9
            GO TO 100
C
C-----------------------------------------------------
C           .OLD_LR: use "old" LR residue program.
C           - "new" code differs only in the number
C             of evaluations of the transition moments
C             and, in particular, eta and ksi vectors.
C-----------------------------------------------------
C
12          CONTINUE
               OLDLRS = .TRUE.
            GO TO 100
C
C--------------------------------------------------------------------
C           .BOTHLR: use both the OLDLR and new codes (debug option).
C--------------------------------------------------------------------
C
13          CONTINUE
               BOTHLRS = .TRUE.
            GO TO 100
C
C-----------------------------------------------------
C           .NEW_LR: use "new" LR residue program.
C           - "new" code differs only in the number
C             of evaluations of the transition moments
C             and, in particular, eta and ksi vectors.
C-----------------------------------------------------
C
14          CONTINUE
               OLDLRS = .FALSE.
            GO TO 100
C
C----------------------------------------------------------------
C           .ECD   : calculate length and velocity gauge rotatory
C                    strengths.
C----------------------------------------------------------------
C
15          CONTINUE
              IF (.NOT.ROTLEN) THEN
                 ROTLEN = .TRUE.
                 IF (NLRSOP+3 .GT. MXLRSO) THEN
                   WRITE(LUPRI,'(2(/A,I5))')
     &             ' NO. OF OPERATOR DOUBLES SPECIFIED  : ',NLRSOP+3,
     &             ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSO
                   CALL QUIT('TOO MANY OPERATOR DOUBLES IN CCLRSD.')
                 END IF
                 IDIP(1) = INDPRPCC('XDIPLEN ')
                 IDIP(2) = INDPRPCC('YDIPLEN ')
                 IDIP(3) = INDPRPCC('ZDIPLEN ')
                 IANG(1) = INDPRPCC('XANGMOM ')
                 IANG(2) = INDPRPCC('YANGMOM ')
                 IANG(3) = INDPRPCC('ZANGMOM ')
                 DO IDXAB=1,3
                   IDX = NLRSOP + IDXAB
                   IALRSOP(IDX) = IDIP(IDXAB)
                   IBLRSOP(IDX) = IANG(IDXAB)
                 END DO
                 NLRSOP = NLRSOP + 3
              END IF
              IF (.NOT.ROTVEL) THEN
                 ROTVEL = .TRUE.
                 IF (NLRSOP+3 .GT. MXLRSO) THEN
                   WRITE(LUPRI,'(2(/A,I5))')
     &             ' NO. OF OPERATOR DOUBLES SPECIFIED  : ',NLRSOP+3,
     &             ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSO
                   CALL QUIT('TOO MANY OPERATOR DOUBLES IN CCLRSD.')
                 END IF
                 IDIP(1) = INDPRPCC('XDIPVEL ')
                 IDIP(2) = INDPRPCC('YDIPVEL ')
                 IDIP(3) = INDPRPCC('ZDIPVEL ')
                 IANG(1) = INDPRPCC('XANGMOM ')
                 IANG(2) = INDPRPCC('YANGMOM ')
                 IANG(3) = INDPRPCC('ZANGMOM ')
                 DO IDXAB=1,3
                   IDX = NLRSOP + IDXAB
                   IALRSOP(IDX) = IDIP(IDXAB)
                   IBLRSOP(IDX) = IANG(IDXAB)
                 END DO
                 NLRSOP = NLRSOP + 3
              END IF
            GO TO 100
C
C----------------------------------------------------------------
C           .OECD  : calculate length and velocity gauge rotatory
C                    strength tensors.
C----------------------------------------------------------------
C
16          CONTINUE
               IF (.NOT.RTNLEN) THEN
                  RTNLEN = .TRUE.
                  ROTLEN = .TRUE.
                  IF (NLRSOP+27 .GT. MXLRSO) THEN
                    WRITE(LUPRI,'(2(/A,I5))')
     &              ' NO. OF OPERATOR DOUBLES SPECIFIED  : ',NLRSOP+27,
     &              ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSO
                    CALL QUIT('TOO MANY OPERATOR DOUBLES IN CCLRSD.')
                  END IF
                  IDIP(1) = INDPRPCC('XDIPLEN ')
                  IDIP(2) = INDPRPCC('YDIPLEN ')
                  IDIP(3) = INDPRPCC('ZDIPLEN ')
                  IQUA(1) = INDPRPCC('XXSECMOM')
                  IQUA(2) = INDPRPCC('XYSECMOM')
                  IQUA(3) = INDPRPCC('XZSECMOM')
                  IQUA(4) = INDPRPCC('YYSECMOM')
                  IQUA(5) = INDPRPCC('YZSECMOM')
                  IQUA(6) = INDPRPCC('ZZSECMOM')
                  IANG(1) = INDPRPCC('XANGMOM ')
                  IANG(2) = INDPRPCC('YANGMOM ')
                  IANG(3) = INDPRPCC('ZANGMOM ')
                  DO IDXA=1,3
                     DO IDXB=1,6
                       IDX = NLRSOP + (IDXA-1)*6+IDXB
                       IALRSOP(IDX) = IDIP(IDXA)
                       IBLRSOP(IDX) = IQUA(IDXB)
                     END DO
                  END DO
                  NLRSOP = NLRSOP + 18
                  DO IDXA=1,3
                     DO IDXB=1,3
                       IDX = NLRSOP + (IDXA-1)*3+IDXB
                       IALRSOP(IDX) = IDIP(IDXA)
                       IBLRSOP(IDX) = IANG(IDXB)
                     END DO
                  END DO
                  NLRSOP = NLRSOP + 9
               END IF
               IF (.NOT.RTNVEL) THEN
                  RTNVEL = .TRUE.
                  ROTVEL = .TRUE.
                  IF (NLRSOP+27 .GT. MXLRSO) THEN
                    WRITE(LUPRI,'(2(/A,I5))')
     &              ' NO. OF OPERATOR DOUBLES SPECIFIED  : ',NLRSOP+27,
     &              ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSO
                    CALL QUIT('TOO MANY OPERATOR DOUBLES IN CCLRSD.')
                  END IF
                  IDIP(1) = INDPRPCC('XDIPVEL ')
                  IDIP(2) = INDPRPCC('YDIPVEL ')
                  IDIP(3) = INDPRPCC('ZDIPVEL ')
                  IQUA(1) = INDPRPCC('XXROTSTR')
                  IQUA(2) = INDPRPCC('XYROTSTR')
                  IQUA(3) = INDPRPCC('XZROTSTR')
                  IQUA(4) = INDPRPCC('YYROTSTR')
                  IQUA(5) = INDPRPCC('YZROTSTR')
                  IQUA(6) = INDPRPCC('ZZROTSTR')
                  IANG(1) = INDPRPCC('XANGMOM ')
                  IANG(2) = INDPRPCC('YANGMOM ')
                  IANG(3) = INDPRPCC('ZANGMOM ')
                  DO IDXA=1,3
                     DO IDXB=1,6
                       IDX = NLRSOP + (IDXA-1)*6+IDXB
                       IALRSOP(IDX) = IDIP(IDXA)
                       IBLRSOP(IDX) = IQUA(IDXB)
                     END DO
                  END DO
                  NLRSOP = NLRSOP + 18
                  DO IDXA=1,3
                     DO IDXB=1,3
                       IDX = NLRSOP + (IDXA-1)*3+IDXB
                       IALRSOP(IDX) = IDIP(IDXA)
                       IBLRSOP(IDX) = IANG(IDXB)
                     END DO
                  END DO
                  NLRSOP = NLRSOP + 9
               END IF
            GO TO 100
C
C------------------------------------
C           .XXXXXX: unused keywords.
C------------------------------------
C
17          CONTINUE
18          CONTINUE
19          CONTINUE
20          CONTINUE
            GO TO 100

          ELSE
           WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF

      END IF

200   CONTINUE
*---------------------------------------------------------------------*
* check, if input consistent.
*---------------------------------------------------------------------*
C
        IF (SELLRS.AND.(NSELRS.EQ.0)) WRITE(LUPRI,'(/A)') '@ INFO: '//
     &     '(*CCLRSD input is strange - no states is requested.)'
        IF (NLRSOP .EQ.0) WRITE(LUPRI,'(/A)') '@ INFO: '//
     &     '(*CCLRSD input ignored, because no operators requested.)'
C
C---------------------------------------------------------------------
C     Finally if we are to calculate anything at all, put CCLRSD true.
C---------------------------------------------------------------------
C
      CCLRSD  = (NLRSOP.GT.0)
C
      RETURN
      END
c/* deck cc_opainp */
      SUBROUTINE CC_OPAINP(WORD)
C---------------------------------------------------------------------*
C
C    Purpose: Read input for absorption strenghts
C             WORD='*CCOPA '  ground  to ex. state one-photon transit.
C             WORD='*CCTPA '  ground  to ex. state two-photon transit.
C             WORD='*CCXOPA'  excited to ex. state one-photon transit.
C
C    Christof Haettig, Dec 2002 / Oct 2003
C
C=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <ccsections.h>
#include <ccrspprp.h>
#include <ccopainf.h>
#include <cctpainf.h>
#include <ccxopainf.h>

* local parameters:
      CHARACTER SECNAM*(9)
      PARAMETER (SECNAM='CC_OPAINP')

      INTEGER NTABLE
      PARAMETER (NTABLE = 15)

* variables:
      LOGICAL SETGSTOPA, SETGSTTPA, SETXSTOPA
      SAVE SETGSTOPA, SETGSTTPA, SETXSTOPA

      CHARACTER WORD*(7)
      CHARACTER LABEL*(8), LABHELP*(80), LABELA*(8), LABELB*(8)
      CHARACTER TABLE(NTABLE)*(8)

      LOGICAL GSTOPA, GSTTPA, XSTOPA
      INTEGER IXSYM, IXSTATE, IXSYM2, IXSTATE2, IJUMP
      INTEGER INDPRPCC

#if defined (SYS_CRAY)
      REAL SMFREQ
#else
      DOUBLE PRECISION SMFREQ
#endif

* data:
      DATA SETGSTOPA /.FALSE./
      DATA SETGSTTPA /.FALSE./
      DATA SETXSTOPA /.FALSE./
      DATA TABLE / '.SELEXC','.NO2N+1','.OPERAT','.DIPLEN','.DIPVEL',
     *             '.ANGMOM','.HALFFR','.PRINT ','.USE X2','.USE O2',
     *             '.SELSTA','.STATES','.TRANSI','.SECMOM','.ROTSTR'/

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCOPA ') THEN
         GSTOPA = .TRUE.
         GSTTPA = .FALSE.
         XSTOPA = .FALSE.
         IF (SETGSTOPA) RETURN
         SETGSTOPA = .TRUE.
      ELSE IF (WORD(1:7) .EQ. '*CCXOPA') THEN
         GSTOPA = .FALSE.
         GSTTPA = .FALSE.
         XSTOPA = .TRUE.
         IF (SETXSTOPA) RETURN
         SETXSTOPA = .TRUE.
      ELSE IF (WORD(1:7) .EQ. '*CCTPA ') THEN
         GSTTPA = .TRUE.
         GSTOPA = .FALSE.
         XSTOPA = .FALSE.
         IF (SETGSTTPA) RETURN
         SETGSTTPA = .TRUE.
      ELSE
         CALL QUIT('CC_OPAINP called for wrong section:'//WORD(1:7))
      END IF

*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*
      IF (GSTOPA) THEN
        SELLRS = .FALSE.
        LRS2N1 = .TRUE. 
        NSELRS = 0
        NLRSOP = 0
      ELSE IF (XSTOPA) THEN
        SELQR2  = .FALSE.
        QR22N1  = .TRUE.
        NSEQR2  = 0 
        NQR2OP  = 0 
      ELSE IF (GSTTPA) THEN
        NSMSEL   = 0
        NSMOPER  = 0
        IPRSM    = 0
        HALFFR   = .FALSE.
        SELSMST  = .FALSE.
        LTPA_USE_X2 = .FALSE.
        LTPA_USE_O2 = .FALSE.
      END IF
 
*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*

100   CONTINUE

! get new input line:
        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
          CALL UPCASE(WORD)
        END DO

        IF (WORD(1:1) .EQ. '.') THEN

c         table look up:
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO

c         jump to the appropiate input section:
          IF (IJUMP .LE. NTABLE) THEN
            GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_OPAINP.')

C           ---------------------------------------------------------
C           .SELEXC, .SELSTA, .STATES, .TRANSI: 
C            select excited states / transitions
C           ---------------------------------------------------------
1           CONTINUE
11          CONTINUE
12          CONTINUE
13          CONTINUE

             IF (GSTOPA) THEN
              ! ground to excited state one-photon transition:
              !   READ IXSYM, IXSTATE    
              !   IXSYM   : symmetry class
              !   IXSTATE : state number within symmetry class
              SELLRS = .TRUE.
              READ (LUCMD,'(A80)') LABHELP
              DO WHILE (LABHELP(1:1).NE.'.' .AND. LABHELP(1:1).NE.'*')
                IF (LABHELP(1:1).NE.'!') THEN
                  READ(LABHELP,*) IXSYM,IXSTATE
                  IF (NSELRS.LT.MXLRSST) THEN
                   NSELRS = NSELRS + 1
                   ISELRSYM(NSELRS) = IXSYM
                   ISELRSTA(NSELRS) = IXSTATE
                  ELSE
                   NWARN = NWARN + 1
                   WRITE(LUPRI,'(/2A,I5//A,2I5/)')
     &              '@ WARNING: NO. OF STATES SPECIFIED',
     &              ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSST,
     &              '@ IGNORE STATE',IXSYM,IXSTATE 
                  END IF
                END IF
                READ (LUCMD,'(A80)') LABHELP
              END DO
             ELSE IF (XSTOPA) THEN
              ! excited to excited state one-photon transition:
              !   READ IXSYM, IXSTATE, IXSYM2, IXSTATE2
              !   IXSYM, IXSYM2     : symmetry classes
              !   IXSTATE, IXSTATE2 : state numbers within sym. classes
              SELQR2 = .TRUE.
              READ (LUCMD,'(A80)') LABHELP
              DO WHILE (LABHELP(1:1).NE.'.' .AND. LABHELP(1:1).NE.'*')
                IF (LABHELP(1:1).NE.'!') THEN
                  READ(LABHELP,*) IXSYM,IXSTATE,IXSYM2,IXSTATE2
                  IF (NSEQR2.LT.MXQR2ST) THEN
                   NSEQR2 = NSEQR2 + 1
                   ISEQR2SYM(NSEQR2,1) = IXSYM
                   ISEQR2STA(NSEQR2,1) = IXSTATE
                   ISEQR2SYM(NSEQR2,2) = IXSYM2
                   ISEQR2STA(NSEQR2,2) = IXSTATE2
                  ELSE
                   NWARN = NWARN + 1
                   WRITE(LUPRI,'(/2A,I5//A,2I5,I10,I5/)')
     &              '@ WARNING: NO. OF STATE PAIRS SPECIFIED',
     &              ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQR2ST,
     &              '@ IGNORE STATE PAIR',
     &              IXSYM,IXSTATE,IXSYM2,IXSTATE2
                  END IF
                END IF
                READ (LUCMD,'(A80)') LABHELP
              END DO
             ELSE IF (GSTTPA) THEN
              ! ground to excited state two-photon transition:
              !   READ IXSYM, IXSTATE, SMFREQ
              !   IXSYM   : symmetry class
              !   IXSTATE : state number within symmetry class
              !   SMFREQ  : photon energies associated with 2. operators
              SELSMST =.TRUE. 
              READ (LUCMD,'(A70)') LABHELP
              DO WHILE (LABHELP(1:1).NE.'.' .AND. LABHELP(1:1).NE.'*')
                IF (LABHELP(1:1).NE.'!') THEN
                  READ(LABHELP,*) IXSYM,IXSTATE,SMFREQ
                  IF (NSMSEL.LT.MXSMSEL) THEN
                    NSMSEL = NSMSEL + 1
                    ISMSEL(NSMSEL,1) = IXSYM
                    ISMSEL(NSMSEL,2) = IXSTATE
                    BSMFR(NSMSEL)    = SMFREQ
                  ELSE
                    WRITE(LUPRI,'(/A,I5)')
     &               ' NO. OF STATES SPECIFIED'//
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ' ,MXSMSEL 
                    CALL QUIT('TOO MANY STATES SPECIFIED BY .SELSTA')
                  END IF
                END IF
                READ (LUCMD,'(A70)') LABHELP
              END DO
             END IF

             BACKSPACE(LUCMD)
            GO TO 100

C           -----------------------------
C           .NO2N+1: do NOT use 2n+1 rule
C           -----------------------------
2           CONTINUE
              IF (GSTOPA) LRS2N1 = .FALSE.
              IF (XSTOPA) QR22N1 = .FALSE.
              IF (GSTTPA) CONTINUE
            GO TO 100
 
C           ------------------------
C           .OPERAT: operator labels 
C           ------------------------
3           CONTINUE
             IF ( GSTOPA .OR. XSTOPA ) THEN
              READ (LUCMD,'(A)') LABEL
              DO WHILE (LABEL(1:1).NE.'.' .AND. LABEL(1:1).NE.'*')
                IF (LABEL(1:1).NE.'!') THEN

                 IF (GSTOPA) THEN
                  IF (NLRSOP.LT.MXLRSO) THEN
                    NLRSOP = NLRSOP + 1
                    ILRSOP(NLRSOP) = INDPRPCC(LABEL)
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)')
     &               ' NO. OF OPERATORS SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSO
                    CALL QUIT('TOO MANY OPERATORS IN CC_OPAINP.')
                  END IF
                 ELSE IF (XSTOPA) THEN
                  IF (NQR2OP.LT.MXQR2O) THEN
                    NQR2OP = NQR2OP + 1
                    IQR2OP(NQR2OP) = INDPRPCC(LABEL)
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)')
     &               ' NO. OF OPERATORS SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQR2O
                    CALL QUIT('TOO MANY OPERATORS IN CC_OPAINP.')
                  END IF
                 END IF

                END IF
                READ (LUCMD,'(A)') LABEL
              END DO
             ELSE IF (GSTTPA) THEN
              READ (LUCMD,'(2A)') LABELA, LABELB
              DO WHILE (LABELA(1:1).NE.'.' .AND. LABELA(1:1).NE.'*')
                IF (LABELA(1:1).NE.'!') THEN
                  IF (NSMOPER.LT.MXSMOP) THEN
                    NSMOPER = NSMOPER + 1
                    IASMOP(NSMOPER) = INDPRPCC(LABELA)
                    IBSMOP(NSMOPER) = INDPRPCC(LABELB)
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)') 
     &               ' NO. OF OPERATOR PAIRS SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXSMOP
                    CALL QUIT('TOO MANY OPERATOR PAIRS IN CC_OPAINP.')
                  END IF
                END IF
                READ (LUCMD,'(2A)') LABELA, LABELB
              END DO
             ELSE
              CALL QUIT('Error in CC_OPAINP.')
             END IF
             BACKSPACE(LUCMD)
            GO TO 100
 
C           -----------------------------------------------------
C           .DIPLEN: calculate complete dipole transition vectors
C                    in length gauge 
C           -----------------------------------------------------
4           CONTINUE
            IF (GSTOPA) THEN
              CALL CC_PUT1OP(ILRSOP,NLRSOP,MXLRSO,'DIPLEN','CC_OPAINP')
            ELSE IF (XSTOPA) THEN
              CALL CC_PUT1OP(IQR2OP,NQR2OP,MXQR2O,'DIPLEN','CC_OPAINP')
            ELSE IF (GSTTPA) THEN
              CALL CC_PUT2OP(IASMOP,IBSMOP,NSMOPER,MXSMOP,
     &                       'DIPLEN','CC_OPAINP')
            END IF
            GO TO 100

C           -----------------------------------------------------
C           .DIPVEL: calculate complete dipole transition vectors
C                    in velocity gauge 
C           -----------------------------------------------------
5           CONTINUE
            IF (GSTOPA) THEN
              CALL CC_PUT1OP(ILRSOP,NLRSOP,MXLRSO,'DIPVEL','CC_OPAINP')
            ELSE IF (XSTOPA) THEN
              CALL CC_PUT1OP(IQR2OP,NQR2OP,MXQR2O,'DIPVEL','CC_OPAINP')
            ELSE IF (GSTTPA) THEN
              CALL CC_PUT2OP(IASMOP,IBSMOP,NSMOPER,MXSMOP,
     &                       'DIPVEL','CC_OPAINP')
            END IF
            GO TO 100

C           ------------------------------------------------------
C           .ANGMOM: calculate complete magnetic dipole transition 
C                    vectors and if possible rotatory strenghts
C           ------------------------------------------------------
6           CONTINUE
            IF (GSTOPA) THEN
              CALL CC_PUT1OP(ILRSOP,NLRSOP,MXLRSO,'ANGMOM','CC_OPAINP')
            ELSE IF (XSTOPA) THEN
              CALL CC_PUT1OP(IQR2OP,NQR2OP,MXQR2O,'ANGMOM','CC_OPAINP')
            ELSE IF (GSTTPA) THEN
              CALL CC_PUT2OP(IASMOP,IBSMOP,NSMOPER,MXSMOP,
     &                       'ANGMOM','CC_OPAINP')
            END IF
            GO TO 100

C           ------------------------------------------------
C           .HALFFR : impose condition of equal frequencies
C                      for the two lasers 
C           ------------------------------------------------
7           CONTINUE
              IF (GSTTPA) THEN
                HALFFR =.TRUE.
              ELSE
                WRITE(LUPRI,*) 'No .HALFFR keyword in section ',WORD
                WRITE(LUPRI,*) 'input will be ignored...'
              END IF
            GO TO 100

C           ------------
C           .PRINT 
C           ------------
8           CONTINUE
              IF (GSTTPA) THEN
                READ (LUCMD,*) IPRSM
              ELSE
                WRITE(LUPRI,*) 'No .PRINT keyword in section ',WORD
                WRITE(LUPRI,*) 'input will be ignored...'
              END IF
            GO TO 100

C           ------------
C           .USE X2
C           ------------
9           CONTINUE
              IF (GSTTPA) THEN
                LTPA_USE_X2 = .TRUE.
              ELSE
                WRITE(LUPRI,*) 'No .USE X2 keyword in section ',WORD
                WRITE(LUPRI,*) 'input will be ignored...'
              END IF
            GO TO 100

C           ------------
C           .USE O2
C           ------------
10          CONTINUE
              IF (GSTTPA) THEN
                LTPA_USE_O2 = .TRUE.
              ELSE
                WRITE(LUPRI,*) 'No .USE O2 keyword in section ',WORD
                WRITE(LUPRI,*) 'input will be ignored...'
              END IF
            GO TO 100

C           ------------------------------------------------------------
C           .SECMOM: calculate complete length gauge electric quadrupole
C                    transition vectors and if possible rotatory
C                    strength tensors.
C           ------------------------------------------------------------
14          CONTINUE
            IF (GSTOPA) THEN
              CALL CC_PUT1OP(ILRSOP,NLRSOP,MXLRSO,'SECMOM','CC_OPAINP')
            ELSE IF (XSTOPA) THEN
              CALL CC_PUT1OP(IQR2OP,NQR2OP,MXQR2O,'SECMOM','CC_OPAINP')
            ELSE IF (GSTTPA) THEN
              WRITE(LUPRI,*) 'No .SECMOM keyword in section ',WORD
              WRITE(LUPRI,*) 'input will be ignored...'
            END IF
            GO TO 100


C           ------------------------------------------------------
C           .ROTSTR: calculate complete velocity gauge electric
C                    quadrupole transition vectors and if possible
C                    rotatory strength tensors.
C           ------------------------------------------------------
15          CONTINUE
            IF (GSTOPA) THEN
              CALL CC_PUT1OP(ILRSOP,NLRSOP,MXLRSO,'ROTSTR','CC_OPAINP')
            ELSE IF (XSTOPA) THEN
              CALL CC_PUT1OP(IQR2OP,NQR2OP,MXQR2O,'ROTSTR','CC_OPAINP')
            ELSE IF (GSTTPA) THEN
              WRITE(LUPRI,*) 'No .ROTSTR keyword in section ',WORD
              WRITE(LUPRI,*) 'input will be ignored...'
            END IF
            GO TO 100

C
          ELSE
           WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF


200   CONTINUE
*---------------------------------------------------------------------*
* warning if for GSTTPA both .SELST AND .SELHLF are specified
*---------------------------------------------------------------------*
      IF (GSTTPA .AND. SELSMST .AND. HALFFR) THEN
         WRITE (LUPRI,*)
     &        ' WARNING: BOTH .SELST and .HALFFR are specified'
         WRITE (LUPRI,*) ' .HALFFR is used to obtain frequences'
      END IF

*----------------------------------------------------------------------*
* check, if any operator labels specified:
* if not, use default: dipole length and velocity, angular momentum, and
*                      electric dipole length and velocity.
*----------------------------------------------------------------------*
      IF (GSTOPA .AND. NLRSOP.EQ.0) THEN

        IF (NLRSOP+3 .LE. MXLRSO) 
     &    CALL CC_PUT1OP(ILRSOP,NLRSOP,MXLRSO,'DIPLEN','CC_OPAINP')
        IF (NLRSOP+3 .LE. MXLRSO) 
     &    CALL CC_PUT1OP(ILRSOP,NLRSOP,MXLRSO,'DIPVEL','CC_OPAINP')
        IF (NLRSOP+3 .LE. MXLRSO) 
     &    CALL CC_PUT1OP(ILRSOP,NLRSOP,MXLRSO,'ANGMOM','CC_OPAINP')
        IF (NLRSOP+6 .LE. MXLRSO) 
     &    CALL CC_PUT1OP(ILRSOP,NLRSOP,MXLRSO,'SECMOM','CC_OPAINP')
        IF (NLRSOP+6 .LE. MXLRSO) 
     &    CALL CC_PUT1OP(ILRSOP,NLRSOP,MXLRSO,'ROTSTR','CC_OPAINP')

      ELSE IF (XSTOPA .AND. NQR2OP.EQ.0) THEN

        IF (NQR2OP+3 .LE. MXQR2O) 
     &    CALL CC_PUT1OP(IQR2OP,NQR2OP,MXQR2O,'DIPLEN','CC_OPAINP')
        IF (NQR2OP+3 .LE. MXQR2O) 
     &    CALL CC_PUT1OP(IQR2OP,NQR2OP,MXQR2O,'DIPVEL','CC_OPAINP')
        IF (NQR2OP+3 .LE. MXQR2O) 
     &    CALL CC_PUT1OP(IQR2OP,NQR2OP,MXQR2O,'ANGMOM','CC_OPAINP')
        IF (NQR2OP+6 .LE. MXQR2O) 
     &    CALL CC_PUT1OP(IQR2OP,NQR2OP,MXQR2O,'SECMOM','CC_OPAINP')
        IF (NQR2OP+6 .LE. MXQR2O) 
     &    CALL CC_PUT1OP(IQR2OP,NQR2OP,MXQR2O,'ROTSTR','CC_OPAINP')


      ELSE IF (GSTTPA .AND. NSMOPER.EQ.0) THEN

        IF (NSMOPER+9 .LE. MXSMOP) 
     &    CALL CC_PUT2OP(IASMOP,IBSMOP,NSMOPER,MXSMOP,
     &                       'DIPLEN','CC_OPAINP')
        
      END IF

*---------------------------------------------------------------------*
* set CCOPA flag and return:
*---------------------------------------------------------------------*
      IF (GSTOPA) CCOPA  = .TRUE.
      IF (GSTTPA) CCTPA  = .TRUE.
      IF (XSTOPA) CCXOPA = .TRUE.

      RETURN
      END
C=====================================================================*
C                    END OF SUBROUTINE CC_OPAINP
C=====================================================================*
      SUBROUTINE CC_NODINP(WORD,INIT_ONLY)
C---------------------------------------------------------------------*
C
C    Purpose: read flags for different CC3 noddy code options
C
C    Christof Haettig, Jan 2003
C
C=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <ccnoddy.h>

* local parameters:
      CHARACTER SECNAM*(9)
      PARAMETER (SECNAM='CC_NODINP')

      INTEGER NTABLE
      PARAMETER (NTABLE = 20)

* variables:
      LOGICAL SET
      SAVE SET

      LOGICAL INIT_ONLY
      CHARACTER WORD*(7)
      CHARACTER LABEL*(8), LABHELP*(80)
      CHARACTER TABLE(NTABLE)*(8)
      INTEGER IJUMP

* data:
      DATA SET /.FALSE./
      DATA TABLE /'.XI    ','.XIDEN ','.ETA   ','.ETADEN','.FMAT  ',
     *            '.FNOALT','.OVLP  ','.OMEGA ','.LHTR  ','.RHTR  ',
     *            '.FOPDEN','.FINDIF','.FAMAT ','.GMAT  ','.BMAT  ',
     *            '.AAMAT ','.HMAT  ','.FADEN ','.XXXXXX','.XXXXXX'/

*---------------------------------------------------------------------*
* set defaults:
*---------------------------------------------------------------------*
      NODDY_INIT  = .FALSE.

      NODDY_OMEGA = .FALSE.
      NODDY_RHTR  = .FALSE.
      NODDY_LHTR  = .FALSE.
      NODDY_DEN   = .FALSE.
      NODDY_BMAT  = .FALSE.
      NODDY_FMAT  = .FALSE.
      NODDY_GMAT  = .FALSE.
      NODDY_HMAT  = .FALSE.

      NODDY_XI    = .FALSE.
      NODDY_ETA   = .FALSE. 
      NODDY_AAMAT = .FALSE.
      NODDY_FAMAT = .FALSE.

      NODDY_XI_ALTER  = .FALSE.
      NODDY_ETA_ALTER = .FALSE. 
      NODDY_FA_ALTER  = .FALSE.

      CCSDT_F_ALTER = .TRUE.

      NODDY_OVLP = .FALSE.
 
      IF (INIT_ONLY) RETURN

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .NE. '*NODDY') CALL 
     &  QUIT('CC_NODINP was call for wrong input section:'//WORD(1:7))

      IF (SET) RETURN
      SET = .TRUE.

      NODDY_INIT  = .TRUE. ! triggers precalculation of integrals etc.

*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
100   CONTINUE
! get new input line:
        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
          CALL UPCASE(WORD)
        END DO
        IF (WORD(1:1) .EQ. '.') THEN

c         table look up:
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO

c         jump to the appropiate input section:
          IF (IJUMP .LE. NTABLE) THEN
            GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), 
     &            IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_NODINP.')

C           ---------------------------------------
C           .XI    : use noddy code for Xksi vector
C           ---------------------------------------
1           CONTINUE
              NODDY_XI = .TRUE.
            GO TO 100

C           ------------------------------------
C           .XIDEN : use noddy Xksi density code
C           ------------------------------------
2           CONTINUE
              NODDY_XI_ALTER = .TRUE.
            GO TO 100
 
C           ---------------------------------------
C           .ETA   : use noddy code for Xksi vector
C           ---------------------------------------
3           CONTINUE
              NODDY_ETA = .TRUE.
            GO TO 100
 
C           ---------------------------------------
C           .ETADEN: use noddy for Eta density code
C           ---------------------------------------
4           CONTINUE
              NODDY_ETA_ALTER = .TRUE.
            GO TO 100

C           ------------------------------------------------------
C           .FMAT  : use noddy version for F matrix transformation
C           ------------------------------------------------------
5           CONTINUE
              NODDY_FMAT = .TRUE.
            GO TO 100

C           ------------------------------------------------------
C           .FNOALT: don't use alternative noddy code for F matrix
C                    which does triples as B matrix contraction
C           ------------------------------------------------------
6           CONTINUE
              NODDY_FMAT    = .TRUE.
              CCSDT_F_ALTER = .FALSE.
            GO TO 100

C           -------------------------------------------------------
C           .OVLP  : use noddy code for (LE|RE) and (RE|RE) overlap
C           -------------------------------------------------------
7           CONTINUE
              NODDY_OVLP = .TRUE.
            GO TO 100

C           -------------------------------------------
C           .OMEGA : use noddy code for vector function
C           -------------------------------------------
8           CONTINUE
              NODDY_OMEGA = .TRUE.
            GO TO 100

C           --------------------------------------------------------
C           .LHTR  : use noddy code for jacobian left transformation
C           --------------------------------------------------------
9           CONTINUE
              NODDY_LHTR  = .TRUE.
            GO TO 100

C           ---------------------------------------------------------
C           .RHTR  : use noddy code for jacobian right transformation
C           ---------------------------------------------------------
10          CONTINUE
              NODDY_RHTR  = .TRUE.
            GO TO 100
C
C           ------------------------------------------------
C           .FOPDEN: use noddy code for ground state density
C           ------------------------------------------------
11          CONTINUE
              NODDY_DEN  = .TRUE.
            GO TO 100
C
C           --------------------------------------------------------
C           .FINDIF: set flags appropriate for CC3 finite difference
C                    calculations
C           --------------------------------------------------------
12          CONTINUE
              NODDY_OMEGA  = .TRUE.
              NODDY_RHTR   = .TRUE.
              NODDY_LHTR   = .TRUE.
              NODDY_DEN    = .TRUE.
              NODDY_OVLP   = .TRUE.
              NODDY_ETA    = .TRUE.
              NODDY_XI     = .TRUE.
              NODDY_FMAT   = .TRUE.
              NODDY_FAMAT  = .TRUE.
              NODDY_GMAT   = .TRUE.
              NODDY_HMAT   = .TRUE.
              NODDY_BMAT   = .TRUE.
              NODDY_AAMAT  = .TRUE.

              NODDY_XI_ALTER  = .FALSE.
              NODDY_ETA_ALTER = .FALSE. 
              NODDY_FA_ALTER  = .FALSE.
            GO TO 100
C
C           ------------------------------------------------
C           .FAMAT: use noddy code for F{A} matrix 
C           ------------------------------------------------
13          CONTINUE
              NODDY_FAMAT = .TRUE.
            GO TO 100
C
C           ------------------------------------------------
C           .GMAT  : use noddy code for G matrix
C           ------------------------------------------------
14          CONTINUE
              NODDY_GMAT = .TRUE.
            GO TO 100
C
C           ------------------------------------------------
C           .BMAT  : use noddy code for B matrix
C           ------------------------------------------------
15          CONTINUE
              NODDY_BMAT = .TRUE.
            GO TO 100
C
C           ------------------------------------------------
C           .AAMAT : use noddy code for A{A} matrix
C           ------------------------------------------------
16          CONTINUE
              NODDY_AAMAT = .TRUE.
            GO TO 100
C
C           ------------------------------------------------
C           .HMAT  : use noddy code for H matrix
C           ------------------------------------------------
17          CONTINUE
              NODDY_HMAT = .TRUE.
            GO TO 100
C
C           ------------------------------------------------
C           .FADEN : use noddy code for F{A} densities
C           ------------------------------------------------
18          CONTINUE
              NODDY_FA_ALTER  = .TRUE.
            GO TO 100
C
C           ------------------------------------------------
C           .XXXXXX: unused
C           ------------------------------------------------
19          CONTINUE
            GO TO 100
C
C           ------------------------------------------------
C           .XXXXXX: unused
C           ------------------------------------------------
20          CONTINUE
            GO TO 100
C
          ELSE
           WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF

200   CONTINUE
*---------------------------------------------------------------------*
* Check flags for consistency and print some output:
*---------------------------------------------------------------------*
      IF (NODDY_XI_ALTER)  NODDY_XI    = .TRUE.
      IF (NODDY_ETA_ALTER) NODDY_ETA   = .TRUE.
      IF (NODDY_FA_ALTER)  NODDY_FAMAT = .TRUE.

      WRITE(LUPRI,*) 'The Triples section in the response will use '//
     &               'the following modules:'

      WRITE(LUPRI,'(/a)') '  in CCRHSN:'
      IF (NODDY_OMEGA) THEN
       WRITE(LUPRI,*) 'vector function calculation: CCSD_TRIPLE'
      ELSE
       WRITE(LUPRI,*) 'vector function calculation: CC3_OMEG'
      END IF

      WRITE(LUPRI,'(/a)') '  in CC_RHTR:'
      IF (NODDY_RHTR) THEN
       WRITE(LUPRI,*) '    A right transformation : CC_RHTR_NODDY'
      ELSE
       WRITE(LUPRI,*) '    A right transformation : CC3_OMEG'
      END IF

      WRITE(LUPRI,'(/a)') '  in CC_LHTR:'
      IF (NODDY_LHTR) THEN
       WRITE(LUPRI,*) '    A left transformation  : CC_LHTR_NODDY'
      ELSE
       WRITE(LUPRI,*) '    A left transformation  : CC3_T3/L3_LHTR'
      END IF

      WRITE(LUPRI,'(/a)') '  in CC_FOP:'
      IF (NODDY_DEN) THEN
       WRITE(LUPRI,*) '   for one-electron density: CCSDT_XI_CONT_NODDY'
      ELSE
       WRITE(LUPRI,*) '   for one-electron density: CCSDPT_DENS2'
      END IF

      WRITE(LUPRI,'(/a)') '  in CCEXNORM:'
      IF (NODDY_OVLP) THEN
       WRITE(LUPRI,*) '   for (LE|RE) overlapp    : CCOVLPT_NODDY '
      ELSE
       WRITE(LUPRI,*) '   for (LE|RE) overlapp    : CC3_LR_OVLP'
      END IF


      WRITE(LUPRI,'(/a)') '  in CC_XIETA:'

      IF (NODDY_XI) THEN
       WRITE(LUPRI,*) '    xi vector calculation: CCSDT_XI_NODDY'
       IF (NODDY_XI_ALTER) THEN 
         WRITE(LUPRI,*) '    xi contraction: CCSDT_XI_DEN_NODDY'
       ELSE
         WRITE(LUPRI,*) '    xi contraction: CCSDT_XI_NODDY'
       END IF
      ELSE
       WRITE(LUPRI,*) '    xi vector calculation: CC3_XI'
       WRITE(LUPRI,*) '    xi contraction: CC3_XI_DEN'
      END IF

      IF (NODDY_ETA) THEN
       WRITE(LUPRI,*)'    eta vector calculation: CCSDT_ETA_NODDY'
       WRITE(LUPRI,*)'    L A{O} transformation : CCSDT_ETA_NODDY'
       IF (NODDY_ETA_ALTER) THEN 
        WRITE(LUPRI,*)'    eta    contraction: CCSDT_ETA_DEN'
        WRITE(LUPRI,*)'    L A{O} contraction: CCSDT_A_DEN_NODDY'
       ELSE
        WRITE(LUPRI,*)'    eta    contraction: CCSDT_ETA_NODDY'
        WRITE(LUPRI,*)'    L A{O} contraction: CCSDT_ETA_NODDY'
       END IF
      ELSE
       WRITE(LUPRI,*) '    eta vector calculation: CC3_ETASD'
       WRITE(LUPRI,*) '    L A{O} transformation : CC3_ETASD'
       WRITE(LUPRI,*) '    eta    contraction: CCSDT_ETA_DEN'
       WRITE(LUPRI,*) '    L A{O} contraction: CCSDT_ETA_DEN'
      END IF

      WRITE(LUPRI,'(/a)') '  in CC_FMAT:'
      IF (NODDY_FMAT) THEN
       WRITE(LUPRI,*) '    F matrix transformation: CCSDT_FMAT_NODDY'
       WRITE(LUPRI,*) '    F matrix contraction   : CCSDT_FMAT_NODDY'//
     &                ' and CCSDT_FBC_NODDY'
      ELSE
       WRITE(LUPRI,*) '    F matrix transformation: CC3_FMAT'//
     &  ' and CC3_FT3B and CC3_FMATSD'
       WRITE(LUPRI,*) '    F matrix contraction   : CC3_FMAT'//
     &  ' and CCSDT_FBMAT'
      END IF

      WRITE(LUPRI,'(/a)') '  in CCQR_FADRV/CC_FAMAT:'
      IF (NODDY_FAMAT) THEN
       WRITE(LUPRI,*) '    F{A} matrix transform.  : CCSDT_FAMAT_NODDY'
       IF (NODDY_FA_ALTER) THEN
        WRITE(LUPRI,*) ' F{A} matrix contraction : CCSDT_FA_DEN/noddy'
       ELSE
        WRITE(LUPRI,*) '    F{A} matrix contraction : CCSDT_FAMAT_NODDY'
       END IF 
      ELSE
       WRITE(LUPRI,*) '    F{A} matrix transform.  : CCSDT_FAMAT_NODDY'
       WRITE(LUPRI,*) '    F{A} matrix contraction : CCSDT_FA_DEN'
      END IF

      IF (NODDY_GMAT) THEN
       WRITE(LUPRI,*)'    G matrix calculation: CCSDT_GMAT_NODDY'
      ELSE
       WRITE(LUPRI,*)'    G matrix calculation: CC3_GMAT'
      END IF

      IF (NODDY_BMAT) THEN
       WRITE(LUPRI,*)'    B matrix calculation: CCSDT_BMAT_NODDY'
      ELSE
       WRITE(LUPRI,*)'    B matrix calculation: CC3_BMAT'
      END IF

      IF (NODDY_AAMAT) THEN
       WRITE(LUPRI,*)'    A{A} matrix calculation: CCSDT_AAMAT_NODDY'
      ELSE
       WRITE(LUPRI,*)'    A{A} matrix calculation: CC3_AAMAT'
      END IF

      IF (NODDY_HMAT) THEN
       WRITE(LUPRI,*)'    H matrix calculation: CCSDT_HMAT_NODDY'
      ELSE
       WRITE(LUPRI,*)'    H matrix calculation: CC3_HMAT'
      END IF

      RETURN
      END
C=====================================================================*
C                    END OF SUBROUTINE CC_NODINP
C=====================================================================*
c/* deck cc_qr2rinp */
C=====================================================================*
       SUBROUTINE CC_QR2RINP(WORD)
C---------------------------------------------------------------------*
C
C    Purpose: Read input for CC excitec state calculations.
C
C    if (WORD .eq '*CCQR2R ') read & process input and set defaults,
C    else set only defaults
C
C    Ove Christiansen April 1997
C
C=====================================================================*
#include <implicit.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccsdsym.h>
#include <cclr.h>
#include <cclres.h>
#include <leinf.h>
#include <cclrinf.h>
#include <ccrspprp.h>
#include <ccexci.h>
#include <ccqr2r.h>

* local parameters:
      CHARACTER SECNAM*(10)
      PARAMETER (SECNAM='CC_QR2RINP')

      INTEGER NTABLE
      PARAMETER (NTABLE = 5)

* variables:
      LOGICAL SET
      SAVE SET

      CHARACTER WORD*(7)
      CHARACTER LABELA*(8),LABELB*(8),LABHELP*70
      CHARACTER TABLE(NTABLE)*(8)

      INTEGER IJUMP,IDIP(3)
* data:
      DATA SET /.FALSE./
      DATA TABLE /'.DIPOLE','.NO2N+1','.OPERAT','.SELEXC','.DIPVEL'/

*--------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (SET) RETURN
      SET = .TRUE.

*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*
C
      QR22N1 = .TRUE. 
      SELQR2 = .FALSE.
      XOSCST = .FALSE.
      XVELST = .FALSE.
C
      NSEQR2   = 0
      NQR2OP   = 0
C
C     Other initializations
C

      ICHANG = 0

*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCQR2R') THEN

100   CONTINUE

* get new input line:
        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
          CALL UPCASE(WORD)
        END DO

        IF (WORD(1:1) .EQ. '.') THEN

c         table look up:
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO

c         jump to the appropiate input section:
          IF (IJUMP .LE. NTABLE) THEN
            ICHANG = ICHANG + 1
            GOTO (1,2,3,4,5), IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_QR2RINP.')

C
C-------------------------------------------------
C           Calculate dipole oscillator strengths.
C-------------------------------------------------
C
1           CONTINUE
              IF (NQR2OP+9 .GT. MXQR2O) THEN
                WRITE(LUPRI,'(2(/A,I5))')
     &          ' NO. OF OPERATOR DOUBLES SPECIFIED  : ',NQR2OP+9,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQR2O
                CALL QUIT('TOO MANY OPERATOR DOUBLES IN CCQR2R.')
              END IF
              IDIP(1) = INDPRPCC('XDIPLEN ')
              IDIP(2) = INDPRPCC('YDIPLEN ')
              IDIP(3) = INDPRPCC('ZDIPLEN ')
              DO IDXA=1,3
                 DO IDXB=1,3
                   IDX = NQR2OP + (IDXA-1)*3+IDXB
                   IAQR2OP(IDX) = IDIP(IDXA)
                   IBQR2OP(IDX) = IDIP(IDXB)
                 END DO
              END DO
              NQR2OP = NQR2OP + 9
              XOSCST = .TRUE. 
            GO TO 100
C
C-------------------------------------------------------------------
C           Use 2n+1 rule expression for transition matrix elements.
C-------------------------------------------------------------------
C
2           CONTINUE
              QR22N1 = .FALSE. 
            GO TO 100
C
C---------------------------
C           Input OPERATors.
C---------------------------
C
3           CONTINUE
              READ (LUCMD,'(2A)') LABELA, LABELB
              DO WHILE (LABELA(1:1).NE.'.' .AND. LABELA(1:1).NE.'*')
                IF (LABELA(1:1).NE.'!') THEN
                  IF (NQR2OP.LT.MXQR2O) THEN
                    NQR2OP = NQR2OP + 1
                    IAQR2OP(NQR2OP) = INDPRPCC(LABELA)
                    IBQR2OP(NQR2OP) = INDPRPCC(LABELB)
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)')
     &               ' NO. OF OPERATOR DOUBLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQR2O
                    CALL QUIT('TOO MANY OPERATOR DOUBLES IN CCQR2R.')
                  END IF
                END IF
                READ (LUCMD,'(2A)') LABELA, LABELB
              END DO
              BACKSPACE(LUCMD)
            GO TO 100
C
C-------------------------
C           Select states.
C-------------------------
C
4           CONTINUE
              SELQR2 =.TRUE.
              READ (LUCMD,'(A70)') LABHELP
              DO WHILE (LABHELP(1:1).NE.'.' .AND. LABHELP(1:1).NE.'*')
                IF (LABHELP(1:1).NE.'!') THEN
                  READ(LABHELP,*) IXSYM,IXST,IXSYM2,IXST2
                  IF (NSEQR2.LT.MXQR2ST) THEN
                    NSEQR2 = NSEQR2 + 1
                    ISEQR2(NSEQR2,1) = IXSYM
                    ISEQR2(NSEQR2,2) = IXST
                    ISEQR2(NSEQR2,3) = IXSYM2
                    ISEQR2(NSEQR2,4) = IXST2
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)')
     &               ' NO. OF STATES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQR2ST
                    CALL QUIT('TOO MANY STATES IN CCQR2R.')
                  END IF
                END IF
                READ (LUCMD,'(A70)') LABHELP
              END DO
              BACKSPACE(LUCMD)
            GO TO 100
C
C-------------------------------------------------
C           Calculate dipole oscillator strengths.
C-------------------------------------------------
C
5           CONTINUE
              IF (NQR2OP+9 .GT. MXQR2O) THEN
                WRITE(LUPRI,'(2(/A,I5))')
     &          ' NO. OF OPERATOR DOUBLES SPECIFIED  : ',NQR2OP+9,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQR2O
                CALL QUIT('TOO MANY OPERATOR DOUBLES IN CCQR2R.')
              END IF
              IDIP(1) = INDPRPCC('XDIPVEL ')
              IDIP(2) = INDPRPCC('YDIPVEL ')
              IDIP(3) = INDPRPCC('ZDIPVEL ')
              DO IDXA=1,3
                 DO IDXB=1,3
                   IDX = NQR2OP + (IDXA-1)*3+IDXB
                   IAQR2OP(IDX) = IDIP(IDXA)
                   IBQR2OP(IDX) = IDIP(IDXB)
                 END DO
              END DO
              NQR2OP = NQR2OP + 9
              XVELST = .TRUE. 
            GO TO 100

          ELSE
           WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF

      END IF

200   CONTINUE
*---------------------------------------------------------------------*
* check, if input consistent.
*---------------------------------------------------------------------*
C
        IF (SELQR2.AND.(NSEQR2 .EQ.0)) WRITE(LUPRI,'(/A)') '@ INFO: '//
     &     '(*CCQR2R input is strange - no states is requested.)'
        IF (NQR2OP .EQ.0) WRITE(LUPRI,'(/A)') '@ INFO: '//
     &     '(*CCQR2R input ignored, because no operators requested.)'
C
C---------------------------------------------------------------------
C     Finally if we are to calculate anything at all, put CCQR2R true.
C---------------------------------------------------------------------
C
      CCQR2R  = (NQR2OP.GT.0)
C
      RETURN
      END
c/* deck cc_grin */
C=====================================================================*
       SUBROUTINE CC_GRIN(WORD,MSYM)
C---------------------------------------------------------------------*
C
C    Purpose: Read input for CC gradients: ground or excited state
C             walk.
C
C    if (WORD .eq '*CCGR   ') read & process input and set defaults,
C    else set only defaults
C
C    Ove Christiansen august-1997
C
C=====================================================================*
#include <implicit.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccsdsym.h>
#include <cclr.h>
#include <leinf.h>
#include <cclrinf.h>
#include <ccrspprp.h>
#include <ccexci.h>
#include <ccgr.h>
#include <ccfdgeo.h>

* local parameters:
      CHARACTER SECNAM*(7)
      PARAMETER (SECNAM='CC_GRIN')

      INTEGER NTABLE
      PARAMETER (NTABLE = 3)

* variables:
      LOGICAL SET
      SAVE SET

      CHARACTER WORD*(7)
      CHARACTER TABLE(NTABLE)*(8)

      INTEGER IJUMP
* data:
      DATA SET /.FALSE./
      DATA TABLE /'.XSTSYM','.XSTNUM','.NUMGD '/
*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (SET) RETURN
      SET = .TRUE.

*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*
C
      IXSTSY = 0
      IXSTAT = 0
      NUMGD  = .FALSE.
C
C     Other initializations
C
      ICHANG = 0

*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCGR  ') THEN

100   CONTINUE

* get new input line:
        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
          CALL UPCASE(WORD)
        END DO
        IF (WORD(1:1) .EQ. '.') THEN
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO
          IF (IJUMP .LE. NTABLE) THEN
            ICHANG = ICHANG + 1
            GOTO (1,2,3), IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_LRRINP.')
C
C-----------------------------------------
C           Readin excited state symmetry.
C-----------------------------------------
C
1           CONTINUE
               READ (LUCMD,*) IXSTSY
            GO TO 100
C
C---------------------------------------
C           Readin excited state number.
C---------------------------------------
C
2           CONTINUE
               READ (LUCMD,*) IXSTAT 
            GO TO 100
C
C-----------------------------------------------------------------------
C           Numerical differentiation and no analytical derivative calc.
C-----------------------------------------------------------------------
C
3           CONTINUE
               NUMGD = .TRUE.
            GO TO 100
C
          ELSE
           WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF

      END IF

200   CONTINUE
*---------------------------------------------------------------------*
* check, if input consistent.
*---------------------------------------------------------------------*
C
C---------------------------------------------------------------------
C     Finally if we are to calculate anything at all, put CCGR true.
C     Presently this means if numgd then calculate.
C---------------------------------------------------------------------
C
      CCGR   = NUMGD
C
      RETURN
      END
c/* deck cc_exgrin */
C=====================================================================*
       SUBROUTINE CC_EXGRIN(WORD,MSYM)
C---------------------------------------------------------------------*
C
C    Purpose: Read input for CC excited state calculations of 
C             first-order properties.
C
C    if (WORD .eq '*CCEXGR ') read & process input and set defaults,
C    else set only defaults
C
C    Ove Christiansen 4-2-1997
C
C=====================================================================*
#include <implicit.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccsdsym.h>
#include <cclr.h>
#include <leinf.h>
#include <cclrinf.h>
#include <ccrspprp.h>
#include <ccexci.h>
#include <ccexgr.h>

* local parameters:
      CHARACTER SECNAM*(9)
      PARAMETER (SECNAM='CC_EXGRIN')

      INTEGER NTABLE
      PARAMETER (NTABLE = 11)

* variables:
      LOGICAL SET
      SAVE SET

      CHARACTER WORD*(7)
      CHARACTER LABEL*(8), LABHELP*(70)
      CHARACTER TABLE(NTABLE)*(8)

      INTEGER IJUMP, INDPRPCC
* data:
      DATA SET /.FALSE./
      DATA TABLE /'.DIPOLE','.QUADRU','.NQCC  ','.OPERAT','XXXXXXX',
     *            'XXXXXXX','.ALLONE','.RELCOR','.SECMOM','.SELXST',
     *            '.SELEXC'/

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (SET) RETURN
      SET = .TRUE.

*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*
C
      SELXGR = .FALSE.
      SELXST = .FALSE.
      ALLEXE = .FALSE.
      XDIPMO = .FALSE.
      XQUADR = .FALSE.
      XNQCC  = .FALSE.
      XRELCO = .FALSE.
      XSECMO = .FALSE.
      NAXGRO = 0
      CCEXGR = .FALSE.
C     DNSDRV = .FALSE.
C
C     Other initializations
C
      ICHANG = 0

*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCEXGR') THEN

100   CONTINUE

* get new input line:
        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
          CALL UPCASE(WORD)
        END DO
        IF (WORD(1:1) .EQ. '.') THEN
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO
          IF (IJUMP .LE. NTABLE) THEN
            ICHANG = ICHANG + 1
            GOTO (1,2,3,4,5,6,7,8,9,10,11), IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_LRRINP.')
C
C-----------------------------------
C           Calculate dipole moment.
C-----------------------------------
C

1           CONTINUE
               XDIPMO = .TRUE.
            GO TO 100

C
C----------------------------------------
C           Calculate Quadrupole moments.
C----------------------------------------
C

2           CONTINUE
               XQUADR = .TRUE.
            GO TO 100

C
C----------------------------------------------
C           Calculate electric field gradients.
C----------------------------------------------
C

3           CONTINUE
               XNQCC   = .TRUE.
            GO TO 100
C
C----------------------------------------------
C           .OPERAT : General operator section.
C----------------------------------------------
C

4           CONTINUE
              READ (LUCMD,'(A)') LABEL
              DO WHILE ((LABEL(1:1).NE.'.' ).AND.(LABEL(1:1).NE.'*'))
                IF (LABEL(1:1).NE.'!') THEN
                  IF (NAXGRO .LT.MXGROP) THEN
                    NAXGRO  = NAXGRO  + 1
                    IAXGRO(NAXGRO) = INDPRPCC(LABEL)
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)') 
     &               ' NO. OF OPERATORS SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXGROP
                    CALL QUIT('TOO MANY OPERATORS IN CCEXGR.')
                  END IF
                END IF
                READ (LUCMD,'(3A)') LABEL
              END DO
              BACKSPACE(LUCMD)
            GO TO 100
C
C-----------------------------------------
C-----------------------------------------
C
5           CONTINUE
            GO TO 100
C
C---------------------------------------
C---------------------------------------
C
6           CONTINUE
            GO TO 100

C
C---------------------------------------------------------
C           Calculate all standard first order properties.
C---------------------------------------------------------
C
7           CONTINUE
               XDIPMO = .TRUE.
               XQUADR = .TRUE.
               XNQCC  = .TRUE.
               XRELCO = .TRUE.
               XSECMO = .TRUE.
            GO TO 100
C
C------------------------------------
C           Relativistic corrections.
C------------------------------------
C
8           CONTINUE
               XRELCO = .TRUE.
            GO TO 100
C
C--------------------------------
C           Second order moments.
C--------------------------------
C
9           CONTINUE
               XSECMO = .TRUE.
            GO TO 100
C
C---------------------------------------------------------------------
C           Select excited state for first order property calculation.
C---------------------------------------------------------------------
C
10          CONTINUE
               SELXST = .TRUE.
            GO TO 100
C
C---------------------------------------------------------------------
C           Select excited state for first order property calculation.
C---------------------------------------------------------------------
C
11          CONTINUE
              SELXGR = .TRUE.
              READ (LUCMD,'(A70)') LABHELP
              DO WHILE (LABHELP(1:1).NE.'.' .AND. LABHELP(1:1).NE.'*')
                IF (LABHELP(1:1).NE.'!') THEN
                  READ(LABHELP,*) IXSYM,IXST
                  IF (NSEXGR.LT.MXXGST) THEN
                    NSEXGR = NSEXGR + 1
                    ISEXGR(NSEXGR,1) = IXSYM
                    ISEXGR(NSEXGR,2) = IXST
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)')
     &               ' NO. OF STATES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXXGST
                    CALL QUIT('TOO MANY STATES IN CCEXGR.')
                  END IF
                END IF
                READ (LUCMD,'(A70)') LABHELP
              END DO
              BACKSPACE(LUCMD)
            GO TO 100
C
          ELSE
           WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF

      END IF

200   CONTINUE
*---------------------------------------------------------------------*
* check, if input consistent.
*---------------------------------------------------------------------*
C
C---------------------------------------------------------------------
C     Finally if we are to calculate anything at all, put CCEXGR true.
C---------------------------------------------------------------------
C
      CCEXGR = (XDIPMO.OR.XQUADR.OR.XNQCC.OR.XSECMO
     *          .OR.XRELCO.OR.(NAXGRO.GT.0))
C
      RETURN
      END
C---------------------------------------------------------------------*
c /* deck cc_fopinp */
C=====================================================================*
       SUBROUTINE CC_FOPINP(WORD)
C---------------------------------------------------------------------*
C
C  Purpose: read input for CC first order properties; 
C           directs calculation of dipole moments, quadrupole moments,
C           electric field gradients, etc.
C
C  if (WORD .eq '*CCFOP  ') read & process input and set defaults, 
C  else set only defaults 
C 
C  Asger Halkier & Ove Christiansen Oct. 1996/Mar. 1997(RELCOR&APROP)
C  Asger Halkier primo Nov. 1999: relativistic 2-electron Darwin term.
C  Asger Halkier ultimo Nov. 1999: First-order Direct Perturbation
C                                  Theory (DPT) energy corrections.
C
C=====================================================================*
#include <implicit.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccsdsym.h>
#include <cclr.h>
#include <ccfop.h>
#include <cclrinf.h>
#include <ccrspprp.h>

* local parameters:
      CHARACTER SECNAM*(9)
      PARAMETER (SECNAM='CC_FOPINP')

      INTEGER NTABLE
      PARAMETER (NTABLE = 14)

* variables:
      LOGICAL SET
      SAVE SET

      CHARACTER WORD*(7)
      CHARACTER LABEL*(8)
      CHARACTER TABLE(NTABLE)*(8)

      INTEGER IJUMP

* external function:
      INTEGER INDPRPCC

* data:
      DATA SET /.FALSE./
      DATA TABLE /'.DIPMOM','.QUADRU','.NQCC  ','.TSTDEN','.ALLONE',
     *            '.NONREL','.RELCOR','.OPERAT','.SECMOM','.2ELDAR',
     *            '.DPTECO','.BPH2OO','.BPH2SS','.CRONLY'/

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (SET) RETURN
      SET = .TRUE.

*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*
C
      DIPMOM = .FALSE.
      QUADRU = .FALSE.
      NQCC   = .FALSE.
      TSTDEN = .FALSE.
      SECMOM = .FALSE.
      RELCOR = .FALSE.
      RELORB = .TRUE.
      NAFOP  = 0
      DAR2EL = .FALSE.
      DPTECO = .FALSE.
      BP2EOO = .FALSE.
      CORRONLY = .FALSE.
C
*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCFOP ') THEN

100   CONTINUE
       
* get new input line:
        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
          CALL UPCASE(WORD)
        END DO

        IF (WORD(1:1) .EQ. '.') THEN

c         table look up:
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO
     
c         jump to the appropiate input section:
          IF (IJUMP .LE. NTABLE) THEN
            ICHANG = ICHANG + 1
            GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14), IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_FOPINP.')
      
C
C-----------------------------------
C           Calculate dipole moment.
C-----------------------------------
C

1           CONTINUE
               DIPMOM = .TRUE.
CCN         Added for CC-R12:
               IDUM = INDPRPCC('XDIPLEN ')
               IDUM = INDPRPCC('YDIPLEN ')
               IDUM = INDPRPCC('ZDIPLEN ')
            GO TO 100

C
C----------------------------------------
C           Calculate Quadrupole moments.
C----------------------------------------
C

2           CONTINUE
               QUADRU = .TRUE.
CCN         Added for CC-R12:
               IDUM = INDPRPCC('XXTHETA ')
               IDUM = INDPRPCC('XYTHETA ')
               IDUM = INDPRPCC('XZTHETA ')
               IDUM = INDPRPCC('YYTHETA ')
               IDUM = INDPRPCC('YZTHETA ')
               IDUM = INDPRPCC('ZZTHETA ')               
            GO TO 100

C
C----------------------------------------------
C           Calculate electric field gradients.
C----------------------------------------------
C

3           CONTINUE
               NQCC   = .TRUE.
            GO TO 100

C
C--------------------------
C           Test densities.
C--------------------------
C

4           CONTINUE
              TSTDEN = .TRUE.
            GO TO 100

C
C----------------------------------------------------------------------
C           Calculate all standard first-order one-electron properties.
C----------------------------------------------------------------------
C

5           CONTINUE
               DIPMOM = .TRUE.
               QUADRU = .TRUE.
               NQCC   = .TRUE.
               RELCOR = .TRUE.
               SECMOM = .TRUE.
            GO TO 100
C
C---------------------------------
C           No orbital relaxation.
C---------------------------------
C

6           CONTINUE
               RELORB = .FALSE.
            GO TO 100

C
C-------------------------------------------------
C           Relativistic one-electron corrections.
C-------------------------------------------------
C

7          CONTINUE
               RELCOR = .TRUE.
            GO TO 100
C
C-----------------------------------------------------
C           Arbitrary Number of One electron operator.
C-----------------------------------------------------
C

8           CONTINUE
              READ (LUCMD,'(A)') LABEL
              DO WHILE ((LABEL(1:1).NE.'.' ).AND.(LABEL(1:1).NE.'*'))
                IF (LABEL(1:1).NE.'!') THEN
                  IF (NAFOP .LT.MAFOP) THEN
                    NAFOP  = NAFOP + 1
                    IAFOP(NAFOP) = INDPRPCC(LABEL)
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)')
     &               ' NO. OF OPERATORS SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MAFOP 
                    CALL QUIT('TOO MANY OPERATORS IN CCFOP .')
                  END IF
                END IF
                READ (LUCMD,'(3A)') LABEL
              END DO
              BACKSPACE (LUCMD)
            GO TO 100  
C
C--------------------------------
C           Second order moments.
C--------------------------------
C

9          CONTINUE
               SECMOM = .TRUE.
CCN         Added for CC-R12:
               IDUM = INDPRPCC('XXSECMOM')
               IDUM = INDPRPCC('XYSECMOM')
               IDUM = INDPRPCC('XZSECMOM')
               IDUM = INDPRPCC('YYSECMOM')
               IDUM = INDPRPCC('YZSECMOM')
               IDUM = INDPRPCC('ZZSECMOM')
            GO TO 100
C
C-------------------------------------------------
C           Relativistic two-electron Darwin term.
C-------------------------------------------------
C

10         CONTINUE
               DAR2EL = .TRUE.
            GO TO 100
C
C-------------------------------------------------
C           Relativistic DPT
C-------------------------------------------------
C

11         CONTINUE
               DPTECO = .TRUE.
            GO TO 100
C
C
C-------------------------------------------------
C           Breit-Pauli Orbit-Orbit
C-------------------------------------------------
C

12         CONTINUE
               BP2EOO = .TRUE.
            GO TO 100
C
C
C
C-------------------------------------------------
C           Breit-Pauli Spin-Spin = -2 Darwin2E
C-------------------------------------------------
C

13         CONTINUE
               DAR2EL = .TRUE.
           GO TO 100
C
14         CONTINUE
               !removes Hartree-Fock part of densities 
               !yields correlation only contribution to
               !FOP properties
               CORRONLY  = .TRUE.
            GO TO 100
          ELSE
           WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A,/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF

      END IF

200   CONTINUE
C
C--------------------------------------------------------------------
C     Finally if we are to calculate anything at all, put CCFOP true.
C--------------------------------------------------------------------
C
      CCFOP   = (DIPMOM.OR.QUADRU.OR.NQCC.OR.TSTDEN.OR.RELCOR.OR.
     *           SECMOM.OR.DAR2EL.OR.DPTECO.OR.BP2EOO.OR.(NAFOP.GT.0))
C
      IF (CCFOP) RSPIM = .TRUE.
C
      RETURN
      END
C--------------------------------------------------------------------
c /* deck cc_lrinp */
C=====================================================================*
       SUBROUTINE CC_LRINP(WORD)
C---------------------------------------------------------------------*
C
C    Purpose: read input for CC linear response, in particular
C             dynamic polarizabilities
C
C    if (WORD .eq '*CCLR  ') read & process input and set defaults, 
C    else set only defaults 
C 
C    Christof Haettig and Ove Christiansen October 1996
C    Relaxed/Unrelaxed options introduced in Nov' 1998, Ch. Haettig
C
C=====================================================================*
C#if defined (IMPLICIT_NONE)
C      IMPLICIT NONE  
C#else
#  include <implicit.h>
C#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccsdsym.h>
#include <cclrinf.h>
#include <ccrspprp.h>
#include <mxcent.h>
#include <nuclei.h>
#include <codata.h>
Cholesky
#include <maxorb.h>
#include <ccdeco.h>
Cholesky
CTOCD  
#include <ctocdcc.h>
CTOCD  


* local parameters:
      CHARACTER SECNAM*(8)
      PARAMETER (SECNAM='CC_LRINP')

      INTEGER NTABLE
      PARAMETER (NTABLE = 25)

#if defined (SYS_CRAY)
      REAL ZERO, TOLFRQ
#else
      DOUBLE PRECISION ZERO, TOLFRQ
#endif
      PARAMETER (ZERO = 0.0d00)
      PARAMETER (TOLFRQ = 1.0D-09)

* variables:
      LOGICAL SET
      SAVE SET

      CHARACTER WORD*(7), LINE*(80)
      CHARACTER*8 LABELA,LABELB
      CHARACTER*8 LABDIP(3), LABDPV(3), LABANG(3)
      CHARACTER TABLE(NTABLE)*(8)

      LOGICAL LRELAX, LRELAS, LOCSTAT
      INTEGER IDX, IJUMP, IDIP(3), IGRA(MXCOOR)

!     LOGICAL EXCLRL
      INTEGER IGNCHO(4)

* external function:
      INTEGER INDPRPCC

* data:
      DATA LABDIP /'XDIPLEN ','YDIPLEN ','ZDIPLEN '/
      DATA LABDPV /'XDIPVEL ','YDIPVEL ','ZDIPVEL '/
      DATA LABANG /'XANGMOM ','YANGMOM ','ZANGMOM '/

      DATA SET /.FALSE./
      DATA TABLE /'.RELAXE','.UNRELA','.FREQUE','.DIPOLE','.ALLDSP',
     *            '.OLD_LR','.ASYMSD','.DISPCF','.OPERAT','.AVERAG',
     *            '.PRINT ','.STATIC','.DIPGRA','.OR LEN','.OR VEL', 
     *            '.OR    ','.OR MVE','.ORGANL','.ORIGIN','.WAVELE', 
     *            '.INCLRL','.EXCLRL','.CTOSUS','.CTOSHI','.XXXXXX'/ 

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (SET) RETURN
      SET = .TRUE.

*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*
      CAUCHY     = .FALSE.
      NLRDISP    = 0 
      ALLLRDSPCF = .FALSE.
C
      ALPHA_ISO   = .FALSE.
      ALPHA_ANI   = .FALSE.         
      OFFALPHA(1) = -1
      OFFALPHA(2) = -1
C 
      NBLRFR = 0
      NLROP  = 0
      NDIPFR = 0
      NORGIN = 0
      IPRSOP = IPRINT
      DIPPOL = .FALSE.
      ORLEN  = .FALSE.
      ORVEL  = .FALSE.
      ORMVE  = .FALSE.
      ORGANL = .FALSE.
C
Cmodvel
C
!     INCLRL = .FALSE.
!     EXCLRL = .FALSE.
C
Cmodvel
C
      ASYMSD = .FALSE.
      LRELAX = .FALSE.
      DIPGRA = .FALSE.
      ICHANG = 0

      OLDLR  = .FALSE.

      CALL IZERO(IGNCHO,4)
      CALL DZERO(ORGIN,MORGIN+1)

      LOCSTAT = .FALSE.

CTOCD
      CTOSHI = .FALSE.
      CTOSUS = .FALSE.
CTOCD

*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCLR  ') THEN

100   CONTINUE
       
* get new input line:
        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
          CALL UPCASE(WORD)
        END DO

        IF (WORD(1:1) .EQ. '.') THEN

c         table look up:
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO
     
c         jump to the appropriate input section:
          IF (IJUMP .LE. NTABLE) THEN
            ICHANG = ICHANG + 1
            GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &            21,22,23,24,25), IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_LRINP.')
C
C-----------------------
C           .RELAXEd
C-----------------------
C
1           CONTINUE
               LRELAX    = .TRUE.
               IF (CHOINT) THEN
                  IGNCHO(4) = 1
                  LRELAX = .FALSE.
               ENDIF
            GO TO 100
C
C-----------------------
C           .UNRELAxed
C-----------------------
C
2           CONTINUE
               LRELAX = .FALSE.
            GO TO 100
C
C---------------------
C           .FREQUEncy
C---------------------
C
3           CONTINUE
              READ (LUCMD,*) NRDFR
              NFTOT = NRDFR + NBLRFR
              IF (NFTOT .GT. MBLRFR) THEN
                WRITE(LUPRI,'(3(/A,I5),/)')
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',NFTOT,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MBLRFR,
     &          ' THE NUMBER IS RESET TO THE MAXIMUM : ',MBLRFR
                NFTOT = MBLRFR
                NRDFR = NFTOT - NBLRFR
              END IF 
              READ (LUCMD,*) (BLRFR(NBLRFR+I),I=1,NRDFR)
              NBLRFR = NBLRFR + NRDFR
            GO TO 100
C
C-----------------------------------
C           DIPole POLarizabilities.
C-----------------------------------
C
4           CONTINUE
              DIPPOL =.TRUE. 
              CALL CC_LRINPREQ(LABDIP,LABDIP,3,3,.FALSE.,LRELAX)
            GO TO 100
C
C           --------------------------------------------------------
C           .ALLDSP : do not skip odd/even dispersion coefficients
C                     or real/imaginary properties
C           --------------------------------------------------------
5           CONTINUE
              ALLLRDSPCF = .TRUE.
            GO TO 100

C
C           -----------------------------------
C           .OLD_LR : use old LR code^
C           -----------------------------------
C
6           CONTINUE
              OLDLR = .TRUE.
              IF (CHOINT) THEN
                 IGNCHO(3) = 1
                 OLDLR = .FALSE.
              ENDIF
            GO TO 100

C
C           -------------------------------------------------
C           Use asymmetric form for linear response function.
C           (Does not obey 2n+2 rule for multipliers but only
C            response to Y is needed.)
C           -------------------------------------------------
C
7           CONTINUE
              ASYMSD =.TRUE.
              ASYMSD =.TRUE.
              IF (CHOINT) THEN
                 IGNCHO(1) = 1
                 ASYMSD = .FALSE.
              ENDIF
            GO TO 100
C
C           ---------------------------------
C           .DISPCF : dispersion coefficients        
C           ---------------------------------
8           CONTINUE
              CAUCHY = .TRUE.
              READ (LUCMD,*) NLRDISP
              IF (NLRDISP.LT.0) THEN
                CALL QUIT('NLRDISP < 0 not allowed '//
     &                'for .DISPCF in *CCLR')
              END IF
              IF (CHOINT) THEN
                 IGNCHO(2) = 1
              ENDIF
            GO TO 100

C           -------------------------------------------
C           .OPERAT : Operator set for Linear response.
C           -------------------------------------------
 
9           CONTINUE
              READ (LUCMD,'(2A)') LABELA, LABELB
              DO WHILE (LABELA(1:1).NE.'.' .AND. LABELA(1:1).NE.'*')
                IF (LABELA(1:1).NE.'!') THEN
                  IF (NLROP.LT.MXLROP) THEN
                    CALL CC_LRINPREQ(LABELA,LABELB,1,1,.TRUE.,LRELAX)
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)')
     &               ' NO. OF OPERATOR DOUBLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLROP
                    CALL QUIT('TOO MANY OPERATOR DOUBLES IN CCLR.')
                  END IF
                END IF
                READ (LUCMD,'(2A)') LABELA, LABELB
              END DO 
              BACKSPACE(LUCMD)
            GO TO 100
C
C           -------------------------------------------------------
C           .AVERAG : calculate averaged tensor components
C                     implemented: alpha_{iso}, alpha_{ani}
C           -------------------------------------------------------
10          CONTINUE
              READ(LUCMD,'(A)') LINE
              IF (LINE(1:9).EQ.'ALPHA_ISO') THEN
                ALPHA_ISO = .TRUE.
              ELSE IF (LINE(:9).EQ.'ALPHA_ANI') THEN
                ALPHA_ISO = .TRUE.
                ALPHA_ANI = .TRUE.
              ELSE
                WRITE(LUPRI,'(/4A/A/)')
     &           '@ LABEL "',LINE(1:5),'" UNKNOWN FOR .AVERAG KEYWORD',
     &           'IN *CCLR SECTION.','@ INPUT IS IGNORED...'
              END IF

              READ(LUCMD,'(A)') LINE
              CSYM = 'GENERI'
              IF (LINE(1:6).EQ.'ATOMIC') THEN
                CSYM = 'ATOMIC'  ! an atom
              ELSE IF (LINE(1:6).EQ.'SPHTOP') THEN
                CSYM = 'SPHTOP'  ! spherical top
              ELSE IF (LINE(1:6).EQ.'LINEAR') THEN
                CSYM = 'LINEAR'  ! linear molecule
              ELSE IF (LINE(1:6).EQ.'XYDEGN') THEN
                CSYM = 'LINEAR'  ! linear molecule
              ELSE IF (LINE(1:5).EQ.'GENER') THEN
                CSYM = 'GENERI'  ! use generic point group symmetry
              ELSE
                WRITE (LUPRI,*)
     *                'WARNING: unknown symmetry input in *CCLR:'
                WRITE (LUPRI,*) LINE
                WRITE (LUPRI,*)'WARNING: input line ignored...'
              END IF
          
              IF (ALPHA_ISO .OR. ALPHA_ANI) THEN
                IDIP(1) = INDPRPCC('XDIPLEN ')
                IDIP(2) = INDPRPCC('YDIPLEN ')
                IDIP(3) = INDPRPCC('ZDIPLEN ')
                DO IDX = 1, 2
                  IALROP(NLROP+1)           = IDIP(3)   !cmp 1: alph_zz
                  IBLROP(NLROP+1)           = IDIP(3)

                  IALROP(NLROP+2+(IDX-1)*3) = IDIP(IDX) !cmp 2: alph_xx
                  IBLROP(NLROP+2+(IDX-1)*3) = IDIP(IDX) !cmp 5: alph_yy

                  IALROP(NLROP+3+(IDX-1)*3) = IDIP(IDX) !cmp 3: alph_xz
                  IBLROP(NLROP+3+(IDX-1)*3) = IDIP(3)   !cmp 6: alph_yz

                  IALROP(NLROP+4)           = IDIP(1)   !cmp 4: alph_xy
                  IBLROP(NLROP+4)           = IDIP(2)
                END DO
                DO IDX = 1, 6
                    LALORX(NLROP+IDX) = LRELAX
                    LBLORX(NLROP+IDX) = LRELAX
                END DO
                IF (     LRELAX) OFFALPHA(1) = NLROP
                IF (.NOT.LRELAX) OFFALPHA(2) = NLROP
                IF      (CSYM(1:6).EQ.'ATOMIC') THEN
                  NLROP = NLROP + 1
                ELSE IF (CSYM(1:6).EQ.'SPHTOP') THEN
                  NLROP = NLROP + 1
                ELSE IF (CSYM(1:6).EQ.'LINEAR') THEN
                  NLROP = NLROP + 3
                ELSE IF (CSYM(1:6).EQ.'XYDEGN') THEN
                  NLROP = NLROP + 4
                ELSE
                  NLROP = NLROP + 6
                END IF
              END IF
            GO TO 100                                            
C
C           --------------------------------------------------
C           .PRINT set print level for linear response output:
C           --------------------------------------------------
C
11           CONTINUE
               READ (LUCMD,*) IPRSOP
             GO TO 100
C
C-----------------------
C           .STATIC
C-----------------------
C
12          CONTINUE
              IF (.NOT. LOCSTAT) THEN
                 IF (NBLRFR .GE. MBLRFR) THEN
                   WRITE(LUPRI,'(3(/A,I5),/)')
     &            '@ NUMBER OF FREQUENCIES SPECIFIED         : ',NBLRFR,
     &            '@ IS ALREADY EQUAL TO THE MAXIMUM ALLOWED : ',MBLRFR,
     &            '@ THE .STATIC KEYWORD UNDER *CCLR WILL BE IGNORED...'
                 ELSE
                   LOCSTAT = .TRUE.
                   NBLRFR = NBLRFR + 1
                   BLRFR(NBLRFR) = 0.0D0
                 END IF
              ENDIF
            GO TO 100
C
C           -------------------------------------------------
C           .DIPGRA: Dipole gradients and Cioslowski charges.
C           -------------------------------------------------
C
13          CONTINUE
               DIPGRA = .TRUE.
               NDIP   = 3
               NCOOR  = 3*NUCDEP
               NTOT   = NDIP*NCOOR
               IF (NLROP+NTOT .GT. MXLROP) THEN
                  WRITE(LUPRI,'(2(/A,I5))')
     &            ' NO. OF OPERATOR DOUBLES SPECIFIED  : ',NLROP+NTOT,
     &            ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLROP
                  CALL QUIT('TOO MANY OPERATOR DOUBLES IN CCLR')
               ENDIF
               IF (NCOOR .GT. 999) THEN
                  WRITE(LUPRI,*) 'Too many centers in CCLR'
                  WRITE(LUPRI,*)
     &            'Unable to construct labels for gradient ints'
                  CALL QUIT('Too many centers in CCLR')
               ENDIF
               IF (NCOOR .GT. MXCOOR) THEN
                  WRITE(LUPRI,*) 'IGRA dimension error in CC_LRINP:'
                  WRITE(LUPRI,*) ' NCOOR: ',NCOOR
                  WRITE(LUPRI,*) 'MXCOOR: ',MXCOOR
                  CALL QUIT('Error in CC_LRINP')
               ENDIF
               DO I = 1,NCOOR
                  WRITE(LABELA,'(A5,I3)') '1DHAM',I
                  DO J = 6,8
                     IF (LABELA(J:J) .EQ. ' ') LABELA(J:J) = '0'
                  ENDDO
                  IGRA(I) = INDPRPCC(LABELA)
               ENDDO
               IDIP(1) = INDPRPCC('XDIPLEN ')
               IDIP(2) = INDPRPCC('YDIPLEN ')
               IDIP(3) = INDPRPCC('ZDIPLEN ')
               DO IDXB = 1,NDIP
                  DO IDXA = 1,NCOOR
                     IDX = NLROP + NCOOR*(IDXB - 1) + IDXA
                     IALROP(IDX) = IGRA(IDXA)
                     IBLROP(IDX) = IDIP(IDXB)
                     LALORX(IDX) = .TRUE.  ! Force orb. relax. for grad.
                     LBLORX(IDX) = .TRUE.  ! Force orb. relax. for dip.
                  ENDDO
               ENDDO
               NLROP = NLROP + NTOT
            GO TO 100
C
C------------------------------------------------------
C           '.OR LEN': Optical Rotation - LENgth gauge.
C------------------------------------------------------
C
14          CONTINUE
              IF (.NOT. ORLEN) THEN
                 ORLEN =.TRUE.
                 CALL CC_LRINPREQ(LABDIP,LABANG,3,3,.TRUE.,LRELAX)
              ENDIF
            GO TO 100
C
C--------------------------------------------------------
C           '.OR VEL': Optical Rotation - VELocity gauge.
C--------------------------------------------------------
C
15          CONTINUE
              IF (.NOT. ORVEL) THEN
                 ORVEL  = .TRUE.
                 CALL CC_LRINPREQ(LABDPV,LABANG,3,3,.TRUE.,LRELAX)
              ENDIF
 
C----------------------------------------------------
C           '.OR    ': same as '.OR MVE' + '.OR LEN'.
C----------------------------------------------------
C
16          CONTINUE
               IF (.NOT. ORMVE) THEN
                  ORMVE = .TRUE.
                  CALL CC_LRINPREQ(LABDPV,LABANG,3,3,.TRUE.,LRELAX)
                  IF (.NOT. LOCSTAT) THEN
                     IF (NBLRFR .GE. MBLRFR) THEN
                        WRITE(LUPRI,'(2(/A,I5))')
     &             ' NUMBER OF FREQUENCIES SPECIFIED         : ',NBLRFR,
     &             ' IS ALREADY EQUAL TO THE MAXIMUM ALLOWED : ',MBLRFR
                        CALL QUIT('Request for .STATIC under .OR MVE '
     &                            //'failed in CCLR.')
                     ELSE
                        LOCSTAT = .TRUE.
                        NBLRFR  = NBLRFR + 1
                        BLRFR(NBLRFR) = 0.0D0
                     END IF
                  END IF
               END IF
               IF (.NOT. ORLEN) THEN
                  ORLEN =.TRUE.
                  CALL CC_LRINPREQ(LABDIP,LABANG,3,3,.TRUE.,LRELAX)
               ENDIF
            GO TO 100
C
C-------------------------------------------------------------------
C           '.OR MVE': Opt. Rot., modified velocity gauge.
C                      I.e. correct for unphysical static component.
C-------------------------------------------------------------------
C
17          CONTINUE
               IF (.NOT. ORMVE) THEN
                  ORMVE = .TRUE.
                  CALL CC_LRINPREQ(LABDPV,LABANG,3,3,.TRUE.,LRELAX)
                  IF (.NOT. LOCSTAT) THEN
                     IF (NBLRFR .GE. MBLRFR) THEN
                        WRITE(LUPRI,'(2(/A,I5))')
     &             ' NUMBER OF FREQUENCIES SPECIFIED         : ',NBLRFR,
     &             ' IS ALREADY EQUAL TO THE MAXIMUM ALLOWED : ',MBLRFR
                        CALL QUIT('Request for .STATIC under .OR MVE '
     &                            //'failed in CCLR.')
                     ELSE
                        LOCSTAT = .TRUE.
                        NBLRFR  = NBLRFR + 1
                        BLRFR(NBLRFR) = 0.0D0
                     END IF
                  END IF
               END IF
            GO TO 100
C
C------------------------------------------------------------------------
C           '.ORGANL': Calculate OR LEN origin dependence (Delta-vector).
C------------------------------------------------------------------------
C
18          CONTINUE
              IF (.NOT. ORGANL) THEN
                 ORGANL = .TRUE.
                 CALL CC_LRINPREQ(LABDIP,LABDPV,3,3,.FALSE.,LRELAX)
              ENDIF
            GO TO 100
C
C------------------------------------------------------------------------
C           '.ORIGIN': Additional origins for evaluating OR length gauge.
C                      Implies '.ORGANL' and '.OR LEN'.
C------------------------------------------------------------------------
C
19          CONTINUE
              READ(LUCMD,*) NORGIN
              NORGSV = NORGIN
              IF (NORGIN .GT. MORGIN) THEN
                WRITE(LUPRI,'(3(/A,I5))')
     &          ' NUMBER OF OR ORIGINS  SPECIFIED    : ',NORGIN,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MORGIN,
     &          ' THE NUMBER IS RESET TO THE MAXIMUM : ',MORGIN
                NORGIN = MORGIN
              ENDIF
              DO J = 1,NORGIN
                 READ(LUCMD,*) (ORGIN(I,J), I=1,3)
              ENDDO
              DO J = NORGIN+1,NORGSV
                 READ(LUCMD,*) SCR1,SCR2,SCR3
              ENDDO
              IF (.NOT. ORGANL) THEN
                 ORGANL = .TRUE.
                 CALL CC_LRINPREQ(LABDIP,LABDPV,3,3,.FALSE.,LRELAX)
              ENDIF
              IF (.NOT. ORLEN) THEN
                 ORLEN = .TRUE.
                 CALL CC_LRINPREQ(LABDIP,LABANG,3,3,.TRUE.,LRELAX)
              ENDIF
            GO TO 100
C
C-----------------------------------------------------------------------
C           '.WAVELE': Wavelengths in nm (instead of frequencies in au).
C-----------------------------------------------------------------------
C
20          CONTINUE
              READ (LUCMD,*) NWAVEL
              NFTOT = NWAVEL + NBLRFR
              IF (NFTOT .GT. MBLRFR) THEN
                WRITE(LUPRI,'(3(/A,I5))')
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',NFTOT,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MBLRFR,
     &          ' THE NUMBER IS RESET TO THE MAXIMUM : ',MBLRFR
                NFTOT  = MBLRFR
                NWAVEL = NFTOT - NBLRFR
              END IF
              READ (LUCMD,*) (BLRFR(NBLRFR+I),I=1,NWAVEL)
              DO I = 1,NWAVEL
                 IF (DABS(BLRFR(NBLRFR+I)) .LE. 1.0D-7) THEN
                    WRITE(LUPRI,'(/A,I5,A,1P,D22.15,A)')
     &              'Wavelength number',I,' too small: ',
     &              BLRFR(NBLRFR+I),' nm'
                    WRITE(LUPRI,'(A/)')
     &              'Input frequency (in au) instead (.FREQUE keyword).'
                    CALL QUIT('Input wavelength too small in '//SECNAM)
                 ENDIF
                 XWAV = BLRFR(NBLRFR+I)
                 BLRFR(NBLRFR+I) = XTNM/XWAV
              ENDDO
              NBLRFR = NBLRFR + NWAVEL
            GO TO 100
C
C---------------------------------------------------------
C           '.INCLRL': Include commutator terms in OR VEL.
C---------------------------------------------------------
C
21          CONTINUE
              call quit('.INCLRL not implemented in this version')
!             INCLRL = .TRUE.
            GO TO 100
C
C---------------------------------------------------------
C           '.EXCLRL': Exclude commutator terms in OR VEL.
C---------------------------------------------------------
C
22          CONTINUE
              call quit('.EXCLRL not implemented in this version')
!             EXCLRL = .TRUE.
            GO TO 100
C
C-------------------------------------------
C           '.CTOSUS': CTOCD susceptibility.
C-------------------------------------------
C
23          CONTINUE
              CTOSUS = .TRUE.
            GO TO 100
C
C--------------------------------------
C           '.CTOSHI': CTOCD shielding.
C--------------------------------------
C
24          CONTINUE
              CTOSHI = .TRUE.
            GO TO 100
C
C           -------------------------------------------------
C           .XXXXXXX unused keywords
C           -------------------------------------------------
C
25          CONTINUE
            GO TO 100
C
          ELSE
           WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF

      END IF

200   CONTINUE
C
C------------------------------------------------
C     Do some checking for Cholesky calculations:
C------------------------------------------------
C
      IGNSUM = 0
      DO I = 1,4
         IGNSUM = IGNSUM + ABS(IGNCHO(I))
      ENDDO
      IF (CHOINT .AND. (IGNSUM.NE.0)) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(A,A)')
     &   SECNAM,
     &   ': WARNING: *CCLR options not implemented for Cholesky job:'
         IF (IGNCHO(1) .NE. 0) THEN
            WRITE(LUPRI,'(A)')
     &      '.ASYMSD ignored (default 2n+2 expression will be used).'
         ENDIF
         IF (IGNCHO(2) .NE. 0) THEN
            WRITE(LUPRI,'(A)')
     &      'FATAL ERROR: No Cauchy for Cholesky job! (.DISPCF option)'
         ENDIF
         IF (IGNCHO(3) .NE. 0) THEN
            WRITE(LUPRI,'(A)')
     &      '.OLD_LR ignored (default code will be used).'
         ENDIF
         IF (IGNCHO(4) .NE. 0) THEN
            WRITE(LUPRI,'(A)')
     &      '.RELAXE ignored (unrelaxed properties will be calculated).'
         ENDIF
         WRITE(LUPRI,*)
         IF (IGNCHO(2) .NE. 0) CALL QUIT('Error in '//SECNAM)
      ENDIF
C
C--------------------------------------
C     Include operator pairs for CTOCD.
C--------------------------------------
C
      CTOMAG = CTOSUS .OR. CTOSHI
      IF (CTOMAG) CALL CC_CTOMAG
C
C------------------------------------------------------------
C     Check if commutator terms are to be included in OR VEL:
C------------------------------------------------------------
C
!     IF (EXCLRL) INCLRL = .FALSE.
C
C-----------------------------------------------
C check, if operators and frequencies specified:
C-----------------------------------------------
C
      IF (NBLRFR.EQ.0 .AND. (.NOT.CAUCHY)) THEN
         NBLRFR   = 1
         BLRFR(1) = 0.0D0
      ENDIF
C
      IF (ICHANG .NE. 0) THEN 
        IF (NLROP .EQ.0) WRITE(LUPRI,'(/A)')
     &     '(*CCLR   input ignored, because no operators requested.)'
      END IF
C
C----------------------------
C     Make wa frequency list. 
C----------------------------
      DO IFREQ = 1, NBLRFR
        ALRFR(IFREQ) = - BLRFR(IFREQ) 
      END DO
C
C-------------------------------------------------------------------
C     Finally if we are to solve for anything at all, put CCLR true.
C-------------------------------------------------------------------
C
      CCLR   = (NLROP.GT.0)
      IF (CCLR) RSPIM = .TRUE.
C
      RETURN
      END
*---------------------------------------------------------------------*
C  /* Deck cc_lrinpreq */
      SUBROUTINE CC_LRINPREQ(LABELA,LABELB,NA,NB,LDIAGO,LRELAX)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Request linear response calculation of the tensor
C              <<LABELA(i),LABELB(j)>> for i = 1,NA and j=1,NB.
C
C              If LDIAGO: request diagonal only (NA=NB only!).
C
C              LRELAX is the flag that will be associated with
C              each perturbation operator for relaxation.
C
#include <implicit.h>
      CHARACTER*8 LABELA(NA), LABELB(NB)
      LOGICAL     LDIAGO, LRELAX
#include <cclrinf.h>
#include <priunit.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CC_LRINPREQ')

      INTEGER IOPA(MXLROP), IOPB(MXLROP)

      IF ((NA.GT.MXLROP) .OR. (NB.GT.MXLROP)) THEN
         WRITE(LUPRI,'(//A,A,A)')
     &   ' Too many operators in ',SECNAM,':'
         WRITE(LUPRI,'(A,I10/A,I10)')
     &   ' NA =',NA,' NB =',NB
         CALL QUIT('Too many operators in '//SECNAM)
      ELSE IF ((NA.LE.0) .OR. (NB.LE.0)) THEN
         RETURN
      ENDIF

      IF (LDIAGO) THEN

         IF (NA .NE. NB) THEN
            WRITE(LUPRI,'(//A,A,A/A,I10/A,I10/A)')
     &      ' Error in ',SECNAM,':',
     &      ' NA =',NA,' NB =',NB,
     &      ' NA = NB must hold for LDIAGO option.'
            CALL QUIT('NA != NB in '//SECNAM)
         ENDIF

         IF (NLROP+NA .GT. MXLROP) THEN
             WRITE(LUPRI,'(2(/A,I5))')
     &       ' NO. OF OPERATOR DOUBLES SPECIFIED  : ',NLROP+NA,
     &       ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLROP
             CALL QUIT('TOO MANY OPERATOR DOUBLES IN '//SECNAM)
         ENDIF

         DO I = 1,NA
            IOPA(I) = INDPRPCC(LABELA(I))
            IOPB(I) = INDPRPCC(LABELB(I))
         ENDDO

         DO IDXAB = 1,NA
            IDX = NLROP + IDXAB
            IALROP(IDX) = IOPA(IDXAB)
            IBLROP(IDX) = IOPB(IDXAB)
            LALORX(IDX) = LRELAX
            LBLORX(IDX) = LRELAX
         ENDDO

         NLROP = NLROP + NA

      ELSE

         NTOT = NA*NB

         IF (NLROP+NTOT .GT. MXLROP) THEN
             WRITE(LUPRI,'(2(/A,I5))')
     &       ' NO. OF OPERATOR DOUBLES SPECIFIED  : ',NLROP+NTOT,
     &       ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLROP
             CALL QUIT('TOO MANY OPERATOR DOUBLES IN '//SECNAM)
         ENDIF

         DO I = 1,NA
            IOPA(I) = INDPRPCC(LABELA(I))
         ENDDO
         DO I = 1,NB
            IOPB(I) = INDPRPCC(LABELB(I))
         ENDDO

         DO IDXA=1,NA
            DO IDXB=1,NB
               IDX = NLROP + (IDXA - 1)*NB + IDXB
               IALROP(IDX) = IOPA(IDXA)
               IBLROP(IDX) = IOPB(IDXB)
               LALORX(IDX) = LRELAX
               LBLORX(IDX) = LRELAX
            ENDDO
         ENDDO

         NLROP = NLROP + NTOT

      ENDIF

      RETURN
      END
c /* deck cc_qrinp */
*=====================================================================*
       SUBROUTINE CC_QRINP(WORD)
*---------------------------------------------------------------------*
*
*    Purpose: read input for CC dynamic first hyperpolarizabilities
*             and dispersion coefficients
*
*    if (WORD .eq '*CCQR  ') read & process input and set defaults, 
*    else set only defaults 
*
*    Written by Christof Haettig, October 1996, modified December '96
*    dispersion coefficients, October 1997 (Christof Haettig)
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccqrinf.h>

* local parameters:
      CHARACTER MSGDBG*(18)
      PARAMETER (MSGDBG='[debug] CC_QRINP> ')
      CHARACTER SECNAM*(8)
      PARAMETER (SECNAM='CC_QRINP')

      INTEGER NTABLE
      PARAMETER (NTABLE = 20)

#if defined (SYS_CRAY)
      REAL ZERO
#else
      DOUBLE PRECISION ZERO
#endif
      PARAMETER (ZERO = 0.0d00)


* variables:
      LOGICAL SET
      SAVE SET

      CHARACTER WORD*(7), LINE*(80)
      CHARACTER*8 LABELA, LABELB, LABELC
      CHARACTER TABLE(NTABLE)*(7)

      LOGICAL LALRX, LBLRX, LCLRX, LRELAX
      INTEGER IDX, IJUMP, K, M, N
      INTEGER MFREQ
      INTEGER IFREQ, ICA, ICB, ICC, IDXA, IDXB, IDXC, IDIP(3)

      DATA SET /.FALSE./

      DATA TABLE /'.OPERAT','.DIPOLE','.PRINT ','.EXPCOF','.AVERAG',
     &            '.MIXFRE','.SHGFRE','.ORFREQ','.EOPEFR','.STATIC',
     &            '.DISPCF','.ALLDSP','.XYDEGE','.NOBMAT','.USE R2',
     &            '.RELAXE','.UNRELA','.USE AA','.AVANEW','.XXXXXX' /

      INTEGER INDPRPCC

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (SET) RETURN
      SET = .TRUE.

*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*

      NQROPER = 0
      NQRFREQ = 0
      NQRDISP = 0
      NQRDSPE = 0
      NQRDSPO = 0

      CCQR         = .FALSE.
      BETA_AVERAGE = .FALSE.
      XY_DEGENERAT = .FALSE.
      USEBTRAN     = .TRUE.
      USE_R2       = .FALSE.
      USE_AAMAT    = .FALSE.
      ALLDSPCF     = .FALSE.
      LALRX        = .FALSE.
      LBLRX        = .FALSE.
      LCLRX        = .FALSE.
      LRELAX       = .FALSE.    
      LAVANEW      = .FALSE.

      IPRQHYP = 0

      ICHANG = 0

*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCQR  ') THEN

100   CONTINUE
       
* get new input line:
        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
          CALL UPCASE(WORD)
        END DO

        IF (WORD(1:1) .EQ. '.') THEN
C         WRITE (LUPRI,*) WORD
C         CALL FLSHFO(LUPRI)

c         table look up:
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO
     
c         jump to the appropiate input section:
          IF (IJUMP .LE. NTABLE) THEN
            ICHANG = ICHANG + 1
            GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), 
     &           IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_QRINP.')
      
C           ------------------------------------------
C           .OPERAT : triples of operator lables A,B,C
C           ------------------------------------------
1           CONTINUE
              READ (LUCMD,'(3A)') LABELA, LABELB, LABELC
              DO WHILE (LABELA(1:1).NE.'.' .AND. LABELA(1:1).NE.'*')
                IF      (LABELA(1:1).EQ.'!') THEN
                  CONTINUE
                ELSE IF (LABELA(1:1).EQ.'(') THEN
                    LALRX = .FALSE.
                    LBLRX = .FALSE.
                    LCLRX = .FALSE.
                    IF (LABELA(1:7).EQ.'(RELAX)') LALRX = .TRUE.
                    IF (LABELB(1:7).EQ.'(RELAX)') LBLRX = .TRUE.
                    IF (LABELC(1:7).EQ.'(RELAX)') LCLRX = .TRUE.
                    IF (LALRX .OR. LBLRX .OR. LCLRX) THEN
                      KEEPAOTWO = MAX(KEEPAOTWO,1)
                    END IF
                ELSE
                  IF (NQROPER.LT.MXQROP) THEN
                    NQROPER = NQROPER + 1
                    IAQROP(NQROPER) = INDPRPCC(LABELA)
                    IBQROP(NQROPER) = INDPRPCC(LABELB)
                    ICQROP(NQROPER) = INDPRPCC(LABELC)
                    LAQLRX(NQROPER) = LALRX
                    LBQLRX(NQROPER) = LBLRX
                    LCQLRX(NQROPER) = LCLRX
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)')
     &               ' NO. OF OPERATOR TRIPLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQROP
                    CALL QUIT('TOO MANY OPERATOR TRIPLES IN CCQR.')
                  END IF
                END IF
                READ (LUCMD,'(3A)') LABELA, LABELB, LABELC
              END DO
              BACKSPACE(LUCMD)
            GO TO 100

C           -------------------------------------------------------
C           .DIPOL : calculate complete dipole-dipole-dipole tensor
C           -------------------------------------------------------
2           CONTINUE
              IF (NQROPER+27 .GT. MXQROP) THEN
                WRITE(LUPRI,'(2(/A,I5))') 
     &          ' NO. OF OPERATOR TRIPLES SPECIFIED  : ',NQROPER+27,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQROP
                CALL QUIT('TOO MANY OPERATOR TRIPLES IN CCQR.')
              END IF
              IDIP(1) = INDPRPCC('XDIPLEN ')
              IDIP(2) = INDPRPCC('YDIPLEN ')
              IDIP(3) = INDPRPCC('ZDIPLEN ')
              DO IDXA=1,3
              DO IDXB=1,3
              DO IDXC=1,3
                IDX = NQROPER + (IDXA-1)*9+(IDXB-1)*3+IDXC
                IAQROP(IDX) = IDIP(IDXA)
                IBQROP(IDX) = IDIP(IDXB)
                ICQROP(IDX) = IDIP(IDXC)
                LAQLRX(IDX) = LRELAX
                LBQLRX(IDX) = LRELAX
                LCQLRX(IDX) = LRELAX      
              END DO
              END DO
              END DO
              NQROPER = NQROPER + 27
            GO TO 100

C           ------------
C           .PRINT
C           ------------
3           CONTINUE
              READ (LUCMD,*) IPRQHYP
            GO TO 100

C           -----------------------------------------------------------
C           .EXPCOF : coefficients for the expansion of
C                     <<A;B,C>>_{w_B,w_C} in the frequenies w_B and w_C
C           -----------------------------------------------------------
4           CONTINUE
              READ (LUCMD,'(A)') LINE
              DO WHILE (LINE(1:1).NE.'.' .AND. LINE(1:1).NE.'*')
                IF (LINE(1:1).NE.'!') THEN
                  IF (NQRDISP.LT.MXQRDISP) THEN
                    READ(LINE,*) ICA, ICB, ICC
                    IF (ICA.LT.0 .OR. ICB.LT.0 .OR. ICC.LT.0) THEN
                      NWARN = NWARN + 1
                      WRITE(LUPRI,'(/2A/A)')
     &                 '@ WARNING: NEGATIVE EXPANSION COEFFICIENTS NOT',
     &                 ' AVAILABLE FOR FIRST HYPERPOLARIZABILITIES.',
     &                 '@ WARNING: INPUT LINE IGNORED...'
                    ELSE
                      NQRDISP = NQRDISP + 1
                      IQCAUA(NQRDISP) = ICA
                      IQCAUB(NQRDISP) = ICB
                      IQCAUC(NQRDISP) = ICC
                    END IF
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)') 
     &               ' NO. OF EXPANSION COEFFICIENTS ',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQRDISP
                    CALL QUIT('TOO MANY EXPANSION COEFFICIENTS IN CCQR')
                  END IF
                END IF
                READ (LUCMD,'(A)') LINE
              END DO
              BACKSPACE(LUCMD)
            GO TO 100

C           -------------------------------------------------------
C           .AVERAG : calculate averaged tensor components
C                     implemented: beta_{||}, beta_{_|_}, beta_{ms}
C           -------------------------------------------------------
5           CONTINUE
              READ (LUCMD,'(A)') LINE
              IF (LINE(1:8).EQ.'HYPERPOL') THEN
                IF (NQROPER.NE.0) THEN
                  NWARN = NWARN + 1
                  WRITE(LUPRI,'(/2A/A/)') 
     &             '@ WARNING: INPUT FOR .DIPOL OR .OPERATOR OPTIONS',
     &             ' BEFORE THE .AVERAG OPTION',
     &             '@ IN *CCQR SECTION WILL BE IGNORED.'
                  NQROPER = 0
                END IF
                IDIP(1) = INDPRPCC('XDIPLEN ')
                IDIP(2) = INDPRPCC('YDIPLEN ')
                IDIP(3) = INDPRPCC('ZDIPLEN ')
                DO IDX=1,2
                  IAQROP(1)     = IDIP(3)    ! component 1: beta_{zzz}
                  IBQROP(1)     = IDIP(3)
                  ICQROP(1)     = IDIP(3)
                  LAQLRX(1)     = LRELAX
                  LBQLRX(1)     = LRELAX
                  LCQLRX(1)     = LRELAX   

                  IAQROP(2+(IDX-1)*3) = IDIP(3)   ! comp. 2: beta_{zxx}
                  IBQROP(2+(IDX-1)*3) = IDIP(IDX) ! comp. 5: beta_{zyy}
                  ICQROP(2+(IDX-1)*3) = IDIP(IDX)
                  LAQLRX(2+(IDX-1)*3) = LRELAX
                  LBQLRX(2+(IDX-1)*3) = LRELAX
                  LCQLRX(2+(IDX-1)*3) = LRELAX     

                  IAQROP(3+(IDX-1)*3) = IDIP(IDX) ! comp. 3: beta_{xzx}
                  IBQROP(3+(IDX-1)*3) = IDIP(3)   ! comp. 6: beta_{yzy}
                  ICQROP(3+(IDX-1)*3) = IDIP(IDX)
                  LAQLRX(3+(IDX-1)*3) = LRELAX
                  LBQLRX(3+(IDX-1)*3) = LRELAX
                  LCQLRX(3+(IDX-1)*3) = LRELAX     

                  IAQROP(4+(IDX-1)*3) = IDIP(IDX) ! comp. 4: beta_{xxz}
                  IBQROP(4+(IDX-1)*3) = IDIP(IDX) ! comp. 7: beta_{yyz}
                  ICQROP(4+(IDX-1)*3) = IDIP(3)
                  LAQLRX(4+(IDX-1)*3) = LRELAX
                  LBQLRX(4+(IDX-1)*3) = LRELAX
                  LCQLRX(4+(IDX-1)*3) = LRELAX     
                END DO
                NQROPER = 7
                BETA_AVERAGE = .TRUE.
                IF (XY_DEGENERAT) THEN
                  ! forget beta_{zyy}, beta_{yzy}, beta_{yyz}
                  NQROPER = 4
                END IF
              ELSE
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/4A/A/)') '@ WARNING: Label "',
     &           LINE(1:8),'" unknown for .AVERAG keyword',
     &           'in *CCQR section.','@ WARNING: INPUT IS IGNORED...'
              ENDIF 
            GO TO 100


C           ---------------------------
C           .MIXFRE : mixed frequencies
C                     wb, wc, wa=-wb-wc
C           ---------------------------
6           CONTINUE
              READ (LUCMD,*) MFREQ
              IF (NQRFREQ+MFREQ .GT. MXQRFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5),/)') '@ WARNING:',
     &          '@ NUMBER OF FREQUENCIES SPECIFIED    : ',NQRFREQ+MFREQ,
     &          '@ IS GREATER THAN THE ALLOWED NUMBER : ',MXQRFR,
     &          '@ THE NUMBER IS RESET TO THE MAXIMUM : ',MXQRFR
                MFREQ = MXQRFR-NQRFREQ
              END IF
              READ (LUCMD,*) (BQRFR(IDX),IDX=NQRFREQ+1,NQRFREQ+MFREQ)
              READ (LUCMD,*) (CQRFR(IDX),IDX=NQRFREQ+1,NQRFREQ+MFREQ)
              NQRFREQ = NQRFREQ + MFREQ
            GO TO 100

C           ------------------------------------------------
C           .SHGFRE : second harmonic generation frequencies 
C                     wb, wc = wb, wa = -2wb
C           ------------------------------------------------
7           CONTINUE
              READ (LUCMD,*) MFREQ
              IF (NQRFREQ+MFREQ .GT. MXQRFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5),/)') '@ WARNING:',
     &          '@ NUMBER OF FREQUENCIES SPECIFIED    : ',NQRFREQ+MFREQ,
     &          '@ IS GREATER THAN THE ALLOWED NUMBER : ',MXQRFR,
     &          '@ THE NUMBER IS RESET TO THE MAXIMUM : ',MXQRFR
                MFREQ = MXQRFR-NQRFREQ
              END IF
              READ (LUCMD,*) (BQRFR(IDX),IDX=NQRFREQ+1,NQRFREQ+MFREQ)
              DO IDX = NQRFREQ+1, NQRFREQ+MFREQ
                CQRFR(IDX) = BQRFR(IDX)
              END DO
              NQRFREQ = NQRFREQ + MFREQ
            GO TO 100

C           ------------------------------------------------
C           .ORFREQ : optical rectification frequencies
C                     wb, wc = -wb, wa = 0
C           ------------------------------------------------
8           CONTINUE
              READ (LUCMD,*) MFREQ
              IF (NQRFREQ+MFREQ .GT. MXQRFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5),/)') '@ WARNING:',
     &          '@ NUMBER OF FREQUENCIES SPECIFIED    : ',NQRFREQ+MFREQ,
     &          '@ IS GREATER THAN THE ALLOWED NUMBER : ',MXQRFR,
     &          '@ THE NUMBER IS RESET TO THE MAXIMUM : ',MXQRFR
                MFREQ = MXQRFR-NQRFREQ
              END IF
              READ (LUCMD,*) (BQRFR(IDX),IDX=NQRFREQ+1,NQRFREQ+MFREQ)
              DO IDX = NQRFREQ+1, NQRFREQ+MFREQ
                CQRFR(IDX) = -BQRFR(IDX)
              END DO
              NQRFREQ = NQRFREQ + MFREQ
            GO TO 100

C           ------------------------------------------------
C           .EOPEFR : second harmonic generation frequencies 
C                     wb, wc = 0 , wa = -wb
C           ------------------------------------------------
9           CONTINUE
              READ (LUCMD,*) MFREQ
              IF (NQRFREQ+MFREQ .GT. MXQRFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5),/)') '@ WARNING:',
     &          '@ NUMBER OF FREQUENCIES SPECIFIED    : ',NQRFREQ+MFREQ,
     &          '@ IS GREATER THAN THE ALLOWED NUMBER : ',MXQRFR,
     &          '@ THE NUMBER IS RESET TO THE MAXIMUM : ',MXQRFR
                MFREQ = MXQRFR-NQRFREQ
              END IF
              READ (LUCMD,*) (BQRFR(IDX),IDX=NQRFREQ+1,NQRFREQ+MFREQ)
              DO IDX = NQRFREQ+1, NQRFREQ+MFREQ
                CQRFR(IDX) = ZERO
              END DO
              NQRFREQ = NQRFREQ + MFREQ
            GO TO 100

C           ---------------------------------------------------
C           .STATIC : add wb = wc = wa = zero to frequency list
C           ---------------------------------------------------
10          CONTINUE
              IF (NQRFREQ+1 .GT. MXQRFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5),/)') '@ WARNING:',
     &          '@ NUMBER OF FREQUENCIES SPECIFIED    : ',NQRFREQ+1,
     &          '@ IS GREATER THAN THE ALLOWED NUMBER : ',MXQRFR,
     &          '@ INPUT OPTION .STATIC WILL BE IGNORED.'
              ELSE
                NQRFREQ = NQRFREQ + 1
                BQRFR(NQRFREQ) = ZERO
                CQRFR(NQRFREQ) = ZERO
              END IF
            GO TO 100

C           ---------------------------------
C           .DISPCF : dispersion coefficients
C           ---------------------------------
11          CONTINUE
              READ (LUCMD,*) NQRDSPE
              !WRITE (LUPRI,*) 'NQRDSPE = ',NQRDSPE
              IF (NQRDISP.NE.0) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/2A)') 
     &            '@ WARNING: INPUT FOR .EXPCOF OPTION BEFORE .DISPCF',
     &            ' IN *CCQR SECTION WILL BE IGNORED.'
                NQRDISP = 0
              END IF
              DO K = 0, NQRDSPE
                ! WRITE (LUPRI,*) 'NQRDSPE,K = ',NQRDSPE,K
                IF ((NQRDISP+(K+2)*(K+1)/2).LE.MXQRDISP) THEN
                  DO M = 0, K, 1
                    DO N = 0, M, 1
                      NQRDISP = NQRDISP + 1
                      IQCAUA(NQRDISP) = K-M
                      IQCAUB(NQRDISP) = M-N
                      IQCAUC(NQRDISP) = N
                    END DO
                  END DO
                ELSE
                  NWARN = NWARN + 1
                  WRITE(LUPRI,'(/A/2A,I5/)') '@ WARNING:',
     &            '@ NO. OF DISPERSION COEFFICIENTS NEEDED',
     &            ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQRDISP
                  WRITE(LUPRI,'(/A,I3,A)')
     &            '@ DISPERSION COEFFICIENTS OF ORDER',K,' ARE IGNORED.'
                  NQRDSPE = MIN(NQRDSPE,K-1)
                END IF
              END DO
            GO TO 100

C           --------------------------------------------------------
C           .ALLDSP : do not skip odd/even dispersion coefficients
C                     or real/imaginary properties
C           --------------------------------------------------------
12          CONTINUE
              ALLDSPCF = .TRUE.
            GO TO 100

C           ----------------------------------------------
C           .XYDEGE : assume X and Y directions degenerate
C           ----------------------------------------------
13          CONTINUE
              XY_DEGENERAT = .TRUE.
              IF ( BETA_AVERAGE .AND. NQROPER.EQ.7 ) THEN
                ! forget beta_{zyy}, beta_{yzy}, beta_{yyz}
                NQROPER = 4
              END IF
            GO TO 100

C           ---------------------------------------------------------
C           .NOBMAT : don't use B matrix transformation but F matrix
C                     (usually less efficient, because less symmetry)
C           ---------------------------------------------------------
14          CONTINUE
              USEBTRAN = .FALSE.
            GO TO 100

C           -----------------------------------------------------------
C           .USE R2 : use second-order response/Cauchy vectors R2/CR2 
C                     instead first-order left L1/LC vectors times 
C                     B matrix transf. and eta vectors
C                     (test option, computational advantages only in
C                      very rare cases...)
C           -----------------------------------------------------------
15          CONTINUE
              USEBTRAN = .FALSE.
              USE_R2   = .TRUE.
            GO TO 100

C           ----------------------------------------------------------
C           .RELAXE : switch to relaxed modus for all three operators:
C           ----------------------------------------------------------
16          CONTINUE
              ! LRELAX    = .TRUE.
              ! KEEPAOTWO = MAX(KEEPAOTWO,1)
              WRITE (LUPRI,*) 
     *            '.RELAXE keyword in *CCQR section is disabled.'
            GO TO 100

C           ------------------------------------------------------------
C           .UNRELA : switch to unrelaxed modus for all three operators:
C           ------------------------------------------------------------
17          CONTINUE
              LRELAX = .FALSE.
            GO TO 100

C           ------------------------------------------------------------
C           .USE AA : Use A{O} transformation instead of Eta{O} vectors:
C           ------------------------------------------------------------
18          CONTINUE
              USE_AAMAT = .TRUE.
            GO TO 100
C
C           ------------------------------------------------------------
C           .AVANEW: Calculates:
C           beta_i = 1/3 Sum_j=x,y,z [ B_ijj + B_jji + B_jij ]
C           |beta_i*mu_i| for (i=x,y,z) (mu is the dipole moment)
C           <beta>=1/6( B_xyz - B_xzy + B_yzx - B_yxz + B_zxy - B_zyx )
C           ------------------------------------------------------------
19          CONTINUE
              LAVANEW = .TRUE.
            GO TO 100 
C           -----------------------
C           .XXXXXX : unused labels
C           -----------------------
20          CONTINUE
              WRITE (LUPRI,*) 'unused .XXXXXX label... ignored'
            GO TO 100

          ELSE
           WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF

      END IF

200   CONTINUE

*---------------------------------------------------------------------*
* check, if any triples of operator labels specified:
* if not, use default: complete dipole-dipole-dipole tensor
*---------------------------------------------------------------------*
      IF (NQROPER .EQ. 0) THEN 
          IF (NQROPER+27 .GT. MXQROP) THEN
            WRITE(LUPRI,'(2(/A,I5))') 
     &      ' NO. OF OPERATOR TRIPLES SPECIFIED  : ',NQROPER+27,
     &      ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQROP
            CALL QUIT('TOO MANY OPERATOR TRIPLES IN CCQR.')
          END IF
          IDIP(1) = INDPRPCC('XDIPLEN ')
          IDIP(2) = INDPRPCC('YDIPLEN ')
          IDIP(3) = INDPRPCC('ZDIPLEN ')
          DO IDXA=1,3
          DO IDXB=1,3
          DO IDXC=1,3
            IDX = NQROPER + (IDXA-1)*9+(IDXB-1)*3+IDXC
            IAQROP(IDX) = IDIP(IDXA)
            IBQROP(IDX) = IDIP(IDXB)
            ICQROP(IDX) = IDIP(IDXC)
          END DO
          END DO
          END DO
          NQROPER = NQROPER + 27
      END IF

*---------------------------------------------------------------------*
* check, if frequencies specified; if not, use default: static
*---------------------------------------------------------------------*
      IF (NQRFREQ .EQ. 0) THEN 
        NQRFREQ = NQRFREQ + 1
        BQRFR(NQRFREQ) = ZERO
        CQRFR(NQRFREQ) = ZERO
      END IF

*---------------------------------------------------------------------*
* add list with wa frequencies:
*---------------------------------------------------------------------*
      DO IFREQ = 1, NQRFREQ
        AQRFR(IFREQ) = - ( BQRFR(IFREQ) + CQRFR(IFREQ) )
      END DO

*---------------------------------------------------------------------*
* set CCQR flags:
*---------------------------------------------------------------------*
      CCQR  = .TRUE.

      RETURN
      END
*---------------------------------------------------------------------*
c /* deck cc_crinp */
*=====================================================================*
       SUBROUTINE CC_CRINP(WORD)
*---------------------------------------------------------------------*
*
*    Purpose: read input for CC dynamic second hyperpolarizabilities
*
*    if (WORD .eq '*CCCR  ') read & process input and set defaults, 
*    else set only defaults 
*
*    Written by Christof Haettig, October 1996, modified Februar '97
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <cccrinf.h>
C#include <ccrspprp.h>

* local parameters:
      CHARACTER MSGDBG*(18)
      PARAMETER (MSGDBG='[debug] CC_CRINP> ')
      CHARACTER SECNAM*(8)
      PARAMETER (SECNAM='CC_CRINP')

      INTEGER NTABLE
      PARAMETER (NTABLE = 20)

#if defined (SYS_CRAY)
      REAL ZERO
#else
      DOUBLE PRECISION ZERO
#endif
      PARAMETER (ZERO = 0.0d00)


* variables:
      LOGICAL SET
      SAVE SET

      CHARACTER*(7)  WORD
      CHARACTER*(80) LINE
      CHARACTER*(7)  TABLE(NTABLE)
      CHARACTER*(8)  LABELA, LABELB, LABELC, LABELD

      INTEGER IDX, IJUMP, IFREQ, IDIP(3), IDXA, IDXB, IDXC, IDXD
      INTEGER MFREQ, K, L, M, N, ICAUA, ICAUB, ICAUC, ICAUD

      DATA SET /.FALSE./
      DATA TABLE /'.OPERAT','.DIPOLE','.PRINT ','.STATIC','.MIXFRE',
     &            '.THGFRE','.ESHGFR','.DFWMFR','.DCKERR','.USECHI',
     &            '.USEXKS','.EXPCOF','.AVERAG','.DISPCF','.ODDISP',
     &            '.NO2NP1','.L2 BCD','.L2 BC ','.XXXXXX','.XXXXXX'/

      INTEGER INDPRPCC

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (SET) RETURN
      SET = .TRUE.

*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*
      NCROPER = 0
      NCRFREQ = 0
      NCRDISP = 0
      NCRDSPE = -1
      NCRDSPO = -1

      CCCR      = .FALSE.
      GAMMA_PAR = .FALSE.
      GAMMA_ORT = .FALSE.
      CSYM      = 'GENERI'

      L_USE_CHI2   = .FALSE.
      L_USE_XKS3   = .FALSE.
      NO_2NP1_RULE = .FALSE.
      USE_L2BC     = .FALSE.
      USE_LBCD     = .FALSE.

      IPRCHYP = IPRINT

*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCCR  ') THEN

100   CONTINUE
       
* get new input line:
        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
          CALL UPCASE(WORD)
        END DO

        IF (WORD(1:1) .EQ. '.') THEN

c         table look up:
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO
     
c         jump to the appropiate input section:
          IF (IJUMP .LE. NTABLE) THEN
            ICHANG = ICHANG + 1
            GOTO ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &            11,12,13,14,15,16,17,18,19,20), IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_CRINP.')
      
C           ------------
C           .OPERAT
C           ------------
1           CONTINUE
              READ (LUCMD,'(4A)') LABELA,LABELB,LABELC,LABELD
              DO WHILE (LABELA(1:1).NE.'.' .AND. LABELA(1:1).NE.'*')
                IF (LABELA(1:1).NE.'!') THEN
                  IF (NCROPER.LT.MXCROP) THEN
                    NCROPER = NCROPER + 1
                    IACROP(NCROPER) = INDPRPCC(LABELA)
                    IBCROP(NCROPER) = INDPRPCC(LABELB)
                    ICCROP(NCROPER) = INDPRPCC(LABELC)
                    IDCROP(NCROPER) = INDPRPCC(LABELD)
                  ELSE
                    WRITE(LUPRI,'(/2A,I5/)')
     &               ' NO. OF OPERATOR QUADRUPLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXCROP
                    CALL QUIT('TOO MANY OPERATOR QUADRUPLES IN CCCR.')
                  END IF
                END IF
                READ (LUCMD,'(4A)') LABELA,LABELB,LABELC,LABELD
              END DO
              BACKSPACE(LUCMD)
            GO TO 100

C           -------------------------------------------------------
C           .DIPOLE: calculate complete dipole-dipole-dipole-dipole 
C           -------------------------------------------------------
2           CONTINUE
              IF (NCROPER+81 .GT. MXCROP) THEN
                WRITE(LUPRI,'(/2A,I5/)')
     &           ' NO. OF OPERATOR QUADRUPLES SPECIFIED',
     &           ' IS GREATER THAN THE ALLOWED NUMBER : ',MXCROP
                CALL QUIT('TOO MANY OPERATOR QUADRUPLES IN CCCR.')
              END IF
              IDIP(1) = INDPRPCC('XDIPLEN ')
              IDIP(2) = INDPRPCC('YDIPLEN ')
              IDIP(3) = INDPRPCC('ZDIPLEN ')
              DO IDXA=1,3
              DO IDXB=1,3
              DO IDXC=1,3
              DO IDXD=1,3
                IDX = NCROPER + (IDXA-1)*27+(IDXB-1)*9+(IDXC-1)*3+IDXD
                IACROP(IDX) = IDIP(IDXA)
                IBCROP(IDX) = IDIP(IDXB)
                ICCROP(IDX) = IDIP(IDXC)
                IDCROP(IDX) = IDIP(IDXD)
              END DO
              END DO
              END DO
              END DO
              NCROPER = NCROPER + 81
            GO TO 100

C           ------------
C           .PRINT
C           ------------
3           CONTINUE
              READ (LUCMD,*) IPRCHYP
            GO TO 100

C           ------------
C           .STATIC
C           ------------
4           CONTINUE
              IF (NCRFREQ+1 .GT. MXCRFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5)/)') '@ WARNING:',
     &          '@ NUMBER OF FREQUENCIES SPECIFIED    : ',NCRFREQ+1,
     &          '@ IS GREATER THAN THE ALLOWED NUMBER : ',MXCRFR,
     &          '@ INPUT OPTION .STATIC WILL BE IGNORED.'
              ELSE
                NCRFREQ = NCRFREQ + 1
                BCRFR(NCRFREQ) = ZERO
                CCRFR(NCRFREQ) = ZERO
                DCRFR(NCRFREQ) = ZERO
              END IF
            GO TO 100

C           ------------------------------------------------
C           .MIXFRE : mixed frequency input:
C                     read wb, wc, wd  --->  wa = -wb-wc-wd
C           ------------------------------------------------
5           CONTINUE
              READ (LUCMD,*) MFREQ
              IF (NCRFREQ+MFREQ .GT. MXCRFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5)/)') '@ WARNING:',
     &          '@ NUMBER OF FREQUENCIES SPECIFIED    : ',NCRFREQ+MFREQ,
     &          '@ IS GREATER THAN THE ALLOWED NUMBER : ',MXCRFR,
     &          '@ THE NUMBER IS RESET TO THE MAXIMUM : ',MXCRFR
                MFREQ = MXCRFR-NCRFREQ
              END IF
              READ (LUCMD,*) (BCRFR(IDX),IDX=NCRFREQ+1,NCRFREQ+MFREQ)
              READ (LUCMD,*) (CCRFR(IDX),IDX=NCRFREQ+1,NCRFREQ+MFREQ)
              READ (LUCMD,*) (DCRFR(IDX),IDX=NCRFREQ+1,NCRFREQ+MFREQ)
              NCRFREQ = NCRFREQ + MFREQ
            GO TO 100

C           -----------------------------------------------
C           .THGFRE : third harmonic generation frequencies
C                     read wb --> wc=wb, wd=wb, wa= -3wb
C           -----------------------------------------------
6           CONTINUE
              READ (LUCMD,*) MFREQ
              IF (NCRFREQ+MFREQ .GT. MXCRFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5)/)') '@ WARNING:',
     &          '@ NUMBER OF FREQUENCIES SPECIFIED    : ',NCRFREQ+MFREQ,
     &          '@ IS GREATER THAN THE ALLOWED NUMBER : ',MXCRFR,
     &          '@ THE NUMBER IS RESET TO THE MAXIMUM : ',MXCRFR
                MFREQ = MXCRFR-NCRFREQ
              END IF
              READ (LUCMD,*) (BCRFR(IDX),IDX=NCRFREQ+1,NCRFREQ+MFREQ)
              DO IDX = NCRFREQ+1, NCRFREQ+MFREQ
                CCRFR(IDX) = BCRFR(IDX)
                DCRFR(IDX) = BCRFR(IDX)
              END DO
              NCRFREQ = NCRFREQ + MFREQ
            GO TO 100


C           -----------------------------------------------------------
C           .ESHGFR : electric field induced second harmonic generation
C                     read wb --> wc=wb, wd=0, wa= -2wb
C           -----------------------------------------------------------
7           CONTINUE
              READ (LUCMD,*) MFREQ
              IF (NCRFREQ+MFREQ .GT. MXCRFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5)/)') '@ WARNING:',
     &          '@ NUMBER OF FREQUENCIES SPECIFIED    : ',NCRFREQ+MFREQ,
     &          '@ IS GREATER THAN THE ALLOWED NUMBER : ',MXCRFR,
     &          '@ THE NUMBER IS RESET TO THE MAXIMUM : ',MXCRFR
                MFREQ = MXCRFR-NCRFREQ
              END IF
              READ (LUCMD,*) (BCRFR(IDX),IDX=NCRFREQ+1,NCRFREQ+MFREQ)
              DO IDX = NCRFREQ+1, NCRFREQ+MFREQ
                CCRFR(IDX) = BCRFR(IDX)
                DCRFR(IDX) = ZERO
              END DO
              NCRFREQ = NCRFREQ + MFREQ
            GO TO 100


C           -----------------------------------------------------------
C           .DFWMFR : degenerate four wave mixing
C                     read wb --> wc=+wb, wd=-wb, wa= -wb
C           -----------------------------------------------------------
8           CONTINUE
              READ (LUCMD,*) MFREQ
              IF (NCRFREQ+MFREQ .GT. MXCRFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5)/)') '@ WARNING:',
     &          '@ NUMBER OF FREQUENCIES SPECIFIED    : ',NCRFREQ+MFREQ,
     &          '@ IS GREATER THAN THE ALLOWED NUMBER : ',MXCRFR,
     &          '@ THE NUMBER IS RESET TO THE MAXIMUM : ',MXCRFR
                MFREQ = MXCRFR-NCRFREQ
              END IF
              READ (LUCMD,*) (BCRFR(IDX),IDX=NCRFREQ+1,NCRFREQ+MFREQ)
              DO IDX = NCRFREQ+1, NCRFREQ+MFREQ
                CCRFR(IDX) = +BCRFR(IDX)
                DCRFR(IDX) = -BCRFR(IDX)
              END DO
              NCRFREQ = NCRFREQ + MFREQ
            GO TO 100

C           -----------------------------------------------------------
C           .DCKERR : dc Kerr effect, also optical Kerr effect (OKE)
C                     read wd --> wc=wd=0,  wa= -wd
C           -----------------------------------------------------------
9           CONTINUE
              READ (LUCMD,*) MFREQ
              IF (NCRFREQ+MFREQ .GT. MXCRFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5)/)') '@ WARNING:',
     &          '@ NUMBER OF FREQUENCIES SPECIFIED    : ',NCRFREQ+MFREQ,
     &          '@ IS GREATER THAN THE ALLOWED NUMBER : ',MXCRFR,
     &          '@ THE NUMBER IS RESET TO THE MAXIMUM : ',MXCRFR
                MFREQ = MXCRFR-NCRFREQ
              END IF
              READ (LUCMD,*) (DCRFR(IDX),IDX=NCRFREQ+1,NCRFREQ+MFREQ)
              DO IDX = NCRFREQ+1, NCRFREQ+MFREQ
                BCRFR(IDX) = ZERO
                CCRFR(IDX) = ZERO
              END DO
              NCRFREQ = NCRFREQ + MFREQ
            GO TO 100


C           -------------------------------------------------------
C           .USECHI : use second-order chi vectors as intermediates
C                     (test option)
C           -------------------------------------------------------
10          CONTINUE
              L_USE_CHI2 = .TRUE.
              IF (L_USE_XKS3) THEN
                L_USE_XKS3 = .FALSE.
                WRITE(LUPRI,*) '.USECHI and .USEXKS are incompatible'
                WRITE(LUPRI,*) 'in the *CCCR section...'
                WRITE(LUPRI,*) '.USEXKS is switched off'
              END IF
              IF (USE_LBCD) THEN
                USE_LBCD = .FALSE.
                WRITE(LUPRI,*) '.L2 BCD and .USECHI are incompatible'
                WRITE(LUPRI,*) 'in the *CCCR section...'
                WRITE(LUPRI,*) '.L2 BCD is switched off'
              END IF
              IF (USE_L2BC) THEN
                USE_L2BC = .FALSE.
                WRITE(LUPRI,*) '.L2 BC  and .USECHI are incompatible'
                WRITE(LUPRI,*) 'in the *CCCR section...'
                WRITE(LUPRI,*) '.L2 BC  is switched off'
              END IF
            GO TO 100

C           -------------------------------------------------------
C           .USEXKS : use third-order xksi vectors as intermediates
C                     (test option)
C           -------------------------------------------------------
11          CONTINUE
              L_USE_XKS3 = .TRUE.
              IF (L_USE_CHI2) THEN
                L_USE_CHI2 = .FALSE.
                WRITE(LUPRI,*) '.USECHI and .USEXKS are incompatible'
                WRITE(LUPRI,*) 'in the *CCCR section...'
                WRITE(LUPRI,*) '.USECHI is switched off'
              END IF
              IF (USE_LBCD) THEN
                USE_LBCD = .FALSE.
                WRITE(LUPRI,*) '.L2 BCD and .USEXKS are incompatible'
                WRITE(LUPRI,*) 'in the *CCCR section...'
                WRITE(LUPRI,*) '.L2 BCD is switched off'
              END IF
              IF (USE_L2BC) THEN
                USE_L2BC = .FALSE.
                WRITE(LUPRI,*) '.L2 BC  and .USEXKS are incompatible'
                WRITE(LUPRI,*) 'in the *CCCR section...'
                WRITE(LUPRI,*) '.L2 BC  is switched off'
              END IF
            GO TO 100

C           -----------------------------------------------------------
C           .EXPCOF : coefficients for the expansion of
C             <<A;B,C,D>>_{w_B,w_C,w_D} in the frequenies w_B, w_C, w_D
C           -----------------------------------------------------------
12          CONTINUE
              READ (LUCMD,'(A)') LINE
              DO WHILE (LINE(1:1).NE.'.' .AND. LINE(1:1).NE.'*')
                IF (LINE(1:1).NE.'!') THEN
                  IF (NCRDISP.LT.MXCRDISP) THEN
                    READ(LINE,*) ICAUA, ICAUB, ICAUC, ICAUD
                    IF (ICAUA.LT.0 .OR. ICAUB.LT.0 .OR. 
     &                  ICAUC.LT.0 .OR. ICAUD.LT.0       ) THEN
                      NWARN = NWARN + 1
                      WRITE(LUPRI,'(/A/2A,/A)') '@ WARNING:',
     &                 '@ NEGATIVE EXPANSION COEFFICIENTS NOT',
     &                 ' AVAILABLE FOR SECOND HYPERPOLARIZABILITIES.',
     &                 '@ INPUT LINE IGNORED...'
                    ELSE
                      NCRDISP = NCRDISP + 1
                      ICCAUA(NCRDISP) = ICAUA
                      ICCAUB(NCRDISP) = ICAUB
                      ICCAUC(NCRDISP) = ICAUC
                      ICCAUD(NCRDISP) = ICAUD
                    END IF
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)') 
     &               ' NO. OF EXPANSION COEFFICIENTS ',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXCRDISP
                    CALL QUIT('TOO MANY EXPANSION COEFFICIENTS IN CCCR')
                  END IF
                END IF
                READ (LUCMD,'(A)') LINE
              END DO
              BACKSPACE(LUCMD)
            GO TO 100


C           ------------------------------------------------
C           .AVERAG : calculate averaged tensor components
C           ------------------------------------------------
13          CONTINUE

*           first line: type of property:
              READ (LUCMD,'(A)') LINE

              IF (LINE(1:9).EQ.'GAMMA_PAR') THEN
                GAMMA_PAR = .TRUE.
              ELSE IF (LINE(1:9).EQ.'GAMMA_ISO') THEN
                GAMMA_PAR = .TRUE.
                GAMMA_ORT = .TRUE.
              END IF

              IF (GAMMA_PAR .OR. GAMMA_ORT) THEN

*               second line: symmetry:
                READ (LUCMD,'(A)') LINE
                CSYM = 'GENERI'
                IF (LINE(1:6).EQ.'ATOMIC') THEN
                  CSYM = 'ATOMIC'  ! an atom
                ELSE IF (LINE(1:6).EQ.'SPHTOP') THEN
                  CSYM = 'SPHTOP'  ! spherical top
                ELSE IF (LINE(1:6).EQ.'LINEAR') THEN
                  CSYM = 'LINEAR'  ! linear molecule
                ELSE IF (LINE(1:5).EQ.'GENER') THEN
                  CSYM = 'GENERI'  ! use generic point group symmetry
                ELSE
                  WRITE (LUPRI,*)
     &                   'WARNING: unknown symmetry input in *CCCR:'
                  WRITE (LUPRI,*) LINE
                  WRITE (LUPRI,*)'WARNING: input line ignored...'
                END IF

                IF (NCROPER.NE.0) THEN
                  NWARN = NWARN + 1
                  WRITE(LUPRI,'(/2A/A/)') 
     &             '@ WARNING: INPUT FOR .DIPOL OR .OPERATOR OPTIONS',
     &             ' BEFORE THE .AVERAG OPTION',
     &             '@ IN *CCCR SECTION WILL BE IGNORED.'
                  NCROPER = 0
                END IF

*               set operators quadruples for gamma components:
                IDIP(1) = INDPRPCC('XDIPLEN ')
                IDIP(2) = INDPRPCC('YDIPLEN ')
                IDIP(3) = INDPRPCC('ZDIPLEN ')
                DO IDX=1,3
                  IF (IDX.EQ.1) THEN
                    IDXA = 1 ! X \  XXZZ 
                    IDXB = 3 ! Z /     + permutations
                    IDXC = 3 ! Z -  ZZZZ 
                  ELSE IF (IDX.EQ.2) THEN
                    IDXA = 2 ! Y \  YYZZ
                    IDXB = 3 ! Z /     + permutations
                    IDXC = 1 ! X -  XXXX
                  ELSE IF (IDX.EQ.3) THEN
                    IDXA = 1 ! X \  XXYY
                    IDXB = 2 ! Y /     + permutations
                  IDXC = 2 ! Y -  YYYY
                  ELSE
                    CALL QUIT('Error in CC_CRINP.')
                  END IF

*                 note that the order is very important!
                  IACROP(1+(IDX-1)*7) = IDIP(IDXC) !  1.: gamma_{zzzz}
                  IBCROP(1+(IDX-1)*7) = IDIP(IDXC) !  8.: gamma_{xxxx}
                  ICCROP(1+(IDX-1)*7) = IDIP(IDXC) ! 15.: gamma_{yyyy}
                  IDCROP(1+(IDX-1)*7) = IDIP(IDXC)

                  IACROP(2+(IDX-1)*7) = IDIP(IDXB) !  2.: gamma_{zxxz}
                  IBCROP(2+(IDX-1)*7) = IDIP(IDXA) !  9.: gamma_{zyyz}
                  ICCROP(2+(IDX-1)*7) = IDIP(IDXA) ! 16.: gamma_{yxxy}
                  IDCROP(2+(IDX-1)*7) = IDIP(IDXB)

                  IACROP(3+(IDX-1)*7) = IDIP(IDXA) !  3.: gamma_{xxzz}
                  IBCROP(3+(IDX-1)*7) = IDIP(IDXA) ! 10.: gamma_{yyzz}
                  ICCROP(3+(IDX-1)*7) = IDIP(IDXB) ! 17.: gamma_{xxyy}
                  IDCROP(3+(IDX-1)*7) = IDIP(IDXB)

                  IACROP(4+(IDX-1)*7) = IDIP(IDXA) !  4.: gamma_{xzxz}
                  IBCROP(4+(IDX-1)*7) = IDIP(IDXB) ! 11.: gamma_{yzyz}
                  ICCROP(4+(IDX-1)*7) = IDIP(IDXA) ! 18.: gamma_{xyxy}
                  IDCROP(4+(IDX-1)*7) = IDIP(IDXB)

                  IACROP(5+(IDX-1)*7) = IDIP(IDXA) !  5.: gamma_{xzzx}
                  IBCROP(5+(IDX-1)*7) = IDIP(IDXB) ! 12.: gamma_{yzzy}
                  ICCROP(5+(IDX-1)*7) = IDIP(IDXB) ! 19.: gamma_{xyyx}
                  IDCROP(5+(IDX-1)*7) = IDIP(IDXA)

                  IACROP(6+(IDX-1)*7) = IDIP(IDXB) !  6.: gamma_{zzxx}
                  IBCROP(6+(IDX-1)*7) = IDIP(IDXB) ! 13.: gamma_{zzyy}
                  ICCROP(6+(IDX-1)*7) = IDIP(IDXA) ! 20.: gamma_{yyxx}
                  IDCROP(6+(IDX-1)*7) = IDIP(IDXA)

                  IACROP(7+(IDX-1)*7) = IDIP(IDXB) !  7.: gamma_{zxzx}
                  IBCROP(7+(IDX-1)*7) = IDIP(IDXA) ! 13.: gamma_{zyzy}
                  ICCROP(7+(IDX-1)*7) = IDIP(IDXB) ! 21.: gamma_{yxyx}
                  IDCROP(7+(IDX-1)*7) = IDIP(IDXA)
                END DO

                NCROPER = 21
                IF (CSYM(1:6).EQ.'ATOMIC') THEN
                  IF (GAMMA_PAR) NCROPER = 1
                  IF (GAMMA_ORT) NCROPER = 3
                ELSE IF (CSYM(1:6).EQ.'SPHTOP') THEN
                  IF (GAMMA_PAR) NCROPER = 4
                  IF (GAMMA_ORT) NCROPER = 4
                ELSE IF (CSYM(1:6).EQ.'LINEAR') THEN
                  IF (GAMMA_PAR) NCROPER = 8
                  IF (GAMMA_ORT) THEN
                    NCROPER = 10
                    IACROP(9)  = IDIP(1) ! 9.: gamma_{xyyx}
                    IBCROP(9)  = IDIP(2) 
                    ICCROP(9)  = IDIP(2)  
                    IDCROP(9)  = IDIP(1)
                    IACROP(10) = IDIP(1) ! 10.: gamma_{xxyy}
                    IBCROP(10) = IDIP(1) 
                    ICCROP(10) = IDIP(2)  
                    IDCROP(10) = IDIP(2)
                  END IF
                END IF
              END IF 
            GO TO 100


C           ----------------------------------------
C           .DISPCF : (even) dispersion coefficients
C                     for real response functions
C           ----------------------------------------
14          CONTINUE
              READ (LUCMD,*) NCRDSPE
              IF (NCRDISP.NE.0) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/2A)') 
     &            '@ WARNING: INPUT FOR .EXPCOF OPTION BEFORE .DISPCF',
     &            ' IN *CCCR SECTION WILL BE IGNORED.'
                NCRDISP = 0
              END IF
              DO L = 0, NCRDSPE, 2
                IF ((NCRDISP+(L+3)*(L+2)*(L+1)/6).LE.MXCRDISP) THEN
                  DO K = 0, L, 1
                    DO M = 0, K, 1
                      DO N = 0, M, 1
                        NCRDISP = NCRDISP + 1
                        ICCAUA(NCRDISP) = L-K
                        ICCAUB(NCRDISP) = K-M
                        ICCAUC(NCRDISP) = M-N
                        ICCAUD(NCRDISP) = N
                      END DO
                    END DO
                  END DO
                ELSE
                  NWARN = NWARN + 1
                  WRITE(LUPRI,'(/A/2A,I5/)') '@ WARNING:',
     &             '@ NO. OF DISPERSION COEFFICIENTS NEEDED',
     &             ' IS GREATER THAN THE ALLOWED NUMBER : ',MXCRDISP
                  WRITE(LUPRI,'(/A,I2,A)')
     &             '@ DISPERSION COEFFICIENTS OF ORDER',L,' ARE IGNORED'
                  WRITE(LUPRI,'(/2A,I5)') '@ FOR NEXT ORDER INCREASE ',
     &             'MXCRDISP TO:', (NCRDISP+(L+3)*(L+2)*(L+1)/6)
                  NCRDSPE = L-2
                END IF
              END DO
            GO TO 100

C           ------------------------------------------
C           .ODDISP : (odd) dispersion coefficients
C                     for imaginary response functions
C           ------------------------------------------
15          CONTINUE
              WRITE (LUPRI,*)
     &           '.ODDISP option not yet implemented in CCCR.'
            GO TO 100

C           -----------------------------------------------------------
C           .NO2NP1: switch off 2n+1/2n+2 rule for 2.-order Cauchy vec.
C           -----------------------------------------------------------
16          CONTINUE
              NO_2NP1_RULE = .TRUE.
            GO TO 100

C           -----------------------------------------------------------
C           .L2BCD : use L2(BC), L2(BD), L2(CD) vectors instead of 
C                        R2(AD), R2(AC), R2(AB) for freq.-dep. resp.
C           -----------------------------------------------------------
17          CONTINUE
              USE_LBCD = .TRUE.
              IF (L_USE_XKS3) THEN
                L_USE_XKS3 = .FALSE.
                WRITE(LUPRI,*) '.L2 BCD and .USEXKS are incompatible'
                WRITE(LUPRI,*) 'in the *CCCR section...'
                WRITE(LUPRI,*) '.USEXKS is switched off'
              END IF
              IF (L_USE_CHI2) THEN
                L_USE_CHI2 = .FALSE.
                WRITE(LUPRI,*) '.L2 BCD and .USECHI are incompatible'
                WRITE(LUPRI,*) 'in the *CCCR section...'
                WRITE(LUPRI,*) '.USECHI is switched off'
              END IF
            GO TO 100

C           -----------------------------------------------------------
C           .L2BC  : use L2(BC) instead of R2(AD) for freq.-dep. resp.
C           -----------------------------------------------------------
18          CONTINUE
              USE_L2BC = .TRUE.
              IF (L_USE_XKS3) THEN
                L_USE_XKS3 = .FALSE.
                WRITE(LUPRI,*) '.L2 BC  and .USEXKS are incompatible'
                WRITE(LUPRI,*) 'in the *CCCR section...'
                WRITE(LUPRI,*) '.USEXKS is switched off'
              END IF
              IF (L_USE_CHI2) THEN
                L_USE_CHI2 = .FALSE.
                WRITE(LUPRI,*) '.L2 BC  and .USECHI are incompatible'
                WRITE(LUPRI,*) 'in the *CCCR section...'
                WRITE(LUPRI,*) '.USECHI is switched off'
              END IF
            GO TO 100

C           -------------
C           unused labels
C           -------------
19          CONTINUE
20          CONTINUE

          ELSE
           WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF

      END IF

200   CONTINUE
*---------------------------------------------------------------------*
* check, if any quadruples of operator labels specified:
* if not, use default: complete dipole-dipole-dipole-dipole tensor
*---------------------------------------------------------------------*
      IF (NCROPER .EQ. 0) THEN
          IF (NCROPER+81 .GT. MXCROP) THEN
            WRITE(LUPRI,'(2(/A,I5))')
     &      ' NO. OF OPERATOR QUADRUPLES SPECIFIED  : ',NCROPER+81,
     &      ' IS GREATER THAN THE ALLOWED NUMBER : ',MXCROP
            CALL QUIT('TOO MANY OPERATOR QUADRUPLES IN CCQR.')
          END IF
          IDIP(1) = INDPRPCC('XDIPLEN ')
          IDIP(2) = INDPRPCC('YDIPLEN ')
          IDIP(3) = INDPRPCC('ZDIPLEN ')
          DO IDXA=1,3
          DO IDXB=1,3
          DO IDXC=1,3
          DO IDXD=1,3
            IDX = NCROPER + (IDXA-1)*27+(IDXB-1)*9+(IDXC-1)*3+IDXD
            IACROP(IDX) = IDIP(IDXA)
            IBCROP(IDX) = IDIP(IDXB)
            ICCROP(IDX) = IDIP(IDXC)
            IDCROP(IDX) = IDIP(IDXD)
          END DO
          END DO
          END DO
          END DO
          NCROPER = NCROPER + 81
      END IF

*---------------------------------------------------------------------*
* check, if frequencies or dispersion coefficients specified; 
* if not, use default: static hyperpolarizabilities
*---------------------------------------------------------------------*
      IF (NCRFREQ.EQ.0 .AND. NCRDISP.EQ.0) THEN
        NCRFREQ = NCRFREQ + 1
        BCRFR(NCRFREQ) = ZERO
        CCRFR(NCRFREQ) = ZERO
        DCRFR(NCRFREQ) = ZERO
      END IF

*---------------------------------------------------------------------*
* add list with wa frequencies:
*---------------------------------------------------------------------*
      DO IFREQ = 1, NCRFREQ
        ACRFR(IFREQ) = - (BCRFR(IFREQ) + CCRFR(IFREQ) + DCRFR(IFREQ))
      END DO

*---------------------------------------------------------------------*
* set CCCR flags:
*---------------------------------------------------------------------*
      CCCR  = .TRUE.

      RETURN
      END
*=====================================================================*
c /* deck CC_4RINP */
*=====================================================================*
       SUBROUTINE CC_4RINP(WORD)
*---------------------------------------------------------------------*
*
*    Purpose: read input for CC dynamic third hyperpolarizabilities
*             (the quartic response function)
*
*    Written by Christof Haettig, April 1997
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <cc4rinf.h>

* local parameters:
      CHARACTER MSGDBG*(18)
      PARAMETER (MSGDBG='[debug] CC_4RINP> ')
      CHARACTER SECNAM*(8)
      PARAMETER (SECNAM='CC_4RINP')

      INTEGER NTABLE
      PARAMETER (NTABLE = 10)

#if defined (SYS_CRAY)
      REAL ZERO
#else
      DOUBLE PRECISION ZERO
#endif
      PARAMETER (ZERO = 0.0d00)


* variables:
      LOGICAL SET
      SAVE SET

      CHARACTER*(7) WORD
      CHARACTER*(8) LABELA, LABELB, LABELC, LABELD, LABELE
      CHARACTER*(7) TABLE(NTABLE)

      INTEGER IDX, IJUMP, IFREQ, IDIP(3)
      INTEGER IDXA, IDXB, IDXC, IDXD, IDXE
      INTEGER MFREQ

      DATA SET /.FALSE./
      DATA TABLE /'.OPERAT','.DIPOLE','.PRINT ','.STATIC','.MIXFRE',
     &            '.4HGFRE','.USECHI','.XXXXXX','.XXXXXX','.XXXXXX'/

      INTEGER INDPRPCC

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (SET) RETURN
      SET = .TRUE.

*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*
      N4ROPER = 0
      N4RFREQ = 0

      L_USE_CHI3 = .FALSE.

      CC4R      = .FALSE.

      IPR4HYP = IPRINT

*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CC4R  ') THEN

100   CONTINUE
       
* get new input line:
        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
          CALL UPCASE(WORD)
        END DO

        IF (WORD(1:1) .EQ. '.') THEN

c         table look up:
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO
     
c         jump to the appropiate input section:
          IF (IJUMP .LE. NTABLE) THEN
            ICHANG = ICHANG + 1
            GOTO (1,2,3,4,5,6,7,8,9,10), IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_4RINP.')
      
C           ------------
C           .OPERAT
C           ------------
1           CONTINUE
              READ (LUCMD,'(5A)') LABELA,LABELB,LABELC,LABELD,LABELE
              DO WHILE (LABELA(1:1).NE.'.' .AND. LABELA(1:1).NE.'*')
                IF (LABELA(1:1).NE.'!') THEN
                  IF (N4ROPER.LT.MX4ROP) THEN
                    N4ROPER = N4ROPER + 1
                    IA4ROP(N4ROPER) = INDPRPCC(LABELA)
                    IB4ROP(N4ROPER) = INDPRPCC(LABELB)
                    IC4ROP(N4ROPER) = INDPRPCC(LABELC)
                    ID4ROP(N4ROPER) = INDPRPCC(LABELD)
                    IE4ROP(N4ROPER) = INDPRPCC(LABELE)
                  ELSE
                    WRITE(LUPRI,'(/2A,I5/)')
     &               ' NO. OF OPERATOR QUADRUPLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MX4ROP
                    CALL QUIT('TOO MANY OPERATOR QUINTUPLES IN CC4R.')
                  END IF
                END IF
                READ (LUCMD,'(5A)') LABELA,LABELB,LABELC,LABELD,LABELE
              END DO
              BACKSPACE(LUCMD)
            GO TO 100

C           -----------------------------------------------------------
C           .DIPOLE: calculate complete dipole^5 tensor (243 elements!)
C           -----------------------------------------------------------
2           CONTINUE
              IF (N4ROPER+243 .GT. MX4ROP) THEN
                WRITE(LUPRI,'(/2A,I5/)')
     &           ' NO. OF OPERATOR QUINTUPLES SPECIFIED',
     &           ' IS GREATER THAN THE ALLOWED NUMBER : ',MX4ROP
                CALL QUIT('TOO MANY OPERATOR QUINTUPLES IN CC4R.')
              END IF
              IDIP(1) = INDPRPCC('XDIPLEN ')
              IDIP(2) = INDPRPCC('YDIPLEN ')
              IDIP(3) = INDPRPCC('ZDIPLEN ')
              DO IDXA=1,3
              DO IDXB=1,3
              DO IDXC=1,3
              DO IDXD=1,3
              DO IDXE=1,3
                IDX = N4ROPER + (IDXA-1)*81 + (IDXB-1)*27 +
     &                            (IDXC-1)*9 + (IDXD-1)*3 + IDXE
                IA4ROP(IDX) = IDIP(IDXA)
                IB4ROP(IDX) = IDIP(IDXB)
                IC4ROP(IDX) = IDIP(IDXC)
                ID4ROP(IDX) = IDIP(IDXD)
                IE4ROP(IDX) = IDIP(IDXE)
              END DO
              END DO
              END DO
              END DO
              END DO
              N4ROPER = N4ROPER + 243
            GO TO 100

C           ------------
C           .PRINT
C           ------------
3           CONTINUE
              READ (LUCMD,*) IPR4HYP
            GO TO 100

C           ------------
C           .STATIC
C           ------------
4           CONTINUE
              IF (N4RFREQ+1 .GT. MX4RFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5)/)') '@ WARNING:',
     &          '@ NUMBER OF FREQUENCIES SPECIFIED    : ',N4RFREQ+1,
     &          '@ IS GREATER THAN THE ALLOWED NUMBER : ',MX4RFR,
     &          '@ INPUT OPTION STATIC WILL BE IGNORED.'
              ELSE
                N4RFREQ = N4RFREQ + 1
                B4RFR(N4RFREQ) = ZERO
                C4RFR(N4RFREQ) = ZERO
                D4RFR(N4RFREQ) = ZERO
                E4RFR(N4RFREQ) = ZERO
              END IF
            GO TO 100

C           -------------------------------------------------------
C           .MIXFRE : mixed frequency input:
C                     read wb, wc, wd, we  --->  wa = -wb-wc-wd-we
C           -------------------------------------------------------
5           CONTINUE
              READ (LUCMD,*) MFREQ
              IF (N4RFREQ+MFREQ .GT. MX4RFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5)/)') '@ WARNING:',
     &          '@ NUMBER OF FREQUENCIES SPECIFIED    : ',N4RFREQ+MFREQ,
     &          '@ IS GREATER THAN THE ALLOWED NUMBER : ',MX4RFR,
     &          '@ THE NUMBER IS RESET TO THE MAXIMUM : ',MX4RFR
                MFREQ = MX4RFR-N4RFREQ
              END IF
              READ (LUCMD,*) (B4RFR(IDX),IDX=N4RFREQ+1,N4RFREQ+MFREQ)
              READ (LUCMD,*) (C4RFR(IDX),IDX=N4RFREQ+1,N4RFREQ+MFREQ)
              READ (LUCMD,*) (D4RFR(IDX),IDX=N4RFREQ+1,N4RFREQ+MFREQ)
              READ (LUCMD,*) (E4RFR(IDX),IDX=N4RFREQ+1,N4RFREQ+MFREQ)
              N4RFREQ = N4RFREQ + MFREQ
            GO TO 100

C           ----------------------------------------------------
C           .4HGFRE : fourth harmonic generation frequencies
C                     read wb --> wc=wb, wd=wb, we=wb, wa= -4wb
C           ----------------------------------------------------
6           CONTINUE
              READ (LUCMD,*) MFREQ
              IF (N4RFREQ+MFREQ .GT. MX4RFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5)/)') '@ WARNING:',
     &          '@ NUMBER OF FREQUENCIES SPECIFIED    : ',N4RFREQ+MFREQ,
     &          '@ IS GREATER THAN THE ALLOWED NUMBER : ',MX4RFR,
     &          '@ THE NUMBER IS RESET TO THE MAXIMUM : ',MX4RFR
                MFREQ = MX4RFR-N4RFREQ
              END IF
              READ (LUCMD,*) (B4RFR(IDX),IDX=N4RFREQ+1,N4RFREQ+MFREQ)
              DO IDX = N4RFREQ+1, N4RFREQ+MFREQ
                C4RFR(IDX) = B4RFR(IDX)
                D4RFR(IDX) = B4RFR(IDX)
                E4RFR(IDX) = B4RFR(IDX)
              END DO
              N4RFREQ = N4RFREQ + MFREQ
            GO TO 100


C           -------------------------------------------------------
C           .USECHI : use second-order chi vectors as intermediates
C                     (test option)
C           -------------------------------------------------------
7           CONTINUE
              L_USE_CHI3 = .TRUE.
            GO TO 100

C           -------------
C           unused labels
C           -------------
8           CONTINUE
9           CONTINUE
10          CONTINUE
              WRITE (LUPRI,*) 'unused .XXXXXX label... ignored'
            GO TO 100


          ELSE
           WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF

      END IF

200   CONTINUE
*---------------------------------------------------------------------*
* check, if any quintuples of operator labels specified:
* if not, use default: complete dipole^5 tensor
*---------------------------------------------------------------------*
      IF (N4ROPER .EQ. 0) THEN
          IF (N4ROPER+243 .GT. MX4ROP) THEN
            WRITE(LUPRI,'(2(/A,I5))')
     &      ' NO. OF OPERATOR QUINTUPLES SPECIFIED  : ',N4ROPER+243,
     &      ' IS GREATER THAN THE ALLOWED NUMBER : ',MX4ROP
            CALL QUIT('TOO MANY OPERATOR QUINTUPLES IN CC4R.')
          END IF
          IDIP(1) = INDPRPCC('XDIPLEN ')
          IDIP(2) = INDPRPCC('YDIPLEN ')
          IDIP(3) = INDPRPCC('ZDIPLEN ')
          DO IDXA=1,3
          DO IDXB=1,3
          DO IDXC=1,3
          DO IDXD=1,3
          DO IDXE=1,3
            IDX = N4ROPER + (IDXA-1)*81 + (IDXB-1)*27 +
     &                       (IDXC-1)*9 + (IDXD-1)*3 + IDXE
            IA4ROP(IDX) = IDIP(IDXA)
            IB4ROP(IDX) = IDIP(IDXB)
            IC4ROP(IDX) = IDIP(IDXC)
            ID4ROP(IDX) = IDIP(IDXD)
            IE4ROP(IDX) = IDIP(IDXE)
          END DO
          END DO
          END DO
          END DO
          END DO
          N4ROPER = N4ROPER + 243
      END IF

*---------------------------------------------------------------------*
* check, if frequencies specified; if not, use default: static
*---------------------------------------------------------------------*
      IF (N4RFREQ .EQ. 0) THEN
        N4RFREQ = N4RFREQ + 1
        B4RFR(N4RFREQ) = ZERO
        C4RFR(N4RFREQ) = ZERO
        D4RFR(N4RFREQ) = ZERO
        E4RFR(N4RFREQ) = ZERO
      END IF

*---------------------------------------------------------------------*
* add list with wa frequencies:
*---------------------------------------------------------------------*
      DO IFREQ = 1, N4RFREQ
        A4RFR(IFREQ) = - (B4RFR(IFREQ) + C4RFR(IFREQ) 
     &                    + D4RFR(IFREQ) + E4RFR(IFREQ))
      END DO

*---------------------------------------------------------------------*
* set CC4R flags:
*---------------------------------------------------------------------*
      CC4R  = .TRUE.

      RETURN
      END
*=====================================================================*
*=====================================================================*
c /* deck CC_5RINP */
*=====================================================================*
       SUBROUTINE CC_5RINP(WORD)
*---------------------------------------------------------------------*
*
*    Purpose: read input for CC dynamic fourth hyperpolarizabilities
*             (the pentic response function)
*
*    Written by Christof Haettig, Maj 1997
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <cc5rinf.h>
#include <cc5perm.h>

* local parameters:
      CHARACTER MSGDBG*(18)
      PARAMETER (MSGDBG='[debug] CC_5RINP> ')
      CHARACTER SECNAM*(8)
      PARAMETER (SECNAM='CC_5RINP')

      INTEGER NTABLE
      PARAMETER (NTABLE = 10)

#if defined (SYS_CRAY)
      REAL ZERO
#else
      DOUBLE PRECISION ZERO
#endif
      PARAMETER (ZERO = 0.0d00)


* variables:
      LOGICAL SET
      SAVE SET

      CHARACTER*(7) WORD
      CHARACTER*(8) LABEL(6)
      CHARACTER*(7) TABLE(NTABLE)

      INTEGER IDX, IJUMP, IFREQ, IDIP(3)
      INTEGER IDXA, IDXB, IDXC, IDXD, IDXE, IDXF
      INTEGER MFREQ

      DATA SET /.FALSE./
      DATA TABLE /'.OPERAT','.DIPOLE','.PRINT ','.STATIC','.MIXFRE',
     &            '.5HGFRE','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/

      INTEGER INDPRPCC

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (SET) RETURN
      SET = .TRUE.

*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*
      N5ROPER = 0
      N5RFREQ = 0

      CC5R      = .FALSE.

      IPR5HYP = IPRINT

*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CC5R  ') THEN

100   CONTINUE
       
* get new input line:
        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
          CALL UPCASE(WORD)
        END DO

        IF (WORD(1:1) .EQ. '.') THEN

c         table look up:
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO
     
c         jump to the appropiate input section:
          IF (IJUMP .LE. NTABLE) THEN
            ICHANG = ICHANG + 1
            GOTO (1,2,3,4,5,6,7,8,9,10), IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_5RINP.')
      
C           ------------
C           .OPERAT
C           ------------
1           CONTINUE
              READ (LUCMD,'(6A)') (LABEL(IDX),IDX=1,6)
              DO WHILE (LABEL(1)(1:1).NE.'.'.AND.LABEL(1)(1:1).NE.'*')
                IF (LABEL(1)(1:1).NE.'!') THEN
                  IF (N5ROPER.LT.MX5ROP) THEN
                    N5ROPER = N5ROPER + 1
                    DO IDX = 1, 6
                      I5ROP(N5ROPER,IDX) = INDPRPCC(LABEL(IDX))
                    END DO
                    WRITE (LUPRI,*) 'CC_5RINP>',N5ROPER,LABEL
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)')
     &               ' NO. OF OPERATOR QUADRUPLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MX5ROP
                    CALL QUIT('TOO MANY OPERATOR QUINTUPLES IN CC5R.')
                  END IF
                END IF
                READ (LUCMD,'(6A)') (LABEL(IDX),IDX=1,6)
              END DO
              BACKSPACE(LUCMD)
            GO TO 100

C           -----------------------------------------------------------
C           .DIPOLE: calculate complete dipole^6 tensor (729 elements!)
C           -----------------------------------------------------------
2           CONTINUE
              IF (N5ROPER+729 .GT. MX5ROP) THEN
                WRITE(LUPRI,'(/2A,I5)')
     &           ' NO. OF OPERATOR HEXTUPLES SPECIFIED',
     &           ' IS GREATER THAN THE ALLOWED NUMBER : ',MX5ROP
                CALL QUIT('TOO MANY OPERATOR HEXTUPLES IN CC5R.')
              END IF
              IDIP(1) = INDPRPCC('XDIPLEN ')
              IDIP(2) = INDPRPCC('YDIPLEN ')
              IDIP(3) = INDPRPCC('ZDIPLEN ')
              DO IDXA=1,3
              DO IDXB=1,3
              DO IDXC=1,3
              DO IDXD=1,3
              DO IDXE=1,3
              DO IDXF=1,3
                IDX = N5ROPER + (IDXA-1)*243 + (IDXB-1)*81 +
     &               (IDXC-1)*27 + (IDXD-1)*9 + (IDXE-1)*3 + IDXF
                I5ROP(IDX,A) = IDIP(IDXA)
                I5ROP(IDX,B) = IDIP(IDXB)
                I5ROP(IDX,C) = IDIP(IDXC)
                I5ROP(IDX,D) = IDIP(IDXD)
                I5ROP(IDX,E) = IDIP(IDXE)
                I5ROP(IDX,F) = IDIP(IDXF)
C               WRITE (LUPRI,'(8i5)'), IDX, IDIP(IDXA),IDIP(IDXB),IDIP(IDXC)
C    &                              IDIP(IDXD),IDIP(IDXE),IDIP(IDXF)
              END DO
              END DO
              END DO
              END DO
              END DO
              END DO
              N5ROPER = N5ROPER + 729
            GO TO 100

C           ------------
C           .PRINT
C           ------------
3           CONTINUE
              READ (LUCMD,*) IPR5HYP
            GO TO 100

C           ------------
C           .STATIC
C           ------------
4           CONTINUE
              IF (N5RFREQ+1 .GT. MX5RFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5)/)') '@ WARNING:',
     &          '@ NUMBER OF FREQUENCIES SPECIFIED    : ',N5RFREQ+1,
     &          '@ IS GREATER THAN THE ALLOWED NUMBER : ',MX5RFR,
     &          '@ INPUT OPTION .STATIC WILL BE IGNORED.'
              ELSE
                N5RFREQ = N5RFREQ + 1
                FREQ5(N5RFREQ,B) = ZERO
                FREQ5(N5RFREQ,C) = ZERO
                FREQ5(N5RFREQ,D) = ZERO
                FREQ5(N5RFREQ,E) = ZERO
                FREQ5(N5RFREQ,F) = ZERO
              END IF
            GO TO 100

C           -------------------------------------------------------
C           .MIXFRE : mixed frequency input:
C               read wb, wc, wd, we, wf  --->  wa = -wb-wc-wd-we-wf
C           -------------------------------------------------------
5           CONTINUE
              READ (LUCMD,*) MFREQ
              IF (N5RFREQ+MFREQ .GT. MX5RFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5)/)') '@ WARNING:',
     &          '@ NUMBER OF FREQUENCIES SPECIFIED    : ',N5RFREQ+MFREQ,
     &          '@ IS GREATER THAN THE ALLOWED NUMBER : ',MX5RFR,
     &          '@ THE NUMBER IS RESET TO THE MAXIMUM : ',MX5RFR
                MFREQ = MX5RFR-N5RFREQ
              END IF
              READ (LUCMD,*) (FREQ5(IDX,B),IDX=N5RFREQ+1,N5RFREQ+MFREQ)
              READ (LUCMD,*) (FREQ5(IDX,C),IDX=N5RFREQ+1,N5RFREQ+MFREQ)
              READ (LUCMD,*) (FREQ5(IDX,D),IDX=N5RFREQ+1,N5RFREQ+MFREQ)
              READ (LUCMD,*) (FREQ5(IDX,E),IDX=N5RFREQ+1,N5RFREQ+MFREQ)
              READ (LUCMD,*) (FREQ5(IDX,F),IDX=N5RFREQ+1,N5RFREQ+MFREQ)
              N5RFREQ = N5RFREQ + MFREQ
            GO TO 100

C           ----------------------------------------------------
C           .5HGFRE : fourth harmonic generation frequencies
C                read wb --> wc=wb, wd=wb, we=wb, wf=wb, wa= -5wb
C           ----------------------------------------------------
6           CONTINUE
              READ (LUCMD,*) MFREQ
              IF (N5RFREQ+MFREQ .GT. MX5RFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5)/)') '@ WARNING:',
     &          '@ NUMBER OF FREQUENCIES SPECIFIED    : ',N5RFREQ+MFREQ,
     &          '@ IS GREATER THAN THE ALLOWED NUMBER : ',MX5RFR,
     &          '@ THE NUMBER IS RESET TO THE MAXIMUM : ',MX5RFR
                MFREQ = MX5RFR-N5RFREQ
              END IF
              READ (LUCMD,*) (FREQ5(IDX,B),IDX=N5RFREQ+1,N5RFREQ+MFREQ)
              DO IDX = N5RFREQ+1, N5RFREQ+MFREQ
                FREQ5(IDX,C) = FREQ5(IDX,B)
                FREQ5(IDX,D) = FREQ5(IDX,B)
                FREQ5(IDX,E) = FREQ5(IDX,B)
                FREQ5(IDX,F) = FREQ5(IDX,B)
              END DO
              N5RFREQ = N5RFREQ + MFREQ
            GO TO 100


C           -------------
C           unused labels
C           -------------
7           CONTINUE
8           CONTINUE
9           CONTINUE
10          CONTINUE
              WRITE (LUPRI,*) 'unused .XXXXXX label... ignored'
            GO TO 100


          ELSE
           WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF

      END IF

200   CONTINUE
*---------------------------------------------------------------------*
* check, if any quintuples of operator labels specified:
* if not, use default: complete dipole^6 tensor
*---------------------------------------------------------------------*
      IF (N5ROPER .EQ. 0) THEN
          IF (N5ROPER+729 .GT. MX5ROP) THEN
            WRITE(LUPRI,'(2(/A,I5))')
     &      ' NO. OF OPERATOR QUINTUPLES SPECIFIED  : ',N5ROPER+729,
     &      ' IS GREATER THAN THE ALLOWED NUMBER : ',MX5ROP
            CALL QUIT('TOO MANY OPERATOR QUINTUPLES IN CC5R.')
          END IF
          IDIP(1) = INDPRPCC('XDIPLEN ')
          IDIP(2) = INDPRPCC('YDIPLEN ')
          IDIP(3) = INDPRPCC('ZDIPLEN ')
          DO IDXA=1,3
          DO IDXB=1,3
          DO IDXC=1,3
          DO IDXD=1,3
          DO IDXE=1,3
          DO IDXF=1,3
              IDX = N5ROPER + (IDXA-1)*243 + (IDXB-1)*81 +
     &             (IDXC-1)*27 + (IDXD-1)*9 + (IDXE-1)*3 + IDXF
            I5ROP(IDX,A) = IDIP(IDXA)
            I5ROP(IDX,B) = IDIP(IDXB)
            I5ROP(IDX,C) = IDIP(IDXC)
            I5ROP(IDX,D) = IDIP(IDXD)
            I5ROP(IDX,E) = IDIP(IDXE)
            I5ROP(IDX,F) = IDIP(IDXF)
          END DO
          END DO
          END DO
          END DO
          END DO
          END DO
          N5ROPER = N5ROPER + 729
      END IF

*---------------------------------------------------------------------*
* check, if frequencies specified; if not, use default: static
*---------------------------------------------------------------------*
      IF (N5RFREQ .EQ. 0) THEN
        N5RFREQ = N5RFREQ + 1
        FREQ5(N5RFREQ,B) = ZERO
        FREQ5(N5RFREQ,C) = ZERO
        FREQ5(N5RFREQ,D) = ZERO
        FREQ5(N5RFREQ,E) = ZERO
        FREQ5(N5RFREQ,F) = ZERO
      END IF

*---------------------------------------------------------------------*
* add list with wa frequencies:
*---------------------------------------------------------------------*
      DO IFREQ = 1, N5RFREQ
        FREQ5(IFREQ,A) = - (FREQ5(IFREQ,B) + FREQ5(IFREQ,C) 
     &            + FREQ5(IFREQ,D) + FREQ5(IFREQ,E) + FREQ5(IFREQ,F))
      END DO

*---------------------------------------------------------------------*
* set CC5R flags:
*---------------------------------------------------------------------*
      CC5R  = .TRUE.

      RETURN
      END
*=====================================================================*
*---------------------------------------------------------------------*
C  /* Deck indprpcc */
      INTEGER FUNCTION INDPRPCC(NEWLBL)
C
#include <ccrspprp.h>
#include <priunit.h>
C
      CHARACTER*8 NEWLBL
      INTEGER I

      DO 100 I = 1,NPRLBL
         IF ( NEWLBL.EQ.PRPLBL(I) ) THEN
            INDPRPCC = I
            RETURN
         END IF
 100  CONTINUE

      NPRLBL = NPRLBL + 1

      IF (NPRLBL.GT.MAXLBL) THEN
         WRITE(LUPRI,'(A/A,I5,A,I5)')
     *   '@ NUMBER OF SPECIFIED PROPERTIES EXCEED THE MAXIMUM ALLOWED',
     *   '@ MAXPRP =',MAXLBL,' NPRLBL= ',NPRLBL
         CALL QUIT(' INDPRPCC: TOO MANY PROPERTIES SPECIFIED')
      END IF

      PRPLBL(NPRLBL) = NEWLBL
      INDPRPCC       = NPRLBL

      RETURN
      END
*---------------------------------------------------------------------*
      SUBROUTINE CC_PUT1OP(INDOP,NOP,MAXOP,OPERATOR,ROUTINE)
C
#include <priunit.h>
      CHARACTER*(*) ROUTINE
      CHARACTER*(*) OPERATOR
      CHARACTER*80 MESSAGE
      LOGICAL FAILED
      INTEGER MAXOP, NEWOP
      INTEGER INDOP(MAXOP)

      FAILED = .FALSE.
      IF (OPERATOR(1:6).EQ.'DIPLEN') THEN
        NEWOP = 3
        IF ( (NOP+NEWOP) .GT. MAXOP ) THEN
          FAILED = .TRUE.
        ELSE
          INDOP(NOP+1) = INDPRPCC('XDIPLEN ')
          INDOP(NOP+2) = INDPRPCC('YDIPLEN ')
          INDOP(NOP+3) = INDPRPCC('ZDIPLEN ')
        END IF
      ELSE IF (OPERATOR(1:6).EQ.'DIPVEL') THEN
        NEWOP = 3
        IF ( (NOP+NEWOP) .GT. MAXOP ) THEN
          FAILED = .TRUE.
        ELSE
          INDOP(NOP+1) = INDPRPCC('XDIPVEL ')
          INDOP(NOP+2) = INDPRPCC('YDIPVEL ')
          INDOP(NOP+3) = INDPRPCC('ZDIPVEL ')
        END IF
      ELSE IF (OPERATOR(1:6).EQ.'ANGMOM') THEN
        NEWOP = 3
        IF ( (NOP+NEWOP) .GT. MAXOP ) THEN
          FAILED = .TRUE.
        ELSE
          INDOP(NOP+1) = INDPRPCC('XANGMOM ')
          INDOP(NOP+2) = INDPRPCC('YANGMOM ')
          INDOP(NOP+3) = INDPRPCC('ZANGMOM ')
        END IF
      ELSE IF (OPERATOR(1:6).EQ.'SECMOM') THEN
        NEWOP = 6
        IF ( (NOP+NEWOP) .GT. MAXOP ) THEN
          FAILED = .TRUE.
        ELSE
          INDOP(NOP+1) = INDPRPCC('XXSECMOM')
          INDOP(NOP+2) = INDPRPCC('XYSECMOM')
          INDOP(NOP+3) = INDPRPCC('XZSECMOM')
          INDOP(NOP+4) = INDPRPCC('YYSECMOM')
          INDOP(NOP+5) = INDPRPCC('YZSECMOM')
          INDOP(NOP+6) = INDPRPCC('ZZSECMOM')
        END IF
      ELSE IF (OPERATOR(1:6).EQ.'ROTSTR') THEN
        NEWOP = 6
        IF ( (NOP+NEWOP) .GT. MAXOP ) THEN
          FAILED = .TRUE.
        ELSE
          INDOP(NOP+1) = INDPRPCC('XXROTSTR')
          INDOP(NOP+2) = INDPRPCC('XYROTSTR')
          INDOP(NOP+3) = INDPRPCC('XZROTSTR')
          INDOP(NOP+4) = INDPRPCC('YYROTSTR')
          INDOP(NOP+5) = INDPRPCC('YZROTSTR')
          INDOP(NOP+6) = INDPRPCC('ZZROTSTR')
        END IF
      ELSE
        CALL QUIT('Unknown OPERATOR in CC_PUT1OP')
      END IF

      IF (FAILED) THEN
        WRITE(MESSAGE,'(3a)')
     &  'TOO MANY OPERATORS IN ',ROUTINE(1:LEN(ROUTINE)),'.'
        WRITE(LUPRI,'(2(/A,I5))')
     &  ' NO. OF OPERATORS SPECIFIED         : ',NOP+NEWOP,
     &  ' IS GREATER THAN THE ALLOWED NUMBER : ',MAXOP
        CALL QUIT(MESSAGE)
      ELSE
        NOP = NOP + NEWOP
      END IF

      RETURN
      END
*---------------------------------------------------------------------*
      SUBROUTINE CC_PUT2OP(INDOP1,INDOP2,NOP,MAXOP,OPERATOR,ROUTINE)
C
#include <priunit.h>
      CHARACTER*(*) ROUTINE
      CHARACTER*(*) OPERATOR
      CHARACTER*80 MESSAGE
      LOGICAL FAILED
      INTEGER MAXOP, NEWOP
      INTEGER INDOP1(MAXOP), INDOP2(MAXOP), IOP(10)

      FAILED = .FALSE.
      IF (OPERATOR(1:6).EQ.'DIPLEN') THEN
        NEWOP = 3
        IF ( (NOP+NEWOP*NEWOP) .GT. MAXOP ) THEN
          FAILED = .TRUE.
        ELSE
          IOP(1) = INDPRPCC('XDIPLEN ')
          IOP(2) = INDPRPCC('YDIPLEN ')
          IOP(3) = INDPRPCC('ZDIPLEN ')
        END IF
      ELSE IF (OPERATOR(1:6).EQ.'DIPVEL') THEN
        NEWOP = 3
        IF ( (NOP+NEWOP*NEWOP) .GT. MAXOP ) THEN
          FAILED = .TRUE.
        ELSE
          IOP(1) = INDPRPCC('XDIPVEL ')
          IOP(2) = INDPRPCC('YDIPVEL ')
          IOP(3) = INDPRPCC('ZDIPVEL ')
        END IF
      ELSE IF (OPERATOR(1:6).EQ.'ANGMOM') THEN
        NEWOP = 3
        IF ( (NOP+NEWOP*NEWOP) .GT. MAXOP ) THEN
          FAILED = .TRUE.
        ELSE
          IOP(1) = INDPRPCC('XANGMOM ')
          IOP(2) = INDPRPCC('YANGMOM ')
          IOP(3) = INDPRPCC('ZANGMOM ')
        END IF
      ELSE
        CALL QUIT('Unknown OPERATOR in CC_PUT2OP')
      END IF

      IF (FAILED) THEN
        WRITE(MESSAGE,'(3a)')
     &  'TOO MANY OPERATORS IN ',ROUTINE(1:LEN(ROUTINE)),'.'
        WRITE(LUPRI,'(2(/A,I5))')
     &  ' NO. OF OPERATORS SPECIFIED         : ',NOP+NEWOP,
     &  ' IS GREATER THAN THE ALLOWED NUMBER : ',MAXOP
        CALL QUIT(MESSAGE)
      ELSE
        DO IDX1 = 1, NEWOP
          DO IDX2 = 1, NEWOP
            IDX12 = NOP + (IDX1 - 1)*NEWOP + IDX2
            INDOP1(IDX12) = IOP(IDX1)
            INDOP2(IDX12) = IOP(IDX2)
          END DO
        END DO
        NOP = NOP + NEWOP*NEWOP
      END IF

      RETURN
      END
*---------------------------------------------------------------------*
c /* deck cc_exlrinp */
*=====================================================================*
       SUBROUTINE CC_EXLRINP(WORD)
*---------------------------------------------------------------------*
*
* Purpose: read input for coupled cluster excited state linear response
*          calculation of frequency-dependent second-order properties
*          (excited state response functions and two-photon transition
*           moments between two excited states)
*
* Written by Christof Haettig, July 1997
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccexlrinf.h>

* local parameters:
      CHARACTER MSGDBG*(20)
      PARAMETER (MSGDBG='[debug] CC_EXLRINP> ')
      CHARACTER SECNAM*(10)
      PARAMETER (SECNAM='CC_EXLRINP')

      INTEGER NTABLE
      PARAMETER (NTABLE = 12)

#if defined (SYS_CRAY)
      REAL ZERO
#else
      DOUBLE PRECISION ZERO
#endif
      PARAMETER (ZERO = 0.0d00)


* variables:
      LOGICAL SET
      SAVE SET

      CHARACTER WORD*(7), LABHELP*(80)
      CHARACTER*8 LABELA, LABELB
      CHARACTER TABLE(NTABLE)*(7)

      INTEGER IDX, IJUMP, ISYMS(2), IDXS(2), ISTART, IEND
      INTEGER MFREQ
      INTEGER IDXA, IDXB, IDIP(3)

      DATA SET /.FALSE./

      DATA TABLE /'.OPERAT','.DIPOLE','.SELSTA','.PRINT ','.ALLSTA',
     &            '.HALFFR','.USELEF','.FREQ  ','.FREQUE','.STATIC',
     &            '.USE O2','.NOPROJ'/

      INTEGER INDPRPCC

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (SET) RETURN
      SET = .TRUE.

*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*

      NEXLROPER  = 0
      NEXLRFREQ  = 0
      NEXLRST    = 0
      ALLSTATES  = .FALSE.
      HALFFR     = .FALSE.
      USE_EL1    = .FALSE.
      USE_O2     = .FALSE.
      NOPROJ     = .FALSE.

      CCEXLR = .FALSE.

      IPREXLR = 0

      ICHANG = 0

*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCEXLR') THEN

100   CONTINUE
       
* get new input line:
        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
          CALL UPCASE(WORD)
        END DO

        IF (WORD(1:1) .EQ. '.') THEN

*         table look up:
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO
     
*         jump to the appropiate input section:
          IF (IJUMP .LE. NTABLE) THEN
            ICHANG = ICHANG + 1
            GOTO (1,2,3,4,5,6,7,8,9,10,11,12), IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_EXLRINP.')
      
C           --------------------------------------
C           .OPERAT: pair of operator lables (A,B)
C           --------------------------------------
1           CONTINUE
              READ (LUCMD,'(2A)') LABELA, LABELB
              DO WHILE (LABELA(1:1).NE.'.' .AND. LABELA(1:1).NE.'*')
                IF (LABELA(1:1).NE.'!') THEN
                  IF (NEXLROPER.LT.MXEXLROP) THEN
                    NEXLROPER = NEXLROPER + 1
                    IAEXLROP(NEXLROPER) = INDPRPCC(LABELA)
                    IBEXLROP(NEXLROPER) = INDPRPCC(LABELB)
                  ELSE
                    WRITE(LUPRI,'(/2A,I5/)') 
     &               ' NO. OF OPERATOR PAIRS SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXEXLROP
                    CALL QUIT('TOO MANY OPERATOR PAIRS IN CCEXLR.')
                  END IF
                END IF
                READ (LUCMD,'(3A)') LABELA, LABELB
              END DO
              BACKSPACE(LUCMD)
            GO TO 100

C           ------------------------------------------------
C           .DIPOLE: calculate complete dipole-dipole tensor
C           ------------------------------------------------
2           CONTINUE
              IF (NEXLROPER+9 .GT. MXEXLROP) THEN
                WRITE(LUPRI,'(2(/A,I5))') 
     &          ' NO. OF OPERATOR PAIRS SPECIFIED  : ',NEXLROPER+9,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXEXLROP
                CALL QUIT('TOO MANY OPERATOR PAIRS IN CCEXLR.')
              END IF
              IDIP(1) = INDPRPCC('XDIPLEN ')
              IDIP(2) = INDPRPCC('YDIPLEN ')
              IDIP(3) = INDPRPCC('ZDIPLEN ')
              DO IDXA=1,3
              DO IDXB=1,3
                IDX = NEXLROPER + (IDXA-1)*3+IDXB
                IAEXLROP(IDX) = IDIP(IDXA)
                IBEXLROP(IDX) = IDIP(IDXB)
              END DO
              END DO
              NEXLROPER = NEXLROPER + 9
            GO TO 100

C           ------------------------------
C           .SELSTA: select excited states
C           ------------------------------
3           CONTINUE
            READ (LUCMD,'(A80)') LABHELP
            DO WHILE(LABHELP(1:1).NE.'.' .AND. LABHELP(1:1).NE.'*')
              IF (LABHELP(1:1).NE.'!') THEN
                READ(LABHELP,*) ISYMS(1), IDXS(1),ISYMS(2), IDXS(2)
                IF (NEXLRST .LT. MXEXLRST) THEN
                  NEXLRST = NEXLRST + 1
                  IELRSYM(NEXLRST,1) = ISYMS(1)
                  IELRSTA(NEXLRST,1) = IDXS(1)
                  IELRSYM(NEXLRST,2) = ISYMS(2)
                  IELRSTA(NEXLRST,2) = IDXS(2)
                ELSE
                  NWARN = NWARN + 1
                  WRITE(LUPRI,'(/A/2A,I5)') '@ WARNING:',
     &             '@ NO. OF PAIRS OF STATES SPECIFIED',
     &             ' IS GREATER THAN THE ALLOWED NUMBER : ',MXEXLRST
                  WRITE(LUPRI,'(A,2I5/)') '@ IGNORE STATE',ISYMS,IDXS
                END IF
              END IF
              READ (LUCMD,'(A80)') LABHELP
            END DO
            BACKSPACE (LUCMD)
            GO TO 100


C           ------------
C           .PRINT
C           ------------
4           CONTINUE
              READ (LUCMD,*) IPREXLR
            GO TO 100

C           ------------------------------------------------------
C           .ALLSTA: calculate polarizabilities for all states
C           (default, if .SELSTA is not used)
C           ------------------------------------------------------
5           CONTINUE
              ALLSTATES = .TRUE.
            GO TO 100

C           --------------------------------------------------------
C           .HALFFR : use half the excitation energy as frequency
C                     for two-photon transition moments
C                     Note, that .HALFFR is incompatible with a user-
C                     specified frequency list
C                     for polarizabilities .HALFFR is equivalent
C                     to the .STATIC keyword (because the `excitation
C                     energy' is zero)
C           --------------------------------------------------------
6           CONTINUE
               HALFFR = .TRUE.
               IF (NEXLRFREQ.NE.0) THEN
                 NWARN = NWARN + 1
                 WRITE(LUPRI,'(/2a/)')
     &           '@ WARNING: in *CCEXLR on one of the Keywords',
     &           ' .HALFFR and .FREQ/FREQUE',
     &           '         can be specified...',
     &           ' .FREQ/.FREQUE input will be ignored.'
               END IF
               NEXLRFREQ  = 1
               BEXLRFR(1) = ZERO
            GO TO 100

C           -----------------------
C           .USELEF : use left excited state response vectors
C                     (default is to use right excited state responses)
C           -----------------------
7           CONTINUE
              USE_EL1 = .TRUE.
            GO TO 100


C           ------------------------------------------------
C           .FREQ  : external field frequency: wb, wa = -wb
C           .FREQUE: identical, keept for convenience
C           ------------------------------------------------
8           CONTINUE
9           CONTINUE
              READ (LUCMD,*) MFREQ
              IF (NEXLRFREQ+MFREQ .GT. MXEXLRFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5)/)') '@ WARNING:',
     &         '@ NUMBER OF FREQUENCIES SPECIFIED    :',NEXLRFREQ+MFREQ,
     &         '@ IS GREATER THAN THE ALLOWED NUMBER :',MXEXLRFR,
     &         '@ THE NUMBER IS RESET TO THE MAXIMUM :',MXEXLRFR
                MFREQ = MXEXLRFR-NEXLRFREQ
              END IF
              ISTART = NEXLRFREQ+1
              IEND   = NEXLRFREQ+MFREQ
              READ (LUCMD,*) (BEXLRFR(IDX),IDX=ISTART,IEND)
              IF (NEXLRFREQ.NE.0) WRITE (LUPRI,*)
     &           'CC_EXLRINP> ', BEXLRFR(NEXLRFREQ), NEXLRFREQ
              NEXLRFREQ = NEXLRFREQ+MFREQ
              IF (HALFFR .AND. MFREQ.GT.0) THEN
                WRITE(LUPRI,'(/2a/)')
     &          '@ WARNING: in *CCEXLR on one of the Keywords',
     &          ' .HALFFR and .FREQ/FREQUE',
     &          '         can be specified...',
     &          ' option .HALFFR will be ignored.'
              END IF
            GO TO 100

C           ---------------------------------------------------
C           .STATIC : add wb = wa = zero to frequency list
C           ---------------------------------------------------
10          CONTINUE
              IF (NEXLRFREQ+1 .GT. MXEXLRFR) THEN
                NWARN = NWARN + 1
                WRITE(LUPRI,'(/A,3(/A,I5)/)') '@ WARNING:',
     &          '@ NUMBER OF FREQUENCIES SPECIFIED    : ',NEXLRFREQ+1,
     &          '@ IS GREATER THAN THE ALLOWED NUMBER : ',MXEXLRFR,
     &          '@ INPUT OPTION .STATIC WILL BE IGNORED.'
              ELSE
                NEXLRFREQ = NEXLRFREQ + 1
                BEXLRFR(NEXLRFREQ) = ZERO
              END IF
            GO TO 100

C           -----------------------------------------------------------
C           .USE O2 : use rhs vectors for second-order amplitude 
C                     response (might save some time at the 
C                     CCS/CC2/CCSD levels if combined with other 
C                     properties, but is not (yet) implemented for CC3)
C           -----------------------------------------------------------
11          CONTINUE
              USE_O2 = .TRUE.
            GO TO 100

C           -----------------------------------------------------------
C           .NOPROJ: switch off projection onto the orthogonal 
C                    complement for the EL1/ER1 equations for
C                    excited state polarizabilities
C                    (Note that this will cause numerical problems in 
C                     the static limit)
C           -----------------------------------------------------------
12          CONTINUE
              NOPROJ = .TRUE.
            GO TO 100

          ELSE
           WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF

      END IF

200   CONTINUE

*---------------------------------------------------------------------*
* check, if any pairs of operator labels specified:
* if not, use default: complete dipole-dipole tensor
*---------------------------------------------------------------------*
      IF (NEXLROPER .EQ. 0) THEN 
        IF (NEXLROPER+9 .GT. MXEXLROP) THEN
          WRITE(LUPRI,'(2(/A,I5))') 
     &    ' NO. OF OPERATOR PAIRS SPECIFIED  : ',NEXLROPER+9,
     &    ' IS GREATER THAN THE ALLOWED NUMBER : ',MXEXLROP
          CALL QUIT('TOO MANY OPERATOR PAIRS IN CCEXLR.')
        END IF
        IDIP(1) = INDPRPCC('XDIPLEN ')
        IDIP(2) = INDPRPCC('YDIPLEN ')
        DO IDXA=1,3
        DO IDXB=1,3
          IDX = NEXLROPER + (IDXA-1)*3+IDXB
          IAEXLROP(IDX) = IDIP(IDXA)
          IBEXLROP(IDX) = IDIP(IDXB)
        END DO
        END DO
        NEXLROPER = NEXLROPER + 9
      END IF

*---------------------------------------------------------------------*
* check, if frequencies specified; if not, use the default: 
* static polarizabilities and two-photon at half the excitation energy
*---------------------------------------------------------------------*
      IF (NEXLRFREQ .EQ. 0) THEN 
        NEXLRFREQ = NEXLRFREQ + 1
        BEXLRFR(NEXLRFREQ) = ZERO
        HALFFR = .TRUE.
      END IF

*---------------------------------------------------------------------*
* check, if states specificied, if not, use default: all states
*---------------------------------------------------------------------*
      IF (NEXLRST .EQ. 0) ALLSTATES = .TRUE.

*---------------------------------------------------------------------*
* set CCEXLR flags:
*---------------------------------------------------------------------*
      CCEXLR  = .TRUE.

      RETURN
      END
*---------------------------------------------------------------------*
       SUBROUTINE CC_TMINP(WORD)
*---------------------------------------------------------------------*
*
*    Purpose: read input for CC third moment
*             three photon is a special case
*
*    if (WORD .eq '*CCTM  ') read & process input and set defaults, 
*    else set only defaults 
*
*=====================================================================*
C#if defined (IMPLICIT_NONE)
C      IMPLICIT NONE  
C#else
#  include <implicit.h>
C#endif
#include <priunit.h>
#include <cctm.h>
#include <cctminf.h>
#include <ccsdinp.h>
#include <ccsections.h>

* local parameters:
      CHARACTER SECNAM*(8)
      PARAMETER (SECNAM='CC_TMINP')

      INTEGER NTABLE
      PARAMETER (NTABLE = 10)

#if defined (SYS_CRAY)
      REAL ZERO
#else
      DOUBLE PRECISION ZERO
#endif

      PARAMETER (ZERO = 0.0d00)


* variables:
      LOGICAL SET
      SAVE SET

      CHARACTER WORD*(7)
      CHARACTER*8 LABELA, LABELB, LABELC
      CHARACTER*8 LABELD, LABELE, LABELF
      CHARACTER*70 LABHELP
      CHARACTER TABLE(NTABLE)*(7)

#if defined (SYS_CRAY)
      REAL FREQB, FREQC
#else
      DOUBLE PRECISION  FREQB, FREQC
#endif


      INTEGER IDX, IJUMP
      INTEGER  IDXA, IDXB, IDXC, IDXD, IDXE, IDXF, IDIP(3)
      INTEGER  IXSYM , IXST
      DATA SET /.FALSE./

      DATA TABLE /'.OPERAT','.DIPOLE','.PRINT ','.XXXXXX','.XXXXXX',
     &            '.SELSTA','.THIRDF','.XXXXXX','.XXXXXX','.XXXXXX'/

      INTEGER INDPRPCC

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (SET) RETURN
      SET = .TRUE.

*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*
      NTMSEL   = 0

      NTMOPER  = 0

      CCTM = .FALSE.

      IPRTM = 0

      ICHANG = 0
    
      THIRDFR = .FALSE.
 
      SELTMST  = .FALSE.
 

*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCTM  ') THEN

100   CONTINUE
       
* get new input line:
        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
          CALL UPCASE(WORD)
        END DO

        IF (WORD(1:1) .EQ. '.') THEN
C         WRITE (LUPRI,*) WORD
C         CALL FLSHFO(LUPRI)

c         table look up:
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO
     
c         jump to the appropiate input section:
          IF (IJUMP .LE. NTABLE) THEN
            ICHANG = ICHANG + 1
            GOTO (1,2,3,4,5,6,7,8,9,10), IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_TMINP.')
      
C           -------------------------------------------------
C           .OPERAT : hexuples of operator lables A,B,C,D,E,F
C           -------------------------------------------------
1           CONTINUE
              READ (LUCMD,'(6A)') LABELA, LABELB, LABELC,
     &                            LABELD, LABELE, LABELF
              DO WHILE (LABELA(1:1).NE.'.' .AND. LABELA(1:1).NE.'*')
                IF (LABELA(1:1).NE.'!') THEN
                  IF (NTMOPER.LT.MXTMOP) THEN
                    NTMOPER = NTMOPER + 1
                    IATMOP(NTMOPER) = INDPRPCC(LABELA)
                    IBTMOP(NTMOPER) = INDPRPCC(LABELB)
                    ICTMOP(NTMOPER) = INDPRPCC(LABELC)
                    IDTMOP(NTMOPER) = INDPRPCC(LABELD)
                    IETMOP(NTMOPER) = INDPRPCC(LABELE)
                    IFTMOP(NTMOPER) = INDPRPCC(LABELF)
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)') 
     &               ' NO. OF OPERATOR QUADRUPLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXTMOP
                    CALL QUIT('TOO MANY OPERATOR IN CCTM.')
                  END IF
                END IF
                READ (LUCMD,'(6A)') LABELA, LABELB, LABELC,
     &                              LABELD, LABELE, LABELF
              END DO
              BACKSPACE(LUCMD)
            GO TO 100

C           -------------------------------------------------------
C           .DIPOL : calculate full dipole-dipole-dipole
C                                  -dipole-dipole-dipole tensor
C           -------------------------------------------------------
2           CONTINUE
              IF (NTMOPER+729 .GT. MXTMOP) THEN
                WRITE(LUPRI,'(2(/A,I6))') 
     &     ' NO. OF OPERATOR QUADRUPLES SPECIFIED  : ',NTMOPER+729,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXTMOP
                CALL QUIT('TOO MANY OPERATOR QUADRUPLES IN CCTM.')

              END IF
              IDIP(1) = INDPRPCC('XDIPLEN ')
              IDIP(2) = INDPRPCC('YDIPLEN ')
              IDIP(3) = INDPRPCC('ZDIPLEN ')
              DO IDXA=1,3
              DO IDXB=1,3
              DO IDXC=1,3
              DO IDXD=1,3
              DO IDXE=1,3
              DO IDXF=1,3
                IDX = NTMOPER + (IDXA-1)*243+(IDXB-1)*81+(IDXC-1)*27+
     &                          (IDXD-1)*9  +(IDXE-1)*3 + IDXF
                IATMOP(IDX) = IDIP(IDXA)
                IBTMOP(IDX) = IDIP(IDXB)
                ICTMOP(IDX) = IDIP(IDXC)
                IDTMOP(IDX) = IDIP(IDXD)
                IETMOP(IDX) = IDIP(IDXE)
                IFTMOP(IDX) = IDIP(IDXF)
              END DO
              END DO
              END DO
              END DO
              END DO
              END DO
              NTMOPER = NTMOPER + 729
            GO TO 100

C           ------------
C           .PRINT 
C           ------------
3           CONTINUE
              READ (LUCMD,*) IPRTM
            GO TO 100

C           -----------------------
C           .XXXXXX : unused labels
C           -----------------------
4           CONTINUE
5           CONTINUE
              WRITE (LUPRI,*) 'unused .XXXXXX label... ignored'
            GO TO 100
C
C-------------------------
C           Select states.
C-------------------------
C   .SELSTAtes    Select states and frequencies
C                 frequences are overwritten if .THIRDFr are specified
C
6           CONTINUE
              SELTMST =.TRUE. 
              READ (LUCMD,'(A70)') LABHELP
              DO WHILE (LABHELP(1:1).NE.'.' .AND. LABHELP(1:1).NE.'*')
                IF (LABHELP(1:1).NE.'!') THEN
                  READ(LABHELP,*) IXSYM,IXST,FREQB,FREQC
                  IF (NTMSEL.LT.MXTMSEL) THEN
                    NTMSEL = NTMSEL + 1
                    ITMSEL(NTMSEL,1) = IXSYM
                    ITMSEL(NTMSEL,2) = IXST
                    BTMFR(NTMSEL)    = FREQB
                    CTMFR(NTMSEL)    = FREQC
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)')
     &               ' NO. OF STATES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ' ,MXTMSEL 
                    CALL QUIT('TOO MANY STATES SPECIFIED BY .SELST')
                  END IF
                END IF
                READ (LUCMD,'(A70)') LABHELP
              END DO
              BACKSPACE(LUCMD)
            GO TO 100
C
C           ------------------------------------------------
C           .THIRDF : impose condition of equal frequencies
C                      for the two lasers 
C           ------------------------------------------------
7           CONTINUE
               THIRDFR =.TRUE.
            GO TO 100

C           ------------------------------------------------
8           CONTINUE
              WRITE (LUPRI,*) 'unused .XXXXXX label... ignored'
            GO TO 100
C           ------------------------------------------------
9           CONTINUE
              WRITE (LUPRI,*) 'unused .XXXXXX label... ignored'
            GO TO 100
C           _______________________________________________
10          CONTINUE
              WRITE (LUPRI,*) 'unused .XXXXXX label... ignored'
            GO TO 100


          ELSE
           WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF

      END IF

200   CONTINUE

*---------------------------------------------------------------------*
* warning if both .SELST AND .THIRDFr is specified 
* 
      IF (SELTMST.AND.THIRDFR) THEN
         WRITE (LUPRI,*)
     &        ' WARNING: BOTH .SELST and .THIRDFr are specified'
         WRITE (LUPRI,*) ' .THIRDFr is used to obtain frequencies'
      END IF
*---------------------------------------------------------------------*
* check, if any sixtuple of operator labels specified:
* if not, use default: complete dipole tensor
*---------------------------------------------------------------------* 
      IF (NTMOPER .EQ. 0) THEN 
          IF (NTMOPER+729 .GT. MXTMOP) THEN
            WRITE(LUPRI,'(2(/A,I5))') 
     &      ' NO. OF OPERATOR SIXTUPLES SPECIFIED  : ',NTMOPER+729,
     &      ' IS GREATER THAN THE ALLOWED NUMBER : ',MXTMOP
            CALL QUIT('TOO MANY OPERATOR SIXTUPLES IN CCTM.')
          END IF
          IDIP(1) = INDPRPCC('XDIPLEN ')
          IDIP(2) = INDPRPCC('YDIPLEN ')
          IDIP(3) = INDPRPCC('ZDIPLEN ')
          DO IDXA=1,3
          DO IDXB=1,3
          DO IDXC=1,3
          DO IDXD=1,3
          DO IDXE=1,3
          DO IDXF=1,3
            IDX = NTMOPER + (IDXA-1)*243+(IDXB-1)*81+(IDXC-1)*27+
     &                      (IDXD-1)*9  +(IDXE-1)*3 + IDXF
            IATMOP(IDX) = IDIP(IDXA)
            IBTMOP(IDX) = IDIP(IDXB)
            ICTMOP(IDX) = IDIP(IDXC)
            IDTMOP(IDX) = IDIP(IDXC)
            IETMOP(IDX) = IDIP(IDXE)
            IFTMOP(IDX) = IDIP(IDXF)
          END DO
          END DO
          END DO
          END DO
          END DO
          END DO
          NTMOPER = NTMOPER + 729
      END IF

*---------------------------------------------------------------------*
* check, if frequencies are specified; if not, use default: .THIRDFR
*---------------------------------------------------------------------*
      IF ( .NOT. SELTMST ) THEN
         IF ( .NOT. THIRDFR ) THIRDFR = .TRUE.
         NINFO = NINFO + 1
         WRITE(LUPRI,'(/2A)')
     &      '@ INFO: NO FREQUENCIES SPECIFIED IN SECOND MOMENT CALC',
     &      ' DEFAULT  .THIRDFr USED '
      END IF
*---------------------------------------------------------------------*
* set CCTM flags:
*---------------------------------------------------------------------*
      WRITE (LUPRI,*) ' CCTM set to .TRUE.'
      CCTM  = .TRUE.

      RETURN
      END
*======================================================================*
       SUBROUTINE CC_MCDINP(WORD)
*----------------------------------------------------------------------*
*    Purpose: read input for CC magnetic circular dichroism
*
*    if (WORD .eq '*CCMCD  ') read & process input and set defaults, 
*    else set only defaults 
*
*    Use A,B for second order moment, C for first order moment
*
*    Sonia Coriani and Poul Joergensen (fall 1997)
*    Relaxed/PDBS operators, Sonia Coriani (february 2000)
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccmcdinf.h>
#include <ccsdinp.h>
#include <ccsections.h>

* local parameters:
      CHARACTER SECNAM*(9)
      PARAMETER (SECNAM='CC_MCDINP')
      CHARACTER*(19) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_MCDINP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER NTABLE
      PARAMETER (NTABLE = 10)

#if defined (SYS_CRAY)
      REAL ZERO
#else
      DOUBLE PRECISION ZERO
#endif

      PARAMETER (ZERO = 0.0d00)

* variables:
      LOGICAL SET
      SAVE SET

      CHARACTER WORD*(7)
      CHARACTER*8 LABELA, LABELB, LABELC
      CHARACTER*70 LABHELP
      CHARACTER TABLE(NTABLE)*(7)

      LOGICAL LARLX, LBRLX, LCRLX, LRELAX
      INTEGER IJUMP, IJ, ITOT
      INTEGER IDA(6), IDB(6), IDC(6), IDIP(3), IANG(3)
      INTEGER IXSYM , IXST
* data
      DATA SET /.FALSE./
      DATA TABLE /'.OPERAT','.MCD   ','.MCDLAO','.PRINT ','.NO2N+1',
     &            '.SELSTA','.RELAXE','.UNRELA','.USEPL1','.XXXXXX'/
      DATA IDA / 1, 2, 2, 3, 3, 1 /
      DATA IDB / 2, 1, 3, 2, 1, 3 /
      DATA IDC / 3, 3, 1, 1, 2, 2 /
* external function:
      INTEGER INDPRPCC

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (SET) RETURN
      SET = .TRUE.
*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*
      CCMCD     = .FALSE.
      NMCDST    = 0                 ! # MCD (final) states
      NMCDOPER  = 0                 ! # MCD triples
      SELMCDST  = .FALSE.           ! Select MCD fin. state (default)
      LUSE2N1   = .TRUE.            ! 2N+1 rule (Mbar^f vects in LR, default)
      LUSEPL1   = .FALSE.           ! debug use of Left transformed vectors
      IPRMCD    = 0                 ! Print level (default)

      LARLX     = .FALSE.           !Relaxed A operator
      LBRLX     = .FALSE.           !Relaxed B operator
      LCRLX     = .FALSE.           !Relaxed C operator
      LRELAX    = .FALSE.           !Relaxation 

C      RELORB1   = .FALSE.           !orbital relaxation vectors
*
      ICHANG    = 0                   
*---------------------------------------------------------------------*
*     Read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCMCD  ') THEN

100   CONTINUE
       
* get new input line:

        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
          CALL UPCASE(WORD)
        END DO

        IF (WORD(1:1) .EQ. '.') THEN
*         table look up:
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO
*         jump to the appropriate input section:
          IF (IJUMP .LE. NTABLE) THEN
            ICHANG = ICHANG + 1
            GOTO (1,2,3,4,5,6,7,8,9,10), IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_MCDINP.')
      
*           -----------------------------------------------------------
*           .OPERAT :  manually select triples of operator labels A,B,C
*                      A,B  for second order moments
*                      C    for first  order moment 
*           -----------------------------------------------------------
1           CONTINUE
              READ (LUCMD,'(3A)') LABELA, LABELB, LABELC
              DO WHILE (LABELA(1:1).NE.'.' .AND. LABELA(1:1).NE.'*')
                IF (LABELA(1:1).EQ.'!') THEN
                   CONTINUE
                ELSE IF (LABELA(1:1).EQ.'(') THEN
                   LARLX = .FALSE.
                   LBRLX = .FALSE.
                   LCRLX = .FALSE.
                   IF (LABELA(1:7).EQ.'(RELAX)') LARLX = .TRUE.
                   IF (LABELB(1:7).EQ.'(RELAX)') LBRLX = .TRUE.
                   IF (LABELC(1:7).EQ.'(RELAX)') LCRLX = .TRUE.
                   IF (LARLX .OR. LBRLX .OR. LCRLX) THEN
                      KEEPAOTWO = MAX(KEEPAOTWO,1)
C                      RELORB1   = .TRUE.
                   END IF                                    
                ELSE 
                   IF (NMCDOPER.LT.MXMCDOP) THEN
                      NMCDOPER = NMCDOPER + 1
                      IAMCDOP(NMCDOPER) = INDPRPCC(LABELA)
                      IBMCDOP(NMCDOPER) = INDPRPCC(LABELB)
                      ICMCDOP(NMCDOPER) = INDPRPCC(LABELC)
                      LAMCDRX(NMCDOPER) = LARLX
                      LBMCDRX(NMCDOPER) = LBRLX
                      LCMCDRX(NMCDOPER) = LCRLX     
                   ELSE
                      WRITE(LUPRI,'(/2A,I5)') 
     &               ' NO. OF OPERATOR TRIPLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXMCDOP
                       CALL QUIT('TOO MANY OPERATOR-TRIPLETS IN CCMCD.')
                   END IF
                END IF
                READ (LUCMD,'(3A)') LABELA, LABELB, LABELC
              END DO
              BACKSPACE(LUCMD)
            GO TO 100
*           -------------------------------------------------------
*           .MCD : calculate full tensor (r x L) * r = 6 components
*                  all operators UNRELAXED
*           -------------------------------------------------------
2           CONTINUE
              IF (NMCDOPER+6 .GT. MXMCDOP) THEN
                WRITE(LUPRI,'(2(/A,I5))') 
     &          ' NO. OF OPERATOR TRIPLES SPECIFIED  : ',NMCDOPER+6,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXMCDOP
                CALL QUIT('TOO MANY OPERATOR TRIPLES IN CCMCD.')

              END IF
              IDIP(1) = INDPRPCC('XDIPLEN ')
              IDIP(2) = INDPRPCC('YDIPLEN ')
              IDIP(3) = INDPRPCC('ZDIPLEN ')
              IANG(1) = INDPRPCC('XANGMOM ')
              IANG(2) = INDPRPCC('YANGMOM ')
              IANG(3) = INDPRPCC('ZANGMOM ')
              DO IJ = 1,6
                IAMCDOP(IJ+NMCDOPER) = IDIP(IDA(IJ))
                IBMCDOP(IJ+NMCDOPER) = IANG(IDB(IJ))
                ICMCDOP(IJ+NMCDOPER) = IDIP(IDC(IJ))
                LAMCDRX(IJ+NMCDOPER) = LRELAX
                LBMCDRX(IJ+NMCDOPER) = LRELAX
                LCMCDRX(IJ+NMCDOPER) = LRELAX       
              END DO
              NMCDOPER = NMCDOPER + 6
            GO TO 100
*           -------------------------------------------------------
*           .MCDLAO : calculate full tensor (r x L) * r = 6 compnts
*                     L operator is dh/dB
*             UNFINISHED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*           -------------------------------------------------------
3           CONTINUE
              IF (NMCDOPER+6 .GT. MXMCDOP) THEN
                WRITE(LUPRI,'(2(/A,I5))')
     &          ' NO. OF OPERATOR TRIPLES SPECIFIED  : ',NMCDOPER+6,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXMCDOP
                CALL QUIT('TOO MANY OPERATOR TRIPLES IN CCMCD.')

              END IF
              IDIP(1) = INDPRPCC('XDIPLEN ')
              IDIP(2) = INDPRPCC('YDIPLEN ')
              IDIP(3) = INDPRPCC('ZDIPLEN ')
              IANG(1) = INDPRPCC('dh/dBX  ')
              IANG(2) = INDPRPCC('dh/dBY  ')
              IANG(3) = INDPRPCC('dh/dBZ  ')
              DO IJ = 1,6
                IAMCDOP(IJ+NMCDOPER) = IDIP(IDA(IJ))
                IBMCDOP(IJ+NMCDOPER) = IANG(IDB(IJ))
                ICMCDOP(IJ+NMCDOPER) = IDIP(IDC(IJ))
                LAMCDRX(IJ+NMCDOPER) = LRELAX
                LBMCDRX(IJ+NMCDOPER) = LRELAX
                LCMCDRX(IJ+NMCDOPER) = LRELAX
              END DO
              NMCDOPER = NMCDOPER + 6
            GO TO 100
*           -------------------------------------------------------
*           .PRINT : set desired print level (default = 0)
*           -------------------------------------------------------
4           CONTINUE
              READ (LUCMD,*) IPRMCD
            GO TO 100
*           ------------------------------------------------------
*           .NO2N+1 : don't use the 2N+1 rule, ie don't use Mbar^f
*                     for the calculation of the one-photon moment 
*                     for the C operator
*           ------------------------------------------------------
5           CONTINUE
              LUSE2N1 = .FALSE.
              NWARN = NWARN + 1
              WRITE(LUPRI,'(2(/A))')
     &             '@ WARNING MCD: NO2N+1 not yet carried through',
     &             '             LUSE2N1 is reset to TRUE !!!!!'
              LUSE2N1 = .TRUE.
            GO TO 100
*           ---------------------------------------------------------------
*           .SELSTA : Select (final) states (Bfrequency zero by default)
*                     Specify then symmetry (IXSYM) and state number (IXST)
*                     of the state(s) we wish to calculate the transition 
*                     moments (one line with IXSYM,IXST for each state)
*           ---------------------------------------------------------------
6           CONTINUE

              SELMCDST =.TRUE.             
              READ (LUCMD,'(A70)') LABHELP        !read buffer line from input
              DO WHILE (LABHELP(1:1).NE.'.' .AND. LABHELP(1:1).NE.'*')
                IF (LABHELP(1:1).NE.'!') THEN
                  !read sym/number fin.state (from buffer line)
                  READ(LABHELP,*) IXSYM,IXST
                  IF (NMCDST.LT.MXMCDST) THEN
                    NMCDST = NMCDST + 1         !count how many
                    !put state-sym in array IMCDSTSY(*)
                    IMCDSTSY(NMCDST) = IXSYM
                    !put state-nr  in array IMCDSTNR(*)
                    IMCDSTNR(NMCDST) = IXST
                  ELSE
                    WRITE(LUPRI,'(/2A,I5)')
     &               ' NO. OF STATES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ' ,MXMCDST 
                    CALL QUIT(
     &                'TOO MANY STATES SPECIFIED BY .SELSTA in MCD')
                  END IF
                END IF
                READ (LUCMD,'(A70)') LABHELP
              END DO
              BACKSPACE(LUCMD)
            GO TO 100
*           ----------------------------------------------------------
*           .RELAXE : switch to relaxed modus for all three operators:
*           ----------------------------------------------------------
7           CONTINUE
              ! LRELAX    = .TRUE.
              ! KEEPAOTWO = MAX(KEEPAOTWO,1)
              WRITE (LUPRI,*) 
     &            '.RELAXE keyword in *CCMCD section is disabled.'
            GO TO 100
*           ------------------------------------------------------------
*           .UNRELA : switch to unrelaxed modus for all three operators:
*           ------------------------------------------------------------
8           CONTINUE
              LRELAX = .FALSE.
            GO TO 100                         
*           -----------------------------------------------------
*           .USEPL1 : use left transformed contributions (debug)
*           -----------------------------------------------------
9           CONTINUE
              LUSEPL1 = .TRUE.
              WRITE (LUPRI,*) SECNAM, 
     &            ': Use PL1 and left A transformations'
            GO TO 100
*           -----------------------
*           .XXXXXX : unused labels
*           -----------------------
10          CONTINUE
              WRITE (LUPRI,*) SECNAM,': unused .XXXXXX label... ignored'
            GO TO 100
*
          ELSE
            WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
            CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
            CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF

      END IF

200   CONTINUE

*---------------------------------------------------------------------*
* check if any triple of operator labels has been specified:
* if not, use default: complete unrelaxed 
* {dipole x angmom * dipole} tensor
*---------------------------------------------------------------------*
      IF (NMCDOPER .EQ. 0) THEN 
         IF (NMCDOPER+6 .GT. MXMCDOP) THEN
            WRITE(LUPRI,'(2(/A,I5))') 
     &      ' NO. OF OPERATOR TRIPLES SPECIFIED  : ',NMCDOPER+6,
     &      ' IS GREATER THAN THE ALLOWED NUMBER : ',MXMCDOP
            CALL QUIT('TOO MANY OPERATOR TRIPLES IN CCMCD.')
         END IF
         IDIP(1) = INDPRPCC('XDIPLEN ')
         IDIP(2) = INDPRPCC('YDIPLEN ')
         IDIP(3) = INDPRPCC('ZDIPLEN ')
         IANG(1) = INDPRPCC('XANGMOM ')
         IANG(2) = INDPRPCC('YANGMOM ')
         IANG(3) = INDPRPCC('ZANGMOM ')
         DO ITOT=1,6
            IAMCDOP(ITOT+NMCDOPER) = IDIP(IDA(ITOT))
            IBMCDOP(ITOT+NMCDOPER) = IANG(IDB(ITOT))
            ICMCDOP(ITOT+NMCDOPER) = IDIP(IDC(ITOT))
            LAMCDRX(ITOT+NMCDOPER) = LRELAX
            LBMCDRX(ITOT+NMCDOPER) = LRELAX
            LCMCDRX(ITOT+NMCDOPER) = LRELAX       
         END DO
         NMCDOPER = NMCDOPER + 6
      END IF
*---------------------------------------------------------------------*
* set CCMCD = TRUE if we are to calculate anything at all
*---------------------------------------------------------------------*
      CCMCD  = .TRUE.
*---------------------------------------------------------------------*
      RETURN
      END
*---------------------------------------------------------------------*
*=====================================================================*
c /* deck cc_slvinp */
*=====================================================================*
       SUBROUTINE CC_SLVINP(WORD)
C---------------------------------------------------------------------*
C
C    Purpose: read input for CC solvent calculations.
C
C    if (WORD .eq '*CCSLV ') read & process input and set defaults,
C    else set only defaults
C
C    SLV98,OC
C    Ove Christiansen April 1998
C
C=====================================================================*
#include <implicit.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <ccslvinf.h>
#include <qm3.h>

      CHARACTER SECNAM*(9)
      PARAMETER (SECNAM='CC_SLVINP')
      INTEGER NTABLE
      PARAMETER (NTABLE = 15)
 
      LOGICAL SET
      SAVE SET

      CHARACTER WORD*(7)
      CHARACTER TABLE(NTABLE)*(8)


      DATA SET /.FALSE./
      DATA TABLE /'.SOLVAT','.MXSLIT','.ETOLSL','.TTOLSL','.LTOLSL',
     *            '.PTSOLV','.CCMM','.DISCEX','.REPTST','.RELMOM',
     *            '.SLOTH '  ,'.MXINIT','.SKIPNC','.XXXXXX','.XXXXXX'/

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*

      IF (SET) RETURN
      SET = .TRUE.

*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*

      ICHANG   =  0
      IXCCSLIT =  0
      MXCCSLIT = 10
      CVGESOL  = 1.0D-07
      CVGTSOL  = 1.0D-07
      CVGLSOL  = 1.0D-07
      PTSOLV   = .FALSE.
      CCMM     = .FALSE.
      DISCEX   = .FALSE.
      ECCCU    = 0.0D0
      XTNCCCU  = 0.0D0
      XLNCCCU  = 0.0D0
      MXTINIT  = 200
      MXLINIT  = 200
      LOITER = .FALSE.
      REPTST = .FALSE.
      NREPMT = 0
      RELMOM = .FALSE.
      SLOTH = .FALSE.
      SKIPNC = .FALSE.

*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*

      IF (WORD(1:7) .EQ. '*CCSLV ') THEN


100   CONTINUE

        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
          CALL UPCASE(WORD)
        END DO

        IF (WORD(1:1) .EQ. '.') THEN

          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO

          IF (IJUMP .LE. NTABLE) THEN
            ICHANG = ICHANG + 1
            GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_SLVINP.')
C
C-----------------------
C
C-----------------------
C
1           CONTINUE
            READ(LUCMD,'(I5)') NCCSLV
            DO ISLV=1,NCCSLV
              READ(LUCMD,*) LMAXCC(ISLV),RCAVCC(ISLV),
     *                      EPSTCC(ISLV),EPOPCC(ISLV)
              IF (LMAXCC(ISLV).GT.MAXCCL) THEN
                 WRITE(LUPRI,*) 'Maximum Lmax in CC is ',MAXCCL
                 CALL QUIT('Too large LMAX in CC_SLVINP')
              ENDIF
            ENDDO
            GO TO 100
C
C-----------------------
C
C-----------------------
C
2           CONTINUE
              READ(LUCMD,*) MXCCSLIT
            GO TO 100
C
C-----------------------
C
C-----------------------
C
3           CONTINUE
               READ(LUCMD,*) CVGESOL
            GO TO 100
C
C-----------------------
C
C-----------------------
C
4           CONTINUE
               READ(LUCMD,*) CVGTSOL
            GO TO 100
C
C-----------------------
C
C-----------------------
C
5           CONTINUE
               READ(LUCMD,*) CVGLSOL
            GO TO 100
C
C-----------------------
C
C-----------------------
C
6           CONTINUE
               PTSOLV = .TRUE.
            GO TO 100
C
C-----------------------
C
C-----------------------
C
7           CONTINUE 
               CCMM = .TRUE.
            GO TO 100
C
C-----------------------
C
C-----------------------
C
8           CONTINUE
               DISCEX = .TRUE.
            GO TO 100
C
C-----------------------
C
C-----------------------
C
9           CONTINUE
               READ(LUCMD,*) NREPMT 
               REPTST = .TRUE.
            GO TO 100
C
C-----------------------
C
C-----------------------
C
10          CONTINUE
               RELMOM = .TRUE.
            GO TO 100
C
C-----------------------
C
C-----------------------
C
11          CONTINUE
              SLOTH = .TRUE.
            GO TO 100
C
C-----------------------
C
C-----------------------
C
12          CONTINUE
              READ(LUCMD,*) MXTINIT, MXLINIT
              LOITER = .TRUE.
            GO TO 100
C
C-----------------------
C
C-----------------------
C
13           CONTINUE
               SKIPNC = .TRUE.
            GO TO 100
C
C-----------------------
C
C-----------------------
C
14          CONTINUE
            GO TO 100
C
C-----------------------
C
C-----------------------
C
15           CONTINUE
            GO TO 100
C
          ELSE
           WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')

        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF

      END IF

200   CONTINUE
C
C-------------------------------------------------------------------
C     Finally if we have any solvents  put CCSLV true.
C-------------------------------------------------------------------
C
      CCSLV  = (ICHANG.GT.0)
      IF (CCSLV) RSPIM = .TRUE.
      IF (CC2) NONHF = .TRUE.
C
      RETURN
      END

c/* deck cc_r12in */
      SUBROUTINE CC_R12IN(WORD)
C     Purpose: Read input for R12 calculations.
C     Written by Wim Klopper (University of Karlsruhe, 22 November 2002).
#include <implicit.h>
#include <priunit.h>
#include <r12int.h>
CCN
#include <maxorb.h>
#include <infinp.h>
CCN
      LOGICAL SET
      CHARACTER SECNAM*(8)
      PARAMETER (SECNAM='CC_R12IN')
      PARAMETER (NTABLE = 22, D0 = 0.0D0)
      CHARACTER WORD*(7)
      CHARACTER TABLE(NTABLE)*(7)
      CHARACTER*120 CC2LAB
      DATA TABLE /'.NO HYB','.NO A  ',".NO A' ",'.NO B  ',
     &            '.NO RXR','.R12THR','.SVDTHR','.R12XXL',
     &            '.R12DIA','.R12SVD','.R12LEV','.R12RST',
     &            '.BASSCL','.NO 1  ','.NO 2  ','.R12PRP',
     &            '.CABS  ',".NO B' ",'.NO 3  ','.CC2   ',
     &            '.NATVIR','.CCVABK'/
      DATA SET/.FALSE./
      R12LEV = D0
      IF (SET) RETURN
      SET    = .TRUE.
      IF (WORD(1:4) .EQ. '*R12') THEN
  100   CONTINUE
        READ (LUCMD,'(A7)') WORD
        CALL UPCASE(WORD)
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
          CALL UPCASE(WORD)
        END DO
        IF (WORD(1:1) .EQ. '.') THEN
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO
          IF (IJUMP .LE. NTABLE) THEN
            GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &            21,22), IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_R12IN.')
    1       CONTINUE
               R12HYB = .FALSE.
            GO TO 100
    2       CONTINUE
               R12NOA = .TRUE.
            GO TO 100
    3       CONTINUE
               R12NOP = .TRUE.
            GO TO 100
    4       CONTINUE
               R12NOB = .TRUE.
            GO TO 100
    5       CONTINUE
               NORXR = .TRUE.
            GO TO 100
    6       CONTINUE
               READ (LUCMD, *) VCLTHR
            GO TO 100
    7       CONTINUE
               READ (LUCMD, *) SVDTHR
            GO TO 100
    8       CONTINUE
               R12XXL = .TRUE.
            GO TO 100
    9       CONTINUE
               R12DIA = .TRUE.
               R12SVD = .FALSE.
            GO TO 100
   10       CONTINUE
               R12SVD = .TRUE.
               R12DIA = .FALSE.
            GO TO 100
   11       CONTINUE
               READ (LUCMD, *) R12LEV
            GO TO 100
   12       CONTINUE
               R12RST = .TRUE.
            GO TO 100
   13       CONTINUE
               READ (LUCMD, *) BRASCL, KETSCL
            GO TO 100
   14       CONTINUE
               NOTONE = .TRUE.
            GO TO 100
   15       CONTINUE
               NOTTWO = .TRUE.
            GO TO 100
   16       CONTINUE
               R12PRP = .TRUE.
               IANCC2 = 1
               IF (R12NOB) IAPCC2 = 1 
               IF (R12XXL) IAPCC2 = 2
celena
               IF (R12PRP .AND. .NOT. NOTTWO) THEN
                  NOTTWO = .TRUE.
                  NWARN = NWARN + 1
                  write(lupri,'(/A/A)') '@ WARNING',
     &             '@ Sorry, calculation of R12 corrections to '//
     &             'first order properties using '//
     &             'Ansatz 2 not implemented. '//
     &             'Ansatz 2 will be ignored.' 
                  write(lupri,*)
                ENDIF
celena
            GO TO 100
   17       CONTINUE
               R12CBS = .TRUE.
            GO TO 100
   18       CONTINUE
               NOBP = .TRUE.
            GO TO 100
   19       CONTINUE
               NOTTRE = .TRUE.
            GO TO 100
   20       CONTINUE
            DO I = 1, 120
              CC2LAB(I:I) = ' '
            ENDDO
            READ (LUCMD,'(A)') CC2LAB
            DO I = 1, 120
              IF (CC2LAB(I:I) .NE. ' ') THEN
                READ (CC2LAB(I:I),'(I1)',ERR=300) IANCC2
                GOTO 201 
              END IF
            ENDDO
            GOTO 300
  201       CONTINUE
            DO I = 120, 1, -1
              IF (CC2LAB(I:I) .NE. ' ') THEN
                IF (CC2LAB(I:I) .EQ. 'A') THEN
                  IAPCC2 = 1
                  GOTO 100
                ELSE IF (CC2LAB(I:I) .EQ. 'B') THEN
                  IAPCC2 = 2
                  GOTO 100
                ELSE
                  GOTO 300
                END IF
              END IF
            END DO
            GO TO 300
   21       CONTINUE ! .NATVIR
               NATVIR = .TRUE.
               R12NOA = .TRUE.
               R12NOP = .TRUE.
               NOTTWO = .TRUE.
               NOTTRE = .TRUE.
            GO TO 100  
   22       CONTINUE ! .CCVABKL
               !use V^(alpha beta)_(kl) intermediate
               USEVABKL = .TRUE.
            GO TO 100      
          ELSE
           WRITE (LUPRI,'(/5A/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/5A/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
          CALL QUIT('Illegal prompt in '//SECNAM//'.')
        ELSE IF (WORD(1:1) .EQ.'*') THEN
          BACKSPACE (LUCMD)
          GO TO 200
        END IF
      END IF
  200 CONTINUE

      NORXR = NORXR .OR. R12HYB

      IF (IANCC2 .NE. 0) THEN
         IF (IANCC2 .EQ. 1) THEN
            NOTTWO = .TRUE.
            NOTTRE = .TRUE.
            IF (IAPCC2 .EQ. 1) THEN
               IF (NATVIR) THEN 
                  R12NOA = .FALSE.
                  R12NOP = .FALSE.
                  R12NOB = .TRUE.
                  NORXR  = .TRUE.
               END IF
               IAP = 2
            ELSE IF (IAPCC2 .EQ. 2) THEN
               IF (NATVIR) THEN 
                  R12NOA = .TRUE. 
                  R12NOB = .FALSE.
               END IF
               IF (NORXR) THEN
                  IAP = 5
               ELSE
                  IAP = 7
               END IF
            ELSE
               GOTO 300
            END IF
         ELSE IF (IANCC2 .EQ. 2) THEN
            IF (NATVIR) 
     *         CALL QUIT('Sorry, NATVIR for Ansatz 2 not implemented')
            NOTONE = .TRUE.
            NOTTRE = .TRUE.
            IF (IAPCC2 .NE. 2) GOTO 300 
            IF (NORXR) THEN
               IAP = 8
            ELSE
               IAP = 10
            END IF
         ELSE IF (IANCC2 .EQ. 3) THEN
            IF (NATVIR) 
     *         CALL QUIT('Sorry, NATVIR for Ansatz 3 not implemented')
            NOTONE = .TRUE.
            NOTTWO = .TRUE.
            IF (IAPCC2 .NE. 2) GOTO 300 
            IF (NORXR) THEN
               IAP = 8
            ELSE
               IAP = 10
            END IF
         ELSE
            GOTO 300
         END IF
         IAPCC2 = IAP
      END IF
CCN
      DIRFCK = .TRUE.
      write(lupri,*)
      write(lupri,*)'>>>> Detected R12 input:'
      write(lupri,*)'Direct Fock matrix formation switched on!'
CCN
      write(lupri,*)'Scale factors for CC excitations manifolds:'
      write(lupri,*)'for bra states (projection manifold):',brascl
      write(lupri,*)'for ket states (operator manifold)  :',ketscl

      RETURN
  300 write (lupri,*) 'WRONG CC2LAB: ',CC2LAB
      CALL QUIT('WRONG CC2LAB')
      END
C
C  /* Deck cc_chodbinp */
      SUBROUTINE CC_CHODBINP(WORD)
C
C     Thomas Bondo Pedersen, May 2002.
C
C     Purpose: Read input for CC Cholesky debug input section.
C
#include <implicit.h>
      CHARACTER*7 WORD
#include <priunit.h>
#include <chodbg.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CC_CHODBINP')

      PARAMETER (NTABLE = 4)

      LOGICAL SET
      SAVE SET

      CHARACTER*7 TABLE(NTABLE)

      DATA SET /.FALSE./
      DATA TABLE /'.DBIAJB','.STIAJB','.DBAOIN','.STAOIN'/

C     Test SET.
C     ---------

      IF (SET) RETURN
      SET = .TRUE.

C     Initializations and defaults.
C     -----------------------------

      DBAOIN = .FALSE.
      STAOIN = .FALSE.
      DBIAJB = .FALSE.
      STIAJB = .FALSE.
      ICHANG  = 0

C     Process input section.
C     ----------------------

      IF (WORD(1:7) .EQ. '*CHODBG') THEN

  100    CONTINUE

C           Read new input line.
C           --------------------

            READ(LUCMD,'(A7)') WORD
            DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
               READ (LUCMD,'(A7)') WORD
            END DO

            IF (WORD(1:1) .EQ. '.') THEN

               IJUMP = 1
               DO WHILE ((IJUMP.LE.NTABLE) .AND. (TABLE(IJUMP).NE.WORD))
                  IJUMP = IJUMP + 1
               END DO

               IF (IJUMP .LE. NTABLE) THEN

                  ICHANG = ICHANG + 1
                  GOTO (1,2,3,4), IJUMP

                  CALL QUIT
     &            ('Illegal address in computed GOTO in '//SECNAM)

    1             CONTINUE
C                    '.DBIAJB'
C                    Test Cholesky (ia|jb) integrals.
                     DBIAJB = .TRUE.
                  GOTO 100   

    2             CONTINUE
C                    '.STIAJB'
C                    Stop after (ia|jb) test.
                     STIAJB = .TRUE.
                  GOTO 100   

    3             CONTINUE
C                    '.DBAOIN'
C                    Test Cholesky AO integrals.
                     DBAOIN = .TRUE.
                  GOTO 100

    4             CONTINUE
C                    '.STAOIN'
C                    Stop after AO integral test.
                     STAOIN = .TRUE.
                  GOTO 100

               ELSE

                  WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
     &                '" not recognized in ',SECNAM,'.'
                  CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',
     &                       LUPRI)
                  CALL QUIT('Illegal Keyword in '//SECNAM)

               ENDIF

            ELSE IF (WORD(1:1) .NE. '*') THEN

               WRITE (LUPRI,'(/5A,/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
               CALL QUIT('Illegal prompt in '//SECNAM)

            ELSE IF (WORD(1:1) .EQ.'*') THEN

               BACKSPACE (LUCMD)
               GO TO 200

            ENDIF

      ENDIF

  200 CONTINUE

C     Finally, set overall Cholesky debug flag.
C     -----------------------------------------

      CHODBG = DBIAJB .OR. DBAOIN

      RETURN
      END
C  /* Deck cc_chomp2inp */
      SUBROUTINE CC_CHOMP2INP(WORD)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose: Read input for Cholesky based MP2 calculation.
C
#include <implicit.h>
      CHARACTER*7 WORD
#include <priunit.h>
#include <chomp2.h>
Casm
#include <chomp2_b.h>
Casm

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CC_CHOMP2INP')

      PARAMETER (NTABLE = 20)

      LOGICAL SET
      SAVE SET

      CHARACTER*8 TABLE(NTABLE)

      DATA SET /.FALSE./
      DATA TABLE /'.NOCHOM','.THRMP2','.SPAMP2','.MXDECM','.NCHORD',
     &            '.MP2SAV','.SKIPTR','.SKIPCH','.CHOMO ','.ALGORI',
     &            '.SPRMP2','.SCRMP2','.SPLITM','.ZERO  ','.RSTMP2',
     &            '.OLDEN2','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/

C     Test SET.
C     ---------

      IF (SET) RETURN
      SET = .TRUE.

C     Initializations and defaults.
C     Negative values of THRMP2 and SPAMP2 prompt the use of
C     corresponding AO decomposition values.
C     ------------------------------------------------------

      CALL CC_CHOMP2INIT

C     Process input section.
C     ----------------------

      IF (WORD(1:7) .EQ. '*CHOMP2') THEN

  100    CONTINUE

C           Read new input line.
C           --------------------

            READ(LUCMD,'(A7)') WORD
            DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
               READ (LUCMD,'(A7)') WORD
            END DO

            IF (WORD(1:1) .EQ. '.') THEN

               IJUMP = 1
               DO WHILE ((IJUMP.LE.NTABLE) .AND. (TABLE(IJUMP).NE.WORD))
                  IJUMP = IJUMP + 1
               END DO

               IF (IJUMP .LE. NTABLE) THEN

                  ICHANG = ICHANG + 1
                  GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     &                  16,17,18,19,20), IJUMP

                  CALL QUIT
     &            ('Illegal address in computed GOTO in '//SECNAM)

    1             CONTINUE
C                 '.NOCHOM'
C                 Do not decompose (ai|bj).
                  CHOMO = .FALSE.
                  GOTO 100

    2             CONTINUE
C                 '.THRMP2'
C                 Threshold for (ai|bj) decomposition.
                  READ(LUCMD,*) THRMP2
                  GOTO 100

    3             CONTINUE
C                 '.SPAMP2'
C                 Span factor for (ai|bj) decomposition.
                  READ(LUCMD,*) SPAMP2
                  GOTO 100

    4             CONTINUE
C                 '.MXDECM'
C                 Max. qualified diagonals in (ai|bj) decomposition.
                  READ(LUCMD,*) MXDECM
COLD              IF (MXDECM .GT. MAXMOD) THEN
COLD                 WRITE(LUPRI,'(//,5X,A,A,I10)')
COLD &               SECNAM,': MXDECM too large. Max. allowed: ',MAXMOD
COLD                 CALL QUIT('Input error in '//SECNAM)
COLD              ELSE IF (MXDECM .LE. 0) THEN
                  IF (MXDECM .LE. 0) THEN
                     WRITE(LUPRI,'(5X,A,A,I10)')
     &               SECNAM,': MXDECM must be positive!'
                     CALL QUIT('Input error in '//SECNAM)
                  ENDIF
                  GOTO 100

    5             CONTINUE
C                 '.NCHORD'
C                 Max. prev. vectors in (ai|bj) decomposition.
                  READ(LUCMD,*) NCHORD
                  IF (NCHORD .LE. 0) THEN
                     WRITE(LUPRI,'(5X,A,A,I10)')
     &               SECNAM,': NCHORD must be positive!'
                     CALL QUIT('Input error in '//SECNAM)
                  ENDIF
                  GOTO 100

    6             CONTINUE
C                 '.MP2SAV'
C                 Save MP2 amplitudes on disk.
                  MP2SAV = .TRUE.
                  GOTO 100

    7             CONTINUE
C                 '.SKIPTR'
C                 Skip MO transformation; use old vectors.
                  SKIPTR = .TRUE.
                  GOTO 100

    8             CONTINUE
C                 '.SKIPCH'
C                 Skip (ai|bj) decompositon; read info from disk
                  SKIPCH = .TRUE.
                  GOTO 100

    9             CONTINUE
C                 '.CHOMO '
C                 (ai|bj) decompositon
                  CHOMO = .TRUE.
                  GOTO 100

   10             CONTINUE
C                 '.ALGORI'
C                 algorithm:
C                 <1: decided by MP2 routine (default, IALMP2=0)
C                  1: force storage of full-square (ia|jb) in core
C                  2: batch over one virtual index
C                  3: batch over two virtual indices
C                 >3: same as <1.
                  READ(LUCMD,*) IALMP2
                  GOTO 100

   11             CONTINUE
C                 '.SPRMP2'
C                 Use sparse representation of Cholesky vector.
                  SPRMP2 = .TRUE.
                  GOTO 100

   12             CONTINUE
C                 '.SCRMP2'
C                 Screening threshold for sparse representation.
                  READ(LUCMD,*) SCRMP2
                  GOTO 100

   13             CONTINUE
C                 '.SPLITM'
C                 Weight factor for Cholesky part in memory split for
C                 virtual batch algorithms.
                  READ(LUCMD,*) TMP
                  IF (TMP .GT. 0.0D0) SPLITM = TMP
                  GOTO 100

   14             CONTINUE
C                 '.ZERO  '
C                 Threshold for diagonal zeroing in decompositions.
                  READ(LUCMD,*) THZMP2
                  GOTO 100

   15             CONTINUE
C                 '.RSTMP2'
C                 Restart MP2
                  RSTMP2 = .TRUE.
                  READ(LUCMD,*)  IFSYMB, IFVIRB
                  GOTO 100

   16             CONTINUE
C                 '.OLDEN2'
C                 Old MP2 energy (if desided)
                  OLDKNO = .TRUE.
                  READ(LUCMD,*)  OLDEN2
                  GOTO 100

   17             CONTINUE
   18             CONTINUE
   19             CONTINUE
   20             CONTINUE
C                 '.XXXXXX'
C                 Not used
                  GOTO 100

               ELSE

                  WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
     &                '" not recognized in ',SECNAM,'.'
                  CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',
     &                       LUPRI)
                  CALL QUIT('Illegal Keyword in '//SECNAM)

               ENDIF

            ELSE IF (WORD(1:1) .NE. '*') THEN

               WRITE (LUPRI,'(/5A,/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
               CALL QUIT('Illegal prompt in '//SECNAM)

            ELSE IF (WORD(1:1) .EQ.'*') THEN

               BACKSPACE (LUCMD)
               GO TO 200

            ENDIF

      ENDIF

  200 CONTINUE


      RETURN
      END
C  /* Deck cc_chocc2inp */
      SUBROUTINE CC_CHOCC2INP(WORD)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Read input for Cholesky based CC2 calculation.
C
#include <implicit.h>
      CHARACTER*7 WORD
#include <priunit.h>
#include <chocc2.h>

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CC_CHOCC2INP')

      PARAMETER (NTABLE = 15)

      LOGICAL SET
      SAVE SET

      CHARACTER*8 TABLE(NTABLE)

      DATA SET /.FALSE./
      DATA TABLE /'.CHOMO ','.THRCC2','.SPACC2','.MXDECM','.NCHORD',
     &            '.XXXXXX','.CHOT2 ','.NOCHOM','.ALGORI','.THRCCC',
     &            '.SPACCC','.SPLITM','.ZERO  ','.XXXXXX','.XXXXXX'/

C     Test SET.
C     ---------

      IF (SET) RETURN
      SET = .TRUE.

C     Initializations and defaults.
C     Negative values of THRCC2 and SPACC2 prompt the use of
C     corresponding AO decomposition values.
C     ------------------------------------------------------

      CALL CC_CHOCC2INIT

C     Process input section.
C     ----------------------

      IF (WORD(1:7) .EQ. '*CHOCC2') THEN

  100    CONTINUE

C           Read new input line.
C           --------------------

            READ(LUCMD,'(A7)') WORD
            DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
               READ (LUCMD,'(A7)') WORD
            END DO

            IF (WORD(1:1) .EQ. '.') THEN

               IJUMP = 1
               DO WHILE ((IJUMP.LE.NTABLE) .AND. (TABLE(IJUMP).NE.WORD))
                  IJUMP = IJUMP + 1
               END DO

               IF (IJUMP .LE. NTABLE) THEN

                  ICHANG = ICHANG + 1
                  GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), IJUMP

                  CALL QUIT
     &            ('Illegal address in computed GOTO in '//SECNAM)

    1             CONTINUE
C                 '.CHOMO '
C                 Decompose (ai|bj) (=> CHOT2 = .FALSE.).
                  CHOT2  = .FALSE.
                  CHOMO2 = .TRUE.
                  GOTO 100

    2             CONTINUE
C                 '.THRCC2'
C                 Threshold for decomposition.
                  READ(LUCMD,*) THRCC2
                  GOTO 100

    3             CONTINUE
C                 '.SPACC2'
C                 Span factor for decomposition.
                  READ(LUCMD,*) SPACC2
                  GOTO 100

    4             CONTINUE
C                 '.MXDECM'
C                 Max. qualified diagonals in decomposition.
                  READ(LUCMD,*) MXDEC2
                  IF (MXDEC2 .LE. 0) THEN
                     WRITE(LUPRI,'(5X,A,A,I10)')
     &               SECNAM,': .MXDECM input must be positive!'
                     CALL QUIT('Input error in '//SECNAM)
                  ENDIF
                  GOTO 100

    5             CONTINUE
C                 '.NCHORD'
C                 Max. prev. vectors in decomposition.
                  READ(LUCMD,*) NCHRD2
                  IF (NCHRD2 .LE. 0) THEN
                     WRITE(LUPRI,'(5X,A,A,I10)')
     &               SECNAM,': .NCHORD input must be positive!'
                     CALL QUIT('Input error in '//SECNAM)
                  ENDIF
                  GOTO 100

    6             CONTINUE
C                 '.XXXXXX'
                  GOTO 100

    7             CONTINUE
C                 '.CHOT2 '
C                 Decompose CC2 T2 amplitudes. (=> CHOMO2 = .FALSE.)
                  CHOT2  = .TRUE.
                  CHOMO2 = .FALSE.
                  GOTO 100

    8             CONTINUE
C                 '.NOCHOM'
C                 No decompositions in CC2 section.
                  CHOT2  = .FALSE.
                  CHOMO2 = .FALSE.
                  GOTO 100

    9             CONTINUE
C                 '.ALGORI'
C                 Set algorithm (=1 for single virtual batch, =2 for double)
                  READ(LUCMD,*) IALGO
                  IF (IALGO .LE. 1) THEN
                     IALCC2 = 1
                  ELSE
                     IALCC2 = 2
                  ENDIF
                  GOTO 100

   10             CONTINUE
C                 '.THRCCC'
C                 Threshold to use in amplitude decomposition for
C                 response intermediates and right-hand sides.
                  READ(LUCMD,*) THRCCC
                  GOTO 100

   11             CONTINUE
C                 '.SPACCC'
C                 Span factor to use in amplitude decomposition for
C                 response intermediates and right-hand sides.
                  READ(LUCMD,*) SPACCC
                  GOTO 100

   12             CONTINUE
C                 '.SPLITM'
C                 Weight factor for Cholesky part in memory split for
C                 batching over virtuals.
                  READ(LUCMD,*) TMP
                  IF (TMP .GT. 0.0D0) SPLITC = TMP
                  GOTO 100

   13             CONTINUE
C                 '.ZERO  '
C                 Threshold for diagonal zeroing in decompositions.
                  READ(LUCMD,*) THZCC2
                  GOTO 100

   14             CONTINUE
C                 '.XXXXXX'
C                 Not used.
                  GOTO 100

   15             CONTINUE
C                 '.XXXXXX'
C                 Not used.
                  GOTO 100

               ELSE

                  WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
     &                '" not recognized in ',SECNAM,'.'
                  CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',
     &                       LUPRI)
                  CALL QUIT('Illegal Keyword in '//SECNAM)

               ENDIF

            ELSE IF (WORD(1:1) .NE. '*') THEN

               WRITE (LUPRI,'(/5A,/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
               CALL QUIT('Illegal prompt in '//SECNAM)

            ELSE IF (WORD(1:1) .EQ.'*') THEN

               BACKSPACE (LUCMD)
               GO TO 200

            ENDIF

      ENDIF

  200 CONTINUE

      RETURN
      END
C  /* Deck cc_chomp2init */
      SUBROUTINE CC_CHOMP2INIT
C
C     Thomas Bondo Pedersen, October 2002.
C
C     Initialize chomp2.h
#include <implicit.h>
#include <chomp2.h>

      MP2SAV = .FALSE.
      CHOMO  = .FALSE.
      SKIPTR = .FALSE.
      SKIPCH = .FALSE.
      THRMP2 = -1.0D8
      SPAMP2 = -1.0D8
      MXDECM = 50
      NCHORD = 200
      THZMP2 = -1.0D8
      IALMP2 = 0
      SPRMP2 = .FALSE.
      SCRMP2 = -1.0D8
      SPLITM = 1.0D0

      RETURN
      END
C  /* Deck cc_chocc2init */
      SUBROUTINE CC_CHOCC2INIT
C
C     Thomas Bondo Pedersen, October 2002.
C
C     Initialize chocc2.h
#include <implicit.h>
#include <chocc2.h>
Casm
#include <chomp2_b.h>
C
      LOGICAL SET
      SAVE SET
      DATA SET /.FALSE./

C
      IF (SET) THEN
         RETURN
      ELSE
         SET = .TRUE.
      END IF
Casm
      IALCC2 = 2
      CHOT2  = .FALSE.
      CHOMO2 = .FALSE.
      CHOT2C = .FALSE.
      DSKETA = .FALSE.
      DSKFY2 = .FALSE.
      THRCC2 = -1.0D8
      SPACC2 = -1.0D8
      MXDEC2 = 50
      NCHRD2 = 200
      THZCC2 = -1.0D8
      THRCCC = THRCC2
      SPACCC = SPACC2
      MXDECC = MXDEC2
      NCHRDC = NCHRD2
      SPLITC = 1.0D0

Casm
      RSTMP2 = .FALSE.
      OLDKNO = .FALSE.
      OLDEN2 = 0.0D0
Casm
      RETURN
      END
C  /* Deck cc_choptinit */
      SUBROUTINE CC_CHOPTINIT
C
C     TBP, JLC, BFR, AS, and HK,  May 2003.
C
C     Purpose: Set defaults for Cholesky denominator CCSD(T) program.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0)
#include <cc_cho.h>

      MXCHVE = MIN(MAXCHO,10)
      THRCHO = -1.0D10
Casm
C
C     Virtual part
C
      RSTVIR = .FALSE.
      IFVISY = 1
      IFVIOR = 1
C
C     Occupied part
C
      RSTH  = .FALSE.
      RSTH1 = .FALSE.
      RSTF1 = .FALSE.
      RSTC1 = .FALSE.
      RSTC2 = .FALSE.
C
C     Files
C
      SKIVI1 = .FALSE.
      SKIVI2 = .FALSE.
C
C     Previous values
C
      UKNE4V = .TRUE.
      UKNE5V = .TRUE.
      UKNE4O = .TRUE.
      UKNE5O = .TRUE.
C
      OLD4V = ZERO
      OLD5V = ZERO
      OLD4O = ZERO
      OLD5O = ZERO
Casm
      RETURN
      END
C  /* Deck cc_choptinp */
      SUBROUTINE CC_CHOPTINP(WORD)
C
C     TBP, JLC, BFR, AS, and HK,  May 2003.
C
C     Purpose: Process input for changing defaults for the Cholesky denominator
C              CCSD(T) program.
C
#include <implicit.h>
      CHARACTER*7 WORD
#include <priunit.h>
#include <cc_cho.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CC_CHOPTINP')

      PARAMETER (NTABLE = 15)

      LOGICAL SET
      SAVE SET

      CHARACTER*8 TABLE(NTABLE)

      DATA SET /.FALSE./
      DATA TABLE /'.MXCHVE','.THRCHO','.RSTVIR','.RSTH  ','.RSTH1 ',
     &            '.RSTF1 ','.RSTC1 ','.RSTC2 ','.SKIVI1','.SKIVI2',
     &            '.OLD4V ','.OLD5V ','.OLD4O ','.OLD5O ','.XXXXXX'/

C     Test SET.
C     ---------

      IF (SET) RETURN
      SET = .TRUE.

C     Set defaults.
C     -------------

      CALL CC_CHOPTINIT

C     Process input section.
C     ----------------------

      IF (WORD(1:7) .EQ. '*CHO(T)') THEN

  100    CONTINUE

C           Read new input line.
C           --------------------

            READ(LUCMD,'(A7)') WORD
            DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
               READ (LUCMD,'(A7)') WORD
            END DO

            IF (WORD(1:1) .EQ. '.') THEN

               IJUMP = 1
               DO WHILE ((IJUMP.LE.NTABLE) .AND. (TABLE(IJUMP).NE.WORD))
                  IJUMP = IJUMP + 1
               END DO

               IF (IJUMP .LE. NTABLE) THEN

                  ICHANG = ICHANG + 1
                  GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), IJUMP

                  CALL QUIT
     &            ('Illegal address in computed GOTO in '//SECNAM)

    1             CONTINUE
C                 '.MXCHVE'
C                 Maximum number of Cholesky vectors used.
                  NCHDEF = MXCHVE
                  READ(LUCMD,*) MXCHVE
                  IF (MXCHVE .LE. 0) THEN
                     WRITE(LUPRI,'(/,1X,A,A,A,I10,A,/,1X,A,I10,/)')
     &               SECNAM,': NOTICE: Number of Cholesky vectors ',
     &               'specified (',MXCHVE,')',
     &               'is reset to the default value ',NCHDEF
                     MXCHVE = NCHDEF
                  ELSE IF (MXCHVE .GT. MAXCHO) THEN
                     WRITE(LUPRI,'(/,1X,A,A,A,I10,A,/,1X,A,I10,/)')
     &               SECNAM,': NOTICE: Number of Cholesky vectors ',
     &               'specified (',MXCHVE,')',
     &               'is reset to the maximum value ',MAXCHO
                     MXCHVE = MAXCHO
                  ENDIF
                  GOTO 100

    2             CONTINUE
C                 '.THRCHO'
C                 Threshold for skipping remaining Cholesky vectors in each
C                 term.
                  READ(LUCMD,*) THRCHO
                  GO TO 100

    3             CONTINUE
                     RSTVIR = .TRUE.
                     READ(LUCMD,*) IFVISY,IFVIOR
                  GO TO 100

    4             CONTINUE
                     RSTH = .TRUE.
                  GO TO 100
C
    5             CONTINUE
                     RSTH1 = .TRUE.
                  GO TO 100
C
    6             CONTINUE
                     RSTF1 = .TRUE.
                  GO TO 100
C
    7             CONTINUE
                     RSTC1 = .TRUE.
                  GO TO 100
C
    8             CONTINUE
                     RSTC2 = .TRUE.
                  GO TO 100
C
    9             CONTINUE
                     SKIVI1 = .TRUE.
                  GO TO 100
C
   10             CONTINUE
                     SKIVI2 = .TRUE.
                  GO TO 100

   11             CONTINUE
                     UKNE4V = .FALSE.
                     READ(LUCMD,*) OLD4V
                  GO TO 100

   12             CONTINUE
                     UKNE5V = .FALSE.
                     READ(LUCMD,*) OLD5V
                  GO TO 100

   13             CONTINUE
                     UKNE4O = .FALSE.
                     READ(LUCMD,*) OLD4O
                  GO TO 100

   14             CONTINUE
                     UKNE5O = .FALSE.
                     READ(LUCMD,*) OLD5O
                  GO TO 100

   15             CONTINUE
C                    Not used
                  GO TO 100

               ELSE

                  WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
     &                '" not recognized in ',SECNAM,'.'
                  CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',
     &                       LUPRI)
                  CALL QUIT('Illegal Keyword in '//SECNAM)

               ENDIF

            ELSE IF (WORD(1:1) .NE. '*') THEN

               WRITE (LUPRI,'(/5A,/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
               CALL QUIT('Illegal prompt in '//SECNAM)

            ELSE IF (WORD(1:1) .EQ.'*') THEN

               BACKSPACE (LUCMD)
               GO TO 200

            ENDIF

      ENDIF

  200 CONTINUE

      RETURN
C
C
      END
C  /* Deck cc_ctomag */
      SUBROUTINE CC_CTOMAG
C
C     asm September 2005
C
C     Purpose: Set up operator list for CTOCD calculations
C
C
#include <implicit.h>
#include <priunit.h>
#include <cclrinf.h>
#include <ctocdcc.h>
C
      PARAMETER (MAXOPR = 10 * MXLROP)
C
      CHARACTER*8 RECORD(4), STARS, LABEL, LABELA, LABELB
      PARAMETER (STARS = '********')
C
      CHARACTER*8 LSTLBL(MAXOPR)
      INTEGER     SYMLBL(MAXOPR)
C
      LOGICAL SET,LF
      SAVE  SET
      DATA  SET /.FALSE./
      DATA  LF /.FALSE./
C
C
      IF (SET) RETURN
      SET = .TRUE.
C
      LUPROP = -1
      CALL GPOPEN(LUPROP,'AOPROPER','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
      REWIND(LUPROP)
C
C     Read labels in AOPROPER and sort
C
      NOPER = 0
  100 CONTINUE
         READ(LUPROP, END=200, ERR=300) RECORD
         IF (RECORD(1) .NE. STARS) THEN
            GOTO 100
         ELSE
            LABEL = RECORD(4)
            IF ((LABEL .EQ. 'HUCKOVLP') .OR. (LABEL .EQ. 'HUCKEL  ')
     &          .OR. (LABEL .EQ. 'OVERLAP ')) GOTO 100
C
            NOPER = NOPER + 1
            IF (NOPER .GT. MAXOPR) STOP 
     &                             'Too many label found by CC_CTOMAG'
C
            READ(RECORD(2),'(I1)') SYMLBL(NOPER)
            LSTLBL(NOPER) = LABEL
         END IF
      GOTO  100
C
  200 CONTINUE
C
C     Select pairs of operators to compute
C
      DO I = 1,NOPER
C
         LABELA = LSTLBL(I)
         ISYMA  = SYMLBL(I)
C
         IF (LABELA(2:7) .EQ. 'DIPVEL') THEN     ! Most of ctocd properties
C
            DO J = 1,NOPER
C
               LABELB = LSTLBL(J)
               ISYMB  = SYMLBL(J)
C
               IF (ISYMA .EQ. ISYMB) THEN        ! Otherwise, sop is zero
C
                  IF (LABELB(3:6) .EQ. 'RANG') THEN          !Dia suscep
                     IF (CTOSUS)
     &                  CALL CC_LRINPREQ(LABELA,LABELB,1,1,LF,LF)
                  ELSE IF (LABELB(4:7) .EQ. 'RPSO') THEN     !Dia shield
                     IF (CTOSHI)
     &                  CALL CC_LRINPREQ(LABELA,LABELB,1,1,LF,LF)
                  ELSE IF (LABELB(1:3) .EQ. 'PSO') THEN      !Shift shield
                     IF (CTOSHI)
     &                  CALL CC_LRINPREQ(LABELA,LABELB,1,1,LF,LF)
                  END IF
C
               END IF
C
            END DO
C
         ELSE IF (CTOSUS .AND. (LABELA(2:7) .EQ. 'ANGMOM')) THEN   !Para suscep
C
            DO J = 1,NOPER
C
               LABELB = LSTLBL(J)
               ISYMB  = SYMLBL(J)
C
               IF (ISYMA .EQ. ISYMB) THEN        ! Otherwise, sop is zero
C
                  IF (LABELB(2:7) .EQ. 'ANGMOM') THEN    
                     CALL CC_LRINPREQ(LABELA,LABELB,1,1,LF,LF)
                  END IF
C
               END IF
C
            END DO
C
         ELSE IF (CTOSHI .AND. (LABELA(1:3) .EQ. 'PSO')) THEN      !Para shield
C
            DO J = 1,NOPER
C
               LABELB = LSTLBL(J)
               ISYMB  = SYMLBL(J)
C
               IF (ISYMA .EQ. ISYMB) THEN        ! Otherwise, sop is zero
C
                  IF (LABELB(2:7) .EQ. 'ANGMOM') THEN
                     CALL CC_LRINPREQ(LABELA,LABELB,1,1,LF,LF)
                  END IF
C
               END IF
C
            END DO
C
         END IF
C
      END DO
      GOTO 999
C
  300 STOP 'Error when reading in CC_CTOMAG'
C
  999 CONTINUE
C
      RETURN
      END
C  /* Deck cc_actinp */
      SUBROUTINE CC_ACTINP(WORD,MSYM)
C
C     Alfredo Sanchez de Meras. May 2008
C
C     Purpose: Read input for CC Active section
C
#include <implicit.h>
      CHARACTER*7 WORD
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <center.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_ACTINP')
C
      PARAMETER (NTABLE = 15)
      CHARACTER*7 TABLE(NTABLE)
      DATA TABLE /'.ATOMIC','.BOXDEF','.ACTFRE','.DIFADD','.NBOEXP',
     &            '.SELDIR','.THACOC','.THACVI','.SELDIA','.ORDER ',
     &            '.FULDEC','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
C
      LOGICAL SET, CHKACT
      SAVE SET
      DATA SET /.FALSE./
      DATA CHKACT /.FALSE./
C
C
      IF (SET) RETURN
      SET = .TRUE.
C
C     Initializations
C
      ATOMIC = .FALSE.
      ACTFRE = .FALSE.
      DIFADD = .FALSE.
      NBOEXP = .FALSE.
      SELDIR = .FALSE.
      ACTIVE = .TRUE.
C
      THACOC = 1.0D-4
      THACVI = 1.0D-4
C
ctmp  DO I = 1,NUCIND
ctmp     IORDEC(I) = I
ctmp  END DO
C
      ICHANG  = 0
C
C     Process input section.
C
      IF (WORD(1:7) .EQ. '*ACTIVE') THEN
C
  100    CONTINUE
C
C           Read new input line.
C
            READ(LUCMD,'(A7)') WORD
            DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
               READ (LUCMD,'(A7)') WORD
            END DO
C
            IF (WORD(1:1) .EQ. '.') THEN
C
               IJUMP = 1
               DO WHILE ((IJUMP.LE.NTABLE) .AND. (TABLE(IJUMP).NE.WORD))
                  IJUMP = IJUMP + 1
               END DO
C
               IF (IJUMP .LE. NTABLE) THEN
C
                  ICHANG = ICHANG + 1
                  GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), IJUMP
C
                  CALL QUIT
     &            ('Illegal address in computed GOTO in '//SECNAM)
C
    1             CONTINUE
C                    '.ATOMIC'
C                    Select active atoms (no boxes)
                     ATOMIC = .TRUE.
                     CHKACT = .TRUE.
                     READ(LUCMD,*) NACINP
                     IF (NACINP .GT. MXACAT) THEN
                        WRITE(LUPRI,*) 'ERROR: From center.h, ' ,
     &                                 'maximum number of active ',
     &                                 'atoms is',MXACAT
                        STOP 'Too many active atoms'
                     END IF
                     READ(LUCMD,*) (LACINP(I), I=1,NACINP)
C
                     DO I = 1,NACINP
                        DO J = 1,I-1
                           IF (LACINP(J) .EQ. LACINP(I)) THEN
                              WRITE(LUPRI,*) 'ERROR : ',
     &                                       'One atom declared ',
     &                                       'twice as active'
                              STOP 'One atom is hyperactive'
                           ELSE IF (LACINP(J) .GT. LACINP(I)) THEN
                              ITMP = LACINP(J)          ! Not needed,
                              LACINP(J) = LACINP(I)     ! but output
                              LACINP(I) = ITMP          ! looks nicer
                           END IF
                        END DO
                     END DO
                  GOTO 100   
C
    2             CONTINUE
C                    '.BOXDEF'
C                    Boxes definition.           
                     CHKACT = .TRUE.
C_to_do              copy & paste from mkinp.f
                  GOTO 100   
C
    3             CONTINUE
C                    '.ACTFRE'
C                    Freeze orbitals in active atomic space
                     ACTFRE = .TRUE.
                     READ(LUCMD,*) NACTFR
                  GOTO 100
C
    4             CONTINUE
C                    '.DIFADD'
C                    Include in active space selected (diffuse) basis
                     DIFADD = .TRUE.
                     CHKACT = .TRUE.
                     DO ISYM = 1,MSYM
                        READ(LUCMD,*) NEXTBS(ISYM)
                        READ(LUCMD,*) (IEXTBS(I,ISYM), I=1,NEXTBS(ISYM))
                     END DO
                  GOTO 100
C
    5             CONTINUE
C                    '.NBOEXP'
C                    n-body interactions among boxes
                     NBOEXP = .TRUE.
C_to_do              copy & paste from mkinp.f
                  GOTO 100
    6             CONTINUE
C                    '.SELDIR'
C                    Decompose only selected centers
                     SELDIR = .TRUE.
                  GOTO 100
    7             CONTINUE
C                    '.THACOC'
C                    Threshold for decomposition of active occupied block
                     READ(LUCMD,*) THACOC
                  GOTO 100
    8             CONTINUE
C                    '.THACVI'
C                    Threshold for decomposition of virtual occupied block
                     READ(LUCMD,*) THACVI
                  GOTO 100
    9             CONTINUE
C                    '.SELDIA'
C                    Select which diagonals to decompose in one atom
                     READ(LUCMD,*) IAT1,NSLTMP
                     NSELD(IAT1) = NSLTMP
                     READ(LUCMD,*) (ISELD(I,IAT1),I=1,NSELD(IAT1))
                     DIASEL(IAT1) = .TRUE.
                  GOTO 100
   10             CONTINUE
C                    '.ORDER'
C                    Order to decompose atoms
ctmp                 READ(LUCMD,*) (IORDEC(I), I=1,NUCIND)
                  GOTO 100
   11             CONTINUE             
C                    '.FULDEC'
C                    Select all the atoms
                     FULDEC = .TRUE.
                     SELDIR = .TRUE.
                  GOTO 100
C                    '.XXXXXX'
   12             CONTINUE             
   13             CONTINUE             
   14             CONTINUE             
   15             CONTINUE             
                  GOTO 100

C
               ELSE
                  WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
     &                '" not recognized in ',SECNAM,'.'
                  CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',
     &                       LUPRI)
                  CALL QUIT('Illegal Keyword in '//SECNAM)
               END IF
C
            ELSE IF (WORD(1:1) .NE. '*') THEN
               WRITE (LUPRI,'(/5A,/)') 'PROMPT "',WORD,
     &         '" NOT RECOGNIZED IN ',SECNAM,'.'
               CALL QUIT('Illegal prompt in '//SECNAM)
            ELSE IF (WORD(1:1) .EQ.'*') THEN
               BACKSPACE (LUCMD)
               GO TO 200
            END IF
C
      END IF
C
  200 CONTINUE
C
C     Check that input makes sense
C
      WRITE(LUPRI,'(//,A)') '       ----------------------------'
      WRITE(LUPRI,'(A)') '          Info from Selact input'
      WRITE(LUPRI,'(A,/)') '       ----------------------------'
C
      IF (ACTIVE .AND. (.NOT.  CHKACT)) THEN
         WRITE(LUPRI,'(//,A,/,3A,/,A,//)') '>>> ERROR :',
     &                '>>> *ACTIVE was specified, but no definition',
     &                ' of active atoms, boxes or basis',
     &                ' functions was given',
     &                '    Program will stop'
         CALL QUIT(SECNAM//' called but no selection data given')
      END IF
C
      IF (ATOMIC) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(I5,A)') NACINP, ' centers declared as active :'
         WRITE(LUPRI,'(14I5)') (LACINP(I), I=1,NACINP)
      END IF
C
      IF (SELDIR) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Density decomposition only in active space'
      ELSE
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Decompose in atom-by-atom basis and ',
     &                  'select afterwards'
      END IF
C
      WRITE(LUPRI,*)
      WRITE(LUPRI,*)
C
      RETURN
      END
