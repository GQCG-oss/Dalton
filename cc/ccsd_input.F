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
#include <single.h>
C
C  /* Deck cc_inp */
      SUBROUTINE CC_INP(WORD,IREST,MSYM)
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
C     Keywords have been added for R12 method (WK/UniKA/04-11-2002).
      PARAMETER (NTABLE = 98)
      LOGICAL   SET, NEWDEF, SIRFF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER*(80) LINE
C
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
C
      SAVE SET
      DATA TABLE /'.SKIP  ','.PRINT ','.DIRECT','.RESTAR','.CC3   ',
     &            '*CCEXCI','*CCLRSD','.CCSTST','.NSYM  ','.MAXRED',
     &            '*CCEXGR','.R1SKIP','.L1SKIP','.RESKIP','.LESKIP',
     &            '.F1SKIP','.MAX IT','.XXXXXX','.E0SKIP','.L0SKIP',
     &            '.LISKIP','.CC2   ','.MP2   ','.CC(2) ','.CC1B  ',
     &            '.CC(T) ','.CC(3) ','.CCS   ','.FCORE ','.FSECON',
     &            '.CCD   ','.CC1A  ','.CIS   ','.THRENR','.NOCCIT',
     &            '.IMSKIP','.M1SKIP','.FRSKIP','.MINSCR','.MINMEM',
     &            '.B0SKIP','.NEWCAU','.NOT2TC','*CCGR  ','.FROEXP',
     &            '.FROIMP','.MXDIIS','.CCSD  ','.CCR(A)','.CCR(B)',
     &            '*CCFOP ','.XXXXXX','.CCR(3)','.CCR(T)','.FIELD ',
     &            '.DEBUG ','*CCQR2R','.HERDIR','.BUFLEN','*CCLR  ',
     &            '*CCEXLR','.NSIMLE','.THRLEQ','.MXLRV ','*CCTM  ',
     &            '*CCSM  ','*CC5R  ','*CC4R  ','*CCQR  ','*CCCR  ',
     &            '.O2SKIP','.R2SKIP','.X2SKIP','.F2SKIP','.L2SKIP',
     &            '*CCMCD ','.ANAAOD','.PACK  ','.CONNEC','.THRLDP',
     &            '.RCSKIP','.FCSKIP','.LCSKIP','.CO2SKI','.CX2SKI',
     &            '.CR2SKI','.CF2SKI','.CL2SKI','*DERIVA','.N2SKIP',
     &            '.BRSKIP','.FREEZE','*CCSLV ','*R12   ','*R12 IN',
     &            '.PAIRS ','.ETAPTI','.DKABAR'/
      DATA SET/.FALSE./
C
      IF (SET) RETURN
      SET = .TRUE.
C
C     Initialize /CCSDINP/ ,/CCLR / and /CCSDSYM/
C
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
      CCSM   = .FALSE.
      CCTM   = .FALSE.
      CCMCD  = .FALSE.
      CCLRSD = .FALSE.
      CCQR2R = .FALSE.
      CCEXLR = .FALSE.
      CCDERI = .FALSE.

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
      B0SKIP = .FALSE.
      N2SKIP = .FALSE.
      BRSKIP = .FALSE.
      ETADSC = .FALSE.

      FREEZE = .FALSE.
      CCSLV  = .FALSE.
      NCCSLV = 0

      HERDIR = .FALSE.

      ANAAOD   = .FALSE.
      MAXRED   = 200
      MXLRV    = MAXRED - 1
      MXDIIS   = 8
      MAXITE   = 40
      THRLEQ   = 1.0D-05
      THRENR   = 1.0D-08
      THRLDPHF = 1.0D-13
      NFIELD = 0
      NSIMLE = 0
      NEWCAU = .FALSE.
      LBUF   = 250000
C
      KEEPAOTWO = 0
      KEEPAOIN  = .FALSE.
      LPACKINT  = .FALSE.
      THRPCKINT = 1.0D-15
C
      CONNECTION = 'SYMMETR'
C
      MSYMS = MSYM 
C----------------------------------------------------------------
C Initialize additional labels for CCSD(T)-gradient optimization
C and FOP
C----------------------------------------------------------------

      ETACCPT = .FALSE.
      DIRKAPB = .FALSE.
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
C     (for the time being, on non-direct integrals available)
C
      DIRGRD = .FALSE.
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
      CALL DZERO(EFIELD,MXFELT)
C
      ITEST  = 0
C
      MINSCR = .TRUE.
      MINMEM = .FALSE.
      CCPAIR = .FALSE.
ckr      R12CAL = .TRUE.
ckr      R12NOA = .FALSE.
ckr      R12NOP = .FALSE.
ckr      R12NOB = .FALSE.
ckr      R12HYB = .TRUE.
ckr      NORXR  = .FALSE.
      VCLTHR =  0D0
      SVDTHR =  1D-15
ckr      R12XXL = .FALSE.
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
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD

            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
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
     *                92,93,94,94,95,96,97), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                GO TO 100
               END IF
               IF (WORD(1:1) .EQ. '*') THEN
cch
c                either '*END OF' for '**CC   ' section or a sirius
c                keyword. in the former case read next input line
cch
                 IF (WORD.EQ.'*END OF' .AND. WORD1.EQ.'**CC   ') THEN
                   READ (LUCMD, '(A7)') WORD
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
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, '(I5)') IPRINT
               GO TO 100
    3          CONTINUE
                  DIRECT = .TRUE.
               GO TO 100
    4          CONTINUE
                  CCRSTR = .TRUE.
               GO TO 100
    5          CONTINUE
                  CC3    = .TRUE.
                  CCSDT  = .TRUE.
               GO TO 100
    6          CONTINUE
                  CALL CC_EXCINP(WORD,MSYM)
               GO TO 100
    7          CONTINUE
                  CALL CC_LRSINP(WORD,MSYM)
               GO TO 100
    8          CONTINUE
                  CCSTST = .TRUE.
               GO TO 100
    9          CONTINUE
                  READ (LUCMD, '(I5)') MSYM2 
                  IF (((MSYM.LT.8).AND.(MSYM.GT.0))
     *               .AND.(MSYM2.NE.MSYM)) THEN
                    CALL QUIT(' Symmetry mismatch in input')
                  ELSE
                     MSYM  = MSYM2
                  ENDIF 
               GO TO 100
   10          CONTINUE
                  READ (LUCMD, '(I5)') MAXRED
               GO TO 100
   11          CONTINUE
                  CALL CC_EXGRIN(WORD,MSYM)
               GO TO 100
   12          CONTINUE
                  R1SKIP = .TRUE. 
               GO TO 100
   13          CONTINUE
                  L1SKIP = .TRUE. 
               GO TO 100
   14          CONTINUE
                  RESKIP = .TRUE. 
               GO TO 100
   15          CONTINUE
                  LESKIP = .TRUE. 
               GO TO 100
   16          CONTINUE
                  F1SKIP = .TRUE.
               GO TO 100
   17          CONTINUE
                  READ (LUCMD, '(I5)') MAXITE
               GO TO 100
   18          CONTINUE
               GO TO 100
   19          CONTINUE
                  E0SKIP = .TRUE.
               GO TO 100
   20          CONTINUE
                  L0SKIP = .TRUE.
               GO TO 100
   21          CONTINUE
                  LISKIP = .TRUE.
               GO TO 100
   22          CONTINUE
                  CC2   = .TRUE.
               GO TO 100
   23          CONTINUE
                  MP2   = .TRUE.
               GO TO 100
   24          CONTINUE
                  CCP2  = .TRUE.
               GO TO 100
   25          CONTINUE
                  CCSDT = .TRUE.
                  CC1B  = .TRUE.
               GO TO 100
   26          CONTINUE
                  CCPT  = .TRUE.
               GO TO 100
   27          CONTINUE
                  CCP3  = .TRUE.
               GO TO 100
   28          CONTINUE
                  CCS   = .TRUE.
               GO TO 100
   29          CONTINUE
                  LCOR  = .TRUE.
                  READ (LUCMD,*) (ICOR(ISYM),ISYM=1,MSYM)
               GO TO 100
   30          CONTINUE
                  LSEC  = .TRUE.
                  READ (LUCMD,*) (ISEC(ISYM),ISYM=1,MSYM)
               GO TO 100
   31          CONTINUE
                  CCD = .TRUE.
               GO TO 100
   32          CONTINUE
                  CCSDT = .TRUE.
                  CC1A  = .TRUE.
               GO TO 100
   33          CONTINUE
                  CIS = .TRUE.
               GO TO 100
   34          CONTINUE
                  READ (LUCMD, *) THRENR
               GO TO 100
   35          CONTINUE
                  NOCCIT = .TRUE.
               GO TO 100
   36          CONTINUE
                  IMSKIP = .TRUE.
               GO TO 100
   37          CONTINUE
                  M1SKIP = .TRUE.
               GO TO 100
   38          CONTINUE
                  FRSKIP = .TRUE.
               GO TO 100
   39          CONTINUE
                  READ (LUCMD, *) MINSCR
                  ITEST = ITEST + 1
               GO TO 100
   40          CONTINUE
                  READ (LUCMD, *) MINMEM
                  ITEST = ITEST + 1
               GO TO 100
   41          CONTINUE
                  B0SKIP = .TRUE.
               GO TO 100
   42          CONTINUE
                  NEWCAU = .TRUE.
               GO TO 100
   43          CONTINUE
                  T2TCOR = .FALSE.
               GO TO 100
   44          CONTINUE
                  CALL CC_GRIN(WORD,MSYM)
               GO TO 100
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
               GO TO 100
   46          CONTINUE
                  FROIMP = .TRUE.
                  IF (FROEXP) FROEXP = .FALSE.
                  IF (FREEZE) CALL QUIT(' Only one of FREEZE - FROEXP')
                  READ(LUCMD,*) (NRHFFR(I),I=1,MSYM)
                  READ(LUCMD,*) (NVIRFR(I),I=1,MSYM)
               GO TO 100
   47          CONTINUE
                  READ(LUCMD,*) MXDIIS
               GO TO 100
   48          CONTINUE
                  CCSD   = .TRUE.
               GO TO 100
   49          CONTINUE
                  CCR1A  = .TRUE.
               GO TO 100
   50          CONTINUE
                  CCR1B  = .TRUE.
               GO TO 100
   51          CONTINUE
                  CALL CC_FOPINP(WORD)
               GO TO 100
   52          CONTINUE
               GO TO 100
   53          CONTINUE
                  CCR3   = .TRUE.
                  CCSD   = .TRUE.
               GO TO 100
   54          CONTINUE
                  CCRT   = .TRUE.
               GO TO 100
   55          CONTINUE
                  NFIELD = NFIELD + 1
                  IF (NFIELD .LE. MXFELT) THEN
                     READ(LUCMD,*) EFIELD(NFIELD)
                     READ(LUCMD,'(1X,A8)') LFIELD(NFIELD)
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
               GO TO 100
   56          CONTINUE
                  DEBUG = .TRUE.
               GO TO 100
   57          CONTINUE
                  CALL CC_QR2RINP(WORD)
               GO TO 100
   58          CONTINUE
                  HERDIR = .TRUE.
               GO TO 100
   59          CONTINUE
                  READ (LUCMD, *) LBUF
               GO TO 100
   60          CONTINUE
                  CALL CC_LRINP(WORD)
               GO TO 100
   61          CONTINUE
                  CALL CC_EXLRINP(WORD)
               GO TO 100
   62          CONTINUE
                  READ (LUCMD, *) NSIMLE
               GO TO 100
   63          CONTINUE
                  READ (LUCMD, *) THRLEQ
               GO TO 100
   64          CONTINUE
                  READ (LUCMD, *) MXLRV
               GO TO 100
   65          CONTINUE
                  CALL CC_TMINP(WORD)
               GO TO 100
   66          CONTINUE
                  CALL CC_SMINP(WORD)
               GO TO 100
   67          CONTINUE
                  CALL CC_5RINP(WORD)
               GO TO 100
   68          CONTINUE
                  CALL CC_4RINP(WORD)
               GO TO 100
   69          CONTINUE
                  CALL CC_QRINP(WORD)
               GO TO 100
   70          CONTINUE
                  CALL CC_CRINP(WORD)
               GO TO 100
   71          CONTINUE
                  O2SKIP = .TRUE.
               GO TO 100
   72          CONTINUE
                  R2SKIP = .TRUE.
               GO TO 100
   73          CONTINUE
                  X2SKIP = .TRUE.
               GO TO 100
   74          CONTINUE
                  F2SKIP = .TRUE.
               GO TO 100
   75          CONTINUE
                  L2SKIP = .TRUE.
               GO TO 100
   76          CONTINUE
                 CALL CC_MCDINP(WORD)
               GO TO 100
   77          CONTINUE
                  ANAAOD = .TRUE.
               GO TO 100
   78          CONTINUE
C                '.PACK  ' 
                  LPACKINT = .TRUE.
                  READ (LUCMD, *) THRPCKINT
               GO TO 100
   79          CONTINUE
C                '.CONNEC' 
c                 CONNECTION = 'SYMMETR' / 'NATURAL' 
                  READ (LUCMD, '(A7)') CONNECTION
               GO TO 100
   80          CONTINUE
C                '.THRLDP' 
                  READ (LUCMD, *) THRLDPHF
               GO TO 100
   81          CONTINUE
C              '.RCSKIP' 
                 RCSKIP = .TRUE.
               GO TO 100
   82          CONTINUE
C              '.FCSKIP' 
                 FCSKIP = .TRUE.
               GO TO 100
   83          CONTINUE
C              '.LCSKIP' 
                 LCSKIP = .TRUE.
               GO TO 100
   84          CONTINUE
C              '.CO2SKI' 
                 CO2SKIP = .TRUE.
               GO TO 100
   85          CONTINUE
C              '.CX2SKI' 
                 CX2SKIP = .TRUE.
               GO TO 100
   86          CONTINUE
C              '.CR2SKI' 
                 CR2SKIP = .TRUE.
               GO TO 100
   87          CONTINUE
C              '.CF2SKI' 
                 CF2SKIP = .TRUE.
               GO TO 100
   88          CONTINUE
C              '.CL2SKI' 
                 CL2SKIP = .TRUE.
               GO TO 100
   89          CONTINUE
C                '*DERIVA' 
                 CCDERI = .TRUE.
                 RELORB = .TRUE.
               GO TO 100
   90          CONTINUE
C                '.N2SKIP'
                  N2SKIP = .TRUE.
               GO TO 100
   91          CONTINUE
C                '.BRSKIP'
                  BRSKIP = .TRUE.
               GO TO 100
   92          CONTINUE
C                '.FREEZE'
                  FREEZE = .TRUE.
                  IF (FROIMP.OR.FROEXP)
     *            CALL QUIT(' Only one of FREEZE - FROEXP - FROIMP')
                  READ(LUCMD,*) NFC,NFV
                  FROIMP = .TRUE.
               GO TO 100
   93          CONTINUE
C                '*CCSLV '
                  CALL CC_SLVINP(WORD)
               GO TO 100
   94          CONTINUE
C                '*R12   ' OR '*R12 IN'
                 CALL CC_R12IN(WORD)
               GO TO 100
   95          CONTINUE
C                '.PAIRS '
                 CCPAIR = .TRUE.
               GO TO 100
   96          CONTINUE
C                '.ETAPTI', extra integrals for CCSD(T) geopt
                  ETACCPT = .TRUE.
               GO TO 100
   97          CONTINUE
C                '.DKABAR', direct KappaBar calculation in nondir CC
                  DIRKAPB = .TRUE.
               GO TO 100
            ELSE
               WRITE (LUPRI,'(/3A,/)') ' Prompt "',WORD,
     *            '" not recognized in CC2INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in CC2INP.')
            END IF
      END IF
  300 CONTINUE
C
      MSYM = MSYMS
      IF (ITEST .EQ. 0 ) THEN
         MINSCR = .TRUE.
         IF (DIRECT) MINSCR = .FALSE.
         MINMEM = .FALSE.
         IF (DIRECT) MINMEM = .FALSE.
      ENDIF
      IF ( .NOT. MINSCR ) MINMEM = .TRUE.
      RETURN
      END
C  /* Deck CC_PRTI */
      SUBROUTINE CC_PRTI
C
C     K.Ruud, Jan.-96. Split from CC_INP in order to place CC output more 
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
         CALL HEADER('Changes of defaults for CC:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' -CCSD skipped in this run.'
         ELSE
C
            IF (DIRECT) WRITE (LUPRI,'(/A,/A)')
     *         ' -This is a direct atomic orbital integral based',
     *         '  calculation for coupled cluster wave function'
C
            IF (CCSDT) WRITE (LUPRI,'(/A)')
     *         ' -Iterative triple excitations included '
C
            IF (FROIMP) WRITE(LUPRI,'(/A)')
     *         ' -Implicit frozen core calculation'
C
            IF (FROEXP) WRITE(LUPRI,'(/A)')
     *         ' -Explicit frozen core calculation'
C
            IF (.NOT. T2TCOR) WRITE(LUPRI,'(/A)')
     *         ' -Transposed t2-amplitudes not hold in core'
C
            IF (NFIELD.GT.0) THEN
               WRITE(LUPRI,'(A)') ' -Calculation with finite fields:'
               DO IF = 1, NFIELD
                  WRITE(LUPRI,'(1X,A,F10.6,A,1X,A8)')
     *                       'Field strenths: ',EFIELD(IF),
     *                       'Field Label: ',LFIELD(IF)
               END DO
               CALL FLDINP
               IF (CCSDT.AND.NONHF) THEN
                  IF (.NOT. CC3) THEN
                    WRITE(LUPRI,*) 'No triples unrelaxed FF possibility'
                    CALL QUIT('No triples unrelaxed FF possibility')
                  ENDIF
               ENDIF
               IF (CCS.AND.NONHF) THEN
                  WRITE(LUPRI,*) 'No CCS unrelaxed FF possibility '
                  WRITE(LUPRI,*) 'Use instead CC2 with CCSTST option.'
                  CALL QUIT('No CCS unrelaxed FF possibility')
               ENDIF                                             
               IF (NONHF .AND. RELORB) THEN
                  WRITE(LUPRI,*) 'Inconsistency: Non HF reference and '
     *                        //'relaxed derivative requested'
                  CALL QUIT('Inconsistency: in FF '//
     &                 'and relaxed derivative')
               ENDIF
            END IF 
C
            IF (LHTR.AND.((CCLRSD).OR.
     *         (CCR3.OR.CCRT.OR.CCR1A.OR.CC1B.OR.CC1A.OR.CC1B)))
     *      THEN
               WRITE(LUPRI,*) 'Input inconsistent due to LHTR '
               CALL QUIT('Do not use LHTR for this '//
     &              'type of calculation ')
            ENDIF
C
C
            IF ( DEBUG ) WRITE(LUPRI,'(/A)')
     *         ' Debug printout activated '
            IF ( CCEXCI ) WRITE(LUPRI,'(/A)')
     *         ' Excitation energies calculated'
            IF ( CCLRSD ) WRITE(LUPRI,'(/A)')
     *         ' Residues of the linear response functions calculated'
            IF ( JACTST ) WRITE(LUPRI,'(/A)')
     *         ' Jacobian tested agains finite difference Jacobian'
            IF ( JACEXP ) WRITE(LUPRI,'(/A)')
     *         ' Jacobian constructed explicit'
            IF ( FDEXCI ) WRITE(LUPRI,'(/A)')
     *         ' Excitation energies of finite diff. Jacobian calc.'
            IF ( CCLR  ) WRITE(LUPRI,'(/A)')
     *         ' Linear response properties calculated'
            IF ( CAUCHY) WRITE(LUPRI,'(/A)')
     *         ' Dispersion coefficients for linear response calc.'
            IF ( CCQR  ) WRITE(LUPRI,'(/A)')
     *         ' Quadratic response properties calculated'
            IF ( CCCR  ) WRITE(LUPRI,'(/A)')
     *         ' Cubic response properties calculated'
            IF ( OSCSTR) WRITE(LUPRI,'(/A)')
     *         ' Oscillator strengths calculated'
            IF ( CCQR2R) WRITE(LUPRI,'(/A)')
     *        ' Transition strengths between excited states calculated.'
            IF (CCEXGR) WRITE(LUPRI,'(/A)')
     *         ' Excited state properties calculated '
            IF ( CCMCD  ) WRITE(LUPRI,'(/A)')
     *         ' Magnetic circular dichroism B calculated'
            IF (DIPMOM) WRITE(LUPRI,'(/A)')
     *         ' Dipole moment calculated'
            IF (QUADRU) WRITE(LUPRI,'(/A)')
     *         ' Traceless quadrupole moment calculated'
            IF (NQCC) WRITE(LUPRI,'(/A)')
     *         ' Electric field gradient calculated'
            IF (RELCOR) WRITE(LUPRI,'(/A)')
     *         ' Relativistic corrections to energy calculated'
            IF (SECMOM) WRITE(LUPRI,'(/A)')
     *         ' Electronic second moment of charge calculated'
            IF (DAR2EL) WRITE(LUPRI,'(/A)')
     *         ' Relativistic two-electron Darwin term calculated'
            IF (DPTECO) WRITE(LUPRI,'(/A)')
     *         ' First-order DPT energy corrections calculated'

C
         END IF
         WRITE (LUPRI,'(A)') '  '
      END IF
C
      RETURN
      END
c /* deck fldinp */
      SUBROUTINE FLDINP
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
         WRITE(LUPRI,'(/,A,3F14.8)') 'Electric field:  ',
     *                           (ELEFLD(I),I=1,3)
         WRITE(LUPRI,'(A,F14.8)') 'Total norm:      ',ELFLNR
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
      WRITE(LUPRI,'(A,3F14.8)') 'Center of charge:',XCQ,YCQ,ZCQ
      WRITE(LUPRI,'(2A,F14.8)') 'Charge correction to interaction ',
     *                      'with electric field:', CORR
      WRITE(LUPRI,*)
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
C    Purpose: Read input for CC excitec state calculations.
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
      PARAMETER (NTABLE = 19)

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
     *            '.OMEINP','.STVEC ','.STOLD ','.XXXXXX','.THREXC',
     *            '.CCSPIC','.CC2PIC','.CCSDPI','.MARGIN'/

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
      XMARGIN = 1.0
C
      CALL IZERO(NCCEXCI,3*8)
      CALL IZERO(NOMINP,3*8)
C
C     Other initializations
C
      ECURR  = 0.0D0
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
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
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
            GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
     *              IJUMP
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
C
C--------------------
C           .XXXXXX : unused label
C--------------------
C
14          CONTINUE
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
          ELSE
           WRITE (LUPRI,'(/,5A,/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/,5A,/)') 'PROMPT "',WORD,
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
      IF (NEXCI .EQ. 0) OSCSTR = .FALSE.
      IF (NEXCI .EQ. 0) WRITE(LUPRI,'(/1X,A)')
     *'No excitation energy requested even though CCEXCI is set-right?'
C                                                          

      RETURN
      END
C=====================================================================*
c/* deck cc_lrsinp */
       SUBROUTINE CC_LRSINP(WORD,MSYM)
C---------------------------------------------------------------------*
C
C    Purpose: Read input for CC excitec state calculations.
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
      PARAMETER (NTABLE = 7)

* variables:
      LOGICAL SET
      SAVE SET

      CHARACTER WORD*(7)
      CHARACTER LABELA*(8),LABELB*(8),LABHELP*70
      CHARACTER TABLE(NTABLE)*(8)

      INTEGER IJUMP,IDIP(3)
* data:
      DATA SET /.FALSE./
      DATA TABLE /'.DIPOLE','.XXXXXX','.XXXXXX','.NO2N+1','.OPERAT',
     *            '.SELEXC','.DIPVEL'/

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
      LRS2N1 = .TRUE. 
      SELLRS = .FALSE.
C
      NSELRS   = 0
      NLRSOP   = 0
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
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
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
            GOTO (1,2,3,4,5,6,7), IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_LRRINP.')

C
C-------------------------------------------------
C           Calculate dipole oscillator strengths.
C-------------------------------------------------
C
1           CONTINUE
              OSCSTR = .TRUE.
              IF (NLRSOP+9 .GT. MXLRSO) THEN
                WRITE(LUPRI,'(2(/,A,I5),/)')
     &          ' NB. OF OPERATOR DOUBLES SPECIFIED  : ',NLRSOP+9,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSO
                CALL QUIT('TO MANY OPERATOR DOUBLES IN CCLRSD.')
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
C-----------------------------------
C           
C-----------------------------------
C
2           CONTINUE
            GO TO 100
C
C-----------------------------------
C           
C-----------------------------------
C
3           CONTINUE
            GO TO 100
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
                    WRITE(LUPRI,'(/2A,I5,/)')
     &               ' NB. OF OPERATOR DOUBLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSO
                    CALL QUIT('TO MANY OPERATOR DOUBLES IN CCLRSD.')
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
                    WRITE(LUPRI,'(/2A,I5,/)')
     &               ' NB. OF STATES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSST
                    CALL QUIT('TO MANY STATES IN CCLRS.')
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
7           CONTINUE
              VELSTR = .TRUE.
              IF (NLRSOP+9 .GT. MXLRSO) THEN
                WRITE(LUPRI,'(2(/,A,I5),/)')
     &          ' NB. OF OPERATOR DOUBLES SPECIFIED  : ',NLRSOP+9,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLRSO
                CALL QUIT('TO MANY OPERATOR DOUBLES IN CCLRSD.')
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

          ELSE
           WRITE (LUPRI,'(/,5A,/)') ' Prompt "',WORD,
     &             '" not recognized in ',SECNAM,'.'
           CALL PRTAB(NTABLE,TABLE,SECNAM//' input keywords',LUPRI)
           CALL QUIT('Illegal Keyword in '//SECNAM//'.')
          END IF

        ELSE IF (WORD(1:1) .NE. '*') THEN
          WRITE (LUPRI,'(/,5A,/)') 'PROMPT "',WORD,
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
        IF (SELLRS.AND.(NSELRS.EQ.0)) WRITE(LUPRI,'(/A)')
     &     '(*CCLRSD input is strange - no states is requested.)'
        IF (NLRSOP .EQ.0) WRITE(LUPRI,'(/A)')
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
c/* deck cc_qr2rinp */
C=====================================================================*
       SUBROUTINE CC_QR2RINP(WORD,MSYM)
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
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
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
                WRITE(LUPRI,'(2(/A,I5),/)')
     &          ' NB. OF OPERATOR DOUBLES SPECIFIED  : ',NQR2OP+9,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQR2O
                CALL QUIT('TO MANY OPERATOR DOUBLES IN CCQR2R.')
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
                    WRITE(LUPRI,'(/2A,I5,/)')
     &               ' NB. OF OPERATOR DOUBLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQR2O
                    CALL QUIT('TO MANY OPERATOR DOUBLES IN CCQR2R.')
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
                    WRITE(LUPRI,'(/2A,I5,/)')
     &               ' NB. OF STATES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQR2ST
                    CALL QUIT('TO MANY STATES IN CCQR2R.')
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
                WRITE(LUPRI,'(2(/A,I5),/)')
     &          ' NB. OF OPERATOR DOUBLES SPECIFIED  : ',NQR2OP+9,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQR2O
                CALL QUIT('TO MANY OPERATOR DOUBLES IN CCQR2R.')
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
           WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
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
*---------------------------------------------------------------------*
* check, if input consistent.
*---------------------------------------------------------------------*
C
        IF (SELQR2.AND.(NSEQR2 .EQ.0)) WRITE(LUPRI,'(/A)')
     &     '(*CCQR2R input is strange - no states is requested.)'
        IF (NQR2OP .EQ.0) WRITE(LUPRI,'(/A)')
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
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
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
           WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
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
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
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
                    WRITE(LUPRI,'(/2A,I5,/)') 
     &               ' NB. OF OPERATORS SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXGROP
                    CALL QUIT('TO MANY OPERATORS IN CCEXGR.')
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
                    WRITE(LUPRI,'(/2A,I5,/)')
     &               ' NB. OF STATES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXXGRST
                    CALL QUIT('TO MANY STATES IN CCEXGR.')
                  END IF
                END IF
                READ (LUCMD,'(A70)') LABHELP
              END DO
              BACKSPACE(LUCMD)
            GO TO 100
C
          ELSE
           WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
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
      PARAMETER (NTABLE = 11)

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
     *            '.DPTECO'/

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
C
*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCFOP ') THEN

100   CONTINUE
       
* get new input line:
        READ (LUCMD,'(A7)') WORD
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
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
            GOTO (1,2,3,4,5,6,7,8,9,10,11), IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_FOPINP.')
      
C
C-----------------------------------
C           Calculate dipole moment.
C-----------------------------------
C

1           CONTINUE
               DIPMOM = .TRUE.
            GO TO 100

C
C----------------------------------------
C           Calculate Quadrupole moments.
C----------------------------------------
C

2           CONTINUE
               QUADRU = .TRUE.
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
                    WRITE(LUPRI,'(/2A,I5,/)')
     &               ' NB. OF OPERATORS SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MAFOP 
                    CALL QUIT('TO MANY OPERATORS IN CCFOP .')
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
C           Relativistic two-electron Darwin term.
C-------------------------------------------------
C

11         CONTINUE
               DPTECO = .TRUE.
            GO TO 100
C
          ELSE
           WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
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
     *           SECMOM.OR.DAR2EL.OR.DPTECO.OR.(NAFOP.GT.0))
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

* local parameters:
      CHARACTER SECNAM*(8)
      PARAMETER (SECNAM='CC_LRINP')

      INTEGER NTABLE
      PARAMETER (NTABLE = 15)

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
      CHARACTER TABLE(NTABLE)*(8)

      LOGICAL LRELAX
      INTEGER IDX, IJUMP, IDIP(3)

* external function:
      INTEGER INDPRPCC

* data:
      DATA SET /.FALSE./
      DATA TABLE /'.RELAXE','.UNRELA','.FREQUE','.DIPOLE','.ALLDSP',
     *            '.OLD_LR','.ASYMSD','.DISPCF','.OPERAT','.AVERAG',
     *            '.PRINT ','.STATIC','.XXXXXX','.XXXXXX','.XXXXXX'/ 

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
      IPRSOP = IPRINT
      DIPPOL = .FALSE.
      ASYMSD = .FALSE.
      LRELAX = .FALSE.
      ICHANG = 0

      OLDLR  = .FALSE.

*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCLR  ') THEN

100   CONTINUE
       
* get new input line:
        READ (LUCMD,'(A7)') WORD
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
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
            GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), IJUMP
            CALL QUIT('Illegal address in computed GOTO in CC_LRINP.')
C
C-----------------------
C           .RELAXEd
C-----------------------
C
1           CONTINUE
               LRELAX    = .TRUE.
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
C-----------------------
C           .FREQUE
C-----------------------
C
3           CONTINUE
              READ (LUCMD,*) NBLRFR
              IF (NBLRFR .GT. MBLRFR) THEN
                WRITE(LUPRI,'(3(/A,I5),/)') 
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',NBLRFR,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MBLRFR,
     &          ' THE NUMBER IS RESET TO THE MAXIMUM : ',MBLRFR
                NBLRFR = MBLRFR
              END IF
              READ (LUCMD,*) (BLRFR(I),I=1,NBLRFR)
            GO TO 100
C
C-----------------------------------
C           Dipole polarizabilities.
C-----------------------------------
C
4           CONTINUE
              DIPPOL =.TRUE. 
              IF (NLROP+9 .GT. MXLROP) THEN
                WRITE(LUPRI,'(2(/A,I5),/)') 
     &          ' NB. OF OPERATOR DOUBLES SPECIFIED  : ',NLROP+9,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLROP
                CALL QUIT('TO MANY OPERATOR DOUBLES IN CCLR.')
              END IF
              IDIP(1) = INDPRPCC('XDIPLEN ')
              IDIP(2) = INDPRPCC('YDIPLEN ')
              IDIP(3) = INDPRPCC('ZDIPLEN ')
              DO IDXA=1,3
                 DO IDXB=1,3
                   IDX = NLROP + (IDXA-1)*3+IDXB
                   IALROP(IDX) = IDIP(IDXA)
                   IBLROP(IDX) = IDIP(IDXB)
                   LALORX(IDX) = LRELAX
                   LBLORX(IDX) = LRELAX
                 END DO
              END DO
              NLROP = NLROP + 9 
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
            GO TO 100

C           -------------------------------------------
C           .OPERAT : Operator set for Linear response.
C           -------------------------------------------
 
9           CONTINUE
              READ (LUCMD,'(2A)') LABELA, LABELB
              DO WHILE (LABELA(1:1).NE.'.' .AND. LABELA(1:1).NE.'*')
                IF (LABELA(1:1).NE.'!') THEN
                  IF (NLROP.LT.MXLROP) THEN
                    NLROP = NLROP + 1
                    IALROP(NLROP) = INDPRPCC(LABELA)
                    IBLROP(NLROP) = INDPRPCC(LABELB)
                    LALORX(NLROP) = LRELAX
                    LBLORX(NLROP) = LRELAX
                  ELSE
                    WRITE(LUPRI,'(/2A,I5,/)')
     &               ' NB. OF OPERATOR DOUBLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXLROP
                    CALL QUIT('TO MANY OPERATOR DOUBLES IN CCLR.')
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
                WRITE(LUPRI,'(/4A,/4A,/)')
     &           'LABEL "',LINE(1:5),'" UNKNOWN FOR .AVERAG KEYWORD',
     &           'IN *CCLR SECTION.','INPUT IS IGNORED...'
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
              IF (NBLRFR .GE. MBLRFR) THEN
                WRITE(LUPRI,'(3(/A,I5),/)') 
     &          ' NUMBER OF FREQUENCIES SPECIFIED         : ',NBLRFR,
     &          ' IS ALREADY EQUAL TO THE MAXIMUM ALLOWED : ',MBLRFR,
     &          ' THE .STATIC KEYWORD UNDER *CCLR WILL BE IGNORED...'
              ELSE
                NBLRFR = NBLRFR + 1
                BLRFR(NBLRFR) = 0.0D0
              END IF
            GO TO 100
C                                                                     
C           -------------------------------------------------
C           .XXXXXXX unused keywords
C           -------------------------------------------------
C
13           CONTINUE
14           CONTINUE
15           CONTINUE
            GO TO 100                 
C
          ELSE
           WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
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
     &            '.RELAXE','.UNRELA','.XXXXXX','.XXXXXX','.XXXXXX' /

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
      ALLDSPCF     = .FALSE.
      LALRX        = .FALSE.
      LBLRX        = .FALSE.
      LCLRX        = .FALSE.
      LRELAX       = .FALSE.    

      IPRQHYP = 0

      ICHANG = 0

*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCQR  ') THEN

100   CONTINUE
       
* get new input line:
        READ (LUCMD,'(A7)') WORD
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
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
                    WRITE(LUPRI,'(/2A,I5,/)')
     &               ' NB. OF OPERATOR TRIPLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQROP
                    CALL QUIT('TO MANY OPERATOR TRIPLES IN CCQR.')
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
                WRITE(LUPRI,'(2(/A,I5),/)') 
     &          ' NB. OF OPERATOR TRIPLES SPECIFIED  : ',NQROPER+27,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQROP
                CALL QUIT('TO MANY OPERATOR TRIPLES IN CCQR.')
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
                      WRITE(LUPRI,'(/2A,/A)')
     &                 ' NEGATIVE EXPANSION COEFFICIENTS NOT',
     &                 ' AVAILABLE FOR FIRST HYPERPOLARIZABILITIES.',
     &                 ' INPUT LINE IGNORED...'
                    ELSE
                      NQRDISP = NQRDISP + 1
                      IQCAUA(NQRDISP) = ICA
                      IQCAUB(NQRDISP) = ICB
                      IQCAUC(NQRDISP) = ICC
                    END IF
                  ELSE
                    WRITE(LUPRI,'(/2A,I5,/)') 
     &               ' NB. OF EXPANSION COEFFICIENTS ',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQRDISP
                    CALL QUIT('TO MANY EXPANSION COEFFICIENTS IN CCQR.')
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
                  WRITE(LUPRI,'(/2A,/1A,/)') 
     &             ' WARNING: INPUT FOR .DIPOL OR .OPERATOR OPTIONS',
     &             ' BEFORE THE .AVERAG OPTION',
     &             ' IN *CCQR SECTION WILL BE IGNORED.'
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
                WRITE(LUPRI,'(/4A,/A,/)')
     &           'LABEL "',LINE(1:8),'" UNKOWN FOR .AVERAG KEYWORD',
     &           'IN *CCQR SECTION.','INPUT IS IGNORED...'
              ENDIF 
            GO TO 100


C           ---------------------------
C           .MIXFRE : mixed frequencies
C                     wb, wc, wa=-wb-wc
C           ---------------------------
6           CONTINUE
              READ (LUCMD,*) MFREQ
              IF (NQRFREQ+MFREQ .GT. MXQRFR) THEN
                WRITE(LUPRI,'(3(/A,I5),/)') 
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',NQRFREQ+MFREQ,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQRFR,
     &          ' THE NUMBER IS RESET TO THE MAXIMUM : ',MXQRFR
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
                WRITE(LUPRI,'(3(/A,I5),/)') 
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',NQRFREQ+MFREQ,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQRFR,
     &          ' THE NUMBER IS RESET TO THE MAXIMUM : ',MXQRFR
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
                WRITE(LUPRI,'(3(/A,I5),/)') 
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',NQRFREQ+MFREQ,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQRFR,
     &          ' THE NUMBER IS RESET TO THE MAXIMUM : ',MXQRFR
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
                WRITE(LUPRI,'(3(/A,I5),/)') 
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',NQRFREQ+MFREQ,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQRFR,
     &          ' THE NUMBER IS RESET TO THE MAXIMUM : ',MXQRFR
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
                WRITE(LUPRI,'(3(/A,I5)/)') 
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',NQRFREQ+1,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQRFR,
     &          ' INPUT OPTION STATIC WILL BE IGNORED.'
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
                WRITE(LUPRI,'(/2A)') 
     &            ' WARNING: INPUT FOR .EXPCOF OPTION BEFORE .DISPCF',
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
                  WRITE(LUPRI,'(/2A,I5,/)') 
     &             ' NB. OF DISPERSION COEFFICIENTS NEEDED',
     &             ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQRDISP
                  WRITE(LUPRI,'(/A,I2,A)')
     &             'DISPERSION COEFFICIENTS OF ORDER',K,' ARE IGNORED.'
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

C           -----------------------
C           .XXXXXX : unused labels
C           -----------------------
18          CONTINUE
19          CONTINUE
20          CONTINUE
              WRITE (LUPRI,*) 'unused .XXXXXX label... ignored'
            GO TO 100

          ELSE
           WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
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

*---------------------------------------------------------------------*
* check, if any triples of operator labels specified:
* if not, use default: complete dipole-dipole-dipole tensor
*---------------------------------------------------------------------*
      IF (NQROPER .EQ. 0) THEN 
          IF (NQROPER+27 .GT. MXQROP) THEN
            WRITE(LUPRI,'(2(/A,I5),/)') 
     &      ' NB. OF OPERATOR TRIPLES SPECIFIED  : ',NQROPER+27,
     &      ' IS GREATER THAN THE ALLOWED NUMBER : ',MXQROP
            CALL QUIT('TO MANY OPERATOR TRIPLES IN CCQR.')
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
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
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
                    WRITE(LUPRI,'(/2A,I5,/)')
     &               ' NB. OF OPERATOR QUADRUPLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXCROP
                    CALL QUIT('TO MANY OPERATOR QUADRUPLES IN CCCR.')
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
                WRITE(LUPRI,'(/2A,I5,/)')
     &           ' NB. OF OPERATOR QUADRUPLES SPECIFIED',
     &           ' IS GREATER THAN THE ALLOWED NUMBER : ',MXCROP
                CALL QUIT('TO MANY OPERATOR QUADRUPLES IN CCCR.')
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
                WRITE(LUPRI,'(3(/A,I5)/)')
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',NCRFREQ+1,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXCRFR,
     &          ' INPUT OPTION STATIC WILL BE IGNORED.'
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
                WRITE(LUPRI,'(3(/A,I5)/)') 
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',NCRFREQ+MFREQ,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXCRFR,
     &          ' THE NUMBER IS RESET TO THE MAXIMUM : ',MXCRFR
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
                WRITE(LUPRI,'(3(/A,I5)/)')
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',NCRFREQ+MFREQ,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXCRFR,
     &          ' THE NUMBER IS RESET TO THE MAXIMUM : ',MXCRFR
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
                WRITE(LUPRI,'(3(/A,I5)/)')
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',NCRFREQ+MFREQ,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXCRFR,
     &          ' THE NUMBER IS RESET TO THE MAXIMUM : ',MXCRFR
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
                WRITE(LUPRI,'(3(/A,I5)/)')
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',NCRFREQ+MFREQ,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXCRFR,
     &          ' THE NUMBER IS RESET TO THE MAXIMUM : ',MXCRFR
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
                WRITE(LUPRI,'(3(/A,I5)/)')
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',NCRFREQ+MFREQ,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXCRFR,
     &          ' THE NUMBER IS RESET TO THE MAXIMUM : ',MXCRFR
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
                      WRITE(LUPRI,'(/2A,/A)')
     &                 ' NEGATIVE EXPANSION COEFFICIENTS NOT',
     &                 ' AVAILABLE FOR SECOND HYPERPOLARIZABILITIES.',
     &                 ' INPUT LINE IGNORED...'
                    ELSE
                      NCRDISP = NCRDISP + 1
                      ICCAUA(NCRDISP) = ICAUA
                      ICCAUB(NCRDISP) = ICAUB
                      ICCAUC(NCRDISP) = ICAUC
                      ICCAUD(NCRDISP) = ICAUD
                    END IF
                  ELSE
                    WRITE(LUPRI,'(/2A,I5,/)') 
     &               ' NB. OF EXPANSION COEFFICIENTS ',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXCRDISP
                    CALL QUIT('TO MANY EXPANSION COEFFICIENTS IN CCCR.')
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
                  WRITE(LUPRI,'(/2A,/1A,/)') 
     &             ' WARNING: INPUT FOR .DIPOL OR .OPERATOR OPTIONS',
     &             ' BEFORE THE .AVERAG OPTION',
     &             ' IN *CCCR SECTION WILL BE IGNORED.'
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
                WRITE(LUPRI,'(/2A)') 
     &            ' WARNING: INPUT FOR .EXPCOF OPTION BEFORE .DISPCF',
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
                  WRITE(LUPRI,'(/2A,I5,/)') 
     &             ' NB. OF DISPERSION COEFFICIENTS NEEDED',
     &             ' IS GREATER THAN THE ALLOWED NUMBER : ',MXCRDISP
                  WRITE(LUPRI,'(/A,I2,A)')
     &             'DISPERSION COEFFICIENTS OF ORDER',L,' ARE IGNORED.'
                  WRITE(LUPRI,'(/2A,I5)') 'FOR NEXT ORDER INCREASE ',
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
           WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
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
*---------------------------------------------------------------------*
* check, if any quadruples of operator labels specified:
* if not, use default: complete dipole-dipole-dipole-dipole tensor
*---------------------------------------------------------------------*
      IF (NCROPER .EQ. 0) THEN
          IF (NCROPER+81 .GT. MXCROP) THEN
            WRITE(LUPRI,'(2(/A,I5)/)')
     &      ' NB. OF OPERATOR QUADRUPLES SPECIFIED  : ',NCROPER+81,
     &      ' IS GREATER THAN THE ALLOWED NUMBER : ',MXCROP
            CALL QUIT('TO MANY OPERATOR QUADRUPLES IN CCQR.')
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
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
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
                    WRITE(LUPRI,'(/2A,I5,/)')
     &               ' NB. OF OPERATOR QUADRUPLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MX4ROP
                    CALL QUIT('TO MANY OPERATOR QUINTUPLES IN CC4R.')
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
                WRITE(LUPRI,'(/2A,I5,/)')
     &           ' NB. OF OPERATOR QUINTUPLES SPECIFIED',
     &           ' IS GREATER THAN THE ALLOWED NUMBER : ',MX4ROP
                CALL QUIT('TO MANY OPERATOR QUINTUPLES IN CC4R.')
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
                WRITE(LUPRI,'(3(/A,I5)/)')
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',N4RFREQ+1,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MX4RFR,
     &          ' INPUT OPTION STATIC WILL BE IGNORED.'
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
                WRITE(LUPRI,'(3(/A,I5)/)') 
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',N4RFREQ+MFREQ,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MX4RFR,
     &          ' THE NUMBER IS RESET TO THE MAXIMUM : ',MX4RFR
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
                WRITE(LUPRI,'(3(/A,I5)/)')
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',N4RFREQ+MFREQ,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MX4RFR,
     &          ' THE NUMBER IS RESET TO THE MAXIMUM : ',MX4RFR
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
           WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
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
*---------------------------------------------------------------------*
* check, if any quintuples of operator labels specified:
* if not, use default: complete dipole^5 tensor
*---------------------------------------------------------------------*
      IF (N4ROPER .EQ. 0) THEN
          IF (N4ROPER+243 .GT. MX4ROP) THEN
            WRITE(LUPRI,'(2(/A,I5)/)')
     &      ' NB. OF OPERATOR QUINTUPLES SPECIFIED  : ',N4ROPER+243,
     &      ' IS GREATER THAN THE ALLOWED NUMBER : ',MX4ROP
            CALL QUIT('TO MANY OPERATOR QUINTUPLES IN CC4R.')
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
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
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
                    WRITE(LUPRI,'(/2A,I5,/)')
     &               ' NB. OF OPERATOR QUADRUPLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MX5ROP
                    CALL QUIT('TO MANY OPERATOR QUINTUPLES IN CC5R.')
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
                WRITE(LUPRI,'(/2A,I5,/)')
     &           ' NB. OF OPERATOR HEXTUPLES SPECIFIED',
     &           ' IS GREATER THAN THE ALLOWED NUMBER : ',MX5ROP
                CALL QUIT('TO MANY OPERATOR HEXTUPLES IN CC5R.')
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
                WRITE(LUPRI,'(3(/A,I5)/)')
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',N5RFREQ+1,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MX5RFR,
     &          ' INPUT OPTION STATIC WILL BE IGNORED.'
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
                WRITE(LUPRI,'(3(/A,I5)/)') 
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',N5RFREQ+MFREQ,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MX5RFR,
     &          ' THE NUMBER IS RESET TO THE MAXIMUM : ',MX5RFR
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
                WRITE(LUPRI,'(3(/A,I5)/)')
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',N5RFREQ+MFREQ,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MX5RFR,
     &          ' THE NUMBER IS RESET TO THE MAXIMUM : ',MX5RFR
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
           WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
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
*---------------------------------------------------------------------*
* check, if any quintuples of operator labels specified:
* if not, use default: complete dipole^6 tensor
*---------------------------------------------------------------------*
      IF (N5ROPER .EQ. 0) THEN
          IF (N5ROPER+729 .GT. MX5ROP) THEN
            WRITE(LUPRI,'(2(/A,I5)/)')
     &      ' NB. OF OPERATOR QUINTUPLES SPECIFIED  : ',N5ROPER+729,
     &      ' IS GREATER THAN THE ALLOWED NUMBER : ',MX5ROP
            CALL QUIT('TO MANY OPERATOR QUINTUPLES IN CC5R.')
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
         WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *   '@ NUMBER OF SPECIFIED PROPERTIES EXCEED THE MAXIMUM ALLOWED',
     *   '@ MAXPRP =',MAXLBL,' NPRLBL= ',NPRLBL
         CALL QUIT(' INDPRPCC: TOO MANY PROPERTIES SPECIFIED')
      END IF

      PRPLBL(NPRLBL) = NEWLBL
      INDPRPCC       = NPRLBL

      RETURN
      END
*---------------------------------------------------------------------*
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

      CHARACTER WORD*(7), LABHELP*(80)
      CHARACTER*8 LABELA, LABELB
      CHARACTER TABLE(NTABLE)*(7)

      INTEGER IDX, IJUMP, ISYMS(2), IDXS(2), ISTART, IEND
      INTEGER MFREQ
      INTEGER IDXA, IDXB, IDIP(3)

      DATA SET /.FALSE./

      DATA TABLE /'.OPERAT','.DIPOLE','.SELSTA','.PRINT ','.ALLSTA',
     &            '.HALFFR','.USELEF','.FREQ  ','.FREQUE','.STATIC'/

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
      WARNEXLR1  = .FALSE.
      WARNEXLR2  = .FALSE.

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
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
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
            GOTO (1,2,3,4,5,6,7,8,9,10), IJUMP
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
                    WRITE(LUPRI,'(/2A,I5,/)') 
     &               ' NB. OF OPERATOR PAIRS SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXEXLROP
                    CALL QUIT('TO MANY OPERATOR PAIRS IN CCEXLR.')
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
                WRITE(LUPRI,'(2(/A,I5)/)') 
     &          ' NB. OF OPERATOR PAIRS SPECIFIED  : ',NEXLROPER+9,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXEXLROP
                CALL QUIT('TO MANY OPERATOR PAIRS IN CCEXLR.')
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
                  WRITE(LUPRI,'(/2A,I5)')
     &             ' NB. OF PAIRS OF STATES SPECIFIED',
     &             ' IS GREATER THAN THE ALLOWED NUMBER : ',MXEXLRST
                  WRITE(LUPRI,'(A,2I5,/)') 'IGNORE STATE',ISYMS,IDXS
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
               IF (NEXLRFREQ.EQ.0) THEN
                 NEXLRFREQ  = 1
                 BEXLRFR(1) = ZERO
               ELSE IF (NEXLRFREQ.EQ.1 .AND. BEXLRFR(1).EQ.ZERO) THEN
                 CONTINUE
                 WARNEXLR1 = .TRUE.
               ELSE IF (NEXLRFREQ.GT.1) THEN
                 NEXLRFREQ  = 1
                 BEXLRFR(1) = ZERO
                 WARNEXLR2 = .TRUE.
               END IF
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
                WRITE(LUPRI,'(3(/A,I5)/)') 
     &          'NUMBER OF FREQUENCIES SPECIFIED    :',NEXLRFREQ+MFREQ,
     &          'IS GREATER THAN THE ALLOWED NUMBER :',MXEXLRFR,
     &          'THE NUMBER IS RESET TO THE MAXIMUM :',MXEXLRFR
                MFREQ = MXEXLRFR-NEXLRFREQ
              END IF
              ISTART = NEXLRFREQ+1
              IEND   = NEXLRFREQ+MFREQ
              READ (LUCMD,*) (BEXLRFR(IDX),IDX=ISTART,IEND)
              WRITE (LUPRI,*) 'CC_EXLRINP> ', BEXLRFR(NEXLRFREQ),
     &             NEXLRFREQ
              NEXLRFREQ = NEXLRFREQ+MFREQ
            GO TO 100

C           ---------------------------------------------------
C           .STATIC : add wb = wa = zero to frequency list
C           ---------------------------------------------------
10          CONTINUE
              IF (NEXLRFREQ+1 .GT. MXEXLRFR) THEN
                WRITE(LUPRI,'(3(/A,I5)/)') 
     &          ' NUMBER OF FREQUENCIES SPECIFIED    : ',NEXLRFREQ+1,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXEXLRFR,
     &          ' INPUT OPTION STATIC WILL BE IGNORED.'
              ELSE
                NEXLRFREQ = NEXLRFREQ + 1
                BEXLRFR(NEXLRFREQ) = ZERO
              END IF
            GO TO 100


          ELSE
           WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
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

*---------------------------------------------------------------------*
* check, if any pairs of operator labels specified:
* if not, use default: complete dipole-dipole tensor
*---------------------------------------------------------------------*
      IF (NEXLROPER .EQ. 0) THEN 
        IF (NEXLROPER+9 .GT. MXEXLROP) THEN
          WRITE(LUPRI,'(2(/A,I5)/)') 
     &    ' NB. OF OPERATOR PAIRS SPECIFIED  : ',NEXLROPER+9,
     &    ' IS GREATER THAN THE ALLOWED NUMBER : ',MXEXLROP
          CALL QUIT('TO MANY OPERATOR PAIRS IN CCEXLR.')
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
*---------------------------------------------------------------------*
       SUBROUTINE CC_SMINP(WORD)
*---------------------------------------------------------------------*
*
*    Purpose: read input for CC second moment
*             two photon is a special case
*
*    if (WORD .eq '*CCSM  ') read & process input and set defaults, 
*    else set only defaults 
*
*=====================================================================*
C#if defined (IMPLICIT_NONE)
C      IMPLICIT NONE  
C#else
#  include <implicit.h>
C#endif
#include <priunit.h>
#include <ccsm.h>
#include <ccsminf.h>
#include <ccsdinp.h>
#include <ccsections.h>

* local parameters:
      CHARACTER SECNAM*(8)
      PARAMETER (SECNAM='CC_SMINP')

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
      CHARACTER*8 LABELA, LABELB, LABELC, LABELD
      CHARACTER*70 LABHELP
      CHARACTER TABLE(NTABLE)*(7)

#if defined (SYS_CRAY)
      REAL SMFREQ
#else
      DOUBLE PRECISION  SMFREQ
#endif


      INTEGER IDX, IJUMP
      INTEGER  IDXA, IDXB, IDXC, IDXD, IDIP(3)
      INTEGER  IXSYM , IXST
      DATA SET /.FALSE./

      DATA TABLE /'.OPERAT','.DIPOLE','.PRINT ','.XXXXXX','.XXXXXX',
     &            '.SELSTA','.HALFFR','.XXXXXX','.XXXXXX','.XXXXXX'/

      INTEGER INDPRPCC

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (SET) RETURN
      SET = .TRUE.

*---------------------------------------------------------------------*
* initializations & defaults:
*---------------------------------------------------------------------*
      NSMSEL   = 0

      NSMOPER  = 0

      CCSM = .FALSE.

      IPRSM = 0

      ICHANG = 0
    
      HALFFR = .FALSE.
 
      SELSMST  = .FALSE.
 

*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*
      IF (WORD(1:7) .EQ. '*CCSM  ') THEN

100   CONTINUE
       
* get new input line:
        READ (LUCMD,'(A7)') WORD
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
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
            CALL QUIT('Illegal address in computed GOTO in CC_SMINP.')
      
C           ------------------------------------------
C           .OPERAT : quadruples of operator lables A,B,C,D
C           ------------------------------------------
1           CONTINUE
              READ (LUCMD,'(4A)') LABELA, LABELB, LABELC, LABELD 		
              DO WHILE (LABELA(1:1).NE.'.' .AND. LABELA(1:1).NE.'*')
                IF (LABELA(1:1).NE.'!') THEN
                  IF (NSMOPER.LT.MXSMOP) THEN
                    NSMOPER = NSMOPER + 1
                    IASMOP(NSMOPER) = INDPRPCC(LABELA)
                    IBSMOP(NSMOPER) = INDPRPCC(LABELB)
                    ICSMOP(NSMOPER) = INDPRPCC(LABELC)
                    IDSMOP(NSMOPER) = INDPRPCC(LABELD)
                  ELSE
                    WRITE(LUPRI,'(/2A,I5,/)') 
     &               ' NB. OF OPERATOR QUADRUPLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXSMOP
                    CALL QUIT('TO MANY OPERATOR IN CCSM.')
                  END IF
                END IF
                READ (LUCMD,'(4A)') LABELA, LABELB, LABELC, LABELD
              END DO
              BACKSPACE(LUCMD)
            GO TO 100

C           -------------------------------------------------------
C           .DIPOL : calculate full dipole-dipole-dipole-dipole tensor
C           -------------------------------------------------------
2           CONTINUE
              IF (NSMOPER+81 .GT. MXSMOP) THEN
                WRITE(LUPRI,'(2(/A,I5)/)') 
     &          ' NB. OF OPERATOR QUADRUPLES SPECIFIED  : ',NSMOPER+81,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXSMOP
                CALL QUIT('TO MANY OPERATOR QUADRUPLES IN CCSM.')

              END IF
              IDIP(1) = INDPRPCC('XDIPLEN ')
              IDIP(2) = INDPRPCC('YDIPLEN ')
              IDIP(3) = INDPRPCC('ZDIPLEN ')
              DO IDXA=1,3
              DO IDXB=1,3
              DO IDXC=1,3
              DO IDXD=1,3
                IDX = NSMOPER + (IDXA-1)*27+(IDXB-1)*9+(IDXC-1)*3+IDXD
                IASMOP(IDX) = IDIP(IDXA)
                IBSMOP(IDX) = IDIP(IDXB)
                ICSMOP(IDX) = IDIP(IDXC)
                IDSMOP(IDX) = IDIP(IDXD)
              END DO
              END DO
              END DO
              END DO
              NSMOPER = NSMOPER + 81
            GO TO 100

C           ------------
C           .PRINT 
C           ------------
3           CONTINUE
              READ (LUCMD,*) IPRSM
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
C   .SELSTAtes    Select states and frequences
C                 frequencies are overwritten if .SELHALF is specified
C
6           CONTINUE
              SELSMST =.TRUE. 
              READ (LUCMD,'(A70)') LABHELP
              DO WHILE (LABHELP(1:1).NE.'.' .AND. LABHELP(1:1).NE.'*')
                IF (LABHELP(1:1).NE.'!') THEN
                  READ(LABHELP,*) IXSYM,IXST,SMFREQ
                  IF (NSMSEL.LT.MXSMSEL) THEN
                    NSMSEL = NSMSEL + 1
                    ISMSEL(NSMSEL,1) = IXSYM
                    ISMSEL(NSMSEL,2) = IXST
                    BSMFR(NSMSEL)   = SMFREQ
                  ELSE
                    WRITE(LUPRI,'(/2A,I5,/)')
     &               ' NB. OF STATES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ' ,MXSMSEL 
                    CALL QUIT('TO MANY STATES SPECIFIED BY .SELST')
                  END IF
                END IF
                READ (LUCMD,'(A70)') LABHELP
              END DO
              BACKSPACE(LUCMD)
            GO TO 100
C
C           ------------------------------------------------
C           .HALFFR : impose condition of equal frequencies
C                      for the two lasers 
C           ------------------------------------------------
7           CONTINUE
               HALFFR =.TRUE.
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
           WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
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

*---------------------------------------------------------------------*
* warning if both .SELST AND .SELHLF is specified
* 
      IF (SELSMST.AND.HALFFR) THEN
         WRITE (LUPRI,*)
     &        ' WARNING: BOTH .SELST and .HALFFR are specified'
         WRITE (LUPRI,*) ' .HALFFR is used to obtain frequences'
      END IF
*---------------------------------------------------------------------*
* check, if any quadruples of operator labels specified:
* if not, use default: complete dipole-dipole-dipole-dipole tensor
*---------------------------------------------------------------------*
      IF (NSMOPER .EQ. 0) THEN 
          IF (NSMOPER+81 .GT. MXSMOP) THEN
            WRITE(LUPRI,'(2(/A,I5)/)') 
     &      ' NB. OF OPERATOR QUADRUPLES SPECIFIED  : ',NSMOPER+81,
     &      ' IS GREATER THAN THE ALLOWED NUMBER : ',MXSMOP
            CALL QUIT('TO MANY OPERATOR TRIPLES IN CCQR.')
          END IF
          IDIP(1) = INDPRPCC('XDIPLEN ')
          IDIP(2) = INDPRPCC('YDIPLEN ')
          IDIP(3) = INDPRPCC('ZDIPLEN ')
          DO IDXA=1,3
          DO IDXB=1,3
          DO IDXC=1,3
          DO IDXD=1,3
            IDX = NSMOPER + (IDXA-1)*27+(IDXB-1)*9+(IDXC-1)*3+IDXD
            IASMOP(IDX) = IDIP(IDXA)
            IBSMOP(IDX) = IDIP(IDXB)
            ICSMOP(IDX) = IDIP(IDXC)
            IDSMOP(IDX) = IDIP(IDXC)
          END DO
          END DO
          END DO
          END DO
          NSMOPER = NSMOPER + 81
      END IF

*---------------------------------------------------------------------*
* check, if frequencies are specified; if not, use default: .HALFFR
*---------------------------------------------------------------------*
      IF ( .NOT. SELSMST ) THEN
         IF ( .NOT. HALFFR ) HALFFR = .TRUE.
         WRITE(LUPRI,'(/2A)')
     &      ' *WARNING: NO FREQUENCIES SPECIFIED IN SECOND MOMENT CALC',
     &      ' DEFAULT  .HALFFR  USED '
      END IF
*---------------------------------------------------------------------*
* set CCSM flags:
*---------------------------------------------------------------------*
      WRITE (LUPRI,*) 'SET CCSM TO TRUE'
      CCSM  = .TRUE.

      RETURN
      END
*---------------------------------------------------------------------*

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
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
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
                    WRITE(LUPRI,'(/2A,I5,/)') 
     &               ' NB. OF OPERATOR QUADRUPLES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ',MXTMOP
                    CALL QUIT('TO MANY OPERATOR IN CCTM.')
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
              IF (NTMOPER+243 .GT. MXTMOP) THEN
                WRITE(LUPRI,'(2(/A,I5)/)') 
     &     ' NB. OF OPERATOR QUADRUPLES SPECIFIED  : ',NTMOPER+729,
     &          ' IS GREATER THAN THE ALLOWED NUMBER : ',MXTMOP
                CALL QUIT('TO MANY OPERATOR QUADRUPLES IN CCTM.')

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
                    WRITE(LUPRI,'(/2A,I5,/)')
     &               ' NB. OF STATES SPECIFIED',
     &               ' IS GREATER THAN THE ALLOWED NUMBER : ' ,MXTMSEL 
                    CALL QUIT('TO MANY STATES SPECIFIED BY .SELST')
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
           WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
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
            WRITE(LUPRI,'(2(/A,I5)/)') 
     &      ' NB. OF OPERATOR SIXTUPLES SPECIFIED  : ',NTMOPER+729,
     &      ' IS GREATER THAN THE ALLOWED NUMBER : ',MXTMOP
            CALL QUIT('TO MANY OPERATOR SIXTUPLES IN CCTM.')
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
         WRITE(LUPRI,'(/2A)')
     &      ' *WARNING: NO FREQUENCIES SPECIFIED IN SECOND MOMENT CALC',
     &      ' DEFAULT  .THIRDFr USED '
      END IF
*---------------------------------------------------------------------*
* set CCTM flags:
*---------------------------------------------------------------------*
      WRITE (LUPRI,*) 'SET CCTM TO TRUE'
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
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
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
                      WRITE(LUPRI,'(/2A,I5,/)') 
     &               ' NB. OF OPERATOR TRIPLES SPECIFIED',
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
                WRITE(LUPRI,'(2(/A,I5)/)') 
     &          ' NB. OF OPERATOR TRIPLES SPECIFIED  : ',NMCDOPER+6,
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
                WRITE(LUPRI,'(2(/A,I5)/)')
     &          ' NB. OF OPERATOR TRIPLES SPECIFIED  : ',NMCDOPER+6,
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
              WRITE(LUPRI,'(2(/A))')
     &             'Warning MCD: NO2N+1 not yet carried through',
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
                    WRITE(LUPRI,'(/2A,I5,/)')
     &               ' NB. OF STATES SPECIFIED',
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
            WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
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

*---------------------------------------------------------------------*
* check if any triple of operator labels has been specified:
* if not, use default: complete unrelaxed 
* {dipole x angmom * dipole} tensor
*---------------------------------------------------------------------*
      IF (NMCDOPER .EQ. 0) THEN 
         IF (NMCDOPER+6 .GT. MXMCDOP) THEN
            WRITE(LUPRI,'(2(/A,I5)/)') 
     &      ' NB. OF OPERATOR TRIPLES SPECIFIED  : ',NMCDOPER+6,
     &      ' IS GREATER THAN THE ALLOWED NUMBER : ',MXMCDOP
            CALL QUIT('TO MANY OPERATOR TRIPLES IN CCMCD.')
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
#include <ccslvinf.h>

      CHARACTER SECNAM*(9)
      PARAMETER (SECNAM='CC_SLVINP')
      INTEGER NTABLE
      PARAMETER (NTABLE = 9)
 
      LOGICAL SET
      SAVE SET

      CHARACTER WORD*(7)
      CHARACTER TABLE(NTABLE)*(8)


      DATA SET /.FALSE./
      DATA TABLE /'.SOLVAT','.MXSLIT','.ETOLSL','.TTOLSL','.LTOLSL',
     *            '.PTSOLV','.XXXXXX','.XXXXXX','.XXXXXX'/

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
      ECCCU    = 0.0D0 
      XTNCCCU  = 0.0D0
      XLNCCCU  = 0.0D0

*---------------------------------------------------------------------*
* read input:
*---------------------------------------------------------------------*

      IF (WORD(1:7) .EQ. '*CCSLV ') THEN


100   CONTINUE

        READ (LUCMD,'(A7)') WORD
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
        END DO

        IF (WORD(1:1) .EQ. '.') THEN

          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO

          IF (IJUMP .LE. NTABLE) THEN
            ICHANG = ICHANG + 1
            GOTO (1,2,3,4,5,6,7,8,9), IJUMP
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
            GO TO 100
C
C-----------------------
C
C-----------------------
C
8           CONTINUE
            GO TO 100
C
C-----------------------
C
C-----------------------
C
9           CONTINUE
            GO TO 100

          ELSE
           WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
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
C-------------------------------------------------------------------
C     Finally if we have any solvents  put CCSLV true.
C-------------------------------------------------------------------
C
      CCSLV  = (ICHANG.GT.0)
      IF (CCSLV) RSPIM = .TRUE.
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
      LOGICAL SET
      CHARACTER SECNAM*(8)
      PARAMETER (SECNAM='CC_R12IN')
      PARAMETER (NTABLE = 8)
      CHARACTER WORD*(7)
      CHARACTER TABLE(NTABLE)*(7)
      DATA TABLE /'.NO HYB','.NO A  ','.NO A'' ','.NO B  ',
     &            '.NO RXR','.R12THR','.SVDTHR','.R12XXL'/
      DATA SET/.FALSE./
      IF (SET) RETURN
      SET    = .TRUE.
      IF (WORD(1:4) .EQ. '*R12') THEN
  100   CONTINUE
        READ (LUCMD,'(A7)') WORD
        DO WHILE ( WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#' )
          READ (LUCMD,'(A7)') WORD
        END DO
        IF (WORD(1:1) .EQ. '.') THEN
          IJUMP = 1
          DO WHILE ( IJUMP .LE. NTABLE .AND. TABLE(IJUMP) .NE. WORD)
            IJUMP = IJUMP + 1
          END DO
          IF (IJUMP .LE. NTABLE) THEN
            GOTO (1,2,3,4,5,6,7,8), IJUMP
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
          ELSE
           WRITE (LUPRI,'(/5A,/)') ' Prompt "',WORD,
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
  200 CONTINUE
      RETURN
      END
