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
C  FILE: sirius/sirmp2.F
C  ORIGINAL AUTHORS: Hans Joergen Aa. Jensen and Poul Joergensen (1987).
C  Reference: H. J. Aa. Jensen, P. Jorgensen, H. Agren, and J. Olsen,
C             J. Chem. Phys. 88, 3834 (1988); 89, 5354 (1988)
C
#ifdef OLD_REV_LOG
!===========================================================================
!NOTE (961114/hjaaj): if solvent is implemented in MP2 then
!   IF (DOMP2) FLAG(16) = .FALSE.
!in sirctl.F should be removed.
!/* Comdeck mp2log */
!951130-hjaaj
!removed IADR2 from comdeck INFMP2; moved calculation of IADR2 to new
!  subroutine PHADR2, to be called from RESPONSE.
!941115-hjaaj - PHPACK: rewritten code completely
!   (old version used 10% of CPU time on Cray for SOPPA!)
!941005-hjaaj
!FCKIFR: new routine for "not frozen" index arrays
!MP2FCK: use NRHF instead of IFRMP2 so MP2SET can be called after FCKFRO
!--- Major overhaul of SOPPA code and other code, in particular FCKFRO.
!December 1993-mjp (Martin J. Packer) --- Revised module for SOPPA
!===========================================================================
#endif
C  /* Deck mp2inp */
      SUBROUTINE MP2INP(WORD,INPERR,ALLOPT)
C
C 17/9 - 1990 Hans Agren, Revised Sep 2007-Nov 2009, Hans Joergen Aa. Jensen
C
#include <implicit.h>
      PARAMETER (NTABLE = 10)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7
      LOGICAL ALLOPT
C
C Used from common:
C
C INFORB : NSYM
C INFMP2 : IPRMP2, NFRMP2(*)
C
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infmp2.h>
C
      DATA TABLE /'.MP2 FR','.PRINT ','.NO OV ','.TDASOP','.SC_SRM',
     &            '.SRINTS','.SCSMP2','.SOSMP2','.MP2 SC','.LEVELS'/
CMAERKE 900918: MP2SAV not implemented, should it be?
C
C     ***** PROCESS INPUT FOR MP2INP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108,109,110), II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            GO TO 2000
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT(' ILLEGAL KEYWORD IN MP2INP ')
C
C *** Option 1 >MP2 FROZEN<   Frozen MP2 orbitals
  101 READ(LUCMD,*) (NFRMP2(I),I=1,NSYM)
      GO TO 100
C
C *** Option 2 >PRINT <  General print level in MP2
  102 READ(LUCMD,*) IPRMP2
      GO TO 100
C
C *** Option 3 >NO OV <  Do not add occ-virt mp2 part to density matrix
  103 MP2_NO_OCCVIR = .TRUE.
      GO TO 100
C
C *** Option 4 >TDASOP<  Diagonal 2. order 2p-2h correction in SOPPA
  104 MP2_TDA = .TRUE.
      GO TO 100
C *** Option 5 >SC_SRM<  self-consistent short-range contributions included
#ifdef MOD_SRDFT
  105 SRMP2_SELFCONSISTENT = .TRUE.
      IPRMP2 = MAX(IPRMP2,1)
#else
  105 CONTINUE
      call quit(
     &   '.SC_SRMP2 option not included in this non-srdft version')
#endif 
      GO TO 100
C *** Option 6 >SRINTS<  calculate MP2 energy from AOSR2INT integrals
#ifdef MOD_SRDFT
  106 SRMP2_SRINTS = .TRUE.
      IPRMP2 = MAX(IPRMP2,1)
#else
  106 CONTINUE
      call quit(
     &   '.SRINTS  option not included in this non-srdft version')
#endif 
      GO TO 100
C *** Option 7 >SCSMP2< Grimme's spin-component scaled MP2
C                       p_S = 1.2, p_T = 1/3
  107 MP2_SCS = .TRUE.
      IF (MP2_SCALED) CALL QUIT(
     &   'More than one of .SCSMP2, .SOSPMP2, .MP2 SCALED specified.')
      P_S = 1.2D0
      P_T = (1.0D0/3.0D0)
      MP2_SCALED = .TRUE.
      GO TO 100
C *** Option 8 >SOSMP2< Head-Gordon's scaled opposite spin MP2
C                       p_S = 1.3, p_T = 0
  108 MP2_SOS = .TRUE.
      IF (MP2_SCALED) CALL QUIT(
     &   'More than one of .SCSMP2, .SOSPMP2, .MP2 SCALED specified.')
      P_S = 1.3D0
      P_T = 0.0D0
      MP2_SCALED = .TRUE.
      GO TO 100
C *** Option 9 >MP2 SCALED< Scaled MP2 calculation
  109 READ(LUCMD,*) P_S, P_T
      IF (MP2_SCALED) CALL QUIT(
     &   'More than one of .SCSMP2, .SOSPMP2, .MP2 SCALED specified.')
      MP2_SCALED = .TRUE.
      GO TO 100
C *** Option 10 >LEVELSHIFT< Level-shift of  MP2 denominator
  110 READ(LUCMD,*) MP2_LSHIFT
      GO TO 100
C ==========================================================================
 2000 CONTINUE
      IF (MP2_SCALED) THEN
         MP2_SCALEFAC(1) = P_S + P_T
         MP2_SCALEFAC(2) =     - P_T
         IF (MP2_SCS) THEN
            WRITE(LUPRI,'(A)') " S. Grimme's SCS MP2 requested."
         ELSE IF (MP2_SOS) THEN
            WRITE(LUPRI,'(A)') " M. Head-Gordon's SOS MP2 requested."
         ELSE
            WRITE(LUPRI,'(A,/A,2F14.6)') " User-scaled MP2 requested."
            IF (P_T .EQ. D0) MP2_SOS = .TRUE.
C           ... SOS-MP2 but not Head-Gordon's original SOS-MP2
         END IF
         WRITE(LUPRI,'(A,2F14.6)')
     &      ' Scale factors p_S and p_T:',P_S, P_T
      END IF
      IF (MP2_LSHIFT .NE. 0.0D0) WRITE(LUPRI,'(/A,F14.6)')
     &      ' Level-shift of MP2 denominator: ',MP2_LSHIFT
      IF (MP2_NO_OCCVIR) WRITE(LUPRI,'(/A)')
     &  ' ".NO OV BLOCK": Occ-virt block of MP2 density matrix ignored'
C
      RETURN
      END
C  /* Deck mp2ctl */
      SUBROUTINE MP2CTL(WRK,LFRSAV)
C
C MP2 routines written by Poul Joergensen and Hans Joergen Aa. Jensen
C June 1987.
C Modified Oct 1989 to include frozen orbitals.
C 011293 mjp Correlation coefficients saved to LUSIFC if SOPPA follows.
C
C PURPOSES:
C  1)  CALCULATE ONE ELECTRON DENSITY MATRIX CORRECT THROUGH SECOND ORDER
C  1A) CALCULATE SECOND ORDER ENERGY
C  2)  DETERMINE THE NATURAL ORBITALS
C  3)  SAVE INFORMATION FOR SOPPA
C
#include <implicit.h>
C
      DIMENSION WRK(LFRSAV)
C
      LOGICAL FCKFRO, MP2FRZ
      CHARACTER AO2INTFILE_LABEL_SAVE*8
      PARAMETER ( D0 = 0.0D0 )
C
C Used from common blocks:
C     GNRINF: DOSRIN,SRINTS,CHIVAL
C     EXEINF: FTRCTL
C     INFORB: NORBT,NRHF,?
C     SCBRHF: NFRRHF()
C     INFMP2: NPHSYM, MP2FRO
C     INFINP: FLAG(25)
C     INFIND: ISX, ISMO
C     INFTRA: USEDRC, AO2INTFILE_LABEL
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <gnrinf.h>
#include <exeinf.h>
#include <infopt.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <infmp2.h>
#include <scbrhf.h>
#include <infinp.h>
#include <infind.h>
#include <inftra.h>
#include <dftcom.h>
C
      CALL QENTER('MP2CTL')
      CALL GETTIM(TIMMP2,DUMMY)
      IF (IPRMP2 .GE. 0) WRITE (LUPRI,'(//A//A/A)')
     * ' >>>>> Output from SIRIUS MP2 module <<<<<'
     *,' Reference: H.J.Aa.Jensen, P.Jorgensen, H.Agren, and J.Olsen,'
     *,'            J. Chem. Phys. 88, 3834 (1988); 89, 5354 (1988)'
C
C
      KFRSAV= 1
      KFREE = KFRSAV
      LFREE = LFRSAV
      CALL MEMGET('REAL',KORBEN,NORBT ,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDONE ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KCMO  ,NCMOT ,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KIADR1,MAXVIR*MAXOCC,WRK,KFREE,LFREE)
C
C     READ IN MO'S AND CHECK FOR CANONICAL HARTREE FOCK ORBITALS
C     Use DONE as scratch for FC in MP2FCK.
C
      KFC = KDONE
      CALL MP2FCK(FCKFRO,EMY,WRK(KCMO),WRK(KFC),WRK(KORBEN),
     &            WRK(KFREE),LFREE)
C     CALL MP2FCK(FCKFRO,EMY,CMO,FC,ORBEN,SCRA,LSCRA)
C
C
C     DETERMINE NUMBER OF PARTICLES AND HOLES AND SET UP AN
C     INDEX ARRAY FOR PARTICLES AND HOLES
C
C     Variables MP2FRO, IPHORD, IFRMP2, INDXPH, NP, NH placed in INFMP2
C
C     MP2SET called after MP2FCK because if (FCKFRO) then
C     the number of frozen MP2 orbitals may have been changed.
C
      CALL MP2SET(WRK(KIADR1))
      MP2FRZ = MP2FRO
C     Disable "MP2FRO" if result not to be printed
      IF (IPRMP2 .LE. 0) MP2FRO = .FALSE.
      CALL FLSHFO(LUPRI)
      IF (SRMP2_SRINTS) THEN
         IPRMP2 = MAX(0,IPRMP2)
C        ... make sure MP2 energy with short-range integrals is printed!
      END IF
C
C     if SRMP2_SRINTS 
C     1)  calculate "EMP2" from AO2SRINT integrals (i.e. short-range integrals),
C     2)  set SRMP2_SRINTS false, go back to label 10, and redo the 
C         standard srMP2 "EMP2" from AOTWOINT integrals (i.e. long-range integrals)
      ISRINTS_LOOP = 0
   10 CONTINUE
      IF (SRMP2_SRINTS) THEN
         ISRINTS_LOOP = ISRINTS_LOOP + 1
         IF (ISRINTS_LOOP .EQ. 1) THEN
            WRITE(LUPRI,'(/A)')
     &      ' MP2 energy expression with short-range integrals'
            AO2INTFILE_LABEL_SAVE = AO2INTFILE_LABEL
            AO2INTFILE_LABEL = 'AOSR2INT'
         ELSE IF (ISRINTS_LOOP .EQ. 2) THEN
            WRITE(LUPRI,'(/A)')
     &      ' MP2 energy expression with long-range integrals'
            AO2INTFILE_LABEL_SAVE = AO2INTFILE_LABEL
            AO2INTFILE_LABEL = 'AOTWOINT'
         ELSE IF (ISRINTS_LOOP .EQ. 3) THEN
            WRITE(LUPRI,'(/A)') ' MP2 energy expression with '//
     &      'full (short-range + long-range) integrals'
            AO2INTFILE_LABEL_SAVE = AO2INTFILE_LABEL
            AO2INTFILE_LABEL = 'AOTWOINT'
            CHIVAL = -1.0D0
            DOSRIN = .FALSE.
            SRINTS = .FALSE.
            IF (LUINTA .LE. 0) CALL GPOPEN(LUINTA,AO2INTFILE_LABEL,
     &         'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
            CALL GPCLOSE(LUINTA,'DELETE')
         ELSE
            WRITE(LUPRI,'(/A)')
     &      ' Finished .SRINTS MP2 calculations for srDFT, goodbye!'
            CALL QUIT('No error, finished .SRINTS MP2 calculations')
         END IF
         
         FTRCTL = .TRUE.
C        ... force new sort of AO integrals
      END IF
C
cmjp If response calculation then make space for the correlation coefficients.
      IF (FLAG(25) .AND. .NOT.SRMP2_SRINTS) THEN
         CALL MEMGET('REAL',KKAPPA,NPHTOT(1),WRK,KFREE,LFREE)
         CALL DZERO(WRK(KKAPPA),NPHTOT(1))
         JTRLVL = -10
         USEDRC = .TRUE.
      ELSE
         CALL MEMGET('REAL',KKAPPA,0,WRK,KFREE,LFREE)
         IF (MP2_NO_OCCVIR .OR. SRMP2_SRINTS) THEN
C           2. order integral transformation (ir/js)
            JTRLVL = 4
         ELSE
C           3. order integral transformation (ir/st)
C           (needed for occ-virt block of density matrix)
            JTRLVL = 5
         END IF
      ENDIF
C
C
C     PERFORM INTEGRAL TRANSFORMATION
C     JTRLVL =  4 : code for 2. order transformation for MP2 with MP2_NO_OCCVIR.
C     JTRLVL =  5 : code for 3. order transformation for MP2 density matrix.
C     JTRLVL = 10 : code for full transformation if SOPPA response follows.
C     JTRLVL negative: delete AO integral file after transformation.
C
      CALL TRACTL(JTRLVL,WRK(KCMO),WRK(KFREE),LFREE)
C     CALL TRACTL(ITRLVL,CMO,WRK,LFREE)
C
C
      IF (MP2FRO) THEN
         CALL MEMGET('REAL',KDONEF,N2ORBX,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KCMO1 ,NCMOT ,WRK,KFREE,LFREE)
         CALL MEMGET('INTE',KIFRHF,NORBT ,WRK,KFREE,LFREE)
         CALL FCKIFR(WRK(KIFRHF))
      ELSE
         CALL MEMGET('REAL',KDONEF,0,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KCMO1 ,0,WRK,KFREE,LFREE)
         CALL MEMGET('INTE',KIFRHF,0,WRK,KFREE,LFREE)
      END IF
      CALL MEMGET('REAL',KH2M  ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2D  ,N2ORBX,WRK,KFREE,LFREE)
C
C     DETERMINE NUMBER OF PARTICLE-HOLE DIRAC DISTRIBUTIONS WHICH CAN BE
C     KEPT IN CORE IN ONE LOAD
C
C     We assume max 4*N2ORBX additional work space needed :
      IF (MP2_SOS) THEN
         NPSTP = NP
      ELSE
         NPSTP = (LFREE-4*N2ORBX)/(NP*NH*NH)
      END IF
C     Make load per batch even :
      IF (NPSTP .GT. 0) THEN
         NBATCH = (NP-1)/NPSTP  + 1
         NPSTP  = (NP-1)/NBATCH + 1
         NBATCH = (NP-1)/NPSTP  + 1
      ELSE
         NBATCH = 0
      END IF
      IF ( IPRMP2 .GE. 6 .OR. NPSTP .LE. 0) THEN
         WRITE (LUPRI,'(3(/A,I8))')
     &      ' Number of PARTICLE-HOLE distributions ',NP*NH,
     &      ' Number of BATCHES                     ',NBATCH,
     &      ' Number of P-H distributions each BATCH',NPSTP*NH
       IF ( NPSTP .LE. 0 ) THEN
         WRITE (LUPRI,'(/A)')
     &      ' ERROR: Insufficient work memory for one p-h distribution.'
         CALL ERRWRK('MP2CTL Dirac p-h distributions',(-NP*NH*NH),LFREE)
       END IF
      END IF
      IF (MP2_SOS) THEN
         CALL MEMGET('REAL',KH2PHD,0,WRK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KH2PHD,NPSTP*NH*NH*NP,WRK,KFREE,LFREE)
      END IF
C
      EMP2   = D0
      EMP2FL = D0
C
      CALL DZERO(WRK(KDONE),N2ORBX)
      IF (MP2FRO) THEN
         CALL DCOPY(NCMOT,WRK(KCMO),1,WRK(KCMO1),1)
         CALL DZERO(WRK(KDONEF),N2ORBX)
      END IF
C
      IPLOW  = 1
C
C     ... REPEAT UNTIL IPLOW .GT. NP
C
 100  CONTINUE
         IPHIGH  =  MIN(NP,IPLOW+NPSTP-1)
C
C        Retrieve DIRAC particle-hole integral distributions
C        for this load, if NOT SOS-MP2
C
         IF (.NOT. MP2_SOS)
     &   CALL MP2PHD(WRK(KH2PHD),WRK(KH2M),IPLOW,IPHIGH,WRK,KFREE,LFREE)
C        CALL MP2PHD(H2PHD,H2CD,IPLOW,IPHIGH,WRK,KFREE,LFREE)
C
C        Determine second order contributions to one electron
C        density matrix
C
         CALL MP2MUL(IPLOW,IPHIGH, WRK(KDONE),WRK(KKAPPA),WRK(KDONEF),
     *               WRK(KORBEN),WRK(KH2PHD),WRK(KH2M),WRK(KH2D),
     *               WRK(KIFRHF),EMP2,EMP2FL,WRK(KIADR1),WRK,
     *               KFREE,LFREE)
C        CALL MP2MUL(IPLOW,IPHIGH,DONE,COEMP2,DONEFL,ORBEN,H2PHD,
C    *               H2M,H2D,IFRRHF,EMP2,EMP2FL,IADR1,WRK,KFRSAV,LFRSAV)
         IPLOW   =  IPHIGH + 1
      IF (IPLOW.LE.NP) GO TO 100
C     ... END REPEAT
C
C     Modify DONE to obtain correct one electron density matrix
C     ( factors are defined in MP2DEN ) and add zeroth order
C     contributions
C
      IF ( IPRMP2 .GE. 20 ) THEN
         WRITE(LUPRI,'(/A)')
     *     ' MP2 one-electron density matrix, not modified with factors'
         CALL OUTPUT(WRK(KDONE),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      CALL MP2FAC(WRK(KDONE),WRK(KORBEN))
      IF ( IPRMP2 .GE. 10 ) THEN
         WRITE(LUPRI,'(/A)')
     *     ' Final MP2 one-electron density matrix'
         CALL OUTPUT(WRK(KDONE),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      IF (IPRMP2 .GE. 0) THEN
         IF (MP2_NO_OCCVIR) WRITE(LUPRI,'(/A)')
     &      ' * Occ-virt block of MP2 density matrix ignored.'
         IF (MP2FRZ) THEN
            WRITE (LUPRI,'(/A,8I5)')
     *      ' * MP2 frozen orbitals per symmetry :',(NFRMP2(I),I=1,NSYM)
            WRITE (LUPRI,'(/A)')
     *      " * Total energy and NO's calculated with frozen orbitals."
         END IF
         IF (MP2_SCALED) WRITE (LUPRI,'(/A,2(/A,F14.8))')
     *      ' * MP2 calculation was scaled using the factors',
     *      '     p_S :',MP2_SCALEFAC(1)+MP2_SCALEFAC(2),
     *      '     p_T :',-MP2_SCALEFAC(2)
         IF (MP2_LSHIFT .NE. 0.0D0) WRITE (LUPRI,'(/A,F14.8)')
     *      ' * All 2p-2h orbital energy differences have been '//
     *      'shifted with : ',MP2SHF

         IF (SRMP2_SRINTS) THEN
            if (isrints_loop .eq. 1) then
               WRITE (LUPRI,'(//A,F25.10)')
     *          '   MP2 contribution from short-range integrals :',EMP2
               EMP2_SR = EMP2
            ELSE IF (ISRINTS_LOOP .EQ. 2) THEN
               WRITE (LUPRI,'(//A,F25.10)')
     *          '   MP2 contribution from  long-range integrals :',EMP2
               EMP2_LR = EMP2
            ELSE IF (ISRINTS_LOOP .EQ. 3) THEN
               EMP2_SRLR = EMP2 - EMP2_SR - EMP2_LR
               WRITE (LUPRI,'(//A,F25.10)') '   MP2 contribution from'
     &           //'       sr-lr integrals :',EMP2_SRLR
               EMP2_SRLR = EMP2_SRLR + EMP2_SR
               WRITE (LUPRI,'(//A,F25.10)') '   MP2 contribution from'
     &           //' sr-sr+sr-lr integrals :',EMP2_SRLR
               WRITE (LUPRI,'(//A,F25.10)')
     *          '   MP2 contribution from        full integrals :',EMP2
            ELSE
               call quit('You should never see this output!')
            END IF
         ELSE IF (ADDSRI) THEN
            WRITE (LUPRI,'(/,3(/A,F25.10))')
     *       ' Short-range Hartree-Fock total energy         :',EMY,
     *       ' + MP2 contribution from long-range integrals  :',EMP2,
     *       ' = short-range MP2 second order energy         :',EMY+EMP2
         ELSE IF (.NOT.DODFT) THEN
            IF (MP2_SCALED) THEN
               WRITE (LUPRI,'(/2(/A,F32.10),(/A,F25.10))')
     *       ' Hartree-Fock total energy :',EMY,
     *       ' + Scaled MP2 contribution :',EMP2,
     *       ' = Scaled MP2 second order energy :',EMY+EMP2
            ELSE
               WRITE (LUPRI,'(/,3(/A,F25.10))')
     *       ' Hartree-Fock total energy :',EMY,
     *       ' + MP2 contribution        :',EMP2,
     *       ' = MP2 second order energy :',EMY+EMP2
            END IF
         ELSE ! DFT
            WRITE (LUPRI,'(/,(/A,F25.10),(/A,F25.10,A,F8.5,A,F15.10),
     *        (/A,F25.10))')
     *       ' Kohn-Sham contribution    :',EMY,
     *       ' + scaled PT2 contribution :',EMP2,' *',WDFTMP,
     *       ' = ',EMP2*WDFTMP,
     *       ' = KS-PT2 energy           :',EMY+WDFTMP*EMP2
            EMCSCF = EMY + WDFTMP*EMP2
         END IF
      END IF
C
C
      IF (SRMP2_SRINTS) GO TO 20
C     ... skip to calculate energy with full integrals if
C         calculating MP2 correlation energy with short-range integrals
C
      ECORR = EMY + EMP2
C
C     Write the correlation coefficients and second order one-density to
C     the interface file LUSIFC
C
      IF (FLAG(25)) THEN
         CALL MP2IFC(WRK(KDONE),WRK(KKAPPA))
      END IF
C
C Manu July 07 : For the MP2-SRDFT calculation, we want to write the
C occupation numbers of the MP2 natural orbitals including the
C long-range corrections only (the short-range contributions are
C calculated by MP2SRDFTDEN) for analysis purposes in case SRMP2_SELFCONSISTENT is
C TRUE (which means that we also want to calculate the self-consistency
C effects). Thus we must save CMO before calling natorb for the
C first time otherwise we use in MP2SRDFTDEN natural orbitals to
C calculate the F operator instead of HF-SRDFT orbitals
C (which leads to symmetry dependent results ...)
C
      IF (DOHFSRDFT .AND. SRMP2_SELFCONSISTENT) THEN
         CALL MEMGET('REAL',KCMO2,NCMOT,WRK,KFREE,LFREE)
         CALL DCOPY(NCMOT,WRK(KCMO),1,WRK(KCMO2),1)
      END IF
C
      IF (DOHFSRDFT) WRITE(LUPRI,'(4(/A))')
     &     ' *******************************************************'
     &    ,' MP2-SRDFT natural orbitals: Short-range self-consistent'
     &    ,' contributions are NOT taken into account.'
     &    ,' *******************************************************'
C     
      CALL NATORB(IPRMP2,WRK(KDONE),WRK(KCMO),WRK(KFREE),LFREE)
C
#ifdef MOD_SRDFT
C     ======================================================
C     Manu June 07 : If MP2SRDFT calculation ... 
C
      IF (DOHFSRDFT .AND. SRMP2_SELFCONSISTENT) THEN
        WRITE(LUPRI,'(4(/A))')
     &     ' *****************************************************'
     &    ,' MP2-SRDFT natural orbitals: Short-range contributions'
     &    ,' are now taken into account self-consistently.'
     &    ,' *****************************************************'
C
C        Manu July 07 : we do not want to keep the "long-range MP2 " natural
C        orbitals. Thus CMO is reset to its original value (MO coeffs of
C        the HF-SRDFT orbitals)   
C
         CALL DCOPY(NCMOT,WRK(KCMO2),1,WRK(KCMO),1)
         CALL MEMREL('Releasing mem. used for CMO2',
     &             WRK,KFRSAV,KCMO2,KFREE,LFREE)  
         
C
C        Substract the zeroth order contribution to the density matrix
C
         CALL ADDZEROTH(WRK(KDONE),.FALSE.) 
C             ADDZEROTH(DONE,DOADD) 
C
C        Calculate the second order SRDFT contributions to the density
C        matrix
C
         CALL MP2SRDFTDEN(NFRMP2,WRK(KDONE),WRK(KCMO),
     &                    WRK(KORBEN),WRK(KFREE),LFREE)
C             MP2SRDFTDEN(NFROZ,DONE,CMO,ORBEN,WRK,LFRSAV)
C
C        Add the zeroth order contribution to the density matrix
C
         CALL ADDZEROTH(WRK(KDONE),.TRUE.) 
C
C        Get the SRMP2_SELFCONSISTENT natural orbitals and occ. numbers
C
         CALL NATORB(IPRMP2,WRK(KDONE),WRK(KCMO),WRK(KFREE),LFREE)
C
C
      END IF
C     END IF
C
#endif /*  MOD_SRDFT  */
C     Save MP2 orbitals for MCSCF/CI if not SOPPA
C
      IF (.NOT.FLAG(25)) THEN
         CALL NEWORB('MP2SAVE ',WRK(KCMO),.TRUE.)
      END IF
C
   20 CONTINUE
C
      IF (MP2FRO) THEN
C     If (mp2fro) then also calculate full MP2 energy and NO occupations
C     If iprmp2 .le. 0 then mp2fro will be false
         CALL MP2FAC(WRK(KDONEF),WRK(KORBEN))
         IF ( IPRMP2 .GE. 10 ) THEN
            IF ( FCKFRO ) THEN
               WRITE(LUPRI,'(/A)')
     &       ' MP2 one-electron density matrix with RHF frozen orbitals'
            ELSE
               WRITE(LUPRI,'(/A)')
     &       ' MP2 one-electron density matrix with no frozen orbitals'
            END IF
            CALL OUTPUT(WRK(KDONEF),1,NORBT,1,NORBT,NORBT,NORBT,
     &         1,LUPRI)
         END IF
C
         IF (FCKFRO) THEN
            WRITE (LUPRI,'(/A,8I5//A)')
     &      ' RHF frozen orbitals per symmetry :',(NFRRHF(I),I=1,NSYM),
     &      ' Total MP2 energy and NO occupations calculated with'//
     &      ' frozen orbitals from RHF.'
         ELSE
            WRITE (LUPRI,'(/A)')
     *      ' Full MP2 energy and NO occupations obtained without'//
     *      ' frozen orbitals.'
         END IF
         IF (MP2_SCALED) WRITE (LUPRI,'(/A,2(/A,F15.10))')
     *      ' MP2 calculation was scaled using the factors',
     *      '     p_S :',MP2_SCALEFAC(1)+MP2_SCALEFAC(2),
     *      '     p_T :',-MP2_SCALEFAC(2)
         IF (MP2_LSHIFT .NE. 0.0D0) WRITE (LUPRI,'(/A,F14.8)')
     *      ' * All 2p-2h orbital energy differences have been '//
     *      'shifted with : ',MP2SHF
         IF (SRMP2_SRINTS) THEN
            IF (ISRINTS_LOOP .EQ. 1) THEN
               WRITE (LUPRI,'(//A,F25.10)')
     *         '   MP2 contribution from short-range integrals :',EMP2FL
               EMP2FL_SR = EMP2FL
            ELSE IF (ISRINTS_LOOP .EQ. 2) THEN
               WRITE (LUPRI,'(//A,F25.10)')
     *         '   MP2 contribution from  long-range integrals :',EMP2FL
               EMP2FL_LR = EMP2FL
            ELSE IF (ISRINTS_LOOP .EQ. 3) THEN
               EMP2FL_SRLR = EMP2FL - EMP2FL_SR - EMP2FL_LR
               WRITE (LUPRI,'(//A,F25.10)') '   MP2 contribution from'
     &           //'       sr-lr integrals :',EMP2FL_SRLR
               EMP2FL_SRLR = EMP2FL_SRLR + EMP2FL_SR
               WRITE (LUPRI,'(//A,F25.10)') '   MP2 contribution from'
     &           //' sr-sr+sr-lr integrals :',EMP2FL_SRLR
               WRITE (LUPRI,'(//A,F25.10)')
     *         '   MP2 contribution from        full integrals :',EMP2FL
            ELSE
               call quit('You should never see this output!')
            END IF
C
C           now reset name of 2-electron integral file
            AO2INTFILE_LABEL = AO2INTFILE_LABEL_SAVE
            FTRCTL = .TRUE.
            GO TO 10
C           ... and go back and calculate the MP2 contributions from
C               next class of range-separated integrals /hjaaj Oct 09
         ELSE IF (ADDSRI) THEN
           WRITE (LUPRI,'(/A/3(/A,F25.10))')
     *     ' NOTE that the orbitals saved on disk are those obtained'//
     *     ' with frozen orbitals (see above).',
     *     ' Short-range Hartree-Fock total energy         :',EMY,
     *     ' + MP2 contribution from long-range integrals  :',EMP2FL,
     *     ' = short-range MP2 second order energy         :',EMY+EMP2FL
         ELSE 
            IF (MP2_SCALED) THEN
               WRITE (LUPRI,'(/A/,3(/A,F25.10))')
     *      ' NOTE that the orbitals saved on disk are those obtained'//
     *      ' with frozen orbitals (see above).',
     *      ' Hartree-Fock total energy        :',EMY,
     *      ' + Scaled MP2 contribution        :',EMP2FL,
     *      ' = Scaled MP2 second order energy :',EMY+EMP2FL
            ELSE
              WRITE (LUPRI,'(/A/3(/A,F25.10))')
     *     ' NOTE that the orbitals saved on disk are those obtained'//
     *     ' with frozen orbitals (see above).',
     *     ' Hartree-Fock total energy :',EMY,
     *     ' + MP2 contribution        :',EMP2FL,
     *     ' = MP2 second order energy :',EMY+EMP2FL
            END IF
         END IF
C
C        Manu July 07 : we copy (as done previously) CMO1 into CMO2 since we need
C        the HF-SRDFT orbitals in MP2SRDFTDEN
C
         IF (DOHFSRDFT .AND. SRMP2_SELFCONSISTENT) THEN
            CALL MEMGET('REAL',KCMO2,NCMOT,WRK,KFREE,LFREE)
            CALL DCOPY(NCMOT,WRK(KCMO1),1,WRK(KCMO2),1)
         END IF
C
         IF (DOHFSRDFT) WRITE(LUPRI,'(4(/A))')
     &     ' *******************************************************'
     &    ,' MP2-SRDFT natural orbitals: Short-range self-consistent'
     &    ,' contributions are NOT taken into account.'
     &    ,' *******************************************************'
C     
         CALL NATORB(IPRMP2,WRK(KDONEF),WRK(KCMO1),WRK(KFREE),LFREE)
C
C        Manu June 07 : If MP2SRDFT calculation ... 
C
         IF (DOHFSRDFT .AND. SRMP2_SELFCONSISTENT) THEN
#ifndef MOD_SRDFT
         call quit('srdft not implemented in this version')
#else
            WRITE(LUPRI,'(4(/A))')
     &       ' *****************************************************'
     &      ,' MP2-SRDFT natural orbitals: Short-range contributions'
     &      ,' are now taken into account self-consistently.'
     &      ,' *****************************************************'
C
C           Manu July 07 : we do not want to keep the "long-range MP2 " natural
C           orbitals with frozen orbitals. Thus CMO1 is reset to its original value
C           (MO coeffs of the HF-SRDFT orbitals)   
C
            CALL DCOPY(NCMOT,WRK(KCMO2),1,WRK(KCMO1),1)
            CALL MEMREL('Releasing mem. used for CMO2',
     &                  WRK,KFRSAV,KCMO2,KFREE,LFREE)  
C
C           Substract the zeroth order contribution to the density matrix
C
            CALL ADDZEROTH(WRK(KDONEF),.FALSE.) 
C
C           Calculate the second order SRDFT contributions to the density
C           matrix
C
            CALL MP2SRDFTDEN(NFRRHF,WRK(KDONEF),WRK(KCMO1),
     &                       WRK(KORBEN),WRK(KFREE),LFREE)
C
C           Add the zeroth order contribution to the density matrix
C
            CALL ADDZEROTH(WRK(KDONEF),.TRUE.) 
C
C           Get the natural orbitals and occ. numbers
C
            CALL NATORB(IPRMP2,WRK(KDONEF),WRK(KCMO1),WRK(KFREE),LFREE)
C
#endif
C
         END IF
C
C     Here is the end of IF (MP2FRO .AND. .NOT.FCKFRO)
C
      ELSE IF (SRMP2_SRINTS) THEN
C
C         now reset name of 2-electron integral file
         AO2INTFILE_LABEL = AO2INTFILE_LABEL_SAVE
         FTRCTL = .TRUE.
         GO TO 10
C        ... and go back and calculate the MP2 contributions from
C            long-range integrals /hjaaj Oct 09
      END IF
C     END IF (MP2FRO .AND. .NOT.FCKFRO)
C
      CALL GETTIM(TIMEND,DUMMY)
      TIMMP2 = TIMEND - TIMMP2
      WRITE (LUSTAT,1860) TIMMP2
      IF (IPRMP2 .GE. 2) WRITE (LUPRI,1860) TIMMP2
 1860 FORMAT (/' Time used for MP2 natural orbitals :',F12.3,
     *        ' CPU seconds.')
      CALL MEMREL('MP2CTL',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL FLSHFO(LUSTAT)
      CALL FLSHFO(LUPRI)
      CALL QEXIT('MP2CTL')
      RETURN
      END
C  /* Deck natorb */
      SUBROUTINE NATORB(IPRNO,DONE,CMO,WRK,LFRSAV)
C
C PURPOSE: DETERMINE NATURAL ORBITALS
C
C INPUT  : DONE IS ONE ELECTRON DENSITY MATRIX
C
C
#include <implicit.h>
      DIMENSION DONE(NORBT,NORBT),CMO(*),WRK(*)
C
      PARAMETER ( D0 = 0.0D0 )
C
C Used from common blocks:
C
C INFIND : ISSMO
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
#include <inftap.h>
#include <infpri.h>
C
      CALL QENTER('NATORB')
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LFRSAV
      CALL MEMGET('REAL',KORBNO,N2ORBT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KORBEN,NORBT ,WRK,KFREE,LFREE)
      CALL MEMGET('WORK',KWRK1 ,LWRK1 ,WRK,KFREE,LFREE)
      IJ     = 0
      DO 150 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
         IORBI = IORB(ISYM)
         DO 160 I = 1,NORBI
            DO 170 J = 1,I
               IJ = IJ + 1
               WRK(KORBNO-1+IJ) = DONE(IORBI+J,IORBI+I)
 170        CONTINUE
 160     CONTINUE
 150  CONTINUE
C
      DO 200 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
         IORBI = IORB(ISYM)
      IF (NORBI.EQ.0) GO TO 200
         NBASI = NBAS(ISYM)
         ISSYM = IIORB(ISYM)
         ICSYM = ICMO(ISYM)
C
         CALL JACO(WRK(KORBNO+ISSYM),CMO(1+ICSYM),NORBI,NORBI,
     *             NBASI,WRK(KWRK1),WRK(KWRK1+NORBI))
C        CALL JACO(F,VEC,NB,NMAX,NROWV,BIG,JBIG)
C
         OCCSUM = D0
         DO 175 I=1,NORBI
            II = ISSYM + (I*I+I)/2
            WRK(KORBEN-1+IORBI+I)=WRK(KORBNO-1+II)
            OCCSUM = OCCSUM + WRK(KORBEN-1+IORBI + I)
  175    CONTINUE
#if !defined (VAR_MP2SUPSYM)
         CALL ORDER2(CMO(1+ICSYM),WRK(KORBEN+IORBI),NORBI,NBASI)
#else
C     HJ : implement SUPSYM in MP2 ????? 940823 hjaaj
C          Done, Aug. 2004 /hjaaj
         CALL ORD2SS(CMO(1+ICSYM),WRK(KORBEN+IORBI),
     &               ISSMO(1+IORBI),NORBI,NBASI)
#endif
         IF (IPRNO .GE. 1) THEN
            OCCRHF = 2 * NRHF(ISYM)
            WRITE (LUPRI,1840) ISYM,OCCSUM,OCCRHF,OCCSUM-OCCRHF
#if !defined (VAR_MP2SUPSYM)
            WRITE (LUPRI,1850) (WRK(KORBEN-1+IORBI+I),I=1,NORBI)
#else
            IF (.NOT.SUPSYM) THEN
               WRITE (LUPRI,1850) (WRK(KORBEN-1+IORBI+I),I=1,NORBI)
            ELSE
               WRITE (LUPRI,1851)
     &            (WRK(KORBEN-1+IORBI+I),ISSMO(IORBI+I),I=1,NORBI)
            END IF
#endif
         END IF
  200 CONTINUE
#if defined (VAR_MP2SUPSYM)
      IF (SUPSYM) CALL AVEORD()
C     ... remake ISSORD() as ISSMO() may have changed in ORD2SS
#endif
C
 1840 FORMAT(/' Natural orbital occupation numbers, symmetry',I2
     &       /' Sum =',F15.8,'; RHF =',F15.8,'; Difference =',F15.8)
 1850 FORMAT(/,(4X,5F15.8))
#if defined (VAR_MP2SUPSYM)
 1851 FORMAT(' (supersymmetry in parenthesis)'//,(5(F12.8,' (',I2,')')))
#endif
C
      CALL MEMREL('NATORB',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('NATORB')
      RETURN
C     ... end of NATORB
      END
C  /* Deck mp2den */
      SUBROUTINE MP2DEN(NFROZ,IFROZ,IAH,IMP,IAHMP,ISYMAM,H2PHD,DONE,
     *                  COEMP2,IADR1,ORBEN,H2M,H2D,EMP2,KFRO)
C
C 18-Oct-1989 PJ
C THE PURPOSE OF THIS ROUTINE IS THE SAME AS OLD MP2DEN EXCEPT THAT
C SOME OF THE ORBITALS IN THE MP2 CALCULATIONS MAY BE FROZEN
C AND THEREFORE WILL NOT CONTRIBUTE TO THE MP2 ENERGY AND DENSITY
C MATRIX.  NFROZ(NSYM) specifies the number of frozen orbitals
C per symmetry and IFROZ(NORBT) is an index array telling if
C a specific orbital is frozen or not.
C 011293 mjp KFRO indicates if this is the first or second call
C            to MP2DEN during a frozen calculation.
C
C CALCULATE SECOND ORDER PERTURBATION CONTRIBUTION TO ONE ELECTRON
C DENSITY MATRIX WHICH ORIGINATE FROM <**,MA> AND (**,MA)
C
C  A,B,C,D OCCUPIED AND M,N,P,Q UNOCCUPIED HF ORBITAL INDICES
C
C    D(M,A) =  2/(E(A)-E(M)) *
C                ( SUM(P,Q,C) (CP,MQ)*[QA,CP] -
C                  SUM(P,C,D) (CP,DA)*[MD,CP] )
C
C    D(A,B) = -2 * SUM(P,Q,C) (1/(E(A)+E(C)-E(P)-E(Q))*(AQ,CP)*[QB,CP]
C
C    D(M,N) =  2 * SUM(P,C,D) (1/(E(C)+E(D)-E(P)-E(N))*(DN,CP)*[MD,CP]
C
C      WHERE
C
C   [MA,BN] = (1/(E(A)+E(B)-E(M)-E(N)) * (2*(MA,BN)-(MB,AN))
C           = (1/(E(A)+E(B)-E(M)-E(N)) * (2*(MA,BN)-<MA,BN>)
C
C ( FACTORS BETWEEN = AND * ARE MULTIPLIED ONTO D(*,*) IN MP2FAC )
C
C INPUT:
C  MULLIKEN INTEGRAL DISTRIBUTION (**,MA) = H2M(*,*)
C  DIRAC PARTICLE HOLE INTEGRAL DISTRIBUTION H2PHD(*,*,IAHMP)
C
C OUTPUT:
C  SECOND ORDER CONTRIBUTION ADDED IN DONE(*,*)
C
#include <implicit.h>
C
      DIMENSION NFROZ(8),IFROZ(NORBT),H2PHD(NP,NH,*)
      DIMENSION DONE(NORBT,*),ORBEN(*),H2M(NORBT,*),H2D(NORBT,*)
      DIMENSION COEMP2(*)
C
      PARAMETER ( D0 = 0.0D0 , D2 = 2.0D0 )
C
C Used from common blocks :
C    INFINP : FLAG(25)
C
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
#include <infinp.h>
#include <infmp2.h>
      DIMENSION IADR1(MAXVIR,MAXOCC)
C
      SAVE    TVALUE
      DATA    TVALUE/D0/
C
      CALL QENTER('MP2DEN')
C
C CONSTRUCT   (2*(MA,BN)-(MB,AN))/(E(A)+E(B)-E(M)-E(N)) =
C             (2*(MA,BN)-<MA,BN>)/(E(A)+E(B)-E(M)-E(N))
C
      JMP  = IPHORD(NH+IMP)
      JAH  = IPHORD(IAH)
      ENMA = MP2_LSHIFT - ORBEN(JMP) + ORBEN(JAH)
      IF (IPRMP2 .GE. 30) THEN
         WRITE (LUPRI,'(//A/,2(/A,2I5))')
     *      ' >>> TEST OUTPUT FROM MP2DEN',
     *      ' particle indices JMP, IMP     : ',JMP,IMP,
     *      ' hole     indices JAH, IAH     : ',JAH,IAH
         IF (IPRMP2 .GE. 40) WRITE (LUPRI,'(A,2I5,3(/A,1P,D11.4))')
     *      ' IAHMP index, symmetry(A M)    : ',IAHMP,ISYMAM,
     *      ' EMP2 before this distribution :', EMP2,
     *      ' Orbital energy of MP orbital  :', ORBEN(JMP),
     *      ' Orbital energy of AH orbital  :', ORBEN(JAH)
      END IF
      EMP2BF = EMP2
      DO 100 IH = 1,NH
         IBH = IPHORD(IH)
      IF ( IFROZ(IBH) .NE. 0 ) GO TO 100
         DO 200 IP = 1,NP
            INP = IPHORD(NH+IP)
         IF ( IFROZ(INP) .NE. 0 ) GO TO 200
            ENMANB = ENMA - ORBEN(INP) + ORBEN(IBH)
C    Routine modified to include scaled mp2 calculations
            IF (MP2_SOS) THEN
               H2D(INP,IBH) = MP2_SCALEFAC(1)*H2M(INP,IBH) / ENMANB
            ELSE IF (MP2_SCALED) THEN
               H2D(INP,IBH) = ( MP2_SCALEFAC(1)*H2M(INP,IBH)
     *                        + MP2_SCALEFAC(2)*H2PHD(IP,IH,IAHMP))
     *                        / ENMANB
            ELSE
               H2D(INP,IBH) = ( D2*H2M(INP,IBH) - H2PHD(IP,IH,IAHMP))
     *                        / ENMANB
            END IF
            IF (IPRMP2 .GE. 30) THEN
               TVALUE = TVALUE + H2D(INP,IBH)*H2D(INP,IBH)
               IF (IPRMP2 .GE. 50 .AND.
     *             H2D(INP,IBH).NE.D0 .AND. H2M(INP,IBH).NE.D0) THEN
                  WRITE (LUPRI,'(/A,4I5,/A,1P,3D15.6)')
     *            ' IH,IBH,IP,INP         :',IH,IBH,IP,INP,
     *            ' ENMANB, a(MA,NB), H2M :',
     *            ENMANB,H2D(INP,IBH),H2M(INP,IBH)
               END IF
            END IF
C
C Accumulate MP2 energy :
C
            EMP2 = EMP2 + H2M(INP,IBH) * H2D(INP,IBH)
C
C PARTICLE HOLE PART OF MULLIKEN DISTRIBUTIONS ARE DIVIDED WITH ORBITAL
C ENERGY DIFFERENCES TO BE ABLE TO CARRY OUT MATRIX MULTIPLICATIONS.
C (THE MODIFIED MATRIX ELEMENTS ENTER ONLY IN DOUBLES CONTRIBUTIONS,
C  NOT IN SINGLES)
C
            H2M(INP,IBH) = H2M(INP,IBH)/ENMANB
            H2M(IBH,INP) = H2M(INP,IBH)
 200     CONTINUE
 100  CONTINUE
cmjp Pack the correlation coefficients into COEMP2 if SOPPA follow
      IF (FLAG(25) .AND. KFRO .GE. 0) THEN
         CALL PHPACK(H2D,COEMP2(IADR1(IMP,IAH)+1),ISYMAM,1,-1)
      ENDIF
C
C
C  [*,*] = H2D(*,*)
C  (*,*) = H2M(*,*)
C
C ADD CONTRIBUTION TO DONE (FACTORS ARE TAKEN CARE OF LATER)
C
C   DONE(A,M) = (A,D)*[M,D]
C   DONE(M,N) = (N,D)*[M,D]
C
C TODO: if (MP2_NO_OCCVIR) then occ-virt in DONE not needed in DGEMM calls below /hjaaj Apr09
C
      DO 1000 ISSYM = 1,NSYM

         NORBI = NORB(ISSYM)
      IF (NORBI .EQ. 0) GOTO 1000
         IOFFI = IORB(ISSYM) + 1
         NFRMPI= NFROZ(ISSYM)
         NRHFI = MAX(NRHF(ISSYM), NFRMPI)
C        ... This special code makes it possible to exclude low-lying
C            virtual MO's in the MP2 energy /hjaaj Oct 09
         NVIRI = NORBI - NRHFI

         JSYM  = MULD2H(ISYMAM,ISSYM)
         NORBJ = NORB(JSYM)
      IF (NORBJ .EQ. 0) GOTO 1000
         IOFFJ = IORB(JSYM) + 1
         NFRMPJ= NFROZ(JSYM)
         NRHFJ = MAX(NRHF(JSYM),NFRMPJ)
C        ... See comment above for NRHFI
         NVIRJ = NORBJ - NRHFJ
C
C   DONE(A,B) = (A,Q)*[Q,B]
C   DONE(M,B) = (M,Q)*[Q,B]
C        I,I     I,J   J,I
C
         IF (NRHFI .GT. NFRMPI)
     &   CALL DGEMM('N','N',NORBI-NFRMPI,NRHFI-NFRMPI,NVIRJ,1.D0,
     &              H2M(IOFFI+NFRMPI,IOFFJ+NRHFJ),  NORBT,
     &              H2D(IOFFJ+NRHFJ,IOFFI+NFRMPI),  NORBT,1.D0,
     &              DONE(IOFFI+NFRMPI,IOFFI+NFRMPI),NORBT)
C
C   DONE(A,M) = (A,D)*[M,D]
C   DONE(M,N) = (N,D)*[M,D]
C        I,I     I,J   I,J
C
         IF (NRHFJ .GT. NFRMPJ)
     &   CALL DGEMM('N','T',NORBI-NFRMPI,NVIRI,NRHFJ-NFRMPJ,1.D0,
     &              H2M(IOFFI+NFRMPI,IOFFJ+NFRMPJ),NORBT,
     &              H2D(IOFFI+NRHFI,IOFFJ+NFRMPJ), NORBT,1.D0,
     &              DONE(IOFFI+NFRMPI,IOFFI+NRHFI),NORBT)
 1000 CONTINUE
      IF (IPRMP2 .GE. 30) THEN
         WRITE (LUPRI,'(1P,2(/A,2D11.4))')
     *   ' del(EMP2), EMP2 after this distribution :',EMP2-EMP2BF,EMP2,
     *   ' T = sum (a(MA,NB)**2) after this distr. :',TVALUE
      END IF
      CALL QEXIT('MP2DEN')
      RETURN
C     ... return from MP2DEN
      END
C  /* Deck mp2fac */
      SUBROUTINE MP2FAC(DONE,ORBEN)
C
C Modify DONE to obtain correct one electron density matrix
C ( Factors are defined in MP2DEN )
C
#include <implicit.h>
C
      DIMENSION DONE(NORBT,*),ORBEN(*)
C
      PARAMETER (D0 = 0.0D0 , D2 = 2.0D0 ,TOLELE = 1.0D-4 )
C
#include <priunit.h>
#include <maxorb.h>
#include <inforb.h>
#include <infmp2.h>
C
      CALL QENTER('MP2FAC')
C
C     Scale with factor 2 given in formulas in comments in MP2DEN
C
      CALL DSCAL(N2ORBX,D2,DONE,1)
      ELERMV = D0
      ELERCV = D0
      DO 1500 ISYM = 1,NSYM
         IORBI = IORB(ISYM)
         NRHFI = NRHF(ISYM)
         NORBI = NORB(ISYM)
         DO 1600 I = (IORBI+1),(IORBI+NRHFI)
            DO J = (IORBI+1),(IORBI+NRHFI)
               DONE(J,I) = -DONE(J,I)
            END DO
            IF (MP2_NO_OCCVIR) THEN
               DO J = (IORBI+1+NRHFI),(IORBI+NORBI)
                  DONE(J,I) = D0
                  DONE(I,J) = D0
               END DO
            ELSE
               DO J = (IORBI+1+NRHFI),(IORBI+NORBI)
                  DONE(J,I) = (DONE(J,I)-DONE(I,J))/(ORBEN(I)-ORBEN(J))
                  DONE(I,J) = DONE(J,I)
               END DO
            END IF
            ELERMV = ELERMV - DONE(I,I)
C
C           Add zero'th order contribution
C
            DONE(I,I) = D2 + DONE(I,I)
 1600    CONTINUE
         DO 1900 IUNOCC = (IORBI+1+NRHFI),(IORBI+NORBI)
            ELERCV = ELERCV + DONE(IUNOCC,IUNOCC)
 1900    CONTINUE
 1500 CONTINUE
C
      IF ( ABS(-ELERMV+ELERCV).GT.TOLELE ) THEN
         NWARN = NWARN + 1
         WRITE(LUPRI,'(/A/A,F12.6,A)') ' *******WARNING*********'//
     *   ' MP2 CALCULATION CARRIED OUT INCORRECTLY ',
     *   ' TRACE OF ONE ELECTRON DENSITY MATRIX GIVE ',
     *   (-ELERMV+ELERCV),' ELECTRONS IN ADDITION TO SPECIFIED NUMBER'
         I = 1
      ELSE
         I = 0
      END IF
      IF (ELERMV .LT. 0 .OR. ELERMV .GT. 2*NRHFT) THEN
         NWARN = NWARN + 1
         WRITE(LUPRI,'(/A/A,I5/A,F12.6/A/A)')
     &   ' *******WARNING from MP2*********',
     &   ' Number of electrons in Hartree-Fock orbitals :',2*NRHFT,
     &   ' Number of electrons removed from HF orbitals :',ELERMV,
     &   ' Probable error: this molecule is not suitable for MP2',
     &   ' because E(LUMO) - E(HOMO) small or negative.'
         I = 1
      ELSE
         I = 0
      END IF
      IF (IPRMP2 .GE. 2 .OR. I .EQ. 1) THEN
         WRITE(LUPRI,'(/A,F11.6,A)')
     *      ' MP2 move',ELERCV,' electrons to unoccupied HF orbitals'
         IF (IPRMP2 .GE. 6 .OR. I .EQ. 1)
     &      WRITE(LUPRI,'(/A,F11.6,A/A,1P,D10.2,A)')
     *      ' MP2 move',ELERMV,' electrons from occupied HF orbitals',
     *      ' MP2 change in trace of one electron density matrix :',
     *      -ELERMV+ELERCV,' electrons.'
      END IF
C
      CALL QEXIT('MP2FAC')
      RETURN
C
C     ... END OF MP2FAC
C
      END
C  /* Deck mp2fck */
      SUBROUTINE MP2FCK(FCKFRO,EMY,CMO,FC,ORBEN,SCRA,LSCRA)
C
C Purpose:
C  CONSTRUCT INACTIVE FOCK MATRIX AND CHECK IF ORBITALS ARE
C  CANONICAL HARTREE-FOCK ORBITALS
C
C Input:
C  CMO; initial molecular orbitals used to build Fock matrix,
C       assumed to be orthonormal.
C
C Output:
C  EMY; Fock energy
C  CMO; molecular orbitals diagonalizing Fock matrix
C
C Scratch:
C  FC; the inactive Fock matrix and scratch area for overlap matrix
C  SCRA; general scratch area
C
#include <implicit.h>
      DIMENSION CMO(*),FC(*),ORBEN(*),SCRA(*)
      LOGICAL   FCKFRO, FOUND
C
C
      PARAMETER (VALCON_MAX = 1.D-4, VALCON_WARN = 1.D-5,
     &           D0 = 0.0D0, DBIG = 1.D+12, GAPMIN = 0.1D0)
#include <dummy.h>
C
C  Used from common blocks:
C     GNRINF : EMBEDDING
C     INFINP : POTNUC
C     INFORB : NRHF(*), NRHFT, ...
C     SCBRHF : IOPRHF,NFRRHF(*)
C     INFMP2 : NFRMP2(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <gnrinf.h>
#include <infinp.h>
#include <inforb.h>
#include <scbrhf.h>
#include <infmp2.h>
#include <infind.h>
#include <inftap.h>
#include <infpri.h>
C
      IJ_FOCK(I,J) = (J-1)*NORBI + I
C
      CALL QENTER('MP2FCK')

      IF (IOPRHF .GT. 0 .OR. NASHT .GT. 0) THEN
         WRITE (LUPRI,'(//A)')
     &   ' FATAL ERROR in MP2 module: MP2 is not implemented for'//
     &   ' open-shell Hartree-Fock'
         CALL QUIT('ERROR: MP2 is not implemented for open-shell RHF')
      END IF
C
      IF (IPRMP2 .GE. 3) WRITE (LUPRI,'(/A//A,I5/A,8I5)')
     &   ' Check that orbitals are canonical HARTREE-FOCK orbitals',
     &   '    Number of electrons :',2*NRHFT,
     &   '    Orbital occupations :',(NRHF(I),I=1,NSYM)
C
      IF (LSCRA .LT. N2BASX) CALL ERRWRK('MP2FCK',-N2BASX,LSCRA)
C
      CALL READMO(CMO,9,SCRA,LSCRA)

C
C     Step 1 : Construct or read the inactive Fock matrix
C
      IF ( EMBEDDING ) THEN ! read EMY and Fock matrix from SIRFC
         WRITE (LUPRI,'(/A)')
     &   ' MP2 for embedded QM system, reading Fock matrix from SIRIFC'
         KCMO      = 1
         KCMO_diff = KCMO  + NCMOT
         KWRK      = KCMO_diff  + NCMOT
         LWRK      = LSCRA - KWRK
         CALL RD_SIRIFC('CMO',FOUND,SCRA(KCMO),SCRA(KWRK),LWRK)
         CALL DCOPY(NCMOT,CMO,1,SCRA(KCMO_diff),1)
         CALL DAXPY(NCMOT,-1.0D0,CMO,1,SCRA(KCMO_diff),1)
         KCMO_diff_MAX = KCMO_diff - 1 + IDAMAX(NCMOT,SCRA(KCMO_diff),1)
         IF ( abs(SCRA(KCMO_diff_MAX)) .GT. 1.0D-8) THEN
            WRITE(LUPRI,'(//A/A)')
     &      'FATAL ERROR: CMO on SIRIFC not same as CMO from SIRIUS.RST'
     &     ,'             The difference matrix is :'
            !CALL PRORB(SCRA(KCMO_diff),.FALSE.,LUPRI)
            ! problem: PRORB only prints coefficients .gt. 0.0001, thus
            ! call output instead
            DO ISYM = 1,NSYM
               WRITE(LUPRI,'(/A,I5)')
     &         ' CMO(SIRIFC) - CMO(SIRIUS.RST), symmetry',ISYM
               JCMO_diff = KCMO_diff + ICMO(ISYM)
               CALL OUTPUT(SCRA(JCMO_diff),1,NBAS(ISYM),1,NORB(ISYM),
     &            NBAS(ISYM),NORB(ISYM),-1,LUPRI)
            END DO
            WRITE(LUPRI,'(//A/A)')
     &      'FATAL ERROR: CMO on SIRIFC not same as CMO from SIRIUS.RST'
     &     ,'             CMO matrix from SIRIFC :'
            CALL PRORB(SCRA(KCMO),.FALSE.,LUPRI)
            WRITE(LUPRI,'(//A/A)')
     &      'FATAL ERROR: CMO on SIRIFC not same as CMO from SIRIUS.RST'
     &     ,'             CMO matrix from SIRIUS.RST :'
            CALL PRORB(CMO       ,.FALSE.,LUPRI)
            CALL QUIT('CMO on SIRIFC not same as CMO from SIRIUS.RST')
         END IF
         CALL RD_SIRIFC('FC for MP2',FOUND,FC,SCRA,LSCRA)
         IF (.NOT. FOUND) THEN
            CALL QUIT('Info for "FC for MP2" not found on SIRIFC')
         END IF
      ELSE ! construct Fock matrix
         WRITE (LUPRI,'(/A)') 'Generating Fock matrix'
         EMY = D0
         CALL FCKMAT(.TRUE.,DUMMY,CMO,EMY,FC,DUMMY,SCRA,LSCRA)
C        CALL FCKMAT(ONLYFC,DV,CMO,EMY,FC,FV,WRK,LFREE)
      END IF
 
      IF (IPRMP2 .GE. 3) WRITE (LUPRI,1730) EMY
      EMY = EMY + POTNUC
      IF (IPRMP2 .GE. 3) WRITE (LUPRI,1732) EMY
C
C     Step 2 : Check if converged (using that gradient is
C              proportional to occ-unocc Fock matrix elements).
C              Check if NFRMP2(i) .ge. NFRRHF(i).
C
      FCKFRO = .FALSE.
      DO 299 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
      IF (NORBI.EQ.0) GO TO 299
         NFRRHF_I = NFRRHF(ISYM)
         IF (NFRRHF_I  .GT. 0) FCKFRO = .TRUE.
         IF (NFRRHF_I  .GT. NFRMP2(ISYM)) THEN
            NWARN = NWARN + 1
            WRITE (LUPRI,'(3(/A,I5))')
     +      ' ***** WARNING from MP2FCK for symmetry',ISYM,
     +      ' Number of frozen MP2 orbitals this symmetry',NFRMP2(ISYM),
     +      ' is reset to number of frozen RHF orbitals  ',NFRRHF_I 
            NFRMP2(ISYM) = NFRRHF_I 
         ELSE IF (NFRMP2(ISYM) .GT. NRHF(ISYM)) THEN
            WRITE (LUPRI,'(3(/A,I5),/A/A)')
     +      ' INFO from MP2FCK for symmetry',ISYM,
     +      ' INFO Number of frozen MP2 orbitals         ',NFRMP2(ISYM),
     +      ' INFO is greater than number of RHF orbitals',NRHF(ISYM),
     +      ' INFO The thus frozen virtual orbitals are excluded.',
     +      ' INFO We hope you know what you are doing!!!'
         END IF
         NRHFI = NRHF(ISYM)
         IORBI = IORB(ISYM)
         ISSYM = IIORB(ISYM)
         CALL DSPTGE(NORBI,FC(ISSYM+1),SCRA)
         VALCON  = D0
         IVALCON = 0
         DO I = NFRRHF_I+1,NRHFI
C
C           To check if this occupied orbital is converged we find
C           the 2-norm of the coupling to the unoccupied orbitals
            I1 = IJ_FOCK(NRHFI+1,I)
            NORBI_MP2 = NORBI - NRHFI
            VALCONI = DNRM2(NORBI_MP2,SCRA(I1),1)
            IF (VALCONI .GT. VALCON) THEN
               VALCON  = VALCONI
               IVALCON = I
            END IF
         END DO
C        ... VALCON now max euclidean 2-norm of off-diagonal elements for any row
         IF ( VALCON .GT. VALCON_MAX ) THEN
            WRITE(LUPRI,'(/A//A,F12.8,A,I6,A,I3//A)') 'ERROR: '//
     +  ' Molecular orbitals are NOT optimized SCF orbitals',
     +  ' MAX 2-norm of off-diagonal Fock matrix elements is',VALCON,
     +  ' for MO no.',IVALCON,' in symmetry',ISYM,
     +  ' This Fock matrix row is:'
            CALL OUTPUT(SCRA(IJ_FOCK(IVALCON,1)), 1,1, NRHFI+1,NORBI,
     +                  NORBI,NORBI, 1,LUPRI)
            CALL QUIT('MP2FCK: MOs are NOT optimized SCF orbitals')
         ELSE IF ( VALCON .GT. VALCON_WARN ) THEN
            WRITE(LUPRI,'(/A//A,F12.8,A,I6,A,I3//A)') 'WARNING: '//
     +  ' Molecular orbitals are not fully optimized SCF orbitals',
     +  ' MAX 2-norm of off-diagonal Fock matrix elements is',VALCON,
     +  ' for MO no.',IVALCON,' in symmetry',ISYM,
     +  ' This Fock matrix row is:'
            CALL OUTPUT(SCRA(IJ_FOCK(IVALCON,1)), 1,1, NRHFI+1,NORBI,
     +                  NORBI,NORBI, 1,LUPRI)
         ELSE IF (IPRMP2 .GE. 4) THEN
            WRITE(LUPRI,'(/A,I3,A,1P,D10.2)') 'Symmetry',ISYM,
     +  ' MAX 2-norm of a row of off-diagonal Fock matrix elements is',
     +      VALCON
         END IF
  299 CONTINUE
C
C     Step 3 : Make sure orbitals are canonical orbitals
C              (based on code from sircan.F:SIRCNO /hjaaj oct 09)
C              Extract orbital energies.
C
      CALL DZERO(ORBEN,NORBT)
      EHOMO  = -DBIG
      ELUMO  =  DBIG
      DO 399 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
      IF (NORBI.EQ.0) GO TO 399
         NFRRHF_I = NFRRHF(ISYM)
         NRHFI = NRHF(ISYM)
         IORBI = IORB(ISYM)
         NBASI = NBAS(ISYM)
         JFC   = 1 + IIORB(ISYM)
         JCMO  = 1 + ICMO(ISYM)
C
C        Zero couplings between inactive-secondary blocks
C        in Fock matrix and between frozen orbitals and all other.
C
         IF (NFRRHF_I .GT. 0) THEN
C           zero off-diagonal frozen-frozen elements
            DO I = 2,NFRRHF_I
               JSTA = JFC + IROW(I)
               JEND = JSTA + I - 2
               DO J = JSTA,JEND
                  SCRA(J) = D0
               END DO
            END DO
         END IF
         IF (NRHFI .GT. 0) THEN
C           zero inactive-secondary block
C           Note that any frozen orbitals are also "inactive".
            CALL DSPZERO(NORBI,FC(JFC),1,NRHFI,NRHFI+1,NORBI)
C           call dspzero(n,asp, nrsta,nrend, ncsta,ncend)
         END IF
C
C        generate canonical orbitals
         KBIG  = 1
         KJBIG = KBIG + NORBI
         CALL JACO(FC(JFC),CMO(JCMO),NORBI,NORBI,NBASI,
     *             SCRA(KBIG),SCRA(KJBIG))
C        CALL JACO (F,VEC,NB,NMAX,NROWV,BIG,JBIG)
C
         DO I = 1,NORBI
            ORBEN_I = FC(JFC-1 + IROW(I+1))
            ORBEN(IORBI + I) = ORBEN_I
            IF (I .LE. NRHFI) THEN
               IF (EHOMO .LT. ORBEN_I) THEN
                  EHOMO = ORBEN_I
                  ISYM_HOMO = ISYM
               END IF
            ELSE
               IF (ELUMO .GT. ORBEN_I) THEN
                  ELUMO = ORBEN_I
                  ISYM_LUMO = ISYM
               END IF
            END IF
         END DO

         IF (IPRMP2 .GE. 3) THEN
            WRITE (LUPRI,1740) ISYM, NRHFI
            IF (NFRRHF_I .GT. 0) WRITE (LUPRI,'(A,I3,A)')
     +        ' Note: the first',NFRRHF_I ,' orbitals are frozen in SCF'
            WRITE (LUPRI,1750) (ORBEN(IORBI+I),I=1,NORBI)
         END IF
C
  399 CONTINUE
C
      IF (ELUMO-EHOMO .LT. GAPMIN) THEN
         NWARN = NWARN + 1
         WRITE(LUPRI,'(/A/A)')
     &   ' *******WARNING from MP2*********'//
     &   ' This molecule is probably not suitable for MP2',
     &   ' *******WARNING from MP2*********'//
     &   ' because E(LUMO) - E(HOMO) small or negative.'
         I = 1
      ELSE
         I = 0
      END IF
      IF (IPRMP2 .GE. 1 .OR. I .EQ. 1) THEN
         WRITE (LUPRI,'(2(/A,F15.8,A,I2,A)/A/A,F15.8)')
     &      '    E(LUMO) :',ELUMO,'   (in symmetry',ISYM_LUMO,')',
     &      '  - E(HOMO) :',EHOMO,'   (in symmetry',ISYM_HOMO,')',
     &      '  --------------------------',
     &      '    gap     :',ELUMO-EHOMO
      END IF
      CALL QEXIT('MP2FCK')
      RETURN
C
 1730 FORMAT(/' Hartree-Fock electronic energy:',F25.12)
 1732 FORMAT( ' Hartree-Fock total      energy:',F25.12)
 1740 FORMAT(/' Hartree-Fock orbital energies, symmetry',I2,
     &         ', number of occupied SCF orbitals is',I5/)
 1750 FORMAT(4X,5F15.8)
C
C
C *** end of subroutine MP2FCK
C
      END
C  /* Deck MP2IFC */
        SUBROUTINE MP2IFC(DONE,COEMP2)
C
C This routine prepares the SOPPA interface for RESPONSE storing
C the correlation coefficients and the second order one-density
C matrix on the end of the LUSIFC file.
C    MJP 181093
C
#include <implicit.h>
#include <priunit.h>
        DIMENSION DONE(*), COEMP2(*)
C Used from common blocks:
C    INFORB: N2ORBT, NRHFT, NVIRT
C    INFMP2: NPHMAX
C    INFTAP: LUSIFC
C
#include <maxorb.h>
#include <dummy.h>
#include <inforb.h>
#include <infmp2.h>
#include <inftap.h>
      CHARACTER*8 LAB123(3), LABEOD, LABELN, LMP2INFO
      DATA LAB123/'********','********','********'/
      DATA LABEOD /'EODATA  '/
      DATA LMP2INFO /'MP2INFO '/
      CALL GETDAT(LAB123(2),LAB123(3))
C     ... place date in LAB123(2) and time in LAB123(3)
C
C Open the FNSIFC file
      IF (LUSIFC .LT. 0) CALL GPOPEN(LUSIFC,FNSIFC,
     &   'UNKNOWN',' ','UNFORMATTED',IDUMMY,.FALSE.)
C Append MP2->SOPPA info at the end of the file
      REWIND (LUSIFC)
      CALL MOLLAB(LABEOD,LUSIFC,LUERR)
      BACKSPACE (LUSIFC)
C Write "MP2INFO " label
      WRITE (LUSIFC) LAB123,LMP2INFO
C Write the correlation coefficients
      CALL WRITT (LUSIFC,NPHTOT(1),COEMP2)
C Write the second'order one-electron density
      CALL WRITT (LUSIFC,N2ORBX,DONE)
C Write "End Of DATA" label
      WRITE (LUSIFC) LAB123,LABEOD
C Close the FNSIFC file
      CALL GPCLOSE(LUSIFC,'KEEP')
      RETURN
      END
C  /* Deck mp2mul */
      SUBROUTINE MP2MUL(IPLOW,IPHIGH,DONE,COEMP2,DONEFL,ORBEN,H2PHD,
     *                  H2M,H2D,IFRRHF,EMP2,EMP2FL,IADR1,
     *                  WRK,KFRSAV,LFRSAV)
C
C WRITTEN APRIL-1-87
C   Revision 891017-hjaaj: use NXTH2M
C 011293 mjp call to MP2DEN altered since IPHPORD etc. now in /INFMP2/
C
C PURPOSE:
C     READ IN INTEGRAL DISTRIBUTIONS TO BE USED TO CALCULATE
C     SECOND ORDER CONTRIBUTIONS TO THE ONE ELECTRON DENSITY MATRIX
C
C Scratch:
C     WRK(kfrsav:kfrsav-1+lfrsav)
C
#include <implicit.h>
C
      DIMENSION DONEFL(NORBT,NORBT),DONE(NORBT,NORBT),ORBEN(NORBT)
      DIMENSION H2PHD(NP,NH,*),H2M(NORBT,NORBT),H2D(NORBT,NORBT)
      DIMENSION WRK(*), COEMP2(*), IFRRHF(*)
C
C -- local constants
C
      PARAMETER ( D2 = 2.0D0 )
C
C Used from common blocks:
C   INFORB : N2ORBX, ...
C   SCBRHF : NFRRHF(*)
C   INFMP2 : IPHORD, NP, NH, IFRMP2, MP2FRO
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
#include <scbrhf.h>
#include <infmp2.h>
      DIMENSION IADR1(MAXVIR,MAXOCC)
C     
      LOGICAL CHOLE,DHOLE, L_MP2_TDA

      DIMENSION NEEDTP(6)
      DATA NEEDTP /0,1,1,1,1,0/
C     ... only particle-hole distributions needed
C         assumed: NRHF(isym) .ge. NISH(isym)
C              and NVIR(isym) .le. NSSH(isym)
C
      CALL QENTER('MP2MUL')
      KFREE = KFRSAV
      LFREE = LFRSAV
      IF (MP2FRO) THEN
         CALL MEMGET('REAL',KH2M1 ,N2ORBX,WRK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KH2M1 ,0,WRK,KFREE,LFREE)
      END IF
      CALL DZERO(H2D,N2ORBX)
      L_MP2_TDA = MP2_TDA .AND. IPLOW .EQ. 1
      IF (L_MP2_TDA) THEN
         NEEDTP(1) = 1
         NEEDTP(6) = 1
         CALL MEMGET('REAL',KJ,N2ORBX,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KK,N2ORBX,WRK,KFREE,LFREE)
         CALL DZERO(WRK(KJ),N2ORBX)
         CALL DZERO(WRK(KK),N2ORBX)
!        write(0,*) 'hjaaj MP2MUL: allocated RJ and RK, KFREE',KFREE
      ELSE
         NEEDTP(1) = 0
         NEEDTP(6) = 0
      END IF
C
C ****************************************************************
C     Loop over Mulliken distributions allowed in NEEDTP(6)
C
      IDIST = 0
  100 CALL NXTH2M(IC,ID,H2M,NEEDTP,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 800
C     ... if idist .lt. 0 then no more distributions
C
         CHOLE = ( INDXPH(IC) .GT. 0 )
         DHOLE = ( INDXPH(ID) .GT. 0 )
C
C        Check if we need this distribution:
C        HOLE-HOLE AND PARTICLE-PARTICLE  DISTRIBUTIONS NOT NEEDED
C
      IF (L_MP2_TDA) THEN
         IF (IC.EQ.ID) THEN
            DO IA = 1,NORBT
               WRK(KJ+(IC-1)*NORBT+(IA-1)) = H2M(IA,IA)
               WRK(KJ+(IA-1)*NORBT+(IC-1)) = H2M(IA,IA)
            END DO
            WRK(KK+(IC-1)*NORBT+(IC-1)) = H2M(IC,IC)
         ELSE
            WRK(KK+(IC-1)*NORBT+(ID-1)) = H2M(IC,ID)
            WRK(KK+(ID-1)*NORBT+(IC-1)) = H2M(IC,ID)
         END IF
      END IF
      IF ( CHOLE .EQV. DHOLE ) GO TO 100
C
         IF (DHOLE) THEN
            ISWAP = IC
            IC    = ID
            ID    = ISWAP
         END IF
C
C        IC IS A HOLE INDEX
C
         IDPA = - INDXPH(ID)
      IF ( IDPA .LT. IPLOW .OR. IDPA .GT. IPHIGH ) GO TO 100
C     ... we do not need this distribution this time
C
         ICHO  = INDXPH(IC)
         ICHDP = ( IDPA - IPLOW ) * NH + ICHO
C
         ISYMC = ISMO(IC)
         ISYMD = ISMO(ID)
         ISYMCD= MULD2H(ISYMC,ISYMD)
C
C        MULLIKEN DISTRIBUTION (**,CD) IN H2M(*,*)
C        DIRAC PARTICLE-HOLE DISTRIBUTION IN H2PHD(*,*,ICHPH) C<D
C
C        ADD SECOND ORDER CONTRIBUTIONS FROM INTEGRAL DISTRIBUTIONS CD
C        TO ONE ELECTRON DENSITY MATRIX; copy H2M to WRK(KH2M1) because
C        H2M is modified in MP2DEN and it will be needed below if
C        MP2FRO true.
C
         IF (MP2FRO) CALL DCOPY(N2ORBX,H2M,1,WRK(KH2M1),1)
         IF ( IFRMP2(IC).EQ.0 ) THEN
C        If ( IC not frozen   ) then
            CALL MP2DEN(NFRMP2,IFRMP2,ICHO,IDPA,ICHDP,ISYMCD,H2PHD,DONE,
     *                  COEMP2,IADR1,ORBEN,H2M,H2D,EMP2,0)
C           CALL MP2DEN(NFROZ,IFROZ,ICHO,IDPA,ICHDP,ISYMCD,H2PHD,DONE,
C     *                 COEMP2,IADR1,ORBEN,H2M),H2D,EMP2,KFRO)
         END IF
C
C        IF MP2FRO, that is, if some orbitals are frozen above,
C        then also collect full MP2 density matrix and energy.
C        IFRRHF is defined to correspond to only RHF frozen orbitals.
C
Cmjp  If the last argument of MP2DEN < 0, then we will form the energy
Cmjp  in MP2DEN assuming no frozen orbitals. If the argument is
Cmjp  >= 0 then we take account of any frozen orbitals.
C
         IF ( MP2FRO ) THEN
         IF ( IFRRHF(IC) .EQ. 0 ) THEN
            CALL MP2DEN(NFRRHF,IFRRHF,ICHO,IDPA,ICHDP,ISYMCD,H2PHD,
     *                  DONEFL,COEMP2,IADR1,ORBEN,WRK(KH2M1),H2D,
     *                  EMP2FL,-1)
         END IF
         END IF
C
C        Go get next Mulliken distribution
C
      GO TO 100
C
C     arrive at 800 when finished with all Mulliken distributions
C
  800 CONTINUE
C
      IF (L_MP2_TDA) THEN
!        write(0,*) 'hjaaj MP2MUL: writing RJ and RK to LUJK=99'
         LUJK = 99
         REWIND (LUJK)
         CALL WRITT (LUJK,N2ORBX,WRK(KJ))
         CALL WRITT (LUJK,N2ORBX,WRK(KK))
         REWIND (LUJK)
      END IF
C
      CALL MEMREL('MP2MUL',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('MP2MUL')
      RETURN
      END
C  /* Deck fckifr */
      SUBROUTINE FCKIFR(IFRRHF)
C
C  941005-hjaaj
C  Set IFRRHF(NORBT) to 1 if orbital is frozen in RHF
C  otherwise to 0
C
#include <implicit.h>
      DIMENSION IFRRHF(NORBT)
C
C Used from common blocks:
C   INFORB : NSYM, NORBT, IORB(8)
C   SCBRHF : NFRRHF(8)
C
#include <inforb.h>
#include <scbrhf.h>
C
      CALL IZERO(IFRRHF,NORBT)
      DO 100 ISYM = 1,NSYM
         IF (NFRRHF(ISYM) .GT. 0) THEN
            IORBI = IORB(ISYM)
            DO 200 I = IORBI+1,IORBI+NFRRHF(ISYM)
               IFRRHF(I) = 1
  200       CONTINUE
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck mp2phd */
      SUBROUTINE MP2PHD(H2PHD,H2CD,IPLOW,IPHIGH,WRK,KFRSAV,LFRSAV)
C
C (retrieve Particle-Hole integrals in Dirac format)
C Written  18. May 1987 by Hans Jorgen Aa. Jensen.
C Revisions:
C   891017-hjaaj: use NXTH2M
C
C Purpose:
C   For MP-2 density matrix we need particle-hole integrals in Dirac
C   format.
C   Sort Mulliken integrals (i a / j b) = ( A B / C D )
C   to   Dirac    format    <i b / a j> = < A D / B C >
C   (in H2PHD(i,b;a,j) ; PH for particle-hole, D for Dirac)
C
C Scratch:
C   WRK(kfrsav:kfrsav-1+lfrsav)
C   H2CD(NORBT,NORBT)
C
#include <implicit.h>
      DIMENSION H2PHD(NP,NH,*), H2CD(NORBT,NORBT)
      DIMENSION WRK(*)
C
C Used from common blocks:
C   INFORB : NSYM,NRHF(*),NORB(*),IORB(*),?
C   INFIND : ISMO(*),?
C   INFMP2 : NP,NH,?
C   INFPRI : ?
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include <infmp2.h>
#include <infpri.h>
C
      LOGICAL AHOLE, BHOLE, CHOLE, DHOLE

      DIMENSION NEEDTP(6)
      DATA NEEDTP /0,1,1,1,1,0/
C     ... only particle-hole distributions needed
C         assumed: NRHF(isym) .ge. NISH(isym)
C              and NVIR(isym) .le. NSSH(isym)
C
      CALL QENTER('MP2PHD')
C
C     >>> zero H2PHD
C
      LENGTH = NP * NH * NH * (IPHIGH - IPLOW + 1)
      CALL DZERO(H2PHD,LENGTH)
C
C
C ****************************************************************
C     Loop over Mulliken distributions allowed in NEEDTP(6)
C
      KFREE = KFRSAV
      LFREE = LFRSAV
      IDIST = 0
  100 CALL NXTH2M(IC,ID,H2CD,NEEDTP,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 800
C     ... if idist .lt. 0 then no more distributions
C
         ICSYM = ISMO(IC)
         IDSYM = ISMO(ID)
         CHOLE = ( INDXPH(IC) .GT. 0 )
         DHOLE = ( INDXPH(ID) .GT. 0 )
C
C        See if we need this distribution
C        (We only need particle-hole distributions)
C
      IF (CHOLE .EQV. DHOLE) GO TO 100
C     ... we do not need this distribution
         IF (DHOLE) THEN
            ISWAP = IC
            IC    = ID
            ID    = ISWAP
         END IF
C        ... now c is hole, d is particle
C            Do we want this distribution this time ?
         IDPA   = - INDXPH(ID)
      IF ( IDPA .LT. IPLOW .OR. IDPA .GT. IPHIGH ) GO TO 100
C     ... No, we don't
         ICHO   = INDXPH(IC)
         IDPOFF = (IDPA - IPLOW) * NH
         ICDSYM = MULD2H(ICSYM,IDSYM)
C
C        loop over ( A B / C D ) of ( p h / h p ) type
C          sort as < A C / B D > in / B D > distribution.
C                ( < A C / B D > is thus < p h / h p > type )
C
         DO 300 IA = 1,NORBT
            AHOLE = ( INDXPH(IA) .GT. 0 )
            IF (AHOLE) THEN
               IAH   =   INDXPH(IA)
               IASYM = ISMO(IA)
               IBSYM = MULD2H(IASYM, ICDSYM)
               IBST  = IORB(IBSYM) + NRHF(IBSYM) + 1
               IBEND = IORB(IBSYM) + NORB(IBSYM)
               DO 200 IB = IBST,IBEND
                  IBP = - INDXPH(IB)
                  H2PHD(IBP,ICHO,IDPOFF+IAH) = H2CD(IB,IA)
  200          CONTINUE
            END IF
  300    CONTINUE
C
C        Go get next Mulliken distribution
C
      GO TO 100
C
C     arrive at 800 when finished with all Mulliken distributions
C
  800 CONTINUE
C
C *** end of subroutine MP2PHD
C
      CALL MEMREL('MP2PHD',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('MP2PHD')
      RETURN
      END
C  /* Deck mp2set */
      SUBROUTINE MP2SET(IADR1)
C
C 19-May-1987 hjaaj
C 20-MAY-1987 PJ     INSERTED IPHORD
Cold  SUBROUTINE MP2SET(MP2FRO,INDXPH,IPHORD,IFRMP2,NP,NH)
C 011293 mjp MP2FRO,IPHORD,INDXPH,IFRMP2,NP,NH now put into /INFMP2/
C
C Calculate NP, number of particles,
C           NH, number of holes, and
C           INDXPH(NORBT), particle or hole index for each orbital.
C           IPHORD(NORBT), GIVE SYMMETRY ORDERED INDEX FOR PARTICLE OR
C                          HOLE INDEX
C           IFRMP2(NORBT), .EQ.1 FOR FROZEN MP2 ORBITAL ELSE .EQ.0
C
#include <implicit.h>
C
C Used from common blocks:
C  INFORB : NSYM, NRHF(NSYM), NVIR(NSYM), NFRMP2(NSYM), ...
C  INFMP2 : IPHORD, MP2FRO, INDXPH, IFRMP2,NP,NH
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infmp2.h>
      DIMENSION IADR1(MAXVIR,MAXOCC)
#include <infind.h>
C
      CALL QENTER('MP2SET')
C
      NH = NRHFT
      NP = NORBT - NRHFT
C
      IPH = 0
      IH  = 0
      IP  = 0
      IORBT = 0
      CALL IZERO(IFRMP2,NORBT)
      IMP2 = 0
      DO 300 ISYM = 1,NSYM
         NRHFI = NRHF(ISYM)
         NVIRI = NVIR(ISYM)
         NFRMPI= NFRMP2(ISYM)
         DO 100 I = 1,NRHFI
            IPH = IPH + 1
            IH  = IH  + 1
            IPHORD(IH)  = IPH
            INDXPH(IPH) = IH
            IF (I .LE. NFRMPI) THEN
               IMP2        = IMP2 + 1
               IFRMP2(IPH) = 1
            END IF
  100    CONTINUE
         DO 200 I = 1,NVIRI
            IPH = IPH + 1
            IP  = IP  + 1
            IPHORD(NH+IP) = IPH
            INDXPH(IPH)   = -IP
            IF (NRHFI+I .LE. NFRMPI) THEN
C           ... code for freezing low-lying virtuals /hjaaj Oct 09
               IMP2        = IMP2 + 1
               IFRMP2(IPH) = 1
            END IF
  200    CONTINUE
  300 CONTINUE
C
C     Any frozen orbital in MP2?
C
      MP2FRO = (IMP2 .GT. 0)
C
C Calculate the number of ph excitations of each symmetry, ignoring
C orbitals frozen in MP2.
C
         CALL IZERO(NPHSYM,NSYM)
         DO I = 1,NH
            II = IPHORD(I)
            IF (IFRMP2(II).EQ.1) GOTO 500
            DO J = 1,NP
               JA = IPHORD(NH+J)
               KSYM = MULD2H(ISMO(II),ISMO(JA))
               NPHSYM(KSYM) = NPHSYM(KSYM) + 1
            ENDDO
 500     CONTINUE
         ENDDO
         NPHMAX = NPHSYM(1)
         DO I = 2,NSYM
            NPHMAX = MAX(NPHMAX,NPHSYM(I))
         ENDDO
C
C Determine ph offset array IADR1(VIRT,OCC) for symmetry 1
C
         CALL IZERO(NPHTOT,NSYM)
         CALL PHADR1(IADR1,1)
C
      CALL QEXIT('MP2SET')
      RETURN
      END
C  /* Deck phpack */
      SUBROUTINE PHPACK (XK, XKP, JPHSYM, IWAY, IPHLAB)
C
C Pack or unpack ph-array for MP2/SOPPA.
C
C If IWAY.LT.0 then unpack the linear array XKP into XK. Else pack the
C   the array XK into XKP. 
C If IPHLAB>0, pack XK(OCC,VIRT); IPHLAB<0 C pack XK(VIRT,OCC),
C   and IPHLAB=0, pack both (XK(OCC,VIRT) and XK(VIRT,OCC))
C
#include <implicit.h>
      DIMENSION XK(NORBT,NORBT), XKP(*)
C
C Used from common blocks:
C   INFORB: NSYM,NORBT,MULD2H, N2ORBX , NOCC(), ...
C   INFMP2: NFRMP2(*)
C
#include <maxorb.h>
#include <inforb.h>
#include <infmp2.h>
C
      IF (IWAY .LT. 0) CALL DZERO(XK,N2ORBX)
C
      IA = 0
      DO ISYM = 1,NSYM
         JSYM = MULD2H(ISYM,JPHSYM)
         IF (NOCC(ISYM) .EQ. 0) GO TO 600
         IIST = IORB(ISYM) + NFRMP2(ISYM) + 1
         IIEND= IORB(ISYM) + NOCC(ISYM)
         JAST = IORB(JSYM) + NOCC(JSYM) + 1
         JAEND= IORB(JSYM) + NOCC(JSYM) + NSSH(JSYM)
         IF (IWAY .GE. 0) THEN
            IF (IPHLAB .GE. 0) THEN
               DO II = IIST,IIEND
                  DO JA = JAST, JAEND
                     IA = IA + 1
                     XKP(IA) = XK(II,JA)
                  END DO
               END DO
            ELSE
               DO II = IIST,IIEND
                  DO JA = JAST, JAEND
                     IA = IA + 1
                     XKP(IA) = XK(JA,II)
                  END DO
               END DO
            END IF
         ELSE
            IF (IPHLAB .GE. 0) THEN
               DO II = IIST,IIEND
                  DO JA = JAST,JAEND
                     IA = IA + 1
                     XK(II,JA) = XKP(IA)
                  END DO
               END DO
               IF (IPHLAB .EQ. 0) THEN
                  DO II=IIST,IIEND
                     DO JA=JAST,JAEND
                        XK(JA,II) = XK(II,JA)
                     ENDDO
                  ENDDO
               ENDIF
            ELSE
               DO II = IIST,IIEND
                  DO JA = JAST,JAEND
                     IA = IA + 1
                     XK(JA,II) = XKP(IA)
                  END DO
               END DO
            END IF
         END IF
 600  CONTINUE
      END DO
      RETURN
      END
C  /* Deck phadr1 */
      SUBROUTINE PHADR1(IADR1,KSYMOP)
C
C Written by Martin Packer 090394.
C
C IADR1(A,I) is an address C array giving the offsets to
C   the set of ph excitations associated with A,I.
C
#include <implicit.h>
C
C
#include <maxorb.h>
#include <maxash.h>
#include <infind.h>
#include <inforb.h>
#include <infmp2.h>
C
      INTEGER A, AA
      DIMENSION IADR1(MAXVIR,MAXOCC)
C
C
      CALL ISETVC(IADR1,-100 000 000,MXCORB*MAXOCC)
C     call isetvc(ivector,ivalue,ndim)
      IAB = 0
      NPHTOT(KSYMOP) = 0
      DO II = 1,NH
        I = IPHORD(II)
        IF (IFRMP2(I) .EQ. 0) THEN
           DO AA = 1,NP
              A = IPHORD(NH+AA)
              IASYM = MULD2H(ISMO(I),ISMO(A))
              JBSYM = MULD2H(IASYM,KSYMOP)
              NPHTOT(KSYMOP) = NPHTOT(KSYMOP) + NPHSYM(JBSYM)
              IADR1(AA,II) = IAB
              IAB = IAB + NPHSYM(JBSYM)
           ENDDO
         END IF
      ENDDO
      RETURN
      END
C  /* Deck phadr2 */
      SUBROUTINE PHADR2(IADR2,NCOL,NROW)
C
C 30-Nov-1995 Hans Joergen Aa. Jensen.
C
C The address of the MP2 first order 2p2h coefficients
C   kappa(ia,jb) er IADR1(B,J) + IADR2(A,I)
C
C Based on code written by Erik Kragh Dalskov fall 1994.
C
#include <implicit.h>
      DIMENSION IADR2(NCOL,NROW)
C
#include <maxorb.h>
#include <maxash.h>
#include <infind.h>
#include <inforb.h>
#include <infmp2.h>
C
      INTEGER A,AA,ASYM,A1OFF,AST
C
      IF (NOCCT .GT. NROW .OR. NSSHT .GT. NCOL) THEN
         CALL QUIT('prog.error: dimension of IADR2 in PHADR2 too small')
      END IF
C
      CALL ISETVC(IADR2,-100 000 000,NCOL*NROW)
C     call isetvc(ivector,ivalue,ndim)
C
      DO IASYM = 1,NSYM
         IA = 0
         DO I = 1,NH
            II = IPHORD(I)
         IF (IFRMP2(II).EQ.1) GOTO 500
            ISYM = ISMO(II)
            ASYM = MULD2H(ISYM,IASYM)
            AST  = IORB(ASYM) + NOCC(ASYM) + 1
            A1OFF= -INDXPH(AST) - 1
            DO A = 1,NSSH(ASYM)
               IADR2(A1OFF+A,I) = IA+A
            ENDDO
            IA = IA + NSSH(ASYM)
 500     CONTINUE
         ENDDO
      ENDDO
      RETURN
      END
C --- end of sirius/sirmp2.F ---
