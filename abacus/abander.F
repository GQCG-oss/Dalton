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
C  FILE: abacus/abander.F
C
C  /* Deck nmdinp */
      SUBROUTINE NMDINP(WORD,IDRPRI)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NDIR=3,NTABLE = 21)
      CHARACTER PROMPT*1, WORD*7, GRPTMP*15, TABLE(NTABLE)*7,
     &          TABDIR(NDIR)*7, WORD1*7
#include <numder.h>
#include <fcsym.h>
#include <cbinum.h>
#include <abainf.h>
#include <cbiwlk.h>
#include <cbivib.h>
      LOGICAL NEWDEF
C
      DATA TABDIR/'*PROPAV','*XXXXXX','*VIBANA'/
C
      DATA TABLE /'.DORDR ', '.SYMMET', '.SDRTST', '.RESTRT', '.DRYRUN',
     *            '.XXXXXX', '.NORMAL', '.PRECAL', '.REUSE ', '.XXXXXX',
     *            '.VIBANA', '.TEST N', '.DISPLA', '.PROPER', '.PRINT ',
     *            '.MANUAL', '.HARMON', '.SPECTR', '.MIDAS ', '.THRMID',
     *            '.MINOUT'/
C
C     *** Initializing variables for *VIBANA and *HARMON. ***
      CALL NVBINI
      CALL VIBINI
C
C
C     ************************************************
C     **** Finding the analytical differentiation ****
C     **** order of the  energy for the           ****
C     **** wavefunction used.                     ****
C     ************************************************
C
      CALL FNDANA(NAORDR)
C
      ICHANG = 0
      WORD1 = WORD
 100  CONTINUE
      READ (LUCMD, '(A7)') WORD
      CALL UPCASE(WORD)
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 100
      ELSE IF (PROMPT .EQ. '.') THEN
         ICHANG = ICHANG + 1
         DO 200 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
     *                20,21),I
            END IF
 200     CONTINUE
         IF (WORD .EQ. '.OPTION') THEN
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            GO TO 100
         END IF
         WRITE (LUPRI,'(/3A/)') ' Keyword "',WORD,
     *        '" not recognized for '//WORD1
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword for '//WORD1)
 1       CONTINUE
            READ (LUCMD,'(2I4)') NMORDR, NAORDR
            GOTO 100
 2       CONTINUE
            READ (LUCMD,'(A)') GRPTMP
            IJ = 0
            DO II = 1, 15
               IF (GRPTMP(II:II).NE.' ') THEN
                  IJ = IJ + 1
                  FCLASS(IJ:IJ) = GRPTMP(II:II)
               END IF
            END DO
            GOTO 100
 3       CONTINUE
            SDRTST = .TRUE.
            GOTO 100
 4       CONTINUE
            RESTRT = .TRUE.
            GOTO 100
 5       CONTINUE
            DRYRUN = .TRUE.
            READ (LUCMD,*)  NMREDU
            READ (LUCMD,*) (KDRYRN(II),II=1,NMREDU)
            GOTO 100
 6       CONTINUE
            GOTO 100
 7       CONTINUE
            NRMCRD = .TRUE.
            GOTO 100
 8       CONTINUE
            PREHES = .TRUE.
            GOTO 100
 9       CONTINUE
            REUHES = .TRUE.
            GOTO 100
 10      CONTINUE
            GOTO 100
 11      CONTINUE
            NUMVIB = .TRUE.
            GOTO 100
 12      CONTINUE
            NRMCRD = .TRUE.
            HTEST  = .TRUE.
            GOTO 100
 13      CONTINUE
            READ (LUCMD, *) DISPLC
            GOTO 100
 14      CONTINUE
            NPRPDR = .TRUE.
            READ (LUCMD, *) NMRDRP, NARDRP
            IF (NMORDR.EQ.0) NAORDR = 0
            NMORDR = MAX(NMRDRP,NMORDR)
            GOTO 100
 15      CONTINUE
            READ (LUCMD, *) IDRPRI
            GOTO 100
 16      CONTINUE
            MANUAL = .TRUE.
            GOTO 100
 17      CONTINUE
            HARMON = .TRUE.
            VIB    = .TRUE.
            MAXDIF = 2
            GOTO 100
 18      CONTINUE
            SPECTR = .TRUE.
            GOTO 100
 19      CONTINUE
            MIDAS  = .TRUE.
            GOTO 100
 20      CONTINUE
            READ (LUCMD, *) XTHR
            THRMID = ABS(XTHR)
 21      CONTINUE
            MINOUT = .TRUE.
            GOTO 100
      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 300
      ELSE
         WRITE (LUPRI,'(/4A/)') ' Prompt "',WORD,
     *        '" not recognized for ',WORD1
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt for '//WORD1)
      END IF
C
C     *** Print section. ***
C
 300  CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for '//WORD1//':',0)

         IF (NUMVIB) THEN
            WRITE (LUPRI,'(/5X,A)') 'A vibrational analysis is done.'
            WRITE (LUPRI,'(5X,A)')  'Which/how is specified in ' //
     &           '*VIBANA (and **PROPERTIES).'
         ELSE
            WRITE (LUPRI,'(/5X,A,I4/5X,A,I4,A)')
     &         'Numerical derivatives calculated to order', NMORDR,
     &         'using analytical', NAORDR, '. derivatives'
         END IF
C
         WRITE (LUPRI,'(5x,A)') 'Group used for force constants ' //
     &                           GRPTMP
         WRITE (LUPRI,'(5x,A,F6.4)') 'Step size used: ', DISPLC
         IF (SDRTST) THEN
            WRITE (LUPRI,'(/5X,A)') 'Comparison of numerical Hessian '
     &          // 'with analytical Hessian is performed'
         END IF
         IF (DRYRUN) THEN
            WRITE (LUPRI,'(/5X,A)') 'Numerical derivatives will be' //
     &           'conducted as a dry run.'
            WRITE (LUPRI,'(5X,A)') 'No actual derivatives will be' //
     &           'calculated.'
            WRITE (LUPRI,'(A,I5)') 'Number of redundant coordinates:',
     &           NMREDU
         END IF
C
         IF (RESTRT) THEN
            WRITE (LUPRI,'(/5X,A)') 'This is a restart of an old run.'
         END IF
C
         IF (NRMCRD) THEN
            WRITE (LUPRI,'(/5X,A)') 'Normal coordinates will be found.'
            WRITE (LUPRI,'(5X,A)') 'Energy and property derivatives' //
     &           'will be with respect to these coordinates.'
         END IF
C
         IF (PREHES) THEN
            WRITE (LUPRI,'(/5X,A)') 'A precalculated hessian will be' //
     &           ' used to find normal coordinates.'
         END IF
C
         IF (REUHES) THEN
            WRITE (LUPRI,'(/5X,A)') 'Hessian (if specified elsewhere)'
     &           // ' will be saved for future work.',
     &      'Hessian will be saved on the file "DALTON.HES"'
         END IF
C
         IF (NRMCRD.AND.HTEST) THEN
            WRITE (LUPRI,'(/5X,A)') 'A test of the normal coordinates '
     &           // 'will be done.'
         END IF
C
         IF (MANUAL) THEN
            WRITE (LUPRI,'(/A)')
     &         'The mol file will be printed for geometries.'
         END IF
C
         IF (NPRPDR) THEN
            WRITE (LUPRI,'(/5X,A)')
     &         'Property derivatives will be calculated.'
            IF (NMORDR.GT.0) THEN
               WRITE (LUPRI,'(5X,A,I4)') 'Order of the' //
     &              'differentiation is equal to: ', NMORDR
            END IF
         END IF
         IF (MIDAS) THEN 
            WRITE (LUPRI,'(/5X,A)') 
     &                  'Operator file for MidasCpp interface '//
     &                  'will be generated.'
            WRITE (LUPRI,'(5x,A,A,E24.10)')
     &             'Threshold for term coefficient relative to largest',
     &             ' harmonic term: ',THRMID
         ENDIF
      END IF
C
C
C     *** Different * sections. ***
C
 400  CONTINUE
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 400
      ELSE IF (PROMPT .EQ. '*') THEN
         DO 500 I = 1, NDIR
            IF (WORD .EQ. TABDIR(I)) THEN
               GO TO
     *         (101,102,103), I
            END IF
 500     CONTINUE
         IF (WORD(1:2) .EQ. '**') GO TO 600
         WRITE (LUPRI,'(/,3A,/)') ' Directory ',WORD,' nonexistent.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal directory in ABAINP.')
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal or',
     *                        ' out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in ABAINP, error in prompt.')
      END IF
 101  CONTINUE
        CALL NVBINP(WORD)
        GOTO 400
 102  CONTINUE
        GOTO 400
 103  CONTINUE
        CALL VIBINP(WORD)
        GOTO 400
C
 600  CONTINUE
C
      RETURN
      END
C
C
C     /* Deck fndana */
      SUBROUTINE FNDANA(NAORDR)
C     *******************************************************
C     **** Subroutine that keeps track of the analytical ****
C     **** for a given wave-function. NAORDR gives the   ****
C     **** order of the analytical derivative.           ****
C     *******************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
C
      PARAMETER (D0 = 0.0D0)
#include <gnrinf.h>
#include <inforb.h>
#include <dftcom.h>
#include <ecpinf.h>
#include <ccsdinp.h>
      LOGICAL MCSCF, DMP2, DOHFSRDFT, DOROHF, DOCISRDFT, DOMCSRDFT, 
     &        WAVTP, DPCM
C
C     *** Workaround to many common variables with the same name ***
C     *** MCSCF is set to DOMC (in infinp.h), and DMP2 is set    ***
C     *** to DOMP2 (in infinp.h through the logical function     ***
C     *** WAVTP.                                                 ***
C
      MCSCF      = WAVTP('MCSCF')
      DMP2       = WAVTP('MP2')
      DPCM       = WAVTP('PCM')
      DOROHF     = WAVTP('ROHF')
#ifdef MOD_SRDFT
      DOHFSRDFT  = WAVTP('HFSRDFT')
      DOCISRDFT  = WAVTP('CISRDFT')
      DOMCSRDFT  = WAVTP('MCSRDFT')
#else
      DOHFSRDFT  = .FALSE.
      DOCISRDFT  = .FALSE.
      DOMCSRDFT  = .FALSE.
#endif
C
      IF (CCSDT .OR. CIS   .OR. CC1A .OR. CC1B .OR. MCC2 .OR. CCP2 .OR.
     &    CCD   .OR. CC3   .OR. CCP3 .OR. CCRT .OR. CCR3 .OR.
     &    CCR1A .OR. CCR1B .OR. CCT  .OR.
     &    DOHFSRDFT .OR. DOCISRDFT .OR. DOMCSRDFT .OR.
     &    (DMP2.AND..NOT.MCSCF) .OR. ECP .OR. DKTRAN) THEN
         NAORDR = 0
      ELSE IF (CCSD .OR. CCS .OR. CC2 .OR. MP2 .OR.
     &         CCPT .OR. (DFTRUN .AND. (NASHT.GT.0)) .OR.
     &        (DFTRUN .AND. (NSYM.GT.1.OR.HFXMU.NE.D0)) .OR.
     &        DPCM .OR. DOROHF) THEN
         NAORDR = 1
      ELSE
         NAORDR = 2
      END IF
C
      WRITE (LUPRI,'(/5X,A,I4)') 'Order of analytical ' //
     &          'energy-derivatives available:', NAORDR
      WRITE (LUPRI,'(5X,A,I4/)') 'This will be the default.'
C
      RETURN
      END
C
C     /* Deck wavtp */
      LOGICAL FUNCTION WAVTP(STRING)
C     **************************************************************
C     *** Workaround to many common variables with the same name ***
C     *** MCSCF is set to DOMC (in infinp.h) through the logical ***
C     *** function WAVTP.                                        ***
C     **************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include "pcmlog.h"
#include <infinp.h>
      CHARACTER*(*) STRING
C
      IF (STRING.EQ.'MCSCF') THEN
         WAVTP = DOMC
      ELSE IF (STRING.EQ.'MP2') THEN
         WAVTP = DOMP2
      ELSE IF (STRING.EQ.'ROHF') THEN
         WAVTP = HSROHF
      ELSE IF (STRING .EQ. 'PCM') THEN
         WAVTP = PCM
#ifdef MOD_SRDFT
      ELSE IF (STRING.EQ.'HFSRDFT') THEN
         WAVTP = DOHFSRDFT
      ELSE IF (STRING.EQ.'CISRDFT') THEN
         WAVTP = DOCISRDFT
      ELSE IF (STRING.EQ.'MCSRDFT') THEN
         WAVTP = DOMCSRDFT
#endif
      ELSE
         WRITE (LUPRI,'(/2A)')
     &   'Undefined string in WAVTP :', STRING
         CALL QUIT('Wrong string in WAVTP')
      END IF
!     write(lupri,'(2A,T30,A,L10)')
!    &   'WAVE FUNCTION TYPE ', STRING,' : ', WAVTP
C
      RETURN
      END
C
C
C     /* Deck nmdini */
      SUBROUTINE NMDINI(IPRINT)
C
C     Initialize /NUMDER/, /FCSYM/ and some /ABAINF/
C
#include <implicit.h>
#include <mxcent.h>
#include <numder.h>
#include <fcsym.h>
#include <cbinum.h>
#include <abainf.h>
#include <cbiwlk.h>
C
C     Print variable.
      IPRINT = 0
C
C     /CBINUM/
      NRMCRD = .FALSE.
      PGMTST = .FALSE.
      HTEST  = .FALSE.
      PREHES = .FALSE.
      REUHES = .FALSE.
      ANALZ1 = .FALSE.
      NUMVIB = .FALSE.
      NPRPDR = .FALSE.
      HARMON = .FALSE.
      SPECTR = .FALSE.
      MIDAS  = .FALSE.
      MINOUT = .FALSE.
      THRMID = 1.0D-15
C
C     /NUMDER/
      NMORDR = 0
      NAORDR = 0
      NMDPRP = 0
      NMRDRP = 0
      NARDRP = 0
      NMPINI = 0
      NWPROP = .FALSE.
      FSTPRP = .FALSE.
      NOMOVE = .FALSE.
      NUMELC = .FALSE.
      CMPARE = .FALSE.
      SDRTST = .FALSE.
      DRYRUN = .FALSE.
      FRSTNM = .FALSE.
      PRPVIB = .FALSE.
      MANUAL = .FALSE.
C
C     /FCSYM/
      FCLASS = 'C1             '
      MROTAX = .FALSE.
      VPLANE = .FALSE.
      HPLANE = .FALSE.
      ROTAX2 = .FALSE.
      DPLANE = .FALSE.
      ICNTR  = .FALSE.
      ROTARE = .FALSE.
      SEPDEG = .FALSE.
C
C     /ABAINF/
      VIB    = .FALSE.
C
C     /CBIWLK/
      DISPLC = 1.0D-2
C
      RETURN
      END
C
C
C     /* Deck numdrv */
      SUBROUTINE NUMDRV(WORK,LWORK,IPRINT,WRKDLM)
C
C     Driver routine for numerical differentiation
C

#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <infpar.h>
C
c#if defined (VAR_MPI)
c      INCLUDE 'mpif.h'
c      LOGICAL FINISH
c#endif
#include <cbirea.h>
#include <cbiwlk.h>
#include <cbinum.h>
#include <trkoor.h>
#include <nuclei.h>
#include <symmet.h>
#include <numder.h>
#include <molinp.h>
#include <fcsym.h>
#include <gnrinf.h>
#include <huckel.h>
      DIMENSION WORK(LWORK)
      LOGICAL MOLECU
      CHARACTER WORD*7
C
      CALL QENTER('NUMDRV')
C
C     feb 11 - hjaaj
C     cut down on hermit and abacus output during numerical
C     differentiation
C
      IPRUSR_orig = IPRUSR
      IPREAD_orig = IPREAD
      IF (USRIPR) THEN
C        if user has asked for higher print level, no change
         IPRUSR_reduced = IPRUSR
         IPREAD_reduced = IPREAD
      ELSE
         IPRUSR_reduced = -2
         IPREAD_reduced = -2
      END IF
      IPRUSR = IPRUSR_reduced
      IPREAD = IPREAD_reduced
C
C 100  CONTINUE
C
      MOLECU = .TRUE.
      NCOOR  = 3*NUCDEP
C
      NDERIV = 0
      NDIME  = 1
      NINNER = 1
      IF (NAORDR.GE.1) NINNER = NCOOR*NINNER
      IF (NAORDR.GE.2) NINNER =(NCOOR+1)*NINNER/2
C
      MAXADR = 1
      DO 200 J = NMORDR, 1, -2
      DO 200 I = 1, J
         MAXADR = MAXADR + 1
 200  CONTINUE
C
      DO 300 IORDR = 3, NMORDR+NAORDR
         IKDRV = 1
         DO 400 IIORDR =1, IORDR
            IKDRV = IKDRV*(NCOOR+IIORDR-1)/IIORDR
 400     CONTINUE
         NDERIV = NDERIV + IKDRV
 300  CONTINUE
      if (nderiv .le. 0) then
         write(lupri,*) 'hjaaj: nderiv .le. 0',nderiv
         write(lupri,*) 'hjaaj: nmordr,naordr,ncoor',nmordr,naordr,ncoor
      end if
C
      DO 500 IORDR = 1, NMORDR
         IKDIME = 1
         DO 600 IIORDR = 1, IORDR
            IKDIME = IKDIME*(NCOOR+1-IIORDR)/IIORDR
 600     CONTINUE
         NDIME = NDIME + IKDIME*2**IORDR
 500  CONTINUE
      IF (NMORDR .GE. 3) THEN
         NDIME = NDIME + 2*NCOOR
      END IF
      IF (NMORDR .GE. 4) THEN
         NDIME = NDIME + 4*NCOOR*(NCOOR-1)
      END IF
      IF (NMORDR .GE. 5) THEN
         NDIME = NDIME + 2*NCOOR
      END IF
C
      KTEST  = 2
Chjaaj-Oct07: KTEST a simple device to check if arrays which
C     supposedly not are used, are used anyway.
      KDERIV = KTEST  + 1
      KFUVAL = KDERIV + NDERIV
      KCOOR  = KFUVAL + NDIME*NINNER
      KCSTAR = KCOOR  + 3*NCOOR
      KSYMCO = KCSTAR + 3*NCOOR
      KTRNRC = KSYMCO +   NCOOR**2
      KTRMSS = KTRNRC +   NCOOR**2
      KDKIN  = KTRMSS +   NCOOR
      KFREQ  = KDKIN  +   NCOOR
      KRNNRM = KFREQ  +   NCOOR
      KLAST  = KRNNRM +   NCOOR
C
C     *** Memory needed for test on Hessian. ***
C
      IF (SDRTST) THEN
         KTSTGD = KLAST
         KTSTSD = KTSTGD + NCOOR
         KLAST  = KTSTSD + NCOOR**2
      ELSE
         KTSTGD = KTEST
         KTSTSD = KTEST
      END IF
      IF (PGMTST) THEN
         LTHTST = 2**NMORDR
         KENTST = KLAST
         KLAST  = KENTST + LTHTST
      ELSE
         KENTST = KTEST
      END IF
C
C     *** Memory for screening of rendundant force constants. ***
C
      LDPMTX = 0
      IF (NMORDR .GE. 4) THEN
         LDPMTX = LDPMTX + (NCOOR*(NCOOR+1)*(NCOOR+2)*(NCOOR+3))/24
      END IF
      IF (NMORDR .GE. 3) THEN
         LDPMTX = LDPMTX + (NCOOR*(NCOOR+1)*(NCOOR+2))/6
      END IF
      NSTRDR = NMORDR + 1
C     ... hjaaj Dec 07: used for allocation, and NMORDR+1 is sometimes referenced
      IF (NMORDR .GE. 2) THEN
         IFRSTD = 2**NMORDR
         LDPMTX = LDPMTX + (NCOOR*(NCOOR+1))/2
         KDPMTX = KLAST
chj      KDCOEF = KDPMTX +   IFRSTD*NMORDR*LDPMTX
         KDCOEF = KDPMTX +   IFRSTD*NSTRDR*LDPMTX
         KNIDPC = KDCOEF +   IFRSTD       *LDPMTX
         KLAST  = KNIDPC +                 LDPMTX
      ELSE
         KDPMTX = KTEST
         KDCOEF = KTEST
         KNIDPC = KTEST
      END IF
C
C     *** Space for backup of MOLECULE.INP ***
C
      KMBKLN = KLAST
      KLAST  = KMBKLN + 80*NMLINE
C
C     *** Space for backup of isotopes if abacus is run. ***
C
      IF (NAORDR.GT.0) THEN
         KISOTP = KLAST
         KLAST  = KISOTP + NUCDEP
      ELSE
         KISOTP = KTEST
      END IF
C
C     *** Symmetry initialization and symmetry ***
C     ***       related memory allocation      ***
C
      CALL FCSINI
      KGRIRP = KLAST
      KCHRCT = KGRIRP + NGORDR*NGVERT
      KICRIR = KCHRCT + NGORDR*NCVERT
      KLAST  = KICRIR + 2*NCOOR
C
      LWRK1 = LWORK - KLAST + 1
C
Chjaaj-Oct07: KTEST is a simple device to check if arrays which
C     supposedly not are used, are used anyway.
      WORK(KTEST) = -999.9D0
      CALL NUMDR1(WORK(KDERIV),WORK(KFUVAL),WORK(KCOOR),WORK(KCSTAR),
     &            WORK(KSYMCO),WORK(KDCOEF),WORK(KTSTGD),WORK(KTSTSD),
     &            WORK(KENTST),WORK(KGRIRP),WORK(KCHRCT),WORK(KTRNRC),
     &            WORK(KTRMSS),WORK(KDKIN) ,WORK(KFREQ), WORK(KRNNRM),
     &            WORK(KLAST) ,WORK(KDPMTX),WORK(KNIDPC),WORK(KICRIR),
     &            WORK(KISOTP),WORK(KMBKLN),LWRK1,NDERIV,NDIME,
     &            NINNER,MAXADR,LTHTST,LDPMTX,IFRSTD,IPRINT,WRKDLM)
      IF (WORK(KTEST) .NE. -999.9D0) THEN
         CALL QUIT('WORK(KTEST) has been modified!')
      END IF
C
C
c#if defined (VAR_MPI)
cC
cC     We let the slaves wait for the Master to tell them whether to pick up
cC     a new geometry or to end this calculation
cC
c      IF (MYNUM .GT. 0) THEN
c         CALL MPI_BCAST(NTASK,1,my_MPI_INTEGER,NODEID(MASTER),
c     &                 MPI_COMM_WORLD,IERR)
c         IF (NTASK .EQ. 1) THEN
c            CALL PARION
c            RDINPC = .FALSE.
c            CALL READIN(WORK,LWORK,.FALSE.)
c            GOTO 100
c         ELSE IF (NTASK .EQ. 0) THEN
c            CALL MPI_BCAST(FINISH,1,MPI_LOGICAL,NODEID(MASTER),
c     &                  MPI_COMM_WORLD,IERR)
c            CALL MPI_FINALIZE(IERR)
c#if defined (SYS_CRAY)
c            CALL INFO = ISHELL('rm -f $SCRATCHDIR/*')
c#else
c            CALL SYSTEM('rm -f $SCRATCHDIR/*')
c#endif
c            STOP '*** End of DALTON calculation ***'
c         ELSE
c            WRITE (LUPRI,'(/A)') 'Unknown message received by slave'
c            CALL QUIT('Slave received unknown message from master')
c         END IF
c      END IF
c#endif
C
C     *** No more numerical derivative. ***
C
      NMWALK = .FALSE.
      IPRUSR = IPRUSR_orig
      IPREAD = IPREAD_orig
C
      CALL QEXIT('NUMDRV')
      RETURN
      END
C
C  /* Deck numdr1 */
      SUBROUTINE NUMDR1(DERIV,FUVAL,COOR,CSTART,SYMCOR,DCOEFF,TSTGDR,
     &                  TSTSDR,ENTST,GRIREP,CHRCTR,TRNCCR,TRAMSS,DKIN,
     &                  FREQ,RNNORM,WORK,KDPMTX,NMIDPC,ICRIRP,ISOTMP,
     &                  MBKLIN,LWORK,NDERIV,NDIME,NINTIN,MAXADR,
     &                  LTHTST,LDPMTX,IFRSTD,IPRINT,WRKDLM)
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <infpar.h>
      PARAMETER (D0 = 0.0D0)
#include <cbiwlk.h>
#include <cbinum.h>
#include <nuclei.h>
#include <symmet.h>
#include <exeinf.h>
#include <abainf.h>
#include <taymol.h>
#include <trkoor.h>
#include <numder.h>
#include <prpndr.h>
#include <past.h>
#include <gnrinf.h>
#include <inftap.h>
#include <molinp.h>
#include <fcsym.h>
      LOGICAL EXHER, EXSIR, EXABA, GTCOMP, RSTDON, SYMDET, NPRBKP,
     &        NSPNBK
      CHARACTER*80 MBKLIN, PRTEXT
      CHARACTER*8 ANDER
      CHARACTER*6 TXT
      DIMENSION MBKLIN(NMLINE)
      DIMENSION DERIV(NDERIV), FUVAL(NINTIN,NDIME), CSTART(3*NCOOR),
     &          COOR(3*NCOOR), SYMCOR(NCOOR,NCOOR),
     &          DCOEFF(LDPMTX,IFRSTD), TSTGDR(NCOOR), ENTST(LTHTST),
     &          TSTSDR(NCOOR,NCOOR), GRIREP(NGORDR,NGVERT),
     &          CHRCTR(NGORDR,NCVERT), TRNCCR(NCOOR,NCOOR),
     &          TRAMSS(NCOOR), DKIN(NCOOR), FREQ(NCOOR), RNNORM(NCOOR),
     &          ISOTMP(NATOMS), WORK(LWORK)
      DIMENSION ICRIRP(NCOOR,2), KDPMTX(LDPMTX,NSTRDR,IFRSTD),
     &          NMIDPC(LDPMTX)
#include <ibtfun.h>
C
C     ******************************
C     *** Restart initialization ***
C     ******************************
C
      RSTDON = .FALSE.
C
C     *****************************************
C     *** Backing up symmetry for later use ***
C     *****************************************
C
      CALL BKSMNM
c      CALL DALCHG(DUMMY,IDUMMY,IDUMMY,IPRINT,IDUMMY,IDUMMY,.TRUE.)
C
C     *************************************
C     ***Backup of original MOLECULE.INP***
C     ***   To finish off correctly.    ***
C     *************************************
C
      NMBKLN = NMLINE
      DO 100 IMLINE = 1, NMLINE
         MBKLIN(IMLINE) = MLINE(IMLINE)
 100  CONTINUE
C
C     *********************************************
C     *** Backup of isotopes, if abacus is run. ***
C     *********************************************
C
      IF (NAORDR.GT.0) THEN
         CALL ICOPY(NATOMS,ISOTOP,1,ISOTMP,1)
      END IF
C
C     ****************************************
C     ***Unrolling the symmetry coordinates***
C     ***  In order to take proper steps   ***
C     ****************************************
C
      ICOOR = 0
      IATOM = 0
      DO 200 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         DO 300 IOP = 0, MAXOPR
            IF (IBTAND(IOP,MULCNT) .EQ. 0) THEN
               IATOM = IATOM + 1
               DO 400 I = 1, 3
                  ICOOR = ICOOR + 1
                  CSTART(ICOOR) =
     &                 PT(IBTAND(ISYMAX(I,1),IOP))*CORD(I,ICENT)
 400           CONTINUE
            END IF
 300     CONTINUE
 200  CONTINUE
C
C     *** Restart ***
C
      LURSTR = -1
      CALL GPOPEN(LURSTR,'RSTRT.FC','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LURSTR)
      IF (NPRPDR) THEN
         LUNDPR = -1
         CALL GPOPEN(LUNDPR,'PROPERTY.NDER','UNKNOWN',' ','FORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND(LUNDPR)
      END IF
C
      IF (RESTRT) THEN
C
C        *** Restart, find which round it ended  ***
C        *** in this run of nmder. Reread fuval  ***
C        *** values from file.                   ***
C
         KEND = 0
         IDIMAX = 0
         IDIMIN = 2
         CALL RERSTR(FUVAL,SYMCOR,DUMMY,DUMMY,ICRIRP,NDIME,NINTIN,
     &               KEND,IDIMAX,IDIMIN,LURSTR,IPRINT,RSTDON)
C
C        *** Restart for property derivatives. ***
C
         IF (NPRPDR.AND..NOT.NRMCRD) THEN
            KSTRT = 2
            LWRK1 = LWORK - KSTRT
            CALL PRPRER(WORK(KSTRT),IDIMAX,IDIMIN,LURSTR,LWRK1)
         END IF
      ELSE
C
C        *** Open restart file. A zero first in the file ***
C        *** means that the calculation ended here.      ***
C
         WRITE (LURSTR,'(I2)') 0
      END IF
C
      NUMCAL = 0
      SYMDET = .TRUE.
      NDCOOR = NCOOR
      IF (NRMCRD) THEN
         NTMPDR = NMORDR
         NMORDR = 2 - NAORDR
         IF ((NAORDR+NMORDR).LT.2) NMORDR = NTMPDR
      END IF
      MAXINR = 2**NMORDR
C
      IF (.NOT.(PREHES.AND.NRMCRD)) THEN
         IF (PREHES) WRITE (LUPRI,'(/A)') '   Not able to use a' //
     &           ' precalculated hessian, since normal coordinates' //
     &           ' are not specified.'
C
C        *** If normal coordinates then no property derivatives ***
C        *** should be calculated at this time.                 ***
C
         IF (NRMCRD) THEN
C
C           *** First time through. ***
C
            FRSTNM = .TRUE.
C
C           *** If normal coordinates then no  ***
C           *** property derivatives should be ***
C           *** calculated at this time.       ***
C
            NPRBKP = NPRPDR
            NPRPDR = .FALSE.
            NSPNBK = NSPNSP
            NSPNSP = .FALSE.
C
         END IF
C
         KIADRS = 2
         KINDST = KIADRS + MAXADR
         KNPRTN = KINDST + NMORDR
         KINDTM = KNPRTN + NMORDR
         NTORDR = NMORDR
C        ... NTORDR is used for DIMENSION in NMDER /hjaaj
         KIDCMP = KINDTM + MAXINR
         KIEQVG = KIDCMP + NCOOR
         KICIN  = KIEQVG + 2*NMORDR
         KIRPID = KICIN  +   NMORDR
         KEGRAD = KIRPID +   NMORDR
         KEHESS = KEGRAD +   MXCOOR
c#if defined (VAR_MPI)
c         KFTVAL = KEHESS +   MXCOOR**2
c         KLAST  = KFTVAL +   NDIME
c#else
         KLAST  = KEHESS +   MXCOOR**2
c#endif
C
         LWRK  = LWORK - KLAST + 1
C
         CALL NMDER(DERIV,FUVAL,COOR,CSTART,SYMCOR,ENTST,DCOEFF,GRIREP,
     &        CHRCTR,WORK(KEGRAD),WORK(KEHESS),WORK(KLAST),ICRIRP,
     &        KDPMTX,NMIDPC,WORK(KIADRS),WORK(KINDST),WORK(KINDTM),
     &        WORK(KNPRTN),WORK(KIDCMP),WORK(KIEQVG),WORK(KICIN),
     &        WORK(KIRPID),MBKLIN,NMBKLN,NDERIV,LWRK,NDIME,NINTIN,
     &        MAXADR,LDPMTX,IFRSTD,MAXINR,LTHTST,IDIMAX,IDIMIN,
     &        LURSTR,IPRINT,WRKDLM,
c#if defined (VAR_MPI)
c     &        WORK(KFTVAL),SYMDET,RSTDON)
c#else
     &        SYMDET,RSTDON)
c#endif
C
C        *** Restoring property derivatives. ***
C
         IF (NRMCRD) THEN
            NPRPDR = NPRBKP
            NSPNSP = NSPNBK
         END IF
C
C        *******************************************
C        *** Restore isotopes, if abacus is run. ***
C        *******************************************
C
         IF (NAORDR.GT.0) THEN
            CALL ICOPY(NATOMS,ISOTMP,1,ISOTOP,1)
         END IF

      ELSE
         WRITE (LUPRI,'(A)') '                                '
         WRITE (LUPRI,'(A)') 'Reading in precalculated hessian'
         WRITE (LUPRI,'(A)') '                                '
C
         KHSSIN = 1
         KLAST  = KHSSIN + NCOOR**2
         LWRK   = LWORK - KLAST + 1
         CALL RDHESS(SYMCOR,CSTART,GRIREP,CHRCTR,WORK(KHSSIN),
     &               WORK(KLAST),ICRIRP,LWRK,IPRINT,SYMDET)
      END IF
C
C     *** Resetting some variables. ***
C
      IF (NRMCRD) THEN
         NMORDR = NTMPDR
      END IF
C
C     *** Close the restart file ***
C
      CALL GPCLOSE(LURSTR,'KEEP')
C
C     *** Printing the derivatives ***
C
      IF (MYNUM.EQ.0) THEN
         LTEXT  = 8
         PRTEXT(1:8) = 'symmetry'
         IF (NRMCRD) THEN
            NPRRDR = 2
         ELSE
            NPRRDR = NMORDR+NAORDR
         END IF
         NDIMT  = NCOOR*(NCOOR+1)*(NCOOR+2)/6
         NDIMF  = NCOOR*(NCOOR+1)*(NCOOR+2)*(NCOOR+3)/24
         KTDRS  = 1
         KFDRS  = KTDRS + NDIMT
C
         KTTMPD = 2
         KFTMPD = KTTMPD + NCOOR**3
         KLAST  = KFTMPD + NCOOR**4
         LWRK   = LWORK - KLAST + 1
C
         CALL PRDERV(DERIV(KTDRS),DERIV(KFDRS),TSTGDR,TSTSDR,SYMCOR,
     &       CSTART,WORK(KTTMPD),WORK(KFTMPD),RNNORM,WORK(KLAST),ICRIRP,
     &       LWRK,NPRRDR,NDIMT,NDIMF,LTEXT,IPRINT,PRTEXT)
C
         IF (SDRTST) THEN
            KTMPGD = 2
            KTMPHS = KTMPGD + MXCOOR
            KBKPGD = KTMPHS + NCOOR**2
            KBKPHS = KBKPGD + NCOOR
            KLAST  = KBKPHS + NCOOR**2
            LWRK = LWORK - KLAST + 1
            IF (KLAST.GT.LWORK) CALL QUIT('Memory exceeded in SDERTT')
            CALL SDERTT(TSTSDR,TSTGDR,SYMCOR,WORK(KTMPGD),WORK(KTMPHS),
     &           WORK(KBKPGD),WORK(KBKPHS),WORK(KLAST),LWRK,
     &           WRKDLM,IPRINT)
         END IF
C
C     *** Reevaluate restart parameter ***
C
         IF (RESTRT.AND.RSTDON) THEN
            RESTRT = .FALSE.
            RSTDON = .FALSE.
         END IF
C
C
         IF ((NAORDR+NMORDR).GT.1) THEN
            IF (.NOT.DRYRUN) THEN
               IF (HARMON) THEN
                  KBKPHS = 2
                  KTRAMT = KBKPHS + MXCOOR**2
                  KTMPHS = KTRAMT + NCOOR**2
                  KLAST  = KTMPHS + NCOOR**2
                  LWRK = LWORK - KLAST + 1
                  CALL HARMAN(SYMCOR,WORK(KBKPHS),WORK(KTRAMT),
     &                        WORK(KTMPHS),WORK(KLAST),NCOOR,LWRK,
     &                        IPRINT)
               END IF
C
               IF (PREHES.AND..NOT.NRMCRD) THEN
                  KHSSIN = 1
                  KLAST  = KHSSIN + NCOOR**2
                  LWRK   = LWORK - KLAST + 1
                  CALL RDHESS(SYMCOR,CSTART,GRIREP,CHRCTR,WORK(KHSSIN),
     &                        WORK(KLAST),ICRIRP,LWRK,IPRINT,SYMDET)
               END IF
C
C              *** If we are not doing things in normal coordinates ***
C              *** we need to save the symmetry coordinates.        ***
C
               IF (.NOT. NRMCRD) THEN
                  KSYCAR = 2
                  KEIGNV = KSYCAR + NCOOR**2
                  CALL DCOPY(NCOOR**2,SYMCOR,1,WORK(KSYCAR),1)
               ELSE
                  KEIGNV = 2
               END IF
C
               IF (.NOT.RESTRT) THEN
                  KEGNVC = KEIGNV + NCOOR
                  KHSMWT = KEGNVC + NCOOR**2
                  KMT1TP = KHSMWT + NCOOR*(NCOOR+1)/2
                  KMT2TP = KMT1TP + NCOOR**2
                  KAMASS = KMT2TP + NCOOR**2
                  KHTSTM = KAMASS + NATOMS
                  KNATYP = KHTSTM + NCOOR**2
                  KNMSSP = KNATYP + NATOMS
                  KCRTMP = KNMSSP + NCOOR
                  KLAST  = KCRTMP + NCOOR
C
                  LWRK   = LWORK  - KLAST + 1
                  IF (KLAST.GT.LWORK)
     &                        CALL QUIT('Memory exceeded in MKNRMC')
                  CALL MKNRMC(SYMCOR,CSTART,TRNCCR,TRAMSS,WORK(KEIGNV),
     &                 WORK(KEGNVC),WORK(KHSMWT),WORK(KMT1TP),
     &                 WORK(KMT2TP),WORK(KAMASS),DKIN,WORK(KHTSTM),FREQ,
     &                 RNNORM,WORK(KCRTMP),WORK(KLAST),ICRIRP,
     &                 WORK(KNATYP),WORK(KNMSSP),LWRK,IPRINT)
               END IF
C
C              *** Debug printing. ***
C
               IF (.NOT.NRMCRD.AND.(IPRINT.GE.50)) THEN
                  NDIMT  = NCOOR*(NCOOR+1)*(NCOOR+2)/6
                  NDIMF  = NCOOR*(NCOOR+1)*(NCOOR+2)*(NCOOR+3)/24
                  KTDRS  = 1
                  KFDRS  = KTDRS + NDIMT
C
                  KSYCAR = 2
                  KHSNRM = KSYCAR + NCOOR**2
                  KCRTNM = KHSNRM + NCOOR**2
                  KTNRMD = KCRTNM + NCOOR**2
                  KFNRMD = KTNRMD + NCOOR**3
                  KLAST  = KFNRMD + NCOOR**4
                  LWRK   = LWORK - KLAST + 1
C
                  CALL TRAFRC(DERIV(KTDRS),DERIV(KFDRS),WORK(KHSNRM),
     &                        SYMCOR,WORK(KCRTNM),WORK(KSYCAR),
     &                        WORK(KTNRMD),WORK(KFNRMD),WORK(KLAST),
     &                        NCOOR,NDIMF,NDIMT,LWRK,IPRINT)
               END IF
            ELSE
               CALL DRNRMC(SYMCOR,ICRIRP,IPRINT)
            END IF
         END IF
      END IF
C
      IF (NRMCRD .AND. ((NAORDR+NMORDR).GT.2)) THEN
         FRSTNM = .FALSE.
C
C        *** Restart ***
C
         LURSTR = -1
         CALL GPOPEN(LURSTR,'RSTRT.FC','UNKNOWN',' ','FORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LURSTR)
C
         IF (RESTRT) THEN
C
C           *** Restart, find which round it ended  ***
C           *** in this run of nmder. Reread fuval  ***
C           *** values and normal coordinates from  ***
C           *** file.                               ***
C
            KEND   = 1
            IDIMAX = 0
            IDIMIN = 2
            CALL RERSTR(FUVAL,SYMCOR,RNNORM,FREQ,ICRIRP,NDIME,NINTIN,
     &                  KEND,IDIMAX,IDIMIN,LURSTR,IPRINT,RSTDON)
C
C           *** Restart for property derivatives. ***
C
            IF (NPRPDR) THEN
               KSTRT = 2
               LWRK1 = LWORK - KSTRT
               CALL PRPRER(WORK(KSTRT),IDIMAX,IDIMIN,LURSTR,LWRK1)
            END IF
C
C           *** Writing to spectro file if requested. ***
C
            IF (SPECTR) THEN
               NTIME = 1
               IF (NRMCRD) THEN
                  TXT  = 'normal'
               ELSE
                  TXT  = 'cartes'
               END IF
               CALL WRISPC(FREQ,RNNORM,DUMMY,DUMMY,TXT,NCOOR,NDCOOR,
     &                     NTIME,IPRINT)
            END IF
            IF (MIDAS) THEN
               NTIME = 1
               IF (NRMCRD) THEN
                  TXT  = 'normal'
               ELSE
                  TXT  = 'cartes'
               END IF
               IF (NRMCRD) CALL WRIMOP(FREQ,RNNORM,DUMMY,DUMMY,TXT,
     &                                 NCOOR,NDCOOR,NTIME,IPRINT)
            END IF
         ELSE
C
C           *** Open restart file. A 1 first in the file ***
C           *** means that the calculation ended here.   ***
C           *** Normal coordinates are also written to   ***
C           *** file for restart purposes.               ***
C
            WRITE (LURSTR,'(I2)') 1
            CALL WRICOR(SYMCOR,RNNORM,FREQ,ICRIRP,LURSTR,IPRINT)
         END IF
C
         MAXINR = 2**NMORDR
         NTORDR = NMORDR+NAORDR
C
         KIADRS = 2
         KINDST = KIADRS + MAXADR
         KNPRTN = KINDST + NTORDR
         KINDTM = KNPRTN + NTORDR
         KIDCMP = KINDTM + MAXINR
         KIEQVG = KIDCMP + NCOOR
         KICIN  = KIEQVG + 2*NMORDR
         KIRPID = KICIN  +   NMORDR
         KEGRAD = KIRPID +   NMORDR
         KEHESS = KEGRAD +   MXCOOR
c#if defined (VAR_MPI)
c         KFTVAL = KEHESS +   MXCOOR**2
c         KLAST  = KFTVAL +   NDIME
c#else
         KLAST  = KEHESS +   MXCOOR**2
c#endif

C
         LWRK   = LWORK - KLAST + 1
C
         CALL NMDER(DERIV,FUVAL,COOR,CSTART,SYMCOR,ENTST,DCOEFF,GRIREP,
     &              CHRCTR,WORK(KEGRAD),WORK(KEHESS),WORK(KLAST),ICRIRP,
     &              KDPMTX,NMIDPC,WORK(KIADRS),WORK(KINDST),
     &              WORK(KINDTM),WORK(KNPRTN),WORK(KIDCMP),WORK(KIEQVG),
     &              WORK(KICIN),WORK(KIRPID),MBKLIN,NMBKLN,NDERIV,LWRK,
     &              NDIME,NINTIN,MAXADR,LDPMTX,IFRSTD,MAXINR,
     &              LTHTST,IDIMAX,IDIMIN,LURSTR,IPRINT,WRKDLM,
c#if defined (VAR_MPI)
c     &              WORK(KFTVAL),SYMDET,RSTDON)
c#else
     &              SYMDET,RSTDON)
c#endif
C
C        *** Close the restart file ***
C
         CALL GPCLOSE(LURSTR,'KEEP')
C
C        *** Reevaluate restart parameter ***
C
         IF (RESTRT.AND.RSTDON) THEN
            RESTRT = .FALSE.
         END IF
C
         IF (MYNUM.EQ.0) THEN
            LTEXT  = 6
            PRTEXT(1:6) = 'normal'
C
            NPRRDR = NAORDR + NMORDR
            NDIMT  = NDCOOR*(NDCOOR+1)*(NDCOOR+2)/6
            NDIMF  = NDCOOR*(NDCOOR+1)*(NDCOOR+2)*(NCOOR+3)/24
C
            KTDRS  = 1
            KFDRS  = KTDRS + NDIMT
C
            KTTMPD = 2
            KFTMPD = KTTMPD + NCOOR**3
            KLAST  = KFTMPD + NCOOR**4
            LWRK   = LWORK - KLAST + 1
C
            CALL PRDERV(DERIV(KTDRS),DERIV(KFDRS),TSTGDR,TSTSDR,SYMCOR,
     &           CSTART,WORK(KTTMPD),WORK(KFTMPD),RNNORM,WORK(KLAST),
     &           ICRIRP,LWRK,NPRRDR,NDIMT,NDIMF,LTEXT,IPRINT,PRTEXT)
C
C           *** Isotope analysis. ***
C
c            KTPGRD = 2
c            KTPHES = KTPGRD + NCOOR
c            KTPMSS = KTPHES + NCOOR**2
c            KTPTD1 = KTPMSS + NCOOR
c            KTPTD2 = KTPTD1 + NCOOR**3
c            KLAST  = KTPTD2 + NCOOR**3
c            LWRK   = LWORK -KLAST + 1
c            IF (KLAST.GT.LWORK)CALL QUIT('Memory exceeded in NRMISO')
c            CALL NRMISO(DERIV(KTDRS),SYMCOR,DKIN,TRNCCR,TRAMSS,
c     &                  WORK(KTPGRD),WORK(KTPHES),WORK(KTPMSS),
c     &                  WORK(KTPTD1),WORK(KTPTD2),CSTART,WORK(KLAST),
c     &                  NDIMT,LWRK)
         END IF
      END IF
C
      IF (MYNUM.EQ.0) THEN
         IF (NAORDR.EQ.0) THEN
            ANDER = 'energy  '
         ELSE IF (NAORDR .EQ. 1) THEN
            ANDER = 'gradient'
         ELSE
            ANDER = 'hessian '
         END IF
C
         IF (DRYRUN) THEN
            WRITE (LUPRI,'(//5X,A,I10)') 'Number of ' // ANDER //
     &              ' calculations needed:', NUMCAL
         ELSE
            WRITE (LUPRI,'(//A,I10)') '@>>> Number of ' // ANDER //
     &        ' calculations done:', NUMCAL
         END IF

         CALL TITLER('Numerical derivatives have now been calculated.',
     &      '*',118)

C
         IF (ANALZ1) THEN
            CALL PRIPRP
            CALL NVBDRV(DERIV,SYMCOR,FREQ,RNNORM,CSTART,WORK,LWORK,
     &                  NDERIV,IPRINT)
         END IF
      END IF
C
      RETURN
      END
C
C  /* Deck nmder */
      SUBROUTINE NMDER(DERIV,FUVAL,COOR,CSTART,SYMCOR,ENTST,DCOEFF,
     &                 GRIREP,CHRCTR,EGRAD,EHESS,WORK,ICRIRP,KDPMTX,
     &                 NMIDPC,IADRSS,INDSTP,INDTMP,NPRTNR,IDCOMP,
     &                 IEQVGM,ICIN,IRPIND,MBKLIN,NMBKLN,NDERIV,LWORK,
     &                 NDIME,NINTIN,MAXADR,LDPMTX,IFRSTD,MAXINR,
     &                 LTHTST,IDIMAX,IDIMIN,LURSTR,IPRINT,WRKDLM,
c#if defined (VAR_MPI)
c     &                 FTVAL,SYMDET,RSTDON)
c#else
     &                 SYMDET,RSTDON)
c#endif
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D1 = 1.0D0, DM1 = -1.0D0, D0 = 0.0D0, DMAX=1.0D-15)
c#if defined (VAR_MPI)
c      INCLUDE 'mpif.h'
c      DIMENSION my_STATUS(MPI_STATUS_SIZE)
c      DIMENSION FTVAL(NINTIN,NDIME)
c#endif
c#if defined (VAR_MPI2)
cC
cC     pario.h will no longer be needed as an include file when locking RMA
cC     operations become available.
cC
c#include <dummy.h>
c#include <pario.h>
c#endif
#include <infpar.h>
#include <inforb.h>
#include <cbiexc.h>
C
#include <trkoor.h>
#include <nuclei.h>
#include <numder.h>
#include <molinp.h>
#include <cbiwlk.h>
#include <cbinum.h>
#include <fcsym.h>
#include <abainf.h>
#include <taymol.h>
#include <pvibav.h>
#include <prpc.h>
#include <gnrinf.h>
C
      LOGICAL CALCMP, CLNRGY, LASTE, NOSYM, TOTSM, SYMDET, RSTDON,
     &        ALRCAL, PRTNR, CPRPBK
      CHARACTER*80 MBKLIN
      DIMENSION MBKLIN(NMBKLN)
      DIMENSION FUVAL(NINTIN,NDIME), COOR(NCOOR), CSTART(NCOOR),
     &          DERIV(NDERIV), DCOEFF(LDPMTX,IFRSTD),
     &          SYMCOR(NCOOR,NCOOR), GRIREP(NGORDR,NGVERT),
     &          CHRCTR(NGORDR,NCVERT), EGRAD(MXCOOR),
     &          EHESS(MXCOOR,MXCOOR), ENTST(LTHTST), WORK(LWORK)
      DIMENSION INDSTP(NTORDR), INDTMP(NTORDR), ICIN(NMORDR),
     &          IDCOMP(NCOOR), IADRSS(MAXADR), IRPIND(NMORDR),
     &          ICRIRP(NCOOR,2), KDPMTX(LDPMTX,NSTRDR,IFRSTD),
     &          NMIDPC(LDPMTX), IEQVGM(NMORDR,2), NPRTNR(MAXINR)
C
      CALL QENTER('NMDER')
c#if defined (VAR_MPI2)
c      LUNMCL = -9056
c#endif
C
C     *** Finite field header. ***
C
      CALL TITLER('Numerical derivatives.','*',118)
      CALL HEADER('Derivatives calculated:',0)
      WRITE (LUPRI,'(A,I3)')
     &     '           Derivatives calculated to order', NMORDR + NAORDR
      WRITE (LUPRI,'(A,I3,A)')
     &     '           Analytical derivatives from energies to ',
     &     NAORDR, ' order.'
      WRITE (LUPRI,'(A,I3,A,I3,A)')
     &     '         ', NMORDR, '. numerical derivatives from',NAORDR,
     &     '. order analytical derivatives'
C
C     *** Symmetry adapted coordinates ***
C
      IF (SYMDET) THEN
         CALL GRPCHR(CSTART,SYMCOR,GRIREP,CHRCTR,WORK,ICRIRP,LWORK,
     &               IPRINT)
         SYMDET = .FALSE.
      END IF
C
C     *** Finding force constants that are dependent on each-other ***
C
      NLDPMX = 0
      KDIM = IFRSTD*NSTRDR*LDPMTX
      CALL IZERO(KDPMTX,KDIM)
      CALL FSDCST(SYMCOR,GRIREP,DCOEFF,WORK,KDPMTX,NMIDPC,ICRIRP,LDPMTX,
     &            IFRSTD,NLDPMX,LWORK,IPRINT)
C
C     *** Memory allocations for future use ***
C
      NTYPE  = 3
      LASTE  = .FALSE.
      NOSYM  = .FALSE.
      CALL UPCASE(MLINE(NMLAU))
      IPOS = INDEX(MLINE(NMLAU),'ATO')
      IF (IPOS .EQ. 0) THEN
         IF (MLINE(NMLAU)(10:10).EQ.'0') NOSYM = .TRUE.
      ELSE
         IPOS = INDEX(MLINE(NMLAU),'NOS')
         IF (IPOS .NE. 0) NOSYM = .TRUE.
      END IF
C
      KIDTMP = 1
      KIRPDG = KIDTMP + NMORDR
      KIRPST = KIRPDG + NMORDR
C
      ITYPE = 1
C
      IDIME = 2
      IF (PGMTST) EMAX = D0
C
C     *** Order for derivatives ***
C
      DO 100 IORDR  = 1, NMORDR
         IHORDR = INT((IORDR+1)/2)
C
C        *** IHORDR -> The maximum order in one direction for this ***
C        ***           numerical derivative.                       ***
C
         DO 200 IMXRDR = 1, IHORDR, 1
C
            ITYPE         = ITYPE + 1
            IADRSS(ITYPE) = IDIME - 1
C
C           ***   The first component IX1 has always the largest order,   ***
C           ***          and are then independent of the others           ***
C           *** The order of the other components are not larger then one ***
C
            IMINCR = 1
            IF (IMXRDR .EQ. 1) IMINCR = IORDR
            IRSRDR    = IORDR - (2*IMXRDR-1)
            DO 300 IX1    = IMINCR, NDCOOR
C
C              *** Starting values for the component-vector. ***
C
               INDSTP(1) = IX1
               DO 400 IC = IRSRDR+1, 2, -1
                  INDSTP(IC) = IRSRDR+2-IC
 400           CONTINUE
               IF (IRSRDR .GT. 0) INDSTP(IRSRDR+1) = INDSTP(IRSRDR+1)-1
C
               NSTP = 1
               IF (IMXRDR .EQ. 1) THEN
                  DO 500 I = 1, IRSRDR
                     NSTP = NSTP*(IX1-I)/I
 500              CONTINUE
               ELSE
                  DO 600 I = 1, IRSRDR
                     NSTP = NSTP*(NDCOOR-I+1)/I
 600              CONTINUE
               END IF
C
C              *** NSTP -> Number of components for this IORDR, IMXRDR and IX1***
C
               DO 700 ISTP = 1, NSTP
C
C                 *** Finding the other components. ***
C
                  CALCMP = .TRUE.
                  IF (IMXRDR .EQ. 1) THEN
                     DO 800 IC = IRSRDR+1, 2, -1
                        IF (INDSTP(IC) .LT. (INDSTP(IC-1)-1)) THEN
                           INDSTP(IC) = INDSTP(IC) + 1
                           DO 900 I = IC+1, IRSRDR+1
                              INDSTP(I) = (IRSRDR+2) - I
 900                       CONTINUE
                           GOTO 1300
                        END IF
 800                 CONTINUE
                  ELSE
                     DO 1000 IC = IRSRDR+1, 2, -1
                        IF (IC .EQ. 2) THEN
                           INDSTP(2) = INDSTP(2) + 1
                           DO 1100 I = 3, IRSRDR+1
                              INDSTP(I) = 1
 1100                      CONTINUE
                           DO 1150 ICN = 2, IRSRDR+1
                              IF (INDSTP(ICN) .EQ. IX1) CALCMP = .FALSE.
 1150                      CONTINUE
                           GOTO 1300
                        ELSE IF (INDSTP(IC) .LT. (INDSTP(IC-1)-1)) THEN
                           INDSTP(IC) = INDSTP(IC) + 1
                           DO 1200 I = IC+1, IRSRDR
                              INDSTP(I) = 1
 1200                      CONTINUE
                           DO 1250 ICN = 2, IRSRDR+1
                              IF (INDSTP(ICN) .EQ. IX1) CALCMP = .FALSE.
 1250                      CONTINUE
                           GOTO 1300
                        END IF
 1000                CONTINUE
                  END IF
C
 1300             CONTINUE
C
C                 *** Have we calculated this function-value before? ***
C
                  IF (CALCMP) THEN
C
C                    *** IDCOMP(INDSTP(IC)) -> gives the length of the ***
C                    *** steps we need to do in INDSTP(IC) direction.  ***
C
                     CALL IZERO(IDCOMP,NDCOOR)
                     IDCOMP(INDSTP(1)) = IMXRDR
                     DO 1400 IC = 2, IRSRDR+1
                        IDCOMP(INDSTP(IC)) = IDCOMP(INDSTP(IC)) + 1
 1400                CONTINUE
C
C                    *** NINNER -> Number of different steps needed ***
C                    *** If numerical derivatives from energy is    ***
C                    *** calculated, we need to check whether the   ***
C                    *** steps are all totally symmetric.           ***
C
                     IF (NAORDR .EQ. 0) THEN
                        IJ = 1
                        TOTSM = .FALSE.
                        DO IRDR = 1, IRSRDR+1
                           IJ = IJ*ICRIRP(INDSTP(IRDR),1)
                        END DO
                        IF (IJ.EQ.1) TOTSM = .TRUE.
                     END IF
C
                     NINNER = 2**(IRSRDR+1)
C
                     NMPRTN = 0
                     CALL IZERO(NPRTNR,MAXINR)
C
                     DO 1500 IINNER = 1, NINNER
C
C                    *** Initialize alrcal. ***
C
                        ALRCAL = .FALSE.
C
C                       *** Finding the appropriate step-possibility. ***
C
                        IC   = 0
                        IMOD = 1
                        IDIV = 1
                        DO 1600 I = 1, IRSRDR+1
                           ICIN(I) = MOD(INT((IINNER-1)/IDIV),2) + 1
                           IDIV = IDIV*2
 1600                   CONTINUE
C
C                       *** Making the appropriate step, and  ***
C                       ***      get the function value.      ***
C
                        CALL GTNPNT(FUVAL,GRIREP,SYMCOR,EGRAD,EHESS,
     &                       COOR,CSTART,WORK,WRKDLM,INDSTP,ICRIRP,
     &                       NPRTNR,ICIN,KDPMTX,IRPIND,IDCOMP,LDPMTX,
     &                       IFRSTD,NLDPMX,IORDR,IRSRDR,IINNER,NMPRTN,
     &                       NDIME,MAXINR,LWORK,NMBKLN,MBKLIN,IDIME,
     &                       NINTIN,IDIMAX,IDIMIN,LURSTR,
c#if defined (VAR_MPI)
c     &                       FTVAL,RSTDON,PRTNR,ALRCAL,CLNRGY,LASTE)
c#else
     &                       RSTDON,PRTNR,ALRCAL,CLNRGY,LASTE)
c#endif
C
 1500                CONTINUE
                     IDIME = IDIME + NMPRTN
                  END IF
 700           CONTINUE
C
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
C
C     *** Returning to the original geometry, only for master. ***
C
      IF (MYNUM .EQ. 0) THEN
         IORDR = 0
         IDIME = 1
         CLNRGY = .TRUE.
         MINLIM = 1
         IADRSS(1) = 0
         LASTE = .TRUE.
         IF (NAORDR.EQ.1) MINLIM = NCOOR
         CALL GTNPNT(FUVAL,GRIREP,SYMCOR,EGRAD,EHESS,COOR,CSTART,WORK,
     &        WRKDLM,INDSTP,ICRIRP,NPRTNR,ICIN,KDPMTX,IRPIND,IDCOMP,
     &        LDPMTX,IFRSTD,NLDPMX,0,IRSRDR,IINNER,NMPRTN,NDIME,
     &        MAXINR,LWORK,NMBKLN,MBKLIN,IDIME,NINTIN,IDIMAX,IDIMIN,
c#if defined (VAR_MPI)
c     &        LURSTR,FTVAL,RSTDON,PRTNR,ALRCAL,CLNRGY,LASTE)
c#else
     &        LURSTR,RSTDON,PRTNR,ALRCAL,CLNRGY,LASTE)
c#endif
C
c#if defined (VAR_MPI)
c      ELSE
c         DO IDIME = 1, NDIME
c         DO INTIN = 1, NINTIN
c            FUVAL(INTIN,IDIME) = D0
c         END DO
c         END DO
c#endif
      END IF
C
c#if defined (VAR_MPI)
cC
cC     *** If paralell calculation all energies are ***
cC     *** collected into one array.                ***
cC
c#if defined (VAR_MPI2)
cC
cC     However, in the case of "simulated" MPI2 behaviour (RMA operations),
cC     there can occur a "glitch" in the NFS lock file, and points may happen
cC     to be calculated on several processors.
cC
cC     We collect results from one processor at a time, checking for double
cC     counting
cC
c      IF (MYNUM .EQ. 0) THEN
c         CALL DCOPY(NDIME,FTVAL,1,FUVAL,1)
c         DO IWHO = 1, NODTOT
c            CALL MPI_RECV(NWHO,1,my_MPI_INTEGER,MPI_ANY_SOURCE,65,
c     &                    MPI_COMM_WORLD,ISTAT,IERR)
c            CALL MPI_RECV(FTVAL,NDIME,MPI_DOUBLE_PRECISION,NWHO,65,
c     &                    MPI_COMM_WORLD,ISTAT,IERR)
c            DO IPOS = 1, NDIME
cC
cC     Molecular energies ought to be negative
cC
c               IF (.NOT. (FUVAL(IPOS) .LT. D0))
c     &              FUVAL(IPOS) = FTVAL(IPOS)
c            END DO
c         END DO
c      ELSE
c         CALL MPI_SEND(MYNUM,1,my_MPI_INTEGER,NODEID(MASTER),
c     &                 65,MPI_COMM_WORLD,IERR)
c         CALL MPI_SEND(FTVAL,NDIME,MPI_DOUBLE_PRECISION,NODEID(MASTER),
c     &                 65,MPI_COMM_WORLD,IERR)
c      END IF
c#else
c      CALL MPI_REDUCE(FTVAL,FUVAL,NINTIN*NDIME,MPI_DOUBLE_PRECISION,
c     &                MPI_SUM,0,MPI_COMM_WORLD,IERR)
c#endif
c#endif
C
      IF (MYNUM.EQ.0) THEN
C
C        ********************************************
C        *** Preliminary constants to derivatives ***
C        ********************************************
C
         MXCOEF = INT(NMORDR/2) + 1
C
C        ****************************************************
C        *** Calculating force field. Property derivative ***
C        *** needs to be reset for call for NMNDER.       ***
C        ****************************************************
C
         CPRPBK = CNMPRP
         CNMPRP = .FALSE.
C
C        ***********************************
C        *** Calculating the derivatives ***
C        ***********************************
C
         NFINNR = 1
         KCOEF  = 1
         KIMAX  = KCOEF  + (2*MXCOEF+1)*(NMORDR+1)
         KIMIN  = KIMAX  +               NMORDR
         KICNT  = KIMIN  +               NMORDR
         KNCVAL = KICNT  +               NTYPE
         KIDDCP = KNCVAL +               NCOOR
         KLAST  = KIDDCP +               NCOOR
         LWRK1  = LWORK - KLAST
         IF (LWRK1.LT.1) CALL QUIT('Memory exceeded in NMNDER')
         CALL NMNDER(DERIV,WORK(KCOEF),FUVAL,GRIREP,WORK(KLAST),IADRSS,
     &        KDPMTX,ICRIRP,INDSTP,INDTMP,IDCOMP,WORK(KIMAX),
     &        WORK(KIMIN),WORK(KICNT),WORK(KNCVAL),WORK(KIDDCP),MXCOEF,
     &        NMORDR,NDIME,NTYPE,NDERIV,NINTIN,LDPMTX,IFRSTD,NLDPMX,
     &        LWRK1,.TRUE.)
C
C        **************************************
C        *** Resetting property derivative. ***
C        **************************************
C
         CNMPRP = CPRPBK
C
C        **********************************************************
C        *** Assigning values to the dependent force constants. ***
C        **********************************************************
C
         IF (NAORDR.EQ.0) THEN
            CALL ADDPFC(DERIV,DCOEFF,KDPMTX,NMIDPC,LDPMTX,IFRSTD,
     &                  NDERIV,NLDPMX,IPRINT)
         END IF
C
C        *****************************************
C        *** Calculating property derivatives. ***
C        *****************************************
C
         IF ((NPRPDR).AND.((.NOT.NRMCRD).OR.
     &                          (NRMCRD.AND..NOT.FRSTNM))) THEN
C
C           *** Workaround to avoid common commonblock ***
C           *** variables.                             ***
C
            CALL STPPVR
C
            KCOEF  = 1
            KIMAX  = KCOEF  + (2*MXCOEF+1)*(NMORDR+1)
            KIMIN  = KIMAX  +               NMRDRP
            KICNT  = KIMIN  +               NMRDRP
            KNCVAL = KICNT  +               NTYPE
            KIDDCP = KNCVAL +               NCOOR
            KLAST  = KIDDCP +               NCOOR
            NPPDER = NDCOOR
            IF (NMRDRP.GE.2) NPPDER = NPPDER + NDCOOR*(NDCOOR+1)/2
            IF (NMRDRP.GE.3) NPPDER = NPPDER
     &                              + NDCOOR*(NDCOOR+1)*(NDCOOR+2)/6
            IF (NMRDRP.GE.4) NPPDER = NPPDER
     &                      + NDCOOR*(NDCOOR+1)*(NDCOOR+2)*(NDCOOR+3)/24
C
C
            IF (DOCCSD) THEN
C
C           *** Derivatives of cc-properties. ***
C
               KCCPRP = KLAST
               KDCCPR = KCCPRP + NMPCAL*NPRPC
               KLAST  = KDCCPR + NPPDER*NPRPC
            ELSE
               IF (SPNSPN) THEN
                  KSPNSP = KLAST
                  KDSPSP = KSPNSP + 6*NMPCAL*NCOOR**2
                  KLAST  = KDSPSP + 6*NPPDER*NCOOR**2
               END IF
               IF (DODIPS) THEN
                  KTRLEN = KLAST
                  KDRTRL = KTRLEN + 3*NMPCAL*NSYM*MXNEXI
                  KEXEFV = KDRTRL + 3*NPPDER*NSYM*MXNEXI
                  KLAST  = KEXEFV +   NMPCAL*NSYM*MXNEXI
               END IF
            END IF
C
            LWRK1  = LWORK - KLAST
            IF (LWRK1.LT.1) CALL QUIT('Memory exceeded in PRPDER')
            CALL PRPDER(SYMCOR,WORK(KDSPSP),WORK(KCOEF),WORK(KSPNSP),
     &           WORK(KTRLEN),WORK(KDRTRL),WORK(KEXEFV),WORK(KCCPRP),
     &           WORK(KDCCPR),GRIREP,WORK(KLAST),IADRSS,KDPMTX,
     &           ICRIRP,INDSTP,IDCOMP,WORK(KIMAX),WORK(KIMIN),
     &           WORK(KICNT),WORK(KNCVAL),WORK(KIDDCP),MXCOEF,NTYPE,
     &           NPPDER,LDPMTX,IFRSTD,NLDPMX,MXNEXI,NSYM,LWRK1,
     &           IPRINT)
         END IF
C
C        *******************
C        *** Test print. ***
C        *******************
C
         IF (PGMTST) THEN
            WRITE (LUPRI,'(A)') '                                      '
            WRITE (LUPRI,'(A)')
     &        'Test "equal energy for partner geometries" is complete.'
            WRITE (LUPRI,'(A,F20.12)') 'Maximum error in energy is:',
     &             EMAX
            WRITE (LUPRI,'(A,F20.12)') 'Relative error: '
            DO IRDR = 1, NMORDR
               WRITE (LUPRI,'(I2,A,F14.8)') IRDR, '. derivative: ',
     &              EMAX/(DISPLC**(DBLE(IRDR)))
            END DO
         END IF
      END IF
C
      CALL QEXIT('NMDER')
      RETURN
      END
C
C  /* Deck gtnpnt */
      SUBROUTINE GTNPNT(FUVAL,GRIREP,SYMCOR,EGRAD,EHESS,COOR,CSTART,
     &                  WORK,WRKDLM,INDSTP,ICRIRP,NPRTNR,ICIN,KDPMTX,
     &                  IRPIND,IDCOMP,LDPMTX,IFRSTD,NLDPMX,IORDR,IRSRDR,
     &                  IINNER,NMPRTN,NDIME,MAXINR,LWORK,NMBKLN,MBKLIN,
     &                  IDIME,NINTIN,IDIMAX,IDIMIN,LURSTR,
c#if defined (VAR_MPI)
c     &                  FTVAL,RSTDON,PRTNR,ALRCAL,CLNRGY,LASTE)
c#else
     &                  RSTDON,PRTNR,ALRCAL,CLNRGY,LASTE)
c#endif
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (D0=0.0D0)
c#if defined (VAR_MPI)
c      INCLUDE 'mpif.h'
c      DIMENSION FTVAL(NINTIN,NDIME)
c#endif
c#if defined (VAR_MPI2)
cC
cC     pario.h will no longer be needed as an include file when locking RMA
cC     operations become available.
cC
c#include <dummy.h>
c#include <pario.h>
c#endif
#include <abainf.h>
#include <cbinum.h>
#include <cbiwlk.h>
#include <optinf.h>
#include <trkoor.h>
#include <infpar.h>
#include <numder.h>
#include <fcsym.h>
#include <taymol.h>
#include <past.h>
#include <pvibav.h>
#include <gnrinf.h>
#include <nuclei.h>
c
#include <huckel.h>
      LOGICAL RUNPNT, CLNRGY, RSTDON, PRTNR, EXSIR, EXHER, EXABA, EXESG,
     &        ALRCAL, PRPCAL, FNDKEY
      CHARACTER*80 MBKLIN
      DIMENSION MBKLIN(NMBKLN)
      DIMENSION FUVAL(NINTIN,NDIME), COOR (NCOOR), SYMCOR(NCOOR,NCOOR),
     &          CSTART(NCOOR), EGRAD(MXCOOR), EHESS(MXCOOR,MXCOOR),
     &          GRIREP(NGORDR,NGVERT), WORK(LWORK)
      DIMENSION INDSTP(NTORDR), ICRIRP(NCOOR,2), NPRTNR(MAXINR),
     &          ICIN(NMORDR), KDPMTX(LDPMTX,NSTRDR,IFRSTD),
     &          IRPIND(NMORDR), IDCOMP(NCOOR)
C
      CALL QENTER('GTNPNT')
C
C     *** Symmetry initilization. ***
C
      PRTNR = .FALSE.
C
      KDIM = 3*NCOOR
      CALL DCOPY(KDIM,CSTART,1,COOR,1)
C
C     *** If property derivatives are calculated at this ***
C     *** geometry, som variables needs to be set.       ***
C
      IF ((.NOT.NRMCRD).OR.(NRMCRD.AND..NOT.FRSTNM)) THEN
         IF ((PRPVIB).AND.(IORDR.LE.1)) THEN
            NMPCAL = NMPCAL + 1
            CNMPRP = .TRUE.
         ELSE IF ((.NOT.PRPVIB).AND.(IORDR.LE.NMRDRP)
     &                         .AND.(NMRDRP.GT.0)) THEN
            NMPCAL = NMPCAL + 1
            CNMPRP = .TRUE.
         ELSE
            CNMPRP = .FALSE.
         END IF
      END IF
C
C     *** Making the appropriate step, if any. ***
C
      IF (IDIME.NE.1) THEN
         DO 1700 IC = 1, IRSRDR+1
         DO 1700 IMXN = 1, IDCOMP(INDSTP(IC))
            CALL STPCOR(COOR,COOR,SYMCOR,DISPLC,NCOOR,ICIN(IC),
     &                  INDSTP(IC),IPRINT)
 1700    CONTINUE
C
C        *** Symmetry of derivatives calculated ***
C
         CLNRGY = .FALSE.
         KIDTMP = 1
         KIDDBT = KIDTMP + NMORDR
         KIRPDG = KIDDBT + NMORDR
         KIRPST = KIRPDG + NMORDR
         KLAST  = KIRPST + NMORDR
         LWRK   = LWORK  - KLAST + 1
         CALL FCSCRN(GRIREP,WORK(KLAST),KDPMTX,INDSTP,ICRIRP,IRPIND,
     &               WORK(KIDTMP),WORK(KIDDBT),WORK(KIRPDG),
     &               WORK(KIRPST),NPRTNR,LWRK,NLDPMX,LDPMTX,IFRSTD,
     &               IORDR,IRSRDR,MAXINR,IINNER,NMPRTN,IPRINT,CLNRGY,
     &               PRTNR,ALRCAL,.FALSE.)
C
C        *** Test print. ***
         IF (IPRINT.GT.10) THEN
            WRITE (LUPRI,'(A, 12I5)') 'Component: ',
     &                                         (INDSTP(I),I=1,IRSRDR+1)
         END IF
      END IF
C
C     *** Information needed for .MANUAL keyword. ***
C
      IF (MANUAL) THEN
         WRITE (LUPRI,'(A)') '                                    '
         WRITE (LUPRI,'(5X,A,I5)') 'Manual geometry calculated', IDIME
         IF (PRTNR) THEN
            WRITE (LUPRI,'(5X,A,I5)') 'Partner geometry',
     &           IDIME - IINNER + NPRTNR(NMPRTN)
         END IF
      END IF
C
C     *** Calculate the energy, gradient or hessian. If ***
C     *** this is a paralell job, we need to find the   ***
C     *** proper processor.                             ***
C
      IF (RUNPNT(CLNRGY,IRSRDR+1,IDIME)) THEN
C
C        *** Another calculation. ***
C
         NUMCAL = NUMCAL + 1
C
c#if defined (VAR_MPI)
c#if defined (VAR_MPI2)
cC
cC     The following code should be replaced with RMA operations following
cC     the MPI-2 standard and as exemplified in for example Fig.6-9 in
cC     "Using MPI-2" by Gropp, Lusk, and Thakur, and WIN_LOCK becomes
cC     available in an MPI-2 implementation
cC
c 199     CONTINUE
c         OPEN(UNIT=99,FILE=WRKDIR(1:LENWRK)//'LOCK',STATUS='NEW',
c     &        FORM='FORMATTED',ERR=199)
c         CALL GPOPEN(LUNMCL,WRKDIR(1:LENWRK)//'NUMCAL','OLD',' ',
c     &        'FORMATTED',IDUMMY,.FALSE.)
c         READ (LUNMCL,'(I5)') NUMCL
c         IF (NUMCL .EQ. NUMCAL) THEN
c            REWIND (LUNMCL)
c            WRITE (LUNMCL,'(I5)') NUMCL + 1
c            CALL GPCLOSE(LUNMCL,'KEEP')
c            CLOSE(UNIT=99,STATUS='DELETE',ERR=107)
c 107        CONTINUE
c#else
c         IF (MYNUM.EQ.MOD(NUMCAL,(NODTOT+1))) THEN
c#endif
c#endif
C
C           *** Header print ***
C
            CALL HEADER('Next numerical derivative component',0)
            IF (IDIME.NE.1) THEN
               WRITE (LUPRI,'(12I5)')  (INDSTP(I),I=1,IRSRDR+1)
            ELSE
               WRITE (LUPRI,'(5X,A)') 'Starting geometry.'
            END IF
C
            IF ((.NOT.DRYRUN).AND.((.NOT.RESTRT).OR.(RESTRT.AND.RSTDON
     &           .AND.((IDIME.GT.IDIMAX).OR.(IDIME.LT.IDIMIN))))) THEN
C
C              *** Reducing symmetry in the DALTON.INP file. ***
C
c               CALL DALCHG(INDSTP,ICRIRP,IRSRDR,IPRINT,NCOOR,NMORDR,
c     &                     .FALSE.)
C
C              *** Update MOLECULE.INP file. ***
C
               KNATOM = 1
               KNLINE = KNATOM + NUCIND
               KNATIN = KNLINE + NUCIND
               CALL UPCOMB(COOR,WORK(KNATOM),WORK(KNLINE),WORK(KNATIN),
     &                     MBKLIN,NMBKLN,LASTE)
C
C              *** Reset necessary variables ***
C
               CALL NDER_RESET(EXHER,EXSIR,EXABA)
cdj
               EXESG = FNDKEY('*ESG   ')
C
C              *** Find the energy, gradient or hessian ***
C
               IF (NAORDR .EQ. 0) THEN
                  CALL GTNRGY(EXHER,EXSIR,EXABA,EXESG,
     $                 WORK,LWORK,WRKDLM)
c#if defined (VAR_MPI)
c                  FTVAL(1,IDIME) = ENERGY
c#else
                  FUVAL(1,IDIME) = ENERGY
c#endif
                  WRITE(LURSTR,'(2I8,F24.16)') 1, IDIME, ENERGY
                  IF (PRTNR) THEN
                     IMDIME = IDIME - IINNER + NPRTNR(NMPRTN)
                     FUVAL(1,IMDIME) = ENERGY
                     WRITE(LURSTR,'(2I8,F24.16)') 1, IMDIME, ENERGY
                  END IF
                  CALL FLSHFO(LURSTR)
C
                  IDIME = IDIME + 1
C
C                 ****************************************************
C                 *** If derivatives of properties are calculated. ***
C                 ****************************************************
C
                  IF (CNMPRP) THEN
C
C                    ************************************************
C                    *** Calculating properties for this geometry ***
C                    ************************************************
C
                     PASEXC = .FALSE.
                     DOWALK = .FALSE.
                     WRINDX = .TRUE.
                     LUSUPM = -1
                     WORK(1) = WRKDLM
                     CALL ABAINP('**PROPE',WORK(2),LWORK)
                     CALL EXEABA(WORK(1),LWORK-1,WRKDLM)
C
                     KTRAMT = 1
                     KCRTPR = KTRAMT + NCOOR**2
                     CALL TRMTOC(WORK(KTRAMT),COOR,WORK(KCRTPR),NCOOR,
     &                           IPRINT)
                  END IF
C
               ELSE IF (NAORDR .EQ. 1) THEN
C                 ******************************************************
C                 *** If gradient is used, property calculations are ***
C                 *** run through GTGRAD.                            ***
C                 ******************************************************
C
                  CALL GTGRAD(EGRAD,EXHER,EXSIR,EXABA,WORK,LWORK,
     &                        WRKDLM)
C
                  KSEGRD = 1
                  KLAST  = KSEGRD + NCOOR
                  LWRK1  = LWORK  - KLAST + 1
                  IF (KLAST.GT.LWORK)
     &                           CALL QUIT('Memory exceeded in TRFCGD')
                  CALL TRFCGD(EGRAD,SYMCOR,COOR,WORK(KSEGRD),
     &                        WORK(KLAST),NCOOR,NDCOOR,LWRK1,IPRINT)
C
                  DO ICOOR = 1, NDCOOR
c#if defined (VAR_MPI)
c                     FTVAL(ICOOR,IDIME) = EGRAD(ICOOR)
c#else
                     FUVAL(ICOOR,IDIME) = EGRAD(ICOOR)
c#endif
                     WRITE(LURSTR,'(2I8,F24.16)') ICOOR, IDIME,
     &                                           EGRAD(ICOOR)
                  END DO
                  IDIME = IDIME + 1
                  CALL FLSHFO(LURSTR)
               ELSE IF (NAORDR .EQ. 2) THEN
C                 ******************************************************
C                 *** If hessian is used, property calculations are  ***
C                 *** run through GTHESS.                            ***
C                 ******************************************************
C
                  KAHESS = 1
                  KLAST  = KAHESS + NCOOR**2
                  LWRK   = LWORK - KLAST +1
                  CALL GTHESS(EGRAD,EHESS,WORK(KAHESS),EXHER,EXSIR,
     &                        EXABA,WORK(KLAST),LWRK,WRKDLM)
C
                  KSEHSS = 1
                  KLAST  = KSEHSS + NCOOR**2
                  LWRK1  = LWORK  - KLAST + 1
                  IF (KLAST.GT.LWORK)
     &                           CALL QUIT('Memory exceeded in TRFCHS')
                  CALL TRFCHS(EHESS,SYMCOR,COOR,WORK(KSEHSS),
     &                        WORK(KLAST),NCOOR,NDCOOR,LWRK1,IPRINT)
C
                  ICOOR12 = 0
                  DO ICOOR2 = 1, NDCOOR
                  DO ICOOR1 = 1, ICOOR2
                     ICOOR12 = ICOOR12 + 1
c#if defined (VAR_MPI)
c                     FTVAL(ICOOR12,IDIME) = EHESS(ICOOR1,ICOOR2)
c#else
                     FUVAL(ICOOR12,IDIME) = EHESS(ICOOR1,ICOOR2)
c#endif
                     WRITE(LURSTR,'(2I8,F24.16)') ICOOR12, IDIME,
     &                                           EHESS(ICOOR1,ICOOR2)
                  END DO
                  END DO
                  IDIME = IDIME + 1
                  CALL FLSHFO(LURSTR)
               END IF
            ELSE
               IDIME = IDIME + 1
            END IF
c#if defined (VAR_MPI)
c         ELSE
c#if defined (VAR_MPI2)
c            CALL GPCLOSE(LUNMCL,'KEEP')
c            CLOSE(UNIT=99,STATUS='DELETE',ERR=108)
c 108        CONTINUE
c#endif
c            DO INTIN = 1, NINTIN
c               FUVAL(INTIN,IDIME) = D0
c            END DO
c            IDIME = IDIME + 1
c         END IF
c#endif
      ELSE IF (.NOT.ALRCAL) THEN
         IF ((.NOT.RESTRT).OR.(RESTRT.AND.
     &              RSTDON.AND.(IDIME.GT.IDIMAX))) THEN
            DO 2100 INTIN = 1, NINTIN
               FUVAL(INTIN,IDIME) = D0
               WRITE(LURSTR,'(2I8,F24.16)') INTIN, IDIME, D0
 2100       CONTINUE
            CALL FLSHFO(LURSTR)
         END IF
         IDIME = IDIME + 1
      END IF
C
      CALL QEXIT('GTNPNT')
      RETURN
      END
C
C  /* Deck upcomb */
      SUBROUTINE UPCOMB(COOR,NATOM,NLINE,NATIND,MBKLIN,NMBKLN,LASTE)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
C
      PARAMETER (D100 = 100.0D0, THRSH = 1.0D-12)
#include <cbirea.h>
#include <molinp.h>
#include <nuclei.h>
#include <trkoor.h>
#include <symmet.h>
#include <numder.h>
#include <inftap.h>
      LOGICAL BIG, USED, LASTE, DOCART, DOOWN, AUTOSY, NOSYM, ADDSYM,
     &        NEWINP, NEWATO
      CHARACTER*6 CHR
      CHARACTER*4 NAME
      CHARACTER*80 MBKLIN, BSNM
      CHARACTER*11 TMPTXT
      CHARACTER*1 KASYM(3,3), ID3, CRT
      DIMENSION MBKLIN(NMBKLN), JCO1(MXAQN)
      DIMENSION COOR(3,NCOOR/3), NATOM(NUCIND), NLINE(NUCIND),
     &          NATIND(NUCIND)
#include <ibtfun.h>
C
C     Updates geometry in common block
C
      IF (LASTE) THEN
         NMLINE = NMBKLN
         DO 100 IBKLIN = 1, NMBKLN
            MLINE(IBKLIN) = MBKLIN(IBKLIN)
 100     CONTINUE
      ELSE
         NADD   = 0
         IATOM  = 0
         NCLAST = 0
         CALL UPCASE(MLINE(NMLAU))
         IPOS = INDEX(MLINE(NMLAU),'ATO')
         IF (IPOS .EQ. 0) THEN
            IF (MLINE(NMLAU)(10:10).EQ.'0') THEN
               MLINE(NMLAU)(20:20) = ' '
            ELSE
               MLINE(NMLAU)(10:20) = '           '
            END IF
         ELSE
            CALL LINE4(MLINE(NMLAU),NONTYP,NSYMOP,CRT,KCHARG,THRS,
     &                 ADDSYM,KASYM,ID3,DOCART,DOOWN)
            AUTOSY = .TRUE.
            NOSYM = .FALSE.
            ID3 = ' '
            CALL LINE4W(MLINE(NMLAU),NONTYP,NSYMOP,KCHARG,THRS,AUTOSY,
     &                  NOSYM,KASYM,ID3,DOCART,DOOWN)
         END IF
         DO 200 ICENT = 1, NUCIND
            ISYM   = 0
            NRLINE = NCLINE(ICENT)
            NC     = NCLINE(ICENT)
            MULCNT = ISTBNU(ICENT)
            IF (NC .NE. 0) THEN
              READ (MLINE(NC),9100) NAME
              IPOS = INDEX(MLINE(NC),'Isotope=')
              IF (IPOS .NE. 0) THEN
                 READ (MLINE(NC)(IPOS:),'(A11)') TMPTXT
              ELSE
                 TMPTXT = '           '
              END IF
              DO 300 IOP = 0, MAXOPR
                  IF (IBTAND(IOP,MULCNT) .EQ. 0) THEN
                     IATOM = IATOM + 1
                     CRX = COOR(1,IATOM)
                     CRY = COOR(2,IATOM)
                     CRZ = COOR(3,IATOM)
                     BIG = (ABS(CRX) .GE. D100 .OR.
     *                      ABS(CRY) .GE. D100 .OR.
     *                      ABS(CRZ) .GE. D100)
                     IF (ISYM .GT. 0) THEN
                        DO 400 I = NMLINE, NC+1, -1
                           MLINE(I+1) = MLINE(I)
 400                    CONTINUE
                        DO 450 IC2 = ICENT+1, NUCIND
                           NCLINE(IC2) = NCLINE(IC2) + 1
 450                    CONTINUE
                        NRLINE = NRLINE + 1
                        NMLINE = NMLINE + 1
                        NC     = NC     + 1
                     END IF
                     IF (BIG) THEN
                        WRITE (MLINE(NC),9200) NAME,CRX,CRY,CRZ,TMPTXT
                     ELSE
                        WRITE (MLINE(NC),9300) NAME,CRX,CRY,CRZ,TMPTXT
                     END IF
                     ISYM = ISYM + 1
                  END IF
 300           CONTINUE
            END IF
 200     CONTINUE
C
C        Do a count of each type of atom....
C
         CALL UPCASE(MLINE(NCLINE(1)-1))
         NEWINP = (INDEX(MLINE(NCLINE(1) - 1),'CHA') .NE. 0)
         KCENT2 = 0
         ICENT1 = 1
 500     CONTINUE
         ICENT1 = ICENT1 + KCENT2
         DO 700 ICENT2 = ICENT1+1, NUCIND
            IF (NEWINP) THEN
               CALL UPCASE(MLINE(NCLINE(ICENT2)-1))
               NEWATO = (INDEX(MLINE(NCLINE(ICENT2) - 1),'CHA') .NE. 0)
            ELSE
               READ (MLINE(NCLINE(ICENT2)-1),'(A)') CHR
               NEWATO = (CHR .EQ. '      ')
            END IF
            IF (NEWATO) THEN
               NAT = 0
               DO 800 IC = ICENT1, ICENT2-1
                  NAT = NAT + NUCDEG(IC)
 800           CONTINUE
               IF (NEWINP) THEN
                  IF (MLINE(1)(1:5) .EQ. 'BASIS') BASIS  = .TRUE.
                  IF (MLINE(1)(1:5) .EQ. 'ATOMB') ATOMBA = .TRUE.
                  CALL UPCASE(MLINE(NCLINE(ICENT1)-1))
                  CALL LINE5R(MLINE(NCLINE(ICENT1) - 1),Q1,NONT1,MBSI1,
     &                 IQM1,JCO1,MXAQN,BASIS,ATOMBA,LMULBS,BSNM,
     &                 RADIUS_PCM, ALPHA_PCM)
                  CALL LINE5W(MLINE(NCLINE(ICENT1) - 1),Q1,NAT,MBSI1,
     &                 BASIS,ATOMBA,LMULBS,BSNM,IQM1,JCO1,MXAQN,
     &                 RADIUS_PCM, ALPHA_PCM)
               ELSE
                  WRITE (MLINE(NCLINE(ICENT1)-1)(13:15),'(I3)') NAT
               END IF
               KCENT2 = ICENT2-ICENT1
               GOTO 500
            ELSE
               IF (ICENT2 .EQ. NUCIND) THEN
                  KCENT1 = ICENT1
                  GOTO 900
               END IF
            END IF
 700     CONTINUE
         KCENT1 = NUCIND
 900     CONTINUE
         NAT = 0
         DO 1100 IC = KCENT1, NUCIND
            NAT = NAT + NUCDEG(IC)
 1100    CONTINUE
         IF (NEWINP) THEN
            CALL UPCASE(MLINE(NCLINE(KCENT1)-1))
            CALL LINE5R(MLINE(NCLINE(KCENT1) - 1),Q1,NONT1,MBSI1,
     &           IQM1,JCO1,MXAQN,BASIS,ATOMBA,LMULBS,BSNM,
     &           RADIUS_PCM, ALPHA_PCM)
            CALL LINE5W(MLINE(NCLINE(KCENT1) - 1),Q1,NAT,MBSI1,
     &           BASIS,ATOMBA,LMULBS,BSNM,IQM1,JCO1,MXAQN,
     &           RADIUS_PCM, ALPHA_PCM)
            BASIS  = .FALSE.
            ATOMBA = .FALSE.
         ELSE
            WRITE (MLINE(NCLINE(KCENT1)-1)(13:15),'(I3)') NAT
         END IF
      END IF
C
C     Punch MOLECULE input with updated coordinates to LUMOL
C     (And DALTON.OUT if manual is set.)
C
      IF (MANUAL) THEN
         WRITE (LUPRI,'(5X,A)') 'Molecular geometry as requested:'
         WRITE (LUPRI,'(5X,A,I5)') 'Number of lines printed', NMLINE
      END IF
      CALL GPOPEN(LUMOL,'MOLECULE.INP','OLD',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND (LUMOL)
      DO 1300 IMLINE = 1,NMLINE
         WRITE (LUMOL,'(A)') MLINE(IMLINE)
         IF (MANUAL) THEN
            WRITE (LUPRI,'(A)') MLINE(IMLINE)
         END IF
 1300 CONTINUE
      CALL GPCLOSE(LUMOL,'KEEP')
C
 9100 FORMAT (A4)
 9200 FORMAT (A4,3F20.10,1X,A11)
 9300 FORMAT (A4,3F20.15,1X,A11)
C
      RETURN
      END
C
C  /*Deck stpcor*/
      SUBROUTINE STPCOR(COOR,CSTART,SYMCOR,DISPLC,NCOOR,KPM,KSCOOR,
     &                  IPRINT)
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION COOR(NCOOR), CSTART(NCOOR), SYMCOR(NCOOR,NCOOR)
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('CSTART in STPCOR',1)
         CALL OUTPUT(CSTART,1,1,1,NCOOR,1,NCOOR,1,LUPRI)
      END IF
C
      IF (KPM.EQ.1) THEN
        FAC =  DISPLC
      ELSE
        FAC = -DISPLC
      END IF

      DO 100 ICOOR = 1, NCOOR
         COOR(ICOOR) = CSTART(ICOOR) + FAC*SYMCOR(ICOOR,KSCOOR)
 100  CONTINUE
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('COOR in STPCOR',1)
         CALL OUTPUT(COOR,1,1,1,NCOOR,1,NCOOR,1,LUPRI)
      END IF
C
      RETURN
      END
C
C  /* Deck numres */
      SUBROUTINE NDER_RESET(EXHER,EXSIR,EXABA)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
C
#include <ccorb.h>
#include <optinf.h>
#include <symmet.h>
#include <nuclei.h>
#include <gnrinf.h>
#include <huckel.h>
#include <trkoor.h>
#include <cbiwlk.h>
#include <past.h>
#include <abainf.h>
#include <cbinum.h>
#include <numder.h>
      LOGICAL EXHER,EXSIR,EXABA,EX
C
C     This routine resets a few variables, to be able to calculate
C     energy again for a new geometry (and symmetry).
C
      EXHER  = .FALSE.
      EXSIR  = .FALSE.
      EXABA  = .FALSE.
      RDINPC = .FALSE.
      RDMLIN = .FALSE.
C
C     *** unset ABA variable to false ***
C
      CALL ABUSET
C
C     *** If there are possibilities for  new symmetry. ***
C
      IF (((NMORDR+NAORDR).GT.1).OR.(NPRPDR)) THEN
         NEWSYM = .TRUE.
         DOHUCKEL = .TRUE.
      END IF
C
C     *** For higher order derivatives. ****
C
      IF ((NMORDR+NAORDR).GT.1) THEN
         HRINPC = .FALSE.
         KEEPHE = .FALSE.
         RSTARR = .TRUE.
         DOWALK = .FALSE.
         BRKSYM = .FALSE.
         ITRBRK = ITRNMR
         INDOLD = INDTOT
         GECONV = .FALSE.
         CALL IZERO(NUCNUM, MXCENT*8)
         CALL IZERO(NCRREP, 16)
         CALL IZERO(IPTCNT, MXCENT*48)
         CALL IZERO(NAXREP, 16)
         CALL IZERO(INDHES, 8)
C
C        *** For analytical hessians **
C
         ITRNMR = 1
         NCRTOT = NCOOR
         NCART  = NCOOR
         DO I = 0, 7
            DOREPW(I) = .TRUE.
        END DO
      END IF
C
      IF (NPRPDR) THEN
         IF (SPNSPN) THEN
            PASTRP = .FALSE.
         END IF
      END IF
C
C     *** Initialization related to doing CC. ***
C
Ctbp  IF (DOCCSD) THEN
Ctbp     DO ISYM  = 1, 8
Ctbp     DO IXFRO = 1, MAXFRO
Ctbp        FRORHF(IXFRO,ISYM) = .FALSE.
Ctbp     END DO
Ctbp     END DO
Ctbp  END IF
C
C     *** For spin-spin couplings. ***
C
      CALL GPINQ('RSPVEC','EXIST',EX)
      IF (EX) THEN
         LURSP = -1
         CALL GPOPEN(LURSP,'RSPVEC','OLD',
     &               ' ','UNFORMATTED',IDUMMY,.FALSE.)
         CALL GPCLOSE(LURSP,'DELETE')
      END IF
cC
c      CALL GPINQ('RSPRST.E2C','EXIST',EX)
c      IF (EX) THEN
c         CALL GPOPEN(LURSP,'RSPRST.E2C','OLD',' ','UNFORMATTED',IDUMMY,
c     &               .FALSE.)
c         CALL GPCLOSE(LURSP,'DELETE')
c      END IF
cC
c      CALL GPINQ('AOPROPER','EXIST',EX)
c      IF (EX) THEN
c         CALL GPOPEN(LURSP,'AOPROPER','OLD',' ','UNFORMATTED',IDUMMY,
c     &               .FALSE.)
c         CALL GPCLOSE(LURSP,'DELETE')
c      END IF
cC
c      CALL GPINQ('ABACUS.RESTART','EXIST',EX)
c      IF (EX) THEN
c         CALL GPOPEN(LURSP,'ABACUS.RESTART','OLD',' ','UNFORMATTED',
c     &               IDUMMY,.FALSE.)
c         CALL GPCLOSE(LURSP,'DELETE')
c      END IF
cC
c      CALL GPINQ('ABAENR.RST','EXIST',EX)
c      IF (EX) THEN
c         CALL GPOPEN(LURSP,'ABAENR.RST','OLD',' ','UNFORMATTED',
c     &               IDUMMY,.FALSE.)
c         CALL GPCLOSE(LURSP,'DELETE')
c      END IF
cC
c      CALL GPINQ('ABACUS.GDT','EXIST',EX)
c      IF (EX) THEN
c         CALL GPOPEN(LURSP,'ABACUS.GDT','OLD',' ','UNFORMATTED',
c     &               IDUMMY,.FALSE.)
c         CALL GPCLOSE(LURSP,'DELETE')
c      END IF
cC
c      CALL GPINQ('ABACUS.RDT','EXIST',EX)
c      IF (EX) THEN
c         CALL GPOPEN(LURSP,'ABACUS.RDT','OLD',' ','UNFORMATTED',
c     &               IDUMMY,.FALSE.)
c         CALL GPCLOSE(LURSP,'DELETE')
c      END IF
cC
c      CALL GPINQ('DRCINT','EXIST',EX)
c      IF (EX) THEN
c         CALL GPOPEN(LURSP,'DRCINT','OLD',' ','UNFORMATTED',
c     &               IDUMMY,.FALSE.)
c         CALL GPCLOSE(LURSP,'DELETE')
c      END IF
C
      RETURN
      END
C
C
C     /* Deck nmpres */
c      SUBROUTINE NMPRES
c#include <implicit.h>
c#include <priunit.h>
c
C
C
      SUBROUTINE PRNHST(KSCOO1,KSCOO2,DISPLC,I,J,NC)
C
C     Prints out headers to the different steps in finite field.
C
#include <implicit.h>
#include <priunit.h>
      CHARACTER SIGNS*4
C
      SIGNS = '    '
C
      CALL TITLER('Finite field geometry','*',118)
C
      IF (NC .EQ. 1) THEN
         WRITE (LUPRI,'(A,I4)') '        Distortions in coordinate:',
     &                           KSCOO2
      ELSE
         WRITE (LUPRI,'(A,2I4)') '        Distortions in coordinates:',
     &                           KSCOO1, KSCOO2
      END IF
C
      WRITE (LUPRI,'(/A,F12.8/)') '        Length of displacement:',
     &   DISPLC
C
      IF (J .EQ. 1) THEN
         SIGNS(2:2) = '+'
      ELSE
         SIGNS(2:2) = '-'
      END IF
      IF (NC .NE. 1) THEN
         IF (I .EQ. 1) THEN
            SIGNS(4:4) = '+'
         ELSE
            SIGNS(4:4) = '-'
         END IF
      END IF
      WRITE (LUPRI,'(A24,A4//)') '        Step directions:', SIGNS
C
      RETURN
      END
C /* Deck addnuc */
c      SUBROUTINE ADDNUC(WORK,LWORK)
c#include <implicit.h>
c#include <priunit.h>
c#include <mxcent.h>
cC
c      DIMENSION WORK(LWORK)
cC
c#include <energy.h>
c#include <taymol.h>
c#include <trkoor.h>
c#include <cbiwlk.h>
cC
cC     Calculate nuclear contributions
cC
c      NDIM = 3*MXCOOR*MXCOOR
c      IF (LWORK .LT. NDIM) CALL STOPIT('ADDNUC',' ',LWORK,NDIM)
c      CALL NUCREP(WORK,WORK(MXCOOR*MXCOOR+1),WORK(2*MXCOOR*MXCOOR+1))
C
C     Add nuclear contributions to gradient and Hessian
C
c      IF (NUMGRD) THEN
c         DO 100 I = 1, NCOOR
c            GRDMOL(I) = GRDMOL(I) + GRADNN(I)
c  100    CONTINUE
c      END IF
C     IF (NUMHES) THEN
C        DO 200 I = 1, NCOOR
C        DO 200 J = 1, NCOOR
C           HESMOL(J,I) = HESMOL(J,I) + WORK((I-1)*MXCOOR + J)
C 200    CONTINUE
C     END IF
c      RETURN
c      END
C /* Deck GTCOMP */
      LOGICAL FUNCTION GTCOMP()
#include <implicit.h>
#include <cbitro.h>
      GTCOMP = COMPAR
      RETURN
      END
C
C
C /* Deck NMCOEF */
      SUBROUTINE NMCOEF(COEFF,TCOEFF,WORK,MXCOEF,NMNMDR,LWORK)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D0, D05 = 0.5D0)
C
      DIMENSION COEFF (-MXCOEF:MXCOEF,0:NMNMDR),
     &          TCOEFF(-NMNMDR:NMNMDR,0:NMNMDR), WORK(LWORK)

C
      KDIM1 = (2*MXCOEF+1)*(NMNMDR+1)
      KDIM2 = (2*NMNMDR+1)*(NMNMDR+1)
      CALL DZERO(COEFF ,KDIM1)
      CALL DZERO(TCOEFF,KDIM2)
      COEFF(0,0) = D1
C
C     *** Temporary coeffesients used to generate coefficients ***
C     ***            for even-numbered derivatives.            ***
C
      IF (NMNMDR .GT. 1) THEN
         NCOR = 1
         TCOEFF(0,0) = D1
         DO 100 IDR  = 1, NMNMDR
            DO 200 ICOR = -IDR+1, IDR-1
               TCOEFF(ICOR,IDR) = TCOEFF(ICOR-1,IDR-1)
     &                          - TCOEFF(ICOR+1,IDR-1)
 200        CONTINUE
C
            TCOEFF(-IDR,IDR) = (-D1)**IDR
            TCOEFF( IDR,IDR) =   D1
C
 100     CONTINUE
C
C        *** Coefficients for even-numbered derivatives ***
C
         COEFF(0,0) = D1
         NEVEN = INT(NMNMDR/2)
         DO 300 IEVEN = 1, NEVEN
            IDR = 2*IEVEN
            DO 400 IECOR = -IEVEN,IEVEN
               ICOR = 2*IECOR
               COEFF(IECOR,IDR) = TCOEFF(ICOR,IDR)
 400        CONTINUE
 300     CONTINUE
      END IF
C
C     *** Coefficients for odd-numbered derivatives ***
C
      NODD = INT((NMNMDR+1)/2)
      DO 500 IODD = 1, NODD
         IDR = 2*IODD - 1
        DO 600 IOCOR = -IODD+1, IODD-1
            IF (IOCOR .NE. 0) THEN
               COEFF(IOCOR,IDR) = D05*(COEFF(IOCOR-1,IDR-1)
     &                               - COEFF(IOCOR+1,IDR-1))
            END IF
 600     CONTINUE
C
         COEFF(-IODD,IDR) = -D05
         COEFF( IODD,IDR) =  D05
C
 500  CONTINUE
C
      RETURN
      END
C
C
C /* Deck NMNDER */
      SUBROUTINE NMNDER(DERIV,COEFF,FUNVAL,GRIREP,WORK,IADRSS,KDPMTX,
     &                  ICRIRP,INDSTP,INDTMP,IDCOMP,IMAX,IMIN,ICNT,
     &                  NCVAL,IDDCMP,MXCOEF,NORDR,NDIME,NTYPE,NDERIV,
     &                  NFINNR,LDPMTX,IFRSTD,NLDPMX,LWORK,FCCAL)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (DMIN = 1.0D-12, D1=1.0D0, D0=0.0D0)
#include <taymol.h>
#include <trkoor.h>
#include <cbiwlk.h>
#include <cbinum.h>
#include <numder.h>
#include <fcsym.h>
#include <dummy.h>
      LOGICAL CLFVAL, FCCAL, DIAGON
      DIMENSION COEFF(-MXCOEF:MXCOEF,0:NORDR), DERIV(NDERIV),
     &          FUNVAL(NFINNR,NDIME), GRIREP(NGORDR,NGVERT), WORK(LWORK)
      DIMENSION ICNT(NTYPE), IADRSS(NTYPE), IMAX(NORDR), IMIN(NORDR),
     &          INDSTP(NTORDR), INDTMP(NTORDR), IDCOMP(NCOOR),
     &          IDDCMP(NCOOR), NCVAL(NCOOR),
     &          KDPMTX(LDPMTX,NSTRDR,IFRSTD), ICRIRP(NCOOR,2)
C
      CALL QENTER('NMNDER')
C
      IF (FCCAL) THEN
         KDIM = MXCOOR**2
         IF (NAORDR .LT. 1) THEN
            CALL DZERO(GRDMOL,MXCOOR)
         END IF
         IF (NAORDR .LT.2 ) THEN
            CALL DZERO(HESMOL,  KDIM)
         END IF
      END IF
      IF (NDERIV .GT. 0) CALL DZERO(DERIV,NDERIV)
C
C     ***************************************************
C     *** This subroutine calculates the coefficients ***
C     ***      for the numerical differentiation.     ***
C     ***************************************************
C
      KTCOEF = 1
      KLAST  = KTCOEF + (2*NORDR+1)*(NORDR+1)
      LWRK = LWORK - KLAST + 1
      CALL NMCOEF(COEFF,WORK(KTCOEF),WORK(KLAST),MXCOEF,NORDR,LWORK)
C
      IDERIV = 0
      DO 100 IORDR = 1, NORDR
C
         DIAGON = (.NOT.FCCAL ).AND.(IORDR.EQ.NORDR).AND.
     &                  PRPVIB .AND.((NARDRP+NMRDRP).EQ.2)

C        ***********************************************
C        *** Special code for gradients and hessians ***
C        *** due to special memory places.           ***
C        ***********************************************
C
         IF (((IORDR+NAORDR).LE.2).AND.FCCAL) THEN

            IF ((IORDR+NAORDR).EQ.1) THEN
C
C              ****************************************
C              *** Numerical gradient from energies ***
C              ****************************************
C
               IMAX(1)  =  1
               IMIN(1)  = -1
               IDIME    =  1
               HINV     =  D1/DISPLC
               DO 200 ICOOR = 1, NCOOR
               DO 200 I   = IMAX(1), IMIN(1), -2
                  IF (COEFF(I,1)**2 .GT. DMIN) THEN
                     IDIME = IDIME + 1
                     GRDMOL(ICOOR) = GRDMOL(ICOOR)
     &                             + COEFF(I,1)*FUNVAL(1,IDIME)*HINV
                  END IF
 200           CONTINUE
            ELSE
C
C              ******************************
C              *** Numerical Hessian from ***
C              ***  analytical gradients  ***
C              ******************************
C
               IF (NAORDR .EQ. 1) THEN
                  IDIME    = 2
                  HINV     = D1/DISPLC
                  DO 260 ICOOR2 = 1, NCOOR
                     DO 270 ICOOR1 = 1, NCOOR
                        HESMOL(ICOOR1,ICOOR2) =
     &                         (COEFF( 1,1)*FUNVAL(ICOOR1,IDIME      )
     &                        + COEFF(-1,1)*FUNVAL(ICOOR1,IDIME+1))*HINV
 270                 CONTINUE
                     IDIME = IDIME + 2
 260              CONTINUE
               ELSE
C
C              ******************************
C              *** Numerical Hessian from ***
C              ***        energies.       ***
C              ******************************
C
                  CALL IZERO(ICNT,NTYPE)
C
                  IMAX(1)  =  1
                  IMIN(1)  = -1
                  HINV    =  D1/(DISPLC**2)
                  DO 300 IX2 = 1, NDCOOR
                  DO 300 IX1 = 1, IX2
C
                     CALL IZERO(INDSTP,NMORDR)
                     INDSTP(1) = IX2
                     INDSTP(2) = IX1
C
C                    *** Checking whether this component should ***
C                    *** not be calculated, due to symmetry.    ***
C
                     KIDTMP = 1
                     KIDDBT = KIDTMP + NORDR
                     KIRPDG = KIDDBT + NORDR
                     KIRPST = KIRPDG + NORDR
                     KLAST  = KIRPST + NORDR
                     LWRK   = LWORK - KLAST + 1
                     CALL FCSCRN(GRIREP,WORK(KLAST),KDPMTX,INDSTP,
     &                    ICRIRP,IDUMMY,WORK(KIDTMP),WORK(KIDDBT),
     &                    WORK(KIRPDG),WORK(KIRPST),IDUMMY,LWORK,NLDPMX,
     &                    LDPMTX,IFRSTD,IORDR,IORDR-1,IDUMMY,IDUMMY,
     &                    IDUMMY,IPRINT,CLFVAL,.FALSE.,.FALSE.,.TRUE.)
C
                     MX = IMIN(1)
                     IF (IX1 .EQ. IX2) MX = IMAX(1)
                     DO 400 I2 = IMAX(1), IMIN(1),-1
                     DO 400 I1 = IMAX(1), MX, -1
                        IF (IX1 .EQ. IX2) THEN
                           IF (I2.EQ.0) THEN
                              ITYPE       = 1
                              ICNT(ITYPE) = 0
                              BCOEFF      = COEFF(0,2)
                           ELSE
                              ITYPE       = 2
                              BCOEFF      = COEFF(I2,2)
                           END IF
                        ELSE
                           ITYPE       = 3
                           BCOEFF      = COEFF(I1,1)*COEFF(I2,1)
                        END IF
C
                        IF (BCOEFF**2 .GT. DMIN) THEN
                           ICNT(ITYPE) = ICNT(ITYPE) + 1
                           IF (CLFVAL) THEN
                              HESMOL(IX2,IX1) = HESMOL(IX2,IX1)
     &                             + BCOEFF*HINV
     &                             * FUNVAL(1,IADRSS(ITYPE)+ICNT(ITYPE))
                           END IF
                        END IF
 400                 CONTINUE
 300              CONTINUE
               END IF
            END IF
         ELSE IF ((ANALZ1).AND.((NMORDR+NAORDR).EQ.3).AND.NRMCRD.AND.
     &            (FCCAL ).AND. (IORDR.EQ.2).AND.(NAORDR.EQ.1)) THEN
C
C           ****************************************************************
C           *** Special case for cases were only vib. average of         ***
C           *** properties in ANALZ1 method from gradients in normal     ***
C           *** coordinates.                                             ***
C           *** ANALZ1            -> Use ANALZ1 method.                  ***
C           *** NMORDR+NAORDR = 3 -> Calculate up to cubic force field.  ***
C           ***                      When ANALZ1 is also used, then only ***
C           ***                      parts of the cubic force field is   ***
C           ***                      calculated.                         ***
C           *** NRMCRD            -> Normal coordinates.                 ***
C           *** NAORDR = 1        -> Analytical gradients                ***
C           ****************************************************************
C
            ISORDR = IORDR+1
            NINNR2 = 0
C
            NSTP = 1
            DO I = 1, ISORDR
               NSTP = NSTP*(NDCOOR+I-1)/I
            END DO
            POWER  = DBLE(IORDR)
            DIVDIS = D1/(DISPLC**POWER)
C
            CALL IZERO(INDSTP,NMORDR+NAORDR)
C
C           *** NSTP -> Number of components in the ***
C           ***      numerical differentiation      ***
C
            DO ISTP  = 1, NSTP
C
C              *** Finding which component this is ***
C
               DO IC = ISORDR, 1, -1
                  IF (IC .EQ. 1) THEN
                     INDSTP(1) = INDSTP(1) + 1
                     DO I = 2, ISORDR
                        INDSTP(I) = 1
                     END DO
                     GOTO 500
                  ELSE IF (INDSTP(IC) .LE. INDSTP(IC-1)-1) THEN
                     DO I = IC+1, ISORDR
                        INDSTP(I) = 1
                     END DO
                     INDSTP(IC) = INDSTP(IC) + 1
                     GOTO 500
                  END IF
               END DO
 500           CONTINUE
C
C              *** First indices in INDSTP are kept for the ***
C              *** analytical derivative.                   ***
C
               KITMP = 1

               CALL SRTINS(INDSTP,INDTMP,WORK(KITMP),NTORDR,IPRINT)

C              *** Checking whether this component should ***
C              *** not be calculated, due to symmetry.    ***
C
               CLFVAL = .TRUE.
c               IF (FCCAL) THEN
c                  KIDTMP = 1
c                  KIDDBT = KIDTMP + NORDR
c                  KIRPDG = KIDDBT + NORDR
c                  KIRPST = KIRPDG + NORDR
c                  KLAST  = KIRPST + NORDR
c                  LWRK   = LWORK - KLAST + 1
c                  CALL FCSCRN(GRIREP,WORK(KLAST),KDPMTX,INDTMP,ICRIRP,
c     &                        IDUMMY,WORK(KIDTMP),WORK(KIDDBT),
c     &                        WORK(KIRPDG),WORK(KIRPST),IDUMMY,LWORK,
c     &                        NLDPMX,LDPMTX,IFRSTD,IORDR,IORDR-1,IDUMMY,
c     &                        IDUMMY,IDUMMY,IPRINT,CLFVAL,.FALSE.,
c     &                        .FALSE.,.TRUE.)
c               END IF
C
               IDERIV = IDERIV + 1
C
               IF (CLFVAL) THEN
C
                  CALL IZERO(NCVAL,NDCOOR)
                  DO IC = 1, IORDR
                     NCVAL(INDTMP(IC)) = NCVAL(INDTMP(IC)) + 1
                  END DO
C
C                 *******************************************************
C                 *** IDCOMP -> Maks steporder to get the derivative  ***
C                 *** IDDCMP -> Counting array, maks to min steporder ***
C                 *** NTTYPE -> Number of function values needed for  ***
C                 ***           1 component                           ***
C                 *******************************************************
C
                  NTTYPE = 1
                  CALL IZERO(IDCOMP,NDCOOR)
                  CALL IZERO(IDDCMP,NDCOOR)
                  DO IC = 1, NDCOOR
                     IF (NCVAL(IC) .NE. 0) THEN
                        IDCOMP(IC) =   INT((NCVAL(IC)+1)/2)
                        NTTYPE     =   NTTYPE*(2*IDCOMP(IC) + 1)
                     END IF
                  END DO
                  DO I=1,NDCOOR
                     IDDCMP(I) = IDCOMP(I)
                  END DO
C
                  DO ITTYPE = 1, NTTYPE
C
C                    *** Finding the right indices to identify ***
C                    ***        the right function value       ***
C
                     DO IC = 1, NDCOOR
                        IF ((IDDCMP(IC) .GT. -IDCOMP(IC))
     &                              .AND. (ITTYPE .NE. 1)) THEN
                           IDDCMP(IC) = IDDCMP(IC) - 1
                           DO ICT = 1, IC-1
                              IDDCMP(ICT) = IDCOMP(ICT)
                           END DO
                           GOTO 600
                        END IF
                     END DO
 600                 CONTINUE
C
C                    *** Calculate the coefficient for this ***
C                    *** function value                     ***
C
                     BCOEFF = D1
                     NUMCOF = 0
                     DO IC = 1, NDCOOR
                        IF (NCVAL(IC) .NE. 0) THEN
                           BCOEFF = BCOEFF*COEFF(IDDCMP(IC),NCVAL(IC))
                           NUMCOF = NUMCOF + 1
                        END IF
                     END DO
                     IF (NUMCOF .EQ. 0) BCOEFF = D0
C
C                    *** Does the function value contribute? ***
C
                     IF (BCOEFF**2 .GT. DMIN) THEN
C
C                       ************************************************
C                       *** This subroutine finds the adress for the ***
C                       ***      function value, from the indices    ***
C                       ***   NEIND - The adress in the FUVAL-array  ***
C                       ************************************************
C
                        KITCMP = 1
                        CALL GTEIND(IADRSS,IDDCMP,NCVAL,WORK(KITCMP),
     &                              NEIND,ITTYPE,NORDR,IORDR)
C
C                       *** The derivative is calculated. ***
C
                        DERIV(IDERIV) = DERIV(IDERIV)
     &                                + BCOEFF*FUNVAL(INDTMP(3),NEIND)
C
                     END IF
                  END DO
                  DERIV(IDERIV) = DERIV(IDERIV)*DIVDIS
               END IF
            END DO
         ELSE
C
C          *********************************************
C          *** Numerical N'th derivative from NAORDR ***
C          ***         analytival derivative         ***
C          *********************************************
C
            NINNR2 = 0
C
            NSTP = 1
            DO I = 1, IORDR
               NSTP = NSTP*(NDCOOR+I-1)/I
            END DO
            POWER  = DBLE(IORDR)
            DIVDIS = D1/(DISPLC**POWER)
C
            CALL IZERO(INDSTP,NMORDR)
C
C           *** NSTP -> Number of components in the ***
C           ***      numerical differentiation      ***
C
            DO ISTP  = 1, NSTP
C
C              *** Finding which component this is ***
C
               DO IC = IORDR, 1, -1
                  IF (IC .EQ. 1) THEN
                     INDSTP(1) = INDSTP(1) + 1
                     DO I = 2, IORDR
                        INDSTP(I) = 1
                     END DO
                     GOTO 700
                  ELSE IF (INDSTP(IC) .LE. INDSTP(IC-1)-1) THEN
                     DO I = IC+1, IORDR
                        INDSTP(I) = 1
                     END DO
                     INDSTP(IC) = INDSTP(IC) + 1
                     GOTO 700
                  END IF
               END DO
 700           CONTINUE
C
C              *** Checking whether this component should ***
C              *** not be calculated, due to symmetry.    ***
C
               IF (FCCAL) THEN
                  KIDTMP = 1
                  KIDDBT = KIDTMP + NORDR
                  KIRPDG = KIDDBT + NORDR
                  KIRPST = KIRPDG + NORDR
                  KLAST  = KIRPST + NORDR
                  LWRK   = LWORK - KLAST + 1
                  CALL FCSCRN(GRIREP,WORK(KLAST),KDPMTX,INDSTP,ICRIRP,
     &                        IDUMMY,WORK(KIDTMP),WORK(KIDDBT),
     &                        WORK(KIRPDG),WORK(KIRPST),IDUMMY,LWORK,
     &                        NLDPMX,LDPMTX,IFRSTD,IORDR,IORDR-1,IDUMMY,
     &                        IDUMMY,IDUMMY,IPRINT,CLFVAL,.FALSE.,
     &                        .FALSE.,.TRUE.)
               ELSE
                  CLFVAL = .TRUE.
               END IF
C
C              *** In some calculations we only need diagonal. ***
C              *** derivatives.                                ***
C
               IF ((IORDR.GT.1).AND.DIAGON) THEN
                  CLFVAL = (CLFVAL).AND.
     &                     (INDSTP(IORDR-1).EQ.INDSTP(IORDR))
               END IF
C
C              *** The number of innermost elements. ***
C
               IF (FCCAL) THEN
                  NINNR2 = 1
                  IF (NAORDR .GE. 1) THEN
                     NINNR2 = INDSTP(IORDR)
                  END IF
                  IF (NAORDR .GE. 2) THEN
                     NINNR2 = NINNR2*(INDSTP(IORDR)+1)/2
                  END IF
               ELSE
                  NINNR2 = NFINNR
               END IF
C
               IF (CLFVAL) THEN
C
                  CALL IZERO(NCVAL,NDCOOR)
                  DO IC = 1, IORDR
                     NCVAL(INDSTP(IC)) = NCVAL(INDSTP(IC)) + 1
                  END DO
C
C                 *******************************************************
C                 *** IDCOMP -> Maks steporder to get the derivative  ***
C                 *** IDDCMP -> Counting array, maks to min steporder ***
C                 *** NTTYPE -> Number of function values needed for  ***
C                 ***           1 component                           ***
C                 *******************************************************
C
                  NTTYPE = 1
                  CALL IZERO(IDCOMP,NDCOOR)
                  CALL IZERO(IDDCMP,NDCOOR)
                  DO IC = 1, NDCOOR
                     IF (NCVAL(IC) .NE. 0) THEN
                        IDCOMP(IC) =   INT((NCVAL(IC)+1)/2)
                        NTTYPE     =   NTTYPE*(2*IDCOMP(IC) + 1)
                     END IF
                  END DO
                  DO I=1,NDCOOR
                     IDDCMP(I) = IDCOMP(I)
                  END DO
C
                  DO ITTYPE = 1, NTTYPE
C
C                    *** Finding the right indices to identify ***
C                    ***        the right function value       ***
C
                     DO IC = 1, NDCOOR
                        IF ((IDDCMP(IC) .GT. -IDCOMP(IC))
     &                              .AND. (ITTYPE .NE. 1)) THEN
                           IDDCMP(IC) = IDDCMP(IC) - 1
                           DO ICT = 1, IC-1
                              IDDCMP(ICT) = IDCOMP(ICT)
                           END DO
                           GOTO 800
                        END IF
                     END DO
 800                 CONTINUE
C
C                    *** Calculate the coefficient for this ***
C                    *** function value                     ***
C
                     BCOEFF = D1
                     NUMCOF = 0
                     DO IC = 1, NDCOOR
                        IF (NCVAL(IC) .NE. 0) THEN
                           BCOEFF = BCOEFF*COEFF(IDDCMP(IC),NCVAL(IC))
                           NUMCOF = NUMCOF + 1
                        END IF
                     END DO
                     IF (NUMCOF .EQ. 0) BCOEFF = D0
C
C                    *** Does the function value contribute? ***
C
                     IF (BCOEFF**2 .GT. DMIN) THEN
C
C                       ************************************************
C                       *** This subroutine finds the adress for the ***
C                       ***      function value, from the indices    ***
C                       ***   NEIND - The adress in the FUVAL-array  ***
C                       ************************************************
C
                        KITCMP = 1
c                        NINTIN = 1
c                        IF (FCCAL) NINTIN = NINNR2
                        CALL GTEIND(IADRSS,IDDCMP,NCVAL,WORK(KITCMP),
     &                              NEIND,ITTYPE,NORDR,IORDR)
C
C                       *** The derivative is calculated. ***
C
                        DO INNER = 1, NINNR2
                           ID = IDERIV + INNER
                           DERIV(ID) = DERIV(ID)
     &                               + BCOEFF*FUNVAL(INNER,NEIND)
                        END DO
C
                     END IF
                  END DO
C
C                 *** The derivative is correctly scaled. ***
C
                  SCLFCK = D1
                  DO INNER = 1, NINNR2
                     ID = IDERIV + INNER
                     DERIV(ID) = DERIV(ID)*DIVDIS*SCLFCK
                  END DO
               END IF
               IDERIV = IDERIV + NINNR2
            END DO
         END IF
C
 100  CONTINUE
C
      CALL QEXIT('NMNDER')
      RETURN
      END
C
C
C   /*Deck gteind*/
      SUBROUTINE GTEIND(IADRSS,INDCMP,NCVAL,ITCMP,NEIND,NMTYPE,NORDR,
     &                  IORDR)
C
C     *************************************************************************
C     *** This routine finds the adress of the function-value (in numerical ***
C     ***                differentiation) and returns it.                   ***
C     *** NEIND  -> The adress of the function value.                       ***
C     *** IADRSS -> Adress of the start of the steporder                    ***
C     *** NSTP   -> Number of function-component within the order           ***
C     *** NDISP  -> Number value within the component.                      ***
C     *************************************************************************
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <trkoor.h>
#include <numder.h>
      DIMENSION IADRSS(NMTYPE), INDCMP(NCOOR), ITCMP(NORDR),
     &          NCVAL(NCOOR)
C
      CALL IZERO(ITCMP,NORDR)
C
C     *** Preliminary components ***
C
      IRSRDR = 0
      NMX = 0
      ITOT = 0
      ITWOTT = 0
      DO 100 I = NDCOOR, 1, -1
         IF (INDCMP(I) .NE. 0) THEN
            IRSRDR = IRSRDR + 1
            ITOT = ITOT + ABS(INDCMP(I))
            ITWOTT = ITWOTT + (2*ABS(INDCMP(I))-1)
            ITCMP(IRSRDR) = I
            IF (ABS(INDCMP(I)) .GT. NMX) THEN
               NMX   = ABS(INDCMP(I))
               NICMX = I
            END IF
         END IF
 100  CONTINUE
      NRSRDR = IRSRDR
C
      ITYPE = 1
      DO 300 IRDR  = 1, ITWOTT-1
         IHORDR = INT((IRDR+1)/2)
         DO 400 IMXRDR = 1, IHORDR
            ITYPE = ITYPE + 1
 400     CONTINUE
 300  CONTINUE
      NTYPE  = ITYPE + NMX
C
      IF (NMX .EQ. 0) THEN
         ISTP = 1
C
      ELSE IF (NMX .EQ. 1) THEN
         ISTP = 0
         DO 500 IRS1 = 1, NRSRDR - 1
            ITISTP = 1
            DO 600 I  = 1, NRSRDR-IRS1+1
               ITISTP = ITISTP*(ITCMP(IRS1)-I)/I
 600        CONTINUE
            ISTP = ISTP + ITISTP
 500     CONTINUE
         ISTP = ISTP + ITCMP(NRSRDR)
C
C
C
      ELSE
         IF (ITOT .EQ. NMX) THEN
            ISTP = ITCMP(1)
         ELSE
            IF (ITCMP(1) .NE. NICMX) THEN
               ITMP1 = ITCMP(1)
               ITCMP(1) = NICMX
               DO 650 I = 2, NRSRDR
                  ITMP2    = ITCMP(I)
                  ITCMP(I) = ITMP1
                  ITMP1    = ITMP2
 650           CONTINUE
            END IF
C
            ISTP = 0
            ISTP = (ITCMP(1)-1)*(NDCOOR-1) + ITCMP(2)
            IF (NICMX .LT. ITCMP(2)) ISTP = ISTP - 1
C
            DO 700 IRS1 = 3, NRSRDR-1
               ITISTP = 1
               DO 800 I  = 1, NRSRDR-IRS1+1
                  ITISTP = ITISTP*(ITCMP(IRS1)-I)/I
 800           CONTINUE
               ISTP = ISTP + ITISTP
 700        CONTINUE
            IF (NRSRDR .GT. 2) ISTP = ISTP + ITCMP(NRSRDR)
         END IF
      END IF
      NINNER = 2**NRSRDR
      NSTP = (ISTP-1)*NINNER
C
      NDISP = 1
      IORD = NRSRDR-1
      DO 900 IRS = NRSRDR, 2, -1
         IDISP = 0
         IF (INDCMP(ITCMP(IRS)) .LT. 0) THEN
            IDISP = 2**(IORD)
         END IF
         IORD = IORD - 1
         NDISP = NDISP + IDISP
 900  CONTINUE
      IF (ITCMP(1) .NE. 0) THEN
         IF (INDCMP(ITCMP(1)) .LE. 0) NDISP = NDISP + 1
      END IF
C
      NEIND = IADRSS(NTYPE) + (NSTP + NDISP)
C
      RETURN
      END
C
C    /*Deck pritdr*/
      SUBROUTINE PRITDR(TMPTDR,SYMCOR,TDER,SYMTDR,NTMPDM,NUMCOR,LTXT,
     &                  IPRINT,PRWHLE,TEXT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
      PARAMETER (KCOL=6)
#include <trkoor.h>
#include <symmet.h>
#include <numder.h>
#include <fcsym.h>
#include <cbinum.h>
      CHARACTER*80 TEXT
      LOGICAL PRWHLE
      DIMENSION TDER(NCOOR,NCOOR,NCOOR), TMPTDR(NTMPDM),
     &          SYMTDR(NCOOR,NCOOR,NCOOR), SYMCOR(NCOOR,NCOOR)
C
C     *****************************************
C     *** Assigning values to proper places ***
C     ***  according to permutational sym.  ***
C     *****************************************
C
      ITMP = 0
      DO 100 K = 1, NUMCOR
      DO 100 J = 1, K
      DO 100 I = 1, J
         ITMP = ITMP + 1
         SYMTDR(I,J,K) = TMPTDR(ITMP)
         SYMTDR(I,K,J) = TMPTDR(ITMP)
         SYMTDR(J,I,K) = TMPTDR(ITMP)
         SYMTDR(J,K,I) = TMPTDR(ITMP)
         SYMTDR(K,I,J) = TMPTDR(ITMP)
         SYMTDR(K,J,I) = TMPTDR(ITMP)
 100  CONTINUE
C
C     **************************************
C     *** Printing symmetric coordinates ***
C     ***    Cartesian if no symmetry    ***
C     **************************************
C
      IF (PRWHLE) THEN
         IF (.NOT.MINOUT) THEN
            CALL PRTDER(SYMTDR,NCOOR,NUMCOR,TEXT,LTXT,IPRINT)
         ELSE
            WRITE (LUPRI,'(A/)') 
     *         " Output of third derivative suppressed"
         ENDIF
      ELSE
         CALL HEADER('Diagonal of cubic force field, F(I,J,J)',-1)
C
         ISTRT = 1
         LAST  = MIN(NDCOOR,KCOL)
         KCOOR = NDCOOR
         NCOL  = NDCOOR/KCOL
         IF (MOD(NDCOOR,KCOL).NE.0) NCOL = NCOL + 1
C
         DO ICOL = 1, NCOL
            DO ICOOR = 1, NDCOOR
               WRITE (LUPRI,'(5X,6F12.6)')
     &                        (SYMTDR(ICOOR,I,I),I=ISTRT,LAST)
            END DO
            WRITE (LUPRI,'(A)') '                                  '
            ISTRT = ISTRT + KCOL
            LAST  = MIN(NDCOOR,KCOL+LAST)
         END DO
      END IF
C
C     ************************************
C     *** Transformation to cartesian  ***
C     ***  coordinstes, and printing.  ***
C     ************************************
C
      IF ((FCLASS(1:3) .NE. 'C1 ').AND.(TEXT(1:6).NE.'normal')) THEN
         LTXT = 9
         TEXT(1:9) = 'cartesian'
C
         CALL TRATDR(SYMCOR,SYMTDR,TDER,NCOOR,NCOOR,NCOOR,TEXT,LTXT,
     &               IPRINT)
C
         IF (.NOT.MINOUT)
     &            CALL PRTDER(TDER,NCOOR,NUMCOR,TEXT,LTXT,IPRINT)
      ELSE
         CALL DCOPY(NCOOR**3,SYMTDR,1,TDER,1)
      END IF
C
      RETURN
      END

C
C    /*Deck tratdr*/
      SUBROUTINE TRATDR(TRCOOR,CR1TDR,CR2TDR,NMCOR1,NMCOR2,NCOOR,TEXT,
     &                  LTXT,IPRINT)
C     ***********************************************************
C     *** Transforming a third derivative into another set of ***
C     *** coordinates. Tracor is the transformation matrix.   ***
C     ***********************************************************
#include <implicit.h>
#include <priunit.h>
      CHARACTER*80 TEXT
      DIMENSION CR1TDR(NCOOR,NCOOR,NCOOR), CR2TDR(NCOOR,NCOOR,NCOOR),
     &          TRCOOR(NCOOR,NCOOR)
C
      CALL DZERO(CR2TDR,NCOOR**3)
C
      DO 100 ICR1C3 = 1, NMCOR1
      DO 100 ICR1C2 = 1, NMCOR1
      DO 100 ICR1C1 = 1, NMCOR1
      DO 100 ICR2C1 = 1, NMCOR2
         CR2TDR(ICR2C1,ICR1C2,ICR1C3) = CR2TDR(ICR2C1,ICR1C2,ICR1C3)
     &              + TRCOOR(ICR2C1,ICR1C1)*CR1TDR(ICR1C1,ICR1C2,ICR1C3)
 100  CONTINUE
C
      CALL DZERO(CR1TDR,NCOOR**3)
      DO 200 ICR1C3 = 1, NMCOR1
      DO 200 ICR1C2 = 1, NMCOR1
      DO 200 ICR2C2 = 1, NMCOR2
      DO 200 ICR2C1 = 1, NMCOR2
         CR1TDR(ICR2C1,ICR2C2,ICR1C3) = CR1TDR(ICR2C1,ICR2C2,ICR1C3)
     &              + TRCOOR(ICR2C2,ICR1C2)*CR2TDR(ICR2C1,ICR1C2,ICR1C3)
 200  CONTINUE
C
      CALL DZERO(CR2TDR,NCOOR**3)
      DO 300 ICR1C3 = 1, NMCOR1
      DO 300 ICR2C3 = 1, NMCOR2
      DO 300 ICR2C2 = 1, NMCOR2
      DO 300 ICR2C1 = 1, NMCOR2
         CR2TDR(ICR2C1,ICR2C2,ICR2C3) = CR2TDR(ICR2C1,ICR2C2,ICR2C3)
     &              + TRCOOR(ICR2C3,ICR1C3)*CR1TDR(ICR2C1,ICR2C2,ICR1C3)
 300  CONTINUE
C
      IF (IPRINT .GT.6) THEN
         WRITE (LUPRI,'(A)') 'Transformation tensor: '
         CALL PRTRMA(TRCOOR,NCOOR,NCOOR,NMCOR2,NMCOR1,LUPRI)
         CALL PRTDER(CR2TDR,NCOOR,NMCOR2,TEXT,LTXT,IPRINT)
      END IF
C
      RETURN
      END
C
C    /*Deck prifdr*/
      SUBROUTINE PRIFDR(TMPFDR,SYMCOR,FDER,SYMFDR,NTMPDM,NUMCOR,LTXT,
     &                  IPRINT,TEXT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
C
#include <trkoor.h>
#include <symmet.h>
#include <numder.h>
#include <fcsym.h>
#include <cbinum.h>
      CHARACTER*80 TEXT
      DIMENSION FDER(NCOOR,NCOOR,NCOOR,NCOOR), TMPFDR(NTMPDM),
     &          SYMFDR(NCOOR,NCOOR,NCOOR,NCOOR), SYMCOR(NCOOR,NCOOR)
C
C     *****************************************
C     *** Assigning values to proper places ***
C     ***  according to permutational sym.  ***
C     *****************************************
C
      ITMP = 0
      DO 100 L = 1, NUMCOR
      DO 100 K = 1, L
      DO 100 J = 1, K
      DO 100 I = 1, J
         ITMP = ITMP + 1
C
         SYMFDR(I,J,K,L) = TMPFDR(ITMP)
         SYMFDR(I,J,L,K) = TMPFDR(ITMP)
         SYMFDR(I,K,J,L) = TMPFDR(ITMP)
         SYMFDR(I,K,L,J) = TMPFDR(ITMP)
         SYMFDR(I,L,J,K) = TMPFDR(ITMP)
         SYMFDR(I,L,K,J) = TMPFDR(ITMP)
C
         SYMFDR(J,I,K,L) = TMPFDR(ITMP)
         SYMFDR(J,I,L,K) = TMPFDR(ITMP)
         SYMFDR(J,K,I,L) = TMPFDR(ITMP)
         SYMFDR(J,K,L,I) = TMPFDR(ITMP)
         SYMFDR(J,L,I,K) = TMPFDR(ITMP)
         SYMFDR(J,L,K,I) = TMPFDR(ITMP)
C
         SYMFDR(K,I,J,L) = TMPFDR(ITMP)
         SYMFDR(K,I,L,J) = TMPFDR(ITMP)
         SYMFDR(K,J,I,L) = TMPFDR(ITMP)
         SYMFDR(K,J,L,I) = TMPFDR(ITMP)
         SYMFDR(K,L,I,J) = TMPFDR(ITMP)
         SYMFDR(K,L,J,I) = TMPFDR(ITMP)
C
         SYMFDR(L,I,J,K) = TMPFDR(ITMP)
         SYMFDR(L,I,K,J) = TMPFDR(ITMP)
         SYMFDR(L,J,I,K) = TMPFDR(ITMP)
         SYMFDR(L,J,K,I) = TMPFDR(ITMP)
         SYMFDR(L,K,I,J) = TMPFDR(ITMP)
         SYMFDR(L,K,J,I) = TMPFDR(ITMP)
 100  CONTINUE
C
C     **************************************
C     *** Printing symmetric coordinates ***
C     ***    Cartesian if no symmetry    ***
C     **************************************
C
      IF (.NOT.MINOUT) THEN
         CALL PRFDER(SYMFDR,NCOOR,NUMCOR,TEXT,LTXT,IPRINT)
      ELSE
         WRITE (LUPRI,'(A/)') 
     *      " Output of fourth derivative suppressed"
      ENDIF
C
C     ************************************
C     *** Transformation to cartesian  ***
C     ***  coordinstes, and printing.  ***
C     ************************************
C
      IF ((FCLASS(1:3).NE.'C1 ').AND.(TEXT(1:6).NE.'normal')) THEN
C
C        *** Coordinate transformation. ***
C
         CALL TRAFDR(SYMCOR,SYMFDR,FDER,NUMCOR,NUMCOR,NCOOR,TEXT,
     &               LTXT,IPRINT)
C
C        *** Printing in cartesian coordinates. ***
C
         IF (.NOT.MINOUT) 
     &      CALL PRFDER(FDER,NCOOR,NUMCOR,'cartesian',9,IPRINT)
      ELSE
         CALL DCOPY(NCOOR**4,SYMFDR,1,FDER,1)
      END IF
C
      RETURN
      END
C
C    /*Deck trafdr*/
      SUBROUTINE TRAFDR(TRCOOR,CR1FDR,CR2FDR,NMCOR1,NMCOR2,NCOOR,TEXT,
     &                  LTXT,IPRINT)
C     **********************************************************
C     *** Transforming quartic force field to another set of ***
C     *** coordinates                                        ***
C     **********************************************************
#include <implicit.h>
#include <priunit.h>
      CHARACTER*80 TEXT
      DIMENSION CR1FDR(NCOOR,NCOOR,NCOOR,NCOOR),
     &          CR2FDR(NCOOR,NCOOR,NCOOR,NCOOR), TRCOOR(NCOOR,NCOOR)
C
      KDIM = NCOOR**4
C
      CALL DZERO(CR2FDR,KDIM)
      DO 100 ICR1C4 = 1, NMCOR1
      DO 100 ICR1C3 = 1, NMCOR1
      DO 100 ICR1C2 = 1, NMCOR1
      DO 100 ICR1C1 = 1, NMCOR1
      DO 100 ICR2C1 = 1, NMCOR2
         CR2FDR(ICR2C1,ICR1C2,ICR1C3,ICR1C4) =
     &         CR2FDR(ICR2C1,ICR1C2,ICR1C3,ICR1C4)
     &       + TRCOOR(ICR2C1,ICR1C1)*CR1FDR(ICR1C1,ICR1C2,ICR1C3,ICR1C4)
 100  CONTINUE
C
      CALL DZERO(CR1FDR,KDIM)
      DO 200 ICR1C4 = 1, NMCOR1
      DO 200 ICR1C3 = 1, NMCOR1
      DO 200 ICR1C2 = 1, NMCOR1
      DO 200 ICR2C2 = 1, NMCOR2
      DO 200 ICR2C1 = 1, NMCOR2
         CR1FDR(ICR2C1,ICR2C2,ICR1C3,ICR1C4) =
     &         CR1FDR(ICR2C1,ICR2C2,ICR1C3,ICR1C4)
     &       + TRCOOR(ICR2C2,ICR1C2)*CR2FDR(ICR2C1,ICR1C2,ICR1C3,ICR1C4)
 200  CONTINUE
C
      CALL DZERO(CR2FDR,KDIM)
      DO 300 ICR1C4 = 1, NMCOR1
      DO 300 ICR1C3 = 1, NMCOR1
      DO 300 ICR2C3 = 1, NMCOR2
      DO 300 ICR2C2 = 1, NMCOR2
      DO 300 ICR2C1 = 1, NMCOR2
         CR2FDR(ICR2C1,ICR2C2,ICR2C3,ICR1C4) =
     &         CR2FDR(ICR2C1,ICR2C2,ICR2C3,ICR1C4)
     &       + TRCOOR(ICR2C3,ICR1C3)*CR1FDR(ICR2C1,ICR2C2,ICR1C3,ICR1C4)
 300  CONTINUE
C
      CALL DZERO(CR1FDR,KDIM)
      DO 400 ICR1C4 = 1, NMCOR1
      DO 400 ICR2C4 = 1, NMCOR2
      DO 400 ICR2C3 = 1, NMCOR2
      DO 400 ICR2C2 = 1, NMCOR2
      DO 400 ICR2C1 = 1, NMCOR2
         CR1FDR(ICR2C1,ICR2C2,ICR2C3,ICR2C4) =
     &        CR1FDR(ICR2C1,ICR2C2,ICR2C3,ICR2C4)
     &       + TRCOOR(ICR2C4,ICR1C4)*CR2FDR(ICR2C1,ICR2C2,ICR2C3,ICR1C4)
 400  CONTINUE
C
      CALL DCOPY(KDIM,CR1FDR,1,CR2FDR,1)
C
      IF (IPRINT .GT. 7) THEN
         CALL PRFDER(CR2FDR,NCOOR,NMCOR2,TEXT,LTXT,IPRINT)
      END IF
C
      RETURN
      END
C
C
C
      SUBROUTINE HARMAN(SYMCOR,BKPHES,TRAMAT,TMPHES,WORK,NCOOR,LWORK,
     &                  IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <taymol.h>
      DIMENSION BKPHES(MXCOOR,MXCOOR), SYMCOR(NCOOR,NCOOR),
     &          TMPHES(NCOOR ,NCOOR ), TRAMAT(NCOOR,NCOOR),
     &          WORK(LWORK)
C
C     *** Backing up hessian for further use after ***
C     *** harmonic analysis.                       ***
C
      KDIM = MXCOOR**2
      CALL DCOPY(KDIM,HESMOL,1,BKPHES,1)
C
C     *** Transforming the hessian matrix to cartesian ***
C     *** coordinates.                                 ***
C
      CALL OTRTEN(HESMOL,SYMCOR,TMPHES,MXCOOR,NCOOR,NCOOR,IPRINT,'N',
     &            'T')
C
C     *** Transforming the hessian matrix to dalton ***
C     *** symmetry coordinates.                     ***
C
C     *** Transformation matrix. ***
C
      ITYPE = 1
      KTEST = 1
      CALL TRACOR(TRAMAT,WORK(KTEST),ITYPE,NCOOR,IPRINT)
C
C     *** Transformation. ***
C
      CALL OTRTEN(HESMOL,TRAMAT,TMPHES,MXCOOR,NCOOR,NCOOR,IPRINT,'N',
     &            'T')
C
C     *** Run harmonic analysis. ***
C
      KSTART = 1
      CALL VIBCTL(WORK(1),LWORK)
C
C     *** Restoring old hessian. ***
C
      KDIM = MXCOOR**2
      CALL DCOPY(KDIM,BKPHES,1,HESMOL,1)
C
      RETURN
      END
C
C
C     /* Deck mknrmc */
      SUBROUTINE MKNRMC(SYMCOR,CSTART,TRNCCR,TRAMSS,EIGNVL,EGNVCT,
     &                  HESMWT,TM1TMP,TM2TMP,AMASS,DKIN,HTESTM,FREQ,
     &                  RNNORM,CORTMP,WORK,ICRIRP,NATTYP,NMSYSP,
     &                  LWORK,IPRINT)
************************************************************
*** Makes normal coordinates from the molecular hessian, ***
*** and writes out the  harmonic frequensies according   ***
***             general symmetry species.                ***
************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <codata.h>
      PARAMETER (DMTHR = 2.0D-8, D0 = 0.0D0, D1 = 1.0D0)
#include <taymol.h>
#include <trkoor.h>
#include <nuclei.h>
#include <cbinum.h>
#include <numder.h>
#include <abainf.h>
#include <dummy.h>
      CHARACTER*6 TXT
      DIMENSION SYMCOR(NCOOR,NCOOR), EIGNVL(NCOOR), EGNVCT(NCOOR,NCOOR),
     &          HESMWT(NCOOR*(NCOOR+1)/2), AMASS(NATOMS), CSTART(NCOOR),
     &          DKIN(NCOOR), TM1TMP(NCOOR,NCOOR), TRAMSS(NCOOR),
     &          TM2TMP(NCOOR,NCOOR), TRNCCR(NCOOR,NCOOR),
     &          FREQ(NCOOR), HTESTM(NCOOR,NCOOR), RNNORM(NCOOR),
     &          CORTMP(NCOOR), WORK(LWORK), ENORMN(NCOOR)
      DIMENSION ICRIRP(NCOOR,2), NATTYP(NATOMS), NMSYSP(NCOOR)
      INTEGER BEGIN
C
C     *** Calculating center of mass, and mass of each center ***
C
      CALL DCOPY(NCOOR,CSTART,1,CORTMP,1)
      CALL CMMASS(CORTMP,AMASS,NATTYP,WORK,IPRINT)
C
C     *** Diagonal sqrt(mass)^(-1/2) matrix ***
C
      DO 200 IC = 1, NCOOR
         DKIN(IC) = D1/SQRT(XFAMU*AMASS((IC+2)/3))
 200  CONTINUE
C
C     *** The (mass)^(-1/2) matrix for symmetry coordinates. ***
C
      DO 300 IC2 = 1, NCOOR
      DO 300 IC1 = 1, NCOOR
         TM1TMP(IC1,IC2) = DKIN(IC1)*SYMCOR(IC1,IC2)
 300  CONTINUE
C
C     *** TM2TMP is the (mass)^(-1/2) matrix. ***
C
      KDIM = NCOOR**2
      CALL DZERO(TM2TMP,KDIM)
      DO 400 IC3 = 1, NCOOR
      DO 400 IC2 = 1, NCOOR
      DO 400 IC1 = 1, NCOOR
         TM2TMP(IC1,IC3) = TM2TMP(IC1,IC3)
     &                   + SYMCOR(IC2,IC1)*TM1TMP(IC2,IC3)
 400  CONTINUE
C
C     *** Test if TM2TMP is a diagonal matrix (if there are ***
C     ***              different isotopes)                  ***
C
      IF (HTEST) THEN
         DO 500 IC2 = 1, NCOOR
         DO 500 IC1 = 1, NCOOR
            IF ((IC1 .NE. IC2).AND.(ABS(TM2TMP(IC1,IC2)).GT.DMTHR))
     &            CALL QUIT('Diagonal mass test failed. Off-diagonal' //
     &                      'elements present.')
 500     CONTINUE
      END IF
C
C     *** Mass transformation matrix ***
C
      DO 600 IC = 1, NCOOR
         TRAMSS(IC) = D1/TM2TMP(IC,IC)
 600  CONTINUE
C
C     *** Calculating the mass-weighted Hessian ***
C
      CALL DZERO(TM1TMP,KDIM)
      DO 700 IC3 = 1, NCOOR
      DO 700 IC2 = 1, NCOOR
      DO 700 IC1 = 1, NCOOR
         TM1TMP(IC1,IC3) = TM1TMP(IC1,IC3)
     &                   + TM2TMP(IC1,IC2)*HESMOL(IC2,IC3)
 700  CONTINUE
C
C     *** HESMWT is the mass-weighted hessian. ***
C
      IC12 = 0
      CALL DZERO(HESMWT,NCOOR*(NCOOR+1)/2)
      DO 800 IC2 = 1, NCOOR
      DO 800 IC1 = 1, IC2
         IC12 = IC12 + 1
         HESMWT(IC12) = HESMWT(IC12) + TM1TMP(IC1,IC2)*TM2TMP(IC2,IC2)
 800  CONTINUE
C
C     *** Test to check if mass-weighted hessian is symmetric. ***
C
      IF (HTEST) THEN
         CALL DZERO(HTESTM,KDIM)
         DO 900 IC3 = 1, NCOOR
         DO 900 IC2 = 1, NCOOR
         DO 900 IC1 = 1, NCOOR
            HTESTM(IC1,IC3) = HTESTM(IC1,IC3)
     &                      + TM1TMP(IC1,IC2)*TM2TMP(IC2,IC3)
 900     CONTINUE
C
         IF (.NOT. RESTRT) THEN
            DO 1100 IC2 = 1, NCOOR
            DO 1100 IC1 = 1, IC2
               IF (ABS(HTESTM(IC1,IC2)-HTESTM(IC2,IC1)).GT.DMTHR)
     &              CALL QUIT('Mass-weighted hessian is not symmetric.')
 1100       CONTINUE
         END IF
C
         CALL HEADER('Mass weighted Hessian in symmetry coordinates',-1)
         NUMTIM = (NCOOR-1)/6 + 1
         DO 1200 ITIM = 1, NUMTIM
            ISTART =     6*(ITIM-1) + 1
            IEND   = MIN(6* ITIM  ,NCOOR)
            DO 1300 IC1 = 1, NCOOR
               WRITE(LUPRI,'(6F17.14)')(HTESTM(IC1,IC2),IC2=ISTART,IEND)
 1300       CONTINUE
            WRITE (LUPRI,'(A)') '                                  '
 1200    CONTINUE
      END IF
C
C     *** Diagonalizing the mass weighted Hessian. ***
C
      KWRK  = 1
      KIWRK = KWRK + NCOOR
      CALL DZERO(EGNVCT,NCOOR**2)
      CALL DUNIT(EGNVCT,NCOOR)
      CALL JACO(HESMWT,EGNVCT,NCOOR,NCOOR,NCOOR,WORK(KWRK),WORK(KIWRK))

C
C     *** Storing the transformation matrix for later isotope-studies. ****
C
      DO 1400 IC2 = 1, NCOOR
      DO 1400 IC1 = 1, NCOOR
         TRNCCR(IC1,IC2) = EGNVCT(IC1,IC2)
 1400 CONTINUE
C
C     *** Mass-weighting the normal coordinates. ***
C
      KDIM = NCOOR**2
      CALL DCOPY(KDIM,EGNVCT,1,TM1TMP,1)
      CALL DZERO(EGNVCT,KDIM)
      DO 1500 IC3 = 1, NCOOR
      DO 1500 IC2 = 1, NCOOR
      DO 1500 IC1 = 1, NCOOR
         EGNVCT(IC1,IC3) = EGNVCT(IC1,IC3)
     &                   + TM2TMP(IC1,IC2)*TM1TMP(IC2,IC3)
 1500 CONTINUE
C
C     *** Normalizing the normal and transformation coordinates. ***
C
      CALL DZERO(RNNORM,NCOOR)
      DO 1600 IC2 = 1, NCOOR
         RLENGT2 = D0
         DO 1700 IC1 = 1, NCOOR
            RNNORM(IC2) = RNNORM(IC2) + EGNVCT(IC1,IC2)**2
            RLENGT2     = RLENGT2     + TRNCCR(IC1,IC2)**2
 1700    CONTINUE
         RNNORM(IC2) = SQRT(RNNORM(IC2))
C
         DRINV1 = D1/RNNORM(IC2)
         DRINV2 = D1/SQRT(RLENGT2)
         DO 1800 IC1 = 1, NCOOR
            EGNVCT(IC1,IC2) = EGNVCT(IC1,IC2)*DRINV1
            TRNCCR(IC1,IC2) = TRNCCR(IC1,IC2)*DRINV2
 1800    CONTINUE
 1600 CONTINUE
C
C     *** Removing the redundant normal coordinates, and setting some ***
C     *** common variables according to this                          ***
C
      IC12 = 0
      IFREQ = 0
      NUMZRO = 0
      DO 1900 IC = 1, NCOOR
         IC12 = IC12 + IC
         IF (ABS(HESMWT(IC12)).GT.DMTHR) THEN
            IFREQ = IFREQ + 1
            FREQ(IFREQ) = SQRT(ABS(HESMWT(IC12)))
         ELSE
            NUMZRO = NUMZRO + 1
            DO 2100 IC2 = IC-NUMZRO+1, NCOOR-1
               ICRIRP(IC2,1) = ICRIRP(IC2+1,1)
               ICRIRP(IC2,2) = ICRIRP(IC2+1,2)
               DO 2200 IC1 =  1, NCOOR
                  RNNORM(    IC2) = RNNORM(    IC2+1)
                  EGNVCT(IC1,IC2) = EGNVCT(IC1,IC2+1)
                  TRNCCR(IC1,IC2) = TRNCCR(IC1,IC2+1)
 2200          CONTINUE
C
 2100       CONTINUE
C
         END IF
 1900 CONTINUE
      NDCOOR = NCOOR - NUMZRO
C
C     *** Calculating the normal coordinates in cartesian coordinates ***
C
      CALL DZERO(TM1TMP,KDIM)
      DO 2300 IC3 = 1, NCOOR-NUMZRO
      DO 2300 IC2 = 1, NCOOR
      DO 2300 IC1 = 1, NCOOR
         TM1TMP(IC1,IC3) = TM1TMP(IC1,IC3)
     &                   + SYMCOR(IC1,IC2)*EGNVCT(IC2,IC3)
 2300 CONTINUE
C
C     *** Using these coordinates in future differentiations. ***
C
      DO 2400 IC2 = 1, NCOOR-NUMZRO
      DO 2400 IC1 = 1, NCOOR
         SYMCOR(IC1,IC2) = TM1TMP(IC1,IC2)
 2400 CONTINUE
C
C     *** Printing the frequencies. ***
C
      WRITE (LUPRI, '(10X,A,I3)')
     &                  'Number of modes with zero frequency: ', NUMZRO
      CALL HEADER('Vibrational frequencies harmonic approximation: ',-1)
      WRITE (LUPRI,'(20X,A)') ' Mode     cm-1       hartrees  '
      DO 2500 IFREQ = 1, NCOOR-NUMZRO
         WRITE (LUPRI,'(20X,I4,F12.2,F12.6)') ICRIRP(IFREQ,1),
     &                                XTKAYS*FREQ(IFREQ), FREQ(IFREQ)
 2500 CONTINUE
C
C     *** Printing the cartesian components of the normal coordinates ***
C
      NONZRO = NCOOR-NUMZRO
      CALL HEADER('Normal coordinates:  ',0)
      NUMTIM = (NONZRO-1)/6 + 1
      DO ITIM = 1, NUMTIM
         ISTART =     6*(ITIM-1)+1
         IEND   = MIN(6*ITIM,NONZRO)
         WRITE (LUPRI,'(I11,8I13)') (ICRIRP(I,1),I=ISTART,IEND)
         DO ICOOR = 1, NCOOR
            WRITE (LUPRI,'(6F13.7)') (SYMCOR(ICOOR,I),I=ISTART,IEND)
         END DO
         WRITE (LUPRI,'(A)') '                                        '
      END DO

C     We punch out harmonic freqs and normal coordinates on file
C     DALTON.NOR for use in Midas Vibrational calculations

      LUNOR = -1
      CALL GPOPEN(LUNOR,'DALTON.NOR','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      WRITE(LUNOR,'(A)') 'Harmonic Freqs. in cm^-1'
      DO 911 IMODE = 1, NONZRO
        WRITE(LUNOR,'(E23.16)') XTKAYS*FREQ(IMODE)
 911  CONTINUE
      WRITE(LUNOR,*)


C     The normal coordinates in SYMCOR is normalized to one in cartesian (x) space. We
C     would like the normal coordinates to be normalized to one in q-space, q = sqrt(m)*x

      DO 912 IMODE = 1,NONZRO
        ENORM2 = 0.0D0
        DO 913 I = 1, NCOOR
          ENORM2 = ENORM2 + (SYMCOR(I,IMODE)**2)/(DKIN(I)*DKIN(I)*XFAMU)
 913    CONTINUE
      ENORMN(IMODE) = D1/SQRT(ENORM2)
 912  CONTINUE

      WRITE(LUNOR,'(A)') 'Normal Coordinates'
      DO 914 IMODE = 1, NONZRO
        WRITE(LUNOR,8041) (ENORMN(IMODE)*SYMCOR(I,IMODE),I=1,NCOOR)
        WRITE(LUNOR,*)
 914  CONTINUE

      WRITE(LUNOR,'(A)') 'Norm of Vectors'
      DO 915 IMODE = 1, NONZRO
        ENORM2 = 0.0D0
        DO 916 I = 1, NCOOR
          ENORM2 = ENORM2 + (ENORMN(IMODE)*SYMCOR(I,IMODE))**2
 916    CONTINUE
      WRITE(LUNOR,'(E23.16)') SQRT(ENORM2)
 915  CONTINUE

      CALL GPCLOSE(LUNOR,'KEEP')

C     
C     *** Writing to spectro file if requested. ***
C
      IF (SPECTR) THEN
         NTIME = 1
         IF (NRMCRD) THEN
            TXT  = 'normal'
         ELSE
            TXT  = 'cartes'
         END IF
         CALL WRISPC(FREQ,RNNORM,DUMMY,DUMMY,TXT,NCOOR,NDCOOR,NTIME,
     &               IPRINT)
      END IF
      IF (MIDAS) THEN
         NTIME = 1
         IF (NRMCRD) THEN
            TXT  = 'normal'
         ELSE
            TXT  = 'cartes'
         END IF
         IF (NRMCRD) CALL WRIMOP(FREQ,RNNORM,DUMMY,DUMMY,TXT,NCOOR,
     &                           NDCOOR,NTIME,IPRINT)
      END IF
C
C     *** Test printing. ***
C
      IF (IPRINT .GE. 20) THEN
         CALL HEADER('Eigenvectors of the symmetry adapted hessian',-1)
         DO IC1 = 1, NCOOR
            WRITE (LUPRI,'(10X,12F9.4)') (EGNVCT(IC1,IC2),IC2=1,NCOOR)
         END DO
C
         CALL HEADER ('Atomic masses used',-1)
         WRITE (LUPRI,'(3X,A,10I9)') 'Atom number:', (I,I=1,NATOMS)
         WRITE (LUPRI,'(17X, 10F9.4 )') (AMASS(I), I = 1, NATOMS)
C
         CALL HEADER('Diagonal elements of (sqrt(mass))^-1 matrix',-1)
         WRITE (LUPRI,'(24F9.4)') (DKIN(IC), IC=1,NCOOR)
C
         CALL HEADER('Transformed sqrt(mass)^-1 matrix',-1)
         DO IC1 = 1, NCOOR
            WRITE (LUPRI,'(24F9.7)') (TM2TMP(IC1,IC2),IC2= 1, NCOOR)
         END DO
C
C        *** Mass-weighted hessian. ***
C
         CALL HEADER('Diagonalized mass weighted hessian.',-1)
C
         IJ = 0
         DO J = 1, NCOOR
         DO I = 1, J
            IJ = IJ + 1
            TM2TMP(I,J) = HESMWT(IJ)
            TM2TMP(J,I) = HESMWT(IJ)
         END DO
         END DO
C
         BEGIN = 1
         KCOL  = 6
         LAST  = MIN(NCOOR,KCOL)
         NCOL  = NCOOR/KCOL
         IF (MOD(NCOOR,KCOL).NE.0) NCOL = NCOL + 1
C
         DO ICOL = 1, NCOL
            WRITE (LUPRI,1000) (ICRIRP(I,1),I = BEGIN,LAST)
C
            DO ICOOR = BEGIN, NCOOR
               WRITE (LUPRI,2000) ICRIRP(ICOOR,1),
     &                        (TM2TMP(ICOOR,I),I=BEGIN,MIN(LAST,ICOOR))
            END DO
            WRITE (LUPRI,'(A)') '                                 '
            BEGIN = BEGIN + KCOL
            LAST  = MIN(NCOOR,KCOL+LAST)
         END DO
 1000    FORMAT (8X,6(3X,I4,5X),(3X,I4,5X))
 2000    FORMAT (1X,I4,2X,6F12.6)
 8041    FORMAT(1P,3E23.16)
      END IF
C
      RETURN
      END
C
C
C     /*Deck trafrc*/
      SUBROUTINE TRAFRC(TDER,FDER,HESNRM,CORNRM,CRTNRM,SYCART,TNRMDR,
     &                  FNRMDR,WORK,NCOOR,NDIMF,NDIMT,LWORK,IPRINT)
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
#include <maxorb.h>
C
#include <infpar.h>
#include <taymol.h>
#include <numder.h>
#include <cbinum.h>
      CHARACTER*80 TEXT
      LOGICAL PRWHLE
      DIMENSION FDER(NDIMF), TDER(NDIMT), HESNRM(NCOOR,NCOOR),
     &          CORNRM(NCOOR,NCOOR), CRTNRM(NCOOR,NCOOR),
     &          SYCART(NCOOR,NCOOR), TNRMDR(NCOOR,NCOOR,NCOOR),
     &          FNRMDR(NCOOR,NCOOR,NCOOR,NCOOR), WORK(LWORK)
C
C     *** Transformation of hessians. ***
C
      IF (NAORDR+NMORDR.GE.2) THEN
         DO 100 J = 1, NCOOR
         DO 100 I = 1, J
            HESNRM(I,J) = HESMOL(I,J)
            HESNRM(J,I) = HESMOL(I,J)
 100     CONTINUE
C
C        *** Transforming hessian to cartesian coordinates. ***
C
         KTMPHS = 1
         CALL OTRTEN(HESNRM,SYCART,WORK(KTMPHS),NCOOR,NCOOR,NCOOR,
     &        IPRINT,'N','T')
C
C        *** Then to normal coordinates. ***
C
         CALL OTRTEN(HESNRM,CORNRM,WORK(KTMPHS),NCOOR,NCOOR,NDCOOR,
     &               IPRINT,'T','N')
C
C        *** Printing hessian. ***
C
         CALL HEADER('Hessian in normal coordinates',-1)
         CALL PRTRMA(HESNRM,NCOOR,NCOOR,NDCOOR,NDCOOR,LUPRI)
      END IF
C
C     *** Setting up the transformation (cartesian -> normal) ***
C     *** matrix needed for cubic and quartic force fields.   ***
C
      IF (NAORDR+NMORDR.GE.3) THEN
         DO 200 J = 1, NCOOR
         DO 200 I = 1, NCOOR
            CRTNRM(I,J) = CORNRM(J,I)
 200     CONTINUE
      END IF
C
C     *** Transformation of cubic force field. ***
C
      IF (NAORDR+NMORDR.GE.3) THEN
C
C        *** Transformation to cartesian coordinates and  ***
C        *** printing. Force field is returned in TNRMDR. ***
C
         PRWHLE = .TRUE.
         KSYMTD = 1
         LTXT = 8
         TEXT(1:8) = 'Symmetry'
         CALL PRITDR(TDER,SYCART,TNRMDR,WORK(KSYMTD),NDIMT,NCOOR,LTXT,
     &               IPRINT,PRWHLE,TEXT)
C
C        *** Transforming to normal coordinates. ***
C
         LTXT = 6
         TEXT(1:6) = 'normal'
C
         KTMPTD = 1
         CALL DCOPY(NCOOR**3,TNRMDR,1,WORK(KTMPTD),1)
         CALL TRATDR(CRTNRM,WORK(KTMPTD),TNRMDR,NCOOR,NDCOOR,NCOOR,
     &               TEXT,LTXT,IPRINT)
C
         IF (.NOT.MINOUT)
     &            CALL PRTDER(TNRMDR,NCOOR,NDCOOR,TEXT,LTXT,IPRINT)
      END IF
C
C     *** Transformation of quartic force field. ***
C
      IF (NAORDR+NMORDR.GE.4) THEN
C
C        *** Transformation to cartesian coordinates and  ***
C        *** printing. Force field is returned in FNRMDR. ***
C
         KSYMTD = 1
         LTXT = 8
         TEXT(1:8) = 'Symmetry'
         CALL PRIFDR(FDER,SYCART,FNRMDR,WORK(KSYMTD),NDIMF,NCOOR,LTXT,
     &               IPRINT,TEXT)
C
C        *** Transforming to normal coordinates. ***
C
         LTXT = 6
         TEXT(1:6) = 'normal'
C
         KTMPFD = 1
         CALL DCOPY(NCOOR**4,FNRMDR,1,WORK(KTMPFD),1)
         CALL TRAFDR(CRTNRM,WORK(KTMPFD),FNRMDR,NCOOR,NDCOOR,NCOOR,
     &               TEXT,LTXT,IPRINT)
C
C        *** Printing. ***
C
         IF (.NOT.MINOUT)
     &      CALL PRFDER(FNRMDR,NCOOR,NDCOOR,TEXT,LTXT,IPRINT)
      END IF
C
      RETURN
      END
C
C
C     /* Deck prderv*/
      SUBROUTINE PRTDER(TDER,NDIM,NCOR,TEXT,LTXT,IPRINT)
C     **********************************************************
C     *** Printing of third derivatives in TEXT coordinates. ***
C     **********************************************************
#include <implicit.h>
#include <priunit.h>
      CHARACTER*80 TEXT
      DIMENSION TDER(NDIM,NDIM,NDIM)

C
C     *** Header print. ***
C
      CALL HEADER('Third derivative of energy in ' // TEXT(1:LTXT)
     &            // ' coordinates',-1)
C
C     *** Printing of force field. ***
C
      IF (MOD(NCOR,6).EQ.0) THEN
         NLCMAX = NCOR/6
      ELSE
         NLCMAX = INT(NCOR/6)+1
      END IF
C
      DO 100 ICOL2 = 1, NCOR
         WRITE (LUPRI,'(A,I5)') '      Column number', ICOL2
         WRITE (LUPRI,'(A)')    '      ------------------'
         INLC = 0
         DO 200 INLCMX = 1, NLCMAX
            INLC2 = 6*(INLCMX-1) + 1
            INLC  = MIN(INLC+6,NCOR)
            DO 300 ICOL1 = 1, NCOR
               WRITE (LUPRI,'(A,6F10.6)') '   ',
     &                           (TDER(I,ICOL1,ICOL2), I=INLC2, INLC)
 300        CONTINUE
            WRITE (LUPRI,'(A)') '                              '
 200     CONTINUE
 100  CONTINUE
C
      RETURN
      END
C
C
C     /* Deck prfder*/
      SUBROUTINE PRFDER(FDER,NDIM,NCOR,TEXT,LTXT,IPRINT)
C     **********************************************************
C     *** Printing of fourth derivative in TEXT coordinates. ***
C     **********************************************************
#include <implicit.h>
#include <priunit.h>
      CHARACTER*80 TEXT
      DIMENSION FDER(NDIM,NDIM,NDIM,NDIM)
C
C     *** Header print. ***
C
      CALL HEADER('Fourth derivative of energy in ' // TEXT(1:LTXT)
     &            // ' coordinates',-1)
C
C     *** Printing of derivative. ***
C
      IF (MOD(NCOR,6).EQ.0) THEN
         NLCMAX = NCOR/6
      ELSE
         NLCMAX = INT(NCOR/6)+1
      END IF
C
      DO 100 ICOL3 = 1, NCOR
         WRITE (LUPRI,'(A,I4)') '      The fourth dimension', ICOL3
         WRITE (LUPRI,'(A/)') '      ------------------------'
         DO 200 ICOL2 = 1, NCOR
            WRITE (LUPRI,'(A,I4)') '        The third Dimension', ICOL2
            INLC = 0
            DO 300 INLCMX = 1, NLCMAX
               INLC2 = 6*(INLCMX-1) + 1
               INLC  = MIN(INLC+6,NCOR)
               DO 400 ICOL1 = 1, NCOR
                  WRITE (LUPRI,'(A,6F10.6)') '   ',
     &                 (FDER(I,ICOL1,ICOL2,ICOL3), I=INLC2, INLC)
 400           CONTINUE
               WRITE (LUPRI,'(A)') '                              '
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
C
      RETURN
      END
C
C
C     /*Deck prderv*/
      SUBROUTINE PRDERV(TDER,FDER,TSTGDR,TSTSDR,SYMCOR,CSTART,TTMPDR,
     &                  FTMPDR,RNNORM,WORK,ICRIRP,LWORK,NPRRDR,NDIMT,
     &                  NDIMF,LTXT,IPRINT,TEXT)
C     ***************************************************************
C     **** This routine prints out the derivatives of the energy ****
C     **** to NPRRDR order. These are done in 'TEXT' coordinates ****
C     ***************************************************************
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <dummy.h>
C
#include <symmet.h>
#include <nuclei.h>
#include <trkoor.h>
#include <taymol.h>
#include <cbiwlk.h>
#include <cbinum.h>
#include <numder.h>
#include <pvibav.h>
      LOGICAL CPRPBK, PRWHLE
      CHARACTER*80 TEXT
      DIMENSION TDER(NDIMT), FDER(NDIMF), TSTGDR(NCOOR), CSTART(NCOOR),
     &          TSTSDR(NCOOR,NCOOR), SYMCOR(NCOOR,NCOOR),
     &          TTMPDR(NCOOR,NCOOR,NCOOR),
     &          RNNORM(NCOOR),
     &          FTMPDR(NCOOR,NCOOR,NCOOR,NCOOR), WORK(LWORK)
      DIMENSION ICRIRP(NCOOR,2)
C
C     *** We are finished calculating properties. But are ***
C     *** backing it up for later use.                    ***
C
      CPRPBK = CNMPRP
      CNMPRP = .FALSE.
C
C     *** Print gradient ***
C
      IF (NPRRDR.GT.0) THEN
C
         CALL HEADER('Numerical gradient in ' // TEXT(1:LTXT) //
     &               ' coordinates',-1)
C
         IF (NAORDR.GE.1) THEN
            KCSTRA = 1
            KSCTRA = KCSTRA + NCOOR**2
            KEGRAD = KSCTRA + MXCOOR
            KSEGRD = KEGRAD + MXCOOR
            KLAST  = KSEGRD + NCOOR
            LWRK1  = LWORK  - KLAST + 1
            IF (KLAST.GT.LWORK) CALL QUIT('Memory exceeded in TRFCGD')
            CALL TRAGRD(GRDMOL,WORK(KEGRAD),WORK(KCSTRA),WORK(KSCTRA),
     &                  NCRREP(0,1),NCOOR)
            CALL TRFCGD(WORK(KEGRAD),SYMCOR,CSTART,WORK(KSEGRD),
     &                  WORK(KLAST),NCOOR,NDCOOR,LWRK1,IPRINT)
            CALL DCOPY(MXCOOR,WORK(KEGRAD),1,GRDMOL,1)
         END IF
C
         KCGRAD = 1
         KWRK   = KCGRAD + MXCOOR
         LWRK   = LWORK - KWRK
         LNEED  = 2*MXCOOR*MXCOOR
         IF ((LWRK - LNEED) .LT. 0)
     &      CALL STOPIT('PRDERV','GSPGRD',LWRK,LNEED)
         CALL GSPGRD(SYMCOR,WORK(KCGRAD),WORK(KWRK),LWRK,ICRIRP,LTXT,
     &               IPRINT,TEXT)
C
C        *** If comparing with the analytical gradient. ***
C
         IF (SDRTST) THEN
            DO 200 IC1 = 1, NCOOR
               TSTGDR(IC1) = WORK(KCGRAD-1+IC1)
 200        CONTINUE
         END IF
      END IF
C
C     *** Print hessian ***
C
      IF (NPRRDR.GT.1) THEN
C
         CALL HEADER('Numerical hessian in ' // TEXT(1:LTXT) //
     &               ' coordinates',-1)
C
         IF ((NAORDR.GE.2).AND..NOT.PREHES) THEN
            KCSTRA = 1
            KSCTRA = KCSTRA + NCOOR**2
            KEHESS = KSCTRA + NCOOR**2
            KSEHSS = KEHESS + MXCOOR**2
            KLAST  = KSEHSS + NCOOR**2
            LWRK1  = LWORK  - KLAST + 1
            IF (KLAST.GT.LWORK) CALL QUIT('Memory exceeded in TRFCHS')
            CALL TRAHES(HESMOL,MXCOOR,WORK(KEHESS),WORK(KCSTRA),
     &                  WORK(KSCTRA),MXCOOR,NCOOR,1)
            CALL TRFCHS(WORK(KEHESS),SYMCOR,CSTART,WORK(KSEHSS),
     &                  WORK(KLAST),NCOOR,NDCOOR,LWRK1,IPRINT)

            CALL DCOPY(MXCOOR**2,WORK(KEHESS),1,HESMOL,1)
         END IF
C
         KSMCIN = 1
         KCHES1 = KSMCIN + NCOOR**2
         KCHES2 = KCHES1 + NCOOR**2
         KLAST  = KCHES2 + NCOOR**2
         LWRK   = LWORK  - KLAST
         CALL GSPHES(SYMCOR,WORK(KSMCIN),WORK(KCHES1),WORK(KCHES2),
     &               WORK(KLAST),ICRIRP,NDCOOR,LWRK,LTXT,IPRINT,TEXT)
C
C        *** If comparing with the analytical hessian. ***
C
         IF (SDRTST) THEN
            IC12 = 0
            DO 300 IC2 = 1, NCOOR
            DO 300 IC1 = 1, NCOOR
               IC12 = IC12 + 1
               TSTSDR(IC1,IC2) = WORK(KCHES2-1+IC12)
 300        CONTINUE
         END IF
      END IF
C
C     *** Print third derivative of energy ***
C
      IF (NPRRDR.GT.2) THEN
C
         PRWHLE = .NOT.(ANALZ1.AND.NRMCRD.AND.((NMORDR+NAORDR).EQ.3))
C
         KTDER  = 1
         KSYMTD = KTDER  + NCOOR**3
         KLAST  = KSYMTD + NCOOR**3
         CALL HEADER('Numerical third derivative of energy in ' //
     &                TEXT(1:LTXT) // ' coordinates',-1)
         CALL PRITDR(TDER,SYMCOR,TTMPDR,WORK(KSYMTD),NDIMT,NDCOOR,LTXT,
     &               IPRINT,PRWHLE,TEXT)
      END IF
C
C     *** Print fourth derivative of energy ***
C
      IF (NPRRDR.GT.3) THEN
         KFDER  = 1
         KSYMFD = KFDER  + NCOOR**4
         KSCTRA = KSYMFD + NCOOR**4
         KCSTRA = KSCTRA + NCOOR**2
         KLAST  = KCSTRA + NCOOR**2
         CALL HEADER('Numerical fourth derivative of energy in ' //
     &                TEXT(1:LTXT) // ' coordinates',-1)
         CALL PRIFDR(FDER,SYMCOR,FTMPDR,WORK(KSYMFD),NDIMF,NDCOOR,LTXT,
     &               IPRINT,TEXT)
      END IF
C
C     *** Writing to spectro file if requested. ***
C
      IF (SPECTR) THEN
         NTIME = 2
         CALL WRISPC(DUMMY,DUMMY,TTMPDR,FTMPDR,TEXT(1:6),NCOOR,NDCOOR,
     &               NTIME,IPRINT)
      END IF
      IF (MIDAS) THEN
         NTIME = 2
         CALL WRIMOP(DUMMY,RNNORM,TTMPDR,FTMPDR,TEXT(1:6),NCOOR,NDCOOR,
     &               NTIME,IPRINT)
      END IF
C
      IF ((IPRINT .GT. 20).AND.(NPRRDR.EQ.2)) THEN
         CALL HEADER ('Copy of Hessian for test',-1)
         DO IC1 = 1, NCOOR
            WRITE (LUPRI,'(12F12.8)') (TSTSDR(IC1,IC2),IC2=1,NCOOR)
         END DO
      END IF
C
C     *** Restoring CNMPRP. ***
C
      CNMPRP = CPRPBK
C
      RETURN
      END
C
C
C
C     /*Deck gspgrd*/
      SUBROUTINE GSPGRD(SYMCOR,CGRAD,WORK,LWORK,ICRIRP,LTXT,IPRINT,TEXT)
#include <implicit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <priunit.h>
C
#include <symmet.h>
#include <nuclei.h>
#include <trkoor.h>
#include <taymol.h>
#include <cbiwlk.h>
#include <numder.h>
      CHARACTER TEXT*80
      DIMENSION SYMCOR(NCOOR,NCOOR), CGRAD(NCOOR), ICRIRP(NCOOR,2),
     &          WORK(LWORK)
C
      CALL HEADER('Gradient in ' // TEXT(1:LTXT) // ' coordinates',-1)
C
      DO 100 ICOOR = 1, NDCOOR
         IF (ICRIRP(ICOOR,1).EQ.1) THEN
            WRITE (LUPRI,'(A,F16.8)') '                               ',
     &                                GRDMOL(ICOOR)
         END IF
 100  CONTINUE
C
      CALL DZERO(CGRAD,NCOOR)
      DO 200 ICOOR2 = 1, NDCOOR
      DO 200 ICOOR1 = 1, NCOOR
         CGRAD(ICOOR1) = CGRAD(ICOOR1)
     &                 + SYMCOR(ICOOR1,ICOOR2)*GRDMOL(ICOOR2)
 200  CONTINUE
C
      CALL HEADER('Gradient in cartesian coordinates',-1)
C
      IOFF = 0
      DO 300 ICENT = 1, NUCDEP
         WRITE (LUPRI,'(1X,A6,F17.10,2F24.10)') NAMDEP(ICENT),
     &                                         (CGRAD(IOFF+J), J=1,3)
         IOFF = IOFF + 3
 300  CONTINUE
C
C     *** Transform to symmetry basis used in Dalton in ***
C     *** case of geometry optimization                 ***
C
      IF (MAXREP .GT. 0) THEN
         KCSTRA = 1
         KSCTRA = KCSTRA + MXCOOR*MXCOOR
         CALL TRACOR(WORK(KCSTRA),WORK(KSCTRA),1,MXCOOR,IPRINT)
         CALL TRACTS(CGRAD,3*NUCDEP,WORK(KCSTRA))
         CALL DCOPY(3*NUCDEP,CGRAD,1,GRDMOL,1)
C
         CALL HEADER('Gradient in Dalton symmetry coordinates',-1)
C
         DO 202 I = 1, NCRREP(0,1)
            WRITE (LUPRI,'(25X,A6,F17.10)') NAMEX(IPTCOR(I,1)),GRDMOL(I)
 202     CONTINUE
      END IF
C
C     *** If testing ***
C
      IF (SDRTST) THEN
         DO IC1 = 1, NCOOR
            GRDMOL(IC1) = CGRAD(IC1)
         END DO
      END IF
C
C     *** Print ***
C
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Symcor matrix in GSPGRD',-1)
         DO 400 I = 1, NCOOR
            WRITE (LUPRI,'(24F12.7)') (SYMCOR(I,J),J=1,NCOOR)
 400     CONTINUE
      END IF
C
      RETURN
      END
C
C
C     /*Deck gsphes*/
      SUBROUTINE GSPHES(SYMCOR,SMCINV,CHESS1,CHESS2,WORK,ICRIRP,
     &                  NDCOOR,LWORK,LTXT,IPRINT,TEXT)
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
      PARAMETER (KCOL=6)
#include <nuclei.h>
#include <trkoor.h>
#include <taymol.h>
#include <cbiwlk.h>
#include <cbinum.h>
      INTEGER BEGIN, LAST
      LOGICAL HESEXS
      CHARACTER TEXT*80
      DIMENSION SYMCOR(NCOOR,NCOOR), SMCINV(NCOOR,NCOOR),
     &          CHESS1(NCOOR,NCOOR), CHESS2(NCOOR,NCOOR),
     &          WORK(LWORK), ICRIRP(NCOOR,    2)
C
C     *** Print hessian in symmetric coordinates ***
C
      CALL HEADER('Hessian in ' // TEXT(1:LTXT) // ' coordinates',-1)
C
      BEGIN = 1
      LAST  = MIN(NDCOOR,KCOL)
      KCOOR = NDCOOR
      NCOL  = NDCOOR/KCOL
      IF (MOD(NDCOOR,KCOL).NE.0) NCOL = NCOL + 1
C
      DO 100 ICOL = 1, NCOL
         WRITE (LUPRI,1000) (ICRIRP(I,1),I = BEGIN,LAST)
C
         DO 200 ICOOR = BEGIN, NDCOOR
            WRITE (LUPRI,2000) ICRIRP(ICOOR,1),
     &                        (HESMOL(ICOOR,I),I=BEGIN,MIN(LAST,ICOOR))
 200     CONTINUE
         WRITE (LUPRI,'(A)') '                                       '
         BEGIN = BEGIN + KCOL
         LAST  = MIN(NDCOOR,KCOL+LAST)
 100  CONTINUE
C
      DO 300 J = 1, NDCOOR
      DO 300 I = 1, J
         HESMOL(I,J) = HESMOL(J,I)
 300  CONTINUE
C
C     *** Transpose of coordinate transformation matrix ***
C
      IF (TEXT(1:6) .NE. 'normal') THEN
         DO 350 J = 1, NDCOOR
         DO 350 I = 1, NCOOR
            SMCINV(J,I) = SYMCOR(I,J)
 350     CONTINUE
C
C     *** Transform to cartesian hessian ***
C
         KDIM = NCOOR**2
         CALL DZERO(CHESS1,KDIM)
         DO 400 K = 1, NCOOR
         DO 400 J = 1, NDCOOR
         DO 400 I = 1, NDCOOR
            CHESS1(I,K) = CHESS1(I,K) + HESMOL(I,J)*SMCINV(J,K)
 400     CONTINUE
C
         KDIM = NCOOR**2
         CALL DZERO(CHESS2,KDIM)
         DO 500 K = 1, NCOOR
         DO 500 J = 1, NDCOOR
         DO 500 I = 1, NCOOR
            CHESS2(I,K) = CHESS2(I,K) + SYMCOR(I,J)*CHESS1(J,K)
 500     CONTINUE
C
C        *** Print cartesian hessian ***
C
         CALL HEADER('Cartesian Hessian in GSPHES',-1)
         CALL PR2DER(CHESS2,NCOOR,NCOOR,LUPRI)
C
C        *** Print to file if we are going to reuse the Hessian. ***
C
         IF (REUHES) THEN
            INQUIRE(FILE='DALTON.HES',EXIST=HESEXS)
C
C           *** No hessian specified, we can safely write to file. ***
            IF (.NOT. HESEXS) THEN
C
C              *** Open hessian file. ***
               LUHES = -1
               CALL GPOPEN(LUHES,'DALTON.HES','NEW',' ','FORMATTED',
     &                     IDUMMY,.FALSE.)
C
C              *** Checking if this is going to be used with SPECTRO.***
C
               IF (SPECTR) THEN
                  NTIMES = NCOOR/3
                  DO ICOOR2 = 1, NCOOR
                     DO ITIMES = 1, NTIMES
                        ISTART = 3*(ITIMES-1) + 1
                        WRITE (LUHES,'(3F22.12)')
     &                    (CHESS2(ICOOR1,ICOOR2),ICOOR1=ISTART,ISTART+2)
                     END DO
                  END DO
               ELSE
C
C              *** Printing necessary pre-hessian information. ***
                  WRITE(LUHES,'(A)') 'CARTESIAN HESSIAN'
                  WRITE(LUHES,*) NCOOR
                  WRITE(LUHES,'(A)') '                                '
C
                  DO 800 ICOOR2 = 1, NCOOR
                     DO 900 ICOOR1 = 1, NCOOR
                        WRITE (LUHES,'(F22.12)') CHESS2(ICOOR1,ICOOR2)
 900                 CONTINUE
                     WRITE (LUHES,'(A)') '                            '
 800              CONTINUE
               END IF
            ELSE
               WRITE (LUPRI,'(//A/A//)')
     &            'Hessian file "DALTON.HES" already exists.' //
     &              ' This file will NOT be overwritten.',
     &            'Please restart the calculation without this file.'
               CALL QUIT('"DALTON.HES" already exists. See output.')
            END IF
         END IF
C
C        *** Print ***
C
         IF (IPRINT .GT. 20) THEN
C
            KDIM = NCOOR**2
            CALL DZERO(CHESS1,KDIM)
C
            CALL HEADER('Symcor matrix',-1)
            DO 1100 I = 1, NCOOR
               WRITE (LUPRI,'(24F9.6)') (SYMCOR(I,J),J=1,NDCOOR)
 1100       CONTINUE
            WRITE (LUPRI,'(A)') '                                  '
C
            CALL HEADER('Inverse of symcor matrix',-1)
            DO 1200 I = 1, NDCOOR
               WRITE (LUPRI,'(24F9.6)') (SMCINV(I,J),J=1,NCOOR)
 1200       CONTINUE
            WRITE (LUPRI,'(A)') '                                  '
C
            DO 1300 K = 1, NCOOR
            DO 1300 J = 1, NDCOOR
            DO 1300 I = 1, NCOOR
               CHESS1(I,K) = CHESS1(I,K) + SYMCOR(I,J)*SMCINV(J,K)
 1300       CONTINUE
C
            CALL HEADER('Should be unit matrix',-1)
            DO 1400 J = 1, NCOOR
               WRITE (LUPRI,'(24F9.6)') (CHESS1(I,J),I=1,NCOOR)
 1400       CONTINUE
         END IF
      END IF
C
 1000 FORMAT (8X,6(3X,I4,5X),(3X,I4,5X))
 2000 FORMAT (1X,I4,2X,6F12.6)
      RETURN
      END
C
C
C     /* Deck tsths1*/
      SUBROUTINE TSTHS1(SYMCOR,HESMWT,EGNVCT,TM1TMP,TM2TMP,DKIN,WORK,
     &                  ICRIRP,LWORK)
************************************************************
*** Tests the molecular hessian by makeing the cartesian ***
*** hessian, diagonalizing the mass-weighted hessian and ***
*** prints the eigenvalues, and harmonic frequancies     ***
************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <codata.h>
      PARAMETER (DMTHR = 1.0D-9, D1 = 1.0D0)
#include <taymol.h>
#include <trkoor.h>
#include <nuclei.h>
#include <cbinum.h>
      DIMENSION SYMCOR(NCOOR,NCOOR), EGNVCT(NCOOR,NCOOR),
     &          HESMWT(NCOOR*(NCOOR+1)/2), DKIN(NCOOR,NCOOR),
     &          TM1TMP(NCOOR,NCOOR),TM2TMP(NCOOR,NCOOR),
     &          WORK(LWORK)
C
      KDIM = NCOOR**2
      CALL DZERO(TM1TMP,KDIM)
      DO 100 IC3 = 1, NCOOR
      DO 100 IC2 = 1, NCOOR
      DO 100 IC1 = 1, NCOOR
         TM1TMP(IC1,IC3) = TM1TMP(IC1,IC3)
     &                   + HESMOL(IC1,IC2)*SYMCOR(IC3,IC2)
 100  CONTINUE
C
      CALL DZERO(TM2TMP,KDIM)
      DO 200 IC3 = 1, NCOOR
      DO 200 IC2 = 1, NCOOR
      DO 200 IC1 = 1, NCOOR
         TM2TMP(IC1,IC3) = TM2TMP(IC1,IC3)
     &                   + SYMCOR(IC1,IC2)*TM1TMP(IC2,IC3)
 200  CONTINUE
C
      CALL HEADER('Cartesian hessian',-1)
      DO 300 IC1 = 1, NCOOR
         WRITE (LUPRI,'(12F8.5)') (HESMOL(IC1,IC2),IC2=1,NCOOR)
 300  CONTINUE
C
      CALL DZERO(TM1TMP,KDIM)
      DO 400 IC3 = 1, NCOOR
      DO 400 IC2 = 1, NCOOR
      DO 400 IC1 = 1, NCOOR
         TM1TMP(IC1,IC3) = TM1TMP(IC1,IC3)
     &                   + TM2TMP(IC1,IC2)*DKIN(IC2,IC3)
 400  CONTINUE
C
      CALL DZERO(TM2TMP,KDIM)
      DO 500 IC3 = 1, NCOOR
      DO 500 IC2 = 1, NCOOR
      DO 500 IC1 = 1, NCOOR
         TM2TMP(IC1,IC3) = TM2TMP(IC1,IC3)
     &                   + DKIN(IC1,IC2)*TM1TMP(IC2,IC3)
 500  CONTINUE
C
      IC12 = 0
      DO 600 IC2 = 1, NCOOR
      DO 600 IC1 = 1, IC2
         IC12 = IC12 + 1
         HESMWT(IC12) = TM2TMP(IC1,IC2)
 600  CONTINUE
C
      CALL HEADER('Mass-weighted hessian',-1)
      DO 700 IC1 = 1, NCOOR
         WRITE (LUPRI,'(12F8.5)') (TM2TMP(IC1,IC2),IC2=1,IC1)
 700  CONTINUE
C
      KWRK  = 1
      KIWRK = KWRK + NCOOR
      CALL DUNIT(EGNVCT,NCOOR)
      CALL JACO(HESMWT,EGNVCT,NCOOR,NCOOR,NCOOR,WORK(KWRK),WORK(KIWRK))
C
      CALL HEADER('Diagonalized hessian',-1)
      DO 800 IC1 = 1, NCOOR
         ISTART = (IC1*(IC1-1))/2 + 1
         IEND   = (IC1*(IC1+1))/2
         WRITE (LUPRI,'(12F8.5)') (HESMWT(IC12),IC22=ISTART,IEND)
 800  CONTINUE
C
      RETURN
      END
C
C
C     /*Deck sdertt*/
      SUBROUTINE SDERTT(TSTSDR,TSTGDR,SYMCOR,TMPGRD,TMPHES,BKPGRD,
     &                  BKPHES,WORK,LWORK,WRKDLM,IPRINT)
C     *******************************************************************
C     *** This routine tests the numerical derivatives with available ***
C     ***                   analytical derivatives.                   ***
C     *** NOTE: TMPGRD has dimension MXCOOR due to old code.          ***
C     *******************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
      PARAMETER (D0 = 0.0D0)
#include <numder.h>
#include <trkoor.h>
#include <taymol.h>
#include <symmet.h>
#include <abainf.h>
#include <exeinf.h>
#include <gnrinf.h>
#include <past.h>
#include <inftap.h>
      DIMENSION TSTSDR(NCOOR,NCOOR), TSTGDR(NCOOR ),
     &          TMPHES(NCOOR,NCOOR), TMPGRD(MXCOOR),
     &          BKPHES(NCOOR,NCOOR), BKPGRD(NCOOR ),
     &          SYMCOR(NCOOR,NCOOR), WORK  (LWORK )
C
      DO 100 IC2 = 1, NCOOR
         BKPGRD(IC2) = GRDMOL(IC2)
         DO 200 IC1 = 1, NCOOR
            BKPHES(IC1,IC2) = HESMOL(IC1,IC2)
 200     CONTINUE
 100  CONTINUE
C
      DO 300 IC2 = 1, NCOOR
      DO 300 IC1 = 1, NCOOR
         TSTGDR(IC1) = TSTGDR(IC1) + SYMCOR(IC1,IC2)*GRDMOL(IC2)
 300  CONTINUE
C
      CALL DGEMM('N','N',NCOOR,NCOOR,NCOOR,1.D0,
     &           SYMCOR,NCOOR,
     &           HESMOL,MXCOOR,0.D0,
     &           TMPHES,NCOOR)
C
      CALL DGEMM('N','T',NCOOR,NCOOR,NCOOR,1.D0,
     &           TMPHES,NCOOR,
     &           SYMCOR,NCOOR,0.D0,
     &           TSTSDR,NCOOR)
C
      KDIM = MXCOOR**2
      CALL DZERO(GRDMOL,MXCOOR)
      CALL DZERO(HESMOL,KDIM  )
C
      MOLGRD = .TRUE.
      MOLHES = .TRUE.
      PASEXC = .FALSE.
      RNABAC = .TRUE.
      WRINDX = .TRUE.
      FTRONV = .TRUE.
      DOWALK = .FALSE.
      LUSUPM = -1
      WORK(1) = WRKDLM
      CALL ABAINP('**PROPE',WORK(2),LWORK)
      CALL EXEABA(WORK(1),LWORK-1,WRKDLM)
C
C     *** Transforming the analytical derivatives into ***
C     ***               cartesian basis.               ***
C
      IF (MAXREP.GT.0) THEN
         KSCTR = 1
         KCSTR = KSCTR + NCOOR**2
         KLAST = KCSTR + NCOOR**2
         IF (KLAST.GT.LWORK) CALL QUIT('Memory exceeded inside SDERTT')
         CALL TRAGRD(GRDMOL,TMPGRD,WORK(KCSTR),WORK(KSCTR),NCRREP(0,1),
     &               NCOOR)
         CALL TRAHES(HESMOL,MXCOOR,TMPHES,WORK(KCSTR),WORK(KSCTR),NCOOR,
     &               NCOOR,1)
      ELSE
         DO 400 IC2 = 1, NCOOR
            TMPGRD(IC2) = GRDMOL(IC2)
            DO 500 IC1 = 1, NCOOR
               TMPHES(IC1,IC2) = HESMOL(IC1,IC2)
 500        CONTINUE
 400     CONTINUE
      END IF
C
      RMAXGD = D0
      DO 600 J = 1, NCOOR
         RGRDJ = ABS(TMPGRD(J)-TSTGDR(J))
         IF (RGRDJ .GT. RMAXGD) THEN
            RMAXGD = RGRDJ
            NMG    = J
         END IF
 600  CONTINUE
      CALL HEADER('Comparison of numerical and analytical gradients',-1)
      WRITE (LUPRI,'(//A,E13.7,A,I4)')  'Largest difference ', RMAXGD,
     &                                  ' for element:', NMG
      WRITE (LUPRI,'(A,2E13.7//)') 'The values of these elements are: ',
     &                             TMPGRD(NMG),TSTGDR(NMG)

      RLRGST = D0
      DO 700 J = 1, NCOOR
      DO 700 I = 1, J
         RINTMD = (TMPHES(I,J)-TSTSDR(I,J))**2
         IF ( RINTMD .GT. RLRGST) THEN
            RLRGST = RINTMD
            NMI    = I
            NMJ    = J
            HVALC  = TMPHES(I,J)
            HVALN  = TSTSDR(I,J)
         END IF
 700  CONTINUE
C
      CALL HEADER('Comparison of numerical and analytical Hessians',-1)
      WRITE (LUPRI,'(//A,E13.7,A,2I4)') 'Largest difference ',
     &                 SQRT(RLRGST), ' for elements:', NMI, NMJ
      WRITE (LUPRI,'(A,2E13.7//)') 'The values of these elements are: ',
     &                              HVALC, HVALN
C
      DO 800 IC2 = 1, NCOOR
         GRDMOL(IC2) = BKPGRD(IC2)
         DO 900 IC1 = 1, NCOOR
            HESMOL(IC1,IC2) = BKPHES(IC1,IC2)
 900     CONTINUE
 800  CONTINUE
C
      RETURN
      END
C
C
C     /* Deck drnrmc*/
      SUBROUTINE DRNRMC(SYMCOR,ICRIRP,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <trkoor.h>
#include <numder.h>
#include <fcsym.h>
      LOGICAL FOUND
      DIMENSION SYMCOR(NCOOR,NCOOR), ICRIRP(NCOOR,2)
C
      DO 100 IC1 = 1, NMREDU
         FOUND = .FALSE.
         IIREP = KDRYRN(IC1)
C
         DO 200 IC2 = 1, NDCOOR
            IF ((ICRIRP(IC2,1).EQ.IIREP).AND.(.NOT.FOUND)) THEN
               FOUND = .TRUE.
               DO 300 II  = 1,     2
               DO 300 IC3 = IC2+1, NDCOOR
                  ICRIRP(IC3-1,II) = ICRIRP(IC3,II)
 300           CONTINUE
               NDCOOR = NDCOOR - 1
C
               IF (IIREP.GT.N1DIME) THEN
                  FOUND = .FALSE.
                  DO 400 IC3 = IC2-1, NDCOOR
                     IF ((ICRIRP(IC3,1).EQ.IIREP).AND.
     &                   (ICRIRP(IC3,2).EQ.    1).AND.(.NOT.FOUND)) THEN
                        FOUND = .TRUE.
                        DO 500 II  =     1      2
                        DO 500 IC4 = IC3+1, NDCOOR
                           ICRIRP(IC4-1,II) = ICRIRP(IC4,II)
 500                    CONTINUE
                     END IF
 400              CONTINUE
                  NDCOOR = NDCOOR - 1
               END IF
            END IF
 200     CONTINUE
 100  CONTINUE
C
      IF (IPRINT.GT.20) THEN
         WRITE (LUPRI,'(5X,A)') 'Removed translational and rotational'//
     &         'redundencies.'
         WRITE (LUPRI,'(A)') '                                    '
         WRITE (LUPRI,'(5X,A)') 'Symmetry of coordinates left:'
         WRITE (LUPRI,'(5X,24I4)') (ICRIRP(II,1),II=1,NDCOOR)
         WRITE (LUPRI,'(5X,24I4)') (ICRIRP(II,2),II=1,NDCOOR)
      END IF
      RETURN
      END
C
C
C     /* Deck wricor*/
      SUBROUTINE WRICOR(SYMCOR,RNNORM,FREQ,ICRIRP,LURSTR,IPRINT)
C     *************************************************
C     *** This is a routine that writes out normal  ***
C     *** coordinates to file, in case of a restart.***
C     *************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <trkoor.h>
#include <numder.h>
      DIMENSION SYMCOR(NCOOR,NCOOR), RNNORM(NCOOR), FREQ(NCOOR)
      DIMENSION ICRIRP(NCOOR,2)
C
      WRITE (LURSTR,'(I8)') NDCOOR
C
C     *** Writing normal coordinates. ***
C
      DO 100 IC2 = 1, NDCOOR
         WRITE(LURSTR,'(2I5)') (ICRIRP(IC2,I),I=1,2)
         DO 200 IC1 = 1, NCOOR
            WRITE(LURSTR,'(F24.16)') SYMCOR(IC1,IC2)
 200     CONTINUE
 100  CONTINUE
C
C     *** Norm of the non-normalized normal coordinates. ***
C
      WRITE(LURSTR,'(A)') 'Norm'
      DO 300 IC = 1, NDCOOR
         WRITE(LURSTR,'(F24.16)') RNNORM(IC)
 300  CONTINUE
C
C     *** Frequencies. ***
C
      WRITE(LURSTR,'(A)') 'Freq'
      DO 400 IC = 1, NDCOOR
         WRITE(LURSTR,'(F24.16)') FREQ(IC)
 400  CONTINUE
      CALL FLSHFO(LURSTR)
C
      RETURN
      END
C
C
C     /*Deck rerstr*/
      SUBROUTINE RERSTR(FUVAL,SYMCOR,RNNORM,FREQ,ICRIRP,NDIME,NINTIN,
     &                  KEND,IDIMAX,IDIMIN,LURSTR,IPRINT,RSTDON)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <trkoor.h>
#include <numder.h>
      LOGICAL RSTDON
      DIMENSION SYMCOR(NCOOR,NCOOR), FUVAL(NINTIN,NDIME), RNNORM(NCOOR),
     &          FREQ(NCOOR)
      DIMENSION ICRIRP(NCOOR,2)
C
      CALL DZERO(FUVAL,NINTIN*NDIME)
C
      READ(LURSTR,*) II
C
C     *** The program ended the second time around. ***
C     *** Need to read in additional information.   ***
C
      IF ((II.EQ.1) .AND. (II.EQ.KEND)) THEN
         CALL RDHDRS(SYMCOR,RNNORM,FREQ,ICRIRP,NCOOR,NDCOOR,LURSTR)
      END IF

      IF (((II.EQ.0) .AND. (II.EQ.KEND)).OR.
     &    ((II.EQ.1) .AND. (II.EQ.KEND))) THEN
         RSTDON = .TRUE.
C
C        *** Reading the function values ***
C
 100     CONTINUE
         READ(LURSTR,FMT=*,IOSTAT=IOS) IINTIN, IDIME, ENERGY
         IF (IOS.GE.0) THEN
            IDIMAX = MAX(IDIMAX,IDIME)
            IDIMIN = MIN(IDIMIN,IDIME)
            FUVAL(IINTIN,IDIME) = ENERGY
            GOTO 100
         END IF
      END IF
C
      RETURN
      END
C
C
C     /*Deck prprer*/
      SUBROUTINE PRPRER(WORK,IDIMAX,IDIMIN,LURSTR,LWORK)
C     *************************************************
C     *** Restart routine for property derivatives. ***
C     *************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      CHARACTER*9 PRPTXT
      DIMENSION WORK(LWORK)
#include <numder.h>
#include <trkoor.h>
#include <dummy.h>
C
      KNMPRP = 0
      INMTCL = 0
 100  CONTINUE
      READ (LUNDPR,FMT='(A9)',IOSTAT=IOS) PRPTXT
      IF (IOS.GE.0) THEN
         READ (LUNDPR,'(4I7)') NDIM1, NDIM2, NDIM3, KNMCLC
C
C        *** Figures out which property to read in, and reads it. ***
C
         KGRBG = 1
         CALL CHPRRD(WORK(KGRBG),WORK(KGRBG),WORK(KGRBG),WORK(KGRBG),
     &               WORK(KGRBG),NDIM1,NDIM2,NDIM3,PRPTXT,IDUMMY)
C
C        *** Another property has been calculated. ***
C
         INMTCL = INMTCL + 1
C
C        *** Calculating number of properties per geometry. ***
C
         IF (KNMCLC .EQ. 1) THEN
            KNMPRP = KNMPRP + 1
         END IF
C
         GOTO 100
      END IF
C
      IF (KNMCLC.EQ.1) THEN
C
C        *** Only some properties for first geometry has been      ***
C        *** calculated. Nothing to save, continue from beginning. ***
C
         REWIND(LUNDPR)
         IDIMAX = IDIMAX - 1
      ELSE
C
C        *** If not all properties were written for the ***
C        *** last geometry. We need to make sure that we***
C        *** are at the end of a geometry.              ***
C
         IF (KNMCLC*KNMPRP.GT.INMTCL) THEN
C
C           *** We can only use the restart for the previous ***
C           *** geometry.                                    ***
C
            IF (KNMCLC.EQ.IDIMAX) THEN
C              *** Original geometry. ***
               IDIMIN = 2
            ELSE
               IDIMAX = IDIMAX -1
            END IF
            KNMCLC = KNMCLC - 1
C
C           *** Positioning the property file. ***
C
            REWIND(LUNDPR)
            DO INMCLC = 1, KNMCLC
            DO INMPRP = 1, KNMPRP
               READ (LUNDPR,FMT='(A9)',IOSTAT=IOS) PRPTXT
               READ (LUNDPR,'(4I7)') NDIM1, NDIM2, NDIM3, KGRB
C
C              *** Figures out which property to read in, and ***
C              *** reads it.                                  ***
C
               KGRB = 1
               CALL CHPRRD(WORK(KGRB),WORK(KGRB),WORK(KGRB),WORK(KGRB),
     &                     WORK(KGRB),NDIM1,NDIM2,NDIM3,PRPTXT,IDUMMY)
            END DO
            END DO
C
C           *** Positioning the RSTRT.FC file. ***
C
            REWIND(LURSTR)
            READ(LURSTR,*) II
            IF (II.EQ.1) THEN
               KGRB = 1
               CALL RDHDRS(WORK(KGRB),WORK(KGRB),WORK(KGRB),
     &                     WORK(KGRB),NCOOR,IDUMMY,LURSTR)
            END IF
C
            DO ID = IDIMIN, IDIMAX
               READ(LURSTR,FMT=*,IOSTAT=IOS) KGRBG1, KGRBG2, GARBAG
            END DO
         END IF
      END IF
C
C     *** Finally setting number of calculations done. ***
C
      NMDPRP = KNMCLC*KNMPRP
C
      RETURN
      END
C
C
C     /* Deck rdhdrs */
      SUBROUTINE RDHDRS(SYMCOR,RNNORM,FREQ,ICRIRP,NCOOR,NDCOOR,LURSTR)
C     *********************************************************
C     *** Subroutine that reads in header of force constant ***
C     *** restart routine.                                  ***
C     *********************************************************
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION SYMCOR(NCOOR,NCOOR), RNNORM(NCOOR), FREQ(NCOOR)
      DIMENSION ICRIRP(NCOOR,2)

C
C     *** Number of normal coordinates ***
C
      READ(LURSTR,FMT='(I8)') NDCOOR
C
C     *** The normal coordinates. ***
C
      DO IC2 = 1, NDCOOR
         READ(LURSTR,'(2I5)') (ICRIRP(IC2,I),I=1,2)
         DO IC1 = 1, NCOOR
            READ(LURSTR,'(F24.16)') SYMCOR(IC1,IC2)
         END DO
      END DO
C
C     *** Norm of the non-normalized normal coordinates. ***
C
      READ(LURSTR,FMT='(A)')
      DO IC = 1, NDCOOR
         READ(LURSTR,'(F24.16)') RNNORM(IC)
      END DO
C
C     *** Frequencies. ***
C
      READ(LURSTR,FMT='(A)')
      DO IC = 1, NDCOOR
         READ(LURSTR,'(F24.16)') FREQ(IC)
      END DO
C
      RETURN
      END
C
C     /*Deck nrmiso*/
      SUBROUTINE NRMISO(TDER,SYMCOR,DKIN,TRNCCR,TRAMSS,TMPGRD,TMPHES,
     &                  TMPMSS,TMPTD1,TMPTD2,CSTART,WORK,NDIMT,LWORK,
     &                  IPRINT)
C     ***************************************************************
C     *** This routine takes the force constants (with respect to ***
C     *** the most normal masses) and finds the force constants   ***
C     *** with respect to other isotopes).                        ***
C     ***************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <trkoor.h>
#include <numder.h>
#include <taymol.h>
      LOGICAL HESEXS
      INTEGER BEGIN
      DIMENSION TDER(NDIMT), SYMCOR(NCOOR,NCOOR), TRNCCR(NCOOR,NCOOR),
     &          TRAMSS(NCOOR), DKIN(NCOOR), TMPGRD(NCOOR),
     &          TMPHES(NCOOR,NCOOR), TMPMSS(NCOOR), CSTART(NCOOR),
     &          TMPTD1(NCOOR,NCOOR,NCOOR), TMPTD2(NCOOR,NCOOR,NCOOR)
C
      KDIMH = NCOOR**2
      KDIMT = NCOOR**3
C
C     *** Making the third derivative ready for transformation. ***
C
      NCOUNT = 0
      DO 100 IC3 = 1, NCOOR
      DO 100 IC2 = 1, IC3
      DO 100 IC1 = 1, IC2
         NCOUNT = NCOUNT + 1
         TMPTD1(IC1,IC2,IC3) =  TDER(NCOUNT)
         TMPTD1(IC1,IC3,IC2) =  TDER(NCOUNT)
         TMPTD1(IC2,IC1,IC3) =  TDER(NCOUNT)
         TMPTD1(IC2,IC3,IC1) =  TDER(NCOUNT)
         TMPTD1(IC3,IC2,IC1) =  TDER(NCOUNT)
         TMPTD1(IC3,IC1,IC2) =  TDER(NCOUNT)
 100  CONTINUE
C
C     *** Averaged masses needed to be used to transform ***
C     ***          to mass-weighted coordinates.         ***
C
      CALL DZERO(TMPMSS,NDCOOR)
      DO 200 IC2 = 1, NDCOOR
      DO 200 IC1 = 1, NCOOR
         TMPMSS(IC2) = TMPMSS(IC2) + (SYMCOR(IC1,IC2)/DKIN(IC1))**2
 200  CONTINUE
C
C     *** Mass weigting (with the average masses) the force constants ***
C     *** to transform to mass-weighted coordinates.                  ***
C
C     *** Gradient. ***
      DO 300 IC1 = 1, NDCOOR
         GRDMOL(IC1) = GRDMOL(IC1)/SQRT(TMPMSS(IC1))
 300  CONTINUE
C     *** Hessian. ***
      DO 400 IC1 = 1, NDCOOR
         HESMOL(IC1,IC1) = HESMOL(IC1,IC1)/TMPMSS(IC1)
 400  CONTINUE
C     *** Third derivative ***
      DO 500 IC3 = 1, NDCOOR
      DO 500 IC2 = 1, NDCOOR
      DO 500 IC1 = 1, NDCOOR
         TMPTD1(IC1,IC2,IC3) =
     &   TMPTD1(IC1,IC2,IC3)/(SQRT(TMPMSS(IC1)*TMPMSS(IC2)*TMPMSS(IC3)))
 500  CONTINUE
C
C     *** Transforming the force constants back to mass weighted ***
C     ***                cartesian coordinates.                  ***
C
C     *** Gradient ***
      CALL DZERO(TMPGRD,NCOOR)
      DO 600 IC2 = 1, NDCOOR
      DO 600 IC1 = 1, NCOOR
         TMPGRD(IC1) = TMPGRD(IC1) + TRNCCR(IC1,IC2)*GRDMOL(IC2)
 600  CONTINUE
C     *** Hessian ***
      CALL DZERO(TMPHES,KDIMH)
      DO 700 IC3 = 1, NCOOR
      DO 700 IC2 = 1, NDCOOR
      DO 700 IC1 = 1, NDCOOR
         TMPHES(IC1,IC3) = TMPHES(IC1,IC3)
     &                   + HESMOL(IC1,IC2)*TRNCCR(IC3,IC2)
 700  CONTINUE
      CALL DZERO(HESMOL,MXCOOR**2)
      DO 800 IC3 = 1, NCOOR
      DO 800 IC2 = 1, NDCOOR
      DO 800 IC1 = 1, NCOOR
         HESMOL(IC1,IC3) = HESMOL(IC1,IC3)
     &                   + TRNCCR(IC1,IC2)*TMPHES(IC2,IC3)
 800  CONTINUE
C     *** Third derivative ***
      CALL DZERO(TMPTD2,KDIMT)
      DO 900 IC4 = 1, NDCOOR
      DO 900 IC3 = 1, NDCOOR
      DO 900 IC2 = 1, NDCOOR
      DO 900 IC1 = 1, NCOOR
         TMPTD2(IC1,IC3,IC4) = TMPTD2(IC1,IC3,IC4)
     &                       + TRNCCR(IC1,IC2)*TMPTD1(IC2,IC3,IC4)
 900  CONTINUE
      CALL DZERO(TMPTD1,KDIMT)
      DO 1000 IC4 = 1, NDCOOR
      DO 1000 IC3 = 1, NDCOOR
      DO 1000 IC2 = 1, NCOOR
      DO 1000 IC1 = 1, NCOOR
         TMPTD1(IC1,IC2,IC4) = TMPTD1(IC1,IC2,IC4)
     &                       + TRNCCR(IC2,IC3)*TMPTD2(IC1,IC3,IC4)
 1000 CONTINUE
      CALL DZERO(TMPTD2,KDIMT)
      DO 1100 IC4 = 1, NDCOOR
      DO 1100 IC3 = 1, NCOOR
      DO 1100 IC2 = 1, NCOOR
      DO 1100 IC1 = 1, NCOOR
         TMPTD2(IC1,IC2,IC3) = TMPTD2(IC1,IC2,IC3)
     &                       + TRNCCR(IC3,IC4)*TMPTD1(IC1,IC2,IC4)
 1100 CONTINUE
C
C     *** Mass transformation. ***
C
C     *** Gradient ***
      CALL DZERO(GRDMOL,MXCOOR)
      DO 1200 IC1 = 1, NCOOR
         GRDMOL(IC1) = GRDMOL(IC1) + TRAMSS(IC1)*TMPGRD(IC1)
 1200 CONTINUE
C     *** Hessian ***
      CALL DZERO(TMPHES,KDIMH)
      DO 1300 IC2 = 1, NCOOR
      DO 1300 IC1 = 1, NCOOR
         TMPHES(IC1,IC2) = TMPHES(IC1,IC2)+ HESMOL(IC1,IC2)*TRAMSS(IC2)
 1300 CONTINUE
      CALL DZERO(HESMOL,MXCOOR**2)
      DO 1400 IC2 = 1, NCOOR
      DO 1400 IC1 = 1, NCOOR
         HESMOL(IC1,IC2) = HESMOL(IC1,IC2)+ TRAMSS(IC1)*TMPHES(IC1,IC2)
 1400 CONTINUE
C     *** Third derivative. ***
      CALL DZERO(TMPTD1,KDIMT)
      DO 1500 IC3 = 1, NCOOR
      DO 1500 IC2 = 1, NCOOR
      DO 1500 IC1 = 1, NCOOR
         TMPTD1(IC1,IC2,IC3) = TMPTD1(IC1,IC2,IC3)
     &                       + TRAMSS(IC1)*TMPTD2(IC1,IC2,IC3)
 1500 CONTINUE
      CALL DZERO(TMPTD2,KDIMT)
      DO 1600 IC3 = 1, NCOOR
      DO 1600 IC2 = 1, NCOOR
      DO 1600 IC1 = 1, NCOOR
         TMPTD2(IC1,IC2,IC3) = TMPTD2(IC1,IC2,IC3)
     &                       + TRAMSS(IC2)*TMPTD1(IC1,IC2,IC3)
 1600 CONTINUE
      CALL DZERO(TMPTD1,KDIMT)
      DO 1700 IC3 = 1, NCOOR
      DO 1700 IC2 = 1, NCOOR
      DO 1700 IC1 = 1, NCOOR
         TMPTD1(IC1,IC2,IC3) = TMPTD1(IC1,IC2,IC3)
     &                       + TRAMSS(IC3)*TMPTD2(IC1,IC2,IC3)
 1700 CONTINUE
C
C     *** Temporary code, please remove ***
C
      INQUIRE(FILE='DALTON.HES',EXIST=HESEXS)
C
C     *** No hessian specified, we can safely write to file. ***
      IF (.NOT. HESEXS) THEN
C
C        *** Open hessian file. ***
         LUHES = -1
         CALL GPOPEN(LUHES,'DALTON.HES','NEW',' ','FORMATTED',IDUMMY,
     &               .FALSE.)
C
C        *** Printing necessary pre-hessian information. ***
         WRITE(LUHES,'(A)') 'CARTESIAN HESSIAN'
         WRITE(LUHES,*) NCOOR
         WRITE(LUHES,'(A)') '                                   '
C
         DO 1800 ICOOR2 = 1, NCOOR
            DO 1900 ICOOR1 = 1, NCOOR
               WRITE (LUHES,'(F22.12)') HESMOL(ICOOR1,ICOOR2)
 1900       CONTINUE
            WRITE (LUHES,'(A)') '                               '
 1800    CONTINUE
      ELSE
         WRITE (LUPRI,'(//A/A//)')
     &      'Hessian file "DALTON.HES" already exists.' //
     &        ' This file will NOT be overwritten.',
     &      'Please restart the calculation without this file.'
            CALL QUIT('"DALTON.HES" already exists. See output.')
      END IF
C
C     *** Print section. ***
C
C     *** Gradient. ***
      WRITE (LUPRI,'(A)') 'Gradient in symmetry coordinates.'
      DO IC =1, NCOOR
         WRITE (LUPRI,'(F12.7)') GRDMOL(IC)
      END DO
C
C     *** Hessian. ***
C
      BEGIN = 1
      KCOL = 6
      LAST  = MIN(NCOOR,KCOL)
      KCOOR = NCOOR
      NCOL  = NCOOR/KCOL
      IF (MOD(NCOOR,KCOL).NE.0) NCOL = NCOL + 1
C
      DO 101 ICOL = 1, NCOL
         WRITE (LUPRI,1001) (I,I = BEGIN,LAST)
C
         DO 201 ICOOR = BEGIN, NCOOR
            WRITE (LUPRI,2001) ICOOR,
     &                        (HESMOL(ICOOR,I),I=BEGIN,MIN(LAST,ICOOR))
 201     CONTINUE
         WRITE (LUPRI,'(A)') '                                       '
         BEGIN = BEGIN + KCOL
         LAST  = MIN(NCOOR,KCOL+LAST)
 101  CONTINUE
 1001 FORMAT (8X,6(3X,I4,5X),(3X,I4,5X))
 2001 FORMAT (1X,I4,2X,6F12.6)
C
C     *** Qubic force field. ***
C
      CALL HEADER('Third derivative in symmetry coordinates.',-1)
C
      IF (MOD(NCOOR,6).EQ.0) THEN
         NLCMAX = NCOOR/6
      ELSE
         NLCMAX = INT(NCOOR/6)+1
      END IF
C
      DO 202 ICOL2 = 1, NCOOR
         WRITE (LUPRI,'(A,I3)') '      Coloumn number', ICOL2
         WRITE (LUPRI,'(A)') '      -----------------'
         INLC = 0
         DO 402 INLCMX = 1, NLCMAX
            INLC2 = 6*(INLCMX-1) + 1
            INLC  = MIN(INLC+6,NCOOR)
            DO 302 ICOL1 = 1, NCOOR
               WRITE (LUPRI,'(A,6F10.6)') '   ',
     &                       (TMPTD1(I,ICOL1,ICOL2),I=INLC2,INLC)
 302        CONTINUE
            WRITE (LUPRI,'(A)') '                              '
 402     CONTINUE
 202  CONTINUE
C
      IF (IPRINT .GT. 22) THEN
         WRITE (LUPRI,'(A)') 'Inverse of averaged mass'
         WRITE (LUPRI,'(24F10.7)') (TMPMSS(I),I=1,NDCOOR)
      END IF
C
      RETURN
      END
C
C
C     /*Deck dalchg*/
      SUBROUTINE DALCHG(INDSTP,ICRIRP,IRSRDR,IPRINT,NCOOR,NTORDR,FIRST)
C     ***********************************************************
C     *** Routine that reduces the symmetry in the DALTON.INP ***
C     *** file, according to the distortions.                 ***
C     ***********************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <molinp.h>
#include <fcsym.h>
#include <ccorb.h>
      LOGICAL FIRST
      CHARACTER*80 WORD(KMLINE)
      DIMENSION ICRIRP(NCOOR,2), INDSTP(NTORDR)
C
      CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
C
      ILINE = 0
      REWIND (LUCMD,IOSTAT=IOS)
 100  READ (LUCMD,'(A)',ERR=2000) WORD(ILINE+1)
      CALL UPCASE(WORD(ILINE+1))
         ILINE = ILINE + 1
         IF (INDEX(WORD(ILINE),'*END OF').GT.0) GOTO 200
      GOTO 100
C
 200  CONTINUE
      DO 300 I = 1, ILINE
         IF (WORD(I)(1:7) .EQ. '.NSYM  ') THEN
            IF (FIRST) THEN
               READ (WORD(I+1),*) NSMBKP
            ELSE
               NSYM = NSMBKP
               DO 400 J = 1, IRSRDR+1
                  IF ((ICRIRP(INDSTP(J),1).NE.1).AND.(NSYM.GT.1)) THEN
                     NSYM = NSYM/2
                  END IF
 400           CONTINUE
               WRITE (WORD(I+1),'(I4)') NSYM
            END IF
         END IF
 300  CONTINUE
C
      REWIND (LUCMD,IOSTAT=IOS)
      DO 500 I = 1, ILINE
         WRITE (LUCMD,'(A)') WORD(I)
 500  CONTINUE
      CALL GPCLOSE(LUCMD,'KEEP')
C
      RETURN
 2000 CONTINUE
      CALL QUIT('There are problems in correcting the .NSYM parameter')
      END
C
C
C     /* Deck rdhess */
      SUBROUTINE RDHESS(SYMCOR,CSTART,GRIREP,CHRCTR,HESSIN,WORK,ICRIRP,
     &                  LWORK,IPRINT,SYMDET)
C     **********************************************************
C     **** Subroutine that reads in a precalculated hessian ****
C     **** uses this for further work in normal coordinates.****
C     **********************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <fcsym.h>
#include <trkoor.h>
#include <taymol.h>
      LOGICAL SYMDET, SYMADA
      DIMENSION SYMCOR(NCOOR ,NCOOR ), CSTART(NCOOR        ),
     &          GRIREP(NGORDR,NGVERT), CHRCTR(NGORDR,NCVERT),
     &          HESSIN(NCOOR ,NCOOR),  WORK  (LWORK)
      DIMENSION ICRIRP(NCOOR,2)
C
C     *** Initializing. ***
      SYMADA = .FALSE.
C
C     *** Making symmetry adapted coordinates. ***
C
      CALL GRPCHR(CSTART,SYMCOR,GRIREP,CHRCTR,WORK,ICRIRP,LWORK,IPRINT)
C
C     *** Declaring that the symmetry of the system is determined. ***
      SYMDET = .FALSE.
C
C     *** Reading in the hessian ***
C
      CALL RDFHES(WORK,LWORK,IPRINT,SYMADA)
C
      IF (.NOT. SYMADA) THEN
C
C        *** The hessian is cartesian coordinates. Transform it to ***
C        *** Symmetry adapted coordinates.                         ***
C
         KTMPHS = 1
         KLAST  = KTMPHS + NCOOR**2
         LWRK   = LWORK - KLAST
         CALL TRGHES(HESSIN,SYMCOR,WORK(KTMPHS),WORK(KLAST),NCOOR,LWRK,
     &               IPRINT,'symmetry ')
C
         DO 100 ICOOR2 = 1, NCOOR
         DO 100 ICOOR1 = 1, NCOOR
            HESMOL(ICOOR1,ICOOR2) = HESSIN(ICOOR1,ICOOR2)
 100     CONTINUE

      END IF
C
C     *** Print! ***
C
      iprint = 8
      IF (IPRINT .GT. 7) THEN
         CALL HEADER ('Final hessian from RDHESS',0)
C
         N = 0
         IF (MOD(NCOOR,6).NE.0) N = 1
         NCOL = NCOOR/6 + N
         NSTART = 1
         NEND   = MIN(NCOOR,6)
         DO I = 1, NCOL
            DO ICOOR1 = 1, NCOOR
               WRITE (LUPRI,'(6F14.6)')
     &            (HESMOL(ICOOR1,ICOOR2), ICOOR2= NSTART, NEND)
            END DO
            NSTART = NEND + 1
            NEND   = MIN(NCOOR,NEND+6)
            WRITE (LUPRI,'(A)') '                              '
            WRITE (LUPRI,'(A)') '                              '
         END DO
      END IF
      iprint = 0
C
      RETURN
      END
C
C     /* Deck rdfhes */
      SUBROUTINE RDFHES(WORK,LWORK,IPRINT,SYMADA)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <taymol.h>
#include <nuclei.h>
#include <trkoor.h>
      LOGICAL SYMADA, HESEXS
      CHARACTER*5 HSMINF
      DIMENSION WORK(LWORK)
C
      INQUIRE(FILE='DALTON.HES',EXIST=HESEXS)
C
C     *** No hessian specified. ***
      IF (.NOT. HESEXS) CALL QUIT('Unable to open the file DALTON.HES.')
C
C     *** Open hessian file. ***
      LUHES = -1
      CALL GPOPEN(LUHES,'DALTON.HES','OLD',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
C
C     *** Specified hessian in symmetry coordinates? ***
      READ(LUHES,'(A5)') HSMINF
      IF (HSMINF .EQ. 'SYMME') SYMADA = .TRUE.
C
C     *** Check if the speciefied dimensions match those from MOLECULE.INP. ***
      READ(LUHES,*) IDIM
      IF (IDIM .NE. 3*NUCDEP) CALL QUIT('Dimensions for specified ' //
     &     'Hessian does not match those found from the molecule-file.')
      READ(LUHES,*)

C
C     *** Read the hessian from file. ***
C
      DO 100 ICOOR2 = 1, NCOOR
C
         DO 200 ICOOR1 = 1, NCOOR
            READ(LUHES,*) HESMOL(ICOOR1,ICOOR2)
 200     CONTINUE
C
         READ(LUHES,*)
 100  CONTINUE
C
      CALL GPCLOSE(LUHES,'KEEP')
C
      CALL HEADER ('Molecular Hessian read from file.', 0)
      WRITE (LUPRI,'(A)') '                              '
      CALL HEADER ('Molecular Hessian', -1)
      call flshfo(lupri)
C
      N = 0
      IF (MOD(NCOOR,6).NE.0) N = 1
      NCOL = NCOOR/6 + N
      NSTART = 1
      NEND   = MIN(NCOOR,6)
      DO I = 1, NCOL
         DO ICOOR1 = 1, NCOOR
            WRITE (LUPRI,'(6F14.6)')
     &            (HESMOL(ICOOR1,ICOOR2), ICOOR2= NSTART, NEND)
         END DO
         NSTART = NEND + 1
         NEND   = MIN(NCOOR,NEND+6)
         WRITE (LUPRI,'(A)') '                              '
         WRITE (LUPRI,'(A)') '                              '
      END DO
C
      CALL FLSHFO(LUPRI)
C
      RETURN
      END
C
C     /* Deck trahes */
      SUBROUTINE TRGHES(HESSIN,SYMCOR,TMPHES,WORK,NCOOR,LWORK,IPRINT,
     &                  TYPE)
C     **********************************************************
C     *** Transforming hessian in cartesian coordinates, to ****
C     *** symmetry coordinate basis in SYMCOR.              ****
C     **********************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <taymol.h>
      CHARACTER*9 TYPE
      DIMENSION HESSIN(NCOOR,NCOOR), SYMCOR(NCOOR,NCOOR),
     &          TMPHES(NCOOR,NCOOR), WORK(LWORK)
C
C     *** Cartesian to symmetric transformation. ***
      IF (TYPE .EQ. 'symmetry ') THEN
C
        KDIM = NCOOR**2
C
         CALL DZERO(TMPHES,KDIM)
         DO 100 ICOOR3 = 1, NCOOR
         DO 100 ICOOR2 = 1, NCOOR
         DO 100 ICOOR1 = 1, NCOOR
            TMPHES(ICOOR1,ICOOR3) = TMPHES(ICOOR1,ICOOR3)
     &                     + SYMCOR(ICOOR2,ICOOR1)*HESMOL(ICOOR2,ICOOR3)
 100     CONTINUE
C
         CALL DZERO(HESSIN,KDIM)
         DO 200 ICOOR3 = 1, NCOOR
         DO 200 ICOOR2 = 1, NCOOR
         DO 200 ICOOR1 = 1, NCOOR
            HESSIN(ICOOR1,ICOOR3) = HESSIN(ICOOR1,ICOOR3)
     &                     + TMPHES(ICOOR1,ICOOR2)*SYMCOR(ICOOR2,ICOOR3)
 200     CONTINUE
      END IF
C
C     *** Print ***
C
      IF (IPRINT .GT. 7) THEN
         CALL HEADER ('Molecular Hessian in ' // TYPE //
     &                'coordinates, from TRAHES.', 0)
C
         N = 0
         IF (MOD(NCOOR,6).NE.0) N = 1
         NCOL = NCOOR/6 + N
         NSTART = 1
         NEND   = MIN(NCOOR,6)
         DO I = 1, NCOL
            DO ICOOR1 = 1, NCOOR
               WRITE (LUPRI,'(6F14.6)')
     &            (HESSIN(ICOOR1,ICOOR2), ICOOR2= NSTART, NEND)
            END DO
            NSTART = NEND + 1
            NEND   = MIN(NCOOR,NEND+6)
            WRITE (LUPRI,'(A)') '                              '
            WRITE (LUPRI,'(A)') '                              '
         END DO
      END IF
C
      RETURN
      END
C
C     /* Deck priprp */
      SUBROUTINE PRIPRP
C     ****************************************
C     *** Print routine for analyzing part ***
C     *** of the numerical derivatives.    ***
C     ****************************************
#include <implicit.h>
#include <priunit.h>
#include <cbinum.h>
#include <prpndr.h>
C
      CALL HEADER('Analysis using the numerical derivatives',0)
      WRITE (LUPRI,'(/A/)') ' Properties that are analyzed: '
C
      IF (NUMVIB) WRITE (LUPRI,'(A)')
     &   ' - Frequency analysis and efficient geometry.'
      IF (NSPNSP) WRITE (LUPRI,'(A)')
     &   ' - Vibrational average of spin-spin coupling constants'
C
      RETURN
      END
C
C
C     /*Deck stppvr*/
      SUBROUTINE STPPVR
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
C
#include <inforb.h>
#include <cbiexc.h>
#include <pvibav.h>
#include <symmet.h>
#include <ccexcinf.h>
#include <gnrinf.h>
C
      DODIPS = DIPSTR
C
C
      NTOTEX = 0
      DO ISYM = 1, MAXREP+1
         IF (DOCCSD) THEN
            NTOTEX = NTOTEX + (NCCEXCI(ISYM,1))
            NEXCTB(ISYM) = NCCEXCI(ISYM,1)
         ELSE
            NTOTEX = NTOTEX + (NEXCIT(ISYM))
            NEXCTB(ISYM) = NEXCIT(ISYM)
         END IF
      END DO
      EXCIT = (NTOTEX.NE.0)
C
      RETURN
      END
C
C
C     /*Deck prpder*/
      SUBROUTINE PRPDER(SYMCOR,SPSPDR,COEFF,SPSPFV,TRLNFV,TRLNDR,EXENFV,
     &           CCPRFV,CCPRDR,GRIREP,WORK,IADRSS,KDPMTX,ICRIRP,INDSTP,
     &           IDCOMP,IMAX,IMIN,ICNT,NCVAL,IDDCMP,MXCOEF,NTYPE,NPPDER,
     &           LDPMTX,IFRSTD,NLDPMX,MXNEXI,NSYM,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (DMIN = 1.0D-12, D1=1.0D0, D0=0.0D0)
#include <abainf.h>
#include <trkoor.h>
#include <numder.h>
#include <fcsym.h>
#include <cbinum.h>
#include <pvibav.h>
#include <prpc.h>
      LOGICAL PRIVAL,CCPRP
      DIMENSION COEFF(-MXCOEF:MXCOEF,0:NMRDRP), SYMCOR(NCOOR,NCOOR),
     &          SPSPDR(NCOOR ,NCOOR,6,NPPDER), GRIREP(NGORDR,NGVERT),
     &          SPSPFV(NCOOR,NCOOR,6,NMPCAL),EXENFV(NSYM,MXNEXI,NMPCAL),
     &          TRLNFV(3,NSYM,MXNEXI,NMPCAL), CCPRFV(NPRPC,NMPCAL),
     &          TRLNDR(3,NSYM,MXNEXI,NPPDER), CCPRDR(NPRPC,NPPDER),
     &          WORK(LWORK)
      DIMENSION ICNT(NTYPE), IADRSS(NTYPE), IMAX(NMRDRP), IMIN(NMRDRP),
     &          INDSTP(NTORDR), INDTMP(NTORDR), IDCOMP(NCOOR),
     &          IDDCMP(NCOOR), NCVAL(NCOOR),
     &          KDPMTX(LDPMTX,NSTRDR,IFRSTD), ICRIRP(NCOOR,2)
C
C     *************************************
C     *** Reading properties from file. ***
C     *************************************
C
      KTRAMA = 1
      KNREDS = KTRAMA +   NMPCAL*NCOOR**2
      KEXTMP = KNREDS +   NSYM**2
      KTRTMP = KEXTMP +   NSYM*NTOTEX
      KEXERF = KTRTMP + 3*NSYM*NTOTEX
      KLAST  = KEXERF +   NSYM*MXNEXI
      LWRK1  = LWORK - KLAST + 1
      CALL NDRDPP(SPSPFV,TRLNFV,EXENFV,CCPRFV,WORK(KTRAMA),SYMCOR,
     &            WORK(KEXTMP),WORK(KTRTMP),WORK(KEXERF),
     &            WORK(KLAST),WORK(KNREDS),LWRK1,IPRINT,CCPRP)
      IF (NPRPDR) CALL GPCLOSE(LUNDPR,'DELETE')
c      KTRAMA = 1
c      KLAST  = KTRAMA + NMPCAL*NCOOR**2
c      LWRK1  = LWORK - KLAST + 1
c      CALL NDRDPP(SPSPFV,TRLNFV,EXENFV,WORK(KTRAMA),SYMCOR,WORK(KLAST),
c     &            LWRK1,IPRINT,DODIPS)
c      IF (NPRPDR) CALL GPCLOSE(LUNDPR,'DELETE')
C
C     ******************************************
C     *** Calculating numerical derivatives. ***
C     ******************************************
C
C     *** For cc-properties. ***
C
      IF (CCPRP) THEN
         NFINNR = NPRPC
         CALL NMNDER(CCPRDR,COEFF,CCPRFV,GRIREP,WORK,IADRSS,KDPMTX,
     &               ICRIRP,INDSTP,INDTMP,IDCOMP,IMAX,IMIN,ICNT,NCVAL,
     &               IDDCMP,MXCOEF,NMRDRP,NMPCAL,NTYPE,NPPDER*NFINNR,
     &               NFINNR,LDPMTX,IFRSTD,NLDPMX,LWORK,.FALSE.)
C
C        *** Double check sign on excited states property derivatives. ***
C
         KTPCCD = 1
         KTPCCF = KTPCCD + 2*NPRPC*NPPDER
         KLAST  = KTPCCF +   NPRPC*NMPCAL
         LWRK1  = LWORK  - KLAST + 1
         CALL CHK1DR(CCPRFV,CCPRDR,WORK(KTPCCD),WORK(KTPCCF),COEFF,
     &               GRIREP,WORK(KLAST),ICNT,IADRSS,IMAX,IMIN,
     &               INDSTP,INDTMP,IDCOMP,IDDCMP,NCVAL,KDPMTX,ICRIRP,
     &               NPPDER,MXCOEF,NTYPE,NFINNR,LDPMTX,IFRSTD,
     &               NLDPMX,LWRK1,IPRINT)
      END IF
C
C     *** For spin-spin. ***
C
      IF (SPNSPN) THEN
         NFINNR = 6*NCOOR**2
         CALL NMNDER(SPSPDR,COEFF,SPSPFV,GRIREP,WORK,IADRSS,KDPMTX,
     &               ICRIRP,INDSTP,INDTMP,IDCOMP,IMAX,IMIN,ICNT,NCVAL,
     &               IDDCMP,MXCOEF,NMRDRP,NMPCAL,NTYPE,NPPDER*NFINNR,
     &               NFINNR,LDPMTX,IFRSTD,NLDPMX,LWORK,.FALSE.)
      END IF
C
C     *** For transition dipole moments. ***
C
      IF (DODIPS) THEN
         NFINNR = 3*NSYM*MXNEXI
         CALL NMNDER(TRLNDR,COEFF,TRLNFV,GRIREP,WORK,IADRSS,KDPMTX,
     &               ICRIRP,INDSTP,INDTMP,IDCOMP,IMAX,IMIN,ICNT,NCVAL,
     &               IDDCMP,MXCOEF,NMRDRP,NMPCAL,NTYPE,NPPDER*NFINNR,
     &               NFINNR,LDPMTX,IFRSTD,NLDPMX,LWORK,.FALSE.)
      END IF
C
C     ***********************************************
C     *** Write necessary results to file for the ***
C     *** for the vibrational analysis.           ***
C     ***********************************************
C
      IF (PRPVIB) THEN
         NMDPRP = 0
         CALL NDWTPP(SPSPFV,SPSPDR,NPPDER,IPRINT)
      END IF
C
C     ***********************************
C     *** Test print or result print. ***
C     ***********************************
C
      IF ((IPRINT.GE.20).OR.(NPRPDR.AND..NOT.PRPVIB)) THEN
C
C       *** For cc-properties. ***
C
         IF (CCPRP) THEN
            IDERV = 0
            CALL TITLER('Derivatives OF CC-properties.','*',118)
            DO IORDR = 1, NMRDRP
               IF (IORDR.EQ.1) THEN
                  CALL HEADER('1. numerical derivative',0)
                  DO IC = 1, 2
C
                     IF (IC.EQ.2) THEN
                        KCDVAL = 1
                        CALL T1PRSC(CCPRDR,WORK(KCDVAL),SYMCOR,NPRPC,
     &                              NPPDER,IPRINT)
                        KSTART = KCDVAL
                     END IF
C
                     DO ICOOR = 1, NDCOOR
                        IF (IC.EQ.1) THEN
                           IDERV = IDERV + 1
                           WRITE (LUPRI,'(5X,A,I4)')
     &                        'Derivative with respect to ' //
     &                        'symmetry coordinate', ICOOR
                        ELSE
                           WRITE (LUPRI,'(5X,A,I4)')
     &                        'Derivative with respect to ' //
     &                        'cartesian coordinate', ICOOR
                        END IF
C
                        IF (IC.EQ.1) THEN
C
C                          *** Update value for printout in ***
C                          *** symmetry coordinates.        ***
C
                           LUPRPCO = -1
                           CALL GPOPEN(LUPRPCO,'CC_PRPC_O','UNKNOWN',
     &                          ' ','FORMATTED',IDUMMY,.FALSE.)
                           CALL PRPRPC(LUPRPCO,2,CCPRDR(1,IDERV),NPRMI)
                           CALL GPCLOSE(LUPRPCO,'KEEP')
                        ELSE
C
C                          *** Update value for printout in ***
C                          *** cartesian coordinates.       ***
C
                           LUPRPCO = -1
                           CALL GPOPEN(LUPRPCO,'CC_PRPC_O','UNKNOWN',
     &                          ' ','FORMATTED',IDUMMY,.FALSE.)
                           CALL PRPRPC(LUPRPCO,2,WORK(KSTART),NPRMI)
                           CALL GPCLOSE(LUPRPCO,'KEEP')
C
C                          *** Update value for printout. ***
C
                           KSTART = KSTART + NPRPC
                        END IF
                        WRITE (LUPRI,'(/)')
                     END DO
                  END DO
               ELSE IF (IORDR.EQ.2) THEN
C
                  CALL HEADER('2. numerical derivative',0)
C
                  DO ICOOR2 = 1, NDCOOR
                  DO ICOOR1 = 1, ICOOR2
C
                     IDERV = IDERV + 1
C
                     IF (PRPVIB.AND.((NARDRP+NMRDRP).EQ.2)) THEN
                        PRIVAL = ICOOR1.EQ.ICOOR2
                     ELSE
                        PRIVAL = .TRUE.
                     END IF
C
                     IF (PRIVAL) THEN
                        WRITE (LUPRI,'(5X,A,I4,A,I4)')
     &                       'Derivative with respect to coordinate',
     &                       ICOOR2, ' and ', ICOOR1
C
                        LUPRPCO = -1
                        CALL GPOPEN(LUPRPCO,'CC_PRPC_O','UNKNOWN',
     &                              ' ','FORMATTED',IDUMMY,.FALSE.)
                        CALL PRPRPC(LUPRPCO,2,CCPRDR(1,IDERV),NPRMI)
                        CALL GPCLOSE(LUPRPCO,'KEEP')
                     END IF
                  END DO
                  END DO
               END IF
            END DO
         END IF
C
C        *** For spin-spin ***
C
         IF (SPNSPN) THEN
            IDERV = 0
            CALL TITLER('Spin-spin derivatives.','*',118)
            DO IORDR = 1, NMRDRP
               IF (IORDR.EQ.1) THEN
                  CALL HEADER('1. numerical derivative',0)
                  DO ICOOR = 1, NDCOOR
C
                     IDERV = IDERV + 1
C
                     WRITE (LUPRI,'(5X,A,I4)')
     &                    'Derivative with respect to coordinate',
     &                    ICOOR
C
                     CALL PRSPSP(SPSPDR(1,1,1,IDERV),NCOOR,NCOOR,LUPRI)
                  END DO
               ELSE IF (IORDR.EQ.2) THEN
C
                  CALL HEADER('2. numerical derivative',0)
C
                  DO ICOOR2 = 1, NDCOOR
                  DO ICOOR1 = 1, ICOOR2
C
                     IDERV = IDERV + 1
C
                     IF (PRPVIB.AND.((NARDRP+NMRDRP).EQ.2)) THEN
                        PRIVAL = ICOOR1.EQ.ICOOR2
                     ELSE
                        PRIVAL = .TRUE.
                     END IF
C
                     IF (PRIVAL) THEN
                        WRITE (LUPRI,'(5X,A,I4,A,I4)')
     &                       'Derivative with respect to coordinate',
     &                       ICOOR2, ' and ', ICOOR1
C
                        CALL PRSPSP(SPSPDR(1,1,1,IDERV),NCOOR,NCOOR,
     &                              LUPRI)
                     END IF
                  END DO
                  END DO
               END IF
            END DO
         END IF
C
C        *** For transition moments ***
C
         IF (DODIPS) THEN
            IDERV = 0
            CALL TITLER('Transition moment derivatives.','*',118)
            DO IORDR = 1, NMRDRP
               IF (IORDR.EQ.1) THEN
                  CALL HEADER('1. numerical derivative',0)
                  WRITE (LUPRI,'(5X,A)') 'Excitation energies are' //
     &                  ' shown for original geometry.'
                  WRITE (LUPRI,'(5X,A)') '                        '
                  DO ICOOR = 1, NDCOOR
C
                     IDERV = IDERV + 1
C
                     WRITE (LUPRI,'(5X,A,I4)')
     &                    'Derivative with respect to coordinate',
     &                    ICOOR
C
                     CALL PRDPTR(TRLNDR(1,1,1,IDERV),EXENFV(1,1,1),NSYM,
     &                           LUPRI)
                  END DO
               ELSE IF (IORDR.EQ.2) THEN
C
                  CALL HEADER('2. numerical derivative',0)
C
                  DO ICOOR2 = 1, NDCOOR
                  DO ICOOR1 = 1, ICOOR2
C
                     IDERV = IDERV + 1
C
                     IF ((PRPVIB).AND.((NARDRP+NMRDRP).EQ.2)) THEN
                        PRIVAL = ICOOR1.EQ.ICOOR2
                     ELSE
                        PRIVAL = .TRUE.
                     END IF
C
                     IF (PRIVAL) THEN
                        WRITE (LUPRI,'(5X,A,I4,A,I4)')
     &                       'Derivative with respect to coordinate',
     &                       ICOOR2, ' and ', ICOOR1
C
                        CALL PRDPTR(TRLNDR(1,1,1,IDERV),EXENFV(1,1,1),
     &                              NSYM,LUPRI)
                     END IF
                  END DO
                  END DO
               END IF
            END DO
         END IF
      END IF
C
      RETURN
      END
C
C
C     /* Deck prspsp */
      SUBROUTINE PRSPSP(SPSPFV,NDIM1,NDIM2,LPRIUN)
#include <implicit.h>
C
      CHARACTER*26 CNTRIB(6)
      DIMENSION SPSPFV(NDIM1,NDIM2,6)
C
C     *** Different contributions. ***
C
      DATA CNTRIB /'Total spin-spin-coupling. ',
     &             'DSO-contribution.         ',
     &             'PSO-contribution.         ',
     &             'SD-contribution.          ',
     &             'FC-contribution.          ',
     &             'Spin dipole Fermi contact.'/
C
      NTCOL = NDIM1/3 + 1
      IF (MOD(NDIM1,3).EQ.0) NTCOL = NDIM1/3
C
C     *** Printing the cntributions. ***
C
      DO IDIM3 = 1, 6
         KDIM = 0
         CALL HEADER(CNTRIB(IDIM3),-1)
         DO ITCOL = 1, NTCOL
            DO IDIM2 = 1, NDIM2
               WRITE (LPRIUN,'(3F24.16)')
     &              (SPSPFV(IDIM1,IDIM2,IDIM3),
     &                                 IDIM1 = KDIM+1,MIN(KDIM+3,NDIM1))
            END DO
            WRITE (LPRIUN,'(A)') '                          '
            KDIM = KDIM + 3
         END DO
         WRITE (LPRIUN,'(A)') '                          '
      END DO
C
      RETURN
      END
C
C
C     /* Deck prtrma */
      SUBROUTINE PRTRMA(TRAMAT,NDIMT1,NDIMT2,NDIMP1,NDIMP2,LPRIUN)
C     *****************************************************************
C     *** Subroutine that prints a two dimensional matrix (TRAMAT). ***
C     *****************************************************************
#include <implicit.h>
C
      DIMENSION TRAMAT(NDIMT1,NDIMT2)
C
      NTCOL = NDIMP2/6 + 1
      IF (MOD(NDIMP2,6).EQ.0) NTCOL = NDIMP2/6
C
C     *** Printing transformation matrix. ***
C
      KDIM = 0
      DO ITCOL = 1, NTCOL
         DO IDIM2 = 1, NDIMP1
            WRITE(LPRIUN,'(6F10.4)')
     &           (TRAMAT(IDIM2,IDIM1),IDIM1 = KDIM+1,MIN(KDIM+6,NDIMP2))
         END DO
         WRITE (LPRIUN,'(A)') '                          '
         KDIM = KDIM + 6
      END DO
      WRITE (LPRIUN,'(A)') '                          '
C
      RETURN
      END
C
C
C     /* Deck trfcgd */
      SUBROUTINE PRDPTR(TRLEN,EXENG,NSYM,LUPRI)
C     ***************************************************
C     *** Subroutine that prints the dipole transition***
C     *** moments.                                    ***
C     ***************************************************
#include <implicit.h>
C
#include <cbiexc.h>
      DIMENSION TRLEN(3,NSYM,MXNEXI), EXENG(NSYM,MXNEXI)
C
      CALL HEADER('Electric transition dipole moments (in a.u.)',15)
      WRITE (LUPRI,'(1X,A,A,2(/,1X,A))')
     &   ' Sym.   Mode    Frequency  ',
     &   '                   Length       ',
     &   'ex. st.  No.      (au)             x            y         ' //
     &   '    z   ',
     &   '----------------------------------------------------------' //
     &   '---------'
      DO 200 ISYM = 1, NSYM
         DO 100 IEXVAL = 1,NEXCIT(ISYM)
            WRITE (LUPRI,'(2X,I2,6X,I3,1X,F12.6,2X,3F13.5)')
     &            ISYM, IEXVAL, EXENG(ISYM,IEXVAL),
     &            TRLEN(1,ISYM,IEXVAL),  TRLEN(2,ISYM,IEXVAL),
     &            TRLEN(3,ISYM,IEXVAL)
 100     CONTINUE
 200  CONTINUE
      WRITE (LUPRI,'()')
      WRITE (LUPRI,'()')
      WRITE (LUPRI,'()')
      WRITE (LUPRI,'()')
C
      RETURN
      END
C
C
C     /* Deck trfcgd */
      SUBROUTINE TRFCGD(EGRAD,SYMCOR,COOR,SEGRAD,WORK,NCOOR1,NCOOR2,
     &                  LWORK,IPRINT)
C     ************************************************************
C     **** Subroutine that transforms a gradien in a set of   ****
C     **** cartesian coordinates, via the permutation of the  ****
C     **** atoms used by the numerical differentiation, to the****
C     **** symmetry adapted coordinates used by the numerical ****
C     **** differentiation scheme.                            ****
C     ************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      DIMENSION SYMCOR(NCOOR1,NCOOR1), EGRAD(MXCOOR), SEGRAD(NCOOR1),
     &          COOR(NCOOR1), WORK(LWORK)
C
C     *** Transforming into the "old" set of cartesian ***
C     *** coordinates.                                 ***
C
      KCRPRG = 1
      KTRAMT = KCRPRG + NCOOR1
      CALL TROCGD(EGRAD,COOR,SEGRAD,WORK(KCRPRG),WORK(KTRAMT),NCOOR1,
     &            IPRINT)
C
C     *** Transforming into symmetry coordinates. ***
C
      CALL TRSFCG(EGRAD,SYMCOR,SEGRAD,NCOOR1,NCOOR2,IPRINT)
C
      RETURN
      END
C
C
C     /* Deck trocgd*/
      SUBROUTINE TROCGD(EGRAD,COOR,TMPGRD,CRTPRG,TRAMAT,NCOOR,IPRINT)
C     ***********************************************************
C     **** Subroutine that transforms the gradient in the set****
C     **** of cartesian coordinates, to another set of       ****
C     **** cartesian coordinates (stored in coor).           ****
C     ***********************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D1 = 1.0D0, D0 = 0.0D0)
C
      DIMENSION EGRAD (MXCOOR), COOR  (NCOOR      ), CRTPRG(NCOOR),
     &          TMPGRD(NCOOR ), TRAMAT(NCOOR,NCOOR)
C
C     *** Constructing the transformation matrix. ***
C
      CALL TRMTOC(TRAMAT,COOR,CRTPRG,NCOOR,IPRINT)
C
C     *** Transforming the gradient matrix. ***
C
      CALL DGEMM('N','N',NCOOR,1,NCOOR,D1,TRAMAT,NCOOR,EGRAD,MXCOOR,
     &           D0,TMPGRD,NCOOR)
C
      CALL DCOPY(NCOOR,TMPGRD,1,EGRAD,1)
C
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Test-printing of gradient in new cart. coor.',0)
C
         DO ICOOR = 1, NCOOR
            WRITE (LUPRI,'(5X,F12.6)') EGRAD(ICOOR)
         END DO
      END IF
C
      RETURN
      END
C
C
C     /* Deck trsfcg */
      SUBROUTINE TRSFCG(EGRAD,SYMCOR,SEGRAD,NCOOR1,NCOOR2,IPRINT)
C     ***************************************************************
C     **** Subroutine that transforms a gradient in the cartesian****
C     **** coordinates used by the numerical differentiation, to ****
C     **** the symmetry adapted coordinates used by the numerical****
C     **** differentiation scheme.                               ****
C     ***************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D1 = 1.0D0, D0 = 0.0D0)
C
      DIMENSION SYMCOR(NCOOR1,NCOOR1), EGRAD(MXCOOR), SEGRAD(NCOOR1)
C
      CALL DGEMM('T','N',NCOOR2,1,NCOOR1,D1,SYMCOR,NCOOR1,EGRAD,MXCOOR,
     &           D0,SEGRAD,NCOOR1)
C
      CALL DCOPY(NCOOR2,SEGRAD,1,EGRAD,1)
C
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Test-printing of gradient in sym. coordinates',0)
         DO ICOOR1 = 1, NCOOR2
            WRITE (LUPRI,*) SEGRAD(ICOOR1)
         END DO
      END IF
C
      RETURN
      END
C
C
C     /* Deck trfchs */
      SUBROUTINE TRFCHS(EHESS,SYMCOR,COOR,SEHESS,WORK,NCOOR1,NCOOR2,
     &                  LWORK,IPRINT)
C     ************************************************************
C     **** Subroutine that transforms a hessian in a set of   ****
C     **** cartesian coordinates, via the permutation of the  ****
C     **** atoms used by the numerical differentiation, to the****
C     **** symmetry adapted coordinates used by the numerical ****
C     **** differentiation scheme.                            ****
C     ************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D1 = 1.0D0, D0 = 0.0D0, KCOL = 6)
C
      DIMENSION SYMCOR(NCOOR1,NCOOR1), EHESS(MXCOOR,MXCOOR),
     &          COOR(NCOOR1), SEHESS(NCOOR1,NCOOR1), WORK(LWORK)
C
C     *** Transforming into the "old" set of cartesian ***
C     *** coordinates.                                 ***
C
      KCRPRG = 1
      KTRAMT = KCRPRG + NCOOR1
      CALL TROCHS(EHESS,COOR,SEHESS,WORK(KCRPRG),WORK(KTRAMT),NCOOR1,
     &            IPRINT)
C
C     *** Transforming into symmetry coordinates. ***
C
      CALL TRSFC2(EHESS,SYMCOR,SEHESS,NCOOR1,NCOOR2,MXCOOR,IPRINT)
C
      RETURN
      END
C
C
C     /* Deck trochs*/
      SUBROUTINE TROCHS(EHESS,COOR,TMPHES,CRTPRG,TRAMAT,NCOOR,IPRINT)
C     ***********************************************************
C     **** Subroutine that transforms the hessian in the set ****
C     **** of cartesian coordinates, to another set of       ****
C     **** cartesian coordinates (stored in coor).           ****
C     ***********************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D1 = 1.0D0, D0 = 0.0D0, KCOL = 6)
C
      INTEGER BEGIN
      DIMENSION EHESS (MXCOOR,MXCOOR), COOR(NCOOR), CRTPRG(NCOOR),
     &          TMPHES(NCOOR ,NCOOR ), TRAMAT(NCOOR ,NCOOR )
C
C     *** Constructing the transformation matrix. ***
C
      CALL TRMTOC(TRAMAT,COOR,CRTPRG,NCOOR,IPRINT)
C
C     *** Transforming the hessian matrix. ***
C
      CALL OTRTEN(EHESS,TRAMAT,TMPHES,MXCOOR,NCOOR,NCOOR,IPRINT,'N','T')
C

C
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Test-printing of hessian in new cart. coor.',0)
         BEGIN = 1
         LAST  = MIN(NCOOR,KCOL)
         KCOOR = NCOOR
         NCOL  = INT(DBLE(NCOOR)/DBLE(KCOL))
         IF (MOD(NCOOR,KCOL).NE.0) NCOL = NCOL + 1
C
         DO ICOL = 1, NCOL
C
            DO ICOOR = BEGIN, NCOOR
               WRITE (LUPRI,'(2X,6F12.6)')
     &                        (EHESS(ICOOR,I),I=BEGIN,MIN(LAST,ICOOR))
            END DO
            WRITE (LUPRI,'(A)') '                                '
            BEGIN = BEGIN + KCOL
            LAST  = MIN(NCOOR,KCOL+LAST)
         END DO

      END IF
      RETURN
      END
C
C
C     /* Deck trsfc2 */
      SUBROUTINE TRSFC2(SCNDER,SYMCOR,SSCNDR,NCOOR1,NCOOR2,NSCNDR,
     &                  IPRINT)
C     **************************************************************
C     **** Subroutine that transforms a secon derivetive in     ****
C     **** cartesian coordinates used by the numerical          ****
C     **** differentiation, to the symmetry adapted coordinates ****
C     **** used by the numerical differentiation scheme.        ****
C     **************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D0, D0 = 0.0D0, KCOL = 6)
C
      INTEGER BEGIN, LAST
      DIMENSION SYMCOR(NCOOR1,NCOOR1), SCNDER(NSCNDR,NSCNDR),
     &          SSCNDR(NCOOR1,NCOOR1)
C
      CALL DGEMM('T','N',NCOOR2,NCOOR1,NCOOR1,D1,SYMCOR,NCOOR1,
     &           SCNDER,NSCNDR,D0,SSCNDR,NCOOR1)
c      d = 0.0d0
c      do i = 1, ncoor1
c         d = d + SCNDER(1,i)*symcor(i,3)
c         write (lupri,*) SCNDER(1,i), symcor(i,3)
c      end do
c      write (lupri,*) sscndr(3,1),d
c      stop ' '
C
      CALL DGEMM('N','N',NCOOR2,NCOOR2,NCOOR1,D1,SSCNDR,NCOOR1,
     &           SYMCOR,NCOOR1,D0,SCNDER,NSCNDR)
C
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Test-printing of hessian in sym. coordinates',0)
         BEGIN = 1
         LAST  = MIN(NCOOR1,KCOL)
         KCOOR = NCOOR2
         NCOL = NCOOR2/KCOL
         IF (MOD(NCOOR2,KCOL).NE.0) NCOL = NCOL + 1
C
         DO ICOL = 1, NCOL
C
            DO ICOOR = BEGIN, NCOOR2
               WRITE (LUPRI,1000)
     &                        (SCNDER(ICOOR,I),I=BEGIN,MIN(LAST,ICOOR))
            END DO
            WRITE (LUPRI,'(A)') '                                '
            BEGIN = BEGIN + KCOL
            LAST  = MIN(NCOOR2,KCOL+LAST)
         END DO
      END IF
C
      RETURN
 1000 FORMAT (2X,6F12.6)
      END
C
C
C     /* Deck trmtoc */
      SUBROUTINE TRMTOC(TRAMAT,COOR,CRTPRG,NCOOR,IPRINT)
C     **************************************************************
C     **** Subroutine that constructs the transformation matrix ****
C     **** to transform gradient/hessian back to original       ****
C     **** set of coordinates, defined in COOR.                 ****
C     **************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
      PARAMETER (D1 = 1.0D0, DMTHR=1.0D-4)
#include <symmet.h>
#include <nuclei.h>
#include <pvibav.h>
#include <numder.h>
#include <cbinum.h>
       LOGICAL FOUND
       CHARACTER*9 PRPTXT
       DIMENSION TRAMAT(NCOOR,NCOOR), COOR(NCOOR), CRTPRG(NCOOR)
#include <ibtfun.h>
C
C     *** Finding the cartesian coordinates used by the ***
C     *** program at the moment.                        ***
C
      ICOOR = 0
      IATOM = 0
      DO ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         DO IOP = 0, MAXOPR
            IF (IBTAND(IOP,MULCNT) .EQ. 0) THEN
               IATOM = IATOM + 1
               DO I = 1, 3
                  ICOOR = ICOOR + 1
                  CRTPRG(ICOOR) =
     &                 PT(IBTAND(ISYMAX(I,1),IOP))*CORD(I,ICENT)
               END DO
            END IF
         END DO
      END DO
C
C     *** Constructing the transformation matrix by comparing them ***
C     *** to the old set of cartesian coordinates.                 ***
C
      CALL DZERO(TRAMAT,NCOOR**2)
      DO IATOM1 = 1, NATOMS
         ICS1 = 3*(IATOM1-1)
         DO IATOM2 = 1, NATOMS
            ICS2 = 3*(IATOM2-1)
C
            FOUND = .TRUE.
            DO IC = 1, 3
               FOUND = FOUND .AND.
     &                     ((COOR(ICS1+IC)-CRTPRG(ICS2+IC))**2.LT.DMTHR)
            END DO
C
            IF (FOUND) THEN
               DO IC = 1, 3
                  TRAMAT(ICS1+IC,ICS2+IC) = D1
               END DO
            END IF
         END DO
      END DO
C
C     *** If property derivative is calculated, we need to save ***
C     *** the transformation matrix in the property-file.       ***
C
      IF (CNMPRP) THEN
         NDIM3  = 1
         PRPTXT = 'CART-TRAN'
         CALL WRAVFL(TRAMAT,NCOOR,NCOOR,NDIM3,PRPTXT,IPRINT)
      END IF
C
C     *** Test print ***
C
      IF (IPRINT .GT. 50) THEN
         WRITE (LUPRI,'(A)')
     &                  'The nuclear coordinates used by the program:'
         DO ICOOR = 1, NCOOR
            WRITE (LUPRI,'(5X,2F10.5,I4)') CRTPRG(ICOOR)
         END DO
         WRITE (LUPRI,'(A)') '                                       '
         WRITE (LUPRI,'(A)')
     &              'Transforming to using these nuclear coordinates.'
         DO ICOOR = 1, NCOOR
            WRITE (LUPRI,'(5X,F10.5)') COOR(ICOOR)
         END DO
C
         WRITE (LUPRI,'(5X,A)') '                       '
         WRITE (LUPRI,'(5X,A)') 'Transformation matrix: '
         WRITE (LUPRI,'(5X,A)') '                       '
         CALL PRTRMA(TRAMAT,NCOOR,NCOOR,NCOOR,NCOOR,LUPRI)
      END IF
C
      RETURN
      END
C
C
C     /* Deck bksmnm */
      SUBROUTINE BKSMNM
C     ****************************************************
C     *** This routine takes care of symmetry odds and ***
C     *** ends connected to frozen core orbitals in    ***
C     *** distorted symmetry.                          ***
C     ****************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
C
#include <symmet.h>
#include <nmbksym.h>
#include <ccorb.h>
#include <numder.h>
#include <cbinum.h>
      CHARACTER*8 WORD
C
C     *** Backing up symmetry. ***
C
      CALL ICOPY(64    ,IXVAL      ,1,IXVALB,1)
      CALL ICOPY(8     ,JSOP       ,1,JSOPB ,1)
      CALL ICOPY(8     ,NRHFFR     ,1,NRHFRB,1)
      MAXRPB = MAXREP
C
C     *** Making sure that there are no complicating issues, so ***
C     *** that the molecule should not rotate freely.           ***
C
      IF (.NOT.NOMOVE) THEN
         IF (.NOT.((NAORDR.EQ.0).AND.(NMORDR.NE.1).AND.(.NOT.NPRPDR)))
     &                                                  NOMOVE = .TRUE.
      END IF
C
C     *** Different circumstances where molecule ***
C     *** is not allowed to rotate.              ***
C
C
C     *** If symmetry is reported in DALTON.INP ***
C     *** using .NSYM, no rotation is allowed.  ***
C
      CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      ILINE = 0
      REWIND (LUCMD,IOSTAT=IOS)
 100  CONTINUE
      READ (LUCMD,'(A)') WORD
      CALL UPCASE(WORD)
      IF (WORD .EQ. '.NSYM  ') NOMOVE = .TRUE.
      IF((WORD .NE. '*END OF ') .AND.
     &   (WORD .NE. '**END OF')) GOTO 100
      CALL GPCLOSE(LUCMD,'KEEP')
C
C     *** Numerical derivatives of properties, ***
C     *** no rotation is allowed.              ***
C
      IF (NMDPRP.GT.0) NOMOVE = .TRUE.
C
      RETURN
      END
C
C
C     /* Deck fndexs */
      SUBROUTINE FNDEXS(WORD,IPRINT)
C     ******************************************************
C     *** Subroutine that sorts out the symmetry of the  ***
C     *** excited states, and assign them to a new irrep ***
C     *** in the distorted geometry.                     ***
C     ******************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <molinp.h>
C
#include <symmet.h>
#include <nmbksym.h>
#include <ccorb.h>
#include <pgroup.h>
      LOGICAL SAMIRP
      CHARACTER*80 WORD(KMLINE)
C
C     *** Reducing symmetry of frozen orbitals. ***
C
      CALL SDCEIP(NRHFFR,NRHFRB,WORD,'.FROINP')
C
C     *** Print. ***
C
      IF (IPRINT.GT.0) THEN
         WRITE (LUPRI,'(/A)')
     &         'Symmetries of frozen core orbitals in reduced symmetry:'
         WRITE (LUPRI,'(2X,8A4)') (REP(I), I=0,MAXREP)
         WRITE (LUPRI,'(8I4)') (NRHFFR(I), I=1,MAXREP+1)
         WRITE (LUPRI,'(/)')
      END IF
C
      RETURN
      END
C
C
C     /* Deck sdceip */
      SUBROUTINE SDCEIP(NCURNT,NBCKUP,WORD,SWORD)
C     ************************************************************
C     *** Subroutine that subduces properties in input file    ***
C     *** from original symmetry into broken symmetry.         ***
C     *** Original symmetry needs to be backed up in nmbksym.h.***
C     ************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <molinp.h>
C
#include <symmet.h>
#include <nmbksym.h>
#include <pgroup.h>
      LOGICAL SAMIRP
      DIMENSION NCURNT(8), NBCKUP(8)
      CHARACTER*7 SWORD
      CHARACTER*80 WORD(KMLINE)
C
      CALL STCCSM(NSYM)
C
      CALL IZERO(NCURNT,8)
C
C     *** Sorting irep's, and assign them to the appropriate ***
C     *** place in NCURNT.                                   ***
C
      DO 100 IREP = 0, MAXREP
         DO 200 IRPOLD = 0, MAXRPB
            IF (NBCKUP(IRPOLD+1) .GT. 0) THEN
               SAMIRP = .TRUE.
               DO 300 ISYOP2 = 0, MAXREP
               DO 300 ISYOP1 = 0, MAXRPB
                  IF (SYMOP(JSOP(ISYOP2)).EQ.SYMOP(JSOPB(ISYOP1))) THEN
                     IF (IXVALB(JSOPB(ISYOP1),IRPOLD)  .NE.
     &                    IXVAL(JSOP (ISYOP2),IREP )) SAMIRP = .FALSE.
                  END IF
 300           CONTINUE
               IF (SAMIRP) THEN
                  NCURNT(IREP+1) = NCURNT(IREP+1) + NBCKUP(IRPOLD+1)
               END IF
            END IF
 200     CONTINUE
 100  CONTINUE
C
C     *** Writing results to DALTON.INP ***
C
      CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
C
C     *** Updating file.***
C
      ILINE = 0
      REWIND (LUCMD,IOSTAT=IOS)
 400  CONTINUE
      ILINE = ILINE + 1
      READ (LUCMD,'(A)') WORD(ILINE)
      CALL UPCASE(WORD(ILINE))
      IF (WORD(ILINE) .EQ. '.NSYM  ') THEN
         WRITE (WORD(ILINE+1),'(I3)')  NSYM
      ELSE IF (WORD(ILINE) .EQ. SWORD) THEN
         WRITE (WORD(ILINE+1),'(8I4)') (NCURNT(IREP),IREP=1,NSYM)
      END IF
      IF (.NOT.(WORD(ILINE)(1:6) .EQ. '*END O' .OR.
     &          WORD(ILINE)(1:6) .EQ. '**END ')) GOTO 400
C
      REWIND(LUCMD,IOSTAT=IOS)
      DO I = 1, ILINE
         WRITE (LUCMD,'(A)') WORD(I)
      END DO
C
C     *** Closing DALTON.INP. ***
C
      CALL GPCLOSE(LUCMD,'KEEP')
C
      RETURN
      END
C
C
C     /* Deck stccsm */
      SUBROUTINE STCCSM(NSYMCC)
#include <implicit.h>
#include <priunit.h>
C
#include <inforb.h>
C
      NSYMCC = NSYM
C
      RETURN
      END
C
C
C     /* Deck wrispc */
      SUBROUTINE WRISPC(FREQ,RNNORM,QUBIC,QUARTC,TXT,NCOOR,NDCOOR,NTIME,
     &                  IPRINT)
**************************************************************
*** Writes necessary information to DALTON.SPC in order to ***
*** run the dal2spectro.pl script.                         ***
**************************************************************
#include <implicit.h>
#include <priunit.h>
      CHARACTER*6 TXT
      DIMENSION FREQ(NCOOR), RNNORM(NCOOR), QUBIC(NCOOR,NCOOR,NCOOR),
     &          QUARTC(NCOOR,NCOOR,NCOOR,NCOOR)
C
C     *** Open DALTON.SPC file. ***
C
      LSPECT = 0
      CALL GPOPEN(LSPECT,'DALTON.SPC','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
C
C     *** NTIME = 1 -> write the frequencies and norm of ***
C     ***              the normal coordinates.           ***
C     *** NTIME = 2 -> write the cubic and quartic force ***
C     ***              field.                            ***
C
      IF (NTIME .EQ. 1) THEN
C
C        *** Forwarding to the end of the file if necesary. ***
C
         IF (TXT(1:6).EQ.'cartes') THEN
            KTOT = (NCOOR**3 + 1) + (NCOOR**4 + 1) + 1
            DO I = 1, KTOT
               READ(LSPECT,*)
            END DO
         ELSE
            WRITE (LSPECT,*) TXT(1:6)
         END IF
C
C        *** Writing frequencies. ***
C
         WRITE (LSPECT,*) 'Frequencies'
         DO I = 1, NDCOOR
            WRITE (LSPECT,'(F18.10)') FREQ(I)
         END DO
C
C        *** Writing norm of normal coordinates. ***
C
         WRITE (LSPECT,*) 'Coordinate norm'
         DO I = 1, NDCOOR
            WRITE (LSPECT,'(F18.10)') RNNORM(I)
         END DO
      ELSE IF (NTIME .EQ. 2) THEN
C
C        *** Forwarding to the end of the file if necesary. ***
C
         IF (TXT(1:6).EQ.'normal') THEN
            KTOT = 2*NDCOOR + 3
            DO I = 1, KTOT
               READ(LSPECT,*)
            END DO
         ELSE
            WRITE (LSPECT,*) TXT(1:6)
         END IF
C
C        *** Writing cubic force field. ***
C
         WRITE (LSPECT,*) 'Cubic force field'
         DO 100 K = 1, NDCOOR
         DO 100 J = 1, NDCOOR
         DO 100 I = 1, NDCOOR
            WRITE (LSPECT,*) QUBIC(I,J,K), I, J, K
 100     CONTINUE
C
C        *** Writing quartic force field. ***
C
         WRITE (LSPECT,*) 'Quartic force field'
         DO 200 L = 1, NDCOOR
         DO 200 K = 1, NDCOOR
         DO 200 J = 1, NDCOOR
         DO 200 I = 1, NDCOOR
            WRITE (LSPECT,*) QUARTC(I,J,K,L), I, J, K, L
 200     CONTINUE
      END IF
C
C     *** CLOSING FILE. ***
C
      CALL GPCLOSE(LSPECT,'KEEP')
C
      RETURN
      END
C
C
C     /* Deck runpnt */
      LOGICAL FUNCTION RUNPNT(CLNRGY,IWIDTH,IDIME)
C     *********************************************************
C     *** Subroutine that checks if this is a               ***
C     *** point we need to calculate.                       ***
C     *** There are several criteria:                       ***
C     *** CLNRGY = .FALSE. -> need not to calculate because ***
C     ***                     of symmetry.                  ***
C     *** (ANALZ1 = .TRUE.) & (NMORDR=IWIDTH=3) -> Need only***
C     ***                     the diagonal cubic force and  ***
C     ***                     this point contributes to     ***
C     ***                     F(I,J,K), I ne J ne K.        ***
C     *********************************************************
#include <implicit.h>
#include <priunit.h>
      LOGICAL RNPNT1, CLNRGY
#include <numder.h>
#include <cbinum.h>
C
C     *** Original geometry is always calculated. ***
C
      IF (IDIME.EQ.1) THEN
         RNPNT1 = .TRUE.
      ELSE
C
C        *** Initializing. ***
C
         RNPNT1 = CLNRGY
C
C        *** Is this a ANALZ1 vibrational average. ***
C
         IF (RNPNT1.AND.ANALZ1.AND.NRMCRD.AND.(NMORDR.EQ.IWIDTH)) THEN
C
C           *** NMORDR+NAORDR = 3 -> The forcefield we need to do in    ***
C           ***                      ANALZ1.                            ***
C           *** NMORDR-NAORDR > 1 -> Don't need this point for property ***
C           ***                      derivatives.                       ***
C
            IF ((NAORDR.LT.2) .AND. ((NMORDR+NAORDR).EQ.3)) THEN
               RNPNT1 = .FALSE.
            END IF
         END IF
      END IF
C
      RUNPNT = RNPNT1
C
      RETURN
      END
C
C
C     /* Deck srtins*/
      SUBROUTINE SRTINS(INDSTP,INDTMP,ITMP,NORDR,IPRINT)
C     *****************************************************************
C     *** Subroutine that sorts three indices, where two indices    ***
C     *** are equal, so that the one index not equal the two others ***
C     *** are put first. The indices are returned in INDSTP.        ***
C     *****************************************************************
#include <implicit.h>
#include <priunit.h>
      LOGICAL EQUAL, FOUND
      DIMENSION INDSTP(NORDR), INDTMP(NORDR), ITMP(N ORDR)
C
      FOUND = .FALSE.
C
      DO J = 1, 3
         IF (.NOT.FOUND) THEN
            EQUAL = .FALSE.
C
            DO I = 1, 3
               IF (I.NE.J) THEN
                  IF (INDSTP(J).EQ.INDSTP(I)) THEN
                     EQUAL = .TRUE.
                  END IF
               END IF
            END DO
C
            IF (.NOT.EQUAL) THEN
               IJ = 0
               ITMP(3) = INDSTP(J)
               DO I = 1, 3
                  IF (I.NE.J) THEN
                     IJ = IJ + 1
                     ITMP(IJ) = INDSTP(I)
                  END IF
               END DO
               FOUND = .TRUE.
               CALL ICOPY(3,ITMP,1,INDTMP,1)
            END IF
         END IF
      END DO
C
      IF (.NOT.FOUND) THEN
         CALL ICOPY(3,INDSTP,1,INDTMP,1)
      END IF
C
      RETURN
      END
C
C
C     /* Deck prexce */
      SUBROUTINE PREXCE(EXENG,NSYM,LUPRI)
C     **************************************************
C     *** Subroutine that prints excitation energies ***
C     **************************************************
#include <implicit.h>
#include <codata.h>
#include <cbiexc.h>
      DIMENSION EXENG(NSYM,MXNEXI)
C
      IF (EXCTRP) THEN
         CALL HEADER('Triplet electronic excitation energies',15)
      ELSE
         CALL HEADER ('Singlet electronic excitation energies',15)
      END IF
C
      WRITE (LUPRI,'(14X,A,/,14X,A,/,14X,A)')
     &     ' Sym.   Mode   Frequency    Frequency',
     &     'ex. st.  No.      (au)          (eV)',
     &     '---------------------------------------'
      DO 15 ISYM = 1, NSYM
         DO 14 IEXVAL = 1, NEXCIT(ISYM)
            WRITE (LUPRI,'(16X,I2,6X,I3,2F12.6)')
     &           ISYM,IEXVAL,EXENG(ISYM,IEXVAL),
     &           EXENG(ISYM,IEXVAL)*XTEV
 14      CONTINUE
 15   CONTINUE
C
      WRITE (LUPRI,'(//)')
C
      RETURN
      END
C
C
C     /* Deck chksgn */
      SUBROUTINE CHKSGN(TRLNFV,IPRINT)
C     *************************************************************
C     *** Subroutine that checks if the phase of the components ***
C     *** of the transition dipole moment is correct.           ***
C     *************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
      PARAMETER (THRSH=1.0D-8)
#include <inforb.h>
#include <cbiexc.h>
#include <pvibav.h>
#include <symmet.h>
#include <numder.h>
      DIMENSION TRLNFV(3,NSYM,MXNEXI,NMPCAL)
C
      DO 100 IMPCAL = 1, NMPCAL
         CALL NMCOMP(NFIRST,IMPCAL,NMRDRP,IPRINT)
C
         DO 200 ISYM   = 1, MAXREP+1
         DO 200 IEXVAL = 1, NEXCTB(ISYM)
         DO 200 IC = 1, 3
            TADD  =(TRLNFV(IC,ISYM,IEXVAL,IMPCAL)
     &            + TRLNFV(IC,ISYM,IEXVAL,NFIRST))**2
            TSUB  =(TRLNFV(IC,ISYM,IEXVAL,IMPCAL)
     &            - TRLNFV(IC,ISYM,IEXVAL,NFIRST))**2
C
C           *** Checking if there has been a change of sign. ***
C
            IF (ABS(TRLNFV(IC,ISYM,IEXVAL,1)).GT.THRSH) THEN
               IF (TADD.LT.TSUB) THEN
                  TRLNFV(IC,ISYM,IEXVAL,IMPCAL) =
     &                                   - TRLNFV(IC,ISYM,IEXVAL,IMPCAL)
               END IF
            ELSE
               IF (TSUB.LT.TADD) THEN
                  TRLNFV(IC,ISYM,IEXVAL,IMPCAL) =
     &                                   - TRLNFV(IC,ISYM,IEXVAL,IMPCAL)
               END IF
            END IF
 200     CONTINUE
 100  CONTINUE
C
      RETURN
      END
C
C
C     /* Deck chkccs */
      SUBROUTINE CHKCCS(TRLNFV,IPRINT)
C     *************************************************************
C     *** Subroutine that checks if the phase of the components ***
C     *** of the cc transition dipole moment is correct.        ***
C     *************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (THRESH=1.0D-10)
#include <numder.h>
#include <prpc.h>
      DIMENSION TRLNFV(NPRPC,NMPCAL)
      CHARACTER LABEL*10, LABX*8, LABY*8, LABZ*8, LABU*8
C
      LUPRPCO = -1
      CALL GPOPEN(LUPRPCO,'CC_PRPC_O','UNKNOWN',' ','FORMATTED',
     *            IDUMMY,.FALSE.)
C
      REWIND(LUPRPCO)
      DO IPRPC = 1, NPRPC
C
C           Read in info on property
C
         READ(LUPRPCO,
     *         '(I5,I3,I4,1X,A10,1X,E22.16,4(1X,A8),3(1X,E22.16),3I4)')
     *         IPRPC2,ISYMIN,NORD,LABEL,PROP,
     *         LABX,LABY,LABZ,LABU,FRQY,FRQZ,FRQU,ISYEX,ISPEX,IEX
         IF (IPRPC.NE.IPRPC2) CALL QUIT( 'Strange stuff in CHKCCS')
C
         IF (NORD .EQ.-1) THEN
C
            TRLNFV(IPRPC,1) = ABS(TRLNFV(IPRPC,1))
C
            DO IMPCAL = 2, NMPCAL, 2

               TADD=(TRLNFV(IPRPC,IMPCAL)+TRLNFV(IPRPC,IMPCAL+1))**2
               TSUB=(TRLNFV(IPRPC,IMPCAL)-TRLNFV(IPRPC,IMPCAL+1))**2
C
               IF (ISYMIN .NE.1) THEN
                  IF (TSUB.LT.TADD) THEN
                     TRLNFV(IPRPC,IMPCAL) = -TRLNFV(IPRPC,IMPCAL)
                  END IF
               ELSE
                  IF (ABS(TRLNFV(IPRPC,1)).LT.THRESH) THEN
                     IF (TSUB.LT.TADD) THEN
                        TRLNFV(IPRPC,IMPCAL) = -TRLNFV(IPRPC,IMPCAL)
                     END IF
                  ELSE
                     IF (TADD.LT.TSUB) THEN
                        TRLNFV(IPRPC,IMPCAL) = -TRLNFV(IPRPC,IMPCAL)
                     END IF
                  END IF
               END IF
            END DO
         END IF
      END DO
      CALL GPCLOSE(LUPRPCO,'KEEP')
C
      RETURN
      END
C
C
C     /* Deck nmcomp */
      SUBROUTINE NMCOMP(NFIRST,IMPCAL,NORDR,IPRINT)
C     ***********************************************************
C     *** Subroutine that recognizes the first function value ***
C     *** of a derivative and returns this in NFIRST.         ***
C     ***********************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <trkoor.h>
      LOGICAL DONE
C
C     *** Init. ***
C
      DONE = .FALSE.
C
      IDSTRT = 0
      DO IORDR = 1, NORDR
C
         IADD = 1
         DO I = 1, IORDR
            IADD = IADD*2*(NCOOR-I+1)/I
         END DO
         IDSTRT = IDSTRT + IADD
C
C        *** Test that this is a derivative of order iordr. ***
C        *** We need to remove the first point since this   ***
C        *** does not contribute to any derivative.         ***
C
         IF ((IMPCAL-1.LE.IDSTRT).AND.(.NOT.DONE)) THEN
            DONE = .TRUE.
C
C           *** Number of calculation points for this ***
C           *** derivative                            ***
C
            IPNTS = 2**IORDR
C
C           *** This is point number: ***
C
            IDNUM = IMPCAL-2 - (IDSTRT-IADD)
            KPNT = MOD(IDNUM,IPNTS)
C
C           *** First point for this derivative is: ***
C
            NFIRST = IMPCAL - KPNT
C
         END IF
      END DO
C
      RETURN
      END
C
C
C     /* Deck cke1dr */
      SUBROUTINE CHK1DR(CCPRFV,CCPRDR,TMPCCD,TMPCCF,COEFF,GRIREP,
     &                  WORK,ICNT,IADRSS,IMAX,IMIN,INDSTP,INDTMP,
     &                  IDCOMP,IDDCMP,NCVAL,KDPMTX,ICRIRP,NPPDER,MXCOEF,
     &                  NTYPE,NFINNR,LDPMTX,IFRSTD,NLDPMX,LWORK,IPRINT)
C     **********************************************************
C     *** Subroutine that calculates two other possible      ***
C     *** derivatives for the first derivative of transition ***
C     *** dipole moment. This is to check if there is a sign ***
C     *** problem in the molecular system.                   ***
C     **********************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <prpc.h>
#include <numder.h>
#include <fcsym.h>
#include <trkoor.h>
      DIMENSION CCPRFV(NPRPC,NMPCAL), CCPRDR(NPRPC,NPPDER  ),
     &          TMPCCF(NPRPC,NMPCAL), TMPCCD(NPRPC,NPPDER,2),
     &          COEFF(-MXCOEF:MXCOEF,0:NMRDRP),GRIREP(NGORDR,NGVERT),
     &          WORK(LWORK)
      DIMENSION ICNT(NTYPE), IADRSS(NTYPE), IMAX(NMRDRP), IMIN(NMRDRP),
     &          INDSTP(NTORDR), INDTMP(NTORDR), IDCOMP(NCOOR),
     &          IDDCMP(NCOOR),NCVAL(NCOOR),KDPMTX(LDPMTX,NSTRDR,IFRSTD),
     &          ICRIRP(NCOOR,2)
C
      CHARACTER LABEL*10, LABX*8, LABY*8, LABZ*8, LABU*8
C
C     *** First option, reversing the sign on the last value. ***
C     *** df/dx = f(+)-(-f(-)).                               ***
C
C     *** Assigning new function values. ***
C
      IMPCAL = 1
      DO ICOOR = 1, NCOOR
         IMPCAL = IMPCAL + 1
         DO IPRPC = 1, NPRPC
            TMPCCF(IPRPC,IMPCAL) = CCPRFV(IPRPC,IMPCAL)
         END DO
         IMPCAL = IMPCAL + 1
         DO IPRPC = 1, NPRPC
            TMPCCF(IPRPC,IMPCAL) = -CCPRFV(IPRPC,IMPCAL)
         END DO
      END DO
C
C     *** Finding the new derivative. ***
C
      NFINNR = NPRPC
      CALL NMNDER(TMPCCD(1,1,1),COEFF,TMPCCF,GRIREP,WORK,IADRSS,KDPMTX,
     &            ICRIRP,INDSTP,INDTMP,IDCOMP,IMAX,IMIN,ICNT,NCVAL,
     &            IDDCMP,MXCOEF,1,NMPCAL,NTYPE,NPPDER*NFINNR,NFINNR,
     &            LDPMTX,IFRSTD,NLDPMX,LWORK,.FALSE.)
C
C     *** Second option, two, point formula ***
C     *** df/dx = f(+)-f(0).                ***
C
C     *** Assigning new function values. ***
C
      IMPCAL = 1
      DO ICOOR = 1, NCOOR
         IMPCAL = IMPCAL + 1
         DO IPRPC = 1, NPRPC
            TMPCCF(IPRPC,IMPCAL) = CCPRFV(IPRPC,IMPCAL)
         END DO
         IMPCAL = IMPCAL + 1
         DO IPRPC = 1, NPRPC
            TMPCCF(IPRPC,IMPCAL) = CCPRFV(IPRPC,1)
         END DO
      END DO
C
C     *** Finding the new derivative. ***
C
      NFINNR = NPRPC
      CALL NMNDER(TMPCCD(1,1,2),COEFF,TMPCCF,GRIREP,WORK,IADRSS,
     &            KDPMTX,ICRIRP,INDSTP,INDTMP,IDCOMP,IMAX,IMIN,ICNT,
     &            NCVAL,
     &            IDDCMP,MXCOEF,1,NMPCAL,NTYPE,NPPDER*NFINNR,
     &            NFINNR,LDPMTX,IFRSTD,NLDPMX,LWORK,.FALSE.)
C
      LUPRPCO = -1
      CALL GPOPEN(LUPRPCO,'CC_PRPC_O','UNKNOWN',' ','FORMATTED',
     *            IDUMMY,.FALSE.)
C
      DO ICOOR = 1, NCOOR
         CALL AROUND ('Checking derivative with respect to new' //
     &                ' coordinate')
         WRITE (LUPRI,'(A,I4)') 'Coordinate number', ICOOR
         REWIND(LUPRPCO)
         WRITE (LUPRI,'(36X,A)') ' Best guess      Second choice' //
     &           ' Two point formula'
         DO IPRPC = 1, NPRPC
            READ(LUPRPCO,
     *         '(I5,I3,I4,1X,A10,1X,E22.16,4(1X,A8),3(1X,E22.16),3I4)')
     *         IPRPC2,ISYMIN,NORD,LABEL,PROP,
     *         LABX,LABY,LABZ,LABU,FRQY,FRQZ,FRQU,ISYEX,ISPEX,IEX
            IF (IPRPC.NE.IPRPC2) CALL QUIT( 'Strange stuff in CHK1DR')
            IF (NORD .EQ.-1) THEN
               WRITE(LUPRI,
     &       '(I2,A,A8,A3,F9.6,A,3X,F14.7,3X,F14.7,3X,F14.7)')
     &              ISYMIN,' |<O|',LABX,'|i(',FRQY,')>|',
     &              CCPRDR(IPRPC,ICOOR), TMPCCD(IPRPC,ICOOR,1),
     &              2.0D0*TMPCCD(IPRPC,ICOOR,2)
            END IF
         END DO
      END DO
C
      CALL GPCLOSE(LUPRPCO,'KEEP')
C
      RETURN
      END
C
C
C     /* Deck trprsc*/
      SUBROUTINE T1PRSC(DVAL,CDVAL,SYMCOR,NDIM1,NDERV,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D1=1.0D0, D0=0.0D0)
#include <trkoor.h>
      DIMENSION  DVAL(NDIM1,NDERV), SYMCOR(NCOOR,NCOOR),
     &          CDVAL(NDIM1,NCOOR)
C
      CALL DZERO(CDVAL,NDIM1*NCOOR)
      DO ICOOR2 = 1, NCOOR
      DO ICOOR1 = 1, NCOOR
      DO IINNER = 1, NDIM1
         CDVAL(IINNER,ICOOR1) = CDVAL(IINNER,ICOOR1)
     &                 + SYMCOR(ICOOR1,ICOOR2)*DVAL(IINNER,ICOOR2)
      END DO
      END DO
      END DO
C
      RETURN
      END

C
C
C     /* Deck wrimop */
      SUBROUTINE WRIMOP(FREQ,RNNORM,QUBIC,QUARTC,TXT,NCOOR,NDCOOR,NTIME,
     &                  IPRINT)
**************************************************************
*** Writes necessary information to DALTON.MOP :           ***
*** An operator file to be read in by MidasCpp             ***
*** Ove Christiansen based on Torgeirs WRISPC
*** At this stage it makes only sense with norm. coord.
**************************************************************
#include <implicit.h>
#include <priunit.h>
#include <cbinum.h>
      CHARACTER*6 TXT
      DIMENSION FREQ(NCOOR), RNNORM(NCOOR), QUBIC(NCOOR,NCOOR,NCOOR),
     &          QUARTC(NCOOR,NCOOR,NCOOR,NCOOR)
      LOGICAL USESYM
      SAVE FRQLAR
      USESYM = .TRUE.

      NTOT  = 0
      NWRIT = 0
      THRTRM = 0.0D0
C
C     *** Open DALTON.MOP file. ***
C         Only for normal coordinates right now.
C
      IF (TXT(1:6).EQ.'normal') THEN
C         IF ((NTIME.EQ.2).OR.(NTIME.EQ.1 .AND.
      
         WRITE(LUPRI,'(A)') " WRITE TO MIDAS INTERFACE FILE, "
     &           //"DALTON.MOP "
         LMIDAS = 0
         CALL GPOPEN(LMIDAS,'DALTON.MOP','UNKNOWN',' ','FORMATTED',
     &               IDUMMY,.FALSE.)
      ELSE
         RETURN
      ENDIF
C
C     *** NTIME = 1 -> write the frequencies 
C     *** NTIME = 2 -> write the cubic and quartic force field. 
C
C
      FAC = 0.5D0
      IF (NTIME .EQ. 1) THEN
C
         IF (TXT(1:6).EQ.'normal') THEN
            WRITE (LMIDAS,*) "DALTON_FOR_MIDAS "
C
C           *** Writing frequencies. ***
C
C        WRITE (LMIDAS,*) 'Frequencies'
            DO I = 1, NDCOOR
               WRITE (LMIDAS,'(E23.16,2(1X,I5))')
     *               FAC*FREQ(I)*FREQ(I),I,I
c               WRITE (LMIDAS,'(F18.10)') FREQ(I)
            END DO
C
C           Find largest frequency.
C
            FRQLAR=0.0D0
            DO I = 1, NDCOOR
               IF (ABS(FREQ(I)).GT.FRQLAR) FRQLAR = FREQ(I)
            END DO
C
C         ELSE
C            RETURN
         END IF
C
C
      ELSE IF (NTIME .EQ. 2) THEN
C
C
C        *** Forwarding to the end of the file if necesary. ***

         IF (TXT(1:6).EQ.'normal') THEN
            KTOT = NDCOOR + 1
            DO I = 1, KTOT
               READ(LMIDAS,*)
            END DO
         ELSE
            RETURN
         END IF
C
C        Prepare screening and count of significant terms
C        Threshold: do not write out things that are
C        THRMID times smaller than the larges frequency.
C
         NTOT   = NDCOOR
         NWRIT  = NDCOOR
         THRTRM = FAC*FRQLAR*FRQLAR*THRMID
C        WRITE(LUPRI,'(A,E23.16)') " FRQLAR " ,FRQLAR
C        WRITE(LUPRI,'(A,E23.16)') " FAC    " ,FAC
C        WRITE(LUPRI,'(A,E23.16)') " THRMID " ,THRMID
C        WRITE(LUPRI,'(A,E23.16)') " THRTRM " ,THRTRM  
C
C        *** Writing cubic force field. ***
C        Note scaling with norms and symmetry factors!
C
C         WRITE (LMIDAS,*) 'Cubic force field'
         SUM_QUANT = 0.0D0
         IF (.NOT.USESYM) THEN
            FAC=1.0D0/6.0D0
            DO 100 K = 1, NDCOOR
            DO 100 J = 1, NDCOOR
            DO 100 I = 1, NDCOOR
               QUANT =  FAC*QUBIC(I,J,K)*RNNORM(I)*RNNORM(J)*RNNORM(K)
               NTOT=NTOT+1
               IF (ABS(QUANT).GE.THRTRM) THEN

                  SUM_QUANT=SUM_QUANT+ABS(QUANT)
                  WRITE (LMIDAS,'(E23.16,3(1X,I5))')
     *              QUANT, I, J, K
                  NWRIT=NWRIT+1
               ENDIF
 100        CONTINUE
C
C           *** Writing quartic force field. ***
C
C           WRITE (LMIDAS,*) 'Quartic force field'
            FAC=1.0D0/24.0D0
            DO 200 L = 1, NDCOOR
            DO 200 K = 1, NDCOOR
            DO 200 J = 1, NDCOOR
            DO 200 I = 1, NDCOOR
               QUANT = FAC*QUARTC(I,J,K,L)*RNNORM(I)*RNNORM(J)*
     *           RNNORM(K)*RNNORM(L)
               NTOT=NTOT+1
               IF (ABS(QUANT).GE.THRTRM) THEN
                  SUM_QUANT=SUM_QUANT+ABS(QUANT)
                  WRITE (LMIDAS,'(E23.16,4(1X,I5))')
     *            QUANT,I, J, K, L
                  NWRIT=NWRIT+1
               ENDIF
 200        CONTINUE
         ELSE
C one mode
            FAC=1.0D0/6.0D0
            DO 110 I = 1, NDCOOR
               QUANT =  FAC*QUBIC(I,I,I)*RNNORM(I)*RNNORM(I)*RNNORM(I)
               NTOT=NTOT+1
               IF (ABS(QUANT).GE.THRTRM) THEN
                  SUM_QUANT=SUM_QUANT+ABS(QUANT)
                  WRITE (LMIDAS,'(E23.16,3(1X,I5))')
     *               QUANT, I, I, I
                  NWRIT=NWRIT+1
               ENDIF
 110        CONTINUE
            FAC=1.0D0/24.0D0
            DO 210 I = 1, NDCOOR
               NTOT=NTOT+1
               QUANT = FAC*QUARTC(I,I,I,I)*RNNORM(I)*RNNORM(I)*
     *              RNNORM(I)*RNNORM(I)
               IF (ABS(QUANT).GE.THRTRM) THEN
                  SUM_QUANT=SUM_QUANT+ABS(QUANT)
                  WRITE (LMIDAS,'(E23.16,4(1X,I5))')
     *              QUANT,I,I,I,I
                  NWRIT=NWRIT+1
               ENDIF
 210        CONTINUE
C two mode coupling
            FAC=1.0D0/2.0D0
            DO 120 J = 1, NDCOOR
            DO 120 I = 1, NDCOOR
               IF (I.NE.J) THEN
                  NTOT=NTOT+1
                  QUANT = FAC*QUBIC(I,J,J)*RNNORM(I)
     *                       *RNNORM(J)*RNNORM(J)
                  IF (ABS(QUANT).GE.THRTRM) THEN
                     SUM_QUANT=SUM_QUANT+ABS(QUANT)
                     WRITE (LMIDAS,'(E23.16,3(1X,I5))')
     *                 QUANT, I, J, J

                     NWRIT=NWRIT+1
                  ENDIF
               ENDIF
 120        CONTINUE
            FAC=1.0D0/4.0D0
            DO 220 J = 1, NDCOOR
            DO 220 I = 1, J-1
               NTOT=NTOT+1
               QUANT = FAC*QUARTC(I,I,J,J)*RNNORM(I)*RNNORM(I)*
     *           RNNORM(J)*RNNORM(J)
               IF (ABS(QUANT).GE.THRTRM) THEN
                  SUM_QUANT=SUM_QUANT+ABS(QUANT)
                  WRITE (LMIDAS,'(E23.16,4(1X,I5))')
     *               QUANT, I, I, J, J
                  NWRIT=NWRIT+1
               ENDIF
 220        CONTINUE
            FAC=1.0D0/6.0D0
            DO 221 J = 1, NDCOOR
            DO 221 I = 1, NDCOOR
               IF (I.NE.J) THEN
                  NTOT=NTOT+1
                  QUANT = FAC*QUARTC(I,J,J,J)*RNNORM(I)*RNNORM(J)*
     *                    RNNORM(J)*RNNORM(J)
                  IF (ABS(QUANT).GE.THRTRM) THEN
                     SUM_QUANT=SUM_QUANT+ABS(QUANT)
                     WRITE (LMIDAS,'(E23.16,4(1X,I5))') QUANT,I,J,J,J
                     NWRIT=NWRIT+1
                  ENDIF
               ENDIF
 221        CONTINUE
C three mode coupling
            FAC=1.0D0
            DO 130 K = 1, NDCOOR
            DO 130 J = 1, K-1
            DO 130 I = 1, J-1
               NTOT=NTOT+1
               QUANT = FAC*QUBIC(I,J,K)*RNNORM(I)*RNNORM(J)*RNNORM(K)
               IF (ABS(QUANT).GE.THRTRM) THEN
                  SUM_QUANT=SUM_QUANT+ABS(QUANT)
                  WRITE (LMIDAS,'(E23.16,3(1X,I5))') QUANT,I,J,K
                  NWRIT=NWRIT+1
               ENDIF
 130        CONTINUE
            FAC=1.0D0/2.0D0
            DO 230 K = 1, NDCOOR
            DO 230 J = 1, K-1
            DO 230 I = 1, NDCOOR
               IF ((I.NE.J).AND.(I.NE.K)) THEN
                  NTOT=NTOT+1
                  QUANT = FAC*QUARTC(I,I,J,K)*RNNORM(I)*RNNORM(I)*
     *               RNNORM(J)*RNNORM(K)
                  IF (ABS(QUANT).GE.THRTRM) THEN
                     SUM_QUANT=SUM_QUANT+ABS(QUANT)
                     WRITE (LMIDAS,'(E23.16,4(1X,I5))')
     *               QUANT, I, I, J, K
                     NWRIT=NWRIT+1
                  ENDIF
               ENDIF
 230        CONTINUE
C
C four mode coupling
            FAC=1.0D0
            DO 240 L = 1, NDCOOR
            DO 240 K = 1, L-1
            DO 240 J = 1, K-1
            DO 240 I = 1, J-1
               NTOT=NTOT+1
               QUANT = FAC*QUARTC(I,J,K,L)*RNNORM(I)*RNNORM(J)*
     *           RNNORM(K)*RNNORM(L)
               IF (ABS(QUANT).GE.THRTRM) THEN
                  SUM_QUANT=SUM_QUANT+ABS(QUANT)
                  WRITE (LMIDAS,'(E23.16,4(1X,I5))')
     *             QUANT,I, J, K, L
                  NWRIT=NWRIT+1
               END IF
 240        CONTINUE
         END IF
C
C Count terms if permutation symmetry was not used
C
         NUNCON=0
         DO 250 L = 1, NDCOOR
            DO 260 K = 1, NDCOOR
               DO 270 J = 1, NDCOOR
                  DO 280 I = 1, NDCOOR
                     NUNCON=NUNCON+1
 280              ENDDO
                  NUNCON=NUNCON+1
 270           ENDDO
 260        ENDDO
            NUNCON=NUNCON+1
 250     ENDDO
         WRITE (LUPRI,'(/,A,3(/,A,I8))')
     *    " Force field has been written to Midas Operator File " ,
     *    " Number of terms without use of perm sym =    " , NUNCON,
     *    " Number of terms in total using pert sym =    " , NTOT,
     *    " Number of signficant terms written      =    " , NWRIT
         WRITE (LUPRI,'(A,E20.13,/,A)')
     *    " Only terms with coefficients greater than    " ,THRTRM ,
     *    " is written to operator file "
         WRITE (LUPRI,'(A,E20.13,/,A,/)')
     *    " Sum of absolute values of coefficients:      " ,SUM_QUANT,
     *    " for the anharmonic part                      "
      END IF
C
C     *** CLOSING FILE. ***
C
      CALL GPCLOSE(LMIDAS,'KEEP')
C
      RETURN
      END
C --- end of abander.F ---
