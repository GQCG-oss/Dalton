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
C  /* Deck sopinp */
      SUBROUTINE SOPINP(WORD)
C
C     Erik K. Dalskov, Stephan P.A. Sauer 19.2.1996
C     Stephan P. A. Sauer 22.10.2003: SOPPA(CCSD) included
C     Stephan P. A. Sauer 14.03.2006: AO-integral SOPPA included
C
C     This subroutine processes additional input options for SOPPA
C     and SOPPA(CCSD) calculations in abacus. This includes the option
C     to carry out Higher RPA calculations and to calculate the W4 term
C     explicitly.
C
C     It also includes the options to carry out atomic integral direct
C     RPA, RPA(D), SOPPA and SOPPA(CCSD) calculations and several
C     options concerning these calculations.
C     
C     
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 12)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      LOGICAL LNEXC2
#include <abainf.h>
#include "inforb.h"
#include <infrsp.h>
#include <soppinf.h>
C
      DATA TABLE /'.HIRPA ', '.SOPW4 ', '.DIRECT', '.AORPA ', '.DCRPA ',
     &            '.AOSOP ', '.AOSOC ', '.SOPCHK', '.AOTEST', '.NSAVMX',
     &            '.NEXCI2', '.THREX2'/
C     
C------------------
C     Set defaults.
C------------------
C
      SOPW4  = .FALSE.
      HIRPA  = .FALSE.
      AORPA  = .FALSE.
      DCRPA  = .FALSE.
      AOSOP  = .FALSE.
      AOSOC  = .FALSE.
      SOPCHK = .FALSE.
      AOTEST = .FALSE.
C
      NSAVMX = 3
      NSAVMXORIG = 3
      THREX2 = 1.0D-02
      LNEXC2 = .FALSE.
C
      DO 50 I = 1,8
         NEXCI2(I) = 0
   50 CONTINUE
C
      ICHANG = 0
      WORD1 = WORD
  100 CONTINUE
      READ (LUCMD, '(A7)') WORD
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GOTO 100
      ELSE IF (PROMPT .EQ. '.') THEN
         ICHANG = ICHANG + 1
         DO 200 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GOTO (1,2,3,4,5,6,7,8,9,10,11,12), I
            END IF
  200    CONTINUE
         IF (WORD .EQ. '.OPTION') THEN
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            GOTO 100
         END IF
         WRITE (LUPRI,'(/,3A,/)') 'Keyword "', WORD,
     &       '" not recognized in SOPPA.'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword in SOPPA')
    1    CONTINUE
            HIRPA  = .TRUE.
            GOTO 100
    2    CONTINUE
            SOPW4  = .TRUE.
            GOTO 100
    3    CONTINUE
            IF (CCPPA) THEN
               AOSOC = .TRUE.
            ELSE 
               AOSOP = .TRUE.
            ENDIF
            GOTO 100
    4    CONTINUE
            AORPA  = .TRUE.
            GOTO 100
    5    CONTINUE
            DCRPA  = .TRUE.
            AORPA  = .TRUE.
            GOTO 100
    6    CONTINUE
            AOSOP  = .TRUE.
            GOTO 100
    7    CONTINUE
            AOSOC  = .TRUE.
            GOTO 100
    8    CONTINUE
            SOPCHK = .TRUE.
            GOTO 100
    9    CONTINUE
            AOTEST = .TRUE.
            GOTO 100
   10    CONTINUE
            READ (LUCMD,*) NSAVMX
            NSAVMXORIG = NSAVMX
            GOTO 100
   11    CONTINUE
            LNEXC2 = .TRUE.
C           NSYM is set to 8 since NSYM has not yet been set !!
            READ (LUCMD,*) (NEXCI2(I),I=1,NSYM)
            GOTO 100
   12    CONTINUE
            READ (LUCMD,*) THREX2
            GOTO 100
      ELSE IF (PROMPT .EQ. '*') THEN
         GOTO 300
      ELSE
         WRITE (LUPRI,'(/,3A,/)') 'Prompt "', WORD,
     &       '" not recognized in SOPPA.'
      END IF
      CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords', LUPRI)
      CALL QUIT('Illegal prompt in SOPPA')
  300 CONTINUE
C
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults in SOPPA:',0)
         IF (HIRPA)  WRITE(LUPRI,'(A)')
     &       ' HRPA : Higher RPA Polarization Propagator Approximation'
         IF (SOPW4)  WRITE(LUPRI,'(A)')
     &       ' The W4 term of SOPPA is calculated explicitly'
         IF (AORPA)  WRITE(LUPRI,'(A)')
     &       ' AO integral driven RPA'
         IF (DCRPA)  WRITE(LUPRI,'(A/,A,A,/,A)')
     &       ' AO integral driven Double corrected RPA :',
     &       '    (Ref.: O. Christiansen, K. L. Bak, H. Koch and',
     &       ' S. P. A. Sauer',
     &       '     Chem. Phys. Lett. 284, 47-62 (1998))'
         IF (AOSOP)  WRITE(LUPRI,'(A/,A,A,/,A)')
     &       ' AO integral driven SOPPA :',
     &       '    (Ref.: K. L. Bak, H. Koch, J. Oddershede,',
     &       ' O. Christiansen and S. P. A. Sauer',
     &       '     J. Chem. Phys. 112, 4173-4185 (2000))'
         IF (AOSOC)  WRITE(LUPRI,'(A,/,A,A,/,A,A)')
     &       ' AO integral driven SOPPA(CCSD)',
     &       '    (Ref.: H. H. Falden, K. R. Falster-Hansen, ',
     &       'K. L. Bak, S. Rettrup and S. P. A. Sauer,',
     &       '     J. Phys. Chem. A 113, 11995-12012 (2009))'
         IF (SOPCHK) WRITE(LUPRI,'(A)')
     &       ' Explicit construction of E[2] and S[2] for SOPPA'
         IF (AOTEST) WRITE(LUPRI,'(1X,A,A,/,A)') 'Test',
     &       ' orthonormality of trial vectors and compare linear',
     &       ' transformed trial vectors in AO driven SOPPA.'
         IF (NSAVMX .NE. 3) WRITE(LUPRI,'(A,A,I3)') 
     &       ' The maximum number of trial vectors for each ',
     &       'excitation is changed to : ',NSAVMX
         IF (LNEXC2) WRITE(LUPRI,'(A,8I5,A,/,A,E9.2)')
     &       ' The last : ',(NEXCI2(I),I=1,8),' excitation energies',
     &       '    are only converged to : ',THREX2
      END IF
C     
      RETURN
      END
