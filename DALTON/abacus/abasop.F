!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
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
C     RPA, HRPA, RPA(D), SOPPA and SOPPA(CCSD) calculations and several
C     options concerning these calculations.
C     
C   
        use so_info, only: AOSOP, AOSOC, AORPA, AOHRP, AOCC2, DCRPA,
     &                     DCHRP, SDCHR
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 18)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      LOGICAL LNEXC2
#include "abainf.h"
#include "inforb.h"
#include "infrsp.h"
#include "soppinf.h"
C
      CALL QENTER('SOPINP')
C
      DATA TABLE /'.HIRPA ', '.SOPW4 ', '.DIRECT', '.AORPA ', '.RPA(D)',
     &            '.HRPAD', '.AOSOP ', '.AOSOC ', '.SOPCHK', '.AOTEST',
     &            '.NSAVMX', '.NEXCI2', '.THREX2', '.AOHRP ', '.AOCC2 ',
     &            '.HRPA(D', '.DCRPA ', '.SHRPA('/
C
C------------------
C     Set defaults.
C------------------
C
      SOPW4  = .FALSE.
      HIRPA  = .FALSE.
      SOPCHK = .FALSE.
      AOTEST = .FALSE.
C
      AOSOP = .FALSE.
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
               GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18), I
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
    1    CONTINUE !.HIRPA
            HIRPA  = .TRUE.
            GOTO 100
    2    CONTINUE !.SOPW4
            SOPW4  = .TRUE.
            GOTO 100
    3    CONTINUE !.DIRECT
            IF (CCPPA) THEN
               AOSOC = .TRUE.
            ELSE 
               AOSOP = .TRUE.
            ENDIF
            GOTO 100
    4    CONTINUE !.AORPA
            AORPA  = .TRUE.
            GOTO 100
    5    CONTINUE !.RPA(D)
            DCRPA  = .TRUE.
            AORPA  = .TRUE.
            GOTO 100
    6    CONTINUE !.HRPAD
            DCHRP = .TRUE.
            AOHRP = .TRUE.
            GOTO 100
    7    CONTINUE !.AOSOP
            AOSOP  = .TRUE.
            GOTO 100
    8    CONTINUE !.AOSOC
            AOSOC  = .TRUE.
            GOTO 100
    9    CONTINUE !.SOPCHK
            SOPCHK = .TRUE.
            GOTO 100
   10    CONTINUE !.AOTEST
            AOTEST = .TRUE.
            GOTO 100
   11    CONTINUE !.NSAVMX
            READ (LUCMD,*) NSAVMX
            NSAVMXORIG = NSAVMX
            GOTO 100
   12    CONTINUE !.NEXCI2
            LNEXC2 = .TRUE.
C           NSYM is set to 8 since NSYM has not yet been set !!
            READ (LUCMD,*) (NEXCI2(I),I=1,NSYM)
            GOTO 100
   13    CONTINUE !.THREX2
            READ (LUCMD,*) THREX2
            GOTO 100
   14    CONTINUE !.AOHRP
            AOHRP  = .TRUE.
            GOTO 100
   15    CONTINUE !.AOCC2
            AOCC2  = .TRUE.
            GOTO 100
   16    CONTINUE !.HRPA(D(keep old keyword just in case)
            DCHRP  = .TRUE.
            AOHRP  = .TRUE.
            GOTO 100
   17    CONTINUE !.DCRPA (keep old keyword just in case)
            DCRPA  = .TRUE.
            AORPA  = .TRUE.
            GOTO 100
   18    CONTINUE !.SHRPA(
            SDCHR  = .TRUE.
            AOHRP  = .TRUE.
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
         IF (AOHRP)  WRITE(LUPRI,'(A)')
     &       ' AO integral driven Higher RPA (HRPA) :'
         IF (AOSOP)  WRITE(LUPRI,'(A/,A,A,/,A)')
     &       ' AO integral driven SOPPA :',
     &       '    (Ref.: K. L. Bak, H. Koch, J. Oddershede,',
     &       ' O. Christiansen and S. P. A. Sauer',
     &       '     J. Chem. Phys. 112, 4173-4185 (2000))'
CPi
         IF (AOCC2)  WRITE(LUPRI,'(A)')
     &       ' AO integral driven SOPPA(CC2) :'
         IF (DCHRP)  WRITE(LUPRI,'(A)')
     &       ' AO integral driven Double corrected HRPA :'
         IF (SDCHR)  WRITE(LUPRI,'(A)')
     &       ' AO integral driven scaled HRPA(D) :'
Cend-Pi
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
         IF (LNEXC2) WRITE(LUPRI,'(A,8I5,A/A,1P,E9.2)')
     &       ' The last : ',(NEXCI2(I),I=1,8),' excitation energies',
     &       '    are only converged to : ',THREX2
      END IF
C     
      CALL QEXIT('SOPINP')
C
      RETURN
      END
