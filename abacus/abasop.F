C
C...   Copyright (c) 1997 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dalton, an ab initio electronic structure program, Release 1.0
C...   (1997), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, T. Andersen, K. L. Bak, V. Bakken,
C...   O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
C...   H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi, H. Koch,
C...   K. V. Mikkelsen, P. Norman, M. J. Packer, T.Saue,
C...   P. R. Taylor, and O. Vahtras"
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
C  /* Deck sopinp */
      SUBROUTINE SOPINP(WORD)
C
C     Erik K. Dalskov, Stephan P.A. Sauer 19.2.1996
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 7)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <abainf.h>
#include <infrsp.h>
C
      DATA TABLE /'.(CCSD)', 'XXXXXXX', '.SOPW4 ', '.HIRPA ','XXXXXX',
     &            'XXXXXXX', 'XXXXXXX'/
C     
C-----------------
C     Set default.
C-----------------
C
      CCPPA  = .FALSE.
      SOPW4  = .FALSE.
      HIRPA  = .FALSE.
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
               GOTO (1,2,3,4,5,6,7), I
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
            CCPPA  = .TRUE.
            GOTO 100
    2    CONTINUE
            GOTO 100
    3    CONTINUE
            SOPW4  = .TRUE.
            GOTO 100
    4    CONTINUE
            HIRPA  = .TRUE.
            GOTO 100
    5    CONTINUE
            GOTO 100
    6    CONTINUE
            GOTO 100
    7    CONTINUE
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
         IF (CCPPA) WRITE(LUPRI,'(A,A,/,14X,A,/,A,A)')
     &       ' SOPPA(CCSD) :',
     &       '    Second Order Polarization Propagator Approximation',
     &       '    with Coupled Cluster Singles and Doubles Amplitudes',
     &       '    (Ref.: S.P.A. Sauer,',
     &       ' J. Phys. B: At. Mol. Opt. Phys. 30, 3773-3780 (1997))'
         IF (HIRPA)  WRITE(LUPRI,'(A)')
     &       ' HRPA : Higher RPA Polarization Propagator Approximation'
      END IF
C     
      RETURN
      END
