C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
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
C  /* Deck sopinp */
      SUBROUTINE SOPINP(WORD)
C
C     Erik K. Dalskov, Stephan P.A. Sauer 19.2.1996
C     Stephan P. A. Sauer 22.10.2003: SOPPA(CCSD) included
C
C     This subroutine processes additional input options for SOPPA
C     and SOPPA(CCSD) calculations in abacus. This includes the option
C     to carry out Higher RPA calculations and to calculate the W4 term
C     explicitly.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 2)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <abainf.h>
#include <infrsp.h>
C
      DATA TABLE /'.HIRPA', '.SOPW4 '/
C     
C-----------------
C     Set default.
C-----------------
C
      SOPW4  = .FALSE.
      HIRPA  = .FALSE.
C
      ICHANG = 0
      WORD1 = WORD
  100 CONTINUE
      READ (LUCMD, '(A7)') WORD
      CALL UPCASE(WORD)
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GOTO 100
      ELSE IF (PROMPT .EQ. '.') THEN
         ICHANG = ICHANG + 1
         DO 200 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GOTO (1,2), I
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
      END IF
C     
      RETURN
      END
