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
C  /* Deck dftinp */
      SUBROUTINE dftINP(WORD,INPERR,ALLOPT)
#include <implicit.h>
#include <dftinf.h>
#include <dftcom.h>
#include <priunit.h>
      PARAMETER ( NTABLE = 10 )
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7, LINE*80
      LOGICAL ALLOPT
      DATA TABLE /'.PRINT ','.DFTELS','.DFTTHR','.RADINT','.ANGINT',
     &            '.ANGMIN','.NOPRUN','.GRID T','.DFTAC ','.CARPAR'/
C
      WORD1 = WORD
C
 100  CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10), I
                  END IF
 200           CONTINUE
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *              '" not recognized in DFTINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in DFTINP.')
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
              WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in DFTINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 CALL QUIT('Illegal keyword in DFTINP.')
            END IF
c     keyword handling section
c           >PRINT < keyword
 1    CONTINUE
      GO TO 100  
c     >.ELECTRONS< check keyword
 2    CONTINUE
      READ(LUCMD,*) DFTELS
      GO TO 100  
c     >.DFTTHR< check keyword
 3    CONTINUE
      READ (LUCMD,*) DFTHR0, DFTHRL, DFTHRI
      GO TO 100  
c     >.RADINT< check keyword
 4    CONTINUE
      READ (LUCMD,*) RADINT
      GO TO 100  
c     >.ANGINT< check keyword
 5    CONTINUE
      READ(LUCMD,*) ANGINT
      GO TO 100  
c     >.ANGMIN< check keyword
 6    CONTINUE
      GO TO 100  
c     >.NOPRUN< check keyword
 7    CONTINUE
      NOPRUN = .TRUE.
      GO TO 100  
c     >.GRID T< grid type keyword
 8    CONTINUE
      READ(LUCMD,'(A80)') LINE
      CALL DFTGRIDINPUT(LINE)
      GO TO 100  
c     >.DFTAC < asymptotic correction keyword
 9    CONTINUE
      DFTASC = .TRUE.
      DFTPOT = .TRUE.
      READ (LUCMD,*) DFTIPT, DFTBR1, DFTBR2
      GO TO 100  
c     >.CARPAR< cartesian grid parameters.
 10   CONTINUE
      READ(LUCMD,'(A80)') LINE
      CALL DFTCARTESIANINPUT(LINE)
      GO TO 100  

 300  CONTINUE

      CALL HEADER('Settings for DFT calculation:',0)
c     code goes here
#if 0
         IF (DFTHRS) THEN
            WRITE (IWUNIT,'(4X,A,17X,3D12.2)')
     &        ' Redefined thresholds:', DFTHR0, DFTHRL, DFTHRI 
         ELSE
            WRITE (IWUNIT,'(4X,A,17X,3D12.2)')
     &        ' Default thresholds:  ', DFTHR0, DFTHRL, DFTHRI 
         END IF
         IF (DFTHRS) THEN
            WRITE (IWUNIT,'(4X,A,17X,3D12.2)')
     &        ' Redefined threshold for number of electrons: ', DFTELS 
         ELSE
            WRITE (IWUNIT,'(4X,A,17X,3D12.2)')
     &        ' Default threshold for number of electrons:   ', DFTELS 
         END IF
         WRITE (IWUNIT,'(4X,A,E10.2,I4)')
     $        ' DFT radial quadrature accuracy/ang. expansion order:',
     &        RADINT,ANGINT
         IF (DFTPOT) WRITE (IWUNIT,'(/,A)')
     &        ' Kohn-Sham matrix calculated from potential.'
         IF (DFTASC) THEN
            WRITE (IWUNIT,'(/,4X,A)') ' Asymptotic correction applied'
            WRITE (IWUNIT,'(10X,A,F10.5)')
     &        ' Ionization potential used: ', DFTIPT      
            WRITE (IWUNIT,'(10X,A,2F10.5)')
     &        ' Radius factors:            ', DFTBR1, DFTBR2      
         END IF 
         IF (DFTHES) THEN
            WRITE (IWUNIT,'(/,4X,A)') 
     &        ' DFT electronic Hessian constructed.'
         END IF
#endif
      END
