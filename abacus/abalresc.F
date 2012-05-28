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
C  /* Deck lrscinp */
C This module read nesserary input for LRESC calculations.
C jim-gesc : May-2012
C
C
      SUBROUTINE LRSCINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <lrescinf.h>
#include <nuclei.h>
#include <abainf.h>
      PARAMETER (NTABLE = 11, D0 = 0.0D0)
      LOGICAL NEWDEF 
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION IPOINT(MXCENT)
C#include <dorps.h>
#include <spnout.h>
C jim  aca van las opctiones para LRESC , aun no las defino todas
      DATA TABLE /'.GAUGEO','.PRINT ','.SELECT','.xINEAR','.xUADRA',
     &            '.PARA1S','.PARA1T','.PARA3S','.PARA3T','.DIAM0S',
     &            '.DIAM1S'/
C
      write(lupri,*) ' @LRSCINP to read dal file :' , WORD
      NEWDEF = (WORD .EQ. '*LROPTS')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
 100     CONTINUE
            READ (LUCMD,'(A7)') WORD
            write(lupri,*) '  @LRSCINP : reading : ' , WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GOTO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GOTO (1,2,3,4,5,6,7,8,9,10,11), I
                  END IF
 200           CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GOTO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "', WORD,
     &               '" not recognized in LRESC.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Input keywords in SPNINP')
C   READ GAUGEO : DEFINIR AKLGUNA VAR PARA GAUGEO
 1             CONTINUE
                    write(lupri,*)' ..aca leo gaugeo ', WORD
                    READ (LUCMD, *) (LRGAUG(IS), IS = 1, 3)
                    GAUCHANG =.TRUE.
                    ICHANG = ICHANG + 1
               GOTO 100
C   SET PRINT LEVEL 
 2             CONTINUE
                  write(lupri,*)' ..aca leo print levels '
                  READ (LUCMD,*) JIMPRT
c                  IF (JIMPRT.EQ.0) THEN write(lupri,*)'hacer algo  '  
                  ICHANG = ICHANG + 1
               GOTO 100
C   SET WHICH NUC TO DO LRESC 
 3             CONTINUE
                  write(lupri,*)' ..aca leo SELECT NATOM '
                  READ (LUCMD,*) NATOM
              CALL HEADER('Remember to set GAUGEO on INTEGRAL SECTION'//
     &            ' it is important GAUGEO placed on selected NATOM',1)
                  ICHANG = ICHANG + 1
               GOTO 100
C   SET WHICH NUC TO DO LRESC 
 4             CONTINUE
c                  do algo 
               GOTO 100
 5             CONTINUE
c                  do algo 
               GOTO 100
 6             CONTINUE
                  SIGMAP1S = .TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
               GOTO 100
 7             CONTINUE
                  SIGMAP1T = .TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
               GOTO 100
 8             CONTINUE
                  SIGMAP3S = .TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
               GOTO 100
 9             CONTINUE
                  SIGMAP3T = .TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
               GOTO 100
 10            CONTINUE
                  SIGMAD0S = .TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
               GOTO 100
 11            CONTINUE
                  SIGMAD1S=.TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
               GOTO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GOTO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     &               '" not recognized in LRESC'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT ('Illegal keyword in LRESC')
            END IF
      END IF
 300  CONTINUE
c jim-dbg ICHANG > 0 cambios en default variables
c jim-dbg ICHANG < 0 no cambios

      IF ((SIGMAD0S).AND.(SIGMAD1S).AND.(SIGMAP1S).AND.(SIGMAP1T).AND.
     &    (SIGMAP3S).AND.(SIGMAP3T)) THEN 
            ICHANG = ICHANG -6 
            LRESCALL=.TRUE.
            Write(lupri,*) 'You asked for all LRESC corrections.'
      ENDIF
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for LRESC :',0)
         IF (GAUCHANG) THEN
            WRITE (LUPRI,*) 'cambiaste origen de gauge .. ojo  '
            write(lupri,*) ' a ver ', GAUCHANG
         END IF
         IF (JIMPRT.NE.0) THEN
            WRITE (LUPRI,'(A,I5)')
     &       ' Print level in lresc output routine now is : ', JIMPRT
         END IF
         IF (.NOT. LRESCALL) WRITE (LUPRI,'(A)')
     &      ' Not LRESC for all corrections, will be calculated. ' //
     &      'Only the following will be done, as asked for : '
         IF (SIGMAP1S)  WRITE (LUPRI,'(A)')
     &   '   Paramagnetic first order singlet corrections to shielding'
         IF (SIGMAP1T) WRITE (LUPRI,'(A)')
     &   '   Paramagnetic first order triplet corrections to shielding'
         IF (SIGMAP3S) WRITE (LUPRI,'(A)')
     &   '   Paramagnetic third order singlet corrections to shielding'
         IF (SIGMAP3T) WRITE (LUPRI,'(A)')
     &   '   Paramagnetic third order triplet corrections to shielding'
         IF (SIGMAD0S) WRITE (LUPRI,*)
     &   '   Diamagnetic zero order singlet corrections to shielding'
         IF (SIGMAD1S) WRITE (LUPRI,*)
     &   '   Diamagnetic first order singlet corrections to shielding'
         IF (NUCSPI .GT. 0)
     &      WRITE (LUPRI,'(A)')' The spin-spin couplings of selected '//
     &           'nuclei with all the others will be calculated'
      END IF
      RETURN
      END
c
C  /* Deck lrscini */
      SUBROUTINE LRSCINI 
C
C     Initialize /LRSCINI/
C
#include <implicit.h>
#include <lrescinf.h>
C#include <mxcent.h>
C#include <cbiqr.h
#include <abainf.h>
c jim-dbg : sacar esto de abajo
#include <priunit.h>
C
      Write(lupri,*)'  @LRSCINI for LRESC defaults' 
      JIMPRT   = 0
      NATOM    = 1 
      LRGAUG(1) = 0.000
      LRGAUG(2) = 0.000
      LRGAUG(3) = 0.000
      GAUCHANG = .FALSE.
      SIGMAP1S = .FALSE. 
      SIGMAP1T = .FALSE. 
      SIGMAD1S = .FALSE. 
      SIGMAD0S = .FALSE. 
      SIGMAP3S = .FALSE. 
      SIGMAP3T = .FALSE. 
      LRESCALL = .TRUE.
cx      CALL DZERO (QCFREQ,NFMAX)
C
      RETURN
      END
c    ---------------------------------------

