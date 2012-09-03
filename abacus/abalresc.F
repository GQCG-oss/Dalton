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
      DATA TABLE /'.GAUGEO','.PRINT ','.SELECT','.xXXXXx','.xXXXXx',
     &            '.PARA1S','.PARA1T','.PARA3S','.PARA3T','.DIAM0S',
     &            '.DIAM1S'/
C
      NEWDEF = (WORD .EQ. '*LROPTS')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
 100     CONTINUE
            READ (LUCMD,'(A7)') WORD
cv            write(lupri,*) ' @LRSCINP : reading : ' , WORD
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
               CALL QUIT('Input keywords in LRESC Input, LRSCINP')
C   READ GAUGEO : DEFINIR ALGUNA VAR PARA GAUGEO
 1             CONTINUE
cv                    write(lupri,*)' ..aca leo gaugeo ', WORD
                    READ (LUCMD,*) (LRGAUG(IS), IS = 1, 3)
                    GAUCHANG =.TRUE.
                    ICHANG = ICHANG + 1
               GOTO 100
C   SET PRINT LEVEL 
 2             CONTINUE
cv                  write(lupri,*)' ..aca leo print level '
                  READ (LUCMD,*) JIMPRT
                  ICHANG = ICHANG + 1
               GOTO 100
C   SET WHICH NUC TO DO LRESC 
 3             CONTINUE
cv                  write(lupri,*)' ..aca leo SELECT LRATOM '
                  READ (LUCMD,*) LRATOM
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
            ICHANG = ICHANG - 6 
            LRESCALL=.TRUE.
            Write(lupri,*) 'You asked for all LRESC corrections.'
      ENDIF
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for LRESC :',0)
         IF (GAUCHANG) THEN
           write(lupri,*)' - Remember to set GAUGEO on INTEGRAL ' //
     &     ' SECTION. GAUGEO MUST be placed on selected NUCLEUS'
           WRITE (LUPRI,*) ' - cambiaste origen de gauge .. ojo  '
         END IF
         IF (JIMPRT.NE.0) THEN
            WRITE (LUPRI,'(A,I5)')
     &       ' - Print level in lresc output routine now is : ', JIMPRT
         END IF
         IF (.NOT. LRESCALL) WRITE (LUPRI,'(A)')
     &      ' -  Not LRESC for all corrections. ' //
     &      'Only the following will be done : '
         IF (SIGMAP1S)  WRITE (LUPRI,'(A)')
     &   '  **Paramagnetic first order singlet corrections to shielding'
         IF (SIGMAP1T) WRITE (LUPRI,'(A)')
     &   '  **Paramagnetic first order triplet corrections to shielding'
         IF (SIGMAP3S) WRITE (LUPRI,'(A)')
     &   '  **Paramagnetic third order singlet corrections to shielding'
         IF (SIGMAP3T) WRITE (LUPRI,'(A)')
     &   '  **Paramagnetic third order triplet corrections to shielding'
         IF (SIGMAD0S) WRITE (LUPRI,*)
     &   '  **Diamagnetic zero order singlet corrections to shielding'
         IF (SIGMAD1S) WRITE (LUPRI,*)
     &   '  **Diamagnetic first order singlet corrections to shielding'
         IF (NUCSPI .GT. 0)
     &      WRITE (LUPRI,'(A)')' XXXXXXXXTn couplings of selected '//
     &           'nuclei with all the others will be calculated'
       WRITE(LUPRI,*)' ' 
      END IF
      RETURN
      END
C -----------------------------------------------------------------------
C
C  /* Deck lrscini */
      SUBROUTINE LRSCINI 
C
C     Initialize /LRSCINI/
C
#include <implicit.h>
#include <lrescinf.h>
#include <mxcent.h>
#include <cbiqr.h>
#include <abainf.h>
c jim-dbg : sacar esto de abajo
#include <priunit.h>
C
      JIMPRT   = 0
      LRATOM    = 1 
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
C
C     Initialize /CUADRA/
C
      IPRINT = IPRDEF
      IPRQR  = IPRINT
      SKIP   = .FALSE.
      CUT    = .FALSE.
      OOTV   = .FALSE.
      THRESH = 1.D-04
      MAXITE = 60
      MXRM   = 400
      MXPHP  = 0
cs    LABAPP = 0 ?????????
      NABAPP = 0
      LBFREQ = 1
      LCFREQ = 1
      CALL DZERO (QBFREQ,NFMAX)
      CALL DZERO (QCFREQ,NFMAX) 
C
      RETURN
      END
c    ---------------------------------------
C  /* Deck lrscdrv */
      SUBROUTINE LRSCDRV(WORK,LWORK)
C
#include <implicit.h>
#include <dummy.h>
#include <mxcent.h>
#include <priunit.h>
#include <wrkrsp.h>
#include <lrescinf.h>
C

cv      write(lupri,*)'@ LRSCDRV : '

cv      IF (.NOT. LRESCALL) THEN
cv         WRITE (LUPRI,'(A)')  '   LRESC ALL set : NOT '
cv      END IF
C 
C First Order SINLGET 
C
      IF (SIGMAP1S)  THEN
         WRITE (LUPRI,'(A)')'   Calling to SIGMA1S for : PsoKin'
         WRITE (LUPRI,'(A)')'           NOT FOR NOW, some bugs inside'
c         CALL LINEARS('PSOK',WORK,LWORK) ! NOOOO ANDA
         WRITE (LUPRI,'(A)')'   Calling to SIGMA1S for : Lkin'
         CALL LINEARS('LKIN',WORK,LWORK) ! ANDA
      END IF

      IF (SIGMAD1S) THEN
         WRITE (LUPRI,'(A)')'  Calling to SIGMA1S for : Diamagnetic'// 
     &         ' Singlet First Order Shielding corrections'
         CALL LINEARS('D1S ',WORK,LWORK) ! ANDA
      END IF
C 
C First Order TRIPLET
C
      IF (SIGMAP1T) THEN
         WRITE (LUPRI,'(A)')'  Calling to SIGMA1T for : Paramagnetic'//
     &         ' Triplet First Order Shielding corrections'
         CALL LINEART('FCZK',WORK,LWORK) !ANDA
         CALL LINEART('SDZK',WORK,LWORK) !ANDA
         CALL LINEART('FCBS',WORK,LWORK) !ANDA 
         CALL LINEART('SDBS',WORK,LWORK) !ANDA
      END IF
C 
C Third Order SINLGET 
C
      IF (SIGMAP3S) THEN
         WRITE (LUPRI,'(A)')'   Calling to SIGMA3S for : Paramagnetic'//
     &         ' Singlet Second Order Shielding corrections'
cv         CALL CUADRAS('LPSO')
      END IF
C 
C Third Order TRIPLET 
C
      IF (SIGMAP3T) THEN
         WRITE (LUPRI,'(A)')'   Calling to SIGMA3T for : '
cv         CALL CUADRAS('LFCO')
      END IF

      IF (SIGMAD0S) THEN
         WRITE (LUPRI,'(A)')'   Calling to SIGMA0S for : '
      END IF


      RETURN
      END
C...
c    ---------------------------------------
c    ---------------------------------------

