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
C jim - May-2012
C
C
      SUBROUTINE LRSCINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <abainf.h>
      PARAMETER (NTABLE = 11, D0 = 0.0D0)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION IPOINT(MXCENT)
#include <dorps.h>
#include <spnout.h>
C jim  aca van las opctiones para LRESC , aun no las defino
      DATA TABLE /'.JIMOPT','.PRINT ','.SELECT','.LINEAR','.QUADRA',
     &            '.GESC  ','.NOPSO ','.SD+FC ','.SDxFC','.ISOTOP',
     &            '.COUPLI'/
C
      write(lupri,*) ' LRESC OPTS :  called with : ' , WORD
      NEWDEF = (WORD .EQ. '*LROPTS')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
 100     CONTINUE
            READ (LUCMD,'(A7)') WORD
            write(lupri,*) ' LRSCOPTS : read : ' , WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GOTO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
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
     &               '" not recognized in SPIN-S.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Input keywords in SPNINP')
 1             CONTINUE
                    write(lupri,*)' i have just read ', WORD
cx                  READ (LUCMD,*) ABUND
cx                  IF (ABUND .EQ. 1.0D0) ICHANG = ICHANG - 1
               GOTO 100
 2             CONTINUE
                  READ (LUCMD,*) ISPPRI
                  IF (ISPPRI .EQ. 0) ICHANG = ICHANG - 1
               GOTO 100
 3             CONTINUE
                  DO 1030 I = 1, MXCENT
                     DOPERT(I,2) = .FALSE.
 1030             CONTINUE
                  READ (LUCMD,*) NPERT
                  READ (LUCMD, *) (IPOINT(I),I=1,NPERT)
                  CALL HEADER('Only the spin-spin couplings between'//
     &               ' the following nuclei will be calculated:',1)
                  WRITE (LUPRI,*) ' ',(IPOINT(I), I = 1 , NPERT)
                  WRITE (LUPRI,'()')
                  DO 1031 I = 1, NPERT
                     DOPERT(IPOINT(I),2) = .TRUE.
 1031             CONTINUE
                  DOSELE = .TRUE.
               GOTO 100
 4             CONTINUE
                  DODSO = .FALSE.
               GOTO 100
 5             CONTINUE
                  DOSD  = .FALSE.
               GOTO 100
 6             CONTINUE
                  DOFC  = .FALSE.
               GOTO 100
 7             CONTINUE
                  DOPSO = .FALSE.
               GOTO 100
 8             CONTINUE
                  DOSDFC= .TRUE.
                  DOSD  = .FALSE.
                  DOFC  = .FALSE.
               GOTO 100
 9             CONTINUE
                  ANISON = .TRUE.
               GOTO 100
 10            CONTINUE
                  SPNISO = .TRUE.
                  READ (LUCMD, *) (ISOTPS(IS), IS = 1, NATOMS)
                  ABUND = D0
               GOTO 100
 11            CONTINUE
                  READ (LUCMD, *)  NUCSPI
                  READ (LUCMD, *) (IPOINT(IS), IS = 1, NUCSPI)
                  IF (NUCSPI .GT. 0) THEN
                     DO IS = 1, MXCENT
                        NCSPNI(IS) = .FALSE.
                     END DO
                  END IF
                  DO IS = 1, NUCSPI
                     IF (IPOINT(IS) .LE. 0 .OR.
     &                   IPOINT(IS) .GT. NUCIND) THEN
                        WRITE (LUPRI,*) 'Non-existent '//
     &                     'symmetry-dependent nucleus specified for'//
     &                     ' .COUPLING NUCLEUS'
                        CALL QUIT('Inconsistent input in *SPIN-SPIN'//
     &                            'module')
                     ELSE
                        NCSPNI(IPOINT(IS)) = .TRUE.
                     END IF
                  END DO
               GOTO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GOTO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     &               '" not recognized in SPIN-S'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT ('Illegal keyword in SPIN-S')
            END IF
      END IF
 300  CONTINUE
c jim ICHANG > 0 cambios en default variables
c jim ICHANG < 0 no cambios

      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for SPIN-S:',0)
         IF (ABUND .NE. 1.0D0) THEN
            WRITE (LUPRI,'(A,F8.3,A)')
     &            ' Spin-spin couplings printed for atoms with '//
     &            'abundance greater than :', ABUND,' %'
         END IF
         IF (ISPPRI .NE. 0) THEN
            WRITE (LUPRI,'(A,I5)')
     &            ' Print level in spin-spin output routine: ', ISPPRI
         END IF
         IF (.NOT. DODSO) WRITE (LUPRI,'(A)')
     &      ' No diamagnetic spin-orbit contribution calculated'
         IF (.NOT. DOSD)  WRITE (LUPRI,'(A)')
     &      ' No spin-dipole contribution calculated'
         IF (.NOT. DOFC) WRITE (LUPRI,'(A)')
     &      ' No Fermi contact contribution calculated'
         IF (.NOT. DOPSO) WRITE (LUPRI,'(A)')
     &      ' No paramagnetic spin-orbit contribution calculated'
         IF (DOSDFC) WRITE (LUPRI,'(A)')
     &      ' Only the sum of SD and FC will be calculated'
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
C     Initialize /QRSRES/
C     Initialize /LRSCINI/
C
#include <implicit.h>
#include <mxcent.h>
#include <cbiqr.h>
#include <abainf.h>
c jim sacar esto de abajo
#include <priunit.h>
C
cs    iprdef= deve essere il default (applicato per entrambi i print level)
cs
      Write(lupri,*)' jim At LRSCINI ' 
      IPRINT = IPRDEF
      IPRQR  = IPRINT
      SKIP   = .FALSE.
        CUT    = .FALSE.
cx      OOTV   = .FALSE.
cx      THRESH = 1.D-04
cx      MAXITE = 60
cx      MXRM   = 400
cx      MXPHP  = 0
cs    LABAPP = 0 ?????????
cx      NABAPP = 0
cx      LBFREQ = 1
cx      LCFREQ = 1
cx      CALL DZERO (QBFREQ,NFMAX)
cx      CALL DZERO (QCFREQ,NFMAX)
C
      RETURN
      END
c    ---------------------------------------
C  /* Deck lrscdrv */
cs      SUBROUTINE LRSCDRV(WORK,LWORK)
cs jim-lresc
cs we should direct from here which lresc claculation will be done
cs  linear, quadratic, gesc depending on what was read on INPUT file  


